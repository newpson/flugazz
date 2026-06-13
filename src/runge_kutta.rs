use std::ops::{Add, Mul};

pub trait Linear: Clone + Add<Self, Output = Self> + Mul<f64, Output = Self> {}

pub trait System: Clone {
    type State: Linear;

    /// Should current state be terminated
    #[allow(unused)]
    fn should_terminate(&self, time: f64, state: &Self::State) -> bool {
        false
    }

    /// Should current state be transformed into multiple new states
    #[allow(unused)]
    fn should_branch(&self, time: f64, state: &Self::State) -> Option<Vec<Self::State>> {
        None
    }

    /// Make one step of integration
    #[allow(unused)]
    fn integrate(&self, time: f64, state: &Self::State) -> Self::State;

    /// Change some properties after integration step if needed
    #[allow(unused)]
    fn post_integrate(&mut self, time: f64, previous_state: &Self::State, new_state: &mut Self::State);
}

#[derive(Debug)]
pub struct StateSequence<TSystem: System> {
    pub time_begin: f64,
    pub time_end: Option<f64>,
    pub states: Vec<TSystem::State>,
}

pub struct Rk4<'r, TSystem: System> {
    /// Implements the set of methods that are used during the integration process and solves a specific problem (model)
    system: &'r mut TSystem,
    /// Current time
    time: f64,
    /// The time (including) when the integration process will be stopped
    time_end: f64,
    /// The difference between two adjacent time points.
    time_step: f64,
    /// Stores all the state sequences (dead ones have `time_end` set).
    storage: Vec<StateSequence<TSystem>>,
    /// Stores the indicies of the state sequences in storage that are currently alive (`time_end` is not set).
    alive: Vec<usize>,
}

impl<'r, TSystem: System> Rk4<'r, TSystem> {
    /// Let's start with a single `process` at moment `t_begin` and run this process with step `t_step` using `system` till `t_end`
    pub fn new(
        system: &'r mut TSystem,
        states: &[TSystem::State],
        time_begin: f64,
        time_end: f64,
        time_step: f64,
    ) -> Self {
        debug_assert!(time_step > 0.0 && time_end >= time_begin);
        let mut storage = Vec::with_capacity(states.len());
        let mut alive = Vec::with_capacity(states.len());
        for (i, state) in states.iter().enumerate() {
            storage.push(StateSequence {
                time_begin: time_begin,
                time_end: None,
                states: vec![state.clone()],
            });
            alive.push(i);
        }
        Rk4 {
            system: system,
            time: time_begin,
            time_end: time_end,
            time_step: time_step,
            storage: storage,
            alive: alive,
        }
    }

    pub fn integrate(&mut self) -> &[StateSequence<TSystem>] {
        while self.time <= self.time_end {
            let mut num_removed: usize = 0;
            let mut i = 0;
            while self.alive.len() > num_removed && i < self.alive.len() - num_removed {
                let last_state = self.storage[self.alive[i]].states.last().unwrap();

                if self.system.should_terminate(self.time, last_state) {
                    self.kill(i);
                    num_removed += 1;
                    continue;
                }

                if let Some(new_states) = self.system.should_branch(self.time, last_state) {
                    for new_state in new_states {
                        self.spawn(new_state);
                    }
                    self.kill(i);
                    num_removed += 1;
                    continue;
                }

                let h = self.time_step;
                let h2 = h / 2.0;
                let h6 = h / 6.0;
                // FIXME: clone() hell;
                let k1 = self.system.integrate(self.time, last_state);
                let k2 = self
                    .system
                    .integrate(self.time + h2, &(last_state.clone() + k1.clone() * h2));
                let k3 = self
                    .system
                    .integrate(self.time + h2, &(last_state.clone() + k2.clone() * h2));
                let k4 = self
                    .system
                    .integrate(self.time + h, &(last_state.clone() + k3.clone() * h));

                let mut new_state = last_state.clone() + (k1 + k2 * 2.0 + k3 * 2.0 + k4) * h6;
                self.system.post_integrate(self.time, last_state, &mut new_state);
                self.storage[self.alive[i]].states.push(new_state);

                i += 1;
            }
            self.time += self.time_step;
            if self.time > self.time_end && self.time - self.time_end < self.time_step / 2.0 {
                self.time = self.time_end;
            }
        }
        self.storage.as_slice()
    }

    fn spawn(&mut self, state: TSystem::State) {
        self.alive.push(self.storage.len());
        self.storage.push(StateSequence {
            time_begin: self.time,
            time_end: None,
            states: vec![state],
        });
    }

    fn kill(&mut self, i_alive: usize) {
        let i_storage = self.alive.swap_remove(i_alive);
        self.storage[i_storage].time_end = Some(self.time);
    }
}
