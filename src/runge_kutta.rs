use std::ops::{Add, Mul};

pub trait Vector: Clone + Add<Self, Output = Self> + Mul<f64, Output = Self> {}

pub trait System {
    type State: Vector;
    type Shared: Default;

    /// Should current state be terminated
    #[allow(unused)]
    fn should_terminate(&self, time: f64, state: &Self::State, shared: &mut Self::Shared) -> bool {
        false
    }

    /// Should current state be transformed into multiple new states
    #[allow(unused)]
    fn should_branch(&self, time: f64, state: &Self::State, shared: &mut Self::Shared) -> Option<Vec<Self::State>> {
        None
    }

    /// Make one step of integration
    #[allow(unused)]
    fn integrate(&self, time: f64, state: &Self::State, shared: &mut Self::Shared) -> Self::State;

    /// Change some properties after integration step if needed
    #[allow(unused)]
    fn post_integrate(&self, time: f64, previous_state: &Self::State, new_state: &mut Self::State, shared: &mut Self::Shared);
}

#[derive(Debug)]
pub struct StateSequence<TSystem: System>  {
    time_begin: f64,
    time_end: Option<f64>,
    states: Vec<TSystem::State>
}

impl<TSystem: System> StateSequence<TSystem> {
    pub fn is_done(&self) -> bool {
        self.time_end.is_some()
    }

    pub fn time_begin(&self) -> f64 {
        self.time_begin
    }

    pub fn time_end(&self) -> Option<f64> {
        self.time_end
    }

    pub fn states(&self) -> &Vec<TSystem::State> {
        &self.states
    }
}

pub struct Rk4<TSystem: System> {
    system: TSystem,
    time: f64,
    time_end: f64,
    time_step: f64,
    /// `(id, state)`; where `id` is index in the `sequence_storage`
    states: Vec<(usize, TSystem::State)>,
    sequence_storage: Vec<StateSequence<TSystem>>,
}

impl<TSystem: System> Rk4<TSystem>
{
    /// Let's start with a single `process` at moment `t_begin` and run this process with step `t_step` using `system` till `t_end`
    pub fn new(
        system: TSystem,
        states: &[TSystem::State],
        time_begin: f64, time_end: f64, time_step: f64
    ) -> Self {
        debug_assert!(time_step > 0.0 && time_end >= time_begin);
        let mut storage = Vec::with_capacity(states.len());
        for state in states {
            storage.push(StateSequence {
                time_begin: time_begin,
                time_end: None,
                states: vec![state.clone()]
            })
        }
        Rk4 {
            system: system,
            time: time_begin,
            time_end: time_end,
            time_step: time_step,
            states: states.iter().cloned().enumerate().collect(),
            sequence_storage: storage
        }
    }

    pub fn integrate(&mut self) -> &Vec<StateSequence<TSystem>> {
        while self.time <= self.time_end {
            let mut num_removed: usize = 0;
            let mut i = 0;
            let mut shared = TSystem::Shared::default();
            
            while i < self.states.len() - num_removed {
                // check removal conditions
                if self.system.should_terminate(self.time, &self.states[i].1, &mut shared) {
                    self.remove(i);
                    num_removed += 1;
                    continue;
                }

                // check branch conditions (when single process splits into multiple)
                if let Some(new_states) = self.system.should_branch(self.time, &self.states[i].1, &mut shared) {
                    // add new states
                    for new_state in new_states.iter() {
                        self.append(new_state);
                    }
                    // and remove the state from which new states branched off
                    self.remove(i);
                    num_removed += 1;
                    continue;
                }

                // and proceed only after all removal operations
                let state = &self.states[i].1;
                let h = self.time_step;
                let h2 = h/2.0;
                let h6 = h/6.0;
                // FIXME: (may be) clone() hell; Mul and Add for references
                let k1 = self.system.integrate(self.time, state, &mut shared);
                let k2 = self.system.integrate(self.time + h2, &(state.clone() + k1.clone() * h2), &mut shared);
                let k3 = self.system.integrate(self.time + h2, &(state.clone() + k2.clone() * h2), &mut shared);
                let k4 = self.system.integrate(self.time + h, &(state.clone() + k3.clone() * h), &mut shared);

                let mut new_state = state.clone() + (k1 + k2*2.0 + k3*2.0 + k4) * h6;
                self.system.post_integrate(self.time, &self.states[i].1, &mut new_state, &mut shared);
                self.states[i].1 = new_state;
                self.sequence_storage[self.states[i].0].states.push(self.states[i].1.clone());

                i += 1;
            }
            self.time += self.time_step;
            if self.time > self.time_end && self.time - self.time_end < self.time_step {
                self.time = self.time_end;
            }
        }
        &self.sequence_storage
    }

    fn append(&mut self, state: &TSystem::State) {
        // add new state into active states
        self.states.push((self.sequence_storage.len(), state.clone()));
        // and create new state sequence
        self.sequence_storage.push(StateSequence { time_begin: self.time, time_end: None, states: vec![state.clone()] });
    }

    fn remove(&mut self, i: usize) {
        // remove from active states
        let j = self.states.swap_remove(i).0;
        // and mark as dead in the storage
        self.sequence_storage[j].time_end = Some(self.time);
    }
}