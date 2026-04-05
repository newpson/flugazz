use std::ops::{Add, Mul};
use nalgebra::Vector2;
use crate::runge_kutta::System;

#[derive(Debug)]
pub struct LiquidDropProblem {
    /// Плотность газа, мкг/мм^3
    gas_density: f64,

    /// Динамическая вязкость газа, Па*мс
    gas_viscosity: f64,

    /// Скорость однонаправленного потока газа, мм/мс
    gas_speed: Vector2<f64>,

    /// Градиент давления однонаправленного потока газа, Па/мм
    gas_pressure_grad: Vector2<f64>,

    /// Плотность жидкости, мкг/мм^3
    fluid_density: f64,

    /// Коэффициент поверхностного натяжения, мкН/мм
    fluid_surface_tension: f64,

    /// Постоянное число Нуссельта, число
    nusselt: f64,

    /// Критическое значение числа Вебера, число
    weber_critical: f64,

    /// Массовый поток, число из интервала (0,1)
    mass_flow: f64,

    /// Коэффициент сопротивления газовой среды, число
    drag_coefficient: f64,
}

impl LiquidDropProblem {
    pub fn new(gas_density: f64, gas_viscosity: f64, gas_speed: &Vector2<f64>, gas_pressure_grad: &Vector2<f64>,
               fluid_density: f64, fluid_surface_tension: f64, nusselt: f64, weber_critical: f64, mass_flow: f64, c: f64) -> Self
    {
        LiquidDropProblem {
            gas_density: gas_density,
            gas_viscosity: gas_viscosity,
            gas_speed: gas_speed.clone(),
            gas_pressure_grad: gas_pressure_grad.clone(),
            fluid_density: fluid_density,
            fluid_surface_tension: fluid_surface_tension,
            nusselt: nusselt,
            weber_critical: weber_critical,
            mass_flow: mass_flow,
            drag_coefficient: c,
        }
    }
}

trait LiquidDrop<'a>: crate::runge_kutta::Vector
{
    fn position(&'a self) -> &'a Vector2<f64>;
    fn speed(&'a self) -> &'a Vector2<f64>;
    fn diameter3(&'a self) -> f64;
    fn accumulated_stress(&'a self) -> f64;
}

#[derive(Clone, Debug)]
pub struct LiquidDropState {
    /// Радиус-вектор положения капли, мм
    position: Vector2<f64>,

    /// Вектор скорости капли, мм/мс
    speed: Vector2<f64>,

    /// Куб диаметра капли, мм^3
    diameter3: f64,

    /// Накопленное напряжение капли, кГц
    accumulated_stress: f64,
}

impl LiquidDropState {
    pub fn new_with_stress(position: &Vector2<f64>, speed: &Vector2<f64>, diameter3: f64, accumulated_stress: f64) -> Self {
        LiquidDropState {
            position: position.clone(),
            speed: speed.clone(),
            diameter3: diameter3,
            accumulated_stress: accumulated_stress
        }
    }

    pub fn new(position: &Vector2<f64>, speed: &Vector2<f64>, diameter3: f64) -> Self {
        LiquidDropState::new_with_stress(position, speed, diameter3, 0.0)
    }
}

impl<'a> LiquidDrop<'a> for LiquidDropState {
    fn position(&'a self) -> &'a Vector2<f64> {
        &self.position
    }

    fn speed(&'a self) -> &'a Vector2<f64> {
        &self.speed
    }

    fn diameter3(&'a self) -> f64 {
        self.diameter3
    }

    fn accumulated_stress(&'a self) -> f64 {
        self.accumulated_stress
    }
}

impl Mul<f64> for LiquidDropState {
    type Output = Self;
    fn mul(self, rhs: f64) -> Self::Output {
        LiquidDropState {
            position: self.position * rhs,
            speed: self.speed * rhs,
            diameter3: self.diameter3 * rhs,
            accumulated_stress: self.accumulated_stress * rhs,
        }
    }
}

impl Add<LiquidDropState> for LiquidDropState {
    type Output = Self;
    fn add(self, rhs: LiquidDropState) -> Self::Output {
        LiquidDropState {
            position: self.position + rhs.position,
            speed: self.speed + rhs.speed,
            diameter3: self.diameter3 + rhs.diameter3,
            accumulated_stress: self.accumulated_stress + rhs.accumulated_stress,
        }
    }
}

impl crate::runge_kutta::Vector for LiquidDropState {}

/// Skew the vector through its normal forming `angle` between `v` and returned vector
fn skew_transform(v: &Vector2<f64>, angle: f64) -> Vector2<f64> {
    let tan = f64::tan(angle);
    Vector2::new(
        v.x - v.y * tan,
        v.x * tan + v.y
    )
}

impl<'a> System for LiquidDropProblem {
    type State = LiquidDropState;

    fn should_terminate(&self, time: f64, current: &Self::State) -> bool {
        // TODO:
        false
    }

    fn should_branch(&self, _time: f64, state: &Self::State) -> Option<Vec<Self::State>> {
        const BRANCHING_ANGLE: f64 = f64::to_radians(30.0);
        if state.accumulated_stress >= 1.0 {
            let new_diameter3 = state.diameter3 / 2.0;
            return Some(vec![
                LiquidDropState::new(
                    &state.position,
                    &skew_transform(&state.speed, BRANCHING_ANGLE),
                    new_diameter3,
                ),
                LiquidDropState::new(
                    &state.position,
                    &skew_transform(&state.speed, -BRANCHING_ANGLE),
                    new_diameter3
                )]
            );
        }
        None
    }

    fn integrate(&self, _time: f64, state: &Self::State) -> Self::State {
        let evaporation_rate = self.mass_flow / (1.0 - self.mass_flow).powf(0.75);
        let diameter = state.diameter3().powf(1.0/3.0);
        let speed_difference = self.gas_speed - state.speed();
        let relaxation_time = 1.43 * diameter * (self.fluid_density / self.gas_density).sqrt() / (self.gas_speed - state.speed()).magnitude();

        LiquidDropState {
            position: state.speed().clone(),
            speed: 0.75 * self.drag_coefficient * self.gas_density / self.fluid_density / diameter * speed_difference.magnitude() * speed_difference - 1.0/self.gas_density * self.gas_pressure_grad,
            diameter3: -evaporation_rate * self.gas_viscosity * self.nusselt * diameter / self.fluid_density,
            accumulated_stress: 1.0 / relaxation_time
        }
    }

    fn post_integrate(&self, _time: f64, previous_state: &Self::State, new_state: &mut Self::State) {
        let previous_diameter = previous_state.diameter3().powf(1.0/3.0);
        let previous_speed_difference = self.gas_speed - previous_state.speed();
        let weber = self.gas_density * previous_diameter * previous_speed_difference.magnitude_squared() / self.fluid_surface_tension;
        if weber <= self.weber_critical {
            new_state.accumulated_stress = 0.0;
        }
    }
}
