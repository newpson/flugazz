use std::ops::{Add, Mul};
use nalgebra::Vector2;
// use ndarray::Array2;
use crate::runge_kutta::System;

#[derive(Debug)]
pub struct LiquidDropProblem {
    // area: &'a [Vector2<f64>],

    /// I have no f idea what it is, just a numeric coefficient
    c: f64,

    /// kg/m3
    gas_density: f64,

    /// kg/m3
    liquid_density: f64,

    /// В идеальном газе коэффициент кинематической вязкости ν = η/mn совпадает с коэффициентом диффузии.
    /// 
    /// Вязкость газа.
    mu: f64,

    /// Число Нуссельта
    nu: f64,

    // cells: PhaseGrid

    /// Константа при вычислении WL
    sigma: f64,

    /// Критическое значение WL
    stress_moment_critical: f64,
}

impl LiquidDropProblem {
    pub fn new(gas_density: f64, liquid_density: f64, stress_moment_critical: f64,
               c: f64, mu: f64, nu: f64, sigma: f64) -> Self
    {
        LiquidDropProblem {
            // area: &area,
            c: c,
            gas_density: gas_density,
            liquid_density: liquid_density,
            mu: mu,
            nu: nu,
            // cells: grid,
            sigma: sigma,
            stress_moment_critical: stress_moment_critical
        }
    }
}

trait LiquidDrop<'a>: crate::runge_kutta::Vector
{
    fn position(&'a self) -> &'a Vector2<f64>;
    fn speed(&'a self) -> &'a Vector2<f64>;
    fn diameter3(&'a self) -> f64;
}

#[derive(Clone, Debug)]
pub struct LiquidDropState {
    position: Vector2<f64>,
    speed: Vector2<f64>,
    diameter3: f64,

    accumulated_stress: f64,
}

impl LiquidDropState {
    pub fn new_with_stress(position: &Vector2<f64>,
                       speed: &Vector2<f64>,
                       diameter3: f64,
                       accumulated_stress: f64) -> Self
    {
        LiquidDropState {
            position: position.clone(),
            speed: speed.clone(),
            diameter3: diameter3,
            accumulated_stress: accumulated_stress
        }
    }

    pub fn new(position: &Vector2<f64>,
           speed: &Vector2<f64>,
           diameter3: f64) -> Self
    {
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

const gas_speed: Vector2<f64> = Vector2::<f64>::new(3.0, 10.0);
const gas_pressure_grad: Vector2<f64> = Vector2::<f64>::new(-8.0, 7.0);

impl<'a> System for LiquidDropProblem
{
    type State = LiquidDropState;
    fn integrate(&self, time: f64, state: &Self::State) -> Self::State {
        let alpha = 0.75 * self.c * self.gas_density / self.liquid_density / f64::powf(state.diameter3(), 1.0/3.0);
        // let cell = self.cells.extract_cell(&state.position);
        // // println!("Extracted cell: {cell:?}");
        let drop_mass = 1.0/6.0 * std::f64::consts::PI * state.diameter3() * self.liquid_density;
        // // println!("current drop_mass: {drop_mass}");

        let mut accumulated_stress = 1.0 / (1.43 * f64::powf(state.diameter3, 1.0/3.0) * f64::sqrt(self.liquid_density / self.gas_density) / (gas_speed - state.speed()).magnitude());
        let stress_moment = self.gas_density * (gas_speed - state.speed()).magnitude_squared() * f64::powf(state.diameter3, 1.0/3.0) / self.sigma;
        if stress_moment < self.stress_moment_critical {
            println!("zeroing");
            // make stress drop to zero
            accumulated_stress = -state.accumulated_stress;
        }

        LiquidDropState {
            position: state.speed().clone(),
            speed: alpha * (gas_speed - state.speed()).magnitude() * (gas_speed - state.speed()) - 1.0/self.gas_density * gas_pressure_grad,
            diameter3: -drop_mass / f64::powf(1.0 - drop_mass, 0.75) * self.mu * self.nu / self.liquid_density,
            accumulated_stress: accumulated_stress
        }
    }

    fn should_terminate(&self, time: f64, current: &Self::State) -> bool {
        false
    }

    fn should_branch(&self, time: f64, state: &Self::State) -> Option<Vec<Self::State>> {
        const SIN60: f64 = 0.8660254038;
        const COS60: f64 = 0.5;
        const ROTATE_CCW: nalgebra::Matrix2<f64> = nalgebra::Matrix2::new(
            COS60, -SIN60,
            SIN60, COS60
        );
        const ROTATE_CW: nalgebra::Matrix2<f64> = nalgebra::Matrix2::new(
            COS60, SIN60,
            -SIN60, COS60
        );
        let stress_moment = self.gas_density * (gas_speed - state.speed()).magnitude_squared() * f64::powf(state.diameter3, 1.0/3.0) / self.sigma;
        if stress_moment > self.stress_moment_critical && state.accumulated_stress >= 1.0 {
            let state_ccw_speed = ROTATE_CCW * state.speed;
            let state_cw_speed = ROTATE_CW * state.speed;
            let new_diameter3 = state.diameter3 / 2.0;
            return Some(vec![
                LiquidDropState::new(
                    &state.position,
                    &state_ccw_speed, // NOTE: pi/3
                    new_diameter3,
                ),
                LiquidDropState::new(
                    &state.position,
                    &state_cw_speed, // NOTE: -pi/3
                    new_diameter3
                )]
            );
        }
        None
    }
}
