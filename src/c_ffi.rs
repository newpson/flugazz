use nalgebra::{Point2, Vector2};
use std::{
    mem::{self, ManuallyDrop},
    ptr, slice,
};

use crate::{
    liquid_drop::{LiquidDropProblem, LiquidDropState, PhaseCell, PhaseGrid},
    runge_kutta::{Rk4, StateSequence},
};

#[repr(C)]
#[allow(non_camel_case_types)]
pub struct vec2 {
    pub x: f64,
    pub y: f64,
}

#[repr(C)]
#[allow(non_camel_case_types)]
pub struct uvec2 {
    pub x: usize,
    pub y: usize,
}

#[repr(C)]
pub struct LiquidDropStateSequence {
    pub time_begin: f64,
    pub time_end: f64,
    pub is_alive: bool,
    pub states_data: *const LiquidDropState,
    pub states_len: usize,
}

type Rk4Result = Vec<LiquidDropStateSequence>;

impl<'a> From<&StateSequence<LiquidDropProblem<'a>>> for LiquidDropStateSequence {
    fn from(value: &StateSequence<LiquidDropProblem<'a>>) -> Self {
        Self {
            time_begin: value.time_begin(),
            time_end: value.time_end().unwrap_or(0.0),
            is_alive: value.is_alive(),
            states_data: value.states().as_ptr(),
            states_len: value.states().len(),
        }
    }
}

#[unsafe(no_mangle)]
pub unsafe extern "C" fn PhaseGrid_new<'a>(
    polygon_data: *const vec2,
    polygon_len: usize,
    cell_size: vec2,
) -> *mut PhaseGrid<'a> {
    let polygon_data_transmuted: *const Point2<f64> = unsafe { mem::transmute(polygon_data) };
    let polygon = unsafe { slice::from_raw_parts(polygon_data_transmuted, polygon_len) };
    let phase_grid_boxed = Box::new(PhaseGrid::new(
        polygon,
        Vector2::new(cell_size.x, cell_size.y),
    ));
    ManuallyDrop::new(phase_grid_boxed).as_mut()
}

#[unsafe(no_mangle)]
pub unsafe extern "C" fn PhaseGrid_get<'a>(
    phase_grid: *mut PhaseGrid<'a>,
    row: usize,
    col: usize,
) -> *mut PhaseCell {
    unsafe { ptr::from_mut(&mut phase_grid.as_mut().unwrap()[row][col]) }
}

#[unsafe(no_mangle)]
pub unsafe extern "C" fn PhaseGrid_debug_print<'a>(phase_grid: *const PhaseGrid<'a>) {
    unsafe {
        println!("{:?}", phase_grid.as_ref().unwrap());
    }
}

#[unsafe(no_mangle)]
pub unsafe extern "C" fn PhaseGrid_destroy<'a>(phase_grid: *mut PhaseGrid<'a>) {
    unsafe {
        let _ = Box::from_raw(phase_grid);
    }
}

#[unsafe(no_mangle)]
pub unsafe extern "C" fn Rk4_new<'a>(
    gas_density: f64,
    gas_viscosity: f64,
    fluid_density: f64,
    fluid_surface_tension: f64,
    nusselt: f64,
    weber_critical: f64,
    mass_flow: f64,
    drag_coefficient: f64,
    phase_grid: *mut PhaseGrid<'a>,
    initial_states_data: *const LiquidDropState,
    initial_states_len: usize,
    time_begin: f64,
    time_end: f64,
    time_step: f64,
) -> *mut Rk4<LiquidDropProblem<'a>> {
    let phase_grid_ref = unsafe { phase_grid.as_mut().unwrap() };
    let problem = LiquidDropProblem::new(
        gas_density,
        gas_viscosity,
        phase_grid_ref,
        fluid_density,
        fluid_surface_tension,
        nusselt,
        weber_critical,
        mass_flow,
        drag_coefficient,
    );
    let initial_states_data_transmuted: *const LiquidDropState =
        unsafe { mem::transmute(initial_states_data) };
    let initial_states =
        unsafe { slice::from_raw_parts(initial_states_data_transmuted, initial_states_len) };
    let rk4_boxed = Box::new(Rk4::new(
        problem,
        initial_states,
        time_begin,
        time_end,
        time_step,
    ));
    ManuallyDrop::new(rk4_boxed).as_mut()
}

#[unsafe(no_mangle)]
pub unsafe extern "C" fn Rk4_integrate<'a>(
    rk4: *mut Rk4<LiquidDropProblem<'a>>,
) -> *mut Rk4Result {
    let rk4 = unsafe { rk4.as_mut().unwrap() };
    let result = rk4.integrate();
    let mut state_sequences_reprc: Rk4Result = Vec::with_capacity(result.len());
    for sequence in result.iter() {
        state_sequences_reprc.push(LiquidDropStateSequence::from(sequence))
    }
    ManuallyDrop::new(Box::new(state_sequences_reprc)).as_mut()
}

#[unsafe(no_mangle)]
pub unsafe extern "C" fn Rk4_destroy<'a>(rk4: *mut Rk4<LiquidDropProblem<'a>>) {
    unsafe {
        let _ = Box::from_raw(rk4);
    }
}

#[unsafe(no_mangle)]
pub unsafe extern "C" fn Rk4Result_data(
    rk4result: *const Rk4Result,
) -> *const LiquidDropStateSequence {
    unsafe { rk4result.as_ref().unwrap().as_ptr() }
}

#[unsafe(no_mangle)]
pub unsafe extern "C" fn Rk4Result_len(rk4result: *const Rk4Result) -> usize {
    unsafe { rk4result.as_ref().unwrap().len() }
}

#[unsafe(no_mangle)]
pub unsafe extern "C" fn Rk4Result_destroy(rk4result: *mut Rk4Result) {
    unsafe {
        let _ = Box::from_raw(rk4result);
    }
}
