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
type Rk4System = LiquidDropProblem;


impl From<&StateSequence<LiquidDropProblem>> for LiquidDropStateSequence {
    fn from(sequence: &StateSequence<LiquidDropProblem>) -> Self {
        Self {
            time_begin: sequence.time_begin,
            time_end: sequence.time_end.unwrap_or(0.0),
            is_alive: sequence.time_end.is_none(),
            states_data: sequence.states.as_ptr(),
            states_len: sequence.states.len(),
        }
    }
}

#[unsafe(no_mangle)]
pub unsafe extern "C" fn Rk4System_new(
    gas_density: f64,
    gas_viscosity: f64,
    fluid_density: f64,
    fluid_surface_tension: f64,
    nusselt: f64,
    weber_critical: f64,
    mass_flow: f64,
    drag_coefficient: f64,
    phase_grid: *const PhaseGrid,
) -> *mut Rk4System {
    let phase_grid_ref = unsafe { phase_grid.as_ref().unwrap() };
    let problem_boxed = Box::new(LiquidDropProblem::new(
        gas_density,
        gas_viscosity,
        phase_grid_ref,
        fluid_density,
        fluid_surface_tension,
        nusselt,
        weber_critical,
        mass_flow,
        drag_coefficient,
    ));
    ManuallyDrop::new(problem_boxed).as_mut()
}

#[unsafe(no_mangle)]
pub unsafe extern "C" fn Rk4System_grid(rk4system: *mut Rk4System) -> *mut PhaseGrid {
    let rk4system_ref = unsafe { rk4system.as_mut().unwrap() };
    std::ptr::from_mut(&mut rk4system_ref.phase_grid)
}

#[unsafe(no_mangle)]
pub unsafe extern "C" fn Rk4System_destroy(rk4system: *mut Rk4System) {
    unsafe {
        let _ = Box::from_raw(rk4system);
    }
}

#[unsafe(no_mangle)]
pub unsafe extern "C" fn PhaseGrid_new(
    polygon_data: *const vec2,
    polygon_len: usize,
    cell_size: vec2,
) -> *mut PhaseGrid {
    let polygon_data_transmuted: *const Point2<f64> = unsafe { mem::transmute(polygon_data) };
    let polygon = unsafe { slice::from_raw_parts(polygon_data_transmuted, polygon_len) };
    let phase_grid_boxed = Box::new(PhaseGrid::new(
        polygon,
        Vector2::new(cell_size.x, cell_size.y),
    ));
    ManuallyDrop::new(phase_grid_boxed).as_mut()
}

#[unsafe(no_mangle)]
pub unsafe extern "C" fn PhaseGrid_size(phase_grid: *const PhaseGrid) -> uvec2 {
    let grid_size = unsafe { phase_grid.as_ref().unwrap().grid_size() };
    uvec2 {
        x: grid_size.x,
        y: grid_size.y,
    }
}

#[unsafe(no_mangle)]
pub unsafe extern "C" fn PhaseGrid_origin(phase_grid: *const PhaseGrid) -> vec2 {
    let grid_origin = unsafe { phase_grid.as_ref().unwrap().grid_origin() };
    vec2 {
        x: grid_origin.x,
        y: grid_origin.y,
    }
}

#[unsafe(no_mangle)]
pub unsafe extern "C" fn PhaseGrid_get(
    phase_grid: *const PhaseGrid,
    row: usize,
    col: usize,
) -> *const PhaseCell {
    unsafe { ptr::from_ref(&phase_grid.as_ref().unwrap()[row][col]) }
}

#[unsafe(no_mangle)]
pub unsafe extern "C" fn PhaseGrid_get_mut(
    phase_grid: *mut PhaseGrid,
    row: usize,
    col: usize,
) -> *mut PhaseCell {
    unsafe { ptr::from_mut(&mut phase_grid.as_mut().unwrap()[row][col]) }
}

#[unsafe(no_mangle)]
pub unsafe extern "C" fn PhaseGrid_debug_print(phase_grid: *const PhaseGrid) {
    unsafe {
        println!("{:?}", phase_grid.as_ref().unwrap());
    }
}

#[unsafe(no_mangle)]
pub unsafe extern "C" fn PhaseGrid_destroy(phase_grid: *mut PhaseGrid) {
    unsafe {
        let _ = Box::from_raw(phase_grid);
    }
}

#[unsafe(no_mangle)]
pub unsafe extern "C" fn Rk4_new<'r>(
    system: *mut Rk4System,
    initial_states_data: *const LiquidDropState,
    initial_states_len: usize,
    time_begin: f64,
    time_end: f64,
    time_step: f64,
) -> *mut Rk4<'r, LiquidDropProblem> {
    let system_ref = unsafe { system.as_mut().unwrap() };
    let initial_states_data_transmuted: *const LiquidDropState =
        unsafe { mem::transmute(initial_states_data) };
    let initial_states =
        unsafe { slice::from_raw_parts(initial_states_data_transmuted, initial_states_len) };
    let rk4_boxed = Box::new(Rk4::new(
        system_ref,
        initial_states,
        time_begin,
        time_end,
        time_step,
    ));
    ManuallyDrop::new(rk4_boxed).as_mut()
}

#[unsafe(no_mangle)]
pub unsafe extern "C" fn Rk4_integrate(rk4: *mut Rk4<LiquidDropProblem>) -> *mut Rk4Result {
    let rk4 = unsafe { rk4.as_mut().unwrap() };
    let result = rk4.integrate();
    let mut state_sequences_reprc: Rk4Result = Vec::with_capacity(result.len());
    for sequence in result.iter() {
        state_sequences_reprc.push(LiquidDropStateSequence::from(sequence))
    }
    ManuallyDrop::new(Box::new(state_sequences_reprc)).as_mut()
}

#[unsafe(no_mangle)]
pub unsafe extern "C" fn Rk4_destroy(rk4: *mut Rk4<LiquidDropProblem>) {
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
