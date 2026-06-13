#pragma once
#include <cstddef>

struct PhaseGrid;
struct Rk4;
struct Rk4Result;
struct Rk4System;

struct vec2 {
    double x;
    double y;
};

struct uvec2 {
    std::size_t x;
    std::size_t y;
};

// feature = "c_compatible"
struct PhaseCell {
    vec2 gas_speed;
    vec2 gas_pressure_grad;
    double fluid_mass;
};

// feature = "c_compatible"
struct LiquidDropState {
    vec2 position;
    vec2 speed;
    double diameter3;
    double accumulated_stress;
};

struct LiquidDropStateSequence {
    double time_begin;
    double time_end;
    bool is_alive;
    LiquidDropState *states_data;
    std::size_t states_len;
};

extern "C" {
PhaseGrid *PhaseGrid_new(
    const vec2 *const polygon_data,
    const std::size_t polygon_len,
    const vec2 cell_size);
uvec2 PhaseGrid_size(const PhaseGrid *const phase_grid);
vec2 PhaseGrid_origin(const PhaseGrid *const phase_grid);
const PhaseCell *PhaseGrid_get(const PhaseGrid *const phase_grid, const size_t row, const size_t col);
PhaseCell *PhaseGrid_get_mut(PhaseGrid *const phase_grid, const size_t row, const size_t col);
void PhaseGrid_destroy(PhaseGrid *phase_grid);
void PhaseGrid_debug_print(const PhaseGrid *const phase_grid);

Rk4System *Rk4System_new(
    const double gas_density,
    const double gas_viscosity,
    const double fluid_density,
    const double fluid_surface_tension,
    const double nusselt,
    const double weber_critical,
    const double mass_flow,
    const double drag_coefficient,
    PhaseGrid *const phase_grid);
PhaseGrid *Rk4System_grid(Rk4System *system);
void Rk4System_destroy(Rk4System *system);

Rk4 *Rk4_new(
    Rk4System *const system,
    const LiquidDropState *const initial_states_data,
    const std::size_t initial_states_len,
    const double time_begin,
    const double time_end,
    const double time_step);
Rk4Result *Rk4_integrate(Rk4 *const rk4);
void Rk4_destroy(Rk4 *rk4);

const LiquidDropStateSequence *Rk4Result_data(const Rk4Result *const rk4result);
size_t Rk4Result_len(const Rk4Result *const rk4result);
void Rk4Result_destroy(Rk4Result *result);
}
