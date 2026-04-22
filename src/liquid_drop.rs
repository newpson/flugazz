use crate::approx;
use crate::runge_kutta::System;
use nalgebra::Vector2;
use nalgebra::geometry::Point2;
use std::ops::{Add, Mul};
use std::ops::Index;

#[derive(Debug)]
pub struct PhaseCell {
    /// Скорость газа, мм/мс
    gas_speed: Vector2<f64>,

    /// Градиент давления газа, Па/мм
    gas_pressure_grad: Vector2<f64>,

    /// Масса пара, мкг
    fluid_mass: f64,
}

impl PhaseCell {
    pub fn new_empty() -> Self {
        Self {
            gas_speed: Vector2::new(0.0, 0.0),
            gas_pressure_grad: Vector2::new(0.0, 0.0),
            fluid_mass: 0.0,
        }
    }
}

#[derive(Debug)]
pub struct PhaseGrid {
    area_polygon: Vec<Point2<f64>>,
    grid_origin: Point2<f64>,
    grid_size: Vector2<usize>,
    cell_size: Vector2<f64>,
    cells: Vec<PhaseCell>,
}

// FIXME: Eq
#[derive(PartialEq, Eq, Debug)]
pub enum Location {
    OUTSIDE,
    INSIDE,
    EDGE,
}

impl PhaseGrid {
    /// Get bounding box of the polygon.
    /// Returns the corners of the bounding box `(bottom_left, top_right)`.
    pub fn get_bounds(polygon: &[Point2<f64>]) -> (Point2<f64>, Point2<f64>) {
        let mut bottom_left = polygon[0];
        let mut top_right = polygon[0];
        for point in polygon.iter() {
            if point.x < bottom_left.x {
                bottom_left.x = point.x
            }
            if point.y < bottom_left.y {
                bottom_left.y = point.y
            }
            if point.x > top_right.x {
                top_right.x = point.x
            }
            if point.y > top_right.y {
                top_right.y = point.y
            }
        }
        (bottom_left, top_right)
    }

    pub fn new(polygon: &[Point2<f64>], cell_size: Vector2<f64>) -> PhaseGrid {
        let (bottom_left, top_right) = Self::get_bounds(polygon);
        let grid_size_float = approx::ceil_vec(top_right - bottom_left).component_div(&cell_size);
        let grid_size = Vector2::new(grid_size_float.x as usize, grid_size_float.y as usize);
        let mut cells = Vec::with_capacity(grid_size.x * grid_size.y);
        for _ in 0..cells.capacity() {
            cells.push(PhaseCell::new_empty());
        }
        PhaseGrid {
            area_polygon: polygon.to_vec(),
            grid_origin: bottom_left,
            grid_size: grid_size,
            cell_size: cell_size,
            cells: cells,
        }
    }

    pub fn sample<'a>(&'a self, point: Point2<f64>) -> Option<&'a PhaseCell> {
        // TODO: sampling from multiple points (linear interpolation insted of nearest)
        let cool = approx::floor_vec((point - self.grid_origin).component_div(&self.cell_size));
        if cool.x.is_sign_negative() || cool.y.is_sign_negative() {
            return None;
        }
        let integer = Vector2::new(cool.x as usize, cool.y as usize);
        if integer.x >= self.grid_size.x || integer.y >= self.grid_size.y {
            return None;
        }
        Some(&self[integer.y][integer.x])
    }

    // Check whether the target point is inside, on the edge or outside the polygon
    pub fn locate(target: &Point2<f64>, polygon: &[Point2<f64>]) -> Location {
        let mut revolutions: i32 = 0;
        let mut half_revolutions: i32 = 0;

        for points_pair in polygon.windows(2) {
            if points_pair[0] == points_pair[1] {
                continue;
            }

            // move `target` to (0, 0)
            let current = points_pair[0] - target;
            let next = points_pair[1] - target;
            let ycomp = current.y * next.y;
            if ycomp <= 0.0 {
                let u = &current;
                let v = next - current;
                let t = -u.component_div(&v);

                // Check if `target` lies on the edge of `polygon`
                // HACK: t component will be NaN if corresponding u and v components will be 0
                if t.x.is_nan() && t.y >= 0.0 && t.y <= 1.0
                    || t.y.is_nan() && t.x >= 0.0 && t.x <= 1.0
                    || t.x == t.y
                {
                    return Location::EDGE;
                }

                // Count winding number of `polygon` around `target`
                let x = u.x + v.x * t.y;
                // HACK: same as above, current.y != next.y
                if !t.y.is_nan() && x > 0.0 {
                    let delta = if next.y > current.y { 1 } else { -1 };
                    if ycomp == 0.0 {
                        half_revolutions += delta;
                    } else {
                        // ycomp < 0
                        revolutions += delta;
                    }
                }
            }
        }

        // NOTE: non-zero rule for winding number
        if revolutions + half_revolutions / 2 + half_revolutions % 2 != 0 {
            return Location::INSIDE;
        }

        Location::OUTSIDE
    }

    // Check whether the target point is inside, on the edge or outside the self.polygon
    pub fn locate_self(&self, target: &Point2<f64>) -> Location {
        Self::locate(target, self.area_polygon.as_slice())
    }
}

#[cfg(test)]
mod tests {
    // TODO: test_locate with more complicated polygon
    use super::*;
    const AREA_POLYGON: [Point2<f64>; 5] = [
        Point2::new(0.0, 0.0),
        Point2::new(0.0, 5.0),
        Point2::new(5.0, 5.0),
        Point2::new(5.0, 0.0),
        Point2::new(0.0, 0.0),
    ];

    const AREA_WIDE_POLYGON: [Point2<f64>; 5] = [
        Point2::new(0.0, 0.0),
        Point2::new(0.0, 5.0),
        Point2::new(10.0, 5.0),
        Point2::new(10.0, 0.0),
        Point2::new(0.0, 0.0),
    ];

    #[test]
    fn test_sample_negative() {
        let phase_grid: PhaseGrid = PhaseGrid::new(
            &AREA_POLYGON,
            Vector2::new(1.0, 1.0),
        );
        assert!(phase_grid.sample(Point2::new(-1.0, -1.0)).is_none());
    }

    #[test]
    fn test_sample_zero() {
        let phase_grid: PhaseGrid = PhaseGrid::new(
            &AREA_POLYGON,
            Vector2::new(1.0, 1.0),
        );
        assert!(phase_grid.sample(Point2::new(0.0, 0.0)).is_some());
    }

    #[test]
    fn test_sample_outofbounds() {
        let phase_grid: PhaseGrid = PhaseGrid::new(
            &AREA_POLYGON,
            Vector2::new(1.0, 1.0),
        );
        assert!(phase_grid.sample(Point2::new(4.99, 5.0)).is_none());
    }

    #[test]
    fn test_sample_wide() {
        let phase_grid: PhaseGrid = PhaseGrid::new(
            &AREA_WIDE_POLYGON,
            Vector2::new(2.0, 1.0),
        );
        assert!(phase_grid.sample(Point2::new(9.99, 4.99)).is_some());
    }

    #[test]
    fn test_sample_wide_outofbounds() {
        let phase_grid: PhaseGrid = PhaseGrid::new(
            &AREA_WIDE_POLYGON,
            Vector2::new(2.0, 1.0),
        );
        assert!(phase_grid.sample(Point2::new(10.0, 0.0)).is_none());
    }

    #[test]
    fn test_locate_edge() {
        assert_eq!(PhaseGrid::locate(&Point2::new(0.0, 0.0), &AREA_WIDE_POLYGON), Location::EDGE);
        assert_eq!(PhaseGrid::locate(&Point2::new(0.0, 3.0), &AREA_WIDE_POLYGON), Location::EDGE);
        assert_eq!(PhaseGrid::locate(&Point2::new(0.01, 5.0), &AREA_WIDE_POLYGON), Location::EDGE);
        assert_eq!(PhaseGrid::locate(&Point2::new(3.0, 5.0), &AREA_WIDE_POLYGON), Location::EDGE);
        assert_eq!(PhaseGrid::locate(&Point2::new(9.99, 5.0), &AREA_WIDE_POLYGON), Location::EDGE);
        assert_eq!(PhaseGrid::locate(&Point2::new(10.0, 5.0), &AREA_WIDE_POLYGON), Location::EDGE);
        assert_eq!(PhaseGrid::locate(&Point2::new(10.0, 4.99), &AREA_WIDE_POLYGON), Location::EDGE);
        assert_eq!(PhaseGrid::locate(&Point2::new(10.0, 0.01), &AREA_WIDE_POLYGON), Location::EDGE);
    }

    #[test]
    fn test_locate_inside() {
        assert_eq!(PhaseGrid::locate(&Point2::new(0.01, 0.01), &AREA_WIDE_POLYGON), Location::INSIDE);
        assert_eq!(PhaseGrid::locate(&Point2::new(0.01, 4.99), &AREA_WIDE_POLYGON), Location::INSIDE);
        assert_eq!(PhaseGrid::locate(&Point2::new(9.99, 0.01), &AREA_WIDE_POLYGON), Location::INSIDE);
        assert_eq!(PhaseGrid::locate(&Point2::new(9.99, 4.99), &AREA_WIDE_POLYGON), Location::INSIDE);
        assert_eq!(PhaseGrid::locate(&Point2::new(2.0, 3.0), &AREA_WIDE_POLYGON), Location::INSIDE);
    }

    #[test]
    fn test_locate_outside() {
        assert_eq!(PhaseGrid::locate(&Point2::new(-0.01, 0.0), &AREA_WIDE_POLYGON), Location::OUTSIDE);
        assert_eq!(PhaseGrid::locate(&Point2::new(-0.01, 5.0), &AREA_WIDE_POLYGON), Location::OUTSIDE);
        assert_eq!(PhaseGrid::locate(&Point2::new(10.01, 0.0), &AREA_WIDE_POLYGON), Location::OUTSIDE);
        assert_eq!(PhaseGrid::locate(&Point2::new(10.01, 5.0), &AREA_WIDE_POLYGON), Location::OUTSIDE);
        assert_eq!(PhaseGrid::locate(&Point2::new(0.0, -0.01), &AREA_WIDE_POLYGON), Location::OUTSIDE);
        assert_eq!(PhaseGrid::locate(&Point2::new(10.0, -0.01), &AREA_WIDE_POLYGON), Location::OUTSIDE);
        assert_eq!(PhaseGrid::locate(&Point2::new(0.0, 5.01), &AREA_WIDE_POLYGON), Location::OUTSIDE);
        assert_eq!(PhaseGrid::locate(&Point2::new(10.0, 5.01), &AREA_WIDE_POLYGON), Location::OUTSIDE);
        assert_eq!(PhaseGrid::locate(&Point2::new(100.0, 200.0), &AREA_WIDE_POLYGON), Location::OUTSIDE);
        assert_eq!(PhaseGrid::locate(&Point2::new(100.0, -200.0), &AREA_WIDE_POLYGON), Location::OUTSIDE);
        assert_eq!(PhaseGrid::locate(&Point2::new(-100.0, 200.0), &AREA_WIDE_POLYGON), Location::OUTSIDE);
        assert_eq!(PhaseGrid::locate(&Point2::new(-100.0, -200.0), &AREA_WIDE_POLYGON), Location::OUTSIDE);
    }
}

impl Index<usize> for PhaseGrid {
    type Output = [PhaseCell];
    fn index(&self, index: usize) -> &Self::Output {
        debug_assert!(index < self.grid_size.y);
        let start = index * self.grid_size.x;
        &self.cells[start..start + self.grid_size.x]
    }
}

#[derive(Debug)]
pub struct LiquidDropProblem {
    /// Плотность газа, мкг/мм^3
    gas_density: f64,

    /// Динамическая вязкость газа, Па*мс
    gas_viscosity: f64,

    /// Сетка газовых ячеек
    phase_grid: PhaseGrid,

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
    pub fn new(
        gas_density: f64,
        gas_viscosity: f64,
        phase_grid: PhaseGrid,
        fluid_density: f64,
        fluid_surface_tension: f64,
        nusselt: f64,
        weber_critical: f64,
        mass_flow: f64,
        c: f64,
    ) -> Self {
        LiquidDropProblem {
            gas_density: gas_density,
            gas_viscosity: gas_viscosity,
            phase_grid: phase_grid,
            fluid_density: fluid_density,
            fluid_surface_tension: fluid_surface_tension,
            nusselt: nusselt,
            weber_critical: weber_critical,
            mass_flow: mass_flow,
            drag_coefficient: c,
        }
    }
}

trait LiquidDrop<'a>: crate::runge_kutta::Linear {
    fn position(&'a self) -> &'a Point2<f64>;
    fn speed(&'a self) -> &'a Vector2<f64>;
    fn diameter3(&'a self) -> f64;
    fn accumulated_stress(&'a self) -> f64;
}

#[derive(Clone, Debug)]
pub struct LiquidDropState {
    /// Радиус-вектор положения капли, мм
    position: Point2<f64>,

    /// Вектор скорости капли, мм/мс
    speed: Vector2<f64>,

    /// Куб диаметра капли, мм^3
    diameter3: f64,

    /// Накопленное напряжение капли, кГц
    accumulated_stress: f64,
}

impl LiquidDropState {
    pub fn new_with_stress(
        position: &Point2<f64>,
        speed: &Vector2<f64>,
        diameter3: f64,
        accumulated_stress: f64,
    ) -> Self {
        LiquidDropState {
            position: position.clone(),
            speed: speed.clone(),
            diameter3: diameter3,
            accumulated_stress: accumulated_stress,
        }
    }

    pub fn new(position: &Point2<f64>, speed: &Vector2<f64>, diameter3: f64) -> Self {
        LiquidDropState::new_with_stress(position, speed, diameter3, 0.0)
    }
}

impl<'a> LiquidDrop<'a> for LiquidDropState {
    fn position(&'a self) -> &'a Point2<f64> {
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
            position: Point2::new(self.position.x + rhs.position.x, self.position.y + rhs.position.y),
            speed: self.speed + rhs.speed,
            diameter3: self.diameter3 + rhs.diameter3,
            accumulated_stress: self.accumulated_stress + rhs.accumulated_stress,
        }
    }
}

impl crate::runge_kutta::Linear for LiquidDropState {}

/// Skew the vector through its normal forming `angle` between `v` and returned vector
fn skew_transform(v: &Vector2<f64>, angle: f64) -> Vector2<f64> {
    let tan = f64::tan(angle);
    Vector2::new(v.x - v.y * tan, v.x * tan + v.y)
}

impl<'a> System for LiquidDropProblem {
    type State = LiquidDropState;

    fn should_terminate(&self, time: f64, current: &Self::State) -> bool {
        match self.phase_grid.locate_self(&current.position) {
            Location::EDGE | Location::OUTSIDE => true,
            Location::INSIDE => false
        }
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
                    new_diameter3,
                ),
            ]);
        }
        None
    }

    fn integrate(&self, _time: f64, state: &Self::State) -> Self::State {
        let evaporation_rate = self.mass_flow / (1.0 - self.mass_flow).powf(0.75);
        let diameter = state.diameter3().powf(1.0 / 3.0);
        let cell = self.phase_grid.sample(state.position).unwrap();
        let speed_difference = cell.gas_speed - state.speed();
        let relaxation_time = 1.43 * diameter * (self.fluid_density / self.gas_density).sqrt()
            / (cell.gas_speed - state.speed()).magnitude();

        LiquidDropState {
            position: Point2::new(state.speed().x, state.speed().y),
            speed: 0.75 * self.drag_coefficient * self.gas_density / self.fluid_density / diameter
                * speed_difference.magnitude()
                * speed_difference
                - 1.0 / self.gas_density * cell.gas_pressure_grad,
            diameter3: -evaporation_rate * self.gas_viscosity * self.nusselt * diameter
                / self.fluid_density,
            accumulated_stress: 1.0 / relaxation_time,
        }
    }

    fn post_integrate(
        &self,
        _time: f64,
        previous_state: &Self::State,
        new_state: &mut Self::State,
    ) {
        let previous_diameter = previous_state.diameter3().powf(1.0 / 3.0);
        let cell = self.phase_grid.sample(previous_state.position).unwrap();
        let previous_speed_difference = cell.gas_speed - previous_state.speed();
        let weber =
            self.gas_density * previous_diameter * previous_speed_difference.magnitude_squared()
                / self.fluid_surface_tension;
        if weber <= self.weber_critical {
            new_state.accumulated_stress = 0.0;
        }
        // TODO: fill the cells with the evaporated fluid mass uniformly
    }
}
