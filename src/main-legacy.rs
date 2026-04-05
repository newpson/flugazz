use float_cmp::approx_eq;
use nalgebra::{Vector2};
// use ode_solvers::{Vector2, Vector3, Vector5};
use ndarray::Array2;
use ode_solvers::{System, Rk4};
use std::convert::{From, Into};

#[derive(Debug)]
enum PointPolygonLocation {
    OUTSIDE,
    INSIDE,
    EDGE,
}

fn floor_approx(value: f64) -> f64 {
    let floored = f64::floor(value);
    let ceiled = f64::ceil(value);
    if approx_eq!(f64, value, ceiled, ulps = 5) {
        return ceiled;
    }
    floored
}

fn ceil_approx(value: f64) -> f64 {
    let floored = f64::floor(value);
    let ceiled = f64::ceil(value);
    if approx_eq!(f64, value, floored, ulps = 5) {
        return floored;
    }
    ceiled
}

fn calc_point_polygon_location(target: &Vector2<f64>, polygon: &[Vector2<f64>]) -> PointPolygonLocation {
    let mut w: i32 = 0;
    let mut hw: i32 = 0;

    for points_pair in polygon.windows(2) {
        if points_pair[0] == points_pair[1] {
            continue;
        }
        
        // move `target` to (0, 0)
        let current = &points_pair[0] - target;
        let next = &points_pair[1] - target;
        let ycomp = current.y * next.y;
        if ycomp <= 0.0 {
            let u = &current;
            let v = next - current;
            let t = -u.component_div(&v);

            // Check if `target` lies on the edge of `polygon`
            // HACK: t component will be NaN if corresponding u and v components will be 0
            if t.x.is_nan() && t.y >= 0.0 && t.y <= 1.0 ||
               t.y.is_nan() && t.x >= 0.0 && t.x <= 1.0 ||
               t.x == t.y {
                return PointPolygonLocation::EDGE;
            }

            // Count winding number of `polygon` around `target`
            let x = u.x + v.x * t.y;
            // HACK: same as above, current.y != next.y
            if !t.y.is_nan() && x > 0.0 {
                let delta = if next.y > current.y { 1 } else { -1 };
                if ycomp == 0.0 {
                    hw += delta;
                } else { // ycomp < 0
                    w += delta;
                }
            }
        }
    }

    // NOTE: non-zero rule for winding number
    if w + hw/2 + hw%2 != 0 {
        return PointPolygonLocation::INSIDE;
    }
 
    PointPolygonLocation::OUTSIDE
}

#[derive(Debug, Clone)]
struct PhaseCell {
    /// m/s
    gas_speed: Vector2<f64>,
    /// direction of Pa
    gas_pressure_grad: Vector2<f64>,
}

#[derive(Debug)]
struct PhaseGrid {
    origin: Vector2<f64>,
    step: Vector2<f64>,
    size: Vector2<usize>,
    cells: Array2<PhaseCell>
}

impl PhaseGrid {
    fn floor_approx_vector(vector: &Vector2<f64>) -> Vector2<f64> {
        Vector2::new(
            floor_approx(vector.x),
            floor_approx(vector.y)
        )
    }

    fn vector_to_grid(vector: &Vector2<f64>, step: &Vector2<f64>) -> Vector2<usize> {
        let float_vector = PhaseGrid::floor_approx_vector(&vector.component_div(step));
        Vector2::<usize>::new(float_vector.x as usize, float_vector.y as usize)
    }

    fn calc_coords(&self, world_coords: &Vector2<f64>) -> Vector2<usize> {
        PhaseGrid::vector_to_grid(world_coords, &self.step)
    }

    // FIXME: return reference
    fn extract_cell(&self, world_coords: &Vector2<f64>) -> PhaseCell {
        let coords = self.calc_coords(&(world_coords - self.origin));
        // println!("Trying to access coords {coords:?}");
        self.cells[[coords.x, coords.y]].clone()
    }

    /// 1. Определить левую нижнюю точку (min(x), min(y));
    /// 2. Вычесть из всех точек эту точку, чтобы расположить область в I четверти;
    /// 3. Использовать floor_approx(P/h) для определения целочисленных положительных координат точки P;
    fn new_from_polygon(polygon: &[Vector2<f64>], step: &Vector2<f64>) -> PhaseGrid {
        // bottom left corner
        let mut min = polygon[0];
        // top right corner
        let mut max = polygon[0];

        for point in polygon {
            if point.x < min.x { min.x = point.x }
            if point.y < min.y { min.y = point.y }
            if point.x > max.x { max.x = point.x }
            if point.y > max.y { max.y = point.y }
        }

        let size = PhaseGrid::vector_to_grid(&(max - min), &step) + Vector2::<usize>::new(1, 1);

        // FIXME:
        let data = Array2::from_elem(
            (size.y, size.x),
            PhaseCell {
                gas_speed: Vector2::new(10.0, 0.0),
                gas_pressure_grad: Vector2::new(-3.0, 0.0),
            });
        PhaseGrid {
            origin: min,
            step: *step,
            size: size,
            cells: data,
        }
    }
}


#[derive(Debug)]
struct LiquidDropState {
    position: Vector2<f64>,
    speed: Vector2<f64>,
    diameter3: f64,
    
    // fragmentation state
    // accumulated_stress: f64,
    // stress: f64,
    // stress_moment: f64,
}

// some shit for compatibility with ode_solvers
impl Into<ode_solvers::Vector5<f64>> for LiquidDropState {
    fn into(self) -> ode_solvers::Vector5<f64> {
        ode_solvers::Vector5::new(
            self.position.x, self.position.y,
            self.speed.x, self.speed.y,
            self.diameter3)
    }
}

// some shit for compatibility with ode_solvers
impl From<ode_solvers::Vector5<f64>> for LiquidDropState {
    fn from(vector: ode_solvers::Vector5<f64>) -> LiquidDropState {
        LiquidDropState {
            position: Vector2::new(vector[0], vector[1]),
            speed: Vector2::new(vector[2], vector[3]),
            diameter3: vector[4]
        }
    }
}

struct LiquidDropProblem<'a> {
    area: &'a [Vector2<f64>],

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

    // Число Нуссельта
    nu: f64,

    cells: PhaseGrid
}

impl<'a> System<f64, ode_solvers::Vector5<f64>> for LiquidDropProblem<'a> {
    fn system(&self, t: f64, u: &ode_solvers::Vector5<f64>, du: &mut ode_solvers::Vector5<f64>) {
        let state: LiquidDropState = (*u).into();
        let alpha = 0.75 * self.c * self.gas_density / self.liquid_density / f64::powf(state.diameter3, 1.0/3.0);
        let cell = self.cells.extract_cell(&state.position);
        // println!("Extracted cell: {cell:?}");
        let drop_mass = 1.0/6.0 * std::f64::consts::PI * state.diameter3 * self.liquid_density;
        // println!("current drop_mass: {drop_mass}");

        let dspeedx = alpha * (cell.gas_speed - state.speed).magnitude() * (cell.gas_speed.x - state.speed.x) - 1.0/self.gas_density * cell.gas_pressure_grad.x;
        let dspeedy = alpha * (cell.gas_speed - state.speed).magnitude() * (cell.gas_speed.y - state.speed.y) - 1.0/self.gas_density * cell.gas_pressure_grad.y;
        
        let ddiameter3 = -drop_mass / f64::powf(1.0 - drop_mass, 0.75) * self.mu * self.nu / self.liquid_density;
        // println!("diameter3={}, drop_mass={drop_mass}", state.diameter3)1;
        // let ddiameter3 = 0.0;

        let dx = state.speed.x;
        let dy = state.speed.y;

        println!("Current state: {state:?}");
        
        *du = LiquidDropState {
            position: Vector2::new(dx, dy),
            speed: Vector2::new(dspeedx, dspeedy),
            diameter3: ddiameter3,
        }.into()
    }

    fn solout(&mut self, t: f64, u: &ode_solvers::Vector5<f64>, du: &ode_solvers::Vector5<f64>) -> bool {
        let state: LiquidDropState = (*u).into();
        match calc_point_polygon_location(&state.position, self.area) {
            PointPolygonLocation::EDGE => return false,
            PointPolygonLocation::INSIDE => return false,
            PointPolygonLocation::OUTSIDE => {
                println!("Outside: {state:?}");
                return true
            }
        }
    }
}

// fn fragment(drop: &LiquidDropState) -> Option<(LiquidDropState, LiquidDropState)> {
//     const stress_moment_critical: f64 = 1.0;
//     if drop.stress_moment > stress_moment_critical {
//         if drop.accumulated_stress >= 1.0 {
            
//         }
//     }
//     None
// }

fn main() {
    let area = [
        Vector2::new(-3.0, -2.0),
        Vector2::new(-4.0,  9.0),
        Vector2::new( 7.0,  6.0),
        Vector2::new( 5.0, -2.5),
        Vector2::new(-3.0, -2.0)];
    let grid = PhaseGrid::new_from_polygon(&area, &Vector2::new(0.5, 0.5));

    let system = LiquidDropProblem {
        area: &area,
        c: 1.0,
        gas_density: 0.6,
        liquid_density: 1000.0,
        // mu: 1001.596855,
        mu: 1.596855,
        nu: 1.0,
        cells: grid
    };

    let initial_state = LiquidDropState {
        position: Vector2::new(0.0, 2.0),
        speed: Vector2::new(3.0, 0.5),
        diameter3: 0.001
    };

    let mut stepper= Rk4::new(
        system, 0.0, initial_state.into(), 1.0, 0.1);

    let result = stepper.integrate().unwrap();
    let things = stepper.y_out();

    println!("Printing path");
    for state_vec in things {
        let state: LiquidDropState = (*state_vec).into();
        println!("{state:?}");
    }
}
