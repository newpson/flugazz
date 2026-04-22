use drop::liquid_drop::*;
use drop::runge_kutta::*;
use nalgebra::Vector2;
use nalgebra::geometry::Point2;

fn main() {
    let polygon = [
        Point2::new(0.0, 0.0),
        Point2::new(0.0, 30.0),
        Point2::new(30.0, 30.0),
        Point2::new(30.0, 0.0),
        Point2::new(0.0, 0.0),
    ];
    let phase_grid = PhaseGrid::new(&polygon, Vector2::new(10.0, 10.0));
    // Капля воды летит в воздухе при нормальном давлении и температуре 20 градусов Цельсия
    let system = LiquidDropProblem::new(
        1.2041,
        1.8e-05,
        &phase_grid,
        1000.0,
        72.8,
        500.0,
        10.0,
        0.99,
        1.0,
    );

    let initial_state = LiquidDropState::new(
        &Point2::new(0.5, 1.0),
        &Vector2::new(1.0, 1.0),
        // diameter: 0.1mm
        f64::powf(0.1, 3.0),
    );

    let mut stepper = Rk4::new(system, &[initial_state.clone()], 0.0, 1.0, 0.1);
    let result = stepper.integrate();

    println!("Общее число капель: {}", result.len());
    for (i, sequence) in result.iter().enumerate() {
        println!("Капля #{}", i + 1);
        println!("Число состояний: {}", sequence.states().len());
        println!("Начало: {}", sequence.time_begin());
        println!("Конец: {:?}", sequence.time_end());
        for state in sequence.states() {
            println!("  {state:?}");
        }
    }
}
