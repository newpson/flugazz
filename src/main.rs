use nalgebra::Vector2;

mod runge_kutta;
use runge_kutta::Rk4;

mod liquid_drop;
use liquid_drop::{LiquidDropProblem, LiquidDropState};

fn main() {
    // Капля воды летит в воздухе с температурой 20 градусов Цельсия
    let system = LiquidDropProblem::new(
        1.2041, 1.8e-05, &Vector2::new(40.0, -80.0), &Vector2::new(0.0, 0.0),
        1000.0, 72.8, 1.0,
        10.0, 0.9, 1.0
    );

    let initial_state = LiquidDropState::new(
        &Vector2::new(0.0, 30.0),
        &Vector2::new(10.0, 10.0),
        // diameter=0.1mm
        f64::powf(0.1, 3.0));

    let mut stepper = Rk4::new(system, &[initial_state.clone()], 0.0, 1.0, 0.1);
    let result = stepper.integrate();

    println!("Общее число капель: {}", result.len());
    for (i, sequence) in result.iter().enumerate() {
        println!("Капля #{}", i+1);
        println!("Число состояний: {}", sequence.states().len());
        println!("Начало: {}", sequence.time_begin());
        println!("Конец: {:?}", sequence.time_end());
        for state in sequence.states() {
            println!("  {state:?}");
        }
    }
}
