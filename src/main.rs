use nalgebra::Vector2;

mod runge_kutta;
use runge_kutta::Rk4;

mod liquid_drop;
use liquid_drop::{LiquidDropProblem, LiquidDropState};

fn main() {
    // let area = [
    //     Vector2::new(-3.0, -2.0),
    //     Vector2::new(-4.0,  9.0),
    //     Vector2::new( 7.0,  6.0),
    //     Vector2::new( 5.0, -2.5),
    //     Vector2::new(-3.0, -2.0)];
    // let grid = PhaseGrid::new_from_polygon(&area, &Vector2::new(0.5, 0.5));

    let system = LiquidDropProblem::new(
        // area: &area,
        0.6, 1000.0, 0.005,
        1.0, 1.596855, 1.0, 1.0,
        // cells: grid,
    );

    let initial_state = LiquidDropState::new(
        &Vector2::new(0.0, 2.0),
        &Vector2::new(1.0, 0.0),
        0.000000001);

    let mut stepper = Rk4::new(system, &[initial_state.clone()], 0.0, 2.0, 0.1);

    let result = stepper.integrate();

    println!("Общее число капель: {}", result.len());
    for (i, sequence) in result.iter().enumerate() {
        println!("Капля #{}", i+1);
        println!("Число состояний: {}", sequence.states().len());
        println!("Начало: {}", sequence.time_begin());
        println!("Конец: {:?}", sequence.time_end());
        for state in sequence.states() {
            println!("  {state:?}");
            // println!("  {},{}", state.position.x, state.position.y);
        }
    }
}

