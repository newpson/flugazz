use nalgebra::Vector2;

pub fn eq_eps(a: f64, b: f64, eps: f64) -> bool {
    (b - a).abs() < eps
}

pub fn eq(a: f64, b: f64) -> bool {
    const EPS: f64 = 1e-6;
    eq_eps(a, b, EPS)
}

pub fn ceil(value: f64) -> f64 {
    let floored = f64::floor(value);
    let ceiled = f64::ceil(value);
    if eq(value, floored) {
        return floored;
    }
    ceiled
}

pub fn floor(value: f64) -> f64 {
    let floored = f64::floor(value);
    let ceiled = f64::ceil(value);
    if eq(value, ceiled) {
        return ceiled;
    }
    floored
}

pub fn ceil_vec(vec: Vector2<f64>) -> Vector2<f64> {
    Vector2::new(ceil(vec.x), ceil(vec.y))
}

pub fn floor_vec(vec: Vector2<f64>) -> Vector2<f64> {
    Vector2::new(floor(vec.x), floor(vec.y))
}
