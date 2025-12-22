use float_cmp::approx_eq;

fn ceil(value: f64) -> f64 {
    let floored = f64::floor(value);
    let ceiled = f64::ceil(value);
    if approx_eq!(f64, value, floored, ulps = 5) {
        return floored;
    }
    ceiled
}

fn floor(value: f64) -> f64 {
    let floored = f64::floor(value);
    let ceiled = f64::ceil(value);
    if approx_eq!(f64, value, ceiled, ulps = 5) {
        return ceiled;
    }
    floored
}