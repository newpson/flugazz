#[derive(Debug)]
enum PointPolygonLocation {
    OUTSIDE,
    INSIDE,
    EDGE,
}

fn point_polygon_location(target: &Vector2<f64>, polygon: &[Vector2<f64>]) -> PointPolygonLocation {
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