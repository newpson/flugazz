use nalgebra::{Point2};

type Point = Point2<f64>;

#[derive(Debug)]
enum PointPolygonLocation {
    OUTSIDE,
    INSIDE,
    EDGE,
}

fn calc_point_polygon_location(target: &Point, polygon: &Vec<Point>) -> PointPolygonLocation {
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

fn main() {
    let polygon = vec![
            Point::new(0.0, 0.0),
            Point::new(0.0, 3.0),
            Point::new(3.0, 2.0),
            Point::new(2.0, 0.0),
            Point::new(0.0, 0.0)];
    let step = 0.1;
    let mut point = Point::new(0.0, 0.0);
    for y in (-5..36).rev() {
        point.y = f64::from(y) * step;
        for x in -5..36 {
            point.x = f64::from(x) * step;
            let char = match calc_point_polygon_location(&point, &polygon) {
                PointPolygonLocation::EDGE => '.',
                PointPolygonLocation::INSIDE => '@',
                PointPolygonLocation::OUTSIDE => '_'
            };
            print!("{char}");
        }
        println!("");
    }
}
