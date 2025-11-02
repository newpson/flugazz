use nalgebra::{self, Point2};

type Point = Point2<f64>;

#[derive(Debug)]
enum PointPolygonLocation {
    OUTSIDE,
    INSIDE,
    EDGE,
}

fn calc_point_polygon_location(target: &Point, polygon: &Vec<Point>) -> PointPolygonLocation {
    let mut w: isize = 0;
    let mut hw: isize = 0;
    for points_pair in polygon.windows(2) {
        let current = &points_pair[0];
        let next = &points_pair[1];

        if current == next {
            continue;
        }

        let v = next - current;
        let tx = (target.x - current.x) / v.x;
        let ty = (target.y - current.y) / v.y;
        if v.x == 0.0 {
            if target.x == current.x && ty >= 0.0 && ty <= 1.0 {
                return PointPolygonLocation::EDGE;
            }
        } else if v.y == 0.0 {
            if target.y == current.y && tx >= 0.0 && tx <= 1.0 {
                return PointPolygonLocation::EDGE;
            }
        } else if tx == ty && tx >= 0.0 && tx <= 1.0 {
            return PointPolygonLocation::EDGE;
        }
        
        if current.x >= target.x || next.x >= target.x {
            if current.y < target.y {
                if next.y > target.y {
                    w -= 1;
                } else if next.y == target.y {
                    hw -= 1;
                }
            } else if current.y > target.y {
                if next.y < target.y {
                    w += 1;
                } else if next.y == target.y {
                    hw += 1;
                }
            } else { // current.y == target.y
                if next.y > target.y {
                    hw -= 1;
                } else if next.y < target.y {
                    hw += 1;
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
    let wnumber = calc_point_polygon_location(
        &Point::new(1.5, 2.4),
        &vec![
            Point::new(0.0, 0.0),
            Point::new(0.0, 3.0),
            Point::new(3.0, 2.0),
            Point::new(3.0, 0.0),
            Point::new(0.0, 0.0)]);
    println!("Winding number: {wnumber:?}");
}
