use std::ops::{self, Sub};

#[derive(Debug, Copy, Clone)]
pub struct Point {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Point {
    pub fn dot(&self, other: &Point) -> f64 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    pub fn distance(&self, other: &Point) -> f64 {
        let ans =
            (self.x - other.x).powi(2) + (self.y - other.y).powi(2) + (self.z - other.z).powi(2);
        ans.sqrt()
    }

    pub fn cos(&self, p1: &Point, p2: &Point) -> f64 {
        let d1 = self.distance(p1);
        let d2 = self.distance(p2);
        let disp_1 = self.sub(*p1);
        let disp_2 = self.sub(*p2);
        let dot_val = disp_1.dot(&disp_2);
        let ans = dot_val / (d1 * d2);
        ans
    }

    pub fn from_vec(vec: Vec<f64>) -> Point {
        Point {
            x: vec[0],
            y: vec[1],
            z: vec[2],
        }
    }
}

// Operator overloading
impl ops::Add<Point> for Point {
    type Output = Point;

    fn add(self, other: Point) -> Point {
        Point {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        }
    }
}

impl ops::Sub<Point> for Point {
    type Output = Point;

    fn sub(self, other: Point) -> Point {
        Point {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
        }
    }
}

impl ops::Neg for Point {
    type Output = Point;

    fn neg(self) -> Point {
        Point {
            x: -self.x,
            y: -self.y,
            z: -self.z,
        }
    }
}
