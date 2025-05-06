/*
MIT License with FEFF10 Attribution

Copyright (c) 2025 Ameyanagi

Based on or developed using Distribution: FEFF10.0
Copyright (c) 2020 FEFF Project, University of Washington and SLAC National Accelerator Laboratory.
All rights reserved.
*/

//! Vector3D type for representing 3D positions and directions

use std::fmt;
use std::ops::{Add, Sub};

/// Represents a 3D vector for positions and other spatial quantities
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Vector3D {
    /// X coordinate
    pub x: f64,
    /// Y coordinate
    pub y: f64,
    /// Z coordinate
    pub z: f64,
}

impl Vector3D {
    /// Create a new 3D vector
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    /// Create a new vector at the origin
    pub fn origin() -> Self {
        Self::new(0.0, 0.0, 0.0)
    }

    /// Calculate the distance to another vector
    pub fn distance(&self, other: &Self) -> f64 {
        let dx = self.x - other.x;
        let dy = self.y - other.y;
        let dz = self.z - other.z;
        (dx * dx + dy * dy + dz * dz).sqrt()
    }

    /// Calculate the length (magnitude) of the vector
    pub fn length(&self) -> f64 {
        (self.x * self.x + self.y * self.y + self.z * self.z).sqrt()
    }

    /// Calculate the dot product with another vector
    pub fn dot(&self, other: &Self) -> f64 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    /// Calculate the cross product with another vector
    pub fn cross(&self, other: &Self) -> Self {
        Self {
            x: self.y * other.z - self.z * other.y,
            y: self.z * other.x - self.x * other.z,
            z: self.x * other.y - self.y * other.x,
        }
    }

    /// Normalize the vector to unit length
    pub fn normalize(&self) -> Self {
        let len = self.length();
        if len > 1e-10 {
            Self {
                x: self.x / len,
                y: self.y / len,
                z: self.z / len,
            }
        } else {
            Self::origin()
        }
    }
}

impl fmt::Display for Vector3D {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "({:.6}, {:.6}, {:.6})", self.x, self.y, self.z)
    }
}

impl Add for Vector3D {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Self {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        }
    }
}

impl Sub for Vector3D {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        Self {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_vector_operations() {
        let v1 = Vector3D::new(1.0, 2.0, 3.0);
        let v2 = Vector3D::new(4.0, 5.0, 6.0);

        // Test distance
        assert_relative_eq!(v1.distance(&v2), 5.196152, epsilon = 1e-6);

        // Test length
        assert_relative_eq!(v1.length(), 3.741657, epsilon = 1e-6);

        // Test dot product
        assert_relative_eq!(v1.dot(&v2), 32.0, epsilon = 1e-6);

        // Test cross product
        let cross = v1.cross(&v2);
        assert_relative_eq!(cross.x, -3.0, epsilon = 1e-6);
        assert_relative_eq!(cross.y, 6.0, epsilon = 1e-6);
        assert_relative_eq!(cross.z, -3.0, epsilon = 1e-6);

        // Test normalize
        let norm = v1.normalize();
        assert_relative_eq!(norm.length(), 1.0, epsilon = 1e-6);
    }
}
