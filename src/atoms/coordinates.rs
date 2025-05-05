/*
MIT License with FEFF10 Attribution

Copyright (c) 2025 Ameyanagi

Based on or developed using Distribution: FEFF10.0
Copyright (c) 2020 FEFF Project, University of Washington and SLAC National Accelerator Laboratory.
All rights reserved.
*/

//! Coordinate system conversions for atomic positions

use super::errors::{AtomError, Result};
use super::vector::Vector3D;
use std::f64::consts::PI;

/// Represents a coordinate system for atomic positions
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum CoordinateSystem {
    /// Cartesian coordinates (x, y, z)
    Cartesian,
    /// Spherical coordinates (r, θ, φ)
    Spherical,
    /// Cylindrical coordinates (ρ, φ, z)
    Cylindrical,
}

/// Convert from one coordinate system to another
pub fn convert(
    coords: &[f64; 3],
    from_system: CoordinateSystem,
    to_system: CoordinateSystem,
) -> Result<[f64; 3]> {
    // First convert to Cartesian
    let cartesian = match from_system {
        CoordinateSystem::Cartesian => [coords[0], coords[1], coords[2]],
        CoordinateSystem::Spherical => {
            let r = coords[0];
            let theta = coords[1] * PI / 180.0; // Convert to radians
            let phi = coords[2] * PI / 180.0;
            [
                r * theta.sin() * phi.cos(),
                r * theta.sin() * phi.sin(),
                r * theta.cos(),
            ]
        }
        CoordinateSystem::Cylindrical => {
            let rho = coords[0];
            let phi = coords[1] * PI / 180.0; // Convert to radians
            let z = coords[2];
            [rho * phi.cos(), rho * phi.sin(), z]
        }
    };

    // Then convert from Cartesian to the target system
    match to_system {
        CoordinateSystem::Cartesian => Ok(cartesian),
        CoordinateSystem::Spherical => {
            let x = cartesian[0];
            let y = cartesian[1];
            let z = cartesian[2];
            let r = (x * x + y * y + z * z).sqrt();

            if r < 1e-10 {
                return Err(AtomError::CoordinateError(
                    "Cannot convert to spherical coordinates at origin".to_string(),
                ));
            }

            let theta = (z / r).acos() * 180.0 / PI; // To degrees
            let phi = y.atan2(x) * 180.0 / PI; // To degrees
            Ok([r, theta, phi])
        }
        CoordinateSystem::Cylindrical => {
            let x = cartesian[0];
            let y = cartesian[1];
            let z = cartesian[2];
            let rho = (x * x + y * y).sqrt();
            let phi = y.atan2(x) * 180.0 / PI; // To degrees
            Ok([rho, phi, z])
        }
    }
}

/// Create a Vector3D from coordinates in a specific system
pub fn vector_from_coords(coords: &[f64; 3], system: CoordinateSystem) -> Result<Vector3D> {
    let cartesian = convert(coords, system, CoordinateSystem::Cartesian)?;
    Ok(Vector3D::new(cartesian[0], cartesian[1], cartesian[2]))
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_coordinate_conversions() {
        // Cartesian to spherical
        let cart = [1.0, 1.0, 1.0];
        let spherical = convert(
            &cart,
            CoordinateSystem::Cartesian,
            CoordinateSystem::Spherical,
        )
        .unwrap();

        assert_relative_eq!(spherical[0], 1.732051, epsilon = 1e-6); // r
        assert_relative_eq!(spherical[1], 54.735610, epsilon = 1e-6); // theta (degrees)
        assert_relative_eq!(spherical[2], 45.0, epsilon = 1e-6); // phi (degrees)

        // Spherical back to Cartesian
        let cart_again = convert(
            &spherical,
            CoordinateSystem::Spherical,
            CoordinateSystem::Cartesian,
        )
        .unwrap();

        assert_relative_eq!(cart_again[0], cart[0], epsilon = 1e-6);
        assert_relative_eq!(cart_again[1], cart[1], epsilon = 1e-6);
        assert_relative_eq!(cart_again[2], cart[2], epsilon = 1e-6);
    }
}
