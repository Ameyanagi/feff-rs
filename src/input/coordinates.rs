/*
MIT License with FEFF10 Attribution

Copyright (c) 2025 Ameyanagi

Based on or developed using Distribution: FEFF10.0
Copyright (c) 2020 FEFF Project, University of Washington and SLAC National Accelerator Laboratory.
All rights reserved.
*/

//! Coordinate system definition for FEFF input

/// Coordinate system enumeration for atomic positions
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum CoordinateSystem {
    /// Cartesian (x, y, z)
    Cartesian,
    /// Spherical (r, theta, phi)
    Spherical,
    /// Cylindrical (rho, phi, z)
    Cylindrical,
}
