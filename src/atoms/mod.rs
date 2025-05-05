/*
MIT License with FEFF10 Attribution

Copyright (c) 2025 Ameyanagi

Based on or developed using Distribution: FEFF10.0
Copyright (c) 2020 FEFF Project, University of Washington and SLAC National Accelerator Laboratory.
All rights reserved.
*/

//! Atomic data and calculations module
//!
//! This module provides atomic data and calculations for FEFF.
//! The implementation follows FEFF10's approach for representing
//! atoms, potential types, and atomic structures used in calculations.

pub mod atom;
pub mod coordinates;
pub mod database;
pub mod errors;
pub mod potential;
pub mod structure;
pub mod vector;

// Re-export commonly used types and functions for convenience
pub use atom::Atom;
pub use coordinates::CoordinateSystem;
pub use errors::{AtomError, Result};
pub use potential::PotentialType;
pub use structure::AtomicStructure;
pub use vector::Vector3D;
