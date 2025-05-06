/*
MIT License with FEFF10 Attribution

Copyright (c) 2025 Ameyanagi

Based on or developed using Distribution: FEFF10.0
Copyright (c) 2020 FEFF Project, University of Washington and SLAC National Accelerator Laboratory.
All rights reserved.
*/

//! Potential calculation module
//!
//! This module handles calculation of atomic potentials for FEFF.
//! The implementation follows the muffin-tin potential approach used in FEFF10,
//! where the space is divided into non-overlapping spherical regions centered on atoms.

pub mod electron_config;
mod errors;
mod exchange_correlation;
mod muffin_tin;

pub use errors::{PotentialError, Result};
pub use exchange_correlation::ExchangeCorrelationType;
pub use muffin_tin::{GridType, MuffinTinPotential, MuffinTinPotentialResult};

#[cfg(test)]
mod tests {
    use super::*;
    use crate::atoms::{Atom, AtomicStructure, PotentialType, Vector3D};

    #[test]
    fn test_muffin_tin_creation() {
        // Create a simple atomic structure with an iron atom
        let mut structure = AtomicStructure::new();
        let pot_fe = PotentialType::new(0, 26).unwrap(); // Iron
        structure.add_potential_type(pot_fe);

        // Add central iron atom
        let fe_idx = structure.add_atom(Atom::new(26, Vector3D::new(0.0, 0.0, 0.0), 0).unwrap());
        structure.set_central_atom(fe_idx).unwrap();

        // Calculate muffin-tin radius
        structure.calculate_muffin_tin_radii().unwrap();

        // Create muffin-tin potential calculator
        let mt_calculator = MuffinTinPotential::new(&structure);

        // Should be able to create the calculator
        assert!(mt_calculator.is_ok());
    }
}
