/*
MIT License with FEFF10 Attribution

Copyright (c) 2025 Ameyanagi

Based on or developed using Distribution: FEFF10.0
Copyright (c) 2020 FEFF Project, University of Washington and SLAC National Accelerator Laboratory.
All rights reserved.
*/

//! Utility functions for FEFF calculations
//!
//! This module provides common utilities used throughout the FEFF code.

use num_complex::Complex64;

/// Physical constants used in FEFF calculations
pub mod constants {
    /// Speed of light in atomic units
    pub const SPEED_OF_LIGHT: f64 = 137.036;
    
    /// Bohr radius in Angstroms
    pub const BOHR_RADIUS: f64 = 0.529177;
    
    /// Rydberg energy in eV
    pub const RYDBERG: f64 = 13.6057;
    
    /// Conversion from eV to Hartree
    pub const EV_TO_HARTREE: f64 = 1.0 / (2.0 * RYDBERG);
    
    /// Conversion from Hartree to eV
    pub const HARTREE_TO_EV: f64 = 2.0 * RYDBERG;
}

/// Convert from Angstroms to Bohr radii
pub fn angstrom_to_bohr(angstrom: f64) -> f64 {
    angstrom / constants::BOHR_RADIUS
}

/// Convert from Bohr radii to Angstroms
pub fn bohr_to_angstrom(bohr: f64) -> f64 {
    bohr * constants::BOHR_RADIUS
}

/// Convert energy from eV to Hartree
pub fn ev_to_hartree(ev: f64) -> f64 {
    ev * constants::EV_TO_HARTREE
}

/// Convert energy from Hartree to eV
pub fn hartree_to_ev(hartree: f64) -> f64 {
    hartree * constants::HARTREE_TO_EV
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    
    #[test]
    fn test_unit_conversions() {
        let angstrom = 1.0;
        let bohr = angstrom_to_bohr(angstrom);
        assert_relative_eq!(bohr_to_angstrom(bohr), angstrom, epsilon = 1e-10);
        
        let ev = 10.0;
        let hartree = ev_to_hartree(ev);
        assert_relative_eq!(hartree_to_ev(hartree), ev, epsilon = 1e-10);
    }
}
