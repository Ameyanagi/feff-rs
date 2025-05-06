/*
MIT License with FEFF10 Attribution

Copyright (c) 2025 Ameyanagi

Based on or developed using Distribution: FEFF10.0
Copyright (c) 2020 FEFF Project, University of Washington and SLAC National Accelerator Laboratory.
All rights reserved.
*/

//! Physical constants used in FEFF calculations

/// Speed of light in atomic units
pub const SPEED_OF_LIGHT: f64 = 137.036;

/// Bohr radius in Angstroms
pub const BOHR_RADIUS: f64 = 0.529177;

/// Conversion from Bohr to Angstroms
pub const BOHR_TO_ANGSTROM: f64 = BOHR_RADIUS;

/// Conversion from Angstroms to Bohr
pub const ANGSTROM_TO_BOHR: f64 = 1.0 / BOHR_RADIUS;

/// Rydberg energy in eV
pub const RYDBERG: f64 = 13.6057;

/// Conversion from eV to Hartree
pub const EV_TO_HARTREE: f64 = 1.0 / (2.0 * RYDBERG);

/// Conversion from Hartree to eV
pub const HARTREE_TO_EV: f64 = 2.0 * RYDBERG;
