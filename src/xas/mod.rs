/*
MIT License with FEFF10 Attribution

Copyright (c) 2025 Ameyanagi

Based on or developed using Distribution: FEFF10.0
Copyright (c) 2020 FEFF Project, University of Washington and SLAC National Accelerator Laboratory.
All rights reserved.
*/

//! X-ray Absorption Spectroscopy (XAS) module
//!
//! This module handles calculation of XAS spectra for FEFF, including XANES and EXAFS.
//! It provides functionality for computing absorption spectra with different theoretical
//! approaches like multiple scattering and self-consistent field calculations.

mod core_hole;

use crate::atoms::Result as AtomResult;
use crate::scattering::ScatteringMatrixResults;

pub use core_hole::{calculate_with_core_hole, CoreHoleConfig, CoreHoleMethod};

/// Calculate XAS spectrum using multiple scattering theory
///
/// # Arguments
///
/// * `scattering` - Scattering matrix results including phase shifts and path operators
/// * `polarization` - Optional polarization vector (for polarized XANES)
///
/// # Returns
///
/// Vector of (energy, mu) pairs representing the XAS spectrum
pub fn calculate_xas_spectrum(
    scattering: &ScatteringMatrixResults,
    _polarization: Option<[f64; 3]>,
) -> AtomResult<Vec<(f64, f64)>> {
    // This is a placeholder that will be implemented later
    // For now, just return a simple approximation based on the
    // imaginary part of the absorbing atom's phase shift

    let energy = scattering.energy;
    let phase_shifts = &scattering.phase_shifts;

    if phase_shifts.is_empty() {
        return Ok(vec![(energy, 0.0)]);
    }

    // Use the phase shift of the absorbing atom (first potential)
    // and the first angular momentum (usually s-wave for K-edge)
    let absorbing_atom_phase = phase_shifts[0][0];

    // The absorption coefficient is related to the imaginary part of the phase shift
    let mu = absorbing_atom_phase.im;

    Ok(vec![(energy, mu)])
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::atoms::{Atom, AtomicStructure, PotentialType, Vector3D};
    use crate::scattering::{
        calculate_phase_shifts_with_method,
        calculate_scattering_matrices_old as calculate_scattering_matrices, PhaseShiftMethod,
    };

    #[test]
    fn test_placeholder() {
        assert!(true);
    }

    #[test]
    fn test_xas_calculation() {
        // Create a simple iron atom structure
        let fe_potential = PotentialType::new(0, 26).unwrap();
        let fe_atom = Atom::new(26, Vector3D::new(0.0, 0.0, 0.0), 0).unwrap();

        let mut structure = AtomicStructure::new();
        structure.add_potential_type(fe_potential);
        let fe_idx = structure.add_atom(fe_atom);
        structure.set_central_atom(fe_idx).unwrap();

        // Calculate muffin-tin radius
        structure.calculate_muffin_tin_radii().unwrap();

        // Calculate phase shifts
        let energy = 100.0; // eV
        let max_l = 3;

        let _phase_shifts = calculate_phase_shifts_with_method(
            &structure,
            energy,
            max_l,
            PhaseShiftMethod::Approximate,
        )
        .unwrap();

        // Calculate scattering matrices
        let scattering = calculate_scattering_matrices(&structure, energy, max_l).unwrap();

        // Calculate XAS spectrum
        let spectrum = calculate_xas_spectrum(&scattering, None).unwrap();

        // Should return at least one point
        assert!(!spectrum.is_empty());

        // Energy should match input
        assert_eq!(spectrum[0].0, energy);

        // Absorption coefficient should be positive
        assert!(spectrum[0].1 > 0.0);
    }
}
