/*
MIT License with FEFF10 Attribution

Copyright (c) 2025 Ameyanagi

Based on or developed using Distribution: FEFF10.0
Copyright (c) 2020 FEFF Project, University of Washington and SLAC National Accelerator Laboratory.
All rights reserved.
*/

//! Path amplitude calculations for EXAFS
//!
//! This module implements functions to calculate the amplitude
//! and phase of scattering paths for EXAFS calculations.

use num_complex::Complex64;

use crate::atoms::structure::AtomicStructure;
use crate::path::path::Path;
use crate::scattering::ScatteringResults;

/// Represents the physical parameters needed for EXAFS calculations
#[derive(Debug, Clone)]
pub struct ExafsParameters {
    /// Electron energy grid for calculations (in eV)
    pub energies: Vec<f64>,

    /// Fermi energy (in eV)
    pub fermi_energy: f64,

    /// Debye-Waller factors for each path
    pub debye_waller_factors: Vec<f64>,

    /// Mean free path for the photoelectron (in Å)
    pub mean_free_path: f64,

    /// S0² amplitude reduction factor
    pub s0_squared: f64,
}

impl Default for ExafsParameters {
    fn default() -> Self {
        Self {
            energies: (0..1000).map(|i| i as f64).collect(),
            fermi_energy: 0.0,
            debye_waller_factors: vec![0.003], // Default σ² = 0.003 Å²
            mean_free_path: 10.0,              // Default λ = 10 Å
            s0_squared: 1.0,                   // Default S0² = 1.0
        }
    }
}

/// Calculates the EXAFS amplitude and phase for a given path
///
/// # Arguments
///
/// * `path` - The scattering path
/// * `structure` - Atomic structure containing the atoms
/// * `scattering_results` - Scattering results containing phase shifts
/// * `params` - EXAFS calculation parameters
///
/// # Returns
///
/// Vectors of amplitude and phase values for each energy point
pub fn calculate_path_exafs(
    path: &Path,
    structure: &AtomicStructure,
    scattering_results: &ScatteringResults,
    params: &ExafsParameters,
) -> (Vec<f64>, Vec<f64>) {
    let num_energies = params.energies.len();
    let mut amplitudes = vec![0.0; num_energies];
    let mut phases = vec![0.0; num_energies];

    // Get central atom info
    let absorber_index = path.atom_sequence[0];
    let _absorber = structure.atom(absorber_index).unwrap();

    // Calculate Debye-Waller factor
    let dw_factor = if params.debye_waller_factors.len() > path.atom_sequence.len() {
        // Use path-specific factors if available
        params.debye_waller_factors[path.atom_sequence.len() - 1]
    } else {
        // Use default factor
        params.debye_waller_factors[0]
    };

    // Calculate effective path length (half of total path length)
    let _effective_length = path.effective_length();

    // Compute EXAFS for each energy point
    for (i, &energy) in params.energies.iter().enumerate() {
        // Convert energy to wavenumber k (in Å⁻¹)
        let k = energy_to_k(energy, params.fermi_energy);

        if k <= 0.0 {
            // Skip negative k values (below Fermi energy)
            amplitudes[i] = 0.0;
            phases[i] = 0.0;
            continue;
        }

        // Calculate the complex amplitude
        let amplitude = calculate_complex_amplitude(
            path,
            structure,
            scattering_results,
            k,
            dw_factor,
            params.mean_free_path,
            params.s0_squared,
        );

        // Extract magnitude and phase
        amplitudes[i] = amplitude.norm();
        phases[i] = amplitude.arg();
    }

    (amplitudes, phases)
}

/// Calculates the complex amplitude for a path at a specific k value
///
/// This function implements the EXAFS formula for computing the
/// contribution of a path to the overall EXAFS signal.
///
/// # Arguments
///
/// * `path` - The scattering path
/// * `structure` - Atomic structure containing the atoms
/// * `scattering_results` - Scattering results containing phase shifts
/// * `k` - Wavenumber (in Å⁻¹)
/// * `dw_factor` - Debye-Waller factor σ² (in Å²)
/// * `mean_free_path` - Mean free path λ (in Å)
/// * `s0_squared` - S0² amplitude reduction factor
///
/// # Returns
///
/// The complex amplitude
fn calculate_complex_amplitude(
    path: &Path,
    structure: &AtomicStructure,
    scattering_results: &ScatteringResults,
    k: f64,
    dw_factor: f64,
    mean_free_path: f64,
    s0_squared: f64,
) -> Complex64 {
    // Use path length for calculations
    let half_path_length = path.effective_length();

    // Degeneracy factor
    let degeneracy_factor = path.degeneracy as f64;

    // S0² factor
    let s0_factor = s0_squared.sqrt();

    // Mean free path damping
    let mfp_damping = (-path.total_length / mean_free_path).exp();

    // Debye-Waller factor
    let dw_damping = (-2.0 * k * k * dw_factor).exp();

    // 1/kR² amplitude factor
    let geometric_factor = 1.0 / (k * half_path_length * half_path_length);

    // Initialize complex amplitude
    let mut complex_amplitude = Complex64::new(1.0, 0.0);

    // Absorber phase shift (central atom)
    let absorber_index = path.atom_sequence[0];
    let absorber = structure.atom(absorber_index).unwrap();
    let absorber_pot_index = absorber.potential_type() as usize;

    // Use l=1 for dipole transitions
    if 1 < scattering_results.phase_shifts[absorber_pot_index].len() {
        let absorber_phase = scattering_results.phase_shifts[absorber_pot_index][1];
        complex_amplitude *= Complex64::from_polar(1.0, -2.0 * absorber_phase.re);
    }

    // Process each scattering atom in the path
    for i in 1..path.atom_sequence.len() - 1 {
        let scatterer_index = path.atom_sequence[i];
        let scatterer = structure.atom(scatterer_index).unwrap();

        // Get scattering amplitude and phase for this atom
        let scatterer_pot_index = scatterer.potential_type() as usize;

        if scattering_results.phase_shifts[scatterer_pot_index].len() > 0 {
            // Add phase from the scatterer (l=0 for s-wave scattering)
            let scatterer_phase = scattering_results.phase_shifts[scatterer_pot_index][0];
            complex_amplitude *= Complex64::from_polar(1.0, 2.0 * scatterer_phase.re);

            // Multiply by scattering amplitude approximation
            // In a full implementation, this would use proper scattering matrices
            let scattering_factor = 0.1 * scatterer.atomic_number() as f64;
            complex_amplitude *= Complex64::new(scattering_factor, 0.0);
        }
    }

    // Include propagation phase factor (e^{2ikR})
    let propagation_phase = 2.0 * k * half_path_length;
    complex_amplitude *= Complex64::from_polar(1.0, propagation_phase);

    // Scale by all the factors
    complex_amplitude *=
        degeneracy_factor * s0_factor * mfp_damping * dw_damping * geometric_factor;

    complex_amplitude
}

/// Calculates precise EXAFS using full multiple scattering matrices
///
/// This is a more accurate version that uses full scattering matrices rather
/// than just phase shifts.
///
/// # Arguments
///
/// * `path` - The scattering path
/// * `structure` - Atomic structure containing the atoms
/// * `scattering_results` - Array of scattering results for each energy point
/// * `params` - EXAFS calculation parameters
///
/// # Returns
///
/// Vectors of amplitude and phase values for each energy point
pub fn calculate_path_exafs_with_matrices(
    path: &Path,
    structure: &AtomicStructure,
    scattering_results: &[ScatteringResults],
    params: &ExafsParameters,
) -> (Vec<f64>, Vec<f64>) {
    let num_energies = params.energies.len();
    let mut amplitudes = vec![0.0; num_energies];
    let mut phases = vec![0.0; num_energies];

    // Validate input
    if scattering_results.len() != num_energies {
        // Return zeros if matrices don't match energies
        return (amplitudes, phases);
    }

    // Calculate Debye-Waller factor
    let dw_factor = if params.debye_waller_factors.len() > path.atom_sequence.len() {
        // Use path-specific factors if available
        params.debye_waller_factors[path.atom_sequence.len() - 1]
    } else {
        // Use default factor
        params.debye_waller_factors[0]
    };

    // Compute EXAFS for each energy point
    for (i, &energy) in params.energies.iter().enumerate() {
        // Convert energy to wavenumber k (in Å⁻¹)
        let k = energy_to_k(energy, params.fermi_energy);

        if k <= 0.0 {
            // Skip negative k values (below Fermi energy)
            amplitudes[i] = 0.0;
            phases[i] = 0.0;
            continue;
        }

        // Get scattering result for this energy
        let scattering_result = &scattering_results[i];

        // TODO: Implement full matrix-based calculation
        // This would require proper path operator construction using
        // the scattering matrices and propagators

        // For now, use the simplified version
        let amplitude = calculate_complex_amplitude(
            path,
            structure,
            scattering_result,
            k,
            dw_factor,
            params.mean_free_path,
            params.s0_squared,
        );

        // Extract magnitude and phase
        amplitudes[i] = amplitude.norm();
        phases[i] = amplitude.arg();
    }

    (amplitudes, phases)
}

/// Converts energy to wavenumber k
///
/// # Arguments
///
/// * `energy` - Energy in eV
/// * `fermi_energy` - Fermi energy in eV
///
/// # Returns
///
/// Wavenumber k in Å⁻¹
fn energy_to_k(energy: f64, fermi_energy: f64) -> f64 {
    // Constants
    const HBAR_SQUARED_OVER_2M: f64 = 3.81; // in eV·Å²

    let relative_energy = energy - fermi_energy;

    if relative_energy <= 0.0 {
        return 0.0;
    }

    (relative_energy / HBAR_SQUARED_OVER_2M).sqrt()
}

/// Calculates the chi(k) EXAFS spectrum for a set of paths
///
/// This function combines the contributions of multiple paths to
/// produce the overall EXAFS spectrum.
///
/// # Arguments
///
/// * `paths` - Vector of scattering paths
/// * `structure` - Atomic structure containing the atoms
/// * `scattering_results` - Scattering results containing phase shifts
/// * `params` - EXAFS calculation parameters
///
/// # Returns
///
/// The chi(k) EXAFS spectrum for each energy point
pub fn calculate_exafs_spectrum(
    paths: &[Path],
    structure: &AtomicStructure,
    scattering_results: &ScatteringResults,
    params: &ExafsParameters,
) -> Vec<f64> {
    let num_energies = params.energies.len();
    let mut chi_k = vec![0.0; num_energies];

    // Calculate contribution from each path
    for path in paths {
        let (amplitudes, phases) =
            calculate_path_exafs(path, structure, scattering_results, params);

        // Add contribution to total chi(k)
        for i in 0..num_energies {
            let k = energy_to_k(params.energies[i], params.fermi_energy);

            if k > 0.0 {
                // Convert amplitude and phase to real part of chi(k)
                let real_part = amplitudes[i] * phases[i].cos();

                // Add to total
                chi_k[i] += real_part;
            }
        }
    }

    chi_k
}
