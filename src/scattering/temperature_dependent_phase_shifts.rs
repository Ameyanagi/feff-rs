/*
MIT License with FEFF10 Attribution

Copyright (c) 2025 Ameyanagi

Based on or developed using Distribution: FEFF10.0
Copyright (c) 2020 FEFF Project, University of Washington and SLAC National Accelerator Laboratory.
All rights reserved.
*/

//! Temperature-dependent phase shifts for XAS calculations
//!
//! This module extends the phase shift calculations to include temperature effects.
//! Temperature impacts phase shifts through several mechanisms:
//!
//! 1. **Thermal expansion**: Lattice parameters increase with temperature, changing
//!    the muffin-tin radius and effective potentials.
//!
//! 2. **Vibrational effects**: Atomic vibrations modify the effective potential
//!    experienced by photoelectrons, particularly at high temperatures.
//!
//! 3. **Anharmonic effects**: Beyond the harmonic approximation, temperature-dependent
//!    anharmonicity can significantly modify phase shifts at high temperatures.
//!
//! Temperature-dependent phase shifts are essential for accurate EXAFS analysis,
//! especially when comparing data collected at different temperatures or when
//! analyzing thermal behavior of materials.

use super::phase_shifts::calculate_phase_shifts;
use crate::atoms::{AtomicStructure, Result as AtomResult};
use crate::utils::thermal::ThermalModel;
use crate::xas::thermal::ThermalParameters;
use num_complex::Complex64;

/// Temperature-dependent phase shift correction factors
///
/// These values represent typical corrections to phase shifts due to temperature effects.
/// Different angular momentum channels (s, p, d, f) respond differently to temperature.
const TEMP_CORRECTION_FACTORS: [f64; 4] = [0.020, 0.015, 0.010, 0.005];

/// Calculate temperature-dependent phase shifts for an atomic structure
///
/// This function extends the standard phase shift calculation by including
/// thermal effects on phase shifts.
///
/// # Arguments
///
/// * `structure` - The atomic structure containing atoms and potentials
/// * `energy` - Energy in eV
/// * `max_l` - Maximum angular momentum to include in calculations
/// * `temperature` - Temperature in Kelvin
///
/// # Returns
///
/// A ScatteringResults object containing temperature-corrected phase shifts and T-matrices
pub fn calculate_temperature_dependent_phase_shifts(
    structure: &AtomicStructure,
    energy: f64,
    max_l: i32,
    temperature: f64,
) -> AtomResult<super::ScatteringResults> {
    // Calculate base phase shifts first
    let mut scattering_results = calculate_phase_shifts(structure, energy, max_l)?;

    // No correction needed at absolute zero
    if temperature <= 0.0 {
        return Ok(scattering_results);
    }

    // Apply temperature corrections to phase shifts
    for potential_idx in 0..scattering_results.phase_shifts.len() {
        let shifts = &mut scattering_results.phase_shifts[potential_idx];

        // Get potential info
        let potential = structure
            .potential_type(potential_idx as usize)
            .ok_or_else(|| {
                crate::atoms::errors::AtomError::CalculationError(format!(
                    "Invalid potential index: {}",
                    potential_idx
                ))
            })?;

        // Calculate thermal correction based on Debye temperature
        let debye_temperature = potential.debye_temperature().unwrap_or(300.0);

        // Apply temperature corrections to each angular momentum channel
        for l in 0..=max_l as usize {
            if l < shifts.len() {
                // Get the appropriate correction factor for this l value (limit to first 4 values)
                let corr_factor = if l < TEMP_CORRECTION_FACTORS.len() {
                    TEMP_CORRECTION_FACTORS[l]
                } else {
                    TEMP_CORRECTION_FACTORS[TEMP_CORRECTION_FACTORS.len() - 1] / (l as f64 - 2.0)
                };

                // Basic temperature scaling based on the Debye model
                // (T/θD) term represents thermal activation relative to Debye temperature
                let temp_ratio = temperature / debye_temperature;
                let temp_factor = if temp_ratio < 1.0 {
                    // Below Debye temperature, less than linear scaling
                    temp_ratio * temp_ratio.sqrt()
                } else {
                    // Above Debye temperature, approximately linear scaling with diminishing returns
                    temp_ratio.sqrt() * (1.0 + 0.5 * (temp_ratio - 1.0).min(2.0))
                };

                // Photoelectron energy also affects the temperature dependence
                // Higher energy photoelectrons are less affected by thermal disorder
                let energy_factor = 100.0 / (100.0 + energy);

                // Calculate the phase shift correction
                let phase_correction = corr_factor * temp_factor * energy_factor;

                // Modify both real and imaginary parts of the phase shift
                // Real part: primarily affected by changes in potential shape
                shifts[l] = Complex64::new(
                    shifts[l].re + phase_correction,
                    // Imaginary part: increased due to inelastic scattering from vibrating atoms
                    shifts[l].im * (1.0 + 0.5 * phase_correction),
                );
            }
        }
    }

    // Recalculate T-matrices with the updated phase shifts
    for potential_idx in 0..scattering_results.phase_shifts.len() {
        let shifts = &scattering_results.phase_shifts[potential_idx];
        let t_matrix = crate::utils::matrix::compute_t_matrix(shifts, max_l).map_err(|e| {
            crate::atoms::errors::AtomError::CalculationError(format!(
                "Failed to compute T-matrix: {}",
                e
            ))
        })?;

        scattering_results.t_matrices[potential_idx] = t_matrix;
    }

    // Set the temperature field
    scattering_results.temperature = Some(temperature);

    Ok(scattering_results)
}

/// Calculate path-specific temperature-dependent phase shifts
///
/// This function calculates phase shift corrections for specific scattering paths,
/// taking into account path-specific thermal parameters.
///
/// # Arguments
///
/// * `phase_shift` - Original phase shift (without temperature effects)
/// * `energy` - Energy in eV
/// * `l` - Angular momentum quantum number
/// * `thermal_model` - Reference to a thermal model for calculating MSD
/// * `temperature` - Temperature in Kelvin
/// * `path_distance` - Distance of the scattering path in Å
///
/// # Returns
///
/// The temperature-corrected phase shift
pub fn calculate_path_thermal_phase_shift(
    phase_shift: Complex64,
    energy: f64,
    l: i32,
    thermal_model: &dyn ThermalModel,
    temperature: f64,
    path_distance: f64,
) -> Complex64 {
    // No correction at absolute zero
    if temperature <= 0.0 {
        return phase_shift;
    }

    // Calculate thermal mean-square displacement
    let msd = thermal_model.mean_square_displacement(temperature);

    // Convert energy to wavenumber (k ≈ sqrt(E/3.81))
    let k = (energy / 3.81).sqrt();

    // Basic thermal phase shift correction based on MSD and wavenumber
    // This follows the FEFF approach with empirical correction factors
    let l_index = l.min(3) as usize;
    let correction_factor = TEMP_CORRECTION_FACTORS[l_index];

    // Calculate anharmonic term (more pronounced at high temperatures)
    // This uses the cubic anharmonicity approximation
    let anharmonic_term = if temperature > 300.0 {
        let reduced_temp = temperature / 300.0;
        0.001 * (reduced_temp - 1.0).powi(2) * path_distance
    } else {
        0.0
    };

    // Calculate thermal phase correction
    // The phase correction includes:
    // 1. Linear term proportional to k²σ² (Debye-Waller factor related)
    // 2. Path distance scaling (longer paths have more significant corrections)
    // 3. Anharmonic term for high temperatures
    let phase_correction =
        correction_factor * k * k * msd * (path_distance / 2.5) + anharmonic_term;

    // Energy-dependent damping of the correction (less effect at high energy)
    let energy_damping = 100.0 / (100.0 + energy);
    let effective_correction = phase_correction * energy_damping;

    // Apply corrections to both real and imaginary parts
    Complex64::new(
        phase_shift.re + effective_correction,
        // Imaginary part increases more significantly due to inelastic effects
        phase_shift.im * (1.0 + 0.5 * effective_correction),
    )
}

/// Calculate thermal phase shifts for an entire scattering path
///
/// This function applies thermal corrections to all phase shifts relevant
/// for a specific scattering path.
///
/// # Arguments
///
/// * `structure` - The atomic structure containing atoms and potentials
/// * `path_indices` - Indices of atoms in the scattering path
/// * `phase_shifts` - Original phase shifts for each potential type
/// * `energy` - Energy in eV
/// * `max_l` - Maximum angular momentum
/// * `thermal_params` - Thermal parameters for the calculation
///
/// # Returns
///
/// Vector of temperature-corrected phase shifts for each potential type in the path
pub fn apply_thermal_corrections_to_path(
    structure: &AtomicStructure,
    path_indices: &[usize],
    phase_shifts: &[Vec<Complex64>],
    energy: f64,
    max_l: i32,
    thermal_params: &ThermalParameters,
) -> Vec<Vec<Complex64>> {
    // If path is empty or only contains one atom, return original shifts
    if path_indices.len() < 2 {
        return phase_shifts.to_vec();
    }

    // Create a copy of the phase shifts to modify
    let mut corrected_shifts = phase_shifts.to_vec();

    // For each atom in the path, apply thermal corrections
    for i in 0..path_indices.len() {
        // Get atom and its potential type
        if let Some(atom) = structure.atom(path_indices[i]) {
            let pot_type = atom.potential_type();

            // Skip if potential type is invalid
            if pot_type as usize >= corrected_shifts.len() {
                continue;
            }

            // Calculate path distance (use average of connected legs)
            let mut path_distance = 0.0;
            let mut count = 0;

            // Consider legs connecting to this atom
            if i > 0 {
                if let (Some(prev_atom), Some(curr_atom)) = (
                    structure.atom(path_indices[i - 1]),
                    structure.atom(path_indices[i]),
                ) {
                    let prev_pos = prev_atom.position();
                    let curr_pos = curr_atom.position();
                    path_distance += (curr_pos.clone() - prev_pos.clone()).length();
                    count += 1;
                }
            }

            if i < path_indices.len() - 1 {
                if let (Some(curr_atom), Some(next_atom)) = (
                    structure.atom(path_indices[i]),
                    structure.atom(path_indices[i + 1]),
                ) {
                    let curr_pos = curr_atom.position();
                    let next_pos = next_atom.position();
                    path_distance += (next_pos.clone() - curr_pos.clone()).length();
                    count += 1;
                }
            }

            // Calculate average leg length
            if count > 0 {
                path_distance /= count as f64;
            } else {
                // Default to a typical first-shell distance if we can't calculate
                path_distance = 2.5;
            }

            // Create thermal model for this atom
            let atom_z = atom.atomic_number();

            // Check for pair-specific parameters
            let reduced_mass = if i > 0 {
                if let Some(prev_atom) = structure.atom(path_indices[i - 1]) {
                    // Calculate reduced mass if we have adjacent atoms
                    let prev_z = prev_atom.atomic_number();
                    calculate_reduced_mass(prev_z as f64, atom_z as f64)
                } else {
                    // Default to atom's mass if we can't calculate reduced mass
                    atom_z as f64
                }
            } else {
                atom_z as f64
            };

            // Create thermal model
            let thermal_model = thermal_params.create_model(reduced_mass, Some(path_distance));

            // Apply thermal corrections to each angular momentum channel
            for l in 0..=max_l as usize {
                if l < corrected_shifts[pot_type as usize].len() {
                    let original_shift = corrected_shifts[pot_type as usize][l];
                    corrected_shifts[pot_type as usize][l] = calculate_path_thermal_phase_shift(
                        original_shift,
                        energy,
                        l as i32,
                        thermal_model.as_ref(),
                        thermal_params.temperature,
                        path_distance,
                    );
                }
            }
        }
    }

    corrected_shifts
}

/// Calculate reduced mass for a pair of atoms
///
/// The reduced mass is used in thermal vibration calculations and is defined as:
/// μ = (m₁*m₂)/(m₁+m₂)
///
/// # Arguments
///
/// * `mass1` - Mass of first atom in atomic mass units
/// * `mass2` - Mass of second atom in atomic mass units
///
/// # Returns
///
/// The reduced mass in atomic mass units
fn calculate_reduced_mass(mass1: f64, mass2: f64) -> f64 {
    (mass1 * mass2) / (mass1 + mass2)
}

#[cfg(test)]
mod tests {
    use super::*;
    // These are used in other tests in this module
    #[allow(unused_imports)]
    use crate::atoms::{Atom, PotentialType, Vector3D};
    use crate::utils::thermal::{CorrelatedDebyeModel, DebyeModel};

    #[test]
    fn test_temperature_corrections() {
        // Create a simple thermal model
        let debye_model = DebyeModel::new(300.0, 55.845); // Iron

        // Test phase shift at two temperatures
        let base_phase = Complex64::new(0.3, 0.1);
        let energy = 100.0;
        let l = 1;
        let path_distance = 2.5;

        // Room temperature
        let room_temp_phase = calculate_path_thermal_phase_shift(
            base_phase,
            energy,
            l,
            &debye_model,
            300.0,
            path_distance,
        );

        // High temperature
        let high_temp_phase = calculate_path_thermal_phase_shift(
            base_phase,
            energy,
            l,
            &debye_model,
            600.0,
            path_distance,
        );

        // Higher temperature should increase both real and imaginary parts
        assert!(high_temp_phase.re > room_temp_phase.re);
        assert!(high_temp_phase.im > room_temp_phase.im);

        // Both should be larger than the base phase
        assert!(room_temp_phase.re > base_phase.re);
        assert!(room_temp_phase.im > base_phase.im);
    }

    #[test]
    fn test_path_distance_effect() {
        // Create a thermal model
        let debye_model = DebyeModel::new(300.0, 55.845); // Iron

        // Test phase shift for different path distances
        let base_phase = Complex64::new(0.3, 0.1);
        let energy = 100.0;
        let l = 1;
        let temperature = 300.0;

        // Short path
        let short_path_phase = calculate_path_thermal_phase_shift(
            base_phase,
            energy,
            l,
            &debye_model,
            temperature,
            2.0,
        );

        // Long path
        let long_path_phase = calculate_path_thermal_phase_shift(
            base_phase,
            energy,
            l,
            &debye_model,
            temperature,
            4.0,
        );

        // Longer paths should have larger phase corrections
        assert!(long_path_phase.re > short_path_phase.re);
    }

    #[test]
    fn test_angular_momentum_dependence() {
        // Create a thermal model
        let debye_model = DebyeModel::new(300.0, 55.845); // Iron

        // Test phase shift for different angular momentum channels
        let base_phase = Complex64::new(0.3, 0.1);
        let energy = 100.0;
        let temperature = 300.0;
        let path_distance = 2.5;

        // s-orbital (l=0)
        let s_orbital_phase = calculate_path_thermal_phase_shift(
            base_phase,
            energy,
            0,
            &debye_model,
            temperature,
            path_distance,
        );

        // p-orbital (l=1)
        let p_orbital_phase = calculate_path_thermal_phase_shift(
            base_phase,
            energy,
            1,
            &debye_model,
            temperature,
            path_distance,
        );

        // d-orbital (l=2)
        let d_orbital_phase = calculate_path_thermal_phase_shift(
            base_phase,
            energy,
            2,
            &debye_model,
            temperature,
            path_distance,
        );

        // Lower angular momentum should have larger thermal effects
        assert!(s_orbital_phase.re > p_orbital_phase.re);
        assert!(p_orbital_phase.re > d_orbital_phase.re);
    }

    #[test]
    fn test_reduced_mass_calculation() {
        // Equal masses
        let equal_mass = calculate_reduced_mass(12.0, 12.0); // Carbon-Carbon
        assert!((equal_mass - 6.0).abs() < 1e-10); // μ = m/2 for equal masses

        // Different masses
        let fe_o_mass = calculate_reduced_mass(55.845, 16.0); // Iron-Oxygen
        let expected = (55.845 * 16.0) / (55.845 + 16.0);
        assert!((fe_o_mass - expected).abs() < 1e-10);
    }

    #[test]
    fn test_anharmonic_effects() {
        // Create a thermal model
        let debye_model = DebyeModel::new(300.0, 55.845); // Iron

        // Test phase shift at high temperature where anharmonic effects matter
        let base_phase = Complex64::new(0.3, 0.1);
        let energy = 100.0;
        let l = 1;
        let path_distance = 2.5;

        // Below anharmonic threshold
        let low_temp_phase = calculate_path_thermal_phase_shift(
            base_phase,
            energy,
            l,
            &debye_model,
            300.0,
            path_distance,
        );

        // Above anharmonic threshold
        let high_temp_phase = calculate_path_thermal_phase_shift(
            base_phase,
            energy,
            l,
            &debye_model,
            600.0,
            path_distance,
        );

        // Very high temperature with anharmonic effects
        let very_high_temp_phase = calculate_path_thermal_phase_shift(
            base_phase,
            energy,
            l,
            &debye_model,
            900.0,
            path_distance,
        );

        // Anharmonicity means phase shift doesn't scale linearly with temperature
        let low_to_mid_diff = high_temp_phase.re - low_temp_phase.re;
        let mid_to_high_diff = very_high_temp_phase.re - high_temp_phase.re;

        // The difference should be larger between mid and high due to anharmonicity
        assert!(mid_to_high_diff > low_to_mid_diff);
    }

    #[test]
    fn test_correlation_effects() {
        // Create a correlated thermal model
        let correlated_model = CorrelatedDebyeModel::new(300.0, 55.845, 2.5); // Iron, first shell

        // Create a standard thermal model
        let standard_model = DebyeModel::new(300.0, 55.845);

        // Test phase shift calculation
        let base_phase = Complex64::new(0.3, 0.1);
        let energy = 100.0;
        let l = 1;
        let temperature = 300.0;
        let path_distance = 2.5;

        // Calculate with correlation
        let correlated_phase = calculate_path_thermal_phase_shift(
            base_phase,
            energy,
            l,
            &correlated_model,
            temperature,
            path_distance,
        );

        // Calculate without correlation
        let uncorrelated_phase = calculate_path_thermal_phase_shift(
            base_phase,
            energy,
            l,
            &standard_model,
            temperature,
            path_distance,
        );

        // Correlation reduces thermal effects, so the correlated phase shift
        // should be closer to the original phase
        let correlated_diff = (correlated_phase.re - base_phase.re).abs();
        let uncorrelated_diff = (uncorrelated_phase.re - base_phase.re).abs();

        assert!(correlated_diff < uncorrelated_diff);
    }
}
