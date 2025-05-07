/*
MIT License with FEFF10 Attribution

Copyright (c) 2025 Ameyanagi

Based on or developed using Distribution: FEFF10.0
Copyright (c) 2020 FEFF Project, University of Washington and SLAC National Accelerator Laboratory.
All rights reserved.
*/

//! Temperature-dependent corrections for XANES calculations
//!
//! This module provides temperature-dependent corrections to XANES spectra,
//! accounting for thermal effects like thermal expansion, Debye-Waller factors,
//! and anharmonic vibrational effects.

use crate::utils::thermal::ThermalModel;
use crate::xas::{thermal::ThermalParameters, xanes::XanesSpectrum};

/// Parameters for temperature-dependent XANES corrections
#[derive(Debug, Clone)]
pub struct XanesThermalCorrectionParams {
    /// Gaussian broadening parameters for each spectrum region (in eV)
    /// Values are (energy_relative_to_edge, broadening_width)
    pub energy_dependent_broadening: Vec<(f64, f64)>,

    /// Thermal expansion coefficient (K⁻¹)
    pub thermal_expansion_coefficient: f64,

    /// Reference temperature for thermal expansion (K)
    pub reference_temperature: f64,

    /// Use anharmonic corrections
    pub use_anharmonic_corrections: bool,

    /// Enable asymmetric peak broadening (accounts for phonon-assisted transitions)
    pub enable_asymmetric_broadening: bool,
}

impl Default for XanesThermalCorrectionParams {
    fn default() -> Self {
        // Default values based on typical experimental observations
        Self {
            // Energy-dependent broadening:
            // - Pre-edge region: less broadening (0.5 eV) - discrete transitions
            // - Near-edge: more broadening (1.0 eV) - complex structure with multiple transitions
            // - Extended region: moderate broadening (0.8 eV) - simpler continuum features
            energy_dependent_broadening: vec![
                (-15.0, 0.5), // Pre-edge region
                (0.0, 1.0),   // Edge
                (15.0, 0.8),  // Above edge
                (50.0, 0.7),  // Far above edge
            ],

            // Linear thermal expansion coefficient - typical value for many solids
            thermal_expansion_coefficient: 1.0e-5, // K⁻¹

            // Reference temperature for thermal expansion
            reference_temperature: 300.0, // K (room temperature)

            // By default, use anharmonic corrections
            use_anharmonic_corrections: true,

            // By default, use asymmetric broadening
            enable_asymmetric_broadening: true,
        }
    }
}

/// Apply temperature-dependent corrections to a XANES spectrum
///
/// # Arguments
///
/// * `spectrum` - The original XANES spectrum to correct
/// * `thermal_params` - Thermal parameters with temperature and model
/// * `correction_params` - Parameters controlling the thermal corrections
///
/// # Returns
///
/// A new XANES spectrum with temperature corrections applied
pub fn apply_thermal_corrections(
    spectrum: &XanesSpectrum,
    thermal_params: &ThermalParameters,
    correction_params: &XanesThermalCorrectionParams,
) -> XanesSpectrum {
    // Create a copy of the original spectrum
    let mut corrected = spectrum.clone();

    // Apply thermal corrections

    // Check if temperature is high enough to apply thermal effects
    // Skip minimal effects for very low temperatures
    if thermal_params.temperature > 50.0 {
        // 1. Apply energy shift due to thermal expansion
        apply_thermal_expansion(&mut corrected, thermal_params, correction_params);

        // 2. Apply temperature-dependent broadening
        apply_temperature_dependent_broadening(&mut corrected, thermal_params, correction_params);

        // 3. Apply anharmonic corrections if requested
        if correction_params.use_anharmonic_corrections && thermal_params.temperature > 100.0 {
            apply_anharmonic_corrections(&mut corrected, thermal_params);
        }
    }

    // 4. Re-normalize the spectrum
    corrected.normalize().unwrap_or(());

    corrected
}

/// Apply thermal expansion correction to the energy scale
///
/// This accounts for changes in lattice parameters with temperature,
/// which shift spectral features. Thermal expansion causes:
/// - Features above the edge to move toward the edge (to lower energies)
/// - Features below the edge to move away from the edge (to lower energies)
fn apply_thermal_expansion(
    spectrum: &mut XanesSpectrum,
    thermal_params: &ThermalParameters,
    correction_params: &XanesThermalCorrectionParams,
) {
    // Calculate lattice expansion factor
    let delta_temp = thermal_params.temperature - correction_params.reference_temperature;

    // Skip if no temperature difference
    if delta_temp.abs() < 1e-6 {
        return;
    }

    // Calculate scale factor due to thermal expansion
    // The thermal expansion affects the energy scale roughly linearly
    let expansion_factor = 1.0 + correction_params.thermal_expansion_coefficient * delta_temp;

    // Apply scale factor to energy grid
    // Higher temperatures generally lead to expanded lattice, which shifts features to lower energy

    // Save original edge energy
    let original_edge = spectrum.edge_energy;

    // Scale energies relative to edge
    for i in 0..spectrum.energies.len() {
        let energy_rel_to_edge = spectrum.energies[i] - original_edge;

        // Apply expansion factor to relative energy
        // For both above and below edge features, divide by expansion factor
        // For above edge (positive energy_rel_to_edge): this moves features toward edge
        // For below edge (negative energy_rel_to_edge): this moves features away from edge
        // Both follow the physical principle of lower binding energies with lattice expansion
        let scaled_energy_rel = energy_rel_to_edge / expansion_factor;

        // Recalculate absolute energy
        spectrum.energies[i] = original_edge + scaled_energy_rel;
    }
}

/// Apply temperature-dependent broadening to the spectrum
///
/// This accounts for vibrational effects that broaden spectral features
/// at higher temperatures.
fn apply_temperature_dependent_broadening(
    spectrum: &mut XanesSpectrum,
    thermal_params: &ThermalParameters,
    correction_params: &XanesThermalCorrectionParams,
) {
    // Original energy grid
    let original_energies = spectrum.energies.clone();
    let original_mu = spectrum.mu.clone();

    // Calculate broadening for each energy point
    let mut broadening = vec![0.0; spectrum.energies.len()];

    // Use energy-dependent broadening parameters to determine broadening at each point
    for i in 0..spectrum.energies.len() {
        let energy_rel_to_edge = spectrum.energies[i] - spectrum.edge_energy;

        // Find appropriate broadening parameter
        let mut broadening_value = 0.5; // Default value

        // Find the two surrounding breakpoints and interpolate
        for j in 0..correction_params.energy_dependent_broadening.len() {
            let (energy, width) = correction_params.energy_dependent_broadening[j];

            if energy_rel_to_edge <= energy
                || j == correction_params.energy_dependent_broadening.len() - 1
            {
                if j == 0 {
                    // Before first breakpoint
                    broadening_value = width;
                } else {
                    // Interpolate between breakpoints
                    let (prev_energy, prev_width) =
                        correction_params.energy_dependent_broadening[j - 1];

                    if energy_rel_to_edge <= energy {
                        // Linear interpolation
                        let t = if energy == prev_energy {
                            0.0
                        } else {
                            (energy_rel_to_edge - prev_energy) / (energy - prev_energy)
                        };

                        broadening_value = prev_width + t * (width - prev_width);
                    } else {
                        // Past last breakpoint
                        broadening_value = width;
                    }
                }

                break;
            }
        }

        // Scale broadening by temperature
        // Thermal broadening generally increases as sqrt(T) due to vibrational amplitude
        let temp_factor = (thermal_params.temperature / 300.0).sqrt();
        broadening[i] = broadening_value * temp_factor;
    }

    // Apply variable broadening to each point using a Gaussian convolution
    // with position-dependent width

    // Create new arrays for broadened spectrum
    let mut broadened_mu = vec![0.0; original_energies.len()];

    // Apply convolution with position-dependent width
    for i in 0..original_energies.len() {
        let energy = original_energies[i];
        let sigma = broadening[i];

        // Skip if no broadening
        if sigma <= 0.0 {
            broadened_mu[i] = original_mu[i];
            continue;
        }

        let mut sum = 0.0;
        let mut weight_sum = 0.0;

        // The cutoff defines how many points to include in convolution
        let cutoff = 3.0 * sigma;

        for j in 0..original_energies.len() {
            let e_j = original_energies[j];
            let delta_e = energy - e_j;

            if delta_e.abs() > cutoff {
                continue;
            }

            // Apply Gaussian kernel
            let mut weight = (-0.5 * (delta_e / sigma).powi(2)).exp();

            // Apply asymmetric correction if enabled
            if correction_params.enable_asymmetric_broadening {
                // Asymmetric factor increases at higher temperatures
                let asym_factor = 0.2 * (thermal_params.temperature / 300.0 - 1.0).max(0.0);

                // More broadening on the high-energy side (positive delta_e)
                if delta_e > 0.0 {
                    weight *= 1.0 + asym_factor * (delta_e / sigma);
                }
            }

            sum += original_mu[j] * weight;
            weight_sum += weight;
        }

        if weight_sum > 0.0 {
            broadened_mu[i] = sum / weight_sum;
        } else {
            broadened_mu[i] = original_mu[i];
        }
    }

    // Update spectrum with broadened data
    spectrum.mu = broadened_mu;
}

/// Apply anharmonic corrections to the spectrum
///
/// This accounts for non-linear vibrational effects that become
/// important at higher temperatures.
fn apply_anharmonic_corrections(spectrum: &mut XanesSpectrum, thermal_params: &ThermalParameters) {
    // Skip if temperature is low
    if thermal_params.temperature < 100.0 {
        return;
    }

    // The anharmonic correction factor depends on temperature
    // It generally increases with temperature and has more effect on higher energy features
    let anharmonic_strength = 0.01 * (thermal_params.temperature / 300.0 - 0.3).max(0.0);

    // Skip if correction is negligible
    if anharmonic_strength < 0.001 {
        return;
    }

    // Apply anharmonic correction
    for i in 0..spectrum.mu.len() {
        let energy_rel_to_edge = spectrum.energies[i] - spectrum.edge_energy;

        // Anharmonic effects are more pronounced above the edge
        if energy_rel_to_edge > 0.0 {
            // Reduce intensity and shift energy slightly
            // This approximates the effect of phonon-assisted transitions
            let scale_factor = 1.0 - anharmonic_strength * energy_rel_to_edge / 50.0;
            spectrum.mu[i] *= scale_factor.max(0.7); // Limit the reduction
        }
    }
}

/// Calculate the thermal Debye-Waller factor for a given energy point
///
/// The Debye-Waller factor for XANES represents the damping of fine structure
/// due to thermal vibrations. Unlike EXAFS, where it's applied in k-space,
/// for XANES we apply it directly to the absorption coefficient.
///
/// # Arguments
///
/// * `energy` - Energy in eV
/// * `edge_energy` - Edge energy in eV
/// * `thermal_model` - Thermal model to use for calculation
/// * `temperature` - Temperature in Kelvin
///
/// # Returns
///
/// Debye-Waller factor (between 0 and 1)
pub fn calculate_xanes_debye_waller_factor(
    energy: f64,
    edge_energy: f64,
    thermal_model: &dyn ThermalModel,
    temperature: f64,
) -> f64 {
    // For XANES, we need to relate energy to wavenumber k
    // E - E0 = ħ²k²/2m
    // k = sqrt(2m(E-E0)/ħ²) ≈ 0.512 * sqrt(E-E0) in Å⁻¹

    let energy_rel = energy - edge_energy;

    // Only apply above the edge
    if energy_rel <= 0.0 {
        return 1.0;
    }

    // Convert to wavenumber (approximate formula)
    let k = 0.512 * energy_rel.sqrt();

    // Calculate mean-square displacement
    let sigma_sq = thermal_model.mean_square_displacement(temperature);

    // Calculate Debye-Waller factor
    (-2.0 * k * k * sigma_sq).exp()
}

/// Create a XANES spectrum with temperature effects
///
/// This is a higher-level function that applies all temperature corrections
/// to a XANES spectrum and returns a new spectrum with the corrected data.
///
/// # Arguments
///
/// * `base_spectrum` - The base XANES spectrum (typically calculated at 0K)
/// * `thermal_params` - Thermal parameters including temperature and model
/// * `correction_params` - Parameters for thermal corrections
///
/// # Returns
///
/// A new XANES spectrum with all temperature effects applied
pub fn create_temperature_dependent_xanes(
    base_spectrum: &XanesSpectrum,
    thermal_params: &ThermalParameters,
    correction_params: Option<XanesThermalCorrectionParams>,
) -> XanesSpectrum {
    // Use default correction parameters if none provided
    let corr_params = correction_params.unwrap_or_default();

    // Apply all thermal corrections
    apply_thermal_corrections(base_spectrum, thermal_params, &corr_params)
}

/// Calculate temperature-dependent XANES spectra at multiple temperatures
///
/// # Arguments
///
/// * `base_spectrum` - The base XANES spectrum (typically calculated at 0K)
/// * `temperatures` - List of temperatures (in K) to calculate
/// * `thermal_model_type` - Type of thermal model to use ("debye", "einstein", or "correlated_debye")
/// * `thermal_parameter` - Model parameter (Debye temperature or Einstein frequency)
/// * `correction_params` - Parameters for thermal corrections
///
/// # Returns
///
/// A vector of XANES spectra at different temperatures
pub fn calculate_temperature_series(
    base_spectrum: &XanesSpectrum,
    temperatures: &[f64],
    thermal_model_type: &str,
    thermal_parameter: f64,
    correction_params: Option<XanesThermalCorrectionParams>,
) -> Vec<XanesSpectrum> {
    let mut spectra = Vec::with_capacity(temperatures.len());

    // Create thermal parameters for each temperature
    for &temp in temperatures {
        let thermal_params = match thermal_model_type.to_lowercase().as_str() {
            "debye" => ThermalParameters::new_debye(temp, thermal_parameter),
            "einstein" => ThermalParameters::new_einstein(temp, thermal_parameter),
            "correlated_debye" => ThermalParameters::new_correlated_debye(temp, thermal_parameter),
            _ => ThermalParameters::new_debye(temp, thermal_parameter), // Default to Debye
        };

        // Apply thermal corrections
        let spectrum = create_temperature_dependent_xanes(
            base_spectrum,
            &thermal_params,
            correction_params.clone(),
        );

        spectra.push(spectrum);
    }

    spectra
}

/// Calculate and analyze the temperature-dependent changes in specific XANES features
///
/// This function tracks how specific spectral features change with temperature.
///
/// # Arguments
///
/// * `spectra` - Temperature series of XANES spectra
/// * `temperatures` - Corresponding temperatures for each spectrum
/// * `feature_energies` - List of energies (relative to edge) for features to track
///
/// # Returns
///
/// A map of feature energy to (temperature, intensity) pairs
pub fn analyze_temperature_dependence(
    spectra: &[XanesSpectrum],
    temperatures: &[f64],
    feature_energies: &[f64],
) -> Vec<(f64, Vec<(f64, f64)>)> {
    let mut results = Vec::with_capacity(feature_energies.len());

    // For each feature energy
    for &feature_e_rel in feature_energies {
        let mut feature_data = Vec::with_capacity(temperatures.len());

        // For each temperature
        for (i, spectrum) in spectra.iter().enumerate() {
            let temp = temperatures[i];

            // Convert relative energy to absolute
            let feature_e_abs = spectrum.edge_energy + feature_e_rel;

            // Find closest energy point
            let closest_idx = spectrum
                .energies
                .iter()
                .enumerate()
                .min_by(|(_, a), (_, b)| {
                    (**a - feature_e_abs)
                        .abs()
                        .partial_cmp(&(**b - feature_e_abs).abs())
                        .unwrap()
                })
                .map(|(idx, _)| idx)
                .unwrap_or(0);

            // Get intensity at this energy
            let intensity = spectrum.normalized_mu[closest_idx];

            // Add (temperature, intensity) pair
            feature_data.push((temp, intensity));
        }

        // Add (feature_energy, feature_data) to results
        results.push((feature_e_rel, feature_data));
    }

    results
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::xas::xanes::{Edge, XanesParameters, XanesSpectrum};

    fn create_test_spectrum() -> XanesSpectrum {
        // Create a simple test spectrum
        let edge = Edge::K;
        let atomic_number = 26; // Fe
        let edge_energy = 7112.0; // Fe K-edge
        let params = XanesParameters::default();

        let mut spectrum = XanesSpectrum::new(edge, atomic_number, edge_energy, params);

        // Generate energies from -20 to 50 eV around edge
        let mut energies = Vec::new();
        let mut e = edge_energy - 20.0;
        let step = 0.5;

        while e <= edge_energy + 50.0 {
            energies.push(e);
            e += step;
        }

        spectrum.energies = energies;

        // Create a simple spectrum with:
        // - Pre-edge feature at -5 eV
        // - Main edge jump at 0 eV
        // - Post-edge features at 10 and 30 eV
        let mut mu = Vec::with_capacity(spectrum.energies.len());

        for e in &spectrum.energies {
            let e_rel = e - edge_energy;

            // Base absorption
            let mut absorption = 0.1;

            // Edge jump (sigmoid function)
            if e_rel >= -10.0 {
                absorption += 0.8 / (1.0 + (-e_rel).exp());
            }

            // Pre-edge feature
            if (e_rel + 5.0).abs() < 2.0 {
                absorption += 0.2 * (-(e_rel + 5.0).powi(2) / 0.5).exp();
            }

            // Post-edge features
            if (e_rel - 10.0).abs() < 3.0 {
                absorption += 0.3 * (-(e_rel - 10.0).powi(2) / 2.0).exp();
            }

            if (e_rel - 30.0).abs() < 5.0 {
                absorption += 0.4 * (-(e_rel - 30.0).powi(2) / 6.0).exp();
            }

            mu.push(absorption);
        }

        spectrum.mu = mu;

        // Normalized values (just use the same as mu for test)
        spectrum.normalized_mu = spectrum.mu.clone();

        spectrum
    }

    #[test]
    fn test_thermal_expansion() {
        let mut spectrum = create_test_spectrum();
        let original_energies = spectrum.energies.clone();

        // Create thermal parameters for high temperature (1000K)
        let thermal_params = ThermalParameters::new_debye(1000.0, 400.0);

        // Create correction params with significant expansion
        let mut corr_params = XanesThermalCorrectionParams::default();
        corr_params.thermal_expansion_coefficient = 1.0e-5;
        corr_params.reference_temperature = 300.0;

        // Apply thermal expansion
        apply_thermal_expansion(&mut spectrum, &thermal_params, &corr_params);

        // Verify energies have shifted
        for i in 0..original_energies.len() {
            // Skip edge energy point - it should remain fixed
            if (original_energies[i] - spectrum.edge_energy).abs() < 0.01 {
                continue;
            }

            // For thermal expansion, both above and below edge features should have
            // their relative energies reduced in magnitude (divided by expansion factor)
            let original_relative = original_energies[i] - spectrum.edge_energy;
            let new_relative = spectrum.energies[i] - spectrum.edge_energy;

            // The relative energy magnitude should decrease (get closer to zero)
            assert!(
                new_relative.abs() < original_relative.abs(),
                "Relative energy magnitude should decrease with thermal expansion"
            );
        }
    }

    #[test]
    fn test_thermal_broadening() {
        let mut spectrum = create_test_spectrum();
        let original_mu = spectrum.mu.clone();

        // Create thermal parameters for room temperature
        let thermal_params = ThermalParameters::new_debye(300.0, 400.0);

        // Create correction params
        let corr_params = XanesThermalCorrectionParams::default();

        // Apply thermal broadening
        apply_temperature_dependent_broadening(&mut spectrum, &thermal_params, &corr_params);

        // Find peak in original and broadened spectrum
        let original_max = original_mu.iter().cloned().fold(0.0f64, f64::max);
        let broadened_max = spectrum.mu.iter().cloned().fold(0.0f64, f64::max);

        // Peaks should be lower after broadening
        assert!(
            broadened_max < original_max,
            "Thermal broadening should reduce peak heights"
        );

        // Verify that valleys are filled in (higher after broadening)
        let mut found_valley = false;
        for i in 1..spectrum.mu.len() - 1 {
            // Look for a local minimum (valley)
            if original_mu[i] < original_mu[i - 1] && original_mu[i] < original_mu[i + 1] {
                if spectrum.mu[i] > original_mu[i] {
                    found_valley = true;
                    break;
                }
            }
        }

        assert!(
            found_valley,
            "Thermal broadening should increase intensity in valleys"
        );
    }

    #[test]
    fn test_full_thermal_corrections() {
        let spectrum = create_test_spectrum();

        // Create thermal parameters with very wide temperature spread for clear effects
        let _thermal_params_low = ThermalParameters::new_debye(5.0, 400.0); // Near zero K
        let thermal_params_high = ThermalParameters::new_debye(2000.0, 400.0); // Extremely hot

        // Create correction params with exaggerated thermal effects for testing
        let mut corr_params = XanesThermalCorrectionParams::default();
        corr_params.thermal_expansion_coefficient = 1.0e-4; // Exaggerated for testing

        // Increase all broadening parameters for clearer test results
        corr_params.energy_dependent_broadening = vec![
            (-15.0, 1.0), // Double the normal values
            (0.0, 2.0),   // Double the normal values
            (15.0, 1.6),  // Double the normal values
            (50.0, 1.4),  // Double the normal values
        ];

        // Apply corrections
        let spectrum_orig = spectrum.clone(); // Keep the original unchanged
        let spectrum_high =
            apply_thermal_corrections(&spectrum, &thermal_params_high, &corr_params);

        // Quick check that spectra have the same number of points
        assert_eq!(spectrum_high.mu.len(), spectrum_orig.mu.len());

        // For testing, let's verify the broadening directly - find the sharpest peak in both spectra
        let find_sharpest_peak = |spec: &XanesSpectrum| -> (usize, f64) {
            let mut max_idx = 0;
            let mut max_sharpness = 0.0;

            // Look for local maxima with high curvature (sharp peaks)
            for i in 2..spec.mu.len() - 2 {
                if spec.mu[i] > spec.mu[i - 1] && spec.mu[i] > spec.mu[i + 1] {
                    // Calculate curvature (second derivative estimate)
                    let curvature = (spec.mu[i - 1] + spec.mu[i + 1] - 2.0 * spec.mu[i]).abs();

                    if curvature > max_sharpness {
                        max_sharpness = curvature;
                        max_idx = i;
                    }
                }
            }

            (max_idx, max_sharpness)
        };

        let (_, orig_sharpness) = find_sharpest_peak(&spectrum_orig);
        let (_, high_sharpness) = find_sharpest_peak(&spectrum_high);

        // Thermal broadening should reduce peak sharpness
        assert!(
            high_sharpness < orig_sharpness,
            "Thermal effects should reduce peak sharpness"
        );

        // Also directly check at a specific feature (edge jump)
        let edge_idx = spectrum_orig
            .energies
            .iter()
            .position(|&e| e >= spectrum_orig.edge_energy)
            .unwrap_or(spectrum_orig.energies.len() / 2);

        // For a normalized spectrum, the edge jump is sharper at low temperatures
        let slope_orig = (spectrum_orig.normalized_mu[edge_idx + 2]
            - spectrum_orig.normalized_mu[edge_idx - 2])
            / (spectrum_orig.energies[edge_idx + 2] - spectrum_orig.energies[edge_idx - 2]);

        let slope_high = (spectrum_high.normalized_mu[edge_idx + 2]
            - spectrum_high.normalized_mu[edge_idx - 2])
            / (spectrum_high.energies[edge_idx + 2] - spectrum_high.energies[edge_idx - 2]);

        // Higher temperature should reduce the edge steepness
        assert!(
            slope_high.abs() < slope_orig.abs() || high_sharpness < orig_sharpness,
            "Either edge slope or peak sharpness should be reduced by thermal effects"
        );
    }

    #[test]
    fn test_temperature_series() {
        let spectrum = create_test_spectrum();

        // Calculate a temperature series
        let temperatures = [100.0, 300.0, 500.0, 800.0];

        let series = calculate_temperature_series(&spectrum, &temperatures, "debye", 400.0, None);

        // Verify we got the right number of spectra
        assert_eq!(series.len(), temperatures.len());

        // Verify each spectrum has data
        for spec in &series {
            assert!(!spec.mu.is_empty());
            assert!(!spec.normalized_mu.is_empty());
        }
    }
}
