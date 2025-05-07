/*
MIT License with FEFF10 Attribution

Copyright (c) 2025 Ameyanagi

Based on or developed using Distribution: FEFF10.0
Copyright (c) 2020 FEFF Project, University of Washington and SLAC National Accelerator Laboratory.
All rights reserved.
*/

use feff_rs::atoms::{Atom, AtomicStructure, PotentialType, Vector3D};
use feff_rs::utils::thermal::ThermalModel;
use feff_rs::xas::{
    apply_thermal_corrections, calculate_temperature_series, thermal::ThermalParameters, Edge,
    XanesParameters, XanesSpectrum, XanesThermalCorrectionParams,
};

// Create a simple mock-up of a XANES spectrum for testing
fn create_test_spectrum() -> XanesSpectrum {
    let edge = Edge::K;
    let atomic_number = 26; // Fe
    let edge_energy = 7112.0; // Fe K-edge
    let params = XanesParameters::default();

    let mut spectrum = XanesSpectrum::new(edge, atomic_number, edge_energy, params);

    // Generate energy grid from -20 to 50 eV around edge
    let mut energies = Vec::new();
    let e_start = edge_energy - 20.0;
    let e_step = 0.5;
    let e_end = edge_energy + 50.0;

    let mut e = e_start;
    while e <= e_end {
        energies.push(e);
        e += e_step;
    }

    spectrum.energies = energies;

    // Create a mock spectrum with:
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
fn test_basic_thermal_corrections() {
    // Create a test spectrum
    let spectrum = create_test_spectrum();

    // Create thermal parameters for room temperature
    let thermal_params = ThermalParameters::new_debye(300.0, 315.0);

    // Create correction parameters
    let correction_params = XanesThermalCorrectionParams::default();

    // Apply thermal corrections
    let corrected = apply_thermal_corrections(&spectrum, &thermal_params, &correction_params);

    // Ensure corrected spectrum has expected dimensions
    assert_eq!(corrected.mu.len(), spectrum.mu.len());
    assert_eq!(corrected.normalized_mu.len(), spectrum.normalized_mu.len());

    // Thermal effects should reduce sharp features
    let mut found_broadening = false;
    for i in 0..spectrum.mu.len() {
        // Look for local maxima in the original spectrum
        if i > 0
            && i < spectrum.mu.len() - 1
            && spectrum.mu[i] > spectrum.mu[i - 1]
            && spectrum.mu[i] > spectrum.mu[i + 1]
        {
            // Local peak should be lower in the broadened spectrum
            if corrected.mu[i] < spectrum.mu[i] {
                found_broadening = true;
                break;
            }
        }
    }

    assert!(
        found_broadening,
        "Thermal corrections should reduce peak intensities"
    );
}

#[test]
fn test_temperature_dependence() {
    // Create a test spectrum
    let spectrum = create_test_spectrum();

    // Calculate spectra at very different temperatures to ensure clear effects
    let temperatures = [50.0, 2000.0]; // Use extreme temperatures for clear testing

    // Create custom correction params with exaggerated thermal effects
    let mut corr_params = XanesThermalCorrectionParams::default();
    corr_params.thermal_expansion_coefficient = 1.0e-4; // Exaggerated for testing

    // Increase all broadening parameters for clearer test results
    corr_params.energy_dependent_broadening = vec![
        (-15.0, 1.0), // Double the normal values
        (0.0, 2.0),   // Double the normal values
        (15.0, 1.6),  // Double the normal values
        (50.0, 1.4),  // Double the normal values
    ];

    let temp_series = calculate_temperature_series(
        &spectrum,
        &temperatures,
        "debye",
        315.0, // Debye temperature for Fe
        Some(corr_params),
    );

    // Should have one spectrum per temperature
    assert_eq!(temp_series.len(), temperatures.len());

    // For testing, evaluate the overall spectrum smoothness
    // Higher temperatures cause more broadening, which reduces overall "roughness"
    let calculate_roughness = |spec: &XanesSpectrum| -> f64 {
        let mut roughness = 0.0;

        // Sum up absolute differences between adjacent points
        for i in 1..spec.normalized_mu.len() {
            roughness += (spec.normalized_mu[i] - spec.normalized_mu[i - 1]).abs();
        }

        roughness
    };

    let low_temp_roughness = calculate_roughness(&temp_series[0]);
    let high_temp_roughness = calculate_roughness(&temp_series[1]);

    // Higher temperature should reduce the overall "roughness" (jaggedness)
    assert!(
        high_temp_roughness < low_temp_roughness,
        "Higher temperature should smooth the spectrum (reduce roughness)"
    );
}

#[test]
fn test_thermal_expansion_effects() {
    let spectrum = create_test_spectrum();

    // Create thermal parameters for high temperature with exaggerated expansion
    let thermal_params = ThermalParameters::new_debye(1000.0, 315.0);
    let mut corr_params = XanesThermalCorrectionParams::default();
    corr_params.thermal_expansion_coefficient = 5.0e-5; // Exaggerated for testing

    // Make other thermal effects minimal to focus only on expansion
    corr_params.energy_dependent_broadening = vec![
        (-20.0, 0.01), // Very little broadening
        (50.0, 0.01),  // Very little broadening
    ];
    corr_params.use_anharmonic_corrections = false;

    // Apply thermal corrections using the public API
    let expansion_only = apply_thermal_corrections(&spectrum, &thermal_params, &corr_params);

    // Calculate expected behavior for thermal expansion
    let edge_energy = spectrum.edge_energy;
    let delta_temp = thermal_params.temperature - corr_params.reference_temperature;
    let expansion_factor = 1.0 + corr_params.thermal_expansion_coefficient * delta_temp;

    // Print debug info
    println!(
        "Thermal expansion test: ΔT={}, factor={}",
        delta_temp, expansion_factor
    );

    // Test a few specific points
    for i in 0..spectrum.energies.len() {
        let orig_rel = spectrum.energies[i] - edge_energy;
        let scaled_rel = expansion_only.energies[i] - edge_energy;

        // Skip points very close to edge
        if orig_rel.abs() < 0.1 {
            continue;
        }

        // Print some sample points for debugging
        if i % 20 == 0 {
            println!(
                "Rel energy: orig={:.4}, scaled={:.4}, ratio={:.4}",
                orig_rel,
                scaled_rel,
                scaled_rel / orig_rel
            );
        }

        // Verify energy was scaled correctly - for any energy (above or below edge),
        // the relative energy should be divided by the expansion factor
        let expected_rel = orig_rel / expansion_factor;
        assert!(
            (scaled_rel - expected_rel).abs() < 1e-5, // Increased tolerance due to potential numerical issues
            "Energy scaling incorrect at point {}: expected rel {}, got {}",
            i,
            expected_rel,
            scaled_rel
        );

        // For above-edge points, relative energy should decrease in magnitude
        if orig_rel > 0.0 {
            assert!(
                scaled_rel < orig_rel,
                "Above-edge energy should decrease from {} to {}",
                orig_rel,
                scaled_rel
            );
        }
        // For below-edge points, relative energy should also decrease in magnitude
        // but since it's negative, the new value will be higher (less negative)
        else if orig_rel < 0.0 {
            assert!(
                scaled_rel > orig_rel,
                "Below-edge energy should increase from {} to {}",
                orig_rel,
                scaled_rel
            );
        }
    }
}

#[test]
fn test_different_thermal_models() {
    let spectrum = create_test_spectrum();

    // Create thermal parameters for different models
    let debye_params = ThermalParameters::new_debye(300.0, 315.0);
    let einstein_params = ThermalParameters::new_einstein(300.0, 25.0);
    let corr_debye_params = ThermalParameters::new_correlated_debye(300.0, 315.0);

    // Apply thermal corrections with each model
    let corr_params = XanesThermalCorrectionParams::default();

    let debye_spectrum = apply_thermal_corrections(&spectrum, &debye_params, &corr_params);
    let einstein_spectrum = apply_thermal_corrections(&spectrum, &einstein_params, &corr_params);
    let corr_debye_spectrum =
        apply_thermal_corrections(&spectrum, &corr_debye_params, &corr_params);

    // All models should produce valid results
    assert_eq!(debye_spectrum.mu.len(), spectrum.mu.len());
    assert_eq!(einstein_spectrum.mu.len(), spectrum.mu.len());
    assert_eq!(corr_debye_spectrum.mu.len(), spectrum.mu.len());

    // Einstein model typically produces more broadening than Debye at same temperature
    // Find a peak to compare
    for i in 10..spectrum.mu.len() - 10 {
        if spectrum.mu[i] > spectrum.mu[i - 1] && spectrum.mu[i] > spectrum.mu[i + 1] {
            // At a peak, check if Einstein model typically gives more broadening
            let debye_ratio = debye_spectrum.mu[i] / spectrum.mu[i];
            let einstein_ratio = einstein_spectrum.mu[i] / spectrum.mu[i];

            assert!(
                einstein_ratio <= debye_ratio,
                "Einstein model should typically produce more or equal broadening than Debye"
            );
            break;
        }
    }
}
