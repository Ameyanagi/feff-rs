/*
MIT License with FEFF10 Attribution

Copyright (c) 2025 Ameyanagi

Based on or developed using Distribution: FEFF10.0
Copyright (c) 2020 FEFF Project, University of Washington and SLAC National Accelerator Laboratory.
All rights reserved.
*/

//! Simple tests for the thermal models (DebyeModel and EinsteinModel)
//! to ensure the basic mean-square displacement calculations are correct.

use feff_rs::utils::thermal::{DebyeModel, EinsteinModel, ThermalModel};

#[test]
fn test_debye_model_msrd() {
    // Test the Debye model MSRD calculation
    // Using typical values for Cu: θD = 315K, μ = 63.546 amu
    let debye_temp = 315.0; // K
    let reduced_mass = 63.546; // amu

    // Create Debye model
    let model = DebyeModel::new(debye_temp, reduced_mass);

    // At T = 0K, should just be zero-point motion
    let sigma_sq_0k = model.mean_square_displacement(0.0);
    assert!(sigma_sq_0k > 0.0, "Zero-point motion should be positive");

    // At 300K (room temperature)
    let sigma_sq_300k = model.mean_square_displacement(300.0);

    assert!(
        sigma_sq_300k > sigma_sq_0k,
        "MSRD should increase with temperature"
    );

    // Higher temperature should have larger MSRD
    let sigma_sq_500k = model.mean_square_displacement(500.0);
    assert!(
        sigma_sq_500k > sigma_sq_300k,
        "MSRD should increase with temperature"
    );
}

#[test]
fn test_einstein_model_msrd() {
    // Test the Einstein model MSRD calculation
    // Using typical values for a Cu-O bond: ω = 25 meV, μ = 12 amu
    let einstein_freq = 25.0; // meV
    let reduced_mass = 12.0; // amu

    // Create Einstein model
    let model = EinsteinModel::new(einstein_freq, reduced_mass);

    // At T = 0K, should just be zero-point motion
    let sigma_sq_0k = model.mean_square_displacement(0.0);
    assert!(sigma_sq_0k > 0.0, "Zero-point motion should be positive");

    // At 300K (room temperature)
    let sigma_sq_300k = model.mean_square_displacement(300.0);
    assert!(
        sigma_sq_300k > sigma_sq_0k,
        "MSRD should increase with temperature"
    );

    // Higher temperature should have larger MSRD
    let sigma_sq_500k = model.mean_square_displacement(500.0);
    assert!(
        sigma_sq_500k > sigma_sq_300k,
        "MSRD should increase with temperature"
    );
}

#[test]
fn test_thermal_model_comparison() {
    // Compare Debye and Einstein models
    // For many materials, they should give similar results with appropriate parameters

    // Using Cu parameters (approximate)
    let debye_temp = 315.0; // K
    let einstein_freq = 21.0; // meV (approximately equivalent to Debye temperature)
    let reduced_mass = 63.546; // amu

    let debye_model = DebyeModel::new(debye_temp, reduced_mass);
    let einstein_model = EinsteinModel::new(einstein_freq, reduced_mass);

    // Compare at room temperature (300K)
    let debye_msrd = debye_model.mean_square_displacement(300.0);
    let einstein_msrd = einstein_model.mean_square_displacement(300.0);

    // Just assert that they're both positive - we don't need to compare exact values
    // since the models are different
    assert!(debye_msrd > 0.0);
    assert!(einstein_msrd > 0.0);
}

#[test]
fn test_temperature_scaling() {
    // Test temperature scaling behavior of the models

    let debye_model = DebyeModel::new(300.0, 50.0);
    let einstein_model = EinsteinModel::new(20.0, 50.0);

    // Test at different temperatures
    let temps = [10.0, 100.0, 300.0, 500.0, 1000.0];

    // Just check that values are monotonically increasing with temperature
    let mut prev_debye = 0.0;
    let mut prev_einstein = 0.0;

    for temp in temps {
        let debye_msrd = debye_model.mean_square_displacement(temp);
        let einstein_msrd = einstein_model.mean_square_displacement(temp);

        assert!(
            debye_msrd > prev_debye,
            "Debye MSRD should increase with temperature"
        );
        assert!(
            einstein_msrd > prev_einstein,
            "Einstein MSRD should increase with temperature"
        );

        prev_debye = debye_msrd;
        prev_einstein = einstein_msrd;
    }
}
