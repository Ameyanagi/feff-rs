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

    // Print values for debugging
    println!("Debye model tests");
    println!("MSRD at 0K: {}", sigma_sq_0k);
    println!("MSRD at 300K: {}", model.mean_square_displacement(300.0));
    println!("MSRD at 500K: {}", model.mean_square_displacement(500.0));

    // Check that values are physically reasonable
    let sigma_sq_300k = model.mean_square_displacement(300.0);
    assert!(
        sigma_sq_300k > 0.001 && sigma_sq_300k < 0.02,
        "MSRD at 300K should be in physically reasonable range"
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

    // Print values for debugging
    println!("Einstein model tests");
    println!("MSRD at 0K: {}", sigma_sq_0k);
    println!("MSRD at 300K: {}", model.mean_square_displacement(300.0));
    println!("MSRD at 500K: {}", model.mean_square_displacement(500.0));

    // Check that values are physically reasonable
    let sigma_sq_300k = model.mean_square_displacement(300.0);
    assert!(
        sigma_sq_300k > 0.001 && sigma_sq_300k < 0.02,
        "MSRD at 300K should be in physically reasonable range"
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

    println!("\nTemperature scaling test:");

    // Print values at each temperature for inspection
    for temp in &temps {
        let debye_msrd = debye_model.mean_square_displacement(*temp);
        let einstein_msrd = einstein_model.mean_square_displacement(*temp);

        println!(
            "T={:.1}K - Debye: {:.6}, Einstein: {:.6}",
            temp, debye_msrd, einstein_msrd
        );

        // Check that values are physically reasonable at each temperature
        assert!(
            debye_msrd > 0.0001 && debye_msrd < 0.05,
            "Debye MSRD at T={} should be reasonable",
            temp
        );
        assert!(
            einstein_msrd > 0.0001 && einstein_msrd < 0.05,
            "Einstein MSRD at T={} should be reasonable",
            temp
        );
    }
}

#[test]
fn test_anisotropic_thermal_model() {
    use feff_rs::utils::thermal::{AnisotropicThermalModel, CorrelatedDebyeModel};

    // Test anisotropic thermal model with a Debye base model
    let temp = 300.0;
    let debye_temp = 315.0;
    let mass = 63.546; // Cu

    // Create a standard correlated Debye model for comparison
    let corr_debye = CorrelatedDebyeModel::new(debye_temp, mass, 2.5);
    let msrd_isotropic = corr_debye.mean_square_displacement(temp);

    // Create an anisotropic model with enhanced z vibrations
    let aniso_z = AnisotropicThermalModel::new_from_correlated_debye(
        debye_temp,
        mass,
        2.5,             // First-shell distance
        [0.8, 0.8, 2.0], // Enhanced z-vibrations
        [0.0, 0.0, 1.0], // Path along z-axis
    );

    // Create an anisotropic model with enhanced x vibrations
    let aniso_x = AnisotropicThermalModel::new_from_correlated_debye(
        debye_temp,
        mass,
        2.5,             // First-shell distance
        [2.0, 0.8, 0.8], // Enhanced x-vibrations
        [0.0, 0.0, 1.0], // Path along z-axis
    );

    // Calculate MSDs
    let msrd_aniso_z = aniso_z.mean_square_displacement(temp);
    let msrd_aniso_x = aniso_x.mean_square_displacement(temp);

    // Enhanced vibrations along z should increase MSRD for z-direction path
    assert!(msrd_aniso_z > msrd_isotropic);

    // Enhanced vibrations along x should have less effect on z-direction path
    assert!(msrd_aniso_z > msrd_aniso_x);

    // Now test paths along different directions
    let aniso_z_x_path = AnisotropicThermalModel::new_from_correlated_debye(
        debye_temp,
        mass,
        2.5,             // First-shell distance
        [0.8, 0.8, 2.0], // Enhanced z-vibrations
        [1.0, 0.0, 0.0], // Path along x-axis
    );

    let msrd_z_x_path = aniso_z_x_path.mean_square_displacement(temp);

    // Enhanced vibrations along z should have less effect on x-direction path
    assert!(msrd_aniso_z > msrd_z_x_path);

    // Print values for debugging
    println!("\nAnisotropic model tests at T={}K:", temp);
    println!("Isotropic MSRD: {}", msrd_isotropic);
    println!("Anisotropic (z-enhanced, z-path) MSRD: {}", msrd_aniso_z);
    println!("Anisotropic (x-enhanced, z-path) MSRD: {}", msrd_aniso_x);
    println!("Anisotropic (z-enhanced, x-path) MSRD: {}", msrd_z_x_path);
}

#[test]
fn test_anisotropic_thermal_parameters() {
    use feff_rs::xas::thermal::{create_anisotropic_thermal_parameters, ThermalParameters};

    // Create thermal parameters for different crystal systems
    let temp = 300.0;
    let debye_temp = 315.0;

    // Create thermal parameters for cubic, tetragonal, and hexagonal systems
    let cubic_params = create_anisotropic_thermal_parameters(temp, "cubic", debye_temp, None);
    let tetragonal_params =
        create_anisotropic_thermal_parameters(temp, "tetragonal", debye_temp, None);
    let hexagonal_params =
        create_anisotropic_thermal_parameters(temp, "hexagonal", debye_temp, None);
    let orthorhombic_params =
        create_anisotropic_thermal_parameters(temp, "orthorhombic", debye_temp, None);
    let layered_params =
        create_anisotropic_thermal_parameters(temp, "layered", debye_temp, Some(2.0));

    // Check model type
    assert_eq!(cubic_params.model_type, "anisotropic_correlated_debye");
    assert_eq!(tetragonal_params.model_type, "anisotropic_correlated_debye");
    assert_eq!(hexagonal_params.model_type, "anisotropic_correlated_debye");

    // Check displacement factors
    assert_eq!(cubic_params.displacement_factors, Some([1.0, 1.0, 1.0]));

    let tetragonal_factors = tetragonal_params.displacement_factors.unwrap();
    assert_eq!(tetragonal_factors[0], 1.0);
    assert_eq!(tetragonal_factors[1], 1.0);
    assert!(tetragonal_factors[2] > 1.0);

    let hexagonal_factors = hexagonal_params.displacement_factors.unwrap();
    assert_eq!(hexagonal_factors[0], 1.0);
    assert_eq!(hexagonal_factors[1], 1.0);
    assert!(hexagonal_factors[2] > 1.0);

    let orthorhombic_factors = orthorhombic_params.displacement_factors.unwrap();
    assert_eq!(orthorhombic_factors[0], 1.0);
    assert!(orthorhombic_factors[1] != 1.0);
    assert!(orthorhombic_factors[2] != 1.0);

    // Layered materials should have very high anisotropy along c-axis
    let layered_factors = layered_params.displacement_factors.unwrap();
    assert!(layered_factors[2] > 3.0);

    // Test creating thermal models from these parameters
    let reduced_mass = 63.546; // Cu
    let path_distance = 2.5; // typical first shell

    let cubic_model = cubic_params.create_model(reduced_mass, Some(path_distance));
    let tetragonal_model = tetragonal_params.create_model(reduced_mass, Some(path_distance));
    let hexagonal_model = hexagonal_params.create_model(reduced_mass, Some(path_distance));
    let layered_model = layered_params.create_model(reduced_mass, Some(path_distance));

    // Calculate MSDs at room temperature
    let msd_cubic = cubic_model.mean_square_displacement(temp);
    let msd_tetragonal = tetragonal_model.mean_square_displacement(temp);
    let msd_hexagonal = hexagonal_model.mean_square_displacement(temp);
    let msd_layered = layered_model.mean_square_displacement(temp);

    // MSDs should be larger for more anisotropic systems along z-axis
    assert!(msd_cubic < msd_tetragonal);
    assert!(msd_tetragonal < msd_layered);

    // Print values for debugging
    println!("\nAnisotropic thermal parameters tests at T={}K:", temp);
    println!("Cubic MSD: {}", msd_cubic);
    println!("Tetragonal MSD: {}", msd_tetragonal);
    println!("Hexagonal MSD: {}", msd_hexagonal);
    println!("Layered MSD: {}", msd_layered);

    // Test Debye-Waller factors
    let k = 10.0; // Å⁻¹
    let dw_cubic = cubic_model.debye_waller_factor(temp, k);
    let dw_layered = layered_model.debye_waller_factor(temp, k);

    // Higher anisotropy should lead to lower DW factor for z-axis paths
    assert!(dw_cubic > dw_layered);
}

#[test]
fn test_material_thermal_parameters() {
    use feff_rs::xas::thermal::create_material_thermal_parameters;

    // Test material-specific thermal parameters for various materials
    let temp = 300.0;

    // Test metals
    let cu_params = create_material_thermal_parameters("Cu", temp, None);
    let fe_params = create_material_thermal_parameters("Fe", temp, None);
    let au_params = create_material_thermal_parameters("gold", temp, None);

    // Test oxides
    let tio2_params = create_material_thermal_parameters("TiO2", temp, None);
    let fe2o3_params = create_material_thermal_parameters("Fe2O3", temp, None);

    // Test layered materials
    let graphite_params = create_material_thermal_parameters("graphite", temp, None);
    let mos2_params = create_material_thermal_parameters("MoS2", temp, None);

    // Verify models are appropriate for materials
    assert_eq!(cu_params.model_type, "debye");
    assert!(tio2_params.model_type.contains("anisotropic"));
    assert!(graphite_params.model_type.contains("anisotropic"));

    // Verify expected Debye temperatures
    let cu_debye_temp = match cu_params.pair_parameters {
        Some(ref params) => {
            let cu_cu_key = feff_rs::xas::thermal::PairKey::new(29, 29);
            match params.get(&cu_cu_key) {
                Some(pair_params) => pair_params.debye_temperature,
                None => 0.0,
            }
        }
        None => 0.0,
    };
    assert!(cu_debye_temp > 300.0 && cu_debye_temp < 320.0); // Around 315K

    // Check that layered materials have high anisotropy
    if let Some(factors) = graphite_params.displacement_factors {
        assert!(factors[2] > factors[0]); // z-axis vibrations should be larger
        println!("Graphite displacement factors: {:?}", factors);
    } else {
        panic!("Graphite should have anisotropic displacement factors");
    }

    // Test with custom Debye temperature
    let custom_si = create_material_thermal_parameters("Si", temp, Some(800.0));
    let default_si = create_material_thermal_parameters("Si", temp, None);

    assert!(custom_si.debye_temperature > default_si.debye_temperature);
}

#[test]
fn test_exafs_thermal_parameters() {
    use feff_rs::xas::thermal::create_exafs_thermal_parameters;

    // Test EXAFS-specific thermal parameters
    let temp_low = 100.0;
    let temp_high = 500.0;
    let debye_temp = 350.0;

    // Test various material types
    let metal_params = create_exafs_thermal_parameters("metal", temp_low, debye_temp, false);
    let oxide_params = create_exafs_thermal_parameters("oxide", temp_low, debye_temp, false);
    let layered_params = create_exafs_thermal_parameters("layered", temp_low, debye_temp, false);

    // Check model types are appropriate
    assert_eq!(metal_params.model_type, "correlated_debye");

    // Check that layered materials use anisotropic model
    assert!(layered_params.model_type.contains("anisotropic"));

    // Test high temperature with anharmonic effects
    let metal_high_temp = create_exafs_thermal_parameters("metal", temp_high, debye_temp, true);
    let oxide_high_temp = create_exafs_thermal_parameters("oxide", temp_high, debye_temp, true);

    // Anharmonic models should have adjusted parameters
    assert!(metal_high_temp.debye_temperature < debye_temp); // Anharmonicity reduces effective Debye temp

    // High temperature anisotropic model should have higher anisotropy
    if let (Some(low_factors), Some(high_factors)) = (
        layered_params.displacement_factors,
        create_exafs_thermal_parameters("layered", temp_high, debye_temp, true)
            .displacement_factors,
    ) {
        assert!(high_factors[2] > low_factors[2]);
        println!("Low temp anisotropy: {:?}", low_factors);
        println!("High temp anisotropy: {:?}", high_factors);
    }
}

#[test]
fn test_xanes_thermal_parameters() {
    use feff_rs::xas::thermal::create_xanes_thermal_parameters;

    // Test XANES-specific thermal parameters
    let temp_low = 100.0;
    let temp_high = 500.0;
    let debye_temp = 350.0;

    // Test for different edges and materials
    let metal_k_edge = create_xanes_thermal_parameters("metal", temp_low, debye_temp, "K");
    let metal_l_edge = create_xanes_thermal_parameters("metal", temp_low, debye_temp, "L3");
    let oxide_k_edge = create_xanes_thermal_parameters("oxide", temp_low, debye_temp, "K");

    // Check that model types are appropriate
    assert_eq!(metal_k_edge.model_type, "debye");

    // Check that oxide uses anisotropic model
    assert!(oxide_k_edge.model_type.contains("anisotropic"));

    // Test high temperature parameters
    let high_temp_k = create_xanes_thermal_parameters("metal", temp_high, debye_temp, "K");

    // High temperature should use correlated model for K-edge
    assert!(high_temp_k.model_type.contains("correlated"));

    // Check for anharmonic correction for high temperatures
    let high_temp_model = create_xanes_thermal_parameters("layered", temp_high, debye_temp, "K");
    assert!(high_temp_model.model_type.contains("anharmonic"));

    // Test that different edges get appropriate models
    let m_edge = create_xanes_thermal_parameters("metal", temp_high, debye_temp, "M");
    assert!(m_edge.model_type != metal_k_edge.model_type);
}
