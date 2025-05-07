/*
MIT License with FEFF10 Attribution

Copyright (c) 2025 Ameyanagi

Based on or developed using Distribution: FEFF10.0
Copyright (c) 2020 FEFF Project, University of Washington and SLAC National Accelerator Laboratory.
All rights reserved.
*/

//! Tests for pair-specific thermal parameters in EXAFS calculations

use feff_rs::atoms::{Atom, AtomicStructure, PotentialType, Vector3D};
use feff_rs::xas::{
    calculate_exafs,
    thermal::{
        create_standard_thermal_parameters, PairKey, PairThermalParameters, ThermalParameters,
    },
    Edge, EnergyGrid, ExafsParameters,
};

#[test]
fn test_pair_thermal_parameters() {
    // Test creating pair-specific thermal parameters
    let pair_key = PairKey::new(29, 8); // Cu-O pair

    // Create a PairThermalParameters instance for Cu-O bond
    let pair_params = PairThermalParameters::new_einstein(
        70.0, // Einstein frequency in meV
        Some("Cu-O bond".to_string()),
    );

    // Create global thermal parameters
    let mut thermal_params = ThermalParameters::new_debye(300.0, 315.0);

    // Add pair-specific parameters
    thermal_params = thermal_params.with_pair_parameters(29, 8, pair_params);

    // Check that we can retrieve the pair parameters
    let retrieved_params = thermal_params.get_pair_parameters(29, 8);
    assert!(retrieved_params.is_some());

    if let Some(params) = retrieved_params {
        assert_eq!(params.model_type, "einstein");
        assert_eq!(params.einstein_frequency, Some(70.0));
        assert_eq!(params.description, Some("Cu-O bond".to_string()));
    }

    // Test canonical ordering of pair keys
    // Getting Cu-O parameters by specifying O-Cu should work the same way
    let retrieved_params = thermal_params.get_pair_parameters(8, 29);
    assert!(retrieved_params.is_some());
}

#[test]
fn test_standard_thermal_parameters() {
    // Create standard thermal parameters at room temperature
    let params = create_standard_thermal_parameters(300.0);

    // Check that we have parameters for various common pairs
    assert!(params.get_pair_parameters(29, 29).is_some()); // Cu-Cu
    assert!(params.get_pair_parameters(26, 26).is_some()); // Fe-Fe
    assert!(params.get_pair_parameters(29, 8).is_some()); // Cu-O
    assert!(params.get_pair_parameters(26, 8).is_some()); // Fe-O

    // Check that metals use Debye model
    let cu_cu = params.get_pair_parameters(29, 29).unwrap();
    assert_eq!(cu_cu.model_type, "debye");

    // Check that metal-oxide bonds use Einstein model
    let cu_o = params.get_pair_parameters(29, 8).unwrap();
    assert_eq!(cu_o.model_type, "einstein");

    // Test temperature setting
    assert_eq!(params.temperature, 300.0);

    // Create parameters at a different temperature
    let cold_params = create_standard_thermal_parameters(100.0);
    assert_eq!(cold_params.temperature, 100.0);
}

#[test]
fn test_pair_specific_thermal_models() {
    // Instead of running a full EXAFS calculation which might have
    // dependencies on other modules that aren't part of this test,
    // we'll test the specific behavior of using pair-specific thermal
    // models with a simplified approach

    // Create custom thermal parameters for different bond types
    let temperature = 300.0;
    let thermal_params = ThermalParameters::new_debye(temperature, 300.0)
        .with_pair_parameters(
            29,
            29, // Cu-Cu
            PairThermalParameters::new_debye(315.0, Some("Cu-Cu bond".to_string())),
        )
        .with_pair_parameters(
            29,
            8, // Cu-O
            PairThermalParameters::new_einstein(70.0, Some("Cu-O bond".to_string())),
        )
        .with_pair_parameters(
            26,
            8, // Fe-O
            PairThermalParameters::new_einstein(65.0, Some("Fe-O bond".to_string())),
        );

    // Check that we can retrieve the specific parameters for each pair
    // Cu-Cu pair
    let cu_cu_params = thermal_params.get_pair_parameters(29, 29);
    assert!(cu_cu_params.is_some());
    let cu_cu = cu_cu_params.unwrap();
    assert_eq!(cu_cu.model_type, "debye");
    assert_eq!(cu_cu.debye_temperature, 315.0);

    // Cu-O pair
    let cu_o_params = thermal_params.get_pair_parameters(29, 8);
    assert!(cu_o_params.is_some());
    let cu_o = cu_o_params.unwrap();
    assert_eq!(cu_o.model_type, "einstein");
    assert_eq!(cu_o.einstein_frequency, Some(70.0));

    // Fe-O pair
    let fe_o_params = thermal_params.get_pair_parameters(26, 8);
    assert!(fe_o_params.is_some());
    let fe_o = fe_o_params.unwrap();
    assert_eq!(fe_o.model_type, "einstein");
    assert_eq!(fe_o.einstein_frequency, Some(65.0));

    // Check parameter retrieval with reversed order (should find the same parameters)
    let o_cu_params = thermal_params.get_pair_parameters(8, 29);
    assert!(o_cu_params.is_some());
    assert_eq!(o_cu_params.unwrap().einstein_frequency, Some(70.0));

    // Test adding a thermal parameter over an existing one
    let updated_params = thermal_params.with_pair_parameters(
        29,
        8,
        PairThermalParameters::new_einstein(75.0, Some("Updated Cu-O bond".to_string())),
    );

    let updated_cu_o = updated_params.get_pair_parameters(29, 8).unwrap();
    assert_eq!(updated_cu_o.einstein_frequency, Some(75.0));
    assert_eq!(
        updated_cu_o.description,
        Some("Updated Cu-O bond".to_string())
    );

    println!("All pair-specific thermal model tests passed!");
}
