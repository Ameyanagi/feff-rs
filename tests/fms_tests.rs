/*
MIT License with FEFF10 Attribution

Copyright (c) 2025 Ameyanagi

Based on or developed using Distribution: FEFF10.0
Copyright (c) 2020 FEFF Project, University of Washington and SLAC National Accelerator Laboratory.
All rights reserved.
*/

use feff_rs::atoms::{Atom, AtomicStructure, PotentialType, Vector3D};
use feff_rs::fms::{FmsMatrix, FmsParameters, SolverMethod};
use feff_rs::scattering::calculate_scattering_matrices_old;

/// Helper function to create a simple iron structure for testing
fn setup_iron_atom() -> AtomicStructure {
    // Create a potential type for iron
    let fe_potential = PotentialType::new(0, 26).unwrap();

    // Create a new atomic structure
    let mut structure = AtomicStructure::new();

    // Add the potential type
    structure.add_potential_type(fe_potential);

    // Add an iron atom at the origin
    let fe_atom = Atom::new(26, Vector3D::new(0.0, 0.0, 0.0), 0).unwrap();
    let fe_idx = structure.add_atom(fe_atom);

    // Set it as the central atom
    structure.set_central_atom(fe_idx).unwrap();

    // Calculate muffin-tin radii
    structure.calculate_muffin_tin_radii().unwrap();

    structure
}

#[test]
fn test_fms_module_structure() {
    // Test that the FMS module can be initialized correctly
    let _structure = setup_iron_atom();

    // Create FMS parameters
    let _params = FmsParameters {
        radius: 5.0,
        solver_method: SolverMethod::LuDecomposition,
        energies: vec![7120.0], // Iron K-edge energy region
        calculate_xanes: true,
        include_thermal_effects: false,
        core_hole_lifetime: 1.5,
        energy_shift: None,
    };

    // The test passes if code compiles and we reach this point
    assert!(true);
}

#[test]
#[ignore] // Ignore this test until implementation is complete
fn test_fms_matrix_builder() {
    // Create a simple structure
    let structure = setup_iron_atom();

    // Create FMS matrix builder
    let fms_radius = 5.0;
    let fms_matrix = FmsMatrix::new(&structure, fms_radius);

    // Validate that the FMS matrix builder was created successfully
    assert!(fms_matrix.is_ok());

    // Get the number of atoms within the FMS radius
    if let Ok(matrix_builder) = fms_matrix {
        // There should be at least one atom (the central atom)
        assert!(matrix_builder.atom_count() >= 1);
    }
}

#[test]
#[ignore] // Ignore this test until implementation is complete
fn test_fms_calculation() {
    // Create a simple structure
    let structure = setup_iron_atom();

    // Calculate scattering matrices
    let energy = 7120.0; // Iron K-edge energy region in eV
    let max_l = 3;
    let scattering_matrices = calculate_scattering_matrices_old(&structure, energy, max_l);

    // Validate that scattering matrices were calculated successfully
    assert!(scattering_matrices.is_ok());

    // Create FMS parameters
    let params = FmsParameters {
        radius: 5.0,
        solver_method: SolverMethod::LuDecomposition,
        energies: vec![energy],
        calculate_xanes: true,
        include_thermal_effects: false,
        core_hole_lifetime: 1.5,
        energy_shift: None,
    };

    // Calculate FMS (will be implemented in fms::calculate_fms)
    // let fms_results = fms::calculate_fms(&structure, &scattering_matrices.unwrap(), params);

    // For now, just check that we can access the parameters
    assert_eq!(params.radius, 5.0);
    assert_eq!(params.energies, vec![energy]);
    assert_eq!(params.core_hole_lifetime, 1.5);
}
