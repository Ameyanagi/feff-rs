/*
MIT License with FEFF10 Attribution

Copyright (c) 2025 Ameyanagi

Based on or developed using Distribution: FEFF10.0
Copyright (c) 2020 FEFF Project, University of Washington and SLAC National Accelerator Laboratory.
All rights reserved.
*/

//! Tests for the XANES calculator module

use feff_rs::atoms::{Atom, AtomicStructure, Vector3D};
use feff_rs::fms::XanesCalculator;
use ndarray::Array2;
use num_complex::Complex64;

#[test]
fn test_xanes_calculator_creation() {
    // Create a simple structure with a single atom
    let mut structure = AtomicStructure::new();

    // Add a central atom (Fe)
    let central_atom = Atom::new(26, Vector3D::new(0.0, 0.0, 0.0), 0).unwrap();
    let central_idx = structure.add_atom(central_atom);
    structure.set_central_atom(central_idx).unwrap();

    // Create a XANES calculator
    let core_hole_lifetime = 1.5; // eV
    let calculator = XanesCalculator::new(&structure, core_hole_lifetime);

    // Verify properties
    assert_eq!(calculator.get_core_hole_lifetime(), 1.5);
    assert_eq!(calculator.get_polarization(), [1.0, 0.0, 0.0]); // Default x-polarization
    assert!(!calculator.get_quadrupole_included()); // Quadrupole not included by default
}

#[test]
fn test_dipole_matrix_element() {
    // Create a simple structure
    let mut structure = AtomicStructure::new();
    let central_atom = Atom::new(26, Vector3D::new(0.0, 0.0, 0.0), 0).unwrap();
    let central_idx = structure.add_atom(central_atom);
    structure.set_central_atom(central_idx).unwrap();

    // Create a XANES calculator and configure it
    let mut x_pol = XanesCalculator::new(&structure, 1.0);
    x_pol.set_polarization([1.0, 0.0, 0.0]);
    x_pol.set_initial_state(1, 0); // K-edge (1s)

    let mut y_pol = XanesCalculator::new(&structure, 1.0);
    y_pol.set_polarization([0.0, 1.0, 0.0]);
    y_pol.set_initial_state(1, 0); // K-edge (1s)

    let mut z_pol = XanesCalculator::new(&structure, 1.0);
    z_pol.set_polarization([0.0, 0.0, 1.0]);
    z_pol.set_initial_state(1, 0); // K-edge (1s)

    // Test s -> p transitions with different polarizations
    let dipole_x_m1 = x_pol.calculate_dipole_matrix_element(1, -1);
    let dipole_x_0 = x_pol.calculate_dipole_matrix_element(1, 0);
    let dipole_x_1 = x_pol.calculate_dipole_matrix_element(1, 1);

    let dipole_z_m1 = z_pol.calculate_dipole_matrix_element(1, -1);
    let dipole_z_0 = z_pol.calculate_dipole_matrix_element(1, 0);
    let dipole_z_1 = z_pol.calculate_dipole_matrix_element(1, 1);

    // Verify that x-polarization couples to m=±1 states, z-polarization to m=0
    assert!(dipole_x_m1.norm() > 0.0);
    assert!(dipole_x_1.norm() > 0.0);
    assert_eq!(dipole_x_0, Complex64::new(0.0, 0.0)); // No z component

    assert_eq!(dipole_z_m1, Complex64::new(0.0, 0.0)); // No x,y components
    assert!(dipole_z_0.norm() > 0.0);
    assert_eq!(dipole_z_1, Complex64::new(0.0, 0.0)); // No x,y components

    // Verify selection rules - only l±1 transitions allowed for dipole
    assert_eq!(
        x_pol.calculate_dipole_matrix_element(0, 0),
        Complex64::new(0.0, 0.0)
    ); // Forbidden: s -> s
    assert_eq!(
        x_pol.calculate_dipole_matrix_element(2, 0),
        Complex64::new(0.0, 0.0)
    ); // Forbidden direct: s -> d
}

#[test]
fn test_quadrupole_matrix_element() {
    // Create a simple structure
    let mut structure = AtomicStructure::new();
    let central_atom = Atom::new(26, Vector3D::new(0.0, 0.0, 0.0), 0).unwrap();
    let central_idx = structure.add_atom(central_atom);
    structure.set_central_atom(central_idx).unwrap();

    // Create a XANES calculator and configure it
    let mut calculator = XanesCalculator::new(&structure, 1.0);
    calculator.set_initial_state(1, 0); // K-edge (1s)
    calculator.set_include_quadrupole(true);

    // Verify selection rules - only l±2 transitions allowed for quadrupole
    assert_eq!(
        calculator.calculate_quadrupole_matrix_element(1, 0),
        Complex64::new(0.0, 0.0)
    ); // Forbidden: s -> p
    assert!(calculator.calculate_quadrupole_matrix_element(2, 0).norm() > 0.0); // Allowed: s -> d
    assert_eq!(
        calculator.calculate_quadrupole_matrix_element(3, 0),
        Complex64::new(0.0, 0.0)
    ); // Forbidden: s -> f
}

#[test]
fn test_xanes_calculation() {
    // Create a simple structure
    let mut structure = AtomicStructure::new();

    // Add a central atom (Fe)
    let central_atom = Atom::new(26, Vector3D::new(0.0, 0.0, 0.0), 0).unwrap();
    let central_idx = structure.add_atom(central_atom);
    structure.set_central_atom(central_idx).unwrap();

    // Add some surrounding atoms to form an octahedral environment
    let octahedral_distance = 2.0; // Angstroms

    // Add 6 oxygen atoms in octahedral coordination
    structure.add_atom(Atom::new(8, Vector3D::new(octahedral_distance, 0.0, 0.0), 1).unwrap());
    structure.add_atom(Atom::new(8, Vector3D::new(-octahedral_distance, 0.0, 0.0), 1).unwrap());
    structure.add_atom(Atom::new(8, Vector3D::new(0.0, octahedral_distance, 0.0), 1).unwrap());
    structure.add_atom(Atom::new(8, Vector3D::new(0.0, -octahedral_distance, 0.0), 1).unwrap());
    structure.add_atom(Atom::new(8, Vector3D::new(0.0, 0.0, octahedral_distance), 1).unwrap());
    structure.add_atom(Atom::new(8, Vector3D::new(0.0, 0.0, -octahedral_distance), 1).unwrap());

    // Create a mock path operator for testing
    // In a real calculation, this would come from the FMS solver
    let l_max = 3;
    let l_size = (l_max + 1) * (l_max + 1);
    let atom_count = structure.atom_count();
    let matrix_size = atom_count * l_size;

    let mut path_operator = Array2::<Complex64>::zeros((matrix_size, matrix_size));

    // Fill the path operator with some realistic values
    // For the test, we'll use a simplified diagonal path operator
    for i in 0..matrix_size {
        path_operator[(i, i)] = Complex64::new(0.0, 1.0); // Diagonal elements
    }

    // Create the XANES calculator
    let mut calculator = XanesCalculator::new(&structure, 1.0);
    calculator.set_fermi_energy(0.0);
    calculator.set_energy_shift(0.0);
    calculator.set_max_l(l_max as usize);

    // Calculate XANES at different energies
    let edge_energy = calculator.calculate_edge_energy(26).unwrap(); // Edge energy for Fe

    // Calculate below edge
    let below_edge = calculator
        .calculate_xanes(edge_energy - 10.0, &path_operator)
        .unwrap();

    // Calculate at edge
    let at_edge = calculator
        .calculate_xanes(edge_energy, &path_operator)
        .unwrap();

    // Calculate above edge
    let above_edge = calculator
        .calculate_xanes(edge_energy + 10.0, &path_operator)
        .unwrap();

    // Verify that there's no absorption below edge
    assert_eq!(below_edge, 0.0);

    // Verify that absorption is larger at and above edge
    assert!(at_edge > 0.0);
    assert!(above_edge > 0.0);

    // Test with and without quadrupole
    let mut with_quad_calc = XanesCalculator::new(&structure, 1.0);
    with_quad_calc.set_fermi_energy(0.0);
    with_quad_calc.set_energy_shift(0.0);
    with_quad_calc.set_max_l(l_max as usize);
    with_quad_calc.set_include_quadrupole(true);
    let with_quadrupole = with_quad_calc
        .calculate_xanes(edge_energy + 20.0, &path_operator)
        .unwrap();

    let mut without_quad_calc = XanesCalculator::new(&structure, 1.0);
    without_quad_calc.set_fermi_energy(0.0);
    without_quad_calc.set_energy_shift(0.0);
    without_quad_calc.set_max_l(l_max as usize);
    without_quad_calc.set_include_quadrupole(false);
    let without_quadrupole = without_quad_calc
        .calculate_xanes(edge_energy + 20.0, &path_operator)
        .unwrap();

    // Verify that including quadrupole transitions doesn't decrease absorption
    // This is a very simple test - in a real system, the difference would depend on the structure
    assert!(with_quadrupole >= without_quadrupole);
}
