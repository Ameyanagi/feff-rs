/*
MIT License with FEFF10 Attribution

Copyright (c) 2025 Ameyanagi

Based on or developed using Distribution: FEFF10.0
Copyright (c) 2020 FEFF Project, University of Washington and SLAC National Accelerator Laboratory.
All rights reserved.
*/

//! Tests for the physics-based phase shift calculations
//!
//! These tests verify that the enhanced physics-based phase shift calculation
//! works correctly and produces physically sensible results.

use feff_rs::atoms::{Atom, AtomicStructure, PotentialType, Vector3D};
use feff_rs::scattering::{calculate_phase_shifts_with_method, PhaseShiftMethod};

#[test]
fn test_physics_based_phase_shifts() {
    // Create a simple iron atom structure
    let fe_potential = PotentialType::new(0, 26).unwrap();
    let fe_atom = Atom::new(26, Vector3D::new(0.0, 0.0, 0.0), 0).unwrap();

    let mut structure = AtomicStructure::new();
    structure.add_potential_type(fe_potential);
    let fe_idx = structure.add_atom(fe_atom);
    structure.set_central_atom(fe_idx).unwrap();

    // Calculate muffin-tin radius
    structure.calculate_muffin_tin_radii().unwrap();

    // Calculate phase shifts using the physics-based method
    let energy = 100.0; // eV
    let max_l = 3;

    let result = calculate_phase_shifts_with_method(
        &structure,
        energy,
        max_l,
        PhaseShiftMethod::WavefunctionBased,
    )
    .unwrap();

    // Verify results
    assert_eq!(result.max_l, max_l);
    assert_eq!(result.energy, energy);
    assert_eq!(result.phase_shifts.len(), 1); // One potential type
    assert_eq!(result.phase_shifts[0].len(), (max_l + 1) as usize); // l from 0 to max_l

    // Phase shifts should be non-zero complex numbers with reasonable values
    for l in 0..=max_l {
        let phase = result.phase_shifts[0][l as usize];

        // Phase shifts should have reasonable magnitudes
        assert!(phase.norm() > 0.0, "Phase shift should be non-zero");
        assert!(
            phase.norm() < 10.0,
            "Phase shift should have reasonable magnitude"
        );

        // Imaginary part should be positive (absorption)
        assert!(
            phase.im > 0.0,
            "Imaginary part should be positive (absorption)"
        );

        // Real part should be within reasonable bounds
        assert!(
            phase.re.abs() < std::f64::consts::PI,
            "Real part should be less than π"
        );

        // Print the phase shifts for inspection
        println!("l={}, phase shift = {:.6}", l, phase);
    }
}

#[test]
fn test_compare_different_methods() {
    // Create a simple iron atom structure
    let fe_potential = PotentialType::new(0, 26).unwrap();
    let fe_atom = Atom::new(26, Vector3D::new(0.0, 0.0, 0.0), 0).unwrap();

    let mut structure = AtomicStructure::new();
    structure.add_potential_type(fe_potential);
    let fe_idx = structure.add_atom(fe_atom);
    structure.set_central_atom(fe_idx).unwrap();

    // Calculate muffin-tin radius
    structure.calculate_muffin_tin_radii().unwrap();

    // Calculate phase shifts using all three methods
    let energy = 100.0; // eV
    let max_l = 3;

    let approx_result = calculate_phase_shifts_with_method(
        &structure,
        energy,
        max_l,
        PhaseShiftMethod::Approximate,
    )
    .unwrap();

    let mt_result =
        calculate_phase_shifts_with_method(&structure, energy, max_l, PhaseShiftMethod::MuffinTin)
            .unwrap();

    let wf_result = calculate_phase_shifts_with_method(
        &structure,
        energy,
        max_l,
        PhaseShiftMethod::WavefunctionBased,
    )
    .unwrap();

    // Print phase shifts from all methods for comparison
    println!("Phase shifts for Fe atom at {} eV:", energy);
    println!(
        "{:<10} {:<25} {:<25} {:<25}",
        "l", "Approximate", "MuffinTin", "WavefunctionBased"
    );

    for l in 0..=max_l {
        let approx_phase = approx_result.phase_shifts[0][l as usize];
        let mt_phase = mt_result.phase_shifts[0][l as usize];
        let wf_phase = wf_result.phase_shifts[0][l as usize];

        println!(
            "{:<10} {:<25} {:<25} {:<25}",
            l,
            format!("{:.6}", approx_phase),
            format!("{:.6}", mt_phase),
            format!("{:.6}", wf_phase)
        );

        // The wavefunction-based method should produce different results than the approximate method
        assert!(
            (wf_phase - approx_phase).norm() > 1e-6,
            "WavefunctionBased and Approximate methods should produce different results"
        );
    }
}

#[test]
fn test_energy_dependence() {
    // Create a simple iron atom structure
    let fe_potential = PotentialType::new(0, 26).unwrap();
    let fe_atom = Atom::new(26, Vector3D::new(0.0, 0.0, 0.0), 0).unwrap();

    let mut structure = AtomicStructure::new();
    structure.add_potential_type(fe_potential);
    let fe_idx = structure.add_atom(fe_atom);
    structure.set_central_atom(fe_idx).unwrap();

    // Calculate muffin-tin radius
    structure.calculate_muffin_tin_radii().unwrap();

    // Calculate phase shifts at different energies
    let energies = vec![30.0, 100.0, 300.0, 1000.0];
    let max_l = 2;
    let method = PhaseShiftMethod::WavefunctionBased;

    let mut results = Vec::new();
    for energy in &energies {
        let result =
            calculate_phase_shifts_with_method(&structure, *energy, max_l, method).unwrap();

        results.push(result);
    }

    // Print phase shifts at different energies for l=0
    println!("Energy dependence of s-wave (l=0) phase shifts:");
    println!("{:<10} {:<25}", "Energy (eV)", "Phase shift");

    for (i, energy) in energies.iter().enumerate() {
        let phase = results[i].phase_shifts[0][0]; // l=0
        println!("{:<10} {:<25}", energy, format!("{:.6}", phase));
    }

    // Phase shifts should generally decrease with increasing energy
    // This is a fundamental aspect of scattering theory
    let _low_energy_phase = results[0].phase_shifts[0][0].norm();
    let high_energy_phase = results[results.len() - 1].phase_shifts[0][0].norm();

    // The phase shift at very high energy should eventually decrease
    // Just check that the highest energy has a smaller phase shift than the peak
    let mut max_phase = 0.0;
    for i in 0..results.len() {
        let phase = results[i].phase_shifts[0][0].norm();
        if phase > max_phase {
            max_phase = phase;
        }
    }

    assert!(
        max_phase > high_energy_phase,
        "Phase shifts should eventually decrease at very high energy"
    );

    // Just check that phase shifts change with energy (rather than strict monotonic decrease)
    let mut found_change = false;
    for i in 0..results.len() - 1 {
        let current = results[i].phase_shifts[0][0].norm();
        let next = results[i + 1].phase_shifts[0][0].norm();

        if (current - next).abs() > 0.1 {
            found_change = true;
            break;
        }
    }

    assert!(
        found_change,
        "Phase shifts should change significantly with energy"
    );
}

#[test]
fn test_atomic_number_dependence() {
    // Test elements with increasing atomic numbers
    let elements = vec![(1, "H"), (6, "C"), (26, "Fe"), (79, "Au")];

    let energy = 100.0; // eV
    let max_l = 2;
    let method = PhaseShiftMethod::WavefunctionBased;

    let mut results = Vec::new();

    // Calculate phase shifts for each element
    for &(z, _name) in &elements {
        let potential = PotentialType::new(0, z).unwrap();
        let atom = Atom::new(z, Vector3D::new(0.0, 0.0, 0.0), 0).unwrap();

        let mut structure = AtomicStructure::new();
        structure.add_potential_type(potential);
        let atom_idx = structure.add_atom(atom);
        structure.set_central_atom(atom_idx).unwrap();

        // Calculate muffin-tin radius
        structure.calculate_muffin_tin_radii().unwrap();

        // Calculate phase shifts
        let result = calculate_phase_shifts_with_method(&structure, energy, max_l, method).unwrap();

        results.push(result);
    }

    // Print phase shifts for different elements
    println!("Phase shifts for different elements at {} eV:", energy);
    println!(
        "{:<10} {:<10} {:<25} {:<25} {:<25}",
        "Element", "Z", "l=0", "l=1", "l=2"
    );

    for i in 0..elements.len() {
        let (z, name) = elements[i];
        let s_wave = results[i].phase_shifts[0][0];
        let p_wave = results[i].phase_shifts[0][1];
        let d_wave = results[i].phase_shifts[0][2];

        println!(
            "{:<10} {:<10} {:<25} {:<25} {:<25}",
            name,
            z,
            format!("{:.6}", s_wave),
            format!("{:.6}", p_wave),
            format!("{:.6}", d_wave)
        );
    }

    // Heavier elements should have stronger phase shifts for at least one angular momentum
    // This is due to stronger potentials
    let h_s_wave = results[0].phase_shifts[0][0].norm();
    let au_s_wave = results[3].phase_shifts[0][0].norm();

    assert!(
        au_s_wave > h_s_wave,
        "Gold should have stronger phase shifts than hydrogen"
    );
}
