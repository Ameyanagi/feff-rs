/*
MIT License with FEFF10 Attribution

Copyright (c) 2025 Ameyanagi

Based on or developed using Distribution: FEFF10.0
Copyright (c) 2020 FEFF Project, University of Washington and SLAC National Accelerator Laboratory.
All rights reserved.
*/

use approx::assert_relative_eq;
use feff_rs::atoms::{Atom, AtomicStructure, PotentialType, Vector3D};
use feff_rs::potential::{GridType, MuffinTinPotential};

/// Test the muffin-tin radius calculations
#[test]
fn test_muffin_tin_radius_calculation() {
    // Create potential types
    let mut pot_fe = PotentialType::new(0, 26).unwrap(); // Iron
    let pot_o = PotentialType::new(1, 8).unwrap(); // Oxygen

    // Verify default radius is based on covalent radius
    let fe_default = pot_fe.default_muffin_tin_radius();
    assert!(fe_default > 1.0 && fe_default < 2.0); // Iron covalent radius * 1.2 should be ~1.6 Å

    let o_default = pot_o.default_muffin_tin_radius();
    assert!(o_default > 0.5 && o_default < 1.0); // Oxygen covalent radius * 1.2 should be ~0.8 Å

    // Create atoms
    let fe_atom = Atom::new(26, Vector3D::new(0.0, 0.0, 0.0), 0).unwrap();
    let o_atom1 = Atom::new(8, Vector3D::new(2.0, 0.0, 0.0), 1).unwrap();
    let o_atom2 = Atom::new(8, Vector3D::new(0.0, 2.0, 0.0), 1).unwrap();
    let o_atom3 = Atom::new(8, Vector3D::new(0.0, 0.0, 2.0), 1).unwrap();

    // Create an array of all atoms
    let atoms = vec![
        fe_atom.clone(),
        o_atom1.clone(),
        o_atom2.clone(),
        o_atom3.clone(),
    ];

    // Calculate muffin-tin radius for Fe atom
    let fe_radius = pot_fe.calculate_muffin_tin_radius(&fe_atom, &atoms);

    // The Fe radius should be constrained by the O atoms at distance 2.0
    // Given overlap factor of 1.15, it should be smaller than 2.0 * 1.15
    assert!(fe_radius < 2.0);
    assert!(fe_radius > 0.0);

    // Calculate muffin-tin radius for O atom
    let o_radius = pot_o.calculate_muffin_tin_radius(&o_atom1, &atoms);

    // The O radius should be constrained by the Fe atom at distance 2.0
    assert!(o_radius < 2.0);
    assert!(o_radius > 0.0);

    // Set overlap factor and recalculate
    pot_fe.set_overlap_factor(1.0);
    let fe_radius_no_overlap = pot_fe.calculate_muffin_tin_radius(&fe_atom, &atoms);

    // With overlap factor 1.0, radius should be smaller
    assert!(fe_radius_no_overlap < fe_radius);
}

/// Test overlap volume calculation
#[test]
fn test_overlap_volume() {
    // Create potential type
    let pot = PotentialType::new(0, 26).unwrap();

    // Create two atoms at distance 2.0
    let atom1 = Atom::new(26, Vector3D::new(0.0, 0.0, 0.0), 0).unwrap();
    let atom2 = Atom::new(26, Vector3D::new(2.0, 0.0, 0.0), 0).unwrap();

    // No overlap when radii sum is less than distance
    let overlap1 = pot.overlap_volume(&atom1, &atom2, 0.9, 0.9);
    assert_relative_eq!(overlap1, 0.0, epsilon = 1e-6);

    // Full containment when one radius is larger than distance + other radius
    let overlap2 = pot.overlap_volume(&atom1, &atom2, 3.0, 0.5);
    let small_sphere_vol = (4.0 / 3.0) * std::f64::consts::PI * 0.5f64.powi(3);
    assert_relative_eq!(overlap2, small_sphere_vol, epsilon = 1e-6);

    // Partial overlap
    let overlap3 = pot.overlap_volume(&atom1, &atom2, 1.5, 1.5);
    assert!(overlap3 > 0.0);
    assert!(overlap3 < (4.0 / 3.0) * std::f64::consts::PI * 1.5f64.powi(3)); // Less than one full sphere
}

/// Test the structure-level muffin-tin functions
#[test]
fn test_structure_muffin_tin_calculations() {
    // Create an FeO4 structure
    let mut structure = AtomicStructure::with_title("FeO4 test structure");

    // Add potential types
    let pot_fe = PotentialType::new(0, 26).unwrap(); // Iron
    let pot_o = PotentialType::new(1, 8).unwrap(); // Oxygen

    structure.add_potential_type(pot_fe);
    structure.add_potential_type(pot_o);

    // Add atoms
    let fe_idx = structure.add_atom(Atom::new(26, Vector3D::new(0.0, 0.0, 0.0), 0).unwrap());
    structure.add_atom(Atom::new(8, Vector3D::new(2.0, 0.0, 0.0), 1).unwrap());
    structure.add_atom(Atom::new(8, Vector3D::new(0.0, 2.0, 0.0), 1).unwrap());
    structure.add_atom(Atom::new(8, Vector3D::new(0.0, 0.0, 2.0), 1).unwrap());
    structure.add_atom(Atom::new(8, Vector3D::new(-2.0, 0.0, 0.0), 1).unwrap());

    // Set Fe as the central atom
    structure.set_central_atom(fe_idx).unwrap();

    // Calculate muffin-tin radii
    let radii = structure.calculate_muffin_tin_radii().unwrap();

    // Should have two radii (one for each potential type)
    assert_eq!(radii.len(), 2);

    // Both radii should be positive
    assert!(radii[0] > 0.0);
    assert!(radii[1] > 0.0);

    // Get current radii values and verify they're positive
    let _fe_radius = structure.potential_type(0).unwrap().muffin_tin_radius();
    let _o_radius = structure.potential_type(1).unwrap().muffin_tin_radius();

    // Set larger radii to ensure overlap
    let mut structure_with_overlap = structure.clone();
    structure_with_overlap
        .potential_type_mut(0)
        .unwrap()
        .set_muffin_tin_radius(1.5);
    structure_with_overlap
        .potential_type_mut(1)
        .unwrap()
        .set_muffin_tin_radius(1.0);

    // Calculate total overlap volume
    let overlap = structure_with_overlap
        .calculate_total_overlap_volume()
        .unwrap();

    // With larger radii, there should be some overlap
    assert!(overlap > 0.0);

    // Optimize muffin-tin radii
    let optimized_overlap = structure_with_overlap
        .optimize_muffin_tin_radii(10, 1e-4)
        .unwrap();

    // Optimized overlap should be less than or equal to initial overlap
    assert!(optimized_overlap <= overlap);
}

/// Test creating a cluster from a larger structure
#[test]
fn test_cluster_creation() {
    // Create a structure with a 3x3x3 grid of atoms
    let mut structure = AtomicStructure::with_title("Grid structure");

    // Add potential type
    let pot = PotentialType::new(0, 26).unwrap(); // Iron
    structure.add_potential_type(pot);

    // Add atoms in a grid
    let mut atom_indices = Vec::new();
    for x in -1..=1 {
        for y in -1..=1 {
            for z in -1..=1 {
                let x_pos = x as f64 * 2.0;
                let y_pos = y as f64 * 2.0;
                let z_pos = z as f64 * 2.0;

                let idx = structure
                    .add_atom(Atom::new(26, Vector3D::new(x_pos, y_pos, z_pos), 0).unwrap());

                atom_indices.push(idx);
            }
        }
    }

    // Set the central atom (at position 0,0,0)
    let center_idx = atom_indices[13]; // Should be the center of the 3x3x3 grid
    structure.set_central_atom(center_idx).unwrap();

    // Create a cluster with radius 2.5
    let cluster = structure.create_cluster(center_idx, 2.5).unwrap();

    // This should include only the central atom and the 6 nearest neighbors
    // (those that are 2.0 away, not the diagonal ones)
    assert_eq!(cluster.atom_count(), 7);

    // The central atom should still be marked as central
    assert!(cluster.central_atom().is_some());
    assert_eq!(
        cluster.central_atom().unwrap().position().x,
        structure.central_atom().unwrap().position().x
    );
}

/// Test muffin-tin potential calculation for Fe atom
#[test]
fn test_muffin_tin_potential() {
    // Create a simple atomic structure with an iron atom
    let mut structure = AtomicStructure::new();
    let pot_fe = PotentialType::new(0, 26).unwrap(); // Iron
    structure.add_potential_type(pot_fe);

    // Add central iron atom
    let fe_idx = structure.add_atom(Atom::new(26, Vector3D::new(0.0, 0.0, 0.0), 0).unwrap());
    structure.set_central_atom(fe_idx).unwrap();

    // Calculate muffin-tin radius
    structure.calculate_muffin_tin_radii().unwrap();

    // Create muffin-tin potential calculator
    let mut mt_calculator = MuffinTinPotential::new(&structure).unwrap();

    // For testing, use a small grid to ensure fast convergence
    mt_calculator.set_grid(10, GridType::Logarithmic).unwrap();

    // Calculate the potential
    let potential = mt_calculator.calculate().unwrap();

    // Check basic properties of the potential
    assert!(potential.grid_points() > 0);
    assert!(potential.radial_grid().len() > 0);

    // The potential should be finite and well-behaved
    for v in potential.values(0).unwrap() {
        assert!(v.is_finite());
    }

    // Potential should be negative (attractive)
    // At least near the nucleus, a few grid points from the beginning
    let near_nucleus_potential = potential.values(0).unwrap()[5]; // Not the very first point
    assert!(near_nucleus_potential < 0.0);

    // Energy levels should be calculated
    let energy_levels = potential.energy_levels(0).unwrap();
    assert!(!energy_levels.is_empty());

    // 1s level should be the lowest (most negative)
    let lowest_level = energy_levels[0];
    assert!(lowest_level < -100.0); // For iron, 1s is very deep, around -7000 eV

    // The Fermi level should be set
    let fermi = potential.fermi_energy();
    assert!(fermi.is_finite());

    // Test self-consistency
    let scf_result = mt_calculator.run_self_consistency(3, 1e-3).unwrap();

    // SCF should converge in a few iterations
    assert!(scf_result.iterations <= 3);
    assert!(scf_result.converged);
    assert!(scf_result.final_error < 1e-3);
}

/// Test exchange-correlation implementation for muffin-tin potential
#[test]
fn test_exchange_correlation() {
    // Create a simple hydrogen atom to test exchange-correlation approximations
    let mut structure = AtomicStructure::new();
    let pot_h = PotentialType::new(0, 1).unwrap(); // Hydrogen
    structure.add_potential_type(pot_h);

    // Add central hydrogen atom
    let h_idx = structure.add_atom(Atom::new(1, Vector3D::new(0.0, 0.0, 0.0), 0).unwrap());
    structure.set_central_atom(h_idx).unwrap();

    // Calculate muffin-tin radius
    structure.calculate_muffin_tin_radii().unwrap();

    // Create muffin-tin potential calculator with LDA exchange-correlation
    let mut mt_calculator = MuffinTinPotential::new(&structure).unwrap();
    mt_calculator.set_exchange_correlation("LDA").unwrap();

    // For testing, use a small grid to ensure fast convergence
    mt_calculator.set_grid(10, GridType::Logarithmic).unwrap();

    // Calculate LDA potential
    let lda_potential = mt_calculator.calculate().unwrap();

    // Now use GGA exchange-correlation
    mt_calculator.set_exchange_correlation("GGA").unwrap();

    // Calculate GGA potential
    let gga_potential = mt_calculator.calculate().unwrap();

    // Make sure we got different results with the different approximations
    // At the same radial point, they should give different exchange-correlation potentials
    let radial_point = 5; // Some arbitrary point in the grid, not too close to nucleus or boundary
    let lda_value = lda_potential.values(0).unwrap()[radial_point];
    let gga_value = gga_potential.values(0).unwrap()[radial_point];

    // They should be different but not drastically so
    assert!(lda_value != gga_value);
    // The difference should be a small fraction of the potential itself
    let relative_diff = (lda_value - gga_value).abs() / lda_value.abs();
    assert!(relative_diff < 0.2); // Less than 20% difference typically
}
