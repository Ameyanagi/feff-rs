/*
MIT License with FEFF10 Attribution

Copyright (c) 2025 Ameyanagi

Based on or developed using Distribution: FEFF10.0
Copyright (c) 2020 FEFF Project, University of Washington and SLAC National Accelerator Laboratory.
All rights reserved.
*/

use approx::assert_relative_eq;
use feff_rs::atoms::{Atom, AtomicStructure, PotentialType, Vector3D};
use feff_rs::scattering::calculate_scattering_matrices_old as calculate_scattering_matrices;
use ndarray::Array2;
use num_complex::Complex64;

/// Test basic scattering matrix calculation for a single atom
#[test]
fn test_single_atom_scattering() {
    // Create a single iron atom structure
    let fe_potential = PotentialType::new(0, 26).unwrap(); // Iron
    let fe_atom = Atom::new(26, Vector3D::new(0.0, 0.0, 0.0), 0).unwrap();

    let mut structure = AtomicStructure::new();
    structure.add_potential_type(fe_potential);
    let fe_idx = structure.add_atom(fe_atom);
    structure.set_central_atom(fe_idx).unwrap();

    // Calculate muffin-tin radii
    structure.calculate_muffin_tin_radii().unwrap();

    // Calculate scattering matrices at a specific energy
    let energy = 100.0; // eV
    let max_l = 2;

    let result = calculate_scattering_matrices(&structure, energy, max_l).unwrap();

    // Check basic properties of the result
    assert_eq!(result.energy, energy);
    assert_eq!(result.max_l, max_l);

    // For a single atom, G matrix should be zero (no propagation between atoms)
    assert_eq!(result.green_matrix.shape()[0], 9); // (l_max+1)^2 = 9

    // The T matrix should be diagonal with elements determined by the phase shifts
    assert_eq!(result.t_matrices.len(), 1); // One potential type
    assert_eq!(result.t_matrices[0].shape(), [9, 9]); // (l_max+1)^2 x (l_max+1)^2

    // Check that the T matrix is diagonal
    let t_matrix = &result.t_matrices[0];
    for i in 0..9 {
        for j in 0..9 {
            if i == j {
                // Diagonal elements should be non-zero
                assert!(t_matrix[(i, j)].norm() > 0.0);
            } else {
                // Off-diagonal elements should be zero
                assert_relative_eq!(t_matrix[(i, j)].norm(), 0.0, epsilon = 1e-10);
            }
        }
    }
}

/// Test scattering matrices with multiple atoms
#[test]
fn test_multiple_atom_scattering() {
    // Create FeO structure with iron at center and 6 oxygens in octahedral positions
    let fe_potential = PotentialType::new(0, 26).unwrap(); // Iron
    let o_potential = PotentialType::new(1, 8).unwrap(); // Oxygen

    let mut structure = AtomicStructure::new();
    structure.add_potential_type(fe_potential);
    structure.add_potential_type(o_potential);

    // Add Fe atom at center
    let fe_idx = structure.add_atom(Atom::new(26, Vector3D::new(0.0, 0.0, 0.0), 0).unwrap());

    // Add O atoms in octahedral positions (along x, y, z axes)
    let distance = 2.0; // 2 Angstrom Fe-O distance
    structure.add_atom(Atom::new(8, Vector3D::new(distance, 0.0, 0.0), 1).unwrap());
    structure.add_atom(Atom::new(8, Vector3D::new(-distance, 0.0, 0.0), 1).unwrap());
    structure.add_atom(Atom::new(8, Vector3D::new(0.0, distance, 0.0), 1).unwrap());
    structure.add_atom(Atom::new(8, Vector3D::new(0.0, -distance, 0.0), 1).unwrap());
    structure.add_atom(Atom::new(8, Vector3D::new(0.0, 0.0, distance), 1).unwrap());
    structure.add_atom(Atom::new(8, Vector3D::new(0.0, 0.0, -distance), 1).unwrap());

    structure.set_central_atom(fe_idx).unwrap();

    // Calculate muffin-tin radii
    structure.calculate_muffin_tin_radii().unwrap();

    // Calculate scattering matrices
    let energy = 100.0; // eV
    let max_l = 2;

    let result = calculate_scattering_matrices(&structure, energy, max_l).unwrap();

    // Check basic properties
    assert_eq!(result.energy, energy);
    assert_eq!(result.max_l, max_l);
    assert_eq!(result.t_matrices.len(), 2); // Two potential types

    // Calculate the size of the complete matrices
    let l_size = (max_l + 1) * (max_l + 1); // (l_max+1)^2
    let n_atoms = structure.atom_count();
    let total_size = l_size * n_atoms as i32;

    // Green's function matrix should have proper dimensions
    assert_eq!(
        result.green_matrix.shape(),
        [total_size as usize, total_size as usize]
    );

    // Each T matrix should have proper dimensions
    assert_eq!(
        result.t_matrices[0].shape(),
        [l_size as usize, l_size as usize]
    );
    assert_eq!(
        result.t_matrices[1].shape(),
        [l_size as usize, l_size as usize]
    );

    // Global T matrix should have proper dimensions
    assert_eq!(
        result.global_t_matrix.shape(),
        [total_size as usize, total_size as usize]
    );

    // Check that the Green's function matrix has proper structure
    // It should have zero diagonal blocks (no propagation within the same atom)
    let g_matrix = &result.green_matrix;
    for atom_idx in 0..n_atoms {
        for l1 in 0..l_size {
            for l2 in 0..l_size {
                let i = (atom_idx as i32 * l_size + l1) as usize;
                let j = (atom_idx as i32 * l_size + l2) as usize;

                // Diagonal blocks should be zero
                assert_relative_eq!(g_matrix[(i, j)].norm(), 0.0, epsilon = 1e-10);
            }
        }
    }

    // Check non-zero elements in off-diagonal blocks
    let mut has_nonzero = false;
    for atom_i in 0..n_atoms {
        for atom_j in 0..n_atoms {
            if atom_i != atom_j {
                for l1 in 0..l_size {
                    for l2 in 0..l_size {
                        let i = (atom_i as i32 * l_size + l1) as usize;
                        let j = (atom_j as i32 * l_size + l2) as usize;

                        if g_matrix[(i, j)].norm() > 1e-10 {
                            has_nonzero = true;
                        }
                    }
                }
            }
        }
    }

    // There should be some non-zero elements in the Green's function matrix
    assert!(has_nonzero);
}

/// Test path operator calculation
#[test]
fn test_path_operator() {
    // Create a simple two-atom structure
    let fe_potential = PotentialType::new(0, 26).unwrap(); // Iron

    let mut structure = AtomicStructure::new();
    structure.add_potential_type(fe_potential);

    // Add two iron atoms separated by 3 Angstroms
    let fe1_idx = structure.add_atom(Atom::new(26, Vector3D::new(0.0, 0.0, 0.0), 0).unwrap());
    structure.add_atom(Atom::new(26, Vector3D::new(3.0, 0.0, 0.0), 0).unwrap());

    structure.set_central_atom(fe1_idx).unwrap();

    // Calculate muffin-tin radii
    structure.calculate_muffin_tin_radii().unwrap();

    // Calculate scattering matrices
    let energy = 100.0; // eV
    let max_l = 2;

    let result = calculate_scattering_matrices(&structure, energy, max_l).unwrap();

    // Calculate path operator matrix
    let path_operator = result.calculate_path_operator().unwrap();

    // Check dimensions
    let l_size = (max_l + 1) * (max_l + 1); // (l_max+1)^2
    let total_size = l_size * 2; // 2 atoms
    assert_eq!(
        path_operator.shape(),
        [total_size as usize, total_size as usize]
    );

    // Path operator should be different from the global T matrix
    // because it includes multiple scattering effects
    let mut has_difference = false;
    for i in 0..total_size as usize {
        for j in 0..total_size as usize {
            if (path_operator[(i, j)] - result.global_t_matrix[(i, j)]).norm() > 1e-10 {
                has_difference = true;
            }
        }
    }

    assert!(has_difference);
}

/// Test energy dependence of scattering matrices
#[test]
fn test_energy_dependence() {
    // Create a simple iron atom structure
    let fe_potential = PotentialType::new(0, 26).unwrap();
    let fe_atom = Atom::new(26, Vector3D::new(0.0, 0.0, 0.0), 0).unwrap();

    let mut structure = AtomicStructure::new();
    structure.add_potential_type(fe_potential);
    let fe_idx = structure.add_atom(fe_atom);
    structure.set_central_atom(fe_idx).unwrap();

    // Calculate muffin-tin radii
    structure.calculate_muffin_tin_radii().unwrap();

    // Calculate scattering matrices at two different energies
    let energy1 = 50.0; // eV
    let energy2 = 500.0; // eV
    let max_l = 2;

    let result1 = calculate_scattering_matrices(&structure, energy1, max_l).unwrap();
    let result2 = calculate_scattering_matrices(&structure, energy2, max_l).unwrap();

    // Compare the T matrices at different energies
    // The high-energy T matrix should have smaller norm
    let t1_norm = matrix_norm(&result1.t_matrices[0]);
    let t2_norm = matrix_norm(&result2.t_matrices[0]);

    assert!(t2_norm < t1_norm);
}

/// Test l_max convergence
#[test]
fn test_lmax_convergence() {
    // Create a simple iron atom structure
    let fe_potential = PotentialType::new(0, 26).unwrap();
    let fe_atom = Atom::new(26, Vector3D::new(0.0, 0.0, 0.0), 0).unwrap();

    let mut structure = AtomicStructure::new();
    structure.add_potential_type(fe_potential);
    let fe_idx = structure.add_atom(fe_atom);
    structure.set_central_atom(fe_idx).unwrap();

    // Calculate muffin-tin radii
    structure.calculate_muffin_tin_radii().unwrap();

    // Calculate scattering matrices with different l_max values
    let energy = 100.0; // eV
    let l_max1 = 2;
    let l_max2 = 3;

    let result1 = calculate_scattering_matrices(&structure, energy, l_max1).unwrap();
    let result2 = calculate_scattering_matrices(&structure, energy, l_max2).unwrap();

    // Check that higher l_max includes all the same elements as lower l_max
    // plus additional angular momentum channels
    let t1 = &result1.t_matrices[0];
    let t2 = &result2.t_matrices[0];

    // The l_max2 matrix should be larger
    assert!(t2.shape()[0] > t1.shape()[0]);

    // The top-left block of the l_max2 matrix should match the l_max1 matrix
    for i in 0..t1.shape()[0] {
        for j in 0..t1.shape()[1] {
            assert_relative_eq!(t2[(i, j)].re, t1[(i, j)].re, epsilon = 1e-10);
            assert_relative_eq!(t2[(i, j)].im, t1[(i, j)].im, epsilon = 1e-10);
        }
    }
}

/// Calculate the Frobenius norm of a matrix
fn matrix_norm(matrix: &Array2<Complex64>) -> f64 {
    let mut sum_sq = 0.0;

    for i in 0..matrix.shape()[0] {
        for j in 0..matrix.shape()[1] {
            sum_sq += matrix[(i, j)].norm_sqr();
        }
    }

    sum_sq.sqrt()
}
