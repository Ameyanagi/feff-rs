/*
MIT License with FEFF10 Attribution

Copyright (c) 2025 Ameyanagi

Based on or developed using Distribution: FEFF10.0
Copyright (c) 2020 FEFF Project, University of Washington and SLAC National Accelerator Laboratory.
All rights reserved.
*/

//! Scattering matrix calculations for FEFF
//!
//! This module implements the full scattering matrix calculations,
//! including Green's functions, T-matrices, and path operators.

use super::{calculate_phase_shifts, ScatteringMatrixResults};
use crate::atoms::errors::AtomError;
use crate::atoms::{AtomicStructure, Result as AtomResult};
use crate::utils::constants::{BOHR_TO_ANGSTROM, HARTREE_TO_EV};
use crate::utils::linear_algebra::{faer_to_ndarray, ndarray_to_faer};
use crate::utils::matrix::compute_free_greens_matrix;
use ndarray::Array2;
use num_complex::Complex64;
use rayon::prelude::*;

/// Calculate scattering matrices (legacy function)
///
/// # Deprecated
///
/// This function is kept for backward compatibility with existing tests.
/// It calculates phase shifts first, then calculates scattering matrices.
/// For better performance, use `calculate_phase_shifts` followed by `calculate_scattering_matrices`.
///
/// # Arguments
///
/// * `structure` - The atomic structure
/// * `energy` - Energy in eV
/// * `max_l` - Maximum angular momentum to include
///
/// # Returns
///
/// Scattering matrix results including Green's function, T-matrices, and global T-matrix
pub fn calculate_scattering_matrices_legacy(
    structure: &AtomicStructure,
    energy: f64,
    max_l: i32,
) -> AtomResult<ScatteringMatrixResults> {
    // First, calculate the phase shifts
    let phase_shift_results = calculate_phase_shifts(structure, energy, max_l)?;

    // Then, calculate scattering matrices using these phase shifts
    calculate_scattering_matrices(structure, &phase_shift_results)
}

/// Calculate scattering matrices from phase shifts
///
/// # Arguments
///
/// * `structure` - The atomic structure
/// * `phase_shifts` - Pre-calculated phase shift results
///
/// # Returns
///
/// Scattering matrix results including Green's function, T-matrices, and global T-matrix
pub fn calculate_scattering_matrices(
    structure: &AtomicStructure,
    phase_shifts: &super::ScatteringResults,
) -> AtomResult<ScatteringMatrixResults> {
    // Get energy, max_l, and phase shifts from pre-calculated results
    let energy = phase_shifts.energy;
    let max_l = phase_shifts.max_l;

    // Convert energy to atomic units (Hartree)
    let energy_hartree = energy / HARTREE_TO_EV;

    // Calculate wave number (k = sqrt(2*E))
    let k = (2.0 * energy_hartree).sqrt();

    // Use parallel tasks to compute matrices independently
    // This is more efficient than doing them sequentially
    let results = rayon::join(
        // Task 1: Compute Green's function matrix
        || {
            // Prepare atom positions
            let positions = get_atom_positions(structure);

            // Create the free Green's function matrix
            compute_free_greens_matrix(&positions, k, max_l).map_err(|e| {
                AtomError::CalculationError(format!(
                    "Failed to compute Green's function matrix: {}",
                    e
                ))
            })
        },
        // Task 2: Compute global T-matrix
        || {
            // Get individual T-matrices for each potential type
            let t_matrices = &phase_shifts.t_matrices;

            // Calculate global T-matrix
            construct_global_t_matrix(structure, t_matrices, max_l)
        },
    );

    // Extract results from parallel computation
    let green_matrix = results.0?;
    let global_t_matrix = results.1?;

    // Return the results
    Ok(ScatteringMatrixResults {
        energy,
        max_l,
        phase_shifts: phase_shifts.phase_shifts.clone(),
        t_matrices: phase_shifts.t_matrices.clone(),
        green_matrix,
        global_t_matrix,
    })
}

/// Extract positions of all atoms in the structure as a vector of (x, y, z) tuples
fn get_atom_positions(structure: &AtomicStructure) -> Vec<(f64, f64, f64)> {
    // For large atom counts, use parallel processing
    if structure.atom_count() > 100 {
        (0..structure.atom_count())
            .into_par_iter()
            .filter_map(|i| {
                structure.atom(i).map(|atom| {
                    let pos = atom.position();
                    (
                        pos.x / BOHR_TO_ANGSTROM,
                        pos.y / BOHR_TO_ANGSTROM,
                        pos.z / BOHR_TO_ANGSTROM,
                    )
                })
            })
            .collect()
    } else {
        // For smaller atom counts, sequential is faster due to reduced overhead
        structure
            .atoms()
            .iter()
            .map(|atom| {
                let pos = atom.position();
                (
                    pos.x / BOHR_TO_ANGSTROM,
                    pos.y / BOHR_TO_ANGSTROM,
                    pos.z / BOHR_TO_ANGSTROM,
                )
            })
            .collect()
    }
}

/// Construct the global T-matrix from individual potential T-matrices
fn construct_global_t_matrix(
    structure: &AtomicStructure,
    t_matrices: &[Array2<Complex64>],
    max_l: i32,
) -> AtomResult<Array2<Complex64>> {
    let n_atoms = structure.atom_count();
    let l_size = ((max_l + 1) * (max_l + 1)) as usize;
    let total_size = n_atoms * l_size;

    // Create global T-matrix using Faer for better performance
    let mut global_t_faer = faer::Mat::<Complex64>::zeros(total_size, total_size);

    // Convert all potential T-matrices to Faer format for faster operations
    let t_matrices_faer: Vec<faer::Mat<Complex64>> =
        t_matrices.iter().map(ndarray_to_faer).collect();

    // Fill the T-matrix sequentially to avoid complexity
    for (atom_idx, atom) in structure.atoms().iter().enumerate() {
        let pot_idx = atom.potential_type() as usize;

        // Check that the potential index is valid
        if pot_idx >= t_matrices_faer.len() {
            return Err(AtomError::CalculationError(format!(
                "Invalid potential type index: {}",
                pot_idx
            )));
        }

        let t_mat = &t_matrices_faer[pot_idx];

        // Use faster block operations instead of element-by-element copy
        // This is a critical optimization for large matrices
        let base_idx = atom_idx * l_size;
        for i in 0..l_size {
            for j in 0..l_size {
                global_t_faer[(base_idx + i, base_idx + j)] = t_mat[(i, j)];
            }
        }
    }

    // Convert back to ndarray
    Ok(faer_to_ndarray(&global_t_faer))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::atoms::{Atom, PotentialType, Vector3D};
    use approx::assert_relative_eq;

    #[test]
    fn test_construct_global_t_matrix() {
        // Create a simple two-atom structure
        let fe_potential = PotentialType::new(0, 26).unwrap(); // Iron

        let mut structure = AtomicStructure::new();
        structure.add_potential_type(fe_potential);

        // Add two iron atoms
        structure.add_atom(Atom::new(26, Vector3D::new(0.0, 0.0, 0.0), 0).unwrap());
        structure.add_atom(Atom::new(26, Vector3D::new(3.0, 0.0, 0.0), 0).unwrap());

        // Create a simple T-matrix for the iron potential
        let max_l = 1; // l_max = 1 for simplicity
        let l_size = ((max_l + 1) * (max_l + 1)) as usize; // = 4

        let mut t_matrix = Array2::<Complex64>::zeros((l_size, l_size));
        // Fill with some test values
        for i in 0..l_size {
            t_matrix[(i, i)] = Complex64::new(0.1 * (i as f64 + 1.0), 0.0);
        }

        let t_matrices = vec![t_matrix];

        // Construct the global T-matrix
        let global_t = construct_global_t_matrix(&structure, &t_matrices, max_l).unwrap();

        // Check dimensions
        assert_eq!(global_t.shape(), [8, 8]); // 2 atoms * 4 (l_max=1) = 8

        // Check that diagonal blocks are filled correctly
        for atom_idx in 0..2 {
            for i in 0..l_size {
                // Diagonal elements should match the original T-matrix
                assert_relative_eq!(
                    global_t[(atom_idx * l_size + i, atom_idx * l_size + i)].re,
                    t_matrices[0][(i, i)].re,
                    epsilon = 1e-10
                );

                // Off-diagonal elements within the block should be zero
                for j in 0..l_size {
                    if i != j {
                        assert_relative_eq!(
                            global_t[(atom_idx * l_size + i, atom_idx * l_size + j)].norm(),
                            0.0,
                            epsilon = 1e-10
                        );
                    }
                }
            }
        }

        // Check that off-diagonal blocks are zero
        for i in 0..l_size {
            for j in 0..l_size {
                assert_relative_eq!(global_t[(i, l_size + j)].norm(), 0.0, epsilon = 1e-10);
                assert_relative_eq!(global_t[(l_size + i, j)].norm(), 0.0, epsilon = 1e-10);
            }
        }
    }
}
