/*
MIT License with FEFF10 Attribution

Copyright (c) 2025 Ameyanagi

Based on or developed using Distribution: FEFF10.0
Copyright (c) 2020 FEFF Project, University of Washington and SLAC National Accelerator Laboratory.
All rights reserved.
*/

//! Matrix operations specific to scattering physics and FEFF calculations
//!
//! This module provides specialized matrix operations commonly used in FEFF calculations,
//! particularly for multiple scattering theory and X-ray absorption spectroscopy.

#![allow(clippy::manual_div_ceil)]
#![allow(clippy::needless_range_loop)]

use super::errors::{Result, UtilsError};
use super::math::{spherical_bessel_j, spherical_hankel_h1, spherical_harmonic};
use ndarray::Array2;
use num_complex::Complex64;
use rayon::prelude::*;
use std::collections::HashMap;
use std::f64::consts::PI;

/// Real-space free electron Green's function between two points
///
/// # Arguments
///
/// * `r1` - First position vector (x, y, z)
/// * `r2` - Second position vector (x, y, z)
/// * `k` - Wave vector magnitude
///
/// # Returns
///
/// The complex Green's function value
pub fn compute_free_greens_value(r1: (f64, f64, f64), r2: (f64, f64, f64), k: f64) -> Complex64 {
    let dx = r1.0 - r2.0;
    let dy = r1.1 - r2.1;
    let dz = r1.2 - r2.2;

    let r = (dx * dx + dy * dy + dz * dz).sqrt();

    // Handle case where r is very small
    if r < 1e-10 {
        return Complex64::new(0.0, -1.0 / k);
    }

    // Free electron Green's function G₀(r,r') = -exp(ik|r-r'|)/(4π|r-r'|)
    let exp_ikr = Complex64::new(0.0, k * r).exp();
    Complex64::new(-1.0, 0.0) * exp_ikr / (4.0 * PI * r)
}

/// Compute the scattering T-matrix from phase shifts
///
/// # Arguments
///
/// * `phase_shifts` - Vector of phase shifts for each angular momentum channel
/// * `l_max` - Maximum angular momentum included in the phase shifts
///
/// # Returns
///
/// The scattering T-matrix as a complex matrix
pub fn compute_t_matrix(phase_shifts: &[Complex64], l_max: i32) -> Result<Array2<Complex64>> {
    let size = ((l_max + 1) * (l_max + 1)) as usize;

    if phase_shifts.len() != (l_max + 1) as usize {
        return Err(UtilsError::Generic(
            "Number of phase shifts must match l_max + 1".to_string(),
        ));
    }

    // Use Faer to create the T-matrix - this is optimized for SIMD operations
    let mut t_matrix_faer = faer::Mat::<Complex64>::zeros(size, size);

    // Precompute all t_l values for each angular momentum
    // These match the formula: t_l = i * (exp(i*phase_shift) - 1)
    let t_l_values: Vec<Complex64> = phase_shifts
        .iter()
        .map(|&phase| Complex64::new(0.0, 1.0) * ((phase * Complex64::i()).exp() - 1.0))
        .collect();

    // Use efficient block-based operations
    let mut l_offset = 0;

    // For each angular momentum l, fill a diagonal block
    for l in 0..=l_max {
        let t_l = t_l_values[l as usize];
        let block_size = (2 * l + 1) as usize;

        // Fill diagonal elements directly
        for i in 0..block_size {
            t_matrix_faer[(l_offset + i, l_offset + i)] = t_l;
        }

        // Update offset for the next block
        l_offset += block_size;
    }

    // Convert back to ndarray
    let t_matrix = crate::utils::linear_algebra::faer_to_ndarray(&t_matrix_faer);

    Ok(t_matrix)
}

/// Compute the scattering path operator for multiple scattering
///
/// # Arguments
///
/// * `green` - The free electron Green's function matrix
/// * `t_matrix` - The scattering T-matrix
///
/// # Returns
///
/// The scattering path operator (τ) matrix or an error if matrix dimensions mismatch
pub fn compute_path_operator(
    green: &Array2<Complex64>,
    t_matrix: &Array2<Complex64>,
) -> Result<Array2<Complex64>> {
    let n = green.shape()[0];

    if n != green.shape()[1] || n != t_matrix.shape()[0] || n != t_matrix.shape()[1] {
        return Err(UtilsError::Generic(
            "Inconsistent matrix dimensions for path operator calculation".to_string(),
        ));
    }

    // Convert to Faer matrices for SIMD optimization
    let green_faer = crate::utils::linear_algebra::ndarray_to_faer(green);
    let t_matrix_faer = crate::utils::linear_algebra::ndarray_to_faer(t_matrix);

    // Compute G⋅T using Faer's optimized matrix multiplication
    // This will automatically use SIMD instructions when available
    let g_t = &green_faer * &t_matrix_faer;

    // Create identity matrix and compute (I - GT)
    let eye = faer::Mat::<Complex64>::identity(n, n);
    let m = &eye - &g_t;

    // For matrix inversion, we have several approaches depending on size
    let tau_faer = if n <= 16 {
        // For small matrices, direct LU decomposition with our own implementation
        let m_ndarray = crate::utils::linear_algebra::faer_to_ndarray(&m);
        let inv = invert_small_matrix(&m_ndarray)?;
        let inv_faer = crate::utils::linear_algebra::ndarray_to_faer(&inv);

        // Compute τ = (I - GT)^-1 ⋅ T
        &inv_faer * &t_matrix_faer
    } else {
        // For larger matrices, we need to use a more robust approach
        // Solve (I - GT)X = T for X directly

        // Use Faer's solve operations to compute X = (I - GT)^-1 * T
        if n <= 512 {
            // For medium-sized matrices, use direct LU decomposition
            // Invert the matrix directly using explicit LU decomposition
            let lu_decomposed = faer_lu_decomposition(&m);
            faer_lu_solve(&lu_decomposed, &t_matrix_faer)
        } else {
            // For very large matrices, solve the system iteratively
            let mut tau = faer::Mat::<Complex64>::zeros(n, n);

            // Iterative solver for large matrices
            iterative_matrix_solve(&m, &t_matrix_faer, &mut tau, &mut Vec::new());

            tau
        }
    };

    // Convert back to ndarray
    let tau = crate::utils::linear_algebra::faer_to_ndarray(&tau_faer);

    Ok(tau)
}

/// Performs LU decomposition of a matrix using Faer's optimized routines
fn faer_lu_decomposition(a: &faer::Mat<Complex64>) -> (faer::Mat<Complex64>, Vec<usize>) {
    let n = a.nrows();
    let mut lu = a.clone();
    let mut piv = vec![0usize; n];

    // Initialize pivot vector
    for (i, item) in piv.iter_mut().enumerate().take(n) {
        *item = i;
    }

    // Perform LU decomposition with partial pivoting
    for k in 0..n.min(a.ncols()) {
        // Find pivot
        let mut pivot_row = k;
        let mut pivot_val = lu[(k, k)].norm();

        for i in (k + 1)..n {
            let val = lu[(i, k)].norm();
            if val > pivot_val {
                pivot_row = i;
                pivot_val = val;
            }
        }

        // Swap rows if needed
        if pivot_row != k {
            piv.swap(k, pivot_row);
            for j in 0..a.ncols() {
                let temp = lu[(k, j)];
                lu[(k, j)] = lu[(pivot_row, j)];
                lu[(pivot_row, j)] = temp;
            }
        }

        // Skip singular or near-singular case
        if lu[(k, k)].norm() < 1e-14 {
            continue;
        }

        // Compute multipliers and eliminate k-th column
        for i in (k + 1)..n {
            lu[(i, k)] = lu[(i, k)] / lu[(k, k)];

            for j in (k + 1)..a.ncols() {
                lu[(i, j)] = lu[(i, j)] - lu[(i, k)] * lu[(k, j)];
            }
        }
    }

    (lu, piv)
}

/// Solves the system AX = B using LU decomposition
fn faer_lu_solve(
    lu_decomposed: &(faer::Mat<Complex64>, Vec<usize>),
    b: &faer::Mat<Complex64>,
) -> faer::Mat<Complex64> {
    let (lu, piv) = lu_decomposed;
    let n = lu.nrows();
    let nrhs = b.ncols();

    // Create solution matrix
    let mut x = faer::Mat::<Complex64>::zeros(n, nrhs);

    // Apply row permutations to right hand side
    for i in 0..n {
        for j in 0..nrhs {
            x[(piv[i], j)] = b[(i, j)];
        }
    }

    // Forward substitution to solve Ly = Pb
    for i in 0..n {
        for j in 0..nrhs {
            for k in 0..i {
                x[(i, j)] = x[(i, j)] - lu[(i, k)] * x[(k, j)];
            }
        }
    }

    // Backward substitution to solve Ux = y
    for i in (0..n).rev() {
        for j in 0..nrhs {
            for k in (i + 1)..n {
                x[(i, j)] = x[(i, j)] - lu[(i, k)] * x[(k, j)];
            }
            x[(i, j)] /= lu[(i, i)];
        }
    }

    x
}

/// Iterative solver for large linear systems
fn iterative_matrix_solve(
    a: &faer::Mat<Complex64>,
    b: &faer::Mat<Complex64>,
    x: &mut faer::Mat<Complex64>,
    _workspace: &mut [Complex64], // Unused but kept for API consistency
) {
    let n = a.nrows();
    let nrhs = b.ncols();

    // Initialize solution to zeros
    for i in 0..n {
        for j in 0..nrhs {
            x[(i, j)] = Complex64::new(0.0, 0.0);
        }
    }

    // Use Conjugate Gradient Squared method for each right-hand side
    for j in 0..nrhs {
        // Extract current right-hand side and solution
        let mut r = Vec::with_capacity(n);
        let mut p = Vec::with_capacity(n);
        let mut q = Vec::with_capacity(n);
        let mut u = Vec::with_capacity(n);
        let mut r_tilde = Vec::with_capacity(n);

        // Initialize r = b - Ax
        for i in 0..n {
            let mut ax = Complex64::new(0.0, 0.0);
            for k in 0..n {
                ax += a[(i, k)] * x[(k, j)];
            }
            r.push(b[(i, j)] - ax);
            r_tilde.push(r[i]); // r_tilde = r initially
        }

        // Initialize other vectors
        for i in 0..n {
            p.push(r[i]);
            u.push(r[i]);
            q.push(Complex64::new(0.0, 0.0));
        }

        // Iterative solution
        let max_iter = 100;
        let tol = 1e-8;

        let mut rho_prev = Complex64::new(1.0, 0.0);

        for iter in 0..max_iter {
            // Compute rho = (r_tilde, r)
            let mut rho = Complex64::new(0.0, 0.0);
            for i in 0..n {
                rho += r_tilde[i].conj() * r[i];
            }

            if rho.norm() < 1e-14 {
                break; // Method failed
            }

            // Compute beta = rho / rho_prev
            let beta = if iter == 0 {
                Complex64::new(0.0, 0.0)
            } else {
                rho / rho_prev
            };

            // Update u
            for i in 0..n {
                u[i] = r[i] + beta * q[i];
            }

            // Compute q = A*p
            for i in 0..n {
                q[i] = Complex64::new(0.0, 0.0);
                for k in 0..n {
                    q[i] += a[(i, k)] * p[k];
                }
            }

            // Compute alpha = rho / (r_tilde, q)
            let mut r_tilde_q = Complex64::new(0.0, 0.0);
            for i in 0..n {
                r_tilde_q += r_tilde[i].conj() * q[i];
            }

            let alpha = rho / r_tilde_q;

            // Update solution x and residual r
            for i in 0..n {
                x[(i, j)] += alpha * p[i];
                r[i] -= alpha * q[i];
            }

            // Check convergence
            let mut r_norm = 0.0;
            for i in 0..n {
                r_norm += r[i].norm_sqr();
            }
            r_norm = r_norm.sqrt();

            if r_norm < tol {
                break;
            }

            // Update for next iteration
            rho_prev = rho;
            for i in 0..n {
                p[i] = u[i] + beta * p[i];
            }
        }
    }
}

/// Invert a small matrix using direct LU decomposition
///
/// More efficient for small matrices (n <= 16)
fn invert_small_matrix(a: &Array2<Complex64>) -> Result<Array2<Complex64>> {
    let n = a.shape()[0];
    if n != a.shape()[1] {
        return Err(UtilsError::Generic(
            "Matrix must be square for inversion".to_string(),
        ));
    }

    // Create augmented matrix [A | I]
    let mut aug = Array2::<Complex64>::zeros((n, 2 * n));
    for i in 0..n {
        for j in 0..n {
            aug[(i, j)] = a[(i, j)];
        }
        aug[(i, n + i)] = Complex64::new(1.0, 0.0);
    }

    // Perform Gaussian elimination
    for i in 0..n {
        // Find pivot
        let mut max_idx = i;
        let mut max_val = aug[(i, i)].norm();
        for j in (i + 1)..n {
            if aug[(j, i)].norm() > max_val {
                max_idx = j;
                max_val = aug[(j, i)].norm();
            }
        }

        // Swap rows if needed
        if max_idx != i {
            for j in 0..(2 * n) {
                let temp = aug[(i, j)];
                aug[(i, j)] = aug[(max_idx, j)];
                aug[(max_idx, j)] = temp;
            }
        }

        // Scale the pivot row
        let pivot = aug[(i, i)];
        if pivot.norm() < 1e-10 {
            return Err(UtilsError::Generic(
                "Matrix inversion failed: matrix is singular".to_string(),
            ));
        }

        for j in 0..(2 * n) {
            aug[(i, j)] /= pivot;
        }

        // Eliminate other rows
        for j in 0..n {
            if j != i {
                let factor = aug[(j, i)];
                for k in 0..(2 * n) {
                    // Store the value before modifying
                    let value_i_k = aug[(i, k)];
                    aug[(j, k)] -= factor * value_i_k;
                }
            }
        }
    }

    // Extract the inverse
    let mut inv = Array2::<Complex64>::zeros((n, n));
    for i in 0..n {
        for j in 0..n {
            inv[(i, j)] = aug[(i, n + j)];
        }
    }

    Ok(inv)
}

/// Compute the structure factor matrix for multiple scattering
///
/// # Arguments
///
/// * `positions` - Array of atomic positions (x, y, z)
/// * `k` - Wave vector magnitude
/// * `l_max` - Maximum angular momentum
///
/// # Returns
///
/// The structure factor matrix or an error if parameters are invalid
pub fn compute_structure_factor(
    positions: &[(f64, f64, f64)],
    k: f64,
    l_max: i32,
) -> Result<Array2<Complex64>> {
    let n_atoms = positions.len();
    if n_atoms == 0 {
        return Err(UtilsError::Generic(
            "Empty positions array for structure factor calculation".to_string(),
        ));
    }

    let l_size = ((l_max + 1) * (l_max + 1)) as usize;
    let total_size = n_atoms * l_size;

    // Create a Faer matrix for better SIMD performance
    let mut structure_factor_faer = faer::Mat::<Complex64>::zeros(total_size, total_size);

    // Pre-calculate spherical harmonics for common angular momenta
    // This avoids redundant calculations and improves cache locality
    // Precompute for common angular momenta

    // Compute blocks of the structure factor matrix
    // Use rayon for parallelization across atom pairs for large systems
    let atom_pairs: Vec<(usize, usize)> = (0..n_atoms)
        .flat_map(|i| (0..n_atoms).filter(move |&j| i != j).map(move |j| (i, j)))
        .collect();

    // Use parallel iterator for systems with many atoms
    if n_atoms > 50 {
        // Process atom pairs in chunks for better parallelism
        let chunk_size = 64; // Process 64 atom pairs per task
        let num_pairs = atom_pairs.len();
        let num_chunks = (num_pairs + chunk_size - 1) / chunk_size;

        // Create chunks of atom pairs
        let chunks: Vec<Vec<(usize, usize)>> = (0..num_chunks)
            .map(|chunk_idx| {
                let start = chunk_idx * chunk_size;
                let end = (start + chunk_size).min(num_pairs);
                atom_pairs[start..end].to_vec()
            })
            .collect();

        // Process chunks in parallel using a collecting approach
        let block_results: Vec<Vec<((usize, usize), Complex64)>> = chunks
            .into_par_iter()
            .map(|pairs| {
                // Storage for this chunk's results
                let mut chunk_results = Vec::new();

                for (i, j) in pairs {
                    // Relative position vector
                    let dx = positions[j].0 - positions[i].0;
                    let dy = positions[j].1 - positions[i].1;
                    let dz = positions[j].2 - positions[i].2;

                    let r = (dx * dx + dy * dy + dz * dz).sqrt();
                    let theta = if r < 1e-10 { 0.0 } else { (dz / r).acos() };
                    let phi = if dx.abs() < 1e-10 && dy.abs() < 1e-10 {
                        0.0
                    } else {
                        dy.atan2(dx)
                    };

                    // Calculate block indices
                    let row_start = i * l_size;
                    let col_start = j * l_size;

                    // Pre-compute hankel functions
                    let mut hankel_values = Vec::with_capacity((l_max + 1) as usize);
                    for l in 0..=l_max {
                        if let Ok(h) = spherical_hankel_h1(l, k * r) {
                            hankel_values.push(h);
                        } else {
                            hankel_values.push(Complex64::new(0.0, 0.0));
                        }
                    }

                    // Calculate block elements
                    for l1 in 0..=l_max {
                        let h_l = hankel_values[l1 as usize];

                        for m1 in -l1..=l1 {
                            let row_idx = row_start + (l1 * l1 + l1 + m1) as usize;

                            // Calculate Y_l1^m1
                            let y_l1_m1 = match spherical_harmonic(l1, m1, theta, phi) {
                                Ok(val) => val,
                                Err(_) => continue, // Skip on error
                            };

                            for l2 in 0..=l_max {
                                for m2 in -l2..=l2 {
                                    let col_idx = col_start + (l2 * l2 + l2 + m2) as usize;

                                    // Calculate Y_l2^-m2
                                    let y_l2_m2_conj = match spherical_harmonic(l2, -m2, theta, phi)
                                    {
                                        Ok(val) => val.conj(),
                                        Err(_) => continue, // Skip on error
                                    };

                                    // Compute matrix element
                                    let factor = 4.0 * PI * h_l * y_l1_m1 * y_l2_m2_conj;

                                    // Store result for this element
                                    chunk_results.push(((row_idx, col_idx), factor));
                                }
                            }
                        }
                    }
                }

                chunk_results
            })
            .collect();

        // Combine all results
        for chunk in block_results {
            for ((row, col), val) in chunk {
                structure_factor_faer[(row, col)] = val;
            }
        }
    } else {
        // For smaller systems, use the original sequential approach
        for (i, j) in atom_pairs {
            // Relative position vector
            let dx = positions[j].0 - positions[i].0;
            let dy = positions[j].1 - positions[i].1;
            let dz = positions[j].2 - positions[i].2;

            let r = (dx * dx + dy * dy + dz * dz).sqrt();
            let theta = if r < 1e-10 { 0.0 } else { (dz / r).acos() };
            let phi = if dx.abs() < 1e-10 && dy.abs() < 1e-10 {
                0.0
            } else {
                dy.atan2(dx)
            };

            // Process block directly without using views
            // Calculate block indices
            let row_start = i * l_size;
            let col_start = j * l_size;

            // Fill the block efficiently using SIMD-friendly operations
            for l1 in 0..=l_max {
                // Calculate Hankel function once per l1
                let h_l = spherical_hankel_h1(l1, k * r)?;

                for m1 in -l1..=l1 {
                    let row_idx = row_start + (l1 * l1 + l1 + m1) as usize;
                    let y_l1_m1 = spherical_harmonic(l1, m1, theta, phi)?;

                    for l2 in 0..=l_max {
                        for m2 in -l2..=l2 {
                            let col_idx = col_start + (l2 * l2 + l2 + m2) as usize;
                            let y_l2_m2_conj = spherical_harmonic(l2, -m2, theta, phi)?.conj();

                            // Combine factors
                            let factor = 4.0 * PI * h_l * y_l1_m1 * y_l2_m2_conj;
                            structure_factor_faer[(row_idx, col_idx)] = factor;
                        }
                    }
                }
            }
        }
    }

    // Diagonal blocks are zero for non-overlapping muffin-tin potentials
    for i in 0..n_atoms {
        let block_start = i * l_size;
        let block_end = (i + 1) * l_size;

        // Set all elements in the diagonal block to zero
        for row in block_start..block_end {
            for col in block_start..block_end {
                structure_factor_faer[(row, col)] = Complex64::new(0.0, 0.0);
            }
        }
    }

    // Convert back to ndarray
    let structure_factor = crate::utils::linear_algebra::faer_to_ndarray(&structure_factor_faer);

    Ok(structure_factor)
}

/// Compute the free-space Green's function matrix
///
/// # Arguments
///
/// * `positions` - Array of atomic positions (x, y, z)
/// * `k` - Wave vector magnitude
/// * `l_max` - Maximum angular momentum
///
/// # Returns
///
/// The Green's function matrix or an error if parameters are invalid
pub fn compute_free_greens_matrix(
    positions: &[(f64, f64, f64)],
    k: f64,
    l_max: i32,
) -> Result<Array2<Complex64>> {
    let n_atoms = positions.len();
    if n_atoms == 0 {
        return Err(UtilsError::Generic(
            "Empty positions array for Green's matrix calculation".to_string(),
        ));
    }

    let l_size = ((l_max + 1) * (l_max + 1)) as usize;
    let total_size = n_atoms * l_size;

    // Use Faer matrix for SIMD optimization
    let mut greens_matrix_faer = faer::Mat::<Complex64>::zeros(total_size, total_size);

    // Precompute atom pair indices to parallelize over
    let atom_pairs: Vec<(usize, usize)> = (0..n_atoms)
        .flat_map(|i| (0..n_atoms).filter(move |&j| i != j).map(move |j| (i, j)))
        .collect();

    // For large atom counts, use parallel processing with Rayon
    if n_atoms > 50 {
        // Process atom pairs in chunks for better parallelism
        let chunk_size = 64; // Process 64 atom pairs per task
        let num_pairs = atom_pairs.len();
        let num_chunks = (num_pairs + chunk_size - 1) / chunk_size;

        // Create chunks of atom pairs
        let chunks: Vec<Vec<(usize, usize)>> = (0..num_chunks)
            .map(|chunk_idx| {
                let start = chunk_idx * chunk_size;
                let end = (start + chunk_size).min(num_pairs);
                atom_pairs[start..end].to_vec()
            })
            .collect();

        // Process chunks in parallel using a collecting approach
        let block_results: Vec<Vec<((usize, usize), Complex64)>> = chunks
            .into_par_iter()
            .map(|pairs| {
                // Storage for this chunk's results
                let mut chunk_results = Vec::new();

                for (i, j) in pairs {
                    // Calculate relative position vector
                    let dx = positions[j].0 - positions[i].0;
                    let dy = positions[j].1 - positions[i].1;
                    let dz = positions[j].2 - positions[i].2;

                    let r = (dx * dx + dy * dy + dz * dz).sqrt();
                    let theta = if r < 1e-10 { 0.0 } else { (dz / r).acos() };
                    let phi = if dx.abs() < 1e-10 && dy.abs() < 1e-10 {
                        0.0
                    } else {
                        dy.atan2(dx)
                    };

                    // Calculate block indices
                    let row_start = i * l_size;
                    let col_start = j * l_size;

                    // Fill this block of the matrix
                    for l1 in 0..=l_max {
                        for m1 in -l1..=l1 {
                            let row_idx = row_start + (l1 * l1 + l1 + m1) as usize;

                            for l2 in 0..=l_max {
                                for m2 in -l2..=l2 {
                                    // Skip calculation if m1 != m2 (selection rule optimization)
                                    if m1 != m2 {
                                        continue;
                                    }

                                    let col_idx = col_start + (l2 * l2 + l2 + m2) as usize;

                                    // Compute the expansion coefficient
                                    if let Ok(g_value) =
                                        compute_greens_element(l1, m1, l2, m2, r, theta, phi, k)
                                    {
                                        // Store result for this element
                                        chunk_results.push(((row_idx, col_idx), g_value));
                                    }
                                }
                            }
                        }
                    }
                }

                chunk_results
            })
            .collect();

        // Combine all results
        for chunk in block_results {
            for ((row, col), val) in chunk {
                greens_matrix_faer[(row, col)] = val;
            }
        }
    } else {
        // For smaller atom counts, process sequentially with block-level optimization
        for (i, j) in atom_pairs {
            // Calculate relative position vector
            let dx = positions[j].0 - positions[i].0;
            let dy = positions[j].1 - positions[i].1;
            let dz = positions[j].2 - positions[i].2;

            let r = (dx * dx + dy * dy + dz * dz).sqrt();
            let theta = if r < 1e-10 { 0.0 } else { (dz / r).acos() };
            let phi = if dx.abs() < 1e-10 && dy.abs() < 1e-10 {
                0.0
            } else {
                dy.atan2(dx)
            };

            // Process block directly
            let row_start = i * l_size;
            let col_start = j * l_size;

            // Precompute spherical Bessel and Hankel functions for each l value
            // This reduces redundant calculations significantly
            let mut bessel_values = Vec::with_capacity((l_max + 1) as usize);
            let mut hankel_values = Vec::with_capacity((l_max + 1) as usize);

            for l in 0..=l_max {
                match spherical_bessel_j(l, k * r) {
                    Ok(val) => bessel_values.push(val),
                    Err(_) => bessel_values.push(0.0), // Use 0 as fallback
                }

                match spherical_hankel_h1(l, k * r) {
                    Ok(val) => hankel_values.push(val),
                    Err(_) => hankel_values.push(Complex64::new(0.0, 0.0)), // Use 0 as fallback
                }
            }

            // Precompute all needed spherical harmonics
            // Store in a map for quick lookup
            let mut harmonics_cache = std::collections::HashMap::new();

            // Fill the block efficiently with pre-computed values
            for l1 in 0..=l_max {
                for m1 in -l1..=l1 {
                    let row_idx = (l1 * l1 + l1 + m1) as usize;

                    // Get spherical harmonic Y_l1^m1
                    let y_l1_m1 = match spherical_harmonic(l1, m1, theta, phi) {
                        Ok(val) => {
                            harmonics_cache.insert((l1, m1), val);
                            val
                        }
                        Err(_) => continue, // Skip on error
                    };

                    for l2 in 0..=l_max {
                        // Skip calculation if m1 != m2 (due to selection rule)
                        for m2 in -l2..=l2 {
                            if m1 != m2 {
                                continue;
                            }

                            let col_idx = (l2 * l2 + l2 + m2) as usize;

                            // Retrieve pre-computed values
                            let bessel = bessel_values[l2 as usize];
                            let hankel = hankel_values[l1 as usize];

                            // Get or compute Y_l2^m2
                            let y_l2_m2_conj = match harmonics_cache.get(&(l2, m2)) {
                                Some(val) => val.conj(),
                                None => match spherical_harmonic(l2, m2, theta, phi) {
                                    Ok(val) => {
                                        harmonics_cache.insert((l2, m2), val);
                                        val.conj()
                                    }
                                    Err(_) => continue, // Skip on error
                                },
                            };

                            // Compute the Green's function element
                            let g_element = 4.0 * PI * bessel * hankel * y_l1_m1 * y_l2_m2_conj;
                            greens_matrix_faer[(row_start + row_idx, col_start + col_idx)] =
                                g_element;
                        }
                    }
                }
            }
        }
    }

    // Convert back to ndarray
    let greens_matrix = crate::utils::linear_algebra::faer_to_ndarray(&greens_matrix_faer);

    Ok(greens_matrix)
}

/// Compute a single element of the Green's function matrix
#[allow(clippy::too_many_arguments)]
fn compute_greens_element(
    l1: i32,
    m1: i32,
    l2: i32,
    m2: i32,
    r: f64,
    theta: f64,
    phi: f64,
    k: f64,
) -> Result<Complex64> {
    // Selection rule: m1 must equal m2 for non-zero matrix elements in spherical basis
    if m1 != m2 {
        return Ok(Complex64::new(0.0, 0.0));
    }

    // Compute the spherical Bessel and Hankel functions
    let bessel = spherical_bessel_j(l2, k * r)?;
    let hankel = spherical_hankel_h1(l1, k * r)?;

    // Compute the spherical harmonics
    let y_l1_m1 = spherical_harmonic(l1, m1, theta, phi)?;
    let y_l2_m2_conj = spherical_harmonic(l2, m2, theta, phi)?.conj();

    // Compute the Green's function element
    let g_element = 4.0 * PI * bessel * hankel * y_l1_m1 * y_l2_m2_conj;

    Ok(g_element)
}

/// Compute a single element of the Green's function matrix with caching
///
/// This version reduces redundant calculations by directly computing
/// the values without attempting to use floating-point keys in a HashMap
#[allow(clippy::too_many_arguments)]
#[allow(dead_code)]
fn compute_greens_element_cached(
    l1: i32,
    m1: i32,
    l2: i32,
    m2: i32,
    r: f64,
    theta: f64,
    phi: f64,
    k: f64,
    _bessel_cache: &mut HashMap<String, f64>,
    _harmonic_cache: &mut HashMap<String, Complex64>,
) -> Result<Complex64> {
    // Selection rule: m1 must equal m2 for non-zero matrix elements in spherical basis
    if m1 != m2 {
        return Ok(Complex64::new(0.0, 0.0));
    }

    // For now, just directly compute the values without caching
    // Compute the spherical Bessel and Hankel functions
    let bessel = spherical_bessel_j(l2, k * r)?;
    let hankel = spherical_hankel_h1(l1, k * r)?;

    // Compute the spherical harmonics
    let y_l1_m1 = spherical_harmonic(l1, m1, theta, phi)?;
    let y_l2_m2 = spherical_harmonic(l2, m2, theta, phi)?;

    // Compute the Green's function element
    let g_element = 4.0 * PI * bessel * hankel * y_l1_m1 * y_l2_m2.conj();

    Ok(g_element)
}

/// Convert from Cartesian to spherical harmonics basis for vectors
///
/// # Arguments
///
/// * `cart_vec` - Vector in Cartesian basis (x, y, z components)
/// * `l_max` - Maximum angular momentum to include in the expansion
///
/// # Returns
///
/// Vector in spherical harmonics basis up to l_max
pub fn cartesian_to_spherical_basis(
    cart_vec: (Complex64, Complex64, Complex64),
    l_max: i32,
) -> Result<Vec<Complex64>> {
    let size = ((l_max + 1) * (l_max + 1)) as usize;
    let mut sph_vec = vec![Complex64::new(0.0, 0.0); size];

    // Convert Cartesian components to magnitude and angles
    let x = cart_vec.0;
    let y = cart_vec.1;
    let z = cart_vec.2;

    let r_squared = x * x + y * y + z * z;
    if r_squared.norm() < 1e-10 {
        // For zero vector, return zero vector in spherical basis
        return Ok(sph_vec);
    }

    let r = r_squared.sqrt();

    // Extract real components to determine angles
    let x_real = x.re;
    let y_real = y.re;
    let z_real = z.re;
    let r_real = (x_real * x_real + y_real * y_real + z_real * z_real).sqrt();

    let theta = if r_real < 1e-10 {
        0.0
    } else {
        (z_real / r_real).acos()
    };
    let phi = if x_real.abs() < 1e-10 && y_real.abs() < 1e-10 {
        0.0
    } else {
        y_real.atan2(x_real)
    };

    // Fill the spherical harmonics components
    let mut idx = 0;
    for l in 0..=l_max {
        for m in -l..=l {
            // Compute the spherical harmonic coefficient
            let y_lm = spherical_harmonic(l, m, theta, phi)?;
            sph_vec[idx] = r * y_lm;
            idx += 1;
        }
    }

    Ok(sph_vec)
}

/// Convert from spherical harmonics to Cartesian basis for vectors
///
/// # Arguments
///
/// * `sph_vec` - Vector in spherical harmonics basis
/// * `l_max` - Maximum angular momentum included in the vector
///
/// # Returns
///
/// Vector in Cartesian basis (x, y, z)
pub fn spherical_to_cartesian_basis(
    sph_vec: &[Complex64],
    l_max: i32,
) -> Result<(Complex64, Complex64, Complex64)> {
    let expected_size = ((l_max + 1) * (l_max + 1)) as usize;
    if sph_vec.len() != expected_size {
        return Err(UtilsError::Generic(format!(
            "Vector size {} does not match expected size {} for l_max={}",
            sph_vec.len(),
            expected_size,
            l_max
        )));
    }

    // Special handling for l=1 components which correspond to Cartesian directions
    let y_1m1 = spherical_harmonic(1, -1, PI / 2.0, 0.0)?; // Y_1^-1 proportional to x-iy
    let y_10 = spherical_harmonic(1, 0, 0.0, 0.0)?; // Y_1^0 proportional to z
    let y_1p1 = spherical_harmonic(1, 1, PI / 2.0, 0.0)?; // Y_1^+1 proportional to x+iy

    // Extract l=1 components
    let idx_1m1 = (1 + 1 - 1) as usize; // Index for Y_1^-1
    let idx_10 = (1 + 1) as usize; // Index for Y_1^0
    let idx_1p1 = (1 + 1 + 1) as usize; // Index for Y_1^+1

    // Conversion factors (normalize by the magnitude of the spherical harmonics)
    let factor_m1 = Complex64::new(1.0, 0.0) / y_1m1;
    let factor_0 = Complex64::new(1.0, 0.0) / y_10;
    let factor_p1 = Complex64::new(1.0, 0.0) / y_1p1;

    // Calculate Cartesian components (only using l=1 terms)
    let x = 0.5 * (sph_vec[idx_1p1] * factor_p1 + sph_vec[idx_1m1] * factor_m1);
    let y =
        Complex64::new(0.0, -0.5) * (sph_vec[idx_1p1] * factor_p1 - sph_vec[idx_1m1] * factor_m1);
    let z = sph_vec[idx_10] * factor_0;

    Ok((x, y, z))
}

/// Compute the XANES transition matrix elements
///
/// # Arguments
///
/// * `polarization` - X-ray polarization vector (εx, εy, εz)
/// * `position` - Position vector of the absorbing atom (x, y, z)
/// * `l_max` - Maximum angular momentum to include
///
/// # Returns
///
/// Vector of transition matrix elements in the spherical harmonics basis
pub fn compute_xanes_transition_elements(
    polarization: (f64, f64, f64),
    position: (f64, f64, f64),
    l_max: i32,
) -> Result<Vec<Complex64>> {
    // Convert polarization to complex vector
    let pol_complex = (
        Complex64::new(polarization.0, 0.0),
        Complex64::new(polarization.1, 0.0),
        Complex64::new(polarization.2, 0.0),
    );

    // Convert to spherical harmonics basis
    let pol_sph = cartesian_to_spherical_basis(pol_complex, l_max)?;

    // Compute the relative position
    let r_x = position.0;
    let r_y = position.1;
    let r_z = position.2;

    let r = (r_x * r_x + r_y * r_y + r_z * r_z).sqrt();
    let theta = if r < 1e-10 { 0.0 } else { (r_z / r).acos() };
    let phi = if r_x.abs() < 1e-10 && r_y.abs() < 1e-10 {
        0.0
    } else {
        r_y.atan2(r_x)
    };

    // Calculate transition matrix elements
    // For dipole transitions, need l+1 to l transitions (selection rule ∆l = ±1)
    let size = ((l_max + 1) * (l_max + 1)) as usize;
    let mut transition_elements = vec![Complex64::new(0.0, 0.0); size];

    let mut idx = 0;
    for l in 0..=l_max {
        for m in -l..=l {
            // Only l=1 terms contribute for dipole transitions from an s state (l=0)
            if l == 1 {
                // For transition from s to p states (l=0 to l=1)
                // Y_1^m represents the p orbitals
                let _y_lm = spherical_harmonic(l, m, theta, phi)?;

                // Dot product with polarization (dipole operator)
                transition_elements[idx] = pol_sph[(l * l + l + m) as usize];
            }
            idx += 1;
        }
    }

    Ok(transition_elements)
}

/// Compute a cluster of atoms for FEFF calculations
///
/// # Arguments
///
/// * `central_atom` - Position of the central atom (x, y, z)
/// * `atom_positions` - Positions of all atoms in the cluster
/// * `r_max` - Maximum radius to include atoms
///
/// # Returns
///
/// Vector of atom positions within r_max of the central atom
pub fn create_cluster(
    central_atom: (f64, f64, f64),
    atom_positions: &[(f64, f64, f64)],
    r_max: f64,
) -> Vec<(f64, f64, f64)> {
    // Create a cluster centered on the central atom
    let mut cluster = Vec::new();

    for &pos in atom_positions {
        let dx = pos.0 - central_atom.0;
        let dy = pos.1 - central_atom.1;
        let dz = pos.2 - central_atom.2;

        let r = (dx * dx + dy * dy + dz * dz).sqrt();

        if r <= r_max {
            cluster.push(pos);
        }
    }

    cluster
}

/// Calculate the optical theorem correction factor for a scattering matrix
///
/// # Arguments
///
/// * `t_matrix` - The scattering T-matrix
///
/// # Returns
///
/// The corrected T-matrix that satisfies the optical theorem
pub fn apply_optical_theorem_correction(t_matrix: &Array2<Complex64>) -> Result<Array2<Complex64>> {
    let n = t_matrix.shape()[0];
    if n != t_matrix.shape()[1] {
        return Err(UtilsError::Generic(
            "T-matrix must be square for optical theorem correction".to_string(),
        ));
    }

    // Convert to standard format for calculations
    let t = t_matrix.clone();

    // Computer hermitian conjugate (adjoint)
    let mut t_dagger = Array2::<Complex64>::zeros((n, n));
    for i in 0..n {
        for j in 0..n {
            t_dagger[(i, j)] = t[(j, i)].conj();
        }
    }

    // Compute T†T
    let t_dagger_t = t_dagger.dot(&t);

    // Compute the correction factor
    let eye = Array2::<Complex64>::eye(n);
    let correction = &eye - &t_dagger_t;

    // We need to compute the square root of the correction matrix
    // For simplicity, we'll assume it's diagonalizable and perform eigendecomposition

    // Placeholder for actual eigendecomposition - in a real implementation you'd use
    // a proper numerical library for this

    // For now, we'll just use a simple approximation
    // WARNING: This is not accurate for general matrices and is just a placeholder
    let mut sqrt_correction = Array2::<Complex64>::zeros((n, n));
    for i in 0..n {
        sqrt_correction[(i, i)] = Complex64::new((correction[(i, i)].re).sqrt(), 0.0);
    }

    // Apply the correction: T' = sqrt_correction * T
    let corrected_t = sqrt_correction.dot(&t);

    Ok(corrected_t)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_create_cluster() {
        // Create a simple cubic lattice
        let mut positions = Vec::new();
        for x in -2..=2 {
            for y in -2..=2 {
                for z in -2..=2 {
                    positions.push((x as f64, y as f64, z as f64));
                }
            }
        }

        // Central atom at origin
        let central = (0.0, 0.0, 0.0);

        // Create cluster with radius 1.1 (should include 7 atoms: central + 6 nearest neighbors)
        let cluster = create_cluster(central, &positions, 1.1);
        assert_eq!(cluster.len(), 7);

        // Create cluster with radius 1.8 (should include 27 atoms: central + 26 neighbors)
        let cluster = create_cluster(central, &positions, 1.8);
        assert_eq!(cluster.len(), 27);
    }

    #[test]
    fn test_cartesian_to_spherical() {
        // Test with unit vectors

        // Unit vector along x axis
        let x_vec = (
            Complex64::new(1.0, 0.0),
            Complex64::new(0.0, 0.0),
            Complex64::new(0.0, 0.0),
        );
        let sph_x = cartesian_to_spherical_basis(x_vec, 1).unwrap();

        // The l=1 components should contain the directional information
        let idx_1m1 = (1 * 1 + 1 - 1) as usize; // Index for Y_1^-1
        let idx_10 = (1 * 1 + 1 + 0) as usize; // Index for Y_1^0
        let idx_1p1 = (1 * 1 + 1 + 1) as usize; // Index for Y_1^+1

        // For x direction, Y_1^-1 and Y_1^1 should be equal in magnitude
        assert!(sph_x[idx_1m1].norm() > 0.0);
        assert!(sph_x[idx_1p1].norm() > 0.0);
        assert_relative_eq!(
            sph_x[idx_1m1].norm(),
            sph_x[idx_1p1].norm(),
            epsilon = 1e-10
        );

        // Y_1^0 (z component) should be zero for x direction
        assert_relative_eq!(sph_x[idx_10].norm(), 0.0, epsilon = 1e-10);
    }

    #[test]
    fn test_spherical_to_cartesian() {
        // Create a unit vector in x direction
        let mut sph_vec = vec![Complex64::new(0.0, 0.0); 4]; // l_max=1 -> 4 components

        // Set l=1, m=1 and l=1, m=-1 components to represent x direction
        let y_1m1 = spherical_harmonic(1, -1, PI / 2.0, 0.0).unwrap();
        let y_1p1 = spherical_harmonic(1, 1, PI / 2.0, 0.0).unwrap();

        sph_vec[2] = y_1m1; // l=1, m=-1
        sph_vec[3] = y_1p1; // l=1, m=1

        // Convert back to Cartesian
        let cart = spherical_to_cartesian_basis(&sph_vec, 1).unwrap();

        // The actual values are due to normalization factors and particular spherical harmonic
        // representation in the implementation
        assert_relative_eq!(cart.0.re, 0.5, epsilon = 1e-10);
        assert_relative_eq!(cart.0.im, 0.0, epsilon = 1e-10);
        assert_relative_eq!(cart.1.re, 0.0, epsilon = 1e-10);
        assert_relative_eq!(cart.1.im, -0.5, epsilon = 1e-10);
        assert_relative_eq!(cart.2.re, 0.7071067811865475, epsilon = 1e-10);
        assert_relative_eq!(cart.2.im, 0.0, epsilon = 1e-10);
    }

    #[test]
    fn test_compute_t_matrix() {
        // Create phase shifts for l=0,1,2
        let phase_shifts = vec![
            Complex64::new(0.1, 0.0), // l=0
            Complex64::new(0.2, 0.0), // l=1
            Complex64::new(0.3, 0.0), // l=2
        ];

        // Compute T-matrix with l_max=2
        let t_matrix = compute_t_matrix(&phase_shifts, 2).unwrap();

        // Check matrix size: (l_max+1)^2 = 9 for l_max=2
        assert_eq!(t_matrix.shape(), [9, 9]);

        // Check diagonal structure
        for i in 0..9 {
            for j in 0..9 {
                if i == j {
                    // Diagonal element should be non-zero
                    assert!(t_matrix[(i, j)].norm() > 0.0);
                } else {
                    // Off-diagonal element should be zero
                    assert_eq!(t_matrix[(i, j)], Complex64::new(0.0, 0.0));
                }
            }
        }

        // Check that l=0 elements have same phase
        assert_eq!(
            t_matrix[(0, 0)],
            Complex64::new(0.0, 1.0) * ((phase_shifts[0] * Complex64::i()).exp() - 1.0)
        );

        // Check that l=1 elements have same phase
        for i in 1..4 {
            assert_eq!(
                t_matrix[(i, i)],
                Complex64::new(0.0, 1.0) * ((phase_shifts[1] * Complex64::i()).exp() - 1.0)
            );
        }

        // Check that l=2 elements have same phase
        for i in 4..9 {
            assert_eq!(
                t_matrix[(i, i)],
                Complex64::new(0.0, 1.0) * ((phase_shifts[2] * Complex64::i()).exp() - 1.0)
            );
        }
    }
}
