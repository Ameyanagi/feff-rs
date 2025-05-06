/*
MIT License with FEFF10 Attribution

Copyright (c) 2025 Ameyanagi

Based on or developed using Distribution: FEFF10.0
Copyright (c) 2020 FEFF Project, University of Washington and SLAC National Accelerator Laboratory.
All rights reserved.
*/

//! FMS equation solver
//!
//! This module implements different methods for solving the FMS equations.

use super::errors::{FmsError, Result};
use crate::utils::linear_algebra::{faer_to_ndarray, inner_product, ndarray_to_faer, vector_norm};
use faer::{col, Mat};
use ndarray::Array2;
use num_complex::Complex64;
use rayon::prelude::*;

/// Methods for solving the FMS equations
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SolverMethod {
    /// Direct LU decomposition
    LuDecomposition,
    /// Iterative solver using conjugate gradient squared method
    IterativeCGS,
    /// Block diagonal approximation
    BlockDiagonal,
}

/// Solver for FMS equations
///
/// This struct provides methods for solving the FMS equations (I - GT)τ = T,
/// which determine the full multiple scattering path operator τ.
#[derive(Debug)]
pub struct FmsSolver {
    /// Method to use for solving the equations
    method: SolverMethod,
    /// Convergence tolerance for iterative methods
    tolerance: f64,
    /// Maximum number of iterations for iterative methods
    max_iterations: usize,
}

impl FmsSolver {
    /// Create a new FMS solver with the specified method
    ///
    /// # Arguments
    ///
    /// * `method` - Method to use for solving the equations
    ///
    /// # Returns
    ///
    /// A new FMS solver
    pub fn new(method: SolverMethod) -> Self {
        Self {
            method,
            tolerance: 1e-6,
            max_iterations: 100,
        }
    }

    /// Set the convergence tolerance for iterative methods
    ///
    /// # Arguments
    ///
    /// * `tolerance` - Convergence tolerance
    pub fn set_tolerance(&mut self, tolerance: f64) -> &mut Self {
        self.tolerance = tolerance;
        self
    }

    /// Set the maximum number of iterations for iterative methods
    ///
    /// # Arguments
    ///
    /// * `max_iterations` - Maximum number of iterations
    pub fn set_max_iterations(&mut self, max_iterations: usize) -> &mut Self {
        self.max_iterations = max_iterations;
        self
    }

    /// Solve the FMS equations to get the path operator τ
    ///
    /// The FMS equations are (I - GT)τ = T, where G is the free electron
    /// propagator, T is the scattering matrix, and τ is the full multiple
    /// scattering path operator.
    ///
    /// # Arguments
    ///
    /// * `fms_matrix` - The FMS matrix (I - GT)
    ///
    /// # Returns
    ///
    /// The path operator τ
    pub fn solve(&self, fms_matrix: &Array2<Complex64>) -> Result<Array2<Complex64>> {
        match self.method {
            SolverMethod::LuDecomposition => self.solve_lu(fms_matrix),
            SolverMethod::IterativeCGS => self.solve_iterative(fms_matrix),
            SolverMethod::BlockDiagonal => self.solve_block_diagonal(fms_matrix),
        }
    }

    /// Solve the FMS equations using LU decomposition
    ///
    /// This method is most suitable for small to medium-sized matrices
    /// (up to a few thousand atoms).
    ///
    /// # Arguments
    ///
    /// * `fms_matrix` - The FMS matrix (I - GT)
    ///
    /// # Returns
    ///
    /// The path operator τ
    fn solve_lu(&self, fms_matrix: &Array2<Complex64>) -> Result<Array2<Complex64>> {
        // Convert to faer matrix for better performance
        let matrix_faer = ndarray_to_faer(fms_matrix);
        let n = matrix_faer.nrows();

        // For τ = (I - GT)^(-1) T, we need to compute the inverse
        // We'll use LU decomposition, which is efficient for medium-sized matrices

        // Perform manual LU decomposition with optimizations for performance
        // First, clone the matrix
        let mut a = matrix_faer.clone();

        // Create permutation vector (initially identity)
        let mut p = (0..n).collect::<Vec<usize>>();

        // Create the L and U matrices
        let mut l = Mat::<Complex64>::zeros(n, n);
        let mut u = Mat::<Complex64>::zeros(n, n);

        // LU decomposition with partial pivoting
        // For large matrices, split into threshold-based parallel and sequential processing
        if n > 64 {
            // For larger matrices, use blocked LU decomposition to improve cache locality
            let block_size = 32; // Chosen to fit in L1 cache

            for k_block in (0..n).step_by(block_size) {
                let k_end = (k_block + block_size).min(n);

                // Process this diagonal block
                for k in k_block..k_end {
                    // Find pivot within the current column
                    let mut pivot_row = k;
                    let mut pivot_val = a[(k, k)].norm();

                    for i in (k + 1)..n {
                        let val = a[(i, k)].norm();
                        if val > pivot_val {
                            pivot_row = i;
                            pivot_val = val;
                        }
                    }

                    // Swap rows if necessary
                    if pivot_row != k {
                        // Use SIMD-friendly approach for row swapping
                        // Process in chunks for better cache locality
                        for j_chunk in (0..n).step_by(32) {
                            let j_end = (j_chunk + 32).min(n);
                            for j in j_chunk..j_end {
                                let temp = a[(k, j)];
                                a[(k, j)] = a[(pivot_row, j)];
                                a[(pivot_row, j)] = temp;
                            }
                        }

                        // Update permutation
                        p.swap(k, pivot_row);
                    }

                    // Fill in L and U for this row/column
                    for i in k..n {
                        u[(k, i)] = a[(k, i)];
                    }

                    // Use parallelism for large matrices during update step
                    if n > 1000 && (n - k) > 64 {
                        // In parallel, compute updates for L and A
                        let remaining_rows = (k + 1)..n;

                        // Compute L elements and A updates in parallel
                        let updates: Vec<(usize, Vec<(usize, Complex64)>)> = remaining_rows
                            .into_par_iter()
                            .map(|i| {
                                // Calculate the L element
                                let l_element = a[(i, k)] / u[(k, k)];

                                // Compute updates for row i of matrix A
                                let mut row_updates = Vec::new();
                                for j in k..n {
                                    if j >= k {
                                        let update = -l_element * u[(k, j)];
                                        row_updates.push((j, update));
                                    }
                                }

                                (i, row_updates)
                            })
                            .collect();

                        // Apply the updates to L and A
                        for (i, row_updates) in updates {
                            // Update L element
                            l[(i, k)] = a[(i, k)] / u[(k, k)];

                            // Apply updates to A
                            for (j, update) in row_updates {
                                a[(i, j)] += update;
                            }
                        }
                    } else {
                        // Sequential processing for smaller matrices or near the end
                        for i in (k + 1)..n {
                            l[(i, k)] = a[(i, k)] / u[(k, k)];

                            // Process in chunks for better cache locality
                            for j_chunk in (k..n).step_by(32) {
                                let j_end = (j_chunk + 32).min(n);
                                for j in j_chunk..j_end {
                                    a[(i, j)] -= l[(i, k)] * u[(k, j)];
                                }
                            }
                        }
                    }
                }
            }
        } else {
            // For smaller matrices, use sequential processing to avoid overhead
            for k in 0..n {
                // Find pivot
                let mut pivot_row = k;
                let mut pivot_val = a[(k, k)].norm();

                for i in (k + 1)..n {
                    let val = a[(i, k)].norm();
                    if val > pivot_val {
                        pivot_row = i;
                        pivot_val = val;
                    }
                }

                // Swap rows if necessary
                if pivot_row != k {
                    for j in 0..n {
                        let temp = a[(k, j)];
                        a[(k, j)] = a[(pivot_row, j)];
                        a[(pivot_row, j)] = temp;
                    }

                    // Update permutation
                    p.swap(k, pivot_row);
                }

                // Fill in L and U
                for i in k..n {
                    u[(k, i)] = a[(k, i)];
                }

                for i in (k + 1)..n {
                    l[(i, k)] = a[(i, k)] / u[(k, k)];
                    for j in k..n {
                        a[(i, j)] -= l[(i, k)] * u[(k, j)];
                    }
                }
            }
        }

        // Fill in the diagonal of L with 1's
        for i in 0..n {
            l[(i, i)] = Complex64::new(1.0, 0.0);
        }

        // Create identity matrix for right-hand side
        let rhs = Mat::<Complex64>::identity(n, n);

        // Apply permutation to the identity matrix
        let mut p_rhs = Mat::<Complex64>::zeros(n, n);

        // Apply permutation in parallel for large matrices
        if n > 1000 {
            // Compute all permutation updates in parallel
            let updates: Vec<(usize, usize, Complex64)> = (0..n)
                .into_par_iter()
                .flat_map(|i| {
                    let permuted_i = p[i];
                    let mut row_updates = Vec::with_capacity(n);
                    for j in 0..n {
                        row_updates.push((permuted_i, j, rhs[(i, j)]));
                    }
                    row_updates
                })
                .collect();

            // Apply all the updates
            for (i, j, value) in updates {
                p_rhs[(i, j)] = value;
            }
        } else {
            // Sequential processing for smaller matrices
            for i in 0..n {
                for j in 0..n {
                    p_rhs[(p[i], j)] = rhs[(i, j)];
                }
            }
        }

        // Forward substitution to solve L*Y = P*I
        let mut y = Mat::<Complex64>::zeros(n, n);

        // Solve for each right-hand side column in parallel for large matrices
        if n > 500 {
            // We need a different approach for parallel column operations
            // Create a vector to store all column results
            let y_cols: Vec<Vec<Complex64>> = (0..n)
                .into_par_iter()
                .map(|j| {
                    let mut y_col = vec![Complex64::new(0.0, 0.0); n];

                    // Forward substitution for this column
                    for i in 0..n {
                        let mut sum = Complex64::new(0.0, 0.0);
                        for k in 0..i {
                            sum += l[(i, k)] * y_col[k];
                        }
                        y_col[i] = p_rhs[(i, j)] - sum;
                    }

                    // Return the column result
                    y_col
                })
                .collect();

            // Now copy all the results to the y matrix
            for j in 0..n {
                for i in 0..n {
                    y[(i, j)] = y_cols[j][i];
                }
            }
        } else {
            // Sequential processing for smaller matrices
            for j in 0..n {
                for i in 0..n {
                    let mut sum = Complex64::new(0.0, 0.0);
                    for k in 0..i {
                        sum += l[(i, k)] * y[(k, j)];
                    }
                    y[(i, j)] = p_rhs[(i, j)] - sum;
                }
            }
        }

        // Backward substitution to solve U*X = Y
        let mut inverted = Mat::<Complex64>::zeros(n, n);

        // Solve for each right-hand side column in parallel for large matrices
        if n > 500 {
            // We need a different approach for parallel column operations
            // Create a vector to store all column results
            let x_cols: Vec<Vec<Complex64>> = (0..n)
                .into_par_iter()
                .map(|j| {
                    let mut x_col = vec![Complex64::new(0.0, 0.0); n];

                    // Backward substitution for this column
                    for i in (0..n).rev() {
                        let mut sum = Complex64::new(0.0, 0.0);
                        for k in (i + 1)..n {
                            sum += u[(i, k)] * x_col[k];
                        }
                        x_col[i] = (y[(i, j)] - sum) / u[(i, i)];
                    }

                    // Return the column result
                    x_col
                })
                .collect();

            // Now copy all the results to the inverted matrix
            for j in 0..n {
                for i in 0..n {
                    inverted[(i, j)] = x_cols[j][i];
                }
            }
        } else {
            // Sequential processing for smaller matrices
            for j in 0..n {
                for i in (0..n).rev() {
                    let mut sum = Complex64::new(0.0, 0.0);
                    for k in (i + 1)..n {
                        sum += u[(i, k)] * inverted[(k, j)];
                    }
                    inverted[(i, j)] = (y[(i, j)] - sum) / u[(i, i)];
                }
            }
        }

        // Convert back to ndarray
        Ok(faer_to_ndarray(&inverted))
    }

    /// Solve the FMS equations using an iterative method
    ///
    /// This method is suitable for large matrices (thousands of atoms)
    /// where direct decomposition would be too slow or memory-intensive.
    /// It implements the Conjugate Gradient Squared (CGS) method, which is
    /// an iterative method that can converge quickly for well-conditioned systems.
    ///
    /// This implementation uses Rayon for parallel processing of multiple right-hand sides
    /// simultaneously, significantly improving performance for large systems.
    ///
    /// # Arguments
    ///
    /// * `fms_matrix` - The FMS matrix (I - GT)
    ///
    /// # Returns
    ///
    /// The path operator τ
    fn solve_iterative(&self, fms_matrix: &Array2<Complex64>) -> Result<Array2<Complex64>> {
        // Convert to faer matrix for better performance
        let a = ndarray_to_faer(fms_matrix);
        let n = a.nrows();

        // Create a matrix to hold the result
        let mut result = Mat::<Complex64>::zeros(n, n);

        // For large matrices, use Rayon to solve multiple columns in parallel
        // Each column of the result matrix can be calculated independently
        if n > 128 {
            // Solve for each column in parallel
            let column_results: Vec<Result<(usize, col::Col<Complex64>)>> = (0..n)
                .into_par_iter()
                .map(|col_idx| -> Result<(usize, col::Col<Complex64>)> {
                    // Solve the system for this column
                    let column_solution = self.solve_cgs_for_column(&a, col_idx)?;

                    // Return the column index and solution
                    Ok((col_idx, column_solution))
                })
                .collect();

            // Process the column results
            for col_result in column_results {
                match col_result {
                    Ok((col_idx, solution)) => {
                        // Store this column in the result matrix
                        for i in 0..n {
                            result[(i, col_idx)] = solution[i];
                        }
                    }
                    Err(err) => return Err(err),
                }
            }
        } else {
            // For smaller matrices, use sequential processing to avoid overhead
            for col_idx in 0..n {
                // Solve for this column
                let column_solution = self.solve_cgs_for_column(&a, col_idx)?;

                // Store this column in the result matrix
                for i in 0..n {
                    result[(i, col_idx)] = column_solution[i];
                }
            }
        }

        // Convert back to ndarray
        Ok(faer_to_ndarray(&result))
    }

    /// Solve a single column of the FMS matrix using the CGS method
    ///
    /// This helper method solves A*x = b for a single right-hand side vector b.
    /// It extracts the CGS algorithm into a separate function to enable parallel execution.
    ///
    /// # Arguments
    ///
    /// * `a` - The coefficient matrix
    /// * `col_idx` - The index of the column of the identity matrix to use as the right-hand side
    ///
    /// # Returns
    ///
    /// The solution vector x
    fn solve_cgs_for_column(
        &self,
        a: &Mat<Complex64>,
        col_idx: usize,
    ) -> Result<col::Col<Complex64>> {
        let n = a.nrows();

        // Create the right-hand side vector (a column of the identity matrix)
        let mut b = col::Col::<Complex64>::zeros(n);
        b[col_idx] = Complex64::new(1.0, 0.0);

        // Initial guess for the solution
        let mut x = col::Col::<Complex64>::zeros(n);

        // Calculate initial residual r = b - Ax
        let mut r = b.clone();

        // Use blocked matrix-vector multiplication for better cache performance
        let block_size = 32;
        for i_block in (0..n).step_by(block_size) {
            let i_end = (i_block + block_size).min(n);

            for j_block in (0..n).step_by(block_size) {
                let j_end = (j_block + block_size).min(n);

                for i in i_block..i_end {
                    for j in j_block..j_end {
                        r[i] -= a[(i, j)] * x[j];
                    }
                }
            }
        }

        // Choose an arbitrary vector r_tilde (often r_tilde = r)
        let r_tilde = r.clone();

        // Initialize other vectors needed for CGS
        let mut p = r.clone();
        let mut q = col::Col::<Complex64>::zeros(n);
        let mut u = r.clone();

        // Perform iterations
        for iter in 0..self.max_iterations {
            // Check if residual is small enough
            let r_norm = vector_norm(&r);
            if r_norm < self.tolerance {
                break;
            }

            // Calculate rho = (r_tilde, r)
            let rho = inner_product(&r_tilde, &r)?;
            if rho.norm() < 1e-14 {
                return Err(FmsError::IterationFailed(
                    "CGS breakdown (rho ≈ 0)".to_string(),
                ));
            }

            // Calculate q = A*p using blocked approach for better cache performance
            q.fill(Complex64::new(0.0, 0.0));
            for i_block in (0..n).step_by(block_size) {
                let i_end = (i_block + block_size).min(n);

                for j_block in (0..n).step_by(block_size) {
                    let j_end = (j_block + block_size).min(n);

                    for i in i_block..i_end {
                        for j in j_block..j_end {
                            q[i] += a[(i, j)] * p[j];
                        }
                    }
                }
            }

            // Calculate alpha = rho / (r_tilde, q)
            let alpha = rho / inner_product(&r_tilde, &q)?;

            // Update u = r - alpha * q (not u = u + alpha * q as in standard CGS)
            for i in 0..n {
                u[i] = r[i] - alpha * q[i];
            }

            // Calculate A*u using blocked approach for better cache performance
            let mut au = col::Col::<Complex64>::zeros(n);
            for i_block in (0..n).step_by(block_size) {
                let i_end = (i_block + block_size).min(n);

                for j_block in (0..n).step_by(block_size) {
                    let j_end = (j_block + block_size).min(n);

                    for i in i_block..i_end {
                        for j in j_block..j_end {
                            au[i] += a[(i, j)] * u[j];
                        }
                    }
                }
            }

            // Update solution x = x + alpha * (p + u)
            for i in 0..n {
                x[i] += alpha * (p[i] + u[i]);
            }

            // Update residual r = r - alpha * (q + A*u)
            let old_r_norm = r_norm;
            for i in 0..n {
                r[i] -= alpha * (q[i] + au[i]);
            }

            // Check convergence
            let new_r_norm = vector_norm(&r);
            if new_r_norm < self.tolerance {
                break;
            }

            // Check for stagnation (avoid infinite loops)
            if iter > 0 && new_r_norm > old_r_norm * 0.9 {
                // If we're not making progress, try to restart with a new initial guess
                if iter < self.max_iterations / 2 {
                    // Restart the algorithm with the current x
                    // Recalculate the residual r = b - Ax
                    r = b.clone();
                    for i in 0..n {
                        for j in 0..n {
                            r[i] -= a[(i, j)] * x[j];
                        }
                    }

                    // Reset the other vectors
                    p = r.clone();
                    u = r.clone();

                    continue;
                } else {
                    // If we've already tried restarting, report stagnation
                    let stagnation_msg = format!(
                        "CGS stagnation at iteration {}, residual norm: {}",
                        iter, new_r_norm
                    );
                    return Err(FmsError::IterationFailed(stagnation_msg));
                }
            }

            // Calculate beta = (r_tilde, r_new) / (r_tilde, r_old)
            let rho_new = inner_product(&r_tilde, &r)?;
            let beta = rho_new / rho;

            // Update p = r + beta * (p + u)
            for i in 0..n {
                p[i] = r[i] + beta * (p[i] + u[i]);
            }
        }

        Ok(x)
    }

    /// Solve the FMS equations using a block diagonal approximation
    ///
    /// This method is suitable for very large systems where even
    /// iterative methods would be too slow. It divides the matrix into
    /// blocks corresponding to atoms, and solves each block independently.
    /// This approximation is valid when the interaction between distant
    /// atoms is weak.
    ///
    /// This implementation uses Rayon for parallel processing of multiple
    /// diagonal blocks simultaneously, significantly improving performance.
    ///
    /// # Arguments
    ///
    /// * `fms_matrix` - The FMS matrix (I - GT)
    ///
    /// # Returns
    ///
    /// The path operator τ
    fn solve_block_diagonal(&self, fms_matrix: &Array2<Complex64>) -> Result<Array2<Complex64>> {
        let n = fms_matrix.shape()[0];

        // Determine the block size based on angular momentum
        // For l_max = 3, each atom has (l_max+1)^2 = 16 angular momentum components
        let block_size = {
            // Try to infer block size by looking for pattern of matrix elements
            // Typically, blocks correspond to atoms with angular momentum components
            let mut candidate_size = 1;
            for size in 4..=36 {
                // Try common block sizes (4, 9, 16, 25, 36)
                if n % size == 0 {
                    candidate_size = size;
                    break;
                }
            }
            candidate_size
        };

        // If we couldn't determine a good block size, report an error
        if block_size == 1 || n % block_size != 0 {
            return Err(FmsError::InvalidParameter(
                "Could not determine appropriate block size for block diagonal approximation"
                    .to_string(),
            ));
        }

        let num_blocks = n / block_size;

        // Create a matrix for the result
        let mut result = Array2::<Complex64>::zeros((n, n));

        // Process blocks in parallel for better performance if we have enough blocks
        if num_blocks >= 4 {
            // Define a structure to hold the block computation results
            struct BlockResult {
                block_idx: usize,
                block_size: usize,
                inverted_block: Mat<Complex64>,
            }

            // Process blocks in parallel using Rayon
            let block_results: Vec<Result<BlockResult>> = (0..num_blocks)
                .into_par_iter()
                .map(|block_idx| -> Result<BlockResult> {
                    // Extract the current diagonal block
                    let start_idx = block_idx * block_size;

                    // Create a block array for this diagonal block
                    let mut block = Array2::<Complex64>::zeros((block_size, block_size));
                    for i in 0..block_size {
                        for j in 0..block_size {
                            block[(i, j)] = fms_matrix[(start_idx + i, start_idx + j)];
                        }
                    }

                    // Solve this block using LU decomposition
                    let inverted_block = self.invert_block(&block)?;

                    // Return the result for this block
                    Ok(BlockResult {
                        block_idx,
                        block_size,
                        inverted_block,
                    })
                })
                .collect();

            // Process the block results and insert them into the main result matrix
            for block_result in block_results {
                match block_result {
                    Ok(result_data) => {
                        let start_idx = result_data.block_idx * result_data.block_size;
                        let block_result = faer_to_ndarray(&result_data.inverted_block);

                        // Copy the block solution to the appropriate position in the result matrix
                        for i in 0..result_data.block_size {
                            for j in 0..result_data.block_size {
                                result[(start_idx + i, start_idx + j)] = block_result[(i, j)];
                            }
                        }
                    }
                    Err(err) => return Err(err),
                }
            }
        } else {
            // For a small number of blocks, use sequential processing to avoid overhead
            for block_idx in 0..num_blocks {
                // Extract the current diagonal block
                let start_idx = block_idx * block_size;

                let mut block = Array2::<Complex64>::zeros((block_size, block_size));
                for i in 0..block_size {
                    for j in 0..block_size {
                        block[(i, j)] = fms_matrix[(start_idx + i, start_idx + j)];
                    }
                }

                // Solve this block using LU decomposition
                let block_inverted = self.invert_block(&block)?;

                // Convert to ndarray and place in the result matrix
                let block_result = faer_to_ndarray(&block_inverted);

                // Copy the block solution to the appropriate position in the result matrix
                for i in 0..block_size {
                    for j in 0..block_size {
                        result[(start_idx + i, start_idx + j)] = block_result[(i, j)];
                    }
                }
            }
        }

        // For off-diagonal blocks, use an approximation or set to zero
        // In this basic implementation, we set off-diagonal blocks to zero,
        // which means we're neglecting long-range interactions

        Ok(result)
    }

    /// Invert a small block matrix using LU decomposition
    ///
    /// This helper method inverts a small matrix block using LU decomposition.
    /// It is extracted as a separate function to enable parallel execution.
    ///
    /// # Arguments
    ///
    /// * `block` - The block matrix to invert
    ///
    /// # Returns
    ///
    /// The inverted block matrix
    fn invert_block(&self, block: &Array2<Complex64>) -> Result<Mat<Complex64>> {
        // Convert to faer matrix for better performance
        let block_faer = ndarray_to_faer(block);
        let n_block = block_faer.nrows();

        // Perform LU decomposition for this block
        let mut a = block_faer.clone();
        let mut p = (0..n_block).collect::<Vec<usize>>();
        let mut l = Mat::<Complex64>::zeros(n_block, n_block);
        let mut u = Mat::<Complex64>::zeros(n_block, n_block);

        // Standard LU decomposition with partial pivoting and cache optimization
        // We're using a block algorithm with a smaller inner block size to improve cache performance
        let inner_block_size = 4;

        for k_block in (0..n_block).step_by(inner_block_size) {
            let k_end = (k_block + inner_block_size).min(n_block);

            // Process current diagonal block
            for k in k_block..k_end {
                // Find pivot within the current column
                let mut pivot_row = k;
                let mut pivot_val = a[(k, k)].norm();

                for i in (k + 1)..n_block {
                    let val = a[(i, k)].norm();
                    if val > pivot_val {
                        pivot_row = i;
                        pivot_val = val;
                    }
                }

                // Swap rows if necessary
                if pivot_row != k {
                    for j in 0..n_block {
                        let temp = a[(k, j)];
                        a[(k, j)] = a[(pivot_row, j)];
                        a[(pivot_row, j)] = temp;
                    }
                    p.swap(k, pivot_row);
                }

                // Fill in U for this row
                for i in k..n_block {
                    u[(k, i)] = a[(k, i)];
                }

                // Update the current column of L and the remaining submatrix
                for i in (k + 1)..n_block {
                    // Calculate L element
                    l[(i, k)] = a[(i, k)] / u[(k, k)];

                    // Update the row in blocks for better cache performance
                    for j_block in (k..n_block).step_by(inner_block_size) {
                        let j_end = (j_block + inner_block_size).min(n_block);

                        for j in j_block..j_end {
                            a[(i, j)] -= l[(i, k)] * u[(k, j)];
                        }
                    }
                }
            }
        }

        // Fill in the diagonal of L with 1's
        for i in 0..n_block {
            l[(i, i)] = Complex64::new(1.0, 0.0);
        }

        // Now solve for the identity matrix using forward and backward substitution
        let rhs = Mat::<Complex64>::identity(n_block, n_block);

        // Apply permutation to the identity matrix
        let mut p_rhs = Mat::<Complex64>::zeros(n_block, n_block);
        for i in 0..n_block {
            for j in 0..n_block {
                p_rhs[(p[i], j)] = rhs[(i, j)];
            }
        }

        // Forward substitution to solve L*Y = P*I
        let mut y = Mat::<Complex64>::zeros(n_block, n_block);
        // Solve for each right-hand side column
        for j in 0..n_block {
            for i in 0..n_block {
                let mut sum = Complex64::new(0.0, 0.0);
                for k in 0..i {
                    sum += l[(i, k)] * y[(k, j)];
                }
                y[(i, j)] = p_rhs[(i, j)] - sum;
            }
        }

        // Backward substitution to solve U*X = Y
        let mut block_inverted = Mat::<Complex64>::zeros(n_block, n_block);
        // Solve for each right-hand side column
        for j in 0..n_block {
            for i in (0..n_block).rev() {
                let mut sum = Complex64::new(0.0, 0.0);
                for k in (i + 1)..n_block {
                    sum += u[(i, k)] * block_inverted[(k, j)];
                }
                block_inverted[(i, j)] = (y[(i, j)] - sum) / u[(i, i)];
            }
        }

        Ok(block_inverted)
    }
}
