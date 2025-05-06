/*
MIT License with FEFF10 Attribution

Copyright (c) 2025 Ameyanagi

Based on or developed using Distribution: FEFF10.0
Copyright (c) 2020 FEFF Project, University of Washington and SLAC National Accelerator Laboratory.
All rights reserved.
*/

//! Linear algebra utilities using the Faer library
//!
//! This module provides basic linear algebra operations
//! for FEFF calculations using the Faer library.

use super::errors::{Result, UtilsError};
use faer::{col, Mat};
use ndarray::{Array1, Array2};
use num_complex::Complex64;
use rayon::prelude::*;

/// Convert from ndarray::Array2<Complex64> to faer::Mat<Complex64>
///
/// This function efficiently converts an ndarray matrix to a Faer matrix,
/// optimized for SIMD operations. Faer matrices are better suited for
/// high-performance linear algebra operations.
pub fn ndarray_to_faer(array: &Array2<Complex64>) -> Mat<Complex64> {
    let (rows, cols) = array.dim();

    // Create matrix with uninitialized memory for better performance
    let mut result = Mat::<Complex64>::zeros(rows, cols);

    // Use row-major iteration for better cache locality
    if rows * cols > 1000 {
        // For large matrices, use parallel processing
        let chunks = rows.min(16); // Process up to 16 rows in parallel
        let chunk_size = rows.div_ceil(chunks);

        // Process chunks of rows in parallel using Rayon
        let rows_chunks: Vec<_> = (0..chunks)
            .map(|chunk_idx| {
                let start_row = chunk_idx * chunk_size;
                let end_row = (start_row + chunk_size).min(rows);
                (start_row, end_row)
            })
            .collect();

        // Process in parallel with a mutable closure
        let results: Vec<_> = rows_chunks
            .into_par_iter()
            .map(|(start_row, end_row)| {
                // Create a partial result for this chunk
                let mut partial = vec![vec![Complex64::new(0.0, 0.0); cols]; end_row - start_row];

                // Fill the partial result
                for (local_i, i) in (start_row..end_row).enumerate() {
                    for j in 0..cols {
                        partial[local_i][j] = array[(i, j)];
                    }
                }

                (start_row, partial)
            })
            .collect();

        // Combine the partial results
        for (start_row, partial) in results {
            for (local_i, i) in (start_row..(start_row + partial.len())).enumerate() {
                for j in 0..cols {
                    result[(i, j)] = partial[local_i][j];
                }
            }
        }
    } else {
        // For smaller matrices, use sequential processing to avoid overhead
        for i in 0..rows {
            for j in 0..cols {
                result[(i, j)] = array[(i, j)];
            }
        }
    }

    result
}

/// Convert from faer::Mat<Complex64> to ndarray::Array2<Complex64>
///
/// This function efficiently converts a Faer matrix back to an ndarray matrix.
/// It's optimized for copying large blocks of data with good cache locality.
pub fn faer_to_ndarray(matrix: &Mat<Complex64>) -> Array2<Complex64> {
    let rows = matrix.nrows();
    let cols = matrix.ncols();

    // Create output matrix
    let mut result = Array2::<Complex64>::zeros((rows, cols));

    // Use row-major iteration for better cache locality
    if rows * cols > 1000 {
        // For large matrices, use parallel processing
        let chunks = rows.min(16); // Process up to 16 rows in parallel
        let chunk_size = rows.div_ceil(chunks);

        // Process chunks of rows in parallel using Rayon
        let rows_chunks: Vec<_> = (0..chunks)
            .map(|chunk_idx| {
                let start_row = chunk_idx * chunk_size;
                let end_row = (start_row + chunk_size).min(rows);
                (start_row, end_row)
            })
            .collect();

        // Process in parallel with a mutable closure
        let results: Vec<_> = rows_chunks
            .into_par_iter()
            .map(|(start_row, end_row)| {
                // Create a partial result for this chunk
                let mut partial = vec![vec![Complex64::new(0.0, 0.0); cols]; end_row - start_row];

                // Fill the partial result
                for (local_i, i) in (start_row..end_row).enumerate() {
                    for j in 0..cols {
                        partial[local_i][j] = matrix[(i, j)];
                    }
                }

                (start_row, partial)
            })
            .collect();

        // Combine the partial results
        for (start_row, partial) in results {
            for (local_i, i) in (start_row..(start_row + partial.len())).enumerate() {
                for j in 0..cols {
                    result[(i, j)] = partial[local_i][j];
                }
            }
        }
    } else {
        // For smaller matrices, use sequential processing to avoid overhead
        for i in 0..rows {
            for j in 0..cols {
                result[(i, j)] = matrix[(i, j)];
            }
        }
    }

    result
}

/// Create a zero-filled complex matrix with Faer
pub fn zeros(rows: usize, cols: usize) -> Mat<Complex64> {
    Mat::<Complex64>::zeros(rows, cols)
}

/// Create an identity complex matrix with Faer
pub fn eye(size: usize) -> Mat<Complex64> {
    Mat::<Complex64>::identity(size, size)
}

/// Convert from ndarray::Array1<Complex64> to faer::col::Col<Complex64>
pub fn ndarray_to_faer_vector(array: &Array1<Complex64>) -> col::Col<Complex64> {
    let n = array.len();
    let mut result = col::Col::<Complex64>::zeros(n);

    for i in 0..n {
        result[i] = array[i];
    }

    result
}

/// Convert from faer::col::Col<Complex64> to ndarray::Array1<Complex64>
pub fn faer_vector_to_ndarray(vector: &col::Col<Complex64>) -> Array1<Complex64> {
    let n = vector.nrows();
    let mut result = Array1::<Complex64>::zeros(n);

    for i in 0..n {
        result[i] = vector[i];
    }

    result
}

/// Create a zero-filled complex vector with Faer
pub fn zeros_vector(size: usize) -> col::Col<Complex64> {
    col::Col::<Complex64>::zeros(size)
}

/// Compute the dot product of two complex vectors with SIMD optimization
///
/// Uses Faer's optimized SIMD operations when available for better performance.
pub fn dot_product(a: &col::Col<Complex64>, b: &col::Col<Complex64>) -> Result<Complex64> {
    if a.nrows() != b.nrows() {
        return Err(UtilsError::Generic(
            "Vectors must have the same length for dot product".to_string(),
        ));
    }

    let n = a.nrows();

    // For small vectors, use simple sequential calculation
    if n < 16 {
        let mut result = Complex64::new(0.0, 0.0);
        for i in 0..n {
            result += a[i] * b[i];
        }
        return Ok(result);
    }

    // For larger vectors, use Faer's optimized operations
    // Create matrices from vectors for using matrix multiplication
    let a_mat = a.as_ref();
    let b_mat = b.transpose();

    // Use matrix-matrix multiplication, which is optimized for SIMD operations in Faer
    // This performs a * b^T, which gives us the dot product as a 1x1 matrix
    let result_mat = a_mat * b_mat;

    // Extract the single element from the result matrix
    Ok(result_mat[(0, 0)])
}

/// Compute the conjugate dot product (inner product) of two complex vectors
/// with SIMD optimization where possible
///
/// This implementation intelligently adapts to vector size, using SIMD
/// instructions for larger vectors where the overhead is worthwhile.
pub fn inner_product(a: &col::Col<Complex64>, b: &col::Col<Complex64>) -> Result<Complex64> {
    if a.nrows() != b.nrows() {
        return Err(UtilsError::Generic(
            "Vectors must have the same length for inner product".to_string(),
        ));
    }

    let n = a.nrows();

    // For small vectors, use simple sequential calculation
    if n < 16 {
        let mut result = Complex64::new(0.0, 0.0);
        for i in 0..n {
            result += a[i].conj() * b[i];
        }
        return Ok(result);
    }

    // For larger vectors, use SIMD-optimized approach
    // We'll compute this in parallel chunks to benefit from SIMD

    // Create conjugated vector a_conj
    let mut a_conj = col::Col::<Complex64>::zeros(n);

    // Compute conjugate in chunks to leverage cache locality
    // This is especially important for large vectors
    const CHUNK_SIZE: usize = 64; // Typical cache line size

    if n > 1000 {
        // For very large vectors, use parallel processing with Rayon
        let chunks = n.div_ceil(CHUNK_SIZE);

        // Process chunks in parallel
        let partial_sums: Vec<Complex64> = (0..chunks)
            .into_par_iter()
            .map(|chunk| {
                let start = chunk * CHUNK_SIZE;
                let end = (start + CHUNK_SIZE).min(n);

                // Compute partial sum for this chunk
                let mut sum = Complex64::new(0.0, 0.0);
                for i in start..end {
                    sum += a[i].conj() * b[i];
                }

                sum
            })
            .collect();

        // Sum up all partial results
        let result = partial_sums
            .into_iter()
            .fold(Complex64::new(0.0, 0.0), |acc, val| acc + val);

        Ok(result)
    } else {
        // For moderately sized vectors, use Faer's optimized operations
        // First, create a conjugated version of a
        for i in 0..n {
            a_conj[i] = a[i].conj();
        }

        // Convert to matrices for matrix multiplication
        let a_conj_mat = a_conj.as_ref();
        let b_mat = b.transpose();

        // Use matrix multiplication, which is optimized for SIMD
        let result_mat = a_conj_mat * b_mat;

        // Extract the single element
        Ok(result_mat[(0, 0)])
    }
}

/// Compute the vector norm (Euclidean) with SIMD optimization
///
/// This function calculates the Euclidean norm of a complex vector,
/// using SIMD instructions for larger vectors when available.
pub fn vector_norm(v: &col::Col<Complex64>) -> f64 {
    let n = v.nrows();

    // For small vectors, use simple sequential calculation
    if n < 16 {
        let mut sum_squares: f64 = 0.0;
        for item in v.iter().take(n) {
            sum_squares += item.norm_sqr();
        }
        return sum_squares.sqrt();
    }

    // For larger vectors, use Faer's optimized matrix operations
    // Calculate the norm by explicitly computing the dot product

    // For vectors of moderate size, use explicit calculation with chunks for better cache behavior
    if n < 1000 {
        // Create chunks to improve cache locality
        let chunk_size = 64; // Typical cache line size
        let mut sum_squares = 0.0;

        // Process vector in chunks
        for chunk_start in (0..n).step_by(chunk_size) {
            let chunk_end = (chunk_start + chunk_size).min(n);
            let mut chunk_sum = 0.0;

            // Process this chunk
            for i in chunk_start..chunk_end {
                chunk_sum += v[i].norm_sqr();
            }

            sum_squares += chunk_sum;
        }

        return sum_squares.sqrt();
    }

    // For very large vectors, use parallel computation
    let chunk_size = 128;
    let chunks = n.div_ceil(chunk_size);

    // Process chunks in parallel
    let partial_sums: Vec<f64> = (0..chunks)
        .into_par_iter()
        .map(|chunk| {
            let start = chunk * chunk_size;
            let end = (start + chunk_size).min(n);

            // Compute partial sum for this chunk
            let mut sum = 0.0;
            for i in start..end {
                sum += v[i].norm_sqr();
            }

            sum
        })
        .collect();

    // Sum all partial results
    let sum_squares = partial_sums.iter().sum::<f64>();

    // Take the square root
    sum_squares.sqrt()
}

/// Normalize a vector to unit length with SIMD optimization
///
/// This function normalizes a complex vector to unit length,
/// using SIMD acceleration for better performance on larger vectors.
pub fn normalize_vector(v: &col::Col<Complex64>) -> col::Col<Complex64> {
    let n = v.nrows();
    let norm = vector_norm(v);

    // Return original vector if norm is zero to avoid division by zero
    if norm < 1e-14 {
        return v.clone();
    }

    let inv_norm = 1.0 / norm;
    let mut result = col::Col::<Complex64>::zeros(n);

    // For small vectors, use simple sequential division
    if n < 32 {
        for i in 0..n {
            result[i] = v[i] * Complex64::new(inv_norm, 0.0);
        }
        return result;
    }

    // For moderate-sized vectors, process in chunks for better cache behavior
    if n < 1000 {
        let chunk_size = 64; // Typical cache line size

        // Process vector in chunks
        for chunk_start in (0..n).step_by(chunk_size) {
            let chunk_end = (chunk_start + chunk_size).min(n);

            // Process this chunk
            for i in chunk_start..chunk_end {
                result[i] = v[i] * Complex64::new(inv_norm, 0.0);
            }
        }

        return result;
    }

    // For very large vectors, use parallel computation
    let chunk_size = 128;
    let chunks = n.div_ceil(chunk_size);

    // Create chunks for parallel processing
    let chunks_vec: Vec<_> = (0..chunks)
        .map(|chunk| {
            let start = chunk * chunk_size;
            let end = (start + chunk_size).min(n);
            (start, end)
        })
        .collect();

    // Process chunks in parallel
    let results: Vec<_> = chunks_vec
        .into_par_iter()
        .map(|(start, end)| {
            // Create a partial result for this chunk
            let mut partial = vec![Complex64::new(0.0, 0.0); end - start];

            // Fill the partial result
            for (local_i, i) in (start..end).enumerate() {
                partial[local_i] = v[i] * Complex64::new(inv_norm, 0.0);
            }

            (start, partial)
        })
        .collect();

    // Combine the partial results
    for (start, partial) in results {
        for (local_i, i) in (start..(start + partial.len())).enumerate() {
            result[i] = partial[local_i];
        }
    }

    result
}

/// Compute the outer product of two complex vectors: v⊗w* with SIMD optimization
///
/// This function computes the outer product of two complex vectors,
/// using SIMD instructions on compatible hardware for better performance.
pub fn outer_product(v: &col::Col<Complex64>, w: &col::Col<Complex64>) -> Mat<Complex64> {
    let m = v.nrows();
    let n = w.nrows();

    // Special optimization for empty or singleton vectors
    if m == 0 || n == 0 {
        return Mat::<Complex64>::zeros(m, n);
    }

    // Create matrix for result
    let mut result = Mat::<Complex64>::zeros(m, n);

    // For small vectors, use direct computation
    if m * n < 256 {
        for i in 0..m {
            for j in 0..n {
                result[(i, j)] = v[i] * w[j].conj();
            }
        }
        return result;
    }

    // For larger vectors, use optimized block operations with parallelism

    // Create conjugated w vector
    let mut w_conj = col::Col::<Complex64>::zeros(n);
    for j in 0..n {
        w_conj[j] = w[j].conj();
    }

    // Use chunking for better cache behavior
    let row_chunk_size = 32.min(m);
    let col_chunk_size = 32.min(n);

    // Process in chunks
    for i_chunk in (0..m).step_by(row_chunk_size) {
        let i_end = (i_chunk + row_chunk_size).min(m);

        for j_chunk in (0..n).step_by(col_chunk_size) {
            let j_end = (j_chunk + col_chunk_size).min(n);

            // Process this block
            for i in i_chunk..i_end {
                for j in j_chunk..j_end {
                    result[(i, j)] = v[i] * w_conj[j];
                }
            }
        }
    }

    // For very large matrices, use parallel computation for better performance
    if m * n > 10000 {
        // Create a new, larger result matrix
        let mut parallel_result = Mat::<Complex64>::zeros(m, n);

        // Divide work into chunks for parallel processing
        let chunks = 16.min(m);
        let chunk_size = m.div_ceil(chunks);

        // Create chunks for parallel processing
        let chunks_vec: Vec<_> = (0..chunks)
            .map(|chunk| {
                let start_row = chunk * chunk_size;
                let end_row = (start_row + chunk_size).min(m);
                (start_row, end_row)
            })
            .collect();

        // Process chunks in parallel
        let results: Vec<_> = chunks_vec
            .into_par_iter()
            .map(|(start_row, end_row)| {
                // Create a partial result for this chunk
                let mut partial = vec![vec![Complex64::new(0.0, 0.0); n]; end_row - start_row];

                // Fill the partial result
                for (local_i, i) in (start_row..end_row).enumerate() {
                    for j in 0..n {
                        partial[local_i][j] = v[i] * w_conj[j];
                    }
                }

                (start_row, partial)
            })
            .collect();

        // Combine the partial results into the final matrix
        for (start_row, partial) in results {
            for (local_i, i) in (start_row..(start_row + partial.len())).enumerate() {
                for j in 0..n {
                    parallel_result[(i, j)] = partial[local_i][j];
                }
            }
        }

        // Return the parallel result
        return parallel_result;
    }

    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_ndarray_faer_conversion() {
        // Create ndarray matrix
        let mut arr = Array2::<Complex64>::zeros((2, 2));
        arr[(0, 0)] = Complex64::new(1.0, 0.0);
        arr[(0, 1)] = Complex64::new(2.0, 1.0);
        arr[(1, 0)] = Complex64::new(3.0, -1.0);
        arr[(1, 1)] = Complex64::new(4.0, 0.0);

        // Convert to Faer
        let faer_mat = ndarray_to_faer(&arr);

        // Check dimensions and values
        assert_eq!(faer_mat.nrows(), 2);
        assert_eq!(faer_mat.ncols(), 2);
        assert_eq!(faer_mat[(0, 0)].re, 1.0);
        assert_eq!(faer_mat[(0, 0)].im, 0.0);
        assert_eq!(faer_mat[(0, 1)].re, 2.0);
        assert_eq!(faer_mat[(0, 1)].im, 1.0);

        // Convert back to ndarray
        let arr2 = faer_to_ndarray(&faer_mat);

        // Check roundtrip conversion
        assert_eq!(arr[(0, 0)], arr2[(0, 0)]);
        assert_eq!(arr[(0, 1)], arr2[(0, 1)]);
        assert_eq!(arr[(1, 0)], arr2[(1, 0)]);
        assert_eq!(arr[(1, 1)], arr2[(1, 1)]);
    }

    #[test]
    fn test_matrix_creation() {
        // Test zeros
        let zero_mat = zeros(2, 3);
        assert_eq!(zero_mat.nrows(), 2);
        assert_eq!(zero_mat.ncols(), 3);
        for i in 0..2 {
            for j in 0..3 {
                assert_eq!(zero_mat[(i, j)].re, 0.0);
                assert_eq!(zero_mat[(i, j)].im, 0.0);
            }
        }

        // Test identity
        let eye_mat = eye(3);
        assert_eq!(eye_mat.nrows(), 3);
        assert_eq!(eye_mat.ncols(), 3);
        for i in 0..3 {
            for j in 0..3 {
                if i == j {
                    assert_eq!(eye_mat[(i, j)].re, 1.0);
                    assert_eq!(eye_mat[(i, j)].im, 0.0);
                } else {
                    assert_eq!(eye_mat[(i, j)].re, 0.0);
                    assert_eq!(eye_mat[(i, j)].im, 0.0);
                }
            }
        }
    }

    #[test]
    fn test_vector_operations() {
        // Create test vectors
        let mut a = col::Col::<Complex64>::zeros(3);
        a[0] = Complex64::new(1.0, 0.0);
        a[1] = Complex64::new(0.0, 1.0);
        a[2] = Complex64::new(2.0, 0.0);

        let mut b = col::Col::<Complex64>::zeros(3);
        b[0] = Complex64::new(2.0, 0.0);
        b[1] = Complex64::new(0.0, -1.0);
        b[2] = Complex64::new(1.0, 0.0);

        // Test dot product: (1)(2) + (i)(-i) + (2)(1) = 2 + 1 + 2 = 5
        let dot = dot_product(&a, &b).unwrap();
        assert_relative_eq!(dot.re, 5.0, epsilon = 1e-10);
        assert_relative_eq!(dot.im, 0.0, epsilon = 1e-10);

        // Test inner product (conjugate dot product):
        // (1*)(2) + (i*)(-i) + (2*)(1) = 2 + (-i)(-i) + 2 = 2 - 1 + 2 = 3
        let inner = inner_product(&a, &b).unwrap();
        assert_relative_eq!(inner.re, 3.0, epsilon = 1e-10);
        assert_relative_eq!(inner.im, 0.0, epsilon = 1e-10);

        // Test vector norm - |a|² = |1|² + |i|² + |2|² = 1 + 1 + 4 = 6
        let norm_a = vector_norm(&a);
        assert_relative_eq!(norm_a, 6.0_f64.sqrt(), epsilon = 1e-10);
    }

    #[test]
    fn test_outer_product() {
        // Create two vectors
        let mut a = col::Col::<Complex64>::zeros(2);
        a[0] = Complex64::new(1.0, 0.0);
        a[1] = Complex64::new(0.0, 1.0);

        let mut b = col::Col::<Complex64>::zeros(2);
        b[0] = Complex64::new(2.0, 0.0);
        b[1] = Complex64::new(0.0, 2.0);

        // Compute the outer product
        let outer = outer_product(&a, &b);

        // The result should be a 2x2 matrix
        assert_eq!(outer.nrows(), 2);
        assert_eq!(outer.ncols(), 2);

        // Expected result:
        // [1*2, 1*conj(0+2i)]     [2, 0-2i]
        // [(0+i)*2, (0+i)*conj(0+2i)] = [0+2i, 2]
        assert_relative_eq!(outer[(0, 0)].re, 2.0, epsilon = 1e-10);
        assert_relative_eq!(outer[(0, 0)].im, 0.0, epsilon = 1e-10);

        assert_relative_eq!(outer[(0, 1)].re, 0.0, epsilon = 1e-10);
        assert_relative_eq!(outer[(0, 1)].im, -2.0, epsilon = 1e-10);

        assert_relative_eq!(outer[(1, 0)].re, 0.0, epsilon = 1e-10);
        assert_relative_eq!(outer[(1, 0)].im, 2.0, epsilon = 1e-10);

        assert_relative_eq!(outer[(1, 1)].re, 2.0, epsilon = 1e-10);
        assert_relative_eq!(outer[(1, 1)].im, 0.0, epsilon = 1e-10);
    }
}
