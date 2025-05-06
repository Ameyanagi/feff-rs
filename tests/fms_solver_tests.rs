/*
MIT License with FEFF10 Attribution

Copyright (c) 2025 Ameyanagi

Based on or developed using Distribution: FEFF10.0
Copyright (c) 2020 FEFF Project, University of Washington and SLAC National Accelerator Laboratory.
All rights reserved.
*/

//! Tests for the FMS solver module

use feff_rs::fms::SolverMethod;
use ndarray::Array2;
use num_complex::Complex64;

// Create a simple test matrix for solver testing
fn create_test_matrix(size: usize) -> Array2<Complex64> {
    let mut matrix = Array2::<Complex64>::zeros((size, size));

    // Diagonal dominant matrix
    for i in 0..size {
        matrix[(i, i)] = Complex64::new(size as f64, 0.0); // Strong diagonal

        // Off-diagonal elements
        for j in 0..size {
            if i != j {
                // Weaker off-diagonal elements to ensure convergence
                matrix[(i, j)] = Complex64::new(0.1, 0.05) * (-1.0_f64).powi((i + j) as i32);
            }
        }
    }

    matrix
}

#[test]
fn test_lu_solver() {
    // Create a small test matrix (8x8)
    let matrix_size = 8;
    let test_matrix = create_test_matrix(matrix_size);

    // Import the solver for testing
    use feff_rs::fms::FmsSolver;

    // Create solvers with different methods
    let lu_solver = FmsSolver::new(SolverMethod::LuDecomposition);

    // Solve using LU decomposition
    let solution = lu_solver.solve(&test_matrix).unwrap();

    // Verify the solution by checking (A * A^-1) ≈ I
    for i in 0..matrix_size {
        for j in 0..matrix_size {
            let mut sum = Complex64::new(0.0, 0.0);
            for k in 0..matrix_size {
                sum += test_matrix[(i, k)] * solution[(k, j)];
            }

            if i == j {
                // Diagonal elements should be close to 1
                assert!((sum - Complex64::new(1.0, 0.0)).norm() < 1e-10);
            } else {
                // Off-diagonal elements should be close to 0
                assert!(sum.norm() < 1e-10);
            }
        }
    }
}

#[test]
fn test_iterative_solver() {
    // Create a small test matrix (8x8)
    let matrix_size = 8;
    let test_matrix = create_test_matrix(matrix_size);

    // Import the solver for testing
    use feff_rs::fms::FmsSolver;

    // Create iterative solver
    let mut iterative_solver = FmsSolver::new(SolverMethod::IterativeCGS);
    iterative_solver.set_tolerance(1e-8).set_max_iterations(100);

    // Solve using iterative method
    let solution = iterative_solver.solve(&test_matrix).unwrap();

    // Verify the solution by checking (A * A^-1) ≈ I with a looser tolerance
    // Iterative methods may not achieve the same precision as direct methods
    for i in 0..matrix_size {
        for j in 0..matrix_size {
            let mut sum = Complex64::new(0.0, 0.0);
            for k in 0..matrix_size {
                sum += test_matrix[(i, k)] * solution[(k, j)];
            }

            if i == j {
                // Diagonal elements should be close to 1
                assert!((sum - Complex64::new(1.0, 0.0)).norm() < 1e-6);
            } else {
                // Off-diagonal elements should be close to 0
                assert!(sum.norm() < 1e-6);
            }
        }
    }
}

#[test]
fn test_block_diagonal_solver() {
    // Create a test matrix with a block structure (16x16 with 4x4 blocks)
    let block_size = 4;
    let num_blocks = 4;
    let matrix_size = block_size * num_blocks;

    // Initialize with zeros
    let mut test_matrix = Array2::<Complex64>::zeros((matrix_size, matrix_size));

    // Fill diagonal blocks with strong dominant matrices
    for block in 0..num_blocks {
        let start = block * block_size;
        let end = start + block_size;

        // Create a dominant diagonal block
        for i in start..end {
            for j in start..end {
                if i == j {
                    test_matrix[(i, j)] = Complex64::new(block_size as f64 * 2.0, 0.0);
                } else {
                    test_matrix[(i, j)] =
                        Complex64::new(0.1, 0.05) * (-1.0_f64).powi((i + j) as i32);
                }
            }
        }

        // Add weak coupling between blocks (so block-diagonal approximation is reasonable)
        if block > 0 {
            let prev_start = (block - 1) * block_size;
            let prev_end = prev_start + block_size;

            // Weak coupling to previous block
            for i in start..end {
                for j in prev_start..prev_end {
                    test_matrix[(i, j)] = Complex64::new(0.01, 0.005);
                    test_matrix[(j, i)] = Complex64::new(0.01, -0.005);
                }
            }
        }
    }

    // Import the solver for testing
    use feff_rs::fms::FmsSolver;

    // Create block diagonal solver
    let block_solver = FmsSolver::new(SolverMethod::BlockDiagonal);

    // Solve using block diagonal method
    let solution = block_solver.solve(&test_matrix).unwrap();

    // For block diagonal approximation, we verify that the diagonal blocks
    // of (A * A^-1) are close to identity blocks
    for block in 0..num_blocks {
        let start = block * block_size;
        let end = start + block_size;

        for i in start..end {
            for j in start..end {
                let mut sum = Complex64::new(0.0, 0.0);
                for k in start..end {
                    // Only check within the block
                    sum += test_matrix[(i, k)] * solution[(k, j)];
                }

                if i == j {
                    // Diagonal elements should be close to 1
                    assert!(
                        (sum - Complex64::new(1.0, 0.0)).norm() < 1e-6,
                        "Block diagonal approximation failed at ({}, {}): expected 1.0, got {:?}",
                        i,
                        j,
                        sum
                    );
                } else {
                    // Off-diagonal elements should be close to 0
                    assert!(
                        sum.norm() < 1e-6,
                        "Block diagonal approximation failed at ({}, {}): expected 0.0, got {:?}",
                        i,
                        j,
                        sum
                    );
                }
            }
        }
    }
}

#[test]
fn test_solver_method_comparison() {
    // Create a small test matrix (8x8)
    let matrix_size = 8;
    let test_matrix = create_test_matrix(matrix_size);

    // Import the solver for testing
    use feff_rs::fms::FmsSolver;

    // Create solvers with different methods
    let lu_solver = FmsSolver::new(SolverMethod::LuDecomposition);
    let iterative_solver = FmsSolver::new(SolverMethod::IterativeCGS);

    // Solve using both methods
    let lu_solution = lu_solver.solve(&test_matrix).unwrap();
    let iterative_solution = iterative_solver.solve(&test_matrix).unwrap();

    // Compare solutions (they should be close)
    for i in 0..matrix_size {
        for j in 0..matrix_size {
            // Allow for some numerical differences due to different methods
            assert!(
                (lu_solution[(i, j)] - iterative_solution[(i, j)]).norm() < 1e-5,
                "Solutions differ at ({}, {}): LU={:?}, Iterative={:?}",
                i,
                j,
                lu_solution[(i, j)],
                iterative_solution[(i, j)]
            );
        }
    }
}
