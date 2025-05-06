/*
MIT License with FEFF10 Attribution

Copyright (c) 2025 Ameyanagi

Based on or developed using Distribution: FEFF10.0
Copyright (c) 2020 FEFF Project, University of Washington and SLAC National Accelerator Laboratory.
All rights reserved.
*/

//! FMS matrix builder
//!
//! This module implements the FMS matrix construction for XANES calculations.

use super::errors::{FmsError, Result};
use crate::atoms::AtomicStructure;
use crate::scattering::ScatteringMatrixResults;
use crate::utils::linear_algebra::{faer_to_ndarray, ndarray_to_faer};
use ndarray::Array2;
use num_complex::Complex64;
use std::ops::Mul;
use rayon::prelude::*;

/// FMS matrix types
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum FmsMatrixType {
    /// Full FMS matrix (I - GT)
    Full,
    /// Only the free propagator G
    Propagator,
    /// Only the T-matrix
    TMatrix,
}

/// FMS matrix builder
///
/// Constructs the FMS matrices needed for full multiple scattering calculations.
/// This includes filtering atoms within the FMS radius, building the scattering
/// matrices, and constructing the full matrix (I - GT) needed for FMS calculations.
#[derive(Debug)]
pub struct FmsMatrix<'a> {
    /// Reference to the atomic structure
    structure: &'a AtomicStructure,
    /// FMS radius in Angstroms
    radius: f64,
    /// Indices of atoms within the FMS radius
    fms_atoms: Vec<usize>,
}

impl<'a> FmsMatrix<'a> {
    /// Create a new FMS matrix builder
    ///
    /// # Arguments
    ///
    /// * `structure` - The atomic structure
    /// * `radius` - FMS radius in Angstroms
    ///
    /// # Returns
    ///
    /// A new FMS matrix builder
    pub fn new(structure: &'a AtomicStructure, radius: f64) -> Result<Self> {
        // Validate FMS radius
        if radius <= 0.0 {
            return Err(FmsError::InvalidParameter(format!(
                "FMS radius must be positive: {}",
                radius
            )));
        }

        // Get central atom position
        let central_atom = match structure.central_atom() {
            Some(atom) => atom,
            None => {
                return Err(FmsError::InvalidParameter(
                    "No central atom specified in structure".to_string(),
                ));
            }
        };
        let central_pos = *central_atom.position();

        // Find atoms within the FMS radius
        let fms_atoms: Vec<usize> = (0..structure.atom_count())
            .filter(|&idx| {
                if let Some(atom) = structure.atom(idx) {
                    let pos = atom.position();
                    let dx = pos.x - central_pos.x;
                    let dy = pos.y - central_pos.y;
                    let dz = pos.z - central_pos.z;
                    let distance = (dx * dx + dy * dy + dz * dz).sqrt();
                    distance <= radius
                } else {
                    false
                }
            })
            .collect();

        // Check that we have at least one atom
        if fms_atoms.is_empty() {
            return Err(FmsError::RadiusTooSmall(format!(
                "No atoms found within FMS radius {} Å",
                radius
            )));
        }

        Ok(Self {
            structure,
            radius,
            fms_atoms,
        })
    }

    /// Get the number of atoms within the FMS radius
    pub fn atom_count(&self) -> usize {
        self.fms_atoms.len()
    }

    /// Build the FMS matrix
    ///
    /// # Arguments
    ///
    /// * `scattering_matrices` - Pre-calculated scattering matrices
    ///
    /// # Returns
    ///
    /// The FMS matrix (I - GT) for the calculation
    pub fn build_fms_matrix(
        &self,
        scattering_matrices: &ScatteringMatrixResults,
    ) -> Result<Array2<Complex64>> {
        // Extract matrices
        let green_matrix = &scattering_matrices.green_matrix;
        let global_t_matrix = &scattering_matrices.global_t_matrix;
        let max_l = scattering_matrices.max_l;
        let l_size = ((max_l + 1) * (max_l + 1)) as usize;

        // Check dimensions
        let atom_count = self.atom_count();
        let total_size = atom_count * l_size;

        // Filter matrices to include only atoms within the FMS radius
        let filtered_green_matrix = self.filter_matrix(green_matrix, FmsMatrixType::Propagator)?;
        let filtered_t_matrix = self.filter_matrix(global_t_matrix, FmsMatrixType::TMatrix)?;

        // Convert to faer matrices for better performance
        let g_faer = ndarray_to_faer(&filtered_green_matrix);
        let t_faer = ndarray_to_faer(&filtered_t_matrix);

        // Calculate GT product (critical for performance)
        let gt_faer = g_faer.mul(&t_faer);

        // Create identity matrix
        let mut result_faer = faer::Mat::<Complex64>::identity(total_size, total_size);

        // Calculate I - GT in parallel for better performance using a different approach
        // We'll use a parallel iterator to create the matrix directly
        if total_size > 1000 {
            // Use parallel_chunks method for better performance
            let chunk_size = 64;  // Process 64 rows at a time
            
            // Calculate the matrix in parallel using chunks
            let updates: Vec<(usize, usize, Complex64)> = (0..total_size)
                .into_par_iter()
                .flat_map(|i| {
                    let mut row_updates = Vec::with_capacity(total_size);
                    for j in 0..total_size {
                        row_updates.push((i, j, -gt_faer[(i, j)])); // Store the update value
                    }
                    row_updates
                })
                .collect();
            
            // Apply the updates to the result matrix
            for (i, j, value) in updates {
                result_faer[(i, j)] += value; // Add the update value to identity matrix
            }
        } else {
            // For smaller matrices, use sequential processing to avoid overhead
            for i in 0..total_size {
                for j in 0..total_size {
                    result_faer[(i, j)] -= gt_faer[(i, j)];
                }
            }
        }

        // Convert back to ndarray
        Ok(faer_to_ndarray(&result_faer))
    }

    /// Filter a matrix to include only atoms within the FMS radius
    ///
    /// # Arguments
    ///
    /// * `matrix` - The full matrix to filter
    /// * `matrix_type` - Type of the matrix (propagator, T-matrix, or full)
    ///
    /// # Returns
    ///
    /// Filtered matrix
    fn filter_matrix(
        &self,
        matrix: &Array2<Complex64>,
        matrix_type: FmsMatrixType,
    ) -> Result<Array2<Complex64>> {
        // Get dimensions
        let (rows, cols) = matrix.dim();
        
        // Implement proper filtering to include only atoms within the FMS radius
        // This is a placeholder implementation for now
        let filtered_matrix = match matrix_type {
            FmsMatrixType::Propagator => {
                // For Green's function (propagator), we'll filter based on distance
                // between atoms in the future implementation
                matrix.clone()
            },
            FmsMatrixType::TMatrix => {
                // For T-matrix, we may want to apply different filtering criteria
                matrix.clone()
            },
            FmsMatrixType::Full => {
                // For the full matrix, we'll apply combined filtering
                matrix.clone()
            }
        };
        
        // If the matrix is large, use parallel processing for the clone operation
        if rows * cols > 10000 {
            // Create a parallel version with the same dimensions using a similar approach
            let mut parallel_result = Array2::<Complex64>::zeros((rows, cols));
            
            // Use parallel iterators to create a list of updates
            let updates: Vec<(usize, usize, Complex64)> = (0..rows)
                .into_par_iter()
                .flat_map(|i| {
                    let mut row_updates = Vec::with_capacity(cols);
                    for j in 0..cols {
                        row_updates.push((i, j, filtered_matrix[(i, j)]));
                    }
                    row_updates
                })
                .collect();
            
            // Apply all updates to the result matrix
            for (i, j, value) in updates {
                parallel_result[(i, j)] = value;
            }
            
            return Ok(parallel_result);
        }
        
        // For smaller matrices, just return the cloned result
        Ok(filtered_matrix)
    }

    /// Calculate the matrix elements for the full Green's function
    ///
    /// This function calculates the Green's function matrix elements between
    /// all pairs of atoms within the FMS radius. It uses the pre-calculated
    /// free electron Green's function and adds any necessary corrections.
    ///
    /// # Arguments
    ///
    /// * `scattering_matrices` - Pre-calculated scattering matrices
    ///
    /// # Returns
    ///
    /// The Green's function matrix
    fn calculate_greens_matrix(
        &self,
        scattering_matrices: &ScatteringMatrixResults,
    ) -> Result<Array2<Complex64>> {
        // Get the pre-calculated Green's function matrix
        let green_matrix = &scattering_matrices.green_matrix;
        let (rows, cols) = green_matrix.dim();
        
        // In a future implementation, we'll add curved-wave corrections and other
        // physical effects to the free electron Green's function
        
        // For large matrices, use parallel processing
        if rows * cols > 10000 {
            // Create a new result matrix with the same dimensions
            let mut result = Array2::<Complex64>::zeros((rows, cols));
            
            // Use parallel iterators to create a list of updates
            let updates: Vec<(usize, usize, Complex64)> = (0..rows)
                .into_par_iter()
                .flat_map(|i| {
                    let mut row_updates = Vec::with_capacity(cols);
                    for j in 0..cols {
                        // In the future, we'll add physical corrections here
                        // For now, just copy the original values
                        row_updates.push((i, j, green_matrix[(i, j)]));
                    }
                    row_updates
                })
                .collect();
            
            // Apply all updates to the result matrix
            for (i, j, value) in updates {
                result[(i, j)] = value;
            }
            
            return Ok(result);
        }
        
        // For smaller matrices, just clone the original
        Ok(green_matrix.clone())
    }
}