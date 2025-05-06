/*
MIT License with FEFF10 Attribution

Copyright (c) 2025 Ameyanagi

Based on or developed using Distribution: FEFF10.0
Copyright (c) 2020 FEFF Project, University of Washington and SLAC National Accelerator Laboratory.
All rights reserved.
*/

//! Full Multiple Scattering (FMS) module
//!
//! This module handles full multiple scattering calculations for FEFF.
//! It calculates the multiple scattering path operator and XANES spectra
//! by solving the matrix equation (I - GT)^-1 T for the entire cluster.

// Module structure
mod errors;
mod matrix;
mod solver;
mod xanes;

// Public exports
pub use errors::{FmsError, Result};
pub use matrix::FmsMatrix;
pub use solver::{FmsSolver, SolverMethod};
pub use xanes::XanesCalculator;

use crate::atoms::AtomicStructure;
use crate::scattering::ScatteringMatrixResults;
use ndarray::Array2;
use num_complex::Complex64;
use rayon::prelude::*;

/// FMS calculation parameters
#[derive(Debug, Clone)]
pub struct FmsParameters {
    /// FMS radius in Angstroms
    pub radius: f64,
    /// Method to use for solving the FMS equations
    pub solver_method: SolverMethod,
    /// Energy grid for calculation
    pub energies: Vec<f64>,
    /// Whether to calculate full XANES spectrum
    pub calculate_xanes: bool,
    /// Whether to include thermal effects
    pub include_thermal_effects: bool,
    /// Core-hole lifetime broadening in eV
    pub core_hole_lifetime: f64,
    /// Optional energy shift in eV
    pub energy_shift: Option<f64>,
}

impl Default for FmsParameters {
    fn default() -> Self {
        Self {
            radius: 6.0, // Default 6.0 Å FMS radius
            solver_method: SolverMethod::LuDecomposition,
            energies: Vec::new(),
            calculate_xanes: true,
            include_thermal_effects: false,
            core_hole_lifetime: 1.0, // Default 1.0 eV broadening
            energy_shift: None,
        }
    }
}

/// Results from an FMS calculation
#[derive(Debug)]
pub struct FmsResults {
    /// Energy grid in eV
    pub energies: Vec<f64>,
    /// XANES spectrum (absorption coefficient μ(E))
    pub xanes: Option<Vec<f64>>,
    /// Imaginary part of the scattering amplitude
    pub im_scattering_amplitude: Option<Vec<f64>>,
    /// Real part of the scattering amplitude
    pub re_scattering_amplitude: Option<Vec<f64>>,
    /// Number of atoms in the calculation
    pub atom_count: usize,
    /// Calculation parameters
    pub parameters: FmsParameters,
}

/// Perform full multiple scattering calculation
///
/// This is the main entry point for FMS calculations. It takes pre-calculated
/// scattering matrices and parameters, then calculates the FMS path operator
/// and resulting spectra. This implementation uses parallel processing with Rayon
/// to efficiently calculate spectra for multiple energy points simultaneously.
///
/// # Arguments
///
/// * `structure` - The atomic structure
/// * `scattering_matrices` - Pre-calculated scattering matrices
/// * `parameters` - FMS calculation parameters
///
/// # Returns
///
/// FMS calculation results including spectra
pub fn calculate_fms(
    structure: &AtomicStructure,
    scattering_matrices: &ScatteringMatrixResults,
    parameters: FmsParameters,
) -> Result<FmsResults> {
    // Create FMS matrix builder
    let fms_matrix = FmsMatrix::new(structure, parameters.radius)?;

    // Get or generate energy grid
    let energies = if !parameters.energies.is_empty() {
        parameters.energies.clone()
    } else {
        // Generate default energy grid if none provided
        let e0 = scattering_matrices.energy;
        let emin = e0 - 10.0;
        let emax = e0 + 40.0;
        let estep = 0.5;

        // Generate energy grid directly
        let num_points = ((emax - emin) / estep).ceil() as usize + 1;
        (0..num_points).map(|i| emin + i as f64 * estep).collect()
    };

    // Create FMS solver with appropriate method
    let solver = FmsSolver::new(parameters.solver_method);

    // Create XANES calculator if needed
    let calculate_xanes = parameters.calculate_xanes;
    let xanes_calculator = if calculate_xanes {
        Some(XanesCalculator::new(
            structure,
            parameters.core_hole_lifetime,
        ))
    } else {
        None
    };

    // For parallel calculation, we need thread-safe data structures
    // Create a structure to hold results for each energy point
    #[derive(Debug)]
    struct EnergyResult {
        energy: f64,
        xanes_value: Option<f64>,
        amplitude: Complex64,
    }

    // Determine if we should use parallel processing
    // Only parallelize if we have a significant number of energy points
    let results = if energies.len() > 8 {
        // Before using per-energy point calculation, check if we need the XANES spectrum
        // If we do, we can optimize by doing the XANES calculation for all energies at once
        if calculate_xanes && xanes_calculator.is_some() {
            // First, calculate all the path operators for each energy in parallel
            let path_operators: Vec<Result<(f64, Array2<Complex64>, Complex64)>> = energies
                .par_iter()
                .map(|&energy| -> Result<(f64, Array2<Complex64>, Complex64)> {
                    // Build FMS matrix for this energy
                    let fms_mat = fms_matrix.build_fms_matrix(scattering_matrices)?;

                    // Solve FMS equations
                    let path_operator = solver.solve(&fms_mat)?;

                    // Calculate scattering amplitude
                    let amplitude = calculate_scattering_amplitude(structure, &path_operator)?;

                    // Return the energy, path operator, and amplitude
                    Ok((energy, path_operator, amplitude))
                })
                .collect();

            // Check if any calculation failed
            let mut processed_results = Vec::with_capacity(energies.len());
            let mut operators_by_energy = Vec::with_capacity(energies.len());

            for result in path_operators {
                match result {
                    Ok((energy, path_operator, amplitude)) => {
                        processed_results.push(EnergyResult {
                            energy,
                            xanes_value: None, // We'll fill this in later
                            amplitude,
                        });
                        operators_by_energy.push((energy, path_operator));
                    }
                    Err(err) => return Err(err),
                }
            }

            // Now calculate all XANES values at once using the batch method
            // Extract all path operators and calculate XANES for all energies in one go
            if let Some(ref xanes_calc) = xanes_calculator {
                // For each (energy, path_operator) pair, calculate the XANES value
                for (index, (energy, ref path_operator)) in operators_by_energy.iter().enumerate() {
                    let xanes_value = xanes_calc.calculate_xanes(*energy, path_operator)?;
                    processed_results[index].xanes_value = Some(xanes_value);
                }
            }

            processed_results
        } else {
            // No XANES calculation needed, use the regular parallel approach
            // Process energy points in parallel using Rayon
            energies
                .par_iter()
                .map(|&energy| -> Result<EnergyResult> {
                    // Build FMS matrix for this energy
                    let fms_mat = fms_matrix.build_fms_matrix(scattering_matrices)?;

                    // Solve FMS equations
                    let path_operator = solver.solve(&fms_mat)?;

                    // Calculate scattering amplitude
                    let amplitude = calculate_scattering_amplitude(structure, &path_operator)?;

                    // Calculate XANES if requested
                    let xanes_value = if let Some(ref xanes_calc) = xanes_calculator {
                        if calculate_xanes {
                            Some(xanes_calc.calculate_xanes(energy, &path_operator)?)
                        } else {
                            None
                        }
                    } else {
                        None
                    };

                    // Return results for this energy point
                    Ok(EnergyResult {
                        energy,
                        xanes_value,
                        amplitude,
                    })
                })
                .collect::<Result<Vec<_>>>()?
        }
    } else {
        // For small number of energy points, use sequential processing to avoid overhead
        if calculate_xanes && xanes_calculator.is_some() {
            // Even for small number of points, we can still benefit from optimizing the XANES calculation
            let mut operators_by_energy = Vec::with_capacity(energies.len());
            let mut processed_results = Vec::with_capacity(energies.len());

            // First collect all path operators
            for &energy in energies.iter() {
                // Build FMS matrix for this energy
                let fms_mat = fms_matrix.build_fms_matrix(scattering_matrices)?;

                // Solve FMS equations
                let path_operator = solver.solve(&fms_mat)?;

                // Calculate scattering amplitude
                let amplitude = calculate_scattering_amplitude(structure, &path_operator)?;

                // Add to results
                processed_results.push(EnergyResult {
                    energy,
                    xanes_value: None, // We'll fill this in later
                    amplitude,
                });

                // Save the path operator for XANES calculation
                operators_by_energy.push((energy, path_operator));
            }

            // Now calculate XANES for all energies
            if let Some(ref xanes_calc) = xanes_calculator {
                // For each (energy, path_operator) pair, calculate the XANES value
                for (index, (energy, ref path_operator)) in operators_by_energy.iter().enumerate() {
                    let xanes_value = xanes_calc.calculate_xanes(*energy, path_operator)?;
                    processed_results[index].xanes_value = Some(xanes_value);
                }
            }

            processed_results
        } else {
            // Standard sequential processing
            energies
                .iter()
                .map(|&energy| -> Result<EnergyResult> {
                    // Build FMS matrix for this energy
                    let fms_mat = fms_matrix.build_fms_matrix(scattering_matrices)?;

                    // Solve FMS equations
                    let path_operator = solver.solve(&fms_mat)?;

                    // Calculate scattering amplitude
                    let amplitude = calculate_scattering_amplitude(structure, &path_operator)?;

                    // Calculate XANES if requested
                    let xanes_value = if let Some(ref xanes_calc) = xanes_calculator {
                        if calculate_xanes {
                            Some(xanes_calc.calculate_xanes(energy, &path_operator)?)
                        } else {
                            None
                        }
                    } else {
                        None
                    };

                    // Return results for this energy point
                    Ok(EnergyResult {
                        energy,
                        xanes_value,
                        amplitude,
                    })
                })
                .collect::<Result<Vec<_>>>()?
        }
    };

    // Sort results by energy (parallel processing might return results out of order)
    // We need stable sorting to preserve the original order for identical energy points
    let mut sorted_results = results;
    sorted_results.sort_by(|a, b| {
        a.energy
            .partial_cmp(&b.energy)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    // Combine results into the final data structures
    let mut xanes = if calculate_xanes {
        Some(Vec::with_capacity(energies.len()))
    } else {
        None
    };

    let mut im_scattering_amplitude = Some(Vec::with_capacity(energies.len()));
    let mut re_scattering_amplitude = Some(Vec::with_capacity(energies.len()));

    // Extract results from sorted_results
    for result in sorted_results {
        // Store scattering amplitude components
        if let Some(ref mut im_amp) = im_scattering_amplitude {
            im_amp.push(result.amplitude.im);
        }
        if let Some(ref mut re_amp) = re_scattering_amplitude {
            re_amp.push(result.amplitude.re);
        }

        // Store XANES values if calculated
        if let Some(ref mut xanes_spectrum) = xanes {
            if let Some(xanes_value) = result.xanes_value {
                xanes_spectrum.push(xanes_value);
            }
        }
    }

    // Return combined results
    Ok(FmsResults {
        energies,
        xanes,
        im_scattering_amplitude,
        re_scattering_amplitude,
        atom_count: structure.atom_count(),
        parameters,
    })
}

/// Calculate the scattering amplitude from the path operator
///
/// This function computes the physical scattering amplitude by taking the trace
/// of the appropriate blocks of the path operator matrix. This represents the
/// sum over all scattering paths that contribute to the spectrum.
///
/// # Arguments
///
/// * `structure` - The atomic structure
/// * `path_operator` - The calculated path operator matrix
///
/// # Returns
///
/// Complex scattering amplitude
fn calculate_scattering_amplitude(
    structure: &AtomicStructure,
    path_operator: &ndarray::Array2<Complex64>,
) -> Result<Complex64> {
    // Get the central atom index (typically the first atom)
    let central_atom_idx = 0;

    // Determine the angular momentum size based on the path operator matrix dimensions
    let path_op_size = path_operator.shape()[0];
    let atom_count = structure.atom_count();

    if atom_count == 0 {
        return Err(FmsError::CalculationError(
            "No atoms in structure".to_string(),
        ));
    }

    // Calculate l_size (number of angular momentum components per atom)
    let l_size = path_op_size / atom_count;
    if l_size * atom_count != path_op_size {
        return Err(FmsError::DimensionMismatch(format!(
            "Path operator matrix size {} is not divisible by atom count {}",
            path_op_size, atom_count
        )));
    }

    // Calculate max_l from l_size
    // l_size = (max_l + 1)^2, so max_l = sqrt(l_size) - 1
    let max_l = (l_size as f64).sqrt() as usize - 1;

    // Get the starting index for the central atom in the path operator
    let start_idx = central_atom_idx * l_size;

    // Sum over the diagonal elements of the central atom's block
    // Use Rayon for parallel sum if l_size is large enough
    let mut total_amplitude = if l_size > 32 {
        // Parallel summation for large angular momentum basis
        (0..l_size)
            .into_par_iter()
            .map(|i| path_operator[(start_idx + i, start_idx + i)])
            .reduce(|| Complex64::new(0.0, 0.0), |a, b| a + b)
    } else {
        // Sequential summation for smaller basis
        (0..l_size).fold(Complex64::new(0.0, 0.0), |sum, i| {
            sum + path_operator[(start_idx + i, start_idx + i)]
        })
    };

    // Normalize by the number of angular momentum states (2l+1) for each l
    // This is a physical normalization factor based on quantum mechanics
    let mut normalization_factor = 0.0;
    for l in 0..=max_l {
        // Each l contributes (2l+1) states
        normalization_factor += (2 * l + 1) as f64;
    }

    if normalization_factor > 0.0 {
        total_amplitude /= Complex64::new(normalization_factor, 0.0);
    }

    Ok(total_amplitude)
}

#[cfg(test)]
mod tests {
    #[test]
    fn test_placeholder() {
        assert!(true);
    }
}
