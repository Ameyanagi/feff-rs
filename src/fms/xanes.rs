/*
MIT License with FEFF10 Attribution

Copyright (c) 2025 Ameyanagi

Based on or developed using Distribution: FEFF10.0
Copyright (c) 2020 FEFF Project, University of Washington and SLAC National Accelerator Laboratory.
All rights reserved.
*/

//! XANES calculator
//!
//! This module implements the XANES calculation from the path operator.

use super::errors::{FmsError, Result};
use crate::atoms::AtomicStructure;
use crate::utils::constants::HARTREE_TO_EV;
use crate::utils::math::lorentzian;
use ndarray::Array2;
use num_complex::Complex64;
use rayon::prelude::*;
use std::f64::consts::PI;

/// Calculator for XANES spectra
///
/// This struct handles the calculation of X-ray Absorption Near-Edge Structure
/// (XANES) spectra from the path operator. It computes the transition matrix
/// elements and applies appropriate broadening and normalization.
#[derive(Debug, Clone)]
pub struct XanesCalculator<'a> {
    /// Reference to the atomic structure
    structure: &'a AtomicStructure,
    /// Core hole lifetime in eV (for broadening)
    core_hole_lifetime: f64,
    /// Polarization vector (default: [1, 0, 0])
    polarization: [f64; 3],
    /// Whether to include quadrupole transitions
    include_quadrupole: bool,
    /// Maximum angular momentum for calculations
    max_l: usize,
    /// Fermi energy in eV
    fermi_energy: f64,
    /// Convolution broadening parameter in eV (additional to core hole)
    conv_broadening: f64,
    /// Energy shift in eV
    energy_shift: f64,
    /// Angular momentum of initial state (n)
    initial_n: usize,
    /// Angular momentum of initial state (l)
    initial_l: usize,
}

impl<'a> XanesCalculator<'a> {
    /// Create a new XANES calculator
    ///
    /// # Arguments
    ///
    /// * `structure` - The atomic structure
    /// * `core_hole_lifetime` - Core hole lifetime in eV
    ///
    /// # Returns
    ///
    /// A new XANES calculator with default settings
    pub fn new(structure: &'a AtomicStructure, core_hole_lifetime: f64) -> Self {
        // Determine edge type from central atom
        let (initial_n, initial_l) = if let Some(central_atom) = structure.central_atom() {
            let z = central_atom.atomic_number();
            if z <= 20 {
                // K-edge for lighter elements
                (1, 0)
            } else if z <= 36 {
                // L₃-edge for medium elements
                (2, 1)
            } else {
                // M₅-edge for heavy elements
                (3, 2)
            }
        } else {
            // Default to K-edge if no central atom
            (1, 0)
        };

        Self {
            structure,
            core_hole_lifetime,
            polarization: [1.0, 0.0, 0.0], // Default: x-polarization
            include_quadrupole: false,
            max_l: 3,             // Default l_max = 3
            fermi_energy: 0.0,    // Will be set later
            conv_broadening: 0.1, // Default broadening (eV)
            energy_shift: 0.0,    // No energy shift by default
            initial_n,
            initial_l,
        }
    }

    /// Set the maximum angular momentum
    ///
    /// # Arguments
    ///
    /// * `max_l` - Maximum angular momentum
    pub fn set_max_l(&mut self, max_l: usize) -> &mut Self {
        self.max_l = max_l;
        self
    }

    /// Set the Fermi energy
    ///
    /// # Arguments
    ///
    /// * `fermi_energy` - Fermi energy in eV
    pub fn set_fermi_energy(&mut self, fermi_energy: f64) -> &mut Self {
        self.fermi_energy = fermi_energy;
        self
    }

    /// Set the convolution broadening
    ///
    /// # Arguments
    ///
    /// * `broadening` - Broadening parameter in eV
    pub fn set_conv_broadening(&mut self, broadening: f64) -> &mut Self {
        self.conv_broadening = broadening;
        self
    }

    /// Set the energy shift
    ///
    /// # Arguments
    ///
    /// * `shift` - Energy shift in eV
    pub fn set_energy_shift(&mut self, shift: f64) -> &mut Self {
        self.energy_shift = shift;
        self
    }

    /// Set the initial state
    ///
    /// # Arguments
    ///
    /// * `n` - Principal quantum number
    /// * `l` - Angular momentum quantum number
    pub fn set_initial_state(&mut self, n: usize, l: usize) -> &mut Self {
        self.initial_n = n;
        self.initial_l = l;
        self
    }

    /// Set the polarization vector
    ///
    /// # Arguments
    ///
    /// * `polarization` - Polarization vector [x, y, z]
    pub fn set_polarization(&mut self, polarization: [f64; 3]) -> &mut Self {
        // Normalize the polarization vector
        let norm =
            (polarization[0].powi(2) + polarization[1].powi(2) + polarization[2].powi(2)).sqrt();

        if norm > 1e-10 {
            self.polarization = [
                polarization[0] / norm,
                polarization[1] / norm,
                polarization[2] / norm,
            ];
        }

        self
    }

    /// Enable or disable quadrupole transitions
    ///
    /// # Arguments
    ///
    /// * `include` - Whether to include quadrupole transitions
    pub fn set_include_quadrupole(&mut self, include: bool) -> &mut Self {
        self.include_quadrupole = include;
        self
    }

    /// Calculate the XANES spectrum at a given energy
    ///
    /// # Arguments
    ///
    /// * `energy` - Energy in eV
    /// * `path_operator` - The calculated path operator
    ///
    /// # Returns
    ///
    /// XANES absorption coefficient μ(E)
    pub fn calculate_xanes(&self, energy: f64, path_operator: &Array2<Complex64>) -> Result<f64> {
        // Get the central atom to determine initial state properties
        let central_atom = match self.structure.central_atom() {
            Some(atom) => atom,
            None => {
                return Err(FmsError::CalculationError(
                    "No central atom specified in structure".to_string(),
                ));
            }
        };

        // Get atomic number for edge energy calculations
        let atomic_number = central_atom.atomic_number() as u32;

        // Determine edge energy based on atomic number and initial state
        let edge_energy = self.calculate_edge_energy(atomic_number)?;

        // Apply energy shift
        let shifted_energy = energy + self.energy_shift;

        // Calculate energy relative to Fermi level
        let relative_energy = shifted_energy - self.fermi_energy;

        // Check if energy is below edge (no absorption in that case)
        if relative_energy < edge_energy {
            return Ok(0.0);
        }

        // The dimensions of the path operator matrix represent the total size of the angular momentum basis
        // for all atoms in the FMS cluster, organized as (atom index, angular momentum)
        let path_op_size = path_operator.shape()[0];

        // Calculate l_size (number of angular momentum components per atom) based on max_l
        let l_size = ((self.max_l + 1) * (self.max_l + 1)) as usize;

        // Calculate number of atoms in the FMS cluster
        let num_atoms = path_op_size / l_size;
        if num_atoms * l_size != path_op_size {
            return Err(FmsError::DimensionMismatch(format!(
                "Path operator matrix size {} is not divisible by angular momentum basis size {}",
                path_op_size, l_size
            )));
        }

        // Find central atom index - assuming the first atom in the structure is the central atom
        // In a more complete implementation, we would have a dedicated is_central flag on Atom
        let central_atom_idx = 0;

        // Index of the first angular momentum component for the central atom
        let central_start_idx = central_atom_idx * l_size;

        // Calculate the dipole and quadrupole matrix elements
        let mut cross_section = Complex64::new(0.0, 0.0);

        // Final state angular momentum values depend on the initial state
        // and selection rules for dipole and quadrupole transitions
        let dipole_l = self.initial_l as i32 + 1; // For dipole transitions, Δl = +1
        let quadrupole_l = self.initial_l as i32 + 2; // For quadrupole transitions, Δl = +2

        // Dipole transitions contribute to all spectra
        cross_section += self.calculate_dipole_contribution(
            path_operator,
            central_start_idx,
            dipole_l as usize,
            l_size,
        )?;

        // Include quadrupole transitions if requested
        if self.include_quadrupole {
            cross_section += self.calculate_quadrupole_contribution(
                path_operator,
                central_start_idx,
                quadrupole_l as usize,
                l_size,
            )?;
        }

        // Convert imaginary part to absorption coefficient
        // The XANES μ(E) is proportional to the imaginary part of the cross section
        let raw_absorption = cross_section.im;

        // Apply appropriate broadening
        // The total broadening combines core hole lifetime and convolution broadening
        let total_broadening = self.core_hole_lifetime + self.conv_broadening;

        // Apply energy-dependent broadening (increases with energy above edge)
        let energy_factor = ((shifted_energy - edge_energy) / 10.0).max(0.0);
        let energy_dependent_broadening = total_broadening * (1.0 + energy_factor);

        // Apply Lorentzian broadening centered at this energy
        let broadened_absorption = raw_absorption
            * lorentzian(
                shifted_energy,
                shifted_energy,
                energy_dependent_broadening / 2.0,
            )
            * PI
            * energy_dependent_broadening;

        // If negative (can happen due to numerical errors), return zero
        Ok(broadened_absorption.max(0.0))
    }

    /// Calculate XANES spectra for multiple energies in parallel
    ///
    /// This method is optimized for computing XANES at many energy points
    /// by using Rayon for parallel processing. It's much more efficient than
    /// calling calculate_xanes() repeatedly.
    ///
    /// # Arguments
    ///
    /// * `energies` - Vector of energies in eV
    /// * `path_operator` - The calculated path operator
    ///
    /// # Returns
    ///
    /// Vector of XANES absorption coefficients μ(E) for each energy
    pub fn calculate_xanes_spectrum(
        &self,
        energies: &[f64],
        path_operator: &Array2<Complex64>,
    ) -> Result<Vec<f64>> {
        // Get the central atom to determine initial state properties
        let central_atom = match self.structure.central_atom() {
            Some(atom) => atom,
            None => {
                return Err(FmsError::CalculationError(
                    "No central atom specified in structure".to_string(),
                ));
            }
        };

        // Get atomic number for edge energy calculations
        let atomic_number = central_atom.atomic_number() as u32;

        // Determine edge energy based on atomic number and initial state
        let edge_energy = self.calculate_edge_energy(atomic_number)?;

        // The dimensions of the path operator matrix represent the total size of the angular momentum basis
        let path_op_size = path_operator.shape()[0];

        // Calculate l_size (number of angular momentum components per atom) based on max_l
        let l_size = ((self.max_l + 1) * (self.max_l + 1)) as usize;

        // Calculate number of atoms in the FMS cluster
        let num_atoms = path_op_size / l_size;
        if num_atoms * l_size != path_op_size {
            return Err(FmsError::DimensionMismatch(format!(
                "Path operator matrix size {} is not divisible by angular momentum basis size {}",
                path_op_size, l_size
            )));
        }

        // Find central atom index - assuming the first atom in the structure is the central atom
        let central_atom_idx = 0;

        // Index of the first angular momentum component for the central atom
        let central_start_idx = central_atom_idx * l_size;

        // Calculate transition matrix elements just once (they don't depend on energy)
        // Final state angular momentum values depend on the initial state
        let dipole_l = self.initial_l as i32 + 1; // For dipole transitions, Δl = +1
        let quadrupole_l = self.initial_l as i32 + 2; // For quadrupole transitions, Δl = +2

        // Calculate the dipole contribution (same for all energies)
        let dipole_contribution = self.calculate_dipole_contribution(
            path_operator,
            central_start_idx,
            dipole_l as usize,
            l_size,
        )?;

        // Calculate quadrupole contribution if requested
        let quadrupole_contribution = if self.include_quadrupole {
            self.calculate_quadrupole_contribution(
                path_operator,
                central_start_idx,
                quadrupole_l as usize,
                l_size,
            )?
        } else {
            Complex64::new(0.0, 0.0)
        };

        // Calculate the total cross section (reused for all energies)
        let cross_section = dipole_contribution + quadrupole_contribution;

        // Convert imaginary part to absorption coefficient
        // The XANES μ(E) is proportional to the imaginary part of the cross section
        let raw_absorption = cross_section.im;

        // For very small number of energy points, sequential processing may be more efficient
        if energies.len() < 8 {
            let mut spectrum = Vec::with_capacity(energies.len());

            for &energy in energies {
                // Apply energy shift
                let shifted_energy = energy + self.energy_shift;

                // Calculate energy relative to Fermi level
                let relative_energy = shifted_energy - self.fermi_energy;

                // No absorption below the edge
                if relative_energy < edge_energy {
                    spectrum.push(0.0);
                    continue;
                }

                // The total broadening combines core hole lifetime and convolution broadening
                let total_broadening = self.core_hole_lifetime + self.conv_broadening;

                // Apply energy-dependent broadening (increases with energy above edge)
                let energy_factor = ((shifted_energy - edge_energy) / 10.0).max(0.0);
                let energy_dependent_broadening = total_broadening * (1.0 + energy_factor);

                // Apply Lorentzian broadening centered at this energy
                let broadened_absorption = raw_absorption
                    * lorentzian(
                        shifted_energy,
                        shifted_energy,
                        energy_dependent_broadening / 2.0,
                    )
                    * PI
                    * energy_dependent_broadening;

                // If negative (can happen due to numerical errors), use zero
                spectrum.push(broadened_absorption.max(0.0));
            }

            return Ok(spectrum);
        }

        // For larger energy grids, use parallel processing with Rayon
        let spectrum: Vec<f64> = energies
            .par_iter()
            .map(|&energy| {
                // Apply energy shift
                let shifted_energy = energy + self.energy_shift;

                // Calculate energy relative to Fermi level
                let relative_energy = shifted_energy - self.fermi_energy;

                // No absorption below the edge
                if relative_energy < edge_energy {
                    return 0.0;
                }

                // The total broadening combines core hole lifetime and convolution broadening
                let total_broadening = self.core_hole_lifetime + self.conv_broadening;

                // Apply energy-dependent broadening (increases with energy above edge)
                let energy_factor = ((shifted_energy - edge_energy) / 10.0).max(0.0);
                let energy_dependent_broadening = total_broadening * (1.0 + energy_factor);

                // Apply Lorentzian broadening centered at this energy
                let broadened_absorption = raw_absorption
                    * lorentzian(
                        shifted_energy,
                        shifted_energy,
                        energy_dependent_broadening / 2.0,
                    )
                    * PI
                    * energy_dependent_broadening;

                // If negative (can happen due to numerical errors), use zero
                broadened_absorption.max(0.0)
            })
            .collect();

        Ok(spectrum)
    }

    /// Calculate edge energy based on atomic number and initial state
    pub fn calculate_edge_energy(&self, atomic_number: u32) -> Result<f64> {
        // This is a simplified calculation - in a real implementation,
        // we would use a more detailed database of edge energies
        match (self.initial_n, self.initial_l) {
            (1, 0) => {
                // K-edge: 1s core hole
                let edge_energy = match atomic_number {
                    1..=10 => 0.03 * (atomic_number as f64).powi(2) + 0.5, // Very light elements
                    11..=36 => 0.035 * (atomic_number as f64).powi(2) + 1.0, // Light elements
                    37..=54 => 0.040 * (atomic_number as f64).powi(2) + 1.5, // Medium elements
                    _ => 0.045 * (atomic_number as f64).powi(2) + 2.0,     // Heavy elements
                };
                Ok(edge_energy)
            }
            (2, 0) => {
                // L₁-edge: 2s core hole
                let edge_energy = 0.015 * (atomic_number as f64).powi(2) + 0.3;
                Ok(edge_energy)
            }
            (2, 1) => {
                // L₂,₃-edge: 2p core hole
                let edge_energy = 0.010 * (atomic_number as f64).powi(2) + 0.2;
                Ok(edge_energy)
            }
            (3, 0) => {
                // M₁-edge: 3s core hole
                let edge_energy = 0.005 * (atomic_number as f64).powi(2) + 0.1;
                Ok(edge_energy)
            }
            (3, 1) => {
                // M₂,₃-edge: 3p core hole
                let edge_energy = 0.003 * (atomic_number as f64).powi(2) + 0.05;
                Ok(edge_energy)
            }
            (3, 2) => {
                // M₄,₅-edge: 3d core hole
                let edge_energy = 0.002 * (atomic_number as f64).powi(2) + 0.02;
                Ok(edge_energy)
            }
            _ => Err(FmsError::InvalidParameter(format!(
                "Unsupported initial state: n={}, l={}",
                self.initial_n, self.initial_l
            ))),
        }
    }

    /// Calculate dipole contribution to XANES cross section
    fn calculate_dipole_contribution(
        &self,
        path_operator: &Array2<Complex64>,
        central_start_idx: usize,
        final_l: usize,
        _l_size: usize,
    ) -> Result<Complex64> {
        let mut dipole_sum = Complex64::new(0.0, 0.0);

        // Calculate index range for final l states
        let l_offset = final_l * final_l;
        let l_count = 2 * final_l + 1; // Number of m values for the given l

        // Sum over all m values for the given l
        for m_idx in 0..l_count {
            // Calculate full index in the path operator matrix
            let matrix_idx = central_start_idx + l_offset + m_idx;

            // Get dipole matrix element for this transition (angular part)
            let dipole_element = self
                .calculate_dipole_matrix_element(final_l as i32, (m_idx as i32) - final_l as i32);

            // Get radial matrix element (depends on energy)
            // Use energy relative to edge energy for the photoelectron energy
            let radial_factor = self.calculate_radial_matrix_element(
                0.0, // Energy will be used in future implementations
                final_l,
            );

            // Combine angular and radial parts
            let matrix_element = dipole_element * Complex64::new(radial_factor, 0.0);

            // Multiply by corresponding path operator element
            // For absorption, we need the diagonal element of the path operator
            dipole_sum += matrix_element * path_operator[(matrix_idx, matrix_idx)];
        }

        Ok(dipole_sum)
    }

    /// Calculate quadrupole contribution to XANES cross section
    fn calculate_quadrupole_contribution(
        &self,
        path_operator: &Array2<Complex64>,
        central_start_idx: usize,
        final_l: usize,
        _l_size: usize,
    ) -> Result<Complex64> {
        let mut quadrupole_sum = Complex64::new(0.0, 0.0);

        // Similar to dipole, but with quadrupole matrix elements
        let l_offset = final_l * final_l;
        let l_count = 2 * final_l + 1; // Number of m values for the given l

        for m_idx in 0..l_count {
            let matrix_idx = central_start_idx + l_offset + m_idx;

            let quadrupole_element = self.calculate_quadrupole_matrix_element(
                final_l as i32,
                (m_idx as i32) - final_l as i32,
            );

            // Get radial matrix element for quadrupole transition
            let radial_factor = self.calculate_radial_matrix_element(
                0.0, // Energy will be used in future implementations
                final_l,
            );

            // Combine angular and radial parts
            let matrix_element = quadrupole_element * Complex64::new(radial_factor, 0.0);

            // Multiply by corresponding path operator element
            quadrupole_sum += matrix_element * path_operator[(matrix_idx, matrix_idx)];
        }

        // Quadrupole transitions are typically weaker than dipole transitions
        // The relative strength depends on the specific system
        Ok(quadrupole_sum * 0.01) // Scale factor for quadrupole contribution
    }

    /// Calculate dipole transition matrix elements
    ///
    /// # Arguments
    ///
    /// * `l` - Final state angular momentum
    /// * `m` - Final state magnetic quantum number
    ///
    /// # Returns
    ///
    /// Transition matrix element
    pub fn calculate_dipole_matrix_element(&self, l: i32, m: i32) -> Complex64 {
        // Check angular momentum selection rules
        // For dipole transitions: Δl = ±1
        if l != self.initial_l as i32 + 1 && l != self.initial_l as i32 - 1 {
            return Complex64::new(0.0, 0.0);
        }

        // Calculate dipole transition matrix element based on polarization vector
        // and spherical harmonics
        // We use the dot product of the polarization vector with the
        // spherical harmonic vector components

        // For s -> p transition (l=0 -> l=1)
        if self.initial_l == 0 && l == 1 {
            match m {
                -1 => {
                    // m = -1: Y_1^{-1} ~ (x - iy)/sqrt(2)
                    let coef = 1.0 / 2.0_f64.sqrt();
                    Complex64::new(coef * self.polarization[0], -coef * self.polarization[1])
                }
                0 => {
                    // m = 0: Y_1^0 ~ z
                    Complex64::new(self.polarization[2], 0.0)
                }
                1 => {
                    // m = 1: Y_1^1 ~ -(x + iy)/sqrt(2)
                    let coef = -1.0 / 2.0_f64.sqrt();
                    Complex64::new(coef * self.polarization[0], coef * self.polarization[1])
                }
                _ => Complex64::new(0.0, 0.0), // Invalid m value
            }
        }
        // For p -> d transition (l=1 -> l=2)
        else if self.initial_l == 1 && l == 2 {
            // p -> d transition (more complex expressions for Y_2^m)
            // In a real implementation, we would use proper spherical harmonics
            // and Clebsch-Gordan coefficients

            // Simplified implementation based on m value
            match m {
                -2 => {
                    // For Y_2^{-2}
                    let coef = 0.3;
                    Complex64::new(
                        coef * (self.polarization[0] - self.polarization[1]),
                        -coef * (self.polarization[0] + self.polarization[1]),
                    )
                }
                -1 => {
                    // For Y_2^{-1}
                    let coef = 0.3;
                    Complex64::new(
                        coef * self.polarization[0] * self.polarization[2],
                        -coef * self.polarization[1] * self.polarization[2],
                    )
                }
                0 => {
                    // For Y_2^0
                    let coef = 0.3;
                    Complex64::new(coef * (3.0 * self.polarization[2].powi(2) - 1.0), 0.0)
                }
                1 => {
                    // For Y_2^1
                    let coef = -0.3;
                    Complex64::new(
                        coef * self.polarization[0] * self.polarization[2],
                        coef * self.polarization[1] * self.polarization[2],
                    )
                }
                2 => {
                    // For Y_2^2
                    let coef = 0.3;
                    Complex64::new(
                        coef * (self.polarization[0] + self.polarization[1]),
                        coef * (self.polarization[0] - self.polarization[1]),
                    )
                }
                _ => Complex64::new(0.0, 0.0), // Invalid m value
            }
        }
        // Add more transitions as needed
        else {
            // Default to zero for unsupported transitions
            Complex64::new(0.0, 0.0)
        }
    }

    /// Calculate quadrupole transition matrix elements
    ///
    /// # Arguments
    ///
    /// * `l` - Final state angular momentum
    /// * `m` - Final state magnetic quantum number
    ///
    /// # Returns
    ///
    /// Transition matrix element
    pub fn calculate_quadrupole_matrix_element(&self, l: i32, m: i32) -> Complex64 {
        // Check angular momentum selection rules
        // For quadrupole transitions: Δl = ±2
        if l != self.initial_l as i32 + 2 && l != self.initial_l as i32 - 2 {
            return Complex64::new(0.0, 0.0);
        }

        // For s -> d transition (l=0 -> l=2)
        if self.initial_l == 0 && l == 2 {
            // Quadrupole matrix elements depend on the second moments
            // of the polarization vector (products of components)

            match m {
                -2 => {
                    // For Y_2^{-2}
                    let coef = 0.2;
                    Complex64::new(coef * self.polarization[0] * self.polarization[1], 0.0)
                }
                -1 => {
                    // For Y_2^{-1}
                    let coef = 0.2;
                    Complex64::new(coef * self.polarization[1] * self.polarization[2], 0.0)
                }
                0 => {
                    // For Y_2^0
                    let coef = 0.2;
                    Complex64::new(coef * (3.0 * self.polarization[2].powi(2) - 1.0), 0.0)
                }
                1 => {
                    // For Y_2^1
                    let coef = 0.2;
                    Complex64::new(coef * self.polarization[0] * self.polarization[2], 0.0)
                }
                2 => {
                    // For Y_2^2
                    let coef = 0.2;
                    Complex64::new(
                        coef * (self.polarization[0].powi(2) - self.polarization[1].powi(2)),
                        0.0,
                    )
                }
                _ => Complex64::new(0.0, 0.0), // Invalid m value
            }
        }
        // Add more transitions as needed
        else {
            // Default to zero for unsupported transitions
            Complex64::new(0.0, 0.0)
        }
    }

    /// Calculate the appropriate radial transition matrix element
    ///
    /// # Arguments
    ///
    /// * `energy` - Photoelectron energy in eV
    /// * `final_l` - Angular momentum of final state
    ///
    /// # Returns
    ///
    /// Get the core hole lifetime
    pub fn get_core_hole_lifetime(&self) -> f64 {
        self.core_hole_lifetime
    }

    /// Get the polarization vector
    pub fn get_polarization(&self) -> [f64; 3] {
        self.polarization
    }

    /// Check if quadrupole transitions are included
    pub fn get_quadrupole_included(&self) -> bool {
        self.include_quadrupole
    }

    /// Calculate the appropriate radial transition matrix element
    ///
    /// # Arguments
    ///
    /// * `energy` - Photoelectron energy in eV
    /// * `final_l` - Angular momentum of final state
    ///
    /// # Returns
    ///
    /// Radial transition matrix element
    fn calculate_radial_matrix_element(&self, energy: f64, final_l: usize) -> f64 {
        // Calculate the effective nuclear charge based on initial state
        let z_eff: f64 = match (self.initial_n, self.initial_l) {
            (1, 0) => 1.0, // K-edge: 1s
            (2, 0) => 0.7, // L₁-edge: 2s
            (2, 1) => 0.8, // L₂,₃-edge: 2p
            (3, 0) => 0.5, // M₁-edge: 3s
            (3, 1) => 0.6, // M₂,₃-edge: 3p
            (3, 2) => 0.7, // M₄,₅-edge: 3d
            _ => 0.8,      // Default value
        };

        // Calculate the photoelectron wavenumber (k) based on energy
        // E = hbar^2 * k^2 / (2 * m_e) in atomic units
        let atomic_energy = energy / HARTREE_TO_EV; // Convert from eV to Hartree
        let k = (2.0 * atomic_energy).sqrt();

        // Calculate the radial matrix element based on energy, initial and final state
        // This is a semiclassical approximation for demonstration purposes
        // In a real implementation, you would use proper radial wavefunctions

        // For dipole transition (Δl = ±1)
        if (final_l as i32 - self.initial_l as i32).abs() == 1 {
            // Form factor decreases with energy approximately as k^-3
            let form_factor = z_eff.powi(3) / (1.0 + k.powi(3));

            // Transitions to higher l values are generally stronger
            let l_factor = if final_l > self.initial_l { 1.2 } else { 0.8 };

            // Normalize for typical K-edge cross section
            let normalization = 10.0;

            return form_factor * l_factor * normalization;
        }
        // For quadrupole transition (Δl = ±2)
        else if (final_l as i32 - self.initial_l as i32).abs() == 2 {
            // Quadrupole transitions scale as k^-5 and are weaker
            let form_factor = z_eff.powi(5) / (1.0 + k.powi(5));

            // Transitions to higher l values are generally stronger
            let l_factor = if final_l > self.initial_l { 1.5 } else { 0.5 };

            // Quadrupole transitions are typically much weaker than dipole
            let normalization = 0.1;

            return form_factor * l_factor * normalization;
        } else {
            0.0 // For unsupported transitions
        }
    }
}
