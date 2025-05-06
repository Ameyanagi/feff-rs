/*
MIT License with FEFF10 Attribution

Copyright (c) 2025 Ameyanagi

Based on or developed using Distribution: FEFF10.0
Copyright (c) 2020 FEFF Project, University of Washington and SLAC National Accelerator Laboratory.
All rights reserved.
*/

//! Wavefunctions for muffin-tin potentials
//!
//! This module provides data structures and functions for handling wavefunctions
//! in muffin-tin potentials, which are essential for accurate electron density
//! and energy level calculations in FEFF.

use super::errors::{PotentialError, Result};
use crate::utils::constants::HARTREE_TO_EV;

/// Represents a radial wavefunction for a specific orbital
#[derive(Debug, Clone)]
pub struct MuffinTinWavefunction {
    /// Principal quantum number
    pub n: i32,
    /// Angular momentum quantum number
    pub l: i32,
    /// Occupation number
    pub occupation: f64,
    /// Energy eigenvalue in eV
    pub energy: f64,
    /// Radial wavefunction R(r)
    pub wavefunction: Vec<f64>,
}

impl MuffinTinWavefunction {
    /// Create a new wavefunction
    pub fn new(n: i32, l: i32, energy_hartree: f64, wavefunction: Vec<f64>) -> Self {
        Self {
            n,
            l,
            occupation: 0.0,                        // Will be set later
            energy: energy_hartree * HARTREE_TO_EV, // Convert to eV
            wavefunction,
        }
    }

    /// Set the occupation number
    pub fn with_occupation(mut self, occupation: f64) -> Self {
        self.occupation = occupation;
        self
    }

    /// Get the principal quantum number
    pub fn n(&self) -> i32 {
        self.n
    }

    /// Get the angular momentum quantum number
    pub fn l(&self) -> i32 {
        self.l
    }

    /// Get the energy in eV
    pub fn energy(&self) -> f64 {
        self.energy
    }

    /// Get the energy in Hartree
    pub fn energy_hartree(&self) -> f64 {
        self.energy / HARTREE_TO_EV
    }

    /// Get the occupation number
    pub fn occupation(&self) -> f64 {
        self.occupation
    }

    /// Get a reference to the wavefunction
    pub fn wavefunction(&self) -> &[f64] {
        &self.wavefunction
    }

    /// Calculate the probability density |ψ|²
    pub fn probability_density(&self) -> Result<Vec<f64>> {
        let mut density = Vec::with_capacity(self.wavefunction.len());
        for &psi in &self.wavefunction {
            density.push(psi * psi);
        }
        Ok(density)
    }

    /// Calculate the electron density contribution from this orbital
    /// at a given radius r using n(r) = |R(r)|²/(4πr²)
    pub fn electron_density_at(&self, r: f64, psi_value: f64) -> f64 {
        if r < 1e-10 {
            return 0.0; // Avoid division by zero
        }

        // Degeneracy factor: 2 for spin, (2l+1) for m
        let degeneracy = 2.0 * (2.0 * self.l as f64 + 1.0);

        // Electron density = occupation × degeneracy × |R(r)|²/(4πr²)
        self.occupation * degeneracy * psi_value * psi_value / (4.0 * std::f64::consts::PI * r * r)
    }

    /// Normalize the wavefunction to ∫|R(r)|²dr = 1
    pub fn normalize(&mut self, grid: &[f64]) -> Result<()> {
        if grid.len() != self.wavefunction.len() {
            return Err(PotentialError::CalculationError(
                "Grid and wavefunction lengths do not match".to_string(),
            ));
        }

        // Calculate normalization using trapezoidal rule
        let mut norm = 0.0;
        for i in 1..grid.len() {
            let r1 = grid[i - 1];
            let r2 = grid[i];
            let dr = r2 - r1;

            let psi1 = self.wavefunction[i - 1];
            let psi2 = self.wavefunction[i];

            // For radial wavefunction R(r), the normalization is ∫|R(r)|²dr = 1
            // Integrate using trapezoidal rule
            norm += 0.5 * dr * (psi1 * psi1 + psi2 * psi2);
        }

        if norm <= 0.0 {
            return Err(PotentialError::CalculationError(
                "Wavefunction has zero norm".to_string(),
            ));
        }

        // Apply normalization
        let scale = 1.0 / norm.sqrt();
        for psi in &mut self.wavefunction {
            *psi *= scale;
        }

        Ok(())
    }
}

/// Collection of wavefunctions for a specific atom type
#[derive(Debug, Clone)]
pub struct AtomicWavefunctions {
    /// Atomic number
    pub atomic_number: i32,
    /// Wavefunctions for each orbital
    pub wavefunctions: Vec<MuffinTinWavefunction>,
}

impl AtomicWavefunctions {
    /// Create a new empty collection
    pub fn new(atomic_number: i32) -> Self {
        Self {
            atomic_number,
            wavefunctions: Vec::new(),
        }
    }

    /// Add a wavefunction to the collection
    pub fn add_wavefunction(&mut self, wavefunction: MuffinTinWavefunction) {
        self.wavefunctions.push(wavefunction);
    }

    /// Get a reference to the wavefunction for a given (n,l) pair
    pub fn get_wavefunction(&self, n: i32, l: i32) -> Option<&MuffinTinWavefunction> {
        self.wavefunctions.iter().find(|wf| wf.n == n && wf.l == l)
    }

    /// Get a mutable reference to the wavefunction for a given (n,l) pair
    pub fn get_wavefunction_mut(&mut self, n: i32, l: i32) -> Option<&mut MuffinTinWavefunction> {
        self.wavefunctions
            .iter_mut()
            .find(|wf| wf.n == n && wf.l == l)
    }

    /// Calculate the total electron density on a grid
    pub fn calculate_electron_density(&self, grid: &[f64]) -> Result<Vec<f64>> {
        let mut density = vec![0.0; grid.len()];

        for wf in &self.wavefunctions {
            if wf.wavefunction.len() != grid.len() {
                return Err(PotentialError::CalculationError(
                    "Grid and wavefunction lengths do not match".to_string(),
                ));
            }

            // Degeneracy factor: 2 for spin, (2l+1) for m
            let degeneracy = 2.0 * (2.0 * wf.l as f64 + 1.0);

            for (i, &r) in grid.iter().enumerate() {
                if r < 1e-10 {
                    continue; // Skip very small r to avoid division by zero
                }

                let psi = wf.wavefunction[i];

                // Electron density = occupation × degeneracy × |R(r)|²/(4πr²)
                density[i] +=
                    wf.occupation * degeneracy * psi * psi / (4.0 * std::f64::consts::PI * r * r);
            }
        }

        Ok(density)
    }

    /// Calculate the total number of electrons
    pub fn total_electrons(&self) -> f64 {
        self.wavefunctions
            .iter()
            .map(|wf| wf.occupation * 2.0 * (2.0 * wf.l as f64 + 1.0))
            .sum()
    }

    /// Sort wavefunctions by energy (most negative first)
    pub fn sort_by_energy(&mut self) {
        self.wavefunctions.sort_by(|a, b| {
            b.energy
                .partial_cmp(&a.energy)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_wavefunction_creation() {
        let n = 2;
        let l = 1;
        let energy = -0.5; // Hartree
        let wavefunction = vec![0.0, 0.1, 0.2, 0.3, 0.4, 0.3, 0.2, 0.1, 0.0];

        let wf = MuffinTinWavefunction::new(n, l, energy, wavefunction.clone());

        assert_eq!(wf.n, n);
        assert_eq!(wf.l, l);
        assert_relative_eq!(wf.energy, energy * HARTREE_TO_EV);
        assert_eq!(wf.wavefunction, wavefunction);
        assert_relative_eq!(wf.occupation, 0.0);

        // Test with occupation
        let wf = wf.with_occupation(2.0);
        assert_relative_eq!(wf.occupation, 2.0);
    }

    #[test]
    fn test_wavefunction_normalization() {
        // Create a simple hydrogen-like wavefunction on a grid
        let n = 1;
        let l = 0;
        let energy = -0.5; // Hartree

        // Create a radial grid
        let n_points = 100;
        let r_max = 10.0;
        let dr = r_max / (n_points as f64);
        let mut grid = Vec::with_capacity(n_points);
        let mut wavefunction = Vec::with_capacity(n_points);

        for i in 0..n_points {
            let r = (i as f64 + 0.1) * dr; // Avoid r=0
            grid.push(r);

            // 1s wavefunction for hydrogen: R(r) = 2 * exp(-r)
            let psi = 2.0 * (-r).exp();
            wavefunction.push(psi);
        }

        let mut wf = MuffinTinWavefunction::new(n, l, energy, wavefunction);

        // Normalize the wavefunction
        wf.normalize(&grid).unwrap();

        // Calculate the norm after normalization
        let mut norm = 0.0;
        for i in 1..grid.len() {
            let r1 = grid[i - 1];
            let r2 = grid[i];
            let dr = r2 - r1;

            let psi1 = wf.wavefunction[i - 1];
            let psi2 = wf.wavefunction[i];

            norm += 0.5 * dr * (psi1 * psi1 + psi2 * psi2);
        }

        // The norm should be close to 1
        assert_relative_eq!(norm, 1.0, epsilon = 1e-3);
    }

    #[test]
    fn test_electron_density_calculation() {
        // Create atomic wavefunctions for hydrogen
        let mut atomic_wf = AtomicWavefunctions::new(1);

        // Create a radial grid
        let n_points = 100;
        let r_max = 10.0;
        let dr = r_max / (n_points as f64);
        let mut grid = Vec::with_capacity(n_points);
        let mut wavefunction = Vec::with_capacity(n_points);

        for i in 0..n_points {
            let r = (i as f64 + 0.1) * dr; // Avoid r=0
            grid.push(r);

            // 1s wavefunction for hydrogen: R(r) = 2 * exp(-r)
            let psi = 2.0 * (-r).exp();
            wavefunction.push(psi);
        }

        let mut wf = MuffinTinWavefunction::new(1, 0, -0.5, wavefunction);
        wf.normalize(&grid).unwrap();
        wf = wf.with_occupation(1.0); // For hydrogen 1s

        atomic_wf.add_wavefunction(wf);

        // Calculate electron density
        let density = atomic_wf.calculate_electron_density(&grid).unwrap();

        // Check that density is non-negative
        for &d in &density {
            assert!(d >= 0.0);
        }

        // For hydrogen, the electron density peak should be near the nucleus
        let max_density_idx = density
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
            .map(|(idx, _)| idx)
            .unwrap();

        // The maximum density should be within the first 10% of the grid
        assert!(max_density_idx < grid.len() / 10);

        // Calculate total electrons by integrating 4πr²ρ(r)
        let mut total_electrons = 0.0;
        for i in 1..grid.len() {
            let r1 = grid[i - 1];
            let r2 = grid[i];
            let dr = r2 - r1;

            let rho1 = density[i - 1];
            let rho2 = density[i];

            // Volume elements
            let vol1 = 4.0 * std::f64::consts::PI * r1 * r1;
            let vol2 = 4.0 * std::f64::consts::PI * r2 * r2;

            // Integrate using trapezoidal rule
            total_electrons += 0.5 * dr * (rho1 * vol1 + rho2 * vol2);
        }

        // For hydrogen 1s, total electrons should be close to 2.0 due to spin degeneracy
        // (2 electrons for the 1s shell, which has l=0, so degeneracy = 2*(2*0+1) = 2)
        assert_relative_eq!(total_electrons, 2.0, epsilon = 0.1);

        // Check that total_electrons() method gives the same result
        assert_relative_eq!(atomic_wf.total_electrons(), 2.0); // 2 for 1s shell (2 * (2*0 + 1))
    }
}
