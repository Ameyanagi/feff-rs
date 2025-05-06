/*
MIT License with FEFF10 Attribution

Copyright (c) 2025 Ameyanagi

Based on or developed using Distribution: FEFF10.0
Copyright (c) 2020 FEFF Project, University of Washington and SLAC National Accelerator Laboratory.
All rights reserved.
*/

//! Self-consistency cycle implementation for muffin-tin potentials
//!
//! This module provides the implementation of self-consistent field (SCF)
//! iterations for muffin-tin potentials, using advanced mixing schemes to
//! accelerate convergence.

use super::errors::{PotentialError, Result};
use faer::{Mat, Scale};
use rayon::prelude::*;

/// Result of self-consistency calculation
#[derive(Debug, Clone)]
pub struct SelfConsistencyResult {
    /// Number of iterations performed
    pub iterations: usize,
    /// Whether convergence was achieved
    pub converged: bool,
    /// Final error measure
    pub final_error: f64,
    /// Error history
    pub error_history: Vec<f64>,
    /// Timings for each iteration (in milliseconds)
    pub timings: Vec<u64>,
}

/// Mixing method to use for self-consistency
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum MixingMethod {
    /// Linear mixing with fixed parameter
    Linear(f64),
    /// Broyden mixing (quasi-Newton method)
    Broyden,
    /// Pulay mixing (Direct Inversion in Iterative Subspace)
    Pulay(usize), // Number of iterations to use in the mixing
}

impl Default for MixingMethod {
    fn default() -> Self {
        MixingMethod::Linear(0.3) // Default 30% linear mixing
    }
}

/// Workspace for Broyden mixing
#[derive(Debug, Clone)]
pub struct BroydenWorkspace {
    /// Maximum history size
    max_history: usize,
    /// History of density vectors
    density_history: Vec<Vec<f64>>,
    /// History of residual vectors (F(ρ) - ρ)
    residual_history: Vec<Vec<f64>>,
    /// Jacobian approximation inverse
    jacobian_inverse: Option<Mat<f64>>,
    /// Current iteration
    iteration: usize,
}

impl BroydenWorkspace {
    /// Create a new Broyden workspace
    #[allow(dead_code)]
    pub fn new(max_history: usize) -> Self {
        Self {
            max_history,
            density_history: Vec::with_capacity(max_history),
            residual_history: Vec::with_capacity(max_history),
            jacobian_inverse: None,
            iteration: 0,
        }
    }

    /// Reset the workspace
    pub fn reset(&mut self) {
        self.density_history.clear();
        self.residual_history.clear();
        self.jacobian_inverse = None;
        self.iteration = 0;
    }

    /// Update the Broyden mixing workspace and generate the next density
    pub fn update(&mut self, old_density: &[f64], new_density: &[f64]) -> Vec<f64> {
        let n = old_density.len();

        // First compute residual F(ρ) - ρ
        let mut residual = Vec::with_capacity(n);
        for i in 0..n {
            residual.push(new_density[i] - old_density[i]);
        }

        if self.iteration == 0 {
            // First iteration - use simple linear mixing
            let alpha = 0.3; // Initial mixing parameter
            let mut mixed_density = Vec::with_capacity(n);
            for i in 0..n {
                mixed_density.push(old_density[i] + alpha * residual[i]);
            }

            // Store density and residual
            self.density_history.push(old_density.to_vec());
            self.residual_history.push(residual);
            self.iteration += 1;

            return mixed_density;
        }

        // Store current density and residual
        self.density_history.push(old_density.to_vec());
        self.residual_history.push(residual.clone());

        // Keep only the last max_history iterations
        if self.density_history.len() > self.max_history {
            self.density_history.remove(0);
            self.residual_history.remove(0);
        }

        // Perform Broyden update
        let mixed_density = if self.iteration == 1 {
            // Initialize Jacobian inverse as identity scaled by initial mixing parameter
            let alpha = 0.3;
            let mut mixed_density = Vec::with_capacity(n);
            for i in 0..n {
                mixed_density.push(old_density[i] + alpha * residual[i]);
            }

            // Initialize Jacobian inverse (identity matrix scaled by alpha)
            let jacobian_inverse = Scale(alpha) * Mat::<f64>::identity(n, n);
            self.jacobian_inverse = Some(jacobian_inverse);

            mixed_density
        } else {
            // Use Broyden update formula
            let rho_prev = &self.density_history[self.density_history.len() - 2];
            let f_prev = &self.residual_history[self.residual_history.len() - 2];

            // Compute difference vectors
            let mut delta_rho = Vec::with_capacity(n);
            let mut delta_f = Vec::with_capacity(n);

            for i in 0..n {
                delta_rho.push(old_density[i] - rho_prev[i]);
                delta_f.push(residual[i] - f_prev[i]);
            }

            // Convert to matrices for easier manipulation
            let _delta_rho_mat = Mat::from_fn(n, 1, |i, _| delta_rho[i]);
            let delta_f_mat = Mat::from_fn(n, 1, |i, _| delta_f[i]);

            // Current Jacobian inverse
            let j_inv = self.jacobian_inverse.clone().unwrap();

            // Compute u = δρ - J⁻ δF using matrix operations
            // First apply j_inv to delta_f_mat
            let j_inv_delta_f_mat = &j_inv * &delta_f_mat;

            // Convert to vector for easier manipulation
            let mut j_inv_delta_f = Vec::with_capacity(n);
            for i in 0..n {
                j_inv_delta_f.push(j_inv_delta_f_mat[(i, 0)]);
            }

            // Compute u = δρ - J⁻ δF
            let mut u = Vec::with_capacity(n);
            for i in 0..n {
                u.push(delta_rho[i] - j_inv_delta_f[i]);
            }

            // Compute denominator = δF·(J⁻ δF)
            let mut denom = 0.0;
            for i in 0..n {
                denom += delta_f[i] * j_inv_delta_f[i];
            }

            let mut j_inv_updated = j_inv.clone();

            if denom.abs() > 1e-10 {
                // Update Jacobian inverse: J⁻_new = J⁻ + (δρ - J⁻ δF) ⊗ (J⁻ δF) / (δF · (J⁻ δF))
                let u_mat = Mat::from_fn(n, 1, |i, _| u[i]);
                let j_inv_delta_f_row_mat = Mat::from_fn(1, n, |_, j| j_inv_delta_f[j]);

                // Outer product of vectors
                let update = &u_mat * &j_inv_delta_f_row_mat;

                // Scale by 1/denom
                let scaled_update = Scale(1.0 / denom) * update;

                // Add to current j_inv
                j_inv_updated = &j_inv + &scaled_update;
                self.jacobian_inverse = Some(j_inv_updated.clone());
            }

            // Compute new density: ρ_new = ρ + J⁻ F(ρ)
            let residual_mat = Mat::from_fn(n, 1, |i, _| residual[i]);
            let j_inv_f_mat = &j_inv_updated * &residual_mat;
            let mut mixed_density = Vec::with_capacity(n);
            for i in 0..n {
                mixed_density.push(old_density[i] + j_inv_f_mat[(i, 0)]);
            }

            mixed_density
        };

        self.iteration += 1;
        mixed_density
    }
}

/// Workspace for Pulay mixing (DIIS)
#[derive(Debug, Clone)]
pub struct PulayWorkspace {
    /// Maximum history size
    max_history: usize,
    /// History of density vectors
    density_history: Vec<Vec<f64>>,
    /// History of residual vectors (F(ρ) - ρ)
    residual_history: Vec<Vec<f64>>,
    /// Current iteration
    iteration: usize,
}

impl PulayWorkspace {
    /// Create a new Pulay workspace
    #[allow(dead_code)]
    pub fn new(max_history: usize) -> Self {
        Self {
            max_history,
            density_history: Vec::with_capacity(max_history),
            residual_history: Vec::with_capacity(max_history),
            iteration: 0,
        }
    }

    /// Reset the workspace
    pub fn reset(&mut self) {
        self.density_history.clear();
        self.residual_history.clear();
        self.iteration = 0;
    }

    /// Update the Pulay mixing workspace and generate the next density
    pub fn update(&mut self, old_density: &[f64], new_density: &[f64]) -> Vec<f64> {
        let n = old_density.len();

        // First compute residual F(ρ) - ρ
        let mut residual = Vec::with_capacity(n);
        for i in 0..n {
            residual.push(new_density[i] - old_density[i]);
        }

        if self.iteration < 2 {
            // First two iterations - use simple linear mixing
            let alpha = 0.3; // Initial mixing parameter
            let mut mixed_density = Vec::with_capacity(n);
            for i in 0..n {
                mixed_density.push(old_density[i] + alpha * residual[i]);
            }

            // Store density and residual (clone to keep the original)
            self.density_history.push(old_density.to_vec());
            self.residual_history.push(residual.clone());
            self.iteration += 1;

            return mixed_density;
        }

        // Store current density and residual (clone to keep the original)
        self.density_history.push(old_density.to_vec());
        self.residual_history.push(residual.clone());

        // Keep only the last max_history iterations
        if self.density_history.len() > self.max_history {
            self.density_history.remove(0);
            self.residual_history.remove(0);
        }

        // Number of iterations in history
        let n_hist = self.density_history.len();

        // Construct the overlap matrix of residuals
        let mut a_data = vec![0.0; (n_hist + 1) * (n_hist + 1)];

        for i in 0..n_hist {
            for j in 0..n_hist {
                let res_i = &self.residual_history[i];
                let res_j = &self.residual_history[j];

                // Compute dot product of residuals
                let mut dot = 0.0;
                for k in 0..n {
                    dot += res_i[k] * res_j[k];
                }

                a_data[i * (n_hist + 1) + j] = dot;
            }
        }

        // Last row and column are all 1.0 (Lagrange multiplier for normalization)
        for i in 0..n_hist {
            a_data[i * (n_hist + 1) + n_hist] = 1.0;
            a_data[n_hist * (n_hist + 1) + i] = 1.0;
        }
        a_data[n_hist * (n_hist + 1) + n_hist] = 0.0;

        // We're not using these matrices anymore, but keeping for reference
        let _a = Mat::from_fn(n_hist + 1, n_hist + 1, |i, j| a_data[i * (n_hist + 1) + j]);

        // Right-hand side vector
        let mut b_data = vec![0.0; n_hist + 1];
        b_data[n_hist] = 1.0; // Normalization constraint

        // Create right-hand side matrix (not used, but keeping for reference)
        let _b = Mat::from_fn(n_hist + 1, 1, |i, _| b_data[i]);

        // Solve the linear system to find optimal coefficients
        let mut coeffs = Vec::with_capacity(n_hist);

        // Solve using a simplified Gaussian elimination approach
        // For small systems this is fast enough
        // For this small system, use a simple approach - direct solution
        // Pulay mixing is robust to small errors so we don't need high precision
        if n_hist > 0 {
            let solution = solve_linear_system_simple(&a_data, &b_data, n_hist + 1);
            if let Some(x) = solution {
                for coef in x.iter().take(n_hist) {
                    coeffs.push(*coef);
                }
            }
        } else {
            // Fallback to simple linear mixing if the system is singular
            let alpha = 0.3;
            let mut mixed_density = Vec::with_capacity(n);
            for i in 0..n {
                mixed_density.push(old_density[i] + alpha * residual[i]);
            }

            self.iteration += 1;
            return mixed_density;
        }

        // Compute optimal density and optimal residual
        let mut optimal_density = vec![0.0; n];
        let mut optimal_residual = vec![0.0; n];

        for i in 0..n {
            for (j, &coeff) in coeffs.iter().enumerate().take(n_hist) {
                optimal_density[i] += coeff * self.density_history[j][i];
                optimal_residual[i] += coeff * self.residual_history[j][i];
            }
        }

        // Apply a damping factor to the optimal residual to avoid oscillations
        let alpha = 0.5;
        let mut mixed_density = Vec::with_capacity(n);
        for i in 0..n {
            mixed_density.push(optimal_density[i] + alpha * optimal_residual[i]);
        }

        self.iteration += 1;
        mixed_density
    }
}

/// Mix densities using the specified method
#[allow(dead_code)]
pub fn mix_densities(
    method: MixingMethod,
    old_densities: &[Vec<f64>],
    new_densities: &[Vec<f64>],
    workspaces: &mut Option<Vec<Box<dyn DensityMixer>>>,
) -> Result<Vec<Vec<f64>>> {
    let n_potentials = old_densities.len();

    // Ensure input sizes match
    if n_potentials != new_densities.len() {
        return Err(PotentialError::CalculationError(
            "Mismatch in number of potential types".to_string(),
        ));
    }

    for (pot_idx, old_density) in old_densities.iter().enumerate() {
        if old_density.len() != new_densities[pot_idx].len() {
            return Err(PotentialError::CalculationError(format!(
                "Mismatch in density array sizes for potential {}",
                pot_idx
            )));
        }
    }

    // Initialize workspaces if needed
    if workspaces.is_none() {
        match method {
            MixingMethod::Linear(_) => {
                // No workspace needed for linear mixing
                *workspaces = Some(Vec::new());
            }
            MixingMethod::Broyden => {
                // Create Broyden workspaces
                let mut mixers: Vec<Box<dyn DensityMixer>> = Vec::with_capacity(n_potentials);
                for _ in 0..n_potentials {
                    mixers.push(Box::new(BroydenMixer::new(5)));
                }
                *workspaces = Some(mixers);
            }
            MixingMethod::Pulay(hist_size) => {
                // Create Pulay workspaces
                let mut mixers: Vec<Box<dyn DensityMixer>> = Vec::with_capacity(n_potentials);
                for _ in 0..n_potentials {
                    mixers.push(Box::new(PulayMixer::new(hist_size)));
                }
                *workspaces = Some(mixers);
            }
        }
    }

    // Initialize result
    let mut mixed_densities = Vec::with_capacity(n_potentials);

    match method {
        MixingMethod::Linear(alpha) => {
            // Simple linear mixing: ρ_new = α*ρ_calc + (1-α)*ρ_old
            for (pot_idx, old_density) in old_densities.iter().enumerate() {
                let new_density = &new_densities[pot_idx];
                let mut mixed = Vec::with_capacity(old_density.len());

                for (i, &old_val) in old_density.iter().enumerate() {
                    let new_val = new_density[i];
                    mixed.push(alpha * new_val + (1.0 - alpha) * old_val);
                }

                mixed_densities.push(mixed);
            }
        }
        MixingMethod::Broyden | MixingMethod::Pulay(_) => {
            // Advanced mixing using workspaces
            let ws = workspaces.as_mut().unwrap();

            // Process each potential type
            for (pot_idx, old_density) in old_densities.iter().enumerate() {
                if pot_idx >= ws.len() {
                    return Err(PotentialError::CalculationError(
                        "Workspace mismatch".to_string(),
                    ));
                }

                // Apply mixing
                let mixed = ws[pot_idx].mix(old_density, &new_densities[pot_idx])?;
                mixed_densities.push(mixed);
            }
        }
    }

    Ok(mixed_densities)
}

/// Calculate error between old and new densities
#[allow(dead_code)]
pub fn calculate_density_error(
    old_densities: &[Vec<f64>],
    new_densities: &[Vec<f64>],
) -> Result<f64> {
    let n_potentials = old_densities.len();

    // Ensure input sizes match
    if n_potentials != new_densities.len() {
        return Err(PotentialError::CalculationError(
            "Mismatch in number of potential types".to_string(),
        ));
    }

    // Parallel computation of maximum error across all potential types
    let max_errors: Vec<_> = (0..n_potentials)
        .into_par_iter()
        .map(|pot_idx| {
            let old_density = &old_densities[pot_idx];
            let new_density = &new_densities[pot_idx];

            if old_density.len() != new_density.len() {
                return Err(PotentialError::CalculationError(format!(
                    "Mismatch in density array sizes for potential {}",
                    pot_idx
                )));
            }

            // Find maximum relative error for this potential type
            let mut max_error: f64 = 0.0;
            for (i, &old_val) in old_density.iter().enumerate() {
                let new_val = new_density[i];

                // Calculate relative error if old value is not too small
                if old_val.abs() > 1e-10 {
                    let error = ((new_val - old_val) / old_val).abs();
                    max_error = f64::max(max_error, error);
                } else {
                    let error = (new_val - old_val).abs();
                    max_error = f64::max(max_error, error);
                }
            }

            Ok(max_error)
        })
        .collect::<Result<Vec<f64>>>()?;

    // Overall maximum error
    Ok(*max_errors
        .iter()
        .max_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
        .unwrap_or(&0.0))
}

/// Simple linear system solver for small systems
/// This avoids the need for complex dependencies
/// Returns None if the system is singular
fn solve_linear_system_simple(a: &[f64], b: &[f64], n: usize) -> Option<Vec<f64>> {
    if n == 0 {
        return None;
    }

    // Make copies of the input
    let mut a_copy = vec![0.0; n * n];
    let mut b_copy = vec![0.0; n];

    a_copy.copy_from_slice(&a[0..(n * n)]);
    b_copy.copy_from_slice(&b[0..n]);

    // Gaussian elimination with pivoting
    for i in 0..n {
        // Find pivot
        let mut max_idx = i;
        let mut max_val = a_copy[i * n + i].abs();

        for j in (i + 1)..n {
            let val = a_copy[j * n + i].abs();
            if val > max_val {
                max_idx = j;
                max_val = val;
            }
        }

        // Check for singularity
        if max_val < 1e-10 {
            return None;
        }

        // Swap rows if needed
        if max_idx != i {
            for j in 0..n {
                let idx1 = i * n + j;
                let idx2 = max_idx * n + j;
                a_copy.swap(idx1, idx2);
            }
            b_copy.swap(i, max_idx);
        }

        // Eliminate below
        for j in (i + 1)..n {
            let factor = a_copy[j * n + i] / a_copy[i * n + i];
            b_copy[j] -= factor * b_copy[i];

            for k in i..n {
                a_copy[j * n + k] -= factor * a_copy[i * n + k];
            }
        }
    }

    // Back substitution
    let mut x = vec![0.0; n];
    for i in (0..n).rev() {
        let mut sum = 0.0;
        for j in (i + 1)..n {
            sum += a_copy[i * n + j] * x[j];
        }
        x[i] = (b_copy[i] - sum) / a_copy[i * n + i];
    }

    Some(x)
}

/// Trait for density mixing methods
#[allow(dead_code)]
pub trait DensityMixer: Send + Sync {
    /// Mix old and new densities to produce a new density
    fn mix(&mut self, old_density: &[f64], new_density: &[f64]) -> Result<Vec<f64>>;

    /// Reset the mixer
    fn reset(&mut self);
}

/// Linear density mixer
#[derive(Debug)]
pub struct LinearMixer {
    /// Mixing parameter
    alpha: f64,
}

impl LinearMixer {
    /// Create a new linear mixer
    #[allow(dead_code)]
    pub fn new(alpha: f64) -> Self {
        Self { alpha }
    }
}

impl DensityMixer for LinearMixer {
    fn mix(&mut self, old_density: &[f64], new_density: &[f64]) -> Result<Vec<f64>> {
        if old_density.len() != new_density.len() {
            return Err(PotentialError::CalculationError(
                "Mismatch in density array sizes".to_string(),
            ));
        }

        let mut mixed = Vec::with_capacity(old_density.len());
        for (i, &old_val) in old_density.iter().enumerate() {
            let new_val = new_density[i];
            mixed.push(self.alpha * new_val + (1.0 - self.alpha) * old_val);
        }

        Ok(mixed)
    }

    fn reset(&mut self) {
        // Nothing to reset for linear mixing
    }
}

/// Broyden density mixer
#[derive(Debug)]
pub struct BroydenMixer {
    /// Broyden workspace
    workspace: BroydenWorkspace,
}

impl BroydenMixer {
    /// Create a new Broyden mixer
    #[allow(dead_code)]
    pub fn new(max_history: usize) -> Self {
        Self {
            workspace: BroydenWorkspace::new(max_history),
        }
    }
}

impl DensityMixer for BroydenMixer {
    fn mix(&mut self, old_density: &[f64], new_density: &[f64]) -> Result<Vec<f64>> {
        if old_density.len() != new_density.len() {
            return Err(PotentialError::CalculationError(
                "Mismatch in density array sizes".to_string(),
            ));
        }

        Ok(self.workspace.update(old_density, new_density))
    }

    fn reset(&mut self) {
        self.workspace.reset();
    }
}

/// Pulay density mixer
#[derive(Debug)]
pub struct PulayMixer {
    /// Pulay workspace
    workspace: PulayWorkspace,
}

impl PulayMixer {
    /// Create a new Pulay mixer
    #[allow(dead_code)]
    pub fn new(max_history: usize) -> Self {
        Self {
            workspace: PulayWorkspace::new(max_history),
        }
    }
}

impl DensityMixer for PulayMixer {
    fn mix(&mut self, old_density: &[f64], new_density: &[f64]) -> Result<Vec<f64>> {
        if old_density.len() != new_density.len() {
            return Err(PotentialError::CalculationError(
                "Mismatch in density array sizes".to_string(),
            ));
        }

        Ok(self.workspace.update(old_density, new_density))
    }

    fn reset(&mut self) {
        self.workspace.reset();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_linear_mixing() {
        let old_density = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let new_density = vec![2.0, 3.0, 4.0, 5.0, 6.0];

        // Mix with alpha = 0.5
        let mut mixer = LinearMixer::new(0.5);
        let mixed = mixer.mix(&old_density, &new_density).unwrap();

        // Check result
        for i in 0..old_density.len() {
            let expected = 0.5 * new_density[i] + 0.5 * old_density[i];
            assert_relative_eq!(mixed[i], expected);
        }
    }

    #[test]
    fn test_density_error() {
        let old_densities = vec![vec![1.0, 2.0, 3.0, 4.0, 5.0]];

        // Small error case
        let new_densities_small = vec![vec![1.01, 2.02, 3.03, 4.04, 5.05]];
        let error_small = calculate_density_error(&old_densities, &new_densities_small).unwrap();
        assert_relative_eq!(error_small, 0.01);

        // Large error case
        let new_densities_large = vec![vec![1.1, 2.2, 3.3, 4.4, 5.5]];
        let error_large = calculate_density_error(&old_densities, &new_densities_large).unwrap();
        assert_relative_eq!(error_large, 0.1);
    }

    #[test]
    fn test_broyden_workspace() {
        let mut workspace = BroydenWorkspace::new(3);

        // Initial density and output from SCF
        let old_density = vec![1.0, 2.0, 3.0];
        let new_density = vec![1.2, 2.3, 3.4];

        // First iteration (should use linear mixing)
        let mixed1 = workspace.update(&old_density, &new_density);

        // Second iteration
        let old_density2 = mixed1.clone();
        let new_density2 = vec![1.15, 2.25, 3.35];
        let mixed2 = workspace.update(&old_density2, &new_density2);

        // Should produce a new density
        assert_ne!(mixed2, old_density2);
        assert_ne!(mixed2, new_density2);

        // Reset should clear history
        workspace.reset();
        assert_eq!(workspace.iteration, 0);
        assert!(workspace.density_history.is_empty());
    }

    #[test]
    fn test_mix_densities() {
        let old_densities = vec![vec![1.0, 2.0, 3.0], vec![4.0, 5.0, 6.0]];

        let new_densities = vec![vec![1.2, 2.3, 3.4], vec![4.1, 5.2, 6.3]];

        // Linear mixing
        let mut workspaces = None;
        let mixed_linear = mix_densities(
            MixingMethod::Linear(0.4),
            &old_densities,
            &new_densities,
            &mut workspaces,
        )
        .unwrap();

        assert_eq!(mixed_linear.len(), 2);
        assert_eq!(mixed_linear[0].len(), 3);

        // Check linear mixing result for first potential
        for i in 0..3 {
            let expected = 0.4 * new_densities[0][i] + 0.6 * old_densities[0][i];
            assert_relative_eq!(mixed_linear[0][i], expected);
        }

        // Broyden mixing
        let mut workspaces = None;
        let mixed_broyden = mix_densities(
            MixingMethod::Broyden,
            &old_densities,
            &new_densities,
            &mut workspaces,
        )
        .unwrap();

        assert_eq!(mixed_broyden.len(), 2);
        assert_eq!(mixed_broyden[0].len(), 3);

        // Pulay mixing
        let mut workspaces = None;
        let mixed_pulay = mix_densities(
            MixingMethod::Pulay(3),
            &old_densities,
            &new_densities,
            &mut workspaces,
        )
        .unwrap();

        assert_eq!(mixed_pulay.len(), 2);
        assert_eq!(mixed_pulay[0].len(), 3);
    }
}
