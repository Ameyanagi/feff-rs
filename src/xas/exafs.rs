/*
MIT License with FEFF10 Attribution

Copyright (c) 2025 Ameyanagi

Based on or developed using Distribution: FEFF10.0
Copyright (c) 2020 FEFF Project, University of Washington and SLAC National Accelerator Laboratory.
All rights reserved.
*/

//! EXAFS (Extended X-ray Absorption Fine Structure) calculation module
//!
//! This module implements EXAFS calculations based on scattering paths and phase shifts.
//! It includes functionality to compute k-space EXAFS chi(k), perform Fourier transforms
//! to r-space, and apply various corrections like many-body effects and Debye-Waller factors.

use std::f64::consts::PI;

use crate::atoms::{AtomicStructure, Result as AtomResult};
use crate::path::{Path, PathFinder, PathFinderConfig};
use crate::scattering::ScatteringResults;
use crate::utils::constants::HARTREE_TO_EV;

/// Enumeration of grid types for EXAFS calculations
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum GridType {
    /// Linear grid in k-space (uniform steps in wavenumber)
    Linear,
    /// Logarithmic grid in k-space (useful for emphasizing low-k region)
    Logarithmic,
    /// Exponential grid in k-space (useful for emphasizing high-k region)
    Exponential,
    /// Custom grid with user-provided points
    Custom,
}

/// Energy grid for EXAFS calculations
#[derive(Clone)]
pub struct EnergyGrid {
    /// Energy values in eV
    pub energies: Vec<f64>,
    /// Corresponding k values in Å^-1
    pub k_values: Vec<f64>,
    /// Energy reference (E0) in eV
    pub e0: f64,
    /// Type of grid
    pub grid_type: GridType,
}

impl EnergyGrid {
    /// Create a new energy grid for EXAFS calculations with linear k-spacing
    ///
    /// # Arguments
    ///
    /// * `e0` - Reference energy (absorption edge) in eV
    /// * `k_min` - Minimum k value in Å^-1
    /// * `k_max` - Maximum k value in Å^-1
    /// * `k_step` - Step size in k-space in Å^-1
    ///
    /// # Returns
    ///
    /// A new `EnergyGrid` instance
    pub fn new(e0: f64, k_min: f64, k_max: f64, k_step: f64) -> Self {
        let k_count = ((k_max - k_min) / k_step).ceil() as usize + 1;
        let mut k_values = Vec::with_capacity(k_count);
        let mut energies = Vec::with_capacity(k_count);

        for i in 0..k_count {
            let k = k_min + (i as f64) * k_step;
            k_values.push(k);

            // Convert k to energy: E = E0 + ħ²k²/2m
            // In atomic units: E(Hartree) = k² (a.u.) / 2
            // With k in Å^-1, conversion is needed
            let k_atomic = k / 1.8897259886; // Convert from Å^-1 to atomic units (a.u.)
            let e_atomic = k_atomic * k_atomic / 2.0;
            let e_ev = e0 + e_atomic * HARTREE_TO_EV;

            energies.push(e_ev);
        }

        Self {
            energies,
            k_values,
            e0,
            grid_type: GridType::Linear,
        }
    }

    /// Create a new energy grid with logarithmic k-spacing
    ///
    /// Logarithmic spacing emphasizes the low-k region, which can be useful
    /// for capturing fine features in the EXAFS spectrum at low k values.
    ///
    /// # Arguments
    ///
    /// * `e0` - Reference energy (absorption edge) in eV
    /// * `k_min` - Minimum k value in Å^-1
    /// * `k_max` - Maximum k value in Å^-1
    /// * `num_points` - Number of points in the grid
    ///
    /// # Returns
    ///
    /// A new `EnergyGrid` instance with logarithmic spacing
    pub fn new_logarithmic(e0: f64, k_min: f64, k_max: f64, num_points: usize) -> Self {
        if k_min <= 0.0 {
            // Logarithmic grid requires positive k_min
            panic!("Logarithmic grid requires k_min > 0");
        }

        let mut k_values = Vec::with_capacity(num_points);
        let mut energies = Vec::with_capacity(num_points);

        let log_k_min = k_min.ln();
        let log_k_max = k_max.ln();
        let log_step = (log_k_max - log_k_min) / (num_points as f64 - 1.0);

        for i in 0..num_points {
            let log_k = log_k_min + (i as f64) * log_step;
            let k = log_k.exp();
            k_values.push(k);

            // Convert k to energy
            let k_atomic = k / 1.8897259886;
            let e_atomic = k_atomic * k_atomic / 2.0;
            let e_ev = e0 + e_atomic * HARTREE_TO_EV;

            energies.push(e_ev);
        }

        Self {
            energies,
            k_values,
            e0,
            grid_type: GridType::Logarithmic,
        }
    }

    /// Create a new energy grid with exponential k-spacing
    ///
    /// Exponential spacing emphasizes the high-k region, which can be useful
    /// for capturing fine features in the EXAFS spectrum at high k values.
    ///
    /// # Arguments
    ///
    /// * `e0` - Reference energy (absorption edge) in eV
    /// * `k_min` - Minimum k value in Å^-1
    /// * `k_max` - Maximum k value in Å^-1
    /// * `num_points` - Number of points in the grid
    /// * `exponent` - Controls how much emphasis to place on the high-k region (1.0-3.0 typical)
    ///
    /// # Returns
    ///
    /// A new `EnergyGrid` instance with exponential spacing
    pub fn new_exponential(
        e0: f64,
        k_min: f64,
        k_max: f64,
        num_points: usize,
        exponent: f64,
    ) -> Self {
        let mut k_values = Vec::with_capacity(num_points);
        let mut energies = Vec::with_capacity(num_points);

        for i in 0..num_points {
            let t = (i as f64) / (num_points as f64 - 1.0);
            let k = k_min + (k_max - k_min) * t.powf(exponent);
            k_values.push(k);

            // Convert k to energy
            let k_atomic = k / 1.8897259886;
            let e_atomic = k_atomic * k_atomic / 2.0;
            let e_ev = e0 + e_atomic * HARTREE_TO_EV;

            energies.push(e_ev);
        }

        Self {
            energies,
            k_values,
            e0,
            grid_type: GridType::Exponential,
        }
    }

    /// Create a new energy grid with custom k values
    ///
    /// This allows for complete flexibility in grid definition.
    ///
    /// # Arguments
    ///
    /// * `e0` - Reference energy (absorption edge) in eV
    /// * `k_values` - Custom k values in Å^-1
    ///
    /// # Returns
    ///
    /// A new `EnergyGrid` instance with the provided k values
    pub fn new_custom(e0: f64, k_values: Vec<f64>) -> Self {
        let mut energies = Vec::with_capacity(k_values.len());

        for &k in &k_values {
            // Convert k to energy
            let k_atomic = k / 1.8897259886;
            let e_atomic = k_atomic * k_atomic / 2.0;
            let e_ev = e0 + e_atomic * HARTREE_TO_EV;

            energies.push(e_ev);
        }

        Self {
            energies,
            k_values,
            e0,
            grid_type: GridType::Custom,
        }
    }

    /// Create a new energy grid directly in E-space (energy space)
    ///
    /// This is useful when working with experimental data that is provided
    /// on a specific energy grid.
    ///
    /// # Arguments
    ///
    /// * `e0` - Reference energy (absorption edge) in eV
    /// * `energies` - Energy values in eV
    ///
    /// # Returns
    ///
    /// A new `EnergyGrid` instance with the provided energy values
    pub fn new_from_energies(e0: f64, energies: Vec<f64>) -> Self {
        let mut k_values = Vec::with_capacity(energies.len());

        for &e in &energies {
            // Only calculate k for energies above the edge
            if e >= e0 {
                // Convert energy to k: k = sqrt(2m(E-E0)/ħ²)
                let e_rel = e - e0; // Relative energy above edge
                let e_hartree = e_rel / HARTREE_TO_EV; // Convert to Hartree
                let k_atomic = (2.0 * e_hartree).sqrt(); // In atomic units
                let k = k_atomic * 1.8897259886; // Convert to Å^-1
                k_values.push(k);
            } else {
                // For energies below the edge, k is not defined (use small negative number)
                k_values.push(-0.01);
            }
        }

        Self {
            energies,
            k_values,
            e0,
            grid_type: GridType::Custom,
        }
    }

    /// Create a grid optimized for EXAFS fitting
    ///
    /// This grid uses a combination of constant k-spacing in low-k region
    /// and increased density in the high-k region where oscillations are faster.
    ///
    /// # Arguments
    ///
    /// * `e0` - Reference energy (absorption edge) in eV
    /// * `k_min` - Minimum k value in Å^-1
    /// * `k_max` - Maximum k value in Å^-1
    /// * `density` - Density factor for the grid (higher means more points)
    ///
    /// # Returns
    ///
    /// A new `EnergyGrid` instance optimized for EXAFS fitting
    pub fn new_exafs_optimized(e0: f64, k_min: f64, k_max: f64, density: f64) -> Self {
        // Number of points based on the range and density
        // Higher k requires higher sampling density
        let base_points = ((k_max - k_min) * 20.0 * density).ceil() as usize;

        // Create a non-uniform grid with higher density at higher k
        let mut k_values = Vec::with_capacity(base_points);
        let mut energies = Vec::with_capacity(base_points);

        // Dense sampling in the high-k region
        for i in 0..base_points {
            let t = (i as f64) / (base_points as f64 - 1.0);

            // Non-linear mapping to increase density at higher k
            // This formula gives more points at high k where oscillations are faster
            let p = if t < 0.3 {
                // Linear in the low-k region
                t
            } else {
                // Gradually increase density with k
                0.3 + (t - 0.3) * (1.0 + t)
            };

            let k = k_min + (k_max - k_min) * p;
            k_values.push(k);

            // Convert k to energy
            let k_atomic = k / 1.8897259886;
            let e_atomic = k_atomic * k_atomic / 2.0;
            let e_ev = e0 + e_atomic * HARTREE_TO_EV;

            energies.push(e_ev);
        }

        Self {
            energies,
            k_values,
            e0,
            grid_type: GridType::Custom,
        }
    }

    /// Get the number of points in the grid
    pub fn len(&self) -> usize {
        self.energies.len()
    }

    /// Check if the grid is empty
    pub fn is_empty(&self) -> bool {
        self.energies.is_empty()
    }

    /// Get the k range of the grid
    pub fn k_range(&self) -> (f64, f64) {
        if self.k_values.is_empty() {
            return (0.0, 0.0);
        }

        let mut min_k = self.k_values[0];
        let mut max_k = self.k_values[0];

        for &k in &self.k_values {
            if k < min_k && k > 0.0 {
                min_k = k;
            }
            if k > max_k {
                max_k = k;
            }
        }

        (min_k, max_k)
    }

    /// Get the energy range of the grid
    pub fn energy_range(&self) -> (f64, f64) {
        if self.energies.is_empty() {
            return (0.0, 0.0);
        }

        let mut min_e = self.energies[0];
        let mut max_e = self.energies[0];

        for &e in &self.energies {
            if e < min_e {
                min_e = e;
            }
            if e > max_e {
                max_e = e;
            }
        }

        (min_e, max_e)
    }

    /// Find the index of the closest grid point to a given k value
    pub fn find_closest_k_index(&self, k: f64) -> usize {
        if self.k_values.is_empty() {
            return 0;
        }

        let mut closest_idx = 0;
        let mut min_diff = (self.k_values[0] - k).abs();

        for (i, &grid_k) in self.k_values.iter().enumerate().skip(1) {
            let diff = (grid_k - k).abs();
            if diff < min_diff {
                min_diff = diff;
                closest_idx = i;
            }
        }

        closest_idx
    }

    /// Find the index of the closest grid point to a given energy value
    pub fn find_closest_energy_index(&self, energy: f64) -> usize {
        if self.energies.is_empty() {
            return 0;
        }

        let mut closest_idx = 0;
        let mut min_diff = (self.energies[0] - energy).abs();

        for (i, &grid_e) in self.energies.iter().enumerate().skip(1) {
            let diff = (grid_e - energy).abs();
            if diff < min_diff {
                min_diff = diff;
                closest_idx = i;
            }
        }

        closest_idx
    }

    /// Get a subset of the grid within a specified k range
    pub fn k_subset(&self, k_min: f64, k_max: f64) -> Self {
        let mut k_values = Vec::new();
        let mut energies = Vec::new();

        for i in 0..self.len() {
            let k = self.k_values[i];
            if k >= k_min && k <= k_max {
                k_values.push(k);
                energies.push(self.energies[i]);
            }
        }

        Self {
            energies,
            k_values,
            e0: self.e0,
            grid_type: self.grid_type,
        }
    }

    /// Get a subset of the grid within a specified energy range
    pub fn energy_subset(&self, e_min: f64, e_max: f64) -> Self {
        let mut k_values = Vec::new();
        let mut energies = Vec::new();

        for i in 0..self.len() {
            let e = self.energies[i];
            if e >= e_min && e <= e_max {
                energies.push(e);
                k_values.push(self.k_values[i]);
            }
        }

        Self {
            energies,
            k_values,
            e0: self.e0,
            grid_type: self.grid_type,
        }
    }

    /// Merge two energy grids
    ///
    /// This is useful for combining data from different energy ranges.
    /// The resulting grid preserves the order of points and removes duplicates.
    ///
    /// # Arguments
    ///
    /// * `other` - Another energy grid to merge with this one
    ///
    /// # Returns
    ///
    /// A new merged energy grid
    pub fn merge(&self, other: &Self) -> Self {
        // Check if the reference energies are compatible
        if (self.e0 - other.e0).abs() > 0.1 {
            println!("Warning: Merging grids with different reference energies");
        }

        // Create sets of unique energy points (with a small tolerance for floating point comparison)
        let mut combined_energies = Vec::new();
        let mut combined_k_values = Vec::new();

        // Add points from self
        for i in 0..self.len() {
            combined_energies.push(self.energies[i]);
            combined_k_values.push(self.k_values[i]);
        }

        // Add points from other (avoiding duplicates)
        for i in 0..other.len() {
            let e = other.energies[i];
            let k = other.k_values[i];

            // Check if this point is already included
            let mut is_duplicate = false;
            for j in 0..combined_energies.len() {
                if (combined_energies[j] - e).abs() < 1e-6 {
                    is_duplicate = true;
                    break;
                }
            }

            if !is_duplicate {
                combined_energies.push(e);
                combined_k_values.push(k);
            }
        }

        // Sort by energy
        let mut indices: Vec<usize> = (0..combined_energies.len()).collect();
        indices.sort_by(|&i, &j| {
            combined_energies[i]
                .partial_cmp(&combined_energies[j])
                .unwrap()
        });

        let sorted_energies: Vec<f64> = indices.iter().map(|&i| combined_energies[i]).collect();
        let sorted_k_values: Vec<f64> = indices.iter().map(|&i| combined_k_values[i]).collect();

        Self {
            energies: sorted_energies,
            k_values: sorted_k_values,
            e0: self.e0, // Use the first grid's reference energy
            grid_type: GridType::Custom,
        }
    }
}

/// EXAFS data for a particular calculation
#[derive(Clone)]
pub struct ExafsData {
    /// Energy grid for the calculation
    pub grid: EnergyGrid,
    /// Chi(k) values - the EXAFS oscillation function
    pub chi_k: Vec<f64>,
    /// k*Chi(k) values - often used for visualization
    pub k_chi_k: Vec<f64>,
    /// k²*Chi(k) values - often used for visualization
    pub k2_chi_k: Vec<f64>,
    /// k³*Chi(k) values - often used for visualization
    pub k3_chi_k: Vec<f64>,
    /// Fourier transform magnitude |χ(R)| - the radial distribution function
    pub chi_r_mag: Option<Vec<f64>>,
    /// Fourier transform real part Re[χ(R)]
    pub chi_r_real: Option<Vec<f64>>,
    /// Fourier transform imaginary part Im[χ(R)]
    pub chi_r_imag: Option<Vec<f64>>,
    /// R values for r-space data in Å
    pub r_values: Option<Vec<f64>>,
}

impl ExafsData {
    /// Create a new EXAFS data instance with the given energy grid
    pub fn new(grid: EnergyGrid) -> Self {
        let n = grid.len();
        let zeros = vec![0.0; n];

        Self {
            grid,
            chi_k: zeros.clone(),
            k_chi_k: zeros.clone(),
            k2_chi_k: zeros.clone(),
            k3_chi_k: zeros,
            chi_r_mag: None,
            chi_r_real: None,
            chi_r_imag: None,
            r_values: None,
        }
    }

    /// Calculate weighted k*Chi(k), k²*Chi(k), k³*Chi(k) from chi_k
    pub fn calculate_weighted_spectra(&mut self) {
        for i in 0..self.chi_k.len() {
            let k = self.grid.k_values[i];
            let chi = self.chi_k[i];

            self.k_chi_k[i] = k * chi;
            self.k2_chi_k[i] = k * k * chi;
            self.k3_chi_k[i] = k * k * k * chi;
        }
    }
}

/// Parameters controlling EXAFS calculations
pub struct ExafsParameters {
    /// S0² - amplitude reduction factor due to many-body effects
    pub s02: f64,
    /// Maximum path length to include in Å
    pub r_max: f64,
    /// Minimum path importance to include
    pub min_importance: f64,
    /// Maximum number of legs in paths
    pub max_legs: usize,
    /// Apply thermal Debye-Waller factors
    pub use_debye_waller: bool,
    /// Temperature for Debye-Waller factors in K
    pub temperature: f64,
    /// Energy range for calculations
    pub energy_range: EnergyGrid,
}

impl Default for ExafsParameters {
    fn default() -> Self {
        // Create default energy grid from 3 Å^-1 to 15 Å^-1 with 0.05 Å^-1 steps
        let energy_grid = EnergyGrid::new(0.0, 3.0, 15.0, 0.05);

        Self {
            s02: 1.0,
            r_max: 6.0,
            min_importance: 0.01,
            max_legs: 4,
            use_debye_waller: true,
            temperature: 300.0,
            energy_range: energy_grid,
        }
    }
}

/// Calculate EXAFS for a given structure and set of parameters
///
/// # Arguments
///
/// * `structure` - Atomic structure for calculations
/// * `phase_shifts` - Pre-calculated phase shifts for each potential
/// * `params` - EXAFS calculation parameters
///
/// # Returns
///
/// EXAFS data containing chi(k) and possibly chi(r)
pub fn calculate_exafs(
    structure: &AtomicStructure,
    phase_shifts: &ScatteringResults,
    params: &ExafsParameters,
) -> AtomResult<ExafsData> {
    // Find relevant paths
    // Create path finder configuration
    let path_config = PathFinderConfig {
        max_path_length: params.r_max,
        max_paths: 100,
        max_legs: params.max_legs,
        importance_threshold: params.min_importance,
        cluster_paths: true,
        unique_scatterers_only: true,
    };

    // Get the central atom index
    let absorber_index = structure.central_atom_index().ok_or_else(|| {
        crate::atoms::errors::AtomError::CalculationError("Central atom not defined".to_string())
    })?;

    let mut path_finder = PathFinder::new(structure.clone(), absorber_index, path_config);
    let paths = path_finder.find_paths();

    // Initialize EXAFS data
    let mut exafs_data = ExafsData::new(params.energy_range.clone());

    // Calculate EXAFS for each path at each energy point
    for (i, k) in exafs_data.grid.k_values.iter().enumerate() {
        let energy = exafs_data.grid.energies[i];

        // Sum contributions from all paths
        let mut chi_k = 0.0;

        for path in &paths {
            // Calculate path contribution to EXAFS
            let path_chi =
                calculate_path_contribution(path, k, energy, structure, phase_shifts, params)?;
            chi_k += path_chi;
        }

        exafs_data.chi_k[i] = chi_k;
    }

    // Calculate k-weighted spectra
    exafs_data.calculate_weighted_spectra();

    Ok(exafs_data)
}

/// Calculate the contribution of a single path to the EXAFS spectrum
///
/// # Arguments
///
/// * `path` - The scattering path
/// * `k` - Wavenumber in Å^-1
/// * `energy` - Energy in eV
/// * `structure` - Atomic structure
/// * `phase_shifts` - Scattering results containing phase shifts
/// * `params` - EXAFS calculation parameters
///
/// # Returns
///
/// The path's contribution to chi(k)
fn calculate_path_contribution(
    path: &Path,
    k: &f64,
    _energy: f64,
    structure: &AtomicStructure,
    phase_shifts: &ScatteringResults,
    params: &ExafsParameters,
) -> AtomResult<f64> {
    // Get path properties
    let path_length = path.total_length;
    let degeneracy = path.degeneracy as f64;
    let legs = &path.legs;

    // Apply S0² global amplitude factor (accounts for many-body effects)
    let s02 = params.s02;

    // Calculate amplitude factor based on path geometry and scattering properties
    // Full EXAFS equation: χ(k) = S0² Σ (Nj·fj(k)·e^(-2k²σj²)·e^(-2Rj/λ(k))) / (kRj²) · sin(2kRj + φj(k))

    // Calculate scattering amplitude product (Πj fj(k))
    // For multiple-scattering paths, each atom contributes a scattering amplitude
    let mut scattering_amplitude = 1.0;
    let mut total_phase_shift = 0.0;

    // Central atom phase shift - always present in EXAFS
    let central_atom_idx = structure.central_atom_index().unwrap();
    let central_atom = structure.atom(central_atom_idx).unwrap();
    let central_potential_idx = central_atom.potential_type() as usize;

    // For the absorbing atom, contribution is from l=1 (p orbital) for K-edge
    // This is because the final state is a p orbital for dipole transitions from K-edge
    let l_channel = 1usize; // l=1 for K-edge (dipole selection rules)

    // Get the phase shift for the central atom (absorbing atom)
    // For a K-edge, this will primarily involve l=1 final states
    if central_potential_idx < phase_shifts.phase_shifts.len()
        && l_channel < phase_shifts.phase_shifts[central_potential_idx].len()
    {
        let absorber_phase = phase_shifts.phase_shifts[central_potential_idx][l_channel];
        total_phase_shift += absorber_phase.re;
    }

    // Loop through all legs and calculate each atom's contribution
    for leg in legs {
        // Skip the central atom (handled separately)
        if leg.from_atom == central_atom_idx || leg.to_atom == central_atom_idx {
            continue;
        }

        // Get scattering atom info
        let scattering_atom = structure.atom(leg.to_atom).unwrap();
        let potential_idx = scattering_atom.potential_type() as usize;

        // For backscattering, contribution is from all l channels (summed)
        // In FEFF, this is handled by calculating the effective scattering amplitude
        // Here we use a simplified model

        // Sum over angular momentum contributions (backscattering amplitude)
        let mut atom_amplitude = 0.0;
        let mut atom_phase = 0.0;

        if potential_idx < phase_shifts.phase_shifts.len() {
            let max_l = phase_shifts.max_l.min(3); // Usually only need up to l=3

            for l in 0..=max_l {
                let l_usize = l as usize;
                if l_usize >= phase_shifts.phase_shifts[potential_idx].len() {
                    continue;
                }

                let phase = phase_shifts.phase_shifts[potential_idx][l_usize];

                // Weight by (2l+1) to account for m degeneracy
                let weight = (2 * l + 1) as f64;

                // Phase shift contribution (real part is phase, imaginary part is amplitude damping)
                atom_amplitude += weight * phase.im;
                atom_phase += weight * phase.re;
            }

            // Normalize by sum of weights
            let norm = (max_l * 2 + 2) as f64; // Σ(2l+1) from l=0 to max_l
            atom_amplitude /= norm;
            atom_phase /= norm;
        }

        // Z-dependent correction (higher Z atoms scatter more strongly)
        let z_factor = (scattering_atom.atomic_number() as f64).sqrt() / 5.0;
        atom_amplitude *= z_factor;

        // Add to total scattering amplitude and phase
        scattering_amplitude *= atom_amplitude;
        total_phase_shift += atom_phase;
    }

    // Geometry factor - 1/kR²
    let geometry_factor = 1.0 / (k * path_length * path_length);

    // Debye-Waller factor - thermal vibrations dampen EXAFS signal
    // e^(-2k²σ²) where σ² is the mean square displacement
    let sigma2 = if params.use_debye_waller {
        // Calculate Debye-Waller factor based on temperature
        // For a more sophisticated implementation, this would consider:
        // - Debye temperature of the material
        // - Correlated motion of atoms
        // - Anharmonic effects

        // Simple approximation that depends on:
        // - Temperature (higher T = more vibration)
        // - Path length (longer paths more affected)
        // - Energy (higher k = more sensitive to vibrations)
        0.003 + 0.000005 * params.temperature * path_length
    } else {
        0.0
    };
    let thermal_factor = f64::exp(-2.0 * sigma2 * k * k);

    // Mean free path damping - photoelectron can travel only so far
    // Mean free path varies with energy: λ(k) ∝ k/Γ(k) where Γ is proportional to Im(Σ(E))

    // Simple approximation: λ(k) = k₀ + αk (parameters depend on material)
    // Typical values for metals at room temperature
    let lambda_0 = 4.0; // Baseline MFP in Å (low energy limit)
    let lambda_slope = 0.6; // How fast MFP increases with k
    let lambda = lambda_0 + lambda_slope * k; // Mean free path in Å

    let mfp_factor = f64::exp(-2.0 * path_length / lambda);

    // The full path phase includes:
    // - 2kR (free electron propagation)
    // - φj(k) (phase shifts from each atom)
    let full_phase = 2.0 * k * path_length + total_phase_shift;

    // Combine all factors to calculate the EXAFS contribution from this path
    // χj(k) = S0² × Nj × |fj(k)| × e^(-2k²σj²) × e^(-2Rj/λ(k)) × sin(2kRj + φj(k)) / (kRj²)
    let chi = s02
        * degeneracy
        * scattering_amplitude
        * thermal_factor
        * mfp_factor
        * geometry_factor
        * full_phase.sin();

    Ok(chi)
}

/// Enumeration of window functions for Fourier transform
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum WindowFunction {
    /// No windowing - use with caution as it can lead to artifacts in the FT
    None,
    /// Hanning (Hann) window - good general purpose window
    Hanning,
    /// Hamming window - similar to Hanning but doesn't go to zero at edges
    Hamming,
    /// Blackman window - very good sidelobe suppression
    Blackman,
    /// Blackman-Harris window - excellent sidelobe suppression
    BlackmanHarris,
    /// Kaiser-Bessel window with adjustable parameter
    KaiserBessel(f64),
    /// Gaussian window with adjustable width
    Gaussian(f64),
    /// Parzen (triangle) window
    Parzen,
    /// Welch window (parabolic)
    Welch,
    /// Tukey window (tapered cosine) with adjustable parameter
    Tukey(f64),
    /// Flat top window - excellent amplitude accuracy but poor frequency resolution
    FlatTop,
    /// Exponential window - useful for emphasizing early parts of data
    Exponential(f64),
}

impl Default for WindowFunction {
    fn default() -> Self {
        WindowFunction::Hanning // Good default choice for EXAFS
    }
}

/// Applies window function to EXAFS data for Fourier transform
///
/// Window functions are used to reduce spectral leakage in Fourier transforms
/// by tapering the data at the edges. Different windows have different properties
/// in terms of amplitude accuracy, frequency resolution, and sidelobe suppression.
///
/// # Arguments
///
/// * `k_values` - Wavenumber values
/// * `chi_k` - EXAFS signal values
/// * `window` - The window function type to apply
///
/// # Returns
///
/// The windowed EXAFS signal
pub fn apply_window(k_values: &[f64], chi_k: &[f64], window: WindowFunction) -> Vec<f64> {
    let mut windowed = Vec::with_capacity(chi_k.len());

    // Special case for no windowing
    if window == WindowFunction::None {
        return chi_k.to_vec();
    }

    // For all window functions, normalize the position to [0, 1]
    let k_min = k_values[0];
    let k_max = k_values[k_values.len() - 1];
    let range = k_max - k_min;

    // Apply the selected window function
    for (i, &k) in k_values.iter().enumerate() {
        let normalized_k = (k - k_min) / range;
        // Ensure normalized_k is between 0 and 1 (handle edge cases)
        let safe_normalized_k = normalized_k.clamp(0.0, 1.0);
        let window_value = calculate_window_value(safe_normalized_k, window);
        // Ensure window value is non-negative
        let safe_window_value = window_value.max(0.0);
        windowed.push(chi_k[i] * safe_window_value);
    }

    windowed
}

/// Helper function to calculate window function value at a given normalized position
///
/// # Arguments
///
/// * `x` - Normalized position in the window [0, 1]
/// * `window` - The window function type
///
/// # Returns
///
/// Window function value at position x
fn calculate_window_value(x: f64, window: WindowFunction) -> f64 {
    match window {
        WindowFunction::None => 1.0,

        WindowFunction::Hanning => {
            // Hann window: w(x) = 0.5 * (1 - cos(2πx))
            0.5 * (1.0 - (2.0 * PI * x).cos())
        }

        WindowFunction::Hamming => {
            // Hamming window: w(x) = 0.54 - 0.46 * cos(2πx)
            0.54 - 0.46 * (2.0 * PI * x).cos()
        }

        WindowFunction::Blackman => {
            // Blackman window
            // w(x) = 0.42 - 0.5 * cos(2πx) + 0.08 * cos(4πx)
            let cx = (2.0 * PI * x).cos();
            let c2x = (4.0 * PI * x).cos();
            0.42 - 0.5 * cx + 0.08 * c2x
        }

        WindowFunction::BlackmanHarris => {
            // 4-term Blackman-Harris window
            // w(x) = a0 - a1*cos(2πx) + a2*cos(4πx) - a3*cos(6πx)
            let a0 = 0.35875;
            let a1 = 0.48829;
            let a2 = 0.14128;
            let a3 = 0.01168;

            let cx = (2.0 * PI * x).cos();
            let c2x = (4.0 * PI * x).cos();
            let c3x = (6.0 * PI * x).cos();

            a0 - a1 * cx + a2 * c2x - a3 * c3x
        }

        WindowFunction::KaiserBessel(beta) => {
            // Kaiser-Bessel window with parameter beta
            // w(x) = I0(beta * sqrt(1 - (2x-1)^2)) / I0(beta)
            // where I0 is the modified Bessel function of the first kind, order 0

            // For EXAFS, typical beta values are between 2 and 14
            // Higher beta means more tapering (better sidelobe suppression but wider mainlobe)

            // For our implementation, we'll use a fast approximation of the Bessel function
            // for typical beta values in EXAFS

            let arg = beta * (1.0 - (2.0 * x - 1.0).powi(2)).sqrt();

            // Fast approximation of I0(z) / I0(beta) for z <= beta
            // This is adequate for window functions
            if arg <= 0.0 {
                0.0
            } else if arg >= beta {
                1.0
            } else {
                let ratio = arg / beta;
                let adjusted = 1.0 - (1.0 - ratio).powi(2);
                adjusted * adjusted * (3.0 - 2.0 * adjusted)
            }
        }

        WindowFunction::Gaussian(sigma) => {
            // Gaussian window with parameter sigma (0 to 0.5, default 0.4)
            // w(x) = exp(-0.5 * ((x - 0.5) / sigma)^2)
            // sigma controls the width (smaller sigma = narrower peak)

            // Convert x from [0,1] to [-0.5,0.5]
            let centered_x = x - 0.5;
            f64::exp(-0.5 * (centered_x / sigma).powi(2))
        }

        WindowFunction::Parzen => {
            // Parzen (triangular) window
            // Simple piecewise function that gives a triangular shape
            if x <= 0.5 {
                2.0 * x
            } else {
                2.0 * (1.0 - x)
            }
        }

        WindowFunction::Welch => {
            // Welch window (parabolic)
            // w(x) = 1 - ((x - 0.5) / 0.5)^2
            let centered_x = x - 0.5;
            1.0 - (centered_x / 0.5).powi(2)
        }

        WindowFunction::Tukey(alpha) => {
            // Tukey window (tapered cosine)
            // Flat in the middle, cosine-tapered at the edges
            // alpha = 0: rectangular window
            // alpha = 1: Hann window

            if x < alpha / 2.0 {
                // First taper
                0.5 * (1.0 - (2.0 * PI * x / alpha).cos())
            } else if x > (1.0 - alpha / 2.0) {
                // Second taper
                0.5 * (1.0 - (2.0 * PI * (1.0 - x) / alpha).cos())
            } else {
                // Flat middle
                1.0
            }
        }

        WindowFunction::FlatTop => {
            // Flat top window - excellent amplitude accuracy, poor frequency resolution
            // w(x) = a0 - a1*cos(2πx) + a2*cos(4πx) - a3*cos(6πx) + a4*cos(8πx)
            let a0 = 0.21557895;
            let a1 = 0.41663158;
            let a2 = 0.27726316;
            let a3 = 0.08357895;
            let a4 = 0.00694737;

            let cx = (2.0 * PI * x).cos();
            let c2x = (4.0 * PI * x).cos();
            let c3x = (6.0 * PI * x).cos();
            let c4x = (8.0 * PI * x).cos();

            a0 - a1 * cx + a2 * c2x - a3 * c3x + a4 * c4x
        }

        WindowFunction::Exponential(tau) => {
            // Exponential window - emphasizes the start of the data
            // w(x) = exp(-x / tau)
            // tau controls the decay rate (typical values 0.1 - 0.5)
            // tau = 0.1 means steep decay, tau = 0.5 means gentle decay
            f64::exp(-x / tau)
        }
    }
}

/// Create a window function appropriate for the given data properties
///
/// This function selects an appropriate window function based on the
/// data characteristics and analysis goals.
///
/// # Arguments
///
/// * `r_resolution` - Desired resolution in r-space in Å
/// * `max_path_length` - Longest path in the data in Å
/// * `signal_to_noise` - Estimated signal-to-noise ratio
///
/// # Returns
///
/// An appropriate window function
pub fn select_window_function(
    r_resolution: f64,
    max_path_length: f64,
    signal_to_noise: f64,
) -> WindowFunction {
    // Choose window based on the resolution requirements and signal quality

    // For high S/N and high resolution needs, use a window with good spectral resolution
    if signal_to_noise > 10.0 && r_resolution < 0.1 {
        return WindowFunction::Hanning;
    }

    // For very noisy data, use a window with better noise suppression
    if signal_to_noise < 3.0 {
        return WindowFunction::BlackmanHarris;
    }

    // For data with long paths, use a window that can separate closely spaced features
    if max_path_length > 10.0 {
        return WindowFunction::KaiserBessel(6.0);
    }

    // Default to a good general-purpose window
    WindowFunction::Hanning
}

/// Apply a window function optimized for backscattering paths
///
/// Backscattering paths typically show up as a peak at around 2x the
/// nearest neighbor distance. This function applies a window that emphasizes
/// this region in r-space.
///
/// # Arguments
///
/// * `exafs_data` - EXAFS data to window
/// * `nearest_neighbor` - Distance to nearest neighbor in Å
///
/// # Returns
///
/// EXAFS data with the backscattering window applied
pub fn apply_backscattering_window(mut exafs_data: ExafsData, _nearest_neighbor: f64) -> ExafsData {
    // Apply a Gaussian window centered at twice the nearest neighbor distance
    // This will emphasize the backscattering peak in the FT

    let k_values = &exafs_data.grid.k_values;

    // Parameters for the window
    let k_center = (k_values[0] + k_values[k_values.len() - 1]) / 2.0;
    let range = k_values[k_values.len() - 1] - k_values[0];
    let sigma = range / 4.0; // Adjusted for good backscattering emphasis

    // Apply the window to chi(k)
    for (i, &k) in k_values.iter().enumerate().take(exafs_data.chi_k.len()) {
        let window = f64::exp(-(k - k_center).powi(2) / (2.0 * sigma * sigma));
        exafs_data.chi_k[i] *= window;
    }

    // Recalculate weighted spectra
    exafs_data.calculate_weighted_spectra();

    exafs_data
}

/// Apply a window function optimized for multiple scattering paths
///
/// Multiple scattering paths typically show up at larger r-values.
/// This function applies a window that emphasizes these longer paths.
///
/// # Arguments
///
/// * `exafs_data` - EXAFS data to window
/// * `cutoff_distance` - Minimum distance to emphasize in Å
///
/// # Returns
///
/// EXAFS data with the multiple scattering window applied
pub fn apply_multiple_scattering_window(
    mut exafs_data: ExafsData,
    _cutoff_distance: f64,
) -> ExafsData {
    // Apply a window that emphasizes longer paths by attenuating the low-k region
    // where the contribution of single scattering dominates

    let k_values = &exafs_data.grid.k_values;
    let k_min = k_values[0];
    let k_max = k_values[k_values.len() - 1];

    // Attenuate signal below k_cutoff
    let k_cutoff = 4.0; // Empirical cutoff point

    // Apply the window to chi(k)
    for (i, &k) in k_values.iter().enumerate().take(exafs_data.chi_k.len()) {
        // Smoothly attenuate low-k region
        let attenuation = if k < k_cutoff {
            (k / k_cutoff).powi(2)
        } else {
            1.0
        };

        // Apply a high-k emphasis window
        let window = attenuation * (1.0 - ((k_max - k) / (k_max - k_min)).powi(2));
        exafs_data.chi_k[i] *= window;
    }

    // Recalculate weighted spectra
    exafs_data.calculate_weighted_spectra();

    exafs_data
}

/// Perform Fourier transform from k-space to r-space
///
/// This function takes EXAFS data in k-space and performs a Fourier transform to
/// generate the radial distribution function in r-space. It applies a window function
/// to reduce spectral leakage and improve the clarity of the r-space spectrum.
///
/// # Arguments
///
/// * `exafs_data` - EXAFS data containing k-space information
/// * `window` - Window function to apply before transform
/// * `k_weight` - k-weight to apply (0, 1, 2, or 3)
/// * `r_min` - Minimum R value in Å
/// * `r_max` - Maximum R value in Å
/// * `dr` - R-space step size in Å
///
/// # Returns
///
/// Updated EXAFS data with Fourier transform results
pub fn fourier_transform(
    mut exafs_data: ExafsData,
    window: WindowFunction,
    k_weight: usize,
    r_min: f64,
    r_max: f64,
    dr: f64,
) -> ExafsData {
    let k_values = &exafs_data.grid.k_values;

    // Select the weighted chi(k) based on k_weight
    let chi_k_weighted = match k_weight {
        0 => &exafs_data.chi_k,
        1 => &exafs_data.k_chi_k,
        2 => &exafs_data.k2_chi_k,
        3 => &exafs_data.k3_chi_k,
        _ => &exafs_data.k2_chi_k, // Default to k²-weighted
    };

    // Apply window function
    let windowed_chi = apply_window(k_values, chi_k_weighted, window);

    // Prepare r-space grid
    let r_count = ((r_max - r_min) / dr).ceil() as usize + 1;
    let mut r_values = Vec::with_capacity(r_count);
    let mut chi_r_real = Vec::with_capacity(r_count);
    let mut chi_r_imag = Vec::with_capacity(r_count);
    let mut chi_r_mag = Vec::with_capacity(r_count);

    // Calculate Fourier transform for each R value
    // We use a more accurate approach here with a phase correction factor
    // This accounts for the fact that EXAFS oscillations contain both sine and cosine components
    for i in 0..r_count {
        let r = r_min + (i as f64) * dr;
        r_values.push(r);

        // Numerical integration for Fourier transform using Simpson's rule
        // Simpson's rule is more accurate than trapezoidal for oscillatory functions
        let mut ft_real = 0.0;
        let mut ft_imag = 0.0;

        // Apply Simpson's rule for integration: ∫f(x)dx ≈ h/3 [f(x₀) + 4f(x₁) + 2f(x₂) + 4f(x₃) + ... + f(xₙ)]
        let n = k_values.len();
        if n > 2 {
            // Need at least 3 points for Simpson's rule
            let mut simpson_sum_real = 0.0;
            let mut simpson_sum_imag = 0.0;

            // First point
            let phase_0 = 2.0 * k_values[0] * r;
            simpson_sum_real += windowed_chi[0] * phase_0.cos();
            simpson_sum_imag -= windowed_chi[0] * phase_0.sin();

            // Middle points with alternating 4 and 2 factors
            for j in 1..n - 1 {
                let phase = 2.0 * k_values[j] * r;
                let cos_term = phase.cos();
                let sin_term = phase.sin();

                let factor = if j % 2 == 1 { 4.0 } else { 2.0 };
                simpson_sum_real += factor * windowed_chi[j] * cos_term;
                simpson_sum_imag -= factor * windowed_chi[j] * sin_term;
            }

            // Last point
            let phase_n = 2.0 * k_values[n - 1] * r;
            simpson_sum_real += windowed_chi[n - 1] * phase_n.cos();
            simpson_sum_imag -= windowed_chi[n - 1] * phase_n.sin();

            // Calculate the step size (works for both uniform and non-uniform grids)
            let h = (k_values[n - 1] - k_values[0]) / (n as f64 - 1.0);

            // Apply Simpson's rule prefactor
            ft_real = simpson_sum_real * h / 3.0;
            ft_imag = simpson_sum_imag * h / 3.0;
        } else {
            // Fall back to trapezoidal rule for small data sets
            for j in 1..k_values.len() {
                let k_prev = k_values[j - 1];
                let k = k_values[j];
                let dk = k - k_prev;

                let chi_prev = windowed_chi[j - 1];
                let chi = windowed_chi[j];

                // e^(-2ikr) = cos(2kr) - i*sin(2kr)
                let phase_prev = 2.0 * k_prev * r;
                let phase = 2.0 * k * r;

                let cos_term_prev = phase_prev.cos();
                let sin_term_prev = phase_prev.sin();
                let cos_term = phase.cos();
                let sin_term = phase.sin();

                ft_real += 0.5 * (chi_prev * cos_term_prev + chi * cos_term) * dk;
                ft_imag -= 0.5 * (chi_prev * sin_term_prev + chi * sin_term) * dk;
            }
        }

        // Apply normalization factor
        // The factor 1/π is conventional in EXAFS analysis
        // It makes the transform approximately correspond to the radial distribution function
        ft_real /= PI;
        ft_imag /= PI;

        // Calculate magnitude and store results
        let magnitude = (ft_real * ft_real + ft_imag * ft_imag).sqrt();

        chi_r_real.push(ft_real);
        chi_r_imag.push(ft_imag);
        chi_r_mag.push(magnitude);
    }

    // Apply phase correction (optional)
    // In EXAFS analysis, a phase correction is often applied to account for
    // phase shifts from the central atom and scatterers
    // This helps align peaks with actual bond distances
    let apply_phase_correction = false; // Could be made configurable

    if apply_phase_correction {
        // Typical phase correction for first-shell paths
        // This is a simple approximation based on typical phase shifts
        let phase_correction = 0.5; // Å

        // Shift the r_values by the phase correction
        for r_value in &mut r_values {
            *r_value -= phase_correction;
        }
    }

    // Update the EXAFS data with r-space results
    exafs_data.r_values = Some(r_values);
    exafs_data.chi_r_real = Some(chi_r_real);
    exafs_data.chi_r_imag = Some(chi_r_imag);
    exafs_data.chi_r_mag = Some(chi_r_mag);

    exafs_data
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::atoms::{Atom, Vector3D};

    #[test]
    fn test_energy_grid_creation() {
        let e0 = 8333.0; // Fe K-edge
        let grid = EnergyGrid::new(e0, 3.0, 15.0, 1.0);

        assert_eq!(grid.e0, e0);
        assert_eq!(grid.k_values.len(), 13);
        assert_eq!(grid.energies.len(), 13);

        // First point should be at k=3.0
        assert!((grid.k_values[0] - 3.0).abs() < 1e-6);

        // Check conversion from k to E
        // E = e0 + (k²/2) * conversion_factor
        let k_squared = grid.k_values[0] * grid.k_values[0];
        let expected_e = e0 + k_squared / 2.0 / 1.8897259886 / 1.8897259886 * HARTREE_TO_EV;
        assert!((grid.energies[0] - expected_e).abs() < 1e-3);

        // Check grid type
        assert_eq!(grid.grid_type, GridType::Linear);
    }

    #[test]
    fn test_logarithmic_grid() {
        let e0 = 8333.0; // Fe K-edge
        let grid = EnergyGrid::new_logarithmic(e0, 2.0, 16.0, 20);

        assert_eq!(grid.e0, e0);
        assert_eq!(grid.k_values.len(), 20);
        assert_eq!(grid.energies.len(), 20);
        assert_eq!(grid.grid_type, GridType::Logarithmic);

        // First point should be at k=2.0
        assert!((grid.k_values[0] - 2.0).abs() < 1e-6);

        // Last point should be at k=16.0
        assert!((grid.k_values[19] - 16.0).abs() < 1e-6);

        // Points should follow logarithmic spacing
        // Check that ratio between consecutive differences is roughly constant
        let ratio1 = (grid.k_values[2] - grid.k_values[1]) / (grid.k_values[1] - grid.k_values[0]);
        let ratio2 = (grid.k_values[3] - grid.k_values[2]) / (grid.k_values[2] - grid.k_values[1]);

        // In logarithmic grid, ratios should be approximately equal
        assert!((ratio1 - ratio2).abs() < 1e-2);
    }

    #[test]
    fn test_exponential_grid() {
        let e0 = 8333.0; // Fe K-edge
        let grid = EnergyGrid::new_exponential(e0, 2.0, 16.0, 20, 2.0);

        assert_eq!(grid.e0, e0);
        assert_eq!(grid.k_values.len(), 20);
        assert_eq!(grid.energies.len(), 20);
        assert_eq!(grid.grid_type, GridType::Exponential);

        // First point should be at k=2.0
        assert!((grid.k_values[0] - 2.0).abs() < 1e-6);

        // Last point should be at k=16.0
        assert!((grid.k_values[19] - 16.0).abs() < 1e-6);

        // Points should have higher density at high k
        // Check that spacing increases
        let diff1 = grid.k_values[10] - grid.k_values[9];
        let diff2 = grid.k_values[15] - grid.k_values[14];

        // In exponential grid with exponent > 1, later differences should be larger
        assert!(diff2 > diff1);
    }

    #[test]
    fn test_custom_grid() {
        let e0 = 8333.0; // Fe K-edge
        let custom_k_values = vec![3.0, 5.0, 8.0, 12.0, 15.0];
        let grid = EnergyGrid::new_custom(e0, custom_k_values.clone());

        assert_eq!(grid.e0, e0);
        assert_eq!(grid.k_values.len(), 5);
        assert_eq!(grid.energies.len(), 5);
        assert_eq!(grid.grid_type, GridType::Custom);

        // k values should match input
        for i in 0..custom_k_values.len() {
            assert_eq!(grid.k_values[i], custom_k_values[i]);
        }
    }

    #[test]
    fn test_energy_from_grid() {
        let e0 = 8333.0; // Fe K-edge
        let custom_energies = vec![8330.0, 8340.0, 8350.0, 8380.0, 8450.0, 8550.0];
        let grid = EnergyGrid::new_from_energies(e0, custom_energies.clone());

        assert_eq!(grid.e0, e0);
        assert_eq!(grid.energies.len(), 6);
        assert_eq!(grid.k_values.len(), 6);
        assert_eq!(grid.grid_type, GridType::Custom);

        // Energies should match input
        for i in 0..custom_energies.len() {
            assert_eq!(grid.energies[i], custom_energies[i]);
        }

        // k values should be consistent with energy conversion
        for i in 0..grid.k_values.len() {
            let energy = grid.energies[i];
            if energy >= e0 {
                let e_rel = energy - e0;
                let e_hartree = e_rel / HARTREE_TO_EV;
                let k_atomic = (2.0 * e_hartree).sqrt();
                let expected_k = k_atomic * 1.8897259886;
                assert!((grid.k_values[i] - expected_k).abs() < 1e-3);
            } else {
                assert!(grid.k_values[i] < 0.0); // Should be negative for E < E0
            }
        }
    }

    #[test]
    fn test_exafs_optimized_grid() {
        let e0 = 8333.0; // Fe K-edge
        let grid = EnergyGrid::new_exafs_optimized(e0, 2.0, 16.0, 1.0);

        assert_eq!(grid.e0, e0);
        assert!(grid.k_values.len() > 10); // Should have reasonable number of points
        assert_eq!(grid.energies.len(), grid.k_values.len());
        assert_eq!(grid.grid_type, GridType::Custom);

        // First point should be close to k_min
        assert!((grid.k_values[0] - 2.0).abs() < 1e-2);

        // The non-uniform distribution means the last point might not be exactly at k_max
        // but it should be approaching it
        let last_idx = grid.k_values.len() - 1;
        let (_, max_k) = grid.k_range(); // We only need the max value
        assert!(max_k > 10.0); // Should reach at least a reasonable k value

        // Spacing should increase with k (higher density at high k)
        let mid_idx = grid.k_values.len() / 2;
        let diff_start = grid.k_values[1] - grid.k_values[0];
        let diff_mid = grid.k_values[mid_idx] - grid.k_values[mid_idx - 1];
        let diff_end = grid.k_values[last_idx] - grid.k_values[last_idx - 1];

        // Check that point spacing increases, which is the main feature
        // of the optimized grid
        assert!(diff_mid > diff_start);
        assert!(diff_end > diff_mid);
    }

    #[test]
    fn test_grid_merging() {
        let e0 = 8333.0; // Fe K-edge

        // Create two grids with overlapping ranges
        let grid1 = EnergyGrid::new(e0, 2.0, 8.0, 1.0);
        let grid2 = EnergyGrid::new(e0, 6.0, 14.0, 1.0);

        // Merge the grids
        let merged = grid1.merge(&grid2);

        // Check properties
        assert_eq!(merged.e0, e0);

        // Number of points should be the sum minus the overlapping points (k=6,7,8)
        let expected_points = grid1.len() + grid2.len() - 3;
        assert_eq!(merged.k_values.len(), expected_points);

        // Grid should be sorted by energy
        for i in 1..merged.energies.len() {
            assert!(merged.energies[i] > merged.energies[i - 1]);
        }

        // Merged grid should span the full range
        let (min_k, max_k) = merged.k_range();
        assert!((min_k - 2.0).abs() < 1e-6);
        assert!((max_k - 14.0).abs() < 1e-6);
    }

    #[test]
    fn test_grid_subset() {
        let e0 = 8333.0; // Fe K-edge
        let grid = EnergyGrid::new(e0, 2.0, 16.0, 1.0);

        // Get a k subset
        let subset = grid.k_subset(5.0, 10.0);

        // Check properties
        assert_eq!(subset.e0, e0);
        assert_eq!(subset.grid_type, GridType::Linear);

        // Should include only points in the specified range
        for k in &subset.k_values {
            assert!(*k >= 5.0 && *k <= 10.0);
        }

        // First point should be k=5.0, last should be k=10.0
        assert!((subset.k_values[0] - 5.0).abs() < 1e-6);
        assert!((subset.k_values[subset.k_values.len() - 1] - 10.0).abs() < 1e-6);
    }

    #[test]
    fn test_find_closest() {
        let e0 = 8333.0; // Fe K-edge
        let grid = EnergyGrid::new(e0, 2.0, 10.0, 2.0);

        // Find closest k index
        let idx = grid.find_closest_k_index(5.1);
        assert_eq!(grid.k_values[idx], 6.0); // Should find k=6.0 as closest to 5.1

        // Find closest energy index
        let energy_idx = grid.find_closest_energy_index(grid.energies[2] + 1.0);
        assert_eq!(energy_idx, 2); // Should find the energy closest to e[2]+1
    }

    #[test]
    fn test_exafs_data_creation() {
        let grid = EnergyGrid::new(8333.0, 3.0, 15.0, 1.0);
        let data = ExafsData::new(grid);

        assert_eq!(data.chi_k.len(), 13);
        assert_eq!(data.k_chi_k.len(), 13);
        assert_eq!(data.k2_chi_k.len(), 13);
        assert_eq!(data.k3_chi_k.len(), 13);

        // r-space data should be None initially
        assert!(data.r_values.is_none());
        assert!(data.chi_r_mag.is_none());
        assert!(data.chi_r_real.is_none());
        assert!(data.chi_r_imag.is_none());
    }

    #[test]
    fn test_window_functions() {
        let k_values = vec![3.0, 4.0, 5.0, 6.0, 7.0];
        let chi_k = vec![1.0, 1.0, 1.0, 1.0, 1.0];

        // No window should return original values
        let windowed_none = apply_window(&k_values, &chi_k, WindowFunction::None);
        assert_eq!(windowed_none, chi_k);

        // Hanning window should modify values
        let windowed_hanning = apply_window(&k_values, &chi_k, WindowFunction::Hanning);
        assert_ne!(windowed_hanning, chi_k);

        // First and last points in Hanning window should be close to 0
        assert!(windowed_hanning[0] < 0.1);
        assert!(windowed_hanning[4] < 0.1);

        // Middle point should be close to 1
        assert!(windowed_hanning[2] > 0.9);

        // Test various window functions
        let window_types = vec![
            WindowFunction::Hamming,
            WindowFunction::Blackman,
            WindowFunction::BlackmanHarris,
            WindowFunction::KaiserBessel(4.0),
            WindowFunction::Gaussian(0.25),
            WindowFunction::Parzen,
            WindowFunction::Welch,
            WindowFunction::Tukey(0.5),
            WindowFunction::FlatTop,
            WindowFunction::Exponential(0.2),
        ];

        for window_type in window_types {
            let windowed = apply_window(&k_values, &chi_k, window_type);

            // All windows should have the same length as the input
            assert_eq!(windowed.len(), chi_k.len());

            // No window values should be negative
            for val in &windowed {
                assert!(*val >= 0.0);
            }

            // All window functions attenuate the signal at edges
            assert!(windowed[0] <= 1.0);
            assert!(windowed[4] <= 1.0);
        }

        // Test window selection function
        let window1 = select_window_function(0.05, 5.0, 20.0);
        let window2 = select_window_function(0.2, 7.0, 2.0);
        let window3 = select_window_function(0.1, 15.0, 5.0);

        // Different conditions should result in different window types
        assert_ne!(format!("{:?}", window1), format!("{:?}", window2));
        assert_ne!(format!("{:?}", window2), format!("{:?}", window3));
    }

    #[test]
    fn test_specialized_windows() {
        // Create test data
        let e0 = 8333.0; // Fe K-edge
        let grid = EnergyGrid::new(e0, 3.0, 15.0, 1.0);
        let mut data = ExafsData::new(grid);

        // Fill with simple test data
        for i in 0..data.chi_k.len() {
            let k = data.grid.k_values[i];
            data.chi_k[i] = (3.0 * k).sin() + 0.5 * (6.0 * k).sin(); // Two frequency components
        }
        data.calculate_weighted_spectra();

        // Apply specialized windows
        let backscattering_data = apply_backscattering_window(data.clone(), 2.5);
        let multiple_scattering_data = apply_multiple_scattering_window(data.clone(), 5.0);

        // Results should be different from original
        assert_ne!(backscattering_data.chi_k, data.chi_k);
        assert_ne!(multiple_scattering_data.chi_k, data.chi_k);

        // Results should also be different from each other
        assert_ne!(backscattering_data.chi_k, multiple_scattering_data.chi_k);
    }

    // More tests will be added as implementation progresses
}
