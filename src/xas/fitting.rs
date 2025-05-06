/*
MIT License with FEFF10 Attribution

Copyright (c) 2025 Ameyanagi

Based on or developed using Distribution: FEFF10.0
Copyright (c) 2020 FEFF Project, University of Washington and SLAC National Accelerator Laboratory.
All rights reserved.
*/

//! EXAFS fitting utilities module
//!
//! This module provides functions and structures for fitting theoretical EXAFS
//! spectra to experimental data. It includes least-squares fitting algorithms,
//! parameter optimization, and error analysis.

use thiserror::Error;

use crate::atoms::{AtomicStructure, Result as AtomResult};
use crate::path::{Path, PathFinder, PathFinderConfig};
use crate::scattering::ScatteringResults;
use crate::xas::exafs::{
    calculate_exafs, fourier_transform, EnergyGrid, ExafsData, ExafsParameters, WindowFunction,
};

/// Errors that can occur during EXAFS fitting
#[derive(Error, Debug)]
pub enum FittingError {
    /// Error when parameters are out of allowed ranges
    #[error("Parameter {0} value {1} is outside allowed range {2}..{3}")]
    ParameterOutOfRange(String, f64, f64, f64),

    /// Error when fitting diverges or fails to converge
    #[error("Fitting failed to converge after {0} iterations")]
    FittingNotConverged(usize),

    /// Error when input data is invalid
    #[error("Invalid input data: {0}")]
    InvalidData(String),

    /// Error when fitting algorithm encounters a mathematical problem
    #[error("Mathematical error during fitting: {0}")]
    MathError(String),

    /// Error when reference to experimental data is missing
    #[error("No experimental data provided for fitting")]
    NoExperimentalData,
}

/// Type of parameter in EXAFS fitting
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ParameterType {
    /// S0² - amplitude reduction factor
    AmplitudeFactor,
    /// E0 - energy reference shift (eV)
    E0Shift,
    /// ΔR - distance correction (Å)
    DistanceCorrection,
    /// σ² - Debye-Waller factor (Å²)
    DebyeWaller,
    /// Third cumulant - anharmonicity parameter (Å³)
    ThirdCumulant,
    /// Fourth cumulant - anharmonicity parameter (Å⁴)
    FourthCumulant,
    /// Coordination number scale factor
    CoordinationNumber,
    /// Custom parameter for user-defined functions
    Custom,
}

/// Parameter for EXAFS fitting
#[derive(Debug, Clone)]
pub struct FittingParameter {
    /// Name of the parameter
    pub name: String,
    /// Current value of the parameter
    pub value: f64,
    /// Error/uncertainty in the parameter value
    pub error: Option<f64>,
    /// Minimum allowed value
    pub min: f64,
    /// Maximum allowed value
    pub max: f64,
    /// Whether the parameter is allowed to vary during fitting
    pub vary: bool,
    /// Type of parameter
    pub param_type: ParameterType,
    /// Group ID for correlated parameters (if any)
    pub group_id: Option<String>,
    /// Expression for constrained parameters
    pub constraint: Option<String>,
}

impl FittingParameter {
    /// Create a new fitting parameter
    pub fn new(
        name: &str,
        value: f64,
        min: f64,
        max: f64,
        vary: bool,
        param_type: ParameterType,
    ) -> Self {
        Self {
            name: name.to_string(),
            value,
            error: None,
            min,
            max,
            vary,
            param_type,
            group_id: None,
            constraint: None,
        }
    }

    /// Create a new S0² parameter
    pub fn s02(value: f64, vary: bool) -> Self {
        Self::new("s02", value, 0.5, 1.2, vary, ParameterType::AmplitudeFactor)
    }

    /// Create a new E0 shift parameter
    pub fn e0_shift(value: f64, vary: bool) -> Self {
        Self::new("e0_shift", value, -10.0, 10.0, vary, ParameterType::E0Shift)
    }

    /// Create a new distance correction (ΔR) parameter
    pub fn delta_r(value: f64, vary: bool) -> Self {
        Self::new(
            "delta_r",
            value,
            -0.5,
            0.5,
            vary,
            ParameterType::DistanceCorrection,
        )
    }

    /// Create a new Debye-Waller factor (σ²) parameter
    pub fn sigma_squared(value: f64, vary: bool) -> Self {
        Self::new(
            "sigma_squared",
            value,
            0.0,
            0.025,
            vary,
            ParameterType::DebyeWaller,
        )
    }

    /// Create a new third cumulant parameter
    pub fn third_cumulant(value: f64, vary: bool) -> Self {
        Self::new(
            "third_cumulant",
            value,
            -0.001,
            0.001,
            vary,
            ParameterType::ThirdCumulant,
        )
    }

    /// Create a new coordination number parameter
    pub fn coordination_number(value: f64, vary: bool) -> Self {
        Self::new(
            "coordination_number",
            value,
            0.5,
            2.0,
            vary,
            ParameterType::CoordinationNumber,
        )
    }

    /// Check if parameter value is within allowed range
    pub fn is_valid(&self) -> Result<(), FittingError> {
        if self.value < self.min || self.value > self.max {
            return Err(FittingError::ParameterOutOfRange(
                self.name.clone(),
                self.value,
                self.min,
                self.max,
            ));
        }
        Ok(())
    }
}

/// Configuration for path-specific parameters in EXAFS fitting
#[derive(Debug, Clone)]
pub struct PathParameterConfig {
    /// Path filter pattern (e.g., "Fe-O", "*.2.*" for second shell)
    pub path_pattern: String,
    /// Whether this path is included in the fit
    pub include: bool,
    /// ΔR parameter for this path
    pub delta_r: Option<FittingParameter>,
    /// σ² parameter for this path
    pub sigma_squared: Option<FittingParameter>,
    /// Third cumulant parameter for this path
    pub third_cumulant: Option<FittingParameter>,
    /// Fourth cumulant parameter for this path
    pub fourth_cumulant: Option<FittingParameter>,
    /// Coordination number scaling for this path
    pub n_scale: Option<FittingParameter>,
}

/// EXAFS fitting model configuration
#[derive(Clone)]
pub struct FittingModel {
    /// Structure model for EXAFS calculation
    pub structure: AtomicStructure,
    /// Global fitting parameters
    pub global_parameters: Vec<FittingParameter>,
    /// Path-specific parameters
    pub path_parameters: Vec<PathParameterConfig>,
    /// Current fit quality metrics
    pub fit_quality: Option<FitQuality>,
    /// Experimental data to fit against
    pub experimental_data: Option<ExafsData>,
    /// k-range for fitting
    pub k_range: (f64, f64),
    /// k-weight for fitting
    pub k_weight: usize,
    /// Whether to fit in k-space (true) or r-space (false)
    pub fit_in_k_space: bool,
    /// r-range for r-space fitting
    pub r_range: Option<(f64, f64)>,
    /// Window function for Fourier transform
    pub window_function: WindowFunction,
    /// Maximum number of iterations for fitting
    pub max_iterations: usize,
    /// Convergence tolerance
    pub tolerance: f64,
}

impl FittingModel {
    /// Create a new fitting model
    pub fn new(structure: AtomicStructure) -> Self {
        Self {
            structure,
            global_parameters: vec![
                FittingParameter::s02(1.0, true),
                FittingParameter::e0_shift(0.0, true),
            ],
            path_parameters: Vec::new(),
            fit_quality: None,
            experimental_data: None,
            k_range: (3.0, 12.0),
            k_weight: 2,
            fit_in_k_space: false,
            r_range: Some((1.0, 6.0)),
            window_function: WindowFunction::Hanning,
            max_iterations: 50,
            tolerance: 1e-5,
        }
    }

    /// Set the experimental data to fit against
    pub fn set_experimental_data(&mut self, data: ExafsData) {
        self.experimental_data = Some(data);
    }

    /// Add a global parameter to the model
    pub fn add_global_parameter(&mut self, parameter: FittingParameter) {
        self.global_parameters.push(parameter);
    }

    /// Add a path-specific parameter configuration
    pub fn add_path_parameter(&mut self, config: PathParameterConfig) {
        self.path_parameters.push(config);
    }

    /// Find a global parameter by name
    pub fn get_global_parameter(&self, name: &str) -> Option<&FittingParameter> {
        self.global_parameters.iter().find(|p| p.name == name)
    }

    /// Get a mutable reference to a global parameter by name
    pub fn get_global_parameter_mut(&mut self, name: &str) -> Option<&mut FittingParameter> {
        self.global_parameters.iter_mut().find(|p| p.name == name)
    }

    /// Get all parameters that can vary during fitting
    pub fn get_variable_parameters(&self) -> Vec<&FittingParameter> {
        let mut variable = self
            .global_parameters
            .iter()
            .filter(|p| p.vary)
            .collect::<Vec<_>>();

        for path_param in &self.path_parameters {
            if let Some(ref p) = path_param.delta_r {
                if p.vary {
                    variable.push(p);
                }
            }
            if let Some(ref p) = path_param.sigma_squared {
                if p.vary {
                    variable.push(p);
                }
            }
            if let Some(ref p) = path_param.third_cumulant {
                if p.vary {
                    variable.push(p);
                }
            }
            if let Some(ref p) = path_param.fourth_cumulant {
                if p.vary {
                    variable.push(p);
                }
            }
            if let Some(ref p) = path_param.n_scale {
                if p.vary {
                    variable.push(p);
                }
            }
        }

        variable
    }

    /// Check if a path matches a path pattern
    fn path_matches_pattern(&self, path: &Path, pattern: &str) -> bool {
        // Simple pattern matching for now
        // In a more sophisticated implementation, this would handle wildcards and regex

        // Extract the path notation (e.g., "Fe-O-Fe")
        let mut path_notation = String::new();

        for leg in &path.legs {
            let from_atom = self.structure.atom(leg.from_atom).unwrap();
            let to_atom = self.structure.atom(leg.to_atom).unwrap();

            if path_notation.is_empty() {
                // For the first leg, add the starting atom
                path_notation.push_str(&format!("{}", from_atom.atomic_number()));
            }

            // Add the destination atom
            path_notation.push('-');
            path_notation.push_str(&format!("{}", to_atom.atomic_number()));
        }

        // Check if the pattern matches
        if pattern == "*" {
            return true;
        }

        if pattern.contains('*') {
            // Simple wildcard matching
            let parts: Vec<&str> = pattern.split('*').collect();
            let mut remaining = path_notation.as_str();

            for part in parts {
                if part.is_empty() {
                    continue;
                }

                if let Some(idx) = remaining.find(part) {
                    remaining = &remaining[idx + part.len()..];
                } else {
                    return false;
                }
            }

            true
        } else {
            // Exact matching
            path_notation == pattern
        }
    }

    /// Apply the current parameter values to calculate a theoretical EXAFS spectrum
    pub fn calculate_theoretical_spectrum(
        &self,
        phase_shifts: &ScatteringResults,
    ) -> AtomResult<ExafsData> {
        // Get current parameter values
        let s02 = self
            .get_global_parameter("s02")
            .map(|p| p.value)
            .unwrap_or(1.0);

        let e0_shift = self
            .get_global_parameter("e0_shift")
            .map(|p| p.value)
            .unwrap_or(0.0);

        // Create a copy of the structure for modification
        let structure = self.structure.clone();

        // Apply e0 shift to the energy reference
        let experimental_e0 = if let Some(ref exp_data) = self.experimental_data {
            exp_data.grid.e0
        } else {
            0.0
        };

        // Create energy grid with updated e0
        let adjusted_e0 = experimental_e0 + e0_shift;
        let k_min = self.k_range.0;
        let k_max = self.k_range.1;
        let k_step = 0.05; // Default step size

        let energy_grid = EnergyGrid::new(adjusted_e0, k_min, k_max, k_step);

        // Create EXAFS parameters
        let params = ExafsParameters {
            s02,
            r_max: 8.0, // Sufficiently large to include all paths
            min_importance: 0.001,
            max_legs: 4,
            use_debye_waller: true,
            temperature: 300.0,
            energy_range: energy_grid,
        };

        // Find paths
        let absorber_index = structure.central_atom_index().ok_or_else(|| {
            crate::atoms::errors::AtomError::CalculationError(
                "Central atom not defined".to_string(),
            )
        })?;

        let path_config = PathFinderConfig {
            max_path_length: params.r_max,
            max_paths: 100,
            max_legs: params.max_legs,
            importance_threshold: params.min_importance,
            cluster_paths: true,
            unique_scatterers_only: true,
        };

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
                // Check if this path should be included based on path parameters
                let mut include_path = true;
                let mut path_delta_r = 0.0;
                let mut path_sigma_squared = 0.0;
                let mut path_third_cumulant = 0.0;
                let mut path_n_scale = 1.0;

                // Apply path-specific parameters if they exist
                for path_param in &self.path_parameters {
                    if self.path_matches_pattern(path, &path_param.path_pattern) {
                        include_path = path_param.include;

                        if !include_path {
                            break;
                        }

                        // Apply path-specific parameters
                        if let Some(ref param) = path_param.delta_r {
                            path_delta_r = param.value;
                        }

                        if let Some(ref param) = path_param.sigma_squared {
                            path_sigma_squared = param.value;
                        }

                        if let Some(ref param) = path_param.third_cumulant {
                            path_third_cumulant = param.value;
                        }

                        if let Some(ref param) = path_param.n_scale {
                            path_n_scale = param.value;
                        }
                    }
                }

                if !include_path {
                    continue;
                }

                // Calculate path contribution with modified parameters
                let mut modified_path = path.clone();

                // Apply ΔR correction
                modified_path.total_length += path_delta_r;

                // Apply coordination number scaling
                modified_path.degeneracy = (modified_path.degeneracy as f64 * path_n_scale) as u32;

                // Calculate path contribution with our custom implementation
                // (since the original is private)
                let path_chi = calculate_path_contribution_with_fitting(
                    &modified_path,
                    k,
                    energy,
                    &structure,
                    phase_shifts,
                    &params,
                    path_sigma_squared,
                    path_third_cumulant,
                )?;

                chi_k += path_chi;
            }

            exafs_data.chi_k[i] = chi_k;
        }

        // Calculate k-weighted spectra
        exafs_data.calculate_weighted_spectra();

        // If fitting in r-space, perform Fourier transform
        if !self.fit_in_k_space {
            let r_min = self.r_range.unwrap_or((0.0, 8.0)).0;
            let r_max = self.r_range.unwrap_or((0.0, 8.0)).1;
            let dr = 0.02; // Default r-space step

            exafs_data = fourier_transform(
                exafs_data,
                self.window_function,
                self.k_weight,
                r_min,
                r_max,
                dr,
            );
        }

        Ok(exafs_data)
    }

    /// Perform the fitting to match theoretical and experimental EXAFS
    pub fn fit(&mut self, phase_shifts: &ScatteringResults) -> Result<FitQuality, FittingError> {
        // Check if we have experimental data
        let _experimental_data = match &self.experimental_data {
            Some(data) => data,
            None => return Err(FittingError::NoExperimentalData),
        };

        // Check if any parameters can vary
        // Get parameter names and info up front to avoid borrow checker issues
        let parameter_info: Vec<(String, f64, f64, f64)> = self
            .global_parameters
            .iter()
            .filter(|p| p.vary)
            .map(|p| (p.name.clone(), p.value, p.min, p.max))
            .collect();

        // Also collect path parameter info
        let mut path_param_info = Vec::new();
        for path_param in &self.path_parameters {
            if let Some(ref p) = path_param.delta_r {
                if p.vary {
                    path_param_info.push((p.name.clone(), p.value, p.min, p.max));
                }
            }
            if let Some(ref p) = path_param.sigma_squared {
                if p.vary {
                    path_param_info.push((p.name.clone(), p.value, p.min, p.max));
                }
            }
            if let Some(ref p) = path_param.third_cumulant {
                if p.vary {
                    path_param_info.push((p.name.clone(), p.value, p.min, p.max));
                }
            }
            if let Some(ref p) = path_param.n_scale {
                if p.vary {
                    path_param_info.push((p.name.clone(), p.value, p.min, p.max));
                }
            }
        }

        // Combine all parameters
        let all_param_info = [parameter_info, path_param_info].concat();

        if all_param_info.is_empty() {
            return Err(FittingError::InvalidData(
                "No variable parameters specified for fitting".to_string(),
            ));
        }

        // Define the objective function to minimize
        let mut current_quality = self.calculate_fit_quality(phase_shifts)?;
        let mut best_quality = current_quality.clone();
        let mut converged = false;
        let mut iteration = 0;

        // Main optimization loop
        while !converged && iteration < self.max_iterations {
            iteration += 1;
            let mut improvement = false;

            // Try varying each parameter
            for (param_name, current_value, param_min, param_max) in &all_param_info {
                // Try a smaller value
                let smaller = (current_value - current_value * 0.05).max(*param_min);
                self.update_parameter_value(param_name, smaller)?;

                // Calculate quality with smaller value
                let quality_smaller = self.calculate_fit_quality(phase_shifts)?;

                // Try a larger value
                let larger = (current_value + current_value * 0.05).min(*param_max);
                self.update_parameter_value(param_name, larger)?;

                // Calculate quality with larger value
                let quality_larger = self.calculate_fit_quality(phase_shifts)?;

                // Restore original value
                self.update_parameter_value(param_name, *current_value)?;

                // Check if either direction improves the fit
                if quality_smaller.r_factor < current_quality.r_factor {
                    self.update_parameter_value(param_name, smaller)?;
                    current_quality = quality_smaller;
                    improvement = true;
                } else if quality_larger.r_factor < current_quality.r_factor {
                    self.update_parameter_value(param_name, larger)?;
                    current_quality = quality_larger;
                    improvement = true;
                }

                // Update best quality if needed
                if current_quality.r_factor < best_quality.r_factor {
                    best_quality = current_quality.clone();
                }
            }

            // Check for convergence
            if !improvement || current_quality.r_factor < self.tolerance {
                converged = true;
            }
        }

        // If we didn't converge, return an error
        if !converged {
            return Err(FittingError::FittingNotConverged(iteration));
        }

        // Store the fit quality in the model
        self.fit_quality = Some(best_quality.clone());

        Ok(best_quality)
    }

    /// Update a parameter value by name
    fn update_parameter_value(&mut self, name: &str, value: f64) -> Result<(), FittingError> {
        // First check global parameters
        for param in &mut self.global_parameters {
            if param.name == name {
                param.value = value;
                return param.is_valid();
            }
        }

        // Then check path parameters
        for path_param in &mut self.path_parameters {
            if let Some(ref mut param) = path_param.delta_r {
                if param.name == name {
                    param.value = value;
                    return param.is_valid();
                }
            }

            if let Some(ref mut param) = path_param.sigma_squared {
                if param.name == name {
                    param.value = value;
                    return param.is_valid();
                }
            }

            if let Some(ref mut param) = path_param.third_cumulant {
                if param.name == name {
                    param.value = value;
                    return param.is_valid();
                }
            }

            if let Some(ref mut param) = path_param.fourth_cumulant {
                if param.name == name {
                    param.value = value;
                    return param.is_valid();
                }
            }

            if let Some(ref mut param) = path_param.n_scale {
                if param.name == name {
                    param.value = value;
                    return param.is_valid();
                }
            }
        }

        Err(FittingError::InvalidData(format!(
            "Parameter {} not found",
            name
        )))
    }

    /// Calculate the fit quality metrics between theoretical and experimental data
    pub fn calculate_fit_quality(
        &self,
        phase_shifts: &ScatteringResults,
    ) -> Result<FitQuality, FittingError> {
        // Check if we have experimental data
        let experimental_data = match &self.experimental_data {
            Some(data) => data,
            None => return Err(FittingError::NoExperimentalData),
        };

        // Calculate theoretical data
        let theoretical_data = match self.calculate_theoretical_spectrum(phase_shifts) {
            Ok(data) => data,
            Err(e) => return Err(FittingError::MathError(e.to_string())),
        };

        // Prepare data for comparison
        let (exp_x, exp_y, theo_y) = if self.fit_in_k_space {
            // Compare in k-space
            let k_values = &experimental_data.grid.k_values;
            let exp_chi = match self.k_weight {
                0 => &experimental_data.chi_k,
                1 => &experimental_data.k_chi_k,
                2 => &experimental_data.k2_chi_k,
                3 => &experimental_data.k3_chi_k,
                _ => &experimental_data.k2_chi_k, // Default to k²-weighted
            };

            let theo_chi = match self.k_weight {
                0 => &theoretical_data.chi_k,
                1 => &theoretical_data.k_chi_k,
                2 => &theoretical_data.k2_chi_k,
                3 => &theoretical_data.k3_chi_k,
                _ => &theoretical_data.k2_chi_k, // Default to k²-weighted
            };

            // Interpolate if necessary to match k points
            let matched_theo_y =
                if theoretical_data.grid.k_values == experimental_data.grid.k_values {
                    theo_chi.clone()
                } else {
                    // Simple linear interpolation
                    let mut interpolated = Vec::with_capacity(k_values.len());

                    for &k in k_values {
                        if k < theoretical_data.grid.k_values[0]
                            || k > *theoretical_data.grid.k_values.last().unwrap()
                        {
                            interpolated.push(0.0);
                            continue;
                        }

                        // Find the two nearest points in theoretical data
                        let mut idx = 0;
                        while idx < theoretical_data.grid.k_values.len() - 1
                            && theoretical_data.grid.k_values[idx + 1] < k
                        {
                            idx += 1;
                        }

                        if idx == theoretical_data.grid.k_values.len() - 1 {
                            interpolated.push(theo_chi[idx]);
                        } else {
                            let k1 = theoretical_data.grid.k_values[idx];
                            let k2 = theoretical_data.grid.k_values[idx + 1];
                            let y1 = theo_chi[idx];
                            let y2 = theo_chi[idx + 1];

                            let t = (k - k1) / (k2 - k1);
                            let y = y1 + t * (y2 - y1);
                            interpolated.push(y);
                        }
                    }

                    interpolated
                };

            (k_values.clone(), exp_chi.clone(), matched_theo_y)
        } else {
            // Compare in r-space
            let r_min = self.r_range.unwrap_or((0.0, 8.0)).0;
            let r_max = self.r_range.unwrap_or((0.0, 8.0)).1;

            if experimental_data.r_values.is_none() || experimental_data.chi_r_mag.is_none() {
                return Err(FittingError::InvalidData(
                    "Experimental data does not contain r-space information".to_string(),
                ));
            }

            if theoretical_data.r_values.is_none() || theoretical_data.chi_r_mag.is_none() {
                return Err(FittingError::InvalidData(
                    "Theoretical data does not contain r-space information".to_string(),
                ));
            }

            let exp_r = experimental_data.r_values.as_ref().unwrap();
            let exp_chi_r = experimental_data.chi_r_mag.as_ref().unwrap();
            let theo_r = theoretical_data.r_values.as_ref().unwrap();
            let theo_chi_r = theoretical_data.chi_r_mag.as_ref().unwrap();

            // Filter by r-range and interpolate theoretical data to match experimental r points
            let mut filtered_exp_r = Vec::new();
            let mut filtered_exp_chi = Vec::new();
            let mut filtered_theo_chi = Vec::new();

            for (i, &r) in exp_r.iter().enumerate() {
                if r >= r_min && r <= r_max {
                    filtered_exp_r.push(r);
                    filtered_exp_chi.push(exp_chi_r[i]);

                    // Find the matching point in theoretical data using interpolation
                    if r < theo_r[0] || r > *theo_r.last().unwrap() {
                        filtered_theo_chi.push(0.0);
                    } else {
                        // Find the two nearest points
                        let mut idx = 0;
                        while idx < theo_r.len() - 1 && theo_r[idx + 1] < r {
                            idx += 1;
                        }

                        if idx == theo_r.len() - 1 {
                            filtered_theo_chi.push(theo_chi_r[idx]);
                        } else {
                            let r1 = theo_r[idx];
                            let r2 = theo_r[idx + 1];
                            let y1 = theo_chi_r[idx];
                            let y2 = theo_chi_r[idx + 1];

                            let t = (r - r1) / (r2 - r1);
                            let y = y1 + t * (y2 - y1);
                            filtered_theo_chi.push(y);
                        }
                    }
                }
            }

            (filtered_exp_r, filtered_exp_chi, filtered_theo_chi)
        };

        // Calculate fit quality metrics
        if exp_x.is_empty() {
            return Err(FittingError::InvalidData(
                "No data points in fitting range".to_string(),
            ));
        }

        // Calculate R-factor
        let mut sum_squared_diff = 0.0;
        let mut sum_squared_exp = 0.0;

        for i in 0..exp_y.len() {
            let diff = exp_y[i] - theo_y[i];
            sum_squared_diff += diff * diff;
            sum_squared_exp += exp_y[i] * exp_y[i];
        }

        let r_factor = if sum_squared_exp > 0.0 {
            (sum_squared_diff / sum_squared_exp).sqrt()
        } else {
            1.0 // If experimental data is all zeros, r-factor is meaningless
        };

        // Calculate reduced chi-squared
        let n_points = exp_y.len();
        let n_parameters = self.get_variable_parameters().len();
        let degrees_of_freedom = n_points - n_parameters;

        let reduced_chi_squared = if degrees_of_freedom > 0 {
            sum_squared_diff / degrees_of_freedom as f64
        } else {
            f64::INFINITY
        };

        // Calculate Akaike Information Criterion
        let aic =
            n_points as f64 * (sum_squared_diff / n_points as f64).ln() + 2.0 * n_parameters as f64;

        // Calculate Bayesian Information Criterion
        let bic = n_points as f64 * (sum_squared_diff / n_points as f64).ln()
            + n_parameters as f64 * (n_points as f64).ln();

        Ok(FitQuality {
            r_factor,
            reduced_chi_squared,
            aic,
            bic,
            n_points,
            n_parameters,
            residuals: exp_y
                .iter()
                .zip(theo_y.iter())
                .map(|(&e, &t)| e - t)
                .collect(),
        })
    }
}

/// Quality metrics for an EXAFS fit
#[derive(Debug, Clone)]
pub struct FitQuality {
    /// R-factor (measure of misfit, lower is better)
    pub r_factor: f64,
    /// Reduced chi-squared (χ²/ν)
    pub reduced_chi_squared: f64,
    /// Akaike Information Criterion (AIC)
    pub aic: f64,
    /// Bayesian Information Criterion (BIC)
    pub bic: f64,
    /// Number of data points
    pub n_points: usize,
    /// Number of parameters
    pub n_parameters: usize,
    /// Residuals (experimental - theoretical)
    pub residuals: Vec<f64>,
}

impl FitQuality {
    /// Get parameter errors from the covariance matrix
    pub fn parameter_errors(&self, _covariance_matrix: &[Vec<f64>]) -> Vec<f64> {
        // For a basic implementation, just return the diagonal elements
        // of the covariance matrix (the variances)
        _covariance_matrix
            .iter()
            .enumerate()
            .map(|(i, row)| row[i].sqrt())
            .collect()
    }

    /// Calculate the correlation matrix
    pub fn correlation_matrix(&self, covariance_matrix: &[Vec<f64>]) -> Vec<Vec<f64>> {
        let n = covariance_matrix.len();
        let mut correlation = vec![vec![0.0; n]; n];

        for i in 0..n {
            for j in 0..n {
                let cov_ii = covariance_matrix[i][i];
                let cov_jj = covariance_matrix[j][j];

                if cov_ii > 0.0 && cov_jj > 0.0 {
                    correlation[i][j] = covariance_matrix[i][j] / (cov_ii * cov_jj).sqrt();
                } else {
                    correlation[i][j] = 0.0;
                }
            }
        }

        correlation
    }
}

/// Calculate the contribution of a single path to the EXAFS spectrum with cumulant expansion
///
/// This extends the basic path contribution calculation to include higher order cumulants
/// for anharmonic effects and non-Gaussian pair distributions.
///
/// # Arguments
///
/// * `path` - The scattering path
/// * `k` - Wavenumber in Å^-1
/// * `energy` - Energy in eV
/// * `structure` - Atomic structure
/// * `phase_shifts` - Scattering results containing phase shifts
/// * `params` - EXAFS calculation parameters
/// * `sigma_squared` - Debye-Waller factor (σ²) for this path
/// * `third_cumulant` - Third cumulant (C₃) for anharmonicity
///
/// # Returns
///
/// The path's contribution to chi(k)
fn calculate_path_contribution_with_fitting(
    path: &Path,
    k: &f64,
    _energy: f64,
    structure: &AtomicStructure,
    phase_shifts: &ScatteringResults,
    params: &ExafsParameters,
    sigma_squared: f64,
    third_cumulant: f64,
) -> AtomResult<f64> {
    // Get path properties
    let path_length = path.total_length;
    let degeneracy = path.degeneracy as f64;
    let legs = &path.legs;

    // Apply S0² global amplitude factor (accounts for many-body effects)
    let s02 = params.s02;

    // Calculate amplitude factor based on path geometry and scattering properties

    // Calculate scattering amplitude product
    let mut scattering_amplitude = 1.0;
    let mut total_phase_shift = 0.0;

    // Central atom phase shift
    let central_atom_idx = structure.central_atom_index().unwrap();
    let central_atom = structure.atom(central_atom_idx).unwrap();
    let central_potential_idx = central_atom.potential_type() as usize;

    // For the absorbing atom, contribution is from l=1 (p orbital) for K-edge
    let l_channel = 1usize; // l=1 for K-edge (dipole selection rules)

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

                // Phase shift contribution
                atom_amplitude += weight * phase.im;
                atom_phase += weight * phase.re;
            }

            // Normalize by sum of weights
            let norm = (max_l * 2 + 2) as f64; // Σ(2l+1) from l=0 to max_l
            atom_amplitude /= norm;
            atom_phase /= norm;
        }

        // Z-dependent correction
        let z_factor = (scattering_atom.atomic_number() as f64).sqrt() / 5.0;
        atom_amplitude *= z_factor;

        // Add to total
        scattering_amplitude *= atom_amplitude;
        total_phase_shift += atom_phase;
    }

    // Geometry factor - 1/kR²
    let geometry_factor = 1.0 / (k * path_length * path_length);

    // Use provided sigma_squared instead of calculating from temperature
    let thermal_factor = f64::exp(-2.0 * sigma_squared * k * k);

    // Mean free path damping
    let lambda_0 = 4.0; // Baseline MFP in Å
    let lambda_slope = 0.6; // How fast MFP increases with k
    let lambda = lambda_0 + lambda_slope * k; // Mean free path in Å

    let mfp_factor = f64::exp(-2.0 * path_length / lambda);

    // The full path phase includes standard oscillation plus cumulant corrections
    let standard_phase = 2.0 * k * path_length + total_phase_shift;

    // Add third cumulant correction to the phase
    // In the cumulant expansion, C₃ term modifies the phase: -4/3 * k³ * C₃
    let cumulant_phase_correction = -4.0 / 3.0 * k.powi(3) * third_cumulant;
    let full_phase = standard_phase + cumulant_phase_correction;

    // Combine all factors to calculate the EXAFS contribution from this path
    let chi = s02
        * degeneracy
        * scattering_amplitude
        * thermal_factor
        * mfp_factor
        * geometry_factor
        * full_phase.sin();

    Ok(chi)
}

/// Import experimental EXAFS data from file
///
/// Supports various file formats commonly used in EXAFS analysis
///
/// # Arguments
///
/// * `file_path` - Path to the experimental data file
/// * `file_type` - Type of file (e.g., "athena", "ascii", "viper")
/// * `e0` - Reference energy (edge position) in eV
///
/// # Returns
///
/// The imported experimental data as ExafsData
pub fn import_experimental_data(
    _file_path: &str,
    _file_type: &str,
    _e0: f64,
) -> Result<ExafsData, FittingError> {
    // This is a placeholder for a more complete implementation
    // The actual implementation would handle various file formats

    Err(FittingError::InvalidData(
        "Data import not yet implemented".to_string(),
    ))
}

/// Export EXAFS fitting results to file
///
/// # Arguments
///
/// * `model` - The fitted EXAFS model
/// * `file_path` - Path to save results
/// * `format` - Output format (e.g., "json", "csv", "txt")
///
/// # Returns
///
/// Result indicating success or failure
pub fn export_fitting_results(
    _model: &FittingModel,
    _file_path: &str,
    _format: &str,
) -> Result<(), FittingError> {
    // This is a placeholder for a more complete implementation
    // The actual implementation would output the fitting results
    // in the requested format

    Err(FittingError::InvalidData(
        "Results export not yet implemented".to_string(),
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::atoms::{Atom, Vector3D};

    #[test]
    fn test_fitting_parameter_creation() {
        let s02 = FittingParameter::s02(0.9, true);
        assert_eq!(s02.name, "s02");
        assert_eq!(s02.value, 0.9);
        assert_eq!(s02.vary, true);
        assert_eq!(s02.param_type, ParameterType::AmplitudeFactor);

        let dr = FittingParameter::delta_r(0.02, false);
        assert_eq!(dr.name, "delta_r");
        assert_eq!(dr.value, 0.02);
        assert_eq!(dr.vary, false);
        assert_eq!(dr.param_type, ParameterType::DistanceCorrection);
    }

    #[test]
    fn test_parameter_validation() {
        let mut param = FittingParameter::s02(0.9, true);
        assert!(param.is_valid().is_ok());

        // Test out of range
        param.value = 1.5;
        assert!(param.is_valid().is_err());
    }

    #[test]
    fn test_fitting_model_creation() {
        // Create a simple structure
        let fe_potential = crate::atoms::PotentialType::new(0, 26).unwrap();
        let fe_atom = Atom::new(26, Vector3D::new(0.0, 0.0, 0.0), 0).unwrap();

        let mut structure = AtomicStructure::new();
        structure.add_potential_type(fe_potential);
        let fe_idx = structure.add_atom(fe_atom);
        structure.set_central_atom(fe_idx).unwrap();

        let model = FittingModel::new(structure);

        // Check default values
        assert_eq!(model.global_parameters.len(), 2); // s02 and e0_shift
        assert_eq!(model.path_parameters.len(), 0);
        assert_eq!(model.k_range, (3.0, 12.0));
        assert_eq!(model.k_weight, 2);
        assert_eq!(model.fit_in_k_space, false); // Default to r-space
    }

    #[test]
    fn test_get_variable_parameters() {
        // Create a simple structure
        let fe_potential = crate::atoms::PotentialType::new(0, 26).unwrap();
        let fe_atom = Atom::new(26, Vector3D::new(0.0, 0.0, 0.0), 0).unwrap();

        let mut structure = AtomicStructure::new();
        structure.add_potential_type(fe_potential);
        let fe_idx = structure.add_atom(fe_atom);
        structure.set_central_atom(fe_idx).unwrap();

        let mut model = FittingModel::new(structure);

        // By default, s02 and e0_shift are variable
        let variable_params = model.get_variable_parameters();
        assert_eq!(variable_params.len(), 2);

        // Add a fixed parameter
        model.add_global_parameter(FittingParameter::new(
            "fixed_param",
            1.0,
            0.0,
            2.0,
            false,
            ParameterType::Custom,
        ));

        // Should still be 2 variable parameters
        let variable_params = model.get_variable_parameters();
        assert_eq!(variable_params.len(), 2);

        // Add a variable parameter
        model.add_global_parameter(FittingParameter::new(
            "var_param",
            1.0,
            0.0,
            2.0,
            true,
            ParameterType::Custom,
        ));

        // Now 3 variable parameters
        let variable_params = model.get_variable_parameters();
        assert_eq!(variable_params.len(), 3);
    }
}
