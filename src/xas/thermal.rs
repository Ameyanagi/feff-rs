/*
MIT License with FEFF10 Attribution

Copyright (c) 2025 Ameyanagi

Based on or developed using Distribution: FEFF10.0
Copyright (c) 2020 FEFF Project, University of Washington and SLAC National Accelerator Laboratory.
All rights reserved.
*/

//! Thermal parameters for XAS calculations
//!
//! This module defines structures for thermal parameters used in
//! XANES and EXAFS calculations to account for temperature effects.

use serde::{Deserialize, Serialize};

/// Thermal parameters for XAS calculations
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ThermalParameters {
    /// Temperature in Kelvin
    pub temperature: f64,
    /// Type of thermal model: "debye" or "einstein"
    pub model_type: String,
    /// Debye temperature in Kelvin (for Debye model)
    pub debye_temperature: f64,
    /// Einstein frequency in meV (for Einstein model)
    pub einstein_frequency: Option<f64>,
}

impl Default for ThermalParameters {
    fn default() -> Self {
        Self {
            temperature: 300.0, // Room temperature
            model_type: "debye".to_string(),
            debye_temperature: 300.0, // Default Debye temperature
            einstein_frequency: None,
        }
    }
}

impl ThermalParameters {
    /// Create new thermal parameters for a Debye model
    ///
    /// # Arguments
    ///
    /// * `temperature` - Temperature in Kelvin
    /// * `debye_temperature` - Debye temperature in Kelvin
    ///
    /// # Returns
    ///
    /// New thermal parameters
    pub fn new_debye(temperature: f64, debye_temperature: f64) -> Self {
        Self {
            temperature,
            model_type: "debye".to_string(),
            debye_temperature,
            einstein_frequency: None,
        }
    }

    /// Create new thermal parameters for an Einstein model
    ///
    /// # Arguments
    ///
    /// * `temperature` - Temperature in Kelvin
    /// * `einstein_frequency` - Einstein frequency in meV
    ///
    /// # Returns
    ///
    /// New thermal parameters
    pub fn new_einstein(temperature: f64, einstein_frequency: f64) -> Self {
        Self {
            temperature,
            model_type: "einstein".to_string(),
            debye_temperature: 0.0, // Not used for Einstein model
            einstein_frequency: Some(einstein_frequency),
        }
    }
}
