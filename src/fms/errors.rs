/*
MIT License with FEFF10 Attribution

Copyright (c) 2025 Ameyanagi

Based on or developed using Distribution: FEFF10.0
Copyright (c) 2020 FEFF Project, University of Washington and SLAC National Accelerator Laboratory.
All rights reserved.
*/

//! Error types for the FMS module

use thiserror::Error;

/// Result type for FMS operations
pub type Result<T> = std::result::Result<T, FmsError>;

/// FMS-specific errors
#[derive(Error, Debug)]
pub enum FmsError {
    /// Error when matrix operations fail
    #[error("Matrix operation error: {0}")]
    MatrixError(String),

    /// Error when the matrix dimensions are incompatible
    #[error("Matrix dimension mismatch: {0}")]
    DimensionMismatch(String),

    /// Error when the matrix solver fails to converge
    #[error("Solver failed to converge: {0}")]
    ConvergenceError(String),

    /// Error when iterative solver method fails
    #[error("Iterative solver failed: {0}")]
    IterationFailed(String),

    /// Error when parameters are invalid
    #[error("Invalid parameter: {0}")]
    InvalidParameter(String),

    /// Error when the FMS radius is too small
    #[error("FMS radius too small: {0}")]
    RadiusTooSmall(String),

    /// Error from the atoms module
    #[error("Atom error: {0}")]
    AtomError(String),

    /// Error from the scattering module
    #[error("Scattering error: {0}")]
    ScatteringError(String),

    /// General calculation error
    #[error("Calculation error: {0}")]
    CalculationError(String),
}

// Implement conversions from other error types
impl From<crate::atoms::AtomError> for FmsError {
    fn from(err: crate::atoms::AtomError) -> Self {
        FmsError::AtomError(err.to_string())
    }
}

impl From<crate::utils::errors::MathError> for FmsError {
    fn from(err: crate::utils::errors::MathError) -> Self {
        FmsError::CalculationError(format!("Math error: {}", err))
    }
}

impl From<crate::potential::PotentialError> for FmsError {
    fn from(err: crate::potential::PotentialError) -> Self {
        FmsError::CalculationError(format!("Potential error: {}", err))
    }
}

impl From<crate::atoms::Result<String>> for FmsError {
    fn from(result: crate::atoms::Result<String>) -> Self {
        match result {
            Ok(msg) => FmsError::CalculationError(msg),
            Err(err) => FmsError::AtomError(err.to_string()),
        }
    }
}
