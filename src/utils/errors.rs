/*
MIT License with FEFF10 Attribution

Copyright (c) 2025 Ameyanagi

Based on or developed using Distribution: FEFF10.0
Copyright (c) 2020 FEFF Project, University of Washington and SLAC National Accelerator Laboratory.
All rights reserved.
*/

//! Error types for the utils module

use thiserror::Error;

/// Errors that can occur in the utils module
#[derive(Error, Debug)]
pub enum UtilsError {
    /// Generic error with a message
    #[error("Utility error: {0}")]
    Generic(String),

    /// Math-related errors
    #[error("Math error: {0}")]
    Math(String),
}

/// Alias for Math-related errors
pub type MathError = UtilsError;

/// A specialized Result type for utils operations
pub type Result<T> = std::result::Result<T, UtilsError>;
