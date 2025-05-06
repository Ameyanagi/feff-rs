/*
MIT License with FEFF10 Attribution

Copyright (c) 2025 Ameyanagi

Based on or developed using Distribution: FEFF10.0
Copyright (c) 2020 FEFF Project, University of Washington and SLAC National Accelerator Laboratory.
All rights reserved.
*/

//! Error types for potential calculations

use thiserror::Error;

/// Result type for potential calculations
pub type Result<T> = std::result::Result<T, PotentialError>;

/// Error type for potential-related operations
#[derive(Error, Debug)]
pub enum PotentialError {
    /// General error with message
    #[error("{0}")]
    Generic(String),

    /// Invalid atomic structure
    #[error("Invalid atomic structure: {0}")]
    InvalidStructure(String),

    /// Failed to calculate potential
    #[error("Potential calculation failed: {0}")]
    CalculationError(String),

    /// Invalid exchange-correlation functional
    #[error("Invalid exchange-correlation functional: {0}")]
    InvalidExchangeCorrelation(String),

    /// Self-consistency iteration failed to converge
    #[error("Failed to achieve self-consistency: {0}")]
    ConvergenceError(String),

    /// Error in solving Schrödinger equation
    #[error("Failed to solve Schrödinger equation: {0}")]
    SchrodingerError(String),

    /// Propagation of error from atoms module
    #[error("Atom error: {0}")]
    AtomError(#[from] crate::atoms::errors::AtomError),

    /// Propagation of error from utils module
    #[error("Utils error: {0}")]
    UtilsError(#[from] crate::utils::errors::UtilsError),
}
