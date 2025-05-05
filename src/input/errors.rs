/*
MIT License with FEFF10 Attribution

Copyright (c) 2025 Ameyanagi

Based on or developed using Distribution: FEFF10.0
Copyright (c) 2020 FEFF Project, University of Washington and SLAC National Accelerator Laboratory.
All rights reserved.
*/

//! Error types for FEFF input parsing

use std::io;
use thiserror::Error;

/// Errors that can occur during FEFF input parsing
#[derive(Error, Debug)]
pub enum InputError {
    #[error("IO error: {0}")]
    IoError(#[from] io::Error),

    #[error("Parse error: {0}")]
    ParseError(String),

    #[error("Invalid input format: {0}")]
    InvalidFormat(String),

    #[error("Missing required card: {0}")]
    MissingCard(String),

    #[error("Invalid potential: {0}")]
    InvalidPotential(String),

    #[error("Invalid atomic structure: {0}")]
    InvalidStructure(String),
}

/// Result type for input operations
pub type Result<T> = std::result::Result<T, InputError>;
