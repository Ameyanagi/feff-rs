/*
MIT License with FEFF10 Attribution

Copyright (c) 2025 Ameyanagi

Based on or developed using Distribution: FEFF10.0
Copyright (c) 2020 FEFF Project, University of Washington and SLAC National Accelerator Laboratory.
All rights reserved.
*/

//! Error types for the atoms module

/// Error types for the atoms module
#[derive(Debug, thiserror::Error)]
pub enum AtomError {
    #[error("Invalid atomic number: {0}")]
    InvalidAtomicNumber(i32),

    #[error("Invalid potential type: {0}")]
    InvalidPotentialType(i32),

    #[error("Calculation error: {0}")]
    CalculationError(String),

    #[error("Coordinate conversion error: {0}")]
    CoordinateError(String),

    #[error("File error: {0}")]
    FileError(#[from] std::io::Error),

    #[error("Parse error: {0}")]
    ParseError(String),

    #[error("Invalid structure: {0}")]
    InvalidStructure(String),

    #[error("Invalid edge type: {0}")]
    InvalidEdge(String),

    #[error("Potential error: {0}")]
    PotentialError(String),
}

/// Result type for atom operations
pub type Result<T> = std::result::Result<T, AtomError>;
