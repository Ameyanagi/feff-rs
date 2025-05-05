/*
MIT License with FEFF10 Attribution

Copyright (c) 2025 Ameyanagi

Based on or developed using Distribution: FEFF10.0
Copyright (c) 2020 FEFF Project, University of Washington and SLAC National Accelerator Laboratory.
All rights reserved.
*/

//! Configuration for FEFF input parser

use std::path::PathBuf;

/// FEFF input parser configuration
#[derive(Debug, Clone)]
pub struct ParserConfig {
    /// Path to the input file
    pub input_path: PathBuf,
    /// Whether to validate atomic structure against FEFF requirements
    pub validate: bool,
    /// Whether to add hydrogen to under-coordinated atoms
    pub add_hydrogens: bool,
    /// Whether to enable debugging output during parsing
    pub debug: bool,
}

impl Default for ParserConfig {
    fn default() -> Self {
        Self {
            input_path: PathBuf::from("feff.inp"),
            validate: true,
            add_hydrogens: false,
            debug: false,
        }
    }
}
