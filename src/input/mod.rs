/*
MIT License with FEFF10 Attribution

Copyright (c) 2025 Ameyanagi

Based on or developed using Distribution: FEFF10.0
Copyright (c) 2020 FEFF Project, University of Washington and SLAC National Accelerator Laboratory.
All rights reserved.
*/

//! FEFF input parsing and processing
//!
//! This module provides functionality for parsing and handling FEFF input files.
//! It supports all standard FEFF cards and parameters, and follows FEFF10's approach
//! for representing input parameters.

pub mod card;
pub mod config;
pub mod coordinates;
pub mod errors;
pub mod model;
pub mod parameters;
pub mod parser;

// Re-export commonly used types and functions for convenience
pub use card::Card;
pub use config::ParserConfig;
pub use coordinates::CoordinateSystem;
pub use errors::{InputError, Result};
pub use model::FeffInput;
pub use parser::FeffInputParser;

/// Parse a FEFF input file with default configuration
pub fn parse_feff_input<P: AsRef<std::path::Path>>(path: P) -> Result<FeffInput> {
    let config = ParserConfig {
        input_path: path.as_ref().to_path_buf(),
        ..Default::default()
    };

    let mut parser = FeffInputParser::new(config);
    parser.parse::<&std::path::Path>(None)
}

/// Create a new empty FEFF input
pub fn new_feff_input() -> FeffInput {
    FeffInput::new()
}
