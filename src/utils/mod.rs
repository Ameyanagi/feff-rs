/*
MIT License with FEFF10 Attribution

Copyright (c) 2025 Ameyanagi

Based on or developed using Distribution: FEFF10.0
Copyright (c) 2020 FEFF Project, University of Washington and SLAC National Accelerator Laboratory.
All rights reserved.
*/

//! Utility functions for FEFF calculations
//!
//! This module provides common utilities used throughout the FEFF code.

pub mod constants;
pub mod conversions;
pub mod errors;

// Re-export commonly used types and functions for convenience
pub use constants::*;
pub use conversions::*;
pub use errors::{Result, UtilsError};
