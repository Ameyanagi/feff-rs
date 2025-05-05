/*
MIT License with FEFF10 Attribution

Copyright (c) 2025 Ameyanagi

Based on or developed using Distribution: FEFF10.0
Copyright (c) 2020 FEFF Project, University of Washington and SLAC National Accelerator Laboratory.
All rights reserved.
*/

//! Atomic data and calculations module
//!
//! This module provides atomic data and calculations for FEFF.

use num_complex::Complex64;

/// Error types for the atoms module
#[derive(Debug, thiserror::Error)]
pub enum AtomError {
    #[error("Invalid atomic number: {0}")]
    InvalidAtomicNumber(i32),
    
    #[error("Calculation error: {0}")]
    CalculationError(String),
}

/// Result type for atom operations
pub type Result<T> = std::result::Result<T, AtomError>;

/// Represents an atom in the calculation
pub struct Atom {
    // Will be implemented as development progresses
    atomic_number: i32,
    symbol: String,
}

impl Atom {
    /// Create a new atom with the given atomic number
    pub fn new(atomic_number: i32) -> Result<Self> {
        if atomic_number <= 0 || atomic_number > 118 {
            return Err(AtomError::InvalidAtomicNumber(atomic_number));
        }
        
        // Will be expanded with proper symbol lookup
        Ok(Self {
            atomic_number,
            symbol: format!("Element{}", atomic_number),
        })
    }
    
    /// Get the atomic number
    pub fn atomic_number(&self) -> i32 {
        self.atomic_number
    }
    
    /// Get the atomic symbol
    pub fn symbol(&self) -> &str {
        &self.symbol
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_atom_creation() {
        let atom = Atom::new(29).unwrap(); // Copper
        assert_eq!(atom.atomic_number(), 29);
    }
    
    #[test]
    fn test_invalid_atom() {
        assert!(Atom::new(0).is_err());
        assert!(Atom::new(119).is_err());
    }
}
