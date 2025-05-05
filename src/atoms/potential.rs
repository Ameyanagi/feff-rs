/*
MIT License with FEFF10 Attribution

Copyright (c) 2025 Ameyanagi

Based on or developed using Distribution: FEFF10.0
Copyright (c) 2020 FEFF Project, University of Washington and SLAC National Accelerator Laboratory.
All rights reserved.
*/

//! Potential types used in FEFF calculations

use super::database;
use super::errors::{AtomError, Result};
use std::fmt;

/// Represents a potential type used in FEFF calculations
/// Multiple atoms can share the same potential type
#[derive(Debug, Clone)]
pub struct PotentialType {
    /// Index of this potential type
    index: i32,
    /// Atomic number (Z) of the element
    atomic_number: i32,
    /// Element symbol
    symbol: String,
    /// Muffin-tin radius
    muffin_tin_radius: f64,
    /// Overlap factor for potential calculation
    overlap_factor: f64,
    /// Number of atoms with this potential type
    atom_count: usize,
    /// Optional spin configuration
    spin: f64,
    /// Model atom index (which atom is the representative for this potential)
    model_atom_index: Option<usize>,
}

impl PotentialType {
    /// Create a new potential type with default values
    pub fn new(index: i32, atomic_number: i32) -> Result<Self> {
        if atomic_number <= 0 || atomic_number > 118 {
            return Err(AtomError::InvalidAtomicNumber(atomic_number));
        }

        if index < 0 {
            return Err(AtomError::InvalidPotentialType(index));
        }

        // Get the element symbol from the database
        let symbol = match database::element_symbol(atomic_number) {
            Some(s) => s.to_string(),
            None => format!("Element{}", atomic_number),
        };

        Ok(Self {
            index,
            atomic_number,
            symbol,
            muffin_tin_radius: 0.0,
            overlap_factor: 1.15, // FEFF default
            atom_count: 0,
            spin: 0.0,
            model_atom_index: None,
        })
    }

    /// Create a new potential type with specified properties
    pub fn with_properties(
        index: i32,
        atomic_number: i32,
        muffin_tin_radius: f64,
        overlap_factor: f64,
        spin: f64,
    ) -> Result<Self> {
        let mut pot = Self::new(index, atomic_number)?;
        pot.muffin_tin_radius = muffin_tin_radius;
        pot.overlap_factor = overlap_factor;
        pot.spin = spin;

        Ok(pot)
    }

    /// Get the index of this potential type
    pub fn index(&self) -> i32 {
        self.index
    }

    /// Get the atomic number
    pub fn atomic_number(&self) -> i32 {
        self.atomic_number
    }

    /// Get the element symbol
    pub fn symbol(&self) -> &str {
        &self.symbol
    }

    /// Get the muffin-tin radius
    pub fn muffin_tin_radius(&self) -> f64 {
        self.muffin_tin_radius
    }

    /// Set the muffin-tin radius
    pub fn set_muffin_tin_radius(&mut self, radius: f64) {
        self.muffin_tin_radius = radius;
    }

    /// Get the overlap factor
    pub fn overlap_factor(&self) -> f64 {
        self.overlap_factor
    }

    /// Set the overlap factor
    pub fn set_overlap_factor(&mut self, factor: f64) {
        self.overlap_factor = factor;
    }

    /// Get the atom count (number of atoms with this potential)
    pub fn atom_count(&self) -> usize {
        self.atom_count
    }

    /// Increment the atom count
    pub fn increment_atom_count(&mut self) {
        self.atom_count += 1;
    }

    /// Set the atom count explicitly
    pub fn set_atom_count(&mut self, count: usize) {
        self.atom_count = count;
    }

    /// Get the spin configuration
    pub fn spin(&self) -> f64 {
        self.spin
    }

    /// Set the spin configuration
    pub fn set_spin(&mut self, spin: f64) {
        self.spin = spin;
    }

    /// Get the model atom index
    pub fn model_atom_index(&self) -> Option<usize> {
        self.model_atom_index
    }

    /// Set the model atom index
    pub fn set_model_atom_index(&mut self, index: usize) {
        self.model_atom_index = Some(index);
    }

    /// Get the atomic weight from the database
    pub fn atomic_weight(&self) -> Option<f64> {
        database::atomic_weight(self.atomic_number)
    }

    /// Get the covalent radius from the database
    pub fn covalent_radius(&self) -> Option<f64> {
        database::covalent_radius(self.atomic_number)
    }
}

impl fmt::Display for PotentialType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "Potential {} - {} (Z={}), r_mt={:.6}, folp={:.6}, count={}",
            self.index,
            self.symbol,
            self.atomic_number,
            self.muffin_tin_radius,
            self.overlap_factor,
            self.atom_count
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_potential_type_creation() {
        let pot = PotentialType::new(1, 26).unwrap(); // Iron

        assert_eq!(pot.index(), 1);
        assert_eq!(pot.atomic_number(), 26);
        assert_eq!(pot.symbol(), "Fe");
        assert_eq!(pot.atom_count(), 0);
        assert_eq!(pot.spin(), 0.0);
        assert_relative_eq!(pot.overlap_factor(), 1.15, epsilon = 1e-6); // Default
        assert_relative_eq!(pot.muffin_tin_radius(), 0.0, epsilon = 1e-6);
        assert_eq!(pot.model_atom_index(), None);
    }

    #[test]
    fn test_potential_type_with_properties() {
        let pot = PotentialType::with_properties(1, 26, 2.0, 1.2, 1.0).unwrap();

        assert_eq!(pot.index(), 1);
        assert_eq!(pot.atomic_number(), 26);
        assert_eq!(pot.symbol(), "Fe");
        assert_eq!(pot.atom_count(), 0);
        assert_eq!(pot.spin(), 1.0);
        assert_relative_eq!(pot.overlap_factor(), 1.2, epsilon = 1e-6);
        assert_relative_eq!(pot.muffin_tin_radius(), 2.0, epsilon = 1e-6);
    }

    #[test]
    fn test_invalid_potential_type() {
        assert!(PotentialType::new(-1, 26).is_err());
        assert!(PotentialType::new(1, 0).is_err());
        assert!(PotentialType::new(1, 119).is_err());
    }

    #[test]
    fn test_potential_type_methods() {
        let mut pot = PotentialType::new(1, 26).unwrap();

        pot.set_muffin_tin_radius(2.5);
        assert_relative_eq!(pot.muffin_tin_radius(), 2.5, epsilon = 1e-6);

        pot.set_overlap_factor(1.3);
        assert_relative_eq!(pot.overlap_factor(), 1.3, epsilon = 1e-6);

        pot.set_spin(0.5);
        assert_relative_eq!(pot.spin(), 0.5, epsilon = 1e-6);

        pot.increment_atom_count();
        assert_eq!(pot.atom_count(), 1);

        pot.set_atom_count(5);
        assert_eq!(pot.atom_count(), 5);

        pot.set_model_atom_index(3);
        assert_eq!(pot.model_atom_index(), Some(3));
    }
}
