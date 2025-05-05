/*
MIT License with FEFF10 Attribution

Copyright (c) 2025 Ameyanagi

Based on or developed using Distribution: FEFF10.0
Copyright (c) 2020 FEFF Project, University of Washington and SLAC National Accelerator Laboratory.
All rights reserved.
*/

//! Atom representation for FEFF calculations

use super::database;
use super::errors::{AtomError, Result};
use super::vector::Vector3D;
use std::fmt;

/// Represents an atom in the calculation
#[derive(Debug, Clone)]
pub struct Atom {
    /// Atomic number (Z) of the element
    atomic_number: i32,
    /// Atomic symbol (element symbol)
    symbol: String,
    /// Position of the atom in 3D space
    position: Vector3D,
    /// Potential type index for this atom
    potential_type: i32,
    /// Optional spin configuration (1.0 for spin up, -1.0 for spin down, 0.0 for non-magnetic)
    spin: f64,
}

impl Atom {
    /// Create a new atom with the given atomic number and position
    pub fn new(atomic_number: i32, position: Vector3D, potential_type: i32) -> Result<Self> {
        if atomic_number <= 0 || atomic_number > 118 {
            return Err(AtomError::InvalidAtomicNumber(atomic_number));
        }

        if potential_type < 0 {
            return Err(AtomError::InvalidPotentialType(potential_type));
        }

        // Get the element symbol from the database
        let symbol = match database::element_symbol(atomic_number) {
            Some(s) => s.to_string(),
            None => format!("Element{}", atomic_number),
        };

        Ok(Self {
            atomic_number,
            symbol,
            position,
            potential_type,
            spin: 0.0,
        })
    }

    /// Create a new atom with all properties specified
    pub fn with_properties(
        atomic_number: i32,
        position: Vector3D,
        potential_type: i32,
        spin: f64,
    ) -> Result<Self> {
        let mut atom = Self::new(atomic_number, position, potential_type)?;
        atom.spin = spin;
        Ok(atom)
    }

    /// Get the atomic number
    pub fn atomic_number(&self) -> i32 {
        self.atomic_number
    }

    /// Get the atomic symbol
    pub fn symbol(&self) -> &str {
        &self.symbol
    }

    /// Get the atom's position
    pub fn position(&self) -> &Vector3D {
        &self.position
    }

    /// Set the atom's position
    pub fn set_position(&mut self, position: Vector3D) {
        self.position = position;
    }

    /// Get the potential type index
    pub fn potential_type(&self) -> i32 {
        self.potential_type
    }

    /// Set the potential type index
    pub fn set_potential_type(&mut self, potential_type: i32) -> Result<()> {
        if potential_type < 0 {
            return Err(AtomError::InvalidPotentialType(potential_type));
        }
        self.potential_type = potential_type;
        Ok(())
    }

    /// Get the spin configuration
    pub fn spin(&self) -> f64 {
        self.spin
    }

    /// Set the spin configuration
    pub fn set_spin(&mut self, spin: f64) {
        self.spin = spin;
    }

    /// Get the atomic weight from the database
    pub fn atomic_weight(&self) -> Option<f64> {
        database::atomic_weight(self.atomic_number)
    }

    /// Get the covalent radius from the database
    pub fn covalent_radius(&self) -> Option<f64> {
        database::covalent_radius(self.atomic_number)
    }

    /// Calculate the distance to another atom
    pub fn distance_to(&self, other: &Self) -> f64 {
        self.position.distance(&other.position)
    }
}

impl fmt::Display for Atom {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{} (Z={}) at {} (potential type {})",
            self.symbol, self.atomic_number, self.position, self.potential_type
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_atom_creation() {
        let position = Vector3D::new(1.0, 2.0, 3.0);
        let atom = Atom::new(29, position, 1).unwrap(); // Copper

        assert_eq!(atom.atomic_number(), 29);
        assert_eq!(atom.symbol(), "Cu");
        assert_eq!(atom.position(), &position);
        assert_eq!(atom.potential_type(), 1);
        assert_eq!(atom.spin(), 0.0);
    }

    #[test]
    fn test_atom_with_properties() {
        let position = Vector3D::new(1.0, 2.0, 3.0);
        let atom = Atom::with_properties(29, position, 1, 1.0).unwrap(); // Copper with spin

        assert_eq!(atom.atomic_number(), 29);
        assert_eq!(atom.symbol(), "Cu");
        assert_eq!(atom.position(), &position);
        assert_eq!(atom.potential_type(), 1);
        assert_eq!(atom.spin(), 1.0);
    }

    #[test]
    fn test_invalid_atom() {
        let position = Vector3D::new(0.0, 0.0, 0.0);
        assert!(Atom::new(0, position, 1).is_err());
        assert!(Atom::new(119, position, 1).is_err());
        assert!(Atom::new(29, position, -1).is_err());
    }

    #[test]
    fn test_atom_distance() {
        let atom1 = Atom::new(29, Vector3D::new(0.0, 0.0, 0.0), 1).unwrap();
        let atom2 = Atom::new(8, Vector3D::new(3.0, 4.0, 0.0), 2).unwrap();

        assert_eq!(atom1.distance_to(&atom2), 5.0);
    }
}
