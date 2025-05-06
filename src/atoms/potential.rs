/*
MIT License with FEFF10 Attribution

Copyright (c) 2025 Ameyanagi

Based on or developed using Distribution: FEFF10.0
Copyright (c) 2020 FEFF Project, University of Washington and SLAC National Accelerator Laboratory.
All rights reserved.
*/

//! Potential types used in FEFF calculations

use super::atom::Atom;
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

    /// Set the atomic number
    pub fn set_atomic_number(&mut self, z: i32) -> Result<()> {
        if z <= 0 || z > 118 {
            return Err(AtomError::InvalidAtomicNumber(z));
        }

        self.atomic_number = z;
        self.symbol = crate::atoms::database::element_symbol(z)
            .unwrap_or("X")
            .to_string();

        Ok(())
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

    /// Calculate the default muffin-tin radius based on covalent radius
    /// Returns the radius in angstroms
    pub fn default_muffin_tin_radius(&self) -> f64 {
        match self.covalent_radius() {
            Some(radius) => radius * 1.2, // Scale by 1.2 as a default
            None => 1.0,                  // Default fallback value
        }
    }

    /// Calculate the appropriate muffin-tin radius for this potential type
    /// Based on the distances to neighboring atoms in the structure
    pub fn calculate_muffin_tin_radius(&self, atom: &Atom, atoms: &[Atom]) -> f64 {
        let covalent = self.default_muffin_tin_radius();

        // Find the minimum distance to any atom that's not the reference atom
        let mut min_distance = f64::MAX;
        for other in atoms {
            // Skip the atom itself
            if other.position() == atom.position() {
                continue;
            }

            let distance = atom.distance_to(other);

            // Calculate the adjusted distance based on overlap factor
            // The overlap factor determines how much the muffin-tin spheres are allowed to overlap
            let other_cov_radius = match database::covalent_radius(other.atomic_number()) {
                Some(r) => r * 1.2,
                None => 1.0,
            };

            let effective_distance = distance - other_cov_radius * self.overlap_factor;
            if effective_distance < min_distance {
                min_distance = effective_distance;
            }
        }

        // If no neighboring atoms found, return the default radius
        if min_distance == f64::MAX || min_distance <= 0.0 {
            return covalent;
        }

        // The effective muffin-tin radius is the smaller of:
        // 1. The default radius based on covalent radius
        // 2. The minimum distance to a neighbor, adjusted by the overlap factor
        f64::min(covalent, min_distance * self.overlap_factor)
    }

    /// Calculate muffin-tin radii for all atoms of this potential type
    /// This updates the muffin-tin radius for this potential type
    pub fn calculate_optimal_muffin_tin_radius(
        &mut self,
        atoms_of_this_type: &[&Atom],
        all_atoms: &[Atom],
    ) -> Result<f64> {
        if atoms_of_this_type.is_empty() {
            return Err(AtomError::CalculationError(
                "Cannot calculate muffin-tin radius without atoms".to_string(),
            ));
        }

        // Calculate the radius for each atom of this type
        let mut radii = Vec::with_capacity(atoms_of_this_type.len());
        for atom in atoms_of_this_type {
            radii.push(self.calculate_muffin_tin_radius(atom, all_atoms));
        }

        // Use the average radius as the optimal value
        let avg_radius = radii.iter().sum::<f64>() / radii.len() as f64;
        self.muffin_tin_radius = avg_radius;

        Ok(avg_radius)
    }

    /// Calculate the overlap volume between two muffin-tin spheres
    /// Returns the volume in cubic angstroms
    pub fn overlap_volume(&self, atom1: &Atom, atom2: &Atom, radius1: f64, radius2: f64) -> f64 {
        let d = atom1.distance_to(atom2);

        // If spheres don't overlap, return 0
        if d >= radius1 + radius2 {
            return 0.0;
        }

        // If one sphere is contained within the other, return volume of the smaller sphere
        if d <= (radius1 - radius2).abs() {
            let smaller_radius = f64::min(radius1, radius2);
            return (4.0 / 3.0) * std::f64::consts::PI * smaller_radius.powi(3);
        }

        // Calculate the volume of the lens-shaped intersection
        // Formula: V = (pi/12d) * (r1 + r2 - d)^2 * (d^2 + 2d(r1 + r2) - 3(r1 - r2)^2)
        let term1 = (radius1 + radius2 - d).powi(2);
        let term2 = d.powi(2) + 2.0 * d * (radius1 + radius2) - 3.0 * (radius1 - radius2).powi(2);

        (std::f64::consts::PI / (12.0 * d)) * term1 * term2
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
