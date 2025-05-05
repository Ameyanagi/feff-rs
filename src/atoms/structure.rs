/*
MIT License with FEFF10 Attribution

Copyright (c) 2025 Ameyanagi

Based on or developed using Distribution: FEFF10.0
Copyright (c) 2020 FEFF Project, University of Washington and SLAC National Accelerator Laboratory.
All rights reserved.
*/

//! Atomic structure representation for FEFF calculations

use super::atom::Atom;
use super::errors::{AtomError, Result};
use super::potential::PotentialType;

/// AtomicStructure represents a collection of atoms for FEFF calculations
/// This is equivalent to the cluster in FEFF10
#[derive(Debug, Default, Clone)]
pub struct AtomicStructure {
    /// List of atoms in the structure
    atoms: Vec<Atom>,
    /// List of potential types
    potential_types: Vec<PotentialType>,
    /// Central atom index
    central_atom_index: Option<usize>,
    /// Title or description of the structure
    title: String,
}

impl AtomicStructure {
    /// Create a new empty atomic structure
    pub fn new() -> Self {
        Self {
            atoms: Vec::new(),
            potential_types: Vec::new(),
            central_atom_index: None,
            title: String::new(),
        }
    }

    /// Create a new atomic structure with a title
    pub fn with_title(title: &str) -> Self {
        let mut structure = Self::new();
        structure.title = title.to_string();
        structure
    }

    /// Get the title of the atomic structure
    pub fn title(&self) -> &str {
        &self.title
    }

    /// Set the title of the atomic structure
    pub fn set_title(&mut self, title: &str) {
        self.title = title.to_string();
    }

    /// Add an atom to the structure
    pub fn add_atom(&mut self, atom: Atom) -> usize {
        let atom_index = self.atoms.len();

        // Update potential type's atom count if it exists
        let pot_type = atom.potential_type() as usize;
        if pot_type < self.potential_types.len() {
            let pot = &mut self.potential_types[pot_type];
            pot.increment_atom_count();
        }

        self.atoms.push(atom);
        atom_index
    }

    /// Add a potential type to the structure
    pub fn add_potential_type(&mut self, potential: PotentialType) -> usize {
        let pot_index = self.potential_types.len();
        self.potential_types.push(potential);
        pot_index
    }

    /// Set the central atom index
    pub fn set_central_atom(&mut self, index: usize) -> Result<()> {
        if index >= self.atoms.len() {
            return Err(AtomError::CalculationError(format!(
                "Central atom index {} out of range (max {})",
                index,
                self.atoms.len() - 1
            )));
        }
        self.central_atom_index = Some(index);
        Ok(())
    }

    /// Get the central atom index
    pub fn central_atom_index(&self) -> Option<usize> {
        self.central_atom_index
    }

    /// Get a reference to the central atom
    pub fn central_atom(&self) -> Option<&Atom> {
        match self.central_atom_index {
            Some(idx) => self.atoms.get(idx),
            None => None,
        }
    }

    /// Get a mutable reference to the central atom
    pub fn central_atom_mut(&mut self) -> Option<&mut Atom> {
        match self.central_atom_index {
            Some(idx) => self.atoms.get_mut(idx),
            None => None,
        }
    }

    /// Get the number of atoms
    pub fn atom_count(&self) -> usize {
        self.atoms.len()
    }

    /// Get the number of potential types
    pub fn potential_type_count(&self) -> usize {
        self.potential_types.len()
    }

    /// Get a reference to an atom by index
    pub fn atom(&self, index: usize) -> Option<&Atom> {
        self.atoms.get(index)
    }

    /// Get a mutable reference to an atom by index
    pub fn atom_mut(&mut self, index: usize) -> Option<&mut Atom> {
        self.atoms.get_mut(index)
    }

    /// Get a reference to a potential type by index
    pub fn potential_type(&self, index: usize) -> Option<&PotentialType> {
        self.potential_types.get(index)
    }

    /// Get a mutable reference to a potential type by index
    pub fn potential_type_mut(&mut self, index: usize) -> Option<&mut PotentialType> {
        self.potential_types.get_mut(index)
    }

    /// Get a slice of all atoms
    pub fn atoms(&self) -> &[Atom] {
        &self.atoms
    }

    /// Get a slice of all potential types
    pub fn potential_types(&self) -> &[PotentialType] {
        &self.potential_types
    }

    /// Calculate distances between all atoms in the structure
    pub fn calculate_distances(&self) -> Vec<Vec<f64>> {
        let n = self.atoms.len();
        let mut distances = vec![vec![0.0; n]; n];

        for i in 0..n {
            for j in i..n {
                let dist = self.atoms[i].distance_to(&self.atoms[j]);
                distances[i][j] = dist;
                distances[j][i] = dist;
            }
        }

        distances
    }

    /// Find atoms within a given distance of a reference atom
    pub fn atoms_within_distance(
        &self,
        reference_atom_idx: usize,
        max_distance: f64,
    ) -> Vec<usize> {
        let mut result = Vec::new();

        if reference_atom_idx >= self.atoms.len() {
            return result;
        }

        let reference = &self.atoms[reference_atom_idx];

        for (idx, atom) in self.atoms.iter().enumerate() {
            if idx != reference_atom_idx && reference.distance_to(atom) <= max_distance {
                result.push(idx);
            }
        }

        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::atoms::Vector3D;
    use approx::assert_relative_eq;

    #[test]
    fn test_atomic_structure_creation() {
        let structure = AtomicStructure::new();

        assert_eq!(structure.atom_count(), 0);
        assert_eq!(structure.potential_type_count(), 0);
        assert_eq!(structure.central_atom_index(), None);
        assert_eq!(structure.title(), "");
    }

    #[test]
    fn test_atomic_structure_with_title() {
        let structure = AtomicStructure::with_title("Test Structure");

        assert_eq!(structure.title(), "Test Structure");
    }

    #[test]
    fn test_adding_atoms_and_potentials() {
        let mut structure = AtomicStructure::new();

        // Add potential types
        let pot1 = PotentialType::new(0, 26).unwrap(); // Fe
        let pot2 = PotentialType::new(1, 8).unwrap(); // O

        let pot1_idx = structure.add_potential_type(pot1);
        let pot2_idx = structure.add_potential_type(pot2);

        assert_eq!(pot1_idx, 0);
        assert_eq!(pot2_idx, 1);
        assert_eq!(structure.potential_type_count(), 2);

        // Add atoms
        let atom1 = Atom::new(26, Vector3D::new(0.0, 0.0, 0.0), 0).unwrap();
        let atom2 = Atom::new(8, Vector3D::new(2.0, 0.0, 0.0), 1).unwrap();
        let atom3 = Atom::new(8, Vector3D::new(-2.0, 0.0, 0.0), 1).unwrap();

        let atom1_idx = structure.add_atom(atom1);
        let atom2_idx = structure.add_atom(atom2);
        let atom3_idx = structure.add_atom(atom3);

        assert_eq!(atom1_idx, 0);
        assert_eq!(atom2_idx, 1);
        assert_eq!(atom3_idx, 2);
        assert_eq!(structure.atom_count(), 3);

        // Check atom counts in potential types
        assert_eq!(structure.potential_type(0).unwrap().atom_count(), 1);
        assert_eq!(structure.potential_type(1).unwrap().atom_count(), 2);
    }

    #[test]
    fn test_central_atom() {
        let mut structure = AtomicStructure::new();

        // Add potential type
        let pot = PotentialType::new(0, 26).unwrap();
        structure.add_potential_type(pot);

        // Add atoms
        let atom1 = Atom::new(26, Vector3D::new(0.0, 0.0, 0.0), 0).unwrap();
        let atom2 = Atom::new(26, Vector3D::new(2.0, 0.0, 0.0), 0).unwrap();

        structure.add_atom(atom1);
        structure.add_atom(atom2);

        // Set central atom
        assert!(structure.set_central_atom(0).is_ok());
        assert_eq!(structure.central_atom_index(), Some(0));

        // Invalid central atom index
        assert!(structure.set_central_atom(10).is_err());
    }

    #[test]
    fn test_distances_and_neighbor_finding() {
        let mut structure = AtomicStructure::new();

        // Add potential type
        let pot = PotentialType::new(0, 26).unwrap();
        structure.add_potential_type(pot);

        // Add atoms in a line: A --- B --- C
        let atom_a = Atom::new(26, Vector3D::new(0.0, 0.0, 0.0), 0).unwrap();
        let atom_b = Atom::new(26, Vector3D::new(2.0, 0.0, 0.0), 0).unwrap();
        let atom_c = Atom::new(26, Vector3D::new(4.0, 0.0, 0.0), 0).unwrap();

        structure.add_atom(atom_a);
        structure.add_atom(atom_b);
        structure.add_atom(atom_c);

        // Calculate distances
        let distances = structure.calculate_distances();

        // Check distances
        assert_relative_eq!(distances[0][0], 0.0, epsilon = 1e-6); // A to A
        assert_relative_eq!(distances[0][1], 2.0, epsilon = 1e-6); // A to B
        assert_relative_eq!(distances[0][2], 4.0, epsilon = 1e-6); // A to C
        assert_relative_eq!(distances[1][2], 2.0, epsilon = 1e-6); // B to C

        // Find atoms within distance 3 of atom A
        let neighbors = structure.atoms_within_distance(0, 3.0);
        assert_eq!(neighbors.len(), 1);
        assert_eq!(neighbors[0], 1); // only B is within distance 3

        // Find atoms within distance 5 of atom A
        let neighbors = structure.atoms_within_distance(0, 5.0);
        assert_eq!(neighbors.len(), 2);
        assert!(neighbors.contains(&1)); // B
        assert!(neighbors.contains(&2)); // C
    }
}
