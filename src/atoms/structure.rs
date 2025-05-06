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
use super::vector::Vector3D;

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

    /// Replace an existing potential type with a new one
    pub fn replace_potential_type(&mut self, index: usize, pot: PotentialType) -> Result<()> {
        if index >= self.potential_types.len() {
            return Err(AtomError::InvalidPotentialType(index as i32));
        }

        // Ensure the new potential has the same index as the one being replaced
        let current_index = self.potential_types[index].index();
        if pot.index() != current_index {
            return Err(AtomError::InvalidPotentialType(pot.index()));
        }

        self.potential_types[index] = pot;
        Ok(())
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

    /// Calculate the muffin-tin radii for all potential types in the structure
    /// This updates the muffin-tin radius for each potential type
    pub fn calculate_muffin_tin_radii(&mut self) -> Result<Vec<f64>> {
        if self.atoms.is_empty() {
            return Err(AtomError::CalculationError(
                "Cannot calculate muffin-tin radii without atoms".to_string(),
            ));
        }

        let mut radii = Vec::with_capacity(self.potential_types.len());

        // Process each potential type
        for pot_idx in 0..self.potential_types.len() {
            // Create a filtered list of atoms for this potential type
            let atoms_of_this_type: Vec<&Atom> = self
                .atoms
                .iter()
                .filter(|a| a.potential_type() as usize == pot_idx)
                .collect();

            if atoms_of_this_type.is_empty() {
                // If no atoms of this type, use default radius
                let default_radius = self.potential_types[pot_idx].default_muffin_tin_radius();
                self.potential_types[pot_idx].set_muffin_tin_radius(default_radius);
                radii.push(default_radius);
                continue;
            }

            // Calculate optimal radius for this potential type using all atoms
            let radius = self.potential_types[pot_idx]
                .calculate_optimal_muffin_tin_radius(&atoms_of_this_type, &self.atoms)?;

            radii.push(radius);
        }

        Ok(radii)
    }

    /// Create a cluster centered on the specified atom with a given radius
    pub fn create_cluster(&self, center_idx: usize, radius: f64) -> Result<Self> {
        if center_idx >= self.atoms.len() {
            return Err(AtomError::CalculationError(format!(
                "Center atom index {} out of range (max {})",
                center_idx,
                self.atoms.len() - 1
            )));
        }

        let center = *self.atoms[center_idx].position();
        let mut cluster = Self::with_title(&format!("Cluster centered on atom {}", center_idx));

        // Copy the potential types
        for pot in &self.potential_types {
            cluster.add_potential_type(pot.clone());
        }

        // Add atoms within the radius
        let mut center_added = false;
        for (idx, atom) in self.atoms.iter().enumerate() {
            let distance = Vector3D::new(
                atom.position().x - center.x,
                atom.position().y - center.y,
                atom.position().z - center.z,
            )
            .length();

            if distance <= radius {
                let atom_idx = cluster.add_atom(atom.clone());

                // Set the center atom
                if idx == center_idx {
                    cluster.set_central_atom(atom_idx)?;
                    center_added = true;
                }
            }
        }

        if !center_added {
            return Err(AtomError::CalculationError(
                "Center atom was not included in the cluster".to_string(),
            ));
        }

        Ok(cluster)
    }

    /// Calculate the total overlap volume for all muffin-tin spheres
    pub fn calculate_total_overlap_volume(&self) -> Result<f64> {
        if self.potential_types.is_empty() || self.atoms.is_empty() {
            return Err(AtomError::CalculationError(
                "Cannot calculate overlap volume without atoms and potentials".to_string(),
            ));
        }

        let mut total_overlap = 0.0;

        // For each pair of atoms
        for i in 0..self.atoms.len() {
            for j in (i + 1)..self.atoms.len() {
                let atom1 = &self.atoms[i];
                let atom2 = &self.atoms[j];

                let pot_type1 = self
                    .potential_types
                    .get(atom1.potential_type() as usize)
                    .ok_or_else(|| {
                        AtomError::CalculationError(format!(
                            "Invalid potential type index: {}",
                            atom1.potential_type()
                        ))
                    })?;

                let pot_type2 = self
                    .potential_types
                    .get(atom2.potential_type() as usize)
                    .ok_or_else(|| {
                        AtomError::CalculationError(format!(
                            "Invalid potential type index: {}",
                            atom2.potential_type()
                        ))
                    })?;

                let radius1 = pot_type1.muffin_tin_radius();
                let radius2 = pot_type2.muffin_tin_radius();

                // Use the first potential type for calculations
                // The overlap_volume is symmetric, so it doesn't matter which one we use
                let overlap = pot_type1.overlap_volume(atom1, atom2, radius1, radius2);

                total_overlap += overlap;
            }
        }

        Ok(total_overlap)
    }

    /// Optimize the muffin-tin radii to minimize overlap while maintaining coverage
    pub fn optimize_muffin_tin_radii(
        &mut self,
        max_iterations: usize,
        tolerance: f64,
    ) -> Result<f64> {
        // Calculate initial muffin-tin radii
        self.calculate_muffin_tin_radii()?;

        // Get initial overlap volume
        let mut prev_overlap = self.calculate_total_overlap_volume()?;
        let mut overlap_change = f64::MAX;

        // Iteratively optimize
        let mut iteration = 0;
        while iteration < max_iterations && overlap_change > tolerance {
            // Adjust each potential type's muffin-tin radius
            for pot_idx in 0..self.potential_types.len() {
                let current_radius = self.potential_types[pot_idx].muffin_tin_radius();

                // Try reducing the radius to see if it improves overlap
                let reduced_radius = current_radius * 0.95;
                self.potential_types[pot_idx].set_muffin_tin_radius(reduced_radius);

                // Calculate new overlap
                let new_overlap = self.calculate_total_overlap_volume()?;

                // If reducing made things worse, revert and try increasing
                if new_overlap > prev_overlap {
                    self.potential_types[pot_idx].set_muffin_tin_radius(current_radius);

                    let increased_radius = current_radius * 1.05;
                    self.potential_types[pot_idx].set_muffin_tin_radius(increased_radius);

                    let new_overlap = self.calculate_total_overlap_volume()?;

                    // If increasing also made things worse, revert back
                    if new_overlap > prev_overlap {
                        self.potential_types[pot_idx].set_muffin_tin_radius(current_radius);
                    } else {
                        prev_overlap = new_overlap;
                    }
                } else {
                    prev_overlap = new_overlap;
                }
            }

            // Calculate the change in overlap volume
            let current_overlap = self.calculate_total_overlap_volume()?;
            overlap_change = (prev_overlap - current_overlap).abs();
            prev_overlap = current_overlap;

            iteration += 1;
        }

        Ok(prev_overlap)
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
