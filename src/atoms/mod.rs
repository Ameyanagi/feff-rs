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
//! The implementation follows FEFF10's approach for representing
//! atoms, potential types, and atomic structures used in calculations.

use std::fmt;

pub mod database;

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
}

/// Result type for atom operations
pub type Result<T> = std::result::Result<T, AtomError>;

/// Represents a 3D vector for positions and other spatial quantities
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Vector3D {
    /// X coordinate
    pub x: f64,
    /// Y coordinate
    pub y: f64,
    /// Z coordinate
    pub z: f64,
}

impl Vector3D {
    /// Create a new 3D vector
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    /// Create a new vector at the origin
    pub fn origin() -> Self {
        Self::new(0.0, 0.0, 0.0)
    }

    /// Calculate the distance to another vector
    pub fn distance(&self, other: &Self) -> f64 {
        let dx = self.x - other.x;
        let dy = self.y - other.y;
        let dz = self.z - other.z;
        (dx * dx + dy * dy + dz * dz).sqrt()
    }

    /// Calculate the length (magnitude) of the vector
    pub fn length(&self) -> f64 {
        (self.x * self.x + self.y * self.y + self.z * self.z).sqrt()
    }

    /// Calculate the dot product with another vector
    pub fn dot(&self, other: &Self) -> f64 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    /// Calculate the cross product with another vector
    pub fn cross(&self, other: &Self) -> Self {
        Self {
            x: self.y * other.z - self.z * other.y,
            y: self.z * other.x - self.x * other.z,
            z: self.x * other.y - self.y * other.x,
        }
    }

    /// Normalize the vector to unit length
    pub fn normalize(&self) -> Self {
        let len = self.length();
        if len > 1e-10 {
            Self {
                x: self.x / len,
                y: self.y / len,
                z: self.z / len,
            }
        } else {
            Self::origin()
        }
    }
}

impl fmt::Display for Vector3D {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "({:.6}, {:.6}, {:.6})", self.x, self.y, self.z)
    }
}

/// Represents a coordinate system for atomic positions
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum CoordinateSystem {
    /// Cartesian coordinates (x, y, z)
    Cartesian,
    /// Spherical coordinates (r, θ, φ)
    Spherical,
    /// Cylindrical coordinates (ρ, φ, z)
    Cylindrical,
}

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

/// AtomicStructure represents a collection of atoms for FEFF calculations
/// This is equivalent to the cluster in FEFF10
#[derive(Debug, Default)]
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

/// Coordinate conversion functions
pub mod coordinates {
    use super::{AtomError, CoordinateSystem, Result, Vector3D};
    use std::f64::consts::PI;

    /// Convert from one coordinate system to another
    pub fn convert(
        coords: &[f64; 3],
        from_system: CoordinateSystem,
        to_system: CoordinateSystem,
    ) -> Result<[f64; 3]> {
        // First convert to Cartesian
        let cartesian = match from_system {
            CoordinateSystem::Cartesian => [coords[0], coords[1], coords[2]],
            CoordinateSystem::Spherical => {
                let r = coords[0];
                let theta = coords[1] * PI / 180.0; // Convert to radians
                let phi = coords[2] * PI / 180.0;
                [
                    r * theta.sin() * phi.cos(),
                    r * theta.sin() * phi.sin(),
                    r * theta.cos(),
                ]
            }
            CoordinateSystem::Cylindrical => {
                let rho = coords[0];
                let phi = coords[1] * PI / 180.0; // Convert to radians
                let z = coords[2];
                [rho * phi.cos(), rho * phi.sin(), z]
            }
        };

        // Then convert from Cartesian to the target system
        match to_system {
            CoordinateSystem::Cartesian => Ok(cartesian),
            CoordinateSystem::Spherical => {
                let x = cartesian[0];
                let y = cartesian[1];
                let z = cartesian[2];
                let r = (x * x + y * y + z * z).sqrt();

                if r < 1e-10 {
                    return Err(AtomError::CoordinateError(
                        "Cannot convert to spherical coordinates at origin".to_string(),
                    ));
                }

                let theta = (z / r).acos() * 180.0 / PI; // To degrees
                let phi = y.atan2(x) * 180.0 / PI; // To degrees
                Ok([r, theta, phi])
            }
            CoordinateSystem::Cylindrical => {
                let x = cartesian[0];
                let y = cartesian[1];
                let z = cartesian[2];
                let rho = (x * x + y * y).sqrt();
                let phi = y.atan2(x) * 180.0 / PI; // To degrees
                Ok([rho, phi, z])
            }
        }
    }

    /// Create a Vector3D from coordinates in a specific system
    pub fn vector_from_coords(coords: &[f64; 3], system: CoordinateSystem) -> Result<Vector3D> {
        let cartesian = convert(coords, system, CoordinateSystem::Cartesian)?;
        Ok(Vector3D::new(cartesian[0], cartesian[1], cartesian[2]))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

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
    fn test_vector_operations() {
        let v1 = Vector3D::new(1.0, 2.0, 3.0);
        let v2 = Vector3D::new(4.0, 5.0, 6.0);

        // Test distance
        assert_relative_eq!(v1.distance(&v2), 5.196152, epsilon = 1e-6);

        // Test length
        assert_relative_eq!(v1.length(), 3.741657, epsilon = 1e-6);

        // Test dot product
        assert_relative_eq!(v1.dot(&v2), 32.0, epsilon = 1e-6);

        // Test cross product
        let cross = v1.cross(&v2);
        assert_relative_eq!(cross.x, -3.0, epsilon = 1e-6);
        assert_relative_eq!(cross.y, 6.0, epsilon = 1e-6);
        assert_relative_eq!(cross.z, -3.0, epsilon = 1e-6);

        // Test normalize
        let norm = v1.normalize();
        assert_relative_eq!(norm.length(), 1.0, epsilon = 1e-6);
    }

    #[test]
    fn test_coordinate_conversions() {
        // Cartesian to spherical
        let cart = [1.0, 1.0, 1.0];
        let spherical = coordinates::convert(
            &cart,
            CoordinateSystem::Cartesian,
            CoordinateSystem::Spherical,
        )
        .unwrap();

        assert_relative_eq!(spherical[0], 1.732051, epsilon = 1e-6); // r
        assert_relative_eq!(spherical[1], 54.735610, epsilon = 1e-6); // theta (degrees)
        assert_relative_eq!(spherical[2], 45.0, epsilon = 1e-6); // phi (degrees)

        // Spherical back to Cartesian
        let cart_again = coordinates::convert(
            &spherical,
            CoordinateSystem::Spherical,
            CoordinateSystem::Cartesian,
        )
        .unwrap();

        assert_relative_eq!(cart_again[0], cart[0], epsilon = 1e-6);
        assert_relative_eq!(cart_again[1], cart[1], epsilon = 1e-6);
        assert_relative_eq!(cart_again[2], cart[2], epsilon = 1e-6);
    }

    #[test]
    fn test_atom_distance() {
        let atom1 = Atom::new(29, Vector3D::new(0.0, 0.0, 0.0), 1).unwrap();
        let atom2 = Atom::new(8, Vector3D::new(3.0, 4.0, 0.0), 2).unwrap();

        assert_relative_eq!(atom1.distance_to(&atom2), 5.0, epsilon = 1e-6);
    }

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
