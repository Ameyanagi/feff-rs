/*
MIT License with FEFF10 Attribution

Copyright (c) 2025 Ameyanagi

Based on or developed using Distribution: FEFF10.0
Copyright (c) 2020 FEFF Project, University of Washington and SLAC National Accelerator Laboratory.
All rights reserved.
*/

//! Path data structures for FEFF calculations
//!
//! This module defines the data structures used to represent scattering paths
//! for XAFS calculations. These paths describe the journey of a photoelectron
//! from the absorbing atom, through one or more scattering atoms, and back to
//! the absorbing atom.

use std::fmt;
use std::hash::{Hash, Hasher};

use crate::atoms::structure::AtomicStructure;
use crate::atoms::vector::Vector3D;

/// Represents a single leg of a scattering path
///
/// A leg connects two atoms in a scattering path. It contains
/// information about the source and destination atoms, as well as
/// the distance between them and other physical parameters.
#[derive(Debug, Clone, PartialEq)]
pub struct PathLeg {
    /// Index of the source atom in the atomic structure
    pub from_atom: usize,

    /// Index of the destination atom in the atomic structure
    pub to_atom: usize,

    /// The length of the leg (distance between atoms in Å)
    pub length: f64,
}

impl PathLeg {
    /// Creates a new path leg between two atoms
    ///
    /// # Arguments
    ///
    /// * `from_atom` - Index of the source atom
    /// * `to_atom` - Index of the destination atom
    /// * `structure` - Atomic structure containing the atoms
    ///
    /// # Returns
    ///
    /// A new PathLeg with the calculated distance between atoms
    pub fn new(from_atom: usize, to_atom: usize, structure: &AtomicStructure) -> Self {
        let from_pos = structure.atom(from_atom).unwrap().position().clone();
        let to_pos = structure.atom(to_atom).unwrap().position().clone();
        let length = (to_pos - from_pos).length();

        Self {
            from_atom,
            to_atom,
            length,
        }
    }

    /// Returns the vector from the source atom to the destination atom
    ///
    /// # Arguments
    ///
    /// * `structure` - Atomic structure containing the atoms
    ///
    /// # Returns
    ///
    /// A Vector3D representing the direction and distance from source to destination
    pub fn direction_vector(&self, structure: &AtomicStructure) -> Vector3D {
        let from_pos = structure.atom(self.from_atom).unwrap().position().clone();
        let to_pos = structure.atom(self.to_atom).unwrap().position().clone();
        to_pos - from_pos
    }
}

/// Defines the type of a scattering path
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum PathType {
    /// Single scattering path (2 legs, 1 scatterer)
    SingleScattering,

    /// Double scattering path (3 legs, 2 scatterers)
    DoubleScattering,

    /// Triangle path (3 legs, 2 scatterers in a triangular arrangement)
    Triangle,

    /// Multiple scattering path (more than 3 legs)
    MultipleScattering,
}

/// A complete scattering path for XAFS calculations
///
/// A path represents the journey of a photoelectron from the absorbing atom,
/// through one or more scattering atoms, and back to the absorbing atom.
/// Paths form the basis for calculating the EXAFS and XANES spectra.
#[derive(Clone)]
pub struct Path {
    /// The legs that make up this path
    pub legs: Vec<PathLeg>,

    /// The total path length (sum of all leg lengths)
    pub total_length: f64,

    /// The type of the path (single scattering, double scattering, etc.)
    pub path_type: PathType,

    /// The degeneracy of this path (number of equivalent paths)
    pub degeneracy: u32,

    /// The importance factor for this path (used for filtering)
    pub importance: f64,

    /// Sequence of atom indices visited in this path, including the central atom
    pub atom_sequence: Vec<usize>,
}

impl Path {
    /// Creates a new path from a sequence of legs
    ///
    /// # Arguments
    ///
    /// * `legs` - The legs that make up the path
    /// * `absorber_index` - Index of the absorbing atom
    ///
    /// # Returns
    ///
    /// A new Path with calculated properties
    pub fn new(legs: Vec<PathLeg>, absorber_index: usize) -> Self {
        let total_length = legs.iter().map(|leg| leg.length).sum();

        // Determine path type based on number of legs and arrangement
        let path_type = match legs.len() {
            2 => PathType::SingleScattering,
            3 => {
                if legs[0].from_atom == absorber_index && legs[2].to_atom == absorber_index {
                    PathType::DoubleScattering
                } else {
                    PathType::Triangle
                }
            }
            _ => PathType::MultipleScattering,
        };

        // Build atom sequence including the central atom
        let mut atom_sequence = Vec::with_capacity(legs.len() + 1);
        atom_sequence.push(absorber_index);

        for leg in &legs {
            atom_sequence.push(leg.to_atom);
        }

        Self {
            legs,
            total_length,
            path_type,
            degeneracy: 1,
            importance: 0.0,
            atom_sequence,
        }
    }

    /// Creates a single scattering path from the absorber to a scatterer and back
    ///
    /// # Arguments
    ///
    /// * `absorber_index` - Index of the absorbing atom
    /// * `scatterer_index` - Index of the scattering atom
    /// * `structure` - Atomic structure containing the atoms
    ///
    /// # Returns
    ///
    /// A single scattering Path
    pub fn create_single_scattering_path(
        absorber_index: usize,
        scatterer_index: usize,
        structure: &AtomicStructure,
    ) -> Self {
        let outgoing_leg = PathLeg::new(absorber_index, scatterer_index, structure);
        let returning_leg = PathLeg::new(scatterer_index, absorber_index, structure);

        let legs = vec![outgoing_leg, returning_leg];
        Self::new(legs, absorber_index)
    }

    /// Gets the effective path length (half of the total path length)
    ///
    /// This is used in EXAFS calculations where the effective path length
    /// is used to calculate the phase shift.
    ///
    /// # Returns
    ///
    /// The effective path length in Å
    pub fn effective_length(&self) -> f64 {
        self.total_length / 2.0
    }

    /// Gets the maximum leg length in the path
    ///
    /// This can be used as a filtering criterion, as paths with very long
    /// legs have smaller contributions to the EXAFS signal.
    ///
    /// # Returns
    ///
    /// The maximum leg length in Å
    pub fn max_leg_length(&self) -> f64 {
        self.legs.iter().map(|leg| leg.length).fold(0.0, f64::max)
    }

    /// Gets the number of legs in the path
    ///
    /// # Returns
    ///
    /// The number of legs
    pub fn leg_count(&self) -> usize {
        self.legs.len()
    }

    /// Gets the number of scatterers in the path (excluding the absorber)
    ///
    /// For a valid path that starts and ends at the absorber, this is
    /// the number of unique atoms visited excluding the absorber.
    ///
    /// # Returns
    ///
    /// The number of scatterers
    pub fn scatterer_count(&self) -> usize {
        let mut unique_atoms = std::collections::HashSet::new();

        for leg in &self.legs {
            if leg.from_atom != self.atom_sequence[0] {
                unique_atoms.insert(leg.from_atom);
            }
            if leg.to_atom != self.atom_sequence[0] {
                unique_atoms.insert(leg.to_atom);
            }
        }

        unique_atoms.len()
    }

    /// Sets the degeneracy of this path
    ///
    /// # Arguments
    ///
    /// * `degeneracy` - The degeneracy value to set
    pub fn set_degeneracy(&mut self, degeneracy: u32) {
        self.degeneracy = degeneracy;
    }

    /// Sets the importance factor for this path
    ///
    /// # Arguments
    ///
    /// * `importance` - The importance value to set
    pub fn set_importance(&mut self, importance: f64) {
        self.importance = importance;
    }
}

impl fmt::Debug for Path {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Path {{ type: {:?}, length: {:.3} Å, atoms: {:?}, degeneracy: {}, importance: {:.6} }}",
               self.path_type, self.total_length, self.atom_sequence, self.degeneracy, self.importance)
    }
}

impl PartialEq for Path {
    fn eq(&self, other: &Self) -> bool {
        // Two paths are equal if they have the same atom sequence and similar total length
        const LENGTH_TOLERANCE: f64 = 1e-6;

        if self.atom_sequence.len() != other.atom_sequence.len() {
            return false;
        }

        if (self.total_length - other.total_length).abs() > LENGTH_TOLERANCE {
            return false;
        }

        // Compare atom sequences
        self.atom_sequence == other.atom_sequence
    }
}

impl Eq for Path {}

impl Hash for Path {
    fn hash<H: Hasher>(&self, state: &mut H) {
        // Hash based on atom sequence and rounded path length
        // (for robust comparison when looking for path degeneracy)
        let rounded_length = (self.total_length * 1000.0).round() / 1000.0;

        self.atom_sequence.hash(state);
        ((rounded_length * 1000.0) as i64).hash(state);
    }
}
