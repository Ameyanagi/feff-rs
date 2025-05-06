/*
MIT License with FEFF10 Attribution

Copyright (c) 2025 Ameyanagi

Based on or developed using Distribution: FEFF10.0
Copyright (c) 2020 FEFF Project, University of Washington and SLAC National Accelerator Laboratory.
All rights reserved.
*/

//! Path degeneracy calculation
//!
//! This module implements algorithms for calculating path degeneracy and
//! identifying equivalent paths in XAFS calculations.

use std::collections::{HashMap, HashSet};

use crate::atoms::structure::AtomicStructure;
use crate::atoms::vector::Vector3D;
use crate::path::path::{Path, PathType};

/// Checks if two paths are geometrically equivalent
///
/// Two paths are geometrically equivalent if they have the same sequence
/// of leg lengths and scattering angles, even if they involve different atoms.
///
/// # Arguments
///
/// * `path1` - First path to compare
/// * `path2` - Second path to compare
/// * `structure` - Atomic structure containing the atoms
/// * `angle_tolerance` - Tolerance for angle comparison (in radians)
/// * `length_tolerance` - Tolerance for length comparison (in Å)
///
/// # Returns
///
/// `true` if the paths are geometrically equivalent, `false` otherwise
pub fn are_paths_geometrically_equivalent(
    path1: &Path,
    path2: &Path,
    structure: &AtomicStructure,
    angle_tolerance: f64,
    length_tolerance: f64,
) -> bool {
    // Paths must have the same number of legs
    if path1.legs.len() != path2.legs.len() {
        return false;
    }
    
    // Total path lengths must be similar
    if (path1.total_length - path2.total_length).abs() > length_tolerance {
        return false;
    }
    
    // Check individual leg lengths
    for (leg1, leg2) in path1.legs.iter().zip(path2.legs.iter()) {
        if (leg1.length - leg2.length).abs() > length_tolerance {
            return false;
        }
    }
    
    // For paths with 3 or more legs, check scattering angles
    if path1.legs.len() >= 3 {
        // Compare scattering angles at each middle atom
        for i in 1..path1.legs.len() {
            let angle1 = calculate_scattering_angle(path1, i, structure);
            let angle2 = calculate_scattering_angle(path2, i, structure);
            
            if (angle1 - angle2).abs() > angle_tolerance {
                return false;
            }
        }
    }
    
    true
}

/// Calculates the scattering angle at a specific position in a path
///
/// The scattering angle is the angle between the incoming and outgoing
/// legs at a scattering atom. For the absorber, this is the angle between
/// the first and last legs.
///
/// # Arguments
///
/// * `path` - The path containing the scattering atom
/// * `position` - Index in the atom_sequence for the scattering atom
/// * `structure` - Atomic structure containing the atoms
///
/// # Returns
///
/// The scattering angle in radians
fn calculate_scattering_angle(path: &Path, position: usize, structure: &AtomicStructure) -> f64 {
    // Handle boundary cases
    if position == 0 || position >= path.atom_sequence.len() - 1 {
        return 0.0;
    }
    
    let atom_index = path.atom_sequence[position];
    
    // Get the previous and next atoms in the path
    let prev_atom_index = path.atom_sequence[position - 1];
    let next_atom_index = path.atom_sequence[position + 1];
    
    // Get atom positions
    let atom_pos = structure.atom(atom_index).unwrap().position().clone();
    let prev_atom_pos = structure.atom(prev_atom_index).unwrap().position().clone();
    let next_atom_pos = structure.atom(next_atom_index).unwrap().position().clone();
    
    // Calculate vectors for incoming and outgoing legs
    let incoming_vec = (atom_pos - prev_atom_pos).normalize();
    let outgoing_vec = (next_atom_pos - atom_pos).normalize();
    
    // Calculate the angle between the vectors
    let cos_angle = incoming_vec.dot(&outgoing_vec);
    
    // Return the angle in radians, clamped to valid range
    cos_angle.clamp(-1.0, 1.0).acos()
}

/// Groups equivalent paths and calculates degeneracy
///
/// This function identifies equivalent paths in a set of paths and
/// groups them together, calculating the degeneracy of each group.
///
/// # Arguments
///
/// * `paths` - Vector of paths to group
/// * `structure` - Atomic structure containing the atoms
/// * `angle_tolerance` - Tolerance for angle comparison (in radians)
/// * `length_tolerance` - Tolerance for length comparison (in Å)
///
/// # Returns
///
/// A vector of unique paths with degeneracy values set
pub fn calculate_path_degeneracies(
    paths: Vec<Path>,
    structure: &AtomicStructure,
    angle_tolerance: f64,
    length_tolerance: f64,
) -> Vec<Path> {
    // Group paths by length for more efficient comparison
    let mut length_groups: HashMap<i64, Vec<Path>> = HashMap::new();
    
    for path in paths {
        // Round the length to group similar paths
        let rounded_length = (path.total_length * 1000.0).round() as i64;
        length_groups.entry(rounded_length).or_default().push(path);
    }
    
    let mut unique_paths = Vec::new();
    
    // For each length group, find equivalent paths
    for paths in length_groups.values() {
        let mut path_groups: Vec<Vec<Path>> = Vec::new();
        
        'outer: for path in paths {
            // Check if this path is equivalent to any existing group
            for group in &mut path_groups {
                let representative = &group[0];
                
                if are_paths_geometrically_equivalent(
                    representative,
                    path,
                    structure,
                    angle_tolerance,
                    length_tolerance,
                ) {
                    // Path is equivalent to this group, add it
                    group.push(path.clone());
                    continue 'outer;
                }
            }
            
            // Path doesn't match any group, create a new group
            path_groups.push(vec![path.clone()]);
        }
        
        // Create a representative path for each group with degeneracy set
        for group in path_groups {
            let mut representative = group[0].clone();
            representative.set_degeneracy(group.len() as u32);
            unique_paths.push(representative);
        }
    }
    
    unique_paths
}

/// Calculates the effective degeneracy of a path considering symmetry
///
/// This function identifies the degeneracy of a path by considering
/// symmetry operations in the crystal structure.
///
/// # Arguments
///
/// * `path` - The path to analyze
/// * `structure` - Atomic structure containing the atoms
///
/// # Returns
///
/// The calculated degeneracy
pub fn calculate_single_path_degeneracy(path: &Path, structure: &AtomicStructure) -> u32 {
    match path.path_type {
        PathType::SingleScattering => {
            // For single scattering, count atoms of the same type at the same distance
            let scatterer_index = path.atom_sequence[1];
            let scatterer = structure.atom(scatterer_index).unwrap();
            let absorber_index = path.atom_sequence[0];
            let absorber = structure.atom(absorber_index).unwrap();
            
            let mut equivalent_count = 0;
            
            // Find atoms of the same type (atomic number) at approximately the same distance
            for (index, atom) in structure.atoms().iter().enumerate() {
                if index == absorber_index || index == scatterer_index {
                    continue;
                }
                
                if atom.atomic_number() == scatterer.atomic_number() {
                    let distance = atom.distance_to(absorber);
                    
                    // Check if this atom is at approximately the same distance
                    if (distance - path.legs[0].length).abs() < 0.01 {
                        equivalent_count += 1;
                    }
                }
            }
            
            // Add 1 for the original scatterer
            equivalent_count + 1
        },
        PathType::DoubleScattering | PathType::Triangle => {
            // For double scattering and triangle paths, the calculation is more complex
            // This is a simplified version that considers only the atoms in the path
            let mut equivalent_count = 1;
            
            // For each scatterer in the path, find equivalent atoms
            let unique_scatterers: HashSet<usize> = path.atom_sequence
                .iter()
                .skip(1) // Skip absorber
                .copied()
                .collect();
            
            for &scatterer_index in &unique_scatterers {
                let scatterer = structure.atom(scatterer_index).unwrap();
                
                // Count atoms of the same type
                let same_type_count = structure.atoms()
                    .iter()
                    .filter(|atom| atom.atomic_number() == scatterer.atomic_number())
                    .count();
                
                // Adjust the degeneracy
                equivalent_count *= same_type_count as u32;
            }
            
            // This is a simplified calculation; in reality, we would need to consider
            // the exact geometric arrangement of atoms and symmetry operations
            equivalent_count
        },
        PathType::MultipleScattering => {
            // For more complex paths, this would use the full symmetry analysis
            // For now, return 1 as a default
            1
        },
    }
}