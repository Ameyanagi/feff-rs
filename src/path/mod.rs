/*
MIT License with FEFF10 Attribution

Copyright (c) 2025 Ameyanagi

Based on or developed using Distribution: FEFF10.0
Copyright (c) 2020 FEFF Project, University of Washington and SLAC National Accelerator Laboratory.
All rights reserved.
*/

//! Path finding and filtering module
//!
//! This module handles enumeration and selection of scattering paths for FEFF.
//! It provides the necessary functionality to identify and filter important
//! scattering paths, which are essential for EXAFS calculations.
//!
//! The module implements path finding algorithms based on the original FEFF10
//! approach, which uses a heap-based search to efficiently explore the path space.
//! It also includes functionality for calculating path degeneracy and filtering
//! paths based on their importance.

// Re-export public components
mod path;
mod finder;
mod degeneracy;
mod filter;
mod amplitude;
mod visualization;

pub use path::{Path, PathLeg, PathType};
pub use finder::{PathFinder, PathFinderConfig};
pub use degeneracy::{calculate_path_degeneracies, calculate_single_path_degeneracy};
pub use filter::{PathFilterConfig, filter_paths, cluster_paths_by_length, select_optimal_path_set};
pub use amplitude::{ExafsParameters, calculate_path_exafs, calculate_exafs_spectrum};
pub use visualization::{format_path_description, create_path_summary_table, generate_path_xyz, 
                        generate_path_json, export_paths};

#[cfg(test)]
mod tests {
    use super::*;
    use crate::atoms::structure::AtomicStructure;
    use crate::atoms::atom::Atom;
    use crate::atoms::vector::Vector3D;

    /// Creates a simple test structure for path finding tests
    fn create_test_structure() -> AtomicStructure {
        let mut structure = AtomicStructure::new();
        
        // Add a central atom (Fe)
        structure.add_atom(Atom::new(26, Vector3D::new(0.0, 0.0, 0.0), 0).unwrap());
        
        // Add 6 O atoms in octahedral arrangement
        structure.add_atom(Atom::new(8, Vector3D::new(2.0, 0.0, 0.0), 1).unwrap());
        structure.add_atom(Atom::new(8, Vector3D::new(-2.0, 0.0, 0.0), 1).unwrap());
        structure.add_atom(Atom::new(8, Vector3D::new(0.0, 2.0, 0.0), 1).unwrap());
        structure.add_atom(Atom::new(8, Vector3D::new(0.0, -2.0, 0.0), 1).unwrap());
        structure.add_atom(Atom::new(8, Vector3D::new(0.0, 0.0, 2.0), 1).unwrap());
        structure.add_atom(Atom::new(8, Vector3D::new(0.0, 0.0, -2.0), 1).unwrap());
        
        // Add 8 Fe atoms at corners of a cube
        structure.add_atom(Atom::new(26, Vector3D::new(4.0, 4.0, 4.0), 0).unwrap());
        structure.add_atom(Atom::new(26, Vector3D::new(4.0, 4.0, -4.0), 0).unwrap());
        structure.add_atom(Atom::new(26, Vector3D::new(4.0, -4.0, 4.0), 0).unwrap());
        structure.add_atom(Atom::new(26, Vector3D::new(4.0, -4.0, -4.0), 0).unwrap());
        structure.add_atom(Atom::new(26, Vector3D::new(-4.0, 4.0, 4.0), 0).unwrap());
        structure.add_atom(Atom::new(26, Vector3D::new(-4.0, 4.0, -4.0), 0).unwrap());
        structure.add_atom(Atom::new(26, Vector3D::new(-4.0, -4.0, 4.0), 0).unwrap());
        structure.add_atom(Atom::new(26, Vector3D::new(-4.0, -4.0, -4.0), 0).unwrap());
        
        structure
    }

    #[test]
    fn test_path_creation() {
        let structure = create_test_structure();
        
        // Create a single scattering path
        let path = Path::create_single_scattering_path(0, 1, &structure);
        
        assert_eq!(path.legs.len(), 2);
        assert_eq!(path.path_type, PathType::SingleScattering);
        assert_eq!(path.atom_sequence, vec![0, 1, 0]);
        assert!((path.total_length - 4.0).abs() < 1e-10); // 2.0 out + 2.0 back = 4.0
    }

    #[test]
    fn test_path_finder() {
        let structure = create_test_structure();
        let absorber_index = 0;
        
        let config = PathFinderConfig {
            max_path_length: 10.0,
            max_paths: 20,
            max_legs: 4,
            importance_threshold: 0.0, // Accept all paths for testing
            cluster_paths: false,
            unique_scatterers_only: true,
        };
        
        let mut finder = PathFinder::new(structure, absorber_index, config);
        let paths = finder.find_paths();
        
        // We should find at least single scattering paths to the 6 O atoms
        assert!(paths.len() >= 6);
        
        // Check that we have some single scattering paths
        let single_scattering_count = paths.iter()
            .filter(|p| p.path_type == PathType::SingleScattering)
            .count();
        
        assert!(single_scattering_count > 0);
    }
    
    #[test]
    fn test_path_degeneracy() {
        let structure = create_test_structure();
        
        // Create single scattering paths to the 6 O atoms
        let paths: Vec<Path> = (1..7)
            .map(|i| Path::create_single_scattering_path(0, i, &structure))
            .collect();
        
        // Calculate degeneracies
        let degenerate_paths = calculate_path_degeneracies(
            paths,
            &structure,
            0.01, // angle tolerance in radians
            0.01, // length tolerance in Å
        );
        
        // We should have 1 unique path with degeneracy 6 (all O atoms are equidistant)
        assert_eq!(degenerate_paths.len(), 1);
        assert_eq!(degenerate_paths[0].degeneracy, 6);
    }
}
