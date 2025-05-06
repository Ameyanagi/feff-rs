/*
MIT License with FEFF10 Attribution

Copyright (c) 2025 Ameyanagi

Based on or developed using Distribution: FEFF10.0
Copyright (c) 2020 FEFF Project, University of Washington and SLAC National Accelerator Laboratory.
All rights reserved.
*/

//! Tests for path finding and filtering
//!
//! These tests verify the path finding, filtering, and degeneracy
//! calculation functionality in the path module.

use std::collections::HashSet;

use feff_rs::atoms::atom::Atom;
use feff_rs::atoms::structure::AtomicStructure;
use feff_rs::atoms::vector::Vector3D;
use feff_rs::path::{
    calculate_path_degeneracies, filter_paths, Path, PathFilterConfig, PathFinder,
    PathFinderConfig, PathLeg, PathType,
};

/// Creates a test structure for path finding tests
fn create_test_structure() -> AtomicStructure {
    let mut structure = AtomicStructure::new();

    // Add a central Fe atom at the origin
    structure.add_atom(Atom::new(26, Vector3D::new(0.0, 0.0, 0.0), 0).unwrap());

    // Add O atoms in an octahedral arrangement (along x, y, z axes)
    structure.add_atom(Atom::new(8, Vector3D::new(2.0, 0.0, 0.0), 1).unwrap());
    structure.add_atom(Atom::new(8, Vector3D::new(-2.0, 0.0, 0.0), 1).unwrap());
    structure.add_atom(Atom::new(8, Vector3D::new(0.0, 2.0, 0.0), 1).unwrap());
    structure.add_atom(Atom::new(8, Vector3D::new(0.0, -2.0, 0.0), 1).unwrap());
    structure.add_atom(Atom::new(8, Vector3D::new(0.0, 0.0, 2.0), 1).unwrap());
    structure.add_atom(Atom::new(8, Vector3D::new(0.0, 0.0, -2.0), 1).unwrap());

    // Add Fe atoms for second shell scatterers
    structure.add_atom(Atom::new(26, Vector3D::new(4.0, 0.0, 0.0), 0).unwrap());
    structure.add_atom(Atom::new(26, Vector3D::new(0.0, 4.0, 0.0), 0).unwrap());
    structure.add_atom(Atom::new(26, Vector3D::new(0.0, 0.0, 4.0), 0).unwrap());

    // Add atoms for multiple scattering paths
    structure.add_atom(Atom::new(26, Vector3D::new(2.0, 2.0, 0.0), 0).unwrap());
    structure.add_atom(Atom::new(8, Vector3D::new(2.0, 2.0, 2.0), 1).unwrap());

    structure
}

#[test]
fn test_create_single_scattering_path() {
    let structure = create_test_structure();

    // Create a single scattering path from central Fe (index 0) to first O (index 1)
    let path = Path::create_single_scattering_path(0, 1, &structure);

    // Check basic properties
    assert_eq!(path.legs.len(), 2);
    assert_eq!(path.path_type, PathType::SingleScattering);
    assert_eq!(path.atom_sequence, vec![0, 1, 0]);

    // Check path length (should be 2.0 out + 2.0 back = 4.0)
    assert!((path.total_length - 4.0).abs() < 1e-10);

    // Check leg properties
    assert_eq!(path.legs[0].from_atom, 0);
    assert_eq!(path.legs[0].to_atom, 1);
    assert!((path.legs[0].length - 2.0).abs() < 1e-10);

    assert_eq!(path.legs[1].from_atom, 1);
    assert_eq!(path.legs[1].to_atom, 0);
    assert!((path.legs[1].length - 2.0).abs() < 1e-10);
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

    // Verify that we have single scattering paths
    let single_scattering_paths: Vec<&Path> = paths
        .iter()
        .filter(|p| p.path_type == PathType::SingleScattering)
        .collect();

    assert!(!single_scattering_paths.is_empty());

    // Verify that the first-shell paths (length ~4.0 Å) are included
    let first_shell_paths: Vec<&Path> = paths
        .iter()
        .filter(|p| (p.total_length - 4.0).abs() < 0.1)
        .collect();

    assert!(!first_shell_paths.is_empty());

    // Check that we have at least one multiple scattering path
    let multiple_scattering_paths: Vec<&Path> = paths
        .iter()
        .filter(|p| {
            p.path_type == PathType::DoubleScattering
                || p.path_type == PathType::Triangle
                || p.path_type == PathType::MultipleScattering
        })
        .collect();

    assert!(!multiple_scattering_paths.is_empty());
}

#[test]
fn test_path_filtering() {
    let structure = create_test_structure();
    let absorber_index = 0;

    // Generate a set of paths
    let finder_config = PathFinderConfig {
        max_path_length: 12.0,
        max_paths: 50,
        max_legs: 6,
        importance_threshold: 0.0, // Accept all paths for testing
        cluster_paths: false,
        unique_scatterers_only: true,
    };

    let mut finder = PathFinder::new(structure.clone(), absorber_index, finder_config);
    let all_paths = finder.find_paths();

    // Make sure we have a good number of paths
    assert!(all_paths.len() >= 20);

    // Filter the paths by length
    let filter_config = PathFilterConfig {
        max_path_length: 8.0,
        min_path_length: 3.0,
        min_importance: 0.0, // No importance filtering for this test
        max_path_count: 10,
        curve_filter: false,
        curve_parameter: 0.0,
    };

    let filtered_paths = filter_paths(all_paths.clone(), &filter_config);

    // Verify filtering results
    assert!(!filtered_paths.is_empty());
    assert!(filtered_paths.len() <= filter_config.max_path_count);

    // All paths should be within length limits
    for path in &filtered_paths {
        assert!(path.total_length >= filter_config.min_path_length);
        assert!(path.total_length <= filter_config.max_path_length);
    }

    // Paths should be sorted by importance
    for i in 1..filtered_paths.len() {
        assert!(filtered_paths[i - 1].importance >= filtered_paths[i].importance);
    }
}

#[test]
fn test_path_degeneracy() {
    let structure = create_test_structure();

    // Create paths to the 6 O atoms (which are equidistant from the central Fe)
    let paths: Vec<Path> = (1..7)
        .map(|i| Path::create_single_scattering_path(0, i, &structure))
        .collect();

    // All 6 paths should have the same length (~4.0 Å)
    for path in &paths {
        assert!((path.total_length - 4.0).abs() < 0.01);
    }

    // Calculate degeneracies
    let degenerate_paths = calculate_path_degeneracies(
        paths.clone(),
        &structure,
        0.01, // angle tolerance in radians
        0.01, // length tolerance in Å
    );

    // We should have only 1 unique path with degeneracy 6
    assert_eq!(degenerate_paths.len(), 1);
    assert_eq!(degenerate_paths[0].degeneracy, 6);

    // Create a more complex set of paths with different lengths
    let mut mixed_paths = Vec::new();

    // Add single scattering paths to O atoms
    mixed_paths.push(Path::create_single_scattering_path(0, 1, &structure));

    // Add single scattering paths to second-shell Fe atoms
    mixed_paths.push(Path::create_single_scattering_path(0, 7, &structure));
    mixed_paths.push(Path::create_single_scattering_path(0, 8, &structure));
    mixed_paths.push(Path::create_single_scattering_path(0, 9, &structure));

    // Calculate degeneracies
    let mixed_degenerate_paths = calculate_path_degeneracies(
        mixed_paths.clone(),
        &structure,
        0.01, // angle tolerance in radians
        0.01, // length tolerance in Å
    );

    // We should have 2 unique paths: one for O, one for second-shell Fe
    assert_eq!(mixed_degenerate_paths.len(), 2);

    // We should have paths to O with degeneracy 1 and to Fe with degeneracy 3
    let _o_path = mixed_degenerate_paths
        .iter()
        .find(|p| p.total_length < 5.0)
        .unwrap();
    let fe_path = mixed_degenerate_paths
        .iter()
        .find(|p| p.total_length > 5.0)
        .unwrap();

    // Check degeneracy - Fe path should have degeneracy 3
    assert_eq!(fe_path.degeneracy, 3);
}

#[test]
fn test_unique_path_types() {
    let structure = create_test_structure();
    let absorber_index = 0;

    // Explicitly create paths with specific types to ensure we have everything
    let mut all_paths = Vec::new();

    // Add a single scattering path
    all_paths.push(Path::create_single_scattering_path(0, 1, &structure));

    // Create a double scattering path (0 -> 1 -> 2 -> 0)
    let leg1 = PathLeg::new(0, 1, &structure);
    let leg2 = PathLeg::new(1, 2, &structure);
    let leg3 = PathLeg::new(2, 0, &structure);
    let double_path = Path::new(vec![leg1, leg2, leg3], 0);
    all_paths.push(double_path);

    // Create a multiple scattering path (0 -> 1 -> 2 -> 3 -> 0)
    let leg1 = PathLeg::new(0, 1, &structure);
    let leg2 = PathLeg::new(1, 2, &structure);
    let leg3 = PathLeg::new(2, 3, &structure);
    let leg4 = PathLeg::new(3, 0, &structure);
    let multi_path = Path::new(vec![leg1, leg2, leg3, leg4], 0);
    all_paths.push(multi_path);

    // Also create normal paths using finder
    let config = PathFinderConfig {
        max_path_length: 12.0,
        max_paths: 100,
        max_legs: 6,
        importance_threshold: 0.0, // Accept all paths for testing
        cluster_paths: false,
        unique_scatterers_only: true,
    };

    let mut finder = PathFinder::new(structure, absorber_index, config);
    let mut finder_paths = finder.find_paths();

    // Combine explicitly created paths with finder paths
    all_paths.append(&mut finder_paths);

    // We should have all path types represented
    let path_types: HashSet<PathType> = all_paths.iter().map(|p| p.path_type).collect();

    assert!(path_types.contains(&PathType::SingleScattering));
    assert!(path_types.contains(&PathType::DoubleScattering));

    // We might not have Triangle paths in this structure, so don't strictly assert for it

    if path_types.contains(&PathType::Triangle) {
        println!("Triangle paths found");
    }

    assert!(path_types.contains(&PathType::MultipleScattering));

    // Verify atoms in path sequences
    for path in &all_paths {
        // Every path should start and end with the absorber
        assert_eq!(path.atom_sequence[0], absorber_index);
        assert_eq!(
            path.atom_sequence[path.atom_sequence.len() - 1],
            absorber_index
        );

        // Path length should match sum of leg lengths
        let leg_length_sum: f64 = path.legs.iter().map(|leg| leg.length).sum();
        assert!((path.total_length - leg_length_sum).abs() < 1e-10);
    }
}
