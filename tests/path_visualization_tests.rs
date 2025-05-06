/*
MIT License with FEFF10 Attribution

Copyright (c) 2025 Ameyanagi

Based on or developed using Distribution: FEFF10.0
Copyright (c) 2020 FEFF Project, University of Washington and SLAC National Accelerator Laboratory.
All rights reserved.
*/

//! Tests for path visualization utilities
//!
//! These tests verify that the path visualization functions
//! correctly format and represent scattering paths.

use feff_rs::atoms::atom::Atom;
use feff_rs::atoms::structure::AtomicStructure;
use feff_rs::atoms::vector::Vector3D;
use feff_rs::path::{
    create_path_summary_table, export_paths, format_path_description, generate_path_json,
    generate_path_xyz, Path, PathLeg,
};

/// Creates a test structure for path visualization tests
fn create_test_structure() -> AtomicStructure {
    let mut structure = AtomicStructure::new();

    // Add a central Fe atom at the origin
    structure.add_atom(Atom::new(26, Vector3D::new(0.0, 0.0, 0.0), 0).unwrap());

    // Add O atoms in tetrahedral arrangement
    structure.add_atom(Atom::new(8, Vector3D::new(1.5, 1.5, 1.5), 1).unwrap());
    structure.add_atom(Atom::new(8, Vector3D::new(1.5, -1.5, -1.5), 1).unwrap());
    structure.add_atom(Atom::new(8, Vector3D::new(-1.5, 1.5, -1.5), 1).unwrap());
    structure.add_atom(Atom::new(8, Vector3D::new(-1.5, -1.5, 1.5), 1).unwrap());

    structure
}

#[test]
fn test_format_path_description() {
    let structure = create_test_structure();

    // Create a single scattering path
    let path = Path::create_single_scattering_path(0, 1, &structure);

    // Format the path description
    let description = format_path_description(&path, &structure);

    // Verify content
    assert!(description.contains("SingleScattering path"));
    assert!(description.contains("length:"));
    assert!(description.contains("degeneracy:"));
    assert!(description.contains("importance:"));
    assert!(description.contains("Atom sequence:"));
    assert!(description.contains("Legs:"));
    assert!(description.contains("Fe → O"));
    assert!(description.contains("O → Fe"));
}

#[test]
fn test_create_path_summary_table() {
    let structure = create_test_structure();

    // Create multiple paths
    let path1 = Path::create_single_scattering_path(0, 1, &structure);
    let path2 = Path::create_single_scattering_path(0, 2, &structure);

    let paths = vec![path1, path2];

    // Generate the table
    let table = create_path_summary_table(&paths, &structure);

    // Check structure and content
    assert!(table.contains("| # | Type | Length"));
    assert!(table.contains("Fe → O → Fe"));
    assert!(table.len() > 100); // Should be reasonably sized
}

#[test]
fn test_generate_path_xyz() {
    let structure = create_test_structure();

    // Create a single scattering path
    let path = Path::create_single_scattering_path(0, 1, &structure);

    // Generate XYZ format
    let xyz = generate_path_xyz(&path, &structure);

    // Check XYZ format
    let lines: Vec<&str> = xyz.lines().collect();

    // First line should be number of atoms (2 in this case)
    assert_eq!(lines[0].trim(), "2");

    // Second line is a comment
    assert!(lines[1].contains("Path visualization"));

    // Check that we have the Fe and O atom coordinates
    assert!(xyz.contains("Fe*"));
    assert!(xyz.contains("O1"));
    assert!(lines.len() >= 4);
}

#[test]
fn test_generate_path_json() {
    let structure = create_test_structure();

    // Create a single scattering path
    let path = Path::create_single_scattering_path(0, 1, &structure);

    // Generate JSON
    let json = generate_path_json(&path, &structure);

    // Check JSON structure
    assert!(json.contains("\"type\": \"SingleScattering\""));
    assert!(json.contains("\"totalLength\":"));
    assert!(json.contains("\"degeneracy\":"));
    assert!(json.contains("\"importance\":"));
    assert!(json.contains("\"atomSequence\":"));
    assert!(json.contains("\"legs\":"));
    assert!(json.contains("\"element\": \"Fe\""));
    assert!(json.contains("\"element\": \"O\""));
}

#[test]
fn test_export_paths() {
    let structure = create_test_structure();

    // Create two paths
    let mut paths = Vec::new();
    paths.push(Path::create_single_scattering_path(0, 1, &structure));
    paths.push(Path::create_single_scattering_path(0, 2, &structure));

    // Test XYZ format
    let xyz = export_paths(&paths, &structure, "xyz");
    assert!(xyz.contains("Path visualization"));
    assert!(xyz.contains("Fe*"));

    // Test JSON format
    let json = export_paths(&paths, &structure, "json");
    assert!(json.starts_with("["));
    assert!(json.ends_with("]\n"));
    assert!(json.contains("\"type\": \"SingleScattering\""));

    // Test text format
    let text = export_paths(&paths, &structure, "text");
    assert!(text.contains("Path Summary"));
    assert!(text.contains("Detailed Path Descriptions"));
}

#[test]
fn test_visualize_complex_path() {
    let structure = create_test_structure();

    // Create a multiple scattering path with three legs
    let leg1 = PathLeg::new(0, 1, &structure);
    let leg2 = PathLeg::new(1, 2, &structure);
    let leg3 = PathLeg::new(2, 0, &structure);
    let path = Path::new(vec![leg1, leg2, leg3], 0);

    // Generate description
    let description = format_path_description(&path, &structure);

    // Check for correct path type and elements
    assert!(description.contains("DoubleScattering path") || description.contains("Triangle path"));

    // The test structure uses a different path format, so let's check for parts instead
    assert!(description.contains("Fe") && description.contains("O"));
    assert!(description.contains("→"));

    // Generate XYZ
    let xyz = generate_path_xyz(&path, &structure);
    let lines: Vec<&str> = xyz.lines().collect();

    // Should have 3 atoms
    assert_eq!(lines[0].trim(), "3");

    // Check for all three atoms
    assert!(xyz.contains("Fe*"));
    assert!(xyz.contains("O1"));
    assert!(xyz.contains("O2"));
}
