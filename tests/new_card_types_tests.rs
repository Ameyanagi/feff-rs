/*
MIT License with FEFF10 Attribution

Copyright (c) 2025 Ameyanagi

Based on or developed using Distribution: FEFF10.0
Copyright (c) 2020 FEFF Project, University of Washington and SLAC National Accelerator Laboratory.
All rights reserved.
*/

use feff_rs::input::{
    parameters::{BandstructureParams, KmeshParams, RdinpParams},
    parse_feff_input, FeffInput, InputError,
};
use std::fs::File;
use std::io::Write;
use tempfile::tempdir;

/// Test helper to create a temporary test FEFF input file
fn create_test_feff_input(content: &str) -> (tempfile::TempDir, std::path::PathBuf) {
    let dir = tempdir().unwrap();
    let file_path = dir.path().join("feff.inp");
    let mut file = File::create(&file_path).unwrap();
    writeln!(file, "{}", content).unwrap();
    (dir, file_path)
}

#[test]
fn test_parse_rdinp_basic() {
    let content = r#"TITLE Test RDINP card

ATOMS
0 0.0 0.0 0.0 Fe

RDINP some_other_file.inp

POTENTIALS
0 26 Fe
"#;

    let (_dir, file_path) = create_test_feff_input(content);
    let result = parse_feff_input(&file_path);

    assert!(result.is_ok());
    let input = result.unwrap();

    // Check RDINP parameters
    assert!(input.rdinp.is_some());
    let rdinp = input.rdinp.unwrap();
    assert_eq!(rdinp.file_name, "some_other_file.inp");
}

#[test]
fn test_parse_rdinp_empty() {
    let content = r#"TITLE Test RDINP card with empty content

ATOMS
0 0.0 0.0 0.0 Fe

RDINP

POTENTIALS
0 26 Fe
"#;

    let (_dir, file_path) = create_test_feff_input(content);
    let result = parse_feff_input(&file_path);

    // Should fail because RDINP requires a file name
    assert!(result.is_err());
    match result {
        Err(InputError::ParseError(msg)) => {
            assert!(msg.contains("RDINP card requires a file name"));
        }
        _ => panic!("Expected ParseError for empty RDINP card"),
    }
}

#[test]
fn test_write_rdinp() {
    // Create a minimal valid input file
    let content = r#"TITLE Test RDINP write test

ATOMS
0 0.0 0.0 0.0 Fe

POTENTIALS
0 26 Fe
"#;

    let (_dir, file_path) = create_test_feff_input(content);
    let result = parse_feff_input(&file_path);
    assert!(result.is_ok());

    // Get the input and add RDINP parameters
    let mut input = result.unwrap();

    // Add RDINP parameters
    let rdinp_params = RdinpParams {
        file_name: "external_input.inp".to_string(),
    };
    input.rdinp = Some(rdinp_params);

    // Create a temporary directory for the output
    let dir = tempdir().unwrap();
    let output_path = dir.path().join("test_write_rdinp.inp");

    // Write to file
    assert!(input.write(&output_path).is_ok());

    // Read back and verify
    let read_result = parse_feff_input(&output_path);
    assert!(read_result.is_ok());

    let read_input = read_result.unwrap();
    assert!(read_input.rdinp.is_some());

    let read_rdinp = read_input.rdinp.unwrap();
    assert_eq!(read_rdinp.file_name, "external_input.inp");
}

#[test]
fn test_parse_bandstructure_basic() {
    let content = r#"TITLE Test BANDSTRUCTURE card

ATOMS
0 0.0 0.0 0.0 Fe

BANDSTRUCTURE 200 -15.0 15.0 0.2 1 1

POTENTIALS
0 26 Fe
"#;

    let (_dir, file_path) = create_test_feff_input(content);
    let result = parse_feff_input(&file_path);

    assert!(result.is_ok());
    let input = result.unwrap();

    // Check BANDSTRUCTURE parameters
    assert!(input.bandstructure.is_some());
    let bandstructure = input.bandstructure.unwrap();
    assert_eq!(bandstructure.nk, 200);
    assert_eq!(bandstructure.emin, -15.0);
    assert_eq!(bandstructure.emax, 15.0);
    assert_eq!(bandstructure.estep, 0.2);
    assert_eq!(bandstructure.kmesh, 1);
    assert_eq!(bandstructure.symmetry, 1);
}

#[test]
fn test_parse_bandstructure_defaults() {
    let content = r#"TITLE Test BANDSTRUCTURE card with defaults

ATOMS
0 0.0 0.0 0.0 Fe

BANDSTRUCTURE

POTENTIALS
0 26 Fe
"#;

    let (_dir, file_path) = create_test_feff_input(content);
    let result = parse_feff_input(&file_path);

    assert!(result.is_ok());
    let input = result.unwrap();

    // Check BANDSTRUCTURE parameters with defaults
    assert!(input.bandstructure.is_some());
    let bandstructure = input.bandstructure.unwrap();
    assert_eq!(bandstructure.nk, 100);
    assert_eq!(bandstructure.emin, -10.0);
    assert_eq!(bandstructure.emax, 10.0);
    assert_eq!(bandstructure.estep, 0.1);
    assert_eq!(bandstructure.kmesh, 1);
    assert_eq!(bandstructure.symmetry, 1);
}

#[test]
fn test_parse_bandstructure_invalid() {
    let content = r#"TITLE Test invalid BANDSTRUCTURE card

ATOMS
0 0.0 0.0 0.0 Fe

BANDSTRUCTURE 0 -15.0 -20.0 -0.1 3 3

POTENTIALS
0 26 Fe
"#;

    let (_dir, file_path) = create_test_feff_input(content);
    let result = parse_feff_input(&file_path);

    // Should fail due to invalid parameters
    assert!(result.is_err());
    match result {
        Err(InputError::InvalidFormat(msg)) => {
            assert!(
                msg.contains("Invalid number of k-points")
                    || msg.contains("Maximum energy")
                    || msg.contains("Energy step")
                    || msg.contains("Invalid kmesh flag")
                    || msg.contains("Invalid symmetry flag")
            );
        }
        _ => panic!("Expected InvalidFormat error for invalid BANDSTRUCTURE parameters"),
    }
}

#[test]
fn test_write_bandstructure() {
    // Create a minimal valid input file
    let content = r#"TITLE Test BANDSTRUCTURE write test

ATOMS
0 0.0 0.0 0.0 Fe

POTENTIALS
0 26 Fe
"#;

    let (_dir, file_path) = create_test_feff_input(content);
    let result = parse_feff_input(&file_path);
    assert!(result.is_ok());

    // Get the input and add BANDSTRUCTURE parameters
    let mut input = result.unwrap();

    // Add BANDSTRUCTURE parameters
    let bandstructure_params = BandstructureParams {
        nk: 150,
        emin: -12.5,
        emax: 12.5,
        estep: 0.15,
        kmesh: 1,
        symmetry: 1,
    };
    input.bandstructure = Some(bandstructure_params);

    // Create a temporary directory for the output
    let dir = tempdir().unwrap();
    let output_path = dir.path().join("test_write_bandstructure.inp");

    // Write to file
    assert!(input.write(&output_path).is_ok());

    // Read back and verify
    let read_result = parse_feff_input(&output_path);
    assert!(read_result.is_ok());

    let read_input = read_result.unwrap();
    assert!(read_input.bandstructure.is_some());

    let read_bandstructure = read_input.bandstructure.unwrap();
    assert_eq!(read_bandstructure.nk, 150);
    assert_eq!(read_bandstructure.emin, -12.5);
    assert_eq!(read_bandstructure.emax, 12.5);
    assert_eq!(read_bandstructure.estep, 0.15);
    assert_eq!(read_bandstructure.kmesh, 1);
    assert_eq!(read_bandstructure.symmetry, 1);
}

#[test]
fn test_parse_kmesh_basic() {
    let content = r#"TITLE Test KMESH card

ATOMS
0 0.0 0.0 0.0 Fe

KMESH 20 20 20
0.0 0.0 0.0
0.5 0.0 0.0
0.5 0.5 0.0
0.0 0.5 0.0

POTENTIALS
0 26 Fe
"#;

    let (_dir, file_path) = create_test_feff_input(content);
    let result = parse_feff_input(&file_path);

    assert!(result.is_ok());
    let input = result.unwrap();

    // Check KMESH parameters
    assert!(input.kmesh.is_some());
    let kmesh = input.kmesh.unwrap();
    assert_eq!(kmesh.nx, 20);
    assert_eq!(kmesh.ny, 20);
    assert_eq!(kmesh.nz, 20);

    // Check k-points
    assert_eq!(kmesh.kpoints.len(), 4);
    assert_eq!(kmesh.kpoints[0], (0.0, 0.0, 0.0));
    assert_eq!(kmesh.kpoints[1], (0.5, 0.0, 0.0));
    assert_eq!(kmesh.kpoints[2], (0.5, 0.5, 0.0));
    assert_eq!(kmesh.kpoints[3], (0.0, 0.5, 0.0));
}

#[test]
fn test_parse_kmesh_invalid() {
    let content = r#"TITLE Test invalid KMESH card

ATOMS
0 0.0 0.0 0.0 Fe

KMESH -5 10 10

POTENTIALS
0 26 Fe
"#;

    let (_dir, file_path) = create_test_feff_input(content);
    let result = parse_feff_input(&file_path);

    // Should fail due to invalid parameters
    assert!(result.is_err());
    match result {
        Err(InputError::InvalidFormat(msg)) => {
            assert!(msg.contains("Mesh dimensions must be positive"));
        }
        _ => panic!("Expected InvalidFormat error for invalid KMESH parameters"),
    }
}

#[test]
fn test_write_kmesh() {
    // Create a minimal valid input file
    let content = r#"TITLE Test KMESH write test

ATOMS
0 0.0 0.0 0.0 Fe

POTENTIALS
0 26 Fe
"#;

    let (_dir, file_path) = create_test_feff_input(content);
    let result = parse_feff_input(&file_path);
    assert!(result.is_ok());

    // Get the input and add KMESH parameters
    let mut input = result.unwrap();

    // Add KMESH parameters
    let kmesh_params = KmeshParams {
        nx: 15,
        ny: 15,
        nz: 15,
        kpoints: vec![(0.0, 0.0, 0.0), (0.25, 0.25, 0.25), (0.5, 0.5, 0.5)],
    };
    input.kmesh = Some(kmesh_params);

    // Create a temporary directory for the output
    let dir = tempdir().unwrap();
    let output_path = dir.path().join("test_write_kmesh.inp");

    // Write to file
    assert!(input.write(&output_path).is_ok());

    // Read back and verify
    let read_result = parse_feff_input(&output_path);
    assert!(read_result.is_ok());

    let read_input = read_result.unwrap();
    assert!(read_input.kmesh.is_some());

    let read_kmesh = read_input.kmesh.unwrap();
    assert_eq!(read_kmesh.nx, 15);
    assert_eq!(read_kmesh.ny, 15);
    assert_eq!(read_kmesh.nz, 15);
    assert_eq!(read_kmesh.kpoints.len(), 3);
    assert_eq!(read_kmesh.kpoints[0], (0.0, 0.0, 0.0));
    assert_eq!(read_kmesh.kpoints[1], (0.25, 0.25, 0.25));
    assert_eq!(read_kmesh.kpoints[2], (0.5, 0.5, 0.5));
}

#[test]
fn test_combined_bandstructure_kmesh() {
    let content = r#"TITLE Test combined BANDSTRUCTURE and KMESH cards

ATOMS
0 0.0 0.0 0.0 Fe

BANDSTRUCTURE 200 -15.0 15.0 0.2 0 1
KMESH 20 20 20
0.0 0.0 0.0
0.5 0.0 0.0
0.5 0.5 0.0
0.0 0.5 0.0

POTENTIALS
0 26 Fe
"#;

    let (_dir, file_path) = create_test_feff_input(content);
    let result = parse_feff_input(&file_path);

    assert!(result.is_ok());
    let input = result.unwrap();

    // Check both cards are parsed
    assert!(input.bandstructure.is_some());
    assert!(input.kmesh.is_some());

    let bandstructure = input.bandstructure.unwrap();
    let kmesh = input.kmesh.unwrap();

    // Verify BANDSTRUCTURE parameters
    assert_eq!(bandstructure.nk, 200);
    assert_eq!(bandstructure.emin, -15.0);
    assert_eq!(bandstructure.emax, 15.0);
    assert_eq!(bandstructure.estep, 0.2);
    assert_eq!(bandstructure.kmesh, 0); // User-defined k-mesh
    assert_eq!(bandstructure.symmetry, 1);

    // Verify KMESH parameters
    assert_eq!(kmesh.nx, 20);
    assert_eq!(kmesh.ny, 20);
    assert_eq!(kmesh.nz, 20);
    assert_eq!(kmesh.kpoints.len(), 4);
}

#[test]
fn test_write_read_roundtrip_all() {
    // Create an input object with all our new card types
    let mut input = FeffInput::new();

    // Add basic required data
    input.title = Some("Test roundtrip for advanced cards".to_string());

    // We need atoms for a valid FEFF input
    let mut atomic_structure = feff_rs::atoms::AtomicStructure::new();
    let atom = feff_rs::atoms::Atom::new(
        26, // Fe
        feff_rs::atoms::Vector3D::new(0.0, 0.0, 0.0),
        0, // potential type
    )
    .unwrap();
    atomic_structure.add_atom(atom);
    input.atomic_structure = Some(atomic_structure);

    // Add potential for Fe
    let pot_info = feff_rs::input::parameters::PotentialInfo {
        index: 0,
        atomic_number: 26,
        symbol: "Fe".to_string(),
    };
    input.potentials.insert(0, pot_info);

    // Add RDINP parameters
    input.rdinp = Some(RdinpParams {
        file_name: "other_input.inp".to_string(),
    });

    // Add BANDSTRUCTURE parameters
    input.bandstructure = Some(BandstructureParams {
        nk: 150,
        emin: -12.5,
        emax: 12.5,
        estep: 0.15,
        kmesh: 1,
        symmetry: 1,
    });

    // Add KMESH parameters
    input.kmesh = Some(KmeshParams {
        nx: 15,
        ny: 15,
        nz: 15,
        kpoints: vec![(0.0, 0.0, 0.0), (0.25, 0.25, 0.25), (0.5, 0.5, 0.5)],
    });

    // Write to file
    let dir = tempdir().unwrap();
    let file_path = dir.path().join("feff_roundtrip_all.inp");
    let write_result = input.write(&file_path);
    assert!(write_result.is_ok());

    // Read back from file
    let read_result = parse_feff_input(&file_path);
    assert!(read_result.is_ok());
    let read_input = read_result.unwrap();

    // Verify RDINP was preserved
    assert!(read_input.rdinp.is_some());
    let rdinp = read_input.rdinp.unwrap();
    assert_eq!(rdinp.file_name, "other_input.inp");

    // Verify BANDSTRUCTURE was preserved
    assert!(read_input.bandstructure.is_some());
    let bandstructure = read_input.bandstructure.unwrap();
    assert_eq!(bandstructure.nk, 150);
    assert_eq!(bandstructure.emin, -12.5);
    assert_eq!(bandstructure.emax, 12.5);
    assert_eq!(bandstructure.estep, 0.15);
    assert_eq!(bandstructure.kmesh, 1);
    assert_eq!(bandstructure.symmetry, 1);

    // Verify KMESH was preserved
    assert!(read_input.kmesh.is_some());
    let kmesh = read_input.kmesh.unwrap();
    assert_eq!(kmesh.nx, 15);
    assert_eq!(kmesh.ny, 15);
    assert_eq!(kmesh.nz, 15);
    assert_eq!(kmesh.kpoints.len(), 3);
    assert_eq!(kmesh.kpoints[0], (0.0, 0.0, 0.0));
    assert_eq!(kmesh.kpoints[1], (0.25, 0.25, 0.25));
    assert_eq!(kmesh.kpoints[2], (0.5, 0.5, 0.5));
}
