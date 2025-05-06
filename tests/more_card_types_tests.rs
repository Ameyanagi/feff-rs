/*
MIT License with FEFF10 Attribution

Copyright (c) 2025 Ameyanagi

Based on or developed using Distribution: FEFF10.0
Copyright (c) 2020 FEFF Project, University of Washington and SLAC National Accelerator Laboratory.
All rights reserved.
*/

use feff_rs::input::{
    parameters::{CifsParams, DosParams, RestartParams},
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
fn test_parse_restart_basic() {
    let content = r#"TITLE Test RESTART card

ATOMS
0 0.0 0.0 0.0 Fe

RESTART pot

POTENTIALS
0 26 Fe
"#;

    let (_dir, file_path) = create_test_feff_input(content);
    let result = parse_feff_input(&file_path);

    assert!(result.is_ok());
    let input = result.unwrap();

    // Check RESTART parameters
    assert!(input.restart.is_some());
    let restart = input.restart.unwrap();
    assert_eq!(restart.module, "pot");
    assert!(restart.file_name.is_none());
}

#[test]
fn test_parse_restart_with_filename() {
    let content = r#"TITLE Test RESTART card with filename

ATOMS
0 0.0 0.0 0.0 Fe

RESTART pot custom_pot.bin

POTENTIALS
0 26 Fe
"#;

    let (_dir, file_path) = create_test_feff_input(content);
    let result = parse_feff_input(&file_path);

    assert!(result.is_ok());
    let input = result.unwrap();

    // Check RESTART parameters with filename
    assert!(input.restart.is_some());
    let restart = input.restart.unwrap();
    assert_eq!(restart.module, "pot");
    assert_eq!(restart.file_name, Some("custom_pot.bin".to_string()));
}

#[test]
fn test_parse_restart_empty() {
    let content = r#"TITLE Test RESTART card with empty content

ATOMS
0 0.0 0.0 0.0 Fe

RESTART

POTENTIALS
0 26 Fe
"#;

    let (_dir, file_path) = create_test_feff_input(content);
    let result = parse_feff_input(&file_path);

    // Should fail because RESTART requires a module name
    assert!(result.is_err());
    match result {
        Err(InputError::ParseError(msg)) => {
            assert!(msg.contains("RESTART card requires a module name"));
        }
        _ => panic!("Expected ParseError for empty RESTART card"),
    }
}

#[test]
fn test_write_restart() {
    // Create a minimal valid input file
    let content = r#"TITLE Test RESTART write test

ATOMS
0 0.0 0.0 0.0 Fe

POTENTIALS
0 26 Fe
"#;

    let (_dir, file_path) = create_test_feff_input(content);
    let result = parse_feff_input(&file_path);
    assert!(result.is_ok());

    // Get the input and add RESTART parameters
    let mut input = result.unwrap();

    // Add RESTART parameters
    let restart_params = RestartParams {
        module: "fms".to_string(),
        file_name: Some("custom_fms.bin".to_string()),
    };
    input.restart = Some(restart_params);

    // Create a temporary directory for the output
    let dir = tempdir().unwrap();
    let output_path = dir.path().join("test_write_restart.inp");

    // Write to file
    assert!(input.write(&output_path).is_ok());

    // Read back and verify
    let read_result = parse_feff_input(&output_path);
    assert!(read_result.is_ok());

    let read_input = read_result.unwrap();
    assert!(read_input.restart.is_some());

    let read_restart = read_input.restart.unwrap();
    assert_eq!(read_restart.module, "fms");
    assert_eq!(read_restart.file_name, Some("custom_fms.bin".to_string()));
}

#[test]
fn test_parse_dos_basic() {
    let content = r#"TITLE Test DOS card

ATOMS
0 0.0 0.0 0.0 Fe

DOS -25.0 25.0 0.2 0.3

POTENTIALS
0 26 Fe
"#;

    let (_dir, file_path) = create_test_feff_input(content);
    let result = parse_feff_input(&file_path);

    assert!(result.is_ok());
    let input = result.unwrap();

    // Check DOS parameters
    assert!(input.dos.is_some());
    let dos = input.dos.unwrap();
    assert_eq!(dos.emin, -25.0);
    assert_eq!(dos.emax, 25.0);
    assert_eq!(dos.estep, 0.2);
    assert_eq!(dos.gamma, 0.3);
    assert!(dos.params.is_empty());
}

#[test]
fn test_parse_dos_with_additional_params() {
    let content = r#"TITLE Test DOS card with additional parameters

ATOMS
0 0.0 0.0 0.0 Fe

DOS -25.0 25.0 0.2 0.3 1.0 2.0 3.0

POTENTIALS
0 26 Fe
"#;

    let (_dir, file_path) = create_test_feff_input(content);
    let result = parse_feff_input(&file_path);

    assert!(result.is_ok());
    let input = result.unwrap();

    // Check DOS parameters with additional parameters
    assert!(input.dos.is_some());
    let dos = input.dos.unwrap();
    assert_eq!(dos.emin, -25.0);
    assert_eq!(dos.emax, 25.0);
    assert_eq!(dos.estep, 0.2);
    assert_eq!(dos.gamma, 0.3);
    assert_eq!(dos.params, vec![1.0, 2.0, 3.0]);
}

#[test]
fn test_parse_dos_defaults() {
    let content = r#"TITLE Test DOS card with defaults

ATOMS
0 0.0 0.0 0.0 Fe

DOS

POTENTIALS
0 26 Fe
"#;

    let (_dir, file_path) = create_test_feff_input(content);
    let result = parse_feff_input(&file_path);

    assert!(result.is_ok());
    let input = result.unwrap();

    // Check DOS parameters with defaults
    assert!(input.dos.is_some());
    let dos = input.dos.unwrap();
    assert_eq!(dos.emin, -20.0);
    assert_eq!(dos.emax, 20.0);
    assert_eq!(dos.estep, 0.1);
    assert_eq!(dos.gamma, 0.2);
    assert!(dos.params.is_empty());
}

#[test]
fn test_parse_dos_invalid() {
    let content = r#"TITLE Test invalid DOS card

ATOMS
0 0.0 0.0 0.0 Fe

DOS 15.0 10.0 0.1 0.2

POTENTIALS
0 26 Fe
"#;

    let (_dir, file_path) = create_test_feff_input(content);
    let result = parse_feff_input(&file_path);

    // Should fail due to invalid parameters
    assert!(result.is_err());
    match result {
        Err(InputError::InvalidFormat(msg)) => {
            assert!(msg.contains("Maximum energy") && msg.contains("minimum energy"));
        }
        _ => panic!("Expected InvalidFormat error for invalid DOS parameters"),
    }
}

#[test]
fn test_write_dos() {
    // Create a minimal valid input file
    let content = r#"TITLE Test DOS write test

ATOMS
0 0.0 0.0 0.0 Fe

POTENTIALS
0 26 Fe
"#;

    let (_dir, file_path) = create_test_feff_input(content);
    let result = parse_feff_input(&file_path);
    assert!(result.is_ok());

    // Get the input and add DOS parameters
    let mut input = result.unwrap();

    // Add DOS parameters
    let dos_params = DosParams {
        emin: -30.0,
        emax: 30.0,
        estep: 0.25,
        gamma: 0.4,
        params: vec![1.5, 2.5],
    };
    input.dos = Some(dos_params);

    // Create a temporary directory for the output
    let dir = tempdir().unwrap();
    let output_path = dir.path().join("test_write_dos.inp");

    // Write to file
    assert!(input.write(&output_path).is_ok());

    // Read back and verify
    let read_result = parse_feff_input(&output_path);
    assert!(read_result.is_ok());

    let read_input = read_result.unwrap();
    assert!(read_input.dos.is_some());

    let read_dos = read_input.dos.unwrap();
    assert_eq!(read_dos.emin, -30.0);
    assert_eq!(read_dos.emax, 30.0);
    assert_eq!(read_dos.estep, 0.25);
    assert_eq!(read_dos.gamma, 0.4);
    assert_eq!(read_dos.params, vec![1.5, 2.5]);
}

#[test]
fn test_parse_cifs_basic() {
    let content = r#"TITLE Test CIFS card

ATOMS
0 0.0 0.0 0.0 Fe

CIFS structure.cif

POTENTIALS
0 26 Fe
"#;

    let (_dir, file_path) = create_test_feff_input(content);
    let result = parse_feff_input(&file_path);

    assert!(result.is_ok());
    let input = result.unwrap();

    // Check CIFS parameters
    assert!(input.cifs.is_some());
    let cifs = input.cifs.unwrap();
    assert_eq!(cifs.file_name, "structure.cif");
    assert!(cifs.site_index.is_none());
    assert!(cifs.distance_cutoff.is_none());
}

#[test]
fn test_parse_cifs_with_site_index() {
    let content = r#"TITLE Test CIFS card with site index

ATOMS
0 0.0 0.0 0.0 Fe

CIFS structure.cif 2

POTENTIALS
0 26 Fe
"#;

    let (_dir, file_path) = create_test_feff_input(content);
    let result = parse_feff_input(&file_path);

    assert!(result.is_ok());
    let input = result.unwrap();

    // Check CIFS parameters with site index
    assert!(input.cifs.is_some());
    let cifs = input.cifs.unwrap();
    assert_eq!(cifs.file_name, "structure.cif");
    assert_eq!(cifs.site_index, Some(2));
    assert!(cifs.distance_cutoff.is_none());
}

#[test]
fn test_parse_cifs_full() {
    let content = r#"TITLE Test CIFS card with all parameters

ATOMS
0 0.0 0.0 0.0 Fe

CIFS structure.cif 2 7.5

POTENTIALS
0 26 Fe
"#;

    let (_dir, file_path) = create_test_feff_input(content);
    let result = parse_feff_input(&file_path);

    assert!(result.is_ok());
    let input = result.unwrap();

    // Check CIFS parameters with all parameters
    assert!(input.cifs.is_some());
    let cifs = input.cifs.unwrap();
    assert_eq!(cifs.file_name, "structure.cif");
    assert_eq!(cifs.site_index, Some(2));
    assert_eq!(cifs.distance_cutoff, Some(7.5));
}

#[test]
fn test_parse_cifs_empty() {
    let content = r#"TITLE Test CIFS card with empty content

ATOMS
0 0.0 0.0 0.0 Fe

CIFS

POTENTIALS
0 26 Fe
"#;

    let (_dir, file_path) = create_test_feff_input(content);
    let result = parse_feff_input(&file_path);

    // Should fail because CIFS requires a file name
    assert!(result.is_err());
    match result {
        Err(InputError::ParseError(msg)) => {
            assert!(msg.contains("CIFS card requires a file name"));
        }
        _ => panic!("Expected ParseError for empty CIFS card"),
    }
}

#[test]
fn test_parse_cifs_invalid_site() {
    let content = r#"TITLE Test CIFS card with invalid site index

ATOMS
0 0.0 0.0 0.0 Fe

CIFS structure.cif x

POTENTIALS
0 26 Fe
"#;

    let (_dir, file_path) = create_test_feff_input(content);
    let result = parse_feff_input(&file_path);

    // Should fail due to invalid site index
    assert!(result.is_err());
    match result {
        Err(InputError::ParseError(msg)) => {
            assert!(msg.contains("Invalid site index"));
        }
        _ => panic!("Expected ParseError for invalid site index"),
    }
}

#[test]
fn test_write_cifs() {
    // Create a minimal valid input file
    let content = r#"TITLE Test CIFS write test

ATOMS
0 0.0 0.0 0.0 Fe

POTENTIALS
0 26 Fe
"#;

    let (_dir, file_path) = create_test_feff_input(content);
    let result = parse_feff_input(&file_path);
    assert!(result.is_ok());

    // Get the input and add CIFS parameters
    let mut input = result.unwrap();

    // Add CIFS parameters
    let cifs_params = CifsParams {
        file_name: "complex.cif".to_string(),
        site_index: Some(3),
        distance_cutoff: Some(8.2),
    };
    input.cifs = Some(cifs_params);

    // Create a temporary directory for the output
    let dir = tempdir().unwrap();
    let output_path = dir.path().join("test_write_cifs.inp");

    // Write to file
    assert!(input.write(&output_path).is_ok());

    // Read back and verify
    let read_result = parse_feff_input(&output_path);
    assert!(read_result.is_ok());

    let read_input = read_result.unwrap();
    assert!(read_input.cifs.is_some());

    let read_cifs = read_input.cifs.unwrap();
    assert_eq!(read_cifs.file_name, "complex.cif");
    assert_eq!(read_cifs.site_index, Some(3));
    assert_eq!(read_cifs.distance_cutoff, Some(8.2));
}

#[test]
fn test_combined_new_cards() {
    let content = r#"TITLE Test combined new card types

ATOMS
0 0.0 0.0 0.0 Fe

RESTART pot saved_pot.bin
DOS -30.0 30.0 0.2 0.4 1.0 2.0
CIFS structure.cif 2 8.0

POTENTIALS
0 26 Fe
"#;

    let (_dir, file_path) = create_test_feff_input(content);
    let result = parse_feff_input(&file_path);

    assert!(result.is_ok());
    let input = result.unwrap();

    // Check all cards are parsed
    assert!(input.restart.is_some());
    assert!(input.dos.is_some());
    assert!(input.cifs.is_some());

    // Verify RESTART parameters
    let restart = input.restart.unwrap();
    assert_eq!(restart.module, "pot");
    assert_eq!(restart.file_name, Some("saved_pot.bin".to_string()));

    // Verify DOS parameters
    let dos = input.dos.unwrap();
    assert_eq!(dos.emin, -30.0);
    assert_eq!(dos.emax, 30.0);
    assert_eq!(dos.estep, 0.2);
    assert_eq!(dos.gamma, 0.4);
    assert_eq!(dos.params, vec![1.0, 2.0]);

    // Verify CIFS parameters
    let cifs = input.cifs.unwrap();
    assert_eq!(cifs.file_name, "structure.cif");
    assert_eq!(cifs.site_index, Some(2));
    assert_eq!(cifs.distance_cutoff, Some(8.0));
}

#[test]
fn test_roundtrip_all_new_cards() {
    // Create an input object with all new card types
    let mut input = FeffInput::new();

    // Add basic required data
    input.title = Some("Test roundtrip for new card types".to_string());

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

    // Add RESTART parameters
    input.restart = Some(RestartParams {
        module: "pot".to_string(),
        file_name: Some("saved_pot.bin".to_string()),
    });

    // Add DOS parameters
    input.dos = Some(DosParams {
        emin: -30.0,
        emax: 30.0,
        estep: 0.2,
        gamma: 0.4,
        params: vec![1.0, 2.0, 3.0],
    });

    // Add CIFS parameters
    input.cifs = Some(CifsParams {
        file_name: "structure.cif".to_string(),
        site_index: Some(2),
        distance_cutoff: Some(8.0),
    });

    // Write to file
    let dir = tempdir().unwrap();
    let file_path = dir.path().join("feff_roundtrip_new.inp");
    let write_result = input.write(&file_path);
    assert!(write_result.is_ok());

    // Read back from file
    let read_result = parse_feff_input(&file_path);
    assert!(read_result.is_ok());
    let read_input = read_result.unwrap();

    // Verify RESTART was preserved
    assert!(read_input.restart.is_some());
    let restart = read_input.restart.unwrap();
    assert_eq!(restart.module, "pot");
    assert_eq!(restart.file_name, Some("saved_pot.bin".to_string()));

    // Verify DOS was preserved
    assert!(read_input.dos.is_some());
    let dos = read_input.dos.unwrap();
    assert_eq!(dos.emin, -30.0);
    assert_eq!(dos.emax, 30.0);
    assert_eq!(dos.estep, 0.2);
    assert_eq!(dos.gamma, 0.4);
    assert_eq!(dos.params, vec![1.0, 2.0, 3.0]);

    // Verify CIFS was preserved
    assert!(read_input.cifs.is_some());
    let cifs = read_input.cifs.unwrap();
    assert_eq!(cifs.file_name, "structure.cif");
    assert_eq!(cifs.site_index, Some(2));
    assert_eq!(cifs.distance_cutoff, Some(8.0));
}
