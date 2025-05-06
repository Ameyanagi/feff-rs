/*
MIT License with FEFF10 Attribution

Copyright (c) 2025 Ameyanagi

Based on or developed using Distribution: FEFF10.0
Copyright (c) 2020 FEFF Project, University of Washington and SLAC National Accelerator Laboratory.
All rights reserved.
*/

use feff_rs::input::{
    parameters::{MultipoleParams, TdldaParams},
    parse_feff_input, InputError,
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
fn test_parse_tdlda_basic() {
    let content = r#"TITLE Test TDLDA card

ATOMS
0 0.0 0.0 0.0 Fe

TDLDA 2 0 -20.0 30.0 0.1 0.2

POTENTIALS
0 26 Fe
"#;

    let (_dir, file_path) = create_test_feff_input(content);
    let result = parse_feff_input(&file_path);

    assert!(result.is_ok());
    let input = result.unwrap();

    assert!(input.tdlda.is_some());
    let tdlda = input.tdlda.unwrap();

    // Verify parsed parameters
    assert_eq!(tdlda.iscreen, 2);
    assert_eq!(tdlda.icalc, 0);
    assert_eq!(tdlda.elow, -20.0);
    assert_eq!(tdlda.ehigh, 30.0);
    assert_eq!(tdlda.estep, 0.1);
    assert_eq!(tdlda.gamma, 0.2);
}

#[test]
fn test_parse_tdlda_with_defaults() {
    let content = r#"TITLE Test TDLDA card with defaults

ATOMS
0 0.0 0.0 0.0 Fe

TDLDA

POTENTIALS
0 26 Fe
"#;

    let (_dir, file_path) = create_test_feff_input(content);
    let result = parse_feff_input(&file_path);

    assert!(result.is_ok());
    let input = result.unwrap();

    assert!(input.tdlda.is_some());
    let tdlda = input.tdlda.unwrap();

    // Verify default parameters
    assert_eq!(tdlda.iscreen, 2); // Default TDLDA
    assert_eq!(tdlda.icalc, 0); // Default SCF+XAS
    assert_eq!(tdlda.elow, -20.0);
    assert_eq!(tdlda.ehigh, 30.0);
    assert_eq!(tdlda.estep, 0.1);
    assert_eq!(tdlda.gamma, 0.1);
}

#[test]
fn test_parse_tdlda_invalid_params() {
    let content = r#"TITLE Test TDLDA card with invalid parameters

ATOMS
0 0.0 0.0 0.0 Fe

TDLDA 5 0 -20.0 30.0 0.1 0.1

POTENTIALS
0 26 Fe
"#;

    let (_dir, file_path) = create_test_feff_input(content);
    let result = parse_feff_input(&file_path);

    // Should fail validation because iscreen is not in [0, 1, 2]
    assert!(result.is_err());
    match result {
        Err(InputError::InvalidFormat(msg)) => {
            assert!(msg.contains("Invalid iscreen value"));
        }
        _ => panic!("Expected InvalidFormat error for invalid iscreen value"),
    }
}

#[test]
fn test_write_tdlda() {
    // Create a minimal valid input file
    let content = r#"TITLE Test TDLDA write test

ATOMS
0 0.0 0.0 0.0 Fe

POTENTIALS
0 26 Fe
"#;

    let (_dir, file_path) = create_test_feff_input(content);
    let result = parse_feff_input(&file_path);
    assert!(result.is_ok());

    // Get the input and add TDLDA parameters
    let mut input = result.unwrap();

    // Add TDLDA parameters
    let tdlda_params = TdldaParams {
        iscreen: 1,
        icalc: 2,
        elow: -15.0,
        ehigh: 25.0,
        estep: 0.2,
        gamma: 0.3,
    };
    input.tdlda = Some(tdlda_params);

    // Create a temporary directory for the output
    let dir = tempdir().unwrap();
    let output_path = dir.path().join("test_write_tdlda.inp");

    // Write to file
    assert!(input.write(&output_path).is_ok());

    // Read back and verify
    let read_result = parse_feff_input(&output_path);
    assert!(read_result.is_ok());

    let read_input = read_result.unwrap();
    assert!(read_input.tdlda.is_some());

    let read_tdlda = read_input.tdlda.unwrap();
    assert_eq!(read_tdlda.iscreen, 1);
    assert_eq!(read_tdlda.icalc, 2);
    assert_eq!(read_tdlda.elow, -15.0);
    assert_eq!(read_tdlda.ehigh, 25.0);
    assert_eq!(read_tdlda.estep, 0.2);
    assert_eq!(read_tdlda.gamma, 0.3);
}

#[test]
fn test_parse_multipole_basic() {
    let content = r#"TITLE Test MULTIPOLE card

ATOMS
0 0.0 0.0 0.0 Fe

MULTIPOLE 4 2 1

POTENTIALS
0 26 Fe
"#;

    let (_dir, file_path) = create_test_feff_input(content);
    let result = parse_feff_input(&file_path);

    assert!(result.is_ok());
    let input = result.unwrap();

    assert!(input.multipole.is_some());
    let multipole = input.multipole.unwrap();

    // Verify parsed parameters
    assert_eq!(multipole.lmax, 4);
    assert_eq!(multipole.morder, 2);
    assert_eq!(multipole.tensor, 1);
}

#[test]
fn test_parse_multipole_with_defaults() {
    let content = r#"TITLE Test MULTIPOLE card with defaults

ATOMS
0 0.0 0.0 0.0 Fe

MULTIPOLE

POTENTIALS
0 26 Fe
"#;

    let (_dir, file_path) = create_test_feff_input(content);
    let result = parse_feff_input(&file_path);

    assert!(result.is_ok());
    let input = result.unwrap();

    assert!(input.multipole.is_some());
    let multipole = input.multipole.unwrap();

    // Verify default parameters
    assert_eq!(multipole.lmax, 3); // Default maximum l
    assert_eq!(multipole.morder, 2); // Default quadrupole
    assert_eq!(multipole.tensor, 0); // Default tensor off
}

#[test]
fn test_parse_multipole_invalid_params() {
    let content = r#"TITLE Test MULTIPOLE card with invalid parameters

ATOMS
0 0.0 0.0 0.0 Fe

MULTIPOLE 4 5 1

POTENTIALS
0 26 Fe
"#;

    let (_dir, file_path) = create_test_feff_input(content);
    let result = parse_feff_input(&file_path);

    // Should fail validation because morder is not in [1, 2, 3, 4]
    assert!(result.is_err());
    match result {
        Err(InputError::InvalidFormat(msg)) => {
            assert!(msg.contains("Invalid morder value"));
        }
        _ => panic!("Expected InvalidFormat error for invalid morder value"),
    }
}

#[test]
fn test_write_multipole() {
    // Create a minimal valid input file
    let content = r#"TITLE Test MULTIPOLE write test

ATOMS
0 0.0 0.0 0.0 Fe

POTENTIALS
0 26 Fe
"#;

    let (_dir, file_path) = create_test_feff_input(content);
    let result = parse_feff_input(&file_path);
    assert!(result.is_ok());

    // Get the input and add MULTIPOLE parameters
    let mut input = result.unwrap();

    // Add MULTIPOLE parameters
    let multipole_params = MultipoleParams {
        lmax: 4,
        morder: 3,
        tensor: 1,
    };
    input.multipole = Some(multipole_params);

    // Create a temporary directory for the output
    let dir = tempdir().unwrap();
    let output_path = dir.path().join("test_write_multipole.inp");

    // Write to file
    assert!(input.write(&output_path).is_ok());

    // Read back and verify
    let read_result = parse_feff_input(&output_path);
    assert!(read_result.is_ok());

    let read_input = read_result.unwrap();
    assert!(read_input.multipole.is_some());

    let read_multipole = read_input.multipole.unwrap();
    assert_eq!(read_multipole.lmax, 4);
    assert_eq!(read_multipole.morder, 3);
    assert_eq!(read_multipole.tensor, 1);
}

#[test]
fn test_parse_combined() {
    let content = r#"TITLE Test combined TDLDA and MULTIPOLE cards

ATOMS
0 0.0 0.0 0.0 Fe

TDLDA 2 0 -20.0 30.0 0.1 0.2
MULTIPOLE 4 2 1

POTENTIALS
0 26 Fe
"#;

    let (_dir, file_path) = create_test_feff_input(content);
    let result = parse_feff_input(&file_path);

    assert!(result.is_ok());
    let input = result.unwrap();

    // Verify both cards are parsed correctly
    assert!(input.tdlda.is_some());
    assert!(input.multipole.is_some());

    let tdlda = input.tdlda.unwrap();
    let multipole = input.multipole.unwrap();

    // Verify TDLDA parameters
    assert_eq!(tdlda.iscreen, 2);
    assert_eq!(tdlda.gamma, 0.2);

    // Verify MULTIPOLE parameters
    assert_eq!(multipole.lmax, 4);
    assert_eq!(multipole.tensor, 1);
}
