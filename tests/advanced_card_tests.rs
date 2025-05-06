/*
MIT License with FEFF10 Attribution

Copyright (c) 2025 Ameyanagi

Based on or developed using Distribution: FEFF10.0
Copyright (c) 2020 FEFF Project, University of Washington and SLAC National Accelerator Laboratory.
All rights reserved.
*/

use feff_rs::input::{
    parameters::{DimensionsParams, ScreenParams, SpectralParams},
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
fn test_parse_screen_basic() {
    let content = r#"TITLE Test SCREEN card

ATOMS
0 0.0 0.0 0.0 Fe

SCREEN 1 1 1.0 1.0

POTENTIALS
0 26 Fe
"#;

    let (_dir, file_path) = create_test_feff_input(content);
    let result = parse_feff_input(&file_path);

    assert!(result.is_ok());
    let input = result.unwrap();

    assert!(input.screen.is_some());
    let screen = input.screen.unwrap();

    // Verify parsed parameters
    assert_eq!(screen.iself, 1);
    assert_eq!(screen.iscreen, 1);
    assert_eq!(screen.ca1, 1.0);
    assert_eq!(screen.ci1, 1.0);
}

#[test]
fn test_parse_screen_with_defaults() {
    let content = r#"TITLE Test SCREEN card with defaults

ATOMS
0 0.0 0.0 0.0 Fe

SCREEN

POTENTIALS
0 26 Fe
"#;

    let (_dir, file_path) = create_test_feff_input(content);
    let result = parse_feff_input(&file_path);

    assert!(result.is_ok());
    let input = result.unwrap();

    assert!(input.screen.is_some());
    let screen = input.screen.unwrap();

    // Verify default parameters
    assert_eq!(screen.iself, 1); // Default HL scheme
    assert_eq!(screen.iscreen, 1); // Default screened
    assert_eq!(screen.ca1, 1.0);
    assert_eq!(screen.ci1, 1.0);
}

#[test]
fn test_parse_screen_invalid_params() {
    let content = r#"TITLE Test SCREEN card with invalid parameters

ATOMS
0 0.0 0.0 0.0 Fe

SCREEN 3 1 1.0 1.0

POTENTIALS
0 26 Fe
"#;

    let (_dir, file_path) = create_test_feff_input(content);
    let result = parse_feff_input(&file_path);

    // Should fail validation because iself is not in [0, 1, 2]
    assert!(result.is_err());
    match result {
        Err(InputError::InvalidFormat(msg)) => {
            assert!(msg.contains("Invalid iself value"));
        }
        _ => panic!("Expected InvalidFormat error for invalid iself value"),
    }
}

#[test]
fn test_write_screen() {
    // Create a minimal valid input file
    let content = r#"TITLE Test SCREEN write test

ATOMS
0 0.0 0.0 0.0 Fe

POTENTIALS
0 26 Fe
"#;

    let (_dir, file_path) = create_test_feff_input(content);
    let result = parse_feff_input(&file_path);
    assert!(result.is_ok());

    // Get the input and add SCREEN parameters
    let mut input = result.unwrap();

    // Add SCREEN parameters
    let screen_params = ScreenParams {
        iself: 2,
        iscreen: 0,
        ca1: 0.8,
        ci1: 1.2,
    };
    input.screen = Some(screen_params);

    // Create a temporary directory for the output
    let dir = tempdir().unwrap();
    let output_path = dir.path().join("test_write_screen.inp");

    // Write to file
    assert!(input.write(&output_path).is_ok());

    // Read back and verify
    let read_result = parse_feff_input(&output_path);
    assert!(read_result.is_ok());

    let read_input = read_result.unwrap();
    assert!(read_input.screen.is_some());

    let read_screen = read_input.screen.unwrap();
    assert_eq!(read_screen.iself, 2);
    assert_eq!(read_screen.iscreen, 0);
    assert_eq!(read_screen.ca1, 0.8);
    assert_eq!(read_screen.ci1, 1.2);
}

#[test]
fn test_parse_spectral_basic() {
    let content = r#"TITLE Test SPECTRAL card

ATOMS
0 0.0 0.0 0.0 Fe

SPECTRAL 1 1 1 -30.0 40.0 0.2

POTENTIALS
0 26 Fe
"#;

    let (_dir, file_path) = create_test_feff_input(content);
    let result = parse_feff_input(&file_path);

    assert!(result.is_ok());
    let input = result.unwrap();

    assert!(input.spectral.is_some());
    let spectral = input.spectral.unwrap();

    // Verify parsed parameters
    assert_eq!(spectral.ispect, 1);
    assert_eq!(spectral.ispsharp, 1);
    assert_eq!(spectral.isprule, 1);
    assert_eq!(spectral.emin, -30.0);
    assert_eq!(spectral.emax, 40.0);
    assert_eq!(spectral.estep, 0.2);
}

#[test]
fn test_parse_spectral_with_defaults() {
    let content = r#"TITLE Test SPECTRAL card with defaults

ATOMS
0 0.0 0.0 0.0 Fe

SPECTRAL

POTENTIALS
0 26 Fe
"#;

    let (_dir, file_path) = create_test_feff_input(content);
    let result = parse_feff_input(&file_path);

    assert!(result.is_ok());
    let input = result.unwrap();

    assert!(input.spectral.is_some());
    let spectral = input.spectral.unwrap();

    // Verify default parameters
    assert_eq!(spectral.ispect, 1); // Default enabled
    assert_eq!(spectral.ispsharp, 0); // Default no sharpening
    assert_eq!(spectral.isprule, 0);
    assert_eq!(spectral.emin, -20.0);
    assert_eq!(spectral.emax, 20.0);
    assert_eq!(spectral.estep, 0.1);
}

#[test]
fn test_parse_spectral_invalid_params() {
    let content = r#"TITLE Test SPECTRAL card with invalid parameters

ATOMS
0 0.0 0.0 0.0 Fe

SPECTRAL 2 0 0 -20.0 20.0 0.1

POTENTIALS
0 26 Fe
"#;

    let (_dir, file_path) = create_test_feff_input(content);
    let result = parse_feff_input(&file_path);

    // Should fail validation because ispect is not 0 or 1
    assert!(result.is_err());
    match result {
        Err(InputError::InvalidFormat(msg)) => {
            assert!(msg.contains("Invalid ispect value"));
        }
        _ => panic!("Expected InvalidFormat error for invalid ispect value"),
    }
}

#[test]
fn test_write_spectral() {
    // Create a minimal valid input file
    let content = r#"TITLE Test SPECTRAL write test

ATOMS
0 0.0 0.0 0.0 Fe

POTENTIALS
0 26 Fe
"#;

    let (_dir, file_path) = create_test_feff_input(content);
    let result = parse_feff_input(&file_path);
    assert!(result.is_ok());

    // Get the input and add SPECTRAL parameters
    let mut input = result.unwrap();

    // Add SPECTRAL parameters
    let spectral_params = SpectralParams {
        ispect: 1,
        ispsharp: 1,
        isprule: 1,
        emin: -25.0,
        emax: 35.0,
        estep: 0.25,
    };
    input.spectral = Some(spectral_params);

    // Create a temporary directory for the output
    let dir = tempdir().unwrap();
    let output_path = dir.path().join("test_write_spectral.inp");

    // Write to file
    assert!(input.write(&output_path).is_ok());

    // Read back and verify
    let read_result = parse_feff_input(&output_path);
    assert!(read_result.is_ok());

    let read_input = read_result.unwrap();
    assert!(read_input.spectral.is_some());

    let read_spectral = read_input.spectral.unwrap();
    assert_eq!(read_spectral.ispect, 1);
    assert_eq!(read_spectral.ispsharp, 1);
    assert_eq!(read_spectral.isprule, 1);
    assert_eq!(read_spectral.emin, -25.0);
    assert_eq!(read_spectral.emax, 35.0);
    assert_eq!(read_spectral.estep, 0.25);
}

#[test]
fn test_parse_dimensions_basic() {
    let content = r#"TITLE Test DIMENSIONS card

ATOMS
0 0.0 0.0 0.0 Fe

DIMENSIONS 20 250 300 10 600

POTENTIALS
0 26 Fe
"#;

    let (_dir, file_path) = create_test_feff_input(content);
    let result = parse_feff_input(&file_path);

    assert!(result.is_ok());
    let input = result.unwrap();

    assert!(input.dimensions.is_some());
    let dimensions = input.dimensions.unwrap();

    // Verify parsed parameters
    assert_eq!(dimensions.nat, 20);
    assert_eq!(dimensions.nph, 250);
    assert_eq!(dimensions.lx, 300);
    assert_eq!(dimensions.npot, 10);
    assert_eq!(dimensions.nstat, 600);
}

#[test]
fn test_parse_dimensions_with_defaults() {
    let content = r#"TITLE Test DIMENSIONS card with defaults

ATOMS
0 0.0 0.0 0.0 Fe

DIMENSIONS

POTENTIALS
0 26 Fe
"#;

    let (_dir, file_path) = create_test_feff_input(content);
    let result = parse_feff_input(&file_path);

    assert!(result.is_ok());
    let input = result.unwrap();

    assert!(input.dimensions.is_some());
    let dimensions = input.dimensions.unwrap();

    // Verify default parameters
    assert_eq!(dimensions.nat, 15);
    assert_eq!(dimensions.nph, 200);
    assert_eq!(dimensions.lx, 200);
    assert_eq!(dimensions.npot, 8);
    assert_eq!(dimensions.nstat, 500);
}

#[test]
fn test_parse_dimensions_invalid_params() {
    let content = r#"TITLE Test DIMENSIONS card with invalid parameters

ATOMS
0 0.0 0.0 0.0 Fe

DIMENSIONS 20 -10 300 10 600

POTENTIALS
0 26 Fe
"#;

    let (_dir, file_path) = create_test_feff_input(content);
    let result = parse_feff_input(&file_path);

    // Should fail validation because nph is negative
    assert!(result.is_err());
    match result {
        Err(InputError::InvalidFormat(msg)) => {
            assert!(msg.contains("Invalid nph value"));
        }
        _ => panic!("Expected InvalidFormat error for invalid nph value"),
    }
}

#[test]
fn test_write_dimensions() {
    // Create a minimal valid input file
    let content = r#"TITLE Test DIMENSIONS write test

ATOMS
0 0.0 0.0 0.0 Fe

POTENTIALS
0 26 Fe
"#;

    let (_dir, file_path) = create_test_feff_input(content);
    let result = parse_feff_input(&file_path);
    assert!(result.is_ok());

    // Get the input and add DIMENSIONS parameters
    let mut input = result.unwrap();

    // Add DIMENSIONS parameters
    let dimensions_params = DimensionsParams {
        nat: 18,
        nph: 220,
        lx: 250,
        npot: 12,
        nstat: 800,
    };
    input.dimensions = Some(dimensions_params);

    // Create a temporary directory for the output
    let dir = tempdir().unwrap();
    let output_path = dir.path().join("test_write_dimensions.inp");

    // Write to file
    assert!(input.write(&output_path).is_ok());

    // Read back and verify
    let read_result = parse_feff_input(&output_path);
    assert!(read_result.is_ok());

    let read_input = read_result.unwrap();
    assert!(read_input.dimensions.is_some());

    let read_dimensions = read_input.dimensions.unwrap();
    assert_eq!(read_dimensions.nat, 18);
    assert_eq!(read_dimensions.nph, 220);
    assert_eq!(read_dimensions.lx, 250);
    assert_eq!(read_dimensions.npot, 12);
    assert_eq!(read_dimensions.nstat, 800);
}

#[test]
fn test_parse_combined_cards() {
    let content = r#"TITLE Test combined SCREEN, SPECTRAL, and DIMENSIONS cards

ATOMS
0 0.0 0.0 0.0 Fe

SCREEN 2 0 0.8 1.2
SPECTRAL 1 1 0 -25.0 35.0 0.25
DIMENSIONS 18 220 250 12 800

POTENTIALS
0 26 Fe
"#;

    let (_dir, file_path) = create_test_feff_input(content);
    let result = parse_feff_input(&file_path);

    assert!(result.is_ok());
    let input = result.unwrap();

    // Verify all cards are parsed correctly
    assert!(input.screen.is_some());
    assert!(input.spectral.is_some());
    assert!(input.dimensions.is_some());

    let screen = input.screen.unwrap();
    let spectral = input.spectral.unwrap();
    let dimensions = input.dimensions.unwrap();

    // Verify SCREEN parameters
    assert_eq!(screen.iself, 2);
    assert_eq!(screen.iscreen, 0);
    assert_eq!(screen.ca1, 0.8);
    assert_eq!(screen.ci1, 1.2);

    // Verify SPECTRAL parameters
    assert_eq!(spectral.ispect, 1);
    assert_eq!(spectral.ispsharp, 1);
    assert_eq!(spectral.isprule, 0);
    assert_eq!(spectral.emin, -25.0);
    assert_eq!(spectral.emax, 35.0);
    assert_eq!(spectral.estep, 0.25);

    // Verify DIMENSIONS parameters
    assert_eq!(dimensions.nat, 18);
    assert_eq!(dimensions.nph, 220);
    assert_eq!(dimensions.lx, 250);
    assert_eq!(dimensions.npot, 12);
    assert_eq!(dimensions.nstat, 800);
}
