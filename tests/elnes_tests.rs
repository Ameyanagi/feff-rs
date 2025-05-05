use feff_rs::input::{parse_feff_input, InputError};
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
fn test_parse_elnes_card() {
    let content = r#"TITLE Test ELNES card

ATOMS
0 0.0 0.0 0.0 Fe
1 0.0 0.0 2.0 O

ELNES 4.0 0.07 0.0
300
0 1 0
2.4 0.0
5 3
0.0 0.0

POTENTIALS
0 26 Fe
1 8 O
"#;

    let (_dir, file_path) = create_test_feff_input(content);
    let result = parse_feff_input(&file_path);
    assert!(result.is_ok());

    let input = result.unwrap();
    assert!(input.elnes.is_some());
    let elnes = input.elnes.unwrap();

    // Check parsed parameters
    assert_eq!(elnes.xkmax, 4.0);
    assert_eq!(elnes.xkstep, 0.07);
    assert_eq!(elnes.vixan, 0.0);
    assert_eq!(elnes.beam_energy, 300.0);
    assert_eq!(elnes.aver, None);
    assert_eq!(elnes.cross, None);
    assert_eq!(elnes.relat, None);
    assert_eq!(elnes.beam_direction.x, 0.0);
    assert_eq!(elnes.beam_direction.y, 1.0);
    assert_eq!(elnes.beam_direction.z, 0.0);
    assert_eq!(elnes.beta, 2.4);
    assert_eq!(elnes.alpha, 0.0);
    assert_eq!(elnes.nr, 5);
    assert_eq!(elnes.na, 3);
    assert_eq!(elnes.detector_position.x, 0.0);
    assert_eq!(elnes.detector_position.y, 0.0);
}

#[test]
fn test_parse_elnes_card_with_optional_params() {
    let content = r#"TITLE Test ELNES card with optional parameters

ATOMS
0 0.0 0.0 0.0 Fe
1 0.0 0.0 2.0 O

ELNES 5.0 0.08 0.02
300 1.0 0.5 1.0
0 0 1
3.0 0.5
10 5
0.1 0.2

POTENTIALS
0 26 Fe
1 8 O
"#;

    let (_dir, file_path) = create_test_feff_input(content);
    let result = parse_feff_input(&file_path);
    assert!(result.is_ok());

    let input = result.unwrap();
    assert!(input.elnes.is_some());
    let elnes = input.elnes.unwrap();

    // Check parsed parameters including optional ones
    assert_eq!(elnes.xkmax, 5.0);
    assert_eq!(elnes.xkstep, 0.08);
    assert_eq!(elnes.vixan, 0.02);
    assert_eq!(elnes.beam_energy, 300.0);
    assert_eq!(elnes.aver, Some(1.0));
    assert_eq!(elnes.cross, Some(0.5));
    assert_eq!(elnes.relat, Some(1.0));
    assert_eq!(elnes.beam_direction.x, 0.0);
    assert_eq!(elnes.beam_direction.y, 0.0);
    assert_eq!(elnes.beam_direction.z, 1.0);
    assert_eq!(elnes.beta, 3.0);
    assert_eq!(elnes.alpha, 0.5);
    assert_eq!(elnes.nr, 10);
    assert_eq!(elnes.na, 5);
    assert_eq!(elnes.detector_position.x, 0.1);
    assert_eq!(elnes.detector_position.y, 0.2);
}

#[test]
fn test_write_elnes_card() {
    let content = r#"TITLE Test writing ELNES card

ATOMS
0 0.0 0.0 0.0 Fe
1 0.0 0.0 2.0 O

ELNES 4.0 0.07 0.0
300 1.0 0.5 1.0
0 1 0
2.4 0.0
5 3
0.0 0.0

POTENTIALS
0 26 Fe
1 8 O
"#;

    let (_dir, file_path) = create_test_feff_input(content);
    let result = parse_feff_input(&file_path);
    assert!(result.is_ok());

    let input = result.unwrap();

    // Create a new temporary file to write to
    let dir = tempdir().unwrap();
    let output_path = dir.path().join("output.inp");

    // Write the input to the file
    let write_result = input.write(&output_path);
    assert!(write_result.is_ok());

    // Read it back in to verify
    let read_back = parse_feff_input(&output_path);
    assert!(read_back.is_ok());

    let read_input = read_back.unwrap();

    // Verify ELNES card was preserved
    assert!(read_input.elnes.is_some());
    let elnes = read_input.elnes.unwrap();

    // Check key parameters
    assert_eq!(elnes.xkmax, 4.0);
    assert_eq!(elnes.xkstep, 0.07);
    assert_eq!(elnes.vixan, 0.0);
    assert_eq!(elnes.beam_energy, 300.0);
    assert_eq!(elnes.beam_direction.y, 1.0);
    assert_eq!(elnes.beta, 2.4);
    assert_eq!(elnes.nr, 5);
    assert_eq!(elnes.na, 3);
}

#[test]
fn test_elnes_validation_missing_lines() {
    let content = r#"TITLE Test ELNES card validation

ATOMS
0 0.0 0.0 0.0 Fe
1 0.0 0.0 2.0 O

ELNES 4.0 0.07 0.0
300
0 1 0

POTENTIALS
0 26 Fe
1 8 O
"#;

    let (_dir, file_path) = create_test_feff_input(content);
    let result = parse_feff_input(&file_path);
    assert!(result.is_err());

    match result {
        Err(InputError::ParseError(msg)) => {
            assert!(msg.contains("Insufficient content for ELNES card"));
        }
        _ => panic!("Expected ParseError for insufficient ELNES card content"),
    }
}
