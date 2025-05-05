use feff_rs::input::{parameters::OpConsParams, parse_feff_input};
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
fn test_parse_opcons_card() {
    let content = r#"TITLE Test OPCONS card

ATOMS
0 0.0 0.0 0.0 Fe
1 1.0 1.0 1.0 O

OPCONS

POTENTIALS
0 26 Fe
1 8 O
"#;

    let (_dir, file_path) = create_test_feff_input(content);
    let result = parse_feff_input(&file_path);

    assert!(result.is_ok());
    let input = result.unwrap();

    assert!(input.opcons.is_some());
    let opcons = input.opcons.unwrap();

    // Verify OPCONS is enabled
    assert!(opcons.enabled);
}

#[test]
fn test_write_opcons_card() {
    // Create a FeffInput with a minimal structure and OPCONS
    let content = r#"TITLE Test OPCONS write test

ATOMS
0 0.0 0.0 0.0 Fe
1 1.0 1.0 1.0 O

POTENTIALS
0 26 Fe
1 8 O
"#;

    let (_dir, file_path) = create_test_feff_input(content);
    let result = parse_feff_input(&file_path);
    assert!(result.is_ok());

    // Get the input and add OPCONS
    let mut input = result.unwrap();

    // Add OPCONS parameters
    let opcons_params = OpConsParams { enabled: true };
    input.opcons = Some(opcons_params);

    // Create a temporary directory for the output
    let dir = tempdir().unwrap();
    let output_path = dir.path().join("test_write.inp");

    // Write to file
    assert!(input.write(&output_path).is_ok());

    // Read back and verify
    let read_result = parse_feff_input(&output_path);
    assert!(read_result.is_ok());

    let read_input = read_result.unwrap();
    assert!(read_input.opcons.is_some());

    let read_opcons = read_input.opcons.unwrap();
    assert!(read_opcons.enabled);
}

#[test]
fn test_opcons_with_other_cards() {
    let content = r#"TITLE Test OPCONS with other cards

ATOMS
0 0.0 0.0 0.0 Fe
1 1.0 1.0 1.0 O

OPCONS
ELLIPTICITY 0.5 1.0 0.0 0.0

POTENTIALS
0 26 Fe
1 8 O
"#;

    let (_dir, file_path) = create_test_feff_input(content);
    let result = parse_feff_input(&file_path);

    assert!(result.is_ok());
    let input = result.unwrap();

    // Verify both OPCONS and ELLIPTICITY are present
    assert!(input.opcons.is_some());
    assert!(input.ellipticity.is_some());

    let opcons = input.opcons.unwrap();
    let ellipticity = input.ellipticity.unwrap();

    assert!(opcons.enabled);
    assert_eq!(ellipticity.ellipticity, 0.5);
}
