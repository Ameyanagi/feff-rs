use feff_rs::input::parse_feff_input;
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
fn test_parse_nrixs_card_spherical() {
    let content = r#"TITLE Test NRIXS card with spherical averaging

ATOMS
0 0.0 0.0 0.0 Fe
1 0.0 0.0 2.0 O

NRIXS -1 24.0

POTENTIALS
0 26 Fe
1 8 O
"#;

    let (_dir, file_path) = create_test_feff_input(content);
    let result = parse_feff_input(&file_path);
    assert!(result.is_ok());

    let input = result.unwrap();
    assert!(input.nrixs.is_some());
    let nrixs = input.nrixs.unwrap();

    // Check parsed parameters for spherical averaging
    assert_eq!(nrixs.nq, -1);
    assert_eq!(nrixs.qx, 24.0);
    assert_eq!(nrixs.qy, 0.0); // Should be 0 for spherical averaging
    assert_eq!(nrixs.qz, 0.0); // Should be 0 for spherical averaging
    assert_eq!(nrixs.scalar, None);
}

#[test]
fn test_parse_nrixs_card_specific_q() {
    let content = r#"TITLE Test NRIXS card with specific q-vector

ATOMS
0 0.0 0.0 0.0 Fe
1 0.0 0.0 2.0 O

NRIXS 1 5.13 0.0 0.0 1.0

POTENTIALS
0 26 Fe
1 8 O
"#;

    let (_dir, file_path) = create_test_feff_input(content);
    let result = parse_feff_input(&file_path);
    assert!(result.is_ok());

    let input = result.unwrap();
    assert!(input.nrixs.is_some());
    let nrixs = input.nrixs.unwrap();

    // Check parsed parameters for specific q-vector
    assert_eq!(nrixs.nq, 1);
    assert_eq!(nrixs.qx, 5.13);
    assert_eq!(nrixs.qy, 0.0);
    assert_eq!(nrixs.qz, 0.0);
    assert_eq!(nrixs.scalar, Some(1.0));
}

#[test]
fn test_parse_nrixs_card_validation() {
    let content = r#"TITLE Test NRIXS card validation with insufficient parameters

ATOMS
0 0.0 0.0 0.0 Fe
1 0.0 0.0 2.0 O

NRIXS 1

POTENTIALS
0 26 Fe
1 8 O
"#;

    let (_dir, file_path) = create_test_feff_input(content);
    let result = parse_feff_input(&file_path);
    assert!(result.is_err());

    match result {
        Err(err) => {
            let err_str = format!("{}", err);
            assert!(err_str.contains("Insufficient parameters for NRIXS card"));
        }
        _ => panic!("Expected error for insufficient NRIXS parameters"),
    }
}

#[test]
fn test_parse_nrixs_card_missing_q_components() {
    let content = r#"TITLE Test NRIXS card validation with missing q components

ATOMS
0 0.0 0.0 0.0 Fe
1 0.0 0.0 2.0 O

NRIXS 1 5.0

POTENTIALS
0 26 Fe
1 8 O
"#;

    let (_dir, file_path) = create_test_feff_input(content);
    let result = parse_feff_input(&file_path);
    assert!(result.is_err());

    match result {
        Err(err) => {
            let err_str = format!("{}", err);
            assert!(err_str.contains("Missing q-vector components (qy, qz)"));
        }
        _ => panic!("Expected error for missing q-vector components"),
    }
}

#[test]
fn test_write_nrixs_card() {
    let content = r#"TITLE Test writing NRIXS card

ATOMS
0 0.0 0.0 0.0 Fe
1 0.0 0.0 2.0 O

NRIXS -1 24.0

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

    // Verify NRIXS card was preserved
    assert!(read_input.nrixs.is_some());
    let nrixs = read_input.nrixs.unwrap();

    // Check key parameters
    assert_eq!(nrixs.nq, -1);
    assert_eq!(nrixs.qx, 24.0);
}
