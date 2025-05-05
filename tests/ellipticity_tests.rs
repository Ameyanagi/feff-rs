use feff_rs::input::{parameters::EllipticityParams, parse_feff_input, InputError};
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
fn test_parse_ellipticity_basic() {
    let content = r#"TITLE Test ELLIPTICITY card

ATOMS
0 0.0 0.0 0.0 Fe

ELLIPTICITY 0.5 1.0 0.0 0.0

POTENTIALS
0 26 Fe
"#;

    let (_dir, file_path) = create_test_feff_input(content);
    let result = parse_feff_input(&file_path);

    assert!(result.is_ok());
    let input = result.unwrap();

    assert!(input.ellipticity.is_some());
    let ellipticity = input.ellipticity.unwrap();

    assert_eq!(ellipticity.ellipticity, 0.5);
    assert_eq!(ellipticity.beam_direction.x, 1.0);
    assert_eq!(ellipticity.beam_direction.y, 0.0);
    assert_eq!(ellipticity.beam_direction.z, 0.0);
}

#[test]
fn test_parse_ellipticity_circular() {
    let content = r#"TITLE Test ELLIPTICITY card with circular polarization

ATOMS
0 0.0 0.0 0.0 Fe

ELLIPTICITY 1.0 0.0 0.0 1.0

POTENTIALS
0 26 Fe
"#;

    let (_dir, file_path) = create_test_feff_input(content);
    let result = parse_feff_input(&file_path);

    assert!(result.is_ok());
    let input = result.unwrap();

    assert!(input.ellipticity.is_some());
    let ellipticity = input.ellipticity.unwrap();

    assert_eq!(ellipticity.ellipticity, 1.0); // Circular polarization
    assert_eq!(ellipticity.beam_direction.x, 0.0);
    assert_eq!(ellipticity.beam_direction.y, 0.0);
    assert_eq!(ellipticity.beam_direction.z, 1.0);
}

#[test]
fn test_parse_ellipticity_with_polarization() {
    let content = r#"TITLE Test ELLIPTICITY card with POLARIZATION

ATOMS
0 0.0 0.0 0.0 Fe

POLARIZATION 1.0 0.0 0.0
ELLIPTICITY 0.5 0.0 1.0 0.0

POTENTIALS
0 26 Fe
"#;

    let (_dir, file_path) = create_test_feff_input(content);
    let result = parse_feff_input(&file_path);

    assert!(result.is_ok());
    let input = result.unwrap();

    assert!(input.ellipticity.is_some());
    assert!(input.polarization.is_some());

    let ellipticity = input.ellipticity.unwrap();
    let polarization = input.polarization.unwrap();

    // Verify polarization is along x-axis
    assert_eq!(polarization.x, 1.0);
    assert_eq!(polarization.y, 0.0);
    assert_eq!(polarization.z, 0.0);

    // Verify ellipticity params
    assert_eq!(ellipticity.ellipticity, 0.5);

    // Verify beam direction is along y-axis (perpendicular to polarization)
    assert_eq!(ellipticity.beam_direction.x, 0.0);
    assert_eq!(ellipticity.beam_direction.y, 1.0);
    assert_eq!(ellipticity.beam_direction.z, 0.0);
}

#[test]
fn test_parse_ellipticity_invalid_range() {
    let content = r#"TITLE Test ELLIPTICITY card with invalid ellipticity value

ATOMS
0 0.0 0.0 0.0 Fe

ELLIPTICITY 2.0 1.0 0.0 0.0

POTENTIALS
0 26 Fe
"#;

    let (_dir, file_path) = create_test_feff_input(content);
    let result = parse_feff_input(&file_path);

    // Should fail validation because ellipticity is not in [-1, 1]
    assert!(result.is_err());

    match result {
        Err(InputError::InvalidFormat(msg)) => {
            assert!(msg.contains("Ellipticity should be between -1 and 1"));
        }
        _ => panic!("Expected InvalidFormat error for ellipticity out of range"),
    }
}

#[test]
fn test_parse_ellipticity_zero_vector() {
    let content = r#"TITLE Test ELLIPTICITY card with zero beam direction vector

ATOMS
0 0.0 0.0 0.0 Fe

ELLIPTICITY 0.5 0.0 0.0 0.0

POTENTIALS
0 26 Fe
"#;

    let (_dir, file_path) = create_test_feff_input(content);
    let result = parse_feff_input(&file_path);

    // Should fail validation because beam direction is a zero vector
    assert!(result.is_err());

    match result {
        Err(InputError::InvalidFormat(msg)) => {
            assert!(msg.contains("Beam direction vector cannot be a zero vector"));
        }
        _ => panic!("Expected InvalidFormat error for zero beam direction vector"),
    }
}

#[test]
fn test_write_ellipticity() {
    // Instead of creating a minimal structure, let's use an existing valid input file
    let content = r#"TITLE Test ELLIPTICITY write test

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

    // Get the input and add ellipticity
    let mut input = result.unwrap();

    // Add ELLIPTICITY parameters
    let ellipticity_params = EllipticityParams {
        ellipticity: 0.75,
        beam_direction: feff_rs::atoms::Vector3D::new(0.0, 0.0, 1.0),
    };
    input.ellipticity = Some(ellipticity_params);

    // Create a temporary directory for the output
    let dir = tempdir().unwrap();
    let output_path = dir.path().join("test_write.inp");

    // Write to file
    assert!(input.write(&output_path).is_ok());

    // Read back and verify
    let read_result = parse_feff_input(&output_path);
    assert!(read_result.is_ok());

    let read_input = read_result.unwrap();
    assert!(read_input.ellipticity.is_some());

    let read_ellipticity = read_input.ellipticity.unwrap();
    assert_eq!(read_ellipticity.ellipticity, 0.75);
    assert_eq!(read_ellipticity.beam_direction.x, 0.0);
    assert_eq!(read_ellipticity.beam_direction.y, 0.0);
    assert_eq!(read_ellipticity.beam_direction.z, 1.0);
}
