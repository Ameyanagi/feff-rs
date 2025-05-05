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
fn test_edge_type_validation_valid() {
    // Test valid edge types
    let valid_types = ["K", "L1", "L2", "L3", "M1", "N4", "1", "4", "23"];

    for edge_type in valid_types.iter() {
        let content = format!(
            r#"TITLE Test EDGE validation with valid type: {}

ATOMS
0 0.0 0.0 0.0 Fe
1 0.0 0.0 2.0 O

EDGE {}

POTENTIALS
0 26 Fe
1 8 O
"#,
            edge_type, edge_type
        );

        let (_dir, file_path) = create_test_feff_input(&content);
        let result = parse_feff_input(&file_path);
        assert!(result.is_ok(), "Edge type '{}' should be valid", edge_type);

        let input = result.unwrap();
        assert!(input.edge.is_some());
        assert_eq!(input.edge.unwrap().edge_type, *edge_type);
    }
}

#[test]
fn test_edge_type_validation_invalid() {
    // Test invalid edge types
    let invalid_types = ["P", "X1", "L0", "M10", "0", "24", "abc"];

    for edge_type in invalid_types.iter() {
        let content = format!(
            r#"TITLE Test EDGE validation with invalid type: {}

ATOMS
0 0.0 0.0 0.0 Fe
1 0.0 0.0 2.0 O

EDGE {}

POTENTIALS
0 26 Fe
1 8 O
"#,
            edge_type, edge_type
        );

        let (_dir, file_path) = create_test_feff_input(&content);
        let result = parse_feff_input(&file_path);
        assert!(
            result.is_err(),
            "Edge type '{}' should be invalid",
            edge_type
        );

        match result {
            Err(err) => {
                let err_str = format!("{}", err);
                assert!(
                    err_str.contains("Invalid edge type"),
                    "Error message should mention invalid edge type, got: {}",
                    err_str
                );
            }
            _ => panic!("Expected error for invalid edge type '{}'", edge_type),
        }
    }
}

#[test]
fn test_edge_energy_validation() {
    // Test invalid energy values
    let content = r#"TITLE Test EDGE energy validation

ATOMS
0 0.0 0.0 0.0 Fe
1 0.0 0.0 2.0 O

EDGE K -100.0

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
            assert!(
                err_str.contains("Edge energy must be positive"),
                "Error message should mention positive energy, got: {}",
                err_str
            );
        }
        _ => panic!("Expected error for negative edge energy"),
    }
}

#[test]
fn test_polarization_vector_validation() {
    // Test invalid polarization vector (zero vector)
    let content = r#"TITLE Test POLARIZATION vector validation

ATOMS
0 0.0 0.0 0.0 Fe
1 0.0 0.0 2.0 O

POLARIZATION 0.0 0.0 0.0

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
            assert!(
                err_str.contains("Polarization vector cannot be a zero vector"),
                "Error message should mention zero vector, got: {}",
                err_str
            );
        }
        _ => panic!("Expected error for zero polarization vector"),
    }
}

#[test]
fn test_polarization_ellipticity_validation() {
    // Test invalid ellipticity value (outside -1 to 1 range)
    let content = r#"TITLE Test POLARIZATION ellipticity validation

ATOMS
0 0.0 0.0 0.0 Fe
1 0.0 0.0 2.0 O

POLARIZATION 1.0 0.0 0.0 1.5

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
            assert!(
                err_str.contains("Ellipticity should be between -1 and 1"),
                "Error message should mention ellipticity range, got: {}",
                err_str
            );
        }
        _ => panic!("Expected error for invalid ellipticity value"),
    }
}

#[test]
fn test_polarization_normalization() {
    // Test that the polarization vector gets normalized
    let content = r#"TITLE Test POLARIZATION normalization

ATOMS
0 0.0 0.0 0.0 Fe
1 0.0 0.0 2.0 O

POLARIZATION 3.0 4.0 0.0

POTENTIALS
0 26 Fe
1 8 O
"#;

    let (_dir, file_path) = create_test_feff_input(content);
    let result = parse_feff_input(&file_path);
    assert!(result.is_ok());

    let input = result.unwrap();
    assert!(input.polarization.is_some());
    let pol = input.polarization.unwrap();

    // Check that the vector was normalized (3-4-5 triangle)
    assert!((pol.x - 0.6).abs() < 1e-6);
    assert!((pol.y - 0.8).abs() < 1e-6);
    assert!((pol.z - 0.0).abs() < 1e-6);

    // Verify the length is 1.0
    let length = (pol.x * pol.x + pol.y * pol.y + pol.z * pol.z).sqrt();
    assert!((length - 1.0).abs() < 1e-6);
}
