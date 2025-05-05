use feff_rs::input::{parse_feff_input, FeffInputParser, InputError, ParserConfig};
use std::fs::File;
use std::io::Write;
use std::path::Path;
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
fn test_parse_empty_input() {
    let (_dir, file_path) = create_test_feff_input("");
    let result = parse_feff_input(&file_path);
    assert!(result.is_err());
}

#[test]
fn test_parse_minimal_input() {
    let content = r#"TITLE Test minimal input

ATOMS
0 0.0 0.0 0.0 Fe
1 0.0 0.0 2.0 O

POTENTIALS
0 26 Fe
1 8 O
"#;
    let (_dir, file_path) = create_test_feff_input(content);
    let result = parse_feff_input(&file_path);
    assert!(result.is_ok());

    let input = result.unwrap();
    assert_eq!(input.title, Some("Test minimal input".to_string()));
    assert!(input.atomic_structure.is_some());

    let structure = input.atomic_structure.unwrap();
    assert_eq!(structure.atom_count(), 2); // Two atoms total (Fe + O)

    assert_eq!(input.potentials.len(), 2);
    assert_eq!(input.potentials.get(&0).unwrap().atomic_number, 26);
    assert_eq!(input.potentials.get(&0).unwrap().symbol, "Fe".to_string());
    assert_eq!(input.potentials.get(&1).unwrap().atomic_number, 8);
    assert_eq!(input.potentials.get(&1).unwrap().symbol, "O".to_string());
}

#[test]
fn test_parse_with_different_coord_systems() {
    let content = r#"TITLE Test coordinates

ATOMS
* Cartesian coordinates
0 0.0 0.0 0.0 Fe
1 1.0 1.0 1.0 O

ATOMS
* SPHERICAL coordinates
0 0.0 0.0 0.0 Fe
1 1.732 45.0 45.0 O

ATOMS
* CYLINDRICAL coordinates
0 0.0 0.0 0.0 Fe
1 1.414 45.0 1.0 O

POTENTIALS
0 26 Fe
1 8 O
"#;

    let (_dir, file_path) = create_test_feff_input(content);
    let result = parse_feff_input(&file_path);
    assert!(result.is_ok());

    // We should verify the coordinates are all correctly converted
    // but for now we'll just check the parser doesn't error
}

#[test]
fn test_parse_control_card() {
    let content = r#"TITLE Test control card

ATOMS
0 0.0 0.0 0.0 Fe
1 0.0 0.0 2.0 O

CONTROL 1 1 1 1 1 1

POTENTIALS
0 26 Fe
1 8 O
"#;

    let (_dir, file_path) = create_test_feff_input(content);
    let result = parse_feff_input(&file_path);
    assert!(result.is_ok());

    let input = result.unwrap();
    assert!(input.control.is_some());
    let control = input.control.unwrap();
    assert_eq!(control.mpot, 1);
    assert_eq!(control.mphase, 1);
    assert_eq!(control.mfms, 1);
    assert_eq!(control.mpath, 1);
    assert_eq!(control.mfeff, 1);
    assert_eq!(control.mchi, 1);
}

#[test]
fn test_parse_exchange_card() {
    let content = r#"TITLE Test exchange card

ATOMS
0 0.0 0.0 0.0 Fe
1 0.0 0.0 2.0 O

EXCHANGE 0 0.0 0.0 2

POTENTIALS
0 26 Fe
1 8 O
"#;

    let (_dir, file_path) = create_test_feff_input(content);
    let result = parse_feff_input(&file_path);
    assert!(result.is_ok());

    let input = result.unwrap();
    assert!(input.exchange.is_some());
    let exchange = input.exchange.unwrap();
    assert_eq!(exchange.ixc, 0);
    assert_eq!(exchange.vr0, 0.0);
    assert_eq!(exchange.vi0, 0.0);
    assert_eq!(exchange.ixc0, 2);
}

#[test]
fn test_parse_fms_card() {
    let content = r#"TITLE Test FMS card

ATOMS
0 0.0 0.0 0.0 Fe
1 0.0 0.0 2.0 O

FMS 6.0 0

POTENTIALS
0 26 Fe
1 8 O
"#;

    let (_dir, file_path) = create_test_feff_input(content);
    let result = parse_feff_input(&file_path);
    assert!(result.is_ok());

    let input = result.unwrap();
    assert!(input.fms.is_some());
    let fms = input.fms.unwrap();
    assert_eq!(fms.rfms, 6.0);
    assert_eq!(fms.nmultiple, 0);
}

#[test]
fn test_parse_complex_input() {
    let content = r#"TITLE Complex FEFF input test

CONTROL 1 1 1 1 1 1

EXCHANGE 0 0 0 0

SCF 4.0 0

COREHOLE RPA

XANES 8.0 0.07 0.0

FMS 4.0 0

POTENTIALS
0 26 Fe 3 3 0.01
1 8 O 3 3 1.0

ATOMS
* This is the absorbing atom
0 0.00000    0.00000    0.00000    Fe_abs
* These are the nearest neighbors
1 1.95000    1.95000    0.00000    O_1
1 1.95000   -1.95000    0.00000    O_2
1 0.00000    1.95000    1.95000    O_3
1 0.00000    1.95000   -1.95000    O_4
1 0.00000   -1.95000    1.95000    O_5
1 0.00000   -1.95000   -1.95000    O_6
"#;

    let (_dir, file_path) = create_test_feff_input(content);
    let result = parse_feff_input(&file_path);
    assert!(result.is_ok());

    let input = result.unwrap();
    assert_eq!(input.title, Some("Complex FEFF input test".to_string()));

    // Check atomic structure
    assert!(input.atomic_structure.is_some());
    let structure = input.atomic_structure.unwrap();
    assert_eq!(structure.atom_count(), 7); // 1 Fe + 6 O atoms

    // Check potentials
    assert_eq!(input.potentials.len(), 2);
    assert_eq!(input.potentials.get(&0).unwrap().atomic_number, 26);
    assert_eq!(input.potentials.get(&0).unwrap().symbol, "Fe".to_string());
    assert_eq!(input.potentials.get(&1).unwrap().atomic_number, 8);
    assert_eq!(input.potentials.get(&1).unwrap().symbol, "O".to_string());

    // Check other cards
    assert!(input.control.is_some());
    assert!(input.exchange.is_some());
    assert!(input.scf.is_some());
    assert!(input.xanes.is_some());
    assert!(input.fms.is_some());
}

#[test]
fn test_parser_config() {
    let content = r#"TITLE Test with custom config

ATOMS
0 0.0 0.0 0.0 Fe
1 0.0 0.0 2.0 O

POTENTIALS
0 26 Fe
1 8 O
"#;

    let (_dir, file_path) = create_test_feff_input(content);

    // Create a custom config
    let config = ParserConfig {
        input_path: file_path.clone(),
        validate: true,
        add_hydrogens: true,
        debug: true,
    };

    let mut parser = FeffInputParser::new(config);
    let result = parser.parse::<&Path>(None);
    assert!(result.is_ok());
}

#[test]
fn test_invalid_atoms_format() {
    let content = r#"TITLE Invalid atoms format test

ATOMS
0 not_a_number 0.0 0.0 Fe

POTENTIALS
0 26 Fe
"#;

    let (_dir, file_path) = create_test_feff_input(content);
    let result = parse_feff_input(&file_path);
    assert!(result.is_err());

    if let Err(err) = result {
        match err {
            InputError::ParseError(_) => {} // Expected error type
            _ => panic!("Expected ParseError"),
        }
    }
}

#[test]
fn test_invalid_potentials_format() {
    let content = r#"TITLE Invalid potentials format test

ATOMS
0 0.0 0.0 0.0 Fe

POTENTIALS
0 not_a_number Fe
"#;

    let (_dir, file_path) = create_test_feff_input(content);
    let result = parse_feff_input(&file_path);
    assert!(result.is_err());

    if let Err(err) = result {
        match err {
            InputError::ParseError(_) => {} // Expected error type
            _ => panic!("Expected ParseError"),
        }
    }
}

#[test]
fn test_missing_required_cards() {
    // Missing ATOMS card
    let content = r#"TITLE Missing atoms card test

POTENTIALS
0 26 Fe
1 8 O
"#;

    let (_dir, file_path) = create_test_feff_input(content);
    let result = parse_feff_input(&file_path);
    assert!(result.is_err());

    if let Err(err) = result {
        match err {
            InputError::MissingCard(card) => {
                assert_eq!(card, "ATOMS");
            }
            _ => panic!("Expected MissingCard error"),
        }
    }
}

#[test]
fn test_additional_cards() {
    let content = r#"TITLE Test additional cards

ATOMS
0 0.0 0.0 0.0 Fe
1 1.0 1.0 1.0 O

POTENTIALS
0 26 Fe
1 8 O

RPATH 0.05 10.0 4
PRINT 1
CORRECTIONS 0.01 0.05 1
S02 0.9
"#;

    let (_dir, file_path) = create_test_feff_input(content);
    let result = parse_feff_input(&file_path);
    assert!(result.is_ok());

    let input = result.unwrap();

    // Check RPATH card
    assert!(input.rpath.is_some());
    let rpath = input.rpath.unwrap();
    assert_eq!(rpath.pcrit, 0.05);
    assert_eq!(rpath.rmax, 10.0);
    assert_eq!(rpath.nleg, 4);

    // Check PRINT card
    assert!(input.print.is_some());
    assert_eq!(input.print.unwrap().iprint, 1);

    // Check CORRECTIONS card
    assert!(input.corrections.is_some());
    let corr = input.corrections.unwrap();
    assert_eq!(corr.real_correction, 0.01);
    assert_eq!(corr.imag_correction, 0.05);
    assert_eq!(corr.icorr, 1);

    // Check S02 card
    assert!(input.s02.is_some());
    assert_eq!(input.s02.unwrap().s02, 0.9);
}

#[test]
fn test_write_to_file() {
    let content = r#"TITLE Test write functionality

ATOMS
0 0.0 0.0 0.0 Fe
1 1.0 1.0 1.0 O

POTENTIALS
0 26 Fe
1 8 O

CONTROL 1 1 1 1 1 1
EXCHANGE 0 0.0 0.0 0
RPATH 0.05 10.0 4
PRINT 1
CORRECTIONS 0.01 0.05 1
S02 0.9
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

    // Verify the key components were preserved
    assert_eq!(read_input.title, input.title);
    assert!(read_input.atomic_structure.is_some());
    assert_eq!(read_input.potentials.len(), input.potentials.len());
    assert!(read_input.control.is_some());
    assert!(read_input.exchange.is_some());

    // Check that atomic structure was preserved
    let original_structure = input.atomic_structure.unwrap();
    let new_structure = read_input.atomic_structure.unwrap();
    assert_eq!(new_structure.atom_count(), original_structure.atom_count());
}
