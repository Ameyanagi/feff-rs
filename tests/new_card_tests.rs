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
fn test_parse_edge_card() {
    let content = r#"TITLE Test EDGE card

ATOMS
0 0.0 0.0 0.0 Fe
1 0.0 0.0 2.0 O

EDGE K 7112.0

POTENTIALS
0 26 Fe
1 8 O
"#;

    let (_dir, file_path) = create_test_feff_input(content);
    let result = parse_feff_input(&file_path);
    assert!(result.is_ok());

    let input = result.unwrap();
    assert!(input.edge.is_some());
    let edge = input.edge.unwrap();
    assert_eq!(edge.edge_type, "K");
    assert_eq!(edge.energy, Some(7112.0));
}

#[test]
fn test_parse_edge_card_no_energy() {
    let content = r#"TITLE Test EDGE card without energy

ATOMS
0 0.0 0.0 0.0 Fe
1 0.0 0.0 2.0 O

EDGE L3

POTENTIALS
0 26 Fe
1 8 O
"#;

    let (_dir, file_path) = create_test_feff_input(content);
    let result = parse_feff_input(&file_path);
    assert!(result.is_ok());

    let input = result.unwrap();
    assert!(input.edge.is_some());
    let edge = input.edge.unwrap();
    assert_eq!(edge.edge_type, "L3");
    assert_eq!(edge.energy, None);
}

#[test]
fn test_parse_debye_card() {
    let content = r#"TITLE Test DEBYE card

ATOMS
0 0.0 0.0 0.0 Fe
1 0.0 0.0 2.0 O

DEBYE 350.0 0.005 1

POTENTIALS
0 26 Fe
1 8 O
"#;

    let (_dir, file_path) = create_test_feff_input(content);
    let result = parse_feff_input(&file_path);
    assert!(result.is_ok());

    let input = result.unwrap();
    assert!(input.debye.is_some());
    let debye = input.debye.unwrap();
    assert_eq!(debye.temp, 350.0);
    assert_eq!(debye.debye_waller_factor, 0.005);
    assert_eq!(debye.correlated_debye, true);
}

#[test]
fn test_parse_ldos_card() {
    let content = r#"TITLE Test LDOS card

ATOMS
0 0.0 0.0 0.0 Fe
1 0.0 0.0 2.0 O

LDOS -25.0 15.0 0.2

POTENTIALS
0 26 Fe
1 8 O
"#;

    let (_dir, file_path) = create_test_feff_input(content);
    let result = parse_feff_input(&file_path);
    assert!(result.is_ok());

    let input = result.unwrap();
    assert!(input.ldos.is_some());
    let ldos = input.ldos.unwrap();
    assert_eq!(ldos.emin, -25.0);
    assert_eq!(ldos.emax, 15.0);
    assert_eq!(ldos.estep, 0.2);
}

#[test]
fn test_parse_exafs_card() {
    let content = r#"TITLE Test EXAFS card

ATOMS
0 0.0 0.0 0.0 Fe
1 0.0 0.0 2.0 O

EXAFS 10.0 800.0 0.5

POTENTIALS
0 26 Fe
1 8 O
"#;

    let (_dir, file_path) = create_test_feff_input(content);
    let result = parse_feff_input(&file_path);
    assert!(result.is_ok());

    let input = result.unwrap();
    assert!(input.exafs.is_some());
    let exafs = input.exafs.unwrap();
    assert_eq!(exafs.emin, 10.0);
    assert_eq!(exafs.emax, 800.0);
    assert_eq!(exafs.estep, 0.5);
}

#[test]
fn test_parse_danes_card() {
    let content = r#"TITLE Test DANES card

ATOMS
0 0.0 0.0 0.0 Fe
1 0.0 0.0 2.0 O

DANES 6.0 1.0 2.0 3.0

POTENTIALS
0 26 Fe
1 8 O
"#;

    let (_dir, file_path) = create_test_feff_input(content);
    let result = parse_feff_input(&file_path);
    assert!(result.is_ok());

    let input = result.unwrap();
    assert!(input.danes.is_some());
    let danes = input.danes.unwrap();
    assert_eq!(danes.radius, 6.0);
    assert_eq!(danes.parameters, vec![1.0, 2.0, 3.0]);
}

#[test]
fn test_parse_corehole_card() {
    let content = r#"TITLE Test COREHOLE card

ATOMS
0 0.0 0.0 0.0 Fe
1 0.0 0.0 2.0 O

COREHOLE RPA 1.0 0.5

POTENTIALS
0 26 Fe
1 8 O
"#;

    let (_dir, file_path) = create_test_feff_input(content);
    let result = parse_feff_input(&file_path);
    assert!(result.is_ok());

    let input = result.unwrap();
    assert!(input.corehole.is_some());
    let corehole = input.corehole.unwrap();
    assert_eq!(corehole.treatment, "RPA");
    assert_eq!(corehole.params, vec![1.0, 0.5]);
}

#[test]
fn test_parse_polarization_card() {
    let content = r#"TITLE Test POLARIZATION card

ATOMS
0 0.0 0.0 0.0 Fe
1 0.0 0.0 2.0 O

POLARIZATION 0.0 1.0 0.0 0.5

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
    assert_eq!(pol.x, 0.0);
    assert_eq!(pol.y, 1.0);
    assert_eq!(pol.z, 0.0);
    assert_eq!(pol.ellipticity, Some(0.5));
}

#[test]
fn test_parse_polarization_card_no_ellipticity() {
    let content = r#"TITLE Test POLARIZATION card without ellipticity

ATOMS
0 0.0 0.0 0.0 Fe
1 0.0 0.0 2.0 O

POLARIZATION 1.0 0.0 0.0

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
    assert_eq!(pol.x, 1.0);
    assert_eq!(pol.y, 0.0);
    assert_eq!(pol.z, 0.0);
    assert_eq!(pol.ellipticity, None);
}

#[test]
fn test_parse_real_card() {
    let content = r#"TITLE Test REAL card

ATOMS
0 0.0 0.0 0.0 Fe
1 0.0 0.0 2.0 O

REAL 0.03 120

POTENTIALS
0 26 Fe
1 8 O
"#;

    let (_dir, file_path) = create_test_feff_input(content);
    let result = parse_feff_input(&file_path);
    assert!(result.is_ok());

    let input = result.unwrap();
    assert!(input.real_grid.is_some());
    let real = input.real_grid.unwrap();
    assert_eq!(real.spacing, 0.03);
    assert_eq!(real.size, 120);
}

#[test]
fn test_parse_reciprocal_card() {
    let content = r#"TITLE Test RECIPROCAL card

ATOMS
0 0.0 0.0 0.0 Fe
1 0.0 0.0 2.0 O

RECIPROCAL 0.02 150

POTENTIALS
0 26 Fe
1 8 O
"#;

    let (_dir, file_path) = create_test_feff_input(content);
    let result = parse_feff_input(&file_path);
    assert!(result.is_ok());

    let input = result.unwrap();
    assert!(input.reciprocal_grid.is_some());
    let recip = input.reciprocal_grid.unwrap();
    assert_eq!(recip.spacing, 0.02);
    assert_eq!(recip.size, 150);
}

#[test]
fn test_write_new_cards() {
    let content = r#"TITLE Test writing all new cards

ATOMS
0 0.0 0.0 0.0 Fe
1 0.0 0.0 2.0 O

EDGE K 7112.0
DEBYE 350.0 0.005 1
LDOS -25.0 15.0 0.2
EXAFS 10.0 800.0 0.5
DANES 6.0 1.0 2.0 3.0
COREHOLE RPA 1.0 0.5
POLARIZATION 0.0 1.0 0.0 0.5
REAL 0.03 120
RECIPROCAL 0.02 150

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

    // Verify the key components were preserved
    assert_eq!(read_input.title, input.title);

    // Check that all new card types were preserved
    assert!(read_input.edge.is_some());
    assert!(read_input.debye.is_some());
    assert!(read_input.ldos.is_some());
    assert!(read_input.exafs.is_some());
    assert!(read_input.danes.is_some());
    assert!(read_input.corehole.is_some());
    assert!(read_input.polarization.is_some());
    assert!(read_input.real_grid.is_some());
    assert!(read_input.reciprocal_grid.is_some());

    // Check that specific values are preserved
    assert_eq!(read_input.edge.unwrap().edge_type, "K");
    assert_eq!(read_input.debye.unwrap().temp, 350.0);
    assert_eq!(read_input.ldos.unwrap().emin, -25.0);
    assert_eq!(read_input.exafs.unwrap().emax, 800.0);
    assert_eq!(read_input.danes.unwrap().radius, 6.0);
    assert_eq!(read_input.corehole.unwrap().treatment, "RPA");
    assert_eq!(read_input.polarization.unwrap().y, 1.0);
    assert_eq!(read_input.real_grid.unwrap().spacing, 0.03);
    assert_eq!(read_input.reciprocal_grid.unwrap().size, 150);
}
