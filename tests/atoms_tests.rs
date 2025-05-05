/*
MIT License with FEFF10 Attribution

Copyright (c) 2025 Ameyanagi

Based on or developed using Distribution: FEFF10.0
Copyright (c) 2020 FEFF Project, University of Washington and SLAC National Accelerator Laboratory.
All rights reserved.
*/

use feff_rs::atoms::{Atom, Vector3D};

#[test]
fn test_atom_creation() {
    // Test valid atom creation
    let position = Vector3D::new(0.0, 0.0, 0.0);
    let atom = Atom::new(29, position, 0).unwrap(); // Copper
    assert_eq!(atom.atomic_number(), 29);
    assert_eq!(atom.symbol(), "Cu");
}

#[test]
fn test_invalid_atom_creation() {
    let position = Vector3D::new(0.0, 0.0, 0.0);

    // Test invalid atomic numbers
    assert!(Atom::new(0, position, 0).is_err());
    assert!(Atom::new(119, position, 0).is_err());
    assert!(Atom::new(-1, position, 0).is_err());

    // Test invalid potential type
    assert!(Atom::new(29, position, -1).is_err());
}
