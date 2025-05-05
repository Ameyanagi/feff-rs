/*
MIT License with FEFF10 Attribution

Copyright (c) 2025 Ameyanagi

Based on or developed using Distribution: FEFF10.0
Copyright (c) 2020 FEFF Project, University of Washington and SLAC National Accelerator Laboratory.
All rights reserved.
*/

use feff_rs::atoms::Atom;

#[test]
fn test_atom_creation() {
    // Test valid atom creation
    let atom = Atom::new(29).unwrap(); // Copper
    assert_eq!(atom.atomic_number(), 29);

    // Test atom symbol (placeholder implementation for now)
    assert_eq!(atom.symbol(), "Element29");
}

#[test]
fn test_invalid_atom_creation() {
    // Test invalid atomic numbers
    assert!(Atom::new(0).is_err());
    assert!(Atom::new(119).is_err());
    assert!(Atom::new(-1).is_err());
}
