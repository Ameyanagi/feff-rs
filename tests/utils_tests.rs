/*
MIT License with FEFF10 Attribution

Copyright (c) 2025 Ameyanagi

Based on or developed using Distribution: FEFF10.0
Copyright (c) 2020 FEFF Project, University of Washington and SLAC National Accelerator Laboratory.
All rights reserved.
*/

use approx::assert_relative_eq;
use feff_rs::utils::{angstrom_to_bohr, bohr_to_angstrom, constants, ev_to_hartree, hartree_to_ev};

#[test]
fn test_unit_conversions() {
    // Test Angstrom ↔ Bohr conversions
    let angstrom_value = 2.5;
    let bohr_value = angstrom_to_bohr(angstrom_value);
    let converted_back = bohr_to_angstrom(bohr_value);

    assert_relative_eq!(converted_back, angstrom_value, epsilon = 1e-10);
    assert_relative_eq!(
        bohr_value,
        angstrom_value / constants::BOHR_RADIUS,
        epsilon = 1e-10
    );

    // Test eV ↔ Hartree conversions
    let ev_value = 27.211;
    let hartree_value = ev_to_hartree(ev_value);
    let converted_back = hartree_to_ev(hartree_value);

    assert_relative_eq!(converted_back, ev_value, epsilon = 1e-10);
    assert_relative_eq!(
        hartree_value,
        ev_value / (2.0 * constants::RYDBERG),
        epsilon = 1e-10
    );
}
