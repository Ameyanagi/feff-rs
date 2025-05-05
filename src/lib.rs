/*
MIT License with FEFF10 Attribution

Copyright (c) 2025 Ameyanagi

Based on or developed using Distribution: FEFF10.0
Copyright (c) 2020 FEFF Project, University of Washington and SLAC National Accelerator Laboratory.
All rights reserved.
*/

//! # feff-rs
//!
//! A Rust implementation of the FEFF code for calculating X-ray absorption spectra.
//!
//! FEFF is an ab initio multiple-scattering code for calculating excitation spectra
//! and electronic structure. This crate provides a modern, high-performance
//! implementation in Rust.
//!
//! This project is currently under active development.

// Module declarations - these will be implemented incrementally
pub mod atoms;
pub mod cli;
pub mod fms;
pub mod input;
pub mod path;
pub mod potential;
pub mod scattering;
pub mod utils;
pub mod xas;

// Version information
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
pub const AUTHORS: &str = env!("CARGO_PKG_AUTHORS");

/// The main entry point for the FEFF calculation
#[derive(Default)]
pub struct Feff {
    /// The input data for the FEFF calculation
    pub input: input::FeffInput,
    /// The atomic structure for the calculation
    pub atomic_structure: Option<atoms::AtomicStructure>,
}

impl Feff {
    /// Create a new FEFF calculation instance
    pub fn new() -> Self {
        Self::default()
    }

    /// Create a FEFF calculation from input file
    pub fn from_file<P: AsRef<std::path::Path>>(path: P) -> anyhow::Result<Self> {
        let input = input::parse_feff_input(path)?;
        let atomic_structure = input.atomic_structure.clone();

        Ok(Self {
            input,
            atomic_structure,
        })
    }

    /// Run the calculation with the provided input
    pub fn run(&self) -> anyhow::Result<()> {
        // Will be implemented as development progresses
        println!("Running FEFF calculation...");

        // Check if we have valid input
        if self.atomic_structure.is_none() {
            return Err(anyhow::anyhow!(
                "No atomic structure provided for calculation"
            ));
        }

        // Further calculation steps will be implemented as development progresses

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs::File;
    use std::io::Write;
    use tempfile::tempdir;

    #[test]
    fn test_feff_creation() {
        let mut feff = Feff::new();

        // Create a simple atomic structure for testing
        let mut structure = atoms::AtomicStructure::new();

        // Add an absorbing atom
        let absorbing_atom = atoms::Atom::new(
            26, // Fe
            atoms::Vector3D::new(0.0, 0.0, 0.0),
            0, // ipot
        )
        .unwrap();

        // The central/absorbing atom is just the first atom added with index 0
        let absorbing_atom_index = structure.add_atom(absorbing_atom);
        structure.set_central_atom(absorbing_atom_index).unwrap();

        // Add some surrounding atoms
        let atom1 = atoms::Atom::new(
            8, // O
            atoms::Vector3D::new(1.0, 0.0, 0.0),
            1, // ipot
        )
        .unwrap();
        structure.add_atom(atom1);

        feff.atomic_structure = Some(structure);

        assert!(feff.run().is_ok());
    }

    #[test]
    fn test_feff_from_file() {
        // Create a temporary directory for testing
        let dir = tempdir().unwrap();
        let file_path = dir.path().join("feff.inp");

        // Write a test FEFF input file
        let mut file = File::create(&file_path).unwrap();
        writeln!(file, "TITLE Test FEFF input").unwrap();
        writeln!(file, "ATOMS").unwrap();
        writeln!(file, "0 0.0 0.0 0.0 Fe").unwrap();
        writeln!(file, "1 1.0 0.0 0.0 O").unwrap();
        writeln!(file, "POTENTIALS").unwrap();
        writeln!(file, "0 26 Fe").unwrap();
        writeln!(file, "1 8 O").unwrap();

        // Should fail because our Atom::new_with_position_and_potential
        // doesn't actually set the atomic number yet, just the potential type
        // This will be fixed in a future implementation
        let result = Feff::from_file(&file_path);
        assert!(result.is_ok());
    }
}
