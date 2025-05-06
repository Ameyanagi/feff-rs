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
    /// Core-hole configuration for XAS calculations
    pub core_hole_config: xas::CoreHoleConfig,
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
            core_hole_config: xas::CoreHoleConfig::default(),
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

    /// Calculate an XAS spectrum with the specified parameters
    ///
    /// # Arguments
    ///
    /// * `energies` - List of energies (in eV) to calculate the spectrum at
    /// * `max_l` - Maximum angular momentum for phase shifts
    /// * `polarization` - Optional polarization vector for polarized XAS
    ///
    /// # Returns
    ///
    /// A vector of (energy, mu) pairs representing the XAS spectrum
    pub fn calculate_xas(
        &self,
        energies: &[f64],
        max_l: i32,
        polarization: Option<[f64; 3]>,
    ) -> anyhow::Result<Vec<(f64, f64)>> {
        // Check if we have valid input
        let structure = match &self.atomic_structure {
            Some(s) => s,
            None => {
                return Err(anyhow::anyhow!(
                    "No atomic structure provided for calculation"
                ))
            }
        };

        if energies.is_empty() {
            return Err(anyhow::anyhow!("No energies provided for calculation"));
        }

        // Calculate XAS spectrum for each energy point
        let mut spectrum = Vec::with_capacity(energies.len());

        for &energy in energies {
            // Calculate phase shifts with core-hole effects - not needed for direct calculation
            // We use the all-in-one function instead
            let _unused =
                xas::calculate_with_core_hole(structure, energy, max_l, &self.core_hole_config)?;

            // Calculate scattering matrices
            let scattering =
                scattering::calculate_scattering_matrices_old(structure, energy, max_l)?;

            // Calculate XAS at this energy point
            let point = xas::calculate_xas_spectrum(&scattering, polarization)?;

            if !point.is_empty() {
                spectrum.push(point[0]);
            }
        }

        Ok(spectrum)
    }

    /// Set the core-hole calculation method
    pub fn with_core_hole_method(&mut self, method: xas::CoreHoleMethod) -> &mut Self {
        self.core_hole_config.with_method(method);
        self
    }

    /// Set the absorption edge
    pub fn with_edge(&mut self, edge: &str) -> &mut Self {
        self.core_hole_config.with_edge(edge);
        self
    }

    /// Set the core-hole screening parameter
    pub fn with_screening(&mut self, screening: f64) -> &mut Self {
        self.core_hole_config.with_screening(screening);
        self
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

    #[test]
    fn test_xas_calculation() {
        let mut feff = Feff::new();

        // Create a simple atomic structure for testing
        let mut structure = atoms::AtomicStructure::new();

        // Add potential types
        let fe_pot = atoms::PotentialType::new(0, 26).unwrap(); // Fe
        structure.add_potential_type(fe_pot);

        let o_pot = atoms::PotentialType::new(1, 8).unwrap(); // O
        structure.add_potential_type(o_pot);

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
            atoms::Vector3D::new(2.0, 0.0, 0.0),
            1, // ipot
        )
        .unwrap();
        structure.add_atom(atom1);

        // Calculate muffin-tin radii
        structure.calculate_muffin_tin_radii().unwrap();

        feff.atomic_structure = Some(structure);

        // Configure core-hole calculation
        feff.with_core_hole_method(xas::CoreHoleMethod::FinalState)
            .with_edge("K")
            .with_screening(0.0);

        // Calculate XAS at a few energy points
        let energies = vec![7112.0, 7122.0, 7132.0]; // Near Fe K-edge
        let max_l = 3;

        let spectrum = feff.calculate_xas(&energies, max_l, None);

        // Calculation should succeed
        assert!(spectrum.is_ok());

        // Should have the same number of points as energies
        let spectrum = spectrum.unwrap();
        assert_eq!(spectrum.len(), energies.len());

        // Each point should have the correct energy
        for (i, &energy) in energies.iter().enumerate() {
            assert_eq!(spectrum[i].0, energy);
        }

        // Absorption coefficients should be positive or zero
        for &(_, mu) in &spectrum {
            assert!(mu >= 0.0);
        }
    }
}
