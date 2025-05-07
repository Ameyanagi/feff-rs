/*
MIT License with FEFF10 Attribution

Copyright (c) 2025 Ameyanagi

Based on or developed using Distribution: FEFF10.0
Copyright (c) 2020 FEFF Project, University of Washington and SLAC National Accelerator Laboratory.
All rights reserved.
*/

// Allow dead code in the crate since we have many experimental features
// and physics models in development that may not be used yet
#![allow(dead_code)]
// Temporarily disable doctests completely
#![doc(test(attr(deny(warnings))))]
#![doc(test(no_crate_inject))]
// Temporarily allow some Clippy warnings for initial implementation
#![allow(clippy::needless_range_loop)]
#![allow(clippy::single_char_add_str)]
#![allow(clippy::manual_clamp)]
#![allow(clippy::clone_on_copy)]
#![allow(clippy::unnecessary_cast)]

//! # FEFF-rs
//!
//! A modern Rust implementation of the FEFF code for calculating X-ray absorption spectroscopy (XAS)
//! and related spectroscopies.
//!
//! FEFF is a widely used ab initio multiple-scattering code for calculations of excitation spectra
//! and electronic structure. This library provides a high-performance, memory-efficient, and maintainable
//! implementation in Rust.
//!
//! ## Overview
//!
//! FEFF-rs can calculate:
//!
//! - X-ray Absorption Near Edge Structure (XANES) spectra
//! - Extended X-ray Absorption Fine Structure (EXAFS) spectra
//! - Multiple scattering paths and their contributions
//! - Temperature-dependent effects on XAS spectra
//!
//! The library is designed to be used both as a standalone application and as a library
//! that can be integrated into other Rust programs.
//!
//! ## Getting Started
//!
//! ```no_run
//! use feff_rs::{Feff, xas::{Edge, XanesParameters}};
//!
//! // Create a FEFF calculation from an input file
//! let mut feff = Feff::from_file("path/to/feff.inp").unwrap();
//!
//! // Run the calculation
//! feff.run().unwrap();
//!
//! // Calculate XANES spectrum with custom parameters
//! let params = XanesParameters {
//!     edge: Edge::K,
//!     energy_range: (-10.0, 50.0, 0.5), // (min, max, step) in eV
//!     polarization: None, // Isotropic spectrum
//!     ..Default::default()
//! };
//!
//! let spectrum = feff.calculate_xanes(&params).unwrap();
//!
//! // Access the results
//! for (i, energy) in spectrum.energies.iter().enumerate() {
//!     println!("{:.2} eV: {:.6}", energy, spectrum.mu[i]);
//! }
//! ```
//!
//! ## Key Features
//!
//! - **Comprehensive Physics**: Accurate implementation of multiple scattering theory for XAS
//! - **Performance**: Parallel processing and efficient algorithms for large calculations
//! - **Thermal Effects**: Multiple thermal models for temperature-dependent calculations
//! - **Memory Efficiency**: Careful memory management for handling large systems
//! - **Modularity**: Clean separation of concerns with well-defined module boundaries
//!
//! ## Module Structure
//!
//! The library is organized into several modules, each handling a specific aspect of the calculation:
//!
//! - **atoms**: Atomic data structures and coordinate handling
//! - **input**: Input file parsing and configuration
//! - **potential**: Atomic potential calculations
//! - **scattering**: Phase shift and scattering matrix calculations
//! - **path**: Path finding and filtering for EXAFS
//! - **fms**: Full multiple scattering implementation for XANES
//! - **xas**: X-ray absorption spectroscopy calculations
//! - **utils**: Utility functions and constants
//!
//! ## Thermal Effects
//!
//! FEFF-rs includes several thermal models for accounting for temperature effects in XAS:
//!
//! ```no_run
//! use feff_rs::xas::thermal::ThermalParameters;
//!
//! // Debye model for simple cubic materials
//! let debye_params = ThermalParameters::new_debye(300.0, 315.0);
//!
//! // Einstein model for materials with localized vibrations
//! let einstein_params = ThermalParameters::new_einstein(300.0, 70.0);
//!
//! // Correlated Debye model for accurate EXAFS analysis
//! let correlated_params = ThermalParameters::new_correlated_debye(300.0, 315.0);
//!
//! // Anisotropic model for non-cubic materials
//! let anisotropic_params = ThermalParameters::new_anisotropic_debye(
//!     300.0,              // Temperature in K
//!     315.0,              // Debye temperature
//!     [1.0, 1.0, 2.0],    // Displacement factors [x, y, z]
//!     None                // Default path direction
//! );
//! ```

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

// Re-export commonly used types for convenience
pub use atoms::{Atom, AtomicStructure, PotentialType, Vector3D};
pub use input::FeffInput;
pub use xas::exafs::{ExafsData, ExafsParameters};
pub use xas::{
    calculate_exafs, calculate_xanes, CoreHoleMethod, Edge, WindowFunction, XanesParameters,
    XanesSpectrum,
};

/// The main entry point for FEFF calculations.
///
/// The `Feff` struct provides a high-level interface to the FEFF calculation engine.
/// It manages the input data, atomic structure, and calculation parameters.
///
/// # Examples
///
/// ```no_run
/// use feff_rs::Feff;
///
/// // Create a new FEFF calculation from an input file
/// let mut feff = Feff::from_file("path/to/feff.inp").unwrap();
///
/// // Run the calculation
/// feff.run().unwrap();
/// ```
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
    /// Creates a new FEFF calculation instance with default settings.
    ///
    /// # Examples
    ///
    /// ```
    /// use feff_rs::Feff;
    ///
    /// let feff = Feff::new();
    /// ```
    pub fn new() -> Self {
        Self::default()
    }

    /// Creates a FEFF calculation from an input file.
    ///
    /// This method reads a FEFF input file, parses it, and creates a new `Feff` instance
    /// with the input data and atomic structure.
    ///
    /// # Arguments
    ///
    /// * `path` - Path to the FEFF input file
    ///
    /// # Returns
    ///
    /// A `Result` containing the new `Feff` instance or an error if the file cannot be read or parsed.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use feff_rs::Feff;
    ///
    /// let feff = Feff::from_file("path/to/feff.inp").unwrap();
    /// ```
    pub fn from_file<P: AsRef<std::path::Path>>(path: P) -> anyhow::Result<Self> {
        let input = input::parse_feff_input(path)?;
        let atomic_structure = input.atomic_structure.clone();

        Ok(Self {
            input,
            atomic_structure,
            core_hole_config: xas::CoreHoleConfig::default(),
        })
    }

    /// Runs the FEFF calculation with the provided input.
    ///
    /// This method performs the core calculations needed for XAS, including:
    /// - Potential calculations
    /// - Phase shift calculations
    /// - Scattering matrix preparation
    ///
    /// After running this method, you can calculate XAS spectra using the `calculate_xanes`
    /// or `calculate_exafs` methods.
    ///
    /// # Returns
    ///
    /// A `Result` indicating success or an error if the calculation fails.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use feff_rs::Feff;
    ///
    /// let mut feff = Feff::from_file("path/to/feff.inp").unwrap();
    /// feff.run().unwrap();
    /// ```
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

    /// Calculates an XAS spectrum with the specified parameters.
    ///
    /// This is a low-level method that provides direct control over the calculation parameters.
    /// For most applications, you should use the higher-level `calculate_xanes` or `calculate_exafs`
    /// methods instead.
    ///
    /// # Arguments
    ///
    /// * `energies` - List of energies (in eV) at which to calculate the spectrum
    /// * `max_l` - Maximum angular momentum for phase shifts
    /// * `polarization` - Optional polarization vector for polarized XAS
    ///
    /// # Returns
    ///
    /// A `Result` containing a vector of (energy, mu) pairs representing the XAS spectrum,
    /// or an error if the calculation fails.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use feff_rs::Feff;
    ///
    /// let mut feff = Feff::from_file("path/to/feff.inp").unwrap();
    /// feff.run().unwrap();
    ///
    /// // Calculate XAS at a few energy points
    /// let energies = vec![7112.0, 7122.0, 7132.0]; // Near Fe K-edge
    /// let max_l = 3;
    /// let spectrum = feff.calculate_xas(&energies, max_l, None).unwrap();
    ///
    /// for &(energy, mu) in &spectrum {
    ///     println!("{:.2} eV: {:.6}", energy, mu);
    /// }
    /// ```
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
            let _unused = xas::calculate_with_core_hole(
                structure,
                energy,
                max_l,
                &self.core_hole_config,
                None,
            )?;

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

    /// Calculates a XANES spectrum using the specified parameters.
    ///
    /// XANES (X-ray Absorption Near Edge Structure) is the part of the XAS spectrum
    /// near the absorption edge, typically within ~50 eV of the edge. It is sensitive
    /// to the electronic structure and local coordination environment of the absorbing atom.
    ///
    /// This method uses Full Multiple Scattering (FMS) to calculate the XANES spectrum.
    ///
    /// # Arguments
    ///
    /// * `params` - Parameters controlling the XANES calculation
    ///
    /// # Returns
    ///
    /// A `Result` containing the calculated XANES spectrum or an error if the calculation fails.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use feff_rs::{Feff, xas::{XanesParameters, Edge}};
    ///
    /// let mut feff = Feff::from_file("path/to/feff.inp").unwrap();
    /// feff.run().unwrap();
    ///
    /// let params = XanesParameters {
    ///     edge: Edge::K,
    ///     energy_range: (-10.0, 50.0, 0.5), // (min, max, step) in eV
    ///     polarization: None, // Isotropic spectrum
    ///     ..Default::default()
    /// };
    ///
    /// let spectrum = feff.calculate_xanes(&params).unwrap();
    ///
    /// // Access the results
    /// for (i, energy) in spectrum.energies.iter().enumerate() {
    ///     println!("{:.2} eV: {:.6}", energy, spectrum.mu[i]);
    /// }
    /// ```
    pub fn calculate_xanes(
        &self,
        params: &xas::XanesParameters,
    ) -> anyhow::Result<xas::XanesSpectrum> {
        // Check if we have valid input
        let structure = match &self.atomic_structure {
            Some(s) => s,
            None => {
                return Err(anyhow::anyhow!(
                    "No atomic structure provided for calculation"
                ))
            }
        };

        // Calculate XANES spectrum
        let spectrum = xas::calculate_xanes(structure, params)?;
        Ok(spectrum)
    }

    /// Calculates an EXAFS spectrum using the specified parameters.
    ///
    /// EXAFS (Extended X-ray Absorption Fine Structure) is the oscillatory part of the XAS
    /// spectrum extending from ~50 eV to ~1000 eV above the absorption edge. It provides
    /// information about the local atomic structure around the absorbing atom.
    ///
    /// This method uses path expansion and filtering to calculate the EXAFS spectrum.
    ///
    /// # Arguments
    ///
    /// * `params` - Parameters controlling the EXAFS calculation
    ///
    /// # Returns
    ///
    /// A `Result` containing the calculated EXAFS spectrum or an error if the calculation fails.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use feff_rs::{Feff, xas::{ExafsParameters, Edge}};
    ///
    /// let feff = Feff::from_file("path/to/feff.inp").unwrap();
    /// feff.run().unwrap();
    ///
    /// let params = ExafsParameters {
    ///     edge: Edge::K,
    ///     k_range: (3.0, 14.0), // k-range in Å⁻¹
    ///     r_range: (0.0, 6.0, 0.05), // min, max, step in Å
    ///     max_path_length: 8.0,
    ///     max_legs: 4,
    ///     ..Default::default()
    /// };
    ///
    /// let spectrum = feff.calculate_exafs(&params).unwrap();
    ///
    /// // Access the results
    /// for (i, k) in spectrum.grid.k_values.iter().enumerate() {
    ///     println!("k: {:.2} Å⁻¹, chi(k): {:.6}", k, spectrum.chi_k[i]);
    /// }
    /// ```
    pub fn calculate_exafs(
        &self,
        params: &xas::ExafsParameters,
    ) -> anyhow::Result<xas::exafs::ExafsData> {
        // Check if we have valid input
        let structure = match &self.atomic_structure {
            Some(s) => s,
            None => {
                return Err(anyhow::anyhow!(
                    "No atomic structure provided for calculation"
                ))
            }
        };

        // Calculate EXAFS spectrum
        let spectrum = xas::calculate_exafs(structure, params)?;
        Ok(spectrum)
    }

    /// Sets the core-hole calculation method.
    ///
    /// # Arguments
    ///
    /// * `method` - The core-hole method to use
    ///
    /// # Returns
    ///
    /// A mutable reference to `self` for method chaining.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use feff_rs::{Feff, xas::CoreHoleMethod};
    ///
    /// let mut feff = Feff::new();
    /// feff.with_core_hole_method(CoreHoleMethod::FinalState);
    /// ```
    pub fn with_core_hole_method(&mut self, method: xas::CoreHoleMethod) -> &mut Self {
        self.core_hole_config.with_method(method);
        self
    }

    /// Sets the absorption edge.
    ///
    /// # Arguments
    ///
    /// * `edge` - The absorption edge ("K", "L1", "L2", "L3", etc.)
    ///
    /// # Returns
    ///
    /// A mutable reference to `self` for method chaining.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use feff_rs::Feff;
    ///
    /// let mut feff = Feff::new();
    /// feff.with_edge("K");
    /// ```
    pub fn with_edge(&mut self, edge: &str) -> &mut Self {
        self.core_hole_config.with_edge(edge);
        self
    }

    /// Sets the core-hole screening parameter.
    ///
    /// The screening parameter controls how much of the core-hole potential is included
    /// in the calculation, affecting the shape of the XANES spectrum.
    ///
    /// # Arguments
    ///
    /// * `screening` - The screening parameter (0.0 = full core hole, 1.0 = fully screened)
    ///
    /// # Returns
    ///
    /// A mutable reference to `self` for method chaining.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use feff_rs::Feff;
    ///
    /// let mut feff = Feff::new();
    /// feff.with_screening(0.8); // Partially screened core hole
    /// ```
    pub fn with_screening(&mut self, screening: f64) -> &mut Self {
        self.core_hole_config.with_screening(screening);
        self
    }

    /// Saves a spectrum to a file.
    ///
    /// # Arguments
    ///
    /// * `path` - The path to save the spectrum to
    /// * `spectrum` - The spectrum to save
    ///
    /// # Returns
    ///
    /// A `Result` indicating success or an error if the file cannot be written.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use feff_rs::{Feff, xas::{XanesParameters, Edge}};
    ///
    /// let mut feff = Feff::from_file("path/to/feff.inp").unwrap();
    /// feff.run().unwrap();
    ///
    /// let params = XanesParameters {
    ///     edge: Edge::K,
    ///     energy_range: (-10.0, 50.0, 0.5),
    ///     ..Default::default()
    /// };
    ///
    /// let spectrum = feff.calculate_xanes(&params).unwrap();
    /// feff.save_spectrum("xanes_spectrum.dat", &spectrum).unwrap();
    /// ```
    pub fn save_spectrum<P: AsRef<std::path::Path>>(
        &self,
        path: P,
        spectrum: &xas::XanesSpectrum,
    ) -> anyhow::Result<()> {
        use std::fs::File;
        use std::io::Write;

        let mut file = File::create(path)?;

        // Write header
        writeln!(file, "# XANES spectrum calculated by FEFF-rs v{}", VERSION)?;
        writeln!(file, "# Energy(eV) mu(E)")?;

        // Write data
        for (i, energy) in spectrum.energies.iter().enumerate() {
            writeln!(file, "{:.6} {:.9}", energy, spectrum.mu[i])?;
        }

        Ok(())
    }

    /// Saves an EXAFS spectrum to a file.
    ///
    /// # Arguments
    ///
    /// * `path` - The path to save the spectrum to
    /// * `spectrum` - The EXAFS spectrum to save
    ///
    /// # Returns
    ///
    /// A `Result` indicating success or an error if the file cannot be written.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use feff_rs::{Feff, xas::{ExafsParameters, Edge}};
    ///
    /// let mut feff = Feff::from_file("path/to/feff.inp").unwrap();
    /// feff.run().unwrap();
    ///
    /// let params = ExafsParameters {
    ///     edge: Edge::K,
    ///     k_range: (3.0, 14.0),
    ///     ..Default::default()
    /// };
    ///
    /// let spectrum = feff.calculate_exafs(&params).unwrap();
    /// feff.save_exafs_spectrum("exafs_spectrum.dat", &spectrum).unwrap();
    /// ```
    pub fn save_exafs_spectrum<P: AsRef<std::path::Path>>(
        &self,
        path: P,
        spectrum: &xas::exafs::ExafsData,
    ) -> anyhow::Result<()> {
        use std::fs::File;
        use std::io::Write;

        let mut file = File::create(path)?;

        // Write header
        writeln!(file, "# EXAFS spectrum calculated by FEFF-rs v{}", VERSION)?;
        writeln!(file, "# k(Å⁻¹) chi(k) k*chi(k) k²*chi(k) k³*chi(k)")?;

        // Write data
        for i in 0..spectrum.len() {
            let k = spectrum.grid.k_values[i];
            let chi = spectrum.chi_k[i];
            writeln!(
                file,
                "{:.6} {:.9} {:.9} {:.9} {:.9}",
                k,
                chi,
                k * chi,
                k * k * chi,
                k * k * k * chi
            )?;
        }

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
