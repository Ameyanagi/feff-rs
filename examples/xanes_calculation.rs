/*
MIT License with FEFF10 Attribution

Copyright (c) 2025 Ameyanagi

Based on or developed using Distribution: FEFF10.0
Copyright (c) 2020 FEFF Project, University of Washington and SLAC National Accelerator Laboratory.
All rights reserved.
*/

//! Example XANES calculation
//!
//! This example demonstrates how to calculate a XANES spectrum
//! for an iron oxide cluster using the FEFF-rs library.

use feff_rs::atoms::{Atom, AtomicStructure, PotentialType, Vector3D};
use feff_rs::xas::{calculate_xanes, CoreHoleMethod, Edge, XanesParameters};
use std::error::Error;

fn main() -> Result<(), Box<dyn Error>> {
    // Create a simple Fe-O octahedral cluster
    let mut structure = AtomicStructure::new();

    // Add potentials
    let fe_potential = PotentialType::new(0, 26)?; // Iron
    let o_potential = PotentialType::new(1, 8)?; // Oxygen

    structure.add_potential_type(fe_potential);
    structure.add_potential_type(o_potential);

    // Add central Fe atom
    let fe_atom = Atom::new(26, Vector3D::new(0.0, 0.0, 0.0), 0)?;
    let central_idx = structure.add_atom(fe_atom);
    structure.set_central_atom(central_idx)?;

    // Add 6 oxygen atoms in octahedral coordination
    let distance = 2.0; // Fe-O distance in Angstroms

    structure.add_atom(Atom::new(8, Vector3D::new(distance, 0.0, 0.0), 1)?);
    structure.add_atom(Atom::new(8, Vector3D::new(-distance, 0.0, 0.0), 1)?);
    structure.add_atom(Atom::new(8, Vector3D::new(0.0, distance, 0.0), 1)?);
    structure.add_atom(Atom::new(8, Vector3D::new(0.0, -distance, 0.0), 1)?);
    structure.add_atom(Atom::new(8, Vector3D::new(0.0, 0.0, distance), 1)?);
    structure.add_atom(Atom::new(8, Vector3D::new(0.0, 0.0, -distance), 1)?);

    // Calculate muffin-tin radii
    structure.calculate_muffin_tin_radii()?;

    println!(
        "Created Fe-O octahedral cluster with {} atoms",
        structure.atom_count()
    );

    // Set up XANES calculation parameters
    let params = XanesParameters {
        edge: Edge::K,                                // Calculate Fe K-edge
        energy_range: (-20.0, 50.0, 2.0),             // Energy range relative to edge
        fermi_energy: 0.0,                            // Fermi energy in eV
        energy_shift: 0.0,                            // Energy shift in eV
        polarization: None,                           // Isotropic spectrum
        core_hole_lifetime: None,                     // Use default based on element and edge
        gaussian_broadening: 1.0,                     // Additional Gaussian broadening (eV)
        lorentzian_broadening: 0.5,                   // Additional Lorentzian broadening (eV)
        energy_dependent_broadening: 0.1,             // Energy-dependent broadening factor
        thermal_parameters: None,                     // No thermal effects
        include_quadrupole: true,                     // Include quadrupole transitions
        max_l: 3,                                     // Maximum angular momentum
        core_hole_method: CoreHoleMethod::FinalState, // Use final state rule
        core_hole_screening: 0.0,                     // No core-hole screening
    };

    println!("Calculating Fe K-edge XANES spectrum...");

    // Perform XANES calculation
    let spectrum = calculate_xanes(&structure, &params)?;

    println!("XANES calculation complete!");
    println!("Edge energy: {:.2} eV", spectrum.edge_energy);
    println!("Energy points: {}", spectrum.energies.len());

    // Export spectrum to file
    let output_file = "fe_k_edge_xanes.dat";
    spectrum.export_to_file(output_file)?;

    println!("Spectrum saved to {}", output_file);

    // Find the white line using the analyzer
    let analyzer = feff_rs::xas::XanesAnalyzer::new(&spectrum);
    if let Some((energy, intensity)) = analyzer.find_white_line() {
        println!(
            "White line found at {:.2} eV with intensity {:.4}",
            energy, intensity
        );
    }

    Ok(())
}
