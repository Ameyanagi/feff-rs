/*
MIT License with FEFF10 Attribution

Copyright (c) 2025 Ameyanagi

Based on or developed using Distribution: FEFF10.0
Copyright (c) 2020 FEFF Project, University of Washington and SLAC National Accelerator Laboratory.
All rights reserved.
*/

//! Example demonstrating the use of thermal models in XAS calculations
//!
//! This example shows how to create and use different thermal models for EXAFS calculations.

use feff_rs::atoms::{Atom, AtomicStructure, PotentialType, Vector3D};
use feff_rs::utils::thermal::ThermalModel;
use feff_rs::xas::thermal::{create_anisotropic_thermal_parameters, ThermalParameters};
use std::fs::File;
use std::io::Write;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("Thermal effects example");

    // Create a Cu structure
    let structure = create_cu_structure();

    // Print number of atoms in the structure
    println!("Structure has {} atoms", structure.atom_count());

    // Demonstrate calculation of thermal parameters
    demo_thermal_parameters()?;

    println!("All calculations completed successfully!");
    Ok(())
}

/// Demonstrate thermal parameter calculations
fn demo_thermal_parameters() -> Result<(), Box<dyn std::error::Error>> {
    println!("Demonstrating thermal parameter calculations...");

    // Define parameters
    let debye_temp = 315.0; // Debye temperature for Cu in Kelvin
    let einstein_freq = 25.0; // Einstein frequency for Cu in meV
    let reduced_mass = 63.546; // Cu atomic mass in amu
    let path_length = 2.5; // Typical Cu-Cu bond length in Å

    // Create thermal parameters for different models
    let debye_params = ThermalParameters::new_debye(300.0, debye_temp);
    let einstein_params = ThermalParameters::new_einstein(300.0, einstein_freq);
    let corr_debye_params = ThermalParameters::new_correlated_debye(300.0, debye_temp);

    // Create the thermal models from the parameters
    let debye_model = debye_params.create_model(reduced_mass, Some(path_length));
    let einstein_model = einstein_params.create_model(reduced_mass, Some(path_length));
    let corr_debye_model = corr_debye_params.create_model(reduced_mass, Some(path_length));

    // Calculate Mean-Square Displacements (MSDs) at different temperatures
    let mut file = File::create("thermal_models_msd.dat")?;
    writeln!(
        file,
        "# Temperature(K) Debye_MSD Einstein_MSD Correlated_Debye_MSD"
    )?;

    // Calculate MSD at different temperatures from 50K to 900K
    for temp in (50..=900).step_by(50) {
        let debye_msd = debye_model.mean_square_displacement(temp as f64);
        let einstein_msd = einstein_model.mean_square_displacement(temp as f64);
        let corr_debye_msd = corr_debye_model.mean_square_displacement(temp as f64);

        writeln!(
            file,
            "{} {:.6} {:.6} {:.6}",
            temp, debye_msd, einstein_msd, corr_debye_msd
        )?;
    }

    println!("Thermal model MSDs written to thermal_models_msd.dat");

    // Calculate Debye-Waller factors at different k values
    let mut file = File::create("debye_waller_factors.dat")?;
    writeln!(file, "# k(Å⁻¹) Debye_DW Einstein_DW Correlated_DW")?;

    // Calculate for k values from 2 to 15 Å⁻¹
    let temp = 300.0; // Room temperature
    for k in (20..=150).map(|k| k as f64 / 10.0) {
        let debye_dw = debye_model.debye_waller_factor(temp, k);
        let einstein_dw = einstein_model.debye_waller_factor(temp, k);
        let corr_debye_dw = corr_debye_model.debye_waller_factor(temp, k);

        writeln!(
            file,
            "{:.1} {:.6} {:.6} {:.6}",
            k, debye_dw, einstein_dw, corr_debye_dw
        )?;
    }

    println!("Debye-Waller factors written to debye_waller_factors.dat");

    // Demonstrate anisotropic thermal parameters for different crystal systems
    let temp = 300.0;

    // Demonstrate anisotropic thermal parameters for different crystal systems
    let cubic_params = create_anisotropic_thermal_parameters(temp, "cubic", debye_temp, None);
    let tetragonal_params =
        create_anisotropic_thermal_parameters(temp, "tetragonal", debye_temp, Some(1.5));
    let layered_params =
        create_anisotropic_thermal_parameters(temp, "layered", debye_temp, Some(2.5));

    // Create models for different directions and calculate anisotropic MSDs
    let mut file = File::create("anisotropic_thermal.dat")?;
    writeln!(file, "# Direction Cubic_DW Tetragonal_DW Layered_DW")?;

    // Calculate for different directions
    let directions = [
        ([1.0, 0.0, 0.0], "x-axis"),
        ([0.0, 1.0, 0.0], "y-axis"),
        ([0.0, 0.0, 1.0], "z-axis"),
        ([1.0, 1.0, 0.0], "xy-plane"),
        ([0.0, 1.0, 1.0], "yz-plane"),
        ([1.0, 0.0, 1.0], "xz-plane"),
        ([1.0, 1.0, 1.0], "xyz-diagonal"),
    ];

    // Use a specific k value to show anisotropic effects
    let k_value = 8.0; // Typical EXAFS k-value in Å⁻¹

    for (dir, name) in &directions {
        // Set up the models with these directions
        let mut cubic = cubic_params.clone();
        cubic.path_direction = Some(*dir);

        let mut tetragonal = tetragonal_params.clone();
        tetragonal.path_direction = Some(*dir);

        let mut layered = layered_params.clone();
        layered.path_direction = Some(*dir);

        // Create thermal models and calculate Debye-Waller factors
        let cubic_model = cubic.create_model(reduced_mass, Some(path_length));
        let tetragonal_model = tetragonal.create_model(reduced_mass, Some(path_length));
        let layered_model = layered.create_model(reduced_mass, Some(path_length));

        let cubic_dw = cubic_model.debye_waller_factor(temp, k_value);
        let tetragonal_dw = tetragonal_model.debye_waller_factor(temp, k_value);
        let layered_dw = layered_model.debye_waller_factor(temp, k_value);

        writeln!(
            file,
            "{} {:.6} {:.6} {:.6}",
            name, cubic_dw, tetragonal_dw, layered_dw
        )?;
    }

    println!("Anisotropic thermal parameters written to anisotropic_thermal.dat");

    Ok(())
}

/// Create a simple Cu FCC structure for testing
fn create_cu_structure() -> AtomicStructure {
    let mut structure = AtomicStructure::new();

    // Add potential types
    let cu_potential = PotentialType::new(0, 29).unwrap();
    structure.add_potential_type(cu_potential);

    // Add atoms in an FCC arrangement
    let a = 3.6140; // Cu lattice constant in Å

    // Add central atom
    let cu_central = Atom::new(29, Vector3D::new(0.0, 0.0, 0.0), 0).unwrap();
    let central_idx = structure.add_atom(cu_central);

    // Add FCC lattice points
    structure.add_atom(Atom::new(29, Vector3D::new(0.0, 0.0, a), 0).unwrap());
    structure.add_atom(Atom::new(29, Vector3D::new(0.0, a, 0.0), 0).unwrap());
    structure.add_atom(Atom::new(29, Vector3D::new(0.0, a, a), 0).unwrap());
    structure.add_atom(Atom::new(29, Vector3D::new(a, 0.0, 0.0), 0).unwrap());
    structure.add_atom(Atom::new(29, Vector3D::new(a, 0.0, a), 0).unwrap());
    structure.add_atom(Atom::new(29, Vector3D::new(a, a, 0.0), 0).unwrap());
    structure.add_atom(Atom::new(29, Vector3D::new(a, a, a), 0).unwrap());

    // Set central atom
    structure.set_central_atom(central_idx).unwrap();

    structure
}
