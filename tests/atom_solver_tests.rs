/*
MIT License with FEFF10 Attribution

Copyright (c) 2025 Ameyanagi

Based on or developed using Distribution: FEFF10.0
Copyright (c) 2020 FEFF Project, University of Washington and SLAC National Accelerator Laboratory.
All rights reserved.
*/

//! Tests for the atomic solver implementation

use approx::assert_relative_eq;
use feff_rs::potential::{AtomSolver, AtomSolverConfig};
use std::f64::consts::PI;

/// Create a simple hydrogen-like potential for testing
fn create_test_potential(z: i32, n_points: usize, r_max: f64) -> (Vec<f64>, Vec<f64>) {
    let dr = r_max / (n_points as f64);
    let mut grid = Vec::with_capacity(n_points);
    let mut potential = Vec::with_capacity(n_points);

    for i in 0..n_points {
        let r = (i as f64 + 0.001) * dr; // Avoid r=0
        grid.push(r);

        // Simple Coulomb potential: -Z/r
        let v = if r < 1e-10 { -1000.0 } else { -(z as f64) / r };
        potential.push(v);
    }

    (grid, potential)
}

#[test]
fn test_hydrogen_ground_state() {
    // Create potential for hydrogen (Z=1)
    let (grid, potential) = create_test_potential(1, 200, 30.0);

    // Create solver configuration
    let config = AtomSolverConfig {
        atomic_number: 1,
        relativistic: false,
        energy_tolerance: 1e-6,
        max_iterations: 100,
        use_experimental: false,
    };

    // Create solver
    let mut solver = AtomSolver::new(config, grid, potential);

    // Solve for 1s state
    let wavefunction = solver.solve_state(1, 0).unwrap();

    // Check the energy (hydrogen 1s should be -0.5 Hartree)
    assert_relative_eq!(wavefunction.energy, -0.5, epsilon = 0.01);

    // Should have 0 nodes for ground state
    assert_eq!(wavefunction.nodes, 0);

    // Should be normalized
    assert!(wavefunction.normalized);

    // Calculate electron density with occupation 1.0
    let density = wavefunction.calculate_density(1.0).unwrap();

    // Integrate density to get total electrons
    let mut total_electrons = 0.0;
    for i in 1..density.len() {
        let r1 = wavefunction.grid[i - 1];
        let r2 = wavefunction.grid[i];
        let rho1 = density[i - 1];
        let rho2 = density[i];

        // Volume element
        let vol1 = 4.0 * PI * r1 * r1;
        let vol2 = 4.0 * PI * r2 * r2;

        // Integrate using trapezoidal rule
        let dr = r2 - r1;
        total_electrons += 0.5 * dr * (rho1 * vol1 + rho2 * vol2);
    }

    // With occupation 1.0 and l=0, we should have 2 electrons (2 for spin)
    assert_relative_eq!(total_electrons, 2.0, epsilon = 0.1);
}

#[test]
fn test_hydrogen_excited_states() {
    // Create potential for hydrogen (Z=1)
    let (grid, potential) = create_test_potential(1, 300, 40.0);

    // Create solver configuration
    let config = AtomSolverConfig {
        atomic_number: 1,
        relativistic: false,
        energy_tolerance: 1e-6,
        max_iterations: 100,
        use_experimental: false,
    };

    // Create solver
    let mut solver = AtomSolver::new(config, grid, potential);

    // Solve for 2s state
    let wf_2s = solver.solve_state(2, 0).unwrap();

    // Check the energy (hydrogen 2s should be -0.125 Hartree)
    assert_relative_eq!(wf_2s.energy, -0.125, epsilon = 0.01);

    // Should have 1 node for 2s state
    assert_eq!(wf_2s.nodes, 1);

    // Solve for 2p state
    let wf_2p = solver.solve_state(2, 1).unwrap();

    // Check the energy (should be the same as 2s for hydrogen)
    assert_relative_eq!(wf_2p.energy, -0.125, epsilon = 0.01);

    // Should have 0 nodes for 2p state
    assert_eq!(wf_2p.nodes, 0);

    // Solve for 3d state
    let wf_3d = solver.solve_state(3, 2).unwrap();

    // Check the energy (hydrogen 3d should be -0.056 Hartree)
    assert_relative_eq!(wf_3d.energy, -0.056, epsilon = 0.01);

    // Should have 0 nodes for 3d state
    assert_eq!(wf_3d.nodes, 0);
}

#[test]
fn test_helium_energies() {
    // Create potential for helium (Z=2)
    let (grid, potential) = create_test_potential(2, 200, 30.0);

    // Create solver configuration
    let config = AtomSolverConfig {
        atomic_number: 2,
        relativistic: false,
        energy_tolerance: 1e-6,
        max_iterations: 100,
        use_experimental: false,
    };

    // Create solver
    let mut solver = AtomSolver::new(config, grid, potential);

    // Solve for 1s state
    let wf_1s = solver.solve_state(1, 0).unwrap();

    // For pure Coulomb potential, the energies scale with Z²
    // So helium 1s should be approximately -2.0 Hartree
    // (In reality it's more complex due to electron-electron interactions)
    assert_relative_eq!(wf_1s.energy, -2.0, epsilon = 0.1);

    // Now let's add a core hole with Z*=1
    solver.set_core_hole(1.0).unwrap();

    // Solve for 1s state with core hole
    let wf_1s_hole = solver.solve_state(1, 0).unwrap();

    // Energy should be more negative with a core hole
    assert!(wf_1s_hole.energy < wf_1s.energy);

    // The energy should be closer to Z=3 (lithium)
    assert_relative_eq!(wf_1s_hole.energy, -4.5, epsilon = 0.5);

    // Clear the core hole
    solver.clear_core_hole();

    // Solve for 1s state again to verify we're back to original potential
    let wf_1s_again = solver.solve_state(1, 0).unwrap();
    assert_relative_eq!(wf_1s_again.energy, wf_1s.energy, epsilon = 0.001);
}

#[test]
fn test_electronic_configuration() {
    // Test for carbon (Z=6)
    let (grid, potential) = create_test_potential(6, 200, 30.0);

    let config = AtomSolverConfig {
        atomic_number: 6,
        relativistic: false,
        energy_tolerance: 1e-6,
        max_iterations: 100,
        use_experimental: false,
    };

    let solver = AtomSolver::new(config, grid, potential);

    // Get the electronic configuration
    let configuration = solver.electronic_configuration();

    // Carbon should have 1s², 2s², 2p²
    assert_eq!(configuration.len(), 3);

    // Check quantum numbers
    assert_eq!(configuration[0].0, 1); // n for 1s
    assert_eq!(configuration[0].1, 0); // l for 1s

    assert_eq!(configuration[1].0, 2); // n for 2s
    assert_eq!(configuration[1].1, 0); // l for 2s

    assert_eq!(configuration[2].0, 2); // n for 2p
    assert_eq!(configuration[2].1, 1); // l for 2p

    // Check occupations
    assert_relative_eq!(configuration[0].2, 1.0, epsilon = 0.01); // 1s fully occupied
    assert_relative_eq!(configuration[1].2, 1.0, epsilon = 0.01); // 2s fully occupied
    assert_relative_eq!(configuration[2].2, 1.0 / 3.0, epsilon = 0.01); // 2p partially occupied (2 of 6 electrons)

    // Test for iron (Z=26)
    let (grid, potential) = create_test_potential(26, 200, 30.0);

    let config = AtomSolverConfig {
        atomic_number: 26,
        relativistic: false,
        energy_tolerance: 1e-6,
        max_iterations: 100,
        use_experimental: false,
    };

    let solver = AtomSolver::new(config, grid, potential);

    // Get the electronic configuration
    let configuration = solver.electronic_configuration();

    // For now, we're using a simplified filling order for testing
    // Iron should have at least the 1s through 4s shells filled
    assert!(configuration.len() >= 4); // At least through 4s

    // Find the 3d orbital
    let mut found_3d = false;
    for &(n, l, occ) in &configuration {
        if n == 3 && l == 2 {
            found_3d = true;
            // 3d should be partially filled (6 of 10 electrons)
            assert_relative_eq!(occ, 0.6, epsilon = 0.01);
            break;
        }
    }

    assert!(found_3d, "3d orbital not found in iron configuration");
}

#[test]
fn test_multi_electron_density() {
    // Test for beryllium (Z=4) which has 1s²2s²
    let (grid, potential) = create_test_potential(4, 200, 30.0);

    let config = AtomSolverConfig {
        atomic_number: 4,
        relativistic: false,
        energy_tolerance: 1e-6,
        max_iterations: 100,
        use_experimental: false,
    };

    let mut solver = AtomSolver::new(config, grid, potential);

    // Solve for the states
    solver.solve_state(1, 0).unwrap(); // 1s
    solver.solve_state(2, 0).unwrap(); // 2s

    // Set occupations for beryllium
    let occupations = vec![(1, 0, 1.0), (2, 0, 1.0)]; // 1s and 2s fully occupied

    // Calculate the total density
    let density = solver.calculate_density(&occupations).unwrap();

    // Integrate to find total number of electrons
    let mut total_electrons = 0.0;
    for i in 1..density.len() {
        let r1 = solver.grid[i - 1];
        let r2 = solver.grid[i];
        let rho1 = density[i - 1];
        let rho2 = density[i];

        // Volume element
        let vol1 = 4.0 * PI * r1 * r1;
        let vol2 = 4.0 * PI * r2 * r2;

        // Integrate using trapezoidal rule
        let dr = r2 - r1;
        total_electrons += 0.5 * dr * (rho1 * vol1 + rho2 * vol2);
    }

    // Beryllium has 4 electrons
    assert_relative_eq!(total_electrons, 4.0, epsilon = 0.1);

    // Calculate total energy
    let total_energy = solver.calculate_total_energy(&occupations).unwrap();

    // The energy should be negative
    assert!(total_energy < 0.0);
}
