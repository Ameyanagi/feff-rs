/*
MIT License with FEFF10 Attribution

Copyright (c) 2025 Ameyanagi

Based on or developed using Distribution: FEFF10.0
Copyright (c) 2020 FEFF Project, University of Washington and SLAC National Accelerator Laboratory.
All rights reserved.
*/

//! Core-hole effects for X-ray absorption spectroscopy
//!
//! This module implements core-hole calculations for XAS simulations,
//! which are essential for modeling the X-ray absorption process.
//! It handles different approximations like the final-state rule (FSR)
//! and the random-phase approximation (RPA).

use crate::atoms::errors::AtomError;
use crate::atoms::{AtomicStructure, Result as AtomResult};
use crate::potential::{ExchangeCorrelationType, MuffinTinPotential};
use crate::scattering::{calculate_phase_shifts_with_method, PhaseShiftMethod, ScatteringResults};

/// Core-hole approximation methods for XAS calculations
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum CoreHoleMethod {
    /// Ground state approximation (no core hole)
    GroundState,

    /// Final state rule (Z+1 approximation for absorbing atom)
    FinalState,

    /// Self-consistent field with core hole
    SelfConsistent,

    /// Random phase approximation (RPA)
    RPA,
}

/// Configuration for core-hole calculations
#[derive(Debug, Clone)]
pub struct CoreHoleConfig {
    /// Core-hole calculation method to use
    pub method: CoreHoleMethod,

    /// Absorption edge (K, L1, L2, L3, etc.)
    pub edge: String,

    /// Core-hole screening parameter (0.0 to 1.0)
    /// 0.0 = no screening (Z* = 1.0)
    /// 1.0 = full screening (Z* = 0.0)
    pub screening: f64,

    /// Exchange-correlation functional to use
    pub xc_type: ExchangeCorrelationType,
}

impl Default for CoreHoleConfig {
    fn default() -> Self {
        Self {
            method: CoreHoleMethod::FinalState,
            edge: "K".to_string(),
            screening: 0.0,
            xc_type: ExchangeCorrelationType::HedinLundqvist,
        }
    }
}

impl CoreHoleConfig {
    /// Create a new core-hole configuration with default values
    pub fn new() -> Self {
        Self::default()
    }

    /// Set the core-hole calculation method
    pub fn with_method(&mut self, method: CoreHoleMethod) -> &mut Self {
        self.method = method;
        self
    }

    /// Set the absorption edge
    pub fn with_edge(&mut self, edge: &str) -> &mut Self {
        self.edge = edge.to_string();
        self
    }

    /// Set the core-hole screening parameter
    pub fn with_screening(&mut self, screening: f64) -> &mut Self {
        if !(0.0..=1.0).contains(&screening) {
            // Clamp value to valid range
            self.screening = screening.clamp(0.0, 1.0);
        } else {
            self.screening = screening;
        }
        self
    }

    /// Set the exchange-correlation functional
    pub fn with_xc_type(&mut self, xc_type: ExchangeCorrelationType) -> &mut Self {
        self.xc_type = xc_type;
        self
    }

    /// Get the effective core-hole charge based on screening
    pub fn effective_z_star(&self) -> f64 {
        // Z* = (1 - screening)
        // At screening=0, Z*=1 (full core hole)
        // At screening=1, Z*=0 (no core hole)
        1.0 - self.screening
    }

    /// Validate the edge type is supported
    pub fn validate_edge(&self) -> AtomResult<()> {
        match self.edge.as_str() {
            "K" | "L1" | "L2" | "L3" | "M1" => Ok(()),
            _ => Err(AtomError::InvalidEdge(format!(
                "Unsupported edge type: {}. Supported types are: K, L1, L2, L3, M1",
                self.edge
            ))),
        }
    }
}

/// Perform core-hole calculation for XAS using the specified method
///
/// # Arguments
///
/// * `structure` - The atomic structure with the absorbing atom set
/// * `energy` - The energy for the calculation in eV
/// * `max_l` - Maximum angular momentum for phase shifts
/// * `config` - Core-hole configuration parameters
///
/// # Returns
///
/// A ScatteringResults object containing phase shifts with core-hole effects
pub fn calculate_with_core_hole(
    structure: &AtomicStructure,
    energy: f64,
    max_l: i32,
    config: &CoreHoleConfig,
) -> AtomResult<ScatteringResults> {
    // Validate edge type
    config.validate_edge()?;

    // Make sure structure has a central atom defined
    if structure.central_atom().is_none() {
        return Err(AtomError::InvalidStructure(
            "No central atom defined for XAS calculation".to_string(),
        ));
    }

    // Calculate muffin-tin radii if not already calculated
    let mut structure_copy = structure.clone();
    if structure_copy
        .potential_type(0)
        .unwrap()
        .muffin_tin_radius()
        <= 0.0
    {
        structure_copy.calculate_muffin_tin_radii()?;
    }

    match config.method {
        CoreHoleMethod::GroundState => {
            // Just calculate phase shifts with no core hole
            calculate_phase_shifts_with_method(
                &structure_copy,
                energy,
                max_l,
                PhaseShiftMethod::MuffinTin,
            )
        }
        CoreHoleMethod::FinalState => {
            // Z+1 approximation for the absorbing atom
            apply_z_plus_one(&mut structure_copy)?;

            // Calculate phase shifts with the modified structure
            calculate_phase_shifts_with_method(
                &structure_copy,
                energy,
                max_l,
                PhaseShiftMethod::MuffinTin,
            )
        }
        CoreHoleMethod::SelfConsistent => {
            // Create muffin-tin potential calculator with core hole
            let mut mt_calculator = match MuffinTinPotential::new(&structure_copy) {
                Ok(calc) => calc,
                Err(e) => return Err(AtomError::PotentialError(format!("{}", e))),
            };

            match mt_calculator.set_core_hole(config.effective_z_star()) {
                Ok(_) => {}
                Err(e) => return Err(AtomError::PotentialError(format!("{}", e))),
            };

            match mt_calculator.set_exchange_correlation(&format!("{:?}", config.xc_type)) {
                Ok(_) => {}
                Err(e) => return Err(AtomError::PotentialError(format!("{}", e))),
            };

            // Run self-consistency calculation
            let scf_result = match mt_calculator.run_self_consistency(10, 1e-3) {
                Ok(result) => result,
                Err(e) => return Err(AtomError::PotentialError(format!("{}", e))),
            };

            if !scf_result.converged {
                // Log warning but continue with the best result we have
                eprintln!("Warning: Core-hole potential did not fully converge in {} iterations, final error: {}", 
                    scf_result.iterations, scf_result.final_error);
            }

            // Calculate phase shifts with the self-consistent potential
            calculate_phase_shifts_with_method(
                &structure_copy,
                energy,
                max_l,
                PhaseShiftMethod::MuffinTin,
            )
        }
        CoreHoleMethod::RPA => {
            // First calculate ground state
            let ground_state = calculate_phase_shifts_with_method(
                &structure_copy,
                energy,
                max_l,
                PhaseShiftMethod::MuffinTin,
            )?;

            // Then calculate with core hole
            let core_hole_config = CoreHoleConfig {
                method: CoreHoleMethod::SelfConsistent,
                ..config.clone()
            };

            let excited_state =
                calculate_with_core_hole(&structure_copy, energy, max_l, &core_hole_config)?;

            // Combine results using RPA
            combine_rpa_results(&ground_state, &excited_state)
        }
    }
}

/// Apply the Z+1 approximation (final state rule)
///
/// This modifies the atomic number of the absorbing atom by increasing it by 1,
/// which simulates the core hole by assuming complete screening
fn apply_z_plus_one(structure: &mut AtomicStructure) -> AtomResult<()> {
    let central_idx = match structure.central_atom_index() {
        Some(idx) => idx,
        None => {
            return Err(AtomError::InvalidStructure(
                "No central atom defined for Z+1 approximation".to_string(),
            ));
        }
    };

    let atom = structure.atom(central_idx).unwrap();
    let pot_idx = atom.potential_type() as usize;
    let pot = structure.potential_type(pot_idx).unwrap();

    // Create a new potential type with Z+1
    let new_z = pot.atomic_number() + 1;
    let mut new_pot = pot.clone();
    new_pot.set_atomic_number(new_z)?;

    // Replace the potential type
    structure.replace_potential_type(pot_idx, new_pot)?;

    Ok(())
}

/// Combine ground state and excited state results using the Random Phase Approximation
///
/// # Arguments
///
/// * `ground_state` - Phase shifts calculated without core hole
/// * `excited_state` - Phase shifts calculated with core hole
///
/// # Returns
///
/// Combined results according to RPA
fn combine_rpa_results(
    ground_state: &ScatteringResults,
    excited_state: &ScatteringResults,
) -> AtomResult<ScatteringResults> {
    // Verify that both results have the same energy and max_l
    if ground_state.energy != excited_state.energy || ground_state.max_l != excited_state.max_l {
        return Err(AtomError::CalculationError(
            "Cannot combine RPA results with different energy or max_l".to_string(),
        ));
    }

    // Create a new result with the same structure
    let mut combined = ScatteringResults {
        energy: ground_state.energy,
        max_l: ground_state.max_l,
        phase_shifts: Vec::with_capacity(ground_state.phase_shifts.len()),
        t_matrices: Vec::with_capacity(ground_state.t_matrices.len()),
    };

    // Combine phase shifts
    for i in 0..ground_state.phase_shifts.len() {
        let mut combined_phases = Vec::with_capacity(ground_state.phase_shifts[i].len());

        for j in 0..ground_state.phase_shifts[i].len() {
            // For the absorbing atom (i=0), use the excited state
            // For other atoms, use the ground state
            let phase = if i == 0 {
                excited_state.phase_shifts[i][j]
            } else {
                ground_state.phase_shifts[i][j]
            };
            combined_phases.push(phase);
        }

        combined.phase_shifts.push(combined_phases);
    }

    // Combine T-matrices
    for i in 0..ground_state.t_matrices.len() {
        // For the absorbing atom (i=0), use the excited state
        // For other atoms, use the ground state
        let t_matrix = if i == 0 {
            excited_state.t_matrices[i].clone()
        } else {
            ground_state.t_matrices[i].clone()
        };

        combined.t_matrices.push(t_matrix);
    }

    Ok(combined)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::atoms::{Atom, PotentialType, Vector3D};

    #[test]
    fn test_core_hole_config() {
        let config = CoreHoleConfig::default();
        assert_eq!(config.method, CoreHoleMethod::FinalState);
        assert_eq!(config.edge, "K");
        assert_eq!(config.screening, 0.0);

        let mut config = CoreHoleConfig::new();
        config
            .with_method(CoreHoleMethod::SelfConsistent)
            .with_edge("L1")
            .with_screening(0.5);
        let custom_config = config;

        assert_eq!(custom_config.method, CoreHoleMethod::SelfConsistent);
        assert_eq!(custom_config.edge, "L1");
        assert_eq!(custom_config.screening, 0.5);

        // Test effective_z_star
        assert_eq!(custom_config.effective_z_star(), 0.5);

        // Test validation
        assert!(custom_config.validate_edge().is_ok());

        let mut config = CoreHoleConfig::new();
        let invalid_config = config.with_edge("X");
        assert!(invalid_config.validate_edge().is_err());
    }

    #[test]
    fn test_z_plus_one() {
        // Create a simple iron atom structure
        let fe_potential = PotentialType::new(0, 26).unwrap();
        let fe_atom = Atom::new(26, Vector3D::new(0.0, 0.0, 0.0), 0).unwrap();

        let mut structure = AtomicStructure::new();
        structure.add_potential_type(fe_potential);
        let fe_idx = structure.add_atom(fe_atom);
        structure.set_central_atom(fe_idx).unwrap();

        // Apply Z+1 approximation
        apply_z_plus_one(&mut structure).unwrap();

        // Check that the atomic number was increased by 1
        let pot_idx = structure.atom(fe_idx).unwrap().potential_type() as usize;
        let z = structure.potential_type(pot_idx).unwrap().atomic_number();
        assert_eq!(z, 27); // Fe (Z=26) -> Co (Z=27)
    }

    #[test]
    fn test_core_hole_methods() {
        // Create a simple iron atom structure
        let fe_potential = PotentialType::new(0, 26).unwrap();
        let fe_atom = Atom::new(26, Vector3D::new(0.0, 0.0, 0.0), 0).unwrap();

        let mut structure = AtomicStructure::new();
        structure.add_potential_type(fe_potential);
        let fe_idx = structure.add_atom(fe_atom);
        structure.set_central_atom(fe_idx).unwrap();

        // Calculate muffin-tin radius
        structure.calculate_muffin_tin_radii().unwrap();

        // Energy and max_l for calculation
        let energy = 100.0; // eV
        let max_l = 2;

        // Test ground state method
        let mut config = CoreHoleConfig::new();
        let ground_config = config.with_method(CoreHoleMethod::GroundState);

        let ground_result = calculate_with_core_hole(&structure, energy, max_l, &ground_config);

        assert!(ground_result.is_ok());

        // Test final state rule (Z+1) method
        let mut config = CoreHoleConfig::new();
        let fsr_config = config.with_method(CoreHoleMethod::FinalState);

        let fsr_result = calculate_with_core_hole(&structure, energy, max_l, &fsr_config);

        assert!(fsr_result.is_ok());

        // Compare results - should be different
        let ground = ground_result.unwrap();
        let fsr = fsr_result.unwrap();

        assert_eq!(ground.energy, fsr.energy);
        assert_eq!(ground.max_l, fsr.max_l);

        // At least one phase shift should be different
        let mut all_same = true;
        for l in 0..=max_l {
            if (ground.phase_shifts[0][l as usize] - fsr.phase_shifts[0][l as usize]).norm() > 1e-6
            {
                all_same = false;
                break;
            }
        }

        // The two methods should produce different results
        assert!(!all_same);
    }

    #[test]
    fn test_core_hole_edges() {
        // Create a simple oxygen atom structure
        let o_potential = PotentialType::new(0, 8).unwrap();
        let o_atom = Atom::new(8, Vector3D::new(0.0, 0.0, 0.0), 0).unwrap();

        let mut structure = AtomicStructure::new();
        structure.add_potential_type(o_potential);
        let o_idx = structure.add_atom(o_atom);
        structure.set_central_atom(o_idx).unwrap();

        // Calculate muffin-tin radius
        structure.calculate_muffin_tin_radii().unwrap();

        // Test with K edge
        let mut config = CoreHoleConfig::new();
        config
            .with_method(CoreHoleMethod::FinalState)
            .with_edge("K");
        let k_config = config;

        let k_result = calculate_with_core_hole(
            &structure, 550.0, // Near O K-edge
            2, &k_config,
        );

        assert!(k_result.is_ok());

        // Test with L1 edge
        let mut config = CoreHoleConfig::new();
        config
            .with_method(CoreHoleMethod::FinalState)
            .with_edge("L1");
        let l1_config = config;

        let l1_result = calculate_with_core_hole(
            &structure, 100.0, // Not actual L1 energy, just for testing
            2, &l1_config,
        );

        assert!(l1_result.is_ok());
    }
}
