# FEFF-rs Library Architecture and API Reference

This document provides a comprehensive overview of the FEFF-rs library architecture, core modules, and API usage patterns for developers.

## Core Architecture

FEFF-rs is designed as a modular library that follows the physics of X-ray absorption spectroscopy (XAS) calculations. The architecture consists of several layers, with each layer building on the previous one:

```
┌────────────────────────────────────────────┐
│ Top-Level API (Feff struct)                │
└────────────────────────────────────────────┘
                     │
┌────────────────────────────────────────────┐
│ XAS Calculations (XANES, EXAFS)            │
└────────────────────────────────────────────┘
                     │
┌────────────────────┬─────────────────────┐
│ FMS Calculations   │ Path Finding        │
└────────────────────┴─────────────────────┘
                     │
┌────────────────────────────────────────────┐
│ Scattering Calculations                    │
└────────────────────────────────────────────┘
                     │
┌────────────────────────────────────────────┐
│ Potential Calculations                     │
└────────────────────────────────────────────┘
                     │
┌────────────────────┬─────────────────────┐
│ Atomic Structure   │ Input Parsing       │
└────────────────────┴─────────────────────┘
                     │
┌────────────────────────────────────────────┐
│ Utilities and Constants                    │
└────────────────────────────────────────────┘
```

## Module Structure

### 1. Core Data Types (`atoms/`)

The foundation of FEFF-rs starts with the atomic structure representation.

#### Key Types

- `Vector3D`: Three-dimensional vector for atomic positions and directions
- `Atom`: Represents an atom with position, atomic number, and potential index
- `AtomicStructure`: Collection of atoms with their relationships and properties
- `PotentialType`: Defines the potential properties for atoms

#### Example: Creating an Atomic Structure

```rust
use feff_rs::atoms::{Atom, AtomicStructure, PotentialType, Vector3D};

// Create an empty structure
let mut structure = AtomicStructure::new();

// Add potential types
let cu_potential = PotentialType::new(0, 29).unwrap(); // index 0, Cu (Z=29)
structure.add_potential_type(cu_potential);

// Add atoms
let cu_central = Atom::new(29, Vector3D::new(0.0, 0.0, 0.0), 0).unwrap();
let cu_neighbor = Atom::new(29, Vector3D::new(2.55, 0.0, 0.0), 0).unwrap();

let central_idx = structure.add_atom(cu_central);
structure.add_atom(cu_neighbor);

// Set the central atom (absorber)
structure.set_central_atom(central_idx).unwrap();
```

### 2. Input Parsing (`input/`)

This module handles the parsing and processing of FEFF input files.

#### Key Types

- `FeffInput`: Top-level structure representing a FEFF input file
- `Card`: Representation of individual FEFF input cards
- `CoordinateSystem`: Enumeration of coordinate systems (cartesian, spherical, etc.)

#### Example: Parsing an Input File

```rust
use feff_rs::input::FeffInput;
use std::path::Path;

// Parse a FEFF input file
let input = FeffInput::from_file(Path::new("path/to/feff.inp")).unwrap();

// Access parsed data
let atoms = input.get_atoms().unwrap();
let potentials = input.get_potentials().unwrap();
let parameters = input.get_parameters();
```

### 3. Potential Calculations (`potential/`)

Handles atomic potential calculations with various exchange-correlation functionals.

#### Key Types

- `MuffinTin`: Represents the muffin-tin potential
- `AtomSolver`: Solves atomic Schrödinger equation
- `ExchangeCorrelation`: Trait for various exchange-correlation potentials

#### Example: Calculating Potentials

```rust
use feff_rs::potential::{calculate_potentials, AtomSolverParams};

// Calculate potentials for the atomic structure
let params = AtomSolverParams::default();
let potentials = calculate_potentials(&structure, &params).unwrap();
```

### 4. Scattering Calculations (`scattering/`)

Performs phase shift calculations and prepares scattering matrices.

#### Key Types

- `PhaseShifts`: Container for phase shifts
- `ScatteringMatrices`: Collection of scattering matrices
- `PhaseShiftCalculator`: Engine for calculating phase shifts

#### Example: Calculating Phase Shifts

```rust
use feff_rs::scattering::{calculate_phase_shifts, PhaseShiftParams};

// Calculate phase shifts
let params = PhaseShiftParams::default();
let phase_shifts = calculate_phase_shifts(&structure, &potentials, &params).unwrap();
```

### 5. Path Finding (`path/`)

Finds and filters scattering paths for EXAFS calculations.

#### Key Types

- `Path`: Representation of a scattering path
- `PathFinder`: Engine for finding all relevant paths
- `PathFilter`: Filters paths based on importance criteria

#### Example: Finding EXAFS Paths

```rust
use feff_rs::path::{PathFinder, PathParameters};

// Set up path finding parameters
let params = PathParameters {
    max_path_length: 8.0,
    max_legs: 4,
    max_paths: 100,
    min_importance: 0.1,
};

// Find paths
let path_finder = PathFinder::new(&structure);
let paths = path_finder.find_paths(&params).unwrap();

// Access path data
for (i, path) in paths.iter().enumerate() {
    println!("Path {}: deg={}, importance={:.3}", 
             i, path.degeneracy, path.amplitude / paths[0].amplitude);
}
```

### 6. Full Multiple Scattering (`fms/`)

Implements the Full Multiple Scattering (FMS) approach for XANES calculations.

#### Key Types

- `FmsSolver`: Engine for solving FMS equations
- `FmsParameters`: Parameters for FMS calculations
- `FmsResults`: Results of FMS calculations

#### Example: FMS Calculation

```rust
use feff_rs::fms::{FmsSolver, FmsParameters};
use feff_rs::scattering::calculate_scattering_matrices;

// Calculate scattering matrices
let scattering_matrices = calculate_scattering_matrices(
    &structure, &phase_shifts, &energy_grid).unwrap();

// Set up FMS parameters
let params = FmsParameters {
    max_l: 3,
    convergence_threshold: 1e-6,
    ..Default::default()
};

// Solve FMS equations
let fms_solver = FmsSolver::new(&structure, &scattering_matrices, params);
let fms_results = fms_solver.solve().unwrap();
```

### 7. XAS Calculations (`xas/`)

Calculates X-ray absorption spectra (XANES and EXAFS).

#### Key Types

- `XanesParameters`: Parameters for XANES calculations
- `XanesSpectrum`: Results of XANES calculations
- `ExafsParameters`: Parameters for EXAFS calculations
- `ExafsSpectrum`: Results of EXAFS calculations
- `ThermalParameters`: Parameters for modeling thermal effects

#### Example: Calculating XANES Spectrum

```rust
use feff_rs::xas::{calculate_xanes, XanesParameters, Edge};

// Set up XANES parameters
let params = XanesParameters {
    edge: Edge::K,
    energy_range: (-10.0, 50.0, 0.5), // min, max, step in eV
    polarization: None, // isotropic spectrum
    core_hole_method: CoreHoleMethod::FinalState,
    ..Default::default()
};

// Calculate XANES spectrum
let spectrum = calculate_xanes(&structure, &params).unwrap();

// Access results
for (i, energy) in spectrum.energies.iter().enumerate() {
    println!("Energy: {:.2} eV, Absorption: {:.6}", energy, spectrum.mu[i]);
}
```

#### Example: Calculating EXAFS Spectrum

```rust
use feff_rs::xas::{calculate_exafs, ExafsParameters, Edge, WindowFunction};

// Set up EXAFS parameters
let params = ExafsParameters {
    edge: Edge::K,
    k_range: (3.0, 14.0), // k-range in Å⁻¹
    r_range: (0.0, 6.0, 0.05), // min, max, step in Å
    max_path_length: 8.0,
    max_legs: 4,
    max_paths: 100,
    window_function: WindowFunction::Kaiser(2.0),
    ..Default::default()
};

// Calculate EXAFS spectrum
let spectrum = calculate_exafs(&structure, &params).unwrap();

// Access results
for (i, k) in spectrum.k.iter().enumerate() {
    println!("k: {:.2} Å⁻¹, chi(k): {:.6}", k, spectrum.chi_k[i]);
}
```

### 8. Thermal Effects (`xas/thermal.rs` and `utils/thermal.rs`)

Handles temperature effects on XAS spectra with various thermal models.

#### Key Types

- `ThermalParameters`: Configuration for thermal calculations
- `ThermalModel`: Trait for different thermal models
  - `DebyeModel`: Standard thermal model for simple materials
  - `EinsteinModel`: Model for localized vibrations
  - `CorrelatedDebyeModel`: Model with correlation effects
  - `AnisotropicThermalModel`: Model for directional thermal vibrations

#### Example: Using Thermal Models for XANES

```rust
use feff_rs::xas::{calculate_xanes, XanesParameters, Edge};
use feff_rs::xas::thermal::ThermalParameters;

// Create thermal parameters for copper at room temperature
let thermal_params = ThermalParameters::new_debye(300.0, 315.0);

// Set up XANES parameters with thermal effects
let params = XanesParameters {
    edge: Edge::K,
    energy_range: (-10.0, 50.0, 0.5),
    thermal_parameters: Some(thermal_params),
    ..Default::default()
};

// Calculate temperature-dependent XANES
let spectrum = calculate_xanes(&structure, &params).unwrap();
```

#### Example: Anisotropic Thermal Effects

```rust
use feff_rs::xas::thermal::{ThermalParameters, create_anisotropic_thermal_parameters};

// Create thermal parameters for a tetragonal crystal
let thermal_params = create_anisotropic_thermal_parameters(
    300.0,          // Temperature in K
    "tetragonal",   // Crystal system
    400.0,          // Debye temperature
    Some(1.7)       // Anisotropy ratio
);

// Use in XAS calculations
let params = XanesParameters {
    thermal_parameters: Some(thermal_params),
    ..Default::default()
};
```

## Top-Level API (`Feff` struct)

The `Feff` struct in `src/lib.rs` provides a high-level API for working with FEFF calculations.

#### Key Methods

- `new()`: Create a new FEFF calculation instance
- `from_file(path)`: Create a FEFF calculation from an input file
- `run()`: Run the calculation with the provided input
- `calculate_xanes(params)`: Calculate a XANES spectrum
- `calculate_exafs(params)`: Calculate an EXAFS spectrum

#### Example: Complete Workflow with the Feff Struct

```rust
use feff_rs::{Feff, XanesParameters, Edge};

// Create a FEFF calculation from an input file
let mut feff = Feff::from_file("path/to/feff.inp").unwrap();

// Run the initial calculation (potentials, phase shifts)
feff.run().unwrap();

// Calculate XANES with custom parameters
let params = XanesParameters {
    edge: Edge::K,
    energy_range: (-10.0, 50.0, 0.5),
    polarization: None, // isotropic spectrum
    ..Default::default()
};

let spectrum = feff.calculate_xanes(&params).unwrap();

// Save the spectrum to a file
feff.save_spectrum("xanes_spectrum.dat", &spectrum).unwrap();
```

## Performance Considerations

FEFF-rs includes several performance optimizations:

1. **Parallel Processing**: Many calculations use Rayon for parallelization:

```rust
// Example of parallel processing configuration
use rayon::ThreadPoolBuilder;

// Set the number of threads to use
ThreadPoolBuilder::new().num_threads(8).build_global().unwrap();

// Now calculations will use this thread pool
let spectrum = calculate_xanes(&structure, &params).unwrap();
```

2. **Matrix Operations**: FEFF-rs uses the Faer library for high-performance linear algebra:

```rust
use feff_rs::utils::matrix::solve_linear_system;
use faer::Mat;

// Solve a linear system efficiently
let a: Mat<f64> = /* ... */;
let b: Mat<f64> = /* ... */;
let x = solve_linear_system(&a, &b).unwrap();
```

3. **Memory Management**: For large calculations, control memory usage:

```rust
// Example with energy grid chunking for large calculations
use feff_rs::xas::{calculate_xanes_chunked, XanesParameters};

// Calculate XANES in chunks to reduce memory usage
let params = XanesParameters {
    energy_range: (-10.0, 100.0, 0.1), // 1100 energy points
    ..Default::default()
};

let chunk_size = 200; // Process 200 energy points at a time
let spectrum = calculate_xanes_chunked(
    &structure, &params, chunk_size).unwrap();
```

## Error Handling

FEFF-rs uses the `thiserror` crate for library errors and provides specific error types for each module:

```rust
use feff_rs::atoms::errors::AtomError;
use feff_rs::potential::errors::PotentialError;
use feff_rs::scattering::errors::ScatteringError;
use feff_rs::fms::errors::FmsError;
use feff_rs::xas::errors::XasError;

// Example of error handling
fn run_calculation() -> Result<(), Box<dyn std::error::Error>> {
    let result = calculate_xanes(&structure, &params)
        .map_err(|e| match e {
            XasError::InvalidParameters(msg) => {
                // Handle parameter validation errors
                eprintln!("Invalid parameters: {}", msg);
                e
            }
            XasError::CalculationFailed(msg) => {
                // Handle calculation errors
                eprintln!("Calculation failed: {}", msg);
                e
            }
            _ => e
        })?;
    
    Ok(())
}
```

## Advanced Features

### 1. Polarized XAS Calculations

FEFF-rs supports polarized XAS calculations for studying orientation-dependent absorption:

```rust
use feff_rs::xas::{calculate_xanes, XanesParameters};

// Calculate polarized XANES with polarization vector along x-axis
let params = XanesParameters {
    polarization: Some([1.0, 0.0, 0.0]),
    ..Default::default()
};

let spectrum_x = calculate_xanes(&structure, &params).unwrap();

// Calculate polarized XANES with polarization vector along z-axis
let params_z = XanesParameters {
    polarization: Some([0.0, 0.0, 1.0]),
    ..Default::default()
};

let spectrum_z = calculate_xanes(&structure, &params_z).unwrap();
```

### 2. Core-Hole Effects

Different methods for handling the core-hole:

```rust
use feff_rs::xas::{XanesParameters, CoreHoleMethod};

// Final state approximation (fully screened core-hole)
let params_final = XanesParameters {
    core_hole_method: CoreHoleMethod::FinalState,
    ..Default::default()
};

// Self-consistent field approach
let params_scf = XanesParameters {
    core_hole_method: CoreHoleMethod::SelfConsistent,
    core_hole_screening: 0.8, // Partial core-hole screening
    ..Default::default()
};
```

### 3. EXAFS Fitting

FEFF-rs includes tools for EXAFS data analysis and fitting:

```rust
use feff_rs::xas::fitting::{ExafsFitter, ExafsFitParameters};

// Create a fitter for experimental data
let mut fitter = ExafsFitter::new(
    &experimental_k,     // k values from experiment
    &experimental_chi_k, // chi(k) values from experiment
    &structure           // Initial structure model
);

// Set fitting parameters
let fit_params = ExafsFitParameters {
    vary_distances: true,
    vary_debye_waller: true,
    vary_e0: true,
    vary_s02: true,
    ..Default::default()
};

// Perform the fit
let fit_result = fitter.fit(&fit_params).unwrap();

// Access fitted parameters
println!("Fitted S0²: {:.3}", fit_result.s02);
println!("Fitted ΔE0: {:.3} eV", fit_result.e0_shift);
for (i, shell) in fit_result.shells.iter().enumerate() {
    println!("Shell {}: R={:.3} Å, σ²={:.5} Å²", 
             i+1, shell.distance, shell.debye_waller);
}
```

## Best Practices

1. **Memory Management**:
   - Reuse structures when running multiple calculations
   - Use chunked calculations for large energy grids

2. **Performance Optimization**:
   - Use the appropriate level of theory for your needs
   - Limit path lengths and legs appropriately for EXAFS
   - Choose the right `max_l` parameter for FMS calculations

3. **Physics Accuracy**:
   - Validate against experimental data when possible
   - Use appropriate exchange-correlation potentials
   - Include thermal effects for realistic spectra
   - Configure core-hole parameters appropriately

4. **Error Handling**:
   - Always check returned `Result` values
   - Use pattern matching to handle specific error cases
   - Provide useful error messages to users

## References

1. J.J. Rehr, R.C. Albers, Rev. Mod. Phys. 72 (2000) 621-654.
2. A.L. Ankudinov, B. Ravel, J.J. Rehr, S.D. Conradson, Phys. Rev. B 58 (1998) 7565-7576.
3. Y. Joly, Phys. Rev. B 63 (2001) 125120.
4. A.L. Ankudinov, J.J. Rehr, Phys. Rev. B 62 (2000) 2437-2445.
5. J.J. Kas, A.P. Sorini, M.P. Prange, L.W. Cambell, J.A. Soininen, J.J. Rehr, Phys. Rev. B 76 (2007) 195116.