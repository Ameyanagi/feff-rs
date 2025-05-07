# FEFF-rs

A modern Rust implementation of the FEFF code for calculating X-ray absorption spectroscopy (XAS) and related spectroscopies.

[![CI](https://github.com/ameyanagi/feff-rs/actions/workflows/ci.yml/badge.svg)](https://github.com/ameyanagi/feff-rs/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/ameyanagi/feff-rs/branch/main/graph/badge.svg)](https://codecov.io/gh/ameyanagi/feff-rs)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

FEFF-rs provides a high-performance, memory-efficient, and maintainable implementation of the FEFF code in Rust. FEFF is a widely used ab initio multiple-scattering code for calculations of excitation spectra and electronic structure.

The software computes X-ray Absorption Spectra (XAS) including both X-ray Absorption Near Edge Structure (XANES) and Extended X-ray Absorption Fine Structure (EXAFS) from first principles based on multiple scattering theory.

## Features

- **Input Processing**: Parse FEFF input files with support for all standard FEFF cards
- **Atomic Data**: Comprehensive atomic database and structure handling
- **Potential Calculations**: Compute atomic potentials with various exchange-correlation functionals
- **Path Finding**: Find and filter scattering paths for EXAFS calculations
- **Multiple Scattering**: Full multiple scattering (FMS) implementation for accurate XANES
- **XAS Calculations**:
  - XANES: Near-edge spectroscopy with proper treatment of core-hole effects
  - EXAFS: Extended fine structure analysis with accurate phase shifts
- **Thermal Effects**: Model temperature effects on XAS with various thermal models:
  - Debye model for simple cubic materials
  - Einstein model for materials with localized vibrations
  - Correlated Debye model for accurate EXAFS analysis
  - Anisotropic thermal models for non-cubic materials
- **Performance Optimizations**:
  - Parallel processing for multi-core efficiency
  - High-performance linear algebra with the Faer library
  - Memory-efficient algorithms for large systems
- **Modern CLI**: User-friendly command-line interface with clear error messages

## Getting Started

### Installation

```bash
# Clone the repository
git clone https://github.com/ameyanagi/feff-rs.git
cd feff-rs

# Build the project
cargo build --release
```

### Running Examples

```bash
# Run XANES calculation example
cargo run --example xanes_calculation

# Run thermal effects example
cargo run --example thermal_effects_example
```

### Basic Usage

```rust
use feff_rs::{Feff, xas::{Edge, XanesParameters}};

// Create a FEFF calculation from an input file
let mut feff = Feff::from_file("path/to/feff.inp").unwrap();

// Run the calculation
feff.run().unwrap();

// Calculate XANES spectrum with custom parameters
let params = XanesParameters {
    edge: Edge::K,
    energy_range: (-10.0, 50.0, 0.5), // (min, max, step) in eV
    polarization: None, // Isotropic spectrum
    ..Default::default()
};

let spectrum = feff.calculate_xanes(&params).unwrap();

// Access the results
for (i, energy) in spectrum.energies.iter().enumerate() {
    println!("{:.2} eV: {:.6}", energy, spectrum.mu[i]);
}
```

## Thermal Models

FEFF-rs includes several thermal models for accounting for temperature effects in XAS:

### Debye Model

For simple cubic materials and metals with collective lattice vibrations.

```rust
use feff_rs::xas::thermal::ThermalParameters;

// Create thermal parameters for copper at room temperature 
let params = ThermalParameters::new_debye(300.0, 315.0);
```

### Einstein Model

For materials with localized vibrations and stiff bonds.

```rust
use feff_rs::xas::thermal::ThermalParameters;

// Create thermal parameters using Einstein model (70 meV frequency)
let params = ThermalParameters::new_einstein(300.0, 70.0);
```

### Correlated Debye Model

For accurate EXAFS analysis with correlation effects between atomic motions.

```rust
use feff_rs::xas::thermal::ThermalParameters;

// Create thermal parameters using correlated Debye model
let params = ThermalParameters::new_correlated_debye(300.0, 315.0);
```

### Anisotropic Thermal Models

For non-cubic materials with direction-dependent vibrations.

```rust
use feff_rs::xas::thermal::{ThermalParameters, create_anisotropic_thermal_parameters};

// Create parameters for a tetragonal material
let tetragonal_params = create_anisotropic_thermal_parameters(
    300.0,         // Temperature in K
    "tetragonal",  // Crystal system
    315.0,         // Debye temperature
    Some(1.5)      // Anisotropy ratio
);

// Create manually with specific displacement factors
let manual_params = ThermalParameters::new_anisotropic_debye(
    300.0,              // Temperature in K
    315.0,              // Debye temperature
    [1.0, 1.0, 2.0],    // Displacement factors [x, y, z]
    None                // Default path direction
);
```

### Material-Specific Thermal Parameters

```rust
use feff_rs::xas::thermal::create_material_thermal_parameters;

// Create parameters optimized for a specific material
let copper_params = create_material_thermal_parameters(
    "Cu",    // Material name
    300.0,   // Temperature in K
    None     // Use default Debye temperature
);

// Create parameters for a layered material
let graphite_params = create_material_thermal_parameters(
    "graphite",  // Material name
    300.0,       // Temperature in K
    None         // Use default Debye temperature
);
```

## Project Structure

```
feff-rs/
├── src/              # Source code
│   ├── atoms/        # Atomic data and structure
│   │   ├── atom.rs         # Atom implementation
│   │   ├── coordinates.rs  # Coordinate conversions
│   │   ├── database.rs     # Atomic database
│   │   ├── structure.rs    # AtomicStructure implementation
│   │   ├── vector.rs       # Vector3D implementation
│   │   └── potential.rs    # PotentialType implementation
│   ├── input/        # Input file parsing
│   │   ├── card.rs         # FEFF card definitions
│   │   ├── config.rs       # Configuration handling
│   │   ├── model.rs        # Input model
│   │   └── parser.rs       # Input file parser
│   ├── potential/    # Potential calculations
│   │   ├── atom_solver.rs      # Atomic potential solver
│   │   ├── muffin_tin.rs       # Muffin-tin potential
│   │   ├── scf.rs              # Self-consistent field
│   │   └── exchange_correlation.rs # XC potentials
│   ├── scattering/   # Scattering calculations
│   │   ├── phase_shifts.rs           # Phase shift storage
│   │   ├── phase_shift_calculator.rs # Phase shift calculation
│   │   └── scattering_matrices.rs    # Scattering matrices
│   ├── path/         # Path finding algorithms
│   │   ├── finder.rs        # Path finding implementation
│   │   ├── path.rs          # Path representation
│   │   ├── filter.rs        # Path filtering
│   │   └── degeneracy.rs    # Path degeneracy calculation
│   ├── fms/          # Full multiple scattering
│   │   ├── solver.rs        # FMS equation solver
│   │   ├── matrix.rs        # FMS matrix builders
│   │   └── xanes.rs         # XANES from FMS
│   ├── xas/          # X-ray absorption calculations
│   │   ├── xanes.rs         # XANES calculation
│   │   ├── exafs.rs         # EXAFS calculation
│   │   ├── thermal.rs       # Thermal parameters
│   │   ├── core_hole.rs     # Core-hole effects
│   │   └── fitting.rs       # EXAFS fitting
│   └── utils/        # Utility functions
│       ├── constants.rs     # Physical constants
│       ├── math.rs          # Mathematical functions
│       ├── thermal.rs       # Thermal models
│       └── matrix.rs        # Matrix operations
├── examples/         # Example applications
├── tests/            # Tests for various components
└── docs/             # Documentation
```

## Documentation

For detailed documentation:

```bash
# Generate and open documentation
cargo doc --open
```

Check the `docs/` directory for additional information on:
- [Thermal Effects](docs/thermal_effects.md)
- [XANES Calculation](docs/xanes_calculation.md)
- [EXAFS Analysis](docs/exafs_analysis.md)

## Performance

FEFF-rs is designed for high performance with several optimizations:

1. **Parallel Processing**: Uses the Rayon library to parallelize calculations across available CPU cores
2. **Efficient Matrix Operations**: Leverages the Faer library for high-performance linear algebra
3. **Memory Efficiency**: Carefully manages memory allocation for large calculations
4. **SIMD Acceleration**: Takes advantage of SIMD instructions where available

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

Before submitting your contribution, please:
1. Ensure your code follows the Rust standard formatting with `cargo fmt`
2. Run the linter with `cargo clippy -- -D warnings`
3. Run all tests with `cargo test`
4. Add tests for any new functionality

## License

This project is licensed under the MIT License with FEFF10 Attribution - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- Based on or developed using Distribution: FEFF10.0
- Copyright (c) 2020 FEFF Project, University of Washington and SLAC National Accelerator Laboratory