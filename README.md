# FEFF-rs

A modern Rust implementation of the FEFF code for calculating X-ray absorption spectroscopy (XAS) and related spectroscopies.

[![CI](https://github.com/ameyanagi/feff-rs/actions/workflows/ci.yml/badge.svg)](https://github.com/ameyanagi/feff-rs/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/ameyanagi/feff-rs/branch/main/graph/badge.svg)](https://codecov.io/gh/ameyanagi/feff-rs)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

FEFF-rs aims to provide a high-performance, memory-efficient, and maintainable implementation of the FEFF code in Rust. FEFF is a widely used ab initio multiple-scattering code for calculations of excitation spectra and electronic structure.

This project is under active development and is not yet ready for production use.

## Features

- Parse FEFF input files
- Calculate atomic potentials
- Find scattering paths
- Compute X-ray Absorption Near Edge Structure (XANES)
- Compute Extended X-ray Absorption Fine Structure (EXAFS)
- Model temperature effects on XAS with various thermal models
- Support for anisotropic thermal vibrations
- Parallel processing for improved performance
- Modern CLI with good error messages

More features are planned for future releases.

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

## Thermal Models

FEFF-rs includes several thermal models for accounting for temperature effects in XAS:

- **Debye Model**: For simple cubic materials and metals
- **Einstein Model**: For materials with localized vibrations and stiff bonds
- **Correlated Debye Model**: For accurate EXAFS analysis with correlation effects
- **Anisotropic Thermal Models**: For non-cubic materials with direction-dependent vibrations

Example usage:

```rust
use feff_rs::xas::thermal::{ThermalParameters, create_anisotropic_thermal_parameters};

// Create thermal parameters for copper at room temperature 
let params = ThermalParameters::new_debye(300.0, 315.0);

// Or for anisotropic materials like tetragonal structures
let tetragonal_params = create_anisotropic_thermal_parameters(
    300.0,         // Temperature in K
    "tetragonal",  // Crystal system
    315.0,         // Debye temperature
    Some(1.5)      // Anisotropy ratio
);
```

For more details, see the [thermal effects documentation](docs/thermal_effects.md) and the [thermal effects example](docs/thermal_effects_example.md).

## Project Structure

```
feff-rs/
├── src/              # Source code
│   ├── atoms/        # Atomic data and structure
│   ├── input/        # Input file parsing
│   ├── potential/    # Potential calculations
│   ├── scattering/   # Scattering calculations
│   ├── path/         # Path finding algorithms
│   ├── fms/          # Full multiple scattering
│   ├── xas/          # X-ray absorption calculations
│   └── utils/        # Utility functions and constants
├── examples/         # Example applications
├── tests/            # Tests for various components
└── docs/             # Documentation
```

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This project is licensed under the MIT License with FEFF10 Attribution - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- Based on or developed using Distribution: FEFF10.0
- Copyright (c) 2020 FEFF Project, University of Washington and SLAC National Accelerator Laboratory