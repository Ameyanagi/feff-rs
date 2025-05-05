# FEFF-rs

A modern Rust implementation of the FEFF code for calculating X-ray absorption spectroscopy (XAS) and related spectroscopies.

[![CI](https://github.com/ameyanagi/feff-rs/actions/workflows/ci.yml/badge.svg)](https://github.com/ameyanagi/feff-rs/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/ameyanagi/feff-rs/branch/main/graph/badge.svg)](https://codecov.io/gh/ameyanagi/feff-rs)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

FEFF-rs aims to provide a high-performance, memory-efficient, and maintainable implementation of the FEFF code in Rust. FEFF is a widely used ab initio multiple-scattering code for calculations of excitation spectra and electronic structure.

This project is under active development and is not yet ready for production use.

## Features (Planned)

- Parse FEFF input files
- Calculate atomic potentials
- Find scattering paths
- Compute X-ray Absorption Near Edge Structure (XANES)
- Compute Extended X-ray Absorption Fine Structure (EXAFS)
- Support multiple spectroscopies (XANES, EXAFS, EELS, XES, etc.)
- Parallel processing for improved performance
- Modern CLI with good error messages

## Project Structure

```
feff-rs/
├── src/              # Source code
│   ├── atoms/        # Atomic data and calculations
│   ├── input/        # Input file parsing
│   ├── potential/    # Potential calculation
│   ├── scattering/   # Scattering calculations
│   ├── path/         # Path finding and filtering
│   ├── fms/          # Full multiple scattering
│   ├── xas/          # XAS spectrum calculations
│   ├── utils/        # Utility functions
│   └── cli/          # Command-line interface
├── tests/            # Integration tests
├── benches/          # Performance benchmarks
├── docs/             # Documentation
├── plans/            # Development plans
└── .github/          # GitHub workflows
```

## Development

### Getting Started

1. Clone the repository
2. Install pre-commit hooks using `uvx`:
   ```bash
   # Install UV if you haven't already
   curl -LsSf https://astral.sh/uv/install.sh | sh

   # Install the pre-commit hooks
   uvx pre-commit install
   ```

### Building

To build the project, run:

```bash
cargo build
```

For a release build:

```bash
cargo build --release
```

### Testing

Run the tests with:

```bash
cargo test
```

### Benchmarks

Run benchmarks with:

```bash
cargo bench
```

### Code Style

This project follows the Rust standard formatting. To check your code:

```bash
cargo fmt -- --check
```

To automatically format the code:

```bash
cargo fmt
```

### Linting

Run the linter with:

```bash
cargo clippy
```

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

### Pre-commit Hooks

This project uses pre-commit hooks with `uvx pre-commit` to ensure code quality. The hooks run:
- Clippy linting with warnings as errors
- Rustfmt checks
- Cargo check compilation
- Rustfmt.toml validation
- Unit tests
- Check for yanked dependencies

When you commit changes, these checks will run automatically if you've installed the hooks.

## Project Status

This project is in early development. See the [architecture document](plans/architecture.md) for the planned design.

## License

MIT License with FEFF10 Attribution

Copyright (c) 2025 Ameyanagi

This project is based on the FEFF10 code developed by the FEFF Project at the University of Washington and SLAC National Accelerator Laboratory. See the [LICENSE](LICENSE) and [FEFF10_LICENSE](FEFF10_LICENSE) files for details.

## Citation

When using this software, please cite both this implementation and the original FEFF papers:

```
FEFF-rs: A Rust implementation of FEFF
[Citation information to be added]

The original FEFF project:
J.J. Rehr et al., Phys. Rev. B, 80, 115112 (2009)
```

## Acknowledgments

We acknowledge the FEFF Project at the University of Washington and SLAC National Accelerator Laboratory for developing the original FEFF code.