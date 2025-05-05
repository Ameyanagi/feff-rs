# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build/Test Commands
- Build: `cargo build`
- Run: `cargo run`
- Test all: `cargo test`
- Test specific: `cargo test test_name`
- Test specific module: `cargo test module_name`
- Lint: `cargo clippy -- -D warnings`
- Format: `cargo fmt`
- Benchmark: `cargo bench`
- Documentation: `cargo doc --open`

## Development Approach
- **FEFF10 First**: Always consult the original FEFF10 source code before implementing any functionality
- **TDD**: Write failing tests first, then implement code to make them pass
- **Error Handling**: Use anyhow for application errors, thiserror for library errors
- **Parallelism**: Use rayon for data parallelism where applicable

## Development Workflow
1. ALWAYS check the FEFF10 source code to understand the original implementation
2. Create a detailed implementation plan before writing any code
3. Document the current state and intended outcome
4. For each component:
   - Study FEFF10's approach to the problem thoroughly
   - Design the API and data structures to match FEFF10's physics while using idiomatic Rust
   - Write failing tests that validate expected behavior against FEFF10 reference values
   - Implement the code to make tests pass
   - Run tests to verify implementation matches FEFF10 results
   - Document the implementation with thorough rustdoc comments including physics explanation
5. Continuously update README.md with progress and usage instructions

## Documentation Retrieval
- Use context7 for retrieving the latest documentation on Rust, libraries
- Faer documentation: https://docs.rs/faer/latest/faer/
  - High-performance linear algebra library optimized for medium to large dense matrices
  - Provides matrix decompositions (Cholesky, LU, QR, SVD, Eigendecomposition)
  - Supports parallel processing via Rayon and SIMD acceleration
  - Useful for physics calculations involving matrix operations

## Library Dependencies
- anyhow = "1.0.98" - For application error handling
- thiserror = "2.0.12" - For library error type definitions
- clap = { version = "4.5.37", features = ["derive"] } - For command-line argument parsing
- rayon = "1.10.0" - For parallel processing
- num-complex = "0.4.6" - For complex number operations
- faer = "0.22.6" - For high-performance linear algebra operations
- faer-core = "0.17.1" - Core functionality of faer
- faer-ext = { version = "0.6.0", features = ["ndarray"] } - For conversion between faer and ndarray
- serde = { version = "1.0.219", features = ["derive"] } - For serialization/deserialization
- serde_json = "1.0.140" - For JSON handling

## Matrix Types Usage
- Use ndarray or Vec for external APIs and public interfaces
- Use faer matrices for all internal matrix operations and computations
- Use faer-ext for conversions between faer and ndarray at API boundaries
- Take care with version compatibility between faer, faer-core, and faer-ext

## Directory Structure
- **src/**
  - **lib.rs** - Main library entry point and module declarations
  - **input/** - Input file parsing and validation
  - **atoms/** - Atomic data and calculations
  - **potential/** - Potential calculation modules
  - **scattering/** - Scattering calculations
  - **path/** - Path finding and filtering
  - **fms/** - Full multiple scattering
  - **xas/** - X-ray absorption spectroscopy calculations
  - **utils/** - Common utilities and helper functions
  - **cli/** - Command-line interface implementation
- **tests/**
  - Unit and integration tests organized by module
  - Common test fixtures and utilities
- **examples/**
  - Sample applications and usage examples
  - Comparison with FEFF10 reference results
  - Tutorials for common use cases
- **docs/**
  - Additional documentation beyond code comments
  - Physics background and theory
  - Implementation notes and design decisions
- **plans/**
  - Development roadmaps and milestones
  - Feature planning documents
  - Architecture diagrams
- **.states/**
  - Development state tracking
  - Implementation progress records
  - Serialized intermediate calculation states
  - Checkpoint files for long-running calculations

## Code Style Guidelines
- **Formatting**: Follow Rust standard formatting with `cargo fmt`
- **Naming**: Use physics-meaningful names that describe physical quantities
- **Organization**: Split code into modules matching FEFF components (input, potential, scattering, etc.)
- **Documentation**: All public APIs must have rustdoc comments with physics explanation
- **Error Types**: Define specific error types with thiserror for each module
- **Performance**: Optimize numerically intensive operations, use faer efficiently
- **Testing**: Unit test physics correctness against reference values from FEFF10 examples
- **Interfaces**: Create clean abstractions over FEFF10 functionality with idiomatic Rust APIs

## Licensing
- This project uses a hybrid MIT license with FEFF10 attribution requirements
- The code is MIT-licensed with additional conditions for FEFF10 attribution
- Copyright (c) 2025 Ameyanagi
- Ensure all code follows the license requirements by:
  - Including the copyright header in all source files
  - Maintaining attribution to FEFF10.0, University of Washington, and SLAC
  - Providing proper attribution for algorithms and methods
  - Respecting the original FEFF10 license for reference code
- The full FEFF10 license is available in FEFF10_LICENSE
