# FEFF-rs Architecture Design

This document outlines the architecture for FEFF-rs, a Rust implementation of the FEFF code for X-ray absorption spectroscopy calculations.

## Overview

FEFF is a sophisticated code that calculates X-ray Absorption Spectra (XAS) using the real-space multiple-scattering approach. The key components include:

1. Input parsing
2. Atomic potential calculation
3. Scattering path calculation
4. Multiple scattering calculation
5. XAS spectrum generation

Our Rust implementation will follow a modular design that separates these concerns while allowing for efficient computation.

## Core Components

### 1. Input Module

**Purpose**: Parse and validate input files in FEFF format

**Components**:
- `InputParser`: Parses FEFF input files
- `AtomicStructure`: Represents atomic coordinates and potentials
- `CalculationParameters`: Stores configuration for calculations
- `CIFImporter`: Imports crystal information files

**Key Functionality**:
- Parse FEFF input format
- Support CIF crystal structure format
- Validate input parameters
- Generate calculation settings

### 2. Atomic Module

**Purpose**: Handle atomic data and calculations

**Components**:
- `AtomicDatabase`: Data for elements (energy levels, occupations)
- `CoreHole`: Core hole effects calculator
- `ElectronDensity`: Electron density calculator

**Key Functionality**:
- Store/retrieve atomic properties
- Calculate electronic configurations
- Handle core hole effects

### 3. Potential Module

**Purpose**: Calculate atomic potentials

**Components**:
- `MuffinTinPotential`: Muffin-tin potential calculator
- `ExchangeCorrelation`: Exchange-correlation functionals
- `SelfConsistency`: Self-consistent field (SCF) solver
- `PotentialOverlap`: Overlapping potentials handler

**Key Functionality**:
- Calculate muffin-tin potentials
- Support multiple exchange-correlation models
- Handle self-consistency iterations
- Calculate phase shifts for scattering

### 4. Scattering Module

**Purpose**: Calculate scattering matrices and properties

**Components**:
- `PhaseShift`: Phase shift calculator
- `ScatteringMatrix`: Single-atom scattering matrix
- `PropagatorMatrix`: Free-electron propagator
- `KSpaceScattering`: K-space scattering methods

**Key Functionality**:
- Calculate scattering matrices
- Generate propagator matrices
- Handle relativistic effects
- Support energy-dependent calculations

### 5. Path Module

**Purpose**: Find and filter scattering paths

**Components**:
- `PathFinder`: Enumerate possible scattering paths
- `PathFilter`: Filter paths by importance
- `PathCalculator`: Calculate path contributions

**Key Functionality**:
- Generate paths up to specified cutoff
- Filter paths based on importance
- Group similar paths
- Calculate path-specific parameters

### 6. FMS Module (Full Multiple Scattering)

**Purpose**: Perform full multiple scattering calculations

**Components**:
- `FMSMatrix`: Full multiple scattering matrix builder
- `FMSSolver`: Matrix equation solver
- `KSpaceFMS`: K-space full multiple scattering

**Key Functionality**:
- Build FMS matrices efficiently
- Solve matrix equations for XAS
- Support both real-space and k-space approaches
- Parallelize calculations

### 7. XAS Module

**Purpose**: Generate XAS spectra

**Components**:
- `XANESCalculator`: XANES spectrum generator
- `EXAFSCalculator`: EXAFS spectrum generator
- `SpectrumProcessor`: Post-processing for spectra

**Key Functionality**:
- Calculate XANES (near-edge) spectra
- Calculate EXAFS (extended) spectra
- Apply broadening and other corrections
- Generate output files

### 8. Utils Module

**Purpose**: Provide common utilities

**Components**:
- `Math`: Mathematical functions
- `LinearAlgebra`: Matrix operations using faer
- `Parallel`: Parallelization utilities
- `Converters`: Format conversion utilities

**Key Functionality**:
- Matrix operations optimized for scattering
- Complex number operations
- Coordinate transformations
- Serialization/deserialization

### 9. CLI Module

**Purpose**: Command-line interface

**Components**:
- `Commands`: Command implementations
- `Options`: Option parsing
- `Output`: Output handling

**Key Functionality**:
- Parse command-line arguments
- Run calculations based on commands
- Format and output results

## Data Flow

1. Input parsing → Atomic structure and calculation parameters
2. Atomic structure → Potential calculation
3. Potentials → Scattering matrices
4. Atomic structure + Scattering matrices → Path enumeration and filtering
5. Paths + Scattering matrices → XAS calculation
   - Path-by-path approach for EXAFS
   - Full matrix for XANES
6. XAS calculation → Spectrum generation and output

## Key Design Decisions

### 1. Matrix Storage and Computation

- Use `faer` for all internal matrix operations
- Use `ndarray` or `Vec` for external APIs
- Employ `faer-ext` for conversions between matrix types

### 2. Parallelization Strategy

- Use `rayon` for data parallelism
- Parallelize path calculations (embarrassingly parallel)
- Parallelize matrix operations (through faer)
- Support concurrent potential calculations

### 3. Memory Management

- Use memory-efficient representations for sparse matrices
- Implement chunk-based processing for large structures
- Allow configurable memory limits
- Provide checkpointing for long calculations

### 4. Error Handling

- Use `thiserror` for library error types
- Use `anyhow` for application error handling
- Implement graceful degradation where possible
- Provide detailed error messages with physics context

### 5. Configuration and Extensibility

- Make all physical parameters configurable
- Support multiple theory levels
- Design extensible interfaces for new spectroscopies
- Implement plugin architecture for exchange-correlation functionals

## Implementation Priorities

1. Basic input parsing and atomic structure representation
2. Potential calculation with simple models
3. Path enumeration and filtering
4. Simple EXAFS calculation using path expansion
5. FMS implementation for XANES
6. Advanced features (thermal effects, self-consistency, etc.)
7. Performance optimization and parallelization
8. Additional spectroscopies (EELS, XES, etc.)

This architecture provides a solid foundation for implementing FEFF's functionality in Rust while allowing for modern software engineering practices, performance optimizations, and future extensions.