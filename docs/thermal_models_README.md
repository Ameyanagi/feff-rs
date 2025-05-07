# Thermal Models in feff-rs

This README provides an overview of the thermal modeling capabilities in the feff-rs library, including implementation details, usage guidance, and physics background.

## Overview

The thermal models in feff-rs account for the effects of atomic vibrations on X-ray absorption spectroscopy (XAS) measurements. These effects are crucial for accurate modeling of experimental data, especially for:

- Temperature-dependent studies
- Materials with distinct vibrational properties (metals vs. covalent compounds)
- Non-cubic materials with anisotropic thermal vibrations
- Accurate EXAFS analysis where thermal disorder is a major factor

## Implemented Models

### Base Models

1. **Debye Model (`DebyeModel`)**: 
   - Classic model for thermal vibrations in crystalline solids
   - Based on a continuous distribution of vibrational modes up to a cutoff frequency
   - Characterized by the Debye temperature (θD)
   - Well-suited for simple metals and highly symmetric materials

2. **Einstein Model (`EinsteinModel`)**:
   - Treats atoms as independent harmonic oscillators with a single frequency
   - Characterized by the Einstein frequency (ωE)
   - Better for materials with localized vibrations or stiff bonds
   - Good approximation for molecular and covalent systems

3. **Correlated Debye Model (`CorrelatedDebyeModel`)**:
   - Enhanced model that accounts for correlation between atomic vibrations
   - Particularly important for EXAFS, where relative displacements matter
   - Includes path distance-dependent correlation factors
   - More accurate for modeling multiple-scattering paths

### Advanced Model

4. **Anisotropic Thermal Model (`AnisotropicThermalModel`)**:
   - Accounts for direction-dependent thermal vibrations
   - Essential for non-cubic materials (tetragonal, hexagonal, orthorhombic, etc.)
   - Combines any base model with directional displacement factors
   - Projects anisotropic vibrations onto specific scattering paths

## Implementation Details

### ThermalModel Trait

All thermal models implement the `ThermalModel` trait, which defines:

- `mean_square_displacement(&self, temperature: f64) -> f64`: Calculates σ² at a given temperature
- `debye_waller_factor(&self, temperature: f64, k: f64) -> f64`: Calculates exp(-2k²σ²)

### Anisotropic Model Implementation

The anisotropic model is implemented using composition:

1. It wraps a base thermal model (Debye, Einstein, or Correlated Debye)
2. Adds displacement factors [u_x, u_y, u_z] for different crystallographic directions
3. Adds path direction information [d_x, d_y, d_z] as a unit vector
4. Projects the anisotropic thermal ellipsoid onto the path direction
5. Scales the base model's mean-square displacement accordingly

This approach maintains the temperature dependence of the base model while adding directional effects.

## Physics of Thermal Vibrations in XAS

### Debye-Waller Factor

The Debye-Waller factor exp(-2k²σ²) accounts for the dampening of XAS signals due to thermal disorder. The mean-square relative displacement (MSRD) σ² has several important properties:

1. **Temperature dependence**:
   - At high T: σ² ∝ T (classical limit)
   - At low T: σ² approaches non-zero value (zero-point motion)
   - Transition region follows the Debye model with T³ dependence

2. **Path length dependence**:
   - For single-scattering: σ² increases with path length
   - For multiple-scattering: σ² behavior is complex due to correlation effects

3. **Directional dependence** (for anisotropic materials):
   - σ² varies with path orientation relative to crystal axes
   - Thermal ellipsoids describe the directional variation

### Physical Meaning of Parameters

- **Debye temperature (θD)**: Characteristic temperature related to the maximum phonon frequency. Higher values indicate stiffer lattices with smaller thermal vibrations.

- **Einstein frequency (ωE)**: Characteristic frequency of localized vibrations. Relates to bond strength.

- **Displacement factors [u_x, u_y, u_z]**: Relative amplitudes of thermal vibrations along different crystallographic axes (dimensionless ratios).

- **Correlation factor**: Measures how atomic motions are related. Range from 0 (independent motion) to 1 (perfectly correlated).

## Usage Guide

### Creating Basic Thermal Parameters

```rust
use feff_rs::xas::thermal::ThermalParameters;

// Debye model (for metals, simple materials)
let debye_params = ThermalParameters::new_debye(
    300.0,  // Temperature in Kelvin
    315.0   // Debye temperature for Cu
);

// Einstein model (for covalent bonds, molecular materials)
let einstein_params = ThermalParameters::new_einstein(
    300.0,  // Temperature in Kelvin
    65.0    // Einstein frequency in meV (typical for metal-oxide bonds)
);

// Correlated Debye model (for accurate EXAFS analysis)
let correlated_params = ThermalParameters::new_correlated_debye(
    300.0,  // Temperature in Kelvin
    315.0   // Debye temperature
);
```

### Creating Anisotropic Thermal Parameters

```rust
use feff_rs::xas::thermal::{ThermalParameters, create_anisotropic_thermal_parameters};

// Manual specification
let anisotropic_params = ThermalParameters::new_anisotropic_debye(
    300.0,            // Temperature in Kelvin
    315.0,            // Debye temperature
    [1.0, 1.0, 1.8],  // Enhanced vibrations along z-axis
    Some([0.0, 0.0, 1.0]) // Path along z-axis
);

// Using crystal system factory function
let tetragonal_params = create_anisotropic_thermal_parameters(
    300.0,         // Temperature in Kelvin
    "tetragonal",  // Crystal system
    315.0,         // Debye temperature
    Some(1.5)      // Anisotropy ratio
);
```

### Using Thermal Parameters in XAS Calculations

```rust
use feff_rs::xas::exafs::calculate_exafs;

// Calculate EXAFS with thermal effects
let exafs = calculate_exafs(
    &structure,           // Atomic structure
    &paths,               // Scattering paths
    &energies,            // Energy grid
    Some(&thermal_params) // Thermal parameters
)?;
```

## Typical Values for Common Materials

### Debye Temperatures (θD)

| Material | Debye Temperature (K) |
|----------|----------------------|
| Cu       | 315                  |
| Fe       | 470                  |
| Ni       | 375                  |
| Al       | 390                  |
| Si       | 640                  |
| MgO      | 750                  |
| TiO2     | 770                  |

### Einstein Frequencies (ωE)

| Bond Type     | Einstein Frequency (meV) |
|--------------|-----------------------|
| Metal-metal   | 15-25                |
| Metal-oxide   | 40-70                |
| Metal-halide  | 20-40                |
| Si-O          | 100-140              |
| C-C           | 160-180              |
| O-H           | 450-500              |

### Anisotropy Ratios

| Material Type   | Anisotropy Ratio | Description |
|----------------|------------------|-------------|
| Cubic metals    | 1.0              | Isotropic   |
| Tetragonal      | 1.3-1.7          | Enhanced c-axis vibrations |
| Hexagonal       | 1.5-2.0          | Enhanced c-axis vibrations |
| Layered (e.g., graphite) | 2.0-3.0 | Much larger out-of-plane vibrations |
| Chain structures | 1.7-2.5         | Enhanced vibrations perpendicular to chains |

## References

1. Sevillano, E., Meuth, H. & Rehr, J. J. Extended x-ray absorption fine structure Debye-Waller factors. I. Monatomic crystals. Phys. Rev. B 20, 4908–4911 (1979)

2. Beni, G. & Platzman, P. M. Temperature and polarization dependence of extended x-ray absorption fine-structure spectra. Phys. Rev. B 14, 1514–1518 (1976)

3. Fornasini, P. et al. Extended X-ray-absorption fine-structure measurements of copper: Local dynamics, anharmonicity, and thermal expansion. Phys. Rev. B 70, 174301 (2004)

4. Vila, F. D., Rehr, J. J., Rossner, H. H. & Krappe, H. J. Theoretical x-ray absorption Debye-Waller factors. Phys. Rev. B 76, 014301 (2007)