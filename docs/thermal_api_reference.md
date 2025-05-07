# Thermal Effects API Reference

This document provides a detailed API reference for the thermal effects modeling capabilities in FEFF-rs.

## Overview

FEFF-rs provides extensive support for modeling thermal effects in X-ray Absorption Spectroscopy (XAS) through two key modules:

1. `utils::thermal`: Core implementations of thermal models
2. `xas::thermal`: High-level thermal parameters for XAS calculations

## Core Thermal Models (`utils::thermal`)

### `ThermalModel` Trait

The foundation of thermal modeling is the `ThermalModel` trait:

```rust
pub trait ThermalModel {
    /// Calculate the mean-square displacement (σ²) at a given temperature
    fn mean_square_displacement(&self, temperature: f64) -> f64;

    /// Calculate the Debye-Waller factor exp(-2k²σ²) at a given temperature and wavenumber
    fn debye_waller_factor(&self, temperature: f64, k: f64) -> f64;
}
```

This trait is implemented by various thermal models with different physical approaches.

### `DebyeModel`

The standard Debye model for thermal vibrations in solids.

```rust
pub struct DebyeModel {
    /// Debye temperature in Kelvin
    debye_temperature: f64,
    /// Reduced mass in atomic mass units (amu)
    reduced_mass: f64,
}

impl DebyeModel {
    /// Create a new Debye model
    pub fn new(debye_temperature: f64, reduced_mass: f64) -> Self;
}
```

**Physical meaning**: The Debye model treats the crystal lattice as a continuous medium with a phonon density of states proportional to ω² up to a maximum frequency determined by the Debye temperature. It works well for simple cubic systems and metals.

### `EinsteinModel`

The Einstein model for localized vibrations.

```rust
pub struct EinsteinModel {
    /// Einstein frequency in meV
    einstein_frequency: f64,
    /// Reduced mass in atomic mass units (amu)
    reduced_mass: f64,
}

impl EinsteinModel {
    /// Create a new Einstein model
    pub fn new(einstein_frequency: f64, reduced_mass: f64) -> Self;
}
```

**Physical meaning**: The Einstein model treats atoms as independent harmonic oscillators with a single characteristic frequency. It works well for systems with localized vibrations, like covalent bonds and molecular crystals.

### `CorrelatedDebyeModel`

An enhanced Debye model that accounts for correlation between atomic vibrations.

```rust
pub struct CorrelatedDebyeModel {
    /// Debye temperature in Kelvin
    debye_temperature: f64,
    /// Reduced mass in atomic mass units (amu)
    reduced_mass: f64,
    /// Path distance in Å - affects correlation
    path_distance: f64,
    /// Correlation factor (0.0-1.0)
    correlation: f64,
}

impl CorrelatedDebyeModel {
    /// Create a new Correlated Debye model
    pub fn new(debye_temperature: f64, reduced_mass: f64, path_distance: f64) -> Self;

    /// Create a new Correlated Debye model with a specified correlation factor
    pub fn with_correlation(
        debye_temperature: f64,
        reduced_mass: f64,
        path_distance: f64,
        correlation: f64,
    ) -> Self;
}
```

**Physical meaning**: The Correlated Debye model accounts for the fact that atoms connected by chemical bonds tend to move in a coordinated manner. This correlation reduces the effective mean-square relative displacement (MSRD) compared to what would be predicted if atoms moved independently. The correlation is typically stronger for shorter paths (nearest neighbors) and weaker for longer paths.

### `AnisotropicThermalModel`

Advanced model for direction-dependent thermal vibrations.

```rust
pub struct AnisotropicThermalModel {
    /// Base thermal model used for overall magnitude of vibrations
    base_model: Box<dyn ThermalModel>,
    /// Direction-dependent scaling factors for thermal displacements
    displacement_factors: [f64; 3],
    /// Path direction in Cartesian coordinates (unit vector)
    path_direction: [f64; 3],
}

impl AnisotropicThermalModel {
    /// Create a new anisotropic thermal model based on a Debye model
    pub fn new_from_debye(
        debye_temperature: f64,
        reduced_mass: f64,
        displacement_factors: [f64; 3],
        path_direction: [f64; 3],
    ) -> Self;

    /// Create a new anisotropic thermal model based on a Correlated Debye model
    pub fn new_from_correlated_debye(
        debye_temperature: f64,
        reduced_mass: f64,
        path_distance: f64,
        displacement_factors: [f64; 3],
        path_direction: [f64; 3],
    ) -> Self;

    /// Create a new anisotropic thermal model based on an Einstein model
    pub fn new_from_einstein(
        einstein_frequency: f64,
        reduced_mass: f64,
        displacement_factors: [f64; 3],
        path_direction: [f64; 3],
    ) -> Self;
}
```

**Physical meaning**: The anisotropic model accounts for the fact that in non-cubic materials, thermal vibrations have different amplitudes along different crystallographic directions. The model applies direction-dependent scaling factors to a base thermal model (Debye, Einstein, or Correlated Debye) and projects these anisotropic vibrations onto specific scattering paths.

### Factory Functions

The module provides several factory functions to simplify the creation of thermal models:

```rust
/// Create a thermal model from parameters
pub fn create_thermal_model(
    model_type: &str,
    debye_temperature: Option<f64>,
    einstein_frequency: Option<f64>,
    reduced_mass: f64,
) -> Box<dyn ThermalModel>;

/// Create a thermal model with path distance for correlated models
pub fn create_thermal_model_with_path(
    model_type: &str,
    debye_temperature: Option<f64>,
    einstein_frequency: Option<f64>,
    reduced_mass: f64,
    path_distance: f64,
) -> Box<dyn ThermalModel>;

/// Create an anisotropic thermal model with crystal direction information
pub fn create_anisotropic_thermal_model(
    base_model_type: &str,
    debye_temperature: Option<f64>,
    einstein_frequency: Option<f64>,
    reduced_mass: f64,
    path_distance: Option<f64>,
    displacement_factors: [f64; 3],
    path_direction: [f64; 3],
) -> Box<dyn ThermalModel>;
```

## XAS Thermal Parameters (`xas::thermal`)

This module provides high-level structures for configuring thermal effects in XAS calculations.

### `ThermalParameters`

The main structure for thermal configuration in XAS calculations.

```rust
pub struct ThermalParameters {
    /// Temperature in Kelvin
    pub temperature: f64,
    /// Type of thermal model: "debye", "einstein", "correlated_debye", or "anisotropic"
    pub model_type: String,
    /// Debye temperature in Kelvin (for Debye model)
    pub debye_temperature: f64,
    /// Einstein frequency in meV (for Einstein model)
    pub einstein_frequency: Option<f64>,
    /// Pair-specific thermal parameters, keyed by element pair
    pub pair_parameters: Option<HashMap<PairKey, PairThermalParameters>>,
    /// Crystal direction-dependent displacement factors [u_x, u_y, u_z] (for anisotropic model)
    pub displacement_factors: Option<[f64; 3]>,
    /// Path direction in crystal coordinates [d_x, d_y, d_z] (for anisotropic model)
    pub path_direction: Option<[f64; 3]>,
}

impl ThermalParameters {
    /// Create new thermal parameters for a Debye model
    pub fn new_debye(temperature: f64, debye_temperature: f64) -> Self;

    /// Create new thermal parameters for an Einstein model
    pub fn new_einstein(temperature: f64, einstein_frequency: f64) -> Self;

    /// Create new thermal parameters for a Correlated Debye model
    pub fn new_correlated_debye(temperature: f64, debye_temperature: f64) -> Self;

    /// Create new thermal parameters for an anisotropic Debye model
    pub fn new_anisotropic_debye(
        temperature: f64,
        debye_temperature: f64,
        displacement_factors: [f64; 3],
        path_direction: Option<[f64; 3]>,
    ) -> Self;

    /// Create new thermal parameters for an anisotropic Einstein model
    pub fn new_anisotropic_einstein(
        temperature: f64,
        einstein_frequency: f64,
        displacement_factors: [f64; 3],
        path_direction: Option<[f64; 3]>,
    ) -> Self;

    /// Create new thermal parameters for an anisotropic correlated Debye model
    pub fn new_anisotropic_correlated_debye(
        temperature: f64,
        debye_temperature: f64,
        displacement_factors: [f64; 3],
        path_direction: Option<[f64; 3]>,
    ) -> Self;

    /// Add pair-specific thermal parameters for a specific element pair
    pub fn with_pair_parameters(
        self,
        z1: i32,
        z2: i32,
        pair_params: PairThermalParameters,
    ) -> Self;

    /// Get thermal parameters for a specific element pair
    pub fn get_pair_parameters(&self, z1: i32, z2: i32) -> Option<&PairThermalParameters>;

    /// Create a thermal model from these parameters
    pub fn create_model(
        &self,
        reduced_mass: f64,
        path_distance: Option<f64>,
    ) -> Box<dyn crate::utils::thermal::ThermalModel>;
}
```

### `PairThermalParameters` and `PairKey`

Structures for defining bond-specific thermal parameters:

```rust
/// Pair key for identifying specific atom pairs
pub struct PairKey {
    /// Atomic number of the first element (usually the absorber)
    pub z1: u32,
    /// Atomic number of the second element (usually the scatterer)
    pub z2: u32,
}

/// Bond-specific thermal parameters for a specific pair of elements
pub struct PairThermalParameters {
    /// Type of thermal model: "debye" or "einstein"
    pub model_type: String,
    /// Debye temperature in Kelvin (for Debye model)
    pub debye_temperature: f64,
    /// Einstein frequency in meV (for Einstein model)
    pub einstein_frequency: Option<f64>,
    /// Optional description of this bond type
    pub description: Option<String>,
}

impl PairThermalParameters {
    /// Create new pair thermal parameters for a Debye model
    pub fn new_debye(debye_temperature: f64, description: Option<String>) -> Self;

    /// Create new pair thermal parameters for an Einstein model
    pub fn new_einstein(einstein_frequency: f64, description: Option<String>) -> Self;
}
```

### Convenience Factory Functions

Functions for creating thermal parameters for specific materials and scenarios:

```rust
/// Create thermal parameters with common element pairs for typical materials
pub fn create_standard_thermal_parameters(temperature: f64) -> ThermalParameters;

/// Create anisotropic thermal parameters for different crystal systems
pub fn create_anisotropic_thermal_parameters(
    temperature: f64,
    crystal_system: &str,
    debye_temperature: f64,
    anisotropy_ratio: Option<f64>,
) -> ThermalParameters;

/// Create material-specific thermal parameters
pub fn create_material_thermal_parameters(
    material: &str,
    temperature: f64,
    custom_debye_temp: Option<f64>,
) -> ThermalParameters;

/// Create thermal parameters optimized for EXAFS analysis
pub fn create_exafs_thermal_parameters(
    material_type: &str,
    temperature: f64,
    debye_temperature: f64,
    include_anharmonic: bool,
) -> ThermalParameters;

/// Create thermal parameters optimized for XANES analysis
pub fn create_xanes_thermal_parameters(
    material_type: &str,
    temperature: f64,
    debye_temperature: f64,
    edge: &str,
) -> ThermalParameters;
```

## Usage Patterns

### Basic Thermal Parameter Creation

```rust
// For a simple cubic material like copper
let cu_params = ThermalParameters::new_debye(300.0, 315.0);

// For a molecular system with localized vibrations
let covalent_params = ThermalParameters::new_einstein(300.0, 70.0);

// For accurate EXAFS analysis
let exafs_params = ThermalParameters::new_correlated_debye(300.0, 315.0);
```

### Anisotropic Thermal Parameters

```rust
// For a tetragonal material like rutile (TiO2)
let tetragonal_params = create_anisotropic_thermal_parameters(
    300.0,          // Temperature in K
    "tetragonal",   // Crystal system
    450.0,          // Debye temperature for TiO2
    Some(1.5)       // Moderate anisotropy
);

// For a highly anisotropic layered material like graphite
let graphite_params = create_anisotropic_thermal_parameters(
    300.0,          // Temperature in K
    "layered",      // Layered material
    420.0,          // Debye temperature for graphite
    Some(3.0)       // Very high anisotropy
);

// Manually specifying displacement factors
let manual_params = ThermalParameters::new_anisotropic_debye(
    300.0,              // Temperature in K
    420.0,              // Debye temperature
    [1.0, 1.0, 2.5],    // Much larger vibrations along z-axis
    Some([0.0, 0.0, 1.0]) // Path along z-axis
);
```

### Material-Specific Parameters

```rust
// For common materials, use the factory function
let copper_params = create_material_thermal_parameters(
    "Cu",     // Material name
    300.0,    // Temperature in K
    None      // Use default Debye temperature
);

let iron_oxide_params = create_material_thermal_parameters(
    "Fe2O3",  // Material name
    300.0,    // Temperature in K
    None      // Use default Debye temperature
);
```

### Integration with XAS Calculations

```rust
use feff_rs::xas::{calculate_xanes, XanesParameters, Edge};
use feff_rs::xas::thermal::ThermalParameters;

// Create thermal parameters
let thermal_params = ThermalParameters::new_debye(300.0, 315.0);

// Include in XANES parameters
let params = XanesParameters {
    edge: Edge::K,
    energy_range: (-10.0, 50.0, 0.5),
    thermal_parameters: Some(thermal_params),
    ..Default::default()
};

// Calculate XANES with thermal effects
let spectrum = calculate_xanes(&structure, &params).unwrap();
```

### Bond-Specific Thermal Parameters

```rust
use feff_rs::xas::thermal::{ThermalParameters, PairThermalParameters};
use std::collections::HashMap;

// Start with base parameters for the material
let mut params = ThermalParameters::new_debye(300.0, 400.0);

// Add specific parameters for Fe-O bonds
params = params.with_pair_parameters(
    26, // Fe
    8,  // O
    PairThermalParameters::new_einstein(
        65.0, // Einstein frequency
        Some("Fe-O bond".to_string())
    )
);

// Add specific parameters for Fe-Fe bonds
params = params.with_pair_parameters(
    26, // Fe
    26, // Fe
    PairThermalParameters::new_debye(
        470.0, // Debye temperature
        Some("Fe-Fe bond".to_string())
    )
);
```

### Temperature-Dependent Studies

For studying how XAS changes with temperature:

```rust
let temperatures = [100.0, 200.0, 300.0, 400.0, 500.0];
let mut spectra = Vec::new();

for temp in &temperatures {
    // Create thermal parameters for this temperature
    let thermal_params = ThermalParameters::new_debye(*temp, 315.0);
    
    // Set up calculation parameters
    let params = XanesParameters {
        thermal_parameters: Some(thermal_params),
        ..Default::default()
    };
    
    // Calculate spectrum and store
    let spectrum = calculate_xanes(&structure, &params).unwrap();
    spectra.push(spectrum);
}
```

## Physical Parameters Reference

### Debye Temperatures

| Material      | Debye Temperature (K) |
|---------------|----------------------:|
| Cu (copper)   | 315                   |
| Fe (iron)     | 470                   |
| Ni (nickel)   | 375                   |
| Al (aluminum) | 390                   |
| Ag (silver)   | 225                   |
| Au (gold)     | 170                   |
| Pt (platinum) | 230                   |
| Pb (lead)     | 95                    |
| Si (silicon)  | 645                   |
| Ge (germanium)| 360                   |
| C (diamond)   | 2230                  |
| MgO           | 750                   |
| NaCl          | 320                   |
| TiO2 (rutile) | 450                   |
| Fe2O3         | 500                   |
| Graphite      | 420                   |
| MoS2          | 380                   |

### Einstein Frequencies

| Bond Type     | Einstein Frequency (meV) |
|---------------|-------------------------:|
| C-C (diamond) | 160                      |
| Si-O          | 130                      |
| Cu-O          | 70                       |
| Fe-O          | 65                       |
| Ni-O          | 68                       |
| Si-Si         | 60                       |
| Al-O          | 80                       |
| Mo-S          | 60                       |

### Anisotropy Ratios

| Material Type    | Typical Anisotropy Ratio |
|------------------|-------------------------:|
| Cubic            | 1.0 (isotropic)          |
| Tetragonal       | 1.3-1.7                  |
| Hexagonal        | 1.5-2.0                  |
| Orthorhombic     | 1.3-2.0 (multiple axes)  |
| Layered materials| 2.0-4.0                  |
| Chain compounds  | 1.8-3.0                  |

## Best Practices

1. **Choose the appropriate model** for your material:
   - Use the Debye model for simple cubic systems and metals
   - Use the Einstein model for systems with localized vibrations
   - Use the Correlated Debye model for EXAFS analysis
   - Use anisotropic models for non-cubic materials

2. **Use physically reasonable parameters**:
   - For Debye temperatures, refer to literature values
   - For anisotropic materials, ensure the displacement factors match crystal symmetry
   - For path-specific calculations, verify path directions

3. **Check temperature dependence**:
   - Verify that the thermal model produces physically reasonable behavior across temperatures
   - The Debye-Waller factor should increase with temperature
   - Zero-point motion should be present even at T=0K

4. **Consider material-specific factors**:
   - For complex materials, use bond-specific thermal parameters
   - For layered materials, ensure proper modeling of in-plane vs. out-of-plane vibrations
   - For high temperatures, consider anharmonic effects

## References

1. E.D. Crozier, J.J. Rehr, R. Ingalls, "Temperature Dependence of X-Ray Absorption in EXAFS," in X-Ray Absorption: Principles, Applications, Techniques of EXAFS, SEXAFS, and XANES, Wiley, 1988.

2. G. Beni, P.M. Platzman, "Temperature and polarization dependence of extended x-ray absorption fine-structure spectra," Phys. Rev. B 14 (1976) 1514.

3. A.I. Frenkel, J.J. Rehr, "Thermal expansion and x-ray-absorption fine-structure cumulants," Phys. Rev. B 48 (1993) 585-588.

4. N.W. Ashcroft, N.D. Mermin, Solid State Physics, Holt, Rinehart and Winston, New York, 1976.

5. H. Shigematsu, T. Sato, H. Kondo, T. Takeda, "Investigation of the Debye–Waller factor in transition metal alloys and compounds," Journal of Physics F: Metal Physics 7 (1977) 2089-2101.