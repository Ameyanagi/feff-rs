# Thermal Effects in X-ray Absorption Spectroscopy

This document explains how thermal effects are modeled in the feff-rs library and provides guidance for using different thermal models for various applications.

## Overview of Thermal Effects

Temperature affects X-ray absorption spectroscopy (XAS) in several important ways:

1. **Debye-Waller factor damping**: Thermal vibrations cause a reduction in the amplitude of EXAFS oscillations, modeled by the factor exp(-2k²σ²), where σ² is the mean-square relative displacement (MSRD).

2. **Peak broadening**: Thermal effects broaden XANES features, reducing spectral resolution as temperature increases.

3. **Direction-dependent effects**: In anisotropic materials (non-cubic crystals), thermal vibrations have different amplitudes along different crystallographic directions.

4. **Correlation effects**: In multi-atom scattering paths, atomic motions are often correlated, which affects the effective MSRD for longer paths.

## Thermal Models in feff-rs

The feff-rs library provides several thermal models with increasing levels of sophistication:

### 1. Debye Model

The simplest and most widely used model for describing thermal vibrations in crystalline solids.

- **Appropriate for**: Simple cubic materials, metals, and compounds with uniform bonding
- **Key parameter**: Debye temperature (θD), typically 250-500K for most materials
- **Usage example**:

```rust
use feff_rs::xas::thermal::ThermalParameters;

// Create thermal parameters for copper at room temperature
let params = ThermalParameters::new_debye(300.0, 315.0); // 315K is θD for Cu
```

### 2. Einstein Model

Models atoms as independent harmonic oscillators with a characteristic frequency.

- **Appropriate for**: Materials with localized vibrations, stiff bonds, molecular crystals
- **Key parameter**: Einstein frequency (typically 10-100 meV)
- **Usage example**:

```rust
use feff_rs::xas::thermal::ThermalParameters;

// Create thermal parameters for a covalent material (Si-O bond)
let params = ThermalParameters::new_einstein(300.0, 65.0); // 65 meV Einstein frequency
```

### 3. Correlated Debye Model

An enhanced model that accounts for correlation between atomic vibrations.

- **Appropriate for**: Accurate EXAFS analysis, especially for longer paths
- **Key parameters**: Debye temperature, path distance (for correlation calculation)
- **Usage example**:

```rust
use feff_rs::xas::thermal::ThermalParameters;

// Create thermal parameters for iron at room temperature
let params = ThermalParameters::new_correlated_debye(300.0, 470.0); // 470K is θD for Fe
```

### 4. Anisotropic Thermal Models

Advanced models that account for direction-dependent thermal vibrations.

- **Appropriate for**: Non-cubic materials (tetragonal, hexagonal, orthorhombic, etc.)
- **Key parameters**: Base model parameters, displacement factors, path direction
- **Usage examples**:

```rust
use feff_rs::xas::thermal::{ThermalParameters, create_anisotropic_thermal_parameters};

// Manual specification for a tetragonal material
let tetragonal_params = ThermalParameters::new_anisotropic_debye(
    300.0,            // Temperature in Kelvin
    420.0,            // Debye temperature
    [1.0, 1.0, 1.5],  // Larger vibrations along the c-axis (z)
    Some([0.0, 0.0, 1.0]) // Path along z-axis
);

// Using the crystal system factory function
let hexagonal_params = create_anisotropic_thermal_parameters(
    300.0,         // Temperature in Kelvin
    "hexagonal",   // Crystal system
    450.0,         // Debye temperature
    Some(1.8)      // Anisotropy ratio (higher = more anisotropic)
);
```

## Crystal Systems and Thermal Anisotropy

Different crystal systems have characteristic patterns of thermal anisotropy:

| Crystal System | Displacement Factors Pattern | Typical Examples |
|----------------|------------------------------|------------------|
| Cubic          | [1, 1, 1] (isotropic)        | NaCl, Cu, Au     |
| Tetragonal     | [1, 1, r] (r > 1)            | TiO₂, CaSO₄      |
| Hexagonal      | [1, 1, r] (r > 1)            | Zn, graphite     |
| Orthorhombic   | [1, r₁, r₂] (all different)  | MgSO₄, aragonite |
| Monoclinic     | [r₁, r₂, r₃] (all different) | gypsum           |
| Triclinic      | [r₁, r₂, r₃] (all different) | albite           |
| Layered        | [1, 1, r] (r >> 1)           | graphite, MoS₂   |
| Chain          | [r, r, 1] (r > 1)            | polymers         |

The `create_anisotropic_thermal_parameters` function automatically sets appropriate displacement factors based on the crystal system.

## Advanced Use Cases

### Temperature-Dependent Studies

To model how XAS changes with temperature:

```rust
use feff_rs::xas::thermal::ThermalParameters;
use feff_rs::xas::exafs::calculate_exafs_with_thermal_effects;

// Create an array of temperatures
let temperatures = [100.0, 200.0, 300.0, 400.0, 500.0];
let mut spectra = Vec::new();

for temp in temperatures.iter() {
    let thermal_params = ThermalParameters::new_debye(*temp, 315.0);
    let spectrum = calculate_exafs_with_thermal_effects(&structure, &path_list, &thermal_params);
    spectra.push(spectrum);
}
```

### Bond-Specific Thermal Parameters

For materials with different types of bonds that have different vibrational characteristics:

```rust
use feff_rs::xas::thermal::{ThermalParameters, PairThermalParameters, PairKey};
use std::collections::HashMap;

// Start with base parameters
let mut params = ThermalParameters::new_debye(300.0, 350.0);

// Add parameters for specific bonds
let mut pair_params = HashMap::new();

// Metal-oxygen bond (stiffer than metal-metal)
pair_params.insert(
    PairKey::new(28, 8),  // Ni-O pair
    PairThermalParameters::new_einstein(65.0, Some("Ni-O bond".to_string()))
);

// Metal-metal bond
pair_params.insert(
    PairKey::new(28, 28), // Ni-Ni pair
    PairThermalParameters::new_debye(375.0, Some("Ni-Ni bond".to_string()))
);

params.pair_parameters = Some(pair_params);
```

### Anisotropic Effects in Layered Materials

For studying how thermal effects differ for in-plane vs. out-of-plane paths in a layered material:

```rust
use feff_rs::xas::thermal::{ThermalParameters, create_anisotropic_thermal_parameters};

// Create thermal parameters for a layered material like MoS₂
let layered_params = create_anisotropic_thermal_parameters(
    300.0,            // Room temperature
    "layered",        // Layered material type
    450.0,            // Debye temperature
    Some(3.0)         // High anisotropy ratio (3x higher vibrations perpendicular to layers)
);

// For an in-plane Mo-S path
let in_plane_params = layered_params.clone();
in_plane_params.path_direction = Some([1.0, 0.0, 0.0]); // x-direction (in-plane)

// For an out-of-plane Mo-S path
let out_of_plane_params = layered_params.clone();
out_of_plane_params.path_direction = Some([0.0, 0.0, 1.0]); // z-direction (out-of-plane)
```

## Best Practices

1. **Choose the appropriate model** for your material:
   - Use the Debye model for simple cubic systems
   - Use the Einstein model for molecular systems with localized vibrations
   - Use the correlated Debye model for accurate EXAFS analysis
   - Use anisotropic models for non-cubic materials

2. **Use literature values** for Debye temperatures and Einstein frequencies when available.

3. **For anisotropic materials**, be sure to specify the correct crystal system and path directions.

4. **Check temperature dependence** of your results to ensure the thermal model is behaving as expected.

5. **Consider pair-specific parameters** for materials with different types of bonds.

## References

1. E.D. Crozier, J.J. Rehr, R. Ingalls, in: D.C. Koningsberger, R. Prins (Eds.), X-Ray Absorption: Principles, Applications, Techniques of EXAFS, SEXAFS, and XANES, Wiley, New York, 1988, p. 373.

2. G. Beni, P.M. Platzman, Phys. Rev. B 14 (1976) 1514.

3. N.W. Ashcroft, N.D. Mermin, Solid State Physics, Holt, Rinehart and Winston, New York, 1976.

4. E. Sevillano, H. Meuth, J.J. Rehr, Phys. Rev. B 20 (1979) 4908.

5. A.I. Frenkel, J.J. Rehr, Phys. Rev. B 48 (1993) 585.