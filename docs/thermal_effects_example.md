# Thermal Effects Example

This document explains the `thermal_effects_example.rs` example that demonstrates thermal models in the feff-rs library. The example shows how to use various thermal models and their effects on XAS calculations.

## What This Example Demonstrates

The example demonstrates several key features of the thermal modeling capabilities in feff-rs:

1. Creating and comparing different thermal models (Debye, Einstein, Correlated Debye)
2. Calculating Mean-Square Displacements (MSDs) at different temperatures
3. Computing Debye-Waller factors for EXAFS calculations
4. Modeling anisotropic thermal effects for different crystal systems

## Running the Example

You can run the example with the following command:

```
cargo run --example thermal_effects_example
```

The example will create several data files in the current directory:

- `thermal_models_msd.dat`: Shows how MSDs vary with temperature for different models
- `debye_waller_factors.dat`: Shows how Debye-Waller factors vary with wavenumber k
- `anisotropic_thermal.dat`: Shows how thermal effects depend on direction for different crystal systems

## Example Walkthrough

### 1. Creating the Atomic Structure

First, the example creates a simple Cu FCC structure:

```rust
fn create_cu_structure() -> AtomicStructure {
    let mut structure = AtomicStructure::new();
    
    // Add potential types
    let cu_potential = PotentialType::new(0, 29).unwrap();
    structure.add_potential_type(cu_potential);
    
    // Add atoms in an FCC arrangement
    let a = 3.6140; // Cu lattice constant in Å
    
    // Add central atom
    let cu_central = Atom::new(29, Vector3D::new(0.0, 0.0, 0.0), 0).unwrap();
    let central_idx = structure.add_atom(cu_central);
    
    // Add FCC lattice points
    structure.add_atom(Atom::new(29, Vector3D::new(0.0, 0.0, a), 0).unwrap());
    structure.add_atom(Atom::new(29, Vector3D::new(0.0, a, 0.0), 0).unwrap());
    // ... more atoms
    
    // Set central atom
    structure.set_central_atom(central_idx).unwrap();
    
    structure
}
```

### 2. Creating Thermal Parameters and Models

The example creates thermal parameters for different models and then generates thermal models:

```rust
// Create thermal parameters for different models
let debye_params = ThermalParameters::new_debye(300.0, debye_temp);
let einstein_params = ThermalParameters::new_einstein(300.0, einstein_freq);
let corr_debye_params = ThermalParameters::new_correlated_debye(300.0, debye_temp);

// Create the thermal models from the parameters
let debye_model = debye_params.create_model(reduced_mass, Some(path_length));
let einstein_model = einstein_params.create_model(reduced_mass, Some(path_length));
let corr_debye_model = corr_debye_params.create_model(reduced_mass, Some(path_length));
```

### 3. Calculating Mean-Square Displacements

The example calculates MSDs at different temperatures:

```rust
// Calculate MSD at different temperatures from 50K to 900K
for temp in (50..=900).step_by(50) {
    let debye_msd = debye_model.mean_square_displacement(temp as f64);
    let einstein_msd = einstein_model.mean_square_displacement(temp as f64);
    let corr_debye_msd = corr_debye_model.mean_square_displacement(temp as f64);

    writeln!(
        file,
        "{} {:.6} {:.6} {:.6}",
        temp, debye_msd, einstein_msd, corr_debye_msd
    )?;
}
```

### 4. Calculating Debye-Waller Factors

The example calculates Debye-Waller factors for different k values:

```rust
// Calculate for k values from 2 to 15 Å⁻¹
let temp = 300.0; // Room temperature
for k in (20..=150).map(|k| k as f64 / 10.0) {
    let debye_dw = debye_model.debye_waller_factor(temp, k);
    let einstein_dw = einstein_model.debye_waller_factor(temp, k);
    let corr_debye_dw = corr_debye_model.debye_waller_factor(temp, k);

    writeln!(
        file,
        "{:.1} {:.6} {:.6} {:.6}",
        k, debye_dw, einstein_dw, corr_debye_dw
    )?;
}
```

### 5. Modeling Anisotropic Thermal Effects

The example demonstrates anisotropic thermal effects:

```rust
// Demonstrate anisotropic thermal parameters for different crystal systems
let cubic_params = create_anisotropic_thermal_parameters(
    temp, "cubic", debye_temp, None
);
let tetragonal_params = create_anisotropic_thermal_parameters(
    temp, "tetragonal", debye_temp, Some(1.5)
);
let layered_params = create_anisotropic_thermal_parameters(
    temp, "layered", debye_temp, Some(2.5)
);
```

And tests the effects along different directions:

```rust
// Calculate for different directions
let directions = [
    ([1.0, 0.0, 0.0], "x-axis"),
    ([0.0, 1.0, 0.0], "y-axis"),
    ([0.0, 0.0, 1.0], "z-axis"),
    // ... more directions
];

for (dir, name) in &directions {
    // Set up models with these directions
    let mut cubic = cubic_params.clone();
    cubic.path_direction = Some(*dir);
    
    // ... similar for other models
    
    // Calculate Debye-Waller factors
    let cubic_dw = cubic_model.debye_waller_factor(temp, k_value);
    // ... similar for other models
    
    writeln!(file, "{} {:.6} {:.6} {:.6}", name, cubic_dw, tetragonal_dw, layered_dw)?;
}
```

## Understanding the Results

### Mean-Square Displacements

The MSDs in `thermal_models_msd.dat` show:

- MSDs increase with temperature, following the expected physics
- Einstein model typically gives higher MSDs than the Debye model at the same temperature
- Correlated Debye model shows lower values due to correlation effects
- At high temperatures, the relationship becomes more linear (classical limit)

### Debye-Waller Factors

The Debye-Waller factors in `debye_waller_factors.dat` show:

- DW factors decrease exponentially with increasing k
- This explains why EXAFS oscillations are damped at high k values
- Different models produce different damping rates, affecting EXAFS analysis

### Anisotropic Effects

The results in `anisotropic_thermal.dat` show:

- Cubic materials have isotropic behavior (same DW factor in all directions)
- Tetragonal materials show enhanced vibrations along the z-axis
- Layered materials show strongly enhanced vibrations perpendicular to the layers
- These anisotropic effects can significantly impact EXAFS and XANES features in non-cubic materials

## Implementing in Your Own Code

To use thermal models in your XAS calculations:

1. Create thermal parameters with appropriate values for your material
2. If working with non-cubic materials, use anisotropic parameters
3. Create thermal models from these parameters
4. Apply the thermal models in your XAS calculations

Example EXAFS calculation with thermal effects:

```rust
// Create thermal parameters for your material
let thermal_params = ThermalParameters::new_correlated_debye(
    300.0,  // Temperature in K
    315.0   // Debye temperature in K
);

// Set up EXAFS parameters
let mut exafs_params = ExafsParameters {
    edge: Edge::K,
    energy_range: EnergyGrid::new(8979.0, 2.0, 12.0, 0.1),
    // ... other parameters
    thermal_parameters: Some(thermal_params),
};

// Calculate EXAFS with thermal effects
let exafs_data = calculate_exafs(&structure, &exafs_params)?;
```

## Conclusion

This example demonstrates the powerful thermal modeling capabilities in the feff-rs library, which are essential for accurate simulation of XAS spectra, especially for:

- Temperature-dependent studies
- Materials with anisotropic thermal properties
- EXAFS analysis where accurate Debye-Waller factors are crucial

By understanding and applying these thermal models, you can significantly improve the agreement between theoretical calculations and experimental measurements.