# Implementation Plan for Temperature-Dependent Atomic Potentials

## Overview

Temperature affects atomic potentials in several ways that are important for accurate X-ray absorption spectroscopy (XAS) calculations:

1. **Thermal expansion** changes interatomic distances and coordination environments
2. **Thermal disorder** causes smearing of potential features due to atomic vibrations
3. **Electronic structure changes** from temperature-dependent thermal excitation of electrons
4. **Phonon effects** modify the effective potential through electron-phonon coupling

The current feff-rs implementation accounts for thermal effects primarily through the Debye-Waller factor in the scattering paths calculation, but does not directly modify the atomic potentials based on temperature. This implementation plan outlines how to extend the atomic potential calculation to include temperature-dependent effects.

## Implementation Strategy

### 1. Extend `MuffinTinPotential` to Accept Temperature Parameters

First, we'll modify the `MuffinTinPotential` struct to include temperature parameters:

```rust
pub struct MuffinTinPotential<'a> {
    // Existing fields
    
    /// Temperature in Kelvin
    temperature: f64,
    
    /// Thermal model type (Debye, Einstein, etc.)
    thermal_model_type: Option<String>,
    
    /// Thermal model parameter (Debye temperature or Einstein frequency)
    thermal_model_parameter: Option<f64>,
}
```

With corresponding methods:

```rust
impl<'a> MuffinTinPotential<'a> {
    // Existing methods
    
    /// Set the temperature for potential calculation
    pub fn set_temperature(&mut self, temperature: f64) -> &mut Self {
        self.temperature = temperature;
        self
    }
    
    /// Set the thermal model parameters
    pub fn set_thermal_model(&mut self, model_type: &str, parameter: f64) -> Result<&mut Self> {
        // Validate model type
        if !["debye", "einstein", "correlated_debye"].contains(&model_type.to_lowercase().as_str()) {
            return Err(PotentialError::Generic(format!(
                "Invalid thermal model type: {}", model_type
            )));
        }
        
        self.thermal_model_type = Some(model_type.to_lowercase());
        self.thermal_model_parameter = Some(parameter);
        Ok(self)
    }
}
```

### 2. Implement Temperature-Dependent Potential Modifications

#### 2.1 Thermal Smearing of Potentials

We'll implement a method to apply thermal smearing to the calculated potentials based on the mean-square displacement:

```rust
fn apply_thermal_smearing(&self, result: &mut MuffinTinPotentialResult) -> Result<()> {
    // Skip if temperature is too low
    if self.temperature < 50.0 {
        return Ok(());
    }
    
    // Create thermal model based on parameters
    let thermal_model = self.create_thermal_model()?;
    
    // Calculate mean-square displacement at current temperature
    let msd = thermal_model.mean_square_displacement(self.temperature);
    
    // Apply Gaussian smearing to potentials
    for pot_idx in 0..self.structure.potential_type_count() {
        // Get potential values
        let potential = result.potentials[pot_idx].clone(); // Clone to avoid borrowing issues
        let grid = &result.radial_grid;
        
        // Apply Gaussian smearing with width based on mean-square displacement
        // The smearing width is related to the √MSD
        let sigma = msd.sqrt();
        
        // Skip if negligible smearing
        if sigma < 0.01 {
            continue;
        }
        
        // Create convolution kernel
        let kernel_size = (3.0 * sigma / (grid[1] - grid[0])) as usize + 1;
        let mut kernel = vec![0.0; 2 * kernel_size + 1];
        
        for i in 0..kernel.len() {
            let r = (i as isize - kernel_size as isize) as f64 * (grid[1] - grid[0]);
            kernel[i] = (-0.5 * (r / sigma).powi(2)).exp();
        }
        
        // Normalize kernel
        let kernel_sum: f64 = kernel.iter().sum();
        for k in kernel.iter_mut() {
            *k /= kernel_sum;
        }
        
        // Apply convolution
        let mut smoothed_potential = vec![0.0; grid.len()];
        
        for i in 0..grid.len() {
            let mut sum = 0.0;
            let mut weight_sum = 0.0;
            
            for j in 0..kernel.len() {
                let idx = i as isize + (j as isize - kernel_size as isize);
                
                if idx >= 0 && idx < grid.len() as isize {
                    sum += potential[idx as usize] * kernel[j];
                    weight_sum += kernel[j];
                }
            }
            
            if weight_sum > 0.0 {
                smoothed_potential[i] = sum / weight_sum;
            } else {
                smoothed_potential[i] = potential[i];
            }
        }
        
        // Update potential with smoothed version
        result.potentials[pot_idx] = smoothed_potential;
    }
    
    Ok(())
}
```

#### 2.2 Temperature-Dependent Exchange-Correlation

Extend the exchange-correlation calculation to include temperature effects:

```rust
fn calculate_xc_potential(&self, result: &mut MuffinTinPotentialResult) -> Result<()> {
    // Existing implementation...
    
    // For temperature-dependent functionals like Hedin-Lundqvist,
    // use the current temperature rather than assuming zero
    let temperature = self.temperature;
    
    // Process potential types
    for pot_idx in 0..self.structure.potential_type_count() {
        // Existing implementation...
        
        // When calculating XC potential, pass temperature for temperature-dependent functionals
        let xc_potential = match self.xc_type {
            // Existing cases...
            
            ExchangeCorrelationType::HedinLundqvist => {
                // Use temperature in calculation
                let (real_part, _imag_part) = calculate_hedin_lundqvist_potential(rho, 0.0, Some(temperature));
                real_part
            }
            
            // Other cases...
        };
        
        // Existing implementation...
    }
    
    Ok(())
}
```

### 3. Thermal Expansion Effects on Potentials

Implement methods to account for thermal expansion effects on the potential:

```rust
fn apply_thermal_expansion(&self, result: &mut MuffinTinPotentialResult) -> Result<()> {
    // Skip if temperature is too low
    if self.temperature < 50.0 {
        return Ok(());
    }
    
    // Reference temperature (typically room temperature)
    let ref_temp = 300.0;
    
    // Skip if no temperature difference
    if (self.temperature - ref_temp).abs() < 1e-6 {
        return Ok(());
    }
    
    // Get thermal expansion coefficient (typical value ~1e-5 K⁻¹)
    let alpha = 1.0e-5; // Could be material-specific in the future
    
    // Calculate scale factor due to thermal expansion
    let expansion_factor = 1.0 + alpha * (self.temperature - ref_temp);
    
    // Scale the radial grid to account for expansion
    let original_grid = result.radial_grid.clone();
    
    for i in 0..result.radial_grid.len() {
        result.radial_grid[i] = original_grid[i] * expansion_factor;
    }
    
    // Now we need to transform the potentials to match the new grid
    for pot_idx in 0..self.structure.potential_type_count() {
        let original_potential = result.potentials[pot_idx].clone();
        
        // Interpolate potential to new grid
        for i in 0..result.radial_grid.len() {
            if i == 0 {
                // Keep origin point unchanged
                continue;
            }
            
            // Find corresponding point in original grid
            let r = result.radial_grid[i] / expansion_factor;
            
            // Find bracketing points in original grid
            let mut idx_low = 0;
            let mut idx_high = original_grid.len() - 1;
            
            for j in 0..original_grid.len() {
                if original_grid[j] >= r {
                    idx_high = j;
                    if j > 0 {
                        idx_low = j - 1;
                    }
                    break;
                }
            }
            
            // Linear interpolation
            let r_low = original_grid[idx_low];
            let r_high = original_grid[idx_high];
            
            if r_high > r_low {
                let t = (r - r_low) / (r_high - r_low);
                let v_low = original_potential[idx_low];
                let v_high = original_potential[idx_high];
                
                result.potentials[pot_idx][i] = v_low + t * (v_high - v_low);
            }
        }
    }
    
    Ok(())
}
```

### 4. Integrate Changes into Main Calculation Flow

Modify the main calculation method to include these temperature effects:

```rust
pub fn calculate(&self) -> Result<MuffinTinPotentialResult> {
    // Create the radial grid for each potential type
    let mut result = self.create_grids()?;
    
    // Calculate atomic densities (non-overlapping)
    self.calculate_atomic_densities(&mut result)?;
    
    // Calculate Coulomb potential
    self.calculate_coulomb_potential(&mut result)?;
    
    // Add exchange-correlation potential
    self.calculate_xc_potential(&mut result)?;
    
    // Apply temperature effects if temperature is set
    if self.temperature > 0.0 {
        // Apply thermal expansion effects (shifts grid and interpolates potentials)
        self.apply_thermal_expansion(&mut result)?;
        
        // Apply thermal smearing to account for atomic vibrations
        self.apply_thermal_smearing(&mut result)?;
    }
    
    // Calculate energy levels by solving Schrödinger equation
    self.calculate_energy_levels(&mut result)?;
    
    // Calculate Fermi energy
    self.calculate_fermi_energy(&mut result)?;
    
    Ok(result)
}
```

### 5. Extend API to Accept Thermal Parameters from XAS Code

Modify the API for potential calculation in the XAS modules to allow passing thermal parameters:

```rust
pub fn calculate_potential_with_thermal_effects(
    structure: &AtomicStructure,
    thermal_params: Option<&ThermalParameters>,
) -> Result<MuffinTinPotentialResult> {
    // Create muffin-tin potential calculator
    let mut calculator = MuffinTinPotential::new(structure)?;
    
    // Set thermal parameters if provided
    if let Some(params) = thermal_params {
        calculator.set_temperature(params.temperature);
        
        if params.model_type != "none" {
            let param_value = match params.model_type.as_str() {
                "debye" | "correlated_debye" => params.debye_temperature,
                "einstein" => params.einstein_frequency.unwrap_or(25.0),
                _ => params.debye_temperature, // Default to Debye
            };
            
            calculator.set_thermal_model(&params.model_type, param_value)?;
        }
    }
    
    // Calculate the potential with the specified parameters
    calculator.calculate()
}
```

## Testing Plan

1. **Unit Tests**: Create tests that verify the thermal effects on potentials:
   - Test thermal smearing of potentials at different temperatures
   - Verify thermal expansion effects on potentials
   - Test behavior at extreme temperatures

2. **Integration Tests**: Test the end-to-end effects on XAS calculations:
   - Compare XAS calculations with and without temperature-dependent potentials
   - Verify expected spectral changes with temperature variations
   - Test behavior with different thermal models

3. **Validation Tests**: Compare results with FEFF10 reference calculations:
   - Generate reference temperature-dependent XANES spectra with FEFF10
   - Verify that the feff-rs implementation reproduces these results
   - Document any systematic differences

## Implementation Schedule

1. Extend `MuffinTinPotential` with temperature parameters
2. Implement thermal smearing of potentials
3. Implement temperature-dependent exchange-correlation effects
4. Implement thermal expansion effects on potentials
5. Integrate changes into calculation flow
6. Add unit tests for each component
7. Add integration tests
8. Perform validation against FEFF10 reference
9. Update documentation

## Expected Outcomes

1. More physically accurate modeling of temperature effects in XAS calculations
2. Better agreement with experimental data measured at various temperatures
3. Ability to model temperature-dependent spectral features, particularly near the edge
4. Improved handling of anharmonic effects and non-cubic materials