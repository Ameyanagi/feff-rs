/*
MIT License with FEFF10 Attribution

Copyright (c) 2025 Ameyanagi

Based on or developed using Distribution: FEFF10.0
Copyright (c) 2020 FEFF Project, University of Washington and SLAC National Accelerator Laboratory.
All rights reserved.
*/

//! Thermal parameters for XAS calculations
//!
//! This module defines structures for thermal parameters used in
//! XANES and EXAFS calculations to account for temperature effects.
//! It includes the following features:
//!
//! - Core thermal parameters structure for controlling temperature effects
//! - Support for pair-specific thermal parameters to model different bond types
//!   with appropriate vibrational characteristics
//! - Anisotropic thermal modeling for non-cubic materials with directional
//!   thermal vibration parameters
//! - Crystal system-specific factory functions for common materials
//!
//! # Effects of thermal vibrations on XAS
//!
//! Temperature affects XAS spectra in several ways:
//!
//! 1. **EXAFS amplitude reduction**: Thermal disorder decreases the amplitude
//!    of EXAFS oscillations through the Debye-Waller factor exp(-2k²σ²), where
//!    σ² is the mean-square relative displacement (MSRD) of the absorber-scatterer
//!    pair. Higher temperatures lead to larger σ² values and greater damping.
//!
//! 2. **XANES peak broadening**: Thermal vibrations cause broadening of XANES
//!    features, decreasing the sharpness of spectral peaks and reducing overall
//!    amplitude. This effect is particularly important for low-temperature experiments
//!    and theoretical calculations.
//!
//! 3. **Anisotropic effects**: In non-cubic materials, thermal vibrations have
//!    different amplitudes along different crystallographic directions. This causes
//!    direction-dependent thermal dampening, which is crucial for accurate modeling
//!    of materials like layered compounds or chain structures.
//!
//! # Using thermal parameters
//!
//! Different systems require different thermal models:
//!
//! - For simple **cubic materials** (metals, rock salt structures), the standard Debye
//!   model is usually sufficient.
//! - For materials with **stiff bonds** (diamond, covalent compounds), the Einstein
//!   model often works better due to its emphasis on localized vibrations.
//! - For accurate **EXAFS analysis**, the Correlated Debye model provides better
//!   results as it accounts for atomic motion correlation along scattering paths.
//! - For **non-cubic materials**, anisotropic models are essential to capture the
//!   directional nature of thermal vibrations.
//!
//! # Material-specific considerations
//!
//! Different material classes have characteristic thermal behaviors:
//!
//! - **Metals** typically follow the Debye model well, with Debye temperatures
//!   ranging from ~225K (Pb) to ~450K (Fe).
//! - **Oxides** often show anisotropic thermal behavior and may require
//!   correlated models for accurate EXAFS analysis.
//! - **Molecular crystals** are usually better described by Einstein models
//!   due to their localized vibrational modes.
//! - **Layered materials** (like graphite, MoS₂) have strong anisotropy with
//!   much larger vibrations perpendicular to the layers.
//! - **Chain compounds** show anisotropy with enhanced vibrations perpendicular
//!   to the chain direction.

use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// Pair key for identifying specific atom pairs
///
/// Uses atomic numbers (Z) to identify unique element pairs
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct PairKey {
    /// Atomic number of the first element (usually the absorber)
    pub z1: u32,
    /// Atomic number of the second element (usually the scatterer)
    pub z2: u32,
}

impl PairKey {
    /// Create a new pair key
    ///
    /// # Arguments
    ///
    /// * `z1` - Atomic number of the first element
    /// * `z2` - Atomic number of the second element
    ///
    /// # Returns
    ///
    /// A new pair key with elements ordered by atomic number
    pub fn new(z1: u32, z2: u32) -> Self {
        // Store in canonical order (lowest Z first)
        if z1 <= z2 {
            Self { z1, z2 }
        } else {
            Self { z1: z2, z2: z1 }
        }
    }
}

/// Bond-specific thermal parameters for a specific pair of elements
#[derive(Debug, Clone, Serialize, Deserialize)]
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

impl Default for PairThermalParameters {
    fn default() -> Self {
        Self {
            model_type: "debye".to_string(),
            debye_temperature: 300.0,
            einstein_frequency: None,
            description: None,
        }
    }
}

/// Thermal parameters for XAS calculations
#[derive(Debug, Clone, Serialize, Deserialize)]
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
    /// Values represent relative thermal vibration amplitudes along each crystallographic axis
    pub displacement_factors: Option<[f64; 3]>,
    /// Path direction in crystal coordinates [d_x, d_y, d_z] (for anisotropic model)
    /// Used to project displacement factors onto specific scattering paths
    pub path_direction: Option<[f64; 3]>,
}

impl Default for ThermalParameters {
    fn default() -> Self {
        Self {
            temperature: 300.0, // Room temperature
            model_type: "debye".to_string(),
            debye_temperature: 300.0, // Default Debye temperature
            einstein_frequency: None,
            pair_parameters: None,
            displacement_factors: None,
            path_direction: None,
        }
    }
}

impl ThermalParameters {
    /// Create new thermal parameters for a Debye model
    ///
    /// # Arguments
    ///
    /// * `temperature` - Temperature in Kelvin
    /// * `debye_temperature` - Debye temperature in Kelvin
    ///
    /// # Returns
    ///
    /// New thermal parameters
    pub fn new_debye(temperature: f64, debye_temperature: f64) -> Self {
        Self {
            temperature,
            model_type: "debye".to_string(),
            debye_temperature,
            einstein_frequency: None,
            pair_parameters: None,
            displacement_factors: None,
            path_direction: None,
        }
    }

    /// Create new thermal parameters for an Einstein model
    ///
    /// # Arguments
    ///
    /// * `temperature` - Temperature in Kelvin
    /// * `einstein_frequency` - Einstein frequency in meV
    ///
    /// # Returns
    ///
    /// New thermal parameters
    pub fn new_einstein(temperature: f64, einstein_frequency: f64) -> Self {
        Self {
            temperature,
            model_type: "einstein".to_string(),
            debye_temperature: 0.0, // Not used for Einstein model
            einstein_frequency: Some(einstein_frequency),
            pair_parameters: None,
            displacement_factors: None,
            path_direction: None,
        }
    }

    /// Create new thermal parameters for a Correlated Debye model
    ///
    /// # Arguments
    ///
    /// * `temperature` - Temperature in Kelvin
    /// * `debye_temperature` - Debye temperature in Kelvin
    ///
    /// # Returns
    ///
    /// New thermal parameters with a correlated Debye model
    pub fn new_correlated_debye(temperature: f64, debye_temperature: f64) -> Self {
        Self {
            temperature,
            model_type: "correlated_debye".to_string(),
            debye_temperature,
            einstein_frequency: None,
            pair_parameters: None,
            displacement_factors: None,
            path_direction: None,
        }
    }

    /// Create new thermal parameters for an anisotropic Debye model
    ///
    /// This creates a thermal model that accounts for directional thermal vibrations,
    /// which is essential for correctly modeling non-cubic materials like tetragonal,
    /// hexagonal, or orthorhombic crystals where thermal vibrations have different
    /// amplitudes along different crystallographic axes.
    ///
    /// The anisotropic Debye model combines:
    /// - The temperature dependence of the Debye model for the overall vibration magnitude
    /// - Direction-dependent scaling to account for crystallographic anisotropy
    ///
    /// # Physical significance
    ///
    /// - In **layered materials** (graphite, transition metal dichalcogenides),
    ///   out-of-plane vibrations are typically much larger than in-plane vibrations,
    ///   which can be modeled with displacement factors like [1.0, 1.0, 2.0].
    ///   
    /// - In **chain compounds** (polymers, 1D conductors), vibrations perpendicular to
    ///   the chain direction are typically larger than along the chain, which can be
    ///   modeled with displacement factors like [2.0, 2.0, 1.0] if the chain runs along z.
    ///   
    /// - **Path direction** determines how these anisotropic vibrations project onto
    ///   a specific scattering path. For example, a path that runs parallel to the
    ///   z-axis will be primarily affected by z-direction displacement factors.
    ///
    /// # Arguments
    ///
    /// * `temperature` - Temperature in Kelvin
    /// * `debye_temperature` - Debye temperature in Kelvin
    /// * `displacement_factors` - Relative thermal vibration amplitudes [u_x, u_y, u_z]
    /// * `path_direction` - Optional scattering path direction [d_x, d_y, d_z]
    ///
    /// # Returns
    ///
    /// New thermal parameters with an anisotropic Debye model
    pub fn new_anisotropic_debye(
        temperature: f64,
        debye_temperature: f64,
        displacement_factors: [f64; 3],
        path_direction: Option<[f64; 3]>,
    ) -> Self {
        // Use default path direction if none provided
        let direction = path_direction.unwrap_or([0.0, 0.0, 1.0]);

        Self {
            temperature,
            model_type: "anisotropic".to_string(),
            debye_temperature,
            einstein_frequency: None,
            pair_parameters: None,
            displacement_factors: Some(displacement_factors),
            path_direction: Some(direction),
        }
    }

    /// Create new thermal parameters for an anisotropic Einstein model
    ///
    /// This creates a thermal model that accounts for directional thermal vibrations
    /// using the Einstein model as a base.
    ///
    /// # Arguments
    ///
    /// * `temperature` - Temperature in Kelvin
    /// * `einstein_frequency` - Einstein frequency in meV
    /// * `displacement_factors` - Relative thermal vibration amplitudes [u_x, u_y, u_z]
    /// * `path_direction` - Optional scattering path direction [d_x, d_y, d_z]
    ///
    /// # Returns
    ///
    /// New thermal parameters with an anisotropic Einstein model
    pub fn new_anisotropic_einstein(
        temperature: f64,
        einstein_frequency: f64,
        displacement_factors: [f64; 3],
        path_direction: Option<[f64; 3]>,
    ) -> Self {
        // Use default path direction if none provided
        let direction = path_direction.unwrap_or([0.0, 0.0, 1.0]);

        Self {
            temperature,
            model_type: "anisotropic".to_string(),
            debye_temperature: 0.0, // Not used for Einstein-based model
            einstein_frequency: Some(einstein_frequency),
            pair_parameters: None,
            displacement_factors: Some(displacement_factors),
            path_direction: Some(direction),
        }
    }

    /// Create new thermal parameters for an anisotropic correlated Debye model
    ///
    /// This creates a thermal model that accounts for both directional thermal vibrations
    /// and correlation effects between atoms.
    ///
    /// # Arguments
    ///
    /// * `temperature` - Temperature in Kelvin
    /// * `debye_temperature` - Debye temperature in Kelvin
    /// * `displacement_factors` - Relative thermal vibration amplitudes [u_x, u_y, u_z]
    /// * `path_direction` - Optional scattering path direction [d_x, d_y, d_z]
    ///
    /// # Returns
    ///
    /// New thermal parameters with an anisotropic correlated Debye model
    pub fn new_anisotropic_correlated_debye(
        temperature: f64,
        debye_temperature: f64,
        displacement_factors: [f64; 3],
        path_direction: Option<[f64; 3]>,
    ) -> Self {
        // Use default path direction if none provided
        let direction = path_direction.unwrap_or([0.0, 0.0, 1.0]);

        Self {
            temperature,
            model_type: "anisotropic_correlated_debye".to_string(),
            debye_temperature,
            einstein_frequency: None,
            pair_parameters: None,
            displacement_factors: Some(displacement_factors),
            path_direction: Some(direction),
        }
    }

    /// Add pair-specific thermal parameters for a specific element pair
    ///
    /// # Arguments
    ///
    /// * `z1` - Atomic number of the first element
    /// * `z2` - Atomic number of the second element
    /// * `pair_params` - Thermal parameters for this specific pair
    ///
    /// # Returns
    ///
    /// Self for method chaining
    pub fn with_pair_parameters(
        mut self,
        z1: i32,
        z2: i32,
        pair_params: PairThermalParameters,
    ) -> Self {
        // Ensure valid atomic numbers
        if z1 <= 0 || z2 <= 0 {
            return self;
        }

        let pair_key = PairKey::new(z1 as u32, z2 as u32);

        // Initialize the pair parameters map if it doesn't exist
        if self.pair_parameters.is_none() {
            self.pair_parameters = Some(HashMap::new());
        }

        // Add or update the pair parameters
        if let Some(params) = &mut self.pair_parameters {
            params.insert(pair_key, pair_params);
        }

        self
    }

    /// Get thermal parameters for a specific element pair
    ///
    /// # Arguments
    ///
    /// * `z1` - Atomic number of the first element
    /// * `z2` - Atomic number of the second element
    ///
    /// # Returns
    ///
    /// Pair-specific thermal parameters if they exist, or None
    pub fn get_pair_parameters(&self, z1: i32, z2: i32) -> Option<&PairThermalParameters> {
        // Convert i32 to u32, ensuring atomic numbers are positive
        if z1 <= 0 || z2 <= 0 {
            return None;
        }

        let pair_key = PairKey::new(z1 as u32, z2 as u32);

        self.pair_parameters
            .as_ref()
            .and_then(|params| params.get(&pair_key))
    }

    /// Create a thermal model from these parameters
    ///
    /// # Arguments
    ///
    /// * `reduced_mass` - Reduced mass in atomic mass units (amu)
    /// * `path_distance` - Optional path distance in Å (required for correlated models)
    ///
    /// # Returns
    ///
    /// A boxed thermal model configured according to these parameters
    pub fn create_model(
        &self,
        reduced_mass: f64,
        path_distance: Option<f64>,
    ) -> Box<dyn crate::utils::thermal::ThermalModel> {
        use crate::utils::thermal::{
            create_anisotropic_thermal_model, create_thermal_model, create_thermal_model_with_path,
        };

        match self.model_type.as_str() {
            "debye" => {
                create_thermal_model("debye", Some(self.debye_temperature), None, reduced_mass)
            }
            "einstein" => {
                create_thermal_model("einstein", None, self.einstein_frequency, reduced_mass)
            }
            "correlated_debye" => {
                let distance = path_distance.unwrap_or(2.5); // Default to typical first shell
                create_thermal_model_with_path(
                    "correlated_debye",
                    Some(self.debye_temperature),
                    None,
                    reduced_mass,
                    distance,
                )
            }
            "anisotropic" | "anisotropic_correlated_debye" => {
                // For anisotropic models, we need the displacement factors and path direction
                let factors = self.displacement_factors.unwrap_or([1.0, 1.0, 1.0]);
                let direction = self.path_direction.unwrap_or([0.0, 0.0, 1.0]);

                // Choose base model type
                let base_model = if self.model_type == "anisotropic_correlated_debye" {
                    "correlated_debye"
                } else if self.einstein_frequency.is_some() {
                    "einstein"
                } else {
                    "debye"
                };

                // Create anisotropic model
                create_anisotropic_thermal_model(
                    base_model,
                    Some(self.debye_temperature),
                    self.einstein_frequency,
                    reduced_mass,
                    path_distance,
                    factors,
                    direction,
                )
            }
            _ => {
                // Default to Debye if model type not recognized
                create_thermal_model("debye", Some(self.debye_temperature), None, reduced_mass)
            }
        }
    }
}

impl PairThermalParameters {
    /// Create new pair thermal parameters for a Debye model
    ///
    /// # Arguments
    ///
    /// * `debye_temperature` - Debye temperature in Kelvin
    /// * `description` - Optional description of this bond type
    ///
    /// # Returns
    ///
    /// New pair thermal parameters
    pub fn new_debye(debye_temperature: f64, description: Option<String>) -> Self {
        Self {
            model_type: "debye".to_string(),
            debye_temperature,
            einstein_frequency: None,
            description,
        }
    }

    /// Create new pair thermal parameters for an Einstein model
    ///
    /// # Arguments
    ///
    /// * `einstein_frequency` - Einstein frequency in meV
    /// * `description` - Optional description of this bond type
    ///
    /// # Returns
    ///
    /// New pair thermal parameters
    pub fn new_einstein(einstein_frequency: f64, description: Option<String>) -> Self {
        Self {
            model_type: "einstein".to_string(),
            debye_temperature: 0.0, // Not used for Einstein model
            einstein_frequency: Some(einstein_frequency),
            description,
        }
    }
}

/// Create thermal parameters with common element pairs for typical materials
///
/// # Arguments
///
/// * `temperature` - Temperature in Kelvin for the thermal parameters
///
/// # Returns
///
/// ThermalParameters populated with typical values for common bonds
pub fn create_standard_thermal_parameters(temperature: f64) -> ThermalParameters {
    // Start with base parameters (Debye model at specified temperature)
    let mut params = ThermalParameters::new_debye(temperature, 300.0);

    // Add pair-specific parameters for common element pairs
    // These values are based on literature and experiment

    // Metal-metal bonds (use Debye model)
    // Cu-Cu bond
    params = params.with_pair_parameters(
        29,
        29,
        PairThermalParameters::new_debye(315.0, Some("Cu-Cu bond (metal)".to_string())),
    );

    // Fe-Fe bond
    params = params.with_pair_parameters(
        26,
        26,
        PairThermalParameters::new_debye(470.0, Some("Fe-Fe bond (metal)".to_string())),
    );

    // Ni-Ni bond
    params = params.with_pair_parameters(
        28,
        28,
        PairThermalParameters::new_debye(375.0, Some("Ni-Ni bond (metal)".to_string())),
    );

    // Metal-oxide bonds (use Einstein model due to localized vibrations)
    // Cu-O bond
    params = params.with_pair_parameters(
        29,
        8,
        PairThermalParameters::new_einstein(
            70.0, // ~70 meV typical for Cu-O vibrations
            Some("Cu-O bond (metal-oxide)".to_string()),
        ),
    );

    // Fe-O bond
    params = params.with_pair_parameters(
        26,
        8,
        PairThermalParameters::new_einstein(
            65.0, // ~65 meV typical for Fe-O vibrations
            Some("Fe-O bond (metal-oxide)".to_string()),
        ),
    );

    // Ni-O bond
    params = params.with_pair_parameters(
        28,
        8,
        PairThermalParameters::new_einstein(
            68.0, // ~68 meV typical for Ni-O vibrations
            Some("Ni-O bond (metal-oxide)".to_string()),
        ),
    );

    // Covalent bonds (use Einstein model due to stiffness)
    // C-C bond
    params = params.with_pair_parameters(
        6,
        6,
        PairThermalParameters::new_einstein(
            160.0, // Stiff covalent bond
            Some("C-C covalent bond".to_string()),
        ),
    );

    // Si-O bond
    params = params.with_pair_parameters(
        14,
        8,
        PairThermalParameters::new_einstein(
            130.0, // Typical for Si-O in silicates
            Some("Si-O tetrahedral bond".to_string()),
        ),
    );

    params
}

/// Create anisotropic thermal parameters for different crystal systems
///
/// This function creates appropriate anisotropic thermal parameters for
/// various crystal structures based on their symmetry properties. It automatically
/// sets the displacement factors according to typical patterns for each crystal system.
///
/// # Crystal systems and their thermal characteristics
///
/// - **Cubic** (ex: NaCl, Cu, diamond): Equal thermal vibrations in all directions.
///   Displacement factors: [1.0, 1.0, 1.0]
///   
/// - **Tetragonal** (ex: TiO₂-rutile, CaSO₄): Equal vibrations in the a and b
///   directions (x and y), but different along c (z). Often, c-axis vibrations
///   are larger.
///   Displacement factors: [1.0, 1.0, ratio]
///   
/// - **Hexagonal** (ex: Zn, graphite, wurtzite): Similar to tetragonal, with
///   equal a and b, but different c. Typically has more pronounced anisotropy
///   than tetragonal.
///   Displacement factors: [1.0, 1.0, ratio*1.1]
///   
/// - **Orthorhombic** (ex: MgSO₄, aragonite): Three different axes with different
///   thermal vibration amplitudes.
///   Displacement factors: [1.0, ratio*0.8, ratio]
///   
/// - **Monoclinic** (ex: gypsum, clinopyroxenes): Complex anisotropy with the
///   unique axis (usually b) having distinct vibration amplitudes.
///   Displacement factors: [1.0, ratio, ratio*0.9]
///   
/// - **Triclinic** (ex: albite, microcline): Most complex anisotropy with all
///   three axes having different vibration amplitudes.
///   Displacement factors: [1.0, ratio*0.8, ratio*1.2]
///   
/// - **Layered materials** (ex: graphite, MoS₂): Very high anisotropy with
///   much larger vibrations perpendicular to the layers.
///   Displacement factors: [1.0, 1.0, ratio*2.0]
///   
/// - **Chain compounds** (ex: polymers): Anisotropy with larger vibrations
///   perpendicular to the chain direction.
///   Displacement factors: [ratio, ratio, 1.0] (assuming chain along z)
///
/// # Arguments
///
/// * `temperature` - Temperature in Kelvin
/// * `crystal_system` - Crystal system type ("cubic", "tetragonal", "hexagonal", "orthorhombic", etc.)
/// * `debye_temperature` - Debye temperature in Kelvin
/// * `anisotropy_ratio` - Ratio of thermal vibration amplitudes along different axes
///   (default is 1.5 if not specified; higher values create more pronounced anisotropy)
///
/// # Examples
///
/// ```no_run
/// use feff_rs::xas::thermal::create_anisotropic_thermal_parameters;
/// // Create parameters for a tetragonal material at room temperature
/// let _tetragonal_params = create_anisotropic_thermal_parameters(
///     300.0,             // Room temperature
///     "tetragonal",      // Crystal system
///     350.0,             // Debye temperature
///     Some(1.7)          // Higher anisotropy ratio
/// );
///
/// // Create parameters for a layered material at low temperature
/// let _graphite_params = create_anisotropic_thermal_parameters(
///     77.0,              // Liquid nitrogen temperature
///     "layered",         // Layered material
///     420.0,             // Debye temperature for graphite
///     Some(2.5)          // Very high anisotropy ratio
/// );
/// ```
///
/// # Returns
///
/// ThermalParameters configured with appropriate anisotropic settings for the
/// specified crystal system.
pub fn create_anisotropic_thermal_parameters(
    temperature: f64,
    crystal_system: &str,
    debye_temperature: f64,
    anisotropy_ratio: Option<f64>,
) -> ThermalParameters {
    // Default anisotropy ratio (higher means more anisotropy)
    let ratio = anisotropy_ratio.unwrap_or(1.5);

    // Set displacement factors based on crystal system
    let displacement_factors = match crystal_system.to_lowercase().as_str() {
        "cubic" => {
            // Cubic systems have equal thermal vibrations in all directions
            [1.0, 1.0, 1.0]
        }
        "tetragonal" => {
            // Tetragonal systems have equal x,y vibrations but different z
            // Typically, c-axis (z) vibrations are larger in layered materials
            [1.0, 1.0, ratio]
        }
        "hexagonal" => {
            // Hexagonal systems have equal x,y vibrations but different z
            // Similar to tetragonal, but often with more pronounced anisotropy
            [1.0, 1.0, ratio * 1.1]
        }
        "trigonal" | "rhombohedral" => {
            // Trigonal systems have complex anisotropy
            [1.0, 1.0, ratio]
        }
        "orthorhombic" => {
            // Orthorhombic systems have different vibrations along all three axes
            [1.0, ratio * 0.8, ratio]
        }
        "monoclinic" => {
            // Monoclinic systems have more complex anisotropy
            // Typically, the vibrations along the unique axis (usually b) differ
            [1.0, ratio, ratio * 0.9]
        }
        "triclinic" => {
            // Triclinic systems have the most complex anisotropy
            // All three axes have different vibration amplitudes
            [1.0, ratio * 0.8, ratio * 1.2]
        }
        "layered" => {
            // Special case for layered materials with very high anisotropy
            // Examples: graphite, transition metal dichalcogenides
            [1.0, 1.0, ratio * 2.0]
        }
        "chain" => {
            // Special case for chain-like materials
            // Examples: polymers, one-dimensional conductors
            // Assuming chains run along z-axis
            [ratio, ratio, 1.0]
        }
        _ => {
            // Default to isotropic for unknown crystal systems
            [1.0, 1.0, 1.0]
        }
    };

    // Create anisotropic thermal parameters
    // Use correlated model for better physical accuracy
    ThermalParameters::new_anisotropic_correlated_debye(
        temperature,
        debye_temperature,
        displacement_factors,
        None, // Default path direction (will be set when used)
    )
}

/// Create material-specific thermal parameters
///
/// This function provides optimized thermal parameters for specific materials,
/// including appropriate Debye temperatures, Einstein frequencies, and anisotropic
/// factors based on the material type. It serves as a convenient way to get
/// reasonable starting values for common materials without needing to look up
/// parameters manually.
///
/// # Arguments
///
/// * `material` - Material name or type (e.g., "Cu", "Fe2O3", "graphite")
/// * `temperature` - Temperature in Kelvin
/// * `custom_debye_temp` - Optional override for the material's Debye temperature
///
/// # Returns
///
/// ThermalParameters with material-specific values
pub fn create_material_thermal_parameters(
    material: &str,
    temperature: f64,
    custom_debye_temp: Option<f64>,
) -> ThermalParameters {
    // Convert material name to lowercase for case-insensitive matching
    let material_lower = material.to_lowercase();

    // Standard metals
    if material_lower == "cu" || material_lower == "copper" {
        let debye_temp = custom_debye_temp.unwrap_or(315.0);
        let mut params = ThermalParameters::new_debye(temperature, debye_temp);
        params = params.with_pair_parameters(
            29,
            29,
            PairThermalParameters::new_debye(debye_temp, Some("Cu-Cu bond".to_string())),
        );
        return params;
    }

    if material_lower == "fe" || material_lower == "iron" {
        let debye_temp = custom_debye_temp.unwrap_or(470.0);
        let mut params = ThermalParameters::new_debye(temperature, debye_temp);
        params = params.with_pair_parameters(
            26,
            26,
            PairThermalParameters::new_debye(debye_temp, Some("Fe-Fe bond".to_string())),
        );
        return params;
    }

    if material_lower == "ni" || material_lower == "nickel" {
        let debye_temp = custom_debye_temp.unwrap_or(375.0);
        let mut params = ThermalParameters::new_debye(temperature, debye_temp);
        params = params.with_pair_parameters(
            28,
            28,
            PairThermalParameters::new_debye(debye_temp, Some("Ni-Ni bond".to_string())),
        );
        return params;
    }

    if material_lower == "au" || material_lower == "gold" {
        let debye_temp = custom_debye_temp.unwrap_or(170.0);
        let mut params = ThermalParameters::new_debye(temperature, debye_temp);
        params = params.with_pair_parameters(
            79,
            79,
            PairThermalParameters::new_debye(debye_temp, Some("Au-Au bond".to_string())),
        );
        return params;
    }

    if material_lower == "pt" || material_lower == "platinum" {
        let debye_temp = custom_debye_temp.unwrap_or(230.0);
        let mut params = ThermalParameters::new_debye(temperature, debye_temp);
        params = params.with_pair_parameters(
            78,
            78,
            PairThermalParameters::new_debye(debye_temp, Some("Pt-Pt bond".to_string())),
        );
        return params;
    }

    // Silicon and diamond materials
    if material_lower == "si" || material_lower == "silicon" {
        let debye_temp = custom_debye_temp.unwrap_or(645.0);
        let mut params = ThermalParameters::new_debye(temperature, debye_temp);
        params = params.with_pair_parameters(
            14,
            14,
            PairThermalParameters::new_debye(debye_temp, Some("Si-Si covalent bond".to_string())),
        );
        return params;
    }

    if material_lower == "diamond" || material_lower == "c" {
        let debye_temp = custom_debye_temp.unwrap_or(2230.0);
        let mut params = ThermalParameters::new_debye(temperature, debye_temp);
        // For diamond, Einstein model better describes the stiff C-C bonds
        params = params.with_pair_parameters(
            6,
            6,
            PairThermalParameters::new_einstein(160.0, Some("C-C diamond bond".to_string())),
        );
        return params;
    }

    // Oxide materials (anisotropic)
    if material_lower == "tio2" || material_lower.contains("rutile") {
        let debye_temp = custom_debye_temp.unwrap_or(450.0);
        // TiO2 has tetragonal structure
        let mut params = create_anisotropic_thermal_parameters(
            temperature,
            "tetragonal",
            debye_temp,
            Some(1.4), // Moderate anisotropy
        );
        // Add pair parameters for Ti-O bonds
        params = params.with_pair_parameters(
            22,
            8,
            PairThermalParameters::new_einstein(75.0, Some("Ti-O bond".to_string())),
        );
        return params;
    }

    if material_lower == "fe2o3" || material_lower.contains("hematite") {
        let debye_temp = custom_debye_temp.unwrap_or(500.0);
        // Fe2O3 has rhombohedral/trigonal structure
        let mut params = create_anisotropic_thermal_parameters(
            temperature,
            "trigonal",
            debye_temp,
            Some(1.3), // Moderate anisotropy
        );
        // Add pair parameters for Fe-O bonds
        params = params.with_pair_parameters(
            26,
            8,
            PairThermalParameters::new_einstein(65.0, Some("Fe-O bond".to_string())),
        );
        return params;
    }

    // Layered materials (highly anisotropic)
    if material_lower == "graphite" || material_lower == "graphene" {
        let debye_temp = custom_debye_temp.unwrap_or(420.0);
        // Graphite has extremely high anisotropy
        let mut params = create_anisotropic_thermal_parameters(
            temperature,
            "layered",
            debye_temp,
            Some(3.0), // Very high anisotropy
        );
        // Add pair parameters for C-C bonds (in-plane vs interlayer)
        params = params.with_pair_parameters(
            6,
            6,
            PairThermalParameters::new_einstein(180.0, Some("C-C in-plane bond".to_string())),
        );
        return params;
    }

    if material_lower == "mos2" || material_lower.contains("disulfide") {
        let debye_temp = custom_debye_temp.unwrap_or(380.0);
        // MoS2 has hexagonal layered structure
        let mut params = create_anisotropic_thermal_parameters(
            temperature,
            "layered",
            debye_temp,
            Some(2.5), // High anisotropy
        );
        // Add pair parameters for Mo-S bonds
        params = params.with_pair_parameters(
            42,
            16,
            PairThermalParameters::new_einstein(60.0, Some("Mo-S bond".to_string())),
        );
        return params;
    }

    // If no specific material match, use a general approach based on material class

    // Generic material classes
    if material_lower.contains("metal") {
        // Use Debye model for metals
        let debye_temp = custom_debye_temp.unwrap_or(350.0);
        return ThermalParameters::new_debye(temperature, debye_temp);
    }

    if material_lower.contains("oxide") {
        // Use anisotropic model for oxides
        let debye_temp = custom_debye_temp.unwrap_or(450.0);
        return create_anisotropic_thermal_parameters(
            temperature,
            "tetragonal", // Most oxides have tetragonal or similar structure
            debye_temp,
            Some(1.5),
        );
    }

    if material_lower.contains("layer") {
        // Use highly anisotropic model for layered materials
        let debye_temp = custom_debye_temp.unwrap_or(400.0);
        return create_anisotropic_thermal_parameters(
            temperature,
            "layered",
            debye_temp,
            Some(2.5),
        );
    }

    // Default to standard Debye model with user-provided or reasonable Debye temperature
    let debye_temp = custom_debye_temp.unwrap_or(300.0);
    ThermalParameters::new_debye(temperature, debye_temp)
}

/// Create thermal parameters optimized for EXAFS analysis
///
/// EXAFS requires special consideration for thermal effects, particularly:
/// 1. Correlation between atomic motions along scattering paths
/// 2. Anharmonic effects at high temperatures
/// 3. Anisotropic vibrations in non-cubic materials
///
/// This function provides appropriate thermal parameters for EXAFS analysis
/// of different material types.
///
/// # Arguments
///
/// * `material_type` - Material class ("metal", "oxide", "molecular", "layered", etc.)
/// * `temperature` - Temperature in Kelvin
/// * `debye_temperature` - Debye temperature in Kelvin
/// * `include_anharmonic` - Whether to include anharmonic corrections (important above ~500K)
///
/// # Returns
///
/// ThermalParameters optimized for EXAFS analysis
pub fn create_exafs_thermal_parameters(
    material_type: &str,
    temperature: f64,
    debye_temperature: f64,
    include_anharmonic: bool,
) -> ThermalParameters {
    // Lowercase for case-insensitive matching
    let material_lower = material_type.to_lowercase();

    // Adjust Debye temperature for anharmonic effects
    // At high temperatures, the effective Debye temperature decreases
    let effective_debye_temp = if include_anharmonic && temperature > 300.0 {
        // Empirical adjustment: reduce Debye temperature by up to 20% at high temperatures
        let reduction_factor = 1.0 - 0.2 * (temperature - 300.0) / 700.0;
        debye_temperature * reduction_factor.max(0.8)
    } else {
        debye_temperature
    };

    match material_lower.as_str() {
        "metal" | "metallic" => {
            // For metals, use correlated Debye model which works well for EXAFS
            ThermalParameters::new_correlated_debye(temperature, effective_debye_temp)
        }
        "oxide" | "ceramic" => {
            // For oxides, use anisotropic model
            let anisotropy = if include_anharmonic && temperature > 300.0 {
                // Increase anisotropy at high temperatures
                1.5 + 0.5 * (temperature - 300.0) / 700.0
            } else {
                1.5 // Moderate anisotropy at room temperature
            };

            create_anisotropic_thermal_parameters(
                temperature,
                "tetragonal", // Common structure for many oxides
                effective_debye_temp,
                Some(anisotropy),
            )
        }
        "layered" | "2d" => {
            // Highly anisotropic model for layered materials
            let anisotropy = if include_anharmonic && temperature > 300.0 {
                // Increase anisotropy at high temperatures
                2.5 + (temperature - 300.0) / 700.0
            } else {
                2.5 // High anisotropy at room temperature
            };

            create_anisotropic_thermal_parameters(
                temperature,
                "layered",
                effective_debye_temp,
                Some(anisotropy),
            )
        }
        "molecular" | "organic" => {
            // Einstein model works better for molecular crystals
            // Convert Debye temperature to approximate Einstein frequency
            let einstein_freq = debye_temperature * 0.8 * 0.695; // in meV
            ThermalParameters::new_einstein(temperature, einstein_freq)
        }
        "chain" | "polymer" => {
            // Anisotropic model for chain compounds
            let anisotropy = if include_anharmonic && temperature > 300.0 {
                2.0 + 0.5 * (temperature - 300.0) / 700.0
            } else {
                2.0
            };

            create_anisotropic_thermal_parameters(
                temperature,
                "chain",
                effective_debye_temp,
                Some(anisotropy),
            )
        }
        _ => {
            // Default to correlated Debye model for unknown materials
            ThermalParameters::new_correlated_debye(temperature, effective_debye_temp)
        }
    }
}

/// Create thermal parameters optimized for XANES analysis
///
/// XANES has different requirements for thermal modeling compared to EXAFS:
/// 1. Broadening effects dominate over Debye-Waller dampening
/// 2. Core-hole lifetime often outweighs thermal broadening
/// 3. Multiple-scattering pathways require accounting for 3D vibrations
///
/// This function provides thermal parameters optimized for XANES analysis.
///
/// # Arguments
///
/// * `material_type` - Material class ("metal", "oxide", "molecular", "layered", etc.)
/// * `temperature` - Temperature in Kelvin
/// * `debye_temperature` - Debye temperature in Kelvin
/// * `edge` - Absorption edge ("K", "L1", "L2", "L3", "M1", etc.)
///
/// # Returns
///
/// ThermalParameters optimized for XANES analysis
pub fn create_xanes_thermal_parameters(
    material_type: &str,
    temperature: f64,
    debye_temperature: f64,
    edge: &str,
) -> ThermalParameters {
    // Lowercase for case-insensitive matching
    let material_lower = material_type.to_lowercase();
    let edge_lower = edge.to_lowercase();

    // Different edges require different models
    let is_k_edge = edge_lower.starts_with('k');
    let _is_l_edge = edge_lower.starts_with('l');

    // High temperature switch
    let high_temperature = temperature > 400.0;

    // For high temperature XANES, we need to account for anharmonicity
    let model_suffix = if high_temperature { "_anharmonic" } else { "" };

    match material_lower.as_str() {
        "metal" | "metallic" => {
            if high_temperature {
                // For high temperature, use correlated model to account for multiple scattering
                ThermalParameters {
                    temperature,
                    model_type: format!("correlated_debye{}", model_suffix),
                    debye_temperature,
                    ..ThermalParameters::default()
                }
            } else if is_k_edge {
                // K-edge XANES is well described by Debye model at moderate temperatures
                ThermalParameters::new_debye(temperature, debye_temperature)
            } else {
                // L and M edges benefit from correlation effects
                ThermalParameters::new_correlated_debye(temperature, debye_temperature)
            }
        }
        "oxide" | "ceramic" => {
            // Oxides require anisotropic model for accurate XANES
            let anisotropy = if high_temperature { 1.8 } else { 1.5 };

            create_anisotropic_thermal_parameters(
                temperature,
                "tetragonal",
                debye_temperature,
                Some(anisotropy),
            )
        }
        "layered" | "2d" => {
            // Layered materials need strong anisotropy consideration
            let anisotropy = if high_temperature { 3.0 } else { 2.5 };

            // For high temperatures, add anharmonic tag
            let model_type = if high_temperature {
                "anisotropic_anharmonic"
            } else {
                "anisotropic"
            };

            ThermalParameters {
                temperature,
                model_type: model_type.to_string(),
                debye_temperature,
                displacement_factors: Some([1.0, 1.0, anisotropy]),
                ..ThermalParameters::default()
            }
        }
        "molecular" | "organic" => {
            // Einstein model works better for molecular materials
            let einstein_freq = debye_temperature * 0.8 * 0.695; // in meV

            // For XANES of molecular materials, we need to consider correlation
            if is_k_edge && high_temperature {
                ThermalParameters {
                    temperature,
                    model_type: format!("einstein_correlated{}", model_suffix),
                    debye_temperature: 0.0,
                    einstein_frequency: Some(einstein_freq),
                    ..ThermalParameters::default()
                }
            } else {
                ThermalParameters::new_einstein(temperature, einstein_freq)
            }
        }
        _ => {
            // Default case - use standard Debye model
            if high_temperature {
                ThermalParameters {
                    temperature,
                    model_type: format!("debye{}", model_suffix),
                    debye_temperature,
                    ..ThermalParameters::default()
                }
            } else {
                ThermalParameters::new_debye(temperature, debye_temperature)
            }
        }
    }
}
