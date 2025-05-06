/*
MIT License with FEFF10 Attribution

Copyright (c) 2025 Ameyanagi

Based on or developed using Distribution: FEFF10.0
Copyright (c) 2020 FEFF Project, University of Washington and SLAC National Accelerator Laboratory.
All rights reserved.
*/

//! Exchange-correlation functionals for muffin-tin potentials

use super::errors::{PotentialError, Result};

/// Exchange-correlation functional types
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ExchangeCorrelationType {
    /// Local Density Approximation (LDA)
    LDA,
    /// Generalized Gradient Approximation (GGA)
    GGA,
    /// Hedin-Lundqvist complex potential
    HedinLundqvist,
    /// Dirac-Hara complex potential
    DiracHara,
    /// von Barth-Hedin potential
    VonBarth,
}

impl ExchangeCorrelationType {
    /// Create a new exchange-correlation functional type from string
    pub fn from_string(name: &str) -> Result<Self> {
        match name.to_uppercase().as_str() {
            "LDA" => Ok(ExchangeCorrelationType::LDA),
            "GGA" => Ok(ExchangeCorrelationType::GGA),
            "HL" | "HEDIN-LUNDQVIST" | "HEDINLUNDQVIST" => {
                Ok(ExchangeCorrelationType::HedinLundqvist)
            }
            "DH" | "DIRAC-HARA" | "DIRACHARA" => Ok(ExchangeCorrelationType::DiracHara),
            "VBH" | "VON-BARTH" | "VONBARTH" | "VON-BARTH-HEDIN" => {
                Ok(ExchangeCorrelationType::VonBarth)
            }
            _ => Err(PotentialError::InvalidExchangeCorrelation(format!(
                "Unknown exchange-correlation functional: {}",
                name
            ))),
        }
    }

    /// Get a string representation of the exchange-correlation type
    pub fn as_str(&self) -> &'static str {
        match self {
            ExchangeCorrelationType::LDA => "LDA",
            ExchangeCorrelationType::GGA => "GGA",
            ExchangeCorrelationType::HedinLundqvist => "Hedin-Lundqvist",
            ExchangeCorrelationType::DiracHara => "Dirac-Hara",
            ExchangeCorrelationType::VonBarth => "von Barth-Hedin",
        }
    }
}

/// Calculate LDA exchange-correlation potential
///
/// Uses the Perdew-Zunger parameterization of the Ceperley-Alder data
///
/// # Arguments
/// * `density` - electron density in electrons/bohr^3
///
/// # Returns
/// * Exchange-correlation potential in Hartree
pub fn calculate_lda_potential(density: f64) -> f64 {
    // Prevent division by zero or negative density
    if density <= 1e-12 {
        return 0.0;
    }

    // Wigner-Seitz radius (r_s = (3/4π*ρ)^(1/3))
    let rs = (3.0 / (4.0 * std::f64::consts::PI * density)).powf(1.0 / 3.0);

    // Parameters for Perdew-Zunger parameterization
    if rs < 1.0 {
        // High-density formula
        let a = 0.0311;
        let b = -0.0480;
        let c = 0.0020;
        let d = -0.0116;

        // Exchange term
        let exchange = -0.458 / rs;

        // Correlation term
        let log_rs = rs.ln();
        let correlation = a * log_rs + b + c * rs * log_rs + d * rs;

        exchange + correlation
    } else {
        // Low-density formula
        let gamma = -0.1423;
        let beta1 = 1.0529;
        let beta2 = 0.3334;

        // Exchange term
        let exchange = -0.458 / rs;

        // Correlation term
        let correlation = gamma / (1.0 + beta1 * rs.sqrt() + beta2 * rs);

        exchange + correlation
    }
}

/// Calculate GGA exchange-correlation potential
///
/// Uses the Perdew-Burke-Ernzerhof (PBE) functional
///
/// # Arguments
/// * `density` - electron density in electrons/bohr^3
/// * `density_gradient` - gradient of electron density in electrons/bohr^4
///
/// # Returns
/// * Exchange-correlation potential in Hartree
pub fn calculate_gga_potential(density: f64, density_gradient: f64) -> f64 {
    // Start with LDA potential
    let lda_potential = calculate_lda_potential(density);

    // Prevent division by zero or negative density
    if density <= 1e-12 {
        return lda_potential;
    }

    // Calculate the dimensionless gradient s = |∇ρ|/(2(3π²)^(1/3) ρ^(4/3))
    let kf = (3.0 * std::f64::consts::PI.powi(2) * density).powf(1.0 / 3.0);
    let s = density_gradient.abs() / (2.0 * kf * density);

    // Parameters for PBE functional
    let kappa = 0.804;
    let mu = 0.21951;

    // Enhancement factor
    let enhancement = 1.0 + kappa - kappa / (1.0 + mu * s * s / kappa);

    // Add the gradient correction to the LDA potential
    let gradient_correction = lda_potential * (enhancement - 1.0);

    lda_potential + gradient_correction
}

/// Calculate the Hedin-Lundqvist exchange-correlation potential with complex self-energy
///
/// # Arguments
/// * `density` - electron density in electrons/bohr^3
/// * `energy` - electron energy in Hartree
///
/// # Returns
/// * Complex exchange-correlation potential in Hartree
pub fn calculate_hedin_lundqvist_potential(density: f64, energy: f64) -> (f64, f64) {
    // Prevent division by zero or negative density
    if density <= 1e-12 {
        return (0.0, 0.0);
    }

    // Start with real part (LDA potential)
    let real_part = calculate_lda_potential(density);

    // Calculate plasma frequency ω_p = sqrt(4πρ)
    let plasma_freq = (4.0 * std::f64::consts::PI * density).sqrt();

    // Imaginary part (optical potential)
    // Use the Hedin-Lundqvist form: -θ(E-E_F) * ω_p² / (4 * (E-E_F)²)
    // For E < E_F, the imaginary part is zero
    let imag_part = if energy > 0.0 {
        -plasma_freq.powi(2) / (4.0 * energy.powi(2))
    } else {
        0.0
    };

    (real_part, imag_part)
}

/// Calculate the Dirac-Hara exchange-correlation potential with complex self-energy
///
/// # Arguments
/// * `density` - electron density in electrons/bohr^3
/// * `energy` - electron energy in Hartree
///
/// # Returns
/// * Complex exchange-correlation potential in Hartree (real part, imaginary part)
pub fn calculate_dirac_hara_potential(density: f64, energy: f64) -> (f64, f64) {
    // Prevent division by zero or negative density
    if density <= 1e-12 {
        return (0.0, 0.0);
    }

    // Calculate exchange part only (Dirac exchange)
    // Vx = -3/2 * (3ρ/π)^(1/3)
    let exchange = -1.5 * (3.0 * density / std::f64::consts::PI).powf(1.0 / 3.0);

    // Dirac-Hara adds a small imaginary part similar to Hedin-Lundqvist
    // but with a different pre-factor
    let plasma_freq = (4.0 * std::f64::consts::PI * density).sqrt();

    // Imaginary part - weaker than HL
    let imag_part = if energy > 0.0 {
        -plasma_freq.powi(2) / (8.0 * energy.powi(2))
    } else {
        0.0
    };

    (exchange, imag_part)
}

/// Calculate the von Barth-Hedin exchange-correlation potential
///
/// # Arguments
/// * `density` - electron density in electrons/bohr^3
///
/// # Returns
/// * Exchange-correlation potential in Hartree
pub fn calculate_von_barth_potential(density: f64) -> f64 {
    // Prevent division by zero or negative density
    if density <= 1e-12 {
        return 0.0;
    }

    // Wigner-Seitz radius (r_s = (3/4π*ρ)^(1/3))
    let rs = (3.0 / (4.0 * std::f64::consts::PI * density)).powf(1.0 / 3.0);

    // Exchange part (Dirac exchange)
    let exchange = -0.458 / rs;

    // von Barth-Hedin correlation parameters
    let rp = 21.0;
    let cp = 0.045;
    let q = 0.0333 * rp / rs;
    let tan_q = q.atan();
    let x = q.powi(2) + q.powi(4) + q.powi(6) / 3.0;

    // Correlation potential
    let correlation =
        cp * ((x - tan_q) / (q.powi(6) / 3.0 - tan_q + (1.0 / 3.0) * std::f64::consts::PI));

    exchange + correlation
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_exchange_correlation_type() {
        // Test creation from strings
        assert_eq!(
            ExchangeCorrelationType::from_string("LDA").unwrap(),
            ExchangeCorrelationType::LDA
        );
        assert_eq!(
            ExchangeCorrelationType::from_string("gga").unwrap(),
            ExchangeCorrelationType::GGA
        );
        assert_eq!(
            ExchangeCorrelationType::from_string("Hedin-Lundqvist").unwrap(),
            ExchangeCorrelationType::HedinLundqvist
        );
        assert_eq!(
            ExchangeCorrelationType::from_string("Dirac-Hara").unwrap(),
            ExchangeCorrelationType::DiracHara
        );
        assert_eq!(
            ExchangeCorrelationType::from_string("von-barth").unwrap(),
            ExchangeCorrelationType::VonBarth
        );

        // Test invalid string
        assert!(ExchangeCorrelationType::from_string("invalid").is_err());

        // Test string conversion
        assert_eq!(ExchangeCorrelationType::LDA.as_str(), "LDA");
        assert_eq!(ExchangeCorrelationType::GGA.as_str(), "GGA");
        assert_eq!(
            ExchangeCorrelationType::HedinLundqvist.as_str(),
            "Hedin-Lundqvist"
        );
        assert_eq!(ExchangeCorrelationType::DiracHara.as_str(), "Dirac-Hara");
        assert_eq!(
            ExchangeCorrelationType::VonBarth.as_str(),
            "von Barth-Hedin"
        );
    }

    #[test]
    fn test_lda_potential() {
        // Test zero density
        let v_xc_zero = calculate_lda_potential(0.0);
        assert_eq!(v_xc_zero, 0.0);

        // Test metallic density (rs=2)
        let density_rs2 = 3.0 / (4.0 * std::f64::consts::PI * 2.0f64.powi(3));
        let v_xc_rs2 = calculate_lda_potential(density_rs2);

        // Expected value: should be negative and close to -0.4 Hartree
        assert!(v_xc_rs2 < 0.0);
        assert!(v_xc_rs2 > -1.0);

        // Test high density (rs=0.5)
        let density_rs05 = 3.0 / (4.0 * std::f64::consts::PI * 0.5f64.powi(3));
        let v_xc_rs05 = calculate_lda_potential(density_rs05);

        // Higher density should give more negative potential
        assert!(v_xc_rs05 < v_xc_rs2);
    }

    #[test]
    fn test_gga_potential() {
        // Test zero density and gradient
        let v_xc_zero = calculate_gga_potential(0.0, 0.0);
        assert_eq!(v_xc_zero, 0.0);

        // Test with density but zero gradient (should equal LDA)
        let density = 0.01;
        let v_lda = calculate_lda_potential(density);
        let v_gga_no_grad = calculate_gga_potential(density, 0.0);
        assert_relative_eq!(v_lda, v_gga_no_grad, epsilon = 1e-10);

        // Test with gradient
        let gradient = 0.005;
        let v_gga = calculate_gga_potential(density, gradient);

        // GGA potential should differ from LDA
        assert!(v_gga != v_lda);
    }

    #[test]
    fn test_hedin_lundqvist_potential() {
        // Test zero density
        let (re, im) = calculate_hedin_lundqvist_potential(0.0, 1.0);
        assert_eq!(re, 0.0);
        assert_eq!(im, 0.0);

        // Test with density but energy below Fermi level
        let density = 0.01;
        let energy_below_fermi = -0.1;
        let (re_below, im_below) = calculate_hedin_lundqvist_potential(density, energy_below_fermi);

        // Real part should be non-zero
        assert!(re_below != 0.0);
        // Imaginary part should be zero below Fermi
        assert_eq!(im_below, 0.0);

        // Test with energy above Fermi level
        let energy_above_fermi = 0.1;
        let (re_above, im_above) = calculate_hedin_lundqvist_potential(density, energy_above_fermi);

        // Real part should be the same as LDA
        assert_relative_eq!(re_above, calculate_lda_potential(density), epsilon = 1e-10);
        // Imaginary part should be negative above Fermi
        assert!(im_above < 0.0);
    }

    #[test]
    fn test_dirac_hara_potential() {
        // Test zero density
        let (re, im) = calculate_dirac_hara_potential(0.0, 1.0);
        assert_eq!(re, 0.0);
        assert_eq!(im, 0.0);

        // Test with density but energy below Fermi level
        let density = 0.01;
        let energy_below_fermi = -0.1;
        let (re_below, im_below) = calculate_dirac_hara_potential(density, energy_below_fermi);

        // Real part should be non-zero
        assert!(re_below != 0.0);
        // Imaginary part should be zero below Fermi
        assert_eq!(im_below, 0.0);

        // Test with energy above Fermi level
        let energy_above_fermi = 0.1;
        let (re_above, im_above) = calculate_dirac_hara_potential(density, energy_above_fermi);

        // Real part should be negative (exchange only)
        assert!(re_above < 0.0);
        // Imaginary part should be negative above Fermi
        assert!(im_above < 0.0);

        // Dirac-Hara has weaker imaginary part than Hedin-Lundqvist
        let (_, im_hl) = calculate_hedin_lundqvist_potential(density, energy_above_fermi);
        assert!(im_above > im_hl); // Less negative means weaker absorption
    }

    #[test]
    fn test_von_barth_potential() {
        // Test zero density
        let v_xc_zero = calculate_von_barth_potential(0.0);
        assert_eq!(v_xc_zero, 0.0);

        // Test metallic density (rs=2)
        let density_rs2 = 3.0 / (4.0 * std::f64::consts::PI * 2.0f64.powi(3));
        let v_xc_rs2 = calculate_von_barth_potential(density_rs2);

        // Expected value: should be negative
        assert!(v_xc_rs2 < 0.0);

        // Test high density (rs=0.5)
        let density_rs05 = 3.0 / (4.0 * std::f64::consts::PI * 0.5f64.powi(3));
        let v_xc_rs05 = calculate_von_barth_potential(density_rs05);

        // Higher density should give more negative potential
        assert!(v_xc_rs05 < v_xc_rs2);

        // von Barth should have different values from LDA
        assert!(v_xc_rs2 != calculate_lda_potential(density_rs2));
    }
}
