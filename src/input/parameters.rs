/*
MIT License with FEFF10 Attribution

Copyright (c) 2025 Ameyanagi

Based on or developed using Distribution: FEFF10.0
Copyright (c) 2020 FEFF Project, University of Washington and SLAC National Accelerator Laboratory.
All rights reserved.
*/

//! Parameter types for various FEFF input cards

use crate::atoms::Vector3D;

/// Information about a potential type
#[derive(Debug, Clone)]
pub struct PotentialInfo {
    /// Potential index (ipot)
    pub index: i32,
    /// Atomic number (Z)
    pub atomic_number: i32,
    /// Element symbol
    pub symbol: String,
    // Add other potential parameters as needed
}

/// Control parameters
#[derive(Debug, Clone, Default)]
pub struct ControlParams {
    /// Potential flag
    pub mpot: i32,
    /// Phase shift flag
    pub mphase: i32,
    /// Full multiple scattering flag
    pub mfms: i32,
    /// Path finder flag
    pub mpath: i32,
    /// FEFF calculation flag
    pub mfeff: i32,
    /// XAFS calculation flag
    pub mchi: i32,
}

/// Exchange parameters
#[derive(Debug, Clone, Default)]
pub struct ExchangeParams {
    /// Exchange correlation type
    pub ixc: i32,
    /// Constant real part of exchange correlation
    pub vr0: f64,
    /// Constant imaginary part of exchange correlation
    pub vi0: f64,
    /// Exchange correlation type for ground state
    pub ixc0: i32,
}

/// Hole parameters
#[derive(Debug, Clone, Default)]
pub struct HoleParams {
    /// Hole type (1=K, 2=L1, etc.)
    pub hole_type: i32,
    /// S0^2 factor
    pub s02: f64,
}

/// SCF parameters
#[derive(Debug, Clone, Default)]
pub struct ScfParams {
    /// SCF radius
    pub rfms: f64,
    /// Convergence acceleration factor
    pub ca: f64,
}

/// FMS parameters
#[derive(Debug, Clone, Default)]
pub struct FmsParams {
    /// FMS radius
    pub rfms: f64,
    /// Number of multiple scattering paths
    pub nmultiple: i32,
}

/// XANES parameters
#[derive(Debug, Clone, Default)]
pub struct XanesParams {
    /// Grid spacing
    pub rgrid: f64,
    /// Path selection criteria
    pub pcrit: f64,
    /// Energy shift
    pub edge_shift: f64,
}

/// Path parameters for RPATH card
#[derive(Debug, Clone, Default)]
pub struct RPathParams {
    /// Path selection criteria
    pub pcrit: f64,
    /// Maximum half-path length
    pub rmax: f64,
    /// Maximum number of paths
    pub nleg: i32,
}

/// Print parameters for PRINT card
#[derive(Debug, Clone, Default)]
pub struct PrintParams {
    /// Print level
    pub iprint: i32,
}

/// Correction parameters for CORRECTIONS card
#[derive(Debug, Clone, Default)]
pub struct CorrectionsParams {
    /// Real part of atomic correction
    pub real_correction: f64,
    /// Imaginary part of atomic correction
    pub imag_correction: f64,
    /// Correction selection flag
    pub icorr: i32,
}

/// S02 scaling parameter for S02 card
#[derive(Debug, Clone, Default)]
pub struct S02Params {
    /// Overall amplitude reduction factor
    pub s02: f64,
}

/// EDGE card parameters
#[derive(Debug, Clone, Default)]
pub struct EdgeParams {
    /// Edge type (K, L1, L2, L3, etc.)
    pub edge_type: String,
    /// Edge energy in eV (if specified)
    pub energy: Option<f64>,
}

/// DEBYE card parameters for thermal effects
#[derive(Debug, Clone, Default)]
pub struct DebyeParams {
    /// Debye temperature in K
    pub temp: f64,
    /// Debye-Waller factor
    pub debye_waller_factor: f64,
    /// Correlated Debye flag
    pub correlated_debye: bool,
}

/// LDOS card parameters for local density of states calculation
#[derive(Debug, Clone, Default)]
pub struct LdosParams {
    /// Lower energy bound in eV
    pub emin: f64,
    /// Upper energy bound in eV
    pub emax: f64,
    /// Energy step in eV
    pub estep: f64,
}

/// EXAFS card parameters
#[derive(Debug, Clone, Default)]
pub struct ExafsParams {
    /// Lower energy bound in eV
    pub emin: f64,
    /// Upper energy bound in eV
    pub emax: f64,
    /// Energy step in eV
    pub estep: f64,
}

/// DANES card parameters for differential ANES
#[derive(Debug, Clone, Default)]
pub struct DanesParams {
    /// Radius for DANES calculation
    pub radius: f64,
    /// Additional parameters
    pub parameters: Vec<f64>,
}

/// COREHOLE card parameters
#[derive(Debug, Clone, Default)]
pub struct CoreholeParams {
    /// Core hole treatment method (RPA, FSR, etc.)
    pub treatment: String,
    /// Additional parameters
    pub params: Vec<f64>,
}

/// POLARIZATION card parameters
#[derive(Debug, Clone, Default)]
pub struct PolarizationParams {
    /// X component of polarization vector
    pub x: f64,
    /// Y component of polarization vector
    pub y: f64,
    /// Z component of polarization vector
    pub z: f64,
    /// Ellipticity (if specified)
    pub ellipticity: Option<f64>,
}

/// REAL card parameters for real space grid
#[derive(Debug, Clone, Default)]
pub struct RealParams {
    /// Grid spacing
    pub spacing: f64,
    /// Grid size
    pub size: i32,
}

/// RECIPROCAL card parameters for reciprocal space grid
#[derive(Debug, Clone, Default)]
pub struct ReciprocalParams {
    /// Grid spacing
    pub spacing: f64,
    /// Grid size
    pub size: i32,
}

/// ELNES card parameters for electron energy loss near edge structure
#[derive(Debug, Clone, Default)]
pub struct ElnesParams {
    /// Maximum k-value for the calculation (in Å⁻¹)
    pub xkmax: f64,
    /// Step size of the upper part of the k-mesh
    pub xkstep: f64,
    /// Step size of the lower part of the k-mesh
    pub vixan: f64,
    /// Beam energy in keV
    pub beam_energy: f64,
    /// Whether to average over multiple orientations
    pub aver: Option<f64>,
    /// Controls relativistic cross-section calculation
    pub cross: Option<f64>,
    /// Whether to include relativistic effects
    pub relat: Option<f64>,
    /// Beam direction in crystal coordinates
    pub beam_direction: Vector3D,
    /// Collection semi-angle in mrad
    pub beta: f64,
    /// Convergence semi-angle in mrad
    pub alpha: f64,
    /// Number of radial integration points
    pub nr: i32,
    /// Number of angular integration points
    pub na: i32,
    /// Position of the detector (x,y angle in mrad)
    pub detector_position: Vector3D,
}

/// NRIXS card parameters for non-resonant inelastic x-ray scattering
#[derive(Debug, Clone, Default)]
pub struct NrixsParams {
    /// Controls how the momentum transfer is calculated
    pub nq: i32,
    /// X-component of momentum transfer vector or magnitude for spherical averaging
    pub qx: f64,
    /// Y-component of momentum transfer vector (ignored for spherical averaging)
    pub qy: f64,
    /// Z-component of momentum transfer vector (ignored for spherical averaging)
    pub qz: f64,
    /// Optional scalar parameter
    pub scalar: Option<f64>,
}

/// ELLIPTICITY card parameters for elliptical polarization
#[derive(Debug, Clone, Default)]
pub struct EllipticityParams {
    /// Ellipticity value (-1.0 to 1.0)
    pub ellipticity: f64,
    /// Incident beam direction vector
    pub beam_direction: Vector3D,
}

/// OPCONS card parameters for optical constants calculation
#[derive(Debug, Clone, Default)]
pub struct OpConsParams {
    /// Flag indicating that the OPCONS card is present
    pub enabled: bool,
}

/// TDLDA card parameters for time-dependent local density approximation
#[derive(Debug, Clone, Default)]
pub struct TdldaParams {
    /// Screening flag: 0=RPA, 1=TDA, 2=TDLDA
    pub iscreen: i32,
    /// Calculation flag: 0/1=SCF+XAS, 2=store matrix, 3=non-SCF from stored
    pub icalc: i32,
    /// Energy range below absorption edge (in eV)
    pub elow: f64,
    /// Energy range above absorption edge (in eV)
    pub ehigh: f64,
    /// Energy step (in eV)
    pub estep: f64,
    /// Broadening parameter (in eV)
    pub gamma: f64,
}

/// MULTIPOLE card parameters for multipole transitions
#[derive(Debug, Clone, Default)]
pub struct MultipoleParams {
    /// Maximum orbital angular momentum (l) for the final state
    pub lmax: i32,
    /// Highest order multipole calculated (1=dipole, 2=quadrupole, etc.)
    pub morder: i32,
    /// Enable spherical tensor calculation (0=no, 1=yes)
    pub tensor: i32,
}

/// SCREEN card parameters for self-energy corrections to core-hole screening
#[derive(Debug, Clone, Default)]
pub struct ScreenParams {
    /// Self-energy scheme (0=none, 1=HL, 2=DH)
    pub iself: i32,
    /// Core hole screening mode (0=unscreened, 1=screened, 2=partially screened)
    pub iscreen: i32,
    /// Self-energy real part pre-factor
    pub ca1: f64,
    /// Self-energy imaginary part pre-factor
    pub ci1: f64,
}

/// SPECTRAL card parameters for spectral function convolution
#[derive(Debug, Clone, Default)]
pub struct SpectralParams {
    /// Apply spectral convolution (0=no, 1=yes)
    pub ispect: i32,
    /// Sharpen spectra (0=no, 1=yes)
    pub ispsharp: i32,
    /// Parameter for Fermi level determination
    pub isprule: i32,
    /// Spectral energy range below Fermi level in eV
    pub emin: f64,
    /// Spectral energy range above Fermi level in eV  
    pub emax: f64,
    /// Energy grid spacing for spectral function in eV
    pub estep: f64,
}

/// DIMENSIONS card parameters for array dimensioning
#[derive(Debug, Clone, Default)]
pub struct DimensionsParams {
    /// Maximum l quantum number for waves and phase shifts
    pub nat: i32,
    /// Maximum atomic sites
    pub nph: i32,
    /// Maximum number of r-mesh points for potentials
    pub lx: i32,
    /// Maximum number of unique potentials
    pub npot: i32,
    /// Maximum number of paths
    pub nstat: i32,
}

/// RDINP card parameters for reading input from a different file
#[derive(Debug, Clone, Default)]
pub struct RdinpParams {
    /// File name to read input from
    pub file_name: String,
}

/// BANDSTRUCTURE card parameters for band structure calculations
#[derive(Debug, Clone, Default)]
pub struct BandstructureParams {
    /// Number of k-points
    pub nk: i32,
    /// Minimum energy relative to Fermi level (eV)
    pub emin: f64,
    /// Maximum energy relative to Fermi level (eV)
    pub emax: f64,
    /// Energy step (eV)
    pub estep: f64,
    /// k-point mesh flag (0=user-defined, 1=uniform)
    pub kmesh: i32,
    /// Symmetry flag (0=no symmetry, 1=use symmetry)
    pub symmetry: i32,
}

/// KMESH card parameters for k-space mesh in band structure calculations
#[derive(Debug, Clone, Default)]
pub struct KmeshParams {
    /// Number of k-points in x direction
    pub nx: i32,
    /// Number of k-points in y direction
    pub ny: i32,
    /// Number of k-points in z direction
    pub nz: i32,
    /// K-points defined explicitly in content
    pub kpoints: Vec<(f64, f64, f64)>,
}

/// RESTART card parameters for restarting from saved files
#[derive(Debug, Clone, Default)]
pub struct RestartParams {
    /// The module to restart from (e.g., "pot", "phase", "fms", etc.)
    pub module: String,
    /// Additional file name if needed
    pub file_name: Option<String>,
}

/// DOS card parameters for density of states calculations
#[derive(Debug, Clone, Default)]
pub struct DosParams {
    /// Lower energy bound in eV
    pub emin: f64,
    /// Upper energy bound in eV
    pub emax: f64,
    /// Energy step in eV
    pub estep: f64,
    /// Gaussian broadening in eV
    pub gamma: f64,
    /// Additional DOS calculation parameters
    pub params: Vec<f64>,
}

/// CIFS card parameters for crystallographic information
#[derive(Debug, Clone, Default)]
pub struct CifsParams {
    /// File name for CIF (Crystallographic Information File)
    pub file_name: String,
    /// Site number for the absorbing atom
    pub site_index: Option<i32>,
    /// Distance cutoff for cluster construction
    pub distance_cutoff: Option<f64>,
}
