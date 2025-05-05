/*
MIT License with FEFF10 Attribution

Copyright (c) 2025 Ameyanagi

Based on or developed using Distribution: FEFF10.0
Copyright (c) 2020 FEFF Project, University of Washington and SLAC National Accelerator Laboratory.
All rights reserved.
*/

//! FEFF input model representation

use super::card::Card;
use super::errors::{InputError, Result};
use super::parameters::*;

use crate::atoms::AtomicStructure;
use std::collections::HashMap;
use std::fs::File;
use std::io::Write;
use std::path::Path;

/// Main FEFF input data structure
#[derive(Debug, Default)]
pub struct FeffInput {
    /// Title of the calculation
    pub title: Option<String>,

    /// Atomic structure
    pub atomic_structure: Option<AtomicStructure>,

    /// Potentials
    pub potentials: HashMap<i32, PotentialInfo>,

    /// Control parameters
    pub control: Option<ControlParams>,

    /// Exchange parameters
    pub exchange: Option<ExchangeParams>,

    /// Hole parameters
    pub hole: Option<HoleParams>,

    /// SCF parameters
    pub scf: Option<ScfParams>,

    /// FMS parameters
    pub fms: Option<FmsParams>,

    /// XANES parameters
    pub xanes: Option<XanesParams>,

    /// Path parameters
    pub rpath: Option<RPathParams>,

    /// Print parameters
    pub print: Option<PrintParams>,

    /// Correction parameters
    pub corrections: Option<CorrectionsParams>,

    /// S02 scaling factor
    pub s02: Option<S02Params>,

    /// EDGE card parameters
    pub edge: Option<EdgeParams>,

    /// DEBYE card parameters
    pub debye: Option<DebyeParams>,

    /// LDOS card parameters
    pub ldos: Option<LdosParams>,

    /// EXAFS card parameters
    pub exafs: Option<ExafsParams>,

    /// DANES card parameters
    pub danes: Option<DanesParams>,

    /// COREHOLE card parameters
    pub corehole: Option<CoreholeParams>,

    /// POLARIZATION card parameters
    pub polarization: Option<PolarizationParams>,

    /// REAL space grid parameters
    pub real_grid: Option<RealParams>,

    /// RECIPROCAL space grid parameters
    pub reciprocal_grid: Option<ReciprocalParams>,

    /// ELNES card parameters
    pub elnes: Option<ElnesParams>,

    /// NRIXS card parameters
    pub nrixs: Option<NrixsParams>,

    /// ELLIPTICITY card parameters
    pub ellipticity: Option<EllipticityParams>,

    /// OPCONS card parameters
    pub opcons: Option<OpConsParams>,

    /// Unknown cards
    pub unknown_cards: Vec<Card>,
}

impl FeffInput {
    /// Create a new empty FEFF input
    pub fn new() -> Self {
        Self::default()
    }

    /// Validate the input
    pub fn validate(&self) -> Result<()> {
        // Check for required cards
        if self.atomic_structure.is_none() {
            return Err(InputError::MissingCard("ATOMS".to_string()));
        }

        // Validate atomic structure
        if let Some(atomic_structure) = &self.atomic_structure {
            if atomic_structure.atom_count() == 0 {
                return Err(InputError::InvalidStructure(
                    "Atomic structure is empty".to_string(),
                ));
            }
        }

        // Additional validation can be added here

        Ok(())
    }

    /// Write the input to a file
    pub fn write<P: AsRef<Path>>(&self, path: P) -> Result<()> {
        let file = File::create(path).map_err(InputError::IoError)?;
        let mut writer = std::io::BufWriter::new(file);

        // Write TITLE card if available
        if let Some(title) = &self.title {
            writeln!(writer, "TITLE {}", title).map_err(InputError::IoError)?;
            writeln!(writer).map_err(InputError::IoError)?;
        }

        // Write CONTROL card if available
        if let Some(control) = &self.control {
            writeln!(
                writer,
                "CONTROL {} {} {} {} {} {}",
                control.mpot,
                control.mphase,
                control.mfms,
                control.mpath,
                control.mfeff,
                control.mchi
            )
            .map_err(InputError::IoError)?;
            writeln!(writer).map_err(InputError::IoError)?;
        }

        // Write EXCHANGE card if available
        if let Some(exchange) = &self.exchange {
            writeln!(
                writer,
                "EXCHANGE {} {} {} {}",
                exchange.ixc, exchange.vr0, exchange.vi0, exchange.ixc0
            )
            .map_err(InputError::IoError)?;
            writeln!(writer).map_err(InputError::IoError)?;
        }

        // Write HOLE card if available
        if let Some(hole) = &self.hole {
            writeln!(writer, "HOLE {} {}", hole.hole_type, hole.s02)
                .map_err(InputError::IoError)?;
            writeln!(writer).map_err(InputError::IoError)?;
        }

        // Write SCF card if available
        if let Some(scf) = &self.scf {
            writeln!(writer, "SCF {} {}", scf.rfms, scf.ca).map_err(InputError::IoError)?;
            writeln!(writer).map_err(InputError::IoError)?;
        }

        // Write FMS card if available
        if let Some(fms) = &self.fms {
            writeln!(writer, "FMS {} {}", fms.rfms, fms.nmultiple).map_err(InputError::IoError)?;
            writeln!(writer).map_err(InputError::IoError)?;
        }

        // Write XANES card if available
        if let Some(xanes) = &self.xanes {
            writeln!(
                writer,
                "XANES {} {} {}",
                xanes.rgrid, xanes.pcrit, xanes.edge_shift
            )
            .map_err(InputError::IoError)?;
            writeln!(writer).map_err(InputError::IoError)?;
        }

        // Write RPATH card if available
        if let Some(rpath) = &self.rpath {
            writeln!(
                writer,
                "RPATH {} {} {}",
                rpath.pcrit, rpath.rmax, rpath.nleg
            )
            .map_err(InputError::IoError)?;
            writeln!(writer).map_err(InputError::IoError)?;
        }

        // Write PRINT card if available
        if let Some(print) = &self.print {
            writeln!(writer, "PRINT {}", print.iprint).map_err(InputError::IoError)?;
            writeln!(writer).map_err(InputError::IoError)?;
        }

        // Write CORRECTIONS card if available
        if let Some(corrections) = &self.corrections {
            writeln!(
                writer,
                "CORRECTIONS {} {} {}",
                corrections.real_correction, corrections.imag_correction, corrections.icorr
            )
            .map_err(InputError::IoError)?;
            writeln!(writer).map_err(InputError::IoError)?;
        }

        // Write S02 card if available
        if let Some(s02) = &self.s02 {
            writeln!(writer, "S02 {}", s02.s02).map_err(InputError::IoError)?;
            writeln!(writer).map_err(InputError::IoError)?;
        }

        // Write EDGE card if available
        if let Some(edge) = &self.edge {
            let energy_str = if let Some(energy) = edge.energy {
                format!(" {}", energy)
            } else {
                String::new()
            };
            writeln!(writer, "EDGE {}{}", edge.edge_type, energy_str)
                .map_err(InputError::IoError)?;
            writeln!(writer).map_err(InputError::IoError)?;
        }

        // Write DEBYE card if available
        if let Some(debye) = &self.debye {
            writeln!(
                writer,
                "DEBYE {} {} {}",
                debye.temp,
                debye.debye_waller_factor,
                if debye.correlated_debye { 1 } else { 0 }
            )
            .map_err(InputError::IoError)?;
            writeln!(writer).map_err(InputError::IoError)?;
        }

        // Write LDOS card if available
        if let Some(ldos) = &self.ldos {
            writeln!(writer, "LDOS {} {} {}", ldos.emin, ldos.emax, ldos.estep)
                .map_err(InputError::IoError)?;
            writeln!(writer).map_err(InputError::IoError)?;
        }

        // Write EXAFS card if available
        if let Some(exafs) = &self.exafs {
            writeln!(
                writer,
                "EXAFS {} {} {}",
                exafs.emin, exafs.emax, exafs.estep
            )
            .map_err(InputError::IoError)?;
            writeln!(writer).map_err(InputError::IoError)?;
        }

        // Write DANES card if available
        if let Some(danes) = &self.danes {
            let mut line = format!("DANES {}", danes.radius);
            for param in &danes.parameters {
                line.push_str(&format!(" {}", param));
            }
            writeln!(writer, "{}", line).map_err(InputError::IoError)?;
            writeln!(writer).map_err(InputError::IoError)?;
        }

        // Write COREHOLE card if available
        if let Some(corehole) = &self.corehole {
            let mut line = format!("COREHOLE {}", corehole.treatment);
            for param in &corehole.params {
                line.push_str(&format!(" {}", param));
            }
            writeln!(writer, "{}", line).map_err(InputError::IoError)?;
            writeln!(writer).map_err(InputError::IoError)?;
        }

        // Write POLARIZATION card if available
        if let Some(pol) = &self.polarization {
            let ellipticity_str = if let Some(ellipticity) = pol.ellipticity {
                format!(" {}", ellipticity)
            } else {
                String::new()
            };
            writeln!(
                writer,
                "POLARIZATION {} {} {}{}",
                pol.x, pol.y, pol.z, ellipticity_str
            )
            .map_err(InputError::IoError)?;
            writeln!(writer).map_err(InputError::IoError)?;
        }

        // Write REAL card if available
        if let Some(real) = &self.real_grid {
            writeln!(writer, "REAL {} {}", real.spacing, real.size).map_err(InputError::IoError)?;
            writeln!(writer).map_err(InputError::IoError)?;
        }

        // Write RECIPROCAL card if available
        if let Some(recip) = &self.reciprocal_grid {
            writeln!(writer, "RECIPROCAL {} {}", recip.spacing, recip.size)
                .map_err(InputError::IoError)?;
            writeln!(writer).map_err(InputError::IoError)?;
        }

        // Write ELNES card if available
        if let Some(elnes) = &self.elnes {
            // First line: ELNES [xkmax xkstep vixan]
            writeln!(
                writer,
                "ELNES {} {} {}",
                elnes.xkmax, elnes.xkstep, elnes.vixan
            )
            .map_err(InputError::IoError)?;

            // Second line: E [aver [cross [relat]]]
            let mut line = format!("{}", elnes.beam_energy);
            if let Some(aver) = elnes.aver {
                line.push_str(&format!(" {}", aver));
                if let Some(cross) = elnes.cross {
                    line.push_str(&format!(" {}", cross));
                    if let Some(relat) = elnes.relat {
                        line.push_str(&format!(" {}", relat));
                    }
                }
            }
            writeln!(writer, "{}", line).map_err(InputError::IoError)?;

            // Third line: kx ky kz
            let beam_dir = &elnes.beam_direction;
            writeln!(writer, "{} {} {}", beam_dir.x, beam_dir.y, beam_dir.z)
                .map_err(InputError::IoError)?;

            // Fourth line: β α
            writeln!(writer, "{} {}", elnes.beta, elnes.alpha).map_err(InputError::IoError)?;

            // Fifth line: nr na
            writeln!(writer, "{} {}", elnes.nr, elnes.na).map_err(InputError::IoError)?;

            // Sixth line: dx dy
            let detector = &elnes.detector_position;
            writeln!(writer, "{} {}", detector.x, detector.y).map_err(InputError::IoError)?;

            writeln!(writer).map_err(InputError::IoError)?;
        }

        // Write NRIXS card if available
        if let Some(nrixs) = &self.nrixs {
            let mut line = format!("NRIXS {}", nrixs.nq);

            if nrixs.nq < 0 {
                // Spherical averaging - only need magnitude
                line.push_str(&format!(" {}", nrixs.qx));
            } else {
                // Specific q-vector - need all components
                line.push_str(&format!(" {} {} {}", nrixs.qx, nrixs.qy, nrixs.qz));
            }

            // Add optional scalar parameter if present
            if let Some(scalar) = nrixs.scalar {
                line.push_str(&format!(" {}", scalar));
            }

            writeln!(writer, "{}", line).map_err(InputError::IoError)?;
            writeln!(writer).map_err(InputError::IoError)?;
        }

        // Write ELLIPTICITY card if available
        if let Some(ellipticity) = &self.ellipticity {
            let beam_dir = &ellipticity.beam_direction;
            writeln!(
                writer,
                "ELLIPTICITY {} {} {} {}",
                ellipticity.ellipticity, beam_dir.x, beam_dir.y, beam_dir.z
            )
            .map_err(InputError::IoError)?;
            writeln!(writer).map_err(InputError::IoError)?;
        }

        // Write OPCONS card if enabled
        if let Some(opcons) = &self.opcons {
            if opcons.enabled {
                writeln!(writer, "OPCONS").map_err(InputError::IoError)?;
                writeln!(writer).map_err(InputError::IoError)?;
            }
        }

        // Write POTENTIALS card
        if !self.potentials.is_empty() {
            writeln!(writer, "POTENTIALS").map_err(InputError::IoError)?;

            // Sort by potential index
            let mut pot_indices: Vec<i32> = self.potentials.keys().copied().collect();
            pot_indices.sort();

            for ipot in pot_indices {
                if let Some(pot_info) = self.potentials.get(&ipot) {
                    writeln!(
                        writer,
                        "{} {} {}",
                        pot_info.index, pot_info.atomic_number, pot_info.symbol
                    )
                    .map_err(InputError::IoError)?;
                }
            }
            writeln!(writer).map_err(InputError::IoError)?;
        }

        // Write ATOMS card
        if let Some(atomic_structure) = &self.atomic_structure {
            writeln!(writer, "ATOMS").map_err(InputError::IoError)?;
            writeln!(writer, "* Cartesian coordinates").map_err(InputError::IoError)?;

            for atom in atomic_structure.atoms().iter() {
                let position = atom.position();
                let ipot = atom.potential_type();

                // Get element symbol - prefer the atom's own symbol
                let symbol = atom.symbol();

                writeln!(
                    writer,
                    "{} {:.6} {:.6} {:.6} {}",
                    ipot, position.x, position.y, position.z, symbol
                )
                .map_err(InputError::IoError)?;
            }
            writeln!(writer).map_err(InputError::IoError)?;
        }

        // Write any unknown cards
        for card in &self.unknown_cards {
            writeln!(writer, "{}", card.name).map_err(InputError::IoError)?;
            for line in &card.content {
                writeln!(writer, "{}", line).map_err(InputError::IoError)?;
            }
            writeln!(writer).map_err(InputError::IoError)?;
        }

        Ok(())
    }
}
