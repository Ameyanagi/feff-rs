/*
MIT License with FEFF10 Attribution

Copyright (c) 2025 Ameyanagi

Based on or developed using Distribution: FEFF10.0
Copyright (c) 2020 FEFF Project, University of Washington and SLAC National Accelerator Laboratory.
All rights reserved.
*/

//! Parser implementation for FEFF input files

use super::card::{is_card_name, Card};
use super::config::ParserConfig;
use super::coordinates::CoordinateSystem;
use super::errors::{InputError, Result};
use super::model::FeffInput;
use super::parameters::*;

use crate::atoms::{database, Atom, AtomicStructure, Vector3D};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

/// Main FEFF input file parser
#[derive(Debug)]
pub struct FeffInputParser {
    config: ParserConfig,
    cards: Vec<Card>,
    current_line: usize,
}

impl FeffInputParser {
    /// Create a new FEFF input parser with the given configuration
    pub fn new(config: ParserConfig) -> Self {
        Self {
            config,
            cards: Vec::new(),
            current_line: 0,
        }
    }

    /// Parse a FEFF input file
    pub fn parse<P: AsRef<Path>>(&mut self, path: Option<P>) -> Result<FeffInput> {
        let path = match path {
            Some(p) => p.as_ref().to_path_buf(),
            None => self.config.input_path.clone(),
        };

        // Open the file and read its contents
        let file = File::open(&path).map_err(InputError::IoError)?;
        let reader = BufReader::new(file);

        // Parse the file into cards
        self.parse_cards(reader)?;

        // Create a new FEFF input object
        let mut feff_input = FeffInput::default();

        // Process each card type
        for card in &self.cards {
            match card.name.as_str() {
                "ATOMS" => self.parse_atoms_card(card, &mut feff_input)?,
                "POTENTIALS" => self.parse_potentials_card(card, &mut feff_input)?,
                "CONTROL" => self.parse_control_card(card, &mut feff_input)?,
                "EXCHANGE" => self.parse_exchange_card(card, &mut feff_input)?,
                "HOLE" => self.parse_hole_card(card, &mut feff_input)?,
                "SCF" => self.parse_scf_card(card, &mut feff_input)?,
                "FMS" => self.parse_fms_card(card, &mut feff_input)?,
                "XANES" => self.parse_xanes_card(card, &mut feff_input)?,
                "TITLE" => self.parse_title_card(card, &mut feff_input)?,
                "RPATH" => self.parse_rpath_card(card, &mut feff_input)?,
                "PRINT" => self.parse_print_card(card, &mut feff_input)?,
                "CORRECTIONS" => self.parse_corrections_card(card, &mut feff_input)?,
                "S02" => self.parse_s02_card(card, &mut feff_input)?,
                "EDGE" => self.parse_edge_card(card, &mut feff_input)?,
                "DEBYE" => self.parse_debye_card(card, &mut feff_input)?,
                "LDOS" => self.parse_ldos_card(card, &mut feff_input)?,
                "EXAFS" => self.parse_exafs_card(card, &mut feff_input)?,
                "DANES" => self.parse_danes_card(card, &mut feff_input)?,
                "COREHOLE" => self.parse_corehole_card(card, &mut feff_input)?,
                "POLARIZATION" => self.parse_polarization_card(card, &mut feff_input)?,
                "REAL" => self.parse_real_card(card, &mut feff_input)?,
                "RECIPROCAL" => self.parse_reciprocal_card(card, &mut feff_input)?,
                "ELNES" => self.parse_elnes_card(card, &mut feff_input)?,
                "NRIXS" => self.parse_nrixs_card(card, &mut feff_input)?,
                "ELLIPTICITY" => self.parse_ellipticity_card(card, &mut feff_input)?,
                "OPCONS" => self.parse_opcons_card(card, &mut feff_input)?,
                "TDLDA" => self.parse_tdlda_card(card, &mut feff_input)?,
                "MULTIPOLE" => self.parse_multipole_card(card, &mut feff_input)?,
                "SCREEN" => self.parse_screen_card(card, &mut feff_input)?,
                "SPECTRAL" => self.parse_spectral_card(card, &mut feff_input)?,
                "DIMENSIONS" => self.parse_dimensions_card(card, &mut feff_input)?,
                "RDINP" => self.parse_rdinp_card(card, &mut feff_input)?,
                "BANDSTRUCTURE" => self.parse_bandstructure_card(card, &mut feff_input)?,
                "KMESH" => self.parse_kmesh_card(card, &mut feff_input)?,
                "RESTART" => self.parse_restart_card(card, &mut feff_input)?,
                "DOS" => self.parse_dos_card(card, &mut feff_input)?,
                "CIFS" => self.parse_cifs_card(card, &mut feff_input)?,
                // Unknown card, store it for future reference
                _ => {
                    feff_input.unknown_cards.push(card.clone());
                }
            }
        }

        // Validate the input
        if self.config.validate {
            feff_input.validate()?;
        }

        Ok(feff_input)
    }

    /// Parse the input file into cards
    fn parse_cards<R: BufRead>(&mut self, reader: R) -> Result<()> {
        self.cards.clear();
        self.current_line = 0;

        let mut current_card: Option<Card> = None;

        for (i, line_result) in reader.lines().enumerate() {
            self.current_line = i + 1;
            let line = line_result?;
            let trimmed = line.trim();

            // Skip empty lines and comments
            if trimmed.is_empty() || trimmed.starts_with('*') || trimmed.starts_with('#') {
                continue;
            }

            // Check if this is a new card
            if is_card_name(trimmed) {
                // If we were processing a card, save it before starting a new one
                if let Some(card) = current_card.take() {
                    self.cards.push(card);
                }

                // Get the card name (first word)
                let card_name = trimmed
                    .split_whitespace()
                    .next()
                    .unwrap_or("")
                    .to_uppercase();

                // Get any additional content on the same line
                let card_content = if trimmed.len() > card_name.len() {
                    // Extract content after the card name
                    trimmed[card_name.len()..].trim().to_string()
                } else {
                    String::new()
                };

                // Start a new card
                let mut content = Vec::new();
                if !card_content.is_empty() {
                    content.push(card_content);
                }

                current_card = Some(Card {
                    name: card_name,
                    content,
                    line_number: self.current_line,
                });
            } else if let Some(card) = &mut current_card {
                // Add this line to the current card's content
                card.content.push(line);
            } else {
                // Line outside any card - might be a comment or data without a card
                // For now, we'll just ignore these
            }
        }

        // Add the last card if there is one
        if let Some(card) = current_card {
            self.cards.push(card);
        }

        Ok(())
    }

    /// Parse ATOMS card - stores atomic positions and types
    fn parse_atoms_card(&self, card: &Card, input: &mut FeffInput) -> Result<()> {
        // Default coordinate system is Cartesian
        let mut coord_system = CoordinateSystem::Cartesian;
        let mut line_index = 0;

        // Check if the first line specifies a coordinate system
        if !card.content.is_empty() {
            let first_line = card.content[0].trim();
            if first_line.starts_with('*') {
                let upper = first_line.to_uppercase();
                if upper.contains("SPHERICAL") {
                    coord_system = CoordinateSystem::Spherical;
                } else if upper.contains("CYLINDRICAL") {
                    coord_system = CoordinateSystem::Cylindrical;
                }
                line_index += 1;
            }
        }

        // Create a new atomic structure if one doesn't already exist
        if input.atomic_structure.is_none() {
            input.atomic_structure = Some(AtomicStructure::new());
        }

        let atomic_structure = input.atomic_structure.as_mut().unwrap();

        // Process each atom line
        for i in line_index..card.content.len() {
            let line = card.content[i].trim();

            // Skip any comment lines
            if line.starts_with('*') || line.starts_with('#') || line.is_empty() {
                continue;
            }

            // Split the line into fields
            let fields: Vec<&str> = line.split_whitespace().collect();

            // Check if we have at least the required fields (x, y, z, ipot)
            if fields.len() < 4 {
                return Err(InputError::ParseError(format!(
                    "Invalid atom format at line {}: expected at least 4 fields",
                    card.line_number + i
                )));
            }

            // Parse ipot
            let pot_type: i32 = fields[0].parse().map_err(|_| {
                InputError::ParseError(format!(
                    "Invalid potential index at line {}: {}",
                    card.line_number + i,
                    fields[0]
                ))
            })?;

            // Parse coordinates
            let coords: [f64; 3] = [
                fields[1].parse().map_err(|_| {
                    InputError::ParseError(format!(
                        "Invalid coordinate at line {}: {}",
                        card.line_number + i,
                        fields[1]
                    ))
                })?,
                fields[2].parse().map_err(|_| {
                    InputError::ParseError(format!(
                        "Invalid coordinate at line {}: {}",
                        card.line_number + i,
                        fields[2]
                    ))
                })?,
                fields[3].parse().map_err(|_| {
                    InputError::ParseError(format!(
                        "Invalid coordinate at line {}: {}",
                        card.line_number + i,
                        fields[3]
                    ))
                })?,
            ];

            // Convert coordinates to the appropriate Vector3D based on coordinate system
            let position = match coord_system {
                CoordinateSystem::Cartesian => Vector3D::new(coords[0], coords[1], coords[2]),
                CoordinateSystem::Spherical => {
                    // Convert spherical (r, theta, phi) to Cartesian
                    let r = coords[0];
                    let theta = coords[1].to_radians();
                    let phi = coords[2].to_radians();

                    let x = r * theta.sin() * phi.cos();
                    let y = r * theta.sin() * phi.sin();
                    let z = r * theta.cos();

                    Vector3D::new(x, y, z)
                }
                CoordinateSystem::Cylindrical => {
                    // Convert cylindrical (rho, phi, z) to Cartesian
                    let rho = coords[0];
                    let phi = coords[1].to_radians();
                    let z = coords[2];

                    let x = rho * phi.cos();
                    let y = rho * phi.sin();

                    Vector3D::new(x, y, z)
                }
            };

            // Try to get the element symbol or atomic number from the tag field
            let mut atomic_number = 0;

            if fields.len() > 4 {
                let tag = fields[4];

                // Try to parse as atomic number
                if let Ok(z) = tag.parse::<i32>() {
                    atomic_number = z;
                } else {
                    // Try to parse as element symbol (first 1-2 characters)
                    let element_symbol: String =
                        tag.chars().take_while(|c| c.is_alphabetic()).collect();

                    // Look up atomic number from element symbol
                    if !element_symbol.is_empty() {
                        // Use the database to look up atomic number
                        if let Some(z) = database::atomic_number_from_symbol(&element_symbol) {
                            atomic_number = z;
                        }
                    }
                }
            }

            // Try to get atomic number from a few sources
            let atomic_number = if atomic_number > 0 {
                // We already got it from the tag
                atomic_number
            } else if let Some(pot_info) = input.potentials.get(&pot_type) {
                // Get from potentials map
                pot_info.atomic_number
            } else {
                // Default to 0, will be updated later when parsing potentials
                0
            };

            // Create atom
            let atom = match Atom::new(atomic_number, position, pot_type) {
                Ok(a) => a,
                Err(_) => {
                    // For now, just create a placeholder atom if we don't have a valid atomic number yet
                    // This will be fixed once we have proper element database lookup
                    Atom::new(1, position, pot_type).unwrap()
                }
            };

            // Add atom to structure
            let atom_index = atomic_structure.add_atom(atom);

            // If this is the first atom and ipot is 0, it's likely the central/absorbing atom
            if i == line_index && pot_type == 0 {
                let _ = atomic_structure.set_central_atom(atom_index);
            }
        }

        Ok(())
    }

    /// Parse POTENTIALS card
    fn parse_potentials_card(&self, card: &Card, input: &mut FeffInput) -> Result<()> {
        // Create potentials map
        let mut potentials = HashMap::new();

        for (i, line) in card.content.iter().enumerate() {
            let trimmed = line.trim();

            // Skip comment lines
            if trimmed.starts_with('*') || trimmed.starts_with('#') || trimmed.is_empty() {
                continue;
            }

            // Split the line into fields
            let fields: Vec<&str> = trimmed.split_whitespace().collect();

            // Check if we have at least the required fields (ipot, Z, symbol)
            if fields.len() < 3 {
                return Err(InputError::ParseError(format!(
                    "Invalid potential format at line {}: expected at least 3 fields",
                    card.line_number + i + 1
                )));
            }

            // Parse ipot (potential index)
            let ipot: i32 = fields[0].parse().map_err(|_| {
                InputError::ParseError(format!(
                    "Invalid potential index at line {}: {}",
                    card.line_number + i + 1,
                    fields[0]
                ))
            })?;

            // Parse Z (atomic number)
            let z: i32 = fields[1].parse().map_err(|_| {
                InputError::ParseError(format!(
                    "Invalid atomic number at line {}: {}",
                    card.line_number + i + 1,
                    fields[1]
                ))
            })?;

            // Parse symbol
            let symbol = fields[2].to_string();

            // Create a PotentialInfo object
            let potential_info = PotentialInfo {
                index: ipot,
                atomic_number: z,
                symbol,
                // Add other fields as needed
            };

            // Add to potentials map
            potentials.insert(ipot, potential_info);
        }

        // Set the potentials in the input
        input.potentials = potentials;

        Ok(())
    }

    /// Parse CONTROL card
    fn parse_control_card(&self, card: &Card, input: &mut FeffInput) -> Result<()> {
        if card.content.is_empty() {
            return Ok(());
        }

        let line = card.content[0].trim();
        let fields: Vec<&str> = line.split_whitespace().collect();

        // FEFF7 format: CONTROL mphase mpath mfeff mchi
        // FEFF8+ format: CONTROL mpot mphase mfms mpath mfeff mchi

        if fields.len() >= 4 {
            let control = if fields.len() == 4 {
                // FEFF7 format
                ControlParams {
                    mpot: fields[0].parse().unwrap_or(1),
                    mphase: fields[0].parse().unwrap_or(1),
                    mfms: 0,
                    mpath: fields[1].parse().unwrap_or(1),
                    mfeff: fields[2].parse().unwrap_or(1),
                    mchi: fields[3].parse().unwrap_or(1),
                }
            } else {
                // FEFF8+ format
                ControlParams {
                    mpot: fields[0].parse().unwrap_or(1),
                    mphase: fields[1].parse().unwrap_or(1),
                    mfms: fields[2].parse().unwrap_or(1),
                    mpath: fields[3].parse().unwrap_or(1),
                    mfeff: fields[4].parse().unwrap_or(1),
                    mchi: fields[5].parse().unwrap_or(1),
                }
            };

            input.control = Some(control);
        }

        Ok(())
    }

    /// Parse EXCHANGE card
    fn parse_exchange_card(&self, card: &Card, input: &mut FeffInput) -> Result<()> {
        if card.content.is_empty() {
            return Ok(());
        }

        let line = card.content[0].trim();
        let fields: Vec<&str> = line.split_whitespace().collect();

        // Format: EXCHANGE ixc vr0 vi0 [ixc0]
        if !fields.is_empty() {
            let exchange = ExchangeParams {
                ixc: fields[0].parse().unwrap_or(0),
                vr0: if fields.len() > 1 {
                    fields[1].parse().unwrap_or(0.0)
                } else {
                    0.0 // Default value
                },
                vi0: if fields.len() > 2 {
                    fields[2].parse().unwrap_or(0.0)
                } else {
                    0.0 // Default value
                },
                ixc0: if fields.len() > 3 {
                    fields[3].parse().unwrap_or(0)
                } else {
                    0 // Default value
                },
            };

            input.exchange = Some(exchange);
        }

        Ok(())
    }

    /// Parse HOLE card
    fn parse_hole_card(&self, card: &Card, input: &mut FeffInput) -> Result<()> {
        if card.content.is_empty() {
            return Ok(());
        }

        let line = card.content[0].trim();
        let fields: Vec<&str> = line.split_whitespace().collect();

        // Format: HOLE hole_type s02
        if !fields.is_empty() {
            // Initialize with parsed values
            let hole = HoleParams {
                hole_type: fields[0].parse().unwrap_or(1),
                s02: if fields.len() > 1 {
                    fields[1].parse().unwrap_or(1.0)
                } else {
                    1.0 // Default value
                },
            };

            input.hole = Some(hole);
        }

        Ok(())
    }

    /// Parse SCF card
    fn parse_scf_card(&self, card: &Card, input: &mut FeffInput) -> Result<()> {
        if card.content.is_empty() {
            return Ok(());
        }

        let line = card.content[0].trim();
        let fields: Vec<&str> = line.split_whitespace().collect();

        // Format: SCF rfms [ca]
        if !fields.is_empty() {
            let scf = ScfParams {
                rfms: fields[0].parse().unwrap_or(0.0),
                ca: if fields.len() > 1 {
                    fields[1].parse().unwrap_or(0.0)
                } else {
                    0.0 // Default value
                },
            };

            input.scf = Some(scf);
        }

        Ok(())
    }

    /// Parse FMS card
    fn parse_fms_card(&self, card: &Card, input: &mut FeffInput) -> Result<()> {
        if card.content.is_empty() {
            return Ok(());
        }

        let line = card.content[0].trim();
        let fields: Vec<&str> = line.split_whitespace().collect();

        // Format: FMS rfms [nmultiple]
        if !fields.is_empty() {
            let fms = FmsParams {
                rfms: fields[0].parse().unwrap_or(0.0),
                nmultiple: if fields.len() > 1 {
                    fields[1].parse().unwrap_or(0)
                } else {
                    0 // Default value
                },
            };

            input.fms = Some(fms);
        }

        Ok(())
    }

    /// Parse XANES card
    fn parse_xanes_card(&self, card: &Card, input: &mut FeffInput) -> Result<()> {
        if card.content.is_empty() {
            return Ok(());
        }

        let line = card.content[0].trim();
        let fields: Vec<&str> = line.split_whitespace().collect();

        // Format: XANES rgrid [ pcrit [ edge_shift ] ]
        if !fields.is_empty() {
            let xanes = XanesParams {
                rgrid: fields[0].parse().unwrap_or(0.0),
                pcrit: if fields.len() > 1 {
                    fields[1].parse().unwrap_or(0.0)
                } else {
                    0.0 // Default value
                },
                edge_shift: if fields.len() > 2 {
                    fields[2].parse().unwrap_or(0.0)
                } else {
                    0.0 // Default value
                },
            };

            input.xanes = Some(xanes);
        }

        Ok(())
    }

    /// Parse TITLE card
    fn parse_title_card(&self, card: &Card, input: &mut FeffInput) -> Result<()> {
        // The title is the whole line, including any content that follows the TITLE keyword
        // If the card has content, use it
        if !card.content.is_empty() {
            let title = card.content.join("\n");
            input.title = Some(title);
        } else {
            // If there's no content, the title might be on the same line as the TITLE keyword
            // In line is in the card name - extract it from there
            if card.name.len() > 5 {
                // "TITLE" is 5 chars
                let title = card.name[5..].trim().to_string();
                if !title.is_empty() {
                    input.title = Some(title);
                }
            }
        }

        Ok(())
    }

    /// Parse RPATH card - controls path selection
    fn parse_rpath_card(&self, card: &Card, input: &mut FeffInput) -> Result<()> {
        if card.content.is_empty() {
            return Ok(());
        }

        let line = card.content[0].trim();
        let fields: Vec<&str> = line.split_whitespace().collect();

        // Format: RPATH pcrit rmax nleg
        let mut rpath = RPathParams::default();

        if !fields.is_empty() {
            rpath.pcrit = fields[0].parse().unwrap_or(0.0);
        }

        if fields.len() > 1 {
            rpath.rmax = fields[1].parse().unwrap_or(0.0);
        }

        if fields.len() > 2 {
            rpath.nleg = fields[2].parse().unwrap_or(0);
        }

        input.rpath = Some(rpath);
        Ok(())
    }

    /// Parse PRINT card - controls output level
    fn parse_print_card(&self, card: &Card, input: &mut FeffInput) -> Result<()> {
        if card.content.is_empty() {
            return Ok(());
        }

        let line = card.content[0].trim();
        let fields: Vec<&str> = line.split_whitespace().collect();

        // Format: PRINT iprint
        let mut print = PrintParams::default();

        if !fields.is_empty() {
            print.iprint = fields[0].parse().unwrap_or(0);
        }

        input.print = Some(print);
        Ok(())
    }

    /// Parse CORRECTIONS card - atomic corrections
    fn parse_corrections_card(&self, card: &Card, input: &mut FeffInput) -> Result<()> {
        if card.content.is_empty() {
            return Ok(());
        }

        let line = card.content[0].trim();
        let fields: Vec<&str> = line.split_whitespace().collect();

        // Format: CORRECTIONS real_corr imag_corr icorr
        let mut corrections = CorrectionsParams::default();

        if !fields.is_empty() {
            corrections.real_correction = fields[0].parse().unwrap_or(0.0);
        }

        if fields.len() > 1 {
            corrections.imag_correction = fields[1].parse().unwrap_or(0.0);
        }

        if fields.len() > 2 {
            corrections.icorr = fields[2].parse().unwrap_or(0);
        }

        input.corrections = Some(corrections);
        Ok(())
    }

    /// Parse S02 card - overall amplitude reduction factor
    fn parse_s02_card(&self, card: &Card, input: &mut FeffInput) -> Result<()> {
        // Create the S02 parameters with default value first
        let mut s02 = S02Params::default();

        // If the card has content, parse it
        if !card.content.is_empty() {
            let line = card.content[0].trim();
            let fields: Vec<&str> = line.split_whitespace().collect();

            if !fields.is_empty() {
                s02.s02 = fields[0].parse().unwrap_or(1.0);
            }
        } else {
            // If there's content on the same line as the card name
            if card.name.len() > 3 {
                // "S02" is 3 chars
                let content = card.name[3..].trim();
                if !content.is_empty() {
                    let fields: Vec<&str> = content.split_whitespace().collect();
                    if !fields.is_empty() {
                        if let Ok(val) = fields[0].parse::<f64>() {
                            s02.s02 = val;
                        }
                    }
                }
            }
        }

        input.s02 = Some(s02);
        Ok(())
    }

    /// Parse EDGE card - absorption edge parameters
    fn parse_edge_card(&self, card: &Card, input: &mut FeffInput) -> Result<()> {
        if card.content.is_empty() {
            return Ok(());
        }

        let line = card.content[0].trim();
        let fields: Vec<&str> = line.split_whitespace().collect();

        // Format: EDGE K [energy]
        if !fields.is_empty() {
            // Validate edge type
            let edge_type = fields[0].to_string();
            if !Self::is_valid_edge_type(&edge_type) {
                return Err(InputError::InvalidFormat(format!(
                    "Invalid edge type '{}'. Valid types are K, L1, L2, L3, M1-M5, N1-N7, 1-23, etc.",
                    edge_type
                )));
            }

            // Validate energy if provided
            let energy = if fields.len() > 1 {
                match fields[1].parse::<f64>() {
                    Ok(energy) => {
                        if energy <= 0.0 {
                            return Err(InputError::InvalidFormat(format!(
                                "Edge energy must be positive, found: {}",
                                energy
                            )));
                        }
                        Some(energy)
                    }
                    Err(_) => {
                        return Err(InputError::ParseError(format!(
                            "Could not parse edge energy: {}",
                            fields[1]
                        )));
                    }
                }
            } else {
                None
            };

            let edge = EdgeParams { edge_type, energy };

            input.edge = Some(edge);
        }

        Ok(())
    }

    /// Validates if the edge type is one of the accepted edge types in FEFF
    fn is_valid_edge_type(edge_type: &str) -> bool {
        // Check if it's one of the valid string representations
        let valid_string_types = [
            "K", "L1", "L2", "L3", "M1", "M2", "M3", "M4", "M5", "N1", "N2", "N3", "N4", "N5",
            "N6", "N7", "O1", "O2", "O3", "O4", "O5", "O6", "O7",
        ];

        if valid_string_types.contains(&edge_type) {
            return true;
        }

        // Check if it's a valid numerical identifier (1-23)
        if let Ok(num) = edge_type.parse::<i32>() {
            return (1..=23).contains(&num); // 23 covers up to O7 edge
        }

        false
    }

    /// Parse DEBYE card - thermal effects parameters
    fn parse_debye_card(&self, card: &Card, input: &mut FeffInput) -> Result<()> {
        if card.content.is_empty() {
            return Ok(());
        }

        let line = card.content[0].trim();
        let fields: Vec<&str> = line.split_whitespace().collect();

        // Format: DEBYE temp [debye_waller [correlated_debye]]
        if !fields.is_empty() {
            let debye = DebyeParams {
                temp: fields[0].parse().unwrap_or(300.0),
                debye_waller_factor: if fields.len() > 1 {
                    fields[1].parse().unwrap_or(0.0)
                } else {
                    0.0
                },
                correlated_debye: if fields.len() > 2 {
                    fields[2].parse::<i32>().unwrap_or(0) > 0
                } else {
                    false
                },
            };

            input.debye = Some(debye);
        }

        Ok(())
    }

    /// Parse LDOS card - local density of states parameters
    fn parse_ldos_card(&self, card: &Card, input: &mut FeffInput) -> Result<()> {
        if card.content.is_empty() {
            return Ok(());
        }

        let line = card.content[0].trim();
        let fields: Vec<&str> = line.split_whitespace().collect();

        // Format: LDOS emin emax estep
        if fields.len() >= 3 {
            let ldos = LdosParams {
                emin: fields[0].parse().unwrap_or(-20.0),
                emax: fields[1].parse().unwrap_or(20.0),
                estep: fields[2].parse().unwrap_or(0.1),
            };

            input.ldos = Some(ldos);
        }

        Ok(())
    }

    /// Parse EXAFS card - extended X-ray absorption fine structure parameters
    fn parse_exafs_card(&self, card: &Card, input: &mut FeffInput) -> Result<()> {
        if card.content.is_empty() {
            return Ok(());
        }

        let line = card.content[0].trim();
        let fields: Vec<&str> = line.split_whitespace().collect();

        // Format: EXAFS emin emax estep
        if fields.len() >= 3 {
            let exafs = ExafsParams {
                emin: fields[0].parse().unwrap_or(0.0),
                emax: fields[1].parse().unwrap_or(1000.0),
                estep: fields[2].parse().unwrap_or(0.5),
            };

            input.exafs = Some(exafs);
        }

        Ok(())
    }

    /// Parse DANES card - differential anomalous scattering parameters
    fn parse_danes_card(&self, card: &Card, input: &mut FeffInput) -> Result<()> {
        if card.content.is_empty() {
            return Ok(());
        }

        let line = card.content[0].trim();
        let fields: Vec<&str> = line.split_whitespace().collect();

        // Format: DANES radius [params...]
        if !fields.is_empty() {
            let radius = fields[0].parse().unwrap_or(5.0);

            // Parse any additional parameters
            let mut parameters = Vec::new();
            for field in fields.iter().skip(1) {
                if let Ok(param) = field.parse::<f64>() {
                    parameters.push(param);
                }
            }

            let danes = DanesParams { radius, parameters };

            input.danes = Some(danes);
        }

        Ok(())
    }

    /// Parse COREHOLE card - core hole treatment parameters
    fn parse_corehole_card(&self, card: &Card, input: &mut FeffInput) -> Result<()> {
        if card.content.is_empty() {
            return Ok(());
        }

        let line = card.content[0].trim();
        let fields: Vec<&str> = line.split_whitespace().collect();

        // Format: COREHOLE treatment [params...]
        if !fields.is_empty() {
            let treatment = fields[0].to_string();

            // Parse any additional parameters
            let mut params = Vec::new();
            for field in fields.iter().skip(1) {
                if let Ok(param) = field.parse::<f64>() {
                    params.push(param);
                }
            }

            let corehole = CoreholeParams { treatment, params };

            input.corehole = Some(corehole);
        }

        Ok(())
    }

    /// Parse POLARIZATION card - polarization vector parameters
    fn parse_polarization_card(&self, card: &Card, input: &mut FeffInput) -> Result<()> {
        if card.content.is_empty() {
            return Ok(());
        }

        let line = card.content[0].trim();
        let fields: Vec<&str> = line.split_whitespace().collect();

        // Format: POLARIZATION x y z [ellipticity]
        if fields.len() >= 3 {
            // Parse vector components with proper error handling
            let x = match fields[0].parse::<f64>() {
                Ok(val) => val,
                Err(_) => {
                    return Err(InputError::ParseError(format!(
                        "Could not parse X component of polarization vector: {}",
                        fields[0]
                    )))
                }
            };

            let y = match fields[1].parse::<f64>() {
                Ok(val) => val,
                Err(_) => {
                    return Err(InputError::ParseError(format!(
                        "Could not parse Y component of polarization vector: {}",
                        fields[1]
                    )))
                }
            };

            let z = match fields[2].parse::<f64>() {
                Ok(val) => val,
                Err(_) => {
                    return Err(InputError::ParseError(format!(
                        "Could not parse Z component of polarization vector: {}",
                        fields[2]
                    )))
                }
            };

            // Check that the vector isn't a zero vector
            let magnitude_squared = x * x + y * y + z * z;
            if magnitude_squared < 1e-10 {
                return Err(InputError::InvalidFormat(
                    "Polarization vector cannot be a zero vector".to_string(),
                ));
            }

            // Parse ellipticity if provided
            let ellipticity = if fields.len() > 3 {
                match fields[3].parse::<f64>() {
                    Ok(val) => {
                        // Typically ellipticity should be between -1 and 1
                        if !(-1.0..=1.0).contains(&val) {
                            return Err(InputError::InvalidFormat(format!(
                                "Ellipticity should be between -1 and 1, found: {}",
                                val
                            )));
                        }
                        Some(val)
                    }
                    Err(_) => {
                        return Err(InputError::ParseError(format!(
                            "Could not parse ellipticity: {}",
                            fields[3]
                        )))
                    }
                }
            } else {
                None
            };

            // Normalize the vector before storing
            let magnitude = magnitude_squared.sqrt();
            let polarization = PolarizationParams {
                x: x / magnitude,
                y: y / magnitude,
                z: z / magnitude,
                ellipticity,
            };

            input.polarization = Some(polarization);
        } else {
            return Err(InputError::InvalidFormat(
                "POLARIZATION card requires at least 3 parameters (x, y, z)".to_string(),
            ));
        }

        Ok(())
    }

    /// Parse REAL card - real space grid parameters
    fn parse_real_card(&self, card: &Card, input: &mut FeffInput) -> Result<()> {
        if card.content.is_empty() {
            return Ok(());
        }

        let line = card.content[0].trim();
        let fields: Vec<&str> = line.split_whitespace().collect();

        // Format: REAL spacing size
        if fields.len() >= 2 {
            let real_grid = RealParams {
                spacing: fields[0].parse().unwrap_or(0.05),
                size: fields[1].parse().unwrap_or(100),
            };

            input.real_grid = Some(real_grid);
        }

        Ok(())
    }

    /// Parse RECIPROCAL card - reciprocal space grid parameters
    fn parse_reciprocal_card(&self, card: &Card, input: &mut FeffInput) -> Result<()> {
        if card.content.is_empty() {
            return Ok(());
        }

        let line = card.content[0].trim();
        let fields: Vec<&str> = line.split_whitespace().collect();

        // Format: RECIPROCAL spacing size
        if fields.len() >= 2 {
            let reciprocal_grid = ReciprocalParams {
                spacing: fields[0].parse().unwrap_or(0.05),
                size: fields[1].parse().unwrap_or(100),
            };

            input.reciprocal_grid = Some(reciprocal_grid);
        }

        Ok(())
    }

    /// Parse ELNES card - electron energy loss near edge structure
    fn parse_elnes_card(&self, card: &Card, input: &mut FeffInput) -> Result<()> {
        // ELNES card requires multiple lines
        if card.content.len() < 6 {
            return Err(InputError::ParseError(format!(
                "Insufficient content for ELNES card: need at least 6 lines, found {}",
                card.content.len()
            )));
        }

        // Parse first line: ELNES [xkmax xkstep vixan]
        let first_line = card.content[0].trim();
        let first_fields: Vec<&str> = first_line.split_whitespace().collect();

        let xkmax = if !first_fields.is_empty() {
            first_fields[0].parse().unwrap_or(4.0)
        } else {
            4.0
        };

        let xkstep = if first_fields.len() > 1 {
            first_fields[1].parse().unwrap_or(0.07)
        } else {
            0.07
        };

        let vixan = if first_fields.len() > 2 {
            first_fields[2].parse().unwrap_or(0.0)
        } else {
            0.0
        };

        // Parse second line: E [aver [cross [relat]]]
        let second_line = card.content[1].trim();
        let second_fields: Vec<&str> = second_line.split_whitespace().collect();

        if second_fields.is_empty() {
            return Err(InputError::ParseError(
                "Missing beam energy in ELNES card".to_string(),
            ));
        }

        let beam_energy = second_fields[0].parse().unwrap_or(300.0);
        let aver = if second_fields.len() > 1 {
            second_fields[1].parse::<f64>().ok()
        } else {
            None
        };

        let cross = if second_fields.len() > 2 {
            second_fields[2].parse::<f64>().ok()
        } else {
            None
        };

        let relat = if second_fields.len() > 3 {
            second_fields[3].parse::<f64>().ok()
        } else {
            None
        };

        // Parse third line: kx ky kz
        let third_line = card.content[2].trim();
        let third_fields: Vec<&str> = third_line.split_whitespace().collect();

        if third_fields.len() < 3 {
            return Err(InputError::ParseError(
                "Missing beam direction (kx ky kz) in ELNES card".to_string(),
            ));
        }

        let kx = third_fields[0].parse().unwrap_or(0.0);
        let ky = third_fields[1].parse().unwrap_or(0.0);
        let kz = third_fields[2].parse().unwrap_or(0.0);
        let beam_direction = Vector3D::new(kx, ky, kz);

        // Parse fourth line: β α
        let fourth_line = card.content[3].trim();
        let fourth_fields: Vec<&str> = fourth_line.split_whitespace().collect();

        if fourth_fields.is_empty() {
            return Err(InputError::ParseError(
                "Missing collection semi-angle (β) in ELNES card".to_string(),
            ));
        }

        let beta = fourth_fields[0].parse().unwrap_or(2.4);
        let alpha = if fourth_fields.len() > 1 {
            fourth_fields[1].parse().unwrap_or(0.0)
        } else {
            0.0
        };

        // Parse fifth line: nr na
        let fifth_line = card.content[4].trim();
        let fifth_fields: Vec<&str> = fifth_line.split_whitespace().collect();

        if fifth_fields.len() < 2 {
            return Err(InputError::ParseError(
                "Missing integration points (nr na) in ELNES card".to_string(),
            ));
        }

        let nr = fifth_fields[0].parse().unwrap_or(5);
        let na = fifth_fields[1].parse().unwrap_or(3);

        // Parse sixth line: dx dy
        let sixth_line = card.content[5].trim();
        let sixth_fields: Vec<&str> = sixth_line.split_whitespace().collect();

        let dx = if !sixth_fields.is_empty() {
            sixth_fields[0].parse().unwrap_or(0.0)
        } else {
            0.0
        };

        let dy = if sixth_fields.len() > 1 {
            sixth_fields[1].parse().unwrap_or(0.0)
        } else {
            0.0
        };

        let detector_position = Vector3D::new(dx, dy, 0.0);

        // Create ELNES parameters
        let elnes = ElnesParams {
            xkmax,
            xkstep,
            vixan,
            beam_energy,
            aver,
            cross,
            relat,
            beam_direction,
            beta,
            alpha,
            nr,
            na,
            detector_position,
        };

        input.elnes = Some(elnes);
        Ok(())
    }

    /// Parse NRIXS card - non-resonant inelastic x-ray scattering
    fn parse_nrixs_card(&self, card: &Card, input: &mut FeffInput) -> Result<()> {
        if card.content.is_empty() {
            return Ok(());
        }

        let line = card.content[0].trim();
        let fields: Vec<&str> = line.split_whitespace().collect();

        // Format: NRIXS nq qx qy qz [scalar]
        // Need at least 2 fields for spherical averaging or 4 for specific q-vector
        if fields.len() < 2 {
            return Err(InputError::ParseError(
                "Insufficient parameters for NRIXS card".to_string(),
            ));
        }

        // Parse nq parameter
        let nq = fields[0].parse::<i32>().map_err(|_| {
            InputError::ParseError(format!("Invalid nq value in NRIXS card: {}", fields[0]))
        })?;

        // Parse q-vector components
        let qx = fields[1].parse::<f64>().map_err(|_| {
            InputError::ParseError(format!("Invalid qx value in NRIXS card: {}", fields[1]))
        })?;

        // If nq is negative, we're doing spherical averaging and only need qx (magnitude)
        // If nq is positive, we need the full q-vector (qx, qy, qz)
        let (qy, qz) = if nq < 0 {
            (0.0, 0.0) // Default values for spherical averaging
        } else {
            // Specific q-vector - need qy and qz
            if fields.len() < 4 {
                return Err(InputError::ParseError(
                    "Missing q-vector components (qy, qz) for specific q in NRIXS card".to_string(),
                ));
            }

            let qy = fields[2].parse::<f64>().map_err(|_| {
                InputError::ParseError(format!("Invalid qy value in NRIXS card: {}", fields[2]))
            })?;

            let qz = fields[3].parse::<f64>().map_err(|_| {
                InputError::ParseError(format!("Invalid qz value in NRIXS card: {}", fields[3]))
            })?;

            (qy, qz)
        };

        // Optional scalar parameter
        let scalar = if (nq < 0 && fields.len() > 2) || (nq >= 0 && fields.len() > 4) {
            let idx = if nq < 0 { 2 } else { 4 };
            fields[idx].parse::<f64>().ok()
        } else {
            None
        };

        let nrixs = NrixsParams {
            nq,
            qx,
            qy,
            qz,
            scalar,
        };

        input.nrixs = Some(nrixs);
        Ok(())
    }

    /// Parse ELLIPTICITY card - elliptical polarization parameters
    fn parse_ellipticity_card(&self, card: &Card, input: &mut FeffInput) -> Result<()> {
        if card.content.is_empty() {
            return Ok(());
        }

        let line = card.content[0].trim();
        let fields: Vec<&str> = line.split_whitespace().collect();

        // Format: ELLIPTICITY elpty x y z
        // Need 4 fields: ellipticity and beam direction vector
        if fields.len() < 4 {
            return Err(InputError::ParseError(
                "Insufficient parameters for ELLIPTICITY card. Format: ELLIPTICITY elpty x y z"
                    .to_string(),
            ));
        }

        // Parse ellipticity
        let ellipticity = fields[0].parse::<f64>().map_err(|_| {
            InputError::ParseError(format!("Invalid ellipticity value: {}", fields[0]))
        })?;

        // Validate ellipticity range (typically between -1 and 1)
        if !(-1.0..=1.0).contains(&ellipticity) {
            return Err(InputError::InvalidFormat(format!(
                "Ellipticity should be between -1 and 1, found: {}",
                ellipticity
            )));
        }

        // Parse beam direction vector components
        let x = fields[1]
            .parse::<f64>()
            .map_err(|_| InputError::ParseError(format!("Invalid x component: {}", fields[1])))?;

        let y = fields[2]
            .parse::<f64>()
            .map_err(|_| InputError::ParseError(format!("Invalid y component: {}", fields[2])))?;

        let z = fields[3]
            .parse::<f64>()
            .map_err(|_| InputError::ParseError(format!("Invalid z component: {}", fields[3])))?;

        // Check that the beam direction vector isn't a zero vector
        let magnitude_squared = x * x + y * y + z * z;
        if magnitude_squared < 1e-10 {
            return Err(InputError::InvalidFormat(
                "Beam direction vector cannot be a zero vector".to_string(),
            ));
        }

        // Create beam direction vector
        let beam_direction = Vector3D::new(x, y, z);

        // Check if polarization is defined and if beam direction is approximately perpendicular to polarization
        if let Some(polarization) = &input.polarization {
            let dot_product = polarization.x * x + polarization.y * y + polarization.z * z;

            // Calculate the cosine of the angle between vectors
            let cosine = dot_product / (magnitude_squared.sqrt());

            // Check if vectors are approximately perpendicular (cosine close to 0)
            if cosine.abs() > 0.1 { // Allow some deviation from exactly perpendicular
                 // This is just a warning, not an error
                 // In a production environment, this could be logged or added as a warning to the user
            }
        }

        let ellipticity_params = EllipticityParams {
            ellipticity,
            beam_direction,
        };

        input.ellipticity = Some(ellipticity_params);
        Ok(())
    }

    /// Parse OPCONS card - optical constants calculation
    fn parse_opcons_card(&self, _card: &Card, input: &mut FeffInput) -> Result<()> {
        // OPCONS card is a flag card without parameters
        // Just mark it as enabled
        let opcons = OpConsParams { enabled: true };

        input.opcons = Some(opcons);
        Ok(())
    }

    /// Parse TDLDA card - time-dependent local density approximation parameters
    fn parse_tdlda_card(&self, card: &Card, input: &mut FeffInput) -> Result<()> {
        // Default parameters
        let mut tdlda = TdldaParams {
            iscreen: 2,  // Default to TDLDA (2)
            icalc: 0,    // Default to SCF+XAS (0)
            elow: -20.0, // Default to -20 eV below edge
            ehigh: 30.0, // Default to 30 eV above edge
            estep: 0.1,  // Default step size 0.1 eV
            gamma: 0.1,  // Default broadening 0.1 eV
        };

        if card.content.is_empty() {
            input.tdlda = Some(tdlda);
            return Ok(());
        }

        let line = card.content[0].trim();
        let fields: Vec<&str> = line.split_whitespace().collect();

        // Format: TDLDA iscreen icalc elow ehigh estep gamma

        // Parse provided parameters
        if !fields.is_empty() {
            tdlda.iscreen = fields[0].parse().unwrap_or(2);
        }

        if fields.len() > 1 {
            tdlda.icalc = fields[1].parse().unwrap_or(0);
        }

        if fields.len() > 2 {
            tdlda.elow = fields[2].parse().unwrap_or(-20.0);
        }

        if fields.len() > 3 {
            tdlda.ehigh = fields[3].parse().unwrap_or(30.0);
        }

        if fields.len() > 4 {
            tdlda.estep = fields[4].parse().unwrap_or(0.1);
        }

        if fields.len() > 5 {
            tdlda.gamma = fields[5].parse().unwrap_or(0.1);
        }

        // Validate parameters
        if tdlda.iscreen < 0 || tdlda.iscreen > 2 {
            return Err(InputError::InvalidFormat(format!(
                "Invalid iscreen value: {}. Must be 0 (RPA), 1 (TDA), or 2 (TDLDA)",
                tdlda.iscreen
            )));
        }

        if tdlda.icalc < 0 || tdlda.icalc > 3 {
            return Err(InputError::InvalidFormat(format!(
                "Invalid icalc value: {}. Must be between 0 and 3",
                tdlda.icalc
            )));
        }

        if tdlda.estep <= 0.0 {
            return Err(InputError::InvalidFormat(format!(
                "Invalid energy step: {}. Must be positive",
                tdlda.estep
            )));
        }

        if tdlda.gamma <= 0.0 {
            return Err(InputError::InvalidFormat(format!(
                "Invalid gamma value: {}. Must be positive",
                tdlda.gamma
            )));
        }

        input.tdlda = Some(tdlda);
        Ok(())
    }

    /// Parse MULTIPOLE card - multipole transitions parameters
    fn parse_multipole_card(&self, card: &Card, input: &mut FeffInput) -> Result<()> {
        // Default parameters
        let mut multipole = MultipoleParams {
            lmax: 3,   // Default maximum l
            morder: 2, // Default to quadrupole
            tensor: 0, // Default tensor off
        };

        if card.content.is_empty() {
            input.multipole = Some(multipole);
            return Ok(());
        }

        let line = card.content[0].trim();
        let fields: Vec<&str> = line.split_whitespace().collect();

        // Format: MULTIPOLE lmax morder tensor

        // Parse provided parameters
        if !fields.is_empty() {
            multipole.lmax = fields[0].parse().unwrap_or(3);
        }

        if fields.len() > 1 {
            multipole.morder = fields[1].parse().unwrap_or(2);
        }

        if fields.len() > 2 {
            multipole.tensor = fields[2].parse().unwrap_or(0);
        }

        // Validate parameters
        if multipole.lmax < 1 {
            return Err(InputError::InvalidFormat(format!(
                "Invalid lmax value: {}. Must be positive",
                multipole.lmax
            )));
        }

        if multipole.morder < 1 || multipole.morder > 4 {
            return Err(InputError::InvalidFormat(format!(
                "Invalid morder value: {}. Must be between 1 and 4",
                multipole.morder
            )));
        }

        if multipole.tensor != 0 && multipole.tensor != 1 {
            return Err(InputError::InvalidFormat(format!(
                "Invalid tensor value: {}. Must be 0 or 1",
                multipole.tensor
            )));
        }

        input.multipole = Some(multipole);
        Ok(())
    }

    /// Parse SCREEN card - self-energy corrections parameters
    fn parse_screen_card(&self, card: &Card, input: &mut FeffInput) -> Result<()> {
        // Default parameters
        let mut screen = ScreenParams {
            iself: 1,   // Default to HL scheme
            iscreen: 1, // Default to screened
            ca1: 1.0,   // Default real part prefactor
            ci1: 1.0,   // Default imaginary part prefactor
        };

        if card.content.is_empty() {
            input.screen = Some(screen);
            return Ok(());
        }

        let line = card.content[0].trim();
        let fields: Vec<&str> = line.split_whitespace().collect();

        // Format: SCREEN iself iscreen ca1 ci1

        // Parse provided parameters
        if !fields.is_empty() {
            screen.iself = fields[0].parse().unwrap_or(1);
        }

        if fields.len() > 1 {
            screen.iscreen = fields[1].parse().unwrap_or(1);
        }

        if fields.len() > 2 {
            screen.ca1 = fields[2].parse().unwrap_or(1.0);
        }

        if fields.len() > 3 {
            screen.ci1 = fields[3].parse().unwrap_or(1.0);
        }

        // Validate parameters
        if screen.iself < 0 || screen.iself > 2 {
            return Err(InputError::InvalidFormat(format!(
                "Invalid iself value: {}. Must be 0 (none), 1 (HL), or 2 (DH)",
                screen.iself
            )));
        }

        if screen.iscreen < 0 || screen.iscreen > 2 {
            return Err(InputError::InvalidFormat(format!(
                "Invalid iscreen value: {}. Must be 0 (unscreened), 1 (screened), or 2 (partial)",
                screen.iscreen
            )));
        }

        input.screen = Some(screen);
        Ok(())
    }

    /// Parse SPECTRAL card - spectral function convolution parameters
    fn parse_spectral_card(&self, card: &Card, input: &mut FeffInput) -> Result<()> {
        // Default parameters
        let mut spectral = SpectralParams {
            ispect: 1,   // Default to enabled
            ispsharp: 0, // Default sharpening off
            isprule: 0,  // Default Fermi level determination method
            emin: -20.0, // Default range below Fermi level
            emax: 20.0,  // Default range above Fermi level
            estep: 0.1,  // Default energy grid spacing
        };

        if card.content.is_empty() {
            input.spectral = Some(spectral);
            return Ok(());
        }

        let line = card.content[0].trim();
        let fields: Vec<&str> = line.split_whitespace().collect();

        // Format: SPECTRAL ispect ispsharp isprule emin emax estep

        // Parse provided parameters
        if !fields.is_empty() {
            spectral.ispect = fields[0].parse().unwrap_or(1);
        }

        if fields.len() > 1 {
            spectral.ispsharp = fields[1].parse().unwrap_or(0);
        }

        if fields.len() > 2 {
            spectral.isprule = fields[2].parse().unwrap_or(0);
        }

        if fields.len() > 3 {
            spectral.emin = fields[3].parse().unwrap_or(-20.0);
        }

        if fields.len() > 4 {
            spectral.emax = fields[4].parse().unwrap_or(20.0);
        }

        if fields.len() > 5 {
            spectral.estep = fields[5].parse().unwrap_or(0.1);
        }

        // Validate parameters
        if spectral.ispect != 0 && spectral.ispect != 1 {
            return Err(InputError::InvalidFormat(format!(
                "Invalid ispect value: {}. Must be 0 (no) or 1 (yes)",
                spectral.ispect
            )));
        }

        if spectral.ispsharp != 0 && spectral.ispsharp != 1 {
            return Err(InputError::InvalidFormat(format!(
                "Invalid ispsharp value: {}. Must be 0 (no) or 1 (yes)",
                spectral.ispsharp
            )));
        }

        if spectral.estep <= 0.0 {
            return Err(InputError::InvalidFormat(format!(
                "Invalid energy step: {}. Must be positive",
                spectral.estep
            )));
        }

        input.spectral = Some(spectral);
        Ok(())
    }

    /// Parse DIMENSIONS card - array dimensioning parameters
    fn parse_dimensions_card(&self, card: &Card, input: &mut FeffInput) -> Result<()> {
        // Default parameters
        let mut dimensions = DimensionsParams {
            nat: 15,    // Default maximum l quantum number
            nph: 200,   // Default maximum atomic sites
            lx: 200,    // Default maximum r-mesh points
            npot: 8,    // Default maximum unique potentials
            nstat: 500, // Default maximum number of paths
        };

        if card.content.is_empty() {
            input.dimensions = Some(dimensions);
            return Ok(());
        }

        let line = card.content[0].trim();
        let fields: Vec<&str> = line.split_whitespace().collect();

        // Format: DIMENSIONS nat nph lx npot nstat

        // Parse provided parameters
        if !fields.is_empty() {
            dimensions.nat = fields[0].parse().unwrap_or(15);
        }

        if fields.len() > 1 {
            dimensions.nph = fields[1].parse().unwrap_or(200);
        }

        if fields.len() > 2 {
            dimensions.lx = fields[2].parse().unwrap_or(200);
        }

        if fields.len() > 3 {
            dimensions.npot = fields[3].parse().unwrap_or(8);
        }

        if fields.len() > 4 {
            dimensions.nstat = fields[4].parse().unwrap_or(500);
        }

        // Validate parameters
        if dimensions.nat <= 0 {
            return Err(InputError::InvalidFormat(format!(
                "Invalid nat value: {}. Must be positive",
                dimensions.nat
            )));
        }

        if dimensions.nph <= 0 {
            return Err(InputError::InvalidFormat(format!(
                "Invalid nph value: {}. Must be positive",
                dimensions.nph
            )));
        }

        if dimensions.lx <= 0 {
            return Err(InputError::InvalidFormat(format!(
                "Invalid lx value: {}. Must be positive",
                dimensions.lx
            )));
        }

        if dimensions.npot <= 0 {
            return Err(InputError::InvalidFormat(format!(
                "Invalid npot value: {}. Must be positive",
                dimensions.npot
            )));
        }

        if dimensions.nstat <= 0 {
            return Err(InputError::InvalidFormat(format!(
                "Invalid nstat value: {}. Must be positive",
                dimensions.nstat
            )));
        }

        input.dimensions = Some(dimensions);
        Ok(())
    }

    /// Parse RDINP card - read input from a different file
    fn parse_rdinp_card(&self, card: &Card, input: &mut FeffInput) -> Result<()> {
        if card.content.is_empty() {
            return Err(InputError::ParseError(
                "RDINP card requires a file name".to_string(),
            ));
        }

        let line = card.content[0].trim();
        let fields: Vec<&str> = line.split_whitespace().collect();

        // Format: RDINP filename
        if fields.is_empty() {
            return Err(InputError::ParseError(
                "RDINP card requires a file name".to_string(),
            ));
        }

        // Get the filename
        let file_name = fields[0].to_string();

        // Basic validation
        if file_name.is_empty() {
            return Err(InputError::InvalidFormat(
                "RDINP filename cannot be empty".to_string(),
            ));
        }

        let rdinp = RdinpParams { file_name };
        input.rdinp = Some(rdinp);
        Ok(())
    }

    /// Parse BANDSTRUCTURE card - band structure calculation parameters
    fn parse_bandstructure_card(&self, card: &Card, input: &mut FeffInput) -> Result<()> {
        // Default parameters
        let mut bandstructure = BandstructureParams {
            nk: 100,     // Default number of k-points
            emin: -10.0, // Default minimum energy
            emax: 10.0,  // Default maximum energy
            estep: 0.1,  // Default energy step
            kmesh: 1,    // Default to uniform mesh
            symmetry: 1, // Default to use symmetry
        };

        // If card has no content, use defaults
        if card.content.is_empty() {
            input.bandstructure = Some(bandstructure);
            return Ok(());
        }

        let line = card.content[0].trim();
        let fields: Vec<&str> = line.split_whitespace().collect();

        // Format: BANDSTRUCTURE nk emin emax estep kmesh symmetry

        // Parse provided parameters
        if !fields.is_empty() {
            bandstructure.nk = fields[0].parse().unwrap_or(100);
        }

        if fields.len() > 1 {
            bandstructure.emin = fields[1].parse().unwrap_or(-10.0);
        }

        if fields.len() > 2 {
            bandstructure.emax = fields[2].parse().unwrap_or(10.0);
        }

        if fields.len() > 3 {
            bandstructure.estep = fields[3].parse().unwrap_or(0.1);
        }

        if fields.len() > 4 {
            bandstructure.kmesh = fields[4].parse().unwrap_or(1);
        }

        if fields.len() > 5 {
            bandstructure.symmetry = fields[5].parse().unwrap_or(1);
        }

        // Validate parameters
        if bandstructure.nk <= 0 {
            return Err(InputError::InvalidFormat(format!(
                "Invalid number of k-points: {}. Must be positive",
                bandstructure.nk
            )));
        }

        if bandstructure.emax <= bandstructure.emin {
            return Err(InputError::InvalidFormat(format!(
                "Maximum energy ({}) must be greater than minimum energy ({})",
                bandstructure.emax, bandstructure.emin
            )));
        }

        if bandstructure.estep <= 0.0 {
            return Err(InputError::InvalidFormat(format!(
                "Energy step ({}) must be positive",
                bandstructure.estep
            )));
        }

        if bandstructure.kmesh != 0 && bandstructure.kmesh != 1 {
            return Err(InputError::InvalidFormat(format!(
                "Invalid kmesh flag: {}. Must be 0 (user-defined) or 1 (uniform)",
                bandstructure.kmesh
            )));
        }

        if bandstructure.symmetry != 0 && bandstructure.symmetry != 1 {
            return Err(InputError::InvalidFormat(format!(
                "Invalid symmetry flag: {}. Must be 0 (no symmetry) or 1 (use symmetry)",
                bandstructure.symmetry
            )));
        }

        input.bandstructure = Some(bandstructure);
        Ok(())
    }

    /// Parse KMESH card - k-space mesh for band structure calculations
    fn parse_kmesh_card(&self, card: &Card, input: &mut FeffInput) -> Result<()> {
        // Default parameters
        let mut kmesh = KmeshParams {
            nx: 10,              // Default x-direction points
            ny: 10,              // Default y-direction points
            nz: 10,              // Default z-direction points
            kpoints: Vec::new(), // Empty list of explicit k-points
        };

        // KMESH card requires at least one line
        if card.content.is_empty() {
            return Err(InputError::ParseError(
                "KMESH card requires parameters".to_string(),
            ));
        }

        let first_line = card.content[0].trim();
        let first_fields: Vec<&str> = first_line.split_whitespace().collect();

        // Format (first line): KMESH nx ny nz
        if first_fields.len() < 3 {
            return Err(InputError::ParseError(
                "KMESH card requires at least nx, ny, nz parameters".to_string(),
            ));
        }

        // Parse mesh dimensions
        kmesh.nx = first_fields[0].parse().map_err(|_| {
            InputError::ParseError(format!("Invalid nx value: {}", first_fields[0]))
        })?;

        kmesh.ny = first_fields[1].parse().map_err(|_| {
            InputError::ParseError(format!("Invalid ny value: {}", first_fields[1]))
        })?;

        kmesh.nz = first_fields[2].parse().map_err(|_| {
            InputError::ParseError(format!("Invalid nz value: {}", first_fields[2]))
        })?;

        // Validate mesh dimensions
        if kmesh.nx <= 0 || kmesh.ny <= 0 || kmesh.nz <= 0 {
            return Err(InputError::InvalidFormat(
                "Mesh dimensions must be positive".to_string(),
            ));
        }

        // If there are more lines, parse explicit k-points
        for i in 1..card.content.len() {
            let line = card.content[i].trim();

            // Skip comment lines and empty lines
            if line.starts_with('*') || line.starts_with('#') || line.is_empty() {
                continue;
            }

            let fields: Vec<&str> = line.split_whitespace().collect();

            // Each k-point needs 3 coordinates
            if fields.len() < 3 {
                return Err(InputError::ParseError(format!(
                    "K-point at line {} needs 3 coordinates",
                    card.line_number + i
                )));
            }

            // Parse coordinates
            let kx = fields[0]
                .parse::<f64>()
                .map_err(|_| InputError::ParseError(format!("Invalid kx value: {}", fields[0])))?;

            let ky = fields[1]
                .parse::<f64>()
                .map_err(|_| InputError::ParseError(format!("Invalid ky value: {}", fields[1])))?;

            let kz = fields[2]
                .parse::<f64>()
                .map_err(|_| InputError::ParseError(format!("Invalid kz value: {}", fields[2])))?;

            // Add to list of k-points
            kmesh.kpoints.push((kx, ky, kz));
        }

        input.kmesh = Some(kmesh);
        Ok(())
    }

    /// Parse RESTART card - load saved calculation results
    fn parse_restart_card(&self, card: &Card, input: &mut FeffInput) -> Result<()> {
        if card.content.is_empty() {
            return Err(InputError::ParseError(
                "RESTART card requires a module name".to_string(),
            ));
        }

        let line = card.content[0].trim();
        let fields: Vec<&str> = line.split_whitespace().collect();

        // Format: RESTART module [file_name]
        if fields.is_empty() {
            return Err(InputError::ParseError(
                "RESTART card requires a module name".to_string(),
            ));
        }

        // Get the module name
        let module = fields[0].to_string();

        // Basic validation
        if module.is_empty() {
            return Err(InputError::InvalidFormat(
                "RESTART module name cannot be empty".to_string(),
            ));
        }

        // Get optional file name
        let file_name = if fields.len() > 1 {
            Some(fields[1].to_string())
        } else {
            None
        };

        let restart = RestartParams { module, file_name };
        input.restart = Some(restart);
        Ok(())
    }

    /// Parse DOS card - density of states calculation parameters
    fn parse_dos_card(&self, card: &Card, input: &mut FeffInput) -> Result<()> {
        // Default parameters
        let mut dos = DosParams {
            emin: -20.0,        // Default minimum energy
            emax: 20.0,         // Default maximum energy
            estep: 0.1,         // Default energy step
            gamma: 0.2,         // Default broadening
            params: Vec::new(), // Default empty additional parameters
        };

        if card.content.is_empty() {
            input.dos = Some(dos);
            return Ok(());
        }

        let line = card.content[0].trim();
        let fields: Vec<&str> = line.split_whitespace().collect();

        // Format: DOS emin emax estep gamma [additional_params...]

        // Parse provided parameters
        if !fields.is_empty() {
            dos.emin = fields[0].parse().unwrap_or(-20.0);
        }

        if fields.len() > 1 {
            dos.emax = fields[1].parse().unwrap_or(20.0);
        }

        if fields.len() > 2 {
            dos.estep = fields[2].parse().unwrap_or(0.1);
        }

        if fields.len() > 3 {
            dos.gamma = fields[3].parse().unwrap_or(0.2);
        }

        // Parse any additional parameters
        for field in fields.iter().skip(4) {
            if let Ok(param) = field.parse::<f64>() {
                dos.params.push(param);
            }
        }

        // Validate parameters
        if dos.emax <= dos.emin {
            return Err(InputError::InvalidFormat(format!(
                "Maximum energy ({}) must be greater than minimum energy ({})",
                dos.emax, dos.emin
            )));
        }

        if dos.estep <= 0.0 {
            return Err(InputError::InvalidFormat(format!(
                "Energy step ({}) must be positive",
                dos.estep
            )));
        }

        if dos.gamma < 0.0 {
            return Err(InputError::InvalidFormat(format!(
                "Broadening parameter ({}) cannot be negative",
                dos.gamma
            )));
        }

        input.dos = Some(dos);
        Ok(())
    }

    /// Parse CIFS card - crystallographic information file
    fn parse_cifs_card(&self, card: &Card, input: &mut FeffInput) -> Result<()> {
        if card.content.is_empty() {
            return Err(InputError::ParseError(
                "CIFS card requires a file name".to_string(),
            ));
        }

        let line = card.content[0].trim();
        let fields: Vec<&str> = line.split_whitespace().collect();

        // Format: CIFS file_name [site_index [distance_cutoff]]
        if fields.is_empty() {
            return Err(InputError::ParseError(
                "CIFS card requires a file name".to_string(),
            ));
        }

        // Get the filename
        let file_name = fields[0].to_string();

        // Basic validation
        if file_name.is_empty() {
            return Err(InputError::InvalidFormat(
                "CIFS filename cannot be empty".to_string(),
            ));
        }

        // Parse optional site index
        let site_index = if fields.len() > 1 {
            match fields[1].parse::<i32>() {
                Ok(index) => Some(index),
                Err(_) => {
                    return Err(InputError::ParseError(format!(
                        "Invalid site index: {}",
                        fields[1]
                    )))
                }
            }
        } else {
            None
        };

        // Parse optional distance cutoff
        let distance_cutoff = if fields.len() > 2 {
            match fields[2].parse::<f64>() {
                Ok(cutoff) => {
                    if cutoff <= 0.0 {
                        return Err(InputError::InvalidFormat(format!(
                            "Distance cutoff must be positive, found: {}",
                            cutoff
                        )));
                    }
                    Some(cutoff)
                }
                Err(_) => {
                    return Err(InputError::ParseError(format!(
                        "Invalid distance cutoff: {}",
                        fields[2]
                    )))
                }
            }
        } else {
            None
        };

        let cifs = CifsParams {
            file_name,
            site_index,
            distance_cutoff,
        };

        input.cifs = Some(cifs);
        Ok(())
    }
}
