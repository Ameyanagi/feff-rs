use crate::atoms::{database, Atom, AtomicStructure, Vector3D};
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::{Path, PathBuf};
use thiserror::Error;

/// Errors that can occur during FEFF input parsing
#[derive(Error, Debug)]
pub enum InputError {
    #[error("IO error: {0}")]
    IoError(#[from] io::Error),

    #[error("Parse error: {0}")]
    ParseError(String),

    #[error("Invalid input format: {0}")]
    InvalidFormat(String),

    #[error("Missing required card: {0}")]
    MissingCard(String),

    #[error("Invalid potential: {0}")]
    InvalidPotential(String),

    #[error("Invalid atomic structure: {0}")]
    InvalidStructure(String),
}

/// Coordinate system enumeration for atomic positions
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum CoordinateSystem {
    /// Cartesian (x, y, z)
    Cartesian,
    /// Spherical (r, theta, phi)
    Spherical,
    /// Cylindrical (rho, phi, z)
    Cylindrical,
}

/// Represents a FEFF input card with its content
#[derive(Debug, Clone)]
pub struct Card {
    /// Name of the card (e.g., "ATOMS", "POTENTIALS")
    pub name: String,
    /// Content lines of the card
    pub content: Vec<String>,
    /// Line number where the card starts in the input file
    pub line_number: usize,
}

/// FEFF input parser configuration
#[derive(Debug, Clone)]
pub struct ParserConfig {
    /// Path to the input file
    pub input_path: PathBuf,
    /// Whether to validate atomic structure against FEFF requirements
    pub validate: bool,
    /// Whether to add hydrogen to under-coordinated atoms
    pub add_hydrogens: bool,
    /// Whether to enable debugging output during parsing
    pub debug: bool,
}

impl Default for ParserConfig {
    fn default() -> Self {
        Self {
            input_path: PathBuf::from("feff.inp"),
            validate: true,
            add_hydrogens: false,
            debug: false,
        }
    }
}

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
    pub fn parse<P: AsRef<Path>>(&mut self, path: Option<P>) -> Result<FeffInput, InputError> {
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
    fn parse_cards<R: BufRead>(&mut self, reader: R) -> Result<(), InputError> {
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
    fn parse_atoms_card(&self, card: &Card, input: &mut FeffInput) -> Result<(), InputError> {
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
    fn parse_potentials_card(&self, card: &Card, input: &mut FeffInput) -> Result<(), InputError> {
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
    fn parse_control_card(&self, card: &Card, input: &mut FeffInput) -> Result<(), InputError> {
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
    fn parse_exchange_card(&self, card: &Card, input: &mut FeffInput) -> Result<(), InputError> {
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
    fn parse_hole_card(&self, card: &Card, input: &mut FeffInput) -> Result<(), InputError> {
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
    fn parse_scf_card(&self, card: &Card, input: &mut FeffInput) -> Result<(), InputError> {
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
    fn parse_fms_card(&self, card: &Card, input: &mut FeffInput) -> Result<(), InputError> {
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
    fn parse_xanes_card(&self, card: &Card, input: &mut FeffInput) -> Result<(), InputError> {
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
    fn parse_title_card(&self, card: &Card, input: &mut FeffInput) -> Result<(), InputError> {
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
    fn parse_rpath_card(&self, card: &Card, input: &mut FeffInput) -> Result<(), InputError> {
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
    fn parse_print_card(&self, card: &Card, input: &mut FeffInput) -> Result<(), InputError> {
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
    fn parse_corrections_card(&self, card: &Card, input: &mut FeffInput) -> Result<(), InputError> {
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
    fn parse_s02_card(&self, card: &Card, input: &mut FeffInput) -> Result<(), InputError> {
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
    fn parse_edge_card(&self, card: &Card, input: &mut FeffInput) -> Result<(), InputError> {
        if card.content.is_empty() {
            return Ok(());
        }

        let line = card.content[0].trim();
        let fields: Vec<&str> = line.split_whitespace().collect();

        // Format: EDGE K [energy]
        if !fields.is_empty() {
            let edge = EdgeParams {
                edge_type: fields[0].to_string(),
                energy: if fields.len() > 1 {
                    fields[1].parse::<f64>().ok()
                } else {
                    None
                },
            };

            input.edge = Some(edge);
        }

        Ok(())
    }

    /// Parse DEBYE card - thermal effects parameters
    fn parse_debye_card(&self, card: &Card, input: &mut FeffInput) -> Result<(), InputError> {
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
    fn parse_ldos_card(&self, card: &Card, input: &mut FeffInput) -> Result<(), InputError> {
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
    fn parse_exafs_card(&self, card: &Card, input: &mut FeffInput) -> Result<(), InputError> {
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
    fn parse_danes_card(&self, card: &Card, input: &mut FeffInput) -> Result<(), InputError> {
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
    fn parse_corehole_card(&self, card: &Card, input: &mut FeffInput) -> Result<(), InputError> {
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
    fn parse_polarization_card(
        &self,
        card: &Card,
        input: &mut FeffInput,
    ) -> Result<(), InputError> {
        if card.content.is_empty() {
            return Ok(());
        }

        let line = card.content[0].trim();
        let fields: Vec<&str> = line.split_whitespace().collect();

        // Format: POLARIZATION x y z [ellipticity]
        if fields.len() >= 3 {
            let polarization = PolarizationParams {
                x: fields[0].parse().unwrap_or(1.0),
                y: fields[1].parse().unwrap_or(0.0),
                z: fields[2].parse().unwrap_or(0.0),
                ellipticity: if fields.len() > 3 {
                    fields[3].parse::<f64>().ok()
                } else {
                    None
                },
            };

            input.polarization = Some(polarization);
        }

        Ok(())
    }

    /// Parse REAL card - real space grid parameters
    fn parse_real_card(&self, card: &Card, input: &mut FeffInput) -> Result<(), InputError> {
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
    fn parse_reciprocal_card(&self, card: &Card, input: &mut FeffInput) -> Result<(), InputError> {
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
}

/// Helper function to check if a string is a valid FEFF card name
fn is_card_name(s: &str) -> bool {
    // Get first word from string
    let first_word = s.split_whitespace().next().unwrap_or("");

    // Card names are typically all uppercase and 2+ characters
    // Special handling for S02 which is 3 characters but starts with 'S'
    !first_word.is_empty()
        && ((first_word.chars().all(|c| c.is_ascii_uppercase()) && first_word.len() >= 2)
            || first_word == "S02")
}

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

    /// Unknown cards
    pub unknown_cards: Vec<Card>,
}

impl FeffInput {
    /// Create a new empty FEFF input
    pub fn new() -> Self {
        Self::default()
    }

    /// Validate the input
    pub fn validate(&self) -> Result<(), InputError> {
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
    pub fn write<P: AsRef<Path>>(&self, path: P) -> Result<(), InputError> {
        use std::fs::File;
        use std::io::Write;

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

/// Public functions
///
/// Parse a FEFF input file with default configuration
pub fn parse_feff_input<P: AsRef<Path>>(path: P) -> Result<FeffInput, InputError> {
    let config = ParserConfig {
        input_path: path.as_ref().to_path_buf(),
        ..Default::default()
    };

    let mut parser = FeffInputParser::new(config);
    parser.parse::<&Path>(None)
}

/// Create a new empty FEFF input
pub fn new_feff_input() -> FeffInput {
    FeffInput::new()
}
