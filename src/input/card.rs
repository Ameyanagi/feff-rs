/*
MIT License with FEFF10 Attribution

Copyright (c) 2025 Ameyanagi

Based on or developed using Distribution: FEFF10.0
Copyright (c) 2020 FEFF Project, University of Washington and SLAC National Accelerator Laboratory.
All rights reserved.
*/

//! FEFF input card representation

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

/// Helper function to check if a string is a valid FEFF card name
pub fn is_card_name(s: &str) -> bool {
    // Get first word from string
    let first_word = s.split_whitespace().next().unwrap_or("");

    // Card names are typically all uppercase and 2+ characters
    // Special handling for S02 which is 3 characters but starts with 'S'
    !first_word.is_empty()
        && ((first_word.chars().all(|c| c.is_ascii_uppercase()) && first_word.len() >= 2)
            || first_word == "S02")
}
