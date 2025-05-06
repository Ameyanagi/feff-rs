/*
MIT License with FEFF10 Attribution

Copyright (c) 2025 Ameyanagi

Based on or developed using Distribution: FEFF10.0
Copyright (c) 2020 FEFF Project, University of Washington and SLAC National Accelerator Laboratory.
All rights reserved.
*/

//! Path visualization utilities
//!
//! This module provides functions for visualizing scattering paths
//! in various formats for debugging and analysis.

use crate::atoms::structure::AtomicStructure;
use crate::path::path::Path;

/// Formats a path as a human-readable string
///
/// # Arguments
///
/// * `path` - The path to format
/// * `structure` - The atomic structure containing the atoms in the path
///
/// # Returns
///
/// A formatted string describing the path
pub fn format_path_description(path: &Path, structure: &AtomicStructure) -> String {
    let mut description = String::new();

    // Add path type and length
    description.push_str(&format!(
        "{:?} path, length: {:.3} Å, degeneracy: {}, importance: {:.6}\n",
        path.path_type, path.total_length, path.degeneracy, path.importance
    ));

    // Add atom sequence with element symbols
    description.push_str("Atom sequence: ");
    for (i, &atom_idx) in path.atom_sequence.iter().enumerate() {
        let atom = structure.atom(atom_idx).unwrap();

        // Format differently for central atom
        if i == 0 || i == path.atom_sequence.len() - 1 {
            description.push_str(&format!("*{}*", atom.symbol()));
        } else {
            description.push_str(atom.symbol());
        }

        // Add arrow or space
        if i < path.atom_sequence.len() - 1 {
            description.push_str(" → ");
        }
    }

    // Add leg details
    description.push_str("\n\nLegs:\n");
    for (i, leg) in path.legs.iter().enumerate() {
        let from_atom = structure.atom(leg.from_atom).unwrap();
        let to_atom = structure.atom(leg.to_atom).unwrap();

        description.push_str(&format!(
            "  {}. {} → {}: {:.3} Å\n",
            i + 1,
            from_atom.symbol(),
            to_atom.symbol(),
            leg.length
        ));
    }

    description
}

/// Creates a table of path summaries for multiple paths
///
/// # Arguments
///
/// * `paths` - The paths to include in the table
/// * `structure` - The atomic structure containing the atoms in the paths
///
/// # Returns
///
/// A string containing a formatted table of path summaries
pub fn create_path_summary_table(paths: &[Path], structure: &AtomicStructure) -> String {
    let mut table = String::from(
        "| # | Type | Length (Å) | Degeneracy | Importance | Path Sequence |\n\
         |---|------|------------|------------|------------|---------------|\n",
    );

    for (i, path) in paths.iter().enumerate() {
        // Create a simplified path sequence
        let path_sequence = path
            .atom_sequence
            .iter()
            .map(|&idx| structure.atom(idx).unwrap().symbol())
            .collect::<Vec<_>>()
            .join(" → ");

        // Add a row to the table
        table.push_str(&format!(
            "| {} | {:?} | {:.3} | {} | {:.6} | {} |\n",
            i + 1,
            path.path_type,
            path.total_length,
            path.degeneracy,
            path.importance,
            path_sequence
        ));
    }

    table
}

/// Generates XYZ format data for visualizing a path with markers
///
/// The XYZ format is a simple text format that can be read by most
/// molecular visualization software to display atomic structures.
/// This function creates an XYZ representation with the path atoms
/// numbered to show the sequence.
///
/// # Arguments
///
/// * `path` - The path to visualize
/// * `structure` - The atomic structure containing the atoms
///
/// # Returns
///
/// A string in XYZ format representing the path
pub fn generate_path_xyz(path: &Path, structure: &AtomicStructure) -> String {
    let mut xyz = String::new();

    // Count the number of atoms in the path
    let unique_atoms: std::collections::HashSet<usize> =
        path.atom_sequence.iter().copied().collect();
    let atom_count = unique_atoms.len();

    // XYZ header (number of atoms and comment line)
    xyz.push_str(&format!("{}\n", atom_count));
    xyz.push_str(&format!(
        "Path visualization: {:?} path, length: {:.3} Å\n",
        path.path_type, path.total_length
    ));

    // Add each atom in the path
    let mut visited = std::collections::HashSet::new();
    for (i, &atom_idx) in path.atom_sequence.iter().enumerate() {
        // Skip if we've already added this atom (but with a different position in the sequence)
        if visited.contains(&atom_idx) {
            continue;
        }

        let atom = structure.atom(atom_idx).unwrap();
        let position = atom.position();

        // Add sequence number to element symbol for visualization
        let element = if i == 0 {
            // Mark the absorber atom
            format!("{}*", atom.symbol())
        } else {
            format!("{}{}", atom.symbol(), i)
        };

        // Add atom line to XYZ file (element, x, y, z)
        xyz.push_str(&format!(
            "{:<4} {:12.6} {:12.6} {:12.6}\n",
            element, position.x, position.y, position.z
        ));

        visited.insert(atom_idx);
    }

    xyz
}

/// Generates a JSON representation of a path for web visualization
///
/// # Arguments
///
/// * `path` - The path to convert to JSON
/// * `structure` - The atomic structure containing the atoms
///
/// # Returns
///
/// A JSON string representing the path
pub fn generate_path_json(path: &Path, structure: &AtomicStructure) -> String {
    let mut json = String::from("{\n");

    // Add path metadata
    json.push_str(&format!("  \"type\": \"{:?}\",\n", path.path_type));
    json.push_str(&format!("  \"totalLength\": {:.6},\n", path.total_length));
    json.push_str(&format!("  \"degeneracy\": {},\n", path.degeneracy));
    json.push_str(&format!("  \"importance\": {:.6},\n", path.importance));

    // Add atom sequence
    json.push_str("  \"atomSequence\": [\n");
    for (i, &atom_idx) in path.atom_sequence.iter().enumerate() {
        let atom = structure.atom(atom_idx).unwrap();

        json.push_str(&format!(
            "    {{\"index\": {}, \"element\": \"{}\", \"position\": [{:.6}, {:.6}, {:.6}]}}",
            atom_idx,
            atom.symbol(),
            atom.position().x,
            atom.position().y,
            atom.position().z
        ));

        if i < path.atom_sequence.len() - 1 {
            json.push_str(",\n");
        } else {
            json.push_str("\n");
        }
    }
    json.push_str("  ],\n");

    // Add legs
    json.push_str("  \"legs\": [\n");
    for (i, leg) in path.legs.iter().enumerate() {
        json.push_str(&format!(
            "    {{\"from\": {}, \"to\": {}, \"length\": {:.6}}}",
            leg.from_atom, leg.to_atom, leg.length
        ));

        if i < path.legs.len() - 1 {
            json.push_str(",\n");
        } else {
            json.push_str("\n");
        }
    }
    json.push_str("  ]\n");

    json.push_str("}\n");

    json
}

/// Exports multiple paths to a comprehensive visualization file
///
/// # Arguments
///
/// * `paths` - The paths to export
/// * `structure` - The atomic structure containing the atoms
/// * `format` - The output format ("xyz", "json", or "text")
///
/// # Returns
///
/// A string containing the requested visualization format
pub fn export_paths(paths: &[Path], structure: &AtomicStructure, format: &str) -> String {
    match format.to_lowercase().as_str() {
        "xyz" => {
            // For XYZ, we'll create a multiframe XYZ file with one frame per path
            let mut xyz = String::new();

            for path in paths {
                // Generate XYZ for this path and append it
                xyz.push_str(&generate_path_xyz(path, structure));
                // Add a blank line between frames
                xyz.push('\n');
            }

            xyz
        }
        "json" => {
            // For JSON, we'll create an array of path objects
            let mut json = String::from("[\n");

            for (i, path) in paths.iter().enumerate() {
                // Generate JSON for this path (without the outer braces)
                let path_json = generate_path_json(path, structure);

                // Add the path JSON to the array
                json.push_str(&path_json);

                if i < paths.len() - 1 {
                    json.push_str(",\n");
                } else {
                    json.push_str("\n");
                }
            }

            json.push_str("]\n");
            json
        }
        _ => {
            // Default to text format with a summary table
            let mut text = String::from("Path Summary\n============\n\n");
            text.push_str(&create_path_summary_table(paths, structure));

            // Add detailed descriptions for each path
            text.push_str("\n\nDetailed Path Descriptions\n========================\n\n");
            for (i, path) in paths.iter().enumerate() {
                text.push_str(&format!("Path #{}\n", i + 1));
                text.push_str(&format_path_description(path, structure));
                text.push_str("\n---\n\n");
            }

            text
        }
    }
}
