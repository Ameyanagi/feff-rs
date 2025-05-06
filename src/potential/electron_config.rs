/*
MIT License with FEFF10 Attribution

Copyright (c) 2025 Ameyanagi

Based on or developed using Distribution: FEFF10.0
Copyright (c) 2020 FEFF Project, University of Washington and SLAC National Accelerator Laboratory.
All rights reserved.
*/

//! Electron configuration utilities for atomic calculations

/// Shell type in spectroscopic notation
pub enum ShellType {
    S, // l=0
    P, // l=1
    D, // l=2
    F, // l=3
}

/// Convert shell type to angular momentum quantum number
pub fn shell_to_l(shell: ShellType) -> i32 {
    match shell {
        ShellType::S => 0,
        ShellType::P => 1,
        ShellType::D => 2,
        ShellType::F => 3,
    }
}

/// Determine electronic shells for a given element
///
/// # Arguments
///
/// * `atomic_number` - The atomic number Z
///
/// # Returns
///
/// A vector of (n, l) pairs for each shell, ordered by energy (most negative first)
pub fn determine_shells_for_element(atomic_number: i32) -> Vec<(i32, i32)> {
    let mut shells = Vec::new();

    // Core electrons (1s, 2s, 2p) for all elements
    shells.push((1, 0)); // 1s

    if atomic_number > 2 {
        shells.push((2, 0)); // 2s
        shells.push((2, 1)); // 2p
    }

    // 3s, 3p, 3d
    if atomic_number > 10 {
        shells.push((3, 0)); // 3s
        shells.push((3, 1)); // 3p
    }

    if atomic_number > 18 {
        shells.push((3, 2)); // 3d
    }

    // 4s, 4p, 4d, 4f
    if atomic_number > 20 {
        shells.push((4, 0)); // 4s
    }

    if atomic_number > 30 {
        shells.push((4, 1)); // 4p
    }

    if atomic_number > 36 {
        shells.push((4, 2)); // 4d
    }

    if atomic_number > 56 {
        shells.push((4, 3)); // 4f
    }

    // 5s, 5p, 5d
    if atomic_number > 38 {
        shells.push((5, 0)); // 5s
    }

    if atomic_number > 48 {
        shells.push((5, 1)); // 5p
    }

    if atomic_number > 56 {
        shells.push((5, 2)); // 5d
    }

    // 6s, 6p, 6d
    if atomic_number > 56 {
        shells.push((6, 0)); // 6s
    }

    if atomic_number > 80 {
        shells.push((6, 1)); // 6p
    }

    if atomic_number > 88 {
        shells.push((6, 2)); // 6d
    }

    // 7s, 7p
    if atomic_number > 88 {
        shells.push((7, 0)); // 7s
    }

    if atomic_number > 103 {
        shells.push((7, 1)); // 7p
    }

    shells
}

/// Get the electron occupation of a shell in an atom
///
/// # Arguments
///
/// * `atomic_number` - The atomic number Z
/// * `n` - Principal quantum number
/// * `l` - Angular momentum quantum number
///
/// # Returns
///
/// The number of electrons in the shell
pub fn shell_occupation(atomic_number: i32, n: i32, l: i32) -> i32 {
    // Maximum occupation of a shell: 2*(2l+1)
    let max_occupancy = 2 * (2 * l + 1);

    // Use simplified model based on blocks in periodic table
    match (n, l) {
        (1, 0) => i32::min(2, atomic_number), // 1s
        (2, 0) => {
            if atomic_number > 2 {
                i32::min(2, atomic_number - 2)
            } else {
                0
            }
        } // 2s
        (2, 1) => {
            if atomic_number > 4 {
                i32::min(6, atomic_number - 4)
            } else {
                0
            }
        } // 2p
        (3, 0) => {
            if atomic_number > 10 {
                i32::min(2, atomic_number - 10)
            } else {
                0
            }
        } // 3s
        (3, 1) => {
            if atomic_number > 12 {
                i32::min(6, atomic_number - 12)
            } else {
                0
            }
        } // 3p
        // Other shells follow different filling order due to aufbau principle
        // This is a simplified version - real atoms have more complex electron configurations
        _ => {
            // Default to full shell or zero electrons
            if atomic_number > 18 {
                max_occupancy
            } else {
                0
            }
        }
    }
}

/// Calculate the core-hole occupation for a given edge
///
/// # Arguments
///
/// * `edge` - The absorption edge (K, L1, L2, L3, etc.)
/// * `n` - Principal quantum number
/// * `l` - Angular momentum quantum number
///
/// # Returns
///
/// The core-hole occupation (usually one less than the normal occupation)
pub fn core_hole_occupation(edge: &str, n: i32, l: i32) -> i32 {
    match edge {
        "K" => {
            if n == 1 && l == 0 {
                // K-edge: 1s core hole
                1 // instead of 2
            } else {
                2 * (2 * l + 1) // full occupation for all other shells
            }
        }
        "L1" => {
            if n == 2 && l == 0 {
                // L1-edge: 2s core hole
                1 // instead of 2
            } else {
                2 * (2 * l + 1) // full occupation for all other shells
            }
        }
        "L2" | "L3" => {
            if n == 2 && l == 1 {
                // L2,L3-edge: 2p core hole
                5 // instead of 6
            } else {
                2 * (2 * l + 1) // full occupation for all other shells
            }
        }
        // Add other edges as needed
        _ => 2 * (2 * l + 1), // Default to full shell
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_determine_shells() {
        // Test H (Z=1)
        let h_shells = determine_shells_for_element(1);
        assert_eq!(h_shells.len(), 1);
        assert_eq!(h_shells[0], (1, 0)); // 1s

        // Test O (Z=8)
        let o_shells = determine_shells_for_element(8);
        assert_eq!(o_shells.len(), 3);
        assert_eq!(o_shells[0], (1, 0)); // 1s
        assert_eq!(o_shells[1], (2, 0)); // 2s
        assert_eq!(o_shells[2], (2, 1)); // 2p

        // Test Fe (Z=26)
        let fe_shells = determine_shells_for_element(26);
        assert!(fe_shells.len() >= 6); // At least these shells
        assert!(fe_shells.contains(&(1, 0))); // 1s
        assert!(fe_shells.contains(&(2, 0))); // 2s
        assert!(fe_shells.contains(&(2, 1))); // 2p
        assert!(fe_shells.contains(&(3, 0))); // 3s
        assert!(fe_shells.contains(&(3, 1))); // 3p
        assert!(fe_shells.contains(&(3, 2))); // 3d
    }

    #[test]
    fn test_shell_occupation() {
        // Test H (Z=1) - should have 1 electron in 1s
        assert_eq!(shell_occupation(1, 1, 0), 1);

        // Test He (Z=2) - should have 2 electrons in 1s
        assert_eq!(shell_occupation(2, 1, 0), 2);

        // Test Li (Z=3) - should have 2 electrons in 1s, 1 in 2s
        assert_eq!(shell_occupation(3, 1, 0), 2);
        assert_eq!(shell_occupation(3, 2, 0), 1);

        // Test O (Z=8) - should have 2,2,4 electrons in 1s,2s,2p
        assert_eq!(shell_occupation(8, 1, 0), 2);
        assert_eq!(shell_occupation(8, 2, 0), 2);
        assert_eq!(shell_occupation(8, 2, 1), 4);
    }

    #[test]
    fn test_core_hole_occupation() {
        // Test K-edge core hole (1s)
        assert_eq!(core_hole_occupation("K", 1, 0), 1); // reduced by 1
        assert_eq!(core_hole_occupation("K", 2, 0), 2); // unchanged

        // Test L1-edge core hole (2s)
        assert_eq!(core_hole_occupation("L1", 1, 0), 2); // unchanged
        assert_eq!(core_hole_occupation("L1", 2, 0), 1); // reduced by 1

        // Test L2,L3-edge core hole (2p)
        assert_eq!(core_hole_occupation("L2", 2, 1), 5); // reduced by 1
        assert_eq!(core_hole_occupation("L3", 2, 1), 5); // reduced by 1
    }
}
