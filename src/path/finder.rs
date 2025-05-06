/*
MIT License with FEFF10 Attribution

Copyright (c) 2025 Ameyanagi

Based on or developed using Distribution: FEFF10.0
Copyright (c) 2020 FEFF Project, University of Washington and SLAC National Accelerator Laboratory.
All rights reserved.
*/

//! Path finding and enumeration algorithms
//!
//! This module implements algorithms for finding and enumerating scattering paths
//! for XAFS calculations. It includes functionality for building paths from atomic
//! structures and filtering them based on various criteria.

use std::cmp::Ordering;
use std::collections::{BinaryHeap, HashSet};

use crate::atoms::structure::AtomicStructure;
use crate::path::path::{Path, PathLeg, PathType};

/// Represents a path finding parameter configuration
///
/// This configuration controls the path finding process including
/// the maximum path length, the number of paths to keep, and
/// the various criteria for path filtering.
#[derive(Debug, Clone)]
pub struct PathFinderConfig {
    /// Maximum path length to consider (in Å)
    pub max_path_length: f64,

    /// Maximum number of paths to keep
    pub max_paths: usize,

    /// Maximum number of legs in a path
    pub max_legs: usize,

    /// Filter paths with importance below this threshold
    pub importance_threshold: f64,

    /// Cluster paths with similar half path lengths
    pub cluster_paths: bool,

    /// Consider unique scattering atoms only in degeneracy calculation
    pub unique_scatterers_only: bool,
}

impl Default for PathFinderConfig {
    fn default() -> Self {
        Self {
            max_path_length: 10.0,
            max_paths: 100,
            max_legs: 8,
            importance_threshold: 0.001,
            cluster_paths: true,
            unique_scatterers_only: true,
        }
    }
}

/// A wrapper for Path with ordering based on importance for use in priority queues
#[derive(Clone)]
struct RankedPath {
    /// The wrapped path
    path: Path,
}

impl RankedPath {
    /// Creates a new RankedPath from a Path
    fn new(path: Path) -> Self {
        Self { path }
    }
}

impl PartialEq for RankedPath {
    fn eq(&self, other: &Self) -> bool {
        self.path.importance == other.path.importance
    }
}

impl Eq for RankedPath {}

impl PartialOrd for RankedPath {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for RankedPath {
    fn cmp(&self, other: &Self) -> Ordering {
        // Order by importance (higher importance first)
        // This creates a max-heap when used with BinaryHeap
        other
            .path
            .importance
            .partial_cmp(&self.path.importance)
            .unwrap_or(Ordering::Equal)
    }
}

/// Handles generation and filtering of scattering paths
///
/// The PathFinder is responsible for enumerating potential scattering paths,
/// evaluating their importance, and filtering them based on various criteria.
pub struct PathFinder {
    /// Configuration parameters for path finding
    config: PathFinderConfig,

    /// Atomic structure for which paths are being found
    structure: AtomicStructure,

    /// Index of the absorbing atom
    absorber_index: usize,

    /// Cache of nearest neighbors for each atom
    neighbor_cache: Vec<Vec<(usize, f64)>>,
}

impl PathFinder {
    /// Creates a new PathFinder
    ///
    /// # Arguments
    ///
    /// * `structure` - Atomic structure for path finding
    /// * `absorber_index` - Index of the absorbing atom
    /// * `config` - Path finder configuration
    ///
    /// # Returns
    ///
    /// A new PathFinder instance
    pub fn new(
        structure: AtomicStructure,
        absorber_index: usize,
        config: PathFinderConfig,
    ) -> Self {
        let atom_count = structure.atom_count();
        let mut neighbor_cache = Vec::with_capacity(atom_count);

        // Initialize empty neighbor lists for each atom
        for _ in 0..atom_count {
            neighbor_cache.push(Vec::new());
        }

        Self {
            config,
            structure,
            absorber_index,
            neighbor_cache,
        }
    }

    /// Precomputes neighbor atoms for each atom in the structure
    ///
    /// This speeds up path finding by caching distances between atoms.
    ///
    /// # Arguments
    ///
    /// * `max_distance` - Maximum distance to consider for neighbors
    fn precompute_neighbors(&mut self, max_distance: f64) {
        let atom_count = self.structure.atom_count();

        for i in 0..atom_count {
            let mut neighbors = Vec::new();

            for j in 0..atom_count {
                if i == j {
                    continue;
                }

                let pos_i = *self.structure.atom(i).unwrap().position();
                let pos_j = *self.structure.atom(j).unwrap().position();
                let distance = (pos_j - pos_i).length();

                if distance <= max_distance {
                    neighbors.push((j, distance));
                }
            }

            // Sort neighbors by distance (closest first)
            neighbors.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(Ordering::Equal));

            self.neighbor_cache[i] = neighbors;
        }
    }

    /// Finds all scattering paths that satisfy the configuration criteria
    ///
    /// # Returns
    ///
    /// A vector of paths ordered by importance
    pub fn find_paths(&mut self) -> Vec<Path> {
        // Precompute neighbors for efficiency
        self.precompute_neighbors(self.config.max_path_length / 2.0);

        // Generate and score paths
        let mut paths = self.enumerate_paths();

        // Filter and order paths by importance
        paths.sort_by(|a, b| {
            b.importance
                .partial_cmp(&a.importance)
                .unwrap_or(Ordering::Equal)
        });

        // Limit the number of paths
        if paths.len() > self.config.max_paths {
            paths.truncate(self.config.max_paths);
        }

        paths
    }

    /// Calculates the importance factor for a path
    ///
    /// The importance factor determines how much a path contributes to the
    /// EXAFS signal. It depends on the path length, the number of legs,
    /// and the scattering strength of atoms in the path.
    ///
    /// # Arguments
    ///
    /// * `path` - The path to evaluate
    ///
    /// # Returns
    ///
    /// The calculated importance factor
    fn calculate_importance(&self, path: &Path) -> f64 {
        // Base importance is inversely proportional to path length
        let length_factor = 1.0 / path.total_length;

        // Paths with more legs are generally less important
        let leg_factor = match path.path_type {
            PathType::SingleScattering => 1.0,
            PathType::DoubleScattering => 0.5,
            PathType::Triangle => 0.7,
            PathType::MultipleScattering => 0.3 / (path.legs.len() as f64 - 2.0),
        };

        // Scattering strength depends on atomic number (Z)
        // Higher Z atoms scatter more strongly
        let mut scattering_factor = 0.0;
        for leg in &path.legs {
            let atom = self.structure.atom(leg.to_atom).unwrap();
            // Simple Z-dependent scattering approximation
            let z_factor = atom.atomic_number() as f64 / 20.0;
            scattering_factor += z_factor.min(1.0);
        }
        scattering_factor /= path.legs.len() as f64;

        // Combine factors
        let importance = length_factor * leg_factor * scattering_factor;

        // Account for degeneracy
        importance * path.degeneracy as f64
    }

    /// Enumerates all potential scattering paths
    ///
    /// This method uses a heap-based approach to efficiently explore the
    /// path space, similar to the approach used in the original FEFF code.
    ///
    /// # Returns
    ///
    /// A vector of paths that pass the initial filtering criteria
    fn enumerate_paths(&self) -> Vec<Path> {
        let mut path_heap = BinaryHeap::new();
        let mut final_paths = Vec::new();
        let mut visited_paths: HashSet<Vec<usize>> = HashSet::new();

        // First, generate all single scattering paths
        self.generate_single_scattering_paths(&mut path_heap, &mut visited_paths);

        // Process paths from the heap until it's empty or we've found enough paths
        while let Some(ranked_path) = path_heap.pop() {
            let path = ranked_path.path;

            // Add this path to our final set if it passes the threshold
            if path.importance >= self.config.importance_threshold {
                final_paths.push(path.clone());

                if final_paths.len() >= self.config.max_paths {
                    break;
                }
            }

            // Generate extended paths if we haven't reached the max legs
            if path.legs.len() < self.config.max_legs {
                self.extend_path(&path, &mut path_heap, &mut visited_paths);
            }
        }

        final_paths
    }

    /// Generates all single scattering paths and adds them to the path heap
    ///
    /// # Arguments
    ///
    /// * `path_heap` - Heap for ordering paths by importance
    /// * `visited_paths` - Set to track already processed paths
    fn generate_single_scattering_paths(
        &self,
        path_heap: &mut BinaryHeap<RankedPath>,
        visited_paths: &mut HashSet<Vec<usize>>,
    ) {
        // Get neighbor atoms of the absorber
        for &(scatterer_index, distance) in &self.neighbor_cache[self.absorber_index] {
            // Skip if the total path would be too long
            if distance * 2.0 > self.config.max_path_length {
                continue;
            }

            // Create a single scattering path
            let mut path = Path::create_single_scattering_path(
                self.absorber_index,
                scatterer_index,
                &self.structure,
            );

            // Calculate importance
            let importance = self.calculate_importance(&path);
            path.set_importance(importance);

            // Add to visited set
            let path_key = path.atom_sequence.clone();
            visited_paths.insert(path_key);

            // Add to heap
            path_heap.push(RankedPath::new(path));
        }
    }

    /// Extends a path by adding additional legs and adds resulting paths to the heap
    ///
    /// # Arguments
    ///
    /// * `base_path` - Path to extend
    /// * `path_heap` - Heap for ordering paths by importance
    /// * `visited_paths` - Set to track already processed paths
    fn extend_path(
        &self,
        base_path: &Path,
        path_heap: &mut BinaryHeap<RankedPath>,
        visited_paths: &mut HashSet<Vec<usize>>,
    ) {
        // Get the last atom in the path
        let last_atom = base_path.legs.last().unwrap().to_atom;

        // Try extending the path with each neighbor of the last atom
        for &(next_atom, _distance) in &self.neighbor_cache[last_atom] {
            // Skip if this would create a 2-leg loop (except back to absorber)
            if base_path.legs.len() >= 2 {
                let previous_atom = base_path.legs[base_path.legs.len() - 2].from_atom;
                if next_atom == previous_atom && next_atom != self.absorber_index {
                    continue;
                }
            }

            // Check if adding this leg would exceed max path length
            let atom1 = self.structure.atom(last_atom).unwrap();
            let atom2 = self.structure.atom(next_atom).unwrap();
            let new_leg_length = atom1.distance_to(atom2);
            let potential_path_length = base_path.total_length + new_leg_length;

            // If we're not going back to the absorber, we need to add the return leg too
            let potential_total_length = if next_atom != self.absorber_index {
                let next_atom_obj = self.structure.atom(next_atom).unwrap();
                let absorber_atom = self.structure.atom(self.absorber_index).unwrap();
                let return_leg_length = next_atom_obj.distance_to(absorber_atom);
                potential_path_length + return_leg_length
            } else {
                potential_path_length
            };

            if potential_total_length > self.config.max_path_length {
                continue;
            }

            // Create the new leg
            let new_leg = PathLeg::new(last_atom, next_atom, &self.structure);

            // Build the extended path
            let mut new_legs = base_path.legs.clone();
            new_legs.push(new_leg);

            // If we're not back at the absorber, add a final leg to return
            if next_atom != self.absorber_index {
                let final_leg = PathLeg::new(next_atom, self.absorber_index, &self.structure);
                new_legs.push(final_leg);
            }

            // Create the new path
            let mut new_path = Path::new(new_legs, self.absorber_index);

            // Check if we've already visited this path
            let path_key = new_path.atom_sequence.clone();
            if visited_paths.contains(&path_key) {
                continue;
            }

            // Calculate importance
            let importance = self.calculate_importance(&new_path);
            new_path.set_importance(importance);

            // Skip if below threshold
            if importance < self.config.importance_threshold {
                continue;
            }

            // Add to visited set
            visited_paths.insert(path_key);

            // Add to heap
            path_heap.push(RankedPath::new(new_path));
        }
    }
}
