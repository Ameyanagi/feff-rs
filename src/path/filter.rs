/*
MIT License with FEFF10 Attribution

Copyright (c) 2025 Ameyanagi

Based on or developed using Distribution: FEFF10.0
Copyright (c) 2020 FEFF Project, University of Washington and SLAC National Accelerator Laboratory.
All rights reserved.
*/

//! Path filtering mechanisms
//!
//! This module implements filters for selecting important scattering paths
//! based on various criteria such as path length, scattering angle, and amplitude.

use crate::path::path::{Path, PathType};

/// Filter configuration for path selection
///
/// This configuration controls how paths are selected and filtered
/// for use in EXAFS calculations.
#[derive(Debug, Clone)]
pub struct PathFilterConfig {
    /// Maximum path length to consider (in Å)
    pub max_path_length: f64,
    
    /// Minimum path length to consider (in Å)
    pub min_path_length: f64,
    
    /// Minimum importance factor to keep a path
    pub min_importance: f64,
    
    /// Maximum number of paths to keep
    pub max_path_count: usize,
    
    /// Filter curved paths more strictly
    pub curve_filter: bool,
    
    /// Curve filter parameter (higher = stricter filtering)
    pub curve_parameter: f64,
}

impl Default for PathFilterConfig {
    fn default() -> Self {
        Self {
            max_path_length: 10.0,
            min_path_length: 1.0,
            min_importance: 0.001,
            max_path_count: 100,
            curve_filter: true,
            curve_parameter: 2.0,
        }
    }
}

/// Filters a set of paths based on specified criteria
///
/// This function applies various filters to select the most important
/// paths for EXAFS calculations.
///
/// # Arguments
///
/// * `paths` - Vector of paths to filter
/// * `config` - Filter configuration
///
/// # Returns
///
/// A vector of filtered paths
pub fn filter_paths(mut paths: Vec<Path>, config: &PathFilterConfig) -> Vec<Path> {
    // Filter paths by length
    paths.retain(|path| {
        path.total_length >= config.min_path_length && 
        path.total_length <= config.max_path_length
    });
    
    // Apply curve filter for multiple scattering paths
    if config.curve_filter {
        paths.retain(|path| {
            match path.path_type {
                PathType::SingleScattering => true,  // Keep all single scattering paths
                PathType::DoubleScattering | PathType::Triangle | PathType::MultipleScattering => {
                    // Calculate curved path penalty
                    let penalty = curve_penalty(path, config.curve_parameter);
                    
                    // Adjust importance by the penalty factor
                    path.importance / penalty >= config.min_importance
                }
            }
        });
    } else {
        // Simple importance filter
        paths.retain(|path| path.importance >= config.min_importance);
    }
    
    // Sort paths by importance (high to low)
    paths.sort_by(|a, b| b.importance.partial_cmp(&a.importance).unwrap_or(std::cmp::Ordering::Equal));
    
    // Limit the number of paths
    if paths.len() > config.max_path_count {
        paths.truncate(config.max_path_count);
    }
    
    paths
}

/// Calculates a penalty factor for curved paths
///
/// This function applies a penalty to curved paths based on their
/// deviation from a straight line. The penalty increases with the
/// number of legs and the amount of curvature.
///
/// # Arguments
///
/// * `path` - The path to evaluate
/// * `curve_parameter` - Curve sensitivity parameter
///
/// # Returns
///
/// A penalty factor (1.0 = no penalty, higher = more penalty)
fn curve_penalty(path: &Path, curve_parameter: f64) -> f64 {
    match path.path_type {
        PathType::SingleScattering => 1.0,  // No penalty for single scattering
        PathType::DoubleScattering | PathType::Triangle => {
            // Triangle paths generally have slightly higher penalty
            if path.path_type == PathType::Triangle {
                1.0 + curve_parameter * 0.5
            } else {
                // For double scattering, penalty depends on how close to a straight line
                // A perfect double scattering path is a straight line (180° angle)
                // so the penalty is minimal
                1.0 + curve_parameter * 0.2
            }
        },
        PathType::MultipleScattering => {
            // For multiple scattering, penalty increases with the number of legs
            let leg_penalty = (path.legs.len() as f64 - 2.0) * 0.5;
            
            // Overall penalty combines leg count and curvature
            1.0 + curve_parameter * leg_penalty
        }
    }
}

/// Filters paths by clustering similar path lengths
///
/// This function groups paths with similar effective lengths and
/// keeps only the most important path from each group. This helps
/// to avoid oversampling a particular path length region.
///
/// # Arguments
///
/// * `paths` - Vector of paths to filter
/// * `resolution` - Resolution for path length clustering (in Å)
///
/// # Returns
///
/// A vector of filtered paths
pub fn cluster_paths_by_length(mut paths: Vec<Path>, resolution: f64) -> Vec<Path> {
    // Sort paths by length
    paths.sort_by(|a, b| a.total_length.partial_cmp(&b.total_length).unwrap_or(std::cmp::Ordering::Equal));
    
    let mut clustered_paths = Vec::new();
    let mut current_cluster = Vec::new();
    let mut current_length = if !paths.is_empty() { paths[0].total_length } else { 0.0 };
    
    for path in paths {
        // Check if this path belongs to the current cluster
        if (path.total_length - current_length).abs() <= resolution {
            current_cluster.push(path);
        } else {
            // This path starts a new cluster
            if !current_cluster.is_empty() {
                // Find the most important path in the current cluster
                let best_path = current_cluster.iter()
                    .max_by(|a, b| a.importance.partial_cmp(&b.importance).unwrap_or(std::cmp::Ordering::Equal))
                    .unwrap()
                    .clone();
                
                clustered_paths.push(best_path);
            }
            
            // Start a new cluster
            current_cluster = vec![path];
            current_length = current_cluster[0].total_length;
        }
    }
    
    // Process the last cluster
    if !current_cluster.is_empty() {
        let best_path = current_cluster.iter()
            .max_by(|a, b| a.importance.partial_cmp(&b.importance).unwrap_or(std::cmp::Ordering::Equal))
            .unwrap()
            .clone();
        
        clustered_paths.push(best_path);
    }
    
    clustered_paths
}

/// Selects a subset of paths optimized for spectral coverage
///
/// This function selects paths to ensure good coverage of the
/// EXAFS spectrum while keeping the total number of paths manageable.
///
/// # Arguments
///
/// * `paths` - Vector of paths to select from
/// * `max_count` - Maximum number of paths to select
/// * `k_max` - Maximum k-value for EXAFS calculations (in Å⁻¹)
///
/// # Returns
///
/// A vector of selected paths
pub fn select_optimal_path_set(mut paths: Vec<Path>, max_count: usize, k_max: f64) -> Vec<Path> {
    // If we already have fewer paths than the maximum, return all paths
    if paths.len() <= max_count {
        return paths;
    }
    
    // Sort paths by importance
    paths.sort_by(|a, b| b.importance.partial_cmp(&a.importance).unwrap_or(std::cmp::Ordering::Equal));
    
    // Calculate the Nyquist sampling interval for path lengths
    // This is π/k_max in the EXAFS literature
    let nyquist_interval = std::f64::consts::PI / k_max;
    
    // Create bins for path lengths based on the Nyquist interval
    let max_length = paths.iter().map(|p| p.total_length).fold(0.0, f64::max);
    let num_bins = (max_length / nyquist_interval).ceil() as usize;
    
    // Create empty bins
    let mut bins: Vec<Vec<Path>> = (0..num_bins).map(|_| Vec::new()).collect();
    
    // Distribute paths to bins
    for path in paths {
        let bin_index = (path.total_length / nyquist_interval).floor() as usize;
        if bin_index < bins.len() {
            bins[bin_index].push(path);
        } else {
            // This shouldn't happen, but just in case
            bins.last_mut().unwrap().push(path);
        }
    }
    
    // Select paths from each bin
    let mut selected_paths = Vec::new();
    
    // First pass: take the most important path from each bin
    for bin in &mut bins {
        if bin.is_empty() {
            continue;
        }
        
        // Sort the bin by importance
        bin.sort_by(|a, b| b.importance.partial_cmp(&a.importance).unwrap_or(std::cmp::Ordering::Equal));
        
        // Take the most important path
        selected_paths.push(bin[0].clone());
    }
    
    // Second pass: fill remaining slots with the most important remaining paths
    if selected_paths.len() < max_count {
        // Collect all remaining paths
        let mut remaining_paths = Vec::new();
        for bin in &bins {
            for path in bin.iter().skip(1) {
                remaining_paths.push(path.clone());
            }
        }
        
        // Sort by importance
        remaining_paths.sort_by(|a, b| b.importance.partial_cmp(&a.importance).unwrap_or(std::cmp::Ordering::Equal));
        
        // Take as many as needed to reach max_count
        let remaining_slots = max_count - selected_paths.len();
        selected_paths.extend(remaining_paths.into_iter().take(remaining_slots));
    }
    
    // Final sort by importance
    selected_paths.sort_by(|a, b| b.importance.partial_cmp(&a.importance).unwrap_or(std::cmp::Ordering::Equal));
    
    selected_paths
}