/*
MIT License with FEFF10 Attribution

Copyright (c) 2025 Ameyanagi

Based on or developed using Distribution: FEFF10.0
Copyright (c) 2020 FEFF Project, University of Washington and SLAC National Accelerator Laboratory.
All rights reserved.
*/

use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use feff_rs::atoms::{Atom, AtomicStructure, PotentialType, Vector3D};
use feff_rs::potential::MuffinTinPotential;
use feff_rs::scattering::{calculate_phase_shifts, calculate_scattering_matrices_old};
use feff_rs::utils::matrix::{
    compute_free_greens_matrix, compute_path_operator, compute_structure_factor, compute_t_matrix,
};
use feff_rs::utils::{angstrom_to_bohr, bohr_to_angstrom, ev_to_hartree, hartree_to_ev};
use feff_rs::xas::{calculate_with_core_hole, CoreHoleConfig, CoreHoleMethod};
use num_complex::Complex64;

fn unit_conversion_benchmark(c: &mut Criterion) {
    let mut group = c.benchmark_group("Unit Conversions");

    group.bench_function("angstrom_to_bohr", |b| {
        b.iter(|| {
            for i in 0..1000 {
                black_box(angstrom_to_bohr(black_box(i as f64 * 0.1)));
            }
        })
    });

    group.bench_function("bohr_to_angstrom", |b| {
        b.iter(|| {
            for i in 0..1000 {
                black_box(bohr_to_angstrom(black_box(i as f64 * 0.1)));
            }
        })
    });

    group.bench_function("ev_to_hartree", |b| {
        b.iter(|| {
            for i in 0..1000 {
                black_box(ev_to_hartree(black_box(i as f64 * 0.1)));
            }
        })
    });

    group.bench_function("hartree_to_ev", |b| {
        b.iter(|| {
            for i in 0..1000 {
                black_box(hartree_to_ev(black_box(i as f64 * 0.1)));
            }
        })
    });

    group.finish();
}

/// Create a structure of a given size for benchmarking
fn create_benchmark_structure(size: usize) -> AtomicStructure {
    // Create a simple test structure with an iron atom at the center
    // surrounded by oxygen atoms arranged in a grid pattern
    let mut structure = AtomicStructure::new();

    // Add potential types
    let fe_potential = PotentialType::new(0, 26).unwrap(); // Iron
    let o_potential = PotentialType::new(1, 8).unwrap(); // Oxygen

    structure.add_potential_type(fe_potential);
    structure.add_potential_type(o_potential);

    // Add central Fe atom
    let fe_atom = Atom::new(26, Vector3D::new(0.0, 0.0, 0.0), 0).unwrap();
    let fe_idx = structure.add_atom(fe_atom);
    structure.set_central_atom(fe_idx).unwrap();

    // Add oxygen atoms in a grid pattern
    let grid_size = (size as f64).sqrt().ceil() as i32;
    let spacing = 2.0; // Angstrom

    for i in 0..grid_size {
        for j in 0..grid_size {
            if structure.atom_count() >= size {
                break;
            }

            // Skip the center position (already occupied by Fe)
            if i == grid_size / 2 && j == grid_size / 2 {
                continue;
            }

            let x = (i - grid_size / 2) as f64 * spacing;
            let y = (j - grid_size / 2) as f64 * spacing;
            let z = 0.0;

            let o_atom = Atom::new(8, Vector3D::new(x, y, z), 1).unwrap();
            structure.add_atom(o_atom);
        }
    }

    // Calculate muffin-tin radii
    structure.calculate_muffin_tin_radii().unwrap();

    structure
}

fn potential_calculation_benchmark(c: &mut Criterion) {
    let mut group = c.benchmark_group("Potential Calculations");

    // Test with structures of different sizes
    for size in [10, 50, 100] {
        let structure = create_benchmark_structure(size);

        group.bench_with_input(BenchmarkId::new("muffin_tin", size), &structure, |b, s| {
            b.iter(|| {
                let calculator = MuffinTinPotential::new(black_box(s)).unwrap();
                black_box(calculator.calculate().unwrap());
            })
        });
    }

    group.finish();
}

fn phase_shift_benchmark(c: &mut Criterion) {
    let mut group = c.benchmark_group("Phase Shift Calculations");

    // Test with structures of different sizes
    for size in [10, 50, 100] {
        let structure = create_benchmark_structure(size);
        let energy = 100.0; // eV
        let max_l = 3;

        group.bench_with_input(
            BenchmarkId::new("phase_shifts", size),
            &structure,
            |b, s| {
                b.iter(|| {
                    black_box(
                        calculate_phase_shifts(black_box(s), black_box(energy), black_box(max_l))
                            .unwrap(),
                    );
                })
            },
        );
    }

    group.finish();
}

fn scattering_benchmark(c: &mut Criterion) {
    let mut group = c.benchmark_group("Scattering Calculations");

    // Test with structures of different sizes
    for size in [10, 50, 100] {
        let structure = create_benchmark_structure(size);
        let energy = 100.0; // eV
        let max_l = 3;

        group.bench_with_input(
            BenchmarkId::new("scattering_matrices", size),
            &structure,
            |b, s| {
                b.iter(|| {
                    black_box(
                        calculate_scattering_matrices_old(
                            black_box(s),
                            black_box(energy),
                            black_box(max_l),
                        )
                        .unwrap(),
                    );
                })
            },
        );
    }

    group.finish();
}

fn core_hole_benchmark(c: &mut Criterion) {
    let mut group = c.benchmark_group("Core-Hole Calculations");

    // Test with structures of different sizes
    for size in [10, 50] {
        let structure = create_benchmark_structure(size);
        let energy = 100.0; // eV
        let max_l = 3;

        let mut config = CoreHoleConfig::new();
        config
            .with_method(CoreHoleMethod::FinalState)
            .with_edge("K")
            .with_screening(0.0);

        group.bench_with_input(
            BenchmarkId::new("core_hole_fsr", size),
            &structure,
            |b, s| {
                b.iter(|| {
                    black_box(
                        calculate_with_core_hole(
                            black_box(s),
                            black_box(energy),
                            black_box(max_l),
                            black_box(&config),
                        )
                        .unwrap(),
                    );
                })
            },
        );
    }

    group.finish();
}

fn matrix_benchmark(c: &mut Criterion) {
    let mut group = c.benchmark_group("Matrix Operations");

    // Test with matrices of different sizes
    for l_max in [3, 6, 9] {
        // Create phase shifts for testing
        let mut phase_shifts = Vec::with_capacity((l_max + 1) as usize);
        for l in 0..=l_max {
            let phase = Complex64::new(0.1 * l as f64, 0.05 * l as f64);
            phase_shifts.push(phase);
        }

        group.bench_with_input(
            BenchmarkId::new("t_matrix", l_max),
            &phase_shifts,
            |b, shifts| {
                b.iter(|| {
                    black_box(compute_t_matrix(black_box(shifts), black_box(l_max)).unwrap());
                })
            },
        );
    }

    group.finish();
}

fn simd_matrix_benchmark(c: &mut Criterion) {
    let mut group = c.benchmark_group("SIMD Matrix Operations");

    // Test with different atom counts to measure SIMD scaling
    for atom_count in [10, 50, 100] {
        // Create test positions
        let mut positions = Vec::with_capacity(atom_count);
        let grid_size = (atom_count as f64).sqrt().ceil() as i32;
        let spacing = 2.0; // Simulate 2 Å spacing

        // Create a grid of atoms
        for i in 0..grid_size {
            for j in 0..grid_size {
                if positions.len() >= atom_count {
                    break;
                }

                let x = (i - grid_size / 2) as f64 * spacing;
                let y = (j - grid_size / 2) as f64 * spacing;
                let z = 0.0;

                positions.push((x, y, z));
            }
        }

        // Parameters for benchmarks
        let l_max = 3;
        let k = 1.0; // Wave vector

        // Benchmark compute_structure_factor
        group.bench_with_input(
            BenchmarkId::new("structure_factor", atom_count),
            &(positions.clone(), k, l_max),
            |b, (pos, k_val, l_max_val)| {
                b.iter(|| {
                    black_box(
                        compute_structure_factor(
                            black_box(pos),
                            black_box(*k_val),
                            black_box(*l_max_val),
                        )
                        .unwrap(),
                    )
                })
            },
        );

        // Benchmark compute_free_greens_matrix
        group.bench_with_input(
            BenchmarkId::new("greens_matrix", atom_count),
            &(positions.clone(), k, l_max),
            |b, (pos, k_val, l_max_val)| {
                b.iter(|| {
                    black_box(
                        compute_free_greens_matrix(
                            black_box(pos),
                            black_box(*k_val),
                            black_box(*l_max_val),
                        )
                        .unwrap(),
                    )
                })
            },
        );

        // Create matrices for path operator benchmark
        let greens = compute_free_greens_matrix(&positions, k, l_max).unwrap();

        // Create phase shifts for t-matrix
        let mut phase_shifts = Vec::with_capacity((l_max + 1) as usize);
        for l in 0..=l_max {
            let phase = Complex64::new(0.1 * l as f64, 0.05 * l as f64);
            phase_shifts.push(phase);
        }

        let t_matrix = compute_t_matrix(&phase_shifts, l_max).unwrap();

        // Benchmark compute_path_operator
        group.bench_with_input(
            BenchmarkId::new("path_operator", atom_count),
            &(greens.clone(), t_matrix.clone()),
            |b, (g, t)| {
                b.iter(|| black_box(compute_path_operator(black_box(g), black_box(t)).unwrap()))
            },
        );
    }

    group.finish();
}

criterion_group!(
    benches,
    unit_conversion_benchmark,
    potential_calculation_benchmark,
    phase_shift_benchmark,
    scattering_benchmark,
    core_hole_benchmark,
    matrix_benchmark,
    simd_matrix_benchmark
);
criterion_main!(benches);
