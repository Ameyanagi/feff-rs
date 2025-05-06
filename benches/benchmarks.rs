/*
MIT License with FEFF10 Attribution

Copyright (c) 2025 Ameyanagi

Based on or developed using Distribution: FEFF10.0
Copyright (c) 2020 FEFF Project, University of Washington and SLAC National Accelerator Laboratory.
All rights reserved.
*/

use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use feff_rs::atoms::{Atom, AtomicStructure, PotentialType, Vector3D};
use feff_rs::fms::{FmsMatrix, FmsParameters, FmsSolver, SolverMethod, XanesCalculator};
use feff_rs::potential::MuffinTinPotential;
use feff_rs::scattering::{
    calculate_phase_shifts, calculate_scattering_matrices_old, ScatteringMatrixResults,
};
use feff_rs::utils::matrix::{
    compute_free_greens_matrix, compute_path_operator, compute_structure_factor, compute_t_matrix,
};
use feff_rs::utils::{angstrom_to_bohr, bohr_to_angstrom, ev_to_hartree, hartree_to_ev};
use feff_rs::xas::{calculate_with_core_hole, CoreHoleConfig, CoreHoleMethod};
use ndarray::Array2;
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

fn fms_matrix_benchmark(c: &mut Criterion) {
    let mut group = c.benchmark_group("FMS Matrix Construction");

    // Test with structures of different sizes
    for size in [10, 50, 100] {
        let structure = create_benchmark_structure(size);
        let radius = 6.0; // FMS radius in Angstroms

        group.bench_with_input(
            BenchmarkId::new("fms_matrix_construction", size),
            &structure,
            |b, s| {
                b.iter(|| {
                    let fms_matrix = FmsMatrix::new(black_box(s), black_box(radius)).unwrap();
                    black_box(&fms_matrix);
                })
            },
        );
    }

    group.finish();
}

fn fms_solver_benchmark(c: &mut Criterion) {
    let mut group = c.benchmark_group("FMS Solver Methods");

    // Create test matrix and parameters
    let n = 100; // Matrix size
    let mut test_matrix = Array2::<Complex64>::zeros((n, n));

    // Fill matrix with realistic data to simulate FMS matrix
    for i in 0..n {
        for j in 0..n {
            let dist = ((i as f64 - j as f64).powi(2)).sqrt();
            if i == j {
                test_matrix[(i, j)] = Complex64::new(1.0, 0.0);
            } else {
                test_matrix[(i, j)] = Complex64::new(0.1 / (1.0 + dist), 0.01 / (1.0 + dist));
            }
        }
    }

    // Benchmark each solver method
    for &method in &[
        SolverMethod::LuDecomposition,
        SolverMethod::IterativeCGS,
        SolverMethod::BlockDiagonal,
    ] {
        group.bench_with_input(
            BenchmarkId::new("fms_solver", format!("{:?}", method)),
            &(&test_matrix, method),
            |b, (matrix, method)| {
                b.iter(|| {
                    let solver = FmsSolver::new(black_box(*method));
                    black_box(solver.solve(black_box(matrix)).unwrap());
                })
            },
        );
    }

    group.finish();
}

fn xanes_calculator_benchmark(c: &mut Criterion) {
    let mut group = c.benchmark_group("XANES Calculator");

    // Create test structure and path operator for XANES calculation
    let structure = create_benchmark_structure(50);
    let n = 16 * structure.atom_count(); // Size based on l_max=3 and atom count

    // Create a sample path operator matrix
    let mut path_operator = Array2::<Complex64>::zeros((n, n));

    // Fill with realistic data
    for i in 0..n {
        for j in 0..n {
            if i == j {
                path_operator[(i, j)] = Complex64::new(1.0, 0.1);
            } else {
                let dist = ((i as f64 - j as f64).powi(2)).sqrt();
                path_operator[(i, j)] = Complex64::new(0.05 / (1.0 + dist), 0.01 / (1.0 + dist));
            }
        }
    }

    // Create energy grid for spectrum calculation
    let energy_min = 0.0;
    let energy_max = 100.0;
    let energy_step = 1.0;
    let mut energies = Vec::new();
    let mut e = energy_min;
    while e <= energy_max {
        energies.push(e);
        e += energy_step;
    }

    // Benchmark single-point XANES calculation
    group.bench_with_input(
        BenchmarkId::new("xanes_single_point", structure.atom_count()),
        &(&structure, &path_operator, 50.0),
        |b, (structure, path_operator, energy)| {
            b.iter(|| {
                let calculator = XanesCalculator::new(black_box(*structure), 1.0);
                black_box(
                    calculator
                        .calculate_xanes(black_box(*energy), black_box(path_operator))
                        .unwrap(),
                );
            })
        },
    );

    // Benchmark batch XANES calculation using sequential processing
    group.bench_with_input(
        BenchmarkId::new("xanes_spectrum_sequential", energies.len()),
        &(&structure, &path_operator, &energies),
        |b, (structure, path_operator, energies)| {
            b.iter(|| {
                let calculator = XanesCalculator::new(black_box(*structure), 1.0);

                // Calculate each point sequentially
                let mut results = Vec::with_capacity(energies.len());
                for &energy in *energies {
                    let value = calculator
                        .calculate_xanes(black_box(energy), black_box(path_operator))
                        .unwrap();
                    results.push(value);
                }

                black_box(results)
            })
        },
    );

    // Benchmark batch XANES calculation using parallel processing
    group.bench_with_input(
        BenchmarkId::new("xanes_spectrum_parallel", energies.len()),
        &(&structure, &path_operator, &energies),
        |b, (structure, path_operator, energies)| {
            b.iter(|| {
                let calculator = XanesCalculator::new(black_box(*structure), 1.0);

                // Use the batch calculation method
                black_box(
                    calculator
                        .calculate_xanes_spectrum(black_box(energies), black_box(path_operator))
                        .unwrap(),
                )
            })
        },
    );

    group.finish();
}

fn fms_end_to_end_benchmark(c: &mut Criterion) {
    let mut group = c.benchmark_group("FMS End-to-End");

    // Simplified version of scattering matrix results
    struct MockScatteringResults {
        pub green_matrix: Array2<Complex64>,
        pub global_t_matrix: Array2<Complex64>,
        pub max_l: usize,
        pub energy: f64,
    }

    // Test with different structure sizes
    for size in [10, 20] {
        let structure = create_benchmark_structure(size);
        let atom_count = structure.atom_count();
        let max_l = 3;
        let l_size = (max_l + 1) * (max_l + 1);
        let n = atom_count * l_size;

        // Create mock scattering matrix results
        let green_matrix = Array2::<Complex64>::zeros((n, n));
        let global_t_matrix = Array2::<Complex64>::zeros((n, n));

        // Fill with realistic data
        for i in 0..n {
            for j in 0..n {
                if i == j {
                    green_matrix[(i, j)] = Complex64::new(0.0, 0.0);
                    global_t_matrix[(i, j)] = Complex64::new(0.2, 0.1);
                } else {
                    let atom_i = i / l_size;
                    let atom_j = j / l_size;
                    if atom_i == atom_j {
                        green_matrix[(i, j)] = Complex64::new(0.0, 0.0);
                        global_t_matrix[(i, j)] = Complex64::new(0.05, 0.02);
                    } else {
                        let dist = (atom_i as f64 - atom_j as f64).abs() + 1.0;
                        green_matrix[(i, j)] = Complex64::new(0.1 / dist, 0.05 / dist);
                        global_t_matrix[(i, j)] = Complex64::new(0.0, 0.0);
                    }
                }
            }
        }

        let mock_results = ScatteringMatrixResults {
            phase_shifts: Vec::new(), // Not used in the benchmark
            green_matrix,
            global_t_matrix,
            max_l,
            energy: 100.0,
        };

        // Create FMS parameters for benchmark
        let mut parameters = FmsParameters::default();
        parameters.radius = 6.0;
        parameters.solver_method = SolverMethod::LuDecomposition;

        // Add energy grid
        let energy_min = 95.0;
        let energy_max = 105.0;
        let energy_step = 1.0;
        let mut energies = Vec::new();
        let mut e = energy_min;
        while e <= energy_max {
            energies.push(e);
            e += energy_step;
        }
        parameters.energies = energies;

        // Benchmark FMS calculation
        group.bench_with_input(
            BenchmarkId::new("fms_calculation", atom_count),
            &(&structure, &mock_results, parameters.clone()),
            |b, (structure, mock_results, params)| {
                b.iter(|| {
                    black_box(
                        feff_rs::fms::calculate_fms(
                            black_box(*structure),
                            black_box(mock_results),
                            black_box(params.clone()),
                        )
                        .unwrap(),
                    )
                })
            },
        );

        // Set different solver methods for comparison
        let solver_methods = [
            SolverMethod::LuDecomposition,
            SolverMethod::IterativeCGS,
            SolverMethod::BlockDiagonal,
        ];

        for &method in &solver_methods {
            let mut params = parameters.clone();
            params.solver_method = method;

            group.bench_with_input(
                BenchmarkId::new(format!("fms_method_{:?}", method), atom_count),
                &(&structure, &mock_results, params.clone()),
                |b, (structure, mock_results, params)| {
                    b.iter(|| {
                        black_box(
                            feff_rs::fms::calculate_fms(
                                black_box(*structure),
                                black_box(mock_results),
                                black_box(params.clone()),
                            )
                            .unwrap(),
                        )
                    })
                },
            );
        }
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
    simd_matrix_benchmark,
    fms_matrix_benchmark,
    fms_solver_benchmark,
    xanes_calculator_benchmark,
    fms_end_to_end_benchmark
);
criterion_main!(benches);
