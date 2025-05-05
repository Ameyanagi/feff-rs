/*
MIT License with FEFF10 Attribution

Copyright (c) 2025 Ameyanagi

Based on or developed using Distribution: FEFF10.0
Copyright (c) 2020 FEFF Project, University of Washington and SLAC National Accelerator Laboratory.
All rights reserved.
*/

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use feff_rs::utils::{angstrom_to_bohr, bohr_to_angstrom, ev_to_hartree, hartree_to_ev};

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

criterion_group!(benches, unit_conversion_benchmark);
criterion_main!(benches);