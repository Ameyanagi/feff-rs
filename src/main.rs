/*
MIT License with FEFF10 Attribution

Copyright (c) 2025 Ameyanagi

Based on or developed using Distribution: FEFF10.0
Copyright (c) 2020 FEFF Project, University of Washington and SLAC National Accelerator Laboratory.
All rights reserved.
*/

//! Main executable for FEFF-rs

fn main() -> anyhow::Result<()> {
    // Initialize logging
    env_logger::init();

    println!("FEFF-rs v{}", feff_rs::VERSION);
    println!("This is a Rust implementation of the FEFF code for XAS calculations");
    println!("Based on FEFF10 from the FEFF Project at UW and SLAC");
    println!("-----------------------------------------------------------");
    println!("This software is currently under development and not ready for production use.\n");

    let feff = feff_rs::Feff::new();
    feff.run()?;

    Ok(())
}
