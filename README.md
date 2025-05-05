# FEFF-rs

A modern Rust implementation of the FEFF code for calculating X-ray absorption spectroscopy (XAS) and related spectroscopies.

## Overview

FEFF-rs aims to provide a high-performance, memory-efficient, and maintainable implementation of the FEFF code in Rust. FEFF is a widely used ab initio multiple-scattering code for calculations of excitation spectra and electronic structure.

This project is under active development and is not yet ready for production use.

## Features (Planned)

- Parse FEFF input files
- Calculate atomic potentials
- Find scattering paths
- Compute X-ray Absorption Near Edge Structure (XANES)
- Compute Extended X-ray Absorption Fine Structure (EXAFS)
- Support multiple spectroscopies (XANES, EXAFS, EELS, XES, etc.)
- Parallel processing for improved performance
- Modern CLI with good error messages

## Project Status

This project is in early development. See the [architecture document](plans/architecture.md) for the planned design.

## License

MIT License with FEFF10 Attribution

Copyright (c) 2025 Ameyanagi

This project is based on the FEFF10 code developed by the FEFF Project at the University of Washington and SLAC National Accelerator Laboratory. See the [LICENSE](LICENSE) and [FEFF10_LICENSE](FEFF10_LICENSE) files for details.

## Citation

When using this software, please cite both this implementation and the original FEFF papers:

```
FEFF-rs: A Rust implementation of FEFF
[Citation information to be added]

The original FEFF project:
J.J. Rehr et al., Phys. Rev. B, 80, 115112 (2009)
```

## Acknowledgments

We acknowledge the FEFF Project at the University of Washington and SLAC National Accelerator Laboratory for developing the original FEFF code.