# tn_fnds_PBR - Rust Port of tn_fnds_PB UTAU Resampler

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

Rust implementation of the tn_fnds_PB vocal synthesis engine for UTAU, combining Masanori Morise's WORLD vocoder with Zteer's tn_fnds resampler.

## Features

- Rust implementation of the original resampler
- Supports standard tn_fnds flags (g, t, B, A, O, e)
- Updated W flag
- Improved memory safety and performance from Rust implementation

## Missing Features (TODO)

This Rust port is currently missing some functionality from the original implementation:

- `createWaveSpec()` - Wave spectrogram creation (used when B > 50 or O ≠ 0)
- `Opening()` - Vocal tract opening effect (used when O ≠ 0)
- `rebuildWave()` - Wave reconstruction (used when O ≠ 0) 
- `breath2()` - Breathiness effect (used when B > 50)

These features will be implemented in future versions.

## Requirements

- Rust 1.60+ (https://www.rust-lang.org/)
- .frq files (preferably generated with fresamp)
- Input WAV files in UTAU-compatible format

## Installation

```sh
git clone https://github.com/Arcuslok/tn_fnds_PBR.git
cd tn_fnds_PBR
cargo build --release
