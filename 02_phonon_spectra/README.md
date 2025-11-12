# Phonon Spectra Calculations for γ-Ga2O3

This directory contains all input and output files for phonon spectra calculations of gamma-Ga2O3, including phonon dispersion relations for Ga2, Ga3, and Ga4 structures.

## Key File Descriptions

- **Phonon dispersion data** - omega2.out, x_omega.dat, y_omega.dat are data files related to phonon spectra
- **Input files** - kpoints.in, basis.in, run.in, and model.xyz are input files for phonon calculations

## Calculation Methods

Phonon spectra calculations are based on the trained neuroevolution potential (NEP) file `nep.txt` from the `01_nep_training` directory.

## Plotting Scripts

The plotting scripts used for phonon spectra visualization in this work are based on the code from the [tang070205/tools](https://github.com/tang070205/tools/blob/main) repository and have been modified according to the specific requirements of this study.
