## Key File Descriptions

Each directional folder contains:
- **`run.in`** - Input parameter file for GPUMD simulations
- **`kappa.out`** - Output file containing thermal conductivity results

## Calculation Methods

Thermal conductivity calculations are based on the homogeneous non-equilibrium molecular dynamics method, using the trained potential file `nep.txt` from the `01_nep_training` directory. The results are averaged over multiple independent runs.

## Results Analysis

By comparing the thermal conductivity of different sites (Ga2/Ga3/Ga4) and different crystallographic directions ([010]/[100]/[001]), we can analyze the thermal transport anisotropy of γ-Ga2O3 and the influence of different atomic sites on thermal conductivity.
