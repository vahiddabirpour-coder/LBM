# LBM: Natural Convection in a Cavity with Porous Media

Fortran code based on the Lattice Boltzmann Method (LBM) for natural convection in a square cavity, with a clear fluid, a porous medium, or a bi-disperse porous medium.

## Model
- D2Q9 lattice for flow and D2Q9 lattice for temperature (BGK collision)
- Boussinesq buoyancy with Guo forcing; Zou-He boundary conditions
- Porous medium: generalized Darcy-Brinkman-Forchheimer model (Guo & Zhao, 2002)
- Bi-disperse geometry: porous blocks separated by clear-fluid gaps
- Optional run with temperature-dependent viscosity

## Files
| File | Purpose |
|---|---|
| `Source1.f90` | Main program, input reading, initialization, streaming and collision |
| `BC.f90` | Velocity and temperature boundary conditions |
| `result.f90` | Nusselt number, convergence checks, Tecplot output |
| `geomakerf90.f90` | Bi-disperse porous geometry |
| `Input.txt` | Physical parameters (Pr, Ra, Da, geometry) and run number |

## Build and run
```
gfortran -O2 -std=legacy -w -ffree-line-length-none -o lbm Source1.f90 BC.f90 result.f90 geomakerf90.f90
printf '1\n81,81\n1\n\n' | ./lbm
```
Inputs requested at start: run type (1: natural convection, 2: temperature-dependent viscosity), grid size `nx,ny`, and medium type (1: clear fluid, 2: porous, 3: bi-disperse).

## Output
Files are tagged with the run number set in `Input.txt`:
- `Results_Run001_*.txt`: run information and final average Nusselt number
- `Nu_Time_Run001.txt`: Nusselt number versus time
- `Tecplot_Run001_*.DAT`: velocity, pressure and temperature fields (Tecplot format)

## Status
Work in progress. Validation against benchmark results (de Vahl Davis, 1983) is ongoing.
git add README.md
git commit -m "Update README"
git push