# Spacecraft ADCS — Course Project                                                                                                      
**NASA Starling Mission Analysis**                                                                                                      
*16.S897 Spacecraft Attitude Determination and Control, MIT Spring 2026*                                                                
*Aleksander Garbuz · Israel Garcia* 

Analysis of NASA's Starling 6U CubeSat (BCT XB6 bus, SSO at 480 km).

## Structure

Each homework is a self-contained folder:

```
report1/   — spacecraft model, orbital dynamics, Euler's equation
report2/   — safe-mode gyrostat, sensor models
report3/   — MEKF, gyro bias estimation
report4/   — actuators, disturbance torques, attitude control, eigen-axis slew
sources/   — shared Julia files, HW specs
import/    — development/scratch notebooks and Julia modules
```

Each folder contains a Julia notebook, LaTeX source, and `figs/`.

## Running

Notebooks use a Julia kernel

Run notebooks from inside their folder so `include()` and `savefig()` paths resolve.

## Compiling

```bash
cd report4
pdflatex report4.tex && pdflatex report4.tex
```
