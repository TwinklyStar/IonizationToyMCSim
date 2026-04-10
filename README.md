# IonizationToyMCSim

A Monte Carlo simulation of muonium laser ionization. It models two-step photoionization of muonium atoms (Mu = μ⁺e⁻) by solving the **Optical Bloch Equations (OBE)** for a 3-level quantum system driven by 122 nm and 355 nm laser pulses. Output is written in CERN ROOT format.

---

## Table of Contents

- [Prerequisites](#prerequisites)
- [Building](#building)
- [Running a Simulation](#running-a-simulation)
- [Macro File Reference](#macro-file-reference)
- [Input Data Format](#input-data-format)
- [Output](#output)
- [Example Analysis](#example-analysis)
- [Physics Overview](#physics-overview)
- [Known Limitations](#known-limitations)
- [Version History](#version-history)

---

## Prerequisites

The following libraries must be installed before building:

| Library | Version | Purpose |
|---|---|---|
| [Eigen3](https://eigen.tuxfamily.org) | any recent | Linear algebra |
| [Boost](https://www.boost.org) | ≥ 1.86.0 | ODE solver (`odeint`), filesystem |
| [ROOT](https://root.cern) | ≥ 6.x | Output file format, histogram/tree I/O |
| C++ compiler | C++17 | e.g. GCC 9+ or Clang 10+ |
| CMake | ≥ 3.26 | Build system |

### Installing on macOS (Homebrew)

```bash
brew install cmake eigen boost root
```

### Installing on Linux (Ubuntu/Debian)

```bash
sudo apt install cmake libeigen3-dev libboost-all-dev
# ROOT must be installed manually from https://root.cern/install/
```

---

## Building

```bash
git clone <repository-url>
cd LaserToyMC

mkdir build && cd build
cmake ..
make -j$(nproc)
```

The executable `IonizationToyMCSim` will be created inside `build/`.

### CMake options

If CMake cannot find Boost automatically (common on systems with non-standard install paths), set the path explicitly:

```bash
cmake .. -DBOOST_ROOT=/path/to/boost
```

---

## Running a Simulation

The program takes a single argument: the path to a **macro file** (`.mac`) that defines all simulation parameters.

```bash
# From the build directory:
./IonizationToyMCSim ../run/ioni_test.mac
```

> **Note:** Paths inside the macro file (e.g. `MuInputFile`, `OutputFile`) are resolved relative to the **working directory** where you launch the executable, not the macro file's location. It is therefore convenient to run from inside the `run/` directory:

```bash
cd run/
../build/IonizationToyMCSim ioni_test.mac
```

Output is written to the path specified by `OutputFile` in the macro. The output directory is not created automatically, so create it first:

```bash
mkdir -p run/data
```

The provided `ioni_test.mac` writes to `data/g2edmIoni_test.root` (relative to `run/`).

---

## Macro File Reference

A macro file is a plain text file. Lines beginning with `#` are comments. Below is a full list of supported commands with their arguments and units.

### Laser Configuration

```
AddLaser122  <E[J]> <FWHM[ns]> <t_peak[ns]> <linewidth[GHz]> \
             <σ_x[mm]> <σ_y[mm]> \
             <off_x[mm]> <off_y[mm]> <off_z[mm]> \
             <yaw[deg]> <pitch[deg]> <roll[deg]> \
             <detuning[GHz]>

AddLaser355  <E[J]> <FWHM[ns]> <t_peak[ns]> <linewidth[GHz]> \
             <σ_x[mm]> <σ_y[mm]> \
             <off_x[mm]> <off_y[mm]> <off_z[mm]> \
             <yaw[deg]> <pitch[deg]> <roll[deg]>
```

Multiple `AddLaser122` and `AddLaser355` commands can be used to add several laser beams. Parameters:

| Parameter | Description |
|---|---|
| `E` | Total pulse energy in Joules |
| `FWHM` | Pulse duration (full width at half maximum) in ns |
| `t_peak` | Time of peak intensity in ns |
| `linewidth` | Laser linewidth (1σ) in GHz |
| `σ_x`, `σ_y` | Gaussian beam radii in mm |
| `off_x/y/z` | Beam center offset from origin in mm |
| `yaw/pitch/roll` | Beam direction Euler angles (x-y-z order) in degrees |
| `detuning` | Frequency detuning from resonance in GHz (122 nm only) |

### Simulation Control

```
SetRunTime      <duration[ns]>      # Total simulation time window
RandomSeed      <integer>           # RNG seed for reproducibility
EventN          <N> | max           # Number of muonium events to simulate
SetDopplerShift <shift[rad/ns]>     # Fix Doppler shift to a constant value (overrides v·k)
```

### Muonium Source

```
InputMuPar   on | off              # on: read from file; off: Monte Carlo sampling
MuInputFile  <path/to/file.dat>    # Path to muonium input data file (see format below)
```

### Output Control

All per-timestep array branches can be toggled individually to reduce file size:

```
RootOutput  t            on | off   # Time array
RootOutput  RabiFreq     on | off
RootOutput  EField       on | off
RootOutput  Intensity122 on | off
RootOutput  Intensity355 on | off
RootOutput  GammaIon     on | off
RootOutput  rho_gg       on | off   # Ground state population
RootOutput  rho_ee       on | off   # Excited state population
RootOutput  rho_ge_r     on | off   # Coherence (real part)
RootOutput  rho_ge_i     on | off   # Coherence (imaginary part)
RootOutput  rho_ion      on | off   # Ionized state population

OutputFile  <path/to/output.root>
```

### Example Macro

```
# 122 nm: 13.5 µJ, 2 ns FWHM, peak at 5 ns, 80 GHz linewidth, 4×1 mm beam
AddLaser122  13.5e-6  2  5  80  4  1  0  0  2  0  0  0  0
# 355 nm: 8 mJ, same timing and geometry
AddLaser355  8e-3     2  5  80  4  1  0  0  2  0  0  0

SetRunTime     10
RandomSeed     999
InputMuPar     on
MuInputFile    ../datasets/test1k.dat
EventN         1000

RootOutput  t         on
RootOutput  rho_ion   on

OutputFile  data/output.root
```

---

## Input Data Format

When `InputMuPar on` is set, the program reads muonium initial conditions from a plain-text file.

```
<number_of_events>
x1  y1  z1  vx1  vy1  vz1
x2  y2  z2  vx2  vy2  vz2
...
```

- Positions in **mm**
- Velocities in **m/s**
- One event per line, whitespace-separated

An example dataset is provided at `datasets/test1k.dat` (1000 events).

When `InputMuPar off`, velocities are sampled from a Maxwell–Boltzmann distribution at 322 K and positions are sampled uniformly.

---

## Output

The output ROOT file contains a TTree named `obe` with one entry per simulated muonium event.

**Per-event scalars (always present):**

| Branch | Type | Description |
|---|---|---|
| `EventID` | `Int_t` | Event index |
| `x`, `y`, `z` | `Double_t` | Initial position [mm] |
| `vx`, `vy`, `vz` | `Double_t` | Initial velocity [m/s] |
| `DoppFreq` | `Double_t` | Peak Doppler frequency [GHz] |
| `PeakIntensity122` | `Double_t` | Peak 122 nm intensity [W/mm²] |
| `PeakIntensity355` | `Double_t` | Peak 355 nm intensity [W/mm²] |
| `Step_n` | `Int_t` | Number of ODE time steps taken |
| `LastRho_gg` | `Double_t` | Ground state population at end of simulation |
| `LastRho_ee` | `Double_t` | Excited state population at end of simulation |
| `LastRho_ion` | `Double_t` | Ionized state population at end of simulation |
| `IfIonized` | `Int_t` | Ionization flag (1 = ionized, 0 = not) |
| `IoniTime` | `Double_t` | MC-sampled ionization time [ns]; -1 if not ionized |

**Per-timestep arrays (present only if enabled via `RootOutput`):**

| Branch | Description |
|---|---|
| `t` | Time [ns] |
| `EField` | Electric field amplitude [V/mm] |
| `RabiFreq` | Rabi frequency [GHz] |
| `Intensity122` | 122 nm laser intensity [W/mm²] |
| `Intensity355` | 355 nm laser intensity [W/mm²] |
| `GammaIon` | Ionization rate γ_ion [GHz] |
| `rho_gg` | Ground state population |
| `rho_ee` | Excited state population |
| `rho_ge_r` | Coherence ρ_ge, real part |
| `rho_ge_i` | Coherence ρ_ge, imaginary part |
| `rho_ion` | Ionized state population |

To inspect the output interactively:

```bash
root -l run/data/g2edmIoni_test.root
# In the ROOT prompt:
new TBrowser   # GUI file/tree browser
```

---

## Example Analysis

An example ROOT macro is provided at `run/ana/example_analysis.C`. It demonstrates how to open the output file, connect all tree branches, and loop over events. A brief ionization summary is printed at the end. The section marked `ADD YOUR ANALYSIS HERE` is where you add your own code.

**Quick start** — run `ioni_test.mac` first, then the analysis:

```bash
# 1. Build
mkdir build && cd build && cmake .. && make -j$(nproc) && cd ..

# 2. Create the output directory and run the test simulation
mkdir -p run/data
cd run
../build/IonizationToyMCSim ioni_test.mac

# 3. Run the example analysis
cd ana
root -l example_analysis.C
```

The macro accepts an optional file path argument if your output file has a different name:

```bash
root -l 'example_analysis.C("../data/my_output.root")'
```

---

## Physics Overview

The simulation evolves a 4-level density matrix under the OBE for a two-photon ionization scheme:

```
|g⟩ ──── 122 nm ────▶ |e⟩ ──── 355 nm ────▶ |ion⟩
         (Rabi Ω)           (rate γ_ion)
```

The equations (in GHz / ns units) are:

```
dρ_gg/dt  =  γ₁·ρ_ee + Im(Ω·ρ_ge)
dρ_ee/dt  = −(γ₁ + γ_ion)·ρ_ee − Im(Ω·ρ_ge)
dρ_ge/dt  = −(γ₂ + iΔ + γ_ion/2)·ρ_ge + i·Ω/2·(ρ_ee − ρ_gg)
dρ_ion/dt =  γ_ion·ρ_ee
```

Key parameters:
- **γ₁ = 0.627 GHz** — Einstein A coefficient of the 1S–2P transition
- **γ₂** — decoherence rate (includes spontaneous emission and laser linewidth)
- **Ω** — Rabi frequency, proportional to the 122 nm electric field amplitude
- **Δ** — detuning = Doppler shift (v·k) + manual detuning
- **γ_ion** — ionization rate from 355 nm, proportional to its intensity (cross-section: 1.26×10⁻¹⁷ cm²)

Integration uses an adaptive Runge–Kutta Cash–Karp (4/5) stepper from Boost `odeint`.

---

## Known Limitations

- **Doppler shift**: computed using only the first 122 nm laser's wave vector. If multiple 122 nm lasers are configured, the Doppler contributions from lasers 2, 3, … are ignored. The 355 nm laser's Doppler effect is always ignored.
- **OBE approximation**: the rotating wave approximation (RWA) is assumed throughout.

---

## Version History

### v5
Refined README with full build instructions, macro reference, output branch documentation, and physics overview.  
Add `run/ana/example_analysis.C`: a working ROOT macro that reads the simulation output, loops the event tree, and prints per-event results and an ionization probability summary.

### v4
Enable multiple 122 nm and 355 nm lasers. Multiple `AddLaser122`/`AddLaser355` commands can be used in a single macro.  
**Important limitation**: the Doppler shift is calculated using only the first laser (OBE solver constraint).  
Add `Intensity122`/`Intensity355` branches to output.  
Enable setting random seed and simulation end time via macro:
```
SetRunTime  10       # simulation duration in ns
RandomSeed  999
```

### v3
Add laser beam angle control via x-y-z Euler angles (yaw, pitch, roll).  
Simulation parameters are now read from macro files instead of being hardcoded. Example: `run/ioni_test.mac`.

### v2
Read muonium position/velocity distribution from an input file.  
Add Monte Carlo sampling to determine ionized muon based on ρ_ion.
