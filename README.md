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

### Measured Transverse Profile (optional)

By default each laser's transverse intensity is the analytic Gaussian `exp(-2(x²/σ_x² + y²/σ_y²))`. A laser can instead use a measured 2D transverse profile:

```
SetLaser122Profile <path/to/profile.root>   # applies to the most-recently-added 122 nm laser
SetLaser355Profile <path/to/profile.root>   # applies to the most-recently-added 355 nm laser
```

- The file must contain a `TH2D` named `h_profile`: a transverse **density** in mm⁻² (unit integral, `Σ bin·Δx·Δy = 1`) on laser-frame axes — x = σ_x direction (= target Y, vertical), y = σ_y direction (= target Z) — centred on its own intensity centroid.
- Substitution: `I_spatial = E / (√(2π)·τ[s]) · h_profile(x, y)` (bilinear `TH2::Interpolate`); outside the histogram the intensity is exactly 0. For 122 nm the E-field amplitude is `√(2·η·I_spatial)`, η = 376.7303134 Ω, so field and intensity stay consistent.
- Path set ⇒ profile; path absent ⇒ analytic Gaussian (the default, unchanged).
- `σ_x`/`σ_y` on the `AddLaser*` line are **ignored** for a profiled laser (a warning is printed). Everything else — energy, FWHM, `t_peak`, offsets, Euler angles, detuning — keeps its meaning, so the offsets slide the sampling point across the profile.
- To align a 355 nm profile relative to a 122 nm one, put the measured relative shift `centroid_camera_mm(122) − centroid_camera_mm(355)` into the 355 laser's `offset_y`/`offset_z` (helper: `FromClaude/realprofile_validation/relative_offset.py`).

Example macros: `run/example_realprofile.mac` (commented walk-through with the measured profiles and the measured 122↔355 offset), `run/ioni_test_realprofile.mac` (measured profiles, same beam as `ioni_test.mac`) and `run/ioni_test_gaussianprofile.mac` (ideal Gaussian profile files, same beam as `ioni_test.mac`). Profile files: `datasets/122_profile_test.root` and `datasets/355_profile_test.root` (measured); `datasets/122_profile_gaussian.root` and `datasets/355_profile_gaussian.root` (synthetic, generated by `FromClaude/realprofile_validation/make_gaussian_profile.py`).

### Laser Parameter Jitter (shot-to-shot fluctuation)

Any laser parameter can optionally be resampled every event from a Gaussian distribution instead of staying fixed, to model realistic shot-to-shot laser fluctuation:

```
AddLaser122Sigma <same 13 params as AddLaser122, as Gaussian σ>
AddLaser355Sigma <same 12 params as AddLaser355, as Gaussian σ>

LaserJitter on | off    # default off; when on, resample every laser parameter each event
                        # from N(macro value, σ), using the RandomSeed RNG
```

`AddLaser122Sigma`/`AddLaser355Sigma` set the standard deviation for the *most-recently-added* laser of that wavelength, in the same order/units as `AddLaser122`/`AddLaser355`. A parameter left at σ=0 (the default) stays exactly fixed at its macro value every event, even when `LaserJitter` is on. `LaserJitter` defaults to `off`, so existing macros that never call it behave exactly as before (bit-identical output).

### Simulation Control

```
SetRunTime      <duration[ns]>      # Total simulation time window
RandomSeed      <integer>           # RNG seed for reproducibility (also drives LaserJitter)
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

RootOutput  LaserPars122 on | off   # Record the realized (nominal or jittered) 122nm laser parameters per event
RootOutput  LaserPars355 on | off   # Record the realized (nominal or jittered) 355nm laser parameters per event

OutputFile  <path/to/output.root>
```

`LaserPars122`/`LaserPars355` default to `off`. Each supports only a single laser of that wavelength — the simulation throws if more than one `AddLaser122`/`AddLaser355` is configured while the corresponding output is on.

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

A version with shot-to-shot laser jitter enabled (see `run/ioni_test_jitter.mac` for the full example):

```
AddLaser122       13.5e-6  2  5  80  4  1  0  0  2  0  0  0  0
AddLaser122Sigma  1.35e-6  0  0.2  0  0.2  0.05  0  0  0  0  0  0  0   # 10% energy, 5% beam-size, 0.2ns timing jitter
LaserJitter       on

RootOutput  LaserPars122  on   # record the realized energy/sigma_x/sigma_y/... per event
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

An example dataset is provided at `datasets/test1k.dat` (2700 events, scattered). For a smooth spatial map, use the regular grid `datasets/MuGrid_z0-8mm_y-40to40mm_10k.dat` (100 × 100 points over z ∈ [0, 8] mm, y ∈ [−40, 40] mm, x = 0, v = 0), generated by `FromClaude/realprofile_validation/make_intensity_check_grid.py`, and used by `run/ioni_grid_realprofile.mac` / `run/ioni_grid_gaussianprofile.mac`.

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
| `PeakIntensity122` | `Double_t` | Peak 122 nm intensity over the ODE steps (running max) [W/cm²] |
| `PeakIntensity355` | `Double_t` | Peak 355 nm intensity over the ODE steps (running max) [W/cm²] |
| `Step_n` | `Int_t` | Number of ODE time steps taken |
| `LastRho_gg` | `Double_t` | Ground state population at end of simulation |
| `LastRho_ee` | `Double_t` | Excited state population at end of simulation |
| `LastRho_ion` | `Double_t` | Ionized state population at end of simulation |
| `IfIonized` | `Int_t` | Ionization flag (1 = ionized, 0 = not) |
| `IoniTime` | `Double_t` | MC-sampled ionization time [ns]; -1 if not ionized |

**Per-event laser parameter snapshot (present only if `RootOutput LaserPars122`/`LaserPars355 on`):**

The realized value of each laser parameter for that event — equal to the macro's nominal value every event unless `LaserJitter on` is set, in which case it reflects that event's Gaussian-sampled draw. `Laser122_*` has 13 branches (adds `Laser122_Detuning`), `Laser355_*` has the same 12 without detuning:

| Branch | Type | Description |
|---|---|---|
| `Laser122_Energy` / `Laser355_Energy` | `Double_t` | Pulse energy [J] |
| `Laser122_Linewidth` / `Laser355_Linewidth` | `Double_t` | Linewidth [GHz] |
| `Laser122_PeakTime` / `Laser355_PeakTime` | `Double_t` | Peak time [ns] |
| `Laser122_SigmaX` / `Laser355_SigmaX` | `Double_t` | Beam radius σx [mm] |
| `Laser122_SigmaY` / `Laser355_SigmaY` | `Double_t` | Beam radius σy [mm] |
| `Laser122_Tau` / `Laser355_Tau` | `Double_t` | Pulse time constant (0.4247×FWHM) [ns] |
| `Laser122_OffsetX/Y/Z` / `Laser355_OffsetX/Y/Z` | `Double_t` | Beam center offset [mm] |
| `Laser122_Yaw/Pitch/Roll` / `Laser355_Yaw/Pitch/Roll` | `Double_t` | Beam orientation Euler angles [deg] |
| `Laser122_Detuning` | `Double_t` | Frequency detuning [GHz] (122 nm only) |

**Per-timestep arrays (present only if enabled via `RootOutput`):**

| Branch | Description |
|---|---|
| `t` | Time [ns] |
| `EField` | Electric field amplitude [V/mm] |
| `RabiFreq` | Rabi frequency [GHz] |
| `Intensity122` | 122 nm laser intensity [W/cm²] |
| `Intensity355` | 355 nm laser intensity [W/cm²] |
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

### Peak Intensity Map: `run/ana/draw_peakI_map.C`

Plots the per-event peak 122 nm and 355 nm intensity (`PeakIntensity122`/`PeakIntensity355`) as a 2D map for both wavelengths side by side: horizontal = target Z, vertical = target Y, colour = the raw **linear** peak intensity in W/cm². It is meant for a run over the dense 100 × 100 grid (`datasets/MuGrid_z0-8mm_y-40to40mm_10k.dat`), which gives a smooth map; works for any profile (measured or synthetic) run over that grid.

```bash
cd run/ana
root -l draw_peakI_map.C                                              # default: ../data/g2edmIoni_grid_realprofile.root
root -l 'draw_peakI_map.C("../data/g2edmIoni_grid_gaussianprofile.root")'
```

Writes `plots/<input stem>_peakI_map.png`. Requires the `obe` tree of the run.

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
- **`RootOutput LaserPars122`/`LaserPars355`**: only supports a single laser of that wavelength; the simulation throws if more than one `AddLaser122`/`AddLaser355` is configured while that output is enabled.
- **Measured profiles** (`SetLaser*Profile`): the transverse shape is static (no per-event shape change); `σ_x`/`σ_y` of a profiled laser are not used; the intensity is exactly zero outside the profile's histogram range (a hard cutoff, not a smooth tail); the 122 nm E-field is derived from the intensity, so its absolute scale follows the profile normalization and the pulse energy.

---

## Version History

### v7
Add a measured transverse laser profile per wavelength, replacing the analytic Gaussian when set:
```
SetLaser122Profile <file.root>   # TH2D "h_profile", unit-integral density in mm⁻²
SetLaser355Profile <file.root>
```
See [Measured Transverse Profile](#measured-transverse-profile-optional). Example macros: `run/ioni_test_realprofile.mac`, `run/ioni_test_gaussianprofile.mac`, and the dense-grid versions `run/ioni_grid_*profile.mac`. Analysis: `run/ana/draw_peakI_map.C`.
Also corrected the documented units of `PeakIntensity122/355` and `Intensity122/355` to W/cm² (what the code has always produced).
Without any `SetLaser*Profile` line, output is unchanged from v6.

### v6
Add per-event Gaussian jitter of laser parameters, to model shot-to-shot laser fluctuation:
```
AddLaser122Sigma  <same 13 params as AddLaser122, as Gaussian σ>
AddLaser355Sigma  <same 12 params as AddLaser355, as Gaussian σ>
LaserJitter       on | off    # default off
```
Add `RootOutput LaserPars122`/`LaserPars355 on | off` to record each event's realized laser parameters (`Laser122_Energy`, `Laser122_SigmaX`, ... — see [Output](#output)). Both default off; only a single laser per wavelength is supported.
Example macro: `run/ioni_test_jitter.mac`.
Feature is fully opt-in: with `LaserJitter` left off (the default), output is bit-identical to v5.

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
