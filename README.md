# IonizationToyMCSim

## v2
Read muonium distribution from input file.  
Add Monte Carlo to sample ionized muon based on $\rho_{ion}$.
Parameters in `main.cpp` are final designed laser values.

## v3
Add laser angles by x-y-z order euler angle  
Now input parameters from macro files. Example: `run/ioni_test.mac`

Execute the program:
```
IonizationToyMCSim ioni_test.mac
```

## v4
Enable multiple 122nm and 355nm laser. You can simply add laser in the macro.  
**IMPORTANT LIMITATION**: The doppler shift is calculated only use the first laser! (Due to OBE limitation)   
Add Intensity branch, both for 122 and 355 nm.  
Enable set random seed and simulation end time through macro:
```
# Set simulation end time (in ns)
SetRunTime 10
```
```
# Set random seed
RandomSeed 999
```