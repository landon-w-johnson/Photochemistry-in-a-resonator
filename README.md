# Landon's As-of-yet Untitled RT-TDDFT Photochemistry Tool

### The Goal

This is a real-time time dependent density functional theory (RT-TDDFT) package that is designed as an add-on to the Vienna Ab-initio Simulation Package (VASP). The goal is to incorporate explicit electron-photon coupling into VASP calculations so that molecular dynamics (MD) trajectories can be run with electronic occupations being updated "on-the-fly." Further, this is to be done in a scalable fashion, allowing for the simulation of trajectories on the order of picoseconds.

### The Basic Idea

The general methodology here is to piggy-back off of VASP and use time-dependent perturbation theory (TDPT) to update the electronic occupations in between nuclear steps. This is done by running VASP for a single nuclear time step and using its outputs (namely the **unperturbed** electronic eigenstates and eigenenergies) in TDPT to numerically solve for the change in electronic state occupations due to a perturbing electric field.

### The Current Status

This project is still in a validation stage. The major assumptions at present are

1) The system is inside a Fabry-Perot resonator that is perfectly tuned to the excitation energy between the user-specified bands (even if the molecular resonant frequency is changing). This allows all other electric field frequencies to be ignored as they will die off very quickly in such a resonator.
2) The electric field is treated classically. This is for the sake of debugging and validation against the analytical solution of Rabi oscillations, which do not require quantized photons. The photons will be quantized eventually.
3) The molecule can absorb energy from the E-field, but it cannot release energy back into the E-field. This will likely be fixed at the same time that the quantized treatment of photons is implemented.

This package is able to numerically capture the analytical solution of Rabi oscillations *for frozen nuclei* in multiple molecular systems.

This package runs successfully for molecular systems in which the ions are allowed to move, but there are no longer any analytical solutions available for comparison in this regime. As such, there is no guarantee that the results are trustworthy or realistic for mobile nuclei.

Using bands that are degenerate in energy will cause logical errors in this package, as the character of the degenerate bands will switch (seemingly at random) back and forth during a MD trajectory. This, in turn, causes errors in the calculations of the transition dipole moments, which are critical to the TDPT framework. This will be corrected if I can figure out how to do it.

## How To Use It

First and foremost, *please* contact me at <landon.w.johnson@ndsu.edu> if you want to get this package up and running. It's still a long way from being as user-friendly as I'd like it to be.

As this package requires VASP to do anything in the first place, it is designed with the intent of being easy to learn and use for VASP users, i.e. the user specifies the instructions for this package via tags in the `PHOTCAR` file in an identical fashion to how the user specifies the instructions for VASP via tags in the `INCAR` file.

In order to run this package, the user will need the files necessary for VASP itself:
- `INCAR`
- `POSCAR`
- `POTCAR`
- `KPOINTS` (only gamma point calculations at present)
After optimizing the geometry and assuming that you have the afforementioned files for that optimized geometry, you will also need
- `WAVECAR`
- `firstTimeStep.f`
- `firstTimeStep.exe` (the `.f` file is compiled to this automatically in my scripts)
- `updateWAVECAR.f`
- `updateWAVECAR.exe` (the `.f` file is compiled to this automatically in my scripts)
- `twoLevelSystem.sh` (this is designed for use on one of two specific computer clusters. You will need to ensure that you've appropriately commented/uncommented lines that end in `# CCAST` or `# NERSC`. You will likely need to modify such lines appropriately for any other computer cluster that you run this on.)
- a job submission script for the computer cluster you're using (`job.pbs` (PBS script for CCAST) and `submit-NERSC.sh` (slurm script for NERSC))

Ensure that your job submission script runs `twoLevelSystem.sh`! This is very important as this script controls the flow of the overarching package. The package will save important data (some of which it will re-use) to `dataFile.txt` in your working directory. `twoLevelSystem.sh` is designed to direct *all* standard outputs to `outputFile.txt`. This allows the user to check on the status of the job mid-run and debug as needed while allowing the job to continue if desired.

### Submitting Batch Jobs

`makeTrajectoryGrid.sh` is a script that will take the contents its working directory and copy them into subdirectories that it makes for different initial conditions. This script has cluster-specific lines (e.g. lines that end in `# CCAST` or `# NERSC`) that will need to be (un)commented or modified appropriately for the computer cluster being used. Also note that the line that writes the value for `omega` (the frequency of E-field amplitude oscillations in atomic units) may need to be modified. Ensure that there are the correct number of preceeding 0s for the value of `omega` that it writes. The arithmetic that calculates the value for `omega` is integer-based and does not track the leading 0s after the decimal point.

### Other Auxiliary Scipts

The other scripts in `AuxiliaryScripts` are mostly for plotting data. The MATLAB scripts will likely work the best. Some notable exceptions are:
- `wavecar2text.f`: reads `WAVECAR` and converts it from binary to text. This is for debugging purposes. Everything in `WAVECAR` is a floating-point real number except for the plane wave expansion coefficients that VASP uses to represent the bands, which are saved as complex floating-point numbers. So depending on what you want to read from `WAVECAR`, you'll need to swap which data type the script uses.
- `rerun_Debug.sh`: I use this for debugging changes to `firstTimeStep.f` and `updateWAVECAR.f`. It just recompiles and reruns the relevant script quickly so you can see what happens.
- `prepTemp&Freq.sh`: gathers time-dependent internuclear distance (specifically for H2+), temperature readings of VASP, and molecular resonant frequency.
- `patchLog.txt`: This is where I keep track of the things I've been doing with this package. I might change my archaic ways now that I'm learning how to `git` gud.

#### Miscellaneous

Things that show up in file names:
- H2+: It's for an H2+ molecule
- SiC: It's for silicon carbide
- Rabi: This file is intended for recreation of Rabi oscillations with *frozen nuclei*
- Phot0K: This file is intended for mobile nuclei, but starts with no thermal motion.
- *N*x*N*: It's a *N*x*N* supercell.
- DiVac: It has a divacancy defect (V_Si V_C)
- OC: It's off-center (low-dimensional materials directly in the center of the simulation cell cause problems. I believe this is because **every** plane wave that is used to construct the electronic bands has either a node or an anti-node right where we're looking, which causes some instability.)
- Opt: The geometry has already been optimized
