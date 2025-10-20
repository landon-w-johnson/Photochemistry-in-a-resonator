# Constraints

1. If the user can use VASP (even if it's just barely), they can use this. Therefore, this program shall only require a `PHTOCAR` file (or whatever "more appropriate" name) that is equivalent to the `INCAR` file for VASP.

1. This program shall not require anything that is not explicitly required to run VASP. Therefore it should be written in compliance with F2008. [Here are VASP's technical requirements for installation.](https://www.vasp.at/wiki/Installing_VASP.6.X.X)

1. The goal is to make this a high-throughput program, so computational efficiency matters more than prettiness, unless it makes things difficult for the user.




# Pseudocode (assuming a single Fortran file)

```
parse PHOTCAR
instantiate necessary data structures and output files
run first timestep(?) (I previously had issues trying to use the updateWAVECAR scheme for the first timestep since it started with 0 occupation for the excited electronic state. If we can find a way to use the general updateWAVECAR scheme without it slowing things down, then let's do it. Updating the occupancies comes before VASP since it is assumed that the user already has an optimized geometry or a properly heated system.)
    read WAVECAR info
    perform perturbation theory
    write new occupancies into WAVECAR
    output relevant data into output files for plotting and post-processing
loop
    run VASP
    run updateWAVECAR (or whatever we call it)
        read WAVECAR info
        perform perturbation theory
        write new occupancies into WAVECAR
        output relevant data into output files for plotting and post-processing
```