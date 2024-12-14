# Instructions for computing electron-phonon dynamics using rt-BTE with Perturbo

## Software 
Quantum Espresso can be found [here](https://www.quantum-espresso.org).

For SUNDIALS installation, refer to SUNDIALS [tutorials](https://computing.llnl.gov/projects/sundials).

For TDEP installation, refer to [TDEP] (https://github.com/tdep-developers).

After installing Quantum Espresso 7.2, PERTURBO can be installed inside the folder. The code is in the folder ../code
Linking SUNDIALS and TDEP in this version of the code is used by specifying compiler flag in the make.sys file.

Refer to PERTURBO [tutorial](https://perturbo-code.github.io) for Perturbo download, installation, and usage with TDEP.

## Usage
First run DFT and DFPT calculations using QE.

Then use TDEP to compute 2nd- and 3rd-order force constants, copy _infile.forceconstant_ and _infile.forceconstant_thirdorder_ and _infile.ucposcar_ into folder `qe2pert`. Create input _qe2pert.in_ with `tdep = .true.`, and run `qe2pert.x` to generate the _prefix_epr.h5_ file for PERTURBO calculations.

For running dynamics of electronic population, set `dyna_phph = .false.`, and refer to [ultrafast dynamics](https://perturbo-code.github.io/mydoc_dynamics.html) tutorials for PERTURBO. 

To compute dynamic phonon populations with ph-ph interaction, set `dyna_phph = .true.`, and corresponding parameters for smearing and cutoff but with the suffix `_ph`.

`divide_grid = .true.` allows a different set of momentum grids for electrons and phonons (electon grid size is twice denser)

For the solver interface, one can set the following parameters:
`solver = 'rk4'`, `'euler'` or `'sundials'`

For `solver = 'sundials'`, one can set `sun_method` to be `mri` (MRI) or `ark` (representing ERK). 

Set `hs` for slow time step in fs in the MRI mode, `retol` for relative tolerance, and `abstol_c` for absolute tolerance for all carrier populations, and `abstol_ph` for that of phonons.

After dynamics is computed. One will arrive at _prefix_cynda.h5_ file, containing all the electron and phonon populations as a function of time. Please refer to the online tutorials on the PERTURBO website for details. The post-processing using `calc_mode = dynamics-pp` will genearte _prefix_popu.h5_, which contains the populations as a function of energy.

