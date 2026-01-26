In this repository all the programs and scripts used to perform the calculation reported in "A deep learning strategy to calibrate heteroatomic interactions in metal alloys" are present.

The ‘MC’ folder contains the code for MC simulations and the files needed to run it. To compile, simply run the command “make all” and the executable “Unit_cell” will be generated. When run, it will perform a Monte Carlo simulation. To run “Unit_cell”, three other files must be present in the same folder:
  - 'leggi.in': defines the MC simulation setup. Not all lines in this file are used by this programme (lines 3 and 7 are not used)
  - 'potential.in': in this file, the potential parameters are defined. In line 3, “lattice_type_El1” and “lattice_type_El2” are replaced with the lattice of the two pure metals (hcp, bcc, fcc) and the values of the relative parameters are replaced.
  - 'cell_in_#.in': are the files containing the initial configurations to be simulated. There is one file for each configuration. The # must be replaced with the composition of the system in terms of the first element, as defined in line 6 of “leggi.in”.

The “Unit_cell” programme can then be started with “./Unit_cell #”, where the # in this case is an index for enumerating the simulations.

