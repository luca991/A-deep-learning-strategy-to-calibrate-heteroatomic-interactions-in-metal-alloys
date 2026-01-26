# A deep learning strategy to calibrate heteroatomic interactions in metal alloys

In this repository all the programs and scripts used to perform the calculation reported in "A deep learning strategy to calibrate heteroatomic interactions in metal alloys" are present.

## Monte Carlo simulations
The *MC* folder contains the code for Monte Carlo simulations and the files needed to run it. To compile, simply run the command `make all` and the executable *Unit_cell* will be generated. When run, it will perform a Monte Carlo simulation. To run *Unit_cell*, three other files must be present in the same folder:
  - *leggi.in*: defines the MC simulation setup. Not all lines in this file are used by this programme (lines 3 and 7 are not used)
  - *potential.in*: in this file, the potential parameters are defined. In line 3, *lattice_type_El1* and *lattice_type_El2* are replaced with the lattice of the two pure metals (hcp, bcc, fcc) and the values of the relative parameters are replaced.
  - *cell_in_#.in*: are the files containing the initial configurations to be simulated. There is one file for each configuration. The # must be replaced with the composition of the system in terms of the first element, as defined in line 6 of *leggi.in*.

The *Unit_cell* programme can then be started with `./Unit_cell #`, where the # in this case is an index for enumerating the simulations.

After the simulations are complete, if you want to calculate the mixing enthalpy, you can use the Python code *calc_mix_enthalpy.py*, which, starting from the results of *Unit_cell*, calculates the mixing enthalpy for each simulated composition. In this case, the two pure systems must also be simulated.

## Dissolution energy
The *Impurity_energy* folder contains the code for energy of dissolution calculation and the files needed to run it. Run the command `make all` to compile the programme and the executable *impure* will be generated.
For this programme too, files *leggi.in* and *potential.in*, as defined for MC simulations, must be located in the same folder where *impure* is executed. *impure* calculates the two dissolution energies and saves them in a file.

## Dataset
In *Dataset* folder the dataset used to train the Neural Networks is present with the name *data_ML_article.dat*. This file contains, for each simulated system, the potential parameters, mixing enthalpies, lattice parameters and dissolution energies.

## Neural Network training
*NN_training* contains the Python scripts for training the neural networks. Each programme trains the neural network of a physical quantity. 
To run these programmes, you need to have the dataset (located in folder *Dataset* and named *data_ML_article.dat*) and the file with the extreme values for each quantity (*min_max_values.dat*) in the same folder.
