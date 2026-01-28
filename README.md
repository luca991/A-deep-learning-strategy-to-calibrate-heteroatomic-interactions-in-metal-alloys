# A deep learning strategy to calibrate heteroatomic interactions in metal alloys

In this repository all the programs and scripts used to perform the calculation reported in "A deep learning strategy to calibrate heteroatomic interactions in metal alloys" are present.

## Monte Carlo simulations
The *MC* folder contains the code for Monte Carlo simulations and the files needed to run it. To compile, simply run the command `make all` and the executable *Unit_cell* will be generated. When run, it will perform a Monte Carlo simulation. To run *Unit_cell*, three other files must be present in the same folder:
  - *leggi.in*: defines the MC simulation setup. Not all lines in this file are used by this programme (lines 3 and 7 are not used). *leggi_example.in* is an example of how to fill *leggi.in*
  - *potential.in*: in this file, the potential parameters are defined. The parameters related to the first element (El1) are named with '_El1El1', similarly for the second element (EL2) with '_El2El2', while the parameters related to the mixed interactions are named with '_El1El2'. In line 3, *lattice_type_El1* and *lattice_type_El2* are replaced with the lattice of the two pure metals (hcp, bcc, fcc) and the values of the relative parameters are replaced. *potential_example.in* is an example of how to fill *potential.in*
  - *cell_in_#.in*: are the files containing the initial configurations to be simulated. There is one file for each configuration. The # must be replaced with the composition of the system in terms of the first element, as defined in line 6 of *leggi.in*. These files have to be in the `.xyz` format, as follow:

    > N  
    > Lattice=" L_xx L_xy L_xz L_yx L_yy L_yz L_zx L_zy L_zz " Properties=species:S:1:pos:R:3  
    > El	x y z  
    > El	x y z
    > El	x y z ... (N times)

With N the number of atoms, L_** the box size and the list of atomic positions with their atomic species

The *Unit_cell* programme can then be started with `./Unit_cell #`, where the # in this case is an index for enumerating the simulations.

After the simulations are complete, if you want to calculate the mixing enthalpy, you can use the Python code *calc_mix_enthalpy.py*, which, starting from the results of *Unit_cell*, calculates the mixing enthalpy for each simulated composition. In this case, the two pure systems must also be simulated.

## Dissolution energy
The *Impurity_energy* folder contains the code for energy of dissolution calculation and the files needed to run it. Run the command `make all` to compile the programme and the executable *impure* will be generated.
For this programme too, files *leggi.in* and *potential.in*, as defined for MC simulations, must be located in the same folder where *impure* is executed. *impure* calculates the two dissolution energies and saves them in a file.

## Dataset
In *Dataset* folder the dataset used to train the Neural Networks is present with the name *data_ML_article.dat*. This file contains, for each fictitious system:

  - A_00, xi_00, p_00, q_00, r0_00: potential parameters related to the first element (_El1El1 parameters in *leggi.in*)
  - A_11, xi_11, p_11, q_11, r0_11: potential parameters related to the second element (_El2El2 parameters in *leggi.in*)
  - A_01, xi_01, p_01, q_01, r0_01: potential parameters related to the mixed interaction (_El1El2 parameters in *leggi.in*)
  - LS_0, LS_1: lattice structures of El1 and El2 respectively
  - T: simulation temperature
  - H_mix_#: mixing enthalpy at differen composition of the alloy (# = 5, 50, 95)
  - latt_par_#: lattice parameter at differen composition of the alloy (# = 5, 50, 95) 

, the potential parameters, mixing enthalpies, lattice parameters and dissolution energies.

## Neural Network training
*NN_training* contains the Python scripts for training the neural networks. Each programme trains the neural network of a physical quantity. 
To run these programmes, you need to have the dataset (located in folder *Dataset* and named *data_ML_article.dat*) and the file with the extreme values for each quantity (*min_max_values.dat*) in the same folder.
Folder *NN_training* also contains the NNs trained by us that were used in the next step. These are identified by the name of the quantity they predict and have the extension *.h5*.

## Minimization algorithms
The *Minimization* folder contains the algorithms that follow local minimisation (see section 3.3 of the article). There are two programmes, one for systems with positive mixing enthalpy (*minimization_positive_H_mix.py*) and another for systems with negative mixing enthalpy (*minimization_negative_H_mix.py*). 
In order to run, they must be located in the same folder as the neural networks and the *min_max_values.dat* file to be used, and the *Pure_elements.dat* and *Byn_prop_resc.dat* files containing the data for pure metals and alloys respectively must also be present. 
To run the programmes, launch the `python minimization_positive_H_mix.py El1El2` (or `python minimization_negative_H_mix.py El1El2`) command, where El1 and El2 are the symbols of the metals that form the alloy and must be in the order used in the *Byn_prop_resc.dat* file.
