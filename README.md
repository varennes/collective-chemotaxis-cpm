# Collectice-Chemotaxis-CPM

Source code for simulations of collective chemotaxis using the cellular Potts Model (CPM) as the computational implementation. In these simulations biological cells sense a chemical gradient in the environment using intercellular communication and attempt to migrate in the direction of increasing chemical concentration based on their measurements.

The resulting migratory behavior of cell clusters of various sizes has been studied and the results are described in [**this paper**](http://arxiv.org/abs/1605.00712).

## Compilation

To compile the code you need a FORTRAN compiler like [gfrotran](https://gcc.gnu.org/fortran/). The following command will compile the code and output the executable `a.out` used to run simulations.

`$ gfortran -c mod1* mod2* mod3*;rm *.o;gfortran cpm_main.f90 mod*`

## Content

The file `cpm_main.f90` is the main program file which simulates the time evolution of the cell cluster and outputs simulation results. The various module files (all module filnames start with the string `mod`) contain functions and subroutines which are called by the main program file.

Functions and subroutines within the module files are annotated in order to explain their use. Module files are organized in the following manner.
- `mod1util.f90`: Initialization functions and commonly used functions.
- `mod2sense.f90`: All functions relevant to creating the chemical concentration profile and measuring the chemical gradient.
- `mod3goal.f90`: All functions pertaining to the evaluation of the energy and bias term for the CPM.
- `mod3polar.f90`: All function relevant to assigning individual cell polarization vectors.
- `mod3sc.f90`: Functions pertaining to the CPM implementation.
- `mod3wrt.f90`: Functions that print out simulation information to various output files.

## Further Compilation Details

Files that are modules need to be compiled first by themselves before being used in compilation of the main program file. All module files start with the string `mod` in their filenames.
Module files can be individually compiled using the following command:

`$ gfortran -c my_module.f90`

Once all modules are compiled, proceed to compile the main program using the following command.

`$ gfortran -o a.out main.f90 my_module_1.f90 my_module_2.f90 `

For debugging, the [following flags are helpful](http://stackoverflow.com/questions/3676322/what-flags-do-you-set-for-your-gfortran-debugger-compiler-to-catch-faulty-code).
