# Collectice-Chemotaxis-CPM

Source code for simulations of collective chemotaxis using the cellular Potts Model (CPM) as the computational implementation.

## Compilation

To compile the code you need a FORTRAN compiler like [gfrotran](https://gcc.gnu.org/fortran/). The following command will compile the code and output the executable `a.out` used to run simulations.

`gfortran -c mod1* mod2* mod3*;rm *.o;gfortran cpm_main.f90 mod*`

## Compilation Details

Files that are modules need to be compiled first by themselves before being used in compilation of the main program file. All module files start with the string `mod` in their filenames.
Module files can be individually compiled using the following command:
`$ gfortran -c my_module.f90`

Once all modules are compiled, proceed to compile the main program using the following command.
`$ gfortran -o a.out main.f90 my_module_1.f90 my_module_2.f90 `

For debugging, the [following flags are helpful](http://stackoverflow.com/questions/3676322/what-flags-do-you-set-for-your-gfortran-debugger-compiler-to-catch-faulty-code).
