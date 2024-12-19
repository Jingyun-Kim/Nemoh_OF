## 2024. 12. 18.
Add TestCases/VolturnUS-S/

==========================================================================================
# Nemoh_OF
## Discription
This repository contains modified source code from Nemoh, including adjustments to its Compile and QTF Module, to obtain WAMIT format Hydrodata(.hst, .1, .3) used in the HydroDyn Module of OpenFAST

## Reference
Nemoh gitlab - <https://gitlab.com/lheea/Nemoh>
	
==========================================================================================
## Original Nemoh README.md

NEMOH is a Boundary Element Methods (BEM) code dedicated to the computation of first order wave loads on offshore structures (added mass, radiation damping, diffraction forces). It has been developed by researchers at Ecole Centrale de Nantes for more than 30 years. It is still used in many of our research projects. Typical uses include the estimation of dynamic response of floating structures or the performance assessment of wave energy converters.

Copyright © 2022 Ecole Centrale de Nantes

This program has been developed at [Ecole Centrale de Nantes](http://www.ec-nantes.fr) by
G. Delhommeau, A. Babarit, J. Singh, P. Guével, J.C. Daubisse, R. Kurnia.

## Pre-built binaries and documentation

The full documentation is available [here](https://lheea.gitlab.io/Nemoh). Please do not hesitate to report any errors or suggest improvements!

Pre-built binaries (executables) of all the Nemoh programs for Linux and Windows (x64) can be found in the `Releases` section of the project for each version (starting from v3.0.0).

Alternatively, you can compile Nemoh by yourself using the following instructions (of by following [the documentation](https://lheea.gitlab.io/Nemoh)).

## Compilation

Compile all Nemoh executables using CMake (from the root of the repository):

```shell
cmake -S. -Bbuild
cmake --build build
```

The resulting executables are in the `bin` directory.

To compile only one of the executables, use the `--target` option of CMake.
The available targets for Nemoh are:
- `preProc`
- `solver`
- `postProc`
- `mesh`
- `hydrosCal`

The available targets (executables) for Nemoh QTF are:
- `QTFpreProc`
- `QTFsolver`
- `QTFpostProc`

The choice of the compiler is left to CMake, but can be overridden by setting the `CMAKE_Fortran_COMPILER` at the configuration step, e.g.:

```shell
cmake -S. -Bbuild -DCMAKE_Fortran_COMPILER=gfortran
```

Compilation was tested with GNU's `gfortran` and Intel's `ifort`.

## Testing

After building, the fast tests can be run from the `build` directory:

```shell
ctest -V -j <N_concurrent>
```

Where `<N_concurrent>` is the number of simultaneous workers (processes).

The tests can be restricted using their labels and the `-L` option  of `ctest`:

```shell
ctest -V -j <N_concurrent> -L <label>
```

Where label is one of the following:
- `NEMOH1`: only the non-QTF test cases
- `PREPROC`: only the pre-processing operations
- `SOLVER`: only the solving operations (depend on the pre-processing tests)
- `POSTPROC`: only the post-processing operations (depend on the pre-processing and solving tests)
- `NEMOH2`: only the QTF test cases
- `QTF`: only the computation of the QTF (depend on the prior non-QTF Nemoh computation)

Tests with unsatisfied requirements will fail.

In order to run the full test cases (matching reference data), Makefiles are provided in the `TestCases` directory and in each of the test cases' directories.

## Bugs, improvements and questions

Please use the GitLab Issues to report any bugs you might have encountered, to suggests improvements to the code and documentation and to ask questions.

External contributions are welcome, but we suggest that you contact us beforehand to avoid conflicts with ongoing developments.
