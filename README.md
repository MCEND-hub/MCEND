![License: MIT](https://img.shields.io/github/license/MCEND-hub/MCEND)
![GitHub Workflow Status](https://img.shields.io/github/actions/workflow/status/MCEND-hub/MCEND/CI.yml?branch=main)
![Language](https://img.shields.io/github/languages/top/MCEND-hub/MCEND)
![Version](https://img.shields.io/github/v/tag/MCEND-hub/MCEND-tools)

# MCEND
MCEND (multi-configuration electron-nuclear dynamics) is a quantum dynamics code for simultaneous propagation of electrons and nuclei, written in Fortran90.

## Requirements
MCEND requires Blas, Lapack and FFTW3 libraries; the tools require python3 and numpy, plus additional python libraries. The documentation requires doxygen. You will need a Fortran compiler to compile the source code. Intel Fortran is recommended as the resulting binary is computationally more efficient.

On ubuntu, for an Intel Fortran version of 2023.2.0:
```
sudo apt-get install intel-oneapi-compiler-fortran-2023.2.0 intel-oneapi-mkl-2023.2.0 libfftw3-dev
```

## Compilation
To compile the source code, type
```
make
```
in the MCEND directory.

There are also different options available: Compile with Intel classic `ifort` compiler
```
make classic
```
Compile with gnu Fortran
```
make gfortran
```
All of these can also be used with options for debugging and no optimization (otherwise, `-O3` optimization is used). Just combine the compiler option with the `debug` option, for example for intel classic:
```
make classic_debug
```
Alternatively, `make debug` (for `ifx`) and `make debug_gfortran` (for `gfortran`) are available.
The source code compiles and runs on linux and MacOS, it has not been tested on Windows.

## MCEND basis library

In order to run MCEND, you need to install the [basis library](https://github.com/MCEND-hub/MCEND-library). Either clone the library directly into a directory of your choice, or use it through the provided [git submodule](https://git-scm.com/book/en/v2/Git-Tools-Submodules). In order to check out the files of the submodule, either clone MCEND directly using the `--recursive` flag, or type
```
git submodule update --init --recursive
```
in your terminal within the MCEND folder. You should then find all the basis set files inside `MCEND-library/basis_library/`.

To provide the directory of the basis set files, you can also set a symlink that points to the folder.

For the ease of use, set an `MCEND_BASIS_LIBRARY` environment variable to point to the MCEND basis set library. For example, in bash:
```
export MCEND_BASIS_LIBRARY="/User/username/MCEND-libary/basis_library
```
where you need to adapt the path accordingly.

## Developer's documentation
The source code documentation can be build on your local machine using
```
cd doc
doxygen Doxyfile
```
Open the resulting doc/html/index.html with your browser. 
