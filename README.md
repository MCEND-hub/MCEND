![License: MIT](https://img.shields.io/github/license/iulusoy/MCEND)
![GitHub Workflow Status](https://img.shields.io/github/actions/workflow/status/iulusoy/MCEND/CI.yml?branch=main)
![Language](https://img.shields.io/github/languages/top/iulusoy/MCEND)
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
Alternatively, `make debug` (for `ifx`) and `make debug_gfortran (for `gfortran`) are available.
The source code compiles and runs on linux and MacOS, it has not been tested on Windows.

## Developer's documentation
The source code documentation can be build on your local machine using
```
cd doc
doxygen Doxyfile
```
Open the resulting doc/html/index.html with your browser. 
