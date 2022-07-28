# MCEND-v2.0

Issues:
- writeout should be moved to specified module and unified
- not all files should be written at all times, there should be more fine-grained options as each additional write slows down the program
- do not reintroduce variables that are already set in parent modules
- don't leave loops in there with commented out printout
- use an IDE for coding, this will provide you with a linter
- use standard formatting, ie only "" quotes, and spacing, so that one can search through the code easily instead of having to try all different options of writing the phrase
- if you want to write to different files for spinorbital and standard propagation, at least do not create all files for all runs, then only the ones that you are going to write into
- if you comment in the code, keep the comments updated. Check also the doxygen documentation, this will put all the comments together and you can see which ones are outdated. ("Any comment can become a lie in the future.")
## Requirements
MCEND requires Blas, Lapack and FFTW3 libraries; the tools require python3 and numpy, plus additional python libraries. The documentation requires doxygen.

## Documentation
The source code documentation can be build on your local machine using
```
cd doc
doxygen Doxyfile
```
Open the resulting doc/html/index.html with your browser. Please extend the documentation to any new functionalities, and fill in missing ones.
