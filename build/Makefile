FC = ifx
FC_classic = ifort
MKLLIB := ${MKLLIB}
# FFTWLIB := ${FFTWLIB}
ifeq ($(MKLLIB),)
  MKLLIB = ${MKLROOT}/lib/intel64
  $(info    MKLLIB taken from MKLROOT $(MKLLIB))
endif

UNAME_S := $(shell uname -s)
# Conditional flags for macOS
ifeq ($(UNAME_S),Darwin)
  # Add FFTW3 library path
  FFTW_LIB_PATH := $(shell brew --prefix fftw)/lib
  FFTWLIB += -L$(FFTW_LIB_PATH) -lfftw3
  $(info    FFTWLIB set to $(FFTWLIB))
endif

ifeq ($(FFTWLIB),)
  FFTWLIB = -lfftw3
  $(info    FFTWLIB not set, using lfftw3)
endif

FFLAGS = -O3 -fopenmp -c -traceback -pg -integer-size 64 -real-size 64 -extend-source -fpe:0 -fpconstant
LINK_FLAGS = -qopenmp -L${MKLLIB} -traceback ${FFTWLIB} -liomp5 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core
FFLAGS_classic = -fopenmp -c -fast -traceback -pg -i8 -r8 -extend-source -fpe0 -fpconstant
LINK_FLAGS_classic = -qopenmp -L${MKLLIB} -traceback ${FFTWLIB} -liomp5 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core

BINDIR=../bin
SRCDIR=../src
OBJ = read-input.o generic_data.o akima_interpol.o utils.o combinations.o init.o analysis.o propa.o mcend-2.1.o

all: $(BINDIR)/mcend

$(BINDIR)/mcend: $(OBJ)
	$(FC) -o $(BINDIR)/mcend $(OBJ) $(LINK_FLAGS)

debug: FFLAGS = -O0 -fopenmp -c -traceback -pg -integer-size 64 -real-size 64 -extend-source -fpe:0 -fpconstant -warn all # -check all -fp-stack-check -fsanitize=memory
# debug: LINK_FLAGS += -fsanitize=memory
debug: $(BINDIR)/mcend

classic: FC = $(FC_classic)
classic: FFLAGS = -O3 $(FFLAGS_classic)
classic: LINK_FLAGS = $(LINK_FLAGS_classic)
classic: $(BINDIR)/mcend

debug_classic: FC = $(FC_classic)
debug_classic: FFLAGS = $(FFLAGS_classic) -O0 -i8 -warn all -check all -traceback -fp-stack-check -g
debug_classic: LINK_FLAGS = $(LINK_FLAGS_classic)
debug_classic: $(BINDIR)/mcend

gfort: FC = gfortran
gfort: FFLAGS = -O3 -c -g
gfort: LINK_FLAGS = -llapack -lblas ${FFTWLIB} -fopenmp -g
gfort: $(BINDIR)/mcend

debug_gfort: FC = gfortran
debug_gfort: LINK_FLAGS = -llapack -lblas ${FFTWLIB} -fopenmp -g
#debug_gfort: FFLAGS += -O0 -c -g -p -ffree-line-length-0 -fcheck=all -ffpe-trap=zero,overflow,underflow -finit-real=nan -Wline-truncation -Wall -fbacktrace -fbounds-check
debug_gfort: FFLAGS += -O0 -c -g -ffree-line-length-0 -fcheck=all -Wline-truncation -Wall -fbacktrace -fbounds-check
debug_gfort: $(BINDIR)/mcend

clean:
	rm -f *.o *.mod *.inc *.a

%.o:	$(SRCDIR)/%.f90
	$(FC) $(FFLAGS)  $<

