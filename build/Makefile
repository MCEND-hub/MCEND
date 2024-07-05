FC = ifx
FC_classic = ifort
FFLAGS = -O3 -fopenmp -c -traceback -pg -integer-size 64 -real-size 64 -extend-source -fpe:0 -fpconstant 
LINK_FLAGS = -L${MKLROOT}/lib/intel64 -traceback -lfftw3 -liomp5 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core
FFLAGS_classic = -fopenmp -c -fast -traceback -pg -i8 -r8 -extend-source -fpe0 -fpconstant
LINK_FLAGS_classic = -qopenmp -L${MKLROOT}/lib/intel64 -traceback -Wl,--start-group -lfftw3 -liomp5 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core

BINDIR=../bin
SRCDIR=../src
OBJ = read-input.o generic_data.o akima_interpol.o utils.o combinations.o init.o analysis.o propa.o mcend-2.1.o

all: $(BINDIR)/mcend

$(BINDIR)/mcend: $(OBJ)
	$(FC) -o $(BINDIR)/mcend $(OBJ) $(LINK_FLAGS)

# debug: FFLAGS = -O0 -warn all -check all -traceback -fp-stack-check -g -fpe:1 -fopenmp -c -traceback -pg -integer-size 64 -real-size 64 -extend-source -fpconstant
debug: FFLAGS = -O0 -fopenmp -c -warn all -traceback -pg -integer-size 64 -real-size 64 -extend-source -fpe:0 -fpconstant
# debug: LINK_FLAGS += -check all
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
gfort: FFLAGS = -O3 -c -pg -g
#gfort: FFLAGS = -O1 -c -pg -g
gfort: LINK_FLAGS = -llapack -lblas -lfftw3 -fopenmp -pg -g
gfort: $(BINDIR)/mcend

debug_gfort: FC = gfortran
debug_gfort: FFLAGS = -O0 -c
debug_gfort: LINK_FLAGS = -llapack -lblas -lfftw3 -fopenmp
#debug_gfort: FFLAGS += -O0 -c -g -p -ffree-line-length-0 -fcheck=all -ffpe-trap=zero,overflow,underflow -finit-real=nan -Wline-truncation -Wall -fbacktrace -fbounds-check
debug_gfort: FFLAGS += -O0 -c -g -p -ffree-line-length-0 -fcheck=all -Wline-truncation -Wall -fbacktrace -fbounds-check
debug_gfort: $(BINDIR)/mcend

clean:
	rm -f *.o *.mod *.inc *.a

%.o:	$(SRCDIR)/%.f90
	$(FC) $(FFLAGS)  $<

