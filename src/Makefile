# The compiler
FC = gfortran
# flags for debugging or for maximum performance, comment as necessary
FCFLAGS = -g -fbounds-check
#FCFLAGS = -O2
# flags forall (e.g. look for system .mod files, required in gfortran)
FCFLAGS += -I/usr/include -I/usr/lib64/gfortran/modules

# libraries needed for linking, unused in the examples
LDFLAGS = -llog4fortran -lsim_base -lsim_vol7d -lsim_volgrid6d

# List of executables to be built within the package
PROGRAMS = prodsim_vg6d_tcorr prodsim_create_projgrib prodsim_v7d_fast_read

# "make" builds all
all: $(PROGRAMS)

# General rule for building prog from prog.o; $^ (GNU extension) is
# used in order to list additional object files on which the
# executable depends
%: %.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

# General rules for building prog.o from prog.f90 or prog.F90; $< is
# used in order to list only the first prerequisite (the source file)
# and not the additional prerequisites such as module or include files
%.o: %.f90
	$(FC) $(FCFLAGS) -c $<

%.o: %.F90
	$(FC) $(FCFLAGS) -c $<

# Utility targets
.PHONY: clean veryclean

clean:
	rm -f *.o *.mod *.MOD

veryclean: clean
	rm -f *~ $(PROGRAMS)
