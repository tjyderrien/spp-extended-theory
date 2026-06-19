CC = gcc
CCFLAGS = -O2
FC =  ifort
FCFLAGS = 
LD = $(FC)
LDFLAGS = $(FCFLAGS)
LIB = -lgsl -lgslcblas

OBJ = wrap_gsl.o interspectus.o coeffs.o sum_rule.o  call_parser.o  parser.o read_variable_length.o io.o

all: $(OBJ)
	$(LD) $(LDFLAGS) $(LIB) -o interspectus.x $^

clean:
	rm -f $(OBJ) *.a *.mod *.o *~ interspectus.x

# inference rules

.SUFFIXES: .mod .f90

.f90.o:
	$(FC) $(FCFLAGS) $(FREEFORM) -c $<

.c.o:
	$(CC) $(CCFLAGS) -c $<


