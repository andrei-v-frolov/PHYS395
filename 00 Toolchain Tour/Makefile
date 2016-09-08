# Example Makefile (run 'make' to build stuff)

# compiler configuration (Intel Fortran)
#FC = ifort
#FFLAGS = -O3 -ipo -xHOST -heap-arrays 256 -r8 -pc80
#LDFLAGS = -static-intel

# compiler configuration (GNU Fortran)
FC = gfortran
FFLAGS = -O3 -fdefault-real-8

# dependencies (identation in these rules HAS to be a TAB, not spaces!)
all: fibonacci fibonacci.dat fibonacci.pdf

fibonacci.dat: fibonacci
	./fibonacci > $@

fibonacci.pdf: fibonacci.dat plot.gpl
	gnuplot plot.gpl

# generic compilation rules; you can leave them alone
%: %.f90
	$(FC) $(FFLAGS) $^ -o $@ $(LDFLAGS)
