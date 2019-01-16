! ortho.f90  -  orthogonal polynomials primer
! compile with: gfortran -O3 -fdefault-real-8 ortho.f90

program ortho
implicit none

! number of points to print
integer, parameter :: n = 1000

real x(n), y(n)
integer i, k

forall (i=1:n) x(i) = (2*i-n-1.0)/(n-1)
y = sin(10.0*x)

do k = 0,50
	! output k-th order basis function
	do i = 1,n
                write (*,*) x(i), ChebyshevT(x(i), k)
	end do
	
	! two empty lines separate plots in gnuplot
	write (*,*) ""
	write (*,*) ""
end do

contains

! Chebyshev polynomial T_n(x)
function ChebyshevT(x, n)
	real ChebyshevT, x; integer n
	
	ChebyshevT = cos(n*acos(x))
end function

end
