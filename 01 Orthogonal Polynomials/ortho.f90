! ortho.f90  -  orthogonal polynomials primer
! compile with: gfortran -O3 -fdefault-real-8 ortho.f90

program ortho
implicit none

! number of points to print
integer, parameter :: n = 1000

real x(n), y(n)
integer i, k

forall (i=1:n) x(i) = (2*i-n-1.0)/(n-1)

do k = 0,50
	! output k-th order basis function
	do i = 1,n
                write (*,*) x(i), ChebyshevT(x(i), k), LegendreP(x(i), k)
	end do
	
	! two empty lines separate plots in gnuplot
	write (*,*) ""
	write (*,*) ""
end do

contains

! Chebyshev polynomial T_n(x)
elemental function ChebyshevT(x, n)
	real ChebyshevT, x; integer n
	intent(in) x, n
	
	ChebyshevT = cos(n*acos(x))
end function

! Legendre polynomial P_n(x)
pure function LegendreP(x, n)
	real LegendreP, x, P(0:n); integer n, i
	intent(in) x, n
	
	select case(n)
		case(0); LegendreP = 1.0
		case(1); LegendreP = x
		case default
			! recursion formula for P_n(x)
			P(0) = 1.0; P(1) = x
			
			do i = 2,n
				P(i) = ((2*i-1)*x*P(i-1) - (i-1)*P(i-2))/i
			end do
			
			LegendreP = P(n)
	end select
end function

end
