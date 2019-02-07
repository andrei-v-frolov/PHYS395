! descent-2d.f90  -  steepest descent in n dimensions
! compile with: gfortran -O3 -fdefault-real-8 descent-2d.f90

program descent
implicit none

! dimnesionality of the problem
integer, parameter :: n = 2

! tilt of the potential bottom
real, parameter :: mu = 1.0

! value of parameter determines effective damping
! this is very important to get right [more on this later]
real, parameter :: dt = 0.001

real x(n)
integer i

! initial guess
x = 0.0; x(1) = 2.0; x(2) = 1.0

! move in the downward direction
do i = 1,10000
	x = x - df(x)*dt
	write (*,*) x, f(x)
end do

contains

! function to find a minimum of...
pure function f(x); intent(in) x
	real f, x(n)
	f = (sum(x*x) - 1.0)**2 + mu*x(1)
end function

! derivative of a function to find a minimum of
pure function df(x); intent(in) x
	real df(n), x(n)
	df = 4.0*(sum(x*x) - 1.0)*x; df(1) = df(1) + mu
end function

! second derivative of a function to find a minimum of
pure function ddf(x); intent(in) x
	real ddf(n,n), x(n); integer i, j
	
	forall (i=1:n,j=1:n) ddf(i,j) = 8.0*x(i)*x(j)
	forall (i=1:n) ddf(i,i) = ddf(i,i) + 4.0*(sum(x*x) - 1.0)
end function


end program