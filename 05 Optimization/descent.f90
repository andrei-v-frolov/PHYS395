! descent.f90  -  steepest descent in 1D
! compile with: gfortran -O3 -fdefault-real-8 descent.f90

program descent
implicit none

real x
integer i

! value of parameter determines effective damping
! this is very important to get right [more on this later]
real, parameter :: dt = 0.1

! initial guess
x = 1.0

! move in the downward direction
do i = 1,1000
	x = x - df(x)*dt
	write (*,*) x, f(x)
end do

contains

! function to find a minimum of...
pure function f(x); intent(in) x
	real f, x
	f = cos(x) + 1.1
end function

! derivative of a function to find a minimum of
pure function df(x); intent(in) x
	real df, x
	df = -sin(x)
end function

end program