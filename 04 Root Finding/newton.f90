! newton.f90  -  find a root by Newton's method
! compile with: gfortran -O3 -fdefault-real-8 newton.f90

program newton
implicit none

real x
integer i

! initial guess
x = 1.0

! Newton's iteration
do i = 1,8
	x = x - f(x)/df(x)
	write (*,*) x, f(x)
end do

contains

! function to find a root of...
pure function f(x); intent(in) x
	real f, x
	f = cos(x)
end function

! derivative of a function to find a root of
pure function df(x); intent(in) x
	real df, x
	df = -sin(x)
end function

end program