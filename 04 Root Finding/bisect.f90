! bisect.f90  -  find a (bracketed) root by bisection
! compile with: gfortran -O3 -fdefault-real-8 bisect.f90

program bisect
implicit none

! tolerance for interval length
real, parameter :: epsilon = 1.0e-14

real a, b, c, fa, fb, fc

! initial interval
a = 0.0; fa = f(a)
b = 2.5; fb = f(b)

if (fa*fb > 0.0) stop "The root is not bracketed, bailing out..."

! bisect the interval
do while (abs(b-a) > epsilon)
	c = (a+b)/2.0; fc = f(c); if (fc == 0.0) exit
	
	write (*,*) c, fc
	
	if (fa*fc < 0.0) then; b = c; fb = fc; end if
	if (fc*fb < 0.0) then; a = c; fa = fc; end if
end do

contains

! function to find a root of...
pure function f(x); intent(in) x
	real f, x
	f = cos(x)
end function

end program