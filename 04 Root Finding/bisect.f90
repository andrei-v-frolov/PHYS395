program bisect
implicit none

real, parameter :: epsilon = 1.0e-14

real a, b, c, fa, fb, fc

a = 0.0; fa = f(a)
b = 2.5; fb = f(b)

if (fa*fb > 0.0) stop "The root is not bracketed, bailing out..."

do while (abs(b-a) > epsilon)
	c = (a+b)/2.0; fc = f(c)
	
	write (*,*) c, fc
	
	if (fa*fc < 0.0) then
		b = c; fb = fc
	else
		a = c; fa = fc
	end if
end do

contains

! function to find a root of
function f(x)
	real f, x
	f = cos(x)
end function

end program