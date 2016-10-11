program bracketing; implicit none

real x

x = bracket(0.5, 1.5)

write (*,*) x

contains

function bracket(a0, c0)
	real a0, c0, bracket
	real a, b, c, d, fa, fb, fc, fd
	real, parameter :: eps = 1.0e-16
	real, parameter :: alpha = 1.5 - sqrt(5.0)/2.0
	
	a = a0; fa = f(a)
	b = a0 + alpha*(c0-a0); fb = f(b) ! adjust!
	c = c0; fc = f(c)
	
	!if ((fb-fa)*(fc-fb) > 0) call abort
	if (fb >= fa .or. fb >= fc) call abort
	
	write (*,*) d, fd
	
	do while (abs(c-a) > eps)
		d = (b+c)/2.0; fd = f(d) ! subdivide longest interval - bc *or* ab!
		write (*,*) d, fd
		
		! test if abd or bdc bracket the minimum
		if ((fb-fa)*(fd-fb) < 0.0) then; c = d; fc = fd; end if
		if ((fd-fb)*(fc-fd) < 0.0) then; a = b; fa = fb; b = d; fb = fd; end if
	end do
	
	bracket = b
end function

elemental function f(x)
	real f, x; intent(in) x
	
	f = (x*x - 1.0)**2 + x
end function

end program
