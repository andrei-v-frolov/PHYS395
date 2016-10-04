program bisection; implicit none

real x

x = bisect(0.0, 5.0)

write (*,*) x

contains

function bisect(a0, b0)
	real a0, b0, bisect
	real a, b, c, fa, fb, fc
	real, parameter :: eps = 1.0e-16
	
	a = a0; fa = f(a)
	b = b0; fb = f(b)
	
	if (fa*fb > 0) call abort
	
	do while (abs(b-a) > eps)
		c = (a+b)/2.0; fc = f(c); if (fc == 0.0) exit
		if (fa*fc < 0.0) then; b = c; fb = fc; end if
		if (fc*fb < 0.0) then; a = c; fa = fc; end if
	end do
	
	bisect = c
end function

elemental function f(x)
	real f, x; intent(in) x
	
	f = x*x - 2.0*x - 3.0
end function

end program
