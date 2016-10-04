program newtons; implicit none

real x

x = newton(10.0)

!write (*,*) x

contains

function newton(a)
	real a, x, y, newton
	
	integer i
	
	x = a
	
	do i = 1,8
		y = f(x)
		write (*,*) x, y
		x = x - y/df(x)
	end do
end function

elemental function f(x)
	real f, x; intent(in) x
	
	f = x*x - 2.0*x - 3.0
end function

elemental function df(x)
	real df, x; intent(in) x
	
	df = 2.0*x - 2.0
end function


end program
