program descent; implicit none

real, parameter :: mu = 0.01
real, parameter :: lambda = 1.0e-1

real fx, x(2), grad(2)
integer i

x = [0.0,3.0]

do i = 1,10000
	write (*,*) x, f(x)
	grad = df(x)
	x = x - lambda*grad/sqrt(sum(grad**2))
end do

contains

function f(x)
	real x(2), f
	
	associate ( x => x(1), y => x(2) )
		f = (x*x + y*y - 1.0)**2 + mu*x
	end associate
end function

function df(x)
	real x(2), df(2)
	
	associate ( x => x(1), y => x(2) )
		df(1) = 4.0*(x*x + y*y - 1.0) * x + mu
		df(2) = 4.0*(x*x + y*y - 1.0) * y
	end associate
end function


end program