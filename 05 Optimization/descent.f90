program descent; implicit none

real, parameter :: mu = 0.01
real, parameter :: lambda = 0.5

real fx, x(2), grad(2)
integer i

x = [0.0,3.0]

do i = 1,100
	write (*,*) x, f(x)
	grad = df(x)
	x = x - lambda * matmul(inverse(ddf(x)), grad)
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

function ddf(x)
	real x(2), ddf(2,2)
	
	associate ( x => x(1), y => x(2) )
		ddf(1,1) = 12.0*x*x + 4.0*y*y - 4.0
		ddf(2,1) = 8.0*x*y
		ddf(1,2) = 8.0*x*y
		ddf(2,2) = 12.0*y*y + 4.0*x*x - 4.0
	end associate
end function

function inverse(M)
	real M(2,2), inverse(2,2), det
	
	det = M(1,1)*M(2,2) - M(1,2)*M(2,1)
	
	inverse(1,1) =  M(2,2)/det
	inverse(1,2) = -M(2,1)/det
	inverse(2,1) = -M(1,2)/det
	inverse(2,2) =  M(1,1)/det
end function

end program