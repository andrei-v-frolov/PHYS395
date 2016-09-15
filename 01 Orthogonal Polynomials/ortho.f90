program ortho; implicit none

integer, parameter :: n = 15, k = 3
real, parameter :: pi = 3.141592653589793238462643383279502884197169399375Q0


real x(n), y(n)
integer i

do i = 1,n
	x(i) = (2*i-n-1.0)/(n-1.0)
	y(i) = f(x(i))
	
	write (*,*) x(i), y(i), sin(2*pi*k*x(i)), ChebyshevT(k, x(i))
end do

contains

pure function f(x)
	real f, x; intent(in) x
	
	f = 1.0/(1.0 + 10.0*x*x)
end function

pure function ChebyshevT(k, x)
	real ChebyshevT, x; integer k; intent(in) k, x
	
	ChebyshevT = cos(k*acos(x))
end function


end program