program random; implicit none

integer, parameter :: n = 1500
real, parameter :: a = 15.0, b = -3.0

real, parameter :: pi = 3.141592653589793238462643383279502884197169399375Q0

real x(n), y(n)

integer i

call random_number(x)
call random_number(y)

y = tan(pi * (y-0.5))

do i = 1,n
	write (*,*) x(i), a*x(i) + b + y(i)
end do

end program