! random number generator demo

program random
implicit none

integer, parameter :: n = 100
real, parameter :: pi = 3.141592653589793238462643383279502884197169399375105821Q0

real U(n), V(n), A(n), X(n)
integer i

call random_seed()
call random_number(U)
call random_number(V)
call random_number(X)

A = sqrt(-2.0*log(U)) * cos(2*pi*V)

do i = 1,n
	write (*,*) X(i), (0.5*X(i) + 1.0) + 0.1*A(i)
end do

end program