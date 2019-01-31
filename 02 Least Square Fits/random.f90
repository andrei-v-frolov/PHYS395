! random number generator demo
! compile with: gfortran -O3 -fdefault-real-8 random.f90

program random
implicit none

integer, parameter :: n = 100
real, parameter :: pi = 3.141592653589793238462643383279502884197169399375105821Q0

real U(n), V(n), A(n), B(n), C(n), X(n), Y(n)
integer i

! initialize random sequence
call random_seed()

! uniformly-distributed random numbers
call random_number(U)
call random_number(V)
call random_number(X)
call random_number(Y)

! Gaussian-distributed random numbers
A = sqrt(-2.0*log(U)) * cos(2*pi*V)
B = sqrt(-2.0*log(U)) * sin(2*pi*V)

! Cauchy-distributed random number
C = tan(pi*(Y-0.5))

! noisy data model
Y = (0.5*X + 1.0) + 0.1*A

do i = 1,n
	write (*,*) X(i), Y(i)
end do

end program