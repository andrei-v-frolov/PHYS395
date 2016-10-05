! random.f90  -  random number generator demo
! compile with: gfortran -O3 -fdefault-real-8 random.f90

program random; implicit none

integer, parameter :: n = 1500
real, parameter :: a = 15.0, b = -3.0

real, parameter :: pi = 3.141592653589793238462643383279502884197169399375Q0

real x(n), y(n)

integer i, seed(16)

! initialize random number generator (use urandom on clusters!)
open (333, file="/dev/random", action='read', access='stream', form='unformatted')
read (333) seed; call random_seed(PUT=seed); close (333)

! random number can be seeded with defaults; gfortran will always produce the same sequence
! call random_seed()

! harvest random numbers
call random_number(x)
call random_number(y)

! shape random number distribution
y = tan(pi * (y-0.5))

do i = 1,n
	write (*,*) x(i), a*x(i) + b + y(i)
end do

end program