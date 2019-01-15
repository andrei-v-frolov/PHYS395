! sample program illustrating floating point arithmetics pitfalls
! compile with: gfortran -O3 overflow.f90 -o overflow

program overflow
implicit none

real a, b
integer i

integer, parameter :: n = 500000000

a = 0.0; b = 1.0

do i=1,500000000
	a = a + 1.0
end do

b = b + 5.0e-8

write (*,*) a, a == n
write (*,*) b, b == 1.0

end