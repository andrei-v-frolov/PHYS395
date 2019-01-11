! compile with: gfortran -O3 -fdefault-real-8 hello.f90 -o hello
! run with:     ./hello > DATA

program hello
implicit none

! number of points to print
integer, parameter :: n = 100
real x(0:n), y(0:n)
integer i

forall (i=0:n) x(i) = i/real(n)
y = sin(10.0*x)

do i = 0,n
    write (*,*) x(i), y(i)
end do

end
