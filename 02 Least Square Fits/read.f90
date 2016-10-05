! read.f90  -  read data from a file demo
! compile with: gfortran -O3 -fdefault-real-8 read.f90 

program read; implicit none

integer, parameter :: n = 1500

real x(n), y(n)

integer i, m, status

! read data from stdin, at most n points
do i = 1,n
	read (*,*,iostat=status) x(i), y(i)
	if (status < 0) exit
end do

! check actual data extent
m = i-1; if (m == n) write (0,*) "Read data extent was truncated, recompile with larger n"

! write data back to stdout
do i = 1,m
	write (*,*) x(i), y(i)
end do

end program
