program read; implicit none

integer, parameter :: n = 1500

real x(n), y(n)

integer i, status

! read beam from stdin
do i = 1,n
	read (*,*,iostat=status) x(i), y(i)
	if (status < 0) exit
end do

do i = 1,n
	write (*,*) x(i), y(i)
end do


end program
