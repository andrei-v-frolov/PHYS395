! argv.f90  -  command line arguments in Fortran 2003
! compile with: gfortran -O3 argv.f90
! run with: ./a.out 1 2.0 three '(4,5)'

program argv
implicit none

integer i, n, status

! buffer to store arguments in
character(len=8000) :: arg

! parsed argument values
integer k
real x
complex z

! get total number of argument
n = command_argument_count()

do i = 1,n
	call get_command_argument(i, arg)
	write (*,*) trim(arg)
	
	read (arg,*,iostat=status) k
	if (status == 0) write (*,*) "Parses as integer:", k
	
	read (arg,*,iostat=status) x
	if (status == 0) write (*,*) "Parses as real number:", x

	read (arg,*,iostat=status) z
	if (status == 0) write (*,*) "Parses as complex number:", z
end do

end program