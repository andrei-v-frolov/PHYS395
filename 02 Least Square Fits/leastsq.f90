! general linear least squares fit
! compile with gfortran -O3 -fdefault-real-8 leastsq.f90 -llapack

program leastsq
implicit none

! number of basis functions
integer, parameter :: n = 3

real x, y, A(n,n), B(n), U(n)

integer i, j
integer :: status = 0

! initialize accumulators
A = 0.0
U = 0.0

! read data from standard input
do while (status == 0)
	read (*,*,iostat=status) x, y
	if (status < 0) exit
	
	! evaluate basis functions at point x
	call evalb(x, B)
	
	! accumulate least square problem matrices
	forall (i=1:n,j=1:n) A(i,j) = A(i,j) + B(i)*B(j); U = U + y*B
end do

! output accumulated matrices
write (*,*) A
write (*,*) U

! solve for best fit parameters
call lsolve(n, A, U)

! output best fit parameters
write (*,*) U

contains

! basis functions we are fitting
subroutine evalb(x, B)
	real x, B(n)
	
	B = [1.0, x, x*x]
end subroutine

! solve A.x = B using LAPACK xGESV
! A gets destroyed, answer is returned in B
subroutine lsolve(n, A, B)
	integer n, pivot(n), status; real A(n,n), B(n)
	
	! initialize status to all clear
	status = 0
	
	! call appropriate version of xGESV
	select case (kind(A))
		case(4); call sgesv(n, 1, A, n, pivot, B, n, status)
		case(8); call dgesv(n, 1, A, n, pivot, B, n, status)
		case default; call abort
	end select
	
	! abort at first sign of trouble
	if (status /= 0) stop "singular matrix in lsolve()" 
end subroutine

end program