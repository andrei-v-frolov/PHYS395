! linear.f90  -  solve linear system of equations using LAPACK
! compile with: gfortran -O3 -fdefault-real-8 linear.f90

program test
implicit none

real A(3,3), B(3)

! matrix to invert
A(1,:) = [1.0, 2.0, 3.0]
A(2,:) = [3.0, 2.0, 1.0]
A(3,:) = [0.0, 1.0, 0.0]

B = [1.0, 2.0, 3.0]

write (*,*) A

call lsolve(3, A, B)

write (*,*) B

contains

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

end