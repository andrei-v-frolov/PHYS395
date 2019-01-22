! gaussj.f90  -  starting point for Gauss-Jordan elimination
! compile with: gfortran -O3 -fdefault-real-8 gaussj.f90

program test
implicit none

real A(3,3), B(3)

! matrix to invert
A(1,:) = [1.0, 2.0, 3.0]
A(2,:) = [3.0, 2.0, 1.0]
A(3,:) = [0.0, 1.0, 0.0]

B = [1.0, 2.0, 3.0]

write (*,*) A

call gaussj(3, A, B)

write (*,*) B

contains

! solve A.x = B using Gauss-Jordan elimination
! A gets destroyed, answer is returned in B
subroutine gaussj(n, A, B)
	integer n; real A(n,n), B(n)
end subroutine

end