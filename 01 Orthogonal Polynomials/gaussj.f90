! gaussj.f90  -  starting point for Gauss-Jordan elimination
! compile with: gfortran -O3 -fdefault-real-8 gaussj.f90

! the supplied code does neither work nor compiles, you should fix it yourself!

program test; implicit none

integer, parameter :: n = 2

real A(n,n), B(n)

A(1,1) = 1.0
A(1,2) = 7.0
A(2,1) = 5.0
A(2,2) = 1.0

B(1) = 1.0
B(2) = 3.0

call gaussj(n, A, B)

write (*,*) B

contains
	
! solve A.x = B using Gauss-Jordan elimination
! A gets destroyed, answer is returned in B
subroutine gaussj(n, A, B)
	integer n; real A(n,n), B(n)
	
	integer i, j
	
	! get the matrix into upper-diagonal form first
	do i = 1,n
		! find a pivot element somewhere...
		
		do j = i+1,n
			! eliminate a row along the lines of...
			A(j,:) = A(j,:) - A(j,i)/A(i,i) * A(i,:)
			B(j,:) = B(j,:) - A(j,i)/A(i,i) * B(i,:)
		end do
	end do
	
	! then solve the equations by back-substitution
	! ...
end subroutine

end program