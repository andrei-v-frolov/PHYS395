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

write (*,*) maxloc(A(2:,2:))

!call gaussj(3, A, B)

!write (*,*) B

contains

! solve A.x = B using Gauss-Jordan elimination
! A gets destroyed, answer is returned in B
recursive subroutine gaussj(n, A, B)
	integer n, i, j; real A(n,n), B(n)
	
	! pivot indices
	integer row(n), col(n), pivot(2)
	
	forall (i=1:n) row(i) = i
	forall (j=1:n) col(j) = j
	
	! access data like A(row(i),col(j))
	
	! find location of maximal value
	pivot = maxloc(abs(A))
	
	! swap rows/columns
	row([1,pivot(1)]) = row([pivot(1),1])
	col([1,pivot(2)]) = col([pivot(2),1])
	
	! get rid of first column
	do i = 2,n
		A(i,:) = A(i,:) - (A(i,1)/A(1,1)) * A(1,:)
		B(i) = B(i) - (A(i,1)/A(1,1)) * B(1)
	end do
	
	! if you want to iterate by recursion
	call gaussj(n-1, A(2:,2:), B(2:))
	
	! otherwise, do loop on k
end subroutine

end