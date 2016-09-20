program test; implicit none

integer, parameter :: n = 2
real A(n,n), B(n)

integer i

B = 0.0

A = 0.0

A(1,1) = 1.0
A(1,2) = 7.0
A(2,1) = 5.0
A(2,2) = 1.0
!A(3,1) = 7.0
!A(3,3) = 1

B(1) = 1.0
B(2) = 3.0

call lsolve(n, A, B)
write (*,*) B


!do i = 1,n
!	write (*,*) A(i,:)
!end do

contains
	
! solve A.x = B using Gauss-Jordan elimination
! A gets destroyed, answer is returned in B
subroutine gaussj(n, A, B)
	integer n; real A(n,n), B(n)
	
	integer i, j
	
	! ... actually solve it! ...
	do i = 1,n
		do j = i+1,n
			A(j,:) = A(j,:) - A(j,i)/A(i,i) * A(i,:)
			!B(j,:) = B(j,:) - A(j,i)/A(i,i) * B(i,:)
			write (*,*) i, j, A(j,:)
		end do
	end do
end subroutine

! solve A.x = B using LAPACK canned routine
! A gets destroyed, answer is returned in B
subroutine lsolve(n, A, B)
	integer n; real A(n,n), B(n)
	
	integer status, pivot(n)
	
        ! find static solution by direct inversion
        status = 0; select case (kind(A))
                case(4); call sgesv(n, 1, A, n, pivot, B, n, status)
                case(8); call dgesv(n, 1, A, n, pivot, B, n, status)
                case default; call abort
        end select
        
        ! bail at first sign of trouble
        if (status /= 0) call abort
end subroutine

end program