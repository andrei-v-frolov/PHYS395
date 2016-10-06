! regression.f90  -  fit y = a*x + b to (x,y) data read from stdin
! compile with: gfortran -O3 -fdefault-real-8 regression.f90

program regression; implicit none

real, parameter :: pi = 3.141592653589793238462643383279502884197169399375Q0

integer, parameter :: n = 4

real x, y
real A(n,n), B(n), C(n)
integer i, j, status

A = 0.0; B = 0.0

do
	! read data from stdin
	read (*,*,iostat=status) x, y
	C = basis(n, x)
	
	! accumulate moments
	B = B + y*C
	forall (i=1:n,j=1:n) A(i,j) = A(i,j) + C(i)*C(j)
	
	! exit read loop if end of file reached
	if (status < 0) exit
end do

! best fit values
call msolve(n, A, B, -1.0)

write (*,*) B

contains

! evaluate basis functions b_i(x) at point x
pure function basis(n, x)
	integer k, n; real x, basis(n); intent(in) n, x
	
	! Chebyshev polynomials
	forall (k=1:n) basis(k) = cos(2.0*pi*(k-1)*x)
end function

! minimize |A.x - B|^2 using LAPACK canned SVD routine (A gets destroyed, the answer is returned in B)
! rcond determines the effective rank of A as described in LAPACK docs. Pass -1.0 for machine precision
subroutine msolve(n, A, B, rcond)
	integer n; real A(n,n), B(n), S(n), W(6*n), rcond
	
	integer rank, status
	
        ! find solution by singular value decomposition
        status = 0; select case (kind(A))
                case(4); call sgelss(n, n, 1, A, n, B, n, S, rcond, rank, W, 6*n, status)
                case(8); call dgelss(n, n, 1, A, n, B, n, S, rcond, rank, W, 6*n, status)
                case default; call abort
        end select
        
        ! bail at first sign of trouble
        if (status /= 0) call abort
end subroutine

end program
