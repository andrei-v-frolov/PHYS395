! poly.f90  -  showcase for polynomial approximation techniques
! compile with: gfortran -O3 -fdefault-real-8 poly.f90 -llapack

program expansion; implicit none

integer, parameter :: n = 15, m = 3001

real, parameter :: pi = 3.141592653589793238462643383279502884197169399375Q0

real x(n)	! function is evaluated at these (grid) points
real v(n)	! values of the function f(x) at these points
real B(n,n)	! B(i,j) contains value of j-th basis function at point i

integer i

! initialize grid and evaluate function
x = lobatto(n); v = f(x)

! evaluate the basis on the grid points
do i = 1,n; B(i,:) = basis(n, x(i)); end do

! find coefficients of polynomial expansion
call lsolve(n, B, v)

! dump polynomial expansion at a finer grid
call dump(m)

contains

! function to be approximated
elemental function f(x)
	real f, x; intent(in) x
	
	f = 1.0/(1.0 + 10.0*x*x)
end function

! evaluate basis functions b_i(x) at point x
pure function basis(n, x)
	integer k, n; real x, basis(n); intent(in) n, x
	
	! Chebyshev polynomials
	forall (k=1:n) basis(k) = cos((k-1)*acos(x))
end function

! uniformly spaced grid of n points covering interval [-1,1]
pure function uniform(n)
	integer i, n; real uniform(n); intent(in) n
	
	! uniform grid
	forall (i=1:n) uniform(i) = (2*i-n-1.0)/(n-1.0)
end function

! Gauss-Lobatto grid of n points covering interval [-1,1]
pure function lobatto(n)
	integer i, n; real lobatto(n); intent(in) n
	
	! Gauss-Lobatto grid for Chebyshev polynomials
	forall (i=1:n) lobatto(i) = cos(pi*(n-i)/(n-1.0))
end function

! solve A.x = B using LAPACK canned routine (A gets destroyed, the answer is returned in B)
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

! evaluate polynomial expansion
subroutine dump(m)
	integer i, m; real x(m)
	
	x = uniform(m)
	
	do i = 1,m
		write (*,*) x(i), sum(v*basis(n,x(i)))
	end do
end subroutine

end program
