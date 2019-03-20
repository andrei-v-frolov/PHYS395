! spectral.f90  -  relaxation solver for BVP using spectral basis
! compile with: gfortran -O3 -fdefault-real-8 spectral.f90 -llapack

program spectral
implicit none

! order of the spectral scheme
integer, parameter :: n = 15

! scale of compactification
real, parameter :: ell = 1.0

! this is what you think it is...
real, parameter :: pi = 3.1415926535897932384626433832795028842Q0

! collocation grids
real, dimension(n) :: x, theta, phi

! second derivative operator [d^2/dx^2 phi = matmul(L,phi)]
real, dimension (n,n) :: L

integer k

! initialize spectral operators
call initg(); call initl()

! initial guess
phi = cos(theta); call dump(phi)

! try to relax using Newton's iteration
do k = 1,64
	phi = phi + lsolve(phi); call dump(phi)
end do

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! potential in the problem we are trying to solve
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! potential
elemental function V(phi); intent(in) phi
	real V, phi
	
	V = (phi*phi - 1.0)**2/4.0
end function

! potential derivative
elemental function DV(phi); intent(in) phi
	real DV, phi
	
	DV = (phi*phi - 1.0)*phi
end function

! potential second derivative
elemental function DDV(phi); intent(in) phi
	real DDV, phi
	
	DDV = 3.0*phi*phi - 1.0
end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! spectral grid, basis functions and derivative operators
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! initialize the collocation grid
subroutine initg()
	integer i
	
	forall (i=1:n) theta(i) = pi*(n-i+0.5)/n; x = ell/tan(theta)
end subroutine

! evaluate rational Chebyshev basis on a grid theta
subroutine evalb(n, pts, theta, Tn, Tnx, Tnxx)
        integer n, pts; real, dimension(pts), intent(in) :: theta
        real, dimension(pts), intent(out), optional :: Tn, Tnx, Tnxx
        
        ! Chebyshev basis and its derivatives
        if (present(Tn))   Tn = cos(n*theta)
        if (present(Tnx))  Tnx = n * sin(n*theta) * sin(theta)**2/ell
        if (present(Tnxx)) Tnxx = -n * (n*cos(n*theta)*sin(theta) + 2.0*sin(n*theta)*cos(theta)) * sin(theta)**3/ell**2
end subroutine evalb

! initialize linear spectral derivative operator
subroutine initl()
	integer i, pivot(n), status; real A(n,n), B(n,n)
	
	! evaluate basis and differential operator values on collocation grid
	do i = 1,n
		call evalb(i-1, n, theta, Tn=A(i,:), Tnxx=B(i,:))
	end do
	
        ! find linear operator matrix
        status = 0; select case (kind(A))
                case(4); call sgesv(n, n, A, n, pivot, B, n, status)
                case(8); call dgesv(n, n, A, n, pivot, B, n, status)
                case default; call abort
        end select
        
        ! bail at first sign of trouble
        if (status /= 0) call abort
        
        L = transpose(B)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! residual and Newton itertation solver
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! residual for equation d^2/dx^2 phi = V'(phi)
function residual(phi)
	real phi(n), residual(n)
	
	residual = matmul(L,phi) - DV(phi)
end function

! Newton's iteration improving solution:
! [d^2/dx^2 - V''(phi)] delta = - residual(phi)
function lsolve(phi)
	real phi(n), lsolve(n), A(n,n), B(n)
	integer i, pivot(n), status
	
	! linear improvement
	A = L; forall (i=1:n) A(i,i) = L(i,i) - DDV(phi(i))
	B = -residual(phi)
	
	! find linear operator matrix
        status = 0; select case (kind(A))
                case(4); call sgesv(n, 1, A, n, pivot, B, n, status)
                case(8); call dgesv(n, 1, A, n, pivot, B, n, status)
                case default; call abort
        end select
        
        ! bail at first sign of trouble
        if (status /= 0) call abort
        
        lsolve = B
end function

! dump the solution and its residual
subroutine dump(phi)
	real phi(n), delta(n); integer i
	
	delta = residual(phi)
	
	do i = 1,n
		write (*,'(3g24.16)') x(i), phi(i), delta(i)
	end do
	
	write (*,*) ""; write (*,*) ""
end subroutine

end program