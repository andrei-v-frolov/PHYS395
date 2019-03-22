! rayleigh.f90  -  Rayleigh iteration solver for eigenvalue problem
! compile with: gfortran -O3 -fdefault-real-8 rayleigh.f90 -llapack

program rayleigh
implicit none

! order of the spectral scheme
integer, parameter :: n = 150

! scale of compactification
real, parameter :: ell = 1.0

! this is what you think it is...
real, parameter :: pi = 3.1415926535897932384626433832795028842Q0

! collocation grids
real, dimension(n) :: x, theta, psi

! second derivative operator
real, dimension (n,n) :: L, H

real lambda
integer i, k

! initialize spectral operators
call initg(); call initl()

! Hamiltonian in static Schrödinger equation
H = -L/2.0; forall (i=1:n) H(i,i) = -L(i,i)/2.0 + V(x(i))

! initial guess of the wavefunction
psi = 1.0/(1.0+x*x); call dump(psi)

! try to relax using Rayleigh's iteration
do k = 1,64
	lambda = dot_product(psi,matmul(H,psi))/dot_product(psi,psi)
	psi = lsolve(psi,lambda); psi = psi/sqrt(dot_product(psi,psi))
	!write (*,*) lambda
	call dump(psi)
end do

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! potential in the problem we are trying to solve
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! simple harmonic oscillator potential
elemental function V(x); intent(in) x
	real V, x
	
	V = x*x/2.0
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
! Rayleigh itertation solver
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Rayleigh's iteration solving eigenvalue problem:
! [-1/2 d^2/dx^2 + V(x)] psi(x) = lambda psi(x)
function lsolve(psi, lambda)
	real lambda, psi(n), lsolve(n), A(n,n), B(n)
	integer i, pivot(n), status
	
	! linear improvement inflating eigenvector
	A = H; forall (i=1:n) A(i,i) = H(i,i) - lambda
	B = psi
	
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
subroutine dump(psi)
	real psi(n), delta(n); integer i
	
	delta = matmul(H,psi) - lambda*psi
	
	do i = 1,n
		write (*,'(3g24.16)') x(i), psi(i), delta(i)
	end do
	
	write (*,*) ""; write (*,*) ""
end subroutine

end program