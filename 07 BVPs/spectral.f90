program spectral; implicit none

real, parameter :: pi = 3.141592653589793238462643383279502884197169399375Q0

integer, parameter :: nn = 150

real x(nn), phi(nn), L(nn,nn)

integer i

! initialize grid and Laplacian operator matric
call initg(x)
call initl()

! initialize initial solution guess
call evalb(3,x,phi); call dump()

! try to relax using Newton's iteration
do i = 1,50
	phi = phi - lsolve(DDV(phi), matmul(L,phi) - DV(phi)); call dump()
end do

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! derivative of the potential
elemental function DV(phi)
    real DV, phi; intent(in) phi
    
    DV = (phi**2 - 1.0) * phi
end function

! second derivative of the potential
elemental function DDV(phi)
    real DDV, phi; intent(in) phi
    
    DDV = 3.0*phi**2 - 1.0
end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! initialize collocation grid
subroutine initg(x)
	integer i; real x(nn)
	
	forall(i=1:nn) x(i) = atanh(cos(pi*(nn-i+0.5)/nn))
end subroutine

! evaluate basis function and its derivatives on a grid
subroutine evalb(n, x, Tn, Tnxx)
	integer n; real, dimension(nn) :: x, theta, Tn, Tnxx; optional Tnxx
	
	theta = n*acos(tanh(x)); Tn = cos(theta)
	if (present(Tnxx)) Tnxx = -n*(sinh(x)*sin(theta) + n*cos(theta))/cosh(x)**2
end subroutine

! initialize Laplacian operator matrix
subroutine initl()
        integer i, pivot(nn), status; real A(nn,nn), B(nn,nn)
        
        ! evaluate basis and differential operator values on a grid
        do i = 1,nn; call evalb(i-1, x, A(i,:), B(i,:)); end do
        
        ! find linear differentiation matrix
        status = 0; select case (kind(A))
                case(4); call sgesv(nn, nn, A, nn, pivot, B, nn, status)
                case(8); call dgesv(nn, nn, A, nn, pivot, B, nn, status)
                case default; call abort
        end select
        
        ! bail at first sign of trouble
        if (status /= 0) call abort
        
        ! to evaluate Laplacian of function f, simply do matmul(L,f)
        L = transpose(B)
end subroutine initl

! solve linear Laplace problem [L - m_eff^2] phi = RHS
function lsolve(m2eff, rhs)
        real lsolve(nn), m2eff(nn), rhs(nn), A(nn,nn), B(nn,1)
        integer i, pivot(nn), status
        
        ! set up Laplace equation for massive field
        A = L; do i = 1,nn
                A(i,i) = A(i,i) - m2eff(i)
                B(i,1) = rhs(i)
        end do
        
        ! find static solution by direct inversion
        status = 0; select case (kind(A))
                case(4); call sgesv(nn, 1, A, nn, pivot, B, nn, status)
                case(8); call dgesv(nn, 1, A, nn, pivot, B, nn, status)
                case default; call abort
        end select
        
        ! bail at first sign of trouble
        if (status /= 0) call abort
        
        ! return solution
        lsolve = B(:,1)
end function lsolve

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! dump the the solution
subroutine dump()
	real lapl(nn), residual(nn)
	integer i
	
	lapl = matmul(L, phi)
	residual = lapl - DV(phi)
	
	do i = 1,nn
		write (*,*) x(i), phi(i), lapl(i), residual(i)
	end do
	write (*,*) ""; write (*,*) ""
end subroutine

end program