program spectral; implicit none

real, parameter :: pi = 3.141592653589793238462643383279502884197169399375Q0

integer, parameter :: nn = 150

real x(nn), phi(nn), V(nn), L(nn,nn), H(nn,nn), lambda

integer i

! initialize grid and Laplacian operator matric
call initg(x); V = x*x/2.0
call initl();  H = -L/2.0; forall (i = 1:nn) H(i,i) = H(i,i) + V(i)

! initialize initial solution guess
call evalb(2,x,phi); !call dump()

! try to relax using Newton's iteration
do i = 1,50
	lambda = dot_product(phi, matmul(H,phi))/dot_product(phi,phi)
	write (*,*) lambda
	phi = lsolve(lambda, phi)
	phi = phi/sqrt(dot_product(phi,phi))
	!call dump()
end do

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! initialize collocation grid
subroutine initg(x)
	integer i; real x(nn)
	
	forall(i=1:nn) x(i) = atanh(cos(pi*(nn-i+0.5)/nn))
end subroutine

! evaluate basis function and its derivatives on a grid
subroutine evalb(n, x, Tn, Tnxx)
	integer n; real, dimension(nn) :: x, theta, Tn, Tnxx; optional Tnxx
	
	if (mod(n,2) == 0) then
		theta = n*acos(tanh(x)); Tn = cos(theta) - 1.0
		if (present(Tnxx)) Tnxx = -n*(sinh(x)*sin(theta) + n*cos(theta))/cosh(x)**2
	else
		theta = n*acos(tanh(x)); Tn = cos(theta) - tanh(x)
		if (present(Tnxx)) Tnxx = -n*(sinh(x)*sin(theta) + n*cos(theta))/cosh(x)**2 + 2.0*sinh(x)/cosh(x)**3
	end if
end subroutine

! initialize Laplacian operator matrix
subroutine initl()
        integer i, pivot(nn), status; real A(nn,nn), B(nn,nn)
        
        ! evaluate basis and differential operator values on a grid
        do i = 1,nn; call evalb(i+1, x, A(i,:), B(i,:)); end do
        
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

! solve linear Laplace problem H phi = RHS
function lsolve(lambda, rhs)
        real lsolve(nn), lambda, rhs(nn), A(nn,nn), B(nn,1)
        integer i, pivot(nn), status
        
        A = H; forall(i=1:nn) A(i,i) = A(i,i) - lambda; B(:,1) = rhs
        
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
	!residual = lapl - DV(phi)
	
	do i = 1,nn
		write (*,*) x(i), phi(i), lapl(i)
	end do
	write (*,*) ""; write (*,*) ""
end subroutine

end program