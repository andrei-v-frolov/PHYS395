! minimize.f90  -  minimization in n dimensions
! compile with: gfortran -O3 -fdefault-real-8 minimize.f90

program descent
implicit none

! dimnesionality of the problem
integer, parameter :: n = 2

! tilt of the potential bottom
real, parameter :: mu = 0.1

! value of parameter determines effective damping
! this is very important to get right [more on this later]
real, parameter :: dt = 1.0

real x(n)
integer i

! initial guess
x = 0.0; x(1) = 2.0; x(2) = 1.0

! move in the downward direction
do i = 1,20
	! steepest descent: go down local gradient
	! x = x - df(x)*dt
	
	! Newton's method: go to projected minimum location
	! x = x - matmul(invert2x2(ddf(x)), df(x))*dt
	
	! Fletcher's regularization: modify Hessian matrix
	! this is the basic building block in Levenberg–Marquardt algorithm
	x = x - matmul(invert2x2(regular(ddf(x), 0.1)), df(x))*dt
	
	write (*,*) x, f(x)
end do

contains

! function to find a minimum of...
pure function f(x); intent(in) x
	real f, x(n)
	f = (sum(x*x) - 1.0)**2 + mu*x(1)
end function

! derivative of a function to find a minimum of
pure function df(x); intent(in) x
	real df(n), x(n)
	df = 4.0*(sum(x*x) - 1.0)*x; df(1) = df(1) + mu
end function

! second derivative of a function to find a minimum of
pure function ddf(x); intent(in) x
	real ddf(n,n), x(n); integer i, j
	
	forall (i=1:n,j=1:n) ddf(i,j) = 8.0*x(i)*x(j)
	forall (i=1:n) ddf(i,i) = ddf(i,i) + 4.0*(sum(x*x) - 1.0)
end function

! regularize matrix for inversion
pure function regular(M, lambda); intent(in) M, lambda
	real regular(n,n), M(n,n), lambda; integer i
	
	regular = M; forall (i=1:n) regular(i,i) = (1.0+lambda)*M(i,i)
end function

! invert 2x2 matrix
pure function invert2x2(M); intent(in) M
	real invert2x2(2,2), M(2,2)
	
	invert2x2 = reshape([M(2,2), -M(2,1), -M(1,2), M(1,1)]/(M(1,1)*M(2,2) - M(1,2)*M(2,1)), [2,2])
end function

end program