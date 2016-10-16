! read.f90  -  non-linear least square fit demo
! compile with: gfortran -O3 -fdefault-real-8 read.f90 

program leastsq; implicit none

real, parameter :: pi = 3.141592653589793238462643383279502884197169399375Q0

!real, parameter :: lambda = 0.5

integer, parameter :: n = 10000

real x(n), y(n), alpha(2), H(2,2)

integer i, j, m, status

! read data from stdin, at most n points
do i = 1,n
	read (*,*,iostat=status) x(i), y(i)
	if (status < 0) exit
end do

! check actual data extent
m = i-1; if (m == n) write (0,*) "Read data extent was truncated, recompile with larger n"

! initial parameter guess
alpha = [1,0]

! optimization loop
do i = 1,100
	write (*,*) alpha, chi2(alpha)
	H = ddchi2(alpha)
	!forall (j=1:2) H(j,j) = (1.0+lambda)*H(j,j) + 1.0e-6
	alpha = alpha - 0.5 * matmul(inverse(H), dchi2(alpha))
end do

contains

! model function and its derivatives
elemental function f(x,a,b)
	real x,a,b,f
	intent(in) x,a,b
	
	f = exp(a*cos(2.0*pi*x)) + b
end function

elemental function fa(x,a,b)
	real x,a,b,fa
	intent(in) x,a,b
	
	fa = exp(a*cos(2.0*pi*x)) * cos(2.0*pi*x)
end function

elemental function fb(x,a,b)
	real x,a,b,fb
	intent(in) x,a,b
	
	fb = 1.0
end function

! cost function and its derivatives
function chi2(alpha)
	real alpha(2), chi2
	
	associate( a => alpha(1), b => alpha(2) )
	chi2 = sum((y(1:m) - f(x(1:m),a,b))**2)
	end associate
end function

function dchi2(alpha)
	real alpha(2), dchi2(2)
	
	associate( a => alpha(1), b => alpha(2) )
	dchi2(1) = -2.0*sum((y(1:m) - f(x(1:m),a,b)) * fa(x(1:m),a,b))
	dchi2(2) = -2.0*sum((y(1:m) - f(x(1:m),a,b)) * fb(x(1:m),a,b))
	end associate
end function

function ddchi2(alpha)
	real alpha(2), ddchi2(2,2)
	
	associate( a => alpha(1), b => alpha(2) )
	ddchi2(1,1) = 2.0*sum(fa(x(1:m),a,b)*fa(x(1:m),a,b))
	ddchi2(1,2) = 2.0*sum(fa(x(1:m),a,b)*fb(x(1:m),a,b))
	ddchi2(2,1) = ddchi2(1,2)
	ddchi2(2,2) = 2.0*sum(fb(x(1:m),a,b)*fb(x(1:m),a,b))
	end associate
end function

! invert 2x2 matrix M
function inverse(M)
	real M(2,2), inverse(2,2), det
	
	det = M(1,1)*M(2,2) - M(1,2)*M(2,1)
	
	inverse(1,1) =  M(2,2)/det
	inverse(1,2) = -M(2,1)/det
	inverse(2,1) = -M(1,2)/det
	inverse(2,2) =  M(1,1)/det
end function

end program
