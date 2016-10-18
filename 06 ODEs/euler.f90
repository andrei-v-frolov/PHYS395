program euler; implicit none

real, parameter :: m = 1.0
real, parameter :: k = 1.0
real, parameter :: dt = 0.1

real x, v, E0

integer i

x = 1.0; v = 0.0; E0 = E(x,v)

do i = 1,1000
	write (*,*) (i-1)*dt, x, v, (E(x,v)-E0)/E0
	
	x = x + v*(dt/2.0)
	v = v + F(x,v)/m * dt
	x = x + v*(dt/2.0)
end do

contains

! force
function F(x,v)
	real x, v, F
	
	F = - k*x
end function

! energy
function E(x,v)
	real x, v, E
	
	E = (m/2.0) * v*v + (k/2.0) * x*x
end function

end program