program euler
implicit none

real, parameter :: omega = 1.0
real, parameter :: dt = 0.1/omega

real x, v, E

integer i

x = 1.0; v = 0.0

do i = 1,1000
	x = x + v*(dt/2.0)
	v = v - (omega**2 * x) * dt
	x = x + v*(dt/2.0)
	
	E = v**2/2.0 + omega**2 * x**2/2.0
	
	write (*,'(4g24.16)') i*dt, x, v, E
end do

end