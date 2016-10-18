program rungekutta; implicit none

integer, parameter :: n = 2

real, parameter :: m = 1.0
real, parameter :: k = 1.0
real, parameter :: dt = 0.1/2.0

real y(n), t, E0

integer i

t = 0.0; y(1) = 1.0; y(2) = 0.0; E0 = E(y,t)

do i = 1,1000
	write (*,*) t, y, (E(y,t)-E0)/E0
	call rk4(y, t, dt)
end do

contains

function F(y,t)
	real t, y(n), F(n)
	
	associate (x => y(1), v => y(2))
		F(1) = v
		F(2) = -(k/m)*x
	end associate
end function

function E(y,t)
	real t, y(n), E
	
	associate (x => y(1), v => y(2))
		E = (m/2.0) * v*v + (k/2.0) * x*x
	end associate
end function

subroutine rk4(y, t, dt)
	real y(n), t, dt
	real f1(n), f2(n), f3(n), f4(n)
	
	f1 = F(y,t)
	f2 = F(y + f1*(dt/2.0), t + (dt/2.0))
	f3 = F(y + f2*(dt/2.0), t + (dt/2.0))
	f4 = F(y + f3*dt, t + dt)
	
	y = y + (f1 + 2.0*f2 + 2.0*f3 + f4)*(dt/6.0)
	t = t + dt
end subroutine

end program