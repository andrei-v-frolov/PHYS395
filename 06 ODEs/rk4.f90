! rk4.f90  -  4th order Runge-Kutta method demo
! compile with: gfortran -O3 -fdefault-real-8 rk4.f90

program euler
implicit none

! dynamical system dimension
integer, parameter :: n = 2

real, parameter :: omega = 1.0, beta = 0.001
real, parameter :: dt = 0.1/omega

! state vector (which packs all my variables)
real y(2), E0

integer i

! initial conditions
y = [3.1415, 0.1]; E0 = E(y)

! time evolution
do i = 1,10000
	!y = y + F(y)*dt
	call rk4(y, dt)
	
	write (*,'(4g24.16)') i*dt, y, (E(y)-E0)/E0
end do

contains

! 4th order Runge-Kutta step
subroutine rk4(y, dt)
	real y(n), dt
	real k1(n), k2(n), k3(n), k4(n)
	
	k1 = F(y)*dt
	k2 = F(y + k1/2.0)*dt
	k3 = F(y + k2/2.0)*dt
	k4 = F(y + k3)*dt
	
	y = y + (k1 + 2.0*k2 + 2.0*k3 + k4)/6.0
end subroutine

! force
pure function F(y); intent(in) y
	real y(n), F(n)
	
	associate(x => y(1), v => y(2))
	F(1) = v
	F(2) = - omega**2 * sin(x) - beta*v
	end associate
end function

! energy
pure function E(y); intent(in) y
	real y(n), E
	
	associate(x => y(1), v => y(2))
	E = v**2/2.0 - omega**2 * cos(x) 
	end associate
end function

end
