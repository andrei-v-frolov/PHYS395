! euler.f90  -  Euler and Verlet methods demo
! compile with: gfortran -O3 -fdefault-real-8 euler.f90

program euler
implicit none

real, parameter :: omega = 1.0
real, parameter :: dt = 0.1/omega

real x, v, E0

integer i

! initial conditions
x = 1.0; v = 0.0; E0 = E(x,v)

! time evolution
do i = 1,1000
	x = x + v*(dt/2.0)
	v = v + F(x,v)*dt
	x = x + v*(dt/2.0)
	
	write (*,'(4g24.16)') i*dt, x, v, (E(x,v)-E0)/E0
end do

contains

! force
pure function F(x,v); intent(in) x, v
	real x, v, F
	
	F = - omega**2 * x
end function

! energy
pure function E(x,v); intent(in) x, v
	real x, v, E
	
	E = v**2/2.0 + omega**2 * x**2/2.0
end function

end