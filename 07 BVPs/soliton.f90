! soliton.f90  -  solve soliton BVP using shooting method
! compile with: gfortran -O3 -fdefault-real-8 soliton.f90 

program soliton
implicit none

! iterations for bisection cycle
integer, parameter :: iterations = 1000

real a, b, c, fa, fb, fc
integer i

!fa = integrate(0.70710678119832193, 10.0, 0.01); stop "Done"

! initial interval
a = 0.5; fa = f(a)
b = 1.5; fb = f(b)

if (fa*fb > 0.0) stop "The root is not bracketed, bailing out..."

! bisect the interval
do i = 1,iterations
	c = (a+b)/2.0; fc = f(c); if (fc == 0.0) exit
	
	! Ridder's variation on the basic method
	c = c + (c-a)*sign(1.0,fa-fb)*fc/sqrt(fc*fc-fa*fb); fc = f(c)
	
	write (*,*) c, fc
	
	if (fa*fc < 0.0) then; b = c; fb = fc; end if
	if (fc*fb < 0.0) then; a = c; fa = fc; end if
end do


contains

! potential
pure function V(phi); intent(in) phi
	real V, phi
	
	V = (phi*phi - 1.0)**2/4.0
end function

! potential derivative
pure function DV(phi); intent(in) phi
	real DV, phi
	
	DV = (phi*phi - 1.0)*phi
end function

! total energy of the dynamical system
pure function E(u); intent(in) u
        real u(2), E
        
        associate( phi => u(1), pi => u(2) )
        E = pi**2/2.0 - V(phi)
        end associate
end function

! evaluate derivatives
subroutine evalf(u, dudt)
        real u(2), dudt(2)
        
        associate( phi => u(1), pi => u(2) )
        dudt(1) = pi; dudt(2) = DV(phi)
        end associate
end subroutine evalf

! 10th order implicit Gauss-Legendre integrator
subroutine gl10(y, dt)
        integer, parameter :: s = 5, n = 2
        real y(n), g(n,s), dt; integer i, k
        
        ! Butcher tableau for 8th order Gauss-Legendre method
        real, parameter :: a(s,s) = reshape((/ &
                  0.5923172126404727187856601017997934066Q-1, -1.9570364359076037492643214050884060018Q-2, &
                  1.1254400818642955552716244215090748773Q-2, -0.5593793660812184876817721964475928216Q-2, &
                  1.5881129678659985393652424705934162371Q-3,  1.2815100567004528349616684832951382219Q-1, &
                  1.1965716762484161701032287870890954823Q-1, -2.4592114619642200389318251686004016630Q-2, &
                  1.0318280670683357408953945056355839486Q-2, -2.7689943987696030442826307588795957613Q-3, &
                  1.1377628800422460252874127381536557686Q-1,  2.6000465168064151859240589518757397939Q-1, &
                  1.4222222222222222222222222222222222222Q-1, -2.0690316430958284571760137769754882933Q-2, &
                  4.6871545238699412283907465445931044619Q-3,  1.2123243692686414680141465111883827708Q-1, &
                  2.2899605457899987661169181236146325697Q-1,  3.0903655906408664483376269613044846112Q-1, &
                  1.1965716762484161701032287870890954823Q-1, -0.9687563141950739739034827969555140871Q-2, &
                  1.1687532956022854521776677788936526508Q-1,  2.4490812891049541889746347938229502468Q-1, &
                  2.7319004362580148889172820022935369566Q-1,  2.5888469960875927151328897146870315648Q-1, &
                  0.5923172126404727187856601017997934066Q-1/), (/s,s/))
        real, parameter ::   b(s) = [ &
                  1.1846344252809454375713202035995868132Q-1,  2.3931433524968323402064575741781909646Q-1, &
                  2.8444444444444444444444444444444444444Q-1,  2.3931433524968323402064575741781909646Q-1, &
                  1.1846344252809454375713202035995868132Q-1]
        
        ! iterate trial steps
        g = 0.0; do k = 1,16
                g = matmul(g,a)
                do i = 1,s
                        call evalf(y + g(:,i)*dt, g(:,i))
                end do
        end do
        
        ! update the solution
        y = y + matmul(g,b)*dt
end subroutine gl10

! integrate equations of motion for a given amount of time
function integrate(pi, t, dt)
	real pi, t, dt, integrate
	real u(2), E0; integer i, n
	
	! start from zero at a given gradient
	u = [0.0, pi]; E0 = E(u)
	
	! number of time steps needed
	n = floor(t/dt)
	
	do i = 1,n
		call gl10(u, dt); if (abs(u(1)) > 100.0) exit
		write (*,'(4g24.16)') i*dt, u(1), u(2), E(u)-E0
	end do
	
	call gl10(u,t-n*dt)
	
	! return state at time t
	integrate = u(1)
end function

! function to find a root of...
function f(x); intent(in) x
	real f, x
	f = integrate(x, 10.0, 0.01) - 1.0
end function


end program