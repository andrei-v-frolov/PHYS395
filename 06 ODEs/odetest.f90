! $Id$
! ODE integration methods tested on a simple anharmonic oscillator
! 
! [compile with: gfortran -O3 -fdefault-real-8 odetest.f90 -o odetest]
! [run as: ./odetest > LOG] [plot in gnuplot: plot 'LOG' using 1:2 with lines ]

program odetest; implicit none

real, parameter :: dt = 7.416298709205487673735401388781040185Q0/2**7

integer l; real :: y(2) = (/ 1.0, 0.0 /)

do l = 1,2**10
	! output time, position, velocity, error in energy
        write (*,'(4g24.16)') (l-1)*dt, y, y(1)**4/4.0 + y(2)**2/2.0 - 0.25
        
        ! step forward with your choice of numerical scheme
	!call rk4(y, 2, (l-1)*dt, dt, y)
	!call si(8, y, dt)
        call gl6(y, dt)
end do


contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! explicit symplectic integrators; these methods work for separable Hamiltonians
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! split Hamiltonian evolution step
subroutine si2(y, dt)
	real y(2), dt
	
	y(1) = y(1) + y(2) * (dt/2.0)
	y(2) = y(2) - y(1)**3 * dt
	y(1) = y(1) + y(2) * (dt/2.0)
end subroutine si2

! we can increase the order by taking a magical sequence of time steps
! see Yoshida's paper: http://dx.doi.org/10.1016/0375-9601(90)90092-3
subroutine si6(y, dt)
	real y(2), dt; integer i
	real, parameter :: w(0:3) = (/ &
		 1.31518632068391121888424972823886251Q0, &
		-1.17767998417887100694641568096431573Q0, &
		 0.235573213359358133684793182978534602Q0, &
		 0.784513610477557263819497633866349876Q0 /)
	
	do i = -3,3; call si2(y, w(abs(i))*dt); end do
end subroutine si6

! here's a general trick which bumps the order by two
! k-th order recursive symplectic integrator (k should be even)
recursive subroutine si(k, y, dt)
        real y(2), dt, gamma, w1, w0; integer k
        
        select case (k)
                case (2); call si2(y, dt)
                case (6); call si6(y, dt)
                case default
                        gamma = 1.0/(k-1); w1 = 1.0/(2.0 - 2.0**gamma); w0 = 1.0 - 2.0*w1
                        call si(k-2, y, w1*dt); call si(k-2, y, w0*dt); call si(k-2, y, w1*dt)
        end select
end subroutine si


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! implicit Gauss-Legendre methods; symplectic with arbitrary Hamiltonian, A-stable
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! evaluate derivatives
subroutine evalf(y, dydx)
        real y(2), dydx(2)
        
        dydx(1) = y(2)
        dydx(2) = -y(1)**3
end subroutine evalf

! 4th order implicit Gauss-Legendre integrator
subroutine gl4(y, dt)
        integer, parameter :: s = 2, n = 2
        real y(n), g(n,s), dt; integer i, k
        
        ! Butcher tableau for 4th order Gauss-Legendre method
        real, parameter :: a(s,s) = reshape((/ 0.25, 0.25 - 0.5/sqrt(3.0), 0.25 + 0.5/sqrt(3.0), 0.25 /), (/s,s/))
        real, parameter ::   b(s) = (/ 0.5, 0.5 /)
        
        ! iterate trial steps
        g = 0.0; do k = 1,16
                g = matmul(g,a)
                do i = 1,s
                        call evalf(y + g(:,i)*dt, g(:,i))
                end do
        end do
        
        ! update the solution
        y = y + matmul(g,b)*dt
end subroutine gl4

! 6th order implicit Gauss-Legendre integrator
subroutine gl6(y, dt)
        integer, parameter :: s = 3, n = 2
        real y(n), g(n,s), dt; integer i, k
        
        ! Butcher tableau for 6th order Gauss-Legendre method
        real, parameter :: a(s,s) = reshape((/ &
                5.0/36.0, 2.0/9.0 - 1.0/sqrt(15.0), 5.0/36.0 - 0.5/sqrt(15.0), &
                5.0/36.0 + sqrt(15.0)/24.0, 2.0/9.0, 5.0/36.0 - sqrt(15.0)/24.0, &
                5.0/36.0 + 0.5/sqrt(15.0), 2.0/9.0 + 1.0/sqrt(15.0), 5.0/36.0 /), (/s,s/))
        real, parameter ::   b(s) = (/ 5.0/18.0, 4.0/9.0, 5.0/18.0/)
        
        ! iterate trial steps
        g = 0.0; do k = 1,16
                g = matmul(g,a)
                do i = 1,s
                        call evalf(y + g(:,i)*dt, g(:,i))
                end do
        end do
        
        ! update the solution
        y = y + matmul(g,b)*dt
end subroutine gl6

! 8th order implicit Gauss-Legendre integrator
subroutine gl8(y, dt)
        integer, parameter :: s = 4, n = 2
        real y(n), g(n,s), dt; integer i, k
        
        ! Butcher tableau for 8th order Gauss-Legendre method
        real, parameter :: a(s,s) = reshape((/ &
                 0.869637112843634643432659873054998518Q-1, -0.266041800849987933133851304769531093Q-1, &
                 0.126274626894047245150568805746180936Q-1, -0.355514968579568315691098184956958860Q-2, &
                 0.188118117499868071650685545087171160Q0,   0.163036288715636535656734012694500148Q0,  &
                -0.278804286024708952241511064189974107Q-1,  0.673550059453815551539866908570375889Q-2, &
                 0.167191921974188773171133305525295945Q0,   0.353953006033743966537619131807997707Q0,  &
                 0.163036288715636535656734012694500148Q0,  -0.141906949311411429641535704761714564Q-1, &
                 0.177482572254522611843442956460569292Q0,   0.313445114741868346798411144814382203Q0,  &
                 0.352676757516271864626853155865953406Q0,   0.869637112843634643432659873054998518Q-1 /), (/s,s/))
        real, parameter ::   b(s) = (/ &
                 0.173927422568726928686531974610999704Q0,   0.326072577431273071313468025389000296Q0,  &
                 0.326072577431273071313468025389000296Q0,   0.173927422568726928686531974610999704Q0  /)
        
        ! iterate trial steps
        g = 0.0; do k = 1,16
                g = matmul(g,a)
                do i = 1,s
                        call evalf(y + g(:,i)*dt, g(:,i))
                end do
        end do
        
        ! update the solution
        y = y + matmul(g,b)*dt
end subroutine gl8

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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! explicit 4-th order Runge-Kutta method; workhorse and baseline for comparison
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! field derivatives on initial surface
subroutine derivs(x, y, dydx)
	real x, y(2), dydx(2)
	
	dydx(1) = y(2)
	dydx(2) = -y(1)**3
end subroutine derivs

! 4th order Runge-Kutta integrator with fixed step size
subroutine rk4(y, n, x, h, yout)
        integer n
        real y(n), x, h, yout(n)
        
        real hh, xh, yt(n), dydx(n), dyt(n), dym(n)
        
        hh = h/2.0
        xh = x + hh
        
        call derivs(x, y, dydx)
        yt = y + hh*dydx
        
        call derivs(xh, yt, dyt)
        yt = y + hh*dyt
        
        call derivs(xh, yt, dym)
        yt = y + h*dym
        dym = dyt + dym
        
        call derivs(x+h, yt, dyt)
        yout = y + h*(dydx + dyt + 2.0*dym)/6.0
end subroutine rk4

end
