! $Id$
! ODE integration methods tested on a simple anharmonic oscillator
! 
! [compile with: gfortran -O3 -fdefault-real-8 odetest-gfortran.f90 -o odetest]
! [run as: ./odetest > LOG] [plot in gnuplot: plot 'LOG' using 1:2 with lines ]

program odetest; implicit none

real, parameter :: dt = 1.0e-2

integer l; real :: y(2) = (/ 0.0, 0.5 /)

do l = 1,2**10
	! output time, position, velocity, error in energy
        write (*,'(4g24.16)') (l-1)*dt, y, y(1)**4/4.0 + y(2)**2/2.0 - 0.25
        
        ! step forward with your choice of numerical scheme
        call gl8(y, dt)
end do


contains

elemental function F(phi)
    real F, phi
    intent(in) phi
    
    F = (phi**2 - 1.0) * phi
end function


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! implicit Gauss-Legendre methods; symplectic with arbitrary Hamiltonian, A-stable
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! evaluate derivatives
subroutine evalf(y, dydx)
        real y(2), dydx(2)
        
        dydx(1) = y(2)
        dydx(2) = F(y(1))
end subroutine evalf

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

end