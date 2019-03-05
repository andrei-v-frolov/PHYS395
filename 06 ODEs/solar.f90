! Direct summation N-body code by Anatoly Klypin (NMSU)
! http://astronomy.nmsu.edu/aklypin/CosSim/index.html

! compile with 'gfortran -O3 solar.f90' or 'ifort -fast solar.f90'

!-------------------------------------------------------------------------------

module Nbody            ! not the best code, but a very short one
implicit none

integer, parameter      :: D = 3                ! spatial dimensions
integer, parameter      :: N = 10               ! number of particles

real*8, dimension(N)    :: M                    ! masses of each particle
real*8, dimension(D,N)  :: X, V, A              ! position, velocity, acceleration

real*8, parameter       :: eps2 = 1.0D-8**2     ! potential softening
real*8                  :: Ekin, Epot, Etot     ! energy of the system

contains

subroutine Accel        !-------------------------------------------------------
    integer i, j
    
    A = 0.0
    
    !$omp parallel do
    do i=1,N
    do j=1,N
        A(:,i) = A(:,i) + M(j)*(X(:,j)-X(:,i))/sqrt( sum( (X(:,j)-X(:,i))**2 + eps2 )**3 )
    end do
    end do
end subroutine Accel

subroutine Energy       !-------------------------------------------------------
    integer i, j
    
    Epot = 0.0
    Ekin = 0.5*sum(M*sum(V**2,1))
    
    do i=1,N-1
    do j=i+1,N
        Epot = Epot - M(i)*M(j)/sqrt( sum( (X(:,j)-X(:,i))**2 + eps2 ) )
    end do
    end do
    
    Etot = Ekin + Epot
end subroutine Energy

end module Nbody

!-------------------------------------------------------------------------------

module SolarSystem      ! Keplerian orbits of Solar system planets (JPL data)
implicit none

! algebraic constants
real*8, parameter :: radian = 57.29577951308232087679815481410517
real*8, parameter :: twopi = 6.283185307179586476925286766559006

! astrodynamic constants [http://ssd.jpl.nasa.gov/?constants]
real*8, parameter :: k = 0.01720209895   ! Gaussian gravitational constant, (AU^3/d^2)½
real*8, parameter :: GM = (k*365.25)**2  ! Heliocentric gravitational constant G*M_sun, (AU^3/y^2)

! mass ratios M(sun)/M(planet system)
real*8, parameter, dimension(9) :: mratio = (/ 6023600.0, 408523.71, 328900.56, &
                     3098708.0, 1047.3486, 3497.898, 22902.98, 19412.24, 1.35e8 /)

! Keplerian elements [http://ssd.jpl.nasa.gov/?planet_pos]
integer, parameter :: aau = 1, ecc = 2, inc = 3, lon = 4, per = 5, asc = 6
real*8, parameter, dimension(6,2,9) :: J2000_elements = reshape((/ &
  0.38709927,  0.20563593,  7.00497902,    252.25032350,  77.45779628,  48.33076593, &
  0.00000037,  0.00001906, -0.00594749, 149472.67411175,   0.16047689,  -0.12534081, &
  0.72333566,  0.00677672,  3.39467605,    181.97909950, 131.60246718,  76.67984255, &
  0.00000390, -0.00004107, -0.00078890,  58517.81538729,   0.00268329,  -0.27769418, &
  1.00000261,  0.01671123, -0.00001531,    100.46457166, 102.93768193,   0.0,        &
  0.00000562, -0.00004392, -0.01294668,  35999.37244981,   0.32327364,   0.0,        &
  1.52371034,  0.09339410,  1.84969142,     -4.55343205, -23.94362959,  49.55953891, &
  0.00001847,  0.00007882, -0.00813131,  19140.30268499,   0.44441088,  -0.29257343, &
  5.20288700,  0.04838624,  1.30439695,     34.39644051,  14.72847983, 100.47390909, &
 -0.00011607, -0.00013253, -0.00183714,   3034.74612775,   0.21252668,   0.20469106, &
  9.53667594,  0.05386179,  2.48599187,     49.95424423,  92.59887831, 113.66242448, &
 -0.00125060, -0.00050991,  0.00193609,   1222.49362201,  -0.41897216,  -0.28867794, &
 19.18916464,  0.04725744,  0.77263783,    313.23810451, 170.95427630,  74.01692503, &
 -0.00196176, -0.00004397, -0.00242939,    428.48202785,   0.40805281,   0.04240589, &
 30.06992276,  0.00859048,  1.77004347,    -55.12002969,  44.96476227, 131.78422574, &
  0.00026291,  0.00005105,  0.00035372,    218.45945325,  -0.32241464,  -0.00508664, &
 39.48211675,  0.24882730, 17.14001206,    238.92903833, 224.06891629, 110.30393684, &
 -0.00031596,  0.00005170,  0.00004818,    145.20780515,  -0.04062942,  -0.01183482/), (/6,2,9/))

contains

! calculate planet position coordinates in J2000 ecliptic plane at time t
! in Julian years past J2000.0, approximation accurate 1800AD - 2050AD
subroutine ecliptic_coord(x, t)
    real*8 x(3,9), r(2,9), t, elements(6,9), omega(9), M(9), E(9); integer i
    
    ! current elements
    elements = J2000_elements(:,1,:) + (t/100.0) * J2000_elements(:,2,:)
    
    ! convert angular elements to radians
    elements(3:6,:) = elements(3:6,:)/radian
    
    ! argument of perihelion and the mean anomaly
    omega = elements(per,:) - elements(asc,:)
    M = elements(lon,:) - elements(per,:)
    M = M - twopi*floor(M/twopi+0.5)
    
    ! eccentric anomaly E, satisfies M = E - e sin(M), found by Newton's method
    E = M + elements(ecc,:)*sin(M)
    do i = 1,16
        E = E + (M - (E - elements(ecc,:)*sin(M)))/(1.0 - elements(ecc,:)*cos(E))
    end do
    
    ! heliocentric coordinates
    r(1,:) = elements(aau,:) * (cos(E) - elements(ecc,:))
    r(2,:) = elements(aau,:) * sqrt(1.0 - elements(ecc,:)**2) * sin(E)
    
    ! rotate to J2000 ecliptic plane
    x(1,:) = (cos(omega)*cos(elements(asc,:)) - sin(omega)*sin(elements(asc,:))*cos(elements(inc,:))) * r(1,:) - &
             (sin(omega)*cos(elements(asc,:)) + cos(omega)*sin(elements(asc,:))*cos(elements(inc,:))) * r(2,:)
    x(2,:) = (cos(omega)*sin(elements(asc,:)) + sin(omega)*cos(elements(asc,:))*cos(elements(inc,:))) * r(1,:) - &
             (sin(omega)*sin(elements(asc,:)) - cos(omega)*cos(elements(asc,:))*cos(elements(inc,:))) * r(2,:)
    x(3,:) = sin(omega)*sin(elements(inc,:)) * r(1,:) + cos(omega)*sin(elements(inc,:)) * r(2,:)
end subroutine ecliptic_coord

subroutine ecliptic_velocity(v, t)
    real*8, dimension(3,9) :: v, a, b, dx1, dx2, dx3
    real*8, parameter :: dt = 1.0D-2/365.25; real*8 t
    
    ! sixth order finite difference of the ecliptic_coord
    call ecliptic_coord(a,t+0.5*dt); call ecliptic_coord(b,t-0.5*dt); dx1 = (a-b)/dt
    call ecliptic_coord(a,t+1.0*dt); call ecliptic_coord(b,t-1.0*dt); dx2 = (a-b)/dt
    call ecliptic_coord(a,t+1.5*dt); call ecliptic_coord(b,t-1.5*dt); dx3 = (a-b)/dt
    
    v = 1.5*dx1 - 0.3*dx2 + dx3/30.0 
end subroutine ecliptic_velocity

end module SolarSystem

!-------------------------------------------------------------------------------

program SimpleSolarSystem     ! dynamical Solar system model
    use Nbody; use SolarSystem; implicit none
    
    real*8, parameter       :: dt = 1.0/365.25  ! simulation step
    real*8, parameter       :: t_end = 10.0     ! simulation end
    
    real*8 :: t, E0, Vcm(3)
    integer i
    
    ! initial conditions at Julian year 2000.0
    M(1) = 1.0; M(2:10) = 1.0/mratio
    X(:,1) = 0.0; call ecliptic_coord(X(:,2:10), t)
    V(:,1) = 0.0; call ecliptic_velocity(V(:,2:10), t)
    
    ! zero out center of mass velocity
    do i = 1,3
        Vcm(i) = sum(M*V(i,:))/sum(M)
        V(i,:) = V(i,:) - Vcm(i)
    end do
    
    ! initial energy of configuration
    t = 0.0; call Energy; E0 = Ekin + GM*Epot
    
    ! time evolution loop
    do while (t < t_end)
        ! dump calculated results once in a while
        if (mod(int(t/dt),7)==0) then
            ! output time, kinetic and potential energy, and fractional total error
            call Energy; write (*,'(4g20.11)') t, Ekin, Epot, (Ekin + GM*Epot - E0)/abs(E0)
            
            ! dump planet positions relative to the Sun
            do i=1,9; write (10+i,'(4g20.11)') t, X(:,i+1)-X(:,1); end do
        end if
        
        ! advance one time step
        call si6(dt); t = t + dt
    end do
contains

! leapfrog (aka Verlet) method (second order, symplectic)
subroutine leapfrog(dt)
    real(8) dt
    
    X = X + V*(dt/2.0)
    call Accel
    V = V + A*(GM*dt) ! remember, G is not unity anymore...
    X = X + V*(dt/2.0)
end subroutine leapfrog

! we can increase the order by taking a magical sequence of time steps
! see Yoshida's paper: http://dx.doi.org/10.1016/0375-9601(90)90092-3
subroutine si6(dt)
    real dt; integer i
    real, parameter :: w(0:3) = (/ &
         1.31518632068391121888424972823886251Q0, &
        -1.17767998417887100694641568096431573Q0, &
         0.235573213359358133684793182978534602Q0, &
         0.784513610477557263819497633866349876Q0 /)
    
    do i = -3,3; call leapfrog(w(abs(i))*dt); end do
end subroutine si6


end program SimpleSolarSystem
