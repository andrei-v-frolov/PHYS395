! Direct summation N-body code by Anatoly Klypin (NMSU)
! http://astronomy.nmsu.edu/aklypin/CosSim/index.html

! compile with 'gfortran -O3 -fdefault-real-8 -fopenmp nbody.f90' or 'ifort -fast -r8 -openmp nbody.f90'
! this is actually a parallelized code, and it will use all of your CPUs efficiently
! omit OpenMP flags for single-threaded version (if you have a few hundred particles)

module Nbody            ! not the best code, but a very short one
implicit none

integer, parameter      :: D = 3                ! spatial dimensions
integer, parameter      :: N = 100              ! number of particles

real*8, dimension(N)    :: M                    ! masses of each particle
real*8, dimension(D,N)  :: X, V, A              ! position, velocity, acceleration

real*8, parameter       :: eps2 = 1.0D-3**2     ! potential softening
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

program SimpleNBody     !-------------------------------------------------------
    use Nbody; implicit none
    
    real*8, parameter       :: dt = 1.0D-4      ! simulation step
    real*8, parameter       :: t_end = 100.0      ! simulation end
    
    real*8 :: t, E0
    
    ! all particles are of equal mass
    M = 1.0/N
    
    ! initial conditions are random, uniformly distributed
    call random_seed()
    call random_number(X); X = (X-0.5)*2.0
    call random_number(V); V = (V-0.5)/3.0
    
    ! initial energy of configuration
    t = 0.0; call Energy; E0 = Etot
    
    ! time evolution loop
    do while (t < t_end)
        ! dump calculated results once in a while
        if (mod(int(t/dt),100)==0) then
            ! output time, kinetic and potential energy, and fractional total error
            call Energy; write (*,'(4g20.11)') t, Ekin, Epot, (Etot-E0)/abs(E0)
            call DumpX
        end if
        
        ! take a single time step forward (using Euler's method)
        X = X + V*(dt/2.0)
        call Accel
        V = V + A*dt
        X = X + V*(dt/2.0)
        t = t + dt
    end do

contains

subroutine DumpX        !-------------------------------------------------------
    integer i
    
    ! output time and positions of all the particles
    do i=1,N
        write (11,'(4g20.11)') t, X(:,i)
    end do
    
    ! write two empty lines to separate different time steps
    write (11,'(a)') "", ""
end subroutine DumpX

end program SimpleNBody
