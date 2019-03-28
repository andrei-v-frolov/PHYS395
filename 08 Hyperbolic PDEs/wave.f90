! wave.f90  -  one dimension wave equation using leapfrog finite difference scheme
! compile with: gfortran -O3 -fdefault-real-8 wave.f90

program wave
implicit none

integer, parameter :: fields = 1	! total number of fields to evolve

integer, parameter ::  n = 2**10	! total number of spatial points
integer, parameter :: tt = 2**15	! total number of time steps
integer, parameter :: nt = 2**2		! output solution every nt steps
integer, parameter :: nx = 2**2		! output solution every nx points

real, parameter :: alpha = 1.0		! parameter setting dt^2/dx^2
real, parameter :: dx = 1.0/(n-1)	! spacing between grid points
real, parameter :: dt = sqrt(alpha)*dx	! time step size

real smp(fields,n,3)

integer l

! initialize the problem
call init(smp(:,:,1),smp(:,:,2)); call dump(0.0, smp(:,:,2))

! time evolution loop
do l=1,tt,3
	call step(smp(:,:,1), smp(:,:,2), smp(:,:,3)); if (mod(l+0,nt) == 0) call dump((l+0)*dt, smp(:,:,3))
	call step(smp(:,:,2), smp(:,:,3), smp(:,:,1)); if (mod(l+1,nt) == 0) call dump((l+1)*dt, smp(:,:,1))
	call step(smp(:,:,3), smp(:,:,1), smp(:,:,2)); if (mod(l+2,nt) == 0) call dump((l+2)*dt, smp(:,:,2))
end do

contains

! initial field profile
subroutine init(prev, now)
	real prev(fields,n), now(fields,n); integer i
	
	forall (i=1:n) now(:,i) = exp( -1000.0*((i-n/2)*dx)**2 )
	
	prev = now
end subroutine

! step solution forward
subroutine step(prev, now, next)
	real prev(fields,n), now(fields,n), next(fields,n); integer i
	
	! use equations of motion for interior of the grid
	do i=2,n-1
		next(:,i) = alpha*(now(:,i-1)+now(:,i+1)) + 2.0*(1.0-alpha)*now(:,i) - prev(:,i)
	end do
	
	! use boundary conditions at boundary points
	next(:,1) = 0.0			! Dirichlet boundary conditions
	next(:,n) = next(:,n-1)		! Neumann boundary conditions
end subroutine

! dump solution at timestep t
subroutine dump(t, phi)
	real t, phi(fields,n); integer i
	
	do i=1,n,nx
		write (*,'(8g32.16e3)') t, (i-1)*dx, phi(:,i)
	end do
	
	write (*,*) ""; write (*,*) ""
end subroutine

end program