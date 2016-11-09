program wave; implicit none

integer, parameter :: n = 2**12		! total points in the grid
integer, parameter :: nx = 2**10	! points you want output
integer, parameter :: tt = n		! total number of time steps
integer, parameter :: nt = 2**3		! output every nt-th time step

real, parameter :: alpha = 0.9
real, parameter :: dx = 1.0/n		! spatial gridpoint separation
real, parameter :: dt = alpha*dx	! timestep to take (total time is t = tt*dt)

real, parameter :: lambda = 1000.0

integer, parameter :: fields = 2

real samples(fields,n,3)

integer l

call init(samples(:,:,1), samples(:,:,2))

! evolution loop reusing slices in a ring buffer
do l = 1,tt
	if (mod(l-1,nt) == 0) call dump((l-1)*dt, samples(:,:,2))
	call step(samples(:,:,1), samples(:,:,2), samples(:,:,3))
	call step(samples(:,:,2), samples(:,:,3), samples(:,:,1))
	call step(samples(:,:,3), samples(:,:,1), samples(:,:,2))
end do

contains

! scalar field potential derivative
elemental function DV(phi)
	real DV, phi; intent(in) phi
	
	DV = lambda*(phi**2 - 1.0)*phi
end function

! initial field profiles
subroutine init(dn, hr)
	real, dimension(fields,n) :: dn, hr
	real, parameter :: v = 0.1, x1 = 0.3, x2 = 0.7
	real, parameter :: w = sqrt(2.0/lambda*(1.0- v*v))
	
	integer i
	
	forall (i=1:n) dn(:,i) = tanh(((i-1)*dx - x1)/w) - tanh(((i-1)*dx - x2)/w) - 1.0
	forall (i=1:n) hr(:,i) = tanh(((i-1)*dx - (x1 + v*dt))/w) - tanh(((i-1)*dx - (x2 - v*dt))/w) - 1.0
end subroutine

! step forward in time
subroutine step(dn, hr, up)
	real, dimension(fields,n) :: dn, hr, up
	
	integer i
	real, parameter :: c1 = alpha**2, c0 = 2.0*(1.0 - alpha**2)
	
	do i = 2,n-1
		up(:,i) = c1 * (hr(:,i+1) + hr(:,i-1)) + c0 * hr(:,i) - dn(:,i) - DV(hr(:,i))*(dt**2)
	end do
	
	up(:,1) = hr(:,1)	! Dirichlet boundary condition
	up(:,n) = hr(:,n)	! Dirichlet boundary condition
end subroutine

! initial field profiles
subroutine dump(t, hr)
	real t, hr(fields,n)
	integer i, j
	
	! output fields on sub-sampled grid
	do j = 1,nx; i = j*n/nx
		write (*,*) t, i*dx, hr(:,i)
	end do
	
	! two empty lines to animate stuff
	write (*,*) ""; write (*,*) ""
end subroutine

end program
