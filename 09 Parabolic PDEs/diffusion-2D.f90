program wave; implicit none

integer, parameter :: n = 2**8		! total points in the grid
integer, parameter :: nx = 2**6		! points you want output
integer, parameter :: tt = 2**14	! total number of time steps
integer, parameter :: nt = 2**5		! output every nt-th time step

real, parameter :: alpha = 0.25
real, parameter :: dx = 1.0/(n-1)	! spatial gridpoint separation
real, parameter :: dt = alpha*dx**2	! timestep to take (total time is t = tt*dt)

integer, parameter :: fields = 1

real samples(fields,n,n,2)

integer l

call init(samples(:,:,:,1))

! evolution loop reusing slices in a ring buffer
do l = 1,tt
	if (mod(l-1,nt) == 0) call dump((l-1)*dt, samples(:,:,:,1))
	call step(samples(:,:,:,1), samples(:,:,:,2))
	call step(samples(:,:,:,2), samples(:,:,:,1))
end do

contains

! source term in the equation's RHS
elemental function rho(x,y)
	real x, y, rho; intent(in) x, y
	
	rho = 0.0
	if ((x-0.5)**2 + (y-0.5)**2 < 0.4**2) rho = 1.0
end function

! initial field profiles
subroutine init(hr)
	real, dimension(fields,n,n) :: hr
	
	hr = 0.0
end subroutine

! step forward in time
subroutine step(hr, up)
	real, dimension(fields,n,n) :: hr, up
	
	integer i, j
	real, parameter :: c1 = alpha, c0 = 1.0 - 4.0*alpha
	
	!$omp parallel do
	do j = 2,n-1; do i = 2,n-1
		up(:,i,j) = c1 * (hr(:,i+1,j) + hr(:,i-1,j) + hr(:,i,j+1) + hr(:,i,j-1)) + c0 * hr(:,i,j) - rho((i-1)*dx,(j-1)*dx)*dt
	end do; end do
	
	up(:,1,:) = 0.0		! Dirichlet boundary condition
	up(:,n,:) = 0.0		! Dirichlet boundary condition
	up(:,:,1) = 0.0		! Dirichlet boundary condition
	up(:,:,n) = 0.0		! Dirichlet boundary condition
end subroutine

! initial field profiles
subroutine dump(t, hr)
	real t, hr(fields,n,n)
	integer i, j, ii, jj
	
	! output fields on sub-sampled grid
	do j = 1,nx; jj = (j-1)*n/nx
	do i = 1,nx; ii = (i-1)*n/nx
		write (*,*) t, ii*dx, jj*dx, hr(:,ii+1,jj+1)
	end do; end do
	
	! two empty lines to animate stuff
	write (*,*) ""; write (*,*) ""
end subroutine

end program