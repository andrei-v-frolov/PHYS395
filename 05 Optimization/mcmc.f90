! mcmc.f90  -  maximum likelihood finder using Markov Chain Monte Carlo
! compile with: gfortran -O3 -fdefault-real-8 mcmc.f90 && a.out > DATA
! animate on-screen: gnuplot movie-screen.gpl

program mcmc
implicit none

! dimnesionality of the problem
integer, parameter :: n = 2

! tilt of the potential bottom
real, parameter :: mu = 3.0

! number of chains to run simultaneously
integer, parameter :: chains = 100

! size of the random step to take
real, parameter :: sigma = 0.1

! number of random step to take
integer, parameter :: steps = 2000

! this is what you think it is...
real, parameter :: pi = 3.141592653589793238462643383279502884197169399375105821Q0

! parameters, proposed update, and likelihood acceptance
real x(n,chains), y(n,chains), A(chains)

integer i, l

! initialize random sequence generator
call random_seed()

! initial parameter values are random and unform
call random_number(x); x = 4.0*x - 2.0

call dump(x)

! random walk
do l = 1,steps
	call random_number(A)
	call step(x, y, sigma)
	
	! accept the proposed step according to Metropolis–Hastings criterium
	do i = 1,chains
		if (f(y(:,i))/f(x(:,i)) >= A(i)) x(:,i) = y(:,i)
	end do
	
	call dump(x)
end do

contains

! likelihood to find a maximum of...
pure function f(x); intent(in) x
	real f, x(n)
	f = exp(-((sum(x*x) - 1.0)**2 + mu*x(1)))
end function

! do random step with Gaussian displacement
subroutine step(x, y, sigma)
	real x(n,chains), y(n,chains), sigma, U(chains), V(chains)
	
	call random_number(U)
	call random_number(V)
	
	y(1,:) = x(1,:) + sigma*sqrt(-2.0*log(U)) * cos(2*pi*V)
	y(2,:) = x(2,:) + sigma*sqrt(-2.0*log(U)) * sin(2*pi*V)
end subroutine

! print current parameter values
subroutine dump(x)
	real x(n,chains); integer i
	
	do i = 1,chains
		write (*,*) x(:,i), f(x(:,i))
	end do
	
	write (*,*) ""
	write (*,*) ""
end subroutine

end program