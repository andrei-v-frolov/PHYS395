! FFTW demo, compile with
! gfortran -O3 -fdefault-real-8 -I/opt/local/include fft.f90 -L/opt/local/lib -lfftw3

program fftdemo; implicit none

include "fftw3.f"

integer, parameter :: N = 2**10
real, parameter :: pi = 3.141592653589793238462643383279502884197169399375Q0

complex(8) in, ft, out
dimension in(N), ft(N), out(N)
integer*8 plan

integer i, k

!in = (0,0)
!in(1:N/2) = (1,0)
!forall (i=1:N) in(i) = exp(2.0*sin(2*pi*(i-1)/N))
forall (i=1:N) in(i) = exp( -1000.0*((i-N/2)*(1.0/N))**2 )

call dfftw_plan_dft_1d(plan,N,in,ft,FFTW_FORWARD,FFTW_ESTIMATE)
call dfftw_execute_dft(plan, in, ft)
call dfftw_destroy_plan(plan)

do k = 1,N/2
	ft(k) = (0,1)*(2.0*pi)*(k-1) * ft(k)		! positive frequencies
	ft(n+1-k) = -(0,1)*(2.0*pi)*k * ft(n+1-k)	! negative frequensies
end do

call dfftw_plan_dft_1d(plan,N,ft,out,FFTW_BACKWARD,FFTW_ESTIMATE)
call dfftw_execute_dft(plan, ft, out)
call dfftw_destroy_plan(plan)

do i = 1,N
	write (*,*) abs(in(i)), abs(ft(i)), real(out(i)/N)
end do

end program