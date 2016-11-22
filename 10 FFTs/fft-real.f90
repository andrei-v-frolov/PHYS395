! FFTW demo, compile with
! gfortran -O3 -fdefault-real-8 -I/opt/local/include fft.f90 -L/opt/local/lib -lfftw3

program fftdemo; implicit none

include "fftw3.f"

integer, parameter :: N = 2**22, NN = N/2 + 1
real, parameter :: pi = 3.141592653589793238462643383279502884197169399375Q0

real(8) in(N), out(N)
complex(8) fft(NN)
integer*8 plan

integer i, k

in = (0,0)

!in(1:N/2) = (1,0)
forall (i=1:N) in(i) = exp(2.0*sin(2*pi*(i-1)/N))


call dfftw_plan_dft_r2c_1d(plan,N,in,fft,FFTW_ESTIMATE)
call dfftw_execute_dft_r2c(plan, in, fft)
call dfftw_destroy_plan(plan)

do k = 1,NN
	fft(k) = (0,1)*(2.0*pi)*(k-1) * fft(k)		! positive frequencies
end do


call dfftw_plan_dft_c2r_1d(plan,N,fft,out,FFTW_ESTIMATE)
call dfftw_execute_dft_c2r(plan, fft, out)
call dfftw_destroy_plan(plan)

do i = 1,N
	!write (*,*) in(i), out(i)/N
end do

end program