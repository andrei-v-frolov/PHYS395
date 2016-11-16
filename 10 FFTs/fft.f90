! FFTW demo, compile with
! gfortran -O3 -fdefault-real-8 -I/opt/local/include fft.f90 -L/opt/local/lib -lfftw3

program fftdemo; implicit none

include "fftw3.f"

integer, parameter :: N = 2**10

complex(8) in, out, inv
dimension in(N), out(N), inv(N)
integer*8 plan

integer i

in = (0,0)

in(1:N/2) = (1,0)

call dfftw_plan_dft_1d(plan,N,in,out,FFTW_FORWARD,FFTW_ESTIMATE)
call dfftw_execute_dft(plan, in, out)
call dfftw_destroy_plan(plan)

do i = 1,N/2
	out(i) = (0,1)*(i-1) * out(i)		! positive frequencies
	out(n+1-i) = (0,1)*i * out(n+1-i)	! negative frequensies
end do



call dfftw_plan_dft_1d(plan,N,out,inv,FFTW_BACKWARD,FFTW_ESTIMATE)
call dfftw_execute_dft(plan, out, inv)
call dfftw_destroy_plan(plan)

do i = 1,N
	write (*,*) abs(in(i)), abs(out(i)), abs(inv(i))/N
end do

end program