! FFTW demo, compile with
! gfortran -O3 -fdefault-real-8 -I/opt/local/include fft.f90 -L/opt/local/lib -lfftw3

program fftdemo; implicit none

include "fftw3.f"

integer, parameter :: N = 2**22, NN = N/2 + 1
real, parameter :: pi = 3.141592653589793238462643383279502884197169399375Q0, dk = pi

real(8) in(N), out(N), dst(N)
integer*8 plan

integer i, k

in = (0,0)

in(3*N/8:5*N/8) = 1.0


call dfftw_plan_r2r_1d(plan,N,in,dst,FFTW_RODFT00,FFTW_ESTIMATE)
call dfftw_execute_r2r(plan, in, dst)
call dfftw_destroy_plan(plan)

do k = 1,N
	dst(k) = -dst(k)/(k*dk)**2	! inverse Laplacian
end do

call dfftw_plan_r2r_1d(plan,N,dst,out,FFTW_RODFT00,FFTW_ESTIMATE)
call dfftw_execute_r2r(plan, dst, out)
call dfftw_destroy_plan(plan)

do i = 1,N,1024
	write (*,*) in(i), out(i)/N/2
end do

end program