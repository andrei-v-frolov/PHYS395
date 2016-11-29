! Hello world in OpenMP environment, compile with
! gfortran -fopenmp -fdefault-real-8 hello.f90 

program hello
use omp_lib
implicit none

integer, parameter :: n = 2**18
integer i, k, threads
real A(n), S(8)

S = 0.0

forall (i=1:n) A(i) = i

! we are in the master thread when we start
write (*,*) "Hello, world!", omp_get_num_threads(), omp_get_thread_num()

! parallel pragma will run execution block in all threads
!$omp parallel
	write (*,*) "Hello, world in parallel!", omp_get_num_threads(), omp_get_thread_num()
	if (mod(omp_get_thread_num(),2) == 0) then
		write (*,*) "I am even!"
	else
		write (*,*) "I am odd!"
	end if		
!$omp end parallel

! reduction using parallel do construct
!$omp parallel do
do k = 1,8
	do i = k,n,8
		S(k) = S(k) + A(i)
	end do
	
	write (*,*) "Doing partial sum", k, "in thread", omp_get_thread_num(), "result is", S(k)
end do
!$omp end parallel do

write (*,*) "Total sum is", sum(S), "should be", n*(n+1.0)/2.0

end program