program fibonacci; implicit none

integer, parameter :: n = 50

integer(8) i, F(n)

F(1) = 1; F(2) = 1

do i = 3,n
	F(i) = F(i-1) + F(i-2)
	write (*,*) F(i)
end do

end program