program bisect
implicit none

real x
integer i

x = 1.0

do i = 1,8
	x = x - f(x)/df(x)
	write (*,*) x, f(x)
end do

contains

! function to find a root of
function f(x)
	real f, x
	f = cos(x)
end function

! derivative of a function to find a root of
function df(x)
	real df, x
	df = -sin(x)
end function

end program