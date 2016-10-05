! regression.f90  -  fit y = a*x + b to (x,y) data read from stdin
! compile with: gfortran -O3 -fdefault-real-8 regression.f90

program regression; implicit none

real sxx, sxy, sy, sx, s
real x, y, a, b
integer status

sxx = 0; sxy = 0; sx = 0; sy = 0; s = 0

do
	! read data from stdin
	read (*,*,iostat=status) x, y
	
	! accumulate moments
	sxx = sxx + x*x
	sxy = sxy + x*y
	sx = sx + x
	sy = sy + y
	s = s + 1.0
	
	! exit read loop if end of file reached
	if (status < 0) exit
end do

! best fit values
a = (s*sxy - sx*sy)/(s*sxx - sx*sx)
b = -(sx*sxy - sxx*sy)/(s*sxx - sx*sx)

write (*,*) a, b

end program
