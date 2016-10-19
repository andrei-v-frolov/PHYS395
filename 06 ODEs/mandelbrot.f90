! mandelbrot.f90 - 2D image output demo
! gfortran -O3 mandelbrot.f90 -lcfitsio

program mandelbrot; implicit none

integer, parameter :: nx = 1024, ny = 1024, iterations = 1000
real(8), parameter :: xrange(2) = [-2.5,1.0], yrange(2) = [-1.0,1.0]

real(4), allocatable :: image(:,:,:)

real x, y
complex z, c
integer i, j, k

allocate(image(1,nx,ny))

do i = 1,nx; x = xrange(1) + (xrange(2) - xrange(1))*(i-1)/(nx-1)
	do j = 1,ny; y = yrange(1) + (yrange(2) - yrange(1))*(j-1)/(ny-1)
		c = cmplx(x,y); z = (0,0)
		
		do k = 1,iterations
			z = z*z + c
			if (abs(z) > 2.0) exit
		end do
		
		image(1,i,j) = k
	end do
end do

call write2fits('mandelbrot.fit', image, xrange, yrange, ['N'])

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! write array data into FITS file as sequence of image extensions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine write2fits(file, array, xx, yy, vars, coords)
        character(len=*) file, vars(:), coords
        real(4) array(:,:,:); real(8) xx(2), yy(2)
        optional xx, yy, vars, coords
        
        integer i, j, status, unit
        integer :: hdus, naxis = 2, n(2), npix
        integer :: bitpix = -32, group = 1, blocksize = -1
        
        ! data dimansions
        hdus = size(array,1)
        n(1) = size(array,2)
        n(2) = size(array,3)
        npix = n(1)*n(2)
        
        ! delete file if it already exists
        open(unit=1234, iostat=status, file=file, status='old')
        if (status == 0) close(1234, status='delete'); status = 0
        
        ! initialize FITS file
        call ftgiou(unit, status)
        call ftinit(unit, file, blocksize, status)
        
        ! write image extensions
        do i = 1,hdus
                call ftiimg(unit, bitpix, naxis, n, status)
                call ftppre(unit, group, 1, npix, array(i,:,:), status)
                
                if (present(vars)) then
                        if (present(coords)) then
                                call ftpkys(unit, 'EXTNAME', vars(i)//coords, 'variable stored in extension', status)
                        else
                                call ftpkys(unit, 'EXTNAME', vars(i), 'variable stored in extension', status)
                        end if
                end if
                if (present(xx)) then
                        call ftpkyj(unit, 'CRPIX1', 1, 'x-axis origin pixel', status)
                        call ftpkyd(unit, 'CRVAL1', xx(1), 14, 'x-axis origin coordinate', status)
                        call ftpkyd(unit, 'CDELT1', (xx(2)-xx(1))/n(1), 14, 'x-axis increment', status)
                end if
                if (present(yy)) then
                        call ftpkyj(unit, 'CRPIX2', 1, 'y-axis origin pixel', status)
                        call ftpkyd(unit, 'CRVAL2', yy(1), 14, 'y-axis origin coordinate', status)
                        call ftpkyd(unit, 'CDELT2', (yy(2)-yy(1))/n(2), 14, 'y-axis increment', status)
                end if
        end do
        
        ! clean up
        call ftclos(unit, status)
        call ftfiou(unit, status)
end subroutine write2fits

end program