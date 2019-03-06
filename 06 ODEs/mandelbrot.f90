! mandelbrot.f90  -  Mandelbrot fractal renderer (a demo for FITS format output)
! compile with: gfortran -O3 -fdefault-real-8 mandelbrot.f90 -lcfitsio [-lcurl]

program mandelbrot
implicit none

! image size, number of iterations, and scan bounds
integer, parameter :: nx = 3 * 2**11, ny = 2**12, iterations = 2**10
real(8), parameter :: xx(2) = [-2.0, 1.0], yy(2) = [-1.0, 1.0]

! data array
real(4) data(1,nx,ny)

! temporary variables
integer i, j, k
complex z, c
real x, y

! Mandelbrot iterator
do i = 1,nx; x = xx(1) + (xx(2)-xx(1))*(i-1)/(nx-1)
        do j = 1,ny; y = yy(1) + (yy(2)-yy(1))*(j-1)/(ny-1)
                z = 0.0; c = complex(x,y)
                
                do k = 0,iterations
                        if (abs(z) > 2.0) exit
                        z = z*z + c
                end do
                
                data(1,i,j) = k
                
                ! make color map continious
                if (k < iterations) data(1,i,j) = k+1 - log(log(abs(z))/log(2.0))/log(2.0)
        end do
end do

! write out image to file
call write2fits('data.fit', log(1+data), xx, yy, ['iterations'], '(x,y)')

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
                                call ftpkys(unit, 'EXTNAME', trim(vars(i))//coords, 'variable stored in extension', status)
                        else
                                call ftpkys(unit, 'EXTNAME', trim(vars(i)), 'variable stored in extension', status)
                        end if
                end if
                if (present(xx)) then
                        call ftpkyj(unit, 'CRPIX1', 1, 'x-axis origin pixel', status)
                        call ftpkyd(unit, 'CRVAL1', xx(1), 14, 'x-axis origin coordinate', status)
                        call ftpkyd(unit, 'CDELT1', (xx(2)-xx(1))/(n(1)-1), 14, 'x-axis increment', status)
                end if
                if (present(yy)) then
                        call ftpkyj(unit, 'CRPIX2', 1, 'y-axis origin pixel', status)
                        call ftpkyd(unit, 'CRVAL2', yy(1), 14, 'y-axis origin coordinate', status)
                        call ftpkyd(unit, 'CDELT2', (yy(2)-yy(1))/(n(2)-1), 14, 'y-axis increment', status)
                end if
        end do
        
        ! clean up
        call ftclos(unit, status)
        call ftfiou(unit, status)
end subroutine write2fits

end program
