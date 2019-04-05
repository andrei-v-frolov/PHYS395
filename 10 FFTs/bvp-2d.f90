! FFTW demo, compile with
! gfortran -O3 -fdefault-real-8 -I/opt/local/include bvp-2d.f90 -L /opt/healpix/lib -lcfitsio -L/opt/local/lib -lfftw3_threads -lfftw3

program fftdemo; implicit none

include "fftw3.f"

integer, parameter :: N = 2**11
real, parameter :: pi = 3.141592653589793238462643383279502884197169399375Q0
real, parameter :: dx = 1.0/(N+1), dk = pi

real(8) in(N,N), out(N,N), dst(N,N)
integer*8 plan

real(4) image(2,N,N)

integer i, j, kx, ky, status

forall (j=1:N,i=1:N) in(i,j) = rho(i*dx - 0.15, j*dx - 0.15)

call dfftw_init_threads(status)
call dfftw_plan_with_nthreads(4)

call dfftw_plan_r2r_2d(plan,N,N,in,dst,FFTW_RODFT00,FFTW_RODFT00,FFTW_ESTIMATE)
call dfftw_execute_r2r(plan, in, dst)
call dfftw_destroy_plan(plan)

forall (ky=1:N,kx=1:N) dst(kx,ky) = -dst(kx,ky)/((kx*dk)**2 + (ky*dk)**2)	! inverse Laplacian

call dfftw_plan_r2r_2d(plan,N,N,dst,out,FFTW_RODFT00,FFTW_RODFT00,FFTW_ESTIMATE)
call dfftw_execute_r2r(plan, dst, out)
call dfftw_destroy_plan(plan)

call dfftw_cleanup_threads()

image(1,:,:) = out/4/(N+1)**2
image(2,:,:) = in

call write2fits('potential.fit', image, [0.0,1.0], [0.0,1.0], ['phi', 'rho'])

!do j = 1,N,16
!	do i = 1,N,16
!		write (*,*) i*dx, j*dx, in(i,j), out(i,j)/(4*(N+1)**2)
!	end do
!	write (*,*) ""
!end do

contains

pure function rho(x,y)
	real x, y, rho; intent(in) x, y
	
	rho = 0.0; if (x*x+y*y < 0.1**2) rho = 1.0
end function

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