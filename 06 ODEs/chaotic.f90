! chaotic.f90  -  simple chaotic evolution demo
! compile with: gfortran -O3 -fdefault-real-8 -fopenmp chaotic.f90 -lcfitsio [-lcurl]

program chaotic
implicit none

! coupling parameter in the potential
real, parameter :: g2 = 2.0

! image size and scan bounds
integer, parameter :: nx = 3 * 2**5, ny = 2**6
real(8), parameter :: xx(2) = [-1.5, 1.5], yy(2) = [-1.0, 1.0]

! data array (allocatable to avoid problems with system resources)
real(4), allocatable :: data(:,:,:)

! temporary variables
integer i, j
real dt

allocate(data(4,nx,ny))

! timestep resolving fastest timescale in the scan
dt = 0.01/sqrt((1.0+g2)*maxval(yy**2) + g2*maxval(xx*2))

! phase space scan
!$omp parallel do
do j = 1,ny
do i = 1,nx
                data(:,i,j) = integrate(xx(1) + (xx(2)-xx(1))*(i-1)/(nx-1), yy(1) + (yy(2)-yy(1))*(j-1)/(ny-1), 100.0, dt)
end do
end do

! write out image to file
call write2fits('data.fit', data, xx, yy, ['x','y','v','w'], '(x0,y0)')

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! dynamical system and Gauss-Legendre integrator
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! potential for 2D particle motion
pure function V(x,y); intent(in) x, y
	real V, x, y
	
	V = y**4/4.0 + (g2/2.0) * (x*y)**2
end function

! total energy of the dynamical system
pure function E(u); intent(in) u
        real u(4), E
        
        associate( x => u(1), y => u(2) )
        E = u(3)**2/2.0 + u(4)**2/2.0 + V(x,y)
        end associate
end function

! evaluate derivatives
subroutine evalf(u, dudt)
        real u(4), dudt(4)
        
        associate( x => u(1), y => u(2), v => u(3), w => u(4) )
        dudt(1) = v
        dudt(2) = w
        dudt(3) = - (g2*y*y) * x
        dudt(4) = - (y*y + g2*x*x) * y
        end associate
end subroutine evalf

! 10th order implicit Gauss-Legendre integrator
subroutine gl10(y, dt)
        integer, parameter :: s = 5, n = 4
        real y(n), g(n,s), dt; integer i, k
        
        ! Butcher tableau for 8th order Gauss-Legendre method
        real, parameter :: a(s,s) = reshape((/ &
                  0.5923172126404727187856601017997934066Q-1, -1.9570364359076037492643214050884060018Q-2, &
                  1.1254400818642955552716244215090748773Q-2, -0.5593793660812184876817721964475928216Q-2, &
                  1.5881129678659985393652424705934162371Q-3,  1.2815100567004528349616684832951382219Q-1, &
                  1.1965716762484161701032287870890954823Q-1, -2.4592114619642200389318251686004016630Q-2, &
                  1.0318280670683357408953945056355839486Q-2, -2.7689943987696030442826307588795957613Q-3, &
                  1.1377628800422460252874127381536557686Q-1,  2.6000465168064151859240589518757397939Q-1, &
                  1.4222222222222222222222222222222222222Q-1, -2.0690316430958284571760137769754882933Q-2, &
                  4.6871545238699412283907465445931044619Q-3,  1.2123243692686414680141465111883827708Q-1, &
                  2.2899605457899987661169181236146325697Q-1,  3.0903655906408664483376269613044846112Q-1, &
                  1.1965716762484161701032287870890954823Q-1, -0.9687563141950739739034827969555140871Q-2, &
                  1.1687532956022854521776677788936526508Q-1,  2.4490812891049541889746347938229502468Q-1, &
                  2.7319004362580148889172820022935369566Q-1,  2.5888469960875927151328897146870315648Q-1, &
                  0.5923172126404727187856601017997934066Q-1/), (/s,s/))
        real, parameter ::   b(s) = [ &
                  1.1846344252809454375713202035995868132Q-1,  2.3931433524968323402064575741781909646Q-1, &
                  2.8444444444444444444444444444444444444Q-1,  2.3931433524968323402064575741781909646Q-1, &
                  1.1846344252809454375713202035995868132Q-1]
        
        ! iterate trial steps
        g = 0.0; do k = 1,16
                g = matmul(g,a)
                do i = 1,s
                        call evalf(y + g(:,i)*dt, g(:,i))
                end do
        end do
        
        ! update the solution
        y = y + matmul(g,b)*dt
end subroutine gl10

! integrate equations of motion for a given amount of time
function integrate(x, y, t, dt)
	real x, y, t, dt, integrate(4)
	real u(4), E0; integer i, n
	
	! start from a given position at rest
	u = [x, y, 0.0, 0.0]; E0 = V(x,y)
	
	! number of time steps needed
	n = floor(t/dt)
	
	do i = 1,n
		call gl10(u, dt)
		!write (*,'(4g24.16)') i*dt, u(1), u(2), E(u)-E0
	end do
	
	call gl10(u,t-n*dt)
	
	! return state at time t
	integrate = u
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
