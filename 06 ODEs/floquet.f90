! $Id$
! ODE integration methods tested on a simple anharmonic oscillator
! 
! [compile with: gfortran -O3 -fdefault-real-8 odetest-gfortran.f90 -o odetest]
! [run as: ./odetest > LOG] [plot in gnuplot: plot 'LOG' using 1:2 with lines ]

program odetest; implicit none

integer, parameter :: n = 2**7 ! timesteps per period

real, parameter :: P = 6.2831853071795864769252867665590057684Q0 ! period
!real, parameter :: P = 7.416298709205487673735401388781040185Q0 ! period

real, parameter :: dt = P/n

integer, parameter :: na = 2**10, nq = 2**10
real(8), parameter :: arange(2) = [-10.0,10.0], qrange(2) = [-10.0,10.0]

real(4), allocatable :: image(:,:,:)

integer i,j; real a, q; complex mu

allocate(image(2,na,nq))

do i = 1,na; a = arange(1) + (arange(2) - arange(1))*(i-1)/(na-1)
    do j = 1,nq; q = qrange(1) + (qrange(2) - qrange(1))*(j-1)/(nq-1)
                mu = floquet()
                !write (*,'(4g24.16)') q, a, real(mu), aimag(mu)
                image(1,i,j) = real(mu)
                image(2,i,j) = aimag(mu)
        end do
        !write (*,'(a)') ""
end do

call write2fits('mathieu.fit', image, arange, qrange, ['ReMu', 'ImMu'], '(a,q)')

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! implicit Gauss-Legendre methods; symplectic with arbitrary Hamiltonian, A-stable
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


function floquet()
        complex floquet, trace; real u(5); integer l
        
        u = (/ 1.0, 0.0, 0.0, 1.0, 0.0 /)
        
        do l = 1,n
            call gl8(u, dt)
        end do
        
        trace = complex(u(1)+u(4),0.0)
        floquet = acosh(trace/2.0)/P
end function floquet

! evaluate derivatives
subroutine evalf(u, dudt)
        real u(5), dudt(5)
        
        ! principal matrix
        !  [ u(1) u(3) ]
        !  [ u(2) u(4) ]
        
        associate(x1 => u(1), v1 => u(2), x2 => u(3), v2 => u(4), t => u(5))
                dudt(1) = v1
                dudt(2) = -(a - 2.0*q*cos(2*t))*x1
                dudt(3) = v2
                dudt(4) = -(a - 2.0*q*cos(2*t))*x2
                dudt(5) = 1.0
        end associate
end subroutine evalf

! 8th order implicit Gauss-Legendre integrator
subroutine gl8(u, dt)
        integer, parameter :: s = 4, n = 5
        real u(n), g(n,s), dt; integer i, k
        
        ! Butcher tableau for 8th order Gauss-Legendre method
        real, parameter :: a(s,s) = reshape((/ &
                 0.869637112843634643432659873054998518Q-1, -0.266041800849987933133851304769531093Q-1, &
                 0.126274626894047245150568805746180936Q-1, -0.355514968579568315691098184956958860Q-2, &
                 0.188118117499868071650685545087171160Q0,   0.163036288715636535656734012694500148Q0,  &
                -0.278804286024708952241511064189974107Q-1,  0.673550059453815551539866908570375889Q-2, &
                 0.167191921974188773171133305525295945Q0,   0.353953006033743966537619131807997707Q0,  &
                 0.163036288715636535656734012694500148Q0,  -0.141906949311411429641535704761714564Q-1, &
                 0.177482572254522611843442956460569292Q0,   0.313445114741868346798411144814382203Q0,  &
                 0.352676757516271864626853155865953406Q0,   0.869637112843634643432659873054998518Q-1 /), (/s,s/))
        real, parameter ::   b(s) = (/ &
                 0.173927422568726928686531974610999704Q0,   0.326072577431273071313468025389000296Q0,  &
                 0.326072577431273071313468025389000296Q0,   0.173927422568726928686531974610999704Q0  /)
        
        ! iterate trial steps
        g = 0.0; do k = 1,16
                g = matmul(g,a)
                do i = 1,s
                        call evalf(u + g(:,i)*dt, g(:,i))
                end do
        end do
        
        ! update the solution
        u = u + matmul(g,b)*dt
end subroutine gl8

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

end