module setprob_module

    ! horizontal forces:
    integer :: xforce_m
    real(kind=8), allocatable, dimension(:)  :: xforce_t, xforce_accel
    real(kind=8) :: xsalt,usalt
    integer, parameter :: iunit_salt = 49
    integer, parameter :: iunit_energy = 51

    save


contains

    subroutine setprob

        implicit none
        integer :: i

        !open(unit=7, file='Radial_accel_R001.txt', status='old', form='formatted')
        !open(unit=7, file='Radial_accel_R032_45.txt', status='old', form='formatted')
        !open(unit=7, file='Radial_accel_R032_10_T40-500s.txt', status='old', form='formatted')
        !open(unit=7, file='Radial_accel_Filter-period_100-500s_R032_Chicxulub_100_-Eq-Mw11.txt', status='old', form='formatted')
        !open(unit=7, file='Radial_accel_Filter-period_40-500s_R032_Chicxulub_100_-Eq-Mw11.txt', status='old', form='formatted')
        open(unit=7, file='Radial_accel_R032_HD60s_source-depth_10km_filter_T40-500s.txt', status='old', form='formatted')
        read(7,*) xforce_m
        write(6,*) '+++ xforce_m = ',xforce_m
        allocate(xforce_t(xforce_m), xforce_accel(xforce_m))

        do i=1,xforce_m
            read(7,*) xforce_t(i), xforce_accel(i)
            !xforce_t(i) = xforce_t(i) - 400.d0 ! shift in time to match prior
            xforce_t(i) = xforce_t(i)
            !xforce_accel(i) = 5.d0 * xforce_accel(i)   ! amplify !
        enddo

        close(7)

        open(unit=iunit_salt, file='salt.txt', status='unknown', form='formatted')
        xsalt = 0.d0  ! initial location of river/saltwater boundary
        usalt = 0.d0

        open(unit=iunit_energy, file='energy.txt', status='unknown', form='formatted')
        write(iunit_energy,*) '# t, PE, KE, total energy'

    end subroutine setprob

end module setprob_module
