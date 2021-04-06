module driver
!use parameters
implicit none
contains
subroutine model_driver(SNOWOP,   &           ! OPTION   SNOW option
                        ETOP,     &           ! OPTION   Percolation option
                        DRAINOP,  &           ! OPTION   Snow option
                        SFROP,    &           ! OPTION   Soil Moisture option
                        BASEOP,   &           ! OPTION   Surface Runoff option
                        SMOP,     &           ! OPTION   Surface Runoff option
                        ntimes,   &           ! FORCING   Number of model timesteps
                        PET,      &           ! FORCING   Potential Evapotranspiration
                        jday,     &           ! FORCING   Day of Year
                        tair,     &           ! FORCING   Air Temperature, length N
                        precip,   &           ! FORCING   Precipitation, length N
                        dt,        &          ! PARAMETER   SNOW
                        elevation, &          ! PARAMETER   SNOW
                        t_base,   &           ! PARAMETER   SNOW01
                        t_power,  &           ! PARAMETER   SNOW01
                        t_melt,   &           ! PARAMETER   SNOW01
                        rvs ,     &           ! PARAMETER   SNOW17
                        uadj,     &           ! PARAMETER   SNOW17
                        mbase,    &           ! PARAMETER   SNOW17
                        mfmax,    &           ! PARAMETER   SNOW17
                        mfmin,    &           ! PARAMETER   SNOW17
                        tipm,     &           ! PARAMETER   SNOW17
                        nmf,      &           ! PARAMETER   SNOW17
                        plwhc,    &           ! PARAMETER   SNOW17
                        pxtemp,   &           ! PARAMETER   SNOW17
                        pxtemp1,  &           ! PARAMETER   SNOW17
                        pxtemp2,  &           ! PARAMETER   SNOW17
                        sm1max,   &           ! PARAMETER   SNOW17
                        sm2max,   &           ! PARAMETER   SNOW17
                        beta,     &           ! PARAMETER   PERCOLATION
                        ku,       &           ! PARAMETER   PERCOLATION
                        c,        &           ! PARAMETER   PERCOLATION
                        sm1Fmax,  &           ! OPTIONAL PARAMETER   PERCOLATION
                        psi,      &           ! OPTIONAL PARAMETER   PERCOLATION
                        alpha,    &           ! OPTIONAL PARAMETER   PERCOLATION
                        ks,       &           ! PARAMETER            BASEFLOW
                        lambda,   &           ! OPTIONAL PARAMETER   BASEFLOW
                        lowercasen, &         ! OPTIONAL PARAMETER   BASEFLOW
                        sweVecOutput)  ! OUTPUT


    ! use soilwater
    ! use percolation
    use snowmodule17
    use evaporation
    use percolation
    use baseflow


    ! MODEL RUN OPTIONS
    integer, intent(in) :: SNOWOP
    integer, intent(in) :: ETOP
    integer, intent(in) :: DRAINOP
    integer, intent(in) :: SFROP
    integer, intent(in) :: BASEOP
    integer, intent(in) :: SMOP

    ! MODEL FORCINGS
    integer, intent(in) :: ntimes
    real, intent(in), dimension(0:ntimes)  :: jday
    real, intent(in), dimension(0:ntimes)  :: PET
    real, intent(in), dimension(0:ntimes)  :: tair
    real, intent(in), dimension(0:ntimes)  :: precip

    ! SNOW PARAMETERS
    real, intent(in) :: dt        ! timestep (hours)
    real, intent(in) :: elevation ! tens of meters
    real, intent(in) :: t_base
    real, intent(in) :: t_power
    real, intent(in) :: t_melt
    real, intent(in) :: rvs
    real, intent(in) :: uadj
    real, intent(in) :: mbase
    real, intent(in) :: mfmax
    real, intent(in) :: mfmin
    real, intent(in) :: tipm
    real, intent(in) :: nmf
    real, intent(in) :: plwhc
    real, intent(in) :: pxtemp
    real, intent(in) :: pxtemp1
    real, intent(in) :: pxtemp2

    ! SOIL PARAMETERS
    real, intent(in) :: sm1max
    real, intent(in) :: sm2max

    ! SURFACE RUNOFF PARAMETERS
    real, intent(in) :: beta

    ! PERCOLATION PARAMETERS
    real, intent(in) :: ku                ! Percolation option A,?
    integer, intent(in) :: c              ! Percolation option A,? exponent
    real, intent(in), optional :: sm1Fmax  ! Percolation option C
    real, intent(in), optional :: psi     ! Percolation option C
    real, intent(in), optional :: alpha   ! Percolation option C

    ! BASEFLOW PARAMETERS
    real, intent(in) :: ks
    real, intent(in), optional :: lambda
    integer, intent(in), optional :: lowercasen

    ! OUTPUT
    real, intent(out), dimension(0:ntimes) :: sweVecOutput


    ! INTERNAL
    integer:: i                     ! timestep
    real :: tair_i                  ! forcing-- tair at time i
    real :: precip_i                ! forcing -- precip at time i
    real :: PET_i                   ! forcing -- potential ET at time i
    real :: sm1                     ! state (sm in top layer)
    real :: sm2                     ! state (sm in bottom layer)
    real :: sm1F                    ! state (free water in sm1?)
    real :: E                       ! flux from sm1 --> atmos. (evaporation)
    real :: q12                     ! flux from sm1 --> sm2
    real :: qb                      ! flux from sm2 --> channel (baseflow )
    real :: q0                      ! parameter (percolation) -- baseflow at soil moisture saturation (constant)


    real, dimension(0:ntimes) :: meltVec  ! output from snow model (either option)
    real, dimension(0:ntimes) :: sweVec   ! output from snow model (either option)


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!           INITIAL CONDITIONS                   !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    q0 = 100. ! not sure how to calc this yet



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!           BEGIN MAIN MODEL EXECUTION           !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    !! NOTE: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! After calling the snowmodel, the "meltVec" and "sweVec"
    !! Populated with melt and swe values for the entire model run
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    select case(SNOWOP)
        case(0)
            call snow17driver(ntimes, &
                              jday, &
                              precip, &
                              tair, &
                              elevation, &
                              dt, &
                              rvs, &
                              uadj, &
                              mbase, &
                              mfmax, &
                              mfmin, &
                              tipm, &
                              nmf, &
                              plwhc, &
                              pxtemp, &
                              pxtemp1, &
                              pxtemp2, &
                              meltVec, &
                              sweVec)
        case(1)
            print*, "NOT IMPLEMENTED"
    end select

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Begin time stepping (1st order explicit Runge-Kutta)
    ! ET, soilwater, flow must be called simultaneously
    ! since they all depend on soil water values
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! do i=0,ntimes
    !     print*, i
    !     ! !!!!!!!!!!!!!!!!!!!!!!!!!
    !     ! !! COMPUTE EVAPORATION !!
    !     ! !!!!!!!!!!!!!!!!!!!!!!!!!
    !     ! select case(ETOP)
    !     !     ! several different case options
    !     !     case(0)
    !     !         ! use ET option "A" from the evaporation.f90 module
    !     !         E = ETa(PET(i), sm1, sm1max)  ! this is just a function, not a subroutine

    !     !     case(1)
    !     !         print*, 'NOT IMPLEMENTED, ETOP=',ETOP
    !     ! end select

    !     ! !!!!!!!!!!!!!!!!!!!!!!!!!
    !     ! !! COMPUTE PERCOLATION !!
    !     ! !!!!!!!!!!!!!!!!!!!!!!!!!
    !     ! select case(DRAINOP)
    !     !     !
    !     !     case(0)
    !     !         q12 = q12a(sm1,  sm1max, ku, c)

    !     !     case(1)
    !     !         print*, 'DRAINOP NOT IMPLEMENTED 1'

    !     !     case(2)
    !     !         ! use option C from percolation.f90 module
    !     !         q12 = q12c(sm1F, sm1Fmax, q0, sm2, sm2max, psi, alpha)

    !     ! end select

    !     ! !!!!!!!!!!!!!!!!!!!!!!!!!
    !     ! !! COMPUTE BASEFLOW    !!
    !     ! !!!!!!!!!!!!!!!!!!!!!!!!!
    !     ! select case(BASEOP)
    !     !     !
    !     !     case(0)
    !     !         qb = Qba(sm2, ks)

    !     !     case(1)
    !     !         qb = Qbb(sm2, sm2max, ks, lowercasen)
    !     ! end select

    !     !!!!!!!!!!!!!!!!!!!!!!!!!!
    !     !! UPDATE SOIL MOISTURE !!
    !     !!!!!!!!!!!!!!!!!!!!!!!!!!

    !     ! select case(SMOP)
    !     !     case(0)

    !     ! end select



    ! end do



    ! There is now a vector of melt...

    ! if (ETOP == 0) then
    !     print*, "ETOP 0"
    ! else if (ETOP == 1) then
    !     print*, "ETOP 1"
    ! else
    !     print*, "NOT IMPLEMENTED"
    ! end if




!  print*, params%nmf

end subroutine model_driver
end module driver





! ! Unit test !!!!!
! program test
! use driver
! implicit none






! end program test




