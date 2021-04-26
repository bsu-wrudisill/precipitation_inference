module driver
!use parameters
implicit none
contains
subroutine model_driver(SNOWOP,     &         ! OPTION   SNOW option
                        ETOP,       &         ! OPTION   Percolation option
                        DRAINOP,    &         ! OPTION   Snow option
                        SATOP,      &         ! OPTION   Saturated Area Option
                        BASEOP,     &         ! OPTION   Surface Runoff option
                        SMOP1,      &         ! OPTION   Soil Water Option -- Top
                        SMOP2,      &         ! OPTION   Soil Water Option -- Bottom
                        ntimes,     &         ! FORCING   Number of model timesteps
                        PET,        &         ! FORCING   Potential Evapotranspiration
                        jday,       &         ! FORCING   Day of Year
                        tair,       &         ! FORCING   Air Temperature, length N
                        precip,     &         ! FORCING   Precipitation, length N
                        dt,         &         ! PARAMETER   SNOW
                        elevation,  &         ! PARAMETER   SNOW
                        t_base,     &         ! PARAMETER   SNOW01
                        t_power,    &         ! PARAMETER   SNOW01
                        t_melt,     &         ! PARAMETER   SNOW01
                        rvs ,       &         ! PARAMETER   SNOW17
                        uadj,       &         ! PARAMETER   SNOW17
                        mbase,      &         ! PARAMETER   SNOW17
                        mfmax,      &         ! PARAMETER   SNOW17
                        mfmin,      &         ! PARAMETER   SNOW17
                        tipm,       &         ! PARAMETER   SNOW17
                        nmf,        &         ! PARAMETER   SNOW17
                        plwhc,      &         ! PARAMETER   SNOW17
                        pxtemp,     &         ! PARAMETER   SNOW17
                        pxtemp1,    &         ! PARAMETER   SNOW17
                        pxtemp2,    &         ! PARAMETER   SNOW17
                        sm1max,     &         ! PARAMETER   SNOW17
                        sm2max,     &         ! PARAMETER   SNOW17
                        ku,         &         ! PARAMETER   PERCOLATION
                        c,          &         ! PARAMETER   PERCOLATION
                        sm1Fmax,    &         ! PARAMETER   PERCOLATION  --- OPTIONAL
                        psi,        &         ! PARAMETER   PERCOLATION  --- OPTIONAL
                        alpha,      &         ! PARAMETER   PERCOLATION  --- OPTIONAL
                        ks,         &         ! PARAMETER   BASEFLOW
                        lambda,     &         ! PARAMETER   BASEFLOW  --- OPTIONAL
                        lowercasen, &         ! PARAMETER   BASEFLOW  --- OPTIONAL
                        beta,       &         ! PARAMETER   SFROFF
                        sweVecOutput, &       ! OUTPUT
                        qVecOutput)          ! OUTPUT


    ! use soilwater
    ! use percolation
    use snowmodule17
    use evaporation
    use percolation
    use baseflow
    use surfacewater
    use overflow
    use soilwater

    ! MODEL RUN OPTIONS
    integer, intent(in) :: SNOWOP
    integer, intent(in) :: ETOP
    integer, intent(in) :: DRAINOP
    integer, intent(in) :: SATOP
    integer, intent(in) :: BASEOP
    integer, intent(in) :: SMOP1
    integer, intent(in) :: SMOP2

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

    ! SFROFF PARAMETERS
    real, intent(in) :: beta

    ! OUTPUT
    real, intent(out), dimension(0:ntimes) :: sweVecOutput
    real, intent(out), dimension(0:ntimes) :: qVecOutput


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
    real :: Ac                      ! Saturated area
    real :: qsx                     ! Overland flow
    real :: qin                     ! rain + melt
    real :: qif                     ! "interflow" (?)
    real :: qufof                   ! saturation excess flow

    real, dimension(0:ntimes) :: meltVec  ! output from snow model (either option)
    real, dimension(0:ntimes) :: sweVec   ! output from snow model (either option)


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!           INITIAL CONDITIONS                   !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    q0 = 30. ! not sure how to calc this yet
    qif = 0   ! interflow (?) is always zero...
    sm1 = 100.
    sm2 = 100.


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!           BEGIN MAIN MODEL EXECUTION           !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    !! NOTE: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! After calling the snowmodel, the "meltVec" and "sweVec"
    !! Populated with melt and swe values for the entire model run
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    select case(SNOWOP)
        case(0)
            call snow17driver(ntimes, jday, precip, tair, &
                              nlayers, opg, dz, stat_elev, &
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

    do i=0,ntimes
        ! get the snowmelt and the rain ...
        qin = meltVec(i) ! + rain(?)

        ! COMPUTE EVAPORATION
        select case(ETOP)
            ! several different case options
            case(0);  E = ETa(PET(i), sm1, sm1max)  ! this is just a function, not a subroutine
        end select


        ! COMPUTE PERCOLATION
        select case(DRAINOP)
            case(0); q12 = q12a(sm1,  sm1max, ku, c)
            case(2); q12 = q12c(sm1F, sm1Fmax, q0, sm2, sm2max, psi, alpha)
        end select


        ! COMPUTE BASEFLOW
        select case(BASEOP)
            case(0);  qb = Qba(sm2, ks)
            case(1);  qb = Qbb(sm2, sm2max, ks, lowercasen)
        end select

        ! COMPUTE SURFACE RUNOFF  !
        ! Select the method for computing the saturated area
        select case(SATOP)
            case(0); Ac = satarea_a(sm1, sm1max, beta)
        end select

        !Now compute the overland flow
        qsx = Ac*qin

        ! UPDATE SOIL MOISTURE TOP LAYER
        select case(SMOP1)
            case(0); call s1_topvic(qin, E, qsx, qif, q12, qufof, sm1)
        end select

        ! UPDATE SOIL MOISTURE BOTTOM LAYER
        select case(SMOP2)
            case(0); call s2_topprms(q12, qb, sm2)
        end select

        ! Compute the total q
        qVecOutput(i) = qb + qsx
    end do



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




