module driver
!use parameters
implicit none

contains

! Numerical functions

!!!! F2PY does not allow for allocatable function inputs!!!
!!!! major limitation..... !!!
subroutine convolve(nb, nx, ny, bb, xx, yy)
    integer, intent(in) :: nb             !# number of coefficients in filter
    integer, intent(in) :: nx             !# number of coefficients in input
    integer, intent(in) :: ny             !# number of coefficients in output will be nx+nb-1
    real, intent(in), dimension(nb) :: bb         !# filter coefficients
    real, intent(in), dimension(nx) :: xx         !# input trace
    real, intent(out), dimension(ny) :: yy

    !internal
    integer :: ib, ix, iy
    !integer, dimension(ny) :: yy

    yy = 0.0
    !ny = nx + nb -1
    do ib = 1, nb
        do ix = 1, nx
            yy(ix+ib-1) = yy( ix+ib-1) + xx(ix) * bb(ib)
        end do
    end do
end subroutine convolve


! This is the gamma function kernel used to convolve
! runoff signal to produce streamfow
real function ht(t, k, n) !=3.5, N=4)
    real, intent(in) :: t, k, N
    ht = (t/k)**(N-1) * EXP(-t/k)/k*GAMMA(N)
end function ht

! Helper functions
function PET(L, t, etpar)
    real :: PET, L, t, etpar
    real :: sat_vap_pres, abs_humidity
    ! Potential evapotranspiration calculator
    ! L: Length of day
    ! rho_v_sat: staturated absolute humidity
    ! etpar: parameter
    sat_vap_pres = (.6112)*exp((17.67*t)/(t + 243.5))
    abs_humidity = (2165 * sat_vap_pres)/(t + 273.15)
    PET=L*abs_humidity*etpar
end function PET



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
                        nlayers,    &         ! PARAMETER   SNOW17
                        rvs,        &         ! PARAMETER   SNOW17
                        opg_method, &         ! PARAMETER   SNOW17
                        dz,         &         ! PARAMETER   SNOW17
                        dt,         &         ! PARAMETER   SNOW17
                        opg,        &         ! PARAMETER   SNOW17
                        bias,       &         ! PARAMETER   SNOW17
                        uadj,       &         ! PARAMETER   SNOW17
                        mbase,      &         ! PARAMETER   SNOW17
                        mfmax,      &         ! PARAMETER   SNOW17
                        mfmin,      &         ! PARAMETER   SNOW17
                        tipm,       &         ! PARAMETER   SNOW17
                        nmf,        &         ! PARAMETER   SNOW17
                        plwhc,      &         ! PARAMETER   SNOW17
                        pxtemp,     &         ! PARAMETER   SNOW17
                        pxtemp1,    &         ! PARAMETER   SNOW17
                        pxtemp2,    &         ! PARAMETER  SNOW17
                        sm1max,     &         ! PARAMETER   ?
                        sm2max,     &         ! PARAMETER   ?
                        ku,         &         ! PARAMETER   PERCOLATION
                        c,          &         ! PARAMETER   PERCOLATION
                        sm1Fmax,    &         ! PARAMETER   PERCOLATION  --- OPTIONAL
                        psi,        &         ! PARAMETER   PERCOLATION  --- OPTIONAL
                        alpha,      &         ! PARAMETER   PERCOLATION  --- OPTIONAL
                        ks,         &         ! PARAMETER   BASEFLOW
                        lam,     &         ! PARAMETER   BASEFLOW  --- OPTIONAL
                        lowercasen, &         ! PARAMETER   BASEFLOW  --- OPTIONAL
                        beta,       &         ! PARAMETER   SFROFF
                        Nr, &
                        kr, &
                        qVecOutput, chanVecOutput,sweVecOutput)          ! OUTPUT





    use snowmodule17
    ! use percolation
    ! use evaporation
    ! use baseflow
    ! use surfacewater
    ! use overflow
    ! use soilwater

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
    real, intent(in), dimension(ntimes) :: PET
    real, intent(in), dimension(ntimes) :: jday
    real, intent(in), dimension(ntimes) :: tair
    real, intent(in), dimension(ntimes) :: precip ! change this later. for testing purposes

    ! SNOW17 PARAMETERS
    integer, intent(in) :: nlayers
    integer, intent(in) :: rvs
    integer, intent(in) :: opg_method
    real, intent(in), dimension(nlayers)  :: dz
    real, intent(in) :: dt
    real, intent(in) :: opg
    real, intent(in) :: bias
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

    !SOIL PARAMETERS
    real, intent(in) :: sm1max
    real, intent(in) :: sm2max

    ! PERCOLATION PARAMETERS
    real,    intent(in) :: ku              ! Percolation option A,?
    real, intent(in) :: c                  ! Percolation option A,? exponent
    real, intent(in), optional :: sm1Fmax  ! Percolation option C
    real, intent(in), optional :: psi      ! Percolation option C
    real, intent(in), optional :: alpha    ! Percolation option C

    ! BASEFLOW PARAMETERS
    real, intent(in) :: ks
    real, intent(in), optional :: lam
    real, intent(in), optional :: lowercasen

    ! SFROFF PARAMETERS
    real, intent(in) :: beta


    ! ROUTING PARAMETERS
    real, intent(in) :: Nr
    real, intent(in) :: kr

    ! OUTPUT
    real, intent(out), dimension(ntimes) :: qVecOutput
    real, intent(out), dimension(ntimes) :: chanVecOutput
    real, intent(out), dimension(nlayers,ntimes) :: sweVecOutput

    ! INTERNAL
    integer:: i                     ! timestep
    integer :: ny
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
    real :: dsm1dt
    real :: dsm2dt
    real, dimension(ntimes) :: htv ! routing kernel
    real, dimension(ntimes*2-1) :: cnvrt

    ! SNOW MODEL INPUTS
    real, dimension(nlayers,ntimes) :: outflowVec
    real, dimension(nlayers,ntimes) :: sweVec
    real, dimension(nlayers,ntimes) :: rainVec
    real, dimension(nlayers,ntimes) :: ptotVec

    ! Collapse Melt into the single step
    real, dimension(ntimes) :: outflowVecTotal

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!           INITIAL CONDITIONS                   !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    q0 = 30. ! not sure how to calc this yet
    qif = 0   ! interflow (?) is always zero...
    sm1 = 100.
    sm2 = 700.
    qsx = 0.0


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!           BEGIN MAIN MODEL EXECUTION           !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    !! NOTE: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! After calling the snowmodel, the "meltVec" and "sweVec"
    !! Populated with melt and swe values for the entire model run
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    select case(SNOWOP)
        case(0)

        call snow17driver(ntimes, jday, precip, tair, & ! INPUTS
                          nlayers, dz, dt, rvs,                & ! INPUTS
                          OPG_METHOD, opg, bias, &               ! parameters
                          uadj, &                                ! parameters
                          mbase, &                               ! parameters
                          mfmax, &                               ! parameters
                          mfmin, &                               ! parameters
                          tipm, &                                ! parameters
                          nmf, &                                 ! parameters
                          plwhc, &                               ! parameters
                          pxtemp, &                              ! parameters
                          pxtemp1, &                             ! parameters
                          pxtemp2, &                             ! parameters
                          outflowVec, &                          ! OUTPUT
                          sweVec, &                              !  OUTPUT
                          rainVec, &
                          ptotVec)                               !  OUTPUT

        ! print*, "SNOWOP NOT IMPLEMENTED"
        sweVecOutput = sweVec
        case(1)
            print*, "SNOWOP NOT IMPLEMENTED"
    end select

    ! Compute the total outflow for all snow layers...
    do i=1,ntimes
        outflowVecTotal(i) = SUM(outflowVec(:,i))
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Begin time stepping (1st order explicit Runge-Kutta)
    ! ET, soilwater, flow must be called simultaneously
    ! since they all depend on soil water values
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do i=1,ntimes
        ! get the snowmelt and the rain ...
        qin = outflowVecTotal(i)/3.0 ! + rain(?)

        ! compute PET
        E = PET(i)*min(sm1,sm1max)/sm1max

        ! Compute Percolation
        q12 = ku * (sm1/sm1max)**c

        ! Compute baseflow
        qb = MIN(sm2, ks*(sm2/sm2max)**lowercasen)

        ! Compute saturated area
        Ac = 1 - (1 - sm1/sm1max) ** beta
        print*, Ac

        ! Compute surface runoff
        qsx = Ac*qin

        ! Update soil moisture in both layers
        ! top layer change
        sm1 = sm1 + (qin - qsx) - E - q12 - qif

        ! bottom layer change
        sm2 = sm2 + q12 - qb

        ! store streamflow
        qVecOutput(i) = qb + qsx

    end do

    ! ROUTING MODULE
    ! now convert qVecOutput to streamflow...
    do i=1,ntimes
        htv(i) = ht(REAL(i),Nr,Kr)
    end do

    ! normalize it by 1 -- so as to not add/subtract flow
    htv = htv/SUM(htv)
    ny = ntimes*2-1
    call convolve(ntimes, ntimes, ny, htv, qVecOutput, cnvrt)

    ! the convolution adds on a bunch of points at hte end that we don't need
    do i=1,ntimes
        chanVecOutput(i) = cnvrt(i)
    end do


    ! call convolve(nb, nx, ny, bb, xx, yy)


!        print*, sm1, sm2
        ! ! COMPUTE EVAPORATION
        ! select case(ETOP)
        !     ! several different case options
        !     case(0);  E = ETa(PET(i), sm1, sm1max)  ! this is just a function, not a subroutine
        ! end select


        ! ! COMPUTE PERCOLATION
        ! select case(DRAINOP)
        !     case(0); q12 = q12a(sm1,  sm1max, ku, c)
        !     case(2); q12 = q12c(sm1F, sm1Fmax, q0, sm2, sm2max, psi, alpha)
        ! end select


        ! ! COMPUTE BASEFLOW
        ! select case(BASEOP)
        !     case(0);  qb = Qba(sm2, ks)
        !     case(1);  qb = Qbb(sm2, sm2max, ks, lowercasen)
        ! end select

        ! ! COMPUTE SURFACE RUNOFF  !
        ! ! Select the method for computing the saturated area
        ! select case(SATOP)
        !     case(0); Ac = satarea_a(sm1, sm1max, beta)
        ! end select

        ! !Now compute the overland flow
        ! qsx = Ac*qin

        ! ! UPDATE SOIL MOISTURE TOP LAYER
        ! select case(SMOP1)
        !     case(0); call s1_topvic(qin, E, qsx, qif, q12, qufof, sm1)
        ! end select

        ! ! UPDATE SOIL MOISTURE BOTTOM LAYER
        ! select case(SMOP2)
        !     case(0); call s2_topprms(q12, qb, sm2)
        ! end select





!  print*, params%nmf

end subroutine model_driver


end module driver





