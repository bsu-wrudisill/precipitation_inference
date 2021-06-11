module Hymod
!use parameters
implicit none

contains

! Helper functions
function PET(L, t)
    real :: PET, L, t
    real :: sat_vap_pres, abs_humidity
    ! Potential evapotranspiration calculator
    ! L: Length of day
    ! rho_v_sat: staturated absolute humidity
    ! etpar: parameter
    sat_vap_pres = (.6112)*exp((17.67*t)/(t + 243.5))
    abs_humidity = (2165 * sat_vap_pres)/(t + 273.15)
    PET = .55 * (L/12.)**2 * abs_humidity
end function PET



subroutine nash(k, inp, xbeg, N, xend, out)
    real,    intent(in)                 :: k         ! flow rate limiter
    real,    intent(in)                 :: inp       ! input value -- note dimension
    real,    intent(in),  dimension(N)  :: xbeg      ! beginning of ... tank
    integer, intent(in)                 :: N         !
    real,    intent(out), dimension(N)  :: xend      ! end of ... tank
    real,    intent(out)                :: out

    ! internal
    real, dimension(N) :: oo
    integer :: i

    ! oo = 0.
    ! xend = 0.

    do i=1,N
        oo(i) = k*xbeg(i)
        xend(i) = xbeg(i) - oo(i)

        if (i == 1) then
            xend(i) = xend(i) + inp
        else
            xend(i) = xend(i) + oo(i-1)
        end if
    out = oo(N)

    end do
end subroutine nash


subroutine Pdm01(huz, b, hbeg, pp, PET, &   ! inputs
                 ov, ET, hend, cend)            ! outputs
    !input
    real, intent(in) :: huz
    real, intent(in) :: b
    real, intent(in) :: hbeg
    real, intent(in) :: pp
    real, intent(in) :: PET       ! parameters

    ! output
    real, intent(out) :: ov
    real, intent(out) :: ET
    real, intent(out) :: hend
    real, intent(out) :: cend

    ! internal
    real :: bv
    real :: cpar
    real :: cbeg
    real :: ov1
    real :: ov2
    real :: ppinf
    real :: hint
    real :: cint


    ! Begin
    bv = LOG(1.0-MIN(1.9999,b)/2.0)/LOG(0.5)
    cpar = huz/(1+bv)
    cbeg = cpar*(1-(1-hbeg/huz)**(1+bv))

    ! overland flow...
    ov2 = MAX(pp+hbeg-huz, 0.)

    ! Infiltraion
    ppinf = pp-ov2

    ! something
    hint = MIN((ppinf+hbeg), huz)
    cint = cpar*(1.-(1.-hint/huz)**(1.+bv))
    ov1 = MAX(ppinf+cbeg-cint, 0.)

    ! these are the outputs
    ov = ov1+ov2
    ET = MIN(PET, cint)
    cend = cint-ET
    hend = huz*(1.-(1.-cend/cpar)**(1./(1.+bv)))

end subroutine Pdm01



subroutine Hymod01(QINvec, &
	               PETvec, &
				   ntimes, &
                   Nq, &
                   Kq, &
                   Ks, &
                   Alp,&
                   Huz,&
                   B, &
                   q)

    real, intent(in), dimension(ntimes) :: QINvec ! change this later. for testing purposes
    real, intent(in), dimension(ntimes) :: PETvec ! change this later. for testing purposes
    integer, intent(in) :: ntimes
    integer, intent(in) :: Nq   ! parameters
    real, intent(in) :: Kq      ! parameters
    real, intent(in) :: Ks      ! parameters
    real, intent(in) :: Alp     ! parameters
    real, intent(in) :: Huz     ! parameters
    real, intent(in) :: B       ! parameters

    ! output
    real, dimension(ntimes), intent(out) :: q

    ! internal states/fluxes
    real, dimension(ntimes) :: xhuz
    real, dimension(ntimes) :: xs
    real, dimension(ntimes, Nq) :: xq
    real, dimension(ntimes) :: et
    real, dimension(ntimes) :: ov
    real, dimension(ntimes) :: xcuz
    real, dimension(ntimes) :: qq
    real, dimension(ntimes) :: qs
    real :: alpha_ov, alpha_ov_m1

    ! internal counters
    integer :: i

    ! initialize arrays
    xhuz = 0.
    xs = 1.0
    xq = 1.0
    et = 0.
    ov = 0.
    xcuz = 0.
    qq = 0.
    qs = 0.

    ! Begin
    do i=1,ntimes

        call Pdm01(Huz, B, xhuz(i), QINvec(i), PETvec(i), & ! inputs
                   ov(i), et(i), xhuz(i), xcuz(i))          ! outputs

        ! call the nash cascade thing
        alpha_ov = Alp * ov(i)
        alpha_ov_m1 = (1.- Alp) * ov(i)

        ! do the quick flow

!nash(k, inp, xbeg, N, xend, out)

        call nash(kq, alpha_ov, xq(i, :), Nq, xq(i, :), qq(i))

        ! now the slow flow
        call nash(ks, alpha_ov_m1, xs(i), 1, xs(i), qs(i))

        ! now...
        if (i < ntimes) then
            xhuz(i+1) = xhuz(i)
            xq(i+1,:) = xq(i, :)
            xs(i+1) = xs(i)
        end if

        ! compute the total
        q(i) = qs(i) + qq(i)
	end do


end subroutine Hymod01


subroutine hymod_driver(ntimes,     &         ! FORCING   Number of model timesteps
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
                        Nq,         &         ! PARAMETER   SNOW17
                        Kq,         &         ! PARAMETER   SNOW17
                        Ks,         &         ! PARAMETER   SNOW17
                        Alp,        &         ! PARAMETER   SNOW17
                        Huz,        &         ! PARAMETER   SNOW17
                        B,          &         ! PARAMETER  SNOW17
                        q) ! OUTPUT


    ! use the snowmod code
    use snowmodule17

    ! stuff
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
    ! output
    real, intent(out), dimension(ntimes) :: q

    ! hymod params
    integer, intent(in) :: Nq   ! parameters
    real, intent(in) :: Kq      ! parameters
    real, intent(in) :: Ks      ! parameters
    real, intent(in) :: Alp     ! parameters
    real, intent(in) :: Huz     ! parameters
    real, intent(in) :: B       ! parameters
    ! internal
    real, dimension(nlayers,ntimes) :: outflowVec
    real, dimension(nlayers,ntimes) :: sweVec
    real, dimension(nlayers,ntimes) :: rainVec
    real, dimension(nlayers,ntimes) :: ptotVec
    real, dimension(ntimes) :: outflowVec_adjusted
    real :: scale

    ! BEGIN
    ! call snow17
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



    ! Do some flow adjusting ...
    scale = (SUM(outflowVec)/real(nlayers))/SUM(outflowVec(3,:))
    print*, scale
    outflowVec_adjusted = outflowVec(3,:)*scale

    ! Now call the hymod model ...
    call Hymod01(outflowVec_adjusted, &
                 PET, &
                 ntimes, &
                 Nq, &
                 Kq, &
                 Ks, &
                 Alp,&
                 Huz,&
                 B,  &
                 q)



end subroutine hymod_driver
end module Hymod
