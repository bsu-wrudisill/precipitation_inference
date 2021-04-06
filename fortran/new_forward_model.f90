module fast_foward_model
implicit none




contains


!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Function definitions
!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

real function Drainage(Wu, S1max)
    real, intent(in) :: Wu, S1max
    if(Wu.gt.S1max) then
        Drainage = Wu - S1max
    else
        Drainage = 0
    end if
end function Drainage


real function PET(L, t, etpar)
    real, intent(in) :: L, t, etpar
    real :: sat_vap_pres, abs_humidity
    ! Potential evapotranspiration calculator
    ! L: Length of day
    ! rho_v_sat: staturated absolute humidity
    ! etpar: parameter
    sat_vap_pres = (.6112)*exp((17.67*t)/(t + 243.5))
    abs_humidity = (2165 * sat_vap_pres)/(t + 273.15)
    PET=L*abs_humidity*etpar
end function PET


! Snow Melting Functions
real function DDPar(T30, t_base, t_power)
    real, intent(in) :: T30, t_base, t_power

    ! T30: previous 30day mean temp
    ! t_base: parameters
    ! t_power: parameters

    if (T30.gt.t_base) then
        DDPar = (T30 - t_base)**t_power
    else
        DDPar = t_base
    end if
end function DDPar


! Snow melt function
real function Ms(Td, Snow, t_melt, ddpar)
    real, intent(in) :: Td, Snow, t_melt, ddpar
    real :: melt

    ! # Compute the snowmelt
    ! # Td: temperature
    ! # Snow: amount in the bucket
    ! # ddpar: parameter (but not constant)
    ! # t_melt: parameter

    if (Td.gt.t_melt) then
        melt = (Td-t_melt)*ddpar

        if (melt.gt.Snow) then
            Ms = Snow

        else
            Ms = melt
        end if
    else
        Ms = 0
    end if
end function Ms


! Baseflow functions
real function Baseflow01(frtgw, S2)
    real, intent(in) :: frtgw   ! parameter
    real, intent(in) :: S2      ! storage in the lower soil layer
    real, intent(out) :: Qb
    Qb = frtgw * S2
end function Baseflow01



real function Baseflow01(S2, n)
    real, intent(in) :: frtgw   ! parameter
    real, intent(in) :: S2      ! storage in the lower soil layer
    real, intent(out) :: Qb
    Qb = frtgw * S2
end function Baseflow01



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! real function routing(N, Q, alpha, beta)
!   integer, intent(in) :: N                ! Len of Q
!   real, intent(in), dimension(0:N) :: Q   ! Input timeseries of Q
!   real, intent(in) :: alpha, beta         ! Parameters
!   real, intent(out) :: Qr                 ! 'Routed' Q




subroutine Snow(k, iPd, iTd, iT30, iSnow)
integer, intent(in)                      :: k               ! num of snow layers
real,    intent(in),    dimension(0:k-1) :: iPd, iTd, iT30  ! daily values of t, p, t30
real,    intent(inout), dimension(0:k-1) :: Snow            ! snow pack





end subroutine Snow






subroutine one_step(Nls, Pd, Td, T30, L, &
                    S1, S2, Snow, Q, E, PET_val, &
                    frtdir, frtgw, S1max, etpar, t_snow, t_melt, t_base, t_power)


  implicit none
  integer, intent(in)  :: Nls
  real, intent(in)     :: L, frtdir, frtgw, S1max, etpar, t_snow, t_melt, t_base, t_power    !Parameters
  real, intent(in), dimension(0:Nls-1)    :: Pd, Td, T30                                     !Forcings
  real, intent(inout), dimension(0:Nls-1) :: Snow
  real, intent(inout)  :: S1, S2, E, PET_val
  real, intent(out) :: Q

  ! Dummy variables
  real :: Prk, Psk, Mdk, Prsum, Mdsum, Tdavg, Qd, Qb, Wi, D, ddpar_value, dS1dt, dS2dt
  integer :: k

 ! Holders for layer calcs
  Mdsum = 0
  Prsum = 0
  Tdavg = 0

  do k=0,Nls-1
    ! determine precip phase
    if (Td(k) > t_snow) then
        Prk = Pd(k)
    else
        Prk = 0
    end if
    ! create snow and rain values
    Psk = Pd(k) - Prk

    ! Compute melt factors for each layer
    ddpar_value = DDPar(T30(k), t_base, t_power)

    ! Compute melt
    Mdk = Ms(Td(k), Snow(k), t_melt, ddpar_value)

    ! Gather up the total melt
    Mdsum = Mdk + Mdsum
    Prsum = Prsum + Prk

    ! melt/accumulate the snowpack
    Snow(k) = Snow(k) + Psk - Mdk ! last snow, new snow, snow melt
    Tdavg = Tdavg + Td(k)

  end do

  !
  Tdavg = Tdavg/Nls


  ! Compute the direct runoff (Qd)
  Qd = frtdir * (Prsum + Mdsum)

  ! Compute the Soil Water input
  Wi = (1-frtdir) * (Prsum + Mdsum)

  ! Compute drainage; 1st layer of soil --> 2nd layer of soil
  D = Drainage(S1, S1max)

  ! Compute baseflow (bottom layer to runoff)
  Qb = Baseflow(S2, frtgw)

  ! Compute PET
  ! First compute saturation specific humidity
  PET_val = PET(L, Tdavg, etpar)

  ! Compute AET
  E = PET_val*(S1/S1max)

  ! Compute change in the water balance of top layer
  dS1dt = Wi - D - E

  ! Compute change in water balance at bottom layer
  dS2dt = D - Qb

  ! update values
  ! soil moisture
  S2 = S2 + dS2dt
  S1 = S1 + dS1dt

  ! Compute discharge
  Q = Qb + Qd

end subroutine one_step




end module fast_foward_model
!end program testing





