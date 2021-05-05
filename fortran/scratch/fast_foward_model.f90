module fast_foward_model
!program testing
implicit none
!integer, parameter :: grid_size = 100




contains


!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Function definitions
!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

real function Drainage(Wu, smcap)
    real, intent(in) :: Wu, smcap
    if(Wu.gt.smcap) then
        Drainage = Wu - smcap
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine subr1(Q)
  implicit none
  real, intent(inout) :: Q
  real :: Wu, Wb, Snow

  ! Setup the initial conditions
   Wu = 1.0
   Wb = 1.0
   Snow = 50.

!   call one_step(0.0, 10.0, 5.0, 8., Wu, Wb, Snow, Q, .1, .1, 200., .01, 0., 0., 1., 1.)
   ! print*, Q
end subroutine subr1


!Pd, Td, T30, L, Q, Wu, Wb, Snow, frtdir, frtgw, smcap, etpar, t_snow, t_melt, t_base, t_power)

subroutine one_step(Nls, Pd, Td, T30, L, Wu, Wb, Snow, Q, frtdir, frtgw, smcap, etpar, t_snow, t_melt, t_base, t_power)
  implicit none
  integer, intent(in)  :: Nls
  real, intent(in)     :: L, frtdir, frtgw, smcap, etpar, t_snow, t_melt, t_base, t_power    !Parameters
  real, intent(in), dimension(0:Nls-1)    :: Pd, Td, T30                                    !Forcings
  real, intent(inout), dimension(0:Nls-1) :: Snow
  real, intent(inout)  :: Wu, Wb
  real, intent(out) :: Q

  ! Dummy variables
  real :: Prk, Psk, Mdk, Prsum, Mdsum, Tdavg, Qd, Qb, Wi, D, E, ddpar_value, pet_value, dWudt, dWbdt
  integer :: k

 ! Holders for layer calcs
  Mdsum = 0
  Prsum = 0
  Tdavg = 0

  ! loop through the snow layers
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
  D = Drainage(Wu, smcap)

  ! Compute baseflow (bottom layer to runoff)
  Qb = Wb * frtgw

  ! Compute PET
  ! First compute saturation specific humidity
  pet_value = PET(L, Tdavg, etpar)

  ! Compute AET
  E = pet_value*(Wu/smcap)

  ! Compute change in the water balance of top layer
  dWudt = Wi - D - E

  ! Compute change in water balance at bottom layer
  dWbdt = D - Qb

  ! update values
  ! soil moisture
  Wb = Wb + dWbdt
  Wu = Wu + dWudt

  ! Compute discharge
  Q = Qb + Qd

end subroutine one_step

! subroutine vectest(N,a)
!   implicit none
!   integer, intent(inout) :: N
!   real, intent(out) :: a(N), K
!   integer :: i
!   do i = 1,N
!      a(i) = 1.0
!   end do
!   K = 4
! end subroutine vectest



subroutine fwd(N, Ndz, Nls, PdVec, TdVec, LVec, QVec, dz, frtdir, frtgw, smcap, etpar, t_snow, t_melt, t_base, t_power, bias, opg)!, opg, mult)
  implicit none
  integer, intent(inout) :: N, Ndz, Nls   ! length of the Forcings
  real, intent(in), dimension(0:N) :: PdVec, TdVec, LVec                               !Forcings
  real, intent(in), dimension(0:Ndz) :: dz                                             !Forcings
  real, intent(out), dimension(0:N):: QVec                                             !Q
  real, intent(in) :: frtdir, frtgw, smcap, etpar, t_snow, t_melt, t_base, t_power, bias, opg

  ! internal ?
  real, dimension(0:Nls-1) :: Pd, Td, T30, Snow
  real :: dPd, dTd, dT30, L, Wu, Wb, Q
  integer :: i, maxi, j, k, ll, llr, llr_0, llc, llc_0, li


  Wb = 0.0
  Wu = 0.0
  Q = 0.0
  Snow = 0.0

  ! Get the layer lenthg
  ll = Ndz/Nls          ! The 'layer length'
  llr = mod(Ndz, Nls)   ! The 'layer length' remainder

  !--------------------------------
  ! Loop through the time dimension
  !--------------------------------
  do i = 1,N
    ! get length of day..
    L = LVec(i)

    ! Setup the layer counters/indices ...
    llc = 0
    li = 0
    llr_0 = llr

    ! Initialize the arrays for forcings ...
    Pd = 0
    Td = 0
    T30 = 0

    ! Compute the 30-day temperature average
    dT30 = 0.0
    maxi = max(i-30, 0)

    do j = maxi,i
       dT30 = dT30 + TdVec(j)
    end do


    !----------------------------------------
    ! Loop through the model layers dimension
    !----------------------------------------
    do while (llc < Ndz)

      ! zero out the remainder term
      llc_0 = llc+ll+llr_0

      ! these get reset for each layer. they are scalars
      dTd = 0.0
      dPd = 0.0


      ! --------------------------
      ! Loop through the dz points
      ! --------------------------
      do k=llc, llc_0  ! get the starting allnd ending indices
         dTd = dTd + dz(k) * -0.0065

         ! Add bias only if there is precipitation
         if (PdVec(i).gt.0.0) then
           ! Also  make sture that it's never non-negative...
           dPd = dPd + PdVec(i) * dz(k) * opg + bias
         end if

      ! End looping through the dz points...
      end do

      ! compute average adjustment
      dTd = dTd/(llc_0-llc)
      dPd = dPd/(llc_0-llc)

      ! Note the dimensions here...
      Td(li) = TdVec(i) + dTd
      T30(li) = TdVec(i) + dTd
      Pd(li) = (PdVec(i) + dPd)/Nls


      li = li + 1
      llc = llc_0
      llr_0 = 0
    !----------------------------------------
    ! END Loop through the model layers dimension
    !----------------------------------------
    end do


    call one_step(Nls, Pd, Td, T30, L, Wu, Wb, Snow, Q, frtdir, frtgw, smcap, etpar, t_snow, t_melt, t_base, t_power)
    QVec(i) = Q


  end do

end subroutine fwd


end module fast_foward_model
!end program testing





