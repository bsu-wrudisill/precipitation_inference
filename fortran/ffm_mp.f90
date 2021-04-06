module fast_foward_model
use hydrofunctions
use snowmodule17
use snowmodule01
implicit none
!integer, parameter :: grid_size = 100




contains


!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Function definitions
!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine one_step(Nls,&
                    Pd,&
                    Td,&
                    T30,&
                    L,&
                    Wu,&
                    Wb,&
                    Snow,&
                    Q,&
                    frtdir,&
                    frtgw,&
                    smcap,&
                    etpar,&
                    t_snow,&
                    t_melt,&
                    t_base,&
                    t_power)
  implicit none
  integer, intent(in)  :: Nls

  ! Model Parameters
  real, intent(in)     :: frtdir
  real, intent(in)     :: frtgw
  real, intent(in)     :: smcap
  real, intent(in)     :: etpar
  real, intent(in)     :: t_snow
  real, intent(in)     :: t_melt
  real, intent(in)     :: t_base
  real, intent(in)     :: t_power

  ! Forcing variables
  real, intent(in)                     :: L
  real, intent(in), dimension(0:Nls-1) :: Pd
  real, intent(in), dimension(0:Nls-1) :: Td
  real, intent(in), dimension(0:Nls-1) :: T30

  ! Model States
  real, intent(inout), dimension(0:Nls-1) :: Snow
  real, intent(inout)  :: Wu
  real, intent(inout)  :: Wb
  real, intent(out) :: Q

  ! Dummy variables
  real :: Prk,
  real :: Psk
  real :: Mdk
  real :: Prsum
  real :: Mdsum
  real :: Tdavg
  real :: Qd
  real :: Qb
  real :: Wi
  real :: D
  real :: E
  real :: ddpar_value
  real :: pet_value
  real :: dWudt
  real :: dWbdt
  integer :: k

  ! Holders for layer calcs
  Mdsum = 0
  Prsum = 0
  Tdavg = 0

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!                 ~~~BEGIN~~~                    !!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! loop through the snow layers
  do k=0,Nls-1

    ! call snow01(precip, tair, tair30, pxtemp, t_melt, t_base, t_power, Snow)

    ! ! determine precip phase
    ! if (Td(k) > t_snow) then
    !     Prk = Pd(k)
    ! else
    !     Prk = 0
    ! end if
    ! ! create snow and rain values
    ! Psk = Pd(k) - Prk

    ! ! Compute melt factors for each layer
    ! ddpar_value = DDPar(T30(k), t_base, t_power)

    ! ! Compute melt
    ! Mdk = Ms(Td(k), Snow(k), t_melt, ddpar_value)

    ! ! Gather up the total melt
    ! Mdsum = Mdk + Mdsum
    ! Prsum = Prsum + Prk

    ! ! melt/accumulate the snowpack
    ! Snow(k) = Snow(k) + Psk - Mdk ! last snow, new snow, snow melt
    ! Tdavg = Tdavg + Td(k)
    print*, k

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




subroutine fwd(N
               Ndz,&
               Nls,&
               PdVec,&
               TdVec,&
               LVec,&
               QVec,&
               dz,&
               frtdir,&
               frtgw,&
               smcap,&
               etpar,&
               t_snow,&
               t_melt,&
               t_base,&
               t_power,&
               bias,&
               opg)!, opg, mult)&

  implicit none
  integer, intent(inout) :: N                ! Number of timesteps
  integer, intent(inout) :: Ndz              ! Length of the elevation (dz) vector
  integer, intent(inout) :: Nls              ! Number of layers
  real, intent(in), dimension(0:N) :: PdVec  ! Precipitation vector
  real, intent(in), dimension(0:N) :: TdVec  ! Temperature vector
  real, intent(in), dimension(0:N) :: LVec   ! Day Length vector (for PET)
  real, intent(in), dimension(0:Ndz) :: dz   ! Elevation vector (length Ndz)
  real, intent(out), dimension(0:N):: QVec   ! Channel discharge vector

  real, intent(in) :: frtdir
  real, intent(in) :: frtgw
  real, intent(in) :: smcap
  real, intent(in) :: etpar
  real, intent(in) :: t_snow
  real, intent(in) :: t_melt
  real, intent(in) :: t_base
  real, intent(in) :: t_power
  real, intent(in) :: bias
  real, intent(in) :: opg

  ! internal ?
  real, dimension(0:Nls-1) :: Pd
  real, dimension(0:Nls-1) :: Td
  real, dimension(0:Nls-1) :: T30
  real, dimension(0:Nls-1) :: Snow
  real :: dPd
  real :: dTd
  real :: dT30
  real :: L
  real :: Wu
  real :: Wb
  real :: Q
  integer :: i
  integer :: maxi
  integer :: j
  integer :: k
  integer :: ll
  integer :: llr
  integer :: llr_0
  integer :: llc
  integer :: llc_0
  integer :: li
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!                 ~~~BEGIN~~~                    !!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
           dPd = max(0.0, dPd + PdVec(i) * dz(k) * opg + bias)
         end if

      ! End looping through the dz points...
      end do

      ! compute average adjustment
      dTd = dTd/llc_0
      dPd = dPd/llc_0

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





