module snowmodule17
implicit none
contains


!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 1) Snow-17 accumulation and ablation model.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Adapted from
! This version of Snow-17 is intended for use at a point location.
! Based on Anderson (2006) and Mark Raleigh's matlab code.
! Primary Citations:
! 1.  Anderson, E. A. (1973), National Weather Service River Forecast System
!     Snow   Accumulation   and   Ablation   Model,   NOAA   Tech.   Memo.   NWS
!     HYDro-17, 217 pp., U.S. Dep. of Commer., Silver Spring, Md.
! 2.  Anderson, E. A. (1976), A point energy and mass balance model of a snow
!     cover, NOAA Tech. Rep. 19, 150 pp., U.S. Dep. of Commer., Silver Spring, Md.
!
! Python code Written by Joe Hamman April, 2013
!
! Translated into Fortran (lol why) by William Rudisill, March 2021
!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

real function melt_function(jday,  &
	                          dt,    &
	                          mfmax, &
	                          mfmin)
    implicit none

    ! input
    real, intent(in) :: jday
    real, intent(in) :: dt
    real, intent(in) :: mfmax
    real, intent(in) :: mfmin

    ! internal
    integer :: n_mar21
    integer :: days
    real :: sv
    real :: av
    ! BEGIN
    ! ------

    ! a few parameters
    days = 365
    n_mar21 = jday - 80
    sv = (0.5 * SIN((n_mar21 * 2 * 3.14159265359) / days)) + 0.5

    ! assumes lat < 45
    av = 1.0

    ! do something different for northern latitudes... NOT IMPLEMENTED
    ! av = <code here>
    ! ....
    melt_function = (dt / 6) * ((sv * av * (mfmax - mfmin)) + mfmin)
end function melt_function



! Begin subroutines
subroutine snow17(jday, &
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
                  swe, &
                  ait, &
                  w_qx, &
                  w_q, &
                  w_i, &
                  deficit, &
                  melt_outflow)

    implicit none

    ! input
    real, intent(in) :: jday
    real, intent(in) :: precip
    real, intent(in) :: tair
    real, intent(in) :: elevation
    real, intent(in) :: dt
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

    ! inouts (snow states. these can change)
    real, intent(inout) :: swe
    real, intent(inout) :: ait
    real, intent(inout) :: w_qx
    real, intent(inout) :: w_q
    real, intent(inout) :: w_i
    real, intent(inout) :: deficit
    real, intent(out) :: melt_outflow


    !!!!! Declare f2py outputs... inout will not return anything in python
    !!!!! ----------------------------------------------------------------
    !f2py intent(in,out) :: swe
    !f2py intent(in,out) :: ait
    !f2py intent(in,out) :: w_qx
    !f2py intent(in,out) :: w_q
    !f2py intent(in,out) :: w_i
    !f2py intent(in,out) :: deficit

    ! internal
    real, dimension(2) :: transitionx
    real, dimension(2) :: transitiony
    real :: stefan        ! constants
    real :: p_atm
    real :: mf
    real :: tipm_dt
    real :: pn
    real :: rain
    real :: delta_hd_snow
    real :: fracrain
    real :: fracsnow
    real :: t_snow_new
    real :: t_rain
    real :: e_sat
    real :: delta_hd_t
    real :: m_ros
    real :: m_ros1
    real :: m_ros2
    real :: m_ros3
    real :: m_nr
    real :: qw
    !f2py intent(in,out) :: e_sat


    ! BEGIN
    ! ------


    ! constant
    !stefan = 6.12 * (10 ** (-10))
    stefan=6.12E-10

    ! atmospheric pressure (mb) where elevation is in HUNDREDS of meters
    ! (this is incorrectly stated in the manual)
    p_atm = 33.86 * (29.9 - (0.335 * elevation / 100) + (0.00022 * ((elevation / 100) ** 2.4)))

    ! what's happening here?
    tipm_dt = 1.0 - ((1.0 - tipm) ** (dt / 6))

    ! Loop over time dimension
    mf = melt_function(jday, dt, mfmax, mfmin)

    !!!air temperature at this time step (deg C)


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  Select the 'rain versus snow' method !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Option 0
    if (rvs == 0) then
        if (tair <= pxtemp) then
            ! then the air temperature is cold enough for snow to occur
            fracsnow = 1.0
        else
            fracsnow = 0.0
        end if

    ! Option 1
    else if  (rvs == 1) then
         if (tair <= pxtemp1) then
            fracsnow = 1.0
         else if (tair >= pxtemp2) then
             fracsnow = 0.0
         else
             ! Linear interpolate between 0 and 1
             fracsnow  = (tair - pxtemp1)/(pxtemp2 - pxtemp1)

         end if

    ! Option 2
    else if (rvs == 2) then
         fracsnow = 1.0
    else
    !     raise ValueError('Invalid rain vs snow option')
    end if

   ! now move on...
   fracrain = 1.0 - fracsnow

   ! snow
   pn = precip * fracsnow !* scf -- get rid of this part. undercatch factor

   ! not sure what these are (yet)
   w_i =  w_i + pn
   melt_outflow = 0.0

   ! amount of precip (mm) that is rain during this time step
   rain = fracrain * precip

   ! Temperature and Heat deficit from new Snow
   ! The new snow temperature will be, at a minimum, zero
   if (tair < 0.0) then
       t_snow_new = tair
       ! delta_hd_snow = change in the heat deficit due to snowfall (mm)
       delta_hd_snow = - (t_snow_new * pn) / (80 / 0.5)
       t_rain = pxtemp
   else
       t_snow_new = 0.0
       delta_hd_snow = 0.0
       t_rain = tair
   end if


  ! Antecedent temperature Index
  if (pn > (1.5 * dt)) then
      ait = t_snow_new
  else
      ait = ait + tipm_dt * (tair - ait)
  end if

  if (ait > 0) then
      ait = 0
  end if

  ! heat deficit
  delta_hd_t = nmf * (dt / 6.0) * ((mf) / mfmax) * (ait - t_snow_new)

  ! saturated vapor pressure at tair (mb)
  e_sat = 2.7489 * (10 ** 8) * EXP((-4278.63 / (tair + 242.792)))
  ! Rain-on-snow melt_outflow
  if (rain  > (0.25 * dt)) then
    !melt_outflow (mm) during rain-on-snow periods is:
    m_ros1 = MAX(stefan * dt * (((tair + 273.0) ** 4) - (273.0 ** 4)), 0.0)
    m_ros2 = MAX(0.0125 * rain * t_rain, 0.0)
    m_ros3 = MAX(8.5 * uadj * dt / 6.0 * ((0.9 * e_sat) - 6.11) +  (0.00057 * p_atm * tair), 0.0)
    m_ros = m_ros1 + m_ros2 + m_ros3
  else
    m_ros = 0.0
  end if

  ! Non-Rain melt_outflow
  if (rain <= (0.25 * dt).and.(tair > mbase)) then
      ! melt_outflow during non-rain periods is:
      m_nr = (mf * (tair - mbase)) + (0.0125 * rain * t_rain)
  else
      m_nr = 0.0
  end if


  ! Ripeness of the snow cover
  melt_outflow = m_ros + m_nr
  if (melt_outflow <= 0) then
      melt_outflow = 0.0
  end if

  if (melt_outflow < w_i) then
      w_i = w_i - melt_outflow
  else
      melt_outflow = w_i + w_q
      w_i = 0.0
  end if

  !qw = liquid water available melt_outflowed/rained at the snow surface (mm)
  qw = melt_outflow + rain
  ! w_qx = liquid water capacity (mm)
  w_qx = plwhc * w_i
  ! deficit = heat deficit (mm)
  deficit = deficit + delta_hd_snow + delta_hd_t

  ! limits of heat deficit
  if (deficit < 0) then
      deficit = 0.0
  else if (deficit > (0.33 * w_i)) then
      deficit = 0.33 * w_i
  end if

  !!!
  ! Snow cover is ripe when both (deficit=0) & (w_q = w_qx)
  if (w_i > 0.0) then
      if ((qw + w_q) > ((deficit * (1 + plwhc)) + w_qx)) then
          ! # THEN the snow is RIPE
          ! # Excess liquid water (mm)
          melt_outflow = qw + w_q - w_qx - (deficit * (1 + plwhc))
          ! # fills liquid water capacity
          w_q = w_qx
          ! # w_i increases because water refreezes as heat deficit is
          ! # decreased
          w_i = w_i + deficit
          deficit = 0.0
      else if ((qw >= deficit).and.( (qw + w_q) <= (deficit * (1 + plwhc) + w_qx))) then
          ! # THEN the snow is NOT yet ripe, but ice is being melt_outflowed
          melt_outflow = 0.0
          w_q = w_q + qw - deficit
          ! # w_i increases because water refreezes as heat deficit is
          ! # decreased
          w_i = w_i + deficit
          deficit = 0.0

      else
          ! # (qw < deficit) %elseif ((qw + w_q) < deficit):
          ! # THEN the snow is NOT yet ripe
          melt_outflow = 0.0
          ! # w_i increases because water refreezes as heat deficit is
          ! # decreased
          w_i = w_i + qw
          deficit = deficit - qw

      end if
      ! now update swe
      swe = w_i + w_q
  ! if there is no water input (w_i)
  else
      melt_outflow = qw
      swe = 0
  end if

!  deficit...
  if (deficit == 0) then
      ait = 0
  end if

  ! End of model execution
  !model_swe = swe  ! total swe (mm) at this time step
  !melt_outflow = e


end subroutine snow17



subroutine snow17driver(ntimes, jdayVec, precipVec, tairVec, &  ! INPUTS
                        nlayers, opg, dz,           & ! INPUTS
                        dt,   &                                 ! parameters
                        rvs, &                                  ! parameters
                        uadj, &                                 ! parameters
                        mbase, &                                ! parameters
                        mfmax, &                                ! parameters
                        mfmin, &                                ! parameters
                        tipm, &                                 ! parameters
                        nmf, &                                  ! parameters
                        plwhc, &                                ! parameters
                        pxtemp, &                               ! parameters
                        pxtemp1, &                              ! parameters
                        pxtemp2, &                              ! parameters
                        melt_outflowVec, &                      ! OUTPUT
                        sweVec)                                 !  OUTPUT

                        implicit none

                        ! INPUT forcings
                        integer, intent(in) :: ntimes        ! Number of timesteps..
                        real, intent(in), dimension(ntimes) :: jdayVec   ! julian day
                        real, intent(in), dimension(ntimes) :: precipVec ! precipitation
                        real, intent(in), dimension(ntimes) :: tairVec   ! air temperature

                        ! INPUT paramters  -- preciptiation adjustment
                        real, intent(in) :: opg
                        integer, intent(in) :: nlayers       ! Number of model layers
                        real, intent(in), dimension(nlayers) :: dz       ! vector of *mean* elevtaion differences for each nlayers


                        ! INPUT paramters -- Snow 17
                        real, intent(in) :: dt
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

                        ! OUTPUTS
                        real, intent(out), dimension(nlayers,ntimes) :: melt_outflowVec
                        real, intent(out), dimension(nlayers,ntimes) :: sweVec

                        !internal

                        ! states
                        real, parameter :: t_lapse = -.0065  ! 6.5 c/m
                        real, parameter :: stat_elev = 35.  ! elevation in HUNDREDS of meters. Approx elevation of the Butte Snotel
                        real  :: swe
                        real  :: ait
                        real  :: w_qx
                        real  :: w_q
                        real  :: w_i
                        real  :: deficit
                        real  :: melt
                        real  :: elevation
                        real :: precip_il ! precipitation in the lth layer
                        real :: temp_il
                        ! misc
                        integer :: i   ! for time loop
                        integer :: l   ! for layer loop

                        ! loop through layers
                        do l=1,nlayers

                          ! Set Initial conditions for the layer
                          ait = 0.
                          swe = 0.
                          ait = 0.
                          w_qx = 0.
                          w_q = 0.
                          w_i = 0.
                          deficit = 0.
                          melt = 0

                          ! midpoint elevation ...
                          elevation = stat_elev + dz(l)/100 ! elevation is in units of 100s of meters for whatever reason...

                          ! loop through timesteps
                          do i=1,ntimes

                            ! adjust the precipitaiton
                            precip_il = precipVec(i) + opg * dz(l) * precipVec(i)
                            temp_il = tairVec(i) + t_lapse * dz(l)

                            call snow17(jdayVec(i), &
                                        precip_il, &
                                        temp_il, &
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
                                        swe, &
                                        ait, &
                                        w_qx, &
                                        w_q, &
                                        w_i, &
                                        deficit, &
                                        melt)

                            ! store the data...
                            melt_outflowVec(l,i) = swe
                            sweVec(l,i) = melt
                        end do
                      end do


end subroutine snow17driver



end module snowmodule17





