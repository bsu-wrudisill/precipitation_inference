module fast_foward_model
implicit none
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


real function Ms(Td, Snow, tmelt, ddpar)
    real, intent(in) :: Td, Snow, tmelt, ddpar
    real :: melt

    ! # Compute the snowmelt
    ! # Td: temperature
    ! # Snow: amount in the bucket
    ! # ddpar: parameter (but not constant)
    ! # tmelt: parameter

    if (Td.gt.tmelt) then
        melt = (Td-tmelt)*ddpar

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

   call one_step(0.0, 10.0, 5.0, 8., Wu, Wb, Snow, Q, .1, .1, 200., .01, 0., 0., 1., 1., 1.)
   ! print*, Q
end subroutine subr1


!Pd, Td, T30, L, Q, Wu, Wb, Snow, frtdir, frtgw, smcap, etpar, tmelt, t_snow, t_melt, t_base, t_power)

subroutine one_step(Pd, Td, T30, L, Wu, Wb, Snow, Q, frtdir, frtgw, smcap, etpar, tmelt, t_snow, t_melt, t_base, t_power)
   implicit none
   real, intent(in) :: Pd, Td, T30, L                                                         !Forcings
   real, intent(in) :: frtdir, frtgw, smcap, etpar, tmelt, t_snow, t_melt, t_base, t_power    !Parameters
   real, intent(inout) :: Wu, Wb, Snow
   real, intent(out) :: Q

   ! internal ?
   real :: Pr, Ps, Md, Qd, Qb, Wi, D, E, ddpar_value, pet_value, dWudt, dWbdt


   if (Td > t_snow) then
       Pr = Pd
   else
       Pr = 0
   end if

    ! ! # compute daily snow accumulation
    Ps = Pd - Pr

    ! ! # Compute the snowmelt
    ! ! # compute ddpar first
    ddpar_value = DDPar(T30, t_base, t_power)
    Md = Ms(Td, Snow, tmelt, ddpar_value)

    ! ! # Compute the direct runoff (Qd)
    Qd = frtdir * (Pr + Md)

    ! ! # Compute the Soil Water input
    Wi = (1-frtdir) * (Pr + Md)

    ! ! # Compute drainage; 1st layer of soil --> 2nd layer of soil
    D = Drainage(Wu, smcap)

    ! ! # Compute baseflow (bottom layer to runoff)
    Qb = Wb * frtgw

    ! Compute PET
    ! First compute saturation specific humidity
    pet_value = PET(L, Td, etpar)

    ! Compute AET
    E = pet_value*(Wu/smcap)

    ! Compute change in the water balance of top layer
    dWudt = Wi - D - E

    ! ! # Compute change in water balance at bottom layer
    dWbdt = D - Qb

    ! ! # update values
    ! ! # soil moisture
    Wb = Wb + dWbdt
    Wu = Wu + dWudt

    ! ! Compute snow
    Snow = Snow + Ps - Md

    ! ! Compute discharge
    Q = Qb + Qd

end subroutine one_step


end module fast_foward_model

