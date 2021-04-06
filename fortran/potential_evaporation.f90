! Different ways of estimating PET...
module potential_evaportation
implicit none


! PET (a)-- Hamon Equation
real function PETa(L,    &
                   t,    &
                   etpar)

    ! inputs
    real, intent(in) :: L      ! length of day
    real, intent(in) :: t      ! temperature
    real, intent(in) :: etpar  ! scaling parameter

    ! internal
    real :: sat_vap_pres   ! saturation vapor pressure
    real :: abs_humidity   ! staturated absolute humidity

    sat_vap_pres = (.6112) * exp((17.67 * t)/(t + 243.5))
    abs_humidity = (2165 * sat_vap_pres)/(t + 273.15)
    PET=L*abs_humidity*etpar
end function PETa



! PET (b) Makkink equation
real function PETb(K,     &
                   s,     &
                   gamma, &
                   etpar)

    ! inputs
    real, intent(in) :: K     ! Incoming solar radiation (w/m2)
    real, intent(in) :: s     ! slope of saturation water vapor temperature curve at temperature (t)
    real, intent(in) :: gamma ! psychometric constant
    real, intent(in) :: etpar

    PETb = .65 * etpar * (s/(s+gamma)) * K

end function PETb



end module potential_evaportation