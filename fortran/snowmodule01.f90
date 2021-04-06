module snowmodule01
implicit none
contains


real function DDPar(T30,&
                    t_base,&
                    t_power)

    real, intent(in) :: T30
    real, intent(in) :: t_base
    real, intent(in) :: t_power

    ! T30: previous 30day mean temp
    ! t_base: parameters
    ! t_power: parameters

    if (T30.gt.t_base) then
        DDPar = (T30 - t_base)**t_power
    else
        DDPar = t_base
    end if
end function DDPar



real function Ms(tair,&
                 swe,&
                 t_melt,&
                 ddpar)

    real, intent(in) :: tair   ! Air temp
    real, intent(in) :: swe    ! Snow Wat Equiv
    real, intent(in) :: t_melt ! (PARAM) melt temperature
    real, intent(in) :: ddpar  ! value

    !internal
    real :: melt

    ! # Compute the snowmelt
    ! # tair: temperature
    ! # swe: amount in the bucket
    ! # ddpar: parameter (but not constant)
    ! # t_melt: parameter

    if (tair.gt.t_melt) then
        melt = (tair-t_melt)*ddpar

        if (melt.gt.swe) then
            Ms = swe

        else
            Ms = melt
        end if
    else
        Ms = 0
    end if
end function Ms



subroutine snow01(precip, &
                  tair, &
                  tair30, &
                  pxtemp, &
                  t_melt, &
                  t_base, &
                  t_power, &
                  swe, &
                  melt_outflow)

    ! forcings
    real, intent(in)     :: precip
    real, intent(in)     :: tair
    real, intent(in)     :: tair30

    ! parameters
    real, intent(in)     :: pxtemp
    real, intent(in)     :: t_melt
    real, intent(in)     :: t_base
    real, intent(in)     :: t_power

    ! states/fluxes
    real, intent(inout)  :: swe
    real, intent(out)    :: melt_outflow
    !f2py intent(in,out) :: swe
    !f2py intent(out) :: melt_outflow

    ! internal
    real :: pr
    real :: ps
    real :: ddpar_value
    real :: Mdk
    real :: Mdsum
    real :: Prsum
    real :: tairavg
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!                 ~~~BEGIN~~~                    !!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! all or nothing snow partitioning
    if (tair > pxtemp) then
        pr = precip
    else
        pr = 0
    end if

    ! create snow and rain values
    ps = precip - pr

    ! Compute melt factors for each layer
    ddpar_value = DDPar(tair30, t_base, t_power)

    ! Compute melt
    melt_outflow = Ms(tair, swe, t_melt, ddpar_value)

    ! melt/accumulate the snowpack
    swe = swe + ps - melt_outflow ! last snow, new snow, snow melt

end subroutine snow01

end module snowmodule01

