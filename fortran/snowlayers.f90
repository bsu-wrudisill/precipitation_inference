module snowlayers
implicit none
contains

subroutine layerer(nlayers,       &
				   ntimes,        &
				   precipitation, &
				   temperature,   &
				   dz,            &
				   opg,           &
				   bias,          &
                   elevation,     &
                   dt,            &
                   rvs,           &
                   uadj,          &
                   mbase,         &
                   mfmax,         &
                   mfmin,         &
                   tipm,          &
                   nmf,           &
                   plwhc,         &
                   pxtemp,        &
                   pxtemp1,       &
                   pxtemp2,       &
                   snow)


    ! these are the two different options...
    use snowmodule17
    use snowmodule01

    ! Snow "layer" parameters
    integer, intent(in) :: nlayers
    integer, intent(in) :: ntimes

    ! FORCGING variables
    real, intent(in), dimension(0:ntimes) :: precipitation
    real, intent(in), dimension(0:ntimes) :: temperature
    real, intent(in), dimension(0:nlayers) :: dz
    real, intent(in) :: opg
    real, intent(in) :: bias

    ! Snow17 paramters
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

    ! Outputs
    real, dimension(ntimes:nlayers) :: snow


    ! Internal
    integer :: l


    do l=0,nlayers
    call snow17driver(ntimes, &
                      jday, &
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
                      meltVec, &
                      sweVec)


end subroutine layerer



end module snowlayers