module soilwater
implicit none
contains
! Soil State Equations

! TOP LAYER
subroutine s1_topvic(qin,   &
				     E,     &
				     qsx,   &
				     qif,   &
				     q12,   &
				     qufof, &
				     sm1)



    implicit none
    ! inputs-- forcings
    real, intent(in) :: qin      ! water input rain + melt
    real, intent(in) :: E        ! evaporation
    real, intent(in) :: qsx      ! overland flow
    real, intent(in) :: qif      ! interflow
    real, intent(in) :: q12      ! interflow (layer 1 --> layer 2)
    real, intent(in) :: qufof    ! overflow from free storage in top layer

    ! output
    real, intent(inout) :: sm1

    ! internal
    real :: dsm1dt   ! Change in the top layer soil water storage

    !!!!!!! BEGIN  !!!!!!!!!

    ! Equation X in Clark. change in top soil moisture layer w.r.t time.
    dsm1dt = (qin - qsx) - E - q12 - qif - qufof

    ! Update soild moisture and return
    sm1 = sm1 + dsm1dt
end subroutine s1_topvic




! BOTTOM LAYER
subroutine s2_topprms(q12,  &
                      qb,   &
                      sm2)

    implicit none
    real, intent(in)    :: q12      ! incoming water from layer 1
    real, intent(in)    :: qb       ! baseflow
    real, intent(inout) :: sm2     ! soil water storage in the 2nd layer

    ! internal
    real :: dsm2dt

    !!!!!! BEGIN !!!!!!
    sm2 = sm2 + dsm2dt
end subroutine s2_topprms



end module soilwater
