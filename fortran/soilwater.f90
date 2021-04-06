module soilwater
implicit none
contains



! Soil State Equations
subroutine s1_topvic(prain,  &
				     pmelt,  &
				     E,      &
				     qsx,    &
				     qif,    &
				     q12,    &
				     qufof,  &
				     sm1)



    implicit none
    ! inputs-- forcings
    real, intent(in) :: prain    ! rain
    real, intent(in) :: pmelt    ! snowmelt
    real, intent(in) :: E        ! evaporation
    real, intent(in) :: qsx      ! overland flow
    real, intent(in) :: qif      ! interflow
    real, intent(in) :: q12      ! interflow (layer 1 --> layer 2)
    real, intent(in) :: qufof    ! overflow from free storage in top layer


    ! output
    real, intent(inout) :: sm1

    ! internal
    real :: dsm1dt   ! Change in the top layer soil water storage

    !!!!!!!!!!!!!!!!!!!!!!!!!!BEGIN!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    ! Equation X in Clark. change in top soil moisture layer w.r.t time.
    dsm1dt = (prain + pmelt - qsx) - E - q12 - qif - qufof

    ! Update soild moisture and return
    sm1 = sm1 + dsm1dt
end subroutine s1_topvic




end module soilwater


!!!!!!! UNIT TESTING !!!!!!!
! program test
! use soilwater
! implicit none
! real :: sm
! integer :: N
! integer :: i
! real :: q12

! sm=110.

! N=300
! do i=1,N
!     call s1_topvic(2.0,  1.0,  .12, .3,  180., 1., 1.,sm)
!     print*, sm
! end do

! q12 = q12a(100.,110.,.05,1)
! print*, q12
! end program test

