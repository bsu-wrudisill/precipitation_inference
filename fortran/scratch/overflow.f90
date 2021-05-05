module overflow
implicit none
contains

! ET FUNCTIONS
real function Qufof_a(sm1,  &
	                sm1max, &
	                p, &
	                qsx, &
	                w)
    ! ET parameterization (a). ET comes only from the top soil moisture layer
    real, intent(in) :: sm1
    real, intent(in) :: sm1max
    real, intent(in) :: p
    real, intent(in) :: qsx
    real, intent(in) :: w
    real, parameter :: e=5

    !internal
    real :: scale

    scale  = 1./(1 + EXP((sm1 - sm1max - w*e)/e))
    Qufof_a = (p - qsx)*scale

end function Qufof_a

end module overflow
