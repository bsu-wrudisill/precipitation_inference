module evaporation
implicit none
contains

! ET FUNCTIONS
real function ETa(PET, &
	              sm1,  &
	              sm1max)
    ! ET parameterization (a). ET comes only from the top soil moisture layer
    real, intent(in) :: PET
    real, intent(in) :: sm1
    real, intent(in) :: sm1max
    ETa = PET*min(sm1,sm1max)/sm1max
end function ETa


real function ETb(PET, &
	              sm1,  &
	              sm1max)
    ! ET parameterization (a). ET comes only from the top soil moisture layer
    real, intent(in) :: PET
    real, intent(in) :: sm1
    real, intent(in) :: sm1max
    ETb = PET*min(sm1,sm1max)/sm1max
end function ETb



end module evaporation


