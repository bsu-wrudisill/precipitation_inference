module percolation
implicit none
contains

! Option A -- VIC
real function q12a(sm1,   &
	               sm1max,&
	               ku,   &
	               c)
    ! percolation function #1
    real, intent(in) :: sm1
    real, intent(in) :: sm1max
    real, intent(in) :: ku
    integer, intent(in) :: c
    q12a = ku * (sm1/sm1max)**c
end function q12a


! Option C -- Sacramento Model
real function q12c(sm1F,   &
	               sm1Fmax,&
	               q0,    &
	 			   sm2,    &
	 			   sm2max, &
	               psi,   &
	               alpha)

    ! percolation function #1
    real, intent(in) :: sm1F     ! "free" soil moisture storage in first layer
    real, intent(in) :: sm1Fmax  ! "free" soil moisture storage in first layer MAXIMUM
    real, intent(in) :: q0       ! basefolw at saturation
    real, intent(in) :: sm2      ! soil moisture storate in 2nd layer
    real, intent(in) :: sm2max   ! soil moisture storage in 2nd layer MAXIMUM
    real, intent(in) :: psi      ! parameter
    real, intent(in) :: alpha    ! parameter

    ! internal
    real :: dlz


    dlz = 1. + alpha * (sm2/sm2max)**psi
    q12c = q0 * dlz * (sm1F/sm1Fmax)
end function q12c

end module percolation