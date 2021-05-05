module baseflow
implicit none
contains


! Baseflow option a --- (super simple)
real function Qba(sm2, &
	              ks)
    real, intent(in) :: sm2
    real, intent(in) :: ks


    Qba = ks*sm2
end function Qba

! Baseflow option b --- (name)
real function Qbb(sm2,    &
	              sm2max, &
	              ks,    &
	              n)
    real, intent(in) :: sm2
    real, intent(in) :: sm2max
    real, intent(in) :: ks
    integer, intent(in) :: n
    ! done

    Qbb = ks*(sm2/sm2max)**n

end function Qbb

! Baseflow option c -- Topmodel
real function Qbc(sm2,     &
				  sm2max,  &
	              ks,     &
	              lambda, &
	              n)
    real, intent(in) :: sm2
    real, intent(in) :: sm2max
    real, intent(in) :: lambda
    real, intent(in) :: ks
    integer, intent(in) :: n

    ! internal
    real :: m

    ! done
    m = sm2max/n
    Qbc = (ks*m)/(lambda*n)*(sm2/sm2max)**n

end function Qbc

end module baseflow





! program test
! use baseflow
! implicit none
! real :: sm2
! real :: sm2max
! real :: lambda
! real :: ks
! integer :: n

! sm2=100.
! sm2max=200.
! lambda = .3
! ks=.03
! n = 2


! print*, Qba(sm2, ks)
! print*, Qbb(sm2, sm2max, ks, n)
! print*, Qbc(sm2, sm2max, ks, lambda, n)
! end program test
