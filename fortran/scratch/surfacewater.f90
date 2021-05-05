module surfacewater
implicit none
contains


! staturated area formulation (a) -- prms (?)
real function satarea_a(sm1,    &
                        sm1max, &
                        beta)  result(Ac)

    real, intent(in) :: sm1
    real, intent(in) :: sm1max
    real, intent(in) :: beta

    Ac = 1 - (1 - sm1/sm1max) **beta
end function satarea_a


! saturated area formulation (b) -- topmodel?
real function satarea_b(sm1,    &
                        sm1max, &
                        beta)  result(Ac)

    real, intent(in) :: sm1
    real, intent(in) :: sm1max
    real, intent(in) :: beta

    Ac = 1 - (1 - sm1/sm1max) **beta
end function satarea_b

! END
end module surfacewater




! !!!!!!!!!!!!!!!!

! program test
! use surfacerunoff
! implicit none
! real :: sm1
! real :: sm1max
! real :: beta

! sm1=10.
! sm1max=200.
! beta=.5

! print*, satarea_a(sm1, sm1max, beta)


! end program test