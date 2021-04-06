module surfacerunoff
implicit none
contains


! staturated area formulation (a) -- prms (?)
real function satarea_a(s1,    &
                        s1max, &
                        beta)  result(Ac)

    real, intent(in) :: s1
    real, intent(in) :: s1max
    real, intent(in) :: beta

    Ac = 1 - (1 - s1/s1max) **beta
end function satarea_a


! saturated area formulation (b) -- topmodel?
real function satarea_b(s1,    &
                        s1max, &
                        beta)  result(Ac)

    real, intent(in) :: s1
    real, intent(in) :: s1max
    real, intent(in) :: beta

    Ac = 1 - (1 - s1/s1max) **beta
end function satarea_b

! END
end module surfacerunoff




! !!!!!!!!!!!!!!!!

! program test
! use surfacerunoff
! implicit none
! real :: s1
! real :: s1max
! real :: beta

! s1=10.
! s1max=200.
! beta=.5

! print*, satarea_a(s1, s1max, beta)


! end program test