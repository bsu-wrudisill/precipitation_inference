module test
contains

subroutine subr1(a,b,c)
implicit none
real, intent(in) :: a,b
real, intent(out) :: c
c = b+a
end subroutine subr1


end module test
