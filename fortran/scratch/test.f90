program test
integer, parameter :: n=2
integer, parameter :: m=3
real, dimension(m,n) :: foo
real, dimension(m,n) :: output
integer :: i,j



do i=1,m
    do j=1,n
        foo(i,j) = 1
! bar = (/3., 6., 8./, /2., 3., 8./)
	end do
end do
print*, foo




end program test