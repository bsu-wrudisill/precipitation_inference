
! module class_Circle
!   implicit none
!   private
!   real :: pi = 3.1415926535897931d0 ! Class-wide private constant

!   type, public :: Circle
!      real :: radius
!    contains
!      procedure :: area => circle_area
!      procedure :: print => circle_print
!   end type Circle
! contains
!   function circle_area(this) result(area)
!     class(Circle), intent(in) :: this
!     real :: area
!     area = pi * this%radius**2
!   end function circle_area

!   subroutine circle_print(this)
!     class(Circle), intent(in) :: this
!     real :: area
!     area = this%area()  ! Call the type-bound function
!     print *, 'Circle: r = ', this%radius, ' area = ', area
!   end subroutine circle_print
! end module class_Circle


! program circle_test
!   use class_Circle
!   implicit none

!   type(Circle) :: c     ! Declare a variable of type Circle.
!   c = Circle(1.5)       ! Use the implicit constructor, radius = 1.5.
!   call c%print          ! Call the type-bound subroutine
! end program circle_test


!! CLASS !!
module class_snow_model
  implicit none
  private
  type, public :: snow_states
    ! These are all of the possible snow parameters that the multiple
    ! snow models might use.
    real ::


  contains
    procedure :: snow17 => snow17fx
!    procedure :: snow01 => snow01fx
  end type snow_states

! these are the functions that the method contains
contains

subroutine snow17fx(this, newval)
     class(snow_states), intent(inout) :: this
     real, intent(in) :: newval

end subroutine snow17fx
end module class_snow_model




!! CLASS !!
module class_hydro_model
  implicit none
  private
  type, public :: hydro_parameters
    real :: frtdir, frtgw, smcap, etpar, t_snow, t_melt, t_base, t_power, bias, opg
  end type hydro_parameters


  type, public :: hydro_states
    real :: wu, wb, PET
  end type hydro_parameters


  type, public :: hydro_fluxes
    real :: Q, ET
  end type hydro_parameters


end module hydro_model







!! RUN THE TEST
program circle_test
  use class_snow_model
  implicit none

  type(snow_states) :: hp     ! Declare a variable of type Circle.
  real :: foo
  hp = snow_states

  print*, hp%frtdir
  call hp%snow17(20.)
  print*, hp%frtdir

  !hp%frtdir = 12.0 ! this is allowed. we can reassign

end program circle_test



