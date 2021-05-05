module forcing
   type                           ! Define new type 'pipe', which
     real diameter                    ! is made up of two reals, an
     real flowrate                    ! integer, and a character.
     integer length
     character(len=10) :: flowtype
   end type pipe
end module forcing
