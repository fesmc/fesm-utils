module precision

  implicit none

  ! Floating point section
  
  integer, parameter :: dp  = kind(1.d0)
  integer, parameter :: sp  = kind(1.0)

  ! Set working precision
  integer, parameter :: wp = dp

end module precision
