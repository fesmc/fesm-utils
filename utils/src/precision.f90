module precision

  implicit none

  ! Floating point section
  
  integer, parameter :: dp  = kind(1.d0)
  integer, parameter :: sp  = kind(1.0)

  ! Set working precision
  integer, parameter :: wp = sp

  ! Define default missing value 
  real(wp), parameter :: mv = -9999.0_wp 

  ! Machine tolerance
  real(wp), parameter :: TOL = 1e-8_wp
  real(wp), parameter :: TOL_UNDERFLOW = real(1e-15,wp)
  
end module precision
