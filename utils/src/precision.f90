module precision
  ! Floating-point kinds for fesm-utils. Domain values (missing value,
  ! tolerances, math constants) live in the `constants` module.

  implicit none

  integer, parameter :: dp = kind(1.d0)
  integer, parameter :: sp = kind(1.0)

  ! Working precision
  integer, parameter :: wp = sp

end module precision
