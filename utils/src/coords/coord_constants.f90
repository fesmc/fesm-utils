module coord_constants
    ! Domain constants for the coords module set (math constants, missing-value
    ! and error sentinels). Precision kinds are NOT defined or re-exported here:
    ! source `dp`/`sp` directly from the fesm-utils `precision` module.

    use precision, only: dp, sp

    implicit none

    private :: dp, sp   ! used here for kind specification only; not re-exported

    ! Working precision of the coords library (pinned to dp; see coords-design.md).
    integer,  parameter :: prec = dp

    ! Missing value and aliases
    real(dp), parameter :: MISSING_VALUE_DEFAULT = -9999.0_dp
    real(dp), parameter :: mv = MISSING_VALUE_DEFAULT

    ! Error distance (very large) and error index
    real(dp), parameter :: ERR_DIST = 1E8_dp
    integer,  parameter :: ERR_IND  = -1

    ! Mathematical constants
    real(dp), parameter :: pi = 2._dp*acos(0._dp)
    real(dp), parameter :: degrees_to_radians = pi / 180._dp
    real(dp), parameter :: radians_to_degrees = 180._dp / pi

    public :: prec, MISSING_VALUE_DEFAULT, mv, ERR_DIST, ERR_IND, &
              pi, degrees_to_radians, radians_to_degrees

end module coord_constants
