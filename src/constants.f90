module constants
    ! Package-wide domain constants for fesm-utils: missing value, error
    ! sentinels, mathematical constants, and tolerances. Precision kinds are
    ! NOT defined here -- source `wp`/`sp`/`dp` directly from the `precision`
    ! module. Values default to working precision (`wp`); a double-precision
    ! `mv_dp` is provided for the double-only internals (e.g. coords geometry
    ! and weights).

    use precision, only: wp, dp

    implicit none

    private :: wp, dp   ! kinds used here for declarations only; not re-exported

    ! Missing value (working precision), plus a double-precision form.
    real(wp), parameter :: mv    = -9999.0_wp
    real(dp), parameter :: mv_dp = -9999.0_dp

    ! Error distance (very large) and error index
    real(dp), parameter :: ERR_DIST = 1E8_dp
    integer,  parameter :: ERR_IND  = -1

    ! Mathematical constants
    real(dp), parameter :: pi = 2._dp*acos(0._dp)
    real(dp), parameter :: degrees_to_radians = pi / 180._dp
    real(dp), parameter :: radians_to_degrees = 180._dp / pi

    ! Machine tolerances
    real(wp), parameter :: TOL = 1e-8_wp
    real(wp), parameter :: TOL_UNDERFLOW = real(1e-15, wp)

    public :: mv, mv_dp, ERR_DIST, ERR_IND, &
              pi, degrees_to_radians, radians_to_degrees, TOL, TOL_UNDERFLOW

end module constants
