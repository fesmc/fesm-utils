module fesm_units
    ! Conversion factors between units: x_to = scale*x_from + offset.
    !
    ! fesm_fields recognises equal units spelled differently
    ! (units_canonical); this module adds the factors that fesm_fields leaves
    ! to adapters, in its own module so that model adapters can use it too.
    ! Small and explicit, no udunits. Classes (every row also inverts, and
    ! any two units of one class convert into each other):
    !
    !   class          units                                   canonical      factor / offset
    !   mass flux      kg m-2 {s,d,yr}-1, mm [w.e.] {s,d,yr}-1  kg m-2 yr-1    x sec_year/{1, 86400, sec_year}
    !                  m w.e. yr-1                                             x 1000
    !                  m i.e. yr-1, m ice yr-1                                 x rho_ice (absent: error)
    !   temperature    degC                                    K              + 273.15
    !   heat flux      mW m-2                                  W m-2          x 1e-3
    !   length         km                                      m              x 1000
    !   velocity       m {s,d}-1                               m yr-1         x sec_year/{1, 86400}
    !   dimensionless  %                                       1              x 0.01
    !   salinity       g kg-1                                  psu            1
    !
    ! Units are matched after units_canonical (lower case, no blanks, `*`
    ! or `^`, aliases resolved), so "mm w.e. a-1", "mm/d", "kg m**-2 s**-1"
    ! and "degrees_Celsius" are recognised. Rates are written `<q> <t>-1`
    ! or `<q>/<t>`, with t in s | sec | second, d | day, yr | a | year.
    ! Units with the same canonical key convert with scale = 1, offset = 0
    ! (identity, also for units this module does not know).
    !
    ! sec_year: the run's year length [s] when the caller passes it, else
    ! 365.25 d (as FesmData products are built).

    use precision,   only : dp
    use fesm_log,    only : fesm_error
    use fesm_fields, only : units_canonical

    implicit none

    private

    public :: units_conv

    real(dp), parameter :: SEC_DAY  = 86400.0_dp
    real(dp), parameter :: SEC_YEAR_DEFAULT = 365.25_dp*SEC_DAY

    ! Unit classes
    integer, parameter :: C_NONE = 0, C_MASS_FLUX = 1, C_TEMP = 2, C_HEAT_FLUX = 3, &
                          C_LENGTH = 4, C_VELOCITY = 5, C_DIMLESS = 6, C_SALINITY = 7

    ! Factors that depend on the arguments
    integer, parameter :: F_ONE = 0, F_RHO_ICE = 1

contains

    subroutine units_conv(from, to, scale, offset, rho_ice, sec_year, stat)
        ! scale, offset such that x_to = scale*x_from + offset. Stops when
        ! no conversion is known (different classes, unknown units, or ice
        ! equivalent without rho_ice); with `stat` present, sets stat = 1
        ! (scale = 1, offset = 0) and returns instead.

        character(len=*), intent(in)  :: from, to
        real(dp),         intent(out) :: scale, offset
        real(dp),         intent(in),  optional :: rho_ice    ! [kg m-3], enables m i.e. yr-1
        real(dp),         intent(in),  optional :: sec_year   ! [s], default 365.25 d
        integer,          intent(out), optional :: stat

        integer  :: cf, ct, ff, ft
        real(dp) :: mf, mt, tf, tt, bf, bt, sy
        character(len=256) :: kf, kt
        character(len=:), allocatable :: why

        if (present(stat)) stat = 0
        scale  = 1.0_dp
        offset = 0.0_dp

        kf = units_canonical(from)
        kt = units_canonical(to)
        if (kf == kt) return

        sy = SEC_YEAR_DEFAULT
        if (present(sec_year)) sy = sec_year

        call classify(kf, cf, mf, tf, bf, ff, sy)
        call classify(kt, ct, mt, tt, bt, ft, sy)

        why = ""
        if (cf == C_NONE) why = "'"//trim(from)//"' is not a known unit"
        if (ct == C_NONE .and. len(why) == 0) why = "'"//trim(to)//"' is not a known unit"
        if (len(why) == 0 .and. cf /= ct) why = "the units are of different quantities"
        if (len(why) == 0 .and. (ff == F_RHO_ICE .or. ft == F_RHO_ICE) .and. .not. present(rho_ice)) &
            why = "ice equivalent units need rho_ice"
        if (len(why) > 0) then
            call fesm_error("units_conv", "no conversion from '"//trim(from)//"' to '"//trim(to)//"'", &
                "from = '"//trim(from)//"' (key '"//trim(kf)//"')"//new_line("a")// &
                "to   = '"//trim(to)//"' (key '"//trim(kt)//"')"//new_line("a")//why, stat)
            return
        end if

        if (ff == F_RHO_ICE) mf = mf*rho_ice
        if (ft == F_RHO_ICE) mt = mt*rho_ice

        ! Rates: x_c = x*m*(sec_year/t); the year length cancels between two
        ! rates. Other classes have t = 1 and an offset. The products are
        ! exact for the factors of the table, so scale is rounded once.
        scale  = (mf*tt)/(mt*tf)
        offset = (bf - bt)/mt

    end subroutine units_conv


    subroutine classify(key, cls, m, t, b, f, sec_year)
        ! Class of a units key and its relation to the canonical units of the
        ! class: x_canonical = x*m*(sec_year/t) + b for rates (t in s), else
        ! x*m + b (t = 1). f = F_RHO_ICE: m is to be multiplied by rho_ice.

        character(len=*), intent(in)  :: key
        integer,          intent(out) :: cls, f
        real(dp),         intent(out) :: m, t, b
        real(dp),         intent(in)  :: sec_year

        character(len=256) :: q
        real(dp) :: ts

        cls = C_NONE
        f   = F_ONE
        m   = 1.0_dp
        t   = 1.0_dp
        b   = 0.0_dp

        select case(trim(key))
            case("k", "degk", "degreesk", "degreek")
                cls = C_TEMP
            case("degc", "degreesc", "degreec", "degrees_c", "degree_c", "deg_c", "c", "celsius", &
                 "degreescelsius", "degreecelsius", "degrees_celsius", "degree_celsius", &
                 "degreescelcius", "degreecelcius", "degrees_celcius", "degree_celcius")
                cls = C_TEMP
                b   = 273.15_dp
            case("wm-2")
                cls = C_HEAT_FLUX
            case("mwm-2", "mw/m2")
                cls = C_HEAT_FLUX
                m   = 1.0e-3_dp
            case("m")
                cls = C_LENGTH
            case("km")
                cls = C_LENGTH
                m   = 1000.0_dp
            case("1")
                cls = C_DIMLESS
            case("%", "percent")
                cls = C_DIMLESS
                m   = 0.01_dp
            case("psu", "gkg-1", "g/kg")
                cls = C_SALINITY
        end select
        if (cls /= C_NONE) return

        ! Rates: "<q><t>-1" or "<q>/<t>"
        call split_rate(key, q, ts, sec_year)
        if (ts <= 0.0_dp) return
        t = ts

        select case(trim(q))
            case("kgm-2", "kg/m2", "mm", "mmw.e.", "mmwe", "mmw.e")
                cls = C_MASS_FLUX
            case("mw.e.", "mwe", "mw.e")
                cls = C_MASS_FLUX
                m   = 1000.0_dp
            case("mi.e.", "mie", "mi.e", "mice")
                cls = C_MASS_FLUX
                f   = F_RHO_ICE
            case("m")
                cls = C_VELOCITY
            case default
                t = 1.0_dp
        end select

    end subroutine classify


    subroutine split_rate(key, q, t, sec_year)
        ! Split a rate key into its quantity q and the length t [s] of its
        ! time unit; t = 0 if key is not a rate.

        character(len=*), intent(in)  :: key
        character(len=*), intent(out) :: q
        real(dp),         intent(out) :: t
        real(dp),         intent(in)  :: sec_year

        character(len=8), parameter :: TOK(8) = [character(len=8) :: &
            "second", "year", "sec", "day", "yr", "s", "d", "a"]
        character(len=256) :: k
        integer :: n, i, p, lt

        q = ""
        t = 0.0_dp
        k = key
        n = len_trim(k)

        p = index(k(1:n), "/", back=.true.)
        if (p > 1 .and. p < n) then
            ! "<q>/<t>"
            q = k(1:p-1)
            t = tunit(k(p+1:n))
            if (t <= 0.0_dp) q = ""
            return
        end if

        if (n < 3) return
        if (k(n-1:n) /= "-1") return
        n = n - 2
        do i = 1, size(TOK)
            lt = len_trim(TOK(i))
            if (n > lt) then
                if (k(n-lt+1:n) == TOK(i)(1:lt)) then
                    q = k(1:n-lt)
                    t = tunit(TOK(i))
                    return
                end if
            end if
        end do

    contains

        real(dp) function tunit(s)
            character(len=*), intent(in) :: s
            select case(trim(s))
                case("s", "sec", "second")
                    tunit = 1.0_dp
                case("d", "day")
                    tunit = SEC_DAY
                case("yr", "a", "year")
                    tunit = sec_year
                case default
                    tunit = 0.0_dp
            end select
        end function tunit

    end subroutine split_rate

end module fesm_units
