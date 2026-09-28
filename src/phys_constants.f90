module phys_constants
    ! Shared physical constants for coupled FESM programs.
    !
    ! This module provides the SHAPE of a physical-constants record and the
    ! machinery to fill, derive, check, record and compare it. It deliberately
    ! defines NO numeric values: the values are the program's policy, not the
    ! library's. A program supplies them either from a namelist
    ! (`phys_const_load`) or from its own parameters (`phys_const_set`), then
    ! hands the resulting object to every component that needs a constant.
    !
    ! A reference set of Earth values is shipped as DATA in
    ! `par/phys_const_earth.nml`. Pass it as `defaults_file` to
    ! `phys_const_load` to have it act as the fallback under a program's own
    ! (possibly partial) namelist group; nothing here reads it implicitly.
    !
    ! Components are `real(dp)` regardless of any program's working precision,
    ! so one build of this library serves an `sp` program and a `dp` program at
    ! the same time. Consumers narrow once, where they copy into their own
    ! parameter struct:
    !
    !     mshlf%par%rho_ice = real(cnst%rho_ice, wp)
    !
    ! There are no defaults anywhere: every primitive starts at
    ! PHYS_CONST_UNSET, and `phys_const_require` / `phys_const_get` turn a
    ! forgotten value into a loud failure at init rather than a silently wrong
    ! number at runtime.
    !
    ! Time is NOT part of the record. Which year length a program uses is its
    ! calendar choice, so no field or derived quantity here depends on one. The
    ! named year and day lengths below are a different thing: definitions of the
    ! standard conventions, offered so that a program can name the convention it
    ! picked instead of scattering bare literals. Choosing among them stays the
    ! program's business.

    use precision, only: sp, dp
    use nml
    use ncio, only: nc_write_attr

    implicit none

    private

    ! Sentinel for "this value was never supplied". All physical quantities
    ! carried here are positive, so a negative sentinel cannot be mistaken for
    ! a legitimate value.
    real(dp), parameter :: PHYS_CONST_UNSET = -9999.0_dp

    ! =========================================================================
    ! Calendar conventions
    !
    ! Definitions, not policy: each is exactly `n_days * sec_day` for the stated
    ! day count, so naming one records which convention a program adopted. They
    ! are deliberately not fields of phys_const_class and nothing derived uses
    ! them -- a program that needs a rate conversion multiplies by whichever it
    ! chose, at the point where it converts.
    !
    ! The four values in current use across these models are all legitimate and
    ! distinct conventions rather than roundings of one another, which is exactly
    ! why they are worth naming.
    ! =========================================================================

    real(dp), parameter :: sec_day  = 86400.0_dp        ! [s] nominal day
    real(dp), parameter :: sec_hour =  3600.0_dp        ! [s]
    real(dp), parameter :: sec_min  =    60.0_dp        ! [s]

    ! Year lengths, as an exact number of nominal days.
    real(dp), parameter :: sec_year_360d      = 360.0_dp      * sec_day  ! 31104000
    real(dp), parameter :: sec_year_365d      = 365.0_dp      * sec_day  ! 31536000
    real(dp), parameter :: sec_year_366d      = 366.0_dp      * sec_day  ! 31622400
    real(dp), parameter :: sec_year_julian    = 365.25_dp     * sec_day  ! 31557600
    real(dp), parameter :: sec_year_gregorian = 365.2425_dp   * sec_day  ! 31556952
    real(dp), parameter :: sec_year_tropical  = 365.2422_dp   * sec_day  ! 31556926.08

    ! Day lengths for calendars that divide a full astronomical year into a
    ! fixed number of equal, longer-than-nominal days. A 360-day calendar of
    ! this kind has days of ~87658 s, not 86400 s.
    real(dp), parameter :: sec_day_360d_tropical = sec_year_tropical / 360.0_dp
    real(dp), parameter :: sec_day_365d_tropical = sec_year_tropical / 365.0_dp

    ! =========================================================================
    ! Constant names, for phys_const_get
    !
    ! Fortran has no reflection, so a consumer names the quantity it wants
    ! rather than reaching into the record. The name is the same word used for
    ! the field, for the namelist key and in every message, so there is one
    ! vocabulary to learn and it is the one already visible in
    ! par/phys_const_earth.nml. Matching ignores case and trailing blanks.
    !
    ! This list is what phys_const_get prints when a name is not recognised, and
    ! test/phys_constants resolves every entry in it against the shipped
    ! reference set -- so a name listed here without a matching entry in
    ! phys_const_lookup, or a quantity the reference set forgets to supply,
    ! fails the test rather than a run.
    ! =========================================================================

    integer, parameter :: PHYS_NAME_LEN = 20

    ! The first PHYS_N_PRIMITIVE entries are the quantities a program supplies;
    ! the rest are computed by phys_const_derive. The split lets a failure report
    ! which primitives are missing, which is what a reader needs when a derived
    ! quantity comes out unset.
    integer, parameter :: PHYS_N_PRIMITIVE = 12

    character(len=PHYS_NAME_LEN), parameter :: PHYS_CONST_NAMES(22) = [ &
        ! Primitives
        character(len=PHYS_NAME_LEN) :: &
        "g                   ", "T0                  ", &
        "rho_ice             ", "rho_w               ", &
        "rho_sw              ", "rho_asth            ", &
        "L_ice               ", "cp_ice              ", &
        "cp_w                ", "cp_ocn              ", &
        "T_pmp_beta          ", "area_seasurf        ", &
        ! Derived
        "rho_ice_sw          ", "rho_sw_ice          ", &
        "conv_we_ie          ", "conv_ie_we          ", &
        "conv_mmawe_maie     ", "conv_m3_Gt          ", &
        "conv_km3_Gt         ", "conv_millionkm3_Gt  ", &
        "conv_km3_sle        ", "omega_melt          " ]

    type :: phys_const_class

        logical            :: initialized = .false.  ! set by _load / _set
        character(len=64)  :: label       = ""       ! e.g. "Earth", "EISMINT", "climber-x"
        character(len=512) :: source      = ""       ! file path, or "<program>:<module>"

        ! --- Primitives: supplied by the program -------------------------------
        real(dp) :: g            = PHYS_CONST_UNSET  ! [m s-2]     Gravitational acceleration
        real(dp) :: T0           = PHYS_CONST_UNSET  ! [K]         Reference freezing temperature
        real(dp) :: rho_ice      = PHYS_CONST_UNSET  ! [kg m-3]    Density of ice
        real(dp) :: rho_w        = PHYS_CONST_UNSET  ! [kg m-3]    Density of pure water
        real(dp) :: rho_sw       = PHYS_CONST_UNSET  ! [kg m-3]    Density of seawater
        real(dp) :: rho_asth     = PHYS_CONST_UNSET  ! [kg m-3]    Density of the asthenosphere
        real(dp) :: L_ice        = PHYS_CONST_UNSET  ! [J kg-1]    Latent heat of fusion, ice/water
        real(dp) :: cp_ice       = PHYS_CONST_UNSET  ! [J kg-1 K-1] Specific heat capacity of ice
        real(dp) :: cp_w         = PHYS_CONST_UNSET  ! [J kg-1 K-1] Specific heat capacity of pure water
        real(dp) :: cp_ocn       = PHYS_CONST_UNSET  ! [J kg-1 K-1] Specific heat capacity, ocean mixed layer
        real(dp) :: T_pmp_beta   = PHYS_CONST_UNSET  ! [K Pa-1]    Melting-point pressure slope
        real(dp) :: area_seasurf = PHYS_CONST_UNSET  ! [km2]       Global sea-surface area

        ! --- Derived: computed by phys_const_derive, never supplied ------------
        ! Each is left UNSET when its inputs are.
        real(dp) :: rho_ice_sw   = PHYS_CONST_UNSET  ! [1] rho_ice/rho_sw, floating-ice draft fraction
        real(dp) :: rho_sw_ice   = PHYS_CONST_UNSET  ! [1] rho_sw/rho_ice
        real(dp) :: conv_we_ie   = PHYS_CONST_UNSET  ! [1] water equiv. => ice equiv.
        real(dp) :: conv_ie_we   = PHYS_CONST_UNSET  ! [1] ice equiv.   => water equiv.
        real(dp) :: conv_mmawe_maie = PHYS_CONST_UNSET ! [m a-1 / mm a-1] mm a-1 w.e. => m a-1 i.e.
        real(dp) :: conv_m3_Gt   = PHYS_CONST_UNSET  ! [Gt m-3]  m3 ice => Gt
        real(dp) :: conv_km3_Gt  = PHYS_CONST_UNSET  ! [Gt km-3] km3 ice => Gt
        real(dp) :: conv_millionkm3_Gt = PHYS_CONST_UNSET ! [Gt / 1e6 km3]
        real(dp) :: conv_km3_sle = PHYS_CONST_UNSET  ! [m km-3]  km3 ice above flotation => m s.l.e.
        real(dp) :: omega_melt   = PHYS_CONST_UNSET  ! [K-1]     rho_sw*cp_ocn/(rho_ice*L_ice)

    end type phys_const_class

    ! Fetch one constant into a consumer's own working precision, checking that
    ! it was supplied. Resolved on the kind of the destination, so a consumer
    ! built at sp and one built at dp both call it the same way.
    interface phys_const_get
        module procedure phys_const_get_sp
        module procedure phys_const_get_dp
    end interface phys_const_get

    public :: phys_const_class
    public :: PHYS_CONST_UNSET
    public :: phys_const_isset
    public :: phys_const_load
    public :: phys_const_set
    public :: phys_const_derive
    public :: phys_const_require
    public :: phys_const_get
    public :: phys_const_write
    public :: phys_const_log
    public :: phys_const_to_nc
    public :: phys_const_compare

    ! Calendar conventions
    public :: sec_min, sec_hour, sec_day
    public :: sec_year_360d, sec_year_365d, sec_year_366d
    public :: sec_year_julian, sec_year_gregorian, sec_year_tropical
    public :: sec_day_360d_tropical, sec_day_365d_tropical

    ! The recognised constant names, for phys_const_get
    public :: PHYS_CONST_NAMES, PHYS_NAME_LEN

contains

    elemental function phys_const_isset(value) result(isset)
        ! .true. when `value` holds a supplied number rather than the sentinel.
        !
        ! Compared with a tolerance rather than for exact equality, so that a
        ! value which has been through a text round-trip (phys_const_write then
        ! phys_const_load) is still recognised. The sentinel is far from any
        ! physical constant this record carries, so the tolerance cannot mask a
        ! real value.

        real(dp), intent(IN) :: value
        logical :: isset

        isset = (abs(value - PHYS_CONST_UNSET) > 1.0e-6_dp)

    end function phys_const_isset


    subroutine phys_const_load(c, filename, group, defaults_file, defaults_group, label)
        ! Fill `c` from namelist group `group` of `filename`, then derive.
        !
        ! With `defaults_file` present, each value is taken from that file (in
        ! `defaults_group`, or `group` when not given) and then overridden by
        ! `filename` if the parameter appears there. This is how a program uses
        ! the shipped `par/phys_const_earth.nml` as its reference set while
        ! overriding only what it means to change. Without `defaults_file`,
        ! every parameter must be present in `filename`.

        type(phys_const_class), intent(OUT) :: c
        character(len=*), intent(IN) :: filename
        character(len=*), intent(IN) :: group
        character(len=*), intent(IN), optional :: defaults_file
        character(len=*), intent(IN), optional :: defaults_group
        character(len=*), intent(IN), optional :: label

        character(len=256) :: def_group

        def_group = trim(group)
        if (present(defaults_group)) def_group = trim(defaults_group)

        if (present(defaults_file)) then

            call nml_read(filename, group, "g",            c%g,            &
                          defaults_file=defaults_file, defaults_group=def_group)
            call nml_read(filename, group, "T0",           c%T0,           &
                          defaults_file=defaults_file, defaults_group=def_group)
            call nml_read(filename, group, "rho_ice",      c%rho_ice,      &
                          defaults_file=defaults_file, defaults_group=def_group)
            call nml_read(filename, group, "rho_w",        c%rho_w,        &
                          defaults_file=defaults_file, defaults_group=def_group)
            call nml_read(filename, group, "rho_sw",       c%rho_sw,       &
                          defaults_file=defaults_file, defaults_group=def_group)
            call nml_read(filename, group, "rho_asth",     c%rho_asth,     &
                          defaults_file=defaults_file, defaults_group=def_group)
            call nml_read(filename, group, "L_ice",        c%L_ice,        &
                          defaults_file=defaults_file, defaults_group=def_group)
            call nml_read(filename, group, "cp_ice",       c%cp_ice,       &
                          defaults_file=defaults_file, defaults_group=def_group)
            call nml_read(filename, group, "cp_w",         c%cp_w,         &
                          defaults_file=defaults_file, defaults_group=def_group)
            call nml_read(filename, group, "cp_ocn",       c%cp_ocn,       &
                          defaults_file=defaults_file, defaults_group=def_group)
            call nml_read(filename, group, "T_pmp_beta",   c%T_pmp_beta,   &
                          defaults_file=defaults_file, defaults_group=def_group)
            call nml_read(filename, group, "area_seasurf", c%area_seasurf, &
                          defaults_file=defaults_file, defaults_group=def_group)

        else

            call nml_read(filename, group, "g",            c%g)
            call nml_read(filename, group, "T0",           c%T0)
            call nml_read(filename, group, "rho_ice",      c%rho_ice)
            call nml_read(filename, group, "rho_w",        c%rho_w)
            call nml_read(filename, group, "rho_sw",       c%rho_sw)
            call nml_read(filename, group, "rho_asth",     c%rho_asth)
            call nml_read(filename, group, "L_ice",        c%L_ice)
            call nml_read(filename, group, "cp_ice",       c%cp_ice)
            call nml_read(filename, group, "cp_w",         c%cp_w)
            call nml_read(filename, group, "cp_ocn",       c%cp_ocn)
            call nml_read(filename, group, "T_pmp_beta",   c%T_pmp_beta)
            call nml_read(filename, group, "area_seasurf", c%area_seasurf)

        end if

        c%label  = trim(group)
        if (present(label)) c%label = trim(label)
        c%source = trim(filename)

        c%initialized = .true.

        call phys_const_derive(c)

        return

    end subroutine phys_const_load


    subroutine phys_const_set(c, label, source, g, T0, rho_ice, rho_w, rho_sw, &
                              rho_asth, L_ice, cp_ice, cp_w, cp_ocn,           &
                              T_pmp_beta, area_seasurf)
        ! Fill `c` from values the program already holds (e.g. its own
        ! parameters), then derive. Every quantity is optional: supply what the
        ! program has, and anything omitted stays UNSET so that a consumer
        ! depending on it fails loudly via `phys_const_get` rather than
        ! silently using a wrong number.

        type(phys_const_class), intent(OUT) :: c
        character(len=*), intent(IN), optional :: label
        character(len=*), intent(IN), optional :: source
        real(dp), intent(IN), optional :: g, T0
        real(dp), intent(IN), optional :: rho_ice, rho_w, rho_sw, rho_asth
        real(dp), intent(IN), optional :: L_ice, cp_ice, cp_w, cp_ocn
        real(dp), intent(IN), optional :: T_pmp_beta, area_seasurf

        if (present(label))  c%label  = trim(label)
        if (present(source)) c%source = trim(source)

        if (present(g))            c%g            = g
        if (present(T0))           c%T0           = T0
        if (present(rho_ice))      c%rho_ice      = rho_ice
        if (present(rho_w))        c%rho_w        = rho_w
        if (present(rho_sw))       c%rho_sw       = rho_sw
        if (present(rho_asth))     c%rho_asth     = rho_asth
        if (present(L_ice))        c%L_ice        = L_ice
        if (present(cp_ice))       c%cp_ice       = cp_ice
        if (present(cp_w))         c%cp_w         = cp_w
        if (present(cp_ocn))       c%cp_ocn       = cp_ocn
        if (present(T_pmp_beta))   c%T_pmp_beta   = T_pmp_beta
        if (present(area_seasurf)) c%area_seasurf = area_seasurf

        c%initialized = .true.

        call phys_const_derive(c)

        return

    end subroutine phys_const_set


    subroutine phys_const_derive(c)
        ! Recompute every derived quantity from the primitives. Each is left
        ! UNSET unless all of its inputs are set, so an incomplete record
        ! cannot produce a plausible-looking derived number.
        !
        ! This routine is the single definition of these relations. Any module
        ! that needs one of them takes it from here rather than recomputing it,
        ! which is what keeps the ratios consistent with the densities.

        type(phys_const_class), intent(INOUT) :: c

        c%rho_ice_sw = PHYS_CONST_UNSET
        c%rho_sw_ice = PHYS_CONST_UNSET
        c%conv_we_ie = PHYS_CONST_UNSET
        c%conv_ie_we = PHYS_CONST_UNSET
        c%conv_mmawe_maie = PHYS_CONST_UNSET
        c%conv_m3_Gt  = PHYS_CONST_UNSET
        c%conv_km3_Gt = PHYS_CONST_UNSET
        c%conv_millionkm3_Gt = PHYS_CONST_UNSET
        c%conv_km3_sle = PHYS_CONST_UNSET
        c%omega_melt   = PHYS_CONST_UNSET

        ! Ice/seawater density ratio: the draft fraction of floating ice.
        if (phys_const_isset(c%rho_ice) .and. phys_const_isset(c%rho_sw)) then
            c%rho_ice_sw = c%rho_ice / c%rho_sw
            c%rho_sw_ice = c%rho_sw / c%rho_ice
        end if

        ! Water-equivalent <=> ice-equivalent thickness.
        if (phys_const_isset(c%rho_ice) .and. phys_const_isset(c%rho_w)) then
            c%conv_we_ie = c%rho_w / c%rho_ice
            c%conv_ie_we = c%rho_ice / c%rho_w

            ! [mm a-1 w.e.] => [m a-1 i.e.]
            c%conv_mmawe_maie = 1.0e-3_dp * c%conv_we_ie
        end if

        ! Ice volume => mass.
        if (phys_const_isset(c%rho_ice)) then
            c%conv_m3_Gt         = c%rho_ice * 1.0e-12_dp        ! [kg m-3] * [Gt / 1e12 kg]
            c%conv_km3_Gt        = 1.0e9_dp  * c%conv_m3_Gt      ! [1e9 m3 / km3]
            c%conv_millionkm3_Gt = 1.0e15_dp * c%conv_m3_Gt      ! [1e6 km3] * [1e9 m3 / km3]
        end if

        ! Ice volume above flotation => sea-level equivalent: convert to an
        ! equivalent volume of liquid water and spread it over the ocean.
        if (phys_const_isset(c%rho_ice) .and. phys_const_isset(c%rho_w) .and. &
            phys_const_isset(c%area_seasurf)) then
            c%conv_km3_sle = (c%rho_ice / c%rho_w) * 1.0e9_dp / (c%area_seasurf * 1.0e6_dp)
        end if

        ! Thermal-forcing melt scaling.
        if (phys_const_isset(c%rho_sw) .and. phys_const_isset(c%cp_ocn) .and. &
            phys_const_isset(c%rho_ice) .and. phys_const_isset(c%L_ice)) then
            c%omega_melt = (c%rho_sw * c%cp_ocn) / (c%rho_ice * c%L_ice)
        end if

        return

    end subroutine phys_const_derive


    subroutine phys_const_require(c, caller)
        ! Assert that `c` was actually filled. Call this first in any routine
        ! that receives a phys_const_class, so that an unfilled record stops
        ! the program at init instead of propagating sentinel values.

        type(phys_const_class), intent(IN) :: c
        character(len=*), intent(IN) :: caller

        if (.not. c%initialized) then
            write(*,*) ""
            write(*,*) "phys_const_require:: Error: physical constants not initialized."
            write(*,*) "  Required by: ", trim(caller)
            write(*,*) "  Fill the record with phys_const_load or phys_const_set"
            write(*,*) "  before passing it to any component."
            error stop "Program stopped."
        end if

        return

    end subroutine phys_const_require


    subroutine phys_const_get_sp(c, name, out)
        ! Copy one constant into a consumer's own working precision, refusing to
        ! hand over a value that was never supplied.
        !
        !   call phys_const_get(cnst, "cp_ocn", mshlf%par%cp_ocn)
        !
        ! Fusing the lookup, the check and the narrowing is the point: a consumer
        ! cannot take a constant without asserting that it exists, it never
        ! reaches into the record itself, and the conversion from the dp
        ! interchange record happens here rather than at every assignment.
        !
        ! `name` is the field's own name, which is also its namelist key. Case
        ! and trailing blanks are ignored.

        type(phys_const_class), intent(IN) :: c
        character(len=*), intent(IN) :: name
        real(sp), intent(OUT) :: out

        real(dp) :: value

        call phys_const_lookup(c, name, value)
        out = real(value, sp)

        return

    end subroutine phys_const_get_sp


    subroutine phys_const_get_dp(c, name, out)
        ! Double-precision form of phys_const_get; see phys_const_get_sp.

        type(phys_const_class), intent(IN) :: c
        character(len=*), intent(IN) :: name
        real(dp), intent(OUT) :: out

        call phys_const_lookup(c, name, out)

        return

    end subroutine phys_const_get_dp


    subroutine phys_const_lookup(c, name, value)
        ! Resolve a name to its value, refusing a name that is not recognised or
        ! a quantity that was never supplied.

        type(phys_const_class), intent(IN) :: c
        character(len=*), intent(IN) :: name
        real(dp), intent(OUT) :: value

        integer :: i
        logical :: known
        character(len=256) :: missing

        call phys_const_resolve(c, name, value, known)

        if (.not. known) then
            write(*,*) ""
            write(*,*) "phys_const_get:: Error: no constant named '"//trim(name)//"'."
            write(*,*) "  Either the name is misspelled, or it is a field of"
            write(*,*) "  phys_const_class with no entry in phys_const_resolve -- in which"
            write(*,*) "  case add it there, beside the others, and to PHYS_CONST_NAMES."
            write(*,*) "  Recognised names:"
            do i = 1, size(PHYS_CONST_NAMES)
                write(*,*) "    "//trim(PHYS_CONST_NAMES(i))
            end do
            error stop "Program stopped."
        end if

        if (.not. phys_const_isset(value)) then
            missing = unset_primitives(c)
            write(*,*) ""
            write(*,*) "phys_const_get:: Error: physical constant '"//trim(name)// &
                       "' is not defined (UNSET)."
            write(*,*) "  Constants in use: '"//trim(c%label)//"' from '"//trim(c%source)//"'"
            write(*,*) "  Primitives not supplied: "//trim(missing)
            write(*,*) "  A derived constant is UNSET whenever any primitive it is computed"
            write(*,*) "  from is. Define the missing primitives in that parameter file, or"
            write(*,*) "  in the program code that filled the record, and check the run log."
            error stop "Program stopped."
        end if

        return

    end subroutine phys_const_lookup


    function unset_primitives(c) result(list)
        ! Comma-separated names of the primitives this record never received.
        ! This is the actionable part of a lookup failure: a derived quantity is
        ! unset precisely because one of these is.

        type(phys_const_class), intent(IN) :: c
        character(len=256) :: list

        integer  :: i
        real(dp) :: value
        logical  :: known

        list = ""

        do i = 1, PHYS_N_PRIMITIVE
            call phys_const_resolve(c, PHYS_CONST_NAMES(i), value, known)
            if (known .and. .not. phys_const_isset(value)) then
                if (len_trim(list) == 0) then
                    list = trim(PHYS_CONST_NAMES(i))
                else
                    list = trim(list)//", "//trim(PHYS_CONST_NAMES(i))
                end if
            end if
        end do

        if (len_trim(list) == 0) list = "(none -- every primitive was supplied)"

        return

    end function unset_primitives


    pure subroutine phys_const_resolve(c, name, value, known)
        ! Map a name onto the record's field, without judging the result. This is
        ! the one such map; it sits beside the type on purpose, so that a field
        ! and its entry here are read together. `known` is .false. for a name
        ! that has no entry, which lets a caller probe a field without failing.

        type(phys_const_class), intent(IN) :: c
        character(len=*), intent(IN) :: name
        real(dp), intent(OUT) :: value
        logical,  intent(OUT) :: known

        known = .true.
        value = PHYS_CONST_UNSET

        select case(lower(name))

            ! Primitives
            case("g");                  value = c%g
            case("t0");                 value = c%T0
            case("rho_ice");            value = c%rho_ice
            case("rho_w");              value = c%rho_w
            case("rho_sw");             value = c%rho_sw
            case("rho_asth");           value = c%rho_asth
            case("l_ice");              value = c%L_ice
            case("cp_ice");             value = c%cp_ice
            case("cp_w");               value = c%cp_w
            case("cp_ocn");             value = c%cp_ocn
            case("t_pmp_beta");         value = c%T_pmp_beta
            case("area_seasurf");       value = c%area_seasurf

            ! Derived
            case("rho_ice_sw");         value = c%rho_ice_sw
            case("rho_sw_ice");         value = c%rho_sw_ice
            case("conv_we_ie");         value = c%conv_we_ie
            case("conv_ie_we");         value = c%conv_ie_we
            case("conv_mmawe_maie");    value = c%conv_mmawe_maie
            case("conv_m3_gt");         value = c%conv_m3_Gt
            case("conv_km3_gt");        value = c%conv_km3_Gt
            case("conv_millionkm3_gt"); value = c%conv_millionkm3_Gt
            case("conv_km3_sle");       value = c%conv_km3_sle
            case("omega_melt");         value = c%omega_melt

            case DEFAULT
                known = .false.

        end select

        return

    end subroutine phys_const_resolve


    pure function lower(str) result(out)
        ! Lowercase a name for matching, so that callers need not remember
        ! whether a field is written T0, L_ice or rho_ice.

        character(len=*), intent(IN) :: str
        character(len=len_trim(adjustl(str))) :: out

        integer :: i, code

        out = trim(adjustl(str))

        do i = 1, len(out)
            code = iachar(out(i:i))
            if (code >= iachar("A") .and. code <= iachar("Z")) then
                out(i:i) = achar(code - iachar("A") + iachar("a"))
            end if
        end do

        return

    end function lower


    subroutine phys_const_log(c, unit)
        ! Print the record. Used for the run log, so that the constants a
        ! simulation actually ran with are visible in its output.

        type(phys_const_class), intent(IN) :: c
        integer, intent(IN), optional :: unit

        integer :: io

        io = 6
        if (present(unit)) io = unit

        write(io,*) ""
        write(io,*) "phys_constants:: label  = ", trim(c%label)
        write(io,*) "phys_constants:: source = ", trim(c%source)
        if (.not. c%initialized) then
            write(io,*) "    *** NOT INITIALIZED ***"
            return
        end if
        write(io,*) "    g                  = ", c%g
        write(io,*) "    T0                 = ", c%T0
        write(io,*) "    rho_ice            = ", c%rho_ice
        write(io,*) "    rho_w              = ", c%rho_w
        write(io,*) "    rho_sw             = ", c%rho_sw
        write(io,*) "    rho_asth           = ", c%rho_asth
        write(io,*) "    L_ice              = ", c%L_ice
        write(io,*) "    cp_ice             = ", c%cp_ice
        write(io,*) "    cp_w               = ", c%cp_w
        write(io,*) "    cp_ocn             = ", c%cp_ocn
        write(io,*) "    T_pmp_beta         = ", c%T_pmp_beta
        write(io,*) "    area_seasurf       = ", c%area_seasurf
        write(io,*) "  derived:"
        write(io,*) "    rho_ice_sw         = ", c%rho_ice_sw
        write(io,*) "    rho_sw_ice         = ", c%rho_sw_ice
        write(io,*) "    conv_we_ie         = ", c%conv_we_ie
        write(io,*) "    conv_ie_we         = ", c%conv_ie_we
        write(io,*) "    conv_mmawe_maie    = ", c%conv_mmawe_maie
        write(io,*) "    conv_m3_Gt         = ", c%conv_m3_Gt
        write(io,*) "    conv_km3_Gt        = ", c%conv_km3_Gt
        write(io,*) "    conv_millionkm3_Gt = ", c%conv_millionkm3_Gt
        write(io,*) "    conv_km3_sle       = ", c%conv_km3_sle
        write(io,*) "    omega_melt         = ", c%omega_melt
        write(io,*) ""

        return

    end subroutine phys_const_log


    subroutine phys_const_write(c, filename, group)
        ! Write the record as a namelist group. A program calls this once per
        ! run, into its output directory, so that the constants are archived
        ! with the results rather than only in a shared input tree that may
        ! change afterwards. The file it writes is a valid input to
        ! `phys_const_load`.

        type(phys_const_class), intent(IN) :: c
        character(len=*), intent(IN) :: filename
        character(len=*), intent(IN), optional :: group

        integer :: io
        character(len=64) :: grp

        grp = "phys_const"
        if (present(group)) grp = trim(group)

        open(newunit=io, file=trim(filename), status="replace", action="write")

        write(io,"(a)")    "! Physical constants used by this simulation."
        write(io,"(a,a)")  "! label  : ", trim(c%label)
        write(io,"(a,a)")  "! source : ", trim(c%source)
        write(io,"(a)")    "! Written by phys_const_write; valid input to phys_const_load."
        write(io,"(a)")    ""
        write(io,"(a)")    "&"//trim(grp)
        call write_par(io, "g",            c%g,            "[m s-2]      Gravitational acceleration")
        call write_par(io, "T0",           c%T0,           "[K]          Reference freezing temperature")
        call write_par(io, "rho_ice",      c%rho_ice,      "[kg m-3]     Density of ice")
        call write_par(io, "rho_w",        c%rho_w,        "[kg m-3]     Density of pure water")
        call write_par(io, "rho_sw",       c%rho_sw,       "[kg m-3]     Density of seawater")
        call write_par(io, "rho_asth",     c%rho_asth,     "[kg m-3]     Density of the asthenosphere")
        call write_par(io, "L_ice",        c%L_ice,        "[J kg-1]     Latent heat of fusion")
        call write_par(io, "cp_ice",       c%cp_ice,       "[J kg-1 K-1] Heat capacity of ice")
        call write_par(io, "cp_w",         c%cp_w,         "[J kg-1 K-1] Heat capacity of pure water")
        call write_par(io, "cp_ocn",       c%cp_ocn,       "[J kg-1 K-1] Heat capacity, ocean mixed layer")
        call write_par(io, "T_pmp_beta",   c%T_pmp_beta,   "[K Pa-1]     Melting-point pressure slope")
        call write_par(io, "area_seasurf", c%area_seasurf, "[km2]        Global sea-surface area")
        write(io,"(a)")    "/"

        close(io)

        return

    end subroutine phys_const_write


    subroutine write_par(io, name, value, comment)
        ! One namelist line, or a commented-out line when the value is unset.

        integer, intent(IN) :: io
        character(len=*), intent(IN) :: name
        real(dp), intent(IN) :: value
        character(len=*), intent(IN) :: comment

        if (phys_const_isset(value)) then
            write(io,"(4x,a,t20,a,es16.8,4x,a,a)") trim(name), "= ", value, "! ", trim(comment)
        else
            write(io,"(4x,a,a,t22,a,a,a)") "! ", trim(name), "= <unset>    ! ", trim(comment), &
                                           "  [not supplied by this program]"
        end if

        return

    end subroutine write_par


    subroutine phys_const_to_nc(c, filename)
        ! Attach the record to a NetCDF file as global attributes, so that an
        ! output file carries the constants it was produced with.

        type(phys_const_class), intent(IN) :: c
        character(len=*), intent(IN) :: filename

        call nc_write_attr(filename, "phys_const_label",  trim(c%label))
        call nc_write_attr(filename, "phys_const_source", trim(c%source))

        call put_attr(filename, "g",            c%g)
        call put_attr(filename, "T0",           c%T0)
        call put_attr(filename, "rho_ice",      c%rho_ice)
        call put_attr(filename, "rho_w",        c%rho_w)
        call put_attr(filename, "rho_sw",       c%rho_sw)
        call put_attr(filename, "rho_asth",     c%rho_asth)
        call put_attr(filename, "L_ice",        c%L_ice)
        call put_attr(filename, "cp_ice",       c%cp_ice)
        call put_attr(filename, "cp_w",         c%cp_w)
        call put_attr(filename, "cp_ocn",       c%cp_ocn)
        call put_attr(filename, "T_pmp_beta",   c%T_pmp_beta)
        call put_attr(filename, "area_seasurf", c%area_seasurf)

        return

    end subroutine phys_const_to_nc


    subroutine put_attr(filename, name, value)
        ! One global attribute, skipped when the value is unset.

        character(len=*), intent(IN) :: filename
        character(len=*), intent(IN) :: name
        real(dp), intent(IN) :: value

        if (phys_const_isset(value)) then
            call nc_write_attr(filename, "phys_const_"//trim(name), value)
        end if

        return

    end subroutine put_attr


    subroutine phys_const_compare(c1, c2, label1, label2, tol, strict, n_diff)
        ! Report where two records disagree.
        !
        ! This is for coupled configurations in which a component insists on
        ! its own constants (a GIA model calibrated to a benchmark density, a
        ! vendored ice-sheet model reading its own namelist). Comparing makes
        ! the disagreement an explicit, logged decision instead of a silent
        ! inconsistency in the coupled mass budget. With `strict`, a difference
        ! is an error rather than a warning.

        type(phys_const_class), intent(IN) :: c1, c2
        character(len=*), intent(IN), optional :: label1, label2
        real(dp), intent(IN), optional :: tol
        logical,  intent(IN), optional :: strict
        integer,  intent(OUT), optional :: n_diff

        character(len=64) :: nm1, nm2
        real(dp) :: rtol
        logical  :: is_strict
        integer  :: ndiff

        nm1 = trim(c1%label)
        nm2 = trim(c2%label)
        if (present(label1)) nm1 = trim(label1)
        if (present(label2)) nm2 = trim(label2)

        rtol = 1.0e-6_dp
        if (present(tol)) rtol = tol

        is_strict = .false.
        if (present(strict)) is_strict = strict

        ndiff = 0

        write(*,*) ""
        write(*,*) "phys_const_compare:: ", trim(nm1), " vs ", trim(nm2), &
                   "  (relative tolerance ", rtol, ")"

        call cmp("g",            c1%g,            c2%g)
        call cmp("T0",           c1%T0,           c2%T0)
        call cmp("rho_ice",      c1%rho_ice,      c2%rho_ice)
        call cmp("rho_w",        c1%rho_w,        c2%rho_w)
        call cmp("rho_sw",       c1%rho_sw,       c2%rho_sw)
        call cmp("rho_asth",     c1%rho_asth,     c2%rho_asth)
        call cmp("L_ice",        c1%L_ice,        c2%L_ice)
        call cmp("cp_ice",       c1%cp_ice,       c2%cp_ice)
        call cmp("cp_w",         c1%cp_w,         c2%cp_w)
        call cmp("cp_ocn",       c1%cp_ocn,       c2%cp_ocn)
        call cmp("T_pmp_beta",   c1%T_pmp_beta,   c2%T_pmp_beta)
        call cmp("area_seasurf", c1%area_seasurf, c2%area_seasurf)

        if (ndiff == 0) then
            write(*,*) "    no differences."
        else
            write(*,*) "    ", ndiff, " difference(s)."
            if (is_strict) then
                write(*,*) ""
                write(*,*) "phys_const_compare:: Error: constants differ between ", &
                           trim(nm1), " and ", trim(nm2), "."
                error stop "Program stopped."
            end if
        end if
        write(*,*) ""

        if (present(n_diff)) n_diff = ndiff

        return

    contains

        subroutine cmp(name, v1, v2)
            ! Compare one quantity. Skipped unless both records supplied it:
            ! a value one program does not carry is not a disagreement.

            character(len=*), intent(IN) :: name
            real(dp), intent(IN) :: v1, v2

            real(dp) :: denom, reldiff

            if (.not. (phys_const_isset(v1) .and. phys_const_isset(v2))) return

            denom = max(abs(v1), abs(v2))
            if (denom < tiny(denom)) return

            reldiff = abs(v1 - v2) / denom

            if (reldiff > rtol) then
                ndiff = ndiff + 1
                write(*,"(4x,a,t20,es16.8,3x,es16.8,4x,a,f10.4,a)") &
                    trim(name), v1, v2, "(", 100.0_dp*reldiff, " %)"
            end if

            return

        end subroutine cmp

    end subroutine phys_const_compare

end module phys_constants
