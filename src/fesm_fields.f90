module fesm_fields
    ! Containers for fields that cross a component boundary.
    !
    ! A field_state_t is one side of an exchange: the exports or the imports
    ! of a component instance. It holds a list of field_t, each a `dp` buffer
    ! v(nx,ny,nz) with the metadata an exchange layer needs: units, grid tag,
    ! validity interval [t0,t1], time semantics (INSTANT, MEAN, ACCUM),
    ! conservation flag (INTENSIVE, EXTENSIVE) and an optional valid mask.
    ! Rank is metadata: scalars are v(1,1,1) with grid = "", 2-D fields have
    ! nz = 1, rank-3 fields name their third axis in `zdim`.
    !
    ! This module knows no field names, does no remapping and no unit
    ! conversion. The vocabulary, the wiring and the lags live in the model
    ! (kryos); the exchange there combines this module with `coupler`. Units
    ! are only recognised as equal when spelled differently (units_canonical).
    !
    ! Life cycle of a state, all through free procedures:
    !
    !     call fields_init(ex, "ice_grl:exports")
    !     call fields_add(ex, "smb", "kg m-2 yr-1", FIELD_MEAN, FIELD_EXTENSIVE, h_smb, &
    !                     grid="grl-16km", dims=[nx,ny])
    !     call fields_lock(ex)                         ! allocates every buffer once
    !     call fields_put(ex, h_smb, smb, t0=ta, t1=tb) ! accumulates (MEAN)
    !     call fields_finalize_interval(ex, tb)        ! publishes the window mean
    !     call fields_get(ex, h_smb, smb_out)          ! the published value
    !
    ! Rules that make the exchange exact and the failures loud:
    !   * Components are public (Fortran has no protected components): read
    !     st%f(h) freely, write only through fields_put.
    !   * Buffers are allocated in fields_lock only; put, get,
    !     finalize_interval and restart_read never allocate them.
    !   * Accumulators reset in fields_finalize_interval, never in fields_get
    !     (several consumers read one published value). get takes the state
    !     intent(in).
    !   * Accumulation windows must be contiguous to the bit (times come from
    !     the clock as n*tick); a gap or an overlap is an error.
    !   * Every misuse stops through fesm_error with the owner and the field
    !     in the message; there are no stat arguments.

    use precision, only : sp, dp
    use fesm_log,  only : LOG_INFO, log_msg, fesm_error, to_str
    use ncio,      only : nc_create, nc_open, nc_close, nc_write_dim, nc_write, nc_write_attr, &
                          nc_read, nc_read_attr, nc_exists_var, nc_exists_attr, nc_dims

    implicit none

    private

    ! Time semantics of a field's value
    integer, parameter, public :: FIELD_INSTANT = 1     ! value at t1 (t0 = t1)
    integer, parameter, public :: FIELD_MEAN    = 2     ! time mean over [t0, t1]
    integer, parameter, public :: FIELD_ACCUM   = 3     ! time integral over [t0, t1]

    ! Conservation flag (selects the default remap method in the exchange)
    integer, parameter, public :: FIELD_INTENSIVE = 1   ! per-point quantity (T_srf, z_srf); bilin
    integer, parameter, public :: FIELD_EXTENSIVE = 2   ! integral to conserve (H_ice, smb); con

    integer, parameter, public :: FIELD_NAME_LEN  = 32
    integer, parameter, public :: FIELD_UNITS_LEN = 32
    integer, parameter, public :: FIELD_GRID_LEN  = 64
    integer, parameter, public :: FIELD_OWNER_LEN = 64

    ! Initial capacity of a state's field list; it doubles while declaring.
    integer, parameter :: CAPACITY0 = 8

    type, public :: field_t
        ! --- Metadata: set by fields_add, fixed after fields_lock ------------
        character(len=FIELD_NAME_LEN)  :: name  = ""   ! e.g. "smb"; unique within a state
        character(len=FIELD_UNITS_LEN) :: units = ""   ! as declared, "kg m-2 yr-1"
        character(len=FIELD_UNITS_LEN) :: ukey  = ""   ! units_canonical(units)
        character(len=FIELD_GRID_LEN)  :: grid  = ""   ! coupler grid tag; "" for rank 0
        character(len=FIELD_NAME_LEN)  :: zdim  = ""   ! name of axis 3, rank 3 only ("month", "depth")
        integer :: rank = 0                            ! 0 | 2 | 3
        integer :: nx = 1, ny = 1, nz = 1              ! rank 0: 1,1,1; rank 2: nz = 1
        integer :: tsem = FIELD_INSTANT
        integer :: cons = FIELD_INTENSIVE
        logical :: masked = .false.                    ! carries a valid mask

        ! --- Published value: what fields_get returns ------------------------
        logical  :: published = .false.
        real(dp) :: t0 = 0.0_dp, t1 = 0.0_dp           ! validity interval [yr]
        real(dp), allocatable :: v(:,:,:)              ! (nx,ny,nz), owned
        logical,  allocatable :: valid(:,:)            ! (nx,ny) iff masked

        ! --- Open accumulation window: MEAN and ACCUM only --------------------
        real(dp), allocatable :: acc(:,:,:)            ! sum of var*(t1-t0) over puts
        logical,  allocatable :: acc_valid(:,:)        ! AND of the put masks (iff masked)
        real(dp) :: wsum = 0.0_dp                      ! sum of weights (t1-t0) [yr]
        real(dp) :: ta0 = 0.0_dp, ta1 = 0.0_dp         ! window; empty while wsum == 0
    end type field_t

    type, public :: field_state_t
        character(len=FIELD_OWNER_LEN) :: owner = ""   ! "ice_grl:exports"; in every message
        integer :: n = 0                               ! number of declared fields
        logical :: locked = .false.
        type(field_t), allocatable :: f(:)             ! f(h), h = 1..n in declaration order;
                                                       ! size(f) == n after fields_lock
    end type field_state_t

    ! Copy in. Specifics {0d,2d,3d} x {dp,sp}; t1 is always dp.
    interface fields_put
        module procedure fields_put_0d_dp, fields_put_0d_sp
        module procedure fields_put_2d_dp, fields_put_2d_sp
        module procedure fields_put_3d_dp, fields_put_3d_sp
    end interface fields_put

    ! Copy out the published value. Specifics as fields_put.
    interface fields_get
        module procedure fields_get_0d_dp, fields_get_0d_sp
        module procedure fields_get_2d_dp, fields_get_2d_sp
        module procedure fields_get_3d_dp, fields_get_3d_sp
    end interface fields_get

    public :: fields_init, fields_add, fields_lock, fields_handle, fields_has
    public :: fields_put, fields_get, fields_finalize_interval
    public :: fields_describe, fields_restart_write, fields_restart_read
    public :: fields_check_link
    public :: units_key, units_canonical, units_equivalent

    ! Unit identities (factor 1, offset 0), as written; matched after
    ! units_key. Units that need a factor or an offset (m i.e., degC, s-1
    ! rates) are deliberately absent: conversions belong in the adapters.
    integer, parameter :: UNITS_LEN = 16
    integer, parameter :: N_ALIAS = 21
    character(len=UNITS_LEN), parameter :: UNITS_ALIAS(2, N_ALIAS) = reshape([ &
        character(len=UNITS_LEN) :: &
        ! canonical       accepted spelling
        "kg m-2 yr-1",    "kg m-2 a-1",       &
        "kg m-2 yr-1",    "kg/m2/yr",         &
        "kg m-2 yr-1",    "mm w.e. yr-1",     &
        "kg m-2 yr-1",    "mm we yr-1",       &
        "kg m-2 yr-1",    "mm w.e./yr",       &
        "kg m-2 yr-1",    "mm/yr w.e.",       &
        "kg m-2",         "kg/m2",            &
        "kg m-2",         "mm w.e.",          &
        "kg m-2",         "mm we",            &
        "m",              "meter",            &
        "m",              "meters",           &
        "m",              "metre",            &
        "m",              "metres",           &
        "m yr-1",         "m a-1",            &
        "m yr-1",         "m/yr",             &
        "m yr-1",         "m/a",              &
        "K",              "kelvin",           &
        "W m-2",          "W/m2",             &
        "W m-2",          "W m**-2",          &
        "1",              "-",                &
        "1",              "fraction"          ], [2, N_ALIAS])

contains

    ! =========================================================================
    ! Declaration
    ! =========================================================================

    subroutine fields_init(st, owner)
        ! Reset `st` to an empty, unlocked state. `owner` names the instance
        ! and role ("ice_grl:exports") and appears in every message.

        type(field_state_t), intent(out) :: st
        character(len=*),    intent(in)  :: owner

        if (len_trim(owner) > FIELD_OWNER_LEN) then
            call fesm_error("fesm_fields_init", "invalid declaration", &
                "owner = '"//trim(owner)//"'"//new_line("a")// &
                "owner longer than FIELD_OWNER_LEN = "//to_str(FIELD_OWNER_LEN))
        end if

        st%owner  = owner
        st%n      = 0
        st%locked = .false.
        allocate(st%f(CAPACITY0))

    end subroutine fields_init


    subroutine fields_add(st, name, units, tsem, cons, h, grid, dims, zdim, masked)
        ! Declare one field and return its handle (= st%n after the call).
        ! dims absent: scalar; dims = [nx,ny]: rank 2; [nx,ny,nz]: rank 3.
        ! `grid` is required iff dims is present; `zdim` (default "z") only
        ! for rank 3; `masked` (default .false.) not for rank 0. The handle is
        ! an argument, not a function result, so that a call with side
        ! effects never sits in an expression.

        type(field_state_t), intent(inout) :: st
        character(len=*),    intent(in)    :: name, units
        integer,             intent(in)    :: tsem, cons   ! FIELD_INSTANT|MEAN|ACCUM, FIELD_INTENSIVE|EXTENSIVE
        integer,             intent(out)   :: h
        character(len=*), intent(in), optional :: grid
        integer,          intent(in), optional :: dims(:)
        character(len=*), intent(in), optional :: zdim
        logical,          intent(in), optional :: masked

        character(len=*), parameter :: proc = "fesm_fields_add"
        type(field_t), allocatable :: tmp(:)
        type(field_t) :: f
        character(len=256) :: key
        integer :: i
        logical :: has_grid

        h = 0

        if (st%locked) then
            call fesm_error(proc, "state is locked", ctx_name(st, name)//new_line("a")// &
                "every field must be added before fields_lock")
        end if

        ! --- Names and strings ------------------------------------------------
        if (len_trim(name) == 0 .or. len_trim(name) > FIELD_NAME_LEN) then
            call bad_decl(st, name, "name = '"//trim(name)//"' (1 to "// &
                to_str(FIELD_NAME_LEN)//" characters)")
        end if
        do i = 1, st%n
            if (trim(st%f(i)%name) == trim(name)) then
                call fesm_error(proc, "duplicate field name", ctx_name(st, name)//new_line("a")// &
                    "existing handle = "//to_str(i))
            end if
        end do
        if (len_trim(units) == 0 .or. len_trim(units) > FIELD_UNITS_LEN) then
            call bad_decl(st, name, "units = '"//trim(units)//"' (1 to "// &
                to_str(FIELD_UNITS_LEN)//" characters; '1' for dimensionless)")
        end if

        ! --- Semantics ----------------------------------------------------------
        if (tsem /= FIELD_INSTANT .and. tsem /= FIELD_MEAN .and. tsem /= FIELD_ACCUM) then
            call bad_decl(st, name, "tsem = "//to_str(tsem)// &
                " (FIELD_INSTANT, FIELD_MEAN or FIELD_ACCUM)")
        end if
        if (cons /= FIELD_INTENSIVE .and. cons /= FIELD_EXTENSIVE) then
            call bad_decl(st, name, "cons = "//to_str(cons)// &
                " (FIELD_INTENSIVE or FIELD_EXTENSIVE)")
        end if

        ! --- Shape --------------------------------------------------------------
        f%rank = 0
        if (present(dims)) then
            if (size(dims) /= 2 .and. size(dims) /= 3) then
                call bad_decl(st, name, "dims = ["//to_str(dims)//"] (size 2 or 3)")
            end if
            if (any(dims <= 0)) then
                call bad_decl(st, name, "dims = ["//to_str(dims)//"] (every extent > 0)")
            end if
            f%rank = size(dims)
            f%nx = dims(1)
            f%ny = dims(2)
            if (f%rank == 3) f%nz = dims(3)
        end if

        has_grid = .false.
        if (present(grid)) has_grid = (len_trim(grid) > 0)
        if (f%rank >= 2 .and. .not. has_grid) then
            call bad_decl(st, name, "grid missing (required for rank "//to_str(f%rank)//")")
        end if
        if (f%rank == 0 .and. has_grid) then
            call bad_decl(st, name, "grid = '"//trim(grid)//"' given for a scalar")
        end if
        if (has_grid) then
            if (len_trim(grid) > FIELD_GRID_LEN) then
                call bad_decl(st, name, "grid = '"//trim(grid)//"' (at most "// &
                    to_str(FIELD_GRID_LEN)//" characters)")
            end if
            f%grid = grid
        end if

        if (present(zdim)) then
            if (f%rank /= 3) then
                call bad_decl(st, name, "zdim = '"//trim(zdim)//"' given for rank "//to_str(f%rank))
            end if
            if (len_trim(zdim) == 0 .or. len_trim(zdim) > FIELD_NAME_LEN) then
                call bad_decl(st, name, "zdim = '"//trim(zdim)//"' (1 to "// &
                    to_str(FIELD_NAME_LEN)//" characters)")
            end if
            f%zdim = zdim
        else if (f%rank == 3) then
            f%zdim = "z"
        end if

        if (present(masked)) then
            if (masked .and. f%rank == 0) then
                call bad_decl(st, name, "masked = T given for a scalar")
            end if
            f%masked = masked
        end if

        f%name  = name
        f%units = units
        key     = units_canonical(units)     ! short: at most the length of units or of a canonical key
        f%ukey  = key(1:FIELD_UNITS_LEN)
        f%tsem  = tsem
        f%cons  = cons

        ! --- Append (no buffers yet: growth copies metadata only) ------------
        if (.not. allocated(st%f)) allocate(st%f(CAPACITY0))
        if (st%n == size(st%f)) then
            allocate(tmp(max(CAPACITY0, 2*size(st%f))))
            tmp(1:st%n) = st%f(1:st%n)
            call move_alloc(tmp, st%f)
        end if

        st%n = st%n + 1
        st%f(st%n) = f
        h = st%n

    end subroutine fields_add


    subroutine fields_lock(st)
        ! End the declaration: trim the list to n fields and allocate every
        ! buffer once (v; acc for MEAN and ACCUM; masks for masked fields).
        ! v and acc start at zero, acc_valid at .true. (an empty window).

        type(field_state_t), intent(inout) :: st

        type(field_t), allocatable :: tmp(:)
        integer :: h

        if (st%locked) then
            call fesm_error("fesm_fields_lock", "state is locked", &
                "owner = '"//trim(st%owner)//"'"//new_line("a")//"fields_lock called twice")
        end if

        if (.not. allocated(st%f)) allocate(st%f(0))
        if (size(st%f) /= st%n) then
            allocate(tmp(st%n))
            tmp = st%f(1:st%n)
            call move_alloc(tmp, st%f)
        end if

        do h = 1, st%n
            associate (f => st%f(h))
                allocate(f%v(f%nx, f%ny, f%nz))
                f%v = 0.0_dp
                if (f%masked) then
                    allocate(f%valid(f%nx, f%ny))
                    f%valid = .false.
                end if
                if (f%tsem /= FIELD_INSTANT) then
                    allocate(f%acc(f%nx, f%ny, f%nz))
                    f%acc = 0.0_dp
                    if (f%masked) then
                        allocate(f%acc_valid(f%nx, f%ny))
                        f%acc_valid = .true.
                    end if
                end if
            end associate
        end do

        st%locked = .true.

    end subroutine fields_lock


    integer function fields_handle(st, name) result(h)
        ! Handle of the field `name`; stops if there is none.

        type(field_state_t), intent(in) :: st
        character(len=*),    intent(in) :: name

        h = find_field(st, name)
        if (h == 0) then
            call fesm_error("fesm_fields_handle", "no field named '"//trim(name)//"'", &
                "owner = '"//trim(st%owner)//"'"//new_line("a")// &
                "declared: "//declared_names(st))
        end if

    end function fields_handle


    logical function fields_has(st, name)
        ! Is a field `name` declared? Never stops.

        type(field_state_t), intent(in) :: st
        character(len=*),    intent(in) :: name

        fields_has = (find_field(st, name) > 0)

    end function fields_has

    ! =========================================================================
    ! fields_put: copy in
    ! =========================================================================
    !
    ! INSTANT:     v = var, [t0,t1] as given (t0 defaults to t1), published.
    ! MEAN, ACCUM: acc += var*(t1-t0), wsum += t1-t0, acc_valid AND= valid;
    !              requires t0 < t1 and, unless the window is empty,
    !              t0 == ta1 exactly. v is untouched until
    !              fields_finalize_interval.
    ! A rank-3 var with size(var,3) == 1 is accepted for a rank-2 field. sp
    ! input is widened with real(var, dp). Values at valid cells must be
    ! finite; values at invalid cells are stored as given (INSTANT) or not
    ! accumulated (MEAN, ACCUM) and carry no meaning.

    subroutine fields_put_0d_dp(st, h, var, t1, t0)
        type(field_state_t), intent(inout) :: st
        integer,  intent(in) :: h
        real(dp), intent(in) :: var
        real(dp), intent(in) :: t1                  ! end of validity [yr]
        real(dp), intent(in), optional :: t0        ! start; default t1 (INSTANT only)

        real(dp) :: ta, buf(1,1,1)
        integer  :: nbad, ibad(3)

        call put_check(st, h, 0, [1,1,1], .false., [0,0], t1, t0, ta)
        buf = var
        call put_copy_dp(st%f(h), 1, 1, 1, buf, ta, t1, nbad=nbad, ibad=ibad)
        call put_finish(st, h, ta, t1, nbad, ibad)

    end subroutine fields_put_0d_dp

    subroutine fields_put_0d_sp(st, h, var, t1, t0)
        type(field_state_t), intent(inout) :: st
        integer,  intent(in) :: h
        real(sp), intent(in) :: var
        real(dp), intent(in) :: t1
        real(dp), intent(in), optional :: t0

        real(dp) :: ta
        real(sp) :: buf(1,1,1)
        integer  :: nbad, ibad(3)

        call put_check(st, h, 0, [1,1,1], .false., [0,0], t1, t0, ta)
        buf = var
        call put_copy_sp(st%f(h), 1, 1, 1, buf, ta, t1, nbad=nbad, ibad=ibad)
        call put_finish(st, h, ta, t1, nbad, ibad)

    end subroutine fields_put_0d_sp

    subroutine fields_put_2d_dp(st, h, var, t1, t0, valid)
        type(field_state_t), intent(inout) :: st
        integer,  intent(in) :: h
        real(dp), intent(in) :: var(:,:)
        real(dp), intent(in) :: t1
        real(dp), intent(in), optional :: t0
        logical,  intent(in), optional :: valid(:,:)   ! required iff the field is masked

        real(dp) :: ta
        integer  :: nbad, ibad(3), vs(2)

        vs = 0
        if (present(valid)) vs = shape(valid)
        call put_check(st, h, 2, [size(var,1), size(var,2), 1], present(valid), vs, t1, t0, ta)
        call put_copy_dp(st%f(h), size(var,1), size(var,2), 1, var, ta, t1, valid, nbad, ibad)
        call put_finish(st, h, ta, t1, nbad, ibad)

    end subroutine fields_put_2d_dp

    subroutine fields_put_2d_sp(st, h, var, t1, t0, valid)
        type(field_state_t), intent(inout) :: st
        integer,  intent(in) :: h
        real(sp), intent(in) :: var(:,:)
        real(dp), intent(in) :: t1
        real(dp), intent(in), optional :: t0
        logical,  intent(in), optional :: valid(:,:)

        real(dp) :: ta
        integer  :: nbad, ibad(3), vs(2)

        vs = 0
        if (present(valid)) vs = shape(valid)
        call put_check(st, h, 2, [size(var,1), size(var,2), 1], present(valid), vs, t1, t0, ta)
        call put_copy_sp(st%f(h), size(var,1), size(var,2), 1, var, ta, t1, valid, nbad, ibad)
        call put_finish(st, h, ta, t1, nbad, ibad)

    end subroutine fields_put_2d_sp

    subroutine fields_put_3d_dp(st, h, var, t1, t0, valid)
        type(field_state_t), intent(inout) :: st
        integer,  intent(in) :: h
        real(dp), intent(in) :: var(:,:,:)
        real(dp), intent(in) :: t1
        real(dp), intent(in), optional :: t0
        logical,  intent(in), optional :: valid(:,:)

        real(dp) :: ta
        integer  :: nbad, ibad(3), vs(2)

        vs = 0
        if (present(valid)) vs = shape(valid)
        call put_check(st, h, 3, shape(var), present(valid), vs, t1, t0, ta)
        call put_copy_dp(st%f(h), size(var,1), size(var,2), size(var,3), var, ta, t1, valid, nbad, ibad)
        call put_finish(st, h, ta, t1, nbad, ibad)

    end subroutine fields_put_3d_dp

    subroutine fields_put_3d_sp(st, h, var, t1, t0, valid)
        type(field_state_t), intent(inout) :: st
        integer,  intent(in) :: h
        real(sp), intent(in) :: var(:,:,:)
        real(dp), intent(in) :: t1
        real(dp), intent(in), optional :: t0
        logical,  intent(in), optional :: valid(:,:)

        real(dp) :: ta
        integer  :: nbad, ibad(3), vs(2)

        vs = 0
        if (present(valid)) vs = shape(valid)
        call put_check(st, h, 3, shape(var), present(valid), vs, t1, t0, ta)
        call put_copy_sp(st%f(h), size(var,1), size(var,2), size(var,3), var, ta, t1, valid, nbad, ibad)
        call put_finish(st, h, ta, t1, nbad, ibad)

    end subroutine fields_put_3d_sp


    subroutine put_check(st, h, rk, shp, has_valid, vshp, t1, t0, ta)
        ! Every check of a put except finiteness (fused with the copy).
        ! Returns the start of the interval, ta.

        type(field_state_t), intent(in) :: st
        integer,  intent(in)  :: h, rk, shp(3), vshp(2)
        logical,  intent(in)  :: has_valid
        real(dp), intent(in)  :: t1
        real(dp), intent(in), optional :: t0
        real(dp), intent(out) :: ta

        character(len=*), parameter :: proc = "fesm_fields_put"

        call check_access(st, h, proc)
        call check_shape(st, h, proc, rk, shp)

        associate (f => st%f(h))

            if (f%masked .and. .not. has_valid) then
                call fesm_error(proc, "valid mask required", ctx(st, h))
            end if
            if (has_valid .and. .not. f%masked) then
                call fesm_error(proc, "field has no mask", ctx(st, h))
            end if
            if (has_valid) then
                if (vshp(1) /= f%nx .or. vshp(2) /= f%ny) then
                    call fesm_error(proc, "shape mismatch", ctx(st, h)//new_line("a")// &
                        "declared mask [nx,ny] = ["//to_str([f%nx, f%ny])//"]"//new_line("a")// &
                        "given mask shape      = ["//to_str(vshp)//"]")
                end if
            end if

            if (f%tsem == FIELD_INSTANT) then
                ta = t1
                if (present(t0)) ta = t0
                if (ta > t1) then
                    call fesm_error(proc, "invalid interval", ctx(st, h)//new_line("a")// &
                        "t0 = "//t_str(ta)//new_line("a")//"t1 = "//t_str(t1)//" (need t0 <= t1)")
                end if
            else
                if (.not. present(t0)) then
                    call fesm_error(proc, "invalid interval", ctx(st, h)//new_line("a")// &
                        "t0 absent (required for MEAN and ACCUM)"//new_line("a")// &
                        "t1 = "//t_str(t1))
                end if
                ta = t0
                if (.not. (ta < t1)) then
                    call fesm_error(proc, "invalid interval", ctx(st, h)//new_line("a")// &
                        "t0 = "//t_str(ta)//new_line("a")//"t1 = "//t_str(t1)//" (need t0 < t1)")
                end if
                if (f%wsum > 0.0_dp .and. (ta < f%ta1 .or. ta > f%ta1)) then
                    call fesm_error(proc, "gap or overlap in accumulation", ctx(st, h)//new_line("a")// &
                        "window = ["//t_str(f%ta0)//", "//t_str(f%ta1)//"]"//new_line("a")// &
                        "given  = ["//t_str(ta)//", "//t_str(t1)//"]"//new_line("a")// &
                        merge("t0 > ta1: a gap  ", "t0 < ta1: overlap", ta > f%ta1))
                end if
            end if

        end associate

    end subroutine put_check


    subroutine put_copy_dp(f, nx, ny, nz, var, t0, t1, valid, nbad, ibad)
        ! The copy (INSTANT) or accumulation (MEAN, ACCUM) with the
        ! finiteness check fused in. nbad: non-finite values at valid cells,
        ! ibad: the first of them. Finite is abs(x) <= huge(x): NaN fails
        ! (a comparison with NaN is false), +-Inf exceed huge. Not
        ! ieee_is_finite, which ifx calls out of line for every value.

        type(field_t), intent(inout) :: f
        integer,  intent(in)  :: nx, ny, nz
        real(dp), intent(in)  :: var(nx,ny,nz)
        real(dp), intent(in)  :: t0, t1
        logical,  intent(in), optional :: valid(nx,ny)
        integer,  intent(out) :: nbad, ibad(3)

        real(dp) :: x, dt
        integer  :: i, j, k

        nbad = 0
        ibad = 0

        if (f%tsem == FIELD_INSTANT) then
            do k = 1, nz
            do j = 1, ny
            do i = 1, nx
                x = var(i,j,k)
                f%v(i,j,k) = x
                if (.not. (abs(x) <= huge(x))) then
                    if (present(valid)) then
                        if (.not. valid(i,j)) cycle
                    end if
                    nbad = nbad + 1
                    if (nbad == 1) ibad = [i, j, k]
                end if
            end do
            end do
            end do
            if (present(valid)) f%valid = valid
        else
            dt = t1 - t0
            do k = 1, nz
            do j = 1, ny
            do i = 1, nx
                if (present(valid)) then
                    if (.not. valid(i,j)) cycle
                end if
                x = var(i,j,k)
                if (.not. (abs(x) <= huge(x))) then
                    nbad = nbad + 1
                    if (nbad == 1) ibad = [i, j, k]
                    cycle
                end if
                f%acc(i,j,k) = f%acc(i,j,k) + x*dt
            end do
            end do
            end do
            if (present(valid)) f%acc_valid = f%acc_valid .and. valid
        end if

    end subroutine put_copy_dp

    subroutine put_copy_sp(f, nx, ny, nz, var, t0, t1, valid, nbad, ibad)
        ! As put_copy_dp, widening each value with real(var, dp) (every finite
        ! sp value is finite in dp, Inf and NaN stay so).

        type(field_t), intent(inout) :: f
        integer,  intent(in)  :: nx, ny, nz
        real(sp), intent(in)  :: var(nx,ny,nz)
        real(dp), intent(in)  :: t0, t1
        logical,  intent(in), optional :: valid(nx,ny)
        integer,  intent(out) :: nbad, ibad(3)

        real(dp) :: x, dt
        integer  :: i, j, k

        nbad = 0
        ibad = 0

        if (f%tsem == FIELD_INSTANT) then
            do k = 1, nz
            do j = 1, ny
            do i = 1, nx
                x = real(var(i,j,k), dp)
                f%v(i,j,k) = x
                if (.not. (abs(x) <= huge(x))) then
                    if (present(valid)) then
                        if (.not. valid(i,j)) cycle
                    end if
                    nbad = nbad + 1
                    if (nbad == 1) ibad = [i, j, k]
                end if
            end do
            end do
            end do
            if (present(valid)) f%valid = valid
        else
            dt = t1 - t0
            do k = 1, nz
            do j = 1, ny
            do i = 1, nx
                if (present(valid)) then
                    if (.not. valid(i,j)) cycle
                end if
                x = real(var(i,j,k), dp)
                if (.not. (abs(x) <= huge(x))) then
                    nbad = nbad + 1
                    if (nbad == 1) ibad = [i, j, k]
                    cycle
                end if
                f%acc(i,j,k) = f%acc(i,j,k) + x*dt
            end do
            end do
            end do
            if (present(valid)) f%acc_valid = f%acc_valid .and. valid
        end if

    end subroutine put_copy_sp


    subroutine put_finish(st, h, ta, t1, nbad, ibad)
        ! Stop on non-finite values, else record the interval: INSTANT
        ! publishes, MEAN/ACCUM extend the window.

        type(field_state_t), intent(inout) :: st
        integer,  intent(in) :: h, nbad, ibad(3)
        real(dp), intent(in) :: ta, t1

        if (nbad > 0) then
            call fesm_error("fesm_fields_put", "non-finite value", ctx(st, h)//new_line("a")// &
                "count = "//to_str(nbad)//" (at valid cells)"//new_line("a")// &
                "first at (i,j,k) = ("//to_str(ibad)//")")
        end if

        associate (f => st%f(h))
            if (f%tsem == FIELD_INSTANT) then
                f%t0 = ta
                f%t1 = t1
                f%published = .true.
            else
                if (.not. (f%wsum > 0.0_dp)) f%ta0 = ta
                f%ta1  = t1
                f%wsum = f%wsum + (t1 - ta)
            end if
        end associate

    end subroutine put_finish

    ! =========================================================================
    ! fields_get: copy out the published value
    ! =========================================================================
    !
    ! The open window is invisible. dp to sp narrows with real(v, sp) (no
    ! range check). Stops if nothing was published. A rank-3 var with
    ! size(var,3) == 1 is accepted for a rank-2 field, as in fields_put.

    subroutine fields_get_0d_dp(st, h, var, t0, t1)
        type(field_state_t), intent(in) :: st
        integer,  intent(in)  :: h
        real(dp), intent(out) :: var
        real(dp), intent(out), optional :: t0, t1   ! validity of the value returned

        call get_check(st, h, 0, [1,1,1], .false., [0,0])
        var = st%f(h)%v(1,1,1)
        if (present(t0)) t0 = st%f(h)%t0
        if (present(t1)) t1 = st%f(h)%t1

    end subroutine fields_get_0d_dp

    subroutine fields_get_0d_sp(st, h, var, t0, t1)
        type(field_state_t), intent(in) :: st
        integer,  intent(in)  :: h
        real(sp), intent(out) :: var
        real(dp), intent(out), optional :: t0, t1

        call get_check(st, h, 0, [1,1,1], .false., [0,0])
        var = real(st%f(h)%v(1,1,1), sp)
        if (present(t0)) t0 = st%f(h)%t0
        if (present(t1)) t1 = st%f(h)%t1

    end subroutine fields_get_0d_sp

    subroutine fields_get_2d_dp(st, h, var, t0, t1, valid)
        type(field_state_t), intent(in) :: st
        integer,  intent(in)  :: h
        real(dp), intent(out) :: var(:,:)           ! caller's array; shape must match
        real(dp), intent(out), optional :: t0, t1
        logical,  intent(out), optional :: valid(:,:)

        integer :: vs(2)

        vs = 0
        if (present(valid)) vs = shape(valid)
        call get_check(st, h, 2, [size(var,1), size(var,2), 1], present(valid), vs)
        call get_copy_dp(st%f(h), size(var,1), size(var,2), 1, var)
        if (present(valid)) valid = st%f(h)%valid
        if (present(t0)) t0 = st%f(h)%t0
        if (present(t1)) t1 = st%f(h)%t1

    end subroutine fields_get_2d_dp

    subroutine fields_get_2d_sp(st, h, var, t0, t1, valid)
        type(field_state_t), intent(in) :: st
        integer,  intent(in)  :: h
        real(sp), intent(out) :: var(:,:)
        real(dp), intent(out), optional :: t0, t1
        logical,  intent(out), optional :: valid(:,:)

        integer :: vs(2)

        vs = 0
        if (present(valid)) vs = shape(valid)
        call get_check(st, h, 2, [size(var,1), size(var,2), 1], present(valid), vs)
        call get_copy_sp(st%f(h), size(var,1), size(var,2), 1, var)
        if (present(valid)) valid = st%f(h)%valid
        if (present(t0)) t0 = st%f(h)%t0
        if (present(t1)) t1 = st%f(h)%t1

    end subroutine fields_get_2d_sp

    subroutine fields_get_3d_dp(st, h, var, t0, t1, valid)
        type(field_state_t), intent(in) :: st
        integer,  intent(in)  :: h
        real(dp), intent(out) :: var(:,:,:)
        real(dp), intent(out), optional :: t0, t1
        logical,  intent(out), optional :: valid(:,:)

        integer :: vs(2)

        vs = 0
        if (present(valid)) vs = shape(valid)
        call get_check(st, h, 3, shape(var), present(valid), vs)
        call get_copy_dp(st%f(h), size(var,1), size(var,2), size(var,3), var)
        if (present(valid)) valid = st%f(h)%valid
        if (present(t0)) t0 = st%f(h)%t0
        if (present(t1)) t1 = st%f(h)%t1

    end subroutine fields_get_3d_dp

    subroutine fields_get_3d_sp(st, h, var, t0, t1, valid)
        type(field_state_t), intent(in) :: st
        integer,  intent(in)  :: h
        real(sp), intent(out) :: var(:,:,:)
        real(dp), intent(out), optional :: t0, t1
        logical,  intent(out), optional :: valid(:,:)

        integer :: vs(2)

        vs = 0
        if (present(valid)) vs = shape(valid)
        call get_check(st, h, 3, shape(var), present(valid), vs)
        call get_copy_sp(st%f(h), size(var,1), size(var,2), size(var,3), var)
        if (present(valid)) valid = st%f(h)%valid
        if (present(t0)) t0 = st%f(h)%t0
        if (present(t1)) t1 = st%f(h)%t1

    end subroutine fields_get_3d_sp


    subroutine get_check(st, h, rk, shp, has_valid, vshp)
        type(field_state_t), intent(in) :: st
        integer, intent(in) :: h, rk, shp(3), vshp(2)
        logical, intent(in) :: has_valid

        character(len=*), parameter :: proc = "fesm_fields_get"

        call check_access(st, h, proc)

        associate (f => st%f(h))
            if (.not. f%published) then
                call fesm_error(proc, "field read before any value was published", ctx(st, h))
            end if
            call check_shape(st, h, proc, rk, shp)
            if (has_valid .and. .not. f%masked) then
                call fesm_error(proc, "field has no mask", ctx(st, h))
            end if
            if (has_valid) then
                if (vshp(1) /= f%nx .or. vshp(2) /= f%ny) then
                    call fesm_error(proc, "shape mismatch", ctx(st, h)//new_line("a")// &
                        "declared mask [nx,ny] = ["//to_str([f%nx, f%ny])//"]"//new_line("a")// &
                        "given mask shape      = ["//to_str(vshp)//"]")
                end if
            end if
        end associate

    end subroutine get_check

    subroutine get_copy_dp(f, nx, ny, nz, var)
        type(field_t), intent(in)  :: f
        integer,       intent(in)  :: nx, ny, nz
        real(dp),      intent(out) :: var(nx,ny,nz)

        var = f%v

    end subroutine get_copy_dp

    subroutine get_copy_sp(f, nx, ny, nz, var)
        type(field_t), intent(in)  :: f
        integer,       intent(in)  :: nx, ny, nz
        real(sp),      intent(out) :: var(nx,ny,nz)

        var = real(f%v, sp)

    end subroutine get_copy_sp

    ! =========================================================================
    ! fields_finalize_interval
    ! =========================================================================

    subroutine fields_finalize_interval(st, t)
        ! Close the accumulation window of every MEAN and ACCUM field at `t`:
        ! v = acc/wsum (MEAN) or acc (ACCUM), valid = acc_valid,
        ! [t0,t1] = [ta0,ta1], published; then the window is cleared. Each
        ! window must be non-empty and end exactly at t. INSTANT fields are
        ! untouched.

        type(field_state_t), intent(inout) :: st
        real(dp),            intent(in)    :: t      ! end of the interval being closed

        character(len=*), parameter :: proc = "fesm_fields_finalize_interval"
        integer :: h

        if (.not. st%locked) then
            call fesm_error(proc, "state is not locked", "owner = '"//trim(st%owner)//"'")
        end if

        do h = 1, st%n
            associate (f => st%f(h))
                if (f%tsem == FIELD_INSTANT) cycle

                if (.not. (f%wsum > 0.0_dp) .or. f%ta1 < t .or. f%ta1 > t) then
                    call fesm_error(proc, "accumulation does not reach t", ctx(st, h)//new_line("a")// &
                        "ta1  = "//t_str(f%ta1)//new_line("a")// &
                        "t    = "//t_str(t)//new_line("a")// &
                        "wsum = "//t_str(f%wsum))
                end if

                if (f%tsem == FIELD_MEAN) then
                    f%v = f%acc / f%wsum
                else
                    f%v = f%acc
                end if
                if (f%masked) f%valid = f%acc_valid
                f%t0 = f%ta0
                f%t1 = f%ta1
                f%published = .true.

                f%acc  = 0.0_dp
                if (f%masked) f%acc_valid = .true.
                f%wsum = 0.0_dp
                f%ta0  = 0.0_dp
                f%ta1  = 0.0_dp
            end associate
        end do

    end subroutine fields_finalize_interval

    ! =========================================================================
    ! Units
    ! =========================================================================

    pure function units_key(units) result(key)
        ! Units in lower case without blanks and exponent markers ("*", "^"):
        ! "kg m**-2 d**-1" => "kgm-2d-1".

        character(len=*), intent(in) :: units
        character(len=256) :: key

        integer :: i, n, c

        key = ""
        n = 0
        do i = 1, len_trim(units)
            if (units(i:i) == " " .or. units(i:i) == "*" .or. units(i:i) == "^") cycle
            c = iachar(units(i:i))
            if (c >= iachar("A") .and. c <= iachar("Z")) c = c + 32
            n = n + 1
            key(n:n) = achar(c)
        end do

    end function units_key


    pure function units_canonical(units) result(key)
        ! units_key of the canonical spelling: an accepted alias (see
        ! UNITS_ALIAS) is replaced by its canonical form, anything else is
        ! its own units_key. "mm w.e./yr" => "kgm-2yr-1".

        character(len=*), intent(in) :: units
        character(len=256) :: key

        integer :: i

        key = units_key(units)
        do i = 1, N_ALIAS
            if (key == units_key(UNITS_ALIAS(2,i))) then
                key = units_key(UNITS_ALIAS(1,i))
                return
            end if
        end do

    end function units_canonical


    pure logical function units_equivalent(a, b)
        ! Are `a` and `b` the same units, possibly spelled differently?

        character(len=*), intent(in) :: a, b

        units_equivalent = (units_canonical(a) == units_canonical(b))

    end function units_equivalent


    subroutine fields_check_link(src, hs, dst, hd)
        ! Stop unless field hs of src can feed field hd of dst: equivalent
        ! units, same rank and nz, same conservation flag, and time
        ! semantics allowed by
        !
        !     source \ destination   INSTANT   MEAN   ACCUM
        !     INSTANT                ok        -      -
        !     MEAN                   ok        ok     -
        !     ACCUM                  -         -      ok
        !
        ! The horizontal shape and the grid may differ (the exchange remaps).

        type(field_state_t), intent(in) :: src, dst
        integer,             intent(in) :: hs, hd

        character(len=*), parameter :: proc = "fesm_fields_check_link"
        character(len=:), allocatable :: detail
        logical :: ok

        call check_handle(src, hs, proc)
        call check_handle(dst, hd, proc)

        associate (fs => src%f(hs), fd => dst%f(hd))

            detail = link_side("source     ", src, hs)//new_line("a")// &
                     link_side("destination", dst, hd)

            if (fs%ukey /= fd%ukey) then
                call fesm_error(proc, "units differ and no equivalence is known", detail)
            end if
            if (fs%rank /= fd%rank) then
                call fesm_error(proc, "rank differs", detail)
            end if
            if (fs%nz /= fd%nz) then
                call fesm_error(proc, "nz differs", detail)
            end if
            if (fs%cons /= fd%cons) then
                call fesm_error(proc, "conservation flag differs", detail)
            end if
            ok = (fs%tsem == fd%tsem) .or. (fs%tsem == FIELD_MEAN .and. fd%tsem == FIELD_INSTANT)
            if (.not. ok) then
                call fesm_error(proc, "time semantics incompatible", detail)
            end if

        end associate

    end subroutine fields_check_link


    function link_side(label, st, h) result(s)
        character(len=*),    intent(in) :: label
        type(field_state_t), intent(in) :: st
        integer,             intent(in) :: h
        character(len=:), allocatable :: s

        associate (f => st%f(h))
            s = label//": owner = '"//trim(st%owner)//"', field = '"//trim(f%name)// &
                "', units = '"//trim(f%units)//"' (key '"//trim(f%ukey)//"'), "// &
                "shape = "//shape_str(f)//", "//trim(tsem_name(f%tsem))//", "// &
                trim(cons_name(f%cons))
        end associate

    end function link_side

    ! =========================================================================
    ! fields_describe
    ! =========================================================================

    subroutine fields_describe(st, level)
        ! A header ("[ice_grl:exports] 14 fields, 6.9 MB") and one line per
        ! field: handle, name, units, grid, shape, tsem, cons, validity (or
        ! "unpublished"), open window, mask. Through log_msg at `level`
        ! (default LOG_INFO).

        type(field_state_t), intent(in) :: st
        integer, intent(in), optional :: level

        character(len=*), parameter :: who = "fesm_fields"
        character(len=:), allocatable :: line
        character(len=32) :: buf
        integer  :: lvl, h
        real(dp) :: bytes, nval

        lvl = LOG_INFO
        if (present(level)) lvl = level

        bytes = 0.0_dp
        do h = 1, st%n
            associate (f => st%f(h))
                nval  = real(f%nx, dp)*real(f%ny, dp)*real(f%nz, dp)
                bytes = bytes + 8.0_dp*nval
                if (f%tsem /= FIELD_INSTANT) bytes = bytes + 8.0_dp*nval
                if (f%masked) then
                    nval  = real(f%nx, dp)*real(f%ny, dp)*real(storage_size(.true.)/8, dp)
                    bytes = bytes + nval
                    if (f%tsem /= FIELD_INSTANT) bytes = bytes + nval
                end if
            end associate
        end do
        write(buf, "(f12.1)") bytes/1.0e6_dp

        line = "["//trim(st%owner)//"] "//to_str(st%n)//" fields, "//trim(adjustl(buf))//" MB"
        if (.not. st%locked) line = line//" (not locked)"
        call log_msg(lvl, who, line)

        do h = 1, st%n
            associate (f => st%f(h))
                line = "["//trim(st%owner)//"] "//pad(to_str(h), 3)//" "// &
                       pad(trim(f%name), 12)//" "//pad("["//trim(f%units)//"]", 16)//" "// &
                       pad(trim(f%grid), 12)//" "//pad(shape_str(f), 18)//" "// &
                       pad(trim(tsem_name(f%tsem)), 7)//" "//pad(trim(cons_name(f%cons)), 9)//" "
                if (f%published) then
                    line = line//"["//to_str(f%t0)//", "//to_str(f%t1)//"]"
                else
                    line = line//"unpublished"
                end if
                if (f%wsum > 0.0_dp) then
                    line = line//", window ["//to_str(f%ta0)//", "//to_str(f%ta1)//"]"
                end if
                if (f%masked) line = line//", masked"
                call log_msg(lvl, who, line)
            end associate
        end do

    end subroutine fields_describe

    ! =========================================================================
    ! Restart
    ! =========================================================================
    !
    ! One NetCDF file per state, owned entirely by this module. Per field:
    ! `<name>` (v, double), `<name>_acc` (MEAN, ACCUM), `<name>_valid` and
    ! `<name>_acc_valid` (int 0/1, masked fields); attributes on `<name>`:
    ! units, tsem, cons, grid (rank >= 2), zdim (rank 3), published, t0, t1,
    ! ta0, ta1, wsum. Dimensions `<grid>_x`, `<grid>_y`, `<zdim>`; scalars on
    ! the length-1 dimension `one`. Global attributes: fesm_fields_owner,
    ! fesm_fields_n and fesm_fields_names (the declared names, comma
    ! separated), against which fields_restart_read detects a renamed,
    ! dropped or added field. Each call opens the file once and does all its
    ! reads or writes through that handle (every open and close is a file
    ! system round trip, slow on a parallel file system).

    subroutine fields_restart_write(st, filename)
        ! Write the whole state to `filename` (created or overwritten).

        type(field_state_t), intent(in) :: st
        character(len=*),    intent(in) :: filename

        character(len=*), parameter :: proc = "fesm_fields_restart_write"
        character(len=FIELD_GRID_LEN+2), allocatable :: dname(:)
        integer, allocatable :: dlen(:)
        character(len=FIELD_GRID_LEN+2) :: dn(3)
        integer :: h, nd, k, q, ncid

        if (.not. st%locked) then
            call fesm_error(proc, "state is not locked", "owner = '"//trim(st%owner)//"'"// &
                new_line("a")//"file  = "//trim(filename))
        end if

        call check_restart_names(st, filename)

        ! Dimensions, each once; a name used with two lengths is an error.
        allocate(dname(3*st%n), dlen(3*st%n))
        nd = 0
        do h = 1, st%n
            call restart_dims(st%f(h), dn)
            do k = 1, max(1, st%f(h)%rank)
                q = find_dim(dname(1:nd), dn(k))
                if (q > 0) then
                    if (dlen(q) /= dim_len(st%f(h), k)) then
                        call fesm_error(proc, "restart dimension length conflict", ctx(st, h)// &
                            new_line("a")//"dimension '"//trim(dn(k))//"' has length "// &
                            to_str(dlen(q))//" for an earlier field, "// &
                            to_str(dim_len(st%f(h), k))//" here"//new_line("a")// &
                            "file  = "//trim(filename))
                    end if
                else
                    nd = nd + 1
                    dname(nd) = dn(k)
                    dlen(nd)  = dim_len(st%f(h), k)
                end if
            end do
        end do

        call nc_create(filename)
        call nc_open(filename, ncid, writable=.true.)
        call nc_write_attr(filename, ncid, "fesm_fields_owner", trim(st%owner))
        call nc_write_attr(filename, ncid, "fesm_fields_n", st%n)
        if (st%n > 0) call nc_write_attr(filename, ncid, "fesm_fields_names", declared_names(st, ","))

        do q = 1, nd
            call nc_write_dim(filename, trim(dname(q)), 1.0_dp, dx=1.0_dp, nx=dlen(q), ncid=ncid)
        end do

        do h = 1, st%n
            associate (f => st%f(h))
                call restart_dims(f, dn)
                call write_real(filename, ncid, trim(f%name), f, f%v, dn)
                if (f%tsem /= FIELD_INSTANT) call write_real(filename, ncid, trim(f%name)//"_acc", f, f%acc, dn)
                if (f%masked) then
                    call nc_write(filename, trim(f%name)//"_valid", f%valid, &
                                  dim1=trim(dn(1)), dim2=trim(dn(2)), ncid=ncid)
                    if (f%tsem /= FIELD_INSTANT) &
                        call nc_write(filename, trim(f%name)//"_acc_valid", f%acc_valid, &
                                      dim1=trim(dn(1)), dim2=trim(dn(2)), ncid=ncid)
                end if

                call nc_write_attr(filename, ncid, trim(f%name), "units", trim(f%units))
                call nc_write_attr(filename, ncid, trim(f%name), "tsem",  f%tsem)
                call nc_write_attr(filename, ncid, trim(f%name), "cons",  f%cons)
                if (f%rank >= 2) call nc_write_attr(filename, ncid, trim(f%name), "grid", trim(f%grid))
                if (f%rank == 3) call nc_write_attr(filename, ncid, trim(f%name), "zdim", trim(f%zdim))
                call nc_write_attr(filename, ncid, trim(f%name), "published", merge(1, 0, f%published))
                call nc_write_attr(filename, ncid, trim(f%name), "t0",   f%t0)
                call nc_write_attr(filename, ncid, trim(f%name), "t1",   f%t1)
                call nc_write_attr(filename, ncid, trim(f%name), "ta0",  f%ta0)
                call nc_write_attr(filename, ncid, trim(f%name), "ta1",  f%ta1)
                call nc_write_attr(filename, ncid, trim(f%name), "wsum", f%wsum)
            end associate
        end do

        call nc_close(ncid)

    end subroutine fields_restart_write


    subroutine fields_restart_read(st, filename)
        ! Restore a state written by fields_restart_write into a locked state
        ! with the same declaration. Every declared field must be in the file
        ! with the same shape, units key, tsem, cons and mask, and the file
        ! must hold no other field; then values, masks, windows and validity
        ! are restored, so the state is bit-identical to the written one
        ! (up to ncio's read conventions: |x| < 1e-30 reads as 0, NaN as
        ! -9999; both only matter at invalid cells or for subnormal values).

        type(field_state_t), intent(inout) :: st
        character(len=*),    intent(in)    :: filename

        character(len=*), parameter :: proc = "fesm_fields_restart_read"
        character(len=:), allocatable :: names_file
        character(len=4096) :: buf
        character(len=256)  :: ustr
        integer, allocatable :: dl(:)
        integer :: h, nfile, ival, ncid
        logical :: exists

        if (.not. st%locked) then
            call fesm_error(proc, "state is not locked", "owner = '"//trim(st%owner)//"'"// &
                new_line("a")//"file  = "//trim(filename))
        end if

        inquire(file=trim(filename), exist=exists)
        if (.not. exists) then
            call fesm_error(proc, "restart file not found", "owner = '"//trim(st%owner)//"'"// &
                new_line("a")//"file  = "//trim(filename))
        end if

        call nc_open(filename, ncid, writable=.false.)

        ! --- The set of fields ------------------------------------------------
        if (.not. nc_exists_attr(filename, ncid, "fesm_fields_n")) then
            call mismatch(0, "global attribute fesm_fields_n", "a fesm_fields restart file", "absent")
        end if
        call nc_read_attr(filename, ncid, "fesm_fields_n", nfile)
        names_file = ""
        if (nfile > 0) then
            if (.not. nc_exists_attr(filename, ncid, "fesm_fields_names")) then
                call mismatch(0, "global attribute fesm_fields_names", "present", "absent")
            end if
            buf = ""
            call nc_read_attr(filename, ncid, "fesm_fields_names", buf)
            names_file = trim(buf)
        end if
        if (nfile /= st%n) then
            call mismatch(0, "number of fields", to_str(st%n)//" ("//declared_names(st)//")", &
                to_str(nfile)//" ("//names_file//")")
        end if
        do h = 1, st%n
            if (index(","//names_file//",", ","//trim(st%f(h)%name)//",") == 0) then
                call mismatch(h, "fields", declared_names(st), names_file)
            end if
        end do

        ! --- Every field: present, same shape and metadata --------------------
        do h = 1, st%n
            associate (f => st%f(h))

                if (.not. nc_exists_var(filename, trim(f%name), ncid=ncid)) then
                    call mismatch(h, "variable '"//trim(f%name)//"'", "present", "absent")
                end if

                call nc_dims(filename, trim(f%name), dims=dl, ncid=ncid)
                if (.not. same_shape(f, dl)) then
                    call mismatch(h, "shape", shape_str(f), "["//to_str(dl)//"]")
                end if

                ustr = ""
                if (nc_exists_attr(filename, ncid, trim(f%name), "units")) &
                    call nc_read_attr(filename, ncid, trim(f%name), "units", ustr)
                if (units_canonical(ustr) /= f%ukey) then
                    call mismatch(h, "units key", "'"//trim(f%ukey)//"' ('"//trim(f%units)//"')", &
                        "'"//trim(units_canonical(ustr))//"' ('"//trim(ustr)//"')")
                end if

                call read_int_attr(h, "tsem", ival)
                if (ival /= f%tsem) then
                    call mismatch(h, "tsem", trim(tsem_name(f%tsem)), trim(tsem_name(ival)))
                end if
                call read_int_attr(h, "cons", ival)
                if (ival /= f%cons) then
                    call mismatch(h, "cons", trim(cons_name(f%cons)), trim(cons_name(ival)))
                end if

                call expect_var(h, trim(f%name)//"_acc",       f%tsem /= FIELD_INSTANT)
                call expect_var(h, trim(f%name)//"_valid",     f%masked)
                call expect_var(h, trim(f%name)//"_acc_valid", f%masked .and. f%tsem /= FIELD_INSTANT)

            end associate
        end do

        ! --- Restore ----------------------------------------------------------
        do h = 1, st%n
            associate (f => st%f(h))

                call read_real(filename, ncid, trim(f%name), f, f%v)
                if (f%tsem /= FIELD_INSTANT) call read_real(filename, ncid, trim(f%name)//"_acc", f, f%acc)
                if (f%masked) then
                    call nc_read(filename, trim(f%name)//"_valid", f%valid, ncid=ncid)
                    if (f%tsem /= FIELD_INSTANT) &
                        call nc_read(filename, trim(f%name)//"_acc_valid", f%acc_valid, ncid=ncid)
                end if

                call read_int_attr(h, "published", ival)
                f%published = (ival /= 0)
                call nc_read_attr(filename, ncid, trim(f%name), "t0",   f%t0)
                call nc_read_attr(filename, ncid, trim(f%name), "t1",   f%t1)
                call nc_read_attr(filename, ncid, trim(f%name), "ta0",  f%ta0)
                call nc_read_attr(filename, ncid, trim(f%name), "ta1",  f%ta1)
                call nc_read_attr(filename, ncid, trim(f%name), "wsum", f%wsum)

            end associate
        end do

        call nc_close(ncid)

    contains

        subroutine mismatch(hh, what, expected, found)
            integer,          intent(in) :: hh
            character(len=*), intent(in) :: what, expected, found

            character(len=:), allocatable :: who

            if (hh > 0) then
                who = ctx(st, hh)
            else
                who = "owner = '"//trim(st%owner)//"'"
            end if
            call fesm_error(proc, "restart file does not match the declaration", &
                who//new_line("a")//"file     = "//trim(filename)//new_line("a")// &
                "item     = "//what//new_line("a")// &
                "expected = "//expected//new_line("a")// &
                "found    = "//found)

        end subroutine mismatch

        subroutine read_int_attr(hh, att, val)
            integer,          intent(in)  :: hh
            character(len=*), intent(in)  :: att
            integer,          intent(out) :: val

            if (.not. nc_exists_attr(filename, ncid, trim(st%f(hh)%name), att)) then
                call mismatch(hh, "attribute '"//att//"'", "present", "absent")
            end if
            call nc_read_attr(filename, ncid, trim(st%f(hh)%name), att, val)

        end subroutine read_int_attr

        subroutine expect_var(hh, vname, wanted)
            integer,          intent(in) :: hh
            character(len=*), intent(in) :: vname
            logical,          intent(in) :: wanted

            logical :: found

            found = nc_exists_var(filename, vname, ncid=ncid)
            if (found .neqv. wanted) then
                call mismatch(hh, "variable '"//vname//"'", &
                    merge("present", "absent ", wanted), merge("present", "absent ", found))
            end if

        end subroutine expect_var

    end subroutine fields_restart_read


    subroutine check_restart_names(st, filename)
        ! The variable and dimension names of the file must be distinct
        ! (a field "a_acc" next to a MEAN field "a" would collide).

        type(field_state_t), intent(in) :: st
        character(len=*),    intent(in) :: filename

        character(len=FIELD_GRID_LEN+12), allocatable :: nm(:)
        character(len=FIELD_GRID_LEN+2),  allocatable :: dd(:)
        character(len=FIELD_GRID_LEN+2) :: dn(3)
        integer :: h, k, m, md, i, j

        ! Variables of every field, then each distinct dimension once.
        allocate(nm(7*st%n), dd(3*st%n))
        m  = 0
        md = 0
        do h = 1, st%n
            associate (f => st%f(h))
                call add(trim(f%name))
                if (f%tsem /= FIELD_INSTANT) call add(trim(f%name)//"_acc")
                if (f%masked) call add(trim(f%name)//"_valid")
                if (f%masked .and. f%tsem /= FIELD_INSTANT) call add(trim(f%name)//"_acc_valid")
                call restart_dims(f, dn)
                do k = 1, max(1, f%rank)
                    if (find_dim(dd(1:md), dn(k)) == 0) then
                        md = md + 1
                        dd(md) = dn(k)
                    end if
                end do
            end associate
        end do
        do k = 1, md
            call add(trim(dd(k)))
        end do

        do i = 1, m
            do j = i+1, m
                if (nm(i) == nm(j)) then
                    call fesm_error("fesm_fields_restart_write", "restart variable name collision", &
                        "owner = '"//trim(st%owner)//"'"//new_line("a")// &
                        "name  = '"//trim(nm(i))//"' (a field, a companion _acc/_valid "// &
                        "variable or a dimension)"//new_line("a")// &
                        "file  = "//trim(filename))
                end if
            end do
        end do

    contains

        subroutine add(s)
            character(len=*), intent(in) :: s
            m = m + 1
            nm(m) = s
        end subroutine add

    end subroutine check_restart_names


    subroutine restart_dims(f, dn)
        ! Dimension names of a field in the restart file.
        type(field_t), intent(in) :: f
        character(len=FIELD_GRID_LEN+2), intent(out) :: dn(3)

        dn = ""
        if (f%rank == 0) then
            dn(1) = "one"
        else
            dn(1) = trim(f%grid)//"_x"
            dn(2) = trim(f%grid)//"_y"
            if (f%rank == 3) dn(3) = f%zdim
        end if

    end subroutine restart_dims

    integer function dim_len(f, k)
        type(field_t), intent(in) :: f
        integer,       intent(in) :: k

        select case (k)
            case (1);     dim_len = f%nx
            case (2);     dim_len = f%ny
            case default; dim_len = f%nz
        end select

    end function dim_len

    integer function find_dim(names, name)
        character(len=*), intent(in) :: names(:), name
        integer :: i

        find_dim = 0
        do i = 1, size(names)
            if (names(i) == name) then
                find_dim = i
                return
            end if
        end do

    end function find_dim

    logical function same_shape(f, dl)
        type(field_t), intent(in) :: f
        integer,       intent(in) :: dl(:)

        select case (f%rank)
            case (0)
                same_shape = (size(dl) == 1)
                if (same_shape) same_shape = (dl(1) == 1)
            case (2)
                same_shape = (size(dl) == 2)
                if (same_shape) same_shape = (dl(1) == f%nx .and. dl(2) == f%ny)
            case default
                same_shape = (size(dl) == 3)
                if (same_shape) same_shape = (dl(1) == f%nx .and. dl(2) == f%ny .and. dl(3) == f%nz)
        end select

    end function same_shape

    subroutine write_real(filename, ncid, vname, f, a, dn)
        ! One (nx,ny,nz) buffer of f, written with the rank of the field.
        character(len=*), intent(in) :: filename, vname
        integer,          intent(in) :: ncid
        type(field_t),    intent(in) :: f
        real(dp),         intent(in) :: a(:,:,:)
        character(len=*), intent(in) :: dn(3)

        select case (f%rank)
            case (0)
                call nc_write(filename, vname, a(:,1,1), dim1=trim(dn(1)), ncid=ncid)
            case (2)
                call nc_write(filename, vname, a(:,:,1), dim1=trim(dn(1)), dim2=trim(dn(2)), ncid=ncid)
            case default
                call nc_write(filename, vname, a, dim1=trim(dn(1)), dim2=trim(dn(2)), dim3=trim(dn(3)), ncid=ncid)
        end select

    end subroutine write_real

    subroutine read_real(filename, ncid, vname, f, a)
        ! Read into an existing (nx,ny,nz) buffer, with the rank of the field.
        character(len=*), intent(in)    :: filename, vname
        integer,          intent(in)    :: ncid
        type(field_t),    intent(in)    :: f
        real(dp),         intent(inout) :: a(:,:,:)

        select case (f%rank)
            case (0)
                call nc_read(filename, vname, a(:,1,1), ncid=ncid)
            case (2)
                call nc_read(filename, vname, a(:,:,1), ncid=ncid)
            case default
                call nc_read(filename, vname, a, ncid=ncid)
        end select

    end subroutine read_real

    ! =========================================================================
    ! Internal helpers
    ! =========================================================================

    integer function find_field(st, name)
        ! Index of `name` (exact match after trim), 0 if absent.
        type(field_state_t), intent(in) :: st
        character(len=*),    intent(in) :: name
        integer :: i

        find_field = 0
        do i = 1, st%n
            if (trim(st%f(i)%name) == trim(name)) then
                find_field = i
                return
            end if
        end do

    end function find_field

    subroutine check_handle(st, h, proc)
        type(field_state_t), intent(in) :: st
        integer,             intent(in) :: h
        character(len=*),    intent(in) :: proc

        if (h < 1 .or. h > st%n) then
            call fesm_error(proc, "invalid handle", "owner = '"//trim(st%owner)//"'"//new_line("a")// &
                "h = "//to_str(h)//", n = "//to_str(st%n)//new_line("a")// &
                "declared: "//declared_names(st))
        end if

    end subroutine check_handle

    subroutine check_access(st, h, proc)
        ! The state is locked and h is one of its handles.
        type(field_state_t), intent(in) :: st
        integer,             intent(in) :: h
        character(len=*),    intent(in) :: proc

        if (.not. st%locked) then
            call fesm_error(proc, "state is not locked", "owner = '"//trim(st%owner)//"'"// &
                new_line("a")//"h = "//to_str(h)//new_line("a")//"declared: "//declared_names(st))
        end if
        call check_handle(st, h, proc)

    end subroutine check_access

    subroutine check_shape(st, h, proc, rk, shp)
        ! The rank and shape of a put/get argument match the declaration. A
        ! rank-3 argument with a third extent of 1 is accepted for a rank-2
        ! field.
        type(field_state_t), intent(in) :: st
        integer,             intent(in) :: h, rk, shp(3)
        character(len=*),    intent(in) :: proc

        logical :: ok

        associate (f => st%f(h))
            ok = (rk == f%rank .or. (rk == 3 .and. f%rank == 2)) .and. &
                 shp(1) == f%nx .and. shp(2) == f%ny .and. shp(3) == f%nz
            if (.not. ok) then
                call fesm_error(proc, "shape mismatch", ctx(st, h)//new_line("a")// &
                    "declared rank "//to_str(f%rank)//", [nx,ny,nz] = ["// &
                    to_str([f%nx, f%ny, f%nz])//"]"//new_line("a")// &
                    "given    rank "//to_str(rk)//", shape      = ["//to_str(shp(1:max(1,rk)))//"]")
            end if
        end associate

    end subroutine check_shape

    subroutine bad_decl(st, name, what)
        type(field_state_t), intent(in) :: st
        character(len=*),    intent(in) :: name, what

        call fesm_error("fesm_fields_add", "invalid declaration", ctx_name(st, name)//new_line("a")//what)

    end subroutine bad_decl

    function ctx(st, h) result(s)
        ! "owner = '...'" and "field = '...' (h = ...)" on two lines.
        type(field_state_t), intent(in) :: st
        integer,             intent(in) :: h
        character(len=:), allocatable :: s

        s = "owner = '"//trim(st%owner)//"'"//new_line("a")// &
            "field = '"//trim(st%f(h)%name)//"' (h = "//to_str(h)//")"

    end function ctx

    function ctx_name(st, name) result(s)
        type(field_state_t), intent(in) :: st
        character(len=*),    intent(in) :: name
        character(len=:), allocatable :: s

        s = "owner = '"//trim(st%owner)//"'"//new_line("a")//"field = '"//trim(name)//"'"

    end function ctx_name

    function declared_names(st, sep) result(s)
        ! The declared names, joined by `sep` (default ", ").
        type(field_state_t), intent(in) :: st
        character(len=*), intent(in), optional :: sep
        character(len=:), allocatable :: s, sp_
        integer :: i

        sp_ = ", "
        if (present(sep)) sp_ = sep
        s = ""
        do i = 1, st%n
            if (i > 1) s = s//sp_
            s = s//trim(st%f(i)%name)
        end do
        if (st%n == 0 .and. .not. present(sep)) s = "(none)"

    end function declared_names

    function shape_str(f) result(s)
        type(field_t), intent(in) :: f
        character(len=:), allocatable :: s

        select case (f%rank)
            case (0)
                s = "scalar"
            case (2)
                s = to_str(f%nx)//"x"//to_str(f%ny)
            case default
                s = to_str(f%nx)//"x"//to_str(f%ny)//"x"//to_str(f%nz)//" ("//trim(f%zdim)//")"
        end select

    end function shape_str

    function tsem_name(tsem) result(s)
        integer, intent(in) :: tsem
        character(len=16) :: s

        select case (tsem)
            case (FIELD_INSTANT); s = "INSTANT"
            case (FIELD_MEAN);    s = "MEAN"
            case (FIELD_ACCUM);   s = "ACCUM"
            case default;         s = "tsem="//to_str(tsem)
        end select

    end function tsem_name

    function cons_name(cons) result(s)
        integer, intent(in) :: cons
        character(len=16) :: s

        select case (cons)
            case (FIELD_INTENSIVE); s = "intensive"
            case (FIELD_EXTENSIVE); s = "extensive"
            case default;           s = "cons="//to_str(cons)
        end select

    end function cons_name

    function t_str(t) result(s)
        ! A time at full precision (17 significant digits), for interval errors.
        real(dp), intent(in) :: t
        character(len=:), allocatable :: s
        character(len=32) :: buf

        write(buf, "(es25.17e3)") t
        s = trim(adjustl(buf))

    end function t_str

    function pad(s, n) result(out)
        ! `s` padded with blanks to at least n characters.
        character(len=*), intent(in) :: s
        integer,          intent(in) :: n
        character(len=max(n, len(s))) :: out

        out = s

    end function pad

end module fesm_fields
