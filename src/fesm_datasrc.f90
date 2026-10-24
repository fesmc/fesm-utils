module fesm_datasrc
    ! Fields served from files into a field_state_t.
    !
    ! A datasrc_t serves a list of named fields of a caller's state (exports
    ! of a "data" component): each name is bound, in the caller's instance
    ! group, to a field group (file, variable, mode) found in the par file or
    ! a library file; the group is read with varslice (2-D, 3-D) or series
    ! (0-D), with the time axis and the units of the file; the field is
    ! remapped once at init onto the instance grid, converted to the units
    ! declared in the state (or init stops), evaluated at a time or over an
    ! interval, optionally as an anomaly against a reference, and put with
    ! its validity. It knows no field names and does no physics.
    !
    !     call datasrc_init(ds, par_file, "ice_bg", exports, names, grid, defaults_file, &
    !                       var_defs=var_defs, time_ref=t_ref, time_span=t_span)
    !     call datasrc_describe(ds)
    !     call datasrc_update(ds, t0, t1)     ! or datasrc_update(ds, t) for INSTANT only
    !     call datasrc_fill(ds, exports)
    !
    ! Field groups are read with the schema [datasrc_field] of defaults_file
    ! (par/datasrc_defaults.nml), so a group holds only what differs; every
    ! key of a group must be in the schema. Modes: const (value), static and
    ! clim (read at init), transient (read at each new update). Time methods
    ! (time_method = auto): INSTANT fields interp, MEAN/ACCUM interval_mean.
    ! Outside the axis: out_of_range = error (at init against time_span and
    ! at every update) or hold. Units: units(1) (or the attribute) is the
    ! file's, units(2) (or units(1)) that after scaling; units_conv
    ! (fesm_units) converts to the declared units (ACCUM: "<declared> yr-1",
    ! the state integrates the mean rate). Anomalies: anomaly_ref = slice
    ! (the field at ref_time of the same data) or group (a static group of
    ! the same quantity); converted before differencing.
    !
    ! Errors stop through fesm_error naming the instance, field, group and
    ! file; with `stat` present, datasrc_init and datasrc_update set stat = 1
    ! and return instead (for tests). Errors raised inside varslice and nml
    ! (a malformed file or group) always stop.

    use, intrinsic :: iso_fortran_env, only : int64
    use precision,   only : dp
    use constants,   only : mv
    use fesm_log,    only : LOG_INFO, log_msg, fesm_error, to_str
    use fesm_fields, only : field_state_t, FIELD_INSTANT, FIELD_ACCUM, FIELD_EXTENSIVE, &
                            fields_has, fields_handle, fields_put, units_canonical
    use fesm_units,  only : units_conv
    use nml,         only : nml_find_group, nml_validate, nml_read, nml_has_param, nml_replace
    use ncio,        only : nc_exists_attr, nc_read_attr
    use coordinates, only : grid_class
    use varslice,    only : varslice_class, varslice_init_nml, varslice_update, varslice_get, &
                            varslice_end, varslice_nsub, varslice_time_decode, get_matching_files
    use series,      only : series_class, series_load

    implicit none

    private

    character(len=*), parameter :: DEFAULTS_GROUP = "datasrc_field"
    integer, parameter :: LEN_FILE = 1024

    type :: datasrc_field_t
        ! One bound field (private)
        character(len=32)        :: name  = ""         ! export name (the caller's)
        character(len=128)       :: group = ""         ! field group
        character(len=LEN_FILE)  :: gfile = ""         ! file the group was found in
        character(len=LEN_FILE)  :: file  = ""         ! data file (first of the set)
        integer :: h = 0, rank = 0                     ! state handle; rank 0 | 2 | 3
        integer :: nx = 1, ny = 1, nz = 1
        integer :: tsem = FIELD_INSTANT
        logical :: masked = .false.
        character(len=16) :: mode = "", method = "", out_of_range = "", anomaly_ref = ""
        character(len=128) :: units_data = "", units_decl = "" ! after scaling -> declared (rate for ACCUM)
        real(dp) :: scale = 1.0_dp, offset = 0.0_dp   ! after scaling -> declared
        logical  :: identity = .true.                 ! no arithmetic
        real(dp) :: scaling(2) = [1.0_dp, 0.0_dp]     ! series only (varslice applies its own)
        real(dp) :: value = 0.0_dp                    ! const
        real(dp) :: trange(2) = 0.0_dp                ! coverage of the file [clock yr]
        integer  :: nsub = 1                          ! sub-annual steps of a transient axis
        character(len=32)  :: remap = "none"
        character(len=256) :: grid_src = "", grid_tgt = ""
        type(varslice_class)  :: vs                   ! rank 2, 3
        type(series_class)    :: ser                  ! rank 0
        real(dp), allocatable :: ref(:,:,:)           ! anomaly reference, declared units
        logical,  allocatable :: ref_ok(:,:,:)
        real(dp), allocatable :: buf(:,:,:)           ! current value (nx,ny,nz); rank 0: (1,1,1)
        logical,  allocatable :: valid(:,:)           ! masked fields only
    end type datasrc_field_t

    type, public :: datasrc_t
        character(len=128) :: group    = ""           ! instance group
        real(dp)           :: time_ref = 0.0_dp       ! calendar year of the caller's time zero
        real(dp)           :: t0 = 0.0_dp, t1 = 0.0_dp ! last update
        logical            :: updated  = .false.
        type(datasrc_field_t), allocatable :: f(:)    ! bound fields only
    end type datasrc_t

    interface datasrc_update
        module procedure datasrc_update_t, datasrc_update_iv
    end interface datasrc_update

    interface datasrc_get
        module procedure datasrc_get_0d, datasrc_get_2d, datasrc_get_3d
    end interface datasrc_get

    public :: datasrc_init, datasrc_update, datasrc_fill, datasrc_has
    public :: datasrc_describe, datasrc_get, datasrc_end

contains

    ! =========================================================================
    ! Init
    ! =========================================================================

    subroutine datasrc_init(ds, par_file, group, state, names, grid, defaults_file, &
                            required, var_defs, domain, subs, time_ref, time_span, &
                            rho_ice, sec_year, stat)
        ! Bind every name of `names` that the instance group `group` of
        ! par_file binds (key <name> = "<field group>"; "" or absent:
        ! unbound), read static, clim and const fields, and check units,
        ! shape, grid and coverage. Names in `required` (default: all) must
        ! be bound.

        type(datasrc_t),     intent(out) :: ds
        character(len=*),    intent(in)  :: par_file, group     ! the instance group
        type(field_state_t), intent(in)  :: state               ! declared (locked) exports
        character(len=*),    intent(in)  :: names(:)            ! fields this backend may serve
        type(grid_class),    intent(in)  :: grid                ! instance grid: {grid_name}, remap target
        character(len=*),    intent(in)  :: defaults_file       ! holds [datasrc_field]
        character(len=*),    intent(in), optional :: required(:)    ! must be bound (default: all)
        character(len=*),    intent(in), optional :: var_defs(:)    ! library files after par_file
        character(len=*),    intent(in), optional :: domain, subs(:,:)
        real(dp),            intent(in), optional :: time_ref       ! default 0: calendar years
        real(dp),            intent(in), optional :: time_span(2)   ! run interval: coverage checked
        real(dp),            intent(in), optional :: rho_ice        ! enables "m ice yr-1"
        real(dp),            intent(in), optional :: sec_year       ! default 365.25 d
        integer,             intent(out), optional :: stat          ! tests: errors return

        character(len=*), parameter :: proc = "datasrc_init"
        character(len=LEN_FILE), allocatable :: files(:)
        type(datasrc_field_t),   allocatable :: tmp(:)
        character(len=128) :: key
        logical :: req
        integer :: i, n

        if (present(stat)) stat = 0

        allocate(ds%f(0))                ! an empty source if init fails with stat
        ds%group    = group
        ds%time_ref = 0.0_dp
        if (present(time_ref)) ds%time_ref = time_ref
        ds%updated  = .false.

        ! Lookup order: the par file, then the library
        if (present(var_defs)) then
            allocate(files(1+size(var_defs)))
            files(2:) = var_defs
        else
            allocate(files(1))
        end if
        files(1) = par_file

        if (present(required)) then
            do i = 1, size(required)
                if (.not. any(names == required(i))) then
                    call fesm_error(proc, "a required field is not among the names served", &
                        "instance = ["//trim(group)//"]"//new_line("a")// &
                        "field    = "//trim(required(i))//new_line("a")// &
                        "names    = "//join(names), stat)
                    return
                end if
            end do
        end if

        allocate(tmp(size(names)))
        n = 0
        do i = 1, size(names)

            if (.not. fields_has(state, names(i))) then
                call fesm_error(proc, "field is not declared in the state", &
                    "instance = ["//trim(group)//"]"//new_line("a")// &
                    "field    = "//trim(names(i))//new_line("a")// &
                    "state    = "//trim(state%owner), stat)
                return
            end if

            req = .true.
            if (present(required)) req = any(required == names(i))

            key = ""
            if (nml_has_param(par_file, group, names(i))) call nml_read(par_file, group, names(i), key)
            if (len_trim(key) == 0) then
                if (req) then
                    call fesm_error(proc, "field `"//trim(names(i))//"` of ["//trim(group)// &
                        "] is not bound: set `"//trim(names(i))//" = ""<group>""`", &
                        "par file = "//trim(par_file), stat)
                    return
                end if
                cycle
            end if

            n = n + 1
            call bind_field(tmp(n), names(i), key, files, ds, state, grid, defaults_file, &
                            domain, subs, time_span, rho_ice, sec_year, stat)
            if (failed(stat)) return

        end do

        deallocate(ds%f)
        allocate(ds%f(n))
        do i = 1, n
            ds%f(i) = tmp(i)
        end do

    end subroutine datasrc_init


    recursive subroutine bind_field(fd, name, key, files, ds, state, grid, defaults_file, &
                                    domain, subs, time_span, rho_ice, sec_year, stat, force_static)
        ! Bind field `name` to field group `key`: find and validate the group,
        ! read its keys, open the reader, derive the unit conversion, check
        ! shape, grid and coverage, build the anomaly reference, and read
        ! const, static and clim fields. force_static: a reference group.

        type(datasrc_field_t), intent(inout) :: fd
        character(len=*),      intent(in)    :: name, key, files(:)
        type(datasrc_t),       intent(in)    :: ds
        type(field_state_t),   intent(in)    :: state
        type(grid_class),      intent(in)    :: grid
        character(len=*),      intent(in)    :: defaults_file
        character(len=*),      intent(in), optional :: domain, subs(:,:)
        real(dp),              intent(in), optional :: time_span(2), rho_ice, sec_year
        integer,               intent(inout), optional :: stat
        logical,               intent(in), optional :: force_static

        character(len=*), parameter :: proc = "datasrc_init"
        character(len=256) :: filename, vname, gsrc_key, ref_group, time_units
        character(len=128) :: units(2), u_file, u_attr, target
        character(len=32)  :: remap_key
        character(len=16)  :: mode_key, time_method
        real(dp) :: scaling(2), value, ref_time
        integer  :: h, nb, st
        logical  :: is_ref
        type(datasrc_field_t) :: rf
        character(len=:), allocatable :: ctx

        is_ref = .false.
        if (present(force_static)) is_ref = force_static

        h = fields_handle(state, name)
        fd%name   = name
        fd%group  = key
        fd%h      = h
        fd%rank   = state%f(h)%rank
        fd%nx     = state%f(h)%nx
        fd%ny     = state%f(h)%ny
        fd%nz     = state%f(h)%nz
        fd%tsem   = state%f(h)%tsem
        fd%masked = state%f(h)%masked
        fd%grid_tgt = grid%name

        ctx = "instance = ["//trim(ds%group)//"]"//new_line("a")// &
              "field    = "//trim(name)//new_line("a")// &
              "group    = "//trim(key)

        ! --- The group: lookup and keys ----------------------------------------
        fd%gfile = nml_find_group(files, key, what="["//trim(ds%group)//"] "//trim(name), stat=st)
        if (st /= 0) then
            call fesm_error(proc, "field group not found (see above)", ctx, stat)
            return
        end if
        ctx = ctx//new_line("a")//"in       = "//trim(fd%gfile)

        call nml_validate(fd%gfile, defaults_file, key, DEFAULTS_GROUP, halt_on_error=.false., nbad=nb)
        if (nb > 0) then
            call fesm_error(proc, "unknown key(s) in the field group (see above)", &
                ctx//new_line("a")//"schema   = ["//DEFAULTS_GROUP//"] of "//trim(defaults_file), stat)
            return
        end if

        mode_key = "auto"; time_method = "auto"; units = "auto"; ref_group = ""
        filename = ""; vname = ""; remap_key = "auto"; gsrc_key = ""; time_units = ""
        scaling = [1.0_dp, 0.0_dp]
        value = 0.0_dp; ref_time = 0.0_dp
        fd%out_of_range = "error"; fd%anomaly_ref = "none"
        call rd_str(key, "mode", mode_key)
        call rd_str(key, "time_method", time_method)
        call rd_str(key, "out_of_range", fd%out_of_range)
        call rd_str(key, "anomaly_ref", fd%anomaly_ref)
        call rd_str(key, "ref_group", ref_group)
        call rd_str(key, "filename", filename)
        call rd_str(key, "name", vname)
        call rd_str(key, "remap", remap_key)
        call rd_str(key, "grid_src", gsrc_key)
        call rd_str(key, "time_units", time_units)
        call nml_read(fd%gfile, key, "units",    units,    defaults_file=defaults_file, defaults_group=DEFAULTS_GROUP)
        call nml_read(fd%gfile, key, "scaling",  scaling,  defaults_file=defaults_file, defaults_group=DEFAULTS_GROUP)
        call nml_read(fd%gfile, key, "value",    value,    defaults_file=defaults_file, defaults_group=DEFAULTS_GROUP)
        call nml_read(fd%gfile, key, "ref_time", ref_time, defaults_file=defaults_file, defaults_group=DEFAULTS_GROUP)

        if (.not. one_of(mode_key, [character(len=16) :: "auto", "const", "static", "clim", "transient"])) then
            call fesm_error(proc, "mode must be auto, const, static, clim or transient", &
                ctx//new_line("a")//"mode     = "//trim(mode_key), stat)
            return
        end if
        if (.not. one_of(time_method, [character(len=16) :: "auto", "interp", "mean", "exact"])) then
            call fesm_error(proc, "time_method must be auto, interp, mean or exact", &
                ctx//new_line("a")//"time_method = "//trim(time_method), stat)
            return
        end if
        if (.not. one_of(fd%out_of_range, [character(len=16) :: "error", "hold"])) then
            call fesm_error(proc, "out_of_range must be error or hold", &
                ctx//new_line("a")//"out_of_range = "//trim(fd%out_of_range), stat)
            return
        end if
        if (.not. one_of(fd%anomaly_ref, [character(len=16) :: "none", "slice", "group"])) then
            call fesm_error(proc, "anomaly_ref must be none, slice or group", &
                ctx//new_line("a")//"anomaly_ref = "//trim(fd%anomaly_ref), stat)
            return
        end if
        if (is_ref) then
            mode_key = "static"
            fd%anomaly_ref = "none"
        end if

        ! --- Shape and grid of the declaration --------------------------------
        if (fd%rank >= 2) then
            if (trim(state%f(h)%grid) /= trim(grid%name)) then
                call fesm_error(proc, "the field is declared on another grid than the instance grid", &
                    ctx//new_line("a")//"declared on = "//trim(state%f(h)%grid)//new_line("a")// &
                    "instance grid = "//trim(grid%name), stat)
                return
            end if
            if (fd%nx /= grid%G%nx .or. fd%ny /= grid%G%ny) then
                call fesm_error(proc, "the declared shape is not that of the instance grid", &
                    ctx//new_line("a")//"declared [nx,ny] = ["//to_str([fd%nx, fd%ny])//"]"//new_line("a")// &
                    "grid "//trim(grid%name)//" = ["//to_str([grid%G%nx, grid%G%ny])//"]", stat)
                return
            end if
        end if

        allocate(fd%buf(fd%nx, fd%ny, fd%nz))
        fd%buf = 0.0_dp
        if (fd%masked) then
            allocate(fd%valid(fd%nx, fd%ny))
            fd%valid = .true.
        end if

        ! --- const --------------------------------------------------------------
        if (trim(mode_key) == "const") then
            if (trim(fd%anomaly_ref) /= "none") then
                call fesm_error(proc, "a const field has no anomaly reference", ctx, stat)
                return
            end if
            fd%mode   = "const"
            fd%method = "none"
            fd%value  = value
            fd%buf    = value
            fd%units_data = state%f(h)%units
            fd%units_decl = state%f(h)%units
            fd%trange = [-huge(1.0_dp), huge(1.0_dp)]
            return
        end if

        if (fd%rank == 1) then
            call fesm_error(proc, "rank 1 is not supported", ctx, stat)
            return
        end if

        ! --- Open the reader -----------------------------------------------------
        if (fd%rank == 0) then
            call open_series(stat)
        else
            call open_varslice()
        end if
        if (failed(stat)) return
        ctx = ctx//new_line("a")//"file     = "//trim(fd%file)

        ! --- Units -----------------------------------------------------------------
        if (trim(units(1)) /= "auto") then
            u_file = units(1)
            if (len_trim(u_attr) > 0 .and. units_canonical(u_attr) /= units_canonical(u_file)) then
                call fesm_error(proc, "units(1) of the group and the units attribute of the file differ", &
                    ctx//new_line("a")//"units(1)  = '"//trim(u_file)//"'"//new_line("a")// &
                    "attribute = '"//trim(u_attr)//"'", stat)
                return
            end if
        else
            u_file = u_attr
            if (len_trim(u_file) == 0) then
                call fesm_error(proc, "no units: the variable has no units attribute, set units(1)", ctx, stat)
                return
            end if
        end if
        fd%units_data = u_file
        if (trim(units(2)) /= "auto") fd%units_data = units(2)

        target = state%f(h)%units
        if (fd%tsem == FIELD_ACCUM) target = trim(target)//" yr-1"
        fd%units_decl = target
        call units_conv(fd%units_data, target, fd%scale, fd%offset, rho_ice=rho_ice, &
                        sec_year=sec_year, stat=st)
        if (st /= 0) then
            call fesm_error(proc, "no conversion from '"//trim(fd%units_data)//"' to '"//trim(target)// &
                "' for `"//trim(name)//"` of ["//trim(ds%group)//"]", ctx, stat)
            return
        end if
        fd%identity = (units_canonical(fd%units_data) == units_canonical(target))

        ! --- Time method -------------------------------------------------------------
        select case(trim(fd%mode))
            case("static")
                fd%method = "none"
            case("clim")
                if (trim(time_method) == "mean") then
                    fd%method = "mean"
                else if (trim(time_method) == "auto") then
                    fd%method = "cycle"
                else
                    call fesm_error(proc, "a climatology is served as its cycle (time_method = auto) "// &
                        "or its annual mean (time_method = mean)", &
                        ctx//new_line("a")//"time_method = "//trim(time_method), stat)
                    return
                end if
            case("transient")
                select case(trim(time_method))
                    case("auto")
                        fd%method = "interp"
                        if (fd%tsem /= FIELD_INSTANT) fd%method = "interval_mean"
                    case("mean")
                        fd%method = "interval_mean"
                    case default
                        fd%method = time_method
                end select
        end select

        call check_shape(stat)
        if (failed(stat)) return

        ! --- Coverage -------------------------------------------------------------------
        if (trim(fd%mode) == "transient" .and. trim(fd%out_of_range) == "error" .and. present(time_span)) then
            if (time_span(1) < fd%trange(1) .or. time_span(2) > fd%trange(2)) then
                call fesm_error(proc, "the file does not cover the run", &
                    ctx//new_line("a")//"run      = ["//to_str(time_span(1))//", "//to_str(time_span(2))//"]"// &
                    new_line("a")//"axis     = ["//to_str(fd%trange(1))//", "//to_str(fd%trange(2))//"]"// &
                    new_line("a")//"(set out_of_range = ""hold"" to hold the first or last sample)", stat)
                return
            end if
        end if

        ! --- Anomaly reference ------------------------------------------------------------
        select case(trim(fd%anomaly_ref))
            case("slice")
                if (trim(fd%mode) /= "transient") then
                    call fesm_error(proc, "anomaly_ref = ""slice"" needs a transient field", ctx, stat)
                    return
                end if
                if (ref_time < fd%trange(1) .or. ref_time > fd%trange(2)) then
                    call fesm_error(proc, "ref_time is outside the time axis of the file", &
                        ctx//new_line("a")//"ref_time = "//to_str(ref_time)//new_line("a")// &
                        "axis     = ["//to_str(fd%trange(1))//", "//to_str(fd%trange(2))//"]", stat)
                    return
                end if
                allocate(fd%ref(fd%nx,fd%ny,fd%nz), fd%ref_ok(fd%nx,fd%ny,fd%nz))
                call field_value(fd, ref_time, ref_time, "interp", fd%ref, fd%ref_ok)
            case("group")
                if (len_trim(ref_group) == 0) then
                    call fesm_error(proc, "anomaly_ref = ""group"" needs ref_group", ctx, stat)
                    return
                end if
                call bind_field(rf, name, ref_group, files, ds, state, grid, defaults_file, &
                                domain, subs, time_span, rho_ice, sec_year, stat, force_static=.true.)
                if (failed(stat)) return
                allocate(fd%ref(fd%nx,fd%ny,fd%nz), fd%ref_ok(fd%nx,fd%ny,fd%nz))
                call field_value(rf, 0.0_dp, 0.0_dp, "none", fd%ref, fd%ref_ok)
                call varslice_end(rf%vs)
        end select

        ! --- static and clim: read once ---------------------------------------------------
        if (trim(fd%mode) /= "transient") then
            call evaluate(fd, ds%group, 0.0_dp, 0.0_dp, stat)
            if (failed(stat)) return
        end if

    contains

        subroutine rd_str(grp, k, v)
            character(len=*), intent(in)    :: grp, k
            character(len=*), intent(inout) :: v
            call nml_read(fd%gfile, grp, k, v, defaults_file=defaults_file, defaults_group=DEFAULTS_GROUP)
        end subroutine rd_str

        subroutine open_varslice()
            ! varslice with the schema, the resolved remap and the instance time_ref
            character(len=LEN_FILE), allocatable :: fl(:)
            character(len=256) :: gsrc, fpath
            character(len=16)  :: vmode

            fpath = subst(filename)
            call get_matching_files(fl, trim(fpath))
            fd%file = fl(1)

            ! remap = auto: from a known grid of the file other than the instance
            ! grid, conservative for extensive fields
            gsrc = subst(gsrc_key)
            if (len_trim(gsrc) == 0) then
                if (nc_exists_attr(fd%file, "grid_name")) call nc_read_attr(fd%file, "grid_name", gsrc)
            end if
            fd%grid_src = gsrc
            fd%remap = remap_key
            if (trim(remap_key) == "auto") then
                fd%remap = "none"
                if (len_trim(gsrc) > 0 .and. trim(gsrc) /= trim(grid%name)) then
                    fd%remap = "bilin"
                    if (state%f(h)%cons == FIELD_EXTENSIVE) fd%remap = "con"
                end if
            end if

            vmode = mode_key
            call varslice_init_nml(fd%vs, fd%gfile, key, domain=domain, grid_name=grid%name, &
                                   verbose=.false., subs=subs, grid=grid, &
                                   defaults_file=defaults_file, defaults_group=DEFAULTS_GROUP, &
                                   mode=vmode, remap=fd%remap, grid_src=gsrc, time_ref=ds%time_ref)
            fd%mode = fd%vs%par%mode
            u_attr  = fd%vs%par%units_file

            if (trim(fd%mode) == "transient") then
                fd%nsub = varslice_nsub(fd%vs)
                if (allocated(fd%vs%time_bnds)) then
                    fd%trange = [fd%vs%time_bnds(1,1), fd%vs%time_bnds(2,size(fd%vs%time))]
                else
                    fd%trange = [fd%vs%time(1), fd%vs%time(size(fd%vs%time))]
                end if
                if (fd%nsub > 1 .and. fd%rank == 3 .and. fd%nz == fd%nsub) then
                    ! Served as the cycle of the year of t: whole years
                    fd%trange = [real(floor(fd%vs%time(1)),dp), real(floor(fd%vs%time(size(fd%vs%time)))+1,dp)]
                end if
            end if
        end subroutine open_varslice

        subroutine open_series(stat)
            ! series (ascii or netCDF) with its axis decoded (time_units, else
            ! the attributes of the time variable) to clock years
            integer, intent(inout), optional :: stat
            character(len=256) :: fpath, tu, cal
            real(dp), allocatable :: y(:)
            integer :: k, st2
            logical :: is_nc

            if (trim(mode_key) /= "auto" .and. trim(mode_key) /= "transient") then
                call fesm_error(proc, "a scalar field from a file is a time series (mode = transient; "// &
                    "a constant: mode = ""const"")", ctx//new_line("a")//"mode = "//trim(mode_key), stat)
                return
            end if
            if (trim(fd%anomaly_ref) == "group") then
                call fesm_error(proc, "anomaly_ref = ""group"" is not available for a scalar field", ctx, stat)
                return
            end if

            fpath = subst(filename)
            fd%file = fpath
            k = index(fpath, ".", back=.true.)
            is_nc = .false.
            if (k > 0) is_nc = one_of(fpath(k+1:), [character(len=16) :: "nc", "nc4", "cdf", "netcdf"])

            u_attr = ""
            cal = ""
            tu  = time_units
            if (is_nc) then
                call series_load(fd%ser, trim(fpath), varname=trim(vname), time_name="time")
                if (nc_exists_attr(fpath, trim(vname), "units")) call nc_read_attr(fpath, trim(vname), "units", u_attr)
                if (len_trim(tu) == 0) then
                    if (nc_exists_attr(fpath, "time", "units")) call nc_read_attr(fpath, "time", "units", tu)
                end if
                if (nc_exists_attr(fpath, "time", "calendar")) call nc_read_attr(fpath, "time", "calendar", cal)
            else
                call series_load(fd%ser, trim(fpath), fmt="ascii")
                if (trim(units(1)) == "auto") then
                    call fesm_error(proc, "an ascii series has no units attribute: set units(1)", ctx, stat)
                    return
                end if
            end if
            if (len_trim(tu) == 0) then
                call fesm_error(proc, "the time axis of the series has no units: set time_units", ctx, stat)
                return
            end if
            if (fd%ser%nc /= 1) then
                call fesm_error(proc, "a scalar field needs a series with one channel", &
                    ctx//new_line("a")//"channels = "//to_str(fd%ser%nc), stat)
                return
            end if

            if (trim(tu) /= "clock") then
                call varslice_time_decode(fd%ser%time, tu, cal, y, stat=st2)
                if (st2 /= 0) then
                    call fesm_error(proc, "the time axis of the series cannot be decoded (see above)", ctx, stat)
                    return
                end if
                fd%ser%time = y - ds%time_ref
            end if
            do k = 2, fd%ser%nt
                if (.not. (fd%ser%time(k) > fd%ser%time(k-1))) then
                    call fesm_error(proc, "the time axis of the series is not strictly increasing", &
                        ctx//new_line("a")//"index = "//to_str(k), stat)
                    return
                end if
            end do

            fd%mode    = "transient"
            fd%scaling = scaling
            fd%trange  = [fd%ser%time(1), fd%ser%time(fd%ser%nt)]

        end subroutine open_series

        subroutine check_shape(stat)
            ! The data as served against the declared shape
            integer, intent(inout), optional :: stat
            integer :: nd, nsp, nexp
            character(len=:), allocatable :: what

            if (fd%rank == 0) return

            nd  = fd%vs%par%ndim
            nsp = nd
            if (trim(fd%mode) /= "static") nsp = nd - 1      ! the last dimension is time or the cycle
            what = ""
            nexp = 1
            select case(trim(fd%mode))
                case("static")
                    if (nsp == 3) nexp = fd%vs%dim(3)
                    if (nsp /= fd%rank) what = "a static field of "//to_str(nsp)//" dimensions"
                case("clim")
                    if (nsp == 2) then
                        if (trim(fd%method) == "mean") then
                            if (fd%rank /= 2) what = "the annual mean of a 2-D climatology is rank 2"
                        else
                            nexp = fd%vs%dim(3)
                            if (fd%rank /= 3) what = "the cycle of a 2-D climatology is rank 3 (nz = steps)"
                        end if
                    else if (nsp == 3) then
                        nexp = fd%vs%dim(3)
                        if (trim(fd%method) /= "mean") then
                            what = "a 3-D climatology is served only as its annual mean (time_method = ""mean"")"
                        else if (fd%rank /= 3) then
                            what = "the annual mean of a 3-D climatology is rank 3"
                        end if
                    else
                        what = "a climatology of "//to_str(nsp)//" dimensions"
                    end if
                case("transient")
                    if (nsp == 2) then
                        if (fd%rank == 3) then
                            nexp = fd%nsub
                            if (fd%nsub <= 1) what = "rank 3 from a 2-D transient field needs a sub-annual axis"
                        end if
                    else if (nsp == 3) then
                        nexp = fd%vs%dim(3)
                        if (fd%rank /= 3) what = "a 3-D transient field is rank 3"
                    else
                        what = "a transient field of "//to_str(nsp)//" dimensions"
                    end if
            end select
            if (len(what) == 0 .and. fd%rank == 3 .and. fd%nz /= nexp) &
                what = "nz = "//to_str(fd%nz)//", the data have "//to_str(nexp)
            if (len(what) > 0) then
                call fesm_error(proc, "shape mismatch", ctx//new_line("a")// &
                    "declared rank "//to_str(fd%rank)//" (nz = "//to_str(fd%nz)//"), mode "// &
                    trim(fd%mode)//new_line("a")//what, stat)
            end if
        end subroutine check_shape

        function subst(path) result(p)
            character(len=*), intent(in) :: path
            character(len=256) :: p
            integer :: k
            p = path
            if (present(domain)) call nml_replace(p, "{domain}", trim(domain))
            call nml_replace(p, "{grid_name}", trim(grid%name))
            if (present(subs)) then
                do k = 1, size(subs,1)
                    call nml_replace(p, "{"//trim(subs(k,1))//"}", trim(subs(k,2)))
                end do
            end if
        end function subst

    end subroutine bind_field

    ! =========================================================================
    ! Update and fill
    ! =========================================================================

    subroutine datasrc_update_t(ds, t, stat)
        ! INSTANT fields at t; stops if a MEAN or ACCUM field is bound.
        ! A no-op for the time of the last update.

        type(datasrc_t), intent(inout) :: ds
        real(dp),        intent(in)    :: t
        integer,         intent(out), optional :: stat

        integer :: i

        if (present(stat)) stat = 0

        do i = 1, size(ds%f)
            if (ds%f(i)%tsem /= FIELD_INSTANT) then
                call fesm_error("datasrc_update", "a MEAN or ACCUM field needs datasrc_update(ds, t0, t1)", &
                    "instance = ["//trim(ds%group)//"]"//new_line("a")// &
                    "field    = "//trim(ds%f(i)%name), stat)
                return
            end if
        end do

        call update_all(ds, t, t, stat)

    end subroutine datasrc_update_t


    subroutine datasrc_update_iv(ds, t0, t1, stat)
        ! INSTANT fields at t1; MEAN and ACCUM fields: the mean (rate) over
        ! [t0,t1], t0 < t1. A no-op for the interval of the last update.

        type(datasrc_t), intent(inout) :: ds
        real(dp),        intent(in)    :: t0, t1
        integer,         intent(out), optional :: stat

        integer :: i

        if (present(stat)) stat = 0

        if (t0 > t1) then
            call fesm_error("datasrc_update", "invalid interval (t0 > t1)", &
                "instance = ["//trim(ds%group)//"]"//new_line("a")// &
                "t0 = "//to_str(t0)//", t1 = "//to_str(t1), stat)
            return
        end if
        do i = 1, size(ds%f)
            if (ds%f(i)%tsem /= FIELD_INSTANT .and. .not. (t0 < t1)) then
                call fesm_error("datasrc_update", "a MEAN or ACCUM field needs t0 < t1", &
                    "instance = ["//trim(ds%group)//"]"//new_line("a")// &
                    "field    = "//trim(ds%f(i)%name)//new_line("a")// &
                    "t0 = "//to_str(t0)//", t1 = "//to_str(t1), stat)
                return
            end if
        end do

        call update_all(ds, t0, t1, stat)

    end subroutine datasrc_update_iv


    subroutine update_all(ds, t0, t1, stat)

        type(datasrc_t), intent(inout) :: ds
        real(dp),        intent(in)    :: t0, t1
        integer,         intent(inout), optional :: stat

        integer :: i

        if (ds%updated .and. same(t0, ds%t0) .and. same(t1, ds%t1)) return

        do i = 1, size(ds%f)
            if (trim(ds%f(i)%mode) /= "transient") cycle
            call evaluate(ds%f(i), ds%group, t0, t1, stat)
            if (failed(stat)) return
        end do

        ds%t0 = t0
        ds%t1 = t1
        ds%updated = .true.

    end subroutine update_all


    subroutine evaluate(fd, inst, t0, t1, stat)
        ! buf (declared units, minus the reference) and valid of field fd for
        ! the interval [t0,t1] (INSTANT: at t1); coverage and missing values
        ! checked.

        type(datasrc_field_t), intent(inout) :: fd
        character(len=*),      intent(in)    :: inst
        real(dp),              intent(in)    :: t0, t1
        integer,               intent(inout), optional :: stat

        logical, allocatable :: ok(:,:,:)
        real(dp) :: ta, tb
        character(len=16) :: method
        character(len=:), allocatable :: ctx
        integer :: k

        ctx = "instance = ["//trim(inst)//"]"//new_line("a")// &
              "field    = "//trim(fd%name)//new_line("a")// &
              "group    = "//trim(fd%group)//new_line("a")// &
              "file     = "//trim(fd%file)

        method = fd%method
        ta = t0
        tb = t1
        if (trim(fd%mode) == "transient") then
            ! The time(s) at which the data are needed
            if (fd%tsem == FIELD_INSTANT) then
                ta = t1
            else if (trim(method) /= "interval_mean") then
                ta = 0.5_dp*(t0 + t1)
                tb = ta
            end if
            if (trim(method) /= "interval_mean") ta = tb
            if (trim(fd%out_of_range) == "error" .and. (ta < fd%trange(1) .or. tb > fd%trange(2))) then
                call fesm_error("datasrc_update", "time outside the axis of the file", &
                    ctx//new_line("a")//"time     = ["//to_str(ta)//", "//to_str(tb)//"]"//new_line("a")// &
                    "axis     = ["//to_str(fd%trange(1))//", "//to_str(fd%trange(2))//"]", stat)
                return
            end if
        end if

        allocate(ok(fd%nx, fd%ny, fd%nz))
        call field_value(fd, ta, tb, method, fd%buf, ok, stat)
        if (failed(stat)) return

        if (allocated(fd%ref)) then
            ok = ok .and. fd%ref_ok
            where (ok)
                fd%buf = fd%buf - fd%ref
            end where
        end if

        if (.not. all(ok)) then
            if (fd%masked) then
                do k = 1, fd%nz
                    where (.not. ok(:,:,k)) fd%buf(:,:,k) = 0.0_dp
                end do
                fd%valid = all(ok, dim=3)
            else
                call fesm_error("datasrc_update", "missing values in a field without a mask", &
                    ctx//new_line("a")//"time     = ["//to_str(ta)//", "//to_str(tb)//"]"//new_line("a")// &
                    "missing  = "//to_str(count(.not. ok))//" of "//to_str(size(ok))//new_line("a")// &
                    "(set fill = ""nearest"" or ""const"" in the group, or declare the field masked)", stat)
                return
            end if
        else if (fd%masked) then
            fd%valid = .true.
        end if

    end subroutine evaluate


    subroutine field_value(fd, ta, tb, method, out, ok, stat)
        ! The value of the data over [ta,tb] (a time: ta = tb) in declared
        ! units, without reference; ok = not missing.

        type(datasrc_field_t), intent(inout) :: fd
        real(dp),              intent(in)    :: ta, tb
        character(len=*),      intent(in)    :: method
        real(dp),              intent(out)   :: out(:,:,:)
        logical,               intent(out)   :: ok(:,:,:)
        integer,               intent(inout), optional :: stat

        real(dp) :: v
        integer  :: rep

        ok = .true.

        select case(trim(fd%mode))

            case("const")
                out = fd%value
                return

            case("static")
                call varslice_update(fd%vs)

            case("clim")
                if (trim(method) == "mean") then
                    call varslice_update(fd%vs, [0.0_dp, 0.0_dp], method="range_mean", rep=1)
                else
                    call varslice_update(fd%vs)
                end if

            case("transient")
                if (fd%rank == 0) then
                    select case(trim(method))
                        case("interval_mean")
                            v = ser_mean(fd, ta, tb)
                        case("exact")
                            v = ser_exact(fd, ta, stat)
                            if (failed(stat)) return
                        case default
                            v = ser_interp(fd, ta)
                    end select
                    out = (v*fd%scaling(1) + fd%scaling(2))
                    if (.not. fd%identity) out = out*fd%scale + fd%offset
                    return
                end if

                rep = 1
                if (fd%rank == 3 .and. fd%nsub > 1 .and. fd%nz == fd%nsub .and. fd%vs%par%ndim == 3) rep = fd%nsub
                select case(trim(method))
                    case("interval_mean")
                        call varslice_update(fd%vs, [ta, tb], method="interval_mean", rep=rep)
                    case("exact")
                        call varslice_update(fd%vs, [ta], method="exact")
                        if (fd%vs%idx(1) <= 0) then
                            call fesm_error("datasrc_update", "no sample of the file at this time "// &
                                "(time_method = ""exact"")", &
                                "field = "//trim(fd%name)//new_line("a")//"group = "//trim(fd%group)// &
                                new_line("a")//"time  = "//to_str(ta), stat)
                            return
                        end if
                    case default
                        ! interp, holding the first or last sample outside the axis
                        if (rep > 1) then
                            call varslice_update(fd%vs, [real(floor(ta), dp)], method="extrap", rep=rep)
                        else
                            call varslice_update(fd%vs, [ta], method="extrap")
                        end if
                end select

        end select

        ! dp copy at the declared rank
        if (fd%rank == 2) then
            call varslice_get(fd%vs, out(:,:,1))
        else
            call varslice_get(fd%vs, out)
        end if
        ok = .not. is_mv(out)
        if (.not. fd%identity) then
            where (ok) out = out*fd%scale + fd%offset
        end if

    end subroutine field_value


    subroutine datasrc_fill(ds, state)
        ! Put every bound field into state with the interval of the last
        ! update (INSTANT: [t1,t1]), static ones included.

        type(datasrc_t),     intent(in)    :: ds
        type(field_state_t), intent(inout) :: state

        integer  :: i
        real(dp) :: t0

        if (.not. ds%updated) then
            call fesm_error("datasrc_fill", "fill before the first update", &
                "instance = ["//trim(ds%group)//"]")
        end if

        do i = 1, size(ds%f)
            associate (fd => ds%f(i))
                t0 = ds%t0
                if (fd%tsem == FIELD_INSTANT) t0 = ds%t1
                select case(fd%rank)
                    case(0)
                        call fields_put(state, fd%h, fd%buf(1,1,1), t1=ds%t1, t0=t0)
                    case(2)
                        if (fd%masked) then
                            call fields_put(state, fd%h, fd%buf(:,:,1), t1=ds%t1, t0=t0, valid=fd%valid)
                        else
                            call fields_put(state, fd%h, fd%buf(:,:,1), t1=ds%t1, t0=t0)
                        end if
                    case(3)
                        if (fd%masked) then
                            call fields_put(state, fd%h, fd%buf, t1=ds%t1, t0=t0, valid=fd%valid)
                        else
                            call fields_put(state, fd%h, fd%buf, t1=ds%t1, t0=t0)
                        end if
                end select
            end associate
        end do

    end subroutine datasrc_fill

    ! =========================================================================
    ! Queries
    ! =========================================================================

    logical function datasrc_has(ds, name)
        ! Is `name` bound?

        type(datasrc_t),  intent(in) :: ds
        character(len=*), intent(in) :: name

        datasrc_has = (find(ds, name) > 0)

    end function datasrc_has


    subroutine datasrc_describe(ds, level)
        ! One line per bound field: group, file, mode, method, units ->
        ! declared (scale, offset), coverage, grid_src -> grid. Through
        ! log_msg at `level` (default LOG_INFO).

        type(datasrc_t), intent(in) :: ds
        integer, intent(in), optional :: level

        integer :: lvl, i
        character(len=:), allocatable :: line, conv

        lvl = LOG_INFO
        if (present(level)) lvl = level

        call log_msg(lvl, "fesm_datasrc", "["//trim(ds%group)//"] "//to_str(size(ds%f))//" bound fields")
        do i = 1, size(ds%f)
            associate (fd => ds%f(i))
                if (fd%identity) then
                    conv = "identity"
                else
                    conv = "x "//to_str(fd%scale)//" + "//to_str(fd%offset)
                end if
                line = "["//trim(ds%group)//"] "//trim(fd%name)//" <- &"//trim(fd%group)// &
                       " ("//trim(fd%gfile)//"): "//trim(fd%mode)
                if (trim(fd%mode) == "const") then
                    line = line//" "//to_str(fd%value)//" ["//trim(fd%units_decl)//"]"
                else
                    if (trim(fd%mode) /= "static") line = line//" ("//trim(fd%method)//")"
                    line = line//", "//trim(fd%file)//", '"//trim(fd%units_data)// &
                           "' -> '"//trim(fd%units_decl)//"' ("//conv//")"
                    if (trim(fd%mode) == "transient") line = line//", coverage ["//to_str(fd%trange(1))// &
                        ", "//to_str(fd%trange(2))//"] ("//trim(fd%out_of_range)//")"
                    if (fd%rank >= 2) then
                        if (trim(fd%remap) == "none") then
                            line = line//", on "//trim(fd%grid_tgt)
                        else
                            line = line//", "//trim(fd%grid_src)//" -> "//trim(fd%grid_tgt)// &
                                   " ("//trim(fd%remap)//")"
                        end if
                    end if
                    if (trim(fd%anomaly_ref) /= "none") line = line//", anomaly ("//trim(fd%anomaly_ref)//")"
                end if
                call log_msg(lvl, "fesm_datasrc", line)
            end associate
        end do

    end subroutine datasrc_describe


    subroutine datasrc_get_0d(ds, name, var)
        ! dp copy of the current value (tests, backend-internal use)
        type(datasrc_t),  intent(in)  :: ds
        character(len=*), intent(in)  :: name
        real(dp),         intent(out) :: var
        integer :: i
        i = get_check(ds, name, 0, [1,1,1])
        var = ds%f(i)%buf(1,1,1)
    end subroutine datasrc_get_0d

    subroutine datasrc_get_2d(ds, name, var)
        type(datasrc_t),  intent(in)  :: ds
        character(len=*), intent(in)  :: name
        real(dp),         intent(out) :: var(:,:)
        integer :: i
        i = get_check(ds, name, 2, [size(var,1), size(var,2), 1])
        var = ds%f(i)%buf(:,:,1)
    end subroutine datasrc_get_2d

    subroutine datasrc_get_3d(ds, name, var)
        type(datasrc_t),  intent(in)  :: ds
        character(len=*), intent(in)  :: name
        real(dp),         intent(out) :: var(:,:,:)
        integer :: i
        i = get_check(ds, name, 3, shape(var))
        var = ds%f(i)%buf
    end subroutine datasrc_get_3d


    subroutine datasrc_end(ds)

        type(datasrc_t), intent(inout) :: ds

        integer :: i

        if (allocated(ds%f)) then
            do i = 1, size(ds%f)
                call varslice_end(ds%f(i)%vs)
            end do
            deallocate(ds%f)
        end if
        ds%updated = .false.

    end subroutine datasrc_end

    ! =========================================================================
    ! Helpers
    ! =========================================================================

    real(dp) function ser_interp(fd, t) result(v)
        ! Linear interpolation of the series at t in dp, held outside
        type(datasrc_field_t), intent(in) :: fd
        real(dp),              intent(in) :: t
        integer :: n, k
        real(dp) :: lam
        n = fd%ser%nt
        if (t <= fd%ser%time(1)) then
            v = real(fd%ser%var(1,1), dp)
        else if (t >= fd%ser%time(n)) then
            v = real(fd%ser%var(n,1), dp)
        else
            k = 1
            do while (fd%ser%time(k+1) <= t)
                k = k + 1
            end do
            lam = (t - fd%ser%time(k))/(fd%ser%time(k+1) - fd%ser%time(k))
            v = (1.0_dp - lam)*real(fd%ser%var(k,1), dp) + lam*real(fd%ser%var(k+1,1), dp)
        end if
    end function ser_interp

    real(dp) function ser_mean(fd, ta, tb) result(v)
        ! Mean over [ta,tb] of the linear interpolant (held outside), exact:
        ! trapezoids between the interval ends and the samples inside it
        type(datasrc_field_t), intent(in) :: fd
        real(dp),              intent(in) :: ta, tb
        real(dp) :: x0, x1
        integer  :: k
        if (.not. (tb > ta)) then
            v = ser_interp(fd, ta)
            return
        end if
        v  = 0.0_dp
        x0 = ta
        do k = 1, fd%ser%nt + 1
            if (k <= fd%ser%nt) then
                if (fd%ser%time(k) <= x0) cycle
                x1 = min(fd%ser%time(k), tb)
            else
                x1 = tb
            end if
            v  = v + 0.5_dp*(x1 - x0)*(ser_interp(fd, x0) + ser_interp(fd, x1))
            x0 = x1
            if (x0 >= tb) exit
        end do
        v = v/(tb - ta)
    end function ser_mean

    real(dp) function ser_exact(fd, t, stat) result(v)
        type(datasrc_field_t), intent(in) :: fd
        real(dp),              intent(in) :: t
        integer,               intent(inout), optional :: stat
        integer :: k
        v = 0.0_dp
        do k = 1, fd%ser%nt
            if (abs(fd%ser%time(k) - t) <= 1e-9_dp*max(1.0_dp, abs(t))) then
                v = real(fd%ser%var(k,1), dp)
                return
            end if
        end do
        call fesm_error("datasrc_update", "no sample of the series at this time (time_method = ""exact"")", &
            "field = "//trim(fd%name)//new_line("a")//"group = "//trim(fd%group)//new_line("a")// &
            "time  = "//to_str(t), stat)
    end function ser_exact

    integer function get_check(ds, name, rk, shp) result(i)
        type(datasrc_t),  intent(in) :: ds
        character(len=*), intent(in) :: name
        integer,          intent(in) :: rk, shp(3)
        integer :: r
        i = find(ds, name)
        if (i == 0) then
            call fesm_error("datasrc_get", "field is not bound", &
                "instance = ["//trim(ds%group)//"]"//new_line("a")//"field    = "//trim(name))
        end if
        r = ds%f(i)%rank
        if (.not. (r == rk .or. (r == 2 .and. rk == 3 .and. shp(3) == 1)) .or. &
            shp(1) /= ds%f(i)%nx .or. shp(2) /= ds%f(i)%ny .or. (rk == 3 .and. shp(3) /= ds%f(i)%nz)) then
            call fesm_error("datasrc_get", "shape mismatch", &
                "instance = ["//trim(ds%group)//"]"//new_line("a")//"field    = "//trim(name)// &
                new_line("a")//"declared = rank "//to_str(r)//", ["// &
                to_str([ds%f(i)%nx, ds%f(i)%ny, ds%f(i)%nz])//"]"//new_line("a")// &
                "given    = rank "//to_str(rk)//", ["//to_str(shp)//"]")
        end if
    end function get_check

    integer function find(ds, name) result(i)
        type(datasrc_t),  intent(in) :: ds
        character(len=*), intent(in) :: name
        integer :: k
        i = 0
        if (.not. allocated(ds%f)) return
        do k = 1, size(ds%f)
            if (trim(ds%f(k)%name) == trim(name)) then
                i = k
                return
            end if
        end do
    end function find

    elemental logical function is_mv(x)
        ! The missing value of varslice (read as sp, widened exactly)
        real(dp), intent(in) :: x
        is_mv = abs(x - real(mv, dp)) < 1.0e-6_dp*abs(real(mv, dp))
    end function is_mv

    logical function same(a, b)
        ! Equal to the bit (times from the clock are n*tick)
        real(dp), intent(in) :: a, b
        same = (transfer(a, 0_int64) == transfer(b, 0_int64))
    end function same

    logical function failed(stat)
        integer, intent(in), optional :: stat
        failed = .false.
        if (present(stat)) failed = (stat /= 0)
    end function failed

    logical function one_of(s, list)
        character(len=*), intent(in) :: s, list(:)
        integer :: k
        one_of = .false.
        do k = 1, size(list)
            if (trim(s) == trim(list(k))) one_of = .true.
        end do
    end function one_of

    function join(list) result(s)
        character(len=*), intent(in) :: list(:)
        character(len=:), allocatable :: s
        integer :: k
        s = ""
        do k = 1, size(list)
            if (k > 1) s = s//", "
            s = s//trim(list(k))
        end do
    end function join

end module fesm_datasrc
