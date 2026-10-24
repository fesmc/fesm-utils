program test_bind
    ! fesm_datasrc binding: the lookup order (par file over the library);
    ! an unbound required field stops, an unbound optional one is skipped;
    ! unknown key, missing group, undeclared name stop; a clim_db group is
    ! accepted unchanged; a MEAN field with datasrc_update(ds, t) stops.

    use precision
    use ncio
    use coordinates,  only : grid_class, grid_init
    use fesm_fields
    use fesm_datasrc

    implicit none

    character(len=*), parameter :: fpar = "test_ds_bind_par.nml", flib = "test_ds_bind_lib.nml", &
                                   fdef = "par/datasrc_defaults.nml"
    integer, parameter :: nx = 3, ny = 2
    type(grid_class)    :: grid
    type(field_state_t) :: st
    type(datasrc_t)     :: ds
    real(dp) :: x, v2(nx,ny)
    real(wp) :: sf(nx,ny)
    integer  :: h, stat, nfail

    nfail = 0

    call grid_init(grid, name="TEST-DS-BIND", mtype="cartesian", units="kilometers", &
                   x0=0.0_dp, dx=10.0_dp, nx=nx, y0=0.0_dp, dy=10.0_dp, ny=ny)

    sf = 1000.0_wp
    call nc_create("test_ds_bind.nc")
    call nc_write_dim("test_ds_bind.nc", "xc", x=0.0_dp, dx=10.0_dp, nx=nx)
    call nc_write_dim("test_ds_bind.nc", "yc", x=0.0_dp, dx=10.0_dp, nx=ny)
    call nc_write("test_ds_bind.nc", "sf", sf, dim1="xc", dim2="yc", units="kg m-2 yr-1")

    call write_nml()

    call fields_init(st, "bind_test:exports")
    call fields_add(st, "a",   "m",           FIELD_INSTANT, FIELD_INTENSIVE, h)
    call fields_add(st, "b",   "m",           FIELD_INSTANT, FIELD_INTENSIVE, h)
    call fields_add(st, "c",   "m",           FIELD_INSTANT, FIELD_INTENSIVE, h)
    call fields_add(st, "sf",  "kg m-2 yr-1", FIELD_INSTANT, FIELD_EXTENSIVE, h, grid="TEST-DS-BIND", dims=[nx,ny])
    call fields_add(st, "m",   "m",           FIELD_MEAN,    FIELD_INTENSIVE, h)
    call fields_lock(st)

    ! --- lookup order ---------------------------------------------------------------
    call datasrc_init(ds, fpar, "inst", st, [character(len=2) :: "a", "b", "c"], grid, fdef, &
                      required=[character(len=2) :: "a", "b"], var_defs=[character(len=64) :: "", flib], stat=stat)
    call check("init", stat .eq. 0)
    call datasrc_update(ds, 0.0_dp)
    call datasrc_get(ds, "a", x)
    call check("group in the par file and the library: the par file's", x .eq. 1.0_dp)
    call datasrc_get(ds, "b", x)
    call check("group only in the library", x .eq. 3.0_dp)
    call check("unbound optional field: not bound", .not. datasrc_has(ds, "c"))
    call datasrc_end(ds)

    ! --- a clim_db group, unchanged (static, kg m-2 yr-1 -> mm d-1 by scaling) -------
    call datasrc_init(ds, fpar, "db", st, ["sf"], grid, fdef, var_defs=[flib], stat=stat)
    call check("clim_db group accepted", stat .eq. 0)
    call datasrc_update(ds, 0.0_dp)
    call datasrc_get(ds, "sf", v2)
    call check("clim_db group: scaling, then mm d-1 -> kg m-2 yr-1", &
               maxval(abs(v2 - real(1000.0_wp*2.73785e-3_wp, dp)*365.25_dp)) .lt. 1e-2_dp)
    call datasrc_end(ds)

    ! --- errors (stat) ------------------------------------------------------------------
    call datasrc_init(ds, fpar, "inst", st, [character(len=2) :: "a", "b", "c"], grid, fdef, &
                      var_defs=[flib], stat=stat)
    call check("unbound required field: stat /= 0", stat .ne. 0)
    call datasrc_init(ds, fpar, "unknown_key", st, ["a"], grid, fdef, var_defs=[flib], stat=stat)
    call check("unknown key in the field group: stat /= 0", stat .ne. 0)
    call datasrc_init(ds, fpar, "missing_group", st, ["a"], grid, fdef, var_defs=[flib], stat=stat)
    call check("field group in no file: stat /= 0", stat .ne. 0)
    call datasrc_init(ds, fpar, "inst", st, ["zz"], grid, fdef, var_defs=[flib], stat=stat)
    call check("name not declared in the state: stat /= 0", stat .ne. 0)
    call datasrc_init(ds, fpar, "inst", st, ["a"], grid, fdef, required=["b"], var_defs=[flib], stat=stat)
    call check("required name not among the names: stat /= 0", stat .ne. 0)

    call datasrc_init(ds, fpar, "mean", st, ["m"], grid, fdef, var_defs=[flib], stat=stat)
    call check("init MEAN scalar const", stat .eq. 0)
    call datasrc_update(ds, 10.0_dp, stat=stat)
    call check("MEAN field with update(t): stat /= 0", stat .ne. 0)
    call datasrc_update(ds, 0.0_dp, 10.0_dp, stat=stat)
    call datasrc_fill(ds, st)
    call fields_finalize_interval(st, 10.0_dp)
    call fields_get(st, fields_handle(st, "m"), x)
    call check("MEAN const over [t0,t1]", stat .eq. 0 .and. x .eq. 7.0_dp)
    call datasrc_end(ds)

    if (nfail .gt. 0) then
        write(*,*) "test_bind: ", nfail, " FAILED"
        error stop 1
    end if
    write(*,*) "test_bind: all passed"

contains

    subroutine write_nml()
        integer :: u
        open(newunit=u, file=fpar, status="replace", action="write")
        write(u,"(a)") "&inst"
        write(u,"(a)") '    a = "g_both"'
        write(u,"(a)") '    b = "g_lib"'
        write(u,"(a)") '    c = ""'
        write(u,"(a)") "/"
        write(u,"(a)") "&db"
        write(u,"(a)") '    sf = "mar_sf"'
        write(u,"(a)") "/"
        write(u,"(a)") "&unknown_key"
        write(u,"(a)") '    a = "g_typo"'
        write(u,"(a)") "/"
        write(u,"(a)") "&missing_group"
        write(u,"(a)") '    a = "g_nowhere"'
        write(u,"(a)") "/"
        write(u,"(a)") "&mean"
        write(u,"(a)") '    m = "g_seven"'
        write(u,"(a)") "/"
        write(u,"(a)") "&g_both"
        write(u,"(a)") '    mode  = "const"'
        write(u,"(a)") '    value = 1.0'
        write(u,"(a)") "/"
        close(u)

        open(newunit=u, file=flib, status="replace", action="write")
        write(u,"(a)") "&g_both"
        write(u,"(a)") '    mode  = "const"'
        write(u,"(a)") '    value = 2.0'
        write(u,"(a)") "/"
        write(u,"(a)") "&g_lib"
        write(u,"(a)") '    mode  = "const"'
        write(u,"(a)") '    value = 3.0'
        write(u,"(a)") "/"
        write(u,"(a)") "&g_seven"
        write(u,"(a)") '    mode  = "const"'
        write(u,"(a)") '    value = 7.0'
        write(u,"(a)") "/"
        write(u,"(a)") "&g_typo"
        write(u,"(a)") '    mode   = "const"'
        write(u,"(a)") '    valeu  = 2.0'
        write(u,"(a)") "/"
        write(u,"(a)") "! as in yelmox clim_db.nml"
        write(u,"(a)") "&mar_sf"
        write(u,"(a)") '    filename = "test_ds_bind.nc"'
        write(u,"(a)") '    name     = "sf"'
        write(u,"(a)") '    units    = "kg m-2 yr-1" "mm d-1"'
        write(u,"(a)") '    scaling  = 2.73785e-3 0.0'
        write(u,"(a)") '    time     = 0 0 0 0 1'
        write(u,"(a)") "/"
        close(u)
    end subroutine write_nml

    subroutine check(label, ok)
        character(len=*), intent(in) :: label
        logical,          intent(in) :: ok
        write(*,"(a,a)") merge("  ok   ", "  FAIL ", ok), label
        if (.not. ok) nfail = nfail + 1
    end subroutine check

end program test_bind
