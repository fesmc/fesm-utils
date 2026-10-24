program test_transient
    ! fesm_datasrc, transient fields: interp exact at samples and linear
    ! between; interval_mean equal to the analytic mean for linear data and
    ! to the overlap-weighted mean with time_bnds; units converted (mm d-1 ->
    ! kg m-2 yr-1; an ACCUM field fed its mean rate); hold vs error;
    ! coverage stop at init; a repeated update reads nothing; time_ref.

    use precision
    use ncio
    use coordinates,  only : grid_class, grid_init
    use fesm_fields
    use fesm_datasrc

    implicit none

    character(len=*), parameter :: fpar = "test_ds_trans.nml", fdef = "par/datasrc_defaults.nml"
    integer, parameter :: nx = 3, ny = 2, nt = 5
    real(dp), parameter :: tol = 1e-4_dp
    type(grid_class)    :: grid
    type(field_state_t) :: st
    type(datasrc_t)     :: ds
    real(dp) :: tt(nt), tb(2,nt), v2(nx,ny), t0, t1, x
    real(wp) :: v(nx,ny,nt)
    integer  :: h, k, stat, nfail

    nfail = 0

    call grid_init(grid, name="TEST-DS-TR", mtype="cartesian", units="kilometers", &
                   x0=0.0_dp, dx=10.0_dp, nx=nx, y0=0.0_dp, dy=10.0_dp, ny=ny)

    ! Calendar years 1950 + [0, 10, 30, 60, 100] ("yr CE"); value = 10 + 2*(y-1950) + i
    tt = 1950.0_dp + [0.0_dp, 10.0_dp, 30.0_dp, 60.0_dp, 100.0_dp]
    call field_lin(v, tt)
    call ncfile("test_ds_trans_lin.nc", v, tt, "m", .false.)
    ! Same samples as cell means with bounds [1945,1955], [1955,1975], ...
    tb(1,:) = [1945.0_dp, 1955.0_dp, 1975.0_dp, 2000.0_dp, 2030.0_dp]
    tb(2,:) = [1955.0_dp, 1975.0_dp, 2000.0_dp, 2030.0_dp, 2070.0_dp]
    call ncfile("test_ds_trans_bnds.nc", v, tt, "mm d-1", .true.)

    call write_nml()

    call fields_init(st, "tr_test:exports")
    call fields_add(st, "H",    "m",           FIELD_INSTANT, FIELD_EXTENSIVE, h, grid="TEST-DS-TR", dims=[nx,ny])
    call fields_add(st, "Hm",   "m",           FIELD_MEAN,    FIELD_EXTENSIVE, h, grid="TEST-DS-TR", dims=[nx,ny])
    call fields_add(st, "smb",  "kg m-2 yr-1", FIELD_MEAN,    FIELD_EXTENSIVE, h, grid="TEST-DS-TR", dims=[nx,ny])
    call fields_add(st, "macc", "kg m-2",      FIELD_ACCUM,   FIELD_EXTENSIVE, h, grid="TEST-DS-TR", dims=[nx,ny])
    call fields_add(st, "Hh",   "m",           FIELD_INSTANT, FIELD_EXTENSIVE, h, grid="TEST-DS-TR", dims=[nx,ny])
    call fields_lock(st)

    ! --- INSTANT: interp ---------------------------------------------------------
    call datasrc_init(ds, fpar, "inst", st, [character(len=4) :: "H", "Hh"], grid, fdef, &
                      time_ref=1950.0_dp, time_span=[0.0_dp, 100.0_dp], stat=stat)
    call check("init (time_ref = 1950, span covered)", stat .eq. 0)
    call datasrc_describe(ds)
    call datasrc_update(ds, 30.0_dp)
    call datasrc_get(ds, "H", v2)
    call check("interp exact at a sample", maxval(abs(v2 - lin(30.0_dp))) .lt. tol)
    call datasrc_update(ds, 45.0_dp)
    call datasrc_get(ds, "H", v2)
    call check("interp linear between samples", maxval(abs(v2 - lin(45.0_dp))) .lt. tol)
    call datasrc_fill(ds, st)
    call fields_get(st, fields_handle(st, "H"), v2, t0, t1)
    call check("INSTANT validity [t1,t1]", t0 .eq. 45.0_dp .and. t1 .eq. 45.0_dp)

    ! hold vs error outside the axis
    call datasrc_update(ds, 0.0_dp, 120.0_dp, stat=stat)
    call check("error outside the axis: stat /= 0", stat .ne. 0)
    call datasrc_end(ds)
    call datasrc_init(ds, fpar, "hold", st, ["Hh"], grid, fdef, time_ref=1950.0_dp, stat=stat)
    call datasrc_update(ds, 150.0_dp)
    call datasrc_get(ds, "Hh", v2)
    call check("hold after the axis: last sample", maxval(abs(v2 - lin(100.0_dp))) .lt. tol)
    call datasrc_update(ds, -20.0_dp)
    call datasrc_get(ds, "Hh", v2)
    call check("hold before the axis: first sample", maxval(abs(v2 - lin(0.0_dp))) .lt. tol)
    call datasrc_end(ds)

    ! coverage stop at init
    call datasrc_init(ds, fpar, "inst", st, [character(len=4) :: "H", "Hh"], grid, fdef, &
                      time_ref=1950.0_dp, time_span=[-10.0_dp, 100.0_dp], stat=stat)
    call check("coverage stop at init: stat /= 0", stat .ne. 0)

    ! --- MEAN and ACCUM: interval_mean --------------------------------------------
    call datasrc_init(ds, fpar, "mean", st, [character(len=4) :: "Hm", "smb", "macc"], grid, fdef, &
                      time_ref=1950.0_dp, stat=stat)
    call check("init MEAN/ACCUM", stat .eq. 0)
    call datasrc_describe(ds)
    call datasrc_update(ds, 20.0_dp, 50.0_dp)
    call datasrc_get(ds, "Hm", v2)
    call check("interval_mean = analytic mean (linear data)", maxval(abs(v2 - lin(35.0_dp))) .lt. tol)
    ! smb: cells over bounds, mm d-1 -> kg m-2 yr-1 (x 365.25):
    ! [20,50] covers [20,25] of cell 2 (5 yr), [25,50] of cell 3 (25 yr)
    call datasrc_get(ds, "smb", v2)
    x = (5.0_dp*lin1(10.0_dp) + 25.0_dp*lin1(30.0_dp))/30.0_dp
    call check("interval_mean with time_bnds, converted", abs(v2(1,1) - 365.25_dp*x) .lt. 365.25_dp*tol)
    call datasrc_fill(ds, st)
    call fields_finalize_interval(st, 50.0_dp)
    call fields_get(st, fields_handle(st, "macc"), v2, t0, t1)
    call check("ACCUM: the state integrates the mean rate", abs(v2(1,1) - 30.0_dp*365.25_dp*x) .lt. 30.0_dp*365.25_dp*tol &
               .and. t0 .eq. 20.0_dp .and. t1 .eq. 50.0_dp)

    ! a repeated update reads nothing: the files may be gone
    call execute_command_line("mv test_ds_trans_lin.nc test_ds_trans_lin.nc.away; "// &
                              "mv test_ds_trans_bnds.nc test_ds_trans_bnds.nc.away")
    call datasrc_update(ds, 20.0_dp, 50.0_dp, stat=stat)
    call datasrc_get(ds, "Hm", v2)
    call check("repeated update reads nothing", stat .eq. 0 .and. maxval(abs(v2 - lin(35.0_dp))) .lt. tol)
    call execute_command_line("mv test_ds_trans_lin.nc.away test_ds_trans_lin.nc; "// &
                              "mv test_ds_trans_bnds.nc.away test_ds_trans_bnds.nc")

    ! MEAN with update(t) stops
    call datasrc_update(ds, 60.0_dp, stat=stat)
    call check("MEAN field with update(t): stat /= 0", stat .ne. 0)
    call datasrc_update(ds, 60.0_dp, 60.0_dp, stat=stat)
    call check("MEAN field with t0 = t1: stat /= 0", stat .ne. 0)
    call datasrc_end(ds)

    if (nfail .gt. 0) then
        write(*,*) "test_transient: ", nfail, " FAILED"
        error stop 1
    end if
    write(*,*) "test_transient: all passed"

contains

    function lin(t) result(f)
        ! The field at clock time t
        real(dp), intent(in) :: t
        real(dp) :: f(nx,ny)
        integer :: i
        do i = 1, nx
            f(i,:) = 10.0_dp + 2.0_dp*t + i
        end do
    end function lin

    real(dp) function lin1(t)
        real(dp), intent(in) :: t
        lin1 = 10.0_dp + 2.0_dp*t + 1.0_dp
    end function lin1

    subroutine field_lin(v, t)
        real(wp), intent(out) :: v(:,:,:)
        real(dp), intent(in)  :: t(:)
        integer :: k
        do k = 1, size(t)
            v(:,:,k) = real(lin(t(k) - 1950.0_dp), wp)
        end do
    end subroutine field_lin

    subroutine ncfile(fname, v, t, units, bnds)
        character(len=*), intent(in) :: fname, units
        real(wp),         intent(in) :: v(:,:,:)
        real(dp),         intent(in) :: t(:)
        logical,          intent(in) :: bnds
        call nc_create(fname)
        call nc_write_dim(fname, "xc", x=0.0_dp, dx=10.0_dp, nx=nx)
        call nc_write_dim(fname, "yc", x=0.0_dp, dx=10.0_dp, nx=ny)
        call nc_write_dim(fname, "time", x=t, units="yr CE", unlimited=.false.)
        call nc_write(fname, "v", v, dim1="xc", dim2="yc", dim3="time", units=units)
        if (bnds) then
            call nc_write_dim(fname, "nv", x=1.0_dp, dx=1.0_dp, nx=2)
            call nc_write(fname, "time_bnds", tb, dim1="nv", dim2="time")
            call nc_write_attr(fname, "time", "bounds", "time_bnds")
        end if
    end subroutine ncfile

    subroutine write_nml()
        integer :: u
        open(newunit=u, file=fpar, status="replace", action="write")
        write(u,"(a)") "&inst"
        write(u,"(a)") '    H  = "tr_lin"'
        write(u,"(a)") '    Hh = "tr_lin"'
        write(u,"(a)") "/"
        write(u,"(a)") "&hold"
        write(u,"(a)") '    Hh = "tr_lin_hold"'
        write(u,"(a)") "/"
        write(u,"(a)") "&mean"
        write(u,"(a)") '    Hm   = "tr_lin"'
        write(u,"(a)") '    smb  = "tr_bnds"'
        write(u,"(a)") '    macc = "tr_bnds"'
        write(u,"(a)") "/"
        write(u,"(a)") "&tr_lin"
        write(u,"(a)") '    filename = "test_ds_trans_lin.nc"'
        write(u,"(a)") '    name     = "v"'
        write(u,"(a)") "/"
        write(u,"(a)") "&tr_lin_hold"
        write(u,"(a)") '    filename     = "test_ds_trans_lin.nc"'
        write(u,"(a)") '    name         = "v"'
        write(u,"(a)") '    out_of_range = "hold"'
        write(u,"(a)") "/"
        write(u,"(a)") "&tr_bnds"
        write(u,"(a)") '    filename = "test_ds_trans_bnds.nc"'
        write(u,"(a)") '    name     = "v"'
        write(u,"(a)") "/"
        close(u)
    end subroutine write_nml

    subroutine check(label, ok)
        character(len=*), intent(in) :: label
        logical,          intent(in) :: ok
        write(*,"(a,a)") merge("  ok   ", "  FAIL ", ok), label
        if (.not. ok) nfail = nfail + 1
    end subroutine check

end program test_transient
