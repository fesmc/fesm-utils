program test_varslice_legacy
    ! The legacy path of varslice (varslice_init_nml without the new optional
    ! arguments, on groups with a `time` key, as in yelmox clim_db.nml) gives
    ! byte-identical vs%var and vs%time across the datasrc preparatory steps
    ! (V1-V5). The groups and files are written here; every case is reduced
    ! to its shape and a hash of the bit patterns of var and time, compared
    ! with test/datasrc/test_varslice_legacy.ref.
    !
    !     ./bin/test_varslice_legacy.x          compare with the reference
    !     ./bin/test_varslice_legacy.x record   (re)write the reference
    !
    ! The reference was recorded with varslice before V1 (fesm-utils 1725385).

    use, intrinsic :: iso_fortran_env, only : int32, int64
    use precision
    use constants,   only : mv
    use ncio
    use coordinates, only : grid_class, grid_init
    use varslice,    only : varslice_class, varslice_init_nml, varslice_update, &
                            varslice_sub_mean, varslice_end

    implicit none

    character(len=*), parameter :: fnml = "test_ds_legacy.nml"
    character(len=*), parameter :: fref = "test/datasrc/test_varslice_legacy.ref"
    integer, parameter :: nx = 5, ny = 4, nr = 64

    type(varslice_class) :: vs
    type(grid_class) :: grid_tgt
    character(len=160) :: rec(nr), ref(nr), arg
    real(wp), allocatable :: sm(:,:,:)
    integer :: n, nref, i, u, ios, nfail
    logical :: record

    call get_command_argument(1, arg)
    record = (trim(arg) .eq. "record")

    call write_files()
    call write_nml()

    n = 0

    ! --- transient 2-D, regular axis (ice history) --------------------------
    call varslice_init_nml(vs, fnml, "leg_tr", verbose=.false.)
    call add("tr init")
    call varslice_update(vs, [-3000.0_dp], method="exact");              call add("tr exact")
    call varslice_update(vs, [-2300.0_dp], method="interp");             call add("tr interp")
    call varslice_update(vs, [-5000.0_dp], method="interp");             call add("tr interp edge")
    call varslice_update(vs, [-6000.0_dp], method="extrap");             call add("tr extrap lo")
    call varslice_update(vs, [  500.0_dp], method="extrap");             call add("tr extrap hi")
    call varslice_update(vs, [  500.0_dp], method="interp");             call add("tr interp out")
    call varslice_update(vs, [-4000.0_dp,-2000.0_dp], method="range");   call add("tr range")
    call varslice_update(vs, [-5000.0_dp,0.0_dp], method="range_mean", rep=1); call add("tr range_mean")
    call varslice_update(vs, [-4000.0_dp,-1000.0_dp], method="range_min", rep=1); call add("tr range_min")
    call varslice_update(vs, [-4000.0_dp,-1000.0_dp], method="range_max", rep=1); call add("tr range_max")
    call varslice_update(vs, [-4000.0_dp,-1000.0_dp], method="range_sd",  rep=1); call add("tr range_sd")
    call varslice_update(vs, [-3000.0_dp,-1000.0_dp], method="range_sum", rep=1); call add("tr range_sum")
    call varslice_end(vs)

    call varslice_init_nml(vs, fnml, "leg_tr_sc", verbose=.false.)
    call varslice_update(vs, [-1500.0_dp], method="interp");             call add("tr_sc interp")
    call varslice_end(vs)

    ! --- transient monthly (sub-annual axis) --------------------------------
    call varslice_init_nml(vs, fnml, "leg_mon", verbose=.false.)
    call add("mon init")
    call varslice_update(vs, [2001.0_dp], method="interp", rep=12);      call add("mon interp rep12")
    call varslice_update(vs, [2003.0_dp], method="extrap", rep=12);      call add("mon extrap rep12")
    call varslice_update(vs, [2000.0_dp,2002.0_dp], method="range_mean", rep=12); call add("mon range_mean rep12")
    call varslice_update(vs, [2000.0_dp,2001.0_dp], method="range_mean", rep=1);  call add("mon range_mean rep1")
    call varslice_update(vs, [2001.0_dp,2001.0_dp], method="range", rep=12);      call add("mon range rep12")
    call varslice_end(vs)

    ! --- climatology (cycle without years) ----------------------------------
    call varslice_init_nml(vs, fnml, "leg_clim", verbose=.false.)
    call varslice_update(vs);                                            call add("clim cycle")
    sm = varslice_sub_mean(vs)
    call add_arr("clim sub_mean", sm)
    call varslice_update(vs, [0.0_dp,0.0_dp], method="range_mean", rep=1); call add("clim collapsed")
    call varslice_end(vs)

    ! --- static with scaling and fill ---------------------------------------
    call varslice_init_nml(vs, fnml, "leg_static", verbose=.false.)
    call varslice_update(vs);                                            call add("static nearest")
    call varslice_update(vs, fill="none");                               call add("static none")
    call varslice_end(vs)
    call varslice_init_nml(vs, fnml, "leg_static_c", verbose=.false.)
    call varslice_update(vs);                                            call add("static const")
    call varslice_end(vs)

    ! --- multifile ----------------------------------------------------------
    call varslice_init_nml(vs, fnml, "leg_mf", verbose=.false.)
    call add("mf init")
    call varslice_update(vs, [4.0_dp], method="exact");                  call add("mf exact")
    call varslice_update(vs, [2.0_dp,5.0_dp], method="range");           call add("mf range")
    call varslice_update(vs, [3.5_dp], method="interp");                 call add("mf interp")
    call varslice_end(vs)

    ! --- online remap from a lon-lat file (snapshot style) ------------------
    call grid_init(grid_tgt, name="TEST-LL5", mtype="latlon", units="degrees", &
                   x0=-10.0_dp, dx=5.0_dp, nx=nx, y0=40.0_dp, dy=5.0_dp, ny=ny)
    call varslice_init_nml(vs, fnml, "leg_remap", verbose=.false., grid=grid_tgt)
    call varslice_update(vs);                                            call add("remap bilin")
    call varslice_end(vs)

    ! --- record or compare ---------------------------------------------------
    if (record) then
        open(newunit=u, file=fref, status="replace", action="write")
        write(u,"(a)") "# test_varslice_legacy reference: label | shape | hash(var) | hash(time)"
        do i = 1, n
            write(u,"(a)") trim(rec(i))
        end do
        close(u)
        write(*,*) "test_varslice_legacy: reference written to "//fref//" (", n, " cases)"
        stop
    end if

    open(newunit=u, file=fref, status="old", action="read", iostat=ios)
    if (ios .ne. 0) then
        write(*,*) "test_varslice_legacy: cannot open the reference "//fref
        error stop 1
    end if
    read(u,"(a)")    ! header
    nref = 0
    do
        read(u,"(a)",iostat=ios) arg
        if (ios .ne. 0) exit
        nref = nref + 1
        ref(nref) = arg
    end do
    close(u)

    nfail = 0
    if (nref .ne. n) then
        write(*,*) "  FAIL number of cases: ", n, " reference: ", nref
        nfail = nfail + 1
    end if
    do i = 1, min(n, nref)
        if (trim(rec(i)) .eq. trim(ref(i))) then
            write(*,"(a)") "  ok   "//trim(rec(i))
        else
            write(*,"(a)") "  FAIL "//trim(rec(i))
            write(*,"(a)") "  ref  "//trim(ref(i))
            nfail = nfail + 1
        end if
    end do

    if (nfail .gt. 0) then
        write(*,*) "test_varslice_legacy: ", nfail, " FAILED"
        error stop 1
    end if
    write(*,*) "test_varslice_legacy: all passed"

contains

    subroutine add(label)
        ! One record: label, shape of var, hash of var, hash of time.
        character(len=*), intent(in) :: label
        character(len=64) :: shp
        n = n + 1
        if (allocated(vs%var)) then
            write(shp,"(i0,'x',i0,'x',i0,'x',i0)") shape(vs%var)
            rec(n) = label//" | "//trim(shp)//" | "//hash_sp(reshape(vs%var,[size(vs%var)]))// &
                     " | "//hash_dp(vs%time)
        else
            rec(n) = label//" | - | - | "//hash_dp(vs%time)
        end if
    end subroutine add

    subroutine add_arr(label, a)
        character(len=*), intent(in) :: label
        real(wp),         intent(in) :: a(:,:,:)
        character(len=64) :: shp
        n = n + 1
        write(shp,"(i0,'x',i0,'x',i0)") shape(a)
        rec(n) = label//" | "//trim(shp)//" | "//hash_sp(reshape(a,[size(a)]))//" | -"
    end subroutine add_arr

    function hash_sp(x) result(s)
        ! Two polynomial hashes (mod two primes) of the 32-bit patterns.
        real(sp), intent(in) :: x(:)
        character(len=17) :: s
        integer(int64) :: h1, h2, w
        integer :: i
        h1 = 7; h2 = 11
        do i = 1, size(x)
            w  = iand(int(transfer(x(i), 0_int32), int64), int(z'FFFFFFFF', int64))
            h1 = mod(h1*1000003_int64 + w, 2147483647_int64)
            h2 = mod(h2*1000033_int64 + w, 2147483629_int64)
        end do
        write(s,"(z8.8,'-',z8.8)") h1, h2
    end function hash_sp

    function hash_dp(x) result(s)
        real(dp), intent(in) :: x(:)
        character(len=17) :: s
        integer(int64) :: h1, h2, w, b
        integer :: i, k
        h1 = 7; h2 = 11
        do i = 1, size(x)
            b = transfer(x(i), 0_int64)
            do k = 0, 1
                w  = iand(ishft(b, -32*k), int(z'FFFFFFFF', int64))
                h1 = mod(h1*1000003_int64 + w, 2147483647_int64)
                h2 = mod(h2*1000033_int64 + w, 2147483629_int64)
            end do
        end do
        write(s,"(z8.8,'-',z8.8)") h1, h2
    end function hash_dp

    pure real(sp) function val(i, j, k)
        integer, intent(in) :: i, j, k
        val = real(100.0_dp + 10.0_dp*i + j + 0.37_dp*k*sin(real(i+j,dp)) + 0.01_dp*k*k, sp)
    end function val

    subroutine write_files()
        real(sp) :: v3(nx,ny,36), v2(nx,ny), vll(8,6)
        integer  :: i, j, k

        ! transient, 6 steps
        do k = 1, 6; do j = 1, ny; do i = 1, nx
            v3(i,j,k) = val(i,j,k)
        end do; end do; end do
        call nc_create("test_ds_legacy_tr.nc")
        call nc_write_dim("test_ds_legacy_tr.nc", "xc", x=1.0_dp, dx=1.0_dp, nx=nx)
        call nc_write_dim("test_ds_legacy_tr.nc", "yc", x=1.0_dp, dx=1.0_dp, nx=ny)
        call nc_write_dim("test_ds_legacy_tr.nc", "time", x=0.0_dp, dx=1.0_dp, nx=6)
        call nc_write("test_ds_legacy_tr.nc", "H_ice", v3(:,:,1:6), dim1="xc", dim2="yc", dim3="time", &
                      units="m")

        ! monthly, 3 years
        do k = 1, 36; do j = 1, ny; do i = 1, nx
            v3(i,j,k) = val(i,j,k) + real(mod(k-1,12),sp)
        end do; end do; end do
        call nc_create("test_ds_legacy_mon.nc")
        call nc_write_dim("test_ds_legacy_mon.nc", "xc", x=1.0_dp, dx=1.0_dp, nx=nx)
        call nc_write_dim("test_ds_legacy_mon.nc", "yc", x=1.0_dp, dx=1.0_dp, nx=ny)
        call nc_write_dim("test_ds_legacy_mon.nc", "time", x=1.0_dp, dx=1.0_dp, nx=36)
        call nc_write("test_ds_legacy_mon.nc", "tas", v3, dim1="xc", dim2="yc", dim3="time", units="K")

        ! climatology, 12 months, one cell missing
        do k = 1, 12; do j = 1, ny; do i = 1, nx
            v3(i,j,k) = val(i,j,k) - 100.0_sp
        end do; end do; end do
        v3(2,3,1:12) = mv
        call nc_create("test_ds_legacy_clim.nc")
        call nc_write_dim("test_ds_legacy_clim.nc", "xc", x=1.0_dp, dx=1.0_dp, nx=nx)
        call nc_write_dim("test_ds_legacy_clim.nc", "yc", x=1.0_dp, dx=1.0_dp, nx=ny)
        call nc_write_dim("test_ds_legacy_clim.nc", "month", x=1.0_dp, dx=1.0_dp, nx=12)
        call nc_write("test_ds_legacy_clim.nc", "to", v3(:,:,1:12), dim1="xc", dim2="yc", dim3="month", &
                      units="degC", missing_value=mv)

        ! static with missing cells
        do j = 1, ny; do i = 1, nx
            v2(i,j) = val(i,j,3)*100.0_sp
        end do; end do
        v2(1,1) = mv
        v2(4,2) = mv
        call nc_create("test_ds_legacy_static.nc")
        call nc_write_dim("test_ds_legacy_static.nc", "xc", x=1.0_dp, dx=1.0_dp, nx=nx)
        call nc_write_dim("test_ds_legacy_static.nc", "yc", x=1.0_dp, dx=1.0_dp, nx=ny)
        call nc_write("test_ds_legacy_static.nc", "sf", v2, dim1="xc", dim2="yc", units="kg m-2 yr-1", &
                      missing_value=mv)

        ! multifile: 3 + 3 steps
        do k = 1, 6; do j = 1, ny; do i = 1, nx
            v3(i,j,k) = val(i,j,k) + 7.0_sp
        end do; end do; end do
        call nc_create("test_ds_legacy_mf_01.nc")
        call nc_write_dim("test_ds_legacy_mf_01.nc", "xc", x=1.0_dp, dx=1.0_dp, nx=nx)
        call nc_write_dim("test_ds_legacy_mf_01.nc", "yc", x=1.0_dp, dx=1.0_dp, nx=ny)
        call nc_write_dim("test_ds_legacy_mf_01.nc", "time", x=1.0_dp, dx=1.0_dp, nx=3)
        call nc_write("test_ds_legacy_mf_01.nc", "v", v3(:,:,1:3), dim1="xc", dim2="yc", dim3="time")
        call nc_create("test_ds_legacy_mf_02.nc")
        call nc_write_dim("test_ds_legacy_mf_02.nc", "xc", x=1.0_dp, dx=1.0_dp, nx=nx)
        call nc_write_dim("test_ds_legacy_mf_02.nc", "yc", x=1.0_dp, dx=1.0_dp, nx=ny)
        call nc_write_dim("test_ds_legacy_mf_02.nc", "time", x=4.0_dp, dx=1.0_dp, nx=3)
        call nc_write("test_ds_legacy_mf_02.nc", "v", v3(:,:,4:6), dim1="xc", dim2="yc", dim3="time")

        ! lon-lat source of the remap case, 8 x 6 cells of 3 degrees
        do j = 1, 6; do i = 1, 8
            vll(i,j) = val(i,j,1)
        end do; end do
        call nc_create("test_ds_legacy_ll.nc")
        call nc_write_dim("test_ds_legacy_ll.nc", "lon", x=-12.0_dp, dx=3.0_dp, nx=8)
        call nc_write_dim("test_ds_legacy_ll.nc", "lat", x=38.0_dp, dx=3.0_dp, nx=6)
        call nc_write("test_ds_legacy_ll.nc", "so", vll, dim1="lon", dim2="lat", units="PSU")

    end subroutine write_files

    subroutine write_nml()
        integer :: u
        open(newunit=u, file=fnml, status="replace", action="write")
        write(u,"(a)") "! test_varslice_legacy: clim_db-style groups (written by the test)"
        call grp(u, "leg_tr", "test_ds_legacy_tr.nc", "H_ice", '"m" "m"', "1.0 0.0", "1 -5000 0 1000 1")
        call grp(u, "leg_tr_sc", "test_ds_legacy_tr.nc", "H_ice", '"degC" "K"', "1.0 273.15", "1 -5000 0 1000 1")
        call grp(u, "leg_mon", "test_ds_legacy_mon.nc", "tas", '"K" "K"', "1.0 0.0", "1 2000 2002 1 12")
        call grp(u, "leg_clim", "test_ds_legacy_clim.nc", "to", '"degC" "K"', "1.0 273.15", "0 0 0 0 12")
        call grp(u, "leg_static", "test_ds_legacy_static.nc", "sf", '"kg m-2 yr-1" "mm d-1"', &
                 "2.73785e-3 0.0", "0 1961 1990 0 1", extra='    fill     = "nearest"')
        call grp(u, "leg_static_c", "test_ds_legacy_static.nc", "sf", '"kg m-2 yr-1" "mm d-1"', &
                 "2.73785e-3 0.0", "0 1961 1990 0 1", extra='    fill     = "const"'//new_line("a")// &
                 '    fill_value = 5.0')
        call grp(u, "leg_mf", "test_ds_legacy_mf_*.nc", "v", '"m" "m"', "1.0 0.0", "1 1 6 1 1")
        call grp(u, "leg_remap", "test_ds_legacy_ll.nc", "so", '"PSU" "PSU"', "1.0 0.0", "0 0 0 0 1", &
                 extra='    remap    = "bilin"'//new_line("a")//'    grid_src = "TEST-LL8"'// &
                       new_line("a")//'    fill     = "nearest"')
        close(u)
    end subroutine write_nml

    subroutine grp(u, name, file, var, units, scaling, time, extra)
        integer,          intent(in) :: u
        character(len=*), intent(in) :: name, file, var, units, scaling, time
        character(len=*), intent(in), optional :: extra
        write(u,"(a)") ""
        write(u,"(a)") "&"//name
        write(u,"(a)") '    filename = "'//file//'"'
        write(u,"(a)") '    name     = "'//var//'"'
        write(u,"(a)") '    units    = '//units
        write(u,"(a)") '    scaling  = '//scaling
        write(u,"(a)") '    time     = '//time
        if (present(extra)) write(u,"(a)") extra
        write(u,"(a)") "/"
    end subroutine grp

end program test_varslice_legacy
