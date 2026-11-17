program test_varslice_legacy
    ! The legacy path of varslice (varslice_init_nml without the new optional
    ! arguments, on groups with a `time` key, as in yelmox clim_db.nml) keeps
    ! vs%var and vs%time across the datasrc steps. The groups and files are
    ! written here; every case is reduced to its shape, the counts of valid and
    ! missing values, and min, max, sum and an index-weighted sum of var and of
    ! time, compared with test/datasrc/test_varslice_legacy.ref: counts exactly,
    ! the values to rtol_var (sp) and rtol_time (dp) of the field's magnitude.
    ! The index weight makes a permuted field fail.
    !
    ! A tolerance, not bit patterns: the last bits of a mean depend on the
    ! platform and the flags (fused multiply-add on arm64, -ffast-math or the
    ! value-safe flags), so a bit hash only held on the machine that recorded it.
    !
    !     ./bin/test_varslice_legacy.x          compare with the reference
    !     ./bin/test_varslice_legacy.x record   (re)write the reference
    !
    ! The reference was recorded at v1.5 with gfortran -O2 on arm64, where every
    ! case is bit-identical to the hashes of varslice before V1 (fesm-utils 1725385).

    use precision
    use constants,   only : mv
    use ncio
    use coordinates, only : grid_class, grid_init
    use varslice,    only : varslice_class, varslice_init_nml, varslice_update, &
                            varslice_sub_mean, varslice_end

    implicit none

    type case_t
        character(len=40) :: label = ""
        character(len=24) :: shp   = "-"
        logical  :: has_var  = .false.
        integer  :: nv = 0, nmv = 0      ! valid and missing values of var
        real(dp) :: v(4) = 0.0_dp        ! var: min, max, sum, index-weighted sum (valid values)
        logical  :: has_time = .false.
        integer  :: nt = 0
        real(dp) :: t(4) = 0.0_dp        ! time: min, max, sum, index-weighted sum
    end type case_t

    character(len=*), parameter :: fnml = "test_ds_legacy.nml"
    character(len=*), parameter :: fref = "test/datasrc/test_varslice_legacy.ref"
    integer,  parameter :: nx = 5, ny = 4, nr = 64
    real(dp), parameter :: rtol_var  = 1.0e-6_dp    ! ~8 ulp of sp
    real(dp), parameter :: rtol_time = 1.0e-12_dp

    type(varslice_class) :: vs
    type(grid_class) :: grid_tgt
    type(case_t) :: rec(nr), ref(nr)
    character(len=400) :: line
    character(len=160) :: arg
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
        write(u,"(a)") "# test_varslice_legacy reference: label | shape | var: nvalid nmissing min max sum wsum"// &
                       " | time: n min max sum wsum"
        do i = 1, n
            write(u,"(a)") trim(case_line(rec(i)))
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
        read(u,"(a)",iostat=ios) line
        if (ios .ne. 0) exit
        nref = nref + 1
        ref(nref) = case_parse(line)
    end do
    close(u)

    nfail = 0
    if (nref .ne. n) then
        write(*,*) "  FAIL number of cases: ", n, " reference: ", nref
        nfail = nfail + 1
    end if
    do i = 1, min(n, nref)
        if (case_match(rec(i), ref(i))) then
            write(*,"(a)") "  ok   "//trim(rec(i)%label)//" | "//trim(rec(i)%shp)
        else
            write(*,"(a)") "  FAIL "//trim(case_line(rec(i)))
            write(*,"(a)") "  ref  "//trim(case_line(ref(i)))
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
        ! One case: label, shape and statistics of vs%var, statistics of vs%time.
        character(len=*), intent(in) :: label
        integer :: nm
        n = n + 1
        rec(n)%label = label
        if (allocated(vs%var)) then
            write(rec(n)%shp,"(i0,'x',i0,'x',i0,'x',i0)") shape(vs%var)
            rec(n)%has_var = .true.
            call stats(real(reshape(vs%var,[size(vs%var)]),dp), &
                       reshape(vs%var .ne. mv,[size(vs%var)]), rec(n)%nv, rec(n)%nmv, rec(n)%v)
        end if
        rec(n)%has_time = .true.
        call stats(vs%time, spread(.true., 1, size(vs%time)), rec(n)%nt, nm, rec(n)%t)
    end subroutine add

    subroutine add_arr(label, a)
        character(len=*), intent(in) :: label
        real(wp),         intent(in) :: a(:,:,:)
        n = n + 1
        rec(n)%label = label
        write(rec(n)%shp,"(i0,'x',i0,'x',i0)") shape(a)
        rec(n)%has_var = .true.
        call stats(real(reshape(a,[size(a)]),dp), reshape(a .ne. mv,[size(a)]), &
                   rec(n)%nv, rec(n)%nmv, rec(n)%v)
    end subroutine add_arr

    subroutine stats(x, valid, nv, nmv, v)
        ! min, max, sum and index-weighted sum (weight i/size) of the valid values.
        real(dp), intent(in)  :: x(:)
        logical,  intent(in)  :: valid(:)
        integer,  intent(out) :: nv, nmv
        real(dp), intent(out) :: v(4)
        integer :: k
        nv  = count(valid)
        nmv = size(x) - nv
        v   = 0.0_dp
        if (nv .eq. 0) return
        v(1) = minval(x, mask=valid)
        v(2) = maxval(x, mask=valid)
        do k = 1, size(x)
            if (.not. valid(k)) cycle
            v(3) = v(3) + x(k)
            v(4) = v(4) + x(k)*real(k,dp)/real(size(x),dp)
        end do
    end subroutine stats

    logical function stats_close(a, b, n, rtol)
        ! a ~ b to rtol of the field's magnitude; the sums to rtol of n times it.
        real(dp), intent(in) :: a(4), b(4), rtol
        integer,  intent(in) :: n
        real(dp) :: tol
        tol = rtol*max(abs(b(1)), abs(b(2)), tiny(1.0_dp))
        stats_close = all(abs(a(1:2) - b(1:2)) .le. tol) .and. &
                      all(abs(a(3:4) - b(3:4)) .le. tol*real(max(n,1),dp))
    end function stats_close

    logical function case_match(a, b)
        type(case_t), intent(in) :: a, b
        case_match = trim(a%label) .eq. trim(b%label) .and. trim(a%shp) .eq. trim(b%shp) .and. &
                     (a%has_var .eqv. b%has_var) .and. (a%has_time .eqv. b%has_time)
        if (.not. case_match) return
        if (a%has_var) case_match = a%nv .eq. b%nv .and. a%nmv .eq. b%nmv .and. &
                                    stats_close(a%v, b%v, b%nv, rtol_var)
        if (.not. case_match) return
        if (a%has_time) case_match = a%nt .eq. b%nt .and. stats_close(a%t, b%t, b%nt, rtol_time)
    end function case_match

    function case_line(c) result(s)
        type(case_t), intent(in) :: c
        character(len=400) :: s
        character(len=200) :: sv, st
        sv = "-"
        st = "-"
        if (c%has_var)  write(sv,"(i0,1x,i0,4(1x,es24.16e3))") c%nv, c%nmv, c%v
        if (c%has_time) write(st,"(i0,4(1x,es24.16e3))") c%nt, c%t
        s = trim(c%label)//" | "//trim(c%shp)//" | "//trim(sv)//" | "//trim(st)
    end function case_line

    function case_parse(s) result(c)
        character(len=*), intent(in) :: s
        type(case_t) :: c
        character(len=400) :: f(4)
        integer :: k, p, q
        p = 1
        do k = 1, 3
            q = index(s(p:), "|")
            f(k) = adjustl(s(p:p+q-2))
            p = p + q
        end do
        f(4) = adjustl(s(p:))
        c%label = trim(f(1))
        c%shp   = trim(f(2))
        c%has_var = trim(f(3)) .ne. "-"
        if (c%has_var) read(f(3),*) c%nv, c%nmv, c%v
        c%has_time = trim(f(4)) .ne. "-"
        if (c%has_time) read(f(4),*) c%nt, c%t
    end function case_parse

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
