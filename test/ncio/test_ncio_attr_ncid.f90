program test_ncio_attr_ncid
    ! Attribute routines on an open file (ncid forms of nc_write_attr,
    ! nc_read_attr, nc_exists_attr, nc_size_attr; ncid of nc_exists_var and
    ! nc_dims):
    !   * every attribute type, variable and global scope, scalar and vector,
    !     written through an ncid and read back through the filename forms
    !     (which open the file themselves) and through an ncid;
    !   * attributes written by the filename forms read through an ncid;
    !   * an attribute written twice keeps the last value;
    !   * the file opened once can mix nc_write_dim, nc_write and attributes;
    !   * exists/size/dims through an ncid agree with the filename forms.
    ! The file is written by the test itself, then removed.

    use ncio
    use netcdf

    implicit none

    character(len=*), parameter :: fnm = "test_ncio_attr_ncid.nc"
    integer, parameter :: nx = 5
    integer :: nfail, ncid, u

    double precision :: v(nx), dv(3), dr(3)
    real(4)          :: sv(3), sr(3)
    integer          :: iv(3), ir(3)
    character(len=64) :: cs
    integer :: i0
    real(4) :: s0
    double precision :: d0
    integer, allocatable :: dl(:)

    nfail = 0
    v  = [(0.5d0*i0, i0 = 1, nx)]
    iv = [3, -7, 1000000]
    sv = [1.5, -2.25, 3.0e10]
    dv = [1.0d0/3.0d0, -2.0d300, 4.9d-300]

    write(*,*) "==================================================="
    write(*,*) " ncio attribute ncid test suite"
    write(*,*) "==================================================="

    ! --- One open file: dimension, variable, every attribute ----------------
    call nc_create(fnm)
    call nc_open(fnm, ncid, writable=.true.)
    call nc_write_dim(fnm, "x", 1.0d0, dx=1.0d0, nx=nx, ncid=ncid)
    call nc_write(fnm, "v", v, dim1="x", ncid=ncid)

    call nc_write_attr(fnm, ncid, "v", "a_char", "some text")
    call nc_write_attr(fnm, ncid, "v", "a_int",  -42)
    call nc_write_attr(fnm, ncid, "v", "a_sp",   0.1)
    call nc_write_attr(fnm, ncid, "v", "a_dp",   0.1d0)
    call nc_write_attr(fnm, ncid, "v", "a_int1", iv)
    call nc_write_attr(fnm, ncid, "v", "a_sp1",  sv)
    call nc_write_attr(fnm, ncid, "v", "a_dp1",  dv)
    call nc_write_attr(fnm, ncid, "g_char", "global text")
    call nc_write_attr(fnm, ncid, "g_int",  17)
    call nc_write_attr(fnm, ncid, "g_sp",   0.2)
    call nc_write_attr(fnm, ncid, "g_dp",   0.2d0)
    call nc_write_attr(fnm, ncid, "g_int1", iv)
    call nc_write_attr(fnm, ncid, "g_sp1",  sv)
    call nc_write_attr(fnm, ncid, "g_dp1",  dv)

    ! Overwritten: the last value stays
    call nc_write_attr(fnm, ncid, "v", "a_int", 43)
    call nc_write_attr(fnm, ncid, "g_char", "replaced")

    ! Data written after attributes still works (data mode kept)
    call nc_write(fnm, "w", 2.0d0*v, dim1="x", ncid=ncid)
    call nc_write_attr(fnm, ncid, "w", "units", "m")
    call nc_close(ncid)
    call check(ncid == 0, "nc_close returns ncid = 0")

    ! --- Read back through the filename forms -------------------------------
    cs = ""; call nc_read_attr(fnm, "v", "a_char", cs)
    call check(trim(cs) == "some text", "variable char")
    call nc_read_attr(fnm, "v", "a_int", i0)
    call check(i0 == 43, "variable int, overwritten")
    call nc_read_attr(fnm, "v", "a_sp", s0)
    call check(s0 == 0.1, "variable sp")
    call nc_read_attr(fnm, "v", "a_dp", d0)
    call check(d0 == 0.1d0, "variable dp")
    call nc_read_attr(fnm, "v", "a_int1", ir)
    call check(all(ir == iv), "variable int 1D")
    call nc_read_attr(fnm, "v", "a_sp1", sr)
    call check(all(sr == sv), "variable sp 1D")
    call nc_read_attr(fnm, "v", "a_dp1", dr)
    call check(all(dr == dv), "variable dp 1D")
    cs = ""; call nc_read_attr(fnm, "g_char", cs)
    call check(trim(cs) == "replaced", "global char, overwritten")
    call nc_read_attr(fnm, "g_int", i0)
    call check(i0 == 17, "global int")
    call nc_read_attr(fnm, "g_sp", s0)
    call check(s0 == 0.2, "global sp")
    call nc_read_attr(fnm, "g_dp", d0)
    call check(d0 == 0.2d0, "global dp")
    call nc_read_attr(fnm, "g_int1", ir)
    call check(all(ir == iv), "global int 1D")
    call nc_read_attr(fnm, "g_sp1", sr)
    call check(all(sr == sv), "global sp 1D")
    call nc_read_attr(fnm, "g_dp1", dr)
    call check(all(dr == dv), "global dp 1D")
    cs = ""; call nc_read_attr(fnm, "w", "units", cs)
    call check(trim(cs) == "m", "attribute after data written through the ncid")

    ! Data written through the same ncid is intact and read back
    call check_data("v", v)
    call check_data("w", 2.0d0*v)

    ! --- Read back through an ncid (read-only) ------------------------------
    call nc_open(fnm, ncid, writable=.false.)
    cs = ""; call nc_read_attr(fnm, ncid, "v", "a_char", cs)
    call check(trim(cs) == "some text", "ncid variable char")
    call nc_read_attr(fnm, ncid, "v", "a_int", i0)
    call check(i0 == 43, "ncid variable int")
    call nc_read_attr(fnm, ncid, "v", "a_sp", s0)
    call check(s0 == 0.1, "ncid variable sp")
    call nc_read_attr(fnm, ncid, "v", "a_dp", d0)
    call check(d0 == 0.1d0, "ncid variable dp")
    call nc_read_attr(fnm, ncid, "v", "a_int1", ir)
    call check(all(ir == iv), "ncid variable int 1D")
    call nc_read_attr(fnm, ncid, "v", "a_sp1", sr)
    call check(all(sr == sv), "ncid variable sp 1D")
    call nc_read_attr(fnm, ncid, "v", "a_dp1", dr)
    call check(all(dr == dv), "ncid variable dp 1D")
    cs = ""; call nc_read_attr(fnm, ncid, "g_char", cs)
    call check(trim(cs) == "replaced", "ncid global char")
    call nc_read_attr(fnm, ncid, "g_int", i0)
    call check(i0 == 17, "ncid global int")
    call nc_read_attr(fnm, ncid, "g_sp", s0)
    call check(s0 == 0.2, "ncid global sp")
    call nc_read_attr(fnm, ncid, "g_dp", d0)
    call check(d0 == 0.2d0, "ncid global dp")
    call nc_read_attr(fnm, ncid, "g_int1", ir)
    call check(all(ir == iv), "ncid global int 1D")
    call nc_read_attr(fnm, ncid, "g_sp1", sr)
    call check(all(sr == sv), "ncid global sp 1D")
    call nc_read_attr(fnm, ncid, "g_dp1", dr)
    call check(all(dr == dv), "ncid global dp 1D")

    ! exists / size / dims / variable through the ncid agree with the filename forms
    call check(nc_exists_attr(fnm, ncid, "g_int") .eqv. nc_exists_attr(fnm, "g_int"), "exists global, present")
    call check(nc_exists_attr(fnm, ncid, "g_none") .eqv. nc_exists_attr(fnm, "g_none"), "exists global, absent")
    call check(nc_exists_attr(fnm, ncid, "v", "a_dp1"), "ncid exists variable, present")
    call check(.not. nc_exists_attr(fnm, ncid, "v", "a_none"), "ncid exists variable, absent")
    call check(nc_size_attr(fnm, ncid, "g_dp1") == 3, "ncid size global")
    call check(nc_size_attr(fnm, ncid, "v", "a_sp1") == nc_size_attr(fnm, "v", "a_sp1"), "ncid size variable")
    call check(nc_size_attr(fnm, ncid, "v", "a_char") == len("some text"), "ncid size char")
    call check(nc_exists_var(fnm, "v", ncid=ncid), "ncid exists var, present")
    call check(.not. nc_exists_var(fnm, "nope", ncid=ncid), "ncid exists var, absent")
    call nc_dims(fnm, "v", dims=dl, ncid=ncid)
    call check(size(dl) == 1 .and. dl(1) == nx, "ncid dims")
    call nc_close(ncid)

    ! --- Written by the filename forms, read through an ncid ----------------
    call nc_write_attr(fnm, "v", "b_int", 5)
    call nc_write_attr(fnm, "b_dp", 2.5d0)
    call nc_open(fnm, ncid, writable=.false.)
    call nc_read_attr(fnm, ncid, "v", "b_int", i0)
    call check(i0 == 5, "filename write, ncid read: variable")
    call nc_read_attr(fnm, ncid, "b_dp", d0)
    call check(d0 == 2.5d0, "filename write, ncid read: global")
    call nc_close(ncid)

    open(newunit=u, file=fnm, status="old")
    close(u, status="delete")

    write(*,"(a)") ""
    if (nfail == 0) then
        write(*,*) "All ncio attribute ncid tests passed."
    else
        write(*,*) "ncio attribute ncid tests FAILED: ", nfail
        stop 1
    end if

contains

    subroutine check(cond, label)
        logical,          intent(in) :: cond
        character(len=*), intent(in) :: label
        if (.not. cond) then
            write(*,"(a,a)") "  FAIL  ", label
            nfail = nfail + 1
        end if
    end subroutine check

    subroutine check_data(name, ref)
        character(len=*), intent(in) :: name
        double precision, intent(in) :: ref(:)
        double precision :: r(nx)
        call nc_read(fnm, name, r)
        call check(all(r == ref), "data "//name)
    end subroutine check_data

end program test_ncio_attr_ncid
