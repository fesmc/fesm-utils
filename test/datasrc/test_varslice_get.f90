program test_varslice_get
    ! varslice_get (varslice V4): dp and sp copies of the slice at its
    ! natural rank, 0-D to 3-D, at a trailing index; dp widening is exact.

    use precision
    use ncio
    use constants, only : mv
    use varslice,  only : varslice_class, varslice_init_arg, varslice_update, varslice_get, varslice_end

    implicit none

    integer, parameter :: nx = 4, ny = 3, nz = 2, nt = 5
    type(varslice_class) :: vs
    real(wp) :: v4(nx,ny,nz,nt), v3(nx,ny,nt), v1(nt)
    real(dp) :: d0, d1(nx), d2(nx,ny), d3(nx,ny,nz)
    real(sp) :: s0, s2(nx,ny), s3(nx,ny,nz)
    integer :: i, j, k, t, nfail

    nfail = 0

    do t = 1, nt; do k = 1, nz; do j = 1, ny; do i = 1, nx
        v4(i,j,k,t) = real(1000*t + 100*k + 10*j + i, wp) + 0.1_wp
    end do; end do; end do; end do
    v3 = v4(:,:,1,:)
    v3(2,2,3) = mv
    do t = 1, nt
        v1(t) = real(t, wp)/3.0_wp
    end do

    call nc_create("test_ds_get.nc")
    call nc_write_dim("test_ds_get.nc", "xc", x=1.0_dp, dx=1.0_dp, nx=nx)
    call nc_write_dim("test_ds_get.nc", "yc", x=1.0_dp, dx=1.0_dp, nx=ny)
    call nc_write_dim("test_ds_get.nc", "zc", x=1.0_dp, dx=1.0_dp, nx=nz)
    call nc_write_dim("test_ds_get.nc", "time", x=1.0_dp, dx=1.0_dp, nx=nt, unlimited=.false.)
    call nc_write("test_ds_get.nc", "v4", v4, dim1="xc", dim2="yc", dim3="zc", dim4="time")
    call nc_write("test_ds_get.nc", "v3", v3, dim1="xc", dim2="yc", dim3="time", missing_value=mv)
    call nc_write("test_ds_get.nc", "v1", v1, dim1="time")

    ! 3-D + time
    call varslice_init_arg(vs, "test_ds_get.nc", "v4", "1", "1", with_time=.true., &
                           time_par=[1.0_dp, real(nt,dp), 1.0_dp, 0.0_dp])
    call varslice_update(vs, [2.0_dp, 4.0_dp], method="range")
    call varslice_get(vs, d3, k=2)
    call check("3-D dp at k=2", all(d3 .eq. real(v4(:,:,:,3), dp)))
    call varslice_get(vs, s3)
    call check("3-D sp at k=1", all(s3 .eq. v4(:,:,:,2)))
    call varslice_end(vs)

    ! 2-D + time
    call varslice_init_arg(vs, "test_ds_get.nc", "v3", "1", "1", with_time=.true., &
                           time_par=[1.0_dp, real(nt,dp), 1.0_dp, 0.0_dp])
    call varslice_update(vs, [3.0_dp], method="exact")
    call varslice_get(vs, d2)
    call check("2-D dp", all(d2 .eq. real(v3(:,:,3), dp)))
    call check("2-D dp: missing value as mv", d2(2,2) .eq. real(mv, dp))
    call varslice_get(vs, s2)
    call check("2-D sp", all(s2 .eq. v3(:,:,3)))
    call varslice_get(vs, d1, k=1)
    call check("1-D dp (column of the first level)", all(d1 .eq. real(v3(:,1,3), dp)))
    call varslice_end(vs)

    ! 0-D + time
    call varslice_init_arg(vs, "test_ds_get.nc", "v1", "1", "1", with_time=.true., &
                           time_par=[1.0_dp, real(nt,dp), 1.0_dp, 0.0_dp])
    call varslice_update(vs, [1.0_dp, 5.0_dp], method="range")
    call varslice_get(vs, d0, k=4)
    call check("0-D dp at k=4 (exact widening)", d0 .eq. real(v1(4), dp))
    call varslice_get(vs, s0)
    call check("0-D sp at k=1", s0 .eq. v1(1))
    call varslice_end(vs)

    if (nfail .gt. 0) then
        write(*,*) "test_varslice_get: ", nfail, " FAILED"
        error stop 1
    end if
    write(*,*) "test_varslice_get: all passed"

contains

    subroutine check(label, ok)
        character(len=*), intent(in) :: label
        logical,          intent(in) :: ok
        write(*,"(a,a)") merge("  ok   ", "  FAIL ", ok), label
        if (.not. ok) nfail = nfail + 1
    end subroutine check

end program test_varslice_get
