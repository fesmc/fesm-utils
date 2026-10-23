program test_ncio_precision
    ! Exercise the precision of nc_write (prec_abs, prec_rel: bit rounding):
    !   * the error of each value read back is within the precision asked;
    !   * lossless without prec_abs/prec_rel;
    !   * missing values, zeros, NaN and Inf are kept (NaN read with nf90,
    !     as nc_read replaces it by its missing value);
    !   * slices written with start/count are rounded each;
    !   * the size of the files, rounded vs lossless;
    !   * nc_set_lossless: no rounding and no precision_* attributes when on;
    !   * nc_write_t (ncio_transpose) passes prec_abs through.
    ! All files are written by the test itself, then removed.

    use ncio
    use ncio_transpose
    use netcdf
    use, intrinsic :: ieee_arithmetic

    implicit none

    integer, parameter :: nx = 300, ny = 300
    double precision, parameter :: mv = -9999.d0
    integer :: nfail

    double precision :: Hd(nx,ny)
    real(4)          :: Hf(nx,ny)

    nfail = 0

    write(*,*) "==================================================="
    write(*,*) " ncio precision test suite"
    write(*,*) "==================================================="
    write(*,"(a)") ""

    call make_field(Hd)
    Hf = real(Hd,4)

    call test_float("lossless",   0.d0,  0.d0)
    call test_float("abs 0.01",   1d-2,  0.d0)
    call test_float("rel 0.1%",   0.d0,  1d-3)
    call test_float("abs or rel", 1d-1,  1d-3)
    call test_double("lossless",  0.d0,  0.d0)
    call test_double("abs 0.01",  1d-2,  0.d0)
    call test_double("rel 1e-6",  0.d0,  1d-6)
    call test_slices()
    call test_nonfinite()
    call test_lossless_switch()
    call test_transpose()

    write(*,"(a)") ""
    if (nfail .eq. 0) then
        write(*,*) "All ncio precision tests passed."
    else
        write(*,*) "ncio precision tests FAILED: ", nfail
        stop 1
    end if

contains

    ! Ice-thickness-like field (m): 0 to ~3000 with small-scale noise and a
    ! few missing values
    subroutine make_field(H)
        double precision, intent(OUT) :: H(:,:)
        integer :: i, j
        double precision :: r
        do j = 1, ny
            do i = 1, nx
                r = sqrt(dble(i-nx/2)**2 + dble(j-ny/2)**2) / dble(nx/2)
                H(i,j) = max(0.d0, 3000.d0*sqrt(max(0.d0, 1.d0-r**2))) &
                         + 3.d0*sin(0.37d0*i)*cos(0.53d0*j)
                if (r .ge. 1.d0) H(i,j) = 0.d0
            end do
        end do
        H(1:5,1) = mv
    end subroutine make_field

    subroutine test_float(label, pabs, prel)
        character(len=*), intent(IN) :: label
        double precision, intent(IN) :: pabs, prel
        character(len=*), parameter :: fnm = "test_ncio_precision.nc"
        real(4) :: R(nx,ny)

        call new_file(fnm)
        if (pabs .gt. 0.d0 .and. prel .gt. 0.d0) then
            call nc_write(fnm, "H", Hf, dim1="x", dim2="y", missing_value=real(mv,4), &
                          prec_abs=real(pabs,4), prec_rel=real(prel,4))
        else if (pabs .gt. 0.d0) then
            call nc_write(fnm, "H", Hf, dim1="x", dim2="y", missing_value=real(mv,4), prec_abs=real(pabs,4))
        else if (prel .gt. 0.d0) then
            call nc_write(fnm, "H", Hf, dim1="x", dim2="y", missing_value=real(mv,4), prec_rel=real(prel,4))
        else
            call nc_write(fnm, "H", Hf, dim1="x", dim2="y", missing_value=real(mv,4))
        end if
        call nc_read(fnm, "H", R)

        ! The precision asked is that of the float values: compare in float
        call check("float  "//label, dble(Hf), dble(R), real(pabs,4)*1.d0, real(prel,4)*1.d0, fnm)
        call delete_file(fnm)
    end subroutine test_float

    subroutine test_double(label, pabs, prel)
        character(len=*), intent(IN) :: label
        double precision, intent(IN) :: pabs, prel
        character(len=*), parameter :: fnm = "test_ncio_precision.nc"
        double precision :: R(nx,ny)

        call new_file(fnm)
        if (pabs .gt. 0.d0 .and. prel .gt. 0.d0) then
            call nc_write(fnm, "H", Hd, dim1="x", dim2="y", missing_value=mv, prec_abs=pabs, prec_rel=prel)
        else if (pabs .gt. 0.d0) then
            call nc_write(fnm, "H", Hd, dim1="x", dim2="y", missing_value=mv, prec_abs=pabs)
        else if (prel .gt. 0.d0) then
            call nc_write(fnm, "H", Hd, dim1="x", dim2="y", missing_value=mv, prec_rel=prel)
        else
            call nc_write(fnm, "H", Hd, dim1="x", dim2="y", missing_value=mv)
        end if
        call nc_read(fnm, "H", R)

        call check("double "//label, Hd, R, pabs, prel, fnm)
        call delete_file(fnm)
    end subroutine test_double

    ! Two time slices of different magnitude, written one at a time
    subroutine test_slices()
        character(len=*), parameter :: fnm = "test_ncio_precision.nc"
        real(4) :: R(nx,ny), F(nx,ny)
        integer :: n

        call new_file(fnm, nt=2)
        do n = 1, 2
            F = Hf * real(10**(n-1),4)
            where (Hf .eq. real(mv,4)) F = real(mv,4)
            call nc_write(fnm, "H", F, dim1="x", dim2="y", dim3="time", start=[1,1,n], count=[nx,ny,1], &
                          missing_value=real(mv,4), prec_abs=0.01)
        end do
        do n = 1, 2
            F = Hf * real(10**(n-1),4)
            where (Hf .eq. real(mv,4)) F = real(mv,4)
            call nc_read(fnm, "H", R, start=[1,1,n], count=[nx,ny,1])
            call check("float  slice "//char(48+n)//"  ", dble(F), dble(R), real(0.01,4)*1.d0, 0.d0, fnm)
        end do
        call delete_file(fnm)
    end subroutine test_slices

    ! NaN and Inf are kept, the values beside them rounded
    subroutine test_nonfinite()
        character(len=*), parameter :: fnm = "test_ncio_precision.nc"
        double precision :: A(4)
        integer :: nbad, ncid, varid, stat

        A = [ieee_value(1.d0, ieee_quiet_nan), ieee_value(1.d0, ieee_positive_inf), &
             ieee_value(1.d0, ieee_negative_inf), 1234.5678d0]
        call nc_create(fnm)
        call nc_write_dim(fnm, "n", x=1.d0, dx=1.d0, nx=4)
        call nc_write(fnm, "A", A, dim1="n", prec_abs=1d-2)

        A = 0.d0
        stat = nf90_open(fnm, nf90_nowrite, ncid)
        stat = nf90_inq_varid(ncid, "A", varid)
        stat = nf90_get_var(ncid, varid, A)
        stat = nf90_close(ncid)
        call delete_file(fnm)

        nbad = 0
        if (.not. ieee_is_nan(A(1))) nbad = nbad + 1
        if (.not. (A(2) .gt. huge(1.d0))) nbad = nbad + 1
        if (.not. (A(3) .lt. -huge(1.d0))) nbad = nbad + 1
        if (abs(A(4) - 1234.5678d0) .gt. 1d-2 .or. A(4) .eq. 1234.5678d0) nbad = nbad + 1
        write(*,"(a,i0)") "  nan/inf kept, finite rounded:  bad = ", nbad
        if (nbad .gt. 0) nfail = nfail + 1
    end subroutine test_nonfinite

    ! nc_set_lossless(.true.): prec_abs ignored, values bit-identical and no
    ! precision_abs attribute; nc_set_lossless(.false.): rounded again
    subroutine test_lossless_switch()
        character(len=*), parameter :: fnm = "test_ncio_precision.nc"
        double precision :: R(nx,ny)
        integer :: nbad

        call nc_set_lossless(.TRUE.)
        call new_file(fnm)
        call nc_write(fnm, "H", Hd, dim1="x", dim2="y", missing_value=mv, prec_abs=1d-2)
        call nc_read(fnm, "H", R)
        nbad = count(R .ne. Hd)
        if (has_att(fnm, "H", "precision_abs")) nbad = nbad + 1
        call delete_file(fnm)
        write(*,"(a,i0)") "  lossless on:  identical, no precision_abs:  bad = ", nbad
        if (nbad .gt. 0) nfail = nfail + 1

        call nc_set_lossless(.FALSE.)
        call new_file(fnm)
        call nc_write(fnm, "H", Hd, dim1="x", dim2="y", missing_value=mv, prec_abs=1d-2)
        call nc_read(fnm, "H", R)
        nbad = 0
        if (all(R .eq. Hd)) nbad = nbad + 1
        if (.not. has_att(fnm, "H", "precision_abs")) nbad = nbad + 1
        write(*,"(a,i0)") "  lossless off: rounded, precision_abs:      bad = ", nbad
        if (nbad .gt. 0) nfail = nfail + 1
        call check("double lossless off", Hd, R, 1d-2, 0.d0, fnm)
        call delete_file(fnm)
    end subroutine test_lossless_switch

    ! nc_write_t writes [j,i] as [i,j], rounded to prec_abs
    subroutine test_transpose()
        character(len=*), parameter :: fnm = "test_ncio_precision.nc"
        double precision :: R(nx,ny), T(ny,nx)

        T = transpose(Hd)
        call new_file(fnm)
        call nc_write_t(fnm, "H", T, dim1="x", dim2="y", missing_value=mv, prec_abs=1d-2)
        call nc_read(fnm, "H", R)
        if (.not. has_att(fnm, "H", "precision_abs")) nfail = nfail + 1
        call check("double transpose  ", Hd, R, 1d-2, 0.d0, fnm)
        call delete_file(fnm)
    end subroutine test_transpose

    logical function has_att(fnm, vnm, anm)
        character(len=*), intent(IN) :: fnm, vnm, anm
        integer :: ncid, varid, stat
        stat = nf90_open(fnm, nf90_nowrite, ncid)
        stat = nf90_inq_varid(ncid, vnm, varid)
        has_att = (nf90_inquire_attribute(ncid, varid, anm) .eq. nf90_noerr)
        stat = nf90_close(ncid)
    end function has_att

    ! Every value within max(pabs, prel*|x|) (exact if both 0); missing and
    ! zero kept
    subroutine check(label, X, R, pabs, prel, fnm)
        character(len=*), intent(IN) :: label, fnm
        double precision, intent(IN) :: X(:,:), R(:,:), pabs, prel
        integer :: i, j, nbad
        integer(8) :: fsize
        double precision :: err, tol, errmax

        nbad = 0
        errmax = 0.d0
        do j = 1, size(X,2)
            do i = 1, size(X,1)
                if (X(i,j) .eq. mv .or. X(i,j) .eq. 0.d0) then
                    if (R(i,j) .ne. X(i,j)) nbad = nbad + 1
                    cycle
                end if
                err = abs(R(i,j) - X(i,j))
                tol = max(pabs, prel*abs(X(i,j)))
                if (err .gt. tol) nbad = nbad + 1
                errmax = max(errmax, err/max(tol, tiny(1.d0)))
            end do
        end do

        inquire(file=fnm, size=fsize)
        write(*,"(a,a,i9,a,f6.3,a,i0)") "  ", label, fsize, " bytes, max err/tol = ", &
            merge(0.d0, errmax, pabs .eq. 0.d0 .and. prel .eq. 0.d0), ", bad = ", nbad
        if (nbad .gt. 0) nfail = nfail + 1
    end subroutine check

    subroutine new_file(fnm, nt)
        character(len=*), intent(IN) :: fnm
        integer, intent(IN), optional :: nt
        call nc_create(fnm)
        call nc_write_dim(fnm, "x", x=1.d0, dx=1.d0, nx=nx)
        call nc_write_dim(fnm, "y", x=1.d0, dx=1.d0, nx=ny)
        if (present(nt)) call nc_write_dim(fnm, "time", x=1.d0, dx=1.d0, nx=nt, unlimited=.TRUE.)
    end subroutine new_file

    subroutine delete_file(fnm)
        character(len=*), intent(IN) :: fnm
        integer :: u
        open(newunit=u, file=fnm, status="old")
        close(u, status="delete")
    end subroutine delete_file

end program test_ncio_precision
