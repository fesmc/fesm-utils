program test_ncio_tiny
    ! nc_read returns exactly what is on file: values below 1e-30 are not
    ! flushed to zero (they were, until the underflow flush was removed from the
    ! reader). Tiny values, written losslessly, are read back bit-identical:
    !   * double: 1e-35, 3e-120, 7e-300, with and without missing_value, 1-D/2-D/3-D;
    !   * float:  1e-35, 1e-38 (subnormal), 1e-30 and 1e-20;
    !   * a slice written and read with start/count;
    !   * zeros, and values of ordinary size, are unchanged beside them.
    ! Not covered: -0.0 on file is read as +0.0 (the identity scale/offset step).
    ! The file is written by the test itself, then removed.

    use ncio

    implicit none

    character(len=*), parameter :: fnm = "test_ncio_tiny.nc"
    integer, parameter :: dp = kind(1.d0), sp = kind(1.0)
    integer, parameter :: i8 = selected_int_kind(18)
    integer :: nfail

    nfail = 0

    write(*,*) "==================================================="
    write(*,*) " ncio tiny-value (no flush on read) test"
    write(*,*) "==================================================="

    call test_double()
    call test_float()
    call test_slice()

    write(*,"(a)") ""
    if (nfail .eq. 0) then
        write(*,*) "All ncio tiny-value tests passed."
    else
        write(*,*) "ncio tiny-value tests FAILED: ", nfail
        stop 1
    end if

contains

    subroutine new_file(nz)
        integer, intent(IN) :: nz
        integer :: stat
        open(unit=99, file=fnm, status="old", iostat=stat)
        if (stat .eq. 0) close(99, status="delete")
        call nc_create(fnm)
        call nc_write_dim(fnm, "x", x=1.d0, dx=1.d0, nx=4)
        call nc_write_dim(fnm, "y", x=1.d0, dx=1.d0, nx=3)
        if (nz .gt. 0) call nc_write_dim(fnm, "z", x=1.d0, dx=1.d0, nx=nz)
    end subroutine new_file

    subroutine delete_file()
        open(unit=99, file=fnm, status="old")
        close(99, status="delete")
    end subroutine delete_file

    subroutine report(label, nbad)
        character(len=*), intent(IN) :: label
        integer, intent(IN) :: nbad
        write(*,"(a,a,i0)") "  ", label, nbad
        if (nbad .gt. 0) nfail = nfail + 1
    end subroutine report

    integer function ndiff_dp(a, b)
        real(dp), intent(IN) :: a(:), b(:)
        integer :: i
        ndiff_dp = 0
        do i = 1, size(a)
            if (transfer(a(i), 1_i8) .ne. transfer(b(i), 1_i8)) ndiff_dp = ndiff_dp + 1
        end do
    end function ndiff_dp

    integer function ndiff_sp(a, b)
        real(sp), intent(IN) :: a(:), b(:)
        integer :: i
        ndiff_sp = 0
        do i = 1, size(a)
            if (transfer(a(i), 1) .ne. transfer(b(i), 1)) ndiff_sp = ndiff_sp + 1
        end do
    end function ndiff_sp

    subroutine test_double()
        real(dp) :: a1(4), a2(4,3), a3(4,3,2), r1(4), r2(4,3), r3(4,3,2), r0, s0
        integer :: k

        a1 = [1.d-35, 3.d-120, 7.d-300, 0.d0]
        a2 = reshape([(a1(mod(k,4)+1) * (1.d0 + 0.125d0*mod(k,3)), k = 0, 11)], [4,3])
        a2(2,2) = 12.5d0
        a2(3,3) = -2.d-40
        a3(:,:,1) = a2
        a3(:,:,2) = 0.d0 - a2          ! not -a2: keeps the zeros +0 (a read gives -0 as +0)
        s0 = 3.d-120

        call new_file(2)
        call nc_write(fnm, "a1", a1, dim1="x")
        call nc_write(fnm, "a2", a2, dim1="x", dim2="y")
        call nc_write(fnm, "a3", a3, dim1="x", dim2="y", dim3="z")
        call nc_write(fnm, "s0", s0)
        call nc_write(fnm, "a2m", a2, dim1="x", dim2="y", missing_value=-9999.d0)

        call nc_read(fnm, "a1", r1)
        call report("double 1-D bit-identical:             bad = ", ndiff_dp(a1, r1))
        call nc_read(fnm, "a2", r2)
        call report("double 2-D bit-identical:             bad = ", ndiff_dp(reshape(a2,[12]), reshape(r2,[12])))
        call nc_read(fnm, "a3", r3)
        call report("double 3-D bit-identical:             bad = ", ndiff_dp(reshape(a3,[24]), reshape(r3,[24])))
        call nc_read(fnm, "s0", r0)
        call report("double scalar bit-identical:          bad = ", ndiff_dp([s0], [r0]))
        r2 = 0.d0
        call nc_read(fnm, "a2m", r2, missing_value=-9999.d0)
        call report("double with missing_value identical:  bad = ", ndiff_dp(reshape(a2,[12]), reshape(r2,[12])))
        call delete_file()
    end subroutine test_double

    subroutine test_float()
        real(sp) :: a1(4), a2(4,3), r1(4), r2(4,3), r0, s0
        integer :: k

        a1 = [1.e-35, 1.e-38, 1.e-30, 1.e-20]
        a2 = reshape([(a1(mod(k,4)+1) * (1.0 + 0.125*mod(k,3)), k = 0, 11)], [4,3])
        a2(2,2) = 12.5
        a2(3,3) = 0.0
        s0 = 1.e-38

        call new_file(0)
        call nc_write(fnm, "a1", a1, dim1="x")
        call nc_write(fnm, "a2", a2, dim1="x", dim2="y")
        call nc_write(fnm, "s0", s0)

        call nc_read(fnm, "a1", r1)
        call report("float  1-D bit-identical:             bad = ", ndiff_sp(a1, r1))
        call nc_read(fnm, "a2", r2)
        call report("float  2-D bit-identical:             bad = ", ndiff_sp(reshape(a2,[12]), reshape(r2,[12])))
        call nc_read(fnm, "s0", r0)
        call report("float  scalar bit-identical:          bad = ", ndiff_sp([s0], [r0]))
        call delete_file()
    end subroutine test_float

    ! Slices of a 3-D variable, written and read one at a time
    subroutine test_slice()
        real(dp) :: a(4,3), r(4,3)
        integer :: n, nbad

        call new_file(2)
        nbad = 0
        do n = 1, 2
            a = 7.d-300 * real(n, dp)
            a(1,1) = 1.d-35
            a(2,2) = 3.d-120 * real(n, dp)
            call nc_write(fnm, "a", a, dim1="x", dim2="y", dim3="z", start=[1,1,n], count=[4,3,1])
        end do
        do n = 1, 2
            a = 7.d-300 * real(n, dp)
            a(1,1) = 1.d-35
            a(2,2) = 3.d-120 * real(n, dp)
            call nc_read(fnm, "a", r, start=[1,1,n], count=[4,3,1])
            nbad = nbad + ndiff_dp(reshape(a,[12]), reshape(r,[12]))
        end do
        call report("double slices bit-identical:          bad = ", nbad)
        call delete_file()
    end subroutine test_slice

end program test_ncio_tiny
