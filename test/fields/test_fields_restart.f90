program test_fields_restart
    ! fesm_fields restart: a state holding a published INSTANT field (masked),
    ! a published MEAN field with a new window already open (rank 3, and a
    ! masked rank 2), an ACCUM window not yet published, an unpublished
    ! INSTANT field and a scalar is written, then read into a fresh state with
    ! the same declaration: every component is bit-identical, and closing the
    ! open windows gives the same result in both states. A changed
    ! declaration (units, shape, extra field) stops: abort/fields_abort_restart_*.
    ! The layout of the file (variables, dimension order, attributes) is
    ! checked with ncio, and a second write over the same file reads back the
    ! same.

    use precision, only : dp
    use ncio, only : nc_dims, nc_exists_var, nc_exists_attr, nc_read_attr, nc_read
    use fesm_fields

    implicit none

    character(len=*), parameter :: fname = "test_fields_restart.nc"
    type(field_state_t) :: a, b, c
    integer  :: h, nfail, i, j, k
    real(dp) :: x2(4,3), x3(4,3,2)
    logical  :: m2(4,3)
    character(len=32), allocatable :: dnames(:)
    character(len=64) :: str
    integer  :: ival
    real(dp) :: dval
    integer  :: vm(4,3)

    nfail = 0

    call declare(a)
    call declare(b)
    call declare(c)

    do j = 1, 3
        do i = 1, 4
            x2(i,j) = sin(real(i + 7*j, dp)) * 1.0e3_dp
            m2(i,j) = mod(i + j, 3) /= 0
            do k = 1, 2
                x3(i,j,k) = cos(real(i*j + 11*k, dp)) / 3.0_dp
            end do
        end do
    end do

    ! Published INSTANT values.
    call fields_put(a, 1, 0.1_dp*acos(-1.0_dp), t1=100.0_dp)
    call fields_put(a, 2, x2, t1=100.0_dp, valid=m2)

    ! MEAN: one published interval [99, 100], then a window [100, 100.25] open.
    call fields_put(a, 3, x3,       t0=99.0_dp,  t1=99.5_dp)
    call fields_put(a, 3, 2.0_dp*x3, t0=99.5_dp, t1=100.0_dp)
    call fields_put(a, 4, x2, t0=99.0_dp, t1=100.0_dp, valid=m2)
    call fields_put(a, 5, x2/7.0_dp, t0=99.0_dp, t1=100.0_dp)
    call fields_finalize_interval(a, 100.0_dp)
    call fields_put(a, 3, x3 + 1.0_dp, t0=100.0_dp, t1=100.25_dp)
    call fields_put(a, 4, x2 - 1.0_dp, t0=100.0_dp, t1=100.25_dp, valid=.not. m2)
    call fields_put(a, 5, x2/3.0_dp,   t0=100.0_dp, t1=100.25_dp)

    call fields_restart_write(a, fname)
    call fields_restart_read(b, fname)

    ! The layout of the file.
    str = ""
    call nc_read_attr(fname, "fesm_fields_owner", str)
    call check("global owner", trim(str) == "test:restart")
    call nc_read_attr(fname, "fesm_fields_n", ival)
    call check("global n", ival == 6)
    str = ""
    call nc_read_attr(fname, "fesm_fields_names", str)
    call check("global names", trim(str) == "bsl,H,tas,smb,melt,T_b")
    call nc_dims(fname, "tas", names=dnames)
    call check("tas dimensions (g1_x, g1_y, month)", size(dnames) == 3 .and. trim(dnames(1)) == "g1_x" .and. &
               trim(dnames(2)) == "g1_y" .and. trim(dnames(3)) == "month")
    call nc_dims(fname, "bsl", names=dnames)
    call check("scalar on the dimension one", size(dnames) == 1 .and. trim(dnames(1)) == "one")
    call check("companions: tas_acc, smb_valid, smb_acc_valid, H_valid", nc_exists_var(fname, "tas_acc") .and. &
               nc_exists_var(fname, "smb_valid") .and. nc_exists_var(fname, "smb_acc_valid") .and. &
               nc_exists_var(fname, "H_valid"))
    call check("no companions for INSTANT unmasked: bsl_acc, H_acc, bsl_valid, melt_valid", &
               .not. (nc_exists_var(fname, "bsl_acc") .or. nc_exists_var(fname, "H_acc") .or. &
                      nc_exists_var(fname, "bsl_valid") .or. nc_exists_var(fname, "melt_valid")))
    call nc_read_attr(fname, "tas", "tsem", ival)
    call check("tas:tsem = MEAN", ival == FIELD_MEAN)
    str = ""
    call nc_read_attr(fname, "tas", "units", str)
    call check("tas:units", trim(str) == "K")
    str = ""
    call nc_read_attr(fname, "tas", "grid", str)
    call check("tas:grid", trim(str) == "g1")
    str = ""
    call nc_read_attr(fname, "tas", "zdim", str)
    call check("tas:zdim", trim(str) == "month")
    call nc_read_attr(fname, "tas", "wsum", dval)
    call check("tas:wsum", dval == 0.25_dp)
    call check("no grid, zdim on a scalar; no zdim on rank 2", .not. (nc_exists_attr(fname, "bsl", "grid") .or. &
               nc_exists_attr(fname, "bsl", "zdim") .or. nc_exists_attr(fname, "H", "zdim")))
    call nc_read(fname, "H_valid", vm)
    call check("H_valid is the mask as 0/1", all(vm == merge(1, 0, m2)))

    do h = 1, a%n
        associate (fa => a%f(h), fb => b%f(h))
            call check(trim(fa%name)//": published, t0, t1", (fa%published .eqv. fb%published) .and. &
                       fa%t0 == fb%t0 .and. fa%t1 == fb%t1)
            call check(trim(fa%name)//": v bit-identical", all(fa%v == fb%v))
            if (allocated(fa%acc)) then
                call check(trim(fa%name)//": acc bit-identical", all(fa%acc == fb%acc))
                call check(trim(fa%name)//": window and wsum", fa%ta0 == fb%ta0 .and. fa%ta1 == fb%ta1 .and. &
                           fa%wsum == fb%wsum)
            end if
            if (allocated(fa%valid)) then
                call check(trim(fa%name)//": valid identical", all(fa%valid .eqv. fb%valid))
            end if
            if (allocated(fa%acc_valid)) then
                call check(trim(fa%name)//": acc_valid identical", all(fa%acc_valid .eqv. fb%acc_valid))
            end if
        end associate
    end do
    call check("open MEAN window survives", b%f(3)%wsum == 0.25_dp .and. b%f(3)%ta0 == 100.0_dp)

    ! A second write over the same file leaves nothing of the first.
    call fields_restart_write(a, fname)
    call fields_restart_read(c, fname)
    call check("second write over the file reads back the same", &
               all(a%f(3)%v == c%f(3)%v) .and. all(a%f(3)%acc == c%f(3)%acc) .and. &
               all(a%f(4)%valid .eqv. c%f(4)%valid) .and. a%f(3)%wsum == c%f(3)%wsum)
    call check("unpublished field stays unpublished", .not. b%f(6)%published)

    ! Closing the open windows gives the same result in both states.
    call fields_put(a, 3, x3, t0=100.25_dp, t1=101.0_dp)
    call fields_put(b, 3, x3, t0=100.25_dp, t1=101.0_dp)
    call fields_put(a, 4, x2, t0=100.25_dp, t1=101.0_dp, valid=m2)
    call fields_put(b, 4, x2, t0=100.25_dp, t1=101.0_dp, valid=m2)
    call fields_put(a, 5, x2, t0=100.25_dp, t1=101.0_dp)
    call fields_put(b, 5, x2, t0=100.25_dp, t1=101.0_dp)
    call fields_finalize_interval(a, 101.0_dp)
    call fields_finalize_interval(b, 101.0_dp)
    call check("continued MEAN identical", all(a%f(3)%v == b%f(3)%v) .and. all(a%f(4)%v == b%f(4)%v))
    call check("continued mask identical", all(a%f(4)%valid .eqv. b%f(4)%valid))
    call check("continued ACCUM identical", all(a%f(5)%v == b%f(5)%v))

    call execute_command_line("rm -f "//fname)

    if (nfail > 0) then
        write(*,*) "test_fields_restart: ", nfail, " check(s) failed."
        stop 1
    end if
    write(*,*) "PASS: test_fields_restart"

contains

    subroutine declare(st)
        type(field_state_t), intent(out) :: st
        integer :: hh

        call fields_init(st, "test:restart")
        call fields_add(st, "bsl",  "m",           FIELD_INSTANT, FIELD_INTENSIVE, hh)
        call fields_add(st, "H",    "m",           FIELD_INSTANT, FIELD_EXTENSIVE, hh, grid="g1", dims=[4,3], &
                        masked=.true.)
        call fields_add(st, "tas",  "K",           FIELD_MEAN,    FIELD_INTENSIVE, hh, grid="g1", dims=[4,3,2], &
                        zdim="month")
        call fields_add(st, "smb",  "kg m-2 yr-1", FIELD_MEAN,    FIELD_EXTENSIVE, hh, grid="g1", dims=[4,3], &
                        masked=.true.)
        call fields_add(st, "melt", "kg m-2",      FIELD_ACCUM,   FIELD_EXTENSIVE, hh, grid="g1", dims=[4,3])
        call fields_add(st, "T_b",  "K",           FIELD_INSTANT, FIELD_INTENSIVE, hh, grid="g2", dims=[2,2])
        call fields_lock(st)

    end subroutine declare

    subroutine check(label, ok)
        character(len=*), intent(in) :: label
        logical,          intent(in) :: ok
        if (ok) then
            write(*,"(a)") "  ok    "//label
        else
            write(*,"(a)") "  FAIL  "//label
            nfail = nfail + 1
        end if
    end subroutine check

end program test_fields_restart
