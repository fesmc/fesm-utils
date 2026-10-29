program test_var_io
    ! Exercise the variable tables of variable_io and var_io_write:
    !   * 5-column (old) and 6-column (prec) tables, all prec syntaxes, and
    !     stray non-row lines between the rows (rows stored contiguously);
    !   * var_io_write of float/double/integer, 1D-3D, over time records,
    !     with the dims/missing_value/grid_mapping overrides;
    !   * precision_abs/precision_rel attributes as in the table, none
    !     with var_io_lossless, and the lossless values bit-identical;
    !   * a bad prec string aborts with a non-zero exit status (checked in a
    !     child process: this program run with --bad-prec).
    ! All files are written by the test itself, then removed.

    use ncio
    use variable_io

    implicit none

    integer, parameter :: nx = 40, ny = 30, nz = 5, ncol = 7
    double precision, parameter :: mv = -9999.d0
    character(len=*), parameter :: ftab5 = "test_var_io_5col.md"
    character(len=*), parameter :: ftab6 = "test_var_io_6col.md"
    character(len=*), parameter :: fbad  = "test_var_io_badprec.md"
    character(len=*), parameter :: fbadlog = "test_var_io_badprec.log"
    character(len=*), parameter :: fnc   = "test_var_io.nc"
    character(len=*), parameter :: fnc_ll = "test_var_io_lossless.nc"
    integer :: nfail
    character(len=256) :: arg

    type(var_io_type), allocatable :: t5(:), t6(:)

    nfail = 0

    ! Child mode: load a table with a bad prec entry (must abort, never return)
    if (command_argument_count() >= 1) then
        call get_command_argument(1, arg)
        if (trim(arg) .eq. "--bad-prec") call child_bad_prec()
    end if

    write(*,*) "==================================================="
    write(*,*) " variable_io test suite"
    write(*,*) "==================================================="
    write(*,"(a)") ""

    call write_tables()

    call load_var_io_table(t5,ftab5)
    call load_var_io_table(t6,ftab6)
    write(*,"(a)") ""

    call test_tables()
    call test_write()
    call test_bad_prec()

    call delete_file(ftab5)
    call delete_file(ftab6)

    write(*,"(a)") ""
    if (nfail .eq. 0) then
        write(*,*) "All variable_io tests passed."
    else
        write(*,*) "variable_io tests FAILED: ", nfail
        stop 1
    end if

contains

    subroutine write_tables()
        integer :: u

        open(newunit=u, file=ftab5, status="replace", action="write")
        write(u,"(a)") "# old table, 5 columns"
        write(u,"(a)") ""
        write(u,"(a)") "| id | variable | dimensions | units | long_name     |"
        write(u,"(a)") "|----|----------|------------|-------|---------------|"
        write(u,"(a)") "|  1 | H_ice    | xc, yc     | m     | Ice thickness |"
        write(u,"(a)") ""
        write(u,"(a)") "A stray note between the rows."
        write(u,"(a)") "|  2 | col      | column     | m     | Column field  |"
        close(u)

        open(newunit=u, file=ftab6, status="replace", action="write")
        write(u,"(a)") "# new table, 6 columns"
        write(u,"(a)") ""
        write(u,"(a)") "| id | variable | dimensions   | units | long_name          | prec       |"
        write(u,"(a)") "|----|----------|--------------|-------|--------------------|------------|"
        write(u,"(a)") "|  1 | H_ice    | xc, yc       | m     | Ice thickness      | 1e-2       |"
        write(u,"(a)") "|  2 | uxy_s    | xc, yc       | m/yr  | Surface velocity   | 1e-3r      |"
        write(u,"(a)") "## A stray heading"
        write(u,"(a)") "|  3 | T_ice    | xc, yc, zeta | K     | Ice temperature    | 1e-2,1e-4r |"
        write(u,"(a)") "|  4 | mask_bed | xc, yc       | -     | Bed mask           | -          |"
        write(u,"(a)") ""
        write(u,"(a)") "|  5 | smb_col  | column       | m/yr  | Column mass balance| 1e-4r, 1e-3|"
        close(u)
    end subroutine write_tables

    subroutine child_bad_prec()
        ! Not reached past load_var_io_table when the library aborts as it should.
        type(var_io_type), allocatable :: t(:)
        integer :: u
        open(newunit=u, file=fbad, status="replace", action="write")
        write(u,"(a)") "| id | variable | dimensions | units | long_name     | prec  |"
        write(u,"(a)") "|----|----------|------------|-------|---------------|-------|"
        write(u,"(a)") "|  1 | H_ice    | xc, yc     | m     | Ice thickness | 1e-3x |"
        close(u)
        call load_var_io_table(t, fbad)
        write(*,*) "child_bad_prec: the table was loaded (no abort)"
        stop 0
    end subroutine child_bad_prec

    subroutine test_bad_prec()
        ! A bad prec string must stop the program with a non-zero exit status
        ! and say what is wrong (error_unit, collected with the output).
        character(len=512) :: self, line
        integer :: es, cs, u, ios, nbad
        logical :: msg

        nbad = 0
        call get_command_argument(0, self)
        call execute_command_line('"'//trim(self)//'" --bad-prec > '//fbadlog//' 2>&1', &
                                  exitstat=es, cmdstat=cs)
        if (cs .ne. 0 .or. es .eq. 0) nbad = nbad + 1

        msg = .false.
        open(newunit=u, file=fbadlog, status="old", action="read", iostat=ios)
        if (ios .eq. 0) then
            do
                read(u,"(a)",iostat=ios) line
                if (ios .ne. 0) exit
                if (index(line, "bad prec entry") .gt. 0) msg = .true.
            end do
            close(u)
        end if
        if (.not. msg) nbad = nbad + 1

        call delete_file(fbad)
        call delete_file(fbadlog)
        write(*,"(a,i0,a,i0,a)") "  bad prec string aborts:        bad = ", nbad, "  (child exit status ", es, ")"
        nfail = nfail + nbad
    end subroutine test_bad_prec

    subroutine test_tables()
        type(var_io_type) :: vl
        integer :: nbad

        nbad = 0

        ! 5 columns: rows contiguous across the stray lines, lossless
        if (size(t5) .ne. 2) nbad = nbad + 1
        if (size(t5) .eq. 2) then
            if (trim(t5(1)%varname) .ne. "H_ice" .or. trim(t5(2)%varname) .ne. "col") nbad = nbad + 1
            if (t5(2)%id .ne. 2 .or. t5(2)%ndims .ne. 1 .or. trim(t5(2)%dims(1)) .ne. "column") nbad = nbad + 1
            if (any([t5%prec_abs, t5%prec_rel] .ne. 0.d0)) nbad = nbad + 1
        end if
        write(*,"(a,i0)") "  5-column table:                bad = ", nbad
        if (nbad .gt. 0) nfail = nfail + 1

        ! 6 columns: all prec syntaxes
        nbad = 0
        if (size(t6) .ne. 5) nbad = nbad + 1
        if (size(t6) .eq. 5) then
            if (trim(t6(3)%varname) .ne. "T_ice" .or. t6(3)%ndims .ne. 3) nbad = nbad + 1
            if (trim(t6(3)%long_name) .ne. "Ice temperature") nbad = nbad + 1
            call check_prec(t6(1), 1d-2, 0.d0, nbad)
            call check_prec(t6(2), 0.d0, 1d-3, nbad)
            call check_prec(t6(3), 1d-2, 1d-4, nbad)
            call check_prec(t6(4), 0.d0, 0.d0, nbad)
            call check_prec(t6(5), 1d-3, 1d-4, nbad)
            vl = var_io_lossless(t6(3))
            call check_prec(vl, 0.d0, 0.d0, nbad)
            if (trim(vl%varname) .ne. "T_ice" .or. vl%ndims .ne. 3) nbad = nbad + 1
        end if
        write(*,"(a,i0)") "  6-column table, prec syntaxes: bad = ", nbad
        if (nbad .gt. 0) nfail = nfail + 1

    end subroutine test_tables

    subroutine check_prec(v, pabs, prel, nbad)
        type(var_io_type), intent(IN) :: v
        double precision,  intent(IN) :: pabs, prel
        integer, intent(INOUT) :: nbad
        if (v%prec_abs .ne. pabs .or. v%prec_rel .ne. prel) then
            write(*,*) "    wrong prec: ", trim(v%varname), v%prec_abs, v%prec_rel
            nbad = nbad + 1
        end if
    end subroutine check_prec

    subroutine test_write()
        double precision :: H(nx,ny), H1(nx,ny), Hll(nx,ny)
        real(4)  :: U(nx,ny), T(nx,ny,nz), T1(nx,ny,nz)
        integer  :: M(nx,ny), M1(nx,ny)
        real(4)  :: C(ncol), C1(ncol)
        integer  :: i, j, k, n, nbad
        double precision :: tol

        do j = 1, ny
            do i = 1, nx
                H(i,j) = 1000.d0 + 500.d0*sin(0.31d0*i)*cos(0.17d0*j) + 1d-5*i
                U(i,j) = real(10.d0**(4.d0*dble(i)/nx - 1.d0) * (1.d0 + 0.3d0*sin(0.7d0*j)), 4)
                M(i,j) = mod(i+j,5)
                do k = 1, nz
                    T(i,j,k) = real(250.d0 + 0.37d0*i - 0.11d0*j + 3.3d0*k, 4)
                end do
            end do
        end do
        do i = 1, ncol
            C(i) = real(0.123456789d0*i, 4)
        end do
        C(ncol) = real(mv,4)

        ! Rounded, as in the table, over two time records
        call new_file(fnc)
        do n = 1, 2
            call var_io_write(fnc, t6(1), H + dble(n), n)
            call var_io_write(fnc, t6(2), U, n)
            call var_io_write(fnc, t6(3), T, n)
            call var_io_write(fnc, t6(4), M, n)
            call var_io_write(fnc, t6(5), C, n, dims=["column","time  "], &
                              missing_value=real(mv,4), grid_mapping="")
        end do

        ! Lossless (var_io_lossless) of the same table variables
        call new_file(fnc_ll)
        call var_io_write(fnc_ll, var_io_lossless(t6(1)), H, 1)
        call var_io_write(fnc_ll, var_io_lossless(t6(3)), T, 1)

        ! Values within the precision of the table, written per record
        nbad = 0
        do n = 1, 2
            call nc_read(fnc, "H_ice", H1, start=[1,1,n], count=[nx,ny,1])
            if (maxval(abs(H1 - (H + dble(n)))) .gt. 1d-2) nbad = nbad + 1
            if (all(H1 .eq. H + dble(n))) nbad = nbad + 1                  ! not rounded
        end do
        call nc_read(fnc, "uxy_s", H1, start=[1,1,2], count=[nx,ny,1])
        if (any(abs(H1 - U) .gt. 1d-3*abs(U))) nbad = nbad + 1
        call nc_read(fnc, "T_ice", T1, start=[1,1,1,2], count=[nx,ny,nz,1])
        tol = 1d-2
        if (any(abs(T1 - T) .gt. max(tol, 1d-4*abs(T)))) nbad = nbad + 1
        call nc_read(fnc, "mask_bed", M1, start=[1,1,2], count=[nx,ny,1])
        if (any(M1 .ne. M)) nbad = nbad + 1
        call nc_read(fnc, "smb_col", C1, start=[1,2], count=[ncol,1])
        if (C1(ncol) .ne. real(mv,4)) nbad = nbad + 1
        if (any(abs(C1(1:ncol-1) - C(1:ncol-1)) .gt. max(1d-3, 1d-4*abs(C(1:ncol-1))))) nbad = nbad + 1
        if (nc_size(fnc, "time") .ne. 2) nbad = nbad + 1
        write(*,"(a,i0)") "  rounded values within prec:    bad = ", nbad
        if (nbad .gt. 0) nfail = nfail + 1

        ! Attributes: as in the table, none for lossless and integer data
        nbad = 0
        call check_attr(fnc, "H_ice", 1d-2, 0.d0, nbad)
        call check_attr(fnc, "uxy_s", 0.d0, 1d-3, nbad)
        call check_attr(fnc, "T_ice", 1d-2, 1d-4, nbad)
        call check_attr(fnc, "mask_bed", 0.d0, 0.d0, nbad)
        call check_attr(fnc, "smb_col", 1d-3, 1d-4, nbad)
        call check_attr(fnc_ll, "H_ice", 0.d0, 0.d0, nbad)
        call check_attr(fnc_ll, "T_ice", 0.d0, 0.d0, nbad)
        write(*,"(a,i0)") "  precision attributes:          bad = ", nbad
        if (nbad .gt. 0) nfail = nfail + 1

        ! Lossless: bit-identical
        nbad = 0
        call nc_read(fnc_ll, "H_ice", Hll, start=[1,1,1], count=[nx,ny,1])
        if (any(Hll .ne. H)) nbad = nbad + 1
        call nc_read(fnc_ll, "T_ice", T1, start=[1,1,1,1], count=[nx,ny,nz,1])
        if (any(T1 .ne. T)) nbad = nbad + 1
        write(*,"(a,i0)") "  lossless bit-identical:        bad = ", nbad
        if (nbad .gt. 0) nfail = nfail + 1

        call delete_file(fnc)
        call delete_file(fnc_ll)

    end subroutine test_write

    ! precision_abs/precision_rel present with the value given (float
    ! precision: the value written is rounded to real(4)), absent if 0
    subroutine check_attr(fnm, varname, pabs, prel, nbad)
        character(len=*), intent(IN) :: fnm, varname
        double precision, intent(IN) :: pabs, prel
        integer, intent(INOUT) :: nbad
        call check_attr_one(fnm, varname, "precision_abs", pabs, nbad)
        call check_attr_one(fnm, varname, "precision_rel", prel, nbad)
    end subroutine check_attr

    subroutine check_attr_one(fnm, varname, attname, val, nbad)
        character(len=*), intent(IN) :: fnm, varname, attname
        double precision, intent(IN) :: val
        integer, intent(INOUT) :: nbad
        double precision :: a
        if (val .eq. 0.d0) then
            if (nc_exists_attr(fnm, varname, attname)) then
                write(*,*) "    unexpected attribute: ", trim(varname), ":", trim(attname)
                nbad = nbad + 1
            end if
        else
            if (.not. nc_exists_attr(fnm, varname, attname)) then
                write(*,*) "    missing attribute: ", trim(varname), ":", trim(attname)
                nbad = nbad + 1
            else
                call nc_read_attr(fnm, varname, attname, a)
                if (abs(a - val) .gt. 1d-6*val) then
                    write(*,*) "    wrong attribute: ", trim(varname), ":", trim(attname), a
                    nbad = nbad + 1
                end if
            end if
        end if
    end subroutine check_attr_one

    subroutine new_file(fnm)
        character(len=*), intent(IN) :: fnm
        call nc_create(fnm)
        call nc_write_dim(fnm, "xc", x=1.d0, dx=1.d0, nx=nx)
        call nc_write_dim(fnm, "yc", x=1.d0, dx=1.d0, nx=ny)
        call nc_write_dim(fnm, "zeta", x=1.d0, dx=1.d0, nx=nz)
        call nc_write_dim(fnm, "column", x=1.d0, dx=1.d0, nx=ncol)
        call nc_write_dim(fnm, "time", x=1.d0, dx=1.d0, nx=1, unlimited=.TRUE.)
    end subroutine new_file

    subroutine delete_file(fnm)
        character(len=*), intent(IN) :: fnm
        integer :: u
        open(newunit=u, file=fnm, status="old")
        close(u, status="delete")
    end subroutine delete_file

end program test_var_io
