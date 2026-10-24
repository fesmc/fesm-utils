module grid_signature
    ! Grid identity for the map cache (docs/dev/coupler-global-design.md, 4).
    !
    ! A grid signature is a one-line text record of everything that defines the
    ! cells of a grid: the coordinate system (mtype, units, planet a/f, lon180,
    ! xy_conv, projection parameters), the axis sizes, and the first / last /
    ! sum of the axes and of the cell bounds. It is written into each new map
    ! cache file (global attributes coords_src_sig / coords_dst_sig) and
    ! compared on load.
    !
    ! Signatures are compared with tolerance (1e-9 of a cell for axis values,
    ! 1e-9 relative for the coordinate-system scalars), not as exact strings: a
    ! grid read back from grid_<name>.txt (axes written f0.12) differs from the
    ! in-memory grid in the last bits, and an exact hash would raise false
    ! mismatches.
    !
    ! grid_hash8 is a short tag for file names only: FNV-1a 64 over the
    ! signature with reals rounded to 6 significant digits, first 8 hex
    ! characters.

    use, intrinsic :: iso_fortran_env, only : int64
    use precision,   only : dp
    use coordinates, only : grid_class

    implicit none
    private

    character(len=*), parameter :: SIG_TAG = "coords-sig-1"
    integer,          parameter :: NREAL   = 21      ! number of real fields

    public :: grid_signature_str, grid_signature_match, grid_hash8

contains

    function grid_signature_str(grid) result(sig)
        ! The signature record of a grid.
        type(grid_class), intent(in) :: grid
        character(len=:), allocatable :: sig
        real(dp) :: r(NREAL)
        character(len=32) :: buf
        character(len=1024) :: line
        integer :: i

        call sig_reals(grid, r)
        write(line, "(a,i0,a,i0,a,l1)") SIG_TAG//" mtype="//trim(grid%cs%mtype)// &
            " units="//trim(grid%cs%units)//" nx=", grid%G%nx, " ny=", grid%G%ny, &
            " lon180=", grid%cs%is_lon180
        sig = trim(line)//" r="
        do i = 1, NREAL
            write(buf, "(es24.16e3)") r(i)
            sig = sig//trim(adjustl(buf))
            if (i < NREAL) sig = sig//","
        end do
    end function grid_signature_str

    logical function grid_signature_match(sig1, sig2, why) result(ok)
        ! Do two signature records describe the same grid (with tolerance)?
        ! `why` names the first differing field.
        character(len=*), intent(in) :: sig1, sig2
        character(len=*), intent(out), optional :: why
        character(len=256) :: t1, t2
        integer  :: nx1, ny1, nx2, ny2
        real(dp) :: r1(NREAL), r2(NREAL), cx, cy, tol
        integer  :: i
        logical  :: ok1, ok2

        ok = .false.
        if (present(why)) why = ""
        call sig_parse(sig1, t1, nx1, ny1, r1, ok1)
        call sig_parse(sig2, t2, nx2, ny2, r2, ok2)
        if (.not. (ok1 .and. ok2)) then
            if (present(why)) why = "unreadable signature"
            return
        end if
        if (trim(t1) /= trim(t2)) then
            if (present(why)) why = "coordinate system (mtype/units/lon180)"
            return
        end if
        if (nx1 /= nx2 .or. ny1 /= ny2) then
            if (present(why)) why = "axis sizes"
            return
        end if

        ! coordinate-system scalars: relative tolerance
        do i = 1, 9
            if (abs(r1(i) - r2(i)) > 1.0e-9_dp*max(1.0_dp, abs(r1(i)), abs(r2(i)))) then
                if (present(why)) why = "coordinate-system parameter "//field_name(i)
                return
            end if
        end do

        ! axis values: 1e-9 of a cell (sums: per summed value, plus rounding)
        cx = cell_size(r1(10), r1(11), r1(16), r1(17), nx1)
        cy = cell_size(r1(13), r1(14), r1(19), r1(20), ny1)
        do i = 10, NREAL
            select case (i)
                case (10, 11, 16, 17); tol = 1.0e-9_dp*cx
                case (13, 14, 19, 20); tol = 1.0e-9_dp*cy
                case (12)
                    tol = real(nx1,dp)*(1.0e-9_dp*cx + 4.0_dp*epsilon(1.0_dp)*max(abs(r1(10)), abs(r1(11))))
                case (18)
                    tol = real(nx1+1,dp)*(1.0e-9_dp*cx + 4.0_dp*epsilon(1.0_dp)*max(abs(r1(16)), abs(r1(17))))
                case (15)
                    tol = real(ny1,dp)*(1.0e-9_dp*cy + 4.0_dp*epsilon(1.0_dp)*max(abs(r1(13)), abs(r1(14))))
                case default   ! 21
                    tol = real(ny1+1,dp)*(1.0e-9_dp*cy + 4.0_dp*epsilon(1.0_dp)*max(abs(r1(19)), abs(r1(20))))
            end select
            if (abs(r1(i) - r2(i)) > tol) then
                if (present(why)) why = "axis "//field_name(i)
                return
            end if
        end do
        ok = .true.
    end function grid_signature_match

    function grid_hash8(grid) result(h8)
        ! Short hash tag of a grid for file names: FNV-1a 64 over the signature
        ! with reals rounded to 6 significant digits; the top 32 bits in hex.
        type(grid_class), intent(in) :: grid
        character(len=8) :: h8
        real(dp) :: r(NREAL)
        character(len=:), allocatable :: s
        character(len=32) :: buf
        character(len=512) :: line
        integer(int64) :: hi, lo, p, c
        integer(int64), parameter :: M32 = 4294967295_int64     ! 2^32 - 1
        integer :: i

        call sig_reals(grid, r)
        write(line, "(a,i0,a,i0,a,l1)") trim(grid%cs%mtype)//"|"//trim(grid%cs%units)//"|", &
            grid%G%nx, "|", grid%G%ny, "|", grid%cs%is_lon180
        s = trim(line)
        do i = 1, NREAL
            if (r(i) == 0.0_dp) then
                buf = "0"
            else
                write(buf, "(es13.5e3)") r(i)
            end if
            s = s//"|"//trim(adjustl(buf))
        end do

        ! FNV-1a 64 on (hi, lo) 32-bit halves held in int64 (no overflow):
        ! offset basis 0xcbf29ce484222325; prime 0x100000001b3 = 2^40 + 435.
        hi = 3421674724_int64
        lo = 2216829733_int64
        do i = 1, len(s)
            lo = ieor(lo, int(ichar(s(i:i)), int64))
            p  = lo*435_int64                                   ! < 2^41
            c  = shiftr(p, 32)
            hi = iand(hi*435_int64 + c + iand(shiftl(lo, 8), M32), M32)
            lo = iand(p, M32)
        end do
        write(h8, "(z8.8)") hi
        do i = 1, 8
            if (h8(i:i) >= "A" .and. h8(i:i) <= "F") h8(i:i) = achar(iachar(h8(i:i)) + 32)
        end do
    end function grid_hash8

    ! ----- internals ---------------------------------------------------------------

    subroutine sig_reals(grid, r)
        ! The real fields of the signature, in a fixed order:
        !  1 a  2 f  3 xy_conv  4 lambda  5 phi  6 alpha  7 x_e  8 y_n  9 k0
        ! 10 x(1) 11 x(nx) 12 sum(x)   13 y(1) 14 y(ny) 15 sum(y)
        ! 16 x_bnds(0) 17 x_bnds(nx) 18 sum(x_bnds)
        ! 19 y_bnds(0) 20 y_bnds(ny) 21 sum(y_bnds)
        type(grid_class), intent(in)  :: grid
        real(dp),         intent(out) :: r(NREAL)
        integer :: nx, ny
        nx = grid%G%nx; ny = grid%G%ny
        r(1) = grid%cs%planet%a
        r(2) = grid%cs%planet%f
        r(3) = grid%cs%xy_conv
        r(4) = grid%cs%proj%lambda
        r(5) = grid%cs%proj%phi
        r(6) = grid%cs%proj%alpha
        r(7) = grid%cs%proj%x_e
        r(8) = grid%cs%proj%y_n
        r(9) = grid%cs%proj%k0
        r(10) = grid%G%x(1);  r(11) = grid%G%x(nx); r(12) = sum(grid%G%x)
        r(13) = grid%G%y(1);  r(14) = grid%G%y(ny); r(15) = sum(grid%G%y)
        r(16) = grid%G%x_bnds(0); r(17) = grid%G%x_bnds(nx); r(18) = sum(grid%G%x_bnds)
        r(19) = grid%G%y_bnds(0); r(20) = grid%G%y_bnds(ny); r(21) = sum(grid%G%y_bnds)
    end subroutine sig_reals

    subroutine sig_parse(sig, text, nx, ny, r, ok)
        ! Split a signature record into its text key (tag, mtype, units,
        ! lon180), the axis sizes and the real fields.
        character(len=*), intent(in)  :: sig
        character(len=*), intent(out) :: text
        integer,          intent(out) :: nx, ny
        real(dp),         intent(out) :: r(NREAL)
        logical,          intent(out) :: ok
        integer :: inx, iny, il, ir, ios
        character(len=1) :: l180

        ok = .false.; nx = 0; ny = 0; r = 0.0_dp; text = ""
        if (index(sig, SIG_TAG) /= 1) return
        inx = index(sig, " nx="); iny = index(sig, " ny="); il = index(sig, " lon180=")
        ir  = index(sig, " r=")
        if (inx == 0 .or. iny == 0 .or. il == 0 .or. ir == 0) return
        read(sig(inx+4:iny-1), *, iostat=ios) nx;  if (ios /= 0) return
        read(sig(iny+4:il-1),  *, iostat=ios) ny;  if (ios /= 0) return
        l180 = sig(il+8:il+8)
        read(sig(ir+3:), *, iostat=ios) r;        if (ios /= 0) return
        text = sig(1:inx-1)//" lon180="//l180
        ok = .true.
    end subroutine sig_parse

    real(dp) function cell_size(first, last, b0, bn, n) result(c)
        ! Mean cell size of an axis (from its bounds for a single cell).
        real(dp), intent(in) :: first, last, b0, bn
        integer,  intent(in) :: n
        if (n > 1) then
            c = abs(last - first)/real(n-1, dp)
        else
            c = abs(bn - b0)
        end if
        if (c == 0.0_dp) c = 1.0_dp
    end function cell_size

    function field_name(i) result(nm)
        integer, intent(in) :: i
        character(len=16) :: nm
        character(len=16), parameter :: names(NREAL) = [character(len=16) :: &
            "a", "f", "xy_conv", "lambda", "phi", "alpha", "x_e", "y_n", "k0", &
            "x(1)", "x(nx)", "sum(x)", "y(1)", "y(ny)", "sum(y)", &
            "x_bnds(0)", "x_bnds(nx)", "sum(x_bnds)", "y_bnds(0)", "y_bnds(ny)", "sum(y_bnds)"]
        nm = names(i)
    end function field_name

end module grid_signature
