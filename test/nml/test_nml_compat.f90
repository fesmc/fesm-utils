program test_nml_compat
    ! Read every key of a key list (lines `file group key`, from nml_keys) as
    ! a string vector with nml_read and print one line per key:
    !     file group key : n : v1 | v2 | ...
    ! (n = index of the last non-empty element), or `file group key : absent`
    ! when nml_has_param does not find it. Uses only the nml API of v1.4, so
    ! test/nml/nml_compat.sh builds it against a reference commit as well and
    ! diffs the two outputs.

    use nml, only : nml_read, nml_has_param, nml_set_verbose

    implicit none

    integer, parameter :: NMAX = 128
    character(len=2048) :: list, line, fname, group, key
    character(len=256)  :: vec(NMAX)
    character(len=:), allocatable :: out
    integer :: u, ios, n, i

    if (command_argument_count() /= 1) then
        write(*,"(a)") "usage: test_nml_compat.x <key list: lines 'file group key'>"
        error stop 2
    end if
    call get_command_argument(1, list)

    call nml_set_verbose(.false.)

    open(newunit=u, file=trim(list), status="old", action="read")
    do
        read(u, "(a)", iostat=ios) line
        if (ios /= 0) exit
        if (len_trim(line) == 0) cycle
        call split3(line, fname, group, key)
        out = trim(fname)//" "//trim(group)//" "//trim(key)//" : "
        if (.not. nml_has_param(trim(fname), trim(group), trim(key))) then
            write(*,"(a)") out//"absent"
            cycle
        end if
        vec = ""
        call nml_read(trim(fname), trim(group), trim(key), vec)
        n = 0
        do i = 1, NMAX
            if (len_trim(vec(i)) > 0) n = i
        end do
        out = out//itoa(n)//" :"
        do i = 1, n
            if (i > 1) out = out//" |"
            out = out//" "//trim(vec(i))
        end do
        write(*,"(a)") out
    end do
    close(u)

contains

    subroutine split3(s, a, b, c)
        ! The three blank-separated tokens of s (list-directed input would
        ! stop at the '/' of a path).
        character(len=*), intent(in)  :: s
        character(len=*), intent(out) :: a, b, c
        character(len=len(s)) :: r
        integer :: p
        r = adjustl(s)
        p = index(r, " ")
        a = r(1:p-1)
        r = adjustl(r(p:))
        p = index(r, " ")
        b = r(1:p-1)
        c = adjustl(r(p:))
    end subroutine split3

    function itoa(k) result(s)
        integer, intent(in) :: k
        character(len=:), allocatable :: s
        character(len=16) :: buf
        write(buf, "(i0)") k
        s = trim(buf)
    end function itoa

end program test_nml_compat
