
module nml 

    use, intrinsic :: iso_fortran_env, only : input_unit, output_unit, error_unit

    implicit none 

    logical :: VERBOSE = .TRUE.            !!  should freshly read namelist be printed to screen?
    logical :: ERROR_NO_PARAM = .TRUE.      !! Should error be thrown if parameter isn't found? 

    integer, parameter :: io_unit_err = error_unit

    ! Maximum number of vector elements to echo in nml_print_*_vector.
    ! Longer vectors are truncated with a "... (n total)" note to avoid
    ! building huge strings that overflow the fixed-length print line.
    integer, parameter :: nml_print_vector_max = 50

    interface nml_read 
        module procedure nml_read_string, nml_read_double, nml_read_float 
        module procedure nml_read_integer, nml_read_logical
        module procedure nml_read_string_vector, nml_read_double_vector
        module procedure nml_read_float_vector, nml_read_integer_vector
        module procedure nml_read_logical_vector
    end interface 

    interface nml_write
        module procedure nml_write_string 
    end interface 

    interface nml_print 
        module procedure nml_print_string, nml_print_double, nml_print_float  
        module procedure nml_print_integer, nml_print_logical
        module procedure nml_print_string_vector, nml_print_double_vector
        module procedure nml_print_float_vector, nml_print_integer_vector
        module procedure nml_print_logical_vector
    end interface 

    private
    public :: nml_read, nml_write
    public :: nml_print
    public :: nml_set_verbose
    public :: nml_replace
    public :: nml_validate

contains

    subroutine nml_set_verbose(verb)

        implicit none 

        logical :: verb 

        VERBOSE = verb 

        write(*,*) "nml_set_verbose:: ", VERBOSE 

        return 

    end subroutine nml_set_verbose

    subroutine nml_replace(s,text,rep,outs)
        ! Adapted from FUNCTION Replace_Text:
        ! http://fortranwiki.org/fortran/show/String_Functions
        CHARACTER(len=*)           :: s,text,rep
        CHARACTER(len=*), optional :: outs
        INTEGER :: i, nt, nr

        character(len=LEN(s)+100) :: tmps  ! Temp string to hold output 

        tmps = s ; nt = LEN_TRIM(text) ; nr = LEN_TRIM(rep)
        DO
           i = INDEX(tmps,text(:nt)) ; IF (i == 0) EXIT
           tmps = tmps(:i-1) // rep(:nr) // tmps(i+nt:)
        END DO

        if (present(outs)) then
            outs = trim(tmps)
        else
            s = trim(tmps)
        end if

        return

    end subroutine nml_replace

    ! =============================================================
    !
    ! nml parameter reading functions
    !
    ! =============================================================

    ! This is the basic nml reading subroutine.
    ! All interfaces use this to read the parameter, then it
    ! is converted to the correct type.
    !
    ! When defaults_file is present, the parameter MUST exist in the
    ! defaults file (the schema). The value is taken from defaults_file
    ! and then overlaid by filename when present. defaults_group lets the
    ! caller request a different group name in defaults (used for yelmo's
    ! group aliasing, e.g. user group "ydyn_north" -> defaults "ydyn").
    subroutine nml_read_internal(filename,group,name,value,comment, &
                                 defaults_file,defaults_group)

        implicit none

        character(len=*), intent(INOUT) :: value
        character(len=*), intent(IN)    :: filename, group, name
        character(len=*), intent(INOUT), optional :: comment
        character(len=*), intent(IN),    optional :: defaults_file
        character(len=*), intent(IN),    optional :: defaults_group

        character(len=256)  :: def_group_actual
        character(len=1024) :: value_def, value_user
        character(len=1024) :: comment_def, comment_user
        logical :: found_def, found_user

        comment_def  = ""
        comment_user = ""

        if (present(defaults_file)) then
            def_group_actual = group
            if (present(defaults_group)) def_group_actual = trim(defaults_group)

            ! Read default value (must exist)
            call nml_find_param(defaults_file,trim(def_group_actual),name, &
                                value_def,found_def,comment_def)
            if (.not. found_def) then
                write(io_unit_err,*) ""
                write(io_unit_err,*) "nml:: Error: parameter not declared in defaults file."
                write(io_unit_err,*) "  Defaults file: ", trim(defaults_file)
                write(io_unit_err,*) "  Group (def):   ", trim(def_group_actual)
                write(io_unit_err,*) "  Parameter:     ", trim(name)
                write(io_unit_err,*) "  (Requested via user file: ", trim(filename), &
                                     ", group: ", trim(group), ")"
                error stop "Program stopped."
            end if

            value = trim(value_def)
            if (present(comment)) comment = trim(comment_def)

            ! Overlay with value from user file (optional)
            call nml_find_param(filename,group,name,value_user,found_user,comment_user)
            if (found_user) then
                value = trim(value_user)
                if (present(comment) .and. trim(comment_user) /= "") &
                    comment = trim(comment_user)
            end if

        else
            ! Legacy path: scan only the user file.
            call nml_find_param(filename,group,name,value_user,found_user,comment_user)
            if (found_user) then
                value = trim(value_user)
                if (present(comment)) comment = trim(comment_user)
            else if (ERROR_NO_PARAM) then
                write(io_unit_err,*) ""
                write(io_unit_err,*) "nml:: Error: parameter not found."
                write(io_unit_err,*) "Filename:  ", trim(filename)
                write(io_unit_err,*) "Group:     ", trim(group)
                write(io_unit_err,*) "Parameter: ", trim(name)
                error stop "Program stopped."
            end if
        end if

        return

    end subroutine nml_read_internal

    ! Scan filename for the first (target_group, target_name) match, returning
    ! its value (and comment). Supports both classic &group ... / blocks and
    ! flat 'group.param = value' lines (single- or multi-parameter per line).
    subroutine nml_find_param(filename,target_group,target_name,value,found,comment)

        implicit none

        character(len=*), intent(IN)  :: filename, target_group, target_name
        character(len=*), intent(OUT) :: value
        logical,          intent(OUT) :: found
        character(len=*), intent(OUT), optional :: comment

        integer :: io, file_unit
        integer :: iostat, l, ltype
        character(len=1000) :: line, name1, value1, comment1
        character(len=256)  :: group_now
        logical :: ingroup, is_flat
        integer :: i, n_flat
        integer, parameter :: MAX_FLAT = 32
        character(len=128) :: flat_groups(MAX_FLAT), flat_names(MAX_FLAT)
        character(len=1024) :: flat_values(MAX_FLAT)

        value = ""
        if (present(comment)) comment = ""
        found = .false.

        ! Open the file, or get unit if already open
        inquire(file=trim(filename),NUMBER=file_unit)
        if (file_unit .eq. -1) then
            open(unit=newunit(io),file=trim(filename),status="old",iostat=iostat)
            if (iostat /= 0) then
                write(io_unit_err,*) "nml:: namelist file could not be opened: "//trim(filename)
                error stop
            end if
        else
            io = file_unit
            rewind(io)
        end if

        ingroup   = .false.
        group_now = ""

        do l = 1, 5000
            read(io,"(a1000)",iostat=iostat) line
            if (iostat /= 0) exit

            ! Detect flat-format line first
            call parse_flat_line(line, is_flat, flat_groups, flat_names, flat_values, n_flat)
            if (is_flat) then
                do i = 1, n_flat
                    if (trim(flat_groups(i)) == trim(target_group) .and. &
                        trim(flat_names(i))  == trim(target_name)) then
                        value = trim(flat_values(i))
                        if (present(comment)) comment = ""
                        found = .true.
                        exit
                    end if
                end do
                if (found) exit
                cycle
            end if

            ! Classic-format parse
            call parse_line(line, ltype, name1, value1, comment1)

            if (ingroup .and. ltype == 3) then
                if (trim(group_now) == trim(target_group) .and. &
                    trim(name1)     == trim(target_name)) then
                    value = trim(value1)
                    if (present(comment)) comment = trim(comment1)
                    found = .true.
                    exit
                end if
            end if

            if (ltype == 1) then
                ingroup   = .true.
                group_now = trim(name1)
            end if
            if (ltype == 2) then
                ingroup   = .false.
                group_now = ""
            end if

            if (l .eq. 5000) then
                write(*,*) "nml:: Warning: maximum nml length of 5000 lines reached."
            end if
        end do

        if (file_unit .eq. -1) close(io)

        return

    end subroutine nml_find_param

    ! Detect whether a line is in flat format (group.param=value ...) and if
    ! so, extract all (group, name, value) tokens it contains. A line is flat
    ! when it does not start with & or /, contains an '=' sign outside of any
    ! quoted string, and the substring before the first such '=' contains a
    ! '.' outside of quotes.
    !
    ! Quote-aware: any '=', '!', whitespace, or '.' that lies inside a
    ! "..." or '...' span is ignored for tokenisation. This allows
    ! string values like  g.method = "polar mode"  to round-trip cleanly,
    ! including when multiple parameters share one line.
    subroutine parse_flat_line(line, is_flat, groups, names, values, n)

        implicit none

        character(len=*), intent(IN)  :: line
        logical,          intent(OUT) :: is_flat
        character(len=*), intent(OUT) :: groups(:), names(:)
        character(len=*), intent(OUT) :: values(:)
        integer,          intent(OUT) :: n

        character(len=len(line)) :: work
        character(len=len(line)) :: key, val
        logical :: in_q(len(line))
        integer :: q_eq, q_dot
        integer :: p, len_w, j, q_space_back
        integer :: ndot, i, k

        is_flat = .false.
        n = 0
        groups(:) = ""
        names(:)  = ""
        values(:) = ""

        work = line
        call build_quote_mask(work, len(work), in_q)

        ! Quote-aware comment strip
        len_w = len_trim(work)
        do i = 1, len_w
            if (work(i:i) == "!" .and. .not. in_q(i)) then
                work(i:) = " "
                exit
            end if
        end do
        len_w = len_trim(work)
        if (len_w == 0) return

        ! Skip leading whitespace
        p = 1
        do while (p <= len_w)
            if (work(p:p) /= ' ' .and. ichar(work(p:p)) /= 9) exit
            p = p + 1
        end do
        if (p > len_w) return
        if (work(p:p) == "&" .or. work(p:p) == "/") return

        ! Must contain unquoted '=' with an unquoted '.' before it.
        q_eq = find_unquoted_char(work, in_q, "=", p, len_w)
        if (q_eq == 0) return
        q_dot = find_unquoted_char(work, in_q, ".", p, q_eq - 1)
        if (q_dot == 0) return

        is_flat = .true.

        do while (p <= len_w .and. n < size(groups))
            q_eq = find_unquoted_char(work, in_q, "=", p, len_w)
            if (q_eq == 0) exit
            key = trim(adjustl(work(p:q_eq-1)))

            ! Validate key: exactly one '.', no quote chars
            ndot = 0
            do k = 1, len_trim(key)
                if (key(k:k) == '.') ndot = ndot + 1
            end do
            if (ndot /= 1) exit

            ! Where does this value end?
            j = find_unquoted_char(work, in_q, "=", q_eq + 1, len_w)
            if (j == 0) then
                val = trim(adjustl(work(q_eq+1:len_w)))
                n = n + 1
                call split_group_dot_name(key, groups(n), names(n))
                values(n) = val
                call remove_quotes_comma(values(n))
                return
            end if

            ! The next key is the last unquoted-whitespace-separated token
            ! before j. Walk left from j-1, skip trailing unquoted whitespace,
            ! then continue past the next-key token until the unquoted
            ! whitespace that bounds it.
            q_dot = j - 1
            do while (q_dot >= q_eq + 1)
                if (in_q(q_dot)) exit
                if (work(q_dot:q_dot) /= ' ' .and. ichar(work(q_dot:q_dot)) /= 9) exit
                q_dot = q_dot - 1
            end do
            if (q_dot < q_eq + 1) exit  ! no key between two '='

            q_space_back = q_eq
            do i = q_dot - 1, q_eq + 1, -1
                if (in_q(i)) cycle
                if (work(i:i) == ' ' .or. ichar(work(i:i)) == 9) then
                    q_space_back = i
                    exit
                end if
            end do

            val = trim(adjustl(work(q_eq+1:q_space_back)))
            n = n + 1
            call split_group_dot_name(key, groups(n), names(n))
            values(n) = val
            call remove_quotes_comma(values(n))

            p = q_space_back + 1
            do while (p <= len_w)
                if (in_q(p)) exit
                if (work(p:p) /= ' ' .and. ichar(work(p:p)) /= 9) exit
                p = p + 1
            end do
        end do

        return

    end subroutine parse_flat_line

    ! Fill in_q(i)=.TRUE. for every position inside a "..." or '...' span,
    ! including the quote characters themselves. Quotes do not nest; an
    ! opening " ignores subsequent ' until the matching " (and vice versa).
    subroutine build_quote_mask(s, len_s, in_q)

        implicit none

        character(len=*), intent(IN)  :: s
        integer,          intent(IN)  :: len_s
        logical,          intent(OUT) :: in_q(:)

        integer :: i
        logical :: in_dq, in_sq

        in_q(:) = .false.
        in_dq = .false.
        in_sq = .false.

        do i = 1, len_s
            if (in_dq) then
                in_q(i) = .true.
                if (s(i:i) == '"') in_dq = .false.
            else if (in_sq) then
                in_q(i) = .true.
                if (s(i:i) == "'") in_sq = .false.
            else
                if (s(i:i) == '"') then
                    in_q(i) = .true.
                    in_dq = .true.
                else if (s(i:i) == "'") then
                    in_q(i) = .true.
                    in_sq = .true.
                end if
            end if
        end do

        return

    end subroutine build_quote_mask

    ! Find the first occurrence of character c in s(start_pos:end_pos)
    ! whose position is not inside a quoted span (per in_q). Returns 0 if
    ! not found.
    function find_unquoted_char(s, in_q, c, start_pos, end_pos) result(pos)

        implicit none

        character(len=*), intent(IN) :: s
        logical,          intent(IN) :: in_q(:)
        character(len=1), intent(IN) :: c
        integer,          intent(IN) :: start_pos, end_pos
        integer :: pos

        integer :: i

        pos = 0
        do i = start_pos, end_pos
            if (.not. in_q(i) .and. s(i:i) == c) then
                pos = i
                return
            end if
        end do

        return

    end function find_unquoted_char

    subroutine count_char(s, c, n)
        character(len=*), intent(IN) :: s
        character(len=1), intent(IN) :: c
        integer,          intent(OUT) :: n
        integer :: i
        n = 0
        do i = 1, len_trim(s)
            if (s(i:i) == c) n = n + 1
        end do
        return
    end subroutine count_char

    subroutine split_group_dot_name(key, group_out, name_out)
        character(len=*), intent(IN)  :: key
        character(len=*), intent(OUT) :: group_out, name_out
        integer :: q
        q = index(key, ".")
        group_out = trim(adjustl(key(1:q-1)))
        name_out  = trim(adjustl(key(q+1:)))
        return
    end subroutine split_group_dot_name

    ! Walk filename collecting all (group, name) pairs whose group matches
    ! target_group (or all if target_group is empty). Used by nml_validate.
    subroutine nml_collect_group_params(filename, target_group, names, n)

        implicit none

        character(len=*), intent(IN)  :: filename, target_group
        character(len=*), intent(OUT) :: names(:)
        integer,          intent(OUT) :: n

        integer :: io, file_unit, iostat, l, ltype
        character(len=1000) :: line, name1, value1, comment1
        character(len=256)  :: group_now
        logical :: ingroup, is_flat
        integer :: i, n_flat
        integer, parameter :: MAX_FLAT = 32
        character(len=128) :: flat_groups(MAX_FLAT), flat_names(MAX_FLAT)
        character(len=1024) :: flat_values(MAX_FLAT)

        n = 0
        names(:) = ""

        inquire(file=trim(filename),NUMBER=file_unit)
        if (file_unit .eq. -1) then
            open(unit=newunit(io),file=trim(filename),status="old",iostat=iostat)
            if (iostat /= 0) then
                write(io_unit_err,*) "nml:: namelist file could not be opened: "//trim(filename)
                error stop
            end if
        else
            io = file_unit
            rewind(io)
        end if

        ingroup   = .false.
        group_now = ""

        do l = 1, 5000
            read(io,"(a1000)",iostat=iostat) line
            if (iostat /= 0) exit

            call parse_flat_line(line, is_flat, flat_groups, flat_names, flat_values, n_flat)
            if (is_flat) then
                do i = 1, n_flat
                    if (trim(flat_groups(i)) == trim(target_group)) then
                        if (n < size(names)) then
                            n = n + 1
                            names(n) = trim(flat_names(i))
                        end if
                    end if
                end do
                cycle
            end if

            call parse_line(line, ltype, name1, value1, comment1)

            if (ingroup .and. ltype == 3 .and. trim(group_now) == trim(target_group)) then
                if (n < size(names)) then
                    n = n + 1
                    names(n) = trim(name1)
                end if
            end if

            if (ltype == 1) then
                ingroup   = .true.
                group_now = trim(name1)
            end if
            if (ltype == 2) then
                ingroup   = .false.
                group_now = ""
            end if
        end do

        if (file_unit .eq. -1) close(io)

        return

    end subroutine nml_collect_group_params

    ! Validate that every parameter set under group in filename is also
    ! declared (under defaults_group or, if not present, group) in
    ! defaults_file. Reports each unknown parameter; halts at the end if
    ! halt_on_error (default .TRUE.) and any were found.
    subroutine nml_validate(filename, defaults_file, group, defaults_group, halt_on_error)

        implicit none

        character(len=*), intent(IN) :: filename, defaults_file, group
        character(len=*), intent(IN), optional :: defaults_group
        logical,          intent(IN), optional :: halt_on_error

        integer, parameter :: MAX_PARAMS = 2000
        character(len=128) :: names(MAX_PARAMS)
        character(len=1024) :: dummy_value
        character(len=256) :: def_group_actual
        integer :: n_user, i
        logical :: halt, any_err, found

        halt = .true.
        if (present(halt_on_error)) halt = halt_on_error

        def_group_actual = group
        if (present(defaults_group)) def_group_actual = trim(defaults_group)

        call nml_collect_group_params(filename, group, names, n_user)

        any_err = .false.
        do i = 1, n_user
            call nml_find_param(defaults_file, trim(def_group_actual), trim(names(i)), &
                                dummy_value, found)
            if (.not. found) then
                if (.not. any_err) then
                    write(io_unit_err,*) ""
                    write(io_unit_err,*) "nml_validate:: Unknown parameter(s) in user file."
                    write(io_unit_err,*) "  user file:    ", trim(filename)
                    write(io_unit_err,*) "  defaults:     ", trim(defaults_file)
                    write(io_unit_err,*) "  user group:   ", trim(group)
                    write(io_unit_err,*) "  def  group:   ", trim(def_group_actual)
                end if
                write(io_unit_err,*) "  unknown:      ", trim(names(i))
                any_err = .true.
            end if
        end do

        if (any_err .and. halt) error stop "Program stopped."

        return

    end subroutine nml_validate

    subroutine nml_read_string(filename,group,name,value,comment,init,defaults_file,defaults_group)

        implicit none 

        character(len=*), intent(INOUT) :: value 
        character(len=*), intent(IN)    :: filename, group, name 
        character(len=*), intent(INOUT), optional :: comment 
        logical, optional :: init
        character(len=*), intent(IN), optional :: defaults_file, defaults_group
        logical :: init_var 
        character(len=256) :: value_str 

        init_var = .FALSE. 
        if (present(init)) init_var = init 
        if (init_var) value = "" 

        ! First find parameter value as a string 
        value_str = ""
        call nml_read_internal(filename,group,name,value_str,comment,defaults_file,defaults_group)

        if (value_str /= "") value = trim(value_str)

        if (VERBOSE) call nml_print(name,value,comment)  ! Check

        return 

    end subroutine nml_read_string

    subroutine nml_read_double(filename,group,name,value,comment,init,defaults_file,defaults_group)

        implicit none 

        double precision, intent(INOUT) :: value 
        character(len=*), intent(IN)    :: filename, group, name 
        character(len=*), intent(INOUT), optional :: comment   
        logical, optional :: init
        character(len=*), intent(IN), optional :: defaults_file, defaults_group
        logical :: init_var 
        character(len=256) :: value_str 

        init_var = .FALSE. 
        if (present(init)) init_var = init 
        if (init_var) value = 0.d0 

        ! First find parameter value as a string 
        value_str = ""
        call nml_read_internal(filename,group,name,value_str,comment,defaults_file,defaults_group)
        if (value_str /= "") value = string_to_double(value_str)

        if (VERBOSE) call nml_print(name,value,comment)  ! Check

        return 

    end subroutine nml_read_double 

    subroutine nml_read_float(filename,group,name,value,comment,init,defaults_file,defaults_group)

        implicit none 

        real(4), intent(INOUT) :: value 
        character(len=*), intent(IN)    :: filename, group, name 
        character(len=*), intent(INOUT), optional :: comment  
        logical, optional :: init
        character(len=*), intent(IN), optional :: defaults_file, defaults_group
        logical :: init_var  
        character(len=256) :: value_str 

        init_var = .FALSE. 
        if (present(init)) init_var = init 
        if (init_var) value = 0.0 

        ! First find parameter value as a string 
        value_str = ""
        call nml_read_internal(filename,group,name,value_str,comment,defaults_file,defaults_group)
        if (value_str /= "") value = real(string_to_double(value_str))

        if (VERBOSE) call nml_print(name,value,comment)  ! Check

        return 

    end subroutine nml_read_float 

    subroutine nml_read_integer(filename,group,name,value,comment,init,defaults_file,defaults_group)

        implicit none 

        integer, intent(INOUT) :: value 
        character(len=*), intent(IN)    :: filename, group, name 
        character(len=*), intent(INOUT), optional :: comment 
        logical, optional :: init
        character(len=*), intent(IN), optional :: defaults_file, defaults_group
        logical :: init_var   
        character(len=256) :: value_str 

        init_var = .FALSE. 
        if (present(init)) init_var = init 
        if (init_var) value = 0 

        ! First find parameter value as a string 
        value_str = ""
        call nml_read_internal(filename,group,name,value_str,comment,defaults_file,defaults_group)

        if (value_str /= "") then 
            value = nint(string_to_double(value_str))
        end if 

        if (VERBOSE) call nml_print(name,value,comment)  ! Check

        return 

    end subroutine nml_read_integer

    subroutine nml_read_logical(filename,group,name,value,comment,init,defaults_file,defaults_group)

        implicit none 

        logical, intent(INOUT) :: value 
        character(len=*), intent(IN)    :: filename, group, name 
        character(len=*), intent(INOUT), optional :: comment 
        logical, optional :: init
        character(len=*), intent(IN), optional :: defaults_file, defaults_group
        logical :: init_var   
        character(len=256) :: value_str 

        init_var = .FALSE. 
        if (present(init)) init_var = init 
        if (init_var) value = .FALSE. 

        ! First find parameter value as a string 
        value_str = ""
        call nml_read_internal(filename,group,name,value_str,comment,defaults_file,defaults_group)
        if (value_str /= "") value = string_to_logical(value_str)

        if (VERBOSE) call nml_print(name,value,comment)  ! Check

        return 

    end subroutine nml_read_logical 

    !! Vectors 

    subroutine nml_read_string_vector(filename,group,name,value,comment,init,defaults_file,defaults_group)

        implicit none 

        character(len=*), intent(INOUT) :: value(:) 
        character(len=*), intent(IN)    :: filename, group, name 
        character(len=*), intent(INOUT), optional :: comment 
        logical, optional :: init
        character(len=*), intent(IN), optional :: defaults_file, defaults_group
        logical :: init_var  
        character(len=256) :: value_str 

        init_var = .FALSE. 
        if (present(init)) init_var = init 
        if (init_var) value(:) = "" 

        ! First find parameter value as a string 
        value_str = ""
        call nml_read_internal(filename,group,name,value_str,comment,defaults_file,defaults_group)

        if (value_str /= "") call string_to_vector(value_str,value)

        if (VERBOSE) call nml_print(name,value,comment)  ! Check

        return 

    end subroutine nml_read_string_vector 

    subroutine nml_read_double_vector(filename,group,name,value,comment,init,defaults_file,defaults_group)

        implicit none 

        double precision, intent(INOUT) :: value(:) 
        character(len=*), intent(IN)    :: filename, group, name 
        character(len=*), intent(INOUT), optional :: comment   
        character(len=256) :: value_str, value_str_vec(size(value))
        logical, optional :: init
        character(len=*), intent(IN), optional :: defaults_file, defaults_group
        logical :: init_var  
        integer :: q 

        init_var = .FALSE. 
        if (present(init)) init_var = init 
        if (init_var) value(:) = 0.d0

        ! First find parameter value as a string 
        value_str = ""
        call nml_read_internal(filename,group,name,value_str,comment,defaults_file,defaults_group)

        if (value_str /= "") then
            call string_to_vector(value_str,value_str_vec)
            do q = 1, size(value)
                if (trim(value_str_vec(q)) /= "") then 
                    value(q) = string_to_double(trim(adjustl(value_str_vec(q))))
                end if 

            end do 
        end if 

        if (VERBOSE) call nml_print(name,value,comment)  ! Check

        return 

    end subroutine nml_read_double_vector 

    subroutine nml_read_float_vector(filename,group,name,value,comment,init,defaults_file,defaults_group)

        implicit none 

        real(4), intent(INOUT) :: value(:) 
        character(len=*), intent(IN)    :: filename, group, name 
        character(len=*), intent(INOUT), optional :: comment   
        character(len=256) :: value_str, value_str_vec(size(value)) 
        logical, optional :: init
        character(len=*), intent(IN), optional :: defaults_file, defaults_group
        logical :: init_var
        integer :: q 

        init_var = .FALSE. 
        if (present(init)) init_var = init 
        if (init_var) value(:) = 0.0

        ! First find parameter value as a string 
        value_str = ""
        call nml_read_internal(filename,group,name,value_str,comment,defaults_file,defaults_group)

        if (value_str /= "") then
            call string_to_vector(value_str,value_str_vec)
            do q = 1, size(value)
                if (trim(value_str_vec(q)) /= "") then 
                    value(q) = real(string_to_double(trim(adjustl(value_str_vec(q)))))
                end if 

            end do 
        end if 

        if (VERBOSE) call nml_print(name,value,comment)  ! Check

        return 

    end subroutine nml_read_float_vector 

    subroutine nml_read_integer_vector(filename,group,name,value,comment,init,defaults_file,defaults_group)

        implicit none 

        integer, intent(INOUT) :: value(:) 
        character(len=*), intent(IN)    :: filename, group, name 
        character(len=*), intent(INOUT), optional :: comment   
        character(len=256) :: value_str, value_str_vec(size(value))
        logical, optional :: init
        character(len=*), intent(IN), optional :: defaults_file, defaults_group
        logical :: init_var 
        integer :: q 

        init_var = .FALSE. 
        if (present(init)) init_var = init 
        if (init_var) value(:) = 0

        ! First find parameter value as a string 
        value_str = ""
        call nml_read_internal(filename,group,name,value_str,comment,defaults_file,defaults_group)

        if (value_str /= "") then
            call string_to_vector(value_str,value_str_vec)
            do q = 1, size(value)
                if (trim(value_str_vec(q)) /= "") then 
                    value(q) = nint(string_to_double(trim(adjustl(value_str_vec(q)))))
                end if 

            end do 
        end if 

        if (VERBOSE) call nml_print(name,value,comment)  ! Check

        return 

    end subroutine nml_read_integer_vector

    subroutine nml_read_logical_vector(filename,group,name,value,comment,init,defaults_file,defaults_group)

        implicit none 

        logical, intent(INOUT) :: value(:) 
        character(len=*), intent(IN)    :: filename, group, name 
        character(len=*), intent(INOUT), optional :: comment   
        character(len=256) :: value_str, value_str_vec(size(value))
        logical, optional :: init
        character(len=*), intent(IN), optional :: defaults_file, defaults_group
        logical :: init_var 
        integer :: q 

        init_var = .FALSE. 
        if (present(init)) init_var = init 
        if (init_var) value(:) = .FALSE.

        ! First find parameter value as a string 
        value_str = ""
        call nml_read_internal(filename,group,name,value_str,comment,defaults_file,defaults_group)

        if (value_str /= "") then
            call string_to_vector(value_str,value_str_vec)
            do q = 1, size(value)
                if (trim(value_str_vec(q)) /= "") then 
                    value(q) = string_to_logical(trim(adjustl(value_str_vec(q))))
                end if 

            end do 
        end if 

        if (VERBOSE) call nml_print(name,value,comment)  ! Check

        return 

    end subroutine nml_read_logical_vector

    ! =============================================================
    !
    ! nml parameter writing functions
    !
    ! =============================================================

    ! This is the basic nml writing subroutine
    ! All interfaces use this to write the parameter to file, after
    ! converting the parameter to a string
    subroutine nml_write_string(filename,group,name,value,comment)

        implicit none 

        character(len=*), intent(IN) :: filename, group, name, value 
        character(len=*), intent(IN), optional :: comment  


        return 

    end subroutine nml_write_string 

    ! =============================================================
    !
    ! nml line printing functions
    !
    ! =============================================================

    ! This is the basic routine for printing a parameter to a formatted line
    ! All other interfaces use this routine after converting to a string.
    subroutine nml_print_string(name,value,comment,io,no_quotes)

        implicit none 
        character(len=*), intent(in) :: name, value 
        character(len=*), intent(in), optional :: comment 
        integer, intent(in), optional :: io 
        integer :: io_val 
        logical, intent(in), optional :: no_quotes
        logical :: no_quotes_loc
        character(len=1000) :: line
        character(len=500)  :: comment1 
        character(len=len(value))  :: val_repr

        io_val = 6 
        if (present(io)) io_val = io 
        val_repr = value

        no_quotes_loc = .false.
        if (present(no_quotes)) no_quotes_loc = no_quotes
        if (.not.no_quotes_loc) val_repr = '"'//trim(value)//'"'

        comment1 = "" 
        if (present(comment)) comment1 = "   "//trim(comment)

        write(line,"(a)") "    "//trim(name)//" = "//trim(val_repr)//trim(comment1)
        write(io_val,*) trim(line)

        return 

    end subroutine nml_print_string 

    subroutine nml_print_double(name,value,comment,io)

        implicit none 
        double precision :: value
        character(len=*) :: name 
        character(len=*), optional :: comment
        integer, optional :: io 
        character(len=500) :: value_str  

        write(value_str,*) value 
        if (VERBOSE) call nml_print_string(name,value_str,comment,io,no_quotes=.true.)

        return 

    end subroutine nml_print_double
    
    subroutine nml_print_float(name,value,comment,io)

        implicit none 
        real(4) :: value
        character(len=*) :: name 
        character(len=*), optional :: comment
        integer, optional :: io 
        character(len=500) :: value_str  

        write(value_str,*) value 
        if (VERBOSE) call nml_print_string(name,value_str,comment,io)

        return 

    end subroutine nml_print_float
    
    subroutine nml_print_integer(name,value,comment,io)

        implicit none 
        integer :: value
        character(len=*) :: name 
        character(len=*), optional :: comment
        integer, optional :: io 
        character(len=500) :: value_str  

        write(value_str,*) value 
        if (VERBOSE) call nml_print_string(name,value_str,comment,io,no_quotes=.true.)

        return 

    end subroutine nml_print_integer

    subroutine nml_print_logical(name,value,comment,io)

        implicit none 
        logical :: value
        character(len=*) :: name 
        character(len=*), optional :: comment
        integer, optional :: io 
        character(len=500) :: value_str  

        value_str = "F"
        if (value) value_str = "T" 
        if (VERBOSE) call nml_print_string(name,value_str,comment,io,no_quotes=.true.)

        return 

    end subroutine nml_print_logical
    
    !! Vectors

    ! Note appended to a truncated vector print, e.g. " ... (100000 total)".
    ! Returns "" when the vector is short enough to print in full.
    function nml_print_vector_note(n) result(note)

        implicit none
        integer, intent(in) :: n
        character(len=:), allocatable :: note
        character(len=32) :: tmp

        if (n .gt. nml_print_vector_max) then
            write(tmp,"(i0)") n
            note = " ... ("//trim(adjustl(tmp))//" total)"
        else
            note = ""
        end if

        return

    end function nml_print_vector_note

    subroutine nml_print_string_vector(name,value,comment,io)

        implicit none
        character(len=*) :: value(:)
        character(len=*) :: name
        character(len=*), optional :: comment
        integer, optional :: io
        character(len=:), allocatable :: value_str
        character(len=16) :: tmp
        integer :: q

        value_str = '"'//trim(value(1))//'"'
        do q = 2, min(size(value),nml_print_vector_max)
            value_str = value_str // " " // trim(adjustl(value(q)))
        end do
        value_str = value_str // nml_print_vector_note(size(value))

        if (VERBOSE) call nml_print_string(name,value_str,comment,io,no_quotes=.true.)

        return

    end subroutine nml_print_string_vector
    
    subroutine nml_print_double_vector(name,value,comment,io)

        implicit none 
        double precision :: value(:)
        character(len=*) :: name 
        character(len=*), optional :: comment
        integer, optional :: io 
        character(len=:), allocatable :: value_str
        character(len=16) :: tmp
        integer :: q 

        value_str = ""
        do q = 1, min(size(value),nml_print_vector_max)
            write(tmp,"(g12.3)") value(q)
            value_str = value_str // " " // trim(adjustl(tmp))
        end do
        value_str = value_str // nml_print_vector_note(size(value))

        if (VERBOSE) call nml_print_string(name,value_str,comment,io,no_quotes=.true.)

        return

    end subroutine nml_print_double_vector
    
    subroutine nml_print_float_vector(name,value,comment,io)

        implicit none 
        real(4) :: value(:)
        character(len=*) :: name 
        character(len=*), optional :: comment
        integer, optional :: io 
        character(len=:), allocatable :: value_str
        character(len=16) :: tmp
        integer :: q 

        value_str = ""
        do q = 1, min(size(value),nml_print_vector_max)
            write(tmp,"(g12.3)") value(q)
            value_str = value_str // " " // trim(adjustl(tmp))
        end do
        value_str = value_str // nml_print_vector_note(size(value))

        if (VERBOSE) call nml_print_string(name,value_str,comment,io,no_quotes=.true.)

        return

    end subroutine nml_print_float_vector
    
    subroutine nml_print_integer_vector(name,value,comment,io)

        implicit none
        integer :: value(:)
        character(len=*) :: name
        character(len=*), optional :: comment
        integer, optional :: io
        character(len=:), allocatable :: value_str
        character(len=16) :: tmp
        integer :: q

        value_str = ""
        do q = 1, min(size(value),nml_print_vector_max)
            write(tmp,"(i12)") value(q)
            value_str = value_str // " " // trim(adjustl(tmp))
        end do
        value_str = value_str // nml_print_vector_note(size(value))

        if (VERBOSE) call nml_print_string(name,value_str,comment,io,no_quotes=.true.)

        return

    end subroutine nml_print_integer_vector

    subroutine nml_print_logical_vector(name,value,comment,io)

        implicit none 
        logical :: value(:)
        character(len=*) :: name 
        character(len=*), optional :: comment
        integer, optional :: io 
        character(len=:), allocatable :: value_str
        character(len=16) :: tmp 
        integer :: q 

        value_str = ""
        do q = 1, min(size(value),nml_print_vector_max)
            if (value(q)) then
                value_str = value_str // " " // "T"
            else
                value_str = value_str // " " // "F"
            end if
        end do
        value_str = value_str // nml_print_vector_note(size(value))

        if (VERBOSE) call nml_print_string(name,value_str,comment,io,no_quotes=.true.)

        return 

    end subroutine nml_print_logical_vector
    

    ! =============================================================
    !
    ! Type conversion functions
    !
    ! =============================================================

    function string_to_double(string) result(value)

        implicit none 

        character(len=*), intent(IN) :: string 
        double precision :: value 

        character(len=256) :: tmpstr 
        integer :: stat, n
        double precision :: x 

        tmpstr = trim(adjustl(string))
        n      = len_trim(tmpstr)

        read(tmpstr(1:n),*,IOSTAT=stat) x

        value = 0
        if (stat .eq. 0) then 
            value = x 
        else
            n = len_trim(tmpstr)-1
            READ(tmpstr(1:n),*,IOSTAT=stat) x
            if (stat .ne. 0) then 
                write(*,*) "nml:: Error converting string to number!"
                write(*,*) "|",trim(tmpstr),"|",n,stat,x
                error stop
            else
                value = x 
            end if 
        end if 

        return 

    end function string_to_double

    function string_to_logical(string) result(value)

        implicit none 

        character(len=*), intent(IN) :: string 
        logical :: value 

        character(len=256) :: tmpstr 
        integer :: stat, n
        double precision :: x 

        tmpstr = trim(adjustl(string))
        
        select case(trim(tmpstr))
            case("T","True","TRUE","true",".TRUE.", ".true.")
                value = .TRUE. 
            case("F","False","FALSE","false",".FALSE.", ".false.")
                value = .FALSE. 
            case DEFAULT
                write(*,*) "nml:: Error reading logical parameter."
                error stop
        end select  

        return 

    end function string_to_logical

    subroutine string_to_vector(string,value)

        implicit none 

        character(len=*), intent(IN) :: string 
        character(len=*) :: value(:)
        character(len=256) :: tmpvec(size(value))
        character(len=256) :: tmpstr, fmt 
        integer :: stat, n, q, q1, q2, j 

        tmpstr = trim(adjustl(string))
        n      = len_trim(tmpstr)+2

        tmpvec(:) = "" 

        q1 = 1 
        do q = 1, size(tmpvec)
            q2 = index(tmpstr(q1:n)," ") + q1
            if (q2 .gt. q1 .and. q2 .le. n) then 
                tmpvec(q) = tmpstr(q1:q2-1)
                q1 = q2

                ! Make sure gaps of more than one space are properly handled
                do j = 1, 1000
                    if (tmpstr(q1:q1) == " ") q1 = q1+1
                    if (q1 .ge. n) exit 
                end do 

!                 ! Eliminate quotes
!                 q2 = len_trim(tmpvec(q))
!                 if (tmpvec(q)(1:1) == '"') tmpvec(q) = trim(adjustl(tmpvec(q)(2:q2)))
!                 q2 = len_trim(tmpvec(q))
!                 if (tmpvec(q)(q2:q2) == '"') tmpvec(q) = trim(tmpvec(q)(1:q2-1))
                ! Remove quotes around string if they exist 
                call remove_quotes_comma(tmpvec(q))
            
            end if 
        end do 
        
        value = tmpvec 

        return 

    end subroutine string_to_vector

    ! =============================================================
    !
    ! Helper functions
    !
    ! =============================================================

    ! This function parses a namelist line into pieces, determining
    ! whether it is a blank line (-2), comment (-1), group name (1)
    ! end-of-group (2), or a parameter line (3)
    subroutine parse_line(line,linetype,name,value,comment)

        implicit none 
        character(len=*), intent(IN)    :: line
        character(len=*), intent(INOUT) :: name, value, comment 
        integer :: linetype

        character(len=1001) :: line1 
        integer :: q, q1, q2

        name     = ""
        value    = ""
        comment  = "" 

        line1 = trim(adjustl(line))

        if (trim(line1) == "") then         ! Blank line
            linetype = -2 
            
        else if (line1(1:1) == "!") then    ! Comment line 
            linetype = -1 
            comment = trim(line1)

        else if (line1(1:1) == "&") then    ! Group name 
            linetype = 1 
            q = len_trim(line1)
            name     = line1(2:q)

        else if (line1(1:1) == "/") then    ! End of group 
            linetype = 2 

        else   ! Line must contain parameter to read
            linetype = 3

            q = index(line1,"=")
            if (q == 0) then 
                write(*,*) "nml:: Error reading namelist file."
                write(*,*) "No '=' found on parameter line."
                write(*,*) "line: ", trim(line1)
                error stop
            end if 

            name = trim(adjustl(line1(1:q-1)))

            q1 = index(line1,"!")
            q2 = len_trim(line1)

            if (q1 > 0) then 
                comment = trim(adjustl(line1(q1:q2)))
                value   = trim(adjustl(line1(q+1:q1-1)))
            else
                value   = trim(adjustl(line1(q+1:q2)))
            end if 

            ! Remove quotes around string, and final line comma, if they exist
            call remove_quotes_comma(value)

        end if 

        return 

    end subroutine parse_line 

    subroutine remove_quotes_comma(string)

        implicit none 
        character(len=*), intent(INOUT) :: string 
        integer :: i, n 

!         ! Eliminate quotes
!         n = len_trim(string)
!         if (n == 1 .and. trim(string) == '"') then 
!             string = ""
!         else if (n > 0) then 
!             if (string(1:1) == '"') string = trim(adjustl(string(2:n)))
!             n = len_trim(string)
!             if (n > 1  .and. string(n:n) == '"') string = trim(string(1:n-1))
!             if (n == 1 .and. string(n:n) == '"') string = ""
            
!         end if 

        ! Eliminate quotes
        n = len_trim(string)
        do i = 1,n 
            if (string(i:i) == '"' .or. string(i:i) == "'") string(i:i) = " "
        end do 
        string = trim(adjustl(string))

        ! Remove final comma too
        n = len_trim(string)
        if (n > 0) then 
            if (string(n:n) == ",") string(n:n) = " "
            string = trim(adjustl(string))
        end if 
        
        return 

    end subroutine remove_quotes_comma



    integer function newunit(unit)
        ! This is a simple function to search for an available unit.
        ! LUN_MIN and LUN_MAX define the range of possible LUNs to check.
        ! The UNIT value is returned by the function, and also by the optional
        ! argument. This allows the function to be used directly in an OPEN
        ! statement, and optionally save the result in a local variable.
        ! If no units are available, -1 is returned.
        ! Modified from: http://fortranwiki.org/fortran/show/newunit

        implicit none 

        integer, intent(out), optional :: unit

        ! Local variables
        integer, parameter :: LUN_MIN=10, LUN_MAX=1000
        logical :: opened
        integer :: lun

        newunit=-1

        ! Search for new unit
        do lun=LUN_MIN,LUN_MAX
            inquire(unit=lun,opened=opened)
            if (.not. opened) then
                newunit=lun
                exit
            end if
        end do

        ! Assign new unit if desired
        if (present(unit)) unit=newunit

        return 

    end function newunit

end module nml 

