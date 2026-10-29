module variable_io

    use, intrinsic :: iso_fortran_env, only : input_unit, output_unit, error_unit
    use ncio

    implicit none

    type var_io_type
        integer :: id
        character(len=32)  :: varname
        character(len=32)  :: dimnames
        character(len=32)  :: units
        character(len=256) :: long_name

        integer :: ndims
        character(len=32), allocatable :: dims(:)

        ! Precision of the values written (ncio bit rounding), 0 for lossless:
        ! absolute (units of the field) and relative to the magnitude
        double precision :: prec_abs = 0.d0
        double precision :: prec_rel = 0.d0
    end type

    interface var_io_write
        module procedure var_io_write_float_1D, var_io_write_float_2D, var_io_write_float_3D
        module procedure var_io_write_double_1D, var_io_write_double_2D, var_io_write_double_3D
        module procedure var_io_write_int_1D, var_io_write_int_2D, var_io_write_int_3D
    end interface

    private
    public :: var_io_type
    public :: find_var_io_in_file
    public :: find_var_io_in_table
    public :: load_var_io_table
    public :: var_io_print
    public :: var_io_lossless
    public :: var_io_write

contains

    subroutine find_var_io_in_file(var,varname,filename)

        implicit none

        type(var_io_type), intent(OUT) :: var
        character(len=*),  intent(IN)  :: varname
        character(len=*),  intent(IN)  :: filename

        ! Local variables 
        type(var_io_type), allocatable :: var_table(:)
        
        ! Load the variable table
        call load_var_io_table(var_table,filename)

        ! Find the current variable in the table
        call find_var_io_in_table(var,varname,var_table)

        return

    end subroutine find_var_io_in_file

    subroutine find_var_io_in_table(var,varname,var_table,with_error)

        implicit none
        
        type(var_io_type), intent(OUT) :: var
        character(len=*),  intent(IN)  :: varname 
        type(var_io_type), intent(IN)  :: var_table(:)
        logical,           intent(IN), optional :: with_error

        ! Local variables
        integer :: n 
        logical :: found

        found = .FALSE.
        var%varname = "none"

        do n = 1, size(var_table)
            if (trim(var_table(n)%varname) .eq. trim(varname)) then
                var = var_table(n)
                found = .TRUE.
                exit
            end if
        end do

        if (.not. found .and. present(with_error)) then
            if (with_error) then
                write(error_unit,*) "find_var_io_in_table:: Error: variable not found in variable table."
                write(error_unit,*) "varname = ", trim(varname)
                error stop 1
            end if
        end if

        return

    end subroutine find_var_io_in_table


    subroutine load_var_io_table(var_table,filename)

        implicit none

        type(var_io_type), allocatable, intent(OUT) :: var_table(:)
        character(len=*), intent(IN) :: filename

        ! Local variables
        integer, parameter :: nmax = 1000
        type(var_io_type) :: vt(nmax)       ! Max 1000 variables
        integer :: n, ntot
        character(len=1024) :: line
        integer :: io, check

        ! Open the file for reading
        open(newunit=io, file=trim(filename), status='old', action='read')

        ! Skip the header (header finishes with first line starting with "---" )
        do n = 1, 50
            read(unit=io, fmt='(a)') line
            if (line(1:3) .eq. "|--") exit
        end do

        ! Read the table rows, stored contiguously in vt(1:ntot)
        ! (other lines after the header are skipped)
        ntot = 0
        do
            read(unit=io, fmt='(a)',iostat=check) line

            ! Exit at end of file
            if (check .eq. -1) exit

            if (line(1:2) .eq. "| ") then

                if (ntot .eq. nmax) then
                    write(error_unit,*) "load_var_table:: Error: maximum number of variables read in, but &
                    &end of the table was not reached. Increase value of nmax in this routine to be safe."
                    write(error_unit,*) "nmax = ", nmax
                    write(error_unit,*) "filename = ", trim(filename)
                    write(error_unit,*) "Last line read in: "
                    write(error_unit,*) trim(line)
                    error stop 1
                end if

                ntot = ntot+1

                ! Parse the line
                call parse_line_to_variable(vt(ntot),line)

                ! Parse the dimensions
                call parse_dims(vt(ntot)%dims,vt(ntot)%dimnames)
                vt(ntot)%ndims = size(vt(ntot)%dims)
                call var_io_print(vt(ntot))
            end if

        end do

        ! Close the file
        close(io)

        if (allocated(var_table)) deallocate(var_table)
        allocate(var_table(ntot))
        var_table = vt(1:ntot)

        return

    end subroutine load_var_io_table

    subroutine parse_line_to_variable(var,line)
        ! Parse one table row:
        !   | id | variable | dimensions | units | long_name |
        ! with an optional 6th column `prec` (see parse_prec). Without it,
        ! the variable is written lossless (prec_abs = prec_rel = 0).

        implicit none

        type(var_io_type), intent(INOUT) :: var
        character(len=*),  intent(IN)    :: line

        ! Local variables
        integer, parameter :: ncol_max = 6
        character(len=1024) :: temp_line
        character(len=256)  :: col(ncol_max)
        integer :: pos, ncol
        integer :: ios

        ! Initialize
        col(:) = ""

        ! Remove leading spaces and the leading '|'
        temp_line = adjustl(line)
        temp_line = temp_line(2:)

        ! Split the columns between the '|'
        ncol = 0
        do
            pos = index(temp_line, '|')
            if (pos .eq. 0) exit
            ncol = ncol + 1
            if (ncol .gt. ncol_max) exit
            col(ncol) = adjustl(temp_line(:pos-1))
            temp_line = temp_line(pos+1:)
        end do

        if (ncol .ne. 5 .and. ncol .ne. 6) then
            write(error_unit,*) "variable_io:: parse_line_to_variable:: Error: a table row must have &
            &5 columns (id, variable, dimensions, units, long_name) or 6 (+ prec)."
            write(error_unit,*) "line = ", trim(line)
            error stop 1
        end if

        ! Extract id
        if (.not. trim(col(1)) .eq. "") then
            ! Convert to integer
            read(col(1),"(i20)",iostat=ios) var%id
            if (ios .ne. 0) then
                write(error_unit,*) "variable_io:: parse_line_to_variable:: Error: id must be an integer."
                write(error_unit,*) "id = ", trim(col(1))
                error stop 1
            end if
        else
            var%id = 0
        end if

        var%varname   = col(2)
        var%dimnames  = col(3)
        var%units     = col(4)
        var%long_name = col(5)

        ! Extract precision (lossless without the prec column)
        var%prec_abs = 0.d0
        var%prec_rel = 0.d0
        if (ncol .eq. 6) call parse_prec(var%prec_abs,var%prec_rel,col(6),line)

        return

    end subroutine parse_line_to_variable

    subroutine parse_prec(prec_abs,prec_rel,str,line)
        ! Parse the `prec` column of a table row:
        !   -            lossless (also an empty entry)
        !   1e-3         absolute precision, in the units of the field
        !   1e-4r        relative precision
        !   1e-3,1e-4r   both (ncio keeps the larger error of the two)

        implicit none

        double precision, intent(OUT) :: prec_abs
        double precision, intent(OUT) :: prec_rel
        character(len=*), intent(IN)  :: str
        character(len=*), intent(IN)  :: line         ! Full row, for error messages

        ! Local variables
        character(len=256) :: temp_line
        character(len=64)  :: tok
        integer :: pos, k, nc, ios
        logical :: is_rel, has_abs, has_rel
        double precision :: val

        prec_abs = 0.d0
        prec_rel = 0.d0

        temp_line = adjustl(str)
        if (trim(temp_line) .eq. "" .or. trim(temp_line) .eq. "-") return

        has_abs = .FALSE.
        has_rel = .FALSE.

        do k = 1, 3

            ! Next comma-separated entry
            pos = index(temp_line, ',')
            if (pos .eq. 0) then
                tok = adjustl(temp_line)
            else
                tok = adjustl(temp_line(:pos-1))
                temp_line = temp_line(pos+1:)
            end if

            if (k .eq. 3) call parse_prec_error("at most two entries, one absolute and one relative.")

            ! Relative if suffixed by 'r'
            nc = len_trim(tok)
            is_rel = .FALSE.
            if (nc .gt. 0) is_rel = (tok(nc:nc) .eq. "r")
            if (is_rel) nc = nc - 1

            ! Read the value: a plain non-negative number
            if (nc .eq. 0) call parse_prec_error("empty entry.")
            if (verify(tok(1:nc), "0123456789.eEdD+-") .ne. 0) &
                call parse_prec_error("entry is not a number: "//trim(tok))
            read(tok(1:nc),*,iostat=ios) val
            if (ios .ne. 0) call parse_prec_error("entry is not a number: "//trim(tok))
            if (val .lt. 0.d0) call parse_prec_error("entry must be >= 0: "//trim(tok))

            if (is_rel) then
                if (has_rel) call parse_prec_error("relative precision given twice.")
                prec_rel = val
                has_rel  = .TRUE.
            else
                if (has_abs) call parse_prec_error("absolute precision given twice.")
                prec_abs = val
                has_abs  = .TRUE.
            end if

            if (pos .eq. 0) exit

        end do

        return

    contains

        subroutine parse_prec_error(msg)

            implicit none

            character(len=*), intent(IN) :: msg

            write(error_unit,*) "variable_io:: parse_prec:: Error: bad prec entry, "//trim(msg)
            write(error_unit,*) "  syntax: - (lossless), 1e-3 (absolute), 1e-4r (relative), 1e-3,1e-4r (both)"
            write(error_unit,*) "  prec = ", trim(str)
            write(error_unit,*) "  line = ", trim(line)
            error stop 1

        end subroutine parse_prec_error

    end subroutine parse_prec


    subroutine parse_dims(dims,dimnames)

        implicit none

        character(len=*), allocatable, intent(OUT) :: dims(:)
        character(len=*), intent(IN)  :: dimnames

        ! Local variables
        character(len=200) :: temp_line
        integer :: pos1, pos2, length
        character(len=32) :: d(10)
        integer :: n, ndims 

        ! Initialize
        d(:) = ""

        ! Remove leading/trailing spaces
        temp_line = adjustl(trim(dimnames))

            
        ndims = 0

        do n = 1, size(d)

            pos1 = index(temp_line, ',')

            if (pos1 .eq. 0) then
                ! Only one dimension currently available
                ndims = n
                d(ndims) = adjustl(trim(temp_line))
                exit

            else
                ndims = n
                d(ndims) = adjustl(trim(temp_line(1:pos1-1)))
                temp_line = temp_line(pos1+1:)
            end if

        end do

        if (allocated(dims)) deallocate(dims)
        allocate(dims(ndims))
        dims = d(1:ndims)

        return

    end subroutine parse_dims

    subroutine var_io_print(var)

        implicit none

        type(var_io_type), intent(IN) :: var

        write(*,"(i4,3a20,2x,a50)") var%id, trim(var%varname), trim(var%dimnames), trim(var%units), trim(var%long_name)

        return

    end subroutine var_io_print

    function var_io_lossless(v) result(vl)
        ! Copy of v with the precision reset to lossless, to write a
        ! table variable without rounding (e.g. restart files)

        implicit none

        type(var_io_type), intent(IN) :: v
        type(var_io_type) :: vl

        vl = v
        vl%prec_abs = 0.d0
        vl%prec_rel = 0.d0

        return

    end function var_io_lossless

    ! === var_io_write: write one table variable for one time record ===
    !
    ! call var_io_write(filename,v,dat,n,[ncid],[dims],[missing_value],[grid_mapping])
    !
    ! Writes dat to record n of variable v%varname with the dimensions
    ! v%dims + "time" (or the full list `dims`, incl. the time dimension),
    ! the units and long_name of v, and, for float/double data, the
    ! precision of v (prec_abs, prec_rel; integer data is never rounded).

    subroutine var_io_write_float_1D(filename,v,dat,n,ncid,dims,missing_value,grid_mapping)

        implicit none

        character(len=*),  intent(IN) :: filename
        type(var_io_type), intent(IN) :: v
        real(4), contiguous, intent(IN) :: dat(:)
        integer,           intent(IN) :: n
        integer,           intent(IN), optional :: ncid
        character(len=*),  intent(IN), optional :: dims(:)
        real(4), intent(IN), optional :: missing_value
        character(len=*),  intent(IN), optional :: grid_mapping

        ! Local variables
        character(len=32), allocatable :: dims_out(:)
        integer, allocatable :: start(:), count(:)

        call var_io_write_dims(dims_out,start,count,v,shape(dat),n,dims)

        call nc_write(filename,trim(v%varname),dat,dims=dims_out,start=start,count=count, &
                      units=trim(v%units),long_name=trim(v%long_name), &
                      missing_value=missing_value,grid_mapping=grid_mapping,ncid=ncid, &
                      prec_abs=real(v%prec_abs,4),prec_rel=real(v%prec_rel,4))

        return

    end subroutine var_io_write_float_1D

    subroutine var_io_write_float_2D(filename,v,dat,n,ncid,dims,missing_value,grid_mapping)

        implicit none

        character(len=*),  intent(IN) :: filename
        type(var_io_type), intent(IN) :: v
        real(4), contiguous, intent(IN) :: dat(:,:)
        integer,           intent(IN) :: n
        integer,           intent(IN), optional :: ncid
        character(len=*),  intent(IN), optional :: dims(:)
        real(4), intent(IN), optional :: missing_value
        character(len=*),  intent(IN), optional :: grid_mapping

        ! Local variables
        character(len=32), allocatable :: dims_out(:)
        integer, allocatable :: start(:), count(:)

        call var_io_write_dims(dims_out,start,count,v,shape(dat),n,dims)

        call nc_write(filename,trim(v%varname),dat,dims=dims_out,start=start,count=count, &
                      units=trim(v%units),long_name=trim(v%long_name), &
                      missing_value=missing_value,grid_mapping=grid_mapping,ncid=ncid, &
                      prec_abs=real(v%prec_abs,4),prec_rel=real(v%prec_rel,4))

        return

    end subroutine var_io_write_float_2D

    subroutine var_io_write_float_3D(filename,v,dat,n,ncid,dims,missing_value,grid_mapping)

        implicit none

        character(len=*),  intent(IN) :: filename
        type(var_io_type), intent(IN) :: v
        real(4), contiguous, intent(IN) :: dat(:,:,:)
        integer,           intent(IN) :: n
        integer,           intent(IN), optional :: ncid
        character(len=*),  intent(IN), optional :: dims(:)
        real(4), intent(IN), optional :: missing_value
        character(len=*),  intent(IN), optional :: grid_mapping

        ! Local variables
        character(len=32), allocatable :: dims_out(:)
        integer, allocatable :: start(:), count(:)

        call var_io_write_dims(dims_out,start,count,v,shape(dat),n,dims)

        call nc_write(filename,trim(v%varname),dat,dims=dims_out,start=start,count=count, &
                      units=trim(v%units),long_name=trim(v%long_name), &
                      missing_value=missing_value,grid_mapping=grid_mapping,ncid=ncid, &
                      prec_abs=real(v%prec_abs,4),prec_rel=real(v%prec_rel,4))

        return

    end subroutine var_io_write_float_3D

    subroutine var_io_write_double_1D(filename,v,dat,n,ncid,dims,missing_value,grid_mapping)

        implicit none

        character(len=*),  intent(IN) :: filename
        type(var_io_type), intent(IN) :: v
        real(8), contiguous, intent(IN) :: dat(:)
        integer,           intent(IN) :: n
        integer,           intent(IN), optional :: ncid
        character(len=*),  intent(IN), optional :: dims(:)
        real(8), intent(IN), optional :: missing_value
        character(len=*),  intent(IN), optional :: grid_mapping

        ! Local variables
        character(len=32), allocatable :: dims_out(:)
        integer, allocatable :: start(:), count(:)

        call var_io_write_dims(dims_out,start,count,v,shape(dat),n,dims)

        call nc_write(filename,trim(v%varname),dat,dims=dims_out,start=start,count=count, &
                      units=trim(v%units),long_name=trim(v%long_name), &
                      missing_value=missing_value,grid_mapping=grid_mapping,ncid=ncid, &
                      prec_abs=real(v%prec_abs,8),prec_rel=real(v%prec_rel,8))

        return

    end subroutine var_io_write_double_1D

    subroutine var_io_write_double_2D(filename,v,dat,n,ncid,dims,missing_value,grid_mapping)

        implicit none

        character(len=*),  intent(IN) :: filename
        type(var_io_type), intent(IN) :: v
        real(8), contiguous, intent(IN) :: dat(:,:)
        integer,           intent(IN) :: n
        integer,           intent(IN), optional :: ncid
        character(len=*),  intent(IN), optional :: dims(:)
        real(8), intent(IN), optional :: missing_value
        character(len=*),  intent(IN), optional :: grid_mapping

        ! Local variables
        character(len=32), allocatable :: dims_out(:)
        integer, allocatable :: start(:), count(:)

        call var_io_write_dims(dims_out,start,count,v,shape(dat),n,dims)

        call nc_write(filename,trim(v%varname),dat,dims=dims_out,start=start,count=count, &
                      units=trim(v%units),long_name=trim(v%long_name), &
                      missing_value=missing_value,grid_mapping=grid_mapping,ncid=ncid, &
                      prec_abs=real(v%prec_abs,8),prec_rel=real(v%prec_rel,8))

        return

    end subroutine var_io_write_double_2D

    subroutine var_io_write_double_3D(filename,v,dat,n,ncid,dims,missing_value,grid_mapping)

        implicit none

        character(len=*),  intent(IN) :: filename
        type(var_io_type), intent(IN) :: v
        real(8), contiguous, intent(IN) :: dat(:,:,:)
        integer,           intent(IN) :: n
        integer,           intent(IN), optional :: ncid
        character(len=*),  intent(IN), optional :: dims(:)
        real(8), intent(IN), optional :: missing_value
        character(len=*),  intent(IN), optional :: grid_mapping

        ! Local variables
        character(len=32), allocatable :: dims_out(:)
        integer, allocatable :: start(:), count(:)

        call var_io_write_dims(dims_out,start,count,v,shape(dat),n,dims)

        call nc_write(filename,trim(v%varname),dat,dims=dims_out,start=start,count=count, &
                      units=trim(v%units),long_name=trim(v%long_name), &
                      missing_value=missing_value,grid_mapping=grid_mapping,ncid=ncid, &
                      prec_abs=real(v%prec_abs,8),prec_rel=real(v%prec_rel,8))

        return

    end subroutine var_io_write_double_3D

    subroutine var_io_write_int_1D(filename,v,dat,n,ncid,dims,missing_value,grid_mapping)

        implicit none

        character(len=*),  intent(IN) :: filename
        type(var_io_type), intent(IN) :: v
        integer, contiguous, intent(IN) :: dat(:)
        integer,           intent(IN) :: n
        integer,           intent(IN), optional :: ncid
        character(len=*),  intent(IN), optional :: dims(:)
        integer, intent(IN), optional :: missing_value
        character(len=*),  intent(IN), optional :: grid_mapping

        ! Local variables
        character(len=32), allocatable :: dims_out(:)
        integer, allocatable :: start(:), count(:)

        call var_io_write_dims(dims_out,start,count,v,shape(dat),n,dims)

        call nc_write(filename,trim(v%varname),dat,dims=dims_out,start=start,count=count, &
                      units=trim(v%units),long_name=trim(v%long_name), &
                      missing_value=missing_value,grid_mapping=grid_mapping,ncid=ncid)

        return

    end subroutine var_io_write_int_1D

    subroutine var_io_write_int_2D(filename,v,dat,n,ncid,dims,missing_value,grid_mapping)

        implicit none

        character(len=*),  intent(IN) :: filename
        type(var_io_type), intent(IN) :: v
        integer, contiguous, intent(IN) :: dat(:,:)
        integer,           intent(IN) :: n
        integer,           intent(IN), optional :: ncid
        character(len=*),  intent(IN), optional :: dims(:)
        integer, intent(IN), optional :: missing_value
        character(len=*),  intent(IN), optional :: grid_mapping

        ! Local variables
        character(len=32), allocatable :: dims_out(:)
        integer, allocatable :: start(:), count(:)

        call var_io_write_dims(dims_out,start,count,v,shape(dat),n,dims)

        call nc_write(filename,trim(v%varname),dat,dims=dims_out,start=start,count=count, &
                      units=trim(v%units),long_name=trim(v%long_name), &
                      missing_value=missing_value,grid_mapping=grid_mapping,ncid=ncid)

        return

    end subroutine var_io_write_int_2D

    subroutine var_io_write_int_3D(filename,v,dat,n,ncid,dims,missing_value,grid_mapping)

        implicit none

        character(len=*),  intent(IN) :: filename
        type(var_io_type), intent(IN) :: v
        integer, contiguous, intent(IN) :: dat(:,:,:)
        integer,           intent(IN) :: n
        integer,           intent(IN), optional :: ncid
        character(len=*),  intent(IN), optional :: dims(:)
        integer, intent(IN), optional :: missing_value
        character(len=*),  intent(IN), optional :: grid_mapping

        ! Local variables
        character(len=32), allocatable :: dims_out(:)
        integer, allocatable :: start(:), count(:)

        call var_io_write_dims(dims_out,start,count,v,shape(dat),n,dims)

        call nc_write(filename,trim(v%varname),dat,dims=dims_out,start=start,count=count, &
                      units=trim(v%units),long_name=trim(v%long_name), &
                      missing_value=missing_value,grid_mapping=grid_mapping,ncid=ncid)

        return

    end subroutine var_io_write_int_3D

    subroutine var_io_write_dims(dims_out,start,count,v,shp,n,dims)
        ! Dimensions, start and count of one time record of v, for var_io_write

        implicit none

        character(len=32), allocatable, intent(OUT) :: dims_out(:)
        integer, allocatable, intent(OUT) :: start(:)
        integer, allocatable, intent(OUT) :: count(:)
        type(var_io_type), intent(IN) :: v
        integer,           intent(IN) :: shp(:)         ! Shape of the data
        integer,           intent(IN) :: n              ! Time record
        character(len=*),  intent(IN), optional :: dims(:)

        ! Local variables
        integer :: nd

        nd = size(shp)

        allocate(dims_out(nd+1))

        if (present(dims)) then
            ! Explicit dimensions, incl. time
            if (size(dims) .ne. nd+1) then
                write(error_unit,*) "var_io_write:: Error: dims must hold one dimension per data &
                &dimension plus time."
                write(error_unit,*) "varname = ", trim(v%varname)
                write(error_unit,*) "rank(dat) = ", nd, ", size(dims) = ", size(dims)
                error stop 1
            end if
            dims_out = dims
        else
            ! Dimensions of the table, plus time
            if (v%ndims .ne. nd) then
                write(error_unit,*) "var_io_write:: Error: rank of the data does not match &
                &the dimensions of the variable table."
                write(error_unit,*) "varname = ", trim(v%varname)
                write(error_unit,*) "rank(dat) = ", nd, ", dimensions = ", trim(v%dimnames)
                error stop 1
            end if
            dims_out(1:nd) = v%dims
            dims_out(nd+1) = "time"
        end if

        allocate(start(nd+1),count(nd+1))
        start(1:nd)    = 1
        start(nd+1)    = n
        count(1:nd)    = shp
        count(nd+1)    = 1

        return

    end subroutine var_io_write_dims

end module variable_io
