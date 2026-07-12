module series
    ! Lightweight tabulated time-series reader for scalar forcing/index series.
    !
    ! A `series_class` holds an ordered time axis and one or more value channels
    ! (`nc`), interpolated in time on demand. It is the 0-D/1-channel cousin of
    ! `varslice` (which handles gridded x-y-z-time fields): use `series` for a
    ! scalar (nc=1) or multi-channel (nc>1, e.g. 12 monthly) forcing curve.
    !
    ! Two on-disk formats, auto-detected from the filename extension (override
    ! with the `fmt` argument):
    !   * ASCII  : whitespace-delimited [time  v1 [v2 ...]]; nc = ncols-1.
    !              Blank lines and lines beginning with '#' or '!' are skipped.
    !   * netCDF : a named value variable over a time coordinate. The variable
    !              may be 1-D (nc=1) or 2-D (nc = the non-time dimension). An
    !              optional per-time/per-channel standard deviation is read from
    !              a companion variable (see `sigma_name`).
    !
    ! Interpolation is clamped-linear (endpoints held outside the time range) and
    ! is agnostic to the meaning of `time`: whatever axis the caller passes must
    ! match the file's own axis (absolute year, year BP, elapsed time, ...).

    use precision, only: wp
    use constants, only: mv
    use ncio
    use nml
    use interp1D, only: interp_linear

    implicit none

    type series_class
        character(len=512)    :: filename
        integer               :: nt            ! number of time points
        integer               :: nc            ! number of value channels
        real(wp), allocatable :: time(:)       ! (nt)
        real(wp), allocatable :: var(:,:)      ! (nt,nc)
        real(wp), allocatable :: sigma(:,:)    ! (nt,nc); unallocated if absent
    end type

    private
    public :: series_class
    public :: series_load
    public :: series_init_nml
    public :: series_interp
    public :: series_interp1
    public :: series_interp_sig

contains

    subroutine series_load(ser,filename,varname,time_name,sigma_name,fmt)
        ! Load a time series from `filename`. For netCDF, `varname` selects the
        ! value variable, `time_name` the time coordinate (default "time"), and
        ! `sigma_name` an optional standard-deviation variable (read only if it
        ! exists in the file). All three are ignored for ASCII input.

        type(series_class), intent(INOUT) :: ser
        character(len=*),   intent(IN)    :: filename
        character(len=*),   intent(IN), optional :: varname
        character(len=*),   intent(IN), optional :: time_name
        character(len=*),   intent(IN), optional :: sigma_name
        character(len=*),   intent(IN), optional :: fmt

        character(len=56) :: file_fmt

        ! Reset any previous contents
        if (allocated(ser%time))  deallocate(ser%time)
        if (allocated(ser%var))   deallocate(ser%var)
        if (allocated(ser%sigma)) deallocate(ser%sigma)

        ser%filename = trim(filename)

        ! Determine the format
        if (present(fmt)) then
            file_fmt = trim(fmt)
        else
            file_fmt = detect_format(filename)
        end if

        select case(trim(file_fmt))
            case("ascii")
                call series_load_ascii(ser,filename)
            case("netcdf")
                if (.not. present(varname)) then
                    write(*,*) "series_load:: Error: netCDF input requires 'varname'."
                    write(*,*) "  filename = "//trim(filename)
                    stop
                end if
                call series_load_netcdf(ser,filename,varname,time_name,sigma_name)
            case default
                write(*,*) "series_load:: Error: unknown format '"//trim(file_fmt)//"'."
                stop
        end select

        return

    end subroutine series_load

    subroutine series_init_nml(ser,par_filename,group)
        ! Convenience initializer: read the file location and (netCDF) variable
        ! names from a namelist group, then load. Expected keys:
        !   filename, [varname], [time_name], [sigma_name]
        ! (varname/time_name/sigma_name are only consulted for netCDF files.)

        type(series_class), intent(INOUT) :: ser
        character(len=*),   intent(IN)    :: par_filename
        character(len=*),   intent(IN)    :: group

        character(len=512) :: filename
        character(len=56)  :: varname, time_name, sigma_name

        call nml_read(par_filename,group,"filename",filename)

        if (trim(detect_format(filename)) .eq. "netcdf") then
            call nml_read(par_filename,group,"varname",   varname)
            call nml_read(par_filename,group,"time_name", time_name)
            call nml_read(par_filename,group,"sigma_name",sigma_name)
            call series_load(ser,filename,varname=varname, &
                                time_name=time_name,sigma_name=sigma_name)
        else
            call series_load(ser,filename)
        end if

        return

    end subroutine series_init_nml

    function series_interp(ser,time) result(v)
        ! Interpolate all channels to `time`. Returns a length-nc vector.

        type(series_class), intent(IN) :: ser
        real(wp),           intent(IN) :: time
        real(wp) :: v(ser%nc)

        integer :: c

        do c = 1, ser%nc
            v(c) = interp_linear(ser%time,ser%var(:,c),xout=time)
        end do

        return

    end function series_interp

    function series_interp1(ser,time) result(f)
        ! Scalar convenience for single-channel series (nc==1).

        type(series_class), intent(IN) :: ser
        real(wp),           intent(IN) :: time
        real(wp) :: f

        if (ser%nc .ne. 1) then
            write(*,*) "series_interp1:: Error: series has nc = ", ser%nc, &
                        " (expected 1). Use series_interp for multi-channel series."
            stop
        end if

        f = interp_linear(ser%time,ser%var(:,1),xout=time)

        return

    end function series_interp1

    function series_interp_sig(ser,time) result(s)
        ! Interpolate the per-channel standard deviation to `time`. Returns zeros
        ! if the series carries no sigma information.

        type(series_class), intent(IN) :: ser
        real(wp),           intent(IN) :: time
        real(wp) :: s(ser%nc)

        integer :: c

        if (.not. allocated(ser%sigma)) then
            s = 0.0_wp
            return
        end if

        do c = 1, ser%nc
            s(c) = interp_linear(ser%time,ser%sigma(:,c),xout=time)
        end do

        return

    end function series_interp_sig

    ! === Internals ============================================================

    function detect_format(filename) result(fmt)
        ! Map a filename extension to "ascii" or "netcdf".

        character(len=*), intent(IN) :: filename
        character(len=56) :: fmt

        integer :: idot

        idot = index(filename,".",back=.TRUE.)

        fmt = "ascii"
        if (idot .gt. 0) then
            select case(trim(adjustl(filename(idot+1:))))
                case("nc","nc4","cdf","netcdf")
                    fmt = "netcdf"
                case default
                    fmt = "ascii"
            end select
        end if

        return

    end function detect_format

    subroutine series_load_ascii(ser,filename)
        ! Read a whitespace-delimited [time  v1 [v2 ...]] table. The channel
        ! count is inferred from the first data row.

        type(series_class), intent(INOUT) :: ser
        character(len=*),   intent(IN)    :: filename

        integer :: u, ios, i, ncols
        character(len=4096)   :: line
        real(wp), allocatable :: row(:)

        open(newunit=u,file=filename,status="old",action="read",iostat=ios)
        if (ios .ne. 0) then
            write(*,*) "series_load_ascii:: Error: cannot open file "//trim(filename)
            stop
        end if

        ! Pass 1: count data rows and columns
        ser%nt = 0
        ncols  = 0
        do
            read(u,"(a)",iostat=ios) line
            if (ios .ne. 0) exit
            if (is_blank_or_comment(line) .or. .not. is_numeric_start(line)) cycle
            if (ncols .eq. 0) ncols = count_tokens(line)
            ser%nt = ser%nt + 1
        end do

        if (ncols .lt. 2) then
            write(*,*) "series_load_ascii:: Error: need >=2 columns [time value...] in "//trim(filename)
            stop
        end if
        if (ser%nt .lt. 2) then
            write(*,*) "series_load_ascii:: Error: need >=2 data rows in "//trim(filename)
            stop
        end if

        ser%nc = ncols - 1
        allocate(ser%time(ser%nt))
        allocate(ser%var(ser%nt,ser%nc))
        allocate(row(ncols))

        ! Pass 2: parse
        rewind(u)
        i = 0
        do
            read(u,"(a)",iostat=ios) line
            if (ios .ne. 0) exit
            if (is_blank_or_comment(line) .or. .not. is_numeric_start(line)) cycle
            i = i + 1
            read(line,*,iostat=ios) row(1:ncols)
            if (ios .ne. 0) then
                write(*,*) "series_load_ascii:: Error: could not parse row ", i, &
                            " of "//trim(filename)
                stop
            end if
            ser%time(i)  = row(1)
            ser%var(i,:) = row(2:ncols)
        end do

        close(u)
        deallocate(row)

        return

    end subroutine series_load_ascii

    subroutine series_load_netcdf(ser,filename,varname,time_name,sigma_name)
        ! Read a value variable (1-D or 2-D) over a time coordinate, plus an
        ! optional standard-deviation variable if present in the file.

        type(series_class), intent(INOUT) :: ser
        character(len=*),   intent(IN)    :: filename
        character(len=*),   intent(IN)    :: varname
        character(len=*),   intent(IN), optional :: time_name
        character(len=*),   intent(IN), optional :: sigma_name

        character(len=56) :: tname

        tname = "time"
        if (present(time_name)) then
            if (len_trim(time_name) .gt. 0) tname = trim(time_name)
        end if

        ! Time axis
        ser%nt = nc_size(filename,tname)
        allocate(ser%time(ser%nt))
        call nc_read(filename,tname,ser%time)

        ! Value variable
        call read_var_channels(filename,varname,tname,ser%nt,ser%var)
        ser%nc = size(ser%var,2)

        ! Optional standard deviation (read only if the variable exists)
        if (present(sigma_name)) then
            if (len_trim(sigma_name) .gt. 0) then
                if (nc_exists_var(filename,trim(sigma_name))) then
                    call read_var_channels(filename,sigma_name,tname,ser%nt,ser%sigma)
                    if (size(ser%sigma,2) .ne. ser%nc) then
                        write(*,*) "series_load_netcdf:: Error: sigma channel count /= var channel count in " &
                                    //trim(filename)
                        stop
                    end if
                end if
            end if
        end if

        return

    end subroutine series_load_netcdf

    subroutine read_var_channels(filename,varname,tname,nt,out)
        ! Read a 1-D (time) or 2-D (time x channel, either order) netCDF variable
        ! into a canonical (nt,nc) array with time as the first index.

        character(len=*), intent(IN) :: filename
        character(len=*), intent(IN) :: varname
        character(len=*), intent(IN) :: tname
        integer,          intent(IN) :: nt
        real(wp), allocatable, intent(OUT) :: out(:,:)

        character(len=64), allocatable :: dnames(:)
        integer,           allocatable :: dsizes(:)
        real(wp),          allocatable :: buf(:,:)
        integer :: ndim, it, nc

        call nc_dims(filename,varname,names=dnames,dims=dsizes)
        ndim = size(dsizes)

        select case(ndim)

            case(1)
                nc = 1
                allocate(out(nt,nc))
                call nc_read(filename,varname,out(:,1))

            case(2)
                ! Identify the time axis: by name first, else by matching length
                it = 0
                if (trim(dnames(1)) .eq. trim(tname)) it = 1
                if (trim(dnames(2)) .eq. trim(tname)) it = 2
                if (it .eq. 0) then
                    if (dsizes(1) .eq. nt) it = 1
                    if (dsizes(2) .eq. nt) it = 2
                end if
                if (it .eq. 0) then
                    write(*,*) "read_var_channels:: Error: cannot find time axis of '"// &
                                trim(varname)//"' in "//trim(filename)
                    stop
                end if

                allocate(buf(dsizes(1),dsizes(2)))
                call nc_read(filename,varname,buf)

                if (it .eq. 1) then
                    nc = dsizes(2)
                    allocate(out(nt,nc))
                    out = buf
                else
                    nc = dsizes(1)
                    allocate(out(nt,nc))
                    out = transpose(buf)
                end if
                deallocate(buf)

            case default
                write(*,*) "read_var_channels:: Error: variable '"//trim(varname)// &
                            "' must be 1-D or 2-D in "//trim(filename)
                stop

        end select

        return

    end subroutine read_var_channels

    logical function is_blank_or_comment(line)
        ! .TRUE. for an empty/whitespace line or one starting with '#' or '!'.

        character(len=*), intent(IN) :: line
        character(len=:), allocatable :: s

        s = trim(adjustl(line))
        is_blank_or_comment = (len(s) .eq. 0)
        if (.not. is_blank_or_comment) then
            is_blank_or_comment = (s(1:1) .eq. "#" .or. s(1:1) .eq. "!")
        end if

        return

    end function is_blank_or_comment

    logical function is_numeric_start(line)
        ! .TRUE. if the first whitespace-delimited token parses as a real number.
        ! Used to skip a plain (non-comment) header row, e.g. "time  value", so
        ! ascii series files need not prefix their header with '#'/'!'.

        character(len=*), intent(IN) :: line
        character(len=:), allocatable :: s
        real(wp) :: x
        integer  :: ios

        s = trim(adjustl(line))
        is_numeric_start = .FALSE.
        if (len(s) .eq. 0) return
        read(s,*,iostat=ios) x
        is_numeric_start = (ios .eq. 0)

        return

    end function is_numeric_start

    integer function count_tokens(line)
        ! Count whitespace-delimited tokens (spaces and tabs) in a line.

        character(len=*), intent(IN) :: line

        integer :: i, n
        logical :: in_tok
        character(len=1) :: c
        character(len=1), parameter :: TAB = char(9)

        n      = 0
        in_tok = .FALSE.
        do i = 1, len_trim(line)
            c = line(i:i)
            if (c .eq. " " .or. c .eq. TAB) then
                in_tok = .FALSE.
            else
                if (.not. in_tok) n = n + 1
                in_tok = .TRUE.
            end if
        end do

        count_tokens = n
        return

    end function count_tokens

end module series
