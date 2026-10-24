module varslice

    use, intrinsic :: iso_fortran_env, only : input_unit, output_unit, error_unit
    
    use precision
    use constants, only: mv, TOL
    use ncio
    use nml
    use mapping, only : map_class, map_field, map_init
    use coordinates, only : grid_class, points_class, grid_init, points_init
    use interp2D,    only : fill_nearest_levels
    use fesmdata,    only : fesmdata_grid_desc_exists, fesmdata_grid_read
    use fesm_log,    only : fesm_error, to_str

    implicit none

    logical, parameter :: verbose = .FALSE.

    ! Values with magnitude >= var_lim_crazy are treated as unflagged
    ! fill/garbage and reset to the missing value. This catches NetCDF
    ! default fill values (~9.97e36) that were not tagged via missing_value.
    ! Note: legitimate data with magnitude above this limit would also be
    ! discarded — raise it if a variable can genuinely exceed 1e10.
    real(wp), parameter :: var_lim_crazy = 1.0e10_wp

    type varslice_param_class

        character(len=1024) :: filename
        character(len=1024), allocatable :: filenames(:)
        character(len=56)   :: name 
        character(len=56)   :: units_in
        character(len=56)   :: units_out
        real(wp) :: unit_scale 
        real(wp) :: unit_offset
        logical  :: with_time 
        logical  :: with_time_sub

        ! Online remapping (optional group keys; remap = "none" reads the file
        ! as it is): remap = kernel ("nn", "shepard", "quadrant", "bilinear",
        ! "con") onto the target grid given at init; grid_src = name of the
        ! file's grid (the map cache key).
        character(len=32)   :: remap
        character(len=256)  :: grid_src

        ! Filling of missing values (optional group keys; fill = "none" leaves
        ! them): fill = "nearest" (value of the nearest valid cell, per
        ! horizontal slice) or "const" (fill_value).
        character(len=32)   :: fill = "none"
        real(wp)            :: fill_value = mv

        ! What the variable is, resolved at init: "static" (no time axis),
        ! "clim" (a sub-annual cycle without years, the last dimension) or
        ! "transient" (a time axis). A legacy `time` key fixes it; the group
        ! key or argument `mode` (auto | static | clim | transient) must agree.
        character(len=16)   :: mode = ""

        ! Calendar year of time zero of the caller's clock: an axis decoded
        ! from the file in calendar years y is served as t = y - time_ref.
        real(dp)            :: time_ref = 0.0_dp

        ! Time axis from the file (time(1) < 0, or no time key): the time
        ! coordinate (the variable's last dimension, time_name) of each file,
        ! decoded with its units and calendar attributes, or with these when
        ! set (time_units = "clock": the values are clock years as they stand).
        logical             :: time_from_file = .false.
        character(len=128)  :: time_units = ""
        character(len=32)   :: calendar   = ""
        character(len=56)   :: time_name  = ""

        ! The units attribute of the variable in the (first) file; "" without
        ! one. Metadata: varslice converts nothing (scaling does), a caller
        ! such as fesm_datasrc checks it against the group's units.
        character(len=128)  :: units_file = ""

        ! Internal parameters
        integer  :: ndim 
        real(dp) :: time_par(4) 

    end type

    type varslice_class 

        type(varslice_param_class) :: par 

        ! Parameters defined during update call
        real(dp)          :: time_range(2)
        character(len=56) :: slice_method
        integer           :: range_rep
        character(len=32) :: fill_method = ""
        logical           :: clim_collapsed = .FALSE.   ! climatology read as the mean of its cycle
        
        ! Variable information
        integer,  allocatable :: dim(:)  
        real(wp), allocatable :: x(:) 
        real(wp), allocatable :: y(:)
        real(wp), allocatable :: z(:)  
        real(dp), allocatable :: time(:)
        real(dp), allocatable :: time_sub(:)
        real(dp), allocatable :: time_bnds(:,:)   ! (2,nt): CF bounds of a file axis, if any
        integer, allocatable  :: idx(:)
        integer, allocatable  :: nt_files(:)   ! per-file time length (cached at init)

        real(wp), allocatable :: var(:,:,:,:)

        ! Online remapping: the map from the file's grid (or, for curvilinear
        ! grids, its point set) to the target grid. var, x and y are on the
        ! target grid; dim stays the file's.
        type(map_class) :: map
        logical         :: map_points = .false.     ! source is a point set (2D lon/lat)

    contains

        ! Convenience accessors returning a natural-rank copy of the data,
        ! indexed by the trailing time/level position (default 1). The classic
        ! vs%var(:,:,k,1) indexing continues to work unchanged.
        !   v1(k) => var(:,k,1,1)   v2(k) => var(:,:,k,1)   v3(k) => var(:,:,:,k)
        procedure :: v1 => varslice_v1
        procedure :: v2 => varslice_v2
        procedure :: v3 => varslice_v3

    end type

    interface axis_init
        module procedure axis_init_sp
        module procedure axis_init_dp
    end interface

    ! Copy of the current slice at its natural rank, dp or sp, at trailing
    ! index k (default 1), into the caller's array (shape checked):
    !   0-D var(k,1,1,1)  1-D var(:,k,1,1)  2-D var(:,:,k,1)  3-D var(:,:,:,k)
    ! Storage stays real(wp); missing values come out as mv.
    interface varslice_get
        module procedure varslice_get_0d_dp, varslice_get_0d_sp
        module procedure varslice_get_1d_dp, varslice_get_1d_sp
        module procedure varslice_get_2d_dp, varslice_get_2d_sp
        module procedure varslice_get_3d_dp, varslice_get_3d_sp
    end interface

    private
    public :: varslice_param_class
    public :: varslice_class
    public :: varslice_update
    public :: varslice_nsub
    public :: varslice_sub_mean
    public :: varslice_init_nml
    public :: varslice_init_arg
    public :: varslice_init_data
    public :: varslice_end 
    public :: varslice_time_decode
    public :: varslice_get

    public :: varslice_map_to_grid

    public :: get_matching_files
    
    public :: print_var_range
    
contains

    subroutine varslice_map_to_grid(vs_tgt,vs_src,map,mask_out,stat,method,reset,missing_value, &
                                    mask_pack,fill_method,filt_method,filt_par,verbose)
        ! Given a input vs object on a source grid src,
        ! map the variable to a target grid tgt and store
        ! in new vs object.

        implicit none

        type(varslice_class),  intent(INOUT) :: vs_tgt
        type(varslice_class),  intent(IN)    :: vs_src
        type(map_class),       intent(IN)    :: map
        integer,                intent(OUT), optional :: mask_out(:,:)   ! Mask showing where interpolation was done
        character(len=*),       intent(IN),  optional :: stat            ! Aggregation: mean (default), count, stdev
        character(len=*),       intent(IN),  optional :: method          ! Interpolation-kernel override (distance maps)
        logical,                intent(IN),  optional :: reset           ! Reset var_tgt initially to missing_value?
        real(wp),               intent(IN),  optional :: missing_value   ! Points not included in mapping
        logical,                intent(IN),  optional :: mask_pack(:,:)  ! Mask for where to interpolate
        character(len=*),       intent(IN),  optional :: fill_method     ! Method to fill in remaining missing values
        character(len=*),       intent(IN),  optional :: filt_method     ! Method to use for filtering
        real(wp),               intent(IN),  optional :: filt_par(:)     ! gaussian=[sigma,dx]; poisson=[tol]
        logical,                intent(IN),  optional :: verbose         ! Print information
        
        ! Local variables
        integer :: nx, ny, nz, nt
        integer :: i, j, k, t
        integer :: ndim
        real(wp) :: miss                        ! missing value used for mapping
        logical, allocatable :: mask_int(:,:)   ! per-slice interpolation mask (logical form of mask_out)

        ! Missing value: caller override, else the package default
        miss = mv
        if (present(missing_value)) miss = missing_value

        ! Consistency check: is this at least a 2D field?
        ndim = size(vs_src%dim,1)

        if (ndim .lt. 2 .or. (vs_src%par%with_time .and. ndim .lt. 3)) then
            call fesm_error("varslice_map_to_grid", &
                "mapping can only be done for fields with two spatial dimensions.", &
                "name = "//trim(vs_src%par%name)//new_line("a")// &
                "ndim = "//to_str(ndim)//new_line("a")// &
                "dim  = "//to_str(vs_src%dim))
        end if

        ! Determine size of target var
        nx = map%G%nx
        ny = map%G%ny
        nz = size(vs_src%var,3)         ! z should keep the same dimension
        nt = size(vs_src%var,4)         ! time should keep the same dimension
        
        ! Initialize meta information by copying entire vs object
        vs_tgt = vs_src

        ! Re-allocate target varslice object fields as needed
        deallocate(vs_tgt%x)
        deallocate(vs_tgt%y)
        deallocate(vs_tgt%var)
        
        allocate(vs_tgt%x(nx))
        allocate(vs_tgt%y(ny))
        allocate(vs_tgt%var(nx,ny,nz,nt))

        ! Per-slice interpolation mask (logical); converted to mask_out below
        allocate(mask_int(nx,ny))

        ! Loop over dimensions and remap variable
        do t = 1, nt
        do k = 1, nz
            call map_field(map,trim(vs_tgt%par%name),vs_src%var(:,:,k,t),vs_tgt%var(:,:,k,t), &
                            stat=stat,method=method,missing_value=miss,reset=reset,mask_pack=mask_pack, &
                            fill_method=fill_method,filt_method=filt_method,filt_par=filt_par, &
                            verbose=verbose,mask2=mask_int)
        end do
        end do

        ! Report where interpolation was done (last slice), as an integer 0/1 mask
        if (present(mask_out)) mask_out = merge(1, 0, mask_int)

        ! Done! target field should now be available.

        return

    end subroutine varslice_map_to_grid

    subroutine varslice_remap_init(vs,grid)
        ! Prepare online remapping of the variable from the grid of its file
        ! onto `grid`. A projected grid grid_src with a cdo description
        ! (grid_<grid_src>.txt next to the file, as FesmData writes it, else in
        ! maps/) is read from it; otherwise the file's geometry comes from its
        ! lon/lat: those named in the variable's CF "coordinates" attribute,
        ! else its first two dimensions.
        !   - 1D lon/lat: a regular lon-lat grid; the map is cached in maps/
        !     under par%grid_src.
        !   - 2D lon/lat (curvilinear grids, e.g. tripolar ocean grids): a
        !     lon-lat point set.
        ! Afterwards vs%x, vs%y are the axes of `grid`.

        implicit none

        type(varslice_class), intent(INOUT) :: vs
        type(grid_class),     intent(IN)    :: grid

        ! Local variables
        character(len=1024) :: fname
        character(len=512)  :: coords_att
        character(len=56)   :: lon_name, lat_name, word
        character(len=12), allocatable :: dim_names(:)
        integer,  allocatable :: dims(:)
        real(dp), allocatable :: lon1(:), lat1(:), lon2(:,:), lat2(:,:)
        type(grid_class)   :: grid_src
        type(points_class) :: pts_src
        integer :: nx, ny, ndim_sp, i0, i1

        fname = vs%par%filenames(1)

        ! Number of spatial dimensions: 2 is required (the 3rd, if any, is the
        ! vertical, which is kept as it is)
        ndim_sp = vs%par%ndim
        if (vs%par%with_time) ndim_sp = ndim_sp - 1
        if (ndim_sp .lt. 2) then
            call fesm_error("varslice_remap_init", &
                "remapping needs a field with two horizontal dimensions.", &
                "name = "//trim(vs%par%name)//new_line("a")// &
                "ndim = "//to_str(vs%par%ndim))
        end if
        if (len_trim(vs%par%grid_src) .eq. 0) then
            call fesm_error("varslice_remap_init", &
                "grid_src (the name of the file's grid) must be set with remap.", &
                "name = "//trim(vs%par%name))
        end if

        if (fesmdata_grid_desc_exists(fname,vs%par%grid_src)) then
            ! Projected grid of its cdo description
            call fesmdata_grid_read(grid_src,fname,trim(vs%par%grid_src))
            call map_init(vs%map,grid_src,grid,method=trim(vs%par%remap),fldr="maps")
            vs%map_points = .false.
        else
            call remap_init_lonlat()
        end if

        ! The variable is returned on the target grid
        if (allocated(vs%x)) deallocate(vs%x)
        if (allocated(vs%y)) deallocate(vs%y)
        allocate(vs%x(grid%G%nx), vs%y(grid%G%ny))
        vs%x = grid%G%x
        vs%y = grid%G%y

        return

    contains

        subroutine remap_init_lonlat()
            ! The map from the lon/lat of the file.

            call nc_dims(fname,vs%par%name,dim_names,dims)
            nx = dims(1)
            ny = dims(2)

            coords_att = ""
            if (nc_exists_attr(fname,vs%par%name,"coordinates")) &
                call nc_read_attr(fname,vs%par%name,"coordinates",coords_att)

            ! lon/lat: named by the coordinates attribute (any order), else the
            ! variable's first two dimensions
            lon_name = dim_names(1)
            lat_name = dim_names(2)
            i0 = 1
            do while (i0 .le. len_trim(coords_att))
                i1 = index(coords_att(i0:)," ")
                if (i1 .eq. 0) i1 = len_trim(coords_att(i0:)) + 1
                word = coords_att(i0:i0+i1-2)
                if (len_trim(word) .gt. 0) then
                    if (index(word,"lon") .gt. 0) lon_name = word
                    if (index(word,"lat") .gt. 0) lat_name = word
                end if
                i0 = i0 + i1
            end do

            ! Their rank decides: 2D lon/lat is a curvilinear grid, 1D a regular one
            call nc_dims(fname,lon_name,dim_names,dims)

            if (size(dims) .eq. 2) then
                ! Curvilinear grid: a point set of its 2D lon/lat
                allocate(lon2(nx,ny), lat2(nx,ny))
                call nc_read(fname,lon_name,lon2)
                call nc_read(fname,lat_name,lat2)
                call points_init(pts_src,name=trim(vs%par%grid_src),mtype="latlon",units="degrees", &
                                 x=reshape(lon2,[nx*ny]),y=reshape(lat2,[nx*ny]))
                call map_init(vs%map,pts_src,grid,max_neighbors=10,method=trim(vs%par%remap))
                vs%map_points = .true.
            else
                ! Regular lon-lat grid: its 1D axes
                allocate(lon1(nx), lat1(ny))
                call nc_read(fname,lon_name,lon1)
                call nc_read(fname,lat_name,lat1)
                call grid_init(grid_src,name=trim(vs%par%grid_src),mtype="latlon",units="degrees", &
                               x=lon1,y=lat1)
                call map_init(vs%map,grid_src,grid,method=trim(vs%par%remap),fldr="maps")
                vs%map_points = .false.
            end if

        end subroutine remap_init_lonlat

    end subroutine varslice_remap_init

    subroutine varslice_remap(vs)
        ! Remap vs%var from the file's grid onto the target grid (the map of
        ! varslice_remap_init), slice by slice over the trailing (vertical /
        ! time) dimensions. Missing values are excluded from the kernel; target
        ! points without a valid source value are missing.

        implicit none

        type(varslice_class), intent(INOUT) :: vs

        ! Local variables
        real(wp), allocatable :: src(:,:,:,:)
        integer :: nxs, nys, n3, n4, k, t

        nxs = size(vs%var,1)
        nys = size(vs%var,2)
        n3  = size(vs%var,3)
        n4  = size(vs%var,4)

        call move_alloc(vs%var,src)
        allocate(vs%var(vs%map%G%nx,vs%map%G%ny,n3,n4))

        do t = 1, n4
        do k = 1, n3
            if (vs%map_points) then
                call map_field(vs%map,trim(vs%par%name),reshape(src(:,:,k,t),[nxs*nys]), &
                               vs%var(:,:,k,t),missing_value=mv,reset=.TRUE.)
            else
                call map_field(vs%map,trim(vs%par%name),src(:,:,k,t), &
                               vs%var(:,:,k,t),missing_value=mv,reset=.TRUE.)
            end if
        end do
        end do

        return

    end subroutine varslice_remap

    subroutine varslice_fill(vs,method)
        ! Fill the missing values of the slice: "nearest" gives each missing cell
        ! the value of the nearest valid cell, per horizontal (x,y) slice, and a
        ! level without data the values of the level above (fill_nearest_levels
        ! of interp2D); "const" sets them to par%fill_value.

        implicit none

        type(varslice_class), intent(INOUT) :: vs
        character(len=*),     intent(IN)    :: method

        ! Local variables
        integer :: t, ndim_sp

        select case(trim(method))

            case("none")
                ! Leave the missing values

            case("const")
                where (vs%var .eq. mv) vs%var = vs%par%fill_value

            case("nearest")
                ndim_sp = vs%par%ndim
                if (vs%par%with_time) ndim_sp = ndim_sp - 1
                if (ndim_sp .lt. 2 .and. trim(vs%par%remap) .eq. "none") then
                    call fesm_error("varslice_fill", &
                        "fill = ""nearest"" needs a variable with two horizontal dimensions.", &
                        "variable = "//trim(vs%par%name))
                end if

                do t = 1, size(vs%var,4)
                    if (any(vs%var(:,:,:,t) .eq. mv) .and. any(vs%var(:,:,:,t) .ne. mv)) then
                        call fill_nearest_levels(vs%var(:,:,:,t), mv)
                    end if
                end do

        end select

        return

    end subroutine varslice_fill

    subroutine check_fill_method(method,proc)

        implicit none

        character(len=*), intent(IN) :: method
        character(len=*), intent(IN) :: proc

        select case(trim(method))
            case("none","nearest","const")
                ! ok
            case DEFAULT
                call fesm_error(proc, "fill must be ""none"", ""nearest"" or ""const"".", &
                    "fill = "//trim(method))
        end select

    end subroutine check_fill_method

    subroutine varslice_update(vs,time,method,fill,rep,print_summary)
        ! Routine to update transient climate forcing to match 
        ! current `time`. 

        ! time = a specific time or a range of times

        ! method = ["exact","range","interp","extrap","range_min","range_mean","range_max"]
        ! difference time-slicing methods applied to the data read from the file.
        ! Check routine `get_indices()` to see how data is loaded for each slice_method.
        ! ["exact","range"] produce data identical to the input file, where
        ! "exact" loads the data from the file for the time index matching the desired time
        ! exactly - if it is not available, missing values are returned.
        ! "range" loads all the data within a certain time range.
        ! ["interp","extrap"] both return a time slice by interpolating to the
        ! desired time from the two closest bracketing timesteps available in the file. "extrap"
        ! allows for setting the time slice equal to the first or last timestep available, if
        ! the desired time is out of bounds, while "interp" returns missing values in this case. 
        ! "range_*" methods return one time slice with the method applied to the data within the 
        ! range given by `time`. The range must lie within the time axis of the file, else
        ! an error is raised.
        ! "interval_mean" returns the overlap-weighted mean over [t0,t1] of the file's
        ! piecewise representation: cell means over the CF time bounds (time_bnds) when the
        ! file has them, else the linear interpolant between the samples, with partial
        ! segments at both ends; the part of [t0,t1] outside the axis holds the first or
        ! last sample. With rep = nsub (> 1) each sub-annual step is averaged over the years,
        ! a year's step standing for [year, year+1). A single time t: the value at t.

        ! rep: frequency to apply slice_method over time. If rep=1, then calculation (mean/sd/etc)
        ! will be applied to each time index, returning a field with no time dimension. 
        ! If rep=12, then calculation will be applied to every 12th index, resulting in 
        ! 12 values along dimension.

        ! fill: method to fill the missing values of the slice ("none", "nearest", "const";
        ! see varslice_fill). Overrides the group's fill (par%fill).

        ! The slice is read once: a later call for the same time (range), method, rep
        ! and fill returns it as it is, and data without a time dimension is read only
        ! the first time. Callers must not modify vs%var.

        implicit none

        type(varslice_class), target, intent(INOUT) :: vs
        real(dp),         optional, intent(IN)    :: time(:)        ! [yr] Current time, or time range
        character(len=*), optional, intent(IN)    :: method         ! slice_method (only if with_time==True)
        character(len=*), optional, intent(IN)    :: fill           ! none, nearest, const (how to fill in missing values)
        integer,          optional, intent(IN)    :: rep            ! Only if with_time==True
        logical,          optional, intent(IN)    :: print_summary  ! Print summary of updated variable

        ! Local variables 
        integer :: k, k0, k1, nt
        integer :: nt_tot, nt_rep, nt_major, nt_out
        integer :: nsub
        integer :: n1, n2, i, j, l
        type(varslice_param_class) :: par 
        logical  :: with_time
        logical  :: with_time_sub
        real(dp) :: time_range(2) 
        real(wp), allocatable :: time_wt(:)
        character(len=56) :: slice_method
        character(len=56) :: fill_method
        character(len=56) :: vec_method 
        integer  :: range_rep
        integer,  allocatable :: kk(:)
        real(wp), allocatable, target :: var(:,:,:,:)
        real(wp), allocatable :: clim_mean(:,:,:)
        logical  :: is_clim, collapse
        integer  :: s, nsp
        real(wp), pointer, contiguous :: src2d(:,:), dst2d(:,:)

        ! Define shortcuts
        par = vs%par 
        with_time = par%with_time 

        if (present(rep)) then
            range_rep     = rep
        else
            range_rep     = 1 
        end if

        ! A climatology: a sub-annual cycle without years (time inactive, nsub > 1),
        ! its last dimension the cycle; a range mean with rep=1 collapses it
        is_clim  = (.not. with_time) .and. nint(par%time_par(4)) .gt. 1
        collapse = .FALSE.
        if (is_clim .and. present(method)) &
            collapse = (trim(method) .eq. "range_mean" .and. range_rep .eq. 1)

        if (range_rep .gt. 1) then
            with_time_sub = .TRUE.
        else
            with_time_sub = .FALSE.
        end if

        if (with_time) then 

            if (.not. present(time)) then
                call fesm_error("varslice_update", &
                    "current time or time range must be given as an argument (1D array).")
            end if

            ! Consistency check
            if (size(time,1) .eq. 2) then
                if (time(2) .lt. time(1)) then
                    call fesm_error("varslice_update", &
                        "time(2) should be >= time(1).", &
                        "time = "//to_str(time))
                end if
            end if

        end if 

        slice_method = "exact"
        if (present(method)) slice_method = trim(method)

        ! A range over a variable with a sub-annual axis selects whole years
        ! (all their sub-annual steps) also when it is collapsed to one field
        ! (rep=1), e.g. the annual mean of a monthly series over [y0,y1].
        if (par%with_time_sub .and. index(slice_method,"range") .eq. 1) with_time_sub = .TRUE.

        fill_method = par%fill
        if (present(fill)) fill_method = trim(fill)
        call check_fill_method(fill_method, "varslice_update")

        ! Validate that the requested repetition count is consistent with the
        ! variable's declared sub-annual period (time_par(4), i.e. varslice_nsub).
        ! For the reducing/interpolating methods, `rep` must either be 1 (collapse
        ! the whole range to a single field) or match the file's period (preserve
        ! the sub-annual cycle, e.g. rep=12 for monthly data). This catches
        ! mismatches such as requesting a 12-month climatology from annual data,
        ! or a typo'd rep that would silently reshape the time axis.
        if (with_time) then
            select case(trim(slice_method))
                case("interp","extrap","range_mean","range_sd","range_min","range_max","range_sum", &
                     "interval_mean")
                    nsub = max(1, nint(par%time_par(4)))
                    if (range_rep .ne. 1 .and. range_rep .ne. nsub) then
                        call fesm_error("varslice_update", &
                            "requested rep is inconsistent with the variable's declared "// &
                            "sub-annual period (time_par(4)). rep must be 1 or equal to the period.", &
                            "variable    = "//trim(par%name)//new_line("a")// &
                            "method      = "//trim(slice_method)//new_line("a")// &
                            "rep         = "//to_str(range_rep)//new_line("a")// &
                            "period      = "//to_str(nsub)//new_line("a")// &
                            "time_par(4) = "//to_str(par%time_par(4)))
                    end if
            end select
        end if

        if (present(time)) then

            if (size(time) .eq. 2) then 
                time_range = time
            else 
                time_range(1:2) = time(1) 
            end if 

        else 

            time_range = vs%time_range 

        end if 

        if ( allocated(vs%var) .and. trim(fill_method) .eq. trim(vs%fill_method) .and. &
             ( (.not. with_time .and. (.not. is_clim .or. collapse .eqv. vs%clim_collapsed)) .or. &
               ( trim(slice_method) .eq. trim(vs%slice_method) &
                .and. range_rep .eq. vs%range_rep &
                .and. vs%time_range(1) .eq. time_range(1) &
                .and. vs%time_range(2) .eq. time_range(2) ) ) ) then

            ! Do nothing, the varslice object is already up to date
            ! for the current time, method and fill (data without time: always,
            ! a climatology as long as it is (not) collapsed to its mean alike).

        else

            ! Set parameters in varslice object
            vs%slice_method = trim(slice_method)
            vs%range_rep    = range_rep
            vs%time_range   = time_range
            vs%fill_method  = trim(fill_method)

            ! 2. Read variable and convert units as needed

            ! Make sure var variable is deallocated and ready to be modified
            if (allocated(vs%var)) deallocate(vs%var)
            
            if (.not. with_time) then 
                ! Handle cases that do not have time dimension (simpler)

                select case(par%ndim)

                    case(1)

                        ! Allocate var to the right size 
                        allocate(vs%var(vs%dim(1),1,1,1))

                        ! 1D variable
                        call nc_read(par%filename,par%name,vs%var(:,1,1,1),missing_value=mv)

                    case(2)

                        ! Allocate var to the right size 
                        allocate(vs%var(vs%dim(1),vs%dim(2),1,1))

                        ! 2D variable 
                        call nc_read(par%filename,par%name,vs%var(:,:,1,1),missing_value=mv)

                    case(3)

                        ! Allocate var to the right size 
                        allocate(vs%var(vs%dim(1),vs%dim(2),vs%dim(3),1))

                        ! 3D variable 
                        call nc_read(par%filename,par%name,vs%var(:,:,:,1),missing_value=mv)

                    case(4)

                        ! 3D variable with its sub-annual cycle (a climatology:
                        ! varslice_init_data allows it only with nsub > 1)
                        allocate(vs%var(vs%dim(1),vs%dim(2),vs%dim(3),vs%dim(4)))
                        call nc_read(par%filename,par%name,vs%var,missing_value=mv)

                    case DEFAULT 

                        call fesm_error("varslice_update", &
                            "ndim > 4 with no time dimension is not allowed.", &
                            "ndim = "//to_str(par%ndim))

                end select

                ! A climatology (sub-annual cycle without years) asked for a
                ! range mean with rep=1: the mean over its cycle
                vs%clim_collapsed = collapse
                if (collapse) then
                    clim_mean = varslice_sub_mean(vs)
                    deallocate(vs%var)
                    allocate(vs%var(size(clim_mean,1),size(clim_mean,2),size(clim_mean,3),1))
                    vs%var(:,:,:,1) = clim_mean
                end if

            else if (trim(slice_method) .eq. "interval_mean") then
                ! Overlap-weighted mean over [t0,t1] of the file's piecewise
                ! representation (rep = 1 or nsub)

                call interval_mean_read(vs, time_range, range_rep)

            else 
                ! Cases with a time dimension (more complicated)

                if (trim(slice_method) .eq. "interp" .or. &
                    trim(slice_method) .eq. "extrap") then 
                    ! Update time range for interp/extrap methods 

                    ! Additional consistency check 
                    if (size(time,1) .ne. 1) then
                        call fesm_error("varslice_update", &
                            "to use slice_method=['interp','extrap'], only one time should be "// &
                            "provided as an argument.", &
                            "time = "//to_str(time))
                    end if

                end if

                ! Determine indices of data to load 

                call get_indices(vs%idx,vs%time,vs%time_range,slice_method,with_time_sub)
                
                k0 = minval(vs%idx)
                k1 = maxval(vs%idx)

                ! A reduction over a range that the file does not cover would
                ! return a field of missing values.
                if (k0 .le. 0 .and. index(slice_method,"range_") .eq. 1) then
                    call fesm_error("varslice_update", &
                        "time range is not covered by the time axis of the file.", &
                        "variable   = "//trim(par%name)//new_line("a")// &
                        "filename   = "//trim(par%filename)//new_line("a")// &
                        "method     = "//trim(slice_method)//new_line("a")// &
                        "time_range = "//to_str(vs%time_range)//new_line("a")// &
                        "time axis  = "//to_str([vs%time(1),vs%time(size(vs%time))]))
                end if

                !write(*,*) "idx: ", k0, k1
                !if (k0 .gt. 0) write(*,*) vs%time(k0)
                !if (k1 .gt. 0) write(*,*) vs%time(k1)

                if (k0 .gt. 0 .and. k1 .gt. 0) then 
                    ! Dimension range is available for loading, proceed 

                    ! Get size of time dimension needed for loading
                    nt_tot = k1-k0+1

                    ! Get size of time dimensions for major axis and sub axis
                    nt_rep   = vs%range_rep
                    nt_major = max(nt_tot / nt_rep, 1)

                    if (nt_major .ne. int(real(nt_tot)/real(nt_rep))) then
                        call fesm_error("varslice_update", &
                            "number of major time axis points does not match number of total "// &
                            "points divided by number of sub-time points.", &
                            "nt_rep   = "//to_str(nt_rep)//new_line("a")// &
                            "nt_tot   = "//to_str(nt_tot)//new_line("a")// &
                            "nt_major = "//to_str(nt_major)//new_line("a")// &
                            "int(real(nt_tot)/real(nt_rep)) = "//to_str(int(real(nt_tot)/real(nt_rep))))
                    end if

                    if (verbose) then
                        write(*,*) "nt_tot:   ", nt_tot
                        write(*,*) "nt_rep:   ", nt_rep
                        write(*,*) "nt_major: ", nt_major
                        ! write(*,*) "time: ", vs%time_range, k0, k1 
                        ! write(*,*) "      ", vs%time(k0), vs%time(k1)
                    end if

                    call read_slab(vs, k0, nt_tot, var)

                    ! At this point, the local var variable has been defined 
                    ! with data from the file for the appropriate time indices 

                    ! Next, we need to allocate the vs%var variable to the 
                    ! appropriate size and perform any calculations on the time 
                    ! indices of the local var variable as needed. 

                    ! Handle special case: if only one time is available 
                    ! for interp/extrap methods, then change method 
                    ! to exact 
                    if ( (trim(vs%slice_method) .eq. "interp" .or. & 
                          trim(vs%slice_method) .eq. "extrap") .and. &
                          nt_major .eq. 1) then 
                        ! Same time is given for upper and lower bound

                        slice_method = "exact"
                    end if

                    ! Now, allocate the vs%var variable to the right size

                    select case(trim(slice_method)) 

                        case("exact","range")
                            ! Allocate vs%var to the same size as var 
                            ! and store all values 

                            
                            if (allocated(vs%var)) deallocate(vs%var)
                            allocate(vs%var(size(var,1),size(var,2),size(var,3),size(var,4)))

                            ! Store data in vs%var 
                            vs%var = var 

                        case("interp","extrap","range_mean","range_sd","range_min","range_max","range_sum")
                            ! var should have two major time dimensions to interpolate
                            ! between. Allocate vs%var to appropriate size and 
                            ! perform interpolation 

                            ! Note: slice_method='interp' and 'extrap' use the same method, since 
                            ! the indices determine interpolation weights (ie, 
                            ! for slice_method='interp', if the time of interest lies
                            ! outside of the bounds of the data, then k0=k1=-1 and 
                            ! output data are set to missing values) 

                            select case(trim(slice_method))

                                case("interp","extrap")
                                    
                                    vec_method = "mean"     ! interp methods use the (weighted) mean

                                case("range_mean","range_sd","range_min","range_max","range_sum")

                                    ! Define 'vec_method'
                                    n1 = index(slice_method,"_")
                                    n2 = len_trim(slice_method)
                                    vec_method = slice_method(n1+1:n2)

                            end select

                            ! Size of dimension out is the size of the 
                            ! repitition desired. Ie, range_rep = 1 means 
                            ! to apply mean/sd/etc to all values along dimension
                            ! with the result of calculating 1 value. 
                            ! range_rep = 12 means apply calculation to every 12th 
                            ! value, resulting in 12 values along dimension.

                            nt_out = vs%range_rep

                            select case(trim(slice_method))

                                case("interp","extrap")
                                    
                                    if (nt_tot .ne. 2*nt_out) then
                                        ! Remember that if nt_tot==nt_out, this means that nt_major=1,
                                        ! and so only one time slice was available. So the method was
                                        ! changed to 'exact'.
                                        call fesm_error("varslice_update", &
                                            "something went wrong during interpolation. Exactly 2 major "// &
                                            "time slices should be available to interpolate between. Check!", &
                                            "nt_tot   = "//to_str(nt_tot)//new_line("a")// &
                                            "rep      = "//to_str(nt_out)//new_line("a")// &
                                            "2*rep    = "//to_str(2*nt_out)//new_line("a")// &
                                            "time     = "//to_str(time)//new_line("a")// &
                                            "k0, k1   = "//to_str([k0,k1])//new_line("a")// &
                                            "times    = "//to_str(vs%time(k0:k1)))
                                    end if

                                    ! Calculate time weighting between two extremes
                                    allocate(time_wt(2))
                                    if (with_time_sub) then
                                        time_wt(2) = real(floor(time(1))-floor(vs%time(k0))) / real(floor(vs%time(k1)) - floor(vs%time(k0)))
                                        time_wt(1) = 1.0_wp - time_wt(2) 
                                    else
                                        time_wt(2) = (time(1)-vs%time(k0)) / (vs%time(k1) - vs%time(k0))
                                        time_wt(1) = 1.0_wp - time_wt(2) 
                                    end if

                                case("range_mean","range_sd","range_min","range_max","range_sum")

                                    ! Equal weights to all values on major axis
                                    allocate(time_wt(nt_major))
                                    time_wt = 1.0

                            end select

                            if (minval(time_wt) .lt. 0.0_wp .or. maxval(time_wt) .gt. 1.0_wp) then
                                call fesm_error("varslice_update", &
                                    "interpolation weights are incorrect.", &
                                    "time_wt  = "//to_str(time_wt)//new_line("a")// &
                                    "time     = "//to_str(time)//new_line("a")// &
                                    "time(k0) = "//to_str(vs%time(k0))//new_line("a")// &
                                    "time(k1) = "//to_str(vs%time(k1)))
                            end if
                            
                            ! Make sure that var has at least as many values as we expect 
                            if (nt_out .gt. nt_tot) then
                                call fesm_error("varslice_update", &
                                    "the specified time range does not provide enough data points to "// &
                                    "be consistent with the specified value of range_rep "// &
                                    "(range_rep must be <= nt).", &
                                    "time_range      = "//to_str(vs%time_range)//new_line("a")// &
                                    "nt (time_range) = "//to_str(nt_tot)//new_line("a")// &
                                    "range_rep       = "//to_str(vs%range_rep))
                            end if

                            ! Allocate output at the natural rank (spatial dims
                            ! from the loaded slab, nt_out in the time slot).
                            call alloc_var(vs%var, shape(var), par%ndim, nt_out)

                            ! Collapse every rank to a single (space, time) view:
                            ! nsp is the flattened spatial extent, and the reduction
                            ! (mean/sd/interp/...) is one loop regardless of ndim.
                            nsp = size(var) / nt_tot
                            src2d(1:nsp,1:nt_tot) => var
                            dst2d(1:nsp,1:nt_out) => vs%var

                            do k = 1, nt_out

                                ! Get indices for current repitition
                                call get_rep_indices(kk,i0=k,i1=nt_tot,nrep=vs%range_rep)

                                do s = 1, nsp
                                    ! Calculate the vector value desired (mean,sd,...)
                                    call calc_vec_value(dst2d(s,k),src2d(s,kk),vec_method,mv,wt=time_wt)
                                end do

                            end do


                    end select

                else 
                    ! Dimension range was not found, set variable to missing values 

                    select case(par%ndim)
                        case(1)
                            allocate(vs%var(vs%dim(1),1,1,1))
                        case(2)
                            allocate(vs%var(vs%dim(1),vs%dim(2),1,1))
                        case(3)
                            allocate(vs%var(vs%dim(1),vs%dim(2),vs%dim(3),1))
                        case(4)
                            ! 3-D space + time: one (missing) time slice at natural rank
                            allocate(vs%var(vs%dim(1),vs%dim(2),vs%dim(3),1))
                    end select

                    vs%var = mv

                end if

            end if 

            ! Make sure crazy values have been set to missing (for safety)
            where (abs(vs%var) .ge. var_lim_crazy) vs%var = mv

            ! Apply scaling 
            where (vs%var .ne. mv) 
                vs%var = vs%var*par%unit_scale + par%unit_offset
            end where 

            ! Remap the slice onto the target grid
            if (trim(par%remap) .ne. "none") call varslice_remap(vs)

            ! Fill the missing values
            call varslice_fill(vs, fill_method)

        end if 

        if (present(print_summary)) then
            if (print_summary) then
                call print_var_range(vs%var, trim(slice_method), mv) 
            end if
        end if

        return 

    end subroutine varslice_update

    subroutine read_slab(vs, k0, nt_tot, var)
        ! Read time steps k0 .. k0+nt_tot-1 of the variable (over all files)
        ! into var, at its natural rank with the time in slot ndim.

        implicit none

        type(varslice_class),  intent(IN)  :: vs
        integer,               intent(IN)  :: k0, nt_tot
        real(wp), allocatable, intent(OUT) :: var(:,:,:,:)

        associate(par => vs%par)

            select case(par%ndim)

                case(1)
                    ! 0D (point) variable plus time dimension
                    allocate(var(nt_tot,1,1,1))
                    call nc_read_multifile(par%filenames,par%name,var,missing_value=mv, &
                            start=[k0],count=[nt_tot],nt_files=vs%nt_files)

                case(2)
                    ! 1D variable plus time dimension
                    allocate(var(vs%dim(1),nt_tot,1,1))
                    call nc_read_multifile(par%filenames,par%name,var,missing_value=mv, &
                            start=[1,k0],count=[vs%dim(1),nt_tot],nt_files=vs%nt_files)

                case(3)
                    ! 2D variable plus time dimension
                    allocate(var(vs%dim(1),vs%dim(2),nt_tot,1))
                    call nc_read_multifile(par%filenames,par%name,var,missing_value=mv, &
                            start=[1,1,k0],count=[vs%dim(1),vs%dim(2),nt_tot],nt_files=vs%nt_files)

                case(4)
                    ! 3D variable plus time dimension
                    allocate(var(vs%dim(1),vs%dim(2),vs%dim(3),nt_tot))
                    call nc_read_multifile(par%filenames,par%name,var,missing_value=mv, &
                            start=[1,1,1,k0],count=[vs%dim(1),vs%dim(2),vs%dim(3),nt_tot],nt_files=vs%nt_files)

                case DEFAULT
                    call fesm_error("varslice_update", &
                        "ndim > 4 with time dimension not allowed.", &
                        "ndim = "//to_str(par%ndim))

            end select

        end associate

    end subroutine read_slab

    subroutine interval_mean_read(vs, tr, rep)
        ! vs%var = the overlap-weighted mean over [tr(1),tr(2)] (see
        ! varslice_update), before scaling, remap and fill. rep = 1: one
        ! field; rep = nsub: one per sub-annual step.

        implicit none

        type(varslice_class), target, intent(INOUT) :: vs
        real(dp),             intent(IN)    :: tr(2)
        integer,              intent(IN)    :: rep

        real(dp), allocatable :: w(:), b(:,:), wg(:)
        integer,  allocatable :: grp(:), kk(:)
        real(wp), allocatable, target :: var(:,:,:,:)
        real(wp), allocatable :: wt(:)
        real(wp), pointer, contiguous :: src2d(:,:), dst2d(:,:)
        integer :: nt, k, k0, k1, o, nsp, s, nk

        nt = size(vs%time)
        allocate(w(nt), grp(nt))
        w = 0.0_dp

        if (rep .le. 1) then
            grp = 1
            if (allocated(vs%time_bnds)) then
                call cell_weights(vs%time_bnds, tr(1), tr(2), w)
            else
                call linear_weights(vs%time, tr(1), tr(2), w)
            end if
        else
            ! Per sub-annual step: the years as cells [y, y+1)
            grp = nint(vs%time_sub)
            do o = 1, rep
                kk = pack([(k, k=1,nt)], grp .eq. o)
                nk = size(kk)
                if (nk .eq. 0) cycle
                allocate(b(2,nk), wg(nk))
                b(1,:) = real(floor(vs%time(kk)), dp)
                b(2,:) = b(1,:) + 1.0_dp
                wg = 0.0_dp
                call cell_weights(b, tr(1), tr(2), wg)
                w(kk) = wg
                deallocate(b, wg)
            end do
        end if

        ! The steps with weight, read at once
        k0 = 0
        k1 = 0
        do k = 1, nt
            if (w(k) .gt. 0.0_dp) then
                if (k0 .eq. 0) k0 = k
                k1 = k
            end if
        end do
        if (k0 .eq. 0) then
            call fesm_error("varslice_update", "interval_mean: no time step has weight.", &
                "name  = "//trim(vs%par%name)//new_line("a")// &
                "range = "//to_str(tr))
        end if

        call read_slab(vs, k0, k1-k0+1, var)

        call alloc_var(vs%var, shape(var), vs%par%ndim, max(rep,1))
        nsp = size(var) / (k1-k0+1)
        src2d(1:nsp,1:k1-k0+1) => var
        dst2d(1:nsp,1:max(rep,1)) => vs%var

        do o = 1, max(rep,1)
            kk = pack([(k, k=1,k1-k0+1)], grp(k0:k1) .eq. o .and. w(k0:k1) .gt. 0.0_dp)
            wt = real(w(k0-1+kk), wp)
            do s = 1, nsp
                call calc_vec_value(dst2d(s,o), src2d(s,kk), "mean", mv, wt=wt)
            end do
        end do

        return

    end subroutine interval_mean_read

    subroutine linear_weights(t, t0, t1, w)
        ! Weights of the samples t in the mean over [t0,t1] of their linear
        ! interpolant, held constant outside the axis; t0 = t1: the
        ! interpolation weights at t0.

        implicit none

        real(dp), intent(IN)    :: t(:), t0, t1
        real(dp), intent(INOUT) :: w(:)

        integer  :: n, k
        real(dp) :: a, b, h, lam

        n = size(t)
        if (n .eq. 1) then
            w(1) = 1.0_dp
            return
        end if

        if (t1 .le. t0) then
            if (t0 .le. t(1)) then
                w(1) = 1.0_dp
            else if (t0 .ge. t(n)) then
                w(n) = 1.0_dp
            else
                do k = 1, n-1
                    if (t0 .ge. t(k) .and. t0 .lt. t(k+1)) then
                        lam = (t0 - t(k))/(t(k+1) - t(k))
                        w(k)   = w(k)   + (1.0_dp - lam)
                        w(k+1) = w(k+1) + lam
                        exit
                    end if
                end do
            end if
            return
        end if

        if (t0 .lt. t(1)) w(1) = w(1) + (min(t1, t(1)) - t0)
        if (t1 .gt. t(n)) w(n) = w(n) + (t1 - max(t0, t(n)))

        do k = 1, n-1
            a = max(t0, t(k))
            b = min(t1, t(k+1))
            if (b .gt. a) then
                h   = t(k+1) - t(k)
                lam = (0.5_dp*(a + b) - t(k))/h
                w(k)   = w(k)   + (b - a)*(1.0_dp - lam)
                w(k+1) = w(k+1) + (b - a)*lam
            end if
        end do

    end subroutine linear_weights

    subroutine cell_weights(bnd, t0, t1, w)
        ! Weights of cells [bnd(1,k), bnd(2,k)] (cell means) in the mean over
        ! [t0,t1]: their overlap; the part before the first or after the last
        ! cell holds that cell; t0 = t1, or an interval in a gap between
        ! cells: the cell containing, else nearest to, the midpoint.

        implicit none

        real(dp), intent(IN)    :: bnd(:,:), t0, t1
        real(dp), intent(INOUT) :: w(:)

        integer  :: n, k, kn
        real(dp) :: a, b, tm, d, dn
        logical  :: any_w

        n  = size(bnd,2)
        tm = 0.5_dp*(t0 + t1)
        any_w = .false.

        if (t1 .gt. t0) then
            if (t0 .lt. bnd(1,1)) then
                w(1) = w(1) + (min(t1, bnd(1,1)) - t0)
                any_w = .true.
            end if
            if (t1 .gt. bnd(2,n)) then
                w(n) = w(n) + (t1 - max(t0, bnd(2,n)))
                any_w = .true.
            end if
            do k = 1, n
                a = max(t0, bnd(1,k))
                b = min(t1, bnd(2,k))
                if (b .gt. a) then
                    w(k) = w(k) + (b - a)
                    any_w = .true.
                end if
            end do
        end if

        if (.not. any_w) then
            kn = 1
            dn = huge(1.0_dp)
            do k = 1, n
                if (tm .ge. bnd(1,k) .and. tm .lt. bnd(2,k)) then
                    kn = k
                    exit
                end if
                d = min(abs(tm - bnd(1,k)), abs(tm - bnd(2,k)))
                if (d .lt. dn) then
                    dn = d
                    kn = k
                end if
            end do
            w(kn) = 1.0_dp
        end if

    end subroutine cell_weights

    subroutine get_indices(idx, x, xrange, slice_method, with_sub)
        ! Get indices in range x0 <= x <= x1
        ! This way all decimal values within x0 and x1 are included.
        
        implicit none
        
        integer,  allocatable, intent(INOUT) :: idx(:)          ! Output indices
        real(dp),              intent(IN)    :: x(:)            ! Time array in years
        real(dp),              intent(IN)    :: xrange(:)       ! [Start, End] x (inclusive) or [x_interp]
        character(len=*),      intent(IN)    :: slice_method    ! method = "exact", "interp", "extrap"
        logical,               intent(IN)    :: with_sub        ! Use fractional time unit too

        ! Local variables
        integer :: i, n, nidx
        integer,  allocatable :: ii(:)
        real(dp), allocatable :: xmain(:)
        real(dp), allocatable :: dist(:)
        real(dp) :: dist_min_lo, dist_min_hi
        real(dp) :: x0, x1

        n     = size(x)

        ! Size work arrays to the time axis (no fixed upper bound), and
        ! keep them in double precision to match the incoming time axis x.
        allocate(ii(n))
        allocate(xmain(n))
        allocate(dist(n))

        xmain = 0
        nidx  = 0
        ii    = 0

        if (trim(slice_method) .eq. "interp" .or. trim(slice_method) .eq. "extrap") then
            x0 = xrange(1)
            x1 = xrange(1)
        else
            x0 = xrange(1)
            x1 = xrange(2) 
        end if

        if (with_sub) then
            if ( abs(x0-floor(x0)) .gt. TOL .or. abs(x1-floor(x1)) .gt. TOL) then
                call fesm_error("get_indices", &
                    "when time sub-axis is used, the time range should be specified by whole numbers.", &
                    "xrange = "//to_str(xrange))
            end if
        end if

        if (with_sub) then
            xmain(1:n) = floor(x)
        else
            xmain(1:n) = x
        end if

        select case(trim(slice_method))

            case("exact","range","range_mean","range_sd","range_min","range_max","range_sum")

                if ( x0 .ge. xmain(1)-TOL .and. x1 .le. xmain(n)+TOL) then
                    ! All values within the range are available, proceed
                    
                    do i = 1, n
                        if (xmain(i) >= x0-TOL .and. xmain(i) <= x1+TOL) ii(i) = i
                    end do
                    
                end if 

            case("interp","extrap")

                if (x0 .ge. xmain(1)-TOL .and. x0 .le. xmain(n)+TOL) then
                    ! Interpolation is possible, proceed
                    
                    dist = x0 - xmain(1:n)
                    dist_min_lo = maxval(dist(1:n),mask=dist(1:n).lt.0.0)
                    dist_min_hi = minval(dist(1:n),mask=dist(1:n).ge.0.0)

                    do i = 1, n
                        if (dist(i) .eq. dist_min_lo .or. dist(i) .eq. dist_min_hi) then
                            ii(i) = i
                        end if
                    end do
                
                else if (trim(slice_method) .eq. "extrap") then
                    ! Interp value is above or below bounds, and extrapolation is desired

                    if (x0 .lt. xmain(1)-TOL) then
                        dist = x0 - xmain(1:n)
                        dist_min_lo = maxval(dist(1:n),mask=dist(1:n).lt.0.0)
                        do i = 1, n
                            if (dist(i) .eq. dist_min_lo) then
                                ii(i) = i
                            end if
                        end do
                    else if (x0 .gt. xmain(n)+TOL) then
                        dist = x0 - xmain(1:n)
                        dist_min_hi = minval(dist(1:n),mask=dist(1:n).ge.0.0)
                        do i = 1, n
                            if (dist(i) .eq. dist_min_hi) then
                                ii(i) = i
                            end if
                        end do
                    end if

                end if

            case DEFAULT
                ! Guard against a slice_method that has no index rule here.
                ! Without this, an unhandled method silently yields no
                ! indices and the whole field is returned as missing.
                call fesm_error("get_indices", &
                    "slice_method not recognized.", &
                    "slice_method = "//trim(slice_method))

        end select

        nidx = count(ii .gt. 0)

        if (nidx .gt. 0) then
            if (allocated(idx)) deallocate(idx)
            allocate(idx(nidx)) 
            idx = pack(ii, ii.gt.0)
        else
            if (allocated(idx)) deallocate(idx)
            allocate(idx(1))
            idx(1) = -1
        end if

        if (verbose) then
            write(*,*) "get_indices: ", x0, x1, slice_method, with_sub
            !write(*,*) "x: ", xmain(1:n)
            write(*,*) "nidx: ", nidx
            !write(*,*) "idx: ", pack(ii,ii.gt.0)
            write(*,*) "idx: ", idx
            write(*,*) "x:   ", x(idx)
        end if

        return

    end subroutine get_indices

    subroutine get_rep_indices(ii,i0,i1,nrep)
        ! Given starting value i0 and final value i1, 
        ! and number of values to skip nrep, generate 
        ! a vector of indices ii. 
        ! eg, i0 = 1, i1 = 36, nrep = 12
        ! => ii = [1,13,25]
        ! eg, i0 = 2, i1 = 36, nrep = 12
        ! => ii = [2,14,26]
        
        implicit none 

        integer, allocatable, intent(OUT) :: ii(:) 
        integer, intent(IN) :: i0 
        integer, intent(IN) :: i1 
        integer, intent(IN) :: nrep 

        ! Local variables
        integer :: i, ni, ntot
        integer, allocatable :: jj(:)

        ni = i1-i0+1

        allocate(jj(ni))
        jj = 0

        do i = 1, ni
            jj(i) = i0 + nrep*(i-1) 
            if (jj(i) .gt. i1) then 
                jj(i) = 0
                exit 
            end if 
        end do 

        ntot = count(jj .gt. 0)
        if (allocated(ii)) deallocate(ii)
        allocate(ii(ntot)) 

        ii = jj(1:ntot) 

        return 

    end subroutine get_rep_indices

    subroutine calc_vec_value(val,var,method,mv,wt)

        implicit none 

        real(wp),         intent(OUT) :: val 
        real(wp),         intent(IN)  :: var(:) 
        character(len=*), intent(IN)  :: method 
        real(wp),         intent(IN)  :: mv 
        real(wp),         intent(IN), optional :: wt(:) 

        ! Local variables 
        integer  :: ntot 
        real(wp) :: mean, variance 
        real(wp), allocatable :: wt_now(:) 

        ntot = count(var .ne. mv)
        
        if (ntot .gt. 0) then 
            ! Values available for calculations 

            allocate(wt_now(size(var)))
            if (present(wt)) then 
                wt_now = wt

                ! Safety check
                if (size(wt,1) .ne. size(var,1)) then
                    call fesm_error("calc_vec_value", &
                        "wt vector must be the same length as the var vector.", &
                        "size(wt)  = "//to_str(size(wt,1))//new_line("a")// &
                        "size(var) = "//to_str(size(var,1)))
                end if
            else 
                wt_now = 1.0_wp
            end if 

            where(var .eq. mv) wt_now = 0.0_wp 

            select case(trim(method))

                case("sum","mean")
                    ! Calculate a weighted sum/mean with the right weights for each case 

                    ! Normalize weights to sum to one
                    if (trim(method) .eq. "mean") wt_now = wt_now / sum(wt_now)
                    
                    val = sum(wt_now*var)

                case("sd")
                    ! Calculate the weighted mean and then weighted standard deviation 

                    ! Normalize weights to sum to one
                    wt_now = wt_now / sum(wt_now)

                    if (ntot .ge. 2) then

                        mean     = sum(wt_now*var)
                        variance = real(ntot/(ntot-1),wp) * sum( wt_now*(var-mean)**2 )
                        val      = sqrt(variance)

                    else 

                        val = mv 

                    end if 

                case("min") 

                    val = minval(var,mask=var.ne.mv)

                case("max")

                    val = maxval(var,mask=var.ne.mv)

                case DEFAULT 

                    call fesm_error("calc_vec_value", &
                        "method not recognized.", &
                        "method = "//trim(method))

            end select

        else 
            ! No values available in vector, set to missing value 

            val = mv 

        end if 

        return 

    end subroutine calc_vec_value

    subroutine varslice_init_nml(vs,filename,group,domain,grid_name,verbose,subs,grid, &
                                 defaults_file,defaults_group,mode,remap,grid_src,time_ref)
        ! Routine to load information related to a given 
        ! transient variable, so that it can be processed properly.
        ! grid: the target grid. Required with remap /= "none"; without
        ! remapping the file must be on it (its two horizontal dimensions
        ! are checked against the grid's).
        !
        ! defaults_file: every key is read through the schema group
        ! defaults_group (default: group) of this file, so the group holds
        ! only what differs. Without it (the legacy path) filename, name,
        ! units, scaling and time are required, as always.
        ! mode, remap, grid_src, time_ref override the group's keys.
        ! remap = "auto": remap with "bilin" when the grid of the file
        ! (grid_src, else its global attribute grid_name) is known and is
        ! not `grid`; else read as it is.

        implicit none 

        type(varslice_class),   intent(INOUT) :: vs
        character(len=*),       intent(IN)    :: filename
        character(len=*),       intent(IN)    :: group
        character(len=*),       intent(IN), optional :: domain
        character(len=*),       intent(IN), optional :: grid_name
        logical,                intent(IN), optional :: verbose
        character(len=*),       intent(IN), optional :: subs(:,:)   ! extra {key}->value path substitutions
        type(grid_class),       intent(IN), optional :: grid
        character(len=*),       intent(IN), optional :: defaults_file    ! schema of the group's keys
        character(len=*),       intent(IN), optional :: defaults_group   ! its group (default: group)
        character(len=*),       intent(IN), optional :: mode       ! auto | static | clim | transient
        character(len=*),       intent(IN), optional :: remap      ! auto | none | a map method
        character(len=*),       intent(IN), optional :: grid_src   ! grid of the file
        real(dp),               intent(IN), optional :: time_ref   ! calendar year of time zero
        ! Local variables

        ! First load parameters from nml file
        call varslice_par_load(vs%par,filename,group,domain,grid_name,verbose,subs, &
                               defaults_file,defaults_group,mode,remap,grid_src,time_ref)

        ! remap = "auto": decided by the grid of the file and the target grid
        if (trim(vs%par%remap) .eq. "auto") call remap_auto(vs%par,grid)

        ! Perform remaining init operations 
        call varslice_init_data(vs) 

        if (trim(vs%par%remap) .ne. "none") then
            if (.not. present(grid)) then
                call fesm_error("varslice_init_nml", &
                    "remap is set, but no target grid was given.", &
                    "group = "//trim(group)//new_line("a")// &
                    "remap = "//trim(vs%par%remap))
            end if
            call varslice_remap_init(vs,grid)
        else if (present(grid)) then
            call varslice_check_grid(vs,grid,group)
        end if

        return 

    end subroutine varslice_init_nml

    subroutine remap_auto(par,grid)
        ! remap = "auto": "bilin" from the grid of the file (grid_src, else the
        ! file's global attribute grid_name) when it is known and is not the
        ! target grid; else "none" (read as it is, size checked).

        implicit none

        type(varslice_param_class), intent(INOUT) :: par
        type(grid_class), intent(IN), optional    :: grid

        character(len=256) :: gsrc

        gsrc = par%grid_src
        if (len_trim(gsrc) .eq. 0) then
            if (nc_exists_attr(par%filenames(1),"grid_name")) &
                call nc_read_attr(par%filenames(1),"grid_name",gsrc)
        end if

        par%remap = "none"
        if (present(grid)) then
            if (len_trim(gsrc) .gt. 0 .and. trim(gsrc) .ne. trim(grid%name)) then
                par%remap    = "bilin"
                par%grid_src = gsrc
            end if
        end if

        return

    end subroutine remap_auto

    subroutine varslice_check_grid(vs,grid,group)
        ! A variable read without remapping must be on the target grid: its
        ! first two dimensions (when it has two besides time) must have the
        ! grid's size, else the caller would get a block of a larger field
        ! (e.g. a corner of a finer grid) or an incomplete one.

        implicit none

        type(varslice_class), intent(IN) :: vs
        type(grid_class),     intent(IN) :: grid
        character(len=*),     intent(IN) :: group

        ! Local variables
        integer :: ndim_sp
        character(len=256) :: msg

        ndim_sp = vs%par%ndim
        if (vs%par%with_time) ndim_sp = ndim_sp - 1
        if (ndim_sp < 2) return

        if (vs%dim(1) /= grid%G%nx .or. vs%dim(2) /= grid%G%ny) then
            write(msg,"(a,i0,a,i0,a,i0,a,i0,a)") "file: ", vs%dim(1), " x ", vs%dim(2), &
                                               ", grid "//trim(grid%name)//": ", grid%G%nx, " x ", grid%G%ny, "."
            call fesm_error("varslice_init_nml", &
                "the variable is not on the target grid (set remap and grid_src to remap it).", &
                "group = "//trim(group)//new_line("a")// &
                "file  = "//trim(vs%par%filenames(1))//new_line("a")//trim(msg))
        end if

        return

    end subroutine varslice_check_grid

    subroutine varslice_init_arg(vs,filename,name,units_in,units_out,scale,offset, &
                                                            with_time,time_par)
        ! Routine to load information related to a given 
        ! transient variable, so that it can be processed properly.

        implicit none 

        type(varslice_class), intent(INOUT) :: vs
        character(len=*),      intent(IN)   :: filename
        character(len=*),      intent(IN)   :: name
        character(len=*),      intent(IN)   :: units_in 
        character(len=*),      intent(IN)   :: units_out 
        real(wp), optional,    intent(IN)   :: scale 
        real(wp), optional,    intent(IN)   :: offset 
        logical,  optional,    intent(IN)   :: with_time 
        real(dp), optional,    intent(IN)   :: time_par(:) 

        ! Local variables 
        
        ! Define parameters from subroutine arguments 
        vs%par%filename  = trim(filename) 
        vs%par%name      = trim(name) 
        vs%par%units_in  = trim(units_in) 
        vs%par%units_out = trim(units_out) 
        
        vs%par%unit_scale = 1.0_wp 
        if (present(scale)) vs%par%unit_scale = scale 

        vs%par%unit_offset = 0.0_wp 
        if (present(offset)) vs%par%unit_offset = offset 
        
        vs%par%with_time = .TRUE. 
        if (present(with_time)) vs%par%with_time = with_time 

        vs%par%time_par = [0.0,0.0,0.0,0.0]
        if (present(time_par)) vs%par%time_par(1:size(time_par)) = time_par

        vs%par%remap    = "none"
        vs%par%grid_src = ""

        vs%par%fill       = "none"
        vs%par%fill_value = mv

        vs%par%time_ref = 0.0_dp
        if (vs%par%with_time) then
            vs%par%mode = "transient"
        else if (nint(vs%par%time_par(4)) .gt. 1) then
            vs%par%mode = "clim"
        else
            vs%par%mode = "static"
        end if

        ! Resolve file list and derive internal time parameters
        ! (mirrors the nml path; without this par%filenames and
        !  par%with_time_sub would be left unset)
        call varslice_par_finalize(vs%par)

        ! Perform remaining init operations
        call varslice_init_data(vs)

        return

    end subroutine varslice_init_arg

    subroutine varslice_init_data(vs)

        implicit none 

        type(varslice_class), intent(INOUT) :: vs 

        ! Local variables  
        character(len=12), allocatable :: dim_names(:) 
        logical :: with_time 
        logical :: with_time_sub
        real(dp) :: dt 
        integer  :: i, k 
        integer  :: nt
        integer, allocatable :: dim_now(:)
        character(len=1024)  :: filename_dims
        character(len=:), allocatable :: fnames
        real(dp), allocatable :: t_file(:), tb_file(:,:)

        ! A time axis from the file: resolve the mode and read the axis
        if (vs%par%time_from_file) call file_time_axis(vs%par, t_file, tb_file)

        ! Local shortcut
        with_time       = vs%par%with_time 
        with_time_sub   = vs%par%with_time_sub 

        ! Store the first filename locally, to be used for getting variable dimensions (except time)
        filename_dims = vs%par%filenames(1)

        ! First make sure all data objects are deallocated 
        call varslice_end(vs)

        ! Get information from netcdf file - first file in case their are multiple sources 
        call nc_dims(filename_dims,vs%par%name,dim_names,vs%dim)
        vs%par%ndim = size(vs%dim,1)

        ! Units of the variable in the file (metadata)
        vs%par%units_file = ""
        if (nc_exists_attr(filename_dims,vs%par%name,"units")) &
            call nc_read_attr(filename_dims,vs%par%name,"units",vs%par%units_file)

        ! Determine total time dimension size and cache the per-file time
        ! lengths. nt_files is reused by nc_read_multifile so it need not
        ! re-query the files on every update.
        if (with_time) then
            if (allocated(vs%nt_files)) deallocate(vs%nt_files)
            allocate(vs%nt_files(size(vs%par%filenames,1)))
            do i = 1, size(vs%par%filenames,1)
                call nc_dims(vs%par%filenames(i),vs%par%name,dim_names,dim_now)
                vs%nt_files(i) = dim_now(vs%par%ndim)
            end do
            vs%dim(vs%par%ndim) = sum(vs%nt_files)
        end if

        if (with_time .and. vs%par%time_from_file) then

            ! The axis of the file(s); with a sub-annual cycle, the position
            ! of each step within its year
            vs%time = t_file
            if (allocated(tb_file)) vs%time_bnds = tb_file
            allocate(vs%time_sub(size(vs%time)))
            vs%time_sub = 1.0_dp
            if (with_time_sub) then
                do k = 2, size(vs%time)
                    if (floor(vs%time(k)) .eq. floor(vs%time(k-1))) vs%time_sub(k) = vs%time_sub(k-1) + 1.0_dp
                end do
            end if

        else if (with_time) then

            if (with_time_sub) then
                ! Need to generate a sub-annual axis too

                dt = 1.0 / vs%par%time_par(4)

                call axis_init(vs%time, &
                                x0=vs%par%time_par(1) + dt/2.0, &
                                x1=vs%par%time_par(2) + 1.0 - dt/2.0, &
                                dx=dt)

                ! Store number of sub-value (1:time_par(4)), so 1:12 for monthly values, repeating for all of time.

                if (allocated(vs%time_sub)) deallocate(vs%time_sub)
                allocate(vs%time_sub(size(vs%time)))

                do k = 1, size(vs%time)
                    vs%time_sub(k) = ceiling( (vs%time(k)-floor(vs%time(k)))/dt )
                    !write(*,*) vs%time(k), vs%time_sub(k)
                end do

                ! write(*,*) "Testing axis generation: "
                ! write(*,*) 
                ! do k = 1, size(vs%time)
                !     write(*,*) vs%time(k), floor(vs%time(k)), ceiling( (vs%time(k)-floor(vs%time(k)))/dt )
                ! end do
                ! write(*,*) 
                ! write(*,*) "nt = ", size(vs%time)
                ! stop 

                vs%par%with_time_sub = .TRUE. 

            else
                ! No special sub-annual axis needed.
                ! Note: time_par(4) must be greater than 1 to work
                ! as a fractional amount of time between major units.

                ! Initialize time vector from user parameters 
                call axis_init(vs%time,x0=vs%par%time_par(1), &
                                    x1=vs%par%time_par(2), &
                                    dx=vs%par%time_par(3))

            end if

            ! Check to make sure the netcdf file covers the time vector. A longer
            ! file is allowed: the time vector describes its leading steps (e.g.
            ! a CMIP experiment that one model extended beyond the protocol).
            if (size(vs%time,1) .lt. vs%dim(vs%par%ndim)) then
                write(*,*) "varslice_init_data:: note: ", trim(vs%par%name), ": the file has ", &
                    vs%dim(vs%par%ndim), " time steps; the first ", size(vs%time,1), &
                    " (time_par) are used."
            else if (size(vs%time,1) .gt. vs%dim(vs%par%ndim)) then
                fnames = ""
                if (size(vs%par%filenames,1) .gt. 1) then
                    fnames = "filenames    ="
                    do i = 1, size(vs%par%filenames,1)
                        fnames = fnames//new_line("a")//"               "//trim(vs%par%filenames(i))
                    end do
                    fnames = fnames//new_line("a")
                end if
                call fesm_error("varslice_init_data", &
                    "generated time coordinate is longer than the time dimension in the "// &
                    "netcdf file.", &
                    "time_par    = "//to_str(vs%par%time_par)//new_line("a")// &
                    "size(time)  = "//to_str(size(vs%time,1))//new_line("a")// &
                    "nt (netcdf) = "//to_str(vs%dim(vs%par%ndim))//new_line("a")// &
                    "filename    = "//trim(vs%par%filename)//new_line("a")// &
                    fnames// &
                    "time        = "//to_str(vs%time))
            end if

        end if 

        ! Allocate coordinate variables and load coordinates too. The number
        ! of spatial axes is ndim, minus one when the last dimension is time.
        select case(vs%par%ndim)

            case(1)
                ! Point (0D space), or a 1D-in-space variable with no time.
                ! No spatial axis coordinates to load here.

            case(2)
                if (with_time) then
                    ! 1D space + time
                    call load_or_generate_axis(filename_dims,dim_names(1),vs%dim(1),vs%x)
                else
                    ! 2D space
                    call load_or_generate_axis(filename_dims,dim_names(1),vs%dim(1),vs%x)
                    call load_or_generate_axis(filename_dims,dim_names(2),vs%dim(2),vs%y)
                end if

            case(3)
                if (with_time) then
                    ! 2D space + time
                    call load_or_generate_axis(filename_dims,dim_names(1),vs%dim(1),vs%x)
                    call load_or_generate_axis(filename_dims,dim_names(2),vs%dim(2),vs%y)
                else
                    ! 3D space
                    call load_or_generate_axis(filename_dims,dim_names(1),vs%dim(1),vs%x)
                    call load_or_generate_axis(filename_dims,dim_names(2),vs%dim(2),vs%y)
                    call load_or_generate_axis(filename_dims,dim_names(3),vs%dim(3),vs%z)
                end if

            case(4)
                if (with_time) then
                    ! 3D space + time
                    call load_or_generate_axis(filename_dims,dim_names(1),vs%dim(1),vs%x)
                    call load_or_generate_axis(filename_dims,dim_names(2),vs%dim(2),vs%y)
                    call load_or_generate_axis(filename_dims,dim_names(3),vs%dim(3),vs%z)
                else if (nint(vs%par%time_par(4)) .gt. 1) then
                    ! 3D space + sub-annual cycle without years (a climatology,
                    ! e.g. monthly ocean fields): the last dimension is the cycle
                    call load_or_generate_axis(filename_dims,dim_names(1),vs%dim(1),vs%x)
                    call load_or_generate_axis(filename_dims,dim_names(2),vs%dim(2),vs%y)
                    call load_or_generate_axis(filename_dims,dim_names(3),vs%dim(3),vs%z)
                else
                    call fesm_error("varslice_init_data", &
                        "4D array without time dimension needs a sub-annual period "// &
                        "(time = 0 ... nsub, nsub > 1: the last dimension is the cycle).", &
                        "name = "//trim(vs%par%name))
                end if

            case DEFAULT
                call fesm_error("varslice_init_data", &
                    "ndim > 4 not allowed.", &
                    "ndim = "//to_str(vs%par%ndim))

        end select

        ! Make sure all axis variables exist
        if (.not. allocated(vs%x)) then 
            allocate(vs%x(1))
            vs%x = 0.0_wp 
        end if 
        if (.not. allocated(vs%y)) then 
            allocate(vs%y(1))
            vs%y = 0.0_wp 
        end if 
        if (.not. allocated(vs%z)) then 
            allocate(vs%z(1))
            vs%z = 0.0_wp 
        end if 
        if (.not. allocated(vs%time)) then 
            allocate(vs%time(1))
            vs%time = 0.0_wp 
        end if 
        if (.not. allocated(vs%time_sub)) then 
            allocate(vs%time_sub(1))
            vs%time_sub = 0.0_wp 
        end if 
        if (.not. allocated(vs%idx)) then 
            allocate(vs%idx(1))
            vs%idx = 0 
        end if 

        return  

    end subroutine varslice_init_data

    subroutine alloc_var(var, src_shape, ndim, nout)
        ! (Re)allocate the 4D storage array so that the first ndim-1 axes are
        ! the spatial extents (taken from src_shape, the shape of the loaded
        ! slab), the ndim-th axis has length nout (the time-like output), and
        ! any remaining trailing axes are singleton. This keeps a single code
        ! path for every rank instead of a per-ndim allocate.

        implicit none

        real(wp), allocatable, intent(INOUT) :: var(:,:,:,:)
        integer,               intent(IN)    :: src_shape(4)
        integer,               intent(IN)    :: ndim
        integer,               intent(IN)    :: nout

        integer :: shp(4)

        shp = 1
        if (ndim .ge. 2) shp(1:ndim-1) = src_shape(1:ndim-1)
        shp(ndim) = nout

        if (allocated(var)) deallocate(var)
        allocate(var(shp(1),shp(2),shp(3),shp(4)))

        return

    end subroutine alloc_var

    function varslice_v1(vs, k) result(v)
        ! 1D (vector) copy at trailing index k: var(:,k,1,1)
        class(varslice_class), intent(IN) :: vs
        integer, optional,     intent(IN) :: k
        real(wp), allocatable :: v(:)
        integer :: kk
        kk = 1
        if (present(k)) kk = k
        v = vs%var(:,kk,1,1)
    end function varslice_v1

    function varslice_v2(vs, k) result(v)
        ! 2D (field) copy at trailing index k: var(:,:,k,1)
        class(varslice_class), intent(IN) :: vs
        integer, optional,     intent(IN) :: k
        real(wp), allocatable :: v(:,:)
        integer :: kk
        kk = 1
        if (present(k)) kk = k
        v = vs%var(:,:,kk,1)
    end function varslice_v2

    function varslice_v3(vs, k) result(v)
        ! 3D (volume) copy at trailing index k: var(:,:,:,k)
        class(varslice_class), intent(IN) :: vs
        integer, optional,     intent(IN) :: k
        real(wp), allocatable :: v(:,:,:)
        integer :: kk
        kk = 1
        if (present(k)) kk = k
        v = vs%var(:,:,:,kk)
    end function varslice_v3

    subroutine varslice_get_0d_dp(vs, var, k)
        type(varslice_class), intent(IN)  :: vs
        real(dp),             intent(OUT) :: var
        integer, optional,    intent(IN)  :: k
        integer :: kk
        kk = get_k(vs, k, 1, "0-D")
        var = real(vs%var(kk,1,1,1), dp)
    end subroutine varslice_get_0d_dp

    subroutine varslice_get_1d_dp(vs, var, k)
        type(varslice_class), intent(IN)  :: vs
        real(dp),             intent(OUT) :: var(:)
        integer, optional,    intent(IN)  :: k
        integer :: kk
        kk = get_k(vs, k, 2, "1-D")
        call get_shape(vs, [size(var,1), 1, 1], [size(vs%var,1), 1, 1], "1-D")
        var = real(vs%var(:,kk,1,1), dp)
    end subroutine varslice_get_1d_dp

    subroutine varslice_get_2d_dp(vs, var, k)
        type(varslice_class), intent(IN)  :: vs
        real(dp),             intent(OUT) :: var(:,:)
        integer, optional,    intent(IN)  :: k
        integer :: kk
        kk = get_k(vs, k, 3, "2-D")
        call get_shape(vs, [size(var,1), size(var,2), 1], [size(vs%var,1), size(vs%var,2), 1], "2-D")
        var = real(vs%var(:,:,kk,1), dp)
    end subroutine varslice_get_2d_dp

    subroutine varslice_get_3d_dp(vs, var, k)
        type(varslice_class), intent(IN)  :: vs
        real(dp),             intent(OUT) :: var(:,:,:)
        integer, optional,    intent(IN)  :: k
        integer :: kk
        kk = get_k(vs, k, 4, "3-D")
        call get_shape(vs, shape(var), [size(vs%var,1), size(vs%var,2), size(vs%var,3)], "3-D")
        var = real(vs%var(:,:,:,kk), dp)
    end subroutine varslice_get_3d_dp

    subroutine varslice_get_0d_sp(vs, var, k)
        type(varslice_class), intent(IN)  :: vs
        real(sp),             intent(OUT) :: var
        integer, optional,    intent(IN)  :: k
        integer :: kk
        kk = get_k(vs, k, 1, "0-D")
        var = real(vs%var(kk,1,1,1), sp)
    end subroutine varslice_get_0d_sp

    subroutine varslice_get_1d_sp(vs, var, k)
        type(varslice_class), intent(IN)  :: vs
        real(sp),             intent(OUT) :: var(:)
        integer, optional,    intent(IN)  :: k
        integer :: kk
        kk = get_k(vs, k, 2, "1-D")
        call get_shape(vs, [size(var,1), 1, 1], [size(vs%var,1), 1, 1], "1-D")
        var = real(vs%var(:,kk,1,1), sp)
    end subroutine varslice_get_1d_sp

    subroutine varslice_get_2d_sp(vs, var, k)
        type(varslice_class), intent(IN)  :: vs
        real(sp),             intent(OUT) :: var(:,:)
        integer, optional,    intent(IN)  :: k
        integer :: kk
        kk = get_k(vs, k, 3, "2-D")
        call get_shape(vs, [size(var,1), size(var,2), 1], [size(vs%var,1), size(vs%var,2), 1], "2-D")
        var = real(vs%var(:,:,kk,1), sp)
    end subroutine varslice_get_2d_sp

    subroutine varslice_get_3d_sp(vs, var, k)
        type(varslice_class), intent(IN)  :: vs
        real(sp),             intent(OUT) :: var(:,:,:)
        integer, optional,    intent(IN)  :: k
        integer :: kk
        kk = get_k(vs, k, 4, "3-D")
        call get_shape(vs, shape(var), [size(vs%var,1), size(vs%var,2), size(vs%var,3)], "3-D")
        var = real(vs%var(:,:,:,kk), sp)
    end subroutine varslice_get_3d_sp

    integer function get_k(vs, k, dim, what) result(kk)
        ! The trailing index of varslice_get (default 1), checked against
        ! dimension dim of vs%var.
        type(varslice_class), intent(IN) :: vs
        integer, optional,    intent(IN) :: k
        integer,              intent(IN) :: dim
        character(len=*),     intent(IN) :: what
        if (.not. allocated(vs%var)) then
            call fesm_error("varslice_get", "no slice yet (call varslice_update first).", &
                "name = "//trim(vs%par%name))
        end if
        kk = 1
        if (present(k)) kk = k
        if (kk .lt. 1 .or. kk .gt. size(vs%var,dim)) then
            call fesm_error("varslice_get", "index out of range.", &
                "name  = "//trim(vs%par%name)//new_line("a")// &
                "rank  = "//what//new_line("a")// &
                "k     = "//to_str(kk)//new_line("a")// &
                "shape = "//to_str(shape(vs%var)))
        end if
    end function get_k

    subroutine get_shape(vs, shp, shp_vs, what)
        type(varslice_class), intent(IN) :: vs
        integer,              intent(IN) :: shp(3), shp_vs(3)
        character(len=*),     intent(IN) :: what
        if (any(shp .ne. shp_vs)) then
            call fesm_error("varslice_get", "shape mismatch.", &
                "name   = "//trim(vs%par%name)//new_line("a")// &
                "rank   = "//what//new_line("a")// &
                "slice  = "//to_str(shp_vs)//new_line("a")// &
                "given  = "//to_str(shp))
        end if
    end subroutine get_shape

    subroutine load_or_generate_axis(filename, varname, n, x)
        ! Allocate axis x to length n and populate it: read the coordinate
        ! variable from file if it exists, otherwise generate a 1..n index axis.

        implicit none

        character(len=*),      intent(IN)    :: filename
        character(len=*),      intent(IN)    :: varname
        integer,               intent(IN)    :: n
        real(wp), allocatable, intent(INOUT) :: x(:)

        if (allocated(x)) deallocate(x)
        allocate(x(n))

        if (nc_exists_var(filename, varname)) then
            call nc_read(filename, varname, x)
        else
            call axis_init(x, nx=n)
        end if

        return

    end subroutine load_or_generate_axis

    subroutine file_time_axis(par, t, tb)
        ! A time axis from the file. The mode (auto | static | clim |
        ! transient) is resolved from the variable's last dimension: a time
        ! axis if its coordinate variable exists and is named time, has
        ! axis = "T", standard_name = "time", a decodable units attribute,
        ! the CF climatology attribute (then: clim), or time_units is set.
        ! transient: the axis t of every file, decoded to calendar years
        ! (varslice_time_decode) minus time_ref, concatenated, strictly
        ! increasing; nsub = steps per year when every year holds the same
        ! number of them (> 1), else 1. clim: the last dimension is the
        ! cycle, nsub = its length. tb: the CF bounds of the axis (the time
        ! variable's bounds attribute), decoded alike, when every file has
        ! them.

        implicit none

        type(varslice_param_class), intent(INOUT) :: par
        real(dp), allocatable,      intent(OUT)   :: t(:)
        real(dp), allocatable,      intent(OUT)   :: tb(:,:)

        character(len=56), allocatable :: dn(:)
        integer,  allocatable :: dl(:), cnt(:)
        real(dp), allocatable :: tf(:), yf(:), y(:), bf(:,:), yb(:), tbf(:,:)
        character(len=256) :: units, cal, att, bname
        logical  :: with_bnds
        character(len=:), allocatable :: fnames
        logical  :: has_coord, is_time, is_clim
        integer  :: nd, i, n0, nt, k, stat, nsub, ny0
        real(dp) :: dx

        call nc_dims(par%filenames(1), par%name, dn, dl)
        nd = size(dl)
        par%time_name = dn(nd)

        has_coord = nc_exists_var(par%filenames(1), trim(par%time_name))
        is_time   = .false.
        is_clim   = .false.
        if (has_coord) then
            is_clim = nc_exists_attr(par%filenames(1), trim(par%time_name), "climatology")
            is_time = is_clim .or. trim(par%time_name) .eq. "time" .or. len_trim(par%time_units) .gt. 0
            if (.not. is_time) is_time = attr_is(par%filenames(1), par%time_name, "axis", "t")
            if (.not. is_time) is_time = attr_is(par%filenames(1), par%time_name, "standard_name", "time")
            if (.not. is_time .and. nc_exists_attr(par%filenames(1), trim(par%time_name), "units")) then
                call nc_read_attr(par%filenames(1), trim(par%time_name), "units", att)
                call varslice_time_decode([0.0_dp], att, "", y, stat=stat, quiet=.true.)
                is_time = (stat .eq. 0)
            end if
        end if

        select case(trim(par%mode))
            case("auto")
                par%mode = "static"
                if (is_time) par%mode = "transient"
                if (is_clim) par%mode = "clim"
            case("transient")
                if (.not. has_coord) then
                    call fesm_error("varslice_init_data", &
                        "mode = ""transient"", but the last dimension of the variable has no coordinate variable.", &
                        "name      = "//trim(par%name)//new_line("a")// &
                        "dimension = "//trim(par%time_name)//new_line("a")// &
                        "file      = "//trim(par%filenames(1)))
                end if
        end select

        par%with_time     = .false.
        par%with_time_sub = .false.
        par%time_par      = 0.0_dp

        select case(trim(par%mode))

            case("static")
                return

            case("clim")
                nsub = dl(nd)
                if (nsub .lt. 2) then
                    call fesm_error("varslice_init_data", &
                        "mode = ""clim"" needs a cycle of at least two steps (the last dimension).", &
                        "name      = "//trim(par%name)//new_line("a")// &
                        "dimension = "//trim(par%time_name)//" ("//to_str(nsub)//")"//new_line("a")// &
                        "file      = "//trim(par%filenames(1)))
                end if
                par%time_par(4)   = real(nsub, dp)
                par%with_time_sub = .true.
                return

        end select

        ! transient: decode and concatenate the axis of every file
        allocate(t(0), tbf(2,0))
        with_bnds = .true.
        fnames = ""
        do i = 1, size(par%filenames)
            nt = nc_size(par%filenames(i), trim(par%time_name))
            if (allocated(tf)) deallocate(tf)
            allocate(tf(nt), yf(nt))
            call nc_read(par%filenames(i), trim(par%time_name), tf)

            units = par%time_units
            if (len_trim(units) .eq. 0) then
                if (.not. nc_exists_attr(par%filenames(i), trim(par%time_name), "units")) then
                    call fesm_error("varslice_init_data", &
                        "the time coordinate has no units attribute: set time_units.", &
                        "name = "//trim(par%name)//new_line("a")// &
                        "time = "//trim(par%time_name)//new_line("a")// &
                        "file = "//trim(par%filenames(i)))
                end if
                call nc_read_attr(par%filenames(i), trim(par%time_name), "units", units)
            end if
            cal = par%calendar
            if (len_trim(cal) .eq. 0) then
                if (nc_exists_attr(par%filenames(i), trim(par%time_name), "calendar")) &
                    call nc_read_attr(par%filenames(i), trim(par%time_name), "calendar", cal)
            end if

            if (trim(units) .eq. "clock") then
                yf = tf
            else
                call varslice_time_decode(tf, units, cal, yf, stat=stat)
                if (stat .ne. 0) then
                    call fesm_error("varslice_init_data", &
                        "the time axis of the file cannot be decoded (see above).", &
                        "name     = "//trim(par%name)//new_line("a")// &
                        "time     = "//trim(par%time_name)//new_line("a")// &
                        "units    = "//trim(units)//new_line("a")// &
                        "calendar = "//trim(cal)//new_line("a")// &
                        "file     = "//trim(par%filenames(i)))
                end if
                yf = yf - par%time_ref
            end if

            ! CF bounds, decoded alike
            bname = ""
            if (nc_exists_attr(par%filenames(i), trim(par%time_name), "bounds")) &
                call nc_read_attr(par%filenames(i), trim(par%time_name), "bounds", bname)
            if (with_bnds .and. len_trim(bname) .gt. 0) then
                with_bnds = nc_exists_var(par%filenames(i), trim(bname))
            else
                with_bnds = .false.
            end if
            if (with_bnds) then
                allocate(bf(2,nt))
                call nc_read(par%filenames(i), trim(bname), bf)
                if (trim(units) .eq. "clock") then
                    tbf = reshape([reshape(tbf,[size(tbf)]), reshape(bf,[2*nt])], [2, size(tbf,2)+nt])
                else
                    call varslice_time_decode(reshape(bf,[2*nt]), units, cal, yb)
                    yb = yb - par%time_ref
                    tbf = reshape([reshape(tbf,[size(tbf)]), yb], [2, size(tbf,2)+nt])
                end if
                deallocate(bf)
            end if

            t = [t, yf]
            deallocate(yf)
            fnames = fnames//new_line("a")//"  "//trim(par%filenames(i))//" ["//to_str(nt)//"]"
        end do

        nt = size(t)
        do k = 2, nt
            if (.not. (t(k) .gt. t(k-1))) then
                call fesm_error("varslice_init_data", &
                    "the time axis is not strictly increasing.", &
                    "name  = "//trim(par%name)//new_line("a")// &
                    "index = "//to_str(k)//": "//to_str(t(k-1))//" -> "//to_str(t(k))//new_line("a")// &
                    "files:"//fnames)
            end if
        end do

        ! Steps per year: the same number (> 1) in every year of the axis
        nsub = 1
        if (nt .gt. 1) then
            ny0 = floor(t(1))
            allocate(cnt(floor(t(nt)) - ny0 + 1))
            cnt = 0
            do k = 1, nt
                cnt(floor(t(k)) - ny0 + 1) = cnt(floor(t(k)) - ny0 + 1) + 1
            end do
            n0 = cnt(1)
            if (n0 .gt. 1 .and. all(cnt .eq. n0)) nsub = n0
        end if
        if (nsub .gt. 1 .and. abs(par%time_ref - anint(par%time_ref)) .gt. 0.0_dp) then
            call fesm_error("varslice_init_data", &
                "a sub-annual axis needs time_ref in whole years.", &
                "name     = "//trim(par%name)//new_line("a")// &
                "time_ref = "//to_str(par%time_ref))
        end if

        dx = 0.0_dp
        if (nt .gt. 1) then
            dx = t(2) - t(1)
            if (any(abs((t(2:nt) - t(1:nt-1)) - dx) .gt. 1e-6_dp*max(1.0_dp, abs(dx)))) dx = 0.0_dp
        end if

        par%with_time     = .true.
        par%with_time_sub = (nsub .gt. 1)
        par%time_par      = [t(1), t(nt), dx, real(nsub, dp)]

        if (with_bnds) then
            if (any(tbf(2,:) .lt. tbf(1,:))) then
                call fesm_error("varslice_init_data", &
                    "time bounds with an upper bound below the lower one.", &
                    "name = "//trim(par%name)//new_line("a")//"files:"//fnames)
            end if
            tb = tbf
        end if

    contains

        logical function attr_is(fname, vname, att, val)
            character(len=*), intent(IN) :: fname, vname, att, val
            character(len=256) :: a
            attr_is = .false.
            if (.not. nc_exists_attr(fname, trim(vname), att)) return
            call nc_read_attr(fname, trim(vname), att, a)
            attr_is = (trim(lower(a)) .eq. val)
        end function attr_is

    end subroutine file_time_axis

    subroutine varslice_time_decode(t, units, calendar, y, stat, quiet)
        ! Calendar years y of the values t of a time coordinate with CF units
        ! and calendar (case-insensitive):
        !   {seconds,minutes,hours,days} since <date>: by calendar
        !       standard | gregorian (exact day count; dates before
        !       1582-10-15: error), proleptic_gregorian, noleap | 365_day,
        !       all_leap | 366_day, 360_day ("" = standard);
        !       y = year + day-of-year / days-in-year
        !   months since <date>: 360_day only
        !   years since <date>, common_years since <date>: year of date + t
        !   yr BP, a BP, years BP, ka BP, kyr BP: 1950 - t*{1, 1000}
        !   yr CE, years CE, AD: t
        !   yr, years, ka (no reference): error (set time_units)
        ! Errors stop, or with stat present set stat = 1 and return (quiet:
        ! without a message).

        implicit none

        real(dp),         intent(IN)  :: t(:)
        character(len=*), intent(IN)  :: units, calendar
        real(dp), allocatable, intent(OUT) :: y(:)
        integer,          intent(OUT), optional :: stat
        logical,          intent(IN),  optional :: quiet

        character(len=256) :: u, word, date
        character(len=32)  :: cal
        integer  :: p, yr0, mo0, dy0, k, yy
        real(dp) :: sec0, fac, ref, a, ystart, ylen
        logical  :: ok, loud

        loud = .true.
        if (present(quiet)) loud = .not. quiet
        if (present(stat)) stat = 0
        allocate(y(size(t)))
        y = 0.0_dp

        u   = adjustl(lower(units))
        cal = adjustl(lower(calendar))
        if (len_trim(cal) .eq. 0) cal = "standard"

        p = index(u, " since ")
        if (p .eq. 0) then
            select case(trim(squeeze(u)))
                case("yr bp", "a bp", "years bp")
                    y = 1950.0_dp - t
                case("ka bp", "kyr bp")
                    y = 1950.0_dp - 1000.0_dp*t
                case("yr ce", "years ce", "ad")
                    y = t
                case("yr", "years", "ka", "a", "year", "kyr")
                    call fail("time units without a reference: set time_units", "")
                case default
                    call fail("time units not recognised", "")
            end select
            return
        end if

        word = u(1:p-1)
        date = adjustl(u(p+7:))
        call parse_date(date, yr0, mo0, dy0, sec0, ok)
        if (.not. ok) then
            call fail("reference date not recognised", "date = "//trim(date))
            return
        end if

        select case(trim(word))
            case("years", "year", "common_years", "common_year")
                y = real(yr0, dp) + t
                return
            case("seconds", "second", "secs", "sec", "s")
                fac = 1.0_dp/86400.0_dp
            case("minutes", "minute", "mins", "min")
                fac = 1.0_dp/1440.0_dp
            case("hours", "hour", "hrs", "hr", "h")
                fac = 1.0_dp/24.0_dp
            case("days", "day", "d")
                fac = 1.0_dp
            case("months", "month")
                if (trim(cal) .ne. "360_day") then
                    call fail("months since <date> needs the 360_day calendar", "calendar = "//trim(cal))
                    return
                end if
                fac = 30.0_dp
            case default
                call fail("time units not recognised", "")
                return
        end select

        select case(trim(cal))
            case("standard", "gregorian", "proleptic_gregorian", "noleap", "365_day", &
                 "all_leap", "366_day", "360_day")
            case default
                call fail("calendar not supported", "calendar = "//trim(cal))
                return
        end select

        if (mo0 .lt. 1 .or. mo0 .gt. 12 .or. dy0 .lt. 1 .or. dy0 .gt. month_len(yr0, mo0)) then
            call fail("reference date not valid in the calendar", "date = "//trim(date))
            return
        end if

        ref = day_number(yr0, mo0, dy0) + sec0/86400.0_dp

        do k = 1, size(t)
            a = ref + t(k)*fac
            if ((trim(cal) .eq. "standard" .or. trim(cal) .eq. "gregorian") .and. &
                (a .lt. day_number(1582, 10, 15) .or. ref .lt. day_number(1582, 10, 15))) then
                call fail("dates before 1582-10-15 in the standard calendar (use proleptic_gregorian)", "")
                return
            end if
            ! The year holding day number a
            yy = floor(a/year_mean()) + epoch_year()
            do while (day_number(yy, 1, 1) .gt. a)
                yy = yy - 1
            end do
            do while (day_number(yy+1, 1, 1) .le. a)
                yy = yy + 1
            end do
            ystart = day_number(yy, 1, 1)
            ylen   = day_number(yy+1, 1, 1) - ystart
            y(k)   = real(yy, dp) + (a - ystart)/ylen
        end do

    contains

        subroutine fail(msg, detail)
            character(len=*), intent(IN) :: msg, detail
            character(len=:), allocatable :: d
            d = "units    = "//trim(units)//new_line("a")//"calendar = "//trim(calendar)
            if (len_trim(detail) .gt. 0) d = d//new_line("a")//trim(detail)
            if (present(stat)) then
                stat = 1
                if (.not. loud) return
            end if
            call fesm_error("varslice_time_decode", msg, d, stat)
        end subroutine fail

        logical function is_leap(yy)
            integer, intent(IN) :: yy
            select case(trim(cal))
                case("noleap", "365_day", "360_day")
                    is_leap = .false.
                case("all_leap", "366_day")
                    is_leap = .true.
                case default
                    is_leap = (mod(yy,4) .eq. 0 .and. mod(yy,100) .ne. 0) .or. mod(yy,400) .eq. 0
            end select
        end function is_leap

        integer function month_len(yy, mm)
            integer, intent(IN) :: yy, mm
            integer, parameter :: ML(12) = [31,28,31,30,31,30,31,31,30,31,30,31]
            if (trim(cal) .eq. "360_day") then
                month_len = 30
            else
                month_len = ML(mm)
                if (mm .eq. 2 .and. is_leap(yy)) month_len = 29
            end if
        end function month_len

        real(dp) function day_number(yy, mm, dd)
            ! Days since 0000-01-01 of the calendar (proleptic for the
            ! Gregorian ones).
            integer, intent(IN) :: yy, mm, dd
            integer :: i, ym
            real(dp) :: d
            select case(trim(cal))
                case("360_day")
                    d = 360.0_dp*yy + 30.0_dp*(mm-1) + (dd-1)
                case("noleap", "365_day")
                    d = 365.0_dp*yy
                    do i = 1, mm-1
                        d = d + month_len(yy, i)
                    end do
                    d = d + (dd-1)
                case("all_leap", "366_day")
                    d = 366.0_dp*yy
                    do i = 1, mm-1
                        d = d + month_len(yy, i)
                    end do
                    d = d + (dd-1)
                case default
                    ! Gregorian: leap days of the years before yy
                    ym = yy - 1
                    d  = 365.0_dp*yy + real(floor_div(ym,4) - floor_div(ym,100) + floor_div(ym,400) + 1, dp)
                    do i = 1, mm-1
                        d = d + month_len(yy, i)
                    end do
                    d = d + (dd-1)
            end select
            day_number = d
        end function day_number

        integer function floor_div(a, b)
            integer, intent(IN) :: a, b
            floor_div = floor(real(a, dp)/real(b, dp))
        end function floor_div

        real(dp) function year_mean()
            select case(trim(cal))
                case("360_day")
                    year_mean = 360.0_dp
                case("noleap", "365_day")
                    year_mean = 365.0_dp
                case("all_leap", "366_day")
                    year_mean = 366.0_dp
                case default
                    year_mean = 365.2425_dp
            end select
        end function year_mean

        integer function epoch_year()
            epoch_year = 0
        end function epoch_year

    end subroutine varslice_time_decode

    subroutine parse_date(str, yr, mo, dy, sec, ok)
        ! "Y-M-D[ h:m:s]" (also "T" between date and time; a leading "-"
        ! for negative years; a trailing "z"/"utc" ignored). Missing month
        ! or day: 1.

        implicit none

        character(len=*), intent(IN)  :: str
        integer,          intent(OUT) :: yr, mo, dy
        real(dp),         intent(OUT) :: sec
        logical,          intent(OUT) :: ok

        character(len=256) :: s, d, tm
        integer :: p, q, ios, sgn, hh, mi
        real(dp) :: ss

        yr = 0; mo = 1; dy = 1; sec = 0.0_dp
        ok = .false.

        s = adjustl(str)
        p = index(trim(s), " ")
        q = index(s, "t")
        if (q .gt. 1 .and. (p .eq. 0 .or. q .lt. p)) p = q
        if (p .gt. 0) then
            d  = s(1:p-1)
            tm = adjustl(s(p+1:))
        else
            d  = s
            tm = ""
        end if
        p = index(tm, "z");   if (p .gt. 0) tm = tm(1:p-1)
        p = index(tm, "utc"); if (p .gt. 0) tm = tm(1:p-1)

        sgn = 1
        if (d(1:1) .eq. "-") then
            sgn = -1
            d = d(2:)
        end if

        p = index(d, "-")
        if (p .eq. 0) then
            read(d, *, iostat=ios) yr
            if (ios .ne. 0) return
        else
            read(d(1:p-1), *, iostat=ios) yr
            if (ios .ne. 0) return
            d = d(p+1:)
            p = index(d, "-")
            if (p .eq. 0) then
                read(d, *, iostat=ios) mo
                if (ios .ne. 0) return
            else
                read(d(1:p-1), *, iostat=ios) mo
                if (ios .ne. 0) return
                read(d(p+1:), *, iostat=ios) dy
                if (ios .ne. 0) return
            end if
        end if
        yr = sgn*yr

        if (len_trim(tm) .gt. 0) then
            hh = 0; mi = 0; ss = 0.0_dp
            p = index(tm, ":")
            if (p .eq. 0) then
                read(tm, *, iostat=ios) hh
                if (ios .ne. 0) return
            else
                read(tm(1:p-1), *, iostat=ios) hh
                if (ios .ne. 0) return
                tm = tm(p+1:)
                p = index(tm, ":")
                if (p .eq. 0) then
                    read(tm, *, iostat=ios) mi
                    if (ios .ne. 0) return
                else
                    read(tm(1:p-1), *, iostat=ios) mi
                    if (ios .ne. 0) return
                    read(tm(p+1:), *, iostat=ios) ss
                    if (ios .ne. 0) return
                end if
            end if
            sec = 3600.0_dp*hh + 60.0_dp*mi + ss
        end if

        ok = .true.

    end subroutine parse_date

    pure function lower(s) result(l)
        character(len=*), intent(IN) :: s
        character(len=len(s)) :: l
        integer :: i, c
        l = s
        do i = 1, len(s)
            c = iachar(s(i:i))
            if (c .ge. iachar("A") .and. c .le. iachar("Z")) l(i:i) = achar(c + 32)
        end do
    end function lower

    pure function squeeze(s) result(r)
        ! Runs of blanks reduced to one, leading blanks removed
        character(len=*), intent(IN) :: s
        character(len=len(s)) :: r
        integer :: i, n
        logical :: blank
        r = ""
        n = 0
        blank = .true.
        do i = 1, len_trim(s)
            if (s(i:i) .eq. " ") then
                if (blank) cycle
                blank = .true.
            else
                blank = .false.
            end if
            n = n + 1
            r(n:n) = s(i:i)
        end do
    end function squeeze

    subroutine varslice_end(vs)
        ! Deallocate all variables

        implicit none 

        type(varslice_class), intent(INOUT) :: vs 

        if (allocated(vs%dim))          deallocate(vs%dim)
        if (allocated(vs%x))            deallocate(vs%x)
        if (allocated(vs%y))            deallocate(vs%y)
        if (allocated(vs%z))            deallocate(vs%z)
        if (allocated(vs%time))         deallocate(vs%time)
        if (allocated(vs%time_sub))     deallocate(vs%time_sub)
        if (allocated(vs%time_bnds))    deallocate(vs%time_bnds)
        if (allocated(vs%idx))          deallocate(vs%idx)
        if (allocated(vs%nt_files))     deallocate(vs%nt_files)
        if (allocated(vs%var))          deallocate(vs%var)

        return

    end subroutine varslice_end

    function varslice_nsub(vs) result(nsub)
        ! Number of sub-annual steps per year declared for this variable in its
        ! namelist `time` entry (the 5th value, stored as time_par(4)). Returns
        ! 12 for monthly data and 1 for annual/static data. This is the natural
        ! value of `rep` for slice methods that preserve the sub-annual cycle,
        ! letting callers source it from the data rather than hard-coding it.

        implicit none

        type(varslice_class), intent(IN) :: vs
        integer :: nsub

        nsub = max(1, nint(vs%par%time_par(4)))

        return

    end function varslice_nsub

    function varslice_sub_mean(vs) result(v)
        ! Mean of the current slice over its time-like axis (e.g. the sub-annual
        ! cycle after an update with rep=varslice_nsub(vs)), as an (x,y,z) field.
        ! The time-like axis is axis ndim of the natural-rank storage (var(:,:,:,k)
        ! for 3D space, var(:,:,k,1) for 2D space, ...), so the caller does not
        ! need to know the variable's rank: spatial axes the variable lacks are
        ! singleton, e.g. a 2D field gives (nx,ny,1), the same layout as
        ! var(:,:,:,1) of a rep=1 slice. Missing values are excluded from the
        ! mean (mv where none is available). A variable without time is returned
        ! unchanged.

        implicit none

        type(varslice_class), intent(IN) :: vs
        real(wp), allocatable :: v(:,:,:)

        ! Local variables
        integer :: nd, nt, nsp, s
        integer :: shp(3)
        real(wp), allocatable :: src(:,:), dst(:)

        if (.not. vs%par%with_time .and. varslice_nsub(vs) .le. 1) then
            v = vs%var(:,:,:,1)
            return
        end if

        nd  = vs%par%ndim
        nt  = size(vs%var,nd)
        nsp = size(vs%var) / nt

        ! Spatial extents (axes 1..nd-1), singleton for absent axes
        shp = 1
        do s = 1, nd-1
            shp(s) = size(vs%var,s)
        end do

        ! Axes after nd are singleton, so the storage is (space, time) in memory order
        src = reshape(vs%var, [nsp,nt])

        allocate(dst(nsp))
        do s = 1, nsp
            call calc_vec_value(dst(s),src(s,:),"mean",mv)
        end do

        v = reshape(dst, shp)

        return

    end function varslice_sub_mean

    subroutine varslice_par_load(par,filename,group,domain,grid_name,verbose,subs, &
                                 defaults_file,defaults_group,mode,remap,grid_src,time_ref)

        type(varslice_param_class), intent(OUT) :: par
        character(len=*), intent(IN) :: filename
        character(len=*), intent(IN) :: group
        character(len=*), intent(IN), optional :: domain
        character(len=*), intent(IN), optional :: grid_name
        logical, optional :: verbose
        character(len=*), intent(IN), optional :: subs(:,:)   ! extra {key}->value path substitutions
        character(len=*), intent(IN), optional :: defaults_file, defaults_group
        character(len=*), intent(IN), optional :: mode, remap, grid_src
        real(dp),         intent(IN), optional :: time_ref

        ! Local variables
        logical  :: print_summary
        integer  :: i
        character(len=56) :: units(2)     ! [units_in, units_out]
        real(wp)          :: scaling(2)   ! [unit_scale, unit_offset]
        real(dp)          :: time(5)      ! [active, x0, x1, dx, sub]
        character(len=16) :: mode_key

        print_summary = .TRUE.
        if (present(verbose)) print_summary = verbose

        mode_key = "auto"

        if (present(defaults_file)) then
            ! Every key through the schema
            call par_load_schema()
        else
            ! Legacy path
            call par_load_legacy()
        end if

        ! Argument overrides
        if (present(mode))     mode_key     = mode
        if (present(remap))    par%remap    = remap
        if (present(grid_src)) par%grid_src = grid_src
        par%time_ref = 0.0_dp
        if (present(time_ref)) par%time_ref = time_ref

        ! Unpack into the internal parameter fields
        par%units_in    = units(1)
        par%units_out   = units(2)
        par%unit_scale  = scaling(1)
        par%unit_offset = scaling(2)
        par%with_time   = (time(1) .ge. 0.5_dp)
        par%time_par    = time(2:5)

        call resolve_mode()

        ! Parse filename placeholders ({domain}/{grid_name} and any extra subs).
        ! parse_path skips whichever optionals are absent.
        call parse_path(par%filename,domain,grid_name,subs)
        if (len_trim(par%grid_src) .gt. 0) call parse_path(par%grid_src,domain,grid_name,subs)

        ! Resolve file list and derive internal time parameters
        call varslice_par_finalize(par)

        ! Summary 
        if (print_summary) then  
            write(*,*) "Loading: ", trim(filename), ":: ", trim(group)
            write(*,*) "filename      = ", trim(par%filename)
            if (size(par%filenames,1) .gt. 1) then
                write(*,*) "filenames     = "
                do i = 1, size(par%filenames,1)
                    write(*,*) "                  ", trim(par%filenames(i))
                end do
            end if
            write(*,*) "name          = ", trim(par%name)
            write(*,*) "units_in      = ", trim(par%units_in)
            write(*,*) "units_out     = ", trim(par%units_out)
            write(*,*) "unit_scale    = ", par%unit_scale
            write(*,*) "unit_offset   = ", par%unit_offset
            write(*,*) "with_time     = ", par%with_time
            if (trim(par%remap) .ne. "none") then
                write(*,*) "remap         = ", trim(par%remap)
                write(*,*) "grid_src      = ", trim(par%grid_src)
            end if
            if (trim(par%fill) .ne. "none") then
                write(*,*) "fill          = ", trim(par%fill)
                if (trim(par%fill) .eq. "const") write(*,*) "fill_value    = ", par%fill_value
            end if
            write(*,*) "with_time_sub = ", par%with_time_sub
            if (par%with_time) then
                write(*,*) "time_par    = ", par%time_par
            end if
        end if 

        return

    contains

        subroutine par_load_schema()
            ! Every key from the group, else from the schema group.

            character(len=256) :: dg

            dg = group
            if (present(defaults_group)) dg = defaults_group

            par%filename = ""
            par%name     = ""
            units        = ""
            scaling      = [1.0_wp, 0.0_wp]
            time         = 0.0_dp
            par%remap    = "none"
            par%grid_src = ""
            par%fill     = "none"
            par%fill_value = mv

            call nml_read(filename,group,"filename",  par%filename,   defaults_file=defaults_file, defaults_group=dg)
            call nml_read(filename,group,"name",      par%name,       defaults_file=defaults_file, defaults_group=dg)
            call nml_read(filename,group,"units",     units,          defaults_file=defaults_file, defaults_group=dg)
            call nml_read(filename,group,"scaling",   scaling,        defaults_file=defaults_file, defaults_group=dg)
            call nml_read(filename,group,"time",      time,           defaults_file=defaults_file, defaults_group=dg)
            call nml_read(filename,group,"mode",      mode_key,       defaults_file=defaults_file, defaults_group=dg)
            call nml_read(filename,group,"remap",     par%remap,      defaults_file=defaults_file, defaults_group=dg)
            call nml_read(filename,group,"grid_src",  par%grid_src,   defaults_file=defaults_file, defaults_group=dg)
            call nml_read(filename,group,"fill",      par%fill,       defaults_file=defaults_file, defaults_group=dg)
            call nml_read(filename,group,"fill_value",par%fill_value, defaults_file=defaults_file, defaults_group=dg)
            call nml_read(filename,group,"time_units",par%time_units, defaults_file=defaults_file, defaults_group=dg)
            call nml_read(filename,group,"calendar",  par%calendar,   defaults_file=defaults_file, defaults_group=dg)
            call check_fill_method(par%fill, "varslice_par_load")

            par%time_from_file = (time(1) .lt. 0.0_dp)

        end subroutine par_load_schema

        subroutine par_load_legacy()

        ! Read the condensed namelist keys. Each value bundles related
        ! quantities so the group stays short but fully explicit:
        !   units   = "<in>" "<out>"
        !   scaling = <unit_scale> <unit_offset>
        !   time    = <active> <x0> <x1> <dx> <sub>
        ! `active` (1.0/0.0) flags whether the variable is time-varying, or
        ! whether x0..x1 merely document the period the (static) data represents.
        ! `sub` is optional (trailing) and defaults to 0; every other value is
        ! required, keeping the provenance visible in the file.

            ! Pre-set defaults so an omitted trailing element falls back sensibly
            ! (the keys themselves remain mandatory).
            units   = ""
            scaling = [1.0_wp, 0.0_wp]
            time    = 0.0_dp

            call nml_read(filename,group,"filename", par%filename)
            call nml_read(filename,group,"name",     par%name)
            call nml_read(filename,group,"units",    units)
            call nml_read(filename,group,"scaling",  scaling)

            ! time: the legacy axis; absent or time(1) < 0: the axis of the file
            if (nml_has_param(filename,group,"time")) then
                call nml_read(filename,group,"time", time)
                par%time_from_file = (time(1) .lt. 0.0_dp)
            else
                time = [-1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
                par%time_from_file = .true.
            end if
            if (nml_has_param(filename,group,"time_units")) &
                call nml_read(filename,group,"time_units", par%time_units)
            if (nml_has_param(filename,group,"calendar")) &
                call nml_read(filename,group,"calendar", par%calendar)

            ! Optional: online remapping onto the target grid given at init
            par%remap    = "none"
            par%grid_src = ""
            if (nml_has_param(filename,group,"remap")) then
                call nml_read(filename,group,"remap",    par%remap)
                call nml_read(filename,group,"grid_src", par%grid_src)
            end if

            ! Optional: filling of missing values (fill_value with fill = "const")
            par%fill       = "none"
            par%fill_value = mv
            if (nml_has_param(filename,group,"fill")) then
                call nml_read(filename,group,"fill", par%fill)
                call check_fill_method(par%fill, "varslice_par_load")
                if (trim(par%fill) .eq. "const") then
                    if (.not. nml_has_param(filename,group,"fill_value")) then
                        call fesm_error("varslice_par_load", &
                            "fill = ""const"" needs fill_value.", "group = "//trim(group))
                    end if
                    call nml_read(filename,group,"fill_value", par%fill_value)
                end if
            end if

            ! Optional: what the variable is (must agree with time)
            if (nml_has_param(filename,group,"mode")) call nml_read(filename,group,"mode",mode_key)

        end subroutine par_load_legacy

        subroutine resolve_mode()
            ! The mode implied by the legacy time key; an explicit mode must
            ! agree with it. A time axis from the file: resolved at init
            ! (varslice_init_data), the key is kept.

            character(len=16) :: implied

            if (par%time_from_file) then
                select case(trim(mode_key))
                    case("auto","static","clim","transient")
                        par%mode = mode_key
                    case default
                        call fesm_error("varslice_par_load", &
                            "mode must be auto, static, clim or transient.", &
                            "group = "//trim(group)//new_line("a")// &
                            "name  = "//trim(par%name)//new_line("a")// &
                            "mode  = "//trim(mode_key))
                end select
                par%with_time = .false.
                par%time_par  = 0.0_dp
                return
            end if

            if (par%with_time) then
                implied = "transient"
            else if (nint(par%time_par(4)) .gt. 1) then
                implied = "clim"
            else
                implied = "static"
            end if

            select case(trim(mode_key))
                case("auto")
                    par%mode = implied
                case("static","clim","transient")
                    if (trim(mode_key) .ne. trim(implied)) then
                        call fesm_error("varslice_par_load", &
                            "mode conflicts with the time key.", &
                            "group = "//trim(group)//new_line("a")// &
                            "name  = "//trim(par%name)//new_line("a")// &
                            "mode  = "//trim(mode_key)//new_line("a")// &
                            "time  = "//to_str(time)//" (implies "//trim(implied)//")")
                    end if
                    par%mode = implied
                case default
                    call fesm_error("varslice_par_load", &
                        "mode must be auto, static, clim or transient.", &
                        "group = "//trim(group)//new_line("a")// &
                        "name  = "//trim(par%name)//new_line("a")// &
                        "mode  = "//trim(mode_key))
            end select

        end subroutine resolve_mode

    end subroutine varslice_par_load

    subroutine varslice_par_finalize(par)
        ! Resolve the (possibly wildcard) filename into the sorted list of
        ! matching files and derive the internal time parameters. Shared by
        ! both the nml- and argument-based initialization paths so they stay
        ! consistent.

        implicit none

        type(varslice_param_class), intent(INOUT) :: par

        ! See if multiple files are available (also allocates par%filenames)
        call get_matching_files(par%filenames, par%filename)

        ! For a time-varying field with no step given, collapse to a single
        ! time slice (x1=x0). Skip this for a static field, where x0..x1 are
        ! only documenting the represented period and must be preserved.
        if (par%with_time .and. par%time_par(3) .eq. 0.0) par%time_par(2) = par%time_par(1)

        if (par%time_par(4) .gt. 1.0) then
            par%with_time_sub = .TRUE.
        else
            par%with_time_sub = .FALSE.
        end if

        return

    end subroutine varslice_par_finalize

    subroutine parse_path(path,domain,grid_name,subs)
        ! Substitute placeholder tokens in a path. {domain} and {grid_name} are the
        ! built-in tokens; `subs` supplies any additional {key}->value pairs, with
        ! subs(k,1)=key (no braces) and subs(k,2)=value, e.g. {snapshot}->"lgm".

        implicit none

        character(len=*), intent(INOUT) :: path
        character(len=*), intent(IN), optional :: domain, grid_name
        character(len=*), intent(IN), optional :: subs(:,:)

        integer :: k

        if (present(domain))    call nml_replace(path,"{domain}",   trim(domain))
        if (present(grid_name)) call nml_replace(path,"{grid_name}",trim(grid_name))

        if (present(subs)) then
            do k = 1, size(subs,1)
                call nml_replace(path,"{"//trim(subs(k,1))//"}",trim(subs(k,2)))
            end do
        end if

        return

    end subroutine parse_path
    
    subroutine axis_init_sp(x,x0,x1,dx,nx)

        implicit none 

        real(sp), allocatable, intent(OUT) :: x(:)
        real(sp), optional :: x0
        real(sp), optional :: x1
        real(sp), optional :: dx
        integer,  optional :: nx 

        ! Local variables 
        integer :: i  
        real(sp) :: x0_now
        real(sp) :: x1_now
        real(sp) :: dx_now
        integer  :: nx_now 
        real(sp) :: nx_check

        dx_now = 1.0_sp 
        if (present(dx)) dx_now = dx 

        x0_now = 0.0_sp 
        if (present(x0)) x0_now = x0 
        
        ! Note: if x1 is present, nx is ignored 
        if (present(x1)) then 
            x1_now = x1 
        else if (present(nx)) then 
            x1_now = x0_now + (nx-1)*dx_now 
        else
            call fesm_error("axis_init", &
                "either x1 or nx must be present.")
        end if

        if (allocated(x)) deallocate(x)

        nx_now   = nint((x1_now-x0_now)/dx_now) + 1
        nx_check = (x1_now-x0_now)/dx_now + 1

        if ( abs(nx_now - nx_check) .gt. TOL ) then
            ! Make sure nx is a round number.
            call fesm_error("axis_init", &
                "desired axis bounds [x0,x1] do not divide evenly with dx.", &
                "x0 = "//to_str(x0_now)//new_line("a")// &
                "x1 = "//to_str(x1_now)//new_line("a")// &
                "dx = "//to_str(dx_now)//new_line("a")// &
                "nx = "//to_str((x1_now-x0_now)/dx_now + 1))
        end if
        
        allocate(x(nx_now))

        do i = 1, nx_now 
            x(i) = x0_now + dx_now*real(i-1,wp)
        end do 

        return

    end subroutine axis_init_sp
    
    subroutine axis_init_dp(x,x0,x1,dx,nx)

        implicit none 

        real(dp), allocatable, intent(OUT) :: x(:)
        real(dp), optional :: x0
        real(dp), optional :: x1
        real(dp), optional :: dx
        integer,  optional :: nx 

        ! Local variables 
        integer :: i  
        real(dp) :: x0_now
        real(dp) :: x1_now
        real(dp) :: dx_now
        integer  :: nx_now 
        real(dp) :: nx_check

        dx_now = 1.0_dp 
        if (present(dx)) dx_now = dx 

        x0_now = 0.0_dp 
        if (present(x0)) x0_now = x0 
        
        ! Note: if x1 is present, nx is ignored 
        if (present(x1)) then 
            x1_now = x1 
        else if (present(nx)) then 
            x1_now = x0_now + (nx-1)*dx_now 
        else
            call fesm_error("axis_init", &
                "either x1 or nx must be present.")
        end if

        if (allocated(x)) deallocate(x)

        nx_now   = nint((x1_now-x0_now)/dx_now) + 1
        nx_check = (x1_now-x0_now)/dx_now + 1

        if ( abs(nx_now - nx_check) .gt. TOL ) then
            ! Make sure nx is a round number.
            call fesm_error("axis_init", &
                "desired axis bounds [x0,x1] do not divide evenly with dx.", &
                "x0 = "//to_str(x0_now)//new_line("a")// &
                "x1 = "//to_str(x1_now)//new_line("a")// &
                "dx = "//to_str(dx_now)//new_line("a")// &
                "nx = "//to_str((x1_now-x0_now)/dx_now + 1))
        end if
        
        allocate(x(nx_now))

        do i = 1, nx_now 
            x(i) = x0_now + dx_now*real(i-1,dp)
        end do 

        return

    end subroutine axis_init_dp

    subroutine print_var_range(var,name,missing_value,time)

        implicit none 

        real(wp),         intent(IN) :: var(:,:,:,:) 
        character(len=*), intent(IN) :: name
        real(wp),         intent(IN) :: missing_value 
        real(dp), intent(IN), optional :: time 

        ! Local variables 
        real(wp) :: vmin, vmax 

        if (count(var.ne.missing_value) .gt. 0) then 
            vmin = minval(var,mask=var.ne.missing_value)
            vmax = maxval(var,mask=var.ne.missing_value)
        else 
            vmin = missing_value 
            vmax = missing_value 
        end if 

        if (present(time)) then 
            write(*,"(f10.1,2x,a16,a3,2f14.3)") time, trim(name), ": ", vmin, vmax
        else 
            write(*,"(10x,2x,a16,a3,2f14.3)") trim(name), ": ", vmin, vmax
        end if 

        return 

    end subroutine print_var_range

    subroutine get_matching_files(filenames, pattern)

        use, intrinsic :: iso_c_binding, only : c_int

        implicit none

        character(len=*), allocatable, intent(OUT) :: filenames(:)
        character(len=*), intent(in) :: pattern

        ! Local variables
        character(len=1024) :: command
        character(len=256) :: temp_filename
        integer :: i, ios, unit
        integer :: num_files

        integer, parameter :: max_num_files_allowed = 10000

        ! Bind the POSIX getpid() from libc directly. gfortran offers a
        ! getpid() intrinsic extension, but ifx does not (it lives in the
        ! non-standard IFPORT module); the C interop route is portable.
        interface
            function c_getpid() bind(c, name="getpid")
                import :: c_int
                integer(c_int) :: c_getpid
            end function c_getpid
        end interface

        ! Temporary file to store file names. Make the name unique per
        ! process (PID) so concurrent runs sharing a working directory
        ! (e.g. batch experiments) do not clobber each other's list.
        write(temp_filename,"(a,i0,a)") "filelist_", c_getpid(), ".tmp"

        ! Create command to list files sorted alphabetically
        command = "ls -1 " // trim(pattern) // " | sort > " // trim(temp_filename)
        call execute_command_line(trim(command))

        ! Open the temporary file
        open(newunit=unit, file=trim(temp_filename), status="old", action="read", iostat=ios)
        
        ! Make sure at least one file was found
        if (ios /= 0) then
            call fesm_error("get_matching_files", &
                "temporary file could not be opened.", &
                "temp_filename = "//trim(temp_filename)//new_line("a")// &
                "iostat        = "//to_str(ios))
        end if

        ! First, count the number of files
        num_files = 0
        
        do i = 1, max_num_files_allowed
            read(unit, '(a)', iostat=ios)
            if (ios /= 0) exit
            num_files = num_files + 1
        end do
        rewind(unit)

        ! Make sure at least one file was found
        if (num_files .eq. 0) then
            call fesm_error("get_matching_files", &
                "filename(s) not found.", &
                "filename = "//trim(pattern))
        end if

        ! Allocate the output array
        if (allocated(filenames)) deallocate(filenames)
        allocate(filenames(num_files))

        ! Read the filenames into the array
        do i = 1, num_files
            read(unit, '(a)', iostat=ios) filenames(i)
        end do
        close(unit)

        ! Remove the temporary file
        call execute_command_line("rm " // trim(temp_filename))

        return

    end subroutine get_matching_files


    ! === Extensions to nc_read to handle multiple filenames ====

    ! Should this eventually be incorporated into ncio itself??
    ! It might be hard to do, since it would have to propagate through all routines,
    ! not just nc4_read_internal, but also nc_dims, etc...
    ! Now routine benefits from generic var[4D], but perhaps harder to collapse to 1D reading.

    subroutine nc_read_multifile(filenames,name,var,start,count,nt_files,missing_value)

        implicit none

        character (len=*),  intent(IN)      :: filenames(:)
        character (len=*),  intent(IN)      :: name
        real(wp),           intent(INOUT)   :: var(:,:,:,:)
        integer,            intent(IN)      :: start(:)
        integer,            intent(IN)      :: count(:)
        integer,            intent(IN)      :: nt_files(:)   ! per-file time length (cached)
        real(wp),           intent(IN), optional :: missing_value

        ! Local variables
        integer :: i, k0, nk, j0, nj, t0, nt, num_files
        integer :: ndim
        character(len=:), allocatable :: fnames

        ! Get number of dimensions we are working with
        ndim = size(start,1)

        num_files = size(filenames)

        ! Consistency check
        if (sum(nt_files) .lt. count(ndim)) then
            fnames = "per-file nt:"
            do i = 1, num_files
                fnames = fnames//new_line("a")// &
                    "  "//trim(filenames(i))//"  "//to_str(nt_files(i))
            end do
            call fesm_error("nc_read_multifile", &
                "number of time axis values read in is not sufficient to cover count.", &
                "count        = "//to_str(count)//new_line("a")// &
                "sum(nt_files) = "//to_str(sum(nt_files))//new_line("a")// &
                fnames)
        end if

        ! To do: figure out where index k0 begins within nt_files.
        k0 = start(size(start,1))
        nk = 0
        j0 = 0
        nt = 0
        do i = 1, num_files
            nk = nk + nt_files(i)       ! count maximum available including this file
            !nk = min(nk,count(ndim))       ! Limit maximum to total count in case it is less
            !write(*,*) "k0: ", k0, nk
            if (k0 .le. nk) then 
                j0 = k0 - j0            ! start index for current file
                nj = nk - k0 + 1        ! count available from current file
                nj = min(nj,count(ndim)-nt)  ! limit to values still needed (not total)
                nt = nt + nj            ! store total to be loaded so far
                t0 = nt - nj + 1        ! start index in output array
                !write(*,*) "j: ", j0, nj, k0, nk, t0, nt

                select case(ndim)

                    case(1)
                        call nc_read(filenames(i),name,var(t0:nt,1,1,1),missing_value=missing_value, &
                                                start=[j0],count=[nj])
                    case(2)
                        call nc_read(filenames(i),name,var(:,t0:nt,1,1),missing_value=missing_value, &
                                                start=[start(1),j0],count=[count(1),nj])
                    case(3)
                        call nc_read(filenames(i),name,var(:,:,t0:nt,1),missing_value=missing_value, &
                                                start=[start(1),start(2),j0],count=[count(1),count(2),nj])
                    case(4)
                        call nc_read(filenames(i),name,var(:,:,:,t0:nt),missing_value=missing_value, &
                                                start=[start(1),start(2),start(3),j0],count=[count(1),count(2),count(3),nj])
                end select

                k0 = k0 + nj            ! start index for whole dimension over all files
            end if
            j0 = nk                 ! reset j0 index to end of current total

            if (nt .ge. count(ndim)) exit
        end do

        if (nt .ne. count(ndim)) then
            call fesm_error("nc_read_multifile", &
                "number of time axis values read in does not match the expected total.", &
                "count = "//to_str(count)//new_line("a")// &
                "nk    = "//to_str(nk)//new_line("a")// &
                "nt    = "//to_str(nt))
        end if

        return

    end subroutine nc_read_multifile

end module varslice
