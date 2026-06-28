module ncio_interp
    ! Read a field from a NetCDF file and, when its grid does not match the target
    ! array, remap it onto the target with the coords mapping core (map_class /
    ! map_field). This is the ncio x mapping bridge: it was previously embedded in
    ! mapping_scrip as the nc_read_interp family and applied with the legacy SCRIP
    ! applier; it now lives here and applies with the unified map_field.
    !
    ! nc_read_interp(filename, vnm, var, [var_in], [ncid], [start], [count],
    !                [map], [stat], [fill_method], [filt_method], [filt_par])
    !
    !   - var_in present : take the source field from the argument (used by the 3D
    !                      variants to feed one level at a time) instead of reading.
    !   - source dims == target dims : copy through, no map needed.
    !   - source dims /= target dims : a map_class must be supplied; remap with it.
    !
    ! Variants: dp/sp/integer in 2D and 3D, plus logical 2D. The 3D variants loop
    ! over the third dimension applying the 2D remap (third dim is not interpolated).

    use, intrinsic :: iso_fortran_env, only : error_unit

    use precision,       only : dp, sp
    use constants, only : MV => mv_dp
    use ncio
    use mapping,         only : map_class, map_field

    implicit none

    private

    interface nc_read_interp
        module procedure nc_read_interp_dp_2D
        module procedure nc_read_interp_dp_3D
        module procedure nc_read_interp_sp_2D
        module procedure nc_read_interp_sp_3D
        module procedure nc_read_interp_int_2D
        module procedure nc_read_interp_int_3D
        module procedure nc_read_interp_logical_2D
    end interface

    public :: nc_read_interp

contains

    subroutine nc_read_interp_dp_2D(filename,vnm,var2D,var2D_in,ncid,start,count,map,stat, &
                                        fill_method,filt_method,filt_par)
        ! Read in a 2D field and remap it onto var2D if needed.

        implicit none

        character(len=*),   intent(IN)  :: filename
        character(len=*),   intent(IN)  :: vnm
        real(dp),           intent(OUT) :: var2D(:,:)
        real(dp),optional,  intent(IN)  :: var2D_in(:,:)
        integer, optional,  intent(IN)  :: ncid
        integer, optional,  intent(IN)  :: start(:)
        integer, optional,  intent(IN)  :: count(:)
        type(map_class),  intent(IN),  optional :: map
        character(len=*), intent(IN),  optional :: stat            ! Aggregation: mean (default), count, stdev
        character(len=*), intent(IN),  optional :: fill_method     ! Method to fill in remaining missing values
        character(len=*), intent(IN),  optional :: filt_method     ! Method to use for filtering
        real(dp),         intent(IN),  optional :: filt_par(:)     ! gaussian=[sigma,dx]; poisson=[tol]

        ! Local variables
        integer :: nx, ny
        integer,  allocatable :: dims(:)
        real(dp), allocatable :: var2D_src(:,:)
        character(len=56) :: mapping_method

        if (present(var2D_in)) then
            ! Get source array from argument (useful for handling 3D arrays with this routine)
            nx = size(var2D_in,1)
            ny = size(var2D_in,2)
            allocate(var2D_src(nx,ny))
            var2D_src = var2D_in
        else
            ! Load 2D array from file
            call nc_dims(filename,vnm,dims=dims)
            nx = dims(1)
            ny = dims(2)
            allocate(var2D_src(nx,ny))
            call nc_read(filename,vnm,var2D_src,ncid=ncid,start=start,count=[nx,ny,1],missing_value=real(MV,dp))
        end if

        if (nx .eq. size(var2D,1) .and. ny .eq. size(var2D,2)) then
            ! Assume no interpolation needed, copy variable for output directly
            var2D = var2D_src
        else
            ! Map local source array onto our target array
            mapping_method = "mean"
            if (present(stat)) mapping_method = trim(stat)

            ! Safety check
            if (.not. present(map)) then
                write(error_unit,*) ""
                write(error_unit,*) "nc_read_interp:: Error: map_class object must &
                        &be provided as an argument since array read from file does not &
                        &match the target array size."
                write(error_unit,*) "filename: ", trim(filename)
                write(error_unit,*) "variable: ", trim(vnm)
                write(error_unit,*) "dims in file:          ", nx, ny
                write(error_unit,*) "dims in target object: ", size(var2D,1), size(var2D,2)
                stop
            end if

            ! Perform interpolation
            var2D = real(MV,dp)
            call map_field(map,vnm,var2D_src,var2D,stat=mapping_method,missing_value=real(MV,dp), &
                        fill_method=fill_method,filt_method=filt_method,filt_par=filt_par)
        end if

        return

    end subroutine nc_read_interp_dp_2D

    subroutine nc_read_interp_dp_3D(filename,vnm,var3D,ncid,start,count,map,stat, &
                                        fill_method,filt_method,filt_par)
        ! Read in a 3D field and remap each level onto var3D if needed.

        implicit none

        character(len=*),   intent(IN)  :: filename
        character(len=*),   intent(IN)  :: vnm
        real(dp),           intent(OUT) :: var3D(:,:,:)
        integer, optional,  intent(IN)  :: ncid
        integer, optional,  intent(IN)  :: start(:)
        integer, optional,  intent(IN)  :: count(:)
        type(map_class),  intent(IN),  optional :: map
        character(len=*), intent(IN),  optional :: stat            ! Aggregation: mean (default), count, stdev
        character(len=*), intent(IN),  optional :: fill_method
        character(len=*), intent(IN),  optional :: filt_method
        real(dp),         intent(IN),  optional :: filt_par(:)

        ! Local variables
        integer :: nx, ny, nz, k
        integer,  allocatable :: dims(:)
        real(dp), allocatable :: var3D_src(:,:,:)

        call nc_dims(filename,vnm,dims=dims)
        nx = dims(1)
        ny = dims(2)
        nz = dims(3)

        allocate(var3D_src(nx,ny,nz))

        ! Safety check
        if (nz .ne. size(var3D,3)) then
            write(error_unit,*) ""
            write(error_unit,*) "nc_read_interp_dp_3D:: Error: third dimension of variable in &
                    &input file does not match third dimension of array provided as an argument. &
                    &Interpolation of the third dimension is not yet supported."
            write(error_unit,*) "filename  = ", trim(filename)
            write(error_unit,*) "variable  = ", trim(vnm)
            write(error_unit,*) "nz[file]  = ", nz
            write(error_unit,*) "nz[array] = ", size(var3D,3)
            stop
        end if

        call nc_read(filename,vnm,var3D_src,ncid=ncid,start=start,count=[nx,ny,nz,1],missing_value=real(MV,dp))

        do k = 1, nz
            call nc_read_interp_dp_2D(filename,vnm,var3D(:,:,k),var3D_src(:,:,k),ncid,start,count, &
                                map,stat,fill_method=fill_method,filt_method=filt_method,filt_par=filt_par)
        end do

        return

    end subroutine nc_read_interp_dp_3D

    subroutine nc_read_interp_sp_2D(filename,vnm,var2D,var2D_in,ncid,start,count,map,stat, &
                                        fill_method,filt_method,filt_par)
        ! Read in a 2D field and remap it onto var2D if needed.

        implicit none

        character(len=*),   intent(IN)  :: filename
        character(len=*),   intent(IN)  :: vnm
        real(sp),           intent(OUT) :: var2D(:,:)
        real(sp),optional,  intent(IN)  :: var2D_in(:,:)
        integer, optional,  intent(IN)  :: ncid
        integer, optional,  intent(IN)  :: start(:)
        integer, optional,  intent(IN)  :: count(:)
        type(map_class),  intent(IN),  optional :: map
        character(len=*), intent(IN),  optional :: stat            ! Aggregation: mean (default), count, stdev
        character(len=*), intent(IN),  optional :: fill_method
        character(len=*), intent(IN),  optional :: filt_method
        real(sp),         intent(IN),  optional :: filt_par(:)

        ! Local variables
        integer :: nx, ny
        integer,  allocatable :: dims(:)
        real(sp), allocatable :: var2D_src(:,:)
        character(len=56) :: mapping_method

        if (present(var2D_in)) then
            nx = size(var2D_in,1)
            ny = size(var2D_in,2)
            allocate(var2D_src(nx,ny))
            var2D_src = var2D_in
        else
            call nc_dims(filename,vnm,dims=dims)
            nx = dims(1)
            ny = dims(2)
            allocate(var2D_src(nx,ny))
            call nc_read(filename,vnm,var2D_src,ncid=ncid,start=start,count=[nx,ny,1],missing_value=real(MV,sp))
        end if

        if (nx .eq. size(var2D,1) .and. ny .eq. size(var2D,2)) then
            var2D = var2D_src
        else
            mapping_method = "mean"
            if (present(stat)) mapping_method = trim(stat)

            if (.not. present(map)) then
                write(error_unit,*) ""
                write(error_unit,*) "nc_read_interp:: Error: map_class object must &
                        &be provided as an argument since array read from file does not &
                        &match the target array size."
                write(error_unit,*) "filename: ", trim(filename)
                write(error_unit,*) "variable: ", trim(vnm)
                write(error_unit,*) "dims in file:          ", nx, ny
                write(error_unit,*) "dims in target object: ", size(var2D,1), size(var2D,2)
                stop
            end if

            var2D = real(MV,sp)
            call map_field(map,vnm,var2D_src,var2D,stat=mapping_method,missing_value=real(MV,sp), &
                            fill_method=fill_method,filt_method=filt_method,filt_par=filt_par)
        end if

        return

    end subroutine nc_read_interp_sp_2D

    subroutine nc_read_interp_sp_3D(filename,vnm,var3D,ncid,start,count,map,stat, &
                                        fill_method,filt_method,filt_par)
        ! Read in a 3D field and remap each level onto var3D if needed.

        implicit none

        character(len=*),   intent(IN)  :: filename
        character(len=*),   intent(IN)  :: vnm
        real(sp),           intent(OUT) :: var3D(:,:,:)
        integer, optional,  intent(IN)  :: ncid
        integer, optional,  intent(IN)  :: start(:)
        integer, optional,  intent(IN)  :: count(:)
        type(map_class),  intent(IN),  optional :: map
        character(len=*), intent(IN),  optional :: stat            ! Aggregation: mean (default), count, stdev
        character(len=*), intent(IN),  optional :: fill_method
        character(len=*), intent(IN),  optional :: filt_method
        real(sp),         intent(IN),  optional :: filt_par(:)

        ! Local variables
        integer :: nx, ny, nz, k
        integer,  allocatable :: dims(:)
        real(sp), allocatable :: var3D_src(:,:,:)

        call nc_dims(filename,vnm,dims=dims)
        nx = dims(1)
        ny = dims(2)
        nz = dims(3)

        allocate(var3D_src(nx,ny,nz))

        if (nz .ne. size(var3D,3)) then
            write(error_unit,*) ""
            write(error_unit,*) "nc_read_interp_sp_3D:: Error: third dimension of variable in &
                    &input file does not match third dimension of array provided as an argument. &
                    &Interpolation of the third dimension is not yet supported."
            write(error_unit,*) "filename  = ", trim(filename)
            write(error_unit,*) "variable  = ", trim(vnm)
            write(error_unit,*) "nz[file]  = ", nz
            write(error_unit,*) "nz[array] = ", size(var3D,3)
            stop
        end if

        call nc_read(filename,vnm,var3D_src,ncid=ncid,start=start,count=[nx,ny,nz,1],missing_value=real(MV,sp))

        do k = 1, nz
            call nc_read_interp_sp_2D(filename,vnm,var3D(:,:,k),var3D_src(:,:,k),ncid,start,count,map,stat, &
                                fill_method=fill_method,filt_method=filt_method,filt_par=filt_par)
        end do

        return

    end subroutine nc_read_interp_sp_3D

    subroutine nc_read_interp_int_2D(filename,vnm,var2D,var2D_in,ncid,start,count,map,stat, &
                                        fill_method)
        ! Read in a 2D integer field and remap it onto var2D if needed.

        implicit none

        character(len=*),   intent(IN)  :: filename
        character(len=*),   intent(IN)  :: vnm
        integer,            intent(OUT) :: var2D(:,:)
        integer, optional,  intent(IN)  :: var2D_in(:,:)
        integer, optional,  intent(IN)  :: ncid
        integer, optional,  intent(IN)  :: start(:)
        integer, optional,  intent(IN)  :: count(:)
        type(map_class),  intent(IN),  optional :: map
        character(len=*), intent(IN),  optional :: stat            ! Aggregation: mean (default), count, stdev
        character(len=*), intent(IN),  optional :: fill_method

        ! Local variables
        integer :: nx, ny
        integer,  allocatable :: dims(:)
        integer,  allocatable :: var2D_src(:,:)
        character(len=56) :: mapping_method

        if (present(var2D_in)) then
            nx = size(var2D_in,1)
            ny = size(var2D_in,2)
            allocate(var2D_src(nx,ny))
            var2D_src = var2D_in
        else
            call nc_dims(filename,vnm,dims=dims)
            nx = dims(1)
            ny = dims(2)
            allocate(var2D_src(nx,ny))
            call nc_read(filename,vnm,var2D_src,ncid=ncid,start=start,count=[nx,ny,1],missing_value=int(MV))
        end if

        if (nx .eq. size(var2D,1) .and. ny .eq. size(var2D,2)) then
            var2D = var2D_src
        else
            mapping_method = "mean"
            if (present(stat)) mapping_method = trim(stat)

            if (.not. present(map)) then
                write(error_unit,*) ""
                write(error_unit,*) "nc_read_interp:: Error: map_class object must &
                        &be provided as an argument since array read from file does not &
                        &match the target array size."
                write(error_unit,*) "filename: ", trim(filename)
                write(error_unit,*) "variable: ", trim(vnm)
                write(error_unit,*) "dims in file:          ", nx, ny
                write(error_unit,*) "dims in target object: ", size(var2D,1), size(var2D,2)
                stop
            end if

            var2D = int(MV)
            call map_field(map,vnm,var2D_src,var2D,stat=mapping_method, &
                                    missing_value=int(MV),fill_method=fill_method)
        end if

        return

    end subroutine nc_read_interp_int_2D

    subroutine nc_read_interp_int_3D(filename,vnm,var3D,ncid,start,count,map,stat,fill_method)
        ! Read in a 3D integer field and remap each level onto var3D if needed.

        implicit none

        character(len=*),   intent(IN)  :: filename
        character(len=*),   intent(IN)  :: vnm
        integer,            intent(OUT) :: var3D(:,:,:)
        integer, optional,  intent(IN)  :: ncid
        integer, optional,  intent(IN)  :: start(:)
        integer, optional,  intent(IN)  :: count(:)
        type(map_class),  intent(IN),  optional :: map
        character(len=*), intent(IN),  optional :: stat            ! Aggregation: mean (default), count, stdev
        character(len=*), intent(IN),  optional :: fill_method

        ! Local variables
        integer :: nx, ny, nz, k
        integer,  allocatable :: dims(:)
        integer,  allocatable :: var3D_src(:,:,:)

        call nc_dims(filename,vnm,dims=dims)
        nx = dims(1)
        ny = dims(2)
        nz = dims(3)

        allocate(var3D_src(nx,ny,nz))

        if (nz .ne. size(var3D,3)) then
            write(error_unit,*) ""
            write(error_unit,*) "nc_read_interp_int_3D:: Error: third dimension of variable in &
                    &input file does not match third dimension of array provided as an argument. &
                    &Interpolation of the third dimension is not yet supported."
            write(error_unit,*) "filename  = ", trim(filename)
            write(error_unit,*) "variable  = ", trim(vnm)
            write(error_unit,*) "nz[file]  = ", nz
            write(error_unit,*) "nz[array] = ", size(var3D,3)
            stop
        end if

        call nc_read(filename,vnm,var3D_src,ncid=ncid,start=start,count=[nx,ny,nz,1],missing_value=int(MV))

        do k = 1, nz
            call nc_read_interp_int_2D(filename,vnm,var3D(:,:,k),var3D_src(:,:,k),ncid,start,count,map,stat, &
                                fill_method=fill_method)
        end do

        return

    end subroutine nc_read_interp_int_3D

    subroutine nc_read_interp_logical_2D(filename,vnm,var2D,var2D_in,ncid,start,count,map,stat)
        ! Read in a 2D logical field and remap it onto var2D if needed. The remap
        ! is done on a 0/1 integer image (default method "count") and thresholded.

        implicit none

        character(len=*),   intent(IN)  :: filename
        character(len=*),   intent(IN)  :: vnm
        logical,            intent(OUT) :: var2D(:,:)
        logical, optional,  intent(IN)  :: var2D_in(:,:)
        integer, optional,  intent(IN)  :: ncid
        integer, optional,  intent(IN)  :: start(:)
        integer, optional,  intent(IN)  :: count(:)
        type(map_class),  intent(IN),  optional :: map
        character(len=*), intent(IN),  optional :: stat            ! Aggregation: mean (default), count, stdev

        ! Local variables
        integer :: nx, ny
        integer,  allocatable :: dims(:)
        integer,  allocatable :: var2D_src(:,:)
        integer,  allocatable :: var2D_int(:,:)
        character(len=56) :: mapping_method

        if (present(var2D_in)) then
            nx = size(var2D_in,1)
            ny = size(var2D_in,2)
            allocate(var2D_src(nx,ny))
            var2D_src = 0
            where (var2D_in) var2D_src = 1
        else
            call nc_dims(filename,vnm,dims=dims)
            nx = dims(1)
            ny = dims(2)
            allocate(var2D_src(nx,ny))
            call nc_read(filename,vnm,var2D_src,ncid=ncid,start=start,count=[nx,ny,1],missing_value=int(MV))
        end if

        if (nx .eq. size(var2D,1) .and. ny .eq. size(var2D,2)) then
            var2D = .FALSE.
            where(var2D_src .eq. 1) var2D = .TRUE.
        else
            mapping_method = "count"
            if (present(stat)) mapping_method = trim(stat)

            if (.not. present(map)) then
                write(error_unit,*) ""
                write(error_unit,*) "nc_read_interp:: Error: map_class object must &
                        &be provided as an argument since array read from file does not &
                        &match the target array size."
                write(error_unit,*) "filename: ", trim(filename)
                write(error_unit,*) "variable: ", trim(vnm)
                write(error_unit,*) "dims in file:          ", nx, ny
                write(error_unit,*) "dims in target object: ", size(var2D,1), size(var2D,2)
                stop
            end if

            allocate(var2D_int(size(var2D,1),size(var2D,2)))
            var2D_int = int(MV)
            call map_field(map,vnm,var2D_src,var2D_int,stat=mapping_method, &
                                                missing_value=int(MV),fill_method="nn")

            var2D = .FALSE.
            where(var2D_int .eq. 1) var2D = .TRUE.
        end if

        return

    end subroutine nc_read_interp_logical_2D

end module ncio_interp
