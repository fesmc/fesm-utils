module coupler
    ! A thin multigrid coupling layer over the coords mapping core.
    !
    ! A coupler_class holds two things:
    !   * a registry of named grids (grid_class objects, keyed by name), and
    !   * a lazily-populated cache of grid->grid maps (map_class objects).
    !
    ! Grids are nodes; maps are directed, method-typed edges. remap(cpl, ...)
    ! looks up (or builds+caches on first use) the map for a given
    ! (source, target, method) and applies it, sizing the destination array
    ! itself. Each map is built once, stored once, and shared across all callers
    ! -- the map cache uses fixed-capacity storage + a counter so that adding a
    ! map never reallocates (and so never deep-copies the large weight arrays of
    ! the maps already built).
    !
    ! method selects the interpolation kernel and is baked into the map at build
    ! time, so it is part of the cache key: "con" (conservative/area-weighted,
    ! the default) or a distance kernel (nn / shepard / quadrant / bilin). stat
    ! (mean / count / stdev) is applied when the field is mapped and is NOT part
    ! of the key -- it is forwarded to map_field per call.
    !
    ! gen selects the map generator used when a map's weights must actually be
    ! built: "coords" (in-package, the default) or "cdo" (external cdo call). It
    ! is a build-time-only concern (irrelevant once the map is cached), so it is
    ! NOT part of the cache key. A coupler carries a default gen (set in
    ! coupler_init); a specific map can override it via coupler_prime(gen=),
    ! which is the point at which that map's weights are generated.

    use coords, only : dp, sp, grid_class, map_class, map_init, map_field, &
                       grid_cdo_read_desc

    implicit none
    private

    integer, parameter :: GRID_MAX = 16    ! max registered grids per coupler
    integer, parameter :: MAP_MAX  = 64    ! max cached maps per coupler

    type grid_entry
        character(len=256) :: name = ""
        type(grid_class)   :: grid
    end type

    type coupler_class
        type(grid_entry)   :: grids(GRID_MAX)
        integer            :: ngrids = 0
        type(map_class)    :: maps(MAP_MAX)
        integer            :: nmaps  = 0
        character(len=512) :: map_fldr = "maps"   ! disk cache dir for map_init
        character(len=32)  :: gen = "coords"       ! default map generator; override
                                                   ! per map via coupler_prime(gen=)
    end type

    interface remap
        ! One wrapper per (kind x rank); each sizes var_dst and calls the generic
        ! map_field, which dispatches on kind (dp / sp / int).
        module procedure remap_2D_dp, remap_2D_sp, remap_2D_int
        module procedure remap_3D_dp, remap_3D_sp, remap_3D_int
    end interface

    public :: coupler_class
    public :: coupler_init, coupler_add_grid, coupler_prime
    public :: remap

contains

    subroutine coupler_init(cpl, map_fldr, gen)
        ! Reset a coupler to an empty registry + cache. `gen` sets the coupler-wide
        ! default map generator ("coords" / "cdo"); it can be overridden per map at
        ! prime time via coupler_prime(gen=).
        type(coupler_class), intent(out) :: cpl
        character(len=*), intent(in), optional :: map_fldr
        character(len=*), intent(in), optional :: gen

        cpl%ngrids = 0
        cpl%nmaps  = 0
        cpl%map_fldr = "maps"
        if (present(map_fldr)) cpl%map_fldr = trim(map_fldr)
        cpl%gen = "coords"
        if (present(gen)) cpl%gen = trim(gen)
    end subroutine coupler_init

    subroutine coupler_add_grid(cpl, name, grid)
        ! Register a grid under `name`. Grid identity (and the map cache/disk
        ! filenames) flow through grid%name, so `name` should match grid%name.
        type(coupler_class), intent(inout) :: cpl
        character(len=*),    intent(in)    :: name
        type(grid_class),    intent(in)    :: grid

        if (find_grid(cpl, name) > 0) &
            call coupler_err("coupler_add_grid: grid '"//trim(name)//"' already registered.")
        if (cpl%ngrids >= GRID_MAX) &
            call coupler_err("coupler_add_grid: GRID_MAX exceeded (raise it in coupler.f90).")

        cpl%ngrids = cpl%ngrids + 1
        cpl%grids(cpl%ngrids)%name = trim(name)
        cpl%grids(cpl%ngrids)%grid = grid
        ! Align the grid's own name with the registry name so the map cache keys
        ! (map%name1/name2, taken from grid%name) match the names callers pass.
        cpl%grids(cpl%ngrids)%grid%name = trim(name)
    end subroutine coupler_add_grid

    subroutine coupler_prime(cpl, src, dst, method, gen)
        ! Eagerly build a map up front (fail fast, pay the build/disk-load cost
        ! at init rather than mid-timestep). Same find-or-build path as remap.
        ! `gen` overrides the coupler-wide default generator for this one map --
        ! prime is where a specific map's weights get generated, so it is the
        ! natural override point.
        type(coupler_class), intent(inout) :: cpl
        character(len=*),    intent(in)    :: src, dst
        character(len=*), intent(in), optional :: method, gen

        integer :: im
        im = get_map(cpl, src, dst, method_or_default(method), gen)
    end subroutine coupler_prime

    ! ----- internals ---------------------------------------------------------

    function find_grid(cpl, name) result(ig)
        ! Registry index of `name`, or -1 if not registered.
        type(coupler_class), intent(in) :: cpl
        character(len=*),    intent(in) :: name
        integer :: ig, i

        ig = -1
        do i = 1, cpl%ngrids
            if (trim(cpl%grids(i)%name) == trim(name)) then
                ig = i
                return
            end if
        end do
    end function find_grid

    function get_map(cpl, src, dst, method, gen) result(im)
        ! Find the cached (src->dst, method) map, or build+cache it on miss.
        ! Grids are resolved by name: an in-memory registered grid takes
        ! precedence; otherwise the definition is read from grid_<name>.txt in
        ! the map folder (grid_cdo_read_desc). `gen` selects the generator on a
        ! build miss; when absent the coupler-wide default (cpl%gen) is used.
        type(coupler_class), intent(inout) :: cpl
        character(len=*),    intent(in)    :: src, dst, method
        character(len=*), intent(in), optional :: gen
        integer :: im

        integer :: i
        character(len=32) :: g
        type(grid_class) :: grid_src, grid_dst

        ! Cache hit: registered/disk grids both carry grid%name == the name the
        ! caller passes, so map%name1/name2 can be matched against src/dst directly.
        do i = 1, cpl%nmaps
            if (trim(cpl%maps(i)%name1)  == trim(src)    .and. &
                trim(cpl%maps(i)%name2)  == trim(dst)    .and. &
                trim(cpl%maps(i)%method) == trim(method)) then
                im = i
                return
            end if
        end do

        if (cpl%nmaps >= MAP_MAX) &
            call coupler_err("remap: MAP_MAX exceeded (raise it in coupler.f90).")

        call resolve_grid(cpl, src, grid_src)
        call resolve_grid(cpl, dst, grid_dst)

        g = cpl%gen
        if (present(gen)) g = trim(gen)

        cpl%nmaps = cpl%nmaps + 1
        im = cpl%nmaps
        call map_init(cpl%maps(im), grid_src, grid_dst, &
                      method=trim(method), gen=trim(g), fldr=trim(cpl%map_fldr), load=.true.)
    end function get_map

    subroutine resolve_grid(cpl, name, grid)
        ! Resolve a grid name to a grid_class: use the in-memory registry if the
        ! name is registered, else read grid_<name>.txt from the map folder.
        type(coupler_class), intent(in)  :: cpl
        character(len=*),    intent(in)  :: name
        type(grid_class),    intent(out) :: grid
        integer :: ig

        ig = find_grid(cpl, name)
        if (ig > 0) then
            grid = cpl%grids(ig)%grid
        else
            call grid_cdo_read_desc(grid, trim(name), trim(cpl%map_fldr))
        end if
    end subroutine resolve_grid

    function method_or_default(method) result(mtd)
        character(len=*), intent(in), optional :: method
        character(len=32) :: mtd
        mtd = "con"
        if (present(method)) mtd = trim(method)
    end function method_or_default

    subroutine coupler_err(msg)
        character(len=*), intent(in) :: msg
        write(*,*) "coupler:: "//trim(msg)
        stop 1
    end subroutine coupler_err

    ! ----- remap: 2D ---------------------------------------------------------

    subroutine remap_2D_dp(cpl, var_src, src, var_dst, dst, method, stat)
        type(coupler_class),   intent(inout) :: cpl
        real(dp),              intent(in)    :: var_src(:,:)
        character(len=*),      intent(in)    :: src, dst
        real(dp), allocatable, intent(inout) :: var_dst(:,:)
        character(len=*), intent(in), optional :: method, stat
        integer :: im

        im = get_map(cpl, src, dst, method_or_default(method))
        call alloc_2D_dp(var_dst, cpl%maps(im)%G%nx, cpl%maps(im)%G%ny)
        call map_field(cpl%maps(im), "coupler", var_src, var_dst, stat=stat)
    end subroutine remap_2D_dp

    subroutine remap_2D_sp(cpl, var_src, src, var_dst, dst, method, stat)
        type(coupler_class),   intent(inout) :: cpl
        real(sp),              intent(in)    :: var_src(:,:)
        character(len=*),      intent(in)    :: src, dst
        real(sp), allocatable, intent(inout) :: var_dst(:,:)
        character(len=*), intent(in), optional :: method, stat
        integer :: im

        im = get_map(cpl, src, dst, method_or_default(method))
        call alloc_2D_sp(var_dst, cpl%maps(im)%G%nx, cpl%maps(im)%G%ny)
        call map_field(cpl%maps(im), "coupler", var_src, var_dst, stat=stat)
    end subroutine remap_2D_sp

    subroutine remap_2D_int(cpl, var_src, src, var_dst, dst, method, stat)
        type(coupler_class),  intent(inout) :: cpl
        integer,              intent(in)    :: var_src(:,:)
        character(len=*),     intent(in)    :: src, dst
        integer, allocatable, intent(inout) :: var_dst(:,:)
        character(len=*), intent(in), optional :: method, stat
        integer :: im

        im = get_map(cpl, src, dst, method_or_default(method))
        call alloc_2D_int(var_dst, cpl%maps(im)%G%nx, cpl%maps(im)%G%ny)
        call map_field(cpl%maps(im), "coupler", var_src, var_dst, stat=stat)
    end subroutine remap_2D_int

    ! ----- remap: 3D (loop the 2D map over the trailing dimension) ------------

    subroutine remap_3D_dp(cpl, var_src, src, var_dst, dst, method, stat)
        type(coupler_class),   intent(inout) :: cpl
        real(dp),              intent(in)    :: var_src(:,:,:)
        character(len=*),      intent(in)    :: src, dst
        real(dp), allocatable, intent(inout) :: var_dst(:,:,:)
        character(len=*), intent(in), optional :: method, stat
        integer :: im, nz, k

        im = get_map(cpl, src, dst, method_or_default(method))
        nz = size(var_src, 3)
        call alloc_3D_dp(var_dst, cpl%maps(im)%G%nx, cpl%maps(im)%G%ny, nz)
        do k = 1, nz
            call map_field(cpl%maps(im), "coupler", var_src(:,:,k), var_dst(:,:,k), stat=stat)
        end do
    end subroutine remap_3D_dp

    subroutine remap_3D_sp(cpl, var_src, src, var_dst, dst, method, stat)
        type(coupler_class),   intent(inout) :: cpl
        real(sp),              intent(in)    :: var_src(:,:,:)
        character(len=*),      intent(in)    :: src, dst
        real(sp), allocatable, intent(inout) :: var_dst(:,:,:)
        character(len=*), intent(in), optional :: method, stat
        integer :: im, nz, k

        im = get_map(cpl, src, dst, method_or_default(method))
        nz = size(var_src, 3)
        call alloc_3D_sp(var_dst, cpl%maps(im)%G%nx, cpl%maps(im)%G%ny, nz)
        do k = 1, nz
            call map_field(cpl%maps(im), "coupler", var_src(:,:,k), var_dst(:,:,k), stat=stat)
        end do
    end subroutine remap_3D_sp

    subroutine remap_3D_int(cpl, var_src, src, var_dst, dst, method, stat)
        type(coupler_class),  intent(inout) :: cpl
        integer,              intent(in)    :: var_src(:,:,:)
        character(len=*),     intent(in)    :: src, dst
        integer, allocatable, intent(inout) :: var_dst(:,:,:)
        character(len=*), intent(in), optional :: method, stat
        integer :: im, nz, k

        im = get_map(cpl, src, dst, method_or_default(method))
        nz = size(var_src, 3)
        call alloc_3D_int(var_dst, cpl%maps(im)%G%nx, cpl%maps(im)%G%ny, nz)
        do k = 1, nz
            call map_field(cpl%maps(im), "coupler", var_src(:,:,k), var_dst(:,:,k), stat=stat)
        end do
    end subroutine remap_3D_int

    ! ----- allocate-on-demand helpers (reshape only when needed) -------------

    subroutine alloc_2D_dp(v, nx, ny)
        real(dp), allocatable, intent(inout) :: v(:,:)
        integer, intent(in) :: nx, ny
        if (allocated(v)) then
            if (size(v,1) /= nx .or. size(v,2) /= ny) deallocate(v)
        end if
        if (.not. allocated(v)) allocate(v(nx,ny))
    end subroutine alloc_2D_dp

    subroutine alloc_2D_sp(v, nx, ny)
        real(sp), allocatable, intent(inout) :: v(:,:)
        integer, intent(in) :: nx, ny
        if (allocated(v)) then
            if (size(v,1) /= nx .or. size(v,2) /= ny) deallocate(v)
        end if
        if (.not. allocated(v)) allocate(v(nx,ny))
    end subroutine alloc_2D_sp

    subroutine alloc_2D_int(v, nx, ny)
        integer, allocatable, intent(inout) :: v(:,:)
        integer, intent(in) :: nx, ny
        if (allocated(v)) then
            if (size(v,1) /= nx .or. size(v,2) /= ny) deallocate(v)
        end if
        if (.not. allocated(v)) allocate(v(nx,ny))
    end subroutine alloc_2D_int

    subroutine alloc_3D_dp(v, nx, ny, nz)
        real(dp), allocatable, intent(inout) :: v(:,:,:)
        integer, intent(in) :: nx, ny, nz
        if (allocated(v)) then
            if (size(v,1) /= nx .or. size(v,2) /= ny .or. size(v,3) /= nz) deallocate(v)
        end if
        if (.not. allocated(v)) allocate(v(nx,ny,nz))
    end subroutine alloc_3D_dp

    subroutine alloc_3D_sp(v, nx, ny, nz)
        real(sp), allocatable, intent(inout) :: v(:,:,:)
        integer, intent(in) :: nx, ny, nz
        if (allocated(v)) then
            if (size(v,1) /= nx .or. size(v,2) /= ny .or. size(v,3) /= nz) deallocate(v)
        end if
        if (.not. allocated(v)) allocate(v(nx,ny,nz))
    end subroutine alloc_3D_sp

    subroutine alloc_3D_int(v, nx, ny, nz)
        integer, allocatable, intent(inout) :: v(:,:,:)
        integer, intent(in) :: nx, ny, nz
        if (allocated(v)) then
            if (size(v,1) /= nx .or. size(v,2) /= ny .or. size(v,3) /= nz) deallocate(v)
        end if
        if (.not. allocated(v)) allocate(v(nx,ny,nz))
    end subroutine alloc_3D_int

end module coupler
