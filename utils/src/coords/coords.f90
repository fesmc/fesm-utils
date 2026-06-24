module coords
    ! Umbrella module for the coords library. `use coords` gives access to the
    ! public geometry, projection, planet and geodesic interfaces. Additional
    ! submodules (mapping, conservative) are re-exported here as they land.

    use coord_constants
    use geodesic
    use planet
    use oblimap_projection_module
    use coordinates

    implicit none

    public

end module coords
