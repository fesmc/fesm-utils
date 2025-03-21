program test

    !use mapping_scrip
    use gaussian_quadrature

    type(gq2D_class) :: gq2D
    type(gq3D_class) :: gq3D

    call gq2D_init(gq2D)

    write(*,*)
    write(*,*) "----"
    write(*,*)

    call gq3D_init(gq3D)

end program test
