program test

    call test_derivatives()
    !call test_gaussian_quadrature()

contains

    subroutine test_derivatives()
        ! Test program adapted from: http://www.chides.org/APPM5720/html/FORTRAN_FD_EXAMPLE.html

        use precision
        use derivatives

        implicit none

        real(kind = 8), parameter :: pi = acos(-1.d0)
        integer :: n,i,j
        real(wp) :: h
        real(wp), allocatable :: x(:),f(:),df(:)
        real(dp), allocatable :: df_exact(:)
        logical, allocatable :: m(:)

        do n  = 4,100
            ! Allocate memory for the various arrays
            allocate(x(0:n),f(0:n),df(0:n),df_exact(0:n))
            allocate(m(0:n))

            m = .TRUE. 

            ! Set up the grid.
            h = 2.d0/dble(n)
            do i = 0,n
                x(i) = -1.d0+dble(i)*h
            end do
            ! The function and the exact derivative
            f = exp(cos(pi*x))
            df_exact = -pi*sin(pi*x)*exp(cos(pi*x))

            ! Calculatea derivative
            call calc_dvdx_1D(df,f,h,m,bc="neumann")

            write(*,'(I3,2(ES12.4))') n,h,maxval(abs(df-df_exact))

            ! Deallocate the arrays
            deallocate(x,f,df,df_exact,m)
        end do

    end subroutine test_derivatives


    subroutine test_gaussian_quadrature()

        use precision
        
        !use mapping_scrip
        use gaussian_quadrature

        type(gq2D_class) :: gq2D
        type(gq3D_class) :: gq3D

        real(wp) :: var_qp(4)
        real(wp) :: var(2,2)
        real(wp) :: dx
        real(wp) :: dy 
        integer  :: i, j 
        integer  :: im1, ip1, jm1, jp1

        logical, parameter :: test_2D = .TRUE.
        logical, parameter :: test_3D = .FALSE.
        
        if (test_2D) then
            call gq2D_init(gq2D,verbose=.TRUE.)

            var(1,1) = 2.0
            var(1,2) = 2.0
            var(2,1) = 4.0
            var(2,2) = 4.0

            dx = 1.0
            dy = 1.0 
            
            i = 2
            j = 2
            im1 = i-1
            ip1 = i
            jm1 = j-1
            jp1 = j

            call gq2D_to_nodes(gq2D, var_qp, var, dx, dy, "ab",i,j,im1,ip1,jm1,jp1)

            write(*,*) "var_qp: ", var_qp
        
        end if

        write(*,*)
        write(*,*) "----"
        write(*,*)

        if (test_3D) then
            call gq3D_init(gq3D,verbose=.TRUE.)
        end if

        return

    end subroutine test_gaussian_quadrature

end program test
