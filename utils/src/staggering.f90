module staggering

    use, intrinsic :: iso_fortran_env, only : input_unit, output_unit, error_unit    
    use precision

    implicit none

    private
    
contains

    function calc_magnitude_from_staggered(u,v,f_ice,boundaries) result(umag)
        ! Calculate the centered (aa-nodes) magnitude of a vector 
        ! from the staggered (ac-nodes) components

        implicit none 
        
        real(wp), intent(IN)  :: u(:,:), v(:,:)
        real(wp), intent(IN)  :: f_ice(:,:) 
        real(wp) :: umag(size(u,1),size(u,2)) 
        character(len=*), intent(IN) :: boundaries 

        ! Local variables 
        integer :: i, j, nx, ny 
        integer :: im1, ip1, jm1, jp1
        real(wp) :: unow, vnow 
        real(wp) :: f1, f2, H1, H2 
        
        nx = size(u,1)
        ny = size(u,2) 

        umag = 0.0_wp 

        do j = 1, ny 
        do i = 1, nx 

            ! Get neighbor indices
            call get_neighbor_indices(im1,ip1,jm1,jp1,i,j,nx,ny,boundaries)
            
            if (f_ice(i,j) .eq. 1.0) then 
                unow = 0.5*(u(im1,j)+u(i,j))
                vnow = 0.5*(v(i,jm1)+v(i,j))
                
                if (abs(unow) .lt. TOL_UNDERFLOW) unow = 0.0_wp 
                if (abs(vnow) .lt. TOL_UNDERFLOW) vnow = 0.0_wp 
            else 
                unow = 0.0 
                vnow = 0.0 
            end if 

            umag(i,j) = sqrt(unow*unow+vnow*vnow)

            if (abs(umag(i,j)) .lt. TOL_UNDERFLOW) umag(i,j) = 0.0_wp

        end do 
        end do 

        return

    end function calc_magnitude_from_staggered

    function stagger_ac_aa(u,v) result(umag)
        ! Calculate the centered (aa-node) magnitude of a scalar 
        ! from the staggered (ac-node) components

        implicit none 
        
        real(wp), intent(IN)  :: u(:,:), v(:,:)    ! acx-, acy-nodes 
        real(wp) :: umag(size(u,1),size(u,2))      ! aa-nodes 

        ! Local variables 
        integer :: i, j, nx, ny 
        integer :: im1, ip1, jm1, jp1 

        nx = size(u,1)
        ny = size(u,2) 

        umag = 0.0_wp 

        do j = 1, ny
        do i = 1, nx 
            ! BC: Periodic boundary conditions
            im1 = i-1
            if (im1 == 0) then
                im1 = nx
            end if
            ip1 = i+1
            if (ip1 == nx+1) then
                ip1 = 1
            end if

            jm1 = j-1
            if (jm1 == 0) then
                jm1 = ny
            end if
            jp1 = j+1
            if (jp1 == ny+1) then
                jp1 = 1
            end if

            umag(i,j) = 0.25_wp*(u(i,j)+u(im1,j)+v(i,j)+v(i,jm1))
        end do 
        end do 

        return

    end function stagger_ac_aa
    
    function stagger_aa_ab(u) result(ustag)
        ! Stagger from Aa => Ab
        ! Four point average from corner Aa nodes to central Ab node 

        implicit none 

        real(wp), intent(IN)  :: u(:,:) 
        real(wp) :: ustag(size(u,1),size(u,2)) 

        ! Local variables 
        integer :: i, j, nx, ny  
        integer :: im1, ip1, jm1, jp1 
        
        nx = size(u,1)
        ny = size(u,2) 

        ustag = 0.0_wp 

        do j = 1, ny 
        do i = 1, nx
            ! BC: Periodic boundary conditions
            ip1 = i+1
            if (ip1 == nx+1) then
                ip1 = 1
            end if
            jp1 = j+1
            if (jp1 == ny+1) then
                jp1 = 1
            end if

            ustag(i,j) = 0.25_wp*(u(ip1,jp1)+u(ip1,j)+u(i,jp1)+u(i,j))
        end do 
        end do 

        return

    end function stagger_aa_ab
    
    function stagger_aa_ab_ice(u,H_ice,f_ice) result(ustag)
        ! Stagger from Aa => Ab
        ! Four point average from corner Aa nodes to central Ab node 

        implicit none 

        real(wp), intent(IN)  :: u(:,:) 
        real(wp), intent(IN)  :: H_ice(:,:) 
        real(wp), intent(IN)  :: f_ice(:,:) 
        real(wp) :: ustag(size(u,1),size(u,2)) 

        ! Local variables 
        integer :: i, j, nx, ny, k   
        integer :: im1, ip1, jm1, jp1 

        nx = size(u,1)
        ny = size(u,2) 

        ustag = 0.0_wp 

        do j = 1, ny 
        do i = 1, nx

            ! BC: Periodic boundary conditions
            ip1 = i+1
            if (ip1 == nx+1) then
                ip1 = 1
            end if
            jp1 = j+1
            if (jp1 == ny+1) then
                jp1 = 1
            end if

            k = 0 
            ustag(i,j) = 0.0 
            if (f_ice(i,j) .eq. 0.0) then 
                ustag(i,j) = ustag(i,j) + u(i,j) 
                k = k+1
            end if 

            if (f_ice(ip1,j) .eq. 0.0) then 
                ustag(i,j) = ustag(i,j) + u(ip1,j) 
                k = k+1 
            end if 
            
            if (f_ice(i,jp1) .eq. 0.0) then 
                ustag(i,j) = ustag(i,j) + u(i,jp1) 
                k = k+1 
            end if 
            
            if (f_ice(ip1,jp1) .eq. 0.0) then 
                ustag(i,j) = ustag(i,j) + u(ip1,jp1) 
                k = k+1 
            end if 
            
            if (k .gt. 0) then 
                ustag(i,j) = ustag(i,j) / real(k,wp)
            else 
                ustag(i,j) = 0.25_wp*(u(ip1,jp1)+u(ip1,j)+u(i,jp1)+u(i,j))
            end if 

        end do 
        end do 

        return

    end function stagger_aa_ab_ice
    
    function stagger_ab_aa(u) result(ustag)
        ! Stagger from Ab => Aa
        ! Four point average from corner Ab nodes to central Aa node 

        implicit none 

        real(wp), intent(IN)  :: u(:,:) 
        real(wp) :: ustag(size(u,1),size(u,2)) 

        ! Local variables 
        integer :: i, j, nx, ny  
        integer :: im1, jm1, ip1, jp1

        nx = size(u,1)
        ny = size(u,2) 

        ustag = 0.0_wp 

        do j = 1, ny 
        do i = 1, nx
            ! BC: Periodic boundary conditions
            im1 = i-1
            if (im1 == 0) then
                im1 = nx
            end if
            jm1 = j-1
            if (jm1 == 0) then
                jm1 = ny
            end if

            ustag(i,j) = 0.25_wp*(u(i,j)+u(im1,j)+u(i,jm1)+u(im1,jm1))
        end do 
        end do 

        return

    end function stagger_ab_aa
    
    function stagger_ab_aa_ice(u,H_ice,f_ice) result(ustag)
        ! Stagger from ab => aa
        ! Four point average from corner ab-nodes to central aa-node 

        implicit none 

        real(wp), intent(IN)  :: u(:,:) 
        real(wp), intent(IN)  :: H_ice(:,:) 
        real(wp), intent(IN)  :: f_ice(:,:) 
        real(wp) :: ustag(size(u,1),size(u,2)) 

        ! Local variables 
        integer :: i, j, nx, ny, k   
        integer :: im1, jm1, ip1, jp1 
        real(wp), allocatable ::H_ice_ab(:,:) 

        nx = size(u,1)
        ny = size(u,2) 

        allocate(H_ice_ab(nx,ny))

        H_ice_ab = stagger_aa_ab_ice(H_ice,H_ice,f_ice) 

        ustag = 0.0_wp 

        do j = 1, ny 
        do i = 1, nx

            ! BC: Periodic boundary conditions
            im1 = i-1
            if (im1 == 0) then
                im1 = nx
            end if

            jm1 = j-1
            if (jm1 == 0) then
                jm1 = ny
            end if

            k = 0 
            ustag(i,j) = 0.0 
            if (H_ice_ab(i,j) .gt. 0.0) then 
                ustag(i,j) = ustag(i,j) + u(i,j) 
                k = k+1
            end if 

            if (H_ice_ab(im1,j) .gt. 0.0) then 
                ustag(i,j) = ustag(i,j) + u(im1,j) 
                k = k+1 
            end if 
            
            if (H_ice_ab(i,jm1) .gt. 0.0) then 
                ustag(i,j) = ustag(i,j) + u(i,jm1) 
                k = k+1 
            end if 
            
            if (H_ice_ab(im1,jm1) .gt. 0.0) then 
                ustag(i,j) = ustag(i,j) + u(im1,jm1) 
                k = k+1 
            end if 
            
            if (k .gt. 0) then 
                ustag(i,j) = ustag(i,j) / real(k,wp)
            else 
                ustag(i,j) = 0.25_wp*(u(im1,jm1)+u(im1,j)+u(i,jm1)+u(i,j))
            end if 

        end do 
        end do 

        return

    end function stagger_ab_aa_ice
    
    function stagger_aa_acx(u) result(ustag)
        ! Stagger from Aa => Ac, x-direction 

        implicit none

        real(wp), intent(IN)  :: u(:,:) 
        real(wp) :: ustag(size(u,1),size(u,2)) 

        ! Local variables 
        integer :: i, j, nx, ny  

        nx = size(u,1)
        ny = size(u,2) 

        ustag = 0.0_wp 

        do j = 1, ny 
        do i = 1, nx-1
            ustag(i,j) = 0.5_wp*(u(i,j)+u(i+1,j))
        end do 
        end do 

        return

    end function stagger_aa_acx
    
    function stagger_aa_acy(u) result(ustag)
        ! Stagger from Aa => Ac 

        implicit none 

        real(wp), intent(IN)  :: u(:,:) 
        real(wp) :: ustag(size(u,1),size(u,2)) 

        ! Local variables 
        integer :: i, j, nx, ny  

        nx = size(u,1)
        ny = size(u,2) 

        ustag = 0.0_wp 

        do j = 1, ny-1 
        do i = 1, nx
            ustag(i,j) = 0.5_wp*(u(i,j)+u(i,j+1))
        end do 
        end do 

        return

    end function stagger_aa_acy
    
    function stagger_acx_aa(u) result(ustag)
        ! Stagger from Aa => Ac, x-direction 

        implicit none

        real(wp), intent(IN)  :: u(:,:) 
        real(wp) :: ustag(size(u,1),size(u,2)) 

        ! Local variables 
        integer :: i, j, nx, ny  

        nx = size(u,1)
        ny = size(u,2) 

        ustag = 0.0_wp 

        do j = 1, ny 
        do i = 2, nx
            ustag(i,j) = 0.5_wp*(u(i-1,j)+u(i,j))
        end do 
        end do 

        return

    end function stagger_acx_aa
    
    function stagger_acy_aa(u) result(ustag)
        ! Stagger from Aa => Ac 

        implicit none 

        real(wp), intent(IN)  :: u(:,:) 
        real(wp) :: ustag(size(u,1),size(u,2)) 

        ! Local variables 
        integer :: i, j, nx, ny  

        nx = size(u,1)
        ny = size(u,2) 

        ustag = 0.0_wp 

        do j = 2, ny
        do i = 1, nx
            ustag(i,j) = 0.5_wp*(u(i,j-1)+u(i,j))
        end do 
        end do 

        return

    end function stagger_acy_aa
    
    function stagger_ab_acx(u) result(ustag)
        ! Stagger from Ab => Ac, x-direction 

        implicit none

        real(wp), intent(IN)  :: u(:,:) 
        real(wp) :: ustag(size(u,1),size(u,2)) 

        ! Local variables 
        integer :: i, j, nx, ny  

        nx = size(u,1)
        ny = size(u,2) 

        ustag = 0.0_wp 

        do j = 2, ny 
        do i = 1, nx
            ustag(i,j) = 0.5_wp*(u(i,j)+u(i,j-1))
        end do 
        end do 

        return

    end function stagger_ab_acx
    
    function stagger_ab_acy(u) result(ustag)
        ! Stagger from Ab => Ac 

        implicit none 

        real(wp), intent(IN)  :: u(:,:) 
        real(wp) :: ustag(size(u,1),size(u,2)) 

        ! Local variables 
        integer :: i, j, nx, ny  

        nx = size(u,1)
        ny = size(u,2) 

        ustag = 0.0_wp 

        do j = 1, ny 
        do i = 2, nx
            ustag(i,j) = 0.5_wp*(u(i,j)+u(i-1,j))
        end do 
        end do 

        return

    end function stagger_ab_acy
    
end module staggering
