      module gravity_grid_construction
      ! Pre-evaluate gravity on cell centers and at interfaces
      
      use global_parameters
      use grav_func
      
      contains
      
      subroutine set_gravity_grid
      
      ! Evaluate gravity at cell centers
      Gphi_c = (/ ( phi(r(j)), j = 1-Ng,N+Ng) /)
      
      ! Evaluate gravity at cell interfaces
      Gphi_i = (/ ( phi(r_edg(j)), j = 1-Ng,N+Ng) /)
      
           
      end subroutine set_gravity_grid
      
      ! End of module
      end module gravity_grid_construction
