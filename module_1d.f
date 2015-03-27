
      module param
      integer,parameter :: ndim = 1
      integer,parameter :: nvar = 3
      end module param

      !----------------------------------------------------------------

      module sol_param
      integer :: niter,iorder
      real(8) :: cfl
      end module sol_param

      !----------------------------------------------------------------

      module init_cond
      real(8) :: xdia,t_out
      real(8) :: rl,ul,pl,rr,ur,pr
      end module init_cond

      !----------------------------------------------------------------

      module flow_property
      real(8) :: gam,gm1,gm2,gm3,gp1,gm1g,gp1g,ggm1,gm1i,gami
      end module flow_property

      !----------------------------------------------------------------

      module grid_dimen
      integer :: nx
      real(8) :: x_min,x_max
      end module grid_dimen

      !----------------------------------------------------------------

      module grid_array
      real(8),allocatable,dimension(:) :: znd,zcl
      real(8),allocatable,dimension(:) :: delx
      end module grid_array

      module sol_array
      real(8),allocatable,dimension(:,:) :: q
      end module sol_array


