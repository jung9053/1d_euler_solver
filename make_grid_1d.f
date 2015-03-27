      !----------------------------------------------------------------
      !================================================================
      !----------------------------------------------------------------

      subroutine make_grid_1d(znd)
      use grid_dimen
      implicit none
      real(8) :: znd(nx)
      !local
      integer :: i

      x_min = 0.d0
      x_max = 1.d0


      do i=1,nx
         znd(i) = x_min + (x_max - x_min)/dble(nx-1)*dble(i-1)
      end do
      return
      end subroutine make_grid_1d
