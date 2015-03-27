      !----------------------------------------------------------------
      !================================================================
      !----------------------------------------------------------------

      subroutine grid_metric_1d(znd,delx,zcl)
      use grid_dimen
      implicit none
      real(8) :: znd(nx)
      real(8) :: delx(nx-1),zcl(nx-1)
      !local
      integer :: i
      real(8) :: dx

      dx = (x_max - x_min)/dble(nx-1)
      do i=1,nx-1
         delx(i) = dx
         zcl(i) = znd(i) + 0.5d0*dx
      end do


      return
      end subroutine grid_metric_1d
