      !----------------------------------------------------------------
      !================================================================
      !----------------------------------------------------------------

      subroutine uctime_1d(q,delx,delt,rtime,dtmin)
      use param
      use grid_dimen
      use flow_property
      use sol_param
      implicit none
      real(8) :: dtmin,rtime
      real(8) :: q(nvar,nx-1),delx(nx-1)
      real(8) :: delt(nx-1)
      !local
      integer :: i
      real(8) :: uu,aa,sprad,dt

      dtmin = 1.d10
      do i=1,nx-1
         uu = q(2,i)
         aa = dsqrt(gam*q(3,i)/q(1,i))

         sprad = dabs(uu)+aa

         dt = cfl * delx(i) / sprad
         if(dt.lt.dtmin) dtmin = dt
      end do

      do i=1,nx-1
         delt(i) = dtmin
      end do

      rtime = rtime + dtmin

      return
      end subroutine uctime_1d
