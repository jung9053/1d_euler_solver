      !----------------------------------------------------------------
      !================================================================
      !----------------------------------------------------------------

      subroutine initialize_1d(zcl,q)
      use param
      use init_cond
      use grid_dimen
      implicit none
      real(8) :: zcl(nx-1),q(nvar,nx-1)
      !local
      integer :: i

      do i=1,nx-1
         if(zcl(i).lt.xdia) then
            q(1,i) = rl
            q(2,i) = ul
            q(3,i) = pl
         else
            q(1,i) = rr
            q(2,i) = ur
            q(3,i) = pr
         endif
      end do

!      do i=1,nx-1
!         if(zcl(i).lt.xdia) then
!            q(1,i) = 27./7.
!            q(2,i) = 4.*sqrt(35.)/9.
!            q(3,i) = 31./3.
!         else
!            q(1,i) = 1+0.2*dsin(5.*zcl(i))
!            q(2,i) = 0.
!            q(3,i) = 1.
!         endif
!      end do





      return
      end subroutine initialize_1d
