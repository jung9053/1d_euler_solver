      !----------------------------------------------------------------
      !================================================================
      !----------------------------------------------------------------

      subroutine euler_fwd_1d(q,res,delx,delt)
      use param
      use grid_dimen
      implicit none
      real(8) :: q(nvar,nx-1)
      real(8) :: res(nvar,nx-1)
      real(8) :: delt(nx-1),delx(nx-1)
      !local
      integer :: i,j

      do i=1,nx-1
      do j=1,nvar
         q(j,i) = q(j,i) - delt(i)/delx(i)*res(j,i)
      end do
      end do

      return
      end subroutine euler_fwd_1d

