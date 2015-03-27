      !----------------------------------------------------------------
      !================================================================
      !----------------------------------------------------------------

      subroutine wrt_out_1d(q,zcl)
      use param
      use grid_dimen
      implicit none
      real(8) :: q(nvar,nx-1)
      real(8) :: zcl(nx-1)
      !local
      integer :: i

      open(11,file='output/1d_sol.dat')
      write(11,'(a)') 'variables=x,r,u,p'
      write(11,'(a,i10)') 'zone i=',nx-1
      do i=1,nx-1
         write(11,'(100e25.10)') zcl(i),q(1:nvar,i)
      end do
      close(11)

      return
      end subroutine wrt_out_1d
