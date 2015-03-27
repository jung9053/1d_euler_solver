      !----------------------------------------------------------------
      !================================================================
      !----------------------------------------------------------------

      subroutine calc_inv_flux_1d(qf1,qf2,res)
      use param
      use grid_dimen
      implicit none
      real(8) :: res(nvar,nx-1)
      real(8) :: qf1(nvar,nx),qf2(nvar,nx)
      !..local array
      real(8),allocatable,dimension(:) :: ff
      !local
      integer :: i

      !-> initialize
      res(1:nvar,1:nx-1) = 0.d0
      
      allocate(ff(nvar))

      !----------------\
      !-> boundary face >
      !----------------/
      !-> leftmost boundary
      i = 1
      call calc_roe_flux_1d(qf1(1,i),qf2(1,i),ff)
      res(1:nvar,i) = res(1:nvar,i) - ff(1:nvar)
      !-> rightmost boundary
      i = nx
      call calc_roe_flux_1d(qf1(1,i),qf2(1,i),ff)
      res(1:nvar,i-1) = res(1:nvar,i-1) + ff(1:nvar)

      !----------------\
      !-> interior face >
      !----------------/
      do i=2,nx-1 
         call calc_roe_flux_1d(qf1(1,i),qf2(1,i),ff)
         res(1:nvar,i-1) = res(1:nvar,i-1) + ff(1:nvar)
         res(1:nvar,i  ) = res(1:nvar,i  ) - ff(1:nvar)
      end do

      deallocate(ff)

      return
      end subroutine calc_inv_flux_1d

