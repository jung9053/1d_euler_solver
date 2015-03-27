      !----------------------------------------------------------------
      !================================================================
      !----------------------------------------------------------------

      subroutine calc_res_1d(q,res)
      use param
      use grid_dimen
      implicit none
      real(8) :: q(nvar,nx-1)
      real(8) :: res(nvar,nx-1)
      !..local array
      real(8),allocatable,dimension(:,:) :: qf1,qf2
      !local

      allocate(qf1(nvar,nx),qf2(nvar,nx))

      !---------------------\
      ! calculate face value >
      !---------------------/
      call qface_1d(qf1,qf2,q)

      !----------------------------------\
      ! calculate inviscid numerical flux >
      !----------------------------------/
      call calc_inv_flux_1d(qf1,qf2,res)

      deallocate(qf1,qf2)

      return
      end subroutine calc_res_1d

