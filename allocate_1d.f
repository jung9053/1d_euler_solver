      !----------------------------------------------------------------
      !================================================================
      !----------------------------------------------------------------

      subroutine alloc_sol_array_1d
      use param
      use grid_dimen
      use sol_array
      implicit none

      allocate(q(nvar,nx-1))

      return
      end subroutine alloc_sol_array_1d

      !----------------------------------------------------------------
      !================================================================
      !----------------------------------------------------------------

      subroutine alloc_grid_array_1d
      use grid_dimen
      use grid_array
      implicit none

      allocate(znd(nx))
      allocate(zcl(nx-1))
      allocate(delx(nx-1))

      return
      end subroutine alloc_grid_array_1d

      !----------------------------------------------------------------
      !================================================================
      !----------------------------------------------------------------

      subroutine dealloc_array_1d
      use grid_array
      implicit none

      call dealloc_grid_array_1d
      call dealloc_sol_array_1d

      !----------------------------------------------------------------
      contains
      !----------------------------------------------------------------

      !----------------------------------------------------------------

      subroutine dealloc_grid_array_1d
      use grid_array
      implicit none

      deallocate(znd)
      deallocate(zcl)
      deallocate(delx)

      return
      end subroutine dealloc_grid_array_1d

      !----------------------------------------------------------------

      subroutine dealloc_sol_array_1d
      use sol_array
      implicit none

      deallocate(q)

      return
      end subroutine dealloc_sol_array_1d

      !----------------------------------------------------------------


      end subroutine dealloc_array_1d

