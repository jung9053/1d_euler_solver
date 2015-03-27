

      subroutine euler_1d
      use sol_array
      use grid_array
      implicit none

      !----------------\
      ! read input file >
      !----------------/
      call rinput_1d

      !--------------\
      ! set dimension >
      !--------------/
      call set_dimen_1d

      !--------------------\
      ! allocate grid array >
      !--------------------/
      call alloc_grid_array_1d

      !--------------\
      ! Generate grid >
      !--------------/
      call make_grid_1d(znd)

      !-------------------\
      ! obtain grid metric >
      !-------------------/
      call grid_metric_1d(znd,delx,zcl)

      !-------------\
      ! initialize q >
      !-------------/
      call alloc_sol_array_1d
      call initialize_1d(zcl,q)

      !-----------------\
      ! main iteration q >
      !-----------------/
      call main_iter_1d(q,delx)

      !------------\
      ! write out q >
      !------------/
      call wrt_out_1d(q,zcl)

      !----------------------\
      ! deallocate all arrays >
      !----------------------/
      call dealloc_array_1d

      return
      end subroutine euler_1d

