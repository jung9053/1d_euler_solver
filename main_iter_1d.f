      !----------------------------------------------------------------
      !================================================================
      !----------------------------------------------------------------

      subroutine main_iter_1d(q,delx)
      use param
      use grid_dimen
      use sol_param
      use init_cond
      implicit none
      real(8) :: q(nvar,nx-1),delx(nx-1)
      !..local array
      real(8),allocatable,dimension(:) :: delt
      real(8),allocatable,dimension(:,:) :: res
      !local
      integer :: iter
      real(8) :: rtime,dtmin

      allocate(res(nvar,nx-1))
      allocate(delt(nx-1))

      rtime = 0.d0
      !----------------------------------------------------------------
      do 10 iter=1,niter
      !----------------------------------------------------------------

      !-------------------------\
      ! calculate time step size >
      !-------------------------/
      call uctime_1d(q,delx,delt,rtime,dtmin)

      !-------------------\
      ! calculate residual >
      !-------------------/
      call calc_res_1d(q,res)

      !-------------------------------------\
      !-> convert primitive to conservative  >
      !-------------------------------------/
      call conpq_1d(nx-1,q)

      !-----------------\
      ! time integration >
      !-----------------/
      call euler_fwd_1d(q,res,delx,delt)

      !--------------------------------------\
      !-> convert conservative to primitive   >
      !--------------------------------------/
      call conqp_1d(nx-1,q)

      write(*,'(a,i10,2e25.10)') 'iter,dtmin,rtime=',iter,dtmin,rtime

      if(rtime.gt.t_out) exit

      !----------------------------------------------------------------
   10 continue
      !----------------------------------------------------------------


      deallocate(res,delt)

      return
      end subroutine main_iter_1d

