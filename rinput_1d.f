      !----------------------------------------------------------------
      !================================================================
      !----------------------------------------------------------------

      subroutine rinput_1d
      use sol_param
      use init_cond
      use grid_dimen
      use flow_property
      implicit none

      open(1,file="euler.inp")
      read(1,*) 
      read(1,*) niter,iorder
      read(1,*) 
      read(1,*) cfl
      read(1,*) 
      read(1,*) nx
      read(1,*) 
      read(1,*) xdia,t_out
      read(1,*) 
      read(1,*) rl,ul,pl
      read(1,*) 
      read(1,*) rr,ur,pr
      close(1)

      !-> gamma
      gam = 1.4d0
      gm1  = gam - 1.d0
      gm2  = gam - 2.d0
      gm3  = gam - 3.d0
      gp1  = gam + 1.d0
      gm1g = gm1/gam
      gp1g = gp1/gam
      ggm1 = gam/gm1
      gm1i = 1.d0/gm1
      gami = 1.d0/gam

      return
      end subroutine rinput_1d

