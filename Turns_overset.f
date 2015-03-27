c********1*********2*********3*********4*********5*********6*********7**
      program umturns2d
c                                                                      c
c  the code was cleaned-up primarily by j. sitaraman to f90 type 
c  constructs and  options for oversetting and other good stuff to work
c
c  last modified 09/05/2005 by karthik.
c
c  j.d baeder did previously claim to have cleaned it up in 1991 :)    c
c  to structure, vectorize and streamline to make all options working. c
c                                                                      c
c  many new limiters etc. added in 1996, along with better metrics     c
c                                                                      c
c**********************************************************************c
c
c
c     note: b.c. elements must be configured for each new grid topology
c           (currently c-grid)
c     note: mdim must be .ge. jdim and kdim
c
c         tape1:    grid (input)
c         tape3:    restart file (input)
c         tape5:    input file
c         tape6:    output file of summary data
c         tape7:    output file of run norms
c         tape8:    restart file (output - suitable for plot3d)
c         tape9:    grid (output)
c         tape11:   output file of cl,cd,cmp
c         tape17:   file of alpha(t) (input)
c 
c*********************************************************************** 
c 
c   input variable to airfoil in namelist inputs:
c   
c     iread  = tells whether to read in a restart file from unit 3
c            = 1 => restart
c            = 0 => initial run with no restart
c     jmax   = points in wrap around direction
c     kmax   = points in normal direction
c     jtail1 = j location of beginning airfoil (at tail on underside)
c     half   = symmetrical airfoil?
c            = 0 no symmetry assumed, 
c                        jtail2 = jmax-jtail1+1 & jle =(jmax+1)/2
c            = 1 symmetry assumed, jtail2 = jmax-1 & jle = jmax-1
c   
c   
c     npnorm = output residual to fort.7 every npnorm iterations
c              output force to fort.11 every npnorm its.
c              output spanwise forces to fort.12 every npnorm its. (rotor only)
c     nrest  = write out restart to fort.8 every nrest iterations
c     nsteps = total number of steps at end of this run
c              (typically 200-2000 larger than last run)
c   
c   
c     fsmach = free stream mach number (0 for hover)
c     alfa   = angle of attack for freestream (collective set by grid)
c     rey    = reynolds number
c     invisc = euler or navier-stokes
c            = .false. => navier-stokes
c            = .true. => euler
c     lamin  = is the flow laminar
c            = .true. => laminar
c            = .false. => fully turbulent
c   
c   
c     iunst  = type of unsteady flow
c            = 0 => steady
c            = 1 => unsteady pitching oscillation
c            = 2 => unsteady pitching ramp
c                         (need to input grid motion)
c            = 3 => prescribed pitching
c     ntac   = temporal order of accuracy
c            = 1 => first order
c            = 2 => second order
c            = 3 => third order 
c     itnmax = number of newton iterations per time step
c     dt     = time step size, should be less than 0.10 for time accuracy
c              (typically 50.0 if steady, 0.05 for time-accurate)
c     timeac = how does dt vary in space
c            = 1. => constant time step everywhere
c            otherwise space-varying dt
c     totime = total time for unsteady calculation (reset from restart file)
c   
c   
c     ilhs   = left hand side scheme
c            = 1 => LU-SGS algorithm
c            = 2 => ARC-2D with second and fourth order dissipation
c            = 3 => ARC-2D with upwind
c    iprecon = use low Mach preconditioning or not
c            = .true. => use preconditioning
c            = .false. => no preconditioning
c     Mp    = low Mach preconditioning factor
c     epse   = dissipation for implicit side (usually 0.01)
c     irhsy  = spatial order of accuracy
c            = -1 => 1st order
c            = -2 => 2nd order
c            = -3 => 3rd order
c     ilim   = type of limiting
c            =  0 => no limiting at all (muscl scheme)
c            <  0 => no limiting in k-direction (muscl scheme)
c            >  0 => limiting in both directions
c abs(ilim)  =  1 => differentiable limiter for 3rd order (koren's)
c            =  2 => differentiable limiter for 2nd order (van albada)
c            =  3 => chakravarthy-osher minmod
c            =  4 => sonic-a scheme of hyunh et al.
c            =  5 => sonic extension of chakravarthy-osher minmod
c            =  7 => cubic interpolation with no limiting
c            =  8 => quartic interpolation with no limiting
c            =  9 => quadratic reconstruction with no limiting (pade' scheme)
c    
c
c     jint   = calculate solution on only every jint points in j-direction (1)
c     kint   = calculate solution on only every kint points in k-direction (1)
c    
c     ipert  = use perturbation technique (0)
c     initpe = initialize perturbation technique (1)
c     xvor   = beginning x-location of vortex   (-5.0 chords)
c     yvor   = beginning y-location of vortex   (-.26 chords)
c     rcore  = radius of vortex core (0.05 chords)
c     vorgam = strength of vortex (.20)
c
c     rf     = reduced frequency or time for unsteady motion
c     angmax = maximum change in angle of attack
c 
c************end prologue ********************************************
      use params_global
      use ihc
c*********************************************************************
      implicit none
c*********************************************************************
      ! allocatable arrays

      real, allocatable :: s(:),q(:),qtn(:),qtnm1(:),qnewt(:)
      real, allocatable :: x(:),y(:),tscale(:),bt(:)
      integer, allocatable :: iblank(:)
      real, allocatable :: xx(:),xy(:),yx(:),yy(:)
      real, allocatable :: ug(:),vg(:),turmu(:)
      real, allocatable :: yx0(:),yy0(:),yt0(:)
      real, allocatable :: wgrid(:),wvec(:,:)
      real, allocatable :: pwork(:,:),awork(:),bwork(:),cwork(:)
      real, allocatable :: xbig(:),ybig(:)
      real, allocatable :: xole(:),yole(:),xold(:),yold(:)
      real, allocatable :: vnut(:),vnut0(:)
      integer, allocatable :: jgmx(:),kgmx(:)
      integer, allocatable :: ipointer(:,:)
C asitav
      real, allocatable :: xt2(:),yt2(:)
      real, allocatable :: xg(:),yg(:) !base grid

		real, allocatable :: tau_global(:),utau_global(:)
		integer           :: tau_dim

      ! arrays for overset meshing

      integer              :: Nj,Nk,j,k
      integer,allocatable  :: ndonor(:),nfringe(:),iisptr(:),iieptr(:)
      integer, allocatable :: imesh(:,:,:),idonor(:,:,:)
      integer, allocatable :: ibc(:,:)
      real, allocatable    :: frac(:,:,:)
      real, allocatable    :: xgl(:,:,:,:)
      integer, allocatable :: ibgl(:,:,:)

      ! local scalar variables
      
      integer ig,ig1,igq,igs,igb,igd,idsize,igl,im1,igql
      integer nmesh,jd,kd,im
      integer nrc2,logw,ii,mstop,ihar

      real resmax,resrho,rsum,resold
      real x0,y0,cl,cd,cm,cl2d,cd2d,cm2d
      real clvv2d,cdvv2d,cmvv2d
C asitav
      real cl_tot,cd_tot,cm_tot,chord !,xmin,xmin_sl
      real xle,yle,xte,yte,dtheta,psi_rot
      real, allocatable, dimension(:) :: cla,cda,cma,chorda,xlea,ylea,xtea,ytea

c** first executable statement

      write(6,*) ' '
      write(6,*) ' welcome to maryland overset turns-2d '
      write(6,*) '  this research code should help you with all your',
     <            ' airfoil problems :). '
      write(6,*) '  overset version : 04/13/2005 (jaina)'
      write(6,*)

      nmesh=6
      allocate(jgmx(nmesh),kgmx(nmesh)) ! allocate mesh pointers with 
                                        ! dummy number of meshes to start
c..  initialize data and read inputs

      call read_inputs(jgmx,kgmx,nmesh)

c*************** memory allocation for multiple mesh pointers ***********

      allocate(ipointer(nmesh,5))
      call determine_size(jgmx,kgmx,nq,nv,ipointer,nmesh,
     &     igrd,iqdim,isdim,iadim,mdim)

c***************** memory allocation block for connectivity variables**********

      allocate(nfringe(nmesh),ndonor(nmesh))
      allocate(iisptr(nmesh),iieptr(nmesh))
      allocate(imesh(igrd,2,nmesh),idonor(igrd,2,nmesh),frac(igrd,2,nmesh),
     >         ibc(igrd,nmesh))

c***************** memory allocation block for flow variables**********

      allocate(s(isdim),q(iqdim),qtn(isdim),qtnm1(isdim),qnewt(isdim))
      allocate(x(igrd),y(igrd),tscale(igrd),bt(igrd))
C asitav
      allocate(xt2(igrd),yt2(igrd))
      allocate(xg(igrd),yg(igrd))
C
      allocate(iblank(igrd))
      allocate(xx(igrd),xy(igrd),yx(igrd),yy(igrd))
      allocate(ug(igrd),vg(igrd),turmu(igrd))
      allocate(yx0(iadim),yy0(iadim),yt0(iadim))
      allocate(wgrid(igrd),wvec(mdim,25))
      allocate(pwork(mdim,3),awork(mdim),bwork(mdim),cwork(mdim))
      allocate(xbig(isdim),ybig(isdim),xole(isdim),yole(isdim))
      allocate(xold(isdim),yold(isdim))
      allocate(vnut(igrd),vnut0(igrd))

c***************** memory allocation block for airloads variables**********
C asitav (airload arrays)
      allocate(cla(nmesh),cda(nmesh),cma(nmesh),chorda(nmesh),
     >         xlea(nmesh),ylea(nmesh),xtea(nmesh),ytea(nmesh))

c********* end memory allocation block *******************************

      do im=1,nmesh
         call set_pointers_globals(im,ipointer,ig,ig1,
     &        igq,igs,igb,jd,kd,jgmx,kgmx,nmesh)
         call initia(q(igq),s(igs),x(ig),y(ig),xg(ig),yg(ig),
     &        xt2(ig),yt2(ig), 
     &        xx(ig),xy(ig),yx(ig),
     &        yy(ig),ug(ig),vg(ig),turmu(ig),vnut(ig),
     &        yx0(ig1),yy0(ig1),yt0(ig1),
     &        xbig(igb),ybig(igb),xold(igb),yold(igb),xole(igb),
     &        yole(igb),iblank(ig),im,jd,kd,cl)
         if (ipert .eq. -3) call pert3(q(igq),x(ig),y(ig),ug(ig),vg(ig),
     &        yt0(ig1),yx(ig),yy(ig),jd,kd)
      enddo

C..   read in eddy viscosity
      if (iread.gt.0) then
        if(iturb.eq.1) then
           do im=1,nmesh
              call set_pointers_globals(im,ipointer,ig,ig1,
     &             igq,igs,igb,jd,kd,jgmx,kgmx,nmesh)
              
              call restr_vnu(vnut(ig),jd,kd,3)
           enddo
        elseif (iturb.eq.2) then
           do im=1,nmesh
              call set_pointers_globals(im,ipointer,ig,ig1,
     &             igq,igs,igb,jd,kd,jgmx,kgmx,nmesh)
              
              call restr_komega(q(igq),jd,kd,3)
           enddo
        endif
        read(3) istep0
      endif

      if(iunst.gt.0) then
       theta_col=theta0-angmax*cos(rf*totime
     >                             +phase_pitch*pi/180.)

       if(iunst.gt.0 .and. islat.eq.1) then
        theta_slat = theta_slat0
        dtheta = 0.
        psi_rot = rf_slat*totime
        do ihar=1,nharmSlat
         dtheta=dtheta+ampSlat(ihar)*cos(ihar*psi_rot-phiSlat(ihar))
        enddo
        theta_slat=theta_slat+dtheta!*180./pi
       end if
      end if

         do im=1,nmesh
            call set_pointers_globals(im,ipointer,ig,ig1,
     &           igq,igs,igb,jd,kd,jgmx,kgmx,nmesh)
         call movie(q(igq),x(ig),y(ig),iblank(ig),ug(ig),
     <   vg(ig),jd,kd,1,19+100*im,18+100*im)
         enddo

      call astore(x,y,iblank,q,vnut,iturb,itrans,jgmx,kgmx,nmesh,ipointer,
     &     igrd,iqdim,jd,kd,istep0)

      if(num_grids.gt.1) then

c...perform donor search, fringe point evaluation  and hole-cutting

        Nj=jgmx(1); Nk=kgmx(1)
        do im=2,nmesh
          if(Nj.le.jgmx(im)) Nj=jgmx(im)
          if(Nk.le.kgmx(im)) Nk=kgmx(im)
        end do
        allocate(xgl(2,Nj,Nk,nmesh))
        allocate(ibgl(Nj,Nk,nmesh))
         
c.....collect grids to xgl,ibgl
        do im=1,nmesh
          call set_pointers_globals(im,ipointer,ig,ig1,
     &         igq,igs,igb,jd,kd,jgmx,kgmx,nmesh)
          do k=1,kgmx(im)
            do j=1,jgmx(im)
              xgl(1,j,k,im)=x(ig-1 + jgmx(im)*(k-1) + j)
              xgl(2,j,k,im)=y(ig-1 + jgmx(im)*(k-1) + j)
              ibgl(j,k,im) =iblank(ig-1 + jgmx(im)*(k-1) + j)
            enddo
          enddo
        enddo

c......call connectivity now

        call do_connectihc(xgl,ibgl,jgmx,kgmx,imesh,idonor,frac,
     &       ibc,ndonor,nfringe,iisptr,iieptr,igrd,Nj,Nk,nmesh)
        
c......connectivity info to all meshes
        do im=1,nmesh
          call set_pointers_globals(im,ipointer,ig,ig1,
     &         igq,igs,igb,jd,kd,jgmx,kgmx,nmesh)
          do k=1,kgmx(im)
            do j=1,jgmx(im)
               iblank(ig-1 + jgmx(im)*(k-1) + j)=ibgl(j,k,im)
            enddo
          enddo
        enddo 
         
      endif

		!### modification for generalized law of wall ###
		print*,'****** Allocating variables for wall shear stress *******'
		tau_dim = max(maxval(jgmx,1),maxval(kgmx,1))
		allocate(tau_global(tau_dim),utau_global(tau_dim))
		do j = 1,tau_dim
			tau_global(j) = 1.e-16
			utau_global(j) = 1.e-16
		enddo
		print*,'****** Completed allocating variables for wall shear stress *******'
		!### end modification for generalized law of wall ###
      
c..now perform the flow solution for each iteration or time step 

      nrc2 = nsteps - istep0
      logw=20

!	totime=0.0
      do 10 istep = 1,nrc2

c..update time and move the grid if unsteady problem

        istep0 = istep0 + 1
        if(iunst.ne.0) totime = totime + dt
        
        do ii=1,nmesh
           
           call set_pointers_globals(ii,ipointer,ig,ig1,
     &          igq,igs,igb,jd,kd,jgmx,kgmx,nmesh)
           
           if (iunst.gt.0) then
            if(ii.eq.1.and.islat.eq.1) !only the slat moves
!     >       call rigid_slat_twodof(xg(ig),yg(ig),jd,kd)
     >       call rigid_slat_twodof_nh(xg(ig),yg(ig),jd,kd)

            if(iteflap.eq.1) then
             if( nmesh.eq.1 .or.  (nmesh.eq.2 .and. (.not.bg) .and. ii.eq.2)
     >          .or. (nmesh.eq.2 .and. bg .and. ii.eq.1) 
     >          .or. (nmesh.eq.3 .and. ii.eq.2))
     >       call rigid_flap(xg(ig),yg(ig),jd,kd,.false.)
            end if

            !pure pitching about (c/4)
            if( (ii.lt.nmesh) .or. (ii.eq.nmesh.and.(.not.bg))) then
             call pitch(x(ig),y(ig),xg(ig),yg(ig),xt2(ig),yt2(ig),
     >             xx(ig),xy(ig),ug(ig),
     >             yx(ig),yy(ig),vg(ig),
     >             yx0(ig1),yy0(ig1),yt0(ig1),jd,kd,ii,.false.)
            end if
!              theta_col=alfa + angmax*sin(rf*istep*dt)
!asitav
!	      if (ii.ne.nmesh) then ! do for the blade mesh only
!                call move_new(0,x(ig),y(ig),xx(ig),xy(ig),ug(ig),
!     &             yx(ig),yy(ig),vg(ig),yx0(ig1),yy0(ig1),
!     &             yt0(ig1),jd,kd)
!	      endif
              call metfv(q(igq),x(ig),y(ig),xx(ig),xy(ig),yx(ig),yy(ig),
     &             xbig(igb),ybig(igb),xold(igb),yold(igb),
     &             xole(igb),yole(igb),jd,kd)
 7         endif

           if (itnmax.gt.1) then
             call stqol(q(igq),qtn(igs),qtnm1(igs),qnewt(igs),vnut(ig),
     &                  vnut0(ig),jd,kd)
           endif

           if (ipert.eq.-3) call pert3(q(igq),x(ig),y(ig),ug(ig),vg(ig),
     &          yt0(ig1),yx(ig1),yy(ig1),jd,kd)

           if( mod(istep,npnorm).eq.0)
     <                      write(6,101) istep,theta_col,totime
 101       format(/,' istep,angle,time =',x,i5,2(x,f10.5))

        enddo

c..    recalculate connectivity


        if(num_grids.gt.1.and.iunst.gt.0) then

c...perform donor search, fringe point evaluation  and hole-cutting

c.....collect grids to xgl,ibgl
          do im=1,nmesh
            call set_pointers_globals(im,ipointer,ig,ig1,
     &           igq,igs,igb,jd,kd,jgmx,kgmx,nmesh)
            do k=1,kgmx(im)
              do j=1,jgmx(im)
                xgl(1,j,k,im)=x(ig-1 + jgmx(im)*(k-1) + j)
                xgl(2,j,k,im)=y(ig-1 + jgmx(im)*(k-1) + j)
                ibgl(j,k,im) =iblank(ig-1 + jgmx(im)*(k-1) + j)
              enddo
            enddo
          enddo

c......call connectivity now

          call do_connectihc(xgl,ibgl,jgmx,kgmx,imesh,idonor,frac,
     &         ibc,ndonor,nfringe,iisptr,iieptr,igrd,Nj,Nk,nmesh)
        
c......connectivity info to all meshes
          do im=1,nmesh
            call set_pointers_globals(im,ipointer,ig,ig1,
     &           igq,igs,igb,jd,kd,jgmx,kgmx,nmesh)
            do k=1,kgmx(im)
              do j=1,jgmx(im)
                 iblank(ig-1 + jgmx(im)*(k-1) + j)=ibgl(j,k,im)
              enddo
            enddo
          enddo 
         
        endif

c..perform the step depending on number of newton iterations
c..stop if negative speed of sound

        do 999 itn = 1,itnmax
           
           do im=1,nmesh
              call set_pointers_globals(im,ipointer,ig,ig1,
     &        igq,igs,igb,jd,kd,jgmx,kgmx,nmesh)
           
              if (.not. iprecon) Mp = 1.0

              dualtime = dtpseudo(itn)
              call time_step(q(igq),xx(ig),xy(ig),yx(ig),yy(ig),
     <               ug(ig),vg(ig),tscale(ig),bt(ig),iblank(ig),jd,kd)

              call step(q(igq),qtn(igs),qtnm1(igs),qnewt(igs),s(igs),! q (at time steps)
     &             x(ig),y(ig),iblank(ig),                           ! grid
     &             xx(ig),xy(ig),yx(ig),yy(ig),ug(ig),vg(ig),        ! metrics
     &             yx0(ig1),yy0(ig1),yt0(ig1),                       ! surface metrics
     &             xbig(igb),ybig(igb),                              !
     &             xold(igb),yold(igb),xole(igb),yole(igb),          ! fine mesh 
     &             vnut(ig),vnut0(ig),turmu(ig),                     ! eddy viscosity
     &             tscale(ig),bt(ig),im,jd,kd,resmax,resrho,rsum,cl, ! timestep, dimensions
     &				 tau_dim,tau_global,utau_global)							! wall stress,fric vel	  

              call monitor(mstop,q(igq),jd,kd,im)
              
              if(mstop .gt. 0 ) go to 19
              if(itn.eq.1) resold=resmax
           enddo
           
           if (nmesh.gt.1) then
             call do_interpolations(q,vnut,jgmx,kgmx,ibc,imesh,idonor,frac,
     >              nfringe,ndonor,iisptr,iieptr,igrd,iqdim,nmesh)
           endif

  999   continue
  901   continue


c..   write restart files

        if (mod(istep0,nrest).eq.0.or.istep0.eq.2) then
	  		  print*,'******** WRITING RESTART FILES AT ITER = ',istep
          call astore(x,y,iblank,q,vnut,iturb,itrans,jgmx,kgmx,nmesh,
     <       ipointer,igrd,iqdim,jd,kd,istep0)

			 endif
	
c...call force and moment routine

        do im=1,nmesh
           x0 = x_cm
           y0 = 0.00
           call set_pointers_globals(im,ipointer,ig,ig1,
     &        igq,igs,igb,jd,kd,jgmx,kgmx,nmesh)

           if (bodyflag(im) .or. flatplate) then
              call force2d(logw,jd,kd,x0,y0,x(ig),y(ig),q(igq),
     <             xx(ig),xy(ig),yx(ig),yy(ig),
     <             cl,cd,cm,cl2d,cd2d,cm2d,clvv2d,cdvv2d,cmvv2d,
     <             chord,xle,yle,xte,yte,im,turmu,tau_dim,tau_global,utau_global)

C..asitav     store airload and airfoil data of each element
              cla(im)=cl; cda(im)=cd; cma(im)=cm;
              chorda(im)=chord; 
              xlea(im)=xle; ylea(im)=yle;
              xtea(im)=xte; ytea(im)=yte;
C..

              if( mod(istep,npnorm).eq.0 ) 
     >             write(10+im,610)theta_col,cl,cd,cm,float(istep0),totime
           endif
          
	  if(mod(istep,nmovie).eq.0) then
          if (istep.lt.isin) goto 10
          write(1819,*) istep,theta_col
          call movie(q(igq),x(ig),y(ig),iblank(ig),
     <		ug(ig),vg(ig),jd,kd,0,19+100*im,18+100*im)
	  endif
       enddo

       call get_totairload(nmesh,cla,cda,cma,
     >        chorda,xlea,ylea,xtea,ytea,x0,y0,
     >        cl_tot,cd_tot,cm_tot,alfa,bg,
     >        onlyslat,theta_col,mainqc,mainch)

       if( mod(istep,npnorm).eq.0 )  then
         write(20,610)theta_col,cl_tot,cd_tot,cm_tot,float(istep0),totime
       end if
          
610             FORMAT (6(x,E14.6))
611             FORMAT (7(x,E14.6))
   10 continue
   19 continue

c..finished iterations or time-steps
c..write restart files

		deallocate(tau_global,utau_global)
      write(6,*) 'finished iterations, storing solution files now..'
      call astore(x,y,iblank,q,vnut,iturb,itrans,jgmx,kgmx,nmesh,ipointer,
     &     igrd,iqdim,jd,kd,istep0)

      stop 'umturns2d completed successfully'
      end


c**********************************************************************
      subroutine astore(x,y,iblank,q,vnut,iturb,itrans,jgmx,kgmx,nmesh,
     &     ipointer,igrd,iqdim,jd,kd,istep0)
c
c     make calls to store the solution and eddy viscosity 
c     in multizoned plot3d file format
c**********************************************************************      
      implicit none
c**********************************************************************      

      integer nmesh,igrd,iqdim,jd,kd,im
      real x(igrd),y(igrd),q(iqdim),vnut(igrd)
      integer iturb,itrans
      integer iblank(igrd)
      integer jgmx(nmesh),kgmx(nmesh)
      integer ipointer(nmesh,5)
      integer istep0
      logical lamin
      
      ! local variables

      integer logg,logq,logtur
      integer i,ig,ig1,igq,igs,igb

c***  first executable statement

      logg=9  ! grid     file => fort.9
      logq=8  ! solution file => fort.8
      logtur=10 ! turbulence file => SA.dat or komega.dat 

      rewind(logg)
      rewind(logq)

      if(nmesh.gt.1) write(logg) nmesh
      write(logg) (jgmx(i),kgmx(i),i=1,nmesh)
      if(nmesh.gt.1) write(logq) nmesh
      write(logq) (jgmx(i),kgmx(i),i=1,nmesh)

      ! store solution

      do im=1,nmesh
         call set_pointers_globals(im,ipointer,ig,ig1,
     &        igq,igs,igb,jd,kd,jgmx,kgmx,nmesh)
         call store(q(igq),x(ig),y(ig),iblank(ig),jd,kd,logq,logg)
      enddo

      ! store turbulent quantities

      if(iturb.eq.1) then

        open(unit=logtur,file='SA.dat',status='unknown',
     c                                               form='formatted')
        if (itrans.eq.0) then
          write(logtur,*) "TITLE = ""Turbulence Model solution File"""
          write(logtur,*) "VARIABLES = ""X"" ""Y"" ""IBlank"" ""Vnu_t"""
        else
          write(logtur,*) "TITLE = ""Turbulence and Transition Model
     & solution File"""
          write(logtur,*) "VARIABLES = ""X"" ""Y"" ""IBlank"" ""Vnu_t""
     & ""Intermittency"" ""Re_theta"""
        endif

        do im=1,nmesh
          call set_pointers_globals(im,ipointer,ig,ig1,
     &         igq,igs,igb,jd,kd,jgmx,kgmx,nmesh)
          call storevnu(x(ig),y(ig),iblank(ig),q(igq),vnut(ig),jd,kd,logq,logtur)
        enddo

      elseif(iturb.eq.2) then

        open(unit=logtur,file='komega.dat',status='unknown',
     c                                               form='formatted')
        if (itrans.eq.0) then
          write(logtur,*) "TITLE = ""Turbulence Model solution File"""
          write(logtur,*) "VARIABLES = ""X"" ""Y"" ""IBlank"" ""Vnu_t"" ""K""
     & ""Omega""" 
        else
          write(logtur,*) "TITLE = ""Turbulence and Transition Model
     & solution File"""
          write(logtur,*) "VARIABLES = ""X"" ""Y"" ""IBlank"" ""Vnu_t"" ""K""
     & ""Omega"" ""Intermittency"" ""Re_theta"""
        endif

        do im=1,nmesh
          call set_pointers_globals(im,ipointer,ig,ig1,
     &         igq,igs,igb,jd,kd,jgmx,kgmx,nmesh)
          call storekomega(x(ig),y(ig),iblank(ig),q(igq),vnut(ig),jd,kd,logq,logtur)

        enddo
      endif

      write(logq) istep0

      close(logg)
      close(logq)
      close(logtur)

      return
      end

c***********************************************************************
      subroutine bc(q,x,y,xx,xy,yx,yy,ug,vg,yx0,yy0,yt0,im,jd,kd,bt)
c
c  explicitly update the mesh boundaries
c
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************

      integer im,jd,kd,j,k,js,je,ks,ke,idir,ib
      real q(jd,kd,nq),x(jd,kd),y(jd,kd)
      real xx(jd,kd),xy(jd,kd),yx(jd,kd),yy(jd,kd),bt(jd,kd)
      real ug(jd,kd),vg(jd,kd),yx0(jd),yy0(jd),yt0(jd)

c***  first executable statement

      do ib=1,nbc_all(im)
        js = jbcs_all(ib,im)
        je = jbce_all(ib,im)
        ks = kbcs_all(ib,im)
        ke = kbce_all(ib,im)
        if(js.lt.0) js = jmax+js+1
        if(ks.lt.0) ks = kmax+ks+1
        if(je.lt.0) je = jmax+je+1
        if(ke.lt.0) ke = kmax+ke+1
        idir = ibdir_all(ib,im)

c.. outflow bc - crude implementation (shivaji)
        if (ibtyp_all(ib,im).eq.53) then
          call bcoutflow(q,js,je,ks,ke,idir)

c.. inviscid wall bc (k = 1 only) 
        elseif (ibtyp_all(ib,im).eq.3) then
          call bcwall_generic(q,xx,xy,yx,yy,ug,vg,yx0,yy0,yt0,
     &         jd,kd,js,je,ks,ke,idir,.true.)

c.. viscous wall bc
        elseif (ibtyp_all(ib,im).eq.4) then
          call bcwall_generic(q,xx,xy,yx,yy,ug,vg,yx0,yy0,yt0,
     &         jd,kd,js,je,ks,ke,idir,.false.)

Casitav
Ca..wall bc at k=1,k=-1 
        else if (ibtyp_all(ib,im).eq.40) then
          call bcwall_wind(q,xx,xy,yx,yy,ug,vg,yx0,yy0,yt0,
     &         jd,kd,js,je,ks,ke,idir)

c.. wall bc (specifically for k=1 of the airfoil mesh)
        elseif (ibtyp_all(ib,im).eq.5) then
          call bcwall_universal(q,xx,xy,yx,yy,ug,vg,yx0,yy0,yt0,
     &         jd,kd,js,je,ks,ke,idir)
c.. symmetric bc
        elseif (ibtyp_all(ib,im).eq.11) then
          call bcsym(q,js,je,ks,ke,idir)

c.. periodic bc
        elseif (ibtyp_all(ib,im).eq.22) then
          call bcperiodic(q,js,je,ks,ke,idir)

c.. nishan: periodic bc for channel with pressure gradient to drive the mean flow
        elseif (ibtyp_all(ib,im).eq.23) then
          call bcperiodic_channel(q,js,je,ks,ke,idir)

c.. nishan: inflow bc for channel which reads inflow profile from a datafile
        elseif (ibtyp_all(ib,im).eq.28) then
          call bcinflow_channel(q,js,je,ks,ke,idir)

c.. nishan: sponge outflow bc for channel
        elseif (ibtyp_all(ib,im).eq.32) then
          call bcoutflow_channel(q,js,je,ks,ke,idir)

c.. nishan: this boundary condition that rescales the inlet boundary condition from rescaling station based on LWS rescaling
        elseif (ibtyp_all(ib,im).eq.36) then
          call bcinflow_rescale(q,x,y,js,je,ks,ke,idir)

c.. averaging bc for wake
        elseif (ibtyp_all(ib,im).eq.51) then
          call bcwake(q,x,y,jd,kd,js,je,ks,ke,idir)

c.. averaging bc for wake of O-grid
        elseif (ibtyp_all(ib,im).eq.52) then
          call bcwake_ogrid(q,x,y,jd,kd,js,je,ks,ke,idir)

c.. freesream bc
        elseif (ibtyp_all(ib,im).eq.47) then
          if (.not.iprecon) then
            call bcout(q,x,y,xx,xy,yx,yy,ug,vg,jd,kd,js,je,ks,ke,idir)
          else
             call bcout_trkl(q,x,y,xx,xy,yx,yy,ug,vg,bt,
     &                       jd,kd,js,je,ks,ke,idir)
          endif
        endif
      enddo

      return
      end

c***********************************************************************
      subroutine bcoutflow_channel(q,js,je,ks,ke,idir)

c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************

      real q(jmax,kmax,nq)
      integer js,je,ks,ke,idir

c.. local variables

      integer j,k,jc,jj,jj1,kc,kk,kk1
      integer iadd
      real scale
      
      iadd = sign(1,idir)
     
      print*,'bcoutflow running'
      if(idir.eq.1) then

			print*,'bc outflow not implemented for bcdir = ',idir

      elseif(idir.eq.-1) then
        do j = jmax-2,jmax
        do k = 1,kmax
		     scale = q(j-1,k,nq)/q(j,k,nq)
             q(j,k,1) = q(j-1,k,1)*scale
             q(j,k,2) = q(j-1,k,2)*scale
             q(j,k,3) = q(j-1,k,3)*scale
             q(j,k,4) = q(j-1,k,4)*scale
        enddo
        enddo

      elseif(idir.eq.2) then

			print*,'bc outflow not implemented for bcdir = ',idir

      elseif(idir.eq.-2) then

			print*,'bc outflow not implemented for bcdir = ',idir

      endif

      return
      end

c***********************************************************************
c***********************************************************************
      subroutine bcoutflow(q,js,je,ks,ke,idir)

c.. shivaji
c..outflow bc - just extrapolation based on mach number (poor implementation)
c
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************

      real q(jmax,kmax,nq)
      integer js,je,ks,ke,idir

c.. local variables

      integer j,k,jc,jj,jj1,kc,kk,kk1
      integer iadd,scale
      
      iadd = sign(1,idir)

      if(idir.eq.1) then

			print*,'bc outflow not implemented for bcdir = ',idir

      elseif(idir.eq.-1) then
        do j = js,je
           jc = js - 1
           do k = ks,ke
				 scale = q(jc,k,nq)/q(j,k,nq)
             q(j,k,1) = q(jc,k,1)*scale
             q(j,k,2) = q(jc,k,2)*scale
             q(j,k,3) = q(jc,k,3)*scale
             q(j,k,4) = q(jc,k,4)*scale
           enddo
        enddo

      elseif(idir.eq.2) then

			print*,'bc outflow not implemented for bcdir = ',idir

      elseif(idir.eq.-2) then

			print*,'bc outflow not implemented for bcdir = ',idir

      endif

      return
      end

c***********************************************************************
      subroutine bcperiodic(q,js,je,ks,ke,idir)
c
c..periodic boundary for overlapping planes
c
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************

      real q(jmax,kmax,nq)
      integer js,je,ks,ke,idir

c.. local variables

      integer j,k,jc,jj,jj1,kc,kk,kk1
      integer iadd
      
      iadd = sign(1,idir)

      pi = 4.0*atan(1.0)
      
      if(idir.eq.1) then
        jj  = je - js + 1
        do j = js,je
           jj1 = j - js
           jc = jmax - 2*jj + jj1
 
           do k = ks,ke
             q(j,k,1) = q(jc,k,1)
             q(j,k,2) = q(jc,k,2)
             q(j,k,3) = q(jc,k,3)
             q(j,k,4) = q(jc,k,4)
           enddo
        enddo

      elseif(idir.eq.-1) then
        jj  = je - js + 1
        do j = js,je
           jj1 = je - j
           jc = 1 + 2*jj - jj1
 
           do k = ks,ke
             q(j,k,1) = q(jc,k,1)
             q(j,k,2) = q(jc,k,2)
             q(j,k,3) = q(jc,k,3)
             q(j,k,4) = q(jc,k,4)
           enddo
        enddo

      elseif(idir.eq.2) then
        kk  = ke - ks + 1
        do k = ks,ke
           kk1 = k - ks
           kc = kmax - 2*kk + kk1
 
           do j = js,je
             q(j,k,1) = q(j,kc,1)
             q(j,k,2) = q(j,kc,2)
             q(j,k,3) = q(j,kc,3)
             q(j,k,4) = q(j,kc,4)
           enddo
        enddo

      elseif(idir.eq.-2) then
        kk  = ke - ks + 1
        do k = ks,ke
           kk1 = ke - k
           kc = 1 + 2*kk - kk1 
 
           do j = js,je
             q(j,k,1) = q(j,kc,1)
             q(j,k,2) = q(j,kc,2)
             q(j,k,3) = q(j,kc,3)
             q(j,k,4) = q(j,kc,4)
           enddo
        enddo

      endif

      return
      end

c***********************************************************************
c***********************************************************************
      subroutine bcperiodic_channel(q,js,je,ks,ke,idir)
c
c..periodic boundary for overlapping planes
c
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************

      real q(jmax,kmax,nq)
      integer js,je,ks,ke,idir

c.. local variables

      integer j,k,jc,jj,jj1,kc,kk,kk1
      integer iadd
c.. local variables for channel flow
      real rho1,u1,v1,p1,e1,difac
      real rho2,u2,v2,p2,e2,umag
      
      iadd = sign(1,idir)

      pi = 4.0*atan(1.0)
      
      if(idir.eq.1) then
        jj  = je - js + 1
        do j = js,je
           jj1 = j - js
           jc = jmax - 2*jj + jj1
 
           do k = ks,ke
             q(j,k,1) = q(jc,k,1)
             q(j,k,2) = q(jc,k,2)
             q(j,k,3) = q(jc,k,3)
             q(j,k,4) = q(jc,k,4)
           enddo
        enddo


      elseif(idir.eq.-1) then
        jj  = je - js + 1
        do j = js,je
           jj1 = je - j
           jc = 1 + 2*jj - jj1
           do k = ks,ke
             q(j,k,1) = q(jc,k,1)
             q(j,k,2) = q(jc,k,2)
             q(j,k,3) = q(jc,k,3)
             q(j,k,4) = q(jc,k,4)
           enddo
        enddo

c..    -----------------------------------------------------
c..    Full Stanford channel and Half Stanford channel
       umag = 0.1
       difac = 0.975

       p1 = 1.0/gamma
       rho1 = 1.0
       e1 = p1/(gamma-1) + 0.5*rho1*umag*umag
     
       p2   = p1*difac
       rho2 = rho1*difac
!      rho2 = rho1*(   (p2/p1)**(1/gamma)  )  
       e2 = p2/(gamma-1) + 0.5*rho2*umag*umag  
       
       do k = 1,kmax
	       q(1,k,4)    = e1/q(1,k,nq)
           q(1,k,1)    = rho1/q(1,k,nq)
	       q(jmax,k,4) = e2/q(jmax,k,nq)
           q(jmax,k,1) = rho2/q(jmax,k,nq)
   	   enddo

c..    -------------------------------------------------------
       if (mod(istep0,npnorm).eq.0) print*,q(2,42,2)*q(2,42,nq),q(jmax-1,42,2)*q(jmax-1,42,nq)

      elseif(idir.eq.2) then
		print*,'periodic channel not implemented in direction',idir

      elseif(idir.eq.-2) then
		print*,'periodic channel not implemented in direction',idir
      endif

      return
      end


c***********************************************************************
      subroutine bcinflow_channel(q,js,je,ks,ke,idir)
c
c.. Inflow condition is valid only when js = je
c
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************

      real q(jmax,kmax,nq)
      integer js,je,ks,ke,idir

c.. local variables

      integer j,k,jc,jj,jj1,kc,kk,kk1
      integer iadd
c.. local variables for inflow
      real umean(kmax),vmean(kmax),umag,rho1,u1,v1,e1,p1,dum1(kmax,1),dum2(kmax,1)
      real scale
      
c.. First executable statement
      iadd = sign(1,idir)

      pi = 4.0*atan(1.0)
      
      print*,'bcinflow is running...'    

      if(idir.eq.1) then

         open (unit=35, form='formatted', file='umean_turb.dat')
   	     do k=1,kmax
 		    read(35,"(E24.16)") umean(k)
 		 enddo
 	     close(35)

         umag = maxval(abs(umean))
         p1 = 1.0/gamma
         rho1 = 1.0
         vmean = 0.0

         do j = 1,1
            call myrand(dum1,kmax,1)
            call myrand(dum2,kmax,1)
         do k = 1,kmax
            scale = 1.0/q(j,k,nq)   
            q(j,k,1) = rho1*scale
            u1 = umean(k) + dum1(k,1)*0.1*umean(k)
            v1 = vmean(k) + dum2(k,1)*0.1*vmean(k)
            if(k.eq.1) then
               u1 = 0.0
               v1 = 0.0         
            endif
            q(j,k,2) = rho1*u1*scale
            q(j,k,3) = rho1*v1*scale 
            q(j,k,4) = ( p1/(gamma - 1.0) + 0.5*(u1*u1 + v1*v1) )*scale
         enddo
         enddo

      elseif(idir.eq.-1) then
		print*,'inflow channel not implemented in direction',idir

      elseif(idir.eq.2) then
		print*,'inflow channel not implemented in direction',idir

      elseif(idir.eq.-2) then
		print*,'inflow channel not implemented in direction',idir


      endif

      return
      end


c***********************************************************************
c***********************************************************************
      subroutine bcinflow_rescale(q,x,y,js,je,ks,ke,idir)
c
c
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************

      real q(jmax,kmax,nq)
      integer js,je,ks,ke,idir
      real x(jmax,kmax),y(jmax,kmax)

c.. local variables
      integer j,k,jc,jj,jj1,kc,kk,kk1
      integer iadd
c.. local variables for rescaling
      real umean_inlt(kmax),umean_recy(kmax),vmean_inlt(kmax),vmean_recy(kmax)
      real ufluc_inlt(kmax),ufluc_recy(kmax),vfluc_inlt(kmax),vfluc_recy(kmax)
      real umean_inlt_inner(kmax),umean_inlt_outer(kmax),vmean_inlt_inner(kmax),vmean_inlt_outer(kmax)
      real ufluc_inlt_inner(kmax),ufluc_inlt_outer(kmax),vfluc_inlt_inner(kmax),vfluc_inlt_outer(kmax) 
      real yplus_inlt(kmax),yplus_recy(kmax),eta_inlt(kmax),eta_recy(kmax)
      real utau_inlt,utau_recy,utau_ratio,delta_inlt,delta_recy
      real uimax,urmax,yplusmax,etamax,vnu0,ufreestream
      real rho1,u1,v1,p1,scale
      real dumy,dumu,dumv,aa,bb,ww      
      integer jrecy
c.. First executable statement
      iadd = sign(1,idir)
      pi = 4.0*atan(1.0)      
      print*,'rescaling now...'    
      jrecy = 85
      if(idir.eq.1) then
c......................................................................               		
c.. Reading the mean components from data file and computing fluctuating component
c..        print*,'Reading the mean components from data file and computing fluctuating component'
        vmean_inlt = 0.0
        open (unit=35, form='formatted', file='umean_turb.dat')
        open (unit=36, form='formatted', file='umean_recy.dat')
        open (unit=37, form='formatted', file='vmean_recy.dat')
   	    do k=1,kmax
 		   read(35,"(E24.16)") umean_inlt(k)
 		   read(36,"(E24.16)") umean_recy(k)
 		   read(37,"(E24.16)") vmean_recy(k)
           dumu = q(1,k,2)/q(1,k,1)
           dumv = q(1,k,3)/q(1,k,1)
 		   ufluc_inlt(k) = dumu - umean_inlt(k)
 		   vfluc_inlt(k) = dumv - vmean_inlt(k)
           dumu = q(jrecy,k,2)/q(jrecy,k,1)
           dumv = q(jrecy,k,3)/q(jrecy,k,1)
 		   ufluc_recy(k) = dumu - umean_recy(k)
 		   vfluc_recy(k) = dumv - vmean_recy(k)
 		enddo
 	    close(35)
        close(36)
        close(37)
        ufreestream = umean_inlt(kmax)
        vnu0 = 1.7e-3
c.. calculating the friction velocities
c..        print*,'calculating the friction velocities'
        utau_inlt = sqrt(abs(  vnu0*(   q(1    ,2,2)/q(1    ,2,1) - q(1    ,1,2)/q(1    ,1,1)  )/( y(1,2) - y(1,1) )  ))
        utau_recy = sqrt(abs(  vnu0*(   q(jrecy,2,2)/q(jrecy,2,1) - q(jrecy,1,2)/q(jrecy,1,1)  )/( y(1,2) - y(1,1) )  ))
        utau_ratio = utau_inlt/utau_recy
c.. calculating the wall/inner coordinates
c..        print*,'calculating the wall/inner coordinates'
		do k = 1,kmax
           yplus_inlt(k) = y(1,k)*utau_inlt/vnu0
           yplus_recy(k) = y(1,k)*utau_recy/vnu0
        enddo
c.. calculating the boundary layer thickness
c..        print*,'calculating the boundary layer thickness'
        uimax = maxval(umean_inlt)
        k = 1
        do while (umean_inlt(k).lt.0.99*uimax)
            k = k+1
        enddo
        delta_inlt = y(1,k)

        urmax = maxval(umean_recy)
        k = 1
        do while (umean_recy(k).lt.0.99*urmax)
            k = k+1
        enddo
        delta_recy = y(1,k);
c.. calculating the bl/outer coordinates
c..        print*,'calculating the bl/outer coordinates'
		do k = 1,kmax
           eta_inlt(k) = y(1,k)/delta_inlt
           eta_recy(k) = y(1,k)/delta_recy
        enddo
c.. rescaling
c..        print*,'rescaling as inner and outer layer'
        yplusmax = maxval(yplus_recy)
        etamax = maxval(eta_recy)
        do k = 1,kmax            
c.. rescaling assuming entire layer is inner layer
c..        print*,'rescaling assuming entire layer is inner layer'
           dumy = yplus_inlt(k)
           if(dumy.gt.yplusmax) dumy = yplusmax

           call interp1(kmax,yplus_recy,umean_recy,dumy,dumu)
           umean_inlt_inner(k) = utau_ratio*dumu              
           call interp1(kmax,yplus_recy,vmean_recy,dumy,dumv)
           vmean_inlt_inner(k) = dumv               
           call interp1(kmax,yplus_recy,ufluc_recy,dumy,dumu)
           ufluc_inlt_inner(k) = utau_ratio*dumu              
           call interp1(kmax,yplus_recy,vfluc_recy,dumy,dumv)
           vfluc_inlt_inner(k) = utau_ratio*dumv               
c.. rescaling assuming entire layer is outer layer
c..        print*,'rescaling assuming entire layer is outer layer'
           dumy = eta_inlt(k)               
           if(dumy.gt.etamax) dumy = etamax

           call interp1(kmax,eta_recy,umean_recy,dumy,dumu)
           umean_inlt_outer(k) = utau_ratio*dumu  + (1.0-utau_ratio)*ufreestream
           call interp1(kmax,eta_recy,vmean_recy,dumy,dumv)
           vmean_inlt_outer(k) = dumv
           call interp1(kmax,eta_recy,ufluc_recy,dumy,dumu)
           ufluc_inlt_outer(k) = utau_ratio*dumu
           call interp1(kmax,eta_recy,vfluc_recy,dumy,dumv)
           vfluc_inlt_outer(k) = utau_ratio*dumv            
        enddo
c.. creating the composite profile
c..        print*,'creating the composite profile'
        aa = 0.01
        bb = 0.055
        rho1 = 1.0
        p1 = 1.0/gamma
        do k = 1,kmax
            scale = 1.0/q(1,k,nq)
            dumy = eta_inlt(k)
            dumy = aa*(dumy - bb)/( (1.0-2.0*bb)*dumy  + bb )
            ww = 0.5*(1.0 + tanh(dumy)/tanh(aa) )
            if (dumy.gt.1.0) ww = 1.0            
c..            rho1 = q(jrecy,k,1)*q(jrecy,k,nq)
            u1 = (umean_inlt_inner(k) + ufluc_inlt_inner(k) )*(1.0-ww) 
     &         + (umean_inlt_outer(k) + ufluc_inlt_outer(k) )*ww    
            v1 = (vmean_inlt_inner(k) + vfluc_inlt_inner(k) )*(1.0-ww)  
     &         + (vmean_inlt_outer(k) + vfluc_inlt_outer(k) )*ww 
            q(1,k,1) = rho1*scale
            q(1,k,2) = rho1*u1*scale
            q(1,k,3) = rho1*v1*scale 
            q(1,k,4) = ( p1/(gamma - 1.0) + 0.5*(u1*u1 + v1*v1) )*scale
        enddo

c......................................................................      

      elseif(idir.eq.-1) then
		print*,'inflow channel not implemented in direction',idir

      elseif(idir.eq.2) then
		print*,'inflow channel not implemented in direction',idir

      elseif(idir.eq.-2) then
		print*,'inflow channel not implemented in direction',idir


      endif

      return
      end


c***********************************************************************
c***********************************************************************
      subroutine bcout(q,x,y,xx,xy,yx,yy,ug,vg,jd,kd,js,je,ks,ke,idir)
c
c     Outflow boundary condition 
c     Characteristic extrapolation based on Riemann invariants
c*********************************************************************
      use params_global
c*********************************************************************
      implicit none
c*********************************************************************
      integer jd,kd,js,je,ks,ke,idir
      real  q(jd,kd,nq)
      real  xx(jd,kd), xy(jd,kd), yx(jd,kd), yy(jd,kd),
     $      ug(jd,kd),vg(jd,kd),x(jd,kd),y(jd,kd)

      ! local variables
      real gm1i,gi,foso,xa,ya,rjk,wcyl
      real uvort,vvort,uind,vind
      real rk,rkm1,rkm2,rhoext,uext,vext,eext,pext,snorm
      real yxn,yyn,uinn,uexn,r2,qn,cspe,c2,vel2,vel2ext,velcheck
      real stj,entro,u,v,press,rj,rjp1,rjp2,xxn,xyn,rjm1,rjm2
      integer k,k1,k2,j,j1,j2,iadd,iadir
      
c***  first executable statement

      gm1i = 1./gm1
      gi   = 1./gamma
      foso = 1.0

      iadd = sign(1,idir)
      iadir = abs(idir)

      if(iadir.eq.1) then

        j  = js
        j1 = j  + iadd
        j2 = j1 + iadd
c
        do k = ks,ke
          uind=uinf
          vind=vinf
c
          rj    = 1./q(j,k,nq)
          rjp1  = 1./q(j1,k,1)
          rjp2  = 1./q(j2,k,1)
c..first-order extrapolation
          rhoext= (1.+foso)*q(j1,k,1)*q(j1,k,nq)-foso*q(j2,k,1)*q(j2,k,nq)
          uext  = (1.+foso)*q(j1,k,2)*rjp1-foso*q(j2,k,2)*rjp2
          vext  = (1.+foso)*q(j1,k,3)*rjp1-foso*q(j2,k,3)*rjp2
          eext  = (1.+foso)*q(j1,k,4)*q(j1,k,nq)-foso*q(j2,k,4)*q(j2,k,nq)
          pext  = gm1*(eext - 0.5*rhoext*(uext**2+vext**2))
c..calculate metrics at boundary
          snorm = -1.*iadd/sqrt(xx(j,k)**2+xy(j,k)**2)
          xxn = xx(j,k)*snorm
          xyn = xy(j,k)*snorm
c..calculate riemann invariants
          uinn = (uind-ug(j,k))*xxn + (vind-vg(j,k))*xyn
          r1 = uinn -2.*sqrt(gamma*pinf/rinf)*gm1i
          uexn = (uext-ug(j,k))*xxn + (vext-vg(j,k))*xyn
          r2 = uexn +2.*sqrt(gamma*pext/rhoext)*gm1i
c..calculate normal velocity and speed of sound based on riemann
          qn = 0.5*(r1+r2)
          cspe = (r1-r2)*gm1*0.25
          c2 = cspe**2
c..is flow relatively subsonic or supersonic?
          vel2 = (uind-ug(j,k))**2+(vind-vg(j,k))**2
          vel2ext = (uext-ug(j,k))**2+(vext-vg(j,k))**2
          velcheck = 0.5*(vel2+vel2ext)
c..calculate contributions from interior and exterior
          if(qn.lt. 0) then
c..inflow boundary
            if(velcheck .lt. 1.0) then
c..fix three and extrapolate one (pressure)
              stj = qn - uinn
              entro = rinf**gamma/pinf
              u = uind + stj*xxn
              v = vind + stj*xyn
              q(j,k,1) = (c2*entro*gi)**gm1i
              press = c2*q(j,k,1)*gi
              q(j,k,1) = q(j,k,1)*rj
              q(j,k,2) = q(j,k,1)*u
              q(j,k,3) = q(j,k,1)*v
              q(j,k,4) = press/gm1*rj + 0.5*(q(j,k,2)**2 +
     <                     q(j,k,3)**2)/q(j,k,1)
            else
c..fix four
              q(j,k,1) = rinf*rj
              q(j,k,2) = rinf*uind*rj
              q(j,k,3) = rinf*vind*rj
              press = pinf
              q(j,k,4) = press/gm1*rj + 0.5*(q(j,k,2)**2 +
     <                     q(j,k,3)**2)/q(j,k,1)
            endif
          else
c..outflow boundary
            if(velcheck .lt. 1.0) then
c..prescribe p and extrapolate rho,u,and v
              stj = qn - uexn
              entro = rhoext**gamma/pext
              u = uext + stj*xxn
              v = vext + stj*xyn
              q(j,k,1) = (c2*entro*gi)**gm1i
              press = c2*q(j,k,1)*gi
              q(j,k,1) = q(j,k,1)*rj
              q(j,k,2) = q(j,k,1)*u
              q(j,k,3) = q(j,k,1)*v
              q(j,k,4) = press/gm1*rj + 0.5*(q(j,k,2)**2 +
     <                     q(j,k,3)**2)/q(j,k,1)
            else
c..extrapolate four
              q(j,k,1) = rhoext*rj
              q(j,k,2) = rhoext*uext*rj
              q(j,k,3) = rhoext*vext*rj
              press = pext
              q(j,k,4) = press/gm1*rj + 0.5*(q(j,k,2)**2 +
     <                     q(j,k,3)**2)/q(j,k,1)
            endif
          endif
        enddo

      elseif (iadir.eq.2) then
        k  = ks
        k1 = k  + iadd
        k2 = k1 + iadd
c
        do j = 2,jmax-1
          uind=uinf
          vind=vinf
c
          rk    = 1./q(j,k,nq)
          rkm1  = 1./q(j,k1,1)
          rkm2  = 1./q(j,k2,1)
c..first-order extrapolation
          rhoext= (1.+foso)*q(j,k1,1)*q(j,k1,nq)-foso*q(j,k2,1)*q(j,k2,nq)
          uext  = (1.+foso)*q(j,k1,2)*rkm1-foso*q(j,k2,2)*rkm2
          vext  = (1.+foso)*q(j,k1,3)*rkm1-foso*q(j,k2,3)*rkm2
          eext  = (1.+foso)*q(j,k1,4)*q(j,k1,nq)-foso*q(j,k2,4)*q(j,k2,nq)
          pext  = gm1*(eext - 0.5*rhoext*(uext**2+vext**2))
c..calculate metrics at boundary
          snorm = -1.*iadd/sqrt(yx(j,k)**2+yy(j,k)**2)
          yxn = yx(j,k)*snorm
          yyn = yy(j,k)*snorm
c..calculate riemann invariants
          uinn = (uind-ug(j,k))*yxn + (vind-vg(j,k))*yyn
          r1 = uinn -2.*sqrt(gamma*pinf/rinf)*gm1i
          uexn = (uext-ug(j,k))*yxn + (vext-vg(j,k))*yyn
          r2 = uexn +2.*sqrt(gamma*pext/rhoext)*gm1i
c..calculate normal velocity and speed of sound based on riemann
          qn = 0.5*(r1+r2)
          cspe = (r2-r1)*gm1*0.25
          c2 = cspe**2
c..is flow relatively subsonic or supersonic?
          vel2 = (uind-ug(j,k))**2+(vind-vg(j,k))**2
          vel2ext = (uext-ug(j,k))**2+(vext-vg(j,k))**2
          velcheck = 0.5*(vel2+vel2ext)
c..calculate contributions from interior and exterior
          if(qn .lt. 0) then
c..inflow boundary
            if(velcheck .lt. 1.0) then
c..fix three and extrapolate one (pressure)
              stj = qn - uinn
              entro = rinf**gamma/pinf
              u = uind + stj*yxn
              v = vind + stj*yyn
              q(j,k,1) = (c2*entro*gi)**gm1i
              press = c2*q(j,k,1)*gi
              q(j,k,1) = q(j,k,1)*rk
              q(j,k,2) = q(j,k,1)*u
              q(j,k,3) = q(j,k,1)*v
              q(j,k,4) = press/gm1*rk + 0.5*(q(j,k,2)**2 +
     <                     q(j,k,3)**2)/q(j,k,1)
            else
c..fix four
              q(j,k,1) = rinf*rk
              q(j,k,2) = rinf*uind*rk
              q(j,k,3) = rinf*vind*rk
              press = pinf
              q(j,k,4) = press/gm1*rk + 0.5*(q(j,k,2)**2 +
     <                     q(j,k,3)**2)/q(j,k,1)
            endif
          else
c..outflow boundary
            if(velcheck .lt. 1.0) then
c..prescribe p and extrapolate rho,u,and v
              stj = qn - uexn
              entro = rhoext**gamma/pext
              u = uext + stj*yxn
              v = vext + stj*yyn
              q(j,k,1) = (c2*entro*gi)**gm1i
              press = c2*q(j,k,1)*gi
              q(j,k,1) = q(j,k,1)*rk
              q(j,k,2) = q(j,k,1)*u
              q(j,k,3) = q(j,k,1)*v
              q(j,k,4) = press/gm1*rk + 0.5*(q(j,k,2)**2 +
     <                     q(j,k,3)**2)/q(j,k,1)
            else
c..extrapolate four
              q(j,k,1) = rhoext*rk
              q(j,k,2) = rhoext*uext*rk
              q(j,k,3) = rhoext*vext*rk
              press = pext
              q(j,k,4) = press/gm1*rk + 0.5*(q(j,k,2)**2 +
     <                     q(j,k,3)**2)/q(j,k,1)
            endif
          endif
        enddo

      endif

      return
      end

c***********************************************************************
      subroutine bcout_trkl(q,x,y,xx,xy,yx,yy,ug,vg,bt,jd,kd,
     <                      js,je,ks,ke,idir)
c
c     Outflow boundary condition with Turkel Preconditioning 
c     Characteristic extrapolation based on Riemann invariants
c*********************************************************************
      use params_global
c*********************************************************************
      implicit none
c*********************************************************************
      integer jd,kd,js,je,ks,ke,idir
      real  q(jd,kd,nq)
      real  xx(jd,kd), xy(jd,kd), yx(jd,kd), yy(jd,kd),
     $      ug(jd,kd),vg(jd,kd),x(jd,kd),y(jd,kd),bt(jd,kd)

      ! local variables
      real,allocatable :: qext(:),qbar(:),delq(:),delw(:)
      real smax,foso,xa,ya,rjk,wcyl
      real uvort,vvort,uind,vind
      real rhoext,uext,vext,eext,pext,snorm
      real xxn,xyn,yxn,yyn
      real rho,rrho,uu,vv,e,uv2,cjkl,c2i,ge,qq
      real Z,R,S,Zi,bSq
      integer k,k1,k2,j,j1,j2,iadd,iadir
      
      allocate(qext(4),qbar(4),delq(4),delw(4))

c***  first executable statement
      foso = 1.0

      iadd = sign(1,idir)
      iadir = abs(idir)

      if(iadir.eq.1) then

        j  = js
        j1 = j  + iadd
        j2 = j1 + iadd
c
        do k = ks,ke
          xa = x(j,k)
          ya = y(j,k)
          uind=uinf
          vind=vinf
c..calculate metrics at boundary
          snorm = -1.*iadd/sqrt(xx(j,k)**2+xy(j,k)**2)
          xxn = xx(j,k)*snorm
          xyn = xy(j,k)*snorm
c..first-order extrapolation
          qext(1) = (1.+foso)*q(j1,k,1)*q(j1,k,nq)-foso*q(j2,k,1)*q(j2,k,nq)
          qext(2) = (1.+foso)*q(j1,k,2)*q(j1,k,nq)-foso*q(j2,k,2)*q(j2,k,nq)
          qext(3) = (1.+foso)*q(j1,k,3)*q(j1,k,nq)-foso*q(j2,k,3)*q(j2,k,nq)
          qext(4) = (1.+foso)*q(j1,k,4)*q(j1,k,nq)-foso*q(j2,k,4)*q(j2,k,nq)

          qbar(1) = 0.5*(qext(1)+rinf)
          qbar(2) = 0.5*(qext(2)+rinf*uind)
          qbar(3) = 0.5*(qext(3)+rinf*vind)
          qbar(4) = 0.5*(qext(4)+einf)
       
          delq(1) = 0.5*(rinf-qext(1))
          delq(2) = 0.5*(rinf*uind-qext(2))
          delq(3) = 0.5*(rinf*vind-qext(3))
          delq(4) = 0.5*(einf-qext(4))

          rho = qbar(1)
          rrho = 1./rho
          uu = qbar(2)*rrho
          vv = qbar(3)*rrho
          e = qbar(4)*rrho
          uv2 = 0.5*(uu*uu+vv*vv)
          cjkl = sqrt(ggm1*(e-uv2))
          c2i = 1./(cjkl*cjkl)
          ge = gamma*e - gm1*uv2
          qq = (uu - ug(j,k))*xxn + (vv - vg(j,k))*xyn

          bSq = Mp**2/(bt(j,k)-Mp**2*(bt(j,k)-1))

          Z = sqrt((1.-bSq)*qq*(1.-bSq)*qq+4.*bSq*cjkl*cjkl)
          Zi = 1./Z
          R = 0.5*((1.-bSq)*qq+Z)
          S = 0.5*((1.-bSq)*qq-Z)

c..multiplying by sign(lamda_pa)*inv(X_pa)
          delw(1) = (1.-uv2*gm1*c2i)*sign(1.,qq)*delq(1)
          delw(1) = delw(1) + uu*c2i*gm1*sign(1.,qq)*delq(2)
          delw(1) = delw(1) + vv*c2i*gm1*sign(1.,qq)*delq(3)
          delw(1) = delw(1) - gm1*c2i*sign(1.,qq)*delq(4)
       
          delw(2) = (vv*xxn-uu*xyn)*sign(1.,qq)*delq(1)
          delw(2) = delw(2) + xyn*sign(1.,qq)*delq(2)
          delw(2) = delw(2) - xxn*sign(1.,qq)*delq(3)
       
          delw(3) = -Zi*(uv2*gm1*c2i*S+bSq*qq)*sign(1.,R+bSq*qq)*delq(1)
          delw(3) = delw(3) + (xxn*bSq+S*uu*gm1*c2i)*Zi*
     c                        sign(1.,R+bSq*qq)*delq(2)
          delw(3) = delw(3) + (xyn*bSq+S*vv*gm1*c2i)*Zi*
     c                        sign(1.,R+qq*bSq)*delq(3)
          delw(3) = delw(3) - S*gm1*c2i*Zi*sign(1.,R+qq*bSq)*delq(4)
       
          delw(4) = Zi*(uv2*gm1*c2i*R+qq*bSq)*sign(1.,S+bSq*qq)*delq(1)
          delw(4) = delw(4) - (xxn*bSq+R*uu*gm1*c2i)*Zi*
     c                        sign(1.,S+bSq*qq)*delq(2)
          delw(4) = delw(4) - (xyn*bSq+R*vv*gm1*c2i)*Zi*
     c                        sign(1.,S+qq*bSq)*delq(3)
          delw(4) = delw(4) + R*gm1*c2i*Zi*sign(1.,S+qq*bSq)*delq(4)
       
c.. multiplying by (X_pa)
          delq(1) = delw(1)+delw(3)+delw(4)
          delq(2) = uu*delw(1)+xyn*delw(2)+(uu+xxn*R/bSq)*delw(3)
          delq(2) = delq(2) + (uu+xxn*S/bSq)*delw(4)
          delq(3) = vv*delw(1)-xxn*delw(2)+(vv+xyn*R/bSq)*delw(3)
          delq(3) = delq(3) + (vv+xyn*S/bSq)*delw(4)
          delq(4) = uv2*delw(1)+(uu*xyn-vv*xxn)*delw(2)
          delq(4) = delq(4) + (ge+R*qq/bSq)*delw(3)
          delq(4) = delq(4) + (ge+S*qq/bSq)*delw(4)

c.. qbar + delq
          q(j,k,1) = (qbar(1)-delq(1))/q(j,k,nq)
          q(j,k,2) = (qbar(2)-delq(2))/q(j,k,nq)
          q(j,k,3) = (qbar(3)-delq(3))/q(j,k,nq)
          q(j,k,4) = (qbar(4)-delq(4))/q(j,k,nq)

        enddo
      elseif(iadir.eq.2) then

        k  = ks
        k1 = k  + iadd
        k2 = k1 + iadd
c
        do j = js,je
          xa = x(j,k)
          ya = y(j,k)
          uind=uinf
          vind=vinf

c..calculate metrics at boundary
          snorm = -1.*iadd/sqrt(yx(j,k)**2+yy(j,k)**2)
          yxn = yx(j,k)*snorm
          yyn = yy(j,k)*snorm
c..first-order extrapolation
          qext(1) = (1.+foso)*q(j,k1,1)*q(j,k1,nq)-foso*q(j,k2,1)*q(j,k2,nq)
          qext(2) = (1.+foso)*q(j,k1,2)*q(j,k1,nq)-foso*q(j,k2,2)*q(j,k2,nq)
          qext(3) = (1.+foso)*q(j,k1,3)*q(j,k1,nq)-foso*q(j,k2,3)*q(j,k2,nq)
          qext(4) = (1.+foso)*q(j,k1,4)*q(j,k1,nq)-foso*q(j,k2,4)*q(j,k2,nq)

          qbar(1) = 0.5*(qext(1)+rinf)
          qbar(2) = 0.5*(qext(2)+rinf*uind)
          qbar(3) = 0.5*(qext(3)+rinf*vind)
          qbar(4) = 0.5*(qext(4)+einf)
       
          delq(1) = 0.5*(rinf-qext(1))
          delq(2) = 0.5*(rinf*uind-qext(2))
          delq(3) = 0.5*(rinf*vind-qext(3))
          delq(4) = 0.5*(einf-qext(4))

          rho = qbar(1)
          rrho = 1./rho
          uu = qbar(2)*rrho
          vv = qbar(3)*rrho
          e = qbar(4)*rrho
          uv2 = 0.5*(uu*uu+vv*vv)
          cjkl = sqrt(ggm1*(e-uv2))
          c2i = 1./(cjkl*cjkl)
          ge = gamma*e - gm1*uv2
          qq = (uu-ug(j,k))*yxn + (vv-vg(j,k))*yyn

          bSq = Mp**2/(bt(j,k)-Mp**2*(bt(j,k)-1))

          Z = sqrt((1.-bSq)*qq*(1.-bSq)*qq+4.*bSq*cjkl*cjkl)
          Zi = 1./Z
          R = 0.5*((1.-bSq)*qq+Z)
          S = 0.5*((1.-bSq)*qq-Z)

c..multiplying by sign(lamda_pa)*inv(X_pa)
          delw(1) = (1.-uv2*gm1*c2i)*sign(1.,qq)*delq(1)
          delw(1) = delw(1) + uu*c2i*gm1*sign(1.,qq)*delq(2)
          delw(1) = delw(1) + vv*c2i*gm1*sign(1.,qq)*delq(3)
          delw(1) = delw(1) - gm1*c2i*sign(1.,qq)*delq(4)
       
          delw(2) = (vv*yxn-uu*yyn)*sign(1.,qq)*delq(1)
          delw(2) = delw(2) + yyn*sign(1.,qq)*delq(2)
          delw(2) = delw(2) - yxn*sign(1.,qq)*delq(3)
       
          delw(3) = -Zi*(uv2*gm1*c2i*S+bSq*qq)*sign(1.,R+bSq*qq)*delq(1)
          delw(3) = delw(3) + (yxn*bSq+S*uu*gm1*c2i)*Zi*
     c                        sign(1.,R+bSq*qq)*delq(2)
          delw(3) = delw(3) + (yyn*bSq+S*vv*gm1*c2i)*Zi*
     c                        sign(1.,R+qq*bSq)*delq(3)
          delw(3) = delw(3) - S*gm1*c2i*Zi*sign(1.,R+qq*bSq)*delq(4)
       
          delw(4) = Zi*(uv2*gm1*c2i*R+qq*bSq)*sign(1.,S+bSq*qq)*delq(1)
          delw(4) = delw(4) - (yxn*bSq+R*uu*gm1*c2i)*Zi*
     c                        sign(1.,S+bSq*qq)*delq(2)
          delw(4) = delw(4) - (yyn*bSq+R*vv*gm1*c2i)*Zi*
     c                        sign(1.,S+qq*bSq)*delq(3)
          delw(4) = delw(4) + R*gm1*c2i*Zi*sign(1.,S+qq*bSq)*delq(4)
       
c.. multiplying by (X_pa)
          delq(1) = delw(1)+delw(3)+delw(4)
          delq(2) = uu*delw(1)+yyn*delw(2)+(uu+yxn*R/bSq)*delw(3)
          delq(2) = delq(2) + (uu+yxn*S/bSq)*delw(4)
          delq(3) = vv*delw(1)-yxn*delw(2)+(vv+yyn*R/bSq)*delw(3)
          delq(3) = delq(3) + (vv+yyn*S/bSq)*delw(4)
          delq(4) = uv2*delw(1)+(uu*yyn-vv*yxn)*delw(2)
          delq(4) = delq(4) + (ge+R*qq/bSq)*delw(3)
          delq(4) = delq(4) + (ge+S*qq/bSq)*delw(4)

c.. qbar - delq
          q(j,k,1) = (qbar(1)-delq(1))/q(j,k,nq)
          q(j,k,2) = (qbar(2)-delq(2))/q(j,k,nq)
          q(j,k,3) = (qbar(3)-delq(3))/q(j,k,nq)
          q(j,k,4) = (qbar(4)-delq(4))/q(j,k,nq)

        enddo

      endif

      return
      end

c***********************************************************************
      subroutine bctany( q,xx,xy,yx,yy,ug,vg,jd,kd,js,je,ks,ke,idir )
c
c  inviscid : tangency b.c. at a eta equal constant surface
c  viscous  : no slip condition
c
c  presently, set up for surface at k = kreq
c             direction of freestream = idir   
c
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************
      
      integer jd,kd,ks,ke,idir
      real q(jd,kd,nq)
      real xx(jd,kd), xy(jd,kd), yx(jd,kd), yy(jd,kd)
      real ug(jd,kd), vg(jd,kd)
      integer js,je

      ! local variables
      real,allocatable :: u(:), v(:)
      real,allocatable :: scale(:)

      integer k,k1,k2,k3,ihigh,iadd,iadir
      integer j,j1,j2,j3
      real foso,t,scal,u1,u2,u3
      real v1,v2,v3

      allocate(u(mdim), v(mdim))
      allocate(scale(mdim))
      
c***  first executable statement

      iadd = sign(1,idir)
      iadir = abs(idir)

      if(iadir.eq.1) then
        j  = js
        j1 = j  + iadd
        j2 = j1 + iadd
        j3 = j2 + iadd
        ihigh = 1
        foso = 0.9999

		if( invisc ) then

          do k = ks,ke
            v1 = (q(j1,k,2)*yx(j,k)+q(j1,k,3)*yy(j,k))/q(j1,k,1)
            v2 = (q(j2,k,2)*yx(j,k)+q(j2,k,3)*yy(j,k))/q(j2,k,1)
            v3 = (q(j3,k,2)*yx(j,k)+q(j3,k,3)*yy(j,k))/q(j3,k,1)
            v(k)  = (1.+foso)*v1 - foso*v2
c...higher order bc?
            if(ihigh.eq.1) v(k)  = 3.*v1-3.*v2+v3
            u(k)  = ug(j,k)*xx(j,k)+vg(j,k)*xy(j,k)
		  enddo
		
		else
		  do k = ks,ke
            u(k) = ug(j,k)*xx(j,k) + vg(j,k)*xy(j,k)
            v(k) = ug(j,k)*yx(j,k) + vg(j,k)*yy(j,k)
		  enddo 
		  call uv(js,je,ks,ke,q,xx,xy,yx,yy,u,v,jd,kd,scal,iadir)	
		endif
		

      elseif(iadir.eq.2) then
        k  = ks
        k1 = k + iadd
        k2 = k1 + iadd
        k3 = k2 + iadd
        ihigh = 1
        foso = 0.9999

c..slowly turn on the tangency boundary condition 
c..over 30 time steps

        t = (float(istep0)-1.)/30.
        if(t.gt.1.) t=1.
        scal = (10.-15.*t+6.*t*t)*t**3
        !scal = 1.

        if( invisc ) then

c..for inviscid:
c  extrapolate contravariant velocities to the surface
c  solve for the surface cartesian velocity components

cjdb..note: really we are extrapolating the physical values!

          do j = js,je
            u1 = (q(j,k1,2)*xx(j,k)+q(j,k1,3)*xy(j,k))/q(j,k1,1)
            u2 = (q(j,k2,2)*xx(j,k)+q(j,k2,3)*xy(j,k))/q(j,k2,1)
            u3 = (q(j,k3,2)*xx(j,k)+q(j,k3,3)*xy(j,k))/q(j,k3,1)
            u(j)  = (1.+foso)*u1 - foso*u2
c...higher order bc?
            if(ihigh.eq.1) u(j)  = 3.*u1-3.*u2+u3
            v(j)  = ug(j,k)*yx(j,k)+vg(j,k)*yy(j,k)
          enddo
        else
c
c..for viscous flows the no slip condition at the surface requires
c  setting the contravariant velocities u = 0., & v = 0. for
c  no suction or injection at the surface. this satisfies tangency
c  at the walls  
c..grs
c  solve for the surface cartesian velocity components
c
          do 55 j = js, je
            u(j) = ug(j,k)*xx(j,k) + vg(j,k)*xy(j,k)
            v(j) = ug(j,k)*yx(j,k) + vg(j,k)*yy(j,k)
   55     continue
        endif
        call uv(js,je,ks,ke,q,xx,xy,yx,yy,u,v,jd,kd,scal,iadir)

      endif
c
      return
      end


c***********************************************************************
      subroutine bcwake(q,x,y,jd,kd,js,je,ks,ke,idir)
c
c  for complete cgrid   
c  treatment for in the trailing wake at k = 1
c
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************

      integer jd,kd,js,je,ks,ke,idir
      real q(jd,kd,nq)
      real x(jd,kd), y(jd,kd)

      ! local variables
      
      integer k,k1,k2,ihigh,j,jj,iadd,iadir
      real sc1,sc2,scale1,scale2,sc,qav1,qav2,qav3,qav4
      
c**   first executable statement
      iadd = sign(1,idir)
      iadir = abs(idir)

      if(iadir.eq.1) then
        print*,'idir = ',idir,' is not implemented in bcwake'
      elseif(iadir.eq.2) then

        k  = ks
        k1 = k  + iadd
        k2 = k1 + iadd

c..in the wake, average points above and below

        ihigh = 1
          do 30 j = js, je
            jj = jmax - j + 1
            sc1 = q(j,k1,nq)/q(j,k,nq)
            sc2 = q(j,k2,nq)/q(j,k,nq)
            scale1= q(jj,k1,nq)/q(j,k,nq)
            scale2= q(jj,k2,nq)/q(j,k,nq)
            sc= q(j,k,nq)/q(jj,k,nq)
            qav1 = 0.5*(q(j,k1,1)*sc1+q(jj,k1,1)*scale1)
            qav2 = 0.5*(q(j,k1,2)*sc1+q(jj,k1,2)*scale1)
            qav3 = 0.5*(q(j,k1,3)*sc1+q(jj,k1,3)*scale1)
            qav4 = 0.5*(q(j,k1,4)*sc1+q(jj,k1,4)*scale1)
c...higher order averaging?
            if(ihigh.eq.1) qav1 = (-q(j,k2,1)*sc2+4.*q(j,k1,1)*sc1
     <              +4.*q(jj,k1,1)*scale1-q(jj,k2,1)*scale2)/6.
            if(ihigh.eq.1) qav2 = (-q(j,k2,2)*sc2+4.*q(j,k1,2)*sc1
     <              +4.*q(jj,k1,2)*scale1-q(jj,k2,2)*scale2)/6.
            if(ihigh.eq.1) qav3 = (-q(j,k2,3)*sc2+4.*q(j,k1,3)*sc1
     <              +4.*q(jj,k1,3)*scale1-q(jj,k2,3)*scale2)/6.
            if(ihigh.eq.1) qav4 = (-q(j,k2,4)*sc2+4.*q(j,k1,4)*sc1
     <              +4.*q(jj,k1,4)*scale1-q(jj,k2,4)*scale2)/6.
            q(j,k,1)= qav1
            q(j,k,2)= qav2
            q(j,k,3)= qav3
            q(j,k,4)= qav4
            q(jj,k,1)= qav1*sc
            q(jj,k,2)= qav2*sc
            q(jj,k,3)= qav3*sc
            q(jj,k,4)= qav4*sc
   30   continue
      
      endif
c
      return
      end
c***********************************************************************
      subroutine bcwake_ogrid(q,x,y,jd,kd,js,je,ks,ke,idir)
c
c  for complete ogrid   
c
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************

      integer jd,kd,js,je,ks,ke,idir
      real q(jd,kd,nq)
      real x(jd,kd), y(jd,kd)

      ! local variables
      
      integer j,k,jj,j1,iadd,iadir
      real sc1,sc2,qav1,qav2,qav3,qav4
      
c**   first executable statement
      iadd = sign(1,idir)
      iadir = abs(idir)

      if(iadir.eq.1) then

        j  = js
        j1 = j + iadd
        jj = jmax-j1+1

c..in the wake, average points above and below

        do 30 k = ks, ke
          sc1 = q(j1,k,nq)/q(j,k,nq)
          sc2 = q(jj,k,nq)/q(j,k,nq)
          qav1 = 0.5*(q(j1,k,1)*sc1+q(jj,k,1)*sc2)
          qav2 = 0.5*(q(j1,k,2)*sc1+q(jj,k,2)*sc2)
          qav3 = 0.5*(q(j1,k,3)*sc1+q(jj,k,3)*sc2)
          qav4 = 0.5*(q(j1,k,4)*sc1+q(jj,k,4)*sc2)
          q(j,k,1)= qav1
          q(j,k,2)= qav2
          q(j,k,3)= qav3
          q(j,k,4)= qav4
   30   continue

      elseif (iadir.eq.2) then
        print*,'idir = ',idir,' is not implemented in bcwake_ogrid'
      endif
c
      return
      end

c***********************************************************************
      subroutine bcwall(q,xx,xy,yx,yy,ug,vg,yx0,yy0,yt0,
     &     jd,kd,js,je,ks,ke,idir)
c
c  solid wall boundary condition for a k = 1 solid surface. 
c
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************

      integer jd,kd,js,je,ks,ke,idir
      real q(jd,kd,nq)
      real xx(jd,kd),xy(jd,kd),yx(jd,kd),yy(jd,kd)
      real yx0(jd),yy0(jd),yt0(jd)
      real ug(jd,kd),vg(jd,kd)
      
      ! local variables
      real,allocatable :: p(:,:),a(:),b(:),c(:)

      integer k,k1,k2,k3,ihigh,j,iadd,iadir
      real foso,use,vse,rinv,u,v,rj,us,vs,ue,ve,t,scal,rscal,ajacinv

      allocate(p(jd,kd), a(jd), b(jd), c(jd))

c**** first executable statement
      iadd = sign(1,idir)
      iadir = abs(idir)

      if(iadir.eq.1) then
        print*,'idir = ',idir,' is not implemented in bcwall'
      elseif(iadir.eq.2) then
        k  = ks
        k1 = k  + iadd
        k2 = k1 + iadd
        k3 = k2 + iadd
        ihigh = 0
        foso = 0.9999

c..compute surface velocities

        call bctany( q,xx,xy,yx,yy,ug,vg,jd,kd,js,je,ks,ke,idir )

c..at the trailing edge (dont do if half-plane!)

        if(half.ne.1 .and. invisc) then
          use = 0.5*(q(je,k1,2)/q(je,k1,1)+q(js,k1,2)/q(js,k1,1))
          vse = 0.5*(q(je,k1,3)/q(je,k1,1)+q(js,k1,3)/q(js,k1,1))
c...higher order?
          if(ihigh.eq.1) 
     <      use = (-q(je,k2,2)/q(je,k2,1)+4.*q(je,k1,2)/q(je,k1,1)
     <           +4.*q(js,k1,2)/q(js,k1,1)-q(js,k2,2)/q(js,k2,1))/6.
          if(ihigh.eq.1) 
     <      vse = (-q(je,k2,3)/q(je,k2,1)+4.*q(je,k1,3)/q(je,k1,1)
     <           +4.*q(js,k1,3)/q(js,k1,1)-q(js,k2,3)/q(js,k2,1))/6.
          q(js,k,2) = use*q(js,k,1)
          q(je,k,2) = use*q(je,k,1)
          q(js,k,3) = vse*q(js,k,1)
          q(je,k,3) = vse*q(je,k,1)
        endif
c
        if(half.eq.1) q(js,k,3) = vg(js,k)

c..compute static wall pressure and total energy 

        if(ibcwp.ne.0) then
          call bcwp( a,b,c,p,q,xx,xy,yx,yy,yx0,yy0,yt0,ug,vg,jd,kd,
     <               js,je,ks,ke,idir )
        else
          call bcwp0(p,q,xx,xy,yx,yy,jd,kd,js,je,ks,ke,idir)
        endif

c..extrapolate density to the surface

        if( invisc ) then
          do j = js,je
            rinv = 1./q(j,k,1)
            u = q(j,k,2)*rinv
            v = q(j,k,3)*rinv
            q(j,k,1) = (1.+foso)*q(j,k1,1)*q(j,k1,nq)/q(j,k,nq)
     <                    - foso*q(j,k2,1)*q(j,k2,nq)/q(j,k,nq)
c...higher order?
            if(ihigh.eq.1)
     <        q(j,k,1) = 3.*q(j,k1,1)*q(j,k1,nq)/q(j,k,nq)
     <                 - 3.*q(j,k2,1)*q(j,k2,nq)/q(j,k,nq)
     <                     +q(j,k3,1)*q(j,k3,nq)/q(j,k,nq)
            q(j,k,2) = u*q(j,k,1)
            q(j,k,3) = v*q(j,k,1)
            q(j,k,4) = p(j,k)/(gm1*q(j,k,nq)) +0.5*(q(j,k,2)**2
     $                  +q(j,k,3)**2)/q(j,k,1)
          enddo
        else
          do 12 j = js,je
            rinv = 1./q(j,k,1)
            u = q(j,k,2)*rinv
            v = q(j,k,3)*rinv
            q(j,k,1) = q(j,k1,1)*q(j,k1,nq)/q(j,k,nq)
            q(j,k,2) = u*q(j,k,1)
            q(j,k,3) = v*q(j,k,1)
            q(j,k,4) = p(j,k)/(gm1*q(j,k,nq)) +0.5*(q(j,k,2)**2
     $                  +q(j,k,3)**2)/q(j,k,1)
 12       continue
        endif

c..at the trailing edge (dont do if half-plane!)

        if(half.ne.1) then
          rj = q(js,k,nq)/q(je,k,nq)
          rinv = 1./q(js,k,1)
          us = q(js,k,2)*rinv
          vs = q(js,k,3)*rinv
          rinv = 1./q(je,k,1)
          ue = q(je,k,2)*rinv
          ve = q(je,k,3)*rinv
          q(js,k,1) = 0.5*( q(js,k,1)+q(je,k,1)/rj )
          q(js,k,2) = us*q(js,k,1)
          q(js,k,3) = vs*q(js,k,1)
          q(js,k,4) = p(js,k)/(gm1*q(js,k,nq)) +0.5*(q(js,k,2)**2
     $                 +q(js,k,3)**2)/q(js,k,1)
          q(je,k,1) = q(js,k,1)*rj
          q(je,k,2) = ue*q(je,k,1)
          q(je,k,3) = ve*q(je,k,1)
          q(je,k,4) = p(je,k)/(gm1*q(je,k,nq)) +0.5*(q(je,k,2)**2
     $                 +q(je,k,3)**2)/q(je,k,1)
        endif

c..slowly turn on the wall boundary condition 
c  over 30 time steps
        t = float(istep0)/30.
        if(t.gt.1.) t=1.
        scal = (10.-15.*t+6.*t*t)*t**3
        rscal = 1.-scal
        do j=js,je
          ajacinv = 1./q(j,k,nq)
          q(j,k,1) = q(j,k,1)*scal+rscal*rinf*ajacinv
          q(j,k,2) = q(j,k,2)*scal+rscal*rinf*uinf*ajacinv
          q(j,k,3) = q(j,k,3)*scal+rscal*rinf*vinf*ajacinv
          q(j,k,4) = q(j,k,4)*scal+rscal*einf*ajacinv
        enddo
      
      endif

      return
      end





c***********************************************************************
      subroutine bcwall_universal(q,xx,xy,yx,yy,ug,vg,yx0,yy0,yt0,
     &     jd,kd,js,je,ks,ke,idir)
c
c  solid wall boundary condition for a k = 1 solid surface. 
c
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************

      integer jd,kd,js,je,ks,ke,idir
      real q(jd,kd,nq)
      real xx(jd,kd),xy(jd,kd),yx(jd,kd),yy(jd,kd)
      real yx0(jd),yy0(jd),yt0(jd)
      real ug(jd,kd),vg(jd,kd)
      
      ! local variables
      real,allocatable :: p(:,:),a(:),b(:),c(:)

      integer k,k1,k2,k3,ihigh,iadd,iadir
      integer j,j1,j2,j3
      real foso,use,vse,rinv,u,v,rj,us,vs,ue,ve,t,scal,rscal,ajacinv

      allocate(p(jd,kd), a(jd), b(jd), c(jd))

c**** first executable statement
      iadd = sign(1,idir)
      iadir = abs(idir)

      if(iadir.eq.1) then
        print*,'idir = ',idir,' testing in bcwall'
        j  = js
        j1 = j  + iadd
        j2 = j1 + iadd
        j3 = j2 + iadd
        ihigh = 0
        foso = 0.9999

c..compute surface velocities

        call bctany( q,xx,xy,yx,yy,ug,vg,jd,kd,js,je,ks,ke,idir )

c..at the trailing edge (dont do if half-plane!)

        if(half.ne.1 .and. invisc) then
          print*,'special condition for trailing edge not implemented for vertical wall'
c...higher order?
          if(ihigh.eq.1) print*,'higher order vertical bc not implemented'
        endif
c
        if(half.eq.1) print*,'half plane vertical bc not implemented'

c..compute static wall pressure and total energy 

        if(ibcwp.ne.0) then
          call bcwp( a,b,c,p,q,xx,xy,yx,yy,yx0,yy0,yt0,ug,vg,jd,kd,
     <               js,je,ks,ke,idir )
        else
          call bcwp0(p,q,xx,xy,yx,yy,jd,kd,js,je,ks,ke,idir)
        endif

c..extrapolate density to the surface

        if( invisc ) then
          do k = ks,ke
            rinv = 1./q(j,k,1)
            u = q(j,k,2)*rinv
            v = q(j,k,3)*rinv
            q(j,k,1) = (1.+foso)*q(j1,k,1)*q(j1,k,nq)/q(j,k,nq)
     <                    - foso*q(j2,k,1)*q(j2,k,nq)/q(j,k,nq)
c...higher order?
            if(ihigh.eq.1)
     <        q(j,k,1) = 3.*q(j1,k,1)*q(j1,k,nq)/q(j,k,nq)
     <                 - 3.*q(j2,k,1)*q(j2,k,nq)/q(j,k,nq)
     <                     +q(j3,k,1)*q(j3,k,nq)/q(j,k,nq)
            q(j,k,2) = u*q(j,k,1)
            q(j,k,3) = v*q(j,k,1)
            q(j,k,4) = p(j,k)/(gm1*q(j,k,nq)) +0.5*(q(j,k,2)**2
     $                  +q(j,k,3)**2)/q(j,k,1)
		  enddo
        else
          do 12 k = ks,ke
            rinv = 1./q(j,k,1)
            u = q(j,k,2)*rinv
            v = q(j,k,3)*rinv
            q(j,k,1) = q(j1,k,1)*q(j1,k,nq)/q(j,k,nq)
            q(j,k,2) = u*q(j,k,1)
            q(j,k,3) = v*q(j,k,1)
            q(j,k,4) = p(j,k)/(gm1*q(j,k,nq)) +0.5*(q(j,k,2)**2
     $                  +q(j,k,3)**2)/q(j,k,1)
 12       continue
        endif

c..at the trailing edge (dont do if half-plane!)

        if(half.ne.1) then
			print*,'aab trailing edge calculation not done'
        endif

c..slowly turn on the wall boundary condition 
c  over 30 time steps
        t = float(istep0)/30.
        if(t.gt.1.) t=1.
        scal = (10.-15.*t+6.*t*t)*t**3
        rscal = 1.-scal
        do k=ks,ke
          ajacinv = 1./q(j,k,nq)
          q(j,k,1) = q(j,k,1)*scal+rscal*rinf*ajacinv
          q(j,k,2) = q(j,k,2)*scal+rscal*rinf*uinf*ajacinv
          q(j,k,3) = q(j,k,3)*scal+rscal*rinf*vinf*ajacinv
          q(j,k,4) = q(j,k,4)*scal+rscal*einf*ajacinv
        enddo
     

 

      elseif(iadir.eq.2) then
        k  = ks
        k1 = k  + iadd
        k2 = k1 + iadd
        k3 = k2 + iadd
        ihigh = 0
        foso = 0.9999

c..compute surface velocities

        call bctany( q,xx,xy,yx,yy,ug,vg,jd,kd,js,je,ks,ke,idir )

c..at the trailing edge (dont do if half-plane!)

        if(half.ne.1 .and. invisc) then
          use = 0.5*(q(je,k1,2)/q(je,k1,1)+q(js,k1,2)/q(js,k1,1))
          vse = 0.5*(q(je,k1,3)/q(je,k1,1)+q(js,k1,3)/q(js,k1,1))
c...higher order?
          if(ihigh.eq.1) 
     <      use = (-q(je,k2,2)/q(je,k2,1)+4.*q(je,k1,2)/q(je,k1,1)
     <           +4.*q(js,k1,2)/q(js,k1,1)-q(js,k2,2)/q(js,k2,1))/6.
          if(ihigh.eq.1) 
     <      vse = (-q(je,k2,3)/q(je,k2,1)+4.*q(je,k1,3)/q(je,k1,1)
     <           +4.*q(js,k1,3)/q(js,k1,1)-q(js,k2,3)/q(js,k2,1))/6.
          q(js,k,2) = use*q(js,k,1)
          q(je,k,2) = use*q(je,k,1)
          q(js,k,3) = vse*q(js,k,1)
          q(je,k,3) = vse*q(je,k,1)
        endif
c
        if(half.eq.1) q(js,k,3) = vg(js,k)

c..compute static wall pressure and total energy 

        if(ibcwp.ne.0) then
          call bcwp( a,b,c,p,q,xx,xy,yx,yy,yx0,yy0,yt0,ug,vg,jd,kd,
     <               js,je,ks,ke,idir )
        else
          call bcwp0(p,q,xx,xy,yx,yy,jd,kd,js,je,ks,ke,idir)
        endif

c..extrapolate density to the surface

        if( invisc ) then
          do j = js,je
            rinv = 1./q(j,k,1)
            u = q(j,k,2)*rinv
            v = q(j,k,3)*rinv
            q(j,k,1) = (1.+foso)*q(j,k1,1)*q(j,k1,nq)/q(j,k,nq)
     <                    - foso*q(j,k2,1)*q(j,k2,nq)/q(j,k,nq)
c...higher order?
            if(ihigh.eq.1)
     <        q(j,k,1) = 3.*q(j,k1,1)*q(j,k1,nq)/q(j,k,nq)
     <                 - 3.*q(j,k2,1)*q(j,k2,nq)/q(j,k,nq)
     <                     +q(j,k3,1)*q(j,k3,nq)/q(j,k,nq)
            q(j,k,2) = u*q(j,k,1)
            q(j,k,3) = v*q(j,k,1)
            q(j,k,4) = p(j,k)/(gm1*q(j,k,nq)) +0.5*(q(j,k,2)**2
     $                  +q(j,k,3)**2)/q(j,k,1)
		  enddo
        else
          do j = js,je
            rinv = 1./q(j,k,1)
            u = q(j,k,2)*rinv
            v = q(j,k,3)*rinv
            q(j,k,1) = q(j,k1,1)*q(j,k1,nq)/q(j,k,nq)
            q(j,k,2) = u*q(j,k,1)
            q(j,k,3) = v*q(j,k,1)
            q(j,k,4) = p(j,k)/(gm1*q(j,k,nq)) +0.5*(q(j,k,2)**2
     $                  +q(j,k,3)**2)/q(j,k,1)
		  enddo
        endif

c..at the trailing edge (dont do if half-plane!)

        if(half.ne.1) then
          rj = q(js,k,nq)/q(je,k,nq)
          rinv = 1./q(js,k,1)
          us = q(js,k,2)*rinv
          vs = q(js,k,3)*rinv
          rinv = 1./q(je,k,1)
          ue = q(je,k,2)*rinv
          ve = q(je,k,3)*rinv
          q(js,k,1) = 0.5*( q(js,k,1)+q(je,k,1)/rj )
          q(js,k,2) = us*q(js,k,1)
          q(js,k,3) = vs*q(js,k,1)
          q(js,k,4) = p(js,k)/(gm1*q(js,k,nq)) +0.5*(q(js,k,2)**2
     $                 +q(js,k,3)**2)/q(js,k,1)
          q(je,k,1) = q(js,k,1)*rj
          q(je,k,2) = ue*q(je,k,1)
          q(je,k,3) = ve*q(je,k,1)
          q(je,k,4) = p(je,k)/(gm1*q(je,k,nq)) +0.5*(q(je,k,2)**2
     $                 +q(je,k,3)**2)/q(je,k,1)
        endif

c..slowly turn on the wall boundary condition 
c  over 30 time steps
        t = float(istep0)/30.
        if(t.gt.1.) t=1.
        scal = (10.-15.*t+6.*t*t)*t**3
        rscal = 1.-scal
        do j=js,je
          ajacinv = 1./q(j,k,nq)
          q(j,k,1) = q(j,k,1)*scal+rscal*rinf*ajacinv
          q(j,k,2) = q(j,k,2)*scal+rscal*rinf*uinf*ajacinv
          q(j,k,3) = q(j,k,3)*scal+rscal*rinf*vinf*ajacinv
          q(j,k,4) = q(j,k,4)*scal+rscal*einf*ajacinv
        enddo
      
      endif

      return
      end


c***********************************************************************





c***********************************************************************
      subroutine bcwall_generic(q,xx,xy,yx,yy,ug,vg,yx0,yy0,yt0,
     &     jd,kd,js,je,ks,ke,idir,invsc)
c
c  generic solid wall boundary conditions.
c
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************

      integer jd,kd,js,je,ks,ke,idir
      logical invsc
      real q(jd,kd,nq)
      real xx(jd,kd),xy(jd,kd),yx(jd,kd),yy(jd,kd)
      real yx0(jd),yy0(jd),yt0(jd)
      real ug(jd,kd),vg(jd,kd)
      
      ! local variables
      real,allocatable :: a(:),b(:),c(:)

      integer k,k1,k2,k3,ihigh,j,iadd,iadir
      real foso,rinv,u,v,rj,us,vs,ue,ve,t,scal,rscal,ajacinv
      real p1,p2,press
      logical tmp

      allocate(a(jd), b(jd), c(jd))
c**** first executable statement
      iadd = sign(1,idir)
      iadir = abs(idir)

c...setting wind tunnel walls to be inviscid
      tmp = invisc
      invisc = invsc

      if(iadir.eq.1) then
        print*,'idir = ',idir,' is not implemented in bcwall_generic'
      elseif(iadir.eq.2) then
        k  = ks
        k1 = k  + iadd
        k2 = k1 + iadd
        k3 = k2 + iadd
        foso = 0.!0.9999

c..compute surface velocities

        call bctany( q,xx,xy,yx,yy,ug,vg,jd,kd,js,je,ks,ke,idir )

        k  = ks
        k1 = k  + iadd
        k2 = k1 + iadd

c..extrapolate pressure and density to the surface

        do j = js,je
          rinv = 1./q(j,k,1)
          u = q(j,k,2)*rinv
          v = q(j,k,3)*rinv
          q(j,k,1) = (1.+foso)*q(j,k1,1)*q(j,k1,nq)/q(j,k,nq)
     <                  - foso*q(j,k2,1)*q(j,k2,nq)/q(j,k,nq)
          q(j,k,2) = u*q(j,k,1)
          q(j,k,3) = v*q(j,k,1)

          p1 = gm1*(q(j,k1,4) -0.5*(q(j,k1,2)**2+q(j,k1,3)**2)
     $              /q(j,k1,1))*q(j,k1,nq)
          p2 = gm1*(q(j,k2,4) -0.5*(q(j,k2,2)**2+q(j,k2,3)**2)
     $              /q(j,k2,1))*q(j,k2,nq)

          press = (1.+foso)*p1 - foso*p2

          q(j,k,4) = press/(gm1*q(j,k,nq)) +0.5*(q(j,k,2)**2
     $                +q(j,k,3)**2)/q(j,k,1)
		enddo
      endif

c...resetting the inviscid to original value
      invisc = tmp

      return
      end




c***********************************************************************
      subroutine bcwall_wind(q,xx,xy,yx,yy,ug,vg,yx0,yy0,yt0,
     &     jd,kd,js,je,ks,ke,idir)
c
c  solid wall boundary condition for a k = 1 solid surface. 
c
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************

      integer jd,kd,js,je,ks,ke,idir
      real q(jd,kd,nq)
      real xx(jd,kd),xy(jd,kd),yx(jd,kd),yy(jd,kd)
      real yx0(jd),yy0(jd),yt0(jd)
      real ug(jd,kd),vg(jd,kd)

      ! local variables
      real,allocatable :: p(:,:),a(:),b(:),c(:)

      integer k,k1,k2,k3,ihigh,j,iadd,iadir
      real foso,use,vse,rinv,u,v,rj,us,vs,ue,ve,t,scal,rscal,ajacinv
      real p1,p2,press

      allocate(p(jd,kd), a(jd), b(jd), c(jd))
c**** first executable statement
      iadd = sign(1,idir)
      iadir = abs(idir)

      if(iadir.eq.1) then
        print*,'idir = ',idir,' is not implemented in bcwall_wind'
      elseif(iadir.eq.2) then

       k  = ks!1
       k1 = k  + iadd !k  + 1
       k2 = k1 + iadd !k1 + 1
       k3 = k2 + iadd !k2 + 1
       ihigh = 1
       foso = 0.!0.9999
c
       do j = js,je
         q(j,k,1) = (1.+foso)*q(j,k1,1)*q(j,k1,nq)/q(j,k,nq)
     <                 - foso*q(j,k2,1)*q(j,k2,nq)/q(j,k,nq)
         q(j,k,2) = (1.+foso)*q(j,k1,2)/q(j,k1,1)*q(j,k,1)
     <                 - foso*q(j,k2,2)/q(j,k2,1)*q(j,k,1)
         q(j,k,3) = 0.0
	 p1 = gm1*(q(j,k1,4) -0.5*(q(j,k1,2)**2+q(j,k1,3)**2)
     $              /q(j,k1,1))*q(j,k1,nq)
	 p2 = gm1*(q(j,k2,4) -0.5*(q(j,k2,2)**2+q(j,k2,3)**2)
     $              /q(j,k2,1))*q(j,k2,nq)
	 press = (1.+foso)*p1-foso*p2
         q(j,k,4) = press/(gm1*q(j,k,nq)) +0.5*(q(j,k,2)**2
     $               +q(j,k,3)**2)/q(j,k,1)
       enddo   
      end if

      return
      end


c***********************************************************************
      subroutine bcwp( a,b,rhs,p,q,xx,xy,yx,yy,yx0,yy0,yt0,
     $                 ug,vg,jd,kd,js,je,ks,ke,idir )
c
c  compute the total energy by solving the normal momentum
c  equation for the pressure at the wall. the pressure and
c  and density are used to compute the total energy. the wall
c  is at k = kreq, idir = direction of freestream
c
c  note that the wall pressure is computed only in the interior
c  of the k = constant surface.
c
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************
      
      integer jd,kd,js,je,ks,ke,idir
      real p(jd,kd), q(jd,kd,nq)
      real xx(jd,kd), xy(jd,kd), yx(jd,kd), yy(jd,kd)
      real ug(jd,kd), vg(jd,kd)
      real a(jd), b(jd), rhs(jd)
      real yx0(jd),yy0(jd),yt0(jd)

      ! local variables
      
      real eps,foso,drj1,drj2,drk1,drk2
      real ru,rv,dux,duy,dvy,dvx,exa,dtinv,rzt1,rzt2,rzt3,rztt
      real thet1,exb,eyb,eya,c1,c2,dinv,beta,pav
      integer k,k1,k2,ja,jb,jam1,jbp1,kk,iadd,iadir
      integer j,j1,j2,ka,kb,kam1,kbp1,jj
      
      data eps /1.0e-20/

c***  first executable statement
      iadd = sign(1,idir)
      iadir = abs(idir)

      if(iadir.eq.1) then
        j = js
        j1 = j + iadd
        j2 = j1+ iadd
c
        ka = ks
        kb = ke
c
        kam1 = ka - 1
        kbp1 = kb + 1
c
c..compute pressure field
c
        do jj = js,js+2*iadd,iadd
	        do k = kam1,kbp1
    	      p(jj,k)  = gm1*(q(jj,k,4) -0.5*(q(jj,k,2)**2
     &                   +q(jj,k,3)**2)/q(jj,k,1))*q(jj,k,nq)
	      	enddo
		enddo
c
        foso = 0.9999
c..skip the rhs part ??? why would you want to do that ??

c..update the total energy
c     
        do k = ks,ke
          q(j,k,4) = p(j,k)/(gm1*q(j,k,nq)) +0.5*(q(j,k,2)**2
     $                 +q(j,k,3)**2)/q(j,k,1)
   		enddo

      elseif(iadir.eq.2) then
        k = ks
        k1 = k + iadd
        k2 = k1+ iadd
c
        ja = js
        jb = je
c
        jam1 = ja - 1
        jbp1 = jb + 1
c
c..compute pressure field
c
        do kk = ks,ks+2*iadd,iadd
	        do j = jam1,jbp1
    	      p(j,kk)  = gm1*(q(j,kk,4) -0.5*(q(j,kk,2)**2
     &                   +q(j,kk,3)**2)/q(j,kk,1))*q(j,kk,nq)
	      	enddo
		enddo
c
c..set up rhs in p(j,1) and compute difference coefficients a and b
c
        foso = 0.9999
        do j = ja,jb
          drj1 = 1./q(j+1,k,1)
          drj2 = 1./q(j-1,k,1)
          drk1 = 1./q(j,k1,1)
          drk2 = 1./q(j,k,1)


          ru   = (q(j,k,2)*xx(j,k) +q(j,k,3)*xy(j,k))*q(j,k,nq)
          rv   = (q(j,k,2)*yx(j,k) +q(j,k,3)*yy(j,k))*q(j,k,nq)

c
          dux  = q(j+1,k,2)*drj1 -q(j-1,k,2)*drj2
          dvx  = q(j+1,k,3)*drj1 -q(j-1,k,3)*drj2
          duy  = q(j,k1,2)*drk1 -q(j,k,2)*drk2
          dvy  = q(j,k1,3)*drk1 -q(j,k,3)*drk2

c
c..we have to account for the unsteady motion of the blade
c..the following is the first term on the right hand side of
c..eq. (11) of pulliam and stegers aiaa j. paper, feb. 1980.
c..for quasi-steady rztt = 0.
c..for unsteady flow rztt is non-zero
c
          dtinv = 1.0/dt
          rzt1 = (-ug(j,k)*yx(j,k)-vg(j,k)*yy(j,k)-yt0(j))/dt*q(j,k,1)
          rzt2 = ((yx(j,k) - yx0(j))/dt) * q(j,k,2)
          rzt3 = ((yy(j,k) - yy0(j))/dt) * q(j,k,3)
          rztt = rzt1 + rzt2 + rzt3
          rztt = rztt * q(j,k,nq)
c
          rhs(j)   = (-yx(j,k)*(ru*dux)
     $                -yy(j,k)*(ru*dvx) )*.5 + rztt


        enddo
c
        do j = ja,jb
          c1 = yx(j,k)*xx(j,k)+xy(j,k)*yy(j,k)+eps
          c2 = yx(j,k)**2+yy(j,k)**2+eps
          dinv = -1./((1.+.5*foso)*c2)
          rhs(j) = (rhs(j)+c2*(-(1.+foso)*p(j,2)+.5*foso*p(j,3)))*dinv
          a(j)   = -.5*c1*dinv
        enddo
c
c..solve for pstar
c..set boundary conditions
c
        rhs(ja)  = -a(ja)*( p(jam1,k) ) +rhs(ja)
        rhs(jb)  = a(jb)*( p(jbp1,k) ) +rhs(jb)
c
c..invert xi operator
c     
        do j = ja+1,jb
          beta     = 1.0+a(j)*a(j-1)
          rhs(j) = (rhs(j)-a(j)*rhs(j-1))/beta
          a(j)   = a(j)/beta
        enddo

        do jj  = ja+1,jb
          j = ja + jb -jj
          rhs(j)   = rhs(j) +a(j)*rhs(j+1)
        enddo

        do j = ja,jb
          p(j,k)   = rhs(j)
        enddo
c     
c..now bcs
c     
        if(half.ne.1) then
c..higher order?
c            pav = (-p(js,k2)+4.*p(js,k1)+4.*p(je,k1)-p(je,k2))/6.
            pav = (p(js,k1)+p(je,k1))*.5
            p(js,k) = pav
            p(je,k) = p(js,k)
        else
            p(js,k) = p(js,k1)
            p(je,k) = p(je,k1)
        endif
c
c..update the total energy
c     

        do j = js,je
          q(j,k,4) = p(j,k)/(gm1*q(j,k,nq)) +0.5*(q(j,k,2)**2
     $                 +q(j,k,3)**2)/q(j,k,1)
   		enddo

      endif
c
      return
      end

c***********************************************************************
      subroutine bcwp0( p,q,xx,xy,yx,yy,jd,kd,js,je,ks,ke,idir )
c
c  compute the total energy by using p1 = p2 in lieu of
c  the normal momentum equation for the pressure at the wall. 
c  the pressure and and density are used to compute the total energy. 
c  the present implementation assumes the wall is at k = kreq
c  idir = freestream direction
c
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************

      integer jd,kd,js,je,ks,ke,idir
      real p(jd,kd), q(jd,kd,nq)
      real xx(jd,kd), xy(jd,kd), yx(jd,kd), yy(jd,kd)
      
      ! local variables
      integer k,j,k1,k2,jp,iadd,iadir
      real foso,thet1,exb,eyb,dydy,u,v,uu,vv,duj,dvj,duk,dvk,rhs

c***  first executable statement

      iadd = sign(1,idir)
      iadir = abs(idir)

      if(iadir.eq.1) then
        print*,'idir = ',idir,' is not implemented in bcwp0'
      elseif(iadir.eq.2) then

c..   compute the pressure field at k = 1,3

        do k = ks,ks+2*iadd,iadd
        do j = js,je
          p(j,k) = gm1*(q(j,k,4) -0.5*(q(j,k,2)**2+q(j,k,3)**2)
     $                  /q(j,k,1))*q(j,k,nq)
		enddo
		enddo
c
        k  = ks
        k1 = k  + iadd
        k2 = k1 + iadd
c
c..extrapolate pressure to the surface
c     
        if( invisc ) then
          foso = 1.0
          do 12 j = js,je
            p(j,k) = (1.+foso)*p(j,k1)-foso*p(j,k2)
 12       continue
        else
          do 13 j = js,je
            p(j,k) = p(j,k1)
            if(ichoice.eq.0) then
            jp=j+1
            jm=j-1
          if(j.ge.j_s2.and.j.le.j_e2) then
                  thet1=asin(xy(j,k)/sqrt(xx(j,k)**2+xy(j,k)**2))
                  exb=cos(thet1+angle_b)
                  eyb=sin(thet1+angle_b)
           dydy=(yx(j,k)*yx(j,k)+yy(j,k)*yy(j,k))
          u=q(j,k,2)/q(j,k,1)
          v=q(j,k,3)/q(j,k,1)
              uu  = xx(j,k)*u + xy(j,k)*v
              vv  = yx(j,k)*u + yy(j,k)*v
              duj=0.5*(q(jp,k,2)/q(jp,k,1) - q(jm,k,2)/q(jm,k,1))
              dvj=0.5*(q(jp,k,3)/q(jp,k,1) - q(jm,k,3)/q(jm,k,1))
              duk=q(j,k1,2)/q(j,k1,1)-q(j,k,2)/q(j,k,1)
              dvk=q(j,k1,3)/q(j,k1,1)-q(j,k,3)/q(j,k,1)
              rhs  = q(j,k,1)
     &            *( uu*(yx(j,k)*duj + yy(j,k)*dvj)
     &             + vv*(yx(j,k)*duk + yy(j,k)*dvk))*q(j,k,nq)
             rhs=rhs+vnb*2*pi*f1b*cos(2*pi*f1b*(atime-t1b))*q(j,k,1)*
     &            (exb*yx(j,k)+eyb*yy(j,k))*q(j,k,nq)
           p(j,k) = p(j,k1)+ rhs/dydy

          endif
          endif

 13       continue
        endif
c     
c..now bcs
c     
        if(half.ne.1) then
          p(js,k) = 0.5*(p(js,k1)+p(je,k1))
          p(je,k) = 0.5*(p(js,k1)+p(je,k1))
        else
          p(js,k) = p(js,k1)
          p(je,k) = p(je,k1)
        endif
c
c..update the total energy
c     
        do 51 j = js,je
          q(j,k,4) = p(j,k)/(gm1*q(j,k,nq)) +0.5*(q(j,k,2)**2
     $                 +q(j,k,3)**2)/q(j,k,1)
   51   continue

      endif
c
      return
      end

c***********************************************************************
      subroutine bcsym(q,js,je,ks,ke,idir)
c
c..symmetric bc
c
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************

      real q(jmax,kmax,nq)
      integer js,je,ks,ke,idir

c.. local variables

      integer j,k,j1,k1
      integer iadd,iadir
      real scale

c**** first executable statement

      iadd = sign(1,idir)
      iadir = abs(idir)
      if (iadir.eq.1) then
         j = js
         j1 = j + iadd
         do k = ks,ke
           scale = q(j1,k,nq)/q(j,k,nq)
           q(j,k,1) = q(j1,k,1)*scale
           q(j,k,2) = 0.0
           q(j,k,3) = q(j1,k,3)*scale
           q(j,k,4) = q(j1,k,4)*scale
         enddo
      elseif (iadir.eq.2) then
         k = ks
         k1 = k + iadd
         do j = js,je
           scale = q(j,k1,nq)/q(j,k,nq)
           q(j,k,1) = q(j,k1,1)*scale
           q(j,k,2) = q(j,k1,2)*scale
           q(j,k,3) = 0.0
           q(j,k,4) = q(j,k1,4)*scale
         enddo
      endif
c
      return
      end

c***********************************************************************
      subroutine determine_size(jgmx,kgmx,nq,nv,ipointer,nmesh,
     &     igrd,iqdim,isdim,iadim,mdim)
c
c     Determines the sizes of all arrays for allocation purposes
c     Also sets pointers
c
c     igrd  -> coordinates,metrics  -> ipointer(:,1)
c     iqdim -> q-variables          -> ipointer(:,2)
c     isdim -> flux                 -> ipointer(:,3)
c     iadim -> surface metrics      -> ipointer(:,4)
c     mdim  -> work arrays          -> ipointer(:,5) 
c***********************************************************************
      implicit none
c***********************************************************************

      integer nq,nv,nmesh,igrd,iqdim,isdim,iadim,mdim
      integer jgmx(nmesh),kgmx(nmesh)
      integer ipointer(nmesh,5)

      ! local variables
      
      integer i

c**   first executable statement

      do i=1,5
         ipointer(1,i)=1
      enddo

      do i=2,nmesh
         ipointer(i,1)=ipointer(i-1,1)+jgmx(i-1)*kgmx(i-1)
         ipointer(i,2)=ipointer(i-1,2)+jgmx(i-1)*kgmx(i-1)*nq
         ipointer(i,3)=ipointer(i-1,3)+jgmx(i-1)*kgmx(i-1)*nv
         ipointer(i,4)=ipointer(i-1,4)+jgmx(i-1)
         ipointer(i,5)=ipointer(i-1,5)+(kgmx(i-1)*2+jgmx(i-1))*20
      enddo

      i=nmesh+1
      igrd =ipointer(i-1,1)+jgmx(i-1)*kgmx(i-1)
      iqdim=ipointer(i-1,2)+jgmx(i-1)*kgmx(i-1)*nq
      isdim=ipointer(i-1,3)+jgmx(i-1)*kgmx(i-1)*nv
      iadim=ipointer(i-1,4)+jgmx(i-1)

      mdim=0
      do i=1,nmesh
         if (mdim.lt.jgmx(i)) mdim=jgmx(i)
         if (mdim.lt.kgmx(i)) mdim=kgmx(i)
      enddo
      
      return
      end

c***********************************************************************
      subroutine set_pointers_globals(im,ipointer,ig,ig1,igq,
     &     igs,igb,jd,kd,jgmx,kgmx,nmesh)
c
c     Sets the pointers of the current mesh to the global size parameters
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************
      integer im,ig,ig1,igq,igs,igb,jd,kd,nmesh
      integer ipointer(nmesh,5),jgmx(nmesh),kgmx(nmesh)
      
c***  first executable statement

      ig=ipointer(im,1)
      igq=ipointer(im,2)
      igs=ipointer(im,3)
      igb=igs
      ig1=ipointer(im,4)

      jd=jgmx(im)
      kd=kgmx(im)

      jmax=jd
      kmax=kd

      jm=jmax-1
      km=kmax-1

      jtail1 = jtail(im)
C
      if (half.eq.0) then
        jle = (jmax+1)/2
        jtail2 = jmax - jtail1 + 1
      else
        jle = jmax - 1
        jtail2 = jmax -  1
      endif

		if(flatplate) then
			jtail1 = jle_flatplate
			jtail2 = jmax
		endif

      return
      end

c***********************************************************************
      subroutine eigen( q,xx,xy,yx,yy,ug,vg,jd,kd )
c
c  if dt.lt.0 compute dt from input cnbr
c
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************
      integer jd,kd
      real q(jd,kd,nq)
      real xx(jd,kd), xy(jd,kd), yx(jd,kd), yy(jd,kd)
      real ug(jd,kd), vg(jd,kd)
      
      ! local variables

      real,allocatable :: sig(:)
      real sigmax,dxmin,dtmin,dtmax,u,v,edr,a,uu,vv,sigx,sigy
      integer k,j,jsig,ksig

      allocate(sig(jd))
   
c***  first executable statement

      sigmax = -1.e35
      dxmin  = 1.0
      dtmin  = 1000.0
      dtmax  = 0.0
c
      do 13 k = 2,km
        do j = 2,jm
c
c..the physical velocities u & v and the energy are
c
          u    = q(j,k,2)/q(j,k,1)
          v    = q(j,k,3)/q(j,k,1)
          edr  = q(j,k,4)/q(j,k,1)
c
c..the speed of sound
c
          a    = sqrt( ggm1*(edr -0.5*(u**2+v**2)) )
c
c..the contravariant velocites u & v are
c
          uu   = +(u-ug(j,k))*xx(j,k)+(v-vg(j,k))*xy(j,k)
          vv   = +(u-ug(j,k))*yx(j,k)+(v-vg(j,k))*yy(j,k)
c
          sigx = abs( uu ) +a*sqrt( (xx(j,k)**2+xy(j,k)**2) )
          sigy = abs( vv ) +a*sqrt( (yx(j,k)**2+yy(j,k)**2) )
          sig(j) = amax1( sigx,sigy )
     <             *( 1.0 + 0.002*(1.-timeac)*sqrt(q(j,k,nq)))
     <             /( 1.0 + (1.-timeac)*sqrt(q(j,k,nq)) )
		enddo

c
        do 12 j = 2,jm
          if( sig(j).gt.sigmax ) then
            sigmax = sig(j)
            jsig   = j
            ksig   = k
          end if
          dtmin = amin1( dtmin,h
     <            *( 1.0 + 0.002*(1.-timeac)*sqrt(q(j,k,nq)))
     <            /( 1.0 + (1.-timeac)*sqrt(q(j,k,nq)) ) )
          dtmax = amax1( dtmax,h
     <            *( 1.0 + 0.002*(1.-timeac)*sqrt(q(j,k,nq)))
     <            /( 1.0 + (1.-timeac)*sqrt(q(j,k,nq)) ) )
   12   continue
   13 continue
c
      if( dt.lt.0.0 ) then
        dt = cnbr*dxmin/sigmax
      end if
c
      cnbr = dt*sigmax/dxmin
      h = dt
      hd = .5*dt
c
c     write(6,601)
c     write(6,602) cnbr,dt,sigmax,jsig,ksig
c     write(6,*) ' dtmin,dtmax =',dtmin,dtmax
  601 format( '0',5x,'cnbr, dt, sigmax, jsig, ksig')
  602 format( ' ',5x,3f11.5,2x,2i5,/)
c
      return
      end

c**********************************************************************
      subroutine force2d(logw,jd,kd,x0,y0,x,y,q,xx,xy,yx,yy,
     <                   cl,cd,cm,cl2d,cd2d,cm2d,clvv2d,cdvv2d,cmvv2d,
     <                   chord,xle,yle,xte,yte,im,turmu,tau_dim,tau_global,utau_global)
c
c  this subroutine computes 2-d sectional cl, cd values.
c  modified version of subroutine clcd of arc2d.  vr. 01-15-91
c
c**********************************************************************
      use params_global
c**********************************************************************
      implicit none
c**********************************************************************
      integer jd,kd,logw,im
      real q(jd,kd,nq),x(jd,kd),y(jd,kd)
		real turmu(jd,kd)
      real xx(jd,kd),xy(jd,kd)
      real yx(jd,kd),yy(jd,kd)
      real x0,y0,cl,cd,cm,cl2d,cd2d,cm2d,clvv2d,cdvv2d,cmvv2d
c asitav
      real chord,xle,yle,xte,yte

      ! local variables
      
      real,allocatable :: zz(:)
      real ca,sa,rr,rho,u,v,e,vsq,pp,pp1,cp,cn,cc,cmle,uxi,vxi
      real ueta,veta,xix,xiy,etaxxx,etax,etay
      real tauw,cnv,ccv,cmlev,cfav,cpav,uinf2,amu
      real alngth!,chord
		real rough_fac

      integer k,k2,k3,j1,j2,jtp,j,imovstart,istart,jm1,jp1
      integer numb
      integer ja,jb,iend
      
		! friction variables
		integer tau_dim
		real tau_global(tau_dim),utau_global(tau_dim)

		real u_tau,yplus,uplus


c***  first executable statement

		!print*,'from force calc = ',size(tau_global,1),tau_dim
      allocate(zz(jd))

c..initialize constants and variables

		rough_fac = 0.0
		if(irough) rough_fac = 1.0
      pi  = 4.0 * atan( 1.0 )
      ca  = cos( pi * alf / 180.0 )
      sa  = sin( pi * alf / 180.0 )

      cl2d = 0.
      cd2d = 0.
      cm2d = 0.
      clvv2d = 0.
      cdvv2d = 0.
      cmvv2d = 0.
c
c..set limits of integration
c     
      k  = 1
      k2 = k+1
      k3 = k2+1
      j1 = jtail1
      j2 = jtail2
      jtp = jtail1+1

c
c..compute cp at grid points and store in an array 
c
      do 10 j=j1,j2
        rr=1./q(j,k,1)
        rho=q(j,k,1)*q(j,k,nq)
        u=q(j,k,2)*rr
        v=q(j,k,3)*rr
        e=q(j,k,4)*q(j,k,nq)
cjdb  
        vsq=u*u+v*v
        pp=gm1*(e-rho*vsq/2.)
        pp1=pp/pinf
        cp=2.*(pp1-1.)/(gamma*fsmach**2)
        zz(j) = cp
10    continue

        imovstart=4002
	istart=1
	iend=30000
	numb=1000
c 	if(angnewle.gt.angmax+alf-6.) numb=500

c       if(istep.gt.istart.and.istep.lt.iend.and.
c     &     mod(istep-istart,numb).eq.0) then
        if(istep0.eq.nsteps.or.
     &     mod(istep-istart,numb).eq.0) then
	 write(imovstart+1,*)angnewle
         do j=jtail1,jtail2
          write(imovstart,*) j,x(j,1),-zz(j)
         enddo
       endif

		if(mod(istep0,npnorm).eq.0 .or. istep0.eq.nsteps) then
		if(im==1) then
      open(unit=10,file='cp1.dat',status='unknown')
		elseif(im==2) then
      open(unit=10,file='cp2.dat',status='unknown')
		endif

      do j=j1,j2
        write(10,*) x(j,k),-zz(j)
      enddo
      close(10)
		endif

c..compute normal force coefficient and chord directed force coeff
c..chord taken as one in all cases, modified later

      cn = 0.
      cc = 0.
      cmle = 0.
      do 11 j=jtp,j2
        jm1 = j-1
        cpav = ( zz(j) + zz(jm1))*.5
        cn = cn - cpav*( x(j,k) - x(jm1,k))
        cc = cc + cpav*( y(j,k) - y(jm1,k))
        cmle = cmle + cpav*
     &     (  ( x(j,k)+x(jm1,k))*.5*(x(j,k)-x(jm1,k)) +
     &        ( y(j,k)+y(jm1,k))*.5*(y(j,k)-y(jm1,k)) )
   11 continue

      cn = cn - half*cn
      cmle = cmle - half*cmle
      cc = cc + half*cc
      cl2d = - cc * sa + cn * ca
      cd2d =   cc * ca + cn * sa
      cm2d = cmle - y0*cc + x0*cn

c..viscous coefficent of friction calculation
c..re already has fsmach scaling

      if(.not.invisc) then
        alngth= 1.
        uinf2 = fsmach**2
c
        j = j1
        jp1 = j+1
        amu   = (1.0+rough_fac*turmu(j,k))*rinf*alngth/rey
        uxi = q(jp1,k,2)/q(jp1,k,1)-q(j,k,2)/q(j,k,1)
        vxi = q(jp1,k,3)/q(jp1,k,1)-q(j,k,3)/q(j,k,1)
        ueta= -1.5*q(j,k,2)/q(j,k,1)+2.*q(j,k2,2)/q(j,k2,1)
     *       -.5*q(j,k3,2)/q(j,k3,1)
        veta= -1.5*q(j,k,3)/q(j,k,1)+2.*q(j,k2,3)/q(j,k2,1)
     *       -.5*q(j,k3,3)/q(j,k3,1)
        xix = xx(j,k)
        xiy = xy(j,k)
        etax = yx(j,k)
        etay = yy(j,k)
        tauw= amu*((uxi*xiy+ueta*etay)-(vxi*xix+veta*etax))
        zz(j)= tauw/(.5*rinf*uinf2)

        j = j2
        jm1 = j-1
        amu   = (1.0+rough_fac*turmu(j,k))*rinf*alngth/rey
        uxi = q(j,k,2)/q(j,k,1)-q(jm1,k,2)/q(jm1,k,1)
        vxi = q(j,k,3)/q(j,k,1)-q(jm1,k,3)/q(jm1,k,1)
        ueta= -1.5*q(j,k,2)/q(j,k,1)+2.*q(j,k2,2)/q(j,k2,1)
     *       -.5*q(j,k3,2)/q(j,k3,1)
        veta= -1.5*q(j,k,3)/q(j,k,1)+2.*q(j,k2,3)/q(j,k2,1)
     *       -.5*q(j,k3,3)/q(j,k3,1)
        xix = xx(j,k)
        xiy = xy(j,k)
        etax = yx(j,k)
        etay = yy(j,k)
        tauw= amu*((uxi*xiy+ueta*etay)-(vxi*xix+veta*etax))
        zz(j)= tauw/(.5*rinf*uinf2)

c..set limits

        ja = j1+1
        jb = j2-1
        do 110 j = ja,jb
			 amu   = (1.0+rough_fac*turmu(j,k))*rinf*alngth/rey
          jp1 = j+1
          jm1 = j-1
          uxi=.5*(q(jp1,k,2)/q(jp1,k,1)-q(jm1,k,2)/q(jm1,k,1))
          vxi=.5*(q(jp1,k,3)/q(jp1,k,1)-q(jm1,k,3)/q(jm1,k,1))
          ueta= -1.5*q(j,k,2)/q(j,k,1)+2.*q(j,k2,2)/q(j,k2,1)
     *               -.5*q(j,k3,2)/q(j,k3,1)
          veta= -1.5*q(j,k,3)/q(j,k,1)+2.*q(j,k2,3)/q(j,k2,1)
     *               -.5*q(j,k3,3)/q(j,k3,1)
          xix = xx(j,k)
          xiy = xy(j,k)
          etax = yx(j,k)
          etay = yy(j,k)
          tauw= amu*((uxi*xiy+ueta*etay)-(vxi*xix+veta*etax))
          zz(j)= tauw/(.5*rinf*uinf2)
110     continue

		if(mod(istep0,npnorm).eq.0 .or. istep0.eq.nsteps) then
		if(im==1) then
      open(unit=10,file='cf1.dat',status='unknown')
		elseif(im==2) then
      open(unit=10,file='cf2.dat',status='unknown')
		endif
		
      do j=j1,j2
		  tau_global(j) = 1.e-16	! init tau_global to zero
		  rho = q(j,k,1)*q(j,k,nq)
		  tau_global(j) = zz(j)*0.5*rinf*uinf2*rey
		  utau_global(j) = sqrt(abs(tau_global(j))/rho)
        write(10,*) x(j,k),rey*fsmach*x(j,k),zz(j)
		  !print*,tau_global(j),utau_global(j)
      enddo
      close(10)
		endif
c
c..compute viscous normal and axial forces
c
        cnv = 0.
        ccv = 0.
        cmlev = 0.
        do 111 j=jtp,j2
          jm1 = j-1
          cfav = ( zz(j) + zz(jm1))*.5
          ccv = ccv + cfav*( x(j,k) - x(jm1,k))
          cnv = cnv + cfav*( y(j,k) - y(jm1,k))
          cmlev = cmlev - cfav*
     *     (  (x(j,k)+x(jm1,k))*.5*(y(j,k) -y(jm1,k)) -
     *        (y(j,k)+y(jm1,k))*.5*(x(j,k) -x(jm1,k))   )
111     continue

        cnv = cnv - half*cnv
        cmlev = cmlev - half*cmlev
        ccv = ccv + half*ccv
        clvv2d = - ccv * sa + cnv * ca
        cdvv2d =   ccv * ca + cnv * sa
        cmvv2d = cmlev - y0*ccv + x0*cnv
      endif

c..modified for chord being non-unity

      chord = 1.0 !sqrt((x(jtail1,k)-x(jle,k))**2
!     &       +     (y(jtail1,k)-y(jle,k))**2)
      cl2d=cl2d/chord
      clvv2d=clvv2d/chord
      cd2d=cd2d/chord
      cdvv2d=cdvv2d/chord
      !cm2d=cm2d/chord
      !cmvv2d=cmvv2d/chord
      cm2d=cm2d/chord**2
      cmvv2d=cmvv2d/chord**2

      cl = cl2d+clvv2d
      cd = cd2d+cdvv2d
      cm = cm2d+cmvv2d

C store le and te info
      xle=x(jle,k); yle=y(jle,k)
      xte=x(jtail1,k); yte=y(jtail1,k)

		!print*,'enter j-index'
		!read(6,*)j
		j = 241
		do k = 2,kd
			rr = 1.0/q(j,k,1)
			u = q(j,k,2)*rr
			v = q(j,k,3)*rr
			!write(251,*)y(j,k),u,v
			rho = q(j,1,1)*q(j,1,nq)
			u_tau = sqrt(0.5*rinf*uinf2*abs(zz(j))*rey/rho)
			rho = q(j,k,1)*q(j,k,nq)
			yplus = sqrt((x(j,k)-x(j,1))**2.0+(y(j,k)-y(j,1))**2.0)
			yplus = yplus*u_tau*rho*sqrt(rey)
			uplus = sqrt(rey)*sqrt(u*u+v*v)/u_tau
			!write(999,*)yplus,uplus,log(yplus)/0.41+5.1
		enddo
		!close(251)
		!close(999)

      return
      end

c*********************************************************************
      subroutine grid(x,y,iblank,jd,kd,logg)
c
c  read grid from disk
c
c*********************************************************************
      use params_global
c*********************************************************************
      implicit none
c*********************************************************************
      integer jd,kd,logg
      real x(jd,kd),y(jd,kd)
      integer iblank(jd,kd)
      
      ! local variables
      integer j,k

c***  first executable statement

	if(num_grids.gt.1) then
      
      read(logg)((x(j,k),j=1,jmax),k=1,kmax),
     <     ((y(j,k),j=1,jmax),k=1,kmax),
     <     ((iblank(j,k),j=1,jmax),k=1,kmax)

	else
      read(logg)((x(j,k),j=1,jmax),k=1,kmax),
     <     ((y(j,k),j=1,jmax),k=1,kmax)

	do k=1,kmax
		do j=1,jmax
		iblank(j,k)=1
		enddo
	enddo


      endif
c
      return
      end

c***********************************************************************
      subroutine preilu2d(q,s,jd,kd,js,je,ks,ke,xx,xy,yx,yy,ug,vg,turmu,
     >            tscale,bt)
c
c  calculate the implicit inversion of the lhs
c  this involves two bidiagonal scalar inversions
c
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************
      integer jd,kd,js,je,ks,ke
      real q(jd,kd,nq), s(jd,kd,nv),turmu(jd,kd),tscale(jd,kd),bt(jd,kd)
      real ug(jd,kd), vg(jd,kd)
      real xx(jd,kd), xy(jd,kd), yx(jd,kd), yy(jd,kd)

      ! local variables
      real,allocatable :: d(:,:,:),st(:,:,:),tmp(:)
      real,allocatable :: a(:,:),c(:,:)
      real,allocatable :: uv(:,:),vn(:,:),ge(:,:)
      real,allocatable :: qx(:,:),qy(:,:),cx(:,:),cy(:,:)
      real,allocatable :: b(:,:)
      integer,allocatable :: ms(:),me(:)
      integer j,k,m,i,k1,jj
      integer itgs,j1,n
      real eps2,epsv,dj,uu,vv,uv2,cjkl,rj1,rj2
      real qq1,qqx,rr1,rk1,rk2,qq2,qqy,rr2,vnu,svt
      real ri1,ri2,qq,cc,sp1,sm2,chkx,spec,s1,s2,bSq
      real a2,a4,chky,sv,di,s3,s4
      real scale1,scale2,sav1,sav2,sav3,sav4,scale,term

      allocate(d(jd,kd,4),st(jd,kd,4),tmp(4))
      allocate(a(mdim,4),c(mdim,4))
      allocate(uv(jd,kd),vn(jd,kd),ge(jd,kd))
      allocate(qx(jd,kd),qy(jd,kd),cx(jd,kd),cy(jd,kd))
      allocate(b(mdim,4))
      allocate(ms(mdim*2),me(mdim*2))
c***  first executable statement

      eps2 = epse*2.
      epsv = 1. + eps2

c..set-up for hyper-plane loop

      do 1 m=js+ks,je+ke
        i     = max(js,m-ke)
        ms(m) = i
        me(m) = min(je,m-js)
    1 continue

c
c..store density,u,v in nonconservative variables                  
c
      do 2 k = ks-1, ke+1
      do 2 j = js-1, je+1
        dj       = 1.0 / q(j,k,1)
        q(j,k,2) = q(j,k,2)*dj
        q(j,k,3) = q(j,k,3)*dj
        q(j,k,4) = q(j,k,4)*dj
        q(j,k,1) = q(j,k,1)*q(j,k,nq)
    2 continue
c
c..in the wake, average points above and below
c..for faster convergence?
c
      print*,'Not averaging the wake in ilu3d'
      if(half.ne.1.and.1.eq.0) then
        k  = 1
        k1 = k + 1
        do 11 j = 1, jtail1 - 1
           jj = jmax - j + 1
c
          uu  = q(j,k,2)
          vv  = q(j,k,3)
          uv2 = 0.5*( uu*uu + vv*vv )
          cjkl = sqrt( ggm1*( q(j,k,4) - uv2 ) )
          ge(j,k) = gamma*q(j,k,4) - gm1*uv2
c
          rj1 = xx(j,k)
          rj2 = xy(j,k)
          qq1 = rj1*uu + rj2*vv
          qqx = abs( -rj1*ug(j,k)-rj2*vg(j,k) + qq1 )
          rr1 = sqrt( rj1**2 + rj2**2 )
          rk1 = yx(j,k)
          rk2 = yy(j,k)
          qq2 = rk1*uu + rk2*vv
          qqy = abs( -rk1*ug(j,k)-rk2*vg(j,k) + qq2 )
          rr2 = sqrt( rk1**2 + rk2**2 )
          vnu = 2.0*rr2*rr2*(rmue+turmu(j,k))/(rey*q(j,k,1))
          svt = tscale(j,k)
     <

          bSq = Mp**2/(bt(j,k)-Mp**2*(bt(j,k)-1))

          term=gm1*(bSq-1.)/(cjkl*cjkl)
          d(j,k,1)=1.0/
     &     (1.0+svt*(1.+term*uv2)*((qqx+qqy+cjkl*(rr1+rr2))*epsv+vnu))
          d(j,k,2)=1.0/
     &     (1.0+svt*(1.-term*uu*uu)*((qqx+qqy+cjkl*(rr1+rr2))*epsv+vnu))
          d(j,k,3)=1.0/
     &     (1.0+svt*(1.-term*vv*vv)*((qqx+qqy+cjkl*(rr1+rr2))*epsv+vnu))
          d(j,k,4)=1.0/
     &     (1.0+svt*(1.+term*ge(j,k))*((qqx+qqy+cjkl*(rr1+rr2))*epsv+vnu))

            scale = q(j,k1,nq)/q(j,k,nq)
            scale1= q(jj,k1,nq)/q(j,k,nq)
            scale2= q(j,k,nq)/q(jj,k,nq)
            sav1 = 0.5*(s(j,k1,1)*scale+s(jj,k1,1)*scale1)
            sav2 = 0.5*(s(j,k1,2)*scale+s(jj,k1,2)*scale1)
            sav3 = 0.5*(s(j,k1,3)*scale+s(jj,k1,3)*scale1)
            sav4 = 0.5*(s(j,k1,4)*scale+s(jj,k1,4)*scale1)
            s(j,k,1)= sav1*d(j,k,1)*.95
            s(j,k,2)= sav2*d(j,k,2)*.95
            s(j,k,3)= sav3*d(j,k,3)*.95
            s(j,k,4)= sav4*d(j,k,4)*.95
            s(jj,k,1)= sav1*scale2*d(j,k,1)*.95
            s(jj,k,2)= sav2*scale2*d(j,k,2)*.95
            s(jj,k,3)= sav3*scale2*d(j,k,3)*.95
            s(jj,k,4)= sav4*scale2*d(j,k,4)*.95
   11   continue
      endif
c
c..forward sweep
c..setup d, the diagonal term and and store some arrays
c
      do 111 k = ks-1,ke+1
      do 111 j = js-1,je+1
          uu  = q(j,k,2)
          vv  = q(j,k,3)
          uv2 = 0.5*( uu*uu + vv*vv )
          cjkl = sqrt( ggm1*( q(j,k,4) - uv2 ) )
c
          rj1 = xx(j,k)
          rj2 = xy(j,k)
          qq1 = rj1*uu + rj2*vv
          qqx = abs( -rj1*ug(j,k)-rj2*vg(j,k) + qq1 )
          rr1 = sqrt( rj1**2 + rj2**2 )
          rk1 = yx(j,k)
          rk2 = yy(j,k)
          qq2 = rk1*uu + rk2*vv
          qqy = abs( -rk1*ug(j,k)-rk2*vg(j,k) + qq2 )
          rr2 = sqrt( rk1**2 + rk2**2 )
          vnu = 2.0*rr2*rr2*(rmue+turmu(j,k))/(rey*q(j,k,1))
          svt = tscale(j,k)
          ge(j,k) = gamma*q(j,k,4) - gm1*uv2

          bSq = Mp**2/(bt(j,k)-Mp**2*(bt(j,k)-1))

          term=gm1*(bSq-1.)/(cjkl*cjkl)
          d(j,k,1)=1.0/
     &     (1.0+svt*(1.+term*uv2)*((qqx+qqy+cjkl*(rr1+rr2))*epsv+vnu))
          d(j,k,2)=1.0/
     &     (1.0+svt*(1.-term*uu*uu)*((qqx+qqy+cjkl*(rr1+rr2))*epsv+vnu))
          d(j,k,3)=1.0/
     &     (1.0+svt*(1.-term*vv*vv)*((qqx+qqy+cjkl*(rr1+rr2))*epsv+vnu))
          d(j,k,4)=1.0/
     &     (1.0+svt*(1.+term*ge(j,k))*((qqx+qqy+cjkl*(rr1+rr2))*epsv+vnu))
c
          uv(j,k) = uv2
          qx(j,k) = qq1
          cx(j,k) = cjkl*rr1
          qy(j,k) = qq2
          cy(j,k) = cjkl*rr2
          vn(j,k) = vnu
          
  111 continue
c
      do 1000 itgs=1,1 
c
c..loop on hyper-plane
c
      do 120 m = js+ks,je+ke
c
c..setup a contribution in j-direction
c
      do 121 j = ms(m),me(m)
        k  = m-j
        j1 = j-1
c
        uu  = q(j1,k,2)
        vv  = q(j1,k,3)
        uv2 = uv(j1,k)
        ri1 = xx(j1,k)
        ri2 = xy(j1,k)
        qq  = qx(j1,k)
        cc  = cx(j1,k)
        qqx = -ri1*ug(j1,k)-ri2*vg(j1,k) + qq
        sp1 = abs(qqx) + cc
        sm2 = eps2*sp1
        chkx= 0.5 + sign( 0.5, qqx+cc )
        spec= chkx*( qqx + sp1 ) + sm2
c
        s1 = s(j1,k,1)
        s2 = s(j1,k,2)
        s3 = s(j1,k,3)
        s4 = s(j1,k,4)
c
        a4 = chkx*( ri1*s2 + ri2*s3 - qq*s1 )
        a2 = chkx*gm1*( uv2*s1 - ( uu*s2 + vv*s3) + s4 )
        a(j,1) = a4                  +spec*s1
        a(j,2) = ri1*a2 + uu*a4      +spec*s2
        a(j,3) = ri2*a2 + vv*a4      +spec*s3
        a(j,4) = qq*a2 + ge(j1,k)*a4 +spec*s4
  121 continue
c
c..setup b contribution in k-direction
c
      do 122 j = ms(m),me(m)
        k  = m-j
        k1 = k-1
c
        uu  = q(j,k1,2)
        vv  = q(j,k1,3)
        uv2 = uv(j,k1)
        ri1 = yx(j,k1)
        ri2 = yy(j,k1)
        qq  = qy(j,k1)
        cc  = cy(j,k1)
        vnu = vn(j,k1)
        qqy = -ri1*ug(j,k1)-ri2*vg(j,k1) + qq
        sp1 = abs(qqy) + cc
        sm2 = eps2*sp1
        chky= 0.5 + sign( 0.5, qqy+cc )
        spec= chky*( qqy + sp1 ) + vnu + sm2
c
        s1 = s(j,k1,1)
        s2 = s(j,k1,2)
        s3 = s(j,k1,3)
        s4 = s(j,k1,4)
c
        a4 = chky*( ri1*s2 + ri2*s3 - qq*s1 )
        a2 = chky*gm1*( uv2*s1 - ( uu*s2 + vv*s3) + s4 )
        b(j,1) = a4                  +spec*s1
        b(j,2) = ri1*a2 + uu*a4      +spec*s2
        b(j,3) = ri2*a2 + vv*a4      +spec*s3
        b(j,4) = qq*a2 + ge(j,k1)*a4 +spec*s4
  122 continue
c
c..bi-diagonal inversion
c     
cdir$ ivdep
        do 124 j = ms(m),me(m)
          k  = m-j
          svt = tscale(j,k)
          sv = svt*0.5

          uu=q(j,k,2)
          vv=q(j,k,3)
         
          uv2=0.5*(uu*uu+vv*vv)
          cjkl=sqrt(ggm1*(q(j,k,4)-uv2))
          do n=1,4
             tmp(n)=a(j,n)+b(j,n)
          enddo
         
          a(j,1) =uv2*tmp(1)-uu*tmp(2)-vv*tmp(3)+tmp(4)
          a(j,2) =uu*a(j,1)
          a(j,3) =vv*a(j,1)
          a(j,4) =ge(j,k)*a(j,1)

          bSq = Mp**2/(bt(j,k)-Mp**2*(bt(j,k)-1))

          term=gm1*(bSq-1.)/(cjkl*cjkl)
          do n=1,4
             a(j,n)=a(j,n)*term
             a(j,n)=a(j,n)+tmp(n)
          enddo
         
          do n=1,4
             s(j,k,n)=(s(j,k,n)+sv*(a(j,n)))*d(j,k,n)
          enddo

 124    continue
c
  120 continue
c
c..backward sweep
c..loop on hyper-plane
c
      do 220 m = je+ke,js+ks,-1
c
c..setup a contribution in j-direction
c
      do 221 j = ms(m),me(m)
        k  = m-j
        j1 = j+1
c
        uu  = q(j1,k,2)
        vv  = q(j1,k,3)
        uv2 = uv(j1,k)
        ri1 = xx(j1,k)
        ri2 = xy(j1,k)
        qq  = qx(j1,k)
        cc  = cx(j1,k)
        qqx = -ri1*ug(j1,k)-ri2*vg(j1,k) + qq
        sp1 = abs(qqx) + cc
        sm2 = eps2*sp1
        chkx= 0.5 - sign( 0.5, qqx-cc )
        spec= chkx*( qqx - sp1 ) - sm2
c
        s1 = s(j1,k,1)
        s2 = s(j1,k,2)
        s3 = s(j1,k,3)
        s4 = s(j1,k,4)
c
        a4 = chkx*( ri1*s2 + ri2*s3 - qq*s1 )
        a2 = chkx*gm1*( uv2*s1 - ( uu*s2 + vv*s3 ) + s4 )
        a(j,1) = a4                  +spec*s1
        a(j,2) = ri1*a2 + uu*a4      +spec*s2
        a(j,3) = ri2*a2 + vv*a4      +spec*s3
        a(j,4) = qq*a2 + ge(j1,k)*a4 +spec*s4
  221 continue
c
c..setup b contribution in k-direction
c
      do 222 j = ms(m),me(m)
        k  = m-j
        k1 = k+1
c
        uu  = q(j,k1,2)
        vv  = q(j,k1,3)
        uv2 = uv(j,k1)
        ri1 = yx(j,k1)
        ri2 = yy(j,k1)
        qq  = qy(j,k1)
        cc  = cy(j,k1)
        vnu = vn(j,k1)
        qqy = -ri1*ug(j,k1)-ri2*vg(j,k1) + qq
        sp1 = abs(qqy) + cc
        sm2 = eps2*sp1
        chky= 0.5 - sign( 0.5, qqy-cc )
        spec= chky*( qq - sp1 ) - vnu - sm2
c
        s1 = s(j,k1,1)
        s2 = s(j,k1,2)
        s3 = s(j,k1,3)
        s4 = s(j,k1,4)
c
        a4 = chky*( ri1*s2 + ri2*s3 - qq*s1 )
        a2 = chky*gm1*( uv2*s1 - ( uu*s2 + vv*s3 ) + s4 )
        b(j,1) = a4                  +spec*s1
        b(j,2) = ri1*a2 + uu*a4      +spec*s2
        b(j,3) = ri2*a2 + vv*a4      +spec*s3
        b(j,4) = qq*a2 + ge(j,k1)*a4 +spec*s4
  222 continue
c
c..bi-diagonal inversion
c     
cdir$ ivdep
      do 223 j = ms(m),me(m)
        k  = m-j
        svt = tscale(j,k)
        sv = svt*0.5

        uu=q(j,k,2)
        vv=q(j,k,3)
       
        uv2=0.5*(uu*uu+vv*vv)
        cjkl=sqrt(ggm1*(q(j,k,4)-uv2))
        do n=1,4
           tmp(n)=a(j,n)+b(j,n)
        enddo
       
        a(j,1) =uv2*tmp(1)-uu*tmp(2)-vv*tmp(3)+tmp(4)
        a(j,2) =uu*a(j,1)
        a(j,3) =vv*a(j,1)
        a(j,4) =ge(j,k1)*a(j,1)

        bSq = Mp**2/(bt(j,k)-Mp**2*(bt(j,k)-1))

        term=gm1*(bSq-1)/(cjkl*cjkl)
        do n=1,4
           a(j,n)=a(j,n)*term
           a(j,n)=a(j,n)+tmp(n)
        enddo
        
        do n=1,4
           s(j,k,n)=s(j,k,n)-sv*a(j,n)*d(j,k,n)
        enddo

        !di = d(j,k)
        !s(j,k,1) = s(j,k,1) - sv*( a(j,1)+b(j,1) )*d(j,k,1)
        !s(j,k,2) = s(j,k,2) - sv*( a(j,2)+b(j,2) )*d(j,k,2)
        !s(j,k,3) = s(j,k,3) - sv*( a(j,3)+b(j,3) )*d(j,k,3)
        !s(j,k,4) = s(j,k,4) - sv*( a(j,4)+b(j,4) )*d(j,k,4)
  223 continue
c
  220 continue
 1000 continue
c
c..restore density , u , v in conservative variables             
c                                                                       
      do 4 k = ks-1, ke+1
      do 4 j = js-1, je+1
        q(j,k,1) = q(j,k,1) / q(j,k,nq)
        q(j,k,2) = q(j,k,2)*q(j,k,1)
        q(j,k,3) = q(j,k,3)*q(j,k,1)
        q(j,k,4) = q(j,k,4)*q(j,k,1)
    4 continue
c                                                                       
      return
      end

c***********************************************************************
      subroutine ilu2d(q,s,jd,kd,js,je,ks,ke,xx,xy,yx,yy,ug,vg,turmu,
     >                 tscale)
c
c  calculate the implicit inversion of the lhs
c  this involves two bidiagonal scalar inversions
c
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************
      integer jd,kd,js,je,ks,ke
      real q(jd,kd,nq), s(jd,kd,nv),turmu(jd,kd),tscale(jd,kd)
      real ug(jd,kd), vg(jd,kd)
      real xx(jd,kd), xy(jd,kd), yx(jd,kd), yy(jd,kd)

      ! local variables
      real,allocatable :: d(:,:),st(:,:,:),tmp(:)
      real,allocatable :: a(:,:),c(:,:)
      real,allocatable :: uv(:,:),vn(:,:),ge(:,:)
      real,allocatable :: qx(:,:),qy(:,:),cx(:,:),cy(:,:)
      real,allocatable :: b(:,:)
      integer,allocatable :: ms(:),me(:)
      integer j,k,m,i,k1,jj
      integer itgs,j1
      real eps2,epsv,dj,uu,vv,uv2,cjkl,rj1,rj2
      real qq1,qqx,rr1,rk1,rk2,qq2,qqy,rr2,vnu,svt
      real ri1,ri2,qq,cc,sp1,sm2,chkx,spec,s1,s2
      real a2,a4,chky,sv,di,s3,s4
      real scale1,scale2,sav1,sav2,sav3,sav4,scale

      allocate(d(jd,kd),st(jd,kd,4),tmp(4))
      allocate(a(mdim,4),c(mdim,4))
      allocate(uv(jd,kd),vn(jd,kd),ge(jd,kd))
      allocate(qx(jd,kd),qy(jd,kd),cx(jd,kd),cy(jd,kd))
      allocate(b(mdim,4))
      allocate(ms(mdim*2),me(mdim*2))

c***  first executable statement

      eps2 = epse*2.
      epsv = 1. + eps2

c..set-up for hyper-plane loop

      do 1 m=js+ks,je+ke
        i     = max(js,m-ke)
        ms(m) = i
        me(m) = min(je,m-js)
    1 continue

c
c..store density,u,v in nonconservative variables                  
c
      do 2 k = ks-1, ke+1
      do 2 j = js-1, je+1
        dj       = 1.0 / q(j,k,1)
        q(j,k,2) = q(j,k,2)*dj
        q(j,k,3) = q(j,k,3)*dj
        q(j,k,4) = q(j,k,4)*dj
        q(j,k,1) = q(j,k,1)*q(j,k,nq)
    2 continue
c
c..in the wake, average points above and below
c..for faster convergence?
c
      if(half.ne.1.) then
        k  = 1
        k1 = k + 1
        do 11 j = 1, jtail1 - 1
            jj = jmax - j + 1
c
          uu  = q(j,k,2)
          vv  = q(j,k,3)
          uv2 = 0.5*( uu*uu + vv*vv )
          cjkl = sqrt( ggm1*( q(j,k,4) - uv2 ) )
c
          rj1 = xx(j,k)
          rj2 = xy(j,k)
          qq1 = rj1*uu + rj2*vv
          qqx = abs( -rj1*ug(j,k)-rj2*vg(j,k) + qq1 )
          rr1 = sqrt( rj1**2 + rj2**2 )
          rk1 = yx(j,k)
          rk2 = yy(j,k)
          qq2 = rk1*uu + rk2*vv
          qqy = abs( -rk1*ug(j,k)-rk2*vg(j,k) + qq2 )
          rr2 = sqrt( rk1**2 + rk2**2 )
          vnu = 2.0*rr2*rr2*(rmue+turmu(j,k))/(rey*q(j,k,1))
          svt = tscale(j,k)
     <
          d(j,k)=1.0/
     &     ( 1.0+svt*((qqx+qqy +cjkl*(rr1+rr2))*epsv +vnu) )
            scale = q(j,k1,nq)/q(j,k,nq)
            scale1= q(jj,k1,nq)/q(j,k,nq)
            scale2= q(j,k,nq)/q(jj,k,nq)
            sav1 = 0.5*(s(j,k1,1)*scale+s(jj,k1,1)*scale1)
            sav2 = 0.5*(s(j,k1,2)*scale+s(jj,k1,2)*scale1)
            sav3 = 0.5*(s(j,k1,3)*scale+s(jj,k1,3)*scale1)
            sav4 = 0.5*(s(j,k1,4)*scale+s(jj,k1,4)*scale1)
            s(j,k,1)= sav1*d(j,k)*.95
            s(j,k,2)= sav2*d(j,k)*.95
            s(j,k,3)= sav3*d(j,k)*.95
            s(j,k,4)= sav4*d(j,k)*.95
            s(jj,k,1)= sav1*scale2*d(j,k)*.95
            s(jj,k,2)= sav2*scale2*d(j,k)*.95
            s(jj,k,3)= sav3*scale2*d(j,k)*.95
            s(jj,k,4)= sav4*scale2*d(j,k)*.95
   11   continue
      endif
c
c..forward sweep
c..setup d, the diagonal term and and store some arrays
c
      do 111 k = ks-1,ke+1
      do 111 j = js-1,je+1
          uu  = q(j,k,2)
          vv  = q(j,k,3)
          uv2 = 0.5*( uu*uu + vv*vv )
          cjkl = sqrt( ggm1*( q(j,k,4) - uv2 ) )
c
          rj1 = xx(j,k)
          rj2 = xy(j,k)
          qq1 = rj1*uu + rj2*vv
          qqx = abs( -rj1*ug(j,k)-rj2*vg(j,k) + qq1 )
          rr1 = sqrt( rj1**2 + rj2**2 )
          rk1 = yx(j,k)
          rk2 = yy(j,k)
          qq2 = rk1*uu + rk2*vv
          qqy = abs( -rk1*ug(j,k)-rk2*vg(j,k) + qq2 )
          rr2 = sqrt( rk1**2 + rk2**2 )
          vnu = 2.0*rr2*rr2*(rmue+turmu(j,k))/(rey*q(j,k,1))
          svt = tscale(j,k)
          d(j,k)=1.0/
     &     ( 1.0+svt*((qqx+qqy +cjkl*(rr1+rr2))*epsv +vnu) )
c
          uv(j,k) = uv2
          qx(j,k) = qq1
          cx(j,k) = cjkl*rr1
          qy(j,k) = qq2
          cy(j,k) = cjkl*rr2
          vn(j,k) = vnu
          ge(j,k) = gamma*q(j,k,4) - gm1*uv2
  111 continue
c
      do 1000 itgs=1,1 
c
c..loop on hyper-plane
c
      do 120 m = js+ks,je+ke
c
c..setup a contribution in j-direction
c
      do 121 j = ms(m),me(m)
        k  = m-j
        j1 = j-1
c
        uu  = q(j1,k,2)
        vv  = q(j1,k,3)
        uv2 = uv(j1,k)
        ri1 = xx(j1,k)
        ri2 = xy(j1,k)
        qq  = qx(j1,k)
        cc  = cx(j1,k)
        qqx = -ri1*ug(j1,k)-ri2*vg(j1,k) + qq
        sp1 = abs(qqx) + cc
        sm2 = eps2*sp1
        chkx= 0.5 + sign( 0.5, qqx+cc )
        spec= chkx*( qqx + sp1 ) + sm2
c
        s1 = s(j1,k,1)
        s2 = s(j1,k,2)
        s3 = s(j1,k,3)
        s4 = s(j1,k,4)
c
        a4 = chkx*( ri1*s2 + ri2*s3 - qq*s1 )
        a2 = chkx*gm1*( uv2*s1 - ( uu*s2 + vv*s3) + s4 )
        a(j,1) = a4                  +spec*s1
        a(j,2) = ri1*a2 + uu*a4      +spec*s2
        a(j,3) = ri2*a2 + vv*a4      +spec*s3
        a(j,4) = qq*a2 + ge(j1,k)*a4 +spec*s4
  121 continue
c
c..setup b contribution in k-direction
c
      do 122 j = ms(m),me(m)
        k  = m-j
        k1 = k-1
c
        uu  = q(j,k1,2)
        vv  = q(j,k1,3)
        uv2 = uv(j,k1)
        ri1 = yx(j,k1)
        ri2 = yy(j,k1)
        qq  = qy(j,k1)
        cc  = cy(j,k1)
        vnu = vn(j,k1)
        qqy = -ri1*ug(j,k1)-ri2*vg(j,k1) + qq
        sp1 = abs(qqy) + cc
        sm2 = eps2*sp1
        chky= 0.5 + sign( 0.5, qqy+cc )
        spec= chky*( qqy + sp1 ) + vnu + sm2
c
        s1 = s(j,k1,1)
        s2 = s(j,k1,2)
        s3 = s(j,k1,3)
        s4 = s(j,k1,4)
c
        a4 = chky*( ri1*s2 + ri2*s3 - qq*s1 )
        a2 = chky*gm1*( uv2*s1 - ( uu*s2 + vv*s3) + s4 )
        b(j,1) = a4                  +spec*s1
        b(j,2) = ri1*a2 + uu*a4      +spec*s2
        b(j,3) = ri2*a2 + vv*a4      +spec*s3
        b(j,4) = qq*a2 + ge(j,k1)*a4 +spec*s4
  122 continue
c
c..bi-diagonal inversion
c     
cdir$ ivdep
        do 124 j = ms(m),me(m)
          k  = m-j
          svt = tscale(j,k)
          sv = svt*0.5
          di = d(j,k)
          s(j,k,1) = ( s(j,k,1) + sv*( a(j,1)+b(j,1) ) )*di
          s(j,k,2) = ( s(j,k,2) + sv*( a(j,2)+b(j,2) ) )*di
          s(j,k,3) = ( s(j,k,3) + sv*( a(j,3)+b(j,3) ) )*di
          s(j,k,4) = ( s(j,k,4) + sv*( a(j,4)+b(j,4) ) )*di
 124    continue
c
  120 continue
c
c..backward sweep
c..loop on hyper-plane
c
      do 220 m = je+ke,js+ks,-1
c
c..setup a contribution in j-direction
c
      do 221 j = ms(m),me(m)
        k  = m-j
        j1 = j+1
c
        uu  = q(j1,k,2)
        vv  = q(j1,k,3)
        uv2 = uv(j1,k)
        ri1 = xx(j1,k)
        ri2 = xy(j1,k)
        qq  = qx(j1,k)
        cc  = cx(j1,k)
        qqx = -ri1*ug(j1,k)-ri2*vg(j1,k) + qq
        sp1 = abs(qqx) + cc
        sm2 = eps2*sp1
        chkx= 0.5 - sign( 0.5, qqx-cc )
        spec= chkx*( qqx - sp1 ) - sm2
c
        s1 = s(j1,k,1)
        s2 = s(j1,k,2)
        s3 = s(j1,k,3)
        s4 = s(j1,k,4)
c
        a4 = chkx*( ri1*s2 + ri2*s3 - qq*s1 )
        a2 = chkx*gm1*( uv2*s1 - ( uu*s2 + vv*s3 ) + s4 )
        a(j,1) = a4                  +spec*s1
        a(j,2) = ri1*a2 + uu*a4      +spec*s2
        a(j,3) = ri2*a2 + vv*a4      +spec*s3
        a(j,4) = qq*a2 + ge(j1,k)*a4 +spec*s4
  221 continue
c
c..setup b contribution in k-direction
c
      do 222 j = ms(m),me(m)
        k  = m-j
        k1 = k+1
c
        uu  = q(j,k1,2)
        vv  = q(j,k1,3)
        uv2 = uv(j,k1)
        ri1 = yx(j,k1)
        ri2 = yy(j,k1)
        qq  = qy(j,k1)
        cc  = cy(j,k1)
        vnu = vn(j,k1)
        qqy = -ri1*ug(j,k1)-ri2*vg(j,k1) + qq
        sp1 = abs(qqy) + cc
        sm2 = eps2*sp1
        chky= 0.5 - sign( 0.5, qqy-cc )
        spec= chky*( qq - sp1 ) - vnu - sm2
c
        s1 = s(j,k1,1)
        s2 = s(j,k1,2)
        s3 = s(j,k1,3)
        s4 = s(j,k1,4)
c
        a4 = chky*( ri1*s2 + ri2*s3 - qq*s1 )
        a2 = chky*gm1*( uv2*s1 - ( uu*s2 + vv*s3 ) + s4 )
        b(j,1) = a4                  +spec*s1
        b(j,2) = ri1*a2 + uu*a4      +spec*s2
        b(j,3) = ri2*a2 + vv*a4      +spec*s3
        b(j,4) = qq*a2 + ge(j,k1)*a4 +spec*s4
  222 continue
c
c..bi-diagonal inversion
c     
cdir$ ivdep
      do 223 j = ms(m),me(m)
        k  = m-j
        svt = tscale(j,k)
        sv = svt*0.5
        di = d(j,k)
        s(j,k,1) = s(j,k,1) - sv*( a(j,1)+b(j,1) )*di
        s(j,k,2) = s(j,k,2) - sv*( a(j,2)+b(j,2) )*di
        s(j,k,3) = s(j,k,3) - sv*( a(j,3)+b(j,3) )*di
        s(j,k,4) = s(j,k,4) - sv*( a(j,4)+b(j,4) )*di
  223 continue
c
  220 continue
 1000 continue
c
c..restore density , u , v in conservative variables             
c                                                                       
      do 4 k = ks-1, ke+1
      do 4 j = js-1, je+1
        q(j,k,1) = q(j,k,1) / q(j,k,nq)
        q(j,k,2) = q(j,k,2)*q(j,k,1)
        q(j,k,3) = q(j,k,3)*q(j,k,1)
        q(j,k,4) = q(j,k,4)*q(j,k,1)
    4 continue
c                                                                       
      return
      end

c***********************************************************************
      subroutine arc2d(q,s,jd,kd,js,je,ks,ke,xx,xy,yx,yy,ug,vg,turmu,
     c                                                iblank,tscale,bt)
c
c  calculate the implicit inversion of the lhs
c  this involves two bidiagonal scalar inversions
c
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************
      integer jd,kd,js,je,ks,ke
      integer n
      integer iblank(jd,kd)
      real q(jd,kd,nq), s(jd,kd,nv),turmu(jd,kd),tscale(jd,kd),bt(jd,kd)
      real ug(jd,kd), vg(jd,kd)
      real xx(jd,kd), xy(jd,kd), yx(jd,kd), yy(jd,kd)

      ! local variables

      real,allocatable :: p(:,:),d2px(:,:),d2py(:,:)
      real,allocatable :: uv(:,:),vn(:,:),ge(:,:)
      real,allocatable :: qx(:,:),qy(:,:),cx(:,:),cy(:,:)

      integer j,k
      real dj,uu,vv,uv2,cjkl,rj1,rj2
      real qq1,rr1,rk1,rk2,qq2,rr2,vnu
      real c2i
      real a1,a2,a3,a4,a5,a6

      allocate(p(jd,kd),d2px(jd,kd),d2py(jd,kd))
      allocate(uv(jd,kd),vn(jd,kd),ge(jd,kd))
      allocate(qx(jd,kd),qy(jd,kd),cx(jd,kd),cy(jd,kd))

c..store density,u,v in nonconservative variables                  
c
      do 2 k = ks-1, ke+1
      do 2 j = js-1, je+1
        dj       = 1.0 / q(j,k,1)
        q(j,k,2) = q(j,k,2)*dj
        q(j,k,3) = q(j,k,3)*dj
        q(j,k,4) = q(j,k,4)*dj
        q(j,k,1) = q(j,k,1)*q(j,k,nq)
    2 continue
c
c..multiply by T_si_inverse

      do 111 k = ks-1,ke+1
      do 111 j = js-1,je+1
         uu  = q(j,k,2)
         vv  = q(j,k,3)
         uv2 = 0.5*( uu*uu + vv*vv )
         cjkl   = sqrt( ggm1*( q(j,k,4) - uv2 ) )
         c2i    = 1.0/(cjkl*cjkl)
         p(j,k) = ( q(j,k,4) - uv2 )*gm1*q(j,k,1)
c
         rj1 = xx(j,k)
         rj2 = xy(j,k)
         qq1 = rj1*uu + rj2*vv
         qx(j,k) = qq1 - ( rj1*ug(j,k) + rj2*vg(j,k) )
         rr1 = sqrt( rj1**2 + rj2**2 )
         rj1 = rj1/rr1
         rj2 = rj2/rr1

         rk1 = yx(j,k)
         rk2 = yy(j,k)
         qq2 = rk1*uu + rk2*vv
         qy(j,k) = qq2 - ( rk1*ug(j,k) + rk2*vg(j,k) )
         rr2 = sqrt( rk1**2 + rk2**2 )
         rk1 = rk1/rr2
         rk2 = rk2/rr2

         vnu = (rmue+turmu(j,k))/(rey*q(j,k,1))

         uv(j,k) = uv2
         cx(j,k) = cjkl*rr1
         cy(j,k) = cjkl*rr2
         vn(j,k) = vnu
         ge(j,k) = gamma*q(j,k,4) - gm1*uv2

         a1 = s(j,k,2)*uu + s(j,k,3)*vv - s(j,k,4)
         a1 = a1*gm1*c2i
         a1 = a1+s(j,k,1)*( 1.0 - uv2*gm1*c2i )

         a2 =(rj1*vv-rj2*uu)*s(j,k,1)+rj2*s(j,k,2)-rj1*s(j,k,3)
  
         a3 = uv2*s(j,k,1)-s(j,k,2)*uu-s(j,k,3)*vv+s(j,k,4)
         a3 = a3*gm1*c2i

         a4 = qq1*s(j,k,1)/rr1-rj1*s(j,k,2)-rj2*s(j,k,3)
         a4 = a4

         s(j,k,1) = a1
         s(j,k,2) = a2
         s(j,k,3) = 0.5*(a3-a4/cjkl)
         s(j,k,4) = 0.5*(a3+a4/cjkl)

  111 continue

      do 112 k = ks-1,ke+1
      do 112 j = js,je
         d2px(j,k) = abs( p(j+1,k) - 2.*p(j,k) + p(j-1,k) )
         d2px(j,k) = d2px(j,k)/abs( p(j+1,k) + 2.*p(j,k) + p(j-1,k) )
  112 continue

      j=js-1
      do 113 k = ks-1,ke+1
         d2px(j,k) = d2px(j+1,k)
         d2px(j,k) = 0
  113 continue

      j=je+1
      do 114 k = ks-1,ke+1
         d2px(j,k) = d2px(j-1,k)
         d2px(j,k) = 0
  114 continue

c   finding lhs inverse
      if (ilhs.eq.2) then
         call lhsinv(q,s,qx,cx,d2px,js-1,je+1,ks-1,ke+1,1,xx,xy,
     >               vn,jd,kd,iblank,tscale,bt)
      elseif (ilhs.eq.3) then
         call lhsinv_up(q,s,qx,cx,js-1,je+1,ks-1,ke+1,1,xx,
     >                  xy,vn,jd,kd,iblank,tscale,bt)
      endif
       

c multiply by N_inverse, s(j,k,1) is unchanged

      do 121 k = ks-1,ke+1
      do 121 j = js-1,je+1

         cjkl = sqrt( ggm1*( q(j,k,4) - uv(j,k) ) )
        
         rj1  = xx(j,k)
         rj2  = xy(j,k)
         rr1  = sqrt(rj1*rj1+rj2*rj2)
        
         rk1  = yx(j,k)
         rk2  = yy(j,k)
         rr2  = sqrt(rk1*rk1+rk2*rk2)

         a1 = (rj1*rk1+rj2*rk2)/rr1/rr2
         a2 = (rj1*rk2-rj2*rk1)/rr1/rr2
         a3 = cjkl

         a4 = a1*s(j,k,2)+(cjkl*s(j,k,3)-cjkl*s(j,k,4))*a2
         a5 =-0.5*(a2*s(j,k,2)/cjkl-(a1+1)*s(j,k,3)+(a1-1)*s(j,k,4))
         a6 = 0.5*(a2*s(j,k,2)/cjkl-(a1-1)*s(j,k,3)+(a1+1)*s(j,k,4))


         s(j,k,2)=a4
         s(j,k,3)=a5
         s(j,k,4)=a6
  121 continue

      do 122 k = ks,ke
      do 122 j = js-1,je+1
         d2py(j,k) = abs( p(j,k+1) - 2.*p(j,k) + p(j,k-1) ) 
         d2py(j,k) = d2py(j,k)/abs( p(j,k+1) + 2.*p(j,k) + p(j,k-1) )
  122 continue

      k=ks-1
      do 123 j = js-1,je+1
         d2py(j,k) = d2py(j,k+1)
         d2py(j,k) = 0
  123 continue

      k=ke+1
      do 124 j = js-1,je+1
         d2py(j,k) = d2py(j,k-1)
         d2py(j,k) = 0
  124 continue

c...finding lhs inverse
      if (ilhs.eq.2) then
         call lhsinv(q,s,qy,cy,d2py,ks-1,ke+1,js-1,je+1,2,yx,yy,
     >               vn,jd,kd,iblank,tscale,bt)
      elseif (ilhs.eq.3) then
         call lhsinv_up(q,s,qy,cy,ks-1,ke+1,js-1,je+1,2,yx,
     >                  yy,vn,jd,kd,iblank,tscale,bt)
      endif
  

c...multiply by T_eta
      do 131 k = ks-1,ke+1
      do 131 j = js-1,je+1
         uu   = q(j,k,2)
         vv   = q(j,k,3)
         cjkl = sqrt(ggm1*(q(j,k,4)-uv(j,k)))
         c2i  = 1.0/(cjkl*cjkl)
         
         rk1  = yx(j,k)
         rk2  = yy(j,k)
         qq2  = qy(j,k) + ( rk1*ug(j,k) + rk2*vg(j,k) )
         rr2  = sqrt(rk1*rk1+rk2*rk2)
         rk1  = rk1/rr2
         rk2  = rk2/rr2
         
         a1 = s(j,k,1)+s(j,k,3)+s(j,k,4)

         a2 = uu*s(j,k,1)+rk2*s(j,k,2)+(uu+rk1*cjkl)*s(j,k,3)
         a2 = a2 + (uu-rk1*cjkl)*s(j,k,4)
 
         a3 = vv*s(j,k,1)-rk1*s(j,k,2)+(vv+rk2*cjkl)*s(j,k,3)
         a3 = a3 + (vv-rk2*cjkl)*s(j,k,4)

         a4 = uv(j,k)*s(j,k,1)+(rk2*uu-rk1*vv)*s(j,k,2)
         a4 = a4 + (ge(j,k)+cjkl*qq2/rr2)*s(j,k,3)
         a4 = a4 + (ge(j,k)-cjkl*qq2/rr2)*s(j,k,4)
         
         s(j,k,1) = a1
         s(j,k,2) = a2
         s(j,k,3) = a3
         s(j,k,4) = a4
  131 continue
c
c..restore density , u , v in conservative variables             
c                                                                       
      do 4 k = ks-1, ke+1
      do 4 j = js-1, je+1
        q(j,k,1) = q(j,k,1) / q(j,k,nq)
        q(j,k,2) = q(j,k,2)*q(j,k,1)
        q(j,k,3) = q(j,k,3)*q(j,k,1)
        q(j,k,4) = q(j,k,4)*q(j,k,1)
    4 continue
c                                                                       
      return
      end

c***********************************************************************
      subroutine arc2d_precon(q,s,jd,kd,js,je,ks,ke,xx,xy,yx,yy,ug,vg,turmu,
     c                                                 iblank,tscale,bt)
c
c  calculate the implicit inversion of the lhs
c  this involves two bidiagonal scalar inversions
c
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************
      integer jd,kd,js,je,ks,ke
      integer n
      integer iblank(jd,kd)
      real q(jd,kd,nq), s(jd,kd,nv),turmu(jd,kd),tscale(jd,kd),bt(jd,kd)
      real ug(jd,kd), vg(jd,kd)
      real xx(jd,kd), xy(jd,kd), yx(jd,kd), yy(jd,kd)

      ! local variables

      real,allocatable :: p(:,:),d2px(:,:),d2py(:,:)
      real,allocatable :: uv(:,:),vn(:,:),ge(:,:)
      real,allocatable :: qx(:,:),qy(:,:),cx(:,:),cy(:,:)

      integer j,k
      real dj,uu,vv,uv2,cjkl,rj1,rj2
      real qqmxt,qq1,rr1,rk1,rk2,qq2,rr2,vnu
      real c2i,X,Y,Z,X1,Y1,Z1,bSq
      real a1,a2,a3,a4,a5,a6

      allocate(p(jd,kd),d2px(jd,kd),d2py(jd,kd))
      allocate(uv(jd,kd),vn(jd,kd),ge(jd,kd))
      allocate(qx(jd,kd),qy(jd,kd),cx(jd,kd),cy(jd,kd))

c..store density,u,v in nonconservative variables                  
c
      do 2 k = ks-1, ke+1
      do 2 j = js-1, je+1
        dj       = 1.0 / q(j,k,1)
        q(j,k,2) = q(j,k,2)*dj
        q(j,k,3) = q(j,k,3)*dj
        q(j,k,4) = q(j,k,4)*dj
        q(j,k,1) = q(j,k,1)*q(j,k,nq)
    2 continue
c
c..multiply by T_si_inverse

      do 111 k = ks-1,ke+1
      do 111 j = js-1,je+1
         uu  = q(j,k,2)
         vv  = q(j,k,3)
         uv2 = 0.5*( uu*uu + vv*vv )
         cjkl   = sqrt( ggm1*( q(j,k,4) - uv2 ) )
         c2i    = 1.0/(cjkl*cjkl)
         p(j,k) = ( q(j,k,4) - uv2 )*gm1*q(j,k,1)
c
         rj1 = xx(j,k)
         rj2 = xy(j,k)
         qq1 = rj1*( uu - ug(j,k) ) + rj2*( vv - vg(j,k) )
         rr1 = sqrt( rj1**2 + rj2**2 )
         rj1 = rj1/rr1
         rj2 = rj2/rr1

         rk1 = yx(j,k)
         rk2 = yy(j,k)
         qq2 = rk1*( uu - ug(j,k) ) + rk2*( vv - vg(j,k) )
         rr2 = sqrt( rk1**2 + rk2**2 )
         rk1 = rk1/rr2
         rk2 = rk2/rr2

         vnu = (rmue+turmu(j,k))/(rey*q(j,k,1))

         uv(j,k) = uv2
         qx(j,k) = qq1
         cx(j,k) = cjkl*rr1
         qy(j,k) = qq2
         cy(j,k) = cjkl*rr2
         vn(j,k) = vnu
         ge(j,k) = gamma*q(j,k,4) - gm1*uv2

         bSq = Mp**2/(bt(j,k)-Mp**2*(bt(j,k)-1))

         X = sqrt(((1.-bSq)*qq1)*((1.-bSq)*qq1)
     c       +4.*bSq*cx(j,k)*cx(j,k))/rr1
         Y = 0.5*((1.0 - bSq)*qq1/rr1+X)
         Z = 0.5*((1.0 - bSq)*qq1/rr1-X)

         qqmxt = qq1 + ug(j,k)*rj1*rr1 + vg(j,k)*rj2*rr1

         a1 = s(j,k,2)*uu + s(j,k,3)*vv - s(j,k,4)
         a1 = a1*gm1*c2i
         a1 = a1+s(j,k,1)*( 1.0 - uv2*gm1*c2i )

         a2 =(rj1*vv-rj2*uu)*s(j,k,1)+rj2*s(j,k,2)-rj1*s(j,k,3)
  
         a3 = uv2*s(j,k,1)-s(j,k,2)*uu-s(j,k,3)*vv+s(j,k,4)
         a3 = a3*gm1*c2i

         a4 = qqmxt*s(j,k,1)/rr1-rj1*s(j,k,2)-rj2*s(j,k,3)
         a4 = a4*bSq

         s(j,k,1) = a1
         s(j,k,2) = a2
         s(j,k,3) =-(Z*a3+a4)/X
         s(j,k,4) = (Y*a3+a4)/X

  111 continue

      do 112 k = ks-1,ke+1
      do 112 j = js,je
         d2px(j,k) = abs( p(j+1,k) - 2.*p(j,k) + p(j-1,k) )
         d2px(j,k) = d2px(j,k)/abs( p(j+1,k) + 2.*p(j,k) + p(j-1,k) )
  112 continue

      j=js-1
      do 113 k = ks-1,ke+1
         d2px(j,k) = d2px(j+1,k)
         d2px(j,k) = 0
  113 continue

      j=je+1
      do 114 k = ks-1,ke+1
         d2px(j,k) = d2px(j-1,k)
         d2px(j,k) = 0
  114 continue

c   finding lhs inverse
      if (ilhs.eq.2) then
         call lhsinv(q,s,qx,cx,d2px,js-1,je+1,ks-1,ke+1,1,xx,xy,
     >               vn,jd,kd,iblank,tscale,bt)
      elseif (ilhs.eq.3) then
         call lhsinv_up(q,s,qx,cx,js-1,je+1,ks-1,ke+1,1,xx,
     >                  xy,vn,jd,kd,iblank,tscale,bt)
      endif
       

c multiply by N_inverse, s(j,k,1) is unchanged

      do 121 k = ks-1,ke+1
      do 121 j = js-1,je+1

         cjkl = sqrt( ggm1*( q(j,k,4) - uv(j,k) ) )
        
         rj1  = xx(j,k)
         rj2  = xy(j,k)
         rr1  = sqrt(rj1*rj1+rj2*rj2)
        
         rk1  = yx(j,k)
         rk2  = yy(j,k)
         rr2  = sqrt(rk1*rk1+rk2*rk2)

         bSq = Mp**2/(bt(j,k)-Mp**2*(bt(j,k)-1))

         X = sqrt((( 1.0 - bSq )*qx(j,k))*(( 1.0 - bSq )*qx(j,k))
     c      +4.*bSq*cx(j,k)*cx(j,k))/rr1
         Y = 0.5*(( 1.0 - bSq )*qx(j,k)/rr1+X)
         Z = 0.5*(( 1.0 - bSq )*qx(j,k)/rr1-X)

         X1 = sqrt((( 1.0 - bSq )*qy(j,k))*(( 1.0 - bSq )*qy(j,k))
     c       +4.*bSq*cy(j,k)*cy(j,k))/rr2
         Y1 = 0.5*(( 1.0 - bSq )*qy(j,k)/rr2+X1)
         Z1 = 0.5*(( 1.0 - bSq )*qy(j,k)/rr2-X1)

         a1 = (rj1*rk1+rj2*rk2)/rr1/rr2
         a2 = (rj1*rk2-rj2*rk1)/rr1/rr2
         a3 = cjkl

         a4 = a1*s(j,k,2)+(Y*s(j,k,3)+Z*s(j,k,4))*a2/bSq
         a5 =-(a2*bSq*s(j,k,2)-(Y*a1-Z1)*s(j,k,3)-(Z*a1-Z1)*s(j,k,4))/X1
         a6 = (a2*bSq*s(j,k,2)-(Y*a1-Y1)*s(j,k,3)-(Z*a1-Y1)*s(j,k,4))/X1


         s(j,k,2)=a4
         s(j,k,3)=a5
         s(j,k,4)=a6
  121 continue

      do 122 k = ks,ke
      do 122 j = js-1,je+1
         d2py(j,k) = abs( p(j,k+1) - 2.*p(j,k) + p(j,k-1) ) 
         d2py(j,k) = d2py(j,k)/abs( p(j,k+1) + 2.*p(j,k) + p(j,k-1) )
  122 continue

      k=ks-1
      do 123 j = js-1,je+1
         d2py(j,k) = d2py(j,k+1)
         d2py(j,k) = 0
  123 continue

      k=ke+1
      do 124 j = js-1,je+1
         d2py(j,k) = d2py(j,k-1)
         d2py(j,k) = 0
  124 continue

c...finding lhs inverse
      if (ilhs.eq.2) then
         call lhsinv(q,s,qy,cy,d2py,ks-1,ke+1,js-1,je+1,2,yx,yy,
     >               vn,jd,kd,iblank,tscale,bt)
      elseif (ilhs.eq.3) then
         call lhsinv_up(q,s,qy,cy,ks-1,ke+1,js-1,je+1,2,yx,
     >                  yy,vn,jd,kd,iblank,tscale,bt)
      endif
  

c...multiply by T_eta
      do 131 k = ks-1,ke+1
      do 131 j = js-1,je+1
         uu   = q(j,k,2)
         vv   = q(j,k,3)
         cjkl = sqrt(ggm1*(q(j,k,4)-uv(j,k)))
         c2i  = 1.0/(cjkl*cjkl)
         
         rk1  = yx(j,k)
         rk2  = yy(j,k)
         qq2  = qy(j,k)
         rr2  = sqrt(rk1*rk1+rk2*rk2)
         rk1  = rk1/rr2
         rk2  = rk2/rr2
         
         bSq = Mp**2/(bt(j,k)-Mp**2*(bt(j,k)-1))

         X = sqrt(((1.-bSq)*qq2)*((1.-bSq)*qq2)
     c       +4.*bSq*cy(j,k)*cy(j,k))/rr2
         Y = 0.5*((1.-bSq)*qq2/rr2+X)
         Z = 0.5*((1.-bSq)*qq2/rr2-X)
         
         qqmxt = qq2 + ug(j,k)*rk1*rr2 + vg(j,k)*rk2*rr2

         a1 = s(j,k,1)+s(j,k,3)+s(j,k,4)

         a2 = uu*s(j,k,1)+rk2*s(j,k,2)+(uu+rk1*Y/bSq)*s(j,k,3)
         a2 = a2 + (uu+rk1*Z/bSq)*s(j,k,4)
 
         a3 = vv*s(j,k,1)-rk1*s(j,k,2)+(vv+rk2*Y/bSq)*s(j,k,3)
         a3 = a3 + (vv+rk2*Z/bSq)*s(j,k,4)

         a4 = uv(j,k)*s(j,k,1)+(rk2*uu-rk1*vv)*s(j,k,2)
         a4 = a4 + (ge(j,k)+Y*qqmxt/rr2/bSq)*s(j,k,3)
         a4 = a4 + (ge(j,k)+Z*qqmxt/rr2/bSq)*s(j,k,4)
         
         s(j,k,1) = a1
         s(j,k,2) = a2
         s(j,k,3) = a3
         s(j,k,4) = a4
  131 continue
c
c..restore density , u , v in conservative variables             
c                                                                       
      do 4 k = ks-1, ke+1
      do 4 j = js-1, je+1
        q(j,k,1) = q(j,k,1) / q(j,k,nq)
        q(j,k,2) = q(j,k,2)*q(j,k,1)
        q(j,k,3) = q(j,k,3)*q(j,k,1)
        q(j,k,4) = q(j,k,4)*q(j,k,1)
    4 continue
c                                                                       
      return
      end

c**********************************************************************
      subroutine lhsinv_up(q,s,qn,cn,ms,me,ns,ne,idir,xn,yn,vnu,
     >			   jd,kd,iblank,tscale,bt)

c**********************************************************************
      use params_global
c**********************************************************************
      implicit none
c**********************************************************************
      integer ms,me,ns,ne,idir,jd,kd
      integer iblank(jd,kd)
      real q(jd,kd,nq),s(jd,kd,nv),tscale(jd,kd),bt(jd,kd)
      real qn(jd,kd),cn(jd,kd),xn(jd,kd),yn(jd,kd),vnu(jd,kd)

      !local variables
      integer,allocatable :: iblnk(:)
      real,allocatable :: jcbn(:),tscal(:)
      real,allocatable :: vistrm1(:),vistrm2(:),vistrm3(:),g(:)
      real,allocatable :: a(:),b(:),c(:),d(:),e(:),f(:)
      real,allocatable :: a1(:),b1(:),c1(:),d1(:),e1(:)
      real,allocatable :: diag_plus(:,:),diag_minus(:,:)

      integer i,m,n
      real svt,dis2,dis4,bSq,X
      real c2,c2m,c4m2,c4m,c4,c4p
      real eig1,eig2,eig3,epsval,eps,fac

      allocate(iblnk(me))
      allocate(jcbn(me),tscal(me))
      allocate(vistrm1(me),vistrm2(me),vistrm3(me),g(me))
      allocate(a(me),b(me),c(me),d(me),e(me),f(me))
      allocate(a1(me),b1(me),c1(me),d1(me),e1(me))
      allocate(diag_plus(me,nq),diag_minus(me,nq))

c..some initialization

      epsval  = 0.05
      fac     = 1.05

      do n = ns,ne
         if (idir.eq.1) then
            do m = ms,me
               tscal(m) = tscale(m,n)
               jcbn(m)  = q(m,n,nq)
               iblnk(m) = max(iblank(m,n),0)
               bSq = Mp**2/(bt(m,n)-Mp**2*(bt(m,n)-1))
               X = sqrt(((1.-bSq)*qn(m,n))*((1.-bSq)*qn(m,n))
     c            +4.*bSq*cn(m,n)*cn(m,n))
               eig1 = qn(m,n)
               eig2 = 0.5*((bSq+1.)*qn(m,n)+X)
               eig3 = 0.5*((bSq+1.)*qn(m,n)-X)
               eps  = epsval*sqrt( xn(m,n)**2 + yn(m,n)**2 )
               diag_plus(m,1)  = 0.5*(eig1 + fac*sqrt(eig1**2 + eps**2))
               diag_plus(m,2)  = 0.5*(eig1 + fac*sqrt(eig1**2 + eps**2))
               diag_plus(m,3)  = 0.5*(eig2 + fac*sqrt(eig2**2 + eps**2))
               diag_plus(m,4)  = 0.5*(eig3 + fac*sqrt(eig3**2 + eps**2))
               diag_minus(m,1) = 0.5*(eig1 - fac*sqrt(eig1**2 + eps**2))
               diag_minus(m,2) = 0.5*(eig1 - fac*sqrt(eig1**2 + eps**2))
               diag_minus(m,3) = 0.5*(eig2 - fac*sqrt(eig2**2 + eps**2))
               diag_minus(m,4) = 0.5*(eig3 - fac*sqrt(eig3**2 + eps**2))
            enddo

            do m = ms+1,me-1
               vistrm1(m) = (xn(m,n)*xn(m,n)+yn(m,n)*yn(m,n))/q(m,n,nq)
               vistrm1(m) = vistrm1(m) + (xn(m+1,n)*xn(m+1,n)
     c                                 +yn(m+1,n)*yn(m+1,n))/q(m+1,n,nq)
               vistrm1(m) = vistrm1(m)*q(m,n,nq)*vnu(m,n)

               vistrm3(m) = (xn(m,n)*xn(m,n)+yn(m,n)*yn(m,n))/q(m,n,nq)
               vistrm3(m) = vistrm3(m) + (xn(m-1,n)*xn(m-1,n)
     c                                 +yn(m-1,n)*yn(m-1,n))/q(m-1,n,nq)
               vistrm3(m)= vistrm3(m)*q(m,n,nq)*vnu(m,n)

               vistrm2(m) = 0.5*(vistrm1(m) + vistrm3(m))
            enddo
            m = ms
            vistrm1(m) = 0
            vistrm2(m) = 0
            vistrm3(m) = 0
            m = me
            vistrm1(m) = 0
            vistrm2(m) = 0
            vistrm3(m) = 0

         elseif (idir.eq.2) then

            do m = ms,me
               tscal(m) = tscale(n,m)
               jcbn(m)  = q(n,m,nq)
               iblnk(m) = max(iblank(n,m),0)
               bSq = Mp**2/(bt(n,m)-Mp**2*(bt(n,m)-1))
               X= sqrt(((1.-bSq)*qn(n,m))*((1.-bSq)*qn(n,m))
     c            +4.*bSq*cn(n,m)*cn(n,m))
               eig1 = qn(n,m)
               eig2 = 0.5*((bSq+1.)*qn(n,m)+X)
               eig3 = 0.5*((bSq+1.)*qn(n,m)-X)
               eps  = epsval*sqrt( xn(n,m)**2 + yn(n,m)**2 )
               diag_plus(m,1)  = 0.5*(eig1 + fac*sqrt(eig1**2 + eps**2))
               diag_plus(m,2)  = 0.5*(eig1 + fac*sqrt(eig1**2 + eps**2))
               diag_plus(m,3)  = 0.5*(eig2 + fac*sqrt(eig2**2 + eps**2))
               diag_plus(m,4)  = 0.5*(eig3 + fac*sqrt(eig3**2 + eps**2))
               diag_minus(m,1) = 0.5*(eig1 - fac*sqrt(eig1**2 + eps**2))
               diag_minus(m,2) = 0.5*(eig1 - fac*sqrt(eig1**2 + eps**2))
               diag_minus(m,3) = 0.5*(eig2 - fac*sqrt(eig2**2 + eps**2))
               diag_minus(m,4) = 0.5*(eig3 - fac*sqrt(eig3**2 + eps**2))
            enddo

            do m = ms+1,me-1
               vistrm1(m) = (xn(n,m)*xn(n,m)+yn(n,m)*yn(n,m))/q(n,m,nq)
               vistrm1(m) = vistrm1(m) + (xn(n,m+1)*xn(n,m+1)
     c                                 +yn(n,m+1)*yn(n,m+1))/q(n,m+1,nq)
               vistrm1(m) = vistrm1(m)*q(n,m,nq)*vnu(n,m)
         
               vistrm3(m) = (xn(n,m)*xn(n,m)+yn(n,m)*yn(n,m))/q(n,m,nq)
               vistrm3(m) = vistrm3(m) + (xn(n,m-1)*xn(n,m-1)
     c                                 +yn(n,m-1)*yn(n,m-1))/q(n,m-1,nq)
               vistrm3(m) = vistrm3(m)*q(n,m,nq)*vnu(n,m)

               vistrm2(m) = 0.5*(vistrm1(m)+vistrm3(m))
            enddo
            m = ms
            vistrm1(m) = 0
            vistrm2(m) = 0
            vistrm3(m) = 0
            m = me
            vistrm1(m) = 0
            vistrm2(m) = 0
            vistrm3(m) = 0
 
         endif
 
         do m = ms,me

            a1(m) = 0
            b1(m) = 0
            c1(m) = 0
            d1(m) = 0
            e1(m) = 0

c        viscous terms
            b1(m) = b1(m) - 0.5*vistrm1(m)
            c1(m) = c1(m) + vistrm2(m)
            d1(m) = d1(m) - 0.5*vistrm3(m)

         enddo
c
         do i = 1,nv
	   do m=ms,me
               a(m) = a1(m)
               b(m) = b1(m)
               c(m) = c1(m)
               d(m) = d1(m)
               e(m) = e1(m)

               if (idir.eq.1) then
                  f(m) = s(m,n,i)
               elseif (idir.eq.2) then
                  f(m) = s(n,m,i)
               endif

c           euler terms
               b(m) = b(m) - diag_plus(m,i)
               d(m) = d(m) + diag_minus(m,i)
               c(m) = c(m) + diag_plus(m,i) - diag_minus(m,i)

c           multiply by delta t
               if(m.le.me-1) then
			b(m) = b(m)*tscal(m+1)
	       else
			b(m) = b(m)*tscal(m)
	       endif
               if(m.ge.ms+1) then
			d(m) = d(m)*tscal(m-1)
	       else
			d(m) = d(m)*tscal(m)
	       endif
               c(m) = c(m)*tscal(m) + 1.0
            enddo

            do m = ms+1,me-1
               a(m-ms) = a(m) 
               b(m-ms) = b(m) 
               c(m-ms) = c(m) 
               d(m-ms) = d(m) 
               e(m-ms) = e(m) 
               f(m-ms) = f(m) 
            enddo
            d(1) = 0
            e(1) = 0
            e(2) = 0
            b(me-ms-1) = 0
            a(me-ms-1) = 0
            a(me-ms-2) = 0

            call pentadag(a,b,c,d,e,f,me-ms-1)

            do m=ms+1,me-1
               if (idir.eq.1) then
                  s(m,n,i)=f(m-ms)
               elseif (idir.eq.2) then
                  s(n,m,i)=f(m-ms)
               endif
            enddo
         enddo
      enddo

      return
      end

c**********************************************************************
      subroutine lhsinv(q,s,qn,cn,d2p,ms,me,ns,ne,idir,xn,yn,vnu,jd,kd,
     >                                                iblank,tscale,bt)

c**********************************************************************
      use params_global
c**********************************************************************
      implicit none
c**********************************************************************
      integer ms,me,ns,ne,idir,jd,kd
      integer iblank(jd,kd)
      real q(jd,kd,nq),s(jd,kd,nv),tscale(jd,kd),bt(jd,kd)
      real qn(jd,kd),cn(jd,kd),d2p(jd,kd),xn(jd,kd),yn(jd,kd),vnu(jd,kd)

      !local variables
      integer,allocatable :: iblnk(:) 
      real,allocatable :: jcbn(:),tscal(:)
      real,allocatable :: vistrm1(:),vistrm2(:),vistrm3(:),g(:)
      real,allocatable :: a(:),b(:),c(:),d(:),e(:),f(:)
      real,allocatable :: a1(:),b1(:),c1(:),d1(:),e1(:)
      real,allocatable :: diag(:,:)

      integer i,m,n
      real svt,dis2,dis4,X,bSq
      real c2,c2m,c4m2,c4m,c4,c4p

      allocate(iblnk(me))
      allocate(jcbn(me),tscal(me))
      allocate(vistrm1(me),vistrm2(me),vistrm3(me),g(me))
      allocate(a(me),b(me),c(me),d(me),e(me),f(me))
      allocate(a1(me),b1(me),c1(me),d1(me),e1(me))
      allocate(diag(me,nq))

c..some initialization

      dis2 = 10.0
      dis4 = 0.1

      do n = ns,ne
         if (idir.eq.1) then
            do m = ms,me
               jcbn(m)  = q(m,n,nq)
               tscal(m) = tscale(m,n)
               g(m)     = d2p(m,n)
               iblnk(m) = max(iblank(m,n),0)
               bSq = Mp**2/(bt(m,n)-Mp**2*(bt(m,n)-1))
               X= sqrt(((1.-bSq)*qn(m,n))*((1.-bSq)*qn(m,n))
     c            +4.*bSq*cn(m,n)*cn(m,n))
               diag(m,1) = qn(m,n)
               diag(m,2) = qn(m,n)
               diag(m,3) = 0.5*((bSq+1.)*qn(m,n)+X)
               diag(m,4) = 0.5*((bSq+1.)*qn(m,n)-X)
c              spectral radius divided by jacobian
               diag(m,nq)= 0.5*((bSq+1.)*abs(qn(m,n))+X)/q(m,n,nq)
            enddo

            do m = ms+1,me-1
               vistrm1(m) = (xn(m,n)*xn(m,n)+yn(m,n)*yn(m,n))/q(m,n,nq)
               vistrm1(m) = vistrm1(m) + (xn(m+1,n)*xn(m+1,n)
     c                                 +yn(m+1,n)*yn(m+1,n))/q(m+1,n,nq)
               vistrm1(m) = vistrm1(m)*q(m,n,nq)*vnu(m,n)

               vistrm3(m) = (xn(m,n)*xn(m,n)+yn(m,n)*yn(m,n))/q(m,n,nq)
               vistrm3(m) = vistrm3(m) + (xn(m-1,n)*xn(m-1,n)
     c                                 +yn(m-1,n)*yn(m-1,n))/q(m-1,n,nq)
               vistrm3(m)= vistrm3(m)*q(m,n,nq)*vnu(m,n)

               vistrm2(m) = 0.5*(vistrm1(m) + vistrm3(m))
            enddo
            m = ms
            vistrm1(m) = 0
            vistrm2(m) = 0
            vistrm3(m) = 0
            m = me
            vistrm1(m) = 0
            vistrm2(m) = 0
            vistrm3(m) = 0

         elseif (idir.eq.2) then

            do m = ms,me
               jcbn(m) = q(n,m,nq)
               tscal(m) = tscale(n,m)
               g(m)= d2p(n,m)
               iblnk(m) = max(iblank(n,m),0)
               bSq = Mp**2/(bt(n,m)-Mp**2*(bt(n,m)-1))
               X= sqrt(((1.-bSq)*qn(n,m))*((1.-bSq)*qn(n,m))
     c            +4.*bSq*cn(n,m)*cn(n,m))
               diag(m,1)= qn(n,m)
               diag(m,2)= qn(n,m)
               diag(m,3)= 0.5*((bSq+1.)*qn(n,m)+X)
               diag(m,4)= 0.5*((bSq+1.)*qn(n,m)-X)
               diag(m,nq)= 0.5*((bSq+1.)*abs(qn(n,m))+X)/q(n,m,nq)
            enddo

            do m = ms+1,me-1
               vistrm1(m) = (xn(n,m)*xn(n,m)+yn(n,m)*yn(n,m))/q(n,m,nq)
               vistrm1(m) = vistrm1(m) + (xn(n,m+1)*xn(n,m+1)
     c                                 +yn(n,m+1)*yn(n,m+1))/q(n,m+1,nq)
               vistrm1(m) = vistrm1(m)*q(n,m,nq)*vnu(n,m)
         
               vistrm3(m) = (xn(n,m)*xn(n,m)+yn(n,m)*yn(n,m))/q(n,m,nq)
               vistrm3(m) = vistrm3(m) + (xn(n,m-1)*xn(n,m-1)
     c                                 +yn(n,m-1)*yn(n,m-1))/q(n,m-1,nq)
               vistrm3(m) = vistrm3(m)*q(n,m,nq)*vnu(n,m)

               vistrm2(m) = 0.5*(vistrm1(m)+vistrm3(m))
            enddo
            m = ms
            vistrm1(m) = 0
            vistrm2(m) = 0
            vistrm3(m) = 0
            m = me
            vistrm1(m) = 0
            vistrm2(m) = 0
            vistrm3(m) = 0
 
         endif

 
         do m = ms+3,me-3
c           second order dissipation
            c2   = 0.5*((g(m)+g(m+1)))*dis2
            c2m  = 0.5*((g(m)+g(m-1)))*dis2

            c4p  = (dis4-min(dis4,.5*dis2*(g(m+2)+g(m+1))))
            c4   = (dis4-min(dis4,1.*c2))
            c4m  = (dis4-min(dis4,1.*c2m))
            c4m2 = (dis4-min(dis4,.5*dis2*(g(m-2)+g(m-1))))

            a1(m) = (diag(m+1,nq)+diag(m+2,nq))*c4p
            b1(m) = -((diag(m+2,nq)+diag(m+1,nq))*c4p+
     c               (diag(m+1,nq)+diag(m,nq))*(3.*c4+c2))
            c1(m) = +((diag(m-1,nq)+diag(m,nq))*(3.*c4m+c2m)+
     c               (diag(m+1,nq)+diag(m,nq))*(3.*c4+c2))
            d1(m) = -((diag(m-2,nq)+diag(m-1,nq))*c4m2+
     c               (diag(m,nq)+diag(m-1,nq))*(3.*c4m+c2m))
            e1(m) = (diag(m-1,nq)+diag(m-2,nq))*c4m2
         enddo 

         m    = ms+2
         c2   = 0.5*((g(m)+g(m+1)))*dis2
         c2m  = 0.5*((g(m)+g(m-1)))*dis2
         c4p  = (dis4-min(dis4,.5*dis2*(g(m+2)+g(m+1))))
         c4   = (dis4-1.*min(dis4,c2))
         c4m  = (dis4-min(dis4,1.*c2m))
         c4m2 = (dis4-min(dis4,1.*dis2*(g(m-1))))
         a1(m) = (diag(m+1,nq)+diag(m+2,nq))*c4p
         b1(m) = -((diag(m+2,nq)+diag(m+1,nq))*c4p+
     c            (diag(m+1,nq)+diag(m,nq))*(3.*c4+c2))
         c1(m) = ((diag(m-1,nq)+diag(m,nq))*(3.*c4m+c2m)+
     c           (diag(m+1,nq)+diag(m,nq))*(3.*c4+c2))
         d1(m) = -((2.*diag(m-1,nq))*c4m2+
     c            (diag(m,nq)+diag(m-1,nq))*(3.*c4m+c2m))
         
         m    = ms+1
         c2   = 0.5*((g(m)+g(m+1)))*dis2
         c2m  = ((g(m)))*dis2
         c4p  = (dis4-min(dis4,.5*dis2*(g(m+2)+g(m+1))))
         c4   = (dis4-min(dis4,1.*c2))
         c4m  = (dis4-min(dis4,1.*c2m))
         a1(m) = (diag(m+1,nq)+diag(m+2,nq))*c4p
         b1(m) = -((diag(m+2,nq)+diag(m+1,nq))*c4p+
     c            (diag(m+1,nq)+diag(m,nq))*(3.*c4+c2))
         c1(m) = (2.*(diag(m,nq))*(2.*c4m+c2m)+
     c            (diag(m+1,nq)+diag(m,nq))*(3.*c4+c2))
         m    = me-2
         c2   = 0.5*((g(m)+g(m+1)))*dis2
         c2m  = 0.5*((g(m)+g(m-1)))*dis2
         c4p  = (dis4-min(dis4,1.*dis2*(g(m+1))))
         c4   = (dis4-min(dis4,1.*c2))
         c4m  = (dis4-min(dis4,1.*c2m))
         c4m2 = (dis4-min(dis4,.5*dis2*(g(m-2)+g(m-1))))
         b1(m) = -(2.*(diag(m+1,nq))*c4p+
     c            (diag(m+1,nq)+diag(m,nq))*(3.*c4+c2))
         c1(m) = ((diag(m-1,nq)+diag(m,nq))*(3.*c4m+c2m)+
     c           (diag(m+1,nq)+diag(m,nq))*(2.*c4+c2))
         d1(m) = -((diag(m-2,nq)+diag(m-1,nq))*c4m2+
     c            (diag(m,nq)+diag(m-1,nq))*(3.*c4m+c2m))
         e1(m) = (diag(m-1,nq)+diag(m-2,nq))*c4m2
         
         m    = me-1  
         c2   = ((g(m)))*dis2
         c2m  = 0.5*((g(m)+g(m-1)))*dis2
         c4   = (dis4-min(dis4,1.*c2))
         c4m  = (dis4-min(dis4,1.*c2m))
         c4m2 = (dis4-min(dis4,.5*dis2*(g(m-2)+g(m-1))))
         c1(m) = +((diag(m-1,nq)+diag(m,nq))*(3.*c4m+c2m)+
     c            2.*(diag(m,nq))*(2.*c4+c2))
         d1(m) = -((diag(m-2,nq)+diag(m-1,nq))*c4m2+
     c            (diag(m,nq)+diag(m-1,nq))*(3.*c4m+c2m))
         e1(m) = (diag(m-1,nq)+diag(m-2,nq))*c4m2

         do m = ms,me
            a1(m) = a1(m)*jcbn(m)
            b1(m) = b1(m)*jcbn(m)
            c1(m) = c1(m)*jcbn(m)
            d1(m) = d1(m)*jcbn(m)
            e1(m) = e1(m)*jcbn(m)
c         viscous terms
            b1(m) = b1(m) - 0.5*vistrm1(m)
            c1(m) = c1(m) + vistrm2(m)
            d1(m) = d1(m) - 0.5*vistrm3(m)
         enddo

         do i = 1,nv
            do m = ms,me
               a(m) = a1(m)
               b(m) = b1(m)
               c(m) = c1(m)
               d(m) = d1(m)
               e(m) = e1(m)
               if (idir.eq.1) then
                  f(m)    = s(m,n,i)
               elseif (idir.eq.2) then
                  f(m)    = s(n,m,i)
               endif

c           euler terms

               b(m) = b(m) - diag(m,i)*0.5
               d(m) = d(m) + diag(m,i)*0.5

c         multiply by delta t
               if(m.le.me-2) then
			a(m) = a(m)*tscal(m+2)
               else
			a(m) = a(m)*tscal(m)
	       endif
               if(m.le.me-1) then
			b(m) = b(m)*tscal(m+1)
               else
			b(m) = b(m)*tscal(m)
	       endif
               if(m.ge.ms+1) then
			d(m) = d(m)*tscal(m-1)
               else
			d(m) = d(m)*tscal(m)
	       endif
               if(m.ge.ms+2) then
			e(m) = e(m)*tscal(m-2)
               else
			e(m) = e(m)*tscal(m)
	       endif

c          add identity to matrix
               c(m) = c(m)*tscal(m) + 1.0
            enddo

            do m = ms+1,me-1
               a(m-ms) = a(m) 
               b(m-ms) = b(m) 
               c(m-ms) = c(m) 
               d(m-ms) = d(m) 
               e(m-ms) = e(m) 
               f(m-ms) = f(m) 
            enddo
            d(1) = 0
            e(1) = 0
            e(2) = 0
            b(me-ms-1) = 0
            a(me-ms-1) = 0
            a(me-ms-2) = 0

            call pentadag(a,b,c,d,e,f,me-ms-1)

            do m=ms+1,me-1
               if (idir.eq.1) then
                  s(m,n,i)=f(m-ms)
               elseif (idir.eq.2) then
                  s(n,m,i)=f(m-ms)
               endif
            enddo
         enddo
      enddo 

      return
      end

c**********************************************************************
      subroutine pentadag(a,b,c,d,e,f,N)

c**********************************************************************
      use params_global
c**********************************************************************
      implicit none
c**********************************************************************
      integer N
      real a(N),b(N),c(N),d(N),e(N),f(N)
 
      !local variables
      integer i
      real d0,d1,d2

c     making diagonal 1.
      i = 1
      d0 = 1.0/c(i)
      d(i+1) = d(i+1)*d0
      e(i+2) = e(i+2)*d0
      f(i)   = f(i)*d0
      
c     making diagonal 1. and lower diagonal 0.
      i = 2
      d1 = b(i-1)
      d0 = 1.0/(c(i)-d1*d(i))
      f(i)   = (f(i)-d1*f(i-1))*d0
      d(i+1) = (d(i+1)-d1*e(i+1))*d0
      e(i+2) = e(i+2)*d0
      
c     making diaongal 1. and lower diagonals 0.
      do i = 3,N
         d1 = a(i-2)
         d2 = (b(i-1) - a(i-2)*d(i-1))
         d0 = 1./(c(i)-d2*d(i)-d1*e(i))
         f(i) = (f(i)-d2*f(i-1)-d1*f(i-2))*d0
         if (i.le.N-1) then
            d(i+1) = (d(i+1)-d2*e(i+1))*d0
            if (i.le.n-2) then
               e(i+2) = (e(i+2))*d0
            endif
         endif
      enddo 

c     backward sweep
      i  = N-1
      d1 = 0
      d2 = 0
      do i = N-1,1,-1
         d1 = d(i+1)
         if (i.lt.N-1) then
            d2 = e(i+2)
C..asitav
            f(i) = f(i) - d1*f(i+1) - d2*f(i+2)
         else
            f(i) = f(i) - d1*f(i+1) 
         endif
         !f(i) = f(i) - d1*f(i+1) - d2*f(i+2)
      enddo

      return
      end
c**********************************************************************
      subroutine read_inputs(jgmx,kgmx,nmesh)

c**********************************************************************
      use params_global
c**********************************************************************
      implicit none
c**********************************************************************
      integer nmesh
      integer jgmx(nmesh),kgmx(nmesh)

      ! local variables
      real a,alphamean,cs,ss,tu
      integer istop,i,i1,ii,im,ihar
      logical file_exists

C..Parameter for boundary conditions
      integer :: nbc,ibtyp(25),ibdir(25)
      integer :: jbcs(25),kbcs(25)
      integer :: jbce(25),kbce(25)
      
      namelist/bcinp/ nbc,ibtyp,ibdir,jbcs,jbce,kbcs,kbce

      namelist/flapinp/ pxc,dela,rf_tef,
     &                 nharmflap,ampFlap,phiFlap,theta_f0,theta_finit

      namelist/vdleinp/ pxc_vdle,dela_vdle,rf_vdle,
     &                 nharmvdle,ampvdle,phivdle,theta_vdle0,theta_vdleinit

      namelist/slatinp/ theta_slat0,slatx0,slaty0,slatdx0,slatdy0,rf_slat,
     &                 slatxc,slatyc,nharmslat,ampSlat,phiSlat,
     &                 slatdx0_nh,slatdy0_nh

c***  first executable statement

c..read in the input and write it back out

      read (5, inputs)
      write (6, inputs)

      inquire(file="bc.inp", exist=file_exists)
      
      if (file_exists) then
        open(unit=10,file='bc.inp',status='unknown')
        do im = 1,num_grids
          read(10,bcinp)
          write(6,bcinp)
          nbc_all(im) = nbc
          ibtyp_all(:,im) = ibtyp(:)
          ibdir_all(:,im) = ibdir(:)
          jbcs_all(:,im) = jbcs(:)
          jbce_all(:,im) = jbce(:)
          kbcs_all(:,im) = kbcs(:)
          kbce_all(:,im) = kbce(:)
        enddo
        close(10)
      else
        write(6,*) 'USING DEFAULT BC'
        write(6,*) '    - FIRST MESH --> C-type Airfoil Mesh'
        write(6,*) '    - SECOND MESH --> Rectangular Wind-tunnel Mesh '
      endif 

      pi=4*atan(1.0)

c..read in the TEF inputs and write it back out
      if (iteflap.eq.1) then
        rewind(5)
        read(5,flapinp)
        write (6, flapinp)

        dela=dela*pi/180
        theta_f0 = -theta_f0*pi/180.0
        theta_finit = -theta_finit*pi/180.0
        do ihar=1,nharmflap
           ampFlap(ihar)=-ampFlap(ihar)*pi/180.
           phiFlap(ihar)= phiFlap(ihar)*pi/180.
        enddo
      endif

c..read in the VDLE inputs and write it back out
      if (ivdle.eq.1) then
        rewind(5)
        read(5,vdleinp)
        write (6, vdleinp)

        dela_vdle=dela_vdle*pi/180
        theta_vdle0 = -theta_vdle0*pi/180.0
        theta_vdleinit = -theta_vdleinit*pi/180.0
        do ihar=1,nharmvdle
           ampvdle(ihar)=-ampvdle(ihar)*pi/180.
           phivdle(ihar)= phivdle(ihar)*pi/180.
        enddo
      endif
c
      if (islat.eq.1) then
        rewind(5)
        read(5,slatinp)
        write (6, slatinp)

        theta_slat0 = -theta_slat0*pi/180.0
        do ihar=1,nharmSlat
           ampSlat(ihar)=-ampSlat(ihar)*pi/180.
           phiSlat(ihar)=phiSlat(ihar)*pi/180.
        enddo

      endif

c..calculate a few things based on inputs

      if (invisc) then
        rmue = 0.0
        lamin = .true.
      endif
      if (lamin) iturb = -1
      if (half.eq.0) then
        jle = (jmax+1)/2
        jtail2 = jmax - jtail1 + 1
      else
        jle = jmax - 1
        jtail2 = jmax -  1
      endif

c..write summary of inputs (and check validity!)

      istop = 0
      write(6,*) ' '
      write(6,*) 'here is your input:'

      write(6,*) 'restart'
      if(iread.eq.0) then
        write(6,*) '  initial run w/o restart for this case '
      elseif(iread.eq.1) then
        write(6,*) '  restart run for this case '
      else
        write(6,*) '**invalid iread, must be 0 or 1**'
        istop = istop+1
      endif

      write(6,*) '  this run will last until ',nsteps,' time steps'
      write(6,*) 'flow parameters'
      write(6,*) '  free stream mach number is ',fsmach
      if(alfa.ne.0.) write(6,*) '  flow at ',
     <               alfa,' degrees angle of attack'
      if(invisc) then
        write(6,*) '  inviscid flow'
      else
        write(6,*) '  viscous flow'
        write(6,*) '    reynolds number = ',rey
        if(lamin) then
          write(6,*) '      laminar flow'
        else
          write(6,*) '      turbulent flow'
          if (iturb.eq.0) then
            write(6,*) '     Baldwin-Lomax turbulence model'
          elseif (iturb.eq.1) then
            write(6,*) '     Spalart-Allmaras turbulence model'
          elseif (iturb.eq.2) then
            write(6,*) '     k-omega SST turbulence model'
            if (itrans.eq.1) then
              write(6,*) '       Langter-Mentry transition model'
            endif
          endif
        endif
      endif

      write(6,*) 'time info'
      if(iunst.eq.0) then
        write(6,*) '  steady flow '
      elseif(iunst.eq.1) then
        write(6,*) '  unsteady flow (pitching oscillation)'
      elseif(iunst.eq.2) then
        write(6,*) '  unsteady flow (pitching ramp)'
      elseif(iunst.eq.3) then
        write(6,*) '  unsteady flow (prescribed pitching)' 
      elseif(iunst.eq.4) then
        write(6,*) '  unsteady flow (pitching oscillation about 1/4c)'
      elseif(iunst.eq.5) then
        write(6,*) '  unsteady flow (pitching oscillation about 1/4c
     c               with deformation)'
      else
        write(6,*) '**invalid iunst, must be between 0 and 5**'
        istop = istop+1
      endif
      if(ntac.eq.1) then
        write(6,*) '  1st order in time '
      elseif(ntac.eq.2) then
        write(6,*) '  2nd order in time '
      elseif(ntac.eq.3) then
        write(6,*) '  3rd order in time '
      else
        write(6,*) '**invalid ntac, must be 1 or 2 or 3**'
        istop = istop+1
      endif
      if(itnmax.eq.1) then
        write(6,*) '  no use of newton iterations '
      elseif(itnmax.gt.1) then
        write(6,*) '  using newton iterations with ',
     <        itnmax,' iterations'
      else
        write(6,*) '**invalid itnmax, must be 0 or greater**'
        istop = istop+1
      endif
      if(iunst.gt.0 .and. dt.gt.0.10) then
        write(6,*) '**invalid dt for unsteady flow**'
        istop=istop+1
      elseif(dt.gt.0.) then
        write(6,*) '  dt is ',dt
      else
        write(6,*) '**invalid dt, must be greater than 0'
        istop=istop+1
      endif
      if(timeac.ne.1.0 .and. iunst.gt.0) then
        write(6,*) '**invalid timeac for unsteady flow**'
        istop=istop+1
      elseif(timeac.le.1.0 .and. timeac.ge.0.0) then
        write(6,*) '  timeac for jacobian scaling is ',timeac
      else
        write(6,*) '**invalid timeac, must be between 0 and 1**'
        istop=istop+1
      endif
c     
      write(6,*) 'algorithm stuff'
      if(epse.gt. 0.0 .and. epse.lt. 1.0) then
        write(6,*) '  dissipation for ilu3d is ',epse
      else
        write(6,*) '**invalid epse, must be between 0 and 1**'
        istop=istop+1
      endif
      if(iprecon) then
        write(6,*) '  Preconditioning will be used '
      endif
      if(iallmach) then
        write(6,*) '  All mach correction will be used '
      endif
      if(iprecon.AND.iallmach) then
        write(6,*) '**Both all mach correction and preconditioning cannot be used together**'
        istop=istop+1
      endif
      if(irhsy.eq.-1) then
        write(6,*) '  1st order in space '
      elseif(irhsy.eq.-2) then
        write(6,*) '  2nd order in space '
      elseif(irhsy.eq.-3) then
        write(6,*) '  3rd order in space '
      else
        write(6,*) '**invalid irhsy, must be between -1 and -3**'
        istop=istop+1
      endif
      if(ilim.ge.1) then
        write(6,*) '  limiting in both directions'
      elseif(ilim.le.-1) then
        write(6,*) '  limiting in only the j-direction'
      elseif(ilim.eq.0) then
        write(6,*) '  no limiting'
      endif
      if(abs(ilim).eq.0) then
        write(6,*) '     using muscl scheme'
      elseif(abs(ilim).eq.1) then
        write(6,*) '     using muscl with korens differentiable limiter'
      elseif(abs(ilim).eq.2) then
        write(6,*) '     using muscl with van a. differentiable limiter'
      elseif(abs(ilim).eq.3) then
        write(6,*) '     using muscl with c-o minmod scheme'            
      elseif(abs(ilim).eq.4) then
        write(6,*) '     using sonic-a scheme of hyunh et al.'          
      elseif(abs(ilim).eq.5) then
        write(6,*) '     using sonic extension of c-o minmod scheme'    
      elseif(abs(ilim).eq.7) then
        write(6,*) '     using cubic interpolation with no limiting'    
      elseif(abs(ilim).eq.8) then
        write(6,*) '     using quartic interpolation with no limiting'  
      elseif(abs(ilim).eq.9) then
        write(6,*) '     using quadratic reconstruction (pade) '        
      endif
      if(iread.eq.0) then
        write(6,*) '  totime = ',totime
      else
        write(6,*) '  totime = ',totime,' will be reset from q-file'
      endif
c
      if(iunst.eq.1) then
        write(6,*) 'reduced frequency is ',rf
        write(6,*) 'amplitude of pitching oscillation is (deg) ',angmax
      endif
c     
      if(iunst.eq.2) then
        write(6,*) 'pitch ramp time is ',rf
        write(6,*) 'amplitude of pitch ramp is (deg)',angmax
      endif
c     
      if(jint+kint.gt.2) then
        write(6,*) 'coarsening of grid for this run, use with caution'
        if(jint.ne.1)
     <     write(6,*) ' use every ',jint,' points in j-direction'
        if(kint.ne.1)
     <     write(6,*) ' use every ',kint,' points in k-direction'
      endif
c
      if(istop.gt.0) write(6,*) 'note: ',istop,' errors in input******'
      write(6,*) ' '
c
c...setting pseudo-time step cycle
c
      ii = itnmax + 1
      do i = 1,itnmax
        if (dtpseudo(i).eq.0) then
          ii = i
          exit 
        endif
      enddo

!      if (ii.eq.0) goto 107

      if (ii.eq.1) then
        dtpseudo = 1.
      else
        do i = ii,itnmax
          i1 = mod(i,ii-1)
          if (i1.eq.0) i1 = ii-1
          dtpseudo(i) = dtpseudo(i1)
        enddo
      endif

! 107  continue

c...setting freestream conditions

      h  = dt
      hd = .5*dt
c
      gm1   = gamma -1.0
      ggm1  = gamma*gm1
      cs    = cos( pi*alfa/180.0 )
      ss    = sin( pi*alfa/180.0 )
      theta_col=alfa
      alf   = alfa
c
      einf  = 1.0/ggm1 +0.5*fsmach**2
      pinf  = 1.0/gamma
      rinf  = 1.0
      ainf  = gamma*pinf/rinf
      htinf = gamma*einf/rinf - 0.5*gm1*fsmach**2

      nq = 5
      nv = 4
      nmv = 4

		! SA with transition model
      if (iturb.eq.1.and.itrans.eq.1) then
        nq = nq + 2
        nv = nv + 2
		  vnuinf = 1.e-6
      endif

		! k-omega model
      if (iturb.eq.2) then
        nq = 7
        nv = 6
        if (itrans.eq.1) then
          nq = nq + 2
          nv = nv + 2
        endif
      endif
c
      uinf  = fsmach*cs
      vinf  = fsmach*ss
c
      rey   = rey/fsmach

Casitav 
C For fully-turbulent k-omega model inflow conditions, refer to
C http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/20070035069_2007034647.pdf 
C (Effective Inﬂow Conditions for Turbulence Models in Aerodynamic Calculations,
C  AIAA Journal, Vol. 45, No. 10, 2007, pp. 2544-2553.)
C
C  last working: a) for af, tu = 0.2, eddy viscosity ratio mu_t/mu = 2
C
		if(iturb==2) then

        if (itrans.eq.0) then

			 ! Fully turbulent k-omega model

			 !Recommended values of High-Re applications (full-scale)
          tkeinf = 1e-6*fsmach**2 
          tomegainf = 5*fsmach/rey 

        else

			 ! k-omega w/ transition model
			 ! Based on user specified tuinf & vnuinf
          tkeinf = 1.5*(0.01*tuinf*fsmach)**2
			 tomegainf = rinf*tkeinf/vnuinf

			 ! On-going transition model calibration process
          !tkeinf = 1e-5*fsmach**2 
          !tomegainf = fsmach/rey 

        endif

			 print*,'tuinf, vnuinf = ',tuinf, vnuinf
			 print*,'tkein, omegainf = ',tkeinf, tomegainf

		endif

      if(dt.lt.0.0 .and. iunst.ge.1) dt = abs(dt)*pi/abs(rf)/180.

	 nmesh=num_grids

      if (iread.eq.0) then
         if(num_grids.gt.1) read(1) nmesh
         read(1) (jgmx(i),kgmx(i),i=1,nmesh)
			print*,'grid dimensions = ',jgmx(1),kgmx(1)
      else

         if(num_grids.gt.1) read(3) nmesh
         read(3) (jgmx(i),kgmx(i),i=1,nmesh)
         
         if(num_grids.gt.1) read(4) nmesh
         read(4) (jgmx(i),kgmx(i),i=1,nmesh)

         !for basegrid
         if(num_grids.gt.1) read(1) nmesh
         read(1) (jgmx(i),kgmx(i),i=1,nmesh)
         
      endif

      return
      end


c**********************************************************************
      !subroutine initia(q,s,x,y,xx,xy,yx,yy,ug,vg,turmu,vnut,
      subroutine initia(q,s,x,y,xg,yg,xt2,yt2,xx,xy,yx,yy,ug,vg,turmu,vnut,
     &     yx0,yy0,yt0,xbig,ybig,xold,yold,xole,yole,iblank,im,jd,kd,cl)

c
c  initialize the variables to default values, read inputs, check them
c  and set up the working files.
c
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************
      integer jd,kd,im
      
      real q(jd,kd,nq),s(jd,kd,nv),vnut(jd,kd)
      real x(jd,kd),y(jd,kd)
C asitav
      real xg(jd,kd),yg(jd,kd)
      real xt2(jd,kd),yt2(jd,kd)
C
      real xx(jd,kd),xy(jd,kd),yx(jd,kd),yy(jd,kd)
      real ug(jd,kd),vg(jd,kd)
      real yy0(jd),yx0(jd),yt0(jd)
      real turmu(jd,kd),cl
      integer iblank(jd,kd)
      real xbig(2*jd,2*kd),ybig(2*jd,2*kd)
      real xold(2*jd,2*kd),yold(2*jd,2*kd)
      real xole(2*jd,2*kd),yole(2*jd,2*kd)

      ! local variables

      integer j,k,ib

c**   first executable statement


      pi = 4.0*atan(1.)
      dang = 0.
      cl=0.0
      if(O_GRID.EQ.1) then 
	IBCWP=0
	print*,'OGRID - therefore IBCWP= 0'
      endif
      do 882 j = 1, jd
         do 882 k=1,kd
            turmu(j,k) = 0.
            ug(j,k) = 0.
            vg(j,k) = 0.
            iblank(j,k)=1
882   continue
      
      do 884 j=1,jd
         do 884 k=1,kd
            s(j,k,1) = 0.
            s(j,k,2) = 0.
            s(j,k,3) = 0.
            s(j,k,4) = 0.
884   continue

c
      jm    = jmax - 1
      km    = kmax - 1

c..setup q and x

      istep0  = 0

	print*,'NUM GRIDS=========== ',num_grids

      if (iread.eq.0) then
         call grid(x,y,iblank,jd,kd,1)
         xg = x; yg = y;
      else
         call grid(x,y,iblank,jd,kd,4)
         call grid(xg,yg,iblank,jd,kd,1)
      endif
      
      bodyflag(im) = .FALSE.
      jtail(im) = 1 
      do ib=1,nbc_all(im)
        if (ibtyp_all(ib,im).eq.5) then
          bodyflag(im) = .TRUE.
          jtail(im) = jbcs_all(ib,im)
        endif
      enddo

      jtail1 = jtail(im)
C
      if (half.eq.0) then
        jle = (jmax+1)/2
        jtail2 = jmax - jtail1 + 1
      else
        jle = jmax - 1
        jtail2 = jmax -  1
      endif

		if(flatplate) then
			print*,'enter jtail1 as leading edge j-index of flat plate'
			read(6,*)jle_flatplate
			jtail(im) = jle_flatplate
			jtail2 = jmax
		endif
C
      call qzero( q,jd,kd)

      if(iread.eq.0.or.invisc) then
c         call lamvis(q,vnut,jd,kd)
         do j=1,jmax
            do k=1,kmax
               vnut(j,k)=1.0*vnuinf
            enddo
         enddo
      endif

      if( iread .gt. 0 ) then
        call restr2( q,jd,kd,3 )
      endif

      !initial flap deflection
      if(iteflap.eq.1 .and. abs(theta_finit).lt.30.) then
       if( num_grids.eq.1 .or.  (num_grids.eq.2 .and. (.not.bg) .and. im.eq.2)
     >     .or. (num_grids.eq.2 .and. bg .and. im.eq.1) 
     >     .or. (num_grids.eq.3 .and. im.eq.2))
     >  call rigid_flap(xg,yg,jd,kd,.true.)
      end if

!a
      !initial vdle deflection
      if(ivdle.eq.1 .and. abs(theta_finit).lt.30.) then
       if( num_grids.eq.1 .or.  (num_grids.eq.2 .and. (.not.bg) .and. im.eq.2)
     >     .or. (num_grids.eq.2 .and. bg .and. im.eq.1)
     >     .or. (num_grids.eq.3 .and. im.eq.2))
     >  call rigid_vdle(xg,yg,jd,kd,.true.)
      end if
!a
      if(theta_init.lt.30.0) then
      if( (im.lt.num_grids).or.(im.eq.num_grids.and.(.not.bg))) then
             call pitch(x,y,xg,yg,xt2,yt2,xx,xy,ug,yx,yy,vg,
     >             yx0,yy0,yt0,jd,kd,im,.true.)
      end if
      end if
c     call move_new(1,x,y,xx,xy,ug,yx,yy,vg,yx0,yy0,yt0,jd,kd)

c..compute the metrics and the jacobians using finite volume formula

      call metfv(q,x,y,xx,xy,yx,yy,xbig,ybig,xold,yold,xole,yole,jd,kd)

c..divide q by the jacobian, q/q5 

      call qdivj( q,jd,kd )
 
c..calculate dt if cnbr given or else calculate cnbrmax

      call eigen(q,xx,xy,yx,yy,ug,vg,jd,kd )

      return
      end

c***********************************************************************
      subroutine metfv( q,x,y,xx,xy,yx,yy,xbig,ybig,xold,yold,
     &     xole,yole,jd,kd)
c  finite volume formulation         11/4/87  s.o.
c  compute the metrics and the jacobian for computational space
c  uniform computational space, deltas = 1 < averaged metrics >
c  modified to incorporate refined mesh 1996
c
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************
      integer jd,kd
 
      real q(jd,kd,nq), x(jd,kd),y(jd,kd)
      real xx(jd,kd), xy(jd,kd), yx(jd,kd), yy(jd,kd)
      real xbig(2*jd,2*kd),ybig(2*jd,2*kd)
      real xold(2*jd,2*kd),yold(2*jd,2*kd)
      real xole(2*jd,2*kd),yole(2*jd,2*kd)

      ! local variables

      real,allocatable :: a3(:,:)
      integer,allocatable :: jjp(:),jjr(:),kkp(:),kkr(:)
      real,allocatable :: fline(:),fdotline(:),fintline(:)
      real fac1,fac2,dx2,dy2,qj
      integer j,k,nneg

      allocate(a3(jd,kd))
      allocate(jjp(jd),jjr(jd),kkp(kd),kkr(kd))
      allocate(fline(mdim),fdotline(mdim),fintline(mdim))

c***  first executable statement

      fac1    = 0.5
      fac2    = 0.5

      do 1  j = 1,jmax
       jjp(j) = j + 1
       jjr(j) = j - 1
    1 continue
      jjp(jmax)=jmax
      jjr(1   )=1
      do 2  k = 1,kmax
        kkp(k) = k + 1
        kkr(k) = k - 1
    2 continue
      kkp(kmax)=kmax
      kkr(1   )=1
c..store the old refined mesh for time metrics
      do k=1,kmax*2
      do j=1,jmax*2
        xole(j,k) = xold(j,k)
        yole(j,k) = yold(j,k)
        xold(j,k) = xbig(j,k)
        yold(j,k) = ybig(j,k)
      enddo
      enddo
c..let's define the refined mesh for metric calculations
      do k=1,kmax
      do j=1,jmax
        xbig(2*j-1,2*k-1)=x(j,k)
        ybig(2*j-1,2*k-1)=y(j,k)
      enddo
      enddo
c
      do k=1,kmax
        do j=1,jmax
          fline(j) = y(j,k)
        enddo
        call metslope(fline,fdotline,fintline,1,jmax)
        do j=1,jmax-1
          ybig(2*j,2*k-1) = fintline(j)
        enddo
      enddo
      do k=1,kmax
        do j=1,jmax
          fline(j) = x(j,k)
        enddo
        call metslope(fline,fdotline,fintline,1,jmax)
        do j=1,jmax-1
          xbig(2*j,2*k-1) = fintline(j)
        enddo
      enddo
c
      do j=1,2*jmax-1
        do k=1,kmax
          fline(k) = ybig(j,2*k-1)
        enddo
        call metslope(fline,fdotline,fintline,1,kmax)
        do k=1,kmax-1
          ybig(j,2*k) = fintline(k)
        enddo
      enddo
c
      do j=1,2*jmax-1
        do k=1,kmax
          fline(k) = xbig(j,2*k-1)
        enddo
        call metslope(fline,fdotline,fintline,1,kmax)
        do k=1,kmax-1
          xbig(j,2*k) = fintline(k)
        enddo
      enddo
c..now find the volume near j,k
      do k=1,kmax
      do j=1,jmax
        a3(j,k) = 0.0
      enddo
      enddo
      do k=1,kmax-1
      do j=1,jmax-1
          dx1 = xbig(j*2,k*2)-xbig(j*2-1,k*2-1)
          dy1 = ybig(j*2,k*2)-ybig(j*2-1,k*2-1)
          dx2 = xbig(j*2-1,k*2)-xbig(j*2,k*2-1)
          dy2 = ybig(j*2-1,k*2)-ybig(j*2,k*2-1)
          a3(j,k) = a3(j,k) + 0.5*( dx1*dy2 -dx2*dy1 )
      enddo
      enddo
      do k=2,kmax
      do j=1,jmax-1
          dx1 = xbig(j*2,k*2-1)-xbig(j*2-1,k*2-2)
          dy1 = ybig(j*2,k*2-1)-ybig(j*2-1,k*2-2)
          dx2 = xbig(j*2-1,k*2-1)-xbig(j*2,k*2-2)
          dy2 = ybig(j*2-1,k*2-1)-ybig(j*2,k*2-2)
          a3(j,k) = a3(j,k) + 0.5*( dx1*dy2 -dx2*dy1 )
      enddo
      enddo
      do k=1,kmax-1
      do j=2,jmax
          dx1 = xbig(j*2-1,k*2)-xbig(j*2-2,k*2-1)
          dy1 = ybig(j*2-1,k*2)-ybig(j*2-2,k*2-1)
          dx2 = xbig(j*2-2,k*2)-xbig(j*2-1,k*2-1)
          dy2 = ybig(j*2-2,k*2)-ybig(j*2-1,k*2-1)
          a3(j,k) = a3(j,k) + 0.5*( dx1*dy2 -dx2*dy1 )
      enddo
      enddo
      do k=2,kmax
      do j=2,jmax
          dx1 = xbig(j*2-1,k*2-1)-xbig(j*2-2,k*2-2)
          dy1 = ybig(j*2-1,k*2-1)-ybig(j*2-2,k*2-2)
          dx2 = xbig(j*2-2,k*2-1)-xbig(j*2-1,k*2-2)
          dy2 = ybig(j*2-2,k*2-1)-ybig(j*2-1,k*2-2)
          a3(j,k) = a3(j,k) + 0.5*( dx1*dy2 -dx2*dy1 )
      enddo
      enddo
      do k=1,kmax
      do j=1,jmax
        q(j,k,nq) = 1./a3(j,k)
      enddo
      enddo
      do k=1,kmax,kmax-1
      do j=1,jmax
        q(j,k,nq) = 0.5*q(j,k,nq)
      enddo
      enddo
      do j=1,jmax,jmax-1
      do k=1,kmax
        q(j,k,nq) = 0.5*q(j,k,nq)
      enddo
      enddo
c
c..let's actually use the refined mesh in finite volume like way
c  rather than use slopes form metslope
c
c..xi derivatives
c
      do j=1,jmax
        k=1
        xx(j,k) = (ybig(2*j-1,2*k)-ybig(2*j-1,2*k-1))*2.
        xy(j,k) =-(xbig(2*j-1,2*k)-xbig(2*j-1,2*k-1))*2.
        do k=2,km
          xx(j,k) = (ybig(2*j-1,2*k)-ybig(2*j-1,2*k-2))
          xy(j,k) =-(xbig(2*j-1,2*k)-xbig(2*j-1,2*k-2))
        enddo
        k=kmax
        xx(j,k) = (ybig(2*j-1,2*k-1)-ybig(2*j-1,2*k-2))*2.
        xy(j,k) =-(xbig(2*j-1,2*k-1)-xbig(2*j-1,2*k-2))*2.
      enddo
c
c..eta derivatives
c
      do k=1,kmax
        j=1
        yx(j,k) =-(ybig(2*j,2*k-1)-ybig(2*j-1,2*k-1))*2.
        yy(j,k) = (xbig(2*j,2*k-1)-xbig(2*j-1,2*k-1))*2.
        do j=2,jm
          yx(j,k) =-(ybig(2*j,2*k-1)-ybig(2*j-2,2*k-1))
          yy(j,k) = (xbig(2*j,2*k-1)-xbig(2*j-2,2*k-1))
        enddo
        j=jmax
        yx(j,k) =-(ybig(2*j-1,2*k-1)-ybig(2*j-2,2*k-1))*2.
        yy(j,k) = (xbig(2*j-1,2*k-1)-xbig(2*j-2,2*k-1))*2.
      enddo
c
      do 44 k = 1,kmax
       do 45 j  = 1,jmax
       qj     = q(j,k,nq)
       xx(j,k)= qj*(xx(j,k))
       xy(j,k)= qj*(xy(j,k))
       yx(j,k)= qj*(yx(j,k))
       yy(j,k)= qj*(yy(j,k))
   45  continue
   44 continue
c
c..check for negative jacobians
c
      nneg = 0
        do 910 k = 1, kmax
          do 920 j = 1, jmax
ccray            nneg = cvmgt(nneg+1,nneg,q(j,k,nq).le.0.)
            if( q(j,k,nq).le.0.0 ) then
              nneg = nneg+1
            end if
 920      continue
 910    continue
c
      if(nneg .ne. 0) then
        write(6,*) nneg, ' negative jacobians in block'
        do 74 k = 1,kmax
        do 74 j = 1,jmax
          if( q(j,k,nq).lt.0.0 ) then
            write(6,603) q(j,k,nq), j, k
          end if
   74   continue
      endif
c
  603 format( ' ',10x,'negative jacobian = ',1p,e10.3,1x,'at j,k =',
     $                 2i5,5x)
c
      return
      end


*************************************************************************
      subroutine metslope(f,fdot,fint,is,ie)
c
c m3-quartic scheme of hyunh et al to determine the slope
c                                               and midpoints
c*************************************************************************
      use params_global
c*************************************************************************
      implicit none

      integer is,ie

      real f(mdim),fdot(mdim),fint(mdim)
      real f0(mdim),f1(mdim),f2(mdim),f2m(mdim),quar(mdim)

      ! local variables

      real ammd,fl,fr
      integer i
      real at,ati1,s1,t1,sl,tmax
c*************************************************************************
      ammd(fl,fr) = 0.5*(sign(1.,fl)+sign(1.,fr))*amin1(abs(fl),abs(fr))
c*************************************************************************

c***  first executable statement

c..the m3-quartic interpolation scheme of hyunh follows


c..let's load up a few difference arrays
c
        do i=is,ie
          f0(i) = f(i)
        enddo
c..1st difference at i+1/2
        do i=is,ie-1
          f1(i) = f0(i+1) - f0(i)
        enddo
c..2nd difference at i
        do i=is+1,ie-1
          f2(i)  = f1(i) - f1(i-1)
        enddo
c..extrapolate at the boundaries
        f2(is) = 2.*f2(is+1)-f2(is+2)
        f2(ie) = 2.*f2(ie-1)-f2(ie-2)
c..limit 2nd difference to i+1/2
        do i = is,ie-1
          f2m(i) = ammd(f2(i),f2(i+1))
        enddo
c..quartic slope
      do i=is+2,ie-2
        quar(i) = (+f(i-2)-8.*f(i-1)+8.*f(i+1)-f(i+2))/12.
      enddo
      i=is
      quar(i) = (-22.*f(i)+36.*f(i+1)-18.*f(i+2)+4.*f(i+3))/12.
      i=is+1
      quar(i) = (-2.*f(i-1)-3.*f(i)+6.*f(i+1)-f(i+2))/6.
      i=ie-1
      quar(i) = (2.*f(i+1)+3.*f(i)-6.*f(i-1)+f(i-2))/6.
      i=ie
      quar(i) = (22.*f(i)-36.*f(i-1)+18.*f(i-2)-4.*f(i-3))/12.
c
c..now combine everything to get the slope
c
      do i=is+1,ie-1
c..include limited curvatures to calculate new slopes
        at    = f1(i)   - 0.5*f2m(i)
        ati1  = f1(i-1) + 0.5*f2m(i-1)
c..limit slopes at i
        s1    = ammd(f1(i),f1(i-1))
        t1    = ammd(at,ati1)
c..now find appropriate slope
c       sl = 0.5*(at+ati1)
        sl = quar(i)
        sl = sl+ammd(at-sl,ati1-sl)
        tmax  = sign(1.,t1)*amax1(3.*abs(s1),1.5*abs(t1))
        sl = ammd(sl,tmax)
c..use slope to calculate ql and qr
        fdot(i) = sl
      enddo
c..first-order at boundary?
c
      sl = quar(is)
      sl = ammd(sl,3.*f1(is))
      fdot(is) = sl
c
      sl = quar(ie)
      sl = ammd(sl,3.*f1(ie-1))
      fdot(ie) = sl
c
      do i=is,ie-1
        fint(i) = f(i)+fdot(i)*.5+(3.*f1(i)-2.*fdot(i)-fdot(i+1))*.25
     <           +(fdot(i)+fdot(i+1)-2.*f1(i))*.125
      enddo
c
      return
      end

c***********************************************************************
      subroutine monitor(mstop,q,jd,kd,im)
c
c  this subroutine checks for negative speed of sound
c
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************
      integer jd,kd,im
      real q(jd,kd,nq)

      ! local variables
      integer mstop,nneg,k,j
      real asqmin,rho,qq,asq

c*** first executable statement

c
c..negative speed of sound check
c
      mstop = 0
      asqmin = 0.002
      nneg = 0
       do 910 k = 1, kmax
        do 920 j = 1, jmax
          rho = q(j,k,1)
          qq = q(j,k,2)**2 +q(j,k,3)**2
          asq = ggm1*(q(j,k,4) - .5*qq/rho)/rho
ccray          nneg = cvmgt(nneg+1,nneg,asq.le.asqmin)
          if( asq.le.asqmin ) then
            nneg = nneg+1
          end if
 920    continue
 910   continue
c
      if(nneg .ne. 0) then
        write(6,*) nneg, ' negative speed of sound in mesh: ',im
        do 74 k = 1,kmax
        do 74 j = 1,jmax
          rho = q(j,k,1)
          qq = q(j,k,2)**2 +q(j,k,3)**2
          asq = ggm1*(q(j,k,4) - .5*qq/rho)/rho
          if( asq.le.asqmin ) then
            mstop = 1
            write(6,601) j,k,asq
            return
          end if
   74   continue
      end if
c
 601  format(' j,k,asq from monitor ',2i5,f12.5)
c
      return
      end

c***********************************************************************
      subroutine move(init,x,y,xx,xy,ug,yx,yy,vg,yx0,yy0,yt0,jd,kd)
c
c  this stores the old surface metrics as well as rotates the grid and 
c  updates the metrics.  
c
c  if iunst=1
c      rf = reduced frequency
c      maxang = amplitude of pitching oscillation
c         pos/neg about origin
c  if iunst=2
c      rf = time for pitching ramp
c      maxang = change in pitching angle
c  if iunst=3
c      read in angle for each time step
c  if iunst=4
c      rf = reduced frequency
c      maxang = amplitude of pitching oscillation
c         pos only about orgin
c  if iunst=5
c      interpolation between grids, unsteady but not simply transformation
c      allowed here
c      all metrics get recalculated after adapt subroutine call
c      just keep this so iunst is not zero, steady state
c  if iunst=6
c      step change in grid velocity (modify, grid vel metric, keep
c                                    grid postion metric const)
c      phasestart=total time for pitching motion	
c

c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************
      integer jd,kd
      real xx(jd,kd), xy(jd,kd), yx(jd,kd), yy(jd,kd)
      real ug(jd,kd), vg(jd,kd)
      real x(jd,kd), y(jd,kd)
      real yx0(jd), yy0(jd), yt0(jd)

      ! local variables
      real,allocatable :: xtmpm(:,:),ytmpm(:,:)
      real thnm1,vfact,angnew,angold,xac,yac,dummy,t
      real scalnew,scalold,sn,cn,sno,cno,xold,yold,xnm1,ynm1
      integer j,k,init
      
      allocate(xtmpm(jd,kd),ytmpm(jd,kd))

      thnm1 = 0.0

c for grid interpolation
      if(iunst.eq.5) then
        thetao = 0.0
        thetan = 0.0
      endif

c for steady problems
      if(iunst.eq.0) then
        thetao = 0.0
        thetan = 0.0
      endif

c for ramps
      if(iunst.eq.3) then
        if(init.eq.1) then
          thnm1 = 0.0
          thetao = 0.0
          read(17,*) thetan     
          vfact  = 0.0
        elseif(init.eq.2) then
          thnm1 = thetao
          thetao = thetan 
          read(17,*) thetan
          vfact = 1.0
        else
          thnm1 = thetao
          thetao = thetan
          read(17,*) thetan
          vfact = 1.0
        endif
      endif
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c this is oscillations about origin both positive and negative
      if (iunst.eq.1) then
        angnew = totime*rf
        if(init.eq.1) then
          angold = 0.
          vfact=0.
        else
          angold = (totime-dt)*rf
          vfact=1.
        endif
        thetan = sin(angnew)*(pi/180.)*angmax
        thetao = sin(angold)*(pi/180.)*angmax
      endif
c
c this is oscillations purely positive from origin (le of airfoil, not 1/4c)
      if (iunst.eq.4) then
        angnew = totime*rf
        if(init.eq.1) then
          angold = 0.
          vfact=0.
        else
          angold = (totime-dt)*rf
          vfact=1.
        endif
        thetan = sin(angnew-pi/2.)*(pi/180.)*angmax/2.
     c                                      +angmax/2.*(pi/180.)
        thetao = sin(angold-pi/2.)*(pi/180.)*angmax/2.
     c                                      +angmax/2.*(pi/180.)
        xac=.25
        yac=0.0
      endif
      if (iunst.eq.6) then
        vfact=1.
        thetao = 0.0
        thetan = 0.0
	if (istep.eq.0) then
          read(200) dummy,dummy
          read(200)((xtmpm(j,k),j=1,jmax),k=1,kmax),
     <             ((ytmpm(j,k),j=1,jmax),k=1,kmax)
          close(200)
	endif
      endif

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      if (iunst.eq.2) then
        t = totime/rf
        if(t.gt.1.) t=1.
        if(t.lt.0.) t=0.
        scalnew = (10.-15.*t+6.*t*t)*t**3
        if(init.eq.1) then
          scalold = 0.
          vfact=0.
        else
          t = (totime-dt)/rf
          if(t.gt.1.) t=1.
          if(t.lt.0.) t=0.
          scalold = (10.-15.*t+6.*t*t)*t**3
          vfact=1.
        endif
        thetan = scalnew*(pi/180.)*angmax
        thetao = scalold*(pi/180.)*angmax
      endif
******************************************************************
******************************************************************
******************************************************************
c now calculating metrics of transformation for grid
c
      dang = thetan*180./pi
c
      sn = sin(thetan-thetao)
      cn = cos(thetan-thetao)
      sno = sin(thnm1-thetao)
      cno = cos(thnm1-thetao)
c..store old metrics at surface
      do 450 j = 1,jmax
        yx0(j) = yx(j,1)
        yy0(j) = yy(j,1)
        yt0(j) = -ug(j,1)*yx(j,1)-vg(j,1)*yy(j,1)
  450 continue
c
      do 81 k = 1,kmax
      do 81 j = 1,jmax
c..grid
        xold  = x(j,k)
        yold  = y(j,k)
        if (iunst.ne.4) then    !for rotation about origin
          x(j,k)  = +sn*yold+cn*xold
          y(j,k)  = -sn*xold+cn*yold
        else                    !for rotation about xac
          x(j,k)  = (+sn*(yold-yac)+cn*(xold-xac))+xac
          y(j,k)  = (-sn*(xold-xac)+cn*(yold-yac))+yac
        endif
        xnm1    = +sno*yold+cno*xold
        ynm1    = -sno*xold+cno*yold
c..grid velocities
        ug(j,k) = (x(j,k)-xold)/dt*vfact
        vg(j,k) = (y(j,k)-yold)/dt*vfact
        if(iunst.eq.3 .and. ntac.ge.2 .and. istep.gt.1) then
          ug(j,k) = (1.5*x(j,k)-2.*xold+0.5*xnm1)/dt*vfact
          vg(j,k) = (1.5*y(j,k)-2.*yold+0.5*ynm1)/dt*vfact
        endif
	if (iunst.eq.6) then		!for instataneous change in pitch rate
c	  ug(j,k) = (xtmpm(j,k)-xold)/phasestart*vfact
c         vg(j,k) = (ytmpm(j,k)-yold)/phasestart*vfact
c         ug(j,k) = (xold-xtmpm(j,k))/phasestart*vfact
c         vg(j,k) = (yold-ytmpm(j,k))/phasestart*vfact
c         ug(j,k) = (xtmpm(j,k)-xold)/dt*vfact
c         vg(j,k) = (ytmpm(j,k)-yold)/dt*vfact
	if ((k.lt.60).and.((j.gt.58).and.(j.lt.160))) then
          ug(j,k) = -(xtmpm(j,k)-xold)/.125*vfact*fsmach
          vg(j,k) = -(ytmpm(j,k)-yold)/.125*vfact*fsmach
	else
	  ug(j,k) = 0.
	  vg(j,k) = 0.
	endif
	

	endif
c..xi metrics       
        xold  = xx(j,k)
        yold  = xy(j,k)
        xx(j,k)  = +sn*yold+cn*xold
        xy(j,k)  = -sn*xold+cn*yold
c..eta metrics   
        xold  = yx(j,k)
        yold  = yy(j,k)
        yx(j,k)  = +sn*yold+cn*xold
        yy(j,k)  = -sn*xold+cn*yold
 81   continue
c
      return
      end


c***********************************************************************
      subroutine move_new(init,x,y,xx,xy,ug,yx,yy,vg,yx0,yy0,yt0,jd,kd)
c
c  this stores the old surface metrics as well as rotates the grid and 
c  updates the metrics.  
c
c  if iunst=1
c      interpolation between grids, unsteady but not simply transformation
c      allowed here
c      all metrics get recalculated after adapt subroutine call
c      just keep this so iunst is not zero, steady state

c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************
      integer jd,kd,init
      real xx(jd,kd), xy(jd,kd), yx(jd,kd), yy(jd,kd)
      real ug(jd,kd), vg(jd,kd)
      real x(jd,kd), y(jd,kd)
      real yx0(jd), yy0(jd), yt0(jd)
      ! local variables

      real xac,yac,dtheta,sn,cn,snold,cnold,xnew,ynew,ang,xold,yold
      integer j,k
c***********************************

!     jaina, changing this subroutine, it looked too rough
!     can do better for ug(j,k),vg(j,k) with analytical expressions

      pi = 4.0*atan(1.0)
c     
c     xac=0.25
c     yac=0.0228
      xac=0.25
      yac=0.0
c	print*,'xac, yac= ',xac,yac
      
      dtheta=angmax*(sin(rf*istep*dt)-sin(rf*(istep-1)*dt))
      !dtheta=angmax*dt*rf*cos(rf*totime)
      theta_col=theta_col+dtheta
      sn = sin(dtheta*pi/180.)
      cn = cos(dtheta*pi/180.)
      snold = sin(-1.*dtheta*pi/180.)
      cnold = cos(-1.*dtheta*pi/180.)

c..store old metrics at surface

      do 450 j = 1,jmax
        yx0(j) = yx(j,1)
        yy0(j) = yy(j,1)
        yt0(j) = -ug(j,1)*yx(j,1)-vg(j,1)*yy(j,1)
  450 continue
c
      do 81 k = 1,kmax
      do 81 j = 1,jmax
c..grid

        xnew=(x(j,k)-xac)*cn+(y(j,k)-yac)*sn+xac
        ynew=-(x(j,k)-xac)*sn+(y(j,k)-yac)*cn+yac

	if(ntac.eq.1) then
 	ug(j,k)=(xnew-x(j,k))/dt
   	vg(j,k)=(ynew-y(j,k))/dt
	else
        xold=(x(j,k)-xac)*cnold+(y(j,k)-yac)*snold+xac
        yold=-(x(j,k)-xac)*snold+(y(j,k)-yac)*cnold+yac
 	ug(j,k)=(1.5*xnew-2*x(j,k)+0.5*xold)/dt
 	vg(j,k)=(1.5*ynew-2*y(j,k)+0.5*yold)/dt
	endif

        x(j,k)=xnew
        y(j,k)=ynew

  81	continue

	ang=atan2(y(jle,1)-y(jtail1,1),x(jtail1,1)-x(jle,1))

c	print *,'ang=',ang*180/pi



      return
      end

c***********************************************************************
      subroutine movie(q,x,y,iblank,ug,vg,jd,kd,init,logg,logq)
c
c  this writes data out to file by planes
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************
      integer jd,kd,init,logg,logq
      real x(jd,kd), y(jd,kd), ug(jd,kd), vg(jd,kd)
      integer iblank(jd,kd)
      real q(jd,kd,nq)
      
      ! local variables
      integer j,k
      real ppp(jd,kd)
c***********************************************************************

      if(init.eq.1) then
        write(logg) jmax,kmax,int((nsteps - istep0-isin)/nmovie)
        write(logq) jmax,kmax,int((nsteps - istep0-isin)/nmovie)
        write(logq) fsmach,alf,rey,totime
      else
        do 1 k = 1, kmax
        do 1 j = 1, jmax
          ppp(j,k) = gm1*( q(j,k,4)-0.5*(q(j,k,2)**2+q(j,k,3)**2)
     <                                    /q(j,k,1))*q(j,k,nq)
    1   continue


	
        write(logg) ((x(j,k),j=1,jmax),k=1,kmax),
     <            ((y(j,k),j=1,jmax),k=1,kmax),
     <            ((0.0000001*float(istep0/nmovie),j=1,jmax),k=1,kmax),
     <	          ((iblank(j,k),j=1,jmax),k=1,kmax)	


        write(logq) ((q(j,k,1)*q(j,k,nq),j=1,jmax),k=1,kmax),
     <   ((q(j,k,2)/q(j,k,1),j=1,jmax),k=1,kmax),
     <   ((q(j,k,3)/q(j,k,1),j=1,jmax),k=1,kmax),
     <   ((ppp(j,k)/(q(j,k,1)*q(j,k,nq))**1.4,j=1,jmax),k=1,kmax),
     <   (((ppp(j,k)-1./1.4)/(0.5*fsmach*fsmach),j=1,jmax),k=1,kmax)
      endif

      return
      end

c************************************************************************
      subroutine iflux(f,ql,qr,is,ie,im,th,qt,eps,fmin,fmax,
     <                  ibmin,ibmax,iba)
c
c 2nd order weno scheme with Van-Leer limiter
c*************************************************************************
      use params_global
c*************************************************************************
      implicit none
c*************************************************************************

      integer is,ie,im,ibmin,ibmax,iba(mdim)
      real th,qt,eps
      real f(mdim,nmv),fmin(nmv),fmax(nmv)
      real ql(mdim,nmv),qr(mdim,nmv)

      ! local variables
      integer i,n
      real s1,s2,slope,lim_vanleer

c..this is just 1st order upwind

      if(qt.eq.0.0)then
        do n=1,nmv
          do i=is,ie
            ql(i,n)=f(i,n)
            qr(i,n)=f(i,n)
          enddo
        enddo
      else
        do n=1,nmv
          do i = is+1,ie-1
            s1     = f(i,n)-f(i-1,n)
            s2     = f(i+1,n)-f(i,n)
            !slope  = lim_vanleer(s1,s2)
            slope  = (s1+s2)/2
            ql(i,n)= f(i,n) + 0.5*slope
            qr(i,n)= f(i,n) - 0.5*slope
          enddo

          if(ibmin.eq.2) then
            s1   = f(is,n)-fmin(n)
            s2   = f(is+1,n)-f(is,n)
            slope  = lim_vanleer(s1,s2)
            ql(is,n)= f(is,n) + 0.5*slope
            qr(is,n)= f(is,n) - 0.5*slope
          else
            ql(is,n)= f(is,n)
            qr(is,n)= f(is,n)
          endif

          if(ibmax.eq.2) then
            s1   = f(ie,n)-f(ie-1,n)
            s2   = fmax(n)-f(ie,n)
            slope  = lim_vanleer(s1,s2)
            ql(ie,n)= f(ie,n) + 0.5*slope
            qr(ie,n)= f(ie,n) - 0.5*slope
          else
            ql(ie,n)= f(ie,n)
            qr(ie,n)= f(ie,n)
          endif
        enddo
      endif

      return
      end

c*************************************************************************
      real function lim_minmod(x,y)
c*************************************************************************
      real :: x,y

      if(sign(1.0,x).ne.sign(1.0,y)) then
        lim_minmod = 0.0
      elseif(abs(x).lt.abs(y)) then
        lim_minmod = x
      else
        lim_minmod = y
      endif

      end function lim_minmod

c*************************************************************************
      real function lim_vanleer(x,y)
c*************************************************************************
      real :: x,y
      real :: a,b
      real :: lim_minmod
      a = 0.5*(x+y)
      b = 2.0*lim_minmod(x,y)

      lim_vanleer = lim_minmod(a,b)

      end function lim_vanleer

c*************************************************************************
      subroutine muscld(f,ql,qr,
     <                     is,imax,im,th,qt,eps,fmin,fmax,ibmin,ibmax)
c
c  muscl interpolation for higher order accuracy
c  differentiable limiter for 3rd order accuracy
c
c*************************************************************************
      use params_global
c*************************************************************************
      implicit none
c*************************************************************************

      real f(mdim,nmv)
      real ql(mdim,nmv),qr(mdim,nmv)
      real fmin(nmv),fmax(nmv)
      ! local variables
      real,allocatable :: f2(:,:)
      integer is,imax,im,ibmin,ibmax,n,i
      real th,qt,eps,thm,thp,f2i,f2i1,a1,a2,epsf,f3i,f3qt

      allocate(f2(mdim,nmv))
c***  first executable statement

      if(qt.eq.0.0)then
        do 20 n=1,nmv
        do 20 i=is,imax
          ql(i,n)=f(i,n)
          qr(i,n)=f(i,n)
   20   continue
        return
      else
        thm = 1.0-th
        thp = 1.0+th
        do 1 n=1,nmv
          do 11 i=is,im
            f2(i,n) = f(i+1,n) - f(i,n)
 11       continue
          do 12 i=is+1,im
           f2i    = f2(i,n)
           f2i1   = f2(i-1,n)
           a1     = 3.0*f2i*f2i1
           a2     = 2.0*(f2i-f2i1)**2 + a1
           epsf   = eps*( 0.5+sign( 0.5,-abs(a2) ) )
           f3i    = ( a1 +epsf )/( a2 +epsf )
           f3qt   = qt*f3i
           ql(i,n)= f(i,n)+f3qt*( thm*f2i1 + thp*f2i )
           qr(i,n)= f(i,n)-f3qt*( thp*f2i1 + thm*f2i )
 12     continue
c..first-order at boundary
           qr(is  ,n)= f(is,n)
           ql(is  ,n)= f(is,n)
c..central at boundary?
c         qr(is,n) = 0.5*(f(is,n)+f(is+1,n))
c         ql(is,n) = 0.5*(f(is,n)+f(is+1,n))
          if(ibmin.eq.2) then
            f2i    = f(is+1,n)-f(is,n)
            f2i1   = f(is,n)-fmin(n)
            a1     = 3.0*f2i*f2i1
            a2     = 2.0*(f2i-f2i1)**2 + a1
            epsf   = eps*( 0.5+sign( 0.5,-abs(a2) ) )
c	    print*,'fix the epsilon - look at turns3d'
            f3i    = ( a1 +epsf )/( a2 +epsf )
            f3qt   = qt*f3i
            ql(is,n)= f(is,n)+f3qt*( thm*f2i1 + thp*f2i )
            qr(is,n)= f(is,n)-f3qt*( thp*f2i1 + thm*f2i )
          endif
c..first-order at boundary
           ql(imax,n)= f(imax,n)
           qr(imax,n)= f(imax,n)
c..central at boundary?
c         ql(imax,n) = 0.5*(f(imax,n)+f(im,n))
c         qr(imax,n) = 0.5*(f(imax,n)+f(im,n))
          if(ibmax.eq.2) then
            f2i    = fmax(n)-f(imax,n)
            f2i1   = f(imax,n)-f(im,n)
            a1     = 3.0*f2i*f2i1
            a2     = 2.0*(f2i-f2i1)**2 + a1
            epsf   = eps*( 0.5+sign( 0.5,-abs(a2) ) )
            f3i    = ( a1 +epsf )/( a2 +epsf )
            f3qt   = qt*f3i
            ql(imax,n)= f(imax,n)+f3qt*( thm*f2i1 + thp*f2i )
            qr(imax,n)= f(imax,n)-f3qt*( thp*f2i1 + thm*f2i )
          endif
    1   continue
      endif
c
      return
      end

c************************************************************************
      subroutine muscld_new(f,ql,qr,is,imax,im,th,qt,eps,
     <     fmin,fmax,ibmin,ibmax,iba)
c
c  muscl interpolation for higher order accuracy
c  differentiable limiter for 3rd order accuracy
c
c*************************************************************************
      use params_global
c*************************************************************************
      implicit none
c*************************************************************************
      real f(mdim,nmv)
      real ql(mdim,nmv),qr(mdim,nmv)
      real fmin(nmv),fmax(nmv)
      integer iba(mdim)
      
      ! local variables
      real,allocatable :: f2(:,:)
      integer is,imax,im,ibmin,ibmax,n,i
      real thm,thp,f2i,f2i1,a1,a2,f3i,f3qt
      real th,qt,eps

      allocate(f2(mdim,nmv))

      if(qt.eq.0.0)then
        do 20 n=1,nmv
        do 20 i=is,imax
          ql(i,n)=f(i,n)
          qr(i,n)=f(i,n)
   20   continue
        return
      else
        thm = 1.0-th
        thp = 1.0+th
        do 1 n=1,nmv
          do 11 i=is,im
            f2(i,n) = f(i+1,n) - f(i,n)
 11       continue
          do 12 i=is+1,im
           f2i    = f2(i,n)
           f2i1   = f2(i-1,n)
           a1     = 3.0*(f2i*f2i1+eps)
           a2     = 2.0*(f2i-f2i1)**2 + a1
c           epsf   = eps*( 0.5+sign( 0.5,-abs(a2) ) )
           f3i    = a1/a2
           f3qt   = qt*f3i
     >           *iba(i)*iba(i+1)*iba(i-1)
           ql(i,n)= f(i,n)+f3qt*( thm*f2i1 + thp*f2i )
           qr(i,n)= f(i,n)-f3qt*( thp*f2i1 + thm*f2i )
 12     continue
c..first-order at boundary
           qr(is  ,n)= f(is,n)
           ql(is  ,n)= f(is,n)
c..central at boundary?
c         qr(is,n) = 0.5*(f(is,n)+f(is+1,n))
c         ql(is,n) = 0.5*(f(is,n)+f(is+1,n))
          if(ibmin.eq.2) then
            f2i    = f(is+1,n)-f(is,n)
            f2i1   = f(is,n)-fmin(n)
            a1     = 3.0*(f2i*f2i1+eps)
            a2     = 2.0*(f2i-f2i1)**2 + a1
c           epsf   = eps*( 0.5+sign( 0.5,-abs(a2) ) )
            f3i    = a1/a2
            f3qt   = qt*f3i
            ql(is,n)= f(is,n)+f3qt*( thm*f2i1 + thp*f2i )
            qr(is,n)= f(is,n)-f3qt*( thp*f2i1 + thm*f2i )
          endif
c..first-order at boundary
           ql(imax,n)= f(imax,n)
           qr(imax,n)= f(imax,n)
c..central at boundary?
c         ql(imax,n) = 0.5*(f(imax,n)+f(im,n))
c         qr(imax,n) = 0.5*(f(imax,n)+f(im,n))
          if(ibmax.eq.2) then
            f2i    = fmax(n)-f(imax,n)
            f2i1   = f(imax,n)-f(im,n)
            a1     = 3.0*(f2i*f2i1+eps)
            a2     = 2.0*(f2i-f2i1)**2 + a1
c           epsf   = eps*( 0.5+sign( 0.5,-abs(a2) ) )
            f3i    = a1/a2
            f3qt   = qt*f3i
            ql(imax,n)= f(imax,n)+f3qt*( thm*f2i1 + thp*f2i )
            qr(imax,n)= f(imax,n)-f3qt*( thp*f2i1 + thm*f2i )
          endif
    1   continue
      endif
c
      return
      end


c************************************************************************
      subroutine weno3(f,ql,qr,is,imax,im,th,qt,eps,
     <     fmin,fmax,ibmin,ibmax,iba)
c
c  3rd order WENO interpolation for higher order accuracy
c
c*************************************************************************
      use params_global
c*************************************************************************
      implicit none
c*************************************************************************
      real f(mdim,nmv)
      real ql(mdim,nmv),qr(mdim,nmv)
      real fmin(nmv),fmax(nmv)
      integer iba(mdim)
      
      ! local variables
      real,allocatable :: f2(:,:)
      integer is,imax,im,ibmin,ibmax,n,i
      real f2i,f2i1,a1,a2,w01,w02,w11,w12,f3
      real th,qt,eps,epsw

      allocate(f2(mdim,nmv))

      epsw = 1e-6

      if(qt.eq.0.0)then
        do 20 n=1,nmv
        do 20 i=is,imax
          ql(i,n)=f(i,n)
          qr(i,n)=f(i,n)
   20   continue
        return
      else
        do 1 n=1,nmv
          do 11 i=is,im
            f2(i,n) = f(i+1,n) - f(i,n)
 11       continue
          do 12 i=is+1,im
           f2i    = f2(i,n)
           f2i1   = f2(i-1,n)
           a1     = 1./(f2i1**2+epsw)
           a2     = 1./(f2i**2+epsw)
           w01     = a1/(a1+2*a2) 
           w02     = 1. - w01 
           w11     = 2*a1/(2*a1+a2) 
           w12     = 1. - w11 
           f3   = 0.5*iba(i)*iba(i+1)*iba(i-1)
           ql(i,n)= f(i,n)+f3*( w01*f2i1 + w02*f2i )
           qr(i,n)= f(i,n)-f3*( w11*f2i1 + w12*f2i )
 12     continue
c..first-order at boundary
           qr(is  ,n)= f(is,n)
           ql(is  ,n)= f(is,n)
c..central at boundary?
          if(ibmin.eq.2) then
            f2i    = f(is+1,n)-f(is,n)
            f2i1   = f(is,n)-fmin(n)
            a1     = 1./(f2i1**2+epsw)
            a2     = 1./(f2i**2+epsw)
            w01     = a1/(a1+2*a2) 
            w02     = 1. - w01 
            w11     = 2*a1/(2*a1+a2) 
            w12     = 1. - w11 
            f3   = 0.5
            ql(is,n)= f(is,n)+f3*( w01*f2i1 + w02*f2i )
            qr(is,n)= f(is,n)-f3*( w11*f2i1 + w12*f2i )
          endif
c..first-order at boundary
           ql(imax,n)= f(imax,n)
           qr(imax,n)= f(imax,n)
c..central at boundary?
          if(ibmax.eq.2) then
            f2i    = fmax(n)-f(imax,n)
            f2i1   = f(imax,n)-f(im,n)
            a1     = 1./(f2i1**2+epsw)
            a2     = 1./(f2i**2+epsw)
            w01     = a1/(a1+2*a2) 
            w02     = 1. - w01 
            w11     = 2*a1/(2*a1+a2) 
            w12     = 1. - w11 
            f3   = 0.5
            ql(imax,n)= f(imax,n)+f3*( w01*f2i1 + w02*f2i )
            qr(imax,n)= f(imax,n)-f3*( w11*f2i1 + w12*f2i )
          endif
    1   continue
      endif
c
      return
      end

c**************************************************************************
      subroutine weno(f,ql,qr,is,ie,im,th,qt,eps,fmin,fmax,
     <                  ibmin,ibmax,iba)
c
c 5th order weno scheme
c note: qr(is) and ql(ie) are never used.
c*************************************************************************
      use params_global
c*************************************************************************
      implicit none
c*************************************************************************

      integer is,ie,im,ibmin,ibmax,iba(mdim)
      real th,qt,eps,at1
      real f(mdim,nmv),fmin(nmv),fmax(nmv)
      real ql(mdim,nmv),qr(mdim,nmv)

      ! local variables
      real,allocatable :: f0(:),f1(:),f2(:),f2m(:)
      real,allocatable :: slope(:)

      real dre,duj,vjp1,vjm1,ejp1,ammd,fl,fr
      integer i,n
      real f0bmin,f1bmin,f2bmin,f0bmax,f1bmax,f2bmax,at,s1,ati1,t1
      real weno5,blankv
c*************************************************************************
      ammd(fl,fr) = 0.5*(sign(1.,fl)+sign(1.,fr))*amin1(abs(fl),abs(fr))
c*************************************************************************

      allocate(f0(mdim),f1(mdim),f2(mdim),f2m(mdim))
      allocate(slope(mdim))

c..this is just 1st order upwind
      if(qt.eq.0.0)then
        do n=1,nmv
          do i=is,ie
            ql(i,n)=f(i,n)
            qr(i,n)=f(i,n)
          enddo
        enddo
        return
      else
c
c 5th order weno scheme follows
c
        do 1 n=1,nmv
c
c..let's load up a few difference arrays
c
          do i=is,ie
            f0(i) = f(i,n)
          enddo
c..1st difference at i+1/2
          do i=is,ie-1
            f1(i) = f0(i+1) - f0(i)
          enddo
c..2nd difference at i
          do i=is+1,ie-1
            f2(i)  = f1(i) - f1(i-1)
          enddo
c..extrapolate at the boundaries
          f2(is) = 2.*f2(is+1)-f2(is+2)
          f2(ie) = 2.*f2(ie-1)-f2(ie-2)
c..modify at boundaries, if needed
          if(ibmin.eq.2) then
            f0bmin = fmin(n)
            f1bmin = f0(is) - f0bmin
            f2(is) = f1(is) - f1bmin
            f2bmin = 2.*f2(is)-f2(is+1)
          endif
          if(ibmax.eq.2) then
            f0bmax = fmax(n)
            f1bmax = f0bmax - f0(ie)
            f2(ie) = f1bmax - f1(ie-1)
            f2bmax = 2.*f2(ie)-f2(ie-1)
          endif
c..limit 2nd difference to i+1/2
          do i = is,ie-1
            f2m(i) = ammd(f2(i),f2(i+1))
          enddo
c
c..now combine everything to get ql and qr
c
c
c..first-order at boundary?
c
          at    = f1(is) - 0.5*f2m(is)
          slope(i) = ammd(at,2.*f1(is))
          ql(is,n) = f0(is) + 0.5*slope(i)
c
          if(ibmin.eq.2) then
            i = is
            s1    = ammd(f1(i),f1bmin)
            at    = f1(i)  - 0.5*f2m(i)
            ati1  = f1bmin + 0.5*ammd(f2(i),f2bmin)
            t1    = ammd(at,ati1)
            slope(i)  =  sign(1.,t1)*
     <         amin1(0.5*abs(at+ati1),amax1(2.*abs(s1),abs(t1)))
            ql(i,n) = f0(i) + 0.5*slope(i)
          endif
c
          at1   = f1(ie-1) + 0.5*f2m(ie-1)
          slope(i) = ammd(at1,2.*f1(ie-1))
          qr(ie,n) = f0(ie) - 0.5*slope(i)
c
          if(ibmax.eq.2) then
            i = ie
            s1    = ammd(f1bmax,f1(i-1))
            at    = f1bmax  - 0.5*ammd(f2(i-1),f2bmax)
            ati1  = f1(i-1) + 0.5*f2m(i-1)
            t1    = ammd(at,ati1)
            slope(i)  =  sign(1.,t1)*
     <         amin1(0.5*abs(at+ati1),amax1(2.*abs(s1),abs(t1)))
            qr(i,n) = f0(i) - 0.5*slope(i)
          endif
c
c..sonic-a near the boundary?
c
          do i=is+1,ie-1,ie-is-2
c..include limited curvatures to calculate new slopes
            at    = f1(i)   - 0.5*f2m(i)
            ati1  = f1(i-1) + 0.5*f2m(i-1)
c..limit slopes at i
            s1    = ammd(f1(i),f1(i-1))
            t1    = ammd(at,ati1)
c..now find appropriate slope
            slope(i)  =  sign(1.,t1)*
     <         amin1(0.5*abs(at+ati1),amax1(2.*abs(s1),abs(t1)))
c..use slope to calculate ql and qr
            ql(i,n) = f0(i) + 0.5*slope(i)
            qr(i,n) = f0(i) - 0.5*slope(i)
          enddo
c
c..suresh at interior
c
          do i=is+2,ie-2
          ql(i,n) = weno5(f0(i-2),f0(i-1),f0(i),f0(i+1),f0(i+2))
          qr(i,n) = weno5(f0(i+2),f0(i+1),f0(i),f0(i-1),f0(i-2))
          blankv  = iba(i)*iba(i-1)*iba(i+1)*iba(i-2)*iba(i+2)
	  ql(i,n) = ql(i,n)*blankv+(1.-blankv)*f0(i)
	  qr(i,n) = qr(i,n)*blankv+(1.-blankv)*f0(i)
          enddo
c
c..sonic-a near the boundary?
c
          i=is+2
c..include limited curvatures to calculate new slopes
            at    = f1(i)   - 0.5*f2m(i)
            ati1  = f1(i-1) + 0.5*f2m(i-1)
c..limit slopes at i
            s1    = ammd(f1(i),f1(i-1))
            t1    = ammd(at,ati1)
c..now find appropriate slope
            slope(i)  =  sign(1.,t1)*
     <         amin1(0.5*abs(at+ati1),amax1(2.*abs(s1),abs(t1)))
c..use slope to calculate ql and qr
            ql(i,n) = f0(i) + 0.5*slope(i)
c
c..sonic-a near the boundary?
c
          i=ie-1
c..include limited curvatures to calculate new slopes
            at    = f1(i)   - 0.5*f2m(i)
            ati1  = f1(i-1) + 0.5*f2m(i-1)
c..limit slopes at i
            s1    = ammd(f1(i),f1(i-1))
            t1    = ammd(at,ati1)
c..now find appropriate slope
            slope(i)  =  sign(1.,t1)*
     <         amin1(0.5*abs(at+ati1),amax1(2.*abs(s1),abs(t1)))
c..use slope to calculate ql and qr
            qr(i,n) = f0(i) - 0.5*slope(i)
    1   continue
      endif
c
      return
      end
c*************************************************************************
      function weno5(a,b,c,d,e)
c
c
c*************************************************************************
      implicit none
c*************************************************************************
      real weno5
      real a,b,c,d,e
      real b1,b2,epsw,djm1,ejm1,dj,ej,djp1,ejp1
      real dis0,dis1,dis2,q30,q31,q32,d01,d02,a1ba0,a2ba0
      real w0,w1,w2
      
      b1 = 13./12.
      b2 = 1./6.
      epsw = 1.e-6
      djm1 = a-2.*b+c
      ejm1 = a-4.*b+3.*c
      dj   = b-2.*c+d
      ej   = b-d
      djp1 = c-2.*d+e
      ejp1 = 3.*c-4.*d+e
      dis0 = b1*djm1*djm1+0.25*ejm1*ejm1+epsw
      dis1 = b1*dj*dj+0.25*ej*ej+epsw
      dis2 = b1*djp1*djp1+0.25*ejp1*ejp1+epsw
      q30 = 2.*a-7.*b+11.*c
      q31 = -b+5.*c+2.*d
      q32 = 2.*c+5.*d-e
      d01 = dis0/dis1
      d02 = dis0/dis2
      a1ba0 = 6.*d01
      a2ba0 = 3.*d02
      w0 = 1./(1.+a1ba0+a2ba0)
      w1 = a1ba0*w0
      w2 = 1.-w0-w1
      weno5 = b2*(w0*q30+w1*q31+w2*q32)
      return
      end

c*************************************************************************
      subroutine quad(f,ql,qr,
     <                    is,ie,im,th,qt,eps,fmin,fmax,ibmin,ibmax)
c
c..piecewise quadratic reconstruction with 6th order compact evaluation of
c  nodal derivatives and new monotonicity-preserving constraint.
c
c**************************************************************************
      use params_global
c**************************************************************************
      implicit none
c**************************************************************************
      integer is,ie,im,ibmin,ibmax
      real th,qt,eps
      real f(mdim,nmv),ql(mdim,nmv),qr(mdim,nmv)
      real fmin(nmv),fmax(nmv)

      ! local variables
      real,allocatable :: d1(:),d2(:)
      real,allocatable :: con1(:),con2(:),con3(:)
      real,allocatable :: s1(:),sf(:),sc(:),sb(:)

      real amedian,ammd
      real b,c,d,fl,fr
      integer n,i,is1,is2,ie2,ie1
      real f1i,f1i1,f1i2,f1i3,f2i1,f2i2,f2i3,f2i,f2mi
      real f1bmin,f2bmin,f2mi1,f2mi2,at,sl,ati1,t1,slopes2
      real f1i4,f1bmax,f2bmax,f2mi3,at1,slopee2,sm,cm,si,ci
      real sr,cr,cl
c**************************************************************************
      ammd(fl,fr) = 0.5*(sign(1.,fl)+sign(1.,fr))*amin1(abs(fl),abs(fr))
      amedian(b,c,d)=b+0.5*(sign(1.,c-b)+sign(1.,d-b))*
     &     amin1(abs(c-b),abs(d-b))
c**************************************************************************

      allocate(d1(mdim),d2(mdim))
      allocate(con1(mdim),con2(mdim),con3(mdim))
      allocate(s1(mdim),sf(mdim),sc(mdim),sb(mdim))
c***  first executable statement

      if(qt.eq.0.)then
        do 10 n=1,nmv
        do 10 i=is,ie
          ql(i,n)=f(i,n)
          qr(i,n)=f(i,n)
   10   continue
        return
      else
        is1=is+1
        is2=is+2
        ie2=ie-2
        ie1=ie-1
c
        do 1 n=1,nmv
c
c..slope:
c
        i=is
        f1i=f(i+1,n)-f(i,n)
        f1i1=f(i+2,n)-f(i+1,n)
        f1i2=f(i+3,n)-f(i+2,n)
        f1i3=f(i+4,n)-f(i+3,n)
        f2i1=f1i1-f1i
        f2i2=f1i2-f1i1
        f2i3=f1i3-f1i2
        f2i=2.*f2i1-f2i2
        if(ibmin.eq.2) then
          f1bmin=f(i,n)-fmin(n)
          f2i=f1i-f1bmin
          f2bmin=2.*f2i-f2i1
        endif
        f2mi=ammd(f2i,f2i1)
        f2mi1=ammd(f2i1,f2i2)
        f2mi2=ammd(f2i2,f2i3)
        at=f1i-0.5*f2mi
        d1(i)=ammd(at,2.*f1i)
        if(ibmin.eq.2) then
          sl=ammd(f1i,f1bmin)
          at=f1i-0.5*f2mi
          ati1=f1bmin+0.5*ammd(f2i,f2bmin)
          t1=ammd(at,ati1)
          d1(i)=sign(1.,t1)*amin1(0.5*abs(at+ati1),
     &          amax1(2.*abs(sl),abs(t1)))
        endif
c
c..sonica scheme at the boundary:
c
        ql(i,n)=f(i,n)+0.5*d1(i)
        qr(i,n)=f(i,n)-0.5*d1(i)
c
        i=is1
        at=f1i1-0.5*f2mi1
        ati1=f1i+0.5*f2mi
        sl=ammd(f1i1,f1i)
        t1=ammd(at,ati1)
        d1(i)=sign(1.,t1)*
     &        amin1(0.5*abs(at+ati1),amax1(2.*abs(sl),abs(t1)))
c
        ql(i,n)=f(i,n)+0.5*d1(i)
        qr(i,n)=f(i,n)-0.5*d1(i)
c
        i=is2
        at=f1i2-0.5*f2mi2
        ati1=f1i1+0.5*f2mi1
        sl=ammd(f1i2,f1i1)
        t1=ammd(at,ati1)
        slopes2=sign(1.,t1)*
     &          amin1(0.5*abs(at+ati1),amax1(2.*abs(sl),abs(t1)))
c
        con1(i)=0.25
        con2(i)=1.
        con3(i)=0.75*(f(i+1,n)-f(i-1,n))
        con3(i)=con3(i)-con1(i)*d1(i-1)
c
        do i=is+3,ie-3
          con1(i)=1./3.
          con2(i)=1.
          con3(i)=(f(i+2,n)+28.*f(i+1,n)-28.*f(i-1,n)-f(i-2,n))/36.
        enddo
c
        i=ie
        f1i1=f(i,n)-f(i-1,n)
        f1i2=f(i-1,n)-f(i-2,n)
        f1i3=f(i-2,n)-f(i-3,n)
        f1i4=f(i-3,n)-f(i-4,n)
        f2i1=f1i1-f1i2
        f2i2=f1i2-f1i3
        f2i3=f1i3-f1i4
        f2i=2.*f2i1-f2i2
        if(ibmax.eq.2) then
          f1bmax=fmax(n)-f(i,n)
          f2i=f1bmax-f1i1
          f2bmax=2.*f2i-f2i1
        endif
        f2mi1=ammd(f2i1,f2i)
        f2mi2=ammd(f2i2,f2i1)
        f2mi3=ammd(f2i3,f2i2)
        at1=f1i1+0.5*f2mi1
        d1(i)=ammd(at1,2.*f1i1)
        if(ibmax.eq.2) then
          sl=ammd(f1bmax,f1i1)
          at=f1bmax-0.5*ammd(f2i1,f2bmax)
          ati1=f1i1+0.5*f2mi1
          t1=ammd(at,ati1)
          d1(i)=sign(1.,t1)*amin1(0.5*abs(at+ati1),
     &          amax1(2.*abs(sl),abs(t1)))
        endif
c
        ql(i,n)=f(i,n)+0.5*d1(i)
        qr(i,n)=f(i,n)-0.5*d1(i)
c
        i=ie1
        at=f1i1-0.5*f2mi1
        ati1=f1i2+0.5*f2mi2
        sl=ammd(f1i1,f1i2)
        t1=ammd(at,ati1)
        d1(i)=sign(1.,t1)*
     &        amin1(0.5*abs(at+ati1),amax1(2.*abs(sl),abs(t1)))
c
        ql(i,n)=f(i,n)+0.5*d1(i)
        qr(i,n)=f(i,n)-0.5*d1(i)
c
        i=ie2
        at=f1i2-0.5*f2mi2
        ati1=f1i3+0.5*f2mi3
        sl=ammd(f1i2,f1i3)
        t1=ammd(at,ati1)
        slopee2=sign(1.,t1)*
     &          amin1(0.5*abs(at+ati1),amax1(2.*abs(sl),abs(t1)))
c
        con1(i)=0.25
        con2(i)=1.
        con3(i)=0.75*(f(i+1,n)-f(i-1,n))
        con3(i)=con3(i)-con1(i)*d1(i+1)
        call tridag(con1,con2,con1,con3,d1,is2,ie2)
c
c..curvature:
c
        i=is
        d2(i)=0.
        i=is1
        d2(i)=0.
        i=is2
        con1(i)=0.1
        con2(i)=1.
        con3(i)=1.2*(f(i+1,n)-2.*f(i,n)+f(i-1,n))
        con3(i)=con3(i)-con1(i)*d2(i-1)
        do i=is+3,ie-3
          con1(i)=2./11.
          con2(i)=1.
          con3(i)=(3.*f(i+2,n)+48.*f(i+1,n)-102.*f(i,n)+48.*f(i-1,n)
     &             +3.*f(i-2,n))/44.
        enddo
        i=ie
        d2(i)=0.
        i=ie1
        d2(i)=0.
        i=ie2
        con1(i)=0.1
        con2(i)=1.
        con3(i)=1.2*(f(i+1,n)-2.*f(i,n)+f(i-1,n))
        con3(i)=con3(i)-con1(i)*d2(i+1)
        call tridag(con1,con2,con1,con3,d2,is2,ie2)
c
c..correct the sign of the slope and curvature:
c
        do i=is1,ie
          s1(i)=f(i,n)-f(i-1,n)
        enddo
        do i=is,ie2
          sf(i)=(-f(i+2,n)+4.*f(i+1,n)-3.*f(i,n))/2.
        enddo
        do i=is1,ie1
          sc(i)=(f(i+1,n)-f(i-1,n))/2.
        enddo
        do i=is2,ie
          sb(i)=(f(i-2,n)-4.*f(i-1,n)+3.*f(i,n))/2.
        enddo
c
        do i=is2,ie2
          sm=amedian(sf(i),sc(i),sb(i))
          if(sm.eq.sf(i)) cm=f(i+2,n)-2.*f(i+1,n)+f(i,n)
          if(sm.eq.sc(i)) cm=f(i+1,n)-2.*f(i,n)+f(i-1,n)
          if(sm.eq.sb(i)) cm=f(i,n)-2.*f(i-1,n)+f(i-2,n)
          d1(i)=amedian(d1(i),sm,sc(i))
          if(d1(i).eq.sm) d2(i)=cm
          if(d1(i).eq.sc(i)) d2(i)=f(i+1,n)-2.*f(i,n)+f(i-1,n)
        enddo
c
        do i=is2,ie2
          si=d1(i)
          ci=d2(i)
c..impose monotonicity-preserving constraints:
c..method 1:
c
c         if(abs(d1(i)).gt.2.0*abs(s1(i+1)).or.abs(d1(i)).gt.
c    &       2.0*abs(s1(i))) then
c           if(abs(d1(i)).le.abs(d1(i+1))) then
c             sr=d1(i)
c             cr=s1(i+1)-d1(i)
c           else
c             sr=2.*s1(i+1)-d1(i+1)
c             cr=-2.*s1(i+1)+2.*d1(i+1)
c           endif
c           if(abs(d1(i)).le.abs(d1(i-1))) then
c             sl=d1(i)
c             cl=-s1(i)+d1(i)
c           else
c             sl=2.*s1(i)-d1(i-1)
c             cl=2.*s1(i)-2.*d1(i-1)
c           endif
c           si=ammd(sl,sr)
c           ci=ammd(cl,cr)
c         endif
c
c..method 2:
c
          if(abs(d1(i)).gt.2.0*abs(s1(i+1)).and.abs(d1(i)).gt.
     &       2.0*abs(s1(i))) then
            if(abs(d1(i)).le.abs(d1(i+1))) then
              sr=d1(i)
              cr=2.*(s1(i+1)-d1(i))
            else
              sr=2.*s1(i+1)-d1(i+1)
              cr=2.*(-s1(i+1)+d1(i+1))
            endif
            if(abs(d1(i)).le.abs(d1(i-1))) then
              sl=d1(i)
              cl=2.*(-s1(i)+d1(i))
            else
              sl=2.*s1(i)-d1(i-1)
              cl=2.*(s1(i)-d1(i-1))
            endif
            si=ammd(sl,sr)
            ci=ammd(cl,cr)
          else
            if(abs(d1(i)).gt.2.0*abs(s1(i+1))) then
              if(abs(d1(i)).le.abs(d1(i+1))) then
                sr=d1(i)
                cr=2.*(s1(i+1)-d1(i))
              else
                sr=2.*s1(i+1)-d1(i+1)
                cr=2.*(-s1(i+1)+d1(i+1))
              endif
              si=ammd(d1(i),sr)
              ci=ammd(d2(i),cr)
            endif
            if(abs(d1(i)).gt.2.0*abs(s1(i))) then
              if(abs(d1(i)).le.abs(d1(i-1))) then
                sl=d1(i)
                cl=2.*(-s1(i)+d1(i))
              else
                sl=2.*s1(i)-d1(i-1)
                cl=2.*(s1(i)-d1(i-1))
              endif
              si=ammd(sl,d1(i))
              ci=ammd(cl,d2(i))
            endif
          endif
c..piecewise quadratic reconstruction:
          ql(i,n)=f(i,n)+0.5*si+ci/12.
          qr(i,n)=f(i,n)-0.5*si+ci/12.
        enddo
c
c..consistency with the sonica scheme:
c
        qr(is2,n)=f(is2,n)-0.5*slopes2
        ql(ie2,n)=f(ie2,n)+0.5*slopee2
 1      continue
      endif
c
      return
      end

c**********************************************************************
      subroutine newton(q,qtn,qtnm1,qnewt,s,ka,kb,jd,kd)
c
c  take into account time terms on rhs for newton iterations
c
c**********************************************************************
      use params_global
c**********************************************************************
      implicit none
c**********************************************************************
      integer ka,kb,jd,kd
      real q(jd,kd,nq),qtn(jd,kd,nv),qtnm1(jd,kd,nv)
      real qnewt(jd,kd,nv),s(jd,kd,nv)

      ! local variables
      integer j,k,n
      real oat,tac

c..set time-accuracy
      oat = 1.0
      if (ntac.eq.-2) oat = 0.5
      if (ntac.ge.2 .and. istep.gt.1) oat = 2./3.
      if (ntac.eq.3 .and. istep.gt.2) oat = 6./11.

      tac = 1./(oat*h) 
      do k = ka, kb
      do j = 2, jm
        do n = 1,nv
          s(j,k,n) = s(j,k,n) - (q(j,k,n) - qnewt(j,k,n))*tac
        enddo
      enddo
      enddo
c     
      return
      end

c**********************************************************************
      subroutine time_step(q,xx,xy,yx,yy,ug,vg,tscale,bt,iblank,jd,kd)
c
c     time step calculation
c**********************************************************************
      use params_global
c**********************************************************************
      implicit none
c**********************************************************************
      integer ka,kb,jd,kd
      real q(jd,kd,nq),xx(jd,kd),xy(jd,kd),yx(jd,kd),yy(jd,kd)
      real ug(jd,kd),vg(jd,kd),tscale(jd,kd),bt(jd,kd)
      integer iblank(jd,kd)

      ! local variables
      integer j,k,jmax_j,kmax_j,jmax_k,kmax_k
      real cflj,cflk,eigj,eigk,hnew
      real u,v,e,uu,vv,dxdx,dydy,asq,lam_j,lam_k
      real cflmaxj,cflmaxk,oat

c..set time-accuracy
      oat = 1.0
      if (ntac.eq.-2) oat = 0.5
      if (ntac.ge.2 .and. istep.gt.1) oat = 2./3.
      if (ntac.eq.3 .and. istep.gt.2) oat = 6./11.

c..if time accurate preconditioning use dual time stepping

      if (iprecon.and.timeac.eq.1) then
        iconstant_cfl = 0
        idual = 1
        idual_turb = 1
      endif

	if(iconstant_cfl.eq.0) then
        if(idual.eq.0) then
	do k=1,kmax
	   do j=1,jmax
	      tscale(j,k)  = oat*h*max(iblank(j,k),0)*( 1.0 + 0.002*(1.-timeac)*
     <                    sqrt(q(j,k,nq)))/( 1.0 + (1.-timeac)*sqrt(q(j,k,nq)))
              bt(j,k) = 1.
	   enddo
        enddo
	else
	do k=1,kmax
	   do j=1,jmax
	      tscale(j,k)   = max(iblank(j,k),0)*( 1.0 + 0.002*sqrt(q(j,k,nq)))
     <                 	      /(1.+sqrt(q(j,k,nq)))
              if (iprecon.and.1.eq.0) then
	        tscale(j,k)   = 100.*max(iblank(j,k),0)*
     <                     ( 1.0 + 0.002*sqrt(q(j,k,nq)))/(1.+sqrt(q(j,k,nq)))
                if (tscale(j,k).gt.100) tscale(j,k) = 100.
              endif
              tscale(j,k)=tscale(j,k)*dualtime
              bt(j,k) = 1./(1.+tscale(j,k)/h/oat)
              tscale(j,k)=tscale(j,k)/(1.+tscale(j,k)/h/oat)
	   enddo
        enddo
	endif


	else
c	print*,'CONSTANT CFL MODE'


	do k=1,kmax
	   do j=1,jmax
 	  u=q(j,k,2)/q(j,k,1)
 	  v=q(j,k,3)/q(j,k,1)
 	  e=q(j,k,4)/q(j,k,1)
 	  uu=(u-ug(j,k))*xx(j,k)+(v-vg(j,k))*xy(j,k)
 	  vv=(u-ug(j,k))*yx(j,k)+(v-vg(j,k))*yy(j,k)
 	  dxdx=xx(j,k)**2 + xy(j,k)**2
 	  dydy=yx(j,k)**2 + yy(j,k)**2
 	  asq=ggm1*(e-0.5*(u**2+v**2))
 	  lam_j=abs(uu)+sqrt(asq*dxdx)
 	  lam_k=abs(vv)+sqrt(asq*dydy)
 	  	hnew=dualtime/(lam_j+lam_k)
	      	tscale(j,k)=hnew*max(iblank(j,k),0)           
                tscale(j,k)=tscale(j,k)/(1.+tscale(j,k)/h/oat)
          bt(j,k) = 1. 
           enddo
	enddo


	endif

  	  if(mod(istep,npnorm).eq.1.and.timeac.eq.1.and.1.eq.0) then
	 do k=1,kmax
	   do j=1,jmax
 	  u=q(j,k,2)/q(j,k,1)
 	  v=q(j,k,3)/q(j,k,1)
 	  e=q(j,k,4)/q(j,k,1)
 	  uu=(u-ug(j,k))*xx(j,k)+(v-vg(j,k))*xy(j,k)
 	  vv=(u-ug(j,k))*yx(j,k)+(v-vg(j,k))*yy(j,k)
 	  dxdx=xx(j,k)**2 + xy(j,k)**2
 	  dydy=yx(j,k)**2 + yy(j,k)**2
 	  asq=ggm1*(e-0.5*(u**2+v**2))
 	  lam_j=abs(uu)+sqrt(asq*dxdx)
 	  lam_k=abs(vv)+sqrt(asq*dydy)
 	  cflj=tscale(j,k)*lam_j
 	  cflk=tscale(j,k)*lam_k

 	  if(j.eq.1.and.k.eq.1) then 
 	  cflmaxj=cflj
 	  cflmaxk=cflk
 	  endif

 	  if(cflj.gt.cflmaxj) then
 	  cflmaxj=cflj
 	  jmax_j=j
 	  kmax_j=k
 	  endif

 	  if(cflk.gt.cflmaxk) then
 	  cflmaxk=cflk
 	  jmax_k=j
 	  kmax_k=k
   	  endif
	   enddo
        enddo

         write(6,103)cflmaxj,jmax_j,kmax_j,cflmaxk,jmax_k,kmax_k 
  103    format('cflmax_j,j,k,cflmaxk,j,k =',2(x,e10.4,i4,i4))
   	  endif



      return
      end


c**********************************************************************
      subroutine pert3(q,x,y,ug,vg,yt0,yx,yy,jd,kd)
c
c..used for field velocity approach
c
c         initialize pert, read inputs,
c         and set up the working files.
c
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************
      integer jd,kd
      real x(jd,kd), y(jd,kd)
      real q(jd,kd,nq)
      real ug(jd,kd), vg(jd,kd)
      real yt0(jd)
      real yx(jd,kd), yy(jd,kd)

      ! local variables
      real smax,sculinv,xa,ya,rjk,wcyl
      integer j,k,iscul,iexp
      
c***********************************
c.  ..  cylindrical velocity distribution.. after scully-s model
c
c     note: smax=gama/(2.*pi*ainf*chord) where gama is the strenth of
c     vortex.  smax is therefore the dimensionless strength of vortex.
c
      smax  = fsmach*vorgam/(2.*pi)
c
c..   load flow variables for vortex in to pq array
c
      do j=1,jmax
        yt0(j) = -ug(j,1)*yx(j,1)-vg(j,1)*yy(j,1)
      enddo
c
      k = 1
      iscul=1
      iexp=2*iscul
      sculinv=1./iscul
      do 30 k = 1,kmax
      do 30 j = 1,jmax
        xa = x(j,k) - xvor 
        ya = y(j,k) - yvor
        rjk = sqrt(xa**2+ya**2)
        wcyl = smax*rjk/(rjk**2+rcore**2)
c        wcyl = smax*rjk/(rjk**iexp+rcore**iexp)**sculinv
        ug(j,k)=ug(j,k)-ya/rjk*wcyl
        vg(j,k)=vg(j,k)+xa/rjk*wcyl
 30   continue
c
      return
      end

c***********************************************************************
      subroutine qdivj( q,jd,kd )
c
c  divide q by jacobian
c
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************
      integer jd,kd
      real q(jd,kd,nq)

      ! local variables
      integer j,k

c***  first executable statement

      do 71 k = 1,kmax
      do 71 j = 1,jmax
        q(j,k,1) = q(j,k,1)/q(j,k,nq)
        q(j,k,2) = q(j,k,2)/q(j,k,nq)
        q(j,k,3) = q(j,k,3)/q(j,k,nq)
        q(j,k,4) = q(j,k,4)/q(j,k,nq)
   71 continue

      return
      end

c***********************************************************************
      subroutine qzero( q,jd,kd )
c
c  set initial values of the flow variables to their free-stream
c  values for an impulsive start. otherwise the initial values
c  are read from a restart file.
c
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************
      integer jd,kd
      real q(jd,kd,nq)

      ! local variables
      real ruinf,rvinf,rtkeinf,rtomegainf,tu,re_theta
      integer j,k

		  real re0,re1,re2,re3,re4,re5,re6
		  real tu0,tu1,tu2,tu3,tu4,tu5,tu6
		  real m0,m1,m2,m3,m4,m5
          real m5

     ! local variables for channel initialization
          real umean(kmax),vmean(kmax),umag,rho1,u1,v1,e1,p1,dum1(kmax,1),dum2(kmax,1) 
          integer write_turb_init, write_init_infile,n


c**   first executable statement
      
      write_turb_init = 0

      if (write_turb_init.eq.1) then
         
         print*,'Initializing with turbulent inflow'
          
 			  p1 = 1.0/gamma
         rho1 = 1.0
     
         open (unit=35, form='formatted', file='udomain.dat') 
         open (unit=36, form='formatted', file='vdomain.dat') 
         do j = 1,jmax
            do k = 1,kmax
               read(35,"(E24.16)") u1
               read(36,"(E24.16)") v1
               q(j,k,1) = rho1
               q(j,k,2) = rho1*u1
               q(j,k,3) = rho1*v1 
               q(j,k,4) = p1/(gamma - 1.0) + 0.5*(u1*u1 + v1*v1)   				
            enddo
         enddo
         close(35)         
         close(36)         

c..         open (unit=35, form='formatted', file='umean_turb.dat') 
c..         do k=1,kmax                   
c..             read(35,"(E24.16)") umean(k)   
c..         enddo     
c..         close(35)         
c..
c..         umag = maxval(abs(umean))
c..         p1 = 1.0/gamma
c..         rho1 = 1.0
c..         vmean = 0.0
c..
c..         do j = 1,jmax   
c..            call myrand(dum1,kmax,1)
c..            call myrand(dum2,kmax,1)
c..            do k = 1,kmax
c..               u1 = umean(k) + dum1(k,1)*0.1*umag
c..               v1 = vmean(k) + dum2(k,1)*0.1*umag
c..               if(k.eq.1) then
c..                  u1 = 0.0
c..                  v1 = 0.0         
c..               endif
c..               q(j,k,1) = rho1
c..               q(j,k,2) = rho1*u1
c..               q(j,k,3) = rho1*v1 
c..               q(j,k,4) = p1/(gamma - 1.0) + 0.5*(u1*u1 + v1*v1)
c..            enddo
c..         enddo
      else
         ruinf = rinf*uinf
         rvinf = rinf*vinf
         do k = 1,kmax
            do j = 1,jmax
               q(j,k,1) = rinf
               q(j,k,2) = ruinf
               q(j,k,3) = rvinf
               q(j,k,4) = einf
            enddo
         enddo
      endif
      
      write_init_infile = 0
      if (write_init_infile.eq.1) then
           open (unit=35, form='unformatted', file='initial_solution.p3d',status='replace') 
           write(35) jmax,kmax
           write(35) fsmach,alf,rey,totime
c..           write(35) 0.0,0.0,0.0,0.0
           write(35) (((q(j,k,n),j=1,jmax),k=1,kmax),n=1,4)       
c..           write(35) ((vnut(j,k),j=1,jd),k=1,kd)
           write(35) ((1.0      ,j=1,jmax),k=1,kmax)
           close(35)
      endif

      if (iturb.eq.2) then
        rtkeinf = rinf*tkeinf
        rtomegainf = rinf*tomegainf
        !tuinf = 100.*sqrt(2./3*tkeinf)/fsmach
        do k=1,kmax
        do j=1,jmax
          q(j,k,5) = rtkeinf
          q(j,k,6) = rtomegainf
        enddo
        enddo
      end if

      if (itrans.eq.1) then
        itmcinf = 1. 

		  ! langtry correlations
!        if(tuinf.le.1.3) then
!          re_theta = (1173.51-589.428*tuinf+0.2196/(tuinf*tuinf))
!        else
!          re_theta = 331.5*(tuinf-0.5658)**(-0.671)
!        endif
!
!        retinf = max(re_theta,20.)

		  ! medida correlations
!		  re1 = 1135.0
!		  re2 = 894.0
!		  re3 = 252.0
!		  re4 = 165.0
!		  re5 = 100.0
!		  tu1 = 0.03
!		  tu2 = 0.51
!		  tu3 = 2.0
!		  tu4 = 5.25
!		  tu5 = 6.5
!		  m1 = (re2-re1)/(tu2-tu1)
!		  m2 = (re3-re2)/(tu3-tu2)
!		  m3 = (re4-re3)/(tu4-tu3)
!		  m4 = (re5-re4)/(tu5-tu4)
!
!			if(tuinf <= tu2) then
!			retinf = re1 + m1*(tuinf-tu1)
!			elseif(tuinf>tu2 .AND. tuinf<=tu3) then
!			retinf = re2 + m2 * (tuinf-tu2)
!			elseif(tuinf>tu3 .AND. tuinf<=tu4) then
!			retinf = re3 + m3 * (tuinf-tu3)
!			elseif(tuinf>tu4 .AND. tuinf<=tu5) then
!			retinf = re4 + m4 * (tuinf-tu4)
!			else
!			retinf = re5
!			endif
		Re0 = 1800.0
		Re1 = 1135.0
		Re2 = 894.0
		Re3 = 392.0
		Re4 = 252.0
		Re5 = 165.0
		Re6 = 100.0

		Tu0 = 0.01
		Tu1 = 0.03
		Tu2 = 0.51
		Tu3 = 1.33
		Tu4 = 2.0
		Tu5 = 5.25
		Tu6 = 6.5

		m0 = (Re1-Re0)/(Tu1-Tu0)
		m1 = (Re2-Re1)/(Tu2-Tu1)
		m2 = (Re3-Re2)/(Tu3-Tu2)
		m3 = (Re4-Re3)/(Tu4-Tu3)
		m4 = (Re5-Re4)/(Tu5-Tu4)
		m5 = (Re6-Re5)/(Tu6-Tu5)

		IF(tuinf <= Tu1) THEN
			retinf = Re0 + m0*(tuinf-Tu0)
		ELSEIF(tuinf>Tu1 .AND. tuinf<=Tu2) THEN
			retinf = Re1 + m1 * (tuinf-Tu1)
		ELSEIF(tuinf>Tu2 .AND. tuinf<=Tu3) THEN
			retinf = Re2 + m2 * (tuinf-Tu2)
		ELSEIF(tuinf>Tu3 .AND. tuinf<=Tu4) THEN
			retinf = Re3 + m3 * (tuinf-Tu3)
		ELSEIF(tuinf>Tu4 .AND. tuinf<=Tu5) THEN
			retinf = Re4 + m4 * (tuinf-Tu4)
		ELSEIF(tuinf>Tu5 .AND. tuinf<=Tu6) THEN
			retinf = Re5 + m5 * (tuinf-Tu5)
		ELSE
			retinf = Re6
		ENDIF


        do k=1,kmax
        do j=1,jmax
          q(j,k,nv-1) = itmcinf!*rinf
          q(j,k,nv) = retinf!*rinf
        enddo
        enddo
      endif
c
      return
      end
c
c***********************************************************************
      subroutine restr2( q,jd,kd,logq )
c
c  read the restart file from unit 3 
c
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************
      integer jd,kd,logq
      real q(jd,kd,nq),reypr
      integer j,k
      ! local variables
      integer n
c***  first executable statement

      
      read(logq) fsmach,alf,reypr,totime
      read(logq) (((q(j,k,n),j=1,jmax),k=1,kmax),n=1,4)

      return
      end


c***********************************************************************
      subroutine restr_vnu( vnut,jd,kd,logq )

c     read in eddy viscosity values from restart file
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************
      integer jd,kd,logq,j,k
      real vnut(jd,kd)

      read(logq) ((vnut(j,k),j=1,jd),k=1,kd)

      return
      end

c***********************************************************************
      subroutine restr_komega( q,jd,kd,logq )

c     read in eddy viscosity values from restart file
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************
      integer jd,kd,logq,j,k,n
      real q(jd,kd,nq)

      read(logq) (((q(j,k,n),j=1,jd),k=1,kd),n=5,nv)

      return
      end

c***********************************************************************
      subroutine rhsup( q,s,xx,xy,yx,yy,x,y,xbig,ybig,xold,yold,
     &     xole,yole,iblank,ug,vg,jd,kd,nei,im,bt)
c
c  muscl approach:
c  qt = 0                       1st-order  ; |irhsy| = 1
c  qt = 0.25
c    th = 0     upwind-biased   2nd-order  ; |irhsy| = 2
c    th = 1/3   upwind-biased   3rd-order  ; |irhsy| = 3
c
c
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************
      integer jd,kd,nei,im
      real q(jd,kd,nq), s(jd,kd,nv)
      real xx(jd,kd), xy(jd,kd), yx(jd,kd), yy(jd,kd)
      real x(jd,kd), y(jd,kd)
      real ug(jd,kd), vg(jd,kd)
      real bt(jd,kd)
      real xbig(2*jd,2*kd),ybig(2*jd,2*kd)
      real xold(2*jd,2*kd),yold(2*jd,2*kd)
      real xole(2*jd,2*kd),yole(2*jd,2*kd)
      integer iblank(jd,kd)

      ! local variables
      real,allocatable :: ppp(:,:)
      real,allocatable :: tj(:), xa(:), ya(:)
      real,allocatable :: f(:,:),ql(:,:),qr(:,:)
      real,allocatable :: fmin(:),fmax(:)
      integer,allocatable :: iba(:)
      real,allocatable :: a(:),b(:),c(:),fin(:),fout(:)
      real,allocatable :: bbt(:)

      integer irhs,ilima,k,j,jj
      real th,qt,eps,epsj,epsk,rhoi,dx2,dy2,temp
      integer ibmin,ibmax

      allocate(ppp(jd,kd))
      allocate(tj(mdim), xa(mdim), ya(mdim))
      allocate(f(mdim,nmv),ql(mdim,nmv),qr(mdim,nmv))
      allocate(fmin(nmv),fmax(nmv),iba(mdim))
      allocate(a(mdim),b(mdim),c(mdim),fin(mdim),fout(mdim))
      allocate(bbt(mdim))

c**   first executable statement

      irhs   = iabs(irhsy)
      ilima = abs(ilim)
c     limter = 1
c     if( irhs .eq. 2 .and. limter.eq.1 ) limter = 2
      th     = real( irhs - 2 )/real( irhs )
      qt     = 0.25
      if( irhs .eq. 1 ) qt = 0.0
      eps    = 1.e-6
      epsj   =  (10./real(jm))**3
      epsk   =  (10./real(km))**3
c
      do 1 k = 1, kmax
      do 1 j = 1, jmax
        ppp(j,k) = gm1*( q(j,k,4)-0.5*(q(j,k,2)**2+q(j,k,3)**2)
     <                                                      /q(j,k,1) )
    1 continue
c
c..xi fluxes
c..nonconservative variables
c     
      do 13 k = 2,km
c
        do 15 j = 1,jmax
          rhoi   = 1.0/q(j,k,1)
          f(j,1) = q(j,k,1)*q(j,k,nq)
          f(j,2) = q(j,k,2)*rhoi
          f(j,3) = q(j,k,3)*rhoi
          f(j,4) = ppp(j,k)*q(j,k,nq)
          iba(j) = abs(iblank(j,k))
          bbt(j) = 1.!bt(j,k)
   15   continue
         
c..at boundaries
        ibmin = 1
        ibmax = 1
        if(half.eq.1) then
          ibmax = 2
          fmax(1) = f(jmax-3,1)
          fmax(2) = f(jmax-3,2)
          fmax(3) = -f(jmax-3,3)
          fmax(4) = f(jmax-3,4)
        endif
c..limit
        if(ilima.eq.0)then
          call iflux(f,ql,qr,1,jmax,jm,th,qt,eps,fmin,fmax,ibmin,ibmax,iba)
        endif
        if(ilima.eq.1)then
         call muscld_new(f,ql,qr,1,jmax,jm,th,qt,epsj,fmin,fmax,ibmin,ibmax,iba)
c         call muscld(f,ql,qr,1,jmax,jm,th,qt,eps,fmin,fmax,ibmin,ibmax)
        endif
        if(ilima.eq.2)then
          call quad(f,ql,qr,1,jmax,jm,th,qt,eps,fmin,fmax,ibmin,ibmax)
        endif
        if(ilim.eq.4)then
          call weno(f,ql,qr,1,jmax,jm,th,qt,eps,fmin,fmax,ibmin,ibmax,iba)
        endif
c
c..metric terms
c.do finite-volume like on finer mesh

        if(nei.eq.0) then
         do j = 1,jm
          xa(j) = (yold(2*j,2*k)-yold(2*j,2*k-2))
          ya(j) =-(xold(2*j,2*k)-xold(2*j,2*k-2))
         enddo
        else
         do j = 1,jm
          xa(j) = (ybig(2*j,2*k)-ybig(2*j,2*k-2))
          ya(j) =-(xbig(2*j,2*k)-xbig(2*j,2*k-2))
         enddo
        endif

	if(iunst.gt.0) then
        do j = 1,jm
          dx1 = xbig(j*2,k*2)-xold(j*2,k*2-1)
          dy1 = ybig(j*2,k*2)-yold(j*2,k*2-1)
          dx2 = xold(j*2,k*2)-xbig(j*2,k*2-1)
          dy2 = yold(j*2,k*2)-ybig(j*2,k*2-1)
          temp = -0.5*( dx1*dy2 - dx2*dy1 )/dt
          dx1 = xbig(j*2,k*2-1)-xold(j*2,k*2-2)
          dy1 = ybig(j*2,k*2-1)-yold(j*2,k*2-2)
          dx2 = xold(j*2,k*2-1)-xbig(j*2,k*2-2)
          dy2 = yold(j*2,k*2-1)-ybig(j*2,k*2-2)
          tj(j)= temp-0.5*( dx1*dy2 - dx2*dy1 )/dt
        enddo
        if(ntac.eq.2.and.istep.gt.1) then
         do j = 1,jm
          dx1 = xold(j*2,k*2)-xole(j*2,k*2-1)
          dy1 = yold(j*2,k*2)-yole(j*2,k*2-1)
          dx2 = xole(j*2,k*2)-xold(j*2,k*2-1)
          dy2 = yole(j*2,k*2)-yold(j*2,k*2-1)
          temp = -0.5*( dx1*dy2 - dx2*dy1 )/dt
          dx1 = xold(j*2,k*2-1)-xole(j*2,k*2-2)
          dy1 = yold(j*2,k*2-1)-yole(j*2,k*2-2)
          dx2 = xole(j*2,k*2-1)-xold(j*2,k*2-2)
          dy2 = yole(j*2,k*2-1)-yold(j*2,k*2-2)
          temp = temp-0.5*( dx1*dy2 - dx2*dy1 )/dt
          tj(j)= 1.5*tj(j)-0.5*temp
         enddo
        endif
 	else
 	do j = 1,jm
 	tj(j)=0.0
 	enddo
 	endif

c..compute the generalized numerical flux in roe!'s upwinding
c
        if (.not. iprecon) then
           call roeflx( f,ql,qr,xa,ya,tj,1,jm )
        else
           call roetrklflx( f,ql,qr,xa,ya,tj,1,jm,bbt)
        endif

c
        do 12 j = 2,jm
          s(j,k,1) = s(j,k,1) - ( f(j,1) - f(j-1,1) )
          s(j,k,2) = s(j,k,2) - ( f(j,2) - f(j-1,2) )
          s(j,k,3) = s(j,k,3) - ( f(j,3) - f(j-1,3) )
          s(j,k,4) = s(j,k,4) - ( f(j,4) - f(j-1,4) )
   12   continue

   13 continue
c
c..eta fluxes
c
      do 23 j = 2,jm
c
        do 25 k = 1,kmax
          rhoi   = 1.0/q(j,k,1)
          f(k,1) = q(j,k,1)*q(j,k,nq)
          f(k,2) = q(j,k,2)*rhoi
          f(k,3) = q(j,k,3)*rhoi
          f(k,4) = ppp(j,k)*q(j,k,nq)
          iba(k) = abs(iblank(j,k))
          bbt(k) = 1.!bt(j,k)
   25   continue

c..at boundaries
        ibmin = 1
        ibmax = 1
        if(half.eq.1 .and. j.le.jtail1) then
          ibmin = 2
          fmin(1) = f(2,1)
          fmin(2) = f(2,2)
          fmin(3) = -f(2,3)
          fmin(4) = f(2,4)
        endif
c        if(half.ne.1 .and. 
c     <          (j.le.jtail1 .or. j.ge.jtail2.and.bodyflag(im))) then
        if(half.ne.1 .and. (j.le.jtail1 .or. j.ge.jtail2)) then
			  if(flatplate .or. bodyflag(im)) then
          ibmin = 2
          jj = jmax-j+1
          rhoi   = 1.0/q(jj,2,1)
          fmin(1) = q(jj,2,1)*q(jj,2,nq)
          fmin(2) = q(jj,2,2)*rhoi
          fmin(3) = q(jj,2,3)*rhoi
          fmin(4) = ppp(jj,2)*q(jj,2,nq) 
        endif
		  endif
c..limit
        if(ilima.eq.0)then
          call iflux(f,ql,qr,1,jmax,jm,th,qt,eps,fmin,fmax,ibmin,ibmax,iba)
        endif
        if(ilima.eq.1)then
         call muscld_new(f,ql,qr,1,kmax,km,th,qt,epsk,fmin,fmax,ibmin,ibmax,iba)
c          call muscld(f,ql,qr,1,kmax,km,th,qt,eps,fmin,fmax,ibmin,ibmax)
        endif
        if(ilima.eq.2)then
          call quad(f,ql,qr,1,kmax,km,th,qt,eps,fmin,fmax,ibmin,ibmax)
        endif
	if(ilim.eq.4)then
          call weno(f,ql,qr,1,kmax,km,th,qt,eps,fmin,fmax,ibmin,ibmax,iba)
	endif

c..metric terms
c.do finite-volume like on finer mesh

        if(nei.eq.0) then
         do k = 1,km
          xa(k) =-(yold(2*j,2*k)-yold(2*j-2,2*k))
          ya(k) = (xold(2*j,2*k)-xold(2*j-2,2*k))
         enddo
        else
         do k = 1,km
          xa(k) =-(ybig(2*j,2*k)-ybig(2*j-2,2*k))
          ya(k) = (xbig(2*j,2*k)-xbig(2*j-2,2*k))
         enddo
        endif

 	if(iunst.gt.0) then
        do k = 1,km
          dx1 = xbig(j*2,k*2)-xold(j*2-1,k*2)
          dy1 = ybig(j*2,k*2)-yold(j*2-1,k*2)
          dx2 = xold(j*2,k*2)-xbig(j*2-1,k*2)
          dy2 = yold(j*2,k*2)-ybig(j*2-1,k*2)
          temp  = 0.5*( dx1*dy2 - dx2*dy1 )/dt
          dx1 = xbig(j*2-1,k*2)-xold(j*2-2,k*2)
          dy1 = ybig(j*2-1,k*2)-yold(j*2-2,k*2)
          dx2 = xold(j*2-1,k*2)-xbig(j*2-2,k*2)
          dy2 = yold(j*2-1,k*2)-ybig(j*2-2,k*2)
          tj(k) = temp+0.5*( dx1*dy2 - dx2*dy1 )/dt
        enddo
        if(ntac.eq.2.and.istep.gt.1) then
         do k = 1,km
          dx1 = xold(j*2,k*2)-xole(j*2-1,k*2)
          dy1 = yold(j*2,k*2)-yole(j*2-1,k*2)
          dx2 = xole(j*2,k*2)-xold(j*2-1,k*2)
          dy2 = yole(j*2,k*2)-yold(j*2-1,k*2)
          temp  = 0.5*( dx1*dy2 - dx2*dy1 )/dt
          dx1 = xold(j*2-1,k*2)-xole(j*2-2,k*2)
          dy1 = yold(j*2-1,k*2)-yole(j*2-2,k*2)
          dx2 = xole(j*2-1,k*2)-xold(j*2-2,k*2)
          dy2 = yole(j*2-1,k*2)-yold(j*2-2,k*2)
          temp = temp+0.5*( dx1*dy2 - dx2*dy1 )/dt
          tj(k)= 1.5*tj(k)-0.5*temp
         enddo
        endif
 	else
 	do k = 1,km
 	tj(k)=0.
 	enddo
 	endif
c
c..compute the generalized numerical flux in roe!'s upwinding
c
        if (.not. iprecon) then
           call roeflx( f,ql,qr,xa,ya,tj,1,km )
        else
           call roetrklflx( f,ql,qr,xa,ya,tj,1,km,bbt)
        endif
c
        do 22 k = 2,km
          s(j,k,1) = s(j,k,1) - ( f(k,1) - f(k-1,1) )
          s(j,k,2) = s(j,k,2) - ( f(k,2) - f(k-1,2) )
          s(j,k,3) = s(j,k,3) - ( f(k,3) - f(k-1,3) )
          s(j,k,4) = s(j,k,4) - ( f(k,4) - f(k-1,4) )
   22   continue

   23 continue
c
      return
      end

c***********************************************************************
      subroutine roeflx( f,ql,qr,xa,ya,tj,is,ie )
c
c  compute the generalized numerical flux in roe!'s upwinding
c  by s.o.                          
c  mod by jdb to incorporate smoother entropy check
c
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************

      integer is,ie
      real f(mdim,nmv)
      real tj(mdim), xa(mdim), ya(mdim)
      real ql(mdim,nmv),qr(mdim,nmv)

      ! local variables
      real eps,rlft,ulft,vlft,plft
      real rlfti,rulft,rvlft,uvl,elft,hlft,clft
      real rrht,urht,vrht,prht
      real rrhti,rurht,rvrht,uvr,erht,hrht,crht
      real tklft,tomegalft,tkrht,tomegarht
      real rat,rati,rav,uav,vav,hav,uv,cav,tkav,tomegaav
      real aq1,aq2,aq3,aq4,aq5,aq6,ri1,ri2,ri3,rr2,rr,r0,r2,r3
      real uu,c2,c2i,auu,aupc,aumc,uulft,uurht,upclft,upcrht
      real umclft,umcrht,dauu,dauus,daupc,daumc,daumcs,rcav,aquu
      real daupcs,c2ih,ruuav,b1,b2,b3,b4,b5,b6,b7,b8,b9,aj
      real plar,eplft,eprht,fssub

      integer i,i1
c..   nishan: variables for all mach correction
      integer write_allmach
      real fm

c***  first executable statement

      eps = 1.e-6
      do 11 i = is,ie
c
       i1    = i + 1
       rlft = ql(i,1)
       ulft = ql(i,2)
       vlft = ql(i,3)
       plft = ql(i,4)
       rlfti = 1.0/rlft
       rulft = rlft*ulft
       rvlft = rlft*vlft
       uvl = 0.5*( ulft*ulft + vlft*vlft )
       elft = plft/gm1 + rlft*uvl
       hlft = ( elft + plft )*rlfti
       clft = sqrt( gm1*( hlft - uvl ) )
c
       rrht = qr(i1,1)
       urht = qr(i1,2)
       vrht = qr(i1,3)
       prht = qr(i1,4)
       rrhti = 1.0/rrht
       rurht = rrht*urht
       rvrht = rrht*vrht
       uvr = 0.5*( urht*urht + vrht*vrht )
       erht = prht/gm1 + rrht*uvr
       hrht = ( erht + prht )*rrhti
       crht = sqrt( gm1*( hrht - uvr ) )
c
       rat  = sqrt( rrht*rlfti )
       rati = 1.0/( rat + 1. )
       rav  =   rat*rlft
       uav  = ( rat*urht + ulft )*rati
       vav  = ( rat*vrht + vlft )*rati
       hav  = ( rat*hrht + hlft )*rati
       uv   = 0.5*( uav*uav + vav*vav )
       cav  = sqrt( gm1*( hav - uv ) )
c
       aq1  = rrht - rlft
       aq2  = urht - ulft
       aq3  = vrht - vlft
       aq4  = prht - plft
c
       ri1 = xa(i)
       ri2 = ya(i)
       ri3 = tj(i)
       rr2 = ri1*ri1 + ri2*ri2
       rr  = sqrt( rr2 )
       r0  = 1.0 / rr
       r1  = ri1*r0
       r2  = ri2*r0
       r3  = ri3*r0
c
       uu  = r1*uav + r2*vav + r3
       c2  = cav*cav
       c2i = 1.0/c2
c
       auu   = abs( uu    )
       aupc  = abs( uu+cav )
       aumc  = abs( uu-cav )
c     
       uulft = r1*ulft + r2*vlft + r3
       uurht = r1*urht + r2*vrht + r3
       upclft= uulft + clft
       upcrht= uurht + crht
       umclft= uulft - clft
       umcrht= uurht - crht
c
       dauu = 4.*(uurht-uulft)+eps
       dauus = amax1(dauu,0.0)
ccray       auu = cvmgt(auu**2/dauu+0.25*dauu,auu,auu.le.0.5*dauus)
       if( auu.le.0.5*dauus ) then
         auu = auu**2/dauu+0.25*dauu
       end if
       daupc = 4.*(upcrht-upclft)+eps
       daupcs = amax1(daupc,0.0)
ccray       aupc = cvmgt(aupc**2/daupc+0.25*daupc,aupc,aupc.le.0.5*daupcs)
       if( aupc.le.0.5*daupcs ) then
         aupc = aupc**2/daupc+0.25*daupc
       end if
       daumc = 4.*(umcrht-umclft)+eps
       daumcs = amax1(daumc,0.0)
ccray       aumc = cvmgt(aumc**2/daumc+0.25*daumc,aumc,aumc.le.0.5*daumcs)
       if( aumc.le.0.5*daumcs ) then
         aumc = aumc**2/daumc+0.25*daumc
       end if
c     
       rcav = rav*cav
   
       aquu = uurht - uulft
        
	   if (iallmach) then
          fm=min(     max(  0.0001,sqrt(2.0*uv*c2i)   ),  1. )          
          aquu = fm*aquu
       end if

   
       c2ih = 0.5*c2i
       ruuav= auu*rav
       b1   = auu*( aq1 - c2i*aq4 )
       b2   = c2ih*aupc*( aq4 + rcav*aquu )
       b3   = c2ih*aumc*( aq4 - rcav*aquu )
       b4   = b1 + b2 + b3
       b5   = cav*( b2 - b3 )
       b6   = ruuav*( aq2 - r1*aquu )
       b7   = ruuav*( aq3 - r2*aquu )
c
       aq1 = b4
       aq2 = uav*b4 + r1*b5 + b6
       aq3 = vav*b4 + r2*b5 + b7
       aq4 = hav*b4 + ( uu-r3 )*b5 + uav*b6 + vav*b7 - c2*b1/gm1
c
       aj    = 0.5*rr
       plar  = plft + prht
       eplft = elft + plft
       eprht = erht + prht
       fssub = rr*r3
       fssub = 0.0
       f(i,1) = aj*( rlft*uulft+rrht*uurht-aq1 ) - fssub*rinf
       f(i,2) = aj*( rulft*uulft+rurht*uurht+r1*plar-aq2 )
     <                                          -fssub*rinf*uinf
       f(i,3) = aj*( rvlft*uulft+rvrht*uurht+r2*plar-aq3 )
     <                                          -fssub*rinf*vinf
       f(i,4) = aj*( eplft*uulft+eprht*uurht-r3*plar-aq4 )
     <                                          -fssub*einf

c      fssub = rr*(r1*uinf + r2*vinf + r3)
c      f(i,1) = aj*( rlft*uulft+rrht*uurht-aq1 ) - fssub*rinf
c      f(i,2) = aj*( rulft*uulft+rurht*uurht+r1*plar-aq2 ) 
c    <                            -fssub*rinf*uinf-rr*r1*pinf
c      f(i,3) = aj*( rvlft*uulft+rvrht*uurht+r2*plar-aq3 )
c    <                            -fssub*rinf*vinf-rr*r2*pinf
c      f(i,4) = aj*( eplft*uulft+eprht*uurht-r3*plar-aq4 )
c    <                            -fssub*(einf+pinf)+rr*r3*pinf
   11 continue
c
      return
      end


c***********************************************************************
      subroutine roetrklflx( f,ql,qr,xa,ya,tj,is,ie,b)
c
c  compute the generalized numerical flux in roe!'s upwinding
c  by s.o. and uses turkel preconditioning                         
c  mod by jdb to incorporate smoother entropy check
c
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************

      integer is,ie
      real f(mdim,nmv)
      real tj(mdim), xa(mdim), ya(mdim),b(mdim)
      real ql(mdim,nmv),qr(mdim,nmv)

      ! local variables
      real eps,rlft,ulft,vlft,plft
      real rlfti,rulft,rvlft,uvl,elft,hlft,clft
      real rrht,urht,vrht,prht
      real rrhti,rurht,rvrht,uvr,erht,hrht,crht
      real tklft,tomegalft,tkrht,tomegarht
      real rat,rati,rav,uav,vav,hav,uv,cav,tkav,tomegaav
      real aq1,aq2,aq3,aq4,aq5,aq6,ri1,ri2,ri3,rr2,rr,r0,r2,r3
      real uumxt,uu,c2,c2i,auu,aupc,aumc,uulft,uurht,upclft,upcrht
      real umclft,umcrht,dauu,dauus,daupc,daumc,daumcs,rcav,aquu
      real daupcs,c2ih,ruuav,b1,b2,b3,b4,b5,b6,b7,b8,b9,aj
      real plar,eplft,eprht,fssub
      real R,S,X,bSq

      integer i,i1

c***  first executable statement

      eps = 1.e-6
      do 11 i = is,ie
c
       i1    = i + 1
       rlft = ql(i,1)
       ulft = ql(i,2)
       vlft = ql(i,3)
       plft = ql(i,4)
       rlfti = 1.0/rlft
       rulft = rlft*ulft
       rvlft = rlft*vlft
       uvl = 0.5*( ulft*ulft + vlft*vlft )
       elft = plft/gm1 + rlft*uvl
       hlft = ( elft + plft )*rlfti
       clft = sqrt( gm1*( hlft - uvl ) )
c
       rrht = qr(i1,1)
       urht = qr(i1,2)
       vrht = qr(i1,3)
       prht = qr(i1,4)
       rrhti = 1.0/rrht
       rurht = rrht*urht
       rvrht = rrht*vrht
       uvr = 0.5*( urht*urht + vrht*vrht )
       erht = prht/gm1 + rrht*uvr
       hrht = ( erht + prht )*rrhti
       crht = sqrt( gm1*( hrht - uvr ) )
c
       rat  = sqrt( rrht*rlfti )
       rati = 1.0/( rat + 1. )
       rav  =   rat*rlft
       uav  = ( rat*urht + ulft )*rati
       vav  = ( rat*vrht + vlft )*rati
       hav  = ( rat*hrht + hlft )*rati
       uv   = 0.5*( uav*uav + vav*vav )
       cav  = sqrt( gm1*( hav - uv ) )
c
       aq1  = rrht - rlft
       aq2  = urht - ulft
       aq3  = vrht - vlft
       aq4  = prht - plft
c
       ri1 = xa(i)
       ri2 = ya(i)
       ri3 = tj(i)
       rr2 = ri1*ri1 + ri2*ri2
       rr  = sqrt( rr2 )
       r0  = 1.0 / rr
       r1  = ri1*r0
       r2  = ri2*r0
       r3  = ri3*r0
c
       uumxt  = r1*uav + r2*vav
       uu  = uumxt + r3
       c2  = cav*cav
       c2i = 1.0/c2
c
       bSq = Mp**2/(b(i)-Mp**2*(b(i)-1))

       X = sqrt( (1.-bSq)*uu*(1.-bSq)*uu+4.*bSq*c2 )
       auu   = abs( uu    )
       aupc  = 0.5*abs( (1.+bSq)*uu + X )
       aumc  = 0.5*abs( (1.+bSq)*uu - X )
c     
       uulft = r1*ulft + r2*vlft + r3
       uurht = r1*urht + r2*vrht + r3
       X = sqrt( (1.-bSq)*uulft*((1.-bSq)*uulft)+4.*bSq*clft*clft )
       upclft= 0.5*( (1.+bSq)*uulft + X )
       umclft= 0.5*( (1.+bSq)*uulft - X )
       X = sqrt( (1.-bSq)*uurht*((1.-bSq)*uurht)+4.*bSq*crht*crht )
       upcrht= 0.5*( (1.+bSq)*uurht + X )
       umcrht= 0.5*( (1.+bSq)*uurht - X )
c
       dauu = 4.*(uurht-uulft)+eps
       dauus = amax1(dauu,0.0)
ccray       auu = cvmgt(auu**2/dauu+0.25*dauu,auu,auu.le.0.5*dauus)
       if( auu.le.0.5*dauus ) then
         auu = auu**2/dauu+0.25*dauu
       end if
c
       daupc = 4.*(upcrht-upclft)+eps
       daupcs = amax1(daupc,0.0)
ccray       aupc = cvmgt(aupc**2/daupc+0.25*daupc,aupc,aupc.le.0.5*daupcs)
       if( aupc.le.0.5*daupcs ) then
         aupc = aupc**2/daupc+0.25*daupc
       end if
c
       daumc = 4.*(umcrht-umclft)+eps
       daumcs = amax1(daumc,0.0)
ccray       aumc = cvmgt(aumc**2/daumc+0.25*daumc,aumc,aumc.le.0.5*daumcs)
       if( aumc.le.0.5*daumcs ) then
         aumc = aumc**2/daumc+0.25*daumc
       end if
c     
       rcav = rav*cav
       aquu = uurht - uulft
       c2ih = 0.5*c2i
       ruuav= auu*rav
       X = sqrt( (1.-bSq)*uu*((1.-bSq)*uu)+4.*bSq*c2 )
       R     = 0.5*( (1.-bSq)*uu + X)
       S     = 0.5*( (1.-bSq)*uu - X )
       b1   = auu*( aq1 - c2i*aq4 )
       b2   = aupc*( aq4/R + rav*aquu )/X
       b3   = aumc*(-aq4/S - rav*aquu )/X
       b4   = b1 + b2 + b3
       b5   = ( R*b2 + S*b3 )
       b6   = ruuav*( aq2 - r1*aquu )
       b7   = ruuav*( aq3 - r2*aquu )
c
       aq1 = b4
       aq2 = uav*b4 + r1*b5 + b6
       aq3 = vav*b4 + r2*b5 + b7
       aq4 = hav*b4 + uumxt*b5 + uav*b6 + vav*b7 - c2*b1/gm1
c
       aj    = 0.5*rr
       plar  = plft + prht
       eplft = elft + plft
       eprht = erht + prht
       fssub = rr*r3
       fssub = 0.0
       f(i,1) = aj*( rlft*uulft+rrht*uurht-aq1 ) - fssub*rinf
       f(i,2) = aj*( rulft*uulft+rurht*uurht+r1*plar-aq2 )
     <                                          -fssub*rinf*uinf
       f(i,3) = aj*( rvlft*uulft+rvrht*uurht+r2*plar-aq3 )
     <                                          -fssub*rinf*vinf
       f(i,4) = aj*( eplft*uulft+eprht*uurht-r3*plar-aq4 )
     <                                          -fssub*einf

c      fssub = rr*(r1*uinf + r2*vinf + r3)
c      f(i,1) = aj*( rlft*uulft+rrht*uurht-aq1 ) - fssub*rinf
c      f(i,2) = aj*( rulft*uulft+rurht*uurht+r1*plar-aq2 ) 
c    <                            -fssub*rinf*uinf-rr*r1*pinf
c      f(i,3) = aj*( rvlft*uulft+rvrht*uurht+r2*plar-aq3 )
c    <                            -fssub*rinf*vinf-rr*r2*pinf
c      f(i,4) = aj*( eplft*uulft+eprht*uurht-r3*plar-aq4 )
c    <                            -fssub*(einf+pinf)+rr*r3*pinf
   11 continue
c
      return
      end

c***********************************************************************
      subroutine rhslom(q,s,jd,kd,js,je,ks,ke,bt)
c
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c*********************************************************************c

         real q(jd,kd,nq),s(jd,kd,nv),bt(jd,kd)
         integer jd,kd,js,je,ks,ke

      ! local variables
         real,allocatable :: tmp(:)

         integer j,k,n
         real u,v,ge,aSq,phiSq,bSq

         allocate(tmp(4))

c***  first executable statement

         do j = js,je
            do k = ks,ke
               u = q(j,k,2)/q(j,k,1)
               v = q(j,k,3)/q(j,k,1)
               phiSq = 0.5*( u*u + v*v )
               ge    = q(j,k,4)/q(j,k,1)*gamma-phiSq*gm1
               aSq   = (q(j,k,4)/q(j,k,1)-phiSq)*gm1*gamma
               do n = 1,4 
                 tmp(n)=s(j,k,n)
               enddo
               s(j,k,1) = phiSq*tmp(1)-u*tmp(2)-v*tmp(3)+tmp(4)
               s(j,k,2) = u*s(j,k,1)
               s(j,k,3) = v*s(j,k,1)
               s(j,k,4) = ge*s(j,k,1)
      
               bSq = Mp**2/(bt(j,k)-Mp**2*(bt(j,k)-1))

               do n = 1,4
                  s(j,k,n) = s(j,k,n)*(gm1*(bSq-1.)/aSq)
                  s(j,k,n) = s(j,k,n) + tmp(n)
               enddo

           enddo
        enddo

      return
      end

c***********************************************************************
      subroutine step(q,qtn,qtnm1,qnewt,s,
     &     x,y,iblank,
     &     xx,xy,yx,yy,ug,vg,
     &     yx0,yy0,yt0,
     &     xbig,ybig,xold,yold,xole,yole,
     &     vnut,vnut0,turmu,tscale,bt,
     &     im,jd,kd,resmax,resrho,rsum,cl,
     &     tau_dim,tau_global,utau_global)

c  note that the boundaries are updated explicitly
c
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************
      integer im,jd,kd
      real resmax,resrho,rsum,cl
      real q(jd,kd,nq), s(jd,kd,nv),tscale(jd,kd),bt(jd,kd)
      real vnut(jd,kd), vnut0(jd,kd),turmu(jd,kd)
      real x(jd,kd),y(jd,kd),xx(jd,kd),xy(jd,kd),yx(jd,kd)
      real yy(jd,kd),ug(jd,kd),vg(jd,kd)
      real yx0(jd), yy0(jd), yt0(jd)
      real xbig(2*jd,2*kd),ybig(2*jd,2*kd)
      real xold(2*jd,2*kd),yold(2*jd,2*kd)
      real xole(2*jd,2*kd),yole(2*jd,2*kd)
      real qtn(jd,kd,nv),qtnm1(jd,kd,nv)
      real qnewt(jd,kd,nv)
      integer iblank(jd,kd)

		! friction variables
		integer tau_dim
		real tau_global(tau_dim),utau_global(tau_dim)

      ! local variables
      real,allocatable :: stn(:,:,:)

      integer k,j,n,jmax_j,jmax_k,kmax_j,jd1,kd1,nd1,kmax_k
      real smrs,smrs1,smrs2,smrs3,smrs4,volum
      integer write_source_bodyforce

      allocate(stn(jd,kd,4))

c***  first executable statement

!		print*,'from step = ',size(tau_global,1),tau_dim

      if(ntac.eq.-2.and.itn.eq.1) then
        do k=1,kmax
        do j=1,jmax
        do n=1,nv
          stn(j,k,n)=0.
        enddo
        enddo
        enddo
c
        call rhsup(q,stn,xx,xy,yx,yy,x,y,xbig,ybig,xold,yold,
     &       xole,yole,iblank,ug,vg,jd,kd,0,im,bt)


        if( .not. invisc ) then
          call visrhs(turmu,q,x,y,stn,xx,xy,yx,yy,ug,vg,vnut,
     &          vnut0,tscale,iblank,jd,kd,im,tau_dim,tau_global,utau_global)
        endif
      endif

c..zero s array 
c
      do k = 1,kmax
      do j = 1,jmax
      do n = 1,nv
        s(j,k,n) = 0.
      enddo
      enddo
      enddo

      call bc(q,x,y,xx,xy,yx,yy,ug,vg,yx0,yy0,yt0,im,jd,kd,bt)
c     
c..compute the right hand side and store in s array
c
      call rhsup(q,s,xx,xy,yx,yy,x,y,xbig,ybig,xold,yold,
     &       xole,yole,iblank,ug,vg,jd,kd,1,im,bt)

c.. nishan.. add additional source term in the form of body force
       write_source_bodyforce = 0
       if (write_source_bodyforce == 1) then
           call source_bodyforce(q,s)
       endif
c
c..start newton iteration here at each time step for convergence
c
      if (itnmax.gt.1) call newton(q,qtn,qtnm1,qnewt,s,2,km,jd,kd)
c
c..compute viscous fluxes
c
      if( .not. invisc ) then
        call visrhs( turmu,q,x,y,s,xx,xy,yx,yy,ug,vg,vnut,vnut0,
     >		tscale,iblank,jd,kd,im,tau_dim,tau_global,utau_global)
      endif
c     
c..make sure to account for time accuracy
c     
      if (ntac.eq.-2) then
        do k=1,kmax
        do j=1,jmax
        do n=1,4
          s(j,k,n)=0.5*(s(j,k,n)+stn(j,k,n))
        enddo
        enddo
        enddo
      endif
c
      rsum = 0.
      resrho  = 0.0
      resmax  = 0.0
      jmax_j=2
      jmax_k=2
      kmax_j=2
      jd1=2
      kd1=2
      nd1=1

      kmax_k=2
      do 23 k = 2,km
c
        do 22 n = 1,4
        do 22 j = 2,jm
          smrs = abs(s(j,k,n))*max(iblank(j,k),0)
          if (smrs .gt. resmax) then
            jd1 = j
            kd1 = k
            nd1 = n
            resmax = smrs
          endif
  22    continue
c
        do 24 j = 2,jm
          smrs1 = s(j,k,1)*max(iblank(j,k),0)
          smrs2 = s(j,k,2)*max(iblank(j,k),0)
          smrs3 = s(j,k,3)*max(iblank(j,k),0)
          smrs4 = s(j,k,4)*max(iblank(j,k),0)
          resrho  = resrho + s(j,k,1)**2
          rsum = rsum + smrs1*smrs1 + smrs2*smrs2 + smrs3*smrs3
     <                + smrs4*smrs4
          s(j,k,1) = smrs1*tscale(j,k)
          s(j,k,2) = smrs2*tscale(j,k)
          s(j,k,3) = smrs3*tscale(j,k)
          s(j,k,4) = smrs4*tscale(j,k)
  24    continue
c
  23  continue
c
      volum  = float(jmax*kmax)
      resrho = sqrt(resrho/volum)
      rsum   = sqrt(rsum/volum)
c
      if( mod(istep0,npnorm).eq.0 ) then
         write(70+im,101) istep0,rsum,resrho,resmax,totime,theta_col
  101    format(i7,5(e18.10))
         write(6,102) jd1,kd1,nd1,resmax,resrho,rsum
  102    format('  j,k,n,rmax,l2rho,l2 =',3i4,3(x,e10.4))
      endif
c
c..check on convergence
c
      if( rsum .gt.1000.0 ) then
        write(6,602)  rsum
  602   format(' ',10x,'norm is out of bounds,'
     $         ,1x,'l2norm = ',f16.12,1x,'solution suspended' )
        stop 'norm'
      end if
      if (iprecon) call rhslom(q,s,jd,kd,2,jm,2,km,bt)
c
c..now do the implicit part
c  
      if(ilhs.eq.1) then

         if(.not.iprecon) then
           call ilu2d(q,s,jd,kd,2,jm,2,km,xx,xy,yx,yy,ug,vg,turmu,tscale)
         else
           call preilu2d(q,s,jd,kd,2,jm,2,km,xx,xy,yx,yy,ug,vg,turmu,tscale,bt)
         endif  

      elseif(ilhs.eq.2.or.ilhs.eq.3) then
         if(.not.iprecon) then 
            call arc2d(q,s,jd,kd,2,jm,2,km,xx,xy,yx,yy,ug,vg,turmu,
     c                 iblank,tscale,bt)
         else
            call arc2d_precon(q,s,jd,kd,2,jm,2,km,xx,xy,yx,yy,ug,vg,turmu,
     c                 iblank,tscale,bt)
         endif
      endif

c..update q with corrections
c
cjdb      do 31 k = 2,km
      do 31 k = 2,km
      do 31 j = 2,jm
        q(j,k,1) = q(j,k,1) + s(j,k,1)*max(iblank(j,k),0)
        q(j,k,2) = q(j,k,2) + s(j,k,2)*max(iblank(j,k),0)
        q(j,k,3) = q(j,k,3) + s(j,k,3)*max(iblank(j,k),0)
        q(j,k,4) = q(j,k,4) + s(j,k,4)*max(iblank(j,k),0)
   31 continue

c..update b c

      call bc(q,x,y,xx,xy,yx,yy,ug,vg,yx0,yy0,yt0,im,jd,kd,bt)

      return
      end

c***********************************************************************
      subroutine storevnu(x,y,iblank,q,vnut,jd,kd,logq,logtur)

c      write the turbulent viscosity out
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************
      integer jd,kd,logq,logtur
      real vnut(jd,kd), q(jd,kd,nq), x(jd,kd), y(jd,kd)
      integer iblank(jd,kd)
      ! local variables
      integer j,k

!      real,allocatable :: rho(:,:)
!      allocate(rho(jd,kd))

!      do j = 1,jd
!      do k = 1,kd
!        rho(j,k) = q(j,k,1)*q(j,k,nq)
!      enddo
!      enddo

      write (logq) ((vnut(j,k),j=1,jd),k=1,kd)

      write(logtur,*) "ZONE"
      write(logtur,*) "I = ",jmax,", J = ",kmax
      write(logtur,*) "ZONETYPE = Ordered, DATAPACKING = BLOCK" 
      write(logtur,*) ((x(j,k),j=1,jmax),k=1,kmax)
      write(logtur,*) ((y(j,k),j=1,jmax),k=1,kmax)
      write(logtur,*) ((iblank(j,k),j=1,jmax),k=1,kmax)
      write(logtur,*) ((vnut(j,k),j=1,jmax),k=1,kmax)
      if (itrans.eq.1) then
        write(logtur,*) ((q(j,k,nv-1),j=1,jmax),k=1,kmax)
        write(logtur,*) ((q(j,k,nv),j=1,jmax),k=1,kmax)
      endif

      return
      end

c***********************************************************************
      subroutine storekomega(x,y,iblank,q,vnut,jd,kd,logq,logtur)

c      write the turbulent viscosity out
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************
      integer jd,kd,logq,logtur
      real q(jd,kd,nq), vnut(jd,kd), x(jd,kd), y(jd,kd)
      integer iblank(jd,kd)
      ! local variables
      integer j,k,n

      real,allocatable :: rho(:,:)
      allocate(rho(jd,kd))

      do j = 1,jd
      do k = 1,kd
        rho(j,k) = q(j,k,1)*q(j,k,nq)
      enddo
      enddo

      write(logq) (((q(j,k,n),j=1,jd),k=1,kd),n=5,nv)

      write(logtur,*) "ZONE"
      write(logtur,*) "I = ",jmax,", J = ",kmax
      write(logtur,*) "ZONETYPE = Ordered, DATAPACKING = BLOCK" 
      write(logtur,*) ((x(j,k),j=1,jmax),k=1,kmax)
      write(logtur,*) ((y(j,k),j=1,jmax),k=1,kmax)
      write(logtur,*) ((iblank(j,k),j=1,jmax),k=1,kmax)
      write(logtur,*) ((vnut(j,k),j=1,jmax),k=1,kmax)
      write(logtur,*) ((q(j,k,5)/rho(j,k),j=1,jmax),k=1,kmax)
      write(logtur,*) ((q(j,k,6)/rho(j,k),j=1,jmax),k=1,kmax)
      if (itrans.eq.1) then
        write(logtur,*) ((q(j,k,nv-1),j=1,jmax),k=1,kmax)
        write(logtur,*) ((q(j,k,nv),j=1,jmax),k=1,kmax)
      endif

      return
      end

c***********************************************************************
      subroutine store( q,x,y,iblank,jd,kd,logq,logg)
c
c  write the solution out
c
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************

      integer jd,kd,logq,logg
      real q(jd,kd,nq), x(jd,kd), y(jd,kd)
      integer iblank(jd,kd)

      ! local variables
      integer j,k,n
      real fstip,reypr

      ! local variables.. nishan for writing time wise solution
      integer write_time_wise
      character(len=200) :: wrfile,wrnum,wrfull


c***  first executable statement


      do 11 k = 1,kmax
      do 11 j = 1,jmax
        q(j,k,1) = q(j,k,1)*q(j,k,nq)
        q(j,k,2) = q(j,k,2)*q(j,k,nq)
        q(j,k,3) = q(j,k,3)*q(j,k,nq)
        q(j,k,4) = q(j,k,4)*q(j,k,nq)
   11 continue

      fstip = fsmach
      reypr = rey * fsmach
      write(logq) fstip,alf,reypr,totime
      write(logq) (((q(j,k,n),j=1,jmax),k=1,kmax),n=1,4)

	if(num_grids.gt.1) then
      write(logg)((x(j,k),j=1,jmax),k=1,kmax),
     <           ((y(j,k),j=1,jmax),k=1,kmax),
     <           ((iblank(j,k),j=1,jmax),k=1,kmax)
	else
      write(logg)((x(j,k),j=1,jmax),k=1,kmax),
     <           ((y(j,k),j=1,jmax),k=1,kmax)
	endif

    
      write_time_wise = 0
      if (write_time_wise == 1) then
         write(wrnum,*)istep0
         wrfile = 'fort.8_'
         wrfull = trim(adjustl(wrfile))//trim(adjustl(wrnum))
         open (unit=23, file = trim(adjustl(wrfull)), form = 'unformatted',status='replace')
         write(23) jmax,kmax
         write(23) fstip,alf,reypr,totime
         write(23) (((q(j,k,n),j=1,jmax),k=1,kmax),n=1,4)
         write(23) istep0
         close(23)
       end if

    

c..scale q back with jacobian

      call qdivj( q,jd,kd)

      return
      end

c***********************************************************************
      subroutine stqol(q,qtn,qtnm1,qnewt,vnut,vnut0,jd,kd)
c     
c  store data at previous time steps
c
c*********************************************************************** 
      use params_global
c*********************************************************************** 
      implicit none
c*********************************************************************** 
      integer jd,kd
      real q(jd,kd,nq),qtn(jd,kd,nv),qtnm1(jd,kd,nv)
      real qnewt(jd,kd,nv),vnut(jd,kd),vnut0(jd,kd)

      ! local variables
      integer j,k,n
c**   first executable statement

      do k = 1, kmax

        do j = 1, jmax
          do n = 1,nv
            qnewt(j,k,n) = q(j,k,n)
          enddo
        enddo

        if ( (ntac.eq.2 .and. istep.gt.1) .or. 
     <       (ntac.eq.3 .and. istep.eq.2) ) then
          do j = 1, jmax
            do n = 1,nv
              qnewt(j,k,n) = qnewt(j,k,n)+(q(j,k,n)-qtn(j,k,n))/3.
            enddo
          enddo
        endif

        if (ntac.eq.3 .and. istep.gt.2) then
          do j = 1, jmax 
            do n = 1,nv
              qnewt(j,k,n) = qnewt(j,k,n)+7.*(q(j,k,n)-qtn(j,k,n))/11.-
     <                       2.*(qtn(j,k,n)-qtnm1(j,k,n))/11.
            enddo
          enddo
        endif

        do j = 1, jmax
          do n = 1,nv
            qtnm1(j,k,n) = qtn(j,k,n)
          enddo
        enddo

        do j = 1, jmax
          do n = 1,nv
            qtn(j,k,n) = q(j,k,n)
          enddo
          vnut0(j,k)  = vnut(j,k)
        enddo

      enddo
c
      return
      end

c***********************************************************************
      subroutine tridag(a,b,c,f,z,ni,nl)

c     tri diagonal matrix inversion
c***********************************************************************
      implicit none
c***********************************************************************
      integer jd,ni,nl
      parameter (jd=1001)
      real a(jd),b(jd),c(jd),f(jd),z(jd)

      ! local variables
      real,allocatable :: w(:),g(:)

      integer nipl,j,nd,j1
      real d,rd

      allocate(w(jd),g(jd))
c
      w(ni)=c(ni)/b(ni)
      g(ni)=f(ni)/b(ni)
      nipl=ni+1
      do 10 j=nipl,nl
      d=b(j)-a(j)*w(j-1)
      rd=1.0/d
      w(j)=c(j)*rd
      g(j)=(f(j)-a(j)*g(j-1))*rd
 10   continue
      z(nl)=g(nl)
      nd=nl-ni
      do 20 j1=1,nd
      j=nl-j1
      z(j)=g(j)-w(j)*z(j+1)
 20   continue
      return
      end

c***********************************************************************
      subroutine uv(js,je,ks,ke,q,xx,xy,yx,yy,u,v,jd,kd,scal,iadir)
c
c  solve for the cartesian momemtum components (ru,rv) from
c  the contravariant velocity components (u,v) along lines of j.
c     
c  note: consistent with metfv     
c
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************
      integer jd,kd
      real xx(jd,kd),xy(jd,kd),yx(jd,kd),yy(jd,kd)
      real q(jd,kd,nq)
      real u(mdim), v(mdim),scal
      integer js,je,ks,ke,iadir
      ! local variables
      integer k,j
      real degtorad,radtodeg,angl,thet
      real t11,t12,t21,t22,vol,uwall,vwall,thet1
      real exb,eyb,exa,eya
	  real uw,vw

c***  first executable statement

	do j = js,je
		do k = ks,ke

        	t11   = (yy(j,k))
        	t12   =-(xy(j,k))
c
        	t21   =-(yx(j,k))
        	t22   = (xx(j,k))
c
        	vol=1./(t22*t11-t12*t21)
c
c..the following are the physical plane velocities u & v
			if (iadir.eq.1) then
				uw = u(k)
				vw = v(k)
			elseif (iadir.eq.2) then
				uw = u(j)
				vw = v(j)
			endif
c
        	uwall = (t11*uw +t12*vw)*vol
        	vwall = (t21*uw +t22*vw)*vol
c
c..physical plane velocities multiplied by the density to form q
c
        	q(j,k,2) = scal*uwall*q(j,k,1)+(1.-scal)*q(j,k,2)
        	q(j,k,3) = scal*vwall*q(j,k,1)+(1.-scal)*q(j,k,3)
		
		enddo
	enddo
	        

      return
      end

c***********************************************************************
      subroutine visrhs(turmu,q,x,y,s,xx,xy,yx,yy,ug,vg,vnut,vnut0,
     >                  tscale,iblank,jd,kd,im,tau_dim,tau_global,utau_global)
c
c  compute the viscous rhs 
c
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************
      integer jd,kd,im
      real q(jd,kd,nq), s(jd,kd,nv), turmu(jd,kd), tscale(jd,kd)
      real x(jd,kd),y(jd,kd),ug(jd,kd),vg(jd,kd)
      real xx(jd,kd), xy(jd,kd), yx(jd,kd), yy(jd,kd)
      real vnut(jd,kd),vnut0(jd,kd)
      integer iblank(jd,kd)
c
      ! local variables
      real,allocatable :: f2(:),f3(:),f4(:)
      real,allocatable :: f2v(:),f3v(:),f4v(:),u(:),v(:),e(:)
      real,allocatable :: a1v(:),a2v(:),a4v(:),a5v(:)
      real,allocatable :: b1v(:),b2v(:),b4v(:),b5v(:),d1v(:)
      real,allocatable :: d2v(:),d4v(:),d5v(:),d7v(:)

      real gkpr,prtr,dr,c2b,c2bp,t4,t5,t6,t7,t8
      integer itlns,j,k,km2,km1,k1,k2,jm2,jm1,j1,j2
      real ra,tt,vmue,turm,vnu,gkap,rj,b1,b2,b5,b4
      real du,dv,dei,rjm2,rjm,rjp,rjp2,d1,d2,d4,d5,d7
      real ujp2,ujp1,ujm1,ujm2,dvj
      real ejm2,ejm1,dej,rk,rkm2,rkm,rkp,rkp2,a1,a2,a4,a5
      real ukp2,ukp1,ukm1,ukm2,duk
      real vkp2,vkp1,vkm1,vkm2,dvk,ekp2,ekp1,ekm1,ekm2,dek
      real dre,duj,vjp2,vjp1,vjm1,vjm2,ejp1,ejp2,dv2,fkleb
      real oneby3,oneby6,oneby8,oneby12,oneby24
      integer ifourth

		! friction variables
		integer tau_dim
		real tau_global(tau_dim),utau_global(tau_dim)

      allocate(f2(mdim),f3(mdim),f4(mdim))
      allocate(f2v(mdim),f3v(mdim),f4v(mdim),u(mdim),v(mdim),
     <         e(mdim))
      allocate(a1v(mdim),a2v(mdim),a4v(mdim),a5v(mdim))
      allocate(b1v(mdim),b2v(mdim),b4v(mdim),b5v(mdim),d1v(mdim))
      allocate(d2v(mdim),d4v(mdim),d5v(mdim),d7v(mdim))

c*** first executable statement

!		print*,'from visrhs = ',size(tau_global,1),tau_dim

      if (iturb.eq.0) then
        if (bodyflag(im) .or. flatplate) call vmutur( x,y,q,s,turmu,xx,xy,yx,yy,ug,vg,jd,kd)
      elseif (iturb.eq.1) then
        call vmu_sa( x,y,q,turmu,xx,xy,yx,yy,
     <               ug,vg,jd,kd,vnut,vnut0,tscale,iblank,im,tau_dim,tau_global,utau_global)
      elseif (iturb.eq.2) then
        call komegasst( q,s,turmu,x,y,xx,xy,yx,yy,
     <                  ug,vg,jd,kd,tscale,iblank,im,tau_dim,utau_global)
      endif

      gkpr = gamma/pr
      prtr = pr/0.9
      dre  = .5/rey
      c2b  =198.6/tinf
      c2bp = c2b +1.

	ifourth = 0.
	itlns   = 0.

	if(ifourth.eq.0) then

	if(itlns.eq.1) then



      do 10 j = 2,jm
        do 20 k = 1,kmax
          ra    = 1./q(j,k,1)
          u(k)  = q(j,k,2)*ra
          v(k)  = q(j,k,3)*ra
          e(k)  = q(j,k,4)*ra-.5*(u(k)**2+v(k)**2)
          tt    = ggm1*e(k)
          vmue  = c2bp*tt*sqrt(tt)/(c2b + tt)
          turm  = turmu(j,k)
          vnu   = vmue+turm
          gkap  = vmue+prtr*turm
          rj    = 1./q(j,k,nq)
          b4v(k) = (yx(j,k)**2+yy(j,k)**2)*rj
          b1v(k) = (b4v(k)+yx(j,k)**2/3.*rj)*vnu*dre
          b2v(k) = (b4v(k)+yy(j,k)**2/3.*rj)*vnu*dre
          b5v(k) = (yx(j,k)*yy(j,k)/3.*rj)*vnu*dre
          b4v(k) = b4v(k)*gkpr*gkap*dre
   20   continue
        do 30 k = 1,km
          k1    = k+1

          b1    = b1v(k1)+b1v(k)
          b2    = b2v(k1)+b2v(k)
          b5    = b5v(k1)+b5v(k)
          b4    = b4v(k1)+b4v(k)

          t4    = u(k1)*b1v(k1)+u(k)*b1v(k)
          t5    = v(k1)*b2v(k1)+v(k)*b2v(k)
          t6    = u(k1)*b5v(k1)+u(k)*b5v(k)
          t7    = v(k1)*b5v(k1)+v(k)*b5v(k)
          t8    = b4v(k1)+b4v(k)

          du    = u(k1)-u(k)
          dv    = v(k1)-v(k)
          dei   = e(k1)-e(k)
          f2(k) = b1*du+b5*dv
          f3(k) = b5*du+b2*dv
          f4(k) = (t4+t7)*du+(t5+t6)*dv+t8*dei
   30   continue
c
        do 40 k = 2,km
          s(j,k,2) = s(j,k,2) + (f2(k)-f2(k-1))
          s(j,k,3) = s(j,k,3) + (f3(k)-f3(k-1))
          s(j,k,4) = s(j,k,4) + (f4(k)-f4(k-1))
   40   continue
   10 continue



	else


  	do 110 j = 2,jm
        do 120 k = 1,kmax
          ra    = 1./q(j,k,1)
          u(k)  = q(j,k,2)*ra
          v(k)  = q(j,k,3)*ra
          e(k)  = q(j,k,4)*ra-.5*(u(k)**2+v(k)**2)
          tt    = ggm1*e(k)
          vmue  = c2bp*tt*sqrt(tt)/(c2b + tt)
          turm  = turmu(j,k)
          vnu   = vmue+turm
          gkap  = vmue+prtr*turm
          rj    = 1./q(j,k,nq)
          rjm    = 1./q(j-1,k,nq)
          rjp    = 1./q(j+1,k,nq)

          b4v(k) = (yx(j,k)**2+yy(j,k)**2)*rj
          b1v(k) = (b4v(k)+yx(j,k)**2/3.*rj)*vnu*dre
          b2v(k) = (b4v(k)+yy(j,k)**2/3.*rj)*vnu*dre
          b5v(k) = (yx(j,k)*yy(j,k)/3.*rj)*vnu*dre
          b4v(k) = b4v(k)*gkpr*gkap*dre

        d1v(k) = (4./3.*xx(j-1,k)*yx(j-1,k)+xy(j-1,k)*yy(j-1,k))*rjm
        d1v(k) = d1v(k)+(4./3.*xx(j+1,k)*yx(j+1,k)+xy(j+1,k)*
     1  yy(j+1,k))*rjp
        d1v(k) = vnu*dre*d1v(k)

        d2v(k) = (4./3.*xy(j-1,k)*yy(j-1,k)+xx(j-1,k)*yx(j-1,k))*rjm
        d2v(k) = d2v(k)+(4./3.*xy(j+1,k)*yy(j+1,k)+xx(j+1,k)*
     1  yx(j+1,k))*rjp
        d2v(k) = vnu*dre*d2v(k)

        d4v(k) = (xy(j-1,k)*yy(j-1,k)+xx(j-1,k)*yx(j-1,k))*rjm
        d4v(k) = d4v(k)+(xy(j+1,k)*yy(j+1,k)+xx(j+1,k)*
     1  yx(j+1,k))*rjp
        d4v(k) = d4v(k)*gkpr*gkap*dre

        d5v(k) = (-2./3.*xy(j-1,k)*yx(j-1,k)+xx(j-1,k)*yy(j-1,k))*rjm
        d5v(k) = d5v(k)+(-2./3.*xy(j+1,k)*yx(j+1,k)+xx(j+1,k)*
     1  yy(j+1,k))*rjp
        d5v(k) = vnu*dre*d5v(k)

        d7v(k) = (-2./3.*xx(j-1,k)*yy(j-1,k)+xy(j-1,k)*yx(j-1,k))*rjm
        d7v(k) = d7v(k)+(-2./3.*xx(j+1,k)*yy(j+1,k)+xy(j+1,k)*
     1  yx(j+1,k))*rjp
        d7v(k) = vnu*dre*d7v(k)

120     continue
        do 130 k = 1,km
          k1    = k+1

          b1    = b1v(k1)+b1v(k)
          b2    = b2v(k1)+b2v(k)
          b5    = b5v(k1)+b5v(k)
          b4    = b4v(k1)+b4v(k)
          d1    = d1v(k)
          d2    = d2v(k)
          d4    = d4v(k)
          d5    = d5v(k)
          d7    = d7v(k)
          ujp1  = q(j+1,k,2)/q(j+1,k,1)
          ujm1  = q(j-1,k,2)/q(j-1,k,1)
          duj   = 0.5*(ujp1-ujm1)

          vjp1  = q(j+1,k,3)/q(j+1,k,1)
          vjm1  = q(j-1,k,3)/q(j-1,k,1)
          dvj   = 0.5*(vjp1-vjm1)

          ejp1  = q(j+1,k,4)/q(j+1,k,1)-.5*(ujp1**2+vjp1**2)
          ejm1  = q(j-1,k,4)/q(j-1,k,1)-.5*(ujm1**2+vjm1**2)
          dej   = 0.5*(ejp1-ejm1)

          t4    = u(k1)*b1v(k1)+u(k)*b1v(k)
          t5    = v(k1)*b2v(k1)+v(k)*b2v(k)
          t6    = u(k1)*b5v(k1)+u(k)*b5v(k)
          t7    = v(k1)*b5v(k1)+v(k)*b5v(k)
          t8    = b4v(k1)+b4v(k)


          du    = u(k1)-u(k)
          dv    = v(k1)-v(k)
          dei   = e(k1)-e(k)
          f2(k) = b1*du+b5*dv
          f2v(k)= d1*duj+d5*dvj
          f3(k) = b5*du+b2*dv
          f3v(k)= d7*duj+d2*dvj
          f4(k) = (t4+t7)*du+(t5+t6)*dv+t8*dei
          f4v(k) = 0.5*d1*0.5*(ujp1**2-ujm1**2)
          f4v(k) = f4v(k)+0.5*d2*0.5*(vjp1**2-vjm1**2)
          f4v(k) = f4v(k)+d5*0.5*(ujp1+ujm1)*dvj
          f4v(k) = f4v(k)+d7*0.5*(vjp1+vjm1)*duj
          f4v(k) = f4v(k)+d4*dej
130     continue
c
        do 140 k = 2,km-1
          s(j,k,2) = s(j,k,2)+(f2(k)-f2(k-1))+0.5*(f2v(k+1)-f2v(k-1))
          s(j,k,3) = s(j,k,3)+(f3(k)-f3(k-1))+0.5*(f3v(k+1)-f3v(k-1))
          s(j,k,4) = s(j,k,4)+(f4(k)-f4(k-1))+0.5*(f4v(k+1)-f4v(k-1))
140     continue

          k=km
          s(j,k,2) = s(j,k,2)+(f2(k)-f2(k-1))+0.5*(f2v(k-2)
     1          -4.*f2v(k-1)+3*f2v(k))
          s(j,k,3) = s(j,k,3)+(f3(k)-f3(k-1))+0.5*(f3v(k-2)
     1          -4.*f3v(k-1)+3*f3v(k))
          s(j,k,4) = s(j,k,4)+(f4(k)-f4(k-1))+0.5*(f4v(k-2)
     1          -4.*f4v(k-1)+3.*f4v(k))




110 	continue
c
	   do 210 k = 2,km
        do 220 j = 1,jmax
          ra    = 1./q(j,k,1)
          u(j)  = q(j,k,2)*ra
          v(j)  = q(j,k,3)*ra
          e(j)  = q(j,k,4)*ra-.5*(u(j)**2+v(j)**2)
          tt    = ggm1*e(j)
          vmue  = c2bp*tt*sqrt(tt)/(c2b + tt)
          turm  = turmu(j,k)
          vnu   = vmue+turm
          gkap  = vmue+prtr*turm
          rj    = 1./q(j,k,nq)
          rkm    = 1./q(j,k-1,nq)
          rkp    = 1./q(j,k+1,nq)

          a4v(j) = (xx(j,k)**2+xy(j,k)**2)*rj
          a1v(j) = (a4v(j)+xx(j,k)**2/3.*rj)*vnu*dre
          a2v(j) = (a4v(j)+xy(j,k)**2/3.*rj)*vnu*dre
          a5v(j) = (xx(j,k)*xy(j,k)/3.*rj)*vnu*dre
          a4v(j) = a4v(j)*gkpr*gkap*dre

        d1v(j) = (4./3.*xx(j,k-1)*yx(j,k-1)+xy(j,k-1)*yy(j,k-1))*rkm
        d1v(j) = d1v(j)+(4./3.*xx(j,k+1)*yx(j,k+1)+xy(j,k+1)*
     1  yy(j,k+1))*rkp
        d1v(j) = vnu*dre*d1v(j)

        d2v(j) = (4./3.*xy(j,k-1)*yy(j,k-1)+xx(j,k-1)*yx(j,k-1))*rkm
        d2v(j) = d2v(j)+(4./3.*xy(j,k+1)*yy(j,k+1)+xx(j,k+1)*
     1  yx(j,k+1))*rkp
        d2v(j) = vnu*dre*d2v(j)

        d4v(j) = (xy(j,k-1)*yy(j,k-1)+xx(j,k-1)*yx(j,k-1))*rkm
        d4v(j) = d4v(j)+(xy(j,k+1)*yy(j,k+1)+xx(j,k+1)*
     1  yx(j,k+1))*rkp
        d4v(j) = d4v(j)*gkpr*gkap*dre

        d5v(j) = (-2./3.*xy(j,k-1)*yx(j,k-1)+xx(j,k-1)*yy(j,k-1))*rkm
        d5v(j) = d5v(j)+(-2./3.*xy(j,k+1)*yx(j,k+1)+xx(j,k+1)*
     1  yy(j,k+1))*rkp
        d5v(j) = vnu*dre*d5v(j)

        d7v(j) = (-2./3.*xx(j,k-1)*yy(j,k-1)+xy(j,k-1)*yx(j,k-1))*rkm
        d7v(j) = d7v(j)+(-2./3.*xx(j,k+1)*yy(j,k+1)+xy(j,k+1)*
     1  yx(j,k+1))*rkp
        d7v(j) = vnu*dre*d7v(j)

220     continue

	  do 230 j = 1,jm
          j1    = j+1

          a1    = a1v(j1)+a1v(j)
          a2    = a2v(j1)+a2v(j)
          a5    = a5v(j1)+a5v(j)
          a4    = a4v(j1)+a4v(j)
          d1    = d1v(j)
          d2    = d2v(j)
          d4    = d4v(j)
          d5    = d5v(j)
          d7    = d7v(j)

          ukp1  = q(j,k+1,2)/q(j,k+1,1)
          ukm1  = q(j,k-1,2)/q(j,k-1,1)
          duk   = 0.5*(ukp1-ukm1)

          vkp1  = q(j,k+1,3)/q(j,k+1,1)
          vkm1  = q(j,k-1,3)/q(j,k-1,1)
          dvk   = 0.5*(vkp1-vkm1)

          ekp1  = q(j,k+1,4)/q(j,k+1,1)-.5*(ukp1**2+vkp1**2)
          ekm1  = q(j,k-1,4)/q(j,k-1,1)-.5*(ukm1**2+vkm1**2)
          dek   = 0.5*(ekp1-ekm1)

          t4    = u(j1)*a1v(j1)+u(j)*a1v(j)
          t5    = v(j1)*a2v(j1)+v(j)*a2v(j)
          t6    = u(j1)*a5v(j1)+u(j)*a5v(j)
          t7    = v(j1)*a5v(j1)+v(j)*a5v(j)
          t8    = a4v(j1)+a4v(j)


          du    = u(j1)-u(j)
          dv    = v(j1)-v(j)
          dei   = e(j1)-e(j)
          f2(j) = a1*du+a5*dv
          f2v(j)= d1*duk+d7*dvk
          f3(j) = a5*du+a2*dv
          f3v(j)= d5*duk+d2*dvk
          f4(j) = (t4+t7)*du+(t5+t6)*dv+t8*dei
          f4v(j) = 0.5*d1*0.5*(ukp1**2-ukm1**2)
          f4v(j) = f4v(j)+0.5*d2*0.5*(vkp1**2-vkm1**2)
          f4v(j) = f4v(j)+d5*0.5*(vkp1+vkm1)*duk
          f4v(j) = f4v(j)+d7*0.5*(ukp1+ukm1)*dvk
          f4v(j) = f4v(j)+d4*dek
230     continue
c
        do 240 j = 2,jm-1
          s(j,k,2) = s(j,k,2)+(f2(j)-f2(j-1))+0.5*(f2v(j+1)-f2v(j-1))
          s(j,k,3) = s(j,k,3)+(f3(j)-f3(j-1))+0.5*(f3v(j+1)-f3v(j-1))
          s(j,k,4) = s(j,k,4)+(f4(j)-f4(j-1))+0.5*(f4v(j+1)-f4v(j-1))
240     continue
          j=jm
          s(j,k,2) = s(j,k,2)+(f2(j)-f2(j-1))+0.5*(f2v(j-2)
     1          -4.*f2v(j-1)+3*f2v(j))
          s(j,k,3) = s(j,k,3)+(f3(j)-f3(j-1))+0.5*(f3v(j-2)
     1          -4.*f3v(j-1)+3*f3v(j))
          s(j,k,4) = s(j,k,4)+(f4(j)-f4(j-1))+0.5*(f4v(j-2)
     1          -4.*f4v(j-1)+3*f4v(j))
210     continue

	endif

        elseif (ifourth.eq.1) then

        oneby3 = 1.d0/3.d0
        oneby6 = 1.d0/6.d0
        oneby8 = 1.d0/8.d0
        oneby12 = 1.d0/12.d0
        oneby24 = 1.d0/24.d0

        if(itlns.eq.1) then
        do 310 j = 2,jm
        do 320 k = 1,kmax
          ra    = 1./q(j,k,1)
          u(k)  = q(j,k,2)*ra
          v(k)  = q(j,k,3)*ra
          e(k)  = q(j,k,4)*ra-.5*(u(k)**2+v(k)**2)
          tt    = ggm1*e(k)
          vmue  = c2bp*tt*sqrt(tt)/(c2b + tt)
          turm  = turmu(j,k)
          vnu   = vmue+turm
          gkap  = vmue+prtr*turm
          rj    = 1./q(j,k,nq)
          b4v(k) = (yx(j,k)**2+yy(j,k)**2)*rj
          b1v(k) = (b4v(k)+yx(j,k)**2/3.*rj)*vnu*dre
          b2v(k) = (b4v(k)+yy(j,k)**2/3.*rj)*vnu*dre
          b5v(k) = (yx(j,k)*yy(j,k)/3.*rj)*vnu*dre
          b4v(k) = b4v(k)*gkpr*gkap*dre
  320   continue

        k     = 1
        k1    = k+1
 
        b1    = b1v(k1)+b1v(k)
        b2    = b2v(k1)+b2v(k)
        b5    = b5v(k1)+b5v(k)
        b4    = b4v(k1)+b4v(k)

        t4    = u(k1)*b1v(k1)+u(k)*b1v(k)
        t5    = v(k1)*b2v(k1)+v(k)*b2v(k)
        t6    = u(k1)*b5v(k1)+u(k)*b5v(k)
        t7    = v(k1)*b5v(k1)+v(k)*b5v(k)
        t8    = b4v(k1)+b4v(k)

        du    = u(k1)-u(k)
        dv    = v(k1)-v(k)
        dei   = e(k1)-e(k)

        f2(k) = b1*du+b5*dv
        f3(k) = b5*du+b2*dv
        f4(k) = (t4+t7)*du+(t5+t6)*dv+t8*dei

        do 330 k = 2,km-1
          km1   = k-1
          k1    = k+1
          k2    = k+2

          b1    = (-b1v(km1)+9*b1v(k)+9*b1v(k1)-b1v(k2))*oneby8
          b2    = (-b2v(km1)+9*b2v(k)+9*b2v(k1)-b2v(k2))*oneby8
          b5    = (-b5v(km1)+9*b5v(k)+9*b5v(k1)-b5v(k2))*oneby8
          b4    = (-b4v(km1)+9*b4v(k)+9*b4v(k1)-b4v(k2))*oneby8

          t4    = (-u(km1)*b1v(km1)+9*u(k)*b1v(k)+9*u(k1)*b1v(k1)
     1            -u(k2)*b1v(k2))*oneby8
          t5    = (-v(km1)*b2v(km1)+9*v(k)*b2v(k)+9*v(k1)*b2v(k1)
     1            -v(k2)*b2v(k2))*oneby8
          t6    = (-u(km1)*b5v(km1)+9*u(k)*b5v(k)+9*u(k1)*b5v(k1)
     1            -u(k2)*b5v(k2))*oneby8
          t7    = (-v(km1)*b5v(km1)+9*v(k)*b5v(k)+9*v(k1)*b5v(k1)
     1            -v(k2)*b5v(k2))*oneby8
          t8    = (-b4v(km1)+9*b4v(k)+9*b4v(k1)-b4v(k2))*oneby8

          du    = (u(km1)-27*u(k)+27*u(k1)-u(k2))*oneby24
          dv    = (v(km1)-27*v(k)+27*v(k1)-v(k2))*oneby24
          dei   = (e(km1)-27*e(k)+27*e(k1)-e(k2))*oneby24
          f2(k) = b1*du+b5*dv
          f3(k) = b5*du+b2*dv
          f4(k) = (t4+t7)*du+(t5+t6)*dv+t8*dei
  330   continue
        k     = km
        k1    = km+1
 
        b1    = b1v(K1)+b1v(K)
        b2    = b2v(K1)+b2v(K)
        b5    = b5v(K1)+b5v(K)
        b4    = b4v(K1)+b4v(K)

        t4    = u(k1)*b1v(k1)+u(k)*b1v(k)
        t5    = v(k1)*b2v(k1)+v(k)*b2v(k)
        t6    = u(k1)*b5v(k1)+u(k)*b5v(k)
        t7    = v(k1)*b5v(k1)+v(k)*b5v(k)
        t8    = b4v(k1)+b4v(k)

        du    = u(k1)-u(k)
        dv    = v(k1)-v(k)
        dei   = e(k1)-e(k)

        f2(k) = b1*du+b5*dv
        f3(k) = b5*du+b2*dv
        f4(k) = (t4+t7)*du+(t5+t6)*dv+t8*dei

c
        s(j,2,2) = s(j,2,2) + (f2(2)-f2(1))
        s(j,2,3) = s(j,2,3) + (f3(2)-f3(1))
        s(j,2,4) = s(j,2,4) + (f4(2)-f4(1))
        do 340 k = 3,km-1
          km2   = k-2
          km1   = k-1
          k1    = k+1

          s(j,k,2) = s(j,k,2) + (f2(km2)-27*f2(km1)
     1               +27*f2(k)-f2(k1))*oneby24
          s(j,k,3) = s(j,k,3) + (f3(km2)-27*f3(km1)
     1               +27*f3(k)-f3(k1))*oneby24
          s(j,k,4) = s(j,k,4) + (f4(km2)-27*f4(km1)
     1               +27*f4(k)-f4(k1))*oneby24
  340   continue
        s(j,km,2) = s(j,km,2) + (f2(km)-f2(km-1))
        s(j,km,3) = s(j,km,3) + (f3(km)-f3(km-1))
        s(j,km,4) = s(j,km,4) + (f4(km)-f4(km-1))

  310 continue

        else

        do 410 j = 2,jm
        do 420 k = 1,kmax
          ra    = 1./q(j,k,1)
          u(k)  = q(j,k,2)*ra
          v(k)  = q(j,k,3)*ra
          e(k)  = q(j,k,4)*ra-.5*(u(k)**2+v(k)**2)
          tt    = ggm1*e(k)
          vmue  = c2bp*tt*sqrt(tt)/(c2b + tt)
          turm  = turmu(j,k)
          vnu   = vmue+turm
          gkap  = vmue+prtr*turm
          rj    = 1./q(j,k,nq)
          rjm    = 1./q(j-1,k,nq)
          rjp    = 1./q(j+1,k,nq)

          b4v(k) = (yx(j,k)**2+yy(j,k)**2)*rj
          b1v(k) = (b4v(k)+yx(j,k)**2/3.*rj)*vnu*dre
          b2v(k) = (b4v(k)+yy(j,k)**2/3.*rj)*vnu*dre
          b5v(k) = (yx(j,k)*yy(j,k)/3.*rj)*vnu*dre
          b4v(k) = b4v(k)*gkpr*gkap*dre
          
          d1v(k) = (4./3.*xx(j-1,k)*yx(j-1,k)+xy(j-1,k)*yy(j-1,k))*rjm
          d1v(k) = d1v(k)+(4./3.*xx(j+1,k)*yx(j+1,k)+xy(j+1,k)*
     1             yy(j+1,k))*rjp
          d2v(k) = (4./3.*xy(j-1,k)*yy(j-1,k)+xx(j-1,k)*yx(j-1,k))*rjm
          d2v(k) = d2v(k)+(4./3.*xy(j+1,k)*yy(j+1,k)+xx(j+1,k)*
     1             yx(j+1,k))*rjp
          d4v(k) = (xy(j-1,k)*yy(j-1,k)+xx(j-1,k)*yx(j-1,k))*rjm
          d4v(k) = d4v(k)+(xy(j+1,k)*yy(j+1,k)+xx(j+1,k)*
     1             yx(j+1,k))*rjp
          d5v(k) = (-2./3.*xy(j-1,k)*yx(j-1,k)+xx(j-1,k)*yy(j-1,k))*rjm
          d5v(k) = d5v(k)+(-2./3.*xy(j+1,k)*yx(j+1,k)+xx(j+1,k)*
     1             yy(j+1,k))*rjp
          d7v(k) = (-2./3.*xx(j-1,k)*yy(j-1,k)+xy(j-1,k)*yx(j-1,k))*rjm
          d7v(k) = d7v(k)+(-2./3.*xx(j+1,k)*yy(j+1,k)+xy(j+1,k)*
     1             yx(j+1,k))*rjp

          if (j.gt.2 .and. j.lt.jm) then

            rjm2   = 1./q(j-2,k,nq)
            rjp2   = 1./q(j+2,k,nq)

            d1v(k) = 4*d1v(k)-(4./3.*xx(j-2,k)*yx(j-2,k)+xy(j-2,k)*
     1               yy(j-2,k))*rjm2
            d1v(k) = d1v(k)-(4./3.*xx(j+2,k)*yx(j+2,k)+xy(j+2,k)*
     1               yy(j+2,k))*rjp2
            d1v(k) = d1v(k)*oneby3

            d2v(k) = 4*d2v(k)-(4./3.*xy(j-2,k)*yy(j-2,k)+xx(j-2,k)*
     1               yx(j-2,k))*rjm2
            d2v(k) = d2v(k)-(4./3.*xy(j+2,k)*yy(j+2,k)+xx(j+2,k)*
     1               yx(j+2,k))*rjp2
            d2v(k) = d2v(k)*oneby3

            d4v(k) = 4*d4v(k)-(xy(j-2,k)*yy(j-2,k)+xx(j-2,k)*
     1               yx(j-2,k))*rjm2
            d4v(k) = d4v(k)-(xy(j+2,k)*yy(j+2,k)+xx(j+2,k)*
     1               yx(j+2,k))*rjp2
            d4v(k) = d4v(k)*oneby3

            d5v(k) = 4*d5v(k)-(-2./3.*xy(j-2,k)*yx(j-2,k)+xx(j-2,k)*
     1               yy(j-2,k))*rjm2
            d5v(k) = d5v(k)-(-2./3.*xy(j+2,k)*yx(j+2,k)+xx(j+2,k)*
     1               yy(j+2,k))*rjp2
            d5v(k) = d5v(k)*oneby3

            d7v(k) = 4*d7v(k)-(-2./3.*xx(j-2,k)*yy(j-2,k)+xy(j-2,k)*
     1               yx(j-2,k))*rjm2
            d7v(k) = d7v(k)-(-2./3.*xx(j+2,k)*yy(j+2,k)+xy(j+2,k)*
     1               yx(j+2,k))*rjp2
            d7v(k) = d7v(k)*oneby3

          endif 

          d1v(k) = vnu*dre*d1v(k)
          d2v(k) = vnu*dre*d2v(k)
          d4v(k) = d4v(k)*gkpr*gkap*dre
          d5v(k) = vnu*dre*d5v(k)
          d7v(k) = vnu*dre*d7v(k)

420     continue

        do 430 k = 1,km
          d1    = d1v(k)
          d2    = d2v(k)
          d4    = d4v(k)
          d5    = d5v(k)
          d7    = d7v(k)

          ujp1  = q(j+1,k,2)/q(j+1,k,1)
          ujm1  = q(j-1,k,2)/q(j-1,k,1)

          vjp1  = Q(J+1,K,3)/Q(J+1,K,1)
          vjm1  = Q(J-1,K,3)/Q(J-1,K,1)

          ejp1  = q(j+1,k,4)/q(j+1,k,1)-.5*(ujp1**2+vjp1**2)
          ejm1  = q(j-1,k,4)/q(j-1,k,1)-.5*(ujm1**2+vjm1**2)

          if (j.eq.2 .or. j.eq.jm) then

            duj   = 0.5*(ujp1-ujm1)
            dvj   = 0.5*(vjp1-vjm1)
            dej   = 0.5*(ejp1-ejm1)
            
            f4v(k) = 0.5*d1*0.5*(ujp1**2-ujm1**2)
            f4v(k) = f4v(k)+0.5*d2*0.5*(vjp1**2-vjm1**2)
            f4v(k) = f4v(k)+d5*0.5*(ujp1+ujm1)*dvj
            f4v(k) = f4v(k)+d7*0.5*(vjp1+vjm1)*duj

           
          else

          ujp2  = q(j+2,k,2)/q(j+2,k,1)
          ujm2  = q(j-2,k,2)/q(j-2,k,1)
          duj   = (-ujp2+8*ujp1-8*ujm1+ujm2)*oneby12

          vjp2  = q(j+2,k,3)/q(j+2,k,1)
          vjm2  = q(j-2,k,3)/q(j-2,k,1)
          dvj   = (-vjp2+8*vjp1-8*vjm1+vjm2)*oneby12

          ejp2  = q(j+2,k,4)/q(j+2,k,1)-.5*(ujp2**2+vjp2**2)
          ejm2  = q(j-2,k,4)/q(j-2,k,1)-.5*(ujm2**2+vjm2**2)
          dej   = (-ejp2+8*ejp1-8*ejm1+ejm2)*oneby12

          f4v(k) = 0.5*d1*(-ujp2**2+8*ujp1**2-8*ujm1**2+ujm2**2)*oneby12
          f4v(k) = f4v(k)+0.5*d2*(-vjp2**2+8*vjp1**2
     1            -8*vjm1**2+vjm2**2)*oneby12
          f4v(k) = f4v(k)+d5*oneby6*(-ujm2+4*ujm1+4*ujp1-ujp2)*dvj
          f4v(k) = f4v(k)+d7*oneby6*(-vjm2+4*vjm1+4*vjp1-vjp2)*duj

          endif

          f2v(k) = d1*duj+d5*dvj
          f3v(k) = d7*duj+d2*dvj
          f4v(k) = f4v(k)+d4*dej

430     continue

          k   =  1
          k1    = k+1

          b1    = b1v(k1)+b1v(k)
          b2    = b2v(k1)+b2v(k)
          b5    = b5v(k1)+b5v(k)
          b4    = b4v(k1)+b4v(k)

          t4    = u(k1)*b1v(k1)+u(k)*b1v(k)
          t5    = v(k1)*b2v(k1)+v(k)*b2v(k)
          t6    = u(k1)*b5v(k1)+u(k)*b5v(k)
          t7    = v(k1)*b5v(k1)+v(k)*b5v(k)
          t8    = b4v(k1)+b4v(k)

          du     = u(k1)-u(k)
          dv     = v(k1)-v(k)
          dei    = e(k1)-e(k)
          f2(k)  = b1*du+b5*dv
          f3(k)  = b5*du+b2*dv
          f4(k)  = (t4+t7)*du+(t5+t6)*dv+t8*dei

        do 440 k = 2,km-1
          km1   = k-1
          k1    = k+1
          k2    = k+2

          b1    = (-b1v(km1)+9*b1v(k)+9*b1v(k1)-b1v(k2))*oneby8
          b2    = (-b2v(km1)+9*b2v(k)+9*b2v(k1)-b2v(k2))*oneby8
          b5    = (-b5v(km1)+9*b5v(k)+9*b5v(k1)-b5v(k2))*oneby8
          b4    = (-b4v(km1)+9*b4v(k)+9*b4v(k1)-b4v(k2))*oneby8

          t4    = (-u(km1)*b1v(km1)+9*u(k)*b1v(k)+9*u(k1)*b1v(k1)
     1            -u(k2)*b1v(k2))*oneby8
          t5    = (-v(km1)*b2v(km1)+9*v(k)*b2v(k)+9*v(k1)*b2v(k1)
     1            -v(k2)*b2v(k2))*oneby8
          t6    = (-u(km1)*b5v(km1)+9*u(k)*b5v(k)+9*u(k1)*b5v(k1)
     1            -u(k2)*b5v(k2))*oneby8
          t7    = (-v(km1)*b5v(km1)+9*v(k)*b5v(k)+9*v(k1)*b5v(k1)
     1            -v(k2)*b5v(k2))*oneby8
          t8    = (-b4v(km1)+9*b4v(k)+9*b4v(k1)-b4v(k2))*oneby8

          du    = (u(km1)-27*u(k)+27*u(k1)-u(k2))*oneby24
          dv    = (v(km1)-27*v(k)+27*v(k1)-v(k2))*oneby24
          dei   = (e(km1)-27*e(k)+27*e(k1)-e(k2))*oneby24

          f2(k) = b1*du+b5*dv
          f3(k) = b5*du+b2*dv
          f4(k) = (t4+t7)*du+(t5+t6)*dv+t8*dei
440     continue
        k   =  km
        k1    = k+1

        b1    = b1v(k1)+b1v(k)
        b2    = b2v(k1)+b2v(k)
        b5    = b5v(k1)+b5v(k)
        b4    = b4v(k1)+b4v(k)

        t4    = u(k1)*b1v(k1)+u(k)*b1v(k)
        t5    = v(k1)*b2v(k1)+v(k)*b2v(k)
        t6    = u(k1)*b5v(k1)+u(k)*b5v(k)
        t7    = v(k1)*b5v(k1)+v(k)*b5v(k)
        t8    = b4v(k1)+b4v(k)

        du     = u(k1)-u(k)
        dv     = v(k1)-v(k)
        dei    = e(k1)-e(k)
        f2(k)  = b1*du+b5*dv
        f3(k)  = b5*du+b2*dv
        f4(k)  = (t4+t7)*du+(t5+t6)*dv+t8*dei

c
        k = 2 
        s(j,k,2) = s(j,k,2)+(f2(k)-f2(k-1))+0.5*(f2v(k+1)-f2v(k-1))
        s(j,k,3) = s(j,k,3)+(f3(k)-f3(k-1))+0.5*(f3v(k+1)-f3v(k-1))
        s(j,k,4) = s(j,k,4)+(f4(k)-f4(k-1))+0.5*(f4v(k+1)-f4v(k-1))

        do 450 k = 3,km-2
          km2   = k-2
          km1   = k-1
          k1    = k+1
          k2    = k+2

          s(j,k,2) = s(j,k,2) + (f2(km2)-27*f2(km1)
     1                       +27*f2(k)-f2(k1))*oneby24
     1            +(f2v(km2)-8*f2v(km1)+8*f2v(k1)-f2v(k2))*oneby12
          s(j,k,3) = s(j,k,3) + (f3(km2)-27*f3(km1)
     1                       +27*f3(k)-f3(k1))*oneby24
     1            +(f3v(km2)-8*f3v(km1)+8*f3v(k1)-f3v(k2))*oneby12
          s(j,k,4) = s(j,k,4) + (f4(km2)-27*f4(km1)
     1                       +27*f4(k)-f4(k1))*oneby24
     1            +(f4v(km2)-8*f4v(km1)+8*f4v(k1)-f4v(k2))*oneby12

450     continue

        k = km-1 
        km2   = k-2
        km1   = k-1
        k1    = k+1
        s(j,k,2) = s(j,k,2) + (f2(km2)-27*f2(km1)
     1                      + 27*f2(k)-f2(k1))*oneby24
     1                      + 0.5*(f2v(k1)-f2v(km1))
        s(j,k,3) = s(j,k,3) + (f3(km2)-27*f3(km1)
     1                      + 27*f3(k)-f3(k1))*oneby24
     1                      + 0.5*(f3v(k1)-f3v(km1))
        s(j,k,4) = s(j,k,4) + (f4(km2)-27*f4(km1)
     1                      + 27*f4(k)-f4(k1))*oneby24
     1                      + 0.5*(f4v(k1)-f4v(km1))

        k=km
        s(j,k,2) = s(j,k,2)+(f2(k)-f2(k-1))+0.5*(f2v(k-2)
     1            - 4.*f2v(k-1)+3*f2v(k))
        s(j,k,3) = s(j,k,3)+(f3(k)-f3(k-1))+0.5*(f3v(k-2)
     1            - 4.*f3v(k-1)+3*f3v(k))
        s(j,k,4) = s(j,k,4)+(f4(k)-f4(k-1))+0.5*(f4v(k-2)
     1            - 4.*f4v(k-1)+3.*f4v(k))


410   continue

      do 510 k = 2,km
        do 520 j = 1,jmax
          ra    = 1./q(j,k,1)
          u(j)  = q(j,k,2)*ra
          v(j)  = q(j,k,3)*ra
          e(j)  = q(j,k,4)*ra-.5*(u(j)**2+v(j)**2)
          tt    = ggm1*e(j)
          vmue  = c2bp*tt*sqrt(tt)/(c2b + tt)
          turm  = turmu(j,k)
          vnu   = vmue+turm
          gkap  = vmue+prtr*turm
          rk    = 1./q(j,k,nq)
          rkm    = 1./q(j,k-1,nq)
          rkp    = 1./q(j,k+1,nq)

          a4v(j) = (xx(j,k)**2+xy(j,k)**2)*rk
          a1v(j) = (a4v(j)+xx(j,k)**2/3.*rk)*vnu*dre
          a2v(j) = (a4v(j)+xy(j,k)**2/3.*rk)*vnu*dre
          a5v(j) = (xx(j,k)*xy(j,k)/3.*rk)*vnu*dre
          a4v(j) = a4v(j)*gkpr*gkap*dre

          d1v(j) = (4./3.*xx(j,k-1)*yx(j,k-1)+xy(j,k-1)*yy(j,k-1))*rkm
          d1v(j) = d1v(j)+(4./3.*xx(j,k+1)*yx(j,k+1)+xy(j,k+1)*
     1    yy(j,k+1))*rkp

          d2v(j) = (4./3.*xy(j,k-1)*yy(j,k-1)+xx(j,k-1)*yx(j,k-1))*rkm
          d2v(j) = d2v(j)+(4./3.*xy(j,k+1)*yy(j,k+1)+xx(j,k+1)*
     1    yx(j,k+1))*rkp

          d4v(j) = (xy(j,k-1)*yy(j,k-1)+xx(j,k-1)*yx(j,k-1))*rkm
          d4v(j) = d4v(j)+(xy(j,k+1)*yy(j,k+1)+xx(j,k+1)*
     1    yx(j,k+1))*rkp

          d5v(j) = (-2./3.*xy(j,k-1)*yx(j,k-1)+xx(j,k-1)*yy(j,k-1))*rkm
          d5v(j) = d5v(j)+(-2./3.*xy(j,k+1)*yx(j,k+1)+xx(j,k+1)*
     1    yy(j,k+1))*rkp

          d7v(j) = (-2./3.*xx(j,k-1)*yy(j,k-1)+xy(j,k-1)*yx(j,k-1))*rkm
          d7v(j) = d7v(j)+(-2./3.*xx(j,k+1)*yy(j,k+1)+xy(j,k+1)*
     1    yx(j,k+1))*rkp

          if (k.gt.2 .and. k.lt.km) then

            rkm2   = 1./q(j,k-2,nq)
            rkp2   = 1./q(j,k+2,nq)

            d1v(j) = 4*d1v(j)-(4./3.*xx(j,k-2)*yx(j,k-2)+xy(j,k-2)*
     1              yy(j,k-2))*rkm2
            d1v(j) = d1v(j)-(4./3.*xx(j,k+2)*yx(j,k+2)+xy(j,k+2)*
     1              yy(j,k+2))*rkp2
            d1v(j) = d1v(j)*oneby3

            d2v(j) = 4*d2v(j)-(4./3.*xy(j,k-2)*yy(j,k-2)+xx(j,k-2)*
     1              yx(j,k-2))*rkm2
            d2v(j) = d2v(j)-(4./3.*xy(j,k+2)*yy(j,k+2)+xx(j,k+2)*
     1              yx(j,k+2))*rkp2
            d2v(j) = d2v(j)*oneby3

            d4v(j) = 4*d4v(j)-(xy(j,k-2)*yy(j,k-2)+xx(j,k-2)*
     1              yx(j,k-2))*rkm2
            d4v(j) = d4v(j)-(xy(j,k+2)*yy(j,k+2)+xx(j,k+2)*
     1              yx(j,k+2))*rkp2
            d4v(j) = d4v(j)*oneby3

            d5v(j) = 4*d5v(j)-(-2./3.*xy(j,k-2)*yx(j,k-2)+xx(j,k-2)*
     1              yy(j,k-2))*rkm2
            d5v(j) = d5v(j)-(-2./3.*xy(j,k+2)*yx(j,k+2)+xx(j,k+2)*
     1              yy(j,k+2))*rkp2
            d5v(j) = d5v(j)*oneby3

            d7v(j) = 4*d7v(j)-(-2./3.*xx(j,k-2)*yy(j,k-2)+xy(j,k-2)*
     1              yx(j,k-2))*rkm2
            d7v(j) = d7v(j)-(-2./3.*xx(j,k+2)*yy(j,k+2)+xy(j,k+2)*
     1              yx(j,k+2))*rkp2
            d7v(j) = d7v(j)*oneby3

          endif 

          d1v(j) = vnu*dre*d1v(j)
          d2v(j) = vnu*dre*d2v(j)
          d4v(j) = d4v(j)*gkpr*gkap*dre
          d5v(j) = vnu*dre*d5v(j)
          d7v(j) = vnu*dre*d7v(j)


520     continue

        do 530 j = 1,jm
          d1    = d1v(j)
          d2    = d2v(j)
          d4    = d4v(j)
          d5    = d5v(j)
          d7    = d7v(j)


          ukp1  = q(j,k+1,2)/q(j,k+1,1)
          ukm1  = q(j,k-1,2)/q(j,k-1,1)

          vkp1  = q(j,k+1,3)/q(j,k+1,1)
          vkm1  = q(j,k-1,3)/q(j,k-1,1)

          ekp1  = q(j,k+1,4)/q(j,k+1,1)-.5*(ukp1**2+vkp1**2)
          ekm1  = q(j,k-1,4)/q(j,k-1,1)-.5*(ukm1**2+vkm1**2)

          if (k.eq.2 .or. k.eq.km) then

            duk   = 0.5*(ukp1-ukm1)
            dvk   = 0.5*(vkp1-vkm1)
            dek   = 0.5*(ekp1-ekm1)
            
            f4v(j) = 0.5*d1*0.5*(ukp1**2-ukm1**2)
            f4v(j) = f4v(j)+0.5*d2*0.5*(vkp1**2-vkm1**2)
            f4v(j) = f4v(j)+d5*0.5*(vkp1+vkm1)*duk
            f4v(j) = f4v(j)+d7*0.5*(ukp1+ukm1)*dvk

           
          else

            ukp2  = q(j,k+2,2)/q(j,k+2,1)
            ukm2  = q(j,k-2,2)/q(j,k-2,1)
            duk   = (-ukp2+8*ukp1-8*ukm1+ukm2)*oneby12

            vkp2  = q(j,k+2,3)/q(j,k+2,1)
            vkm2  = q(j,k-2,3)/q(j,k-2,1)
            dvk   = (-vkp2+8*vkp1-8*vkm1+vkm2)*oneby12

            ekp2  = q(j,k+2,4)/q(j,k+2,1)-.5*(ukp2**2+vkp2**2)
            ekm2  = q(j,k-2,4)/q(j,k-2,1)-.5*(ukm2**2+vkm2**2)
            dek   = (-ekp2+8*ekp1-8*ekm1+ekm2)*oneby12

	    f4v(j) = 0.5*d1*(-ukp2**2+8*ukp1**2-8*ukm1**2+ukm2**2)*oneby12
            f4v(j) = f4v(j)+0.5*d2*(-vkp2**2+8*vkp1**2
     1              -8*vkm1**2+vkm2**2)*oneby12
            f4v(j) = f4v(j)+d5*oneby6*(-vkm2+4*vkm1+4*vkp1-vkp2)*duk
            f4v(j) = f4v(j)+d7*oneby6*(-ukm2+4*ukm1+4*ukp1-ukp2)*dvk

          endif

          f2v(j)= d1*duk+d7*dvk
          f3v(j)= d5*duk+d2*dvk
          f4v(j) = f4v(j)+d4*dek

530     continue
          j     = 1
          j1    = j+1

          a1    = a1v(j1)+a1v(j)
          a2    = a2v(j1)+a2v(j)
          a5    = a5v(j1)+a5v(j)
          a4    = a4v(j1)+a4v(j)

          t4    = u(j1)*a1v(j1)+u(j)*a1v(j)
          t5    = v(j1)*a2v(j1)+v(j)*a2v(j)
          t6    = u(j1)*a5v(j1)+u(j)*a5v(j)
          t7    = v(j1)*a5v(j1)+v(j)*a5v(j)
          t8    = a4v(j1)+a4v(j)

          du    = u(j1)-u(j)
          dv    = v(j1)-v(j)
          dei   = e(j1)-e(j)

          f2(j) = a1*du+a5*dv
          f3(j) = a5*du+a2*dv
          f4(j) = (t4+t7)*du+(t5+t6)*dv+t8*dei
        do 540 j = 2,jm-1
          jm1   = j-1
          j1    = j+1
          j2    = j+2

          a1    = (-a1v(jm1)+9*a1v(j)+9*a1v(j1)-a1v(j2))*oneby8
          a2    = (-a2v(jm1)+9*a2v(j)+9*a2v(j1)-a2v(j2))*oneby8
          a5    = (-a5v(jm1)+9*a5v(j)+9*a5v(j1)-a5v(j2))*oneby8
          a4    = (-a4v(jm1)+9*a4v(j)+9*a4v(j1)-a4v(j2))*oneby8

          t4    = (-u(jm1)*a1v(jm1)+9*u(j)*a1v(j)+9*u(j1)*a1v(j1)
     1           -u(j2)*a1v(j2))*oneby8
          t5    = (-v(jm1)*a2v(jm1)+9*v(j)*a2v(j)+9*v(j1)*a2v(j1)
     1            -v(j2)*a2v(j2))*oneby8
          t6    = (-u(jm1)*a5v(jm1)+9*u(j)*a5v(j)+9*u(j1)*a5v(j1)
     1            -u(j2)*a5v(j2))*oneby8
          t7    = (-v(jm1)*a5v(jm1)+9*v(j)*a5v(j)+9*v(j1)*a5v(j1)
     1            -v(j2)*a5v(j2))*oneby8
          t8    = (-a4v(jm1)+9*a4v(j)+9*a4v(j1)-a4v(j2))*oneby8

          du    = (u(jm1)-27*u(j)+27*u(j1)-u(j2))*oneby24
          dv    = (v(jm1)-27*v(j)+27*v(j1)-v(j2))*oneby24
          dei   = (e(jm1)-27*e(j)+27*e(j1)-e(j2))*oneby24

          f2(j) = a1*du+a5*dv
          f3(j) = a5*du+a2*dv
          f4(j) = (t4+t7)*du+(t5+t6)*dv+t8*dei
540     continue
        j     = jm
        j1    = j+1

        a1    = a1v(j1)+a1v(j)
        a2    = a2v(j1)+a2v(j)
        a5    = a5v(j1)+a5v(j)
        a4    = a4v(j1)+a4v(j)

        t4    = u(j1)*a1v(j1)+u(j)*a1v(j)
        t5    = v(j1)*a2v(j1)+v(j)*a2v(j)
        t6    = u(j1)*a5v(j1)+u(j)*a5v(j)
        t7    = v(j1)*a5v(j1)+v(j)*a5v(j)
        t8    = a4v(j1)+a4v(j)

        du    = u(j1)-u(j)
        dv    = v(j1)-v(j)
        dei   = e(j1)-e(j)

        f2(j) = a1*du+a5*dv
        f3(j) = a5*du+a2*dv
        f4(j) = (t4+t7)*du+(t5+t6)*dv+t8*dei


        j = 2
        s(j,k,2) = s(j,k,2)+(f2(j)-f2(j-1))+0.5*(f2v(j+1)-f2v(j-1))
        s(j,k,3) = s(j,k,3)+(f3(j)-f3(j-1))+0.5*(f3v(j+1)-f3v(j-1))
        s(j,k,4) = s(j,k,4)+(f4(j)-f4(j-1))+0.5*(f4v(j+1)-f4v(j-1))
        do 550 j = 3,jm-2
          jm2   = j-2
          jm1   = j-1
          j1    = j+1
          j2    = j+2
          s(j,k,2) = s(j,k,2) + (f2(jm2)-27*f2(jm1)
     1                       +27*f2(j)-f2(j1))*oneby24
     1              +(f2v(jm2)-8*f2v(jm1)+8*f2v(j1)-f2v(j2))*oneby12
          s(j,k,3) = s(j,k,3) + (f3(jm2)-27*f3(jm1)
     1                       +27*f3(j)-f3(j1))*oneby24
     1              +(f3v(jm2)-8*f3v(jm1)+8*f3v(j1)-f3v(j2))*oneby12
          s(j,k,4) = s(j,k,4) + (f4(jm2)-27*f4(jm1)
     1                       +27*f4(j)-f4(j1))*oneby24
     1              +(f4v(jm2)-8*f4v(jm1)+8*f4v(j1)-f4v(j2))*oneby12


550     continue
        j=jm-1
        jm2   = j-2
        jm1   = j-1
        j1    = j+1
        s(j,k,2) = s(j,k,2) + (f2(jm2)-27*f2(jm1)
     1                      + 27*f2(j)-f2(j1))*oneby24
     1                      + 0.5*(f2v(j+1)-f2v(j-1))
        s(j,k,3) = s(j,k,3) + (f3(jm2)-27*f3(jm1)
     1                      + 27*f3(j)-f3(j1))*oneby24
     1                      + 0.5*(f3v(j+1)-f3v(j-1))
        s(j,k,4) = s(j,k,4) + (f4(jm2)-27*f4(jm1)
     1                      + 27*f4(j)-f4(j1))*oneby24
     1                      + 0.5*(f4v(j+1)-f4v(j-1))

        j=jm
        s(j,k,2) = s(j,k,2)+(f2(j)-f2(j-1))+0.5*(f2v(j-2)
     1        -4.*f2v(j-1)+3*f2v(j))
        s(j,k,3) = s(j,k,3)+(f3(j)-f3(j-1))+0.5*(f3v(j-2)
     1        -4.*f3v(j-1)+3*f3v(j))
        s(j,k,4) = s(j,k,4)+(f4(j)-f4(j-1))+0.5*(f4v(j-2)
     1        -4.*f4v(j-1)+3*f4v(j))
510     continue

      endif

      endif


 	return
      end



c***********************************************************************
      subroutine vmutur(x,y,q,s,turmu,xx,xy,yx,yy,ug,vg,jd,kd)
c
c  turbulent eddy viscosity. model is algebraic model of baldwin
c  and lomax
c
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************
      integer jd,kd
      real q(jd,kd,nq), s(jd,kd,nv), turmu(jd,kd)
      real x(jd,kd), y(jd,kd)
      real xx(jd,kd), xy(jd,kd), yx(jd,kd), yy(jd,kd)
      real ug(jd,kd), vg(jd,kd)

      ! local variables
      real,allocatable :: fmax(:),smax(:),vmax(:),vmin(:),grdnrm(:),
     <                    ra(:),snm(:),snp(:),vrtm(:),vrtp(:)

      real,allocatable :: sn(:,:),vort(:,:)
      real,allocatable :: qu(:,:),qv(:,:)

      integer kedge,kedge1,j1,j2,k,j
      real dre,duj,vjp1,vjm1,ejp1,dx2,dy2
      real c1,c2,c3,c4,c5,re,wmu,ux,vx,uy,vy,tx,ty
      real rhomuw,uu,dx,dy,fl,flp,flm,dfm,dfp,dsm,dsp
      real am,bm,si,fli,dvi2,t1,t2,fwake,t3,flkeb,slen,tmi
      real dv2,fkleb

      allocate(fmax(mdim),smax(mdim),vmax(mdim),vmin(mdim),
     <          grdnrm(mdim),ra(mdim),snm(mdim),snp(mdim),
     <          vrtm(mdim),vrtp(mdim))

      allocate(sn(jd,kd),vort(jd,kd))
      allocate(qu(jd,kd),qv(jd,kd))

c***********************************************************************
      data c1,c2,c3,c4,c5/0.4,26.0,0.01688,1.00,0.3/
c***********************************************************************
c
      re = rey
c
c..compute vorticity and normal distances from wall for entire flowfield
c
      kedge   = (3*kmax)/4
      kedge1  = kedge + 1
      j1      = jtail1
      j2      = jtail2
      dx2     = 0.5
      dy2     = 0.5
c
c..zero out turmu
c
      do 12 k = 1,kmax
      do 12 j = 1,jmax
        turmu(j,k) = 0.0
   12 continue
c
c..set wall molecular viscosity
c
      wmu = 1.0
c
cgrs..for compatibility with turbulence model
c
        do 120 k = 1,kedge1+1
        do 120 j = j1-1,j2+1
          qu(j,k)  = q(j,k,2)  - q(j,k,1)*ug(j,k)
          qv(j,k)  = q(j,k,3)  - q(j,k,1)*vg(j,k)
 120    continue
c
        k   = 1
        do 1 j = j1,j2
c
          ux  = (qu(j+1,k)/q(j+1,k,1) - qu(j-1,k)/q(j-1,k,1)) * dx2
          vx  = (qv(j+1,k)/q(j+1,k,1) - qv(j-1,k)/q(j-1,k,1)) * dx2
          uy  = -(3.0*qu(j,k)/q(j,k,1) - 4.0*qu(j,k+1)/q(j,k+1,1) +
     &                qu(j,k+2)/q(j,k+2,1))*dy2
          vy  = -(3.0*qv(j,k)/q(j,k,1) - 4.0*qv(j,k+1)/q(j,k+1,1) +
     &                qv(j,k+2)/q(j,k+2,1))*dy2
c     
          tx  =  xy(j,k)*ux -xx(j,k)*vx +yy(j,k)*uy -yx(j,k)*vy
          ty  =  xx(j,k)*vx -xy(j,k)*ux +yx(j,k)*vy -yy(j,k)*uy
c
          vort(j,k) = sqrt(tx**2 +ty**2)
          sn(j,k)  = 0.0
    1   continue
        k   = 1
        do 2 j = j1,j2
          rhomuw  = q(j,k,1)*q(j,k,nq)/wmu
          ra(j)   = sqrt( re*rhomuw*vort(j,k) )/c2
c
cgrs..for moving blades
cgrs..   uu = sqrt(uwal**2 + vwal**2)/q(j,k,1)
c
          uu = sqrt(qu(j,k)**2 +qv(j,k)**2)/q(j,k,1)
c
          fmax(j)  = 1.e-3
          vmax(j)  = uu
          vmin(j)  = uu
          vrtp(j)  = 0.0
          vrtm(j)  = 0.0
          grdnrm(j) = sqrt(yx(j,k)**2 +yy(j,k)**2)
    2   continue
c
        do 3 k = 2,kedge1
        do 3 j = j1,j2
          ux  = (qu(j+1,k)/q(j+1,k,1) - qu(j-1,k)/q(j-1,k,1)) * dx2
          vx  = (qv(j+1,k)/q(j+1,k,1) - qv(j-1,k)/q(j-1,k,1)) * dx2
          uy  = -(3.0*qu(j,k)/q(j,k,1) - 4.0*qu(j,k+1)/q(j,k+1,1) +
     &                qu(j,k+2)/q(j,k+2,1))*dy2
          vy  = -(3.0*qv(j,k)/q(j,k,1) - 4.0*qv(j,k+1)/q(j,k+1,1) +
     &                qv(j,k+2)/q(j,k+2,1))*dy2
c     
          tx  =  xy(j,k)*ux -xx(j,k)*vx +yy(j,k)*uy -yx(j,k)*vy
          ty  =  xx(j,k)*vx -xy(j,k)*ux +yx(j,k)*vy -yy(j,k)*uy
c
          vort(j,k) = sqrt(tx**2 +ty**2)
c
          dx  = x(j,k) -x(j,1)
          dy  = y(j,k) -y(j,1)
c
          sn(j,k)  = (yx(j,1)*dx +yy(j,1)*dy)/grdnrm(j)
    3   continue
        k   = 2
        do 4 j = j1,j2
          snp(j)   = sn(j,k+1)
          smax(j)  = sn(j,k)
          snm(j)   = sn(j,k-1)
    4   continue
c
c..compute fmax, smax, vmax, vmin
c
        do 11 k = 2,kedge
        do 11 j = j1,j2
c
          fl       = sn(j,k)*vort(j,k)*(1.0 -exp(-ra(j)*sn(j,k)))
          fmax(j)  = amax1(fmax(j),fl)
ccray          smax(j)  = cvmgt(sn(j,k),smax(j),fmax(j).eq.fl)
ccray          snp(j)   = cvmgt(sn(j,k+1),snp(j),fmax(j).eq.fl)
ccray          snm(j)   = cvmgt(sn(j,k-1),snm(j),fmax(j).eq.fl)
ccray          vrtp(j)  = cvmgt(vort(j,k+1),vrtp(j),fmax(j).eq.fl)
ccray          vrtm(j)  = cvmgt(vort(j,k-1),vrtm(j),fmax(j).eq.fl)
          if( fmax(j).eq.fl ) then
            smax(j)  = sn(j,k)
            snp(j)   = sn(j,k+1)
            snm(j)   = sn(j,k-1)
            vrtp(j)  = vort(j,k+1)
            vrtm(j)  = vort(j,k-1)
          end if
c
          uu = sqrt(qu(j,k)**2 +qv(j,k)**2)/q(j,k,1)
c
cgrs..for moving blades
cgrs..    uu = sqrt(uwal**2 + vwal**2)/q(j,k,1)
c
          vmax(j)  = amax1(uu,vmax(j))
          vmin(j)  = amin1(uu,vmin(j))
   11   continue
c
c..interpolation to improve estimate of fmax and smax
c
        do 21 j = j1,j2
          flp = snp(j)*vrtp(j)*(1.0 -exp(-ra(j)*snp(j)))
          flm = snm(j)*vrtm(j)*(1.0 -exp(-ra(j)*snm(j)))
          flp = amax1(flp,1.e-3)
          flm = amax1(flm,1.e-3)
          dfm = fmax(j) -flm
          dfp = flp -fmax(j)
          dsm = smax(j) -snm(j)
          dsp = snp(j) -smax(j)
          am  = (dsp**2*dfm +dsm**2*dfp)/(dsp*dsm*(dsp+dsm))
          bm  = (dsm*dfp -dsp*dfm)/(dsp*dsm*(dsp+dsm))
ccray            si  = cvmgt(smax(j) -0.5*am/(bm+1.e-21),smax(j),bm.lt.0.0 )
ccray            fli = cvmgt(fmax(j) -0.25*am**2/(bm+1.e-21),fmax(j),bm.lt.0.0)
          if( bm.lt.0.0 ) then
            si  = smax(j) -0.5*am/(bm+1.e-21)
            fli = fmax(j) -0.25*am**2/(bm+1.e-21)
          else
            si  = smax(j)
            fli = fmax(j)
          end if
ccray          fmax(j)  = cvmgt(fmax(j),fli,fli.lt.fmax(j)
ccray     *                       .or. si.lt.snm(j) .or. si.gt.snp(j))
ccray          smax(j)  = cvmgt(smax(j),si,fli.lt.fmax(j)
ccray     *                      .or. si.lt.snm(j) .or. si.gt.snp(j))
          if( fli.lt.fmax(j).or. si.lt.snm(j) .or. si.gt.snp(j) ) then
            fmax(j)  = fmax(j)
            smax(j)  = smax(j)
          else
            fmax(j)  = fli
            smax(j)  = si
          end if
   21   continue
c
c..compute outer eddy viscosity
c
        do 31 k = 2,kedge
        do 31 j = j1,j2
          dv2   = (vmax(j) -vmin(j))**2
          t1    = smax(j)*fmax(j)
          t2    = c4*smax(j)*dv2/fmax(j)
          fwake = amin1(t1,t2)
          t3    = (c5*sn(j,k)/smax(j))
          fkleb = 1.0 +5.5*amin1(1.e5,t3)**6
          turmu(j,k) = c3*q(j,k,1)*q(j,k,nq)*fwake/fkleb
   31   continue
c
c..compute inner eddy viscosity and set final eddy viscosity
c
        do 41 k = 2,kedge
        do 41 j = j1,j2
          slen  = c1*sn(j,k)*(1.0 -exp(-ra(j)*sn(j,k)))
          tmi   = q(j,k,1)*q(j,k,nq)*slen**2*vort(j,k)
          turmu(j,k) = amin1(turmu(j,k),tmi)*re
   41   continue

      return
      end


c*************************************************************************
      subroutine do_interpolations(qg,vnutg,jmx,kmx,ibcg,imeshg,idonorg,
     c        fracg,nfringeg,ndonorg,iisptrg,iieptrg,idsize,qsize,nmesh)
c***********************************************************************
      use params_global
c***********************************************************************
      implicit none
c***********************************************************************
      integer idsize,qsize,nmesh
      integer jmx(nmesh),kmx(nmesh)
      integer nfringeg(nmesh),ndonorg(nmesh),iisptrg(nmesh),iieptrg(nmesh)
      real qg(qsize), vnutg(idsize)
      integer imeshg(idsize,2,nmesh),idonorg(idsize,2,nmesh)
      integer ibcg(idsize,nmesh)
      real fracg(idsize,2,nmesh)

c..local variables

      integer bcdim,ii,jj,kk,iip,jjp,kkp,id,j,k,n,is,nf,ndon,im
      integer qptr,vnutptr,nfringe,ndonor,iisptr,iieptr
      real dj0,dj1,dk0,dk1
      real w1,w2,w3,w4
      real,allocatable :: qbc(:,:),vnutbc(:)
      real,allocatable :: q(:,:,:),vnut(:,:)

      bcdim=iieptrg(nmesh)

      allocate(qbc(bcdim,nv),vnutbc(bcdim))

c...LOOP THROUGH ALL THE MESHES AND COLLECT GLOBAL QBC

      qptr = 1
      vnutptr = 1
      do im = 1,nmesh

        ndonor = ndonorg(im)
        iisptr = iisptrg(im)
        iieptr = iieptrg(im)
        jmax = jmx(im); kmax = kmx(im)
        
        allocate(q(jmax,kmax,nq),vnut(jmax,kmax))
         
c.....assign local q,vnut from global values

        do n = 1,nq
          do k=1,kmax
            do j=1,jmax
              q(j,k,n) = qg(qptr-1+jmax*kmax*(n-1) + jmax*(k-1) + j)
            enddo
          enddo
        enddo

        do k=1,kmax
          do j=1,jmax
            vnut(j,k) = vnutg(vnutptr-1+jmax*(k-1) + j)
          enddo
        enddo
       
        qptr = qptr + jmax*kmax*nq
        vnutptr = vnutptr + jmax*kmax
   
        do id=1,ndonor

          ii=idonorg(id,1,im)
          jj=idonorg(id,2,im)

          iip=min(ii+1,jmax)
          jjp=min(jj+1,kmax)

          dj0= 1.-fracg(id,1,im)
          dk0= 1.-fracg(id,2,im)
          dj1=fracg(id,1,im)
          dk1=fracg(id,2,im)
          
          w1   = (dj0*dk0)
          w2   = (dj1*dk0)
          w3   = (dj0*dk1)
          w4   = (dj1*dk1)
          
c.....collect in global pointer qbc from pointer iisptr->iieptr

          do n=1,nmv
             qbc(iisptr-1+id,n)= 
     &                   w1*q(ii ,jj ,n)*q(ii ,jj ,nq) 
     &                 + w2*q(iip,jj ,n)*q(iip,jj ,nq)
     &                 + w3*q(ii ,jjp,n)*q(ii ,jjp,nq) 
     &                 + w4*q(iip,jjp,n)*q(iip,jjp,nq)
          enddo

          do n=nmv+1,nv
             qbc(iisptr-1+id,n)= 
     &                   w1*q(ii ,jj ,n)
     &                 + w2*q(iip,jj ,n)
     &                 + w3*q(ii ,jjp,n)
     &                 + w4*q(iip,jjp,n)
          enddo
          
          vnutbc(iisptr-1+id)= 
     &                 w1*vnut(ii ,jj )
     &               + w2*vnut(iip,jj )
     &               + w3*vnut(ii ,jjp)
     &               + w4*vnut(iip,jjp)
        enddo

        deallocate(q,vnut)
      enddo

c...LOOP THROUGH ALL MESHES AND UPDATE VALUES FROM QBC ARRAY

      qptr = 1
      vnutptr = 1
      do im = 1,nmesh

        nfringe = nfringeg(im)
        jmax = jmx(im); kmax = kmx(im)

        allocate(q(jmax,kmax,nq),vnut(jmax,kmax))

c.....re-assign local q,vnut from global values
        do n=1,nq
          do k=1,kmax
            do j=1,jmax
              q(j,k,n) = qg(qptr-1+jmax*kmax*(n-1) + jmax*(k-1) + j)
            enddo
          enddo
        enddo

        do k=1,kmax
          do j=1,jmax
            vnut(j,k) = vnutg(vnutptr-1+jmax*(k-1) + j)
          enddo
        enddo
       
c.....over write fringe points solution w/ donor global qbc solution

        do id=1,nfringe
          j = ibcg(id,im)
          ii = imeshg(id,1,im)
          jj = imeshg(id,2,im)

          do n = 1,nmv
             q(ii,jj,n) = qbc(j,n)/q(ii,jj,nq)
          enddo

          do n = nmv+1,nv
             q(ii,jj,n) = qbc(j,n)
          enddo

          vnut(ii,jj) = vnutbc(j)
        enddo
       
c.....reassign qbc to global q (containing all Ng meshes)
        do n=1,nv
          do k=1,kmax
            do j=1,jmax
              qg(qptr-1+jmax*kmax*(n-1) + jmax*(k-1) + j) = q(j,k,n)
            enddo
          enddo
        enddo

        do k=1,kmax
          do j=1,jmax
            vnutg(vnutptr-1+jmax*(k-1) + j) = vnut(j,k)
          end do
        end do

        qptr = qptr + jmax*kmax*nq
        vnutptr = vnutptr + jmax*kmax

        deallocate(q,vnut)
      enddo 

      return
      end

c************************************************************************
      subroutine source_bodyforce(q,s)
c        
c  add source term to the rhs
c     
c************************************************************************
     
      use params_global
      implicit none
      
      real q(jmax,kmax,nq), s(jmax,kmax,nv)
      integer k,j
      
      j = 90
      k = 42 
      s(j,k,2) = s(j,k,2) + 1.2078e-5/q(j,k,nq)

c..      do j=2,jmax-1     
c..      do k=2,kmax-1
c..           s(j,k,2) = s(j,k,2) + 2.5273e-9           
c..      enddo
c..      enddo

      return
      end

     
c***********************************************************************

