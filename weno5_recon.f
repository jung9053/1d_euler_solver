      function weno5(a,b,c,d,e,epsw)
      implicit none
      !..global
      real(8) :: a,b,c,d,e,epsw,weno5
      !..local
      real(8) :: b1,b2,djm1,ejm1,dj,ej,djp1,ejp1
      real(8) :: dis0,dis1,dis2,q30,q31,q32,d01,d02,a1ba0,a2ba0
      real(8) :: w0,w1,w2,tau
      real(8) :: wopt0,wopt1,wopt2,wtot,w
      integer :: p,iweight,wexp
      
      p = 1
      epsw = 1.e-6

      ! iweight = -1  ----> WENO-JS without mapping (Jiang and Shu)
      ! iweight =  0  ----> WENO-JS with mapping (Henick et al)
      ! iweight =  1  ----> WENO-Z   (Borges et al)
      ! iweight =  2  ----> WENO-NW  (Yamaleev et al)
      
      ! iweight =  3  ----> no limiting
      ! iweight =  4  ----> ENO3
      iweight = 0

      if(iweight.eq.0) epsw = 1.e-12
      if(iweight.eq.1) epsw = 1.e-20
      if(iweight.eq.2) epsw = 1.e-12

      b1 = 13./12.
      b2 = 1./6.
      
      djm1 = a-2.*b+c
      ejm1 = a-4.*b+3.*c
      dj   = b-2.*c+d
      ej   = b-d
      djp1 = c-2.*d+e
      ejp1 = 3.*c-4.*d+e
      
      dis0 = (b1*djm1*djm1+0.25*ejm1*ejm1)+epsw
      dis1 = (b1*dj*dj+0.25*ej*ej)+epsw
      dis2 = (b1*djp1*djp1+0.25*ejp1*ejp1)+epsw
    
      if(iweight.eq.1) then
        tau = abs(dis2-dis0)
        p   = 1
      endif
      if(iweight.eq.2) then
        tau = 2**2*(b-2.0*c+d)*(b-2.0*c+d)
        p   = 1
      endif
      if((iweight.eq.1).or.(iweight.eq.2)) then
        wexp = 2
        dis0  = (1+(tau/dis0)**wexp)**(-1)
        dis1  = (1+(tau/dis1)**wexp)**(-1)
        dis2  = (1+(tau/dis2)**wexp)**(-1)
      endif

      d01 = dis0/dis1
      d02 = dis0/dis2
      
      wopt0 = 0.1
      wopt1 = 0.6
      wopt2 = 0.3

      a1ba0 = wopt1/wopt0*(abs(d01))**p
      a2ba0 = wopt2/wopt0*(abs(d02))**p


      w0 = 1.0/(1.0+a1ba0+a2ba0)
      w1 = a1ba0*w0
      w2 = a2ba0*w0
      wtot = w0+w1+w2
      w0 = w0/wtot
      w1 = w1/wtot 
      w2 = w2/wtot

      if(iweight.eq.3) then
        w0 = wopt0+0.0*w0
        w1 = wopt1+0.0*w1
        w2 = wopt2+0.0*w2
      endif

      if(iweight.eq.4) then
        if((dis0.lt.dis1).and.(dis0.lt.dis2)) then
          w0 = 1
          w1 = 0
          w2 = 0
        elseif((dis1.lt.dis0).and.(dis0.lt.dis2)) then
          w0 = 0
          w1 = 1
          w2 = 0
        else
          w0 = 0
          w1 = 0
          w2 = 1
        endif
      endif
      

      if(iweight.eq.0) then
        w = w0
        w = w*(wopt0+wopt0*wopt0-w*3.0*wopt0+w*w)/
     .    (wopt0*wopt0+w*(1.0-2*wopt0))
        w0 = w
        w = w1
        w = w*(wopt1+wopt1*wopt1-w*3.0*wopt1+w*w)/
     .     (wopt1*wopt1+w*(1.0-2*wopt1))
        w1 = w
        w = w2
        w = w*(wopt2+wopt2*wopt2-w*3.0*wopt2+w*w)/
     .     (wopt2*wopt2+w*(1.0-2*wopt2))
        w2 = w
        wtot = w0+w1+w2;
        w0 = w0/wtot; w1=w1/wtot; w2=w2/wtot;
      endif


      q30 = 2.*a-7.*b+11.*c
      q31 = -b+5.*c+2.*d
      q32 = 2.*c+5.*d-e
      
      weno5 = b2*(w0*q30+w1*q31+w2*q32)
      
      end function weno5


      !=================================================================
      !=================================================================


      subroutine weno5_recon(q,qf1,qf2)
      use param
      use grid_dimen
      implicit none
      !..global
      real(8) :: q(nvar,nx-1)
      real(8) :: qf1(nvar,nx),qf2(nvar,nx)
      !..local
      real(8) :: epsw,weno5
      !..local array
      real(8),allocatable,dimension(:) :: q2
      integer :: i,j
     
      allocate(q2(nx-1))

      epsw = 1e-15 

      do j=1,nvar

        do i=1,nx-1
          q2(i) = q(j,i)
        enddo
        
        do i=3,nx-3
          qf1(j,i+1) = weno5(q2(i-2),q2(i-1),q2(i),q2(i+1),q2(i+2),epsw)
          qf2(j,i) = weno5(q2(i+2),q2(i+1),q2(i),q2(i-1),q2(i-2),epsw)
        enddo

        qf1(j,2)    = q(j,1)
        qf1(j,3)    = q(j,2)
        qf1(j,nx-1) = q(j,nx-2)
        qf1(j,nx)   = q(j,nx-1)

        qf2(j,1)    = q(j,1)
        qf2(j,2)    = q(j,2)
        qf2(j,nx-2) = q(j,nx-2)
        qf2(j,nx-1) = q(j,nx-1)

      enddo
      
      deallocate(q2)

      return
      end subroutine weno5_recon


