
      
      subroutine crweno5(a,b,c,d,e,r_ele,a_ele,b_ele,c_ele)
      implicit none
      !..global
      real(8) :: a,b,c,d,e,eps
      real(8) :: r_ele,a_ele,b_ele,c_ele
      !..local
      real(8) :: f1,f2,f3,b1,b2,b3,tau,w1_opt,w2_opt,w3_opt
      real(8) :: a1,a2,a3
      real(8) :: w3,w1,w2
    
      eps = 1e-7
            
      ! smoothness indicators
      b1 = (13./12.)*(a-2.0*b+c)*(a-2.0*b+c)+(1.0/4.0)*(a-4.0*b+3.0*c)
     .     *(a-4.0*b+3.0*c)
      b2 = (13./12.)*(b-2.0*c+d)*(b-2.0*c+d)+(1.0/4.0)*(b-d)*(b-d)

      b3 = (13./12.)*(c-2.0*d+e)*(c-2.0*d+e)+(1.0/4.0)*(3.0*c-4.0*d+e)
     .     *(3.0*c-4.0*d+e)

      !..optinal weight
      w1_opt = 0.2
      w2_opt = 0.5
      w3_opt = 0.3
      
      !..limiting
      a1 = w1_opt/((eps+b1)**2.0)
      a2 = w2_opt/((eps+b2)**2.0)
      a3 = w3_opt/((eps+b3)**2.0)

      !..making the weights convex
      w1 = a1 / (a1+a2+a3)
      w2 = a2 / (a1+a2+a3)
      w3 = a3 / (a1+a2+a3)

      ! candidate stencils
      f1 = (b + 5.0*c)/6.0
      f2 = (5.0*c + d)/6.0
      f3 = (c + 5.0*d)/6.0

      !..
      !a_ele = (2.0/3.0)*w1+(1.0/3.0)*w2
      a_ele = (2.0*w1+w2)/3.0
      !b_ele = (1.0/3.0)*w1+(2.0/3.0)*(w2+w3)
      b_ele = (w1+2.0*(w2+w3))/3.0
      c_ele = w3/3.0

      r_ele = w1*b/6.0 +(5.0*(w1+w2)+w3)/6.0*c + (w2+5.0*w3)/6.0*d

      return
      end subroutine crweno5

      !=================================================================
      !=================================================================

      subroutine crweno5_recon(q,qf1,qf2)
      use param
      use grid_dimen
      implicit none
      !..global
      real(8) :: q(nvar,nx-1)
      real(8) :: qf1(nvar,nx),qf2(nvar,nx)
      !..local array
      real(8),allocatable,dimension(:) :: q3
      real(8),allocatable,dimension(:) :: rmat,amat,bmat,cmat
      integer :: i,j
      real(8) :: r_ele,a_ele,b_ele,c_ele
      real(8) :: epsw,weno5

      allocate(rmat(nx-5),amat(nx-5),bmat(nx-5),cmat(nx-5))
      allocate(q3(nx-1))

      ! initializing for left biased
      do j=1,nvar
        do i=1,nx-1
          q3(i) = q(j,i)
        enddo

        do i=4,nx-4
           call crweno5(q3(i-2),q3(i-1),q3(i),q3(i+1),q3(i+2),
     .                 r_ele,a_ele,b_ele,c_ele)
           rmat(i-2) = r_ele
           amat(i-2) = a_ele
           bmat(i-2) = b_ele
           cmat(i-2) = c_ele

        enddo
        
        !.. weno at boundary
        epsw = 1.e-6 
        rmat(1)    =  weno5(q3(1),q3(2),q3(3),q3(4),q3(5),epsw)
        rmat(nx-5) =  weno5(q3(nx-5),q3(nx-4),q3(nx-3),q3(nx-2)
     .               ,q3(nx-1),epsw)

        amat(1)    = 0.0 
        bmat(1)    = 1.0
        cmat(1)    = 0.0

       
        amat(nx-5) = 0.0 
        bmat(nx-5) = 1.0
        cmat(nx-5) = 0.0

        call tri(1,nx-5,amat,bmat,cmat,rmat)
        qf1(j,4:nx-2) = rmat(1:nx-5) 
      
      enddo

      ! initializing for right biased 
      do j=1,nvar
        do i=1,nx-1
          q3(i) = q(j,i)
        enddo
       
        do i=4,nx-4
           call crweno5(q3(i+2),q3(i+1),q3(i),q3(i-1),q3(i-2),
     .                 r_ele,a_ele,b_ele,c_ele)

           rmat(i-2) = r_ele
           amat(i-2) = a_ele
           bmat(i-2) = b_ele
           cmat(i-2) = c_ele
        enddo

        ! weno at boundary
        rmat(1)    =  weno5(q3(5),q3(4),q3(3),q3(2),q3(1),epsw)
        rmat(nx-5) =  weno5(q3(nx-1),q3(nx-2),q3(nx-3),q3(nx-4)
     .               ,q3(nx-5),epsw)
        amat(1)    = 0.0 
        bmat(1)    = 1.0
        cmat(1)    = 0.0

        amat(nx-5) = 0.0 
        bmat(nx-5) = 1.0
        cmat(nx-5) = 0.0

        
        !write(*,*) rmat(:)
        call tri(1,nx-5,amat,bmat,cmat,rmat)
        qf2(j,3:nx-3) = rmat(1:nx-5) 
      enddo

      qf1(j,2)    = q(j,1)
      qf1(j,3)    = q(j,2)
      qf1(j,nx-1) = q(j,nx-2)
      qf1(j,nx)   = q(j,nx-1)

      qf2(j,1)    = q(j,1)
      qf2(j,2)    = q(j,2)
      qf2(j,nx-2) = q(j,nx-2)
      qf2(j,nx-1) = q(j,nx-1)

      deallocate(rmat,amat,bmat,cmat)
      deallocate(q3)

      return
      end subroutine crweno5_recon


