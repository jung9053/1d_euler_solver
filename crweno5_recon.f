
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
      real(8) :: eps
      real(8) :: f1,f2,f3,b1,b2,b3,w1_opt,w2_opt,w3_opt
      real(8) :: a1,a2,a3
      real(8) :: w3,w1,w2
      real(8) :: qm2,qm1,qc,qp1,qp2,qp3

      allocate(rmat(nx-1),amat(nx-1),bmat(nx-1),cmat(nx-1))
      allocate(q3(nx-1))

      amat= 0.0
      bmat= 0.0
      cmat= 0.0
      rmat= 0.0

      ! initializing for left biased
      do j=1,nvar
        do i=1,nx-1
          q3(i) = q(j,i)
        enddo

        !.. weno at boundary
        eps = 1.e-7 
        amat(1)    = 0.0 
        bmat(1)    = 1.0
        cmat(1)    = 0.0
        rmat(1)    = (1.0/3.0)*q3(1)+(5.0/6.0)*q3(2)-(1.0/6.0)*q3(3)
        
        amat(2)    = 0.0 
        bmat(2)    = 1.0
        cmat(2)    = 0.0
        rmat(2)    = -(1.0/6.0)*q3(1)+(5.0/6.0)*q3(2)+(1.0/3.0)*q3(3)
        
        amat(nx-2)    = 0.0 
        bmat(nx-2)    = 1.0
        cmat(nx-2)    = 0.0
        rmat(nx-2)    = (1.0/3.0)*q3(nx-4)-(7.0/6.0)*q3(nx-3)+
     .               (11.0/6.0)*q3(nx-2)
        amat(nx-1) = 0.0 
        bmat(nx-1) = 1.0
        cmat(nx-1) = 0.0
        rmat(nx-1) = 0.0

        do i=3,nx-3
           qm2 = q3(i-2)
           qm1 = q3(i-1)
           qc  = q3(i)
           qp1 = q3(i+1)
           qp1 = q3(i+2)

           f1 = (qm1 + 5.0*qc) / 6.0
           f2 = (5.0*qc + qp1) / 6.0
           f3 = (qc + 5.0*qp1) / 6.0

        ! smoothness indicators
        b1 = (13./12.)*(qm2-2.0*qm1+qc)*(qm2-2.0*qm1+qc)
     .      +(1.0/4.0)*(qm2-4.0*qm1+3.0*qc)*(qm2-4.0*qm1+3.0*qc)
        b2 = (13./12.)*(qm1-2.0*qc+qp1)*(qm1-2.0*qc+qp1)
     .      +(1.0/4.0)*(qm1-qp1)*(qm1-qp1)
        b3 = (13./12.)*(qc-2.0*qp1+qp2)*(qc-2.0*qp1+qp2)
     .      +(1.0/4.0)*(3.0*qc-4.0*qp1+qp2)*(3.0*qc-4.0*qp1+qp2)

        !..optinal weight
        w1_opt = 0.2
        w2_opt = 0.5
        w3_opt = 0.3
        
        !..limiting
        a1 = w1_opt/((eps+b1)**1.0)
        a2 = w2_opt/((eps+b2)**1.0)
        a3 = w3_opt/((eps+b3)**1.0)

        !..making the weights convex
        w1 = a1 / (a1+a2+a3)
        w2 = a2 / (a1+a2+a3)
        w3 = a3 / (a1+a2+a3)

!        !..mapping
!        a1 = w1*(w1_opt+w1_opt*w1_opt - 3.0*w1_opt*w1+w1*w1)/
!     .       (w1_opt*w1_opt+w1*(1.0-2.0*w1_opt))
!        a2 = w2*(w2_opt+w2_opt*w2_opt - 3.0*w2_opt*w2+w2*w2)/
!     .       (w2_opt*w2_opt+w2*(1.0-2.0*w2_opt))
!        a3 = w3*(w3_opt+w3_opt*w3_opt - 3.0*w3_opt*w3+w3*w3)/
!     .       (w3_opt*w3_opt+w3*(1.0-2.0*w3_opt))
!
!        !..making the weights convex
!        w1 = a1 / (a1+a2+a3)
!        w2 = a2 / (a1+a2+a3)
!        w3 = a3 / (a1+a2+a3)
        
        rmat(i) = w1*f1 + w2*f2 + w3*f3
        amat(i) = (2.0/3.0)*w1 + (1.0/3.0)*w2
        bmat(i) = (1.0/3.0)*w1 + (2.0/3.0)*(w2+w3)
        cmat(i) = (1.0/3.0)*w3
        enddo
        
        call tri(1,nx-1,amat,bmat,cmat,rmat)
        qf1(j,2:nx-1) = rmat(1:nx-2) 
      enddo

      ! initializing for right biased 
      amat= 0.0
      bmat= 0.0
      cmat= 0.0
      rmat= 0.0

      do j=1,nvar
        do i=1,nx-1
          q3(i) = q(j,i)
        enddo
         
        !..weno boundary
        amat(2)    = 0.0 
        bmat(2)    = 1.0
        cmat(2)    = 0.0
        rmat(2)    = (1.0/3.0)*q3(4)-(7.0/6.0)*q3(3)+(11.0/6.0)*q3(2)
        
        amat(nx-1)    = 0.0 
        bmat(nx-1)    = 1.0
        cmat(nx-1)    = 0.0
        rmat(nx-1)    = (1.0/3.0)*q3(nx-1)+(5.0/6.0)*q3(nx-2)
     .                  -(1.0/6.0)*q3(nx-3)
        
        amat(nx-2)    = 0.0 
        bmat(nx-2)    = 1.0
        cmat(nx-2)    = 0.0
        rmat(nx-2)    = -(1.0/6.0)*q3(nx-1)+(5.0/6.0)*q3(nx-2)+
     .               (1.0/3.0)*q3(nx-3)
        amat(1) = 0.0 
        bmat(1) = 1.0
        cmat(1) = 0.0
        rmat(1) = 0.0
      
        do i=3,nx-3
           qm1 = q3(i-2)
           qc  = q3(i-1)
           qp1  = q3(i)
           qp2 = q3(i+1)
           qp3 = q3(i+2)

           f1 = (qp2 + 5.0*qp1) / 6.0
           f2 = (5.0*qp1 + qc) / 6.0
           f3 = (qp1 + 5.0*qc) / 6.0

        ! smoothness indicators
        b1 = (13.0/12.0)*(qp3-2.0*qp2+qp1)*(qp3-2.0*qp2+qp1)
     .      +(1.0/4.0)*(qp3-4.0*qp2+3.0*qp1)*(qp3-4.0*qp2+3.0*qp1)
        b2 = (13.0/12.0)*(qp2-2.0*qp1+qc)*(qp2-2.0*qp1+qc)
     .      +(1.0/4.0)*(qp2-qc)*(qp2-qc)
        b3 = (13.0/12.0)*(qp1-2.0*qc+qm1)*(qp1-2.0*qc+qm1)
     .      +(1.0/4.0)*(3.0*qp1-4.0*qc+qm1)*(3.0*qp1-4.0*qc+qm1)

        !..optinal weight
        w1_opt = 0.2
        w2_opt = 0.5
        w3_opt = 0.3
        
        !..limiting
        a1 = w1_opt/((eps+b1)**1.0)
        a2 = w2_opt/((eps+b2)**1.0)
        a3 = w3_opt/((eps+b3)**1.0)

        !..making the weights convex
        w1 = a1 / (a1+a2+a3)
        w2 = a2 / (a1+a2+a3)
        w3 = a3 / (a1+a2+a3)

!        !..mapping
!        a1 = w1*(w1_opt+w1_opt*w1_opt - 3.0*w1_opt*w1+w1*w1)/
!     .       (w1_opt*w1_opt+w1*(1.0-2.0*w1_opt))
!        a2 = w2*(w2_opt+w2_opt*w2_opt - 3.0*w2_opt*w2+w2*w2)/
!     .       (w2_opt*w2_opt+w2*(1.0-2.0*w2_opt))
!        a3 = w3*(w3_opt+w3_opt*w3_opt - 3.0*w3_opt*w3+w3*w3)/
!     .       (w3_opt*w3_opt+w3*(1.0-2.0*w3_opt))
!
!        !..making the weights convex
!        w1 = a1 / (a1+a2+a3)
!        w2 = a2 / (a1+a2+a3)
!        w3 = a3 / (a1+a2+a3)


        rmat(i) = w1*f1 + w2*f2 + w3*f3
        cmat(i) = (2.0/3.0)*w1 + (1.0/3.0)*w2
        bmat(i) = (1.0/3.0)*w1 + (2.0/3.0)*(w2+w3)
        amat(i) = (1.0/3.0)*w3
        enddo
        
        call tri(1,nx-1,amat,bmat,cmat,rmat)

        qf2(j,2:nx-1) = rmat(2:nx-1) 
      enddo

      do j=1,nvar
       qf1(j,nx)   = q(j,nx-1)
       qf2(j,1)   = q(j,1)

      enddo
 
      deallocate(rmat,amat,bmat,cmat)
      deallocate(q3)

      return
      end subroutine crweno5_recon


