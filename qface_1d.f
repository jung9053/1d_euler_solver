      !----------------------------------------------------------------
      !================================================================
      !----------------------------------------------------------------

      subroutine qface_1d(qf1,qf2,q)
      use param
      use grid_dimen
      use sol_param
      implicit none
      real(8) :: q(nvar,nx-1)
      real(8) :: qf1(nvar,nx),qf2(nvar,nx)
      !local
      integer :: i

            
      select case(iorder)       
     

         case(1)
           do i=1,nx-1
              qf1(1:nvar,i+1) = q(1:nvar,i)
              qf2(1:nvar,i)   = q(1:nvar,i)
           end do

         case(2)
            call muscl_recon(q,qf1,qf2)
           
         case(3) 
            call weno5_recon(q,qf1,qf2)

         case(4)
            !write(*,*) "CRWENO is not implemented yet!"
            !stop
            call crweno5_recon(q,qf1,qf2)

         case default
         write(*,*) 'There is no available scheme!'
         stop
      
      end select


      !-> update boundary conditions
      call ubc_1d(qf1(1,1),qf2(1,nx),q(1,1),q(1,nx-1))

      return
      end subroutine qface_1d

      !----------------------------------------------------------------
      !================================================================
      !----------------------------------------------------------------

      subroutine ubc_1d(qf_min,qf_max,qq0,qq1)
      use param
      implicit none
      real(8) :: qf_min(nvar),qf_max(nvar)
      real(8) :: qq0(nvar),qq1(nvar)
      !local
      integer :: i

      do i=1,nvar
         qf_min(i) = qq0(i)
         qf_max(i) = qq1(i)
      end do

      return
      end subroutine ubc_1d
