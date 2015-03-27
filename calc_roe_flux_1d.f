      !----------------------------------------------------------------
      !================================================================
      !----------------------------------------------------------------

      subroutine calc_roe_flux_1d(q1,q2,ff)
      use param
      use flow_property
      implicit none
      real(8) :: q1(nvar),q2(nvar),ff(nvar)
      !..local array
      real(8),allocatable :: dw(:),lambda(:),r(:,:),fl(:),fr(:)
      !local
      integer :: i,j
      real(8) :: e0l,rhol,vl,pl,al,hl
      real(8) :: e0r,rhor,vr,pr,ar,hr
      real(8) :: rt,rho,v,h,a
      real(8) :: drho,dv,dp,da

      allocate(dw(nvar),lambda(nvar),r(nvar,nvar))
      allocate(fl(nvar),fr(nvar))

      e0l = q1(3)*gm1i + 0.5d0*q1(1)*(q1(2)*q1(2))
      rhol = q1(1)
      vl = q1(2)
      pl = q1(3)
      al = dsqrt(gam * pl / rhol)
      hl = (e0l+ pl) / rhol

      e0r = q2(3)*gm1i + 0.5d0*q2(1)*(q2(2)*q2(2))
      rhor = q2(1)
      vr = q2(2)
      pr = q2(3)
      ar = dsqrt(gam * pr / rhor)
      hr = (e0r+ pr) / rhor

      rt = dsqrt(rhor / rhol)
      rho = rt * rhol
      v = (vl + rt * vr) / (1.d0 + rt)
      h = (hl + rt * hr) / (1.d0 + rt)
      a = dsqrt((gam - 1.d0) * (h - 0.5d0 * v * v))

      drho = rhor - rhol
      dv = vr - vl
      dp = pr - pl

      dw(1) = 0.5d0 * (dp - rho * a * dv) / (a * a)
      dw(2) = -1.d0 * (dp / (a * a) - drho)
      dw(3) = 0.5d0 * (dp + rho * a * dv) / (a * a)

      lambda(1) = dabs(v - a)
      lambda(2) = dabs(v)
      lambda(3) = dabs(v + a)

      da = max(0.d0,4.d0 * ((vr - ar) - (vl - al)))
      if(lambda(1).lt.(0.9d-1 * da)) then
              lambda(1) = lambda(1) * lambda(1) / da + 0.25d0 * da
      endif
      da = max(0.d0,4.d0 * ((vr + ar) - (vl + al)))
      if(lambda(3).lt.(0.9d-1 * da)) then
              lambda(3) = lambda(3) * lambda(3) / da + 0.25d0 * da
      endif

      r(1,1) = 1.d0
      r(1,2) = 1.d0
      r(1,3) = 1.d0

      r(2,1) = v - a
      r(2,2) = v
      r(2,3) = v + a

      r(3,1) = h - v * a
      r(3,2) = 0.5d0 * v * v
      r(3,3) = h + v * a

      call euler_flux_1d(q1,fl)
      call euler_flux_1d(q2,fr)
      ff = 0.5d0 * (fl + fr)

      do i = 1,3
      do j = 1,3
      ff(i) = ff(i) - 0.5d0 * lambda(j) * dw(j) * r(i,j)
      enddo
      enddo

      deallocate(dw,lambda,r)
      deallocate(fl,fr)

      return
      end subroutine calc_roe_flux_1d

      !----------------------------------------------------------------
      !================================================================
      !----------------------------------------------------------------

      subroutine euler_flux_1d(w,f)
      use param
      use flow_property
      implicit none
      real(8) :: w(nvar), f(nvar)
      !local
      real(8) :: rho, u, p
      real(8) :: a2

      rho = w(1)
      u = w(2)
      p = w(3)

      a2 = gam * p / rho

      f(1) = rho * u
      f(2) = rho * u * u + p
      f(3) = rho * u * (a2 / gm1 + 0.5d0 * u * u)

      return
      end subroutine euler_flux_1d


