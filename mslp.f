      subroutine mslp ( ps, pmsl, zs, t, sig, ndim, npts, kl )
c------------------------------------------------------------------------------
c
c routine to calculate the mean sea level pressure
c
c   input:
c     ps = sfc.pres.
c     zs = sfc geop in m2/s2
c     t  = temps ( need level 2 temps )
c     sig= sigma values
c     ndim = 1st dim.
c     npts = number of points to do
c     kl = 2nd dim.
c   output:
c     ps = mslp (same units as ps)
c
c------------------------------------------------------------------------------

      parameter ( c=9.806/6.5e-3, conr=c/287., rconr=287./c )

      dimension t(ndim,kl),ps(ndim),pmsl(ndim),zs(ndim)
      dimension sig(kl)

      data ktemp/2/

      con=sig(ktemp)**rconr/c

!     write(6,*)c,conr,rconr,con,ktemp

      do i=1,npts
!       write(6,*)i,ps(i),zs(i),t(i,ktemp)
        pmsl(i) = ps(i)*(1.+con*zs(i)/t(i,ktemp))**conr
      end do ! i=1,npts

      return
      end
