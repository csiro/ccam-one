      subroutine smooth(a,il,jl,value,b,niter,lsm_gbl)

! now assumes actual spval < value

c     routine smooths interior of an array which has undefined points
      real a(il*jl)         ! input and output array
      real value            ! values below this need smoothing
      real b(il*jl)
      real lsm_gbl(il*jl)
      integer niter

      idia=1090
      jdia=540
      iqdia=idia+(jdia-1)*il
      write(6,*)iqdia,idia,jdia,il

      write(6,*)"smooth il,jl,value,niter,iqdia=",
     &                  il,jl,value,niter,iqdia

c******************************************************************************
      do iter=1,niter
c******************************************************************************
        b=a

        im=mod(iter,2)
!       im=0

!       write(6,*)"########################iter,im,a=",iter,im,a(iqdia)

        npts=0

        do j=3,jl-2
        do i=3,il-2
          iq=i+(j-1)*il

          if(lsm_gbl(iq).lt.value)then
            a(iq)=(3.*b(iq-2)+3.*b(iq-1)+3.*b(iq)+3.*b(iq+1)+3.*b(iq+2)
     &            +b(iq-2*il)+b(iq-il)+b(iq+il)+b(iq+2*il))/19.
          endif ! filval

        enddo ! i
        enddo ! j

!       write(6,*)"a(iqdia),b(iqdia),npts="
!    &            ,a(iqdia),b(iqdia),npts

c******************************************************************************
      enddo ! iter
c******************************************************************************

      return
      end
