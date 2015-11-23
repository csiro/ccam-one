      program test

      parameter (il=10,jl=30)

      real a(il,jl)
      real b(il,jl)

      read(*,*) niter
      write(*,*) niter

!     write(6,*)"init"
!     do j = jl,1,-1
!       write(6,'(30f3.0)') (a(i,j),i=1,il)
!     enddo

      do j=1,jl
       do i=1,il
         if(i+j.lt.20) then
           a(i,j)=10
         else
           a(i,j)=0
         endif
       enddo
      enddo

      write(6,*)"initial"
      do j = jl-1,2,-1
        write(6,'(30f5.1)') (a(i,j),i=1,il)
      enddo

      call smooth(a,il,jl,30.,b,niter)

      write(6,*)"final"
      do j = jl,1,-1
        write(6,'(30f5.1)') (a(i,j),i=1,il)
      enddo

      stop
      end
      subroutine smooth(a,il,jl,value,b,niter)

! now assumes actual spval < value

c     routine fills in interior of an array which has undefined points
      real a(il*jl)         ! input and output array
      real value            ! values below this need smoothing
      real b(il*jl)
!     real lsm_gbl(il*jl)
      integer niter
      logical filval

      idia=1090
      jdia=539
      iqdia=idia+(jdia-1)*il
      write(6,*)iqdia,idia,jdia,il

      write(6,*)"smooth il,jl,value,niter,iqdia=",
     &                  il,jl,value,niter,iqdia

      b=a

c******************************************************************************
      do iter = 1,niter
c******************************************************************************

        im=mod(iter,2)
        im=0

!       write(6,*)"########################iter,im,a=",iter,im,a(iqdia)

        npts=0

        do j=1,jl
        do i=1,il
          iq=i+(j-1)*il

          iw=max(1,min(il,i-(1-im)))
          ie=max(1,min(il,i+(1+im)))
          jn=max(1,min(jl,j+(1+im)))
          js=max(1,min(jl,j-(1+im)))

          iw=iw+(j-1)*il
          ie=ie+(j-1)*il
          jn=i+(jn-1)*il
          js=i+(js-1)*il

!         filval=(lsm_gbl(iq).lt.value)
!         filval=filval.and.(lsm_gbl(iw).lt.value)
!         filval=filval.and.(lsm_gbl(ni).lt.value)
!         filval=filval.and.(lsm_gbl(ie).lt.value)
!         filval=filval.and.(lsm_gbl(is).lt.value)
          filval=.true.

!         if(i.eq.idia)then
!           write(6,*)"i,j,iq,filval,lsm_gbl,value="
!    &                ,i,j,iq,filval,lsm_gbl(iq),value
!         endif
!         if(iq.eq.iqdia)then
!           write(6,*)"iq,filval,lsm_gbl,value="
!    &                ,iq,filval,lsm_gbl(iq),value
!         endif

          if(filval)then

            npts=npts+1

!           b(iq)=(4.*a(iq)+a(iw)+a(ni)+a(ie)+a(is))/8.
            b(iq)=(a(iq)+a(iw)+a(ni)+a(ie)+a(is))/5.

!           if(iq.eq.iqdia)then
!             write(6,*)iq,iw,ni,ie,is
!             write(6,*)a(iq),a(iw),a(ni),a(ie),a(is)
!              write(6,*)"iq,b(iq),npts="
!     &                ,iq,b(iq),npts
!           endif

          endif ! filval

        enddo ! i
        enddo ! j

        write(6,*)"a(iqdia),b(iqdia),npts="
     &            ,a(iqdia),b(iqdia),npts

        a=b

c******************************************************************************
      enddo ! iter
c******************************************************************************

      return
      end
