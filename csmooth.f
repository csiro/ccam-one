      subroutine smooth(a,il,jl,value,b,niter,lsm_gbl)

      use indices_m

! now assumes actual spval < value

c     routine fills in interior of an array which has undefined points
      real a(il*jl)         ! input and output array
      real value            ! values below this need smoothing
      real b(il*jl)
      real lsm_gbl(il*jl)
      integer niter
      logical test

      write(6,*)"smooth il,jl,value=",il,jl,value

      id=1
      write(6,*)id,i_w(id),i_n(id),i_e(id),i_s(id)
      write(6,*)id,i_ww(id),i_nn(id),i_ee(id),i_ss(id)

c******************************************************************************
      do iter = 1,niter
c******************************************************************************

        im=mod(iter,2)

        write(6,*)"iter,im=",iter,im

        b=a

        npts=0

        do iq=1,il*jl

          test=(lsm_gbl(iq).gt.value)

          if(im.eq.1)then
            test=test.or.(lsm_gbl(i_ww(iq)).gt.value)
            test=test.or.(lsm_gbl(i_nn(iq)).gt.value)
            test=test.or.(lsm_gbl(i_ee(iq)).gt.value)
            test=test.or.(lsm_gbl(i_ss(iq)).gt.value)
          else
            test=test.or.(lsm_gbl(i_w(iq)).gt.value)
            test=test.or.(lsm_gbl(i_n(iq)).gt.value)
            test=test.or.(lsm_gbl(i_e(iq)).gt.value)
            test=test.or.(lsm_gbl(i_s(iq)).gt.value)
          endif
          write(6,*)"iq,test,lsm_gbl,im=",iq,test,lsm_gbl(iq),im

          if(.not.test)then

            npts=npts+1

            sum=4.*a(iq)/8.
            if(im.eq.1)then
              sum=sum+a(i_ww(iq))/8.
              sum=sum+a(i_nn(iq))/8.
              sum=sum+a(i_ee(iq))/8.
              sum=sum+a(i_ss(iq))/8.
            else
              sum=sum+a(i_w(iq))/8.
              sum=sum+a(i_n(iq))/8.
              sum=sum+a(i_e(iq))/8.
              sum=sum+a(i_s(iq))/8.
            endif ! im

            b(iq)=sum

            write(6,*)"np,iq,a,b=",npts,iq,a(iq),b(iq)

          endif ! test

        enddo ! iq

        write(6,*)"npts=",npts

        a=b

c******************************************************************************
      enddo ! iter
c******************************************************************************

      return
      end
