      subroutine fillzonal(a,il,jl,value,b)

! now assumes actual spval < value

c     routine fillzonal in interior of an array which has undefined points
      real a(il,jl)         ! input and output array
      real value            ! array value denoting undefined
      real b(il,jl)

      write(6,*)"fillzonal il,jl,value=",il,jl,value

      do j=1,jl

        sum=0.
        np=0

        do i=1,il
          if(a(i,j).gt.value)then
            sum=sum+a(i,j)
            np=np+1
          endif
        enddo ! i

        if(np.gt.0)then
          do i=1,il
            if(a(i,j).lt.value)a(i,j)=sum/real(np)
          enddo ! i
!         write(6,*)j,sum,np,sum/real(np)
        endif

      enddo ! j

      return
      end
