      subroutine fillj(a,il,jl,value,b)

! now assumes actual spval < value

c     routine fills in interior of an array which has undefined points
      real a(il,jl)         ! input and output array
      real value            ! array value denoting undefined
      logical mytest
      real b(il,jl)
      dimension in(8), jn(8)   ! specifies neighbours
      real bb(il,jl),cc(il,jl),dd(il,jl)  ! JLM
      data in/-1,-1,-1,0,1,1, 1, 0/
      data jn/-1, 0, 1,1,1,0,-1,-1/

      write(6,*)"fill il,jl,value=",il,jl,value

! start JLM  for further fixing of land points when extending tss over land
! make sure only done for correct "value"      
      dd=a

      mytest=.true.
      niter=0
      nrem=0
      nrems=0
      neighb=0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do while(mytest)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      nrems=nrem
      nrem=0
      niter=niter+1

      b=a

      do j=2,jl-1
      do i=2,il-1

        if(a(i,j).lt.value)then

          neighb=0
          av=0.

          do nbs=1,8
            ni=min(il,max(1,i+in(nbs)))
            nj=min(jl,max(1,j+jn(nbs)))
            if(b(ni,nj).gt.value)then
                neighb=neighb+1
                av=av+b(i+in(nbs),j+jn(nbs))
            endif
          enddo

          if(neighb.gt.0)then
            a(i,j)=av/neighb
          else
            nrem=nrem+1    ! number of remaining points
          endif

        endif

      enddo ! i
      enddo ! j

      mytest=(nrem.gt.0)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      enddo ! mytest
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     fix up any boundary points
      do 7 i=2,il-1
      if(a(i,1).lt.value)a(i,1)=a(i,2)
      if(a(i,jl).lt.value)a(i,jl)=a(i,jl-1)
7     continue
      do 8 j=2,jl-1
      if(a(1,j).lt.value)a(1,j)=a(2,j)
      if(a(il,j).lt.value)a(il,j)=a(il-1,j)
8     continue
      if(a(1,1).lt.value)a(1,1)=a(2,2)
      if(a(il,1).lt.value)a(il,1)=a(il-1,2)
      if(a(1,jl).lt.value)a(1,jl)=a(2,jl-1)
      if(a(il,jl).lt.value)a(il,jl)=a(il-1,jl-1)
           
      write(6,*)"this is section for doing x linear fits"
      bb=a
      cc=a

      do j=1,jl
         itest=0
         do i=1,il
          if(dd(i,j).lt.value)then  ! land point
           if(itest==0)then
             i1=i
             itest=1
           else
             i2=i
           endif
          else  ! sea point
            if(itest==1)then   ! reached end of land segment
              do ii=i1,i2
              bb(ii,j)=a(i1,j)+(ii-i1)*(a(i2,j)-a(i1,j))/max(1,i2-i1)
              enddo
              write(6,*)"i ",i,j,bb((i1+i2)/2,j),i1,i2,a(i2,j),a(i1,j)
              itest=0
            endif
          endif
        enddo ! i
      enddo ! j
c this is section for doing y linear fits
        do i=1,il
         itest=0
         do j=1,jl
         if(dd(i,j).lt.value)then  ! land point
           if(itest==0)then
             j1=j
             itest=1
           else
             j2=j
           endif
          else  ! sea point
            if(itest==1)then   ! reached end of land segment
              do jj=j1,j2
              cc(i,jj)=a(i,j1)+(jj-j1)*(a(i,j2)-a(i,j1))/max(1,j2-j1)
              enddo
              write(6,*)"j ",i,j,cc(i,(j1+j2)/2),j1,j2,a(i,j2),a(i,j1)
              itest=0
            endif
          endif
      enddo
      enddo
      a=.5*(bb+cc)     ! end JLM
      
      return
      end
c******************************************************************************
