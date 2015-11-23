      subroutine fill(a,il,jl,value,b,maxit)

! now assumes actual spval < value

c     routine fills in interior of an array which has undefined points
      real a(il,jl)         ! input and output array
      real value            ! array value denoting undefined
      logical mytest
      real b(il,jl)
      integer maxit
      dimension in(8), jn(8)   ! specifies neighbours
      dimension in6(6), jn6(6)   ! specifies neighbours
      data in/-1,-1,-1,0,1,1, 1, 0/
      data jn/-1, 0, 1,1,1,0,-1,-1/
      data in6/-3,-2,-1, 1, 2, 3/
      data jn6/ 0, 0, 0, 0, 0, 0/

      write(6,*)"fill il,jl,value,maxit=",il,jl,value,maxit

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

      write(6,*)"niter,nrem=",niter,nrem
      mytest=(nrem.gt.0)
      if(maxit.gt.1)mytest=mytest.and.(niter.le.maxit)

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
      return
      end
c******************************************************************************
