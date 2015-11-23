      subroutine outcdf(nameout,varun,ihr,idy,imon,iyr,iarch,mtimer
     &                 ,ofile,ds,il)
     
      integer il,jl,ifull

      common/mapproj/du,tanl,rnml,stl1,stl2

      character rundate*10

      parameter(nihead=54)
      integer nahead(nihead)

      parameter(nrhead=14)
      real ahead(nrhead)

      include 'netcdf.inc'
      character ofile*(*)
      character nameout*(*),varun*(*)

      include 'cdfind.h'

      integer dim(4)
      integer xdim,ydim,zdim,tdim
      character timorg*20
      character grdtim*33
      character*3 month(12)
      data month/'jan','feb','mar','apr','may','jun'
     &          ,'jul','aug','sep','oct','nov','dec'/

      data rundate/"analyses"/

      jl=6*il
      ifull=il*jl

      write(6,*)"outcdf iyr,imon,idy,ihr=",iyr,imon,idy,ihr
      write(6,*)"nameout,varun=",nameout,varun
      write(6,*)"iarch,mtimer=",iarch,mtimer
      write(6,*)"ds,il=",ds,il
      write(6,*)"ofile=",ofile

      dt=0
      ndt=dt
      kdate=iyr*10000+imon*100+idy
      ktime=ihr*100
      write(6,*)"kdate,ktime=",kdate,ktime

      write(6,'("outcdf iarch,ofile=",i10," ",a80)')
     &                  iarch,ofile

c#######################################################################
c netcdf output
c#######################################################################

      if ( iarch.lt.1 ) stop "wrong iarch in outcdf"

      if ( iarch.eq.1 ) then
        write(6,*)'nccre of ',ofile
#ifdef usenc4
	ier=nf_create(ofile,NF_NETCDF4,idnco)
#else
        idnco = nccre(ofile, ncclob, ier)
#endif
        write(6,*)'idnco,ier=',idnco,ier
c Turn off the data filling
        imode = ncsfil(idnco,ncnofill,ier)
        write(6,*)'imode=',imode
c Create dimensions, lon, lat
        xdim = ncddef(idnco, 'longitude', il, ier)
        ydim = ncddef(idnco, 'latitude', jl, ier)
        tdim = ncddef(idnco, 'time',ncunlim,ier)
        write(6,*) "xdim=",xdim," ydim=",ydim
     &           ," zdim=",zdim," tdim=",tdim

c define coords.
        idlon = ncvdef(idnco,'longitude',NCFLOAT,1,xdim,ier)
        call ncaptc(idnco,idlon,'point_spacing',NCCHAR,4,'even',ier)
        call ncaptc(idnco,idlon,'units',NCCHAR,12,'degrees_east',ier)
        idlat = ncvdef(idnco,'latitude',NCFLOAT,1,ydim,ier)
        call ncaptc(idnco,idlat,'point_spacing',NCCHAR,4,'even',ier)
        call ncaptc(idnco,idlat,'units',NCCHAR,13,'degrees_north',ier)
        write(6,*)'idlon,idlat=',idlon,idlat

        write(6,*)'tdim,idnco=',tdim,idnco
        idnt = ncvdef(idnco,'time',NCLONG,1,tdim,ier)
        write(6,*)'idnt=',idnt
        call ncaptc(idnco,idnt,'point_spacing',NCCHAR,4,'even',ier)

        write(6,*)'kdate,ktime,ktau=',kdate,ktime,ktau

        icy=kdate/10000
        icm=max(1,min(12,(kdate-icy*10000)/100))
        icd=max(1,min(31,(kdate-icy*10000-icm*100)))
        ich=ktime/100
        icmi=(ktime-ich*100)
        ics=0
        write(6,*) icy,icm,icd,ich,icmi,ics
        write(timorg,'(i2.2,"-",a3,"-",i4.4,3(":",i2.2))')
     &               icd,month(icm),icy,ich,icmi,ics
        write(6,*)'timorg=',timorg
        call ncaptc(idnco,idnt,'time_origin',NCCHAR,20,timorg,ier)

        write(grdtim,'("minutes since ",i4.4,"-",i2.2,"-",i2.2," ",
     &       2(i2.2,":"),i2.2)') icy,icm,icd,ich,icmi,ics
        write(6,*)'grdtim=',grdtim
        call ncaptc(idnco,idnt,'units',NCCHAR,33,grdtim,ier)

        dim(1) = xdim
        dim(2) = ydim
        dim(3) = zdim
        dim(4) = tdim

c create the attributes of the header record of the file
        nahead=0
        nahead(1)=il         ! needed by cc2hist
        nahead(2)=jl         ! needed by cc2hist
        nahead(3)=1         ! needed by cc2hist
        nahead(12)=mtimer ! should be mtimer
        write(6,*)"nahead=",nahead

        ahead=0.
        ahead(1)=ds
        ahead(5)=du ! rlong0     ! needed by cc2hist
        ahead(6)=tanl ! rlat0      ! needed by cc2hist
        ahead(7)=rnml ! schmidt    ! needed by cc2hist
        write(6,*) "ahead=",ahead

       call ncapt(idnco,ncglobal,'int_header',nclong,nihead,nahead,ier)
       if(ier.ne.0)write(6,*)"ncapt int idnco,ier=",idnco,ier
       call ncapt(idnco,ncglobal,'real_header',ncfloat,nrhead,ahead,ier)
       if(ier.ne.0)write(6,*)"ncapt real idnco,ier=",idnco,ier
       call ncaptc(idnco,ncglobal,'date_header',ncchar,10,rundate,ier)
       if(ier.ne.0)write(6,*)"ncaptc date idnco,ier=",idnco,ier

        idv=ncvdef(idnco,'ds',ncfloat,0,1,ier)
        if(ier.ne.0)write(6,*)"ncvdef ds idnco,ier=",idnco,ier
        idv=ncvdef(idnco,'du',ncfloat,0,1,ier)
        if(ier.ne.0)write(6,*)"ncvdef du idnco,ier=",idnco,ier
        idv=ncvdef(idnco,'rnml',ncfloat,0,1,ier)
        if(ier.ne.0)write(6,*)"ncvdef rnml idnco,ier=",idnco,ier
        idv=ncvdef(idnco,'tanl',ncfloat,0,1,ier)
        if(ier.ne.0)write(6,*)"ncvdef tanl idnco,ier=",idnco,ier
        idv=ncvdef(idnco,'stl1',ncfloat,0,1,ier)
        if(ier.ne.0)write(6,*)"ncvdef stl1 idnco,ier=",idnco,ier
        idv=ncvdef(idnco,'stl2',ncfloat,0,1,ier)
        if(ier.ne.0)write(6,*)"ncvdef stl2 idnco,ier=",idnco,ier
        idv=ncvdef(idnco,'dt',ncfloat,0,1,ier)
        if(ier.ne.0)write(6,*)"ncvdef dt idnco,ier=",idnco,ier
      endif ! ( iarch=1 ) then

      write(6,*)'call openhist for iarch= ',iarch

      call openhist(idnco,iarch,dim,nameout,varun
     &             ,kdate,ktime,mtimer,il,1)

      write(6,*)"done openhist"

      call ncsnc(idnco,ier)
      if(ier.ne.0)write(6,*)"ncsnc idnco,ier=",idnco,ier

      return ! outcdf
      end
c=======================================================================
      subroutine openhist(idnc,iarch,dim,nameout,varun
     &            ,kdate,ktime,mtimer,il,kl)

      use cll_m
      use sigdata_m
      use xyzinfo_m, only : em,f

c this routine creates attributes and writes output

      integer il,jl,kl,ifull

      include 'netcdf.inc'
      include 'cdfind.h'

      character lname*50,expdesc*50
      character nameout*60,varun*60
      integer dim(4)
      integer idim2(3)

      real xpnt(il),ypnt(6*il)

      real tst(6*il*il),tsb(6*il*il)
!       *** qscrn_ave not presently written     
      real aa(6*il*il),bb(6*il*il),cc(6*il*il)
      real cfrac(6*il*il,kl)

      jl=6*il
      ifull=il*jl

!     character*3 mon(12)
!     data mon/'JAN','FEB','MAR','APR','MAY','JUN'
!    &        ,'JUL','AUG','SEP','OCT','NOV','DEC'/

      write(6,*)'openhist iarch,idnc=',iarch,idnc

c if this is the first archive, set up some global attributes

      if(iarch.eq.1) then
        write(6,*)'dim=',dim
        idim2(1)=dim(1)
        idim2(2)=dim(2)
        idim2(3)=dim(4)
        write(6,*)'idim2=',idim2

c       Create global attributes
c       Model run number
        call ncapt(idnc,ncglobal,'nrun',nclong,1,nrun,ier)
        write(6,*)"nrun=",nrun," ier=",ier

c       Experiment description
        expdesc = 'CCAM model run'
        call ncaptc(idnc,ncglobal,'expdesc',ncchar,50,expdesc,ier)
        write(6,*)"expdesc=",expdesc," ier=",ier

        lname = 'year-month-day at start of run'
        idkdate = ncvdef(idnc,'kdate',nclong,1,dim(4),ier)
        call ncaptc(idnc,idkdate,'long_name',ncchar
     &             ,lngstr(lname),lname,ier)

        lname = 'hour-minute at start of run'
        idktime = ncvdef(idnc,'ktime',nclong,1,dim(4),ier)
        call ncaptc(idnc,idktime,'long_name',ncchar
     &             ,lngstr(lname),lname,ier)

        lname = 'timer (hrs)'
        idnter = ncvdef(idnc,'timer',ncfloat,1,dim(4),ier)
        call ncaptc(idnc,idnter,'long_name',ncchar
     &             ,lngstr(lname),lname,ier)

        lname = 'mtimer (mins)'
        idmtimer = ncvdef(idnc,'mtimer',nclong,1,dim(4),ier)
        call ncaptc(idnc,idmtimer,'long_name',ncchar
     &             ,lngstr(lname),lname,ier)

        lname = 'timeg (UTC)'
        idnteg = ncvdef(idnc,'timeg',ncfloat,1,dim(4),ier)
        call ncaptc(idnc,idnteg,'long_name',ncchar
     &             ,lngstr(lname),lname,ier)

        lname = 'number of time steps from start'
        idktau = ncvdef(idnc,'ktau',nclong,1,dim(4),ier)
        call ncaptc(idnc,idktau,'long_name',ncchar
     &             ,lngstr(lname),lname,ier)

        write(6,*)'define attributes of variables'

        lname = 'Surface geopotential'
        call attrib(idnc,idim2,2,'zs',lname,'m',-5.e3,12.e3)
        lname = 'Surface land mask'
        call attrib(idnc,idim2,2,'lsm',lname,'0-1',0.,1.) ! MJT lsmask

c       For time invariant surface fields
        lname = 'clon'
        call attrib(idnc,idim2,2,'clon',lname,'none',-360.,360.)
        lname = 'clat'
        call attrib(idnc,idim2,2,'clat',lname,'none',-90.,90.)

c       For time varying surface fields
        call attrib(idnc,idim2,3,nameout,nameout,varun,-999.,999.)

        write(6,*)'finished defining attributes'
c       Leave define mode
        call ncendf(idnc,ier)
        write(6,*)'leave define mode: ier=',ier

        do i=1,il
         xpnt(i) = float(i)
        end do
        call ncvpt(idnc,idlon,1,il,xpnt,ier)
        write(6,*)"xpnt ",idnc,idlon,1,il,ier
        write(6,*)"xpnt ",xpnt

        do j=1,jl
         ypnt(j) = float(j)
        end do
        call ncvpt(idnc,idlat,1,jl,ypnt,ier)
        write(6,*)"ypnt ",idnc,idlat,1,jl,ier
        write(6,*)"ypnt ",ypnt

        idv = ncvid(idnc,'ds',ier)
        call ncvpt1(idnc,idv,1,ds,ier)
        idv = ncvid(idnc,'tanl',ier)
        call ncvpt1(idnc,idv,1,tanl,ier)
        idv = ncvid(idnc,'rnml',ier)
        call ncvpt1(idnc,idv,1,rnml,ier)
        idv = ncvid(idnc,'du',ier)
        call ncvpt1(idnc,idv,1,du,ier)
        idv = ncvid(idnc,'stl1',ier)
        call ncvpt1(idnc,idv,1,stl1,ier)
        idv = ncvid(idnc,'stl2',ier)
        call ncvpt1(idnc,idv,1,stl2,ier)
        idv = ncvid(idnc,'dt',ier)
        call ncvpt1(idnc,idv,1,dt,ier)
      endif ! iarch.eq.1
!------------------------------------------------------------------      
      ktau=0
      timer=mtimer/60 ! time in hours
      timeg=mod(mtimer/60,24) ! MJT quick fix
      write(6,*)'outcdf processing kdate,ktime,ktau,time,mtimer: ',
     .                           kdate,ktime,ktau,time,mtimer
      write(6,*)'outcdf processing timer,timeg: ',
     .                           timer,mtimeg

c     set time to number of minutes since start 
      idv = ncvid(idnc,'time',ier)
      call ncvpt1(idnc,idv,iarch,mtimer,ier)
      write(6,*)"mtimer(mins)=",mtimer

      idv = ncvid(idnc,'kdate',ier)
      call ncvpt1(idnc,idv,iarch,kdate,ier)
      idv = ncvid(idnc,'ktime',ier)
      call ncvpt1(idnc,idv,iarch,ktime,ier)
      idv = ncvid(idnc,'timer',ier)
      call ncvpt1(idnc,idv,iarch,timer,ier)
      idv = ncvid(idnc,'mtimer',ier)
      call ncvpt1(idnc,idv,iarch,mtimer,ier)
      idv = ncvid(idnc,'timeg',ier)
      call ncvpt1(idnc,idv,iarch,timeg,ier)
      idv = ncvid(idnc,'ktau',ier)
      call ncvpt1(idnc,idv,iarch,ktau,ier)

      write(6,*)'kdate,ktime,ktau=',kdate,ktime,ktau
      write(6,*)'timer,timeg=',timer,timeg

      write(6,*)'now write out variables'

      if(ktau.eq.0)then
!       write time-invariant fields      
        call histwrt3(clon,'clon',idnc,-1,il)
        call histwrt3(clat,'clat',idnc,-1,il)
      endif ! (ktau.eq.0) 

      write(6,*)"zs(m)=",(zs(is+i),i=1,5)
      call prt_pan(zs,il,jl,2,'zs(m)')

      call histwrt3(zs,'zs',idnc,-1,il)   ! always from 13/9/02
      call histwrt3(lsm_m*65.e3,'lsm',idnc,-1,il)

      call histwrt3(ovar,nameout,idnc,iarch,il)

      return ! subroutine openhist(idnc,iarch,dim,sig
      end
c=======================================================================
      subroutine attrib(cdfid,dim,ndim,name,lname,units,xmin,xmax)

      include 'netcdf.inc'

      integer*2 minv, maxv, missval   ! was integer*2
      parameter(minv = -32499, maxv = 32500, missval = -32500)
      integer cdfid, idv, dim(3)
      character name*(*), lname*(*), units*(*)
      real xmin, xmax
      integer lngstr

      if ( xmin .ge. xmax ) then
        idv = ncvdef(cdfid, name, ncfloat, ndim, dim, ier)
      else
        idv = ncvdef(cdfid, name, ncshort, ndim, dim, ier)
        call ncapt(cdfid,idv,'valid_min'    ,ncshort,1,minv,ier)
        call ncapt(cdfid,idv,'valid_max'    ,ncshort,1,maxv,ier)
        call ncapt(cdfid,idv,'missing_value',ncshort,1,missval,ier)
        scalef=(xmax-xmin)/(real(maxv)-real(minv)) ! jlm fix for precision problems
        addoff=xmin-scalef*minv
        call ncapt(cdfid,idv,'add_offset',ncfloat,1,addoff,ier)
        call ncapt(cdfid,idv,'scale_factor',ncfloat,1,scalef,ier)
      endif

      if ( ier.ne.0 ) then
        write(6,*)ier,' Error in variable declaration ', name
        stop
      end if

      call ncaptc(cdfid,idv,'long_name',ncchar,lngstr(lname),lname,ier)
      if(lngstr(units).ne.0)then
        call ncaptc(cdfid,idv,'units',ncchar,lngstr(units),units,ier)
      end if
      call ncaptc(cdfid,idv,'FORTRAN_format',ncchar,5,'G11.4',ier)
      return
      end
c=======================================================================
      function lngstr( string )
      character*(*) string
      ilen = len(string)
c     print*,'string=',string
c     print*,'ilen=',ilen
      do 100 lngstr=ilen,1,-1
        if ( string(lngstr:lngstr) .ne. ' ' ) go to 99
  100 continue
      lngstr = 0
   99 continue
      return
      end
c=======================================================================
      subroutine histwrt3(var,sname,idnc,iarch,il)
c Write 2d+t fields from the savegrid array.

      integer il,jl,ifull

      integer idvar
      integer*2 ipack(il,6*il) ! was integer*2 
      character* (*) sname
c     character*8 sname
      integer*2 minv, maxv, missval ! was integer*2 
      parameter(minv = -32499, maxv = 32500, missval = -32500)

      real var(il,6*il)

      integer, dimension(:), allocatable :: start, count

      jl=6*il
      ifull=il*jl

      write(6,*)"histwrt3 sname=",sname," iarch=",iarch," idnc=",idnc
      write(6,*)"il,jl,kl,ifull",il,jl,kl,ifull

      if ( iarch .gt. 0 ) then
        allocate(start(3),count(3))
        start(1) = 1
        start(2) = 1
        count(1) = il
        count(2) = jl
        start(3) = iarch
        count(3) = 1
      else
        allocate(start(2),count(2))
        start(1) = 1
        start(2) = 1
        count(1) = il
        count(2) = jl
      endif

      write(6,*)"start=",start
      write(6,*)"count=",count

c find variable index
      idvar = ncvid(idnc,sname,ier)
      call ncagt(idnc,idvar,'add_offset',addoff,ier)
      call ncagt(idnc,idvar,'scale_factor',scale_f,ier)

      write(6,*) 'add_offset=',addoff,' scale_factor=',scale_f

      xmin=addoff+scale_f*minv
!     xmax=xmin+scale_f*float(maxv-minv)
      xmax=xmin+scale_f*(real(maxv)-real(minv)) ! jlm fix for precision problems
      write(6,*) "xmin=",xmin," xmax=",xmax

      varn= 1.e29
      varx=-1.e29
      do j=1,jl
        do i=1,il
	  if(var(i,j).lt.varn)then
	     varn=var(i,j)
	     imn=i
	     jmn=j
	  endif
	  if(var(i,j).gt.varx)then
	     varx=var(i,j)
	     imx=i
	     jmx=j
	  endif
          if(xmin.lt.xmax)then
            pvar = max(xmin,min(xmax,var(i,j))) ! limited output variable
            ipack(i,j)=nint((pvar-addoff)/scale_f)
            ipack(i,j)=max(min(ipack(i,j),maxv),minv)
          endif
	end do
      end do

      write(6,*)"pvar=",pvar

      if(xmin.lt.xmax)then
        call ncvpt(idnc, idvar, start, count, ipack, ier)
        if(ier.ne.0)stop "in histwrt3 short ier not zero"
      else
        call ncvpt(idnc, idvar, start, count, var, ier)
        if(ier.ne.0)stop "in histwrt3 float ier not zero"
      endif

      write(6,'("histwrt3:",a7," nt=",i4," n=",f12.4," ij=",2i4
     &          ," x=",f12.4," ij=",2i4)') 
     &          sname,iarch,varn,imn,jmn,varx,imx,jmx

      deallocate(start,count)

      return
      end ! histwrt3
c=======================================================================
      subroutine histwrt4(var,sname,idnc,iarch,il,kl)
c Write 3d+t fields from the savegrid array.

      integer il,jl,kl,ifull

      integer idvar, start(4), count(4)
      integer*2 ipack(il,6*il,kl) ! was integer*2 
      character* (*) sname
c     character*8 sname
      integer*2 minv, maxv, missval ! was integer*2 
      parameter(minv = -32499, maxv = 32500, missval = -32500)

      real var(il,6*il,kl)

      jl=6*il
      ifull=il*jl

      write(6,*)"histwrt4 sname=",sname," iarch=",iarch," idnc=",idnc
      write(6,*)"il,jl,kl,ifull",il,jl,kl,ifull
      start(1) = 1
      start(2) = 1
      start(3) = 1
      start(4) = iarch
      count(1) = il
      count(2) = jl
      count(3) = kl
      count(4) = 1

c find variable index
      idvar = ncvid(idnc,sname,ier)
      call ncagt(idnc,idvar,'add_offset',addoff,ier)
      call ncagt(idnc,idvar,'scale_factor',scale_f,ier)

      xmin=addoff+scale_f*minv
!     xmax=xmin+scale_f*float(maxv-minv)
      xmax=xmin+scale_f*(real(maxv)-real(minv)) ! jlm fix for precision problems

      varn= 1.e29
      varx=-1.e29
	imx=0
	jmx=0
	kmx=0
	imn=0
	jmn=0
	kmn=0
      do k=1,kl
       do j=1,jl
        do i=1,il
          pvar = max(xmin,min(xmax,var(i,j,k))) ! limited output variable
          ipack(i,j,k)=nint((pvar-addoff)/scale_f)
          ipack(i,j,k)=max(min(ipack(i,j,k),maxv),minv)
	   if(var(i,j,k).gt.varx)then
	     varx=var(i,j,k)
	     imx=i
	     jmx=j
	     kmx=k
	   endif
	   if(var(i,j,k).lt.varn)then
	     varn=var(i,j,k)
	     imn=i
	     jmn=j
	     kmn=k
	   endif
	end do
       end do
      end do

      call ncvpt(idnc, idvar, start, count, ipack, ier)

      write(6,'("histwrt4:",a7," nt=",i4," n=",f12.4," ijk=",3i4
     &          ," x=",f12.4," ijk=",3i4)') 
     &          sname,iarch,varn,imn,jmn,kmn,varx,imx,jmx,kmx

      return
      end ! histwrt4
     
      subroutine mtimerget(mtimer,kdate1,ktime1,kdate2,ktime2) ! jlm
!     returns mtimer in minutes, corr. to (kdate2,ktime2) -  (kdate1,ktime1)    
      dimension ndoy(12)   ! days from beginning of year (1st Jan is 0)
      data ndoy/ 0,31,59,90,120,151,181,212,243,273,304,334/
      common/leap_yr/leap  ! 1 to allow leap years
 
      if(leap.ne.0)stop 'leap years not catered for in mtimerget'
!     Set up number of minutes from beginning of year
!     For GCM runs assume year is <1980 (e.g. ~321-460 for 140 year run)
      jyear1=kdate1/10000
      jmonth=(kdate1-jyear1*10000)/100
      jday=kdate1-jyear1*10000-jmonth*100
      jhour=ktime1/100
      jmin=ktime1-jhour*100
      mstart1=1440*(ndoy(jmonth)+jday-1) + 60*jhour + jmin ! mins from start of y

      jyear2=kdate2/10000
      jmonth=(kdate2-jyear2*10000)/100
      jday=kdate2-jyear2*10000-jmonth*100
      jhour=ktime2/100
      jmin=ktime2-jhour*100
      mstart2=1440*(ndoy(jmonth)+jday-1) + 60*jhour + jmin ! mins from start of y

      mtimer=mstart2-mstart1+(jyear2-jyear1)*365*24*60
      return
      end
