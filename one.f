      program one

      use ccinterp
      use cll_m
      use latlong_m
      use sigdata_m

      implicit none

      character*80 inf,zsfil,ofile

      common/mapproj/du,tanl,rnml,stl1,stl2

      include 'netcdf.inc'

      real pi, g, ds, du, tanl
      real rnml, stl1, stl2
      parameter ( pi=3.1415926536 )
      parameter ( g=9.80616 )
      
      include 'cdfind.h'

      integer  il,jl,kl,ifull,ix,iy
      integer inzs, io_out, igd, jgd
      integer id, jd, ntimes, in, ids, ide
      integer jds, jde, iccs, icce, jccs, jcce
      integer, dimension(2) :: ccdim
      
      logical netcdf_zsfil, debug, idv
      integer ierr, ncid, varid
      integer sarch, nsmooth, ifill
      integer ilx, jlx, i, j, iq, ijd
      integer ier, ndims, nvars, ngatts
      integer irecd, lonid, latid
      integer ilatx, ilonx, ilatn, ilonn
      integer narch, ivtim, iarch
      integer ntau, idvar
      integer i1, i2, iyr, imn, idy, ihr, imi
      integer ipval, nlpnts, nopnts
      integer ijgd, kx, ieri, iero, kn
      integer mtimer
      integer icmonth_to_imn
      integer, dimension(3) :: spos, npos

      real ther
      common/lconther/ther
      logical sdiag

      real dst
      real rlong0, rlat0, schmidt
      real rlatx, rlatn, rlonx, rlonn
      real elon, elat, rdays
      real spval, xa, an
      real, dimension(2) :: lonlat
      real, dimension(:,:,:), allocatable :: rlld,xyz,axyz,bxyz      
      real, dimension(:,:), allocatable :: grid
      real, dimension(:), allocatable :: ax,ay,az,bx,by,bz,x,y,z
      real, dimension(:), allocatable :: datan
      real, dimension(:), allocatable :: lsm_gbl
      real, dimension(:), allocatable :: glon, glat

      real, dimension(:), allocatable :: var_m
      real, dimension(:), allocatable :: lsmg_m
      real, dimension(:), allocatable :: zsg

      common/datatype/in_type
      character*1 in_type

      logical olsm_gbl,lnleap
      logical calout
      character*60 timorg
      character*60 namein,nameout,varun,calendar
      character*60 cu
      character*3 cmonth
      character*10 header

      namelist/gnml/inf,ofile,ds,du,tanl,rnml,stl1,stl2,inzs,zsfil
     &             ,io_out,igd,jgd,id,jd,sarch,ntimes
     &             ,namein,nameout,varun
     &     ,in,calout,sdiag,nsmooth,slon,slat
     &     ,ids,ide,jds,jde
     &     ,iccs,icce,jccs,jcce,ifill

      data igd/1/,jgd/1/,id/1/,jd/1/
      data io_out/3/
      data ifill/3/
      data nsmooth/20/
      data namein/'tos'/,varun/'K'/
      data nameout/'tos'/
      data inzs/10/, in/50/
      data sarch/1/,ntimes/999999/
      data debug/.false./,calout/.true./
      data sdiag/.false./
      data ids/1085/,ide/1100/,jds/520/,jde/590/
!     data iccs/100/,icce/129/,jccs/760/,jcce/700/ ! west asia
      data iccs/1/,icce/20/,jccs/1090/,jcce/1040/ ! southern GL

      save
   
      slon=0.
      slat=90.
      lnleap=.false.
      olsm_gbl=.false.

!####################### read namelists ############################
      write(6,*)'read namelist'
      read (5, gnml)
      write(6,nml=gnml)
! read and write namelist input for vidar
!     open  ( 98, file='vidar.nml',status='unknown' )
!     read  ( 98, nml=vi )
!     write ( unit=6, nml=vi)
!####################### read namelist ############################

!####################### read topography data ############################
      write(6,*)'open ',inzs,' zsfil=',zsfil

        write(6,*)"set up cc geometry"

        ierr = nf_open(zsfil,nf_nowrite,ncid)
        netcdf_zsfil = ( ierr==0 )
        if ( netcdf_zsfil ) then
          ierr=nf_get_att_real(ncid,nf_global,'lon0',rlong0)
          ierr=nf_get_att_real(ncid,nf_global,'lat0',rlat0)
          ierr=nf_get_att_real(ncid,nf_global,'schmidt',schmidt)
          ierr=nf_inq_dimid(ncid,'longitude',varid)
          ierr=nf_inq_dimlen(ncid,varid,ilx)
          ierr=nf_inq_dimid(ncid,'latitude',varid)
          ierr=nf_inq_dimlen(ncid,varid,jlx)
        else
          open(unit=inzs,file=zsfil,status='old',form='formatted')
          read(inzs,*)ilx,jlx,rlong0,rlat0,schmidt,ds,header
        end if
        du=rlong0
        tanl=rlat0
        rnml=schmidt
        write(6,*)"gbl mapproj=",ilx,jlx,rlong0,rlat0,schmidt,ds
!        if(ilx.ne.il.or.jlx.ne.jl)
!     &     stop 'wrong topo file supplied (il,jl) for one'

        il=ilx
        jl=6*il
        ifull=il*jl

        call latlongalloc(il)
        call cllalloc(il)
        call sigdataalloc(il)

        allocate(var_m(ifull),lsmg_m(ifull))
        allocate(zsg(ifull))
        allocate(x(ifull),y(ifull),z(ifull))
        allocate(ax(ifull),ay(ifull),az(ifull))
        allocate(bx(ifull),by(ifull),bz(ifull))
        
        ! set-up CC grid
        ccdim(1)=il
        ccdim(2)=jl
        lonlat(1)=rlong0
        lonlat(2)=rlat0
        allocate(rlld(ccdim(1),ccdim(2),2),grid(ccdim(1),ccdim(2)))
        allocate(xyz(ccdim(1),ccdim(2),3),axyz(ccdim(1),ccdim(2),3))
        allocate(bxyz(ccdim(1),ccdim(2),3))        
        Call getcc(rlld,grid,xyz,axyz,bxyz,ccdim,lonlat,schmidt,dst)
        write(6,*)"dst=",dst
        do i=1,ccdim(1)
          do j=1,ccdim(2)
            iq=i+(j-1)*ccdim(1)
            rlong(iq)=rlld(i,j,1)
            rlat(iq)=rlld(i,j,2)
            x(iq)=xyz(i,j,1)
            y(iq)=xyz(i,j,2)
            z(iq)=xyz(i,j,3)
            ax(iq)=axyz(i,j,1)
            ay(iq)=axyz(i,j,2)
            az(iq)=axyz(i,j,3)
            bx(iq)=bxyz(i,j,1)
            by(iq)=bxyz(i,j,2)
            bz(iq)=bxyz(i,j,3)
          end do
        end do
        deallocate(xyz,axyz,bxyz)
        deallocate (rlld,grid)
        !call setxyz

        rlatx=-1.e29
        rlatn= 1.e29
        rlonx=-1.e29
        rlonn= 1.e29

        do iq=1,ifull

c       convert conformal cubic lats & longs to degrees (-90 to 90) & (0 to 360)
c       used in sintp16; N.B. original rlong is -pi to pi
          !rlat(iq)=rlat(iq)*180./pi
          !rlong(iq)=rlong(iq)*180./pi
          if(rlong(iq).lt.0.)rlong(iq)=rlong(iq)+360.
          if(rlat(iq).gt.rlatx)then
            rlatx=rlat(iq)
            ilatx=iq
          endif
          if(rlong(iq).gt.rlonx)then
            rlonx=rlong(iq)
            ilonx=iq
          endif
          if(rlat(iq).lt.rlatn)then
            rlatn=rlat(iq)
            ilatn=iq
          endif
          if(rlong(iq).lt.rlonn)then
            rlonn=rlong(iq)
            ilonn=iq
          endif

        enddo  ! iq loop

        write(6,*)"rlong,rlat(1,1)=",rlong(1),rlat(1)
        write(6,*)"rlong:x,n=",rlonx,ilonx,rlonn,ilonn
        write(6,*)"rlatg:x,n=",rlatx,ilatx,rlatn,ilatn

      write(6,*)'read model grid zsg = g*zs'
      if ( netcdf_zsfil ) then
        spos(1:3) = 1
        npos(1) = il
        npos(2) = jl
        npos(3) = 1
        ierr=nf_inq_varid(ncid,'zs',varid)
        ierr=nf_get_vara_real(ncid,varid,spos,npos,zsi_m)
      else
        read(inzs,*)zsi_m
      end if

      write(6,*)'convert g*zs to zs(m)'
      do iq=1,ifull
        zs(iq)=zsi_m(iq)/g ! convert ascii read in zs*g to zs(m)
      enddo !iq=1,ifull

      write(6,*)'read model grid land-sea mask (0=ocean, 1=land)'
      if ( netcdf_zsfil ) then
        spos(1:3) = 1
        npos(1) = il
        npos(2) = jl
        npos(3) = 1
        ierr=nf_inq_varid(ncid,'lsm',varid)
        ierr=nf_get_vara_real(ncid,varid,spos,npos,lsm_m)
        ierr=nf_close(ncid)
      else
        read(inzs,*)lsm_m
        close(inzs)
      end if

      ijd=id+il*(jd-1)
      write(6,*)"ijd=",ijd," zs(m)=",zs(ijd)," lsm_m=",lsm_m(ijd)
!####################### read topography data ############################

!####################### open input netcdf file ############################
      write(6,*)'inf='
      write(6,*)inf
      idnci = ncopn(inf,ncnowrit,ier)
      write(6,*)'idnci=',idnci
      if(ier.ne.0) then
        write(6,*)' cannot open netCDF file; error code ',ier
        stop
      end if

!####################### get attributes of input netcdf file ############################
      call ncinq(idnci,ndims,nvars,ngatts,irecd,ier)
      write(6,'("ndims,nvars,ngatts,irecd,ier")')
      write(6,'(5i6)') ndims,nvars,ngatts,irecd,ier

c Get dimensions
      write(6,*) "get dim1 idnci=",idnci
c turn OFF fatal netcdf errors
      call ncpopt(0)
      lonid = ncdid(idnci,'lon',ier)
      write(6,*)"lon idnci,lonid,ier=",idnci,lonid,ier
c turn on fatal netcdf errors
c     write(6,*)"NCVERBOS,NCFATAL=",NCVERBOS,NCFATAL
c     call ncpopt(NCVERBOS+NCFATAL)

      if ( ier.eq.0 ) then
        write(6,*)"idnci,lonid=",idnci,lonid
        ier= nf_inq_dimlen(idnci,lonid,ix)
        write(6,*)"input ix,ier=",ix,ier
        latid= ncdid(idnci,'lat',ier)
        ier= nf_inq_dimlen(idnci,latid,iy)
        write(6,*)"input iy,ier=",iy,ier
        allocate(glon(ix),glat(iy))
        ier = nf_inq_varid(idnci,'lon',idv)
! get glon from input dataset
        ier = nf_get_var_real(idnci,idv,glon)
        ier = nf_inq_varid(idnci,'lat',idv)
        ier = nf_get_var_real(idnci,idv,glat)
      else
        write(6,*)"now try longitude"
        lonid = ncdid(idnci,'longitude',ier)
        if ( ier==0 ) then
          write(6,*)"lonid=",lonid," ier=",ier
          ier= nf_inq_dimlen(idnci,lonid,ix)
          write(6,*)"input ix=",ix," ier=",ier
          latid= ncdid(idnci,'latitude',ier)
          ier= nf_inq_dimlen(idnci,latid,iy)
          write(6,*)"input iy=",iy
          allocate(glon(ix),glat(iy))
          ier = nf_inq_varid(idnci,'longitude',idv)
          ier = nf_get_var_real(idnci,idv,glon)
          write(6,*)"glon=",(glon(i),i=1,ix)
          ier = nf_inq_varid(idnci,'latitude',idv)
          ier = nf_get_var_real(idnci,idv,glat)
          write(6,*)"glat=",(glat(i),i=1,iy)
        else
          write(6,*) "now try xt_ocean"
          lonid = ncdid(idnci,'xt_ocean',ier)
          if ( ier==0 ) then
            write(6,*)"lonid=",lonid," ier=",ier
            ier= nf_inq_dimlen(idnci,lonid,ix)
            write(6,*)"input ix=",ix," ier=",ier
            latid= ncdid(idnci,'yt_ocean',ier)
            ier= nf_inq_dimlen(idnci,latid,iy)
            write(6,*)"input iy=",iy
            allocate(glon(ix),glat(iy))
            ier = nf_inq_varid(idnci,'xt_ocean',idv)
            ier = nf_get_var_real(idnci,idv,glon)
            write(6,*)"glon=",(glon(i),i=1,ix)
            ier = nf_inq_varid(idnci,'yt_ocean',idv)
            ier = nf_get_var_real(idnci,idv,glat)
            write(6,*)"glat=",(glat(i),i=1,iy)
          else
            write(6,*) "ERROR: Cannot find valid coordinates in ",
     &                  trim(inf)
            stop
          end if
        end if

      endif ! ( ier .eq. 0 ) then

! find min/max  for input data set
      slon = glon(1)
      elon = glon(ix)
      slat = glat(1)
      elat = glat(iy)
      write(6,*)"==================> slon=",slon," elon=",elon," ix=",ix
      write(6,*)"==================> slat=",slat," elat=",elat," iy=",iy

      write(6,*)"allocate arrays ix,iy=",ix,iy
      allocate(datan(3*ix*iy))
      allocate(lsm_gbl(2*ix*iy))

      ier = nf_inq_dimid(idnci,'time',idnt)
      write(6,*)"ier=",ier," idnt=",idnt

!########## get number of times in input netcdf file ###########

      ier= nf_inq_dimlen(idnci,idnt,narch)
      write(6,*)"ier=",ier," narch=",narch
      narch=min(narch,ntimes)

      ier = nf_inq_varid(idnci,'time',ivtim)
      write(6,*)"ier=",ier," ivtim=",ivtim

      ier = nf_get_att_text(idnci,ivtim,'calendar',calendar) ! PH check calendar
      write(6,*)"ier=",ier," calendar=",calendar
      
      if (trim(calendar) .eq. '365_day') then
      lnleap = .true.
      write(*,*) lnleap
      endif

      ier = nf_get_att_text(idnci,ivtim,'units',timorg) ! MJT quick fix
      write(6,*)"ier=",ier," timorg=",timorg

      if (ier.eq.0) then
        i=index(timorg,'since')
      else
        timorg='hours'
        i=0
      end if

      write(6,*)"i=",i

      if (i.ne.0) then
        i=scan(timorg,' ')-1
        cu=''
        cu(1:i)=timorg(1:i)
        write(6,*)"cu=",cu
        timorg(1:19)=timorg(i+8:i+26)
        write(6,*)"timorg=",timorg
        i1=scan(timorg,'-')-1
        write(6,*)i1
        i2=scan(timorg(i1+2:19),'-')-1
        write(6,*)i2
        read(timorg(1:i1),*) iyr
        write(6,*)iyr
        read(timorg(i1+2:i1+2+i2-1),*) imn
        write(6,*)imn
        read(timorg(i1+2+i2+1:i1+2+i2+2),*) idy
        write(6,*)idy
        i=scan(timorg,' ')-1
        if ( i>0 .and. i<len_trim(timorg) ) then
          read(timorg(i+2:i+3),*) ihr
          write(6,*)ihr
          read(timorg(i+5:i+6),*) imi 
          write(6,*)imi
        else
          i=scan(timorg,'T')-1
          write(6,*)i
          read(timorg(i+2:i+3),*) ihr
          write(6,*)ihr
          read(timorg(i+5:i+6),*) imi
          write(6,*)imi
        end if 
      else
        cu=timorg
        ier = nf_get_att_text(idnci,ivtim,'time_origin',timorg)
        write(6,*)"ier=",ier," timorg=",timorg
        if (ier.ne.0) stop "timorg"
        read(timorg,'(i2)') idy
        read(timorg,'(3x,a3)') cmonth
        write(6,*)"cmonth=",cmonth
        imn = icmonth_to_imn(cmonth)
        write(6,*)"imn=",imn
        read(timorg,'(9x,i2)') iyr
        read(timorg,'(12x,i2)') ihr
        read(timorg,'(15x,i2)') imi
      end if

!     if ( iyr .lt. 10 ) iyr = iyr+2000
!     if ( iyr .lt. 100 ) iyr = iyr+1900

      write(6,'("iyr,imn,idy,ihr,imi=",5i4)')iyr,imn,idy,ihr,imi

      do j=1,iy
       do i=1,ix
         datan(i+(j-1)*ix     )=glon(i)
         datan(i+(j-1)*ix+ix*iy)=glat(j)
       enddo ! i
      enddo ! j


! printout of glon
      do j=1,iy,iy-1
        do i=1,ix,ix-1
          write(6,*)i,j,datan(i+(j-1)*ix)
        enddo
      enddo

       call prt_pan(rlong,il,jl,2,'rlong')
       call prt_pan(rlat ,il,jl,2,'rlat')

      write(6,*)"============= sintp16 clon++++++++++++++++++++++++++++"
      write(6,*)"ix,iy,sdiag,il=",ix,iy,sdiag,il
      call sintp16(datan(1:ix*iy),ix,iy,clon,glon,glat,sdiag,il)

      call prt_pan(clon,il,jl,2,'clon')

! printout of glat
      do j=1,iy,iy-1
        do i=1,ix,ix-1
          write(6,*)i,j,datan(i+(j-1)*ix+ix*iy)
        enddo
      enddo

      write(6,*)"============= sintp16 clat++++++++++++++++++++++++++++"
      call sintp16(datan(1+ix*iy:2*ix*iy),ix,iy,clat,glon,glat,
     &             sdiag,il)

      call prt_pan(clat,il,jl,2,'clat pan2')

      write(6,'("ix,iy,narch=",4i5)')ix,iy,narch

      write(6,*)"++++++++++++++++++++++++++++++++++++++++++++++++++++++"

      write(6,*)' reading variables sarch,narch=',sarch,narch

c***********************************************************************
      do iarch=sarch,narch
       ntau=iarch-sarch+1
c***********************************************************************

       sdiag=.false.

       ier = nf_inq_varid(idnci,'time',ivtim)
       ier = nf_get_var1_real(idnci,ivtim,iarch,rdays)

       write(6,*)"############ iarch,ntau,rdays,cu=",iarch,ntau,rdays,cu
      
       !convert rdays to minutes
       select case(cu) ! MJT quick fix
        case('days')
           mtimer=rdays*24.*60.
      	case('hours')
           mtimer=rdays*60. 
 	case('minutes')
          ! no change	
          mtimer=rdays
	case DEFAULT
	  write(6,*) "cannot convert unknown time unit ",trim(cu)
	  stop
       end select

       write(6,*)"mtimer=",mtimer

       write(6,*)"===========================",namein

       ier = nf_inq_varid(idnci,namein,idvar)
       write(6,*)"ier=",ier," idvar=",idvar
       if ( ier .ne. 0 ) then
         write(6,*)"ier=",ier," idvar=",idvar
         stop
       endif

       write(6,*)"input data has ",namein," data, now read in gbl"

       call ncread_2d(idnci,iarch,idvar,ix,iy,datan(1:ix*iy))

       !call amap ( datan(1:ix*iy), ix, iy, namein//"_gbl", 0., 0. )

       spval=-1.e10
       write(6,*)"spval=",spval

       ipval=2

! set any NAN to spval
! set any values greater than 400 to spval (if also has NAN somewhere in grid)
       if (any(isnan(datan(1:ix*iy)))) then
         where (isnan(datan(1:ix*iy)))
           datan(1:ix*iy)=spval
         end where  
       endif
       
       where (datan(1:ix*iy).gt.400.)
         datan(1:ix*iy)=spval
       else where ( datan(1:ix*iy).lt.-40.)
         datan(1:ix*iy)=spval
       end where

!       if (any(datan(1:ix*iy).gt.400.)) then
!	  write(6,*) "Missing data found in ",namein
!	  where (datan(1:ix*iy).gt.400.)
!	    datan(1:ix*iy)=spval
!	  end where
!       end if

! create lsm mask from spvals in input dataset
! lsm_gbl(1) = 0 over ocean
! lsm_gbl(1) = 1 over land
! lsm_gbl(1*ix*iy) = 1 over ocean
! lsm_gbl(1*ix*iy) = 0 over land
       nlpnts=0
       nopnts=0
       do iq=1,ix*iy
         if ( datan(iq) .lt. .1*spval)then
           olsm_gbl=.true.
           lsm_gbl(iq)=1. ! land
           lsm_gbl(iq+ix*iy)=0. ! land
           nlpnts=nlpnts+1
         else
           lsm_gbl(iq)=0. ! ocean
           lsm_gbl(iq+ix*iy)=1. ! ocean
           nopnts=nopnts+1
         endif
       enddo

       write(6,*)"lsm_gbl ocn"
       write(6,'("glbl",30f5.0)')(real(i),i=ids,ide)
       do j=jde,jds,-1
          write(6,'(i4,30f5.0)')j,(lsm_gbl(i+(j-1)*ix+ix*iy),i=ids,ide)
       enddo

       write(6,*)"input datan"
       write(6,'("glbl",30f5.0)')(real(i),i=ids,ide)
       do j=jde,jds,-1
          write(6,'(i4,30f5.0)')j,(datan(i+(j-1)*ix),i=ids,ide)
       enddo

       write(6,*)"fill in missing values in datan spval=",.1*spval

       if ( ifill.eq.2 )then
         call fillj(datan(1:ix*iy),ix,iy,.1*spval,
     &            datan(1+2*ix*iy:3*ix*iy))
       elseif ( ifill.eq.3 )then
! first fill ocean values over land for two iterations
         call fill(datan(1:ix*iy),ix,iy,.1*spval,
     &            datan(1+2*ix*iy:3*ix*iy),-1)
       else
! first fill ocean values over land for two iterations
         call fill(datan(1:ix*iy),ix,iy,.1*spval,
     &            datan(1+2*ix*iy:3*ix*iy),2)
! now fill in rest of land with zonal mean ocean values
         call fillzonal(datan(1:ix*iy),ix,iy,.1*spval,
     &            datan(1+2*ix*iy:3*ix*iy))
! finally fill any emaining missing values
         call fill(datan(1:ix*iy),ix,iy,.1*spval,
     &            datan(1+2*ix*iy:3*ix*iy),-1)
       endif

! make second copiy for ocean values
       datan(1+ix*iy:2*ix*iy)=datan(1:ix*iy)

       write(6,*)"postfill datan"
       write(6,'("glbl",30f5.0)')(real(i),i=ids,ide)
       do j=jde,jds,-1
          write(6,'(i4,30f5.0)')j,(datan(i+(j-1)*ix),i=ids,ide)
       enddo

       call amap ( datan(1:ix*iy), ix, iy, namein//'_gblfill', 0., 0. )

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       write(6,*)"###################### do we have olsm_gbl=",olsm_gbl
       if ( olsm_gbl ) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       ijgd=igd+ix*(jgd-1)
       write(6,*)"igd,jgd,ijgd=",igd,jgd,ijgd
       ijd=id+il*(jd-1)
       write(6,*)"id,jd,ijd=",id,jd,ijd

       write(6,*)"prepare to interp. tss for sea and land separately"
       write(6,*)"igd,jgd,gtss=",igd,jgd,datan(ijgd)
       write(6,*)"putting only land values into datan"
       write(6,*)"putting only ocean values into datan(+ix*iy)"
!not done since tss already at sea level     write(6,*)"First: reduce tss to sea level"

       if ( nlpnts.eq.0 .and. nopnts.eq.0 ) then
         nlpnts=ix*iy
       endif
       write(6,*)"fill in missing values nlpnts,nopnts=",nlpnts,nopnts

! now interpolate land values to model grid
       call sintp16(datan(1:ix*iy),ix,iy,ovar,glon,glat,
     &                  sdiag,il)   ! land

! Now for ocean values
! smooth land points

       if(nopnts.gt.0)then

         write(6,*)"presmooth opnts"
         write(6,'("glbl",30f5.0)')(real(i),i=ids,ide)
         do j=jde,jds,-1
          write(6,'(i4,30f5.0)')j,(datan(i+(j-1)*ix+ix*iy),i=ids,ide)
         enddo

         write(6,*)"smooth lsm_gbl"
         write(6,'("glbl",30f5.0)')(real(i),i=ids,ide)
         do j=jde,jds,-1
          write(6,'(i4,30f5.0)')j,(lsm_gbl(i+(j-1)*ix+ix*iy),i=ids,ide)
         enddo

         if(nsmooth.gt.0)then
         call smooth(datan(1+ix*iy:2*ix*iy),ix,iy,.5,
     &        datan(1+2*ix*iy:3*ix*iy),nsmooth,lsm_gbl(1+ix*iy:2*ix*iy))
         write(6,*)"postsmooth opnts"
         write(6,'("glbl",30f5.0)')(real(i),i=ids,ide)
         do j=jde,jds,-1
          write(6,'(i4,30f5.0)')j,(datan(i+(j-1)*ix+ix*iy),i=ids,ide)
         enddo
         endif

         write(6,*)"igd,jgd,ogtss=",igd,jgd,datan(ijgd+ix*iy)
         write(6,*)"=========================> now interp. ocean data"

! now interpolate ocean values to model grid
         call sintp16(datan(1+ix*iy:2*ix*iy),ix,iy,var_m,glon,glat,
     &                  sdiag,il)   ! ocean

         !write(6,'("ccam",30f5.0)')(real(i),i=iccs,icce)
         !do j=jccs,jcce,-1
         ! write(6,'(i4,30f5.0)')j,(var_m(i+(j-1)*il),i=iccs,icce)
         !enddo

       endif!(nopnts.gt.0)then

       call prt_pan(ovar  ,il,jl,2,'tss')
       call prt_pan(var_m,il,jl,2,'tsso')

       write(6,*)"id,jd,ltss=",id,jd,ovar(ijd)
       write(6,*)"id,jd,otss=",id,jd,var_m(ijd)

       write(6,*)"now recombine two (land/ocean) fields"
!not done since tss already at sea level     write(6,*)"Also need to recompute tss at zs"

       do j=1,jl
           do i=1,il
            iq=i+(j-1)*il
!not done since tss already at sea level  ovar(i,j)=ovar(i,j)-zs(iq)*.0065
!not done since tss already at sea level  varl_m(i,j)=varl_m(i,j)-zs(iq)*.0065
! remember, land < .5 is an ocean point
            if ( lsm_m(iq) .lt. .5 ) then
               ovar(iq)=var_m(iq)  ! set to ocean interp pnt
            endif
           enddo ! i
       enddo ! j

       write(6,*)"final sst on model grid"
       !write(6,'("ccam",30f5.0)')(real(i),i=iccs,icce)
       !do j=jccs,jcce,-1
       !   write(6,'(i4,30f5.0)')j,(ovar(i+(j-1)*il),i=iccs,icce)
       !enddo

       write(6,*)"final lsm on model grid"
       !write(6,'("ccam",30f5.0)')(real(i),i=iccs,icce)
       !do j=jccs,jcce,-1
       !   write(6,'(i4,30f5.0)')j,(lsm_m(i+(j-1)*il),i=iccs,icce)
       !enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       else!(olsm_gbl)then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         write(6,*)"=========================> now interp. data"
         call sintp16(datan(1:ix*iy),ix,iy,ovar,glon,glat,sdiag,il)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       endif!(olsm_gbl)then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       write(6,*)"id,jd,ovar=",id,jd,ovar(ijd)

       write(6,*)" findxn model ovar"
       call findxn(ovar,ifull,-1.e29,xa,kx,an,kn)

       call prt_pan(ovar,il,jl,2,nameout)

       i=il/2
       j=jl/2
       iq=il/2+(jl/2-1)*il
       write(6,'("ovar=",f12.2)') ovar(iq)

      write(6,*)nameout,varun,ihr,idy,imn,iyr,ntau,mtimer,ofile,ds,il

      call outcdf(nameout,varun,ihr,idy,imn,iyr,ntau,mtimer,ofile,ds,il,
     &           lnleap) ! PH added noleap
!#######################################################################
      enddo ! do ntau=1,narch
!#######################################################################

      deallocate(datan,lsm_gbl)
      deallocate(glon,glat)

      write(6,*)'*********** Finished one ************************'
      
      deallocate(x,y,z,ax,ay,az,bx,by,bz)
      deallocate(var_m,lsmg_m)
      deallocate(zsg)

      write(6,*)"call latlongdealloc"
      call latlongdealloc      
      write(6,*)"call clldealloc"
      call clldealloc
      write(6,*)"call sigdatadealloc"
      call sigdatadealloc

      write(6,*)"idnci,idnco=",idnci,idnco
      ieri= nf_close(idnci)
      iero= nf_close(idnco)
      write(6,*)"ieri,iero=",ieri,iero

      stop
      end ! one
c***************************************************************************
      subroutine ncread_2d(idhist,ntau,idvar,il,jl,var)

      implicit none
!     include 'gblparm.h'
      include 'netcdf.inc'

      integer, intent(in) :: idhist, ntau, il, jl
      integer start(3),count(3)
      integer ier, itype, idvar, i, j, ij

      real var(il*jl), addoff, sf
      real dx, dn

      integer*2, dimension(:), allocatable :: ivar
      real*8, dimension(:), allocatable :: dvar

      character*30 name

      write(6,*)"=============================> ncread_2d idhist=",
     &  idhist
      write(6,*)"ntau=",ntau," idvar=",idvar
      write(6,*)"il=",il," jl=",jl

c read name
      ier = nf_inq_varname(idhist,idvar,name)
      write(6,*)"ier=",ier," name=",name

      if(ier.eq.0)then

!       if(il*jl.gt.nnx*nny)stop "ncread_2d il*jl.gt.nnx*nny"

        start(1) = 1
        start(2) = 1
        start(3) = ntau
        count(1) = il
        count(2) = jl
        count(3) = 1

        write(6,'("start=",3i4)') start
        write(6,'("count=",3i4)') count

        ier = nf_inq_vartype(idhist,idvar,itype)
        write(6,*)"itype=",itype," ier=",ier
        write(6,*)"nf_short=",nf_short," nf_float=",nf_float

        if ( itype .eq. nf_short ) then
           write(6,*)"variable is short"
           allocate(ivar(il*jl))
           call ncvgt(idhist,idvar,start,count,ivar,ier)
           write(6,*)"ivar(1)=",ivar(1)," ier=",ier
           write(6,*)"ivar(il*jl)=",ivar(il*jl)
        else if ( itype .eq. nf_float ) then
           write(6,*)"variable is float"
           call ncvgt(idhist,idvar,start,count,var,ier)
           write(6,*)"var(1)=",var(1)," ier=",ier
           write(6,*)"var(il*jl)=",var(il*jl)
        else if ( itype .eq. nf_double ) then
           write(6,*)"variable is double"
           allocate(dvar(il*jl))
           ier=nf_get_vara_double(idhist,idvar,start,count,dvar)
           var=dvar
           write(6,*)"var(1)=",var(1)," ier=",ier
           write(6,*)"var(il*jl)=",var(il*jl)
        else
           write(6,*)"variable is unknown"
           stop
        endif

c obtain scaling factors and offsets from attributes
        call ncagt(idhist,idvar,'add_offset',addoff,ier)
        if ( ier.ne.0 ) addoff=0.
        write(6,*)"ier=",ier," addoff=",addoff

        call ncagt(idhist,idvar,'scale_factor',sf,ier)
        if ( ier.ne.0 ) sf=1.
        write(6,*)"ier=",ier," scale_factor=",sf

      else!(ier.eq.0)then
c no data found
        do i=1,il*jl
         var(i)=0
        enddo
        sf=0.
        addoff=0.
      endif!(ier.eq.0)then

c unpack data
      dx=-1.e29
      dn= 1.e29
      do j=1,jl
        do i=1,il
          ij=i+(j-1)*il
          if ( itype .eq. nf_short ) then
           if(i.eq.1.and.j.eq.1)
     &      write(6,*)"ivar,sf,addoff=",ivar(ij),sf,addoff
            var(ij) = ivar(ij)*sf + addoff
          else
           if(i.eq.1.and.j.eq.1)
     &      write(6,*)"var,sf,addoff=",var(ij),sf,addoff
            var(ij) = var(ij)*sf + addoff
          endif
          dx=max(dx,var(ij))
          dn=min(dn,var(ij))
        end do
      end do

      write(6,*)"ncread_2d idvar=",idvar," ntau=",ntau
      write(6,*)"ncread_2d dx=",dx," dn=",dn

      if ( itype .eq. nf_short ) then
        deallocate(ivar)
      endif
      if ( itype .eq. nf_double ) then
        deallocate(dvar)
      endif

      return ! ncread_2d
      end
c***************************************************************************
      subroutine ncread_3d(idhist,ntau,idvar,il,jl,kl,var)

!     include 'gblparm.h'
      include 'netcdf.inc'

      integer start(4),count(4)

!     integer*2 ivar(nmax*35)
      integer*2, dimension(:), allocatable :: ivar

      real var(il*jl*kl)
      character*30 name

      write(6,*)"ncread_2d idhist=",idhist
      write(6,*)"ntau=",ntau," idvar=",idvar
      write(6,*)"il=",il," jl=",jl

      ier = nf_inq_varname(idhist,idvar,name)
      write(6,*)"ier=",ier," name=",name

      start(1) = 1
      start(2) = 1
      start(3) = 1
      start(4) = ntau

      count(1) = il
      count(2) = jl
      count(3) = kl
      count(4) = 1

      write(6,'("start=",4i4)') start
      write(6,'("count=",4i4)') count

c read data
      write(6,*)"idhist=",idhist," idvar=",idvar
      ier = nf_inq_vartype(idhist,idvar,itype)
      write(6,*)"ier=",ier," itype=",itype

      if ( itype .eq. nf_short ) then
         write(6,*)"variable is short"
         allocate(ivar(il*jl*kl))
         call ncvgt(idhist,idvar,start,count,ivar,ier)
      else if ( itype .eq. nf_float ) then
         write(6,*)"variable is float"
         call ncvgt(idhist,idvar,start,count,var,ier)
      else
         write(6,*)"variable is unknown"
         stop
      endif

      addoff=0.
      sf=1.
      if ( itype .eq. nf_short ) then
      write(6,*)"obtain scaling factors and offsets from attributes"
      call ncagt(idhist,idvar,'add_offset',addoff,ier)
      if ( ier.ne.0 ) addoff=0.
      write(6,*)"ier=",ier," addoff=",addoff

      call ncagt(idhist,idvar,'scale_factor',sf,ier)
      if ( ier.ne.0 ) sf=1.
      write(6,*)"ier=",ier," sf=",sf
      endif

c unpack data
      dx=-1.e29
      dn= 1.e29
      do k=1,kl
       do j=1,jl
        do i=1,il
          ijk=i+(j-1)*il+(k-1)*il*jl
          if(i.eq.1.and.j.eq.1.and.k.eq.1)
     &       write(6,*)"i,j,k,ijk=",i,j,k,ijk
      	  if ( itype .eq. nf_short ) then
           if(i.eq.1.and.j.eq.1.and.k.eq.1)
     &      write(6,*)"ivar,sf,addoff=",ivar(ijk),sf,addoff
            var(ijk) = ivar(ijk)*sf + addoff
          else
           if(i.eq.1.and.j.eq.1.and.k.eq.1)
     &      write(6,*)"var,sf,addoff=",var(ijk),sf,addoff
            var(ijk) = var(ijk)*sf + addoff
          endif
          if(i.eq.1.and.j.eq.1.and.k.eq.1)
     &      write(6,*)"var=",var(ijk)
          dx=max(dx,var(ijk))
          dn=min(dn,var(ijk))
        end do
       end do
      end do

      write(6,*)"ncread_3d idvar=",idvar," ntau=",ntau
      write(6,*)"ncread_3d dx=",dx," dn=",dn

      if ( itype .eq. nf_short ) then
         deallocate(ivar)
      endif

      return ! ncread_3d
      end
c***********************************************************************
      subroutine filt_nc(var,il,jl,kl)

      real var(il,jl,kl)

      write(6,*) "filt_nc"

!     do k=1,kl
!      do j=1,jl
!       do i=1,il
!         var(i,j,k) = var(i,j,k)
!       end do
!      end do
!     end do

      return ! filt_nc
      end
!***********************************************************************
      function icmonth_to_imn(cmonth)

      integer icmonth_to_imn
      character*(*) cmonth

      write(6,*)"icmonth_to_imn cmonth=",cmonth

      icmonth_to_imn=0
      if ( cmonth.eq.'jan' ) icmonth_to_imn=1
      if ( cmonth.eq.'feb' ) icmonth_to_imn=2
      if ( cmonth.eq.'mar' ) icmonth_to_imn=3
      if ( cmonth.eq.'apr' ) icmonth_to_imn=4
      if ( cmonth.eq.'may' ) icmonth_to_imn=5
      if ( cmonth.eq.'jun' ) icmonth_to_imn=6
      if ( cmonth.eq.'jul' ) icmonth_to_imn=7
      if ( cmonth.eq.'aug' ) icmonth_to_imn=8
      if ( cmonth.eq.'sep' ) icmonth_to_imn=9
      if ( cmonth.eq.'oct' ) icmonth_to_imn=10
      if ( cmonth.eq.'nov' ) icmonth_to_imn=11
      if ( cmonth.eq.'dec' ) icmonth_to_imn=12

      write(6,*)"icmonth_to_imn=",icmonth_to_imn

      return 
      end
!***********************************************************************
