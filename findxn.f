      subroutine findxn ( ar, nval, spp, amax, k1, amin, k2 )
c
c determines maximum (amax) of ar located at k1
c        and minimum (amin) of ar located at k2
c
c ar dimensioned nval
c
c ignores any values less than spp
c
      dimension ar(nval)
c
c preset amax and amin
c
      amax = ar(1)
      amin = ar(1)
      k1   = 1
      k2   = 1
c
c find maximum value and its location
c
      xm=0.
      num=0
      do 100 i=1,nval
        if ( ar(i).lt.spp ) go to 100
        xm=xm+ar(i)
        num=num+1
        if ( ar(i).lt.amax ) go to 100
        k1 = i
        amax = ar(i)
100   continue
      xm=xm/float(num)

c find minimum value and its location

      rms=0.
      do 200 i=1,nval
        if ( ar(i).lt.spp ) go to 200
        rms=rms+(ar(i)-xm)**2
        if ( ar(i).gt.amin ) go to 200
        k2 = i
        amin = ar(i)
 200  continue
      rms=sqrt(rms/float(num))

      write(6,300) xm, rms, amax, k1, amin, k2
 300  format("m=",f12.4," rms=",f12.4," x=",f12.6," i=",i6
     &         ," n=",f12.4," i=",i6)

      return
      end
