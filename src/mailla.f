      subroutine mailla
c     -----------------------------------------------------------------
      include 'param2d.h'
c     -----------------------------------------------------------------
c     local variables definition 
      integer is ,nunit
      integer nslogi, nslog2, nslog3, nslog4, nslog5, nslogt  
      real    xmin, xmax
      real    ymin, ymax 
      nunit   = 14
      open(nunit,file="MESH",status="old")
      rewind(nunit)
      read(nunit, * ) ns,nt
      print*,'number of vertices       : ', ns
      print*,'number of triangles      : ', nt
       do is=1,ns                                   
       read(nunit,*) i,coor(1,is),coor(2,is),logfr(is)
       enddo
       do jt=1,nt
       read(nunit,*) k,nu(1,jt),nu(2,jt),nu(3,jt),i
       enddo
c
      close(nunit)
c
      xmin         =  1.0e+3
      xmax         =-1.0e+3
      ymin         =  1.0e+3
      ymax         =-1.0e+3
      do 702 is=1,ns
         xmin      = min(xmin , coor(1,is))
         xmax      = max(xmax , coor(1,is))
         ymin      = min(ymin , coor(2,is))
         ymax      = max(ymax , coor(2,is))
702   continue
c
      print*,'----------------------------------------------'
      print*,'coordinates extrema   :'
      print*,'xmin = ', xmin
      print*,'xmax = ', xmax
      print*,'ymin = ', ymin
      print*,'ymax = ', ymax
c
c     vertices renumerotation strategy
c
      nslogi       = 0
      nslog2       = 0
      nslog3       = 0
      nslog4       = 0
      nslog5       = 0
      nslog6       = 0
      do 30 is=1,ns
         if (logfr(is) .eq. 0) nslogi  = nslogi+1
c
       if(iaxi.ne.0) then
       coor(2,is)=abs(coor(2,is))
       if(logfr(is).ge.4.and.coor(2,is).lt.1.e-6) logfr(is)=6
       endif
c
         if(logfr(is).eq.2.or.
     1     (logfr(is).eq.3.and.ivis.eq.0)) then
         nslog2  = nslog2+1
         log2(nslog2)=is
         endif
c
         numc(is)=0
         if (logfr(is).eq.3) then
         nslog3 = nslog3+1 
         log3(nslog3)=is
         numc(is)=nslog3
         endif
c
         if (logfr(is) .eq. 4) then
         nslog4  = nslog4+1
         log4(nslog4)=is
         endif
         if (logfr(is) .eq. 5) then
         nslog5  = nslog5+1
         log5(nslog5)=is
         endif
         if (logfr(is) .eq. 6) then
         nslog6 =nslog6+1
         log6(nslog6)=is
         endif
30    continue
      nslogt = nslogi+nslog2+nslog3+nslog4+nslog5+nslog6
      print*,'----------------------------------------------'
      print*,'internal mesh vertices   : ',nslogi
      print*,'slipping vertices        : ',nslog2
      print*,'no-slipping vertices     : ',nslog3
      print*,'upstream   vertices      : ',nslog5
      print*,'downstream vertices      : ',nslog4
      print*,'inflow profile vertices  : ',nslog6
      print*,'total for verification   : ',nslogt
c
56    format(2e12.5)
57    format(4i6)
78    format(12i6)            
c        
      return
      end
