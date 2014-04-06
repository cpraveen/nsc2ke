      subroutine result(dt)
c     -----------------------------------------------------------------
c     saving on file the result of a computation 
c     -----------------------------------------------------------------
      include 'param2d.h'
c     -----------------------------------------------------------------
      DIMENSION Ft(3),COORt(2,3),VAL(100)
c
       isave=2
      if(isave.eq.0) goto 1234
c
c     saving on file SOL_NS
c
      open(10,file="SOL_NS",status="unknown")
      rewind(10)
      do is=1,ns
      write(10,*)(ua(i,is),i=1,4)
      enddo
      write(10,*) kt,t
      close(10)
c
c     saving on file SOL_KE for k-epsilon
c
      if(iturb.ne.0) then
      open(11,file="SOL_KE",status="unknown")
      rewind(11)
      do is=1,ns
      reytot=reylam(is)+reyturb(is)
      write(11,*) ua(5,is),ua(6,is),reytot,reyturb(is)
      enddo
      close(11)
      endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
      open(31,file='WALL.DATA',status='unknown')
      rewind(31)
      xcoef=1./sqrt(reynolds)
      if(ivis.eq.0) nnn = nslog2
      if(ivis.eq.1) nnn = nslog3
      do i=1,nnn
      if(ivis.eq.0) is = log2(i)
      if(ivis.eq.1) is = log3(i)
      ro=ua(1,is)
      reytot=reylam(is)+reyturb(is)
      ux=dx(2,is)
      uy=dy(2,is)
      vx=dx(3,is)
      vy=dy(3,is)
      tx=(dx(4,is)/ro-pres(is)*dx(1,is)/ro**2)/gam1
      ty=(dy(4,is)/ro-pres(is)*dy(1,is)/ro**2)/gam1
      xn=vnox(is)
      yn=vnoy(is)
      cp=-2.0*(pres(is)-pin)
c      cp=pres(is)/pin
      cf=2.0*((ux*xn+uy*yn)*yn-xn*(vx*xn+vy*yn))*reytot
      ch=2.0*gam*((reylam(is)/pr)+(reyturb(is)/prt))
     1            *(tx*xn+ty*yn)
      write(31,*) coor(1,is),cp,cf*xcoef,ch*xcoef
      enddo
      close(31)
c
1234  continue
      print*,'----------------------------------------------'
      print 10,t,dt
10      format( '    physical time     = ',f10.4,
     &        /,'    minimal time step = ',e15.8)
c
      nss          = 0
      xmasu        = 0.0
      xmanu        = 1.0e+8
      pmin         = 1.0e+8
      pmax         = 0.0
      rmin         = 1.0e+8
      rmax         = 0.0
      umin         = 1.0e+8
      umax         =-1.0e+8
      vmin         = 1.0e+8
      vmax         =-1.0e+8
      xkmin        = 1.0e+8
      xkmax        =-1.0e+8
      xemin        = 1.0e+8
      xemax        =-1.0e+8
      tempmi        = 1.0e+8
      tempma        =-1.0e+8
      do 100 is=1,ns
         uz        = ua(2,is)/ua(1,is)
         vz        = ua(3,is)/ua(1,is)
         qz        = sqrt(uz*uz + vz*vz)
         pp        = pres(is)
         temp      = pp/(gam1*ua(1,is))
         cc        = sqrt(gam*pp/ua(1,is))
         rmach     = qz/cc
         xmasu     = max(xmasu, rmach)
         xmanu     = min(xmanu, rmach)
         umin      = min(umin , uz)
         umax      = max(umax , uz)
         vmin      = min(vmin , vz)
         vmax      = max(vmax , vz)
         pmin      = min(pmin , pres(is))
         pmax      = max(pmax , pres(is))
         rmin      = min(rmin , ua(1,is))
         rmax      = max(rmax , ua(1,is))
         tempmi      = min(tempmi , temp)
         tempma      = max(tempma , temp)
         xkmin      = min(xkmin , ua(5,is))
         xkmax      = max(xkmax , ua(5,is))
         xemin      = min(xemin , ua(6,is))
         xemax      = max(xemax , ua(6,is))
         if (rmach .gt. 1.) nss = nss+1
100   continue
c     
      print 998,nss
      print 1000,xmasu,xmanu,pmax,pmin,rmax,rmin,umax,umin,vmax,vmin
      if(iturb.ne.0) print 1001,xkmax,xkmin,xemax,xemin
998   format('    number of supersonic points = ',i4)
1000  format(/
     &     , '    maximum mach number            = ', e15.8,/
     &     , '    minimum mach number            = ', e15.8,/
     &     , '    maximum pressure               = ', e15.8,/
     &     , '    minimum pressure               = ', e15.8,/
     &     , '    maximum density                = ', e15.8,/
     &     , '    minimum density                = ', e15.8,/
     &     , '    maximum x-velocity             = ', e15.8,/
     &     , '    minimum x-velocity             = ', e15.8,/
     &     , '    maximum y-velocity             = ', e15.8,/
     &     , '    minimum y-velocity             = ', e15.8,/)
1001   format(/
     &     , '    maximum k                      = ', e15.8,/
     &     , '    minimum k                      = ', e15.8,/
     &     , '    maximum epsilon                = ', e15.8,/
     &     , '    minimum epsilon                = ', e15.8,/)
c
        if(isave.ne.2) return
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         open(21,file='GNU.MESH',status='unknown')
         open(22,file='GNU.PRES',status='unknown')
         open(23,file='GNU.MACH',status='unknown')
         open(24,file='GNU.TURB',status='unknown')
         open(25,file='GNU.TEMP',status='unknown')
         open(26,file='GNU.VECT',status='unknown')
         rewind(21)
         rewind(22)
         rewind(23)
         rewind(24)
         rewind(25)
         rewind(26)
cccc
c mesh for gnuplot
cccc
         do iseg =1,nseg
         nubo1=nubo(1,iseg)
         nubo2=nubo(2,iseg)
         write(21,*) coor(1,nubo1),coor(2,nubo1)
         write(21,*) coor(1,nubo2),coor(2,nubo2)
         write(21,*) 
        if((logfr(nubo1).ne.0).and.logfr(nubo2).ne.0) then
         write(22,*) coor(1,nubo1),coor(2,nubo1)
         write(22,*) coor(1,nubo2),coor(2,nubo2)
         write(22,*) 
         write(23,*) coor(1,nubo1),coor(2,nubo1)
         write(23,*) coor(1,nubo2),coor(2,nubo2)
         write(23,*) 
         write(24,*) coor(1,nubo1),coor(2,nubo1)
         write(24,*) coor(1,nubo2),coor(2,nubo2)
         write(24,*) 
         write(25,*) coor(1,nubo1),coor(2,nubo1)
         write(25,*) coor(1,nubo2),coor(2,nubo2)
         write(25,*) 
         write(26,*) coor(1,nubo1),coor(2,nubo1)
         write(26,*) coor(1,nubo2),coor(2,nubo2)
         write(26,*) 
        endif
         enddo
         
         niso=30
ccc
c iso-pressure for gnuplot
ccc
         ifile=22
         delro=(pmax-pmin)/niso
         do ii=1,niso+1
          val(ii)=pmin+(ii-1)*delro
         enddo
         do it=1,nt
          do i=1,3
          coort(1,i)=coor(1,nu(i,it))
          coort(2,i)=coor(2,nu(i,it))
          ft(i)=pres(nu(i,it))
          enddo
          call isovat(ifile,ft,coort,niso,val,xx0,yy0,xx1,yy1)
          enddo
ccc
c iso-mach for gnuplot
ccc
         ifile=23
         delro=(xmasu-xmanu)/niso
         do ii=1,niso+1
          val(ii)=xmanu+(ii-1)*delro
         enddo
         do it=1,nt
          do i=1,3
           is=nu(i,it)
          coort(1,i)=coor(1,is)
          coort(2,i)=coor(2,is)
         uz        = ua(2,is)/ua(1,is)
         vz        = ua(3,is)/ua(1,is)
         qz        = sqrt(uz*uz + vz*vz)
         pp        = pres(is)
         cc        = sqrt(gam*pp/ua(1,is))
         rmach     = qz/cc
          ft(i)=rmach
          enddo
          call isovat(ifile,ft,coort,niso,val,xx0,yy0,xx1,yy1)
          enddo
ccc
c iso-temp for gnuplot
ccc
         ifile=25
         delro=(tempma-tempmi)/niso
         do ii=1,niso+1
          val(ii)=tempmi+(ii-1)*delro
         enddo
         do it=1,nt
          do i=1,3
           is=nu(i,it)
          coort(1,i)=coor(1,is)
          coort(2,i)=coor(2,is)
          temp = pres(is)/(gam1*ua(1,is))
          ft(i)=temp
          enddo
          call isovat(ifile,ft,coort,niso,val,xx0,yy0,xx1,yy1)
          enddo
c
cc velocity vector for gnuplot
c
          ifile=26
          deltat=1.e-2
          do is=1,ns
          xx0=coor(1,is)
          yy0=coor(2,is)
          xx1=xx0+ua(2,is)*deltat/ua(1,is)
          yy1=yy0+ua(3,is)*deltat/ua(1,is)
          write(26,*) xx0,yy0
          write(26,*) xx1,yy1
          write(26,*)
          enddo
ccc
c iso-k for gnuplot
ccc
         if(iturb.ne.0) then
         ifile=24
         delro=(xkmax-xkmin)/niso
         do ii=1,niso+1
          val(ii)=xkmin+(ii-1)*delro
         enddo
         do it=1,nt
          do i=1,3
          coort(1,i)=coor(1,nu(i,it))
          coort(2,i)=coor(2,nu(i,it))
          ft(i)=ua(5,nu(i,it))/ua(1,nu(i,it))
          enddo
          call isovat(ifile,ft,coort,niso,val,xx0,yy0,xx1,yy1)
          enddo
          endif
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         close(21)
         close(22)
         close(23)
         close(24)
         close(25)
         close(26)
      
      return
      end
