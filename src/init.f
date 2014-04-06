      subroutine init(kt0, t0)
c     -----------------------------------------------------------------
c     physical solution initialization
c     -----------------------------------------------------------------
      include 'param2d.h'
c     -----------------------------------------------------------------
c
         kt0=0
         t0=0.0
c
         do 200 is=1,ns
            ua(1,is)    = roin
            ua(2,is)    = roin*uxin
            ua(3,is)    = roin*uyin
            ua(4,is)    = ein
            ua(5,is)    = xkin
            ua(6,is)    = xein                              
            reyturb(is) = 0.0
            if(iturb.ne.0) reyturb(is) = xkin*xkin*0.09/xein
200      continue   
c
c
c     reinitialization strategy
c
         if(ncont.eq.1) then
         open(20,file="INIT_NS",status="unknown")
         rewind(20)
         do is=1,ns
c
         read(20,*) (ua(i,is),i=1,4)
c
         enddo
          read(20,*) kt0,t0
         close(20)
         endif
         if(ncont1.eq.1) then
         open(21,file="INIT_KE",status="unknown")
         rewind(21)
         do is=1,ns
         read(21,*) ua(5,is),ua(6,is),reytot,reyturb(is)
         reylam(is)=abs(reytot-reyturb(is))
         enddo
         close(21)
         endif
c
      if (ivis.eq.1.and.ilaw.eq.0) then
         do 9001 i=1,nslog3
            is=log3(i)
            ro=ua(1,is)
            ru=ua(2,is)
            rv=ua(3,is)
            ua(2,is) = 0.0
            ua(3,is) = 0.0
            if(iecc.eq.2) then
              ua(4,is) = ua(1,is)*tbrd2
                          else
              ua(4,is) = ua(4,is)-0.5*(ru**2+rv**2)/ro
              if(ua(4,is).lt.0.) stop 'init ua(4,is)<0'
            endif
            ua(5,is) = 0.0
9001     continue
      endif 
c
      return
      end
