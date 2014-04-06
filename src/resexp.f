       subroutine resexp(ialpha)
c
c explicit resolution
c
       include 'param2d.h'
c
         do 200 is=1,ns
         do 200 ivar=1,nvar
            ce(ivar,is) = 0.0
200      continue
c
         do is=1,ns
            ua(2,is)         = ua(2,is)/ua(1,is)
            ua(3,is)         = ua(3,is)/ua(1,is)
            ua(5,is)         = ua(5,is)/ua(1,is)
            ua(6,is)         = ua(6,is)/ua(1,is)
         enddo

         call viscdt
c
         if(iflux.eq.1) call fluroe
         if(iflux.eq.2) call fluosh
         if(iflux.eq.3) call flucin
c
         if(iturb.ne.0) call keps2d
c
         do 777 is=1,ns
            ua(2,is)         = ua(2,is)*ua(1,is)
            ua(3,is)         = ua(3,is)*ua(1,is)
            ua(5,is)         = ua(5,is)*ua(1,is)
            ua(6,is)         = ua(6,is)*ua(1,is)
777      continue      
c
         if(abs(froud).gt.1.e-3) call gravity
         if(ilaw.ge.1) call vitfrot
         if(iaxi.eq.1) call sa
c
c        boundary conditions treatment
c
         if(iflux.ne.0) call cdl
c
c        updating the physical solution 
c
         do 52 is =1,ns
            usais            = dtl(is)/airsa(is)
         do 52 ivar=1,nvar
            ua(ivar,is)  = un(ivar,is) + alpha(ialpha)*usais*ce(ivar,is)
52       continue
c
         do is=1,ns
            ua(5,is)  = abs(ua(5,is))
            ua(6,is)  = abs(ua(6,is))
         enddo
c
        return
       end
