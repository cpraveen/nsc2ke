      subroutine caldtl(dtmin,dt)
      include 'param2d.h'
c
      cste              = 110.*tinf/tinfd
cc
cc       sutherland   viscosity
cc
         if (ivis .ne. 0) then
         do is=1,ns
         temp = pres(is)/(gam1*ua(1,is))
         xnum=tinf+cste
         xden=temp+cste
         reylam(is)=reynolds*(temp/tinf)**1.5*(xnum/xden)
         enddo
         endif
c
            if(iloc.eq.2) then
            xcoef=2.*gam/pr
            else
            xcoef=0.0
            endif

            dtmin              = 1.0e+12
            do is=1,ns                                       
               reytot=reylam(is)+reyturb(is)
               dsrey= xcoef*reytot/ua(1,is)
               unorm=sqrt(ua(2,is)*ua(2,is) +
     &                    ua(3,is)*ua(3,is))
               sigmax=gam*pres(is)/ua(1,is)
               sigmax=max(1.e-2,sqrt(sigmax)+unorm)
               dtl(is)       = cfl*dthaut(is)*dthaut(is)/
     &                        (dthaut(is)*sigmax + ivis*dsrey)
               dtmin           = min(dtmin, dtl(is))
            enddo
c
            if (iloc .eq. 0) then
               do is=1,ns
                  dtl(is)    = dtmin
               enddo
            endif
c
            dt               = min(dtmin, tmax-t)
            t                = t + dt
c
           return
           end
