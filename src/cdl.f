      subroutine cdl
c     -----------------------------------------------------------------
c     computation of the convective fluxes at  boundaries
c     -----------------------------------------------------------------
      include 'param2d.h'
c     -----------------------------------------------------------------
c     local variables definition
      integer  is  ,  i
      real     cap , capd1, capd2, uq41, uq42, 
     &         u, v, c    , usc  , vp1 , vp3 , vp4, 
     &         fdc1, fdc2 , fdc3 , fdc4, 
     &         fgp1, fgp2 , fgp3 , fgp4
c
c     vertices on the body, slip boundary condition 
c
      do 700 i=1,nslog2
         is=log2(i)
         ce(2,is)  = ce(2,is) - vno(is)*vnox(is)*pres(is)
         ce(3,is)  = ce(3,is) - vno(is)*vnoy(is)*pres(is)
700   continue
c
c solid bodies
c
      if(ivis.eq.1.and.ilaw.eq.0) then
c
      do 100 i=1,nslog3
         is=log3(i)
         ce(2,is)  = 0.0
         ce(3,is)  = 0.0
         ce(5,is)  = 0.0
         ce(6,is)  = 0.0
100   continue
c
      if(iecc.eq.2) then
      do 101 i=1,nslog3
         is=log3(i)
         ce(4,is)  = ce(1,is)*tbrd2
101   continue
      endif
      endif
c
c     vertices on the upstream boundary 
c
      do 850 i=1,nslog5
         is=log5(i)
         cap       = vno(is)
         capd1     = vnox(is)
         capd2     = vnoy(is)
         u         = ua(2,is)/ua(1,is)
         v         = ua(3,is)/ua(1,is)
         c         = sqrt(gam*pres(is)/ua(1,is))
         usc       = 1.0/c
         uq41      = 0.5*(u*u + v*v)*usc + gam4*c
         uq42      = capd1*u + capd2*v
c
c     computation of a-(wi).winf
c
         vp1       = (capd1*u + capd2*v)*cap 
         vp3       = min(vp1 + c*cap, 0.0)
         vp4       = min(vp1 - c*cap, 0.0)
         vp1       = min(vp1        , 0.0)
         fdc1  = vp1*(roin + gam1*usc*usc*(-0.5*(u*u + v*v)*roin   +
     &                                      u*ruxin + v*ruyin - ein))
         fdc2  = vp1*((capd1*v - capd2*u)*roin + 
     &                    capd2*ruxin - capd1*ruyin)
         fdc3  = vp3*(0.5*(-uq42*roin + capd1*ruxin + capd2*ruyin) +
     &                    0.5*gam1*usc*(0.5*(u*u + v*v)*roin - 
     &                                  u*ruxin - v*ruyin + ein))
         fdc4  = vp4*(0.5*(uq42*roin - capd1*ruxin - capd2*ruyin)  + 
     &                    0.5*gam1*usc*(0.5*(u*u + v*v)*roin - 
     &                                  u*ruxin - v*ruyin + ein))
         fgp1  = fdc1 + usc*(fdc3 + fdc4)
         fgp2  = u*fdc1 + capd2*fdc2  + 
     &               (u*usc + capd1)*fdc3 + (u*usc - capd1)*fdc4
         fgp3  = v*fdc1 - capd1*fdc2  + 
     &               (v*usc + capd2)*fdc3 + (v*usc - capd2)*fdc4
         fgp4  = 0.5*(u*u + v*v)*fdc1 + (capd2*u - capd1*v)*fdc2   +
     &               (uq41 + uq42)*fdc3   + (uq41 - uq42)*fdc4
c
         ce(1,is)  = ce(1,is) - fgp1
         ce(2,is)  = ce(2,is) - fgp2
         ce(3,is)  = ce(3,is) - fgp3
         ce(4,is)  = ce(4,is) - fgp4
         ce(5,is)  = ce(5,is) - fgp1*ua(5,is)/ua(1,is)
         ce(6,is)  = ce(6,is) - fgp1*ua(6,is)/ua(1,is)
c
c     computation of a+(wi).wi
c
         vp1   = (capd1*u + capd2*v)*cap 
         vp3   = max(vp1 + c*cap, 0.0)
         vp4   = max(vp1 - c*cap, 0.0)
         vp1   = max(vp1        , 0.0)
         fdc1  = vp1*(ua(1,is) +
     &               gam1*usc*usc*(-0.5*(u*u + v*v)*ua(1,is)    +
     &                    u*ua(2,is) + v*ua(3,is) - ua(4,is)))
         fdc2  = vp1*((capd1*v - capd2*u)*ua(1,is) +
     &                    capd2*ua(2,is) - capd1*ua(3,is))
         fdc3  = vp3*(0.5*(-uq42*ua(1,is) +
     &                capd1*ua(2,is) + capd2*ua(3,is)) +
     &             0.5*gam1*usc*(0.5*(u*u + v*v)*ua(1,is)     -
     &                           u*ua(2,is) - v*ua(3,is)      +
     &                                  ua(4,is)))
         fdc4  = vp4*(0.5*(uq42*ua(1,is)  -
     &                   capd1*ua(2,is) - capd2*ua(3,is))      +
     &              0.5*gam1*usc*(0.5*(u*u + v*v)*ua(1,is)     -
     &                            u*ua(2,is) - v*ua(3,is)      +
     &                                ua(4,is)))
         fgp1  = fdc1 + usc*(fdc3 + fdc4)
         fgp2  = u*fdc1 + capd2*fdc2  +
     &               (u*usc + capd1)*fdc3 + (u*usc - capd1)*fdc4
         fgp3  = v*fdc1 - capd1*fdc2  +
     &               (v*usc + capd2)*fdc3 + (v*usc - capd2)*fdc4
         fgp4  = 0.5*(u*u + v*v)*fdc1 + (capd2*u - capd1*v)*fdc2 +
     &               (uq41 + uq42)*fdc3   + (uq41 - uq42)*fdc4
c
         ce(1,is)  = ce(1,is) - fgp1
         ce(2,is)  = ce(2,is) - fgp2
         ce(3,is)  = ce(3,is) - fgp3
         ce(4,is)  = ce(4,is) - fgp4
         ce(5,is)  = ce(5,is) - fgp1*ua(5,is)/ua(1,is)
         ce(6,is)  = ce(6,is) - fgp1*ua(6,is)/ua(1,is)
850   continue  
c
c     vertices on the downstream boundary 
c
      do 950 i=1,nslog4
         is=log4(i)
         cap       = vno(is)
         capd1     = vnox(is)
         capd2     = vnoy(is)
         u         = ua(2,is)/ua(1,is)
         v         = ua(3,is)/ua(1,is)
         c         = sqrt(gam*pres(is)/ua(1,is))
         usc       = 1.0/c
         uq41      = 0.5*(u*u + v*v)*usc + gam4*c
         uq42      = capd1*u + capd2*v
c
c  constant invariants of riemann for u-c
c
         pstar=pout
         xroout=ua(1,is)*(pstar/pres(is))**(1./gam)
         unorm=sqrt(u*u+v*v)
         unorm=unorm+2./gam1*(c-sqrt(gam*pstar/xroout))
         xruxout=ruxout*unorm*xroout
         xruyout=ruyout*unorm*xroout
         xeout=pstar/gam1+xroout*unorm*unorm*0.5
c
c     computation of a-(wi).winf
c
         vp1       = (capd1*u + capd2*v)*cap 
         vp3       = min(vp1 + c*cap, 0.0)
         vp4       = min(vp1 - c*cap, 0.0)
         vp1       = min(vp1        , 0.0)
         fdc1      = vp1*(xroout + 
     &                    gam1*usc*usc*(-0.5*(u*u + v*v)*xroout +
     &                                  u*xruxout + v*xruyout - xeout))
         fdc2      = vp1*((capd1*v - capd2*u)*xroout + 
     &                    capd2*xruxout - capd1*xruyout)
         fdc3      = vp3*(0.5*(-uq42*xroout  + 
     &                         capd1*xruxout + capd2*xruyout) +
     &                    0.5*gam1*usc*(0.5*(u*u + v*v)*xroout  - 
     &                                  u*xruxout - v*xruyout + xeout))
         fdc4      = vp4*(0.5*(uq42*xroout   - 
     &                         capd1*xruxout - capd2*xruyout) + 
     &                    0.5*gam1*usc*(0.5*(u*u + v*v)*xroout  - 
     &                                  u*xruxout - v*xruyout + xeout))
         fgp1      = fdc1 + usc*(fdc3 + fdc4)
         fgp2      = u*fdc1 + capd2*fdc2  + 
     &               (u*usc + capd1)*fdc3 + (u*usc - capd1)*fdc4
         fgp3      = v*fdc1 - capd1*fdc2  + 
     &               (v*usc + capd2)*fdc3 + (v*usc - capd2)*fdc4
         fgp4    = 0.5*(u*u + v*v)*fdc1 + (capd2*u - capd1*v)*fdc2 +
     &               (uq41 + uq42)*fdc3   + (uq41 - uq42)*fdc4
c
         ce(1,is)  = ce(1,is) - fgp1
         ce(2,is)  = ce(2,is) - fgp2
         ce(3,is)  = ce(3,is) - fgp3
         ce(4,is)  = ce(4,is) - fgp4
         ce(5,is)  = ce(5,is) - fgp1*ua(5,is)/ua(1,is)
         ce(6,is)  = ce(6,is) - fgp1*ua(6,is)/ua(1,is)
c
c     computation of a+(wi).wi
c
         vp1       = (capd1*u + capd2*v)*cap 
         vp3       = max(vp1 + c*cap, 0.0)
         vp4       = max(vp1 - c*cap, 0.0)
         vp1       = max(vp1        , 0.0)
         fdc1      = vp1*(ua(1,is) +
     &                    gam1*usc*usc*(-0.5*(u*u + v*v)*ua(1,is)    +
     &                    u*ua(2,is) + v*ua(3,is) - ua(4,is)))
         fdc2      = vp1*((capd1*v - capd2*u)*ua(1,is) +
     &                    capd2*ua(2,is) - capd1*ua(3,is))
         fdc3      = vp3*(0.5*(-uq42*ua(1,is) +
     &                    capd1*ua(2,is) + capd2*ua(3,is)) +
     &                    0.5*gam1*usc*(0.5*(u*u + v*v)*ua(1,is)     -
     &                                  u*ua(2,is) - v*ua(3,is)      +
     &                                  ua(4,is)))
         fdc4   = vp4*(0.5*(uq42*ua(1,is)  -
     &                         capd1*ua(2,is) - capd2*ua(3,is))      +
     &                    0.5*gam1*usc*(0.5*(u*u + v*v)*ua(1,is)     -
     &                                  u*ua(2,is) - v*ua(3,is)      +
     &                                  ua(4,is)))
         fgp1      = fdc1 + usc*(fdc3 + fdc4)
         fgp2      = u*fdc1 + capd2*fdc2  +
     &               (u*usc + capd1)*fdc3 + (u*usc - capd1)*fdc4
         fgp3      = v*fdc1 - capd1*fdc2  +
     &               (v*usc + capd2)*fdc3 + (v*usc - capd2)*fdc4
         fgp4    = 0.5*(u*u + v*v)*fdc1 + (capd2*u - capd1*v)*fdc2 +
     &               (uq41 + uq42)*fdc3   + (uq41 - uq42)*fdc4
c
         ce(1,is)  = ce(1,is) - fgp1
         ce(2,is)  = ce(2,is) - fgp2
         ce(3,is)  = ce(3,is) - fgp3
         ce(4,is)  = ce(4,is) - fgp4
         ce(5,is)  = ce(5,is) - fgp1*ua(5,is)/ua(1,is)
         ce(6,is)  = ce(6,is) - fgp1*ua(6,is)/ua(1,is)

950   continue  
c      
      do 1000 i=1,nslog6
         is=log6(i)
         ce(1,is)  = 0.0
         ce(2,is)  = 0.0
         ce(3,is)  = 0.0
         ce(4,is)  = 0.0
         ce(5,is)  = 0.0
         ce(6,is)  = 0.0
1000   continue

      return
      end
