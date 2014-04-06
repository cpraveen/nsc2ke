      subroutine keps2d
c     -----------------------------------------------------------------
      include 'param2d.h'
c     -----------------------------------------------------------------
c     local variables definition
      real    dbxx(3)  , dbyy(3), dxt(2), dyt(2),
     &        uph(2,3) 
c
      us6               = 1.0/6.0
      us3               = 1.0/3.0
      us43              = 1.0/4.0/3.0
      gampr             = gam/pr
      gamprt            = gam/prt
      cmu               = 0.09
      cep               = 1./1.3
c
c     loop on global list of triangles
c
      do 1000 jt=1,nt
c
c        scatter operation 
c        getting the physical solutions for the three vertices 
c        of the current triangle 
c
         nubo1          = nu(1,jt)
         nubo2          = nu(2,jt)
         nubo3          = nu(3,jt)
c
         uph(1,1)       = ua(5,nubo1)
         uph(2,1)       = ua(6,nubo1)
c
         uph(1,2)       = ua(5,nubo2)
         uph(2,2)       = ua(6,nubo2)
c
         uph(1,3)       = ua(5,nubo3)
         uph(2,3)       = ua(6,nubo3)
c
c        computation of the p1-gradients
c
         x1             = coor(1,nubo1)
         y1             = coor(2,nubo1)
         x2             = coor(1,nubo2)
         y2             = coor(2,nubo2)
         x3             = coor(1,nubo3)
         y3             = coor(2,nubo3)
c
         dbxx(1)        = y2 - y3
         dbxx(2)        = y3 - y1
         dbxx(3)        = y1 - y2
         dbyy(1)        = x3 - x2
         dbyy(2)        = x1 - x3
         dbyy(3)        = x2 - x1
c                  
            dxt(1)      = uph(1,1)*dbxx(1) +
     &                    uph(1,2)*dbxx(2) +
     &                    uph(1,3)*dbxx(3)
            dyt(1)      = uph(1,1)*dbyy(1) +
     &                    uph(1,2)*dbyy(2) +
     &                    uph(1,3)*dbyy(3)
            dxt(2)      = uph(2,1)*dbxx(1) +
     &                    uph(2,2)*dbxx(2) +
     &                    uph(2,3)*dbxx(3)
            dyt(2)      = uph(2,1)*dbyy(1) +
     &                    uph(2,2)*dbyy(2) +
     &                    uph(2,3)*dbyy(3)
c
       xmlam=us3*(reylam(nubo1)+reylam(nubo2)+reylam(nubo3))
       xmtur=us3*(reyturb(nubo1)+reyturb(nubo2)+reyturb(nubo3))
       xmtot=xmlam+xmtur
c        
         aitt           = 0.25*airta(jt)
c
         ce(5,nubo1)    = ce(5,nubo1) - aitt*xmtot*
     &                    (dbxx(1)*dxt(1) + dbyy(1)*dyt(1))
         ce(6,nubo1)    = ce(6,nubo1) - aitt*(xmlam+cep*xmtur)*
     &                    (dbxx(1)*dxt(2) + dbyy(1)*dyt(2))
c
         ce(5,nubo2)    = ce(5,nubo2) - aitt*xmtot*
     &                    (dbxx(2)*dxt(1) + dbyy(2)*dyt(1))
         ce(6,nubo2)    = ce(6,nubo2) - aitt*(xmlam+cep*xmtur)*
     &                    (dbxx(2)*dxt(2) + dbyy(2)*dyt(2))
c
         ce(5,nubo3)    = ce(5,nubo3) - aitt*xmtot*
     &                    (dbxx(3)*dxt(1) + dbyy(3)*dyt(1))
         ce(6,nubo3)    = ce(6,nubo3) - aitt*(xmlam+cep*xmtur)*
     &                    (dbxx(3)*dxt(2) + dbyy(3)*dyt(2))
c
1000  continue
c
       do 500 nsg=1,nseg
         nubo1     = nubo(1,nsg)
         nubo2     = nubo(2,nsg)
c
         uas51     = ua(5,nubo1)
         uas61     = ua(6,nubo1)
         uas52     = ua(5,nubo2)
         uas62     = ua(6,nubo2)
c
         sgn=-sign(.5,fluro(nsg))
c
         ce(5,nubo1)    = ce(5,nubo1) + fluro(nsg)*
     1                   ( (.5+sgn)*uas51+(.5-sgn)*uas52 )
         ce(6,nubo1)    = ce(6,nubo1) + fluro(nsg)*
     1                   ( (.5+sgn)*uas61+(.5-sgn)*uas62 )
         ce(5,nubo2)    = ce(5,nubo2) - fluro(nsg)*
     1                   ( (.5+sgn)*uas51+(.5-sgn)*uas52 )
         ce(6,nubo2)    = ce(6,nubo2) - fluro(nsg)*
     1                   ( (.5+sgn)*uas61+(.5-sgn)*uas62 )
c
500     continue

         if(ilaw.ne.0) call ke_law
         if(ilaw.eq.0) call ke_two
c
      return
      end
