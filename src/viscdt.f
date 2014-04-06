      subroutine viscdt
c     -----------------------------------------------------------------
c     computation of the hermitian nodal gradients 
c     computation of the viscous fluxes
c     -----------------------------------------------------------------
      include 'param2d.h'
c     -----------------------------------------------------------------
c     local variables definition
      real    dbxx(3)  , dbyy(3), dxt(4), dyt(4),
     &        uph(4,3) , um(2)  , ei(3) 
c
      us6               = 1.0/6.0
      us3               = 1.0/3.0
      us43              = 1.0/4.0/3.0
      gampr             = gam/pr
      gamprt            = gam/prt
      cmu               = 0.09
      cep               = 1./1.3
c
c     initializing the hermitian nodal gradients
c
      do 2 is=1,ns
         dx(1,is)       = 0.0
         dx(2,is)       = 0.0
         dx(3,is)       = 0.0
         dx(4,is)       = 0.0
         dy(1,is)       = 0.0
         dy(2,is)       = 0.0
         dy(3,is)       = 0.0
         dy(4,is)       = 0.0
2     continue
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
         uph(1,1)       = ua(1,nubo1)
         uph(2,1)       = ua(2,nubo1)
         uph(3,1)       = ua(3,nubo1)
         uph(4,1)       = pres(nubo1)
c
         uph(1,2)       = ua(1,nubo2)
         uph(2,2)       = ua(2,nubo2)
         uph(3,2)       = ua(3,nubo2)
         uph(4,2)       = pres(nubo2)
c
         uph(1,3)       = ua(1,nubo3)
         uph(2,3)       = ua(2,nubo3)
         uph(3,3)       = ua(3,nubo3)
         uph(4,3)       = pres(nubo3)
c
c        specific internal energy
c
         ei(1)          = gam4*uph(4,1)/uph(1,1)
         ei(2)          = gam4*uph(4,2)/uph(1,2)
         ei(3)          = gam4*uph(4,3)/uph(1,3)
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
         do 4 k=1,4
            dxt(k)      = uph(k,1)*dbxx(1) +
     &                    uph(k,2)*dbxx(2) +
     &                    uph(k,3)*dbxx(3)
            dyt(k)      = uph(k,1)*dbyy(1) +
     &                    uph(k,2)*dbyy(2) +
     &                    uph(k,3)*dbyy(3)
4        continue
c
         dx(1,nubo1)    = dx(1,nubo1) + dxt(1)
         dx(1,nubo2)    = dx(1,nubo2) + dxt(1)
         dx(1,nubo3)    = dx(1,nubo3) + dxt(1)
         dx(2,nubo1)    = dx(2,nubo1) + dxt(2)
         dx(2,nubo2)    = dx(2,nubo2) + dxt(2)
         dx(2,nubo3)    = dx(2,nubo3) + dxt(2)
         dx(3,nubo1)    = dx(3,nubo1) + dxt(3)
         dx(3,nubo2)    = dx(3,nubo2) + dxt(3)
         dx(3,nubo3)    = dx(3,nubo3) + dxt(3)
         dx(4,nubo1)    = dx(4,nubo1) + dxt(4)
         dx(4,nubo2)    = dx(4,nubo2) + dxt(4)
         dx(4,nubo3)    = dx(4,nubo3) + dxt(4)
c
         dy(1,nubo1)    = dy(1,nubo1) + dyt(1)
         dy(1,nubo2)    = dy(1,nubo2) + dyt(1)
         dy(1,nubo3)    = dy(1,nubo3) + dyt(1)
         dy(2,nubo1)    = dy(2,nubo1) + dyt(2)
         dy(2,nubo2)    = dy(2,nubo2) + dyt(2)
         dy(2,nubo3)    = dy(2,nubo3) + dyt(2)
         dy(3,nubo1)    = dy(3,nubo1) + dyt(3)
         dy(3,nubo2)    = dy(3,nubo2) + dyt(3)
         dy(3,nubo3)    = dy(3,nubo3) + dyt(3)
         dy(4,nubo1)    = dy(4,nubo1) + dyt(4)
         dy(4,nubo2)    = dy(4,nubo2) + dyt(4)
         dy(4,nubo3)    = dy(4,nubo3) + dyt(4)
c
         if (ivis .eq. 0) goto 1000
c
c        computation of the viscous fluxes
c        mean value of the velocity on the current triangle
c
         um(1)          = us3*(uph(2,1) + uph(2,2) + uph(2,3))
         um(2)          = us3*(uph(3,1) + uph(3,2) + uph(3,3))
c
       xmlam=us3*(reylam(nubo1)+reylam(nubo2)+reylam(nubo3))
       xmtur=us3*(reyturb(nubo1)+reyturb(nubo2)+reyturb(nubo3))
       xmtot=xmlam+xmtur
c        
         eix            = dbxx(1)*ei(1) + dbxx(2)*ei(2) +
     &                    dbxx(3)*ei(3)
         eiy            = dbyy(1)*ei(1) + dbyy(2)*ei(2) +
     &                    dbyy(3)*ei(3)

c        deformation tensor components
c
      r3  = (dyt(2) + dxt(3))*xmtot
      r2  = (2.0*us3*(2.0*dxt(2) - dyt(3)))*xmtot
      s3  = (2.0*us3*(2.0*dyt(3) - dxt(2)))*xmtot
      r4  = (um(1)*r2+um(2)*r3)+(gampr*xmlam+gamprt*xmtur)*eix
      s4  = (um(1)*r3+um(2)*s3)+(gampr*xmlam+gamprt*xmtur)*eiy
c
c        gathering of the elementary flux into the global one    
c                           
         aitt           = 0.25*airta(jt)
c
         ce(2,nubo1)    = ce(2,nubo1) - 
     &                    aitt*(dbxx(1)*r2 + dbyy(1)*r3)
         ce(3,nubo1)    = ce(3,nubo1) - 
     &                    aitt*(dbxx(1)*r3 + dbyy(1)*s3)
         ce(4,nubo1)    = ce(4,nubo1) - 
     &                    aitt*(dbxx(1)*r4 + dbyy(1)*s4)
c
         ce(2,nubo2)    = ce(2,nubo2) - 
     &                    aitt*(dbxx(2)*r2 + dbyy(2)*r3)
         ce(3,nubo2)    = ce(3,nubo2) - 
     &                    aitt*(dbxx(2)*r3 + dbyy(2)*s3)
         ce(4,nubo2)    = ce(4,nubo2) - 
     &                    aitt*(dbxx(2)*r4 + dbyy(2)*s4)
c
         ce(2,nubo3)    = ce(2,nubo3) - 
     &                    aitt*(dbxx(3)*r2 + dbyy(3)*r3)
         ce(3,nubo3)    = ce(3,nubo3) - 
     &                    aitt*(dbxx(3)*r3 + dbyy(3)*s3)
         ce(4,nubo3)    = ce(4,nubo3) - 
     &                    aitt*(dbxx(3)*r4 + dbyy(3)*s4)
c
1000  continue
c
c     completing the computation of the nodal gradients 
c
      do 110 is=1,ns
         ais            = us6/airs(is)
         dx(1,is)       = dx(1,is)*ais
         dx(2,is)       = dx(2,is)*ais
         dx(3,is)       = dx(3,is)*ais
         dx(4,is)       = dx(4,is)*ais
         dy(1,is)       = dy(1,is)*ais
         dy(2,is)       = dy(2,is)*ais
         dy(3,is)       = dy(3,is)*ais
         dy(4,is)       = dy(4,is)*ais
110   continue
c
      return
      end
