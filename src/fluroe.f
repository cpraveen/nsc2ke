      subroutine fluroe
c     -----------------------------------------------------------------
c     edgewise convective fluxes computation using roe's
c     approximate riemann solver 
c     the nodal gradients are computed using a beta-combination 
c     of centered and hermitian (half-upwind) gradients
c     -----------------------------------------------------------------
      include 'param2d.h'
c     -----------------------------------------------------------------
c     local variables definition
      integer  nsg    , nubo1  , nubo2
      real invas1, invcnt
      real dpm(4),dpor(4),dpex(4),aux1(4),aux2(4),
     1     gradi(4),gradj(4)

c
      usgam1=1.0/gam1
      gsgam1=gam/gam1     
      beta2 = beta
      beta3 = 0.5*(1.0 - 2.0*beta)
      if(nordre.eq.2) then
      beta2 = 2.0*beta
      beta3 = 1.0 - 2.0*beta
      endif
      e2 = 1.e-16
c
c     loop on global list of edges
c
      do 500 nsg=1,nseg 
c
c        local indexing of the vertices of the current edge
c
         nubo1     = nubo(1,nsg)
         nubo2     = nubo(2,nsg)
c
         aix       = coor(1,nubo2)-coor(1,nubo1)
         aiy       = coor(2,nubo2)-coor(2,nubo1) 
c
c        indirect addressing on vertices physical states
c
         uas11     = ua(1,nubo1)
         uas21     = ua(2,nubo1)
         uas31     = ua(3,nubo1)
         uas41     = pres(nubo1)
         uas12     = ua(1,nubo2)
         uas22     = ua(2,nubo2)
         uas32     = ua(3,nubo2)
         uas42     = pres(nubo2)
c
         if(nordre.eq.1) goto 1234
c
         gradI(1)  = beta2*(aix*dx(1,nubo1) + aiy*dy(1,nubo1)) +
     &               beta3*(uas12 - uas11)
c
         gradI(2)  = beta2*(aix*dx(2,nubo1) + aiy*dy(2,nubo1)) +
     &               beta3*(uas22 - uas21)
c
         gradI(3)  = beta2*(aix*dx(3,nubo1) + aiy*dy(3,nubo1)) +
     &               beta3*(uas32 - uas31)
c
         gradI(4)  = beta2*(aix*dx(4,nubo1) + aiy*dy(4,nubo1)) +
     &               beta3*(uas42 - uas41)
c
         gradJ(1)  = beta2*(aix*dx(1,nubo2) + aiy*dy(1,nubo2)) +
     &               beta3*(uas12 - uas11)
c
         gradJ(2)  = beta2*(aix*dx(2,nubo2) + aiy*dy(2,nubo2)) +
     &               beta3*(uas22 - uas21)
c
         gradJ(3)  = beta2*(aix*dx(3,nubo2) + aiy*dy(3,nubo2)) +
     &               beta3*(uas32 - uas31)
c
         gradJ(4)  = beta2*(aix*dx(4,nubo2) + aiy*dy(4,nubo2)) +
     &               beta3*(uas42 - uas41)
c

         if(nordre.eq.2) then
c
         uas11     = uas11 + 0.5*gradI(1)
         uas21     = uas21 + 0.5*gradI(2)
         uas31     = uas31 + 0.5*gradI(3)
         uas41     = uas41 + 0.5*gradI(4)
c
         uas12     = uas12 - 0.5*gradJ(1)
         uas22     = uas22 - 0.5*gradJ(2)
         uas32     = uas32 - 0.5*gradJ(3)
         uas42     = uas42 - 0.5*gradJ(4)
c
         elseif(nordre.eq.3) then
c
         dpm(1)    =-(uas12 - uas11)
         dpex(1)   =-4.0*gradI(1) - dpm(1)
         aux1(1)   = 0.25*(1.0+  SIGN(1.0, dpex(1)*dpm(1)))
         dpor(1)   =-4.0*gradJ(1) - dpm(1)
         aux2(1)   = 0.25*(1.0 + SIGN(1.0, dpor(1)*dpm(1)))
c
         dpm(2)    =-(uas22 - uas21)
         dpex(2)   =-4.0*gradI(2) - dpm(2)
         aux1(2)   = 0.25*(1.0 + SIGN(1.0, dpex(2)*dpm(2)))
         dpor(2)   =-4.0*gradJ(2) - dpm(2)
         aux2(2)   = 0.25*(1.0 + SIGN(1.0, dpor(2)*dpm(2)))
c
         dpm(3)    =-(uas32 - uas31)
         dpex(3)   =-4.0*gradI(3) - dpm(3)
         aux1(3)   = 0.25*(1.0 + SIGN(1.0, dpex(3)*dpm(3)))
         dpor(3)   =-4.0*gradJ(3) - dpm(3)
         aux2(3)   = 0.25*(1.0 + SIGN(1.0, dpor(3)*dpm(3)))
c
         dpm(4)    =-(uas42 - uas41)
         dpex(4)   =-4.0*gradI(4) - dpm(4)
         aux1(4)   = 0.25*(1.0 + SIGN(1.0, dpex(4)*dpm(4)))
         dpor(4)   =-4.0*gradJ(4) - dpm(4)
         aux2(4)   = 0.25*(1.0 + SIGN(1.0, dpor(4)*dpm(4)))
c
         gradI(1)  = aux1(1)*
     &               ((dpex(1)*dpex(1) + e2)*dpm(1)+
     &                (dpm(1)*dpm(1)   + e2)*dpex(1))/
     &               (dpex(1)*dpex(1)  + dpm(1)*dpm(1) + 2.0*e2)
         gradJ(1)  = aux2(1)*
     &               ((dpor(1)*dpor(1) + e2)*dpm(1)+
     &                (dpm(1)*dpm(1)   + e2)*dpor(1))/
     &               (dpor(1)*dpor(1)  + dpm(1)*dpm(1) + 2.0*e2)
c
         gradI(2)  = aux1(2)*
     &               ((dpex(2)*dpex(2) + e2)*dpm(2)+
     &                (dpm(2)*dpm(2)   + e2)*dpex(2))/
     &               (dpex(2)*dpex(2)  + dpm(2)*dpm(2) + 2.0*e2)
         gradJ(2)  = aux2(2)*
     &               ((dpor(2)*dpor(2) + e2)*dpm(2)+
     &                (dpm(2)*dpm(2)   + e2)*dpor(2))/
     &               (dpor(2)*dpor(2)  + dpm(2)*dpm(2) + 2.0*e2)
         gradI(3)  = aux1(3)*
     &               ((dpex(3)*dpex(3) + e2)*dpm(3)+
     &                (dpm(3)*dpm(3)   + e2)*dpex(3))/
     &               (dpex(3)*dpex(3)  + dpm(3)*dpm(3) + 2.0*e2)
         gradJ(3)  = aux2(3)*
     &               ((dpor(3)*dpor(3) + e2)*dpm(3)+
     &                (dpm(3)*dpm(3)   + e2)*dpor(3))/
     &               (dpor(3)*dpor(3)  + dpm(3)*dpm(3) + 2.0*e2)
c
         gradI(4)  = aux1(4)*
     &               ((dpex(4)*dpex(4) + e2)*dpm(4)+
     &                (dpm(4)*dpm(4)   + e2)*dpex(4))/
     &               (dpex(4)*dpex(4)  + dpm(4)*dpm(4) + 2.0*e2)
         gradJ(4)  = aux2(4)*
     &               ((dpor(4)*dpor(4) + e2)*dpm(4)+
     &                (dpm(4)*dpm(4)   + e2)*dpor(4))/
     &               (dpor(4)*dpor(4)  + dpm(4)*dpm(4) + 2.0*e2)
c
         uas11     = uas11 - gradI(1)
         uas21     = uas21 - gradI(2)
         uas31     = uas31 - gradI(3)
         uas41     = uas41 - gradI(4)
c
         uas12     = uas12 + gradJ(1)
         uas22     = uas22 + gradJ(2)
         uas32     = uas32 + gradJ(3)
         uas42     = uas42 + gradJ(4)
c
         endif
 
c
1234   continue
c
         invas1    = 1.0/uas11
c
         dua1      = sqrt(uas11)
         dua2      = sqrt(uas12)
         dua3      = 1.0/(dua1 + dua2)
c
c        enthalpy
c
         h1        = gsgam1*uas41*invas1 + 
     &               0.5*(uas21*uas21 + uas31*uas31)
         h2        = gsgam1*uas42/uas12  +
     &               0.5*(uas22*uas22 + uas32*uas32)
c
c        roe's mean values computation
c
         ucent     = (uas21*dua1 + uas22*dua2)*dua3
         vcent     = (uas31*dua1 + uas32*dua2)*dua3
         hcent     = (h1*dua1 + h2*dua2)*dua3
         uvnc      = 0.5*(ucent*ucent + vcent*vcent)
         ccent     = sqrt(abs(gam1*(hcent-uvnc)))
         invcnt    = 1.0/ccent
c
         pres1     = uas41 
         uas41     = usgam1*uas41 + 0.5*uas11*
     &               (uas21*uas21 + uas31*uas31)
         uas42     = usgam1*uas42 + 0.5*uas12*
     &               (uas22*uas22 + uas32*uas32)
         uas21     = uas21*uas11
         uas22     = uas22*uas12
         uas31     = uas31*uas11
         uas32     = uas32*uas12
c
         dua1      = uas12 - uas11
         dua2      = uas22 - uas21
         dua3      = uas32 - uas31
         dua4      = uas42 - uas41
c
c        eigenvalues computation
c
         xnn       = vnocl(1,nsg)
         ynn       = vnocl(2,nsg)
         rnn       = vnocl(3,nsg)
c
         vp1       = rnn*(xnn*ucent + ynn*vcent)
         vp3       = vp1 + ccent*rnn
         vp4       = vp1 - ccent*rnn
         uq41      = uvnc*invcnt + usgam1*ccent
         uq42      = xnn*ucent + ynn*vcent
         rc1       = gam1*invcnt
         rc2       = rc1*invcnt
c
c        computation of the centered part of the flux
c
         flu11     = rnn*(xnn*uas21 + ynn*uas31)
         flu21     = xnn*rnn*pres1 + flu11*uas21*invas1
         flu31     = ynn*rnn*pres1 + flu11*uas31*invas1
         flu41     = h1*flu11
c
c        computation of the diffusive part of the flux
c
         fdc1      = max(vp1,0.0)*
     &               (dua1 + rc2*(-uvnc*dua1 + ucent*dua2 + 
     &                            vcent*dua3 - dua4)) 
         fdc2      = max(vp1,0.0)*
     &               ((xnn*vcent - ynn*ucent)*
     &                dua1 + ynn*dua2 - xnn*dua3)
         fdc3      = max(vp3,0.0)*
     &               (0.5*(-uq42*dua1 + xnn*dua2 + ynn*dua3) + 
     &                0.5*rc1*(uvnc*dua1 - ucent*dua2 - vcent*dua3 + 
     &                dua4))
         fdc4      = max(vp4,0.0)*
     &               (0.5*( uq42*dua1 - xnn*dua2 - ynn*dua3) +
     &                0.5*rc1*(uvnc*dua1 - ucent*dua2 - vcent*dua3 + 
     &                dua4))   
c
         duv1      = fdc1 + (fdc3 + fdc4)*invcnt
         duv2      = ucent*fdc1    + ynn*fdc2  +
     &               (ucent*invcnt + xnn)*fdc3 +
     &               (ucent*invcnt - xnn)*fdc4
         duv3      = vcent*fdc1    - xnn*fdc2  +
     &               (vcent*invcnt + ynn)*fdc3 + 
     &               (vcent*invcnt - ynn)*fdc4
         duv4      = uvnc*fdc1 + (ynn*ucent - xnn*vcent)* fdc2 +
     &               (uq41+uq42)*fdc3 + (uq41-uq42)*fdc4
c
c        gathering of the elementary fluxes into the global ones 
c
c
         fluro(nsg)=flu11+duv1
c
         ce(1,nubo1)    = ce(1,nubo1) + flu11 + duv1
         ce(2,nubo1)    = ce(2,nubo1) + flu21 + duv2
         ce(3,nubo1)    = ce(3,nubo1) + flu31 + duv3
         ce(4,nubo1)    = ce(4,nubo1) + flu41 + duv4
c
         ce(1,nubo2)    = ce(1,nubo2) - flu11 - duv1
         ce(2,nubo2)    = ce(2,nubo2) - flu21 - duv2
         ce(3,nubo2)    = ce(3,nubo2) - flu31 - duv3
         ce(4,nubo2)    = ce(4,nubo2) - flu41 - duv4
c
500   continue
c
      return
      end
