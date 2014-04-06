      subroutine flucin
c     -----------------------------------------------------------------
c     edgewise convective fluxes computation using the kinetic boltzmann
c     flux.
c     the nodal gradients are computed using a beta-combination 
c     of centered and hermitian (half-upwind) gradients
c
c     bijan mohammadi, inria-menusin
c     -----------------------------------------------------------------
      include 'param2d.h'
c     -----------------------------------------------------------------
c     local variables definition
      integer  nsg    , nubo1  , nubo2
      real dpm(4),dpor(4),dpex(4),aux1(4),aux2(4),
     1     gradi(4),gradj(4)
c
      gam4m5       = (2.-gam)/gam1
      ra3          = sqrt(3.0)
      alp          = 1./2./ra3
c
      beta2 = beta
      beta3 = 0.5*(1.0 - 2.0*beta)
      if(nordre.eq.2) then
      beta2 = 2.0*beta
      beta3 = 1.0 - 2.0*beta
      endif
      e2           = 1.e-6
c
c     loop on global list of edges
c
      do 500 nsg=1,nseg 
c
         nubo1     = nubo(1,nsg)
         nubo2     = nubo(2,nsg)
c
         xnn       =-vnocl(1,nsg)
         ynn       =-vnocl(2,nsg)
         rnn       = vnocl(3,nsg) 
c
         aix       = coor(1,nubo2)-coor(1,nubo1)
         aiy       = coor(2,nubo2)-coor(2,nubo1) 
c
         uas11     = ua(1,nubo1)
         uas12     = ua(1,nubo2)
         uas21     = ua(2,nubo1)
         uas22     = ua(2,nubo2)
         uas31     = ua(3,nubo1)
         uas32     = ua(3,nubo2)
         uas41     = pres(nubo1)
         uas42     = pres(nubo2)
c
         if(nordre.eq.1)  goto 1234
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
c
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
c rotation
c
         uas210 = uas21
         uas21  = xnn*uas210+ynn*uas31
         uas31  =-ynn*uas210+xnn*uas31
         uas41  = uas41/uas11
         uas220 = uas22
         uas22  = xnn*uas220+ynn*uas32
         uas32  =-ynn*uas220+xnn*uas32
         uas42  = uas42/uas12
c
         rat1 = sqrt(uas41)
         ext1 = min(ra3,max(-ra3,-uas21/rat1))
c
         a01  = alp*(ra3-ext1)
         a11  = alp*(ra3**2-ext1**2)/2.
         a21  = alp*(ra3**3-ext1**3)/3.
         a31  = alp*(ra3**4-ext1**4)/4.
c
         flu11= uas11*(uas21*a01+rat1*a11)
         flu21= uas11*(uas21*uas21*a01+2.*uas21*rat1*a11
     1                 +uas41*a21)
         flu31= uas31*flu11
         flu41= uas11*(uas21*(uas21*uas21+uas31*uas31+uas41)*a01
     1         +rat1*(3.*uas21*uas21+uas31*uas31+uas41)*a11
     1         +3.*uas21*uas41*a21+rat1*uas41*a31)*0.5
     1         +gam4m5*uas41*flu11
c
         rat2 = sqrt(uas42)
         ext2 = min(ra3,max(-ra3,-uas22/rat2))
c
         a01  = alp*(ra3-ext2)
         a11  = alp*(ra3**2-ext2**2)/2.
         a21  = alp*(ra3**3-ext2**3)/3.
         a31  = alp*(ra3**4-ext2**4)/4.
c
         a02  =1.0 -a01
         a12  =-a11 
         a22  =1.0 -a21
         a32  =-a31 
c
         flu12= uas12*(uas22*a02+rat2*a12)
         flu22= uas12*(uas22*uas22*a02+2.*uas22*rat2*a12
     1                 +uas42*a22)
         flu32= uas32*flu12
         flu42= uas12*(uas22*(uas22*uas22+uas32*uas32+uas42)*a02
     1         +rat2*(3.*uas22*uas22+uas32*uas32+uas42)*a12
     1         +3.*uas22*uas42*a22+rat2*uas42*a32)*0.5
     1         +gam4m5*uas42*flu12
c
c rotation inverse
c
         flu210 = flu21
         flu21  = xnn*flu210-ynn*flu31
         flu31  = ynn*flu210+xnn*flu31
         flu220 = flu22
         flu22  = xnn*flu220-ynn*flu32
         flu32  = ynn*flu220+xnn*flu32
c
         flu11=(flu11 + flu12)*rnn
         flu12=(flu21 + flu22)*rnn
         flu13=(flu31 + flu32)*rnn
         flu14=(flu41 + flu42)*rnn
c
         fluro(nsg)=-flu11
c
         ce(1,nubo1) = ce(1,nubo1) - flu11
         ce(2,nubo1) = ce(2,nubo1) - flu12
         ce(3,nubo1) = ce(3,nubo1) - flu13 
         ce(4,nubo1) = ce(4,nubo1) - flu14 
c
         ce(1,nubo2) = ce(1,nubo2) + flu11 
         ce(2,nubo2) = ce(2,nubo2) + flu12 
         ce(3,nubo2) = ce(3,nubo2) + flu13 
         ce(4,nubo2) = ce(4,nubo2) + flu14 
c
500   continue
c
      return
      end
