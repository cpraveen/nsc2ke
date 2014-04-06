      subroutine fluosh
c     ---------------------------------------------------------
c     edgewise convective fluxes computation using osher's
c     approximate riemann solver 
c     the nodal gradients are computed using a beta-combination 
c     of centered and hermitian (half-upwind) gradients
c     bijan mohammadi, INRIA
c     ---------------------------------------------------------
      include 'param2d.h'
c     ---------------------------------------------------------
c        
           real dpm(4),dpor(4),dpex(4),aux1(4),aux2(4),
     1     gradi(4),gradj(4)

      e2       =1.e-16
      beta2 = beta
      beta3 = 0.5*(1.0 - 2.0*beta)
      if(nordre.eq.2) then
      beta2 = 2.0*beta
      beta3 = 1.0 - 2.0*beta
      endif
      gg1      = gam/gam1
      GAMO     = (GAM-1.)*.5
      USG0     = 1./GAMO
      POW      = 1./(2.*GAM)
      COEFF    = GAM1/(GAM+1.)
c
c     loop on global list of edges
c
      do 500 nsg=1,nseg 
c
         xnn       =- vnocl(1,nsg)
         ynn       =- vnocl(2,nsg)
         rnn       =  vnocl(3,nsg) * 0.5
c 
         nubo1     = nubo(1,nsg)
         nubo2     = nubo(2,nsg)
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
      PROD1 = UAS21*xnn+UAS31*ynn   
      C1    = SQRT(GAM*UAS41/UAS11)
      PROD2 = UAS22*xnn+UAS32*ynn
      C2    = SQRT(GAM*UAS42/UAS12)
C
C    -- FLUXOSHER = F(NUBO1)+F(NUBO2)
C 
      FM1 = PROD1*UAS11+PROD2*UAS12
      FM2 = UAS11*UAS21*PROD1+UAS41*xnn+
     &      UAS12*UAS22*PROD2+UAS42*xnn
      FM3 = UAS11*UAS31*PROD1+UAS41*ynn+
     &      UAS12*UAS32*PROD2+UAS42*ynn
      FM4 = (GG1*UAS41+.5*UAS11*
     & (UAS21*UAS21+UAS31*UAS31))*PROD1+
     &      (GG1*UAS42+.5*UAS12*
     & (UAS22*UAS22+UAS32*UAS32))*PROD2
C
C  II/ CALCUL DES QUANTITES AUX POINTS LIMITES POUR LE C.L.D
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      DNUM       = GAMO*(PROD2-PROD1)+C1+C2
      DEN        = ((UAS41/UAS42)**POW)*SQRT(UAS12/UAS11)
      ULD11      = ((DNUM/((1.+1./DEN)*C1))**USG0)*UAS11
      ULD21      = ((DNUM/((1.+   DEN)*C2))**USG0)*UAS12
      ULD4       = ((ULD11/UAS11)**GAM)*UAS41
C
      AUX        = USG0*(C1-SQRT(GAM*ULD4/ULD11))
      ULD12      = UAS21-AUX*xnn
      ULD13      = UAS31-AUX*ynn
      ULD22      = UAS22-(PROD2-PROD1+AUX)*xnn
      ULD23      = UAS32-(PROD2-PROD1+AUX)*ynn
C
C  III/ CALCUL DES FLUX AUX POINTS LIMITES DU C.L.D
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     -- FLULD = FLUX LINEAIREMENT DEGENERE
c
      PRODLD       = ULD12*xnn+ULD13*ynn
c
      FLULD11      = ULD11*PRODLD
      FLULD12      = ULD11*ULD12*PRODLD+ULD4*xnn
      FLULD13      = ULD11*ULD13*PRODLD+ULD4*ynn
      FLULD14      = (gg1*ULD4+.5*(ULD12*ULD12+ulD13*ULD13)*ULD11)
     & *PRODLD
c
      FLULD21      = ULD21*PRODLD
      FLULD22      = ULD21*ULD22*PRODLD+ULD4*xnn
      FLULD23      = ULD21*ULD23*PRODLD+ULD4*ynn
      FLULD24      = (gg1*ULD4+.5*(ULD22*ULD22+ULD23*ULD23)*ULD21)
     & *PRODLD
c
C                           -- FLUOSHER = FLUOSHER - FLULD
c
      AUX     =ABS(PRODLD)
      FM1=FM1-AUX*(ULD21-ULD11)
      FM2=FM2-AUX*(ULD21*ULD22-ULD11*ULD12)
      FM3=FM3-AUX*(ULD21*ULD23-ULD11*ULD13)
      FM4=FM4-AUX*.5*(ULD21*(ULD22*ULD22+ULD23*ULD23)-
     &                ULD11*(ULD12*ULD12+ULD13*ULD13))
C
C  IV/ CALCUL DES QUANTITES AUX POINTS LIMITES POUR LES 2 C.V.N.L.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      AUX        = (PROD2+C2*USG0)/C2
      AX         = .5*(1.+SIGN(1.,AUX))
      UVNL21     = AX*((COEFF*ABS(AUX))**USG0)*UAS12
      UVNL24     = AX*((UVNL21/UAS12)**GAM)*UAS42
      AUX        = USG0*(C2-SQRT(GAM*UVNL24/(UVNL21+(1.-AX))))
      UVNL22     = AX*(UAS22+AUX*xnn)
      UVNL23     = AX*(UAS32+AUX*ynn)
C
      AUX        = (C1*USG0-PROD1)/C1
      AX         = .5*(1.+SIGN(1.,AUX))
      UVNL11     = AX*((COEFF*ABS(AUX))**USG0)*UAS11
      UVNL14     = AX*((UVNL11/UAS11)**GAM)*UAS41
      AUX        = USG0*(C1-SQRT(GAM*UVNL14/(UVNL11+(1.-AX))))
      UVNL12     = AX*(UAS21-AUX*xnn)
      UVNL13     = AX*(UAS31-AUX*ynn)
c                         -- FLUVN = FLUX VRAIMENT NON LINEAIRE
      PROVL1     = UVNL12*xnn+UVNL13*ynn
      FLUVN11    = UVNL11*PROVL1
      FLUVN12    = UVNL11*UVNL12*PROVL1+UVNL14*xnn
      FLUVN13    = UVNL11*UVNL13*PROVL1+UVNL14*ynn
      FLUVN14    = (gg1*UVNL14+.5*(UVNL12*UVNL12+
     &             UVNL13*UVNL13)*UVNL11)*PROVL1
C
      PROVL2     = UVNL22*xnn+UVNL23*ynn
      FLUVN21    = UVNL21*PROVL2
      FLUVN22    = UVNL21*UVNL22*PROVL2+UVNL24*xnn
      FLUVN23    = UVNL21*UVNL23*PROVL2+UVNL24*ynn
      FLUVN24    = (gg1*UVNL24+.5*(UVNL22*UVNL22+
     &             UVNL23*UVNL23)*UVNL21)*PROVL2
C
      SI1        = SIGN(1.,PROD2 -C2)
      SI2        = SIGN(1.,PRODLD-SQRT(GAM*ULD4/ULD21))
      SI3        = SIGN(1.,PRODLD+SQRT(GAM*ULD4/ULD11))
      SI4        = SIGN(1.,PROD1 +C1)
c
C       -- FLUOSHER = FLUOSHER - FLUVN
c
      FM1=FM1-(UAS12*PROD2-FLUVN21)*SI1
     &       +(  FLULD21-FLUVN21)*SI2
     &       -(  FLULD11-FLUVN11)*SI3
     &       +(  UAS11*PROD1-FLUVN11)*SI4
      FM2=FM2-(UAS12*UAS22*PROD2+
     &         UAS42*xnn-FLUVN22)*SI1
     &       +(  FLULD22-FLUVN22)*SI2
     &       -(  FLULD12-FLUVN12)*SI3
     &       +(  UAS11*UAS21*PROD1+
     &           UAS41*xnn-FLUVN12)*SI4
      FM3=FM3-(UAS12*UAS32*PROD2+
     &         UAS42*ynn-FLUVN23)*SI1
     &       +(  FLULD23-FLUVN23)*SI2
     &       -(  FLULD13-FLUVN13)*SI3
     &       +(  UAS11*UAS31*PROD1+
     &           UAS41*ynn-FLUVN13)*SI4
      FM4=FM4-((gg1*UAS42+.5*(UAS22*UAS22+UAS32*UAS32)*
     &        UAS12)*PROD2-FLUVN24)*SI1
     &       +(    FLULD24-FLUVN24)*SI2
     &       -(    FLULD14-FLUVN14)*SI3
     &       +((gg1*UAS41+.5*(UAS21*UAS21+UAS31*UAS31)*
     &        UAS11)*PROD1-FLUVN14)*SI4

C
C 3/ REPORT SUR LES NOEUDS DES SEGMENTS
C======================================
C                              
      fluro(nsg)=-fm1*rnn
c 
      CE(1,NUBO1)=CE(1,NUBO1)-FM1*rnn
      CE(2,NUBO1)=CE(2,NUBO1)-FM2*rnn
      CE(3,NUBO1)=CE(3,NUBO1)-FM3*rnn
      CE(4,NUBO1)=CE(4,NUBO1)-FM4*rnn
C
      CE(1,NUBO2)=CE(1,NUBO2)+FM1*rnn
      CE(2,NUBO2)=CE(2,NUBO2)+FM2*rnn
      CE(3,NUBO2)=CE(3,NUBO2)+FM3*rnn
      CE(4,NUBO2)=CE(4,NUBO2)+FM4*rnn
c
500   continue
c
      return
      end
