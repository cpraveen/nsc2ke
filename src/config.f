      subroutine config
c     -----------------------------------------------------------------
c     reading the application dependent parameters
c     -----------------------------------------------------------------
      include 'param2d.h'
c     -----------------------------------------------------------------
c
c     perferct gaz constant 
c
      gam     = 1.4
      gam1    = gam - 1.0
      gam4    = 1.0/gam1
c
c     prandtl's number 
c
      pr      = 0.72
      prt     = 0.9
c
cccccccccccccccccccccccccccccccccccccccccccc
      open(1,file="DATA",status='old')
      rewind(1)
      read(1, *) iaxi
      read(1, *) ivis
      read(1, *) reynolds
      read(1, *) froud
      read(1, *) xmach
      read(1, *) xpcoef
      read(1, *) iecc
      read(1, *) tinfd
      read(1, *) tbrd1
      read(1, *) tetadeg
      read(1, *) iflux
      read(1, *) nordre
      read(1, *) iloc
      read(1, *) cfl
      read(1, *) ktmax
      read(1, *) ifre
      read(1, *) tmax 
      read(1, *) resf
      read(1, *) ncont
      read(1, *) 
      read(1, *) iturb
      read(1, *) ilaw
      read(1, *) delta
      read(1, *) ncont1
      read(1, *) xtmin,xtmax,ytmin,ytmax
      close(1)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      print *,'*********************************************'
      print *,'        NSC2KE : RELEASE 1.0, 1994           '
      print *,'*********************************************'
      if(iaxi.eq.1) print *,'AXI SYMMETRIC        '
      if(iaxi.eq.0) print *,'2D              '
      if(ivis.eq.1) print *,'Navier-Stokes computation '
      if(ivis.eq.0) print *,'Euler computation      '
      print *,'*********************************************'
      if(iturb.eq.1) then
          print *,'k-epsilon turbulence model      '
      if(ilaw.eq.1) print *,'classical wall-law  technique  '
      if(ilaw.eq.0) print *,'two-layer technique            '
      endif
      print *,'*********************************************'
c
c     mach number, angle of attack
c
      if(iecc.eq.1) print *,'adiabatic walls '
      if(iecc.eq.2) print *,'isothermal walls '
      print*,'free stream mach number  : ', xmach
      print*,'angle of attack          : ', tetadeg
c
      if (ivis .eq. 1) then
      print*,'reynolds number          : ', reynolds
      reynolds = 1./reynolds
      endif
      print*,'Froude number            : ', froud
c
        pin     = 1.0/(gam*xmach**2)
        tinf    = 1.0/(gam*gam1*xmach*xmach)
c
      tbrd2=tinf*tbrd1/tinfd
c
c     free stream uniform solution 
c     external flow around a body
c
      roin    = 1.0
      teta    = tetadeg*3.141592654/180.0
      uxin    = cos(teta)
      uyin    = sin(teta)
      ruxin   = roin*uxin
      ruyin   = roin*uyin
      ein     = 0.5*roin*(uxin*uxin+uyin*uyin)+pin/gam1
      roout   = roin
      uxout   = uxin
      uyout   = uyin
      ruxout  = ruxin
      ruyout  = ruyin
      pout    = pin*xpcoef
      eout    = 0.5*roout*(uxout**2+uyout**2)+pout/gam1
      xkin    = 1.e-5
      xein    = 1.e-5
c
      print*,'courant number           : ', cfl
c
      if(iloc.ne.0) print *,'local time stepping'
c
c     upwinding parameter 
c
      beta=0.33333
      beta1   = 1.0 - beta
c
      if (iflux .eq. 1) print*,'Roe''s scheme'
      if (iflux .eq. 2) print*,'Osher''s scheme'
      if (iflux .eq. 3) print*,'Kinetic''s scheme'
c
c     runge-kutta time integration process parameters 
c
c      irk=3
c      alpha(1)=0.3333
c      alpha(2)=0.5
c      alpha(3)=1.0 
ccccccccccccccccccccccccccccccccc
       irk=4
       alpha(1)=0.11
       alpha(2)=0.2766
       alpha(3)=0.5
       alpha(4)=1.0
cccccccccccccccccccccccccccccccc        
      print *,'*********************************************'
      close(1)
      return
      end
