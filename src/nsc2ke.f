      program nsc2ke
c     -----------------------------------------------------------------
c
c     2D and AXI euler and navier-stokes equations solver 
c
c     explicit multi-steps time integration process
c     upwind schemes and linear interpolation method for 
c     the computation of the convective fluxes using a finite 
c     volume formulation.
c     classical central galerkin p1-finite element method 
c     for the computation of the diffusive fluxes
c
c     k-epsilon turbulence model with two-layer approach or wall laws
c
c     Bijan Mohammadi-Stephane Lanteri
c     INRIA Domaine de Voluceau 78153 Le Chesnay, France
c
c     release 1.0 Jan. 1994
c     Copyright(C) 1994 Bijan Mohammadi-Stephane Lanteri
c
c     Please send bugs and comments to bijan.mohamadi@inria.fr
c
c     -----------------------------------------------------------------
      include 'param2d.h'
c     -----------------------------------------------------------------
c     reading the application dependent parameters
c
      call config
c
c     loading the global mesh data structure
c
      call mailla
c
      call geowall
c
c     computation of the control volume and triangle areas
c
      call aires
c
c     computation of the control volume altitudes
c
      call clhaut
c
c     construction of the mesh edges 
c
      call seg2d
c
      call axigeo
c
c     initializing the physical solution 
c
      call init(kt0,t0)
c
      kt                     = kt0
      t                      = t0
      eps                    = 1.0e-10
c
c     swapping the initial physical solution 
c
      do 2000 is=1,ns
      do 2000 ivar=1,nvar
         un(ivar,is) = ua(ivar,is)
2000  continue
c
c     computation of the initial pressures
c
      call calprc
c
c     saving the initial physical solution 
c
      call result(dt)
c
      open(30,file="RESIDUAL",status="unknown")
      rewind(30)
c
c     physical time loop
c
100   kt                     = kt + 1

      call caldtl(dtmin,dt)
c
c     runge-kutta time integration process
c
      do 51 ialpha=1,irk

          call resexp(ialpha)
c
c        computation of the pressures
c
         call calprc
c
c        computation of the residual
c
         if (ialpha .eq. 1) then
               som           = 0.0
            do 115 is=1,ns
               som           = som + ce(1,is)*ce(1,is)
     1                             + ce(2,is)*ce(2,is)
     1                             + ce(3,is)*ce(3,is)
     1                             + ce(4,is)*ce(4,is)
115         continue
            som             = sqrt(som)/(float(ns))
            if ((kt - kt0) .eq. 1) som0 = som
            write(30,*) kt,som/som0
c
         endif                   
c
51    continue
c
c     swapping of the new physical solution
c
      do 200 is=1,ns
      do 200 ivar=1,nvar
         un(ivar,is)  = ua(ivar,is)
200   continue
c
         if(mod(kt,ifre).eq.0) then
         print *,'  kt  =  ',kt,'  residu  =  ',som/som0
         call result(dt)
         endif
c
c     testing for the end of the physical time loop
c     saving the final physical solution 
c
      if (kt  .eq. ktmax) then
         print*,'end of execution : maximal number of time steps : ',kt
         goto 300
      endif
      if (som .lt. resf) then
         print*,'end of execution : minimal residual : ',som
         goto 300
      endif
      if (abs(t - tmax) .lt. eps) then
         print*,'end of execution : maximal physical time: ',t
         goto 300
      endif
      goto 100
300   continue
      close(30)
         if(mod(kt,ifre).ne.0) then
         print *,'  kt  =  ',kt,'  residu  =  ',som
         call result(dt)
         endif
      stop
      end
