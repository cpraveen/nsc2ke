c     -----------------------------------------------------------------
      implicit real (a-h,o-z)
c     -----------------------------------------------------------------
c     parameter and common definitions for the 2d navier-stokes solver
c     -----------------------------------------------------------------
c     parameters for the global representation of the mesh 
c          nn    : maximum number of vertices
c          nnt   : maximum number of triangles
c          nnsg  : maximum number of edges 
c          nnfr  : maximum number of boundary edges 
      parameter (nn  = 20000,
     &           nnt = 2*nn,
     &           nnsg= 3*nn,  
     &           nnfr= 1500 ,
     &           maxp=20,maxt=20)    
      parameter(nvar=6)
c     -----------------------------------------------------------------
c     maximum number of steps in the runge-kutta time integration 
c     process
      parameter (irkmax= 4)
c     -----------------------------------------------------------------
c     data structures for the global representation of the mesh 
c          ns    : effective number of vertices
c          nt    : effective number of triangles
c          nseg  : effective number of edges 
      common/comge0/ns,nt,nseg,nfr
c     x and y coordinates
      common/comge1/coor(2,nn)
c     control volume and triangle areas
      common/comge2/airs(nn),airt(nnt)
c     edges connectivity table 
      common/comge3/nubo(2,nnsg)
c     boundary references, boundary edges list and 
c     triangles connectivity table
      common/comge4/logfr(nn),nufr(nnfr),nu(3,nnt)
c     control volume normals definition
      common/comge5/vnocl(3,nnsg)
c     boundary normals definition
      common/comge6/vnox(nn),vnoy(nn),vno(nn)
c     control volume altitude
      common/comge9/dthaut(nn)
c     -----------------------------------------------------------------
c     new and old solution 
c          un1, ua1 : density
c          un2, ua2 : horizontal mometum
c          un3, ua3 : vertical mometum
c          un4, ua4 : total energy per unit of volume
c          pres     : pressure                       
c          un5, ua6 : kinetic energy of turbulence k
c          un6, ua7 : epsilon = rate of dissiparion of k
      common/comso1/un(nvar,nn)
      common/comso2/ua(nvar,nn),pres(nn)
c     -----------------------------------------------------------------
c     time step strategy
      common/comdt0/iloc,isave
      common/comdt1/cfl,dtl(nn)
      common/comkt/t,kt,tmax,ktmax,resf,ifre
c     runge-kutta time integration process definition
      common/comrkk/irk,alpha(irkmax)
c     -----------------------------------------------------------------
c     convective and diffusive fluxes
      common/comfl1/ce(nvar,nn)
      common/comfl2/iflux,nordre
c     -----------------------------------------------------------------
c     nodal gradients
      common/comgr0/beta,beta1
      common/comgr1/dx(4,nn),dy(4,nn)
c     -----------------------------------------------------------------
c     free stream solution
      common/cominf/roin,uxin,uyin,pin,roout,uxout,uyout,pout,
     &              ein,ruxin,ruyin,eout,ruxout,ruyout
c     initial solution strategy
      common/comcnt/ncont
c     -----------------------------------------------------------------
c     real gaz physical quantities
      common/comgam/gam,gam1,gam4
c     viscous calculation flag
      common/comvi1/ivis
c     body temperature, viscosity, prandtl and reynols numbers 
      common/comvi2/tbrd2,pr,reynolds,tinf,tinfd,xmach,iecc,froud
c     -----------------------------------------------------------------
c     frontiere
c     -----------------------------------------------------------------
      common/frontier/nslog2,log2(nnfr),nslog3,log3(nnfr),
     1                nslog4,log4(nnfr),nslog5,log5(nnfr),
     1                nslog6,log6(nnfr)
c     -----------------------------------------------------------------
      common/comaxi/iaxi,rrs(nn),airsa(nn),airta(nnt)
c     -----------------------------------------------------------------
c     turbulence 
c     -----------------------------------------------------------------
      common/comvi3/reylam(nn),dist(nn),reyturb(nn),reyt2(nn)
      common/turb1/iturb,delta,ilaw,ncont1,prt,xkin,xein 
      common/turb2/numc(nn),nbt(nnfr),jtb(maxt,nnfr),nbp(nnfr),
     1             ipb(maxp,nnfr),nswall(nn)
      common/turb3/xtmin,xtmax,ytmin,ytmax
      common/turb4/fluro(nnsg),flurob(nn)
c     -----------------------------------------------------------------
