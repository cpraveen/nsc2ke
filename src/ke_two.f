      subroutine ke_two
      include 'param2d.h'
c
      cmu=0.09
      c1=0.1296
      c2=1.83333
      xkap=0.41
      qs3=4./3.
      ds3=2./3.
      xmact0=0.25
      cl=xkap*cmu**(-3./4.)
      if(xmach.lt.1.) xclmu=0.0142
      if(xmach.ge.1.) xclmu=0.0085
c
      do 300 is=1,ns
      ro=ua(1,is)
      y=dist(is)
      iwall=nswall(is)
      x1=coor(1,is)
      x2=coor(2,is)
      xk=ua(5,is)
      xes=ua(6,is)
      enorm=(dx(3,is)+dy(2,is))**2
      divs=dx(2,is)+dy(3,is)
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc compressibility corrections
cccccccccccccccccccccccccccccccccccccccccccccccccccc
      xmactu=sqrt( 2.*ro*xk / (gam*pres(is)) )
      xmactu2=xmactu**2
      xtest=xmactu-xmact0
      heavyside=1.
      if(xtest.le.0.) heavyside=0.
      fzeman=heavyside*(1.-exp(-(xtest/0.66)**2))
      fcomp=fzeman
      xed=xes*(1.+fcomp)
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
      xmu=reylam(iwall)
      yplus=sqrt(xk*ua(1,iwall)*ro)*y/xmu
c
      if(yplus.gt.200..or.y.gt.delta) then
      reyturb(is)=0.0
ccc high reynolds region
      if(x1.lt.xtmin.or.x1.gt.xtmax.or.
     1   x2.lt.ytmin.or.x2.gt.ytmax) goto 300
c
      reyturb(is)=cmu*ro*xk*xk/xes
      sk=enorm*reyturb(is)-ro*xed 
      se=c1*enorm*ro*xk-c2*ro*xed*xed/xk  
      else
ccc near wall region
      reyturb(is)=0.0
      if(logfr(is).ne.0.or.x1.lt.xtmin.or.x1.gt.xtmax.or.
     1   x2.lt.ytmin.or.x2.gt.ytmax) goto 300
      xlmu=cl*y*(1.-exp(-yplus*xclmu))
      xlep=cl*y*(1.-exp(-yplus/(2.0*cl)))
      reyturb(is)=cmu*xlmu*ro*sqrt(xk) 
      un(6,is)=ro*xk**1.5/xlep
      sk=enorm*reyturb(is)-ro*xed
      se=0.0
      ce(6,is)=0.0 
      endif
c
      ce(5,is)=ce(5,is)+sk*airsa(is)
      ce(6,is)=ce(6,is)+se*airsa(is)
c
300   continue
c
      return
      end
