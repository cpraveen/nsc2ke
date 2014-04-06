      subroutine ke_law
      include 'param2d.h'
c
      cmu=0.09
      c1=0.1296
      c2=1.83333
      xmact0=0.25
c
      do 200 is=1,ns       
      x1=coor(1,is)
      x2=coor(2,is)
      reyturb(is)=0.0
      if(x1.lt.xtmin.or.x1.gt.xtmax.or.
     1   x2.lt.ytmin.or.x2.gt.ytmax) goto 200
      ro=ua(1,is)
      xk=ua(5,is)
      xes=ua(6,is)
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc compressibility corrections
cccccccccccccccccccccccccccccccccccccccccccccccccccc
      xmactu=sqrt( 2.*ro*xk / (gam*pres(is)) )
      xmactu2=xmactu**2
      xtest=xmactu-xmact0
      heavyside=1.
      if(xtest.le.0.) heavyside=0.
      fcomp=heavyside*(1.-exp(-(xtest/0.66)**2))
      xed=xes*(1.+fcomp)
cccccccccccccccccccccccccccccccccccccccccccccccccccc
      reyturb(is)=cmu*ro*xk*xk/xes
      enorm=(dx(3,is)+dy(2,is))**2
c
      sk=enorm*reyturb(is)-ro*xed
      se=c1*enorm*ro*xk-c2*ro*xed**2/xk
c
      ce(5,is)=ce(5,is)+airsa(is)*sk
      ce(6,is)=ce(6,is)+airsa(is)*se
c
200   continue
c
      return
      end
