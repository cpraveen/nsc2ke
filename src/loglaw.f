      subroutine loglaw
      include 'param2d.h'
c
c classical wall laws technique
c        weak form
c
      rklaw  =0.419
      claw   =5.445
      nitmax =500
      eps1   =1.e-5           
      cmu    =0.09
c
      do 1000 isc=1,nslog3
c
        is=log3(isc)
        ro=ua(1,is)
        temp=pres(is)/(gam1*ro)
        rey=reylam(is)
        reytot=reylam(is)+reyturb(is)
c     
c on norme le vecteur normal entrant
c
        rnorm=vno(is)
        xn1=-vnox(is)
        xn2=-vnoy(is)
        xt1= xn2
        xt2=-xn1
c
c calcul de la vitesse tangentielle
c
      u=ua(2,is)/ro
      v=ua(3,is)/ro
      utang=xt1*u+xt2*v
c
c methode de newton
c
        uf=sqrt(rey*abs(utang)/delta)
        yplus=ro*delta*uf/rey
c
        if(yplus.lt.5.0) then
           uf=rey*5.0/delta/ro
        endif
c
        if(yplus.gt.11.6) then
c     
           nit=0
500     uplus=log(yplus)/rklaw+claw
        uplusp=uplus+1./rklaw
        uf1=uf+(abs(utang)-uf*uplus)/uplusp
           resuf=abs((uf1-uf)/uf)
           uf=uf1
           if(resuf.gt.eps1.and.nit.lt.nitmax) then
             yplus=ro*delta*abs(uf)/rey
             nit=nit+1
             go to 500
           endif
           if(nit.eq.nitmax) then
           print*,'nit petit,is=',is
           stop
           endif          
           yplus=ro*delta*uf/rey
c
        endif
c
        if(iturb.ne.0) then
        ce(5,is)=0.0
        ce(6,is)=0.0
        un(5,is)=ro*uf*uf/0.3
        un(6,is)=ro*abs(uf*uf*uf)/(rklaw*delta)
        reyturb(is)=cmu*un(5,is)**2/un(6,is)
        endif
cc        print *,yplus,ua(1,is)*uf**2
c
        utgt=1.0
        if(nbp(isc).gt.0) then
        utgt=0.0
        do k=1,nbp(isc)
        isb=ipb(k,isc)
        utgt=utgt+(xt1*ua(2,isb)+xt2*ua(3,isb))/ua(1,isb)
        enddo                            
        utgt=utgt/float(nbp(isc))
        utgt=sign(1.,utgt)
        endif

      conpa=utgt*ro*uf*uf
      ce(2,is)=ce(2,is)-(xt1*conpa-xn1*pres(is))*rnorm
      ce(3,is)=ce(3,is)-(xt2*conpa-xn2*pres(is))*rnorm
      ce(4,is)=ce(4,is)-abs(conpa*utang*rnorm)
1000  continue
      return
      end
