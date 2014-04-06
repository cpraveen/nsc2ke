      subroutine sa
c 
cc axisymmetric corrections
c
      include 'param2d.h'
c
       dtier=2./3.
c
      do is=1,ns
c
      airr=airsa(is)
      air=airs(is)
      r=rrs(is)
      reyt=reylam(is)+reyturb(is)
      p=pres(is)
      u=ua(2,is)/ua(1,is)
      v=ua(3,is)/ua(1,is)
      ux=dx(2,is)
      vx=dx(3,is)
      uy=dy(2,is)
      vy=dy(3,is)
      divu=ux+vy+v/r
c
      ce(2,is)=ce(2,is)-air*ivis*dtier*reyt*vx
c
      xcof=p-ivis*dtier*reyt*(vy-v/r)
      tautt=reyt*(2*v/r-dtier*divu)
      ce(3,is)=ce(3,is)+air*(xcof-tautt)
c
      xcof1=u*vx+v*(ux+2*vy)-v**2/r
      ce(4,is)=ce(4,is)-air*ivis*dtier*reyt*xcof1
c
      if(coor(2,is).lt.1.e-6.and.logfr(is).ne.3) then
      ce(3,is)=0.0
      endif
c
      enddo
c
          return
          end
