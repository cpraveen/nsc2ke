      SUBROUTINE calprc
      include 'param2d.h'
c
       XRHO = 1.E-6
       DO 100 is=1,ns
       if(ua(1,is).lt.xrho) then
              print *,logfr(is),coor(1,is),coor(2,is),pres(is)
              stop 'density < 0 , calprc '
       endif
       pres(is)=gam1*(ua(4,is) - 0.5*(ua(2,is)*ua(2,is)+
     &               ua(3,is)*ua(3,is))/ua(1,is))
        IF(pres(is).lt.xrho) then
              print *,logfr(is),coor(1,is),coor(2,is),pres(is)
              stop 'pressure < 0 calprc '
       endif
100    CONTINUE
C
      RETURN
      END
