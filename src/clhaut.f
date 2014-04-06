      subroutine clhaut
c     -----------------------------------------------------------------
c     computation of the control volume altitudes
c     -----------------------------------------------------------------
      include 'param2d.h'
c     -----------------------------------------------------------------
c     local variables definition
      integer nubo1  , nubo2  , nubo3  , is , jt
      real    dbxx(3), dbyy(3), haut(3), htt, ait, 
     &        x1, x2 , x3, y1 , y2, y3
c
c     initializing the vertex altitudes
c
      do 77 is=1,ns
         dthaut(is)     = 1.0e+12
77    continue
c
      do 1000 jt=1,nt
c
         nubo1          = nu(1,jt)
         nubo2          = nu(2,jt)
         nubo3          = nu(3,jt)
         x1             = coor(1,nubo1)
         y1             = coor(2,nubo1)
         x2             = coor(1,nubo2)
         y2             = coor(2,nubo2)
         x3             = coor(1,nubo3)
         y3             = coor(2,nubo3)
         dbxx(1)        = y2 - y3
         dbxx(2)        = y3 - y1
         dbxx(3)        = y1 - y2
         dbyy(1)        = x3 - x2
         dbyy(2)        = x1 - x3
         dbyy(3)        = x2 - x1
c
         ait            = 1.0/(2.0*airt(jt))
         haut(1)        = ait*(sqrt(dbxx(1)*dbxx(1) + dbyy(1)*dbyy(1)))
         haut(2)        = ait*(sqrt(dbxx(2)*dbxx(2) + dbyy(2)*dbyy(2)))
         haut(3)        = ait*(sqrt(dbxx(3)*dbxx(3) + dbyy(3)*dbyy(3)))
c
         htt            = min(1.0/haut(1), 1.0/haut(2))
         htt            = min(htt        , 1.0/haut(3))

         dthaut(nubo1)  = min(dthaut(nubo1), htt)
         dthaut(nubo2)  = min(dthaut(nubo2), htt)
         dthaut(nubo3)  = min(dthaut(nubo3), htt)
c
1000  continue
c
      return
      end
