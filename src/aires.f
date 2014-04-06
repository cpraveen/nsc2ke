      subroutine aires
c     -----------------------------------------------------------------
c     computation of the triangle and control volume areas
c     -----------------------------------------------------------------    
      include 'param2d.h'
c     -----------------------------------------------------------------
      real    a(2), b(2), c(2)
c
      do 5 is=1,ns
         airs(is)  = 0.0
5     continue
c
      do 10 jt=1,nt
c
         is1       = nu(1,jt)
         is2       = nu(2,jt)
         is3       = nu(3,jt)
c
c        computing the current triangle area
c
         a(1)      = coor(1,is1)
         a(2)      = coor(2,is1)
         b(1)      = coor(1,is2)
         b(2)      = coor(2,is2)
         c(1)      = coor(1,is3)
         c(2)      = coor(2,is3)
         u         = a(2)-c(2)
         v         = c(1)-a(1)
         w         = a(1)*c(2)-a(2)*c(1)
         vol       = u*b(1)+v*b(2)+w
         airt(jt)  = 0.5*abs(vol)
c
c        gathering the current triangle area into the 
c        control volume one 
c
         airs(is1) = airs(is1)+airt(jt)/3.0
         airs(is2) = airs(is2)+airt(jt)/3.0
         airs(is3) = airs(is3)+airt(jt)/3.0
c
10    continue
c
        som_airt=0.
        do it=1,nt
         som_airt=som_airt+airt(it)
        enddo
c
        som_airs=0.
        do is=1,ns
         som_airs=som_airs+airs(is)
        enddo
        print *,'test airt , som airt =',som_airt
        print *,'test airs , som airs =',som_airs
c
      return
      end
