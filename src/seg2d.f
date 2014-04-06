      subroutine seg2d
c     -----------------------------------------------------------------
c     construction of the edges of the triangulation
c     -----------------------------------------------------------------
      include 'param2d.h'
c     -----------------------------------------------------------------
c     local variables definition 
      parameter (nvmax = 15)
      integer is    , jt  , kv    , k 
      integer iseg
      integer is1   , is2   
      integer nuor  , nuex, nub1
      integer jaret(nn,nvmax,2)   , nor(3) , nex(3)
      real    x1,x2 , x3  , y1, y2, y3
      real    esp   , cap
      real    vx    , vy
      real    dbxx(3)    , dbyy(3)
      save
c
      nex(1)                 = 2
      nex(2)                 = 3
      nex(3)                 = 1
      nor(1)                 = 1
      nor(2)                 = 2
      nor(3)                 = 3
c
      do 2 iseg=1,nnsg
         vnocl(1,iseg)       = 0.0
         vnocl(2,iseg)       = 0.0    
         vnocl(3,iseg)       = 0.0    
2     continue
      do 10 is=1,ns
         do 5 kv=1,nvmax
            jaret(is,kv,1)   = 0
            jaret(is,kv,2)   = 0
5        continue
10    continue
c
      esp                    = 1.0/6.0
      nseg                   = 0
c     
      do 1000 jt=1,nt
         is1                 = nu(1,jt)
         is2                 = nu(2,jt)
         is3                 = nu(3,jt)
         x1                  = coor(1,is1)
         y1                  = coor(2,is1)
         x2                  = coor(1,is2)
         y2                  = coor(2,is2)
         x3                  = coor(1,is3)
         y3                  = coor(2,is3)
         dbxx(1)             = y2 - y3
         dbxx(2)             = y3 - y1
         dbxx(3)             = y1 - y2
         dbyy(1)             = x3 - x2
         dbyy(2)             = x1 - x3
         dbyy(3)             = x2 - x1                  
c
         do 500 k=1,3
            is1              = nu(nor(k),jt)
            is2              = nu(nex(k),jt)
            do 100 kv=1,nvmax
               if (jaret(is1,kv,1) .eq. 0) goto 110
               if (jaret(is1,kv,1) .eq. is2) then
                  ig0        = jaret(is1,kv,2)
                  nub1       = iabs(nubo(1,ig0))
                  nubo(1,ig0)=-nubo(1,ig0)
                  nub2       = nubo(2,ig0)
                  if (is1 .eq. nub1 .and. is2 .eq. nub2) then
                     esp     = abs(esp)
                  else
                     if (is1 .eq. nub2 .and. is2 .eq. nub1) then
                        esp  =-abs(esp)
                     else
                        print*,'error in seg2d ',is1,is2,nub1,nub2
                     endif
                  endif
                  goto 200
               endif
100         continue
c
110         do 120 kv=1,nvmax
               if (jaret(is2,kv,1) .eq. 0) goto 130
               if (jaret(is2,kv,1) .eq. is1) then
                  ig0        = jaret(is2,kv,2)
                  nub1       = nubo(1,ig0)
                  nubo(1,ig0)=-nubo(1,ig0)
                  nub2       = nubo(2,ig0)
                  if (is1 .eq. nub1 .and. is2 .eq. nub2) then
                     esp     = abs(esp)
                  else
                     if (is1 .eq. nub2 .and. is2 .eq. nub1) then
                        esp  =-abs(esp)
                     else
                        print*,'error in seg2d ',is1,is2,nub1,nub2
                     endif
                  endif
                  goto 200
               endif
120         continue
c
            kv               = kv + 1
            if (kv .gt. nvmax) then
               print*,'increase nvmax : ', nvmax
               stop 'seg2d'
            endif
130         nseg  = nseg + 1
            if (nseg .gt. nnsg) then
               print*,'increase nnsg  : ', nnsg
               stop 'seg2d'
            endif
c
            jaret(is2,kv,1)  = is1
            jaret(is2,kv,2)  = nseg
            nubo(1,nseg)     =-is1
            nubo(2,nseg)     = is2
            ig0              = nseg
            esp              = abs(esp)
200         continue
c
c           computation of the control volume boundary normals
c
            vnocl(1,ig0)     = vnocl(1,ig0) + 
     &                         esp*(dbxx(nor(k)) - dbxx(nex(k)))
            vnocl(2,ig0)     = vnocl(2,ig0) + 
     &                         esp*(dbyy(nor(k)) - dbyy(nex(k)))
500      continue
1000  continue
c
      do 3000 iseg=1,nseg
         if (abs(vnocl(1,iseg)) .lt. 1.0e-16 .and.
     &       abs(vnocl(2,iseg)) .lt .1.0e-16) then
             vnocl(1,iseg)   = 1.0e-8
             vnocl(2,iseg)   = 1.0e-8
         endif
3000  continue
c
c     normalization of the control volume boundary normals
c     construction of the boundary edges
c
      nfr                    = 0
      do 2000 iseg=1,nseg
         cap                 = sqrt(vnocl(1,iseg)*vnocl(1,iseg) +
     &                              vnocl(2,iseg)*vnocl(2,iseg))
         vnocl(1,iseg)       = vnocl(1,iseg)/cap
         vnocl(2,iseg)       = vnocl(2,iseg)/cap
         vnocl(3,iseg)       = cap
         nub1                = nubo(1,iseg)
         if (nub1 .lt. 0) then
c
c           this is a boundary edge
c
            nfr              = nfr + 1
            nufr(nfr)        = iseg
            nubo(1,iseg)     =-nubo(1,iseg)
         endif
2000  continue
c
      print*,'----------------------------------------------'
      print*,'total number of edges    : ', nseg
      print*,'number of boundary edges : ', nfr
c
c     computation of the boundary vertices normals
c          normales sortantes (tauy, -taux)
c
      do 5000 is=1,ns
         vnox(is)            = 0.0
         vnoy(is)            = 0.0
5000   continue
      do 6000 ifr=1,nfr
         iseg                = nufr(ifr)
         nuor                = nubo(1,iseg)
         nuex                = nubo(2,iseg)
         vx                  = coor(2,nuex) - coor(2,nuor)
         vy                  = coor(1,nuor) - coor(1,nuex)  
            vnox(nuor)       = vnox(nuor) + vx * 0.5
            vnoy(nuor)       = vnoy(nuor) + vy * 0.5
            vnox(nuex)       = vnox(nuex) + vx * 0.5
            vnoy(nuex)       = vnoy(nuex) + vy * 0.5     
6000  continue
c
      do 5550 is=1,ns
         if(logfr(is).ne.0) then
         vno(is)   = sqrt(vnox(is)**2+vnoy(is)**2)
         vnox(is)  = vnox(is)/vno(is)
         vnoy(is)  = vnoy(is)/vno(is)
                            else
         vno(is)   = 0.0
         vnox(is)  = 0.0
         vnoy(is)  = 0.0
         endif
5550  continue
c
      return
      end
