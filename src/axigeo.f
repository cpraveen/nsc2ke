      subroutine axigeo
      include 'param2d.h'
      real vnt(nn)
       
      if(iaxi.eq.0) then
       do is=1,ns
        airsa(is)=airs(is)
        rrs(is)=1.0
       enddo
       do it=1,nt
        airta(it)=1.0/airt(it)
       enddo
       return
      endif
c
       do it=1,nt
         is1=nu(1,it)
         is2=nu(2,it)
         is3=nu(3,it)
         y1=coor(2,is1)
         y2=coor(2,is2)
         y3=coor(2,is3)
         rrt=(y1+y2+y3)/3.0
         rrs(is1)=rrs(is1)+airt(it)*rrt
         rrs(is2)=rrs(is2)+airt(it)*rrt
         rrs(is3)=rrs(is3)+airt(it)*rrt
         airta(it)=rrt/airt(it)
       enddo  
c
       do  is=1,ns
         if(coor(2,is).lt.1.e-6) then
           rrs(is)=rrs(is)/(3.*airs(is))
                                else
           rrs(is)=coor(2,is)
         endif
        airsa(is)=airs(is)*rrs(is)
        vnt(is)=0.
       enddo
c
c redefinition de vno pour l'axi
c
      do ifr=1,nfr
         iseg          = nufr(ifr)
         nuor          = nubo(1,iseg)
         nuex          = nubo(2,iseg)
         r             =(coor(2,nuor)+coor(2,nuex))*0.5
         vnt(nuex)=vnt(nuex)+0.5*r
         vnt(nuor)=vnt(nuor)+0.5*r
      enddo
c
      do is=1,ns
         vno(is)   = vno(is)*vnt(is)
      enddo
c
       do iseg=1,nseg             
          is1=nubo(1,iseg)
          is2=nubo(2,iseg)   
          r=(coor(2,is1)+coor(2,is2))*0.5
          vnocl(3,iseg)=vnocl(3,iseg)*r       
       enddo
c
         return
         end
