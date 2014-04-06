       subroutine geowall
       include 'param2d.h'
c
c maxt  : nombre maximal de tetraedres contenant une normale
c maxp  : nombre maximal de points pour le calcul de la contrainte
c         parietale
c nbt : nombre de tetraedre interceptant la normale en if11
c jtb : numero du tetraedre solution
c nbp : nombre de points pour le calcul de la contrainte en if11
c ipb : numero des points solution
c
      integer deja(nn)
      integer inb,is,is1,is2,jt,ip1(3),ip2(3),nb,isb
      integer k,if11,j
c   
      data ip1 /2,3,1/
      data ip2 /3,1,2/
c      
      do if11=1,nslog3
         nbt(if11)=0
         do inb=1,maxt
            jtb(inb,if11)=-1
         enddo
      enddo
c
      do 100 jt=1,nt
c
         do 50 k=1,3    
c
            is=nu(k,jt)
            is1=nu(ip1(k),jt)
            is2=nu(ip2(k),jt)
c 
c si is est un point adherent
c ...........................
c
            if(numc(is).gt.0) then    
c
            if11=numc(is)
c
               if(nbt(if11).lt.maxt) then
c
c un nouveau tetraedre est solution
c
                  nbt(if11)=nbt(if11)+1
                  jtb(nbt(if11),if11)=jt
c
               endif
         endif
c
50    continue     
100   continue
c 
c *** points pour le calcul de la vitesse de frottement *****
c ...........................................................
c
      do 110 is=1,ns
         deja(is)=-1
110   continue
c
      do 120 if11=1,nslog3
         nbp(if11)=0
120   continue
c
c boucle sur les points adherents      
c
      do 500 if11=1,nslog3
c
         do 400 nb=1,nbt(if11)
c
            do 300 k=1,3
               isb=nu(k,jtb(nb,if11))
c
               if(numc(isb).eq.0.and.deja(isb).lt.0) then
                   nbp(if11)=nbp(if11)+1
                   ipb(nbp(if11),if11)=isb
               endif
c
               if(nbt(if11).gt.1) deja(isb)=1
c
               if(nbt(if11).eq.nb) deja(isb)=-1
c
300         continue
400      continue
c
      do 450 is=1,ns
         deja(is)=-1
450   continue

500   continue
c
       do is=1,ns
       dist(is)=1.e10
       x1=coor(1,is)
       y1=coor(2,is)
       do j=1,nslog3
       js=log3(j)
       x2=coor(1,js)
       y2=coor(2,js)
       dd=sqrt((x1-x2)**2+(y1-y2)**2)
       if(dd.lt.dist(is)) then
        dist(is)=dd
        nswall(is)=js
       endif
       enddo
       enddo
c
      return
      end
