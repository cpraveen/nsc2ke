      subroutine gravity
c 
cc gravity forces
c bijan mohammadi, INRIA
c
      include 'param2d.h'
c 
      do 100 is=1,ns
      air=airsa(is)
      ce(3,is)=ce(3,is)-air*froud*ua(1,is)
100   continue
c
          return
          end
