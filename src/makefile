F77=gfortran
LINK=gfortran -O
FLAGS=-c -O

OBJS=nsc2ke.o config.o mailla.o aires.o seg2d.o\
clhaut.o init.o calprc.o result.o fluosh.o flucin.o\
viscdt.o fluroe.o cdl.o vitfrot.o geowall.o ke_two.o\
sa.o axigeo.o ke_law.o resexp.o caldtl.o keps2d.o\
loglaw.o isovat.o gravity.o

.f.o:
	$(F77) $(FLAGS) $*.f

NSC2KE:$(OBJS)
	$(LINK) $(OBJS) -o NSC2KE

clean:;    rm -f $(OBJS) core
