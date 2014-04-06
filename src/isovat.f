      SUBROUTINE isovat(ifile,F,COOR,NVAL,VAL,x0,y0,x1,y1)
C
C       TRACE DES LIGNES ISOVALEURS D'UNE FONCTION LINEAIRE
C         SUR UN TRIANGLE
C       IOPFEN=1 RESTRICTION A UNE FENETRE
C       IOPFEN=0 PAS DE RESTRICTION
C       IOPFEN=-1 RESTRICTION A L'EXTERIEUR D'UNE FENETRE
C  *************************************************************
C
      DIMENSION F(3),COOR(2,3),VAL(100)
      DIMENSION IP1(3)
C
      epsi=1.e-5
      IP1(1)=2
      IP1(2)=3
      IP1(3)=1
      FF1=F(1)
      FF2=F(2)
      FF3=F(3)
      FFMA=MAX(ABS(FF1),ABS(FF2))
      FFMA=MAX(ffma,ABS(FF3))
      D12=ABS(FF1-FF2)
      D23=ABS(FF2-FF3)
      IF(D12+D23.LT.AMAX1(epsi,epsi*FFMA)) GOTO 1000
C  PAS DE RESTRICTION
C  ******************    
      DO 100 IVAL=1,NVAL
      VAL1=VAL(IVAL)
      ITR=0
      DO 110 K=1,3
      FK=F(K)
      FK1=F(IP1(K))
      FMI=MIN(FK,FK1)
      FMA=MAX(FK,FK1)
      DIF=FMA-FMI
      IF(DIF.LT.epsi) GOTO 110
      EPS=epsi*DIF
      IF(VAL1.LT.FMI-EPS.OR.VAL1.GT.FMA+EPS) GOTO 110
      HH=ABS(FK-VAL1)/DIF
      X=COOR(1,K)+HH*(COOR(1,IP1(K))-COOR(1,K))
      Y=COOR(2,K)+HH*(COOR(2,IP1(K))-COOR(2,K))
      IF(ITR.EQ.0) GOTO 115
      write(ifile,*) x,y
      write(ifile,*) 
      GOTO 100
115   ITR=1
      write(ifile,*) x,y
110   CONTINUE
100   CONTINUE
1000   return      
      END
