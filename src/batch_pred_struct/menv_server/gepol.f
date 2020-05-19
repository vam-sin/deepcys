C
C***********************************************************************************
C	this routine calculated the solvent exposed fraction of each atom, averaged
C	over groups and stores the values in the array DGS
C
C	The routine is based on GEPOL93 written by 
C*** Written by
C
C       J.L. Pascual-Ahuir, E. Silla and I. Tunon
C       Departamento de Quimica Fisica.
C       Facultad de Quimica.
C       Universidad de Valencia.
C       C/Dr.Moliner 50.
C       Burjassot (Valencia) 46100
C       SPAIN
C
C	see
C    D- J. L. Pascual-Ahuir, I Tunon and E. Silla
C       GEPOL: An improved description of molecular surfaces. III. A New
C       algorithm for the computation of the Solvent-Excluding Surface.
C       (explain the algorithm used in GEPOL93)
C       J. Comput. Chem., 15:1127-1138, 1994.
C
C	and earlier references.
C*********************************************************************************
C
      SUBROUTINE GEPOL(NATOM,IOUT,XE,YE,ZE,RE,NG,DGS,LPR,RD)
      IMPLICIT NONE
      INTEGER I,IOUT,IJ,K,IM1,II
      INTEGER*2 NG(*)
      
      LOGICAL GHOST
      LOGICAL LPR
      LOGICAL REDU

      INTEGER*2 IUSE

      INTEGER*4 ISA,ISO,ITO
      INTEGER*4 LT
      INTEGER*4 MC,MV
      INTEGER*4 NATOM,NCOR,NDIV,NP
       
      REAL*8 DGS(*)
      REAL*8 AP
      REAL*8 DVEC
      REAL*8 OFAC
      REAL*8 RD,RE(*),RMIN
      REAL*8 XC1,XE(*),XP
      REAL*8 YC1,YE(*),YP
      REAL*8 ZC1,ZE(*),ZP

      REAL*8 CV
      REAL*8 STOT,VOL,AVE

      CHARACTER*5  KSURF
      
      PARAMETER (MC=1000000,MV=1000000)
      
      COMMON/POLI/CV(32,3),XC1(15360),YC1(15360),ZC1(15360)
      COMMON/CSFE/IUSE(MC)
      COMMON/PUN/ITO(MV),ISO(MV),ISA(MV),XP(MV),YP(MV),ZP(MV),AP(MV)
    
C*****Read input****

C Set defaults and Read General Input
       write(*,*)"natom in gepol",natom
       if (NATOM.LT.3000) THEN
       DO I=1,NATOM
c      WRITE (IOUT,1000) I,XE(I),YE(I),ZE(I),RE(I),DGS(I),NG(I)
       WRITE (71,1000) I,XE(I),YE(I),ZE(I),RE(I),DGS(I),NG(I)
       if (NATOM.GE.5000) return
       ENDDO
       endif
      CALL READIN (KSURF,REDU,LT,RMIN,OFAC,RD,NDIV,DVEC,LPR)

C Read COOR File
      CALL READCOOR(NATOM,GHOST,LPR,RE)

      NCOR=NATOM

C*****Make tesselation for sphere of radius 1.0*********

      CALL TES(LPR)

      CALL DIVIDE(NDIV,LPR)

C*****Beginning  the calculations ******

C
C  Compute the surface

      CALL SUME(NCOR,RD,'SUMA',LPR,RE)
C     stop 1235

      CALL GEOCAV(NCOR,NP,NDIV,GHOST,LPR,XE,YE,ZE,RE)


C
C Compute the area and volume
      CALL VOLARE(NP,LPR,GHOST,STOT,VOL,XE,YE,ZE,RE)
      write(*,*) "NP from volare",np
      CALL PRISPH(NATOM,NCOR,NP,REDU,LPR,DGS)
       if (NATOM.LT.3000) THEN
       DO I=1,NATOM
       if (I.eq.1) then
         write(IOUT,*)'Before averaging'
         write(IOUT,*)'    I         X         Y         Z        R',
     *'    DGS       NG   '
       endif
c        WRITE (IOUT,1000) I,XE(I),YE(I),ZE(I),RE(I),DGS(I),NG(I)
         WRITE (72,1000) I,XE(I),YE(I),ZE(I),RE(I),DGS(I),NG(I)
       if (NATOM.GE.100) return
       ENDDO
       endif
CC
C	calculate group averages
C
      AVE=DGS(1)
      IJ=1
      DO I=2,NATOM! reinitated on 25/11/15
c     DO I=1,NATOM ! changed on 21st Augu 2015| again changed on 25/11/15
      IF (NG(I).EQ.0.AND.I.NE.NATOM) THEN
        AVE=AVE+DGS(I)
      ELSE IF (NG(I).EQ.1.OR.I.EQ.NATOM) THEN
        II=I
        IF (I.EQ.NATOM) THEN
          II=I+1
          AVE=AVE+DGS(I)
        ENDIF
        AVE=AVE/(II-IJ)
        IM1=II-1
        IF (IM1.GT.0) THEN
c         write(*,*)"within gepol, IJ,IM1,II", IJ,IM1,II  !5th Aug 2015
          DO K=IJ,IM1 !changed back on 21st Aug 2015
c         DO K=IJ,IM1+1 !changed on 5th Aug 2015
           DGS(K)=AVE
          ENDDO
        ELSE IF (IM1.EQ.0) THEN
          DGS(IJ)=AVE
        ENDIF
        IJ=I ! changed back on 21st Aug 2015
c       IJ=I+1!changed on 5th Aug 2015
        AVE=DGS(I)
      ENDIF
      ENDDO

      IF (LPR) THEN
       if (NATOM.LT.3000) THEN
       DO I=1,NATOM
       if (I.eq.1) then
         write(IOUT,*)'    I         X         Y         Z        R',
     *'    DGS       NG   '
       endif
c        WRITE (IOUT,1000) I,XE(I),YE(I),ZE(I),RE(I),DGS(I),NG(I)
         WRITE (73,1000) I,XE(I),YE(I),ZE(I),RE(I),DGS(I),NG(I)
       if (NATOM.GE.5000) return
       ENDDO
       endif
       ENDIF
       if (NATOM.LT.3000) THEN
       DO I=1,NATOM
       if (I.eq.1) then
         write(IOUT,*)'after averaging'
         write(IOUT,*)'    I         X         Y         Z        R',
     *'    DGS       NG   '
       endif
         WRITE (IOUT,1000) I,XE(I),YE(I),ZE(I),RE(I),DGS(I),NG(I)
       if (NATOM.GE.100) return
       ENDDO
       endif
 1000 FORMAT(I10,5F10.4,I5)
      RETURN
      END

C
C

      SUBROUTINE READIN (KSURF,REDU,LT,RMIN,OFAC,RD,NDIV,DVEC,LPR)
C     ------------------------------------------------------------------
C     This subroutine reads the general input
C     ------------------------------------------------------------------
      IMPLICIT NONE
      
      LOGICAL ERROR
      LOGICAL LPR
      LOGICAL REDU

      INTEGER*4 LN,LT
      INTEGER*4 NDIV
  
      REAL*8 DVEC
      REAL*8 OFAC
      REAL*8 RD,RMIN

      CHARACTER*5 KEY,KSURF
      
C Set Defaults

      KSURF='ASURF'
      REDU=.FALSE.
C     LPR=.TRUE.
C     RD=1.0D0
      OFAC=0.80
      RMIN=0.50D0
      NDIV=3
      DVEC=1.0D0
      RETURN
      END
C

      SUBROUTINE READCOOR(NATOM,GHOST,LPR,RE)
C     ---------------------------------------------------------------
C     This reads coordinates and radii
C     ---------------------------------------------------------------
      IMPLICIT NONE

      LOGICAL GHOST
      LOGICAL LPR

      INTEGER*2 IUSE

      INTEGER*4 I
      INTEGER*4 MC
      INTEGER*4 NAT0,NATOM,NESFI,NESFP

      REAL*8 RE(*)

      PARAMETER (MC=1000000)
      
      COMMON/CSFE/IUSE(MC)

      NAT0=0
      NESFP=0
      NESFI=0
      GHOST=.FALSE.

      DO  I = 1,NATOM
       IF((RE(I).GT.-0.00001).AND.(RE(I).LT.0.00001)) THEN

          RE(I)=0.0000000D0
          IUSE(I)=1
          NAT0=NAT0+1

       ELSE IF(RE(I).GT.0.00001) THEN

          IUSE(I)=6
          NESFI=NESFI+1        
  
       ELSE IF(RE(I).LT.-0.00001) THEN
       write(*,*)'GHOST!!! ',I

          GHOST=.TRUE.
          IUSE(I)=3
          RE(I)=ABS(RE(I))
          NESFP=NESFP+1

      END IF

      END DO


      RETURN
      END
C


      BLOCK DATA
C     -----------------------------------------------------------------
C     This has the information about the vertices
C     -----------------------------------------------------------------
      IMPLICIT NONE
      INTEGER*4 JVT1,JVT2
      COMMON/PENTA/JVT1(3,60),JVT2(3,4)
      DATA JVT1 /  1,   6,   2,  1,   2,   3,  1,   3,   4,
     &             1,   4,   5,  1,   5,   6,  7,   2,   6,
     &             8,   3,   2,  9,   4,   3, 10,   5,   4,
     &            11,   6,   5,  8,   2,  12,  9,   3,  13,
     &            10,   4,  14, 11,   5,  15,  7,   6,  16,
     &             7,  12,   2,  8,  13,   3,  9,  14,   4,
     &            10,  15,   5, 11,  16,   6,  8,  12,  18,
     &             9,  13,  19, 10,  14,  20, 11,  15,  21,
     &             7,  16,  17,  7,  17,  12,  8,  18,  13,
     &             9,  19,  14, 10,  20,  15, 11,  21,  16,
     &            22,  12,  17, 23,  13,  18, 24,  14,  19,
     &            25,  15,  20, 26,  16,  21, 22,  18,  12,
     &            23,  19,  13, 24,  20,  14, 25,  21,  15,
     &            26,  17,  16, 22,  17,  27, 23,  18,  28,
     &            24,  19,  29, 25,  20,  30, 26,  21,  31,
     &            22,  28,  18, 23,  29,  19, 24,  30,  20,
     &            25,  31,  21, 26,  27,  17, 22,  27,  28,
     &            23,  28,  29, 24,  29,  30, 25,  30,  31,
     &            26,  31,  27, 32,  28,  27, 32,  29,  28,
     &            32,  30,  29, 32,  31,  30, 32,  27,  31 /
      DATA JVT2 /  1,   5,   4,
     &             5,   2,   6,
     &             4,   6,   3,
     &             6,   4,   5 /
      END
C 
      SUBROUTINE TES(LPR)
C     --------------------------------------------------------------------
C     This computes the triangle vertex coordinates for a sphere of radius
C     one, projecting the pentakisdodecahedro onto it.
C     --------------------------------------------------------------------
      IMPLICIT NONE
      
      LOGICAL LPR

      INTEGER*4 I,II
      INTEGER*4 J

      REAL*8 XC1,YC1,ZC1

      REAL*8 CTH,CV
      REAL*8 FI,FIR,FIV
      REAL*8 STH
      REAL*8 TH,THEV
      
      COMMON/POLI/CV(32,3),XC1(15360),YC1(15360),ZC1(15360)
      
      DIMENSION THEV(6),FIV(6)
      
      DATA THEV/0.6523581397843682D0,1.1071487177940905D0,
     $          1.3820857960113345D0,1.7595068575784587D0,
     $          2.0344439357957027D0,2.4892345138054251D0/
      DATA FIV/0.6283185307179586D0,0.0 D0              ,
     $         0.6283185307179586D0,0.0 D0              ,
     $         0.6283185307179586D0,0.0 D0              /
      DATA FIR/1.2566370614359173 D0/

      IF(LPR)WRITE(6,'(/A)')' ===> Start Subroutine TES '

      CV(1,1)=0.D0
      CV(1,2)=0.D0
      CV(1,3)=1.D0
      CV(32,1)=0.D0
      CV(32,2)=0.D0
      CV(32,3)=-1.D0
      II=1
      DO 520 I=1,6
      TH=THEV(I)
      FI=FIV(I)
      CTH=DCOS(TH)
      STH=DSIN(TH)
      DO 521 J=1,5
      FI=FI+FIR
      IF(J.EQ.1) FI=FIV(I)
      II=II+1
      CV(II,1)=STH*DCOS(FI)
      CV(II,2)=STH*DSIN(FI)
      CV(II,3)=CTH
  521 CONTINUE
  520 CONTINUE
      RETURN
      END
C
      SUBROUTINE DIVIDE(NDIV,LPR)
C     ---------------------------------------------------------------
C     This divides the initial 60 spherical triangles to the level
C     indicated by NDIV
C     ---------------------------------------------------------------
      IMPLICIT NONE

      LOGICAL LPR

      INTEGER*4 IJ
      INTEGER*4 J,J2,J3,J4,J5,JVT1,JVT2
      INTEGER*4 NDIV,NV1,NV2,NV21,NV22,NV23,NV3,NV31,NV32,NV33
      INTEGER*4 NV41,NV42,NV43,NV51,NV52,NV53
     
      REAL*8 XC1,YC1,ZC1
   
      REAL*8 CC,CV,CVN2,CVN3,CVN4,CVN5
      REAL*8 FOUR     
      REAL*8 PI
      REAL*8 XV1,XV2,XV3      
      REAL*8 YV1,YV2,YV3     
      REAL*8 ZERO,ZV1,ZV2,ZV3     

      COMMON/POLI/CV(32,3),XC1(15360),YC1(15360),ZC1(15360)
      COMMON/PENTA/JVT1(3,60),JVT2(3,4)

      DIMENSION CVN2(6,3),CVN3(6,3),CVN4(6,3),CVN5(6,3),CC(3)
 
      DATA ZERO/0.0D0/
      DATA PI/3.1415926535897932D0/
      
      IF(LPR)WRITE(6,'(/A)')' ===> Start Subroutine DIVIDE '

      IJ=0
C*****Level 1****************************
      DO 10 J=1,60
      NV1=JVT1(1,J)
      NV2=JVT1(2,J)
      NV3=JVT1(3,J)
      XV1=CV(NV1,1)
      YV1=CV(NV1,2)
      ZV1=CV(NV1,3)
      XV2=CV(NV2,1)
      YV2=CV(NV2,2)
      ZV2=CV(NV2,3)
      XV3=CV(NV3,1)
      YV3=CV(NV3,2)
      ZV3=CV(NV3,3)
      IF(NDIV.GT.1) GO TO 20
      CALL CALCEN(XV1,YV1,ZV1,XV2,YV2,ZV2,XV3,YV3,ZV3,CC)
      IJ=IJ+1
      XC1(IJ)=SNGL(CC(1))
      YC1(IJ)=SNGL(CC(2))
      ZC1(IJ)=SNGL(CC(3))
      GO TO 10
C*****Level 2**********************
   20 CONTINUE
      CVN2(1,1)=XV1
      CVN2(1,2)=YV1
      CVN2(1,3)=ZV1
      CVN2(2,1)=XV2
      CVN2(2,2)=YV2
      CVN2(2,3)=ZV2
      CVN2(3,1)=XV3
      CVN2(3,2)=YV3
      CVN2(3,3)=ZV3
      CALL CALVER(CVN2)
      DO 21 J2=1,4
      NV21=JVT2(1,J2)
      NV22=JVT2(2,J2)
      NV23=JVT2(3,J2)
      XV1=CVN2(NV21,1)
      YV1=CVN2(NV21,2)
      ZV1=CVN2(NV21,3)
      XV2=CVN2(NV22,1)
      YV2=CVN2(NV22,2)
      ZV2=CVN2(NV22,3)
      XV3=CVN2(NV23,1)
      YV3=CVN2(NV23,2)
      ZV3=CVN2(NV23,3)
      IF(NDIV.GT.2) GO TO 30
      CALL CALCEN(XV1,YV1,ZV1,XV2,YV2,ZV2,XV3,YV3,ZV3,CC)
      IJ=IJ+1
      XC1(IJ)=SNGL(CC(1))
      YC1(IJ)=SNGL(CC(2))
      ZC1(IJ)=SNGL(CC(3))
      GO TO 21
C*****Level 3**********************************
   30 CONTINUE
      CVN3(1,1)=XV1
      CVN3(1,2)=YV1
      CVN3(1,3)=ZV1
      CVN3(2,1)=XV2
      CVN3(2,2)=YV2
      CVN3(2,3)=ZV2
      CVN3(3,1)=XV3
      CVN3(3,2)=YV3
      CVN3(3,3)=ZV3
      CALL CALVER(CVN3)
      DO 31 J3=1,4
      NV31=JVT2(1,J3)
      NV32=JVT2(2,J3)
      NV33=JVT2(3,J3)
      XV1=CVN3(NV31,1)
      YV1=CVN3(NV31,2)
      ZV1=CVN3(NV31,3)
      XV2=CVN3(NV32,1)
      YV2=CVN3(NV32,2)
      ZV2=CVN3(NV32,3)
      XV3=CVN3(NV33,1)
      YV3=CVN3(NV33,2)
      ZV3=CVN3(NV33,3)
      IF(NDIV.GT.3) GO TO 40
      CALL CALCEN(XV1,YV1,ZV1,XV2,YV2,ZV2,XV3,YV3,ZV3,CC)
      IJ=IJ+1
      XC1(IJ)=SNGL(CC(1))
      YC1(IJ)=SNGL(CC(2))
      ZC1(IJ)=SNGL(CC(3))
      GO TO 31
C*****Level 4******************************
   40 CONTINUE
      CVN4(1,1)=XV1
      CVN4(1,2)=YV1
      CVN4(1,3)=ZV1
      CVN4(2,1)=XV2
      CVN4(2,2)=YV2
      CVN4(2,3)=ZV2
      CVN4(3,1)=XV3
      CVN4(3,2)=YV3
      CVN4(3,3)=ZV3
      CALL CALVER(CVN4)
      DO 41 J4=1,4
      NV41=JVT2(1,J4)
      NV42=JVT2(2,J4)
      NV43=JVT2(3,J4)
      XV1=CVN4(NV41,1)
      YV1=CVN4(NV41,2)
      ZV1=CVN4(NV41,3)
      XV2=CVN4(NV42,1)
      YV2=CVN4(NV42,2)
      ZV2=CVN4(NV42,3)
      XV3=CVN4(NV43,1)
      YV3=CVN4(NV43,2)
      ZV3=CVN4(NV43,3)
      IF(NDIV.GT.4) GO TO 50
      CALL CALCEN(XV1,YV1,ZV1,XV2,YV2,ZV2,XV3,YV3,ZV3,CC)
      IJ=IJ+1
      XC1(IJ)=SNGL(CC(1))
      YC1(IJ)=SNGL(CC(2))
      ZC1(IJ)=SNGL(CC(3))
      GO TO 41
C*****Level 5*************************************
   50 CONTINUE
      CVN5(1,1)=XV1
      CVN5(1,2)=YV1
      CVN5(1,3)=ZV1
      CVN5(2,1)=XV2
      CVN5(2,2)=YV2
      CVN5(2,3)=ZV2
      CVN5(3,1)=XV3
      CVN5(3,2)=YV3
      CVN5(3,3)=ZV3
      CALL CALVER(CVN5)
      DO 51 J5=1,4
      NV51=JVT2(1,J5)
      NV52=JVT2(2,J5)
      NV53=JVT2(3,J5)
      XV1=CVN5(NV51,1)
      YV1=CVN5(NV51,2)
      ZV1=CVN5(NV51,3)
      XV2=CVN5(NV52,1)
      YV2=CVN5(NV52,2)
      ZV2=CVN5(NV52,3)
      XV3=CVN5(NV53,1)
      YV3=CVN5(NV53,2)
      ZV3=CVN5(NV53,3)
      CALL CALCEN(XV1,YV1,ZV1,XV2,YV2,ZV2,XV3,YV3,ZV3,CC)
      IJ=IJ+1
      XC1(IJ)=SNGL(CC(1))
      YC1(IJ)=SNGL(CC(2))
      ZC1(IJ)=SNGL(CC(3))
   51 CONTINUE
   41 CONTINUE
   31 CONTINUE
   21 CONTINUE
   10 CONTINUE
      RETURN
      END
C
      SUBROUTINE CALVER(CVN)
C     ---------------------------------------------------------------------
C     This divides one triangle into four.
C     ---------------------------------------------------------------------
      IMPLICIT NONE
      REAL*8 CVN,XXX,YYY,ZZZ,RRR,FC
      INTEGER*4 N,N1,N2
      DIMENSION CVN(6,3)
      DO 7 N=1,3
      N2=N+3
      N1=N-1
      IF(N.EQ.1)N1=3
      XXX=(CVN(N,1)+CVN(N1,1))/2.0D0
      YYY=(CVN(N,2)+CVN(N1,2))/2.0D0
      ZZZ=(CVN(N,3)+CVN(N1,3))/2.0D0
      RRR=SQRT(XXX*XXX+YYY*YYY+ZZZ*ZZZ)
      FC=1.0D0/RRR
      CVN(N2,1)=XXX*FC
      CVN(N2,2)=YYY*FC
      CVN(N2,3)=ZZZ*FC
    7 CONTINUE
      RETURN
      END
C
      SUBROUTINE CALCEN(XV1,YV1,ZV1,XV2,YV2,ZV2,XV3,YV3,ZV3,CC)
C     ---------------------------------------------------------------------
C     This computes the center of a spherical triangle.
C     ---------------------------------------------------------------------
      IMPLICIT NONE
      REAL*8 XV1,YV1,ZV1,XV2,YV2,ZV2,XV3,YV3,ZV3,CC,XXX,YYY,ZZZ,RRR,FC
      DIMENSION CC(3)
      XXX=(XV1+XV2+XV3)/3.0D0
      YYY=(YV1+YV2+YV3)/3.0D0
      ZZZ=(ZV1+ZV2+ZV3)/3.0D0
      RRR=SQRT(XXX*XXX+YYY*YYY+ZZZ*ZZZ)
      FC=1.0D0/RRR
      CC(1)=XXX*FC
      CC(2)=YYY*FC
      CC(3)=ZZZ*FC
      RETURN
      END

C
      SUBROUTINE GEOCAV(NCOR,NP,NDIV,GHOST,LPR,XE,YE,ZE,RE)
C     ------------------------------------------------------------------
C     This computes the surface.
C     ------------------------------------------------------------------
      IMPLICIT NONE

      LOGICAL GHOST
      LOGICAL LPR

      INTEGER*2 IUSE

      INTEGER*4 I,IIJ,IJE,ITO,ISA,ISO
      INTEGER*4 J
      INTEGER*4 K 
      INTEGER*4 L1
      INTEGER*4 MC,MV
      INTEGER*4 N3,N4,NCOR,NDIV,NEJCI,NINF,NP,NSUP,NTRIAN,NTS
      INTEGER*4 UN3

      REAL*8 AP,ATP,ATS
      REAL*8 DD,DIJ2
      REAL*8 FC,FNDIV
      REAL*8 PI
      REAL*8 RE(*),REI,RREJ,RRR
      REAL*8 SRE2
      REAL*8 XC1,XE(*),XEI,XP,XPL,XSL,XSM
      REAL*8 YC1,YE(*),YEI,YP,YPL,YSL,YSM
      REAL*8 ZC1,ZE(*),ZEI,ZP,ZPL,ZSL,ZSM

      REAL*8 CV

      PARAMETER (MC=1000000,MV=1000000)
      
      COMMON/POLI/CV(32,3),XC1(15360),YC1(15360),ZC1(15360)
      COMMON/CSFE/IUSE(MC)
      COMMON/PUN/ITO(MV),ISO(MV),ISA(MV),XP(MV),YP(MV),ZP(MV),AP(MV)

      DIMENSION IJE(MC)

      DATA PI/3.141593E0/

      IF(LPR)WRITE(6,'(/A)')' ===> Start Subroutine GEOCAV '

      NP=0

C begin
      NTRIAN=4**(NDIV-1)
      FNDIV=60*NTRIAN

C It selects one sphere
      DO 1 I=1,NCOR
      
      IF(IUSE(I).EQ.6) THEN
      REI=RE(I)
      XEI=XE(I)
      YEI=YE(I)
      ZEI=ZE(I)
      ATS=4.0E0*PI*REI*REI/FNDIV

      IIJ=0
C
C  It determines which spheres are linked to sphere I
      DO J=1,NCOR

      IF(IUSE(J).GE.3) THEN

      IF(I.NE.J)THEN

      DIJ2=(XEI-XE(J))*(XEI-XE(J))+
     &     (YEI-YE(J))*(YEI-YE(J))+
     &     (ZEI-ZE(J))*(ZEI-ZE(J))
      SRE2=(REI+RE(J))*(REI+RE(J))

          IF(DIJ2.LT.SRE2) THEN
            IIJ=IIJ+1
            IJE(IIJ)=J
            NEJCI=IIJ
          END IF

      END IF
      END IF
      END DO

C 
C It selects one main triangle.
      NSUP=0
      UN3=1
      DO 2 J=1,60
      XPL=0.0E0
      YPL=0.0E0
      ZPL=0.0E0
      NTS=0
      NINF=NSUP+1
      NSUP=NINF+NTRIAN-1
C
C It selects one secondary triangle.
      DO 3 K=NINF,NSUP
      XSL=XC1(K)*REI
      YSL=YC1(K)*REI
      ZSL=ZC1(K)*REI
      XSM=XSL+XEI
      YSM=YSL+YEI
      ZSM=ZSL+ZEI
C
C It fixes if the secondary triangle is inside or outside.
      L1=UN3
      DO N3=L1,NEJCI
      UN3=N3
      N4=IJE(N3)
      DD=(XSM-XE(N4))*(XSM-XE(N4))+
     &   (YSM-YE(N4))*(YSM-YE(N4))+
     &   (ZSM-ZE(N4))*(ZSM-ZE(N4))
      RREJ=RE(N4)*RE(N4)
      IF(DD.LT.RREJ) GO TO 3
      END DO

      DO N3=1,L1-1
      UN3=N3
      N4=IJE(N3)
      DD=(XSM-XE(N4))*(XSM-XE(N4))+
     &   (YSM-YE(N4))*(YSM-YE(N4))+
     &   (ZSM-ZE(N4))*(ZSM-ZE(N4))
      RREJ=RE(N4)*RE(N4)
      IF(DD.LT.RREJ) GO TO 3
      END DO
C
C It prepares the coordinates for the main triangle
      XPL=XPL+XSL
      YPL=YPL+YSL
      ZPL=ZPL+ZSL
      NTS=NTS+1
    3 CONTINUE
C
C It reduces the secondary triangles to the main triangle.
      IF(NTS.EQ.0)GO TO 2

      ATP=ATS*NTS
      XPL=XPL/NTS
      YPL=YPL/NTS
      ZPL=ZPL/NTS
      RRR=SQRT(XPL*XPL+YPL*YPL+ZPL*ZPL)
      FC=REI/RRR

      NP=NP+1
      XP(NP)=XPL*FC+XEI
      YP(NP)=YPL*FC+YEI
      ZP(NP)=ZPL*FC+ZEI
      AP(NP)=ATP
      ITO(NP)=J 
      ISO(NP)=I
      ISA(NP)=I

      
    2 CONTINUE

      END IF

    1 CONTINUE

      RETURN
      END
C
      SUBROUTINE VOLARE(NP,LPR,GHOST,STOT,VOL,XE,YE,ZE,RE)
C     ---------------------------------------------------------------------
C     This calculates total area and volume.
C     ---------------------------------------------------------------------

      IMPLICIT NONE    

      LOGICAL LPR,GHOST

      INTEGER*2 IUSE

      INTEGER*4 I,ISA,ISO,ITO
      INTEGER*4 J,JU
      INTEGER*4 NP
      INTEGER*4 MV,MC

      REAL*8 AP
      REAL*8 DP
      REAL*8 RE(*),REJ
      REAL*8 VD,VN
      REAL*8 XE(*),XP,XEJ
      REAL*8 YE(*),YP,YEJ
      REAL*8 ZE(*),ZP,ZEJ

      REAL*8 STOT,VOL

      PARAMETER (MC=1000000,MV=1000000)

      COMMON/CSFE/IUSE(MC)
      COMMON/PUN/ITO(MV),ISO(MV),ISA(MV),XP(MV),YP(MV),ZP(MV),AP(MV)

      IF(LPR)WRITE(6,'(/A)')' ===> Start Subroutine VOLARE '

      STOT=0.0E0
      VOL=0.0E0

      DO I=1,NP
      STOT=STOT+AP(I)
      END DO


      IF(.NOT.GHOST)THEN
      
        JU=0
        DO I=1,NP
        J=ISO(I)

           IF(JU.NE.J)THEN
            JU=J
            XEJ=XE(J)
            YEJ=YE(J)
            ZEJ=ZE(J)
            REJ=RE(J)
           END IF

        VOL=VOL+((XP(I)-XEJ)*XP(I)+
     &           (YP(I)-YEJ)*YP(I)+
     &           (ZP(I)-ZEJ)*ZP(I)) *AP(I)/(3.0E0*REJ)


        END DO

      END IF


      RETURN 
      END
C


      SUBROUTINE SUME(NCOR,RD,OP,LPR,RE)
C     --------------------------------------------------------------------
C     Add the solvent radius to every sphere raddi
C     Or subtract the solvent radius.
C     --------------------------------------------------------------------
      IMPLICIT NONE

      LOGICAL LPR

      INTEGER*2 IUSE

      INTEGER*4 I
      INTEGER*4 MC
      INTEGER*4 NCOR

      REAL*8 RD,RE(*)
      
      CHARACTER*4 OP
 
      PARAMETER (MC=1000000)

      COMMON/CSFE/IUSE(MC)     

      IF(LPR)WRITE(6,'(/A)')' ==> Start subroutine SUM'

      IF(OP.EQ.'SUMA')THEN

        IF(LPR)WRITE(6,'(A)')' Adding the radius of the solvent'
        DO I=1,NCOR
         IF(IUSE(I).NE.1)RE(I)=RE(I)+RD
        END DO

      ELSE IF(OP.EQ.'REST')THEN

        IF(LPR)WRITE(6,'(A)')' Subtracting the radius of the solvent'
        DO I=1,NCOR
         IF(IUSE(I).NE.1)RE(I)=RE(I)-RD
        END DO

      ELSE

        WRITE(6,'(A)')' Variable OP badly defined'
        STOP

      END IF
 
      RETURN
      END


      SUBROUTINE PRISPH(NATOM,NCOR,NP,REDU,LPR,AE)
C     ---------------------------------------------------------------
C     This prints the coordinates, radius, label and area of each
C     sphere
C     ---------------------------------------------------------------
      IMPLICIT NONE

      LOGICAL LPR
      LOGICAL REDU

      INTEGER*2 IUSE

      INTEGER*4 I,ISA,ISO,ITO
      INTEGER*4 J
      INTEGER*4 L,LT
      INTEGER*4 MC,MV
      INTEGER*4 NATOM,NCOR,NP

      REAL*8 AE(*),AP
      REAL*8 XP
      REAL*8 YP
      REAL*8 ZP

      CHARACTER*4 LNEW
      PARAMETER (MC=1000000,MV=1000000)

      COMMON/CSFE/IUSE(MC)
      COMMON/PUN/ITO(MV),ISO(MV),ISA(MV),XP(MV),YP(MV),ZP(MV),AP(MV)


      IF(LPR)WRITE(6,'(/A)')' ===> Start Subroutine PRISPH '

      DO J=1,NCOR
      AE(J)=0.0E0
      END DO

      DO I=1,NP
      J=ISA(I)
      AE(J)=AE(J)+AP(I)
      END DO

      RETURN
      END

