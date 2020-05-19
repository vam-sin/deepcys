      SUBROUTINE FSASA(NRES,IEPS,DEBUG,IDB,INGRP,NATM,
     *VEPS,BFS,NGRP,PCOD)
C
C	calculate solvent accessible surface area of charged groups and
C	fraction relative to group in H2O
C
C	elm, Feb 10, 1997
C
C	DGS Fraction solvent accessible surface area
C	1.-DGS is fraction buried in protein
C	DG2 transfer energy for unit charge
C	ITYPE =2 for titratable group, =1 otherwise
C
      IMPLICIT REAL*8 (A-H,O-Z)
c     parameter (DST=78.4)
c     CHARACTER RES*3,line*80,EPSTYP(2)*4,PCOD*(*)!,temp*4
      CHARACTER PCOD*(*)
c     REAL*8 NTERM,MTOTSA
      integer lres
      LOGICAL DEBUG,INGRP(1),VEPS
      include 'esp.inc'
      include 'chd.inc'
      include 'screen.inc'
c     DATA EPSTYP/'D 1 ','D 2 '/
       write(*,*)'inside FSASA'
      TRF=0
      FRF=0
      DO I=1,NATM
      DG2(I)=0.
      ENDDO
C
      LRESO=0
      ITYP=0
      itypo=ityp
      IATM=0
      LGR=0
      LLGR=0
c     lres=1 !residue number
C
C	get atom SASA's and calculate microenvironment properties
C
       write (6,1101) 
 1101   FORMAT(/8x,'RES NO      RES  GRP NO      FB  ',
     *  '       RFHP ',
     *  '     RFHT      RFHT/RFW')
      IST=1
      write(*,*)"in resfsa, ngrp",ngrp
      DO I=1,NGRP
        lres = ares(ist)
        ISTP=IST+NGRPS(2,I)-1
        ITYP=NGRPS(3,I)
        write(*,*)"Within resfsa ITYP, NGRPS(3,I)",
     *  ITYP, NGRPS(3,I)
        if (lres.eq.lreso) then 
         llgr=llgr+1
        else
         llgr=1
         lreso = lres
        endif
        write(*,*)"NGRPS(3,I),(2,I),(1,I)",NGRPS(3,I),NGRPS(2,I)
     *  ,NGRPS(1,I)
 1100  FORMAT (2(I10,2x),3(I2,3x),2x,F8.4)
        FB=1.-DGS(IST)
        TRF=DGS(IST)*SMENV(LLGR,ITYP)+FB*RKRH(I)
        IF (SMENV(LLGR,ITYP).NE.0) THEN
          write(*,*)"LLGR,ITYP,SMENV(LLGR,ITYP),RKRH(I)",
     *    LLGR,ITYP,SMENV(LLGR,ITYP),RKRH(I)
          FRF=TRF/SMENV(LLGR,ITYP)
        ELSE
          FRF=-1.
        ENDIF
        WRITE (30,1103) PCOD,LRES,RESGRP(LNRES(LRES)),LLGR,
     *  GTYP(LLGR,LNRES(LRES)),FB,
     *  RKRH(I),TRF,FRF
 1103   FORMAT(A4,',  ',I4,',  ',A3,',  ',I6,',  ',
     *  A2,',  ',F6.3,',  ',2(F7.3,',  '),F7.3)
c     IST=ISTP+1
      ist=ist+ngrps(2,i)
      ENDDO
      CLOSE(30)
c     RETURN
 9001 WRITE (*,*) ' end of file on unit 2 in resfsa'
      RETURN
 9002 WRITE (*,*) 'NOTHING IN UNIT 3'
      STOP 2
      END
