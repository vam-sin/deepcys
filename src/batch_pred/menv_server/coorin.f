C
C	READ IN ATOMIC COORDINATES AND FIND MIN AND MAX
C	VALUES FOR UNIT CHARGE MODEL, ASSIGN CHARGES AND
C	COORDINATES
C
      SUBROUTINE COORIN(DEBUG,LATM,IGRP,INPT,IOUT,PATM,PRES,
     *NRES,Q,IPT,FILTP,INGRP,IRES,NPRES,NTRES,Q0,PCOD)
C
C	FILTP input file tyoe, i.e., origin
C	DEBUG extra output
C	THE vdW radii values are squared
C
C       NTATM        = total number of atoms in system. 
C       NTRES        = total number of residues in system
C       NGRD   = total number of grid residues in system.
C       Assume grid points are at the end of the  coordinate file. 
C
C	NGRPS(1,IGRP) = charge type(1 if neutral, 2 if charged)
C       NGRPS(2,IGRP) = # of atoms in group
C       NGRPS(3,IGRP) = residue type, for each group in residue
C       LNRES(IRES) = residue type (0 for residues not part of pK calc)
C
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL DEBUG,INGRP(1),GRDTEST
      INTEGER NGRD,NGRDATMST,tres,treso
      dimension NGRPS3(100)
      CHARACTER LINE*150,RES(100000)*3,PRES(1)*(*),PATM(1)*(*),
     *NRES(1)*(*),FILTP*(*),PCOD*(*),CHMTYPE*4,
     *F1000*120,LINE1(100000)*150
      REAL*8 Q(1),FCH0(9)
      INCLUDE 'esp.inc'
      INCLUDE 'chd.inc'
      DATA FCH0/9*0./,RPR/1.0D0/
c       write(*,*)'entering the FILOPN'
      CALL FILOPN(INPT,IOUT,PCOD)
c       write(*,*)'entering the RDCHG'
C     CALL RDCHG(R,NG,RV,CHMTYPE)
       CALL RDCHG(R,NG,RV,FILTP)
C     CALL RDCHG(CHMTYPE)
C       write(*,*)'entering the RDSENV'
      CALL RDSENV
C       write(*,*)'out of RDSENV'
        I=0
c       IRA=0
c       IRES=1  
        NPATM=0
        NGRD=0
        NGRDATMST=0
        TRESO=0
C
C----------------Reading number of atoms in file-----------------------------
       IF (FILTP.EQ.'CHM') THEN
   30   READ (1,'(A)',END=900,IOSTAT=IERR) LINE
        IF (IERR.NE.0) GO TO 900
          if (line(:1).eq.'*') then
c         write (iout,'(a)') line
          go to 30
          endif
        READ (LINE(1:5),'(I5)') NTATM
       ELSE IF (FILTP.EQ.'PDB') THEN
        LN0=0
        IPDB=0
        Do nn=1,10000000
   40   READ (1,'(A)',END=41,IOSTAT=IERR) LINE
        IF (IERR.NE.0) GO TO 900
c       IF (LINE(:6).EQ.'ATOM  '.OR.LINE(:6).EQ.'HETATM') IPDB=IPDB+1
        IF (LINE(:6).EQ.'ATOM  ') then
         IPDB=IPDB+1
c        Line1(ipdb)=line
c        write(*,'(a150)')line1(ipdb)
        ENDIF
        IF (IPDB.EQ.0) LN0=LN0+1
        enddo
c       GO TO 40
   41   NTATM=IPDB
c       REWIND 1
       ENDIF
c      write(*,*)NTATM
C--------------------End of reading file to retrive total # of atoms-----------
C
C--------------------Initializing some parameters------------------------------
 1013 FORMAT(' atoms in protein =',I5,' total number of atoms =',I5)
      IPPT=IPT
      NEUT=0
      NCHG=0
      IGRP=1
      IRES=1
      IRA=0
      Q0=0.
C-------------------END of initialization--------------------------------------
C
C-------------Reading file to retrive coordinates and other data---------------
c1     Continue
      DO i = 1,1000000
   10  READ (1,'(A)',END=2,IOSTAT=IERR) LINE
       IF (IERR.NE.0) GO TO 900
       IF (FILTP.EQ.'CHM') THEN
          IRA=IRA+1
          READ (LINE,1000) aRES(i),RES(i),ATM(I),(R(J,I),J=1,3),dgs(i)
c         write (15,'(a)') line
c         write(*,*)"i and NG(i)",i,NG(i)
          TRES = ares(i)
c         write(*,*)"writing TRES", TRES, ares(i)
          if(i.eq.1) TRESO = TRES
       ELSE IF (FILTP.EQ.'PDB') THEN
          READ (LINE,1100) ATM(I),RES(i),aRES(i),(R(J,I),J=1,3)
c         write(*,1100) ATM(I),RES(i),aRES(i),(R(J,I),J=1,3)
       ENDIF
C
C---------------Reading number of atoms per residue-------------------------------
        IF (TRES.NE.TRESO) THEN
          NRATM(IRES+1)=IRA-1
          IRA=1
          IRES=IRES+1
          TRESO=TRES
        ENDIF
       NRATM(IRES)=IRA
       NTRES=IRES
c        write(*,*)"IRA",ira
c       write(*,*)"IRES ",irES," nratm(ires)",nratm(ires)
C---------------END of reading number of atoms in group------------------------
C
C---------------Creating new identifiers- for groups--------------------------------------
c       ILST=1
        DO II=1,IGRMX
          If (res(i).eq.resgrp(ii)) go to 22
c         write(*,*)"res(i),atm(i),NG(i)",res(i), atm(i),NG(i)
        ENDDO
22      XG(I)=R(1,I)
        YG(I)=R(2,I)
        ZG(I)=R(3,I)
        IG=NG(I)
c       DO II=1,IGRMX
c       IF (RES(i).EQ.RESGRP(II).AND.INGRP(II)) THEN
c         write(*,'(a3)') "RESGRP(II)",RESGRP(II),"RES",res
c         NGRPS(3,IGRP)=II
c         LNRES(IRES)=II
c         write(*,*)"ngrps",ngrps(3,igrp)
c         GO TO 20
c         ENDIF
c       ENDDO

c       write(*,*)"i,ng(I)",i,ng(i)
        IF (ATM(I)(:3).EQ.'OT1' .or. ATM(I)(:3).EQ.'OXT'
     *  .or.ATM(I)(:2).EQ.'NT') then 
             NTGRP=IGRP
c            write(*,*)"ntgrp--", ntgrp
        ENDIF
c       write(*,*)"IG",IG
        DO IX=1,IGRMX
c           write(*,'(a4)') resgrp(ix)
c        IF (RES(i).EQ.RESGRP(IX).and.ingrp(ix).and.IG.eq.1) THEN
         IF (RES(i).EQ.RESGRP(IX)) THEN
            NGRPS(3,IGRP)=IX
            LNRES(IRES)=IX
c           write(*,*)"ix,i, NGRPS(3,IGRP),LNRES", ix,i,
c    *      NGRPS(3,IGRP),LNRES(IRES)
          GO TO 20
         ENDIF
        ENDDO
20      IF (IG.EQ.1) THEN ! edited on 8th July 2015 by DB
c          NGRPS(3,IGRP+1)=0
c          LNRES(IRES)=0
           NGRAT=I-ILST
c          write(*,*)"i,ilst, NGRAT=I-ILST",i,ilst,NGRAT
c          DO IJ=1,3 !inactivated on 16/11/15
c             RGRP(IJ,IGRP)=RGRP(IJ,IGRP)/NGRAT !inactivated on 16/11/15
c          ENDDO !inactivated on 16/11/15
c      	   NGRPS(1,IGRP)=ITYP
           NGRPS(1,IGRP)=2
           NGRPS(2,IGRP)=NGRAT
c          NGRPS(3,IGRP)=IX
c          write(*,*) "NGRPS(2,IGRP),(3,IGRP), I, ILST, 1",
c    *     NGRPS(2,IGRP),NGRPS(3,IGRP),I,ILST
           IGRP=IGRP+1
           IF (IGRP.GT.IGMX) THEN
             WRITE (*,*) ' TOO MANY GROUPS, EXECUTION TERMINATED,', 
     *      ' NO. OF ATOMS =',I
             CALL EXIT
           ENDIF
c          IGRP=IGRP-1
           ILST=I
c          DO IJ=1,3 !inactivated on 16/11/15
c              RGRP(IJ,IGRP-1)=R(IJ,I) !inactivated on 16/11/15
c          ENDDO !inactivated on 16/11/15
        ELSE IF (IG.EQ.0) THEN
c          DO IJ=1,3! inactivated on 16/11/15
c              RGRP(IJ,IGRP-1)=RGRP(IJ,IGRP-1)+R(IJ,I)! inactivated on 16/11/15
c          ENDDO! inactivated on 16/11/15
        ENDIF
c        DO IJ=1,3
c          RGRP(IJ,IGRP-1)=R(IJ,I)
c        ENDDO
c        ILST=1
c        NGRPS(3,IGRP-1)=IX
        NGRPS3(IGRP-1)=NGRPS(3,IGRP-1)
c       write(*,*)"NGRPS(3,IGRP) out of the loop in coorin",
c    *  NGRPS(3,IGRP-1)
c    * (NGRPS(3,IT),IT=1,IGRP-1)
c    *   (NGRPS3(IT), IT=1,IGRP-1)
         DO II=1,3
           RMAX(II)=R(II,I)
           RMIN(II)=R(II,I)
         ENDDO

            DO II=1,3
            IF (RMAX(II).LT.R(II,I)) RMAX(II)=R(II,I)
            IF (RMIN(II).GT.R(II,I)) RMIN(II)=R(II,I)
            ENDDO
      ENDDO
c     write(*,*)"i at the end of Do loop", i
C---------------END of reading for retriving coordinates-----------------------
C----------------------------END of DO Loop for reading atom and groups--------
C
C--------------Start-calculating properties per group per residue--------------
C
    2 LATM=I-1
C     C-------------Initialization of some parameters--------------------------
      NGRAT=I-ILST+1  !!! NGRAT == Number of atoms per group in a residue
c     DO IJ=1,3! inactivated on 16/11/15
c             RGRP(IJ,IGRP-1)=RGRP(IJ,IGRP-1)/NGRAT! inactivated on 16/11/15
c     ENDDO! inactivated on 16/11/15
c             NGRPS(1,IGRP-1)=ITYP
c     NGRPS(1,IGRP-1)=2 ! commented on 25/11/15
c     NGRPS(2,IGRP-1)=NGRAT ! commented on 25/11/15
c     write(*,*) "NGRPS, I, ILST, IGRP-1",
c    *NGRPS(2,IGRP-1),NGRPS(3,IGRP-1),I,ILST, igrp
      WRITE (IOUT,1010) NTATM,RMIN,RMAX,IGRP-1,NTGRP
      IF (NEUT*NCHG.GT.200000) THEN
        WRITE (*,*) 'NEUT*NCHG exceeds maximum(200000)'
        stop 12345
      ENDIF
C     C--------------End of initialization-------------------------------------
C
C     C--------------Calculating accessible surface area per residue per group-

      DO I=1,NTATM
        DGS(I)=1.0+E0
        DG(I)=1.0+E0
        RVX(I)=RV(I)  !!! RV(I)==
      ENDDO
      CALL GEPOL(NTATM,IOUT,XG,YG,ZG,RVX,NG,DGS,DEBUG,RPR)
c     write(*,*)"DGS from Gepol",DGS
      IS0=1
      IS1=0!changed from 0 to 1 on 25/11/15
      NTRES=ares(LATM) !!! NTRES== Number of residues
C
c     write(*,*)"NTRES,NGRD,LATM",NTRES,NGRD,LATM
C
      do itest=1,igrp-1
c     write(*,*)"itest in coorin outside igrp loop,ngrps(2,itest)",
c    * itest,ngrps(2,itest)
      enddo
      DO Ip=1,NTRES-NGRD
c     DO I=1,NTRES!deactivated on 21/11/15
c     DO Ip=1,igrp-1!reactivated on 21/11/15
c       write(*,*)"WITHIN coorin, NGRPS(3,I)", NGRPS(3,I)
cc     IS1=nratm(I-1)+IS1 ! 21st Aug 2015
        IS1=nratm(Ip)+IS1! deactivated on 21/11/15
c       IS1=ngrps(2,ip)+IS1 !changed on 21st Aug 2015, relevant when considering per group; reactivated on 21/11/15
c       write(*,*)"ngrps(2,ip) at the end, IS1",ngrps(2,ip),IS1
        NGRPS(3,I)=NGRPS3(Ip)
c       write(*,*)"WITHIN coorin, NGRPS(3,I)", NGRPS(3,I)
        DO II=1,NTATM  !!! NTATM== Number of atoms
          RVX(II)=0.
        ENDDO
c       write(*,*)"IS0,IS1",is0,is1
        DO II=IS0,IS1
          RVX(II)=RV(II)
c         write(*,*)"II,RVX(II)",ii,rvx(ii)
        ENDDO
        CALL GEPOL(ntatm,IOUT,XG,YG,ZG,RVX,NG,FR,.FALSE.,RPR)
c       CALL GEPOL(ntatm,IOUT,XG,YG,ZG,RVX,NG,FR,.FALSE.,RPR)
c       if (ip.lt. 4) write(*,*) "FR from Gepol",FR
        DO II=IS0,IS1
          DG(II)=FR(II)
          write(*,*)"ip, ares(ip),II, DG(II)",ip,ares(ip),ii,dg(ii)
        ENDDO
        IS0=IS1+1
      ENDDO
      IF (DEBUG) write (IOUT,*) ' ATM   exps    exps0   expf'
c     DO I=1,igrp-1 !changed from latm to igrp-1 on 25th Nov 2015
      DO I=1,ntatm !changed from NTATM to LATM on 3rd Sep 2015
c       IF (DEBUG) THEN
          EXPF=DGS(I)/DG(I)
          write (*,'("DGS and DG",I5,3f10.4)') I,DGS(I),DG(I),EXPF
c       ENDIF
        DGS(I)=DGS(I)/DG(I)
      ENDDO
      IGRP=IGRP-1
      RETURN
  900 WRITE (IOUT,'(I5)') ' ERROR ON READ, IOERR = ',IERR
      CALL EXIT
c1000 FORMAT(6x,i4,1X,a3,2x,a4,3F10.5)
 1000 FORMAT(6x,i4,1X,a3,2x,a4,3F10.5,12x,f8.6)
c1100 FORMAT(10x,a10,2X,a3,7x,a4,4x,3F20.10)
 1100 FORMAT(12x,A4,1x,A3,2x,i4,4x,3F8.3)
 1001 FORMAT(1X,A4,A3,2X,A4,6X,4F8.3,I5)
 1002 FORMAT(A4,1X,3F15.9,1X,A3,1X,A4,5X,I1,F10.4)
 1010 FORMAT(/' NO OF COORDINATES = ',I10,' RMIN: ',3F8.3/
     *28X,'RMAX: ',3F8.3/
     *'TOTAL GROUPS =',I10/
     *' NO OF GROUPS IN PROTEIN =',I10)
 1011 FORMAT(10X,' COORDINATES,'//(6(I10,I2,2X,F7.3)))
 1012 FORMAT(10X,' GROUP  ITYP   NO OF ATOMS  RES TYPE Grp Center'//
     *(I6,3X,I7,4x,I7,4x,i7,3f10.4))
      END
      