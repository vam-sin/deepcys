      SUBROUTINE CNTRL(INPT,IOUT,IEPS,GI,RCUT,
     *DEBUG,NGRP,NTATM,PHS,PHF,DPH,TEMP,INGRP,IRES,Q0,IDB,
     *ITERM,BFS,ALF,VEPS,RESET,CRCUT,CDAMP,RENV,NROF,PCOD,ISGRD)
C	
C	control program for iterative procedure to calculate average
C	charge distribution of acidic & basic groups and pK's 
C	using a variational procedure
C
C	E.L. MEHLER March 24, 1993
C
C	SCREENING FUNTION TYPE:
C 		CONSTANT = 1
C		LINEAR   = 2
C		WDS      = 3
C		SCL      = 4
C	FOR TYPES 3 AND 4 SEE MEHLER&EICHELE, BIOCHEMISTRY,1984
C	DEFAULT RCUT(1)=15A; RCUT(2)=30A
C	DEFAULT SETTINGS: FILTYP=DIS,SCRN FNCTN=3, I=0., DGRID=1A,EXTENT=4.A
C	eps(WDS for r>20A) = 77.8
C
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER ANS*1,NRES*4,PCOD*(*)
      LOGICAL DEBUG,INGRP(1),VEPS,RESET,LCHG,LEXP
      REAL*8 DR(3),KAPPA,
     *RCUT(2),RCUT0(2),PK0(9),FCH0(9)
      PARAMETER (zro=0.,EAWDS=76.8,EPS=78.4)
      INCLUDE 'esp.inc'
      INCLUDE 'chd.inc'
      INCLUDE 'screen.inc'
      DATA KAPPA/.32915/,
     *SR1,SR2/10.,30./
      DATA PK0/6.3,4.4,4.0,10.0,10.4,12.0,0.,0.,8.3/
      DATA FCH0/1.,-1.,-1.,-1.,1.,1.,-1.,-1.,-1./
      DATA APK0,CPK0/8.0,3.1/
C
      write(*,*)'Inside cntrl.f'
      
c     IF (DEBUG.AND.IDB.GE.2) !commented on 25.11.15
      WRITE (88,1001) (IGRP,NGRPS(1,IGRP),NGRPS(2,IGRP),
     *NGRPS(3,IGRP),IGRP=1,NGRP)
 1001 FORMAT(4I8)
      IF (DEBUG.AND.IDB.GE.1) 
     *WRITE (IOUT,'(5(2I8))') (IIRES,LNRES(IIRES),IIRES=1,IRES)
      JTRES=IRES
      IF (GI.NE.0) GI=SQRT(GI)
      write(*,*)'entering env.f'
      CALL ENV(INPT,IOUT,RENV,DEBUG,NGRP,NTATM,IDB,JTRES-NROF,ISGRD)
      write(*,*)'entering FSASA'
      CALL FSASA(JTRES-NROF,IEPS,DEBUG,IDB,INGRP,NTATM,
     *VEPS,BFS,NGRP,PCOD)
      write(*,*)'out of FSASA'
      return
C     stop 1234
      end
