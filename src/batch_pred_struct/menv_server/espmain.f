C
C	MAIN PROGRAM TO CALULATE THE ELECTROSTATIC POTENTIAL
C
C	E.L. MEHLER  MAY 9, 1990, REVISED NOV 29, 1990
C	Revised to calculate pK only April, 1997
C
C	CALCULATE SELF-CONSISTENT CHARGE DISTRIBUTION OF
C	MACROMOLECULE AND pK OF IONIZABLE GROUPS.
C
C
C	SCREENING FUNCTION 1 = CONSTANT; 
C                          2 = LINEAR; 
C                          3 = WDS;
C	                   4 = SCALED
C			   5 Combined 3 & 4 according to SASA       
C
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER FILTP*3,ANS*1,ATMM(100)*4,RES(100)*3,
     *NRES(100)*4,UNT*2,PCOD*4
      LOGICAL DEBUG,INGRP(1),VEPS,RESET
      REAL*8 RCUT(2),RCUT0(2),Q(100),RMIN0(3),RMAX0(3)
      REAL*4 ti,tim(2),dtime
      PARAMETER (INPT=5,IOUT=6)
      INCLUDE 'esp.inc'
      INCLUDE 'screen.inc'
      DATA DEBUG/.FALSE./,RCUT0/100.,900./,IPT,IPT2/0,0/,Q/100*0./
      DATA VEPS/.FALSE./,reset/.true./   
      DATA IDB/0/
C     DATA NGRES/3,4,3,6,4,4,3,4,3,2,3,3,4,3,6,5,5,4,3,4/
C     DATA SMENV/0.,-18.453,6*0.,-12.832,4*0.,-13.492,8*0.,-8.368,3*0.,
C    *-11.420,5*0.,-11.420,3*0.,84*0./
C
C       RESGRP = code for residues which can be protonated
C       and/or are charged
C	In the program the residue type is stored as the position
C	in array RESGRP
C
C	INGRP = code of ionizable groups included in pK calculation
C	default is GLU, ASP, LYS, ARG and termini
C
C 1 'HIS'	'AM','HS','CM',2*'  ',
C 2 'HSE'	'AM','HS','CM',2*'  ',
C 3 'HSD'	'AM','HS','CM',2*'  ',
C 4 'GLU'	'AM','HP','CO','CM','  ',
C 5 'ASP'	'AM','CO','CM',2*'  ',
C 6 'TYR'	'AM','HP','RS','PH','CM',
C 7 'LYS'	'AM','HP','AS','CM','  ',
C 8 'ARG'	'AM','HP','GS','CM','  ',
C 9 'SER'	'AM','OL','CM',2*'  ',
C10 'THR'	'AM','OL','HP','CM','  ',
C11 'CYS'	'AM','TO','CM',2*'  ',
C12 'GLY'	'AM','CM',3*'  ',
C13 'ALA'	'AM','HP','CM',2*'  ',
C14 'LEU'	'AM','HP','CM',2*'  ',
C15 'ILE'	'AM','HP','CM',2*'  ',
C16 'VAL'	'AM','HP','CM',2*'  ',
C17 'TRP'	'AM','HP','RS','RS','CM',
C18 'GLN'	'AM','HP','AD','CM','  ',
C19 'ASN'	'AM','AD','CM',2*'  ',
C20 'PHE'	'AM','HP','RS','CM','  ',
C21 'PRO'	'AP','HP','CM',2*'  ',
C22 'MET'	'AM','HP','TE','CM','  ',
C23 'HSP'	'AM','HS','CM',2*'  ',
C24 'HSC'	'AM','HS','CM',2*'  ',
C25 'GDP'	5*'  ',
C26 'GTP'	5*'  ',
C27 'GDQ'	5*'  ',
C28 'GTQ'	5*'  ',
C29 'GDR'	5*'  ',
C30 'GTR'	5*'  ',
C31 'OCT'	'HP',4*'  ',
C32 'ALY'	'AM','HP','AS','CM','HP',
C33 'GRD'	'GR',4*'  '/
C34 'HOH'       'WT',4*'  '/

      DATA RESGRP/
     *'HIS','HSE','HSD','GLU','ASP',
     *'TYR','LYS','ARG','SER','THR',
     *'CYS','GLY','ALA','LEU','ILE',
     *'VAL','TRP','GLN','ASN','PHE',
     *'PRO','MET','HSP','HSC','GDP',
     *'GTP','GDQ','GTQ','GDR','GTR',
     *'OCT','ALY','GRD'/ !,'HOH'/
      DATA GTYP/
     *'AM','HS','CM',2*'  ',
     *'AM','HS','CM',2*'  ',
     *'AM','HS','CM',2*'  ',
     *'AM','HP','CO','CM','  ',
     *'AM','CO','CM',2*'  ',
     *'AM','HP','RS','PH','CM',
     *'AM','HP','AS','CM','  ',
     *'AM','HP','GS','CM','  ',
     *'AM','OL','CM',2*'  ',
     *'AM','OL','HP','CM','  ',
     *'AM','TO','CM',2*'  ',
     *'AM','CM',3*'  ',
     *'AM','HP','CM',2*'  ',
     *'AM','HP','CM',2*'  ',
     *'AM','HP','CM',2*'  ',
     *'AM','HP','CM',2*'  ',
     *'AM','HP','RS','RS','CM',
     *'AM','HP','AD','CM','  ',
     *'AM','AD','CM',2*'  ',
     *'AM','HP','RS','CM','  ',
     *'AP','HP','CM',2*'  ',
     *'AM','HP','TE','CM','  ',
     *'AM','HS','CM',2*'  ',
     *'AM','HS','CM',2*'  ',
     *30*'  ',
     *'HP',4*'  ',
     *'AM','HP','AS','CM','HP',
     *'GR',4*'  '/
c    *'WT',4*'  '/
c     DATA INGRP/IGRMX*.TRUE./
      DATA INGRP(1)/.TRUE./
C
C	INPUT
C
C
C	PHS	Starting pH
C	PHF	Final pH
C	DPH	pH steps
C	TEMP	Temperature (K)
C
C
C	Coordinate file type:	CHM Charm input coordinates
C				PDB Protein data bank coordinates (not yet implemented)
C
C	ANS,IDB Debug control swithches
C
      WRITE (IOUT,1001)
 1001 FORMAT('$ENTER COORDINATE FILE TYPE',
     *' (PDB,CHarMm): ')
      READ (INPT,'(A3,I4,A,I1)') FILTP,ISGRD,ANS,IDB
      IF (ANS.EQ.'D') DEBUG=.TRUE.
      RENV=10.
        ti=dtime(tim)
C
C	FETCH COORDINATES AND CHARGES
C
        write(*,*)'enteringin coorin.f'
        CALL COORIN(DEBUG,NATM,NGRP,INPT,IOUT,ATM,RES,
     *  NRES,Q,IPT,FILTP,INGRP,IRES,NPRES,NTRES,Q0,PCOD)
        write(*,*)'exiting coorin.f'
        write(*,*)'entering cntrl.f'
        CALL CNTRL(INPT,IOUT,IEPS,GI,RCUT,DEBUG,
     *  NGRP,NATM,PHS,PHF,DPH,TEMP,INGRP,IRES,Q0,IDB,
     * ITERM,BFS,ALF,VEPS,RESET,CRCUT,CDAMP,RENV,NTRES-NPRES,PCOD,ISGRD)
        write(*,*)'exiting cntrl.f'
      WRITE(IOUT,'(/)')
        ti=dtime(tim)
        write (iout,*) ' elapsed cpu time =',ti,' sec.'
      stop 1234
      CALL EXIT
      END
