C
C	general common and dimensions
C
C     LATMX          = maximum number of atoms
C character
C     IARES(LATMX)    = residue number in characters
C     ATM(LATMX)     = atom name in characters
C     RESGRP(IGRMX)  = residue name
C     GTYP(5,IGRMX)  = (1-5 possible groups per residue, residue number)
C integer
C     NGRPS(1,IGMX)  = charge type (1 if neutral, 2 if charged) 
C     NGRPS(2,IGMX)  = # atoms in group
C     NGRPS(3,IGMX)  = residue type, for each group in residue
C     LNRES(MXRES)   = residue type (0 for residues not part of pK calc
C     NG(LATMX)      = group number for each atom. 
C     NRATM(MXRES)   = first atm number of each residue
C real*8  
C     RGRP(3,IGMX)   = (X,Y,Z of group center) 
C common/data/
C     R(6,LATMX)     = atom properties (X,Y,Z,
C        For nonprotonatable residues charges are read into
C       R(4,I); for titratable residues, charges for charged state are
C       read into R(5,I) & charges for neutral state are read into R(4,I)
C        Hydrophobic constants are read into R(6,I)
C     RMIN(3)        = minimum xyz of system
C     RMAX(3)        = maximum xyz of system
C     RDUM(3)        = legacy variable from pk calc. 
C common/res/
C     RESGRP         = resnames of residues that can be protonated, see espmain.f
C     SMENV(IGRMX,50)= solvent microenvironment data per residue
C     NGRES(50)      = residue group number
C common/dgepol
C     XG(LATMX)      = X coord of group center
C     YG(LATMX)      = Y coord of group center
C     ZG(LATMX)      = Z coord of group center
C     RV(LATMX)      = atom radius for sasa calculation
C     RVX(LATMX)     = RV(LATMX) mostly. used in coorin. check before using.
C     
      PARAMETER (LATMX=1000000,IGMX=1000000,MXRES=1000000,IGRMX=33)
      CHARACTER ATM(LATMX)*4,RESGRP(IGRMX)*3,GTYP(5,IGRMX)*2
      INTEGER*2 NGRPS(4,IGMX),LNRES(MXRES),NG(LATMX),NRATM(MXRES)
      INTEGER::IARES(LATMX) !remodified to integer 10 on Oct 11th 2011
c      INTEGER*8 ARES(LATMX) !modified to kind on Oct 9th 2011
c     #INTEGER (KIND=8) ARES(LATMX) #modified to kind on Oct 9th 2011
c      INTEGER*8 ARES(LATMX) #originally integer length mentioned as 10; 
c     however, for compiler compatibility dimension reduced to 22nd July, 2011
      REAL*8 RGRP(3,IGMX)
      COMMON/DATA/R(6,LATMX),RMIN(3),RMAX(3),RDUM(3),NGRPS,RGRP,IARES,
     *ATM,LNRES,NG
      COMMON/RES/SMENV(IGRMX,50),NGRES(5000),GTYP,RESGRP
      COMMON/DGEPOL/XG(LATMX),YG(LATMX),ZG(LATMX),RV(LATMX),
     *RVX(LATMX),NRATM

