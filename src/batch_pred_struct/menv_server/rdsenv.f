C
C	read in solvent microenvironment data
C
      SUBROUTINE RDSENV
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER LINE*120
      INCLUDE 'esp.inc'
      DO I=1,IGRMX
      NGRES(I)=0
      DO J=1,35
      SMENV(J,I)=0.
      ENDDO
      ENDDO
   1  READ (3,'(a)',END=900) LINE
      IPCO=INDEX(LINE,',')+1
      DO IT=1,IGRMX
      IF (LINE(:3).EQ.RESGRP(IT)(:3)) THEN
    2   IPCN=INDEX(LINE(IPCO:),',')+IPCO-1
        NGRES(IT)=NGRES(IT)+1
        IF (NGRES(IT).GT.30) THEN
          WRITE (*,*) 'error in solvent menv input file; res: ',line(:3)
          stop 321
          endif
        READ (LINE(IPCO:IPCN-1),1000) SMENV(NGRES(IT),IT)
        write(*,*)"within rdsenv, smenv(NGRES(IT),IT)",
     *  smenv(NGRES(IT),IT)
C1000   FORMAT(F<IPCN-IPCO>)
 1000   FORMAT(F6.3)
        IPCO=IPCN+1
        IF (LINE(IPCO:IPCO).EQ.'/') GO TO 1
        GO TO 2
        ENDIF
      ENDDO
  900 RETURN
      END
