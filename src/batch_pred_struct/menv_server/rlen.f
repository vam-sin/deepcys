C
C	return trimmed length of character aray
C
      INTEGER FUNCTION RLEN(CARR)
      CHARACTER CARR*(*),BL*1
      DATA BL/' '/
      N=LEN(CARR)
      DO I=1,N
      IF (CARR(N+1-I:N+1-I).NE.BL) THEN
        RLEN=N+1-I
        RETURN
        ENDIF
      ENDDO
      RLEN=0
      RETURN
      END
