C This subroutine is a part of original HADGEN code
C written by N.M.Sobolevsky
C changed by A.Timofeev at 19/07/2011
      SUBROUTINE INITIS
      COMMON /ABUNUC/ AIST(500),PERC(500),ABUNOR(500),
     *                NUCADR(2,110),IZMAX,L110
      CHARACTER*2 ELEM
C
C --------------------- INITIALIZATION -------------------
      IAD1=0
      IAD2=0
      L110=110
C      
      DO I=1,L110
        DO J=1,2
           NUCADR(J,I)=0
        END DO
      END DO
C
      DO I=1,500
        ABUNOR(I)=0.
      END DO
C ----------------------- DATA INPUT -------------------------
      OPEN(19,FILE='TABNUC.DAT',STATUS='OLD')
      
      DO J=1,L110+1
C Input of chemical symbol (ELEM) and atomic number (IZ) of the next
C nuclide and number (NISOT) of existing stable isitopes 
        READ(19,10)ELEM,IZ,NISOT
   10     FORMAT(A2,I3,1X,I7,60X)
C
        IF(ELEM.EQ.'ZZ' .OR. IZ.EQ.-1 .OR. NISOT.EQ.-1)THEN 
C         Maximal atomic number IZMAX, which is presented in TABNUC.DAT.
            IZMAX=J-1
            GO TO 100
C         End of file
        END IF
C
C Initial and last addresses of data for nuclide IZ 
C in arrays AIST, PERC, ABUNOR
        IAD1=IAD2+1
        IAD2=IAD1+NISOT-1
            NUCADR(1,J)=IAD1
            NUCADR(2,J)=IAD2
C Input of stable isotope mass numbers (AIST) and natural abundonces (PERC,%)       
        READ(19,20) (AIST(K),K=IAD1,IAD2)
   20     FORMAT(6X,15F8.4)
        READ(19,30) (PERC(K),K=IAD1,IAD2)
   30     FORMAT(6X,15F8.4)
C End of input cycle on J
      END DO
  100 CLOSE(19)
      RETURN
      END
c
      SUBROUTINE SETDKY() 
      IMPLICIT REAL*8 (A-H,O-Z) 
C          THIS SUBROUTINE READS IN THE DECAY TABLE FROM TAPE ITDKY 
C          AND SETS UP /DKYTAB/. 
C          QUARK-BASED IDENT CODE 
      COMMON /ATABDATA/ IRES(600),ITYPE(600),BR(600),IMODE(5,600),NLP
      LOOP=0 

      OPEN(18, FILE='ATAB.DAT', STATUS='OLD')
      
200   LOOP=LOOP+1 
      IF(LOOP.GT.600) GO TO 300
      READ(18,*,END=300) IRES22,ITYPE22,BR22,IMODE1,IMODE2,IMODE3
     *,IMODE4,IMODE5
      IRES(LOOP) = IRES22
      IF (IRES22.NE.0)THEN
         ITYPE(LOOP) = ITYPE22
         BR(LOOP) = BR22
         IMODE(1,LOOP) = IMODE1
         IMODE(2,LOOP) = IMODE2
         IMODE(3,LOOP) = IMODE3
         IMODE(4,LOOP) = IMODE4
         IMODE(5,LOOP) = IMODE5
         GO TO 200 
      ENDIF
300   CONTINUE 
      NLP = LOOP
      CLOSE(18)
      RETURN 
      END 
