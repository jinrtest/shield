C     These subroutines are written by Alexander Timofeev, Dubna, 11.07.2011
C     Used in C code to open and close fortran output files
C     These subs are temportary, should be removed in future
      SUBROUTINE OPENFILES
!      OPEN(33,FILE='FOR033.DAT',STATUS='OLD')
!      OPEN(34,FILE='FOR034.DAT',STATUS='replace')
!      OPEN(25,FILE='FOR025.DAT',STATUS='replace')
c      OPEN(18,FILE='ATAB.DAT',STATUS='OLD')       ! NOT NEEDED IN HADGEN
!      OPEN(19,FILE='TABNUC.DAT',STATUS='OLD')      ! NOT NEEDED - implemented in C Code
!      OPEN(19,FILE='RANLUX.DAT', STATUS='replace')
!      OPEN(37,FILE='FOR037',FORM='BINARY')        ! HISTAGT
!      OPEN(37,FILE='FOR037.DAT',FORM='UNFORMATTED')        ! HISTAGT
c      OPEN(38,FILE='FOR038',STATUS='NEW')        ! HADES
!      IF(LSTAR.eq.1 .or. LSTAR.eq.2)
!     *OPEN(36,FILE='FOR036.DAT',STATUS='replace')            ! File for stars
C      IF(LSTAR.eq.1 .or. LSTAR.eq.2)CALL STITLE      ! Title for stars
      RETURN
      END
C
      SUBROUTINE CLOSEFILES
C      CLOSE(33)
!      CLOSE(34)
!      CLOSE(25)
!      CLOSE(18)
!      CLOSE(19)
C===================== HISTAGT ====================
!      CONTINUE!WRITE(37) (((SPECAGT(I,J,K), I=1,126), J=1,30), K=1,180)
!      IF(LSTAR.eq.1 .or. LSTAR.eq.2)CLOSE(36)
!      CLOSE(37)     
      RETURN
      END
C
      SUBROUTINE INITLUX2(IXFIRS,LUXCNT)
      CALL INIT_LUX(IXFIRS, LUXCNT, 34)
      RETURN
      END
