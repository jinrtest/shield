C   These FORTRAN subroutines were written by A.Timofeev, 2012
      SUBROUTINE OPENSH
      OPEN(16,FILE='PASIN.DAT',STATUS='OLD')        ! GEMCA
!      OPEN(17,FILE='FOR017',STATUS='UNKNOWN')       ! GEMCA
      OPEN(22,FILE='FOR022.DAT',STATUS='OLD')       ! ERR
!      OPEN(25,FILE='FOR025',STATUS='UNKNOWN')           ! MEDIUM inmacr
      RETURN
      END
C
      SUBROUTINE CLOSSH
!      CLOSE(25)                                     ! GEMCA
!      CLOSE(17)                                     ! GEMCA
      CLOSE(16)                                     ! ERR
      CLOSE(22)                                     ! MEDIUM inmacr
      RETURN
      END

      SUBROUTINE INITLUX3(IX,LUX)
      CALL INIT_LUX(IX,LUX,24)
      RETURN
      END
