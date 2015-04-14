**********************************************************************
*                                                                    *
*        The following subroutines are written by A.V.Dementyev,     *
*        Institute for Nuclear Research RAS, Moscow,                 *
*        E-mail: AVDEM@AL20.INR.Troitsk.RU                           *
*                                                                    *
*        Corrections of February 1999 by N.M.Sobolevsky              *
*        (see the label "Sobol")                                     *
*        E-mail: Sobolevs@AL20.INR.Troitsk.RU                        *
*                                                                    *
*        Additional corrections of 16.04.99 by Sobolevsky            *
**********************************************************************

      BLOCK DATA GBLCD
      IMPLICIT REAL*8 (A-H,R-Z)
C
      PARAMETER(MAXB=1000,MAXZ=1000)
      PARAMETER(JDIMBD=37*MAXB,JDIMZD=5*MAXZ)
      PARAMETER(MKINDB=10,MKINDS=10)
C----------------------------------------------------------
C MAXB   - MAX NUMBER OF BODYS
C MAXZ   - MAX NUMBER OF ZONES
C JDIMBD - DIMENTION FOR BODY DATA
C JDIMZD - DIMENTION FOR ZONE DATA
C----------------------------------------------------------
C
      COMMON/GDATA0/ NUMBOD,NUMZON
      COMMON/GDATA1/ ISYMB(MKINDB),NBSURF(MKINDB),NSPAR(MKINDS)
      COMMON/GDATA2/ BODYDA(MAXB),BODYDB(JDIMBD)
      INTEGER BODYDA
      COMMON/GDATA3/ ZONEDA(MAXZ),ZONEDB(JDIMZD)
      INTEGER ZONEDA,ZONEDB
      COMMON/GDATA4/ JFLSUF(3,MAXB),BTRZLT(MAXB),
     *                 DIGINF,DIGMIN,DELTA,STPLIM
      COMMON/GDATAP/ X,Y,Z,CX,CY,CZ,XN,YN,ZN,CXN,CYN,CZN,
     *                 NZONO,MEDOLD,NBPO,NSO,
     *                 NZONC,MEDCUR,NBPC,NSCI,NSCO,
     *                 PINSFL,PARTZD
      INTEGER PINSFL,PARTZD
C
      DATA  ISYMB/ 3HSPH, 3HWED, 3HARB,
     *                 3HBOX, 3HRPP, 3HRCC,
     *                 3HREC, 3HTRC, 3HELL,
     *                 3HEND/
C
C-------------------------------------------------------------
C               The number of body surfaces :                |
C-------------------------------------------------------------
C                     1   2   3   4   5   6   7   8   9   10 |
C BODY               SPH WED ARB BOX RPP RCC REC TRC ELL END |
C-------------------------------------------------------------
      DATA NBSURF/  1,  5,  6,  6,  6,  3,  3,  3,  1,  0/
C-------------------------------------------------------------
C               The number of surface parameters :           |
C-------------------------------------------------------------
C                     1   2   3   4   5   6   7   8   9   10 |
C SURFACE                                                    |
C-------------------------------------------------------------
      DATA NSPAR/   1,  1,  1,  4,  4,  1,  2,  2,  3,  0/
C-------------------------------------------------------------
      DATA DIGINF/1.0E+10/
        DATA DIGMIN/0.0000000000001/
      DATA DELTA/0.000001/
      DATA STPLIM/0.0001/
C-------------------------------------------------------------
      END


      SUBROUTINE GNEXTZ(NEXTZ)
C****************************************************************
C*  Next zone number determination                              *
C****************************************************************
	IMPLICIT REAL*8 (A-H,R-Z)
C
	PARAMETER(MAXB=1000,MAXZ=1000)
	PARAMETER(JDIMBD=37*MAXB,JDIMZD=5*MAXZ)
	PARAMETER(MKINDB=10,MKINDS=10)
C
	COMMON/GDATA0/ NUMBOD,NUMZON
	COMMON/GDATA1/ ISYMB(MKINDB),NBSURF(MKINDB),NSPAR(MKINDS)
	COMMON/GDATA2/ BODYDA(MAXB),BODYDB(JDIMBD)
	INTEGER BODYDA
	COMMON/GDATA3/ ZONEDA(MAXZ),ZONEDB(JDIMZD)
	INTEGER ZONEDA,ZONEDB
	COMMON/GDATA4/ JFLSUF(3,MAXB),BTRZLT(MAXB),
     *                 DIGINF,DIGMIN,DELTA,STPLIM
	COMMON/GDATAP/ X,Y,Z,CX,CY,CZ,XN,YN,ZN,CXN,CYN,CZN,
     *                 NZONO,MEDOLD,NBPO,NSO,
     *                 NZONC,MEDCUR,NBPC,NSCI,NSCO,
     *                 PINSFL,PARTZD
	INTEGER PINSFL,PARTZD
C
      PINSFL=10
C
      DO I=1,NUMBOD
        CALL GTBODN(I)
      END DO
C
C   Searching next zone:
C
      DO I=1,NUMZON
        PRIZNZ=1
        PARTZD=1
        JZADDR=ZONEDA(I)
        NZB=ZONEDB(JZADDR)
C
        IF(ZONEDB(JZADDR+2).EQ.0)THEN
          GOTO 2
C       ELSE
C         GOTO 1
        ENDIF
C
    1   DO J=1,NZB
          KINDB=ZONEDB(JZADDR+2+J)
          IF(KINDB*JFLSUF(1,ABS(KINDB)).LT.0)THEN
C           PRIZNZ=0
            GOTO 3
          ENDIF
        END DO
      GOTO 4
C
    2   PRIZNZ=0
        PARTZD=1
        DO J=1,NZB
          KINDB=ZONEDB(JZADDR+2+J)
          IF(KINDB.EQ.0)THEN
            IF(PRIZNZ.EQ.1)THEN
              GOTO 4
            ELSE
              PRIZNZ=1
              PARTZD=PARTZD+1
	    ENDIF
          ELSE
            IF(KINDB*JFLSUF(1,ABS(KINDB)).LT.0)PRIZNZ=0
          ENDIF
        END DO
C
        IF(PRIZNZ.EQ.1)THEN
          GOTO 4
        ELSE
          GOTO 3
        ENDIF
    3 END DO
C
C------------ DEBUG -----------
C       DO I=1,NUMBOD
C        IF(JFLSUF(1,I).GT.0)THEN
C          CALL GCLCER
C          STOP'EXIT from GNEXTZ: undefined zone !'
C        ENDIF
C       ENDDO
C------------------------------
C
C Particle go out from configuration:
C
      NEXTZ=0
      NZONO=NZONC
      NZONC=NEXTZ
      MEDOLD=MEDCUR
      MEDCUR=0
      RETURN
C
C Next zone is found:
C
    4 NEXTZ=I
      NZONO=NZONC
      NZONC=NEXTZ
      MEDOLD=MEDCUR
      MEDCUR=ZONEDB(JZADDR+1)
      NBPO=NBPC
      NSO=NSCO
C
      RETURN
      END


	SUBROUTINE GTBODN(NB)
C****************************************************************
C* Test of body (belong only)                                   *
C****************************************************************
	IMPLICIT REAL*8 (A-H,R-Z)
C
	PARAMETER(MAXB=1000,MAXZ=1000)
	PARAMETER(JDIMBD=37*MAXB,JDIMZD=5*MAXZ)
	PARAMETER(MKINDB=10,MKINDS=10)
C
	COMMON/GDATA0/ NUMBOD,NUMZON
	COMMON/GDATA1/ ISYMB(MKINDB),NBSURF(MKINDB),NSPAR(MKINDS)
	COMMON/GDATA2/ BODYDA(MAXB),BODYDB(JDIMBD)
	INTEGER BODYDA
	COMMON/GDATA3/ ZONEDA(MAXZ),ZONEDB(JDIMZD)
	INTEGER ZONEDA,ZONEDB
	COMMON/GDATA4/ JFLSUF(3,MAXB),BTRZLT(MAXB),
     *                 DIGINF,DIGMIN,DELTA,STPLIM
	COMMON/GDATAP/ X,Y,Z,CX,CY,CZ,XN,YN,ZN,CXN,CYN,CZN,
     *                 NZONO,MEDOLD,NBPO,NSO,
     *                 NZONC,MEDCUR,NBPC,NSCI,NSCO,
     *                 PINSFL,PARTZD
	INTEGER PINSFL,PARTZD
C
	IBADDR=BODYDA(NB)
	ITYPEB=BODYDB(IBADDR)
	NSURF=NBSURF(ITYPEB)
C
C--------------------------------------------
	 IF(ITYPEB.GT.5)THEN
C----------------------
C   going to new system
C
	JP=IBADDR+1
C
	XN=X*BODYDB(JP)+Y*BODYDB(JP+1)+Z*BODYDB(JP+2)-BODYDB(JP+3)
	YN=X*BODYDB(JP+4)+Y*BODYDB(JP+5)+Z*BODYDB(JP+6)-BODYDB(JP+7)
	ZN=X*BODYDB(JP+8)+Y*BODYDB(JP+9)+Z*BODYDB(JP+10)-BODYDB(JP+11)
C
	CXN=CX*BODYDB(JP)+CY*BODYDB(JP+1)+CZ*BODYDB(JP+2)
	CYN=CX*BODYDB(JP+4)+CY*BODYDB(JP+5)+CZ*BODYDB(JP+6)
	CZN=CX*BODYDB(JP+8)+CY*BODYDB(JP+9)+CZ*BODYDB(JP+10)
C
	IBADDR=IBADDR+12
C
	 ELSE
C
	XN=X
	YN=Y
	ZN=Z
	CXN=CX
	CYN=CY
	CZN=CZ
C
	 END IF
C--------------------------------------------
C
	ISTAT=PINSFL
C
	ITSURF=BODYDB(IBADDR+1)
	JPADDR=IBADDR+2
C
	     DO I=1,NSURF
	CALL GSURFN(ITSURF,JPADDR,JSFLAG,NB,I)
C
	  IF(JSFLAG.LT.0)THEN
	JFLSUF(1,NB)=-1
	RETURN
	  END IF
C
	IF(JSFLAG.EQ.0)ISTAT=-10
C
	JPADDR=JPADDR+NSPAR(ITSURF)+2
	ITSURF=BODYDB(JPADDR-1)
C
	    ENDDO
C--------------------------------------------
C
	PINSFL=ISTAT
	JFLSUF(1,NB)=1
C
	RETURN
	END



	SUBROUTINE GSURFN(ITSURF,JPADDR,JSFLAG,NB,NS)
C****************************************************************
C* Test of body surface (belong to body only)                   *
C****************************************************************
	IMPLICIT REAL*8 (A-H,R-Z)
C
	PARAMETER(MAXB=1000,MAXZ=1000)
	PARAMETER(JDIMBD=37*MAXB,JDIMZD=5*MAXZ)
	PARAMETER(MKINDB=10,MKINDS=10)
C
	COMMON/GDATA0/ NUMBOD,NUMZON
	COMMON/GDATA1/ ISYMB(MKINDB),NBSURF(MKINDB),NSPAR(MKINDS)
	COMMON/GDATA2/ BODYDA(MAXB),BODYDB(JDIMBD)
	INTEGER BODYDA
	COMMON/GDATA3/ ZONEDA(MAXZ),ZONEDB(JDIMZD)
	INTEGER ZONEDA,ZONEDB
	COMMON/GDATA4/ JFLSUF(3,MAXB),BTRZLT(MAXB),
     *                 DIGINF,DIGMIN,DELTA,STPLIM
	COMMON/GDATAP/ X,Y,Z,CX,CY,CZ,XN,YN,ZN,CXN,CYN,CZN,
     *                 NZONO,MEDOLD,NBPO,NSO,
     *                 NZONC,MEDCUR,NBPC,NSCI,NSCO,
     *                 PINSFL,PARTZD
	INTEGER PINSFL,PARTZD
C
	SIGN=GSSIGN(ITSURF,JPADDR,XN,YN,ZN)
C
	IF(SIGN.LT.-DIGMIN)THEN
1           JSFLAG=-1
	    RETURN
	END IF
C
	    IF(SIGN.LE.DIGMIN)THEN
	SIGN=GSSIGN(ITSURF,JPADDR,XN+DELTA*CXN,YN+DELTA*CYN,
     *                            ZN+DELTA*CZN)
	 IF(SIGN.LT.-DIGMIN)THEN
	  JSFLAG=-1
	  RETURN
	 ENDIF
C
		   IF(SIGN.LE.DIGMIN)THEN
           IF(PINSFL.LT.0)GOTO 1   ! Sobol: GOTO into IF-ENDIF
         JSFLAG=0                  ! It is not very well
		   ELSE
		 JFLSUF(2,NB)=NS
		 JSFLAG=1
		   END IF
	    ELSE
		 JSFLAG=1
	    END IF
C
	    RETURN
	END


	SUBROUTINE GCURZL(JZONE,S)
C****************************************************************
C* Determination of length to the zone boundary                 *
C****************************************************************
	IMPLICIT REAL*8 (A-H,R-Z)
C
	PARAMETER(MAXB=1000,MAXZ=1000)
	PARAMETER(JDIMBD=37*MAXB,JDIMZD=5*MAXZ)
	PARAMETER(MKINDB=10,MKINDS=10)
C
	COMMON/GDATA0/ NUMBOD,NUMZON
	COMMON/GDATA1/ ISYMB(MKINDB),NBSURF(MKINDB),NSPAR(MKINDS)
	COMMON/GDATA2/ BODYDA(MAXB),BODYDB(JDIMBD)
	INTEGER BODYDA
	COMMON/GDATA3/ ZONEDA(MAXZ),ZONEDB(JDIMZD)
	INTEGER ZONEDA,ZONEDB
	COMMON/GDATA4/ JFLSUF(3,MAXB),BTRZLT(MAXB),
     *                 DIGINF,DIGMIN,DELTA,STPLIM
	COMMON/GDATAP/ X,Y,Z,CX,CY,CZ,XN,YN,ZN,CXN,CYN,CZN,
     *                 NZONO,MEDOLD,NBPO,NSO,
     *                 NZONC,MEDCUR,NBPC,NSCI,NSCO,
     *                 PINSFL,PARTZD
	INTEGER PINSFL,PARTZD
C------------------------------------------------------
C Maybe it is need to check JZONE on zero (=0) or < 0 |
C------------------------------------------------------
C
	PINSFL=10
C
	NZONC=JZONE
	JZADDR=ZONEDA(JZONE)
	NZBODY=ZONEDB(JZADDR)
	MEDCUR=ZONEDB(JZADDR+1)
C
	IF(ZONEDB(JZADDR+2).EQ.-1)THEN
	 GOTO 3
	ELSEIF(ZONEDB(JZADDR+2).EQ.0)THEN
	 GOTO 2
C       ELSE
C        GOTO 1
	ENDIF
C
C-----------------------------------
C only bodys with "+" sign are in the zone definition
C
1              DO I=1,NZBODY
	  NB=ZONEDB(JZADDR+2+I)
C
	  CALL GTBODY(NB)
	       END DO
	  GOTO 4
C
C-----------------------------------
C "OR" is in the zone definition
C
2       PARTCZ=1
	       DO I=1,NZBODY
	  NB=ZONEDB(JZADDR+2+I)
	IF(PARTCZ.EQ.PARTZD)THEN
C
	 IF(NB.GT.0)THEN
	   CALL GTBODY(NB)
	 ELSEIF(NB.LT.0)THEN
	   CALL GTOUTB(ABS(NB))
	 ELSE
	   PARTCZ=PARTCZ+1
	 END IF
C
	ELSE
C
	 IF(NB.EQ.0)THEN
	   PARTCZ=PARTCZ+1
	 ELSE
	   JFLSUF(1,ABS(NB))=-1
	 END IF
C
	END IF
	       END DO
	  GOTO 4
C
C-----------------------------------
C body with "-" sign is in the zone definition
C
3               DO I=1,NZBODY
	  NB=ZONEDB(JZADDR+2+I)
C
	IF(NB.GT.0)THEN
	  CALL GTBODY(NB)
	ELSE
	  CALL GTOUTB(ABS(NB))
	END IF
	       END DO
C         GOTO 4
C
C-----------------------------------
4         S=DIGINF
C
	DO I=1,NZBODY
	  NB=ZONEDB(JZADDR+2+I)
	  NBA=ABS(NB)
	IF(JFLSUF(1,NBA).LT.0 .OR. NBA.EQ.0)GOTO 5
	  T=BTRZLT(NBA)
		      IF(T.LT.S)THEN
			S=T
			NBPC=NBA
			NSCI=JFLSUF(2,NBA)
			NSCO=JFLSUF(3,NBA)
		      END IF
C
5       END DO
C-------------------------------
C Step limitation
C
	IF(S.LT.STPLIM)S=STPLIM
C
	IF(S.LT.DIGINF)THEN
	 RETURN
	ELSE
Cdeb:             STOP'EXIT from GCURZL : particle is not in this zone !'
	 S=0.
	 RETURN
	END IF
C
	END



	SUBROUTINE GTBODY(NB)
C****************************************************************
C* Test of body - point inside the body                         *
C****************************************************************
	IMPLICIT REAL*8 (A-H,R-Z)
C
	PARAMETER(MAXB=1000,MAXZ=1000)
	PARAMETER(JDIMBD=37*MAXB,JDIMZD=5*MAXZ)
	PARAMETER(MKINDB=10,MKINDS=10)
C
	COMMON/GDATA0/ NUMBOD,NUMZON
	COMMON/GDATA1/ ISYMB(MKINDB),NBSURF(MKINDB),NSPAR(MKINDS)
	COMMON/GDATA2/ BODYDA(MAXB),BODYDB(JDIMBD)
	INTEGER BODYDA
	COMMON/GDATA3/ ZONEDA(MAXZ),ZONEDB(JDIMZD)
	INTEGER ZONEDA,ZONEDB
	COMMON/GDATA4/ JFLSUF(3,MAXB),BTRZLT(MAXB),
     *                 DIGINF,DIGMIN,DELTA,STPLIM
	COMMON/GDATAP/ X,Y,Z,CX,CY,CZ,XN,YN,ZN,CXN,CYN,CZN,
     *                 NZONO,MEDOLD,NBPO,NSO,
     *                 NZONC,MEDCUR,NBPC,NSCI,NSCO,
     *                 PINSFL,PARTZD
	INTEGER PINSFL,PARTZD
C
	IBADDR=BODYDA(NB)
	ITYPEB=BODYDB(IBADDR)
	NSURF=NBSURF(ITYPEB)
C
	 IF(ITYPEB.GT.5)THEN
C------------------
C   going to new system
C
	JP=IBADDR+1
C
	XN=X*BODYDB(JP)+Y*BODYDB(JP+1)+Z*BODYDB(JP+2)-BODYDB(JP+3)
	YN=X*BODYDB(JP+4)+Y*BODYDB(JP+5)+Z*BODYDB(JP+6)-BODYDB(JP+7)
	ZN=X*BODYDB(JP+8)+Y*BODYDB(JP+9)+Z*BODYDB(JP+10)-BODYDB(JP+11)
C
	CXN=CX*BODYDB(JP)+CY*BODYDB(JP+1)+CZ*BODYDB(JP+2)
	CYN=CX*BODYDB(JP+4)+CY*BODYDB(JP+5)+CZ*BODYDB(JP+6)
	CZN=CX*BODYDB(JP+8)+CY*BODYDB(JP+9)+CZ*BODYDB(JP+10)
C
	IBADDR=IBADDR+12
C
	 ELSE
C
	XN=X
	YN=Y
	ZN=Z
	CXN=CX
	CYN=CY
	CZN=CZ
C
	 END IF
C
	ISTAT=PINSFL
C
	ITSURF=BODYDB(IBADDR+1)
	JPADDR=IBADDR+2
C
	T=DIGINF
C-------------------------------------------------------
	     DO I=1,NSURF
	CALL GTSURF(ITSURF,JPADDR,10,JSFLAG,TS1)
C
	  IF(JSFLAG.LT.0)THEN
	JFLSUF(1,NB)=-1
	RETURN
	  ELSE
C
	IF(JSFLAG.EQ.0)ISTAT=-1
C
C   IF(TS1.GT.0. .AND. TS1.LT.T)THEN      ! Error of AVDEM. TS1 can
C                                         ! be equal zero EXACTLY!!!
        IF(TS1.ge.0. .AND. TS1.LT.T)THEN  ! Sobol, 02.99
	 T=TS1
c original version  22.01.97
c        NS=I
c for full CG compatibility
	 NS=-I
	END IF
C
	  END IF
C
	JPADDR=JPADDR+NSPAR(ITSURF)+2
	ITSURF=BODYDB(JPADDR-1)
C
	     END DO
C------------------------------------------------------
	PINSFL=ISTAT
	JFLSUF(1,NB)=1
	JFLSUF(3,NB)=NS
	BTRZLT(NB)=T
C
	RETURN
	END



	SUBROUTINE GTOUTB(NB)
C****************************************************************
C* Test of body - point outside the body                        *
C****************************************************************
	IMPLICIT REAL*8 (A-H,R-Z)
C
	PARAMETER(MAXB=1000,MAXZ=1000)
	PARAMETER(JDIMBD=37*MAXB,JDIMZD=5*MAXZ)
	PARAMETER(MKINDB=10,MKINDS=10)
C
	COMMON/GDATA0/ NUMBOD,NUMZON
	COMMON/GDATA1/ ISYMB(MKINDB),NBSURF(MKINDB),NSPAR(MKINDS)
	COMMON/GDATA2/ BODYDA(MAXB),BODYDB(JDIMBD)
	INTEGER BODYDA
	COMMON/GDATA3/ ZONEDA(MAXZ),ZONEDB(JDIMZD)
	INTEGER ZONEDA,ZONEDB
	COMMON/GDATA4/ JFLSUF(3,MAXB),BTRZLT(MAXB),
     *                 DIGINF,DIGMIN,DELTA,STPLIM
	COMMON/GDATAP/ X,Y,Z,CX,CY,CZ,XN,YN,ZN,CXN,CYN,CZN,
     *                 NZONO,MEDOLD,NBPO,NSO,
     *                 NZONC,MEDCUR,NBPC,NSCI,NSCO,
     *                 PINSFL,PARTZD
	INTEGER PINSFL,PARTZD
C
	IBADDR=BODYDA(NB)
	ITYPEB=BODYDB(IBADDR)
	NSURF=NBSURF(ITYPEB)
C
	 IF(ITYPEB.GT.5)THEN
	JP=IBADDR+1
C
	XN=X*BODYDB(JP)+Y*BODYDB(JP+1)+Z*BODYDB(JP+2)-BODYDB(JP+3)
	YN=X*BODYDB(JP+4)+Y*BODYDB(JP+5)+Z*BODYDB(JP+6)-BODYDB(JP+7)
	ZN=X*BODYDB(JP+8)+Y*BODYDB(JP+9)+Z*BODYDB(JP+10)-BODYDB(JP+11)
C
	CXN=CX*BODYDB(JP)+CY*BODYDB(JP+1)+CZ*BODYDB(JP+2)
	CYN=CX*BODYDB(JP+4)+CY*BODYDB(JP+5)+CZ*BODYDB(JP+6)
	CZN=CX*BODYDB(JP+8)+CY*BODYDB(JP+9)+CZ*BODYDB(JP+10)
C
	IBADDR=IBADDR+12
C
	 ELSE
C
	XN=X
	YN=Y
	ZN=Z
	CXN=CX
	CYN=CY
	CZN=CZ
C
	 END IF
C
	ITSURF=BODYDB(IBADDR+1)
	JPADDR=IBADDR+2
C
	T=DIGINF
      DO I=1,NSURF                         !----------begin loop 1
	CALL GTSURF(ITSURF,JPADDR,-10,JSFLAG,TS1)
C       IF(TS1.GT.0.0 .AND. TS1.LT.STPLIM)TS1=STPLIM   ! AVDem 20.07.1998
C                                                      ! Error of AVDEM (?)
        IF(TS1.ge.0.0 .AND. TS1.LT.STPLIM)TS1=STPLIM   ! Sobol, 02.99
C
        IF(TS1.GT.DIGMIN .AND. TS1.LT.DIGINF)THEN
          XS=XN+TS1*CXN
          YS=YN+TS1*CYN
          ZS=ZN+TS1*CZN
          KSURF=BODYDB(IBADDR+1)
          KADDR=IBADDR+2
C
            DO K=1,NSURF                   !---------- begin loop 2
C                IF(ITSURF.EQ.KSURF)GOTO 1   ! Error of AVDEM; C - Sobol !!!
              SIGN=GSSIGN(KSURF,KADDR,XS,YS,ZS)
C                IF(SIGN.LT.-DIGMIN)GOTO 2  ! Sobol 16.04.99
                IF(SIGN.LT.-0.000001)GOTO 2  ! This is very sensitive point,
C                                           ! especially for SPH !
    1         KADDR=KADDR+NSPAR(KSURF)+2
              KSURF=BODYDB(KADDR-1)
            END DO                          !---------- end loop 2
C
C
            IF(T.GT.TS1)THEN
              T=TS1
c             NS=-I  ! original version  22.01.97
              NS=I   ! for full CG compatibility
            END IF
	  END IF
C
    2   JPADDR=JPADDR+NSPAR(ITSURF)+2
	ITSURF=BODYDB(JPADDR-1)
      END DO                                 !---------- end loop 1
C
	   IF(T.LT.DIGINF)THEN
	JFLSUF(1,NB)=1
	JFLSUF(3,NB)=NS
	BTRZLT(NB)=T
	   ELSE
	JFLSUF(1,NB)=-1
	   END IF
C
	RETURN
	END


	SUBROUTINE GTSURF(ITSURF,JPADDR,JBFLAG,JSFLAG,TS1)
C****************************************************************
C* Test of body surface                                         *
C****************************************************************
	IMPLICIT REAL*8 (A-H,R-Z)
C
	PARAMETER(MAXB=1000,MAXZ=1000)
	PARAMETER(JDIMBD=37*MAXB,JDIMZD=5*MAXZ)
	PARAMETER(MKINDB=10,MKINDS=10)
C
	COMMON/GDATA0/ NUMBOD,NUMZON
	COMMON/GDATA1/ ISYMB(MKINDB),NBSURF(MKINDB),NSPAR(MKINDS)
	COMMON/GDATA2/ BODYDA(MAXB),BODYDB(JDIMBD)
	INTEGER BODYDA
	COMMON/GDATA3/ ZONEDA(MAXZ),ZONEDB(JDIMZD)
	INTEGER ZONEDA,ZONEDB
	COMMON/GDATA4/ JFLSUF(3,MAXB),BTRZLT(MAXB),
     *                 DIGINF,DIGMIN,DELTA,STPLIM
	COMMON/GDATAP/ X,Y,Z,CX,CY,CZ,XN,YN,ZN,CXN,CYN,CZN,
     *                 NZONO,MEDOLD,NBPO,NSO,
     *                 NZONC,MEDCUR,NBPC,NSCI,NSCO,
     *                 PINSFL,PARTZD
	INTEGER PINSFL,PARTZD,SIGNFL
C
	     SIGNFL=10
	IF(JBFLAG.LT.0)GOTO 102
C
100     SIGN=GSSIGN(ITSURF,JPADDR,XN,YN,ZN)   ! call from GTBODY
C
	IF(SIGN.LT.-DIGMIN)THEN
101         JSFLAG=-1
	    RETURN
	END IF
C
	IF(SIGN.LE.DIGMIN)THEN
	  IF(PINSFL.LT.0)GOTO 101  ! Sobol: GOTO into IF-ENDIF
	    IF(SIGNFL.LT.0)THEN    ! it is not very well.
	       XN=XN-DELTA*CXN
	       YN=YN-DELTA*CYN
	       ZN=ZN-DELTA*CZN
C
	       JSFLAG=0
	       TS1=DIGINF
	       RETURN
	    ELSE
	       XN=XN+DELTA*CXN
	       YN=YN+DELTA*CYN
	       ZN=ZN+DELTA*CZN
	      SIGNFL=-10
	       GOTO 100
	    END IF
	END IF
C
  102  CONTINUE                               ! call from GTOUTB
            GOTO ( 1, 2, 3, 4, 5, 6, 7, 8, 9, 10),ITSURF
C
1       IF(CXN.EQ.0.)GOTO 104
	 C=BODYDB(JPADDR)
	 A=C-XN
	 TS1=A/CXN
	GOTO 105
C
2       IF(CYN.EQ.0.)GOTO 104
	 C=BODYDB(JPADDR)
	 A=C-YN
	 TS1=A/CYN
	GOTO 105
C
    3 CONTINUE
        IF(CZN.EQ.0.)GOTO 104
	 C=BODYDB(JPADDR)
	 A=C-ZN
	 TS1=A/CZN
	GOTO 105
C
4        A=BODYDB(JPADDR)
	 B=BODYDB(JPADDR+1)
	 C=BODYDB(JPADDR+2)
	 D=BODYDB(JPADDR+3)
	ZNUM=A*CXN+B*CYN+C*CZN
	IF(ZNUM.EQ.0.)GOTO 104
	CHIS=D-A*XN-B*YN-C*ZN
	 TS1=CHIS/ZNUM
	GOTO 105
C
5        X0=XN-BODYDB(JPADDR)
	 Y0=YN-BODYDB(JPADDR+1)
	 Z0=ZN-BODYDB(JPADDR+2)
	 R2=BODYDB(JPADDR+3)
C       A=CXN**2.+CYN**2.+CZN**2.
	A=1.
	B=CXN*X0+CYN*Y0+CZN*Z0
	C=X0**2.+Y0**2.+Z0**2.-R2
	  CALL GSQRTS(A,B,C,NR,TS1,TS2)
	IF(NR.LE.1)GOTO 103
	GOTO 105
C
6        R2=BODYDB(JPADDR)
	A=CXN**2.+CYN**2.
	B=XN*CXN+YN*CYN
	C=XN**2.+YN**2.-R2
	  CALL GSQRTS(A,B,C,NR,TS1,TS2)
	IF(NR.LE.1)GOTO 103
	GOTO 105
C
7        A2=BODYDB(JPADDR)
	 B2=BODYDB(JPADDR+1)
	A=CXN**2./A2+CYN**2./B2
	B=XN*CXN/A2+YN*CYN/B2
	C=XN**2./A2+YN**2./B2-1.
	  CALL GSQRTS(A,B,C,NR,TS1,TS2)
	IF(NR.LE.1)GOTO 103
	GOTO 105
C
    8 CONTINUE
         A1=BODYDB(JPADDR)
	  B2=BODYDB(JPADDR+1)
	A=CXN**2.+CYN**2.-CZN**2./B2
	B=XN*CXN+YN*CYN-(ZN-A1)*CZN/B2
	C=XN**2.+YN**2.-(ZN-A1)**2./B2
	  CALL GSQRTS(A,B,C,NR,TS1,TS2)
	IF(NR.EQ.0)GOTO 103
	GOTO 105
C
9        A2=BODYDB(JPADDR)
	 B2=BODYDB(JPADDR+1)
	 C2=BODYDB(JPADDR+2)
	A=CXN**2./A2+CYN**2./B2+CZN**2./C2
	B=XN*CXN/A2+YN*CYN/B2+ZN*CZN/C2
	C=XN**2./A2+YN**2./B2+ZN**2./C2-1.
	  CALL GSQRTS(A,B,C,NR,TS1,TS2)
	IF(NR.LE.1)GOTO 103
	GOTO 105
C
10      PRINT *,'Surface=',ITSURF
	STOP'EXIT from GTSURF : no such surface !'
C
103     JSFLAG=-1
	TS1=DIGINF
	RETURN
C
104     TS1=DIGINF
C
105     JSFLAG=1
	    IF(SIGNFL.LT.0)THEN
	       XN=XN-DELTA*CXN
	       YN=YN-DELTA*CYN
	       ZN=ZN-DELTA*CZN
	 TS1=TS1+DELTA
	IF(ITSURF.LE.4)TS1=DIGINF
	    END IF
C
	RETURN
	END


	SUBROUTINE GSQRTS(A,B,C,NR,R1,R2)
C****************************************************************
C* Square equation solution                                     *
C****************************************************************
	IMPLICIT REAL*8 (A-H,R-Z)
C
	PARAMETER(MAXB=1000,MAXZ=1000)
C
	COMMON/GDATA4/ JFLSUF(3,MAXB),BTRZLT(MAXB),
     *                 DIGINF,DIGMIN,DELTA,STPLIM
C
        if(A.lt.0.)then     !
          A=-A              !    Error of AVDEM !
          B=-B              ! Added by Sobol 02.99
          C=-C              !
        end if              !
C
	IF(A.EQ.0.)THEN
	 R1=DIGINF
	 NR=2
	RETURN
	ENDIF
C
	D=B**2.-A*C        ! AVDEM
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*   Note of Sobol: If B is about 1, while A and C are very small,     *
*                  such that A*C < 1.0E-14, then D=B**2 EXACTLY!!!    *
*      13.02.99    In this case R1=0.0 EXACTLY!!!, while R2 is some   *
*                  negative value.                                    *
*                  Therefore the case R1=0.0 must be forseen at all   *
*                  subroutines, in particular in GTBODY (see comment) *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
	IF(D.LT.-DIGMIN)THEN   ! -DIGMIN ?
	 NR=0
	 RETURN
	END IF
C
	IF(D.GT.DIGMIN)THEN
	 NR=2
	 R1=(-B-DSQRT(D))/A
	    IF(R1.LT.0.)THEN
	 R1=(-B+DSQRT(D))/A
	    END IF
	ELSE
	 NR=1
	 R1=-B/A
	END IF
C
	RETURN
	END




	REAL*8 FUNCTION GSSIGN(KSURF,KADDR,X,Y,Z)
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
* There is a question. The variable SSIGN gives a measure of deviation  *
* of the point (X,Y,Z) from the surface KSURF. The question is:         *
* Why for some surfaces (namely for KSURF=1,2,3,4,5,6,8) SSIGN is the   *
* absolute deviation while for others (KSURF=7,9) SSIGN is the relative *
* deviation?                                                            *
* Examples:                                                             *
* For sphere (KSURF=5):  SSIGN=X0**2.+Y0**2.+Z0**2.-R2                  *
* For ellipsoid (KSURF=9):  SSIGN=X**2./A2+Y**2./B2+Z**2./C2-1.         *
* As a rule the relative deviation is much more less (for several       *
* orders) then the absolute one. (One can compare numerically SPH       *
* and ELL).                                                             *
*                                    Sobol 16.04.99                     *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
	IMPLICIT REAL*8 (A-H,R-Z)
C
	PARAMETER(MAXB=1000,MAXZ=1000)
	PARAMETER(JDIMBD=37*MAXB,JDIMZD=5*MAXZ)
	PARAMETER(MKINDB=10,MKINDS=10)
C
	COMMON/GDATA2/ BODYDA(MAXB),BODYDB(JDIMBD)
	INTEGER BODYDA
C
	    GOTO ( 1, 2, 3, 4, 5, 6, 7, 8, 9, 10),KSURF
C
C
1        C=BODYDB(KADDR)
	 TRUES=BODYDB(KADDR+1)
	SSIGN=X-C
	GSSIGN=SSIGN*TRUES
	RETURN
C
2        C=BODYDB(KADDR)
	 TRUES=BODYDB(KADDR+1)
	SSIGN=Y-C
	GSSIGN=SSIGN*TRUES
	RETURN
C
3        C=BODYDB(KADDR)
	 TRUES=BODYDB(KADDR+1)
	SSIGN=Z-C
	GSSIGN=SSIGN*TRUES
	RETURN
C
    4 CONTINUE
         A=BODYDB(KADDR)
	 B=BODYDB(KADDR+1)
	 C=BODYDB(KADDR+2)
	 D=BODYDB(KADDR+3)
	 TRUES=BODYDB(KADDR+4)
	SSIGN=A*X+B*Y+C*Z-D
	GSSIGN=SSIGN*TRUES
	RETURN
C
5        X0=X-BODYDB(KADDR)
	 Y0=Y-BODYDB(KADDR+1)
	 Z0=Z-BODYDB(KADDR+2)
	 R2=BODYDB(KADDR+3)
	 TRUES=BODYDB(KADDR+4)
C        SSIGN=X0**2.+Y0**2.+Z0**2.-R2
      SSIGN=X0**2./R2+Y0**2./R2+Z0**2./R2-1.   ! Sobol 16.04.99
        GSSIGN=SSIGN*TRUES
	RETURN
C
6        R2=BODYDB(KADDR)
	 TRUES=BODYDB(KADDR+1)
	SSIGN=X**2.+Y**2.-R2
	GSSIGN=SSIGN*TRUES
	RETURN
C
7        A2=BODYDB(KADDR)
	 B2=BODYDB(KADDR+1)
	 TRUES=BODYDB(KADDR+2)
	SSIGN=X**2./A2+Y**2./B2-1.
	GSSIGN=SSIGN*TRUES
	RETURN
C
    8 CONTINUE
          A1=BODYDB(KADDR)
	  B2=BODYDB(KADDR+1)
	 TRUES=BODYDB(KADDR+2)
        SSIGN=X**2.+Y**2.-(Z-A1)**2./B2
	GSSIGN=SSIGN*TRUES
	RETURN
C
9        A2=BODYDB(KADDR)
	 B2=BODYDB(KADDR+1)
	 C2=BODYDB(KADDR+2)
	 TRUES=BODYDB(KADDR+3)
	SSIGN=X**2./A2+Y**2./B2+Z**2./C2-1.
	GSSIGN=SSIGN*TRUES
	RETURN
C
10      PRINT *,'Surface=',KSURF
	STOP'EXIT from GSSIGN : no such surface !'
C
	END



	SUBROUTINE GEOINI
C****************************************************************
C* Initialization                                               *
C****************************************************************
	IMPLICIT REAL*8 (A-H,R-Z)
C
	PARAMETER(MAXB=1000,MAXZ=1000)
	PARAMETER(JDIMBD=37*MAXB,JDIMZD=5*MAXZ)
	PARAMETER(MKINDB=10,MKINDS=10)
C
	COMMON/GDATA0/ NUMBOD,NUMZON
	COMMON/GDATA1/ ISYMB(MKINDB),NBSURF(MKINDB),NSPAR(MKINDS)
	COMMON/GDATA2/ BODYDA(MAXB),BODYDB(JDIMBD)
	INTEGER BODYDA
	COMMON/GDATA3/ ZONEDA(MAXZ),ZONEDB(JDIMZD)
	INTEGER ZONEDA,ZONEDB
	COMMON/GDATA4/ JFLSUF(3,MAXB),BTRZLT(MAXB),
     *                 DIGINF,DIGMIN,DELTA,STPLIM

      CHARACTER*80 ffor017(1000)
      INTEGER f017In, f017Out
      COMMON/FILEFOR017/ ffor017, f017In, f017Out
C
C-------------------------------------------------------------
C                     1   2   3   4   5   6   7   8   9   10 |
C Body               SPH WED ARB BOX RPP RCC REC TRC ELL END |
C N Surfaces          1   5   6   6   6   3   3   3   1   0  |
C-------------------------------------------------------------
C-------------------------------------------------------------
C Surface type:      Equation:                 N Parameters: |
C             1      X-C=0                                 1 |
C             2      Y-C=0                                 1 |
C             3      Z-C=0                                 1 |
C             4      A*x+B*y+C*z-D=0                       4 |
C             5      (x-x0)**2+(y-y0)**2+(z-z0)**2-R**2=0  4 |
C             6      x**2+y**2-R**2=0                      1 |
C             7      x**2/A2+y**2/B2-1=0                   2 |
C             8      x**2+y**2-(z-A1)**2/B2=0              2 |
C             9      x**2/A2+y**2/B2+z**2/C2-1=0           3 |
C-------------------------------------------------------------
C
	DIMENSION RTMP(30)
	DIMENSION JTMP(20)
	DIMENSION IZSYM(20)
	DIMENSION JTY1(MAXZ),JTY2(MAXZ)          ! 21 may 1998 AVDem
C
	DIMENSION IBCHK(MAXB)
C
      DATA IOR,IBL,IEND/2HOR,3H   ,3HEND/
      DATA IMOV,IROT/3HMOV,3HROT/
C
C=========================================================
      DIGORT=0.0000001   ! The measure of nonorthogonality
C                        ! Introduced by Sobol, 02.99
C=========================================================
C
	DO J=1,MAXB
	   IBCHK(J)=-1
	ENDDO
C
	JIN=16
	JERR=25
	f017Out=0
	f017In=0
C
	READ(JIN,20)JDBG1,JDBG2,(JTMP(I),I=1,15)  ! Input of title line
20      FORMAT(2I5,10X,15A4)
C
	IF(JDBG2.NE.0)THEN
	STPLIM=10.0**DBLE(-ABS(JDBG2))
	ENDIF
C
	NUMBOD=0 ! Number of bodies in the target description
	JBADDR=1 ! Address of the body; (where?)
C       I=0      ! AVDEM
C           Sobol: error of AVDEM!: variable I is used below in other way!
	ISOB=0   ! ISOB instead of I is introduced by Sobol, 02.99
C                ! ISOB is a number of the body: 1- SPH, 2-WED, etc.
C
97      FORMAT(2X,A3,1X,I4,6D10.3)
98      FORMAT(10X,6D10.3)
99      FORMAT(10X,6(I9,1X))
C
  100   CONTINUE  
        CALL GBWRT(ffor017,ISOB,RTMP,JTMP,f017Out) ! Sobol: ISOB instead of I
C                  ! At ISOB=0 GBWRT writes the title line
C
	READ( JIN,97) ITYPE,IALP,(RTMP(I),I=1,6)
	NBT=0
C
101       DO I=1,MKINDB
	IF( ISYMB(I).EQ.ITYPE )then
         ISOB=I                ! The line is added by Sobol, 02.99
         GOTO 102
        end if
	  END DO
C
	IF(ITYPE.EQ.IMOV)THEN
	   CALL GBMOV(ffor017,NBT,RTMP,JTMP,f017In)
	   ISOB=NBT                       ! Sobol: ISOB instead of I
	   GOTO 102
	ELSEIF(ITYPE.EQ.IROT)THEN
	   CALL GBROT(ffor017,NBT,RTMP,JTMP,f017In)
	   ISOB=NBT                       ! Sobol: ISOB instead of I
	   GOTO 102
	ENDIF
C
	 WRITE (JERR,103)ITYPE,ISYMB
103     FORMAT(' ITYPE=',A3,'  Does not equal any of the following: ',
     *       10(A3,2X))
	 PRINT *,'GEOINI: body number ',(NUMBOD+1)
	IF(ITYPE.EQ.IBL)THEN
	 PRINT *,' "   " body type - check body ',NUMBOD
	 PRINT *,' on exceeding of description lines number.'
	 PRINT *,'                    '
	ENDIF
	 STOP'EXIT from GEOINI : invalid body type !'
C
  102    CONTINUE
         NUMBOD=NUMBOD+1
	IF(ISOB.EQ.10)THEN            ! Sobol: ISOB instead of I
	 NUMBOD=NUMBOD-1
	 GOTO 200
	ELSE
	 BODYDA(NUMBOD)=JBADDR
	 BODYDB(JBADDR)=ISOB          ! Sobol: ISOB instead of I
	 JBADDR=JBADDR+1
C
	   IF(NBT.EQ.0)THEN
	      GOTO ( 1, 2, 3, 4, 5, 6, 7, 8, 9),ISOB
	   ELSE                        ! Sobol: ISOB instead of I
	      GOTO ( 1, 502, 503, 504, 5, 506, 507, 508, 509),ISOB
	   ENDIF                       ! Sobol: ISOB instead of I
	ENDIF
C
C----- BODY =  "SPH" -----
1       IF(RTMP(4).LT.0.)THEN
	PRINT *,'SPH: body number ',NUMBOD
	STOP'EXIT from GEOINI: negative radius !'
	ENDIF
C
	BODYDB(JBADDR)=5
	 RTMP(4)=RTMP(4)**2.
	 RTMP(5)=-1.
	DO K=1,5
	 BODYDB(JBADDR+K)=RTMP(K)
	END DO
	 JBADDR=JBADDR+6
C
	 GOTO 100
C
C----- BODY =  "WED" -----
2       READ(JIN,98) (RTMP(I),I=7,12)
C
502     DO I=1,2
	DORT=RTMP(4)*RTMP(4+I*3)+RTMP(5)*RTMP(5+I*3)+RTMP(6)*RTMP(6+I*3)
c        IF(DABS(DORT).GT.DIGMIN)THEN               ! AVDEM
        IF(DABS(DORT).GT.DIGORT)THEN               ! Sobol
	PRINT *,'WED: body number ',NUMBOD,' i=',I
	STOP'EXIT from GEOINI: vector H not orthogonal Ri !'
	ENDIF
	ENDDO
C
	X=RTMP(1)+0.25*(RTMP(4)+RTMP(7)+RTMP(10))
	Y=RTMP(2)+0.25*(RTMP(5)+RTMP(8)+RTMP(11))
	Z=RTMP(3)+0.25*(RTMP(6)+RTMP(9)+RTMP(12))
C
	BODYDB(JBADDR)=4
	CALL GEQPLN(RTMP(1),RTMP(2),RTMP(3),RTMP(4),RTMP(5),RTMP(6),
     * BODYDB(JBADDR+1),BODYDB(JBADDR+2),
     * BODYDB(JBADDR+3),BODYDB(JBADDR+4))
	BODYDB(JBADDR+5)=1.
	BODYDB(JBADDR+5)=GSSIGN(4,JBADDR+1,X,Y,Z)
C
	BODYDB(JBADDR+6)=4
	CALL GEQPLN(RTMP(1)+RTMP(4),RTMP(2)+RTMP(5),RTMP(3)+RTMP(6),
     * RTMP(4),RTMP(5),RTMP(6),
     * BODYDB(JBADDR+7),BODYDB(JBADDR+8),
     * BODYDB(JBADDR+9),BODYDB(JBADDR+10))
	BODYDB(JBADDR+11)=1.
	BODYDB(JBADDR+11)=GSSIGN(4,JBADDR+7,X,Y,Z)
C
	BODYDB(JBADDR+12)=4
	CALL GVEC90(1,RTMP(4),RTMP(5),RTMP(6),RTMP(7),RTMP(8),RTMP(9),
     *   V1,V2,V3)
	CALL GEQPLN(RTMP(1),RTMP(2),RTMP(3),V1,V2,V3,
     * BODYDB(JBADDR+13),BODYDB(JBADDR+14),
     * BODYDB(JBADDR+15),BODYDB(JBADDR+16))
	BODYDB(JBADDR+17)=1.
	BODYDB(JBADDR+17)=GSSIGN(4,JBADDR+13,X,Y,Z)
C
	BODYDB(JBADDR+18)=4
	CALL GVEC90(1,RTMP(4),RTMP(5),RTMP(6),RTMP(10),RTMP(11),RTMP(12),
     *   V1,V2,V3)
	CALL GEQPLN(RTMP(1),RTMP(2),RTMP(3),V1,V2,V3,
     * BODYDB(JBADDR+19),BODYDB(JBADDR+20),
     * BODYDB(JBADDR+21),BODYDB(JBADDR+22))
	BODYDB(JBADDR+23)=1.
	BODYDB(JBADDR+23)=GSSIGN(4,JBADDR+19,X,Y,Z)
C
	VX=RTMP(10)-RTMP(7)
	VY=RTMP(11)-RTMP(8)
	VZ=RTMP(12)-RTMP(9)
C
	BODYDB(JBADDR+24)=4
	CALL GVEC90(1,RTMP(4),RTMP(5),RTMP(6),VX,VY,VZ,V1,V2,V3)
	CALL GEQPLN(RTMP(1)+RTMP(7),RTMP(2)+RTMP(8),RTMP(3)+RTMP(9),
     * V1,V2,V3,
     * BODYDB(JBADDR+25),BODYDB(JBADDR+26),
     * BODYDB(JBADDR+27),BODYDB(JBADDR+28))
	BODYDB(JBADDR+29)=1.
	BODYDB(JBADDR+29)=GSSIGN(4,JBADDR+25,X,Y,Z)
C
	JBADDR=JBADDR+30
	GOTO 100
C
C----- BODY =  "ARB" -----
3       READ(JIN,98) (RTMP(I),I=7,12)
	READ(JIN,98) (RTMP(I),I=13,18)
	READ(JIN,98) (RTMP(I),I=19,24)
	READ(JIN,99) (JTMP(I),I=1,6)
C
503     X=(RTMP(1)+RTMP(4)+RTMP(7)+RTMP(10)+
     *     RTMP(13)+RTMP(16)+RTMP(19)+RTMP(22))/8.
	Y=(RTMP(2)+RTMP(5)+RTMP(8)+RTMP(11)+
     *     RTMP(14)+RTMP(17)+RTMP(20)+RTMP(23))/8.
	Z=(RTMP(3)+RTMP(6)+RTMP(9)+RTMP(12)+
     *     RTMP(15)+RTMP(18)+RTMP(21)+RTMP(24))/8.
C
	DO I=1,6
	IP1=JTMP(I)/1000
	IP2=(JTMP(I)-IP1*1000)/100
	IP3=(JTMP(I)-IP1*1000-IP2*100)/10
	IP4=JTMP(I)-IP1*1000-IP2*100-IP3*10
C
	IP1=(IP1-1)*3+1
	IP2=(IP2-1)*3+1
	IP3=(IP3-1)*3+1
	IP4=(IP4-1)*3+1
C
	CALL GVEC90(1,RTMP(IP1)-RTMP(IP2),RTMP(IP1+1)-RTMP(IP2+1),
     *                RTMP(IP1+2)-RTMP(IP2+2),
     *  RTMP(IP1)-RTMP(IP3),RTMP(IP1+1)-RTMP(IP3+1),
     *  RTMP(IP1+2)-RTMP(IP3+2),
     *   V1,V2,V3)
C
	BODYDB(JBADDR+6*(I-1))=4
	CALL GEQPLN(RTMP(IP1),RTMP(IP1+1),RTMP(IP1+2),V1,V2,V3,
     * BODYDB(JBADDR+6*(I-1)+1),BODYDB(JBADDR+6*(I-1)+2),
     * BODYDB(JBADDR+6*(I-1)+3),BODYDB(JBADDR+6*(I-1)+4))
	BODYDB(JBADDR+6*(I-1)+5)=1.
	BODYDB(JBADDR+6*(I-1)+5)=GSSIGN(4,JBADDR+6*(I-1)+1,X,Y,Z)
C
	CHECKS=GSSIGN(4,JBADDR+6*(I-1)+1,
     *                  RTMP(IP4),RTMP(IP4+1),RTMP(IP4+2))
	IF(DABS(CHECKS).GT.DIGMIN)THEN
	PRINT *,'ARB: body number ',NUMBOD
	PRINT *,'Plane: ',JTMP(I)
	PRINT *,'Check point coordinates.'
	PRINT *,'The last point don`t lie in the plane'
	PRINT *,'that has been drawn through the three'
	PRINT *,'previous points.'
	STOP'EXIT from GEOINI: Invalid coordinates !'
	ENDIF
C
	ENDDO
C
	JBADDR=JBADDR+36
	GOTO 100
C
C----- BODY =  "BOX" -----
4       READ(JIN,98) (RTMP(I),I=7,12)
C
504     DO I=1,2
	DORT=RTMP(4)*RTMP(4+I*3)+RTMP(5)*RTMP(5+I*3)+RTMP(6)*RTMP(6+I*3)
c        IF(DABS(DORT).GT.DIGMIN)THEN              ! AVDEM
        IF(DABS(DORT).GT.DIGORT)THEN              ! Sobol
	PRINT *,'BOX: body number ',NUMBOD,'        i=',I+1
	STOP'EXIT from GEOINI: vector Vec_1 not orthogonal Vec_i !'
	ENDIF
	ENDDO
	DORT=RTMP(7)*RTMP(10)+RTMP(8)*RTMP(11)+RTMP(9)*RTMP(12)
c        IF(DABS(DORT).GT.DIGMIN)THEN              ! AVDEM
        IF(DABS(DORT).GT.DIGORT)THEN              ! Sobol
	PRINT *,'BOX: body number ',NUMBOD
	STOP'EXIT from GEOINI: vector Vec_2 not orthogonal Vec_3 !'
	ENDIF
C
	X=RTMP(1)+0.5*(RTMP(4)+RTMP(7)+RTMP(10))
	Y=RTMP(2)+0.5*(RTMP(5)+RTMP(8)+RTMP(11))
	Z=RTMP(3)+0.5*(RTMP(6)+RTMP(9)+RTMP(12))
C
	  DO I=1,3
	     K=3*I-2
	BODYDB(JBADDR+12*(I-1))=4
	CALL GEQPLN(RTMP(1),RTMP(2),RTMP(3),
     *              RTMP(K+3),RTMP(K+4),RTMP(K+5),
     * BODYDB(JBADDR+12*(I-1)+1),BODYDB(JBADDR+12*(I-1)+2),
     * BODYDB(JBADDR+12*(I-1)+3),BODYDB(JBADDR+12*(I-1)+4))
	BODYDB(JBADDR+12*(I-1)+5)=1.
	BODYDB(JBADDR+12*(I-1)+5)=GSSIGN(4,JBADDR+12*(I-1)+1,X,Y,Z)
C
	BODYDB(JBADDR+12*(I-1)+6)=4
	CALL GEQPLN(RTMP(1)+RTMP(K+3),RTMP(2)+RTMP(K+4),RTMP(3)+RTMP(K+5),
     *              RTMP(K+3),RTMP(K+4),RTMP(K+5),
     * BODYDB(JBADDR+12*(I-1)+7),BODYDB(JBADDR+12*(I-1)+8),
     * BODYDB(JBADDR+12*(I-1)+9),BODYDB(JBADDR+12*(I-1)+10))
	BODYDB(JBADDR+12*(I-1)+11)=1.
	BODYDB(JBADDR+12*(I-1)+11)=GSSIGN(4,JBADDR+12*(I-1)+7,X,Y,Z)
C
	  ENDDO
C
	JBADDR=JBADDR+36
	GOTO 100
C
C----- BODY =  "RPP" -----
5       IF(RTMP(1).GE.RTMP(2))THEN
	PRINT *,'RPP: body number ',NUMBOD
	STOP'EXIT from GEOINI: X1 must be less than X2 !'
	ENDIF
	IF(RTMP(3).GE.RTMP(4))THEN
	PRINT *,'RPP: body number ',NUMBOD
	STOP'EXIT from GEOINI: Y1 must be less than Y2 !'
	ENDIF
	IF(RTMP(5).GE.RTMP(6))THEN
	PRINT *,'RPP: body number ',NUMBOD
	STOP'EXIT from GEOINI: Z1 must be less than Z2 !'
	ENDIF
C
	BODYDB(JBADDR)=1
	BODYDB(JBADDR+1)=RTMP(1)
	BODYDB(JBADDR+2)=1.
C
	BODYDB(JBADDR+3)=1
	BODYDB(JBADDR+4)=RTMP(2)
	BODYDB(JBADDR+5)=-1.
C
	BODYDB(JBADDR+6)=2
	BODYDB(JBADDR+7)=RTMP(3)
	BODYDB(JBADDR+8)=1.
C
	BODYDB(JBADDR+9)=2
	BODYDB(JBADDR+10)=RTMP(4)
	BODYDB(JBADDR+11)=-1.
C
	BODYDB(JBADDR+12)=3
	BODYDB(JBADDR+13)=RTMP(5)
	BODYDB(JBADDR+14)=1.
C
	BODYDB(JBADDR+15)=3
	BODYDB(JBADDR+16)=RTMP(6)
	BODYDB(JBADDR+17)=-1.
C
	JBADDR=JBADDR+18
	GOTO 100
C
C----- BODY =  "RCC" -----
6       READ(JIN,98) RTMP(7)
C
506     IF(RTMP(7).LT.0.)THEN
	PRINT *,'RCC: body number ',NUMBOD
	STOP'EXIT from GEOINI: negative radius !'
	ENDIF
C
	CALL GVEC90(2,RTMP(4),RTMP(5),RTMP(6),A1,B1,C1,A2,B2,C2)
C
	CALL GCOSAX(RTMP(1),RTMP(2),RTMP(3),A1,B1,C1,
     *           BODYDB(JBADDR),BODYDB(JBADDR+1),
     *           BODYDB(JBADDR+2),BODYDB(JBADDR+3),R1,R12)
C
	CALL GCOSAX(RTMP(1),RTMP(2),RTMP(3),A2,B2,C2,
     *           BODYDB(JBADDR+4),BODYDB(JBADDR+5),
     *           BODYDB(JBADDR+6),BODYDB(JBADDR+7),R2,R22)
C
	CALL GCOSAX(RTMP(1),RTMP(2),RTMP(3),RTMP(4),RTMP(5),RTMP(6),
     *           BODYDB(JBADDR+8),BODYDB(JBADDR+9),
     *           BODYDB(JBADDR+10),BODYDB(JBADDR+11),H,H2)
C
	JBADDR=JBADDR+12
	BODYDB(JBADDR)=3
	BODYDB(JBADDR+1)=0.
	BODYDB(JBADDR+2)=1.
	BODYDB(JBADDR+3)=3
	BODYDB(JBADDR+4)=H
	BODYDB(JBADDR+5)=-1.
	BODYDB(JBADDR+6)=6
	BODYDB(JBADDR+7)=RTMP(7)**2.
	BODYDB(JBADDR+8)=-1.
	JBADDR=JBADDR+9
	GOTO 100
C
C----- BODY =  "REC" -----
7       READ(JIN,98) (RTMP(I),I=7,12)
C
507     DO I=1,2
	DORT=RTMP(4)*RTMP(4+I*3)+RTMP(5)*RTMP(5+I*3)+RTMP(6)*RTMP(6+I*3)
c        IF(DABS(DORT).GT.DIGMIN)THEN                ! AVDEM
        IF(DABS(DORT).GT.DIGORT)THEN                ! Sobol
	PRINT *,'REC: body number ',NUMBOD,' i=',I
	STOP'EXIT from GEOINI: vector H not orthogonal Ri !'
	ENDIF
	ENDDO
	DORT=RTMP(7)*RTMP(10)+RTMP(8)*RTMP(11)+RTMP(9)*RTMP(12)
c        IF(DABS(DORT).GT.DIGMIN)THEN                ! AVDEM
        IF(DABS(DORT).GT.DIGORT)THEN                ! Sobol
	PRINT *,'REC: body number ',NUMBOD
	STOP'EXIT from GEOINI: vector R1 not orthogonal R2 !'
	ENDIF
C
	CALL GCOSAX(RTMP(1),RTMP(2),RTMP(3),RTMP(7),RTMP(8),RTMP(9),
     *           BODYDB(JBADDR),BODYDB(JBADDR+1),
     *           BODYDB(JBADDR+2),BODYDB(JBADDR+3),R1,R12)
C
	CALL GCOSAX(RTMP(1),RTMP(2),RTMP(3),RTMP(10),RTMP(11),RTMP(12),
     *           BODYDB(JBADDR+4),BODYDB(JBADDR+5),
     *           BODYDB(JBADDR+6),BODYDB(JBADDR+7),R2,R22)
C
	CALL GCOSAX(RTMP(1),RTMP(2),RTMP(3),RTMP(4),RTMP(5),RTMP(6),
     *           BODYDB(JBADDR+8),BODYDB(JBADDR+9),
     *           BODYDB(JBADDR+10),BODYDB(JBADDR+11),H,H2)
C
	JBADDR=JBADDR+12
	BODYDB(JBADDR)=3
	BODYDB(JBADDR+1)=0.
	BODYDB(JBADDR+2)=1.
	BODYDB(JBADDR+3)=3
	BODYDB(JBADDR+4)=H
	BODYDB(JBADDR+5)=-1.
	BODYDB(JBADDR+6)=7
	BODYDB(JBADDR+7)=R12
	BODYDB(JBADDR+8)=R22
	BODYDB(JBADDR+9)=-1.
	JBADDR=JBADDR+10
	GOTO 100
C
C----- BODY =  "TRC" -----
8       READ(JIN,98) (RTMP(I),I=7,8)
C
508     IF(RTMP(8).GE.RTMP(7))THEN
	PRINT *,'TRC: body number ',NUMBOD
	PRINT *,'Check vector H.'
	STOP'EXIT from GEOINI: R2 must be less than R1 !'
	ENDIF
C
	CALL GVEC90(2,RTMP(4),RTMP(5),RTMP(6),A1,B1,C1,A2,B2,C2)
C
	CALL GCOSAX(RTMP(1),RTMP(2),RTMP(3),A1,B1,C1,
     *           BODYDB(JBADDR),BODYDB(JBADDR+1),
     *           BODYDB(JBADDR+2),BODYDB(JBADDR+3),R1,R12)
C
	CALL GCOSAX(RTMP(1),RTMP(2),RTMP(3),A2,B2,C2,
     *           BODYDB(JBADDR+4),BODYDB(JBADDR+5),
     *           BODYDB(JBADDR+6),BODYDB(JBADDR+7),R2,R22)
C
	CALL GCOSAX(RTMP(1),RTMP(2),RTMP(3),RTMP(4),RTMP(5),RTMP(6),
     *           BODYDB(JBADDR+8),BODYDB(JBADDR+9),
     *           BODYDB(JBADDR+10),BODYDB(JBADDR+11),H,H2)
C
	JBADDR=JBADDR+12
	BODYDB(JBADDR)=3
	BODYDB(JBADDR+1)=0.
	BODYDB(JBADDR+2)=1.
	BODYDB(JBADDR+3)=3
	BODYDB(JBADDR+4)=H
	BODYDB(JBADDR+5)=-1.
	BODYDB(JBADDR+6)=8
	BODYDB(JBADDR+7)=H/(1.-RTMP(8)/RTMP(7))
	BODYDB(JBADDR+8)=BODYDB(JBADDR+7)**2./RTMP(7)**2.
	BODYDB(JBADDR+9)=-1.
	JBADDR=JBADDR+10
	GOTO 100
C
C----- BODY =  "ELL" -----
9       READ(JIN,98) (RTMP(I),I=7,12)
C
509     DO I=1,2
	DORT=RTMP(4)*RTMP(4+I*3)+RTMP(5)*RTMP(5+I*3)+RTMP(6)*RTMP(6+I*3)
c        IF(DABS(DORT).GT.DIGMIN)THEN           ! AVDEM
        IF(DABS(DORT).GT.DIGORT)THEN           ! Sobol
	PRINT *,'ELL: body number ',NUMBOD,'        i=',I+1
	STOP'EXIT from GEOINI: vector Vec_1 not orthogonal Vec_i !'
	ENDIF
	DORT=RTMP(7)*RTMP(10)+RTMP(8)*RTMP(11)+RTMP(9)*RTMP(12)
c        IF(DABS(DORT).GT.DIGMIN)THEN           ! AVDEM
        IF(DABS(DORT).GT.DIGORT)THEN           ! Sobol
	PRINT *,'ELL: body number ',NUMBOD
	STOP'EXIT from GEOINI: vector Vec_2 not orthogonal Vec_3 !'
	ENDIF
	ENDDO
C
C            GCOSAX(X0,Y0,Z0,X,Y,Z,*COSX,*COSY,*COSZ,*DELTAB,*R,*R2)
C
	CALL GCOSAX(RTMP(1),RTMP(2),RTMP(3),RTMP(4),RTMP(5),RTMP(6),
     *         BODYDB(JBADDR),BODYDB(JBADDR+1),
     *         BODYDB(JBADDR+2),BODYDB(JBADDR+3),R,BODYDB(JBADDR+13))
C
	CALL GCOSAX(RTMP(1),RTMP(2),RTMP(3),RTMP(7),RTMP(8),RTMP(9),
     *         BODYDB(JBADDR+4),BODYDB(JBADDR+5),
     *         BODYDB(JBADDR+6),BODYDB(JBADDR+7),R,BODYDB(JBADDR+14))
C
	CALL GCOSAX(RTMP(1),RTMP(2),RTMP(3),RTMP(10),RTMP(11),RTMP(12),
     *         BODYDB(JBADDR+8),BODYDB(JBADDR+9),
     *         BODYDB(JBADDR+10),BODYDB(JBADDR+11),R,BODYDB(JBADDR+15))
C
	BODYDB(JBADDR+12)=9
	BODYDB(JBADDR+16)=-1.
	JBADDR=JBADDR+17
	GOTO 100
C
200     CALL GBWRT(ffor017,-1,RTMP,JTMP,f017Out)
C
	IF(NUMBOD.GT.MAXB)THEN
	   PRINT *,'ERROR - body number limit is exceeded!'
	   PRINT *,'Current body number limit is ',MAXB
	   STOP'EXIT from GEOINI: correct Fortran text.'
	ENDIF
C
C----- END OF BODY DATA ---------------
C
	NUMZON=0
	JZADDR=1
	JZTADD=1
	NUMZB=0
	ZONEDB(3)=1
C
201     FORMAT(2X,A3,I5,9(A2,I5))
202     READ(JIN,201)IALP,NAZ,(IZSYM(I),JTMP(I),I=1,9)
C
	IF( IALP.EQ.IEND)GOTO 299
	IF( IALP.EQ.IBL )GOTO 203
	 ZONEDB(JZTADD)=NUMZB
	 NUMZON=NUMZON+1
	 JZTADD=JZADDR
	 JZADDR=JZADDR+3
	 NUMZB=0
	 ZONEDA(NUMZON)=JZTADD
	 ZONEDB(JZTADD+2)=1
C
203     DO L=1,9
	IF(JTMP(L).NE.0)THEN
	 IF(IZSYM(L).EQ.IOR)THEN
	   NUMZB=NUMZB+1
	   ZONEDB(JZADDR)=0
		  JZADDR=JZADDR+1
	   ZONEDB(JZTADD+2)=0
	 ENDIF
	   NUMZB=NUMZB+1
	   ZONEDB(JZADDR)=JTMP(L)
	   IIIN=ABS(JTMP(L))
	 IF(IIIN.GT.NUMBOD)THEN
	   PRINT *,'Zone - ',NUMZON
	   PRINT *,'Body - ',JTMP(L)
	   STOP'EXIT from GEOINI: undefined body number!'
	 ELSE
	   IBCHK(IIIN)=1
	 ENDIF
	 IF(JTMP(L).LT.0 .AND. ZONEDB(JZTADD+2).NE.0)ZONEDB(JZTADD+2)=-1
	 JZADDR=JZADDR+1
	ELSE
	 GOTO 202
	END IF
	END DO
	 GOTO 202
C
299      ZONEDB(JZTADD)=NUMZB
C
	 DO J=1,NUMBOD
	  IF(IBCHK(J).LT.0)THEN
	     PRINT *,'Body - ',J
	     PRINT *,'WARNING - this body is not used in zone description'
	  ENDIF
	 ENDDO
C
	 IF(NUMZON.GT.MAXZ)THEN
	   PRINT *,'ERROR - zone number limit is exceeded!'
	   PRINT *,'Current zone number limit is ',MAXZ
	   STOP'EXIT from GEOINI: correct FORTRAN text.'
	ENDIF
C
	L=(NUMZON+13)/14
C
301     FORMAT(14I5)
C
	DO I=1,L
	READ(JIN,301)(JTY1(J+(I-1)*14),J=1,14)
	END DO
C
	DO I=1,L
	READ(JIN,301)(JTY2(J+(I-1)*14),J=1,14)
	END DO
C
	  DO I=1,NUMZON
	IF(JTY1(I).EQ.0)THEN
	 PRINT *,'Full number of zone=',NUMZON
	 PRINT *,'Number of zone with media definition=',I-1
	 PRINT *,'                    '
	 STOP'EXIT from GEOINI: undefined zone media !'
	ENDIF
C
	JZADDR=ZONEDA(I)
	ZONEDB(JZADDR+1)=JTY2(I)
	  END DO
C
	 CLOSE(JIN,IOSTAT=IFSTA1,ERR=1000)
	IF(IFSTA1.NE.0)THEN
	 PRINT *,'GEOINI: Some file error was encountered'
	 PRINT *,'Input  file IOSTAT =',IFSTA1
	 PRINT *,'Output file IOSTAT =',IFSTA2 !We have not output file
	 PRINT *,'                    '
	ENDIF
1001     PRINT *,'GEOINI: initialization has completed successfully .'
	RETURN
1000     PRINT *,'GEOINI: Problem with files closing .'
	 PRINT *,'                    '
	GOTO 1001
C
	END



	SUBROUTINE GCOSAX(X0,Y0,Z0,X,Y,Z,COSX,COSY,COSZ,DELTA,R,R2)
C****************************************************************
C*                                                              *
C****************************************************************
	IMPLICIT REAL*8 (A-H,R-Z)
C
	R2=X**2.+Y**2.+Z**2.
	R=DSQRT(R2)
C
	COSX=X/R
	COSY=Y/R
	COSZ=Z/R
C
	DELTA=X0*COSX+Y0*COSY+Z0*COSZ
C
	RETURN
	END



	SUBROUTINE GEQPLN(X,Y,Z,V1,V2,V3,A,B,C,D)
C****************************************************************
C*                                                              *
C****************************************************************
	IMPLICIT REAL*8 (A-H,R-Z)
C
	D=V1*X+V2*Y+V3*Z
C
	 IF(D.EQ.0.)THEN
	A=V1
	B=V2
	C=V3
	 ELSE
	A=V1/D
	B=V2/D
	C=V3/D
	D=1.
	 ENDIF
C
	RETURN
	END



	SUBROUTINE GVEC90(NV,X1,Y1,Z1,X2,Y2,Z2,A,B,C)
C****************************************************************
C*                                                              *
C****************************************************************
	IMPLICIT REAL*8 (A-H,R-Z)
C
	GOTO ( 2, 1),NV
C
1       CONTINUE
C
	IF(X1.EQ.0.)THEN
	 X2=1.
	ELSE
	 X2=-Z1
	END IF
C
	IF(Y1.EQ.0.)THEN
	 Y2=1.
	ELSE
	 Y2=0.
	END IF
C
	IF(Z1.EQ.0.)THEN
	 Z2=1.
	ELSE
	 Z2=X1
	END IF
C
2       CONTINUE
C
	IF(X1.EQ.0. .AND. X2.EQ.0.)THEN
	 A=1.
	 B=0.
	 C=0.
	RETURN
	END IF
C
	IF(Y1.EQ.0. .AND. Y2.EQ.0.)THEN
	 A=0.
	 B=1.
	 C=0.
	RETURN
	END IF
C
	IF(Z1.EQ.0. .AND. Z2.EQ.0.)THEN
	 A=0.
	 B=0.
	 C=1.
	RETURN
	END IF
C
	IF(X1.EQ.0. .AND. Y1.EQ.0.)THEN
	 A=-Y2
	 B=X2
	 C=0.
	RETURN
	END IF
C
	IF(X1.EQ.0. .AND. Z1.EQ.0.)THEN
	 A=-Z2
	 B=0.
	 C=X2
	RETURN
	END IF
C
	IF(Y1.EQ.0. .AND. Z1.EQ.0.)THEN
	 A=0.
	 B=-Z2
	 C=Y2
	RETURN
	END IF
C
	IF(X2.EQ.0. .AND. Y2.EQ.0.)THEN
	 A=-Y1
	 B=X1
	 C=0.
	RETURN
	END IF
C
	IF(X2.EQ.0. .AND. Z2.EQ.0.)THEN
	 A=-Z1
	 B=0.
	 C=X1
	RETURN
	END IF
C
	IF(Y2.EQ.0. .AND. Z2.EQ.0.)THEN
	 A=0.
	 B=-Z1
	 C=Y1
	RETURN
	END IF
C
	IF(X1.NE.0. .AND. X2.NE.0.)THEN
	 IF(Y1/X1.EQ.Y2/X2)THEN
	A=Y1
	B=-X1
	C=0.
	 ELSE
	C=1.
	B=(Z2*X1-Z1*X2)/(Y1*X2-Y2*X1)
	A=-(B*Y1+Z1)/X1
	 END IF
	RETURN
	END IF
C
	IF(Y1.NE.0. .AND. Y2.NE.0.)THEN
	 IF(Z1/Y1.EQ.Z2/Y2)THEN
	A=0.
	B=-Z1
	C=Y1
	 ELSE
	A=1.
	C=(X2*Y1-X1*Y2)/(Z1*Y2-Z2*Y1)
	B=-(C*Z1+X1)/Y1
	 END IF
	RETURN
	END IF
C
	IF(Z1.NE.0. .AND. Z2.NE.0.)THEN
	 IF(X1/Z1.EQ.X2/Z2)THEN
	A=-Z1
	B=0.
	C=X1
	 ELSE
	B=1.
	A=(Y2*Z1-Y1*Z2)/(X1*Z2-X2*Z1)
	C=-(A*X1+Y1)/Z1
	 END IF
	RETURN
	END IF
C
	STOP'EXIT from GVEC90 : vector not found !'
	END




	SUBROUTINE GBMOV(JREAD,NBT,RTMP,JTMP,f017In)
C****************************************************************
C* MOV card interpretation                                      *
C****************************************************************
C
	IMPLICIT REAL*8 (A-H,R-Z)
C
	PARAMETER(MAXB=1000,MAXZ=1000)
	PARAMETER(JDIMBD=37*MAXB,JDIMZD=5*MAXZ)
	PARAMETER(MKINDB=10,MKINDS=10)
C
	COMMON/GDATA0/ NUMBOD,NUMZON
	COMMON/GDATA1/ ISYMB(MKINDB),NBSURF(MKINDB),NSPAR(MKINDS)
	COMMON/GDATA2/ BODYDA(MAXB),BODYDB(JDIMBD)
	INTEGER BODYDA
	CHARACTER*80 JREAD(1000)
      INTEGER f017In
C
	DIMENSION RBDATA(30)
	DIMENSION JBDATA(20)
	DIMENSION RTMP(30)    ! Sobol: (30) instead of (*) to printout
	DIMENSION JTMP(30)    ! Sobol: (30) instead of (*) to printout
C
	NB=INT(RTMP(1))
	NBT=INT(BODYDB(BODYDA(NB)))
	CALL GBREAD(JREAD,NB,RBDATA,JBDATA,f017In)
	CALL GFEND(JREAD,f017In)
C
	X=RTMP(2)
	Y=RTMP(3)
	Z=RTMP(4)
C
	IF(X.EQ.0.0 .AND. Y.EQ.0.0 .AND. Z.EQ.0.0)THEN
	PRINT *,'Body number =',NUMBOD+1
	PRINT *,'MOV on (0.0 , 0.0 , 0.0)'
	STOP'EXIT from GBMOV: can not move!'
	ENDIF
C
	GOTO ( 1, 2, 3, 4, 5, 6, 7, 8, 9),NBT
C
C----- BODY =  "SPH" -----
1       RTMP(1)=RBDATA(1)+X
	RTMP(2)=RBDATA(2)+Y
	RTMP(3)=RBDATA(3)+Z
	RTMP(4)=RBDATA(4)
C
	 RETURN
C
C----- BODY =  "WED" -----
2       RTMP(1)=RBDATA(1)+X
	RTMP(2)=RBDATA(2)+Y
	RTMP(3)=RBDATA(3)+Z
	  DO I=4,12
	RTMP(I)=RBDATA(I)
	  ENDDO
C
	RETURN
C
C----- BODY =  "ARB" -----
3         DO I=1,8
	RTMP(3*(I-1)+1)=RBDATA(3*(I-1)+1)+X
	RTMP(3*(I-1)+2)=RBDATA(3*(I-1)+2)+Y
	RTMP(3*(I-1)+3)=RBDATA(3*(I-1)+3)+Z
	  ENDDO
	  DO I=1,6
	JTMP(I)=JBDATA(I)
	  ENDDO
C
	RETURN
C
C----- BODY =  "BOX" -----
4       RTMP(1)=RBDATA(1)+X
	RTMP(2)=RBDATA(2)+Y
	RTMP(3)=RBDATA(3)+Z
	  DO I=4,12
	RTMP(I)=RBDATA(I)
	  ENDDO
C
	RETURN
C
C----- BODY =  "RPP" -----
5       RTMP(1)=RBDATA(1)+X
	RTMP(2)=RBDATA(2)+X
	RTMP(3)=RBDATA(3)+Y
	RTMP(4)=RBDATA(4)+Y
	RTMP(5)=RBDATA(5)+Z
	RTMP(6)=RBDATA(6)+Z
C
	RETURN
C
C----- BODY =  "RCC" -----
6       RTMP(1)=RBDATA(1)+X
	RTMP(2)=RBDATA(2)+Y
	RTMP(3)=RBDATA(3)+Z
	  DO I=4,7
	RTMP(I)=RBDATA(I)
	  ENDDO
C
	RETURN
C
C----- BODY =  "REC" -----
7       RTMP(1)=RBDATA(1)+X
	RTMP(2)=RBDATA(2)+Y
	RTMP(3)=RBDATA(3)+Z
	  DO I=4,12
	RTMP(I)=RBDATA(I)
	  ENDDO
C
	RETURN
C
C----- BODY =  "TRC" -----
8       RTMP(1)=RBDATA(1)+X
	RTMP(2)=RBDATA(2)+Y
	RTMP(3)=RBDATA(3)+Z
	  DO I=4,8
	RTMP(I)=RBDATA(I)
	  ENDDO
C
	RETURN
C
C----- BODY =  "ELL" -----
9       RTMP(1)=RBDATA(1)+X
	RTMP(2)=RBDATA(2)+Y
	RTMP(3)=RBDATA(3)+Z
	  DO I=4,12
	RTMP(I)=RBDATA(I)
	  ENDDO
C
	RETURN
C
C----- END OF BODY DATA ---------------
C
	END




	SUBROUTINE GBROT(JREAD,NBT,RTMP,JTMP,f017In)
C****************************************************************
C* ROT card interpretation                                      *
C****************************************************************
C
	PRINT *,'SUBROUTINE GBROT is not ready yet - sorry'
	STOP'EXIT from GBROT: ROT card is not ready,use other!'
C
	END




	SUBROUTINE GFEND(JREAD,f017In)
C****************************************************************
C*  Jump to the end of file                                     *
C****************************************************************
       CHARACTER*80 JREAD(1000)
       INTEGER f017In
1       FORMAT(2X,A3)
2       f017In=f017In+1
	READ(JREAD(f017In),1,END=3)ITMP
	GOTO 2
3       RETURN
	END




	SUBROUTINE GBREAD(JREAD,NB,RBDATA,JBDATA,f017In)
C****************************************************************
C* MOV and ROT propotypes reading                               *
C****************************************************************
C
	IMPLICIT REAL*8 (A-H,R-Z)
C
	PARAMETER(MAXB=1000,MAXZ=1000)
	PARAMETER(JDIMBD=37*MAXB,JDIMZD=5*MAXZ)
	PARAMETER(MKINDB=10,MKINDS=10)
C
	COMMON/GDATA0/ NUMBOD,NUMZON
	COMMON/GDATA1/ ISYMB(MKINDB),NBSURF(MKINDB),NSPAR(MKINDS)
	COMMON/GDATA2/ BODYDA(MAXB),BODYDB(JDIMBD)
	INTEGER BODYDA
	COMMON/GDATA3/ ZONEDA(MAXZ),ZONEDB(JDIMZD)
	INTEGER ZONEDA,ZONEDB
	  CHARACTER*80 JREAD(1000)
      INTEGER f017In
C
	DIMENSION RBDATA(*)
	DIMENSION JBDATA(*)
C
	f017In=0 ! REWIND JREAD
C
	NBCUR=0
C
	f017In=f017In+1
	READ(JREAD(f017In),20)JDBG1,JDBG2
20      FORMAT(2I5)
C
97      FORMAT(2X,A3,1X,I4,6D10.3)
98      FORMAT(10X,6D10.3)
99      FORMAT(10X,6(I9,1X))
C
100     IF(NB.EQ.NBCUR)THEN
	  RETURN
	ELSE
	  NBCUR=NBCUR+1
	ENDIF
C
	f017In=f017In+1
	READ(JREAD(f017In),97) ITYPE,NBFIL,(RBDATA(I),I=1,6)
C
101       DO IBT=1,10
	IF( ISYMB(IBT).EQ.ITYPE ) GOTO 102
	  END DO
C
	 STOP'EXIT from GBREAD : invalid body type !'
C
102     GOTO ( 1, 2, 3, 4, 5, 6, 7, 8, 9),IBT
C
C----- BODY =  "SPH" -----
1        GOTO 100
C
C----- BODY =  "WED" -----
2       f017In=f017In+1
	READ(JREAD(f017In),98) (RBDATA(I),I=7,12)
C
	GOTO 100
C
C----- BODY =  "ARB" -----
3       f017In=f017In+1
	READ(JREAD(f017In),98) (RBDATA(I),I=7,12)
	f017In=f017In+1
	READ(JREAD(f017In),98) (RBDATA(I),I=13,18)
	f017In=f017In+1
	READ(JREAD(f017In),98) (RBDATA(I),I=19,24)
	f017In=f017In+1
	READ(JREAD(f017In),99) (JBDATA(I),I=1,6)
C
	GOTO 100
C
C----- BODY =  "BOX" -----
4       f017In=f017In+1
	READ(JREAD(f017In),98) (RBDATA(I),I=7,12)
C
	GOTO 100
C
C----- BODY =  "RPP" -----
5       GOTO 100
C
C----- BODY =  "RCC" -----
6       f017In=f017In+1
	READ(JREAD(f017In),98) RBDATA(7)
C
	GOTO 100
C
C----- BODY =  "REC" -----
7       f017In=f017In+1
	READ(JREAD(f017In),98) (RBDATA(I),I=7,12)
C
	GOTO 100
C
C----- BODY =  "TRC" -----
8       f017In=f017In+1
	READ(JREAD(f017In),98) (RBDATA(I),I=7,8)
C
	GOTO 100
C
C----- BODY =  "ELL" -----
9       f017In=f017In+1
	READ(JREAD(f017In),98) (RBDATA(I),I=7,12)
C
	GOTO 100
C
C----- END OF BODY DATA ---------------
C
	END




	SUBROUTINE GBWRT(JOUT,IBT,RTMP,JTMP,f017Out)
C****************************************************************
C* Output file creation (body data)                             *
C****************************************************************
	IMPLICIT REAL*8 (A-H,R-Z)
C
	PARAMETER(MAXB=1000,MAXZ=1000)
	PARAMETER(JDIMBD=37*MAXB,JDIMZD=5*MAXZ)
	PARAMETER(MKINDB=10,MKINDS=10)
	  CHARACTER*80 JOUT(1000)
	  INTEGER f017Out
C
	COMMON/GDATA0/ NUMBOD,NUMZON
	COMMON/GDATA1/ ISYMB(MKINDB),NBSURF(MKINDB),NSPAR(MKINDS)
	COMMON/GDATA2/ BODYDA(MAXB),BODYDB(JDIMBD)
	INTEGER BODYDA
	COMMON/GDATA3/ ZONEDA(MAXZ),ZONEDB(JDIMZD)
	INTEGER ZONEDA,ZONEDB
	COMMON/GDATA4/ JFLSUF(3,MAXB),BTRZLT(MAXB),
     *                 DIGINF,DIGMIN,DELTA,STPLIM
	COMMON/GDATAP/ X,Y,Z,CX,CY,CZ,XN,YN,ZN,CXN,CYN,CZN,
     *                 NZONO,MEDOLD,NBPO,NSO,
     *                 NZONC,MEDCUR,NBPC,NSCI,NSCO,
     *                 PINSFL,PARTZD
	INTEGER PINSFL,PARTZD
C
	DIMENSION RTMP(30)    ! Sobol: (30) instead of (*) to printout
	DIMENSION JTMP(30)    ! Sobol: (30) instead of (*) to printout
C
	IF(IBT.LT.0)THEN
	f017Out=f017Out+1
	WRITE(JOUT(f017Out),20)'END'
20      FORMAT(2X,A3,70X)
	RETURN
	END IF
C
	IF(IBT.EQ.0)THEN
	JDBG1=0
	JDBG2=0
	f017Out=f017Out+1
	WRITE(JOUT(f017Out),21)JDBG1,JDBG2,(JTMP(I),I=1,15)
21      FORMAT(2I5,10X,15A4)
	RETURN
	END IF
C
97      FORMAT(2X,A3,1X,I4,6F10.3)
98      FORMAT(10X,6F10.3)
99      FORMAT(10X,6(I9,1X))
C
	f017Out=f017Out+1
	WRITE(JOUT(f017Out),97) ISYMB(IBT),NUMBOD,(RTMP(I),I=1,6)
C
	GOTO ( 1, 2, 3, 4, 5, 6, 7, 8, 9),IBT
C
C----- BODY =  "SPH" -----
1        RETURN
C
C----- BODY =  "WED" -----
2       f017Out=f017Out+1
	WRITE(JOUT(f017Out),98) (RTMP(I),I=7,12)
C
	RETURN
C
C----- BODY =  "ARB" -----
3       f017Out=f017Out+1
	WRITE(JOUT(f017Out),98) (RTMP(I),I=7,12)
	f017Out=f017Out+1
	WRITE(JOUT(f017Out),98) (RTMP(I),I=13,18)
	f017Out=f017Out+1
	WRITE(JOUT(f017Out),98) (RTMP(I),I=19,24)
	f017Out=f017Out+1
	WRITE(JOUT(f017Out),99) (JTMP(I),I=1,6)
C
	RETURN
C
C----- BODY =  "BOX" -----
4       f017Out=f017Out+1
	WRITE(JOUT(f017Out),98) (RTMP(I),I=7,12)
C
	RETURN
C
C----- BODY =  "RPP" -----
5       RETURN
C
C----- BODY =  "RCC" -----
6       f017Out=f017Out+1
	WRITE(JOUT(f017Out),98) RTMP(7)
C
	RETURN
C
C----- BODY =  "REC" -----
7       f017Out=f017Out+1
	WRITE(JOUT(f017Out),98) (RTMP(I),I=7,12)
C
	RETURN
C
C----- BODY =  "TRC" -----
8       f017Out=f017Out+1
	WRITE(JOUT(f017Out),98) (RTMP(I),I=7,8)
C
	RETURN
C
C----- BODY =  "ELL" -----
9       f017Out=f017Out+1
	WRITE(JOUT(f017Out),98) (RTMP(I),I=7,12)
C
	RETURN
C
C----- END OF BODY DATA ---------------
C
	END





	SUBROUTINE GZWRT(JOUT,NNN,IALP,NAZ,IZSYM,JTMP,JTY1,JTY2,f017Out)
C****************************************************************
C* Output file creation (zone data)                             *
C****************************************************************
	IMPLICIT REAL*8 (A-H,R-Z)
C
	PARAMETER(MAXB=1000,MAXZ=1000)
	PARAMETER(JDIMBD=37*MAXB,JDIMZD=5*MAXZ)
	PARAMETER(MKINDB=10,MKINDS=10)
C
	COMMON/GDATA0/ NUMBOD,NUMZON
	COMMON/GDATA1/ ISYMB(MKINDB),NBSURF(MKINDB),NSPAR(MKINDS)
	COMMON/GDATA2/ BODYDA(MAXB),BODYDB(JDIMBD)
	INTEGER BODYDA
	COMMON/GDATA3/ ZONEDA(MAXZ),ZONEDB(JDIMZD)
	INTEGER ZONEDA,ZONEDB
	COMMON/GDATA4/ JFLSUF(3,MAXB),BTRZLT(MAXB),
     *                 DIGINF,DIGMIN,DELTA,STPLIM
	COMMON/GDATAP/ X,Y,Z,CX,CY,CZ,XN,YN,ZN,CXN,CYN,CZN,
     *                 NZONO,MEDOLD,NBPO,NSO,
     *                 NZONC,MEDCUR,NBPC,NSCI,NSCO,
     *                 PINSFL,PARTZD
	INTEGER PINSFL,PARTZD
	  CHARACTER*80 JOUT(1000)
      INTEGER f017Out
C
	DATA IOR,IBL,IEND/2HOR,3H   ,3HEND/
	DIMENSION IZSYM(*)
	DIMENSION JTMP(*)
	DIMENSION JTY1(*),JTY2(*)
C
      IF(NNN.EQ.1)THEN
C
1       FORMAT(2X,A3,I5,9(A2,I5))
2       FORMAT(2X,A3,5X,9(A2,I5))
	IF(IALP.EQ.IBL)THEN
	 f017Out=f017Out+1
	WRITE(JOUT(f017Out),1)IALP,(IZSYM(I),JTMP(I),I=1,9)
	ELSE
	 f017Out=f017Out+1
	WRITE(JOUT(f017Out),2)IALP,NAZ,(IZSYM(I),JTMP(I),I=1,9)
	ENDIF
C
	ELSE
C
3       FORMAT(14I5)
	DO I=1,L
	 f017Out=f017Out+1
	WRITE(JOUT(f017Out),3)(JTY1(J+(I-1)*14),J=1,14)
	END DO
C
	DO I=1,L
	 f017Out=f017Out+1
	WRITE(JOUT(f017Out),3)(JTY2(J+(I-1)*14),J=1,14)
	END DO
C
	ENDIF
C
	END




	FUNCTION GZMED(IZONE)
C****************************************************************
C* Determination medium for zone number IZONE                   *
C****************************************************************
C
	PARAMETER(MAXB=1000,MAXZ=1000)
	PARAMETER(JDIMBD=37*MAXB,JDIMZD=5*MAXZ)
C
	COMMON/GDATA3/ ZONEDA(MAXZ),ZONEDB(JDIMZD)
	INTEGER ZONEDA,ZONEDB
C
	JZADDR=ZONEDA(IZONE)
	GZMED=ZONEDB(JZADDR+1)
C
	END



	SUBROUTINE GCLCER
C****************************************************************
C* Debug: data of common print                                  *
C****************************************************************
	IMPLICIT REAL*8 (A-H,R-Z)
C
	PARAMETER(MAXB=1000,MAXZ=1000)
	PARAMETER(JDIMBD=37*MAXB,JDIMZD=5*MAXZ)
	PARAMETER(MKINDB=10,MKINDS=10)
C
	COMMON/GDATA0/ NUMBOD,NUMZON
	COMMON/GDATA1/ ISYMB(MKINDB),NBSURF(MKINDB),NSPAR(MKINDS)
	COMMON/GDATA2/ BODYDA(MAXB),BODYDB(JDIMBD)
	INTEGER BODYDA
	COMMON/GDATA3/ ZONEDA(MAXZ),ZONEDB(JDIMZD)
	INTEGER ZONEDA,ZONEDB
	COMMON/GDATA4/ JFLSUF(3,MAXB),BTRZLT(MAXB),
     *                 DIGINF,DIGMIN,DELTA,STPLIM
	COMMON/GDATAP/ X,Y,Z,CX,CY,CZ,XN,YN,ZN,CXN,CYN,CZN,
     *                 NZONO,MEDOLD,NBPO,NSO,
     *                 NZONC,MEDCUR,NBPC,NSCI,NSCO,
     *                 PINSFL,PARTZD
	INTEGER PINSFL,PARTZD
C
	PRINT *,'Call GCLCER:'
	PRINT *,'                  '
	PRINT *,'NUMBOD,NUMZON=',NUMBOD,NUMZON
C
	PRINT *,'                  '
	PRINT *,'xyz=',X,Y,Z
	PRINT *,'cos=',CX,CY,CZ
	PRINT *,'xyzN=',XN,YN,ZN
	PRINT *,'cosN=',CXN,CYN,CZN
	PRINT *,'                  '
C
	PRINT *,' I   JFLSUF(1,I) JFLSUF(2,I) JFLSUF(3,I) BTRZLT(I)'
	DO I=1,NUMBOD
	PRINT *,I,JFLSUF(1,I),JFLSUF(2,I),JFLSUF(3,I),BTRZLT(I)
	ENDDO
C
	PRINT *,'                  '
	PRINT *,'nzono,nbpo,nso      =',NZONO,NBPO,NSO
	PRINT *,'nzonc,nbpc,nsci,nsco=',NZONC,NBPC,NSCI,NSCO
C
	PRINT *,'                  '
	PRINT *,'pinsfl=',PINSFL
C
	CALL GUSERI
C
	RETURN
	END



	SUBROUTINE GUSERI
C
	PRINT *,'Call GUSERI:'
	PRINT *,'Write your SUBROUTINE GUSERI if need.'
C
	RETURN
	END
