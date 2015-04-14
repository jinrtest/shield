***********************************************************************
*                                                                     *
*        The following subroutines are written by N.M.Sobolevsky,     *
*        Institute for Nuclear Research RAS, Moscow,                  *
*        E-mail: Sobolevs@AL20.INR.Troitsk.RU                         *
*                                                                     *
* The modifications of 16.10.99 to take into account TIME of life of  *
* neutron and transport inside THERMAL enegry group are labelled with:*
*                                                      ! TIME         *
*                                                      ! THERMAL      *
*                                       Sobol 16.10.99                *
***********************************************************************
C
      SUBROUTINE NEUTLO                                              ! GEMCA
C MAIN SUBROUTINE FOR LOENT. IT TAKES THE NEUTRON FROM THE SOURCE ARRAY
C AND CALLS LOENT.
      COMMON /SOURN/ SOURN(13,100000),ISOURN(7,100000),LSN13,LSN7,LSN100   ! GSITR
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      COMMON /WNPN/ WN(11),WN0(11),WNMEM(11),WNI(11),WNF(11),NGRUP   ! GEMCA
C
C NSOR - NUMBER OF NEUTRONS IN SOURCE ARRAY  SOURN
      NSOR=ISOURN(6,1)
          IF(NSOR.EQ.0)RETURN
C
      DO 1 J=2,NSOR+1
          DO 2 K=1,11                                                ! GEMCA
    2     WN(K)=SOURN(K,J)
      CALL LOENT(LABNOR)                                             ! GEMCA
          IF(LABNOR.NE.0)THEN
              CONTINUE!WRITE(25,1001)LABNOR,IXINIT,NUMTRE
 1001         FORMAT(' NEUTLO DIAGNOSTIC: LABNOR=',I2,', IXINIT=',I12,
     *               ', NUMTRE=',I6)
              GO TO 1
          ELSE
              GO TO 1
          END IF
    1 CONTINUE
C
      RETURN
      END



      SUBROUTINE LNMACR
C CALCULATION OF 28-GROUP MACROSCOPICAL TOTAL AND TRANSPORT CROSS SECTIONS
C SIMALN(I,J,K): I=-1,0,1-26 - NO. OF GROUP
C                    J=1 - MACROSCOPICAL TOTAL CROSS SECTION
C                    J=2 - MACROSCOPICAL TRANSPORT CROSS SECTION
C                        K=1-48 - NO. OF MEDIUM  ! CORNEL
C--------------------------------------------------------------------------
      COMMON /MACRIN/ ELMIX(8,24,48),RHOMED(48),NELEMD(2,48),NUMMED  ! CORNEL
      COMMON /BNAB/ EGR(-1:27),SPFU(10,-1:13),VNU(10),SGLN(-1:26,7,110),
     *              SUV(13,-1:12,110)                             ! BNAB28
      COMMON /SIMALN/ SIMALN(-1:26,2,48)      ! CORNEL            ! BNAB28
      COMMON /RENUMB/ IREN(110)
C Modification of 13.04.94
          COMMON /OLN/ OLN
C-------------------------
      COMMON /INFEL/ ZNUC(110),ATWEI(110),DENS(110),RIEV(110),SYMB(110)
      CHARACTER*2 SYMB
c      DIMENSION ITEMP(16)  !  CORNEL
C
C CHECK OF PRESENCE OF NEUTRON DATA IN ARRAY SGLN FOR ALL NUCLIDS IN THE TASK
      DO 1 K=1,NUMMED
      MEDIA=K
      NELEM=NELEMD(2,MEDIA)
          DO 2 J=1,NELEM
          NUCLID=IFIX(ELMIX(1,J,MEDIA))
          IRN=IREN(NUCLID)
              IF(IRN.LE.-1)THEN
c                 COMMENTED BY A.TIMOFEEV
                 CONTINUE!WRITE(25,1001)SYMB(NUCLID),IFIX(ZNUC(NUCLID)),
!     *                                      NINT(ATWEI(NUCLID))
1001             FORMAT(' SORRY, THERE IS NO NEUTRON DATA FOR',
     *                  ' NUCLEUS ',A2,'(',I2,',',I3,') IN ARRAY SGLN')
C---------------- Modification of 13.04.94 ------------------------
C                 STOP
!                  IF(OLN.LT.14.5)THEN                         ! BNAB28
!                      STOP
!                  ELSE
!C                      CONTINUE!WRITE(27,1001)SYMB(NUCLID),IFIX(ZNUC(NUCLID)),
!C     *                                           NINT(ATWEI(NUCLID))
!                  END IF
                      WRITE(6,1001)SYMB(NUCLID),IFIX(ZNUC(NUCLID)),
     *                                           NINT(ATWEI(NUCLID))                  
C------------------------------------------------------------------
              ELSE
                  DO 3 I=-1,26                                ! BNAB28
                  SGLN(I,7,IRN)=SGLN(I,1,IRN)+SGLN(I,3,IRN)+
     *                          SGLN(I,4,IRN)+SGLN(I,5,IRN)
    3             CONTINUE
c                  CALL PRBNAB(NUCLID)
                  GO TO 2
              END IF
    2     CONTINUE
    1 CONTINUE
C
C CALCULATION OF 28-GROUP MACROSCOPICAL TOTAL AND TRANSPORT CROSS SECTION
      DO 11 K=1,NUMMED
      MEDIA=K
      NELEM=NELEMD(2,MEDIA)
          DO 12 J=-1,26                                       ! BNAB28
              SMTOT=0.
              SMTRAN=0.
              DO 13 I=1,NELEM
              NUCLID=IFIX(ELMIX(1,I,MEDIA))
                  IRN=IREN(NUCLID)
C Modification 13.04.94-----------------
                  IF(IRN.LE.-1) GO TO 13
C --------------------------------------
              CONC=ELMIX(2,I,MEDIA)
              SIGFIS=SGLN(J,1,IRN)
              SIGCAP=SGLN(J,3,IRN)
              SIGIN =SGLN(J,4,IRN)
              SIGEL =SGLN(J,5,IRN)
              COSLAB=SGLN(J,6,IRN)
              SMTOT=SMTOT+CONC*(SIGFIS+SIGCAP+SIGIN+SIGEL)
              SMTRAN=SMTRAN+CONC*(SIGFIS+SIGCAP+SIGIN+SIGEL*(1.-COSLAB))
   13         CONTINUE
          SIMALN(J,1,MEDIA)=SMTOT*1000.
          SIMALN(J,2,MEDIA)=SMTRAN*1000.
   12     CONTINUE
   11 CONTINUE
C
C PRINTING OF 28-GROUP MACROSCOPIC CROSS SECTIONS
C      CONTINUE!WRITE(27,20)
C   20 FORMAT(/20X,' 28-GROUP MACROSCOPIC TOTAL NEUTRON CROSS',    ! BNAB28
C     *       ' SECTIONS (1/CM)')
C
C ---------------------- CORNEL ----------------------------
C Printout for MEDIA Nos. 1-12
C      CONTINUE!WRITE(27,21)
C   21 FORMAT(' GR','   Med1   Med2   Med3   Med4   Med5   Med6',
C     *             '   Med7   Med8   Med9  Med10  Med11  Med12')
C      DO J=-1,26                                               ! BNAB28
C        CONTINUE!WRITE(27,30)J,(SIMALN(J,1,K),K=1,12)
C   30   FORMAT(I3,12F7.4)
C      END DO
C
C Printout for MEDIA Nos.13-24, if any
C                   IF(NUMMED.gt.12)THEN
C      CONTINUE!WRITE(27,22)
C   22 FORMAT(/' GR','  Med13  Med14  Med15  Med16  Med17  Med18',
C     *              '  Med19  Med20  Med21  Med22  Med23  Med24')
C      DO J=-1,26                                               ! BNAB28
C        CONTINUE!WRITE(27,30)J,(SIMALN(J,1,K),K=13,24)
C      END DO
C                   END IF
C
C Printout for MEDIA Nos.25-36, if any
C                   IF(NUMMED.gt.24)THEN
C      CONTINUE!WRITE(27,23)
C   23 FORMAT(/' GR','  Med25  Med26  Med27  Med28  Med29  Med30',
C     *              '  Med31  Med32  Med33  Med34  Med35  Med36')
C      DO J=-1,26                                               ! BNAB28
C        CONTINUE!WRITE(27,30)J,(SIMALN(J,1,K),K=25,36)
C      END DO
C                   END IF
C
C Printout for MEDIA Nos.37-48, if any
C                   IF(NUMMED.gt.36)THEN
C      CONTINUE!WRITE(27,24)
C   24 FORMAT(/' GR','  Med37  Med38  Med39  Med40  Med41  Med42',
C     *              '  Med43  Med44  Med45  Med46  Med47  Med48')
C      DO J=-1,26                                               ! BNAB28
C        CONTINUE!WRITE(27,30)J,(SIMALN(J,1,K),K=37,48)
C      END DO
C                   END IF
C
c      CONTINUE!WRITE(27,25)
c   25 FORMAT(/36X,' 28-GROUP MACROSCOPIC TRANSPORT NEUTRON CROSS',  ! BNAB28
c     *       ' SECTIONS (1/CM)')
c      CONTINUE!WRITE(27,22)ITEMP
c      DO 26 J=-1,26                                               ! BNAB28
c   26 CONTINUE!WRITE(27,24)J,(SIMALN(J,2,K),K=1,16)
C ---------------------- end CORNEL ------------------------
C
      RETURN
      END



      SUBROUTINE PRBNAB(NUCLID)
C PRINTING OF 28-GROUP BNAB NEUTRON CONSTANTS FOR NUCLEUS
C NO. IRN=IREN(NUCLID)
      COMMON /BNAB/ EGR(-1:27),SPFU(10,-1:13),VNU(10),SGLN(-1:26,7,110),
     *              SUV(13,-1:12,110)                             ! BNAB28
      COMMON /RENUMB/ IREN(110)
      COMMON /INFEL/ ZNUC(110),ATWEI(110),DENS(110),RIEV(110),SYMB(110)
      CHARACTER*2 SYMB
C
      IRN=IREN(NUCLID)
          IZ=IFIX(ZNUC(NUCLID))
          AT=ATWEI(NUCLID)
C
      CONTINUE!WRITE(27,10)SYMB(NUCLID),IZ,AT
   10 FORMAT(/40X,'BNAB 28-GROUP NEUTRON CONSTANTS FOR ',         ! BNAB28
     *        A2,'(',I3,',',F9.5,')'/)
      CONTINUE!WRITE(27,20)
   20 FORMAT(50X,'MAIN GROUP CONSTANTS')
      CONTINUE!WRITE(27,30)
   30 FORMAT(' GR.',5X,'ENERGY (MEV)',8X,'S-TOT',6X,'S-FIS',5X,
     *       'AVER NU',5X,'S-CAP',6X,'S-IN',7X,'S-EL',6X,'COS EL')
      DO 40 J=-1,25                                               ! BNAB28
   40 CONTINUE!WRITE(27,50)J,EGR(J+1),EGR(J),SGLN(J,7,IRN),(SGLN(J,K,IRN),K=1,6)
   50 FORMAT(1X,I2,2X,E9.3,'-',E9.3,7(2X,F9.4))
      CONTINUE!WRITE(27,60)SGLN(26,7,IRN),(SGLN(26,K,IRN),K=1,6)
   60 FORMAT(' 26',8X,'THERMAL',6X,7(2X,F9.4))
C
      CONTINUE!WRITE(27,70)
   70 FORMAT(/50X,'SIGMA IN (I,I+K)  AT K EQUAL:')
      CONTINUE!WRITE(27,80)
   80 FORMAT('  I',3X,'K= 0',5X,'K= 1',5X,'K= 2',5X,'K= 3',5X,'K= 4',
     *             5X,'K= 5',5X,'K= 6',5X,'K= 7',5X,'K= 8',5X,'K= 9',
     *             5X,'K=10',5X,'K=11',5X,'SUMMA')
C
      DO 90 J=-1,12                                               ! BNAB28
   90 CONTINUE!WRITE(27,100)J,(SUV(K,J,IRN),K=1,13)
  100 FORMAT(1X,I2,13(2X,F7.4))
      CONTINUE!WRITE(27,110)
  110 FORMAT(1X,10('------------'))
C
      RETURN
      END



      SUBROUTINE LOENT(LABNOR)                                       ! GEMCA
C MONTE CARLO NEUTRON TRASPORT CODE. 28-GROUP BNAB NEUTRON CONSTANTS.
C LINEAR-ANIZOTROPIC APPROXIMATION OF ELASTIC SCATTERING.
C INPUT:
C       NEUTRON IN ARRAY WN(11):                                     ! GEMCA
C                        WN(1),WN(2),WN(3) - X0,Y0,Z0 (CM)
C                        WN(4),WN(5),WN(6) - COS(THETA),SIN(PHI),COS(PHI)
C                        WN(7) - KINETIC ENERGY (MEV)
C                        WN(8) - WEIGHT
C                        WN(9) - CURRENT GEOM. ZONE                  ! GEMCA
C                        WN(10) - TIME OF BORN                       ! GEMCA
C                        WN(11) - CURRENT TIME                       ! GEMCA
C OUTPUT:
C        LABNOR.GT.0  IF ABNORMAL RETURN FROM LOENT OCCURES:
C                     LABNOR=1 - RANGE (CG) FAILURE   / not used     ! GEMCA
C                           =2 - NEUTRON MEMORY  REMLN  OVERFILLING
C                           =3 - OTHERS
C        LABNOR=0 IF NORMAL RETURN
C------------------------------------------------------------------------
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      COMMON /WNPN/ WN(11),WN0(11),WNMEM(11),WNI(11),WNF(11),NGRUP   ! GEMCA
      COMMON /REMLN/ REMLN(11,5000),NREMLN,LRL50   ! Koshka_ADS      ! GEMCA
      COMMON /MANLEN/ LNSTAT,LNMED,LNZON,LNCOD,LNUCLD,IRENLN,IPRLN,
     *                NFISN,NSECN
      COMMON /RENUMB/ IREN(110)
      COMMON /BNAB/ EGR(-1:27),SPFU(10,-1:13),VNU(10),SGLN(-1:26,7,110),
     *              SUV(13,-1:12,110)                             ! BNAB28
      COMMON /SPECLN/ SPEFIS(-1:13),SPECIN(12)                    ! BNAB28
      COMMON /OLN/ OLN
C
C
      LABNOR=0
      WN(11)=WN(10)   ! to give born time to current time         ! TIME
C
C REMEMBERING OF SOURCE NEUTRON
        DO J=1,11                                                    ! GEMCA
          WN0(J)=WN(J)                                               ! GEMCA
        ENDDO                                                        ! GEMCA
C
C HISTOGRAMMING OF INPUT NEUTRON SOURCE
      CALL HISTLN(8)
C
C CLEANING OF INNER NEUTRON MEMORY ARRAY  REMLN
        DO J=1,LRL50                                                 ! GEMCA
          DO K=1,11                                                  ! GEMCA
            REMLN(K,J)=0.                                            ! GEMCA
          ENDDO                                                      ! GEMCA
        ENDDO                                                        ! GEMCA
          NREMLN=0
C
C CHECK OF NEUTRON ENERGY (T<TMIN ?)
  111 IF(WN(7).LT.OLN)THEN
          CALL HISTLN(1)
          GO TO 222
      END IF
C
C DETERMINATION OF NEUTRON ENERGY GROUP NUMBER  NGRUP
      NGRUP=NGROUP(WN(7))
C
C SAMPLING OF NEUTRON RANGE
        DO J=1,11                                                    ! GEMCA
          WNI(J)=WN(J)                                               ! GEMCA
        ENDDO                                                        ! GEMCA
              CALL RANGLN(LFAILU)                                    ! GEMCA
      IF(LFAILU.EQ.0)THEN
        DO J=1,11                                                    ! GEMCA
          WNF(J)=WN(J)                                               ! GEMCA
        ENDDO                                                        ! GEMCA
          CALL HISTLN(2)
      ELSE
          LABNOR=1
          RETURN
      END IF
C
C DETERMINATION OF THE NEUTRON STATUS IN THE END POINT OF RANGE
      IF(LNSTAT.EQ.3)GO TO 20
      IF(LNSTAT.EQ.0)GO TO 30
      CONTINUE!WRITE(25,1001)LNSTAT,IXINIT,NUMTRE,WN0
 1001 FORMAT(' LOENT DIAGNOSTIC: UNRECOGNIZED PARTICLE STATUS AFTER',
     *       ' RANGE.'/19X,'LNSTAT=',I3,', IXINIT=',I12,', NUMTRE=',I6/
     *       ' SOURCE NEUTRON:   ',3F7.2,1X,3F7.4,2X,2F8.5/)
      LABNOR=3
      RETURN
C
C NEUTRON ABANDON THE CONFIGURATION
   20 CALL HISTLN(3)
cC--------------- TIME/THERMAL Debug -------------------
c      RRR=SQRT(WN(1)**2+WN(2)**2+WN(3)**2)
c      EEE=WN(7)
c      TTT=WN(11)
c      CONTINUE!WRITE(25,3456)NUMTRE,RRR,EEE,TTT
c 3456 FORMAT(' NUMTRE:',I4,'    R=',F6.2,'  E=',F11.8,'   TIME=',F9.3)
cC--------------- end TIME/THERMAL Debug -------------------
      GO TO 222
C
C SAMPLING OF THE NUCLID WITH WHICH NEUTRON INTERACT
   30 CALL PREPLN
C
C SAMPLING OF THE NEUTRON-NUCLEUS INTERACTION TYPE  IPRLN
      CALL TYPELN
      GO TO(41,52,63,74),IPRLN
C
C FISSION (N,F), IPRLN=1 ----------------------------------------------
   41 CONTINUE
C DETERMINATION OF FISSION NEUTRONS NUMBER  NFISN
      AVERN=SGLN(NGRUP,2,IRENLN)
      AVERI=AINT(AVERN)
          DAV=AVERN-AVERI
      IF(RRAN(IX).LE.DAV)THEN
          NFISN=IFIX(AVERI)+1
      ELSE
          NFISN=IFIX(AVERI)
      END IF
C
C CHECK OF NEUTRON MEMORY OVERFILLING
      IF(NREMLN+NFISN .GT. LRL50)THEN
          CONTINUE!WRITE(25,1002)NFISN,IXINIT,NUMTRE
 1002     FORMAT(' LOENT DIAGNOSTIC: NEUTRON MEMORY OVERFILLING, ',
     *           'NFISN=',I2,', IXINIT=',I12,', NUMTRE=',I6)
          LABNOR=2
          RETURN
      END IF
C
C PREPARATION OF FISSION NEUTRONS SPECTRUM
      CALL PREPFI(AVERN)
C
C SAMPLING OF FISSION NEUTRONS CHARACTERISTICS
      DO 42 J=1,NFISN
              JMEM=NREMLN+J
          REMLN(1,JMEM)=WN(1)
          REMLN(2,JMEM)=WN(2)
          REMLN(3,JMEM)=WN(3)
              CALL UNIFI(SF,CF)
          REMLN(4,JMEM)=2.0*RRAN(IX)-1.0
          REMLN(5,JMEM)=SF
          REMLN(6,JMEM)=CF
              REMLN(7,JMEM)=TNEUT(1)
              REMLN(8,JMEM)=WN(8)
          REMLN(9,JMEM)=WN(9)                                        ! GEMCA
          REMLN(10,JMEM)=WN(11)                                      ! GEMCA
          REMLN(11,JMEM)=WN(11)                                      ! GEMCA
   42 CONTINUE
      NREMLN=NREMLN+NFISN
      CALL HISTLN(4)
      GO TO 333
C
C CAPTURE (N,C), IPRLN=2 --------------------------------------------
   52 CONTINUE
      CALL HISTLN(5)
      GO TO 222
C
C INELASTIC SCATTERING (N,N), (N,2N), IPRLN=_3 -----------------------
   63 CONTINUE
C SAMPLING OF SECONDARY NEUTRONS NUMBER  NSECN
      SIGN2N=SUV(13,NGRUP,IRENLN)-SGLN(NGRUP,4,IRENLN)
          IF(SIGN2N.LE.0.)THEN
              NSECN=1
              GO TO 64
          ELSE
              IF(RRAN(IX).LE.(SIGN2N/SGLN(NGRUP,4,IRENLN)))THEN
                  NSECN=2
                  GO TO 64
              ELSE
                  NSECN=1
                  GO TO 64
              END IF
          END IF
C
C CHECK OF NEUTRON MEMORY OVERFILLING
   64 IF(NREMLN+NSECN .GT. LRL50)THEN
          CONTINUE!WRITE(25,1003)NSECN,IXINIT,NUMTRE
 1003     FORMAT(' LOENT DIAGNOSTIC: NEUTRON MEMORY OVERFILLING, ',
     *           'NSECN=',I2,', IXINIT=',I12,', NUMTRE=',I6)
          LABNOR=2
          RETURN
      END IF
C
C PREPARATION OF SECONDARY NEUTRON SPECTRUM
      CALL PREPIN
C
C SUMPLING OF INELASTIC SCATTERED NEUTRONS CHARACTERISTICS
  165 DO 65 J=1,NSECN         ! 2007-05-12
              JMEM=NREMLN+J
          REMLN(1,JMEM)=WN(1)
          REMLN(2,JMEM)=WN(2)
          REMLN(3,JMEM)=WN(3)
              CALL UNIFI(SF,CF)
          REMLN(4,JMEM)=2.0*RRAN(IX)-1.0
          REMLN(5,JMEM)=SF
          REMLN(6,JMEM)=CF
              REMLN(7,JMEM)=TNEUT(3)
              REMLN(8,JMEM)=WN(8)
          REMLN(9,JMEM)=WN(9)                                        ! GEMCA
          REMLN(10,JMEM)=WN(11)                                      ! GEMCA
          REMLN(11,JMEM)=WN(11)                                      ! GEMCA
   65 CONTINUE
C---------------------------- 2007-05-12 ------------------------------
C The sum of energies of two secondary neutrons in the (n,2n) reaction
C must be less than the energy of primary neutron. 
      if(NSECN.eq.2)then
        if(REMLN(7,NREMLN+1)+REMLN(7,NREMLN+2).ge.WN(7))then
c          CONTINUE!WRITE(25,311)
c  311     format(/' NUMTRE       WN(7) REMLN(7,-1)  REMLN(7,N)     DE')
c          CONTINUE!WRITE(25,312)NUMTRE,WN(7),REMLN(7,NREMLN+1),REMLN(7,NREMLN+2),
c     *                       (WN(7)-REMLN(7,NREMLN+1)-REMLN(7,NREMLN+2))
c  312     format(I7,4E12.5)
          goto 165
        end if
      end if
C---------------------------- end 2007-05-12 --------------------------
      NREMLN=NREMLN+NSECN
      CALL HISTLN(6)
      GO TO 333
C
C ELASTIC SCATTERING (N,N), IPRLN=4 ----------------------------------
   74 CONTINUE
      CALL ELSCLN(CT,SF,CF,TPRIM)
      CALL SUBROT(WN(4),WN(5),WN(6),CT,SF,CF,CTR,SFR,CFR,STR)
          WN(4)=CTR
          WN(5)=SFR
          WN(6)=CFR
      WN(7)=TPRIM
      CALL HISTLN(7)
      GO TO 111
C---------------------------------------------------------------------
C
C NEUTRON MEMORY EMPTINESS CONTROL
  222 IF(NREMLN.le.0)THEN           ! .le., NOT .eq.; Sobol, 15.03.98
          RETURN
      ELSE
          GO TO 333
      END IF
C
C TAKE THE NEXT NEUTRON FROM MEMORY
  333 DO 80 J=1,11                                                   ! GEMCA
      WN(J)=REMLN(J,NREMLN)
          REMLN(J,NREMLN)=0.
              WNMEM(J)=WN(J)
   80 CONTINUE
          NREMLN=NREMLN-1
      GO TO 111
C
      END



      SUBROUTINE RANGLN(LFAILU)                                      ! GEMCA
C SAMPLING OF LOW ENERGY NEUTRON RANGE
C OUTPUT: LFAILU>0 IF RANGE CALCULATION FAILURE OCCURES.
C           COMMON OF Geometry                                       ! GEMCA
        COMMON/GDATAP/ X,Y,Z,CX,CY,CZ,XN,YN,ZN,CXN,CYN,CZN,          ! GEMCA
     *                 NZONO,MEDOLD,NBPO,NSO,                        ! GEMCA
     *                 NZONC,MEDCUR,NBPC,NSCI,NSCO,                  ! GEMCA
     *                 PINSFL,IOLEG                                  ! GEMCA
        INTEGER PINSFL                                               ! GEMCA
        REAL*8 X,Y,Z,CX,CY,CZ,XN,YN,ZN,CXN,CYN,CZN                   ! GEMCA
C           COMMONS OF LOENT
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      COMMON /WNPN/ WN(11),WN0(11),WNMEM(11),WNI(11),WNF(11),NGRUP   ! GEMCA
      COMMON /MANLEN/ LNSTAT,LNMED,LNZON,LNCOD,LNUCLD,IRENLN,IPRLN,
     *                NFISN,NSECN
      COMMON /STEPLN/ STEPLN(10,500),ITEPLN(6,500),ISTEP,LNN50
      COMMON /SIMALN/ SIMALN(-1:26,2,48)      ! CORNEL            ! BNAB28
      REAL*8 R,S                                                     ! GEMCA
C
C
      LFAILU=0
C
C CLEANING OF THE NEUTRON RANGE STEPS ARRAYS  STEPLN,ITEPLN.
      DO 1 J=1,LNN50
          DO 2 I=1,10
    2     STEPLN(I,J)=0.
          DO 3 K=1,6
    3     ITEPLN(K,J)=0
    1 CONTINUE
C
C NEUTRON CHARACTERISTICS BEFORE RANGE
      T0=WN(7)
          X0=WN(1)
          Y0=WN(2)
          Z0=WN(3)
              COST=WN(4)
              SINT=SQRT(1.0-COST**2)
              SINF=WN(5)
              COSF=WN(6)
        IRFIX=WN(9)                                                  ! GEMCA
        NEXTZ=WN(9)                                                  ! GEMCA
C
C DETERMINATION OF INPUT ZONE NUMBER  IR  AND MEDIUM NUMBER  MED
C FOR POINT X0,Y0,Z0.
      X=DBLE(X0)                                                     ! GEMCA
      Y=DBLE(Y0)                                                     ! GEMCA
      Z=DBLE(Z0)                                                     ! GEMCA
      CX=DBLE(SINT*COSF)                                             ! GEMCA
      CY=DBLE(SINT*SINF)                                             ! GEMCA
      CZ=DBLE(COST)                                                  ! GEMCA
C
C Dementyev -------------------
  111 DUMMY=RRAN(IX)
      IF(DUMMY.LE.0.0)GO TO 111
      GAMMA=-ALOG(DUMMY)
C------------------------------
          ISTEP=0
              XI=X0
              YI=Y0
              ZI=Z0
        GO TO 27                                                     ! GEMCA
C
C GNEXTZ give NEXT zone                                              ! GEMCA
C GCURZL determinate the pass lenght  S  into input zone No IZONE,   ! GEMCA
   26 CALL GNEXTZ(NEXTZ)                                             ! GEMCA
C
        IF(NEXTZ.EQ.0 .OR. MEDCUR.EQ.0)THEN                          ! GEMCA
C neutron exit into external vacuum                                  ! GEMCA
              LNSTAT=3
                  LNMED=MEDCUR                                       ! GEMCA
                  LNZON=NEXTZ                                        ! GEMCA
                  LNCOD=IRFIX                                        ! GEMCA
                      WN(1)=SNGL(X)                                  ! GEMCA
                      WN(2)=SNGL(Y)                                  ! GEMCA
                      WN(3)=SNGL(Z)                                  ! GEMCA
                      WN(9)=NEXTZ                                    ! GEMCA
          CALL STETLENPA(WN(8))               ! TLEST
          RETURN
        ELSE
          CONTINUE
C neutron crosses some internal surface (sobol, 02.03.98)
        END IF
C
   27 CALL GCURZL(NEXTZ,S)                                           ! GEMCA
C
      ISTEP=ISTEP+1
C
      IF(MEDCUR.EQ.1000)THEN                                         ! GEMCA
C neutron cross internal vacuum on step No.ISTEP                     ! GEMCA
          R=S                                                        ! GEMCA
            KSTATE=3
                SMTOT=0.
          GO TO 21
      END IF
C
C DETERMINATION OF THE NEUTRON STATUS AS A RESULT OF PASSING
C THROW THE LAYER OF MEDIUM No.MEDCUR (NO.MED ! cg) DEPTH S.
      SMTOT=SIMALN(NGRUP,1,MEDCUR)                                   ! GEMCA
      R=DBLE(GAMMA/SMTOT)                                            ! GEMCA
          IF(R.GT.S)THEN
              R=S                                                    ! GEMCA
                KSTATE=3
                    GAMMA=GAMMA-SMTOT*SNGL(S)                        ! GEMCA
          ELSE
              KSTATE=0
          END IF
C
C COORDINATES OF THE PASS END POINT.
   21  X=X+R*CX                                                      ! GEMCA
       Y=Y+R*CY                                                      ! GEMCA
       Z=Z+R*CZ                                                      ! GEMCA
C Input time culculation !                                           ! GEMCA
C      TIME=SHTIME()                                                 ! GEMCA
      TIMEmk=(SNGL(R)/30000.0)*SQRT(940.0/(2.0*WN(7)))           ! TIME
      WN(11)=WN(11)+TIMEmk                                       ! TIME
C
C FORMATION OF NEUTRON RANGE CHARACTERISTICS ARRAYS ON STEP NO.ISTEP.
C     CHECK OF STEPLN OVERFILLING
      IF(ISTEP.GT.LNN50)THEN
          CONTINUE!WRITE(25,1002)ISTEP,NUMTRE,IXINIT
 1002     FORMAT(' RANGLN DIAGNOSTIC: STEPLN OVERFILLING.',
     *           '  ISTEP=',I3,', NUMTRE=',I6,', IXINIT=',I12)
          LFAILU=1
          RETURN
      END IF
C
      STEPLN( 1,ISTEP)=T0
      STEPLN( 2,ISTEP)=SMTOT
      STEPLN( 3,ISTEP)=SNGL(S)
      STEPLN( 4,ISTEP)=SNGL(R)                                       ! GEMCA
      STEPLN( 5,ISTEP)=XI
      STEPLN( 6,ISTEP)=YI
      STEPLN( 7,ISTEP)=ZI
      STEPLN( 8,ISTEP)=SNGL(X)                                       ! GEMCA
      STEPLN( 9,ISTEP)=SNGL(Y)                                       ! GEMCA
      STEPLN(10,ISTEP)=SNGL(Z)                                       ! GEMCA
              ITEPLN(1,ISTEP)=ISTEP
              ITEPLN(2,ISTEP)=NGRUP
              ITEPLN(3,ISTEP)=MEDCUR                                 ! GEMCA
              ITEPLN(4,ISTEP)=NEXTZ                                  ! GEMCA
              ITEPLN(5,ISTEP)=NBPC                                   ! GEMCA
              ITEPLN(6,ISTEP)=NSCO                                   ! GEMCA
C
C
      IF(KSTATE.EQ.3)GOTO 25                                         ! GEMCA
          IF(KSTATE.EQ.0)THEN
              LNSTAT=0
                  LNMED=MEDCUR                                       ! GEMCA
                  LNZON=NEXTZ                                        ! GEMCA
                  LNCOD=IRFIX                                        ! GEMCA
                      WN(1)=SNGL(X)                                  ! GEMCA
                      WN(2)=SNGL(Y)                                  ! GEMCA
                      WN(3)=SNGL(Z)                                  ! GEMCA
                      WN(9)=NEXTZ                                    ! GEMCA
C                      WN(11)=TIME                                   ! GEMCA
              CALL STETLENPA(WN(8))               ! TLEST
              RETURN
          ELSE
              CONTINUE!WRITE(25,1004)KSTATE,MEDCUR,NUMTRE,IXINIT              ! GEMCA
 1004         FORMAT(' RANGLN DIAGNOSTIC: UNDEFINED NEUTRON STATUS', ! GEMCA
     *               ' AFTER RANGE.'/20X,'KSTATE=',I3,', MEDCUR=',   ! GEMCA
     *               I3,', NUMTRE=',I6,', IXINIT=',I12)              ! GEMCA
                  LFAILU=1                                           ! GEMCA
                  RETURN                                             ! GEMCA
          END IF
C
C neutron cross the medium No.MEDCUR without interaction and entry   ! GEMCA
C into the next input zone.                                          ! GEMCA
   25 IRFIX=NEXTZ                                                    ! GEMCA
          XI=SNGL(X)                                                 ! GEMCA
          YI=SNGL(Y)                                                 ! GEMCA
          ZI=SNGL(Z)                                                 ! GEMCA
      GO TO 26
C
      END



C      INTEGER FUNCTION NGROUP(E)
CC DETERMINATION OF NEUTRON ENERGY GROUP NUMBER IN BNAB 26-GROUP
CC ENERGY MESH.
C      COMMON /BNAB/ EGR(-1:27),SPFU(10,-1:13),VNU(10),SGLN(-1:26,7,110),
C     *              SUV(13,-1:12,110)                             ! BNAB28
CC
C      N1=1        ! It is impossible to modify this FUNCTION for BNAB28
C      N2=27
CC
C    1 NM=(N1+N2)/2
C          IF(E .GT. EGR(NM))THEN
C              N2=NM
C          ELSE
C              N1=NM
C          END IF
C              IF(N1+1 .EQ. N2)THEN
C                  NGROUP=N1
C                  RETURN
C              ELSE
C                  GO TO 1
C              END IF
C      END



      INTEGER FUNCTION NGROUP(E)
C DETERMINATION OF NEUTRON ENERGY GROUP NUMBER IN BNAB 28-GROUP
C ENERGY MESH.
      COMMON /BNAB/ EGR(-1:27),SPFU(10,-1:13),VNU(10),SGLN(-1:26,7,110),
     *              SUV(13,-1:12,110)                             ! BNAB28
C
      DO 1 I=0,27                                                 ! BNAB28
          IF(E .LE. EGR(I))GO TO 1
          NGROUP=I-1
          GO TO 2
    1 CONTINUE
      NGROUP=26                                                   ! THERMAL
    2 RETURN
      END



      SUBROUTINE PREPLN
C SAMPLING OF THE NUCLEUS-TARGET NUMBER  NUCLID  WITH WHICH
C LOW ENERGY NEUTRON INTERACTS.
C INPUT:
C       LNMED - INTERACTION POINT MEDIUM NUMBER
C       NGRUP - NO. OF NEUTRON ENERGY GROUP
C OUTPUT:
C        LNUCLD - NO OF NUCLID, WITH WHICH LOW ENERGY NEUTRON INTERACT
C                 (LNUCLD=NUCLID)
C        IRENLN - NO. OF NUCLEUS IN BNAB: IRENLN=IREN(NUCLID)
C--------------------------------------------------------------------
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      COMMON /WNPN/ WN(11),WN0(11),WNMEM(11),WNI(11),WNF(11),NGRUP   ! GEMCA
      COMMON /MANLEN/ LNSTAT,LNMED,LNZON,LNCOD,LNUCLD,IRENLN,IPRLN,
     *                NFISN,NSECN
      COMMON /RENUMB/ IREN(110)
      COMMON /BNAB/ EGR(-1:27),SPFU(10,-1:13),VNU(10),SGLN(-1:26,7,110),
     *              SUV(13,-1:12,110)                             ! BNAB28
      COMMON /MACRIN/ ELMIX(8,24,48),RHOMED(48),NELEMD(2,48),NUMMED  ! CORNEL
      DIMENSION TEMP(24)  ! CORNEL
C
C SAMPLING OF THE NUCLEUS-TARGET NUMBER NUCLID
      MEDIA=LNMED
      NELEM=NELEMD(2,MEDIA)
          IF(NELEM.EQ.1)THEN
              NUCLID=IFIX(ELMIX(1,1,MEDIA))
              GO TO 11
          END IF
C
      DO 1 I=1,NELEM
      NUCLID=IFIX(ELMIX(1,I,MEDIA))
          IRN=IREN(NUCLID)
      TEMP(I)=ELMIX(2,I,MEDIA)*SGLN(NGRUP,7,IRN)
    1 CONTINUE
C
      DO 2 I=2,NELEM
    2 TEMP(I)=TEMP(I-1)+TEMP(I)
C
      BETHA=RRAN(IX)*TEMP(NELEM)
      DO 3 I=1,NELEM
          IF(BETHA.LE.TEMP(I))THEN
              NUCLID=IFIX(ELMIX(1,I,MEDIA))
              GO TO 11
          ELSE
              GO TO 3
          END IF
    3 CONTINUE
C
C
   11 LNUCLD=NUCLID
      IRENLN=IREN(NUCLID)
      RETURN
      END



      SUBROUTINE TYPELN
C SAMPLING OF LOW ENERGY NEUTRON-NUCLEUS INTERACTION TYPE
C INPUT:
C       IRENLN - NO. OF NUCLEUS IN BNAB
C       NGRUP  - NO. OF ENERGY GROUP
C OUTPUT:
C        IPRLN - TYPE OF INTERACTION  (1-4).
C-----------------------------------------------------------
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      COMMON /WNPN/ WN(11),WN0(11),WNMEM(11),WNI(11),WNF(11),NGRUP   ! GEMCA
      COMMON /MANLEN/ LNSTAT,LNMED,LNZON,LNCOD,LNUCLD,IRENLN,IPRLN,
     *                NFISN,NSECN
      COMMON /BNAB/ EGR(-1:27),SPFU(10,-1:13),VNU(10),SGLN(-1:26,7,110),
     *              SUV(13,-1:12,110)                             ! BNAB28
      DIMENSION TEMP(4)
C
C FISSION CROSS SECTION, IPRLN=1
      TEMP(1)=SGLN(NGRUP,1,IRENLN)
C
C CAPTURE CROSS SECTION, IPRLN=2
      TEMP(2)=SGLN(NGRUP,3,IRENLN)
C
C INELASTIC SCATTERING CROSS SECTION, IPRLN=3
      TEMP(3)=SGLN(NGRUP,4,IRENLN)
C
C ELASTIC SCATTERING CROSS SECTION, IPRLN=4
      TEMP(4)=SGLN(NGRUP,5,IRENLN)
C
      DO 10 J=2,4
   10 TEMP(J)=TEMP(J-1)+TEMP(J)
C
   15 GAMMA=RRAN(IX)                         ! Sobol, 15.04.99
        IF(GAMMA.eq.0.0)GO TO 15            ! Sobol, 15.04.99
      BETHA=GAMMA*TEMP(4)                   ! Sobol, 15.04.99
C
      DO 20 J=1,4
          IF(BETHA.LE.TEMP(J))THEN
              IPRLN=J
              GO TO 30
          ELSE
              GO TO 20
          END IF
   20 CONTINUE
C
   30 RETURN
      END



      SUBROUTINE PREPFI(AVERN)
C PREPARATION OF FISSION NEUTRONS SPECTRUM
C INPUT:
C       AVERN - AVERAGE NUMBER OF FISSION NEUTRONS.
C OUTPUT:
C        SPEFIS(13) - INTEGRATED UNNORMED FISSION SPECTRUM
C-------------------------------------------------------------
      COMMON /BNAB/ EGR(-1:27),SPFU(10,-1:13),VNU(10),SGLN(-1:26,7,110),
     *              SUV(13,-1:12,110)                             ! BNAB28
      COMMON /SPECLN/ SPEFIS(-1:13),SPECIN(12)                    ! BNAB28
C
      IF(AVERN.LE.VNU(1))THEN
          J1=1
          J2=2
          GO TO 10
      END IF
C
      IF(AVERN.GE.VNU(10))THEN
          J1=9
          J2=10
          GO TO 10
      END IF
C
      DO 1 J=2,10
          IF(AVERN.LE.VNU(J))THEN
              J1=J-1
              J2=J
              GO TO 10
          ELSE
              GO TO 1
          END IF
    1 CONTINUE
C
   10 DO 20 K=-1,13                                               ! BNAB28
      X1=VNU(J1)
      X2=VNU(J2)
          Y1=SPFU(J1,K)
          Y2=SPFU(J2,K)
      SPEFIS(K)=Y1+((Y2-Y1)/(X2-X1))*(AVERN-X1)
   20 CONTINUE
C
      DO 30 J=0,13                                                ! BNAB28
   30 SPEFIS(J)=SPEFIS(J-1)+SPEFIS(J)
C
      RETURN
      END



      SUBROUTINE PREPIN
C PREPARATION OF SECONDARY NEUTRONS SPECTRUM IN REACTIONS (N,N'),(N,2N)
C INPUT:
C       IRENLN - NO. OF NUCLEUS IN BNAB
C       NGRUP - INSIDENT NEUTRON GROUP NUMBER
C OUTPUT:
C        SPECIN(12) - INTEGRATED UNNORMED INELASTIC SCATTERING SPECTRUM
C----------------------------------------------------------------------
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      COMMON /WNPN/ WN(11),WN0(11),WNMEM(11),WNI(11),WNF(11),NGRUP   ! GEMCA
      COMMON /MANLEN/ LNSTAT,LNMED,LNZON,LNCOD,LNUCLD,IRENLN,IPRLN,
     *                NFISN,NSECN
      COMMON /BNAB/ EGR(-1:27),SPFU(10,-1:13),VNU(10),SGLN(-1:26,7,110),
     *              SUV(13,-1:12,110)                             ! BNAB28
      COMMON /SPECLN/ SPEFIS(-1:13),SPECIN(12)                    ! BNAB28
C
C
      IF(NGRUP.GT.12)THEN
          CONTINUE!WRITE(25,1001)NGRUP,IXINIT,NUMTRE
 1001     FORMAT(' PREPIN DIAGNOSTIC: ILLEGAL NGRUP=',I3,', IXINIT=',
     *           I12,', NUMTRE=',I6)
          STOP
      END IF
C
      DO 1 J=1,12
    1 SPECIN(J)=SUV(J,NGRUP,IRENLN)
C
      DO 2 J=2,12
    2 SPECIN(J)=SPECIN(J-1)+SPECIN(J)
C
      RETURN
      END



      SUBROUTINE UNIFI(SINF,COSF)
C SAMPLING OF UNIFORM AZIMUTAL ANGLE PHI: SIN(PHI),COS(PHI)
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
C
   10 G1=2.0*RRAN(IX)-1.0
      G2=RRAN(IX)
          G1G1=G1*G1
          G2G2=G2*G2
          RR=G1G1+G2G2
              IF(RR.LE.1.0)THEN
                  COSF=(G1G1-G2G2)/RR
                  SINF=(2.0*G1*G2)/RR
                  GO TO 11
              ELSE
                  GO TO 10
              END IF
   11 RETURN
      END



      FUNCTION TNEUT(IPR)
C SAMPLING OF ENERGY OF FISSION NEUTRONS AND INELASTIC SCATTERED NEUTRONS.
C INPUT:
C       NGRUP - INCIDENT NEUTRON GROUP NUMBER
C       IPR=1 - FISSION NEUTRON
C       IPR=3 - INELASTIC SCATTERED NEUTRON
C OUTPUT:
C        TNEUT - SECONDARY NEUTRON ENERGY (MEV)
C-----------------------------------------------------------------------
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      COMMON /WNPN/ WN(11),WN0(11),WNMEM(11),WNI(11),WNF(11),NGRUP   ! GEMCA
      COMMON /BNAB/ EGR(-1:27),SPFU(10,-1:13),VNU(10),SGLN(-1:26,7,110),
     *              SUV(13,-1:12,110)                             ! BNAB28
      COMMON /SPECLN/ SPEFIS(-1:13),SPECIN(12)                    ! BNAB28
C
      IF(IPR.EQ.1)GO TO 10
      IF(IPR.EQ.3)GO TO 20
      CONTINUE!WRITE(25,1001)IPR
 1001 FORMAT(' TNEUT DIAGNOSTIC: ILLEGAL IPR=',I3)
      STOP
C
   10 BETHA=RRAN(IX)*SPEFIS(13)           ! (n,f)
      DO 11 J=-1,13                                               ! BNAB28
          IF(BETHA .LE. SPEFIS(J))THEN
              NGR=J
              TNEUT=EGR(NGR+1)+(EGR(NGR)-EGR(NGR+1))*RRAN(IX)
              GO TO 30
          ELSE
              GO TO 11
          END IF
   11 CONTINUE
C
   20 BETHA=RRAN(IX)*SPECIN(12)           ! (n,n')
      DO 21 J=1,12
        IF(BETHA.LE.SPECIN(J))THEN
          NGR=NGRUP+J-1
              if(J.eq.1)then    ! correction of 10.04.97
                TNEUT=EGR(NGR+1)+(WN(7)   -EGR(NGR+1))*RRAN(IX)
              else
                TNEUT=EGR(NGR+1)+(EGR(NGR)-EGR(NGR+1))*RRAN(IX)
              end if
          GO TO 30
        ELSE
          GO TO 21
        END IF
   21 CONTINUE
C
   30 RETURN
      END



      SUBROUTINE ELSCLN(COST,SINF,COSF,TPRIM)
C ELASTIC SCATTERING OF LOW ENERGY NEUTRON IN LINEAR-ANIZOTROPICAL
C APPROXIMATION ON NUCLEUS NO.LNUCLD
C INPUT:
C       WN(7) - KINETIC ENERGY OF INSIDENT NEUTRON IN LAB. SYSTEM (MEV).
C               NEUTRON FLY ALONG AXIS Z.
C       LNUCLD (=NUCLID) - NUCLEUS-TARGET
C       IRENLN - NO. OF NUCLEUS-TARGET IN BNAB
C OUTPUT:
C        TPRIM - KINETIC ENERGY OF SCATTERED NEUTRON IN LAB. SYSTEM
C        COST:
C        SINF: - DIRECTION OF SCATTERED NEUTRON IN LAB. SYSTEM
C        COSF:
C-----------------------------------------------------------------------
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      COMMON /WNPN/ WN(11),WN0(11),WNMEM(11),WNI(11),WNF(11),NGRUP   ! GEMCA
      COMMON /MANLEN/ LNSTAT,LNMED,LNZON,LNCOD,LNUCLD,IRENLN,IPRLN,
     *                NFISN,NSECN
      COMMON /BNAB/ EGR(-1:27),SPFU(10,-1:13),VNU(10),SGLN(-1:26,7,110),
     *              SUV(13,-1:12,110)                             ! BNAB28
      COMMON /INFEL/ ZNUC(110),ATWEI(110),DENS(110),RIEV(110),SYMB(110)
      CHARACTER*2 SYMB
C
      IF(LNUCLD.EQ.1)GO TO 20
C
C ELASTIC SCATTERING ON NUCLEUS HEAVIER THEN HIDROGEN H(1,1)
      A=ATWEI(LNUCLD)
      T0=WN(7)
      CTL=SGLN(NGRUP,6,IRENLN)
          GAMMA=RRAN(IX)
C-------------------------COS(THETA-LAB)=0.
      IF(CTL.EQ.0.)THEN
          COST=2.0*GAMMA-1.0
          GO TO 10
      ELSE
          GO TO 1
      END IF
C------------------------- -1/3 < COS(THETA-LAB) < 1/3
    1 ABSCTL=ABS(CTL)
          IF(ABSCTL.LT.0.333333333)THEN
              TEMP=0.333333333/CTL
              COST=-TEMP+(ABSCTL/CTL)*SQRT((TEMP-1.)**2+4.*TEMP*GAMMA)
              GO TO 10
          ELSE
              GO TO 2
          END IF
C------------------------ ABS(COS(TETHA-LAB)) > 1/3
    2 TEMP=(1.-CTL)/(1.+CTL)
          COST=2.0*GAMMA**TEMP-1.0
          GO TO 10
C------------------------UNIFORM AZIMUTAL ANGLE
   10 CONTINUE
      IF(COST .lt. -1.0)COST=-1.0   ! Sobol, 26.10.2004
      IF(COST .gt.  1.0)COST= 1.0   ! Sobol, 26.10.2004
      CALL UNIFI(SINF,COSF)
C------------------------SECONDARY KINETIC ENERGY
                if(T0 .LT. 0.215E-06)then                   ! THERMAL
      TPRIM=T0                                              ! THERMAL
                else                                        ! THERMAL
      TPRIM=T0*((COST+SQRT(COST**2+A**2-1.0))/(A+1.0))**2
                end if                                      ! THERMAL
      RETURN
C
C ELASTIC SCATTERING ON HIDROGEN H(1,1)
C Isotropic scattering in c.m.s (2008-05-02): 
C    Cos(chi)=1-2*gamma; chi is isotropic angle of scattering in c.m.s.
C    theta=chi/2 at m1=m2: theta is angle of scattering of progectile in lab.s.
C    [Cos(chi/2)]**2=(1+Cos(chi))/2=1-gamma=>gamma
C    Cos(chi/2)=Cos(theta)=SQRT(gamma)
   20 COST=SQRT(REAL(RRAN(IX)))     ! Cos(theta)=SQRT(gamma)
      CALL UNIFI(SINF,COSF)
                if(WN(7) .LT. 0.215E-06)then                ! THERMAL
      TPRIM=WN(7)                                           ! THERMAL
                else                                        ! THERMAL
      TPRIM=WN(7)*COST**2
                end if                                      ! THERMAL
      RETURN
C
      END



      SUBROUTINE HISTLN(LPOINT)
C REGISTRATION OF LOW ENERGY NEUTRON TRANSPORT CHARACTERISTICS.
C INPUT:
C       LPOINT - NO. OF POINT IN SUBROUTINE LOENT FROM WHICH
C                CALL HISTLN OCCURES.
C                LPOINT=1 - T < TMIN
C                      =2 - END POINT OF NEUTRON RANGE
C                      =3 - NEUTRON FLY OUT FROM TARGET CONFIGURATION
C                      =4 - FISSION
C                      =5 - CAPTURE
C                      =6 - INELASTIC SCATTERING (N,N'),(N,2N)
C                      =7 - ELASTIC SCATTERING
C                      =8 - INPUT NEUIRON SOURCE
C OUTPUT:
C        ESTIM(I,J), J - NO. OF INPUT ZONE
C              I= 1 - NO. OF INPUT ZONE
C              I= 2 - NO. OF MEDIUM IN THIS INPUT ZONE
C              I= 3 - TOTAL FLUX RANGE ESTIMATION
C              I= 4 - TOTAL FLUX COLLISION ESTIMATION
C              I= 5 - NUMBER OF FISSIONS
C              I= 6 - NUMBER OF CAPTURES
C              I= 7 - SUMMAR ENERGY OF CAPTURED NEUTRONS
C              I= 8 - NUMBER OF (N,N')
C              I= 9 - DELTA E AT (N,N')
C              I=10 - NUMBER OF (N,2N)
C              I=11 - DELTA E AT (N,2N)
C              I=12 - NUMBER OF ELASTIC SCATTERING
C              I=13 - DELTA E AT ELASTIC SCATTERING
C        SOURLN - NEUTRON SOURCE SPECTRUM
C        SPLOUT - SPECTRUM OF OUTGOING NEUTRONS
C        RESTRE - SPECTRUM OF NEUTRONS WITH T<TMIN
C-----------------------------------------------------------------------
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      COMMON /WNPN/ WN(11),WN0(11),WNMEM(11),WNI(11),WNF(11),NGRUP   ! GEMCA
      COMMON /MANLEN/ LNSTAT,LNMED,LNZON,LNCOD,LNUCLD,IRENLN,IPRLN,
     *                NFISN,NSECN
      COMMON /STEPLN/ STEPLN(10,500),ITEPLN(6,500),ISTEP,LNN50
      COMMON /REMLN/ REMLN(11,5000),NREMLN,LRL50   ! Koshka_ADS      ! GEMCA
      COMMON /LNHIST/ ESTIM(13,1000),SOURLN(-1:28),SPLOUT(-1:28),  ! 21.02.97 ! KTH
     *                                      RESTRE(-1:28),LES50   ! BNAB28
      REAL*8 ESTIM   ! Loent28_DP
C
C CHECK OF NUMBER OF ZONES
      IF(LNZON.GT.LES50)THEN
          CONTINUE!WRITE(25,1000)LNZON,IXINIT,NUMTRE
 1000     FORMAT(' HISTLN DIAGNOSTIC: LNZON TOO MACH, LNZON=',I4,
     *           ', IXINIT=',I12,', NUMTRE=',I6)
          STOP
      END IF
C
C CHECK OF LPOINT VALIDITY
      IF(LPOINT.LT.1 .OR. LPOINT.GT.8)THEN
          CONTINUE!WRITE(25,1001)LPOINT
 1001     FORMAT(' HISTLN DIAGNOSTIC: ILLEGAL LPOINT=',I3)
          STOP
      END IF
C
C
      GO TO(10,20,30,40,50,60,70,80),LPOINT
C
C SPECTRUM OF NEUTRONS WITH ENERGY T<TMIN (OLN)
   10 CONTINUE
          IGR=NGROUP(WN(7))
      RESTRE(IGR)=RESTRE(IGR)+WN(8)
      RESTRE(27)=RESTRE(27)+WN(8)
      RESTRE(28)=RESTRE(28)+WN(7)*WN(8)
      RETURN
C
C TOTAL FLUX ESTIMATION ON RANGES AND COLLISIONS
   20 CONTINUE
C         ESTIMATION ON RANGES:
      DO 21 J=1,ISTEP
          NZON=ITEPLN(4,J)
          R   =STEPLN(4,J)
          ESTIM(3,NZON)=ESTIM(3,NZON)+R*WN(8)
   21 CONTINUE
C
C         ESTIMATION ON COLLISINS:
      IF(LNSTAT.EQ.0)THEN
          SMTOT=STEPLN(2,ISTEP)
          ESTIM(4,LNZON)=ESTIM(4,LNZON)+WN(8)/SMTOT
          RETURN
      ELSE
          RETURN
      END IF
C
C REGISTRATIONS OF OUTGOING NEUTRONS
   30 CONTINUE
          IGR=NGROUP(WN(7))
      SPLOUT(IGR)=SPLOUT(IGR)+WN(8)
      SPLOUT(27)=SPLOUT(27)+WN(8)
      SPLOUT(28)=SPLOUT(28)+WN(7)*WN(8)
      RETURN
C
C REGISTRATION OF FISSIONS
   40 CONTINUE
      ESTIM(5,LNZON)=ESTIM(5,LNZON)+WN(8)
      RETURN
C
C REGISTRATION OF NEUTRON CAPTURE
   50 CONTINUE
      ESTIM(6,LNZON)=ESTIM(6,LNZON)+WN(8)
      ESTIM(7,LNZON)=ESTIM(7,LNZON)+WN(7)*WN(8)
      RETURN
C
C REGISTRATION OF INELASTICU INTERACTIONS (N,N'),(N,2N)
   60 CONTINUE
      IF(NSECN.LT.1 .OR. NSECN.GT.2)THEN
          CONTINUE!WRITE(25,1002)NSECN,IXINIT,NUMTRE
 1002     FORMAT(' HISTLN DIAGNOSTIC: ILLEGAL NSECN=',I3,', IXINIT=',
     *           I12,', NUMTRE=',I6)
          STOP
      END IF
C
      IF(NSECN.EQ.1)THEN
          ESTIM(8,LNZON)=ESTIM(8,LNZON)+WN(8)
              DE=WN(7)-REMLN(7,NREMLN)
          ESTIM(9,LNZON)=ESTIM(9,LNZON)+DE*WN(8)
          RETURN
      ELSE
          ESTIM(10,LNZON)=ESTIM(10,LNZON)+WN(8)
              DE=WN(7)-REMLN(7,NREMLN-1)-REMLN(7,NREMLN)
          ESTIM(11,LNZON)=ESTIM(11,LNZON)+DE*WN(8)
          RETURN
      END IF
C
C REGISTRATION OF ELASTIC SCATTERING
   70 CONTINUE
      ESTIM(12,LNZON)=ESTIM(12,LNZON)+WN(8)
          DE=WNF(7)-WN(7)
      ESTIM(13,LNZON)=ESTIM(13,LNZON)+DE*WN(8)
      RETURN
C
C REGISTRATION OF INPUT NEUTRON SOURCE
   80 CONTINUE
          IGR=NGROUP(WN(7))
      SOURLN(IGR)=SOURLN(IGR)+WN(8)
      SOURLN(27)=SOURLN(27)+WN(8)
      SOURLN(28)=SOURLN(28)+WN(7)*WN(8)
      RETURN
C
      END



      SUBROUTINE CLERLN
C CLEANING OF ARRAYS FOR LOENT TOTAL OUTPUT
      COMMON /LNHIST/ ESTIM(13,1000),SOURLN(-1:28),SPLOUT(-1:28),  ! 21.02.97 ! KTH
     *                                      RESTRE(-1:28),LES50   ! BNAB28
      REAL*8 ESTIM   ! Loent28_DP
C
      DO 1 J=-1,28                                                ! BNAB28
          SOURLN(J)=0.
          SPLOUT(J)=0.
          RESTRE(J)=0.
    1 CONTINUE
C
      DO 2 J=1,LES50
          DO 2 K=1,13
    2     ESTIM(K,J)=0.
C
      RETURN
      END



      SUBROUTINE PRINLN                                              ! GEMCA
C PRINTING OF LOENT TOTAL OUTPUT
C           COMMON OF Geometry:                                      ! GEMCA
        COMMON/GDATA0/ NUMBOD,NUMZON                                 ! GEMCA
      COMMON /LNHIST/ ESTIM(13,1000),SOURLN(-1:28),SPLOUT(-1:28),  ! 21.02.97 ! KTH
     *                                      RESTRE(-1:28),LES50   ! BNAB28
      REAL*8 ESTIM   ! Loent28_DP
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
          DIMENSION ES(13)
C
      DO 1 J=1,NUMZON                                                ! GEMCA
          ESTIM(1,J)=FLOAT(J)
          ESTIM(2,J)=GZMED(J)                                        ! GEMCA
    1 CONTINUE
C
C PRINTING NEUTRON SPECTRA
      CONTINUE!WRITE(24,10)INT(SOURLN(27))
   10 FORMAT(/10X,'TOTAL OUTPUT OF LOENT TRANSPORT CODE.'/
     *        10X,'STATISTICS IS ',I10,' SOURCE NEUTRONS.'/)
      CONTINUE!WRITE(24,15)
   15 FORMAT(13X,'NEUTRONS SPECTRA')
      CONTINUE!WRITE(24,20)
   20 FORMAT('  NO.',4X,'NEUTRONS',4X,'OUTGOING',4X,'NEUTRONS'/
     *       ' GRUP',6X,'SOURCE',4X,'SPECTRUM WITH T<TMIN')
C
      DO 30 J=-1,26                                               ! BNAB28
   30 CONTINUE!WRITE(24,40)J,SOURLN(J),SPLOUT(J),RESTRE(J)
   40 FORMAT(I5,3F12.3)
      CONTINUE!WRITE(24,50)SOURLN(27),SPLOUT(27),RESTRE(27)
   50 FORMAT('    N',3F12.2)
      CONTINUE!WRITE(24,60)SOURLN(28),SPLOUT(28),RESTRE(28)
   60 FORMAT(' TSUM',3E12.4)
C
C PRINTING DISTRIBUTIONS INSIDE TARGET
      CONTINUE!WRITE(24,100)
  100 FORMAT(/45X,'DISTRIBUTIONS INSIDE TARGET')
      CONTINUE!WRITE(24,110)
  110 FORMAT(' INP. NO.',7X,'RANGE   COLLISION',9X,'NO.',9X,'NO.  ',
     *       'SUM ENERGY',9X,'NO.',5X,'DELTA E',9X,'NO.',5X,
     *       'DELTA E',9X,'NO.',5X,'DELTA E')
      CONTINUE!WRITE(24,120)
  120 FORMAT(' ZONE MED  ESTIMATION  ESTIMATION     FISSION',
     *       '      CAPTUR OF CAPTUR N      (N,N'')',6X,'(N,N'')',
     *       '      (N,2N)',6X,'(N,2N)     EL.SCAT   EL. SCAT.')
C
      DO 131 K=3,13
  131 ES(K)=0.
C
      DO 130 J=1,NUMZON                                              ! GEMCA
          INT1=INT(ESTIM(1,J))   ! Loent28_DP
          INT2=INT(ESTIM(2,J))   ! Loent28_DP
              DO 132 K=3,13
  132         ES(K)=ES(K)+ESTIM(K,J)
      CONTINUE!WRITE(24,140)INT1,INT2,(ESTIM(K,J),K=3,13)
  140 FORMAT(2I4,1X,2E12.5,F12.1,4(F12.1,E12.5))
  130 CONTINUE
      CONTINUE!WRITE(24,141)(ES(K),K=3,13)
  141 FORMAT(' SUMMA: ',1X,2E12.5,F12.1,4(F12.1,E12.5))
C
      LASTIX=IX
      CONTINUE!WRITE(24,150)IXFIRS,LASTIX
  150 FORMAT(/' RANDOM GENERATOR STATUS: IXFIRS=',I12/26X,'LASTIX=',I12)
C
      RETURN
      END



      SUBROUTINE STETLENPA(WGHTLE)                 ! TLEST
C Analysing of the array STEPLN(10,500) and building up
C of the Track Length Estimation for ALL CROSSED ZONES.
C INPUT: WGHTLE - Statistical weight of the particle
C
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      COMMON /STEPLN/ STEPLN(10,500),ITEPLN(6,500),ISTEP,LNN50
      COMMON /TLEPA/ TLEPA(122,20,1000),NG73,LTLEPA
      REAL*8 TLEPA
C
      DO J=1,ISTEP
        NZCUR=ITEPLN(4,J)
        Ei=STEPLN(1,J)
        SLENG=STEPLN(4,J)
        CALL HIST3DPA(DBLE(Ei),TLEPA,1,NZCUR,DBLE(WGHTLE*SLENG))
      END DO
C
   10 CONTINUE
      RETURN
      END
