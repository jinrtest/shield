
***********************************************************************
*                                                                     *
*        The following subroutines are written by N.M.Sobolevsky,     *
*        Institute for Nuclear Research RAS, Moscow,                  *
*        E-mail: sobolevs@AL20.inr.troitsk.ru                         *
*                                                                     *
***********************************************************************

      FUNCTION RNDMD(N)                                 !
C Presentation of the AGT random generator RNDMD(N)     !
C via SHIELD random generator RAN(IX)                   !
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX    !
      REAL*8 RNDMD                                      !
    1 RNDMD=RRAN(IX)                 ! Sobol 02.11.99    !
        IF(RNDMD.eq.0.0)GO TO 1     ! Sobol 02.11.99    !
      RETURN                                            !
      END                                               !
C                                                       !
      FUNCTION IDINTG(D)                                !
      REAL*8 D,D1                                       !
      D1=1.0001*D                                       !
      IDINTG=IDINT(D1)                                  !
      RETURN                                            !
      END                                               !
C * * * * * * * * * * * * * * * * * * * * * * * * * * * *


      SUBROUTINE GENAGT(INOVER)
C AMELIN-GUDIMA-TONEEV HADRON-NUCLEUS GENERATOR. OUTPUT VARIABLE
C INOVER=1 IF ARRAY SPT WAS EXCEEDED, ELSE INOVER=0.
C
C --------------- COMMONS of AGT-generator -----------------------------
      COMMON /RESULT/ AN1,AN2,ZN1,ZN2,ENEXT1,ENEXT2,PNUCL1(3),PNUCL2(3),
     *                AMNUC1(3),AMNUC2(3)
      REAL*8 AN1,AN2,ZN1,ZN2,ENEXT1,ENEXT2,PNUCL1,PNUCL2,AMNUC1,AMNUC2
C
      COMMON /HCASC/ ANUCL1,ANUCL2,ZNUCL1,ZNUCL2,T0,EPS1,EPS2,VPI,A1,A2,
     *               C1,C2,D1,D2,R0N1,R0N2,TF01,TF02,RM1,RM2
      REAL*8 ANUCL1,ANUCL2,ZNUCL1,ZNUCL2,T0,EPS1,EPS2,VPI,A1,A2,
     *       C1,C2,D1,D2,R0N1,R0N2,TF01,TF02,RM1,RM2
C
      COMMON /MEMAGT/ PMEMO(9,5999),IMEMO(5,5999)
      REAL*8 PMEMO
      COMMON /IDPME/ IDPME(5999)     ! added 13.02.96
C
      COMMON /EXCIT/ TEX1,TEX2,HEX1,HEX2,PEX1,PEX2    ! for PRECO !
      REAL*8 TEX1,TEX2,HEX1,HEX2,PEX1,PEX2
C
      COMMON /BK1003/ UPRECO,APRECO,ZPRECO     ! for PRECO; U in MeV !
      REAL*4 UPRECO,APRECO,ZPRECO
C
      COMMON /COMFR/ ICMS            ! Sobolevsky, 28.12.96
C
      COMMON/ITAPES/ITDKY,ITEVT,ITCOM,ITLIS
C Sobolevsky, 06.04.95 ----------
C      LOGICAL ITDKY,ITEVT,ITCOM,ITLIS               ! ITDKY ?, 26.07.95 !
C ---------------------------------
      DIMENSION I1(5),I2(5)
      REAL*8 P1(9),P2(9)
      REAL*8 DELTA,RN
C
C --------------- COMMONS of SHIELD and GENAGT -------------------
      COMMON /NUMINT/ INTIN,INTEL
      COMMON /INREAC/ COST,SINF,COSF,TINT,WINT,JPART,NUCLID,KSTATE
        COMMON /HIPROJ/ APROJ,ZPROJ
      COMMON /SECHAR/ SPT(6,5000),SNU(10,101),LS6,LS100,LS10,LS11  ! GSITR
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      COMMON /BEFORE/ BUP,BAP,BZP,BPNX,BPNY,BPNZ
      COMMON /DEBUG/ LSTAR,LCASC
      COMMON /ISLERR/ ISLERR   ! HADES1
C
      ITLIS=16                                ! Don`t used !
C Input in c.m. system (it is need after H calculation)
      ICMS=1                       ! Dementyev 24.06.94
C
***************************************
      A1=0.              ! Barz, Feb.98
      A2=0.              ! Barz, Feb.98
      TF01=0.            ! Barz, Feb.98
      TF02=0.            ! Barz, Feb.98
      RM1=0.             ! Barz, Feb.98
      RM2=0.             ! Barz, Feb.98
***************************************
C
C --------------- Cleaning of INOVER and arrays SPT, SNU
        INOVER=0
      DO 1 J=1,LS100
      DO 1 I=1,LS6
    1 SPT(I,J)=0.
        DO 2 J=1,LS11
        DO 2 I=1,LS10
    2     SNU(I,J)=0.
C
C --------------- Cleaning of arrays PMEMO, IMEMO (Sobolevsky, 20.07.95)
      DO J=1,5999
        DO K=1,9
            PMEMO(K,J)=0.
        END DO
        DO L=1,5
            IMEMO(L,J)=0
        END DO
            IDPME(J)=0             ! added 13.02.96
      END DO
C
C Conversion of SHIELD variables into input GENAGT variables
      CALL CONAGT
C
C BRANCHING: hA- or hp-interaction?
      IF(ANUCL2.LT.1.9 .and. anucl1.lt.1.9)GO TO 3    ! ===> hp !
        DELTA=1.3
        RN=0.2
C
C   INTRANUCLEAR PARAMETERS DISTRIBUTION
      IF(ANUCL1.GT.1.1)CALL GHELP(R0N1,ANUCL1,A1,C1,D1,TF01,RM1)
      IF(ANUCL2.GT.1.1)CALL GHELP(R0N2,ANUCL2,A2,C2,D2,TF02,RM2)
cC =========Debug printout of HELP results ============
c      CONTINUE!WRITE(25,214)
c  214 FORMAT('     R0N   ANUCL       A       C       D     TF0      RM')
c      CONTINUE!WRITE(25,215)R0N1,ANUCL1,A1,C1,D1,TF01,RM1
c      CONTINUE!WRITE(25,215)R0N2,ANUCL2,A2,C2,D2,TF02,RM2
c  215 FORMAT(7F8.5)
cC ================end debug=======================
C
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C   Inelastic AGT-interaction                               !
      CALL CASCAN(NEL,RN,DELTA,MV,IRET,INOVER) ! KarSjh     !!   13.02.96
      IF(INOVER.eq.1)RETURN ! KarSjh: Abnormal completion of CASCAN 
C      CALL IMPACT(NEL)   ! HADES
        if(ISLERR.eq.1)return    ! HADES1 
c        if(ISLERR.eq.1)then
c          CONTINUE!WRITE(25,*)'       CASCAN'
c          return
c        end if
      INTEL=INTEL+NEL     ! Number of elastic scattering    !
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C   CHECK OF ARRAY SPT OVERFILLING or Cascan problems (IRET.ne.0)
      IF(MV.GE.LS100-1)THEN
         INOVER=1
         CONTINUE!WRITE(25,1000)MV,IXINIT,NUMTRE
 1000      FORMAT(' GENAGT DIAGNOSTIC: OVERFILLING OF SPT',
     *     'MV=',I4,'  IXINIT=',I12,'  NUMTRE=',I6)
         MV=LS100-1
      END IF
C Check of IRET           ! 13.02.96
      IF(IRET.NE.0)THEN
         INOVER=1
         CONTINUE!WRITE(25,1001)IRET,IXINIT,NUMTRE
 1001      FORMAT(' GENAGT DIAGNOSTIC: CASCAN PROBLEM',
     *     'IRET=',I4,'  IXINIT=',I12,'  NUMTRE=',I6)
      END IF
C
C ============ Debug printout of /MEMAGT/ after AGT ========
      IF(LSTAR.eq.2)THEN
      CONTINUE!WRITE(36,300)NUMTRE
  300   FORMAT(/20X,'************** STAR No.',I4,' **************')
      CALL PRIAGT(MV)
      END IF
C ================== End of debug printout==================
C
cC * * * * * * * * * * * * * * * * * * * * * * * *
cC Coalescence in arrays of /MEMAGT/             !
c      CALL COALES(MV)                           !
cC * * * * * * * * * * * * * * * * * * * * * * * *
C
cC ============ Debug printout of /MEMAGT/ after COALES ========
c      IF(LSTAR.eq.2)CALL PRICOA(MV)
cC ================== End of debug printout==================
C
C  Conversion of AGT-generator secondary particles in /MEMAGT/ into
C  SHIELD secondary particles in /SECHAR/.
      CALL AGTCON(MV,KSTART)
C
C ---------------------- (Sobol 25.07.95) ------------------
C Pre-compaund emission and equilibrium deexcitation
C                  of TARGET-nucleus.
C Debug printout for nucleus-target, 23.12.96 ===============
      IF(LSTAR.eq.2)CONTINUE!WRITE(36,310)
  310 FORMAT(/20X,'NUCLEUS-TARGET:')
C ============== End of debug printout =======================
C
      NSTART=2      ! for the first run of DEEX !
C
C PRE-COMPAUND:
        KSMEM1=MV+1       ! first free line in MEMAGT for PRECO !
        KSMEMO=KSMEM1     ! current free line number !
      ENEXT= SNGL(ENEXT2)
        PNXPRE=SNGL(PNUCL2(1))
        PNYPRE=SNGL(PNUCL2(2))
        PNZPRE=SNGL(PNUCL2(3))
            ATWGHT=SNGL(AN2)
            CHARGE=SNGL(ZN2)
              AMPRE=0.1
              AMFPRE=0.1
              R0PRE=1.3
                  N0=INT(TEX2)
                  P0=SNGL(TEX2-HEX2)
                  H0=SNGL(HEX2)
                  PZ0=SNGL(PEX2)
      CALL PRECO(ENEXT,ATWGHT,CHARGE,PNXPRE,PNYPRE,PNZPRE,
     *           AMPRE,AMFPRE,R0PRE,KSMEMO,N0,P0,H0,PZ0,K1)
      CALL PRESPT(KSMEM1,KSMEMO,KSTART,INOVER,PNXPRE,PNYPRE,PNZPRE)
C
C DEEXITATION:
C  Transformation GeV to MeV before deexitation
      UP= UPRECO
      PNX=PNXPRE*1000.
      PNY=PNYPRE*1000.
      PNZ=PNZPRE*1000.
        AP = APRECO
        ZP = ZPRECO
C Nucleus before deexcitation
            BUP= UP
            BAP= AP
            BZP= ZP
            BPNX=PNX
            BPNY=PNY
            BPNZ=PNZ
      CALL DEEX(UP,AP,ZP,PNX,PNY,PNZ,KSTART,NSTART,INOVER)
C ===========Debug printout of SPTU after DEEX =======
      IF(LSTAR.eq.2)CALL PRSPTU
C ================== End of debug printout====================
C
C Pre-compaund emission and equilibrium deexcitation
C                      of PROJECTILE-nucleus, if any.
      IF(JPART.le.11)GO TO 140         ! projectile is particle      !HION2
C
C Debug printout for nucleus-projectile, 23.12.96 ===============
      IF(LSTAR.eq.2)CONTINUE!WRITE(36,320)
  320 FORMAT(/20X,'NUCLEUS-PROJECTILE:')
C ============== End of debug printout =======================
C Projectile is nucleus
      IF(NSTART.EQ.2)NSTART=3         ! for second run of DEEX

C PRE-COMPAUND:
      KSMEM1=KSMEMO        ! first free line in MEMAGT for PRECO !
      ENEXT= SNGL(ENEXT1)
        PNXPRE=SNGL(PNUCL1(1))
        PNYPRE=SNGL(PNUCL1(2))
        PNZPRE=SNGL(PNUCL1(3))
            ATWGHT=SNGL(AN1)
            CHARGE=SNGL(ZN1)
              AMPRE=0.1
              AMFPRE=0.1
              R0PRE=1.3
                  N0=INT(TEX1)
                  P0=DBLE(TEX1-HEX1)
                  H0=DBLE(HEX1)
                  PZ0=DBLE(PEX1)
      CALL PRECO(ENEXT,ATWGHT,CHARGE,PNXPRE,PNYPRE,PNZPRE,
     *           AMPRE,AMFPRE,R0PRE,KSMEMO,N0,P0,H0,PZ0,K1)
      CALL PRESPT(KSMEM1,KSMEMO,KSTART,INOVER,PNXPRE,PNYPRE,PNZPRE)
C
C DEEXITATION:
C  Transformation GeV to MeV before deexitation
      UP= UPRECO
      PNX=PNXPRE*1000.
      PNY=PNYPRE*1000.
      PNZ=PNZPRE*1000.
        AP = APRECO
        ZP = ZPRECO
C Nucleus before deexcitation
            BUP= UP
            BAP= AP
            BZP= ZP
            BPNX=PNX
            BPNY=PNY
            BPNZ=PNZ
      CALL DEEX(UP,AP,ZP,PNX,PNY,PNZ,KSTART,NSTART,INOVER)
C ===========Debug printout of SPTU after DEEX =======
      IF(LSTAR.eq.2)CALL PRSPTU
C ================== End of debug printout====================
C
C  FIRST LINE OF ARRAY SNU
  140 CONTINUE
      SNU(1,1)=SNGL(ANUCL2)
      SNU(2,1)=SNGL(ZNUCL2)
      SNU(4,1)=2.
C
C Rotation of residual nuclei directions in array SNU
      DO 130 J=2,LS11
        IF(SNU(1,J).EQ.0.)GO TO 131
        PNXF=SNU(5,J)
        PNYF=SNU(6,J)
        PNZF=SNU(7,J)
        CALL DECARS(PNXF,PNYF,PNZF,PNF,COSTN,SINFN,COSFN)
        CALL SUBROT(COST,SINF,COSF,COSTN,SINFN,COSFN,CTR,SFR,CFR,STR)
        SNU(5,J)=PNF*STR*CFR
        SNU(6,J)=PNF*STR*SFR
        SNU(7,J)=PNF*CTR
  130 CONTINUE
  131 CONTINUE
C
C  ROTATION OF THE SECONDARY PARTICLES DIRECTIONS IN ARRAY SPT
      DO 110 J=1,LS100
        IF(SPT(6,J) .EQ. 0.0) GOTO 111
        CALL SUBROT(COST,SINF,COSF,SPT(1,J),SPT(2,J),SPT(3,J),
     *                CTR,SFR,CFR,STR)
        SPT(1,J)=CTR
        SPT(2,J)=SFR
        SPT(3,J)=CFR
  110 CONTINUE
  111 CONTINUE
C
C  PARTICLES TYPE AND WEIGHT APPROPRIATION TO SECONDARY PARTICLES IN SPT
      DO 120 J=1,LS100
        IF(SPT(6,J).EQ.0.0)GOTO 122
        CALL TYPART(IFIX(SPT(5,J)),SPT(6,J),JPR)
        SPT(5,J)=FLOAT(JPR)
        SPT(6,J)=1.
  120 CONTINUE
  122 CONTINUE
C Debug printout of complete STAR (13.02.96)
       IF(LSTAR.eq.1 .or. LSTAR.eq.2)call star(NUMTRE)
C ==================End of debug printout======
      RETURN
C
C
C **************** HADRON - NUCLEON INTERACTION ****************
  3   CALL CNAGTH(P1,I1,P2,I2)
C================= KTH ===================
      IELACO=0
  4   CALL AGTH1(P1,I1,P2,I2)
      IF(IELA.EQ.1)THEN
         IELACO=IELACO+1
         IF(IELACO.le.1000)THEN
           GO TO 4
         ELSE
           INOVER=1
           RETURN
         END IF
      END IF
C================ end KTH ================
      CALL AGTCNH(KSTART,IELA)
C
C ............... Added by Sobolevsky, 28.12.96 ...............
C  ROTATION OF THE SECONDARY PARTICLES DIRECTIONS IN ARRAY SPT
      DO 410 J=1,LS100
        IF(SPT(6,J) .EQ. 0.0)GO TO 411
        CALL SUBROT(COST,SINF,COSF,SPT(1,J),SPT(2,J),SPT(3,J),
     *                CTR,SFR,CFR,STR)
        SPT(1,J)=CTR
        SPT(2,J)=SFR
        SPT(3,J)=CFR
  410 CONTINUE
  411 CONTINUE
C ...................................
C
C  PARTICLES TYPE AND WEIGHT APPROPRIATION TO SECONDARY PARTICLES IN SPT
      DO 420 J=1,LS100
        IF(SPT(6,J) .EQ. 0.0)GO TO 421
        CALL TYPART(IFIX(SPT(5,J)),SPT(6,J),JPR)
        SPT(5,J)=FLOAT(JPR)
        SPT(6,J)=1.
  420 CONTINUE
  421 CONTINUE
C
C Dementyev 24.06.94 ---------------------
C Appropriation of characteristics of nucleus-target (hidrogen)
C and residual nucleus (full desintegration).
C     First line of array SNU
      SNU(1,1)=SNGL(ANUCL2)
      SNU(2,1)=SNGL(ZNUCL2)
      SNU(4,1)=2.
C     Second line of array SNU
      SNU(4,2)=6.
C ----------------------------------------
C ====Debug printout of Amelin's PPTCL and complete STAR (21.02.96)====
      IF(LSTAR.eq.2)CONTINUE!WRITE(36,300)NUMTRE
      IF(LSTAR.eq.2)call prihh
      IF(LSTAR.eq.1 .or. LSTAR.eq.2)call star(NUMTRE)
C ================ end debug =============
C
      RETURN
      END


      SUBROUTINE SIGEOM(SGEOM,APR,ZPR,ATG,ZTG)
C Calculation of geometric cross section of hA- and AA-interactions
C
C Input: Data in COMMONS of AGT and SHIELD
C Output: SGEOM
C         APR,ZPR - Projectile
C         ATG,ZTG - Target
C
C --------------- COMMONS of AGT-generator
      COMMON /HCASC/ ANUCL1,ANUCL2,ZNUCL1,ZNUCL2,T0,EPS1,EPS2,VPI,A1,A2,
     *               C1,C2,D1,D2,R0N1,R0N2,TF01,TF02,RM1,RM2
      REAL*8 ANUCL1,ANUCL2,ZNUCL1,ZNUCL2,T0,EPS1,EPS2,VPI,A1,A2,
     *       C1,C2,D1,D2,R0N1,R0N2,TF01,TF02,RM1,RM2
      COMMON /STIN/ STIN,AMIN
      REAL*8 STIN,AMIN
      COMMON /XBMAX/ XBMAX,IFIB0
      REAL*8 XBMAX
        REAL*8 DELTA,AM,DLAM
        REAL*8 ATWEI1,ATWEI2
C
C --------------- COMMON of SHIELD
      COMMON /INFEL/ ZNUC(110),ATWEI(110),DENS(110),RIEV(110),SYMB(110)
      CHARACTER*2 SYMB
C
C Conversion of SHIELD variables into input GENAGT variables
      CALL CONAGT
        DELTA=1.3
C
C Atomic weights of natural isitope mixture
      IZN1=INT(ZNUCL1)
      IZN2=INT(ZNUCL2)
        ATWEI1=DBLE(ATWEI(IZN1))
        ATWEI2=DBLE(ATWEI(IZN2))
C
C Calculation of intranuclear parameters
      IF(ANUCL1.GT.1.1)CALL GHELP(R0N1,atwei1,A1,C1,D1,TF01,RM1)
      IF(ANUCL2.GT.1.1)CALL GHELP(R0N2,atwei2,A2,C2,D2,TF02,RM2)
C
C Calculation of geometric cross section SGEOM
        AM=0.940
        IF(ANUCL1.LE.1.1)AM=AMIN
        DLAM=1./(5.06*DSQRT(T0*(T0+2.*AM)))
      SGEOM=sngl(31.459*((RM1+RM2+DLAM+DELTA)**2)*XBMAX**2)
C
C Output (Ap,Zp) and (At,Zt)
      APR=SNGL(ATWEI1)
      ZPR=SNGL(ZNUCL1)
        ATG=SNGL(ATWEI2)
        ZTG=SNGL(ZNUCL2)
C
      RETURN
      END



      SUBROUTINE CONAGT
C  Conversion of extranuclear variables JPART and NUCLID into input
C  variables and parameters for AGT-generator
C
C --------------- COMMONS of AGT-generator -----------------------------
      COMMON /HCASC/ ANUCL1,ANUCL2,ZNUCL1,ZNUCL2,T0,EPS1,EPS2,VPI,A1,A2,
     *               C1,C2,D1,D2,R0N1,R0N2,TF01,TF02,RM1,RM2
      REAL*8 ANUCL1,ANUCL2,ZNUCL1,ZNUCL2,T0,EPS1,EPS2,VPI,A1,A2,
     *       C1,C2,D1,D2,R0N1,R0N2,TF01,TF02,RM1,RM2
      COMMON /STIN/ STIN,AMIN
      REAL*8 STIN,AMIN
C
      COMMON /XBMAX/  XBMAX,IFIB0
      COMMON /PIDABS/ PIDABS
      COMMON /TLIMIT/ TLIMIT
      COMMON /MMATUR/ MMES,MBAR
      COMMON /IACT/   IACT
      COMMON /INDDEC/ INDDEC
      COMMON /INDRES/ INDRES     ! INDRES=1, hadron is decayed !
      COMMON /ISOB2/  ISOB2
      COMMON /ISOB3/  ISOB3
      COMMON /INDINT/ INDINT
      COMMON /INTCC/  INTCC
      COMMON /KSYST/  KSYST
      REAL*8 XBMAX,PIDABS,TLIMIT
      REAL*8 MMES,MBAR
C
C --------------- COMMONS of SHIELD -----------------------------------
      COMMON /INFEL/ ZNUC(110),ATWEI(110),DENS(110),RIEV(110),SYMB(110)
      CHARACTER*2 SYMB
      COMMON /INREAC/ COST,SINF,COSF,TINT,WINT,JPART,NUCLID,KSTATE
        COMMON /HIPROJ/ APROJ,ZPROJ
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
        COMMON /LEVERS/ PARLEV(40)    ! Parameters, see Subr. INLEVR
      COMMON /ANTIL/ lantil                                         ! KarSjh
C
C --------------- Nucleus-target ------------------
C      ANUCL2=DBLE(NINT(ATWEI(NUCLID)))
      ANUCL2=DBLE(AMIXIS(NUCLID))
      ZNUCL2=DBLE(ZNUC(NUCLID))
C========================== KarSjh ===================================
C To give Hev*Ion-Projectite (APROJ,ZPROJ) instead of target-nucleus
C in the case of inelastic interaction of Hev*Ion with Hidrogen-target.
C See SUBROUTINE REACT. This don't affect the separate working of HADGEN.
      IF(lantil.EQ.1)THEN
        ANUCL2=DBLE(APROJ)
      END IF
C========================= end KarSjh ================================
C                                Default
      EPS2=DBLE(PARLEV( 5))    ! 0.007 average (for product far from target)
      VPI= DBLE(PARLEV( 6))    ! 0.025
      C2=  DBLE(PARLEV( 7))  ! 0.545 for A>10; for A<10 don't used
      D2=  DBLE(PARLEV( 8))  ! 0.05
      R0N2=DBLE(PARLEV( 9))  ! 1.07 for A>10
      RM2=0.
C
C Parameter R0N for light nuclei (A<=10)
        LA=INT(ANUCL2)
        LZ=INT(ZNUCL2)                              ! Default
      IF(LA.EQ.2 .AND. LZ.EQ.1)R0N2=DBLE(PARLEV(10))  ! 2.31
      IF(LA.EQ.3 .AND. LZ.EQ.1)R0N2=DBLE(PARLEV(11))  ! 2.4
      IF(LA.EQ.3 .AND. LZ.EQ.2)R0N2=DBLE(PARLEV(12))  ! 2.4
      IF(LA.EQ.4 .AND. LZ.EQ.2)R0N2=DBLE(PARLEV(13))  ! 1.60
      IF(LA.EQ.6 .AND. LZ.EQ.3)R0N2=DBLE(PARLEV(14))  ! 2.4
      IF(LA.EQ.7 .AND. LZ.EQ.3)R0N2=DBLE(PARLEV(15))  ! 2.4
      IF(LA.EQ.9 .AND. LZ.EQ.4)R0N2=DBLE(PARLEV(16))  ! 2.48
      IF(LA.EQ.10.AND. LZ.EQ.5)R0N2=DBLE(PARLEV(17))  ! 2.48
C
C General AGT parameters -------------------
      XBMAX=1.         ! =R/Rmax
      IFIB0=1          ! =0 to fix reaction plane?
C                                Default
      PIDABS=DBLE(PARLEV( 1))  ! 0.01 absorbtion via pi-D
      TLIMIT=DBLE(PARLEV( 2))  ! 25. total interaction time in fm/c
      MMES=  DBLE(PARLEV( 3))  ! 0.2 for T0<1 GeV one may use MMES=999.9(?)
      MBAR=  DBLE(PARLEV( 4))  ! 999.999
        IACT=2
        INDDEC=0
        INDRES=1    ! Duplicated in CASCAN !
        ISOB2=1
        ISOB3=1
        INDINT=1
        INTCC=1   ! 0 - NO CCint; 1 - NN&piN included; 2 - all incl.
        KSYST=1   ! 1- Lab. syst.; 2 - eq.vel.syst.; 3 - c.m.s.
C
C --------- Is the projectile a particle or a heavy ion ? -----------!HION2
      T0=DBLE(0.001*TINT)   ! T0 - kin. energy per nucleon
          IF(JPART.ge.21)GO TO 20   ! Heavy Ion: JPART=21-25         !HION2
          IF(JPART.GT.11 .AND. JPART.le.20)GO TO 1000                !HION2
C
C Incident particle ----------------------
      GO TO(1,2,3,4,5,6,7,8,9,10,11),JPART
    1 ANUCL1=1.
      ZNUCL1=0.
        AMIN=0.940
            STIN=0.
              GO TO 12
    2 ANUCL1=1.
      ZNUCL1=1.
        AMIN=0.940
            STIN=0.
              GO TO 12
    3 ANUCL1=0.
      ZNUCL1=-1.
        AMIN=0.140
            STIN=0.
              GO TO 12
    4 ANUCL1=0.
      ZNUCL1=1.
        AMIN=0.140
            STIN=0.
              GO TO 12
    5 GO TO 1000
    6 ANUCL1=-1.
      ZNUCL1=0.
        AMIN=0.94
            STIN=0.
              GO TO 12
    7 ANUCL1=-1.
      ZNUCL1=-1.
        AMIN=0.94
            STIN=0.
              GO TO 12
    8 ANUCL1=0.
      ZNUCL1=-1.
        AMIN=0.495
            STIN=-1.
              GO TO 12
    9 ANUCL1=0.
      ZNUCL1=1.
        AMIN=0.495
            STIN=1.
              GO TO 12
   10 ANUCL1=0
      ZNUCL1=0.
        AMIN=0.495
            STIN=1.
              GO TO 12
   11 ANUCL1=0.
      ZNUCL1=0.
        AMIN=0.495
            STIN=-1.
              GO TO 12
   12 RETURN
C
 1000 CONTINUE!WRITE(25,1001) JPART,IXINIT,NUMTRE
 1001 FORMAT('CONAGT DIAGNOSTIC: INVALID JPART= ',I2,'  IXINIT=',I12,
     * '  NUMTRE=',I4)
      STOP
C
C Incident heavy ion -----------------------
   20 CONTINUE
      ANUCL1=DBLE(APROJ)
      ZNUCL1=DBLE(ZPROJ)
C                                Default
      EPS1=DBLE(PARLEV( 5))    ! 0.007 average (for product far from target)
      VPI= DBLE(PARLEV( 6))    ! 0.025
      C1=  DBLE(PARLEV( 7))  ! 0.545 for A>10; for A<10 don't used
      D1=  DBLE(PARLEV( 8))  ! 0.05
      R0N1=DBLE(PARLEV( 9))  ! 1.07 for A>10
      RM1=0.
C
C Parameter R0N for light nuclei (A<=10)
        LA=INT(ANUCL1)
        LZ=INT(ZNUCL1)                              ! Default
      IF(LA.EQ.2 .AND. LZ.EQ.1)R0N1=DBLE(PARLEV(10))  ! 2.31
      IF(LA.EQ.3 .AND. LZ.EQ.1)R0N1=DBLE(PARLEV(11))  ! 2.4
      IF(LA.EQ.3 .AND. LZ.EQ.2)R0N1=DBLE(PARLEV(12))  ! 2.4
      IF(LA.EQ.4 .AND. LZ.EQ.2)R0N1=DBLE(PARLEV(13))  ! 1.60
      IF(LA.EQ.6 .AND. LZ.EQ.3)R0N1=DBLE(PARLEV(14))  ! 2.4
      IF(LA.EQ.7 .AND. LZ.EQ.3)R0N1=DBLE(PARLEV(15))  ! 2.4
      IF(LA.EQ.9 .AND. LZ.EQ.4)R0N1=DBLE(PARLEV(16))  ! 2.48
      IF(LA.EQ.10.AND. LZ.EQ.5)R0N1=DBLE(PARLEV(17))  ! 2.48
C
      RETURN
      END



      SUBROUTINE AGTCON(MV,KSTART)
C Conversion of AGT-generator secondary particles in /MEMAGT/
C into SHIELD secondary particles in /SECHAR/.
C     INPUT: MV - number of particles in /MEMAGT/
C     OUTPUT: KSTART - the first free line in SPT
C
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      COMMON /SECHAR/ SPT(6,5000),SNU(10,101),LS6,LS100,LS10,LS11  ! GSITR
      COMMON /MEMAGT/ PMEMO(9,5999),IMEMO(5,5999)
      REAL*8 PMEMO
        COMMON /IDPME/ IDPME(5999)     ! added 30.07.95
C
C COMMONs for AGT==>SHIELD conversion protocol, 30.12.96
      COMMON /CONVEN/ ECONV(65),NCONV(65),JCONV(65),UEXCIT,ENKIN
      COMMON /RESULT/ AN1,AN2,ZN1,ZN2,ENEXT1,ENEXT2,PNUCL1(3),PNUCL2(3),
     *                AMNUC1(3),AMNUC2(3)
      REAL*8 AN1,AN2,ZN1,ZN2,ENEXT1,ENEXT2,PNUCL1,PNUCL2,AMNUC1,AMNUC2
C
      KSTART=MV+1
      JJ=0
C
      DO 100 J=1,MV
      JJ=JJ+1
        IBAR=IMEMO(4,J)
        ISTR=IMEMO(3,J)
        ICHR=IMEMO(1,J)
            IDP=IDPME(J)               ! added 30.07.95
        TKMeV=SNGL(PMEMO(8,J))*1000.
C      
C Is the particle No.J hadron or coalesced lightest nucleus?
      IF(IDP.ne.0)THEN 
          CALL IDPJPR(IDP,IPART,TKMeV)     ! hadron
      ELSE
          IPART=0       ! Coalescence (IDP=0): d, t, He3 or alpha
          IF(IBAR.eq.2 .and. ICHR.eq.1)IPART=21
          IF(IBAR.eq.3 .and. ICHR.eq.1)IPART=22
          IF(IBAR.eq.3 .and. ICHR.eq.2)IPART=23
          IF(IBAR.eq.4 .and. ICHR.eq.2)IPART=24
              IF(IPART.eq.0)CONTINUE!WRITE(25,1001)NUMTRE
 1001         FORMAT(' AGTCON-COALES PROBLEMS, NUMTRE=',I8)
      END IF
C
C Is the particle nucleon, pion, antinucleon, kaon, antikaon OR GAMMA?
        if(IPART .ne. 0) go to 20
        JJ=JJ-1
        KSTART=KSTART-1
C
      GO TO 100
C                  To the end of the loop (this is unknown particle)
C
C  The particle is a nucleon, pion, antinucleon, kaon, antikaon OR GAMMA.
   20 CONTINUE
        PX=SNGL(PMEMO(4,J))
        PY=SNGL(PMEMO(5,J))
        PZ=SNGL(PMEMO(6,J))
      CALL DECARS(PX,PY,PZ,PMOD,COSTP,SINFP,COSFP)
        SPT(1,JJ)=COSTP
        SPT(2,JJ)=SINFP
        SPT(3,JJ)=COSFP
        SPT(4,JJ)=SNGL(PMEMO(8,J))*1000.
C -----------------------------------------------------------------
C Appropriation of charges and masses to particles according to IPART. 
C AntiN and antiK0 have got charges 6 and 11 respectively to differ them
C from N and K0. The inverse thansformation of charge and mass to JPART 
C is made by subroutine TYPART.
      IF(IPART.eq.1)THEN      ! Neutron
        SPT(5,JJ)=0.
        SPT(6,JJ)=940.
        GO TO 100
      END IF
      IF(IPART.eq.2)THEN      ! Proton
        SPT(5,JJ)=1.
        SPT(6,JJ)=940.
        GO TO 100
      END IF
      IF(IPART.eq.3)THEN      ! PI -
        SPT(5,JJ)=-1.
        SPT(6,JJ)=140.
        GO TO 100
      END IF
      IF(IPART.eq.4)THEN      ! PI +
        SPT(5,JJ)=1. 
        SPT(6,JJ)=140.
        GO TO 100
      END IF
      IF(IPART.eq.5)THEN      ! PI 0 ?????? TODO
        SPT(5,JJ)=0.
        SPT(6,JJ)=140.
        GO TO 100
      END IF
      IF(IPART.eq.6)THEN      ! Antineutron  
          SPT(5,JJ)=6.           ! Let it is the charge of ANTI N
        SPT(6,JJ)=940.
        GO TO 100
      END IF
      IF(IPART.eq.7)THEN      ! Antiproton
        SPT(5,JJ)=-1.
        SPT(6,JJ)=940.
        GO TO 100
      END IF
      IF(IPART.eq.8)THEN      ! K-
        SPT(5,JJ)=-1.
        SPT(6,JJ)=495.
        GO TO 100
      END IF
      IF(IPART.eq.9)THEN      ! K+
        SPT(5,JJ)=1.
        SPT(6,JJ)=495.
        GO TO 100
      END IF
      IF(IPART.eq.10)THEN     ! K0
        SPT(5,JJ)=0.
        SPT(6,JJ)=495.
        GO TO 100
      END IF
      IF(IPART.eq.11)THEN     ! Anti K0
          SPT(5,JJ)=11.          ! Let it is the charge of ANTI K0
        SPT(6,JJ)=495.
        GO TO 100
      END IF
      IF(IPART.eq.12)THEN     ! Gamma
          SPT(5,JJ)=12.          ! Let it is the charge of GAMMA QUANT !
        SPT(6,JJ)=12.            ! Let it is the mass of GAMMA QUANT !
        GO TO 100
      END IF
      IF(IPART.eq.21)THEN     ! Deuterium
        SPT(5,JJ)=1.
        SPT(6,JJ)=1880.
        GO TO 100
      END IF
      IF(IPART.eq.22)THEN     ! Tritium
        SPT(5,JJ)=1.
        SPT(6,JJ)=2820.
        GO TO 100
      END IF
      IF(IPART.eq.23)THEN     ! He3
        SPT(5,JJ)=2.
        SPT(6,JJ)=2820.
        GO TO 100
      END IF
      IF(IPART.eq.24)THEN     ! alpha
        SPT(5,JJ)=2.
        SPT(6,JJ)=3760.
        GO TO 100
      END IF
C ------------------ end of IPART appropriation -----------
 100  CONTINUE
C
C Excitation energy of residual nuclei after AGT (for conversion protocol)
      UEXCIT=UEXCIT+SNGL(ENEXT1+ENEXT2)*1000.
C
C Kinetic energy of residual nuclei after AGT (for conversion protocol)
      VM1=SNGL(AN1)*0.94  ! Mass
      VM2=SNGL(AN2)*0.94
        P12=SNGL(PNUCL1(1)**2+PNUCL1(2)**2+PNUCL1(3)**2)  ! P**2
        P22=SNGL(PNUCL2(1)**2+PNUCL2(2)**2+PNUCL2(3)**2)
            EKIN1=SQRT(P12+VM1**2)-VM1
            EKIN2=SQRT(P22+VM2**2)-VM2
      ENKIN=ENKIN+(EKIN1+EKIN2)*1000.
C
      RETURN
      END



      SUBROUTINE IDPJPR(IDP,JPART,TKMeV)
C Conversion of the AGT's numbering scheme to SHIELD's one:
C                 Input IDP --> Output JPART
C                 TKMeV - particle energy for diagnostic printout
      CHARACTER*8  PNAME
      COMMON /PNAME/ PNAME(65),IDENTP(65)
      COMMON /CONVEN/ ECONV(65),NCONV(65),JCONV(65),UEXCIT,ENKIN
C
      DIMENSION PARTIC(24)
      CHARACTER*7 PARTIC
      DATA PARTIC /'NEUTRON',' PROTON','  PI-  ','  PI+  ','  PI0  ',
     *             ' ANTI N',' ANTI P','  K-   ','  K+   ','  K0   ',
     *             'ANTI K0',' GAMMA ','ELECTRN','POSITRN','  MU-  ',
     *             '  MU+  ','  NU E ','ANTINUE',' NU MU ','ANTINUM',
     *             'DEUTRON','TRITIUM','HELIUM3',' ALPHA '/
C
C The AGT's numbering scheme:
C   #   PARTICLE   IDP                #   PARTICLE   IDP
C   1   PI+        120               34   AXI*-    -2331
C   2   PI-       -120               35   AXI*0    -1331
C   3   K+         130               36   AOM-     -3331
C   4   K-        -130               37   P         1120
C   5   K0         230               38   N         1220
C   6   AK0       -230               39   S+        1130
C   7   PI0        110               40   S-        2230
C   8   ETA        220               41   S0        1230
C   9   ETAP       330               42   XI-       2330
C  10   RHO+       121               43   XI0       1330
C  11   RHO-      -121               44   L         2130
C  12   K*+        131               45   DL++      1111
C  13   K*-       -131               46   DL+       1121
C  14   K*0        231               47   DL-       2221
C  15   AK*0      -231               48   DL0       1221
C  16   RHO0       111               49   S*+       1131
C  17   OMEG       221               50   S*-       2231
C  18   PHI        331               51   S*0       1231
C  19   AP       -1120               52   XI*-      2331
C  20   AN       -1220               53   XI*0      1331
C  21   AS+      -1130               54   OM-       3331
C  22   AS-      -2230               55   KS          20
C  23   AS0      -1230               56   KL         -20
C  24   AXI-     -2330               57   GM          10
C  25   AXI0     -1330               58   E+          12
C  26   AL       -2130               59   E-         -12
C  27   ADL++    -1111               60   MU+        -14
C  28   ADL+     -1121               61   MU-         14
C  29   ADL-     -2221               62   NUE         11
C  30   ADL0     -1221               63   NUM         13
C  31   AS*+     -1131               64   ANUE       -11
C  32   AS*-     -2231               65   ANUM       -13
C  33   AS*0     -1231
C
C The SHIELD's numbering scheme:
C
C JPART   PARTICLE         JPART   PARTICLE
C     1   NEUTRON             13   ELECTRON
C     2   PROTON              14   POSITRON
C     3   PI-                 15   MUON-
C     4   PI+                 16   MUON+
C     5   PI0                 17   NU E
C     6   ANTI N              18   ANTI NU E
C     7   ANTI P              19   NU MU
C     8   K-                  20   ANTI NU MU
C     9   K+                  21   DEUTERON
C    10   K0                  22   TRITON
C    11   K0 TILDA            23   HE3
C    12   GAMMA               24   ALPHA
C............................................................
      IF(IDP.eq.0)THEN      ! Input check
        JPART=0
        CONTINUE!WRITE(25,1000)
 1000     FORMAT(' IDPJPR DIAGNOSTIC: Input IDP=0')
        RETURN
      END IF
C
C Initial values: NN=NN(IDP) and JPART=0
      do j=1,65
        if(IDP .eq. IDENTP(J))then
            NN=J
            go to 1
        end if
      end do
    1 JPART=0
C
C--------------- Conversion AGT --> SHIELD -------------------
      IF(IDP .EQ.  120)THEN   !  [1]   PI+  
        JPART= 4
            NCONV( 1)=NCONV( 1)+1
            JCONV( 1)=JPART
            ECONV( 1)=ECONV( 1)+TKMeV
        RETURN
      END IF
C
      IF(IDP .EQ. -120)THEN   !  [2]   PI-  
        JPART= 3
            NCONV( 2)=NCONV( 2)+1
            JCONV( 2)=JPART
            ECONV( 2)=ECONV( 2)+TKMeV
        RETURN
      END IF
C
      IF(IDP .EQ.  130)THEN   !  [3]   K+   
        JPART= 9
            NCONV( 3)=NCONV( 3)+1
            JCONV( 3)=JPART
            ECONV( 3)=ECONV( 3)+TKMeV
        RETURN
      END IF
C
      IF(IDP .EQ. -130)THEN   !  [4]   K-   
        JPART= 8
            NCONV( 4)=NCONV( 4)+1
            JCONV( 4)=JPART
            ECONV( 4)=ECONV( 4)+TKMeV
        RETURN
      END IF
C
      IF(IDP .EQ.  230)THEN   !  [5]   K0  
        JPART=10
            NCONV( 5)=NCONV( 5)+1
            JCONV( 5)=JPART
            ECONV( 5)=ECONV( 5)+TKMeV
        RETURN
      END IF
C
      IF(IDP .EQ. -230)THEN   !  [6]   AK0  
        JPART=11
            NCONV( 6)=NCONV( 6)+1
            JCONV( 6)=JPART
            ECONV( 6)=ECONV( 6)+TKMeV
        RETURN
      END IF
C
      IF(IDP .EQ.  110)THEN   !  [7]   PI0  
        JPART= 5
            NCONV( 7)=NCONV( 7)+1
            JCONV( 7)=JPART
            ECONV( 7)=ECONV( 7)+TKMeV
        RETURN
      END IF
C
      IF(IDP .EQ.  220)THEN   !  [8]   ETA   is converted to PI0
        JPART= 5
            NCONV( 8)=NCONV( 8)+1
            JCONV( 8)=JPART
            ECONV( 8)=ECONV( 8)+TKMeV
        GO TO 10
      END IF
C
      IF(IDP .EQ.  330)THEN   !  [9]   ETAP  is converted to PI0
        JPART= 5
            NCONV( 9)=NCONV( 9)+1
            JCONV( 9)=JPART
            ECONV( 9)=ECONV( 9)+TKMeV
        GO TO 10
      END IF
C
      IF(IDP .EQ.  121)THEN   ! [10]   RHO+  is converted to PI+
        JPART= 4
            NCONV(10)=NCONV(10)+1
            JCONV(10)=JPART
            ECONV(10)=ECONV(10)+TKMeV
        GO TO 10
      END IF
C
      IF(IDP .EQ. -121)THEN   ! [11]   RHO-  is converted to PI-
        JPART= 3
            NCONV(11)=NCONV(11)+1
            JCONV(11)=JPART
            ECONV(11)=ECONV(11)+TKMeV
        GO TO 10
      END IF
C
      IF(IDP .EQ.  131)THEN   ! [12]   K*+   is converted to K+
        JPART= 9
            NCONV(12)=NCONV(12)+1
            JCONV(12)=JPART
            ECONV(12)=ECONV(12)+TKMeV
        GO TO 10
      END IF
C
      IF(IDP .EQ. -131)THEN   ! [13]   K*-   is converted to K-
        JPART= 8
            NCONV(13)=NCONV(13)+1
            JCONV(13)=JPART
            ECONV(13)=ECONV(13)+TKMeV
        GO TO 10
      END IF
C
      IF(IDP .EQ.  231)THEN   ! [14]   K*0   is converted to K0
        JPART=10
            NCONV(14)=NCONV(14)+1
            JCONV(14)=JPART
            ECONV(14)=ECONV(14)+TKMeV
        GO TO 10
      END IF
C
      IF(IDP .EQ. -231)THEN   ! [15]   AK*0  is converted to AK0
        JPART=11
            NCONV(15)=NCONV(15)+1
            JCONV(15)=JPART
            ECONV(15)=ECONV(15)+TKMeV
        GO TO 10
      END IF
C
      IF(IDP .EQ.  111)THEN   ! [16]   RHO0  is converted to PI0
        JPART= 5
            NCONV(16)=NCONV(16)+1
            JCONV(16)=JPART
            ECONV(16)=ECONV(16)+TKMeV
        GO TO 10
      END IF
C
      IF(IDP .EQ.  221)THEN   ! [17]   OMEG  is converted to PI0
        JPART= 5
            NCONV(17)=NCONV(17)+1
            JCONV(17)=JPART
            ECONV(17)=ECONV(17)+TKMeV
        GO TO 10
      END IF
C
      IF(IDP .EQ.  331)THEN   ! [18]   PHI   is converted to PI0
        JPART= 5
            NCONV(18)=NCONV(18)+1
            JCONV(18)=JPART
            ECONV(18)=ECONV(18)+TKMeV
        GO TO 10
      END IF
C
      IF(IDP .EQ.-1120)THEN   ! [19]   AP   
        JPART= 7
            NCONV(19)=NCONV(19)+1
            JCONV(19)=JPART
            ECONV(19)=ECONV(19)+TKMeV
        RETURN
      END IF
C
      IF(IDP .EQ.-1220)THEN   ! [20]   AN   
        JPART= 6
            NCONV(20)=NCONV(20)+1
            JCONV(20)=JPART
            ECONV(20)=ECONV(20)+TKMeV
        RETURN
      END IF
C
      IF(IDP .EQ.-1130)THEN   ! [21]   AS+   is converted to AP
        JPART= 7
            NCONV(21)=NCONV(21)+1
            JCONV(21)=JPART
            ECONV(21)=ECONV(21)+TKMeV
        GO TO 10
      END IF
C
      IF(IDP .EQ.-2230)THEN   ! [22]   AS-   is converted to AN
        JPART= 6
            NCONV(22)=NCONV(22)+1
            JCONV(22)=JPART
            ECONV(22)=ECONV(22)+TKMeV
        GO TO 10
      END IF
C
      IF(IDP .EQ.-1230)THEN   ! [23]   AS0   is converted to AN
        JPART= 6
            NCONV(23)=NCONV(23)+1
            JCONV(23)=JPART
            ECONV(23)=ECONV(23)+TKMeV
        GO TO 10
      END IF
C
      IF(IDP .EQ.-2330)THEN   ! [24]   AXI-  is ignored
        JPART= 0
            NCONV(24)=NCONV(24)+1
            JCONV(24)=JPART
            ECONV(24)=ECONV(24)+TKMeV
        GO TO 10
      END IF
C
      IF(IDP .EQ.-1330)THEN   ! [25]   AXI0  is ignored
        JPART= 0
            NCONV(25)=NCONV(25)+1
            JCONV(25)=JPART
            ECONV(25)=ECONV(25)+TKMeV
        GO TO 10
      END IF
C
      IF(IDP .EQ.-2130)THEN   ! [26]   AL    is converted to antiN
        JPART= 6
            NCONV(26)=NCONV(26)+1
            JCONV(26)=JPART
            ECONV(26)=ECONV(26)+TKMeV
        GO TO 10
      END IF
C
      IF(IDP .EQ.-1111)THEN   ! [27]   ADL++ is ignored
        JPART= 0
            NCONV(27)=NCONV(27)+1
            JCONV(27)=JPART
            ECONV(27)=ECONV(27)+TKMeV
        GO TO 10
      END IF
C
      IF(IDP .EQ.-1121)THEN   ! [28]   ADL+  is ignored
        JPART= 0
            NCONV(28)=NCONV(28)+1
            JCONV(28)=JPART
            ECONV(28)=ECONV(28)+TKMeV
        GO TO 10
      END IF
C
      IF(IDP .EQ.-2221)THEN   ! [29]   ADL-  is ignored
        JPART= 0
            NCONV(29)=NCONV(29)+1
            JCONV(29)=JPART
            ECONV(29)=ECONV(29)+TKMeV
        GO TO 10
      END IF
C
      IF(IDP .EQ.-1221)THEN   ! [30]   ADL0  is ignored
        JPART= 0
            NCONV(30)=NCONV(30)+1
            JCONV(30)=JPART
            ECONV(30)=ECONV(30)+TKMeV
        GO TO 10
      END IF
C
      IF(IDP .EQ.-1131)THEN   ! [31]   AS*+  is ignored
        JPART= 0
            NCONV(31)=NCONV(31)+1
            JCONV(31)=JPART
            ECONV(31)=ECONV(31)+TKMeV
        GO TO 10
      END IF
C
      IF(IDP .EQ.-2231)THEN   ! [32]   AS*-  is ignored
        JPART= 0
            NCONV(32)=NCONV(32)+1
            JCONV(32)=JPART
            ECONV(32)=ECONV(32)+TKMeV
        GO TO 10
      END IF
C
      IF(IDP .EQ.-1231)THEN   ! [33]   AS*0  is ignored
        JPART= 0
            NCONV(33)=NCONV(33)+1
            JCONV(33)=JPART
            ECONV(33)=ECONV(33)+TKMeV
        GO TO 10
      END IF
C
      IF(IDP .EQ.-2331)THEN   ! [34]   AXI*- is ignored
        JPART= 0
            NCONV(34)=NCONV(34)+1
            JCONV(34)=JPART
            ECONV(34)=ECONV(34)+TKMeV
        GO TO 10
      END IF
C
      IF(IDP .EQ.-1331)THEN   ! [35]   AXI*0 is ignored
        JPART= 0
            NCONV(35)=NCONV(35)+1
            JCONV(35)=JPART
            ECONV(35)=ECONV(35)+TKMeV
        GO TO 10
      END IF
C
      IF(IDP .EQ.-3331)THEN   ! [36]   AOM-  is ignored
        JPART= 0
            NCONV(36)=NCONV(36)+1
            JCONV(36)=JPART
            ECONV(36)=ECONV(36)+TKMeV
        GO TO 10
      END IF
C
      IF(IDP .EQ. 1120)THEN   ! [37]   P    
        JPART= 2
            NCONV(37)=NCONV(37)+1
            JCONV(37)=JPART
            ECONV(37)=ECONV(37)+TKMeV
        RETURN
      END IF
C
      IF(IDP .EQ. 1220)THEN   ! [38]   N    
        JPART= 1
            NCONV(38)=NCONV(38)+1
            JCONV(38)=JPART
            ECONV(38)=ECONV(38)+TKMeV
        RETURN
      END IF
C
      IF(IDP .EQ. 1130)THEN   ! [39]   S+    is converted to P
        JPART= 2
            NCONV(39)=NCONV(39)+1
            JCONV(39)=JPART
            ECONV(39)=ECONV(39)+TKMeV
        GO TO 10
      END IF
C
      IF(IDP .EQ. 2230)THEN   ! [40]   S-    is converted to N
        JPART= 1
            NCONV(40)=NCONV(40)+1
            JCONV(40)=JPART
            ECONV(40)=ECONV(40)+TKMeV
        GO TO 10
      END IF
C
      IF(IDP .EQ. 1230)THEN   ! [41]   S0    is converted to N
        JPART= 1
            NCONV(41)=NCONV(41)+1
            JCONV(41)=JPART
            ECONV(41)=ECONV(41)+TKMeV
        GO TO 10
      END IF
C
      IF(IDP .EQ. 2330)THEN   ! [42]   XI-   is ignored
        JPART= 0
            NCONV(42)=NCONV(42)+1
            JCONV(42)=JPART
            ECONV(42)=ECONV(42)+TKMeV
        GO TO 10
      END IF
C
      IF(IDP .EQ. 1330)THEN   ! [43]   XI0   is ignored
        JPART= 0
            NCONV(43)=NCONV(43)+1
            JCONV(43)=JPART
            ECONV(43)=ECONV(43)+TKMeV
        GO TO 10
      END IF
C
      IF(IDP .EQ. 2130)THEN   ! [44]   L     is converted to N
        JPART= 1
            NCONV(44)=NCONV(44)+1
            JCONV(44)=JPART
            ECONV(44)=ECONV(44)+TKMeV
        GO TO 10
      END IF
C
      IF(IDP .EQ. 1111)THEN   ! [45]   DL++  is ignored
        JPART= 0
            NCONV(45)=NCONV(45)+1
            JCONV(45)=JPART
            ECONV(45)=ECONV(45)+TKMeV
        GO TO 10
      END IF
C
      IF(IDP .EQ. 1121)THEN   ! [46]   DL+   is ignored
        JPART= 0
            NCONV(46)=NCONV(46)+1
            JCONV(46)=JPART
            ECONV(46)=ECONV(46)+TKMeV
        GO TO 10
      END IF
C
      IF(IDP .EQ. 2221)THEN   ! [47]   DL-   is ignored
        JPART= 0
            NCONV(47)=NCONV(47)+1
            JCONV(47)=JPART
            ECONV(47)=ECONV(47)+TKMeV
        GO TO 10
      END IF
C
      IF(IDP .EQ. 1221)THEN   ! [48]   DL0   is ignored
        JPART= 0
            NCONV(48)=NCONV(48)+1
            JCONV(48)=JPART
            ECONV(48)=ECONV(48)+TKMeV
        GO TO 10
      END IF
C
      IF(IDP .EQ. 1131)THEN   ! [49]   S*+   is ignored
        JPART= 0
            NCONV(49)=NCONV(49)+1
            JCONV(49)=JPART
            ECONV(49)=ECONV(49)+TKMeV
        GO TO 10
      END IF
C
      IF(IDP .EQ. 2231)THEN   ! [50]   S*-   is ignored
        JPART= 0
            NCONV(50)=NCONV(50)+1
            JCONV(50)=JPART
            ECONV(50)=ECONV(50)+TKMeV
        GO TO 10
      END IF
C
      IF(IDP .EQ. 1231)THEN   ! [51]   S*0   is ignored
        JPART= 0
            NCONV(51)=NCONV(51)+1
            JCONV(51)=JPART
            ECONV(51)=ECONV(51)+TKMeV
        GO TO 10
      END IF
C
      IF(IDP .EQ. 2331)THEN   ! [52]   XI*-  is ignored
        JPART= 0
            NCONV(52)=NCONV(52)+1
            JCONV(52)=JPART
            ECONV(52)=ECONV(52)+TKMeV
        GO TO 10
      END IF
C
      IF(IDP .EQ. 1331)THEN   ! [53]   XI*0  is ignored
        JPART= 0
            NCONV(53)=NCONV(53)+1
            JCONV(53)=JPART
            ECONV(53)=ECONV(53)+TKMeV
        GO TO 10
      END IF
C
      IF(IDP .EQ. 3331)THEN   ! [54]   OM-   is ignored
        JPART= 0
            NCONV(54)=NCONV(54)+1
            JCONV(54)=JPART
            ECONV(54)=ECONV(54)+TKMeV
        GO TO 10
      END IF
C
      IF(IDP .EQ.   20)THEN   ! [55]   KS    is ignored
        JPART= 0
            NCONV(55)=NCONV(55)+1
            JCONV(55)=JPART
            ECONV(55)=ECONV(55)+TKMeV
        GO TO 10
      END IF
C
      IF(IDP .EQ.  -20)THEN   ! [56]   KL    is ignored
        JPART= 0
            NCONV(56)=NCONV(56)+1
            JCONV(56)=JPART
            ECONV(56)=ECONV(56)+TKMeV
        GO TO 10
      END IF
C
      IF(IDP .EQ.   10)THEN   ! [57]   GM   
        JPART=12
            NCONV(57)=NCONV(57)+1
            JCONV(57)=JPART
            ECONV(57)=ECONV(57)+TKMeV
        RETURN
      END IF
C
      IF(IDP .EQ.   12)THEN   ! [58]   E+   
        JPART=14
            NCONV(58)=NCONV(58)+1
            JCONV(58)=JPART
            ECONV(58)=ECONV(58)+TKMeV
        RETURN
      END IF
C
      IF(IDP .EQ.  -12)THEN   ! [59]   E-   
        JPART=13
            NCONV(59)=NCONV(59)+1
            JCONV(59)=JPART
            ECONV(59)=ECONV(59)+TKMeV
        RETURN
      END IF
C
      IF(IDP .EQ.  -14)THEN   ! [60]   MU+  
        JPART=16
            NCONV(60)=NCONV(60)+1
            JCONV(60)=JPART
            ECONV(60)=ECONV(60)+TKMeV
        RETURN
      END IF
C
      IF(IDP .EQ.   14)THEN   ! [61]   MU-  
        JPART=15
            NCONV(61)=NCONV(61)+1
            JCONV(61)=JPART
            ECONV(61)=ECONV(61)+TKMeV
        RETURN
      END IF
C
      IF(IDP .EQ.   11)THEN   ! [62]   NUE  
        JPART=17
            NCONV(62)=NCONV(62)+1
            JCONV(62)=JPART
            ECONV(62)=ECONV(62)+TKMeV
        RETURN
      END IF
C
      IF(IDP .EQ.   13)THEN   ! [63]   NUM  
        JPART=19
            NCONV(63)=NCONV(63)+1
            JCONV(63)=JPART
            ECONV(63)=ECONV(63)+TKMeV
        RETURN
      END IF
C
      IF(IDP .EQ.  -11)THEN   ! [64]   ANUE 
        JPART=18
            NCONV(64)=NCONV(64)+1
            JCONV(64)=JPART
            ECONV(64)=ECONV(64)+TKMeV
        RETURN
      END IF
C
      IF(IDP .EQ.  -13)THEN   ! [65]   ANUM 
        JPART=20
            NCONV(65)=NCONV(65)+1
            JCONV(65)=JPART
            ECONV(65)=ECONV(65)+TKMeV
        RETURN
      END IF
C
C ----------------- end of conversion ----------------------
c   10 CONTINUE
   10 RETURN
C ------------------------ Debug printout -----------------
      if(JPART.ne.0)then
C         AGT's particle is converted
        CONTINUE!WRITE(25,11)PNAME(NN),PARTIC(JPART),TKMeV
   11     format(' IDP-->JPART:       ',A8,' IS CONVERTED TO ',A7,
     *           ', E(MeV)=',F10.1)
        return
      else
C         AGT's particle is ignored
        CONTINUE!WRITE(25,21)PNAME(NN),IDP,TKMeV
   21     format(' IDP-->JPART:       ',A8,' IDP=',I6,
     *           ' is ignored, E(MeV)=',F10.1)
        return
      end if
C
      END
      
      
      
      SUBROUTINE PRESPT(KSMEM1,KSMEMO,KSTART,INOVER,PNX,PNY,PNZ)
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C Conversion of PRECO secondary particles in /MEMAGT/ into SHIELD
C secondary particles in /SECHAR/.
C    This subroutine send the precompound particles from /MEMAGT/
C    beginning from line No. KSMEM1 to line No. KSMEMO-1
C    into array SPT beginning from line No. KSTART.
C INPUT: 
C       KSMEM1 - No. of first PRECO particle in /MEMAGT/
C       KSMEMO - No. of first free line in /MEMAGT/
C  i.e. KSMEMO-KSMEM1 is number of PRECO perticles in /MEMAGT/
C       KSTART - No. of first free line in SPT
C       PNX,PNY,PNZ - Nucleus momentum (GeV/c) for star printout only !
C INPUT and OUTPUT: 
C       KSTART - New No. of first free line in SPT
C       INOVER=1 IF ARRAY SPT WAS EXCEEDED, ELSE INOVER=0.
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      COMMON /SECHAR/ SPT(6,5000),SNU(10,101),LS6,LS100,LS10,LS11  ! GSITR
      COMMON /MEMAGT/ PMEMO(9,5999),IMEMO(5,5999)
      REAL*8 PMEMO
      COMMON /BK1003/ UPRECO,APRECO,ZPRECO     ! for star printout only !
      REAL*4 UPRECO,APRECO,ZPRECO              ! UPRECO in MeV !
      COMMON /DEBUG/ LSTAR,LCASC
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
C
      NPREC=KSMEMO-KSMEM1    ! Number of precompound particles
      KFINIS=KSTART+NPREC    ! New No. of first free line in SPT
C
      IF(KFINIS.GE.LS100-1)THEN
         INOVER=1
         CONTINUE!WRITE(25,1000)KFINIS,IXINIT,NUMTRE
 1000      FORMAT(' PRESPT DIAGNOSTIC: OVERFILLING OF SPT',
     *     'KFINIS=',I4,'  IXINIT=',I12,'  NUMTRE=',I6)
         KFINIS=LS100-1
      END IF
C
      DO J=1,NPREC
            JP=KSMEM1+J-1
            JS=KSTART+J-1
        PX=SNGL(PMEMO(4,JP))
        PY=SNGL(PMEMO(5,JP))
        PZ=SNGL(PMEMO(6,JP))
      CALL DECARS(PX,PY,PZ,PMOD,COSTP,SINFP,COSFP)
        SPT(1,JS)=COSTP
        SPT(2,JS)=SINFP
        SPT(3,JS)=COSFP
        SPT(4,JS)=SNGL(PMEMO(8,JP))*1000.
        SPT(5,JS)=FLOAT(IMEMO(1,JP))
        SPT(6,JS)=SNGL(PMEMO(9,JP))*1000.
      END DO
C
cC Debug printout ==============================================
C
      IF(LSTAR.eq.2)THEN
      IF(NPREC.gt.0)CONTINUE!WRITE(36,1)
    1 FORMAT(' MEMAGT after PRECO')
C
      DO J=KSMEM1,KSMEMO-1
      CONTINUE!WRITE(36,2)J,PMEMO(4,J)*1000.,PMEMO(5,J)*1000.,PMEMO(6,J)*1000.,
!     *               PMEMO(8,J)*1000.,PMEMO(9,J)*1000.,  ! GeV ==> MeV
!     *               IMEMO(1,J),IMEMO(3,J),IMEMO(4,J)
    2   FORMAT(I3,1X,3E11.4,E13.4,F11.1,3X,3I2)
      END DO
C
      PX1=PNX*1000.
      PY1=PNY*1000.   ! Momentum: GeV ==> MeV
      PZ1=PNZ*1000.
        A1=APRECO
        TNUCL1=SQRT(PX1**2+PY1**2+PZ1**2+(940.*A1)**2)-(940.*A1)
      IF(NPREC.gt.0)CONTINUE!WRITE(36,3)
    3 FORMAT(1X,49('-'))
      IF(NPREC.gt.0)CONTINUE!WRITE(36,4)PX1,PY1,PZ1,TNUCL1,APRECO,ZPRECO,UPRECO
    4 FORMAT('RES ',3E11.4,E13.4,2X,2F5.1,E11.4)
      END IF
C
C ===================End debug printout========================
C
      KSTART=KFINIS
C
      RETURN
      END
 
 
 
      SUBROUTINE TYPART(ICHARG,RMASS,JPART)
C DETERMINATION OF PARTICLE TYPE JPART ACCORDING TO PARTICLE CHARGE AND MASS.
C INPUT: ICHARG - CHARGE OF PARTICLE
C        RMASS  - MASS OF PARTICLE (MEV)
C OUTPUT: JPART - PARTICLE TYPE
C
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
C
      IF(ABS(RMASS- 140.).LT.10.)GO TO 10  ! pions              
      IF(ABS(RMASS- 495.).LT.10.)GO TO 15  ! kaons
      IF(ABS(RMASS- 940.).LT.10.)GO TO 20  ! nucleons and antinucleons
      IF(ABS(RMASS-1880.).LT.10.)GO TO 30  ! deuterons
      IF(ABS(RMASS-2820.).LT.10.)GO TO 40  ! Tritons and He3
      IF(ABS(RMASS-3760.).LT.10.)GO TO 50  ! alphas
      IF(ABS(RMASS-12.).LT.1. .and. ICHARG.eq.12)GO TO 60  ! GAMMA
      GO TO 1000
C
   10 IF(ICHARG.LT.-1 .OR. ICHARG.GT. 1)GO TO 1000
      IF(ICHARG.EQ.-1)JPART=3
      IF(ICHARG.EQ. 0)JPART=5
      IF(ICHARG.EQ. 1)JPART=4
      RETURN
C
   15 IF(ICHARG.LT.-1 .OR. ICHARG.GT.11)GO TO 1000
      IF(ICHARG.EQ.-1)JPART=8
      IF(ICHARG.EQ. 0)JPART=10
      IF(ICHARG.EQ. 1)JPART=9
      IF(ICHARG.EQ.11)JPART=11
      RETURN
C
   20 IF(ICHARG.LT.-1 .OR. ICHARG.GT. 6)GO TO 1000
      IF(ICHARG.EQ.0) JPART=1
      IF(ICHARG.EQ.1) JPART=2
      IF(ICHARG.EQ.-1)JPART=7
      IF(ICHARG.EQ.6) JPART=6
      RETURN
C
   30 IF(ICHARG.EQ.1)THEN
        JPART=21
        RETURN
      ELSE
        GO TO 1000
      END IF
C
   40 IF(ICHARG.LT.1 .OR. ICHARG.GT.2)GO TO 1000
      IF(ICHARG.EQ.1)JPART=22
      IF(ICHARG.EQ.2)JPART=23
      RETURN      
C
   50 IF(ICHARG.EQ.2)THEN
        JPART=24
        RETURN
      ELSE
        GO TO 1000
      END IF
C
   60 JPART=12   ! Gamma
      RETURN
C
 1000 CONTINUE!WRITE(25,1001)ICHARG,RMASS,IXINIT,NUMTRE
 1001 FORMAT(' TYPART DIAGNOSTIC: INVALID PARTICLE, ICHARG=',I3,
     *'  RMASS=',F7.1/20X,'IXINIT=',I12,'  NUMTRE=',I6,'   STOP')
      STOP
      END
 

 
      SUBROUTINE DECARS(PX,PY,PZ,PMOD,COST,SINF,COSF)
C CALCULATION OF SPHERICAL CUOORDINATES COST,SINF,COSF OF VECTOR PX,PY,PZ.
C INPUT: PX,PY,PZ
C OUTPUT: PMOD,COST,SINF,COSF
      REAL*8 DPX,DPY,DPZ,DP,DCOST,DSINF,DCOSF,DSINT,DTEMP
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
C
      DPX=DBLE(PX)
      DPY=DBLE(PY)
      DPZ=DBLE(PZ)
      DP=DSQRT(DPX**2+DPY**2+DPZ**2)
        IF(DP.LE.0.)GO TO 2
C
      DCOST=DPZ/DP
      DTEMP=1.0-DCOST**2
        IF(DTEMP.LE.0.)GO TO 1
C
      DSINT=DSQRT(DTEMP)
        DTEMP=DP*DSINT
      DSINF=DPY/DTEMP
      DCOSF=DPX/DTEMP
            COST=SNGL(DCOST)
            SINF=SNGL(DSINF)
            COSF=SNGL(DCOSF)
            PMOD=SNGL(DP)
      RETURN
C
    1 IF(DCOST.LT.0)THEN
        COST=-1.
      ELSE
        COST=+1.
      END IF
            SINF=0.707 106 781 187
            COSF=0.707 106 781 187
            PMOD=SNGL(DP)
      RETURN
C
    2 COST=1.-2.*RRAN(IX)
        FI=6.283 185 307 180*RRAN(IX)
      SINF=SIN(FI)
      COSF=COS(FI)
      PMOD=0.
      RETURN
      END
 
 
      SUBROUTINE SUBROT(CT,SF,CF,CTP,SFP,CFP,CTR,SFR,CFR,STR)
C ROTATION OF SECONDARY PARTICLE DIRECTION
C INPUT: CT,SF,CF - DIRECTION OF PRIMARY PARTICLE IN LABORATORY SYSTEM.
C        CTP,SFP,CFP - DIRECTION OF SECONDARY PARTICLE IN SYSTEM,
C                      WHERE PRIMARY PARTICLE FLY ALONG AXIS Z.
C OUTPUT: CTR,SFR,CFR,STR - DIRECTION OF SECONDARY PARTICLE IN INITIAL
C                           LABORATORY SYSTEM.
C--------------------------------------------------------------------
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      REAL*8 DCT,DSF,DCF,DST,DCTP,DSFP,DCFP,DSTP
      REAL*8 DCTR,DSFR,DCFR,DSTR,DTEMP
C
      DCT=DBLE(CT)
      DSF=DBLE(SF)
      DCF=DBLE(CF)
      DST=DSQRT(1.0-DCT**2)
C
      DCTP=DBLE(CTP)
      DSFP=DBLE(SFP)
      DCFP=DBLE(CFP)
      DSTP=DSQRT(1.0-DCTP**2)
C
      DCTR=DCTP*DCT-DSTP*DCFP*DST
        DTEMP=1.0-DCTR**2
        IF(DTEMP.LE.0.)GO TO 1
C
      DSTR=DSQRT(DTEMP)
      DSFR=(DSTP*DCFP*DCT*DSF+DSTP*DSFP*DCF+DCTP*DST*DSF)/DSTR
      DCFR=(DSTP*DCFP*DCT*DCF-DSTP*DSFP*DSF+DCTP*DST*DCF)/DSTR
        CTR=SNGL(DCTR)
        SFR=SNGL(DSFR)
        CFR=SNGL(DCFR)
        STR=SNGL(DSTR)
      RETURN
C
    1 IF(DCTR.LT.0.)THEN
        CTR=-1.
      ELSE
        CTR=+1.
      END IF
        SFR=0.707 106 781 187
        CFR=0.707 106 781 187
        STR=0.
      RETURN
      END
 
 
 
C -------rotor of Dementyev -----------------------------
      SUBROUTINE ROTORD(SAR,SBR,SPSTAR,SPR)    ! ROTORD - 13.02.96
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
C    BLOCK OF ROTATION.
      DIMENSION SAR(3),SBR(3),SPSTAR(3),SPR(3)
      REAL*4 SAR,SBR,SPSTAR,SPR
      DIMENSION AR(3),BR(3),PSTAR(3),PR(3),AN(3)
C
       DO I=1,3
      AR(I)=DBLE(SAR(I))
      BR(I)=DBLE(SBR(I))
      PSTAR(I)=DBLE(SPSTAR(I))
       ENDDO
C
      SP = 0.
       DO IR=1,3
      SP = SP+AR(IR)*BR(IR)
       ENDDO
      AMOD = DSQRT(AR(1)**2+AR(2)**2+AR(3)**2)
      ALPHA1 = SP/AMOD
      BMOD2 = BR(1)**2+BR(2)**2+BR(3)**2
C
      BMAL=BMOD2-ALPHA1**2
      IF(BMAL.LE.0.0)THEN
      PRINT *,'ROTORD: AR,BR,AMOD,BMOD2,ALPHA1,ALPHA1**2,BMAL',      
     *                   AR,BR,AMOD,BMOD2,ALPHA1,ALPHA1**2,BMAL
      ENDIF
      ALPHA2=DSQRT(BMAL)
C
C      ALPHA2 = DSQRT(BMOD2-ALPHA1**2)
      AN(1) = AR(2)*BR(3)-AR(3)*BR(2)
      AN(2) = AR(3)*BR(1)-AR(1)*BR(3)
      AN(3) = AR(1)*BR(2)-AR(2)*BR(1)
      PR(1)=PSTAR(1)*BR(1)/ALPHA2+(PSTAR(3)-ALPHA1*PSTAR(1)/ALPHA2)
     1*AR(1)/AMOD+(PSTAR(2)*AN(1))/(ALPHA2*AMOD)
      PR(2)=PSTAR(1)*BR(2)/ALPHA2+(PSTAR(3)-ALPHA1*PSTAR(1)/ALPHA2)
     1*AR(2)/AMOD+(PSTAR(2)*AN(2))/(ALPHA2*AMOD)
      PR(3)=PSTAR(1)*BR(3)/ALPHA2+(PSTAR(3)-ALPHA1*PSTAR(1)/ALPHA2)
     1*AR(3)/AMOD+(PSTAR(2)*AN(3))/(ALPHA2*AMOD)
C
       DO I=1,3
      SPR(I)=SNGL(PR(I))
       ENDDO
C
      RETURN
      END



c                KOBR=0 -> p+A, KOBR=1 -> A+p     Where is ?
C************************************************************
C ********* BELOW - INTERACTION WITH HYDROGEN ONLY **********
C
      SUBROUTINE CNAGTH(PM1,IM1,PM2,IM2)
C  CONVERSION OF EXTRANUCLEAR VARIABLES JPART AND NUCLID
C  INTO INPUT VARIABLES AND PARAMETERS FOR
C  AMELIN GENERATOR (HADRON+NUCLEON GENERATOR)
C
C            COMMONS OF H+H -  GENERATOR
C---------------------------------------------------------------------
      REAL*8 ENBOU
      COMMON /H1H2/ GH1H2(11)
      COMMON /NOTRE/ NOTRE
      COMMON /NODCAY/ NODCAY,NOETA,NOPI0,NONUNU,NOEVOL,NOHADR,NOKA0
      COMMON /COMMUL/ MULTP
      COMMON /ITAPES/ ITDKY,ITEVT,ITCOM,ITLIS
      COMMON /COMENB/ ENBOU
      COMMON /COMFR/ ICMS
      COMMON /KEYHH/ KEYHH
      COMMON /COMQSE/ QSEE,QVSEE
      COMMON /KEYPLA/ KEYPLA
      COMMON /COMWTI/ WTIME
      LOGICAL QSEE,QVSEE
      LOGICAL MULTP
      LOGICAL NODCAY,NOETA,NOPI0,NONUNU,NOEVOL,NOHARD,NOKA0
      LOGICAL KEYHH,keypla
      LOGICAL GH1H2
      LOGICAL LPRNT
C Sobolevsky, 06.04.95 ----------
C      LOGICAL ITDKY,ITEVT,ITCOM,ITLIS
C -----------------------------
      LOGICAL NOTRE
      LOGICAL WTIME
      REAL*8 PM1(9),PM2(9)
      DIMENSION IM1(5),IM2(5)
C -----------------------------------------------------------------
C
C                 COMMONS OF SHIELD/CG
      COMMON /INFEL/ ZNUC(110),ATWEI(110),DENS(110),RIEV(110),SYMB(110)
      CHARACTER*2 SYMB
      COMMON /INREAC/ COST,SINF,COSF,TINT,WINT,JPART,NUCLID,KSTATE
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
C
       LPRNT=.FALSE.
         NOTRE=.TRUE.
             KEYHH=.TRUE.
               NODCAY=.FALSE.
                   MULTP=.TRUE.
                     KEYPLA=.FALSE.
                         QSEE =.FALSE.
                           QVSEE =.FALSE.
                               WTIME=.FALSE.
                                    ENBOU= 4.4
      ICMS=0   ! Lab. frame
C      ICMS=1   ! CMS frame, Dementyev 24.06.94
C
      DO 110 I=1,11
110   GH1H2(I)=.TRUE.
      GH1H2(4)=.FALSE.
C
C     INCIDENT PARTICLE:
      T0=DBLE(0.001*TINT)
      IF(JPART .LT. 1 .OR. JPART .GT.11) GOTO 1000
      GO TO (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11),JPART
 1    IM1(1)=0
        IM1(3)=0
            IM1(4)=1
              PM1(9)=0.940
                  GO TO 12
 2    IM1(1)=1
        IM1(3)=0
            IM1(4)=1
              PM1(9)=0.940
                 GO TO 12
 3    IM1(1)=-1
        IM1(3)=0
            IM1(4)=0
              PM1(9)=0.140
                  GO TO 12
 4    IM1(1)=1
        IM1(3)=0
            IM1(4)=0
              PM1(9)=0.140
                  GO TO 12
 5    GO TO 1000
 6    IM1(1)=0
        IM1(3)=0
            IM1(4)=-1
              PM1(9)=0.940
                  GO TO 12
 7    IM1(1)=-1
        IM1(3)=0
            IM1(4)=-1
              PM1(9)=0.940
                  GO TO 12
 8    IM1(1)=-1
        IM1(3)=-1
            IM1(4)=0
              PM1(9)=0.495
                  GO TO 12
 9    IM1(1)=1
        IM1(3)=1
            IM1(4)=0
              PM1(9)=0.495
                  GO TO 12
 10   IM1(1)=0
        IM1(3)=1
            IM1(4)=0
              PM1(9)=0.495
                  GO TO 12
 11   IM1(1)=0
        IM1(3)=-1
            IM1(4)=0
              PM1(9)=0.495
                  GO TO 12
 12   IM1(2)=0
        IM1(5)=0
            PM1(1)=0.
              PM1(2)=0.
                  PM1(3)=0.
                    PM1(4)=0.
                        PM1(5)=0.
                          PM1(6)=DSQRT(T0*(T0+2.*PM1(9)))
                              PM1(7)=0.
                                PM1(8)=T0
C
C   TARGET-NUCLEON:
      IM2(1)=1
        IM2(2)=0
            IM2(3)=0
              IM2(4)=1
                  IM2(5)=0
                    PM2(1)=0.
                        PM2(2)=0.
                          PM2(3)=0.
                              PM2(4)=0.
                                PM2(5)=0.
                                    PM2(6)=1.D-12
                                      PM2(7)=0.
                                    PM2(8)=0.
                                PM2(9)=0.940
      CALL SETCON
      CALL SETDKY(LPRNT)
      RETURN
C
 1000 CONTINUE!WRITE(25,1001) JPART,IXINIT,NUMTRE
 1001 FORMAT('CNAGTH DIAGNOSTIC: INVALID JPART= ',I2,'  IXINIT=',I12,
     * '  NUMTRE=',I4)
      STOP
      END



      SUBROUTINE AGTCNH(KSTART,IELA)
C   CONVERSION OF AGT-GENERATOR SECONDARY PARITCLES IN /PPTCL/ INTO
C   SHIELD/CG SECONDARY PARTICLES IN /SECHAR/
C   INPUT: NPTCL - NUMBER OF PARTICLES IN /PPTCL/
C   OUTPUT: KSTART - THE FIRST FREE LINE IN SPT
C           IELA
C
      COMMON /SECHAR/ SPT(6,5000),SNU(10,101),LS6,LS100,LS10,LS11  ! GSITR
C
      COMMON/PARTCL/PPTCL(9,499),NPTCL,IORIG(499),IDENT(499),IDCAY(499)
      CHARACTER*8 PNA1,PNA3602,PNAJ
      COMMON /INDDEC/ INDDEC
      COMMON /ITHEA/ ITHEA(11)
      COMMON /NODCAY/ NODCAY,NOETA,NOPI0,NONUNU,NOEVOL,NOHARD,NOKA0
      LOGICAL NODCAY,NOETA,NOPI0,NONUNU,NOEVOL,NOHARD,NOKA0
      COMMON /COMLID/ PLIDER(499)
      LOGICAL GH1H2
      REAL*8 PPTCL
C Sobolevsky, 06.04.95 -------------
      real*8 plider
C --------------------
      COMMON /H1H2/ GH1H2(11)
C
      NP=0
      IPOINT=-1
      IELA=0     ! KTH (This statement don't affects the results, 24.03.01)
C
      IF(NPTCL.LE.0)GO TO 401
C
C Decay of unstable particles
      DO 400 ND=1,10
        IF(NPTCL.LE.0)GO TO 401
        NDEC=0
        NDEC2=NPTCL
            DO 360 I=1,NDEC2
              IF(NODCAY)GO TO 400
              CALL DECAY(I,IPOINT)
              IF(IPOINT.LT.0)GO TO 360
              NDEC=NDEC+1
                  DO J=1,9
                    PPTCL(J,I)=PPTCL(J,NPTCL)
                  END DO
              IDENT(I)=IDENT(NPTCL)
              IORIG(I)=IORIG(NPTCL)
              IDCAY(I)=IDCAY(NPTCL)
              PLIDER(I)=PLIDER(NPTCL)
              NPTCL=NPTCL-1
  360         CONTINUE
        IF(NDEC.EQ.0)GO TO 362
  400 CONTINUE
C
  401 CONTINUE
      IF(NPTCL.GT.0)GO TO 362
      IELA=1
      RETURN
C
  362 CONTINUE
      DO I=1,11
        IF(GH1H2(I))ITHEA(I)=1        ! GH1H2 - LOGICAL !
      END DO
C
      KSTART=NPTCL+1
      JJ=0
C
      DO 100 J=1,NPTCL
        JJ=JJ+1
C
C ============== Insert of Sobolevsky, 17.02.96 ================
C Conversion for products of hh-collision
C
      TKMeV=SNGL(PPTCL(4,J)-PPTCL(5,J))*1000.
      IDP=IDENT(J)
      CALL IDPJPR(IDP,IPART,TKMeV)
C
C Is the particle nucleon, pion, antinucleon, kaon or antikaon?
        if(IPART .ne. 0) go to 20
        JJ=JJ-1
        KSTART=KSTART-1
C
      GO TO 100
C                  To the end of cycle (this is other particle)
C
C  THE PARTICLE IS A NUCLEON, PION, ANTINUCLEON, KAON, ANTIKAON or gamma
   20 CONTINUE
        PX=SNGL(PPTCL(1,J))
        PY=SNGL(PPTCL(2,J))
        PZ=SNGL(PPTCL(3,J))
      CALL DECARS(PX,PY,PZ,PMOD,COSTP,SINFP,COSFP)
        SPT(1,JJ)=COSTP
        SPT(2,JJ)=SINFP
        SPT(3,JJ)=COSFP
C Dementyev 24.06.94 ---------------------------
C         SPT(4,JJ)=SNGL(PPTCL(8,J))*1000.
        SPT(4,JJ)=SNGL(PPTCL(4,J)-PPTCL(5,J))*1000.
C-----------------------------------------------
C Appropriation of charges and masses to particles according to IPART. 
C AntiN and antiK0 have got charges 6 and 11 respectively to differ them
C from N and K0. The inverse thansformation of charge and mass to JPART 
C is made by subroutine TYPART.
      IF(IPART.eq.1)THEN
        SPT(5,JJ)=0.
        SPT(6,JJ)=940.
        GO TO 100
      END IF
      IF(IPART.eq.2)THEN
        SPT(5,JJ)=1.
        SPT(6,JJ)=940.
        GO TO 100
      END IF
      IF(IPART.eq.3)THEN
        SPT(5,JJ)=-1.
        SPT(6,JJ)=140.
        GO TO 100
      END IF
      IF(IPART.eq.4)THEN
        SPT(5,JJ)=1.
        SPT(6,JJ)=140.
        GO TO 100
      END IF
      IF(IPART.eq.5)THEN
        SPT(5,JJ)=0.
        SPT(6,JJ)=140.
        GO TO 100
      END IF
      IF(IPART.eq.6)THEN
          SPT(5,JJ)=6.       ! Let it is the charge of ANTI N
        SPT(6,JJ)=940.
        GO TO 100
      END IF
      IF(IPART.eq.7)THEN
        SPT(5,JJ)=-1.
        SPT(6,JJ)=940.
        GO TO 100
      END IF
      IF(IPART.eq.8)THEN
        SPT(5,JJ)=-1.
        SPT(6,JJ)=495.
        GO TO 100
      END IF
      IF(IPART.eq.9)THEN
        SPT(5,JJ)=1.
        SPT(6,JJ)=495.
        GO TO 100
      END IF
      IF(IPART.eq.10)THEN
        SPT(5,JJ)=0.
        SPT(6,JJ)=495.
        GO TO 100
      END IF
      IF(IPART.eq.11)THEN
          SPT(5,JJ)=11.      ! Let it is the charge of ANTI K0
        SPT(6,JJ)=495.
        GO TO 100
      END IF
      IF(IPART.eq.12)THEN
          SPT(5,JJ)=12.      ! Let it is the charge of GAMMA QUANT !
        SPT(6,JJ)=12.      ! Let it is the mass of GAMMA QUANT !
        GO TO 100
      END IF
C
C ============== End of insert of Sobolevsky, 17.02.96 ================
C
 100  CONTINUE
C
      RETURN
      END



      SUBROUTINE AGTH1(P1,I1,P2,I2)
C  HADRON+NECLEON INTERACTION,
C  AMELIN-GUDIMA-TONEEV HADRON-NUCLEON GENERATOR. OUTPUT VARIABLE
C
C                COMMONS OF AGT - GENERATOR
      COMMON /INDDEC/ INDDEC
      COMMON /ITHEA/ ITHEA(11)
      CHARACTER*8 PNA1,PNA2,PNAJ
      REAL*8 PPTCL
      COMMON/PARTCL/PPTCL(9,499),NPTCL,IORIG(499),IDENT(499),IDCAY(499)
      REAL*8 P1(9),P2(9)
      DIMENSION I1(5),I2(5)
C Dementyev 24.06.94 -------------------------------
      REAL*8 AM1,PX1,PY1,PZ1,AM2,PX2,PY2,PZ2,SITO,SIAN,SIEL,SIEX
C --------------------------------------------------
C
C                COMMONS OF SHIELD/GG
      COMMON /INREAC/ COST,SINF,COSF,TINT,WINT,JPART,NUCLID,KSTATE
      COMMON /SECHAR/ SPT(6,5000),SNU(10,101),LS6,LS100,LS10,LS11  ! GSITR
        COMMON /HHSIGM/ SITO,SIEL
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
C
C         CLEANING OF ARRAYS SPT,SNU
C       INOVER=0                       ! KarSjh
       DO 1 J=1,LS100
       DO 1 I=1,LS6
   1   SPT(I,J)=0.
         DO 2 J=1,LS11
         DO 2 I=1,LS10
   2       SNU(I,J)=0.
C
C   CONVERSION  OF  SHIELD/CG VARIABLES INTO INPUT QGSM VARIABLES
      CALL PANUN(P1,I1,IK1)
      CALL PANUN(P2,I2,IK2)
        PX1=P1(4)
            PY1=P1(5)
              PZ1=P1(6)
                  PX2=P2(4)
                    PY2=P2(5)
                        PZ2=P2(6)
                          AM1=P1(9)
                              AM2=P2(9)
      CALL IDPANU(ID1,IK1,PNA1)
      CALL IDPANU(ID2,IK2,PNA2)
C
C   INELASTIC INTERACTION
C
      CALL CROSEC(0,ID1,ID2,PX1,PY1,PZ1,AM1,PX2,PY2,PZ2,AM2,SIEL,0)
      CALL CROSEC(1,ID1,ID2,PX1,PY1,PZ1,AM1,PX2,PY2,PZ2,AM2,SITO,0)
      CALL CROSEC(2,ID1,ID2,PX1,PY1,PZ1,AM1,PX2,PY2,PZ2,AM2,SIEX,0)
      CALL CROSEC(3,ID1,ID2,PX1,PY1,PZ1,AM1,PX2,PY2,PZ2,AM2,SIAN,0)
      CALL COLLHH(ID1,AM1,PX1,PY1,PZ1,ID2,AM2,PX2,PY2,PZ2,SITO,SIAN,
     *SIEL)
C
  121 CONTINUE
C
      RETURN
      END



C ------------------- DEBUG PRINTOUT ------------------------
      SUBROUTINE PRIAGT(MV)
C Debug printout of /MEMAGT/ after CASCAN
      COMMON /INREAC/ COST,SINF,COSF,TINT,WINT,JPART,NUCLID,KSTATE
      COMMON /MEMAGT/ PMEMO(9,5999),IMEMO(5,5999)
      REAL*8 PMEMO
      COMMON /IDPME/ IDPME(5999)     ! added 30.07.95
      COMMON /RESULT/ AN1,AN2,ZN1,ZN2,ENEXT1,ENEXT2,PNUCL1(3),PNUCL2(3),
     *                AMNUC1(3),AMNUC2(3)
      REAL*8 AN1,AN2,ZN1,ZN2,ENEXT1,ENEXT2,PNUCL1,PNUCL2,AMNUC1,AMNUC2
C
      CONTINUE!WRITE(36,10)
   10 FORMAT(/' MEMAGT after CASCAN'/
     *4X,'  Px(MeV/c)  Py(MeV/c)  Pz(MeV/c)    Ekin(MeV)  MASS(MeV)',
     *3X,' E S B   IDP')
C
      PXSUM=0.
      PYSUM=0.
      PZSUM=0.
       ESUM=0.
C Particles
      DO J=1,MV
        PXSUM=PXSUM+PMEMO(4,J)*1000.
        PYSUM=PYSUM+PMEMO(5,J)*1000.
        PZSUM=PZSUM+PMEMO(6,J)*1000.
         ESUM= ESUM+PMEMO(8,J)*1000.
       CONTINUE!WRITE(36,20)J,PMEMO(4,J)*1000.,PMEMO(5,J)*1000.,PMEMO(6,J)*1000.,
!     *               PMEMO(8,J)*1000.,PMEMO(9,J)*1000.,
!     *               IMEMO(1,J),IMEMO(3,J),IMEMO(4,J),IDPME(J)
   20   FORMAT(I3,1X,3E11.4,E13.4,F11.1,3X,3I2,I6)
      END DO
C
      CONTINUE!WRITE(36,30)
   30 FORMAT(1X,49('-'))
      CONTINUE!WRITE(36,40)PXSUM,PYSUM,PZSUM,ESUM
   40 FORMAT('Sum:',3E11.4,E13.4)
C Residuel nuclei
      CONTINUE!WRITE(36,50)
   50 FORMAT('     Excited residues:',33X,'A*   Z*    U*(MeV)')
C
C Residue of target nucleus
      PX2=SNGL(PNUCL2(1)*1000.)
      PY2=SNGL(PNUCL2(2)*1000.)
      PZ2=SNGL(PNUCL2(3)*1000.)
        A2=SNGL(AN2)
        TNUCL2=SQRT(PX2**2+PY2**2+PZ2**2+(940.*A2)**2)-(940.*A2)
      CONTINUE!WRITE(36,60)PX2,PY2,PZ2,TNUCL2,AN2,ZN2,ENEXT2*1000.
   60 FORMAT('TRG ',3E11.4,E13.4,2X,2F5.1,E11.4)
C
C Residue of projectile nucleus
      PX1=SNGL(PNUCL1(1)*1000.)
      PY1=SNGL(PNUCL1(2)*1000.)
      PZ1=SNGL(PNUCL1(3)*1000.)
        A1=SNGL(AN1)
        TNUCL1=SQRT(PX1**2+PY1**2+PZ1**2+(940.*A1)**2)-(940.*A1)
      IF(JPART.ge.21)CONTINUE!WRITE(36,70)PX1,PY1,PZ1,TNUCL1,AN1,ZN1,ENEXT1*1000.!HION2
   70 FORMAT('PRJ ',3E11.4,E13.4,2X,2F5.1,E11.4)
C
      RETURN
      END



      SUBROUTINE PRICOA(MV)
C Debug printout of /MEMAGT/ after COALES
      COMMON /MEMAGT/ PMEMO(9,5999),IMEMO(5,5999)
      REAL*8 PMEMO
      COMMON /IDPME/ IDPME(5999)     ! added 30.07.95
      COMMON /RESULT/ AN1,AN2,ZN1,ZN2,ENEXT1,ENEXT2,PNUCL1(3),PNUCL2(3),
     *                AMNUC1(3),AMNUC2(3)
      REAL*8 AN1,AN2,ZN1,ZN2,ENEXT1,ENEXT2,PNUCL1,PNUCL2,AMNUC1,AMNUC2
C
      CONTINUE!WRITE(36,1)
    1 FORMAT(/' MEMAGT after COALES')
C
C Particles
      DO J=1,MV
      CONTINUE!WRITE(36,2)J,PMEMO(4,J)*1000.,PMEMO(5,J)*1000.,PMEMO(6,J)*1000.,
!     *               PMEMO(8,J)*1000.,PMEMO(9,J)*1000.,
!     *               IMEMO(1,J),IMEMO(3,J),IMEMO(4,J),IDPME(J)
    2   FORMAT(I3,1X,3E11.4,E13.4,F11.1,3X,3I2,I6)
      END DO
C
      RETURN
      END



      SUBROUTINE PRSPTU
C Debug printout of SPTU(9,500) after DEEX, 23.12.96.
      COMMON /BLOKC/ SPTU(10,500)        
C
      IF(SPTU(9,1) .gt. 0.0)CONTINUE!WRITE(36,1)
    1 FORMAT(' SPTU after DEEX')
C      
      DO J=1,500
        IF(SPTU(9,J) .LE. 0.0)GO TO 77   ! ==> go out from the loop
            COST=SPTU(4,J)
            SINT=SQRT(1.0-COST*COST)
            SINF=SPTU(5,J)
            COSF=SPTU(6,J)
              TKIN=SPTU(7,J)*1000.    ! kin. energy in MeV
              RESU=SPTU(10,J)
                  AA=SPTU(9,J)*1000.  ! mass in MeV
                  ZZ=SPTU(8,J)
C
        PMOD=SQRT(TKIN*(TKIN+2.0*AA))
          PX=PMOD*SINT*COSF
          PY=PMOD*SINT*SINF
          PZ=PMOD*COST
            IA=NINT(AA/940.)
            IZ=INT(ZZ)
C
C---------Is it a particle or fragment in line No. J?
        IF((IA.EQ.1 .AND. (IZ.EQ.0 .OR. IZ.EQ.1)) .OR.
     *       (IA.EQ.2 .AND. IZ.EQ.1) .OR.
     *       (IA.EQ.3 .AND. (IZ.EQ.1 .OR. IZ.EQ.2)) .OR.
     *       (IA.EQ.4 .AND. IZ.EQ.2)) THEN
C ------------It is particle
            CONTINUE!WRITE(36,10)J,PX,PY,PZ,TKIN,AA,IZ,0,IA
   10         FORMAT(I3,1X,3E11.4,E13.4,F11.1,3X,3I2)
        ELSE
C ------------It is fragment
            CONTINUE!WRITE(36,12)PX,PY,PZ,TKIN,FLOAT(IA),FLOAT(IZ),RESU
   12         FORMAT('RES ',3E11.4,E13.4,2X,2F5.1,E11.4)
        END IF
C
      END DO
C
   77 RETURN
      END



      SUBROUTINE PRIHH
C DEBUG PRINTOUT OF /PARTCL/
      COMMON/PARTCL/PPTCL(9,499),NPTCL,IORIG(499),IDENT(499),IDCAY(499)
      REAL*8 PPTCL
C
      CONTINUE!WRITE(36,1)
    1 FORMAT(' PPTCL after AGTH1'/
     *4X,'  Px(MeV/c)  Py(MeV/c)  Pz(MeV/c)    Ekin(MeV)  MASS(MeV)',
     *' IORG IDNT IDCY')
C
      DO J=1,NPTCL
      CONTINUE!WRITE(36,2)J,pptcl(1,J)*1000.,pptcl(2,J)*1000.,pptcl(3,J)*1000.,
!     *               (pptcl(4,J)-pptcl(5,J))*1000.,pptcl(5,J)*1000.,
!     *               IORIG(J),IDENT(J),IDCAY(J)
    2   FORMAT(I3,1X,3E11.4,E13.4,E11.4,3I5)
      END DO
C
      RETURN
      END



      SUBROUTINE STAR(NEVENT)
C Printout of inelastic interaction star
C       SECONDARY PARTICLES CHARACTERISTICS IN ARRAY SPT
C           SPT(1,J) - COS OF POLAR ANGLE OF SECONDARY PARTICLE
C           SPT(2,J), SPT(3,J) - SIN, COS OF AZIMUTAL ANGLE
C           SPT(4,J) - KINETIC ENERGY OF SECONDARY PARTICLE (MEV)
C           SPT(5,J) - SECONDARY PARTICLE TYPE
C           SPT(6,J) - SECONDARY PARTICLE WEIGHT
C       NUCLEUS CHARACTERISTICS
C           SNU(1,1) - NUCLEUS-TARGET ATOMIC WEIGHT
C           SNU(2,1) - NUCLEUS-TARGET ATOMIC NUMBER
C       Nuclei after deexitation
C           SNU(1,J) - A
C           SNU(2,J) - Z
C           SNU(3,J) - Excitation energy (MeV)
C           SNU(4,J) - KINTER
C           SNU(5,J) ]
C           SNU(6,J) ]   -  Nucleus momentum (MeV/C)
C           SNU(7,J) ]
C-------------------------------------------------------------------------
      COMMON /SECHAR/ SPT(6,5000),SNU(10,101),LS6,LS100,LS10,LS11  ! GSITR
C
      DIMENSION PARTIC(24)
      CHARACTER*7 PARTIC
      DATA PARTIC /'NEUTRON',' PROTON','  PI-  ','  PI+  ','  PI0  ',
     *             ' ANTI N',' ANTI P','  K-   ','  K+   ','  K0   ',
     *             'ANTI K0',' GAMMA ','ELECTRN','POSITRN','  MU-  ',
     *             '  MU+  ','  NU E ','ANTINUE',' NU MU ','ANTINUM',
     *             'DEUTRON','TRITIUM','HELIUM3',' ALPHA '/
C
      CONTINUE!WRITE(36,10)NEVENT 
   10 FORMAT(/ 20X,'COMPLETED STAR No.',I3)   
      CONTINUE!WRITE(36,20)
   20 FORMAT(4X,'  Px(MeV/c)  Py(MeV/c)  Pz(MeV/c)    Ekin(MeV)', 
     *          ' PARTICLE')
C
C Variables for check of conservation laws
      ESUM=0.  ! Summary kinetic energy
      RMSUM=0. ! Summary mass of MESONS only!
      QSUM=0.  ! Summary charge Q
      BSUM=0.  ! Summary barion number
      SSUM=0.  ! Summary strangeness
      RLSUM=0. ! Summary lepton number
      PXSUM=0. ! Summary Px of all particles
      PYSUM=0. ! Summary Py of all particles
      PZSUM=0. ! Summary Pz of all particles
C
C ------------ Loop on all particles in SPT ---------------
      DO 100 J=1,LS100
        IPART=IFIX(SPT(5,J))          ! type of particle
      IF(IPART.EQ.0)GO TO 200         ! go out from the loop
        CT=SPT(1,J)
        ST=SQRT(1.0-CT*CT)
        SF=SPT(2,J)
        CF=SPT(3,J)
            E=SPT(4,J)
            W=SPT(6,J)
C
C         Conservation laws
        IF(IPART.eq. 1)THEN    ! Neutron
            ESUM=ESUM+E
            RMSUM=RMSUM+0.
            QSUM=QSUM+0.
            BSUM=BSUM+1.
            SSUM=SSUM+0.
            RLSUM=RLSUM+0.
              PMOD=SQRT(E*(E+2.0*940.))
              PX=PMOD*ST*CF
              PY=PMOD*ST*SF
              PZ=PMOD*CT
                PXSUM=PXSUM+PX
                PYSUM=PYSUM+PY
                PZSUM=PZSUM+PZ
        END IF
        IF(IPART.eq. 2)THEN    ! Proton
            ESUM=ESUM+E
            RMSUM=RMSUM+0.
            QSUM=QSUM+1.
            BSUM=BSUM+1.
            SSUM=SSUM+0.
            RLSUM=RLSUM+0.
              PMOD=SQRT(E*(E+2.0*940.))
              PX=PMOD*ST*CF
              PY=PMOD*ST*SF
              PZ=PMOD*CT
                PXSUM=PXSUM+PX
                PYSUM=PYSUM+PY
                PZSUM=PZSUM+PZ
        END IF
        IF(IPART.eq. 3)THEN    ! PI-
            ESUM=ESUM+E
            RMSUM=RMSUM+140.
            QSUM=QSUM-1.
            BSUM=BSUM+0.
            SSUM=SSUM+0.
            RLSUM=RLSUM+0.
              PMOD=SQRT(E*(E+2.0*140.))
              PX=PMOD*ST*CF
              PY=PMOD*ST*SF
              PZ=PMOD*CT
                PXSUM=PXSUM+PX
                PYSUM=PYSUM+PY
                PZSUM=PZSUM+PZ
        END IF
        IF(IPART.eq. 4)THEN    ! PI+
            ESUM=ESUM+E
            RMSUM=RMSUM+140.
            QSUM=QSUM+1.
            BSUM=BSUM+0.
            SSUM=SSUM+0.
            RLSUM=RLSUM+0.
              PMOD=SQRT(E*(E+2.0*140.))
              PX=PMOD*ST*CF
              PY=PMOD*ST*SF
              PZ=PMOD*CT
                PXSUM=PXSUM+PX
                PYSUM=PYSUM+PY
                PZSUM=PZSUM+PZ
        END IF
        IF(IPART.eq. 5)THEN    ! PI0
            ESUM=ESUM+E
            RMSUM=RMSUM+140.
            QSUM=QSUM+0.
            BSUM=BSUM+0.
            SSUM=SSUM+0.
            RLSUM=RLSUM+0.
              PMOD=SQRT(E*(E+2.0*140.))
              PX=PMOD*ST*CF
              PY=PMOD*ST*SF
              PZ=PMOD*CT
                PXSUM=PXSUM+PX
                PYSUM=PYSUM+PY
                PZSUM=PZSUM+PZ
        END IF
        IF(IPART.eq. 6)THEN    ! AntiN
            ESUM=ESUM+E
            RMSUM=RMSUM+0.
            QSUM=QSUM+0.
            BSUM=BSUM-1.
            SSUM=SSUM+0.
            RLSUM=RLSUM+0.
              PMOD=SQRT(E*(E+2.0*940.))
              PX=PMOD*ST*CF
              PY=PMOD*ST*SF
              PZ=PMOD*CT
                PXSUM=PXSUM+PX
                PYSUM=PYSUM+PY
                PZSUM=PZSUM+PZ
        END IF
        IF(IPART.eq. 7)THEN    ! AntiP
            ESUM=ESUM+E
            RMSUM=RMSUM+0.
            QSUM=QSUM-1.
            BSUM=BSUM-1.
            SSUM=SSUM+0.
            RLSUM=RLSUM+0.
              PMOD=SQRT(E*(E+2.0*940.))
              PX=PMOD*ST*CF
              PY=PMOD*ST*SF
              PZ=PMOD*CT
                PXSUM=PXSUM+PX
                PYSUM=PYSUM+PY
                PZSUM=PZSUM+PZ
        END IF
        IF(IPART.eq. 8)THEN    ! K-
            ESUM=ESUM+E
            RMSUM=RMSUM+495.
            QSUM=QSUM-1.
            BSUM=BSUM+0.
            SSUM=SSUM-1.
            RLSUM=RLSUM+0.
              PMOD=SQRT(E*(E+2.0*495.))
              PX=PMOD*ST*CF
              PY=PMOD*ST*SF
              PZ=PMOD*CT
                PXSUM=PXSUM+PX
                PYSUM=PYSUM+PY
                PZSUM=PZSUM+PZ
        END IF
        IF(IPART.eq. 9)THEN    ! K+
            ESUM=ESUM+E
            RMSUM=RMSUM+495.
            QSUM=QSUM+1.
            BSUM=BSUM+0.
            SSUM=SSUM+1.
            RLSUM=RLSUM+0.
              PMOD=SQRT(E*(E+2.0*495.))
              PX=PMOD*ST*CF
              PY=PMOD*ST*SF
              PZ=PMOD*CT
                PXSUM=PXSUM+PX
                PYSUM=PYSUM+PY
                PZSUM=PZSUM+PZ
        END IF
        IF(IPART.eq.10)THEN    ! K0
            ESUM=ESUM+E
            RMSUM=RMSUM+495.
            QSUM=QSUM+0.
            BSUM=BSUM+0.
            SSUM=SSUM+1.
            RLSUM=RLSUM+0.
              PMOD=SQRT(E*(E+2.0*495.))
              PX=PMOD*ST*CF
              PY=PMOD*ST*SF
              PZ=PMOD*CT
                PXSUM=PXSUM+PX
                PYSUM=PYSUM+PY
                PZSUM=PZSUM+PZ
        END IF
        IF(IPART.eq.11)THEN    ! AntiK0
            ESUM=ESUM+E
            RMSUM=RMSUM+0.
            QSUM=QSUM+0.
            BSUM=BSUM+0.
            SSUM=SSUM-1.
            RLSUM=RLSUM+0.
              PMOD=SQRT(E*(E+2.0*495.))
              PX=PMOD*ST*CF
              PY=PMOD*ST*SF
              PZ=PMOD*CT
                PXSUM=PXSUM+PX
                PYSUM=PYSUM+PY
                PZSUM=PZSUM+PZ
        END IF
        IF(IPART.eq.12)THEN    ! Gamma
            ESUM=ESUM+E
            RMSUM=RMSUM+0.
            QSUM=QSUM+0.
            BSUM=BSUM+0.
            SSUM=SSUM+0.
            RLSUM=RLSUM+0.
              PMOD=SQRT(E*(E+2.0*0.0))
              PX=PMOD*ST*CF
              PY=PMOD*ST*SF
              PZ=PMOD*CT
                PXSUM=PXSUM+PX
                PYSUM=PYSUM+PY
                PZSUM=PZSUM+PZ
        END IF
        IF(IPART.eq.21)THEN    ! Deuteron
            ESUM=ESUM+E
            RMSUM=RMSUM+0.
            QSUM=QSUM+1.
            BSUM=BSUM+2.
            SSUM=SSUM+0.
            RLSUM=RLSUM+0.
              PMOD=SQRT(E*(E+2.0*1880.))
              PX=PMOD*ST*CF
              PY=PMOD*ST*SF
              PZ=PMOD*CT
                PXSUM=PXSUM+PX
                PYSUM=PYSUM+PY
                PZSUM=PZSUM+PZ
        END IF
        IF(IPART.eq.22)THEN    ! Triton
            ESUM=ESUM+E
            RMSUM=RMSUM+0.
            QSUM=QSUM+1.
            BSUM=BSUM+3.
            SSUM=SSUM+0.
            RLSUM=RLSUM+0.
              PMOD=SQRT(E*(E+2.0*2820.))
              PX=PMOD*ST*CF
              PY=PMOD*ST*SF
              PZ=PMOD*CT
                PXSUM=PXSUM+PX
                PYSUM=PYSUM+PY
                PZSUM=PZSUM+PZ
        END IF
        IF(IPART.eq.23)THEN    ! He3
            ESUM=ESUM+E
            RMSUM=RMSUM+0.
            QSUM=QSUM+2.
            BSUM=BSUM+3.
            SSUM=SSUM+0.
            RLSUM=RLSUM+0.
              PMOD=SQRT(E*(E+2.0*2820.))
              PX=PMOD*ST*CF
              PY=PMOD*ST*SF
              PZ=PMOD*CT
                PXSUM=PXSUM+PX
                PYSUM=PYSUM+PY
                PZSUM=PZSUM+PZ
        END IF
        IF(IPART.eq.24)THEN    ! Alpha
            ESUM=ESUM+E
            RMSUM=RMSUM+0.
            QSUM=QSUM+2.
            BSUM=BSUM+4.
            SSUM=SSUM+0.
            RLSUM=RLSUM+0.
              PMOD=SQRT(E*(E+2.0*3760.))
              PX=PMOD*ST*CF
              PY=PMOD*ST*SF
              PZ=PMOD*CT
                PXSUM=PXSUM+PX
                PYSUM=PYSUM+PY
                PZSUM=PZSUM+PZ
        END IF
C
        CONTINUE!WRITE(36,30)J,PX,PY,PZ,E,PARTIC(IPART)
   30     FORMAT(I3,1X,3E11.4,E13.4,2X,A7)            
  100 CONTINUE
C -------------------- end of SPT loop -----------------
C
  200 CONTINUE
      CONTINUE!WRITE(36,39)
   39 FORMAT(66('-'))
      CONTINUE!WRITE(36,40)PXSUM,PYSUM,PZSUM,ESUM,RMSUM,
!     *            INT(QSUM),INT(BSUM),INT(SSUM),INT(RLSUM)
   40 FORMAT('Sum:',3E11.4,E13.4,'   Q   B   S   L'/
     *      30X,'MMESON=',E13.4,4I4)
C
C RESIDUAL NUCLEI CHARACTERISTICS
      CONTINUE!WRITE(36,50)
   50 FORMAT('     Residual nuclei:',34X,'A*   Z*    U*(MeV)')
C
      CONTINUE!WRITE(36,61)1,SNU(1,1),SNU(2,1)
   61 FORMAT(I3,1X,32X,'Nucleus-target  ',2F5.1)
C
      DO 80 J=2,LS11                       ! Nuclei-products
        KINTER=IFIX(SNU(4,J))
            IF(KINTER.EQ.0)RETURN
            IF(KINTER.EQ.6)THEN
            CONTINUE!WRITE(36,62)J,SNU(1,J),SNU(2,J)
   62           FORMAT(I3,1X,27X,'Full desintegration  ',2F5.1)
            GO TO 80      
            END IF
        EKIN=SQRT((SNU(5,J)**2+SNU(6,J)**2+SNU(7,J)**2)+
     *              (940.0*SNU(1,J))**2)-940.0*SNU(1,J)
        CONTINUE!WRITE(36,70)J,SNU(5,J),SNU(6,J),SNU(7,J),EKIN,
!     *                  SNU(1,J),SNU(2,J),SNU(3,J)
   70     FORMAT(I3,1X,3E11.4,E13.4,2X,2F5.1,E11.4)
   80 CONTINUE
C
      RETURN
      END



      SUBROUTINE STITLE
C Printout of MSDM stars
      COMMON /INREAC/ COST,SINF,COSF,TINT,WINT,JPART,NUCLID,KSTATE
        COMMON /HIPROJ/ APROJ,ZPROJ
      COMMON /ANTLAB/ lanti
      COMMON /INFEL/ ZNUC(110),ATWEI(110),DENS(110),RIEV(110),SYMB(110)
      CHARACTER*2 SYMB
C
      DIMENSION NAME(11)
      CHARACTER*8 NAME
      DATA NAME /' NEUTRON',' PROTON ','PI-MINUS',' PI-PLUS',' PI-ZERO',
     *'ANTINEUT','ANTIPROT','K-MINUS','K-PLUS','K0','ANTIK0'/
C
C----------------------- TITLE PRINTING -----------------------------
      CONTINUE!WRITE(36,1)
    1 FORMAT(' Printout of stars'/
     *       ' of the Hadron-Nucleus and Nucleus-Nucleus'/
     *       ' Generator of the SHIELD Code (MSDM-generator).')
C
C Is projectile Particle or Nucleus?
          IF(JPART.le.11)THEN                                        !HION2
      CONTINUE!WRITE(36,2) NAME(JPART),TINT,lanti
    2 FORMAT(5X,'INCIDENT PARTICLE IS ',A8/
     *       5X,'T0=',F11.1,' MEV.','  ANTILAB=',I1)
        ELSE
            IAPROJ=IFIX(APROJ)
            IZPROJ=IFIX(ZPROJ)
      CONTINUE!WRITE(36,52) SYMB(IZPROJ),IZPROJ,IAPROJ,TINT
   52 FORMAT(5X,'INCIDENT HEAVY ION IS ',A2,'(',I2,',',I3,').'/
     *       5X,'T0=',F11.1,' MEV/A.')
        END IF
C
C Nucleus-target:      
      CONTINUE!WRITE(36,3) SYMB(NUCLID),IFIX(ZNUC(NUCLID)),ATWEI(NUCLID)
    3 FORMAT(5X,'NUCLEUS-TARGET IS ',A2,'(',I2,',',F7.3,').')
C
      RETURN
      END
C---------------------- end of debug printout --------------------



      SUBROUTINE IMPACT(NEL)   ! HADES
      IMPLICIT REAL*8(A-H,O-Z)
C Debug attempt to build up the impact parameter of nucleus-nucleus
C interaction using coordinates of intranuclear nucleons.
      COMMON /CENTER/ XC(2,240),YC(2,240),ZC(2,240),IZ(2,240)
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX    !
      COMMON /HCASC/ ANUCL1,ANUCL2,ZNUCL1,ZNUCL2,T0,EPS1,EPS2,VPI,A1,A2,
     *               C1,C2,D1,D2,R0N1,R0N2,TF01,TF02,RM1,RM2
      COMMON /RIMPACT/ RIMPACT
C
c      CONTINUE!WRITE(25,10)NUMTRE,NEL
   10 FORMAT(/' NUMTRE=',I8,'     NEL=',I8)
C
c      CONTINUE!WRITE(25,20)
   20 FORMAT(' No.       X1       Y1       Z1')
      X1av=0.
      Y1av=0.
      Z1av=0.
      DO J=1,INT(ANUCL1)
        X1av=X1av+XC(1,J)
        Y1av=Y1av+YC(1,J)
        Z1av=Z1av+ZC(1,J)
c        CONTINUE!WRITE(25,100)J,XC(1,J),YC(1,J),ZC(1,J)
  100   FORMAT(I4,3F9.4)
      END DO
c      CONTINUE!WRITE(25,21)X1av,Y1av,Z1av
   21 FORMAT(' Rav',3F9.4)
C
c      CONTINUE!WRITE(25,30)
   30 FORMAT(' No.       X2       Y2       Z2')
      X2av=0.
      Y2av=0.
      Z2av=0.
      DO J=1,INT(ANUCL2)
        X2av=X2av+XC(2,J)
        Y2av=Y2av+YC(2,J)
        Z2av=Z2av+ZC(2,J)
c        CONTINUE!WRITE(25,100)J,XC(2,J),YC(2,J),ZC(2,J)
      END DO
c      CONTINUE!WRITE(25,21)X2av,Y2av,Z2av
C
      RIMPACT=SQRT((X1av-X2av)**2+(Y1av-Y2av)**2)
c      CONTINUE!WRITE(25,40)RIMPACT
   40 FORMAT(' Impact Parameter =',F9.3,' fm')
      RETURN
      END

