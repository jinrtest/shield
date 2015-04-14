      SUBROUTINE RANGE(INOVER)                                       ! GEMCA
C SAMPLING OF PARTICLE RANGE IN COMPLEX MEDIUM WITHOUT STRAGGLING AND
C MULTIPLE SCATTERING.
C     OUTPUT VARIABLE INOVER=1 IF INTERNAL ARRAY  STEPR  WAS EXCEEDED,
C     ELSE INOVER=0
C
C
C                          COMMONS OF GEOMETRY                       ! GEMCA
	COMMON/GDATAP/ X,Y,Z,CX,CY,CZ,XN,YN,ZN,CXN,CYN,CZN,          ! GEMCA
     *                 NZONO,MEDOLD,NBPO,NSO,                        ! GEMCA
     *                 NZONC,MEDCUR,NBPC,NSCI,NSCO,                  ! GEMCA
     *                 PINSFL,IOLEG                                  ! GEMCA
	INTEGER PINSFL                                               ! GEMCA
	REAL*8 X,Y,Z,CX,CY,CZ,XN,YN,ZN,CXN,CYN,CZN                   ! GEMCA
C                           COMMONS OF SHIELD/CG
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      COMMON /WRPR/ WP(15),IWP(7)                                    !HION1
      COMMON /TREE/ TREE(15,400000),ITREE(7,400000),LTR13,LTR7,LTR300 ! GSITR
      COMMON /MANALL/ NTREE,LINE1,LINE2,NBRMAX,NPRMAX,NNEUTR,NGAMMA,
     *                MELSC,MINTER,MCUT,MOUT,MABSOR,MDECAY,NSTRAN,
     *                NEMIN,NEPLS,NNUEL,NANUE,NNUMU,NANUM
      COMMON /STERAN/ STEPR(9,500),LST9,LST50
      COMMON /TMINPR/ TMINUS(25),TMINDF(25),LPR30                    !HION4
      COMMON /HIPROJ/ APROJ,ZPROJ                ! RDSHLD7, 03.01.01
      REAL*8 R,S                                                     ! GEMCA
C
C
      INOVER=0
C
C CLEANING OF THE RANGE STEPS ARRAY  STEPR.
      DO 1 J=1,LST50
      DO 1 K=1,LST9
    1 STEPR(K,J)=0.
C
C SENDING OF THE TRANSPORTING PARTICLE FROM LINE NO.LINE2 OF TREE TO WORK ARRAY
      DO 2 J=1,LTR13
    2 WP(J)=TREE(J,LINE2)
      DO 3 J=1,LTR7
    3 IWP(J)=ITREE(J,LINE2)
c      CONTINUE!WRITE(25,754)IWP(1),WP(7),WP(9)                           !HIONDEB
c  754 format('CALL RANGE: ',I2,F10.2,'  ZONE=',F6.1)            !HIONDEB
C
      JPART=IWP(1)
        IF(JPART.eq.25)THEN                                          !HION4
          APROJ=WP(14)                                               !HION4
          ZPROJ=WP(15)                                               !HION4
        END IF                                                       !HION4
	  T0=WP(7)
	      IF(TMINUS(JPART).GT.0.)THEN
		TMIN=TMINUS(JPART)
	      ELSE
		TMIN=TMINDF(JPART)
	      END IF
		  X0=WP(1)
		  Y0=WP(2)
		  Z0=WP(3)
		      COST=WP(4)
		      SINT=SQRT(1.0-COST**2)
		      SINF=WP(5)
		      COSF=WP(6)
        NEXTZ=INT(WP(9))                                   !HION4   ! GEMCA
	MEDCUR=GZMED(NEXTZ)                                         ! GEMCA
C
C For neutral kaons it is nesessary to determine: K0L or K0S?
      IF(JPART.EQ.10 .OR. JPART.EQ.11)THEN
	  K0L=IFIX(TREE(10,LINE2))
	  K0S=IFIX(TREE(11,LINE2))
	  IF(K0L.NE.0 .OR. K0S.NE.0)THEN
	      CONTINUE
	  ELSE
	      IF(RRAN(IX).LE.0.5)THEN
		  K0L=1
		  K0S=0
	      ELSE
		  K0L=0
		  K0S=1
	      END IF
	  END IF
      ELSE
	  K0L=0
	  K0S=0
      END IF
C
C DETERMINATION OF INPUT ZONE NUMBER IR AND MEDIUM NUMBER MED FOR POINT X0,Y0,Z0
	X=DBLE(X0)                                                   ! GEMCA
	Y=DBLE(Y0)                                                   ! GEMCA
	Z=DBLE(Z0)                                                   ! GEMCA
	CX=DBLE(SINT*COSF)                                           ! GEMCA
	CY=DBLE(SINT*SINF)                                           ! GEMCA
	CZ=DBLE(COST)                                                ! GEMCA
C
C
C FORMATION OF THE FIRST LINE OF RANGE STEPS ARRAY  STEPR
      STEPR(1,1)=0.
      STEPR(3,1)=FLOAT(NEXTZ)                                        ! GEMCA
      STEPR(5,1)=FLOAT(JPART)
      STEPR(6,1)=X0
      STEPR(7,1)=Y0
      STEPR(8,1)=Z0
      STEPR(9,1)=T0
C
C BRANCHING ON PARTICLE TYPE JPART
      IF(JPART.LT.1 .OR. JPART.GT.25)GO TO 1000                      !HION4
      GO TO (  20,  20,  20,  20,  50,  20,  20,  20,  20,  20,
     *         20,1000,1000,1000,  20,  20,1000,1000,1000,1000,
     *         20,  20,  20,  20,  20),JPART                         !HION4
C
C RANGE OF nucleons, charged pions, antinucleons, kaons and muons
   20 CONTINUE                                                       !HION4
      IF(JPART.eq.25)THEN                                            !HION8
        CALL CALMHI(0)                                               !HION4
        CALL ROPTHI(0)                                               !HION8
      END IF                                                         !HION8
      IF(T0.LE.TMIN)GO TO 900
 1320   DUMMY=RRAN(IX)
	IF(DUMMY.LE.0.0)GO TO 1320
      GAMMA=-ALOG(DUMMY)
C
C Sampling of decay lenght for K0S-meson
      IF(K0S.EQ.1)THEN
 1330   DUMMY=RRAN(IX)
	IF(DUMMY.LE.0.0)GO TO 1330
	  GAMMA1=-ALOG(DUMMY)
	  RDECAY=2.6769*(SQRT(T0*(T0+990.))/495.)*GAMMA1
	  RINT1=0.
	  RINT2=0.
      END IF
C
      ISTEP=0
	  TI=T0
	GO TO 27                                                     ! GEMCA
C
c DETERMINATION OF THE PASS LENGHT  S  INTO INPUT ZONE  NEXTZ        ! GEMCA
   26 CALL GNEXTZ(NEXTZ)                                             ! GEMCA
        IF(NEXTZ.EQ.0 .OR. MEDCUR.EQ.0)THEN                          ! GEMCA
C       Particle exit into external vacuum                           ! GEMCA
          KSTATE=3                                                   ! GEMCA
          GO TO 22                                                   ! GEMCA
        ELSE
C       Particle just have crossed some inrernal surface.
C*******Registrate X,Y,Z,TI=TF here if necessary! (sobol, 29.06.97)
          CONTINUE
	END IF                                                       ! GEMCA
C                                                                    ! GEMCA
   27 CALL GCURZL(NEXTZ,S)                                           ! GEMCA
c      CONTINUE!WRITE(25,755)NEXTZ,S                                     !HIONDEB
c  755 format('After GCURZL: NEXTZ=',I4,'  S=',F10.2)           !HIONDEB
C
      ISTEP=ISTEP+1
C
      IF(MEDCUR.EQ.1000)THEN                                         ! GEMCA
C     particle cross internal vacuum on step No.ISTEP                ! GEMCA
	  R=S                                                        ! GEMCA
	      TF=TI                                                  ! GEMCA
		  KSTATE=3                                           ! GEMCA
	  GOTO 21                                                    ! GEMCA
      END IF                                                         ! GEMCA
C
C DETERMINATION OF THE PARTICLE STATUS AS A RESULT OF PASSING THROUGH
C THE LAYER OF MEDIUM No.MEDCUR (NO.MED ! cg) DEPTH S.
      IF(JPART.EQ.1 .OR. JPART.EQ.6 .OR. JPART.EQ.10 .OR. JPART.EQ.11)
     *THEN
	  CALL PASSN(GAMMA,JPART,TI,MEDCUR,S,TF,R,KSTATE)            ! GEMCA
C         (as result we have: R, TF, KSTATE=0 or 3)
      ELSE
cC                                                                    !HIONDEB
c          CONTINUE!WRITE(25,7654)GAMMA,JPART,TI,TF,R,S,KSTATE                 !HIONDEB
c 7654     format(                                                    !HIONDEB
c     *    13X,'   GAMMA JPART        TI        TF         R         S',
c     *    ' KSTATE'/
c     *    'Befor PASSZ: ',F8.6,I6,4F10.2,I7)                         !HIONDEB
cC                                                                    !HIONDEB
          CALL PASSZ(GAMMA,JPART,TI,TMIN,MEDCUR,S,TF,R,KSTATE)       ! GEMCA
C         (as result we have: R, TF, KSTATE=0,3, OR 4)
cC                                                                    !HIONDEB
c          CONTINUE!WRITE(25,7655)GAMMA,JPART,TI,TF,R,S,KSTATE                 !HIONDEB
c 7655     format('After PASSZ: ',F8.6,I6,4F10.2,I7/)                 !HIONDEB
cC                                                                    !HIONDEB
	  IF(KSTATE.EQ.3)GO TO 19
	  IF(KSTATE.EQ.4 .AND. TMINUS(JPART).GT.0.)KSTATE=7
	  IF(KSTATE.EQ.4 .AND. JPART.EQ. 3)KSTATE=2
	  IF(KSTATE.EQ.4 .AND. JPART.EQ. 4)KSTATE=6
	  IF(KSTATE.EQ.4 .AND. JPART.EQ. 7)KSTATE=2
	  IF(KSTATE.EQ.4 .AND. JPART.EQ. 8)KSTATE=2
	  IF(KSTATE.EQ.4 .AND. JPART.EQ. 9)KSTATE=6
	  IF(KSTATE.EQ.4 .AND. JPART.EQ.15)KSTATE=6
	  IF(KSTATE.EQ.4 .AND. JPART.EQ.16)KSTATE=6
   19 CONTINUE
      END IF
C
C COORDINATES OF THE PASS END POINT.
   21 CONTINUE
C If K0S, then competition between decay and other possibilities
      IF(K0S.EQ.1)THEN
	  RINT2=RINT2+R
	  IF(RINT2.GE.RDECAY)THEN
C             decay
	      RPRIM=RDECAY-RINT1
	      RINT1=RINT1+RPRIM
	      R=RPRIM
	      KSTATE=5
	  ELSE
C             nondecay
	      RINT1=RINT1+R
	  END IF
      END IF
C
      IF((K0L.EQ.1 .OR. K0S.EQ.1) .AND. KSTATE.EQ.0)THEN
	  IF(RRAN(IX).LE.0.5)THEN
	      JPART=10
	  ELSE
	      JPART=11
	  END IF
      END IF
C
C move particle on distanse R in direction (cx,cy,cz).              ! GEMCA
	X=X+R*CX                                                    ! GEMCA
	Y=Y+R*CY                                                    ! GEMCA
	Z=Z+R*CZ                                                    ! GEMCA
C
C FORMATION OF THE LINE NO.ISTEP+1 OF ARRAY STEPR.
      IF(ISTEP+1.GT.LST50)THEN
	  CONTINUE!WRITE(25,1002)ISTEP,NUMTRE,IXINIT
 1002    FORMAT(' RANGE DIAGNOSTIC: OVERFILLING OF ARRAY STEPR, ISTEP=',
     *          I3/19X,'NUMTRE=',I6,'   IXINIT=',I12)
	  INOVER=1
	  ISTEP=ISTEP-1
	  GO TO 22
      END IF
C
	      ISTEP1=ISTEP+1
      STEPR(1,ISTEP1)=FLOAT(ISTEP)
      STEPR(2,ISTEP1)=FLOAT(MEDCUR)                                 ! GEMCA
      STEPR(3,ISTEP1)=FLOAT(NEXTZ)                                  ! GEMCA
      STEPR(4,ISTEP1)=SNGL(R)                                       ! GEMCA
      STEPR(5,ISTEP1)=FLOAT(JPART)
      STEPR(6,ISTEP1)=SNGL(X)                                       ! GEMCA
      STEPR(7,ISTEP1)=SNGL(Y)                                       ! GEMCA
      STEPR(8,ISTEP1)=SNGL(Z)                                       ! GEMCA
      STEPR(9,ISTEP1)=TF
C
      IF(KSTATE.EQ.3)GO TO 25                                       ! GEMCA
C
C THE STEP NO.ISTEP IS THE LAST STEP OF PARTICLE RANGE.
   22 STEPR(2,1)=FLOAT(ISTEP)
C
C RECORD OF PARTICLE IN FINAL STATE INTO LINE NO.LINE2+1 OF TREE
	      WP(1)=SNGL(X)                                         ! GEMCA
	      WP(2)=SNGL(Y)                                         ! GEMCA
	      WP(3)=SNGL(Z)                                         ! GEMCA
		  WP(7)=TF
	      WP(9)=NEXTZ                                           ! GEMCA
		      IWP(6)=1
		      IWP(7)=KSTATE
      DO 23 J=1,LTR13
   23 TREE(J,LINE2+1)=WP(J)
      DO 24 J=1,LTR7
   24 ITREE(J,LINE2+1)=IWP(J)
C
C If particle is neutral kaon and he fly out, then hi is K0L or K0S
C-Dementyev      IF((K0L.NE.0 .OR. K0S.NE.0) .AND. KSTATE.EQ.3)THEN
	  TREE(10,LINE2+1)=FLOAT(K0L)
	  TREE(11,LINE2+1)=FLOAT(K0S)
C-Dementyev      END IF
C
      call STETLEPA(WP(8))    ! TLEST
      RETURN
C                              MEDCUR
C PARTICLE CROSS THE MEDIUM NO.MED WITHOUT INTERACTION OR ABSORPTION
C (or decay?, sobol, 29.06.97) AND ENTRY INTO THE NEXT INPUT ZONE.
   25         TI=TF                                                ! GEMCA
      GO TO 26
C--------------- END OF NEUTRON, PROTON AND PION RENGE----------------------
C
C-------------------- PI0: DECAY IN THE BIRTH POINT-----------------------
   50 IF(T0.LE.TMIN)GO TO 900
	      R=0.
		  KSTATE=5
      ISTEP=1
C
C X,Y,Z,T final = X0,Y0,Z0,T0 respectivly.                          ! GEMCA
C FORMATION OF THE LINE NO.ISTEP+1 OF ARRAY STEPR.
	      ISTEP1=ISTEP+1
      STEPR(1,ISTEP1)=FLOAT(ISTEP)
      STEPR(2,ISTEP1)=FLOAT(MEDCUR)                                 ! GEMCA
      STEPR(3,ISTEP1)=FLOAT(NEXTZ)                                  ! GEMCA
      STEPR(4,ISTEP1)=SNGL(R)                                       ! GEMCA
      STEPR(5,ISTEP1)=FLOAT(JPART)
      STEPR(6,ISTEP1)=X0                                            ! GEMCA
      STEPR(7,ISTEP1)=Y0                                            ! GEMCA
      STEPR(8,ISTEP1)=Z0                                            ! GEMCA
      STEPR(9,ISTEP1)=T0                                            ! GEMCA
C
      STEPR(2,1)=FLOAT(ISTEP)
C
C RECORD OF PARTICLE IN FINAL STATE INTO LINE NO.LINE2+1 OF TREE
C WP(*) - do not change                                             ! GEMCA
		      IWP(6)=1
		      IWP(7)=KSTATE
      DO 53 J=1,LTR13
   53 TREE(J,LINE2+1)=WP(J)
      DO 54 J=1,LTR7
   54 ITREE(J,LINE2+1)=IWP(J)
      call STETLEPA(WP(8))    ! TLEST
      RETURN
C--------------------------END OF PI0 -------------------------------------
C
C---------------------- T0 < TMIN ------------------------------
  900         R=0.                                                  ! GEMCA
		  KSTATE=7
      ISTEP=1
C
C X,Y,Z,T final = X0,Y0,Z0,T0 respectivly.                          ! GEMCA
C FORMATION OF THE LINE NO.ISTEP+1 OF ARRAY STEPR.
	      ISTEP1=ISTEP+1
      STEPR(1,ISTEP1)=FLOAT(ISTEP)
      STEPR(2,ISTEP1)=FLOAT(MEDCUR)                                 ! GEMCA
      STEPR(3,ISTEP1)=FLOAT(NEXTZ)                                  ! GEMCA
      STEPR(4,ISTEP1)=SNGL(R)                                       ! GEMCA
      STEPR(5,ISTEP1)=FLOAT(JPART)
      STEPR(6,ISTEP1)=X0                                            ! GEMCA
      STEPR(7,ISTEP1)=Y0                                            ! GEMCA
      STEPR(8,ISTEP1)=Z0                                            ! GEMCA
      STEPR(9,ISTEP1)=T0                                            ! GEMCA
C
      STEPR(2,1)=FLOAT(ISTEP)
C
C RECORD OF PARTICLE IN FINAL STATE INTO LINE NO.LINE2+1 OF TREE
		      IWP(6)=1
		      IWP(7)=KSTATE
      DO 903 J=1,LTR13
  903 TREE(J,LINE2+1)=WP(J)
      DO 904 J=1,LTR7
  904 ITREE(J,LINE2+1)=IWP(J)
      call STETLEPA(WP(8))    ! TLEST
      RETURN
C ------------------------ END OF T0 < TMIN ----------------------
C
C-----------------------------DIAGNOSTIC------------------------------------
 1000 CONTINUE!WRITE(25,1001)JPART,NTREE
 1001 FORMAT(' RANGE DIAGNOSTIC: ABNORMAL JPART=',I3,'  IN NTREE=',I6)
      STOP
      END



      SUBROUTINE PASSZ(GAMMA,JPART,TI,TMIN,MEDIA,S,TF,R,KSTATE)
C PASS  R  UNTIL INTERACTION, ABSORPTION OR FLY OUT OF CHARGED PARTICLE TYPE
C JPART  WITH ENERGY  TI  IN LAYER DEPTH  S  OF MEDIUM NO.MEDIA.
C INPUT AND OUTPUT:
C                 GAMMA=-LOG(RAN)
C INPUT:
C      JPART - PARTICLE TYPE (2-4,7,8,9,15,16)
C      TI    - PARTICLE ENERGY (MEV)
C      TMIN  - PARTICLE CUT OFF ENERGY (MEV)
C      MEDIA - MEDIUM NUMBER (1-16)
C      S     - THICKNESS OF LAYER (CM)
C OUTPUT:
C       R      - PASS IN MEDIUM (CM)
C       TF     - ENERGY AS A RESULT OF PASS (MEV)
C       KSTATE - PARTICLE STATUS AS A RESULT OF PASS: 0 - INTERACTION
C                                                     3 - FLY OUT
C                                                     4 - ABSORPTION
C----------------------------------------------------------------------------
	REAL*8 R,S                                                  ! GEMCA
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      COMMON /RANMAC/ TR(100),RANMAC(100,13,48),OPTMAC(100,13,48),NUPARZ  ! CORNEL !HION
      COMMON /RECJPR/ LPART(25),INVERP(15),LPARTZ(25),INVERZ(13)      !HION
C FUNCTION-OPERATOR FOR LINEAR INTERPOLATION
      FLINT(X,X1,X2,Y1,Y2)=Y1+((Y2-Y1)/(X2-X1))*(X-X1)
C
      IPARTZ=LPARTZ(JPART)
	  IF(IPARTZ.LE.-1)THEN
	      CONTINUE!WRITE(25,1000)JPART
 1000         FORMAT(' PASSZ DIAGNOSTIC: INVALID JPART=',I3)
	      STOP
	  END IF
C
C For internal and external vacuum, Sobol, 25.03.2000
      IF(MEDIA.eq.0 .OR. MEDIA.eq.1000)THEN
        TF=TI
        R=S
        KSTATE=3
        RETURN
      END IF
C
C PASS  R  OF PARTICLE WITH ENERGY  TI  UNTIL ENERGY  TF.
      CALL SEARCH(TI,TR,100,N1,N2)
      RI=   FLINT(TI,TR(N1),TR(N2),RANMAC(N1,IPARTZ,MEDIA),
     *                             RANMAC(N2,IPARTZ,MEDIA))
	  NI1=N1
	  NI2=N2
      CALL SEARCH(TMIN,TR,100,N1,N2)
      RMIN= FLINT(TMIN,TR(N1),TR(N2),RANMAC(N1,IPARTZ,MEDIA),
     *                               RANMAC(N2,IPARTZ,MEDIA))
	  NF1=N1
	  NF2=N2
      R=RI-RMIN
C
      IF(R.GT.S)THEN
	  RIS=RI-S
          CALL SERCH3(RIS,RANMAC,100,13,16,IPARTZ,MEDIA,N1,N2)       !HION
	  TF=   FLINT(RIS,RANMAC(N1,IPARTZ,MEDIA),
     *                    RANMAC(N2,IPARTZ,MEDIA),TR(N1),TR(N2))
	      NF1=N1
	      NF2=N2
	  R=S
      ELSE
	  TF=TMIN
      END IF
C--------------------- R AND TF ARE FOUND ------------------------------
C
C OPTICAL DEPTH  OPTLEN  CORRESPONDING TO THE PASS  R.
      OPTLI=FLINT(TI,TR(NI1),TR(NI2),OPTMAC(NI1,IPARTZ,MEDIA),
     *                               OPTMAC(NI2,IPARTZ,MEDIA))

      OPTLF=FLINT(TF,TR(NF1),TR(NF2),OPTMAC(NF1,IPARTZ,MEDIA),
     *                               OPTMAC(NF2,IPARTZ,MEDIA))
	  OPTLEN=OPTLI-OPTLF
      IF(JPART.EQ.15 .OR. JPART.EQ.16)OPTLEN=0.
C
C
      IF(OPTLEN.GT.GAMMA)GO TO 1
	  IF(R.LT.S)THEN
C             PARTICLE ABSORPTION
	      KSTATE=4
	  ELSE
C             PARTICLE FLY OUT
	      KSTATE=3
		  GAMMA=GAMMA-OPTLEN
	  END IF
	  RETURN
C
C PARTICLE INTERACTION
    1 KSTATE=0
	  OPTLR=OPTLI-GAMMA
      CALL SERCH3(OPTLR,OPTMAC,100,13,16,IPARTZ,MEDIA,N1,N2)          !HION
      TF=FLINT(OPTLR,OPTMAC(N1,IPARTZ,MEDIA),OPTMAC(N2,IPARTZ,MEDIA),
     *                                               TR(N1),TR(N2))
      RF=FLINT(TF,TR(N1),TR(N2),RANMAC(N1,IPARTZ,MEDIA),
     *                          RANMAC(N2,IPARTZ,MEDIA))
	  R=RI-RF
      RETURN
	  END



      SUBROUTINE PASSN(GAMMA,JPART,TI,MEDIA,S,TF,R,KSTATE)
C PASS  R  UNTIL INTERACTION OR FLY OUT OF NEUTRON WITH ENERGY  TI
C IN LAYER DEPTH  S  OF MEDIUM NO.MEDIA.
C INPUT AND OUTPUT:
C                 GAMMA=-LOG(RAN)
C INPUT:
C      JPART = 1,6,10,11
C      TI    - NEUTRON ENERGY (MEV)
C      MEDIA - MEDIUM NUMBER (1-16)
C      S     - THICKNESS OF LAYER (CM)
C OUTPUT:
C       R      - PASS IN MEDIUM (CM)
C       TF=TI  - ENERGY AS A RESULT OF PASS (MEV)
C       KSTATE - NEUTRON STATUS AS A RESULT OF PASS: 0 - INTERACTION
C                                                    3 - FLU OUT
C----------------------------------------------------------------------------
	REAL*8 R,S                                                  ! GEMCA
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
C
C For internal and external vacuum, Sobol, 25.03.2000
      IF(MEDIA.eq.0 .OR. MEDIA.eq.1000)THEN
        TF=TI
        R=S
        KSTATE=3
        RETURN
      END IF
C
      IF(JPART.EQ.10 .OR. JPART.EQ.11)THEN
	  JP1=10
	  JP2=11
	  CALL MACROS(JP1,TI,MEDIA,SECM01)
	  CALL MACROS(JP2,TI,MEDIA,SECM02)
	  SECMAC=0.5*(SECM01+SECM02)
      ELSE
	  CALL MACROS(JPART,TI,MEDIA,SECMAC)
      END IF
C
      R=GAMMA/SECMAC
	  TF=TI
C
      IF(R.GT.S)THEN
	  KSTATE=3
	  R=S
	  GAMMA=GAMMA-SECMAC*S
      ELSE
	  KSTATE=0
      END IF
C
      RETURN
      END



C============================= TLEST ============================
      SUBROUTINE STETLEPA(WGHTLE)                         ! TLEST
C Analysing of the array STEPR(9,500) and building up
C of the Track Length Estimation for ALL CROSSED ZONES
C WGHTLE - Statistical weight of the particle
C
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      COMMON /STERAN/ STEPR(9,500),LST9,LST50
c      REAL*8 STEPR                     ! DPdEdX5
c      REAL*8 Ei,Ef,SLENG,dE,dR,Eaver   ! DPdEdX5
      COMMON /HIPROJ/ APROJ,ZPROJ
      COMMON /EGRIDL/ EGRIDL(119),NG70               ! KI05NG
C ----------------------- KI05#2 --------------------------------
      COMMON /RANMAC/ TR(100),RANMAC(100,13,48),OPTMAC(100,13,48),NUPARZ
c      COMMON /RANMAC/ TR(432),RANMAC(432,13,16),OPTMAC(432,13,16), ! DKFZ8
c     &                        DEDXRD(432,13,16),NUPARZ,LENGR       ! DKFZ8
c      REAL*8 TR,RANMAC,OPTMAC,DEDXRD   ! DPdEdX2                   ! DKFZ8
      DIMENSION RN(100)
c      DIMENSION RN(LENGR)
c      REAL*8 RN
      COMMON /RECJPR/ LPART(25),INVERP(15),LPARTZ(25),INVERZ(13)
C ----------------------- end KI05#2 ----------------------------
C
      LENGR=100
C
      JPART=INT(STEPR(5,1))
      IF(JPART.ge.5 .AND. JPART.le.20)RETURN  ! To reject nonused particles
C
      NSTEP=INT(STEPR(2,1)) ! Number of steps in the array STEPR
      DO J=2,NSTEP+1
        NZCUR=INT(STEPR(3,J))  ! Number of current zone
        MEDIUM=INT(STEPR(2,J)) ! Number of medium current zone   ! KI05#2
        Ei=STEPR(9,J-1)   ! Energy at the enter to zone NZCUR
        Ef=STEPR(9,J)     ! Energy at the exit from zone NZCUR (or absorption)
        SLENG=STEPR(4,J)  ! Pass lenghth in the zone NZCUR
        dE=Ei-Ef          ! Energy deposition in the zone NZCUR
C
C           HISTOGRAMMING FOR CHARHED PARTICLES AND FRAGMENTS
C Checking: if Ei and Ef both are inside the grid EGRIDL(NG70)
C It is necessary to have the correct input for SEARCH.
        IF(Ei .gt. EGRIDL(NG70) .OR. Ef .ge. EGRIDL(NG70))THEN   ! KI05NG
          GOTO 10  ! Energy is aside from energy grid
        END IF
C
C --------------------------- KI05#2 ----------------------------
      IPARTZ=LPARTZ(JPART)
      DO M=1,LENGR   ! Extract the Range to separate array RN
        RN(M)=RANMAC(M,IPARTZ,MEDIUM)
      END DO
C ----------------------- end KI05#2 ----------------------------
C Find the interval [N1Ef,N2Ei] of numbers of elements of the energy grid
C array EGRIDL(NG70) which includes both values Ei and Ef (Remember: Ei>Ef).
        CALL SEARCH(Ef,EGRIDL,NG70,N1,N2)   ! DPdEdX5   ! KI05NG
c        CALL SEARCH(SNGL(Ef),EGRIDL,NG70,N1,N2)   ! DPdEdX5   ! KI05NG
        N1Ef=N1
        N2Ef=N2
        CALL SEARCH(Ei,EGRIDL,NG70,N1,N2)   ! DPdEdX5   ! KI05NG
c        CALL SEARCH(SNGL(Ei),EGRIDL,NG70,N1,N2)   ! DPdEdX5   ! KI05NG
        N1Ei=N1
        N2Ei=N2
C
        IF(N2Ef.le.N1Ei)THEN
C---------Ei and Ef are in different bins of the energy grid EGRIDL(NG70)
          DO K=N1Ef+1,N2Ei
C
            IF(K.eq.N1Ef+1)THEN
*             dR=(SLENG/dE)*(EGRIDL(K)-Ef)        ! Linear interpolation
              dR=FQUINT4(EGRIDL(K),TR,RN,LENGR)-  ! Quadratic interpolation
     *           FQUINT4(Ef,       TR,RN,LENGR)
              Eaver=(EGRIDL(K)+Ef)*0.5
            END IF
C
            IF(K.gt.N1Ef+1 .AND. K.lt.N2Ei)THEN
*             dR=(SLENG/dE)*(EGRIDL(K)-EGRIDL(K-1))  ! Linear interpolation
              dR=FQUINT4(EGRIDL(K),  TR,RN,LENGR)-   ! Quadratic interpolation
     *           FQUINT4(EGRIDL(K-1),TR,RN,LENGR)
              Eaver=(EGRIDL(K)+EGRIDL(K-1))*0.5
            END IF

            IF(K.eq.N2Ei)THEN
*             dR=(SLENG/dE)*(Ei-EGRIDL(K-1))        ! Linear interpolation
              dR=FQUINT4(Ei,         TR,RN,LENGR)-  ! Quadratic interpolation
     *           FQUINT4(EGRIDL(K-1),TR,RN,LENGR)
              Eaver=(Ei+EGRIDL(K-1))*0.5
            END IF
C
            CALL HISTLEPA(JPART,DBLE(Eaver),(DBLE(WGHTLE)*dR),NZCUR)
          END DO
        ELSE
C---------Ei and Ef are in the same bin of the energy grid EGRIDL(NG70)
          DO K=N1Ef+1,N2Ei
            dR=SLENG
            Eaver=(Ei+Ef)*0.5
            CALL HISTLEPA(JPART,DBLE(Eaver),(DBLE(WGHTLE)*dR),NZCUR)
          END DO
        END IF
      END DO
C
   10 CONTINUE
      RETURN
      END



      SUBROUTINE HISTLEPA(JPART,ENER,WGHT,NZCUR)             ! TLEST
C Histogramming of the particles which cross given zone
C INPUT: JPART - Type of particle
C        ENER  - Energy MeV/A
C        WGHT - Statistical weight
C        NZCUR   - No. of zone which is crossed
C
      COMMON /HIPROJ/ APROJ,ZPROJ
      COMMON /BOOKED/ BOOKED(2,10)  ! (Zi,Ai)x10
      COMMON /TLEPA/ TLEPA(122,20,1000),NG73,LTLEPA
      REAL*8 TLEPA
      REAL*8 ENER,WGHT
C
C Fragments n, p, d, t, He3, He4, Li7, Be9, B11, C12 are registered always.
      IF(JPART.eq.1)THEN                                   ! Neutron
        CALL HIST3DPA(ENER,TLEPA, 1,NZCUR,WGHT)
      END IF
C
      IF(JPART.eq.2)THEN                                   ! proton
        CALL HIST3DPA(ENER,TLEPA, 2,NZCUR,WGHT)
      END IF
C
      IF(JPART.eq.21)THEN                                      ! H2
        CALL HIST3DPA(ENER,TLEPA, 3,NZCUR,WGHT)
      END IF
      IF(JPART.eq.25 .AND. APROJ.eq.2. .AND. ZPROJ.eq.1.)THEN  ! H2
        CALL HIST3DPA(ENER,TLEPA, 3,NZCUR,WGHT)
      END IF
C
      IF(JPART.eq.22)THEN                                      ! H3
        CALL HIST3DPA(ENER,TLEPA, 4,NZCUR,WGHT)
      END IF
      IF(JPART.eq.25 .AND. APROJ.eq.3. .AND. ZPROJ.eq.1.)THEN  ! H3
        CALL HIST3DPA(ENER,TLEPA, 4,NZCUR,WGHT)
      END IF
C
      IF(JPART.eq.23)THEN                                      ! He3
        CALL HIST3DPA(ENER,TLEPA, 5,NZCUR,WGHT)
      END IF
      IF(JPART.eq.25 .AND. APROJ.eq.3. .AND. ZPROJ.eq.2.)THEN  ! He3
        CALL HIST3DPA(ENER,TLEPA, 5,NZCUR,WGHT)
      END IF
C
      IF(JPART.eq.24)THEN                                      ! He4
        CALL HIST3DPA(ENER,TLEPA, 6,NZCUR,WGHT)
      END IF
      IF(JPART.eq.25 .AND. APROJ.eq.4. .AND. ZPROJ.eq.2.)THEN  ! He4
        CALL HIST3DPA(ENER,TLEPA, 6,NZCUR,WGHT)
      END IF
C
      IF(JPART.eq.25 .AND. APROJ.eq.7. .AND. ZPROJ.eq.3.)THEN  ! Li7
        CALL HIST3DPA(ENER,TLEPA, 7,NZCUR,WGHT)
      END IF
C
      IF(JPART.eq.25 .AND. APROJ.eq.9. .AND. ZPROJ.eq.4.)THEN  ! Be9
        CALL HIST3DPA(ENER,TLEPA, 8,NZCUR,WGHT)
      END IF
C
      IF(JPART.eq.25 .AND. APROJ.eq.11. .AND. ZPROJ.eq.5.)THEN ! B11
        CALL HIST3DPA(ENER,TLEPA, 9,NZCUR,WGHT)
      END IF
C
      IF(JPART.eq.25 .AND. APROJ.eq.12. .AND. ZPROJ.eq.6.)THEN ! C12
        CALL HIST3DPA(ENER,TLEPA,10,NZCUR,WGHT)
      END IF
C
C The following 10 (ten) fragments are registered according to booking of user
C in the array /BOOKED/ BOOKED(2,10): BOOKED(1,i)=Zi, BOOKED(2,i)=Ai. 
      DO J=1,10
        IF(JPART.eq.25)THEN                  ! Ion
          IF(BOOKED(2,J) .eq. 1000.)THEN
C              Fixed Z and any A
            IF(ZPROJ.eq.BOOKED(1,J))
     *      CALL HIST3DPA(ENER,TLEPA,10+J,NZCUR,WGHT)
          ELSE
C           Fixed (Z,A)
            IF(APROJ.eq.BOOKED(2,J) .AND. ZPROJ.eq.BOOKED(1,J))
     *      CALL HIST3DPA(ENER,TLEPA,10+J,NZCUR,WGHT)
          END IF
        END IF
C
        IF(JPART.eq.3 .OR. JPART.eq.4)THEN   ! Pi+/-
          IF(BOOKED(1,J).eq.-3. .AND. BOOKED(2,J).eq.-3.)
     *    CALL HIST3DPA(ENER,TLEPA,10+J,NZCUR,WGHT)
        END IF
      END DO
C
      RETURN
      END



      SUBROUTINE HIST3DPA(E,SPEC,NP,NZ,W)                  ! TLEST
      IMPLICIT REAL*8 (A-H,O-Z)
C Histogramming of REAL*8 variable E>0 (energy, MeV/u) 
C with statistical weight W according to fixed logatithmic
C grid EGRIDL(NG70) in the limits 0<E<1000 
C (NG73 is equal NG73-4 bins + 4 comulatives).
C in 3-dimensional REAL*8 array SPEC(NG73,20,1000):
C                  NP - No. of fragment (1-20)
C                  NZ - No. of zone (1-1000)
C
        COMMON /EGRIDL/ EGRIDL(119),NG70   ! REAL*4
        REAL*4 EGRIDL
      DIMENSION SPEC(NG70+3,20,1000)
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
C
C Check of input:
      IF(E.le.0.)THEN
        CONTINUE!WRITE(25,1000)E,NUMTRE
 1000   FORMAT(' HIST3DPA DIAGNOSTIC: E<0, E=',E10.3,',  NUMTRE=',I8)
        RETURN
      END IF
C Registration for all energies E>0.
        SPEC(NG70+2,NP,NZ)=SPEC(NG70+2,NP,NZ)+W
        SPEC(NG70+3,NP,NZ)=SPEC(NG70+3,NP,NZ)+W*E
C Registration for energies in the interval 0<E<1000000(1TeV).
      IF(E.lt.1000000.0)THEN
        SPEC(NG70,  NP,NZ)=SPEC(NG70,  NP,NZ)+W
        SPEC(NG70+1,NP,NZ)=SPEC(NG70+1,NP,NZ)+W*E
          CALL SEARCH(SNGL(E),EGRIDL,NG70,N1,N2)   ! REAL*4
        SPEC(N1,NP,NZ)=SPEC(N1,NP,NZ)+W
      END IF
C
      RETURN
      END



      FUNCTION FQUINT4(X,RX,RF,N)      ! TLEST
C FUNCTION FOR QUADRATIC INTERPOLATION
C INPUT: X argument
C        RX(1:N) - array of argument values
C        RF(1:N) - array of function values
C        WARNING!!! N.ge.3
C OUTPUT: FQUINT4
C------------------------------------------
      DIMENSION RX(N),RF(N)
C
      IF(X.le.RX(1))THEN
        N1=1
        N2=2
        N3=3
        GOTO 1
      END IF
C
      IF(X.ge.RX(N-1))THEN   ! It must be N2<N. If N2=N, then N3=N+1!
        N1=N-2                                         ! out of array!
        N2=N-1
        N3=N
        GOTO 1
      END IF
C
      CALL SEARCH(X,RX,N,N1,N2)
        N3=N2+1              ! It must be N2<N. If N2=N, then N3=N+1!
C                                                        out of array!
    1 CONTINUE
      X1=RX(N1)
      X2=RX(N2)
      X3=RX(N3)
      Y1=RF(N1)
      Y2=RF(N2)
      Y3=RF(N3)
C
      FQUINT4=
     &     (((X1-X)*(X2-X))/((X1-X3)*(X2-X3)))*Y3+
     &     (((X2-X)*(X3-X))/((X2-X1)*(X3-X1)))*Y1+
     &     (((X3-X)*(X1-X))/((X3-X2)*(X1-X2)))*Y2
C
      RETURN
      END
