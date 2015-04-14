
***********************************************************************
*                                                                     *
*        The following subroutines are written by N.M.Sobolevsky,     *
*        Institute for Nuclear Research RAS, Moscow,                  *
*        E-mail: sobolevs@AL20.inr.troitsk.ru                         *
*                                                                     *
***********************************************************************

      SUBROUTINE MICROD(JPART,TMEV,NUCLID,SIGTOT,SIGIN)
C       
C  CALCULATION OF PARTICLE-NUCLEUS TOTAL AND INELASTIC CROSS SECTIONS
C  INPUT:
C       JPART - INCIDENT PARTICLE TYPE.   JPART=1,2,3,...,24,25.    !HION
C       TMEV -  INCIDENT PARTICLE/ION KINETIC ENERGY (MEV/A). TMEV>0.
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*   ATTENTION! Incident energy TMEV is the total kinetic energy of the  *
*   projectile in laboratory system in MeV (not in MeV/A!)  19.08.2005  *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C       NUCLID - NUCLEUS-TARGET NUMBER. NUCLID=1,2,....107.
C  OUTPUT:
C       SIGTOT - TOTALL CROSS SECTION (MB)
C       SIGIN - INELASTIC CROSS SECTION (MB)
C INTRINSIC VARIABLES:
C       NREC(110),NUC12(2,110),INVERN(20),LPART(25),INVERP(15),     !HION
C       LPARTZ(25),INVERZ(13) -                                     !HION
C     - ARRAYS FOR RECODING VARIABLES  NUCLID,JPART  INTO  INUCL,IPART
C       AND VICE VERSA (SEE BLOCK DATA CSDATA)       
C ----------------------------------------------------------------------------
C
      COMMON /INFEL/ ZNUC(110),ATWEI(110),DENS(110),RIEV(110),SYMB(110)
      CHARACTER*2 SYMB
      COMMON /RECOD/ NREC(110),NUC12(2,110),INVERN(20)
      COMMON /RECJPR/ LPART(25),INVERP(15),LPARTZ(25),INVERZ(13)     !HION
      COMMON /HIPROJ/ APROJ,ZPROJ       ! Is used if JPART=25        !HION
      COMMON /LEVERS/ PARLEV(40)  ! DKFZf
      REAL*8 TMEVSDP,TLIMITDP,RELATDP
C
C CHECK OF INPUT VARIABLES VALIDITY
C--------------------------------NUCLID
      IF(NUCLID.LT.1 .OR. NUCLID.GT.107)THEN                  !HION
          WRITE(6,1000)NUCLID                                !HION
          CONTINUE!WRITE(25,1000)NUCLID                                !HION
 1000     FORMAT(' MICROD DIAGNOSTIC: INVALID NUCLID=',I6)    !HION
	  STOP 'MICROD1'
      END IF
C--------------------------------TMEV
      IF(TMEV.LT.0.)THEN
	  WRITE(6,1001)TMEV
	  CONTINUE!WRITE(25,1001)TMEV
 1001     FORMAT(' MICROD DIAGNOSTIC: INVALID TMEV=',F12.4)
	  STOP 'MICROD2'
      END IF
C--------------------------------IPART
      IPART=LPART(JPART)
      IF(IPART.LE.-1)THEN
	  WRITE(6,1002)JPART
	  CONTINUE!WRITE(25,1002)JPART
 1002     FORMAT(' MICROD DIAGNOSTIC: PARTICLE TYPE JPART=',I2,
     *           ' DON''T ESTABLISHED')
	  STOP 'MICROD3'
      END IF
C
          PERMUT1=0.   ! DKFZ06
C
C******************** HEAVY ION, 13.07.97 ******************* !HION
      IF(JPART.ge.21 .and. JPART.le.25)THEN
          IF(JPART.eq.21)THEN   ! Deuteron
            A1=2.
            Z1=1.
          END IF
C
          IF(JPART.eq.22)THEN   ! Triton
            A1=3.
            Z1=1.
          END IF
C
          IF(JPART.eq.23)THEN   ! Helium3
            A1=3.
            Z1=2.
          END IF
C
          IF(JPART.eq.24)THEN   ! Alpha (He4)
            A1=4.
            Z1=2.
          END IF
C
          IF(JPART.eq.25)THEN   ! Heavy ion
            A1=APROJ
            Z1=ZPROJ
          END IF
C
          A2=ATWEI(NUCLID)      ! Nucleus-target
          Z2=ZNUC(NUCLID)
C ++++++++++++++++++++++++++++++++ DKFZ06
C Ion-projectile is proton
          if(A1.eq.1.0 .AND. Z1.eq.1.0)then
            IPART=2
            goto 10
          end if
C +++++++++++++++++++++++++++++end DKFZ06
C
C If nucleus impinges hidrogen, let us make PERMUTATION Proj <---> Targ
          PERMUT1=0.
          IF(NUCLID.eq.1)THEN   ! PERMUTATION
            PERMUT1=1.           ! Flag for permutation
C:::::::::::: 19.08.2005 ::::::::::::::::
            TSAV=TMEV     ! Saving of input energy TMEV of projectile (A1,Z1)
            TMEV=TMEV/A1  ! Kinetic energy (MeV) of target H1 in antilab frame
C:::::::::::: end 19.08.2005 ::::::::::::
            JPRSAV=JPART            ! To save input value.
            NCLSAV=NUCLID           ! To save input value.
            IPRSAV=IPART            ! To save input value.
              JPART=2   ! The target H1 temporary serves as proton-projectile
              IPART=LPART(JPART)   ! 19.08.2005, this is the same as IPART=2
C###################### 2004a ###################
              NUCLID=INT(Z1)    ! The projectile temporary serves as a target
              IF(A1.eq.  2.0 .AND. Z1.eq. 1.0)NUCLID=101   ! Deuteron
              IF(A1.eq.  3.0 .AND. Z1.eq. 1.0)NUCLID=102   ! Tritium
              IF(A1.eq.235.0 .AND. Z1.eq.92.0)NUCLID=103   ! U235
              IF(A1.eq.  3.0 .AND. Z1.eq. 2.0)NUCLID=104   ! He3
              IF(A1.eq.  6.0 .AND. Z1.eq. 3.0)NUCLID=105   ! Li6
              IF(A1.eq. 10.0 .AND. Z1.eq. 5.0)NUCLID=106   ! B10
              IF(A1.eq.240.0 .AND. Z1.eq.94.0)NUCLID=107   ! Pu240
C###################### end 2004a ###############
              GOTO 10   ! Go to proton-nucleus interaction
          END IF   ! End of PERMUTATION Proj <---> Targ
C
C Nucleus-nucleus interaction
          TMEVS=TMEV         ! Saved value (MeV/A)   ! DKFZc
C                            ! WARNING! TMEV is in MeV, 19.08.2005
          TLIMIT=MIN(A1,A2)  ! PARAMETER: Limit energy for exponential smoothing
          IF(TMEV.lt.TLIMIT)TMEV=TLIMIT               ! DKFZc
C
C Special particular case: Projectile is Deuteron	and Target is Deuteron
          IF(A2.eq.2.0 .AND. Z2.eq.1.0 .AND.  ! Projectile is Deuteron ! DKFZc
     &       A1.eq.2.0 .AND. Z1.eq.1.0)THEN   ! and Target is Deuteron (1.0,2.0)
C           Calulation for H2-target via He4-target by means
C           of A**2/3-extrapolation: SIG(H2)=SIG(He4)*[(2/4)**2/3]
            SIGTOT=SIGION(1,A1,Z1,4.0,2.0,TMEV)*0.63
            SIGIN =SIGION(2,A1,Z1,4.0,2.0,TMEV)*0.63
C             Correction of error from SIGION for interaction   ! DKFZc
C             of Deuteron with lightest nuclei
              IF(A1.eq.2.0 .AND. Z1.eq.1.0 .AND. SIGIN.ge.SIGTOT)THEN
                SIGIN=SIGTOT
              END IF
            IF(TMEVS.lt.TLIMIT)THEN   ! DKFZc
C             T<TLIMIT - exponential smoothing to zero
              SIGTOT=SIGTOT*((EXP(TMEVS)-1.0)/(EXP(TLIMIT)-1.0))
              SIGIN =SIGIN *((EXP(TMEVS)-1.0)/(EXP(TLIMIT)-1.0))
            END IF
            TMEV=TMEVS   ! DKFZc
            SIGTOT=SIGTOT*PARLEV(39)   ! DKFZf
            SIGIN =SIGIN *PARLEV(39)   ! DKFZf
            RETURN
          END IF
C
          SIGTOT=SIGION(1,A1,Z1,A2,Z2,TMEV)
          SIGIN =SIGION(2,A1,Z1,A2,Z2,TMEV)
C
C         Correction of error from SIGION for interaction of any nuclei
          IF(SIGIN.ge.SIGTOT)SIGIN=SIGTOT   ! 19.08.2005
C
          IF(TMEVS.lt.TLIMIT)THEN   ! DKFZc
C           T<TLIMIT - exponential smoothing to zero
C           DP to avoid the error M6104: Floating point overflow (26.06.2005)
            TMEVSDP=DBLE(TMEVS)
            TLIMITDP=DBLE(TLIMIT)
            RELATDP=(EXP(TMEVSDP)-1.0D0)/(EXP(TLIMITDP)-1.0D0)
            SIGTOT=SIGTOT*SNGL(RELATDP)
            SIGIN =SIGIN *SNGL(RELATDP)
          END IF
          TMEV=TMEVS   ! DKFZc
          SIGTOT=SIGTOT*PARLEV(39)   ! DKFZf
          SIGIN =SIGIN *PARLEV(39)   ! DKFZf
          RETURN
      END IF
C******************** end !HION *********************
C
C Hadron-nucleus interaction
   10 CONTINUE
C--------------------------------INUCL
      INUCL=NREC(NUCLID)
      IF(INUCL.LE.-1)THEN
	  GOTO 200  ! Nucleus is not the reference nucleus
      ELSE
	  GOTO 100  ! Nucleus is the reference nucleus
      END IF
C
C NUCLEUS  NUCLID  IS REFERENCE NUCLEUS
  100 CONTINUE
      CALL MICLIN(IPART,TMEV,INUCL,SIGTOT,SIGIN)
C
C******************** HEAVY ION, 13.07.97 ******************* !HION
C !!! Correction of 01.03.2001: Inverse permutation after label 100
C                               was omitted !
C Inverse permutation:
      IF(PERMUT1.gt.0.5)THEN
          TMEV=TSAV           ! To restore input value 19.08.2005
          JPART=JPRSAV        ! To restore input value
          NUCLID=NCLSAV       ! To restore input value
          IPART=IPRSAV        ! To restore input value
      END IF
C******************** end HION *******************
      SIGTOT=SIGTOT*PARLEV(40)   ! DKFZf
      SIGIN =SIGIN *PARLEV(40)   ! DKFZf
      RETURN
C
C NUCLEUS  NUCLID  IS NOT REFERENCE NUCLEUS.
C INTERPOLATION ON A**2/3 IS REQUIRED.
  200 CONTINUE
C
C--------------------------FIRST REFERENCE NUCLEUS
      INUCL1=NUC12(1,NUCLID)       ! The first reference nucleus for NUCLID
	  NUCLI1=INVERN(INUCL1)      ! The first reference NUCLID1 for NUCLID
	      ATW1=ATWEI(NUCLI1)     ! Atomic weight of NUCLID1
	      A1P23=ATW1**0.6666667  ! Atomic weight of NUCLID1 to the 2/3 power
		    A=ATWEI(NUCLID)      ! Atomic weight of NUCLID
		    AP23=A**0.6666667    ! Atomic weight of NUCLID to the 2/3 power
      CALL MICLIN(IPART,TMEV,INUCL1,SGTOT1,SGIN1)
C
      IF(INUCL1.EQ.20)THEN
C         EXTRAPOLATION ON A**2/3 FOR NUCLEI HEAVIER THEN U238
        SIGTOT=(SGTOT1/A1P23)*AP23
        SIGIN =(SGIN1 /A1P23)*AP23
        SIGTOT=SIGTOT*PARLEV(40)   ! DKFZf
        SIGIN =SIGIN *PARLEV(40)   ! DKFZf
        RETURN
      END IF
C
C--------------------------SECOND REFERENCE NUCLEUS
      INUCL2=NUC12(2,NUCLID)       ! The second reference nucleus for NUCLID
	  NUCLI2=INVERN(INUCL2)      ! The second reference NUCLID2 for NUCLID
	      ATW2=ATWEI(NUCLI2)     ! Atomic weight of NUCLID2
	      A2P23=ATW2**0.6666667  ! Atomic weight of NUCLID2 to the 2/3 power
		    RATIO=(AP23-A1P23)/(A2P23-A1P23)   ! 0<RATIO<1
      CALL MICLIN(IPART,TMEV,INUCL2,SGTOT2,SGIN2)
C
cC-----------Quadratic INTERPOLATION ON X=A**2/3 across two points (OLD)
c	  WTOT1=SGTOT1/A1P23
c	  WIN1 =SGIN1 /A1P23
c      SIGTOT=(WTOT1+(SGTOT2/A2P23-WTOT1)*RATIO)*AP23
c      SIGIN =(WIN1 +( SGIN2/A2P23-WIN1 )*RATIO)*AP23
C
C-----------Linear INTERPOLATION ON X=A**2/3 across two points, DKFZc
      SIGTOT=SGTOT1+(SGTOT2-SGTOT1)*RATIO
      SIGIN = SGIN1+( SGIN2- SGIN1)*RATIO
C
C******************** HEAVY ION, 13.07.97 ******************* !HION
C Inverse permutation:
      IF(PERMUT1.gt.0.5)THEN
          TMEV=TSAV           ! To restore input value 19.08.2005
          JPART=JPRSAV        ! To restore input value
          NUCLID=NCLSAV       ! To restore input value
          IPART=IPRSAV        ! To restore input value
      END IF
C******************** end HION *******************
C
      SIGTOT=SIGTOT*PARLEV(40)   ! DKFZf
      SIGIN =SIGIN *PARLEV(40)   ! DKFZf
      RETURN
      END

 

      SUBROUTINE MICLIN(IPART,TMEV,INUCL,SGTOT,SGIN)   ! DKFZc
C CALCULATION OF TOTAL AND INELASTIC CROSS SECTIONS BY Quadratic 
C interpolation IN ENERGY (with exponential smoothing to zero at
C low energies) IN ARRAY CROSEC(I,J,K,L) FOR REFERENCE NUCLEI
C  INPUT:
C       IPART - PARTICLE NUMBER IPART=LPART(JPART)
C       TMEV -  PARTICLE KINETIC ENERGY (MEV)
C       INUCL - REFERENCE NUCLEUS NUMBER  (1 - 20)
C  OUTPUT:
C       SGTOT - TOTALL CROSS SECTION (MB)
C       SGIN - INELASTIC CROSS SECTION (MB)
C  CROSS SECTIONS ARRAY:
C        CROSEC(I,J,K,L):  I=1-51, ENERGY MESH
C                              J=1 - ENERGY (MEV)
C                              J=2 - TOTALL CROSS SECTION (MB)
C                              J=3 - INELASTIC CROSS SECTION (MB)
C                                  K=1-20, REFERENCE NUCLEI
C                                      L=1-10 - PARTICLE NUMBER IPART
C----------------------------------------------------------------------------
C
      COMMON /CRSSEC/ CROSEC(51,3,20,10),NUPART
      COMMON /INFEL/ ZNUC(110),ATWEI(110),DENS(110),RIEV(110),SYMB(110)
      CHARACTER*2 SYMB
      COMMON /RECOD/ NREC(110),NUC12(2,110),INVERN(20)
C
      NE=IFIX(CROSEC(51,1,INUCL,IPART))  ! Number of points in energy grid
C
      IF(TMEV.le.CROSEC(   1,1,INUCL,IPART))GOTO 100  ! TMEV<Tmin=T1
C
      IF(TMEV.ge.CROSEC(NE-1,1,INUCL,IPART))THEN  ! TMEV>Tmax-1
        N1=NE-2
        N2=NE-1
        N3=NE
        GOTO 10
      END IF
C
      CALL SERCH4(TMEV,NE,INUCL,IPART,N1,N2)   ! TMEV within [Tmin,Tmax-1]
      N3=N2+1          ! This is correct because TMEV within [Tmin,Tmax-1] 
      GOTO 10
C
   10 CONTINUE
      X=TMEV
C     Energy
      X1=CROSEC(N1,1,INUCL,IPART)
      X2=CROSEC(N2,1,INUCL,IPART)
      X3=CROSEC(N3,1,INUCL,IPART)
C     Total cross section
      Y1=CROSEC(N1,2,INUCL,IPART)
      Y2=CROSEC(N2,2,INUCL,IPART)
      Y3=CROSEC(N3,2,INUCL,IPART)
C     Inelastic cross section
      Z1=CROSEC(N1,3,INUCL,IPART)
      Z2=CROSEC(N2,3,INUCL,IPART)
      Z3=CROSEC(N3,3,INUCL,IPART)
C
      SGTOT=
     &     (((X1-X)*(X2-X))/((X1-X3)*(X2-X3)))*Y3+
     &     (((X2-X)*(X3-X))/((X2-X1)*(X3-X1)))*Y1+
     &     (((X3-X)*(X1-X))/((X3-X2)*(X1-X2)))*Y2
      SGIN=
     &     (((X1-X)*(X2-X))/((X1-X3)*(X2-X3)))*Z3+
     &     (((X2-X)*(X3-X))/((X2-X1)*(X3-X1)))*Z1+
     &     (((X3-X)*(X1-X))/((X3-X2)*(X1-X2)))*Z2
      SGTOT=MAX(0.,SGTOT)  ! To exclude negative cross section     ! DKFZc
      SGIN=MAX(0.,SGIN)    ! at the interpolation near a threshold ! DKFZc
      RETURN
C
C T<Tmin - exponential smoothing to zero
  100 CONTINUE
      T1   =CROSEC(1,1,INUCL,IPART)  ! Minimal energy node Tmin of the grid 
      S1TOT=CROSEC(1,2,INUCL,IPART)  ! SIGTOT(Tmin)
      S1IN =CROSEC(1,3,INUCL,IPART)  ! SIGIN(Tmin)
C
      NUCLID=INVERN(INUCL) ! Reference nucleus
      AREF=ATWEI(NUCLID)   ! Mass of the Reference nucleus
      P=1.0/AREF           ! PARAMETER of EXP smoothing   ! DKFZc 
C
      if(T1.le.0.0)then   ! 04.09.2005
        SGTOT=S1TOT
        SGIN =S1IN
      else
        SGTOT=S1TOT*((EXP(P*TMEV)-1.0)/(EXP(P*T1)-1.0))
        SGIN =S1IN *((EXP(P*TMEV)-1.0)/(EXP(P*T1)-1.0))
      end if
C
      RETURN
      END
 
 
 
      SUBROUTINE SERCH4(E,NE,INUCL,IPART,N1,N2)
C  BINARY SEARCH IN ARRANGED (IN INCREASING ORDER) REFERENSE ENERGY ROW
C  OF 4-DIMENSIONAL CROSS SECTION ARRAY  CROSEC(N,1,INUCL,IPART)
C  INPUT:
C        E - CURRENT ENERGY
C        NE - NUMBER OF ENERGY POINTS FOR A GIVEN   INUCL,IPART
C        INUCL - REFERENCE NUCLEUS (1-20)
C        IPART - PARTICLE NUMBER IPART=LPART(JPART)
C  OUTPUT:
C   N1, N2 SUCH, THAT CROSEC(N1,1,INUCL,IPART) < E < CROSEC(N2,1,INUCL,IPART)
C----------------------------------------------------------------------------
C
      COMMON /CRSSEC/ CROSEC(51,3,20,10),NUPART
C
      N1=1
      N2=NE
    1 NM=(N1+N2)/2
      IF(E.LE.CROSEC(NM,1,INUCL,IPART))THEN
	  N2=NM
      ELSE
	  N1=NM
      END IF
      IF(N2.EQ.N1+1)RETURN
      GOTO 1
      END



      BLOCK DATA CSDATA
C PARTICLE-NUCLEUS CROSS SECTIONS DATA TABLES
C
C*****************************************************************************
      COMMON /RECOD/ NREC(110),NUC12(2,110),INVERN(20)
      COMMON /RECJPR/ LPART(25),INVERP(15),LPARTZ(25),INVERZ(13)     !HION
      COMMON /CRSSEC/ CROSEC(51,3,20,10),NUPART
C*****************************************************************************
C
C RECODING FROM NUCLID TO REFERENCE NUCLEUS NUMBER INUCL=NREC(NUCLID)
C FOR CROSS SECTION INTERPOLATION ON A**2/3
C
C NUCLID NUCLEUS NREC      NUCLID NUCLEUS NREC      NUCLID NUCLEUS NREC
C      1    H1      1          33            -          70            -
C    101    H2      2          34            -          71            -
C    102    H3      -          35            -          72            -
C    104   He3      -          36            -          73            -
C      2   HE4      3          37            -          74    W      18
C    105   Li6      -          38            -          75            -
C      3   LI7      4          39            -          76            -
C      4   BE       5          40            -          77            -
C    106    B10     -          41            -          78            -
C      5    B11     -          42   MO      15          79            -
C      6    C       6          43            -          80            -
C      7    N       7          44            -          81            -
C      8    O       8          45            -          82   PB      19
C      9            -          46            -          83            -
C     10            -          47            -          84            -
C     11   NA       9          48   CD      16          85            -
C     12            -          49            -          86            -
C     13   AL      10          50   SN      17          87            -
C     14   SI      11          51            -          88            -
C     15            -          52            -          89            -
C     16            -          53            -          90            -
C     17            -          54            -          91            -
C     18            -          55            -         103    U235    -
C     19            -          56            -          92    U238   20
C     20   CA      12          57            -          93            -
C     21            -          58            -          94   Pu239    -
C     22            -          59            -         107   Pu240    -
C     23            -          60            -          95            -
C     24            -          61            -          96            -
C     25            -          62            -          97            -
C     26   FE      13          63            -          98            -
C     27            -          64            -          99            -
C     28            -          65            -         100            -
C     29   CU      14          66            -         108            -
C     30            -          67            -         109            -
C     31            -          68            -         110            -
C     32            -          69            -
C........................................................................
C
C In fact there are only 18 reference nuclei (not 20), because the data
C for nuclei-targets H(1,2) and Li(3,7) are omitted in the array CROSEC.
C This is the reason for which NREC(3)=-4 and NREC(101)=-2.
C All arrays: NREC(110), NUC12(2,110) and INVERN(20) are correct!
C                                                     DKFZc, 09.05.2005
C..................1....2....3....4....5....6....7....8....9...10.........NUCLID
C                                                                        : 
      DATA NREC/   1,   3,  -4,   5,  -1,   6,   7,   8,  -1,  -1,       : 10
     *             9,  -1,  10,  11,  -1,  -1,  -1,  -1,  -1,  12,       : 20
     *            -1,  -1,  -1,  -1,  -1,  13,  -1,  -1,  14,  -1,       : 30
     *            -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,       : 40
     *            -1,  15,  -1,  -1,  -1,  -1,  -1,  16,  -1,  17,       : 50
     *            -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,       : 60
     *            -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,       : 70
     *            -1,  -1,  -1,  18,  -1,  -1,  -1,  -1,  -1,  -1,       : 80
     *            -1,  19,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,       : 90
     *            -1,  20,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,       :100
     *            -2,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1/       :110
C                                                                        :
C........................................................................:
C
C..................1..2.....1..2.....1..2.....1..2.....1..2...............NUCLID
C                                                                        :
      DATA NUC12/  1, 1,    3, 3,    3, 5,    5, 5,    5, 6,             :  5
     *             6, 6,    7, 7,    8, 8,    8, 9,    8, 9,             : 10
     *             9, 9,    9,10,   10,10,   11,11,   11,12,             : 15
     *            11,12,   11,12,   11,12,   11,12,   12,12,             : 20
     *            12,13,   12,13,   12,13,   12,13,   12,13,             : 25
     *            13,13,   13,14,   13,14,   14,14,   14,15,             : 30
     *            14,15,   14,15,   14,15,   14,15,   14,15,             : 35
     *            14,15,   14,15,   14,15,   14,15,   14,15,             : 40
     *            14,15,   15,15,   15,16,   15,16,   15,16,             : 45
     *            15,16,   15,16,   16,16,   16,17,   17,17,             : 50
     *            17,18,   17,18,   17,18,   17,18,   17,18,             : 55
     *            17,18,   17,18,   17,18,   17,18,   17,18,             : 60
     *            17,18,   17,18,   17,18,   17,18,   17,18,             : 65
     *            17,18,   17,18,   17,18,   17,18,   17,18,             : 70
     *            17,18,   17,18,   17,18,   18,18,   18,19,             : 75
     *            18,19,   18,19,   18,19,   18,19,   18,19,             : 80
     *            18,19,   19,19,   19,20,   19,20,   19,20,             : 85
     *            19,20,   19,20,   19,20,   19,20,   19,20,             : 90
     *            19,20,   20,20,   20,20,   20,20,   20,20,             : 95
     *            20,20,   20,20,   20,20,   20,20,   20,20,             :100
     *             1, 3,    1, 3,   19,20,    1, 3,    3, 5,             :105
     *             5, 6,   20,20,   20,20,   20,20,   20,20/             :110
C                                                                        :
C.........................................................................
C
      DATA INVERN /  1, 101,   2,   3,   4,   6,   7,   8,  11,  13,
     *              14,  20,  26,  29,  42,  48,  50,  74,  82,  92/
C
C*************************************************************************
C      COMMON /RECJPR/ LPART(25),INVERP(15),LPARTZ(25),INVERZ(13)     !HION
C*************************************************************************
C
C RECODING FROM PARTICLE TYPES  JPART  TO THE CURRENT PARTICLE NUMBER
C LPART(JPART)  FOR ARRAYS  CROSEC,SIGMAC  AND TO THE CHARGE PARTICLE 
C NUMBER LPARTZ(JPART) FOR ARRAYS  RANMAC,OPTMAC.
C
C JPART   PARTICLE   LPART LPARTZ      JPART   PARTICLE   LPART LPARTZ
C     1   NEUTRON        1      -         13   ELECTRON       -      -
C     2   PROTON         2      1         14   POSITRON       -      -
C     3   PI-            3      2         15   MUON-          -      7
C     4   PI+            4      3         16   MUON+          -      8
C     5   PI0            -      -         17   NU E           -      -
C     6   ANTI N         5      -         18   ANTI NU E      -      -
C     7   ANTI P         6      4         19   NU MU          -      -
C     8   K-             7      5         20   ANTI NU MU     -      -
C     9   K+             8      6         21   DEUTERON      11      9  !HION
C    10   K0             9      -         22   TRITON        12     10  !HION
C    11   K0 TILDA      10      -         23   HE3           13     11  !HION
C    12   GAMMA          -      -         24   ALPHA         14     12  !HION
C                                         25   HEAVY ION     15     13  !HION
C............................................................
C
      DATA LPART / 1,  2,  3,  4, -1,  5,  6,  7,  8,  9, 10, -1,       !HION
     *            -1, -1, -1, -1, -1, -1, -1, -1, 11, 12, 13, 14, 15/   !HION
C
      DATA INVERP / 1, 2, 3, 4, 6, 7, 8, 9, 10, 11, 21, 22, 23, 24, 25/ !HION
C
      DATA LPARTZ /-1,  1,  2,  3, -1, -1,  4,  5,  6, -1, -1, -1,      !HION
     *             -1, -1,  7,  8, -1, -1, -1, -1,  9, 10, 11, 12, 13/  !HION
C
      DATA INVERZ / 2, 3, 4, 7, 8, 9, 15, 16, 21, 22, 23, 24, 25/       !HION
C
C******************************************************************************
C      COMMON /CRSSEC/ CROSEC(51,3,20,10),NUPART
C******************************************************************************
C  CROSS SECTION ARRAY CROSEC(I,J,K,L):
C          I=1-51, ENERGY MESH
C              J=1 - ENERGY (MEV)
C              J=2 - TOTAL CROSS SECTION (MB)
C              J=3 - INELASTIC CROSS SECTION (MB)
C                  K=1-20, REFERENCE NUCLEI
C                      L=1-10, PARTICLE NUMBER IPART
C ----------------------------------------------------------------------------
C
C *********************** DATA FOR H(1,1), INUCL=1 (FROM INTRANUCLEAR CASCADE)
C                                   NEUTRON
      DATA (CROSEC(I,1,1,1),I=1,51)/
     *      0.0,      1.0,      3.0,      5.0,      7.0,
     *     10.0,     15.0,     20.0,     25.0,     30.0,
     *     40.0,     50.0,     70.0,    100.0,    150.0,
     *    200.0,    250.0,    300.0,    350.0,    400.0,
     *    450.0,    500.0,    600.0,    700.0,    800.0,
     *    900.0,   1000.0,   1100.0,   1300.0,   1500.0,
     *   2000.0,   3000.0,   4000.0,   5000.0,   7000.0,
     *  10000.0,  20000.0,  50000.0, 100000.0, 500000.0,
     *1000000.0,  9*0.,41./  
      DATA (CROSEC(I,2,1,1),I=1,51)/
     *20357.0,17564.1,12546.5, 8286.4, 4783.7,  950.0,  678.8,
     *  480.0,  380.0,  300.0,  200.0,  160.0,  108.0,   74.0,
     *   50.0,   41.0,   36.5,   34.0,   32.5,   32.0,   33.2,
     *   34.2,   35.5,   36.6,   37.4,   38.1,   38.6,   39.0,
     *   40.0,   40.0,   40.5,   41.1,   42.3,   42.3,   42.0,
     *   41.1,   39.6,   38.2,   39.0,   39.0,   39.0 , 10*0./
      DATA (CROSEC(I,3,1,1),I=1,51)/
     *    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,
     *    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,
     *    0.0,    0.0,    0.0,    0.0,    0.0,    0.8,    1.9,
     *    3.4,    8.6,   13.5,   17.3,   19.6,   21.2,   22.5,
     *   24.0,   24.8,   26.5,   29.0,   31.5,   31.6,   31.0,
     *   31.5,   32.3,   33.8,   33.0,   33.0,   33.0,  10*0./
C                                   PROTON
      DATA (CROSEC(I,1,1,2),I=1,51)/
     *      0.0,      1.0,      3.0,      5.0,      7.0,
     *     10.0,     15.0,     20.0,     25.0,     30.0,
     *     40.0,     50.0,     70.0,    100.0,    150.0,
     *    200.0,    250.0,    300.0,    350.0,    400.0,
     *    450.0,    500.0,    600.0,    700.0,    800.0,
     *    900.0,   1000.0,   1100.0,   1300.0,   1500.0,
     *   2000.0,   3000.0,   4000.0,   5000.0,   7000.0,
     *  10000.0,  20000.0,  50000.0, 100000.0, 500000.0,
     *1000000.0,  9*0.,41./
      DATA (CROSEC(I,2,1,2),I=1,51)/
     *17613.0,15114.9,10631.9, 6833.1, 3718.7,  330.0,  227.3,
     *  154.0,  121.0,   96.0,   70.0,   51.0,   38.4,   30.0,
     *   23.6,   22.4,   22.2,   22.6,   23.4,   24.7,   26.8,
     *   29.5,   37.3,   43.8,   47.8,   47.9,   47.3,   47.0,
     *   46.7,   46.0,   45.0,   43.0,   41.2,   40.8,   40.5,
     *   39.8,   39.0,   39.0,   39.0,   39.0,   39.0,  10*0./
      DATA (CROSEC(I,3,1,2),I=1,51)/
     *    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,
     *    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,
     *    0.0,    0.0,    0.0,    0.0,    0.6,    1.8,    2.7,
     *    4.5,   12.3,   18.8,   22.8,   22.9,   22.0,   22.0,
     *   24.7,   26.5,   27.5,   28.0,   27.5,   28.8,   29.5,
     *   30.0,   30.3,   34.4,   32.0,   32.0,   32.0,  10*0./
C                                   PI-
      DATA (CROSEC(I,1,1,3),I=1,51)/
     *      0.0,     20.0,     40.0,     60.0,     80.0,
     *     90.0,    100.0,    120.0,    135.0,    150.0,
     *    165.0,    175.0,    185.0,    195.0,    210.0,
     *    225.0,    250.0,    300.0,    350.0,    400.0,
     *    450.0,    500.0,    550.0,    580.0,    600.0,
     *    620.0,    650.0,    680.0,    700.0,    750.0,
     *    800.0,    840.0,    880.0,    920.0,    960.0,
     *   1000.0,   1100.0,   1200.0,   1400.0,   2000.0,
     *   3000.0,  10000.0,  20000.0, 100000.0,1000000.0,
     *5*0.,45./
      DATA (CROSEC(I,2,1,3),I=1,51)/
     *    6.3,    6.8,    8.0,   11.6,   16.3,   19.6,   23.8,
     *   39.2,   49.0,   57.5,   65.9,   68.5,   67.9,   65.9,
     *   59.7,   52.5,   40.5,   30.8,   25.7,   26.9,   29.0,
     *   32.1,   41.5,   44.6,   45.6,   43.0,   40.1,   38.4,
     *   38.0,   39.0,   44.3,   52.5,   57.2,   55.2,   50.1,
     *   45.7,   38.9,   35.3,   34.3,   33.7,   32.4,   26.5,
     *   25.0,   25.0,   25.0,  6*0./
      DATA (CROSEC(I,3,1,3),I=1,51)/
     *    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,
     *    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,
     *    0.0,    0.0,    0.0,    0.7,    2.1,    4.6,    6.8,
     *    8.6,   13.9,   16.0,   17.0,   16.8,   16.7,   16.8,
     *   17.0,   18.1,   19.5,   20.7,   23.5,   24.3,   23.9,
     *   23.5,   21.9,   21.6,   22.9,   24.1,   25.4,   21.5,
     *   21.0,   21.0,   21.0,  6*0./
C                                   PI+
      DATA (CROSEC(I,1,1,4),I=1,51)/
     *      0.0,     20.0,     30.0,     40.0,     50.0,
     *     60.0,     80.0,     90.0,    100.0,    120.0,
     *    140.0,    160.0,    170.0,    175.0,    182.0,
     *    190.0,    195.0,    200.0,    225.0,    250.0,
     *    300.0,    350.0,    400.0,    450.0,    500.0,
     *    600.0,    650.0,    700.0,    800.0,    900.0,
     *   1000.0,   1100.0,   1200.0,   1300.0,   1400.0,
     *   1500.0,   1600.0,   1800.0,   2000.0,   2500.0,
     *   3000.0,   4000.0,  10000.0,  20000.0, 100000.0,
     *1000000.0,  4*0.,46./
      DATA (CROSEC(I,2,1,4),I=1,51)/
     *    1.9,    3.5,    5.5,    9.0,   14.0,   17.4,   33.4,
     *   45.6,   60.0,  107.4,  146.6,  180.8,  192.0,  195.0,
     *  195.2,  192.6,  189.4,  185.0,  145.0,  113.0,   73.0,
     *   45.0,   33.0,   25.2,   21.6,   15.6,   14.8,   15.2,
     *   19.5,   24.5,   27.6,   32.2,   36.7,   41.0,   39.0,
     *   35.5,   32.3,   31.7,   31.1,   29.9,   28.9,   27.7,
     *   24.9,   23.5,   23.5,   23.5,  5*0./
      DATA (CROSEC(I,3,1,4),I=1,51)/
     *    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,
     *    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,
     *    0.0,    0.1,    0.1,    0.1,    0.2,    0.2,    0.3,
     *    0.6,    1.2,    2.0,    3.0,    4.8,    6.1,    7.5,
     *   10.5,   13.2,   14.1,   17.2,   19.8,   22.0,   22.1,
     *   20.8,   19.7,   20.6,   21.3,   22.6,   23.2,   22.1,
     *   20.0,   19.5,   19.5,   19.5,  5*0./
C                                   ANTINEUTRON
      DATA (CROSEC(I,1,1,5),I=1,51)/
     *      0.0,      1.0,      3.0,      5.0,      7.0,
     *     10.0,     15.0,     20.0,     25.0,     30.0,
     *     40.0,     50.0,     70.0,    100.0,    150.0,
     *    200.0,    250.0,    300.0,    350.0,    400.0,
     *    450.0,    500.0,    600.0,    700.0,    800.0,
     *    900.0,   1000.0,   1100.0,   1300.0,   1500.0,
     *   2000.0,   3000.0,   4000.0,   5000.0,   7000.0,
     *  10000.0,  20000.0,  50000.0, 100000.0, 500000.0,
     *1000000.0,  9*0.,41./  
      DATA (CROSEC(I,2,1,5),I=1,51)/
     *20357.0,17564.1,12546.5, 8286.4, 4783.7,  950.0,  678.8,
     *  480.0,  380.0,  300.0,  200.0,  160.0,  108.0,   74.0,
     *   50.0,   41.0,   36.5,   34.0,   32.5,   32.0,   33.2,
     *   34.2,   35.5,   36.6,   37.4,   38.1,   38.6,   39.0,
     *   40.0,   40.0,   40.5,   41.1,   42.3,   42.3,   42.0,
     *   41.1,   39.6,   38.2,   39.0,   39.0,   39.0 , 10*0./
      DATA (CROSEC(I,3,1,5),I=1,51)/
     *    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,
     *    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,
     *    0.0,    0.0,    0.0,    0.0,    0.0,    0.8,    1.9,
     *    3.4,    8.6,   13.5,   17.3,   19.6,   21.2,   22.5,
     *   24.0,   24.8,   26.5,   29.0,   31.5,   31.6,   31.0,
     *   31.5,   32.3,   33.8,   33.0,   33.0,   33.0,  10*0./
C                                   ANTIPROTON
      DATA (CROSEC(I,1,1,6),I=1,51)/
     *      0.0,      1.0,      3.0,      5.0,      7.0,
     *     10.0,     15.0,     20.0,     25.0,     30.0,
     *     40.0,     50.0,     70.0,    100.0,    150.0,
     *    200.0,    250.0,    300.0,    350.0,    400.0,
     *    450.0,    500.0,    600.0,    700.0,    800.0,
     *    900.0,   1000.0,   1100.0,   1300.0,   1500.0,
     *   2000.0,   3000.0,   4000.0,   5000.0,   7000.0,
     *  10000.0,  20000.0,  50000.0, 100000.0, 500000.0,
     *1000000.0,  9*0.,41./
      DATA (CROSEC(I,2,1,6),I=1,51)/
     *17613.0,15114.9,10631.9, 6833.1, 3718.7,  330.0,  227.3,
     *  154.0,  121.0,   96.0,   70.0,   51.0,   38.4,   30.0,
     *   23.6,   22.4,   22.2,   22.6,   23.4,   24.7,   26.8,
     *   29.5,   37.3,   43.8,   47.8,   47.9,   47.3,   47.0,
     *   46.7,   46.0,   45.0,   43.0,   41.2,   40.8,   40.5,
     *   39.8,   39.0,   39.0,   39.0,   39.0,   39.0,  10*0./
      DATA (CROSEC(I,3,1,6),I=1,51)/
     *    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,
     *    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,
     *    0.0,    0.0,    0.0,    0.0,    0.6,    1.8,    2.7,
     *    4.5,   12.3,   18.8,   22.8,   22.9,   22.0,   22.0,
     *   24.7,   26.5,   27.5,   28.0,   27.5,   28.8,   29.5,
     *   30.0,   30.3,   34.4,   32.0,   32.0,   32.0,  10*0./
C                                   K-
      DATA (CROSEC(I,1,1,7),I=1,51)/
     *      0.0,     20.0,     40.0,     60.0,     80.0,
     *     90.0,    100.0,    120.0,    135.0,    150.0,
     *    165.0,    175.0,    185.0,    195.0,    210.0,
     *    225.0,    250.0,    300.0,    350.0,    400.0,
     *    450.0,    500.0,    550.0,    580.0,    600.0,
     *    620.0,    650.0,    680.0,    700.0,    750.0,
     *    800.0,    840.0,    880.0,    920.0,    960.0,
     *   1000.0,   1100.0,   1200.0,   1400.0,   2000.0,
     *   3000.0,  10000.0,  20000.0, 100000.0,1000000.0,
     *5*0.,45./
      DATA (CROSEC(I,2,1,7),I=1,51)/
     *    6.3,    6.8,    8.0,   11.6,   16.3,   19.6,   23.8,
     *   39.2,   49.0,   57.5,   65.9,   68.5,   67.9,   65.9,
     *   59.7,   52.5,   40.5,   30.8,   25.7,   26.9,   29.0,
     *   32.1,   41.5,   44.6,   45.6,   43.0,   40.1,   38.4,
     *   38.0,   39.0,   44.3,   52.5,   57.2,   55.2,   50.1,
     *   45.7,   38.9,   35.3,   34.3,   33.7,   32.4,   26.5,
     *   25.0,   25.0,   25.0,  6*0./
      DATA (CROSEC(I,3,1,7),I=1,51)/
     *    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,
     *    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,
     *    0.0,    0.0,    0.0,    0.7,    2.1,    4.6,    6.8,
     *    8.6,   13.9,   16.0,   17.0,   16.8,   16.7,   16.8,
     *   17.0,   18.1,   19.5,   20.7,   23.5,   24.3,   23.9,
     *   23.5,   21.9,   21.6,   22.9,   24.1,   25.4,   21.5,
     *   21.0,   21.0,   21.0,  6*0./
C                                   K+
      DATA (CROSEC(I,1,1,8),I=1,51)/
     *      0.0,     20.0,     30.0,     40.0,     50.0,
     *     60.0,     80.0,     90.0,    100.0,    120.0,
     *    140.0,    160.0,    170.0,    175.0,    182.0,
     *    190.0,    195.0,    200.0,    225.0,    250.0,
     *    300.0,    350.0,    400.0,    450.0,    500.0,
     *    600.0,    650.0,    700.0,    800.0,    900.0,
     *   1000.0,   1100.0,   1200.0,   1300.0,   1400.0,
     *   1500.0,   1600.0,   1800.0,   2000.0,   2500.0,
     *   3000.0,   4000.0,  10000.0,  20000.0, 100000.0,
     *1000000.0,  4*0.,46./
      DATA (CROSEC(I,2,1,8),I=1,51)/
     *    1.9,    3.5,    5.5,    9.0,   14.0,   17.4,   33.4,
     *   45.6,   60.0,  107.4,  146.6,  180.8,  192.0,  195.0,
     *  195.2,  192.6,  189.4,  185.0,  145.0,  113.0,   73.0,
     *   45.0,   33.0,   25.2,   21.6,   15.6,   14.8,   15.2,
     *   19.5,   24.5,   27.6,   32.2,   36.7,   41.0,   39.0,
     *   35.5,   32.3,   31.7,   31.1,   29.9,   28.9,   27.7,
     *   24.9,   23.5,   23.5,   23.5,  5*0./
      DATA (CROSEC(I,3,1,8),I=1,51)/
     *    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,
     *    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,
     *    0.0,    0.1,    0.1,    0.1,    0.2,    0.2,    0.3,
     *    0.6,    1.2,    2.0,    3.0,    4.8,    6.1,    7.5,
     *   10.5,   13.2,   14.1,   17.2,   19.8,   22.0,   22.1,
     *   20.8,   19.7,   20.6,   21.3,   22.6,   23.2,   22.1,
     *   20.0,   19.5,   19.5,   19.5,  5*0./
C
C                                   K0
      DATA (CROSEC(I,1,1,9),I=1,51)/
     *      0.0,     20.0,     30.0,     40.0,     50.0,
     *     60.0,     80.0,     90.0,    100.0,    120.0,
     *    140.0,    160.0,    170.0,    175.0,    182.0,
     *    190.0,    195.0,    200.0,    225.0,    250.0,
     *    300.0,    350.0,    400.0,    450.0,    500.0,
     *    600.0,    650.0,    700.0,    800.0,    900.0,
     *   1000.0,   1100.0,   1200.0,   1300.0,   1400.0,
     *   1500.0,   1600.0,   1800.0,   2000.0,   2500.0,
     *   3000.0,   4000.0,  10000.0,  20000.0, 100000.0,
     *1000000.0,  4*0.,46./
      DATA (CROSEC(I,2,1,9),I=1,51)/
     *    1.9,    3.5,    5.5,    9.0,   14.0,   17.4,   33.4,
     *   45.6,   60.0,  107.4,  146.6,  180.8,  192.0,  195.0,
     *  195.2,  192.6,  189.4,  185.0,  145.0,  113.0,   73.0,
     *   45.0,   33.0,   25.2,   21.6,   15.6,   14.8,   15.2,
     *   19.5,   24.5,   27.6,   32.2,   36.7,   41.0,   39.0,
     *   35.5,   32.3,   31.7,   31.1,   29.9,   28.9,   27.7,
     *   24.9,   23.5,   23.5,   23.5,  5*0./
      DATA (CROSEC(I,3,1,9),I=1,51)/
     *    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,
     *    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,
     *    0.0,    0.1,    0.1,    0.1,    0.2,    0.2,    0.3,
     *    0.6,    1.2,    2.0,    3.0,    4.8,    6.1,    7.5,
     *   10.5,   13.2,   14.1,   17.2,   19.8,   22.0,   22.1,
     *   20.8,   19.7,   20.6,   21.3,   22.6,   23.2,   22.1,
     *   20.0,   19.5,   19.5,   19.5,  5*0./
C                            ANTI K0
      DATA (CROSEC(I,1,1,10),I=1,51)/
     *      0.0,     20.0,     40.0,     60.0,     80.0,
     *     90.0,    100.0,    120.0,    135.0,    150.0,
     *    165.0,    175.0,    185.0,    195.0,    210.0,
     *    225.0,    250.0,    300.0,    350.0,    400.0,
     *    450.0,    500.0,    550.0,    580.0,    600.0,
     *    620.0,    650.0,    680.0,    700.0,    750.0,
     *    800.0,    840.0,    880.0,    920.0,    960.0,
     *   1000.0,   1100.0,   1200.0,   1400.0,   2000.0,
     *   3000.0,  10000.0,  20000.0, 100000.0,1000000.0,
     *5*0.,45./
      DATA (CROSEC(I,2,1,10),I=1,51)/
     *    6.3,    6.8,    8.0,   11.6,   16.3,   19.6,   23.8,
     *   39.2,   49.0,   57.5,   65.9,   68.5,   67.9,   65.9,
     *   59.7,   52.5,   40.5,   30.8,   25.7,   26.9,   29.0,
     *   32.1,   41.5,   44.6,   45.6,   43.0,   40.1,   38.4,
     *   38.0,   39.0,   44.3,   52.5,   57.2,   55.2,   50.1,
     *   45.7,   38.9,   35.3,   34.3,   33.7,   32.4,   26.5,
     *   25.0,   25.0,   25.0,  6*0./
      DATA (CROSEC(I,3,1,10),I=1,51)/
     *    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,
     *    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,
     *    0.0,    0.0,    0.0,    0.7,    2.1,    4.6,    6.8,
     *    8.6,   13.9,   16.0,   17.0,   16.8,   16.7,   16.8,
     *   17.0,   18.1,   19.5,   20.7,   23.5,   24.3,   23.9,
     *   23.5,   21.9,   21.6,   22.9,   24.1,   25.4,   21.5,
     *   21.0,   21.0,   21.0,  6*0./
C
C*********************  DATA FOR HE(2,4), INUCL=3 *******************
C                                   NEUTRON
      DATA (CROSEC(I,1,3,1),I=1,51)/ 0.,
     *14.,15.,17.,20.,22.,25.,27.,30.,35.,40.,45.,50.,60.,70.,80.,90.,
     *100.,120.,140.,150.,160.,180.,200.,250.,300.,350.,
     *400.,500.,600.,700.,800.,900.,1000.,1500.,2000.,3000.,
     *5000.,7000.,10000.,20000.,50000.,100000.,500000.,1000000.,
     *5*0.,45./
      DATA (CROSEC(I,2,3,1),I=1,51)/ 2070.,
     * 1090.,1020.,915.,800.,710.,640.,600.,560.,500.,
     *  440., 390.,360.,295.,256.,220.,192.,168.,136.,120.,
     *  116., 114.,110.,107.,104.,106.,108.,110.,120.,
     *  126., 135.,140.,144.,146.,148.,152.,150.,146.,
     *  142., 138.,132.,129.,126.,127.,128.,  6*0./
      DATA (CROSEC(I,3,3,1),I=1,51)/  0.,
c     *   0.,  5., 10., 20., 35., 55., 70., 80., 90.,
c     *  4.2, 4.5, 10., 20., 35., 55., 70., 80., 90.,      ! DKFZc
     *   0., 4.5, 10., 20., 35., 55., 70., 80., 90.,      ! DKFZc1
     * 105.,115.,115.,100., 90., 86., 84., 84., 82., 80.,
     *  80., 80., 80., 79., 78., 80., 84., 88., 94.,
     * 100.,105.,108.,108.,108.,112.,114.,114.,112.,
     * 110.,108.,106.,104.,101.,102.,102.,  6*0./
C                                   PROTON
      DATA (CROSEC(I,1,3,2),I=1,51)/  0.,
     *14.,15.,17.,20.,22.,25.,27.,30.,35.,40.,45.,50.,60.,70.,80.,90.,
     *100.,120.,140.,150.,160.,180.,200.,250.,300.,350.,
     *400.,500.,600.,700.,800.,900.,1000.,1500.,2000.,3000.,
     *5000.,7000.,10000.,20000.,50000.,100000.,500000.,1000000.,
     *5*0.,45./
      DATA (CROSEC(I,2,3,2),I=1,51)/  2070.,
     * 1090.,1020.,915.,800.,710.,640.,600.,560.,500.,
     *  440., 390.,360.,295.,256.,220.,192.,168.,136.,120.,
     *  116., 114.,110.,107.,104.,106.,108.,110.,120.,
     *  126., 135.,140.,144.,146.,148.,152.,150.,146.,
     *  142., 138.,132.,129.,126.,127.,128.,  6*0./
      DATA (CROSEC(I,3,3,2),I=1,51)/  0.,
c     *   0.,  2.,  3., 13., 30., 50., 65., 77., 90.,
c     * 1.68, 1.8,  4., 13., 30., 50., 65., 77., 90.,      ! DKFZc
     *   0., 1.8,  4., 13., 30., 50., 65., 77., 90.,      ! DKFZc1
     * 105.,115.,115.,100., 90., 86., 84., 84., 82., 80.,
     *  80., 80., 80., 79., 78., 80., 84., 88., 94.,
     * 100.,105.,108.,108.,108.,112.,114.,114.,112.,
     * 110.,108.,106.,104.,101.,102.,102.,  6*0./
C                                   PI-
      DATA (CROSEC(I,1,3,3),I=1,51)/
     *20.,40.,60.,80.,100.,120.,130.,140.,150.,160.,170.,180.,190.,
     *200.,220.,240.,260.,280.,300.,350.,400.,450.,500.,550.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  12*0.,38./
      DATA (CROSEC(I,2,3,3),I=1,51)/
     *  40., 70.,108.0,152.,208.0,276.,300.,320.,329.,333.,
     * 332.,328.,322.0,310.,288.0,260.,240.,216.,196.,144.,
     * 125.,112.,108.5,109.,110.5,117.,123.,128.5,135.,110.,
     *  96., 87., 85.0,83.5, 83.5,83.5,83.5, 83.5,  13*0./
      DATA (CROSEC(I,3,3,3),I=1,51)/
     *  18., 38., 62., 98.,136.,176.,190.,200.,209.,212.,
     * 212.,208.,204.,196.,176.,164.,150.,134.,124.,97.5,
     *  90., 85.,82.5,83.5,86.5, 93.,97.5,100.,102., 83.,
     *  77., 75., 74.,72.5,72.5,72.5,72.5,72.5,  13*0./
C                                   PI+
      DATA (CROSEC(I,1,3,4),I=1,51)/
     *20.,40.,60.,80.,100.,120.,130.,140.,150.,160.,170.,180.,190.,
     *200.,220.,240.,260.,280.,300.,350.,400.,450.,500.,550.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  12*0.,38./
      DATA (CROSEC(I,2,3,4),I=1,51)/
     *  40., 70.,108.0,152.,208.0,276.,300.,320.,329.,333.,
     * 332.,328.,322.0,310.,288.0,260.,240.,216.,196.,144.,
     * 125.,112.,108.5,109.,110.5,117.,123.,128.5,135.,110.,
     *  96., 87., 85.0,83.5, 83.5,83.5,83.5, 83.5,  13*0./
      DATA (CROSEC(I,3,3,4),I=1,51)/
     *  18., 38., 62., 98.,136.,176.,190.,200.,209.,212.,
     * 212.,208.,204.,196.,176.,164.,150.,134.,124.,97.5,
     *  90., 85.,82.5,83.5,86.5, 93.,97.5,100.,102., 83.,
     *  77., 75., 74.,72.5,72.5,72.5,72.5,72.5,  13*0./
C                                   ANTINEUTRON
      DATA (CROSEC(I,1,3,5),I=1,51)/ 0.,
     *14.,15.,17.,20.,22.,25.,27.,30.,35.,40.,45.,50.,60.,70.,80.,90.,
     *100.,120.,140.,150.,160.,180.,200.,250.,300.,350.,
     *400.,500.,600.,700.,800.,900.,1000.,1500.,2000.,3000.,
     *5000.,7000.,10000.,20000.,50000.,100000.,500000.,1000000.,
     *5*0.,45./
      DATA (CROSEC(I,2,3,5),I=1,51)/ 2070.,
     * 1090.,1020.,915.,800.,710.,640.,600.,560.,500.,
     *  440., 390.,360.,295.,256.,220.,192.,168.,136.,120.,
     *  116., 114.,110.,107.,104.,106.,108.,110.,120.,
     *  126., 135.,140.,144.,146.,148.,152.,150.,146.,
     *  142., 138.,132.,129.,126.,127.,128.,  6*0./
      DATA (CROSEC(I,3,3,5),I=1,51)/  0.,
     *   0.,  5., 10., 20., 35., 55., 70., 80., 90.,
     * 105.,115.,115.,100., 90., 86., 84., 84., 82., 80.,
     *  80., 80., 80., 79., 78., 80., 84., 88., 94.,
     * 100.,105.,108.,108.,108.,112.,114.,114.,112.,
     * 110.,108.,106.,104.,101.,102.,102.,  6*0./
C                                   ANTIPROTON
      DATA (CROSEC(I,1,3,6),I=1,51)/  0.,
     *14.,15.,17.,20.,22.,25.,27.,30.,35.,40.,45.,50.,60.,70.,80.,90.,
     *100.,120.,140.,150.,160.,180.,200.,250.,300.,350.,
     *400.,500.,600.,700.,800.,900.,1000.,1500.,2000.,3000.,
     *5000.,7000.,10000.,20000.,50000.,100000.,500000.,1000000.,
     *5*0.,45./
      DATA (CROSEC(I,2,3,6),I=1,51)/  2070.,
     * 1090.,1020.,915.,800.,710.,640.,600.,560.,500.,
     *  440., 390.,360.,295.,256.,220.,192.,168.,136.,120.,
     *  116., 114.,110.,107.,104.,106.,108.,110.,120.,
     *  126., 135.,140.,144.,146.,148.,152.,150.,146.,
     *  142., 138.,132.,129.,126.,127.,128.,  6*0./
      DATA (CROSEC(I,3,3,6),I=1,51)/  0.,
     *   0.,  2.,  3., 13., 30., 50., 65., 77., 90.,
     * 105.,115.,115.,100., 90., 86., 84., 84., 82., 80.,
     *  80., 80., 80., 79., 78., 80., 84., 88., 94.,
     * 100.,105.,108.,108.,108.,112.,114.,114.,112.,
     * 110.,108.,106.,104.,101.,102.,102.,  6*0./
C                                   K-
      DATA (CROSEC(I,1,3,7),I=1,51)/
     *20.,40.,60.,80.,100.,120.,130.,140.,150.,160.,170.,180.,190.,
     *200.,220.,240.,260.,280.,300.,350.,400.,450.,500.,550.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  12*0.,38./
      DATA (CROSEC(I,2,3,7),I=1,51)/
     *  40., 70.,108.0,152.,208.0,276.,300.,320.,329.,333.,
     * 332.,328.,322.0,310.,288.0,260.,240.,216.,196.,144.,
     * 125.,112.,108.5,109.,110.5,117.,123.,128.5,135.,110.,
     *  96., 87., 85.0,83.5, 83.5,83.5,83.5, 83.5,  13*0./
      DATA (CROSEC(I,3,3,7),I=1,51)/
     *  18., 38., 62., 98.,136.,176.,190.,200.,209.,212.,
     * 212.,208.,204.,196.,176.,164.,150.,134.,124.,97.5,
     *  90., 85.,82.5,83.5,86.5, 93.,97.5,100.,102., 83.,
     *  77., 75., 74.,72.5,72.5,72.5,72.5,72.5,  13*0./
C                                   K+
      DATA (CROSEC(I,1,3,8),I=1,51)/
     *20.,40.,60.,80.,100.,120.,130.,140.,150.,160.,170.,180.,190.,
     *200.,220.,240.,260.,280.,300.,350.,400.,450.,500.,550.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  12*0.,38./
      DATA (CROSEC(I,2,3,8),I=1,51)/
     *  40., 70.,108.0,152.,208.0,276.,300.,320.,329.,333.,
     * 332.,328.,322.0,310.,288.0,260.,240.,216.,196.,144.,
     * 125.,112.,108.5,109.,110.5,117.,123.,128.5,135.,110.,
     *  96., 87., 85.0,83.5, 83.5,83.5,83.5, 83.5,  13*0./
      DATA (CROSEC(I,3,3,8),I=1,51)/
     *  18., 38., 62., 98.,136.,176.,190.,200.,209.,212.,
     * 212.,208.,204.,196.,176.,164.,150.,134.,124.,97.5,
     *  90., 85.,82.5,83.5,86.5, 93.,97.5,100.,102., 83.,
     *  77., 75., 74.,72.5,72.5,72.5,72.5,72.5,  13*0./
C                                   K0
      DATA (CROSEC(I,1,3,9),I=1,51)/
     *20.,40.,60.,80.,100.,120.,130.,140.,150.,160.,170.,180.,190.,
     *200.,220.,240.,260.,280.,300.,350.,400.,450.,500.,550.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  12*0.,38./
      DATA (CROSEC(I,2,3,9),I=1,51)/
     *  40., 70.,108.0,152.,208.0,276.,300.,320.,329.,333.,
     * 332.,328.,322.0,310.,288.0,260.,240.,216.,196.,144.,
     * 125.,112.,108.5,109.,110.5,117.,123.,128.5,135.,110.,
     *  96., 87., 85.0,83.5, 83.5,83.5,83.5, 83.5,  13*0./
      DATA (CROSEC(I,3,3,9),I=1,51)/
     *  18., 38., 62., 98.,136.,176.,190.,200.,209.,212.,
     * 212.,208.,204.,196.,176.,164.,150.,134.,124.,97.5,
     *  90., 85.,82.5,83.5,86.5, 93.,97.5,100.,102., 83.,
     *  77., 75., 74.,72.5,72.5,72.5,72.5,72.5,  13*0./
C                             ANTI K0
      DATA (CROSEC(I,1,3,10),I=1,51)/
     *20.,40.,60.,80.,100.,120.,130.,140.,150.,160.,170.,180.,190.,
     *200.,220.,240.,260.,280.,300.,350.,400.,450.,500.,550.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  12*0.,38./
      DATA (CROSEC(I,2,3,10),I=1,51)/
     *  40., 70.,108.0,152.,208.0,276.,300.,320.,329.,333.,
     * 332.,328.,322.0,310.,288.0,260.,240.,216.,196.,144.,
     * 125.,112.,108.5,109.,110.5,117.,123.,128.5,135.,110.,
     *  96., 87., 85.0,83.5, 83.5,83.5,83.5, 83.5,  13*0./
      DATA (CROSEC(I,3,3,10),I=1,51)/
     *  18., 38., 62., 98.,136.,176.,190.,200.,209.,212.,
     * 212.,208.,204.,196.,176.,164.,150.,134.,124.,97.5,
     *  90., 85.,82.5,83.5,86.5, 93.,97.5,100.,102., 83.,
     *  77., 75., 74.,72.5,72.5,72.5,72.5,72.5,  13*0./
C
C******************** DATA FOR BE(4,9), INUCL=5 ***********
C                                   NEUTRON
      DATA (CROSEC(I,1,5,1),I=1,51)/
     *14.,15.,17.,20.,22.,25.,27.,30.,35.,40.,45.,50.,60.,70.,80.,90.,
     *100.,120.,140.,150.,160.,180.,200.,250.,300.,350.,
     *400.,500.,600.,700.,800.,900.,1000.,1500.,2000.,3000.,
     *5000.,7000.,10000.,20000.,50000.,100000.,500000.,1000000.,
     *6*0.,44./
      DATA (CROSEC(I,2,5,1),I=1,51)/
     * 1490.,1460.,1400.,1350.,1270.,1200.,1160.,1100.,1000.,
     *  910., 810., 740., 625., 575., 455., 406., 365., 310., 
     *  275., 262., 255., 240., 235., 225., 225., 230., 238.,
     *  252., 270., 282., 288., 290., 294., 303., 303., 300., 292.,
     *  284., 277., 267., 263., 264., 268., 268.,  7*0./
      DATA (CROSEC(I,3,5,1),I=1,51)/
     * 650.,640.,617.,595.,555.,520.,495.,470.,430.,
     * 385.,350.,320.,270.,250.,210.,190.,185.,178.,
     * 175.,175.,175.,175.,175.,170.,170.,172.,176.,
     * 184.,194.,200.,209.,213.,214.,216.,216.,212.,
     * 210.,210.,210.,210.,210.,210.,210.,210., 7*0./
C                                   PROTON
      DATA (CROSEC(I,1,5,2),I=1,51)/
     *14.,15.,17.,20.,22.,25.,27.,30.,35.,40.,45.,50.,60.,70.,80.,90.,
     *100.,120.,140.,150.,160.,180.,200.,250.,300.,350.,
     *400.,500.,600.,700.,800.,900.,1000.,1500.,2000.,3000.,
     *5000.,7000.,10000.,20000.,50000.,100000.,500000.,1000000.,
     *6*0.,44./
      DATA (CROSEC(I,2,5,2),I=1,51)/
     * 1490.,1460.,1400.,1350.,1270.,1200.,1160.,1100.,1000.,
     *  910., 810., 740., 625., 575., 455., 406., 365., 310., 
     *  275., 262., 255., 240., 235., 225., 225., 230., 238.,
     *  252., 270., 282., 288., 290., 294., 303., 303., 300., 292.,
     *  284., 277., 267., 263., 264., 268., 268.,  7*0./
      DATA (CROSEC(I,3,5,2),I=1,51)/
     * 490.,540.,580.,545.,525.,495.,470.,450.,420.,
     * 370.,340.,310.,262.,242.,205.,185.,180.,175.,
     * 172.,175.,175.,175.,175.,170.,170.,172.,176.,
     * 184.,194.,200.,209.,213.,214.,216.,216.,212.,
     * 210.,210.,210.,210.,210.,210.,210.,210.,  7*0./
C                                   PI-
      DATA (CROSEC(I,1,5,3),I=1,51)/
     *20.,40.,60.,80.,100.,120.,130.,140.,150.,160.,170.,180.,190.,
     *200.,220.,240.,260.,280.,300.,350.,400.,450.,500.,550.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  12*0.,38./
      DATA (CROSEC(I,2,5,3),I=1,51)/
     * 150.,210.,294.,396.,520.,600.,623.,635.,642.,640.,
     * 630.,615.,600.,576.,540.,504.,470.,435.,400.,340.,
     * 294.,258.,236.,230.,233.,244.,257.,270.,276.,250.,
     * 230.,215.,205.,194.,188.,186.,186.,186., 13*0./
      DATA (CROSEC(I,3,5,3),I=1,51)/
     *  90.,126.,177.,240.,320.,380.,400.,410.,414.,410.,
     * 400.,387.,371.,360.,333.,312.,285.,260.,237.,216.,
     * 198.,187.,182.,180.,182.,187.,193.,202.,207.,179.,
     * 172.,165.,159.,155.,144.,144.,144.,144.,  13*0./
C                                   PI+
      DATA (CROSEC(I,1,5,4),I=1,51)/
     *20.,40.,60.,80.,100.,120.,130.,140.,150.,160.,170.,180.,190.,
     *200.,220.,240.,260.,280.,300.,350.,400.,450.,500.,550.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  12*0.,38./
      DATA (CROSEC(I,2,5,4),I=1,51)/
     *  96.,150.,222.,320.,430.,514.,545.,565.,574.,574.,
     * 564.,552.,535.,522.,490.,462.,432.,398.,367.,314.,
     * 276.,248.,232.,230.,233.,244.,257.,270.,276.,250.,
     * 230.,215.,205.,194.,188.,186.,186.,186.,  13*0./
      DATA (CROSEC(I,3,5,4),I=1,51)/
     *  60., 95.,142.,194.,262.,319.,345.,361.,364.,364.,
     * 354.,350.,330.,319.,298.,280.,258.,237.,216.,200.,
     * 189.,183.,182.,180.,182.,187.,193.,202.,207.,179.,
     * 172.,165.,159.,155.,144.,144.,144.,144.,  13*0./
C                                   ANTINEUTRON
      DATA (CROSEC(I,1,5,5),I=1,51)/
     *14.,15.,17.,20.,22.,25.,27.,30.,35.,40.,45.,50.,60.,70.,80.,90.,
     *100.,120.,140.,150.,160.,180.,200.,250.,300.,350.,
     *400.,500.,600.,700.,800.,900.,1000.,1500.,2000.,3000.,
     *5000.,7000.,10000.,20000.,50000.,100000.,500000.,1000000.,
     *6*0.,44./
      DATA (CROSEC(I,2,5,5),I=1,51)/
     * 1490.,1460.,1400.,1350.,1270.,1200.,1160.,1100.,1000.,
     *  910., 810., 740., 625., 575., 455., 406., 365., 310., 
     *  275., 262., 255., 240., 235., 225., 225., 230., 238.,
     *  252., 270., 282., 288., 290., 294., 303., 303., 300., 292.,
     *  284., 277., 267., 263., 264., 268., 268.,  7*0./
      DATA (CROSEC(I,3,5,5),I=1,51)/
     * 650.,640.,617.,595.,555.,520.,495.,470.,430.,
     * 385.,350.,320.,270.,250.,210.,190.,185.,178.,
     * 175.,175.,175.,175.,175.,170.,170.,172.,176.,
     * 184.,194.,200.,209.,213.,214.,216.,216.,212.,
     * 210.,210.,210.,210.,210.,210.,210.,210., 7*0./
C                                   ANTIPROTON
      DATA (CROSEC(I,1,5,6),I=1,51)/
     *14.,15.,17.,20.,22.,25.,27.,30.,35.,40.,45.,50.,60.,70.,80.,90.,
     *100.,120.,140.,150.,160.,180.,200.,250.,300.,350.,
     *400.,500.,600.,700.,800.,900.,1000.,1500.,2000.,3000.,
     *5000.,7000.,10000.,20000.,50000.,100000.,500000.,1000000.,
     *6*0.,44./
      DATA (CROSEC(I,2,5,6),I=1,51)/
     * 1490.,1460.,1400.,1350.,1270.,1200.,1160.,1100.,1000.,
     *  910., 810., 740., 625., 575., 455., 406., 365., 310., 
     *  275., 262., 255., 240., 235., 225., 225., 230., 238.,
     *  252., 270., 282., 288., 290., 294., 303., 303., 300., 292.,
     *  284., 277., 267., 263., 264., 268., 268.,  7*0./
      DATA (CROSEC(I,3,5,6),I=1,51)/
     * 490.,540.,580.,545.,525.,495.,470.,450.,420.,
     * 370.,340.,310.,262.,242.,205.,185.,180.,175.,
     * 172.,175.,175.,175.,175.,170.,170.,172.,176.,
     * 184.,194.,200.,209.,213.,214.,216.,216.,212.,
     * 210.,210.,210.,210.,210.,210.,210.,210.,  7*0./
C                                   K-
      DATA (CROSEC(I,1,5,7),I=1,51)/
     *20.,40.,60.,80.,100.,120.,130.,140.,150.,160.,170.,180.,190.,
     *200.,220.,240.,260.,280.,300.,350.,400.,450.,500.,550.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  12*0.,38./
      DATA (CROSEC(I,2,5,7),I=1,51)/
     * 150.,210.,294.,396.,520.,600.,623.,635.,642.,640.,
     * 630.,615.,600.,576.,540.,504.,470.,435.,400.,340.,
     * 294.,258.,236.,230.,233.,244.,257.,270.,276.,250.,
     * 230.,215.,205.,194.,188.,186.,186.,186., 13*0./
      DATA (CROSEC(I,3,5,7),I=1,51)/
     *  90.,126.,177.,240.,320.,380.,400.,410.,414.,410.,
     * 400.,387.,371.,360.,333.,312.,285.,260.,237.,216.,
     * 198.,187.,182.,180.,182.,187.,193.,202.,207.,179.,
     * 172.,165.,159.,155.,144.,144.,144.,144.,  13*0./
C                                   K+
      DATA (CROSEC(I,1,5,8),I=1,51)/
     *20.,40.,60.,80.,100.,120.,130.,140.,150.,160.,170.,180.,190.,
     *200.,220.,240.,260.,280.,300.,350.,400.,450.,500.,550.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  12*0.,38./
      DATA (CROSEC(I,2,5,8),I=1,51)/
     *  96.,150.,222.,320.,430.,514.,545.,565.,574.,574.,
     * 564.,552.,535.,522.,490.,462.,432.,398.,367.,314.,
     * 276.,248.,232.,230.,233.,244.,257.,270.,276.,250.,
     * 230.,215.,205.,194.,188.,186.,186.,186.,  13*0./
      DATA (CROSEC(I,3,5,8),I=1,51)/
     *  60., 95.,142.,194.,262.,319.,345.,361.,364.,364.,
     * 354.,350.,330.,319.,298.,280.,258.,237.,216.,200.,
     * 189.,183.,182.,180.,182.,187.,193.,202.,207.,179.,
     * 172.,165.,159.,155.,144.,144.,144.,144.,  13*0./
C                                   K0
      DATA (CROSEC(I,1,5,9),I=1,51)/
     *20.,40.,60.,80.,100.,120.,130.,140.,150.,160.,170.,180.,190.,
     *200.,220.,240.,260.,280.,300.,350.,400.,450.,500.,550.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  12*0.,38./
      DATA (CROSEC(I,2,5,9),I=1,51)/
     *  96.,150.,222.,320.,430.,514.,545.,565.,574.,574.,
     * 564.,552.,535.,522.,490.,462.,432.,398.,367.,314.,
     * 276.,248.,232.,230.,233.,244.,257.,270.,276.,250.,
     * 230.,215.,205.,194.,188.,186.,186.,186.,  13*0./
      DATA (CROSEC(I,3,5,9),I=1,51)/
     *  60., 95.,142.,194.,262.,319.,345.,361.,364.,364.,
     * 354.,350.,330.,319.,298.,280.,258.,237.,216.,200.,
     * 189.,183.,182.,180.,182.,187.,193.,202.,207.,179.,
     * 172.,165.,159.,155.,144.,144.,144.,144.,  13*0./
C                              ANTI K0
      DATA (CROSEC(I,1,5,10),I=1,51)/
     *20.,40.,60.,80.,100.,120.,130.,140.,150.,160.,170.,180.,190.,
     *200.,220.,240.,260.,280.,300.,350.,400.,450.,500.,550.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  12*0.,38./
      DATA (CROSEC(I,2,5,10),I=1,51)/
     * 150.,210.,294.,396.,520.,600.,623.,635.,642.,640.,
     * 630.,615.,600.,576.,540.,504.,470.,435.,400.,340.,
     * 294.,258.,236.,230.,233.,244.,257.,270.,276.,250.,
     * 230.,215.,205.,194.,188.,186.,186.,186., 13*0./
      DATA (CROSEC(I,3,5,10),I=1,51)/
     *  90.,126.,177.,240.,320.,380.,400.,410.,414.,410.,
     * 400.,387.,371.,360.,333.,312.,285.,260.,237.,216.,
     * 198.,187.,182.,180.,182.,187.,193.,202.,207.,179.,
     * 172.,165.,159.,155.,144.,144.,144.,144.,  13*0./
C
C ********************* DATA FOR C(6,12), INUCL=6 **************
C                                   NEUTRON
      DATA (CROSEC(I,1,6,1),I=1,51)/
     *14.,15.,17.,20.,22.,25.,27.,30.,35.,40.,45.,50.,60.,70.,80.,90.,
     *100.,120.,140.,150.,160.,180.,200.,250.,300.,350.,
     *400.,500.,600.,700.,800.,900.,1000.,1500.,2000.,3000.,
     *5000.,7000.,10000.,20000.,50000.,100000.,500000.,1000000.,
     *6*0.,44./
      DATA (CROSEC(I,2,6,1),I=1,51)/
     * 1300.,1370.,1450.,1455.,1445.,1385.,1345.,1290.,1210.,
C     * 1240.,1370.,1450.,1455.,1445.,1385.,1345.,1290.,1210., Sobol 01.10.99
     * 1110.,1020., 940., 800., 700., 604., 530., 475., 396.,
     *  350., 336., 320., 303., 294., 280., 280., 286., 296.,
     *  314., 330., 344., 356., 360., 364., 384., 388., 384.,
     *  364., 352., 344., 330., 324., 324., 332., 332.,  7*0./
      DATA (CROSEC(I,3,6,1),I=1,51)/
     * 590.,570.,542.,510.,500.,460.,445.,430.,395.,
     * 380.,350.,330.,295.,270.,255.,240.,228.,222.,
     * 216.,216.,210.,210.,210.,208.,210.,214.,216.,
     * 228.,240.,248.,254.,257.,260.,262.,260.,256.,
     * 252.,252.,250.,250.,248.,248.,248.,248., 7*0./
C                                   PROTON
      DATA (CROSEC(I,1,6,2),I=1,51)/
     *14.,15.,17.,20.,22.,25.,27.,30.,35.,40.,45.,50.,60.,70.,80.,90.,
     *100.,120.,140.,150.,160.,180.,200.,250.,300.,350.,
     *400.,500.,600.,700.,800.,900.,1000.,1500.,2000.,3000.,
     *5000.,7000.,10000.,20000.,50000.,100000.,500000.,1000000.,
     *6*0.,44./
      DATA (CROSEC(I,2,6,2),I=1,51)/
     * 1300.,1370.,1450.,1455.,1445.,1385.,1345.,1290.,1210.,
C     * 1240.,1370.,1450.,1455.,1445.,1385.,1345.,1290.,1210., Sobol 01.10.99
     * 1110.,1020., 940., 800., 700., 604., 530., 475., 396.,
     *  350., 336., 320., 303., 294., 280., 280., 286., 296.,
     *  314., 330., 344., 356., 360., 364., 384., 388., 384.,
     *  364., 352., 344., 330., 324., 324., 332., 332.,  7*0./
      DATA (CROSEC(I,3,6,2),I=1,51)/
     * 310.,330.,400.,440.,450.,435.,430.,420.,385.,
     * 370.,340.,320.,288.,263.,249.,234.,222.,216.,
     * 210.,211.,205.,208.,210.,208.,210.,214.,216.,
     * 228.,240.,248.,254.,257.,260.,262.,260.,256.,
     * 252.,252.,250.,250.,248.,248.,248.,248.,  7*0./
C                                   PI-
      DATA (CROSEC(I,1,6,3),I=1,51)/
     *20.,40.,60.,80.,100.,110.,120.,130.,140.,150.,160.,170.,180.,
     *200.,220.,240.,260.,280.,300.,350.,400.,450.,500.,550.,575.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  11*0.,39./
      DATA (CROSEC(I,2,6,3),I=1,51)/
     * 204.,260.,366.,517.,630.,673.,694.,704.,710.,711.,
     * 706.,694.,676.,648.,616.,584.,548.,518.,489.,426., 
     * 376.,342.,323.,310.,312.,313.,319.,333.,342.,348.,
     * 310.,290.,268.,250.,245.,237.,234.,234.,234.,  12*0./
      DATA (CROSEC(I,3,6,3),I=1,51)/
     * 128.,160.,224.,315.,388.,416.,430.,438.,444.,445.,
     * 440.,432.,416.,400.,380.,354.,320.,304.,288.,264.,
     * 246.,240.,233.,232.,233.,234.,238.,246.,252.,256.,
     * 220.,210.,198.,187.,183.,176.,174.,174.,174.,  12*0./
C                                   PI+
      DATA (CROSEC(I,1,6,4),I=1,51)/
     *20.,40.,60.,80.,100.,110.,120.,130.,140.,150.,160.,170.,180.,
     *200.,220.,240.,260.,280.,300.,350.,400.,450.,500.,550.,575.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  11*0.,39./
      DATA (CROSEC(I,2,6,4),I=1,51)/
     * 140.,192.,294.,428.,594.,642.,662.,678.,685.,688.,
     * 684.,672.,656.,630.,598.,567.,533.,504.,474.,416.,
     * 369.,336.,319.,310.,312.,313.,319.,333.,342.,348.,
     * 310.,290.,268.,250.,245.,237.,234.,234.,234.,  12*0./
      DATA (CROSEC(I,3,6,4),I=1,51)/
     *  94.,132.,184.,260.,370.,398.,408.,420.,426.,428.,
     * 424.,416.,400.,386.,366.,340.,308.,294.,280.,257.,
     * 241.,236.,231.,232.,233.,234.,238.,246.,252.,256.,
     * 220.,210.,198.,187.,183.,176.,174.,174.,174.,  12*0./
C                                   ANTINEUTRON
      DATA (CROSEC(I,1,6,5),I=1,51)/
     *14.,15.,17.,20.,22.,25.,27.,30.,35.,40.,45.,50.,60.,70.,80.,90.,
     *100.,120.,140.,150.,160.,180.,200.,250.,300.,350.,
     *400.,500.,600.,700.,800.,900.,1000.,1500.,2000.,3000.,
     *5000.,7000.,10000.,20000.,50000.,100000.,500000.,1000000.,
     *6*0.,44./
      DATA (CROSEC(I,2,6,5),I=1,51)/
     * 1300.,1370.,1450.,1455.,1445.,1385.,1345.,1290.,1210.,
C     * 1240.,1370.,1450.,1455.,1445.,1385.,1345.,1290.,1210., Sobol 01.10.99
     * 1110.,1020., 940., 800., 700., 604., 530., 475., 396.,
     *  350., 336., 320., 303., 294., 280., 280., 286., 296.,
     *  314., 330., 344., 356., 360., 364., 384., 388., 384.,
     *  364., 352., 344., 330., 324., 324., 332., 332.,  7*0./
      DATA (CROSEC(I,3,6,5),I=1,51)/
     * 590.,570.,542.,510.,500.,460.,445.,430.,395.,
     * 380.,350.,330.,295.,270.,255.,240.,228.,222.,
     * 216.,216.,210.,210.,210.,208.,210.,214.,216.,
     * 228.,240.,248.,254.,257.,260.,262.,260.,256.,
     * 252.,252.,250.,250.,248.,248.,248.,248., 7*0./
C                                   ANTIPROTON
      DATA (CROSEC(I,1,6,6),I=1,51)/
     *14.,15.,17.,20.,22.,25.,27.,30.,35.,40.,45.,50.,60.,70.,80.,90.,
     *100.,120.,140.,150.,160.,180.,200.,250.,300.,350.,
     *400.,500.,600.,700.,800.,900.,1000.,1500.,2000.,3000.,
     *5000.,7000.,10000.,20000.,50000.,100000.,500000.,1000000.,
     *6*0.,44./
      DATA (CROSEC(I,2,6,6),I=1,51)/
     * 1300.,1370.,1450.,1455.,1445.,1385.,1345.,1290.,1210.,
C     * 1240.,1370.,1450.,1455.,1445.,1385.,1345.,1290.,1210., Sobol 01.10.99
     * 1110.,1020., 940., 800., 700., 604., 530., 475., 396.,
     *  350., 336., 320., 303., 294., 280., 280., 286., 296.,
     *  314., 330., 344., 356., 360., 364., 384., 388., 384.,
     *  364., 352., 344., 330., 324., 324., 332., 332.,  7*0./
      DATA (CROSEC(I,3,6,6),I=1,51)/
     * 310.,330.,400.,440.,450.,435.,430.,420.,385.,
     * 370.,340.,320.,288.,263.,249.,234.,222.,216.,
     * 210.,211.,205.,208.,210.,208.,210.,214.,216.,
     * 228.,240.,248.,254.,257.,260.,262.,260.,256.,
     * 252.,252.,250.,250.,248.,248.,248.,248.,  7*0./
C                                   K-
      DATA (CROSEC(I,1,6,7),I=1,51)/
     *20.,40.,60.,80.,100.,110.,120.,130.,140.,150.,160.,170.,180.,
     *200.,220.,240.,260.,280.,300.,350.,400.,450.,500.,550.,575.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  11*0.,39./
      DATA (CROSEC(I,2,6,7),I=1,51)/
     * 204.,260.,366.,517.,630.,673.,694.,704.,710.,711.,
     * 706.,694.,676.,648.,616.,584.,548.,518.,489.,426., 
     * 376.,342.,323.,310.,312.,313.,319.,333.,342.,348.,
     * 310.,290.,268.,250.,245.,237.,234.,234.,234.,  12*0./
      DATA (CROSEC(I,3,6,7),I=1,51)/
     * 128.,160.,224.,315.,388.,416.,430.,438.,444.,445.,
     * 440.,432.,416.,400.,380.,354.,320.,304.,288.,264.,
     * 246.,240.,233.,232.,233.,234.,238.,246.,252.,256.,
     * 220.,210.,198.,187.,183.,176.,174.,174.,174.,  12*0./
C                                   K+
      DATA (CROSEC(I,1,6,8),I=1,51)/
     *20.,40.,60.,80.,100.,110.,120.,130.,140.,150.,160.,170.,180.,
     *200.,220.,240.,260.,280.,300.,350.,400.,450.,500.,550.,575.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  11*0.,39./
      DATA (CROSEC(I,2,6,8),I=1,51)/
     * 140.,192.,294.,428.,594.,642.,662.,678.,685.,688.,
     * 684.,672.,656.,630.,598.,567.,533.,504.,474.,416.,
     * 369.,336.,319.,310.,312.,313.,319.,333.,342.,348.,
     * 310.,290.,268.,250.,245.,237.,234.,234.,234.,  12*0./
      DATA (CROSEC(I,3,6,8),I=1,51)/
     *  94.,132.,184.,260.,370.,398.,408.,420.,426.,428.,
     * 424.,416.,400.,386.,366.,340.,308.,294.,280.,257.,
     * 241.,236.,231.,232.,233.,234.,238.,246.,252.,256.,
     * 220.,210.,198.,187.,183.,176.,174.,174.,174.,  12*0./
C                                   K0
      DATA (CROSEC(I,1,6,9),I=1,51)/
     *20.,40.,60.,80.,100.,110.,120.,130.,140.,150.,160.,170.,180.,
     *200.,220.,240.,260.,280.,300.,350.,400.,450.,500.,550.,575.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  11*0.,39./
      DATA (CROSEC(I,2,6,9),I=1,51)/
     * 140.,192.,294.,428.,594.,642.,662.,678.,685.,688.,
     * 684.,672.,656.,630.,598.,567.,533.,504.,474.,416.,
     * 369.,336.,319.,310.,312.,313.,319.,333.,342.,348.,
     * 310.,290.,268.,250.,245.,237.,234.,234.,234.,  12*0./
      DATA (CROSEC(I,3,6,9),I=1,51)/
     *  94.,132.,184.,260.,370.,398.,408.,420.,426.,428.,
     * 424.,416.,400.,386.,366.,340.,308.,294.,280.,257.,
     * 241.,236.,231.,232.,233.,234.,238.,246.,252.,256.,
     * 220.,210.,198.,187.,183.,176.,174.,174.,174.,  12*0./
C                              ANTI K0
      DATA (CROSEC(I,1,6,10),I=1,51)/
     *20.,40.,60.,80.,100.,110.,120.,130.,140.,150.,160.,170.,180.,
     *200.,220.,240.,260.,280.,300.,350.,400.,450.,500.,550.,575.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  11*0.,39./
      DATA (CROSEC(I,2,6,10),I=1,51)/
     * 204.,260.,366.,517.,630.,673.,694.,704.,710.,711.,
     * 706.,694.,676.,648.,616.,584.,548.,518.,489.,426., 
     * 376.,342.,323.,310.,312.,313.,319.,333.,342.,348.,
     * 310.,290.,268.,250.,245.,237.,234.,234.,234.,  12*0./
      DATA (CROSEC(I,3,6,10),I=1,51)/
     * 128.,160.,224.,315.,388.,416.,430.,438.,444.,445.,
     * 440.,432.,416.,400.,380.,354.,320.,304.,288.,264.,
     * 246.,240.,233.,232.,233.,234.,238.,246.,252.,256.,
     * 220.,210.,198.,187.,183.,176.,174.,174.,174.,  12*0./
C
C ********************* DATA FOR N(7,14), INUCL=7 ****************
C                                   NEUTRON
      DATA (CROSEC(I,1,7,1),I=1,51)/
     *14.,15.,17.,20.,22.,25.,27.,30.,35.,40.,45.,50.,60.,70.,80.,90.,
     *100.,120.,140.,150.,160.,180.,200.,250.,300.,350.,
     *400.,500.,600.,700.,800.,900.,1000.,1500.,2000.,3000.,
     *5000.,7000.,10000.,20000.,50000.,100000.,500000.,1000000.,
     *6*0.,44./
      DATA (CROSEC(I,2,7,1),I=1,51)/
     * 1420.,1480.,1537.,1550.,1525.,1500.,1480.,1425.,1340.,
     * 1260.,1175.,1090., 930., 805., 690., 612., 552., 462., 
     *  402., 384., 372., 350., 345., 326., 324., 328., 336.,
     *  356., 372., 388., 400., 408., 415., 430., 435., 432.,
     *  415., 402., 390., 375., 367., 370., 382., 385.,  7*0./
      DATA (CROSEC(I,3,7,1),I=1,51)/
     * 680.,665.,625.,580.,562.,525.,510.,485.,450.,
     * 435.,410.,387.,340.,310.,290.,280.,276.,274.,
     * 260.,258.,254.,247.,245.,240.,240.,244.,250.,
     * 260.,268.,275.,280.,285.,290.,295.,300.,294.,
     * 292.,290.,285.,285.,282.,282.,282.,282.,  7*0./
C                                   PROTON
      DATA (CROSEC(I,1,7,2),I=1,51)/
     *14.,15.,17.,20.,22.,25.,27.,30.,35.,40.,45.,50.,60.,70.,80.,90.,
     *100.,120.,140.,150.,160.,180.,200.,250.,300.,350.,
     *400.,500.,600.,700.,800.,900.,1000.,1500.,2000.,3000.,
     *5000.,7000.,10000.,20000.,50000.,100000.,500000.,1000000.,
     *6*0.,44./
      DATA (CROSEC(I,2,7,2),I=1,51)/
     * 1420.,1480.,1537.,1550.,1525.,1500.,1480.,1425.,1340.,
     * 1260.,1175.,1090., 930., 805., 690., 612., 552., 462., 
     *  402., 384., 372., 350., 345., 326., 324., 328., 336.,
     *  356., 372., 388., 400., 408., 415., 430., 435., 432.,
     *  415., 402., 390., 375., 367., 370., 382., 385.,  7*0./
      DATA (CROSEC(I,3,7,2),I=1,51)/
     * 420.,440.,470.,490.,497.,500.,480.,462.,440.,
     * 425.,400.,377.,333.,303.,284.,274.,270.,268.,
     * 254.,252.,247.,245.,245.,240.,240.,244.,250.,
     * 260.,268.,275.,280.,285.,290.,295.,300.,294.,
     * 292.,290.,285.,285.,282.,282.,282.,282.,  7*0./
C                                   PI-
      DATA (CROSEC(I,1,7,3),I=1,51)/
     *20.,40.,60.,80.,100.,110.,120.,130.,140.,150.,160.,170.,180.,
     *200.,220.,240.,260.,280.,300.,350.,400.,450.,500.,550.,575.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  11*0.,39./
      DATA (CROSEC(I,2,7,3),I=1,51)/
     * 246.,308.,424.,590.,729.,776.,800.,821.,822.,817., 
     * 800.,778.,768.,728.,690.,654.,615.,584.,556.,480., 
     * 430.,393.,373.,367.,368.,370.,375.,388.,390.,397.,
     * 364.,337.,310.,291.,275.,268.,268.,268.,268.,  12*0./
      DATA (CROSEC(I,3,7,3),I=1,51)/
     * 155.,188.,256.,360.,456.,492.,512.,526.,526.,520.,
     * 504.,491.,475.,450.,425.,396.,376.,360.,340.,300.,
     * 282.,270.,265.,265.,266.,268.,273.,280.,288.,288.,
     * 256.,237.,226.,218.,208.,202.,202.,202.,202., 12*0./
C                                   PI+
      DATA (CROSEC(I,1,7,4),I=1,51)/
     *20.,40.,60.,80.,100.,110.,120.,130.,140.,150.,160.,170.,180.,
     *200.,220.,240.,260.,280.,300.,350.,400.,450.,500.,550.,575.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  11*0.,39./
      DATA (CROSEC(I,2,7,4),I=1,51)/
     * 150.,212.,328.,500.,680.,735.,762.,781.,782.,779.,
     * 770.,748.,740.,706.,672.,633.,600.,569.,541.,467.,
     * 419.,385.,368.,364.,366.,368.,375.,388.,390.,397.,
     * 364.,337.,310.,291.,275.,268.,268.,268.,268.,  12*0./
      DATA (CROSEC(I,3,7,4),I=1,51)/
     *  90.,140.,208.,300.,426.,467.,490.,504.,504.,500.,
     * 484.,474.,460.,437.,413.,381.,365.,350.,330.,292.,
     * 276.,267.,263.,264.,265.,267.,273.,280.,288.,288.,
     * 256.,237.,226.,218.,208.,202.,202.,202.,202.,  12*0./
C                                   ANTINEUTRON
      DATA (CROSEC(I,1,7,5),I=1,51)/
     *14.,15.,17.,20.,22.,25.,27.,30.,35.,40.,45.,50.,60.,70.,80.,90.,
     *100.,120.,140.,150.,160.,180.,200.,250.,300.,350.,
     *400.,500.,600.,700.,800.,900.,1000.,1500.,2000.,3000.,
     *5000.,7000.,10000.,20000.,50000.,100000.,500000.,1000000.,
     *6*0.,44./
      DATA (CROSEC(I,2,7,5),I=1,51)/
     * 1420.,1480.,1537.,1550.,1525.,1500.,1480.,1425.,1340.,
     * 1260.,1175.,1090., 930., 805., 690., 612., 552., 462., 
     *  402., 384., 372., 350., 345., 326., 324., 328., 336.,
     *  356., 372., 388., 400., 408., 415., 430., 435., 432.,
     *  415., 402., 390., 375., 367., 370., 382., 385.,  7*0./
      DATA (CROSEC(I,3,7,5),I=1,51)/
     * 680.,665.,625.,580.,562.,525.,510.,485.,450.,
     * 435.,410.,387.,340.,310.,290.,280.,276.,274.,
     * 260.,258.,254.,247.,245.,240.,240.,244.,250.,
     * 260.,268.,275.,280.,285.,290.,295.,300.,294.,
     * 292.,290.,285.,285.,282.,282.,282.,282.,  7*0./
C                                   ANTIPROTON
      DATA (CROSEC(I,1,7,6),I=1,51)/
     *14.,15.,17.,20.,22.,25.,27.,30.,35.,40.,45.,50.,60.,70.,80.,90.,
     *100.,120.,140.,150.,160.,180.,200.,250.,300.,350.,
     *400.,500.,600.,700.,800.,900.,1000.,1500.,2000.,3000.,
     *5000.,7000.,10000.,20000.,50000.,100000.,500000.,1000000.,
     *6*0.,44./
      DATA (CROSEC(I,2,7,6),I=1,51)/
     * 1420.,1480.,1537.,1550.,1525.,1500.,1480.,1425.,1340.,
     * 1260.,1175.,1090., 930., 805., 690., 612., 552., 462., 
     *  402., 384., 372., 350., 345., 326., 324., 328., 336.,
     *  356., 372., 388., 400., 408., 415., 430., 435., 432.,
     *  415., 402., 390., 375., 367., 370., 382., 385.,  7*0./
      DATA (CROSEC(I,3,7,6),I=1,51)/
     * 420.,440.,470.,490.,497.,500.,480.,462.,440.,
     * 425.,400.,377.,333.,303.,284.,274.,270.,268.,
     * 254.,252.,247.,245.,245.,240.,240.,244.,250.,
     * 260.,268.,275.,280.,285.,290.,295.,300.,294.,
     * 292.,290.,285.,285.,282.,282.,282.,282.,  7*0./
C                                   K-
      DATA (CROSEC(I,1,7,7),I=1,51)/
     *20.,40.,60.,80.,100.,110.,120.,130.,140.,150.,160.,170.,180.,
     *200.,220.,240.,260.,280.,300.,350.,400.,450.,500.,550.,575.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  11*0.,39./
      DATA (CROSEC(I,2,7,7),I=1,51)/
     * 246.,308.,424.,590.,729.,776.,800.,821.,822.,817., 
     * 800.,778.,768.,728.,690.,654.,615.,584.,556.,480., 
     * 430.,393.,373.,367.,368.,370.,375.,388.,390.,397.,
     * 364.,337.,310.,291.,275.,268.,268.,268.,268.,  12*0./
      DATA (CROSEC(I,3,7,7),I=1,51)/
     * 155.,188.,256.,360.,456.,492.,512.,526.,526.,520.,
     * 504.,491.,475.,450.,425.,396.,376.,360.,340.,300.,
     * 282.,270.,265.,265.,266.,268.,273.,280.,288.,288.,
     * 256.,237.,226.,218.,208.,202.,202.,202.,202., 12*0./
C                                   K+
      DATA (CROSEC(I,1,7,8),I=1,51)/
     *20.,40.,60.,80.,100.,110.,120.,130.,140.,150.,160.,170.,180.,
     *200.,220.,240.,260.,280.,300.,350.,400.,450.,500.,550.,575.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  11*0.,39./
      DATA (CROSEC(I,2,7,8),I=1,51)/
     * 150.,212.,328.,500.,680.,735.,762.,781.,782.,779.,
     * 770.,748.,740.,706.,672.,633.,600.,569.,541.,467.,
     * 419.,385.,368.,364.,366.,368.,375.,388.,390.,397.,
     * 364.,337.,310.,291.,275.,268.,268.,268.,268.,  12*0./
      DATA (CROSEC(I,3,7,8),I=1,51)/
     *  90.,140.,208.,300.,426.,467.,490.,504.,504.,500.,
     * 484.,474.,460.,437.,413.,381.,365.,350.,330.,292.,
     * 276.,267.,263.,264.,265.,267.,273.,280.,288.,288.,
     * 256.,237.,226.,218.,208.,202.,202.,202.,202.,  12*0./
C                                   K0
      DATA (CROSEC(I,1,7,9),I=1,51)/
     *20.,40.,60.,80.,100.,110.,120.,130.,140.,150.,160.,170.,180.,
     *200.,220.,240.,260.,280.,300.,350.,400.,450.,500.,550.,575.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  11*0.,39./
      DATA (CROSEC(I,2,7,9),I=1,51)/
     * 150.,212.,328.,500.,680.,735.,762.,781.,782.,779.,
     * 770.,748.,740.,706.,672.,633.,600.,569.,541.,467.,
     * 419.,385.,368.,364.,366.,368.,375.,388.,390.,397.,
     * 364.,337.,310.,291.,275.,268.,268.,268.,268.,  12*0./
      DATA (CROSEC(I,3,7,9),I=1,51)/
     *  90.,140.,208.,300.,426.,467.,490.,504.,504.,500.,
     * 484.,474.,460.,437.,413.,381.,365.,350.,330.,292.,
     * 276.,267.,263.,264.,265.,267.,273.,280.,288.,288.,
     * 256.,237.,226.,218.,208.,202.,202.,202.,202.,  12*0./
C                             ANTI K0
      DATA (CROSEC(I,1,7,10),I=1,51)/
     *20.,40.,60.,80.,100.,110.,120.,130.,140.,150.,160.,170.,180.,
     *200.,220.,240.,260.,280.,300.,350.,400.,450.,500.,550.,575.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  11*0.,39./
      DATA (CROSEC(I,2,7,10),I=1,51)/
     * 246.,308.,424.,590.,729.,776.,800.,821.,822.,817., 
     * 800.,778.,768.,728.,690.,654.,615.,584.,556.,480., 
     * 430.,393.,373.,367.,368.,370.,375.,388.,390.,397.,
     * 364.,337.,310.,291.,275.,268.,268.,268.,268.,  12*0./
      DATA (CROSEC(I,3,7,10),I=1,51)/
     * 155.,188.,256.,360.,456.,492.,512.,526.,526.,520.,
     * 504.,491.,475.,450.,425.,396.,376.,360.,340.,300.,
     * 282.,270.,265.,265.,266.,268.,273.,280.,288.,288.,
     * 256.,237.,226.,218.,208.,202.,202.,202.,202., 12*0./
C 
C ********************** DATA FOR O(8,16), INUCL=8 ******************
C                                   NEUTRON
      DATA (CROSEC(I,1,8,1),I=1,51)/
     *14.,15.,17.,20.,22.,25.,27.,30.,35.,40.,45.,50.,60.,70.,80.,90.,
     *100.,120.,140.,150.,160.,180.,200.,250.,300.,350.,
     *400.,500.,600.,700.,800.,900.,1000.,1500.,2000.,3000.,
     *5000.,7000.,10000.,20000.,50000.,100000.,500000.,1000000.,
     *6*0.,44./
      DATA (CROSEC(I,2,8,1),I=1,51)/
     * 1520.,1570.,1630.,1660.,1647.,1623.,1595.,1555.,1475.,
     * 1395.,1290.,1207.,1035., 925., 816., 720., 645., 540., 
     *  462., 438., 415., 392., 378., 362., 361., 381., 390.,
     *  403., 417., 440., 460., 470., 479., 498., 504., 498., 
     *  477., 457., 443., 427., 420., 425., 429., 430.,  7*0./
      DATA (CROSEC(I,3,8,1),I=1,51)/
     * 750.,740.,700.,650.,620.,575.,555.,530.,505.,
     * 462.,435.,420.,375.,345.,320.,310.,300.,293.,
     * 288.,282.,282.,280.,276.,270.,271.,275.,280.,
     * 290.,295.,304.,310.,315.,318.,332.,335.,330.,
     * 323.,320.,317.,315.,315.,315.,315.,315.,  7*0./
C                                   PROTON
      DATA (CROSEC(I,1,8,2),I=1,51)/
     *14.,15.,17.,20.,22.,25.,27.,30.,35.,40.,45.,50.,60.,70.,80.,90.,
     *100.,120.,140.,150.,160.,180.,200.,250.,300.,350.,
     *400.,500.,600.,700.,800.,900.,1000.,1500.,2000.,3000.,
     *5000.,7000.,10000.,20000.,50000.,100000.,500000.,1000000.,
     *6*0.,44./
      DATA (CROSEC(I,2,8,2),I=1,51)/
     * 1520.,1570.,1630.,1660.,1647.,1623.,1595.,1555.,1475.,
     * 1395.,1290.,1207.,1035., 925., 816., 720., 645., 540., 
     *  462., 438., 415., 392., 378., 362., 361., 381., 390.,
     *  403., 417., 440., 460., 470., 479., 498., 504., 498., 
     *  477., 457., 443., 427., 420., 425., 429., 430.,  7*0./
      DATA (CROSEC(I,3,8,2),I=1,51)/
     * 460.,485.,510.,535.,537.,532.,520.,500.,460.,
     * 432.,405.,390.,350.,320.,310.,304.,293.,287.,
     * 283.,279.,279.,278.,276.,270.,271.,275.,280.,
     * 290.,295.,304.,310.,315.,318.,332.,335.,330.,
     * 323.,320.,317.,315.,315.,315.,315.,315.,  7*0./
C                                   PI-
      DATA (CROSEC(I,1,8,3),I=1,51)/
     *20.,40.,60.,80.,100.,120.,140.,160.,180.,
     *200.,220.,250.,300.,350.,400.,450.,500.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  19*0.,31./
      DATA (CROSEC(I,2,8,3),I=1,51)/
     * 280.,360.,500.,685.,812.,861.,870.,865.,835.,800.,
     * 755.,700.,600.,537.,493.,468.,441.,436.,443.,449.,
     * 460.,463.,432.,385.,350.,325.,312.,307.,303.,303.,
     * 303.,  20*0./
      DATA (CROSEC(I,3,8,3),I=1,51)/
     * 190.,207.,300.,420.,500.,540.,550.,542.,520.,490.,
     * 460.,423.,360.,339.,321.,314.,312.,314.,319.,324.,
     * 328.,330.,300.,275.,250.,240.,229.,225.,222.,222.,
     * 222.,  20*0./
C                                   PI+
      DATA (CROSEC(I,1,8,4),I=1,51)/
     *20.,40.,60.,80.,100.,120.,140.,160.,180.,
     *200.,220.,250.,300.,350.,400.,450.,500.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  19*0.,31./
      DATA (CROSEC(I,2,8,4),I=1,51)/
     * 170.,240.,390.,570.,740.,818.,830.,822.,800.,765.,
     * 725.,675.,585.,525.,483.,458.,444.,447.,453.,449.,
     * 460.,463.,432.,385.,350.,325.,312.,307.,303.,303.,
     * 303.,  20*0./
      DATA (CROSEC(I,3,8,4),I=1,51)/
     * 100.,145.,240.,340.,470.,518.,530.,522.,505.,477.,
     * 448.,412.,350.,330.,316.,310.,308.,311.,317.,324.,
     * 328.,330.,300.,275.,250.,240.,229.,225.,222.,222.,
     * 222.,  20*0./
C                                   ANTINEUTRON
      DATA (CROSEC(I,1,8,5),I=1,51)/
     *14.,15.,17.,20.,22.,25.,27.,30.,35.,40.,45.,50.,60.,70.,80.,90.,
     *100.,120.,140.,150.,160.,180.,200.,250.,300.,350.,
     *400.,500.,600.,700.,800.,900.,1000.,1500.,2000.,3000.,
     *5000.,7000.,10000.,20000.,50000.,100000.,500000.,1000000.,
     *6*0.,44./
      DATA (CROSEC(I,2,8,5),I=1,51)/
     * 1520.,1570.,1630.,1660.,1647.,1623.,1595.,1555.,1475.,
     * 1395.,1290.,1207.,1035., 925., 816., 720., 645., 540., 
     *  462., 438., 415., 392., 378., 362., 361., 381., 390.,
     *  403., 417., 440., 460., 470., 479., 498., 504., 498., 
     *  477., 457., 443., 427., 420., 425., 429., 430.,  7*0./
      DATA (CROSEC(I,3,8,5),I=1,51)/
     * 750.,740.,700.,650.,620.,575.,555.,530.,505.,
     * 462.,435.,420.,375.,345.,320.,310.,300.,293.,
     * 288.,282.,282.,280.,276.,270.,271.,275.,280.,
     * 290.,295.,304.,310.,315.,318.,332.,335.,330.,
     * 323.,320.,317.,315.,315.,315.,315.,315.,  7*0./
C                                   ANTIPROTON
      DATA (CROSEC(I,1,8,6),I=1,51)/
     *14.,15.,17.,20.,22.,25.,27.,30.,35.,40.,45.,50.,60.,70.,80.,90.,
     *100.,120.,140.,150.,160.,180.,200.,250.,300.,350.,
     *400.,500.,600.,700.,800.,900.,1000.,1500.,2000.,3000.,
     *5000.,7000.,10000.,20000.,50000.,100000.,500000.,1000000.,
     *6*0.,44./
      DATA (CROSEC(I,2,8,6),I=1,51)/
     * 1520.,1570.,1630.,1660.,1647.,1623.,1595.,1555.,1475.,
     * 1395.,1290.,1207.,1035., 925., 816., 720., 645., 540., 
     *  462., 438., 415., 392., 378., 362., 361., 381., 390.,
     *  403., 417., 440., 460., 470., 479., 498., 504., 498., 
     *  477., 457., 443., 427., 420., 425., 429., 430.,  7*0./
      DATA (CROSEC(I,3,8,6),I=1,51)/
     * 460.,485.,510.,535.,537.,532.,520.,500.,460.,
     * 432.,405.,390.,350.,320.,310.,304.,293.,287.,
     * 283.,279.,279.,278.,276.,270.,271.,275.,280.,
     * 290.,295.,304.,310.,315.,318.,332.,335.,330.,
     * 323.,320.,317.,315.,315.,315.,315.,315.,  7*0./
C                                   K-
      DATA (CROSEC(I,1,8,7),I=1,51)/
     *20.,40.,60.,80.,100.,120.,140.,160.,180.,
     *200.,220.,250.,300.,350.,400.,450.,500.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  19*0.,31./
      DATA (CROSEC(I,2,8,7),I=1,51)/
     * 280.,360.,500.,685.,812.,861.,870.,865.,835.,800.,
     * 755.,700.,600.,537.,493.,468.,441.,436.,443.,449.,
     * 460.,463.,432.,385.,350.,325.,312.,307.,303.,303.,
     * 303.,  20*0./
      DATA (CROSEC(I,3,8,7),I=1,51)/
     * 190.,207.,300.,420.,500.,540.,550.,542.,520.,490.,
     * 460.,423.,360.,339.,321.,314.,312.,314.,319.,324.,
     * 328.,330.,300.,275.,250.,240.,229.,225.,222.,222.,
     * 222.,  20*0./
C                                   K+
      DATA (CROSEC(I,1,8,8),I=1,51)/
     *20.,40.,60.,80.,100.,120.,140.,160.,180.,
     *200.,220.,250.,300.,350.,400.,450.,500.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  19*0.,31./
      DATA (CROSEC(I,2,8,8),I=1,51)/
     * 170.,240.,390.,570.,740.,818.,830.,822.,800.,765.,
     * 725.,675.,585.,525.,483.,458.,444.,447.,453.,449.,
     * 460.,463.,432.,385.,350.,325.,312.,307.,303.,303.,
     * 303.,  20*0./
      DATA (CROSEC(I,3,8,8),I=1,51)/
     * 100.,145.,240.,340.,470.,518.,530.,522.,505.,477.,
     * 448.,412.,350.,330.,316.,310.,308.,311.,317.,324.,
     * 328.,330.,300.,275.,250.,240.,229.,225.,222.,222.,
     * 222.,  20*0./
C                                   K0
      DATA (CROSEC(I,1,8,9),I=1,51)/
     *20.,40.,60.,80.,100.,120.,140.,160.,180.,
     *200.,220.,250.,300.,350.,400.,450.,500.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  19*0.,31./
      DATA (CROSEC(I,2,8,9),I=1,51)/
     * 170.,240.,390.,570.,740.,818.,830.,822.,800.,765.,
     * 725.,675.,585.,525.,483.,458.,444.,447.,453.,449.,
     * 460.,463.,432.,385.,350.,325.,312.,307.,303.,303.,
     * 303.,  20*0./
      DATA (CROSEC(I,3,8,9),I=1,51)/
     * 100.,145.,240.,340.,470.,518.,530.,522.,505.,477.,
     * 448.,412.,350.,330.,316.,310.,308.,311.,317.,324.,
     * 328.,330.,300.,275.,250.,240.,229.,225.,222.,222.,
     * 222.,  20*0./
C                          ANTI K0
      DATA (CROSEC(I,1,8,10),I=1,51)/
     *20.,40.,60.,80.,100.,120.,140.,160.,180.,
     *200.,220.,250.,300.,350.,400.,450.,500.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  19*0.,31./
      DATA (CROSEC(I,2,8,10),I=1,51)/
     * 280.,360.,500.,685.,812.,861.,870.,865.,835.,800.,
     * 755.,700.,600.,537.,493.,468.,441.,436.,443.,449.,
     * 460.,463.,432.,385.,350.,325.,312.,307.,303.,303.,
     * 303.,  20*0./
      DATA (CROSEC(I,3,8,10),I=1,51)/
     * 190.,207.,300.,420.,500.,540.,550.,542.,520.,490.,
     * 460.,423.,360.,339.,321.,314.,312.,314.,319.,324.,
     * 328.,330.,300.,275.,250.,240.,229.,225.,222.,222.,
     * 222.,  20*0./
C
C ******************** DATA FOR NA(11,23), INUCL=9 ***************
C                                   NEUTRON
      DATA (CROSEC(I,1,9,1),I=1,51)/
     *14.,15.,17.,20.,22.,25.,27.,30.,35.,40.,45.,50.,60.,70.,80.,90.,
     *100.,120.,140.,150.,160.,180.,200.,250.,300.,350.,
     *400.,500.,600.,700.,800.,900.,1000.,1500.,2000.,3000.,
     *5000.,7000.,10000.,20000.,50000.,100000.,500000.,1000000.,
     *6*0.,44./
      DATA (CROSEC(I,2,9,1),I=1,51)/
     * 1570.,1620.,1695.,1730.,1750.,1760.,1755.,1740.,1710.,
     * 1643.,1560.,1480.,1343.,1220.,1073., 953., 860., 720.,
     *  618., 582., 546., 522., 504., 484., 492., 500., 512.,
     *  538., 560., 586., 608., 622., 632., 660., 668., 664.,
     *  640., 616., 596., 568., 568., 568., 568., 568.,  7*0./
      DATA (CROSEC(I,3,9,1),I=1,51)/
     * 960.,930.,890.,822.,790.,750.,725.,686.,620.,
     * 600.,575.,540.,497.,450.,414.,390.,380.,372.,
     * 364.,360.,355.,354.,350.,350.,350.,356.,364.,
     * 384.,392.,400.,408.,410.,420.,408.,412.,420.,
     * 411.,409.,407.,403.,400.,400.,400.,400.,  7*0./
C                                   PROTON
      DATA (CROSEC(I,1,9,2),I=1,51)/
     *14.,15.,17.,20.,22.,25.,27.,30.,35.,40.,45.,50.,60.,70.,80.,90.,
     *100.,120.,140.,150.,160.,180.,200.,250.,300.,350.,
     *400.,500.,600.,700.,800.,900.,1000.,1500.,2000.,3000.,
     *5000.,7000.,10000.,20000.,50000.,100000.,500000.,1000000.,
     *6*0.,44./
      DATA (CROSEC(I,2,9,2),I=1,51)/
     * 1570.,1620.,1695.,1730.,1750.,1760.,1755.,1740.,1710.,
     * 1643.,1560.,1480.,1343.,1220.,1073., 953., 860., 720.,
     *  618., 582., 546., 522., 504., 484., 492., 500., 512.,
     *  538., 560., 586., 608., 622., 632., 660., 668., 664.,
     *  640., 616., 596., 568., 568., 568., 568., 568.,  7*0./
      DATA (CROSEC(I,3,9,2),I=1,51)/
     * 600.,617.,660.,675.,680.,680.,670.,650.,575.,
     * 550.,525.,490.,450.,420.,385.,367.,360.,350.,
     * 350.,350.,345.,347.,350.,350.,350.,356.,364.,
     * 384.,392.,400.,408.,410.,420.,408.,412.,420.,
     * 411.,409.,407.,403.,400.,400.,400.,400.,  7*0./
C                                   PI-
      DATA (CROSEC(I,1,9,3),I=1,51)/
     *20.,40.,60.,80.,100.,120.,140.,160.,180.,
     *200.,220.,250.,300.,350.,400.,450.,500.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  19*0.,31./
      DATA (CROSEC(I,2,9,3),I=1,51)/
     * 450.,545.,705.,910.,1020.,1075.,1087.,1080.,1042.,987.,
     * 943.,885.,790.,700., 650., 610., 585., 575., 585.,595.,
     * 600.,610.,556.,524., 494., 458., 445., 429., 427.,427.,
     * 427., 20*0./
      DATA (CROSEC(I,3,9,3),I=1,51)/
     * 275.,315.,413.,545.,620.,660.,670.,662.,630.,593.,
     * 570.,520.,465.,420.,410.,395.,390.,400.,410.,418.,
     * 420.,422.,372.,348.,330.,320.,310.,294.,292.,292.,
     * 292.,  20*0./
C                                   PI+
      DATA (CROSEC(I,1,9,4),I=1,51)/
     *20.,40.,60.,80.,100.,120.,140.,160.,180.,
     *200.,220.,250.,300.,350.,400.,450.,500.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  19*0.,31./
      DATA (CROSEC(I,2,9,4),I=1,51)/
     * 210.,320.,530.,795.,960.,1035.,1050.,1040.,1007.,957.,
     * 918.,865.,773.,685.,636., 598., 575., 565., 578.,590.,
     * 598.,610.,556.,524.,494., 458., 445., 429., 427.,427.,
     * 427.,  20*0./
      DATA (CROSEC(I,3,9,4),I=1,51)/
     * 115.,210.,340.,495.,585.,630.,645.,637.,605.,572.,
     * 550.,505.,455.,410.,401.,388.,383.,393.,405.,414.,
     * 418.,422.,372.,348.,330.,320.,310.,294.,292.,292.,
     * 292.,  20*0./
C                                   ANTINEUTRON
      DATA (CROSEC(I,1,9,5),I=1,51)/
     *14.,15.,17.,20.,22.,25.,27.,30.,35.,40.,45.,50.,60.,70.,80.,90.,
     *100.,120.,140.,150.,160.,180.,200.,250.,300.,350.,
     *400.,500.,600.,700.,800.,900.,1000.,1500.,2000.,3000.,
     *5000.,7000.,10000.,20000.,50000.,100000.,500000.,1000000.,
     *6*0.,44./
      DATA (CROSEC(I,2,9,5),I=1,51)/
     * 1570.,1620.,1695.,1730.,1750.,1760.,1755.,1740.,1710.,
     * 1643.,1560.,1480.,1343.,1220.,1073., 953., 860., 720.,
     *  618., 582., 546., 522., 504., 484., 492., 500., 512.,
     *  538., 560., 586., 608., 622., 632., 660., 668., 664.,
     *  640., 616., 596., 568., 568., 568., 568., 568.,  7*0./
      DATA (CROSEC(I,3,9,5),I=1,51)/
     * 960.,930.,890.,822.,790.,750.,725.,686.,620.,
     * 600.,575.,540.,497.,450.,414.,390.,380.,372.,
     * 364.,360.,355.,354.,350.,350.,350.,356.,364.,
     * 384.,392.,400.,408.,410.,420.,408.,412.,420.,
     * 411.,409.,407.,403.,400.,400.,400.,400.,  7*0./
C                                   ANTIPROTON
      DATA (CROSEC(I,1,9,6),I=1,51)/
     *14.,15.,17.,20.,22.,25.,27.,30.,35.,40.,45.,50.,60.,70.,80.,90.,
     *100.,120.,140.,150.,160.,180.,200.,250.,300.,350.,
     *400.,500.,600.,700.,800.,900.,1000.,1500.,2000.,3000.,
     *5000.,7000.,10000.,20000.,50000.,100000.,500000.,1000000.,
     *6*0.,44./
      DATA (CROSEC(I,2,9,6),I=1,51)/
     * 1570.,1620.,1695.,1730.,1750.,1760.,1755.,1740.,1710.,
     * 1643.,1560.,1480.,1343.,1220.,1073., 953., 860., 720.,
     *  618., 582., 546., 522., 504., 484., 492., 500., 512.,
     *  538., 560., 586., 608., 622., 632., 660., 668., 664.,
     *  640., 616., 596., 568., 568., 568., 568., 568.,  7*0./
      DATA (CROSEC(I,3,9,6),I=1,51)/
     * 600.,617.,660.,675.,680.,680.,670.,650.,575.,
     * 550.,525.,490.,450.,420.,385.,367.,360.,350.,
     * 350.,350.,345.,347.,350.,350.,350.,356.,364.,
     * 384.,392.,400.,408.,410.,420.,408.,412.,420.,
     * 411.,409.,407.,403.,400.,400.,400.,400.,  7*0./
C                                   K-
      DATA (CROSEC(I,1,9,7),I=1,51)/
     *20.,40.,60.,80.,100.,120.,140.,160.,180.,
     *200.,220.,250.,300.,350.,400.,450.,500.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  19*0.,31./
      DATA (CROSEC(I,2,9,7),I=1,51)/
     * 450.,545.,705.,910.,1020.,1075.,1087.,1080.,1042.,987.,
     * 943.,885.,790.,700., 650., 610., 585., 575., 585.,595.,
     * 600.,610.,556.,524., 494., 458., 445., 429., 427.,427.,
     * 427., 20*0./
      DATA (CROSEC(I,3,9,7),I=1,51)/
     * 275.,315.,413.,545.,620.,660.,670.,662.,630.,593.,
     * 570.,520.,465.,420.,410.,395.,390.,400.,410.,418.,
     * 420.,422.,372.,348.,330.,320.,310.,294.,292.,292.,
     * 292.,  20*0./
C                                   K+
      DATA (CROSEC(I,1,9,8),I=1,51)/
     *20.,40.,60.,80.,100.,120.,140.,160.,180.,
     *200.,220.,250.,300.,350.,400.,450.,500.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  19*0.,31./
      DATA (CROSEC(I,2,9,8),I=1,51)/
     * 210.,320.,530.,795.,960.,1035.,1050.,1040.,1007.,957.,
     * 918.,865.,773.,685.,636., 598., 575., 565., 578.,590.,
     * 598.,610.,556.,524.,494., 458., 445., 429., 427.,427.,
     * 427.,  20*0./
      DATA (CROSEC(I,3,9,8),I=1,51)/
     * 115.,210.,340.,495.,585.,630.,645.,637.,605.,572.,
     * 550.,505.,455.,410.,401.,388.,383.,393.,405.,414.,
     * 418.,422.,372.,348.,330.,320.,310.,294.,292.,292.,
     * 292.,  20*0./
C                                   K0
      DATA (CROSEC(I,1,9,9),I=1,51)/
     *20.,40.,60.,80.,100.,120.,140.,160.,180.,
     *200.,220.,250.,300.,350.,400.,450.,500.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  19*0.,31./
      DATA (CROSEC(I,2,9,9),I=1,51)/
     * 210.,320.,530.,795.,960.,1035.,1050.,1040.,1007.,957.,
     * 918.,865.,773.,685.,636., 598., 575., 565., 578.,590.,
     * 598.,610.,556.,524.,494., 458., 445., 429., 427.,427.,
     * 427.,  20*0./
      DATA (CROSEC(I,3,9,9),I=1,51)/
     * 115.,210.,340.,495.,585.,630.,645.,637.,605.,572.,
     * 550.,505.,455.,410.,401.,388.,383.,393.,405.,414.,
     * 418.,422.,372.,348.,330.,320.,310.,294.,292.,292.,
     * 292.,  20*0./
C                            ANTI K0
      DATA (CROSEC(I,1,9,10),I=1,51)/
     *20.,40.,60.,80.,100.,120.,140.,160.,180.,
     *200.,220.,250.,300.,350.,400.,450.,500.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  19*0.,31./
      DATA (CROSEC(I,2,9,10),I=1,51)/
     * 450.,545.,705.,910.,1020.,1075.,1087.,1080.,1042.,987.,
     * 943.,885.,790.,700., 650., 610., 585., 575., 585.,595.,
     * 600.,610.,556.,524., 494., 458., 445., 429., 427.,427.,
     * 427., 20*0./
      DATA (CROSEC(I,3,9,10),I=1,51)/
     * 275.,315.,413.,545.,620.,660.,670.,662.,630.,593.,
     * 570.,520.,465.,420.,410.,395.,390.,400.,410.,418.,
     * 420.,422.,372.,348.,330.,320.,310.,294.,292.,292.,
     * 292.,  20*0./
C       
C ******************** DATA FOR AL(13,27), INUCL=10 ******************
C                                    NEUTRON
      DATA (CROSEC(I,1,10,1),I=1,51)/
     *14.,15.,16.,17.,20.,22.,25.,27.,30.,35.,40.,45.,50.,60.,70.,80.,
     *90.,100.,120.,140.,150.,160.,180.,200.,250.,300.,350.,
     *400.,500.,600.,700.,800.,900.,1000.,1500.,2000.,3000.,
     *5000.,7000.,10000.,20000.,50000.,100000.,500000.,1000000.,
     *5*0.,45./
      DATA (CROSEC(I,2,10,1),I=1,51)/
     * 1735.,1750.,1760.,1795.,1830.,1855.,1885.,1896.,1900.,1870.,
     * 1835.,1785.,1710.,1522.,1350.,1212.,1080., 972., 816.,
     *  720., 678., 642., 600., 567., 558., 560., 578., 592.,
     *  616., 644., 672., 688., 708., 720., 736., 754., 736.,
     *  706., 680., 672., 646., 632., 632., 632., 632.,  6*0./
      DATA (CROSEC(I,3,10,1),I=1,51)/
     * 1000.,990.,975.,950.,905.,875.,825.,800.,762.,690.,
     *  652.,610.,570.,495.,480.,456.,444.,432.,420.,
     *  420.,420.,420.,410.,410.,400.,402.,404.,408.,
     *  424.,438.,448.,450.,454.,456.,472.,480.,466.,
     *  456.,452.,448.,444.,440.,440.,440.,440.,  6*0./ 
C                                    PROTON
      DATA (CROSEC(I,1,10,2),I=1,51)/
     *14.,15.,16.,17.,20.,22.,25.,27.,30.,35.,40.,45.,50.,60.,70.,80.,
     *90.,100.,120.,140.,150.,160.,180.,200.,250.,300.,350.,
     *400.,500.,600.,700.,800.,900.,1000.,1500.,2000.,3000.,
     *5000.,7000.,10000.,20000.,50000.,100000.,500000.,1000000.,
     *5*0.,45./
      DATA (CROSEC(I,2,10,2),I=1,51)/
     * 1735.,1750.,1760.,1795.,1830.,1855.,1885.,1896.,1900.,1870.,
     * 1835.,1785.,1710.,1522.,1350.,1212.,1080., 972., 816.,
     *  720., 678., 642., 600., 567., 558., 560., 578., 592.,
     *  616., 644., 672., 688., 708., 720., 736., 754., 736.,
     *  706., 680., 672., 646., 632., 632., 632., 632.,  6*0./
      DATA (CROSEC(I,3,10,2),I=1,51)/
     * 650.,682.,690.,715.,750.,762.,750.,740.,720.,655.,
     * 617.,575.,540.,470.,455.,432.,420.,408.,400.,
     * 403.,403.,408.,406.,404.,400.,402.,404.,408.,
     * 424.,438.,448.,450.,454.,456.,472.,480.,466.,
     * 456.,452.,448.,444.,440.,440.,440.,440.,  6*0./
C                                    PI-
      DATA (CROSEC(I,1,10,3),I=1,51)/
     *20.,40.,60.,80.,100.,120.,140.,160.,180.,
     *200.,220.,250.,300.,350.,400.,450.,500.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  19*0.,31./
      DATA (CROSEC(I,2,10,3),I=1,51)/
     * 532.,637.,832.,1057.,1207.,1230.,1210.,1174.,1333.,1095.,
     *1038.,970.,890., 807., 750., 710., 675., 665., 670., 673.,
     * 678.,682.,618., 574., 546., 520., 507., 495., 488., 488.,
     * 488.,  20*0./
      DATA (CROSEC(I,3,10,3),I=1,51)/
     * 300.,360.,495.,665.,750.,765.,750.,730.,700.,660.,
     * 615.,570.,520.,490.,470.,450.,448.,450.,450.,452.,
     * 456.,460.,408.,392.,376.,356.,347.,338.,332.,332.,
     * 332.,  20*0./
C                                    PI+
      DATA (CROSEC(I,1,10,4),I=1,51)/
     *20.,40.,60.,80.,100.,120.,140.,160.,180.,
     *200.,220.,250.,300.,350.,400.,450.,500.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  19*0.,31./
      DATA (CROSEC(I,2,10,4),I=1,51)/
     * 225.,350.,615.,945.,1122.,1175.,1157.,1128.,1088.,1045.,
     * 988.,935.,870.,787., 730., 690., 660., 652., 660., 668.,
     * 678.,682.,618., 574., 546., 520., 507., 495., 488., 488.,
     * 488.,  20*0./
      DATA (CROSEC(I,3,10,4),I=1,51)/
     * 120.,238.,390.,610.,712.,735.,720.,703.,655.,635.,
     * 590.,550.,505.,475.,455.,438.,440.,445.,445.,450.,
     * 456.,460.,408.,392.,376.,356.,347.,338.,332.,332.,
     * 332.,  20*0./
C                                    ANTINEUTRON
      DATA (CROSEC(I,1,10,5),I=1,51)/
     *14.,15.,16.,17.,20.,22.,25.,27.,30.,35.,40.,45.,50.,60.,70.,80.,
     *90.,100.,120.,140.,150.,160.,180.,200.,250.,300.,350.,
     *400.,500.,600.,700.,800.,900.,1000.,1500.,2000.,3000.,
     *5000.,7000.,10000.,20000.,50000.,100000.,500000.,1000000.,
     *5*0.,45./
      DATA (CROSEC(I,2,10,5),I=1,51)/
     * 1735.,1750.,1760.,1795.,1830.,1855.,1885.,1896.,1900.,1870.,
     * 1835.,1785.,1710.,1522.,1350.,1212.,1080., 972., 816.,
     *  720., 678., 642., 600., 567., 558., 560., 578., 592.,
     *  616., 644., 672., 688., 708., 720., 736., 754., 736.,
     *  706., 680., 672., 646., 632., 632., 632., 632.,  6*0./
      DATA (CROSEC(I,3,10,5),I=1,51)/
     * 1000.,990.,975.,950.,905.,875.,825.,800.,762.,690.,
     *  652.,610.,570.,495.,480.,456.,444.,432.,420.,
     *  420.,420.,420.,410.,410.,400.,402.,404.,408.,
     *  424.,438.,448.,450.,454.,456.,472.,480.,466.,
     *  456.,452.,448.,444.,440.,440.,440.,440.,  6*0./ 
C                                    ANTIPROTON
      DATA (CROSEC(I,1,10,6),I=1,51)/
     *14.,15.,16.,17.,20.,22.,25.,27.,30.,35.,40.,45.,50.,60.,70.,80.,
     *90.,100.,120.,140.,150.,160.,180.,200.,250.,300.,350.,
     *400.,500.,600.,700.,800.,900.,1000.,1500.,2000.,3000.,
     *5000.,7000.,10000.,20000.,50000.,100000.,500000.,1000000.,
     *5*0.,45./
      DATA (CROSEC(I,2,10,6),I=1,51)/
     * 1735.,1750.,1760.,1795.,1830.,1855.,1885.,1896.,1900.,1870.,
     * 1835.,1785.,1710.,1522.,1350.,1212.,1080., 972., 816.,
     *  720., 678., 642., 600., 567., 558., 560., 578., 592.,
     *  616., 644., 672., 688., 708., 720., 736., 754., 736.,
     *  706., 680., 672., 646., 632., 632., 632., 632.,  6*0./
      DATA (CROSEC(I,3,10,6),I=1,51)/
     * 650.,682.,690.,715.,750.,762.,750.,740.,720.,655.,
     * 617.,575.,540.,470.,455.,432.,420.,408.,400.,
     * 403.,403.,408.,406.,404.,400.,402.,404.,408.,
     * 424.,438.,448.,450.,454.,456.,472.,480.,466.,
     * 456.,452.,448.,444.,440.,440.,440.,440.,  6*0./
C                                    K-
      DATA (CROSEC(I,1,10,7),I=1,51)/
     *20.,40.,60.,80.,100.,120.,140.,160.,180.,
     *200.,220.,250.,300.,350.,400.,450.,500.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  19*0.,31./
      DATA (CROSEC(I,2,10,7),I=1,51)/
     * 532.,637.,832.,1057.,1207.,1230.,1210.,1174.,1333.,1095.,
     *1038.,970.,890., 807., 750., 710., 675., 665., 670., 673.,
     * 678.,682.,618., 574., 546., 520., 507., 495., 488., 488.,
     * 488.,  20*0./
      DATA (CROSEC(I,3,10,7),I=1,51)/
     * 300.,360.,495.,665.,750.,765.,750.,730.,700.,660.,
     * 615.,570.,520.,490.,470.,450.,448.,450.,450.,452.,
     * 456.,460.,408.,392.,376.,356.,347.,338.,332.,332.,
     * 332.,  20*0./
C                                    K+
      DATA (CROSEC(I,1,10,8),I=1,51)/
     *20.,40.,60.,80.,100.,120.,140.,160.,180.,
     *200.,220.,250.,300.,350.,400.,450.,500.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  19*0.,31./
      DATA (CROSEC(I,2,10,8),I=1,51)/
     * 225.,350.,615.,945.,1122.,1175.,1157.,1128.,1088.,1045.,
     * 988.,935.,870.,787., 730., 690., 660., 652., 660., 668.,
     * 678.,682.,618., 574., 546., 520., 507., 495., 488., 488.,
     * 488.,  20*0./
      DATA (CROSEC(I,3,10,8),I=1,51)/
     * 120.,238.,390.,610.,712.,735.,720.,703.,655.,635.,
     * 590.,550.,505.,475.,455.,438.,440.,445.,445.,450.,
     * 456.,460.,408.,392.,376.,356.,347.,338.,332.,332.,
     * 332.,  20*0./
C                                    K0
      DATA (CROSEC(I,1,10,9),I=1,51)/
     *20.,40.,60.,80.,100.,120.,140.,160.,180.,
     *200.,220.,250.,300.,350.,400.,450.,500.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  19*0.,31./
      DATA (CROSEC(I,2,10,9),I=1,51)/
     * 225.,350.,615.,945.,1122.,1175.,1157.,1128.,1088.,1045.,
     * 988.,935.,870.,787., 730., 690., 660., 652., 660., 668.,
     * 678.,682.,618., 574., 546., 520., 507., 495., 488., 488.,
     * 488.,  20*0./
      DATA (CROSEC(I,3,10,9),I=1,51)/
     * 120.,238.,390.,610.,712.,735.,720.,703.,655.,635.,
     * 590.,550.,505.,475.,455.,438.,440.,445.,445.,450.,
     * 456.,460.,408.,392.,376.,356.,347.,338.,332.,332.,
     * 332.,  20*0./
C                             ANTI K0
      DATA (CROSEC(I,1,10,10),I=1,51)/
     *20.,40.,60.,80.,100.,120.,140.,160.,180.,
     *200.,220.,250.,300.,350.,400.,450.,500.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  19*0.,31./
      DATA (CROSEC(I,2,10,10),I=1,51)/
     * 532.,637.,832.,1057.,1207.,1230.,1210.,1174.,1333.,1095.,
     *1038.,970.,890., 807., 750., 710., 675., 665., 670., 673.,
     * 678.,682.,618., 574., 546., 520., 507., 495., 488., 488.,
     * 488.,  20*0./
      DATA (CROSEC(I,3,10,10),I=1,51)/
     * 300.,360.,495.,665.,750.,765.,750.,730.,700.,660.,
     * 615.,570.,520.,490.,470.,450.,448.,450.,450.,452.,
     * 456.,460.,408.,392.,376.,356.,347.,338.,332.,332.,
     * 332.,  20*0./
C
C ********************* DATA FOR SI(14,24), INUCL=11 ******************
C                                    NEUTRON
      DATA (CROSEC(I,1,11,1),I=1,51)/
     *14.,15.,16.,17.,20.,22.,25.,27.,30.,35.,40.,45.,50.,60.,70.,80.,
     *90.,100.,120.,140.,150.,160.,180.,200.,250.,300.,350.,
     *400.,500.,600.,700.,800.,900.,1000.,1500.,2000.,3000.,
     *5000.,7000.,10000.,20000.,50000.,100000.,500000.,1000000.,
     *5*0.,45./
      DATA (CROSEC(I,2,11,1),I=1,51)/
     * 1810.,1833.,1850.,1872.,1920.,1950.,1995.,2020.,2035.,2000.,
     * 1930.,1850.,1760.,1570.,1400.,1255.,1110.,1008., 846.,
     *  742., 696., 671., 623., 588., 584., 584., 602., 618., 
     *  645., 679., 708., 727., 746., 757., 769., 782., 771.,
     *  734., 710., 698., 672., 654., 650., 650., 650.,  6*0./
      DATA (CROSEC(I,3,11,1),I=1,51)/
     * 1060.,1035.,1015.,990.,935.,900.,860.,830.,790.,725.,
     *  665., 630., 600.,520.,504.,486.,470.,456.,444.,
     *  432., 432., 432.,418.,418.,415.,412.,416.,422.,
     *  440., 460., 472.,476.,479.,480.,492.,496.,488.,
     *  472., 472., 464.,460.,452.,448.,448.,448.,  6*0./
C                                    PROTON
      DATA (CROSEC(I,1,11,2),I=1,51)/
     *14.,15.,16.,17.,20.,22.,25.,27.,30.,35.,40.,45.,50.,60.,70.,80.,
     *90.,100.,120.,140.,150.,160.,180.,200.,250.,300.,350.,
     *400.,500.,600.,700.,800.,900.,1000.,1500.,2000.,3000.,
     *5000.,7000.,10000.,20000.,50000.,100000.,500000.,1000000.,
     *5*0.,45./
      DATA (CROSEC(I,2,11,2),I=1,51)/
     * 1810.,1833.,1850.,1872.,1920.,1950.,1995.,2020.,2035.,2000.,
     * 1930.,1850.,1760.,1570.,1400.,1255.,1110.,1008., 846.,
     *  742., 696., 671., 623., 588., 584., 584., 602., 618., 
     *  645., 679., 708., 727., 746., 757., 769., 782., 771.,
     *  734., 710., 698., 672., 654., 650., 650., 650.,  6*0./
      DATA (CROSEC(I,3,11,2),I=1,51)/
     * 670.,700.,725.,750.,780.,780.,770.,757.,735.,690.,
     * 635.,585.,570.,490.,475.,460.,446.,431.,423.,
     * 425.,425.,425.,425.,422.,422.,412.,416.,422.,
     * 440.,460.,472.,476.,479.,480.,492.,496.,488.,
     * 472.,472.,464.,460.,452.,448.,448.,448.,  6*0./
C                                    PI-
      DATA (CROSEC(I,1,11,3),I=1,51)/
     *20.,40.,60.,80.,100.,120.,140.,160.,180.,
     *200.,220.,250.,300.,350.,400.,450.,500.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  19*0.,31./
      DATA (CROSEC(I,2,11,3),I=1,51)/
     * 545.,652.,852.,1082.,1236.,1260.,1239.,1202.,1365.,1121.,
     *1063.,993.,911., 826., 768., 727., 691., 681., 686., 689.,
     * 694.,698.,633., 588., 559., 532., 519., 507., 499., 499.,
     * 499.,  20*0./
      DATA (CROSEC(I,3,11,3),I=1,51)/
     * 307.,368.,507.,681.,768.,783.,768.,747.,717.,676.,
     * 630.,583.,532.,502.,481.,461.,458.,461.,461.,463.,
     * 467.,471.,418.,401.,385.,364.,355.,346.,340.,340.,
     * 340.,  20*0./
C                                    PI+
      DATA (CROSEC(I,1,11,4),I=1,51)/
     *20.,40.,60.,80.,100.,120.,140.,160.,180.,
     *200.,220.,250.,300.,350.,400.,450.,500.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  19*0.,31./
      DATA (CROSEC(I,2,11,4),I=1,51)/
     * 230.,358.,630.,968.,1149.,1203.,1185.,1155.,1114.,1070.,
     *1012.,957.,891.,806., 747., 706., 676., 668., 676., 684.,
     * 694.,698.,633., 588., 559., 532., 519., 507., 499., 499.,
     * 499.,  20*0./
      DATA (CROSEC(I,3,11,4),I=1,51)/
     * 122.,243.,399.,624.,729.,753.,737.,720.,671.,650.,
     * 604.,563.,517.,486.,466.,448.,450.,455.,455.,461.,
     * 467.,471.,418.,401.,385.,364.,355.,346.,340.,340.,
     * 340.,  20*0./
C                                    ANTINEUTRON
      DATA (CROSEC(I,1,11,5),I=1,51)/
     *14.,15.,16.,17.,20.,22.,25.,27.,30.,35.,40.,45.,50.,60.,70.,80.,
     *90.,100.,120.,140.,150.,160.,180.,200.,250.,300.,350.,
     *400.,500.,600.,700.,800.,900.,1000.,1500.,2000.,3000.,
     *5000.,7000.,10000.,20000.,50000.,100000.,500000.,1000000.,
     *5*0.,45./
      DATA (CROSEC(I,2,11,5),I=1,51)/
     * 1810.,1833.,1850.,1872.,1920.,1950.,1995.,2020.,2035.,2000.,
     * 1930.,1850.,1760.,1570.,1400.,1255.,1110.,1008., 846.,
     *  742., 696., 671., 623., 588., 584., 584., 602., 618., 
     *  645., 679., 708., 727., 746., 757., 769., 782., 771.,
     *  734., 710., 698., 672., 654., 650., 650., 650.,  6*0./
      DATA (CROSEC(I,3,11,5),I=1,51)/
     * 1060.,1035.,1015.,990.,935.,900.,860.,830.,790.,725.,
     *  665., 630., 600.,520.,504.,486.,470.,456.,444.,
     *  432., 432., 432.,418.,418.,415.,412.,416.,422.,
     *  440., 460., 472.,476.,479.,480.,492.,496.,488.,
     *  472., 472., 464.,460.,452.,448.,448.,448.,  6*0./
C                                    ANTIPROTON
      DATA (CROSEC(I,1,11,6),I=1,51)/
     *14.,15.,16.,17.,20.,22.,25.,27.,30.,35.,40.,45.,50.,60.,70.,80.,
     *90.,100.,120.,140.,150.,160.,180.,200.,250.,300.,350.,
     *400.,500.,600.,700.,800.,900.,1000.,1500.,2000.,3000.,
     *5000.,7000.,10000.,20000.,50000.,100000.,500000.,1000000.,
     *5*0.,45./
      DATA (CROSEC(I,2,11,6),I=1,51)/
     * 1810.,1833.,1850.,1872.,1920.,1950.,1995.,2020.,2035.,2000.,
     * 1930.,1850.,1760.,1570.,1400.,1255.,1110.,1008., 846.,
     *  742., 696., 671., 623., 588., 584., 584., 602., 618., 
     *  645., 679., 708., 727., 746., 757., 769., 782., 771.,
     *  734., 710., 698., 672., 654., 650., 650., 650.,  6*0./
      DATA (CROSEC(I,3,11,6),I=1,51)/
     * 670.,700.,725.,750.,780.,780.,770.,757.,735.,690.,
     * 635.,585.,570.,490.,475.,460.,446.,431.,423.,
     * 425.,425.,425.,425.,422.,422.,412.,416.,422.,
     * 440.,460.,472.,476.,479.,480.,492.,496.,488.,
     * 472.,472.,464.,460.,452.,448.,448.,448.,  6*0./
C                                    K-
      DATA (CROSEC(I,1,11,7),I=1,51)/
     *20.,40.,60.,80.,100.,120.,140.,160.,180.,
     *200.,220.,250.,300.,350.,400.,450.,500.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  19*0.,31./
      DATA (CROSEC(I,2,11,7),I=1,51)/
     * 545.,652.,852.,1082.,1236.,1260.,1239.,1202.,1365.,1121.,
     *1063.,993.,911., 826., 768., 727., 691., 681., 686., 689.,
     * 694.,698.,633., 588., 559., 532., 519., 507., 499., 499.,
     * 499.,  20*0./
      DATA (CROSEC(I,3,11,7),I=1,51)/
     * 307.,368.,507.,681.,768.,783.,768.,747.,717.,676.,
     * 630.,583.,532.,502.,481.,461.,458.,461.,461.,463.,
     * 467.,471.,418.,401.,385.,364.,355.,346.,340.,340.,
     * 340.,  20*0./
C                                    K+
      DATA (CROSEC(I,1,11,8),I=1,51)/
     *20.,40.,60.,80.,100.,120.,140.,160.,180.,
     *200.,220.,250.,300.,350.,400.,450.,500.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  19*0.,31./
      DATA (CROSEC(I,2,11,8),I=1,51)/
     * 230.,358.,630.,968.,1149.,1203.,1185.,1155.,1114.,1070.,
     *1012.,957.,891.,806., 747., 706., 676., 668., 676., 684.,
     * 694.,698.,633., 588., 559., 532., 519., 507., 499., 499.,
     * 499.,  20*0./
      DATA (CROSEC(I,3,11,8),I=1,51)/
     * 122.,243.,399.,624.,729.,753.,737.,720.,671.,650.,
     * 604.,563.,517.,486.,466.,448.,450.,455.,455.,461.,
     * 467.,471.,418.,401.,385.,364.,355.,346.,340.,340.,
     * 340.,  20*0./
C                                    K0
      DATA (CROSEC(I,1,11,9),I=1,51)/
     *20.,40.,60.,80.,100.,120.,140.,160.,180.,
     *200.,220.,250.,300.,350.,400.,450.,500.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  19*0.,31./
      DATA (CROSEC(I,2,11,9),I=1,51)/
     * 230.,358.,630.,968.,1149.,1203.,1185.,1155.,1114.,1070.,
     *1012.,957.,891.,806., 747., 706., 676., 668., 676., 684.,
     * 694.,698.,633., 588., 559., 532., 519., 507., 499., 499.,
     * 499.,  20*0./
      DATA (CROSEC(I,3,11,9),I=1,51)/
     * 122.,243.,399.,624.,729.,753.,737.,720.,671.,650.,
     * 604.,563.,517.,486.,466.,448.,450.,455.,455.,461.,
     * 467.,471.,418.,401.,385.,364.,355.,346.,340.,340.,
     * 340.,  20*0./
C                           ANTI K0
      DATA (CROSEC(I,1,11,10),I=1,51)/
     *20.,40.,60.,80.,100.,120.,140.,160.,180.,
     *200.,220.,250.,300.,350.,400.,450.,500.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  19*0.,31./
      DATA (CROSEC(I,2,11,10),I=1,51)/
     * 545.,652.,852.,1082.,1236.,1260.,1239.,1202.,1365.,1121.,
     *1063.,993.,911., 826., 768., 727., 691., 681., 686., 689.,
     * 694.,698.,633., 588., 559., 532., 519., 507., 499., 499.,
     * 499.,  20*0./
      DATA (CROSEC(I,3,11,10),I=1,51)/
     * 307.,368.,507.,681.,768.,783.,768.,747.,717.,676.,
     * 630.,583.,532.,502.,481.,461.,458.,461.,461.,463.,
     * 467.,471.,418.,401.,385.,364.,355.,346.,340.,340.,
     * 340.,  20*0./
C
C ****************** DATA FOR CA(20,40), INUCL=12 ********************
C                                    NEUTRON
      DATA (CROSEC(I,1,12,1),I=1,51)/
     *14.,15.,16.,17.,20.,22.,25.,27.,30.,35.,40.,45.,50.,60.,70.,80.,
     *90.,100.,120.,140.,150.,160.,180.,200.,250.,300.,350.,
     *400.,500.,600.,700.,800.,900.,1000.,1500.,2000.,3000.,
     *5000.,7000.,10000.,20000.,50000.,100000.,500000.,1000000.,
     *5*0.,45./
      DATA (CROSEC(I,2,12,1),I=1,51)/
     * 2180.,2130.,2095.,2075.,2115.,2150.,2220.,2250.,2300.,2365.,
     * 2360.,2280.,2180.,2000.,1805.,1650.,1500.,1340.,1140.,
     *  990., 940., 890., 825., 790., 770., 773., 787., 800.,
     *  830., 870., 905., 930., 950., 965., 990.,1002., 990.,
     *  965., 945., 925., 892., 860., 860., 860., 860.,  6*0./
      DATA (CROSEC(I,3,12,1),I=1,51)/
     * 1240.,1225.,1200.,1180.,1125.,1090.,1045.,1020.,980.,925.,
     *  880., 825., 770., 680., 640., 620., 615.,600.,580.,
     *  565., 560., 560., 560., 550., 535., 530.,540.,550.,
     *  570., 595., 610., 615., 620., 622., 629.,630.,620.,
     *  612., 607., 592., 587., 580., 580., 580.,580., 6*0./
C                                    PROTON
      DATA (CROSEC(I,1,12,2),I=1,51)/
     *14.,15.,16.,17.,20.,22.,25.,27.,30.,35.,40.,45.,50.,60.,70.,80.,
     *90.,100.,120.,140.,150.,160.,180.,200.,250.,300.,350.,
     *400.,500.,600.,700.,800.,900.,1000.,1500.,2000.,3000.,
     *5000.,7000.,10000.,20000.,50000.,100000.,500000.,1000000.,
     *5*0.,45./
      DATA (CROSEC(I,2,12,2),I=1,51)/
     * 2180.,2130.,2095.,2075.,2115.,2150.,2220.,2250.,2300.,2365.,
     * 2360.,2280.,2180.,2000.,1805.,1650.,1500.,1340.,1140.,
     *  990., 940., 890., 825., 790., 770., 773., 787., 800.,
     *  830., 870., 905., 930., 950., 965., 990.,1002., 990.,
     *  965., 945., 925., 892., 860., 860., 860., 860.,  6*0./
      DATA (CROSEC(I,3,12,2),I=1,51)/
     * 770.,800.,823.,850.,900.,925.,935.,920.,895.,835.,
     * 800.,750.,715.,640.,605.,590.,588.,573.,555.,
     * 543.,540.,540.,540.,535.,530.,530.,540.,550.,
     * 570.,595.,610.,615.,620.,622.,629.,630.,620.,
     * 612.,607.,592.,587.,580.,580.,580.,580.,  6*0./
C                                    PI-
      DATA (CROSEC(I,1,12,3),I=1,51)/
     *20.,40.,60.,80.,100.,120.,140.,160.,180.,
     *200.,220.,250.,300.,350.,400.,450.,500.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  19*0.,31./
      DATA (CROSEC(I,2,12,3),I=1,51)/
     * 800., 980.,1240.,1460.,1570.,1600.,1580.,1535.,1475.,1425.,
     *1375.,1295.,1200.,1083.,1000., 948., 915., 895., 900., 908.,
     * 915., 922., 855., 795., 740., 705., 682., 660., 660., 660., 
     * 660.,  20*0./
      DATA (CROSEC(I,3,12,3),I=1,51)/
     * 470.,550.,620.,860.,955.,980.,960.,920.,860.,820.,
     * 780.,740.,665.,637.,615.,600.,590.,590.,600.,608.,
     * 610.,615.,550.,525.,510.,488.,470.,450.,450.,450.,
     * 450.,  20*0./
C                                    PI+
      DATA (CROSEC(I,1,12,4),I=1,51)/
     *20.,40.,60.,80.,100.,120.,140.,160.,180.,
     *200.,220.,250.,300.,350.,400.,450.,500.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  19*0.,31./
      DATA (CROSEC(I,2,12,4),I=1,51)/
     * 275., 445., 790.,1195.,1440.,1485.,1475.,1435.,1385.,1335.,  ! DKFZc
     *1295.,1245.,1160.,1050., 970., 923.,895., 877., 887., 897., 
     * 904., 913., 855., 795., 740., 705.,682., 660., 660., 660., 
     * 660.,  20*0./
      DATA (CROSEC(I,3,12,4),I=1,51)/
     * 160.,315.,500.,745.,870.,905.,900.,860.,810.,770.,
     * 740.,710.,640.,617.,595.,585.,575.,575.,590.,600.,
     * 602.,608.,550.,525.,510.,488.,470.,450.,450.,450.,
     * 450.,  20*0./
C                                    ANTINEUTRON
      DATA (CROSEC(I,1,12,5),I=1,51)/
     *14.,15.,16.,17.,20.,22.,25.,27.,30.,35.,40.,45.,50.,60.,70.,80.,
     *90.,100.,120.,140.,150.,160.,180.,200.,250.,300.,350.,
     *400.,500.,600.,700.,800.,900.,1000.,1500.,2000.,3000.,
     *5000.,7000.,10000.,20000.,50000.,100000.,500000.,1000000.,
     *5*0.,45./
      DATA (CROSEC(I,2,12,5),I=1,51)/
     * 2180.,2130.,2095.,2075.,2115.,2150.,2220.,2250.,2300.,2365.,
     * 2360.,2280.,2180.,2000.,1805.,1650.,1500.,1340.,1140.,
     *  990., 940., 890., 825., 790., 770., 773., 787., 800.,
     *  830., 870., 905., 930., 950., 965., 990.,1002., 990.,
     *  965., 945., 925., 892., 860., 860., 860., 860.,  6*0./
      DATA (CROSEC(I,3,12,5),I=1,51)/
     * 1240.,1225.,1200.,1180.,1125.,1090.,1045.,1020.,980.,925.,
     *  880., 825., 770., 680., 640., 620., 615.,600.,580.,
     *  565., 560., 560., 560., 550., 535., 530.,540.,550.,
     *  570., 595., 610., 615., 620., 622., 629.,630.,620.,
     *  612., 607., 592., 587., 580., 580., 580.,580., 6*0./
C                                    ANTIPROTON
      DATA (CROSEC(I,1,12,6),I=1,51)/
     *14.,15.,16.,17.,20.,22.,25.,27.,30.,35.,40.,45.,50.,60.,70.,80.,
     *90.,100.,120.,140.,150.,160.,180.,200.,250.,300.,350.,
     *400.,500.,600.,700.,800.,900.,1000.,1500.,2000.,3000.,
     *5000.,7000.,10000.,20000.,50000.,100000.,500000.,1000000.,
     *5*0.,45./
      DATA (CROSEC(I,2,12,6),I=1,51)/
     * 2180.,2130.,2095.,2075.,2115.,2150.,2220.,2250.,2300.,2365.,
     * 2360.,2280.,2180.,2000.,1805.,1650.,1500.,1340.,1140.,
     *  990., 940., 890., 825., 790., 770., 773., 787., 800.,
     *  830., 870., 905., 930., 950., 965., 990.,1002., 990.,
     *  965., 945., 925., 892., 860., 860., 860., 860.,  6*0./
      DATA (CROSEC(I,3,12,6),I=1,51)/
     * 770.,800.,823.,850.,900.,925.,935.,920.,895.,835.,
     * 800.,750.,715.,640.,605.,590.,588.,573.,555.,
     * 543.,540.,540.,540.,535.,530.,530.,540.,550.,
     * 570.,595.,610.,615.,620.,622.,629.,630.,620.,
     * 612.,607.,592.,587.,580.,580.,580.,580.,  6*0./
C                                    K-
      DATA (CROSEC(I,1,12,7),I=1,51)/
     *20.,40.,60.,80.,100.,120.,140.,160.,180.,
     *200.,220.,250.,300.,350.,400.,450.,500.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  19*0.,31./
      DATA (CROSEC(I,2,12,7),I=1,51)/
     * 800., 980.,1240.,1460.,1570.,1600.,1580.,1535.,1475.,1425.,
     *1375.,1295.,1200.,1083.,1000., 948., 915., 895., 900., 908.,
     * 915., 922., 855., 795., 740., 705., 682., 660., 660., 660., 
     * 660.,  20*0./
      DATA (CROSEC(I,3,12,7),I=1,51)/
     * 470.,550.,620.,860.,955.,980.,960.,920.,860.,820.,
     * 780.,740.,665.,637.,615.,600.,590.,590.,600.,608.,
     * 610.,615.,550.,525.,510.,488.,470.,450.,450.,450.,
     * 450.,  20*0./
C                                    K+
      DATA (CROSEC(I,1,12,8),I=1,51)/
     *20.,40.,60.,80.,100.,120.,140.,160.,180.,
     *200.,220.,250.,300.,350.,400.,450.,500.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  19*0.,31./
      DATA (CROSEC(I,2,12,8),I=1,51)/
     * 275., 445., 790.,1195.,1440.,1485.,475.,1435.,1385.,1335.,
     *1295.,1245.,1160.,1050., 970., 923.,895., 877., 887., 897., 
     * 904., 913., 855., 795., 740., 705.,682., 660., 660., 660., 
     * 660.,  20*0./
      DATA (CROSEC(I,3,12,8),I=1,51)/
     * 160.,315.,500.,745.,870.,905.,900.,860.,810.,770.,
     * 740.,710.,640.,617.,595.,585.,575.,575.,590.,600.,
     * 602.,608.,550.,525.,510.,488.,470.,450.,450.,450.,
     * 450.,  20*0./
C                                    K0
      DATA (CROSEC(I,1,12,9),I=1,51)/
     *20.,40.,60.,80.,100.,120.,140.,160.,180.,
     *200.,220.,250.,300.,350.,400.,450.,500.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  19*0.,31./
      DATA (CROSEC(I,2,12,9),I=1,51)/
     * 275., 445., 790.,1195.,1440.,1485.,475.,1435.,1385.,1335.,
     *1295.,1245.,1160.,1050., 970., 923.,895., 877., 887., 897., 
     * 904., 913., 855., 795., 740., 705.,682., 660., 660., 660., 
     * 660.,  20*0./
      DATA (CROSEC(I,3,12,9),I=1,51)/
     * 160.,315.,500.,745.,870.,905.,900.,860.,810.,770.,
     * 740.,710.,640.,617.,595.,585.,575.,575.,590.,600.,
     * 602.,608.,550.,525.,510.,488.,470.,450.,450.,450.,
     * 450.,  20*0./
C                            ANTI K0
      DATA (CROSEC(I,1,12,10),I=1,51)/
     *20.,40.,60.,80.,100.,120.,140.,160.,180.,
     *200.,220.,250.,300.,350.,400.,450.,500.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  19*0.,31./
      DATA (CROSEC(I,2,12,10),I=1,51)/
     * 800., 980.,1240.,1460.,1570.,1600.,1580.,1535.,1475.,1425.,
     *1375.,1295.,1200.,1083.,1000., 948., 915., 895., 900., 908.,
     * 915., 922., 855., 795., 740., 705., 682., 660., 660., 660., 
     * 660.,  20*0./
      DATA (CROSEC(I,3,12,10),I=1,51)/
     * 470.,550.,620.,860.,955.,980.,960.,920.,860.,820.,
     * 780.,740.,665.,637.,615.,600.,590.,590.,600.,608.,
     * 610.,615.,550.,525.,510.,488.,470.,450.,450.,450.,
     * 450.,  20*0./
C
C ******************** DATA FOR FE(26,56), INUCL=13 *****************
C                                    NEUTRON
      DATA (CROSEC(I,1,13,1),I=1,51)/
     *14.,15.,17.,20.,22.,25.,27.,30.,33.,35.,37.,40.,45.,50.,55.,60.,
     *70.,80.,90.,100.,120.,140.,150.,160.,180.,200.,250.,300.,350.,
     *400.,500.,600.,700.,800.,900.,1000.,1500.,2000.,3000.,
     *5000.,7000.,10000.,20000.,50000.,100000.,500000.,1000000.,
     *3*0.,47./
      DATA (CROSEC(I,2,13,1),I=1,51)/
     * 2580.,2490.,2370.,2282.,2275.,2285.,2320.,2370.,2423.,2445.,
     * 2460.,2485.,2530.,2540.,2517.,2480.,2290.,2110.,1940.,1790.,
     * 1510.,1290.,1220.,1150.,1070.,1030.,1013.,1020.,1030.,1043.,
     * 1075.,1110.,1133.,1163.,1185.,1225.,1252.,1260.,1260.,
     * 1233.,1207.,1185.,1140.,1110.,1110.,1110.,1110.,  4*0./
      DATA (CROSEC(I,3,13,1),I=1,51)/
     * 1440.,1433.,1390.,1325.,1280.,1260.,1215.,1180.,1140.,1110.,
     * 1080.,1040., 990., 955., 920., 885., 835., 800., 780., 765.,
     *  750., 725., 720., 720., 710., 700., 700., 700., 712., 705.,
     *  735., 750., 765., 775., 780., 795., 810., 813., 810.,
     *  784., 757., 743., 735., 720., 720., 720., 720.,  4*0./
C                                    PROTON
      DATA (CROSEC(I,1,13,2),I=1,51)/
     *14.,15.,17.,20.,22.,25.,27.,30.,33.,35.,37.,40.,45.,50.,55.,60.,
     *70.,80.,90.,100.,120.,140.,150.,160.,180.,200.,250.,300.,350.,
     *400.,500.,600.,700.,800.,900.,1000.,1500.,2000.,3000.,
     *5000.,7000.,10000.,20000.,50000.,100000.,500000.,1000000.,
     *3*0.,47./
      DATA (CROSEC(I,2,13,2),I=1,51)/
     * 2580.,2490.,2370.,2282.,2275.,2285.,2320.,2370.,2423.,2445.,
     * 2460.,2485.,2530.,2540.,2517.,2480.,2290.,2110.,1940.,1790.,
     * 1510.,1290.,1220.,1150.,1070.,1030.,1013.,1020.,1030.,1043.,
     * 1075.,1110.,1133.,1163.,1185.,1225.,1252.,1260.,1260.,
     * 1233.,1207.,1185.,1140.,1110.,1110.,1110.,1110.,  4*0./
      DATA (CROSEC(I,3,13,2),I=1,51)/
     * 900.,960.,1070.,1090.,1115.,1120.,1115.,1080.,1045.,1025.,
     *1000.,960., 900., 885., 865., 790., 765., 740., 720., 700.,
     * 697.,697., 697., 697., 695., 690., 688., 690., 712., 705.,
     * 735.,750., 765., 775., 780., 795., 810., 813., 810.,
     * 784.,757., 743., 735., 720., 720., 720., 720.,  4*0./
C                                    PI-
      DATA (CROSEC(I,1,13,3),I=1,51)/
     *20.,40.,60.,80.,100.,120.,140.,160.,180.,
     *200.,220.,250.,300.,350.,400.,450.,500.,550.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  18*0.,32./
      DATA (CROSEC(I,2,13,3),I=1,51)/
     * 1175.,1363.,1670.,1950.,2050.,2040.,1975.,1886.,1834.,1773.,
     * 1720.,1635.,1474.,1380.,1269.,1225.,1182.,1162.,1159.,1162., 
     * 1178.,1190.,1197.,1102.,1035., 975., 945., 925., 905., 905.,
     *  905., 905.,  19*0./
      DATA (CROSEC(I,3,13,3),I=1,51)/
     * 625.,725.,910.,1180.,1275.,1250.,1200.,1150.,1100.,1040.,
     * 995.,925.,825., 810., 780., 760., 745., 740., 740., 740.,
     * 750.,760.,765., 690., 660., 635., 615., 600., 585., 585.,
     * 585.,585.,  19*0./
C                                    PI+
      DATA (CROSEC(I,1,13,4),I=1,51)/
     *20.,40.,60.,80.,100.,120.,140.,160.,180.,
     *200.,220.,250.,300.,350.,400.,450.,500.,550.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  18*0.,32./
      DATA (CROSEC(I,2,13,4),I=1,51)/
     * 330., 575.,1010.,1500.,1837.,1875.,1820.,1751.,1691.,1636.,
     *1590.,1450.,1396.,1305.,1219.,1190.,1148.,1138.,1134.,1144.,
     *1163.,1175.,1183.,1098.,1035., 975., 945., 925., 905., 905.,
     * 905., 905.,  19*0./
      DATA (CROSEC(I,3,13,4),I=1,51)/
     * 210.,410.,707.,1010.,1125.,1150.,1100.,1070.,1010.,960.,
     * 920.,775.,780., 760., 750., 740., 720., 725., 725.,730.,
     * 740.,750.,755., 690., 660., 635., 615., 600., 585., 585.,
     * 585.,585.,  19*0./
C                                    ANTINEUTRON
      DATA (CROSEC(I,1,13,5),I=1,51)/
     *14.,15.,17.,20.,22.,25.,27.,30.,33.,35.,37.,40.,45.,50.,55.,60.,
     *70.,80.,90.,100.,120.,140.,150.,160.,180.,200.,250.,300.,350.,
     *400.,500.,600.,700.,800.,900.,1000.,1500.,2000.,3000.,
     *5000.,7000.,10000.,20000.,50000.,100000.,500000.,1000000.,
     *3*0.,47./
      DATA (CROSEC(I,2,13,5),I=1,51)/
     * 2580.,2490.,2370.,2282.,2275.,2285.,2320.,2370.,2423.,2445.,
     * 2460.,2485.,2530.,2540.,2517.,2480.,2290.,2110.,1940.,1790.,
     * 1510.,1290.,1220.,1150.,1070.,1030.,1013.,1020.,1030.,1043.,
     * 1075.,1110.,1133.,1163.,1185.,1225.,1252.,1260.,1260.,
     * 1233.,1207.,1185.,1140.,1110.,1110.,1110.,1110.,  4*0./
      DATA (CROSEC(I,3,13,5),I=1,51)/
     * 1440.,1433.,1390.,1325.,1280.,1260.,1215.,1180.,1140.,1110.,
     * 1080.,1040., 990., 955., 920., 885., 835., 800., 780., 765.,
     *  750., 725., 720., 720., 710., 700., 700., 700., 712., 705.,
     *  735., 750., 765., 775., 780., 795., 810., 813., 810.,
     *  784., 757., 743., 735., 720., 720., 720., 720.,  4*0./
C                                    ANTIPROTON
      DATA (CROSEC(I,1,13,6),I=1,51)/
     *14.,15.,17.,20.,22.,25.,27.,30.,33.,35.,37.,40.,45.,50.,55.,60.,
     *70.,80.,90.,100.,120.,140.,150.,160.,180.,200.,250.,300.,350.,
     *400.,500.,600.,700.,800.,900.,1000.,1500.,2000.,3000.,
     *5000.,7000.,10000.,20000.,50000.,100000.,500000.,1000000.,
     *3*0.,47./
      DATA (CROSEC(I,2,13,6),I=1,51)/
     * 2580.,2490.,2370.,2282.,2275.,2285.,2320.,2370.,2423.,2445.,
     * 2460.,2485.,2530.,2540.,2517.,2480.,2290.,2110.,1940.,1790.,
     * 1510.,1290.,1220.,1150.,1070.,1030.,1013.,1020.,1030.,1043.,
     * 1075.,1110.,1133.,1163.,1185.,1225.,1252.,1260.,1260.,
     * 1233.,1207.,1185.,1140.,1110.,1110.,1110.,1110.,  4*0./
      DATA (CROSEC(I,3,13,6),I=1,51)/
     * 900.,960.,1070.,1090.,1115.,1120.,1115.,1080.,1045.,1025.,
     *1000.,960., 900., 885., 865., 790., 765., 740., 720., 700.,
     * 697.,697., 697., 697., 695., 690., 688., 690., 712., 705.,
     * 735.,750., 765., 775., 780., 795., 810., 813., 810.,
     * 784.,757., 743., 735., 720., 720., 720., 720.,  4*0./
C                                    K-
      DATA (CROSEC(I,1,13,7),I=1,51)/
     *20.,40.,60.,80.,100.,120.,140.,160.,180.,
     *200.,220.,250.,300.,350.,400.,450.,500.,550.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  18*0.,32./
      DATA (CROSEC(I,2,13,7),I=1,51)/
     * 1175.,1363.,1670.,1950.,2050.,2040.,1975.,1886.,1834.,1773.,
     * 1720.,1635.,1474.,1380.,1269.,1225.,1182.,1162.,1159.,1162., 
     * 1178.,1190.,1197.,1102.,1035., 975., 945., 925., 905., 905.,
     *  905., 905.,  19*0./
      DATA (CROSEC(I,3,13,7),I=1,51)/
     * 625.,725.,910.,1180.,1275.,1250.,1200.,1150.,1100.,1040.,
     * 995.,925.,825., 810., 780., 760., 745., 740., 740., 740.,
     * 750.,760.,765., 690., 660., 635., 615., 600., 585., 585.,
     * 585.,585.,  19*0./
C                                    K+
      DATA (CROSEC(I,1,13,8),I=1,51)/
     *20.,40.,60.,80.,100.,120.,140.,160.,180.,
     *200.,220.,250.,300.,350.,400.,450.,500.,550.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  18*0.,32./
      DATA (CROSEC(I,2,13,8),I=1,51)/
     * 330., 575.,1010.,1500.,1837.,1875.,1820.,1751.,1691.,1636.,
     *1590.,1450.,1396.,1305.,1219.,1190.,1148.,1138.,1134.,1144.,
     *1163.,1175.,1183.,1098.,1035., 975., 945., 925., 905., 905.,
     * 905., 905.,  19*0./
      DATA (CROSEC(I,3,13,8),I=1,51)/
     * 210.,410.,707.,1010.,1125.,1150.,1100.,1070.,1010.,960.,
     * 920.,775.,780., 760., 750., 740., 720., 725., 725.,730.,
     * 740.,750.,755., 690., 660., 635., 615., 600., 585., 585.,
     * 585.,585.,  19*0./
C                                    K0
      DATA (CROSEC(I,1,13,9),I=1,51)/
     *20.,40.,60.,80.,100.,120.,140.,160.,180.,
     *200.,220.,250.,300.,350.,400.,450.,500.,550.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  18*0.,32./
      DATA (CROSEC(I,2,13,9),I=1,51)/
     * 330., 575.,1010.,1500.,1837.,1875.,1820.,1751.,1691.,1636.,
     *1590.,1450.,1396.,1305.,1219.,1190.,1148.,1138.,1134.,1144.,
     *1163.,1175.,1183.,1098.,1035., 975., 945., 925., 905., 905.,
     * 905., 905.,  19*0./
      DATA (CROSEC(I,3,13,9),I=1,51)/
     * 210.,410.,707.,1010.,1125.,1150.,1100.,1070.,1010.,960.,
     * 920.,775.,780., 760., 750., 740., 720., 725., 725.,730.,
     * 740.,750.,755., 690., 660., 635., 615., 600., 585., 585.,
     * 585.,585.,  19*0./
C                            ANTI K0
      DATA (CROSEC(I,1,13,10),I=1,51)/
     *20.,40.,60.,80.,100.,120.,140.,160.,180.,
     *200.,220.,250.,300.,350.,400.,450.,500.,550.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  18*0.,32./
      DATA (CROSEC(I,2,13,10),I=1,51)/
     * 1175.,1363.,1670.,1950.,2050.,2040.,1975.,1886.,1834.,1773.,
     * 1720.,1635.,1474.,1380.,1269.,1225.,1182.,1162.,1159.,1162., 
     * 1178.,1190.,1197.,1102.,1035., 975., 945., 925., 905., 905.,
     *  905., 905.,  19*0./
      DATA (CROSEC(I,3,13,10),I=1,51)/
     * 625.,725.,910.,1180.,1275.,1250.,1200.,1150.,1100.,1040.,
     * 995.,925.,825., 810., 780., 760., 745., 740., 740., 740.,
     * 750.,760.,765., 690., 660., 635., 615., 600., 585., 585.,
     * 585.,585.,  19*0./
C
C ********************* DATA FOR CU(29,64), INUCL=14 *******************
C                                    NEUTRON
      DATA (CROSEC(I,1,14,1),I=1,51)/
     *14.,15.,17.,20.,22.,25.,27.,30.,33.,35.,37.,40.,45.,50.,55.,60.,
     *70.,80.,90.,100.,120.,140.,150.,160.,180.,200.,250.,300.,350.,
     *400.,500.,600.,700.,800.,900.,1000.,1500.,2000.,3000.,
     *5000.,7000.,10000.,20000.,50000.,100000.,500000.,1000000.,
     *3*0.,47./
      DATA (CROSEC(I,2,14,1),I=1,51)/
     * 2920.,2800.,2615.,2480.,2455.,2430.,2440.,2460.,2500.,2530.,
     * 2560.,2615.,2690.,2720.,2700.,2645.,2500.,2320.,2140.,1970.,
     * 1670.,1460.,1380.,1285.,1200.,1160.,1140.,1147.,1163.,1170.,
     * 1200.,1237.,1265.,1285.,1305.,1328.,1375.,1390.,1395.,
     * 1370.,1335.,1315.,1270.,1230.,1230.,1230.,1230.,  4*0./
      DATA (CROSEC(I,3,14,1),I=1,51)/
     * 1540.,1535.,1500.,1445.,1407.,1380.,1330.,1300.,1285.,1270.,
     * 1240.,1190.,1090.,1010., 940., 920., 860., 835., 820., 810., 
     *  800., 780., 775., 770., 760., 760., 758., 765., 765., 770.,
     *  795., 810., 825., 830., 840., 848., 870., 870., 868.,
     *  840., 825., 810., 803., 795., 795., 795., 795.,  4*0./
C                                    PROTON
      DATA (CROSEC(I,1,14,2),I=1,51)/
     *14.,15.,17.,20.,22.,25.,27.,30.,33.,35.,37.,40.,45.,50.,55.,60.,
     *70.,80.,90.,100.,120.,140.,150.,160.,180.,200.,250.,300.,350.,
     *400.,500.,600.,700.,800.,900.,1000.,1500.,2000.,3000.,
     *5000.,7000.,10000.,20000.,50000.,100000.,500000.,1000000.,
     *3*0.,47./
      DATA (CROSEC(I,2,14,2),I=1,51)/
     * 2920.,2800.,2615.,2480.,2455.,2430.,2440.,2460.,2500.,2530.,
     * 2560.,2615.,2690.,2720.,2700.,2645.,2500.,2320.,2140.,1970.,
     * 1670.,1460.,1380.,1285.,1200.,1160.,1140.,1147.,1163.,1170.,
     * 1200.,1237.,1265.,1285.,1305.,1328.,1375.,1390.,1395.,
     * 1370.,1335.,1315.,1270.,1230.,1230.,1230.,1230.,  4*0./
      DATA (CROSEC(I,3,14,2),I=1,51)/
     * 935.,1000.,1060.,1190.,1220.,1250.,1240.,1210.,1150.,1130.,  
     *1115.,1050., 985., 950., 890., 870., 820., 800., 785., 780., 
     * 770., 750., 745., 740., 735., 735., 745., 760., 762., 770.,
     * 795., 810., 825., 830., 840., 848., 870., 870., 868.,
     * 840., 825., 810., 803., 795., 795., 795., 795.,  4*0./
C                                    PI-
      DATA (CROSEC(I,1,14,3),I=1,51)/
     *20.,40.,60.,80.,100.,120.,140.,160.,180.,
     *200.,220.,250.,300.,350.,400.,450.,500.,550.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  18*0.,32./
      DATA (CROSEC(I,2,14,3),I=1,51)/
     * 1400.,1600.,1875.,2088.,2200.,2220.,2175.,2125.,2075.,2012.,
     * 1950.,1855.,1670.,1530.,1430.,1370.,1315.,1315.,1315.,1330.,
     * 1345.,1360.,1365.,1250.,1185.,1128.,1070.,1035.,1010.,1010.,
     * 1010.,1010.,  19*0./
      DATA (CROSEC(I,3,14,3),I=1,51)/
     *  725., 840.,1020.,1200.,1295.,1300.,1267.,1240.,1213.,1175.,
     * 1125.,1042., 950., 900., 860., 840., 830., 832., 835., 840.,
     *  850., 860., 865., 785., 735., 705., 680., 650., 630., 630., 
     *  630., 630.,  19*0./
C                                    PI+
      DATA (CROSEC(I,1,14,4),I=1,51)/
     *20.,40.,60.,80.,100.,120.,140.,160.,180.,
     *200.,220.,250.,300.,350.,400.,450.,500.,550.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  18*0.,32./
      DATA (CROSEC(I,2,14,4),I=1,51)/
     *  355., 605.,1120.,1630.,1940.,2010.,2010.,1980.,1925.,1895.,
     * 1830.,1730.,1585.,1490.,1400.,1340.,1290.,1290.,1290.,1310.,
     * 1330.,1345.,1350.,1240.,1185.,1128.,1070.,1035.,1010.,1010.,
     * 1010.,1010.,  19*0./
      DATA (CROSEC(I,3,14,4),I=1,51)/
     * 230., 425.,780.,1025.,1155.,1190.,1190.,1180.,1125.,1100.,
     *1050.,1000.,900., 870., 835., 815., 810., 812., 815., 825.,
     * 840.,850., 855., 780., 735., 705., 680., 650., 630., 630., 
     * 630., 630.,  19*0./
C                                    ANTINEUTRON
      DATA (CROSEC(I,1,14,5),I=1,51)/
     *14.,15.,17.,20.,22.,25.,27.,30.,33.,35.,37.,40.,45.,50.,55.,60.,
     *70.,80.,90.,100.,120.,140.,150.,160.,180.,200.,250.,300.,350.,
     *400.,500.,600.,700.,800.,900.,1000.,1500.,2000.,3000.,
     *5000.,7000.,10000.,20000.,50000.,100000.,500000.,1000000.,
     *3*0.,47./
      DATA (CROSEC(I,2,14,5),I=1,51)/
     * 2920.,2800.,2615.,2480.,2455.,2430.,2440.,2460.,2500.,2530.,
     * 2560.,2615.,2690.,2720.,2700.,2645.,2500.,2320.,2140.,1970.,
     * 1670.,1460.,1380.,1285.,1200.,1160.,1140.,1147.,1163.,1170.,
     * 1200.,1237.,1265.,1285.,1305.,1328.,1375.,1390.,1395.,
     * 1370.,1335.,1315.,1270.,1230.,1230.,1230.,1230.,  4*0./
      DATA (CROSEC(I,3,14,5),I=1,51)/
     * 1540.,1535.,1500.,1445.,1407.,1380.,1330.,1300.,1285.,1270.,
     * 1240.,1190.,1090.,1010., 940., 920., 860., 835., 820., 810., 
     *  800., 780., 775., 770., 760., 760., 758., 765., 765., 770.,
     *  795., 810., 825., 830., 840., 848., 870., 870., 868.,
     *  840., 825., 810., 803., 795., 795., 795., 795.,  4*0./
C                                    ANTIPROTON
      DATA (CROSEC(I,1,14,6),I=1,51)/
     *14.,15.,17.,20.,22.,25.,27.,30.,33.,35.,37.,40.,45.,50.,55.,60.,
     *70.,80.,90.,100.,120.,140.,150.,160.,180.,200.,250.,300.,350.,
     *400.,500.,600.,700.,800.,900.,1000.,1500.,2000.,3000.,
     *5000.,7000.,10000.,20000.,50000.,100000.,500000.,1000000.,
     *3*0.,47./
      DATA (CROSEC(I,2,14,6),I=1,51)/
     * 2920.,2800.,2615.,2480.,2455.,2430.,2440.,2460.,2500.,2530.,
     * 2560.,2615.,2690.,2720.,2700.,2645.,2500.,2320.,2140.,1970.,
     * 1670.,1460.,1380.,1285.,1200.,1160.,1140.,1147.,1163.,1170.,
     * 1200.,1237.,1265.,1285.,1305.,1328.,1375.,1390.,1395.,
     * 1370.,1335.,1315.,1270.,1230.,1230.,1230.,1230.,  4*0./
      DATA (CROSEC(I,3,14,6),I=1,51)/
     * 935.,1000.,1060.,1190.,1220.,1250.,1240.,1210.,1150.,1130.,  
     *1115.,1050., 985., 950., 890., 870., 820., 800., 785., 780., 
     * 770., 750., 745., 740., 735., 735., 745., 760., 762., 770.,
     * 795., 810., 825., 830., 840., 848., 870., 870., 868.,
     * 840., 825., 810., 803., 795., 795., 795., 795.,  4*0./
C                                    K-
      DATA (CROSEC(I,1,14,7),I=1,51)/
     *20.,40.,60.,80.,100.,120.,140.,160.,180.,
     *200.,220.,250.,300.,350.,400.,450.,500.,550.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  18*0.,32./
      DATA (CROSEC(I,2,14,7),I=1,51)/
     * 1400.,1600.,1875.,2088.,2200.,2220.,2175.,2125.,2075.,2012.,
     * 1950.,1855.,1670.,1530.,1430.,1370.,1315.,1315.,1315.,1330.,
     * 1345.,1360.,1365.,1250.,1185.,1128.,1070.,1035.,1010.,1010.,
     * 1010.,1010.,  19*0./
      DATA (CROSEC(I,3,14,7),I=1,51)/
     *  725., 840.,1020.,1200.,1295.,1300.,1267.,1240.,1213.,1175.,
     * 1125.,1042., 950., 900., 860., 840., 830., 832., 835., 840.,
     *  850., 860., 865., 785., 735., 705., 680., 650., 630., 630., 
     *  630., 630.,  19*0./
C                                    K+
      DATA (CROSEC(I,1,14,8),I=1,51)/
     *20.,40.,60.,80.,100.,120.,140.,160.,180.,
     *200.,220.,250.,300.,350.,400.,450.,500.,550.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  18*0.,32./
      DATA (CROSEC(I,2,14,8),I=1,51)/
     *  355., 605.,1120.,1630.,1940.,2010.,2010.,1980.,1925.,1895.,
     * 1830.,1730.,1585.,1490.,1400.,1340.,1290.,1290.,1290.,1310.,
     * 1330.,1345.,1350.,1240.,1185.,1128.,1070.,1035.,1010.,1010.,
     * 1010.,1010.,  19*0./
      DATA (CROSEC(I,3,14,8),I=1,51)/
     * 230., 425.,780.,1025.,1155.,1190.,1190.,1180.,1125.,1100.,
     *1050.,1000.,900., 870., 835., 815., 810., 812., 815., 825.,
     * 840.,850., 855., 780., 735., 705., 680., 650., 630., 630., 
     * 630., 630.,  19*0./
C                                    K0
      DATA (CROSEC(I,1,14,9),I=1,51)/
     *20.,40.,60.,80.,100.,120.,140.,160.,180.,
     *200.,220.,250.,300.,350.,400.,450.,500.,550.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  18*0.,32./
      DATA (CROSEC(I,2,14,9),I=1,51)/
     *  355., 605.,1120.,1630.,1940.,2010.,2010.,1980.,1925.,1895.,
     * 1830.,1730.,1585.,1490.,1400.,1340.,1290.,1290.,1290.,1310.,
     * 1330.,1345.,1350.,1240.,1185.,1128.,1070.,1035.,1010.,1010.,
     * 1010.,1010.,  19*0./
      DATA (CROSEC(I,3,14,9),I=1,51)/
     * 230., 425.,780.,1025.,1155.,1190.,1190.,1180.,1125.,1100.,
     *1050.,1000.,900., 870., 835., 815., 810., 812., 815., 825.,
     * 840.,850., 855., 780., 735., 705., 680., 650., 630., 630., 
     * 630., 630.,  19*0./
C                            ANTI K0
      DATA (CROSEC(I,1,14,10),I=1,51)/
     *20.,40.,60.,80.,100.,120.,140.,160.,180.,
     *200.,220.,250.,300.,350.,400.,450.,500.,550.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  18*0.,32./
      DATA (CROSEC(I,2,14,10),I=1,51)/
     * 1400.,1600.,1875.,2088.,2200.,2220.,2175.,2125.,2075.,2012.,
     * 1950.,1855.,1670.,1530.,1430.,1370.,1315.,1315.,1315.,1330.,
     * 1345.,1360.,1365.,1250.,1185.,1128.,1070.,1035.,1010.,1010.,
     * 1010.,1010.,  19*0./
      DATA (CROSEC(I,3,14,10),I=1,51)/
     *  725., 840.,1020.,1200.,1295.,1300.,1267.,1240.,1213.,1175.,
     * 1125.,1042., 950., 900., 860., 840., 830., 832., 835., 840.,
     *  850., 860., 865., 785., 735., 705., 680., 650., 630., 630., 
     *  630., 630.,  19*0./
C
C ********************* DATA FOR MO(42,96), INUCL=15 *********************
C                                    NEUTRON
      DATA (CROSEC(I,1,15,1),I=1,51)/
     *14.,15.,17.,20.,22.,25.,27.,30.,33.,35.,37.,40.,45.,50.,55.,60.,
     *70.,80.,90.,100.,120.,140.,150.,160.,180.,200.,250.,300.,350.,
     *400.,500.,600.,700.,800.,900.,1000.,1500.,2000.,3000.,
     *5000.,7000.,10000.,20000.,50000.,100000.,500000.,1000000.,
     *3*0.,47./
      DATA (CROSEC(I,2,15,1),I=1,51)/
     * 4150.,4040.,3800.,3490.,3300.,3060.,2960.,2845.,2785.,2820.,
     * 2850.,2980.,3170.,3230.,3270.,3280.,3225.,3075.,2895.,2710.,
     * 2355.,2060.,1925.,1800.,1630.,1560.,1540.,1550.,1570.,1590.,
     * 1650.,1685.,1715.,1740.,1760.,1780.,1850.,1880.,1858.,
     * 1815.,1790.,1782.,1720.,1690.,1690.,1690.,1690.,  4*0./
      DATA (CROSEC(I,3,15,1),I=1,51)/
     * 1790.,1775.,1740.,1680.,1640.,1580.,1550.,1510.,1460.,1440.,
     * 1418.,1380.,1330.,1280.,1240.,1200.,1155.,1140.,1110.,1110.,
     * 1080.,1065.,1050.,1050.,1025.,1020.,1015.,1020.,1022.,1026.,
     * 1060.,1085.,1100.,1110.,1120.,1127.,1150.,1160.,1140.,
     * 1100.,1085.,1080.,1070.,1070.,1070.,1070.,1070.,  4*0./
C                                    PROTON
      DATA (CROSEC(I,1,15,2),I=1,51)/
     *14.,15.,17.,20.,22.,25.,27.,30.,33.,35.,37.,40.,45.,50.,55.,60.,
     *70.,80.,90.,100.,120.,140.,150.,160.,180.,200.,250.,300.,350.,
     *400.,500.,600.,700.,800.,900.,1000.,1500.,2000.,3000.,
     *5000.,7000.,10000.,20000.,50000.,100000.,500000.,1000000.,
     *3*0.,47./
      DATA (CROSEC(I,2,15,2),I=1,51)/
     * 4150.,4040.,3800.,3490.,3300.,3060.,2960.,2845.,2785.,2820.,
     * 2850.,2980.,3170.,3230.,3270.,3280.,3225.,3075.,2895.,2710.,
     * 2355.,2060.,1925.,1800.,1630.,1560.,1540.,1550.,1570.,1590.,
     * 1650.,1685.,1715.,1740.,1760.,1780.,1850.,1880.,1858.,
     * 1815.,1790.,1782.,1720.,1690.,1690.,1690.,1690.,  4*0./
      DATA (CROSEC(I,3,15,2),I=1,51)/
     * 1025.,1080.,1190.,1380.,1440.,1495.,1475.,1420.,1350.,1310.,
     * 1300.,1290.,1250.,1200.,1170.,1130.,1095.,1060.,1040.,1022.,
     * 1020.,1016.,1016.,1016.,1016.,1012.,1005.,1005.,1005.,1010.,
     * 1060.,1085.,1100.,1110.,1120.,1127.,1150.,1160.,1140.,
     * 1100.,1085.,1080.,1070.,1070.,1070.,1070.,1070.,  4*0./
C                                    PI-
      DATA (CROSEC(I,1,15,3),I=1,51)/
     *20.,40.,50.,60.,70.,80.,90.,100.,120.,140.,160.,180.,
     *200.,220.,250.,300.,350.,400.,450.,500.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  16*0.,34./
      DATA (CROSEC(I,2,15,3),I=1,51)/
     * 2430.,2610.,2710.,2790.,2880.,2940.,2965.,2970.,2970.,2920.,
     * 2840.,2720.,2570.,2500.,2365.,2200.,2050.,1926.,1825.,1768.,
     * 1749.,1750.,1778.,1789.,1808.,1690.,1645.,1530.,1492.,1450.,
     * 1425.,1425.,1425.,1425.,  17*0./
      DATA (CROSEC(I,3,15,3),I=1,51)/
     *  925.,1125.,1250.,1375.,1500.,1600.,1680.,1750.,1770.,1730.,
     * 1660.,1580.,1500.,1450.,1330.,1250.,1190.,1140.,1100.,1075.,
     * 1075.,1070.,1088.,1095.,1110.,1035.,1005., 940., 917., 880.,
     *  860., 860., 860., 860.,  17*0./
C                                    PI+
      DATA (CROSEC(I,1,15,4),I=1,51)/
     *20.,40.,50.,60.,70.,80.,90.,100.,120.,140.,160.,180.,
     *200.,220.,250.,300.,350.,400.,450.,500.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  16*0.,34./
      DATA (CROSEC(I,2,15,4),I=1,51)/
     *  410., 730.,1110.,1530.,1920.,2200.,2385.,2520.,2600.,2630.,
     * 2575.,2470.,2320.,2285.,2185.,2053.,1945.,1852.,1776.,1719.,
     * 1710.,1716.,1746.,1759.,1778.,1675.,1645.,1530.,1492.,1450.,
     * 1425.,1425.,1425.,1425.,  17*0./
      DATA (CROSEC(I,3,15,4),I=1,51)/
     *  271., 540., 825., 975.,1140.,1285.,1400.,1480.,1555.,1580.,
     * 1525.,1470.,1360.,1340.,1255.,1160.,1120.,1085.,1060.,1045.,
     * 1045.,1045.,1065.,1075.,1090.,1025.,1005., 940., 917., 880.,
     *  860., 860., 860., 860.,  17*0./
C                                    ANTINEUTRON
      DATA (CROSEC(I,1,15,5),I=1,51)/
     *14.,15.,17.,20.,22.,25.,27.,30.,33.,35.,37.,40.,45.,50.,55.,60.,
     *70.,80.,90.,100.,120.,140.,150.,160.,180.,200.,250.,300.,350.,
     *400.,500.,600.,700.,800.,900.,1000.,1500.,2000.,3000.,
     *5000.,7000.,10000.,20000.,50000.,100000.,500000.,1000000.,
     *3*0.,47./
      DATA (CROSEC(I,2,15,5),I=1,51)/
     * 4150.,4040.,3800.,3490.,3300.,3060.,2960.,2845.,2785.,2820.,
     * 2850.,2980.,3170.,3230.,3270.,3280.,3225.,3075.,2895.,2710.,
     * 2355.,2060.,1925.,1800.,1630.,1560.,1540.,1550.,1570.,1590.,
     * 1650.,1685.,1715.,1740.,1760.,1780.,1850.,1880.,1858.,
     * 1815.,1790.,1782.,1720.,1690.,1690.,1690.,1690.,  4*0./
      DATA (CROSEC(I,3,15,5),I=1,51)/
     * 1790.,1775.,1740.,1680.,1640.,1580.,1550.,1510.,1460.,1440.,
     * 1418.,1380.,1330.,1280.,1240.,1200.,1155.,1140.,1110.,1110.,
     * 1080.,1065.,1050.,1050.,1025.,1020.,1015.,1020.,1022.,1026.,
     * 1060.,1085.,1100.,1110.,1120.,1127.,1150.,1160.,1140.,
     * 1100.,1085.,1080.,1070.,1070.,1070.,1070.,1070.,  4*0./
C                                    ANTIPROTON
      DATA (CROSEC(I,1,15,6),I=1,51)/
     *14.,15.,17.,20.,22.,25.,27.,30.,33.,35.,37.,40.,45.,50.,55.,60.,
     *70.,80.,90.,100.,120.,140.,150.,160.,180.,200.,250.,300.,350.,
     *400.,500.,600.,700.,800.,900.,1000.,1500.,2000.,3000.,
     *5000.,7000.,10000.,20000.,50000.,100000.,500000.,1000000.,
     *3*0.,47./
      DATA (CROSEC(I,2,15,6),I=1,51)/
     * 4150.,4040.,3800.,3490.,3300.,3060.,2960.,2845.,2785.,2820.,
     * 2850.,2980.,3170.,3230.,3270.,3280.,3225.,3075.,2895.,2710.,
     * 2355.,2060.,1925.,1800.,1630.,1560.,1540.,1550.,1570.,1590.,
     * 1650.,1685.,1715.,1740.,1760.,1780.,1850.,1880.,1858.,
     * 1815.,1790.,1782.,1720.,1690.,1690.,1690.,1690.,  4*0./
      DATA (CROSEC(I,3,15,6),I=1,51)/
     * 1025.,1080.,1190.,1380.,1440.,1495.,1475.,1420.,1350.,1310.,
     * 1300.,1290.,1250.,1200.,1170.,1130.,1095.,1060.,1040.,1022.,
     * 1020.,1016.,1016.,1016.,1016.,1012.,1005.,1005.,1005.,1010.,
     * 1060.,1085.,1100.,1110.,1120.,1127.,1150.,1160.,1140.,
     * 1100.,1085.,1080.,1070.,1070.,1070.,1070.,1070.,  4*0./
C                                    K-
      DATA (CROSEC(I,1,15,7),I=1,51)/
     *20.,40.,50.,60.,70.,80.,90.,100.,120.,140.,160.,180.,
     *200.,220.,250.,300.,350.,400.,450.,500.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  16*0.,34./
      DATA (CROSEC(I,2,15,7),I=1,51)/
     * 2430.,2610.,2710.,2790.,2880.,2940.,2965.,2970.,2970.,2920.,
     * 2840.,2720.,2570.,2500.,2365.,2200.,2050.,1926.,1825.,1768.,
     * 1749.,1750.,1778.,1789.,1808.,1690.,1645.,1530.,1492.,1450.,
     * 1425.,1425.,1425.,1425.,  17*0./
      DATA (CROSEC(I,3,15,7),I=1,51)/
     *  925.,1125.,1250.,1375.,1500.,1600.,1680.,1750.,1770.,1730.,
     * 1660.,1580.,1500.,1450.,1330.,1250.,1190.,1140.,1100.,1075.,
     * 1075.,1070.,1088.,1095.,1110.,1035.,1005., 940., 917., 880.,
     *  860., 860., 860., 860.,  17*0./
C                                    K+
      DATA (CROSEC(I,1,15,8),I=1,51)/
     *20.,40.,50.,60.,70.,80.,90.,100.,120.,140.,160.,180.,
     *200.,220.,250.,300.,350.,400.,450.,500.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  16*0.,34./
      DATA (CROSEC(I,2,15,8),I=1,51)/
     *  410., 730.,1110.,1530.,1920.,2200.,2385.,2520.,2600.,2630.,
     * 2575.,2470.,2320.,2285.,2185.,2053.,1945.,1852.,1776.,1719.,
     * 1710.,1716.,1746.,1759.,1778.,1675.,1645.,1530.,1492.,1450.,
     * 1425.,1425.,1425.,1425.,  17*0./
      DATA (CROSEC(I,3,15,8),I=1,51)/
     *  271., 540., 825., 975.,1140.,1285.,1400.,1480.,1555.,1580.,
     * 1525.,1470.,1360.,1340.,1255.,1160.,1120.,1085.,1060.,1045.,
     * 1045.,1045.,1065.,1075.,1090.,1025.,1005., 940., 917., 880.,
     *  860., 860., 860., 860.,  17*0./
C                                    K0
      DATA (CROSEC(I,1,15,9),I=1,51)/
     *20.,40.,50.,60.,70.,80.,90.,100.,120.,140.,160.,180.,
     *200.,220.,250.,300.,350.,400.,450.,500.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  16*0.,34./
      DATA (CROSEC(I,2,15,9),I=1,51)/
     *  410., 730.,1110.,1530.,1920.,2200.,2385.,2520.,2600.,2630.,
     * 2575.,2470.,2320.,2285.,2185.,2053.,1945.,1852.,1776.,1719.,
     * 1710.,1716.,1746.,1759.,1778.,1675.,1645.,1530.,1492.,1450.,
     * 1425.,1425.,1425.,1425.,  17*0./
      DATA (CROSEC(I,3,15,9),I=1,51)/
     *  271., 540., 825., 975.,1140.,1285.,1400.,1480.,1555.,1580.,
     * 1525.,1470.,1360.,1340.,1255.,1160.,1120.,1085.,1060.,1045.,
     * 1045.,1045.,1065.,1075.,1090.,1025.,1005., 940., 917., 880.,
     *  860., 860., 860., 860.,  17*0./
C                          ANTI K0
      DATA (CROSEC(I,1,15,10),I=1,51)/
     *20.,40.,50.,60.,70.,80.,90.,100.,120.,140.,160.,180.,
     *200.,220.,250.,300.,350.,400.,450.,500.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  16*0.,34./
      DATA (CROSEC(I,2,15,10),I=1,51)/
     * 2430.,2610.,2710.,2790.,2880.,2940.,2965.,2970.,2970.,2920.,
     * 2840.,2720.,2570.,2500.,2365.,2200.,2050.,1926.,1825.,1768.,
     * 1749.,1750.,1778.,1789.,1808.,1690.,1645.,1530.,1492.,1450.,
     * 1425.,1425.,1425.,1425.,  17*0./
      DATA (CROSEC(I,3,15,10),I=1,51)/
     *  925.,1125.,1250.,1375.,1500.,1600.,1680.,1750.,1770.,1730.,
     * 1660.,1580.,1500.,1450.,1330.,1250.,1190.,1140.,1100.,1075.,
     * 1075.,1070.,1088.,1095.,1110.,1035.,1005., 940., 917., 880.,
     *  860., 860., 860., 860.,  17*0./
C
C ********************** DATA FOR CD(48,112), INUCL=16 *******************
C                                    NEUTRON
      DATA (CROSEC(I,1,16,1),I=1,51)/
     *14.,15.,17.,18.,20.,22.,25.,27.,30.,33.,35.,40.,45.,50.,55.,60.,
     *65.,70.,80.,90.,100.,120.,140.,150.,160.,180.,200.,250.,300.,
     *350.,400.,500.,600.,700.,800.,900.,1000.,1500.,2000.,3000.,
     *5000.,7000.,10000.,20000.,50000.,100000.,500000.,1000000.,
     *2*0.,48./
      DATA (CROSEC(I,2,16,1),I=1,51)/
     * 4420.,4280.,4170.,4070.,3860.,3680.,3420.,3280.,3125.,3060.,
     * 3080.,3190.,3350.,3445.,3510.,3540.,3560.,3550.,3460.,3300.,
     * 3030.,2640.,2340.,2190.,2070.,1950.,1770.,1732.,1740.,1760.,
     * 1780.,1832.,1885.,1925.,1945.,1960.,1980.,2070.,2080.,2065.,
     * 2040.,2022.,1980.,1940.,1870.,1870.,1870.,1870.,  3*0./
      DATA (CROSEC(I,3,16,1),I=1,51)/
     * 1920.,1910.,1880.,1860.,1840.,1800.,1760.,1720.,1675.,1630.,
     * 1600.,1520.,1465.,1420.,1390.,1340.,1310.,1280.,1275.,1235.,
     * 1225.,1200.,1170.,1170.,1170.,1165.,1145.,1140.,1140.,1135.,
     * 1160.,1180.,1220.,1240.,1250.,1260.,1265.,1270.,1275.,1250.,
     * 1222.,1222.,1220.,1215.,1190.,1190.,1190.,1190.,  3*0./
C                                    PROTON
      DATA (CROSEC(I,1,16,2),I=1,51)/
     *14.,15.,17.,18.,20.,22.,25.,27.,30.,33.,35.,40.,45.,50.,55.,60.,
     *65.,70.,80.,90.,100.,120.,140.,150.,160.,180.,200.,250.,300.,
     *350.,400.,500.,600.,700.,800.,900.,1000.,1500.,2000.,3000.,
     *5000.,7000.,10000.,20000.,50000.,100000.,500000.,1000000.,
     *2*0.,48./
      DATA (CROSEC(I,2,16,2),I=1,51)/
     * 4420.,4280.,4170.,4070.,3860.,3680.,3420.,3280.,3125.,3060.,
     * 3080.,3190.,3350.,3445.,3510.,3540.,3560.,3550.,3460.,3300.,
     * 3030.,2640.,2340.,2190.,2070.,1950.,1770.,1732.,1740.,1760.,
     * 1780.,1832.,1885.,1925.,1945.,1960.,1980.,2070.,2080.,2065.,
     * 2040.,2022.,1980.,1940.,1870.,1870.,1870.,1870.,  3*0./
      DATA (CROSEC(I,3,16,2),I=1,51)/
     * 1020.,1100.,1225.,1290.,1440.,1520.,1575.,1560.,1518.,1460.,
     * 1420.,1400.,1365.,1340.,1300.,1280.,1260.,1200.,1190.,1160.,
     * 1125.,1125.,1125.,1125.,1125.,1125.,1120.,1120.,1120.,1118.,
     * 1146.,1180.,1220.,1240.,1250.,1260.,1265.,1270.,1275.,1250.,
     * 1222.,1222.,1220.,1215.,1190.,1190.,1190.,1190.,  3*0./
C                                    PI-
      DATA (CROSEC(I,1,16,3),I=1,51)/
     *20.,40.,50.,60.,70.,80.,90.,100.,120.,140.,160.,180.,
     *200.,220.,250.,300.,350.,400.,450.,500.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  16*0.,34./
      DATA (CROSEC(I,2,16,3),I=1,51)/
     * 3060.,3125.,3170.,3220.,3255.,3280.,3290.,3260.,3270.,3200.,
     * 3120.,3080.,3090.,2920.,2810.,2640.,2362.,2230.,2115.,2050.,
     * 2020.,2025.,2040.,2070.,2100.,1900.,1795.,1740.,1675.,1645.,
     * 1625.,1620.,1620.,1620.,  17*0./
      DATA (CROSEC(I,3,16,3),I=1,51)/
     * 1025.,1275.,1440.,1625.,1740.,1800.,1880.,1920.,1980.,1920.,
     * 1850.,1810.,1720.,1650.,1560.,1450.,1330.,1290.,1245.,1210.,
     * 1200.,1200.,1205.,1205.,1230.,1130.,1085.,1060.,1000., 985.,
     *  975., 970., 970., 970.,  17*0./
C                                    PI+
      DATA (CROSEC(I,1,16,4),I=1,51)/
     *20.,40.,50.,60.,70.,80.,90.,100.,120.,140.,160.,180.,
     *200.,220.,250.,300.,350.,400.,450.,500.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  16*0.,34./
      DATA (CROSEC(I,2,16,4),I=1,51)/
     *  455., 780.,1170.,1700.,2120.,2400.,2600.,2720.,2820.,2840.,
     * 2800.,2760.,2720.,2640.,2560.,2450.,2252.,2130.,2035.,1985.,
     * 1970.,1975.,2005.,2035.,2070.,1880.,1785.,1740.,1675.,1645.,
     * 1625.,1620.,1620.,1620.,  17*0./
      DATA (CROSEC(I,3,16,4),I=1,51)/
     *  310., 580., 880.,1060.,1270.,1400.,1530.,1610.,1660.,1680.,
     * 1640.,1600.,1560.,1500.,1430.,1330.,1280.,1230.,1200.,1180.,
     * 1170.,1175.,1180.,1180.,1210.,1120.,1080.,1060.,1000., 985.,
     *  975., 970., 970., 970.,  17*0./
C                                    ANTINEUTRON
      DATA (CROSEC(I,1,16,5),I=1,51)/
     *14.,15.,17.,18.,20.,22.,25.,27.,30.,33.,35.,40.,45.,50.,55.,60.,
     *65.,70.,80.,90.,100.,120.,140.,150.,160.,180.,200.,250.,300.,
     *350.,400.,500.,600.,700.,800.,900.,1000.,1500.,2000.,3000.,
     *5000.,7000.,10000.,20000.,50000.,100000.,500000.,1000000.,
     *2*0.,48./
      DATA (CROSEC(I,2,16,5),I=1,51)/
     * 4420.,4280.,4170.,4070.,3860.,3680.,3420.,3280.,3125.,3060.,
     * 3080.,3190.,3350.,3445.,3510.,3540.,3560.,3550.,3460.,3300.,
     * 3030.,2640.,2340.,2190.,2070.,1950.,1770.,1732.,1740.,1760.,
     * 1780.,1832.,1885.,1925.,1945.,1960.,1980.,2070.,2080.,2065.,
     * 2040.,2022.,1980.,1940.,1870.,1870.,1870.,1870.,  3*0./
      DATA (CROSEC(I,3,16,5),I=1,51)/
     * 1920.,1910.,1880.,1860.,1840.,1800.,1760.,1720.,1675.,1630.,
     * 1600.,1520.,1465.,1420.,1390.,1340.,1310.,1280.,1275.,1235.,
     * 1225.,1200.,1170.,1170.,1170.,1165.,1145.,1140.,1140.,1135.,
     * 1160.,1180.,1220.,1240.,1250.,1260.,1265.,1270.,1275.,1250.,
     * 1222.,1222.,1220.,1215.,1190.,1190.,1190.,1190.,  3*0./
C                                    ANTIPROTON
      DATA (CROSEC(I,1,16,6),I=1,51)/
     *14.,15.,17.,18.,20.,22.,25.,27.,30.,33.,35.,40.,45.,50.,55.,60.,
     *65.,70.,80.,90.,100.,120.,140.,150.,160.,180.,200.,250.,300.,
     *350.,400.,500.,600.,700.,800.,900.,1000.,1500.,2000.,3000.,
     *5000.,7000.,10000.,20000.,50000.,100000.,500000.,1000000.,
     *2*0.,48./
      DATA (CROSEC(I,2,16,6),I=1,51)/
     * 4420.,4280.,4170.,4070.,3860.,3680.,3420.,3280.,3125.,3060.,
     * 3080.,3190.,3350.,3445.,3510.,3540.,3560.,3550.,3460.,3300.,
     * 3030.,2640.,2340.,2190.,2070.,1950.,1770.,1732.,1740.,1760.,
     * 1780.,1832.,1885.,1925.,1945.,1960.,1980.,2070.,2080.,2065.,
     * 2040.,2022.,1980.,1940.,1870.,1870.,1870.,1870.,  3*0./
      DATA (CROSEC(I,3,16,6),I=1,51)/
     * 1020.,1100.,1225.,1290.,1440.,1520.,1575.,1560.,1518.,1460.,
     * 1420.,1400.,1365.,1340.,1300.,1280.,1260.,1200.,1190.,1160.,
     * 1125.,1125.,1125.,1125.,1125.,1125.,1120.,1120.,1120.,1118.,
     * 1146.,1180.,1220.,1240.,1250.,1260.,1265.,1270.,1275.,1250.,
     * 1222.,1222.,1220.,1215.,1190.,1190.,1190.,1190.,  3*0./
C                                    K-
      DATA (CROSEC(I,1,16,7),I=1,51)/
     *20.,40.,50.,60.,70.,80.,90.,100.,120.,140.,160.,180.,
     *200.,220.,250.,300.,350.,400.,450.,500.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  16*0.,34./
      DATA (CROSEC(I,2,16,7),I=1,51)/
     * 3060.,3125.,3170.,3220.,3255.,3280.,3290.,3260.,3270.,3200.,
     * 3120.,3080.,3090.,2920.,2810.,2640.,2362.,2230.,2115.,2050.,
     * 2020.,2025.,2040.,2070.,2100.,1900.,1795.,1740.,1675.,1645.,
     * 1625.,1620.,1620.,1620.,  17*0./
      DATA (CROSEC(I,3,16,7),I=1,51)/
     * 1025.,1275.,1440.,1625.,1740.,1800.,1880.,1920.,1980.,1920.,
     * 1850.,1810.,1720.,1650.,1560.,1450.,1330.,1290.,1245.,1210.,
     * 1200.,1200.,1205.,1205.,1230.,1130.,1085.,1060.,1000., 985.,
     *  975., 970., 970., 970.,  17*0./
C                                    K+
      DATA (CROSEC(I,1,16,8),I=1,51)/
     *20.,40.,50.,60.,70.,80.,90.,100.,120.,140.,160.,180.,
     *200.,220.,250.,300.,350.,400.,450.,500.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  16*0.,34./
      DATA (CROSEC(I,2,16,8),I=1,51)/
     *  455., 780.,1170.,1700.,2120.,2400.,2600.,2720.,2820.,2840.,
     * 2800.,2760.,2720.,2640.,2560.,2450.,2252.,2130.,2035.,1985.,
     * 1970.,1975.,2005.,2035.,2070.,1880.,1785.,1740.,1675.,1645.,
     * 1625.,1620.,1620.,1620.,  17*0./
      DATA (CROSEC(I,3,16,8),I=1,51)/
     *  310., 580., 880.,1060.,1270.,1400.,1530.,1610.,1660.,1680.,
     * 1640.,1600.,1560.,1500.,1430.,1330.,1280.,1230.,1200.,1180.,
     * 1170.,1175.,1180.,1180.,1210.,1120.,1080.,1060.,1000., 985.,
     *  975., 970., 970., 970.,  17*0./
C                                    K0
      DATA (CROSEC(I,1,16,9),I=1,51)/
     *20.,40.,50.,60.,70.,80.,90.,100.,120.,140.,160.,180.,
     *200.,220.,250.,300.,350.,400.,450.,500.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  16*0.,34./
      DATA (CROSEC(I,2,16,9),I=1,51)/
     *  455., 780.,1170.,1700.,2120.,2400.,2600.,2720.,2820.,2840.,
     * 2800.,2760.,2720.,2640.,2560.,2450.,2252.,2130.,2035.,1985.,
     * 1970.,1975.,2005.,2035.,2070.,1880.,1785.,1740.,1675.,1645.,
     * 1625.,1620.,1620.,1620.,  17*0./
      DATA (CROSEC(I,3,16,9),I=1,51)/
     *  310., 580., 880.,1060.,1270.,1400.,1530.,1610.,1660.,1680.,
     * 1640.,1600.,1560.,1500.,1430.,1330.,1280.,1230.,1200.,1180.,
     * 1170.,1175.,1180.,1180.,1210.,1120.,1080.,1060.,1000., 985.,
     *  975., 970., 970., 970.,  17*0./
C                              ANTI K0
      DATA (CROSEC(I,1,16,10),I=1,51)/
     *20.,40.,50.,60.,70.,80.,90.,100.,120.,140.,160.,180.,
     *200.,220.,250.,300.,350.,400.,450.,500.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  16*0.,34./
      DATA (CROSEC(I,2,16,10),I=1,51)/
     * 3060.,3125.,3170.,3220.,3255.,3280.,3290.,3260.,3270.,3200.,
     * 3120.,3080.,3090.,2920.,2810.,2640.,2362.,2230.,2115.,2050.,
     * 2020.,2025.,2040.,2070.,2100.,1900.,1795.,1740.,1675.,1645.,
     * 1625.,1620.,1620.,1620.,  17*0./
      DATA (CROSEC(I,3,16,10),I=1,51)/
     * 1025.,1275.,1440.,1625.,1740.,1800.,1880.,1920.,1980.,1920.,
     * 1850.,1810.,1720.,1650.,1560.,1450.,1330.,1290.,1245.,1210.,
     * 1200.,1200.,1205.,1205.,1230.,1130.,1085.,1060.,1000., 985.,
     *  975., 970., 970., 970.,  17*0./
C
C ********************** DATA FOR SN(50,119), INUCL=17 ***************
C                                    NEUTRON
      DATA (CROSEC(I,1,17,1),I=1,51)/
     *14.,15.,17.,18.,20.,22.,25.,27.,30.,33.,35.,40.,45.,50.,55.,60.,
     *65.,70.,80.,90.,100.,120.,140.,150.,160.,180.,200.,250.,300.,
     *350.,400.,500.,600.,700.,800.,900.,1000.,1500.,2000.,3000.,
     *5000.,7000.,10000.,20000.,50000.,100000.,500000.,1000000.,
     *2*0.,48./
      DATA (CROSEC(I,2,17,1),I=1,51)/
     * 4420.,4400.,4260.,4150.,3980.,3770.,3530.,3370.,3245.,3180.,
     * 3170.,3260.,3400.,3500.,3560.,3610.,3650.,3680.,3580.,3390.,
     * 3190.,2760.,2430.,2295.,2175.,1990.,1880.,1810.,1820.,1840.,
     * 1865.,1940.,1985.,2020.,2040.,2060.,2080.,2160.,2185.,2180.,
     * 2110.,2105.,2080.,2050.,1980.,1980.,1980.,1980.,  3*0./
      DATA (CROSEC(I,3,17,1),I=1,51)/
     * 1945.,1940.,1905.,1890.,1860.,1830.,1780.,1755.,1717.,1680.,
     * 1645.,1570.,1500.,1455.,1410.,1370.,1340.,1320.,1290.,1285.,
     * 1260.,1240.,1235.,1212.,1200.,1200.,1200.,1190.,1190.,1200.,
     * 1210.,1240.,1270.,1285.,1300.,1300.,1310.,1320.,1320.,1290.,
     * 1240.,1240.,1240.,1240.,1240.,1240.,1240.,1240.,  3*0./
C                                    PROTON
      DATA (CROSEC(I,1,17,2),I=1,51)/
     *14.,15.,17.,18.,20.,22.,25.,27.,30.,33.,35.,40.,45.,50.,55.,60.,
     *65.,70.,80.,90.,100.,120.,140.,150.,160.,180.,200.,250.,300.,
     *350.,400.,500.,600.,700.,800.,900.,1000.,1500.,2000.,3000.,
     *5000.,7000.,10000.,20000.,50000.,100000.,500000.,1000000.,
     *2*0.,48./
      DATA (CROSEC(I,2,17,2),I=1,51)/
     * 4420.,4400.,4260.,4150.,3980.,3770.,3530.,3370.,3245.,3180.,
     * 3170.,3260.,3400.,3500.,3560.,3610.,3650.,3680.,3580.,3390.,
     * 3190.,2760.,2430.,2295.,2175.,1990.,1880.,1810.,1820.,1840.,
     * 1865.,1940.,1985.,2020.,2040.,2060.,2080.,2160.,2185.,2180.,
     * 2110.,2105.,2080.,2050.,1980.,1980.,1980.,1980.,  3*0./
      DATA (CROSEC(I,3,17,2),I=1,51)/
     * 1020.,1080.,1270.,1335.,1465.,1505.,1610.,1610.,1550.,1535.,
     * 1500.,1440.,1407.,1370.,1340.,1300.,1285.,1260.,1230.,1215.,
     * 1200.,1180.,1170.,1170.,1165.,1165.,1170.,1165.,1165.,1183.,
     * 1195.,1240.,1270.,1285.,1300.,1300.,1310.,1320.,1320.,1290.,
     * 1240.,1240.,1240.,1240.,1240.,1240.,1240.,1240.,  3*0./
C                                    PI-
      DATA (CROSEC(I,1,17,3),I=1,51)/
     *20.,40.,50.,60.,70.,80.,90.,100.,120.,140.,160.,180.,
     *200.,220.,250.,300.,350.,400.,450.,500.,550.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  15*0.,35./
      DATA (CROSEC(I,2,17,3),I=1,51)/
     * 3000.,3180.,3250.,3300.,3380.,3410.,3470.,3450.,3410.,3350.,
     * 3280.,3200.,3120.,3050.,2900.,2630.,2500.,2325.,2190.,2100.,
     * 2060.,2055.,2055.,2055.,2067.,2085.,2000.,1900.,1835.,1770.,
     * 1720.,1700.,1695.,1695.,1695.,  16*0./
      DATA (CROSEC(I,3,17,3),I=1,51)/
     * 1050.,1350.,1520.,1650.,1800.,1980.,2070.,2120.,2090.,2050.,
     * 1980.,1920.,1830.,1770.,1670.,1500.,1435.,1350.,1300.,1230.,
     * 1220.,1235.,1235.,1235.,1237.,1240.,1160.,1120.,1090.,1065.,
     * 1040.,1020.,1015.,1015.,1015.,  16*0./
C                                    PI+
      DATA (CROSEC(I,1,17,4),I=1,51)/
     *20.,40.,50.,60.,70.,80.,90.,100.,120.,140.,160.,180.,
     *200.,220.,250.,300.,350.,400.,450.,500.,550.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  15*0.,35./
      DATA (CROSEC(I,2,17,4),I=1,51)/
     *  465., 800.,1200.,1760.,2170.,2480.,2730.,2885.,2970.,2980.,
     * 2970.,2890.,2840.,2790.,2620.,2450.,2335.,2205.,2080.,2020.,
     * 2010.,1990.,1990.,2015.,2030.,2045.,1980.,1890.,1835.,1770.,
     * 1720.,1700.,1695.,1695.,1695.,  16*0./
      DATA (CROSEC(I,3,17,4),I=1,51)/
     *  315., 590., 880.,1220.,1460.,1580.,1700.,1770.,1810.,1810.,
     * 1800.,1730.,1680.,1630.,1530.,1400.,1335.,1270.,1210.,1180.,
     * 1190.,1190.,1190.,1205.,1210.,1210.,1150.,1115.,1090.,1065.,
     * 1040.,1020.,1015.,1015.,1015.,  16*0./
C                                    ANTINEUTRON
      DATA (CROSEC(I,1,17,5),I=1,51)/
     *14.,15.,17.,18.,20.,22.,25.,27.,30.,33.,35.,40.,45.,50.,55.,60.,
     *65.,70.,80.,90.,100.,120.,140.,150.,160.,180.,200.,250.,300.,
     *350.,400.,500.,600.,700.,800.,900.,1000.,1500.,2000.,3000.,
     *5000.,7000.,10000.,20000.,50000.,100000.,500000.,1000000.,
     *2*0.,48./
      DATA (CROSEC(I,2,17,5),I=1,51)/
     * 4420.,4400.,4260.,4150.,3980.,3770.,3530.,3370.,3245.,3180.,
     * 3170.,3260.,3400.,3500.,3560.,3610.,3650.,3680.,3580.,3390.,
     * 3190.,2760.,2430.,2295.,2175.,1990.,1880.,1810.,1820.,1840.,
     * 1865.,1940.,1985.,2020.,2040.,2060.,2080.,2160.,2185.,2180.,
     * 2110.,2105.,2080.,2050.,1980.,1980.,1980.,1980.,  3*0./
      DATA (CROSEC(I,3,17,5),I=1,51)/
     * 1945.,1940.,1905.,1890.,1860.,1830.,1780.,1755.,1717.,1680.,
     * 1645.,1570.,1500.,1455.,1410.,1370.,1340.,1320.,1290.,1285.,
     * 1260.,1240.,1235.,1212.,1200.,1200.,1200.,1190.,1190.,1200.,
     * 1210.,1240.,1270.,1285.,1300.,1300.,1310.,1320.,1320.,1290.,
     * 1240.,1240.,1240.,1240.,1240.,1240.,1240.,1240.,  3*0./
C                                    ANTIPROTON
      DATA (CROSEC(I,1,17,6),I=1,51)/
     *14.,15.,17.,18.,20.,22.,25.,27.,30.,33.,35.,40.,45.,50.,55.,60.,
     *65.,70.,80.,90.,100.,120.,140.,150.,160.,180.,200.,250.,300.,
     *350.,400.,500.,600.,700.,800.,900.,1000.,1500.,2000.,3000.,
     *5000.,7000.,10000.,20000.,50000.,100000.,500000.,1000000.,
     *2*0.,48./
      DATA (CROSEC(I,2,17,6),I=1,51)/
     * 4420.,4400.,4260.,4150.,3980.,3770.,3530.,3370.,3245.,3180.,
     * 3170.,3260.,3400.,3500.,3560.,3610.,3650.,3680.,3580.,3390.,
     * 3190.,2760.,2430.,2295.,2175.,1990.,1880.,1810.,1820.,1840.,
     * 1865.,1940.,1985.,2020.,2040.,2060.,2080.,2160.,2185.,2180.,
     * 2110.,2105.,2080.,2050.,1980.,1980.,1980.,1980.,  3*0./
      DATA (CROSEC(I,3,17,6),I=1,51)/
     * 1020.,1080.,1270.,1335.,1465.,1505.,1610.,1610.,1550.,1535.,
     * 1500.,1440.,1407.,1370.,1340.,1300.,1285.,1260.,1230.,1215.,
     * 1200.,1180.,1170.,1170.,1165.,1165.,1170.,1165.,1165.,1183.,
     * 1195.,1240.,1270.,1285.,1300.,1300.,1310.,1320.,1320.,1290.,
     * 1240.,1240.,1240.,1240.,1240.,1240.,1240.,1240.,  3*0./
C                                    K-
      DATA (CROSEC(I,1,17,7),I=1,51)/
     *20.,40.,50.,60.,70.,80.,90.,100.,120.,140.,160.,180.,
     *200.,220.,250.,300.,350.,400.,450.,500.,550.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  15*0.,35./
      DATA (CROSEC(I,2,17,7),I=1,51)/
     * 3000.,3180.,3250.,3300.,3380.,3410.,3470.,3450.,3410.,3350.,
     * 3280.,3200.,3120.,3050.,2900.,2630.,2500.,2325.,2190.,2100.,
     * 2060.,2055.,2055.,2055.,2067.,2085.,2000.,1900.,1835.,1770.,
     * 1720.,1700.,1695.,1695.,1695.,  16*0./
      DATA (CROSEC(I,3,17,7),I=1,51)/
     * 1050.,1350.,1520.,1650.,1800.,1980.,2070.,2120.,2090.,2050.,
     * 1980.,1920.,1830.,1770.,1670.,1500.,1435.,1350.,1300.,1230.,
     * 1220.,1235.,1235.,1235.,1237.,1240.,1160.,1120.,1090.,1065.,
     * 1040.,1020.,1015.,1015.,1015.,  16*0./
C                                    K+
      DATA (CROSEC(I,1,17,8),I=1,51)/
     *20.,40.,50.,60.,70.,80.,90.,100.,120.,140.,160.,180.,
     *200.,220.,250.,300.,350.,400.,450.,500.,550.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  15*0.,35./
      DATA (CROSEC(I,2,17,8),I=1,51)/
     *  465., 800.,1200.,1760.,2170.,2480.,2730.,2885.,2970.,2980.,
     * 2970.,2890.,2840.,2790.,2620.,2450.,2335.,2205.,2080.,2020.,
     * 2010.,1990.,1990.,2015.,2030.,2045.,1980.,1890.,1835.,1770.,
     * 1720.,1700.,1695.,1695.,1695.,  16*0./
      DATA (CROSEC(I,3,17,8),I=1,51)/
     *  315., 590., 880.,1220.,1460.,1580.,1700.,1770.,1810.,1810.,
     * 1800.,1730.,1680.,1630.,1530.,1400.,1335.,1270.,1210.,1180.,
     * 1190.,1190.,1190.,1205.,1210.,1210.,1150.,1115.,1090.,1065.,
     * 1040.,1020.,1015.,1015.,1015.,  16*0./
C                                    K0
      DATA (CROSEC(I,1,17,9),I=1,51)/
     *20.,40.,50.,60.,70.,80.,90.,100.,120.,140.,160.,180.,
     *200.,220.,250.,300.,350.,400.,450.,500.,550.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  15*0.,35./
      DATA (CROSEC(I,2,17,9),I=1,51)/
     *  465., 800.,1200.,1760.,2170.,2480.,2730.,2885.,2970.,2980.,
     * 2970.,2890.,2840.,2790.,2620.,2450.,2335.,2205.,2080.,2020.,
     * 2010.,1990.,1990.,2015.,2030.,2045.,1980.,1890.,1835.,1770.,
     * 1720.,1700.,1695.,1695.,1695.,  16*0./
      DATA (CROSEC(I,3,17,9),I=1,51)/
     *  315., 590., 880.,1220.,1460.,1580.,1700.,1770.,1810.,1810.,
     * 1800.,1730.,1680.,1630.,1530.,1400.,1335.,1270.,1210.,1180.,
     * 1190.,1190.,1190.,1205.,1210.,1210.,1150.,1115.,1090.,1065.,
     * 1040.,1020.,1015.,1015.,1015.,  16*0./
C                               ANTI K0
      DATA (CROSEC(I,1,17,10),I=1,51)/
     *20.,40.,50.,60.,70.,80.,90.,100.,120.,140.,160.,180.,
     *200.,220.,250.,300.,350.,400.,450.,500.,550.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  15*0.,35./
      DATA (CROSEC(I,2,17,10),I=1,51)/
     * 3000.,3180.,3250.,3300.,3380.,3410.,3470.,3450.,3410.,3350.,
     * 3280.,3200.,3120.,3050.,2900.,2630.,2500.,2325.,2190.,2100.,
     * 2060.,2055.,2055.,2055.,2067.,2085.,2000.,1900.,1835.,1770.,
     * 1720.,1700.,1695.,1695.,1695.,  16*0./
      DATA (CROSEC(I,3,17,10),I=1,51)/
     * 1050.,1350.,1520.,1650.,1800.,1980.,2070.,2120.,2090.,2050.,
     * 1980.,1920.,1830.,1770.,1670.,1500.,1435.,1350.,1300.,1230.,
     * 1220.,1235.,1235.,1235.,1237.,1240.,1160.,1120.,1090.,1065.,
     * 1040.,1020.,1015.,1015.,1015.,  16*0./
C
C *********************DATA FOR W(74,184), INUCL=18 ********************
C                                    NEUTRON
      DATA (CROSEC(I,1,18,1),I=1,51)/
     *14.,15.,17.,18.,20.,22.,25.,27.,30.,33.,35.,40.,45.,50.,55.,60.,
     *65.,70.,80.,90.,100.,120.,140.,150.,160.,180.,200.,250.,300.,
     *350.,400.,500.,600.,700.,800.,900.,1000.,1500.,2000.,3000.,
     *5000.,7000.,10000.,20000.,50000.,100000.,500000.,1000000.,
     *2*0.,48./
      DATA (CROSEC(I,2,18,1),I=1,51)/
     * 5320.,5430.,5480.,5450.,5330.,5190.,4960.,4790.,4550.,4340.,
     * 4200.,4070.,4000.,4030.,4125.,4220.,4270.,4390.,4440.,4360.,
     * 4200.,3800.,3380.,3200.,3040.,2790.,2660.,2575.,2575.,2600.,
     * 2640.,2690.,2755.,2790.,2812.,2837.,2850.,2950.,3000.,2970.,
     * 2940.,2910.,2880.,2820.,2730.,2730.,2730.,2730.,  3*0./
      DATA (CROSEC(I,3,18,1),I=1,51)/
     * 2440.,2400.,2370.,2350.,2310.,2270.,2220.,2195.,2150.,2100.,
     * 2070.,2010.,1945.,1900.,1850.,1820.,1780.,1760.,1730.,1720.,
     * 1680.,1680.,1660.,1660.,1650.,1650.,1640.,1640.,1612.,1615.,
     * 1625.,1640.,1700.,1720.,1730.,1740.,1750.,1780.,1780.,1750.,
     * 1740.,1735.,1710.,1695.,1680.,1680.,1680.,1680.,  3*0./
C                                    PROTON
      DATA (CROSEC(I,1,18,2),I=1,51)/
     *14.,15.,17.,18.,20.,22.,25.,27.,30.,33.,35.,40.,45.,50.,55.,60.,
     *65.,70.,80.,90.,100.,120.,140.,150.,160.,180.,200.,250.,300.,
     *350.,400.,500.,600.,700.,800.,900.,1000.,1500.,2000.,3000.,
     *5000.,7000.,10000.,20000.,50000.,100000.,500000.,1000000.,
     *2*0.,48./
      DATA (CROSEC(I,2,18,2),I=1,51)/
     * 5320.,5430.,5480.,5450.,5330.,5190.,4960.,4790.,4550.,4340.,
     * 4200.,4070.,4000.,4030.,4125.,4220.,4270.,4390.,4440.,4360.,
     * 4200.,3800.,3380.,3200.,3040.,2790.,2660.,2575.,2575.,2600.,
     * 2640.,2690.,2755.,2790.,2812.,2837.,2850.,2950.,3000.,2970.,
     * 2940.,2910.,2880.,2820.,2730.,2730.,2730.,2730.,  3*0./
      DATA (CROSEC(I,3,18,2),I=1,51)/
     *  950.,1020.,1240.,1400.,1560.,1670.,1760.,1830.,1850.,1855.,
     * 1870.,1840.,1800.,1770.,1740.,1715.,1680.,1670.,1650.,1620.,
     * 1610.,1600.,1600.,1600.,1600.,1600.,1600.,1595.,1585.,1595.,
     * 1615.,1640.,1700.,1720.,1730.,1740.,1750.,1780.,1780.,1750.,
     * 1740.,1735.,1710.,1695.,1680.,1680.,1680.,1680.,  3*0./
C                                    PI-
      DATA (CROSEC(I,1,18,3),I=1,51)/
     *20.,40.,50.,60.,70.,80.,90.,100.,120.,140.,160.,180.,
     *200.,220.,250.,300.,350.,400.,450.,500.,550.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  15*0.,35./
      DATA (CROSEC(I,2,18,3),I=1,51)/
     * 5200.,5115.,5025.,4975.,4900.,4850.,4780.,4725.,4600.,4490.,
     * 4355.,4255.,4125.,4040.,3830.,3580.,3330.,3110.,2955.,2860.,
     * 2852.,2845.,2885.,2900.,2915.,2940.,2800.,2660.,2580.,2490.,
     * 2460.,2425.,2420.,2420.,2420.,  16*0./
      DATA (CROSEC(I,3,18,3),I=1,51)/
     * 1450.,1850.,2100.,2350.,2550.,2700.,2825.,2900.,2850.,2750.,
     * 2630.,2525.,2400.,2300.,2200.,2070.,1880.,1770.,1715.,1680.,
     * 1680.,1680.,1685.,1690.,1700.,1720.,1635.,1560.,1530.,1460.,
     * 1440.,1410.,1410.,1410.,1410.,  16*0./
C                                    PI+
      DATA (CROSEC(I,1,18,4),I=1,51)/
     *20.,40.,50.,60.,70.,80.,90.,100.,120.,140.,160.,180.,
     *200.,220.,250.,300.,350.,400.,450.,500.,550.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  15*0.,35./
      DATA (CROSEC(I,2,18,4),I=1,51)/
     *  480., 900.,1500.,2350.,3020.,3420.,3650.,3775.,3875.,3830.,
     * 3750.,3700.,3630.,3550.,3550.,3290.,3070.,2890.,2840.,2730.,
     * 2725.,2720.,2770.,2805.,2828.,2865.,2770.,2640.,2570.,2490.,
     * 2460.,2425.,2420.,2420.,2420.,  16*0./
      DATA (CROSEC(I,3,18,4),I=1,51)/
     *  341., 680., 990.,1500.,1850.,2150.,2250.,2300.,2350.,2330.,
     * 2280.,2230.,2200.,2120.,2130.,1900.,1780.,1670.,1635.,1600.,
     * 1602.,1605.,1610.,1615.,1630.,1660.,1620.,1550.,1525.,1460.,
     * 1440.,1410.,1410.,1410.,1410.,  16*0./
C                                    ANTINEUTRON
      DATA (CROSEC(I,1,18,5),I=1,51)/
     *14.,15.,17.,18.,20.,22.,25.,27.,30.,33.,35.,40.,45.,50.,55.,60.,
     *65.,70.,80.,90.,100.,120.,140.,150.,160.,180.,200.,250.,300.,
     *350.,400.,500.,600.,700.,800.,900.,1000.,1500.,2000.,3000.,
     *5000.,7000.,10000.,20000.,50000.,100000.,500000.,1000000.,
     *2*0.,48./
      DATA (CROSEC(I,2,18,5),I=1,51)/
     * 5320.,5430.,5480.,5450.,5330.,5190.,4960.,4790.,4550.,4340.,
     * 4200.,4070.,4000.,4030.,4125.,4220.,4270.,4390.,4440.,4360.,
     * 4200.,3800.,3380.,3200.,3040.,2790.,2660.,2575.,2575.,2600.,
     * 2640.,2690.,2755.,2790.,2812.,2837.,2850.,2950.,3000.,2970.,
     * 2940.,2910.,2880.,2820.,2730.,2730.,2730.,2730.,  3*0./
      DATA (CROSEC(I,3,18,5),I=1,51)/
     * 2440.,2400.,2370.,2350.,2310.,2270.,2220.,2195.,2150.,2100.,
     * 2070.,2010.,1945.,1900.,1850.,1820.,1780.,1760.,1730.,1720.,
     * 1680.,1680.,1660.,1660.,1650.,1650.,1640.,1640.,1612.,1615.,
     * 1625.,1640.,1700.,1720.,1730.,1740.,1750.,1780.,1780.,1750.,
     * 1740.,1735.,1710.,1695.,1680.,1680.,1680.,1680.,  3*0./
C                                    ANTIPROTON
      DATA (CROSEC(I,1,18,6),I=1,51)/
     *14.,15.,17.,18.,20.,22.,25.,27.,30.,33.,35.,40.,45.,50.,55.,60.,
     *65.,70.,80.,90.,100.,120.,140.,150.,160.,180.,200.,250.,300.,
     *350.,400.,500.,600.,700.,800.,900.,1000.,1500.,2000.,3000.,
     *5000.,7000.,10000.,20000.,50000.,100000.,500000.,1000000.,
     *2*0.,48./
      DATA (CROSEC(I,2,18,6),I=1,51)/
     * 5320.,5430.,5480.,5450.,5330.,5190.,4960.,4790.,4550.,4340.,
     * 4200.,4070.,4000.,4030.,4125.,4220.,4270.,4390.,4440.,4360.,
     * 4200.,3800.,3380.,3200.,3040.,2790.,2660.,2575.,2575.,2600.,
     * 2640.,2690.,2755.,2790.,2812.,2837.,2850.,2950.,3000.,2970.,
     * 2940.,2910.,2880.,2820.,2730.,2730.,2730.,2730.,  3*0./
      DATA (CROSEC(I,3,18,6),I=1,51)/
     *  950.,1020.,1240.,1400.,1560.,1670.,1760.,1830.,1850.,1855.,
     * 1870.,1840.,1800.,1770.,1740.,1715.,1680.,1670.,1650.,1620.,
     * 1610.,1600.,1600.,1600.,1600.,1600.,1600.,1595.,1585.,1595.,
     * 1615.,1640.,1700.,1720.,1730.,1740.,1750.,1780.,1780.,1750.,
     * 1740.,1735.,1710.,1695.,1680.,1680.,1680.,1680.,  3*0./
C                                    K-
      DATA (CROSEC(I,1,18,7),I=1,51)/
     *20.,40.,50.,60.,70.,80.,90.,100.,120.,140.,160.,180.,
     *200.,220.,250.,300.,350.,400.,450.,500.,550.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  15*0.,35./
      DATA (CROSEC(I,2,18,7),I=1,51)/
     * 5200.,5115.,5025.,4975.,4900.,4850.,4780.,4725.,4600.,4490.,
     * 4355.,4255.,4125.,4040.,3830.,3580.,3330.,3110.,2955.,2860.,
     * 2852.,2845.,2885.,2900.,2915.,2940.,2800.,2660.,2580.,2490.,
     * 2460.,2425.,2420.,2420.,2420.,  16*0./
      DATA (CROSEC(I,3,18,7),I=1,51)/
     * 1450.,1850.,2100.,2350.,2550.,2700.,2825.,2900.,2850.,2750.,
     * 2630.,2525.,2400.,2300.,2200.,2070.,1880.,1770.,1715.,1680.,
     * 1680.,1680.,1685.,1690.,1700.,1720.,1635.,1560.,1530.,1460.,
     * 1440.,1410.,1410.,1410.,1410.,  16*0./
C                                    K+
      DATA (CROSEC(I,1,18,8),I=1,51)/
     *20.,40.,50.,60.,70.,80.,90.,100.,120.,140.,160.,180.,
     *200.,220.,250.,300.,350.,400.,450.,500.,550.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  15*0.,35./
      DATA (CROSEC(I,2,18,8),I=1,51)/
     *  480., 900.,1500.,2350.,3020.,3420.,3650.,3775.,3875.,3830.,
     * 3750.,3700.,3630.,3550.,3550.,3290.,3070.,2890.,2840.,2730.,
     * 2725.,2720.,2770.,2805.,2828.,2865.,2770.,2640.,2570.,2490.,
     * 2460.,2425.,2420.,2420.,2420.,  16*0./
      DATA (CROSEC(I,3,18,8),I=1,51)/
     *  341., 680., 990.,1500.,1850.,2150.,2250.,2300.,2350.,2330.,
     * 2280.,2230.,2200.,2120.,2130.,1900.,1780.,1670.,1635.,1600.,
     * 1602.,1605.,1610.,1615.,1630.,1660.,1620.,1550.,1525.,1460.,
     * 1440.,1410.,1410.,1410.,1410.,  16*0./
C                                    K0
      DATA (CROSEC(I,1,18,9),I=1,51)/
     *20.,40.,50.,60.,70.,80.,90.,100.,120.,140.,160.,180.,
     *200.,220.,250.,300.,350.,400.,450.,500.,550.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  15*0.,35./
      DATA (CROSEC(I,2,18,9),I=1,51)/
     *  480., 900.,1500.,2350.,3020.,3420.,3650.,3775.,3875.,3830.,
     * 3750.,3700.,3630.,3550.,3550.,3290.,3070.,2890.,2840.,2730.,
     * 2725.,2720.,2770.,2805.,2828.,2865.,2770.,2640.,2570.,2490.,
     * 2460.,2425.,2420.,2420.,2420.,  16*0./
      DATA (CROSEC(I,3,18,9),I=1,51)/
     *  341., 680., 990.,1500.,1850.,2150.,2250.,2300.,2350.,2330.,
     * 2280.,2230.,2200.,2120.,2130.,1900.,1780.,1670.,1635.,1600.,
     * 1602.,1605.,1610.,1615.,1630.,1660.,1620.,1550.,1525.,1460.,
     * 1440.,1410.,1410.,1410.,1410.,  16*0./
C                                ANTI K0
      DATA (CROSEC(I,1,18,10),I=1,51)/
     *20.,40.,50.,60.,70.,80.,90.,100.,120.,140.,160.,180.,
     *200.,220.,250.,300.,350.,400.,450.,500.,550.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  15*0.,35./
      DATA (CROSEC(I,2,18,10),I=1,51)/
     * 5200.,5115.,5025.,4975.,4900.,4850.,4780.,4725.,4600.,4490.,
     * 4355.,4255.,4125.,4040.,3830.,3580.,3330.,3110.,2955.,2860.,
     * 2852.,2845.,2885.,2900.,2915.,2940.,2800.,2660.,2580.,2490.,
     * 2460.,2425.,2420.,2420.,2420.,  16*0./
      DATA (CROSEC(I,3,18,10),I=1,51)/
     * 1450.,1850.,2100.,2350.,2550.,2700.,2825.,2900.,2850.,2750.,
     * 2630.,2525.,2400.,2300.,2200.,2070.,1880.,1770.,1715.,1680.,
     * 1680.,1680.,1685.,1690.,1700.,1720.,1635.,1560.,1530.,1460.,
     * 1440.,1410.,1410.,1410.,1410.,  16*0./
C
C ****************** DATA FOR PB(82,207), INUCL=19 *******************
C                                    NEUTRON
      DATA (CROSEC(I,1,19,1),I=1,51)/
     *14.,15.,17.,19.,20.,22.,25.,27.,30.,35.,40.,45.,50.,55.,60.,70.,
     *80.,90.,100.,120.,140.,150.,160.,180.,200.,250.,300.,350.,
     *400.,500.,600.,700.,800.,900.,1000.,1500.,2000.,3000.,
     *5000.,7000.,10000.,20000.,50000.,100000.,500000.,1000000.,
     *4*0.,46./
      DATA (CROSEC(I,2,19,1),I=1,51)/
     * 5300.,5440.,5720.,5880.,5765.,5745.,5480.,5280.,4970.,4550.,
     * 4390.,4300.,4265.,4325.,4450.,4540.,4740.,4710.,4600.,4100.,
     * 3660.,3480.,3300.,3000.,2890.,2865.,2855.,2850.,2865.,
     * 2920.,2955.,3000.,3030.,3060.,3105.,3240.,3290.,3270.,
     * 3240.,3180.,3090.,3060.,2970.,2970.,2970.,2970.,  5*0./
      DATA (CROSEC(I,3,19,1),I=1,51)/
     * 2580.,2550.,2505.,2462.,2460.,2435.,2380.,2355.,2280.,2180.,
     * 2170.,2130.,2080.,2035.,1980.,1940.,1900.,1870.,1840.,1800.,
     * 1800.,1800.,1780.,1760.,1760.,1740.,1730.,1725.,1740.,
     * 1785.,1815.,1835.,1860.,1890.,1895.,1920.,1920.,1890.,
     * 1850.,1835.,1830.,1830.,1830.,1830.,1830.,1830.,  5*0./
C                                    PROTON
      DATA (CROSEC(I,1,19,2),I=1,51)/
     *14.,15.,17.,19.,20.,22.,25.,27.,30.,35.,40.,45.,50.,55.,60.,70.,
     *80.,90.,100.,120.,140.,150.,160.,180.,200.,250.,300.,350.,
     *400.,500.,600.,700.,800.,900.,1000.,1500.,2000.,3000.,
     *5000.,7000.,10000.,20000.,50000.,100000.,500000.,1000000.,
     *4*0.,46./
      DATA (CROSEC(I,2,19,2),I=1,51)/
     * 5300.,5440.,5720.,5880.,5765.,5745.,5480.,5280.,4970.,4550.,
     * 4390.,4300.,4265.,4325.,4450.,4540.,4740.,4710.,4600.,4100.,
     * 3660.,3480.,3300.,3000.,2890.,2865.,2855.,2850.,2865.,
     * 2920.,2955.,3000.,3030.,3060.,3105.,3240.,3290.,3270.,
     * 3240.,3180.,3090.,3060.,2970.,2970.,2970.,2970.,  5*0./
      DATA (CROSEC(I,3,19,2),I=1,51)/
c     *  900.,1060.,1200.,1420.,1515.,1620.,1750.,1800.,1915.,2030.,
     * 1000.,1060.,1200.,1420.,1515.,1620.,1750.,1800.,1915.,2030., ! 2004a
     * 1960.,1940.,1910.,1860.,1840.,1780.,1770.,1760.,1740.,1720.,
     * 1725.,1740.,1740.,1730.,1720.,1700.,1710.,1720.,1730.,
     * 1740.,1815.,1835.,1860.,1890.,1895.,1920.,1920.,1890.,
     * 1850.,1835.,1830.,1830.,1830.,1830.,1830.,1830.,  5*0./
C                                    PI-
      DATA (CROSEC(I,1,19,3),I=1,51)/
     *20.,40.,50.,60.,70.,80.,90.,100.,120.,140.,160.,180.,
     *200.,220.,250.,300.,350.,400.,450.,500.,550.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  15*0.,35./
      DATA (CROSEC(I,2,19,3),I=1,51)/
     * 5890.,5700.,5610.,5580.,5550.,5480.,5400.,5300.,5100.,4930.,
     * 4750.,4600.,4400.,4280.,4170.,3915.,3650.,3470.,3260.,3150.,
     * 3120.,3070.,3085.,3100.,3120.,3160.,3070.,2930.,2820.,2750.,
     * 2710.,2655.,2640.,2640.,2640.,  16*0./
      DATA (CROSEC(I,3,19,3),I=1,51)/
     * 1575.,2025.,2300.,2575.,2850.,3000.,3115.,3180.,3080.,2940.,
     * 2800.,2670.,2550.,2450.,2370.,2220.,2110.,2000.,1920.,1880.,
     * 1850.,1800.,1805.,1810.,1820.,1840.,1800.,1720.,1640.,1620.,
     * 1570.,1530.,1530.,1530.,1530.,  16*0./
C                                    PI+
      DATA (CROSEC(I,1,19,4),I=1,51)/
     *20.,40.,50.,60.,70.,80.,90.,100.,120.,140.,160.,180.,
     *200.,220.,250.,300.,350.,400.,450.,500.,550.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  15*0.,35./
      DATA (CROSEC(I,2,19,4),I=1,51)/
     *  515., 940.,1500.,2400.,3270.,3750.,4050.,4140.,4260.,4200.,
     * 4080.,3990.,3990.,3810.,3730.,3520.,3370.,3186.,3110.,3010.,
     * 2990.,2985.,3005.,3020.,3040.,3080.,3020.,2905.,2790.,2750.,
     * 2710.,2655.,2640.,2640.,2640.,  16*0./
      DATA (CROSEC(I,3,19,4),I=1,51)/
     *  354., 707.,1040.,1650.,2100.,2400.,2580.,2640.,2650.,2520.,
     * 2410.,2330.,2250.,2190.,2130.,2000.,1930.,1870.,1830.,1790.,
     * 1770.,1765.,1775.,1780.,1790.,1800.,1775.,1710.,1620.,1620.,
     * 1570.,1530.,1530.,1530.,1530.,  16*0./
C                                    ANTINEUTRON
      DATA (CROSEC(I,1,19,5),I=1,51)/
     *14.,15.,17.,19.,20.,22.,25.,27.,30.,35.,40.,45.,50.,55.,60.,70.,
     *80.,90.,100.,120.,140.,150.,160.,180.,200.,250.,300.,350.,
     *400.,500.,600.,700.,800.,900.,1000.,1500.,2000.,3000.,
     *5000.,7000.,10000.,20000.,50000.,100000.,500000.,1000000.,
     *4*0.,46./
      DATA (CROSEC(I,2,19,5),I=1,51)/
     * 5300.,5440.,5720.,5880.,5765.,5745.,5480.,5280.,4970.,4550.,
     * 4390.,4300.,4265.,4325.,4450.,4540.,4740.,4710.,4600.,4100.,
     * 3660.,3480.,3300.,3000.,2890.,2865.,2855.,2850.,2865.,
     * 2920.,2955.,3000.,3030.,3060.,3105.,3240.,3290.,3270.,
     * 3240.,3180.,3090.,3060.,2970.,2970.,2970.,2970.,  5*0./
      DATA (CROSEC(I,3,19,5),I=1,51)/
     * 2580.,2550.,2505.,2462.,2460.,2435.,2380.,2355.,2280.,2180.,
     * 2170.,2130.,2080.,2035.,1980.,1940.,1900.,1870.,1840.,1800.,
     * 1800.,1800.,1780.,1760.,1760.,1740.,1730.,1725.,1740.,
     * 1785.,1815.,1835.,1860.,1890.,1895.,1920.,1920.,1890.,
     * 1850.,1835.,1830.,1830.,1830.,1830.,1830.,1830.,  5*0./
C                                    ANTIPROTON
      DATA (CROSEC(I,1,19,6),I=1,51)/
     *14.,15.,17.,19.,20.,22.,25.,27.,30.,35.,40.,45.,50.,55.,60.,70.,
     *80.,90.,100.,120.,140.,150.,160.,180.,200.,250.,300.,350.,
     *400.,500.,600.,700.,800.,900.,1000.,1500.,2000.,3000.,
     *5000.,7000.,10000.,20000.,50000.,100000.,500000.,1000000.,
     *4*0.,46./
      DATA (CROSEC(I,2,19,6),I=1,51)/
     * 5300.,5440.,5720.,5880.,5765.,5745.,5480.,5280.,4970.,4550.,
     * 4390.,4300.,4265.,4325.,4450.,4540.,4740.,4710.,4600.,4100.,
     * 3660.,3480.,3300.,3000.,2890.,2865.,2855.,2850.,2865.,
     * 2920.,2955.,3000.,3030.,3060.,3105.,3240.,3290.,3270.,
     * 3240.,3180.,3090.,3060.,2970.,2970.,2970.,2970.,  5*0./
      DATA (CROSEC(I,3,19,6),I=1,51)/
     *  900.,1060.,1200.,1420.,1515.,1620.,1750.,1800.,1915.,2030.,
     * 1960.,1940.,1910.,1860.,1840.,1780.,1770.,1760.,1740.,1720.,
     * 1725.,1740.,1740.,1730.,1720.,1700.,1710.,1720.,1730.,
     * 1740.,1815.,1835.,1860.,1890.,1895.,1920.,1920.,1890.,
     * 1850.,1835.,1830.,1830.,1830.,1830.,1830.,1830.,  5*0./
C                                    K-
      DATA (CROSEC(I,1,19,7),I=1,51)/
     *20.,40.,50.,60.,70.,80.,90.,100.,120.,140.,160.,180.,
     *200.,220.,250.,300.,350.,400.,450.,500.,550.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  15*0.,35./
      DATA (CROSEC(I,2,19,7),I=1,51)/
     * 5890.,5700.,5610.,5580.,5550.,5480.,5400.,5300.,5100.,4930.,
     * 4750.,4600.,4400.,4280.,4170.,3915.,3650.,3470.,3260.,3150.,
     * 3120.,3070.,3085.,3100.,3120.,3160.,3070.,2930.,2820.,2750.,
     * 2710.,2655.,2640.,2640.,2640.,  16*0./
      DATA (CROSEC(I,3,19,7),I=1,51)/
     * 1575.,2025.,2300.,2575.,2850.,3000.,3115.,3180.,3080.,2940.,
     * 2800.,2670.,2550.,2450.,2370.,2220.,2110.,2000.,1920.,1880.,
     * 1850.,1800.,1805.,1810.,1820.,1840.,1800.,1720.,1640.,1620.,
     * 1570.,1530.,1530.,1530.,1530.,  16*0./
C                                    K+
      DATA (CROSEC(I,1,19,8),I=1,51)/
     *20.,40.,50.,60.,70.,80.,90.,100.,120.,140.,160.,180.,
     *200.,220.,250.,300.,350.,400.,450.,500.,550.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  15*0.,35./
      DATA (CROSEC(I,2,19,8),I=1,51)/
     *  515., 940.,1500.,2400.,3270.,3750.,4050.,4140.,4260.,4200.,
     * 4080.,3990.,3990.,3810.,3730.,3520.,3370.,3186.,3110.,3010.,
     * 2990.,2985.,3005.,3020.,3040.,3080.,3020.,2905.,2790.,2750.,
     * 2710.,2655.,2640.,2640.,2640.,  16*0./
      DATA (CROSEC(I,3,19,8),I=1,51)/
     *  354., 707.,1040.,1650.,2100.,2400.,2580.,2640.,2650.,2520.,
     * 2410.,2330.,2250.,2190.,2130.,2000.,1930.,1870.,1830.,1790.,
     * 1770.,1765.,1775.,1780.,1790.,1800.,1775.,1710.,1620.,1620.,
     * 1570.,1530.,1530.,1530.,1530.,  16*0./
C                                    K0
      DATA (CROSEC(I,1,19,9),I=1,51)/
     *20.,40.,50.,60.,70.,80.,90.,100.,120.,140.,160.,180.,
     *200.,220.,250.,300.,350.,400.,450.,500.,550.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  15*0.,35./
      DATA (CROSEC(I,2,19,9),I=1,51)/
     *  515., 940.,1500.,2400.,3270.,3750.,4050.,4140.,4260.,4200.,
     * 4080.,3990.,3990.,3810.,3730.,3520.,3370.,3186.,3110.,3010.,
     * 2990.,2985.,3005.,3020.,3040.,3080.,3020.,2905.,2790.,2750.,
     * 2710.,2655.,2640.,2640.,2640.,  16*0./
      DATA (CROSEC(I,3,19,9),I=1,51)/
     *  354., 707.,1040.,1650.,2100.,2400.,2580.,2640.,2650.,2520.,
     * 2410.,2330.,2250.,2190.,2130.,2000.,1930.,1870.,1830.,1790.,
     * 1770.,1765.,1775.,1780.,1790.,1800.,1775.,1710.,1620.,1620.,
     * 1570.,1530.,1530.,1530.,1530.,  16*0./
C                               ANTI K0
      DATA (CROSEC(I,1,19,10),I=1,51)/
     *20.,40.,50.,60.,70.,80.,90.,100.,120.,140.,160.,180.,
     *200.,220.,250.,300.,350.,400.,450.,500.,550.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  15*0.,35./
      DATA (CROSEC(I,2,19,10),I=1,51)/
     * 5890.,5700.,5610.,5580.,5550.,5480.,5400.,5300.,5100.,4930.,
     * 4750.,4600.,4400.,4280.,4170.,3915.,3650.,3470.,3260.,3150.,
     * 3120.,3070.,3085.,3100.,3120.,3160.,3070.,2930.,2820.,2750.,
     * 2710.,2655.,2640.,2640.,2640.,  16*0./
      DATA (CROSEC(I,3,19,10),I=1,51)/
     * 1575.,2025.,2300.,2575.,2850.,3000.,3115.,3180.,3080.,2940.,
     * 2800.,2670.,2550.,2450.,2370.,2220.,2110.,2000.,1920.,1880.,
     * 1850.,1800.,1805.,1810.,1820.,1840.,1800.,1720.,1640.,1620.,
     * 1570.,1530.,1530.,1530.,1530.,  16*0./
C
C ******************** DATA FOR U(92,238), INUCL=20 **********************
C                                    NEUTRON
      DATA (CROSEC(I,1,20,1),I=1,51)/
     *14.,15.,17.,19.,20.,22.,25.,27.,30.,35.,40.,45.,50.,55.,60.,70.,
     *80.,90.,100.,120.,140.,150.,160.,180.,200.,250.,300.,350.,
     *400.,500.,600.,700.,800.,900.,1000.,1500.,2000.,3000.,
     *5000.,7000.,10000.,20000.,50000.,100000.,500000.,1000000.,
     *4*0.,46./
      DATA (CROSEC(I,2,20,1),I=1,51)/
     * 5800.,5940.,6160.,6345.,6360.,6350.,6170.,6020.,5760.,5350.,
     * 4990.,4800.,4710.,4690.,4760.,5040.,5190.,5200.,5080.,4600.,
     * 4120.,3920.,3720.,3420.,3240.,3150.,3160.,3180.,3210.,
     * 3240.,3280.,3350.,3390.,3435.,3480.,3560.,3585.,3580.,
     * 3540.,3500.,3470.,3410.,3335.,3335.,3335.,3335.,  5*0./
      DATA (CROSEC(I,3,20,1),I=1,51)/
     * 2820.,2770.,2700.,2660.,2645.,2620.,2580.,2550.,2515.,2450.,
     * 2390.,2320.,2260.,2225.,2200.,2140.,2080.,2060.,2040.,2000., 
     * 1980.,1965.,1960.,1930.,1920.,1890.,1905.,1920.,
     * 1945.,1970.,1985.,2010.,2040.,2070.,2080.,2090.,2095.,2080.,
     * 2063.,2060.,2050.,2040.,2005.,2005.,2005.,2005.,  5*0./
C                                    PROTON
      DATA (CROSEC(I,1,20,2),I=1,51)/
     *14.,15.,17.,19.,20.,22.,25.,27.,30.,35.,40.,45.,50.,55.,60.,70.,
     *80.,90.,100.,120.,140.,150.,160.,180.,200.,250.,300.,350.,
     *400.,500.,600.,700.,800.,900.,1000.,1500.,2000.,3000.,
     *5000.,7000.,10000.,20000.,50000.,100000.,500000.,1000000.,
     *4*0.,46./
      DATA (CROSEC(I,2,20,2),I=1,51)/
     * 5800.,5940.,6160.,6345.,6360.,6350.,6170.,6020.,5760.,5350.,
     * 4990.,4800.,4710.,4690.,4760.,5040.,5190.,5200.,5080.,4600.,
     * 4120.,3920.,3720.,3420.,3240.,3150.,3160.,3180.,3210.,
     * 3240.,3280.,3350.,3390.,3435.,3480.,3560.,3585.,3580.,
     * 3540.,3500.,3470.,3410.,3335.,3335.,3335.,3335.,  5*0./
      DATA (CROSEC(I,3,20,2),I=1,51)/
     *  800., 900.,1100.,1300.,1410.,1510.,1680.,1800.,2000.,2200.,
     * 2080.,2060.,2025.,2100.,2030.,2030.,2000.,1960.,1960.,1960.,
     * 1940.,1925.,1920.,1905.,1890.,1860.,1880.,1910.,1930.,
     * 1945.,1985.,2010.,2040.,2070.,2080.,2090.,2095.,2080.,
     * 2063.,2060.,2050.,2040.,2005.,2005.,2005.,2005.,  5*0./
C                                    PI-
      DATA (CROSEC(I,1,20,3),I=1,51)/
     *20.,40.,50.,60.,70.,80.,90.,100.,120.,140.,160.,180.,
     *200.,220.,250.,300.,350.,400.,450.,500.,550.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  15*0.,35./
      DATA (CROSEC(I,2,20,3),I=1,51)/
     * 7080.,6830.,6650.,6530.,6400.,6280.,6100.,5840.,5660.,5520.,
     * 5330.,5160.,4990.,4810.,4630.,4323.,4130.,3870.,3700.,3550.,
     * 3490.,3465.,3467.,3475.,3495.,3515.,3440.,3360.,3150.,3040.,
     * 2985.,2955.,2940.,2940.,2940.,  16*0./
      DATA (CROSEC(I,3,20,3),I=1,51)/
     * 1740.,2220.,2500.,2820.,3080.,3300.,3420.,3500.,3420.,3330.,
     * 3200.,3060.,2940.,2850.,2710.,2470.,2380.,2250.,2160.,2080.,
     * 2040.,2045.,2047.,2050.,2055.,2060.,2010.,1980.,1830.,1780.,
     * 1735.,1710.,1700.,1700.,1700.,  16*0./
C                                    PI+
      DATA (CROSEC(I,1,20,4),I=1,51)/
     *20.,40.,50.,60.,70.,80.,90.,100.,120.,140.,160.,180.,
     *200.,220.,250.,300.,350.,400.,450.,500.,550.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  15*0.,35./
      DATA (CROSEC(I,2,20,4),I=1,51)/
     *  485., 960.,1580.,2700.,3550.,4050.,4320.,4420.,4620.,4660.,
     * 4580.,4470.,4350.,4295.,4187.,3938.,3755.,3573.,3450.,3342.,
     * 3310.,3295.,3310.,3330.,3375.,3405.,3350.,3338.,3135.,3040.,
     * 2985.,2955.,2940.,2940.,2940.,  16*0./
      DATA (CROSEC(I,3,20,4),I=1,51)/
     *  361., 720.,1020.,1560.,2100.,2300.,2550.,2700.,2880.,2880.,
     * 2760.,2660.,2550.,2510.,2430.,2270.,2130.,2060.,2000.,1970.,
     * 1950.,1950.,1960.,1960.,1970.,1980.,1950.,1978.,1830.,1780.,
     * 1735.,1710.,1700.,1700.,1700.,  16*0./
C                                    ANTINEUTRON
      DATA (CROSEC(I,1,20,5),I=1,51)/
     *14.,15.,17.,19.,20.,22.,25.,27.,30.,35.,40.,45.,50.,55.,60.,70.,
     *80.,90.,100.,120.,140.,150.,160.,180.,200.,250.,300.,350.,
     *400.,500.,600.,700.,800.,900.,1000.,1500.,2000.,3000.,
     *5000.,7000.,10000.,20000.,50000.,100000.,500000.,1000000.,
     *4*0.,46./
      DATA (CROSEC(I,2,20,5),I=1,51)/
     * 5800.,5940.,6160.,6345.,6360.,6350.,6170.,6020.,5760.,5350.,
     * 4990.,4800.,4710.,4690.,4760.,5040.,5190.,5200.,5080.,4600.,
     * 4120.,3920.,3720.,3420.,3240.,3150.,3160.,3180.,3210.,
     * 3240.,3280.,3350.,3390.,3435.,3480.,3560.,3585.,3580.,
     * 3540.,3500.,3470.,3410.,3335.,3335.,3335.,3335.,  5*0./
      DATA (CROSEC(I,3,20,5),I=1,51)/
     * 2820.,2770.,2700.,2660.,2645.,2620.,2580.,2550.,2515.,2450.,
     * 2390.,2320.,2260.,2225.,2200.,2140.,2080.,2060.,2040.,2000., 
     * 1980.,1965.,1960.,1930.,1920.,1890.,1905.,1920.,
     * 1945.,1970.,1985.,2010.,2040.,2070.,2080.,2090.,2095.,2080.,
     * 2063.,2060.,2050.,2040.,2005.,2005.,2005.,2005.,  5*0./
C                                    ANTIPROTON
      DATA (CROSEC(I,1,20,6),I=1,51)/
     *14.,15.,17.,19.,20.,22.,25.,27.,30.,35.,40.,45.,50.,55.,60.,70.,
     *80.,90.,100.,120.,140.,150.,160.,180.,200.,250.,300.,350.,
     *400.,500.,600.,700.,800.,900.,1000.,1500.,2000.,3000.,
     *5000.,7000.,10000.,20000.,50000.,100000.,500000.,1000000.,
     *4*0.,46./
      DATA (CROSEC(I,2,20,6),I=1,51)/
     * 5800.,5940.,6160.,6345.,6360.,6350.,6170.,6020.,5760.,5350.,
     * 4990.,4800.,4710.,4690.,4760.,5040.,5190.,5200.,5080.,4600.,
     * 4120.,3920.,3720.,3420.,3240.,3150.,3160.,3180.,3210.,
     * 3240.,3280.,3350.,3390.,3435.,3480.,3560.,3585.,3580.,
     * 3540.,3500.,3470.,3410.,3335.,3335.,3335.,3335.,  5*0./
      DATA (CROSEC(I,3,20,6),I=1,51)/
     *  800., 900.,1100.,1300.,1410.,1510.,1680.,1800.,2000.,2200.,
     * 2080.,2060.,2025.,2100.,2030.,2030.,2000.,1960.,1960.,1960.,
     * 1940.,1925.,1920.,1905.,1890.,1860.,1880.,1910.,1930.,
     * 1945.,1985.,2010.,2040.,2070.,2080.,2090.,2095.,2080.,
     * 2063.,2060.,2050.,2040.,2005.,2005.,2005.,2005.,  5*0./
C                                    K-
      DATA (CROSEC(I,1,20,7),I=1,51)/
     *20.,40.,50.,60.,70.,80.,90.,100.,120.,140.,160.,180.,
     *200.,220.,250.,300.,350.,400.,450.,500.,550.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  15*0.,35./
      DATA (CROSEC(I,2,20,7),I=1,51)/
     * 7080.,6830.,6650.,6530.,6400.,6280.,6100.,5840.,5660.,5520.,
     * 5330.,5160.,4990.,4810.,4630.,4323.,4130.,3870.,3700.,3550.,
     * 3490.,3465.,3467.,3475.,3495.,3515.,3440.,3360.,3150.,3040.,
     * 2985.,2955.,2940.,2940.,2940.,  16*0./
      DATA (CROSEC(I,3,20,7),I=1,51)/
     * 1740.,2220.,2500.,2820.,3080.,3300.,3420.,3500.,3420.,3330.,
     * 3200.,3060.,2940.,2850.,2710.,2470.,2380.,2250.,2160.,2080.,
     * 2040.,2045.,2047.,2050.,2055.,2060.,2010.,1980.,1830.,1780.,
     * 1735.,1710.,1700.,1700.,1700.,  16*0./
C                                    K+
      DATA (CROSEC(I,1,20,8),I=1,51)/
     *20.,40.,50.,60.,70.,80.,90.,100.,120.,140.,160.,180.,
     *200.,220.,250.,300.,350.,400.,450.,500.,550.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  15*0.,35./
      DATA (CROSEC(I,2,20,8),I=1,51)/
     *  485., 960.,1580.,2700.,3550.,4050.,4320.,4420.,4620.,4660.,
     * 4580.,4470.,4350.,4295.,4187.,3938.,3755.,3573.,3450.,3342.,
     * 3310.,3295.,3310.,3330.,3375.,3405.,3350.,3338.,3135.,3040.,
     * 2985.,2955.,2940.,2940.,2940.,  16*0./
      DATA (CROSEC(I,3,20,8),I=1,51)/
     *  361., 720.,1020.,1560.,2100.,2300.,2550.,2700.,2880.,2880.,
     * 2760.,2660.,2550.,2510.,2430.,2270.,2130.,2060.,2000.,1970.,
     * 1950.,1950.,1960.,1960.,1970.,1980.,1950.,1978.,1830.,1780.,
     * 1735.,1710.,1700.,1700.,1700.,  16*0./
C                                    K0
      DATA (CROSEC(I,1,20,9),I=1,51)/
     *20.,40.,50.,60.,70.,80.,90.,100.,120.,140.,160.,180.,
     *200.,220.,250.,300.,350.,400.,450.,500.,550.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  15*0.,35./
      DATA (CROSEC(I,2,20,9),I=1,51)/
     *  485., 960.,1580.,2700.,3550.,4050.,4320.,4420.,4620.,4660.,
     * 4580.,4470.,4350.,4295.,4187.,3938.,3755.,3573.,3450.,3342.,
     * 3310.,3295.,3310.,3330.,3375.,3405.,3350.,3338.,3135.,3040.,
     * 2985.,2955.,2940.,2940.,2940.,  16*0./
      DATA (CROSEC(I,3,20,9),I=1,51)/
     *  361., 720.,1020.,1560.,2100.,2300.,2550.,2700.,2880.,2880.,
     * 2760.,2660.,2550.,2510.,2430.,2270.,2130.,2060.,2000.,1970.,
     * 1950.,1950.,1960.,1960.,1970.,1980.,1950.,1978.,1830.,1780.,
     * 1735.,1710.,1700.,1700.,1700.,  16*0./
C                              ANTI K0
      DATA (CROSEC(I,1,20,10),I=1,51)/
     *20.,40.,50.,60.,70.,80.,90.,100.,120.,140.,160.,180.,
     *200.,220.,250.,300.,350.,400.,450.,500.,550.,600.,
     *700.,800.,900.,1000.,2000.,3000.,5000.,10000.,20000.,50000.,
     *100000.,500000.,1000000.,  15*0.,35./
      DATA (CROSEC(I,2,20,10),I=1,51)/
     * 7080.,6830.,6650.,6530.,6400.,6280.,6100.,5840.,5660.,5520.,
     * 5330.,5160.,4990.,4810.,4630.,4323.,4130.,3870.,3700.,3550.,
     * 3490.,3465.,3467.,3475.,3495.,3515.,3440.,3360.,3150.,3040.,
     * 2985.,2955.,2940.,2940.,2940.,  16*0./
      DATA (CROSEC(I,3,20,10),I=1,51)/
     * 1740.,2220.,2500.,2820.,3080.,3300.,3420.,3500.,3420.,3330.,
     * 3200.,3060.,2940.,2850.,2710.,2470.,2380.,2250.,2160.,2080.,
     * 2040.,2045.,2047.,2050.,2055.,2060.,2010.,1980.,1830.,1780.,
     * 1735.,1710.,1700.,1700.,1700.,  16*0./
C
C ************************ END OF DATA ******************************
C
      END



      FUNCTION SIGION(ISS,A1,Z1,A2,Z2,T)
C**********************************************************
C     WRITTEN BY V.S.BARASHENKOV AND A.POLANSKI-JINR-DUBNA
C**********************************************************
C     FOR CALCULATION OF NUCLEUS-NUCLEUS TOTAL(ISS=1)
C     AND INELASTIC(ISS=2) CROSS SECTIONS
C     A1,Z1 - PROJECTILE MASS AND CHARGE NUMBERS (A1>1)
C     A2,Z2 - THE SAME FOR TARGET NUCLEUS(3<A2<240)
C     T - LAB. KINETIC ENERGY OF PROGECTALE (1 MEV/NUCLEON< T <1 TEV/NUCLEON)
C
      COMMON /CX/ CX(38)
      COMMON /FH/ AMP,AMT,AP,AT,B0,R0
C
      IS=3-ISS
      IF(IS.LE.0)RETURN
      IF(A1.lt.1. .OR. A1.gt.240. .OR. A2.lt.3. .OR. A2.gt.240.)
     *GOTO 101
      IF(ABS(Z1).LT.1.0)GOTO 101
      IF(T.LT.1.0)GOTO 101
c Sobolevsky, 13.07.97 --------------------
c     If a projectile is heavier then a target, the diagnostic
c     "ERROR IN INPUT OF PARAMETERS OF FUNCTION SIGION" appears
c     and/or results are wrong. Therefore let us, in this case,
c     to make the PERMUTATION  Projectile <----> Target.
        PERMUT2=0.
      IF(A1.gt.A2)THEN
        PERMUT2=1.
C:::::::::::: 17.08.2005 ::::::::::::::::
        T0=T         ! Saving of input energy T(MeV) of projectile (A1,Z1)
        T=(T/A1)*A2  ! Kinetic energy (MeV) of target (A2,Z2) in antilab frame
C:::::::::::: end 17.08.2005 ::::::::::::
          A0=A1
          Z0=Z1
            A1=A2
            Z1=Z2
              A2=A0
              Z2=Z0
      END IF
c ---------------- end of Sobol -------------
      SIGION=0.
      TP=T/A1
      AP=A1**0.333333
      AT=A2**0.333333
      AMP=A1*930.63
      AMT=A2*930.63
C     PARAMETER FOR CALCULATION OF NUCLEAR RADIUS
      R0=1.4
      IF(ABS(A1-4.) .LT. 0.1) R0=1.3
      B0=1.44*Z1*Z2
      I=1
       IF(IS.EQ.2) I=20
C     SELECTION OF PROGECTALES
C     HEVY ION
      IF(A1.GT.4.1) N=I
C     ALFA,HELION,TRITON
      IF(A1.GT.2.1 .AND. A1.LT.4.1) N=I+6
C     DEUTRON
      IF(A1.LT.2.1) N=I+12
C     HIGH-ENERGY CROSS-SECTION
C     SELECTION OF PROGECTALE ENERGY
      IF(TP.LT.CX(N+1)) K=2
      IF(TP.LT.CX(N+4)) K=5
C     CROSS-SACTION PARAMRTERS
      C=CX(I)
      IF(TP.LT.CX(N+1)) C=CX(N+K)+CX(N+K+1)*LOG10(TP)
      CP=CX(N+5)+CX(N+6)
      IF(TP.LT.10.)GOTO 1
C     HIGH-ENERGY CROSS-SECTION
      SIGION=FHS(IS,T,C)
       if(SIGION.le.0.)SIGION=0.   ! 10.08.2005
      RETURN
C     CALCULATION OF LOW-ENERGY CROSS-SECTION
C     NORMALUSED HIGH-ENERGY CROSS-SECTION
 1     SH10=FHS(IS,10.*A1,CP)
      R0=1.45
      IF(ABS(A1-4.).LT. 0.1) R0=1.4
C     RENORMALUSED COULOMB BARRIER
      B=B0/R0/(AP+AT)
C     LOW-ENERGY CROSS-SECTION
      SIGION=SH10*FC(T,B)/FC(10.*A1,B)
C      IF(SIGION) 101,100,100   ! commented of 10.08.2005
C============== 10.08.2005 ===========
      if(SIGION.le.0.)then
        SIGION=0.
        goto 100
      end if
C============== end 10.08.2005 =======
C  101 WRITE(6,1001)   ! commented of 10.08.2005
  101 continue   ! 10.08.2005
Cc Sobolevsky, 13.07.97 -------- ! write/stop are commented of 10.08.2005
C      CONTINUE!WRITE(25,1001)
C      stop
Cc -----------------------------
 1001 FORMAT(' ERROR IN INPUT OF PARAMETERS OF FUNCTION SIGION')
  100 CONTINUE
c Sobolevsky, 13.07.97 --------------------
c     Inverse PERMUTATION to restore input parameters
      IF(PERMUT2.gt.0.5)THEN
        T=T0     ! 17.08.2005
          A0=A1
          Z0=Z1
            A1=A2
            Z1=Z2
              A2=A0
              Z2=Z0
      END IF
c ---------------- end of Sobol -------------
      RETURN
      END
      
      
      FUNCTION FC(T,B)
      COMMON /FH/AMP,AMT,AP,AT,B0,R0
C     CMS ENERGY
      TC=T*AMT/(AMP+AMT)
      X=(TC-B)/1.2
      IF(X.GT.5)GOTO 1
      D=1.+EXP(X)
      FC=ALOG(D)/TC
      RETURN
  1   FC=X/TC
      RETURN
      END
      
      
      FUNCTION FHS(IS,E,C)
C     CALCULATION OF HIGH-ENERGY TOTAL (IS=2) AND
C     INELASTIC (IS=1) CROSS-SECTIONS

C      E - LAB. KINETIC ENERGY OF PROGECTALE(MEV)
      COMMON/FH/AMP,AMT,AP,AT,B0,R0
C     SQUED PROGECTALE CMS MOMENTUM
      PPC=AMT*AMT*E*(E+2.*AMP)/((AMP+AMT)**2+2.*AMT*E)
C     DE BROGLE WAVE LANGTH
      AL=1.41*140./SQRT(PPC)
      EC=SQRT(PPC+AMP*AMP)-AMP
C     COULOMB BARRIER
      B=B0/R0/(AP+AT+AL)
      FHS=31.416*1.21*(1.-B/EC)*(AP+AT+1.85*AP*AT/(AP+AT)
     *+AL-C)**2*IS
      RETURN
      END
      
      
      BLOCK DATA C
C     NUCLEUS-NUCLEUS CROSS-SECTION PARAMETERS

      COMMON/CX/CX(38)
      DATA  CX/2.07,560.,0.8,0.426,
     *           100.,-2.05,1.9,
     *           200.,0.07,0.87,
     *           20.,-1.55,2.1,
     *           700.,-1.01,1.08,
     *           400.,-0.59,0.94,
     *        2.45,225.,-2.25,2.,
     *           100.,-4.61,3.18,
     *           185.,-3.,2.4,
     *           185.,-3.,2.4,
     *           185.,-4.77,3.18,
     *           185.,-4.77,3.18/
      END


      BLOCK DATA A
      COMMON /TAPEX/ AAI(19),IENER(94),
     * ISIG1(969),ISIG2(969),ISIG3(817),ISIG4(817)
C     AAI- TARGET NUCLEUS MASS NUMBERS
      DATA AAI/
     *  4.00,   9.01,  12.00,  14.00,  16.00,  23.00,  26.98,  40.08,
     *  47.90, 55.85,  63.55,  79.90,  95.94, 112.40, 118.69, 137.34,
     * 183.85, 207.19,238.03/
C     IENER- PROJECTILE KINETIC ENERGIES(MEV)
      DATA IENER/
     *14,15,16,17,18,19,20,22,25,27,30,33,35,37,40,45,50,55,60,65,70,
     *80,90,100,120,140,150,160,180,200,250,300,350,400,500,600,700,800,
     *900,1000,1500,2000,3000,5000,7000,10000,20000,50000,100000,
     *500000,1000000,
     *20,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,
     *220,240,250,260,280,300,350,400,450,500,550,600,700,800,900,1000,
     *2000,3000,5000,10000,20000,50000,100000,500000,1000000/
C      ISIG1- NEUTRON-NUCLEUS CROSS SECTIONS
C      FOUR THE LAST DECIMAL POSITION IN THE NUMBER-INELASTIC CROSS-SECTION.
C      FOUR THE FIRST ONE-TOTAL CROSS SECTION.
      DATA (ISIG1(J),J=1,510)/
     *10900000,10200005, 9670007, 9150010, 8760013, 8380016, 8000020, 
     * 7100035,
     * 6400055, 6000070, 5600080, 5400085, 5000090, 4700097, 4400105, 
     * 3900115,
     * 3600115, 3270107, 2950100, 2750095, 2560090, 2200086, 1920084, 
     * 1680084,
     * 1360082, 1200080, 1160080, 1140080, 1100080, 1070079, 1040078, 
     * 1060080,
     * 1080084, 1100088, 1200094, 1260100, 1350105, 1400108, 1440108, 
     * 1460108,
     * 1480112, 1520114, 1500114, 1460112, 1420110, 1380108, 1320106, 
     * 1290104,
     * 1260101, 1270102, 1280102,
c       2,       1,       0,       0,       0,
C     THE COMMENT c (HERE AND FUTHER) ENDS THE CROSS-SECTION ARRAY.
C     THE FIRST NUMBER (IN COMMENT c) IS TARGET CHARGE NUMBER,
C     THE SECOND NUMBER IS THE TYPE OF CROSS-SECTIONS (1-NEUTRONS,
C     2-PROTONS,3-MESONS -PI,4-MESONS +PI)
     *14900650,14600640,14300628,14000617,13830610,13660603,13500595,
     * 12700555,
     *12000520,11600495,11000470,10500450,10000430, 9550407, 9100385, 
     * 8100350,
     * 7400320, 6820285, 6250270, 6000260, 5750250, 4550210, 4060190, 
     * 3650185,
     * 3100178, 2750175, 2620175, 2550175, 2400175, 2350175, 2250170, 
     * 2250170,
     * 2300172, 2380176, 2520184, 2700194, 2820200, 2880209, 2900213, 
     * 2940214,
     * 3030216, 3030216, 3000212, 2920210, 2840210, 2770210, 2670210, 
     * 2630210,
     * 2640210, 2680210, 2680210,
c       4,       1,       0,       0,       0,
     *12400590,13700570,14100556,14500542,14520531,14540520,14550510,
     * 14450500,
     *13850460,13450445,12900430,12500412,12100395,11600387,11100380,
     * 10200350,
     * 9400330, 8700312, 8000295, 7500282, 7000270, 6040255, 5300240, 
     * 4750228,
     * 3960222, 3500216, 3360216, 3200210, 3030210, 2940210, 2800208, 
     * 2800210,
     * 2860214, 2960216, 3140228, 3300240, 3440248, 3560254, 3600257, 
     * 3640260,
     * 3840262, 3880260, 3840256, 3640252, 3520252, 3440250, 3300250, 
     * 3240248,
     * 3240248, 3320248, 3320248,
c       6,       1,       0,       0,       0,
     *14200680,14800665,15080645,15370625,15410610,15460595,15500580,
     * 15250562,
     *15000525,14800510,14250485,13820467,13400450,13000442,12600435,
     * 11750410,
     *10900387,10100364, 9300340, 8670325, 8050310, 6900290, 6120280, 
     * 5520276,
     * 4620274, 4020260, 3840258, 3720252, 3500247, 3450245, 3260240, 
     * 3240240,
     * 3280244, 3360250, 3560260, 3720268, 3880275, 4000280, 4080285, 
     * 4150290,
     * 4300295, 4350300, 4320294, 4150292, 4020290, 3900285, 3750285, 
     * 3670282,
     * 3700282, 3820282, 3850282,
c       7,       1,       0,       0,       0,
     *15200750,15700740,16000720,16300700,16400683,16500667,16600650,
     * 16470620,
     *16230575,15950555,15550530,15150517,14750505,14350484,13950462,
     * 12900435,
     *12070420,11210397,10350375, 9750360, 9250345, 8160320, 7200310, 
     * 6450300,
     * 5400293, 4620288, 4380282, 4150282, 3920280, 3780276, 3620270, 
     * 3610271,
     * 3810275, 3900280, 4030290, 4170295, 4400304, 4600310, 4700315, 
     * 4790318,
     * 4980332, 5040335, 4980330, 4770323, 4570320, 4430317, 4270315, 
     * 4200315,
     * 4250315, 4290315, 4300315,
c       8,       1,       0,       0,       0,
     *15700960,16200930,16550910,16900890,17030867,17160845,17300822,
     * 17500790,
     *17600750,17550725,17400686,17250653,17100620,16760610,16430600,
     * 15600575,
     *14800540,14110518,13430497,12810473,12200450,10730414, 9530390, 
     * 8600380,
     * 7200372, 6180364, 5820360, 5460355, 5220354, 5040350, 4840350, 
     * 4920350,
     * 5000356, 5120364, 5380384, 5600392, 5860400, 6080408, 6220410, 
     * 6320420,
     * 6600408, 6680412, 6640420, 6400414, 6160409, 5960407, 5680403, 
     * 5680400,
     * 5680400, 5680400, 5680400,
c      11,       1,       0,       0,       0,
     *17351000,17500990,17600975,17950950,18070935,18190920,18300905,
     * 18550875,
     *18850825,18960800,19000762,18850726,18700690,18520671,18350652,
     * 17850610,
     *17100570,16160532,15220495,14380487,13500480,12120456,10800444, 
     * 9720432,
     * 8160420, 7200420, 6780420, 6420420, 6000410, 5670410, 5580400, 
     * 5600402,
     * 5780404, 5920408, 6160424, 6440438, 6720448, 6880450, 7080454, 
     * 7200456,
     * 7360472, 7540480, 7360466, 7060456, 6800452, 6720448, 6460444, 
     * 6320440,
     * 6320440, 6320440, 6320440,
c      13,       1,       0,       0,       0,
     *18101060,18331035,18501015,18720990,18880972,19040954,19200935,
     * 19500900,
     *19950860,20200830,20350790,20170758,20000725,19650695,19300665,
     * 18500630,
     *17600600,16650560,15700520,14850512,14000504,12550486,11100470,
     * 10080456,
     * 8460444, 7420432, 6960432, 6710432, 6230418, 5880418, 5840415, 
     * 5840412,
     * 6020416, 6180422, 6450440, 6790460, 7080472, 7270476, 7460479, 
     * 7570480,
     * 7690492, 7820492, 7710488, 7340472, 7100472, 6980464, 6720460, 
     * 6550452,
     * 6500448, 6500448, 6500448,
c      20,       1,       0,       0,       0,
     *21801240,21301225,20951200,20751180,20621162,20751144,21151125,
     * 21501090,
     *22201045,22501020,23000980,23320952,23650925,23670902,23600880,
     * 22800825,
     *21800770,20900725,20000680,19020660,18050640,16500620,15000615,
     * 13400600,
     *11400580, 9900565, 9400560, 8900560, 8250560, 7900550, 7700536, 
     * 7730530,
     * 7870540, 8000550, 8300570, 8700595, 9050610, 9300615, 9500620, 
     * 9650622,
     * 9900629,10020630, 9900620, 9650612, 9450607, 9250592, 8920587, 
     * 8600580,
     * 8600580, 8600580, 8600580,
c      22,       1,       0,       0,       0,
     *23401270,23001255,22381242,22051235,21971217,21891199,21801180,
     * 22001160,
     *22351180,22671090,23181050,23491027,23801005,24150992,24500980,
     * 24000930,
     *23300880,22450830,21700800,20950780,20200760,18600730,16950715,
     * 15400705,
     *13100685,11500680,10850665,10450660, 9700645, 9400630, 9100620, 
     * 9150620,
     * 9270620, 9330625, 9750650,10100670,10300680,10550687,10700695,
     * 11900700,
     *11200710,11330720,11220712,11000705,10800695,10650690,10400680,
     * 10200675,
     *10120670,10100670,10100670/
c      26,       1,       0,       0,       0,
      DATA (ISIG1(J),J=511,969)/
     *25801440,24901433,24301411,23701390,23411368,23111347,22821325,
     * 22751280,
     *22851260,23201215,23701180,24071145,24451110,24651075,24851040,
     * 25300990,
     *25400955,25170920,24800885,23850860,22900835,21100800,19400780,
     * 17900765,
     *15100750,12900725,12200720,11500720,10700710,10300700,10130700,
     * 10200700,
     *10300712,10430724,10750735,11100750,11330765,11630775,11850780,
     * 12250795,
     *12520810,12600813,12600810,12330784,12070757,11850743,11400735,
     * 11100720,
     *11100720,11100720,11100720,
c      29,       1,       0,       0,       0,
     *29201540,28001535,27071517,26151500,25251463,25701481,24801445,
     * 24551407,
     *24301380,24401330,24601300,24951285,25301270,25721230,26151190,
     * 26901090,
     *27201010,27000940,26450920,25720890,25000860,23200835,21400820,
     * 19700810,
     *16700800,14600780,13800775,12850770,12000760,11600760,11400758,
     * 11470765,
     *11630765,11700770,12000795,12370810,12650825,12850830,13050840,
     * 13280848,
     *13750870,13930872,13900870,13600845,13350825,13150810,12700803,
     * 12300795,
     *12300975,12300975,12300795,
c      35,       1,       0,       0,       0,
     *35501670,33701650,32451625,31201600,30521580,29831560,29151540,
     * 28501500,
     *27421460,27001430,27001400,27201360,27271340,28001270,29301200,
     * 30001145,
     *30201100,29701060,28501027,28121023,27751020,26700990,25050975,
     * 23250960,
     *20000935,17700915,16500900,15900900,14600880,14100870,14000865,
     * 14100865,
     *14300872,14650890,15000920,15200950,15400960,15500970,15600975,
     * 16380965,
     *16501000,16351005,16001000,15800990,15700975,15300965,15000960,
     * 14900955,
     *14850955,14850955,14850955,
c      42,       1,       0,       0,       0,
     *41501790,40401775,39201757,38001740,36971720,35931700,34901680,
     * 33001640,
     *30601580,29601550,28451510,27851460,28201440,28501418,29801380,
     * 31701330,
     *32301280,32551240,32801200,32521177,32251155,30751140,28951110,
     * 27101110,
     *23551080,20601065,19251050,18001050,16301025,15601020,15401015,
     * 15501020,
     *15701022,15901026,16501060,16851085,17151100,17401110,17601120,
     * 17801127,
     *18501150,18801160,18581140,18151100,17901085,17821080,17201070,
     * 16901070,
     *16901070,16901070,16901070,
c      48,       1,       0,       0,       0,
     *44201920,42801910,42251895,41701880,40701860,39651850,38601840,
     * 36801800,
     *34201760,32801720,31251675,30601630,30801600,31351500,31901520,
     * 33501465,
     *34451420,35101390,35401340,35601310,35501280,34601275,33001235,
     * 30301225,
     *26401200,23401170,21901170,20701170,19501165,17701145,17321140,
     * 17401140,
     *17601135,17801160,18321180,18851220,19251240,19451250,19601260,
     * 19801265,
     *20701270,20801275,20651250,20401222,20221222,19801220,19401215,
     * 18701190,
     *18701190,18701190,18701190,
c      50,       1,       0,       0,       0,
     *49802160,50002150,49352135,48702120,47202100,45852080,44502060,
     * 42002025,
     *39301970,37901935,36651880,35601830,35001800,35141762,35251725,
     * 35701670,
     *36401610,37201575,37801540,38501510,39201490,38001460,36701440,
     * 35001400,
     *30801385,27301375,26001360,24801350,22801340,21401345,20751325,
     * 20901325,
     *21001330,21251345,21801370,22301400,22651420,22851430,23151435,
     * 23251440,
     *23951455,24251465,24001445,23601440,23451435,23251430,22851420,
     * 22301405,
     *22251400,22251400,22251400,
c      56,       1,       0,       0,       0,
     *53202440,54302400,54552385,54802370,54502350,53902330,53302310,
     * 51902270,
     *49602220,47902195,45502150,43402100,42002070,41352040,40702010,
     * 40001945,
     *40301900,41251850,42201820,43051790,43901760,44401730,43601720,
     * 42001680,
     *38001680,33801660,32001660,30401650,27901650,26601640,25751640,
     * 25751612,
     *26001615,26401625,26901640,27551700,27901720,28121730,28371740,
     * 28501750,
     *29501780,30001780,29701750,29401740,29101735,28801710,28201695,
     * 27301680,
     *27301680,27301680,27301680,
c      74,       1,       0,       0,       0,
     *53002580,54402550,55802527,57202505,57002483,58802462,57652460,
     * 57452435,
     *54802380,52802355,49702280,47602230,45502180,44702175,43902170,
     * 43002130,
     *42652080,43252035,44501980,44951960,45401940,47401900,47101870,
     * 46001840,
     *41001800,36601800,34801800,33001780,30001760,28901760,28651740,
     * 28551730,
     *28551730,28651740,29201785,29551815,30001835,30301860,30601890,
     * 31051895,
     *32401920,32901920,32701890,32401850,31801835,30901830,30601830,
     * 29701830,
     *29701830,29701830,29701830,
c      82,       1,       0,       0,       0,
     *58002820,59402770,60502735,61602700,62522680,63452660,63602645,
     * 63502620,
     *61702580,60202550,57602515,55552482,53502450,51702420,49902390,
     * 48002320,
     *47102260,46902225,47602200,49002170,50402140,51902080,52002060,
     * 50802040,
     *46002000,41201980,39201965,37201960,34201930,32401920,31501890,
     * 31601905,
     *31801920,32101945,32401970,32801985,33502010,33902040,34352070,
     * 34802080,
     *35602090,35852095,35802080,35402063,35002060,34702050,34102040,
     * 33352005,
     *33352005,33352005,33352005/
c      92,       1,       0,       0,       0,
C      ISIG2- PROTON-NUCLEUS CROSS SECTIONS
C      FOUR THE LAST DECIMAL POSITION IN THE NUMBER-INELASTIC CROSS-SECTION.
C      FOUR THE FIRST ONE-TOTAL CROSS SECTION.
      DATA (ISIG2(J),J=1,510)/
     *10900000,10200002, 9670002, 9150003, 8760006, 8380009, 8000013, 
     * 7100030,
     * 6400050, 6000065, 5600077, 5400083, 5000090, 4700097, 4400105, 
     * 3900115,
     * 3600115, 3270107, 2950100, 2750095, 2560090, 2200086, 1920084, 
     * 1680084,
     * 1360082, 1200080, 1160080, 1140080, 1100080, 1070079, 1040078, 
     * 1060080,
     * 1080084, 1100088, 1200094, 1260100, 1350105, 1400108, 1440108, 
     * 1460108,
     * 1480112, 1520114, 1500114, 1460112, 1420110, 1380108, 1320106, 
     * 1290104,
     * 1260101, 1270102, 1280102,
c       2,       2,       0,       0,       0,
C     THE COMMENT c (HERE AND FUTHER) ENDS THE CROSS-SECTION ARRAY.
C     THE FIRST NUMBER (IN COMMENT c) IS TARGET CHARGE NUMBER,
C     THE SECOND NUMBER IS THE TYPE OF CROSS-SECTIONS (1-NEUTRONS,
C     2-PROTONS,3-MESONS -PI,4-MESONS +PI)
     *14900490,14600540,14300560,14000580,13830568,13660556,13500545,
     * 12700525,
     *12000495,11600470,11000450,10500435,10000420, 9550395, 9100370, 
     * 8100340,
     * 7400310, 6820286, 6250262, 6000252, 5750242, 4550205, 4060185, 
     * 3650180,
     * 3100175, 2750172, 2620175, 2550175, 2400175, 2350175, 2250170, 
     * 2250170,
     * 2300172, 2380176, 2520184, 2700194, 2820200, 2880209, 2900213, 
     * 2940214,
     * 3030216, 3030216, 3000212, 2920210, 2840210, 2770210, 2670210, 
     * 2630210,
     * 2640210, 2680210, 2680210,
c       4,       2,       0,       0,       0,
     *12400310,13700330,14100365,14500400,14520413,14540427,14550440,
     * 14450450,
     *13850435,13450430,12900420,12500402,12100385,11600377,11100370,
     * 10200340,
     * 9400320, 8700304, 8000288, 7500275, 7000263, 6040249, 5300234, 
     * 4750222,
     * 3960216, 3500210, 3360211, 3200205, 3030208, 2940210, 2800208, 
     * 2800210,
     * 2860214, 2960216, 3140228, 3300240, 3440248, 3560254, 3600257, 
     * 3640260,
     * 3840262, 3880260, 3840256, 3640252, 3520252, 3440250, 3300250, 
     * 3240248,
     * 3240248, 3320248, 3320248,
c       6,       2,       0,       0,       0,
     *14200420,14800440,15080455,15370470,15410476,15460483,15500490,
     * 15250497,
     *15000500,14800480,14250462,13820451,13400440,13000432,12600425,
     * 11750400,
     *10900377,10100355, 9300333, 8670318, 8050303, 6900248, 6120274, 
     * 5520270,
     * 4620268, 4020254, 3840252, 3720247, 3500245, 3450245, 3260240, 
     * 3240240,
     * 3280244, 3360250, 3560260, 3720268, 3880275, 4000280, 4080285, 
     * 4150290,
     * 4300295, 4350300, 4320294, 4150292, 4020290, 3900285, 3750285, 
     * 3670282,
     * 3700282, 3820282, 3850282,
c       7,       2,       0,       0,       0,
     *15200460,15700485,16000497,16300510,16400518,16500526,16600535,
     * 16470537,
     *16230532,15950520,15550500,15150480,14750460,14350446,13950432,
     * 12900405,
     *12070390,11210370,10350350, 9750335, 9250320, 8160310, 7200304, 
     * 6450293,
     * 5400287, 4620283, 4380279, 4150279, 3920278, 3780276, 3620270, 
     * 3610271,
     * 3810275, 3900280, 4030290, 4170295, 4400304, 4600310, 4700315, 
     * 4790318,
     * 4980332, 5040335, 4980330, 4770323, 4570320, 4430317, 4270315, 
     * 4200315,
     * 4250315, 4290315, 4300315,
c       8,       2,       0,       0,       0,
     *15700600,16200617,16550639,16900660,17030665,17160670,17300675,
     * 17500680,
     *17600680,17550670,17400650,17250612,17100575,16760562,16430550,
     * 15600525,
     *14800490,14110470,13430450,12810435,12200420,10730385, 9530367, 
     * 8600360,
     * 7200350, 6180350, 5820350, 5460345, 5220347, 5040350, 4840350, 
     * 4920350,
     * 5000356, 5120364, 5380384, 5600392, 5860400, 6080408, 6220410, 
     * 6320420,
     * 6600408, 6680412, 6640420, 6400414, 6160409, 5960407, 5680403, 
     * 5680400,
     * 5680400, 5680400, 5680400,
c      11,       2,       0,       0,       0,
     *17350650,17500682,17600690,17950715,18070727,18190739,18300750,
     * 18550762,
     *18850750,18960740,19000720,18850687,18700655,18520636,18350617,
     * 17850575,
     *17100540,16160505,15220470,14380462,13500455,12120432,10800420, 
     * 9720408,
     * 8160400, 7200403, 6780403, 6420408, 6000406, 5670404, 5580400, 
     * 5600402,
     * 5780404, 5920408, 6160424, 6440438, 6720448, 6880450, 7080454, 
     * 7200456,
     * 7360472, 7540480, 7360466, 7060456, 6800452, 6720448, 6460444, 
     * 6320440,
     * 6320440, 6320440, 6320440,
c      13,       2,       0,       0,       0,
     *18100670,18330700,18500725,18720750,18880760,19040770,19200780,
     * 19500780,
     *19950770,20200757,20350735,20170712,20000690,19650657,19300635,
     * 18500585,
     *17600570,16650430,15700490,14850483,14000475,12550460,11100446,
     * 10080431,
     * 8460423, 7420425, 6960425, 6710425, 6230425, 5880422, 5840422, 
     * 5840412,
     * 6020416, 6180422, 6450440, 6790460, 7080472, 7270476, 7460479, 
     * 7570480,
     * 7690492, 7820492, 7710488, 7340472, 7100472, 6980464, 6720460, 
     * 6550452,
     * 6500448, 6500448, 6500448,
c      20,       2,       0,       0,       0,
     *21800770,21300800,20950823,20750850,20620866,20750883,21150900,
     * 21500925,
     *22200935,22500920,23000895,23320865,23650835,23670817,23600800,
     * 22800750,
     *21800715,20900677,20000640,19020622,18050605,16500590,15000588,
     * 13400573,
     *11400553, 9900540, 9400535, 8900540, 8250543, 7900535, 7700520, 
     * 7730517,
     * 7870530, 8000542, 8300568, 8700595, 9050610, 9300615, 9500620, 
     * 9650622,
     * 9900629,10020630, 9900620, 9650612, 9450607, 9250592, 8920582, 
     * 8600580,
     * 8600580, 8600580, 8600580,
c      22,       2,       0,       0,       0,
     *23400840,23000915,22380950,22050980,21970996,21891012,21801028,
     * 22001030,
     *22351020,22671005,23180975,23490947,23800920,24150902,24500895,
     * 24000860,
     *23300825,22450790,21700760,20950740,20200720,18600690,16950685,
     * 15400680,
     *13100680,11500645,10850640,10450625, 9700615, 9400600, 9100600, 
     * 9150600,
     * 9270602, 9330605, 9750643,10100670,10300680,10550687,10700695,
     * 11900700,
     *11200710,11330720,11220712,11000705,10800695,10650690,10400680,
     * 10200675,
     *10120670,10100670,10100670/
c      26,       2,       0,       0,       0,
      DATA (ISIG2(J),J=511,969)/
     *25800900,24900960,24301015,23701070,23411077,23111083,22821090,
     * 22751115,
     *22851120,23201115,23701080,24071052,24451025,24650992,24850960,
     * 25300900,
     *25400885,25170865,24800795,23850782,22900770,21100740,19400720,
     * 17900710,
     *15100700,12900675,12200680,11500682,10700683,10300685,10130685,
     * 10200688,
     *10300700,10430712,10750725,11100742,11330757,11630770,11850780,
     * 12250795,
     *12520810,12600813,12600810,12330784,12070757,11850743,11400735,
     * 11100720,
     *11100720,11100720,11100720,
c      29,       2,       0,       0,       0,
     *29200935,28001000,27071030,26151060,25251103,25701146,24801190,
     * 24551220,
     *24301250,24401240,24601210,24951170,25301130,25721095,26151060,
     * 26900980,
     *27200930,27000870,26450860,25720830,25000800,23200765,21400760,
     * 19700750,
     *16700750,14600730,13800725,12850720,12000720,11600730,11400733,
     * 11470745,
     *11630747,11700755,12000780,12370800,12650818,12850825,13050835,
     * 13280845,
     *13750870,13930872,13900870,13600845,13350825,13150810,12700803,
     * 12300795,
     *12300975,12300975,12300795,
c      35,       2,       0,       0,       0,
     *35501000,33701070,32451120,31201170,30521210,29831250,29151290,
     * 28501330,
     *27421340,27001330,27001305,27201260,27271240,28001165,29301110,
     * 30001060,
     *30201015,29700970,28500935,28120933,27750932,26700930,25050915,
     * 23250900,
     *20000880,17700870,16500855,15900855,14600840,14100830,14000830,
     * 14100840,
     *14300855,14650880,15000915,15200938,15400948,15500960,15600965,
     * 16380955,
     *16500995,16351005,16001000,15800990,15700975,15300965,15000960,
     * 14900955,
     *14850955,14850955,14850955,
c      42,       2,       0,       0,       0,
     *41501025,40401080,39201135,38001190,36971253,35931316,34901380,
     * 33001440,
     *30601470,29601450,28451410,27851360,28201340,28501318,29801280,
     * 31701230,
     *32301180,32551140,32801100,32521090,32251080,30751065,28951050,
     * 27101050,
     *23551010,20601005,19251000,18001000,16300982,15600980,15400975,
     * 15500985,
     *15701092,15901000,16501037,16851065,17151085,17401095,17601108,
     * 17801117,
     *18501145,18801160,18581140,18151100,17901085,17821080,17201070,
     * 16901070,
     *16901070,16901070,16901070,
c      48,       2,       0,       0,       0,
     *44201020,42801100,42251162,41701225,40701290,39651365,38601440,
     * 36801520,
     *34201575,32801560,31251518,30601460,30801420,31351410,31901400,
     * 33501365,
     *34451320,35101245,35401245,35601215,35501185,34601185,33001155,
     * 30301145,
     *26401120,23401100,21901100,20701110,19501105,17701190,17321185,
     * 17401190,
     *17601100,17801125,18321150,18851195,19251220,19451235,19601245,
     * 19801250,
     *20701260,20801270,20651250,20401222,20221222,19801220,19401215,
     * 18701190,
     *18701190,18701190,18701190,
c      50,       2,       0,       0,       0,
     *49801060,50001170,49351245,48701320,47201380,45851437,44501495,
     * 42001580,
     *39301660,37901680,36651665,35601635,35001610,35141585,35251560,
     * 35701535,
     *36401500,37201460,37801430,38501390,39201375,38001345,36701337,
     * 35001335,
     *30801315,27301305,26001320,24801310,22801285,21401280,20751275,
     * 20901280,
     *21001229,21251310,21801340,22301380,22651400,22851410,23151415,
     * 23251430,
     *23951445,24251458,24001440,23601440,23451435,23251430,22851420,
     * 22301405,
     *22251400,22251400,22251400,
c      56,       2,       0,       0,       0,
     *53200910,54301020,54551130,54801240,54501400,53901480,53301560,
     * 51901670,
     *49601760,47901830,45501850,43401870,42001870,41351855,40701840,
     * 40001780,
     *40301760,41251740,42201700,43051675,43901650,44401610,43601600,
     * 42001600,
     *38001570,33801563,32001560,30401550,27901550,26601545,25751547,
     * 25751555,
     *26001560,26401580,26901600,27551620,27901685,28121700,28371710,
     * 28501725,
     *29501760,30001767,29701745,29401740,29101735,28801710,28201695,
     * 27301680,
     *27301680,27301680,27301680,
c      74,       2,       0,       0,       0,
     *53000800,54401060,55801130,57201200,57001290,58801420,57651515,
     * 57451620,
     *54801750,52801800,49701915,47601972,45502030,44701995,43901960,
     * 43001910,
     *42651850,43251815,44501770,44951750,45401730,47401722,47101718,
     * 46001718,
     *41001680,36601680,34801675,33001675,30001650,28901640,28651640,
     * 28551645,
     *28551645,28651665,29201710,29551750,30001780,30301815,30601850,
     * 31051865,
     *32401900,32901905,32701885,32401850,31801835,30901830,30601830,
     * 29701830,
     *29701830,29701830,29701830,
c      82,       2,       0,       0,       0,
     *58000680,59400900,60501000,61601100,62521200,63451300,63601410,
     * 63501510,
     *61701680,60201650,57601900,55551950,53502200,51701985,49901970,
     * 48001930,
     *47101900,46901880,47601880,49001877,50401875,51901870,52001860,
     * 50801840,
     *46001840,41201815,39201800,37201790,34201790,32401780,31501790,
     * 31601815,
     *31801840,32101860,32401890,32801920,33501950,33901990,34352030,
     * 34802080,
     *35602070,35852078,35802070,35402063,35002060,34702050,34102040,
     * 33352005,
     *33352005,33352005,33352005/
c_____ 92,       2,       0,       0,       0,
C      ISIG3- PI- MESON-NUCLEUS CROSS SECTIONS
C      FOUR THE LAST DECIMAL POSITION IN THE NUMBER-INELASTIC CROSS-SECTION.
C      FOUR THE FIRST ONE-TOTAL CROSS SECTION.
      DATA (ISIG3(J),J=1,430)/
     *  400018,  700038,  890050, 1080062, 1300080, 1520098, 1800117, 
     * 2080136,
     * 2420156, 2760176, 2980188, 3200200, 3290209, 3330212, 3300210, 
     * 3280208,
     * 3220204, 3100196, 2880176, 2600164, 2500157, 2400150, 2160134, 
     * 1960124,
     * 1440097, 1250090, 1120085, 1080082, 1090083, 1100086, 1170093, 
     * 1230097,
     * 1280100, 1350102, 1100083,  960077,  870075,  850074,  830072,  
     * 830072,
     *  830072,  830072,  830072,
c       2,       3,       0,       0,       0,
C     THE COMMENT c (HERE AND FUTHER) ENDS THE CROSS-SECTION ARRAY.
C     THE FIRST NUMBER (IN COMMENT c) IS TARGET CHARGE NUMBER,
C     THE SECOND NUMBER IS THE TYPE OF CROSS-SECTIONS (1-NEUTRONS,
C     2-PROTONS,3-MESONS -PI,4-MESONS +PI)
     * 1500090, 2100126, 2520151, 2940177, 3450208, 3960240, 4580280, 
     * 5200320,
     * 5600350, 6000380, 6230400, 6350410, 6420414, 6400410, 6300400, 
     * 6150387,
     * 6000371, 5760360, 5400333, 5040312, 4870298, 4700285, 4350260, 
     * 4000237,
     * 3400216, 2940198, 2580187, 2360182, 2300180, 2330182, 2440187, 
     * 2570193,
     * 2700202, 2760207, 2500179, 2300172, 2150165, 2050159, 1940155, 
     * 1880154,
     * 1860154, 1860154, 1860154,
c       4,       3,       0,       0,       0,
     * 2040128, 2600160, 3130192, 3660224, 4410269, 5170315, 5730351, 
     * 6300388,
     * 6620409, 6940430, 7040438, 7100444, 7110445, 7060440, 6940432, 
     * 6760416,
     * 6620408, 6480400, 6160380, 5840354, 5660337, 5480320, 5180304, 
     * 4890288,
     * 4260264, 3760246, 3420240, 3230233, 3100232, 3130234, 3190238, 
     * 3330246,
     * 3420252, 3480256, 3100220, 2900210, 2680198, 2500187, 2450183, 
     * 2370176,
     * 2340174, 2340174, 2340174,
c       6,       3,       0,       0,       0,
     * 2460155, 3080188, 3660222, 4240256, 5070308, 5900360, 6600408, 
     * 7290456,
     * 7760492, 8000512, 8210526, 8220526, 8170520, 7970505, 7780491, 
     * 7680475,
     * 7480462, 7280450, 6900425, 6540396, 6340386, 6150376, 5840360, 
     * 5560340,
     * 4800300, 4300282, 3930270, 3730265, 3670265, 3700268, 3750273, 
     * 3880280,
     * 3900288, 3970288, 3640256, 3370237, 3100226, 2910218, 2750208, 
     * 2680202,
     * 2680202, 2680202, 2680020,
c       7,       3,       0,       0,       0,
     * 2800190, 3600207, 4300253, 5000300, 5920360, 6850420, 7480460, 
     * 8120500,
     * 8360520, 8610540, 8650545, 8700550, 8670546, 8650542, 8500531, 
     * 8350520,
     * 8170505, 8000490, 7550460, 7270443, 7000427, 6800415, 6400387, 
     * 6000360,
     * 5370339, 4930321, 4680314, 4410312, 4380313, 4360314, 4430319, 
     * 4490324,
     * 4600328, 4630330, 4320300, 3850275, 3500250, 3250240, 3120229, 
     * 3070225,
     * 3030222, 3030222, 3030222,
c       8,       3,       0,       0,       0,
     * 4500275, 5450315, 6250364, 7050413, 8070479, 9100545, 9650582,
     * 10200620,
     *10470640,10750660,10810665,10870670,10830666,10800662,10610646,
     * 10420630,
     *10190611, 9870593, 9430570, 9040544, 8850520, 8650518, 8280491, 
     * 7900465,
     * 7000420, 6500410, 6100395, 5850390, 5800395, 5750400, 5850410, 
     * 5950418,
     * 6000420, 6100422, 5560372, 5240348, 4940330, 4580320, 4450310, 
     * 4290294,
     * 4270292, 4270292, 4270292,
c      11,       3,       0,       0,       0,
     * 5320300, 6370360, 7340427, 8320495, 9440580,10570665,11320712,
     * 12070750,
     *12180757,12300765,12200757,12100750,11920740,11740730,11530715,
     * 11330700,
     *11140680,10950660,10380615,10010591, 9700570, 9610567, 9270544, 
     * 8900520,
     * 8070490, 7500470, 7100450, 6750448, 6700449, 6650450, 6700450, 
     * 6730452,
     * 6780456, 6820460, 6180408, 5740392, 5460376, 5200356, 5070347, 
     * 4950338,
     * 4880332, 4880332, 4880332,
c      13,       3,       0,       0,       0,
     * 8000470, 9800550,11100585,12400620,13500740,14600860,15150907,
     * 15700955,
     *15850967,16000980,15900970,15800960,15570940,15350920,15050890,
     * 14750860,
     *14500840,14250820,13750780,13270756,12950740,12720725,12360695,
     * 12000665,
     *10830637,10000615, 9480600, 9150590, 9050590, 8950590, 9000600, 
     * 9080608,
     * 9150610, 9220615, 8550550, 7950525, 7400510, 7050488, 6820470, 
     * 6600450,
     * 6600450, 6600450, 6600450,
c      20,       3,       0,       0,       0,
     *11000540,12100620,13400715,14700810,15770890,16850970,17321017,
     * 17801065,
     *17951070,18101075,17921067,17751060,17571047,17401035,17101007,
     * 16800980,
     *16400960,16000940,15400900,14900867,14600825,14400835,13900803,
     * 13400770,
     *12200720,11350690,10850675,10400670,10380668,10350665,10400670,
     * 10500680,
     *10600685,10650690, 9650610, 9400600, 8900585, 8500560, 8200550, 
     * 8050535,
     * 8000530, 8000530, 8000530,
c      22,       3,       0,       0,       0,
     *11750625,13630725,15160817,16700910,18101045,19501180,20001227,
     * 20501275,
     *20451262,20401250,20071225,19751200,19301175,18861150,18601125,
     * 18341100,
     *18031070,17731040,17200995,16660945,16350925,16020905,15380865,
     * 14740825,
     *13800810,12690780,12250760,11820745,11620740,11590740,11620740,
     * 11780750,
     *11900760,11970765,11020690,10350660, 9750635, 9450615, 9250600, 
     * 9050585,
     * 9050585, 9050585, 9050585/
c      26,       3,       0,       0,       0,
      DATA (ISIG3(J),J=431,817)/
     *14000725,16000840,17370930,18751020,19811110,20881200,21441247,
     * 22001295,
     *22101297,22201300,21971283,21751267,21501253,21251240,21001226,
     * 20751213,
     *20431194,20121175,19501125,18921059,18551042,18181022,17440986,
     * 16700950,
     *15300900,14300860,13700840,13150830,13150832,13150835,13300840,
     * 13450850,
     *13600860,13650865,12500785,11850735,11280705,10700680,10350650,
     * 10100630,
     *10100630,10100630,10100630,
c      29,       3,       0,       0,       0,
     *18900810,21150990,22271087,23401185,24401280,25001350,25501410,
     * 25651450,
     *25651475,25651500,25421492,25201485,24901462,24601440,24221417,
     * 23851395,
     *23471362,23101330,22101275,19981224,29851200,19721176,10461128,
     * 19201080,
     *17851035,16801000,16050970,15500960,15500960,15500960,15530960,
     * 15650970,
     *15750975,15900985,14800905,14250873,13550835,13050810,12700780,
     * 12500765,
     *12450763,12450763,12450763,
c      35,       3,       0,       0,       0,
     *24300925,26101125,27001250,27901375,28801500,29401600,29651680,
     * 29701750,
     *29701760,29701770,29451750,29201730,28801695,28401660,27801620,
     * 27201580,
     *26451540,25701500,25001450,23981346,23651330,23321314,22661282,
     * 22001250,
     *20501190,19261140,18251100,17681075,17581075,17491075,17501070,
     * 17781088,
     *17891095,18081110,16901035,16451005,15300940,14920917,14500880,
     * 14250860,
     *14250860,14250860,14250860,
c      42,       3,       0,       0,       0,
     *30001050,31811350,32501420,33001650,33801800,34101980,34702070,
     * 34502120,
     *34302105,34102090,33802070,33502050,33152015,32801980,32401950,
     * 32001920,
     *31601875,31201830,30501770,29521704,29001670,28461636,27381568,
     * 26301500,
     *25001435,23251350,21901300,21001230,20601220,20551235,20551235,
     * 20551235,
     *20671237,20851240,20001160,19001120,18351090,17701065,17201040,
     * 17001020,
     *16951015,16951015,16951015,
c      48,       3,       0,       0,       0,
     *30601025,31251275,31701440,32201625,32551740,32801800,32901880,
     * 32601920,
     *32651950,32701980,32351950,32001920,31601885,31201850,31001830,
     * 30801810,
     *30851765,30901720,29201650,28441585,28101560,27761538,27081494,
     * 26401450,
     *23621330,22301290,21151245,20501210,20351205,20201200,20251200,
     * 20401205,
     *20701205,21001230,19001130,17951085,17401060,16751000,16450985,
     * 16250975,
     *16200970,16200970,16200970,
c      50,       3,       0,       0,       0,
     *36501120,37501520,37651680,37901880,38002040,38002140,38002200,
     * 38002250,
     *37802260,37602270,37202235,36802200,36402180,36002160,35472115,
     * 34952070,
     *34352005,33751960,32801880,31711812,31201775,30731746,29791688,
     * 28851630,
     *27001540,25351470,24201425,23451385,23151375,23051380,23061385,
     * 23201390,
     *23351395,23501400,22351305,21501265,20801235,20101195,19651170,
     * 19401150,
     *19351150,19351150,19351150,
c      56,       3,       0,       0,       0,
     *52001450,50151850,50252100,49752350,49002550,48502700,47802825,
     * 47252900,
     *46622875,46002850,45452800,44902750,44222690,43552630,43052576,
     * 42552525,
     *41902462,41252400,40402300,38702226,38302200,37802174,36802122,
     * 35802070,
     *33301880,31101770,29551715,28601680,28521680,28451680,28851685,
     * 29001690,
     *29151700,29401720,28001635,26601560,25801530,24901460,24601440,
     * 24251410,
     *24201410,24201410,24201410,
c      74,       3,       0,       0,       0,
     *58901575,57002025,56102300,55802575,55502850,54803000,54003115,
     * 53003180,
     *52003130,51003080,50153005,49302930,48402865,47502800,46752735,
     * 46002670,
     *45002610,44002550,42802450,42212400,41702370,41192340,40172280,
     * 39152220,
     *36502110,34702000,32601920,31501880,31201850,30701800,30851805,
     * 31001810,
     *31201820,31601840,30701800,29301720,28201640,27501620,27101570,
     * 26551530,
     *26401530,26401530,26401530,
c      82,       3,       0,       0,       0,
     *70801740,68302220,66502500,65302820,64003080,62803300,61003420,
     * 58403500,
     *57503460,56603420,55903375,55203330,54253265,53303200,52453130,
     * 51603060,
     *50753000,49902940,48102850,46902758,46302710,45682662,44452566,
     * 43232470,
     *41302380,38702250,37002160,35502080,34902040,34652045,34672047,
     * 34752050,
     *34952055,35152060,34402010,33601980,31501830,30401780,29851735,
     * 29551710,
     *29401700,29401700,29401700/
c      92,       3,       0,       0,       0,
C      ISIG4- PI + MESON-NUCLEUS CROSS SECTIONS
C      FOUR THE LAST DECIMAL POSITION IN THE NUMBER-INELASTIC CROSS-SECTION.
C      FOUR THE FIRST ONE-TOTAL CROSS SECTION.
      DATA (ISIG4(J),J=1,430)/
     *  400018,  700038,  890050, 1080062, 1300080, 1520098, 1800117, 
     * 2080136,
     * 2420156, 2760176, 2980188, 3200200, 3290209, 3330212, 3300210, 
     * 3280208,
     * 3220204, 3100196, 2880176, 2600164, 2500157, 2400150, 2160134, 
     * 1960124,
     * 1440097, 1250090, 1120085, 1080082, 1090083, 1100086, 1170093, 
     * 1230097,
     * 1280100, 1350102, 1100083,  960077,  870075,  850074,  830072,  
     * 830072,
     *  830072,  830072,  830072,
c       2,       4,       0,       0,       0,
C     THE COMMENT c (HERE AND FUTHER) ENDS THE CROSS-SECTION ARRAY.
C     THE FIRST NUMBER (IN COMMENT c) IS TARGET CHARGE NUMBER,
C     THE SECOND NUMBER IS THE TYPE OF CROSS-SECTIONS (1-NEUTRONS,
C     2-PROTONS,3-MESONS -PI,4-MESONS +PI)
     *  960060, 1500095, 1860118, 2220142, 2710168, 3200194, 3750228, 
     * 4300262,
     * 4720290, 5140319, 5450345, 5650361, 5740364, 5740364, 5640354, 
     * 5520350,
     * 5350330, 5220319, 4900298, 4620280, 4470269, 4320258, 3980237, 
     * 3670216,
     * 3140200, 2760198, 2480183, 2320182, 2300180, 2330182, 2440187, 
     * 2570193,
     * 2700202, 2760207, 2500179, 2300172, 2150165, 2050159, 1940155, 
     * 1880154,
     * 1860154, 1860154, 1860154,
c       4,       4,       0,       0,       0,
     * 1400094, 1920132, 2430158, 2940184, 3610222, 4280260, 5110315, 
     * 5940370,
     * 6280389, 6620408, 6780420, 6850426, 6880428, 6840424, 6720416, 
     * 6560400,
     * 6430393, 6300386, 5980366, 5670340, 5500324, 5330308, 5040294, 
     * 4740280,
     * 4160257, 3690241, 3360236, 3190231, 3100232, 3130234, 3190238, 
     * 3330246,
     * 3420252, 3480256, 3100220, 2900210, 2680198, 2500187, 2450183, 
     * 2370176,
     * 2340174, 2340174, 2340174,
c       6,       4,       0,       0,       0,
     * 1500090, 2120140, 2700174, 3280208, 4140252, 5000300, 5900363, 
     * 6800426,
     * 7350467, 7620490, 7810504, 7820504, 7790500, 7630487, 7480474, 
     * 7400460,
     * 7230448, 7060437, 6720413, 6330381, 6160373, 6000365, 5690350, 
     * 5410330,
     * 4670292, 4190276, 3850267, 3680263, 3640264, 3680267, 3750273, 
     * 3880280,
     * 3900288, 3970288, 3640256, 3370237, 3100226, 2910218, 2750208, 
     * 2680202,
     * 2680202, 2680202, 2680202,
c       7,       4,       0,       0,       0,
     * 1700100, 2400145, 3150292, 3900240, 4800290, 5700340, 6550405, 
     * 7400470,
     * 7790494, 8180518, 8240524, 8300530, 8260526, 8220522, 8110513, 
     * 8000505,
     * 7820491, 7650477, 7250448, 7000430, 6750412, 6610404, 6230337, 
     * 5850350,
     * 5250330, 4830316, 4580310, 4440308, 4450309, 4470311, 4530317, 
     * 4490324,
     * 4600328, 4630330, 4320300, 3850275, 3500250, 3250240, 3120229, 
     * 3070225,
     * 3030222, 3030222, 3030222,
c       8,       4,       0,       0,       0,
     * 2100115, 3200210, 4250275, 5300340, 6620417, 7950495, 8770540, 
     * 9600585,
     * 9970607,10350630,10420637,10500645,10450641,10400637,10230621,
     * 10070605,
     * 9820588, 9570572, 9180550, 8820526, 8650505, 8460502, 8090478, 
     * 7730455,
     * 6850410, 6360401, 5980388, 5750383, 5700388, 5650393, 5780405, 
     * 5900414,
     * 5980418, 6100422, 5560372, 5240348, 4940330, 4580320, 4450310, 
     * 4290294,
     * 4270292, 4270292, 4270292,
c      11,       4,       0,       0,       0,
     * 2250120, 3500238, 4820314, 6150390, 7800500, 9450610,10330661,
     * 11220712,
     *11480723,11750735,11660727,11570720,11420711,11280703,11080679,
     * 10880655,
     *10660645,10450635, 9880590, 9580569, 9350550, 9290547, 9000526, 
     * 8700505,
     * 7870475, 7300455, 6900438, 6600440, 6560442, 6520445, 6600445, 
     * 6680450,
     * 6780456, 6820460, 6180408, 5740392, 5460376, 5200356, 5070347, 
     * 4950338,
     * 4880332, 4880332, 4880332,
c      13,       4,       0,       0,       0,
     * 2750160, 4450315, 5670407, 7900500, 9920622,11950745,13170807,
     * 14400870,
     *14620887,14850905,14800902,14750900,14550880,14350860,14100835,
     * 13850810,
     *13600790,13350770,12950740,12620724,12450710,12280696,11940668,
     * 11600640,
     *10500617, 9700595, 9230585, 8950575, 8860575, 8770575, 8870590, 
     * 8970600,
     * 9040602, 9130608, 8550550, 7950525, 7400510, 7050488, 6820470, 
     * 6600450,
     * 6600450, 6600450, 6600450,
c      20,       4,       0,       0,       0,
     * 3100193, 5100365, 7100487, 9100610,11150725,13200840,14670902,
     * 16150965,
     *16370982,16601000,16550997,16500995,16320985,16150975,15920955,
     * 15700935,
     *15400907,15100880,14600840,14200815,13800780,13700788,13200761,
     * 12800735,
     *11800690,11000665,10550660,10150655,10130653,10100650,10150655,
     * 10350670,
     *10450675,10550680, 9600610, 9400605, 8900585, 8500560, 8200550, 
     * 8050535,
     * 8000530, 8000530, 8000530,
c      22,       4,       0,       0,       0,
     * 3300210, 5750410, 7920558,10100707,12550858,15001010,16681067,
     * 18371125,
     *18561137,18751150,18451125,18201100,17851085,17511070,17211040,
     * 16911010,
     *16630985,16360960,15900920,14620812,14500775,14400776,14180778,
     * 13960780,
     *13050760,12190750,11900740,11480720,11380725,11340725,11440730,
     * 11630740,
     *11750750,11830755,10980690,10350660, 9750635, 9450615, 9250600, 
     * 9050585,
     * 9050585, 9050585, 9050585/
c      26,       4,       0,       0,       0,
      DATA (ISIG4(J),J=431,817)/
     * 3550230, 6050425, 8370602,11200780,13750904,16301025,17851090,
     * 19401155,
     *19701172,20101190,20101190,20101190,19951185,19801180,19521152,
     * 19251125,
     *19101112,18951100,18301050,17591020,17301000,17010980,16430940,
     * 15850900,
     *14900870,14000835,13400815,12900810,12900812,12900815,13100825,
     * 13300840,
     *13450850,13500855,12400780,11850735,11280705,10700680,10350650,
     * 10100630,
     *10100630,10100630,10100630,
c      29,       4,       0,       0,       0,
     * 3750245, 6650480, 9700680,12750880,17401030,19601140,21001200,
     * 21901275,
     *22301312,22701350,22801355,22901360,22701340,22501320,22201305,
     * 21901290,
     *21551260,21201230,20601215,19721140,19501120,19261100,18781060,
     * 18301020,
     *17100985,16200960,15600940,15100930,15110931,15120932,15150935,
     * 15400955,
     *15550960,15650968,14700900,14250873,13550835,13050810,12700780,
     * 12500765,
     *12450763,12450763,12450763,
c      35,       4,       0,       0,       0,
     * 4100270, 7300540,11300757,15300975,19201140,22001285,23851400,
     * 25201480,
     *25601517,26001555,26151567,26301580,26021552,25751525,25221552,
     * 24701470,
     *23951415,23201360,22851340,22091284,21851255,21571236,21051198,
     * 20531160,
     *19451120,18521085,17761060,17191045,17151045,17101045,17161045,
     * 17461065,
     *17591075,17781090,16751025,16451005,15300940,14920917,14500880,
     * 14250860,
     *14250860,14250860,14250860,
c      42,       4,       0,       0,       0,
     * 4650315, 8000590,12000880,17601220,21701460,24801580,27301700,
     * 28851770,
     *29271990,29701810,29751810,29801810,29751805,29701800,29301765,
     * 28901730,
     *28651705,28401680,27901630,26541556,26201530,25861504,25181452,
     * 24501400,
     *23351335,22051270,20801210,20201180,20101190,19901190,19901190,
     * 20151205,
     *20301210,20451210,19801150,18901115,18351090,17701065,17201040,
     * 17001020,
     *16951015,16951015,16951015,
c      48,       4,       0,       0,       0,
     * 4550310, 7800580,11700880,17001060,21201270,24001400,26001530,
     * 27201610,
     *27701635,28201660,28301670,28401680,28201660,28001640,27801620,
     * 27601600,
     *27401580,27201560,26401500,25871450,25601430,25381410,24941370,
     * 24501330,
     *22521280,21301230,20351200,19851180,19771175,19701170,19751175,
     * 20051180,
     *20351180,20701210,18801120,17851080,17401060,16751000,16450985,
     * 16250975,
     *16200970,16200970,16200970,
c      50,       4,       0,       0,       0,
     * 4720320, 8600640,13700960,19801280,24601550,27401680,29701800,
     * 31001870,
     *31501895,32001920,32101925,32201930,32101910,32001890,31621872,
     * 31251855,
     *30851822,30451790,29751725,28931678,28551650,28111624,27531572,
     * 26851520,
     *25251445,24001390,23001350,22501325,22351330,22201330,22351335,
     * 22651350,
     *22851360,23001365,22151295,21351260,20751230,20101195,19651170,
     * 19401150,
     *19351150,19351150,19351150,
c      56,       4,       0,       0,       0,
     * 9000680, 9000680,15000990,23501500,30201850,34202150,36502250,
     * 37752300,
     *38252325,38752350,38542340,38302330,38302355,37502280,37252255,
     * 37002230,
     *36352215,36302200,35502120,35502126,35502130,34982084,33941992,
     * 32901900,
     *30701780,28901670,28401635,27301600,27251603,27201605,27701610,
     * 28051615,
     *28281630,28651660,27701620,26401550,25701525,24901460,24601440,
     * 24251410,
     *24201410,24201410,24201410,
c      74,       4,       0,       0,       0,
     * 5150348, 9400707,15001040,24001650,32702100,37502400,40502580,
     * 41402640,
     *42002645,42602650,42302585,42002520,41402465,40802410,40352370,
     * 39902330,
     *39452290,39002250,38102190,37722156,37302130,36882104,36042052,
     * 35202000,
     *33701930,31861870,31101830,30101790,29901770,29851765,30051775,
     * 30201780,
     *30401790,30801800,30201775,29051710,27901620,22501620,27101570,
     * 26551530,
     *26401530,26401530,26401530,
c      82,       4,       0,       0,       0,
     * 4850334, 9600720,15801020,27001560,35502100,40502300,43202550,
     * 44202700,
     *45202790,46202880,46402880,46602880,46202820,45802760,45252710,
     * 44702660,
     *44702605,43502550,42952510,42452462,41872430,41372398,40372334,
     * 39382270,
     *37552130,35732060,34502000,33421970,33101950,32951950,33101960,
     * 33301960,
     *33751970,34051980,33501950,33381978,31351830,30401780,29851735,
     * 29551710,
     *29401700,29401700,29401700/
c----- 92,       4,       0,       0,       0,
      END
