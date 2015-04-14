
***********************************************************************
*                                                                     *
*        The following subroutines are written by N.M.Sobolevsky,     *
*        Institute for Nuclear Research RAS, Moscow,                  *
*        E-mail: sobolevs@AL20.inr.troitsk.ru                         *
*                                                                     *
***********************************************************************

      SUBROUTINE DECAY1
C Simulation on phase space of the next decays
C     PI-  -->  MU-,ANTI NU MU (100%) 
C     PI+  -->  MU+,NU MU (100%) 
C     PI0  -->  2 GAMMA (100%) 
C          K-  -->  MU-,ANTI NU MU (63.5%)
C              -->  PI-,PI0 (21.2%)
C              -->  PI-,PI-,PI+ (5.6%)
C              -->  PI-,PI0,PI0 (1.7%)
C              -->  MU-,ANTI NU MU,PI0 (3.2%)
C              -->  E-,ANTI NU E,PI0 (4.8%)
C          K+  -->  MU+,NU MU (63.5%)
C              -->  PI+,PI0 (21.2%)
C              -->  PI+,PI+,PI- (5.6%)
C              -->  PI+,PI0,PI0 (1.7%)
C              -->  MU+,NU MU,PI0 (3.2%)
C              -->  E+,NU E,PI0 (4.8%)
C              K0S  -->  PI+,PI- (69%)
C                   -->  PI0,PI0 (31%)
C              K0L  -->  PI0,PI0,PI0 (22%)
C                   -->  PI+,PI-,PI0 (12.4%)
C                   -->  PI+,MU-,ANTI NU MU (13.5%)
C                   -->  PI-,MU+,NU MU (13.5%)
C                   -->  PI+,E-,ANTI NU E (19.3%)
C                   -->  PI-,E+,NU E (19.3%)
C                   MU-  -->  E-,ANTI NU E,NU MU (100%)
C                   MU+  -->  E+,NU E,ANTI NU MU (100%)
C-----------------------------------------------------------
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      COMMON /INREAC/ COST,SINF,COSF,TINT,WINT,JPART,NUCLID,KSTATE
      COMMON /MANALL/ NTREE,LINE1,LINE2,NBRMAX,NPRMAX,NNEUTR,NGAMMA,
     *                MELSC,MINTER,MCUT,MOUT,MABSOR,MDECAY,NSTRAN,
     *                NEMIN,NEPLS,NNUEL,NANUE,NNUMU,NANUM
      COMMON /TREE/ TREE(15,400000),ITREE(7,400000),LTR13,LTR7,LTR300 ! GSITR
      COMMON /SECHAR/ SPT(6,5000),SNU(10,101),LS6,LS100,LS10,LS11  ! GSITR
C
C Cleaning of the first six lines of arrays SPT,SNU - corrected of 5.4.99
C Cleaning  of arrays SPT,SNU                       ! HION, 5.4.99
      DO 100 J=1,LS100                              ! HION, 5.4.99
      DO 100 I=1,LS6
  100 SPT(I,J)=0.
          DO 200 J=1,LS11                           ! HION, 5.4.99
          DO 200 I=1,LS10
  200     SNU(I,J)=0.
C
C Consignment of decaying particle from TREE to /INREAC/
      COST=TREE(4,LINE2+1)
      SINF=TREE(5,LINE2+1)
      COSF=TREE(6,LINE2+1)
          TINT=TREE(7,LINE2+1)
              JPART=ITREE(1,LINE2+1)
                  KSTATE=ITREE(7,LINE2+1)
C
C Check of particle type validity and KSTATE validity
      IF(JPART.EQ. 3 .OR. JPART.EQ. 4 .OR. JPART.EQ. 5 .OR.
     *   JPART.EQ. 8 .OR. JPART.EQ. 9 .OR. JPART.EQ.10 .OR. 
     *   JPART.EQ.11 .OR. JPART.EQ.15 .OR. JPART.EQ.16)THEN
          CONTINUE
      ELSE
          CONTINUE!WRITE(25,1000)JPART,IXINIT,NUMTRE
 1000     FORMAT(' DECAY1 DIAGNOSTIC: INVALID JPART=',I4,', IXINIT=',
     *           I12,', NUMTRE=',I6)
          STOP
      END IF
C
      IF(KSTATE.EQ.5 .OR. KSTATE.EQ.6)THEN
          CONTINUE
      ELSE
          CONTINUE!WRITE(25,1001)KSTATE,IXINIT,NUMTRE
 1001     FORMAT(' DECAY1 DIAGNOSTIC: INVALID KSTATE=',I2,', IXINIT=',
     *           I12,', NUMTRE=',I6)
          STOP
      END IF
C
C Check of neutral kaons types validity
      IF(JPART.EQ.10 .OR. JPART.EQ.11)THEN
          K0L=IFIX(TREE(10,LINE2+1))
          K0S=IFIX(TREE(11,LINE2+1))
      ELSE
          GO TO 300
      END IF
C
      IF((K0L.EQ.0 .AND. K0S.EQ.0) .OR. (K0L.NE.0 .AND. K0S.NE.0))THEN
          CONTINUE!WRITE(25,1002)K0L,K0S,IXINIT,NUMTRE
 1002     FORMAT(' DECAY1 DIAGNOSTIC: INVALID K0, K0L=',I2,', K0S=',I2,
     *           ', IXINIT=',I12,', NUMTRE=',I6)
          STOP
      ELSE
          CONTINUE
      END IF
C
  300 GO TO(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16),JPART
    1 STOP 'DECAY1-1'
    2 STOP 'DECAY1-2'
    3 CALL DECPIM
      RETURN
    4 CALL DECPIP
      RETURN
    5 CALL DECPI0
      RETURN
    6 STOP 'DECAY1-6'
    7 STOP 'DECAY1-7'
    8 CALL DECKMN
      RETURN
    9 CALL DECKPL
      RETURN
   10 IF(K0L.EQ.1)THEN
          CALL DECK0L
          RETURN
      ELSE
          CALL DECK0S
          RETURN
      END IF
   11 GO TO 10
   12 STOP 'DECAY1-12'
   13 STOP 'DECAY1-13'
   14 STOP 'DECAY1-14'
   15 CALL DECMUM
      RETURN
   16 CALL DECMUP
      RETURN
C
      END
 
 
 
      SUBROUTINE DECPIM
C Decay of PI- on MU- and ANTI NU MU
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      COMMON /INREAC/ COST,SINF,COSF,TINT,WINT,JPART,NUCLID,KSTATE
      COMMON /SECHAR/ SPT(6,5000),SNU(10,101),LS6,LS100,LS10,LS11  ! GSITR
C
C Mass and JPART of particles
      RM0=140.
          RM1 =105.
          VJP1=15.
              RM2 =0.
              VJP2=20.
C
C Did decay occures on fly (KSTATE=5) or at rest (KSTATE=6)?
      IF(KSTATE.EQ.5)THEN
          T0=TINT
      ELSE
          T0=0.
      END IF
C
C Two particle decay
      CALL DECAY2(RM0,T0,COST,SINF,COSF,
     *            RM1,RM2,
     *            T1,COST1,SINF1,COSF1,
     *            T2,COST2,SINF2,COSF2)
C
      SPT(1,1)=COST1
      SPT(2,1)=SINF1
      SPT(3,1)=COSF1
      SPT(4,1)=T1
      SPT(5,1)=VJP1
      SPT(6,1)=1.
          SPT(1,2)=COST2
          SPT(2,2)=SINF2
          SPT(3,2)=COSF2
          SPT(4,2)=T2
          SPT(5,2)=VJP2
          SPT(6,2)=1.
C
      RETURN
      END
 
 
 
      SUBROUTINE DECPIP
C Decay of PI+ on MU+ and NU MU
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      COMMON /INREAC/ COST,SINF,COSF,TINT,WINT,JPART,NUCLID,KSTATE
      COMMON /SECHAR/ SPT(6,5000),SNU(10,101),LS6,LS100,LS10,LS11  ! GSITR
C
C Mass and JPART of particles
      RM0=140.
          RM1 =105.
          VJP1=16.
              RM2 =0.
              VJP2=19.
C
C Did decay occures on fly (KSTATE=5) or at rest (KSTATE=6)?
      IF(KSTATE.EQ.5)THEN
          T0=TINT
      ELSE
          T0=0.
      END IF
C
C Two particle decay
      CALL DECAY2(RM0,T0,COST,SINF,COSF,
     *            RM1,RM2,
     *            T1,COST1,SINF1,COSF1,
     *            T2,COST2,SINF2,COSF2)
C
      SPT(1,1)=COST1
      SPT(2,1)=SINF1
      SPT(3,1)=COSF1
      SPT(4,1)=T1
      SPT(5,1)=VJP1
      SPT(6,1)=1.
          SPT(1,2)=COST2
          SPT(2,2)=SINF2
          SPT(3,2)=COSF2
          SPT(4,2)=T2
          SPT(5,2)=VJP2
          SPT(6,2)=1.
C
      RETURN
      END
 
 
 
      SUBROUTINE DECPI0
C Decay of PI0 on two gamma-quants
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      COMMON /INREAC/ COST,SINF,COSF,TINT,WINT,JPART,NUCLID,KSTATE
      COMMON /SECHAR/ SPT(6,5000),SNU(10,101),LS6,LS100,LS10,LS11  ! GSITR
C
C Mass and JPART of particles
      RM0=140.
          RM1 =0.
          VJP1=12.
              RM2 =0.
              VJP2=12.
C
C Did decay occures on fly (KSTATE=5) or at rest (KSTATE=6)?
      IF(KSTATE.EQ.5)THEN
          T0=TINT
      ELSE
          T0=0.
      END IF
C
C Two particle decay
      CALL DECAY2(RM0,T0,COST,SINF,COSF,
     *            RM1,RM2,
     *            T1,COST1,SINF1,COSF1,
     *            T2,COST2,SINF2,COSF2)
C
      SPT(1,1)=COST1
      SPT(2,1)=SINF1
      SPT(3,1)=COSF1
      SPT(4,1)=T1
      SPT(5,1)=VJP1
      SPT(6,1)=1.
          SPT(1,2)=COST2
          SPT(2,2)=SINF2
          SPT(3,2)=COSF2
          SPT(4,2)=T2
          SPT(5,2)=VJP2
          SPT(6,2)=1.
C
      RETURN
      END
 
 
 
      SUBROUTINE DECKMN
C Decay of K- according to branching ratios
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      COMMON /INREAC/ COST,SINF,COSF,TINT,WINT,JPART,NUCLID,KSTATE
      COMMON /SECHAR/ SPT(6,5000),SNU(10,101),LS6,LS100,LS10,LS11  ! GSITR
          DIMENSION BRANCH(6)
          DATA BRANCH /0.635, 0.847, 0.903, 0.920, 0.952, 1./
C
C Mass  of decaying particle
      RM0=495.
C
C Did decay occures on fly (KSTATE=5) or at rest (KSTATE=6)?
      IF(KSTATE.EQ.5)THEN
          T0=TINT
      ELSE
          T0=0.
      END IF
C
C Sampling of decay mode
      GAMMA=RRAN(IX)
      DO 10 J=1,6
          IF(GAMMA.LE.BRANCH(J))THEN
              MODE=J
              GO TO 20
          ELSE
              GO TO 10
          END IF 
   10 CONTINUE
   20 GO TO(1,2,3,4,5,6),MODE
C Decay of K- on MU- and ANTI NU MU (63.5%)
    1 RM1 =105.
      VJP1=15.
          RM2 =0.
          VJP2=20.
      GO TO 200
C Decay of K- on PI- and PI0 (21.2%)
    2 RM1 =140.
      VJP1=3.
          RM2 =140.
          VJP2=5.
      GO TO 200
C Decay of K- on PI-,PI- and PI+ (5.6%)
    3 RM1 =140.
      VJP1=3.
          RM2 =140.
          VJP2=3.
              RM3 =140.
              VJP3=4.
      GO TO 300
C Decay of K- on PI-,PI0 and PI0 (1.7%)
    4 RM1 =140.
      VJP1=3.
          RM2 =140.
          VJP2=5.
              RM3 =140.
              VJP3=5.
      GO TO 300
C Decay of K- on MU-,ANTI NU MU and PI0 (3.2%)
    5 RM1 =105.
      VJP1=15.
          RM2 =0.
          VJP2=20.
              RM3 =140.
              VJP3=5.
      GO TO 300
C Decay of K- on E-,ANTI NU E and PI0 (4.8%)
    6 RM1 =0.511
      VJP1=13.
          RM2 =0.
          VJP2=18.
              RM3 =140.
              VJP3=5.
      GO TO 300
C
C Two particle decay
  200 CALL DECAY2(RM0,T0,COST,SINF,COSF,
     *            RM1,RM2,
     *            T1,COST1,SINF1,COSF1,
     *            T2,COST2,SINF2,COSF2)
C
      SPT(1,1)=COST1
      SPT(2,1)=SINF1
      SPT(3,1)=COSF1
      SPT(4,1)=T1
      SPT(5,1)=VJP1
      SPT(6,1)=1.
          SPT(1,2)=COST2
          SPT(2,2)=SINF2
          SPT(3,2)=COSF2
          SPT(4,2)=T2
          SPT(5,2)=VJP2
          SPT(6,2)=1.
      RETURN
C
C Three particle decay
  300 CALL DECAY3(RM0,T0,COST,SINF,COSF,
     *            RM1,RM2,RM3,
     *            T1,COST1,SINF1,COSF1,
     *            T2,COST2,SINF2,COSF2,
     *            T3,COST3,SINF3,COSF3)
C
      SPT(1,1)=COST1
      SPT(2,1)=SINF1
      SPT(3,1)=COSF1
      SPT(4,1)=T1
      SPT(5,1)=VJP1
      SPT(6,1)=1.
          SPT(1,2)=COST2
          SPT(2,2)=SINF2
          SPT(3,2)=COSF2
          SPT(4,2)=T2
          SPT(5,2)=VJP2
          SPT(6,2)=1.
              SPT(1,3)=COST3
              SPT(2,3)=SINF3
              SPT(3,3)=COSF3
              SPT(4,3)=T3
              SPT(5,3)=VJP3
              SPT(6,3)=1.
      RETURN
C
      END
 
 
 
      SUBROUTINE DECKPL
C Decay of K+ according to branching ratios
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      COMMON /INREAC/ COST,SINF,COSF,TINT,WINT,JPART,NUCLID,KSTATE
      COMMON /SECHAR/ SPT(6,5000),SNU(10,101),LS6,LS100,LS10,LS11  ! GSITR
          DIMENSION BRANCH(6)
          DATA BRANCH /0.635, 0.847, 0.903, 0.920, 0.952, 1./
C
C Mass  of decaying particle
      RM0=495.
C
C Did decay occures on fly (KSTATE=5) or at rest (KSTATE=6)?
      IF(KSTATE.EQ.5)THEN
          T0=TINT
      ELSE
          T0=0.
      END IF
C
C Sampling of decay mode
      GAMMA=RRAN(IX)
      DO 10 J=1,6
          IF(GAMMA.LE.BRANCH(J))THEN
              MODE=J
              GO TO 20
          ELSE
              GO TO 10
          END IF 
   10 CONTINUE
   20 GO TO(1,2,3,4,5,6),MODE
C Decay of K+ on MU+ and NU MU (63.5%)
    1 RM1 =105.
      VJP1=16.
          RM2 =0.
          VJP2=19.
      GO TO 200
C Decay of K+ on PI+ and PI0 (21.2%)
    2 RM1 =140.
      VJP1=4.
          RM2 =140.
          VJP2=5.
      GO TO 200
C Decay of K+ on PI+,PI+ and PI- (5.6%)
    3 RM1 =140.
      VJP1=4.
          RM2 =140.
          VJP2=4.
              RM3 =140.
              VJP3=3.
      GO TO 300
C Decay of K+ on PI+,PI0 and PI0 (1.7%)
    4 RM1 =140.
      VJP1=4.
          RM2 =140.
          VJP2=5.
              RM3 =140.
              VJP3=5.
      GO TO 300
C Decay of K+ on MU+, NU MU and PI0 (3.2%)
    5 RM1 =105.
      VJP1=16.
          RM2 =0.
          VJP2=19.
              RM3 =140.
              VJP3=5.
      GO TO 300
C Decay of K+ on E+, NU E and PI0 (4.8%)
    6 RM1 =0.511
      VJP1=14.
          RM2 =0.
          VJP2=17.
              RM3 =140.
              VJP3=5.
      GO TO 300
C
C Two particle decay
  200 CALL DECAY2(RM0,T0,COST,SINF,COSF,
     *            RM1,RM2,
     *            T1,COST1,SINF1,COSF1,
     *            T2,COST2,SINF2,COSF2)
C
      SPT(1,1)=COST1
      SPT(2,1)=SINF1
      SPT(3,1)=COSF1
      SPT(4,1)=T1
      SPT(5,1)=VJP1
      SPT(6,1)=1.
          SPT(1,2)=COST2
          SPT(2,2)=SINF2
          SPT(3,2)=COSF2
          SPT(4,2)=T2
          SPT(5,2)=VJP2
          SPT(6,2)=1.
      RETURN
C
C Three particle decay
  300 CALL DECAY3(RM0,T0,COST,SINF,COSF,
     *            RM1,RM2,RM3,
     *            T1,COST1,SINF1,COSF1,
     *            T2,COST2,SINF2,COSF2,
     *            T3,COST3,SINF3,COSF3)
C
      SPT(1,1)=COST1
      SPT(2,1)=SINF1
      SPT(3,1)=COSF1
      SPT(4,1)=T1
      SPT(5,1)=VJP1
      SPT(6,1)=1.
          SPT(1,2)=COST2
          SPT(2,2)=SINF2
          SPT(3,2)=COSF2
          SPT(4,2)=T2
          SPT(5,2)=VJP2
          SPT(6,2)=1.
              SPT(1,3)=COST3
              SPT(2,3)=SINF3
              SPT(3,3)=COSF3
              SPT(4,3)=T3
              SPT(5,3)=VJP3
              SPT(6,3)=1.
      RETURN
C
      END
 
 
 
      SUBROUTINE DECK0S
C Decay of K0S on PI+,PI-(69%)  or on PI0,PI0(31%) 
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      COMMON /INREAC/ COST,SINF,COSF,TINT,WINT,JPART,NUCLID,KSTATE
      COMMON /SECHAR/ SPT(6,5000),SNU(10,101),LS6,LS100,LS10,LS11  ! GSITR
C
C Mass  of particles
      RM0=495.
          RM1 =140.
              RM2 =140.
C
C Sampling of decay mode (JPART of particles-products)
      IF(RRAN(IX).LE.0.69)THEN
          VJP1=4.
          VJP2=3.
      ELSE
          VJP1=5.
          VJP2=5.
      END IF
C
C Did decay occures on fly (KSTATE=5) or at rest (KSTATE=6)?
      IF(KSTATE.EQ.5)THEN
          T0=TINT
      ELSE
          T0=0.
      END IF
C
C Two particle decay
      CALL DECAY2(RM0,T0,COST,SINF,COSF,
     *            RM1,RM2,
     *            T1,COST1,SINF1,COSF1,
     *            T2,COST2,SINF2,COSF2)
C
      SPT(1,1)=COST1
      SPT(2,1)=SINF1
      SPT(3,1)=COSF1
      SPT(4,1)=T1
      SPT(5,1)=VJP1
      SPT(6,1)=1.
          SPT(1,2)=COST2
          SPT(2,2)=SINF2
          SPT(3,2)=COSF2
          SPT(4,2)=T2
          SPT(5,2)=VJP2
          SPT(6,2)=1.
C
      RETURN
      END
 
 
 
      SUBROUTINE DECK0L
C Decay of K0L according to branching ratios
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      COMMON /INREAC/ COST,SINF,COSF,TINT,WINT,JPART,NUCLID,KSTATE
      COMMON /SECHAR/ SPT(6,5000),SNU(10,101),LS6,LS100,LS10,LS11  ! GSITR
          DIMENSION BRANCH(6)
          DATA BRANCH /0.220, 0.344, 0.479, 0.614, 0.807, 1./
C
C Mass  of decaying particle
      RM0=495.
C
C Did decay occures on fly (KSTATE=5) or at rest (KSTATE=6)?
      IF(KSTATE.EQ.5)THEN
          T0=TINT
      ELSE
          T0=0.
      END IF
C
C Sampling of decay mode
      GAMMA=RRAN(IX)
      DO 10 J=1,6
          IF(GAMMA.LE.BRANCH(J))THEN
              MODE=J
              GO TO 20
          ELSE
              GO TO 10
          END IF 
   10 CONTINUE
   20 GO TO(1,2,3,4,5,6),MODE
C Decay of K0L on PI0,PI0 and PI0 (22%)
    1 RM1 =140.
      VJP1=5.
          RM2 =140.
          VJP2=5.
              RM3 =140.
              VJP3=5.
      GO TO 300
C Decay of K0L on PI+,PI- and PI0 (12.4%)
    2 RM1 =140.
      VJP1=4.
          RM2 =140.
          VJP2=3.
              RM3 =140.
              VJP3=5.
      GO TO 300
C Decay of K0L on PI+,MU- and ANTI NU MU (13.5%)
    3 RM1 =140.
      VJP1=4.
          RM2 =105.
          VJP2=15.
              RM3 =0.
              VJP3=20.
      GO TO 300
C Decay of K0L on PI-,MU+ and NU MU (13.5%)
    4 RM1 =140.
      VJP1=3.
          RM2 =105.
          VJP2=16.
              RM3 =0.
              VJP3=19.
      GO TO 300
C Decay of K0L on PI+,E- and ANTI NU E (19.3%)
    5 RM1 =140.
      VJP1=4.
          RM2 =0.511
          VJP2=13.
              RM3 =0.
              VJP3=18.
      GO TO 300
C Decay of K0L on PI-,E+ and NU E (19.3%)
    6 RM1 =140.
      VJP1=3.
          RM2 =0.511
          VJP2=14.
              RM3 =0.
              VJP3=17.
      GO TO 300
C
C Three particle decay
  300 CALL DECAY3(RM0,T0,COST,SINF,COSF,
     *            RM1,RM2,RM3,
     *            T1,COST1,SINF1,COSF1,
     *            T2,COST2,SINF2,COSF2,
     *            T3,COST3,SINF3,COSF3)
C
      SPT(1,1)=COST1
      SPT(2,1)=SINF1
      SPT(3,1)=COSF1
      SPT(4,1)=T1
      SPT(5,1)=VJP1
      SPT(6,1)=1.
          SPT(1,2)=COST2
          SPT(2,2)=SINF2
          SPT(3,2)=COSF2
          SPT(4,2)=T2
          SPT(5,2)=VJP2
          SPT(6,2)=1.
              SPT(1,3)=COST3
              SPT(2,3)=SINF3
              SPT(3,3)=COSF3
              SPT(4,3)=T3
              SPT(5,3)=VJP3
              SPT(6,3)=1.
C
      RETURN
      END
 
 
 
      SUBROUTINE DECMUM
C Decay of MU- on E-, ANTI NU E and NU MU
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      COMMON /INREAC/ COST,SINF,COSF,TINT,WINT,JPART,NUCLID,KSTATE
      COMMON /SECHAR/ SPT(6,5000),SNU(10,101),LS6,LS100,LS10,LS11  ! GSITR
C
C Mass and JPART of particles
      RM0=105.
          RM1 =0.511
          VJP1=13.
              RM2 =0.
              VJP2=18.
                  RM3 =0.
                  VJP3=19.
C
C Did decay occures on fly (KSTATE=5) or at rest (KSTATE=6)?
      IF(KSTATE.EQ.5)THEN
          T0=TINT
      ELSE
          T0=0.
      END IF
C
C Three particle decay
  300 CALL DECAY3(RM0,T0,COST,SINF,COSF,
     *            RM1,RM2,RM3,
     *            T1,COST1,SINF1,COSF1,
     *            T2,COST2,SINF2,COSF2,
     *            T3,COST3,SINF3,COSF3)
C
      SPT(1,1)=COST1
      SPT(2,1)=SINF1
      SPT(3,1)=COSF1
      SPT(4,1)=T1
      SPT(5,1)=VJP1
      SPT(6,1)=1.
          SPT(1,2)=COST2
          SPT(2,2)=SINF2
          SPT(3,2)=COSF2
          SPT(4,2)=T2
          SPT(5,2)=VJP2
          SPT(6,2)=1.
              SPT(1,3)=COST3
              SPT(2,3)=SINF3
              SPT(3,3)=COSF3
              SPT(4,3)=T3
              SPT(5,3)=VJP3
              SPT(6,3)=1.
C
      RETURN
      END
 
 
 
      SUBROUTINE DECMUP
C Decay of MU+ on E+, NU E and ANTI NU MU
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      COMMON /INREAC/ COST,SINF,COSF,TINT,WINT,JPART,NUCLID,KSTATE
      COMMON /SECHAR/ SPT(6,5000),SNU(10,101),LS6,LS100,LS10,LS11  ! GSITR
C
C Mass and JPART of particles
      RM0=105.
          RM1 =0.511
          VJP1=14.
              RM2 =0.
              VJP2=17.
                  RM3 =0.
                  VJP3=20.
C
C Did decay occures on fly (KSTATE=5) or at rest (KSTATE=6)?
      IF(KSTATE.EQ.5)THEN
          T0=TINT
      ELSE
          T0=0.
      END IF
C
C Three particle decay
  300 CALL DECAY3(RM0,T0,COST,SINF,COSF,
     *            RM1,RM2,RM3,
     *            T1,COST1,SINF1,COSF1,
     *            T2,COST2,SINF2,COSF2,
     *            T3,COST3,SINF3,COSF3)
C
      SPT(1,1)=COST1
      SPT(2,1)=SINF1
      SPT(3,1)=COSF1
      SPT(4,1)=T1
      SPT(5,1)=VJP1
      SPT(6,1)=1.
          SPT(1,2)=COST2
          SPT(2,2)=SINF2
          SPT(3,2)=COSF2
          SPT(4,2)=T2
          SPT(5,2)=VJP2
          SPT(6,2)=1.
              SPT(1,3)=COST3
              SPT(2,3)=SINF3
              SPT(3,3)=COSF3
              SPT(4,3)=T3
              SPT(5,3)=VJP3
              SPT(6,3)=1.
C
      RETURN
      END



      SUBROUTINE DECAY2(RM0,T0,COST0,SINF0,COSF0,
     *            RM1,RM2,
     *            T1,COST1,SINF1,COSF1,
     *            T2,COST2,SINF2,COSF2)
C Decay of particle 0 on two particles 1 and 2
C Input: RM0,T0 - mass and kinetic energy (MEV) of particle 0 in lab-system
C        COST0,SINF0,COSF0 - direction of fly of particle 0 in lab-system
C        RM1,RM2 - mass (MEV) of particles-products
C Output:T1,COST1,SINF1,COSF1 - kinetic energy and direction of particle 1 
C                               in lab-system
C        T2,COST2,SINF2,COSF2 - kinetic energy and direction of particle 2 
C                               in lab-system
C------------------------------------------------------------------------
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
C
C Check of mass compatibility
      SM2=RM1+RM2
      IF(SM2.GE.RM0)STOP 'DECAY2'
C
C If decay occures at rest, then appropriate direction of particle 0:
      IF(T0.LE.0.)THEN
          COST0=1.
          SINT0=0.
          SINF0=0.70710678
          COSF0=0.70710678
      ELSE  
          CONTINUE
      END IF
C
C Energy E0 and momentum P0 of partcle 0 in lab-system
      E0=RM0+T0
      P0=SQRT(T0*(T0+2.*RM0))
          SINT0=SQRT(1.-COST0**2)
          P0X=P0*SINT0*COSF0 
          P0Y=P0*SINT0*SINF0 
          P0Z=P0*COST0
C
C Energies and momenta of particles 1 and 2 in CMS
      E1S=(RM0**2+RM1**2-RM2**2)/(2.*RM0)
      E2S=(RM0**2-RM1**2+RM2**2)/(2.*RM0)
      P1S=SQRT((RM0**2+RM1**2-RM2**2)**2-4.*(RM0*RM1)**2)/(2.*RM0)
      P2S=P1S
C
C Isotropical direction of particle 1 in CMS
          CT1S=2.*RRAN(IX)-1.
          ST1S=SQRT(1.-CT1S**2)  
    1 G1=2.*RRAN(IX)-1.
      G2=RRAN(IX)
      GG=G1*G1+G2*G2
          IF(GG.LE.1.)THEN
              SF1S=(2.*G1*G2)/GG
              CF1S=(G1*G1-G2*G2)/GG
          ELSE
              GO TO 1
          END IF
C               (End of isotropical direction)
C
C If decay occures at rest, then transfer to lab-system don't required
      IF(T0.LE.0.)THEN
          E1=E1S
          T1=E1-RM1
            COST1=CT1S
            SINF1=SF1S
            COSF1=CF1S
                  E2=E2S
                  T2=E2-RM2
                    COST2=-COST1
                    SINF2=-SINF1
                    COSF2=-COSF1
          RETURN
      ELSE
          CONTINUE
      END IF
C
C Momenta of particles 1 and 2 in CMS
      P1XS=P1S*ST1S*CF1S
      P1YS=P1S*ST1S*SF1S
      P1ZS=P1S*CT1S
          P2XS=-P1XS
          P2YS=-P1YS
          P2ZS=-P1ZS
C
C Energies and momenta of particles 1 and 2 in lab-system
      E1=(E0*E1S+(P0X*P1XS+P0Y*P1YS+P0Z*P1ZS))/RM0
      T1=E1-RM1
          TEMP1=(E1+E1S)/(E0+RM0)
          P1X=P1XS+P0X*TEMP1 
          P1Y=P1YS+P0Y*TEMP1 
          P1Z=P1ZS+P0Z*TEMP1 
      E2=(E0*E2S+(P0X*P2XS+P0Y*P2YS+P0Z*P2ZS))/RM0
      T2=E2-RM2
          TEMP2=(E2+E2S)/(E0+RM0)
          P2X=P2XS+P0X*TEMP2 
          P2Y=P2YS+P0Y*TEMP2 
          P2Z=P2ZS+P0Z*TEMP2 
C
C Directions of particles 1 and 2 in lab-system
      CALL DECARS(P1X,P1Y,P1Z,P1,COST1,SINF1,COSF1)
      CALL DECARS(P2X,P2Y,P2Z,P2,COST2,SINF2,COSF2)
C
      RETURN
      END



      SUBROUTINE DECAY3(RM0,T0,COST0,SINF0,COSF0,
     *            RM1,RM2,RM3,
     *            T1,COST1,SINF1,COSF1,
     *            T2,COST2,SINF2,COSF2,
     *            T3,COST3,SINF3,COSF3)
C Decay of particle 0 on three particles 1,2 and 3
C Input: RM0,T0 - mass and kinetic energy (MEV) of particle 0 in lab-system
C        COST0,SINF0,COSF0 - direction of fly of particle 0 in lab-system
C        RM1,RM2,RM3 - mass (MEV) of particles-products
C Output:T1,COST1,SINF1,COSF1 - kinetic energy and direction of particle 1 
C                               in lab-system
C        T2,COST2,SINF2,COSF2 - kinetic energy and direction of particle 2 
C                               in lab-system
C        T3,COST3,SINF3,COSF3 - kinetic energy and direction of particle 3 
C                               in lab-system
C------------------------------------------------------------------------
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
C
C Check of mass compatibility
      SM3=RM1+RM2+RM3
      IF(SM3.GE.RM0)STOP 'DECAY3'
C
C Sampling of invariant mass of system of particles 1+2
      RIM2=(RM1+RM2)+(RM0-SM3)*RRAN(IX)
C
C Sampling of characteristics of particle 3 in lab-system
      CALL DECAY2(RM0,T0,COST0,SINF0,COSF0,
     *            RIM2,RM3,
     *            TIM2,CTIM2,SFIM2,CFIM2,
     *            T3,COST3,SINF3,COSF3)
C
C Sampling of characteristics of particles 1 and 2 in lab-system
      CALL DECAY2(RIM2,TIM2,CTIM2,SFIM2,CFIM2,
     *            RM1,RM2,
     *            T1,COST1,SINF1,COSF1,
     *            T2,COST2,SINF2,COSF2)
C
      RETURN
      END
