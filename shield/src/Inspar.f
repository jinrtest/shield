
***********************************************************************
*                                                                     *
*        The following subroutines are written by N.M.Sobolevsky,     *
*        Institute for Nuclear Research RAS, Moscow,                  *
*        E-mail: sobolevs@AL20.inr.troitsk.ru                         *
*                                                                     *
***********************************************************************

      SUBROUTINE INSPAR
C Specification of incident particle
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      COMMON /WRPR/ WP(15),IWP(7)                                    !HION1
      COMMON /WRPRIN/ WP2(15),IWP2(7)                              ! Tim
      COMMON /TREE/ TREE(15,400000),ITREE(7,400000),LTR13,LTR7,LTR300 ! GSITR
      COMMON /MANALL/ NTREE,LINE1,LINE2,NBRMAX,NPRMAX,NNEUTR,NGAMMA,
     *                MELSC,MINTER,MCUT,MOUT,MABSOR,MDECAY,NSTRAN,
     *                NEMIN,NEPLS,NNUEL,NANUE,NNUMU,NANUM
      COMMON /T0JPR0/ TMAX0,JPART0    ! to get from SHIELD input
      COMMON /SOURN/ SOURN(13,100000),ISOURN(7,100000),LSN13,LSN7,LSN100   ! GSITR
          COMMON /HIPRO0/ APRO0,ZPRO0    ! for the input only!       !HION3
C
C ----------------- COMMONs of GEMCA --------------------           ! GEMCA
      COMMON/GDATAP/ X,Y,Z,CX,CY,CZ,XN,YN,ZN,CXN,CYN,CZN,           ! GEMCA
     *               NZONO,MEDOLD,NBPO,NSO,                         ! GEMCA
     *               NZONC,MEDCUR,NBPC,NSCI,NSCO,                   ! GEMCA
     *               PINSFL,IOLEG                                   ! GEMCA
	REAL*8 X,Y,Z,CX,CY,CZ,XN,YN,ZN,CXN,CYN,CZN                  ! GEMCA
	INTEGER PINSFL                                              ! GEMCA
C
C Coordinates (X,Y,Z, cm) of entry point into target:
      WP(1)=WP2(1)  ! Tim, 2012
      WP(2)=WP2(2)  ! Tim, 2012
      WP(3)=WP2(3)  ! Tim, 2012
C
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*C Gaussian distribution of the entry point (X,Y).                        *
*      RTARG=0.6    ! Radius of the target (cm).                          *
*      SIGMA=0.2    ! Standard deviation of X and Y (cm).                 *
*   10 G1=RRAN(IX)   ! Rundom number [0,1]                                 *
*      G2=RRAN(IX)   ! Rundom number [0,1]                                 *
*      WP(1)=SIGMA*SQRT(-2.0*ALOG(G1))*COS(6.2831852*G2)  ! X cm          *
*      WP(2)=SIGMA*SQRT(-2.0*ALOG(G1))*SIN(6.2831852*G2)  ! Y cm          *
*          IF(WP(1)**2+WP(2)**2 .GE. RTARG**2)GO TO 10                    *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C Direction of fly (Cos(theta), Sin(phi), Cos(phi)):
      WP(4)=WP2(4) ! Tim, 2012
      WP(5)=WP2(5) ! Tim, 2012
      WP(6)=WP2(6) ! Tim, 2012
C
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*C Gaussian distribution of the polar angle THETA and                     *
*C uniform distribution of the azimutal angle PHI.                        *
*      SIGMRD=10.              ! Standard deviation of THETA (mrad).      *
*C                             ! (1 mrad = 0.05729578 degrees)            *
*      SIGRAD=SIGMRD*0.001     ! Standard deviation of THETA (rad).       *
*      G3=RRAN(IX)              ! Rundom number [0,1]                      *
*      G4=RRAN(IX)              ! Rundom number [0,1]                      *
*      THETA=SIGRAD*SQRT(-2.0*ALOG(G3))*COS(6.2831852*G4)                 *
*      WP(4)=COS(THETA)                                                   *
*          IF(WP(4) .LE. 0.0)THEN                                         *
*              CONTINUE!WRITE(25,20)                                               *
*   20         FORMAT(' INSPAR DIAGNOSTIC: WP(4)=COS(THETA) .LE. 0.')     *
*              STOP                                                       *
*          END IF                                                         *
*      G5=RRAN(IX)              ! Rundom number [0,1]                      *
*      PHI=6.2831852*G5                                                   *
*      WP(5)=SIN(PHI)                                                     *
*      WP(6)=COS(PHI)                                                     *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C Kinetic energy (MeV):
      WP(7)=TMAX0                  ! SHIELD input: (MeV) or (MeV/A)  !HION3
C Statistic weight:
      WP(8)=1.
C Type of particle:
      IWP(1)=JPART0                   ! SHIELD input
C Primary Nucleus-projectile:                                        !HION3
      IF(JPART0.eq.25)THEN                                           !HION3
        WP(14)=APRO0                                                 !HION3
        WP(15)=ZPRO0                                                 !HION3
      END IF                                                         !HION3
C
C-------------------- ZONE DEFINITION ---------                     ! GEMCA
	X=DBLE(WP(1))                                               ! GEMCA
	Y=DBLE(WP(2))                                               ! GEMCA
	Z=DBLE(WP(3))                                               ! GEMCA
C                                                                   ! GEMCA
	COST=WP(4)                                                  ! GEMCA
	SINT=SQRT(1.0-COST**2.)                                     ! GEMCA
	SINF=WP(5)                                                  ! GEMCA
	COSF=WP(6)                                                  ! GEMCA
C                                                                   ! GEMCA
	CX=DBLE(SINT*COSF)                                          ! GEMCA
	CY=DBLE(SINT*SINF)                                          ! GEMCA
	CZ=DBLE(COST)                                               ! GEMCA
C                                                                   ! GEMCA
	CALL GNEXTZ(NEXTZ)                                          ! GEMCA
        WP(9)=FLOAT(NEXTZ)                                          ! GEMCA
C----------------------------------------
C The following information is permanent:
      IWP(2)=NTREE
      IWP(3)=0
      IWP(4)=0
      IWP(5)=0
      IWP(6)=-1
      IWP(7)=-1
      IF(IWP(1).EQ.10 .OR. IWP(1).EQ.11)WP(10)=1.
C
      DO 1 J=1,LTR13
    1 TREE(J,LINE1)=WP(J)
      DO 2 J=1,LTR7
    2 ITREE(J,LINE1)=IWP(J)
C
C If incident particle is NEUTRON with energy En<14.5 MeV,        ! BNAB28
C send NEUTRON from line No.LINE1 of TREE to line No.2 of SOURN.
      IF(ITREE(1,LINE1).eq.1 .and. TREE(7,LINE1).lt.14.5)THEN     ! BNAB28
          DO K=1,LSN13                                               !HION1
	      SOURN(K,2)   =TREE(K,LINE1)
	      TREE(K,LINE1)=0.
	  END DO
C
          DO K=1,LSN7                                                !HION1
	      ISOURN(K,2)   =ITREE(K,LINE1)
	      ITREE(K,LINE1)=0
	  END DO
C
	  NNEUTR=NNEUTR+1
      END IF
C
      RETURN
      END
