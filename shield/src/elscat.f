
***********************************************************************
*                                                                     *
*        The following subroutines are written by N.M.Sobolevsky,     *
*        Institute for Nuclear Research RAS, Moscow,                  *
*        E-mail: sobolevs@AL20.inr.troitsk.ru                         *
*                                                                     *
***********************************************************************

      SUBROUTINE ELSCAT
C     ELASTIC SCATTERING OF PARTICLE ON NUCLEUS IN LAB-SYSTEM.
C INPUT:
C       COST,SINF,COSF - PARTICLE-PROJECTILE DIRECTION
C       TINT - PARTICLE-PROJECTILE ENERGY (MEV)/(MeV/A)              !HION5
C       WINT - PARTICLE-PROJECTILE WEIGHT
C       JPART - PARTICLE-PROJECTILE TYPE
C       NUCLID - NO. OF NUCLEUS-TARGET
C       KSTATE - INTERACTION TYPE: 1 - ELASTIC SCATTERING
C                                  2 - INELASTIC INTERACTION
C OUTPUT:
C        SCATTERED PARTICLE CHARACTERISTICS - IN FIRST LINE OF ARRAY SPT,
C        RECOIL NUCLEUS - IN SECOND LINE OF ARRAY SNU
C---------------------------------------------------------------------------
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      COMMON /SECHAR/ SPT(6,5000),SNU(10,101),LS6,LS100,LS10,LS11  ! GSITR
      COMMON /INFEL/ ZNUC(110),ATWEI(110),DENS(110),RIEV(110),SYMB(110)
      CHARACTER*2 SYMB
      COMMON /INREAC/ COST,SINF,COSF,TINT,WINT,JPART,NUCLID,KSTATE
          COMMON /HIPROJ/ APROJ,ZPROJ                                !HION5
C
      REAL*8 ATD,RM,T,P2,P,PR,TET,COSTET,COS2,SINTET,RM0,
     *       EP,S,TEMP,EPRIM,E0PRIM,PPRIM,P0PRIM,COSTT0
C
      TINTOT=TINT                                                    !HION5
C
C CHECK OF PARTICLE TYPE VALIDITY
      IF((JPART.GE. 1 .AND. JPART.LE. 4) .OR.
     *   (JPART.GE. 6 .AND. JPART.LE.11) .OR.  
     *    JPART.EQ.15 .OR.  JPART.EQ.16  .OR.                        !HION5
     *   (JPART.GE.21 .AND. JPART.LE.25))GO TO 100                   !HION5
      CONTINUE!WRITE(25,1000)JPART,IXINIT,NUMTRE
 1000 FORMAT(' ELSCAT DIAGNOSTIC: INVALID JPART=',I4,', IXINIT=',I12,
     *       ', NUMTRE=',I6)
      STOP
  100 CONTINUE
C
C CLEANING OF THE FIRST FOUR LINES OF ARRAYS  SPT,SNU
      DO 101 J=1,LS100     ! KarSjh: Clean array SPT completely for any case
      DO 101 I=1,LS6
  101 SPT(I,J)=0.
          DO 102 J=1,LS11  ! KarSjh: Clean array SNU completely for any case
          DO 102 I=1,LS10
  102     SNU(I,J)=0.
C
C NUCLEUS-TARGET ATOMIC WEIGHT AND NUMBER
C      AT=ATWEI(NUCLID)
      AT=AMIXIS(NUCLID)                    ! 18.12.99
      ZT=ZNUC(NUCLID)                          
          ATD=DBLE(AT)
C
C APPROPRIATION OF PARTICLE MASS RM (GEV) ACCORDING TO PARTICLE TYPE JPART
      GO TO(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,                  !HION5
     *      17,18,19,20,21,22,23,24,25),JPART                        !HION5
    1 RM=0.940
      GO TO 110
    2 RM=0.940
      GO TO 110
    3 RM=0.140
      GO TO 110
    4 RM=0.140
      GO TO 110
    5 STOP 'ELSCAT5'
    6 RM=0.940
      GO TO 110
    7 RM=0.940
      GO TO 110
    8 RM=0.495
      GO TO 110
    9 RM=0.495
      GO TO 110
   10 RM=0.495
      GO TO 110
   11 RM=0.495
      GO TO 110
   12 STOP 'ELSCAT12'
   13 STOP 'ELSCAT13'
   14 STOP 'ELSCAT14'
   15 RM=0.105
      GO TO 110
   16 RM=0.105
      GO TO 110
   17 STOP 'ELSCAT17'                                                !HION5
   18 STOP 'ELSCAT18'                                                !HION5
   19 STOP 'ELSCAT19'                                                !HION5
   20 STOP 'ELSCAT20'                                                !HION5
   21 RM=1.880                                                       !HION5
      TINTOT=TINT*2.                                                 !HION5
      GO TO 110                                                      !HION5
   22 RM=2.820                                                       !HION5
      TINTOT=TINT*3.                                                 !HION5
      GO TO 110                                                      !HION5
   23 RM=2.820                                                       !HION5
      TINTOT=TINT*3.                                                 !HION5
      GO TO 110                                                      !HION5
   24 RM=3.760                                                       !HION5
      TINTOT=TINT*4.                                                 !HION5
      GO TO 110                                                      !HION5
   25 RM=0.940*APROJ                                                 !HION5
      TINTOT=TINT*APROJ                                              !HION5
      GO TO 110                                                      !HION5
C
C SAMPLING OF THE POLAR ELASTIC SCATTERING ANGLE TET IN LAB-SYSTEM
  110 T=DBLE(TINTOT*0.001)                                           !HION5
          P2=T*(T+2.*RM)
              P=DSQRT(P2)
      PR=7.098*P*ATD**0.333333333333
      TET=(1.7320508075688773/PR)*DSQRT(-DLOG(DBLE(RRAN(IX))))
          IF(TET.GT.3.1415926535897932)TET=3.1415926535897932
      COSTET=DCOS(TET)
          COS2=COSTET**2
      SINTET=DSQRT(1.-COS2)
C
      RM0=ATD*0.940
          IF(RM.LE.RM0)GO TO 120
              IF(SINTET.LE.RM0/RM)GO TO 120
                  SINTET=RM0/RM
                      COS2=1.-SINTET**2
                  COSTET=DSQRT(COS2)
                  GO TO 120
C
C PARTICLE AND NUCLEUS ENERGIES AND MOMENTA (GEV) IN LAB-SYSTEM
C (SEE BARASHENCOV AND TONEEV, P.86, BALDIN ET AL (1968), P.72)
  120 EP=RM0+RM+T
          S =2.*(RM**2+RM*RM0+RM0*T)
              TEMP=EP**2-P2*COS2
      EPRIM=(S*EP+P*COSTET*DSQRT(S**2-4.*(RM**2)*TEMP))/(2.*TEMP)
C
      E0PRIM=EP-EPRIM
      PPRIM=DSQRT(EPRIM**2-RM**2)
      P0PRIM=DSQRT(E0PRIM**2-RM0**2)
C
      IF(P0PRIM.EQ.0.)THEN
          COSTT0=0.       
      ELSE
          COSTT0=(P-PPRIM*COSTET)/P0PRIM
      END IF
      IF(COSTT0.GT.1.)COSTT0=1.
C
C UNIFORMLY DISTRIBUTED AZIMUTAL ANGLE FI
  130 Q=2.*RRAN(IX)-1.
      R=RRAN(IX)
          QQ=Q*Q
          RR=R*R
              QR=QQ+RR
      IF(QR.GT.1.)GO TO 130
          SINFIN=(2.*Q*R)/QR
          COSFIN=(QQ-RR)/QR
              SINFIP=-SINFIN
              COSFIP=-COSFIN
C
C FIRST LINE OF ARRAY  SNU
      SNU(1,1)=AT
      SNU(2,1)=ZT
      SNU(4,1)=1.
C
C SECOND LINE OF ARRAY SNU
      SNU(1,2)=AT
      SNU(2,2)=ZT
      SNU(4,2)=1.
      CALL SUBROT(COST,SINF,COSF,SNGL(COSTT0),SINFIN,COSFIN,
     *            CTR,SFR,CFR,STR)
      SNU(5,2)=SNGL(P0PRIM)*STR*CFR*1000.
      SNU(6,2)=SNGL(P0PRIM)*STR*SFR*1000.
      SNU(7,2)=SNGL(P0PRIM)*CTR*1000.
C
C
      IF(JPART.eq.25)THEN                                            !HION5
C       SCATTERED NUCLEUS IN THIRD LINE OF ARRAY SNU                 !HION5
        SNU(1,3)=APROJ                                               !HION5
        SNU(2,3)=ZPROJ                                               !HION5
        SNU(4,3)=1.                                                  !HION5
        CALL SUBROT(COST,SINF,COSF,SNGL(COSTET),SINFIP,COSFIP,       !HION5
     *              CTR,SFR,CFR,STR)                                 !HION5
        PSNMOD=SNGL(DSQRT((EPRIM-RM)*(EPRIM+RM)))                    !HION5
        SNU(5,3)=PSNMOD*STR*CFR*1000.                                !HION5
        SNU(6,3)=PSNMOD*STR*SFR*1000.                                !HION5
        SNU(7,3)=PSNMOD*CTR*1000.                                    !HION5
c          CONTINUE!WRITE(25,920)JPART,((EPRIM-RM)*1000.)/APROJ            !HIONDEB
c  920     format('ELSCAT: ',I2,F10.2)                            !HIONDEB
      ELSE                                                           !HION5
C       SCATTERED PARTICLE IN THE FIRST LINE OF ARRAY SPT
        CALL SUBROT(COST,SINF,COSF,SNGL(COSTET),SINFIP,COSFIP,
     *              CTR,SFR,CFR,STR)
        SPT(1,1)=CTR
        SPT(2,1)=SFR
        SPT(3,1)=CFR
        SPT(4,1)=SNGL((EPRIM-RM)*1000.)
        SPT(5,1)=FLOAT(JPART)
        SPT(6,1)=1.
c        IF(JPART.eq.24)THEN                                      !HIONDEB
c          CONTINUE!WRITE(25,921)JPART,SPT(4,1)/4.                         !HIONDEB
c  921     format('ELSCAT: ',I2,F10.2)                            !HIONDEB
c        END IF                                                   !HIONDEB
      END IF                                                         !HION5
C
      RETURN
      END
