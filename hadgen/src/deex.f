
*****************************************************************
*                                                               *
*        The subroutine DEEX is written by N.M.Sobolevsky,      *
*        Institute for Nuclear Research RAS, Moscow,            *
*        E-mail: sobolevs@AL20.inr.troitsk.ru                   *
*                                                               *
*        Other subroutines are written by A.S.Botvina,          *
*        Institute for Nuclear Research RAS, Moscow,            *
*        E-mail: botvina@AL20.inr.troitsk.ru                    *
*                                                               *
*****************************************************************

      SUBROUTINE DEEX(ENEXT,ATWGHT,CHARGE,PNX,PNY,PNZ,KSTART,NSTART,
     *                INOVER)
C Interface program for including of A.S.Botvina complex nuclear deexci-
C tation code into SHIELG hadron-nucleus generator.
C Input: ENEXT - excitation energy (MeV)
C        ATWGHT, CHARGE - A and Z
C        PNX,PNY,PNZ - Momentum of nucleus (MeV)
C Input and Output:
C        KSTART - No. of first free line in array SPT before deexcitation
C        NSTART - The same for SNU
C        INOVER=1 is arrays SPT or SNU were exceeded, else INOVER=0.
C Output: Secondary particles and residual nuclei in arrays SPT and SNU
C         correspondingly
C----------------------------------------------------------------------
      COMMON /BLOKC/ SPTU(10,500)
      COMMON /BLFIS/ FIS,EF,SNT(6,6)
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      COMMON /SECHAR/ SPT(6,5000),SNU(10,101),LS6,LS100,LS10,LS11  ! GSITR
      COMMON /BLMULF/ IMULF
C               COMMONS to set the parameter values 
        COMMON /BLCAN/ ICAN
        COMMON /BLIMA/ IMA
        COMMON /BLPARA/ EPSIL0,FKAPPA,FKACOL
        COMMON /BL0999/ RNCL
        COMMON /BL1100/ AMS,AMFS
        COMMON /LEVERS/ PARLEV(40)    ! Parameters, see Subr. INLEVR
C
      ICAN=INT(PARLEV(23))  ! 0:  6 evaporation channels + fission
C                           ! 1: 32 evaporation channels + fission
      IMA= INT(PARLEV(24))  ! Maximal number of fragments (dflt - 3)
      EPSIL0=  PARLEV(25)   ! 1/a for inner degr. of freedom (dflt - 16.)
      FKACOL=  PARLEV(26)   ! 2.<FKACOL<5. (dflt - 2.)
      RNCL=    PARLEV(27)   ! dflt - 0.
      AMS =    PARLEV(28)   ! dflt - 0.125
      AMFS=    PARLEV(29)   ! dflt - 0.125
C
C Cleaning of the arrays SPTU and SNT
      DO 1 I=1,10
      DO 1 J=1,500
    1 SPTU(I,J)=0.
        DO 2 I=1,6
        DO 2 J=1,6
    2     SNT(I,J)=0.
C
      NSTAR1=NSTART      ! Initial value of NSTART is fixed
C      
C Input for deexcitation code (with transformation MeV to Gev)
      U=ENEXT*0.001
      A=ATWGHT
      Z=CHARGE
        P1=PNX*0.001
        P2=PNY*0.001
        P3=PNZ*0.001
            KSTAR1=1
C
C Check of deexcitation input
      IF(A.LT.0.0 .OR. Z.LT.0.0 .OR. U.LT.0.0)THEN
        CONTINUE!WRITE(25,1000)A,Z,U,IXINIT,NUMTRE
 1000     FORMAT(' DEEX DIAGNOSTIC: INVALID INPUT FOR DEEX, A=',
     *    F5.2,' Z=',F5.2,' U=',F8.1/20X,'IXINIT=',I12,'  NUMTRE=',I6)
        IF(A.LT.0.)A=0.
        IF(Z.LT.0.)Z=0.
        IF(U.LT.0.)U=0.
      END IF
      IF(A.LT.0.1)GO TO 21     ! Full desintegration before DEEX
C
C Condition of Fermi break up
      IF(A .LE. PARLEV(30))THEN  ! Default - 16. BotvinaDKFZ <=16!!!
**          PRINT 12
   12     FORMAT(40X,'Fermi break up')      
        IMULF=0
      ELSE
C         Condition of multifragmentation turning on/off
        IF(U .LE. (PARLEV(31)*A))THEN   ! Default - 0.002
            IMULF=0
        ELSE
            IMULF=1
**            PRINT 10
   10         FORMAT(40X,'Multifragmentation')      
        END IF
      END IF
C
                CALL ZVEZD1(U,A,Z,P1,P2,P3,KSTAR1)
C
C Printout on screen if fission occured
      IF(FIS .GT. 0.0)THEN
**        PRINT 11
   11     FORMAT(50X,'Fission')
      END IF
C
C Sending of secondary deexcitation particles from SPTU to SPT and SNU:
C     to SPT -  particles: n, p, d, t, He3, and He4
C     to SNU -  other fragments and nucleui
C
      DO J=1,500
C---------Is there a particle in line No. J?
        IF(SPTU(9,J) .GT. 0.0)THEN
            IA=NINT(SPTU(9,J)/0.94)
            IZ=NINT(SPTU(8,J))
C
C-------------Is it a particle or fragment in line No. J?
            IF((IA.EQ.1 .AND. (IZ.EQ.0 .OR. IZ.EQ.1)) .OR.
     *           (IA.EQ.2 .AND. IZ.EQ.1) .OR.
     *           (IA.EQ.3 .AND. (IZ.EQ.1 .OR. IZ.EQ.2)) .OR.
     *           (IA.EQ.4 .AND. IZ.EQ.2)) THEN
C
C-----------------It is a particle: n, p, d, t, He3, or He4
              IF(KSTART .GT. LS100)THEN
                 CONTINUE!WRITE(25,1001)KSTART,IXINIT,NUMTRE 
 1001                FORMAT('  Deex: Massive SPT exceeded, KSTART=',
     *                      I4,'  IXINIT=',I12,'  NUMTRE=',I12)
                 GO TO 22
              END IF
C                     Send particle in the line No. KSTART of array SPT
              SPT(1,KSTART)=SPTU(4,J)
              SPT(2,KSTART)=SPTU(5,J)
              SPT(3,KSTART)=SPTU(6,J)
              SPT(4,KSTART)=SPTU(7,J)*1000.
              SPT(5,KSTART)=SPTU(8,J)
              SPT(6,KSTART)=SPTU(9,J)*1000.
                  KSTART=KSTART+1
            ELSE
C-----------------It is a fragment
              IF(NSTART .GT. LS11)THEN
                 CONTINUE!WRITE(25,1002)NSTART,IXINIT,NUMTRE 
 1002                FORMAT('  Deex: Massive SNU exceeded, NSTART=',
     *                      I4,'  IXINIT=',I12,'  NUMTRE=',I12)
                 GO TO 23
              END IF
C                     Send fragment in the line No. NSTART of array SNU
C ******************************* HION9 *********Revised of 6.4.99*******
C 1.Checking of (A,Z) of the fragment and conversion of abnormal lightest 
C   fragments such as (2,0), (2,2), (3,0) etc. into legal fragment He3.
C   He3 is taken as approximately average of H2,H3,He3 and He4.
C 2.Conversion of abnormal fragments with A>4, and  Z=A, Z=A-1, Z=0, Z=1
C   to fragments (A,A/2).
C
                  AFR0=ANINT(SPTU(9,J)/0.94)
                  ZFR0=ANINT(SPTU(8,J))
**      CONTINUE!WRITE(25,*)'DEEX Fragment: AFR0,ZFR0=',AFR0,ZFR0
                  AFR=AFR0   ! Let fragment (AFR0,ZFR0) is OK
                  ZFR=ZFR0   ! Let fragment (AFR0,ZFR0) is OK
C
                  IF(AFR0.le.4.0)THEN ! Any fragment A.le.4 is abnormal.
C                                     ! Legal fragment are H2,H3,He3 and He4 only.
C##################### Conversion for A.le.4 #####################
                    IF(INT(AFR0).eq.1)CONTINUE
                    IF(INT(AFR0).eq.2)THEN ! Convertion to d
                      AFR=2.
                      ZFR=1.
                    END IF
                    IF(INT(AFR0).eq.3)THEN
                      if(RRAN(IX).le.0.5)then  ! Convertion to t
                        AFR=3.
                        ZFR=1.
                      else                    ! Convertion to He3
                        AFR=3.
                        ZFR=2.
                      end if
                    END IF
                    IF(INT(AFR0).eq.4)THEN ! Convertion to alpha
                      AFR=4.
                      ZFR=2.
                    END IF
                    CONTINUE!WRITE(25,1003)INT(AFR0),INT(ZFR0),INT(AFR),INT(ZFR)
 1003 FORMAT(' DEEX DIAGNOSTIC: Fragment (',I1,',',I1,') is converted',
     *       ' into (',I1,',',I1,').')
C##################### end of Conversion for A.le.4 ##############
                  END IF
C
      IF(AFR0.gt.4.0 .AND. ZFR0.ge.(AFR0-1.1))THEN  ! Any fragment (A,A) and
                    AFR=AFR0                        ! (A,A-1) is abnormal.
                    ZFR=ANINT(0.5*AFR0)        ! It is converted to (A,A/2)
                    CONTINUE!WRITE(25,1004)INT(AFR0),INT(ZFR0),INT(AFR),INT(ZFR)
 1004 FORMAT(' DEEX DIAGNOSTIC: Fragment (',I2,',',I2,') is converted',
     *       ' into (',I2,',',I2,')')
                  END IF
C
      IF(AFR0.gt.4.0 .AND. ZFR0.lt.1.1)THEN  ! Any fragment (A,0) and
                    AFR=AFR0                 ! (A,1) is abnormal.
                    ZFR=ANINT(0.5*AFR0)      ! It is converted to (A,A/2)
                    CONTINUE!WRITE(25,1004)INT(AFR0),INT(ZFR0),INT(AFR),INT(ZFR)
                  END IF
C
**      CONTINUE!WRITE(25,*)'DEEX Fragment after checking: AFR,ZFR=',AFR,ZFR
C **************************** end of HION9 *******************************
                  SNU(1,NSTART)=AFR                                  !HION9
                  SNU(2,NSTART)=ZFR                                  !HION9
              SNU(3,NSTART)=SPTU(10,J)
                    IF(FIS .GT. 0.0)THEN
                        SNU(4,NSTART)=3.
                    ELSE            
                        SNU(4,NSTART)=2.
                    END IF
C nonrelativistic PFRAG:      
C                     PFRAG =1000.0*SQRT(2.0*SPTU(9,J)*SPTU(7,J))
C relativistic PFRAG, Sobolevsky, 17.12.96:
                  PFRAG =1000.0*SQRT((SPTU(7,J)+2.0*SPTU(9,J))*
     *                                    SPTU(7,J))
                  STFRAG=SQRT(1.0-SPTU(4,J)**2)
              SNU(5,NSTART)=PFRAG*STFRAG*SPTU(6,J)
              SNU(6,NSTART)=PFRAG*STFRAG*SPTU(5,J)
              SNU(7,NSTART)=PFRAG*SPTU(4,J)
                  NSTART=NSTART+1
            END IF
        ELSE
            GO TO 20
        END IF
      END DO
C
C Input A=0 - full desintegration of nucleus before DEEX
   21 CONTINUE
      SNU(4,NSTAR1)=6.
      RETURN
C
C Array SPT is overfilled during deexitation
   22 CONTINUE
      INOVER=1
      RETURN
C
C Array SNU is overfilled during deexitation
   23 CONTINUE
      INOVER=1
      RETURN
C
C Normal output from DEEX
   20 CONTINUE
      IF(SNU(4,NSTAR1).LE.0.0) SNU(4,NSTAR1)=6.  ! Full desintegration
      RETURN                                     ! of nucleus in DEEX
C
      END


C All following codes were received from Botvina, HMI, via e-mail, 
C 05.04.95. Sobolevsky.
c *************** the version left in GSI on 14.09.1993 *************

      FUNCTION COLHOT(L,RADNCL)
      COMMON /BLHT06/ZJ(68)/BLHT10/ZFJ(68)/BLHT09/AFJ(68)
     */BLHT05/AJ(68)/BLHT03/U,A,Z
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
       IF(L-1)1,1,2
1      COLHOT=0.
       RETURN
2      TEMP1=1.44/RADNCL
       COLHOT=TEMP1*((ZJ(L)*ZFJ(L))/(AJ(L)**.33333+AFJ(L)**.33333))
       IF(COLHOT) 3,3,4
3       COLHOT=0.
4       CONTINUE
       RETURN
       END
      
      
      
      FUNCTION GAMHOT(J,A,U,AM,RADNCL)
       COMMON /BLHT09/AFJ(68)/BLHT15/RJ(68)/BLHT14/GAN(68)
      COMMON /BLHT11/VJ(68) /BLHT05/AJ(68)/BLHT06/ZJ(68)
      COMMON /BLHT10/ZFJ(68) /BLHTEX/EXTF(68)
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      IF(RJ(J).LE.0.OR.U.LE.0.) GO TO 12
      PER=2.*SQRT(AM*A*U)
      AMFR=0.5*AM*(1.+(40.**0.333333)/(AJ(J)**0.333333))
      IF(J.GT.6) PER=PER-2.*SQRT(AMFR*AJ(J)*EXTF(J))
      RN=1.5
      CC=0.2
      IF(J.GT.2.AND.J.LE.6)  CC=0.1
      IF(J.GT.6) CC=    (AJ(J)/AFJ(J))**0.6666667
       IF(J-1)1,2,1
2      ALFA=.76+2.2/AFJ(1)**.33333
       BETA=(2.12/AFJ(1)**.66667-.05)/ALFA
      GO TO 3
1      ALFA=1.+CC
       BETA=0.
       GO TO 3
3      Q1=AM*AFJ(J)
       Q2=AM*AFJ(J)*RJ(J)
       Q3=(GAN(J)*AFJ(J)**.66667)*(ALFA/Q1**2)*(AFJ(J)/(AFJ(J)+AJ(J)))
      Q3=Q3*AJ(J)*(3.1416*RN**2)/(2.*41.5*3.1416**2)
       Q4=(2.*BETA*Q1-3.)/2.+Q2
       Q5=(2.*BETA*Q1-3.)*(SQRT(Q2)-.5)+2.*Q2
      IF(PER-160.) 20,20,21
 20   PEX1=Q4*EXP(-PER)
      GO TO 22
 21   PEX1=0
 22   PP2=PER-2.*SQRT(Q2)
      IF(PP2-160.) 23,23,24
 23   PEX2=Q5*EXP(-PP2)
      GO TO 25
 24   PEX2=0.
 25   GAMHOT=Q3*(PEX1+PEX2)
      IF(J.LE.2) RETURN
      TFORM=AJ(J)*AJ(J)*10.
      GFORM=0.21*940./TFORM
      IF(GAMHOT.GT.GFORM) GAMHOT=GFORM
      RETURN
 12   GAMHOT=0.
      RETURN
      END
       
       
       
       FUNCTION TKIHOT(L,AM)
       COMMON /BLHT09/AFJ(68)/BLHT15/RJ(68)/BLHT11/VJ(68)
       COMMON /BLHT05/AJ(68)
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
       RJL=RJ(L)
       RB=4.*AM*AFJ(L)*RJL
       PP1=SQRT(RB)
 5     B1=RRAN(ix)
       IF(B1.LE.0.OR.B1.GT.1) B1=RRAN(ix)
       IF(PP1-160.) 21,21,22
 21    PEX1=EXP(-PP1)
       GO TO 23
 22    PEX1=0.
 23    RK=1.+(1./PP1)*ALOG(B1+(1.-B1)*PEX1)
       IF(L-1) 1,2,1
2      BETA=(2.12/AFJ(1)**0.66667-0.05)/(0.76+2.2/AFJ(1)**0.33333)
       Q1=1.+BETA/RJL
       Q2=Q1*SQRT(Q1)
       FRK=(((3.*SQRT(3.))/2.)/Q2)*(Q1*RK-RK**3)
       GO TO 3
1     FRK=((3.*SQRT(3.))/2.)*(RK-RK**3)
       GO TO 3
 3     B2=RRAN(ix)
       IF(B2-FRK) 4,4,5
4      TKIHOT=  RJL*(1.-RK**2)+VJ(L)
       RETURN
      END



       SUBROUTINE EVAHOT(ENEXT,ATWGHT,CHARGE,PNX,PNY,PNZ,KHOT)
C   EVAPORATION OF HOT FRAGMENTS
       COMMON /BLHT05/AJ(68)
     */BLHT06/ZJ(68)/BLHT14/GAN(68)/BLHT11/VJ(68)/BLHT15/RJ(68)
      COMMON /BLHT03/U,A,Z/BLHT09/AFJ(68)/BLHT10/ZFJ(68)
      COMMON /BLANGL/ANGL(4) /BLHTEX/EXTF(68)
      COMMON /SMPHOT/ SMPA(100),SMPZ(100),SMPE(100),SMPP(100,3)
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
        COMMON /LEVERS/ PARLEV(40)    ! Parameters, see Subr. INLEVR
      DIMENSION GJ(68),BJ(68),GG(68),VN(3),PN(3),PNL(3),PP(3),PPL(3)
      DIMENSION GGJ(68)
      U=ENEXT*1000.
      A=ATWGHT
      Z=CHARGE
      REMN=940.*A
      PNL(1)=PNX*1000.
      PNL(2)=PNY*1000.
      PNL(3)=PNZ*1000.
      ENL=SQRT(PNL(1)**2+PNL(2)**2+PNL(3)**2+REMN**2)
      VN(1)=PNL(1)/ENL
      VN(2)=PNL(2)/ENL
      VN(3)=PNL(3)/ENL
      AM=0.125
      RADNCL=1.5
      DO 1 III=1,100
      SMPA(III)=0.
      SMPZ(III)=0.
      SMPE(III)=0.
      SMPP(III,1)=0.
      SMPP(III,2)=0.
      SMPP(III,3)=0.
 1    CONTINUE
      IN=68
      DO 20 K=1,100
      UA=U/A
      IF(UA.LE.1.7) GO TO 11
C      IF(A.LE.16.) GO TO 11
      IF(A.LE.PARLEV(30)) GO TO 11   ! BotvinaDKFZ
      DO 4 I=1,IN
      VJ(I)=0.
 4    RJ(I)=-1000.
      CALL DELAM(A,Z,DL1,DSHEL1,BAR1)
      DO 6 I=1,IN
      AFJ(I)=A-AJ(I)
      ZFJ(I)=Z-ZJ(I)
      IF(AFJ(I).LT.ZFJ(I)) GO TO 6
      IF(AFJ(I).LT.AJ(I).OR.ZFJ(I).LT.ZJ(I)) GO TO 6
      VJ(I)=COLHOT(I,RADNCL)
      CALL DELAM(AFJ(I),ZFJ(I),DL2,DSHEL2,BAR2)
      CALL DELAM(AJ(I),ZJ(I),DL3,DSHEL3,BAR3)
      BJ(I)=DL2+DL3-DL1
      RR=U-(BJ(I)+VJ(I))
      EXTF(I)=0.
      IF(RR.GT.0.AND.I.GT.6) EXTF(I)=RR*AJ(I)/A
      RJ(I)=RR-EXTF(I)
 6    CONTINUE
       DO 7 I=1,IN
      GJ(I)=GAMHOT(I,A,U,AM,RADNCL)
7     CONTINUE
      G=0.
      DO 10 I=1,IN
10    G=G+GJ(I)
C ***************
c     IF(G.LE.0.) GO TO 11
c     PRINT 33,U,A,Z
c33   FORMAT(2X,' U=',F7.2,' A=',F5.1,' Z=',F5.1)
c     PRINT 30
c30   FORMAT(10X,'absolute widths:')
c     PRINT 31,(GJ(I),I=1,IN)
c31   FORMAT(2X,' GJ=',10E10.3)
c     DO 32 I=1,IN
c     GGJ(I)=GJ(I)/G
c32   CONTINUE
c     PRINT 35
c35   FORMAT(10X,'relative widths:')
c     PRINT 34,(GGJ(I),I=1,IN)
c34   FORMAT(2X,'GGJ=',10F7.4)
C ***************
      IF(G) 11,11,12
 11   PNX=PNL(1)*0.001
      PNY=PNL(2)*0.001
      PNZ=PNL(3)*0.001
      ENEXT=U*0.001
      ATWGHT=A
      CHARGE=Z
      SMPA(K)=A
      SMPZ(K)=Z
      SMPE(K)=ENEXT
      SMPP(K,1)=PNX
      SMPP(K,2)=PNY
      SMPP(K,3)=PNZ
      KHOT=K
      RETURN
 12   CONTINUE
      DO 13 J=2,IN
13    GJ(J)=GJ(J-1)+GJ(J)
      BB=RRAN(ix)
      B=BB*G
      DO 14 J=1,IN
      IF(B-GJ(J)) 15,14,14
15    LM=J
      GO TO 16
 14   CONTINUE
 16   IF(LM.LT.1.OR.LM.GT.IN) PRINT 300,LM
 300  FORMAT(2X,' ERROR IN EPAHOT - LM=',I5)
      IF(LM.LT.1.OR.LM.GT.IN) GO TO 11
      EP1=TKIHOT(LM,AM)
      EP2=ZJ(LM)
      EP3=940.*AJ(LM)
c ********
c     PRINT 36,LM,AJ(LM),ZJ(LM),EP1,EXTF(LM)
c36   FORMAT(2X,'LM=',I3,' AJ,ZJ=',2F5.1,' EP1=',F8.3,' EXTF=',
c    *F8.3)
c ********
      U=U-BJ(LM)-EP1-EXTF(LM)
      A=AFJ(LM)
      Z=ZFJ(LM)
C  VPM - RELATIVE VELOCITY OF FRAGMENT
      VPM=SQRT((2.*EP1)/(EP3*AFJ(LM)/(AFJ(LM)+AJ(LM))))
      CALL ISANGL
C  IN CMS
      PP(1)=VPM*ANGL(4)*ANGL(3)*EP3/(1.+AJ(LM)/A)
      PP(2)=VPM*ANGL(4)*ANGL(2)*EP3/(1.+AJ(LM)/A)
      PP(3)=VPM*ANGL(1)*EP3/(1.+AJ(LM)/A)
      PN(1)=-PP(1)
      PN(2)=-PP(2)
      PN(3)=-PP(3)
      EP=SQRT(PP(1)**2+PP(2)**2+PP(3)**2+EP3**2)
      EN=SQRT(PN(1)**2+PN(2)**2+PN(3)**2+(940.*A)**2)
      CALL CLPV(PP,VN,PPL,EP)
      CALL CLPV(PN,VN,PNL,EN)
      SMPA(K)=AJ(LM)
      SMPZ(K)=ZJ(LM)
      SMPE(K)=EXTF(LM)*0.001
      SMPP(K,1)=PPL(1)*0.001
      SMPP(K,2)=PPL(2)*0.001
      SMPP(K,3)=PPL(3)*0.001
      ENL=SQRT(PNL(1)**2+PNL(2)**2+PNL(3)**2+(940.*A)**2)
      VN(1)=PNL(1)/ENL
      VN(2)=PNL(2)/ENL
      VN(3)=PNL(3)/ENL
20    CONTINUE
      PRINT 21,U,A,Z
21    FORMAT(35X,37HMASSIVS SMP EXCEEDED AFTER EVAHOT    /40X,2HU=,
     *F10.5,4H  A=,F5.1,4H  Z=,F4.1)
       RETURN
       END
      
      
      
      BLOCK DATA HOT
      COMMON /BLHT05/AJ(68) /BLHT06/ZJ(68) /BLHT14/GAN(68)
      DATA AJ/
     * 1, 1, 2, 3, 3, 4, 6, 6, 7, 7, 8, 9, 9,10,11,11,12,13,13,14,
     *15,15,16,17,18,18,19,20,20,21,22,22,23,24,24,25,26,26,27,28,
     *28,29,30,30,31,32,32,33,33,34,35,35,36,36,36,37,37,38,39,39,
     *40,40,40,41,41,41,42,42/
      DATA ZJ/
     * 0, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7,
     * 7, 8, 8, 8, 8, 9, 9, 9,10,10,10,11,11,11,12,12,12,13,13,13,
     *14,14,14,15,15,15,16,15,16,16,16,17,16,17,18,17,18,18,18,19,
     *18,19,20,18,19,20,19,20/
      DATA GAN/2.,2.,3.,2.,2.,1.,62*1./
      END


c         change ZVEZD1, GAMMA1, E11, EPSIL0 and TC
C  PROGRAM OF MULTIFRAGMENTION (WITH RETURN). SMM1 for Soviet VAX
C------------------------------------------------------------------
      SUBROUTINE ZVEZD1(EE,AA,ZZ,P1,P2,P3,KST)
C MAIN PROGRAM OF MULTIFRAGMENTATION. INPUT: RES.NUCL. - AA,ZZ ; ITS
C EXCITATION - EE (GEV); MOMENTA IN LAB.SYS. - P1,P2,P3 (GEV/C);
C KST - FREE NUMBER IN SPT-MASSIVE (KST=1 - IN THE BEGIN.);
C OUTPUT: ALL PART. AND FRAGM. IN MASSIVE SPT(9,500).
C   *** CHANGED SWITCH.MICR.TO MACR.: INSTEAD SR00 - WW3  ***
      REAL N,MU,K,LT,L,KC
      COMMON /BLOKZV/SRN1(500),SRN2(500),SR1,SR2 /ADDNUC/IZPN,IAPN
      COMMON /BLFINL/A0F,Z0F,E0F,TF,FL,SRF,WW3F
      COMMON /BLOKN/N(500) /BLOKFR/FRA(500),SPSR,MFRAGT /BLIMA/IMA
      COMMON /BLPARA/ EPSIL0,FKAPPA,FKACOL /BLMOM/AMOM,AM(500),IAMO
      COMMON /BL503/WCOMP,WW2,WW3,WW4,WW5,WEVAP,WW4G
      COMMON /BL502/M2,M3,M4,M5,MC2,MC3,MC4,MC5
      COMMON /BLERR/ERRMI,ERRMA,ERRFR /BLMULF/IMULF /BENTRO/SSR
      COMMON /SMPHOT/ SMPA(100),SMPZ(100),SMPE(100),SMPP(100,3)
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
        COMMON /LEVERS/ PARLEV(40)    ! Parameters, see Subr. INLEVR
      DIMENSION PA(500),PAZ(500),AMP(500),ZMP(500)
      DIMENSION PN0(3,500),PNC(3),PNL(3),PN(3),VN(3)
      IAMO=0
      AMOM=0.
      A0=AA
      Z0=ZZ
      IA0=INT(A0+0.5)
C      DO 166 I7=1,IA0
C      AI7=I7
C      IF(I7.LE.4) AM(I7)=0.
C      IF(I7.GT.4) AM(I7)=0.4*2.69*(AI7**1.66667)
C 166  CONTINUE
      SR00=2.60
      IMA4=0
      IF(IMA.EQ.4.AND.A0.LE.110.) IMA4=1
      IF(IMA4.EQ.1) SR00=3.30
      IF(A0.LT.0.9) RETURN
      IREP=0
      E00=(EE/A0)*1000.
      IF(AA-PARLEV(30)) 107,107,108     ! Default - 16.
 107  A=AA
      Z=ZZ
      ERRMI=0
      ERRMA=0
      ERRFR=0
      IF(Z.GT.10.OR.(A-Z).GT.11) GO TO 105
      U=EE
      PN(1)=P1
      PN(2)=P2
      PN(3)=P3
      FRA(IA0)=FRA(IA0)+1.
      SPSR=SPSR+1.
      GO TO 157
 108  CONTINUE
      IF(IMULF.LE.0) GO TO 163
      IF(E00)163,163,165
 163  ERRMI=0
      ERRMA=0
      ERRFR=0
      GO TO 105
 165  CONTINUE
      REMN=0.94*A0
      EWHOLE=SQRT(P1*P1+P2*P2+P3*P3+REMN*REMN)
      VN(1)=P1/EWHOLE
      VN(2)=P2/EWHOLE
      VN(3)=P3/EWHOLE
      K=FKAPPA
      KC=FKACOL
      IF(IMA.EQ.0) WW3=1.
      IF(IMA.EQ.0) GO TO 102
      IF(A0.EQ.A0F.AND.Z0.EQ.Z0F.AND.E00.EQ.E0F) GO TO 60
      ERRMI=0
      ERRMA=0
      ERRFR=0
      E0I=E00
C      IF(IAMO.EQ.1) E0I=E0I+41.5*(AMOM**2)/(2.*AM(IA0)*A0)
CCC      CALL FRATE3(A0,Z0,E0I,K,KC,EPSIL0,TK,SR)
      CALL FRAMIK(A0,Z0,E0I,K,KC,EPSIL0,TK,SR)
      SR1=SR
      DO 201 IJ=1,IA0
 201  SRN1(IJ)=N(IJ)
C *******
C      PRINT 170,SR,TK,SSR
C 170  FORMAT(2X,'AFTER FRAMIK    SR=',F9.4,' TSR=',F9.4,' MEV',
C     *' ENTROP. SSR=',F9.4)
C      PRINT 121,(N(I),I=1,IA0)
C      PRINT 169,WCOMP,WEVAP,WW2,WW3,WW4,WW4G
C 169  FORMAT(2X,' WCOMP=',F8.5,' WEVAP=',F8.5,' WW2=',F8.5,' WW3=',
C     *F8.5,' WW4=',F8.5,' WW4G=',F8.5)
C      PRINT 168,M2,M3,M4,MC4
C 168  FORMAT(5X,' M2=',I10,' M3=',I10,' M4=',I10,' MC4=',I10)
C *******
      GO TO 62
 60   TK=TF
      L=FL
      SR=SRF
      WW3=WW3F
 62   CONTINUE
      IF(IMA4.EQ.1) GO TO 204
      IF(WW3-0.55) 101,101,102
 204  IF(SR-SR00) 101,101,102
 101  CONTINUE
      CALL CHOSAZ(A0,Z0,TK,IP,PA,PAZ)
      GO TO 103
 102  CONTINUE
      MFRAGT=MFRAGT+1
      IF(A0.EQ.A0F.AND.Z0.EQ.Z0F.AND.E00.EQ.E0F) GO TO 61
      CALL FRAGTE(A0,Z0,E00,K,KC,EPSIL0,TK,MU,L,SR)
      SR2=SR
      DO 202 IJ=1,IA0
 202  SRN2(IJ)=N(IJ)
C *******
C      PRINT 171,SR,TK,SSR
C 171  FORMAT(2X,'AFTER FRAGTE    SR=',F9.4,' TK=',F8.3,
C     *' MEV    ENTROP. SSR=',F9.4)
C      PRINT 121,(N(I),I=1,IA0)
C *******
 61   CONTINUE
      IF(IMA4.EQ.1) GO TO 207
      WWW=4.*(WW3-0.55)
      IF(WWW.LT.0.) WWW=0.
      IF(WWW.GT.1.) WWW=1.
      DO 203 IJ=1,IA0
 203  N(IJ)=(1.-WWW)*SRN1(IJ)+WWW*SRN2(IJ)
      SR=(1.-WWW)*SR1+WWW*SR2
 207  CONTINUE
      CALL SELECT(A0,Z0,N,TK,L,IP,PA,PAZ)
 103  CONTINUE
      IPROB=1
      DO 500 I=1,IP
      IA=INT(PA(I)+0.5)
      IZ=INT(PAZ(I)+0.5)
      IF(IA.GT.1.AND.IZ.GE.IA) IPROB=0
 500  IF(IA.GT.1.AND.IZ.LE.0) IPROB=0
      IF(IPROB.EQ.1) GO TO 501
      IF(IMA4.EQ.1) GO TO 205
      IF(WW3-0.55) 101,101,61
 205  IF(SR-SR00) 101,101,61
 501  CONTINUE
      A0F=A0
      Z0F=Z0
      E0F=E00
      TF=TK
      FL=L
      SRF=SR
      WW3F=WW3
C *******
c     PRINT 120,SR,TK
c120  FORMAT(5X,'SR=',F7.3,' TK=',F7.3,' MEV',15X,'N(A)')
c     PRINT 121,(N(I),I=1,IA0)
c121  FORMAT(2X,10F9.4)
c     PRINT 122,IP
c122  FORMAT(5X,' IP=',I3,' FRAGMENTS - PA,PAZ:')
c     PRINT 123,(PA(I),PAZ(I),I=1,IP)
c123  FORMAT(2X,20F6.1)
C *******
      IF(IP-1) 105,105,106
 105  A=AA
      Z=ZZ
      U=EE
      PX=P1
      PY=P2
      PZ=P3
      FRA(IA0)=FRA(IA0)+1.
      SPSR=SPSR+1.
      GO TO 156
 106  CONTINUE
      T=TK
      CALL FINDT(E00,PA,PAZ,IP,T,ETP,ECOLB,ECP,EEX,EB)
      IF(ERRFR.LT.0.5) GO TO 111
C *******
c     PRINT 800,IREP,TK,T,EB,ECOLB
c     PRINT 801,(PA(I),I=1,IP)
c801  FORMAT(2X,29F4.0)
c800  FORMAT(2X,'IREP=',I2,' TK,T=',2F9.5,' EB,ECOLB=',2F10.4,'   PA:')
C *******
      IREP=IREP+1
      IF(IREP.GT.10) GO TO 160
      IF(IMA4.EQ.1) GO TO 206
      IF(WW3-0.55) 101,101,61
 206  IF(SR-SR00) 101,101,61
 111  CONTINUE
C................Sobol, 2.1.97: new insert #1 of Botvina........
c **** including implementation additional number of protons.
c     IZPN=INT(0.5*(79.-ZZ))
c     IAPN=INT(0.5*(197.-AA))
      IF(IZPN.LT.1.OR.IAPN.LT.1) GO TO 411
      IPOLD1=IP+1
      IPNEW=IP+IZPN
      DO 410 J=IPOLD1,IPNEW
      PA(J)=1.
      PAZ(J)=1.
 410  CONTINUE
      IP=IPNEW
      ZZ0=IZPN+ZZ
      AA0=IAPN+AA
      ECOEFC=(0.6*1.44/(1.17*(1.+FKACOL)**0.333333))
      ECOLN=ECOEFC*ZZ0*ZZ0/(AA0**0.333333)
      DO 412 J4=1,IP
      ECOLN=ECOLN-ECOEFC*PAZ(J4)*PAZ(J4)/(PA(J4)**0.333333)
 412  CONTINUE
      ECOLB=ECOLN
 411  CONTINUE
c  *******
C...............Sobol, 2.1.97: End of new insert #1 of Botvina.......
      I=0
      INEUTR=0
      DO 400 J=1,IP
      IF(PAZ(J).LE.0.) GO TO 400
      I=I+1
      AMP(I)=0.94*PA(J)
      ZMP(I)=PAZ(J)
 400  CONTINUE
      IP0=I
      DO 401 J=1,IP
      IF(PAZ(J).GT.0.) GO TO 401
      INEUTR=INEUTR+1
      AMP(IP0+INEUTR)=0.94*PA(J)
      ZMP(IP0+INEUTR)=0.
 401  CONTINUE
C *******
c     PRINT 124,T,IP,INEUTR,IP0,ECOLB,ECP,EEX,EB
c124  FORMAT(2X,'T=',F7.3,' MEV  IP=',I3,' INEUTR=',I3,' IP0=',I3,
c    *' ECOLB=',F9.4,' MEV  ECP=',F9.4,'  EEX=',F9.4,'  EB=',F9.4)
C *******
      IF(IP0.LE.1) IP0=IP
CCC     TP=0.001*(ECOLB+1.5*IP0*T)
CCC     CALL DISIMP(IP0,AMP,PN0,TP)
      ECOLB=0.001*ECOLB
      CALL CULIMP(IP0,A0,Z0,AMP,ZMP,PN0,T,ECOLB)
      IF(IP0.GE.IP) GO TO 104
      TN=0.001*1.5*INEUTR*T
      CALL DISNET(INEUTR,IP0,AMP,PN0,T,TN)
 104  CONTINUE
      SPSR=SPSR+IP
      DO 35 I=1,IP
      DO 36 J=1,3
 36   PNC(J)=PN0(J,I)
C *******
c     EPAT=(PNC(1)**2+PNC(2)**2+PNC(3)**2)/(2.*AMP(I))
c     PRINT 125,IP,AMP(I),(PNC(J),J=1,3),EPAT
c125  FORMAT(2X,'PARTICLE IP=',I3,' AMP=',F7.3,'  PNC(1-3) (GEV/C)=',
c    *3F7.3,' EPAT (GEV)=',F7.3)
C *******
      EN=SQRT(PNC(1)**2+PNC(2)**2+PNC(3)**2+AMP(I)**2)
      CALL CLPV(PNC,VN,PNL,EN)
C *******
c     EPAT=(PNL(1)**2+PNL(2)**2+PNL(3)**2)/(2.*AMP(I))
c     PRINT 126,IP,AMP(I),(PNL(J),J=1,3),EPAT
c126  FORMAT(2X,'B  [.C.K. IP=',I3,' AMP=',F7.3,'  PNL(1-3) (GEV/C)=',
c    *3F7.3,' EPAT (GEV)=',F7.3)
C *******
      A=AMP(I)/0.94
      Z=ZMP(I)
      IA=INT(A+0.5)
      IZ=INT(Z+0.5)
      IF(IA.GT.0.AND.IA.LE.500) FRA(IA)=FRA(IA)+1.
CCC     GO TO 35
      IF(IA-3) 51,51,52
 51   U=0.
      GO TO 53
 52   U=UEVA(A,T)
 53   CONTINUE
C.................Sobol, 2.1.97: new insert #2 of Botvina...........
      PX=PNL(1)
      PY=PNL(2)
      PZ=PNL(3)
      CALL EVAHOT(U,A,Z,PX,PY,PZ,KHT)
      DO 703 I3=1,KHT
      U=SMPE(I3)
      A=SMPA(I3)
      Z=SMPZ(I3)
      PX=SMPP(I3,1)
      PY=SMPP(I3,2)
      PZ=SMPP(I3,3)
      IA3=INT(A+0.5)
      IZ3=INT(Z+0.5)
c *******************
c     PRINT 706,A,Z,U,PX,PY,PZ
c706  FORMAT(2X,'AFTER EVAHOT: A,Z=',2F5.1,' U=',F8.4,
c    *' PX,PY,PZ=',3F8.4)
c *******************
      IF(IZ3.GT.10.OR.(IA3-IZ3).GT.11) GO TO 705
      IF(IA3-INT(PARLEV(30))) 704,704,705        ! Default - 16
 704  PN(1)=PX
      PN(2)=PY
      PN(3)=PZ
      CALL RAZVAL(U,A,Z,PN,KST)
      GO TO 703
 705  CONTINUE
      CALL EVANUC(U,A,Z,PX,PY,PZ,KST)
 703  CONTINUE
 35   CONTINUE
      GO TO 160
 157  CALL RAZVAL(U,A,Z,PN,KST)
      GO TO 160
 156  CALL EVAHOT(U,A,Z,PX,PY,PZ,KHT)
      DO 700 I2=1,KHT
      U=SMPE(I2)
      A=SMPA(I2)
      Z=SMPZ(I2)
      PX=SMPP(I2,1)
      PY=SMPP(I2,2)
      PZ=SMPP(I2,3)
      IA2=INT(A+0.5)
      IZ2=INT(Z+0.5)
c *******************
c     PRINT 706,A,Z,U,PX,PY,PZ
c *******************
      IF(IZ2.GT.10.OR.(IA2-IZ2).GT.11) GO TO 701
      IF(IA2-INT(PARLEV(30))) 702,702,701       ! Default - 16
 702  PN(1)=PX
      PN(2)=PY
      PN(3)=PZ
      CALL RAZVAL(U,A,Z,PN,KST)
      GO TO 700
 701  CALL EVANUC(U,A,Z,PX,PY,PZ,KST)
 700  CONTINUE
C................Sobol, 2.1.97. End of new insert #2 of Botvina.....
      GO TO 160
 160  IF(ERRMI.GT.0.5.OR.ERRMA.GT.0.5) PRINT 161,AA,ZZ,EE,ERRMI,ERRMA
 161  FORMAT(2X,'SUSPECTING EVENT: AA,ZZ,EE=',3F9.4,' ERRMI,ERRMA=',
     *2F9.1)
      IF(ERRFR.GT.0.5) PRINT 162,AA,ZZ,EE,ERRFR
 162  FORMAT(5X,'EVENT EXCLUDED: AA,ZZ,EE=',3F9.4,' ERRFR=',F9.1)
      RETURN
      END

      
      
      BLOCK DATA DEFALT
      COMMON /BLMULF/ IMULF 
      COMMON /BLIMA/ IMA 
      COMMON /BLPARA/ EPSIL0,FKAPPA,FKACOL
c ---- for DEEX01 only; these COMMONs are absent in DEEX -----
        COMMON /RESNU0/ EMIN,EMAX,YMIN,YMAX 
        COMMON /BPLACE/ RFI 
        COMMON /ADDNUC/ IZPN,IAPN
c      
      DATA IMULF /1/
      DATA IMA /3/
      DATA EPSIL0 /16./
      DATA FKAPPA /1./
      DATA FKACOL /2./
c ---- for DEEX01 only -----
        DATA IZPN /0/
        DATA IAPN /0/
        DATA RFI /0./
        DATA EMIN /0./
        DATA EMAX /0./
        DATA YMIN /0./
        DATA YMAX /0./
      END
      
      
      
      SUBROUTINE FRAMIK(AA,ZZ,E00,K,KC,EPSIL0,TSR,MSR)
C HAXO[[EH[E [O [HEP[[[ E00 (MEV/N) M[[[T[[PA[MEHTHO[O COCTO[H[[
C B O[[ACT[ MA[[X [HEP[[[ BO[[[[[EH[[
C   ***********    KBA[[-M[KPOKAHOH[KA    *********
      REAL*4  K,KC,LT,MSR
      COMMON /BL500/W,W2(250),W3(25000),W4(10000) /BLIMA/IMA
      COMMON /BL501/IAF2(250),IAF3(25000),IAF4(10000)
      COMMON /BL502/M2,M3,M4,M5,MC2,MC3,MC4,MC5
      COMMON /BL503/WCOMP,WW2,WW3,WW4,WW5,WEVAP,WW4G
      COMMON /BLDIS/ J(20),IEND,JZ(20),IENDZ
      COMMON /BLOK1/W0,TC,RADNCL,G0 /BLOK0/A0,Z0,V0,XZ0,IA0,E0,E0Q
      COMMON /BLOK3/E11(500),XZ(500),GA(500) /BLOKT/TCON
      COMMON /BLOK2/B0,BT,DBT,LT,CP,PKP013 /BLOKN/SRW(500)
      COMMON /BLERR/ERRMI,ERRMA,ERRFR /BENTRO/SSR
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      ERRMI=0.
      A0=AA
      Z0=ZZ
      IA0=INT(A0+0.5)
      E11(1)=0.
      DO 100 I=2,IA0
      A=I
C 100  E11(I)=EPSIL0*(1.+3./(A-1))
C      IF(I.LE.50) E11(I)=EPSIL0*(1.+0.12*(A**0.33333-5.**0.33333))
C 100  IF(I.GT.50) E11(I)=EPSIL0*(1.+0.12*(50.**0.33333-5.**0.33333))
 100  E11(I)=EPSIL0/(1.+0.0020*(A/25.)**2)
      TCON=SQRT(E00/0.125)
      T=0.
      CALL PARAM(T,K,KC,EPSIL0)
      XZ0=Z0/A0
      A013=A0**0.3333333
      V0=(4.*3.1416/3.)*A0*RADNCL**3
      PKP013=1./(1.+KC)**0.3333333
      E0Q=A0*(-W0+G0*(1.-2*XZ0)**2)+B0*A013*A013
     *+0.6*1.44*Z0*Z0/(RADNCL*A013)
      CP=0.6*1.44*(1.-PKP013)/RADNCL
      DO 1 I=1,IA0
      A=I
      GA(I)=1.
      XZ(I)=XZ0
      SRW(I)=0.
 1    CONTINUE
      XZ(2)=0.5
      XZ(3)=0.5
      XZ(4)=0.5
      GA(2)=3.
      GA(3)=4.
      GA(4)=1.
      GA(1)=4.
      W=0.
      MSR=0.
      TSR=0.
      SSR=0.
      WW2=0.
      WW3=0.
      WW4=0.
      WW5=0.
      WEVAP=0.
      WW4G=0.
      CALL SCOMPA(E00,K,KC,EPSIL0,SCOMP)
      WCOMP=EXP(SCOMP-SCOMP)
      W=W+WCOMP
      WEVAP=WEVAP+WCOMP
      MSR=MSR+WCOMP*1.
      TSR=TSR+TCON*WCOMP
      SSR=SSR+WCOMP*SCOMP
      SRW(IA0)=SRW(IA0)+WCOMP
      M2=0
      M3=0
      M4=0
      M5=0
      MMAX=3
      IF(IMA.EQ.2) MMAX=2
      IF(IMA.EQ.4.AND.IA0.LE.110) MMAX=4
      DO 10 M=2,MMAX
      MC=0
      J(M)=IA0
      M1=M-1
      DO 2 I=1,M1
 2    J(I)=0
 4    CALL DISA(M)
      IF(IEND.EQ.1) GO TO 10
      MC=MC+1
      CALL CALCWE(M,J,SCOMP,E00,K,KC,EPSIL0,PR,ECOLB)
C ***********
C      IF(MC.LE.20) PRINT 101,M,TCON,PR,ECOLB,(J(I),I=1,M)
C 101  FORMAT(2X,'M=',I3,' TCON=',F8.3,' PR=',F8.5,' ECOLB=',F8.3,' J=',
C     *5I4)
C ***********
      W=W+PR
      DO 6 I=1,M
      IA=J(I)
 6    SRW(IA)=SRW(IA)+PR
      MSR=MSR+M*PR
      TSR=TSR+TCON*PR
      IF(PR.GT.0.) SSR=SSR+PR*(SCOMP+ALOG(PR))
      IF(M.EQ.2) GO TO 22
      IF(M.EQ.3) GO TO 23
      IF(M.EQ.4) GO TO 24
      IF(M.EQ.5) GO TO 25
 22   MC2=MC
      WW2=WW2+PR
      IF(J(1).LE.4) WEVAP=WEVAP+PR
      M2=M2+1
      W2(M2)=PR
      IAF2(M2)=J(2)+1000*J(1)
      GO TO 4
 23   MC3=MC
      WW3=WW3+PR
      M3=M3+1
      W3(M3)=PR
      IAF3(M3)=J(3)+1000*J(2)+1000000*J(1)
      GO TO 4
 24   MC4=MC
      WW4=WW4+PR
      M4=M4+1
      IF(J(3).LE.4) WW4G=WW4G+PR
      IF(M4.GT.10000) GO TO 4
      W4(M4)=PR
      IAF4(M4)=J(4)+1000*J(3)+100000*J(2)+10000000*J(1)
      GO TO 4
 25   MC5=MC
      WW5=WW5+PR
      M5=M5+1
      GO TO 4
 10   CONTINUE
      IF(M4.GT.10000) M4=10000
      DO 12 I=1,IA0
 12   SRW(I)=SRW(I)/W
      DO 63 I8=1,M2
 63   W2(I8)=W2(I8)/W
      DO 64 I8=1,M3
 64   W3(I8)=W3(I8)/W
      DO 65 I8=1,M4
 65   W4(I8)=W4(I8)/W
      MSR=MSR/W
      TSR=TSR/W
      SSR=SSR/W
      WCOMP=WCOMP/W
      WW2=WW2/W
      WW3=WW3/W
      WW4=WW4/W
      WW4G=WW4G/W
      WW5=WW5/W
      WEVAP=WEVAP/W
      RETURN
      END
      
      
      
      SUBROUTINE CALCWE(M,J,SCOMP,E00,K,KC,EPSIL0,PR,ECOLB)
C  O[PE[E[EH[E BEPO[THOCT[ KOH[[[[PA[[[ PR=EXP(S(J(1-M),E0)
C       [ K[[OHOBCKO[O [AP[EPA ECOLB (MEV).
C  BXO[:M [PA[MEHTOB B MACC[BE J(1-M),[HEP[.BO[[.-E00 (MEV/N),
C  [APAMETP[ MO[E[[ K,KC,EPSIL0 [ [HTPO[[[ KOM[.[[PA SCOMP ([[[ [EPEC-
C  [ETA PR).
      REAL LT,K,KC
      COMMON /BLOK1/W0,TC,RADNCL,G0 /BLOK0/A0,Z0,V0,XZ0,IA0,E0,E0Q
      COMMON /BLOK3/E11(500),XZ(500),GA(500) /BLOKT/TCON
      COMMON /BLOK2/B0,BT,DBT,LT,CP,PKP013 /BLMOM/AMOM,AM(500),IAMO
      COMMON /BLERR/ERRMI,ERRMA,ERRFR
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      DIMENSION EA(500),ECOLA(500),SAIN(500),J(20)
      E0=E00*A0
      A013=A0**0.333333
      VF=V0*((1.+1.4*(M**0.333333-1)/(RADNCL*A013))**3-1.)
      M1=M-1
      FAKT=1.
      DO 32 I=1,M1
      FFF=1.
      IM=I+1
      DO 33 II=IM,M
      IF(J(I).EQ.J(II)) FFF=FFF+1.
 33   CONTINUE
      FAKT=FAKT*FFF
 32   CONTINUE
      PRFAKT=1.
      PRFAKT=PRFAKT*FAKT
      ECON=0.
      ECOLB=0.
      PRGA=1.
      PRA32=1.
      PAM32=1.
C      SAM=0.
      DO 34 I=1,M
      IA=J(I)
C      SAM=SAM+AM(IA)
C      IF(AM(IA).GT.0.) PAM32=PAM32*(AM(IA)*SQRT(AM(IA)))
      A=IA
      PRGA=PRGA*GA(IA)
      PRA32=PRA32*(A*SQRT(A))
      IF(IA-1) 41,41,42
 41   ECOLA(IA)=CP*XZ(IA)*XZ(IA)
      EA(IA)=ECOLA(IA)
      GO TO 50
 42   IF(IA-4) 43,43,44
 43   ECOLA(IA)=CP*XZ(IA)*XZ(IA)*(A**(5./3.))
      IF(IA.EQ.2) EA(IA)=-2.796+ECOLA(IA)
      IF(IA.EQ.3) EA(IA)=-9.224+ECOLA(IA)
      IF(IA.EQ.4) EA(IA)=-30.11+ECOLA(IA)
      GO TO 50
 44   ECOLA(IA)=CP*XZ(IA)*XZ(IA)*(A**(5./3.))
      EA(IA)=(-W0+G0*(1.-2.*XZ(IA))**2)*A+
     #B0*A**0.666667+ECOLA(IA)
 50   CONTINUE
      ECON=ECON+EA(IA)
      ECOLB=ECOLB+ECOLA(IA)
 34   CONTINUE
C      EMOM=41.5*(AMOM**2)/(2*SAM)
      ECOLB=ECOLB+0.6*1.44*Z0*Z0*PKP013/(RADNCL*A013)
      ECON=ECON+0.6*1.44*Z0*Z0*PKP013/(RADNCL*A013)
C      IF(IAMO.EQ.1) ECON=ECON+EMOM
      DO 51 I=1,M
      IA=J(I)
      A=IA
      ECOLB=ECOLB-0.6*1.44*XZ(IA)*XZ(IA)*(A**(5./3.))/RADNCL
 51   CONTINUE
      IF((E0+E0Q-ECON).LT.0.003) GO TO 56
      T=SQRT(E00/0.125)
      IF(T.LT.0.0012) T=0.0012
      HT=0.5
      K1=0
      K2=0
      ID=0
      ICNT=0
 24   CONTINUE
      IF(ICNT.GT.120) GO TO 71
      CALL PARAM(T,K,KC,EPSIL0)
      ECON=0.
      DO 61 I=1,M
      IA=J(I)
      A=IA
      IF(IA-1) 55,55,52
 55   EA(IA)=ECOLA(IA)
      GO TO 60
 52   IF(IA-4) 53,53,54
 53   IF(IA.EQ.2) EA(IA)=-2.796+ECOLA(IA)
      IF(IA.EQ.3) EA(IA)=-9.224+ECOLA(IA)
      IF(IA.EQ.4) EA(IA)=-30.11+ECOLA(IA)+4.*T*T/E11(4)
      GO TO 60
 54   EA(IA)=(-W0+T*T/E11(IA)+G0*(1.-2.*XZ(IA))**2)*A+
     #(BT-T*DBT)*A**0.666667+ECOLA(IA)
 60   CONTINUE
      ECON=ECON+EA(IA)
 61   CONTINUE
      ECON=ECON+0.6*1.44*Z0*Z0*PKP013/(RADNCL*A013)+1.5*T*(M-1)
C      IF(IAMO.EQ.1) ECON=ECON+EMOM+1.5*T*(M-1)
      D=(E0+E0Q-ECON)/E0
      IF(ABS(D).LT.0.003) GO TO 29
      ICNT=ICNT+1
      H=SIGN(HT,D)
      IF(D) 21,21,22
 21   K1=1
      IF(K1.EQ.1.AND.K2.EQ.1) ID=ID+1
      IF(ID.GT.29) GO TO 71
      T=T+H/2**ID
      IF(T.GE.0.001) GO TO 24
      K2=1
      H=HT
      GO TO 21
 22   K2=1
      IF(K1.EQ.1.AND.K2.EQ.1) ID=ID+1
      IF(ID.GT.29) GO TO 71
      T=T+H/2**ID
      IF(T.GE.0.001) GO TO 24
      K1=1
      H=HT
      GO TO 22
 71   ERRMI=ERRMI+1
C      PRINT 72,ICNT,ID,D,T
C 72   FORMAT(5X,'ERROR IN CALCWE  ICNT=',I3,' ID=',I3,' D=',E12.5,
C     *'  T=',F8.3)
 56   PR=0.
      RETURN
 29   TCON=T
      SCON=0.
      DO 62 I=1,M
      IA=J(I)
      A=IA
      IF(IA-4) 63,63,64
 63   SAIN(IA)=0.
      IF(IA.EQ.4) SAIN(IA)=2*T*4./E11(IA)
      GO TO 65
 64   SAIN(IA)=2*T*A/E11(IA)-DBT*A**0.666667
 65   CONTINUE
      SCON=SCON+SAIN(IA)
 62   CONTINUE
      LT=16.15/SQRT(T)
      STRAN=ALOG(PRA32/PRFAKT)+(M-1)*ALOG(VF/LT**3)+1.5*(M-1)-
     *ALOG(A0*SQRT(A0))
      IF(STRAN.LT.0.) STRAN=0.
      SCON=SCON+ALOG(PRGA)+STRAN
C      SAMOM=1.5*(M-1)-3*(M-1)*ALOG(LT)+ALOG(PAM32)-
C     *ALOG(SAM*SQRT(SAM))
C      IF(SAMOM.LT.0.) SAMOM=0.
C      IF(IAMO.EQ.1) SCON=SCON+SAMOM
      PR=EXP(SCON-SCOMP)
C     IF(M.LE.2) PRINT 100,M,(J(I),I=1,M),ECOLB,T,STRAN,SAMOM,SCON,PR
C100  FORMAT(2X,'M=',I3,' J=',2I3,' ECOLB=',F7.3,' T=',F7.3,' STRAN=',
C    *F7.3,' SAMOM=',F7.3,' SCON=',F7.3,' PR=',E10.3)
      RETURN
      END
      
      
      
      SUBROUTINE DISA(K)
C B[[OP O[HO[O PAC[PE[E[EH[[ A0 H[K[OHOB [O K [PA[MEHTAM C [[ETOM HEPA[-
C [[[[MOCT[. KOH[[[[PA[[[: J(1)...J(K)-CO[EP[[T MACC[ [PA[MEHTOB.
C HA[. [C[OB[E: J(K)=A0,J(K-1)=...J(1)=O -[A[AETC[ [EPE[ DISA.
C [P[[HAK KOH[A IEND=1 O[HA[AET,[TO BCE KOH[[[[PA[[[ [[E [C[EP[AH[.
C [P[[HAK BO[MO[HO[ KOH[[[[PA[[[: IEND=0 .
      COMMON /BLDIS/ J(20),IEND,JZ(20),IENDZ
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      L=0
 1    L=L+1
      IF(L-K) 3,4,4
 4    IEND=1
      RETURN
 3    JL=J(L)
      JK=J(K)
      J(L)=J(L)+1
      J(K)=J(K)-1
      IF(J(L).GT.J(L+1).OR.J(K-1).GT.J(K)) GO TO 2
      IEND=0
      RETURN
 2    J(L)=1
      J(K)=JK+JL-1
      GO TO 1
      END
      
      
      
      SUBROUTINE CHOSAZ(A0,Z0,T,IP,PA,PAZ)
C B[[OP KAHA[A PA[BA[A [[ BCEX KOH[[[[PA[[[ [O M=4 .
C [A[O[HEH[E MACC[BA PA(1-IP) - B [OP[[KE [[[BAH[[ A.
      DIMENSION PA(500),PAZ(500)
      COMMON /BL500/W,W2(250),W3(25000),W4(10000)
      COMMON /BL501/IAF2(250),IAF3(25000),IAF4(10000)
      COMMON /BL502/M2,M3,M4,M5,MC2,MC3,MC4,MC5
      COMMON /BL503/WCOMP,WW2,WW3,WW4,WW5,WEVAP,WW4G
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      WCOWW2=WCOMP+WW2
      WCOWW3=WCOWW2+WW3
      BR=RRAN(ix)
      IP=1
      IF(BR.GT.WCOMP.AND.BR.LE.WCOWW2) IP=2
      IF(BR.GT.WCOWW2.AND.BR.LE.WCOWW3) IP=3
      IF(BR.GT.WCOWW3.AND.WW4.GT.0.) IP=4
      GO TO (1,2,3,4),IP
 1    PA(1)=A0
      PAZ(1)=Z0
      RETURN
 2    WP=WCOMP
      JJ=0
      DO 10 I=1,M2
      WP=WP+W2(I)
      JJ=JJ+1
      IF(WP.GT.BR) GO TO 11
 10   CONTINUE
 11   MCH=JJ
      PA(2)=IAF2(MCH)/1000
      PA(1)=IAF2(MCH)-INT(PA(2))*1000
      GO TO 5
 3    WP=WCOMP+WW2
      JJ=0
      DO 12 I=1,M3
      WP=WP+W3(I)
      JJ=JJ+1
      IF(WP.GT.BR) GO TO 13
 12   CONTINUE
 13   MCH=JJ
      PA(3)=IAF3(MCH)/1000000
      PA(2)=(IAF3(MCH)-INT(PA(3))*1000000)/1000
      PA(1)=IAF3(MCH)-INT(PA(3))*1000000-INT(PA(2))*1000
      GO TO 5
 4    WP=WCOMP+WW2+WW3
      JJ=0
      DO 14 I=1,M4
      JJ=JJ+1
      WP=WP+W4(I)
      IF(WP.GT.BR) GO TO 15
 14   CONTINUE
 15   MCH=JJ
      PA(4)=IAF4(MCH)/10000000
      PA(3)=(IAF4(MCH)-INT(PA(4))*10000000)/100000
      PA(2)=(IAF4(MCH)-INT(PA(4))*10000000-INT(PA(3))*100000)/1000
      PA(1)=IAF4(MCH)-INT(PA(4))*10000000-INT(PA(3))*100000
     *-INT(PA(2))*1000
 5    CONTINUE
      CALL CHOSZ(A0,Z0,T,IP,PA,PAZ)
      RETURN
      END
      
      
      
      SUBROUTINE CHOSZ(A0,Z0,T,IP,PA,PAZ)
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      DIMENSION PA(500),PAZ(500)
C O[PE[E[EH[E [AP[[OB IP [PA[MEHTOB B MACC[BE PAZ(1-IP)
      IZ0=INT(Z0+0.5)
      G0=25.
 5    ISZ=0
      DO 1 I=1,IP
      IA=INT(PA(I)+0.5)
 3    CALL RANNOR(BR1,BR2)
      CC=8.*G0
      ZM=PA(I)*Z0/A0
      IF(PA(I).GT.1.5.AND.PA(I).LT.4.5) ZM=0.5*PA(I)
      DZ=SQRT(PA(I)*T/CC)
      Z=BR1*DZ+ZM
      IZ=INT(Z+0.5)
      IF(IZ.GE.0.AND.IZ.LE.IA) GO TO 4
      Z=BR2*DZ+ZM
      IZ=INT(Z+0.5)
      IF(IZ.GE.0.AND.IZ.LE.IA) GO TO 4
      GO TO 3
 4    PAZ(I)=IZ
      ISZ=ISZ+INT(PAZ(I)+0.5)
 1    CONTINUE
      DELZ=IZ0-ISZ
      IF(ABS(DELZ).GT.1.1) GO TO 5
      PAZ(1)=PAZ(1)+DELZ
      RETURN
      END
      
      
      
      SUBROUTINE PARAM(T,K,KC,EPSIL0)
      COMMON /BLOK1/ W0,TC,RADNCL,G0
      COMMON /BLOK0/ A0,Z0,V0,XZ0,IA0,E0,E0Q
      COMMON /BLOK3/ E11(500),XZ(500),GA(500)
      COMMON /BLOK2/ B0,BT,DBT,LT,CP,PKP013
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      REAL*4 K,KC           ! sobol, 26.02.00
      W0=16.
      B0=18.
      TC=18.
c      TC=1000000.
      G0=25.
      RADNCL=1.17
      IF(T.LT.TC) GO TO 81
      BT=0.
      DBT=0.
      GO TO 82
 81   CONTINUE
      BT=B0*((TC*TC-T*T)/(TC*TC+T*T))**1.25
      DBT=B0*(-5*T*TC*TC/(TC*TC+T*T)**2)*((TC*TC-T*T)/(TC*TC+T*T))**.25
 82   CONTINUE
      IF(EPSIL0.GT.9999.) BT=B0
      IF(EPSIL0.GT.9999.) DBT=0.
      RETURN
      END
      
      
      
      SUBROUTINE SCOMPA(E00,K,KC,EPSIL0,SCOMP)
C  O[PE[E[EH[E [HTPO[[[ KOM[A[H[ COCTO[H[[
      COMMON /BLOK1/W0,TC,RADNCL,G0 /BLOK0/A0,Z0,V0,XZ0,IA0,E0,E0Q
      COMMON /BLOK3/E11(500),XZ(500),GA(500) /BLOKT/TCON
      COMMON /BLOK2/B0,BT,DBT,LT,CP,PKP013 /BLMOM/AMOM,AM(500),IAMO
      COMMON /BLERR/ERRMI,ERRMA,ERRFR
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      REAL K,KC,LT
      E0=E00*A0
C      EMOM=41.5*(AMOM**2)/(2*AM(IA0))
C      IF((E0-EMOM).LT.0.003.AND.IAMO.EQ.1) GO TO 73
      A013=A0**0.333333
      T=SQRT(E00/0.125)
      IF(T.LT.0.0012) T=0.0012
      HT=0.5
      K1=0
      K2=0
      ID=0
      ICNT=0
 24   CONTINUE
      IF(ICNT.GT.120) GO TO 71
      CALL PARAM(T,K,KC,EPSIL0)
      ECON=(-W0+T*T/E11(IA0)+G0*(1.-2.*XZ(IA0))**2)*A0+
     *(BT-T*DBT)*A013*A013+0.6*1.44*Z0*Z0/(RADNCL*A013)
C      IF(IAMO.EQ.1) ECON=ECON+EMOM
      D=(E0+E0Q-ECON)/E0
      IF(ABS(D).LT.0.003) GO TO 29
      ICNT=ICNT+1
      H=SIGN(HT,D)
      IF(D) 21,21,22
 21   K1=1
      IF(K1.EQ.1.AND.K2.EQ.1) ID=ID+1
      IF(ID.GT.29) GO TO 71
      T=T+H/2**ID
      IF(T.GE.0.001) GO TO 24
      K2=1
      H=HT
      GO TO 21
 22   K2=1
      IF(K1.EQ.1.AND.K2.EQ.1) ID=ID+1
      IF(ID.GT.29) GO TO 71
      T=T+H/2**ID
      IF(T.GE.0.001) GO TO 24
      K1=1
      H=HT
      GO TO 22
 71   ERRMI=ERRMI+1
C      PRINT 72,ICNT,ID,D,T
C 72   FORMAT(5X,'ERROR IN SCOMPA  ICNT=',I3,' ID=',I3,' D=',E12.5,
C     *'  T=',F8.3)
 73   SCOMP=0.
      RETURN
 29   TCON=T
      SCOMP=2*T*A0/E11(IA0)-DBT*A013*A013
      RETURN
      END
      
      
      
      SUBROUTINE FRAGTE(AA,ZZ,E00,K,KC,EPSIL0,TK,MU,L,SR)
C HAXO[[EH[E TEM[EPAT[P[ TK [ BCE[O M[[[T[[PA[MEHTHO[O COCTO[H[[ N(A),XZ
C [O [HEP[[[ BO[[[[[EH[[ (HA H[K[OH) E00 [[PA A0,Z0 [ [APAM.CB.O[[EMA K.
      REAL N,MU,K,LT,L,KC
      COMMON /BLOK1/W0,TC,RADNCL,G0 /BLOK0/A0,Z0,V0,XZ0,IA0,E0,E0Q
      COMMON /BLOK2/B0,BT,DBT,LT,CP,PKP013
      COMMON /BLOK3/E11(500),XZ(500),GA(500)
      COMMON /BLERR/ERRMI,ERRMA,ERRFR
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      ERRMA=0.
      A0=AA
      Z0=ZZ
      IA0=INT(A0+0.5)
      SM=(1.+2.31*(E00-3.5))*A0/100.
      IF(SM.LT.2.) SM=2.
      K=(1.+1.44*(SM**0.333333-1.)/(1.17*A0**0.333333))**3-1.
      E11(1)=0.
      DO 100 I=2,IA0
      A=I
C 100  E11(I)=EPSIL0*(1.+3./(A-1))
C      IF(I.LE.50) E11(I)=EPSIL0*(1.+0.12*(A**0.33333-5.**0.33333))
C 100  IF(I.GT.50) E11(I)=EPSIL0*(1.+0.12*(50.**0.33333-5.**0.33333))
 100  E11(I)=EPSIL0/(1.+0.0020*(A/25.)**2)
      T=SQRT(E00/0.12)
      IF(T.LT.0.0012) T=0.0012
      HT=1.0
      ICNT=0
      K1=0
      K2=0
      ID=0
 24   CONTINUE
      IF(ICNT.GT.60) GO TO 71
      CALL FRAG(A0,Z0,T,K,KC,EPSIL0,MU,L,E0,SR)
C *********
C      PRINT 110,E00,E0,T
C 110  FORMAT(2X,'**FRAGTE** E00=',F11.5,' E0=',F11.5,' T=',F11.5)
C *********
      D=(E00-E0)/E00
      IF(ABS(D).LT.0.003) GO TO 29
      ICNT=ICNT+1
      H=SIGN(HT,D)
      IF(D) 21,21,22
 21   K1=1
      IF(K1.EQ.1.AND.K2.EQ.1) ID=ID+1
      IF(ID.GT.30) GO TO 71
      T=T+H/2**ID
      IF(T.GE.0.001) GO TO 24
      K2=1
      H=HT
      GO TO 21
 22   K2=1
      IF(K1.EQ.1.AND.K2.EQ.1) ID=ID+1
      IF(ID.GT.30) GO TO 71
      T=T+H/2**ID
      IF(T.GE.0.001) GO TO 24
      K1=1
      H=HT
      GO TO 22
 29   TK=T
      RETURN
 71   ERRMA=ERRMA+1
C      PRINT 72,ICNT,D,T,ID
C 72   FORMAT(5X,'ERROR IN FRAGTE   ICNT=',I3,' D=',E12.5,' T=',F8.3,
C     *'             ID=',I2)
      TK=T
      RETURN
      END
      
      
      
      SUBROUTINE FINDT(E00,PA,PAZ,IP,TK,ETP,ECOL,ECP,EEX,EB)
C HAXO[[EH[E [O IP B[[PAHH[M [PA[MEHTAM [X TEM[EPAT[P[ T [ [HEP[[[,
C C[[TA[,[TO [O[HA[ [HEP[[[ BO[[[[[EH[[ C[CTEM[ E0*A0 - [[KC[POBAHA.
      REAL N,LT
      COMMON /BLOK1/W0,TC,RADNCL,G0 /BLOK0/A0,Z0,V0,XZ0,IA0,E0,E0Q
      COMMON /BLOK2/B0,BT,DBT,LT,CP,PKP013
      COMMON /BLOK3/E11(500),XZ(500),GA(500)
      COMMON /BLERR/ERRMI,ERRMA,ERRFR
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      DIMENSION PA(500),PAZ(500)
      ERRFR=0.
      A013=A0**0.333333
      E0=E00
      HT=0.5
      ICNT=0
      K1=0
      K2=0
      ID=0
      T=TK
      IF(T.LT.0.0012) T=0.0012
      CALL DELAM(A0,Z0,DLM0,DSH0,BAR)
 24   CONTINUE
      IF(ICNT.GT.120) GO TO 71
      EEX=0.
      EB=-DLM0
      ECP=0.
      ETP=1.5*T*IP
      DO 107 I=1,IP
      CALL DELAM(PA(I),PAZ(I),DLM,DSH,BAR)
      EB=EB+DLM
      PA03=PA(I)**0.333333
      ECP=ECP+(0.6*1.44/RADNCL)*PAZ(I)*PAZ(I)/PA03
      IA=INT(PA(I)+0.1)
      IF(IA.LE.3) GO TO 107
      ESUF=PA(I)*T*T*2.5*B0/(TC*TC*PA03)
      IF(DBT.EQ.0.) ESUF=0.
      EEX=EEX+PA(I)*T*T/E11(IA)+ESUF
      IF(IA.EQ.4) EEX=EEX-ESUF
 107  CONTINUE
      ECOL=(0.6*1.44/RADNCL)*Z0*Z0*PKP013/A013-ECP*PKP013
      E0TP=ETP+EEX+EB+ECOL
      E0TPA0=E0TP/A0
C ********
C      PRINT 110,E0,E0TPA0,T,ETP,EEX,EB,ECOL
C 110  FORMAT(2X,'*FINDT* E0=',F11.5,' E0TP/A0=',F11.5,' T=',F11.5,3X,
C     *'ETP,EEX,EB,ECOL=',4F11.5)
C ********
      D=(E0-E0TP/A0)/E0
      IF(ABS(D).LT.0.003) GO TO 29
      ICNT=ICNT+1
      H=SIGN(HT,D)
      IF(D) 21,21,22
 21   K1=1
      IF(K1.EQ.1.AND.K2.EQ.1) ID=ID+1
      IF(ID.GT.30) GO TO 71
      T=T+H/2**ID
      IF(T.GE.0.001) GO TO 24
      K2=1
      H=HT
      GO TO 21
 22   K2=1
      IF(K1.EQ.1.AND.K2.EQ.1) ID=ID+1
      IF(ID.GT.30) GO TO 71
      T=T+H/2**ID
      IF(T.GE.0.001) GO TO 24
      K1=1
      H=HT
      GO TO 22
 29   TK=T
      RETURN
 71   ERRFR=ERRFR+1
C      PRINT 72,ICNT,ID,D,T,IP
C 72   FORMAT(5X,'ERROR IN FINDT    ICNT=',I3,' ID=',I3,' D=',E12.5,
C     *' T=',F9.3,'  IP=',I3)
      RETURN
      END
      
      
      
      FUNCTION UEVA(A,T)
C O[PE[E[EH[E [HEP[[[ BO[[[[[EH[[ [PA[MEHTOB [EPE[ [X [C[APEH[EM [[[
C [EPM[-PA[BA[OM (UEVA B [[B, TEM[EPAT[PA T B M[B).
      COMMON /BLOK2/B0,BT,DBT,LT,CP,PKP013 /BLOK1/W0,TC,RADNCL,G0
      COMMON /BLOK3/E11(500),XZ(500),GA(500)
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      IA=INT(A+0.5)
      ESUF=A*T*T*2.5*B0/(TC*TC*A**0.333333)
      IF(DBT.EQ.0.) ESUF=0.
      UEVA=A*T*T/E11(IA)+ESUF
      IF(IA.EQ.4) UEVA=A*T*T/E11(IA)
      UEVA=0.001*UEVA
      RETURN
      END
      
      
      
      SUBROUTINE FRAG(AA,ZZ,T,K,KC,EPSIL0,MU,L,EE,SR)
C [PA[MEHTA[[[ [[PA A0,Z0 [P[ [AHHO[ TEM[EPAT[PE T [ [APAM.CB. O[[EMA K.
C [APAMETP[ E11(A) [A[A[TC[ PAHEE B /BLOK2/. HA B[XO[E: X[M.[OT. MU [ L,
C [HEP[[[ E0(M[B/N),[HTPO[[[ S (1/N),CPE[H.[[C[O [PA[MEHTOB SR; MHO[ECTB
C HOCT[ N(A) [ [AP[[OBOE PAC[PE[E[EH[E XZ(A) - HAXO[[TC[ B /BLOK2/.
      REAL N,MU,K,LT,L,KC
      COMMON /BLOK2/B0,BT,DBT,LT,CP,PKP013 /BLOKN/N(500) /BENTRO/SSR
      COMMON /BLOK1/W0,TC,RADNCL,G0 /BLOK0/A0,Z0,V0,XZ0,IA0,E0,E0Q
      COMMON /BLOK3/E11(500),XZ(500),GA(500) /BLOKE/E(500)
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      A0=AA
      Z0=ZZ
      TK=T
      CALL PARAM(T,K,KC,EPSIL0)
      A013=A0**0.33333333
      A023=A013*A013
      R0=RADNCL*A013
      V0=(4.*3.1415926/3.)*R0*R0*R0
      P0=A0/V0
      R=R0*(1.+KC)**0.33333333
      PKP013=R0/R
      XZ0=Z0/A0
      E0Q=A0*(-W0+G0*(1.-2.*XZ0)**2)+B0*A023+0.6*1.44*Z0*Z0/R0
      P=P0/(1.+KC)
      CP=(0.6*1.44/RADNCL)*(1.-PKP013)
      LT=16.15/SQRT(TK)
      CALL XZA(TK,K,MU,L)
      CALL ENERGA(TK,K)
      CALL ENTROP(TK,K,S)
      SSR=S
      E0T=0.
      DO 10 IT=1,IA0
      E0T=E0T+N(IT)*E(IT)
   10 CONTINUE
      E0T=E0T+0.6*1.44*Z0*Z0/R
      E0=(E0T-E0Q)/A0
      EE=E0
      SR=0.
      DO 31 J=1,IA0
      SR= SR+N(J)
 31   CONTINUE
      RETURN
      END
      
      
      
      SUBROUTINE CALCMU(T,K,MU,L)
C B[[[C[EH[E MU=MU(T) [C[O[[[[[ [PABHEH[E  N(1)*1+...+N(A0)*A0=A0
      REAL N,MU,K,LT,L
      COMMON /BLOK1/W0,TC,RADNCL,G0 /BLOK0/A0,Z0,V0,XZ0,IA0,E0,E0Q
      COMMON /BLOK2/B0,BT,DBT,LT,CP,PKP013    /BLOKN/N(500)
      COMMON /BLOK3/E11(500),XZ(500),GA(500)
      COMMON /BLERR/ERRMI,ERRMA,ERRFR
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      HM=1.
      MU=-W0-T*T/E11(5)-L*XZ(5)+G0*(1.-2*XZ(5))**2
     *+(2./3.)*BT/1.71+5.*CP*XZ(5)*XZ(5)*2.92/3.-1.5*T/5.
      KK=0
      K1=0
      K2=0
      ID=0
 24   CONTINUE
      CALL MULTIA(T,K,MU,L)
      A0T=0.
      DO 10 IT=1,IA0
      A0T=A0T+IT*N(IT)
 10   CONTINUE
      KK=KK+1
      D=(A0-A0T)/A0
      IF(ABS(D).LT.0.001) GO TO 25
      IF(KK.GT.60.OR.ID.GT.30) GO TO 61
      H=SIGN(HM,D)
      IF(D) 21,21,22
 21   K1=1
      IF(K1.EQ.1.AND.K2.EQ.1) ID=ID+1
      MU=MU+H/2**ID
      GO TO 24
 22   K2=1
      IF(K1.EQ.1.AND.K2.EQ.1) ID=ID+1
      MU=MU+H/2**ID
      GO TO 24
 25   CONTINUE
      RETURN
   61 ERRMA=ERRMA+1
C      PRINT 62,KK,ID,D,MU,T
C   62 FORMAT(5X,'ERROR IN CALCMU  KK=',I2,' ID=',I2,' D=',E12.5,' MU=',
C     *F8.3,' T=',F8.3)
      RETURN
       END
      
      
      
      SUBROUTINE MULTIA(T,K,MU,L)
C B[[[C[EH[E MHO[ECTBEHHOCTE[ [PA[MEHTOB - N(A)
      REAL N,MU,K,LT,L
      COMMON /BLOK1/W0,TC,RADNCL,G0 /BLOK0/A0,Z0,V0,XZ0,IA0,E0,E0Q
      COMMON /BLOK2/B0,BT,DBT,LT,CP,PKP013   /BLOKN/N(500)
      COMMON /BLOK3/E11(500),XZ(500),GA(500) /BLOKY/YN,YP,Y2,Y3,Y4
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      YN=(2*K*V0/(LT*LT*LT))*EXP(MU/T)
      YP=(2*K*V0/(LT*LT*LT))*EXP((MU+L-CP)/T)
      N(1)=YN+YP
      Y2=3*2*SQRT(2.)*K*V0/(LT**3)
      Y2=Y2*EXP((2.796+2*(MU+L*XZ(2))-CP*XZ(2)*XZ(2)*(2.**(5./3.)))/T)
      Y3=4*3*SQRT(3.)*K*V0/(LT**3)
      Y3=Y3*EXP((9.224+3*(MU+L*XZ(3))-CP*XZ(3)*XZ(3)*(3.**(5./3.)))/T)
      Y4=8*K*V0/(LT**3)
      Y4=Y4*EXP((30.11+4*(MU+L*XZ(4)+T*T/E11(4))-CP*XZ(4)*XZ(4)*
     #(4.**(5./3.)))/T)
      N(2)=Y2
      N(3)=Y3
      N(4)=Y4
      DO 2 I=5,IA0
      A=I
      A23=A**0.66666667
      VNE=(MU+L*XZ(I)+W0+T*T/E11(I)-G0*(1.-2*XZ(I))**2)*A
     *-BT*A23-CP*XZ(I)*XZ(I)*A*A23
      VNE=VNE/T
      IF(VNE.GT.30.0) GO TO 1
      VN=EXP(VNE)
      N(I)=(K*V0*SQRT(A)*A/(LT*LT*LT))*VN
      IF(N(I).LT.1.0E-30) N(I)=0.
      GO TO 2
 1    N(I)=999.
 2    CONTINUE
      RETURN
      END
      
      
      
      SUBROUTINE ENERGA(T,K)
C B[[[C[EH[E [HEP[[[ [PA[MEHTOB  E(A)
      REAL N,MU,K,LT,L
      COMMON /BLOK1/W0,TC,RADNCL,G0 /BLOK0/A0,Z0,V0,XZ0,IA0,E0,E0Q
      COMMON /BLOK2/B0,BT,DBT,LT,CP,PKP013    /BLOKE/E(500)
      COMMON /BLOK3/E11(500),XZ(500),GA(500)
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      E(1)=1.5*T+CP*XZ(1)
      E(2)=-2.796+CP*XZ(2)*XZ(2)*(2.**(5./3.))+1.5*T
      E(3)=-9.224+CP*XZ(3)*XZ(3)*(3.**(5./3.))+1.5*T
      E(4)=-30.11+CP*XZ(4)*XZ(4)*(4.**(5./3.))+1.5*T+4*T*T/E11(4)
      DO 2 I=5,IA0
      A=I
      A23=A**0.66666667
      EV=A*(T*T/E11(I)-W0+G0*(1.-2*XZ(I))**2)
      ES=(BT-T*DBT)*A23
      EC=CP*A23*A*XZ(I)*XZ(I)
      ET=1.5*T
      E(I)=EV+ES+EC+ET
 2    CONTINUE
      RETURN
      END
      
      
      
      SUBROUTINE ENTROP(T,K,S)
C B[[[C[EH[E [HTPO[[[ [[PA [OC[E [PA[MEHTA[[[ - S
      REAL N,MU,K,LT,L
      COMMON /BLOK1/W0,TC,RADNCL,G0 /BLOK0/A0,Z0,V0,XZ0,IA0,E0,E0Q
      COMMON /BLOK2/B0,BT,DBT,LT,CP,PKP013    /BLOKN/N(500)
      COMMON /BLOK3/E11(500),XZ(500),GA(500)  /BLOKY/YN,YP,Y2,Y3,Y4
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      S=0.
      IF(YN.GT.0.) S=YN*(2.5+ALOG(2*K*V0/(LT*LT*LT*YN)))
      IF(YP.GT.0.) S=S+YP*(2.5+ALOG(2*K*V0/(LT*LT*LT*YP)))
      IF(Y2.GT.0.) S=S+Y2*(2.5+ALOG(3*K*V0*2*SQRT(2.)/(Y2*LT**3)))
      IF(Y3.GT.0.) S=S+Y3*(2.5+ALOG(4*K*V0*3*SQRT(3.)/(Y3*LT**3)))
      IF(Y4.GT.0.) S=S+Y4*(2.5+ALOG(8*K*V0/(Y4*LT**3))+8*T/E11(4))
      DO 1 I=5,IA0
      A=I
      IF(N(I).LE.0.) GO TO 1
      SV=A*2*T/E11(I)
      SS=-DBT*A**0.6666667
      ST=2.5+ALOG(K*V0*SQRT(A)*A/(LT*LT*LT*N(I)))
      S=S+(SV+SS+ST)*N(I)
 1    CONTINUE
      RETURN
      END
      
      
      
      SUBROUTINE XZA(T,K,MU,L)
C B[[[C[EH[E XZ(A) - [O[[ Z B [PA[MEHTE A,[C[O[[[[[ [PABHEH[E
C                 XZ(1)*1*N(1)+...+XZ(A0)*A0*N(A0)=Z0  ([P[[[[[EHHO)
      REAL N,MU,K,LT,L
      COMMON /BLOK1/W0,TC,RADNCL,G0 /BLOK0/A0,Z0,V0,XZ0,IA0,E0,E0Q
      COMMON /BLOK2/B0,BT,DBT,LT,CP,PKP013     /BLOKN/N(500)
      COMMON /BLOK3/E11(500),XZ(500),GA(500)   /BLOKY/YN,YP,Y2,Y3,Y4
      COMMON /BLERR/ERRMI,ERRMA,ERRFR
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      DIMENSION C1(500),C2(500)
      K1=0
      K2=0
      KK=0
      ID=0
      HM=1.
      DO 1 I=5,IA0
      A=I
      CC=8.*G0+2.*CP*A**0.66666667
      C1(I)=4.*G0/CC
      C2(I)=1./CC
 1    CONTINUE
      L=XZ0*(8.*G0+2.*CP*A0**0.6666667)-4.*G0
 6    CONTINUE
      XZ(2)=0.5
      XZ(3)=0.5
      XZ(4)=0.5
      DO 2 I=5,IA0
 2    XZ(I)=C1(I)+L*C2(I)
      CALL CALCMU(T,K,MU,L)
      XZ(1)=YP/(YP+YN)
      Z0T=0.
      DO 3 I=1,IA0
 3    Z0T=Z0T+XZ(I)*I*N(I)
      KK=KK+1
      D=(Z0-Z0T)/Z0
      IF(ABS(D).LT.0.002) GO TO 9
      IF(KK.GT.60.OR.ID.GT.30) GO TO 7
      H=SIGN(HM,D)
      IF(D) 4,4,5
 4    K1=1
      IF(K1.EQ.1.AND.K2.EQ.1) ID=ID+1
      L=L+H/2**ID
      GO TO 6
 5    K2=1
      IF(K1.EQ.1.AND.K2.EQ.1) ID=ID+1
      L=L+H/2**ID
      GO TO 6
 9    CONTINUE
      RETURN
 7    ERRMA=ERRMA+1
C      PRINT 8,KK,ID,D,L,MU,T
C 8    FORMAT(5X,'ERROR IN XZA KK=',I2,' ID=',I2,' D=',E12.5,' L=',F8.3,
C     *' MU=',F8.3,' T=',F8.3)
      RETURN
      END
      
      
      
      SUBROUTINE RANNOR(A,B)
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      Y=RRAN(ix)
      Z=RRAN(ix)
      X=6.283185*Z
      A1=SQRT(-2.0*ALOG(Y))
      A=A1*SIN(X)
      B=A1*COS(X)
      RETURN
      END
      
      
      
      SUBROUTINE SELECT(A0,Z0,N,T,L,IP,PA,PAZ)
C C[[[A[HOE PA[[[EH[E [[PA (A0,Z0) HA [PA[MEHT[,B COOTBETCTB[[ C [HK[[[[
C PAC[PE[E[EH[EM N(A) [ HOPMA[[H[M PAC[PE[E[EH[EM [PA[MEHTOB [O Z.
C HA B[XO[E: [[C[O [PA[M.-IP,[X MACC[ -PA(1-IP) /B [OP[[KE [[[BAH[[ A/,
C [ [AP[[[ -PAZ(1-IP).
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      INTEGER NA(500)
      REAL*4 N(500),L
      DIMENSION PA(500),PAZ(500)
      IA0=INT(A0+0.5)
      IZ0=INT(Z0+0.5)
      CALL SELECA(N,IA0,NA,IP)
C  [A[O[HEH[E MACC[BA PA ( [[ MACC[BA NA )
      IS=0
      DO 305 I=1,IP
      PA(I)=0.
 305  PAZ(I)=0.
      DO 306 IA=1,IA0
      IF(NA(IA).EQ.0) GO TO 306
      JIA=NA(IA)
      DO 307 J=1,JIA
 307  PA(IS+J)=IA
      IS=IS+JIA
 306  CONTINUE
C [EPEH[MEPA[[[ MACC[BA PA B [OP[[KE [[[BAH[[ A
      DO 10 J=1,IP
      PAMAX=0.
      DO 11 I=J,IP
      IF(PA(I)-PAMAX) 11,11,12
 12   IM=I
      PAMAX=PA(IM)
 11   CONTINUE
      PA(IM)=PA(J)
      PA(J)=PAMAX
 10   CONTINUE
      CALL SELECZ(T,L,IZ0,IP,PA,PAZ)
      RETURN
      END
      
      
      
      SUBROUTINE SELECZ(T,L,IZ0,IP,PA,PAZ)
      COMMON /BLOK2/B0,BT,DBT,LT,CP,PKP013 /BLOKY/YN,YP,Y2,Y3,Y4
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      REAL*4 PA(500),PAZ(500),L,LT
C O[PE[E[EH[E [AP[[OB IP [PA[MEHTOB: B MACC[BE PAZ(1-IP). MACC[ [PA[MEHT
C [AH[ B MACC[BE PA(1-IP),B [OP[[KE [[[BAH[[ A. [AP[[[ [PA[[C[ [[ HOPMA[
C PAC[PE[E[EH[[ C [EHTPOM ZM [ [[C[EPC[E[ DZ. OCTATOK PO[[[P[[A ( OT[[[[
C C[MM[ PAZ OT Z0 HE [O[EE 0,+1,-1 ) [O[AB[[[C[ B PA(1).
      G0=25.
 5    ISZ=0
      DO 1 I=1,IP
      IA=INT(PA(I)+0.5)
      IF(IA-1) 2,2,3
 2    BR=RRAN(ix)
      PAZ(I)=0.
      IF(BR.GT.(YN/(YN+YP))) PAZ(I)=1.
      ISZ=ISZ+INT(PAZ(I)+0.5)
      GO TO 1
 3    CALL RANNOR(BR1,BR2)
      CC=8.*G0+2.*CP*PA(I)**0.666667
      ZM=PA(I)*(4.*G0+L)/CC
      IF(PA(I).GT.1.5.AND.PA(I).LT.4.5) ZM=0.5*PA(I)
      DZ=SQRT(PA(I)*T/CC)
      Z=BR1*DZ+ZM
      IZ=INT(Z+0.5)
      IF(IZ.GE.0.AND.IZ.LE.IA) GO TO 4
      Z=BR2*DZ+ZM
      IZ=INT(Z+0.5)
      IF(IZ.GE.0.AND.IZ.LE.IA) GO TO 4
      GO TO 3
 4    PAZ(I)=IZ
      ISZ=ISZ+INT(PAZ(I)+0.5)
 1    CONTINUE
      DELZ=IZ0-ISZ
      IF(ABS(DELZ).GT.1.1) GO TO 5
      PAZ(1)=PAZ(1)+DELZ
      RETURN
      END
      
      
      
      SUBROUTINE SELECA(N,IA0,NA,IP)
C C[[[A[H[[ B[[OP [[C[A [PA[MEHTOB IP [ [X MHO[ECTBEHHOCTE[ [O A -NA(A)
C [[ [[C[A H[K[OHOB IA0 [ [HK[[[[BH[X MHO[ECTBEHHOCTE[ -N(A).
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      REAL GG(500),N(500)
      INTEGER NA(500)
      IMAXMA=0
      IMAXMI=0
      GH=0.
      DO 2 I=1,IA0
 2    GH=GH+N(I)
      SQGH=SQRT(GH)+0.5
      GG(1)=N(1)
      DO 3 I=2,IA0
 3    GG(I)=GG(I-1)+N(I)
 4    ISA=0
      DO 1 I=1,IA0
 1    NA(I)=0
      IP=0
 5    BRND=RRAN(ix)*GH
      DO 6 I=1,IA0
      IF(BRND-GG(I)) 7,6,6
 7    IA=I
      GO TO 8
 6    CONTINUE
 8    IP=IP+1
      NA(IA)=NA(IA)+1
      ISA=ISA+IA
      IAK=IA0-ISA
      IF(IAK.EQ.0) GO TO 10
      IF(IAK.GT.IMAXMA) GO TO 5
      IF(IAK.LT.IMAXMI) GO TO 4
      IP=IP+1
      NA(IAK)=NA(IAK)+1
 10   CONTINUE
      IF(ABS(GH-IP).GT.SQGH) GO TO 4
      RETURN
      END
      
      
      
      SUBROUTINE DISNET(INET,IP0,AMP,PN0,T,TN)
C HAXO[[EH[E [HEP[[[ [ [M[[[[COB INET HE[TPOHOB [O [X TEM[EPAT[PE T(MEV)
C [O[HO[ K[H.[HEP[[[ TN(GEV),C[[TA[ PAC[PE[E[EH[E MAKCBE[OBCK[M. AMP(IP0
C -[X MACC[ (B GEV);[M[[[[C[ (GEV/C),B C[MME [A[[[E 0,[O[AB[[[TC[ B MACC
C PN0(3,500) C (IP0+1).
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      DIMENSION AMP(500),PN0(3,500),ANL(3),PX(200),PY(200),PZ(200)
      DIMENSION PI(3),PJ(3),PN(3),AR(3),BR(3)
      IF(INET.LE.0) RETURN
      IP01=IP0+1
      IP0M=IP0+INET
      IF(INET-1) 1,1,2
 1    P=SQRT(2.*AMP(IP01)*TN)
      CALL ISOTR(ANL)
      PN0(1,IP01)=P*ANL(1)
      PN0(2,IP01)=P*ANL(2)
      PN0(3,IP01)=P*ANL(3)
      RETURN
 2    IF(INET-2) 3,3,4
 3    P=SQRT(2.*(AMP(IP01)*AMP(IP0M)/(AMP(IP01)+AMP(IP0M)))*TN)
      CALL ISOTR(ANL)
      PN0(1,IP01)=P*ANL(1)
      PN0(2,IP01)=P*ANL(2)
      PN0(3,IP01)=P*ANL(3)
      PN0(1,IP0M)=-PN0(1,IP01)
      PN0(2,IP0M)=-PN0(2,IP01)
      PN0(3,IP0M)=-PN0(3,IP01)
      RETURN
 4    IP0M2=IP0M-2
      ES=0.
      PSX=0.
      PSY=0.
      PSZ=0.
      FEMT=SQRT(0.5*T)*EXP(-0.5)
      DO 5 I=IP01,IP0M2
 7    E=RRAN(ix)*9.*T
      FE=SQRT(E)*EXP(-E/T)
      FERAND=RRAN(ix)*FEMT
      IF(FERAND-FE) 6,6,7
 6    P=SQRT(2.*E*0.001*AMP(I))
      CALL ISOTR(ANL)
      PX(I)=P*ANL(1)
      PY(I)=P*ANL(2)
      PZ(I)=P*ANL(3)
      ES=ES+E*0.001
      PSX=PSX+PX(I)
      PSY=PSY+PY(I)
      PSZ=PSZ+PZ(I)
 5    CONTINUE
      I1=IP0M2+1
      I2=IP0M2+2
      PPX=-PSX
      PPY=-PSY
      PPZ=-PSZ
      P=SQRT(PPX*PPX+PPY*PPY+PPZ*PPZ)
      EE=TN-ES
      EM=P*P/(2.*(AMP(I1)+AMP(I2)))
      IF(EE.LE.EM) GO TO 4
      H=1.+AMP(I2)/AMP(I1)
      CTM12=H*(1.-2.*AMP(I2)*EE/(P*P))
 11   CT1=1.-2.*RRAN(ix)
      IF(CT1*CT1-CTM12) 11,11,12
 12   IF(CTM12) 13,17,17
 13   IZN=1
      GO TO 15
 17   IF(CT1) 11,14,14
 14   IF(RRAN(ix)-0.5) 16,16,13
 16   IZN=-1
 15   CONTINUE
      P1=(P*CT1+IZN*SQRT(P*P*CT1*CT1-P*P*CTM12))/H
      P2=SQRT(P1*P1+P*P-2.*P1*P*CT1)
      PHI=6.28318*RRAN(ix)
      ST1=SQRT(1.-CT1*CT1)
      CPHI1=COS(PHI)
      SPHI1=SIN(PHI)
      CPHI2=-CPHI1
      SPHI2=-SPHI1
      CT2=(P*P+P2*P2-P1*P1)/(2.*P*P2)
      IF(CT2.GT.-1.AND.CT2.LT.1) GO TO 20
      ST2=0.
      GO TO 21
 20   ST2=SQRT(1.-CT2*CT2)
 21   CONTINUE
      PI(1)=P1*ST1*CPHI1
      PI(2)=P1*ST1*SPHI1
      PI(3)=P1*CT1
      PJ(1)=P2*ST2*CPHI2
      PJ(2)=P2*ST2*SPHI2
      PJ(3)=P2*CT2
      AR(1)=PPX
      AR(2)=PPY
      AR(3)=PPZ
      BR(1)=1.
      BR(2)=0.
      BR(3)=0.
      CALL ROTORD(AR,BR,PI,PN)   ! ROTORD - 13.02.96
      PX(I1)=PN(1)
      PY(I1)=PN(2)
      PZ(I1)=PN(3)
      CALL ROTORD(AR,BR,PJ,PN)   ! ROTORD - 13.02.96
      PX(I2)=PN(1)
      PY(I2)=PN(2)
      PZ(I2)=PN(3)
      PSX=PSX+PX(I1)+PX(I2)
      PSY=PSY+PY(I1)+PY(I2)
      PSZ=PSZ+PZ(I1)+PZ(I2)
      ES=ES+(PX(I1)**2+PY(I1)**2+PZ(I1)**2)/(2.*AMP(I1))
     *+(PX(I2)**2+PY(I2)**2+PZ(I2)**2)/(2.*AMP(I2))
      DO 8 I=IP01,IP0M
      PN0(1,I)=PX(I)
      PN0(2,I)=PY(I)
      PN0(3,I)=PZ(I)
 8    CONTINUE
      RETURN
      END
      
      
      
c      SUBROUTINE ROTOR (AR,BR,PSTAR,PR)
cC--------[[OK [OBOPOTA BEKTOPA PSTAR [[ PA[O[E[ C.[.[. B [CXO[H[[ (PR)
cC--------AR, BR-BEKTOPA B [CXO[HO[ C.[.[.
cC    BLOCK OF ROTATION.
c      DIMENSION AR(3),BR(3),PSTAR(3),PR(3),AN(3)
c      SP = 0.
c      DO 31 IR=1,3
c      SP = SP+AR(IR)*BR(IR)
c   31 CONTINUE
c      AMOD = SQRT(AR(1)**2+AR(2)**2+AR(3)**2)
c      ALPHA1 = SP/AMOD
c      BMOD2 = BR(1)**2+BR(2)**2+BR(3)**2
c      ALPHA2 = SQRT(BMOD2-ALPHA1**2)
c      AN(1) = AR(2)*BR(3)-AR(3)*BR(2)
c      AN(2) = AR(3)*BR(1)-AR(1)*BR(3)
c      AN(3) = AR(1)*BR(2)-AR(2)*BR(1)
c      PR(1)=PSTAR(1)*BR(1)/ALPHA2+(PSTAR(3)-ALPHA1*PSTAR(1)/ALPHA2)
c     1*AR(1)/AMOD+(PSTAR(2)*AN(1))/(ALPHA2*AMOD)
c      PR(2)=PSTAR(1)*BR(2)/ALPHA2+(PSTAR(3)-ALPHA1*PSTAR(1)/ALPHA2)
c     1*AR(2)/AMOD+(PSTAR(2)*AN(2))/(ALPHA2*AMOD)
c      PR(3)=PSTAR(1)*BR(3)/ALPHA2+(PSTAR(3)-ALPHA1*PSTAR(1)/ALPHA2)
c     1*AR(3)/AMOD+(PSTAR(2)*AN(3))/(ALPHA2*AMOD)
c      RETURN
c      END
      
      
      
      SUBROUTINE CULIMP(IP0,A0,ZC0,AM,CH,PN0,TEMP,ECOLB)
      DIMENSION X0(99),Y0(99),Z0(99),VX0(99),VY0(99),VZ0(99),
     *VX(99),VY(99),VZ(99),AM(500),CH(500),VN(99),RS(99),
     *PN0(3,500),PNT(3,500)
      COMMON /BLPAT/PA(500),PAZ(500),IP  /BPLACE/RFI
      COMMON /BLCEN/XC(500),YC(500),ZC(500)  /BLDT/DT,IT
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
c     DT=10.
      DT=2.
      RN=1.17
      RN0=RN
      IF(RFI.GT.0.001) RN0=RFI
      RSYS=2.0*RN*(A0**0.333333)
      IP=IP0
      DO 50 I=1,IP
      PA(I)=AM(I)/0.94
 50   PAZ(I)=CH(I)
      CALL PLACE(RSYS,RN0)
      DO 1 I=1,IP
      X0(I)=XC(I)
      Y0(I)=YC(I)
      Z0(I)=ZC(I)
 1    CONTINUE
      ET=1.5*IP*TEMP*0.001
      CALL DISNET(IP,0,AM,PNT,TEMP,ET)
      DO 40 I=1,IP
      VX0(I)=PNT(1,I)/AM(I)
      VY0(I)=PNT(2,I)/AM(I)
      VZ0(I)=PNT(3,I)/AM(I)
 40   CONTINUE
      CALL CULON(X0,Y0,Z0,VX0,VY0,VZ0,VX,VY,VZ,AM,CH,VN,RS,IP,ECOLB,ET)
      DO 7 I=1,IP
      PN0(1,I)=AM(I)*VX(I)
      PN0(2,I)=AM(I)*VY(I)
      PN0(3,I)=AM(I)*VZ(I)
 7    CONTINUE
      RETURN
      END
      
      
      
      SUBROUTINE CULON(X0,Y0,Z0,VX0,VY0,VZ0,VX,VY,VZ,AM,CH,VN,RS,M,EC,
     *ET)
      DIMENSION X0(99),Y0(99),Z0(99),VX0(99),VY0(99),VZ0(99),
     *X(99),Y(99),Z(99),VX(99),VY(99),VZ(99),AM(500),CH(500),
     *R(99,99),F(99,99),FX(99,99),FY(99,99),FZ(99,99),
     *FSX(99),FSY(99),FSZ(99),AX(99),AY(99),AZ(99),VS(99),VN(99),RS(99),
     *VSX(99),VSY(99),VSZ(99)
      COMMON /BLDT/DT,IT
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      IT=0
      TN=0.0
      TS=0.0
      DO 9 I=1,M
      AX(I)=0.0
      AY(I)=0.0
      AZ(I)=0.0
      VS(I)=SQRT(VX0(I)*VX0(I)+VY0(I)*VY0(I)+VZ0(I)*VZ0(I))
      VX(I)=VX0(I)
      VY(I)=VY0(I)
      VZ(I)=VZ0(I)
      X(I)=X0(I)
      Y(I)=Y0(I)
      Z(I)=Z0(I)
      DO 9 J=1,M
      F(I,J)=0.0
      R(I,J)=0.0
      FX(I,J)=0.0
      FY(I,J)=0.0
      FZ(I,J)=0.0
 9    CONTINUE
      IF(EC.LE.0.) RETURN
 10   CONTINUE
      DO 1 I=1,M
      DO 1 J=1,M
      IF(I.EQ.J) GO TO 1
      R(I,J)=SQRT((X(I)-X(J))*(X(I)-X(J))+(Y(I)-Y(J))*(Y(I)-Y(J))+
     *(Z(I)-Z(J))*(Z(I)-Z(J)))
 1    CONTINUE
      ECULT=0.0
      DO 2 I=1,M
      DO 2 J=1,M
      IF(I.EQ.J) GO TO 2
      F(I,J)=1.44*(CH(I)*CH(J))/(R(I,J)*R(I,J))
      ECULT=ECULT+(1.44*0.5)*(CH(I)*CH(J))/R(I,J)
 2    CONTINUE
      DO 3 I=1,M
      DO 3 J=1,M
      IF(I.EQ.J) GO TO 3
      FX(I,J)=F(I,J)*(X(I)-X(J))/R(I,J)
      FY(I,J)=F(I,J)*(Y(I)-Y(J))/R(I,J)
      FZ(I,J)=F(I,J)*(Z(I)-Z(J))/R(I,J)
 3    CONTINUE
      DO 5 I=1,M
      FSX(I)=0.0
      FSY(I)=0.0
 5    FSZ(I)=0.0
      DO 4 I=1,M
      DO 4 J=1,M
      IF(I.EQ.J) GO TO 4
      FSX(I)=FSX(I)+FX(I,J)
      FSY(I)=FSY(I)+FY(I,J)
      FSZ(I)=FSZ(I)+FZ(I,J)
 4    CONTINUE
C     PRINT 21,((FSX(I),FSY(I),FSZ(I)),I=1,M)
C21   FORMAT(2X,'FXYZ=',9(E11.4,1X))
      DO 6 I=1,M
      AX(I)=FSX(I)/(AM(I)*1000.0)
      AY(I)=FSY(I)/(AM(I)*1000.0)
 6    AZ(I)=FSZ(I)/(AM(I)*1000.0)
      TN=TS+DT
      IT=IT+1
      DO 7 I=1,M
      VSX(I)=VX(I)
      VSY(I)=VY(I)
 7    VSZ(I)=VZ(I)
      DO 8 I=1,M
      VX(I)=VSX(I)+AX(I)*(TN-TS)
      VY(I)=VSY(I)+AY(I)*(TN-TS)
      VZ(I)=VSZ(I)+AZ(I)*(TN-TS)
      VS(I)=SQRT(VSX(I)*VSX(I)+VSY(I)*VSY(I)+VSZ(I)*VSZ(I))
      VN(I)=SQRT(VX(I)*VX(I)+VY(I)*VY(I)+VZ(I)*VZ(I))
      X(I)=X(I)+(VSX(I)+VX(I))*(TN-TS)*0.5
      Y(I)=Y(I)+(VSY(I)+VY(I))*(TN-TS)*0.5
      Z(I)=Z(I)+(VSZ(I)+VZ(I))*(TN-TS)*0.5
 8    RS(I)=SQRT(X(I)*X(I)+Y(I)*Y(I)+Z(I)*Z(I))
C     PRINT 20,(VN(I),I=1,M),ECULT,(RS(I),I=1,M)
C     PRINT 23,IT
c     IF(IT.GT.50) GO TO 11
      IF(IT.GT.50) DT=4.
      IF(IT.GT.75) DT=10.
      IF(IT.GT.100) GO TO 11
      TS=TN
      GO TO 10
 11   CONTINUE
C20   FORMAT(1X,' VN=',3(1X,E11.4),3X,'ECULT=',E11.4,
C    *' RS=',3(1X,E11.4))
C23   FORMAT(2X,'IT=',I5)
      EE=0.
      DO 30 I=1,M
 30   EE=EE+(AM(I)*0.5)*VN(I)*VN(I)
      ETA=SQRT((EC+ET)/EE)
      DO 31 I=1,M
      VX(I)=ETA*VX(I)
      VY(I)=ETA*VY(I)
      VZ(I)=ETA*VZ(I)
 31   CONTINUE
      RETURN
      END
      
      
      
      SUBROUTINE PLACE(RSYS,RN)
      COMMON /BLPAT/ PA(500),PAZ(500),IP
      COMMON /BLCEN/ XC(500),YC(500),ZC(500)
      COMMON /CELIPS/CEL
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
C FOR PROLONGED ELLIPSOID CEL>1
      CEL=1.
C FINDING FRAGMENT POSITIONS IN THE SYSTEM
 1    I=1
      R=(RSYS-RN*PA(I)**0.3333333)*(RRAN(ix))**0.3333333
      CT=1.-2.*RRAN(ix)
      PHI=6.28318*RRAN(ix)
      ST=SQRT(1-CT**2)
      CPHI=COS(PHI)
      SPHI=SIN(PHI)
      XC(1)=R*ST*CPHI
      YC(1)=R*ST*SPHI
C      ZC(1)=R*CT
      ZC(1)=R*CT*CEL
 3    K=I
      KK=0
      I=I+1
 4    CONTINUE
      KK=KK+1
      IF(KK.GT.1000) GO TO 2
      R=(RSYS-RN*PA(I)**0.3333333)*(RRAN(ix))**0.3333333
      CT=1.-2.*RRAN(ix)
      PHI=6.28318*RRAN(ix)
      ST=SQRT(1-CT**2)
      CPHI=COS(PHI)
      SPHI=SIN(PHI)
      XC(I)=R*ST*CPHI
      YC(I)=R*ST*SPHI
C      ZC(I)=R*CT
      ZC(I)=R*CT*CEL
      DO 5 J=1,K
      RR=(XC(I)-XC(J))**2+(YC(I)-YC(J))**2+(ZC(I)-ZC(J))**2
      RMIN=RN*PA(I)**0.3333333+RN*PA(J)**0.3333333
      IF(RR-RMIN*RMIN) 4,5,5
 5    CONTINUE
      IF(I.EQ.IP) RETURN
      GO TO 3
 2    CONTINUE
C2    PRINT 10,KK,I,PA(I),IP
C10   FORMAT(3X,'ERROR IN PLACE - KK=',I4,' I=',I3,' PA(I)=',F5.1,
C    *' IP=',I3)
      GO TO 1
      END
      
      
      
      SUBROUTINE RAZVAL(UP,AP,ZP,PN,KST)
C [EPM[ PA[BA[ [E[KO[O [[PA C A=AP [ Z=ZP. [HEP[[[ BO[[[[[EH[[-UP ([[B).
C HA[. [M[[[[C-PN(3) ([[B/C). KST-HOMEP [A[O[HEH[[ MACC[BA SPT.
C ******** HOBA[ BEPC[[ - [ACT[[[ [A[A[TC[ B TA[[[[E *********
C          ( [[[T[BA[TC[ [ACT[[[ [O A=16 )
      COMMON /LEVERS/ PARLEV(40)               ! DKFZg
      COMMON /ILEVRA/ILEVRA /FKAP/FKAP1,FKAP2  ! DKFZg
      COMMON /BLOK79/SFR(6,99)
      COMMON /SPDM/ MS(12,11),DM(12,11) /BLOKC/SPTU(10,500)
      COMMON /WWW/WK(20),WA(20),WN(100) /BLENN/ENN(100) /SUMN/NK(20),NF
      COMMON /BLJ/J(20),IEND /BLIAN/IAN(100) /BLIZN/IZN(100)
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      DIMENSION PN(3),VN(3),PNC(3),PNL(3),PN0(3,500),AMK(500)
      DIMENSION IPA(20),IPZ(20),WCON(2000),JCON(9000)
      KSTOLD=KST
C      ILEVRA=2   ! BotvinaDKFZ: ILEVRA=1 ground states only
      ILEVRA=INT(PARLEV(32))      ! DKFZg
C      FKAP=1.    ! BotvinaDKFZ: FKAP>1 more fragments
      FKAP1=PARLEV(33) ! Coulomb  ! DKFZg
      FKAP2=PARLEV(34) ! Volume   ! DKFZg
C
      IA0=INT(AP+0.5)
      IZ0=INT(ZP+0.5)
      NF=0
      DO 69 I=1,100
 69   WN(I)=0.
      DO 70 I=1,20
      NK(I)=0
      WA(I)=0.
 70   WK(I)=0.
      IF(AP.LT.1.99) GO TO 37
      PQ=PN(1)**2+PN(2)**2+PN(3)**2
      EAP=SQRT(PQ+(0.94*AP)**2)
      DO 40 I=1,3
 40   VN(I)=PN(I)/EAP
      CALL HELPRA
      DAZ=UP+GMASS(IA0,IZ0)
      SRQ=0.
      K2=0
C#### Instead of these instructions
C      DO 71 K=2,6
C      DO 72 I=1,K
C 72   J(I)=1
C 78   CONTINUE
C#### Write the following ones, BotvinaDKFZ
      K00=IA0                        ! DKFZg
      DO 71 K=2,K00
      J(1)=0
      DO 72 I=2,K
 72   J(I)=1
 78   CONTINUE
C##########################################      
      CALL DISNUM(IA0,IZ0,K)
C     IF(IEND.EQ.1) PRINT 7,(J(I),I=1,K)
C7    FORMAT(2X,'**** IEND=1 ****  J=',12I5)
      IF(IEND.EQ.1) GO TO 71
C     PRINT 3,(J(I),I=1,K)
C3    FORMAT(2X,'OTO[P.KOH[[[. J=',12I5)
      NF=NF+1
      NK(K)=NK(K)+1
      TQ=BCUL(K)
      RMQ=WEIGHT(IA0,K,DAZ,TQ)
      SRQ=SRQ+RMQ
      WCON(NF)=RMQ
      DO 73 I=1,K
      JI=J(I)
      JCON(K2+I)=JI
      IA=IAN(JI)
      IZ=IZN(JI)
      WN(JI)=WN(JI)+RMQ
 73   WA(IA)=WA(IA)+RMQ
      K2=K2+K
      WK(K)=WK(K)+RMQ
      GO TO 78
 71   CONTINUE
      IF(SRQ.LE.0.) GO TO 37
      DO 33 I=1,IA0
      WK(I)=WK(I)/SRQ
 33   WA(I)=WA(I)/SRQ
      DO 34 I=1,100
 34   WN(I)=WN(I)/SRQ
      DO 32 I=1,NF
 32   WCON(I)= WCON(I)/SRQ
              BR=RRAN(ix)
C     PRINT 106,K2,BR
C106  FORMAT(2X,'[[C[O BCEX [ACT[[ B MACC[BE K2=',I5,' C[[[. BR=',F8.5)
      SR=0.
      DO 13 N=1,NF
      NB=N
      SR=SR+WCON(N)
      IF(BR.LT.SR) GO TO 14
 13   CONTINUE
   14 CONTINUE
      NNP=0
      NN=0
C#### Instead of this instruction
C      DO 61 I=2,6
C#### Write the following one, BotvinaDKFZ
      DO 61 I=2,K00                 ! DKFZg
C####################
      NNP=NNP+NK(I)*I
      NN=NN+NK(I)
      IF(NN.GE.NB) GO TO 62
 61   CONTINUE
 62   K=I
      NPK1=(NB-(NN-NK(K))-1)*K
      NPK0=(NNP-NK(K)*K)+NPK1
C     PRINT 101,NN,NB,NPK1,NPK0,K
C101  FORMAT(2X,'NN=',I5,' NB=',I5,' NPK1=',I5,' NPK0=',I5,' K=',I5)
      TN=DAZ
      DO 63 I=1,K
      JI=JCON(NPK0+I)
C     PRINT 105,JI
C105  FORMAT(5X,'JI=',I5)
      IPA(I)=IAN(JI)
      IPZ(I)=IZN(JI)
      AMK(I)=GMASS(IPA(I),IPZ(I))+0.001*ENN(JI)
      TN=TN-AMK(I)
 63   CONTINUE
C     PRINT 102,(IPA(I),IPZ(I),I=1,K)
C102  FORMAT(2X,'IPA,IPZ=',6(2I5,5X))
C     PRINT 103,UP,DAZ,TN,(AMK(I),I=1,K)
C103  FORMAT(2X,'UP,DAZ,TN=',3F8.4,' AMK=',8F8.4)
      CALL DISIMP(K,AMK,PN0,TN)
      DO 35 I=1,K
      DO 36 II=1,3
 36   PNC(II)=PN0(II,I)
      EN=SQRT(PNC(1)**2+PNC(2)**2+PNC(3)**2+AMK(I)**2)
      CALL CLPV(PNC,VN,PNL,EN)
      CALL PINT(PNL,CT,ST,CF,SF,TK,AMK(I))
C     PRINT 104,AMK(I),TK,(PNL(II),II=1,3)
C104  FORMAT(2X,'AMK=',F8.4,' TK=',F8.4,' PNL=',3F8.4)
      SPTU(4,KST)=CT
      SPTU(5,KST)=SF
      SPTU(6,KST)=CF
      SPTU(7,KST)=TK
      SPTU(8,KST)=IPZ(I)
      SPTU(9,KST)=0.94*INT(AMK(I)/0.94+0.5)
      KST=KST+1
 35   CONTINUE
      GO TO 200
C      RETURN
 37   CONTINUE
      IF(AP.LT.0.99) RETURN
      WK(1)=1.
      WA(IA0)=1.
      NF=1
      NK(1)=1
      DO 38 II=1,100
      IA1=IAN(II)
      IZ1=IZN(II)
      IF(IA0.EQ.IA1.AND.IZ0.EQ.IZ1) WN(II)=WN(II)+1.
      IF(IA0.EQ.IA1.AND.IZ0.EQ.IZ1) GO TO 39
 38   CONTINUE
 39   CONTINUE
      CMN=0.94*AP
      CALL PINT(PN,CT,ST,CF,SF,TMN,CMN)
      SPTU(4,KST)=CT
      SPTU(5,KST)=SF
      SPTU(6,KST)=CF
      SPTU(7,KST)=TMN
      SPTU(8,KST)=ZP
      SPTU(9,KST)=CMN
      KST=KST+1
      GO TO 200
C      RETURN
 200  KST1=KST-1
C      PRINT 203,KSTOLD,KST1
C 203  FORMAT(5X,'KSTOLD=',I3,' KST1=',I3,' ***** SPTU (BEFORE) *****')
C      PRINT 204,((SPTU(I1,I2),I1=4,9),I2=KSTOLD,KST1)
C 204  FORMAT(2X,6F10.4)
      CALL POVRAZ(KSTOLD,KST1,KST2)
      KSS=KSTOLD-1
      DO 201 IPOV=1,KST2
      KSS=KSS+1
      DO 201 JPOV=1,6
      JPOVN=JPOV+3
      SPTU(JPOVN,KSS)=SFR(JPOV,IPOV)
 201  CONTINUE
      KST=KSS+1
C      PRINT 205,KST2
C 205  FORMAT(5X,' KST2=',I3,' **** SFR ****')
C      PRINT 204,((SFR(I1,I2),I1=1,6),I2=1,KST2)
C      PRINT 206,KSTOLD,KSS
C 206  FORMAT(5X,' KSTOLD=',I3,' KSS=',I3,' **** SPTU (AFTER) ****')
C      PRINT 204,((SPTU(I1,I2),I1=4,9),I2=KSTOLD,KSS)
      RETURN
      END
      
      
      
      FUNCTION WEIGHT(MA,K,DAZ,TQ)
C O[PE[E[EH[E [P[BE[EHHO[ BEPO[THOCT[ PAC[A[A MA HA K [PA[MEHTOB. HOMEP
C [PA[M. [AH[ B MACC[BE J(I). DAZ-[HEP[[[ MA,BK[[[A[ [HEP[[[ [OKO[. TQ-
C K[[OHOBCK[[ [AP[EP. AKV-[P[BE[EHH[[ O[[EM [[[ O[HO[O [PA[MEHTA. AKM-[P
C BE[EH[E C[[HOB[X [ MACCOB[X [AKTOPOB. RPM-[AKTOP[A[ TO[[ECTBEHH[X [[E-
C HOB. GAM-O[PATHA[ [AMMA-[[HK[[[. TN-K[HET[[ECKA[ [HEP[[[ OCKO[KOB.
      COMMON /SPDM/ MS(12,11),DM(12,11) /VAK/VAK /GAF/GAF(20)
      COMMON /BLJ/J(20),IEND /BLIAN/IAN(100) /BLIZN/IZN(100)
      COMMON /BLISP/ISP(100) /BLENN/ENN(100)
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      BNQ=DAZ
      RMQ=0.
      SPM=1.
      AKM=1.
      AKV=MA*VAK
      DO 40 I=1,K
      JI=J(I)
      IA=IAN(JI)
      IZ=IZN(JI)
      AKM=AKM*IA
      SPM=SPM*ISP(JI)
      BNQ=BNQ-GMASS(IA,IZ)-0.001*ENN(JI)
 40   CONTINUE
      TN=BNQ
      IF(TN.LE.TQ) GO TO 20
      TN=TN-TQ
      AKM=AKM/MA
      AKM=AKM*SQRT(AKM)*SPM
      IF(K.GT.2) GO TO 29
      RMQ=1.1283792*AKV*AKM*SQRT(TN)
      IF(J(1).EQ.J(2)) RMQ=0.5*RMQ
      GO TO 20
 29   CONTINUE
      KI=K-1
      RPM=1.
      VMK=1.
      TEQ=2.71828183*TN/(1.5*K-2.5)
      VTK=AKV*TEQ*SQRT(TEQ)
      DO 21 I=1,KI
      MRS=1
      IK=I+1
      VMK=VMK*VTK
      DO 22 I1=IK,K
      IF(J(I).EQ.J(I1)) MRS=MRS+1
 22   CONTINUE
      RPM=RPM*MRS
 21   CONTINUE
      GAM=GAF(K)
      RMQ=VMK*GAM*AKM/(TEQ*RPM)
 20   WEIGHT=RMQ
      RETURN
      END
      
      
      
      FUNCTION GMASS(MA,MZ)
C O[PE[E[EH[E MACC[ [E[K[X [[EP B [[B.([C[O[[[[[TC[ [[[EPO[H[E E[[H[[[)
C ( OCHOBHOE COCTO[H[E)
      COMMON /SPDM/ MS(12,11),DM(12,11)
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      I2=MZ
      I=MA-I2
      IF(I.LT.0.OR.I.GT.11) GO TO 10
      IF(I2.LT.0.OR.I2.GT.10) GO TO 10
      DMP=DM(I+1,I2+1)
      IF(DMP.GT.90) GO TO 10
      GMASS=0.931437*(I+I2)+0.001*DMP
      RETURN
 10   GMASS=0.939*MA
      RETURN
      END
      
      
      
      FUNCTION BCUL(K)
C  O[PE[E[EH[E K[[OHOBCKO[O [AP[EPA [[[ PA[BA[A HA KOH[[[[PA[[[ J(1-K)
C  [[ K [ACT[[.
      COMMON /FKAP/FKAP1,FKAP2   ! DKFZg
      COMMON /BLJ/J(20),IEND /BLIAN/IAN(100) /BLIZN/IZN(100)
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      DIMENSION IA(20),IZ(20)
      A=0.
      Z=0.
      COEF=(3.*1.44/(5.*1.3))*(1./(1.+FKAP1)**0.333333)   ! DKFZg
      DO 1 I=1,K
      JI=J(I)
      IA(I)=IAN(JI)
      A=A+IA(I)
      IZ(I)=IZN(JI)
      Z=Z+IZ(I)
 1    CONTINUE
      EC=0.
      DO 4 I=1,K
      AP=IA(I)
      ZP=IZ(I)
      EC=EC+ZP*ZP/(AP**0.333333)
 4    CONTINUE
      EC=Z*Z/(A**0.333333)-EC
      BCUL=0.001*COEF*EC
      IF(BCUL.LT.0.) BCUL=0.
      RETURN
      END
      
      
      
      SUBROUTINE HELPRA
C MS(I,J)-MACC[B C[[HOB[X [AKTOPOB [E[K[X [[EP =(2*S+1),[[E I=A-Z+1,
C J=Z+1. DM(I,J)-MACC[B [E[EKTOB MACC B [[[EPO[H[X E[[H[[AX. [[[T[BA[TC[
C [[PA C 0.LE.Z.LE.10 [ Z.LE.A.LE.(11+Z); MO[HO C[[TAT[ PA[BA[ [E[K[X
C [[EP C Z0.LE.10 [ (A0-Z0).LE.11.
C ([[[T[BA[TC[ C[[H[             OCHOBH[X      COCTO[H[[,   [EPETC[ [HEP
C [[[ OCHOBHO[O COCTO[H[[ - [[[ [[EP HA[[HA[ C [[T[[).
C KPOME TO[O,[POBO[[TC[ BC[OMO[ATE[[H[E B[[[C[EH[[ [[[ O[PE[E[EH[[ [AMMA
C [[HK[[[ ([O CT[P[[H[[) [ [P[BE[EHHO[O O[[EMA.
      COMMON /ILEVRA/ILEVRA /FKAP/FKAP1,FKAP2   ! DKFZg
      COMMON /VAK/VAK /GAF/GAF(20) /SPDM/MS(12,11),DM(12,11)
      COMMON /BLIAN1/IAN1(100) /BLIZN1/IZN1(100) /BLENN1/ENN1(100)
      COMMON /BLIAN2/IAN2(100) /BLIZN2/IZN2(100) /BLENN2/ENN2(100)
      COMMON /BLIAN/IAN(100) /BLIZN/IZN(100) /BLENN/ENN(100)
      COMMON /BLISP1/ISP1(100) /BLISP2/ISP2(100) /BLISP/ISP(100)
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      DIMENSION IPA(12),JPZ(12),DMP(12),MSP(12),AMP(12),MD(12,11),
     *MS1(12,11)
      DATA MS1/
     * 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     * 2, 3, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     * 0, 2, 1, 4, 1, 0, 1, 0, 0, 0, 0, 0,
     * 0, 0, 4, 3, 4, 5, 4, 0, 0, 0, 0, 0,
     * 0, 0, 1, 4, 1, 4, 1, 2, 1, 0, 0, 0,
     * 0, 0, 0, 5, 4, 7, 4, 3, 4, 0, 0, 0,
     * 0, 0, 0, 4, 1, 4, 1, 2, 1, 2, 1, 0,
     * 0, 0, 0, 0, 0, 3, 2, 3, 2, 5, 2, 0,
     * 0, 0, 0, 0, 0, 4, 1, 2, 1, 6, 1, 6,
     * 0, 0, 0, 0, 0, 0, 0, 1, 6, 3, 2, 5,
     * 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 1, 4/
      DATA MD/
     *99000,
     * 8071,99999,99999,99999,99999,99999,99999,99999,99999,99999,99999,
     * 7289,
     *13136,14950,25920,33790,99999,99999,99999,99999,99999,99999,99999,
     *99000,
     *14931, 2425,11390,17597,26110,31603,42030,50130,65000,75240,89260,
     *99000,
     *25130,11680,14087,14907,20947,24955,33250,40940,52940,61570,72280,
     *99000,
     *31700,18375,15769, 4942,11348,12607,20179,25072,35720,40720,51210,
     *99000,
     *99999,27940,22920,12416,12052, 8668,13370,16560,24230,29410,37960,
     *99000,
     *99999,35093,28912,15702,10650,00000, 3125, 3020, 9873,13693,17560,
     *99000,99999,99999,99999,25450,17338, 5346, 2864,  102, 5683, 7871,
     *13274,99000,99999,99999,42700,32070,23110, 8010, 2860,-4737, -810,
     * -783, 3332,99000,99999,99999,99999,39700,33380,17610,10693, 1952,
     *  873,-1486,  -16,99000,99999,99999,99999,49400,36400,24110,16470,
     * 5320, 1750,-7040,-5730/
      DO 22 I=1,12
      DO 22 I2=1,11
      MS(I,I2)=MS1(I,I2)
 22   DM(I,I2)=0.001*MD(I,I2)
      VL=1.3/(0.21*SQRT(0.94))
      VOL=FKAP2*VL*VL*VL                     ! DKFZg
      VAK=VOL*SQRT(2./3.14159265)*0.333333
      GAF(1)=0.
      GAF(2)=1./SQRT(3.1416)
      DO 10 K=3,20
      QK=1./(1.5*K-2.5)
      GQ=1.+QK*(1./12.+QK*(1./288.-QK*(139./51840.)))
 10   GAF(K)=SQRT(QK*0.1591549)/GQ
      IF(ILEVRA-1) 31,30,31
 30   DO 32 I=1,100
      IAN(I)=IAN1(I)
      IZN(I)=IZN1(I)
      ISP(I)=ISP1(I)
 32   ENN(I)=ENN1(I)
      RETURN
 31   DO 33 I=1,100
      IAN(I)=IAN2(I)
      IZN(I)=IZN2(I)
      ISP(I)=ISP2(I)
 33   ENN(I)=ENN2(I)
      RETURN
      END
      
      
      
      SUBROUTINE DISNUM(IA0,IZ0,K)
C  HAXO[[EH[E K [ACT[[ B KOH[[[[PA[[[ B B[[E [X HOMEPOB J(1-K) B MAC-
C  C[BAX IAN (B [OP[[.BO[POC.A) [ IZN (B [OP[[.BO[POC.Z [P[ [AHHOM A)
C  BO[MO[HA[ KOH[[[[PA[[[ IEND=0, KOHE[ IEND=1.
C  HA[.KOH[[[.J(1-K)=1 [A[AETC[ [EPE[ DISNUM.
      COMMON /BLIAN/IAN(100) /BLIZN/IZN(100) /BLJ/J(20),IEND
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
 5    L=0
 1    L=L+1
      IF(L.LT.K) GO TO 10
      IF(L.GE.K) GO TO 11
 10   J(L)=J(L)+1
      IF(J(L).GT.J(L+1)) GO TO 2
      ISA=0
      DO 12 I=1,K
      JI=J(I)
 12   ISA=ISA+IAN(JI)
      IF(ISA.GT.IA0) GO TO 2
      IF(ISA.LT.IA0) GO TO 5
 4    ISZ=0
      DO 13 I=1,K
      JI=J(I)
 13   ISZ=ISZ+IZN(JI)
      IF(ISZ.NE.IZ0) GO TO 5
      IEND=0
      RETURN
 2    J(L)=1
      GO TO 1
 11   J(L)=J(L)+1
      IF(J(L).GT.100) GO TO 15
      ISA=0
      DO 14 I=1,K
      JI=J(I)
 14   ISA=ISA+IAN(JI)
      IF(ISA.LT.IA0) GO TO 5
      IF(ISA.EQ.IA0) GO TO 4
 15   IEND=1
      RETURN
      END
      
      
      
      BLOCK DATA DMAS01
C  TA[[[[[,[A[A[[[E MACC[ (IAN),[AP[[[ (IZN),C[[H.MHO[.(ISP),[HEP[.OCH.
C  COCT.(ENN) [ACT[[ .
      COMMON /BLIAN1/IAN1(100) /BLIZN1/IZN1(100) /BLENN1/ENN1(100)
      COMMON /BLIAN2/IAN2(100) /BLIZN2/IZN2(100) /BLENN2/ENN2(100)
      COMMON /BLISP1/ISP1(100) /BLISP2/ISP2(100)
C ******   OCHOBH.+ [CTO[[[B[E BO[[[[[. COCTO[H[[  ******
      DATA IAN2/
     * 1, 1, 2, 3, 3, 4, 5, 5, 6, 6,   6, 7, 7, 7, 7, 8, 8, 8, 9, 9,
     *10,10,10,10,10,10,10,10,10,10,  10,10,
     *                                11,11,11,11,11,11,11,11,11,11,
     *11,11,11,11,11,11,11,11,11,12,  12,12,12,12,12,13,13,13,13,13,
     *14,14,14,14,14,14,14,14,14,14,  14,14,14,14,14,15,15,15,15,15,
     *15,15,15,15,15,15,15,15,15,15,  16,16,16,16,16,16,16,16/
      DATA IZN2/
     * 0, 1, 1, 1, 2, 2, 2, 3, 2, 3,   3, 3, 3, 4, 4, 3, 3, 4, 4, 5,
     * 4, 4, 4, 4, 4, 5, 5, 5, 5, 5,   6, 6,
     *                                 5, 5, 5, 5, 5, 5, 5, 5, 6, 6,
     * 6, 6, 6, 6, 6, 6, 6, 6, 6, 5,   5, 5, 5, 6, 6, 6, 6, 6, 6, 7,
     * 6, 6, 6, 6, 6, 7, 7, 7, 7, 7,   7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
     * 7, 7, 7, 7, 7, 8, 8, 8, 8, 8,   7, 7, 7, 7, 8, 8, 8, 8/
      DATA ISP2/
     * 2, 2, 3, 2, 2, 1, 4, 4, 1, 3,   1, 4, 2, 4, 2, 5, 3, 1, 4, 4,
     * 1, 5, 8, 1, 5, 7, 3, 1, 3, 5,   3, 5,
     *                                 4, 2, 6, 4,10, 6, 4, 6, 4, 2,
     * 6, 4, 2, 8, 6, 4, 4, 6, 8, 3,   5, 5, 4, 1, 5, 2, 2, 4, 6, 2,
     * 1, 3, 8, 6, 5, 3, 1, 3, 1, 5,   3, 7, 3, 7, 5, 2, 8, 4,10, 8,
     * 2, 4,14,14, 8, 2, 8, 4,10, 8,   5, 1, 7, 3, 1, 8, 5, 3/
      DATA ENN2/
     * 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,16.76,16.66, 0.00, 0.00,
     * 3.56, 0.00, 0.48, 0.00, 0.43, 0.00, 0.98, 0.00, 0.00, 0.00,
     * 0.00, 3.37, 5.96, 6.18, 6.26, 0.00, 0.72, 1.74, 2.15, 3.59,
     * 0.00, 3.35,
     * 0.00, 2.13, 4.44, 5.02, 6.76, 7.29, 7.98, 8.56, 0.00, 2.00,
     * 4.32, 4.80, 6.34, 6.48, 6.90, 7.50, 8.10, 8.42, 8.66, 0.00,
     * 0.95, 1.67, 2.65, 0.00, 4.44, 0.00, 3.09, 3.68, 3.85, 0.00,
     * 0.00, 6.09, 6.69, 6.96, 7.34, 0.00, 2.31, 3.95, 4.92, 5.11,
     * 5.69, 5.83, 6.20, 6.44, 7.03, 0.00, 5.28, 6.32, 7.22, 7.57,
     * 8.31, 8.57, 9.15, 9.79,10.00, 0.00, 5.22, 6.18, 6.83, 7.28,
     * 0.00, 0.12, 0.30, 0.40, 0.00, 6.10, 6.92, 7.12/
C  ******** OCHOBH. COCTO[H[[ *******
      DATA IAN1/
     * 1, 1, 2, 3, 3, 4, 5, 5, 6, 6,   7, 7, 8, 8, 9, 9,10,10,10,11,11,
     *12,12,13,13,14,14,15,15,16,16,  17,18,67*0/
      DATA IZN1/
     * 0, 1, 1, 1, 2, 2, 2, 3, 2, 3,   3, 4, 3, 4, 4, 5, 4, 5, 6, 5, 6,
     * 5, 6, 6, 7, 6, 7, 7, 8, 7, 8,   8, 8,67*0/
      DATA ISP1/
     * 2, 2, 3, 2, 2, 1, 4, 4, 1, 3,   4, 4, 5, 1, 4, 4, 1, 7, 3, 4, 4,
     * 3, 1, 2, 2, 1, 3, 2, 2, 5, 1,   6, 1,67*0/
      DATA ENN1/6*0.,16.76,16.66,92*0./
      END
      
      
      
      SUBROUTINE PINT(P,CT,ST,CF,SF,T,CM)
C O[PE[E[EH[E [O MACCE CM [ [M[[[[C[ P(3) [ACT[[[,EE K[H. [HEP[[[ T
C [ [[[OB - CT,ST,CF,SF ,O[PE[E[[[[[X HA[PAB[EH[E [O[ETA.
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      DIMENSION P(3)
      PZ=P(3)**2
      PQ=P(1)**2+P(2)**2+PZ
      IF(PQ.LE.0.) GO TO 11
      CTQ=PZ/PQ
      T=SQRT(PQ+CM**2)-CM
      IF(CTQ.GE.1.) GO TO 10
      CT=SQRT(CTQ)
      IF(P(3).LE.0.) CT=-CT
      ST=SQRT(1.-CTQ)
      PMH=ST*SQRT(PQ)
      CF=P(1)/PMH
      SF=P(2)/PMH
      RETURN
 10   ST=0.
      CT=1.
      SF=0.
      CF=1.
      RETURN
 11   ST=0.
      CT=1.
      SF=0.
      CF=1.
      T=0.
      RETURN
      ENTRY TINP(P,CT,ST,CF,SF,T,CM)
C O[PE[E[EH[E [O MACCE CM,K[H.[HEP[[[ T [ [[[AM [ACT[[[,EE [M[[[[CA P(3)
      PM=SQRT(T*(T+2.*CM))
      P(3)=PM*CT
      PM=PM*ST
      P(1)=PM*CF
      P(2)=PM*SF
      RETURN
      END
      
      
      
      SUBROUTINE CLPV(P1,V,P2,E1)
C O[PE[E[EH[E [M[[[[CA [ACT[[[ P2(3) B [.C.K ,C[[TA[,[TO CTAPA[ C.K.([DE
C [M[[[[C [ACT[[[ P1(3)) [B[[ETC[ CO CKOP.V(3).E1-[O[H.[HEP[.[ACT[[[ B C
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      DIMENSION P1(3),V(3),P2(3)
      SPV=0.
      V2=0.
      DO 10 I=1,3
      V2=V2+V(I)**2
      SPV=SPV+P1(I)*V(I)
 10   CONTINUE
      GAM=1./SQRT(1.-V2)
      TEMP=GAM*(SPV*GAM/(GAM+1.)+E1)
      DO 20 K=1,3
 20   P2(K)=P1(K)+V(K)*TEMP
      RETURN
      END
      
      
      
      SUBROUTINE ISOTR(ANL)
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      DIMENSION ANL(3)
      CT=1.-2.*RRAN(ix)
      ST=SQRT(ABS(1.-CT*CT))
      ANL(3)=CT
      FI=6.2831843*RRAN(ix)
      ANL(1)=COS(FI)*ST
      ANL(2)=SIN(FI)*ST
      RETURN
      END
      
      
      
      SUBROUTINE DISIMP(K,AMK,PN0,TN)
C O[PE[E[EH[E [M[[[[COB K [PA[MEHTOB- PN0(3,1-K) ([[B/C) B C[CTEME [X
C [.[.  AMK(1-K) ([[B) -[X MACC[. TN ([[B) - [X [O[HA[ K[H. [HEP[[[.
C (AHA[O[[[HO [.[.  KO[[[OB[-1970, C.229-239)
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      DIMENSION ANL(3),PNC(3),VRS(3),PN0(3,500),AMK(500)
      SMK=0.
      DO 18 I=1,K
 18   SMK=SMK+AMK(I)
      DO 11 I=1,3
 11   VRS(I)=0.
      TKN=TN
      KL=K-1
      DO 10 L=1,KL
      LK=K-L+1
      AMP=AMK(LK)
      SMK=SMK-AMP
      AMR=SMK
      TPR=TKN
      IF(LK.LT.3) GO TO 20
      TKM=TKN*RNKSI(LK-1)
      TPR=TKN-TKM
 20   PMC=SQRT(2.*((AMP*AMR)/(AMP+AMR))*TPR)
      CALL ISOTR(ANL)
      DO 12 I=1,3
      PNC(I)=PMC*ANL(I)
      PN0(I,LK)=VRS(I)*AMP+PNC(I)
      VRS(I)=VRS(I)-PNC(I)/AMR
 12   CONTINUE
      TKN=TKM
 10   CONTINUE
      DO 21 I=1,3
 21   PN0(I,1)=VRS(I)*AMR
      RETURN
      END
      
      
      
      FUNCTION RNKSI(K)
C HAXO[[EH[E RNKSI(K) O[PE[E[[[[EE K[H.[HEP[[[ K [ACT[[ -T(K) B C[CTEME
C [X [.[. [EPE[ K[H.[HEP[[[ K+1 [ACT[[ -T(K+1) B C[CTEMME [X [.[.
C           (   T(K)=RNKSI(K)*T(K+1)    )
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      CSIM=(3.*K-5.)/(3.*K-4.)
      PEX=1.5*K-2.5
      FCSIM=SQRT(1.-CSIM)*CSIM**PEX
 20   CSI=RRAN(ix)
      FCSI =SQRT(1.-CSI)*CSI**PEX
      RF=FCSIM*RRAN(ix)
      IF(RF.GT.FCSI) GO TO 20
      RNKSI=CSI
      RETURN
      END
      
      
      
      SUBROUTINE POVRAZ(KSTOLD,KST1,KST2)
C  [OBTOPH[[ PA[BA[ HECTA[[[[H[X [PA[MEHTOB. BXO[:KST1 [ACT[[ B MACC[BE
C  SPTU(10,500), B[XO[: KST2 [ACT[[ B MACC[BE SFR(6,99)
      COMMON /BLOKC/SPTU(10,500) /BLOK79/SFR(6,99)
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      DIMENSION P9(3,500),A9(500),PC9(3),PL9(3),V9(3),Z9(500)
      DO 1 I=1,6
      DO 1 J=1,99
 1    SFR(I,J)=0.
      KST2=0
      DO 2 I=KSTOLD,KST1
      CMN=SPTU(9,I)
      ZP=SPTU(8,I)
      TMN=SPTU(7,I)
      CF=SPTU(6,I)
      SF=SPTU(5,I)
      CT=SPTU(4,I)
      IF(CT.LT.1.) ST=SQRT(1.-CT*CT)
      IF(CT.GE.1.) ST=0.
      IA=INT(CMN/0.94+0.5)
      IZ=ZP
      GAM=1.+TMN/CMN
      V=SQRT((GAM**2-1.)/GAM**2)
      V9(3)=V*CT
      V9(2)=V*ST*SF
      V9(1)=V*ST*CF
      IF(IA.EQ.5.AND.IZ.EQ.2) GO TO 21
      IF(IA.EQ.5.AND.IZ.EQ.3) GO TO 22
      IF(IA.EQ.8.AND.IZ.EQ.4) GO TO 23
      IF(IA.EQ.9.AND.IZ.EQ.5) GO TO 24
      KST2=KST2+1
      DO 3 J=1,6
      J1=J+3
 3    SFR(J,KST2)=SPTU(J1,I)
      GO TO 2
 21   K9=2
      A9(1)=3.728173
      Z9(1)=2.
      A9(2)=0.939508
      Z9(2)=0.
      T9=0.001*(16.76+11.39-8.071-2.425)
      GO TO 100
 22   K9=2
      A9(1)=3.728173
      Z9(1)=2.
      A9(2)=0.938726
      Z9(2)=1.
      T9=0.001*(16.66+11.68-7.289-2.425)
      GO TO 100
 23   K9=2
      A9(1)=3.728173
      Z9(1)=2.
      A9(2)=3.728173
      Z9(2)=2.
      T9=0.001*(4.942-2*2.425)
      GO TO 100
 24   K9=3
      A9(1)=3.728173
      Z9(1)=2.
      A9(2)=3.728173
      Z9(2)=2.
      A9(3)=0.938726
      Z9(3)=1.
      T9=0.001*(12.416-7.289-2*2.425)
 100  CONTINUE
      CALL DISIMP(K9,A9,P9,T9)
      DO 101 J=1,K9
      DO 102 II=1,3
 102  PC9(II)=P9(II,J)
      E9=SQRT(PC9(1)**2+PC9(2)**2+PC9(3)**2+A9(J)**2)
      CALL CLPV(PC9,V9,PL9,E9)
      CALL PINT(PL9,C1,S1,C2,S2,T,A9(J))
      KST2=KST2+1
      SFR(1,KST2)=C1
      SFR(2,KST2)=S2
      SFR(3,KST2)=C2
      SFR(4,KST2)=T
      SFR(5,KST2)=Z9(J)
      SFR(6,KST2)=A9(J)
      IF(A9(J).GT.0.5) SFR(6,KST2)=0.94*INT(A9(J)/0.94+0.5)
 101  CONTINUE
 2    CONTINUE
C     PRINT 11,KST1,KST2
C11   FORMAT(25X,'[OC[E POVRAZ               KST1=',I3,' KST2=',I3)
C     DO 12 I=1,KST2
C     PRINT 10,SFR(6,I),SFR(5,I),SFR(4,I),SFR(1,I),SFR(3,I),SFR(2,I)
C10   FORMAT(2X,' A=',F7.3,'  Z=',F7.3,'  T=',F9.5,'  C1=',F9.5,
C    *'       C2,S2=',2F9.5)
C12   CONTINUE
      RETURN
      END
       
       
       
       SUBROUTINE EVANUC(ENEXT,ATWGHT,CHARGE,PNX,PNY,PNZ,KSTART)
       COMMON /BL0999/RNCL/BL1000/AM,AMF,AMCB(35)/AMNUCL/AMNEV(3)
       COMMON /BL1001/T1Y(130)/BL1002/T2XY(200)/BL1005/AJ(35)
     */BL1006/ZJ(35)/BL1014/GAM(35)/BL1011/VJ(35)/BL1015/RJ(35)
      COMMON /BL1003/U,A,Z/BL1009/AFJ(35)/BL1010/ZFJ(35)
      COMMON/BLOKC/SPTU(10,500)/BLANGL/ANGL(4)/BL1100/AMS,AMFS
      COMMON/BLFIS/FIS,EF,SNT(6,6)/BLCAN/ICAN
      COMMON /EXILEV/ ENLEV(20,35),SPLEV(20,35)
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
       DIMENSION GJ(35),BJ(35),GG(35),VN(3),PN(3),PNL(3),PP(3),PPL(3)
      DIMENSION DSH2(35)
C [C[APEH[E [E[K[X [PA[MEHTOB [[ T[[.[[PA [O BA[CKO[[[. [[[T[BAETC[
C COXPAHEH[E [HEP[[[,[M[[[[CA [ MOMEHTA [M[[[[CA(K[ACC[[.). [E[EH[E
C [[[T[BAETC[ KAK BO[MO[H[[ KAHA[ [PEKPA[A[[[[ [C[AP[TE[[H[[ KACKA[.
C PA[MEPHOCT[: [HEP[[[-[[B,[M[[[[C-[[B/C,CKOPOCT[-1/C,MACCA-AT.HOMEP.
C       ******************************
C  ICAN=0 -[[[T[BAETC[ 6 KAHA[OB [ [E[EH[E; ICAN=1 - 32 KAHA[A [ [E[EH[E
C  ICAN=2 - 6 KAHA[OB, [E[EH[E, 34 [ 35 KAHA[; ICAN=3 ([[[ [P.) -BCE KAH
C       ******************************
      ILEVEL=3
      CALL DATEVP(ILEVEL)
      FIS=0.
      SNT(1,1)=ENEXT
      SNT(2,1)=ATWGHT
      SNT(3,1)=CHARGE
      SNT(4,1)=PNX
      SNT(5,1)=PNY
      SNT(6,1)=PNZ
      U=ENEXT*1000.
      A=ATWGHT
      Z=CHARGE
      AMNEV(1)=0.
      AMNEV(2)=0.
      AMNEV(3)=0.
      REMN=940.*A
      PNL(1)=PNX*1000.
      PNL(2)=PNY*1000.
      PNL(3)=PNZ*1000.
      ENL=SQRT(PNL(1)**2+PNL(2)**2+PNL(3)**2+REMN**2)
      VN(1)=PNL(1)/ENL
      VN(2)=PNL(2)/ENL
      VN(3)=PNL(3)/ENL
      DO 20 K=KSTART,500
      IN=35
      DO 4 I=1,IN
      VJ(I)=0.
      AMCB(I)=0.
C      DSH2(I)=0.
 4    RJ(I)=-1000.
      CALL DELAM(A,Z,DL1,DSHEL1,BARIER)
      AM=AMS
C      IF(U.GT.0.) AM=AMS*(1.+(1.-EXP(-0.06*U))*DSHEL1/U)
      IF(Z.GE.89.) AMF=1.04*AMS
      IF(Z.GT.85.AND.Z.LT.89) AMF=AMS*(1.04+0.01*(89.-Z))
      IF(Z.LE.85.) AMF=1.08*AMS
      DO 6 I=1,IN
      IF(I.EQ.33) GO TO 6
      AFJ(I)=A-AJ(I)
      ZFJ(I)=Z-ZJ(I)
      IF(ICAN.EQ.0.AND.I.GT.6) GO TO 6
      IF(ICAN.EQ.1.AND.I.GT.32) GO TO 6
      IF(ICAN.EQ.2.AND.I.GT.6.AND.I.LE.32) GO TO 6
      IF(AFJ(I).LT.ZFJ(I)) GO TO 6
      IF(AFJ(I).LT.AJ(I).OR.ZFJ(I).LT.ZJ(I)) GO TO 6
      JAF=INT(AFJ(I)+0.5)
      JZF=INT(ZFJ(I)+0.5)
      JNF=JAF-JZF
      ODD=11.*(2+2*(JNF/2)-JNF+2*(JZF/2)-JZF)/SQRT(AFJ(I))
      IF(RNCL.GE.0.1) RADNCL=RNCL
      IF(RNCL.LT.0.1)
     *RADNCL=2.173*(1+0.006103*ZJ(I)*ZFJ(I))/(1+0.009443*ZJ(I)*ZFJ(I))
      VJ(I)=COLOMB(I,RADNCL)
      CALL DELAM(AFJ(I),ZFJ(I),DL2,DSHEL2,BAR2)
C      DSH2(I)=DSHEL2
      CALL DELAM(AJ(I),ZJ(I),DL3,DSHEL3,BAR3)
      BJ(I)=DL2+DL3-DL1
c ******** in accordance with DELA8
      IF(A.LE.55.) RJ(I)=U-(BJ(I)+VJ(I))-ODD
      TCOEF=1.-(A-55.)/10.
      IF(A.GT.55.AND.A.LT.65.) RJ(I)=U-(BJ(I)+VJ(I))-ODD*TCOEF
      IF(A.GE.65.) RJ(I)=U-(BJ(I)+VJ(I))
c ********
C     RJ(I)=U-(BJ(I)+VJ(I))-ODD
C     RJ(I)=U-(BJ(I)+VJ(I))
C      IF(RJ(I).GT.0.)
C     *AMCB(I)=AMS*(1.+(1.-EXP(-0.06*RJ(I)))*DSHEL2/RJ(I))
 6    CONTINUE
      EF=BARIER
C      AMCB(33)=AMF
      IF(A.GT.100.) RJ(33)=U-BARIER
       DO 7 I=1,IN
      IF(RNCL.GE.0.1) RADNCL=RNCL
      IF(RNCL.LT.0.1)
     *RADNCL=2.173*(1+0.006103*ZJ(I)*ZFJ(I))/(1+0.009443*ZJ(I)*ZFJ(I))
      AMC=AM
C      AMC=AMCB(I)
      GJ(I)=GAMMAC(I,A,U,AM,AMC,AMF,RADNCL)
      RJI=RJ(I)
      GAMI=GAM(I)
      DO 71 IL=1,20
      IF(SPLEV(IL,I).LT.0.1) GO TO 70
      GAM(I)=SPLEV(IL,I)*AJ(I)
      RJ(I)=RJI-ENLEV(IL,I)
C      IF(RJ(I).GT.0.)
C     *AMC=AMS*(1.+(1.-EXP(-0.06*RJ(I)))*DSH2(I)/RJ(I))
      GJ(I)=GJ(I)+GAMMAC(I,A,U,AM,AMC,AMF,RADNCL)
 71   CONTINUE
 70   RJ(I)=RJI
      GAM(I)=GAMI
7     CONTINUE
      G=0.
      DO 10 I=1,IN
10    G=G+GJ(I)
      IF(G) 11,11,12
 11   PNX=PNL(1)*0.001
      PNY=PNL(2)*0.001
      PNZ=PNL(3)*0.001
      ENEXT=U*0.001
      ATWGHT=A
      CHARGE=Z
      REMN=940.*A
      SPTU(9,K)=REMN*0.001
      SPTU(8,K)=Z
      CALL PINT(PNL,CT,ST,CF,SF,T,REMN)
      SPTU(4,K)=CT
      SPTU(5,K)=SF
      SPTU(6,K)=CF
      SPTU(7,K)=T*0.001
      SPTU(1,K)=1.
      SPTU(3,K)=FIS
C-(Sobolevsky)-res. excitation in MeV
      SPTU(10,K)=U
C------------------------------------
      KSTART=K+1
      RETURN
 12   CONTINUE
      DO 13 J=2,IN
13    GJ(J)=GJ(J-1)+GJ(J)
      BB=RRAN(ix)
      B=BB*G
      DO 14 J=1,IN
      IF(B-GJ(J)) 15,14,14
15    LM=J
      GO TO 16
 14   CONTINUE
 16     IF(LM-33) 18,17,18
 17   FIS=1.
      KSTART=K
      CALL FISION(U,A,Z,PNL(1),PNL(2),PNL(3),KSTART)
      RETURN
 18   AMC=AM
C18    AMC=AMCB(LM)
      EP1=TKIN(LM,AMC,AMF)
      EP2=ZJ(LM)
      EP3=940.*AJ(LM)
C     ALMAX=0.219327*1.2*(AFJ(LM)**0.333333+AJ(LM)**0.333333)
C    +*SQRT((AJ(LM)*AFJ(LM)/(AJ(LM)+AFJ(LM)))*(EP1-VJ(LM)))
C     BR=RRAN(ix)
C     AL=ALMAX*SQRT(BR)
C     BR1=1.-2.*RRAN(ix)
C     BR2=6.283185*RRAN(ix)
C     AMNEV(1)=AMNEV(1)-AL*SQRT(1.-BR1**2)*COS(BR2)
C     AMNEV(2)=AMNEV(2)-AL*SQRT(1.-BR1**2)*SIN(BR2)
C     AMNEV(3)=AMNEV(3)-AL*BR1
      U=U-BJ(LM)-EP1
      A=AFJ(LM)
      Z=ZFJ(LM)
C  VPM-OTHOC[TE[[HA[ CKOPOCT[ [PA[MEHTA
      VPM=SQRT((2.*EP1)/(EP3*AFJ(LM)/(AFJ(LM)+AJ(LM))))
      CALL ISANGL
C  [M[[[[C [PA[MEHTA B C.[.[.
      PP(1)=VPM*ANGL(4)*ANGL(3)*EP3/(1.+AJ(LM)/A)
      PP(2)=VPM*ANGL(4)*ANGL(2)*EP3/(1.+AJ(LM)/A)
      PP(3)=VPM*ANGL(1)*EP3/(1.+AJ(LM)/A)
      PN(1)=-PP(1)
      PN(2)=-PP(2)
      PN(3)=-PP(3)
      EP=SQRT(PP(1)**2+PP(2)**2+PP(3)**2+EP3**2)
      EN=SQRT(PN(1)**2+PN(2)**2+PN(3)**2+(940.*A)**2)
      CALL CLPV(PP,VN,PPL,EP)
      CALL CLPV(PN,VN,PNL,EN)
      ENL=SQRT(PNL(1)**2+PNL(2)**2+PNL(3)**2+(940.*A)**2)
      VN(1)=PNL(1)/ENL
      VN(2)=PNL(2)/ENL
      VN(3)=PNL(3)/ENL
      CALL PINT(PPL,CT,ST,CF,SF,T,EP3)
      SPTU(4,K)=CT
      SPTU(5,K)=SF
      SPTU(6,K)=CF
      SPTU(7,K)=T*0.001
      SPTU(8,K)=EP2
      SPTU(9,K)=EP3*0.001
      SPTU(2,K)=LM
20    CONTINUE
      PRINT 21,U,A,Z
21    FORMAT(35X,37HMASSIV SPT EXCEEDED AFTER EVAPORATION/40X,2HU=,
     *F10.5,4H  A=,F5.1,4H  Z=,F4.1)
       RETURN
       END
      
      
      
      SUBROUTINE DATEVP(I)
C I=1 (I<2) -[[ET OCHOBH[X COCTO[H[[ [C[AP[[[E[C[ [ACT[[[,I=2-[[ET
C [EPB[X BO[[[[[EHH[X COCTO[H[[ ([O [POBHE[,KOTOP[E [C[[CKA[T H[K[OH[),
C I=3 (I>2) -[CPE[HEHH[[ [[ET H[[H[X [POBHE[ [P[ B[COK[X [HEP[[[X BO[[.
      COMMON /BLGAM/ GAM1(35),GAM2(35) /BL1014/ GAM(35)
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      IF(I-2) 1,1,2
 1    DO 3 J=1,35
 3    GAM(J)=GAM1(J)
      IF(I.EQ.2) CALL FILLEV
      GO TO 5
 2    DO 4 J=1,35
 4    GAM(J)=GAM2(J)
 5    RETURN
      END
       
       
       
       FUNCTION COLOMB(L,RADNCL)
       COMMON /BL1006/ZJ(35)/BL1010/ZFJ(35)/BL1009/AFJ(35)
     */BL1005/AJ(35)/BL1003/U,A,Z
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
C B[[[C[EH[E K[[OHOBCKO[O [AP[EPA [C[[CKAH[[ [PA[MEHTA COPTA L (M[B)
       IF(L-1)1,1,2
1      COLOMB=0.
       RETURN
2      TEMP1=1.44/RADNCL
       COLOMB=TEMP1*((ZJ(L)*ZFJ(L))/(AJ(L)**.33333+AFJ(L)**.33333))
       IF(COLOMB) 3,3,4
3       COLOMB=0.
4       CONTINUE
       RETURN
       END
      
      
      
      SUBROUTINE ISANGL
C [[OTPO[H[[ PO[[[P[[ [[[OB B[[ETA [PA[MEHTA B C[CTEME OCTAT.[[PA.
      COMMON /BLANGL/ ANGL(4)
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      ANGL(1)=1.-2.*RRAN(ix)
      ANGL(4)=SQRT(1.-ANGL(1)**2)
      F=2.*3.1415926536*RRAN(ix)
      ANGL(2)=SIN(F)
      ANGL(3)=COS(F)
      RETURN
      END
      
      
      
      FUNCTION GAMMAC(J,A,U,AM,AMC,AMF,RADNCL)
      COMMON /BLSEC/ISEC
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
C ISEC=2 - COOTBETCTB[ET [[ET[ [O[[AP[EPHO[O [C[[CKAH[[ [PA[MEHTOB C A.G
      IF(ISEC.NE.2.OR.J.LE.6) GAMMAC=
     *GAMMA1(J,A,U,AM,AMC,AMF,RADNCL)
      IF(ISEC.EQ.2.AND.J.GT.6) GAMMAC=
     *GAMMA2(J,A,U,AM,AMC,AMF,RADNCL)
      RETURN
      END
      
      
      
      FUNCTION TKIN(L,AMC,AMF)
      COMMON /BLSEC/ISEC
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
C ISEC=2 -COOTBETCTB[ET [[ET[ [O[[AP[EPHO[O [C[[CKAH[[ [PA[MEHTOB C A.GT
      IF(ISEC.NE.2.OR.L.LE.6) TKIN=TKIN1(L,AMC,AMF)
      IF(ISEC.EQ.2.AND.L.GT.6) TKIN=TKIN2(L,AMC,AMF)
      RETURN
      END
      
      
      
      FUNCTION GAMMA1(J,A,U,AM,AMC,AMF,RADNCL)
C  B[[[C[EH[E A[CO[[THO[ [[P[H[ [C[[CKAH[[ [PA[MEHTA COPTA J
C  ( B E[. (1./[P[BE[.[OCT.[[AHKA) M[B )
C  [[P[HA [E[EH[[ B TEX [E E[[H[[AX
       COMMON /BL1009/AFJ(35)/BL1015/RJ(35)/BL1014/GAM(35)
      COMMON /BL1011/VJ(35) /BL1005/AJ(35)/BL1006/ZJ(35)
      COMMON /BL1010/ZFJ(35)
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      IF(RJ(J).LE.0.OR.U.LE.0.) GO TO 12
      PER=2.*SQRT(AM*A*U)
      IF(J.EQ.33) GO TO 10
      RN=1.5
      CC=0.2
      IF(J.GT.2.AND.J.LE.6)  CC=0.1
      IF(J.GT.6) CC=    (AJ(J)/AFJ(J))**0.6666667
       IF(J-1)1,2,1
2      ALFA=.76+2.2/AFJ(1)**.33333
       BETA=(2.12/AFJ(1)**.66667-.05)/ALFA
      GO TO 3
1      ALFA=1.+CC
       BETA=0.
       GO TO 3
3      Q1=AMC*AFJ(J)
       Q2=AMC*AFJ(J)*RJ(J)
       Q3=(GAM(J)*AFJ(J)**.66667)*(ALFA/Q1**2)*(AFJ(J)/(AFJ(J)+AJ(J)))
      Q3=Q3*(3.1416*RN**2)/(2.*41.5*3.1416**2)
       Q4=(2.*BETA*Q1-3.)/2.+Q2
       Q5=(2.*BETA*Q1-3.)*(SQRT(Q2)-.5)+2.*Q2
      IF(PER-160.) 20,20,21
 20   PEX1=Q4*EXP(-PER)
      GO TO 22
 21   PEX1=0
 22   PP2=PER-2.*SQRT(Q2)
      IF(PP2-160.) 23,23,24
 23   PEX2=Q5*EXP(-PP2)
      GO TO 25
 24   PEX2=0.
 25   GAMMA1=Q3*(PEX1+PEX2)
      RETURN
 10   Q1=2.*SQRT(AMF*A*RJ(J))
      Q2=1./(4.*3.1416)
      IF(PER-160.) 30,30,31
 30   PEX1=EXP(-PER)
      GO TO 32
 31   PEX1=0.
 32   PP2=PER-Q1
      IF(PP2-160.) 33,33,34
 33   PEX2=(Q1-1.)*EXP(-PP2)
      GO TO 35
 34   PEX2=0.
 35   GAMMA1=(Q2/(AMF*A))*(PEX1+PEX2)
      TFIS=10000.
      GFIS=0.21*940./TFIS
      IF(GAMMA1.GT.GFIS) GAMMA1=GFIS
      RETURN
 12   GAMMA1=0.
      RETURN
       END
       
       
       
       FUNCTION TKIN1(L,AMC,AMF)
C  PO[[[P[[ K[H.[HEP[[[ B[[ETEB[E[O [PA[MEHTA COPTA L B C.[.[. (M[B)
       COMMON /BL1009/AFJ(35)/BL1015/RJ(35)/BL1011/VJ(35)
       COMMON /BL1005/AJ(35)
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
       RJL=RJ(L)
       RB=4.*AMC*AFJ(L)*RJL
       PP1=SQRT(RB)
 5     B1=RRAN(ix)
       IF(B1.LE.0.OR.B1.GT.1) B1=RRAN(ix)
       IF(PP1-160.) 21,21,22
 21    PEX1=EXP(-PP1)
       GO TO 23
 22    PEX1=0.
 23    RK=1.+(1./PP1)*ALOG(B1+(1.-B1)*PEX1)
       IF(L-1) 1,2,1
2      BETA=(2.12/AFJ(1)**0.66667-0.05)/(0.76+2.2/AFJ(1)**0.33333)
       Q1=1.+BETA/RJL
       Q2=Q1*SQRT(Q1)
       FRK=(((3.*SQRT(3.))/2.)/Q2)*(Q1*RK-RK**3)
       GO TO 3
1     FRK=((3.*SQRT(3.))/2.)*(RK-RK**3)
       GO TO 3
 3     B2=RRAN(ix)
       IF(B2-FRK) 4,4,5
4      TKIN1=   RJL*(1.-RK**2)+VJ(L)
       RETURN
      END
       
       
       
       FUNCTION TKIN2(L,AMC,AMF)
C  PO[[[P[[ K[H.[HEP[[[ B[[ETEB[E[O [PA[MEHTA COPTA L B C.[.[. (M[B)
       COMMON /BL1009/AFJ(35)/BL1015/RJ(35)/BL1011/VJ(35)
       COMMON /BL1005/AJ(35) /BL1010/ZFJ(35)
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
C  [[[T[BAETC[ BO[MO[HOCT[ [O[[AP[EPHO[O [C[[CKAH[[ [PA[MEHTOB.
      DCOL=1.
      DALF=0.869+9.91/ZFJ(L)
      RJM=RJ(L)+VJ(L)
      EX=VJ(L)+DCOL
      Q1=AMC*AFJ(L)
      EY=VJ(L)/5.
      PER=2.*SQRT(Q1*RJM)
      IF(RJ(L).LE.DCOL) GO TO 1
      EM=VJ(L)+(SQRT(1.+4.*Q1*RJ(L))-1.)/(2.*Q1)
      IF(EM-EX) 1,1,2
 2    FM=(EM-VJ(L))*EXP(2.*SQRT(Q1*(RJM-EM))-PER)
      GO TO 3
 1    EM=RJM-Q1/(DALF*DALF)
      IF(EM.GT.EX) EM=EX
      IF(EM.LE.EY) EM=EY
      FM=DCOL*EXP(2.*SQRT(Q1*(RJM-EM))+DALF*(EM-EX)-PER)
 3    ER=EY+(RJM-EY)*RRAN(ix)
      FR=FM*RRAN(ix)
      IF(ER.GT.EX) FTR=(ER-VJ(L))*EXP(2.*SQRT(Q1*(RJM-ER))-PER)
      IF(ER.LE.EX) FTR=DCOL*EXP(2.*SQRT(Q1*(RJM-ER))+DALF*(ER-EX)-PER)
      IF(FR-FTR) 4,4,3
 4    TKIN2=ER
       RETURN
      END
      
      
      
      FUNCTION GAMMA2(J,A,U,AM,AMC,AMF,RADNCL)
C  B[[[C[EH[E A[CO[[THO[ [[P[H[ [C[[CKAH[[ [PA[MEHTA COPTA J
C  ( B E[. (1./[P[BE[.[OCT.[[AHKA) M[B )
C  [[P[HA [E[EH[[ B TEX [E E[[H[[AX
       COMMON /BL1009/AFJ(35)/BL1015/RJ(35)/BL1014/GAM(35)
      COMMON /BL1011/VJ(35) /BL1005/AJ(35)/BL1006/ZJ(35)
      COMMON /BL1010/ZFJ(35)
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
C  [[[T[BAETC[ BO[MO[HOCT[ [O[[AP[EPHO[O [C[[CKAH[[ [PA[MEHTOB.
      IF(U.LE.0.) GO TO 12
      PER=2.*SQRT(AM*A*U)
      IF(J.EQ.33) GO TO 10
      DCOL=1.
      DALF=0.869+9.91/ZFJ(J)
      RJM=RJ(J)+VJ(J)
      EY=VJ(J)/5.
      IF(RJM.LE.EY.OR.J.EQ.1) GO TO 12
       CC=    (AJ(J)/AFJ(J))**0.6666667
       ALFA=1.+CC
       Q1=AMC*AFJ(J)
      BETA=DCOL
       Q3=(GAM(J)*AFJ(J)**.66667)*(ALFA/Q1**2)*(AFJ(J)/(AFJ(J)+AJ(J)))
      Q3=Q3*(3.1416*RADNCL**2)/(2.*41.5*3.1416**2)
      RJ(J)=RJ(J)-DCOL
      GAMMA3=0.
      IF(RJ(J).LE.0.) GO TO 8
       Q2=AMC*AFJ(J)*RJ(J)
       Q4=(2.*BETA*Q1-3.)/2.+Q2
       Q5=(2.*BETA*Q1-3.)*(SQRT(Q2)-.5)+2.*Q2
      GAMMA3=Q3*(Q4*EXP(-PER)+Q5*EXP(2.*SQRT(Q2)-PER))
 8    RJ(J)=RJ(J)+DCOL
      EX=VJ(J)+DCOL
      EM=RJM-Q1/(DALF*DALF)
      IF(EM.GE.EX) GO TO 13
      IF(EM.LT.EX.AND.EM.GT.EY) GO TO 16
      IF(EM.LE.EY) GO TO 17
 13   SQ=SQRT(Q1/(RJM-EX))
      F1X=DALF-SQ
      F2X=SQ/(2.*(RJM-EX))
      IF(F1X.GE.(0.5*F2X)) GO TO 14
      CSI=0.48/SQRT(0.5*F2X)
      F1CSI=DALF-SQRT(Q1/(RJM-EX+CSI))
      GAMMA4=Q3*Q1*Q1*DCOL*SQRT(6.2832/F2X)*
     #EXP(2.*SQRT(Q1*(RJM-EX))-PER-F1CSI*CSI)
      GO TO 15
 14   CSI=0.693/F1X
      SQCSI=SQRT(Q1/(RJM-EX+CSI))
      F2CSI=SQCSI/(2.*(RJM-EX+CSI))
      GAMMA4=Q3*2*Q1*Q1*(DCOL/F1X)*
     #EXP(2.*SQRT(Q1*(RJM-EX))-PER-F2CSI*CSI*CSI/2.)
      GO TO 15
 16   SQ=SQRT(Q1/(RJM-EM))
      F2M=SQ/(2.*(RJM-EM))
      CSI=0.48/SQRT(0.5*F2M)
      F1CSI=DALF-SQRT(Q1/(RJM-EM+CSI))
      GAMMA4=Q3*2*Q1*Q1*DCOL*SQRT(6.2832/F2M)*
     #EXP(DALF*(EM-EX)+2.*SQRT(Q1*(RJM-EM))-PER-F1CSI*CSI)
      GO TO 15
 17   GAMMA4=0.
      GO TO 15
 15   GAMMA2=GAMMA3+GAMMA4
      RETURN
 10   IF(RJ(J).LE.0.) GO TO 12
      Q1=2.*SQRT(AMF*A*RJ(J))
      Q2=1./(4.*3.1416)
      GAMMA2=(Q2/(AMF*A))*((Q1-1.)*EXP(Q1-PER)+EXP(-PER))
      RETURN
 12   GAMMA2=0.
      RETURN
       END
      
      
      
      SUBROUTINE FILLEV
C  [A[O[HEH[E BO[[[[[EHH[X [POBHE[ [E[K[X [PA[MEHTOB
      COMMON /EXILEV/ENLEV(20,35),SPLEV(20,35)
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      ENLEV(1,10)=3.56
      SPLEV(1,10)=1.
      ENLEV(1,11)=0.48
      SPLEV(1,11)=2.
      ENLEV(1,12)=0.98
      SPLEV(1,12)=3.
      ENLEV(1,13)=0.43
      SPLEV(1,13)=2.
      ENLEV(1,16)=3.37
      SPLEV(1,16)=5.
      ENLEV(2,16)=5.96
      SPLEV(2,16)=8.
      ENLEV(3,16)=6.18
      SPLEV(3,16)=1.
      ENLEV(4,16)=6.26
      SPLEV(4,16)=5.
      ENLEV(1,18)=0.72
      SPLEV(1,18)=3.
      ENLEV(2,18)=1.74
      SPLEV(2,18)=1.
      ENLEV(3,18)=2.15
      SPLEV(3,18)=3.
      ENLEV(4,18)=3.59
      SPLEV(4,18)=5.
      ENLEV(1,19)=2.13
      SPLEV(1,19)=2.
      ENLEV(2,19)=4.44
      SPLEV(2,19)=6.
      ENLEV(3,19)=5.02
      SPLEV(3,19)=4.
      ENLEV(4,19)=6.76
      SPLEV(4,19)=10.
      ENLEV(5,19)=7.29
      SPLEV(5,19)=6.
      ENLEV(6,19)=7.98
      SPLEV(6,19)=4.
      ENLEV(7,19)=8.56
      SPLEV(7,19)=6.
      ENLEV(1,20)=0.95
      SPLEV(1,20)=5.
      ENLEV(2,20)=1.67
      SPLEV(2,20)=5.
      ENLEV(3,20)=2.65
      SPLEV(3,20)=4.
      ENLEV(1,21)=2.00
      SPLEV(1,21)=2.
      ENLEV(2,21)=4.32
      SPLEV(2,21)=6.
      ENLEV(3,21)=4.80
      SPLEV(3,21)=4.
      ENLEV(4,21)=6.34
      SPLEV(4,21)=2.
      ENLEV(5,21)=6.48
      SPLEV(5,21)=8.
      ENLEV(6,21)=6.90
      SPLEV(6,21)=6.
      ENLEV(7,21)=7.50
      SPLEV(7,21)=4.
      ENLEV(8,21)=8.10
      SPLEV(8,21)=4.
      ENLEV(9,21)=8.42
      SPLEV(9,21)=6.
      ENLEV(10,21)=8.66
      SPLEV(10,21)=8.
      ENLEV(1,22)=4.44
      SPLEV(1,22)=5.
      ENLEV(1,23)=3.09
      SPLEV(1,23)=2.
      ENLEV(2,23)=3.68
      SPLEV(2,23)=4.
      ENLEV(3,23)=3.85
      SPLEV(3,23)=6.
      ENLEV(1,24)=6.09
      SPLEV(1,24)=3.
      ENLEV(2,24)=6.69
      SPLEV(2,24)=8.
      ENLEV(3,24)=6.96
      SPLEV(3,24)=6.
      ENLEV(4,24)=7.34
      SPLEV(4,24)=5.
      ENLEV(1,26)=2.31
      SPLEV(1,26)=1.
      ENLEV(2,26)=3.95
      SPLEV(2,26)=3.
      ENLEV(3,26)=4.92
      SPLEV(3,26)=1.
      ENLEV(4,26)=5.11
      SPLEV(4,26)=5.
      ENLEV(5,26)=5.69
      SPLEV(5,26)=3.
      ENLEV(6,26)=5.83
      SPLEV(6,26)=7.
      ENLEV(7,26)=6.20
      SPLEV(7,26)=3.
      ENLEV(8,26)=6.44
      SPLEV(8,26)=7.
      ENLEV(9,26)=7.03
      SPLEV(9,26)=5.
      ENLEV(1,27)=5.28
      SPLEV(1,27)=8.
      ENLEV(2,27)=6.32
      SPLEV(2,27)=4.
      ENLEV(3,27)=7.22
      SPLEV(3,27)=10.
      ENLEV(4,27)=7.57
      SPLEV(4,27)=8.
      ENLEV(5,27)=8.31
      SPLEV(5,27)=2.
      ENLEV(6,27)=8.57
      SPLEV(6,27)=4.
      ENLEV(7,27)=9.15
      SPLEV(7,27)=14.
      ENLEV(8,27)=9.79
      SPLEV(8,27)=14.
      ENLEV(9,27)=10.00
      SPLEV(9,27)=8.
      ENLEV(1,28)=0.12
      SPLEV(1,28)=1.
      ENLEV(2,28)=0.30
      SPLEV(2,28)=7.
      ENLEV(3,28)=0.40
      SPLEV(3,28)=3.
      ENLEV(1,29)=5.22
      SPLEV(1,29)=8.
      ENLEV(2,29)=6.18
      SPLEV(2,29)=4.
      ENLEV(3,29)=6.83
      SPLEV(3,29)=10.
      ENLEV(4,29)=7.28
      SPLEV(4,29)=8.
      ENLEV(1,30)=6.10
      SPLEV(1,30)=8.
      ENLEV(2,30)=6.92
      SPLEV(2,30)=5.
      ENLEV(3,30)=7.12
      SPLEV(3,30)=3.
      ENLEV(1,31)=0.87
      SPLEV(1,31)=2.
      ENLEV(2,31)=3.06
      SPLEV(2,31)=2.
      ENLEV(3,31)=3.84
      SPLEV(3,31)=6.
      ENLEV(1,32)=1.98
      SPLEV(1,32)=5.
      ENLEV(2,32)=3.57
      SPLEV(2,32)=10.
      ENLEV(3,32)=3.92
      SPLEV(3,32)=5.
      ENLEV(4,32)=4.46
      SPLEV(4,32)=3.
      ENLEV(5,32)=5.10
      SPLEV(5,32)=7.
      ENLEV(6,32)=5.33
      SPLEV(6,32)=13.
      ENLEV(7,32)=5.53
      SPLEV(7,32)=5.
      ENLEV(8,32)=6.20
      SPLEV(8,32)=3.
      ENLEV(9,32)=6.38
      SPLEV(9,32)=12.
      ENLEV(10,32)=6.88
      SPLEV(10,32)=1.
      RETURN
      END
      
      
      
      SUBROUTINE FISION(ENEXT,ATWGHT,CHARGE,PFX,PFY,PFZ,KSTART)
C  FISSION  (IN - MEV, OUT - GEV)
      COMMON/BLOKC/SPTU(10,500)/BLANGL/ANGL(4)
      COMMON/BLFIS/FIS,EF,SNT(6,6)
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      UT=ENEXT
      AT=ATWGHT
      ZT=CHARGE
      VFX=PFX/(940.*AT)
      VFY=PFY/(940.*AT)
      VFZ=PFZ/(940.*AT)
      SNT(1,2)=ENEXT*0.001
      SNT(2,2)=ATWGHT
      SNT(3,2)=CHARGE
      SNT(4,2)=PFX*0.001
      SNT(5,2)=PFY*0.001
      SNT(6,2)=PFZ*0.001
      JJ=0
 1    CONTINUE
      JJ=JJ+1
      IF(JJ.GT.100) GO TO 10
      AR=AFIS(UT,AT,ZT,EF)
      IA1=INT(AR+0.5)
      ZR=ZFIS(AT,ZT,AR)
      IZ1=INT(ZR+0.5)
      AF1=IA1
      ZF1=IZ1
      AF2=AT-AF1
      ZF2=ZT-ZF1
      IA2=INT(AF2+0.5)
      IZ2=INT(ZF2+0.5)
      EKF=FISKIN(UT,AT,ZT,AF1,AF2,ZF1,ZF2)
      AMAS1=940.*AF1
      AMAS2=940.*AF2
      EK1=AF2*EKF/AT
      EK2=AF1*EKF/AT
      CALL DELAM(AT,ZT,DL0,DSH0,BAR0)
      CALL DELAM(AF1,ZF1,DL1,DSH1,BAR1)
      CALL DELAM(AF2,ZF2,DL2,DSH2,BAR2)
      ET=UT+DL0-DL1-DL2
C      PRINT 2,AT,ZT,AF1,ZF1,AF2,ZF2,UT,DL0,DL1,DL2,ET,EKF
C 2    FORMAT(2X,'AT,ZT=',2F6.1,' AF1,ZF1=',2F6.1,' AF2,ZF2=',2F6.1,
C     *' UT,DL0=',2F6.1,' DL1,DL2=',2F6.1,' ET,EKF=',2F6.1)
      IF(ET.LT.EKF) GO TO 1
      U1=0.001*(ET-EKF)*AF1/AT
      U2=0.001*(ET-EKF)*AF2/AT
      CALL ISANGL
      VM=SQRT(2.*EK1/AMAS1)
      VX1=VM*ANGL(4)*ANGL(3)
      VY1=VM*ANGL(4)*ANGL(2)
      VZ1=VM*ANGL(1)
      PX1=AMAS1*(VX1+VFX)*0.001
      PY1=AMAS1*(VY1+VFY)*0.001
      PZ1=AMAS1*(VZ1+VFZ)*0.001
      TEMP=-1.*(AMAS1/AMAS2)
      VX2=TEMP*VX1
      VY2=TEMP*VY1
      VZ2=TEMP*VZ1
      PX2=AMAS2*(VX2+VFX)*0.001
      PY2=AMAS2*(VY2+VFY)*0.001
      PZ2=AMAS2*(VZ2+VFZ)*0.001
C  BEFORE EVAPORATION
      SNT(1,3)=U1
      SNT(2,3)=AF1
      SNT(3,3)=ZF1
      SNT(4,3)=PX1
      SNT(5,3)=PY1
      SNT(6,3)=PZ1
      SNT(1,4)=U2
      SNT(2,4)=AF2
      SNT(3,4)=ZF2
      SNT(4,4)=PX2
      SNT(5,4)=PY2
      SNT(6,4)=PZ2
C  AFTER EVAPORATION
C      PRINT 3,KSTART
C 3    FORMAT(2X,' BEFORE EVAP.FRAGM. KSTART=',I4)
      CALL EVAPIN(U1,AF1,ZF1,PX1,PY1,PZ1,KSTART)
      SNT(1,5)=U1
      SNT(2,5)=AF1
      SNT(3,5)=ZF1
      SNT(4,5)=PX1
      SNT(5,5)=PY1
      SNT(6,5)=PZ1
C      PRINT 4,KSTART
C 4    FORMAT(2X,'AFT.EVAP.FIRST, BEF.EVAP.SEC.FRAGM. KSTART=',I4)
      CALL EVAPIN(U2,AF2,ZF2,PX2,PY2,PZ2,KSTART)
      SNT(1,6)=U2
      SNT(2,6)=AF2
      SNT(3,6)=ZF2
      SNT(4,6)=PX2
      SNT(5,6)=PY2
      SNT(6,6)=PZ2
C      PRINT 5,KSTART
C 5    FORMAT(2X,' AFTER EVAP.FRAGM. KSTART=',I4)
      RETURN
 10   PRINT 11,ET,EKF,JJ,UT,AT,ZT
 11   FORMAT(2X,'ERROR IN FISION: ET=',F6.1,'.LT.EKF=',F6.1,
     *' JJ=',I4,' UT,AT,ZT=',3F6.1)
      RETURN
      END
      
      
      
      FUNCTION YIELDF(A0,A,F,AS,A1,A2,SIGS,SIG1,SIG2)
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      XSIM=EXP(-0.5*(A-AS)**2/SIGS**2)
      XASIM=EXP(-0.5*(A-A2)**2/SIG2**2)
     *     +EXP(-0.5*(A-(A0-A2))**2/SIG2**2)
     * +0.5*EXP(-0.5*(A-A1)**2/SIG1**2)
     * +0.5*EXP(-0.5*(A-(A0-A1))**2/SIG1**2)
      IF(F.GT.1000.) YIELDF=XSIM
      IF(F.LT.0.001) YIELDF=XASIM
      IF(F.GE.0.001.AND.F.LE.1000) YIELDF=F*XSIM+XASIM
      RETURN
      END
      
      
      
      FUNCTION ZFIS(A,Z,AF)
      COMMON /BLFIS0/F,FH,AS,A1,A2,SIGS,SIG1,SIG2,ZP,SIGZ,EAVR,EMAX
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      IF(AF.GE.134.) DZ=-0.45
      IF(AF.LE.(A-134.)) DZ=0.45
      IF(AF.GT.(A-134.).AND.AF.LT.134.) DZ=-0.45*(AF-0.5*A)/(134-0.5*A)
      ZP=AF*Z/A+DZ
      SIGZ=0.6
 3    CALL RANNOR(X,Y)
      XZ=ZP+SIGZ*X
      IF(XZ.LT.1.OR.XZ.GT.(Z-1).OR.XZ.GT.AF) XZ=ZP+SIGZ*Y
      IF(XZ.LT.1.OR.XZ.GT.(Z-1).OR.XZ.GT.AF) GO TO 3
      ZFIS=INT(XZ+0.5)
      RETURN
      END
       
       
       
       SUBROUTINE EVAPIN(ENEXT,ATWGHT,CHARGE,PNX,PNY,PNZ,KSTART)
C   EVAPORATION WITHOUT FISSION FOR SUBR.-FISION
       COMMON /BL0999/RNCL/BL1000/AM,AMF,AMCB(35)/AMNUCL/AMNEV(3)
       COMMON /BL1001/T1Y(130)/BL1002/T2XY(200)/BL1005/AJ(35)
     */BL1006/ZJ(35)/BL1014/GAM(35)/BL1011/VJ(35)/BL1015/RJ(35)
      COMMON /BL1003/U,A,Z/BL1009/AFJ(35)/BL1010/ZFJ(35)
      COMMON/BLOKC/SPTU(10,500)/BLANGL/ANGL(4)/BL1100/AMS,AMFS
      COMMON/BLFIS/FIS,EF,SNT(6,6) /BLCAN/ICAN
      COMMON /EXILEV/ ENLEV(20,35),SPLEV(20,35)
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
       DIMENSION GJ(35),BJ(35),GG(35),VN(3),PN(3),PNL(3),PP(3),PPL(3)
      DIMENSION DSH2(35)
      FIS=1.
      ILEVEL=3
      CALL DATEVP(ILEVEL)
      U=ENEXT*1000.
      A=ATWGHT
      Z=CHARGE
      AMNEV(1)=0.
      AMNEV(2)=0.
      AMNEV(3)=0.
      REMN=940.*A
      PNL(1)=PNX*1000.
      PNL(2)=PNY*1000.
      PNL(3)=PNZ*1000.
      ENL=SQRT(PNL(1)**2+PNL(2)**2+PNL(3)**2+REMN**2)
      VN(1)=PNL(1)/ENL
      VN(2)=PNL(2)/ENL
      VN(3)=PNL(3)/ENL
      DO 20 K=KSTART,500
      IN=35
      DO 4 I=1,IN
      VJ(I)=0.
      AMCB(I)=0.
C      DSH2(I)=0.
 4    RJ(I)=-1000.
      CALL DELAM(A,Z,DL1,DSHEL1,BAR1)
      AM=AMS
C      IF(U.GT.0.) AM=AMS*(1.+(1.-EXP(-0.06*U))*DSHEL1/U)
      DO 6 I=1,IN
      IF(I.EQ.33) GO TO 6
      AFJ(I)=A-AJ(I)
      ZFJ(I)=Z-ZJ(I)
      IF(ICAN.EQ.0.AND.I.GT.6) GO TO 6
      IF(ICAN.EQ.1.AND.I.GT.32) GO TO 6
      IF(ICAN.EQ.2.AND.I.GT.6.AND.I.LE.32) GO TO 6
      IF(AFJ(I).LT.ZFJ(I)) GO TO 6
      IF(AFJ(I).LT.AJ(I).OR.ZFJ(I).LT.ZJ(I)) GO TO 6
      JAF=INT(AFJ(I)+0.5)
      JZF=INT(ZFJ(I)+0.5)
      JNF=JAF-JZF
      ODD=11.*(2+2*(JNF/2)-JNF+2*(JZF/2)-JZF)/SQRT(AFJ(I))
      IF(RNCL.GE.0.1) RADNCL=RNCL
      IF(RNCL.LT.0.1)
     *RADNCL=2.173*(1+0.006103*ZJ(I)*ZFJ(I))/(1+0.009443*ZJ(I)*ZFJ(I))
      VJ(I)=COLOMB(I,RADNCL)
      CALL DELAM(AFJ(I),ZFJ(I),DL2,DSHEL2,BAR2)
C      DSH2(I)=DSHEL2
      CALL DELAM(AJ(I),ZJ(I),DL3,DSHEL3,BAR3)
      BJ(I)=DL2+DL3-DL1
c ******** in accordance with DELA8
      IF(A.LE.55.) RJ(I)=U-(BJ(I)+VJ(I))-ODD
      TCOEF=1.-(A-55.)/10.
      IF(A.GT.55.AND.A.LT.65.) RJ(I)=U-(BJ(I)+VJ(I))-ODD*TCOEF
      IF(A.GE.65.) RJ(I)=U-(BJ(I)+VJ(I))
c ********
C     RJ(I)=U-(BJ(I)+VJ(I))-ODD
C     RJ(I)=U-(BJ(I)+VJ(I))
C      IF(RJ(I).GT.0.)
C     *AMCB(I)=AMS*(1.+(1.-EXP(-0.06*RJ(I)))*DSHEL2/RJ(I))
 6    CONTINUE
      RJ(33)=-1001.
       DO 7 I=1,IN
      IF(RNCL.GE.0.1) RADNCL=RNCL
      IF(RNCL.LT.0.1)
     *RADNCL=2.173*(1+0.006103*ZJ(I)*ZFJ(I))/(1+0.009443*ZJ(I)*ZFJ(I))
      AMC=AM
C      AMC=AMCB(I)
      GJ(I)=GAMMAC(I,A,U,AM,AMC,AMF,RADNCL)
      RJI=RJ(I)
      GAMI=GAM(I)
      DO 71 IL=1,20
      IF(SPLEV(IL,I).LT.0.1) GO TO 70
      GAM(I)=SPLEV(IL,I)*AJ(I)
      RJ(I)=RJI-ENLEV(IL,I)
C      IF(RJ(I).GT.0.)
C     *AMC=AMS*(1.+(1.-EXP(-0.06*RJ(I)))*DSH2(I)/RJ(I))
      GJ(I)=GJ(I)+GAMMAC(I,A,U,AM,AMC,AMF,RADNCL)
 71   CONTINUE
 70   RJ(I)=RJI
      GAM(I)=GAMI
7     CONTINUE
      G=0.
      DO 10 I=1,IN
10    G=G+GJ(I)
      IF(G) 11,11,12
 11   PNX=PNL(1)*0.001
      PNY=PNL(2)*0.001
      PNZ=PNL(3)*0.001
      ENEXT=U*0.001
      ATWGHT=A
      CHARGE=Z
      REMN=940.*A
      SPTU(9,K)=REMN*0.001
      SPTU(8,K)=Z
      CALL PINT(PNL,CT,ST,CF,SF,T,REMN)
      SPTU(4,K)=CT
      SPTU(5,K)=SF
      SPTU(6,K)=CF
      SPTU(7,K)=T*0.001
      SPTU(1,K)=1.
      SPTU(3,K)=FIS
C-(Sobolevsky)-res. excitation in MeV
      SPTU(10,K)=U
C------------------------------------
      KSTART=K+1
      RETURN
 12   CONTINUE
      DO 13 J=2,IN
13    GJ(J)=GJ(J-1)+GJ(J)
      BB=RRAN(ix)
      B=BB*G
      DO 14 J=1,IN
      IF(B-GJ(J)) 15,14,14
15    LM=J
      GO TO 16
 14   CONTINUE
 16     IF(LM-33) 18,17,18
 17   PRINT 300
 300  FORMAT(2X,' ERROR IN EPAPIN - FALSE FISSION')
      GO TO 11
 18   AMC=AM
C18    AMC=AMCB(LM)
      EP1=TKIN(LM,AMC,AMF)
      EP2=ZJ(LM)
      EP3=940.*AJ(LM)
      U=U-BJ(LM)-EP1
      A=AFJ(LM)
      Z=ZFJ(LM)
C  VPM-OTHOC[TE[[HA[ CKOPOCT[ [PA[MEHTA
      VPM=SQRT((2.*EP1)/(EP3*AFJ(LM)/(AFJ(LM)+AJ(LM))))
      CALL ISANGL
C  [M[[[[C [PA[MEHTA B C.[.[.
      PP(1)=VPM*ANGL(4)*ANGL(3)*EP3/(1.+AJ(LM)/A)
      PP(2)=VPM*ANGL(4)*ANGL(2)*EP3/(1.+AJ(LM)/A)
      PP(3)=VPM*ANGL(1)*EP3/(1.+AJ(LM)/A)
      PN(1)=-PP(1)
      PN(2)=-PP(2)
      PN(3)=-PP(3)
      EP=SQRT(PP(1)**2+PP(2)**2+PP(3)**2+EP3**2)
      EN=SQRT(PN(1)**2+PN(2)**2+PN(3)**2+(940.*A)**2)
      CALL CLPV(PP,VN,PPL,EP)
      CALL CLPV(PN,VN,PNL,EN)
      ENL=SQRT(PNL(1)**2+PNL(2)**2+PNL(3)**2+(940.*A)**2)
      VN(1)=PNL(1)/ENL
      VN(2)=PNL(2)/ENL
      VN(3)=PNL(3)/ENL
      CALL PINT(PPL,CT,ST,CF,SF,T,EP3)
      SPTU(4,K)=CT
      SPTU(5,K)=SF
      SPTU(6,K)=CF
      SPTU(7,K)=T*0.001
      SPTU(8,K)=EP2
      SPTU(9,K)=EP3*0.001
      SPTU(2,K)=LM
20    CONTINUE
      PRINT 21,U,A,Z
21    FORMAT(35X,37HMASSIV SPT EXCEEDED AFTER EVAPORATION/40X,2HU=,
     *F10.5,4H  A=,F5.1,4H  Z=,F4.1)
       RETURN
       END
      
      
      
      FUNCTION AFIS(E,A,Z,EF)
C IMPROVED VERSION TAKING INTO ACCOUNT DEEXCIT. BEFORE FISSION
      COMMON /BLFIS0/F,FH,AS,A1,A2,SIGS,SIG1,SIG2,ZP,SIGZ,EAVR,EMAX
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      A1=134.
      A2=141.
      AS=A/2.
      IZ=INT(Z+0.5)
      IF(A.GE.235.) SIG2=5.6
      IF(A.LT.235.) SIG2=5.6+0.096*(A-235.)
      SIG1=0.5*SIG2
C     IF(E.LE.15.) SIGS=EXP((0.2412*E+0.8252)/2.)
      SIGS=EXP((0.01106*E+4.2772)/2.)
      IF(SIGS.GT.20.) SIGS=20.
      YA=2.*EXP(-0.5*(A2-AS)**2/SIG2**2)+EXP(-0.5*(A1-AS)**2/SIG1**2)
      YS=EXP(-0.5*(AS-(A1+A2)/2.)**2/SIGS**2)
      IF(IZ.GE.90) GO TO 1
      IF(IZ.EQ.89) GO TO 2
      IF(IZ.LE.88.AND.IZ.GE.82) GO TO 3
      IF(IZ.LT.82) GO TO 4
 1    IF(E.LE.16.25) F=EXP(0.5385*E-9.9564)
      IF(E.GT.16.25) F=EXP(0.09197*E-2.7003)
      GO TO 60
 2    F=EXP(0.09197*E-1.0808)
      GO TO 60
 3    X=EF-7.5
      IF(X.LT.0.) X=0.
      F=EXP(0.09197*(E-X)-1.0808)
      T1=1.03*F-YA
      T2=1.-YS*F
      IF(T2.LE.0.0001) T2=0.0001
      IF(T1.LE.0.0001) T1=0.0001
      FH=T1/T2
      IF(A.LT.227.) FH=FH*EXP(0.3*(227-A))
      GO TO 5
 4    FH=1001.
      GO TO 5
 60   T1=1.03*F-YA
      T2=1.-YS*F
      IF(T2.LE.0.0001) T2=0.0001
      IF(T1.LE.0.0001) T1=0.0001
      FH=T1/T2
 5    CONTINUE
      C2A=A2+3.72*SIG2
      C2S=AS+3.72*SIGS
      C2=AMAX1(C2A,C2S)
      IF(FH.GT.1000.) C2=C2S
      IF(FH.LT.0.001) C2=C2A
      C1=A-C2
      IF(C1.LT.30.) C2=A-30.
      IF(C1.LT.30.) C1=30.
      AM1=(AS+A1)/2.
      AM2=(A1+A2)/2.
      XM1=YIELDF(A,AS,FH,AS,A1,A2,SIGS,SIG1,SIG2)
      XM2=YIELDF(A,AM1,FH,AS,A1,A2,SIGS,SIG1,SIG2)
      XM3=YIELDF(A,A1,FH,AS,A1,A2,SIGS,SIG1,SIG2)
      XM4=YIELDF(A,AM2,FH,AS,A1,A2,SIGS,SIG1,SIG2)
      XM5=YIELDF(A,A2,FH,AS,A1,A2,SIGS,SIG1,SIG2)
      PMAX=AMAX1(XM1,XM2,XM3,XM4,XM5)
 6    X=C1+RRAN(ix)*(C2-C1)
      PX=YIELDF(A,X,FH,AS,A1,A2,SIGS,SIG1,SIG2)
      Y=RRAN(ix)
      IF(Y-PX/PMAX)7,7,6
 7    AFIS=INT(X+0.5)
      RETURN
      END
      
      
      
      FUNCTION FISKIN(E,A,Z,AF1,AF2,ZF1,ZF2)
C IMPROVED VERSION TAKING INTO ACCOUNT THE CHANNELS SEPARATELY
      COMMON /BLFIS0/F,FH,AS,A1,A2,SIGS,SIG1,SIG2,ZP,SIGZ,EAVR,EMAX
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      AFM=AF1
      IF(AF2.GT.AF1) AFM=AF2
      IF(AFM.LT.(A/2.)) AFM=A-AFM
      EAVR=0.1071*Z*Z/A**0.33333+22.2
      P1=0.5*EXP(-0.5*(AFM-A1)**2/SIG1**2)
      P2=EXP(-0.5*(AFM-A2)**2/SIG2**2)
      PS=FH*EXP(-0.5*(AFM-AS)**2/SIGS**2)
      PAS=P1+P2
      IF(FH.GT.1000.) PAS=0.
      IF(FH.LT.0.001) PS=0.
      PSY=PS/(PAS+PS)
      ISY=0
      BRND=RRAN(ix)
      IF(BRND.LE.PSY) ISY=1
      IF(ISY.EQ.0) SIGE=10.
      IF(ISY.GE.1) SIGE=8.
      PPAS=2.*0.5*SIG1+2.*SIG2
      PPSY=FH*SIGS
      PPP=PPAS+PPSY
      XAS=PPAS/PPP
      XSY=PPSY/PPP
      EAVRSY=EAVR-12.5*XAS
      EAVRAS=EAVR+12.5*XSY
      IF(ISY.EQ.0) GO TO 10
      IF(ISY.GE.1) GO TO 11
 10   EAV0=EAVRAS
      A11=A1-0.7979*SIG1
      A12=A1+0.7979*SIG1
      A21=A2-0.7979*SIG2
      A22=A2+0.7979*SIG2
      EEF11=EEFIS(A,A11,ISY)
      EEF12=EEFIS(A,A12,ISY)
      EEF21=EEFIS(A,A21,ISY)
      EEF22=EEFIS(A,A22,ISY)
      TEMPL=0.5*SIG1*EEF11+0.5*SIG1*EEF12+SIG2*EEF21+SIG2*EEF22
      EMAX=EAV0*PPAS/TEMPL
      GO TO 12
 11   EAV0=EAVRSY
      AS0=AS+0.7979*SIGS
      EEF00=EEFIS(A,AS0,ISY)
      TEMPL=FH*SIGS*EEF00
      EMAX=EAV0*PPSY/TEMPL
 12   EKAFM=EMAX*EEFIS(A,AFM,ISY)
      JJJJ=0
 3    CALL RANNOR(X,Y)
      JJJJ=JJJJ+1
      IF(JJJJ.GT.100) EKIN=EAVR
      IF(JJJJ.GT.100) GO TO 650
      EKIN=EKAFM+SIGE*X
      IF(EKIN.LT.(EAVR-3.72*SIGE).OR.EKIN.GT.(EAVR+3.72*SIGE))
     *EKIN=EKAFM+SIGE*Y
      IF(EKIN.LT.(EAVR-3.72*SIGE).OR.EKIN.GT.(EAVR+3.72*SIGE))
     *GO TO 3
 650  FISKIN=EKIN
      RETURN
      END
      
      
      
      FUNCTION EEFIS(A,AFM,ISY)
C IMPROVED VERSION TAKING INTO ACCOUNT THE CHANNELS SEPARATELY
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      IF(ISY.GE.1) GO TO 1
      A00=134.
      B1=23.5
      GO TO 2
 1    B1=5.32
      A00=A/2.
 2    TEMP1=10./A
      A0010=A00+10.
      AS=A/2.
      IF(AFM.GE.AS.AND.AFM.LE.A0010) EE=1.-B1*((AFM-A00)/A)**2
      IF(AFM.GT.A0010) EE=1.-B1*TEMP1**2-2.*TEMP1*B1*(AFM-A0010)/A
      EEFIS=EE
      RETURN
      END
      
      
      
      SUBROUTINE DELAM(ANUCL,ZNUCL,DM,SH,BARIER)
      COMMON /FUSR/ BARR,SMASS,SHLL
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      IA=INT(ANUCL+0.5)
      IZ=INT(ZNUCL+0.5)
      IN=IA-IZ
      IF(IZ.LT.0.OR.IN.LT.0) GO TO 8
      ODD=-11.*(1+2*(IN/2)-IN+2*(IZ/2)-IZ)/SQRT(ANUCL)
      IF(IA.LT.65) GO TO 7
      IF(IN.GT.240.OR.IZ.GT.240) GO TO 7
      CM=0.
      SHLL=0.
      CBAR=0.
      NBAR=0
      CALL LYMASM(IZ,IA,CM,CBAR,NBAR)
      IF(NBAR.EQ.1) CBAR=BARR
      DM=CM
      SH=SHLL+ODD
      BARIER=CBAR
c ************* inclusion DELA8
      CALL DELA8(ANUCL,ZNUCL,DM8,SH8)
      DM=DM8
c *************
      RETURN
 7    CONTINUE
      CALL DELA1(ANUCL,ZNUCL,DM1,SH1)
      BAR1=1001.
      IF(IA.GE.65) BAR1=BARIE1(ANUCL,ZNUCL)
      BARIER=BAR1
      DM=DM1
      SH=SH1
c ************* inclusion DELA8
      IF(ANUCL.GT.55.) CALL DELA8(ANUCL,ZNUCL,DM8,SH8)
      TCOEF=1.-(ANUCL-55.)/10.
      IF(ANUCL.GT.55.AND.ANUCL.LT.65.) DM=DM1*TCOEF+DM8*(1.-TCOEF)
c *************
      RETURN
 8    PRINT 9,IZ,IN
 9    FORMAT(2X,'ERROR IN DELAM: IZ=',I5,' IN=',I5)
      RETURN
      END
      
      
      
      SUBROUTINE DELA8(X,Y,DELTA8,DSHEL8)
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      DSHEL8=0.
      I=INT(X+0.5)
      J=INT(Y+0.5)
      L=I-J
      DELTA8=1001.
      IF(I.LE.1.OR.J.LE.0.OR.L.LE.0) RETURN
      X03=X**0.3333333
      X06=X03*X03
      X13=X06*X06
      Y13=Y**1.3333333
      D=(X-2.*Y)/X
      DR=1-0.62025/X06
      ES=(25.8357-44.2355*D*D)*DR*DR*X06
      EC=0.779*(Y*(Y-1.)/X03)*(1.-1.5849/X06+1.2273/X+1.5772/X13)
      EEX=-0.4323*(Y13/X03)*(1.-0.57811/X03-0.14518/X06+0.49597/X)
      EV=X*(-17.0354+31.4506*D*D)
      DELTA8=8.071*X-0.783*Y+EV+ES+EC+EEX
      RETURN
      END
      
      
      
      SUBROUTINE DELA1(X,Y,DELTAM,DSHEL)
C  B[[[C[EH[E (M-A) CO[[ACTHO KAMEPOH[-57 (B [[[EPO[H[X E[[H[[AX)
C XOPO[EE COOTBETCTB[E [P[ Z.LE.8 [ CO[[ACOBAH[E [P[ [O[[[[X Z
C      ([[[ [E[E[ [C[APEH[[ [[ [E[K[X [[EP [ [X CAM[X)
      COMMON /BL1001/T1Y(130) /BL1002/T2XY(200)
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      I=INT(X+0.5)
      J=INT(Y+0.5)
      L=I-J
      DELTAM=1001.
      IF(J.EQ.0) DELTAM=8.071*I
      IF(L.EQ.0) DELTAM=7.289*J
      IF(J.GT.10)  GO TO 11
      IF(I.EQ.1.AND.J.EQ.0) DELTAM=8.071
      IF(J.EQ.0)   GO TO 11
      GO TO (1,2,3,4,5,6,7,8,9,10),J
 1    IF(I.EQ.1.AND.J.EQ.1) DELTAM=7.289
      IF(I.EQ.2.AND.J.EQ.1) DELTAM=13.136
      IF(I.EQ.3.AND.J.EQ.1) DELTAM=14.950
      IF(I.EQ.4.AND.J.EQ.1) DELTAM=25.920
      IF(I.EQ.5.AND.J.EQ.1) DELTAM=33.790
      GO TO 11
 2    IF(I.EQ.3.AND.J.EQ.2) DELTAM=14.931
      IF(I.EQ.4.AND.J.EQ.2) DELTAM=2.425
      IF(I.EQ.5.AND.J.EQ.2) DELTAM=28.150
      IF(I.EQ.6.AND.J.EQ.2) DELTAM=17.593
      IF(I.EQ.7.AND.J.EQ.2) DELTAM=26.111
      IF(I.EQ.8.AND.J.EQ.2) DELTAM=31.599
      GO TO 11
 3    IF(I.EQ.5.AND.J.EQ.3) DELTAM=28.340
      IF(I.EQ.6.AND.J.EQ.3) DELTAM=14.087
      IF(I.EQ.7.AND.J.EQ.3) DELTAM=14.908
      IF(I.EQ.8.AND.J.EQ.3) DELTAM=20.947
      GO TO 11
 4    IF(I.EQ.6.AND.J.EQ.4) DELTAM=18.375
      IF(I.EQ.7.AND.J.EQ.4) DELTAM=15.769
      IF(I.EQ.8.AND.J.EQ.4) DELTAM=4.942
      IF(I.EQ.9.AND.J.EQ.4) DELTAM=11.348
      IF(I.EQ.10.AND.J.EQ.4) DELTAM=12.607
      IF(I.EQ.11.AND.J.EQ.4) DELTAM=20.176
      IF(I.EQ.12.AND.J.EQ.4) DELTAM=25.072
      GO TO 11
 5    IF(I.EQ.9.AND.J.EQ.5) DELTAM=12.415
      IF(I.EQ.10.AND.J.EQ.5) DELTAM=12.052
      IF(I.EQ.11.AND.J.EQ.5) DELTAM=8.668
      IF(I.EQ.12.AND.J.EQ.5) DELTAM=13.370
      IF(I.EQ.13.AND.J.EQ.5) DELTAM=16.562
      IF(I.EQ.14.AND.J.EQ.5) DELTAM=24.230
      GO TO 11
 6    IF(I.EQ.10.AND.J.EQ.6) DELTAM=15.702
      IF(I.EQ.11.AND.J.EQ.6) DELTAM=10.650
      IF(I.EQ.12.AND.J.EQ.6) DELTAM=0.
      IF(I.EQ.13.AND.J.EQ.6) DELTAM=3.125
      IF(I.EQ.14.AND.J.EQ.6) DELTAM=3.020
      GO TO 11
 7    IF(I.EQ.13.AND.J.EQ.7) DELTAM=5.346
      IF(I.EQ.14.AND.J.EQ.7) DELTAM=2.864
      IF(I.EQ.15.AND.J.EQ.7) DELTAM=0.102
      IF(I.EQ.16.AND.J.EQ.7) DELTAM=5.683
      IF(I.EQ.17.AND.J.EQ.7) DELTAM=7.871
      IF(I.EQ.18.AND.J.EQ.7) DELTAM=13.274
      IF(I.EQ.19.AND.J.EQ.7) DELTAM=15.790
      GO TO 11
 8    IF(I.EQ.14.AND.J.EQ.8) DELTAM=8.007
      IF(I.EQ.15.AND.J.EQ.8) DELTAM=2.863
      IF(I.EQ.16.AND.J.EQ.8) DELTAM=-4.737
      IF(I.EQ.17.AND.J.EQ.8) DELTAM=-0.809
      IF(I.EQ.18.AND.J.EQ.8) DELTAM=-0.783
      IF(I.EQ.19.AND.J.EQ.8) DELTAM=3.332
      IF(I.EQ.20.AND.J.EQ.8) DELTAM=3.799
      GO TO 11
 9    IF(I.EQ.17.AND.J.EQ.9) DELTAM=1.952
      IF(I.EQ.18.AND.J.EQ.9) DELTAM=0.873
      IF(I.EQ.19.AND.J.EQ.9) DELTAM=-1.486
      IF(I.EQ.20.AND.J.EQ.9) DELTAM=-0.016
      IF(I.EQ.21.AND.J.EQ.9) DELTAM=-0.046
      GO TO 11
 10   IF(I.EQ.19.AND.J.EQ.10) DELTAM=1.752
      IF(I.EQ.20.AND.J.EQ.10) DELTAM=-7.041
 11   DSHEL=0.
      IF(I.LE.1.OR.J.LE.0.OR.L.LE.0) RETURN
      X03=X**0.3333333
      X06=X03*X03
      X13=X06*X06
      Y13=Y**1.3333333
      D=(X-2.*Y)/X
      DR=1-0.62025/X06
      ES=(25.8357-44.2355*D*D)*DR*DR*X06
      EC=0.779*(Y*(Y-1.)/X03)*(1.-1.5849/X06+1.2273/X+1.5772/X13)
      EEX=-0.4323*(Y13/X03)*(1.-0.57811/X03-0.14518/X06+0.49597/X)
      EV=X*(-17.0354+31.4506*D*D)
      IF(J.LE.0.OR.J.GT.130.OR.L.LE.0.OR.L.GT.200) GO TO 22
      T1=T1Y(J)
      T2=T2XY(L)
      GO TO 23
 22   T1=0.
      T2=0.
 23   CONTINUE
      DELTA0=8.071*X-0.783*Y+EV+ES+EC+EEX+T1+T2
      IF(DELTAM.LT.1000.) GO TO 26
      DSHEL=T1+T2
      DELTAM=DELTA0
      RETURN
 26   DSHEL=DELTAM-DELTA0+T1+T2
      RETURN
      END
      
      
      
      BLOCK DATA DMAS02
C [A[A[TC[ [APAMETP[ OCTAT.[[PA,[C[AP[[[[XC[ [PA[MEHTOB,O[PATHO[O
C CE[EH[[([OCTPOBCK[[-59),O[O[O[.[O[PABK[(KAMEPOH-57)
      COMMON /BL1001/T1Y(130)/BL1002/T2XY(200)/BL1014/GAM(35)
      COMMON /BL0999/RNCL  /BL1100/AMS,AMFS /BLGAM/GAM1(35),GAM2(35)
      COMMON /BL1005/AJ(35) /BL1006/ZJ(35) /BLCAN/ICAN /BLSEC/ISEC
      COMMON /EXILEV/ ENLEV(20,35),SPLEV(20,35)
      DIMENSION TT11(65),TT12(65)
      DIMENSION TT21(50),TT22(50),TT23(50),TT24(50)
      EQUIVALENCE (T1Y(1),TT11(1)),(T1Y(66),TT12(1))
      EQUIVALENCE (T2XY(1),TT21(1)),(T2XY(51),TT22(1)),(T2XY(101),
     *TT23(1)),(T2XY(151),TT24(1))
      DATA GAM2/2.,2.,6.,6.,6.,4.,20.,30.,20.,54.,73.,101.,73.,8.,
     *146.,100.,100.,343.,174.,393.,186.,61.,202.,113.,213.,233.,180.,
     *696.,194.,120.,458.,590.,0.,2790.,216./
      DATA GAM1/2.,2.,6.,6.,6.,4.,20., 6.,20.,18.,28., 40.,28.,8.,
     * 36., 10., 36., 70., 44., 36., 44.,12., 26., 14., 26., 42., 30.,
     * 80., 30., 16.,102., 18.,0.,54.,216./
      DATA GAM/2.,2.,6.,6.,6.,4.,20., 6.,20.,18.,28., 40.,28.,8.,
     * 36., 10., 36., 70., 44., 36., 44.,12., 26., 14., 26., 42., 30.,
     * 80., 30., 16.,102., 18.,0.,54.,216./
      DATA AJ/1.,1.,2.,3.,3.,4.,5.,6.,5.,6.,7.,8.,7.,8.,9.,10.,9.,
     *10.,11.,12.,11.,12.,13.,14.,13.,14.,15.,16.,15.,16.,17.,18.,0.,
     *18.,24./
      DATA ZJ/0.,1.,1.,1.,2.,2.,2.,2.,3.,3.,3.,3.,4.,4.,4.,4.,5.,5.,
     *5.,5.,6.,6.,6.,6.,7.,7.,7.,7.,8.,8.,8.,8.,0.,9.,11./
      DATA RNCL/0.0/,AMS/0.125/,AMFS/0.125/,ICAN/1/,ISEC/1/
      DATA ENLEV/700*0./,SPLEV/700*0./
      DATA TT11/
     *  20.80,  15.80,  21.00,  16.80,  19.80,
     *  16.50,  18.80,  16.50,  18.50,  17.20,
     *  18.26,  15.05,  16.01,  12.04,  13.27,
     *  11.09,  12.17,  10.26,  11.04,   8.41,
     *   9.79,   7.36,   8.15,   5.63,   5.88,
     *   3.17,   3.32,   0.82,   1.83,   0.97,
     *   2.33,   1.27,   2.92,   1.61,   2.91,
     *   1.35,   2.40,   0.89,   1.74,   0.36,
     *   0.95,  -0.65,  -0.04,  -1.73,  -0.96,
     *  -2.87,  -2.05,  -4.05,  -3.40,  -5.72,
     *  -3.75,  -4.13,  -2.42,  -2.85,  -1.01,
     *  -1.33,   0.54,  -0.02,   1.74,   0.75,
     *   2.24,   1.00,   1.98,   0.79,   1.54/
      DATA TT12/
     *   0.39,   1.08,   0.00,   0.78,  -0.35,
     *   0.58,  -0.55,   0.59,  -0.61,   0.59,
     *  -0.35,   0.32,  -0.96,  -0.52,  -2.08,
     *  -2.46,  -3.64,  -1.55,  -0.96,   0.97,
     *   0.88,   2.37,   1.75,   2.72,   1.90,
     *   2.55,   1.46,   1.93,   0.86,   1.17,
     *   0.08,   0.39,  -0.76,  -0.39,  -1.51,
     *  -1.17,  -2.36,  -1.95,  -3.06,  -2.62,
     *  -3.55,  -2.95,  -3.75,  -3.07,  -3.79,
     *  -3.06,  -3.77,  -3.05,  -3.78,  -3.12,
     *  -3.90,  -3.35,  -4.24,  -3.86,  -4.92,
     *  -5.06,  -6.77,  -7.41,  -9.18, -10.16,
     * -11.12,  -9.76,  -9.23,  -7.96,  -7.65/
      DATA TT21/
     *  -8.40, -12.90,  -8.00, -11.90,  -9.20,
     * -12.50, -10.80, -13.60, -11.20, -12.20,
     * -12.81, -15.40, -13.07, -15.80, -13.81,
     * -14.98, -12.63, -13.76, -11.37, -12.38,
     *  -9.23,  -9.65,  -7.64,  -9.17,  -8.05,
     *  -9.72,  -8.87, -10.76,  -8.64,  -8.89,
     *  -6.60,  -7.13,  -4.77,  -5.33,  -3.06,
     *  -3.79,  -1.72,  -2.79,  -0.93,  -2.19,
     *  -0.52,  -1.90,  -0.45,  -2.20,  -1.22,
     *  -3.07,  -2.42,  -4.37,  -3.94,  -6.08/
      DATA TT22/
     *  -4.49,  -4.50,  -3.14,  -2.93,  -1.04,
     *  -1.36,   0.69,   0.21,   2.11,   1.33,
     *   3.29,   2.46,   4.30,   3.32,   4.79,
     *   3.62,   4.97,   3.64,   4.63,   3.07,
     *   4.06,   2.49,   3.30,   1.46,   2.06,
     *   0.51,   0.74,  -1.18,  -1.26,  -3.54,
     *  -3.97,  -5.26,  -4.18,  -3.71,  -2.10,
     *  -1.70,  -0.08,  -0.18,   0.94,   0.27,
     *   1.13,   0.08,   0.91,  -0.31,   0.49,
     *  -0.78,   0.08,  -1.15,  -0.23,  -1.41/
      DATA TT23/
     *  -0.42,  -1.55,  -0.55,  -1.66,  -0.66,
     *  -1.73,  -0.75,  -1.74,  -0.78,  -1.69,
     *  -0.78,  -1.60,  -0.75,  -1.46,  -0.67,
     *  -1.26,  -0.51,  -1.04,  -0.53,  -1.84,
     *  -2.42,  -4.52,  -4.76,  -6.33,  -6.76,
     *  -7.81,  -5.80,  -5.37,  -3.63,  -3.35,
     *  -1.75,  -1.88,  -0.61,  -0.90,   0.09,
     *  -0.32,   0.55,  -0.13,   0.70,  -0.06,
     *   0.49,  -0.20,   0.40,  -0.22,   0.36,
     *  -0.09,   0.58,   0.12,   0.75,   0.15/
      DATA TT24/
     *   0.70,   0.17,   1.11,   0.89,   1.85,
     *   1.62,   2.54,   2.29,   3.20,   2.91,
     *   3.84,   3.53,   4.48,   4.15,   5.12,
     *   4.78,   5.75,   5.39,   6.31,   5.91,
     *   6.87,   6.33,   7.13,   6.61,   7.30,
     *   6.31,   6.27,   4.83,   4.49,   2.85,
     *   2.32,   0.58,  -0.11,  -0.98,   0.81,
     *   1.77,   3.37,   4.13,   5.60,   6.15,
     *   7.29,   7.35,   7.95,   7.67,   8.16,
     *   7.83,   8.31,   8.01,   8.53,   8.27/
       END
      
      
      
      FUNCTION BARIE1(A,Z)
C  B[[[C[EH[E [AP[EPA [E[EH[[ [[PA A,Z (M[B)
C  ( BARASHENKOV,ILJINOV,TONEEV -1972)
      COMMON /BL1001/T1Y(130) /BL1002/T2XY(200)
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      I=INT(A+0.1)
      J=INT(Z+0.1)
      L=I-J
      X=Z*Z/A
      IF(X.LE.33.5) BF0=12.5+4.7*(33.5-X)**0.75
      IF(X.GT.33.5) BF0=12.5-2.7*(X-33.5)**0.666667
      IF((2*(J/2)).LT.J) D=0.
      IF((2*(J/2)).EQ.J) D=-0.5
      IF((2*(L/2)).EQ.L) D=D
      IF((2*(L/2)).LT.L) D=D+1.
      IF(J.LE.0.OR.J.GT.130.OR.L.LE.0.OR.L.GT.200) GO TO 22
      T1=T1Y(J)
      T2=T2XY(L)
      GO TO 23
 22   T1=0.
      T2=0.
 23   CONTINUE
      BARIE1=BF0+D-T1-T2
      IF(BARIE1.LE.0.) BARIE1=0.
      RETURN
      END
      
      
      
      SUBROUTINE LYMASM(IZ,IA,CMASS,CBARR,NOBARR)
C
C     WILLIAM D. MYERS - 6 JULY 1970
C
      COMMON /BLOCX/ IPARQ
      COMMON /FFSS/  ENEX
      COMMON /FUSR/   BARR,SMASS,SHLL
C     COMMON /BLOCXX/ ARQ,F1Q,F1MQ,SUFNUC
C     COMMON /CCC/    WOTNUC,VOLNUC,COULMB,
C    *                A,Z,UN,A1,A2,A3,GGMMA,A3RT2,A3RT,ZSQ,
C    *                ODDEV,SYM,PARMAS,ACOR
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      DIMENSION EM(10),EMP(10),XK(10),Y(2),F(2)
      DATA ZVT/1.6666666666/, ZT/.3333333333/,
     1 ZTT/.6666666667/, SR5/2.2360679775/
      DATA EM/0.00,2.00,8.00,14.00,28.00,50.00,82.00,126.00,
     1        184.00,258.00/
      DATA CAY1/0.0/, CAY2/0.0/, CAY3/2.0/, CAY4/11.0/, CAY5/8.07144/,
     * CAY6/7.28899/,
     *            D/.444/, C/5.8/, SMALC/.325/
C
C     DMASS = REMAINDER AFTER CM - NO SHELL EFFECTS SUBTRACTED
C     SHLL = CALCULATED SHELL EFFECT
C     DIFMAS= DMASS - SHLL
C
C------------------------------
      A1 = 15.4941
      IPARQ=0
C..... IPARQ=0  PARAMETRS MYERS-SWIATECKI
      IF(IPARQ.NE.0) GO TO 121
      A2 = 17.9439
      A3 = 0.7053
      GAMMA = 1.7826
      GO TO 126
 121  CONTINUE
C..... IPARQ=1  PARAMETRS KRAPPE-NIX
      IF(IPARQ.NE.1) GO TO 123
      A2 = 24.70
      A3 = 0.74476032
      GAMMA = 4.0
      GO TO 126
 123  CONTINUE
C..... IPARQ=2  PARAMETRS PAULI-LEDERGERBER
      IF(IPARQ.NE.2) GO TO 124
      A2 = 19.008
      A3 = 0.720
      GAMMA = 2.840
      GO TO 126
 124  CONTINUE
C..... IPARQ=3  BF(T)   PARAMETRS MYERS-SWIATECKI
      ALEVEL=0.1
      AMPAR=ALEVEL*FLOAT(IA)
      TSQ=ENEX/AMPAR
      A2 = 17.9439*(1.-0.0063157*TSQ)
      A3 = 0.7053*(1.-0.001*TSQ)
      GAMMA = 1.7826
 126  CONTINUE
C------------------------------
      IF(IZ.NE.0) GO TO 15
      CMASS=0.0
      RETURN
   15 NOBARR=0
      DO 1 I=1,10
      EMP(I)=EM(I)**ZVT
    1 CONTINUE
      DO 2 I=1,9
      XK(I)=.600*(EMP(I+1)-EMP(I)) /(EM(I+1)-EM(I))
    2 CONTINUE
C
C     FOR DEFINITIONS OF CAY1 AND RZ,SEE UCRL-11980
C
      CAY1=3.28637900*A3**3
      RZ=.86398700/A3
    5 Z= FLOAT(IZ)
      ZSQ=Z**2
      N=IA-IZ
      UN= FLOAT(N)
      A= FLOAT(IA)
      A3RT=A**ZT
      A3RT2=A3RT**2
      A2RT= SQRT(A)
      SYM=((UN-Z)/A)**2
      ACOR=1.00-GAMMA*SYM
      PARMAS=CAY5*UN+CAY6*Z
      VOLNUC=-A1*ACOR*A
      SUFNUC=A2*ACOR*A3RT2
      COULMB=A3*ZSQ/A3RT
      FUZSUR=-CAY1*ZSQ/A
      ODDEV=-(1.00+2.00*(N/2)-UN+2.00*(IZ/2)-Z)/A2RT*CAY4
      WTERM=-CAY2*A3RT2* EXP(-CAY3*SYM)
      WOTNUC=PARMAS+COULMB+FUZSUR+ODDEV+WTERM
      SMASS=WOTNUC+VOLNUC+SUFNUC
      SPW=SUFNUC+WTERM
      C2=SPW/A3RT2
      X=.5*COULMB/SPW
      IF(X.GE.1.00) GO TO 4
C------------------------------
      ARQ=X
      IF(IPARQ.EQ.0) BARR=SUFNUC*XI(X)
      IF(IPARQ.EQ.2) BARR=SUFNUC*XI(X)
      IF(IPARQ.EQ.3) BARR=SUFNUC*XI(X)
      IF(IPARQ.EQ.1) BARR=SUFNUC*XIMOD(X)
      IF(IPARQ.EQ.0) F1Q=XI(X)
      IF(IPARQ.EQ.1) F1MQ=XIMOD(X)
C------------------------------
      GO TO 6
    4 BARR=0.0
    6 Y(1)=UN
      Y(2)=Z
      DO 31 J=1,2
      DO 32 I=1,9
      IF(Y(J)-EM(I+1)) 3,3,32
   32 CONTINUE
      PRINT 332,J
  332 FORMAT('1FAILURE IN LYMASS - Y(',I1,') EXCEEDS LAST MAGIC NO.')
      STOP
    3 F(J)=XK(I)*(Y(J)-EM(I))-.600*(Y(J)**ZVT-EMP(I))
   31 CONTINUE
      S=(2.00/A)**ZTT*(F(1)+F(2))-SMALC*A3RT
      C2D2=C2*D**2
      EE=(C2D2+C2D2)*(1.00-X)
      FF=.425917710*C2D2*D*(1.00+X+X)/A3RT
      SSHELL=C*S
      V=SSHELL/EE
      EPS=1.500*FF/EE
      IF(EE*(1.00-3.00*V).LE.0.00) GO TO 51
      QCALC=0.00
      THETA=0.00
      SHLL=SSHELL
      GO TO 52
C
C       ESTIMATE THETA
C
   51 TO=1.00
C
C       ITERATE TO FIND EQUILIBRIUM THETA
C
  101 DO 725 IPQ=1,10
      TO2=TO**2
C----------------------------------------
C     IF (TO2.GT.170.) PRINT 500, IZ,IA
C 500 FORMAT(1X,'LYMASM',2X,2I5)
C----------------------------------------
C     EXMT2= EXP(-TO2)
      EXMT2= 1.E-20
      IF((ABS(TO2)).LT.30.) EXMT2= EXP(-TO2)
C
      T=TO-(1.00-EPS*TO-V*(3.00-TO2-TO2)*EXMT2) /
     1(-EPS+V*TO*(10.00-4.00*TO2)*EXMT2)
      IF(T.LE.0.00) GO TO 728
      IF( ABS(T-TO) .LT.1.E-4) GO TO 732
      TO=T
  725 CONTINUE
      GO TO 729
  732 T2=T**2
C     EXT2= EXP(-T2)
      EXT2= 1.E-20
      IF((ABS(T2)).LT.30.) EXT2= EXP(-T2)
C
      TEST=EE*(1.00-EPS*(T+T)-V*((4.00*T2-12.00)*T2+3.00)* EXT2)
      IF(TEST.GT.0.00) GO TO 81
  728 TO=.100
      DO 100 I=1,20
      TO2=TO**2
      GL=EE*(1.00-EPS*TO-V*(3.00-TO2-TO2)* EXP(-TO2))
      IF(GL.GT.0.00) GO TO 101
  100 CONTINUE
  729 CMASS=SMASS
      CBARR=0.00
      NOBARR=1
      RETURN
   81 THETA=T
      ALPHA0=D*SR5/A3RT
      ALPHA=ALPHA0*THETA
      SIGMA=ALPHA*(1.00+ALPHA/14.00)
      EXS= EXP(SIGMA+SIGMA)- EXP(-SIGMA)
      QCALC=4.E-3*Z*(RZ*A3RT)**2*EXS
      T2=T**2
      SHLL=T2*(EE-FF*T) + SSHELL*(1.00-T2-T2)* EXP(-T2)
   52 CMASS=SMASS+SHLL
      CBARR=BARR-SHLL
      RETURN
      END
      
      
      
      FUNCTION XI(Z)
C
C     6-POINT LAGRANGE INTERPOLATION
C
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      DIMENSION Y(51)
      DATA Y/.25900,.255200,.250700,.245100,.2400,.23400,.228500,
     1      .22200,.21600,.2100,.20300,.196800,.1900,.18300,.175800,
     2      .1692400,.1620300,.1547800,.147500,.1401900,.1328400,
     3      .1254500,.1180100,.1105200,.1029600,.0953500,.0876800,
     4      .0799900,.0722900,.064600,.0569500,.0493700,.0419300,
     5      .0347600,.0281100,.0223600,.0176200,.0137300,.0105600,
     6      .0079800,.0059100,.0042500,.0029600,.0019700,.0012300,
     7      7.1E-4,3.6E-4,1.5E-4,4.E-5,1.E-5,0.00/
C
C     THE X VALUES ARE EVENLY SPACED - X = 0(.02)1
C
      ZBH=Z*50.00
      M=IFIX(ZBH)
      DEL=ZBH- FLOAT(M)
      M=M+1
      IF(M.LE.51) GO TO 105
      M=51
  100 XI=Y(M)
      RETURN
  105 IF (DEL.LT.1.E-4) GO TO 100
      IF(M.GE.3) GO TO 110
      DEL=DEL- FLOAT(3-M)
      M=3
      GO TO 115
  110 IF(M.LE.48) GO TO 115
      DEL=DEL+ FLOAT(M-48)
      M=48
  115 DM3=DEL-3.00
      PROD=DM3*DEL
      W6=1.00/(1.2E2*DM3)
      DM2=DM3+1.00
      PROD=DM2*PROD
       W5=-1.00/(24.00*DM2)
      DM1=DM2+1.00
      PROD=DM1*PROD
      W4=1.00/(12.00*DM1)
      DP1=DM1+2.00
      PROD=DP1*PROD
      W2=1.00/(24.00*DP1)
      DP2=DP1+1.00
      PROD=DP2*PROD
      W1=-1.00/(1.2E2*DP2)
      W3=-1.00/(12.00*DEL)
      XI=PROD*(W1*Y(M-2)+W2*Y(M-1)+W3*Y(M)+W4*Y(M+1)+W5*Y(M+2)
     1 +W6*Y(M+3))
      RETURN
      END
      
      
      
      FUNCTION XIMOD(Z)
C
C     6-POINT LAGRANGE INTERPOLATION
C      IN MODIFIED LIQUID-DROP FORMULA
C         ( KRAPPE [ NIX --  IAEA-SM-174/12 )
C
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      DIMENSION Y(51)
      DATA  Y/
     1    0.12200, 0.12100, 0.11980, 0.11830, 0.11690, 0.11520, 0.1133,
     2    0.11130, 0.10900, 0.10670, 0.10420, 0.10150, 0.09850, 0.09540,
     3    0.09180, 0.08780, 0.08350, 0.07900, 0.07460, 0.06960, 0.06470,
     4    0.05960, 0.05420, 0.04880, 0.04350, 0.03880, 0.03400, 0.02920,
     5    0.02460, 0.02020, 0.01580, 0.01220, 0.00900, 0.00660, 0.00490,
     6    0.00360, 0.00280, 0.00220, 0.00180, 0.00140, 0.00100, 0.00090,
     7    0.00060, 0.00040, 0.00020, 0.00010, 0.00000, 0.00000, 0.00000,
     8    0.00000, 0.00000/
C
C     THE X VALUES ARE EVENLY SPACED - X = 0(.02)1
C
      ZBH=Z*50.00
      M=IFIX(ZBH)
      DEL=ZBH- FLOAT(M)
      M=M+1
      IF(M.LE.51) GO TO 105
      M=51
  100 XI=Y(M)
      RETURN
  105 IF (DEL.LT.1.E-4) GO TO 100
      IF(M.GE.3) GO TO 110
      DEL=DEL- FLOAT(3-M)
      M=3
      GO TO 115
  110 IF(M.LE.48) GO TO 115
      DEL=DEL+ FLOAT(M-48)
      M=48
  115 DM3=DEL-3.00
      PROD=DM3*DEL
      W6=1.00/(1.2E2*DM3)
      DM2=DM3+1.00
      PROD=DM2*PROD
       W5=-1.00/(24.00*DM2)
      DM1=DM2+1.00
      PROD=DM1*PROD
      W4=1.00/(12.00*DM1)
      DP1=DM1+2.00
      PROD=DP1*PROD
      W2=1.00/(24.00*DP1)
      DP2=DP1+1.00
      PROD=DP2*PROD
      W1=-1.00/(1.2E2*DP2)
      W3=-1.00/(12.00*DEL)
      XIMOD=PROD*(W1*Y(M-2)+W2*Y(M-1)+W3*Y(M)+W4*Y(M+1)+W5*Y(M+2)
     1 +W6*Y(M+3))
      RETURN
      END
