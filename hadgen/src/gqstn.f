
************************************************************************
*                                                                      *
*  The following subroutines are written by K.K.Gudima and V.D.Toneev  *
*  Joint Institute for Nuclear Research, Dubna,                        *
*  E-mail: gudima@acad.moldnet.md                                      *
*          gudima@cc.acad.md                                           *
*  E-mail: toneev@thsun1.jinr.ru                                       *
*                                                                      *
************************************************************************

C     VERSION FROM  27-06-95 01:04pm
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE  STRANT(BU,BV,NU,N1,N2,BPIN,IIN,BPN,IPN,NREP,JPSTR,
     *SIGT,SABS)
      IMPLICIT REAL *8(A-H,O-Z)
      REAL *8 BU,BV(3),BPIN(9),BPN(9),SIGT,SABS
      COMMON/NRAPI/NRAPI
      COMMON/NUCSP/VPR(3),VTA(3),RADP(3),RADT(3),VEV(3),VRE(3),GEV,GRE,
     *VEP(3),VET(3),GEP,GET
      COMMON/INSTR/ INSTR
      COMMON /AMPROJ/ AMPROJ
      COMMON/DEUT/VDN(3),UDN,T1DN,WDN,VND(3),UND,T1ND,WND
      COMMON/FLUC/PDN(9),PND(9)
      COMMON /IDN12/ ID1,ID2
      DIMENSION  V(3),IIN(5),PN(9),IPN(5),VN(3),RL(3),RN(3),PIN(9),
     *PN1(9),PN2(9),IPN1(5),IPN2(5),PS1(3),PS2(3),PL1(3),PL2(3)
C  ******************************************************************
C    JPS  ISRES=0    ISRES=1     INPT                     Y(SP)STR
C                                                         BESIDES 3,4
C      1   K+         ETA        1 - E*DS/D3P (TETAL,PL) 1- M+B
C      2   K0         ETAP       2 - DS/D3P   (TETA*,T*) 2- PI+N
C      3   LAMBDA     RHO0       3 - DS/D/DT  (TETAL,TL) 3- B+B
C      4   SIGMA+     OMEGA      4 - DS/DP/DO (TETAL,PL) 4- N+N
C      5   SIGMA-     PHI        5 -E*DS/D3P  (TETAL,TL) 5- N0+N0
C      6   SIGMA0     (AP)       6 -EdS/dPTdY (PT,Y)     6- M+D
C      7   K-                                            7- D+B
C      8   AP                                            8- B+D
C      9   PI-
C     10   PI0
C     11   PI+
C                     INSTR=1     APPROXIMATED EXP. DATA (WITHOUT FLUCTONS)
C                     INSTR=2     APPROXIMATED EXP. DATA (NEW FERMI MOM.)
C                     INSTR=3     AMELIN'S ACT
C                     INSTR=4     AMELIN'S ACT (BOLTZ. SPECTRA)
C                     INSTR=5     APPROX.WITHOUT FLUCTONS
C                     INSTR=6    *APPR. EXP.DATA(BOLTZ. SPECTRA)
C  ******************************************************************
C     KSIGS=1    B+B and M+B prod. cross section = our approximations
C     KSIGS=2    B+B = appr. ZWERMANN (P.L.183B,p.31,1987),M+B=ours
C     KSIGS=3    B+B = appr. ZWERMANN; M+B=appr. CUGNON (N.P.A422,p.635,1984)
C     KSIGS=4    B+B = appr. ZWERMANN; M+B=appr. CUGNON (P.R.C41,P.1701,1990)
C     KSIGS=5    B+B = appr. ZWERMANN; M+B=appr. Sibirtsev(Z.F.A347,p.191,1994)
C  ******************************************************************
C     FOR THE P+A CASE AND INSTR=6 THE FIRST-CHANCE PRODUCTION IS
C     IN THE Y(SP)STR  ARRAY  WITH INDEX=3
C  ******************************************************************
C     IWSTR=1 W=SK/STOT,  DSIGMA=SUM(W)/NCAS*SINEL
C     IWSTR=2 W=SK     ,  DSIGMA=SUM(W)/NCAS
C  ******************************************************************
C                      KSIGS=1, direct pions productions (our.approx.)
C     FOR JPS=9,10,11
C                      KSIGS > 1, N+N=>D+N; D=PI+N
C  ******************************************************************
C      write( *,*) 'in strant'
      SIG=SIGT
      SAB=SABS
      WDN=0.
      WND=0.
      U=BU
      DO  5  K=1,9
      PDN(K)=0.
      PND(K)=0.
      IF(K.LE.3)  V(K)=BV(K)
      PIN(K)=BPIN(K)
    5 PN(K)=BPN(K)
      IB1=IIN(4)
      IB2=IPN(4)
      IE1=IIN(1)
      IE2=IPN(1)
      AMPROJ=PIN(9)
      SVR=0.
      DO   9  K=1,3
      IF(NRAPI.EQ.2.OR.NRAPI.EQ.3)   GO  TO  6
      VN(K)=0.
      RN(K)=0.
      GO  TO  8
    6 IF(NU.EQ.2)  GO  TO  7
      VN(K)=VPR(K)
      RN(K)=RADP(K)
      GO  TO  8
    7 VN(K)=VTA(K)
      RN(K)=RADT(K)
    8 SVR=SVR+VN(K)*(PIN(K)+PN(K))/2.
    9 CONTINUE
      G=1./DSQRT(1.-VN(1)**2-VN(2)**2-VN(3)**2)
      GG=SVR*G/(G+1.)
      DO  10  K=1,3
      RL(K)=(PIN(K)+PN(K))/2.-GG*VN(K)+RN(K)
   10 CONTINUE
C
      DO  13 KR=1,NREP
      IF(NRAPI.EQ.5)       GO  TO  11
c  redefinition of Fermi momenta
      IF(NRAPI.EQ.1)  THEN
       CALL  PARTN(1,N1,PN1,IPN1)
       CALL  PARTN(2,N2,PN2,IPN2)
       PN12=PN1(4)**2+PN1(5)**2+PN1(6)**2
       PN22=PN2(4)**2+PN2(5)**2+PN2(6)**2
       IF(PN12.LE.0.OR.PN22.LE.0.)  GO  TO  13
       C1=DSQRT((PN1(8)*(PN1(8)+2.*PN1(9)))/PN12)
       C2=DSQRT((PN2(8)*(PN2(8)+2.*PN2(9)))/PN22)
       do  k=1,3
       PS1(k)=PN1(3+k)*C1
       PS2(k)=PN2(3+k)*C2
       enddo
       CALL  KINEMA(PS1,VPR,PL1,CT1,ST1,CF1,SF1,TL1,PN1(9))
       CALL  KINEMA(PS2,VTA,PL2,CT2,ST2,CF2,SF2,TL2,PN2(9))
       do  k=1,3
       PIN(3+k)=PL1(k)
       PN(3+k) =PL2(k)
       enddo
       PIN(8)=DSQRT(PL1(1)**2+PL1(2)**2+PL1(3)**2+0.940**2)-0.940 
       PN(8) =DSQRT(PL2(1)**2+PL2(2)**2+PL2(3)**2+0.940**2)-0.940 
       PIN(9)=0.940
       PN(9) =0.940
      ELSE
       CALL  PARTN(NU,N2,PN2,IPN2)
       do  k=1,3
       PN(3+k)=PN2(3+k)
       enddo
       PN(8)=PN2(8)
       PN(9)=0.94
      ENDIF
      PX1=PIN(4)
      PY1=PIN(5)
      PZ1=PIN(6)
      AM1=PIN(9)
      PX2=PN(4)
      PY2=PN(5)
      PZ2=PN(6)
      AM2=PN(9)
      CALL CROSEC(1,ID1,ID2,PX1,PY1,PZ1,AM1,PX2,PY2,PZ2,AM2,SITO,0)
      SIG=SITO
   11 CONTINUE
      CALL  TINVU(PIN,PN,U,V,T12)
      IF(INSTR.NE.6.OR.NRAPI.EQ.5)  GO  TO  12
      CALL  NFERMI(PIN,PN,NU,N2)
      IF(NRAPI.EQ.1.OR.NRAPI.EQ.3) CALL TINVU(PIN,PND,UND,VND,T1ND)
      IF(NRAPI.EQ.1)               CALL TINVU(PDN, PN,UDN,VDN,T1DN)
      IF(NRAPI.EQ.2)               CALL TINVU(PIN,PND,UND,VND,T1ND)
   12 CONTINUE
      STOT=SIG+SAB
C     WRITE(16,*) 'CONEC in'
      CALL  CONEC(U,T12,V,VN,PIN(9),PN(9),IIN,IPN,RL,NREP,JPSTR,STOT)
C     WRITE( *,*) 'CONEC out'
   13 CONTINUE
C     write( *,*) 'out strant'
      RETURN          
      END
C
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
      SUBROUTINE  POSFL(IB1,IE1,IE2,JP,B,STOT,RL,N)
      IMPLICIT REAL *8(A-H,O-Z)
      REAL *8 B,STOT,RL,MPI,MN
      COMMON/HCASC/ANUCL1,ANUCL2,ZNUCL1,ZNUCL2,T0,EPS1,EPS2,VPI,A1,A2,
     *C1,C2,D1,D2,R0N1,R0N2,TF01,TF02,RM1,RM2
      COMMON/DEUT/VDN(3),UDN,T1DN,WDN,VND(3),UND,T1ND,WND
      COMMON/FLUC/PDN(9),PND(9)
      COMMON/NRAPI/NRAPI
      COMMON /IDN12/ ID1,ID2
      COMMON /AMPROJ/ AMPROJ
      COMMON /NUCOLL/ NUCOLL,MVCOLL(5999)
      COMMON /KSIGS/ KSIGS
      DIMENSION  B(3),RL(3),AM(11),VD(3)
      DATA  AM/2*0.494,1.115,3*1.192,0.494,0.939,3*0.139/
      DATA MN/0.939/,MPI/0.139/
      IF(NRAPI.GT.3)  RETURN
      AMJ=AM(JP)
      GO  TO  (10,10,11,11,11,11,12,121,122,122,122),JP
   10 EMN=AMJ+1.115+IB1*.940
      GO  TO  13
   11 EMN=AMJ+0.494+IB1*.940
      GO  TO  13
   12 EMN=AMJ+0.494+IB1*.940+.940
      GO  TO  13
  121 EMN=AMJ+0.939+IB1*.940+.940
      GO  TO  13
  122 EMN=AMJ+0.939+IB1*.940
      GO  TO  13
   13 DO  22  KD=1,2
      IF(NRAPI.GT.1.AND.KD.EQ.1)  GO  TO  22
      IF(KD-1) 14,14,16
   14 UD=UDN
      IB=1
      WD=WDN
      TD=T1DN
      IE=IE2
      ID=ID2
      DO  15  K=1,3
   15 VD(K)=VDN(K)
      GO  TO  18
   16 UD=UND
      IB=IB1
      WD=WND
      TD=T1ND
      IE=IE1
      ID=ID1
      DO  17  K=1,3
   17 VD(K)=VND(K)
   18 IF(UD.LE.EMN)   GO  TO  22
      IF(JP.GE.9)  GO  TO  192
      IF(JP.EQ.8)  GO  TO  93
      IF(JP.EQ.7)  GO  TO  19
C
      ST=0.
      IF(KSIGS.GE.2.AND.IB.EQ.1) THEN
      ST=SZWER(JP,UD)
      GO  TO  201
      ENDIF
      IF((KSIGS.EQ.3.OR.KSIGS.EQ.4).AND.IB.EQ.0) THEN
      ST=(SCT( 1,1,JP,TD,1)+SCT( 1,0,JP,TD,1)+
     +    SCT(-1,1,JP,TD,1)+SCT(-1,0,JP,TD,1)+
     +    SCT( 0,1,JP,TD,1)+SCT( 0,0,JP,TD,1))/6.
      GO  TO  201
      ENDIF
      IF(KSIGS.EQ.5.AND.IB.EQ.0) THEN
      ST=0.
      SS0=0.494+1.115
      IF(UD.GT.SS0.AND.UD.LT.1.7)   ST=10.*(UD-SS0)
      IF(UD.GT.1.7)                 ST=0.09/(UD-1.6)
      GO  TO  201
      ENDIF
C
      IF(IB.NE.0)  GO  TO  91
      IF(ID.EQ.120.OR.ID.EQ.-120.OR.ID.EQ.110) GO  TO  92
C
C MESON+NUCLEON PRODUCTION CROSS SECTION
      ST=(SST(0,1, 1,1,JP,TD,1)+SST(0,1, 1,0,JP,TD,1)+
     +    SST(0,1,-1,1,JP,TD,1)+SST(0,1,-1,0,JP,TD,1)+
     +    SST(0,1, 0,1,JP,TD,1)+SST(0,1, 0,0,JP,TD,1))/6.
      GO  TO  201
   91 CONTINUE
      IF(IB.NE.1)   RETURN
      IF(ID.EQ.1120.OR.ID.EQ.1220) GO  TO 92
C
C BARYON+NUCLEON PRODUCTION CROSS SECTION
      ST=(SST(1,1,1,1,JP,TD,1)+SST(1,1,1,0,JP,TD,1)+
     +    SST(1,1,0,1,JP,TD,1)+SST(1,1,0,0,JP,TD,1))/4.
      GO  TO  201
   92 CONTINUE
C NUCLEON+NUCLEON AND PION+NUCLEON PRODUCTION CROSS SECTION
      ST=0.25*(SST(IB,1,IE,1,JP,TD,1)+SSDT(IB,1,IE,1,JP,TD,1)+
     +         SST(IB,1,IE,0,JP,TD,1)+SSDT(IB,1,IE,0,JP,TD,1))
      S4=0.25*(SST(IB,1,IE,1,JP,4.D0,1)+SSDT(IB,1,IE,1,JP,4.D0,1)+
     +         SST(IB,1,IE,0,JP,4.D0,1)+SSDT(IB,1,IE,0,JP,4.D0,1))
      IF(TD.GE.4.)  ST=STHE(UD,JP,S4)
      GO  TO  201
   93 CONTINUE
C  ******  CHANNELS FOR  ANTI-P  PRODUCTION
      IF(UD.LE.EMN)  RETURN
      EMN1=EMN-AMJ
      PMAX=DSQRT(((UD+AMJ)**2-EMN1**2)*((UD-AMJ)**2-EMN1**2))/(2.*UD)
C !!!!
C  A. SHOR NP A514(1990) p.717;  DANILIEWICZ PR C42(1990) p.1564
C  SEE: PL B256(1991) p.331
      ST=0.01*(UD-EMN)**1.846
C !!!!
      GO  TO  201
   19 CONTINUE
C
C  ******  CHANNELS FOR  K-   PRODUCTION  G: 20.05.93
      IF(UD.LE.EMN)  RETURN
      EMN1=EMN-AMJ
      PMAX=DSQRT(((UD+AMJ)**2-EMN1**2)*((UD-AMJ)**2-EMN1**2))/(2.*UD)
      IF(IB.NE.0)    GO  TO  191
C  M+FL => K-
      ST=PMAX/13.
      GO  TO  201
  191 CONTINUE
C  B+FL => K-
      ST=SIKMI(PMAX)
      IF(TD.LT.4.)    GO  TO  201
      U4=1.88*(4.+1.88)
      PMA4=DSQRT(((U4+AMJ)**2-EMN1**2)*((U4-AMJ)**2-EMN1**2))/(2.*U4)
      S4=SIKMI(PMA4)
      ST=STHE(UD,JP,S4)
C
      GO  TO  201
  192 CONTINUE
C  ******  CHANNELS FOR  PI-0+   PRODUCTION
      IF(IB.EQ.1)  GO  TO  193
      S1= CROSEG(0,0,1,1,7,TD,MPI,0)
      S2= CROSEG(0,0,1,2,7,TD,MPI,0)
      S3= CROSEG(0,0,1,3,7,TD,MPI,0)
      ST=(S1+S2+S3)/3.
      GO  TO  201
  193 CONTINUE
      S11= CROSEG(0,0,2,1,7,TD,MN,0)
      S12= CROSEG(0,0,2,2,7,TD,MN,0)
      ST=(S11+S12)/2.
  201 CONTINUE
      SK=ST*WD
      IF(SK.LE.0.)  GO  TO  22
c     WRITE ( *,999) UD,EMN,TD,ST,WD
 999  FORMAT(1X,'UD,EMN,TD,ST,WD=',3(F6.3,1X),2(1X,E11.4))
C
C
      EM1=EMN-AMJ+RNDMD(-1)*(UD-EMN)
      E1=(UD**2+AMJ**2-EM1**2)/(2.*UD)
      P1=DSQRT(E1**2-AMJ**2)
      CT1=1.-2.*RNDMD(-1)
      IF(IB.EQ.0)  IND=6
      IF(IB.NE.0)  IND=6+KD
      CALL  POSTR(P1,CT1,AMJ,VD,B,JP,SK,STOT,RL,N,IND)
   21 CONTINUE
   22 CONTINUE
      RETURN
      END
C
C * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
      SUBROUTINE  STRANN(BU,BT12,BV,NU,N1,N2,
     *PIN,IIN,PN,IPN,NREP,JPSTR,SIGT,SABS)
      IMPLICIT REAL *8(A-H,O-Z)
      REAL *8 BU,BT12,BV,PIN,PN,SIGT,SABS
      COMMON/NRAPI/NRAPI
      COMMON/NUCSP/VPR(3),VTA(3),RADP(3),RADT(3),VEV(3),VRE(3),GEV,GRE,
     *VEP(3),VET(3),GEP,GET
      COMMON/INSTR/ INSTR
      COMMON /AMPROJ/ AMPROJ
      DIMENSION  BV(3),V(3),PIN(9),IIN(5),PN(9),IPN(5),
     *VN(3),RL(3),RN(3),
     *PN1(9),PN2(9),IPN1(5),IPN2(5),PS1(3),PS2(3),PL1(3),PL2(3)
      SIG=SIGT
      SAB=SABS
      U=BU
      AMPROJ=PIN(9)
      SVR=0.
      Nw1=N1
      DO  10  K=1,3
      V(K)=BV(K)
      IF(NRAPI.EQ.2.OR.NRAPI.EQ.3)   GO  TO  7
      VN(K)=0.
      RN(K)=0.
      GO  TO  9
    7 IF(NU.EQ.2)  GO  TO  8
      VN(K)=VPR(K)
      RN(K)=RADP(K)
      GO  TO  9
    8 VN(K)=VTA(K)
      RN(K)=RADT(K)
    9 SVR=SVR+VN(K)*PIN(K)
   10 CONTINUE
      G=1./DSQRT(1.-VN(1)**2-VN(2)**2-VN(3)**2)
      GG=SVR*G/(G+1.)
      RL(1)=PIN(1)-GG*VN(1)+RN(1)
      RL(2)=PIN(2)-GG*VN(2)+RN(2)
      RL(3)=PIN(3)-GG*VN(3)+RN(3)
C
      DO  13  KR=1,NREP
      IF(NRAPI.EQ.5)       GO  TO  11
c  redefinition of Fermi momenta
      IF(NRAPI.EQ.1)  THEN
       CALL  PARTN(1,N1,PN1,IPN1)
       CALL  PARTN(2,N2,PN2,IPN2)
       C1=DSQRT((PN1(8)*(PN1(8)+2.*PN1(9)))/
     /      (PN1(4)**2+PN1(5)**2+PN1(6)**2))
       C2=DSQRT((PN2(8)*(PN2(8)+2.*PN2(9)))/
     /      (PN2(4)**2+PN2(5)**2+PN2(6)**2))
       do  k=1,3
       PS1(k)=PN1(3+k)*C1
       PS2(k)=PN2(3+k)*C2
       enddo
       CALL  KINEMA(PS1,VPR,PL1,CT1,ST1,CF1,SF1,TL1,PN1(9))
       CALL  KINEMA(PS2,VTA,PL2,CT2,ST2,CF2,SF2,TL2,PN2(9))
       do  k=1,3
       PIN(3+k)=PL1(k)
       PN(3+k) =PL2(k)
       enddo
       PIN(8)=DSQRT(PL1(1)**2+PL1(2)**2+PL1(3)**2+0.940**2)-0.940 
       PN(8) =DSQRT(PL2(1)**2+PL2(2)**2+PL2(3)**2+0.940**2)-0.940 
       PIN(9)=0.940
       PN(9) =0.940
      ELSE
       CALL  PARTN(NU,N2,PN2,IPN2)
       do  k=1,3
       PN(3+k)=PN2(3+k)
       enddo
       PN(8)=PN2(8)
       PN(9)=0.94
      ENDIF
      PX1=PIN(4)
      PY1=PIN(5)
      PZ1=PIN(6)
      AM1=PIN(9)
      PX2=PN(4)
      PY2=PN(5)
      PZ2=PN(6)
      AM2=PN(9)
      CALL CROSEC(1,ID1,ID2,PX1,PY1,PZ1,AM1,PX2,PY2,PZ2,AM2,SITO,0)
      SIG=SITO
   11 CONTINUE
      CALL  TINVU(PIN,PN,U,V,T12)
      IF(INSTR.NE.2.OR.NRAPI.EQ.5)  GO  TO  12
      IF(NRAPI.EQ.2.OR.NRAPI.EQ.3)  THEN
      CALL  NEWFER(PIN,PN,NU)
      CALL  TINVU(PIN,PN,U,V,T12)
      ENDIF
      CALL  PREPD(U,V,PIN)
   12 CONTINUE
      STOT=SIG+SAB
      CALL CONEC(U,T12,V,VN,PIN(9),PN(9),IIN,IPN,RL,NREP,JPSTR,STOT)
   13 CONTINUE
      RETURN
      END
C
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
      SUBROUTINE NEWFER(PIN,PN,NU)
      IMPLICIT REAL *8(A-H,O-Z)
      REAL *8 PIN,PN
      COMMON/DPFER/DPFER(41)
      DIMENSION  PIN(9),PN(9)
      PF=0.250
      X=PIN(1)
      Y=PIN(2)
      Z=PIN(3)
      TF=TFERMI(X,Y,Z,NU)
      PF=DSQRT(TF*(TF+1.88) )
      D=5.06*6.283185*0.545
      PMAX=2.*PF
   10 DRND=RNDMD(-1)
      P=PMAX*DRND**(1./3.)
      F=(1.+DEXP(-D*PF))/(1.+DEXP(D*(P-PF)))
      IF(RNDMD(-1).GT.F)  GO  TO  10
      IFER=P/0.02+1
      IF(IFER.GT.38) IFER=38
      DPFER(IFER)=DPFER(IFER)+1.
      DPFER(39)=DPFER(39)+P
      DPFER(40)=DPFER(40)+P*P
      DPFER(41)=DPFER(41)+1.
      FI=6.283185*RNDMD(-1)
      CT=1.-2.*RNDMD(-1)
      ST=DSQRT(1.-CT*CT)
      PN(1)=PIN(1)
      PN(2)=PIN(2)
      PN(3)=PIN(3)
      PN(4)=P*ST*DCOS(FI)
      PN(5)=P*ST*DSIN(FI)
      PN(6)=P*CT
      PN(7)=0.
      PN(8)=DSQRT(P**2+0.940**2)-0.940
      PN(9)=0.940
      RETURN
      END
C
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
      SUBROUTINE  PRSTR1(NC,INE,SG)
      IMPLICIT REAL *8(A-H,O-Z)
      REAL *8 SG
      COMMON/INDZAP/INDZAP
      COMMON/INDSTR/JS,NREP,JCA
      COMMON/INPT/INPT
      COMMON/NPRI/NPRI
      COMMON/NPIMI/NPIMI
      COMMON/SPSTR/SPSTR(8,18,41),SPSTR1(8,6,41),DP,DT,DPT,DY,
     *PM,PTM,YMI,YMA
      COMMON/NRESC/RESC(2),STRIN
      COMMON/YSTR/YSTR(8),KSTR(8)
      COMMON/DPFER/DPFER(41)
      COMMON/WSUM/ WSUM
      COMMON /IWSTR/ IWSTR
      DIMENSION  DOM(18),TETF(18)
      CHARACTER*4 TEX1(11),TEX2(11),TEX3(10),TEX4(2),TEX5(9),
     * TEX6(8),TEX7(8),TEX8(8),TEX9(8)
      DIMENSION   C1(18),C2(6)
      DATA TEX1/'INVA','RIAN','T (E','DS/D','P**3',
     *')   ','CROS','S SE','CTIO','NS  ','OF  '/,
     *TEX2/' K+ ',' K0 ','LAMB','SIG+','SIG-','SIG0',' K- ','AP  ',
     *'Pi- ','Pi0 ','Pi+ '/,
     *TEX3/'DIFF','EREN','TIAL','  DI','STRI','BUTI','ONS(','1/ND',
     *'N/DX',') OF'/,
     *TEX4/'FROM','    '/,
     *TEX5/'M+B','PI+N','B+B','N+N ','NNPT','M+D','D+B ','B+D ','ALL '/,
     *TEX6/'DS/D','T/DO','    ','CROS','S SE','CTIO','NS  ','OF  '/,
     *TEX7/'DS/D','P DO','    ','CROS','S SE','CTIO','NS  ','OF  '/,
     *TEX8/'E*DS','/D3P','    ','CROS','S SE','CTIO','NS  ','OF  '/,
     *TEX9/'EdS/','dPTd','Y   ','CROS','S SE','CTIO','NS  ','OF  '/
      DPI=3.141592/18.
      RESC(1)=0.
      RESC(2)=0.
      WSUM=0.
      STRIN=0.
      NPIMI=0
c      DP=0.100
c      IF(JS.LE.2)  DP=DP/2.
c      IF(JS.GE.3.AND.JS.LE.6)  DP=DP*2.
c      IF(JS.EQ.7)  DP=DP/2.
c      IF(JS.EQ.8.and.INPT.EQ.2)  DP=DP/2.
c      IF(JS.EQ.8.and.INPT.NE.2)  DP=DP*2.
c      IF(JS.GE.9)  DP=0.010
c      DT=5.
c      DPT=0.020
c      YMI=-1.
c      DY=0.1
c
      READ(15,*) DP,DT,DPT,YMI,DY 
      WRITE(16,'(1x,''gqstn input:DP,DT,DPT,YMI,DY='',5F7.3)')
     *           DP,DT,DPT,YMI,DY 
c
      DO   M=1,18
      TETF(M)=10.*(M-1)+5.
      DOM(M)=(DCOS((M-1)*DPI)-DCOS(M*DPI))*6.283185
      ENDDO
      PM=36*DP
      PTM=36*DPT
      YMA=YMI+DY*36.
      if(INPT.eq.6) then
      PTM=18*DPT
      DO   M=1,18
        TETF(M)=DPT*(M-1)+.5*DPT
        DOM(M)=DPT*3.146593
      ENDDO
      endif
C
      DO  10  L=1,8
      YSTR(L)=0.
      KSTR(L)=0
      DO  10  N=1,41
      DPFER(N)=0.
      DO  10  M=1,18
      SPSTR(L,M,N)=0.
      IF(M.LE.6)   SPSTR1(L,M,N)=0.
   10 CONTINUE
C
      IF(JCA.GE.0)  GO  TO  11
      RETURN
   11 READ(11)  DCA,DEL,NPIMI,WSUM,RESC,STRIN,YSTR,KSTR,SPSTR,SPSTR1
C     BACKSPACE  11
      REWIND  11
      NC=DCA
      INE=DEL
      WRITE(16,12) DCA
   12 FORMAT(  1X,F7.0,' EVENTS ARE READ FROM FILE 11 '  )
      RETURN
C
C  * * * * * * * * * * * * * * * * * * * * * * * *
C
      ENTRY  PRSTR2(NC,INE,SG)
      IF(INDZAP.EQ.0)   RETURN
      SIN=SG*FLOAT(NC)/FLOAT(NC+INE)
      IF(IWSTR.EQ.1) CY=1./FLOAT(NC)
      IF(IWSTR.EQ.2) CY=1./FLOAT(NC)/SIN
      YMB=YSTR(1)*CY
      YPN=YSTR(2)*CY
      YBB=YSTR(3)*CY
      YNN=YSTR(4)*CY
      Y00=YSTR(5)*CY
      YMD=YSTR(6)*CY
      YDB=YSTR(7)*CY
      YBD=YSTR(8)*CY
      YTOT=YMB+YBB+YMD+YDB+YBD
      SMB=YMB*SIN
      SPN=YPN*SIN
      SBB=YBB*SIN
      SNN=YNN*SIN
      S00=Y00*SIN
      SMD=YMD*SIN
      SDB=YDB*SIN
      SBD=YBD*SIN
      STOT=YTOT*SIN
      DSMB=0.
      IF(KSTR(1).NE.0) DSMB=1./DSQRT(DFLOAT(KSTR(1)))*SMB
      DSPN=0.
      IF(KSTR(2).NE.0) DSPN=1./DSQRT(DFLOAT(KSTR(2)))*SPN
      DSBB=0.
      IF(KSTR(3).NE.0) DSBB=1./DSQRT(DFLOAT(KSTR(3)))*SBB
      DSNN=0.
      IF(KSTR(4).NE.0) DSNN=1./DSQRT(DFLOAT(KSTR(4)))*SNN
      DS00=0.
      IF(KSTR(5).NE.0) DS00=1./DSQRT(DFLOAT(KSTR(5)))*S00
      DSMD=0.
      IF(KSTR(6).NE.0) DSMD=1./DSQRT(DFLOAT(KSTR(6)))*SMD
      DSDB=0.
      IF(KSTR(7).NE.0) DSDB=1./DSQRT(DFLOAT(KSTR(7)))*SDB
      DSBD=0.
      IF(KSTR(8).NE.0) DSBD=1./DSQRT(DFLOAT(KSTR(8)))*SBD
      DSTOT=DSMB+DSBB+DSMD+DSDB+DSBD
      RSC1=0.
      RSC2=0.
      STIN=STRIN/FLOAT(NC*NREP)
      WRITE(16,13) NC,INE,SIN
   13 FORMAT(/1X,'NCAS=',I5,2X,'INTEL=',I5,2X,'INELASTIC CROSS ',
     *'SECTION=',E13.6,'MB'//)
      WRITE(16,112) STIN
  112 FORMAT(1X,'MEAN NUMBER OF STRANGE SOURCES = ',E10.3/)
      IF(STIN.GT.1.E-20)   RSC1=RESC(1)/FLOAT(NC*NREP)/STIN
      IF(STIN.GT.1.E-20)   RSC2=RESC(2)/FLOAT(NC*NREP)/STIN
      WRITE(16,113) RSC1,RSC2
  113 FORMAT(/1X,'MEAN NUMBER OF RESCATTERING = ',E10.3,2X,E10.3//)
      WRITE(16,14) TEX2(JS),
     *YMB,SMB,DSMB,YPN,SPN,DSPN,YBB,SBB,DSBB,YNN,SNN,DSNN,
     *Y00,S00,DS00,YMD,SMD,DSMD,YDB,SDB,DSDB,YBD,SBD,DSBD,
     *YTOT,STOT,DSTOT
c     WRITE( *,14) TEX2(JS),
c    *YMB,SMB,DSMB,YPN,SPN,DSPN,YBB,SBB,DSBB,YNN,SNN,DSNN,
c    *Y00,S00,DS00,YMD,SMD,DSMD,YDB,SDB,DSDB,YBD,SBD,DSBD,
c    *YTOT,STOT,DSTOT
   14 FORMAT(/1X,'MEAN MULTIPLICITY (PROD.CROSS SECTION) OF ',A4/
     *2X,'FROM     M+B=',2(E13.6,2X),'+/-',E13.6/
     *2X,'FROM    PI+N=',2(E13.6,2X),'+/-',E13.6/
     *2X,'FROM     B+B=',2(E13.6,2X),'+/-',E13.6/
     *2X,'FROM     N+N=',2(E13.6,2X),'+/-',E13.6/
     *2X,'FROM NP0+NT0=',2(E13.6,2X),'+/-',E13.6/
     *2X,'FROM   M+"D"=',2(E13.6,2X),'+/-',E13.6/
     *2X,'FROM   "D"+B=',2(E13.6,2X),'+/-',E13.6/
     *2X,'FROM   B+"D"=',2(E13.6,2X),'+/-',E13.6/
     *2X,'       TOTAL=',2(E13.6,2X),'+/-',E13.6//)
      PIMI=FLOAT(NPIMI)/FLOAT(NC)
      WRITE(16,114) PIMI
  114 FORMAT(/5X,'MEAN MULTIPLICITY OF PI-MINUS=',E13.6//)
      IF(INDZAP.EQ.2)  GO  TO  1994
      IF(((NC/NPRI)*NPRI).NE.NC)  GO  TO  124
 1994 CONTINUE
      DO  19  K=1,9
      IF(K.LE.8.AND.KSTR(K).EQ.0)  GO  TO  19
      IF(INPT.LE.2)  WRITE(16,18) TEX1,TEX2(JS),TEX4,TEX5(K)
      IF(INPT.EQ.3)  WRITE(16,28) TEX6,TEX2(JS),TEX4,TEX5(K)
      IF(INPT.EQ.4)  WRITE(16,28) TEX7,TEX2(JS),TEX4,TEX5(K)
      IF(INPT.EQ.5)  WRITE(16,28) TEX8,TEX2(JS),TEX4,TEX5(K)
      IF(INPT.EQ.6)  WRITE(16,28) TEX9,TEX2(JS),TEX4,TEX5(K)
      DO  119  MP=1,2
      M1=(MP-1)*9+1
      M2=(MP-1)*9+9
      IF(INPT.EQ.1)   WRITE(16,118) (TETF(MM),MM=M1,M2)
      IF(INPT.EQ.2)   WRITE(16,128) (TETF(MM),MM=M1,M2)
      IF(INPT.EQ.3)   WRITE(16,138) (TETF(MM),MM=M1,M2)
      IF(INPT.EQ.4)   WRITE(16,148) (TETF(MM),MM=M1,M2)
      IF(INPT.EQ.5)   WRITE(16,158) (TETF(MM),MM=M1,M2)
      IF(INPT.EQ.6)   WRITE(16,159) (TETF(MM),MM=M1,M2)
      DPXX=DP
      if(INPT.eq.6)  DPXX=DY
      DO 119  L=1,41
      DO  15  M=M1,M2
      IF(K.LE.8)  C1(M)=SPSTR(K,M,L)*SIN*CY
      IF(K.EQ.9)  C1(M)=(SPSTR(1,M,L)+SPSTR(3,M,L)+
     +SPSTR(6,M,L)+SPSTR(7,M,L)+SPSTR(8,M,L))*SIN*CY
      C1(M)=C1(M)/DPXX/DOM(M)
   15 CONTINUE
      P=DP*(L-1)+DP/2.
      if(INPT.eq.6) P=YMI+DY*(L-1)+0.5*DY
      IF(L.LE.36)  WRITE(16,16) P,(C1(MM),MM=M1,M2)
      IF(L.GT.36)  WRITE(16,17)   (C1(MM),MM=M1,M2)
   16 FORMAT(3X,F5.3,2X,9(E10.3,2X))
   17 FORMAT(10X,9(E10.3,2X))
   18 FORMAT(//5X,15A4/)
   28 FORMAT(//5X,12A4/)
  118 FORMAT(/3X,'P/TETA',4X,9(F4.0,8X)/)
  128 FORMAT(/2X,'T*/TETA',4X,9(F4.0,8X)/)
  138 FORMAT(/2X,'TL/TETA',4X,9(F4.0,8X)/)
  148 FORMAT(/2X,'PL/TETA',4X,9(F4.0,8X)/)
  158 FORMAT(/2X,'TL/TETA',4X,9(F4.0,8X)/)
  159 FORMAT(/2X,' Y/PT  ',4X,9(F5.3,7X)/)
  119 CONTINUE
   19 CONTINUE
C      
      DO 2424 K=1,9
        IF(K.LE.8.AND.KSTR(K).EQ.0)GO TO 2424
        WRITE(16,20) TEX3,TEX2(JS),TEX4,TEX5(K)
        WRITE(16,120)
   20   FORMAT(//5X,14A4/)
  120   FORMAT(/2X,'EL',4X,'dS/dEL',4x,'dS/dEnn',3X,'cosT',2X,
     *  ' dS/DOl',3x,'dsdO*',4X,'PT',2X,'ds/p_tdp_t',3X,'Y',4X,'ds/dY'/
     *  1X,'MeV',4X,'mb/MeV',4x,'mb/MeV ',9X,
     *  ' mb/sr ',3x,'mb/sr',3X,'MeV/c',1X,'nb/(Mev/c)**2',4X,'mb')
        DO 24 L=1,41
          DO  21  M=1,6
            C2(M)=0.
            IF(K.LE.8)C2(M)=SPSTR1(K,M,L)*SIN*CY
            IF(K.EQ.9)C2(M)=(SPSTR1(1,M,L)+SPSTR1(3,M,L)+
     *                SPSTR1(6,M,L)+SPSTR1(7,M,L)+
     *                SPSTR1(8,M,L))*SIN*CY
   21     CONTINUE
        EL=(DP*(L-1)+DP/2.)*1000.
        IF(JS.GE.9.AND.JS.LE.11)EL=EL+140.
        PT=(DPT*(L-1)+DPT/2.)*1000.
        Y=YMI+(L-1)*DY+DY/2.
        CTS=(L-1)*0.1+0.05-1.
        ELNN=EL
        IF(L.LE.20)WRITE(16,22)EL,C2(1),C2(6),CTS,C2(2),C2(5),
     *                         PT,C2(3),Y,C2(4)
        IF(L.GT.20.AND.L.LE.36)WRITE(16,122)EL,C2(1),C2(6),
     *                                      PT,C2(3),Y,C2(4)
        IF(L.GT.36)WRITE(16,23) C2
   22   FORMAT(0PF6.1,2(1PE10.3),0PF6.2,2(1PE10.3),
     *  0PF5.0,1PE10.3,0PF6.2,1PE10.3)
  122   FORMAT(0PF6.1,2(1PE10.3),26X,
     *  0PF5.0,1PE10.3,0PF6.2,1PE10.3)
   23   FORMAT(5X,    2(1PE10.3),6X    ,2(1PE10.3),
     *  5X,    1PE10.3,6X    ,1PE10.3)
   24   CONTINUE
 2424 CONTINUE
C
  124 CONTINUE
C     IF(INDZAP.EQ.0)   RETURN
      IF(DPFER(41).GT.0.)  WSM=WSUM/DPFER(41)
      WRITE(16,*) 'INDZAP,WSM=', INDZAP,WSM
      DCA=FLOAT(NC)
      DEL=FLOAT(INE)
      WRITE(11) DCA,DEL,NPIMI,WSUM,RESC,STRIN,YSTR,KSTR,SPSTR,SPSTR1
C     BACKSPACE 11
      REWIND  11
   25 FORMAT(//3X,I5,' EVENTS ARE WRITTEN ON FILE 11'//)
      C=DPFER(41)
      DLP=0.020
      DO  125  K=1,40
      PF=FLOAT(K-1)*DLP+DLP/2.
      IF(C.NE.0.)  DPFRM=DPFER(K)/C
      IF(K.EQ.39)  AVRP=DPFRM
      IF(K.EQ.40)  DPFRM=DSQRT(DPFRM-AVRP**2)
      IF(K.LE.38)  DPFRM=DPFRM /DLP
      IF(K.LE.38)  WRITE(16,421) PF,DPFRM
      IF(K.EQ.39)  WRITE(16,422) DPFRM
      IF(K.EQ.40)  WRITE(16,423) DPFRM
  125 CONTINUE
  421 FORMAT(1X,F7.3,1X,E10.3)
  422 FORMAT(1X,'<P>    ',1X,E10.3)
  423 FORMAT(1X,'<DISP> ',1X,E10.3)
      RETURN
      END
C
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
      SUBROUTINE  POSTR(P,CT,AM,V,B,JP,ST,STOT,RL,N,INDEX)
      IMPLICIT REAL *8(A-H,O-Z)
      REAL *8 P,CT,AM,V,B,ST,STOT,RL
      DIMENSION  V(3),B(3),PS(3),PL(3),RL(3)
      COMMON/SPSTR/SPSTR(8,18,41),SPSTR1(8,6,41),DP,DT,DPT,DY,
     *PM,PTM,YMI,YMA
      COMMON/CMNN/BNN,GNN
      COMMON/INPT/INPT
      COMMON/NRAPI/NRAPI
      COMMON/YSTR/YSTR(8),KSTR(8)
      COMMON/BGLAB/BLAB,GLAB,KOBR
      COMMON/KSYST/KSYST
      COMMON/DEUT/VDN(3),UDN,T1DN,WDN,VND(3),UND,T1ND,WND
      COMMON/IPRIMC/IPRIMC
      COMMON /IR1IR2/ IR1,IR2
      COMMON /IWSTR/ IWSTR
C     WRITE(16,*) 'AM,ST,INDEX=', AM,ST,INDEX
      IND=INDEX
      BM2=B(1)**2+B(2)**2+B(3)**2
      SIT=DSQRT(1.-CT*CT)
      FI=6.283185*RNDMD(-1)
      PS(1)=P*SIT*DCOS(FI)
      PS(2)=P*SIT*DSIN(FI)
      PS(3)=P*CT
      CALL  KINEMA(PS,V,PL,CTL,STL,CFL,SFL,TL,AM)
      IF(BM2.LT.1.0E-7)   GO  TO  11
      CALL  KINEMA(PL,B,PS,CTL,STL,CFL,SFL,TL,AM)
      DO  10  K=1,3
      PL(K)=PS(K)
   10 CONTINUE
   11 CONTINUE
      IF(IWSTR.EQ.1) W=ST/STOT/FLOAT(N)
      IF(IWSTR.EQ.2) W=ST/FLOAT(N)
      IF(IND.GT.5)  GO TO 12
      IF(NRAPI.EQ.1)  W=W*(1.-WND-WDN)
      IF(NRAPI.EQ.2)  W=W*(1.-WDN)
      IF(NRAPI.EQ.3)  W=W*(1.-WND)
   12 CONTINUE
      CALL  RESCA1(PL,AM,JP,RL,CTL,STL,CFL,SFL,TL,W,WABS,TAUABS)
C      write( *,*) 'out resca1'
      GO  TO  (1,2,3),KSYST
    1 VSYS=0.
      GO  TO  4
    2 T0=(GLAB-1.)*0.940
      P0=SQRT(T0*(T0+1.88))
      VSYS=P0/(T0+1.88)
      GO  TO  4
    3 VSYS=BNN
      GO  TO  4
    4 GSYS=1./DSQRT(1.-VSYS**2)
      PL(3)=GSYS*(PL(3)+(TL+AM)*VSYS)
      PML=DSQRT(PL(1)**2+PL(2)**2+PL(3)**2)
      CTL=PL(3)/PML
      STL=DSQRT(1.-CTL**2)
      TL=DSQRT(PML**2+AM**2)-AM
      IF(KOBR.EQ.0)  GO  TO  111
      PZL=GLAB*(-PL(3)+(TL+AM)*BLAB)
      PML=DSQRT(PL(1)**2+PL(2)**2+PZL**2)
      CTL=PZL/PML
      STL=DSQRT(1.-CTL**2)
      TL=DSQRT(PML**2+AM**2)-AM
  111 CONTINUE
      PML=DSQRT(TL*(TL+2.*AM))
      EL=TL+AM
      ES=GNN*(EL-PML*CTL*BNN)
      TS=ES-AM
      PMS=DSQRT(TS*(TS+2.*AM))
      PZS=GNN*(PML*CTL-EL*BNN)
      CTS=PZS/PMS
      PT=PML*STL
      W2=W/(DP*1000.)
      W3=W/(2.*3.141592*0.1)
      WY=W/DY
C
      WT=W
      WT=0.
      IF(PT.GT.0.)  WT=W/PT
C
      IF(INPT.EQ.1)  W1=W*EL/PML**2
      IF(INPT.EQ.2)  W1=W*GNN*DABS(PML-EL*BNN*CTL)/PML**2
      IF(INPT.EQ.3)  W1=W
      IF(INPT.EQ.4)  W1=W
      IF(INPT.EQ.5)  W1=W/PML
c      IF(INPT.EQ.6)  W1=WT       !for invar. cross section
      IF(INPT.EQ.6)  W1=W        ! for abundancy
      Y=0.5*DLOG((TL+AM+PML*CTL)/(TL+AM-PML*CTL))
      IF(DABS(CTL).GT.1.)   CTL=DSIGN(1.D0,CTL)
      TETA=DACOS(CTL)*180./3.141592
      if(INPT.eq.6) then
      JT=PT/DPT+1
      else
      JT=TETA/10.+1
      endif
      IF(JT.GT.18)  JT=18
      JJRA=1
      IF(IND.EQ.1.AND.IR1.EQ.0.AND.IR2.EQ.0)   JJRA=2
      IF(IND.EQ.3.AND.IR1.EQ.0.AND.IR2.EQ.0)   JJRA=2
      IF(IND.EQ.3.AND.NRAPI.EQ.1)   JJRA=3
      IF(IND.EQ.3.AND.IPRIMC.EQ.1)  JJRA=3
C     WRITE(16,*) 'JJRA,TL=', JJRA,TL
      DO  15  JJ=1,JJRA
      IND=INDEX+JJ-1
      YSTR(IND)=YSTR(IND)+W
      KSTR(IND)=KSTR(IND)+1
C     IF(DABS(W).GT.1.) WRITE(16,1991) IND,ST,STOT,W,WND,WDN
1991  FORMAT(1X,'IND,ST,STOT,W,WND,WDN=',I2,5(1X,E10.3))
      IF(INPT.EQ.1) CALL DHIST4(PML,0.D0,PM,DP,SPSTR,8,18,41,W1,IND,JT)
      IF(INPT.EQ.2) CALL DHIST4( TS,0.D0,PM,DP,SPSTR,8,18,41,W1,IND,JT)
      IF(INPT.EQ.3) CALL DHIST4( TL,0.D0,PM,DP,SPSTR,8,18,41,W1,IND,JT)
      IF(INPT.EQ.4) CALL DHIST4(PML,0.D0,PM,DP,SPSTR,8,18,41,W1,IND,JT)
      IF(INPT.EQ.5) CALL DHIST4( TL,0.D0,PM,DP,SPSTR,8,18,41,W1,IND,JT)
      IF(INPT.EQ.6) CALL DHIST4(  Y,YMI,YMA,DY,SPSTR,8,18,41,W1,IND,JT)
      CTL1=CTL+1.
      CTS1=CTS+1.
      WT1=WT/DPT
      CALL  DHIST4(TL,0.D0,PM,DP,SPSTR1,8,6,41,W2,IND,1)
      CALL  DHIST4(CTL1,0.D0,2.D0,0.1D0,SPSTR1,8,6,41,W3,IND,2)
      CALL  DHIST4(PT,0.D0,PTM,DPT,SPSTR1,8,6,41,WT1,IND,3)
      CALL  DHIST4( Y,YMI,YMA,DY,SPSTR1,8,6,41,WY,IND,4)
      CALL  DHIST4(CTS1,0.D0,2.D0,0.1D0,SPSTR1,8,6,41,W3,IND,5)
      CALL  DHIST4(TS,0.D0,PM,DP,SPSTR1,8,6,41,W2,IND,6)
   15 CONTINUE
C      write( *,*) 'postr out'
      RETURN
      END
C
C * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
      SUBROUTINE  CONEC(U,T,V,B,AMP,AMT,IIN,IPN,RL ,N,JP,STOT)
      IMPLICIT REAL *8(A-H,O-Z)
      REAL *8 U,T,V,B,AMP,AMT,RL,STOT
      COMMON /KSTST/ KSTST
      COMMON /NRAPI/ NRAPI
      COMMON /INSTR/ INSTR
      COMMON /AMPROJ/ AMPROJ
      COMMON /IDN12/ ID1,ID2
      COMMON /IR1IR2/ IR1,IR2
      COMMON /KSIGS/ KSIGS
      DIMENSION  V(3),B(3),RL(3) ,AM(11),IIN(5),IPN(5),SIJ(4),PIS(3)
      DATA  AM/2*0.494,1.115,3*1.192,0.494,0.939,3*0.139/
      IE1=IIN(1)
      IE2=IPN(1)
      IB1=IIN(4)
      IB2=IPN(4)
      IR1=0
      IF(IIN(5).NE.0)  IR1=1
      IR2=0
      IF(IPN(5).NE.0)  IR2=1
      IBA=IB1+IB2
      IF(NRAPI.EQ.5)  GO  TO  9
      IF(INSTR.LE.2)  GO  TO  8
C  FOR FLUCTONS WITH BOLTZMANN SPECTRA
      IF(INSTR.NE.6)  GO  TO  9
      CALL  POSFL(IB1,IE1,IE2,JP,B,STOT,RL,N)
      GO  TO  9
    8 IF(IBA.EQ.2.OR.IBA.EQ.1)
     *CALL  POSND(IB1,IE1,IE2,JP,B,STOT,RL,N)
    9 CONTINUE
      IF(JP.GE.9)      GO  TO  900
      IF(JP.EQ.8)      GO  TO  800
      AM1=AM(JP)
      AM2=1.115
      IF(JP.GT.2)      AM2=0.494
      IF(JP.EQ.7)      GO  TO  700
      DELTM=AM1+AM2+FLOAT(IB1)*0.940
      IF(U.LE.DELTM)   GO  TO  906
      IF(IBA.NE.2)     GO  TO  110
C
      ST=0.
      UN=U-AMP-AMT+1.880
      IF(KSIGS.GE.2) THEN
      ST=SZWER(JP,UN)
      IF((IR1+IR2).EQ.1)  ST=ST*0.75
      IF((IR1+IR2).EQ.2)  ST=ST*0.50
      SK=ST
      SDR=ST
      SL=ST
      SLDR=ST
      GO  TO  12
      ENDIF
C
      IF(IR1.EQ.0.AND.IR2.EQ.0) GO  TO 91
C
C BARYON+BARYON PRODUCTION CROSS SECTION
C !!!!
      TNN=(U**2-4.*0.940**2)/1.88
      ST=(SST(1,1,1,1,JP,TNN,1)+SST(1,1,1,0,JP,TNN,1)+
     +    SST(1,1,0,1,JP,TNN,1)+SST(1,1,0,0,JP,TNN,1))/4.
      IF((IR1+IR2).EQ.1)  ST=ST*0.75
      IF((IR1+IR2).EQ.2)  ST=ST*0.50
      SK=ST
      SL=ST
      SDR=ST
      SLDR=ST
      IF(SK.LT.1.E-20)     GO  TO  906
      GO  TO  12
   91 IF(T.LT.4.)      GO  TO  92
C
C NUCLEON+NUCLEON PRODUCTION CROSS SECTION FOR T0 > 4. GeV
      ST4N=SST(IB1,IB2,IE1,IE2,JP,4.D0,1)
      ST4D=SSDT(IB1,IB2,IE1,IE2,JP,4.D0,1)
      ST4=ST4N
      IF((RNDMD(-1)*(ST4N+ST4D)).GT.ST4N) ST4=ST4D
      ST=STHE(U,JP,ST4)
      SDR=ST
      SL=ST
      SLDR=ST
      GO  TO  10
  110 CONTINUE
C MESON+BARYON PRODUCTION CROSS SECTION
C
      ST=0.
      IF(KSIGS.EQ.3.OR.KSIGS.EQ.4) THEN
      ST=SCT(IE1,IE2,JP,T,1)
      SL=SCL(IE1,IE2,JP,T)
C
      GO  TO  10
      ENDIF
      IF(KSIGS.EQ.5) THEN
      ST=0.
      SS0=0.494+1.115
      IF(U.GT.SS0.AND.U.LT.1.7)  ST=10.*(U-SS0)
      IF(U.GE.1.7)               ST=0.09/(U-1.6)
      SL=ST
      GO  TO  10
      ENDIF
C
C
      IF((ID1.EQ.120.OR.ID1.EQ.-120.OR.ID1.EQ.110).AND.IR2.EQ.0)
     * GO  TO  92
      UPN=U-AMP-AMT+0.940+0.140
      TPN=(U**2-(0.140+0.940)**2)/1.88
      ST=(SST(0,1, 1,1,JP,TPN,1)+SST(0,1, 1,0,JP,TPN,1)+
     +    SST(0,1,-1,1,JP,TPN,1)+SST(0,1,-1,0,JP,TPN,1)+
     +    SST(0,1, 0,1,JP,TPN,1)+SST(0,1, 0,0,JP,TPN,1))/6.
      SL=ST
      GO  TO  10
C
   92 CONTINUE
C NUCLEON+NUCLEON PRODUCTION CROSS SECTION   FOR T < 4. GeV
C   PION+NUCLEON  PRODUCTION CROSS SECTION
      SDR=SSDT(IB1,IB2,IE1,IE2,JP,T,1)
      SLDR=SSDL(IB1,IB2,IE1,IE2,JP,T)
      ST=SST(IB1,IB2,IE1,IE2,JP,T,1)
      SL=SSL(IB1,IB2,IE1,IE2,JP,T)
   10 SK=ST
      IF(SK.LT.1.E-20)     GO  TO  906
      IF((IB1+IB2).EQ.2)   GO  TO  12
C   * * * * * * * * * * * * * * * * * * * * * * * * * * *
C    M+B=>Y+K,   Y=L,S
      IND=1
C
      AM2=1.115
      IF(JP.GT.2)                AM2=0.494
      IF(JP.LE.2.AND.(U-DELTM-(1.192-1.115)).GT.0.) THEN
C       IF((RNDMD(-1)*ST).GT.SL) AM2=1.192
      IF(ST.NE.SL)             AM2=1.192
      ENDIF
      E1=(U**2+AM1**2-AM2**2)/(2.*U)
      P1=DSQRT(E1**2-AM1**2)
      IF(JP.LE.2)  CT1=COSEX(0,T,0.140D0)
      IF(JP.GT.2)  CT1=-COSEX(0,T,0.140D0)
      CT1=DABS(CT1)
      IF(RNDMD(-1).LE.0.5)  CT1=-CT1
      CALL  POSTR(P1,CT1,AM1,V,B,JP,SK,STOT,RL,N,IND)
   11 CONTINUE
C      WRITE( *,1101) ID1,ID2,IND,T,U,SK,STOT
 1101 FORMAT(1X,'M+B: ID1,ID2,IND,T,U,SK,STOT=',2I5,I2,4(E11.4))
      GO  TO  906
C   * * * * * * * * * * * * * * *  * * * * * * * * *
   12 CONTINUE
C    B+B=>B1+Y+K, Y=L,S;  B,B1=N,DELTA
      IND=3
C
  333 NREP=0
      IF(IR1.EQ.0.AND.IR1.EQ.0) GO  TO 121
      AM3=0.939
      GO  TO  206
  121 CONTINUE
C    N+N=>B+Y+K, Y=L,S;  B=N,DELTA
      AM3=0.939
      IF(JP.GT.2)          GO  TO 204
C  ******  CHANNELS FOR K+ AND K0  PRODUCTION
      AM2=1.115
      AM3=0.939
      IF(SDR.GT.1.E-6)     GO  TO 200
      IF((U-DELTM-(1.192-1.115)).GT.0.) THEN
C       IF((RNDMD(-1)*ST).GT.SL)  AM2=1.192
      IF(ST.NE.SL)  AM2=1.192
      ENDIF
      GO  TO  205
  200 CONTINUE
C     SK=ST+SDR
      SK=ST
      RSR=RNDMD(-1)*(ST+SDR)
      IF((U-DELTM-(1.236-0.939)).GT.0..AND.RSR.GT.SDR) THEN
        AM3=1.236
        SK=SDR
C         IF((RNDMD(-1)*SDR).GT.SLDR)   AM2=1.192
        IF(SDR.NE.SLDR)   AM2=1.192
        GO TO 205
      ELSE
C         IF((RNDMD(-1)*ST).GT.SL)      AM2=1.192
        IF(ST.NE.SL)      AM2=1.192
      ENDIF
      GO  TO   205
C  ******  CHANNELS FOR LAMBDA AND SIGMA PRODUCTION
  204 AM2=0.494
      AM3=0.939
C      SK=ST+SDR
      SK=ST
      RSR=RNDMD(-1)*(ST+SDR)
      IF((U-DELTM-(1.236-0.939)).GT.0..AND.RSR.GT.ST)  THEN
      AM3=1.236
      SK=SDR
      ENDIF
  205 NREP=NREP+1
      IF(U.LE.(AM1+AM2+AM3).AND.NREP.LT.10)  GOTO 333
      IF(U.LE.(AM1+AM2+AM3)) THEN
      WRITE(16,*) 'FROM CONEC:', ID1,ID2,U,AM1,AM2,AM3
c     WRITE( *,*) 'FROM CONEC:', ID1,ID2,U,AM1,AM2,AM3
      GO  TO  906
      ENDIF
      IF(KSTST.NE.1)  GO  TO  207
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C       SAMPLING ACCORDING TO 3-BODY PHASE SPACE VOLUME
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  206 CONTINUE
      IF(U.LE.(AM1+AM2+AM3)) WRITE(16,*) ID1,ID2,U,AM1,AM2,AM3
c     IF(U.LE.(AM1+AM2+AM3)) WRITE( *,*) ID1,ID2,U,AM1,AM2,AM3
      IF(U.LE.(AM1+AM2+AM3)) GO  TO  906
      EM1=(U**2+AM1**2-(AM2+AM3)**2)/(2.*U)
      IF(EM1.LE.AM1)         GO  TO  906
      EM2=(U**2+AM2**2-(AM1+AM3)**2)/(2.*U)
      IF(EM2.LE.AM2)         GO  TO  906
      E1=AM1+(EM1-AM1)*RNDMD(-1)
      E2=AM2+(EM2-AM2)*RNDMD(-1)
      E3=U-E1-E2
      IF(E3.LE.AM3)          GO  TO  206
      FNORM=27.*E1*E2*E3/U**3
      IF(RNDMD(-1).GT.FNORM)  GO  TO  206
      P1=DSQRT(E1**2-AM1**2)
      P2=DSQRT(E2**2-AM2**2)
      P3=DSQRT(E3**2-AM3**2)
      IF(((P1+P2-P3)*(P1-P2+P3)*(P2+P3-P1)).LE.0.)  GO  TO  206
      CT1=1.-2.*RNDMD(-1)
      GO  TO  14
C * * * * * *  * * * * * * * * * * * * * * * * * *  * * * * * *
C       SAMPLING ACCORDING TO EXPERIMENTAL APPROXIMATION
C       ASSUMING THAT SIMILARITY N+N=>N+L+K AND N+N=>N+N+PI
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  207 CONTINUE
      U1=U-AM1-AM2-AM3+2.*0.939+0.140
      TP=U1*U1/2./0.939-2.*0.939
  300 IF(JP.GT.2)   GO  TO  301
      IF(TP.LT.4.)   P1=PMOM(2,TP)
      IF(TP.GE.4.)   P1=PMOM(4,TP)
      E1=DSQRT(P1**2+0.14**2)-0.14+AM1
      IF(TP.LT.4.)   P2=PMOM(1,TP)
      IF(TP.GE.4.)   P2=PMOM(3,TP)
      E2=DSQRT(P2**2+0.94**2)-0.94+AM2
      GO  TO  303
  301 IF(TP.LT.4.)   P1=PMOM(1,TP)
      IF(TP.GE.4.)   P1=PMOM(3,TP)
      E1=DSQRT(P1**2+0.94**2)-0.94+AM1
      IF(TP.LT.4.)   P2=PMOM(2,TP)
      IF(TP.GE.4.)   P2=PMOM(4,TP)
      E2=DSQRT(P2**2+0.14**2)-0.14+AM2
  303 E3=U-E1-E2
      IF(E3.LE.AM3)   GO  TO  300
      P1=DSQRT(E1**2-AM1**2)
      P2=DSQRT(E2**2-AM2**2)
      P3=DSQRT(E3**2-AM3**2)
      IF(T.GE.4.)  GO  TO  304
      IF(((P1+P2-P3)*(P1-P2+P3)*(P2+P3-P1)).LE.0.)  GO  TO  300
  304 IF(JP.LE.2.AND.TP.LT.4.)   CT1=COSTA(15,TP)
      IF(JP.LE.2.AND.TP.GE.4.)   CT1=COSTA(17,TP)
      IF(JP.GT.2.AND.TP.LT.4.)   CT1=COSTA(14,TP)
      IF(JP.GT.2.AND.TP.GE.4.)   CT1=COSTA(16,TP)
      CT1=DABS(CT1)
      IF(RNDMD(-1).LE.0.5)  CT1=-CT1
   14 CONTINUE
C
      IF(AM3.EQ.0.939)  THEN
      CALL  PAULBL(P3,V,IBL)
      IF(IBL.EQ.1)  GO  TO  114
      ENDIF
C
      CALL  POSTR(P1,CT1,AM1,V,B,JP,SK,STOT,RL,N,IND)
  114 CONTINUE
C      WRITE( *,1401) ID1,ID2,IND,T,U,SK,STOT
 1401 FORMAT(1X,'B+B: ID1,ID2,IND,T,U,SK,STOT=',2I5,I2,4(E11.4))
      GO  TO  906
C
C  ******  CHANNELS FOR  K-  PRODUCTION  G: 20.05.93
  700 AM1=0.494
      EMN=(IB1+IB2)*0.939+2.*AM1
      IF(U.LE.EMN)  GO  TO  906
      EMN1=EMN-AM1
      PMAX=DSQRT(((U+AM1)**2-EMN1**2)*((U-AM1)**2-EMN1**2))/(2.*U)
      IF(IBA.EQ.2)    GO  TO  710
      SK=PMAX/13.
      GO  TO  711
  710 SK=SIKMI(PMAX)
      IF(T.LT.4.)    GO  TO  711
      U4=1.88*(4.+1.88)
      PMA4=DSQRT(((U4+AM1)**2-EMN1**2)*((U4-AM1)**2-EMN1**2))/(2.*U4)
      SK4=SIKMI(PMA4)
      ST=STHE(U,7,SK4)
  711 CONTINUE
C
      IF(KSTST.NE.1)       GO  TO  701
      EM1=EMN1+RNDMD(-1)*(U-EMN)
      E1=(U**2+AM1**2-EM1**2)/(2.*U)
      P1=DSQRT(E1**2-AM1**2)
      CT1=1.-2.*RNDMD(-1)
      IF(IBA.EQ.1) IND=1
      IF(IBA.EQ.2) IND=3
      GO  TO  124
C * * * * * * * * * * * * * * * * * * * * * * * * * * *
  701 X=RNDMD(-1)
      P1=PMAX*X
      FB=16.*(X*(1.-X))**2
      IF(RNDMD(-1).GT.FB)  GO  TO  701
      E1=DSQRT(P1**2+AM1**2)
      CT1=1.-2.*RNDMD(-1)
      IF(IBA.EQ.1) IND=1
      IF(IBA.EQ.2) IND=3
  124 CALL  POSTR(P1,CT1,AM1,V,B,JP,SK,STOT,RL,N,IND)
      GO  TO  906
C
C  ******  CHANNELS FOR  ANTI-P  PRODUCTION
  800 X=1.0
      AM1=0.939*X
      EMN=(IB1+IB2)*0.939+2.*AM1
      EMN1=EMN-AM1
      IF(U.LE.EMN)  GO  TO  906
      PMAX=DSQRT(((U+AM1)**2-EMN1**2)*((U-AM1)**2-EMN1**2))/(2.*U)
C !!!!
C  A. SHOR NP A514(1990) p.717;  DANILIEWICZ PR C42(1990) p.1564
C  SEE: PL B256(1991) p.331
      ST=0.01*(U-EMN)**1.846
      SK=ST
C !!!!
C
      EM1=EMN1+RNDMD(-1)*(U-EMN)
      E1=(U**2+AM1**2-EM1**2)/(2.*U)
      P1=DSQRT(E1**2-AM1**2)
      CT1=1.-2.*RNDMD(-1)
      IF(IBA.EQ.1) IND=1
      IF(IBA.EQ.2) IND=3
c
c    additional PAULBL-like should be inserted here !
c
  801 CALL  POSTR(P1,CT1,AM1,V,B,JP,SK,STOT,RL,N,IND)
  802 CONTINUE
      GO  TO  906
  900 CONTINUE
C  ******  CHANNELS FOR  PI-0+  PRODUCTION
      AM1=AM(JP)
      EMN=(IB1+IB2)*0.940+AM1
      IF(U.LE.EMN)  GO  TO  906
      IF((IE1+IE2).EQ.0.AND.JP.EQ.11) GO  TO  906
      IF((IE1+IE2).EQ.2.AND.JP.EQ. 9) GO  TO  906
      IF(KSIGS.EQ.1.OR.(IB1+IB2).NE.2)   GO  TO  902
C
      CALL SARNDT(U,SIJ)
      IF(IE1.EQ.IE2.AND.JP.EQ.10)  SK=SIJ(1)
      IF(IE1.EQ.IE2.AND.JP.NE.10)  SK=SIJ(1)+SIJ(2)
      IF(IE1.NE.IE2.AND.JP.EQ.10)  SK=(SIJ(2)+SIJ(3))/2.
      IF(IE1.NE.IE2.AND.JP.EQ.11)  SK=(SIJ(1)+SIJ(3))/2.
      IF(IE1.NE.IE2.AND.JP.EQ.9)   SK=(SIJ(1)+SIJ(3))/2.
      IF(SK.LE.0.)  GO  TO  906
      CALL  PIIZO(U,PIS,IFIZO)
      IF(IFIZO.EQ.0)                  GO  TO  906
      P1=DSQRT(PIS(1)**2+PIS(2)**2+PIS(3)**2)
      E1=DSQRT(P1**2+AM1**2)
      CT1=PIS(3)/P1
      IND=3
C     WRITE( *,*)  'PI: T,U,SK=',T,U,SK
      GO  TO  905
C
  902 CONTINUE
      KSI=1
      IF(IBA.EQ.2.AND.(IE1+IE2).EQ.1)          KSI=2
      IF(IBA.EQ.1.AND.IE1.EQ.-1.AND.IE2.EQ.1)  KSI=2
      IF(IBA.EQ.1.AND.IE1.EQ. 0.AND.IE2.EQ.1)  KSI=2
      IF(IBA.EQ.1.AND.IE1.EQ. 0.AND.IE2.EQ.0)  KSI=3
      IF(IBA.EQ.1.AND.IE1.EQ. 1.AND.IE2.EQ.0)  KSI=3
      IF(JP.EQ.10) SK= CROSEG(0,0,IBA,KSI,4,T,AMPROJ,0)
      IF(JP.NE.10) SK= CROSEG(0,0,IBA,KSI,5,T,AMPROJ,0)
      IF(SK.LE.0.)  GO  TO  906
      AM2=0.939
      IF(IBA.EQ.1)      AM3=0.139
      IF(IBA.EQ.2)      AM3=0.939
      IF(IBA.EQ.1)      IND=1
      IF(IBA.EQ.2)      IND=3
C      WRITE( *,*)  'PI: IND,IB1,IB2,T,U,SK=',IND,IB1,IB2,T,U,SK
C
  903 CONTINUE
      IF(IBA.EQ.1.AND.T.LT.4.)   P1=PMOM(7,T)
      IF(IBA.EQ.2.AND.T.LT.4.)   P1=PMOM(2,T)
      IF(IBA.EQ.1.AND.T.GE.4.)   P1=PMOM(5,T)
      IF(IBA.EQ.2.AND.T.GE.4.)   P1=PMOM(4,T)
      E1=DSQRT(P1**2+AM1**2)
      IF(T.LT.4.)   P2=PMOM(1,T)
      IF(T.GE.4.)   P2=PMOM(3,T)
      E2=DSQRT(P2**2+AM2**2)
      E3=U-E1-E2
      IF(E3.LE.AM3)   GO  TO  903
      P3=DSQRT(E3**2-AM3**2)
      IF(T.GE.4.)  GO  TO  904
      IF(((P1+P2-P3)*(P1-P2+P3)*(P2+P3-P1)).LE.0.)  GO  TO  903
  904 IF(IBA.EQ.1.AND.T.LT.4.)   CT1=COSTA(18,T)
      IF(IBA.EQ.2.AND.T.LT.4.)   CT1=COSTA(15,T)
      IF(IBA.EQ.1.AND.T.GE.4.)   CT1=COSTA(20,T)
      IF(IBA.EQ.2.AND.T.GE.4.)   CT1=COSTA(17,T)
C
       CALL  PAULBL(P2,V,IBL)
       IF(IBL.EQ.1)  GO  TO  906
       IF(AM3.EQ.0.939)  CALL  PAULBL(P3,V,IBL)
       IF(IBL.EQ.1)  GO  TO  906
C
  905 CONTINUE
      CALL  POSTR(P1,CT1,AM1,V,B,JP,SK,STOT,RL,N,IND)
  906 CONTINUE
C      write( *,*) 'out conec'
      RETURN
      END
C
C * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
      SUBROUTINE  PIIZO(U,PL,IFIZO)
      IMPLICIT REAL *8(A-H,O-Z)
      REAL *8 U,PL(3),MN,MPI,MD,MD1,MD2,MD0
      COMMON/CPIIZO/ PN1(3),PN2(3)
      DIMENSION  VD(3),PS(3),PSN(3)
      DATA MN/0.939/,MPI/0.139/,MD0/1.232/
C      write( *,*) 'in piizo'
      IFIZO=0
      MD1=1.081
      IF(U.LE.(MD1+MN)) GO  TO  999
      MD2=U-MN
      IF(MD2.LE.MD1)    GO  TO  999
      FMD0=1.
      IF(MD2.LE.MD0)  THEN
      CALL  WMD(MD2,TAU2,FMD2)
      ENDIF
  901 MD=MD1+RNDMD(-1)*(MD2-MD1)
      CALL  WMD(MD,TAU,FMD)
      IF(RNDMD(-1).GT.FMD/FMD2)  GO  TO  901
      ED=(U**2+MD**2-MN**2)/2./U
      IF(ED.LE.MD)   THEN
       write( *,*) 'ED<MD', ED,MD
       GO  TO  901
      ENDIF
      PD=DSQRT(ED**2-MD**2)
      CTD=1.-2.*RNDMD(-1)
      FID=6.283185*RNDMD(-1)
      STD=DSQRT(1.-CTD**2)
      VD(1)=PD/ED*STD*DCOS(FID)
      VD(2)=PD/ED*STD*DSIN(FID)
      VD(3)=PD/ED*CTD
      DO  k=1,3
      PN1(k)=-VD(k)*ED
      ENDDO
      ESPI=(MD**2+MPI**2-MN**2)/2./MD
      IF(ESPI.LE.MPI)   THEN
       write( *,*) 'ESPI<MPI', ESPI,MPI
       GO  TO  901
      ENDIF
      PSM=DSQRT(ESPI**2-MPI**2)
      CT=1.-2.*RNDMD(-1)
      FI=6.283185*RNDMD(-1)
      ST=DSQRT(1.-CT**2)
      PS(1)=PSM*ST*DCOS(FI)
      PS(2)=PSM*ST*DSIN(FI)
      PS(3)=PSM*CT
      CALL  KINEMA(PS,VD,PL,CT1,ST1,CF1,SF1,TL1,MPI)
      IFIZO=1
      DO  k=1,3
       PSN(k)=-PS(k)
      ENDDO
      CALL  KINEMA(PSN,VD,PN2,CTN,STN,CFN,SFN,TLN,MN)
  999 CONTINUE
C      write( *,*) 'out piizo'
      RETURN
      END
C
C * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
      SUBROUTINE SARNDT(U,SIJ)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 U,SIJ(4)
      REAL*8 MPI,MN,M0,MDEL0,MNUC0,MSDEL,MSN
      DIMENSION ALFA(4),BETA(4),M0(4),GG(4)
      DATA ALFA /3.772,15.28,146.3,6.030/
      DATA BETA/1.262,0.,0.,1.700/
      DATA M0/1188.,1245.,1472.,1203./
      DATA GG/99.02,137.4,26.49,134.3/
C      write( *,*) 'in sarndt'
      do  i=1,4
      SIJ(i)=0.
      enddo
      S=(U*1000.)**2
      PI=3.1415927
      HC=198.
      PHC=PI*HC**2
      MPI=139.0
      MN=939.0
      MNUC0=1430.
      MDEL0=1220.
      G0D=120.
      G0N=200.
      DS=DSQRT(S)
      IF(DS.LE.(2.*MN+MPI))  RETURN
      P2=S/4.-MN**2
      ZPD=(DS-MN-MDEL0)*(2./G0D)
      ZMD=(MN+MPI-MDEL0)*(2./G0D)
      MSDEL=MDEL0+G0D*DLOG((1.+ZPD**2)/(1.+ZMD**2))/
     *(4.*(DATAN(ZPD)-DATAN(ZMD)))
      ZPN=(DS-MN-MNUC0)*(2./G0N)
      ZMN=(MN+MPI-MNUC0)*(2./G0N)
      MSN=MNUC0+G0N*DLOG((1.+ZPN**2)/(1.+ZMN**2))/
     *(4.*(DATAN(ZPN)-DATAN(ZMN)))
      DO 1 I=1,3
      OM=MSDEL
      IF (I.EQ.3) OM=MSN
C
      IF(DS.LE.(MN+OM))    GO  TO  1
      IF(OM.LE.(MN+MPI))  GO  TO  1
C
      GM2=(M0(I)*GG(I))**2
      PR2=(S-(MN-OM)**2)*(S-(MN+OM)**2)/(4.*S)
      Q2S=(OM**2-(MN-MPI)**2)*(OM**2-(MN+MPI)**2)/(4.*OM**2)
      OP=M0(I)
      Q20=(OP**2-(MN-MPI)**2)*(OP**2-(MN+MPI)**2)/(4.*OP**2)
      S0=(MN+M0(I))**2
      P0=S0/4.-MN**2
      AL=ALFA(I)
      BE=BETA(I)
      SIJ(I)=PHC/(2.*P2)*AL*(PR2/P0)**(BE/2.)*GM2*(Q2S/Q20)**(3./2.)/
     *((OM**2-M0(I)**2)**2+GM2)*10.
   1  CONTINUE
c s(NN=D+pi)
      SPIN=(DS-MN)**2
      PR2=(S-(2.*MN-MPI)**2)*(S-(2.*MN+MPI)**2)/(4.*S)
      S0=(MN+M0(4))**2
      P02=S0/4.-MN**2
      GM2=(M0(4)*GG(4))**2
      AL=ALFA(4)
      BE=BETA(4)
      SIJ(4)=PHC/(2.*P2)*AL*(PR2/P02)**(BE/2.)*GM2/
     *((SPIN-M0(4)**2)**2+GM2)*10.
C      write( *,*) 'out sarndt sij=',SIJ
      RETURN
      END
C                 
C * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
      SUBROUTINE  PAULBL(P,V,IBL)
      IMPLICIT REAL *8(A-H,O-Z)
      REAL *8 P,V(3)
      DIMENSION  PS(3),PL(3)
      IBL=0
      PF=0.250
      CTS=1.-2.*RNDMD(-1)
      FIS=2.*RNDMD(-1)
      STS=DSQRT(1.-CTS**2)
      PS(1)=P*STS*DCOS(FIS)
      PS(2)=P*STS*DSIN(FIS)
      PS(3)=P*CTS
      AM=0.939
      CALL  KINEMA(PL,V,PS,CTL,STL,CFL,SFL,TL,AM)
      PLM=DSQRT(TL*(TL+2.*AM))
      IF(PLM.LE.PF)  IBL=1
      RETURN
      END
C
C * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
      FUNCTION  SZWER(JP,U)
C    B+B==>B+Y+K   CROSS SECTION  (ZWERMANN APPROXIMATION)
      IMPLICIT REAL *8(A-H,O-Z)
      REAL *8 SZWER,U
      DIMENSION AM(8)
      DATA  AM/2*0.494,1.115,3*1.192,0.494,0.939/
      AMJ=AM(JP)
      AMK=0.494
      SZWER=0.
      IF(JP.LE.2) AMX=1.115+0.939
      IF(JP.GT.2) AMX=0.494+0.939
      IF(U.LE.(AMX+AMJ))  RETURN
      PMAX=DSQRT((U**2-(AMX+AMJ)**2)*(U**2-(AMX-AMJ)**2))/(2.*U)
      SZWER=0.049*(PMAX/AMK)**4
      RETURN
      END
C
C * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
      FUNCTION  SST(IB1,IB2,IE1,IE2,JP,T,IS)
C TOTAL PRODUCTION CROSS SECTION OF JP
      IMPLICIT REAL *8(A-H,O-Z)
      REAL *8 SST,T
      SST=0.
      GO  TO  (101,102,103,104,105,106,107),JP
  101 IF((IB1+IB2).EQ.1)  GO  TO  51
      IF((IE1+IE2).EQ.2)  GO  TO  31
      IF((IE1+IE2).EQ.1)  GO  TO  32
      IF((IE1+IE2).EQ.0)  SST=SIG(T,13)
      GO  TO  200
   31 S1=SIG(T,1)
      S2=SIG(T,5)
      S3=SIG(T,7)
      IF(IS.EQ.1) THEN
      RS=RNDMD(-1)*(S1+S2+S3)
                              SST=S1
      IF(RS.GT.S1.AND.RS.LE.(S1+S2))  SST=S2
      IF(RS.GT.(S1+S2))               SST=S3
      ELSE
C !!!
      SST=S1+S2+S3
C !!!
      ENDIF
      GO  TO  200
   32 S1=SIG(T,2)
      S2=SIG(T,9)
      S3=SIG(T,11)
      IF(IS.EQ.1) THEN

      RS=RNDMD(-1)*(S1+S2+S3)
                              SST=S1
      IF(RS.GT.S1.AND.RS.LE.(S1+S2))  SST=S2
      IF(RS.GT.(S1+S2))               SST=S3
      ELSE
C !!!
      SST=S1+S2+S3
C !!!
      ENDIF
      GO  TO  200
   51 IF(IE1.EQ.1.AND.IE2.EQ.1)   SST=SIG(T,15)
      IF(IE1.EQ.-1.AND.IE2.EQ.1)  SST=SIG(T,21)
      IF(IE1.EQ.0.AND.IE2.EQ.1)   GO  TO  21
      IF(IE1.EQ.1.AND.IE2.EQ.0)   GO  TO  22
      IF(IE1.EQ.0.AND.IE2.EQ.0)   SST=SIG(T,27)
      GO  TO  200
   21 S1=SIG(T,23)
      S2=SIG(T,25)
      IF(IS.EQ.1) THEN
      RS=RNDMD(-1)*(S1+S2)
                SST=S1
      IF(RS.GT.S1)  SST=S2
      ELSE
C !!!
      SST=S1+S2
C !!!
      ENDIF
      GO  TO  200
   22 S1=SIG(T,17)
      S2=SIG(T,18)
      IF(IS.EQ.1) THEN
      RS=RNDMD(-1)*(S1+S2)
                SST=S1
      IF(RS.GT.S1)  SST=S2
      ELSE
C !!!
      SST=S1+S2
C !!!
      ENDIF
      GO  TO  200
  102 IF((IB1+IB2).EQ.1)   GO  TO  52
      IF((IE1+IE2).EQ.2)   SST=SIG(T,6)
      IF((IE1+IE2).EQ.1)   GO  TO  33
      IF((IE1+IE2).EQ.0)   GO  TO  34
      GO  TO  200
   33 S1=SIG(T,3)
      S2=SIG(T,8)
      S3=SIG(T,10)
      IF(IS.EQ.1) THEN
      RS=RNDMD(-1)*(S1+S2+S3)
                              SST=S1
      IF(RS.GT.S1.AND.RS.LE.(S1+S2))  SST=S2
      IF(RS.GT.(S1+S2))               SST=S3
      ELSE
C !!!
      SST=S1+S2+S3
C !!!
      ENDIF
      GO  TO  200
   34 S1=SIG(T,4)
      S2=SIG(T,14)
      S3=SIG(T,12)
      IF(IS.EQ.1) THEN
      RS=RNDMD(-1)*(S1+S2+S3)
                              SST=S1
      IF(RS.GT.S1.AND.RS.LE.(S1+S2))  SST=S2
      IF(RS.GT.(S1+S2))               SST=S3
      ELSE
C !!!
      SST=S1+S2+S3
C !!!
      ENDIF
      GO  TO  200
   52 IF(IE1.EQ.-1.AND.IE2.EQ.1)  GO  TO  23
      IF(IE1.EQ.0.AND.IE2.EQ.1)   SST=SIG(T,24)
      IF(IE1.EQ.1.AND.IE2.EQ.0)   SST=SIG(T,16)
      IF(IE1.EQ.-1.AND.IE2.EQ.0)  SST=SIG(T,22)
      IF(IE1.EQ.0.AND.IE2.EQ.0)   GO  TO  24
      GO  TO  200
   23 S1=SIG(T,19)
      S2=SIG(T,20)
      IF(IS.EQ.1) THEN
      RS=RNDMD(-1)*(S1+S2)
                SST=S1
      IF(RS.GT.S1)  SST=S2
      ELSE
C !!!
      SST=S1+S2
C !!!
      ENDIF
      GO  TO  200
   24 S1=SIG(T,28)
      S2=SIG(T,26)
      IF(IS.EQ.1) THEN
      RS=RNDMD(-1)*(S1+S2)
                SST=S1
      IF(RS.GT.S1)  SST=S2
      ELSE
C !!!
      SST=S1+S2
C !!!
      ENDIF
      GO  TO  200
  103 IF((IB1+IB2).EQ.1)   GO  TO  53
      IF((IE1+IE2).EQ.2)   SST=SIG(T,1)
      IF((IE1+IE2).EQ.1)   GO  TO  25
      IF((IE1+IE2).EQ.0)   SST=SIG(T,4)
      GO  TO  200
   25 S1=SIG(T,2)
      S2=SIG(T,3)
      IF(IS.EQ.1) THEN
      RS=RNDMD(-1)*(S1+S2)
                SST=S1
      IF(RS.GT.S1)  SST=S2
      ELSE
C !!!
      SST=S1+S2
C !!!
      ENDIF
      GO  TO  200
   53 IF(IE1.EQ.-1.AND.IE2.EQ.1)  SST=SIG(T,19)
      IF(IE1.EQ.0.AND.IE2.EQ.1)   SST=SIG(T,23)
      IF(IE1.EQ.1.AND.IE2.EQ.0)   SST=SIG(T,18)
      IF(IE1.EQ.0.AND.IE2.EQ.0)   SST=SIG(T,26)
      GO  TO  200
  104 IF((IB1+IB2).EQ.1)   GO  TO  54
      IF((IE1+IE2).EQ.2)   GO  TO  26
      IF((IE1+IE2).EQ.1)   SST=SIG(T,8)
      GO  TO  200
   26 S1=SIG(T,5)
      S2=SIG(T,6)
      IF(IS.EQ.1) THEN
      RS=RNDMD(-1)*(S1+S2)
                SST=S1
      IF(RS.GT.S1)  SST=S2
      ELSE
C !!!
      SST=S1+S2
C !!!
      ENDIF
      GO  TO  200
   54 IF(IE1.EQ.1.AND.IE2.EQ.1)   SST=SIG(T,15)
      IF(IE1.EQ.0.AND.IE2.EQ.1)   SST=SIG(T,24)
      IF(IE1.EQ.1.AND.IE2.EQ.0)   SST=SIG(T,16)
      GO  TO  200
  105 IF((IB1+IB2).EQ.1)   GO  TO  55
      IF((IE1+IE2).EQ.1)   SST=SIG(T,11)
      IF((IE1+IE2).EQ.0)   GO  TO  27
      GO  TO  200
   27 S1=SIG(T,13)
      S2=SIG(T,14)
      IF(IS.EQ.1) THEN
      RS=RNDMD(-1)*(S1+S2)
                SST=S1
      IF(RS.GT.S1)  SST=S2
      ELSE
C !!!
      SST=S1+S2
C !!!
      ENDIF
      GO  TO  200
   55 IF(IE1.EQ.-1.AND.IE2.EQ.1)  SST=SIG(T,21)
      IF(IE1.EQ.-1.AND.IE2.EQ.0)  SST=SIG(T,22)
      IF(IE1.EQ.0.AND.IE2.EQ.0)   SST=SIG(T,27)
      GO  TO  200
  106 IF((IB1+IB2).EQ.1)   GO  TO  56
      IF((IE1+IE2).EQ.2)   SST=SIG(T,7)
      IF((IE1+IE2).EQ.1)   GO  TO  28
      IF((IE1+IE2).EQ.0)   SST=SIG(T,12)
      GO  TO  200
   28 S1=SIG(T,9)
      S2=SIG(T,10)
      IF(IS.EQ.1) THEN
      RS=RNDMD(-1)*(S1+S2)
                SST=S1
      IF(RS.GT.S1)  SST=S2
      ELSE
C !!!
      SST=S1+S2
C !!!
      ENDIF
      GO  TO  200
   56 IF(IE1.EQ.-1.AND.IE2.EQ.1)  SST=SIG(T,20)
      IF(IE1.EQ.0.AND.IE2.EQ.1)   SST=SIG(T,25)
      IF(IE1.EQ.1.AND.IE2.EQ.0)   SST=SIG(T,17)
      IF(IE1.EQ.0.AND.IE2.EQ.0)   SST=SIG(T,28)
      GO  TO  200
  107 CONTINUE
  200 CONTINUE
      RETURN
      END
C
C * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
      FUNCTION   SSL(IB1,IB2,IE1,IE2,JP,T)
C TOTAL PRODUCTION CROSS SECTION OF JP+L(AS A PARTENER)
      IMPLICIT REAL *8(A-H,O-Z)
      REAL *8 SSL,T
      SSL=0.
      GO  TO  (201,202,103,300,300,300,300),JP
  103 SSL=SST(IB1,IB2,IE1,IE2,JP,T,0)
      GO  TO  300
  201 IF((IB1+IB2).EQ.1)  GO  TO  61
      IF((IE1+IE2).EQ.2)   SSL=SIG(T,1)
      IF((IE1+IE2).EQ.1)   SSL=SIG(T,2)
      GO  TO  300
   61 IF(IE1.EQ.0.AND.IE2.EQ.1)   SSL=SIG(T,23)
      IF(IE1.EQ.1.AND.IE2.EQ.0)   SSL=SIG(T,18)
      GO  TO  300
  202 IF((IB1+IB2).EQ.1)   GO  TO  62
      IF((IE1+IE2).EQ.1)   SSL=SIG(T,3)
      IF((IE1+IE2).EQ.0)   SSL=SIG(T,4)
      GO  TO  300
   62 IF(IE1.EQ.-1.AND.IE2.EQ.1)  SSL=SIG(T,19)
      IF(IE1.EQ.0.AND.IE2.EQ.0)   SSL=SIG(T,26)
  300 CONTINUE
      RETURN
      END
C
C * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
      FUNCTION  SIG(T,N)
      IMPLICIT REAL *8(A-H,O-Z)
      REAL *8 SIG,T
C                        N                  REACTION
C                        1*             p+p==>p + L  + K+
C                        2              p+n==>n + L  + K+
C                        3              p+n==>p + L  + K0
C                        4              n+n==>n + L  + K0
C                        5*             p+p==>n + S+ + K+
C                        6*             p+p==>p + S+ + K0
C                        7*             p+p==>p + S0 + K+
C                        8              p+n==>n + S+ + K0
C                        9              p+n==>n + S0 + K+
C                       10              p+n==>p + S0 + K0
C                       11              p+n==>p + S- + K+
C                       12              n+n==>n + S0 + K0
C                       13              n+n==>n + S- + K+
C                       14              n+n==>p + S- + K0
C                       15*             (pi+)+p==>S+ + K+
C                       16              (pi+)+n==>S+ + K0
C                       17              (pi+)+n==>S0 + K+
C                       18              (pi+)+n==>L  + K+
C                       19*             (pi-)+p==>L  + K0
C                       20*             (pi-)+p==>S0 + K0
C                       21*             (pi-)+p==>S- + K+
C                       22              (pi-)+n==>S- + K0
C                       23              (pi0)+p==>L  + K+
C                       24              (pi0)+p==>S+ + K0
C                       25              (pi0)+p==>S0 + K+
C                       26              (pi0)+n==>L  + K0
C                       27              (pi0)+n==>S- + K+
C                       28              (pi0)+n==>S0 + K0
C
C   PP   --> L K+ P
      SN1(X)=FAR(0.132D0,0.709D0,3.432D0,-0.564D0,X)
C   PP   --> S+ K+ N
      SN5(X)=FAR(0.144D0,0.709D0,3.432D0,-0.564D0,X)
C   PP   --> S+ K0 P
      SN6(X)=FAR(0.062D0,0.709D0,3.432D0,-0.564D0,X)
C   PP   --> S0 K+ P
      SN7(X)=FAR(0.058D0,0.709D0,3.432D0,-0.564D0,X)
C   PI+P --> K+ S+
      S15(X)=FAR(0.565D0,1.129D0,0.662D0,-1.539D0,X)
C   PI-P --> K0 L
      S19(X)=FAR(0.183D0,1.375D0,0.130D0,-1.213D0,X)
C   PI-P --> K0 S0
      S20(X)=FAR(0.098D0,1.221D0,0.150D0,-1.073D0,X)
C   PI-P --> K+ S-
      S21(X)=FAR(0.112D0,0.873D0,0.457D0,-1.724D0,X)
C
      TL2=T-0.759
      TS2=T-0.888
      TL3=T-1.579
      TS3=T-1.780

c      TL2=T-0.760
c      TS2=T-0.900
c      TL3=T-1.570
c      TS3=T-1.800
      GO  TO  (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,
     *21,22,23,24,25,26,27,28),N
    1 SIG=SN1(TL3)
      RETURN
    2 SIG=5.*SN1(TL3)
      RETURN
    3 GO  TO  2
    4 SIG=1.*SN1(TL3)
      RETURN
    5 SIG=SN5(TS3)
      RETURN
    6 SIG=SN6(TS3)
      RETURN
    7 SIG=SN7(TS3)
      RETURN
    8 SIG=1./3.*SN5(TS3)-2.*SN7(TS3)+4.*SN6(TS3)
      RETURN
    9 SIG=4./3.*SN5(TS3)-7.*SN7(TS3)+6.*SN6(TS3)
      RETURN
   10 GO  TO  9
   11 SIG=9.*SN5(TS3)-8.*SN7(TS3)+9.*SN6(TS3)
      RETURN
   12 GO  TO  7
   13 SIG=2./3.*SN5(TS3)-4.*SN7(TS3)+3.*SN6(TS3)
   14 GO  TO  5
   15 SIG=S15(TS2)
      RETURN
   16 GO  TO  21
   17 GO  TO  20
   18 GO  TO  19
   19 SIG=S19(TL2)
      RETURN
   20 SIG=S20(TS2)
      RETURN
   21 SIG=S21(TS2)
      RETURN
   22 GO  TO  15
   23 SIG=S19(TL2)/2.
      RETURN
   24 SIG=2./3.*(S21(TS2)+S15(TS2)/3.-S20(TS2)/2.)
      RETURN
   25 SIG=1./2.*(S21(TS2)+S15(TS2)-S20(TS2))
      RETURN
   26 GO  TO  23
   27 GO  TO  20
   28 SIG=1./2.*(S21(TS2)+S15(TS2)-S20(TS2))
      RETURN
      END
C
C * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
      FUNCTION  FAR(A1,A2,A3,A4,X)
      IMPLICIT REAL *8(A-H,O-Z)
      REAL *8 FAR,A1,A2,A3,A4,X
      FAR=0.0
      IF(X.LE.0.)   RETURN
      S=A3**2+X**2
      FAR=A1*(X**A2)*(S**A4)
      RETURN
      END
C
C * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
      FUNCTION SPINET(T)
      IMPLICIT REAL *8(A-H,O-Z)
      REAL*8  SPINET,T,mpi,mn,meta
C      (pi-)+p=>eta+n   CROSS SECTIONS  [mb]
c     T - lab. kinetic energy of pion
      DATA mpi/0.140/,mn/0.940/,meta/0.549/
      ss0=mn+meta
      E=T+mpi
      P=DSQRT(T*(T+2.*mpi))
      ss=DSQRT((E+mn)**2-P**2)
      SPINET=0.
      if(ss.le.ss0) RETURN
      x=ss-ss0
      SPINET=FAR(0.485D0,0.643D0,0.039D0,-0.649D0,x)
      RETURN
      END
C
C * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
      FUNCTION SPPETA(T)
      IMPLICIT REAL *8(A-H,O-Z)
      REAL*8  SPPETA,T,mn,meta
C      p+p=>eta+p+p   CROSS SECTIONS  [mb]
c  De Paoli et al. Pys.Lett. B219(1989),194
c     T - lab. kinetic energy of proton
      DATA mn/0.940/,meta/0.549/,A/0.17/,B/0.253/
      ss0=2.*mn+meta
      E=T+mn
      P=DSQRT(T*(T+2.*mn))
      ss=DSQRT((E+mn)**2-P**2)
      SPPETA=0.
      if(ss.le.ss0) RETURN
      x=ss-ss0
      SPPETA=A*x/(B+x**2)
      RETURN
      END
C
C * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
      FUNCTION SCUGNP(IS,T)
      IMPLICIT REAL *8(A-H,O-Z)
      REAL*8  SCUGNP,T
C      ISOSPIN AVERAGED CROSS SECTIONS
C      from J.Cugnon et al. P.R.41C(1990)1701
C             IS=0  PI+N=>L+K
C             IS=1  PI+N=>S+K
      SCUGNP=0.
      IF(IS.EQ.0) THEN
C S( PI N=>L K )=0.5*S(PI- P=>L K0 )
       SCUGNP=0.5*SIC(T,19)
       GO  TO  10
      ELSE
C S( PI N=>S K )=0.5*(S(PI+ P=>S+ K+)+S(PI- P=>S- K+)+S(PI- P=>S0 K0))
       SCUGNP=0.5*(SIC(T,15)+SIC(T,21)+SIC(T,20))
       GO  TO  10
      ENDIF
   10 RETURN
      END
C
C * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
      FUNCTION  SCT(IE1,IE2,JP,T,IS)
C TOTAL PRODUCTION CROSS SECTION OF JP FROM PI+N=>JP+Y
C      from J.Cugnon et al. P.R.41C(1990)1701
      IMPLICIT REAL *8(A-H,O-Z)
      REAL *8 SCT,T
      SCT=0.
      GO  TO  (101,102,103,104,105,106,107),JP
  101 IF(IE1.EQ.1.AND.IE2.EQ.1)   SCT=SIC(T,15)
      IF(IE1.EQ.-1.AND.IE2.EQ.1)  SCT=SIC(T,21)
      IF(IE1.EQ.0.AND.IE2.EQ.1)   GO  TO  21
      IF(IE1.EQ.1.AND.IE2.EQ.0)   GO  TO  22
      IF(IE1.EQ.0.AND.IE2.EQ.0)   SCT=SIC(T,27)
      GO  TO  200
   21 S1=SIC(T,23)
      S2=SIC(T,25)
      IF(IS.EQ.1) THEN
      RS=RNDMD(-1)*(S1+S2)
                SCT=S1
      IF(RS.GT.S1)  SCT=S2
      ELSE
C !!!
      SCT=S1+S2
C !!!
      ENDIF
      GO  TO  200
   22 S1=SIC(T,17)
      S2=SIC(T,18)
      IF(IS.EQ.1) THEN
      RS=RNDMD(-1)*(S1+S2)
                SCT=S1
      IF(RS.GT.S1)  SCT=S2
      ELSE
C !!!
      SCT=S1+S2
C !!!
      ENDIF
      GO  TO  200
  102 IF(IE1.EQ.-1.AND.IE2.EQ.1)  GO  TO  23
      IF(IE1.EQ.0.AND.IE2.EQ.1)   SCT=SIC(T,24)
      IF(IE1.EQ.1.AND.IE2.EQ.0)   SCT=SIC(T,16)
      IF(IE1.EQ.-1.AND.IE2.EQ.0)  SCT=SIC(T,22)
      IF(IE1.EQ.0.AND.IE2.EQ.0)   GO  TO  24
      GO  TO  200
   23 S1=SIC(T,19)
      S2=SIC(T,20)
      IF(IS.EQ.1) THEN
      RS=RNDMD(-1)*(S1+S2)
                SCT=S1
      IF(RS.GT.S1)  SCT=S2
      ELSE
C !!!
      SCT=S1+S2
C !!!
      ENDIF
      GO  TO  200
   24 S1=SIC(T,28)
      S2=SIC(T,26)
      IF(IS.EQ.1) THEN
      RS=RNDMD(-1)*(S1+S2)
                SCT=S1
      IF(RS.GT.S1)  SCT=S2
      ELSE
C !!!
      SCT=S1+S2
C !!!
      ENDIF
      GO  TO  200
  103 IF(IE1.EQ.-1.AND.IE2.EQ.1)  SCT=SIC(T,19)
      IF(IE1.EQ.0.AND.IE2.EQ.1)   SCT=SIC(T,23)
      IF(IE1.EQ.1.AND.IE2.EQ.0)   SCT=SIC(T,18)
      IF(IE1.EQ.0.AND.IE2.EQ.0)   SCT=SIC(T,26)
      GO  TO  200
  104 IF(IE1.EQ.1.AND.IE2.EQ.1)   SCT=SIC(T,15)
      IF(IE1.EQ.0.AND.IE2.EQ.1)   SCT=SIC(T,24)
      IF(IE1.EQ.1.AND.IE2.EQ.0)   SCT=SIC(T,16)
      GO  TO  200
  105 IF(IE1.EQ.-1.AND.IE2.EQ.1)  SCT=SIC(T,21)
      IF(IE1.EQ.-1.AND.IE2.EQ.0)  SCT=SIC(T,22)
      IF(IE1.EQ.0.AND.IE2.EQ.0)   SCT=SIC(T,27)
      GO  TO  200
  106 IF(IE1.EQ.-1.AND.IE2.EQ.1)  SCT=SIC(T,20)
      IF(IE1.EQ.0.AND.IE2.EQ.1)   SCT=SIC(T,25)
      IF(IE1.EQ.1.AND.IE2.EQ.0)   SCT=SIC(T,17)
      IF(IE1.EQ.0.AND.IE2.EQ.0)   SCT=SIC(T,28)
      GO  TO  200
  107 CONTINUE
  200 CONTINUE
      RETURN
      END
C
C * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
      FUNCTION   SCL(IE1,IE2,JP,T)
C TOTAL PRODUCTION CROSS SECTION OF JP+L(AS A PARTNER)
      IMPLICIT REAL *8(A-H,O-Z)
      REAL *8 SCL,T
      SCL=0.
      GO  TO  (201,202,103,300,300,300,300),JP
  103 SCL=SCT(IE1,IE2,JP,T,0)
      GO  TO  300
  201 IF(IE1.EQ.0.AND.IE2.EQ.1)   SCL=SIC(T,23)
      IF(IE1.EQ.1.AND.IE2.EQ.0)   SCL=SIC(T,18)
      GO  TO  300
  202 IF(IE1.EQ.-1.AND.IE2.EQ.1)  SCL=SIC(T,19)
      IF(IE1.EQ.0.AND.IE2.EQ.0)   SCL=SIC(T,26)
  300 CONTINUE
      RETURN
      END
C
C * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
      FUNCTION  SIC(T,N)
      IMPLICIT REAL *8(A-H,O-Z)
      REAL *8 SIC,T
C
      P=DSQRT(T*(T+0.28))
      N1=N-14
      GO  TO  (15,16,17,18,19,20,21,22,23,24,25,26,27,28),N1
   15 SIC=SC(15,P)
      RETURN
   16 GO  TO  21
   17 GO  TO  20
   18 GO  TO  19
   19 SIC=SC(19,P)
      RETURN
   20 SIC=SC(20,P)
      RETURN
   21 SIC=SC(21,P)
      RETURN
   22 GO  TO  15
   23 SIC=SC(19,P)/2.
      RETURN
   24 SIC=2./3.*(SC(21,P)+SC(15,P)/3.-SC(20,P)/2.)
      RETURN
   25 SIC=1./2.*(SC(21,P)+SC(15,P)-SC(20,P))
      RETURN
   26 GO  TO  23
   27 GO  TO  20
   28 SIC=1./2.*(SC(21,P)+SC(15,P)-SC(20,P))
      RETURN
      END
C
C * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
      FUNCTION SC(N,P)
      IMPLICIT REAL *8(A-H,O-Z)
      COMMON /KSIGS/ KSIGS
      REAL*8  SC,P,MPI,MN,MK,ML,MS
      DATA MPI/0.139/,MN/0.939/,MK/0.494/,ML/1.115/,MS/1.192/
C     APPR. CR.SEC. PI+N=>K+Y from J.Cugnon et al.
C     KSIGS=3  CUGNON (N.P.A422,p.635,1984)
C     KSIGS=4  CUGNON (P.R.C41,P.1701,1990)
      SC=0.
      IF(KSIGS.EQ.4) GO  TO  4
      IF(N.EQ.19) THEN
C PI- P=>L K0
       IF(P.LT.0.9)             SC=0.
       IF(P.GE.0.9.AND.P.LT.1.) SC=0.65*P**(42.00)
       IF(P.GE.1.)              SC=0.65*P**(-1.67)
      ELSE IF(N.EQ.15) THEN
C PI+ P=>S+ K+
       IF(P.LT.1.05)              SC=0.
       IF(P.GE.1.05.AND.P.LT.1.5) SC=0.10*P**(4.80)
       IF(P.GE.1.5)               SC=1.48*P**(-1.85)
      ELSE IF(N.EQ.21) THEN
C PI- P=>S- K+
       IF(P.LT.1.03)              SC=0.
       IF(P.GE.1.03.AND.P.LT.1.5) SC=0.12*P**(1.60)
       IF(P.GE.1.5)               SC=1.12*P**(-3.60)
      ELSE IF(N.EQ.20) THEN
C PI- P=>S0 K0
       IF(P.LT.1.03)              SC=0.
       IF(P.GE.1.03)              SC=0.31*P**(-1.53)
      ELSE
       RETURN
      ENDIF
      RETURN
C
    4 CONTINUE
      SS=DSQRT((DSQRT(P**2+MPI**2)+MN)**2-P**2)
      IF(N.EQ.19) THEN
C PI- P=>L K0
       SS0=MK+ML
       IF(SS.LT.SS0)               SC=0.
       IF(SS.GE.SS0.AND.SS.LT.1.7) SC=0.9*(SS-SS0)/0.091
       IF(SS.GT.1.7)               SC=0.09/(SS-1.6)
      ELSE IF(N.EQ.15) THEN
       SS0=MK+MS
C PI+ P=>S+ K+
       IF(SS.LT.SS0)               SC=0.
       IF(SS.GE.SS0.AND.SS.LT.1.9) SC=0.7*(SS-SS0)/0.218
       IF(SS.GE.1.9)               SC=0.14/(SS-1.7)
      ELSE IF(N.EQ.21.OR.N.EQ.20) THEN
C PI- P=>S- K+  OR PI- P=>S0 K0
       SS0=MK+MS
       IF(SS.LT.SS0)               SC=0.
       IF(SS.GE.SS0)               SC=0.25*(1.-0.75*(SS-1.682))
       IF(SC.LE.0.)                SC=0.
      ELSE
       RETURN
      ENDIF
      RETURN
      END
C
C * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
      FUNCTION STSU(N,SS)
      IMPLICIT REAL *8(A-H,O-Z)
      REAL*8  STSU,SS,MPI,MN,MK,ML,MS
      DATA MK/0.494/,ML/1.115/,MS/1.192/
C     APPR. CR.SEC. PI+Delta=>K+Y
C     from K.Tsushima et al. J.Phys.G:Nucl.Part.Phys.21(1995)33-42
      STSU=0.
      SS0S=MK+MS
      SS0L=MK+ML
      IF(N.EQ.1) THEN
C PI- D++=>S0 K+
       IF(SS.LE.SS0S)   RETURN
       STSU=0.004959*((SS-SS0S)**0.7785)/((SS-1.725)**2+0.008147)
      ELSEIF(N.EQ.2) THEN
C PI0 D0 =>S- K+
       IF(SS.LE.SS0S)   RETURN
       STSU=0.006964*((SS-SS0S)**0.8140)/((SS-1.725)**2+0.007713)
      ELSEIF(N.EQ.3) THEN
C PI+ D0 =>S0 K+
       IF(SS.LE.SS0S)   RETURN
       STSU=0.002053*((SS-SS0S)**0.9853)/((SS-1.725)**2+0.005414)+
     +      0.317900*((SS-SS0S)**0.9025)/((SS-2.675)**2+44.88)
      ELSEIF(N.EQ.4) THEN
C PI+ D- =>S- K+
       IF(SS.LE.SS0S)   RETURN
       STSU=0.017410*((SS-SS0S)**1.2078)/((SS-1.725)**2+0.003777)
      ELSEIF(N.EQ.5) THEN
C PI- D++=>L  K+
       IF(SS.LE.SS0L)   RETURN
       STSU=0.006545*((SS-SS0L)**0.7866)/((SS-1.720)**2+0.004852)
      ELSE
       RETURN
      ENDIF
      RETURN
      END
C
C * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
      SUBROUTINE  SIGSTR(JP,T12,S)
      IMPLICIT REAL *8(A-H,O-Z)
      REAL *8 T12,S
      COMMON/ISRES/ISRES
      S=0.
      IF(ISRES.EQ.0)  THEN
C  S(TOT)  K+ OR K0  ON  NUCLEON
      IF(JP.LE.2)  S=13.
C  S(TOT)  L OR S    ON  NUCLEON
      IF(JP.GT.2)  S=30.
C  S(TOT)    K-      ON  NUCLEON
      IF(JP.EQ.7)  S=40.
C  S(TOT)  ANTI-P    ON  NUCLEON
      IF(JP.EQ.8)  S=94.65*T12**(-0.2943)
      ELSE
C  S(TOT)  ETA,ETAP,RHO,OMEGA,PHI  ON  NUCLEON
      IF(JP.LE.5)  S=20.
C  S(TOT)  ANTI-P    ON  NUCLEON
      IF(JP.EQ.6)  S=94.65*T12**(-0.2943)
      ENDIF
      RETURN
      END
C
C * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
      SUBROUTINE RESCA1(PL,AM,JP,RL0,CTL,STL,CFL,SFL,TL,W,WA,TAU)
      IMPLICIT REAL *8(A-H,O-Z)
      REAL *8 PL,AM,RL0,CTL,STL,CFL,SFL,TL,W,WA,TAU,MMES,MBAR
      COMMON/NUCSP/VPR(3),VTA(3),RADP(3),RADT(3),VEV(3),VRE(3),
     *GEV,GRE,VEP(3),VET(3),GEP,GET
      COMMON/NRESC/RESC(2),STRIN
      COMMON/RESULT/AN1,AN2,ZN1,ZN2,EN1,EN2,PN1(3),PN2(3),AM1(3),AM2(3)
      COMMON/MMATUR/MMES,MBAR
      COMMON/ISRES/ISRES
      DIMENSION  PL(3),RL(3),DT(2),VT(3),B(3),P2(3),
     *P1(3),V12(3),PIN(9),PS1(3),RL0(3)
      WA=0.
      IRESC=1 
      DO  9  K=1,3
      IF(ISRES.EQ.0.AND.JP.EQ.8) THEN
      TMAT=1./(5.06*MMES)*(1.+TL/AM)
      RL(K)=RL0(K)+TMAT*PL(K)/(TL+AM)
      ELSE
      RL(K)=RL0(K)
      ENDIF
    9 CONTINUE
      STRIN=STRIN+1.
      TAU=0.
   10 CONTINUE
      CALL  WTAR1(RL,JP,PL,AM,DT,TAU)
      IF(DT(1).LT.0..AND.DT(2).LT.0.)    GO  TO  15
      IF(DT(1).GT.0..AND.DT(2).LT.0.)    NU=1
      IF(DT(1).LT.0..AND.DT(2).GT.0.)   NU=2
      IF(DT(1).GT.0..AND.DT(2).GT.0..AND.DT(2).GT.DT(1))   NU=1
      IF(DT(1).GT.0..AND.DT(2).GT.0..AND.DT(1).GT.DT(2))   NU=2
      DO  11  K=1,3
      IF(NU.EQ.1)  VT(K)=VPR(K)
      IF(NU.EQ.2)  VT(K)=VTA(K)
      B(K)=-VT(K)
   11 CONTINUE
      PM=0.250*(RNDMD(-1))**0.333333
      CT2=1.-2.*RNDMD(-1)
      ST2=DSQRT(1.-CT2**2)
      FI=6.283185*RNDMD(-1)
      P2(1)=PM*ST2*DCOS(FI)
      P2(2)=PM*ST2*DSIN(FI)
      P2(3)=PM*CT2
      T2=DSQRT(PM**2+0.94**2)-0.94
      CALL  KINEMA(PL,B,P1,CT1,ST1,CF1,SF1,T1,AM)
      DO  12  K=1,3
   12 V12(K)=(P1(K)+P2(K))/(T1+AM+T2+0.94)
      U=(T1+AM+T2+0.94)*DSQRT(1.-V12(1)**2-V12(2)**2-V12(3)**2)
      DTAUL=DT(NU)
      DO  13  K=1,3
   13 RL(K)=RL(K)+DTAUL*PL(K)/(TL+AM)
      TAU=TAU+DTAUL
      CTS=1.-2.*RNDMD(-1)
      FIS=0.
      PIN(4)=P1(1)
      PIN(5)=P1(2)
      PIN(6)=P1(3)
      PIN(9)=AM
      AMN=0.94
      CALL  ABEL(PIN,V12,U,P1,P2,CTS,FIS,AM,AMN)
      CALL  KINEMA(P1,V12,PS1,CTS,STS,CFS,SFS,TS1,AM)
      CALL  KINEMA(PS1,VT,PL,CTL,STL,CFL,SFL,TL,AM)
      RESC(NU)=RESC(NU)+1.
C !!!!    Kostja  05.04.93+24.05.93
      IF(ISRES.EQ.0.AND.JP.EQ.7)  THEN
C   K-   ABSORPTION
      PML=DSQRT(TL*(TL+2.*AM))
      IF(NU.EQ.1) RZA=ZN1/AN1
      IF(NU.EQ.2) RZA=ZN2/AN2
      RND=RNDMD(-1)
      IF(RND.LT.RZA) SIGT=QINTG(PML,33)
      IF(RND.LT.RZA) SIGE=QINTG(PML,34)
      IF(RND.GE.RZA) SIGT=QINTG(PML,31)
      IF(RND.GE.RZA) SIGE=QINTG(PML,32)
      W=W*SIGE/SIGT
      IF(IRESC.EQ.1) WA=1.-SIGE/SIGT
      IRESC=IRESC+1
      ENDIF
C !!!!
      IF(JP.NE.8)  GO  TO  14
      IF(ISRES.EQ.0.AND.JP.EQ.8)  THEN
C ANTI-P ABSORPTION
      T12=(U**2-(AM+0.94)**2)/1.88
      SIGT=94.65*T12**(-0.2943)
      SANN=45.97*T12**(-0.3742)
C !!!!  comment for excluding the antiproton absorption
      W=W*(SIGT-SANN)/SIGT
      IF(IRESC.EQ.1) WA=SANN/SIGT
      IRESC=IRESC+1
C !!!!
      ENDIF
C  ABSORPTION FOR ALL MESONIC RESONANCES:  SIGEL/SIGTOT=0.333333
      IF(ISRES.EQ.1.AND.JP.LE.5)  THEN
      W=W*0.333333
      IF(IRESC.EQ.1) WA=0.666667
      IRESC=IRESC+1
      ENDIF
      GO  TO  10
   14 CONTINUE
      IF(JP.LT.9)  GO  TO  10
C PI-0+ ABSORPTION
      T12=(U**2-(AM+0.94)**2)/1.88
      CALL  SIGPI(JP,T12,SIGT,0)
      CALL  SIGPI(JP,T12,SIGE,1)
      CALL  SIGPI(JP,T12,SIGA,3)
      W=W*SIGE/(SIGT+SIGA)
      IF(IRESC.EQ.1) WA=SIGA/(SIGT+SIGA)
      IRESC=IRESC+1
C       WRITE( *,*) 'T12,SIGT,SIGE,SIGA=',T12,SIGT,SIGE,SIGA
      GO  TO  10
   15 CONTINUE
      RETURN
      END
C
C * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
      SUBROUTINE SIGPI(JP,T,SIG,IS)
      IMPLICIT REAL *8(A-H,O-Z)
      REAL*8 T,SIG
      COMMON/PIDABS/PIDABS
      DIMENSION IEPI(3),KSIP(3),KSIN(3)
      DATA IEPI/-1,0,1/,KSIP/2,2,1/,KSIN/1,3,3/
      JPI=JP-8
      IEP=IEPI(JPI)
      KSI=KSIN(JPI)
      IF(RNDMD(-1).GT.0.5)  KSI=KSIP(JPI)
      SIG=SIGMAG(0,0,1,KSI,IS,T)
c      IF(IS.EQ.3) SIG=SIG*PIDABS
c as compared to free pi+d=pp cr.sec. abs. probab. is higher
c by factor 3, see J.Hufner, Phys.Rep.C21,1(1975))
c and by weight  Z/A*N/A~1/4 of pn pairs: PIDABS=3*1/4=0.75
      IF(IS.EQ.3) SIG=SPIDPP(T)*PIDABS
      RETURN
      END
C
C * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
      FUNCTION SPIDPP(Tpi)
      IMPLICIT REAL *8(A-H,O-Z)
      REAL*8 SPIDPP,Tpi,mpi,md
c  cross section for pi+d=p+p; Tpi is pion lab.kin.ener.
c  Parametrization B.G.Ritche,Phys.Rev. C28(1983)926
      DATA a/-1.2/,b/3.5/,c/7.4/,d/5600./,ER/2136./
      DATA  mpi/139.0/,md/1878./
      TpiM=Tpi*1000.
      E=DSQRT((mpi+md)**2+2.*TpiM*md)
      SPIDPP=a+b/DSQRT(TpiM)+c*10000./((E-ER)**2+d)
      RETURN
      END
C
C * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
      SUBROUTINE  WTAR1(RL,JP,PL,AM,DT,TAU)
      IMPLICIT REAL *8(A-H,O-Z)
      REAL *8 RL,PL,AM,DT,TAU
      COMMON/HCASC/ANUCL1,ANUCL2,ZNUCL1,ZNUCL2,T0,EP1,EP2,
     *VPI,A1,A2,C1,C2,D1,D2,R0N1,R0N2,TF01,TF02,RM1,RM2
      COMMON/RESULT/AN1,AN2,ZN1,ZN2,EN1,EN2,PN1(3),PN2(3),AM1(3),AM2(3)
      COMMON/NUCSP/VPR(3),VTA(3),RADP(3),RADT(3),VEV(3),VRE(3),
     *GEV,GRE,VEP(3),VET(3),GEP,GET
      DIMENSION  RL(3),PL(3),B(3),PS(3),DT(2)
      RO=0.170
      DT(1)=-1.
      IF(AN1.LT.2.)   GO  TO  11
      RN=A1
      IF(ANUCL1.LE.10.)  RN=A1*DSQRT(2.5*ZNUCL1)
      G=1./DSQRT(1.-VPR(1)**2-VPR(2)**2-VPR(3)**2)
      GG=G*G/(G+1.)
      S=(RL(1)-RADP(1))*VPR(1)+(RL(2)-RADP(2))*VPR(2)+(RL(3)-RADP(3))*
     *VPR(3)
      X=RL(1)-RADP(1)+VPR(1)*(S*GG-G*TAU)
      Y=RL(2)-RADP(2)+VPR(2)*(S*GG-G*TAU)
      Z=RL(3)-RADP(3)+VPR(3)*(S*GG-G*TAU)
      R=DSQRT(X**2+Y**2+Z**2)
      IF(R.GE.RN)  GO  TO  11
      DO  10   K=1,3
   10 B(K)=-VPR(K)
      CALL  KINEMA(PL,B,PS,CTS,STS,CFS,SFS,TS,AM)
      P=DSQRT(TS*(TS+2.*AM))
      IF(P.LE.0.001)    GO  TO  11
      IF(JP.EQ.7)  THEN
      RND=RNDMD(-1)
      RZA=ZN1/AN1
      IF(RND.LT.RZA) SIG=QINTG(P,33)
      IF(RND.GE.RZA) SIG=QINTG(P,31)
      ELSE
      IF(JP.LE.6.OR.JP.EQ.8) CALL  SIGSTR(JP,TS,SIG)
      IF(JP.GT.8)  THEN
        CALL  SIGPI(JP,TS,SIGT,0)
        CALL  SIGPI(JP,TS,SIGA,3)
        SIG=SIGT+SIGA
      ENDIF
      ENDIF
      IF(SIG.LT.1.E-10)     GO  TO  11
      BRND=RNDMD(-1)
      CALL DENSND(1,R,RON,ROD,WFL)
      IF(RON.LT.0.0001)      GO  TO  11
      DL=-DLOG(BRND)/(SIG*RON)*ANUCL1/AN1*10.
      R1=DSQRT((X+DL*STS*CFS)**2+(Y+DL*STS*SFS)**2+(Z+DL*CTS)**2)
      IF(R1.GT.RN)    GO TO  11
      DT(1)=DL/P*(TS+AM)*G*(1.-(B(1)*CFS+B(2)*SFS)*STS-B(3)*CTS)
   11 DT(2)=-1.
      IF(AN2.LT.2.)   GO  TO  13
      RN=A2
      IF(ANUCL2.LE.10.)   RN=A2*DSQRT(2.5*ZNUCL2)
      G=1./DSQRT(1.-VTA(1)**2-VTA(2)**2-VTA(3)**2)
      GG=G*G/(G+1.)
      S=(RL(1)-RADT(1))*VTA(1)+(RL(2)-RADT(2))*VTA(2)+(RL(3)-RADT(3))*
     *VTA(3)
      X=RL(1)-RADT(1)+VTA(1)*(S*GG-G*TAU)
      Y=RL(2)-RADT(2)+VTA(2)*(S*GG-G*TAU)
      Z=RL(3)-RADT(3)+VTA(3)*(S*GG-G*TAU)
      R=DSQRT(X**2+Y**2+Z**2)
      IF(R.GE.RN)   GO  TO  13
      DO  12   K=1,3
   12 B(K)=-VTA(K)
      CALL  KINEMA(PL,B,PS,CTS,STS,CFS,SFS,TS,AM)
      P=DSQRT(TS*(TS+2.*AM))
      IF(P.LE.0.001)    GO  TO  13
      IF(JP.EQ.7)  THEN
      RND=RNDMD(-1)
      RZA=ZN1/AN1
      IF(RND.LT.RZA) SIG=QINTG(P,33)
      IF(RND.GE.RZA) SIG=QINTG(P,31)
      ELSE
      IF(JP.LE.6.OR.JP.EQ.8) CALL  SIGSTR(JP,TS,SIG)
      IF(JP.GT.8)  THEN
        CALL  SIGPI(JP,TS,SIGT,0)
        CALL  SIGPI(JP,TS,SIGA,3)
        SIG=SIGT+SIGA
      ENDIF
      ENDIF
      IF(SIG.LT.1.E-10)   GO  TO  13
      CALL DENSND(2,R,RON,ROD,WFL)
      IF(RON.LT.0.0001)      GO  TO  13
      BRND=RNDMD(-1)
      DL=-DLOG(BRND)/(SIG*RON)*ANUCL2/AN2*10.
      R1=DSQRT((X+DL*STS*CFS)**2+(Y+DL*STS*SFS)**2+(Z+DL*CTS)**2)
      IF(R1.GT.RN)   GO  TO  13
      DT(2)=DL/P*(TS+AM)*G*(1.-(B(1)*CFS+B(2)*SFS)*STS-B(3)*CTS)
   13 CONTINUE
      RETURN
      END
C
C * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
      FUNCTION  SSDT(IB1,IB2,IE1,IE2,JP,T,IS)
C  PRODUCTION CROSS SECTION N+N=>DELTA+JP+X
      IMPLICIT REAL *8(A-H,O-Z)
      REAL *8 SSDT,T
      SSDT=0.
      IF((IB1+IB2).NE.2)   RETURN
      IF(IE1+IE2-1)  100,101,102
  100 GO  TO  (11,12,13,14,15,16,17),JP
   11 S1=SIGD(T,5)
      S2=SIGD(T,17)
      S3=SIGD(T,20)
      IF(IS.EQ.1) THEN
      RS=RNDMD(-1)*(S1+S2+S3)
                              SSDT=S1
      IF(RS.GT.S1.AND.RS.LE.(S1+S2))  SSDT=S2
      IF(RS.GT.(S1+S2))               SSDT=S3
      ELSE
C !!!
      SSDT=S1+S2+S3
C !!!
      ENDIF
      RETURN
   12 S1=SIGD(T,6)
      S2=SIGD(T,16)
      S3=SIGD(T,18)
      S4=SIGD(T,19)
      IF(IS.EQ.1) THEN
      RS=RNDMD(-1)*(S1+S2+S3+S4)
                                    SSDT=S1
      IF(RS.GT.S1.AND.RS.LE.(S1+S2))          SSDT=S2
      IF(RS.GT.(S1+S2).AND.RS.LE.(S1+S2+S3))  SSDT=S3
      IF(RS.GT.(S1+S2+S3))                    SSDT=S4
      ELSE
C !!!
      SSDT=S1+S2+S3+S4
C !!!
      ENDIF
      RETURN
   13 S1=SIGD(T,5)
      S2=SIGD(T,6)
      IF(IS.EQ.1) THEN
      RS=RNDMD(-1)*(S1+S2)
                SSDT=S1
      IF(RS.GT.S1)  SSDT=S2
      ELSE
C !!!
      SSDT=S1+S2
C !!!
      ENDIF
      RETURN
   14 SSDT=SIGD(T,19)
      RETURN
   15 S1=SIGD(T,16)
      S2=SIGD(T,17)
      IF(IS.EQ.1) THEN
      RS=RNDMD(-1)*(S1+S2)
                SSDT=S1
      IF(RS.GT.S1)  SSDT=S2
      ELSE
C !!!
      SSDT=S1+S2
C !!!
      ENDIF
      RETURN
   16 S1=SIGD(T,18)
      S2=SIGD(T,20)
      IF(IS.EQ.1) THEN
      RS=RNDMD(-1)*(S1+S2)
                SSDT=S1
      IF(RS.GT.S1)  SSDT=S2
      ELSE
C !!!
      SSDT=S1+S2
C !!!
      ENDIF
   17 CONTINUE
      RETURN
  101 GO  TO  (21,22,23,24,25,26,27),JP
   21 S1=SIGD(T,4)
      S2=SIGD(T,14)
      S3=SIGD(T,15)
      IF(IS.EQ.1) THEN
      RS=RNDMD(-1)*(S1+S2+S3)
                              SSDT=S1
      IF(RS.GT.S1.AND.RS.LE.(S1+S2))  SSDT=S2
      IF(RS.GT.(S1+S2))               SSDT=S3
      ELSE
C !!!
      SSDT=S1+S2+S3
C !!!
      ENDIF
      RETURN
   22 S1=SIGD(T,11)
      S2=SIGD(T,12)
      S3=SIGD(T,13)
      IF(IS.EQ.1) THEN
      RS=RNDMD(-1)*(S1+S2+S3)
                              SSDT=S1
      IF(RS.GT.S1.AND.RS.LE.(S1+S2))  SSDT=S2
      IF(RS.GT.(S1+S2))               SSDT=S3
      ELSE
C !!!
      SSDT=S1+S2+S3
C !!!
      ENDIF
      RETURN
   23 S1=SIGD(T,3)
      S2=SIGD(T,4)
      IF(IS.EQ.1) THEN
      RS=RNDMD(-1)*(S1+S2)
                SSDT=S1
      IF(RS.GT.S1)  SSDT=S2
      ELSE
C !!!
      SSDT=S1+S2
C !!!
      ENDIF
      RETURN
   24 S1=SIGD(T,13)
      S2=SIGD(T,14)
      IF(IS.EQ.1) THEN
      RS=RNDMD(-1)*(S1+S2)
                SSDT=S1
      IF(RS.GT.S1)  SSDT=S2
      ELSE
C !!!
      SSDT=S1+S2
C !!!
      ENDIF
      RETURN
   25 SSDT=SIGD(T,11)
      RETURN
   26 S1=SIGD(T,12)
      S2=SIGD(T,15)
      IF(IS.EQ.1) THEN
      RS=RNDMD(-1)*(S1+S2)
                SSDT=S1
      IF(RS.GT.S1)  SSDT=S2
      ELSE
C !!!
      SSDT=S1+S2
C !!!
      ENDIF
   27 RETURN
  102 GO  TO  (31,32,33,34,35,36,37),JP
   31 S1=SIGD(T,2)
      S2=SIGD(T,8)
      S3=SIGD(T,9)
      S4=SIGD(T,10)
      IF(IS.EQ.1) THEN
      RS=RNDMD(-1)*(S1+S2+S3+S4)
                                    SSDT=S1
      IF(RS.GT.S1.AND.RS.LE.(S1+S2))          SSDT=S2
      IF(RS.GT.(S1+S2).AND.RS.LE.(S1+S2+S3))  SSDT=S3
      IF(RS.GT.(S1+S2+S3))                    SSDT=S4
      ELSE
C !!!
      SSDT=S1+S2+S3+S4
C !!!
      ENDIF
      RETURN
   32 S1=SIGD(T,1)
      S2=SIGD(T,7)
      IF(IS.EQ.1) THEN
      RS=RNDMD(-1)*(S1+S2)
                SSDT=S1
      IF(RS.GT.S1)  SSDT=S2
      ELSE
C !!!
      SSDT=S1+S2
C !!!
      ENDIF
      RETURN
   33 S1=SIGD(T,1)
      S2=SIGD(T,2)
      IF(IS.EQ.1) THEN
      RS=RNDMD(-1)*(S1+S2)
                SSDT=S1
      IF(RS.GT.S1)  SSDT=S2
      ELSE
C !!!
      SSDT=S1+S2
C !!!
      ENDIF
      RETURN
   34 S1=SIGD(T,7)
      S2=SIGD(T,9)
      IF(IS.EQ.1) THEN
      RS=RNDMD(-1)*(S1+S2)
                SSDT=S1
      IF(RS.GT.S1)  SSDT=S2
      ELSE
C !!!
      SSDT=S1+S2
C !!!
      ENDIF
      RETURN
   35 SSDT=SIGD(T,10)
      RETURN
   36 SSDT=SIGD(T,8)
   37 RETURN
      END
C
C * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
      FUNCTION  SSDL(IB1,IB2,IE1,IE2,JP,T)
C  PRODUCTION CROSS SECTION N+N=>L+DELTA+K
      IMPLICIT REAL *8(A-H,O-Z)
      REAL *8 SSDL,T
      SSDL=0.
      IF((IB1+IB2).NE.2)  RETURN
      IF(IE1+IE2-1)  100,101,102
  100 GO  TO  (11,12,13,14,14,14,14),JP
   11 SSDL=SIGD(T,5)
      RETURN
   12 SSDL=SIGD(T,6)
      RETURN
   13 SSDL=SSDT(IB1,IB2,IE1,IE2,JP,T,1)
   14 RETURN
  101 GO  TO  (21,22,13,14,14,14,14),JP
   21 SSDL=SIGD(T,4)
      RETURN
   22 SSDL=SIGD(T,3)
      RETURN
  102 GO  TO  (31,32,13,14,14,14,14),JP
   31 SSDL=SIGD(T,2)
      RETURN
   32 SSDL=SIGD(T,1)
      RETURN
      END
C
C * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
      FUNCTION  SIGD(T,N)
      IMPLICIT REAL *8(A-H,O-Z)
      REAL *8 SIGD,T
C
C   PP  --->  D+ L K+
      SD2(X)=FAR(0.132D0,0.709D0,3.432D0,-0.564D0,X)
C   PP  --->  D+ S+ K0
      SD7(X)=FAR(0.062D0,0.709D0,3.432D0,-0.564D0,X)
C   PP  --->  D+ S0 K+
      SD8(X)=FAR(0.058D0,0.709D0,3.432D0,-0.564D0,X)
C   PP  --->  D0 S+ K+
      SD9(X)=0.25*FAR(0.144D0,0.709D0,3.432D0,-0.564D0,X)
C
      TDL=T-2.432
      TDS=T-2.870
      GO  TO  (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20),N
    1 SIGD=3.*SD2(TDL)
      RETURN
    2 SIGD=SD2(TDL)
      RETURN
    3 SIGD=2.*SD2(TDL)
      RETURN
    4 SIGD=2.*SD2(TDL)
      RETURN
    5 SIGD=3.*SD2(TDL)
      RETURN
    6 SIGD=SD2(TDL)
      RETURN
    7 SIGD=SD7(TDS)
      RETURN
    8 SIGD=SD8(TDS)
      RETURN
    9 SIGD=SD9(TDS)
      RETURN
   10 SIGD=1./2.*(9.*SD7(TDS)-6.*SD8(TDS)+2.*SD9(TDS))
      RETURN
   11 SIGD=3.*SD9(TDS)
      RETURN
   12 SIGD=1./6.*(9.*SD7(TDS)-6.*SD8(TDS)+8.*SD9(TDS))
      RETURN
   13 SIGD=1./6.*(15.*SD7(TDS)-6.*SD8(TDS)+2.*SD9(TDS))
      RETURN
   14 GO  TO  11
   15 GO  TO  12
   16 GO  TO  9
   17 SIGD=1./3.*(9.*SD7(TDS)-12.*SD8(TDS)+8.*SD9(TDS))
      RETURN
   18 GO  TO  8
   19 GO  TO  10
   20 SIGD=1./2.*(9.*SD7(TDS)-12.*SD8(TDS)+8.*SD9(TDS))
      RETURN
      END
C
C * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
      FUNCTION  STHE(U,JP,ST4)
      IMPLICIT REAL *8(A-H,O-Z)
      REAL *8 STHE,U,ST4
      DIMENSION  S19(7)
C   S19    =  TOTAL CROSS SECTION AT T0=19.0 GeV  JP+N
C   FOR  JP=     K+  K0  L  S+  S-  S0  K-
      DATA  S19/1.8,1.8,0.9,0.9,0.9,0.9,0.6/
      U4=3.325
      U19=6.13
      STHE=ST4+(S19(JP)-ST4)*DLOG(U/U4)/DLOG(U19/U4)
      RETURN
      END
C
C * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
      SUBROUTINE  PREPD(U,V,PIN)
      IMPLICIT REAL *8(A-H,O-Z)
      REAL *8 U,V,PIN
      COMMON/NRAPI/NRAPI/DEUT/VDN(3),UDN,T1DN,WDN,VND(3),UND,T1ND,WND
      DIMENSION  V(3),PIN(9)
      WDN=0.
      WND=0.
C     IF(IB1.EQ.0)  RETURN
C
      E1=PIN(8)+PIN(9)
      G=1./DSQRT(1.-V(1)**2-V(2)**2-V(3)**2)
      IF(NRAPI.NE.1)  GO  TO  100
      DO  10  K=1,3
   10 VDN(K)=(PIN(K+3)+V(K)*G*U)/(E1+G*U)
      UDN=(E1+G*U)*DSQRT(1.-VDN(1)**2-VDN(2)**2-VDN(3)**2)
      T1DN=(UDN**2-9.*0.94**2)/(4.*0.94)
      WDN=WFLUC(1,PIN(1),PIN(2),PIN(3))
  100 CONTINUE
      DO  11  K=1,3
   11 VND(K)=(2.*V(K)*G*U-PIN(K+3))/(2.*G*U-E1)
      UND=(2.*G*U-E1)*DSQRT(1.-VND(1)**2-VND(2)**2-VND(3)**2)
      T1ND=(UND**2-9.*0.94**2)/(4.*0.94)
      WND=WFLUC(2,PIN(1),PIN(2),PIN(3))
      RETURN
      END
C
C * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
      SUBROUTINE  POSND(IB1,IE1,IE2,JP,B,STOT,RL,N)
      IMPLICIT REAL *8(A-H,O-Z)
      REAL *8 B,STOT,RL,MPI,MN
      COMMON/DEUT/VDN(3),UDN,T1DN,WDN,VND(3),UND,T1ND,WND
      COMMON/NRAPI/NRAPI
      COMMON /IDN12/ ID1,ID2
      COMMON /AMPROJ/ AMPROJ
      COMMON /KSIGS/ KSIGS
      DIMENSION  B(3),RL(3),AM(11),VD(3)
      DATA  AM/2*0.494,1.115,3*1.192,0.494,0.939,3*0.139/
      DATA MN/0.939/,MPI/0.139/
      IF(NRAPI.GT.3)  RETURN
      AMJ=AM(JP)
      GO  TO  (10,10,11,11,11,11,12,121,122,122,122),JP
C !!! de clarificat EMN==>EMN+0.940 ??
   10 EMN=AMJ+1.115+IB1*.940
      GO  TO  13
   11 EMN=AMJ+0.494+IB1*.940
      GO  TO  13
   12 EMN=AMJ+0.494+IB1*.940+.940
      GO  TO  13
  121 EMN=AMJ+0.939+IB1*.940+.940
      GO  TO  13
  122 EMN=AMJ+0.939+(IB1+1)*.940
      GO  TO  13
   13 DO  22  KD=1,2
      IF(NRAPI.GT.1.AND.KD.EQ.1)  GO  TO  22
      IF(KD-1) 14,14,16
   14 UD=UDN
      IB=1
      WD=WDN
      TD=T1DN
      IE=IE2
      ID=ID2
      DO  15  K=1,3
   15 VD(K)=VDN(K)
      GO  TO  18
C  16 IF(KD.EQ.2.AND.NRAPI.EQ.2)  GO  TO  14
   16 UD=UND
      IB=IB1
      WD=WND
      TD=T1ND
      IE=IE1
      ID=ID1
      DO  17  K=1,3
   17 VD(K)=VND(K)
   18 IF(UD.LE.(EMN+.940))   GO  TO  22
      TNMIN=(EMN**2-(AMPROJ+0.94)**2)/1.88
      TDMIN=((EMN+.940)**2-(AMPROJ+1.88)**2)/3.76
      TD=TD+TNMIN-TDMIN
C     WRITE(16,*) 'UD,TD,ST,TNMIN,TDMIN=', UD,TD,ST,TNMIN,TDMIN
      IF(JP.GE.9)  GO  TO  192
      IF(JP.EQ.8)  GO  TO  93
      IF(JP.EQ.7)  GO  TO  19
C
      ST=0.
      IF(KSIGS.GE.2.AND.IB.EQ.1) THEN
      ST=SZWER(JP,UD)
      GO  TO  201
      ENDIF
      IF((KSIGS.EQ.3.OR.KSIGS.EQ.4).AND.IB.EQ.0) THEN
      ST=(SCT( 1,1,JP,TD,1)+SCT( 1,0,JP,TD,1)+
     +    SCT(-1,1,JP,TD,1)+SCT(-1,0,JP,TD,1)+
     +    SCT( 0,1,JP,TD,1)+SCT( 0,0,JP,TD,1))/6.
      GO  TO  201
      ENDIF
      IF(KSIGS.EQ.5.AND.IB.EQ.0)  THEN
      ST=0.
      SS0=0.494+1.115
      IF(UD.GT.SS0.AND.UD.LT.1.7)  ST=10.*(UD-SS0)
      IF(UD.GE.1.7)                ST=0.09/(UD-1.6)
      GO  TO  201
      ENDIF
C
      IF(IB.NE.0)  GO  TO  91
      IF(ID.EQ.120.OR.ID.EQ.-120.OR.ID.EQ.110) GO  TO  92
C
C MESON+NUCLEON PRODUCTION CROSS SECTION
      ST=(SST(0,1, 1,1,JP,TD,1)+SST(0,1, 1,0,JP,TD,1)+
     +    SST(0,1,-1,1,JP,TD,1)+SST(0,1,-1,0,JP,TD,1)+
     +    SST(0,1, 0,1,JP,TD,1)+SST(0,1, 0,0,JP,TD,1))/6.
      GO  TO  201
   91 CONTINUE
      IF(IB.NE.1)   RETURN
      IF(ID.EQ.1120.OR.ID.EQ.1220) GO  TO 92
C
C BARYON+NUCLEON PRODUCTION CROSS SECTION
      ST=(SST(1,1,1,1,JP,TD,1)+SST(1,1,1,0,JP,TD,1)+
     +    SST(1,1,0,1,JP,TD,1)+SST(1,1,0,0,JP,TD,1))/4.
      GO  TO  201
   92 CONTINUE
C BARYON+NUCLEON AND PION+NUCLEON PRODUCTION CROSS SECTION
      ST=0.5*(SST(IB,1,IE,1,JP,TD,1)+SSDT(IB,1,IE,1,JP,TD,1)+
     +        SST(IB,1,IE,0,JP,TD,1)+SSDT(IB,1,IE,0,JP,TD,1))
      S4=0.5*(SST(IB,1,IE,1,JP,4.D0,1)+SSDT(IB,1,IE,1,JP,4.D0,1)+
     +        SST(IB,1,IE,0,JP,4.D0,1)+SSDT(IB,1,IE,0,JP,4.D0,1))
      IF(TD.GE.4.) ST=STHE(UD,JP,S4)
      GO  TO  201
C  ******  CHANNELS FOR  ANTI-P  PRODUCTION
   93 CONTINUE
      EMN1=EMN-AMJ
      IF(UD.LE.EMN)  RETURN
      PMAX=DSQRT(((UD+AMJ)**2-EMN1**2)*((UD-AMJ)**2-EMN1**2))/(2.*UD)
C !!!!
C  A. SHOR NP A514(1990) p.717;  DANILIEWICZ PR C42(1990) p.1564
C  SEE: PL B256(1991) p.331
      ST=0.01*(UD-EMN)**1.846
C !!!!
      GO  TO  201
   19 CONTINUE
C * * * * * * * * * * * * * * * * * * * * * * * * * * *
C  ******  CHANNELS FOR  K-   PRODUCTION  G: 20.05.93
      IF(UD.LE.EMN)  RETURN
      EMN1=EMN-AMJ
      PMAX=DSQRT(((UD+AMJ)**2-EMN1**2)*((UD-AMJ)**2-EMN1**2))/(2.*UD)
      IF(IB.NE.0)    GO  TO  191
C  M+FL => K-
      ST=PMAX/13.
      GO  TO  201
  191 CONTINUE
C  B+FL => K-
      ST=SIKMI(PMAX)
      IF(TD.LT.4.)    GO  TO  201
      U4=1.88*(4.+1.88)
      PMA4=DSQRT(((U4+AMJ)**2-EMN1**2)*((U4-AMJ)**2-EMN1**2))/(2.*U4)
      S4=SIKMI(PMA4)
      ST=STHE(UD,JP,S4)
C
      GO  TO  201
C * * * * * * * * * * * * * * * * * * * * * * * * * * *
  192 CONTINUE
C  ******  CHANNELS FOR  PI-0+   PRODUCTION
      IF(IB.EQ.1)  GO  TO  193
      S1= CROSEG(0,0,1,1,7,TD,MPI,0)
      S2= CROSEG(0,0,1,2,7,TD,MPI,0)
      S3= CROSEG(0,0,1,3,7,TD,MPI,0)
      ST=(S1+S2+S3)/3.
      GO  TO  201
  193 CONTINUE
      S1= CROSEG(0,0,2,1,7,TD,MN,0)
      S2= CROSEG(0,0,2,2,7,TD,MN,0)
      ST=(S1+S2)/2.
  201 CONTINUE
      SK=ST*WD
      IF(SK.LE.0.)  GO  TO  22
C
      EM1=EMN1+RNDMD(-1)*(UD-EMN)
      E1=(UD**2+AMJ**2-EM1**2)/(2.*UD)
      P1=DSQRT(E1**2-AMJ**2)
      CT1=1.-2.*RNDMD(-1)
      IF(IB.EQ.0)  ID=6
      IF(IB.NE.0)  ID=6+KD
      CALL  POSTR(P1,CT1,AMJ,VD,B,JP,SK,STOT,RL,N,ID)
   21 CONTINUE
   22 CONTINUE
      RETURN
      END
C
C * * * * * * * * * *  G: 20.05.93  * * * * * * * * * *
C
      FUNCTION  SIKMI(PMAX)
      IMPLICIT REAL *8(A-H,O-Z)
      REAL *8 SIKMI,PMAX
C   PRODUCTION CROSS SECTION IN  NN ---> N N K+ K-
C     SIKMI=PMAX/40.   (mb)
      SIKMI=PMAX/40.
      RETURN
      END
C
C * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
      SUBROUTINE  DENSND(NU,R,RON,ROD,W)
      IMPLICIT REAL *8(A-H,O-Z)
      REAL *8 R,RON,ROD,W
      COMMON/HCASC/ANUCL1,ANUCL2,ZNUCL1,ZNUCL2,T0,EPS1,EPS2,VPI,A1,A2,
     *C1,C2,D1,D2,R0N1,R0N2,TF01,TF02,RM1,RM2
      COMMON/RESULT/AN1,AN2,ZN1,ZN2,EN1,EN2,PN1(3),PN2(3),AM1(3),AM2(3)
      COMMON /INSTR/ INSTR
      DIMENSION RFL(6)
c  These parameters are not yet fixed ---------------------
      DATA RFL /2*0.12,2*0.30,2*0.30/
      V2=4.*3.141592/3.*(RFL(INSTR))**3
      GO  TO  (11,12),NU
   11 AN0=ANUCL1
      AN=AN1
      TF0=TF01
      A=A1
      C=C1
      GO  TO  13
   12 AN0=ANUCL2
      AN=AN2
      TF0=TF02
      A=A2
      C=C2
   13 IF(AN0.LE.10.)  GO  TO  14
      TFR=TF0*(((1.+DEXP(-A/C))/(1.+DEXP((R-A)/C)))**0.6666667)
      GO  TO  15
   14 TFR=TF0*DEXP(-(2./3.)*(R/A)**2)
   15 RO=2.*(TFR/0.1985)**1.5
      RON=(AN/AN0)*RO
      W=0.
      ROD=0.
      IF(AN.LT.2.1)  RETURN
      V2R=V2*RON
      IF(V2R.GE.AN)  RETURN
      W=(AN-1.)/2./AN*V2R*
     *(1.-V2R/AN)**(AN-2.)
      ROD=W*RON
C
C
      RETURN
      END
C
C * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
      FUNCTION  WFLUC(K,X,Y,Z)
      IMPLICIT REAL *8(A-H,O-Z)
      REAL *8 WFLUC,X,Y,Z
      COMMON/NRAPI/NRAPI
      COMMON/NUCSP/VPR(3),VTA(3),RADP(3),RADT(3),VEV(3),VRE(3),GEV,GRE,
     *VEP(3),VET(3),GEP,GET
      WFLUC=0.
      IF(NRAPI.EQ.5)  RETURN
      IF(NRAPI.NE.1)  GO  TO  10
      IF(K.EQ.2)      GO  TO   9
      DX=X-RADP(1)
      DY=Y-RADP(2)
      DZ=Z-RADP(3)
      VR=DX*VPR(1)+DY*VPR(2)+DZ*VPR(3)
      GPR=1./DSQRT(1.-VPR(1)**2-VPR(2)**2-VPR(3)**2)
      GG=GPR*GPR/(GPR+1.)
      XN=DX+VPR(1)*VR*GG
      YN=DY+VPR(2)*VR*GG
      ZN=DZ+VPR(3)*VR*GG
      NU=1
      GO  TO  11
    9 DX=X-RADT(1)
      DY=Y-RADT(2)
      DZ=Z-RADT(3)
      VR=DX*VTA(1)+DY*VTA(2)+DZ*VTA(3)
      GTA=1./DSQRT(1.-VTA(1)**2-VTA(2)**2-VTA(3)**2)
      GG=GTA*GTA/(GTA+1.)
      XN=DX+VTA(1)*VR*GG
      YN=DY+VTA(2)*VR*GG
      ZN=DZ+VTA(3)*VR*GG
      NU=2
      GO  TO  11
   10 XN=X
      YN=Y
      ZN=Z
      IF(NRAPI.EQ.2)   NU=1
      IF(NRAPI.EQ.3)   NU=2
   11 R=DSQRT(XN**2+YN**2+ZN**2)
C     IF(R.GT.2.) WRITE(16,*)  'WFLUC=0 FOR NU=', NU
      CALL  DENSND(NU,R,RON,ROD,W)
      WFLUC=W
      RETURN
      END
C
C * * * * * *  * * * * * * * * * * * * * * * * * *  * *
C
      SUBROUTINE  STRANG(BU,BV,NU,N1,N2,BPIN,IIN,BPN,IPN,NREP,JPSTR,
     *SIGT,SABS)
      IMPLICIT REAL *8(A-H,O-Z)
      REAL *8 BU,BV,BPIN,BPN,SIGT,SABS
      COMMON/NRAPI/NRAPI
      COMMON/NUCSP/VPR(3),VTA(3),RADP(3),RADT(3),VEV(3),VRE(3),GEV,GRE,
     *VEP(3),VET(3),GEP,GET
      COMMON/INSTR/ INSTR
      COMMON /AMPROJ/ AMPROJ
      COMMON/DEUT/VDN(3),UDN,T1DN,WDN,VND(3),UND,T1ND,WND
      COMMON/FLUC/PDN(9),PND(9)
      DIMENSION  BV(3),V(3),BPIN(9),IIN(5),BPN(9),PN(9),IPN(5),VN(3),
     *RL(3),RN(3),PIN(9),
     *PN1(9),PN2(9),IPN1(5),IPN2(5),PS1(3),PS2(3),PL1(3),PL2(3)
      SIG=SIGT
      SAB=SABS
      WDN=0.
      WND=0.
      U=BU
      DO  5  K=1,9
      PDN(K)=0.
      PND(K)=0.
      IF(K.LE.3)  V(K)=BV(K)
      PIN(K)=BPIN(K)
    5 PN(K)=BPN(K)
      AMPROJ=PIN(9)
      SVR=0.
      DO  9  K=1,3
      IF(NRAPI.EQ.2.OR.NRAPI.EQ.3)   GO  TO  6
      VN(K)=0.
      RN(K)=0.
      GO  TO  8
    6 IF(NU.EQ.2)  GO  TO  7
      VN(K)=VPR(K)
      RN(K)=RADP(K)
      GO  TO  8
    7 VN(K)=VTA(K)
      RN(K)=RADT(K)
    8 SVR=SVR+VN(K)*(PIN(K)+PN(K))/2.
    9 CONTINUE
      G=1./DSQRT(1.-VN(1)**2-VN(2)**2-VN(3)**2)
      GG=SVR*G/(G+1.)
      DO  10  K=1,3
      RL(K)=(PIN(K)+PN(K))/2.-GG*VN(K)+RN(K)
   10 CONTINUE
C
      DO  13 KR=1,NREP
      IF(NRAPI.EQ.5)       GO  TO  11
c  redefinition of Fermi momenta
      IF(NRAPI.EQ.1)  THEN
       CALL  PARTN(1,N1,PN1,IPN1)
       CALL  PARTN(2,N2,PN2,IPN2)
       C1=DSQRT((PN1(8)*(PN1(8)+2.*PN1(9)))/
     /      (PN1(4)**2+PN1(5)**2+PN1(6)**2))
       C2=DSQRT((PN2(8)*(PN2(8)+2.*PN2(9)))/
     /      (PN2(4)**2+PN2(5)**2+PN2(6)**2))
       do  k=1,3
       PS1(k)=PN1(3+k)*C1
       PS2(k)=PN2(3+k)*C2
       enddo
       CALL  KINEMA(PS1,VPR,PL1,CT1,ST1,CF1,SF1,TL1,PN1(9))
       CALL  KINEMA(PS2,VTA,PL2,CT2,ST2,CF2,SF2,TL2,PN2(9))
       do  k=1,3
       PIN(3+k)=PL1(k)
       PN(3+k) =PL2(k)
       enddo
       PIN(8)=DSQRT(PL1(1)**2+PL1(2)**2+PL1(3)**2+0.940**2)-0.940 
       PN(8) =DSQRT(PL2(1)**2+PL2(2)**2+PL2(3)**2+0.940**2)-0.940 
       PIN(9)=0.940
       PN(9) =0.940
      ELSE
       CALL  PARTN(NU,N2,PN2,IPN2)
       do  k=1,3
       PN(3+k)=PN2(3+k)
       enddo
       PN(8)=PN2(8)
       PN(9)=0.94
      ENDIF
      PX1=PIN(4)
      PY1=PIN(5)
      PZ1=PIN(6)
      AM1=PIN(9)
      PX2=PN(4)
      PY2=PN(5)
      PZ2=PN(6)
      AM2=PN(9)
      CALL CROSEC(1,ID1,ID2,PX1,PY1,PZ1,AM1,PX2,PY2,PZ2,AM2,SITO,0)
      SIG=SITO
   11 CONTINUE
      CALL  TINVU(PIN,PN,U,V,T12)
      IF(INSTR.EQ.3.OR.NRAPI.EQ.5)  GO  TO  12
      CALL  NFERMI(PIN,PN,NU,N2)
      IF(NRAPI.EQ.1.OR.NRAPI.EQ.3) CALL TINVU(PIN,PND,UND,VND,T1ND)
      IF(NRAPI.EQ.1)               CALL TINVU(PDN, PN,UDN,VDN,T1DN)
      IF(NRAPI.EQ.2)               CALL TINVU(PIN,PND,UND,VND,T1ND)
   12 CONTINUE
      STOT=SIG+SAB
      CALL  CONECT(U,V,VN,PIN,IIN,PN,IPN,RL,NREP,JPSTR,STOT)
   13 CONTINUE
      RETURN
      END
C
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
      SUBROUTINE NFERMI(PIN,PN,NU,N2)
      IMPLICIT REAL *8(A-H,O-Z)
      REAL *8 PIN,PN
c      REAL *4 XM,XF,FGAM
      COMMON/DPFER/DPFER(41)
      COMMON/NRAPI/NRAPI
      COMMON/DEUT/VDN(3),UDN,T1DN,WDN,VND(3),UND,T1ND,WND
      COMMON/FLUC/PDN(9),PND(9)
      COMMON/CENTER/ XC(2,240),YC(2,240),ZC(2,240),IZ(2,240)
      COMMON/NUCSP/VPR(3),VTA(3),RADP(3),RADT(3),VEV(3),VRE(3),GEV,GRE,
     *VEP(3),VET(3),GEP,GET
      COMMON/WSUM/ WSUM
      DIMENSION  PIN(9),PN(9),PS(3),PL(3)
C  ***** ONLY LAB SYSTEM  ***** ????
      X=PIN(1)
      Y=PIN(2)
      Z=PIN(3)
      P=DSQRT(PN(4)**2+PN(5)**2+PN(6)**2)
      TF=TFERMI(X,Y,Z,NU)
      PF=DSQRT(TF*(TF+1.88))
      IF(NU.EQ.1)  WDN=WFLUC( 1,X,Y,Z)
      IF(NU.EQ.2)  WND=WFLUC( 2,X,Y,Z)
C  @@@@@@@@@@@@
C     IT IS TO IMPROVE SAMPLING OF FLUCTON MOMENTUM.
C     ONE SHOULD CHOOSE ALSO   P<PF  IN FUNCTION BLTZ .
C     XM = 0.5*(0.55/.110)**2
C     XF = 0.5*(PF/.110)**2
C     FGAM=(GAMMP(1.5,XM)-GAMMP(1.5,XF))/GAMMP(1.5,XM)
C     WDN=WDN*FGAM
C     WND=WND*FGAM
C  @@@@@@@@@@@@
C     DRND=RNDMD(-1)
C     P=PF*DRND**(1./3.)
C     FI=6.283185*RNDMD(-1)
C     CT=1.-2.*RNDMD(-1)
C     ST=DSQRT(1.-CT*CT)
C     PN(1)=PIN(1)
C     PN(2)=PIN(2)
C     PN(3)=PIN(3)
C     PN(4)=P*ST*DCOS(FI)
C     PN(5)=P*ST*DSIN(FI)
C     PN(6)=P*CT
C     PN(7)=0.
C     PN(8)=DSQRT(P**2+0.940**2)-0.940
C     PN(9)=0.940
C
C     IF(NRAPI.EQ.2)  GO  TO  20
      PD=BLTZ(0.550D0,0.115D0,PF)
      FI=6.283185*RNDMD(-1)
      CT=1.-2.*RNDMD(-1)
      ST=DSQRT(1.-CT*CT)
      PND(1)=PIN(1)
      PND(2)=PIN(2)
      PND(3)=PIN(3)
      PND(4)=PD*ST*DCOS(FI)
      PND(5)=PD*ST*DSIN(FI)
      PND(6)=PD*CT
      PND(7)=0.
      PND(8)=DSQRT(PD**2+0.940**2)-0.940
      PND(9)=0.940
               WND=WFLUC( 2,X,Y,Z)
      W=WND
      WSUM=WSUM+W
      IFER=P/0.02+1
      IF(IFER.GT.38) IFER=38
      DPFER(IFER)=DPFER(IFER)+1.-W
      DPFER(39)=DPFER(39)+P*(1.-W)
      DPFER(40)=DPFER(40)+P*P*(1.-W)
      DPFER(41)=DPFER(41)+1.-W
      IFER=PD/0.02+1
C  !!!
c     IF(IFER.GT.0)   GO  TO  10
C  !!!
      IF(IFER.GT.38) IFER=38
      DPFER(IFER)=DPFER(IFER)+W
      DPFER(39)=DPFER(39)+PD*W
      DPFER(40)=DPFER(40)+PD*PD*W
      DPFER(41)=DPFER(41)+W
C
   10 CONTINUE
      IF(NRAPI.NE.1)  RETURN
      X=XC( 1,N2)
      Y=YC( 1,N2)
      Z=ZC( 1,N2)
      R=DSQRT(X**2+Y**2+Z**2)
C     CALL DENSND(1,R,RHON,RHOD,WDN)
      WDN=WFLUC( 1,X,Y,Z)
       PF=.270
      P=BLTZ(.550D0,0.115D0,PF)
C  @@@@@@@@@@@@
C     WDN=WDN*FGAM
C  @@@@@@@@@@@@
      FI=6.283185*RNDMD(-1)
      CT=1.-2.*RNDMD(-1)
      ST=DSQRT(1.-CT*CT)
      PS(1)=P*ST*DCOS(FI)
      PS(2)=P*ST*DSIN(FI)
      PS(3)=P*CT
      CALL  KINEMA(PS,VPR,PL,CTL,STL,CFL,SFL,TL,.94D0)
      PDN(1)=PIN(1)
      PDN(2)=PIN(2)
      PDN(3)=PIN(3)
      PDN(4)=PL(1)
      PDN(5)=PL(2)
      PDN(6)=PL(3)
      PDN(7)=0.
      PDN(8)=TL
      PDN(9)=0.940
      RETURN
C 20  PDN(1)=PIN(1)
C     PDN(2)=PIN(2)
C     PDN(3)=PIN(3)
C     PDN(4)=PD*ST*DCOS(FI)
C     PDN(5)=PD*ST*DSIN(FI)
C     PDN(6)=PD*CT
C     PDN(7)=0.
C     PDN(8)=DSQRT(PD**2+0.940**2)-0.940
C     PDN(9)=0.940
C     RETURN
      END
C
C * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
      SUBROUTINE  STRPR1(NC,INE,SG)
      IMPLICIT REAL *8(A-H,O-Z)
      REAL *8 SG
      COMMON/INDZAP/INDZAP
      COMMON/INDSTR/JS,NREP,JCA
      COMMON/INPT/INPT
      COMMON/NPRI/NPRI
      COMMON/NPIMI/NPIMI
      COMMON/STRSP/STRSP(4,18,41),STRSP1(4,6,41),DP,DT,DPT,DY,
     *PM,PTM,YMI,YMA
      COMMON/NRESC/RESC(2),STRIN
      COMMON/STRY/STRY(4)
      COMMON/ISRES/ISRES
      COMMON/DPFER/DPFER(41)
      COMMON/WSUM/ WSUM
      COMMON /IWSTR/ IWSTR
      DIMENSION  DOM(18),TETF(18)
      CHARACTER*4 TEX1(11),TEX2(11),TEX3(10),TEX4(2),TEX5(5),TEX6(8)
      CHARACTER*4 TEXT,TEX20(6),TEX7(8),TEX8(8)
      DIMENSION   C1(18),C2(6)
      DATA TEX1/'INVA','RIAN','T (E','DS/D','P**3',
     *')   ','CROS','S SE','CTIO','NS  ','OF  '/,
     *TEX2/' K+ ',' K0 ','LAMD','SIG+','SIG-','SIG0',' K- ','AP  ',
     *'Pi- ','Pi0 ','Pi+ '/,
     *TEX3/'DIFF','EREN','TIAL','  DI','STRI','BUTI','ONS(','1/ND',
     *'N/DX',') OF'/,
     *TEX4/'FROM','    '/,
     *TEX5/'PIPI','PI+N','N+N ','NNPT','ALL '/,
     *TEX6/'DS/D','T/DO','    ','CROS','S SE','CTIO','NS  ','OF  '/,
     *TEX7/'DS/D','P DO','    ','CROS','S SE','CTIO','NS  ','OF  '/,
     *TEX8/'E*DS','/D3P','    ','CROS','S SE','CTIO','NS  ','OF  '/,
     *TEX20/' ETA','ETAP','RHO0','OMEG',' PHI','  AP'/
      DPI=3.141592/18.
      RESC(1)=0.
      RESC(2)=0.
      STRIN=0.
      NPIMI=0
      WSUM=0.
      DP=0.100
      IF(JS.LE.2)  DP=DP/2.
      IF(JS.GE.3.AND.JS.LE.6)  DP=DP*2.
      IF(JS.EQ.7)  DP=DP/2.
      IF(JS.EQ.8)  DP=DP*2.
      IF(JS.GE.9)  DP=0.010
      DT=5.
      DPT=0.10
      DY=0.1
      DO   M=1,18
      TETF(M)=10.*(M-1)+5.
      DOM(M)=(DCOS((M-1)*DPI)-DCOS(M*DPI))*6.283185
      ENDDO
      if(INPT.eq.6) then
      DO  M=1,18
        TETF(M)=DPT*(M-1)+0.5*DPT
        DOM(M)=DPT*3.141593
      ENDDO
*        DP=DY
      endif
      PM=3.6
      IF(JS.LE.2)  PM=PM/2.
      IF(JS.GE.3.AND.JS.LE.6)  PM=PM*2.
      IF(JS.EQ.7)  PM=PM/2.
      IF(JS.EQ.8)  PM=PM*2.
      IF(JS.GE.9)  PM=0.350
      PTM=1.8
      YMI=-1.
      YMA=YMI+DY*36.
      DO  10  L=1,4
      STRY(L)=0.
      DO  10  M=1,18
      DO  10  N=1,41
      DPFER(N)=0.
      STRSP(L,M,N)=0.
      IF(M.LE.6)   STRSP1(L,M,N)=0.
   10 CONTINUE
      IF(JCA.GE.0)  GO  TO  11
      RETURN
   11 READ(11)  DCA,DEL,NPIMI,WSUM,RESC,STRIN,STRY,STRSP,STRSP1
C     BACKSPACE  11
      REWIND  11
      NC=DCA
      INE=DEL
      WRITE(16,12) DCA
   12 FORMAT(  1X,F7.0,' EVENTS ARE READ FROM FILE 11'  )
      RETURN
C  * * * * * * * * * * * * * * * * * * * * * * * *
      ENTRY  STRPR2(NC,INE,SG)
      IF(INDZAP.EQ.0)   RETURN
      SIN=SG*FLOAT(NC)/FLOAT(NC+INE)
      WRITE(16,13) NC,INE,SIN
   13 FORMAT(/10X,'NCAS=',I5,2X,'INTEL=',I5,2X,'INELASTIC CROSS',
     *' SECTION=',E13.6,'MB'//)
      IF(IWSTR.EQ.1)  CY=1./FLOAT(NC)
      IF(IWSTR.EQ.2)  CY=1./FLOAT(NC)/SIN
      YPP=STRY(1)*CY
      YPN=STRY(2)*CY
      YNN=STRY(3)*CY
      Y00=STRY(4)*CY
      YTOT=YPP+YPN+YNN
      SPP=YPP*SIN
      SPN=YPN*SIN
      SNN=YNN*SIN
      S00=Y00*SIN
      STOT=YTOT*SIN
      RSC1=0.
      RSC2=0.
      STIN=STRIN/FLOAT(NC*NREP)
      WRITE(16,112) STIN
  112 FORMAT(5X,'MEAN NUMBER OF STRANGE SOURCES = ',E10.3/)
      IF(STIN.GT.1.E-20)   RSC1=RESC(1)/FLOAT(NC*NREP)/STIN
      IF(STIN.GT.1.E-20)   RSC2=RESC(2)/FLOAT(NC*NREP)/STIN
      WRITE(16,113) RSC1,RSC2
  113 FORMAT(/5X,'MEAN NUMBER OF RESCATTERING = ',E10.3,2X,E10.3/)
      TEXT=TEX2(JS)
      IF(ISRES.EQ.1)  TEXT=TEX20(JS)
      WRITE(16,14) TEXT,YPP,SPP,YPN,SPN,YNN,SNN,Y00,S00,YTOT,STOT
c     WRITE( *,14) TEXT,YPP,SPP,YPN,SPN,YNN,SNN,Y00,S00,YTOT,STOT
   14 FORMAT(/1X,'MEAN MULTIPLICITY (PROD.CROSS SECTION) OF ',A4/
     *4X,'FROM PI+PI=',2(E13.6,2X)/
     *5X,'FROM PI+N=',2(E13.6,2X)/6X,   'FROM N+N=',2(E13.6,2X)/
     *2X,'FROM NP0+NT0=',2(E13.6,2X)/9X,'TOTAL=',2(E13.6,2X)/)
      PIMI=FLOAT(NPIMI)/FLOAT(NC)
      WRITE(16,114) PIMI
  114 FORMAT(/5X,'MEAN MULTIPLICITY OF PI-MINUS=',E13.6//)
      IF(INDZAP.EQ.2)  GO  TO  1994
      IF(((NC/NPRI)*NPRI).NE.NC)  GO  TO  124
 1994 CONTINUE
      DO  19  K=1,5
      IF(INPT.LE.2)  WRITE(16,18) TEX1,TEX2(JS),TEX4,TEX5(K)
      IF(INPT.EQ.3)  WRITE(16,28) TEX6,TEX2(JS),TEX4,TEX5(K)
      IF(INPT.EQ.4)  WRITE(16,28) TEX7,TEX2(JS),TEX4,TEX5(K)
      IF(INPT.EQ.5)  WRITE(16,28) TEX8,TEX2(JS),TEX4,TEX5(K)
      DO  119  MP=1,2
      M1=(MP-1)*9+1
      M2=(MP-1)*9+9
      IF(INPT.EQ.1)   WRITE(16,118) (TETF(MM),MM=M1,M2)
      IF(INPT.EQ.2)   WRITE(16,128) (TETF(MM),MM=M1,M2)
      IF(INPT.EQ.3)   WRITE(16,138) (TETF(MM),MM=M1,M2)
      IF(INPT.EQ.4)   WRITE(16,148) (TETF(MM),MM=M1,M2)
      IF(INPT.EQ.5)   WRITE(16,158) (TETF(MM),MM=M1,M2)
      DO 119  L=1,41
      DO  15  M=M1,M2
      IF(K.LE.4)  C1(M)=STRSP(K,M,L)*SIN*CY
      IF(K.EQ.5)
     *  C1(M)=(STRSP(1,M,L)+STRSP(2,M,L)+STRSP(3,M,L))*SIN*CY
      C1(M)=C1(M)/DP/DOM(M)
   15 CONTINUE
      P=DP*(L-1)+DP/2.
      IF(L.LE.36)  WRITE(16,16) P,(C1(MM),MM=M1,M2)
      IF(L.GT.36)  WRITE(16,17)   (C1(MM),MM=M1,M2)
   16 FORMAT(3X,F5.3,2X,9(E10.3,2X))
   17 FORMAT(10X,9(E10.3,2X))
   18 FORMAT(//5X,15A4/)
   28 FORMAT(//5X,12A4/)
  118 FORMAT(/3X,'P/TETA',4X,9(F4.0,8X)/)
  128 FORMAT(/2X,'T*/TETA',4X,9(F4.0,8X)/)
  138 FORMAT(/2X,'TL/TETA',4X,9(F4.0,8X)/)
  148 FORMAT(/2X,'PL/TETA',4X,9(F4.0,8X)/)
  158 FORMAT(/2X,'TL/TETA',4X,9(F4.0,8X)/)
  119 CONTINUE
   19 CONTINUE
      DO  24  K=1,5
      WRITE(16,20) TEX3,TEX2(JS),TEX4,TEX5(K)
      WRITE(16,120)
   20 FORMAT(//5X,14A4/)
  120 FORMAT(/5X,'P',6X,'F(P)',6X,'TETA',4X,'F(TETA)',6X,'PT',6X,'F(PT)'
     *,7X,'Y',7X,'F(Y)',6X,'CT*',5X,'F(CT*)',6X,'T*',6X,'F(T*)'/)
      DO  24  L=1,41
      DO  21  M=1,6
      C2(M)=0.
      IF(K.EQ.1.AND.YPP.GT.1.E-20)  C2(M)=STRSP1(1,M,L)/YPN*CY
      IF(K.EQ.2.AND.YPN.GT.1.E-20)  C2(M)=STRSP1(2,M,L)/YPN*CY
      IF(K.EQ.3.AND.YNN.GT.1.E-20)  C2(M)=STRSP1(3,M,L)/YNN*CY
      IF(K.EQ.4.AND.Y00.GT.1.E-20)  C2(M)=STRSP1(4,M,L)/Y00*CY
      IF(K.EQ.5.AND.YTOT.GT.1.E-20) C2(M)=(STRSP1(1,M,L)+
     +STRSP1(2,M,L)+STRSP1(3,M,L))/YTOT*CY
   21 CONTINUE
      P=DP*(L-1)+DP/2.
      TETA=DT*(L-1)+DT/2.
      PT=DPT*(L-1)+DPT/2.
      Y=YMI+(L-1)*DY+DY/2.
      CTS=(L-1)*0.1+0.05-1.
      TS=(L-1)*0.05+0.025
      IF(L.LE.20) WRITE(16,22) P,C2(1),TETA,C2(2),PT,C2(3),Y,C2(4),
     *CTS,C2(5),TS,C2(6)
      IF(L.GT.20.AND.L.LE.36)
     *WRITE(16,122) P,C2(1),TETA,C2(2),PT,C2(3),
     *Y,C2(4),TS,C2(6)
  122 FORMAT(3X,4(F7.3,1X,E10.3,1X),19X,F7.3,1X,E10.3)
      IF(L.GT.36)  WRITE(16,23) C2
   22 FORMAT(3X,6(F7.3,1X,E10.3,1X))
   23 FORMAT(3X,6( 8X,E10.3,1X))
   24 CONTINUE
C 124 IF(INDZAP.EQ.0)   RETURN
  124 WRITE(16,*) 'INDZAP=', INDZAP
      DCA=NC
      DEL=INE
      WRITE(11) DCA,DEL,NPIMI,WSUM,RESC,STRIN,STRY,STRSP,STRSP1
C     BACKSPACE 11
      REWIND  11
   25 FORMAT(/3X,I5,' EVENTS ARE WRITTEN ON THE FILE 11')
      C=DPFER(41)
      IF(C.LE.0.D0)  RETURN
      WSM=WSUM/C
      WRITE(16,420) WSM
  420 FORMAT(/2X,'EFFECTIVE NUCLEON SPECTRA',5X,'WSM=',E10.3/)
      DLP=0.020
      DO  125  K=1,40
      PF=FLOAT(K-1)*DLP+DLP/2.
      IF(C.NE.0.)  DPFRM=DPFER(K)/C
      IF(K.EQ.39)  AVRP=DPFRM
      IF(K.EQ.40)  DPFRM=DSQRT(DPFRM-AVRP**2)
      IF(K.LE.38)  DPFRM=DPFRM /DLP
      IF(K.LE.38)  WRITE(16,421) PF,DPFRM
      IF(K.EQ.39)  WRITE(16,422) DPFRM
      IF(K.EQ.40)  WRITE(16,423) DPFRM
  125 CONTINUE
  421 FORMAT(1X,F7.3,1X,E10.3)
  422 FORMAT(1X,'<P>    ',1X,E10.3)
  423 FORMAT(1X,'<DISP> ',1X,E10.3)
      RETURN
      END
C
C * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
      FUNCTION BLTZ(PMAX,SIG,PF)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 BLTZ, PMAX, SIG,PF
      XMAX=PMAX**2
      AC=2.*SIG**2
      NT=0
   11 RND1=RNDMD(-1)
      NT=NT+1
      X = -AC*DLOG(1.-RND1*(1.-DEXP(-XMAX/AC)))
      P=DSQRT(X)
      RND2=RNDMD(-1)
      PM = RND2*PMAX
      IF(P.GT.PM) GO  TO 11
C @@@@@@@@@@@@@@@@@@@@
C     IF(P.LT.PF) GO  TO 11
C @@@@@@@@@@@@@@@@@@@@
      BLTZ=P
C     WRITE(16,*) 'NT,P=', NT,P
      END
C
C * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
      SUBROUTINE  CONSTR(P,CT,AM,V,B,JP,ST,STOT,RL ,N,IBA,IDD)
      IMPLICIT REAL *8(A-H,O-Z)
      REAL *8 P,CT,AM,E,V,B,ST,STOT,RL
      DIMENSION  V(3),B(3),PS(3),PL(3),RL(3)
      COMMON/STRSP/STRSP(4,18,41),STRSP1(4,6,41),DP,DT,DPT,DY,
     *PM,PTM,YMI,YMA
      COMMON /CMNN/ BNN,GNN
      COMMON /INPT/ INPT
      COMMON /NRAPI/ NRAPI
      COMMON /STRY/ STRY(4)
      COMMON /BGLAB/ BLAB,GLAB,KOBR
      COMMON /KSYST/ KSYST
      COMMON/DEUT/VDN(3),UDN,T1DN,WDN,VND(3),UND,T1ND,WND
      COMMON/FLUC/PDN(9),PND(9)
      COMMON/IPRIMC/IPRIMC
      COMMON /IWSTR/ IWSTR
      BM2=B(1)**2+B(2)**2+B(3)**2
      SIT=DSQRT(1.-CT*CT)
      FI=6.283185*RNDMD(-1)
      PS(1)=P*SIT*DCOS(FI)
      PS(2)=P*SIT*DSIN(FI)
      PS(3)=P*CT
      CALL  KINEMA(PS,V,PL,CTL,STL,CFL,SFL,TL,AM)
      IF(BM2.LT.1.0E-7)   GO  TO  11
      CALL  KINEMA(PL,B,PS,CTL,STL,CFL,SFL,TL,AM)
      DO  10  K=1,3
      PL(K)=PS(K)
   10 CONTINUE
   11 CONTINUE
      IF(IWSTR.EQ.1) W=ST/STOT/FLOAT(N)
      IF(IWSTR.EQ.2) W=ST/FLOAT(N)
      IF(IDD.EQ.1)  W=W*(1.-WND-WDN)
      IF(IDD.EQ.2)  W=W*WDN
      IF(IDD.EQ.3)  W=W*WND
      CALL  RESCA1(PL,AM,JP,RL,CTL,STL,CFL,SFL,TL,W,WABS,TAUABS)
      GO  TO  (1,2,3),KSYST
    1 VSYS=0.
      GO  TO  4
    2 T0=(GLAB-1.)*0.940
      P0=SQRT(T0*(T0+1.88))
      VSYS=P0/(T0+1.88)
      GO  TO  4
    3 VSYS=BNN
      GO  TO  4
    4 GSYS=1./DSQRT(1.-VSYS**2)
      PL(3)=GSYS*(PL(3)+(TL+AM)*VSYS)
      PML=DSQRT(PL(1)**2+PL(2)**2+PL(3)**2)
      CTL=PL(3)/PML
      STL=DSQRT(1.-CTL**2)
      TL=DSQRT(PML**2+AM**2)-AM
      IF(KOBR.EQ.0)  GO  TO  111
      PZL=GLAB*(-PL(3)+(TL+AM)*BLAB)
      PML=DSQRT(PL(1)**2+PL(2)**2+PZL**2)
      CTL=PZL/PML
      STL=DSQRT(1.-CTL**2)
      TL=DSQRT(PML**2+AM**2)-AM
  111 CONTINUE
      PML=DSQRT(TL*(TL+2.*AM))
      EL=TL+AM
      ES=GNN*(EL-PML*CTL*BNN)
      TS=ES-AM
      PMS=DSQRT(TS*(TS+2.*AM))
      PZS=GNN*(PML*CTL-EL*BNN)
      CTS=PZS/PMS+1.
      PT=PML*STL
      W2=W/(4.*0.05*3.141592*PMS)
C
C     WT=W
      WT=0.
      IF(PT.GT.0.)  WT=W/PT
C
      IF(INPT.EQ.1)  W1=W*EL/PML**2
      IF(INPT.EQ.2)  W1=W*GNN*DABS(PML-EL*BNN*CTL)/PML**2
      IF(INPT.EQ.3)  W1=W
      IF(INPT.EQ.4)  W1=W
      IF(INPT.EQ.5)  W1=W/PML
c      IF(INPT.EQ.6)  W1=WT    ! invar. cross section
      IF(INPT.EQ.6)  W1=W     ! abundancy
      Y=0.5*DLOG((TL+AM+PML*CTL)/(TL+AM-PML*CTL))
      IF(DABS(CTL).GT.1.)   CTL=DSIGN(1.D0,CTL)
      TETA=DACOS(CTL)*180./3.141592
      if(INPT.eq.6) then
      JT=PT/DPT+1
      else
      JT=TETA/10.+1
      endif
      IF(JT.GT.18)  JT=18
      IBS=IBA+1
      STRY(IBS)=STRY(IBS)+W
      IF(INPT.EQ.1) CALL DHIST4(PML,0.D0,PM,DP,STRSP,4,18,41,W1,IBS,JT)
      IF(INPT.EQ.2) CALL DHIST4( TS,0.D0,PM,DP,STRSP,4,18,41,W1,IBS,JT)
      IF(INPT.EQ.3) CALL DHIST4( TL,0.D0,PM,DP,STRSP,4,18,41,W1,IBS,JT)
      IF(INPT.EQ.4) CALL DHIST4(PML,0.D0,PM,DP,STRSP,4,18,41,W1,IBS,JT)
      IF(INPT.EQ.5) CALL DHIST4( TL,0.D0,PM,DP,STRSP,4,18,41,W1,IBS,JT)
      IF(INPT.EQ.6) CALL DHIST4(  Y,YMI,YMA,DY,STRSP,4,18,41,W1,IBS,JT)
      CALL  DHIST4(PML,0.D0,PM,DP,STRSP1,4,6,41,W,IBS,1)
      CALL  DHIST4(TETA,0.D0,180.D0,DT,STRSP1,4,6,41,W,IBS,2)
      CALL  DHIST4(PT,0.D0,PTM,DPT,STRSP1,4,6,41,W,IBS,3)
      CALL  DHIST4( Y,YMI,YMA,DY,STRSP1,4,6,41,W,IBS,4)
      CALL  DHIST4(CTS,0.D0,2.D0,0.1D0,STRSP1,4,6,41,W,IBS,5)
      CALL  DHIST4(TS,0.D0,1.8D0,0.05D0,STRSP1,4,6,41,W2,IBS,6)
      IF(IBA.EQ.2.AND.NRAPI.EQ.1)  GO TO 12
      IF(IBA.EQ.2.AND.IPRIMC.EQ.1)  GO TO 12
      RETURN
   12 IBS=4
      STRY(IBS)=STRY(IBS)+W
      IF(INPT.EQ.1) CALL DHIST4(PML,0.D0,PM,DP,STRSP,4,18,41,W1,IBS,JT)
      IF(INPT.EQ.2) CALL DHIST4( TS,0.D0,PM,DP,STRSP,4,18,41,W1,IBS,JT)
      IF(INPT.EQ.3) CALL DHIST4( TL,0.D0,PM,DP,STRSP,4,18,41,W1,IBS,JT)
      IF(INPT.EQ.4) CALL DHIST4(PML,0.D0,PM,DP,STRSP,4,18,41,W1,IBS,JT)
      IF(INPT.EQ.5) CALL DHIST4( TL,0.D0,PM,DP,STRSP,4,18,41,W1,IBS,JT)
      IF(INPT.EQ.6) CALL DHIST4(  Y,YMI,YMA,DY,STRSP,4,18,41,W1,IBS,JT)
      CALL  DHIST4(PML,0.D0,PM,DP,STRSP1,4,6,41,W,IBS,1)
      CALL  DHIST4(TETA,0.D0,180.D0,DT,STRSP1,4,6,41,W,IBS,2)
      CALL  DHIST4(PT,0.D0,PTM,DPT,STRSP1,4,6,41,WT,IBS,3)
      CALL  DHIST4( Y,YMI,YMA,DY,STRSP1,4,6,41,W,IBS,4)
      CALL  DHIST4(CTS,0.D0,2.D0,0.1D0,STRSP1,4,6,41,W,IBS,5)
      CALL  DHIST4(TS,0.D0,1.8D0,0.05D0,STRSP1,4,6,41,W2,IBS,6)
      RETURN
      END
C
C * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
      SUBROUTINE  CONECT(U,V,B,PIN,IIN,PN,IPN,RL,N,JP,STOT)
      IMPLICIT REAL *8(A-H,O-Z)
      REAL *8 U,V,B,PIN,PN,RL,STOT
      CHARACTER*8 PNA1,PNA2
      COMMON/KSTST/KSTST
      COMMON/NRAPI/NRAPI
      COMMON/INSTR/ INSTR
      COMMON /IDN12/ ID1,ID2
      COMMON/ISRES/ISRES
      COMMON/DEUT/VDN(3),UDN,T1DN,WDN,VND(3),UND,T1ND,WND
      COMMON/FLUC/PDN(9),PND(9)
      DIMENSION  PIN(9),IIN(5),PN(9),IPN(5),AM(11),AMR(6)
      DIMENSION  V(3),B(3),RL(3),PS(3),IDSTR(11),IDSRR(6)
      DATA  AM/2*0.494,1.115,3*1.192,0.494,0.940,3*0.139/
      DATA  AMR/.549,.958,.770,.783,1.020,0.940/
      DATA IDSTR
     */ 130,    230,    2130,   1130,   2230,   1230,  -130, -1120,
     * -120,    110,     120/
      DATA IDSRR
     */ 220,    330,     111,    221,    331,  -1120/
      IBA=IIN(4)+IPN(4)
      CALL  PANUID(ID1,IK1,PNA1)
      CALL  PANUID(ID2,IK2,PNA2)
      IF(ISRES.EQ.0) IDS=IDSTR(JP)
      IF(ISRES.EQ.1) IDS=IDSRR(JP)
      IF(IDS.EQ.ID1.OR.IDS.EQ.ID2)   RETURN
C
C   LIMITATION DUE TO THRESHOLD   !!!!!!
      IF(ISRES.EQ.0)  GO  TO  9
      IF(IBA.EQ.0)              UMIN= 0.140+AMR(JP)
      IF(IBA.NE.0)              UMIN= FLOAT(IBA)*0.940+AMR(JP)
      IF(JP.EQ.6)               UMIN= FLOAT(IBA)*0.940+2.*AMR(JP)
      GO  TO  10
   9  IF(IBA.EQ.0)              UMIN= 2.* AM(JP)
      IF(IBA.NE.0.AND.JP.NE.7)  UMIN=1.115+0.494+FLOAT(IBA-1)*.94
      IF(IBA.NE.0.AND.JP.EQ.7)  UMIN=FLOAT(IBA)*0.940+2.*0.494
      IF(IBA.NE.0.AND.JP.EQ.8)  UMIN=FLOAT(IBA)*0.940+2.*0.940
      IF(IBA.EQ.1.AND.JP.GE.9)  UMIN=0.140+0.140+0.940
      IF(IBA.EQ.2.AND.JP.GE.9)  UMIN=0.140+0.940+0.940
C
   10 IF(U.LE.(UMIN+0.001))   GO TO 12
C      WRITE( *,*) 'IBA,U,UMIN=', IBA,U,UMIN
C
C
      CALL  ACTSTR(PIN,PN,SK,PS,JP,NSTR)
      IF(NSTR.EQ.0)  GO  TO  11
      P1=DSQRT(PS(1)**2+PS(2)**2+PS(3)**2)
      CT1=PS(3)/P1
      AM1=AM(JP)
      E1=DSQRT(P1**2+AM1**2)
      CALL  CONSTR(P1,CT1,AM1,V,B,JP,SK,STOT,RL ,N,IBA,1)
   11 CONTINUE
C
   12 IF(INSTR.NE.4.OR.NRAPI.EQ.5)        RETURN
      IF(UND.LE.(UMIN+0.001))           GO  TO  14
C     IF(NRAPI.EQ.2.OR.NRAPI.GT.3)       GO  TO  14
C     WRITE(16,*) 'IBA,UND,UMIN=', IBA,UND,UMIN
C
      CALL  ACTSTR(PIN,PND,SK,PS,JP,NSTR)
      IF(NSTR.EQ.0)  GO  TO  13
      P1=DSQRT(PS(1)**2+PS(2)**2+PS(3)**2)
      CT1=PS(3)/P1
      AM1=AM(JP)
      E1=DSQRT(P1**2+AM1**2)
      CALL  CONSTR(P1,CT1,AM1,VND,B,JP,SK,STOT,RL ,N,IBA,3)
   13 CONTINUE
C
   14 IF(UDN.LE.(UMIN+0.010))   RETURN
      IF(NRAPI.NE.1)            RETURN
C     WRITE(16,*) 'IBA,UDN,UMIN=', IBA,UDN,UMIN
C
      CALL  ACTSTR(PDN,PN,SK,PS,JP,NSTR)
C     IF(NRAPI.EQ.2)  CALL  ACTSTR(PIN,PDN,SK,PS,JP,NSTR)
      IF(NSTR.EQ.0)   GO  TO  15
      P1=DSQRT(PS(1)**2+PS(2)**2+PS(3)**2)
      CT1=PS(3)/P1
      AM1=AM(JP)
      E1=DSQRT(P1**2+AM1**2)
      CALL  CONSTR(P1,CT1,AM1,VDN,B,JP,SK,STOT,RL ,N,IBA,2)
   15 CONTINUE
      RETURN
      END
C
C * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
      SUBROUTINE  ACTSTR(PIN,PN,SINS,PS,JSTR,NSTR)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL *8 PIN,PN,SINS,MA
      LOGICAL GH1H2
      CHARACTER*8 PNSTR(11),PNA1,PNA2,PNAJ,PNRES(6)
      COMMON/PARTCL/PPTCL(9,499),NPTCL,IORIG(499),IDENT(499)
     *,IDCAY(499)
      COMMON /H1H2/ GH1H2(11)
      COMMON /IDN12/ ID1,ID2
      COMMON /NCASCA/ NCAS,NCPRI
      COMMON /NRAPI/ NRAPI
      COMMON /COMENB/ ENBOU
      COMMON /ISRES/ ISRES
      DIMENSION  PIN(9),PN(9),PS(3),IDSTR(11),IDSRR(6)
      DATA PNSTR
     */'K+   ','K0   ','L    ','S+   ','S-   ','S0   ','K-   ','AP   ',
     *'PI- ','PI0 ','PI+ '/
      DATA PNRES
     */'ETA  ','ETAP ','RHO0 ','OMEGA','PHI  ','AP   '/
      DATA IDSTR
     */ 130,    230,    2130,   1130,   2230,   1230,  -130,   -1120,
     * -120,    110,     120/
      DATA IDSRR
     */ 220,    330,     111,    221,    331,  -1120/
C
      ENBOUT=3.8
      ENBOUS=ENBOU
      NSTR=0
      DO  10  I=1,11
   10 GH1H2(I)=.TRUE.
C   &&&&&&&&  18.01.92 &&&&&&&&&&&&&
C     GH1H2(4)=.FALSE.
C
      PX1=PIN(4)
      PY1=PIN(5)
      PZ1=PIN(6)
      PX2= PN(4)
      PY2= PN(5)
      PZ2= PN(6)
      AM1=PIN(9)
      AM2= PN(9)
      CALL  PANUID(ID1,IK1,PNA1)
      CALL  PANUID(ID2,IK2,PNA2)
      CALL  CROSEC(0,ID1,ID2,PX1,PY1,PZ1,AM1,PX2,PY2,PZ2,AM2,SIEL,0)
      CALL  CROSEC(1,ID1,ID2,PX1,PY1,PZ1,AM1,PX2,PY2,PZ2,AM2,SITO,0)
      CALL  CROSEC(3,ID1,ID2,PX1,PY1,PZ1,AM1,PX2,PY2,PZ2,AM2,SIAN,0)
      SIIN=SITO-SIEL
      IF(SITO.LE.0..OR.SIEL.LE.0.OR.SIIN.LE.0.)    RETURN
C    !!!!!
      ENBOU=ENBOUT
c     IF(NCAS.GE.NCPRI)   WRITE( *,*) 'ACTSTR IN'
C    !!!!!
      CALL COLLHH(ID1,AM1,PX1,PY1,PZ1,ID2,AM2,PX2,PY2,PZ2,
     *SITO,SIAN,SIEL)
c     IF(NCAS.GE.NCPRI)   WRITE( *,*) '      ACTSTR OUT'
      NP=0
C
      IF(NPTCL.LE.0)   GO  TO  15
      IF(ISRES.EQ.1.AND.JSTR.NE.6)   GO  TO  116
      DO 115  ND=1,5
      NDEC=0
      NPTCL1=NPTCL
      DO 114  I=1,NPTCL1
      IDI=IDENT(I)
      IF(IDI.EQ.1111.OR.IDI.EQ.1121.OR.IDI.EQ.1221.OR.IDI.EQ.2221)
     *                 GO  TO  114
      IF(IDI.EQ.110.OR.IDI.EQ.230.OR.IDI.EQ.-230)
     *                 GO  TO  114
      IF(IDI.EQ.1230.OR.IDI.EQ.-1230)
     *                 GO  TO  114
      CALL DECAY(I,IPOINT)
      IF(IPOINT.LT.0)  GO  TO  114
      NDEC=NDEC+1
      DO  113  J=1,9
      PPTCL(J,I)=PPTCL(J,NPTCL)
  113 CONTINUE
      IDENT(I)=IDENT(NPTCL)
      IORIG(I)=IORIG(NPTCL)
      IDCAY(I)=IDCAY(NPTCL)
      NPTCL=NPTCL-1
  114 CONTINUE
      IF(NDEC.EQ.0)     GO  TO  116
  115 CONTINUE
  116 CONTINUE
      DO  13  J=1,NPTCL
      MA=PPTCL(5,J)
      PS(1)=PPTCL(1,J)
      PS(2)=PPTCL(2,J)
      PS(3)=PPTCL(3,J)
      IDJ=IDENT(J)
      IF(ISRES.EQ.1)  GO  TO  117
      IF(JSTR.NE.2.AND.IDJ.EQ.IDSTR(JSTR))        GO  TO  14
      IF(JSTR.EQ.2.AND.IABS(IDJ).EQ.IDSTR(JSTR))  GO  TO  14
      GO  TO  13
  117 IF(IDJ.EQ.IDSRR(JSTR))        GO  TO  14
   13 CONTINUE
      GO  TO  15
   14 NSTR=1
      CALL  PANUID(IDJ,JP,PNAJ)
C    &&&&&&&&&&& 18.01.92 &&&&&&&&&&&&&&&&&
C     SINS=SIIN
      SINS=SITO
c     WRITE( *,1717)  PNA1,PNA2,PNAJ,IDJ,SINS,NRAPI
C     WRITE(16,1717)  PNA1,PNA2,PNAJ,IDJ,SINS,NRAPI
 1717 FORMAT(1X,3A5,'IDJ=',I8,' SINS=',F7.3,' NR=',I3)
   15 ENBOU=ENBOUS
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      FUNCTION GAMMP(A,X)
      IF(X.LT.0..OR.A.LE.0.) WRITE(16,*) 'GAMMP: X OR A LT.0'
      IF(X.LT.A+1.)THEN
      CALL GSER(GAMSER,A,X,GLN)
      GAMMP=GAMSER
      ELSE
      CALL GCF(GAMMCF,A,X,GLN)
      GAMMP=1.-GAMMCF
      ENDIF
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      FUNCTION GAMMLN(XX)
      REAL*8 COF(6),STP,HALF,ONE,FPF,X,TMP,SER
      DATA COF,STP/76.18009173D0,-86.50532033D0,24.01409822D0,
     *    -1.231739516D0,.120858003D-2,-.536382D-5,2.50662827465D0/
      DATA HALF,ONE,FPF/0.5D0,1.0D0,5.5D0/
      X=XX-ONE
      TMP=X+FPF
      TMP=(X+HALF)*LOG(TMP)-TMP
      SER=ONE
      DO 11 J=1,6
      X=X+ONE
      SER=SER+COF(J)/X
11    CONTINUE
      GAMMLN=TMP+LOG(STP*SER)
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE GCF(GAMMCF,A,X,GLN)
      PARAMETER (ITMAX=100,EPS=3.E-7)
      GLN=GAMMLN(A)
      GOLD=0.
      A0=1.
      A1=X
      B0=0.
      B1=1.
      FAC=1.
      DO 11 N=1,ITMAX
      AN=FLOAT(N)
      ANA=AN-A
      A0=(A1+A0*ANA)*FAC
      B0=(B1+B0*ANA)*FAC
      ANF=AN*FAC
      A1=X*A0+ANF*A1
      B1=X*B0+ANF*B1
      IF(A1.NE.0.)THEN
        FAC=1./A1
        G=B1*FAC
        IF(ABS((G-GOLD)/G).LT.EPS)GO TO 1
        GOLD=G
      ENDIF
11    CONTINUE
      WRITE(16,*) 'A too large, ITMAX too small'
1     GAMMCF=EXP(-X+A*ALOG(X)-GLN)*G
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE GSER(GAMSER,A,X,GLN)
      PARAMETER (ITMAX=100,EPS=3.E-7)
      GLN=GAMMLN(A)
      IF(X.LE.0.)THEN
      IF(X.LT.0.) WRITE(16,*) 'GSER: X.LT.0'
      GAMSER=0.
      RETURN
      ENDIF
      AP=A
      SUM=1./A
      DEL=SUM
      DO 11 N=1,ITMAX
      AP=AP+1.
      DEL=DEL*X/AP
      SUM=SUM+DEL
      IF(ABS(DEL).LT.ABS(SUM)*EPS)GO TO 1
11    CONTINUE
      WRITE(16,*) 'A too large, ITMAX too small'
1     GAMSER=SUM*EXP(-X+A*LOG(X)-GLN)
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
