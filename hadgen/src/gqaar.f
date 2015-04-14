
************************************************************************
*                                                                      *
*  The following subroutines are written by K.K.Gudima and V.D.Toneev  *
*  Joint Institute for Nuclear Research, Dubna,                        *
*  E-mail: gudima@acad.moldnet.md                                      *
*          gudima@cc.acad.md                                           *
*  E-mail: toneev@thsun1.jinr.ru                                       *
*                                                                      *
************************************************************************

C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE CASCAN(NEL,RN,DELTA,MV,IRET,INOVER)     ! KarSjh
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL *8 RN,DELTA
      COMMON/HCASC/ANUCL1,ANUCL2,ZNUCL1,ZNUCL2,T0,EPS1,EPS2,VPI,A1,A2,
     *C1,C2,D1,D2,R0N1,R0N2,TF01,TF02,RM1,RM2
      COMMON/RESULT/AN1,AN2,ZN1,ZN2,ENEXT1,ENEXT2,PNUCL1(3),PNUCL2(3),
     *AMNUC1(3),AMNUC2(3)
      COMMON/NUCSP/VPR(3),VTA(3),RADP(3),RADT(3),VEV(3),VRE(3),GEV,GRE
     *,VEP(3),VET(3),GEP,GET
      COMMON/ACTIV/MPA(240),MYP(5999),MYT(5999),MYY(5999)
      COMMON/CENTER/XC(2,240),YC(2,240),ZC(2,240),IZ(2,240)
       COMMON/MEMAGT/PMEMO(9,5999),IMEMO(5,5999)
     */KYPT/ KPT,KYP,KYT,KPT1,KYP1,KYT1,KYY,KYY1
      COMMON/EXCIT/TEX1,TEX2,HEX1,HEX2,PEX1,PEX2
      COMMON/INDSTR/JPSTR,NREP,JG
      COMMON/NRAPI/NRAPI
      COMMON/NUCOLL/ NUCOLL,MVCOLL(5999)
      COMMON/ISHOCK/ISHOCK/ACTIM/TINT
      COMMON/INTCC/INTCC
      COMMON/INTCEN/IPO,INT1(240),INT2(5999),INT3(5999),
     *              INT4(50000),IJPA,IST,JST,NRST
      COMMON/IPAUL/IP
      COMMON/TIMGO/TIMGO,MIN,ITMGO
      COMMON/INSTR/ INSTR
      COMMON/IPRIMC/ IPRIMC
      COMMON /INDRES/ INDRES
      COMMON/TLIMIT/TLIMIT
      COMMON /IDPME/ IDPME(5999)
      COMMON /ISLERR/ ISLERR   ! HADES1
      DIMENSION P1(9),P2(9),IP1(5),IP2(5),V12(3),P3(9),IP3(5),
     *PTEMP(3)
      COMMON/NCASCA/NCAS,NCPRI
C============== KarSjh ===============
      INOVER=0
      ICOTA=0
C============ end KarSjh =============
cc
      NCAS=10      ! 13.02.96
      NCPRI=999    ! 13.02.96
cc
      NUCOLL=0
      NEL=0
      IRET=0
    9 MVU=0
C============== KarSjh ===============
          ICOTA=ICOTA+1
          IF(ICOTA.gt.10000)THEN
            INOVER=1
            RETURN
          END IF
C============ end KarSjh =============
      AN1=ANUCL1
      ZN1=ZNUCL1
      AN2=ANUCL2
      ZN2=ZNUCL2
      ENEXT1=0.
      ENEXT2=0.
      PNUCL1(1)=0.
      PNUCL1(2)=0.
      PNUCL1(3)=0.
      PNUCL2(1)=0.
      PNUCL2(2)=0.
      PNUCL2(3)=0.
      AMNUC1(1)=0.
      AMNUC1(2)=0.
      AMNUC1(3)=0.
      AMNUC2(1)=0.
      AMNUC2(2)=0.
      AMNUC2(3)=0.
      TINT=0.
      MV=0
      KPT=0
      KYP=0
      KYT=0
      KYY=0
      KPT1=0
      KYP1=0
      KYT1=0
      KYY1=0
      TEX1=0.
      TEX2=0.
      PEX1=0.
      PEX2=0.
      HEX1=0.
      HEX2=0.
      DO 10 L=1,240
      MPA(L)=1
   10 CONTINUE
      CALL PINPN(RN,R0X,R0Y,R0Z,MV,DELTA,NEL)
      OBR1=0.
      IF(RM1.GT.0.1) OBR1=0.00144*ZN1/RM1
      OBR2=0.
      IF(RM2.GT.0.1) OBR2=0.00144*ZN2/RM2
      CALL  VINIT(VPR,VTA,ANUCL1,ANUCL2,T0)
      RADP(1)=R0X
      RADP(2)=R0Y
      RADP(3)=R0Z
      RADT(1)=0.
      RADT(2)=0.
      RADT(3)=0.
      IP=-1
      IJPA=0
   11 NA1=AN1+0.1
              NA2=AN2+0.1
      IPRIMC=0
c------------------------------
c       CONTINUE!WRITE(25,*)'CASCAN: MV=',MV
c       write(16,*)'CASCAN: MV=',MV
c       write(16,*)'CASCAN: ISHOCK=',ishock
c       if(mv.le.0)then
c        newmv=3
c       else
c        newmv=mv+2
c       endif
c       call priagt(newmv)
c------------------------------
cc      CALL  STOPGU(TIMGO,ISTOP)
C  *****  GANIL INSERT  *****
c      CALL TIMEX(CPUT)
c      RCPU=CPUT/60.D0
c      ISTOP=0
c      IF(RCPU.GT.TIMGO)   ISTOP=1
C  *****  GANIL INSERT  *****
cc      CALL  STOPP(ISTOPP)
cc      IF(ISTOP.NE.0.OR.ISTOPP.NE.0)   IRET=1
cc      IF(ISTOP.NE.0.OR.ISTOPP.NE.0)   RETURN
cc      IF(ISHOCK.EQ.0)  CALL  UPACOV(MV,MVU,0)
      CALL  DISSIP(IP,T0)
      CALL POINTN(M,P1,P2,IP1,IP2,V12,U12,T12,SIG,SAB,TINT,NR,N1,N2,
     *NA1,NA2,MV,DELTA,R0X,R0Y,R0Z)
      IF(M)   203,203,12 
  203 IF(MV.EQ.0)   GO  TO  205
      DO  204  L=1,MV
      MYP(L)=0
      MYT(L)=0
      MYY(L)=0
  204 CONTINUE
  205 CONTINUE
      IF(ENEXT2.LT.1.E-7)  GO  TO  21
cc      IF(ISHOCK.EQ.0)  CALL  UPACOV(MV,MVU,1)
cc      IF(ISHOCK.EQ.0)  MV=MVU
             GO  TO  20
   12 GO  TO  (13,14,13,112,15),NR
   13 NU=2
         GO TO 15
   14 NU=1
   15 CONTINUE
C
c-Sob      IF(T12.LE.0.)
c-Sob     *  write(16,1992) NR,N1,P1,IP1,N2,P2,IP2
c-Sob 1992 FORMAT(' T12=',E11.4,' NR=',I3/
c-Sob     *' N1=',I5,9(1X,F8.3),4(1X,I2),I10/
c-Sob     *' N2=',I5,9(1X,F8.3),4(1X,I2),I10)
      IF(T12.LE.0.) GO  TO  11
C
      CALL TYPNEW(P1,IP1,P2,IP2,V12,U12,T12,SAB,MV,NP,NABS,
     *P3,IP3,N3,NU,N2,NA1,NA2)
        if(ISLERR.eq.1)return   ! HADES1
c        if(ISLERR.eq.1)then
c          CONTINUE!WRITE(25,*)'      TYPNEW'
c          return
c        end if
      NRAPI=NR
C
C
      IP=0
      IF(NP)11,11,16
   16 GO TO  (17,18,19,112,113),NR
   17 CONTINUE
      CALL PAULI1(P1,P2,IP1,IP2,N1,N2,V12,NP,MV,TINT,IP)
C
      KPT1=KPT1+IP
      KPT=KPT+1
C                         IF(IP.EQ.0) GO TO 11
      GO  TO  111
   18 CONTINUE
      CALL PAULI2(P1,P2,P3,IP1,IP2,IP3,N1,N2,N3,V12,NP,MV,TINT,
     *IP,OBR1)
      KYP1=KYP1+IP
      KYP=KYP+1
C                        IF(IP.EQ.0) GO TO 11
      GO  TO  111
   19 CONTINUE
      CALL PAULI3(P1,P2,P3,IP1,IP2,IP3,N1,N2,N3,V12,NP,MV,TINT,
     *IP,OBR2)
      KYT1=KYT1+IP
      KYT=KYT+1
C                                 IF(IP.EQ.0) GO TO 11
      GO  TO  111
  111 CONTINUE
      IF(ANUCL1.LT.1.1.AND.NUCOLL.LE.1)  IPRIMC=1
      TCR=U12-P1(9)-P2(9)
      CALL TEVOLV(TCR,TINT,1)
      IF(INSTR.LE.2.AND.JPSTR.NE.0.AND.IP1(3).EQ.0)
     *CALL STRANN(U12,T12,V12,NU,N1,N2,P1,IP1,P2,IP2,NREP,JPSTR,SIG,SAB)
      IF(INSTR.GT.2.AND.INSTR.LT.5.AND.JPSTR.NE.0)
     *CALL STRANG(U12,V12,NU,N1,N2,P1,IP1,P2,IP2,NREP,JPSTR,SIG,SAB)
      IF(INSTR.GT.4.AND.JPSTR.NE.0.AND.IP1(3).EQ.0)
     *CALL STRANT(U12,V12,NU,N1,N2,P1,IP1,P2,IP2,NREP,JPSTR,SIG,SAB)
      GO  TO  11
  112 CONTINUE
* 10.11.93 T.   22.03.94
      IF(INDRES.EQ.1.AND.TINT.GT.TLIMIT.AND.
     &  (IDPME(N1).EQ.220.or.IDPME(N1).EQ.221.or.IDPME(N1).EQ.331)) then
      IP=0
      else
C
      CALL  DECREN(N1,NUR,MV,NP)
      CALL  PAULID(MV,N1,NUR,1,NP,IP)
      endif
      GO  TO  11
  113 CALL PAULI4(P1,P2,IP1,IP2,N1,N2,V12,NP,MV,TINT,IP)
      KYY1=KYY1+IP
      KYY=KYY+1
      TCR=U12-P1(9)-P2(9)
      CALL TEVOLV(TCR,TINT,2)
      IF(INSTR.GT.2.AND.INSTR.LT.5.AND.JPSTR.NE.0)
     *CALL STRANG(U12,V12,NU,N1,N2,P1,IP1,P2,IP2,NREP,JPSTR,SIG,SAB)
      IF(INSTR.GT.4.AND.JPSTR.NE.0)
     *CALL STRANT(U12,V12,NU,N1,N2,P1,IP1,P2,IP2,NREP,JPSTR,SIG,SAB)
C
      GO  TO  11
   20 IF(ENEXT2-0.0000001)  21,21,22
   21 NEL=NEL+1
      IF(ISHOCK.EQ.1)  LUP=1
              GO TO 9
   22 DO 23 L=1,3
      AMNUC1(L)=AMNUC1(L)*5.06

C   23 AMNUC2(L)=AMNUC2(L)*5.06  ! CHANGED By A.Timofeev 18/07/2011 18:24
   23 CONTINUE                   ! CHANGED By A.Timofeev 18/07/2011 18:24
      IF(ANUCL1-2.1) 24,26,26
   26 AM1=.94*AN1
      CALL CINEMA(PNUCL1,VPR,PTEMP,CT,ST,CF,SF,TTEMP,AM1)
      PNUCL1(1)=PTEMP(1)
      PNUCL1(2)=PTEMP(2)
      PNUCL1(3)=PTEMP(3)
   24 CONTINUE
      IF(ANUCL2.LT.2.1) GOTO 27
      AM2=.94*AN2
      CALL CINEMA(PNUCL2,VTA,PTEMP,CT,ST,CF,SF,TTEMP,AM2)
      PNUCL2(1)=PTEMP(1)
      PNUCL2(2)=PTEMP(2)
      PNUCL2(3)=PTEMP(3)
   27 CONTINUE
      RETURN
             END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *




C * * * * * * * * * * last correction  16.11.95 16:26 * * * * * * *
C Fermi energies are changed via TF=TF*(AN1/ANUCL1)**(2./3.)
c see PANTN, POTEN, TFERMI; 20.06.95
c anti-baryon potential is added 05.11.95
      SUBROUTINE POINTN(M,P1,P2,IP1,IP2,V1,U1,TR2,SIG,SAB,TINT,NR,N1,
     *N2,NA1,NA2,MV,DELTA,R0X,R0Y,R0Z)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 P1,P2,V1,U1,TR2,SIG,SAB,TINT,DELTA,R0X,R0Y,R0Z
      COMMON/INDINT/INDINT
      COMMON/RESULT/AN1,AN2,ZN1,ZN2,ENEXT1,ENEXT2,PNUCL1(3),PNUCL2(3),
     *AMNUC1(3),AMNUC2(3)
       COMMON/MEMAGT/PMEMO(9,5999),IMEMO(5,5999)
      COMMON/TAUE/TPTE,TYPE,TYTE
      COMMON/NUCSP/VPR(3),VTA(3),RADP(3),RADT(3),VEV(3),VRE(3),GEV,GRE
     *,VEP(3),VET(3),GEP,GET
      COMMON/DTINT/DTINT
      COMMON/HCASC/ANUCL1,ANUCL2,ZNUCL1,ZNUCL2,T00,EPS1,EPS2,
     *VPI,A1,A2,C1,C2,D1,D2,R0N1,R0N2,TF01,TF02,RM1,RM2
      COMMON/INTCEN/IPO,INT1(240),INT2(5999),INT3(5999),
     *              INT4(50000),IJPA,IST,JST,NRST
       COMMON/ISHOCK/ISHOCK
      COMMON/NCASCA/NCAS,NCPRI
      COMMON/PIDABS/PIDABS
      COMMON/INTCC/ INTCC
      COMMON/TAUIJ/ TAUK,TPTS,TYPS,TYTS,TYYS,TIJ4(50000)
      COMMON/TLIMIT/TLIMIT
      COMMON/CSLID/CLIDER(5999)
      COMMON/IACT/ IACT/CVALON/ IVALON
      COMMON /HOLPT/ PH1(3),PH2(3),RH1(3),RH2(3),EHOL1,EHOL2,UP,UT
      COMMON /IDPME/ IDPME(5999)
      COMMON /IDN12/ ID1,ID2
      COMMON /INDRES/ INDRES
      COMMON /BARPOT/ POT
      COMMON /SCOLLD/ SCOLLD
      DIMENSION P1(9),P2(9),IP1(5),IP2(5),V1(3)
     *,PT1(9),PT2(9),IPT1(5),IPT2(5),YP1(9),IYP1(5),YT1(9),IYT1(5)
      TFRO1=TINT
    8 CONTINUE
      R0=DSQRT(R0X**2+R0Y**2+R0Z**2)*1.01
      IPO=0
    9 TPT=-0.1
      TYP=-0.1
      TYT=-0.1
      TAUD=-0.1
      TYY=-0.1
      TAU=-0.1
      TY=-.1
      TFRO2=TINT
c      write( *,1991) NCAS,TINT
1991  FORMAT('+',50X,I6,1X,E10.3)
      RT=DSQRT((RADP(1)-RADT(1))**2+(RADP(2)-RADT(2))**2+
     +(RADP(3)-RADT(3))**2)
      IF(RT-R0)10,10,12
   10 IF(NA1.LT.1.OR.NA2.LT.1)   GO  TO  12
      IF(IPO.NE.1)               GO  TO  11
      IF(TPTS.LE.0.)             GO  TO  12
      IF((TPTS-TAUK).LE.0.0001)  GO  TO  11
      TPT=TPTS-TAUK
      GO  TO  12
   11 CONTINUE
      CALL  RAPID1(NA1,NA2,DELTA,PT1,IPT1,PT2,IPT2,NPT1,NPT2,TPT,DLPT)
   12  IF(MV) 124,124,100
  100 IF(NA1.LE.0)                             GO  TO  101
      IF(IPO.NE.1)                             GO  TO  13
      IF(TYPS.LE.0.)                           GO  TO  101
      IF((TYPS-TAUK).LE.0.0001)                GO  TO  13
      TYP=TYPS-TAUK
      IF(TYP.LT.PMEMO(7,NYP1).AND.CLIDER(NYP1).LT.0.3)
     *                                         GO  TO  13
      IF(TYP.LT.PMEMO(7,NYP1).AND.IVALON.EQ.0)
     *                                         GO  TO  13
      GO  TO  101
   13 CONTINUE
      CALL RAPID2(MV,NA1,DELTA,YP1,IYP1,NYP1,NYP2,TYP,DLYP)
  101 IF(NA2.LT.1)   GO  TO  14
  102 CONTINUE
      IF(IPO.NE.1)                             GO  TO  103
      IF(TYTS.LE.0.)                           GO  TO  14
      IF((TYTS-TAUK).LE.0.0001)                GO  TO  103
      TYT=TYTS-TAUK
      IF(TYT.LT.PMEMO(7,NYT1).AND.CLIDER(NYT1).LT.0.3)
     *                                         GO  TO  103
      IF(TYT.LT.PMEMO(7,NYT1).AND.IVALON.EQ.0)
     *                                         GO  TO  103
      GO  TO  14
  103 CONTINUE
      CALL RAPID3(MV,NA2,DELTA,YT1,IYT1,NYT1,NYT2,TYT,DLYT)
   14 CONTINUE
      CALL  RAPIDD(MV,TAUD,ND)
1992  CONTINUE
      IF(MV.LT.2.OR.INTCC.EQ.0.OR.TINT.GT.TLIMIT) GO  TO 124
      IF(IPO.NE.1)                               GO  TO  114
      IF(TYYS.LE.0.)                             GO  TO  124
      IF((TYY-TAUK).LE.0.0001)                   GO  TO  114
      TYY=TYYS-TAUK
      IF(TYY.LT.PMEMO(7,NYY1).AND.CLIDER(NYY1).LT.0.3)
     *                                            GO  TO 114
      IF(TYY.LT.PMEMO(7,NYY2).AND.CLIDER(NYY2).LT.0.3)
     *                                            GO  TO 114
      IF(TYY.LT.PMEMO(7,NYY1).AND.IVALON.EQ.0)
     *                                            GO  TO 114
      IF(TYY.LT.PMEMO(7,NYY2).AND.IVALON.EQ.0)
     *                                            GO  TO 114
      IF(NYY1.EQ.NYY2)                            GO  TO 114
      GO  TO  124
  114 CONTINUE
      CALL  RAPID4(MV,NYY1,NYY2,TYY)
  124 CONTINUE
      IF(TYP)19,19,15
   15 IF(TYT)16,16,17
   16 TY=TYP
      NR=2
      CLID1=CLIDER(NYP1)
      CLID2=1.
C
                               GO  TO  20
   17 IF(TYP-TYT) 16,16,18
   18 TY=TYT
      NR=3
      CLID1=CLIDER(NYT1)
      CLID2=1.
                              GO TO 20
   19 IF(TYT)20,20,18
   20 IF(TPT)24,24,21
   21 IF(TY)23,23,22
   22 IF(TPT-TY)  23,23,25
   23 TAU=TPT
      NR=1
      CLID1=1.
      CLID2=1.
C
                            GO  TO   26
   24 IF(TY)26,26,25
   25 TAU=TY
C
                            GO  TO   26
   26 IF(TAUD.LT.0.)                GO  TO  126
      IF(TAUD.GT.TAU.AND.TAU.GE.0.) GO  TO  126
      TAU=TAUD
      N1=ND
      NR=4
  126 IF(TYY.LT.0.)                 GO  TO  226
      IF(TYY.GT.TAU.AND.TAU.GE.0.)  GO  TO  226
      TAU=TYY
      N1=NYY1
      N2=NYY2
      NR=5
      CLID1=CLIDER(NYY1)
      CLID2=CLIDER(NYY2)
  226 IF(TAU.LT.0.)                 GO  TO  41
C
c-Sob      IF(NCAS.GE.NCPRI)  THEN
c-Sob      write(16,6000) MV,IPO,NR,TPT,TYP,TYT,TYY,TAUD,TAU
c-Sob 6000 FORMAT(1X,'MV,IPO,NR,TPT,TYP,TYT,TYY,TAUD,TAU=',I4,2I2,6(1X,F6.3))
c-Sob      ENDIF
C
      TAUK=TAU
      IF(MV)29,29,27
   27 DO 28 K=1,MV
      EK=PMEMO(8,K)+PMEMO(9,K)
      PMEMO(1,K)=PMEMO(1,K)+PMEMO(4,K)/EK*TAUK
      PMEMO(2,K)=PMEMO(2,K)+PMEMO(5,K)/EK*TAUK
      PMEMO(3,K)=PMEMO(3,K)+PMEMO(6,K)/EK*TAUK
      PMEMO(7,K)=PMEMO(7,K)-TAUK
      IF(PMEMO(7,K).GT.0.)          GO  TO  128
      PMEMO(7,K)=0.
      CLIDER(K)=1.
  128 CONTINUE
      IF(IMEMO(5,K).EQ.0)           GO  TO  28
      IMEMO(5,K)=IMEMO(5,K)-IDINTG(TAUK*1000.)
      IF(IMEMO(5,K).LE.0)  IMEMO(5,K)=1
   28 CONTINUE
   29 TINT=TINT+TAU
      IF(IJPA.LE.0)                 GO  TO  291
      DO  290  IJ=1,IJPA
  290 TIJ4(IJ)=TIJ4(IJ)-TAUK
  291 CONTINUE
      DO  30  K=1,3
      RADT(K)=RADT(K)+TAU*VTA(K)
      RADP(K)=RADP(K)+TAU*VPR(K)
   30 CONTINUE
      GO  TO (33,35,37,130,230),NR
C  ====> RAPIDD
  130 IF(ISHOCK.EQ.1)  CALL  SHOCKF(TFRO2,TINT,NA1,NA2,MV,0)
      M=1
      GO  TO  135
C  ====> RAPID4
  230 DO 231  K=1,9
      P1(K)=PMEMO(K,N1)
      P2(K)=PMEMO(K,N2)
      IF(K.LE.5)   IP1(K)=IMEMO(K,N1)
      IF(K.LE.5)   IP2(K)=IMEMO(K,N2)
  231 CONTINUE
      ID1=IDPME(N1)
      ID2=IDPME(N2)
      GO  TO  32
C  ====> RAPID1
   33 DO  34  K=1,9
      P1(K)=PT1(K)
      P2(K)=PT2(K)
      IF(K.LE.5)   IP1(K)=IPT1(K)
      IF(K.LE.5)   IP2(K)=IPT2(K)
   34 CONTINUE
      N1=NPT1
      N2=NPT2
      DLK=DLPT
      ID1=1120
      IF(IP1(1).EQ.0) ID1=1220
      ID2=1120
      IF(IP2(1).EQ.0) ID2=1220
      GO  TO   32
C  ====> RAPID2
   35 DO  36  K=1,9
      P1(K)=YP1(K)
      IF(K.LE.5)   IP1(K)=IYP1(K)
   36 CONTINUE
      N1=NYP1
      N2=NYP2
      DLK=DLYP
                 GO  TO  39
C  ====> RAPID3
   37 DO  38  K=1,9
      P1(K)=YT1(K)
      IF(K.LE.5)   IP1(K)=IYT1(K)
   38 CONTINUE
      N1=NYT1
      N2=NYT2
      DLK=DLYT
      CALL  PARTN(2,N2,P2,IP2)
      UT=0.940-P2(9)+P2(8)
      EHOL2=P2(8)
      P2(9)=0.94
      POT=0.
      IF(IP1(3).EQ.0.AND.IP1(5).EQ.0) then
      if(IP1(4).EQ.1)   POT=UT
      if(IP1(4).EQ.0)   POT=VPI
      if(IP1(4).EQ.-1)  POT=POTBAR(2,P2(1),P2(2),P2(3))
      ENDIF
      P1(8)=P1(8)+POT
      C=DSQRT(P1(8)/(P1(8)-POT)*(P1(8)+2.*P1(9))/(P1(8)-POT+2.*P1(9)))
      DO  339 K=1,3
      P1(K+3)=P1(K+3)*C
      PH2(K)=-P2(3+K)
  339 RH2(K)=P2(K)
      ID1=IDPME(N1)
      ID2=1120
      IF(IP2(1).EQ.0) ID2=1220
      GO  TO  32
   39 CONTINUE
      CALL  PARTN(1,N2,P2,IP2)
      UP=0.940-P2(9)+P2(8)
      EHOL1=P2(8)
      P2(9)=0.94
      POT=0.
      IF(IP1(3).EQ.0.AND.IP1(5).EQ.0) then
      if(IP1(4).EQ.1)   POT=UP
      if(IP1(4).EQ.0)   POT=VPI
      if(IP1(4).EQ.-1)  POT=POTBAR(1,P2(1),P2(2),P2(3))
      ENDIF
      P1(8)=P1(8)+POT
      C=DSQRT(P1(8)/(P1(8)-POT)*(P1(8)+2.*P1(9))/(P1(8)-POT+2.*P1(9)))
      DO  239 K=1,3
      P1(K+3)=P1(K+3)*C
      PH1(K)=-P2(3+K)
  239 RH1(K)=P2(K)
      ID1=IDPME(N1)
      ID2=1120
      IF(IP2(1).EQ.0) ID2=1220
      GO  TO  32
   32 CONTINUE
C
      PX1=P1(4)
      PY1=P1(5)
      PZ1=P1(6)
      AM1=P1(9)
      PX2=P2(4)
      PY2=P2(5)
      PZ2=P2(6)
      AM2=P2(9)
      CALL CROSEC(1,ID1,ID2,PX1,PY1,PZ1,AM1,PX2,PY2,PZ2,AM2,SITO,0)
      SIG=SITO
      CALL  SLQEK(L,MS,MQ,KSIN,ME,IP1,IP2)
      CALL  TINVU(P1,P2,U1,V1,TR2)
C  *******  CM ENERGY FOR MDELT1   14.09.92 *********
      SCOLLD=U1*U1
      SAB=0.
      GO  TO  (135,131,132,135,135),NR
  131 IF(AN1.LT.2.1.OR.ZN1.LT.1.1)          GO  TO  135
      GO  TO  133
  132 IF(AN2.LT.2.1.OR.ZN2.LT.1.1)          GO  TO  135
  133 CONTINUE
      IF(ID1.EQ.120.OR.ID1.EQ.-120.OR.ID1.EQ.110)  GO  TO  134
                                       GO  TO  135
  134 CONTINUE
C   !!!  ONLY FOR PION(+-0)  !!!
      SAB=CROSEG(L,MS,MQ,KSIN,3,P1(8),P1(9),IP1(5))
      SAB=SAB*PIDABS
  135 CONTINUE
      IPO=1
      IST=N1
      JST=N2
      NRST=NR
      TPTS=TPT
      TYPS=TYP
      TYTS=TYT
      TYYS=TYY
      IF(NR.EQ.4)        GO  TO  42
      SIGV=SIG
      IF(IVALON.NE.0)    SIGV=CLID1*CLID2*SIG
      IF(ISHOCK.EQ.1)  CALL SHOCKF(TFRO2,TINT,NA1,NA2,MV,0)
      IF(INDINT.EQ.2)    GO  TO  139
      IF(NR.EQ.5)        GO  TO  42
      TEMP1=1.-SIGV/(31.41592*(DLK**2))
      DRND=RNDMD(-1)
      IF(DRND.LT.TEMP1)  GO  TO  9
                   GO  TO  42
  139 TEMP1=SIGV/31.41592
      BI2=BIM2(P1,P2)
      IF(BI2.GT.TEMP1)   GO  TO  9
                   GO  TO  42
   41 M=0
      TSH2=0.
      IF(ISHOCK.EQ.1)  CALL SHOCKF(TINT,TSH2,NA1,NA2,MV,0)
c-Sob      IF(NCAS.GE.NCPRI)  write( *,*) ' M=0'
      RETURN
   42 CONTINUE
      M=1
      DTINT=TINT-TFRO1
      IF(NCAS.LT.NCPRI)  RETURN
c-Sob      write(16,665) IJPA,NR,TAUK,TINT
c-Sob 665  FORMAT(1X,'IJPA,NR,TAUK,TINT=',I5,I2,2(1X,F7.3))
C  *****************
C     IF(IJPA.EQ.0)  GO TO 1665
C     DO 1666 KIJ=1,IJPA
C1666 write(16,*) 'M,INT4(M)=', KIJ,INT4(KIJ)
C1665 CONTINUE
C  *****************
c-Sob      IF(NR.GE.2.AND.NR.LE.4)
c-Sob     *write(16,666) NR,N1,(PMEMO(K,N1),K=1,9),(IMEMO(K,N1),K=1,5)
c-Sob 666  FORMAT(' NR,N1=',I1,I5,9(1X,F10.3),4I2,I8)
c-Sob      IF(NR.EQ.1.OR.NR.EQ.5)
c-Sob     *write(16,667) NR,N1,(P1(K),K=1,9),(IP1(K),K=1,5),
c-Sob     *                 N2,(P2(K),K=1,9),(IP2(K),K=1,5)
c-Sob 667  FORMAT(' NR,N1=',I1,I5,9(1X,F10.3),4I2,I8/
c-Sob     *       '   +N2=',1X,I5,9(1X,F10.3),4I2,I8)
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE RAPID1(NA1,NA2,DELTA,P1,IP1,P2,IP2,N1,N2,TAU,DL1)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL *8 DELTA,P1,P2,TAU,DL1
      COMMON/RCOR/ RCOR
      COMMON/NCASCA/NCAS,NCPRI
      COMMON/INTCEN/IPO,INT1(240),INT2(5999),INT3(5999),
     *              INT4(50000),IJPA,IST,JST,NRST
      COMMON/ACTIV/MPA(240),MYP(5999),MYT(5999),MYY(5999)
      COMMON/CENTER/XC(2,240),YC(2,240),ZC(2,240),IZ(2,240)
      COMMON /HOLPT/ PH1(3),PH2(3),RH1(3),RH2(3),EHOL1,EHOL2,UP,UT
      COMMON/NUCSP/VPR(3),VTA(3),RADP(3),RADT(3),VEV(3),VRE(3),GEV,GRE
     *,VEP(3),VET(3),GEP,GET
      DIMENSION P1(9),P2(9),IP1(5),IP2(5),RI(3),RJ(3),RIJ(3),VIJ(3),
     *PI(3),PJ(3),R1(3),R2(3),PLI(3),PLJ(3)
      DL1=1./(5.06*0.940*DSQRT(GRE*GRE-1.))+DELTA
      GPR=1./DSQRT(DABS(1.-VPR(1)**2-VPR(2)**2-VPR(3)**2))
      GTA=1./DSQRT(DABS(1.-VTA(1)**2-VTA(2)**2-VTA(3)**2))
      DO   8  K=1,3
      PI(K)=0.940*GPR*VPR(K)
      PJ(K)=0.940*GTA*VTA(K)
      VIJ(K)=VPR(K)-VTA(K)
    8 CONTINUE
      VIJ2=VIJ(1)**2+VIJ(2)**2+VIJ(3)**2
      EI=0.940*GPR
      EJ=0.940*GTA
      SIJ=(EI+EJ)**2-(PI(1)+PJ(1))**2-(PI(2)+PJ(2))**2-
     -(PI(3)+PJ(3))**2
      TAU=-0.1
      DO  13  I=1,NA1
      IF(MPA(I).EQ.0)   GO  TO  13
      NI=0
      VRI=XC(1,I)*VPR(1)+YC(1,I)*VPR(2)+ZC(1,I)*VPR(3)
      RI(1)=XC(1,I)-VPR(1)*VRI*GPR/(GPR+1.)+RADP(1)
      RI(2)=YC(1,I)-VPR(2)*VRI*GPR/(GPR+1.)+RADP(2)
      RI(3)=ZC(1,I)-VPR(3)*VRI*GPR/(GPR+1.)+RADP(3)
      DO  12  J=1,NA2
      VRJ=XC(2,J)*VTA(1)+YC(2,J)*VTA(2)+ZC(2,J)*VTA(3)
      RJ(1)=XC(2,J)-VTA(1)*VRJ*GTA/(GTA+1.)+RADT(1)
      RJ(2)=YC(2,J)-VTA(2)*VRJ*GTA/(GTA+1.)+RADT(2)
      RJ(3)=ZC(2,J)-VTA(3)*VRJ*GTA/(GTA+1.)+RADT(3)
      DO  9  K=1,3
      RIJ(K)=RI(K)-RJ(K)
    9 CONTINUE
      TIJ=-(RIJ(1)*VIJ(1)+RIJ(2)*VIJ(2)+RIJ(3)*VIJ(3))/VIJ2
      IF(TIJ.LE.0.)      GO  TO  12
      SP=RIJ(1)*PI(1)+RIJ(2)*PI(2)+RIJ(3)*PI(3)
      ST=RIJ(1)*PJ(1)+RIJ(2)*PJ(2)+RIJ(3)*PJ(3)
      RIJ2=RIJ(1)**2+RIJ(2)**2+RIJ(3)**2
      B2=RIJ2+4.*(SIJ*SP*ST-(0.940*(SP+ST))**2)/
     /SIJ/(SIJ-4.*0.940**2)
      IF(B2.GT.(DL1**2))          GO  TO  12
      IF(IPO.EQ.1.AND.NRST.EQ.1.AND.I.EQ.IST.AND.J.EQ.JST)  GO  TO 12
      NI=NI+1
      IF(TAU.LT.0.)        GO  TO  10
      IF(TAU.LT.TIJ)       GO  TO  12
   10 TAU=TIJ
      N1=I
      N2=J
      DO  11  K=1,3
      R1(K)=RI(K)+TIJ*VPR(K)
      R2(K)=RJ(K)+TIJ*VTA(K)
   11 CONTINUE
   12 CONTINUE
      IF(NI.EQ.0)   MPA(I)=0
   13 CONTINUE
      IF(TAU.LE.0)  RETURN
      CALL  PARTN(1,N1,P1,IP1)
      UP=0.940-P1(9)+P1(8)
      EHOL1=P1(8)
      P12=P1(4)**2+P1(5)**2+P1(6)**2
      C1=1.
      IF(P12.GT.1.D-10) C1=DSQRT((P1(8)*(P1(8)+2.*P1(9)))/P12)
      CALL  PARTN(2,N2,P2,IP2)
      UT=0.940-P2(9)+P2(8)
      EHOL2=P2(8)
      P22=P2(4)**2+P2(5)**2+P2(6)**2
      C2=1.
      IF(P22.GT.1.D-10) C2=DSQRT((P2(8)*(P2(8)+2.*P2(9)))/P22)
      DO  14  K=1,3
      PI(K)=P1(3+K)*C1
      PJ(K)=P2(3+K)*C2
   14 CONTINUE
      CALL  KINEMA(PI,VPR,PLI,CTI,STI,CFI,SFI,TLI,P1(9))
      CALL  KINEMA(PJ,VTA,PLJ,CTJ,STJ,CFJ,SFJ,TLJ,P2(9))
      DO  15  K=1,3
      RH1(K)=P1(K)
      RH2(K)=P2(K)
      P1(K) =R1(K)
      P2(K) =R2(K)
      PH1(K)=-P1(3+K)
      PH2(K)=-P2(3+K)
      P1(3+K)=PLI(K)
      P2(3+K)=PLJ(K)
   15 CONTINUE
      P1(7)=0.
      P2(7)=0.
      P1(8)=DSQRT(PLI(1)**2+PLI(2)**2+PLI(3)**2+0.940**2)-0.940
      P2(8)=DSQRT(PLJ(1)**2+PLJ(2)**2+PLJ(3)**2+0.940**2)-0.940
      P1(9)=0.940
      P2(9)=0.940
      RETURN
             END
C * * * * * * *  * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE RAPID2(MV,NA1,DELTA,P1,IP1,N1,N2,TAU,DL)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL *8 DELTA,P1,TAU,DL
      COMMON/RCOR/ RCOR
      COMMON/NCASCA/NCAS,NCPRI
      COMMON/HCASC/ANUCL1,ANUCL2,ZNUCL1,ZNUCL2,T0,EPS1,EPS2,
     *VPI,A1,A2,C1,C2,D1,D2,R0N1,R0N2,TF01,TF02,RM1,RM2
      COMMON/ACTIV/MPA(240),MYP(5999),MYT(5999),MYY(5999)
       COMMON/MEMAGT/PMEMO(9,5999),IMEMO(5,5999)
      COMMON/CENTER/XC(2,240),YC(2,240),ZC(2,240),IZ(2,240)
      COMMON/CENPAR/NUKC(100)
       COMMON/TAUE/TPTE,TYPE,TYTE
      COMMON/NUCSP/VPR(3),VTA(3),RADP(3),RADT(3),VEV(3),VRE(3),GEV,GRE
     *,VEP(3),VET(3),GEP,GET
      COMMON/INTCEN/IPO,INT1(240),INT2(5999),INT3(5999),
     *              INT4(50000),IJPA,IST,JST,NRST
      COMMON/CSLID/CLIDER(5999)
      COMMON/IACT/ IACT/CVALON/ IVALON
      DIMENSION V0(3),PIN(9),IIN(5),CS(3),C(3),P1(9),IP1(5)
      K = 1
      TAU = -.1
      TAUL=-.1
      V0(1)=-VPR(1)
      V0(2)=-VPR(2)
      V0(3)=-VPR(3)
      GPR=1./DSQRT(1.-VPR(1)**2-VPR(2)**2-VPR(3)**2)
      GG=GPR*GPR/(GPR+1.)
    9 IF(MYP(K)) 26,26,10
   10 VRK=(PMEMO(1,K)-RADP(1))*VPR(1)+(PMEMO(2,K)-RADP(2))*VPR(2)+
     +(PMEMO(3,K)-RADP(3))*VPR(3)
      XK0=PMEMO(1,K)-RADP(1)+VPR(1)*VRK*GG
      YK0=PMEMO(2,K)-RADP(2)+VPR(2)*VRK*GG
      ZK0=PMEMO(3,K)-RADP(3)+VPR(3)*VRK*GG
      CS(1)=PMEMO(4,K)
      CS(2)=PMEMO(5,K)
      CS(3)=PMEMO(6,K)
      CALL CINEMA(CS,V0,C,CT,ST,CF,SF,TL,PMEMO(9,K))
      PIN(4)=C(1)
      PIN(5)=C(2)
      PIN(6)=C(3)
      PIN(7)=PMEMO(7,K)
      PIN(9)=PMEMO(9,K)
      IIN(1)=IMEMO(1,K)
      IIN(2)=IMEMO(2,K)
      IIN(3)=IMEMO(3,K)
      IIN(4)=IMEMO(4,K)
      IIN(5)=IMEMO(5,K)
      CLID=CLIDER(K)
      PIN(1)=XK0
      PIN(2)=YK0
      PIN(3)=ZK0
      VLI=DSQRT(TL*(TL+2.*PIN(9)))/(TL+PIN(9))
      PIN(8)=TL
      DLK=1./(5.06*DSQRT(PIN(8)*(PIN(8)+2.*PIN(9))))+DELTA
      ICE=0
      IF(IPO.EQ.0)   GO  TO  12
      K1=INT2(K)
      GO  TO  16
   12 IF(ICE.NE.0)   GO  TO  13
      CALL  CENUM1(NA1,PIN,IIN,DLK,NC,K1,K2,0,1)
      ICE=1
   13 IF(NC.NE.0)   GO  TO  15
      MYP(K)=0
              GO  TO  26
   15 K1=NUKC(ICE)
      ICE=ICE+1
      NC=NC-1
   16 CONTINUE
      DR=((XC(1,K1)-PIN(1))*CF+(YC(1,K1)-PIN(2))*SF)*ST
     *+(ZC(1,K1)-PIN(3))*CT
      IF(DR.LE.0.)    GO TO  12
      IF(IPO.EQ.1.AND.NRST.EQ.2.AND.K.EQ.IST.AND.K1.EQ.JST) GO TO 12
      PIN(1)=PIN(1)+DR*ST*CF
      PIN(2)=PIN(2)+DR*ST*SF
      PIN(3)=PIN(3)+DR*CT
      DX=PIN(1)-XK0
      DY=PIN(2)-YK0
      DZ=PIN(3)-ZK0
      DRK=DSQRT(DX**2+DY**2+DZ**2)
      TAUK=DRK/VLI
      TAUK1=TAUK*GEP*(1.-VLI*(VEP(1)*ST*CF+VEP(2)*ST*SF+VEP(3)+CT))
      TAUKL=TAUK*GPR*(1.+VLI*(VPR(1)*ST*CF+VPR(2)*ST*SF+VPR(3)*CT))
C    !!!!!
      IF(TAUKL.LT.PIN(7).AND.CLID.LT.0.3)   GO  TO  12
      IF(TAUKL.LT.PIN(7).AND.IVALON.EQ.0)   GO  TO  12
C    !!!!!
      INT2(K)=K1
      IF(TAU)  23,23,22
   22 IF(TAU-TAUK)  26,26,23
   23 TAU=TAUK
      TYPE=TAUK1
      TAUL=TAUKL
      DO 24 L=1,9
      P1(L)=PIN(L)
   24 CONTINUE
      DO 25 L=1,5
      IP1(L)=IIN(L)
   25 CONTINUE
      DL=DLK
      N1=K
      N2=K1
   26 IF(K-MV) 27,28,28
   27 K=K+1
            GO TO 9
   28 CONTINUE
      TAU=TAUL
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE RAPID3(MV,NA2,DELTA,P1,IP1,N1,N2,TAU,DL)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL *8 DELTA,P1,TAU,DL
      COMMON/RCOR/ RCOR
      COMMON/NCASCA/NCAS,NCPRI
      COMMON/HCASC/ANUCL1,ANUCL2,ZNUCL1,ZNUCL2,T0,EPS1,EPS2,
     *VPI,A1,A2,C1,C2,D1,D2,R0N1,R0N2,TF01,TF02,RM1,RM2
      COMMON/NUCSP/VPR(3),VTA(3),RADP(3),RADT(3),VEV(3),VRE(3),GEV,GRE
     *,VEP(3),VET(3),GEP,GET
      COMMON/ACTIV/MPA(240),MYP(5999),MYT(5999),MYY(5999)
       COMMON/MEMAGT/PMEMO(9,5999),IMEMO(5,5999)
      COMMON/CENTER/XC(2,240),YC(2,240),ZC(2,240),IZ(2,240)
      COMMON/CENPAR/NUKC(100)
       COMMON/TAUE/TPTE,TYPE,TYTE
      COMMON/INTCEN/IPO,INT1(240),INT2(5999),INT3(5999),
     *              INT4(50000),IJPA,IST,JST,NRST
      COMMON/CSLID/CLIDER(5999)
      COMMON/IACT/ IACT/CVALON/ IVALON
      DIMENSION V0(3),PIN(9),IIN(5),CS(3),C(3),P1(9),IP1(5)
      V0(1)=-VTA(1)
      V0(2)=-VTA(2)
      V0(3)=-VTA(3)
      GTA=1./DSQRT(1.-VTA(1)**2-VTA(2)**2-VTA(3)**2)
      GG=GTA*GTA/(GTA+1.)
      K = 1
      TAU = -.1
      TAUL=-.1
    9 IF(MYT(K)) 26,26,10
   10 VRK=(PMEMO(1,K)-RADT(1))*VTA(1)+(PMEMO(2,K)-RADT(2))*VTA(2)+
     +(PMEMO(3,K)-RADT(3))*VTA(3)
      XK0=PMEMO(1,K)-RADT(1)+VTA(1)*VRK*GG
      YK0=PMEMO(2,K)-RADT(2)+VTA(2)*VRK*GG
      ZK0=PMEMO(3,K)-RADT(3)+VTA(3)*VRK*GG
      CS(1)=PMEMO(4,K)
      CS(2)=PMEMO(5,K)
      CS(3)=PMEMO(6,K)
      CALL  CINEMA(CS,V0,C,CT,ST,CF,SF,TL,PMEMO(9,K))
      PIN(4)=C(1)
      PIN(5)=C(2)
      PIN(6)=C(3)
      PIN(7)=PMEMO(7,K)
      PIN(9)=PMEMO(9,K)
      PIN(8)=TL
      IIN(1)=IMEMO(1,K)
      IIN(2)=IMEMO(2,K)
      IIN(3)=IMEMO(3,K)
      IIN(4)=IMEMO(4,K)
      IIN(5)=IMEMO(5,K)
      CLID=CLIDER(K)
      PIN(1)=XK0
      PIN(2)=YK0
      PIN(3)=ZK0
      VLI=DSQRT(TL*(TL+2.*PIN(9)))/(TL+PIN(9))
      PIN(8)=TL
      DLK=1./(5.06*DSQRT(PIN(8)*(PIN(8)+2.*PIN(9))))+DELTA
      ICE=0
      IF(IPO.EQ.0)   GO  TO  12
      K1=INT3(K)
      GO  TO  16
   12 IF(ICE.NE.0)   GO  TO  13
      CALL  CENUM1(NA2,PIN,IIN,DLK,NC,K1,K2,0,2)
      ICE=1
   13 IF(NC.NE.0)   GO  TO  15
      MYT(K)=0
               GO  TO  26
   15 K1=NUKC(ICE)
      ICE=ICE+1
      NC=NC-1
   16 CONTINUE
      DR=((XC(2,K1)-PIN(1))*CF+(YC(2,K1)-PIN(2))*SF)*ST+
     *(ZC(2,K1)-PIN(3))*CT
      IF(DR.LE.0.)    GO TO  12
      IF(IPO.EQ.1.AND.NRST.EQ.3.AND.K.EQ.IST.AND.K1.EQ.JST) GO TO 12
      PIN(1)=PIN(1)+DR*ST*CF
      PIN(2)=PIN(2)+DR*ST*SF
      PIN(3)=PIN(3)+DR*CT
      DX=PIN(1)-XK0
      DY=PIN(2)-YK0
      DZ=PIN(3)-ZK0
      DRK=DSQRT(DX**2+DY**2+DZ**2)
      TAUK=DRK/VLI
      TAUK1=TAUK*GET*(1.-VLI*(VET(1)*ST*CF+VET(2)*ST*SF+VET(3)*CT))
      TAUKL=TAUK*GTA*(1.+VLI*(VTA(1)*ST*CF+VTA(2)*ST*SF+VTA(3)*CT))
C    !!!!!
      IF(TAUKL.LT.PIN(7).AND.CLID.LT.0.3)   GO  TO  12
      IF(TAUKL.LT.PIN(7).AND.IVALON.EQ.0)   GO  TO  12
C    !!!!!
      INT3(K)=K1
      IF(TAU)  23,23,22
   22 IF(TAU-TAUK) 26,26,23
   23 TAU=TAUK
      TYTE=TAUK1
      TAUL=TAUKL
      DO 24 L=1,9
      P1(L)=PIN(L)
   24 CONTINUE
      DO 25 L=1,5
      IP1(L)=IIN(L)
   25 CONTINUE
      DL=DLK
      N1=K
      N2=K1
   26 IF(K-MV) 27,28,28
   27 K=K+1
            GO TO 9
   28 CONTINUE
      TAU=TAUL
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE  RAPIDD(MV,TAU,MD)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL *8 TAU
      COMMON /MEMAGT/ PME(9,5999),IME(5,5999)
      COMMON/TLIMIT/TLIMIT
      COMMON/ACTIM/TINT
      COMMON/IDPME/ IDPME(5999)
      COMMON/INDRES/ INDRES
      TAU=-.1
      IF(MV.LE.0)   RETURN
      DO  11  M=1,MV
      IF(IME(5,M).EQ.0)  GO  TO  11
C   @@@@@@@@@@@@@@  FOR ETA-MESONS 22.03.94 @@@@ 15.02.93 @@@@@@@@@@@@
      IF(INDRES.EQ.1.AND.TINT.GT.TLIMIT.AND.
     &(IDPME(M).EQ.220.or.IDPME(M).EQ.221.or.IDPME(M).EQ.331))
     *                                           GO  TO  11
C   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      TAUM=DFLOAT(IME(5,M))/1000.
      IF(TAU.LT.0.)  GO  TO  10
      IF(TAUM.GT.TAU)  GO  TO  11
   10 TAU=TAUM
      MD=M
   11 CONTINUE
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE  RAPID4(MV,N1,N2,TAU)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL *8 TAU
      COMMON/NCASCA/NCAS,NCPRI
      COMMON /MEMAGT/ PME(9,5999),IME(5,5999)
      COMMON/INTCEN/IPO,INT1(240),INT2(5999),INT3(5999),
     *              INT4(50000),IJPA,IST,JST,NRST
      COMMON/TAUIJ/ TAUK,TPTS,TYPS,TYTS,TYYS,TIJ4(50000)
      COMMON/NUCOLL/ NUCOLL,MVCOLL(5999)
      COMMON/ACTIV/MPA(240),MYP(5999),MYT(5999),MYY(5999)
      COMMON/TLIMIT/TLIMIT
      COMMON/ACTIM/TINT
      COMMON/IPAUL/IP
      COMMON/INTCC/INTCC
      COMMON/IACT/ IACT/CVALON/ IVALON
      COMMON/CSLID/CLIDER(5999)
      DIMENSION RIJ(3),VIJ(3),NYY(5999)
      TAU=-.1
      IF(MV.LE.1)                                           RETURN
      IF(IP.EQ.0.OR.IPO.EQ.1)                               GO  TO  15
      DO  8  I=1,MV
    8 NYY(I)=MYY(I)
      DO  13 I=2,MV
      CLIDI=CLIDER(I)
      EI=PME(8,I)+PME(9,I)
      INU=0
      IF( IME(2,I).EQ.0.AND.IME(3,I).EQ.0.AND.IME(4,I).EQ.1.
     *AND.IME(5,I).EQ.0)  INU=1
      JM=I-1
      DO  12  J=1,JM
      IF(MVCOLL(I).EQ.MVCOLL(J))                            GO  TO  12
      IF(MYY(I).EQ.1)                                       GO  TO  9
      IF(MYY(J).EQ.1)                                       GO  TO  9
                                              GO  TO  12
    9 JNU=0
      IF( IME(2,J).EQ.0.AND.IME(3,J).EQ.0.AND.IME(4,J).EQ.1.
     *AND.IME(5,J).EQ.0)  JNU=1
      IF((INU+JNU).LT.1.AND.INTCC.EQ.1)                     GO  TO  12
      CLIDJ=CLIDER(J)
      EJ=PME(8,J)+PME(9,J)
      DO  10  K=1,3
      RIJ(K)=PME(K,I)-PME(K,J)
      VIJ(K)=PME(K+3,I)/EI-PME(K+3,J)/EJ
   10 CONTINUE
      RVIJ=RIJ(1)*VIJ(1)+RIJ(2)*VIJ(2)+RIJ(3)*VIJ(3)
      IF(RVIJ.GE.0.)                                        GO  TO  12
      VIJ2=VIJ(1)**2+VIJ(2)**2+VIJ(3)**2
      TIJ=-RVIJ/VIJ2
      IF(TIJ.LT.0.0000001.OR.(TINT+TIJ).GT.TLIMIT)          GO  TO  12
C   !!!!!
      IF(TIJ.LT.PME(7,I).AND.CLIDI.LT.0.3)                  GO  TO  12
      IF(TIJ.LT.PME(7,I).AND.IVALON.EQ.0)                   GO  TO  12
      IF(TIJ.LT.PME(7,J).AND.CLIDJ.LT.0.3)                  GO  TO  12
      IF(TIJ.LT.PME(7,J).AND.IVALON.EQ.0)                   GO  TO  12
C   !!!!!
      CALL  B2IJ(I,J,B2)
      CALL  SGIJ(I,J,SIG)
      SIGV=SIG
      IF(IVALON.NE.0)    SIGV=CLIDI*CLIDJ*SIG
      IF(B2.GT.(SIGV/31.41592))                             GO  TO  12
      IF(IJPA.LT.50000)                                     GO  TO  11
c-Sob      write(16,100)
c-Sob      write( *,100)
c-Sob  100 FORMAT(1X,'IJPA>50000')
      GO  TO  12
   11 IJPA=IJPA+1
      TIJ4(IJPA)=TIJ
      INT4(IJPA)=10000*I+J
CC
C     IF(IME(4,I).EQ.0.AND.IME(4,J).EQ.0)
C    *CALL NAMIJ(IP,IPO,INT4(IJPA),I,J,IJPA,TIJ)
CC
      NYY(J)=NYY(J)+1
      NYY(I)=NYY(I)+1
   12 CONTINUE
   13 CONTINUE
      DO  14  I=1,MV
      IF(NYY(I).EQ.MYY(I).AND.MYY(I).LE.1)    NYY(I)=0
      IF(NYY(I).EQ.1.AND.MYY(I).EQ.0)         NYY(I)=2
      MYY(I)=NYY(I)
   14 CONTINUE
   15 CONTINUE
      IF(IJPA.LE.0)                                          RETURN
      IJM=1
      IJPA0=IJPA
      DO  19  IJP=1,IJPA0
   16 IF(IJP.GT.IJPA)                                      GO  TO  20
      TIJ=TIJ4(IJP)
      IF(TIJ.LT.0.01)                                      GO  TO  17
      IJ=INT4(IJP)
      I=(IJ+1)/10000
      J=IJ-10000*I
      IF(I.EQ.J)                                           GO  TO  17
c-Sob      IF(I.LT.1.OR.I.GT.5999.OR.J.LT.1.OR.J.GT.5999)
c-Sob     *write(16,101) IJ,I,J
c-Sob  101 FORMAT(' IJ,I,J=',3I10)
      IF(IPO.EQ.1.AND.NRST.EQ.5.AND.(I+J).EQ.(IST+JST).
     *                          AND.(I*J).EQ.(IST*JST))    GO  TO  17
      IF((TINT+TIJ).GT.TLIMIT)                             GO  TO  17
C   !!!!!
      IF(TIJ.LT.PME(7,I).AND.CLIDER(I).LT.0.3)             GO  TO  17
      IF(TIJ.LT.PME(7,I).AND.IVALON.EQ.0)                  GO  TO  17
      IF(TIJ.LT.PME(7,J).AND.CLIDER(J).LT.0.3)             GO  TO  17
      IF(TIJ.LT.PME(7,J).AND.IVALON.EQ.0)                  GO  TO  17
C   !!!!!
      IF(TAU.LT.0.)                                        GO  TO  18
      IF(TIJ.GT.TAU)                                       GO  TO  19
                                             GO  TO  18
   17 TIJ4(IJP)=TIJ4(IJPA)
      INT4(IJP)=INT4(IJPA)
      IJPA=IJPA-1
      IF(IJPA.LE.0)                                        GO  TO  21
                                             GO  TO  16
   18 TAU=TIJ
      N1=I
      N2=J
      IJM=IJP
   19 CONTINUE
   20 CONTINUE
C      IJ=INT4(IJM)
C      INT4(IJM)=INT4(IJPA)
C      TIJ4(IJM)=TIJ4(IJPA)
C      IJPA=IJPA-1
   21 IF(TAU.LT.0.)                                          RETURN
CC
C     IF(IME(4,N1).EQ.0.AND.IME(4,N2).EQ.0)
C    *CALL NAMIJ(IP,IPO,IJ,N1,N2,IJPA,TAU)
CC
      IF(INTCC.GE.2)                                         RETURN
      IF((IME(4,N1)+IME(4,N2)).EQ.2.OR.IME(4,N2).EQ.1)       RETURN
      M=N2
      N2=N1
      N1=M
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE ERAIJ(N)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/INTCEN/IPO,INT1(240),INT2(5999),INT3(5999),
     *              INT4(50000),IJPA,IST,JST,NRST
      COMMON/TAUIJ/ TAUK,TPTS,TYPS,TYTS,TYYS,TIJ4(50000)
      IF(IJPA.LE.0)          RETURN
      IJPA0=IJPA
      DO  12  IJP=1,IJPA0
   10 IF(IJP.GT.IJPA)        RETURN
      IJ=INT4(IJP)
      I=(IJ+1)/10000
      J=IJ-10000*I
      IF(N.EQ.I.OR.N.EQ.J)  GO  TO  11
      GO  TO  12
   11 TIJ4(IJP)=TIJ4(IJPA)
      INT4(IJP)=INT4(IJPA)
      IJPA=IJPA-1
      GO  TO  10
   12 CONTINUE
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE REPIJ(N,M)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/INTCEN/IPO,INT1(240),INT2(5999),INT3(5999),
     *              INT4(50000),IJPA,IST,JST,NRST
      COMMON/TAUIJ/ TAUK,TPTS,TYPS,TYTS,TYYS,TIJ4(50000)
      IF(IJPA.LE.0)          RETURN
      DO  11  IJP=1,IJPA
      IJ=INT4(IJP)
      I=(IJ+1)/10000
      J=IJ-10000*I
      IF(M.EQ.I.OR.M.EQ.J)  GO  TO  10
      GO  TO  11
   10 IF(M.EQ.I)  I=N
      IF(M.EQ.J)  J=N
      IJ=10000*I+J
      INT4(IJP)=IJ
   11 CONTINUE
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE NAMIJ(IP,IPO,IJ,I,J,IJPA,TAU)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL *8 TAU
      CHARACTER*8 PNI,PNJ
      COMMON /MEMAGT/ PME(9,5999),IME(5,5999)
      COMMON /IDPME/ IDPME(5999)
      DIMENSION PI(9),IPI(5),PJ(9),IPJ(5)
      DO  10  K=1,9
      PI(K)=PME(K,I)
      PJ(K)=PME(K,J)
   10 CONTINUE
      DO  11  K=1,5
      IPI(K)=IME(K,I)
      IPJ(K)=IME(K,J)
   11 CONTINUE
      IDI=IDPME(I)
      CALL PANUID(IDI,IKI,PNI)
      IDJ=IDPME(J)
      CALL PANUID(IDJ,IKJ,PNJ)
C     IF(IPI(4).EQ.0.AND.IPJ(4).EQ.0)
C    *write(16,100) IP,IPO,PNI,PNJ,IJ,I,J,IJPA,TAU
c-Sob      IF(IPI(4).EQ.0.AND.IPJ(4).EQ.0)
c-Sob     *write( *,100) IP,IPO,PNI,PNJ,IJ,I,J,IJPA,TAU
c-Sob  100 FORMAT(11X,2I3,1X,2A8,I7,2I5,I6,F8.3)
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE  B2IJ(I,J,B2)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL *8 B2
      COMMON /MEMAGT/ PME(9,5999),IME(5,5999)
      DIMENSION R(3),P(3),RS(3),PS(3)
      EI=PME(8,I)+PME(9,I)
      EJ=PME(8,J)+PME(9,J)
      E=EI+EJ
      PR=0.
      PPI=0.
      PPJ=0.
      DO  10  K=1,3
      R(K)=PME(K,I)-PME(K,J)
      P(K)=PME(K+3,I)+PME(K+3,J)
      PR=PR+P(K)*R(K)
      PPI=PPI+P(K)*PME(K+3,I)
      PPJ=PPJ+P(K)*PME(K+3,J)
   10 CONTINUE
      U=DSQRT(E**2-P(1)**2-P(2)**2-P(3)**2)
      EIS=(U**2+PME(9,I)**2-PME(9,J)**2)/(2.*U)
      EJS=U-EIS
      PRS=0.
      DO  11  K=1,3
      PS(K)=(PME(K+3,I)+P(K)/U*(PPI/(E+U)-EI)-
     -       PME(K+3,J)-P(K)/U*(PPJ/(E+U)-EJ))/2.
      RS(K)=R(K)+P(K)*PR/(E+U)/U
      PRS=PRS+PS(K)*RS(K)
   11 CONTINUE
      PS2=PS(1)**2+PS(2)**2+PS(3)**2
      RS2=RS(1)**2+RS(2)**2+RS(3)**2
      B2=RS2-(PRS**2)/PS2
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      FUNCTION  BIM2(P1,P2)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL *8 BIM2
      DIMENSION R(3),P(3),RS(3),PS(3),P1(9),P2(9)
      EI=P1(8)+P1(9)
      EJ=P2(8)+P2(9)
      E=EI+EJ
      PR=0.
      PPI=0.
      PPJ=0.
      DO  10  K=1,3
      R(K)=P1(K)-P2(K)
      P(K)=P1(K+3)+P2(K+3)
      PR=PR+P(K)*R(K)
      PPI=PPI+P(K)*P1(K+3)
      PPJ=PPJ+P(K)*P2(K+3)
   10 CONTINUE
      U=DSQRT(E**2-P(1)**2-P(2)**2-P(3)**2)
      EIS=(U**2+P1(9)**2-P2(9)**2)/(2.*U)
      EJS=U-EIS
      PRS=0.
      DO  11  K=1,3
      PS(K)=(P1(K+3)+P(K)/U*(PPI/(E+U)-EI)-
     -       P2(K+3)-P(K)/U*(PPJ/(E+U)-EJ))/2.
      PS(K)=(P1(K+3)*(EJ+EJS)-P2(K+3)*(EI+EIS))/(E+U)
      RS(K)=R(K)+P(K)*PR/(E+U)/U
      PRS=PRS+PS(K)*RS(K)
   11 CONTINUE
      PS2=PS(1)**2+PS(2)**2+PS(3)**2
      RS2=RS(1)**2+RS(2)**2+RS(3)**2
      BIM2=RS2-(PRS**2)/PS2
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE  SGIJ(I,J,SG)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 SG
      CHARACTER*8 PNA1,PNA2
      COMMON/MEMAGT/PME(9,5999),IME(5,5999)
      COMMON /IDPME/ IDPME(5999)
      DIMENSION P1(9),P2(9),IP1(5),IP2(5)
      DO  10  K=1,9
      P1(K)=PME(K,I)
      P2(K)=PME(K,J)
      IF(K.LE.5)  IP1(K)=IME(K,I)
      IF(K.LE.5)  IP2(K)=IME(K,J)
   10 CONTINUE
      ID1=IDPME(I)
      CALL PANUID(ID1,IK1,PNA1)
      ID2=IDPME(J)
      CALL PANUID(ID2,IK2,PNA2)
      PX1=P1(4)
      PY1=P1(5)
      PZ1=P1(6)
      AM1=P1(9)
      PX2=P2(4)
      PY2=P2(5)
      PZ2=P2(6)
      AM2=P2(9)
      CALL CROSEC(1,ID1,ID2,PX1,PY1,PZ1,AM1,PX2,PY2,PZ2,AM2,SITO,0)
      SG=SITO
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE PAULI1(P1,P2,IP1,IP2,N1,N2,V,NP,MV,TINT,IP)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL *8 P1,P2,V,TINT
C
      REAL*8 MASN
      COMMON/NCASCA/NCAS,NCPRI
      COMMON/TPROD/TPROD(5999)
      COMMON/PORIG/IORI(3,5999)
      COMMON/IDN12/ID1,ID2
      COMMON/NUCOLL/ NUCOLL,MVCOLL(5999)
      COMMON/CSLID/CLIDER(5999)
      COMMON/RESULT/AN1,AN2,ZN1,ZN2,ENEXT1,ENEXT2,PNUCL1(3),PNUCL2(3),
     *AMNUC1(3),AMNUC2(3)
      COMMON/HCASC/ANUCL1,ANUCL2,ZNUCL1,ZNUCL2,T0,EPS1,EPS2,
     *VPI,A1,A2,C1,C2,D1,D2,R0N1,R0N2,TF01,TF02,RM1,RM2
     */CENTER/XC(2,240),YC(2,240),ZC(2,240),IZ(2,240)
       COMMON/MEMAGT/PMEMO(9,5999),IMEMO(5,5999)
      COMMON/ACTIV/MPA(240),MYP(5999),MYT(5999),MYY(5999)
      COMMON /EXCIT/ TEX1,TEX2,HEX1,HEX2,PEX1,PEX2
      COMMON /HOLPT/ PH1(3),PH2(3),RH1(3),RH2(3),EHOL1,EHOL2,UP,UT
      COMMON/NUCSP/VPR(3),VTA(3),RADP(3),RADT(3),VEV(3),VRE(3),GEV,GRE
     *,VEP(3),VET(3),GEP,GET
      COMMON /RINT/ RINT
      COMMON /IDPME/ IDPME(5999)
      COMMON /SORI/ SORI(5999),SSOR
      DIMENSION P1(9),P2(9),IP1(5),IP2(5),PS(3),V(3),PN1(3),
     *PN2(3),P(3),PL(3),BPR(3),BTA(3),RL(3)
      MASN=0.940
      DO  9  L=1,3
      BPR(L)=-VPR(L)
    9 BTA(L)=-VTA(L)
      TFR1=UP-EPS1
      TFR2=UT-EPS2
      IF(NP-2)10,10,13
   10 DO 11 L=1,9
   11 PMEMO(L,MV+2)=PMEMO(L,MV+3)
      DO 12 L=1,5
   12 IMEMO(L,MV+2)=IMEMO(L,MV+3)
      CLIDER(MV+2)=CLIDER(MV+3)
      IDPME(MV+2)=IDPME(MV+3)
      GO TO 16
   13 DO 14 L=1,9
      TEMP=PMEMO(L,MV+2)
      PMEMO(L,MV+2)=PMEMO(L,MV+3)
   14 PMEMO(L,MV+3)=TEMP
      DO 15 L=1,5
      ITEMP=IMEMO(L,MV+2)
      IMEMO(L,MV+2)=IMEMO(L,MV+3)
   15 IMEMO(L,MV+3)=ITEMP
      TEMC=CLIDER(MV+2)
      CLIDER(MV+2)=CLIDER(MV+3)
      CLIDER(MV+3)=TEMC
      IDTE=IDPME(MV+2)
      IDPME(MV+2)=IDPME(MV+3)
      IDPME(MV+3)=IDTE
   16 DO 18 L=1,2
      IF(IMEMO(3,MV+L).NE.0.OR.IMEMO(5,MV+L).NE.0.OR.IMEMO(2,MV+L).NE.0)
     *GO  TO  18
      PS(1)=PMEMO(4,MV+L)
      PS(2)=PMEMO(5,MV+L)
      PS(3)=PMEMO(6,MV+L)
      CALL KINEMA(PS,V,PL,CTL,STL,CFL,SFL,TL,MASN)
      CALL KINEMA(PL,BPR,PN1,CT1,ST1,CF1,SF1,T1,MASN)
      IFORB=1
      IF(T1.LE.TFR1)  GO  TO  19
      CALL KINEMA(PL,BTA,PN2,CT2,ST2,CF2,SF2,T2,MASN)
      IFORB=2
      IF(T2.LE.TFR2)  GO  TO  19
   18 CONTINUE
      GO  TO  20
   19 IP=0
c-Sob      IF(NCAS.GE.NCPRI) THEN
c-Sob        IF(IFORB.EQ.1) write(16,302) T1,TFR1
c-Sob        IF(IFORB.EQ.2) write(16,303) T2,TFR2
c-Sob  302 FORMAT(1X,'PAULI1 FORBIDED: T1=',F6.3,' < TFR1=',F6.3)
c-Sob  303 FORMAT(1X,'PAULI1 FORBIDED: T2=',F6.3,' < TFR2=',F6.3)
c-Sob      ENDIF
           RETURN
   20 CONTINUE
      NUCOLL=NUCOLL+1
      NA1=AN1+0.1
      NA2=AN2+0.1
   24 E1=IP1(1)
      E2=IP2(1)
      J2=AN2+0.1
      XC(2,N2)=XC(2,J2)
      YC(2,N2)=YC(2,J2)
      ZC(2,N2)=ZC(2,J2)
      IZ(2,N2)=IZ(2,J2)
      ZN2=ZN2-E2
      AN2=AN2-1.
      J1=AN1+0.1
      XC(1,N1)=XC(1,J1)
      YC(1,N1)=YC(1,J1)
      ZC(1,N1)=ZC(1,J1)
      IZ(1,N1)=IZ(1,J1)
      ZN1=ZN1-E1
      AN1=AN1-1.
      MPA(N1)=MPA(J1)
      CALL  RECUL(1,PH1(1),PH1(2),PH1(3),RH1(1),RH1(2),RH1(3))
      CALL  RECUL(2,PH2(1),PH2(2),PH2(3),RH2(1),RH2(2),RH2(3))
      ENEXT1=ENEXT1+TFR1-EHOL1
      ENEXT2=ENEXT2+TFR2-EHOL2
      L=1
      NP1=NP
   25 M=MV+L
      EL=IMEMO(1,M)
      QL=IMEMO(4,M)
      P(1)=PMEMO(4,M)
      P(2)=PMEMO(5,M)
      P(3)=PMEMO(6,M)
      CALL  KINEMA(P,V,PL,CTL,STL,CFL,SFL,TL,PMEMO(9,M))
      PMEMO(4,M)=PL(1)
      PMEMO(5,M)=PL(2)
      PMEMO(6,M)=PL(3)
      PMEMO(8,M)=TL
      CALL  KINEMR(V,M,RL,TAUL)
      CALL  TMATUR(NP,V,M,P1,P2,MV,P,TAUL,TMAT)
      PMEMO(7,M)=TMAT
      TPROD(M)=TINT
      MVCOLL(M)=NUCOLL
      IORI(1,M)=ID1
      IORI(2,M)=ID2
      IORI(3,M)=0
      SORI(M)=SSOR
      IF(IMEMO(5,M).EQ.0)    GO  TO  125
      TAUM0= DFLOAT(IMEMO(5,M))/1000.
      TAUML=TAUM0*(1.+PMEMO(8,M)/PMEMO(9,M))
      IMEMO(5,M)=IDINTG((TAUML+PMEMO(7,M))*1000.)
      IF(IMEMO(5,M).EQ.0)  IMEMO(5,M)=1
  125 CONTINUE
      ELM=PMEMO(8,M)+PMEMO(9,M)
C      PMEMO(1,M)=(P1(1)+P2(1))/2.+RL(1)-PMEMO(4,M)/ELM*TAUL
C      PMEMO(2,M)=(P1(2)+P2(2))/2.+RL(2)-PMEMO(5,M)/ELM*TAUL
C      PMEMO(3,M)=(P1(3)+P2(3))/2.+RL(3)-PMEMO(6,M)/ELM*TAUL
      DR=RINT*RNDMD(-1)**(1./3.)
      FI=6.283185*RNDMD(-1)
      CT=1.-2.*RNDMD(-1)
      ST=DSQRT(1.-CT**2)
      DX=DR*ST*DCOS(FI)
      DY=DR*ST*DSIN(FI)
      DZ=DR*CT
      PMEMO(1,M)=(P1(1)+P2(1))/2.+DX
      PMEMO(2,M)=(P1(2)+P2(2))/2.+DY
      PMEMO(3,M)=(P1(3)+P2(3))/2.+DZ
      MYP(M)=1
      MYT(M)=1
      MYY(M)=1
      IF(IMEMO(2,M).NE.0)  MYP(M)=0
      IF(IMEMO(2,M).NE.0)  MYT(M)=0
      IF(IMEMO(2,M).NE.0)  MYY(M)=0
c-Sob      IF(NCAS.GE.NCPRI)  write(16,301)
c-Sob     * M,(PMEMO(K,M),K=1,9),(IMEMO(K,M),K=1,5),CLIDER(M),IDPME(M)
c-Sob  301 FORMAT(1X,'NEWP',I5,9(1X,F 8.3),2X,4I2,I15,F6.3,I5)
      IF(L-NP1)  36,37,37
   36 L=L+1
            GO TO 25
   37 IP=1
      MV=MV+NP1
      TEX1=TEX1+1.
      HEX1=HEX1+1.
      TEX2=TEX2+1.
      HEX2=HEX2+1.
c-Sob      IF(AN1.LT.ZN1.OR.ZN1.LT.0..OR.AN2.LT.ZN2.OR.ZN2.LT.0.)
c-Sob     *write(16,300) AN1,ZN1,AN2,ZN2
c-Sob  300 FORMAT(2X,'PAULI1',4(2X,F5.0))
      RETURN
             END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE PAULI2(P1,P2,P3,IP1,IP2,IP3,N1,N2,N3,V,NP,MV,TINT,
     *IP,OBR1)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 P1,P2,P3,V,TINT,OBR1
      REAL*8 MASN
      COMMON/NCASCA/NCAS,NCPRI
      COMMON/TPROD/TPROD(5999)
      COMMON/PORIG/IORI(3,5999)
      COMMON/IDN12/ID1,ID2
      COMMON/CSLID/CLIDER(5999)
      COMMON/NUCOLL/ NUCOLL,MVCOLL(5999)
      COMMON/RESULT/AN1,AN2,ZN1,ZN2,ENEXT1,ENEXT2,PNUCL1(3),PNUCL2(3),
     *AMNUC1(3),AMNUC2(3)
      COMMON/HCASC/ANUCL1,ANUCL2,ZNUCL1,ZNUCL2,T0,EPS1,EPS2,
     *VPI,A1,A2,C1,C2,D1,D2,R0N1,R0N2,TF01,TF02,RM1,RM2
     */CENTER/XC(2,240),YC(2,240),ZC(2,240),IZ(2,240)
       COMMON/MEMAGT/PMEMO(9,5999),IMEMO(5,5999)
      COMMON/ACTIV/MPA(240),MYP(5999),MYT(5999),MYY(5999)
      COMMON /EXCIT/ TEX1,TEX2,HEX1,HEX2,PEX1,PEX2
      COMMON/NUCSP/VPR(3),VTA(3),RADP(3),RADT(3),VEV(3),VRE(3),GEV,GRE
     *,VEP(3),VET(3),GEP,GET
      COMMON /RINT/ RINT
      COMMON /HOLPT/ PH1(3),PH2(3),RH1(3),RH2(3),EHOL1,EHOL2,UP,UT
      COMMON /IDPME/ IDPME(5999)
      COMMON /SORI/ SORI(5999),SSOR
      COMMON /BARPOT/ POT
      DIMENSION B(3),P1(9),P2(9),IP1(5),IP2(5),PS(3),V(3),PN1(3),
     *PN2(3),P(3),PL(3),P3(9),IP3(5),RL(3)
      MASN=0.940
      B(1)=VPR(1)
      B(2)=VPR(2)
      B(3)=VPR(3)
      GPR=1./DSQRT(1.-VPR(1)**2-VPR(2)**2-VPR(3)**2)
      R1=DSQRT(P2(1)**2+P2(2)**2+P2(3)**2)
      IABS=0
      TFR1=UP-EPS1
      IF(NP-2)32,10,13
   10 DO 11 L=1,9
   11 PMEMO(L,MV+2)=PMEMO(L,MV+3)
      DO 12 L=1,5
   12 IMEMO(L,MV+2)=IMEMO(L,MV+3)
      IF(IP1(4).EQ.0.AND.(IMEMO(4,MV+1)+IMEMO(4,MV+2)).EQ.2)  IABS=1
      CLIDER(MV+2)=CLIDER(MV+3)
      IDPME(MV+2)=IDPME(MV+3)
      GO TO 16
   13 DO 14 L=1,9
      TEMP=PMEMO(L,MV+2)
      PMEMO(L,MV+2)=PMEMO(L,MV+3)
   14 PMEMO(L,MV+3)=TEMP
      DO 15 L=1,5
      ITEMP=IMEMO(L,MV+2)
      IMEMO(L,MV+2)=IMEMO(L,MV+3)
   15 IMEMO(L,MV+3)=ITEMP
      TEMC=CLIDER(MV+2)
      CLIDER(MV+2)=CLIDER(MV+3)
      CLIDER(MV+3)=TEMC
      IDTE=IDPME(MV+2)
      IDPME(MV+2)=IDPME(MV+3)
      IDPME(MV+3)=IDTE
   16 IF(IMEMO(4,MV+1).NE.1)     GO  TO  17
      IF(IMEMO(3,MV+1).NE.0.OR.IMEMO(5,MV+1).NE.0.OR.IMEMO(2,MV+1).NE.0)
     *GO  TO  17
      PS(1)=PMEMO(4,MV+1)
      PS(2)=PMEMO(5,MV+1)
      PS(3)=PMEMO(6,MV+1)
      CALL CINEMA(PS,V,PN1,CT1,ST1,CF1,SF1,T1,MASN)
      TN=T1
      IF(T1-TFR1)   31,31,17
   17 IF(IMEMO(4,MV+2).NE.1)  GO TO 32
      IF(IMEMO(3,MV+2).NE.0.OR.IMEMO(5,MV+2).NE.0.OR.IMEMO(2,MV+2).NE.0)
     *GO  TO  32
      PS(1)=PMEMO(4,MV+2)
      PS(2)=PMEMO(5,MV+2)
      PS(3)=PMEMO(6,MV+2)
      CALL CINEMA(PS,V,PN2,CT2,ST2,CF2,SF2,T2,MASN)
      TN=T2
      IF(T2-TFR1) 31,31,32
   31 IP=0
c-Sob      IF(NCAS.GE.NCPRI) write(16,302) TN,TFR1
c-Sob  302 FORMAT(1X,'PAULI2 FORBIDED: TN=',F6.3,' < TFR1=',F6.3)
      RETURN
   32 CONTINUE
      NUCOLL=NUCOLL+1
      IBAR=0
      if(IP1(4).EQ.-1.AND.IP1(3).EQ.0.AND.IP1(5).EQ.0) IBAR=1
      E2=IP2(1)
      J1=AN1+0.1
      IF(IABS.EQ.1.AND.J1.EQ.N3)   J1=J1-1
      XC(1,N2)=XC(1,J1)
      YC(1,N2)=YC(1,J1)
      ZC(1,N2)=ZC(1,J1)
      IZ(1,N2)=IZ(1,J1)
      MPA(N2)=MPA(J1)
      ZN1=ZN1-E2
      AN1=AN1-1.
      CALL  RECUL(1,PH1(1),PH1(2),PH1(3),RH1(1),RH1(2),RH1(3))
      ENEXT1=ENEXT1+TFR1-EHOL1
      IF(NP.EQ.1)  GO  TO  51
C **************** KOSTYA 08.04.92 ****************
C      IF(IP1(4))  49,49,51
C   49 IF(IMEMO(4,MV+1)+IMEMO(4,MV+2)-2) 51,50,51
       IF(IABS.NE.1)  GO  TO  51
C *************************************************
   50 ENEXT1=ENEXT1+POTEN(P3,IP3,A1,C1,D1,TF01,VPI,EPS1)-EPS1-P3(8)
      TEX1=TEX1+1.
      HEX1=HEX1+1.
      E3=IP3(1)
      CALL RECUL(1,P3(4),P3(5),P3(6),P3(1),P3(2),P3(3))
      J1=AN1+0.1
      IF(J1.LT.N3)   J1=N3
      XC(1,N3)=XC(1,J1)
      YC(1,N3)=YC(1,J1)
      ZC(1,N3)=ZC(1,J1)
      IZ(1,N3)=IZ(1,J1)
      MPA(N3)=MPA(J1)
      ZN1=ZN1-E3
      AN1=AN1-1.
   51 CONTINUE
      L=1
      NP1=NP
      VPI1=VPI
      NABSN=0
   33 M=MV+L
      EL=IMEMO(1,M)
      QL=IMEMO(4,M)
      P(1)=PMEMO(4,M)
      P(2)=PMEMO(5,M)
      P(3)=PMEMO(6,M)
      CALL CINEMA(P,V,PL,CTL,STL,CFL,SFL,TL,PMEMO(9,M))
      TL0=TL
      IF(IMEMO(4,M).NE.0)   GO  TO  105
      IF(IMEMO(1,M))  103,105,104
  103 IF((ZN1*(AN1-1.)/2.).LT.1.1)   GO  TO  46
      GO  TO  105
  104 IF(((AN1-1.)*(AN1-ZN1)/2.).LT.1.1)   GO  TO  46
  105 CONTINUE
      IF(IMEMO(5,M).NE.0.OR.IMEMO(3,M).NE.0)   GO  TO  46
      IF(IMEMO(4,M).LT.0.OR.IMEMO(2,M).NE.0)   GO  TO  46
      TL0=TL-QL*(TFR1+EPS1)-(1.-QL)*VPI
   35 CUT=EPS1*QL+OBR1*EL*ZN1/ZNUCL1
      IF(EL)  101,102,102
  101 CUT=0.
  102 CONTINUE
      IF(TL0)  36,100,100
  100 CONTINUE
      IF(TL0-CUT) 36,36,46
   36 ENEXT1=ENEXT1+TL0+EPS1*QL+(.14+VPI)*(1.-QL)
      CALL  RECUL(1,PL(1),PL(2),PL(3),P1(1),P1(2),P1(3))
      NABSN=NABSN+IMEMO(4,M)
      AN1=AN1+QL
      ZN1=ZN1+EL
      IF(IMEMO(4,M)-1) 141,37,141
   37 IF(NABSN-1) 38,38,39
   38 J1=AN1+.1
      XC(1,J1)=P2(1)
      YC(1,J1)=P2(2)
      ZC(1,J1)=P2(3)
      GO TO 40
  141 J1=AN1+0.1
      DO  108  K=1,J1
      IF(IMEMO(1,M))  106,41,107
  106 IF(IZ(1,K).NE.1)   GO  TO  108
      IZ(1,K)=0
                GO  TO  41
  107 IF(IZ(1,K).NE.0)   GO  TO  108
      IZ(1,K)=1
                GO  TO  41
  108 CONTINUE
      GO  TO  41
   39 J1=AN1+.1
      XC(1,J1)=P1(1)
      YC(1,J1)=P1(2)
      ZC(1,J1)=P1(3)
   40 IZ(1,J1)=IMEMO(1,M)
      MPA(J1)=1
      TEX1=TEX1+1.
      PEX1=PEX1+EL
   41 IF(L-NP1) 42,45,45
   42 DO 43 K=1,9
   43 PMEMO(K,M)=PMEMO(K,MV+NP1)
      DO 44 K=1,5
   44 IMEMO(K,M)=IMEMO(K,MV+NP1)
      CLIDER(M)=CLIDER(MV+NP1)
      IDPME(M)=IDPME(MV+NP1)
      NP1=NP1-1
              GO TO 33
   45 NP1=NP1-1
              GO TO 48
   46 continue
      if(IMEMO(4,M).EQ.-1.AND.IMEMO(3,M).EQ.0.AND.IMEMO(5,M).EQ.0.
     &AND.IBAR.EQ.1) then
      TP0=TL-POT
      if(TL0.LE.0.)  then
        ENEXT1=ENEXT1+TL+0.940
        AN1=AN1-1.
        ZN1=ZN1+EL
        call RECUL(1,PL(1),PL(2),PL(3),P1(1),P1(2),P1(3))
        go to 41
      endif
      endif
      PM=DSQRT(TL0*(TL0+2.*PMEMO(9,M)))
      PN2(1)=PM*STL*CFL
      PN2(2)=PM*STL*SFL
      PN2(3)=PM*CTL
      CALL CINEMA(PN2,B,PS,CT,ST,CF,SF,TTL,PMEMO(9,M))
      PMEMO(4,M)=PS(1)
      PMEMO(5,M)=PS(2)
      PMEMO(6,M)=PS(3)
      PMEMO(8,M)=TTL
      CALL  KINEMR(B,M,RL,TAUL)
      CALL  TMATUR(NP,V,M,P1,P2,MV,P,TAUL,TMAT)
      PMEMO(7,M)=TMAT
      TPROD(M)=TINT
      MVCOLL(M)=NUCOLL
      IF(IDPME(N1).EQ.IDPME(M).AND.NP.EQ.2) THEN
* rescattering (elastic or charge exchange)  of the N1 particle
      IORI(1,M)=IORI(1,N1)
      IORI(2,M)=IORI(2,N1)
      IORI(3,M)=IORI(3,N1)+1
      SORI(M)=SORI(N1)
      ELSE
* production of a new particle M
      IORI(1,M)=ID1
      IORI(2,M)=ID2
      IORI(3,M)=0
      SORI(M)=SSOR
      ENDIF
      MYP(M)=1
      MYT(M)=1
      MYY(M)=1
      IF(IMEMO(2,M).NE.0)  MYP(M)=0
      IF(IMEMO(2,M).NE.0)  MYT(M)=0
      IF(IMEMO(2,M).NE.0)  MYY(M)=0
      IF(IMEMO(5,M).EQ.0)  GO  TO  146
      TAUM0= DFLOAT(IMEMO(5,M))/1000.
      TAUML=TAUM0*(1.+PMEMO(8,M)/PMEMO(9,M))
      IMEMO(5,M)=IDINTG((TAUML+PMEMO(7,M))*1000.)
      IF(IMEMO(5,M).EQ.0)  IMEMO(5,M)=1
  146 CONTINUE
      X=(P1(1)+P2(1))/2.
      Y=(P1(2)+P2(2))/2.
      Z=(P1(3)+P2(3))/2.
      VRK=X*VPR(1)+Y*VPR(2)+Z*VPR(3)
      DR=RINT*RNDMD(-1)**(1./3.)
      FI=6.283185*RNDMD(-1)
      CT=1.-2.*RNDMD(-1)
      ST=DSQRT(1.-CT**2)
      DX=DR*ST*DCOS(FI)
      DY=DR*ST*DSIN(FI)
      DZ=DR*CT
      PMEMO(1,M)=RADP(1)+X-VPR(1)*VRK*GPR/(GPR+1.)+DX
      PMEMO(2,M)=RADP(2)+Y-VPR(2)*VRK*GPR/(GPR+1.)+DY
      PMEMO(3,M)=RADP(3)+Z-VPR(3)*VRK*GPR/(GPR+1.)+DZ
c-Sob      IF(NCAS.GE.NCPRI)  write(16,301)
c-Sob     * M,(PMEMO(K,M),K=1,9),(IMEMO(K,M),K=1,5),CLIDER(M),IDPME(M)
c-Sob  301 FORMAT(1X,'NEWP',I5,9(1X,F 8.3),2X,4I2,I15,F6.3,I5)
      IF(L-NP1) 47,48,48
   47 L=L+1
            GO TO 33
   48 IP=1
             MV=MV+NP1-1
   52 DO 53 K=1,9
   53 PMEMO(K,N1)=PMEMO(K,MV+1)
      DO 54 K=1,5
   54 IMEMO(K,N1)=IMEMO(K,MV+1)
      MYP(N1)=MYP(MV+1)
      MYT(N1)=MYT(MV+1)
      MYY(N1)=MYY(MV+1)
      TPROD(N1)=TPROD(MV+1)
      MVCOLL(N1)=MVCOLL(MV+1)
      IORI(1,N1)=IORI(1,MV+1)
      IORI(2,N1)=IORI(2,MV+1)
      IORI(3,N1)=IORI(3,MV+1)
      SORI(N1)=SORI(MV+1)
      CLIDER(N1)=CLIDER(MV+1)
      IDPME(N1)=IDPME(MV+1)
      CALL ERAIJ(N1)
      IF(NP1.EQ.0)              CALL REPIJ(N1,MV+1)
      TEX1=TEX1+1.
      HEX1=HEX1+1.
c-Sob      IF(AN1.LT.ZN1.OR.ZN1.LT.0.)  write(16,300) AN1,ZN1
c-Sob  300 FORMAT(2X,'PAULI2',2(2X,F5.0))
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE PAULI3(P1,P2,P3,IP1,IP2,IP3,N1,N2,N3,V,NP,MV,TINT,
     *IP,OBR2)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 P1,P2,P3,V,TINT,OBR2
      REAL*8 MASN
      COMMON/NCASCA/NCAS,NCPRI
      COMMON/TPROD/TPROD(5999)
      COMMON/PORIG/IORI(3,5999)
      COMMON/IDN12/ID1,ID2
      COMMON/CSLID/CLIDER(5999)
      COMMON/NUCOLL/ NUCOLL,MVCOLL(5999)
      COMMON/RESULT/AN1,AN2,ZN1,ZN2,ENEXT1,ENEXT2,PNUCL1(3),PNUCL2(3),
     *AMNUC1(3),AMNUC2(3)
      COMMON/HCASC/ANUCL1,ANUCL2,ZNUCL1,ZNUCL2,T0,EPS1,EPS2,
     *VPI,A1,A2,C1,C2,D1,D2,R0N1,R0N2,TF01,TF02,RM1,RM2
     */CENTER/XC(2,240),YC(2,240),ZC(2,240),IZ(2,240)
       COMMON/MEMAGT/PMEMO(9,5999),IMEMO(5,5999)
      COMMON/ACTIV/MPA(240),MYP(5999),MYT(5999),MYY(5999)
      COMMON /EXCIT/ TEX1,TEX2,HEX1,HEX2,PEX1,PEX2
      COMMON/NUCSP/VPR(3),VTA(3),RADP(3),RADT(3),VEV(3),VRE(3),GEV,GRE
     *,VEP(3),VET(3),GEP,GET
      COMMON /RINT/ RINT
      COMMON /HOLPT/ PH1(3),PH2(3),RH1(3),RH2(3),EHOL1,EHOL2,UP,UT
      COMMON /IDPME/ IDPME(5999)
      COMMON /SORI/ SORI(5999),SSOR
      COMMON /BARPOT/ POT
      DIMENSION B(3),P1(9),P2(9),IP1(5),IP2(5),PS(3),V(3),PN1(3),
     *PN2(3),P(3),PT(3),PL(3),P3(9),IP3(5),RL(3)
      MASN=0.940
      B(1)=VTA(1)
      B(2)=VTA(2)
      B(3)=VTA(3)
      GTA=1./DSQRT(1.-VTA(1)**2-VTA(2)**2-VTA(3)**2)
      R2=DSQRT(P2(1)**2+P2(2)**2+P2(3)**2)
      IABS=0
      TFR2=UT-EPS2
      IF(NP-2)32,10,13
   10 DO 11 L=1,9
   11 PMEMO(L,MV+2)=PMEMO(L,MV+3)
      DO 12 L=1,5
   12 IMEMO(L,MV+2)=IMEMO(L,MV+3)
      IF(IP1(4).EQ.0.AND.(IMEMO(4,MV+1)+IMEMO(4,MV+2)).EQ.2)  IABS=1
      CLIDER(MV+2)=CLIDER(MV+3)
      IDPME(MV+2)=IDPME(MV+3)
      GO TO 16
   13 DO 14 L=1,9
      TEMP=PMEMO(L,MV+2)
      PMEMO(L,MV+2)=PMEMO(L,MV+3)
   14 PMEMO(L,MV+3)=TEMP
      DO 15 L=1,5
      ITEMP=IMEMO(L,MV+2)
      IMEMO(L,MV+2)=IMEMO(L,MV+3)
   15 IMEMO(L,MV+3)=ITEMP
      TEMC=CLIDER(MV+2)
      CLIDER(MV+2)=CLIDER(MV+3)
      CLIDER(MV+3)=TEMC
      IDTE=IDPME(MV+2)
      IDPME(MV+2)=IDPME(MV+3)
      IDPME(MV+3)=IDTE
   16 IF(IMEMO(4,MV+1).NE.1)   GO  TO  17
      IF(IMEMO(3,MV+1).NE.0.OR.IMEMO(5,MV+1).NE.0.OR.IMEMO(2,MV+1).NE.0)
     *GO  TO  17
      PS(1)=PMEMO(4,MV+1)
      PS(2)=PMEMO(5,MV+1)
      PS(3)=PMEMO(6,MV+1)
      CALL CINEMA(PS,V,PN1,CT1,ST1,CF1,SF1,TN,MASN)
      IF(TN-TFR2) 31,31,17
   17 IF(IMEMO(4,MV+2).NE.1)  GO TO 32
      IF(IMEMO(3,MV+2).NE.0.OR.IMEMO(5,MV+2).NE.0.OR.IMEMO(2,MV+2).NE.0)
     *GO  TO  32
      PS(1)=PMEMO(4,MV+2)
      PS(2)=PMEMO(5,MV+2)
      PS(3)=PMEMO(6,MV+2)
      CALL CINEMA(PS,V,PN2,CT2,ST2,CF2,SF2,TN,MASN)
      IF(TN-TFR2) 31,31,32
   31 IP=0
c-Sob      IF(NCAS.GE.NCPRI) write(16,302) TN,TFR2
c-Sob  302 FORMAT(1X,'PAULI3 FORBIDED: TN=',F6.3,' < TFR2=',F6.3)
      RETURN
   32 CONTINUE
      IBAR=0
      if(IP1(4).EQ.-1.AND.IP1(3).EQ.0.AND.IP1(5).EQ.0) IBAR=1
      NUCOLL=NUCOLL+1
      E2=IP2(1)
      J2=AN2+0.1
      IF(IABS.EQ.1.AND.J2.EQ.N3)   J2=J2-1
   25 XC(2,N2)=XC(2,J2)
      YC(2,N2)=YC(2,J2)
      ZC(2,N2)=ZC(2,J2)
      IZ(2,N2)=IZ(2,J2)
      ZN2=ZN2-E2
      AN2=AN2-1.
      CALL  RECUL(2,PH2(1),PH2(2),PH2(3),RH2(1),RH2(2),RH2(3))
      ENEXT2=ENEXT2+TFR2-EHOL2
      IF(NP.EQ.1) GO  TO  51
C **************** KOSTYA 08.04.92 ****************
C      IF(IP1(4))  49,49,51
C   49 IF(IMEMO(4,MV+1)+IMEMO(4,MV+2)-2) 51,50,51
       IF(IABS.NE.1)  GO  TO  51
C *************************************************
   50 ENEXT2=ENEXT2+POTEN(P3,IP3,A2,C2,D2,TF02,VPI,EPS2)-EPS2-P3(8)
      TEX2=TEX2+1.
      HEX2=HEX2+1.
      CALL RECUL(2,P3(4),P3(5),P3(6),P3(1),P3(2),P3(3))
      E3=IP3(1)
      J2=AN2+0.1
      IF(J2.LT.N3)   J2=N3
  132 XC(2,N3)=XC(2,J2)
      YC(2,N3)=YC(2,J2)
      ZC(2,N3)=ZC(2,J2)
      IZ(2,N3)=IZ(2,J2)
      ZN2=ZN2-E3
      AN2=AN2-1.
   51 CONTINUE
      L=1
      NP1=NP
      VPI1=VPI
      NABSN=0
   33 M=MV+L
      EL=IMEMO(1,M)
      QL=IMEMO(4,M)
      P(1)=PMEMO(4,M)
      P(2)=PMEMO(5,M)
      P(3)=PMEMO(6,M)
      CALL CINEMA(P,V,PL,CTL,STL,CFL,SFL,TL,PMEMO(9,M))
      TL0=TL
      IF(IMEMO(4,M).NE.0)   GO  TO  105
      IF(IMEMO(1,M))  103,105,104
  103 IF((ZN2*(AN2-1.)/2.).LT.1.1)   GO  TO  46
      GO  TO  105
  104 IF(((AN2-1.)*(AN2-ZN2)/2.).LT.1.1)   GO  TO  46
  105 CONTINUE
      IF(IMEMO(5,M).NE.0.OR.IMEMO(3,M).NE.0)   GO  TO  46
      IF(IMEMO(4,M).LT.0.OR.IMEMO(2,M).NE.0)   GO  TO  46
      TL0=TL-QL*(TFR2+EPS2)-(1.-QL)*VPI
   35 CUT=EPS2*QL+OBR2*EL*ZN2/ZNUCL2
      IF(EL)  101,102,102
  101 CUT=0.
  102 CONTINUE
      IF(TL0)  36,100,100
  100 CONTINUE
      IF(TL0-CUT) 36,36,46
   36 ENEXT2=ENEXT2+TL0+EPS2*QL+(.14+VPI)*(1.-QL)
      CALL RECUL(2,PL(1),PL(2),PL(3),P1(1),P1(2),P1(3))
      NABSN=NABSN+IMEMO(4,M)
      AN2=AN2+QL
      ZN2=ZN2+EL
      IF(IMEMO(4,M)-1)  141,37,141
   37 IF(NABSN-1) 38,38,39
   38 J2=AN2+.1
      XC(2,J2)=P2(1)
      YC(2,J2)=P2(2)
      ZC(2,J2)=P2(3)
      GO TO 40
  141 J2=AN2+0.1
      DO  108  K=1,J2
      IF(IMEMO(1,M))  106,41,107
  106 IF(IZ(2,K).NE.1)   GO  TO  108
      IZ(2,K)=0
                GO  TO  41
  107 IF(IZ(2,K).NE.0)   GO  TO  108
      IZ(2,K)=1
              GO  TO  41
  108 CONTINUE
      GO  TO  41
   39 J2=AN2+.1
      XC(2,J2)=P1(1)
      YC(2,J2)=P1(2)
      ZC(2,J2)=P1(3)
   40 IZ(2,J2)=IMEMO(1,M)
      TEX2=TEX2+1.
                 PEX2=PEX2+EL
   41 IF(L-NP1) 42,45,45
   42 DO 43 K=1,9
   43 PMEMO(K,M)=PMEMO(K,MV+NP1)
      DO 44 K=1,5
   44 IMEMO(K,M)=IMEMO(K,MV+NP1)
      CLIDER(M)=CLIDER(MV+NP1)
      IDPME(M)=IDPME(MV+NP1)
      NP1=NP1-1
              GO TO 33
   45 NP1=NP1-1
              GO TO 48
   46 CONTINUE
      if(IMEMO(4,M).EQ.-1.AND.IMEMO(3,M).EQ.0.AND.IMEMO(5,M).EQ.0.
     &AND.IBAR.EQ.1) then
      TP0=TL-POT
      if(TL0.LE.0.)  then
        ENEXT2=ENEXT2+TL+0.940
        AN2=AN2-1.
        ZN2=ZN2+EL
        call RECUL(2,PL(1),PL(2),PL(3),P1(1),P1(2),P1(3))
        go to 41
      endif
      endif
      PM0=DSQRT(TL0*(TL0+2.*PMEMO(9,M)))
      DPM=PM0/DSQRT(TL*(TL+2.*PMEMO(9,M)))
      PT(1)=PL(1)*DPM
      PT(2)=PL(2)*DPM
      PT(3)=PL(3)*DPM
      CALL CINEMA(PT,B,PS,CTL,STL,CFL,SFL,TTL,PMEMO(9,M))
      PMEMO(4,M)=PS(1)
      PMEMO(5,M)=PS(2)
      PMEMO(6,M)=PS(3)
      PMEMO(8,M)=TTL
      CALL  KINEMR(B,M,RL,TAUL)
      CALL  TMATUR(NP,V,M,P1,P2,MV,P,TAUL,TMAT)
      PMEMO(7,M)=TMAT
      TPROD(M)=TINT
      MVCOLL(M)=NUCOLL
      IF(IDPME(N1).EQ.IDPME(M).AND.NP.EQ.2) THEN
* rescattering (elastic or charge exchange)  of the N1 particle
       IORI(1,M)=IORI(1,N1)
       IORI(2,M)=IORI(2,N1)
       IORI(3,M)=IORI(3,N1)+1
       SORI(M)=SORI(N1)
      ELSE
* production of a new particle M
       IORI(1,M)=ID1
       IORI(2,M)=ID2
       IORI(3,M)=0
       SORI(M)=SSOR
      ENDIF
      MYP(M)=1
      MYT(M)=1
      MYY(M)=1
      IF(IMEMO(2,M).NE.0)  MYP(M)=0
      IF(IMEMO(2,M).NE.0)  MYT(M)=0
      IF(IMEMO(2,M).NE.0)  MYY(M)=0
      IF(IMEMO(5,M).EQ.0)  GO  TO  146
      TAUM0= DFLOAT(IMEMO(5,M))/1000.
      TAUML=TAUM0*(1.+PMEMO(8,M)/PMEMO(9,M))
      IMEMO(5,M)=IDINTG((TAUML+PMEMO(7,M))*1000.)
      IF(IMEMO(5,M).EQ.0)  IMEMO(5,M)=1
  146 CONTINUE
      X=(P1(1)+P2(1))/2.
      Y=(P1(2)+P2(2))/2.
      Z=(P1(3)+P2(3))/2.
      VRK=X*VTA(1)+Y*VTA(2)+Z*VTA(3)
      DR=RINT*RNDMD(-1)**(1./3.)
      FI=6.283185*RNDMD(-1)
      CT=1.-2.*RNDMD(-1)
      ST=DSQRT(1.-CT**2)
      DX=DR*ST*DCOS(FI)
      DY=DR*ST*DSIN(FI)
      DZ=DR*CT
      PMEMO(1,M)=RADT(1)+X-VTA(1)*VRK*GTA/(GTA+1.)+DX
      PMEMO(2,M)=RADT(2)+Y-VTA(2)*VRK*GTA/(GTA+1.)+DY
      PMEMO(3,M)=RADT(3)+Z-VTA(3)*VRK*GTA/(GTA+1.)+DZ
c-Sob      IF(NCAS.GE.NCPRI)  write(16,301)
c-Sob     * M,(PMEMO(K,M),K=1,9),(IMEMO(K,M),K=1,5),CLIDER(M),IDPME(M)
c-Sob  301 FORMAT(1X,'NEWP',I5,9(1X,F 8.3),2X,4I2,I15,F6.3,I5)
      IF(L-NP1) 47,48,48
   47 L=L+1
            GO TO 33
   48 IP=1
      MV=MV+NP1-1
   52 DO 53 K=1,9
   53 PMEMO(K,N1)=PMEMO(K,MV+1)
      DO 54 K=1,5
   54 IMEMO(K,N1)=IMEMO(K,MV+1)
      MYT(N1)=MYT(MV+1)
      MYY(N1)=MYY(MV+1)
      TPROD(N1)=TPROD(MV+1)
      MVCOLL(N1)=MVCOLL(MV+1)
      IORI(1,N1)=IORI(1,MV+1)
      IORI(2,N1)=IORI(2,MV+1)
      IORI(3,N1)=IORI(3,MV+1)
      SORI(N1)=SORI(MV+1)
      CLIDER(N1)=CLIDER(MV+1)
      IDPME(N1)=IDPME(MV+1)
      CALL ERAIJ(N1)
      IF(NP1.EQ.0)              CALL REPIJ(N1,MV+1)
      TEX2=TEX2+1.
                 HEX2=HEX2+1.
c-Sob      IF(AN2.LT.ZN2.OR.ZN2.LT.0.)  write(16,300) AN2,ZN2
c-Sob  300 FORMAT(2X,'PAULI3',2(2X,F5.0))
      RETURN
             END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE  PAULID(MV,MD,NRE,IFD,NP,IP)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/NCASCA/NCAS,NCPRI
      COMMON/ACTIM/TINT/TPROD/TPROD(5999)
      COMMON/PORIG/IORI(3,5999)
      COMMON/IDN12/ID1,ID2
      COMMON/CSLID/CLIDER(5999)
      COMMON/NUCOLL/ NUCOLL,MVCOLL(5999)
      COMMON/ACTIV/MPA(240),MYP(5999),MYT(5999),MYY(5999)
      COMMON /MEMAGT/ PME(9,5999),IME(5,5999)
      COMMON /HCASC/ AN1,AN2,ZN1,ZN2,T0,EPS1,EPS2,VPI,A1,A2,C1,C2,D1,D2,
     *R0N1,R0N2,TF01,TF02,RM1,RM2
      COMMON /NUCSP/VPR(3),VTA(3),RADP(3),RADT(3),VEV(3),VRE(3),GEV,GRE,
     *VEP(3),VET(3),GEP,GET
      COMMON/IDPME/ IDPME(5999)
      COMMON /SORI/ SORI(5999),SSOR
      DIMENSION  PS(3),VN(3),RN(3),PL(3),PIN(9),IIN(5)
      MVN=0
      RESM=PME(9,MD)
      DO  108  L=1,NP
      DO   8  K=1,3
      PME(K,MV+L)=PME(K,MD)
    8 CONTINUE
      IF(IME(3,MV+L).EQ.0.AND.IME(4,MV+L).EQ.1)  MVN=MV+L
  108 CONTINUE
      IF(IFD.EQ.0)  GO  TO  113
      IF(MVN.EQ.0)  GO  TO  113
      DO  9  K=1,9
    9 PIN(K)=PME(K,MVN)
      DO 10  K=1,5
   10 IIN(K)=IME(K,MVN)
      DO  13  NU=1,2
      IF(NU.EQ.1.AND.AN1.LT.2.1)   GO  TO  13
      IF(NU.EQ.2.AND.AN2.LT.2.1)   GO  TO  13
      DO  12  K=1,3
      PL(K)=PME(K+3,MVN)
      IF(NU.EQ.2)  GO  TO  11
      VN(K)=-VPR(K)
      RN(K)=RADP(K)
      GO  TO  12
   11 VN(K)=-VTA(K)
      RN(K)=RADT(K)
   12 CONTINUE
      GN=1./DSQRT(1.-VN(1)**2-VN(2)**2-VN(3)**2)
      GG=GN*GN/(GN+1.)
      VR=(PME(1,MVN)-RN(1))*VN(1)+
     +   (PME(2,MVN)-RN(2))*VN(2)+
     +   (PME(3,MVN)-RN(3))*VN(3)
      PIN(1)=PME(1,MVN)-RN(1)+VN(1)*VR*GG
      PIN(2)=PME(2,MVN)-RN(2)+VN(2)*VR*GG
      PIN(3)=PME(3,MVN)-RN(3)+VN(3)*VR*GG
      RMOD=DSQRT(PIN(1)**2+PIN(2)**2+PIN(3)**2)
      IF(NU.EQ.1.AND.RMOD.GT.RM1)  GO  TO  13
      IF(NU.EQ.2.AND.RMOD.GT.RM2)  GO  TO  13
      CALL  CINEMA(PL,VN,PS,CTS,STS,CFS,SFS,TS,PIN(9))
      PIN(4)=PS(1)
      PIN(5)=PS(2)
      PIN(6)=PS(3)
      PIN(8)=TS
      IF(NU.EQ.1) TF=POTEN(PIN,IIN,A1,C1,D1,TF01,VPI,EPS1)-EPS1
      IF(NU.EQ.2) TF=POTEN(PIN,IIN,A2,C2,D2,TF02,VPI,EPS2)-EPS2
      IF(TS.GT.TF)   GO  TO  13
      IP=0
      TAU0= TAUN(NRE)
      TAUL=TAU0*(1.+PME(8,MD)/PME(9,MD))
      PME(7,MD)=0.
      IME(5,MD)=IDINTG((TAUL+PME(7,MD))*1000.)
      IF(IME(5,MD).EQ.0)  IME(5,MD)=1
c-Sob      IF(NCAS.GE.NCPRI) write(16,600) TS,TF,NU,IME(5,MD)
c-Sob  600 FORMAT(1X,'PAULID FORBIDED: TN=',F6.3,' < TF=',F6.3,' NU=',I2,1X,
c-Sob     *'IME(5,MD)=',I10)
      RETURN
   13 CONTINUE
  113 IP=1
      NUCOLL=NUCOLL+1
      DO  114  L=1,NP
      PME(7,MV+L)=0.
      TPROD(MV+L)=TINT
* for decay of L particle the IORI(2,...)=0 !
      IORI(1,MV+L)=ID1
      IORI(2,MV+L)=0
      IORI(3,MV+L)=0
      SORI(MV+L)=RESM
      CLIDER(MV+L)=1.
      MVCOLL(MV+L)=NUCOLL
      MYP(MV+L)=IFD
      MYT(MV+L)=IFD
      MYY(MV+L)=IFD
      IF(IME(2,MV+L).NE.0) MYP(MV+L)=0
      IF(IME(2,MV+L).NE.0) MYT(MV+L)=0
      IF(IME(2,MV+L).NE.0) MYY(MV+L)=0
      MVL=MV+L
c-Sob      IF(NCAS.GE.NCPRI)  write(16,301)
c-Sob     * MVL,(PME(K,MV+L),K=1,9),(IME(K,MV+L),K=1,5),IDPME(MV+L)
c-Sob  301 FORMAT(1X,'NEWP',I5,9(1X,F 8.3),2X,4I2,I15,I5)
  114 CONTINUE
      DO  14  K=1,9
   14 PME(K,MD)=PME(K,MV+NP)
      DO  15  K=1,5
   15 IME(K,MD)=IME(K,MV+NP)
      MYP(MD)=IFD
      MYT(MD)=IFD
      MYY(MD)=IFD
      IF(IME(2,MD).NE.0)  MYP(MD)=0
      IF(IME(2,MD).NE.0)  MYT(MD)=0
      IF(IME(2,MD).NE.0)  MYY(MD)=0
      TPROD(MD)=TINT
* for decay of MD particle the IORI(2,...)=0 !
      IORI(1,MD)=ID1
      IORI(2,MD)=0
      IORI(3,MD)=0
      SORI(MD)=RESM
      MVCOLL(MD)=NUCOLL
      CLIDER(MD)=1.
      IDPME(MD)=IDPME(MV+NP)
      CALL ERAIJ(MD)
      MV=MV+NP-1
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE PAULI4(P1,P2,IP1,IP2,N1,N2,V,NP,MV,TINT,IP)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL *8 P1,P2,V,TINT
C
      REAL*8 MASN
      COMMON/HCASC/ANUCL1,ANUCL2,ZNUCL1,ZNUCL2,T0,EPS1,EPS2,
     *VPI,A1,A2,C1,C2,D1,D2,R0N1,R0N2,TF01,TF02,RM1,RM2
      COMMON/NCASCA/NCAS,NCPRI
      COMMON/TPROD/TPROD(5999)
      COMMON/PORIG/IORI(3,5999)
      COMMON/IDN12/ID1,ID2
      COMMON/CSLID/CLIDER(5999)
      COMMON/NUCOLL/ NUCOLL,MVCOLL(5999)
      COMMON/MEMAGT/PME(9,5999),IME(5,5999)
      COMMON/ACTIV/MPA(240),MYP(5999),MYT(5999),MYY(5999)
      COMMON/NUCSP/VPR(3),VTA(3),RADP(3),RADT(3),VEV(3),VRE(3),GEV,GRE
     *,VEP(3),VET(3),GEP,GET
      COMMON /RINT/ RINT
      COMMON /IDPME/ IDPME(5999)
      COMMON /SORI/ SORI(5999),SSOR
      DIMENSION P1(9),P2(9),IP1(5),IP2(5),PS(3),V(3),PN1(3),
     *PN2(3),P(3),PL(3),BPR(3),BTA(3),RL(3)
      MASN=0.940
      DO  9  L=1,3
      BPR(L)=-VPR(L)
    9 BTA(L)=-VTA(L)
      X=(P1(1)+P2(1))/2.
      Y=(P1(2)+P2(2))/2.
      Z=(P1(3)+P2(3))/2.
      IF((IABS(IP1(4))+IABS(IP2(4))).LE.0)  GO  TO  19
      IF(NP-2)19,10,13
   10 DO 11 L=1,9
   11 PME(L,MV+2)=PME(L,MV+3)
      DO 12 L=1,5
   12 IME(L,MV+2)=IME(L,MV+3)
      CLIDER(MV+2)=CLIDER(MV+3)
      IDPME(MV+2)=IDPME(MV+3)
      GO TO 16
   13 DO 14 L=1,9
      TEMP=PME(L,MV+2)
      PME(L,MV+2)=PME(L,MV+3)
   14 PME(L,MV+3)=TEMP
      DO 15 L=1,5
      ITEMP=IME(L,MV+2)
      IME(L,MV+2)=IME(L,MV+3)
   15 IME(L,MV+3)=ITEMP
      TEMC=CLIDER(MV+2)
      CLIDER(MV+2)=CLIDER(MV+3)
      CLIDER(MV+3)=TEMC
      IDTE=IDPME(MV+2)
      IDPME(MV+2)=IDPME(MV+3)
      IDPME(MV+3)=IDTE
   16 DO 17 L=1,2
      IF(IME(3,MV+L).NE.0.OR.IME(5,MV+L).NE.0.OR.IME(4,MV+L).NE.1)
     *GO  TO  17
      PS(1)=PME(4,MV+L)
      PS(2)=PME(5,MV+L)
      PS(3)=PME(6,MV+L)
      CALL KINEMA(PS,V,PL,CTL,STL,CFL,SFL,TL,MASN)
      IF(ANUCL1.LT.2.1)  GO  TO  117
      CALL KINEMA(PL,BPR,PN1,CT1,ST1,CF1,SF1,T1,MASN)
      CALL RPTS(X,Y,Z,XP,YP,ZP,1)
      TFR1=TFERMI(XP,YP,ZP,1)
      IFORB=1
      IF(T1.LE.TFR1)  GO  TO  18
  117 CALL KINEMA(PL,BTA,PN2,CT2,ST2,CF2,SF2,T2,MASN)
      CALL RPTS(X,Y,Z,XT,YT,ZT,2)
      TFR2=TFERMI(XT,YT,ZT,2)
      IFORB=2
      IF(T2.LE.TFR2)  GO  TO  18
   17 CONTINUE
      GO  TO  19
   18 IP=0
c-Sob      IF(NCAS.GE.NCPRI) THEN
c-Sob        IF(IFORB.EQ.1) write(16,302) T1,TFR1
c-Sob        IF(IFORB.EQ.2) write(16,303) T2,TFR2
c-Sob  302 FORMAT(1X,'PAULI4 FORBIDED: T1=',F6.3,' < TFR1=',F6.3)
c-Sob  303 FORMAT(1X,'PAULI4 FORBIDED: T2=',F6.3,' < TFR2=',F6.3)
c-Sob      ENDIF
           RETURN
   19 IP=1
      NUCOLL=NUCOLL+1
      DO  21  L=1,NP
      M=MV+L
      P(1)=PME(4,M)
      P(2)=PME(5,M)
      P(3)=PME(6,M)
      CALL  KINEMA(P,V,PL,CTL,STL,CFL,SFL,TL,PME(9,M))
      PME(4,M)=PL(1)
      PME(5,M)=PL(2)
      PME(6,M)=PL(3)
      PME(8,M)=TL
      CALL  KINEMR(V,M,RL,TAUL)
      CALL  TMATUR(NP,V,M,P1,P2,MV,P,TAUL,TMAT)
      PME(7,M)=TMAT
      TPROD(M)=TINT
      IORI(1,M)=ID1
      IORI(2,M)=ID2
      IORI(3,M)=0
      SORI(M)=SSOR
      IF(IDPME(N1).EQ.IDPME(M).AND.NP.EQ.2) THEN
* for rescattering of the N1 particle
       IORI(1,M)=IORI(1,N1)
       IORI(2,M)=IORI(2,N1)
       IORI(3,M)=IORI(3,N1)+1
       SORI(M)=SORI(N1)
      ENDIF
      IF(IDPME(N2).EQ.IDPME(M).AND.NP.EQ.2) THEN
* for rescattering of the N2 particle
       IORI(1,M)=IORI(1,N2)
       IORI(2,M)=IORI(2,N2)
       IORI(3,M)=IORI(3,N2)+1
       SORI(M)=SORI(N2)
      ENDIF
      MVCOLL(M)=NUCOLL
      MYP(M)=1
      MYT(M)=1
      MYY(M)=1
      IF(IME(2,M).NE.0)    MYP(M)=0
      IF(IME(2,M).NE.0)    MYT(M)=0
      IF(IME(2,M).NE.0)    MYY(M)=0
      IF(IME(5,M).EQ.0)    GO  TO  20
      TAUM0= DFLOAT(IME(5,M))/1000.
      TAUML=TAUM0*(1.+PME(8,M)/PME(9,M))
      IME(5,M)=IDINTG((TAUML+PME(7,M))*1000.)
      IF(IME(5,M).EQ.0)  IME(5,M)=1
   20 CONTINUE
      DR=RINT*RNDMD(-1)**(1./3.)
      FI=6.283185*RNDMD(-1)
      CT=1.-2.*RNDMD(-1)
      ST=DSQRT(1.-CT**2)
      DX=DR*ST*DCOS(FI)
      DY=DR*ST*DSIN(FI)
      DZ=DR*CT
      PME(1,M)=(P1(1)+P2(1))/2.+DX
      PME(2,M)=(P1(2)+P2(2))/2.+DY
      PME(3,M)=(P1(3)+P2(3))/2.+DZ
c-Sob      IF(NCAS.GE.NCPRI)  write(16,301)
c-Sob     * M,(PME(K,M),K=1,9),(IME(K,M),K=1,5),CLIDER(M),IDPME(M)
c-Sob  301 FORMAT(1X,'NEWP',I5,9(1X,F 8.3),2X,4I2,I15,F6.3,I5)
   21 CONTINUE
      M1=MV+NP-1
      M2=MV+NP
      DO  22  K=1,9
      PME(K,N1)=PME(K,M1)
   22 PME(K,N2)=PME(K,M2)
      MYP(N1)=MYP(M1)
      MYT(N1)=MYT(M1)
      MYY(N1)=MYY(M1)
      TPROD(N1)=TPROD(M1)
      IORI(1,N1)=IORI(1,M1)
      IORI(2,N1)=IORI(2,M1)
      IORI(3,N1)=IORI(3,M1)
      SORI(N1)=SORI(M1)
      MVCOLL(N1)=MVCOLL(M1)
      CLIDER(N1)=CLIDER(M1)
      IDPME(N1)=IDPME(M1)
      DO  23  K=1,5
      IME(K,N1)=IME(K,M1)
   23 IME(K,N2)=IME(K,M2)
      MYP(N2)=MYP(M2)
      MYT(N2)=MYT(M2)
      MYY(N2)=MYY(M2)
      TPROD(N2)=TPROD(M2)
      IORI(1,N2)=IORI(1,M2)
      IORI(2,N2)=IORI(2,M2)
      IORI(3,N2)=IORI(3,M2)
      SORI(N2)=SORI(M2)
      MVCOLL(N2)=MVCOLL(M2)
      CLIDER(N2)=CLIDER(M2)
      IDPME(N2)=IDPME(M2)
      CALL ERAIJ(N1)
      CALL ERAIJ(N2)
      IF(NP.EQ.1) CALL REPIJ(N1,M1)
      MV=MV+NP-2
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE  UPACOV(MV,MVU,IU)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*4 UP1
      COMMON/RESULT/AN1,AN2,ZN1,ZN2,ENEXT1,ENEXT2,PNUCL1(3),PNUCL2(3),
     *AMNUC1(3),AMNUC2(3)
      COMMON/ACTIM/TINT/TPROD/TPROD(5999)
      COMMON/PORIG/IORI(3,5999)
      COMMON/NUCOLL/ NUCOLL,MVCOLL(5999)
      COMMON/CSLID/CLIDER(5999)
      COMMON/BGLAB/BLAB,GLAB,KOBR
      COMMON/UPAC/UP1(66000),IPER
      COMMON/MEMAGT/PM(9,5999),IM(5,5999)
      COMMON/ACTIV/MPA(240),MYP(5999),MYT(5999),MYY(5999)
      COMMON/NPIMI/NPIMI
      COMMON /NCASCA/ NCAS,NCPRI
      COMMON/INTCC/INTCC
      COMMON /IDPME/ IDPME(5999)
      COMMON/INDRES/INDRES
      COMMON/TLIMIT/TLIMIT
      COMMON /SORI/ SORI(5999),SSOR
      MV1=MV
             IF(MV.EQ.0)   RETURN
      M=1
    9 CONTINUE
      IF(IU.EQ.1)   GO  TO  10
      IF(AN1.LT.1.) MYP(M)=0
      IF(AN2.LT.1.) MYT(M)=0
      IF(MYP(M).NE.0.OR.MYT(M).NE.0)                 GO TO 20
C   @@@@@@@@@@@@@@@ FOR ETA MESONS 22.03.94 @@@@@ 15.02.93 @@@@@@@@@@@
      IF(IM(5,M).NE.0.AND.INDRES.EQ.1.AND.
     &(IDPME(M).EQ.220.or.IDPME(M).EQ. 221.or.IDPME(M).EQ.331))
C     IF(IM(5,M).NE.0.AND.INDRES.EQ.1.AND.TINT.GT.TLIMIT)
     *                                               GO TO 10
C   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      IF(IM(5,M).NE.0)                               GO TO 20
      IF(INTCC.NE.0.AND.IM(2,M).EQ.0)                GO TO 20
   10 CONTINUE
      I=MVU+1
      IF((11*I).LT.66000)      GO  TO  201
      IPER=1
c-Sob      write(16,'(15X,''ARRAY  UP1 IS EXCEEDED IN CASCAN''/)')
c-Sob      write( *,'(15X,''ARRAY  UP1 IS EXCEEDED IN CASCAN''/)')
                        RETURN
  201 CONTINUE
      ENM=PM(8,M)+PM(9,M)
      UP1(11*I-10)=DFLOAT(IDPME(M))
      UP1(11*I- 9)=DFLOAT(IORI(1,M))
       ICODE=IORI(3,M)*10**4 + IABS(IORI(2,M))
       ICODE=ISIGN(ICODE,IORI(2,M))
      UP1(11*I- 8)=DFLOAT(ICODE)
      UP1(11*I- 7)=PM(9,M)
      IF(IM(2,M).NE.0)  UP1(11*I- 7)=DFLOAT(IM(2,M))
      UP1(11*I- 6)=PM(1,M)-PM(4,M)/ENM*(TINT-TPROD(M))
      UP1(11*I- 5)=PM(2,M)-PM(5,M)/ENM*(TINT-TPROD(M))
      UP1(11*I- 4)=PM(3,M)-PM(6,M)/ENM*(TINT-TPROD(M))
      UP1(11*I- 3)=PM(4,M)
      UP1(11*I- 2)=PM(5,M)
      UP1(11*I- 1)=PM(6,M)
      IF(KOBR.EQ.1)
     *UP1(11*I- 1)=-GLAB*(PM(6,M)-BLAB*(PM(8,M)+PM(9,M)))
*      UP1(11*I   )=TPROD(M)    !!!  production time
      UP1(11*I   )=SORI(M)
      IF(IM(4,M).EQ.0.AND.IM(1,M).EQ.-1.AND.IM(3,M).EQ.0)
     &  NPIMI=NPIMI+1
      MVU=MVU+1
      MV1=MV1-1
c-Sob      IF(NCAS.GE.NCPRI)
c-Sob     *write(16,202) M,MVU,MV1,(PM(J,M),J=1,9),(IM(JJ,M),JJ=1,5)
c-Sob  202 FORMAT(1X,'UPACOV: M=',I5,2X,'MVU=',I5,2X,'MV1=',I5/
c-Sob     *1X,'==>',9(1X,F 8.3),2X,4I2,I15,'<==')
      IF(M.GT.MV1)   GO  TO  21
      DO  18  J=1,9
   18 PM(J,M)=PM(J,MV1+1)
      DO  19  J=1,5
   19 IM(J,M)=IM(J,MV1+1)
      MYP(M)=MYP(MV1+1)
      MYT(M)=MYT(MV1+1)
      MYY(M)=MYY(MV1+1)
      TPROD(M)=TPROD(MV1+1)
      IORI(1,M)=IORI(1,MV1+1)
      IORI(2,M)=IORI(2,MV1+1)
      IORI(3,M)=IORI(3,MV1+1)
      SORI(M)=SORI(MV1+1)
      MVCOLL(M)=MVCOLL(MV1+1)
      CLIDER(M)=CLIDER(MV1+1)
      IDPME(M)=IDPME(MV1+1)
      CALL REPIJ(M,MV1+1)
      GO  TO  9
   20 IF(M.EQ.MV1)   GO  TO  21
      M=M+1
            GO  TO  9
   21 MV=MV1
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE RPTS(X,Y,Z,X1,Y1,Z1,NU)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL *8 X,Y,Z,X1,Y1,Z1
      COMMON/NUCSP/VPR(3),VTA(3),RADP(3),RADT(3),VEV(3),VRE(3),GEV,GRE
     *,VEP(3),VET(3),GEP,GET
      DIMENSION V(3),R(3)
      DO  11  K=1,3
      IF(NU.EQ.2)  GO  TO  10
      V(K)=VPR(K)
      R(K)=RADP(K)
      GO  TO  11
   10 V(K)=VTA(K)
      R(K)=RADT(K)
   11 CONTINUE
      VR=(X-R(1))*V(1)+(Y-R(2))*V(2)+(Z-R(3))*V(3)
      G=1./DSQRT(1.-V(1)**2-V(2)**2-V(3)**2)
      GG=G*G/(G+1.)
      X1=X-R(1)+V(1)*VR*GG
      Y1=Y-R(2)+V(2)*VR*GG
      Z1=Z-R(3)+V(3)*VR*GG
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      FUNCTION  TFERMI(X,Y,Z,NU)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 TFERMI,X,Y,Z
      COMMON/HCASC/ANUCL1,ANUCL2,ZNUCL1,ZNUCL2,T0,EPS1,EPS2,VPI,A1,A2,
     *C1,C2,D1,D2,R0N1,R0N2,TF01,TF02,RM1,RM2
      COMMON/RESULT/AN1,AN2,ZN1,ZN2,ENEXT1,ENEXT2,PNUCL1(3),PNUCL2(3),
     *AMNUC1(3),AMNUC2(3)
      TFERMI=0.
      IF(NU.EQ.1.AND.ANUCL1.LT.2.1)  RETURN
      IF(NU.EQ.2.AND.ANUCL2.LT.2.1)  RETURN
      IF(NU.EQ.1)  THEN
      AN=ANUCL1
      A=A1
      C=C1
      D=D1
      TF0=TF01
      RM=RM1
      Ares=AN1
      ELSE
      AN=ANUCL2
      A=A2
      C=C2
      D=D2
      TF0=TF02
      RM=RM2
      Ares=AN2
      ENDIF
      R=DSQRT(X**2+Y**2+Z**2)
      IF((R/RM).GT.1.5)  RETURN
      IF(AN.le.10.)  THEN
      TFERMI=TF0*DEXP(-(2./3.)*(R/A)**2) *(Ares/AN)**0.6666667
      ELSE
      TFERMI=TF0*( (1.+DEXP(-A/C))/(1.+DEXP((R-A)/C))
     &  *(Ares/AN) )**0.6666667            ! 09.04.95
      ENDIF
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE  DISSIP(IP,T0)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL *8 T0
      COMMON/NUCSP/VPR(3),VTA(3),RADP(3),RADT(3),VEV(3),VRE(3),GEV,GRE
     *,VEP(3),VET(3),GEP,GET
      COMMON/RESULT/AN1,AN2,ZN1,ZN2,ENEXT1,ENEXT2,PNUCL1(3),PNUCL2(3),
     *AMNUC1(3),AMNUC2(3)
      DIMENSION PPL(3),PTL(3)
      IF(IP.EQ.-1)  GO TO  14
      IF(T0.LT.5.0.AND.AN1.GT.1.1) CALL  COTRAN
      IF(IP)  14,14,9
    9 IF(AN1.LT.0.5)   GO  TO  11
      AMP=0.940*AN1
      CALL  CINEMA(PNUCL1,VPR,PPL,CT,ST,CF,SF,TLP,AMP)
      EP=TLP+AMP
      VPR(1)=PPL(1)/EP
      VPR(2)=PPL(2)/EP
      VPR(3)=PPL(3)/EP
      DO  10  K=1,3
   10 PNUCL1(K)=0.0
   11 IF(AN2.LT.0.5)  GO  TO  14
      AMT=0.940*AN2
      CALL  CINEMA(PNUCL2,VTA,PTL,CT,ST,CF,SF,TLT,AMT)
      ET=TLT+AMT
      VTA(1)=PTL(1)/ET
      VTA(2)=PTL(2)/ET
      VTA(3)=PTL(3)/ET
      DO  12  K=1,3
   12 PNUCL2(K)=0.0
   14 GP=1./DSQRT(1.-VPR(1)**2-VPR(2)**2-VPR(3)**2)
      GT=1./DSQRT(1.-VTA(1)**2-VTA(2)**2-VTA(3)**2)
      VTP=VPR(1)*VTA(1)+VPR(2)*VTA(2)+VPR(3)*VTA(3)
      DO  15  K=1,3
      VRE(K)=(VPR(K)+GT*VTA(K)*(VTP*GT/(GT+1.)-1.))/(1.-VTP)/GT
   15 VEV(K)=(VPR(K)*GP+VTA(K)*GT)/(GT+GP)
      GEV=1./DSQRT(1.-VEV(1)**2-VEV(2)**2-VEV(3)**2)
      GRE=1./DSQRT(1.-VRE(1)**2-VRE(2)**2-VRE(3)**2)
      SET=VEV(1)*VTA(1)+VEV(2)*VTA(2)+VEV(3)*VTA(3)
      SEP=VEV(1)*VPR(1)+VEV(2)*VPR(2)+VEV(3)*VPR(3)
      DO 115 K=1,3
      VET(K)=(VEV(K)+GT*VTA(K)*(SET*GT/(GT+1.)-1.))/(1.-SET)/GT
      VEP(K)=(VEV(K)+GP*VPR(K)*(SEP*GP/(GP+1.)-1.))/(1.-SEP)/GP
  115 CONTINUE
      GET=GEV*GT*(1.-SET)
      GEP=GEV*GP*(1.-SEP)
   16 CONTINUE
  200 FORMAT(5X,'VPR,VTA,GP,GT',2(3(F7.4,2X),5X),2(F8.4,2X))
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE  CENUM1(NA,P,IP,DELTA,NC,K1,K2,IK,NU)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL *8 P,DELTA
      COMMON/CENTER/XC(2,240),YC(2,240),ZC(2,240),IZ(2,240)
      COMMON/CENPAR/NUKC(100)
      COMMON/HCASC/ANUCL1,ANUCL2,ZNUCL1,ZNUCL2,T0,EPS1,EPS2,
     *VPI,A1,A2,C1,C2,D1,D2,R0N1,R0N2,TF01,TF02,RM1,RM2
      DIMENSION   P(9),RIN1(3),RC1(3),IP(5)
     *,ZKC(100)
      PM=DSQRT(P(4)**2+P(5)**2+P(6)**2)
      IF(IK)  9,9,16
    9 CONTINUE
      NC=0
      K=1
      K2=0
      DELTA2=DELTA**2
      R=RM2
      IF(NU.EQ.1)   R=RM1
      RIN1(3)=(P(1)*P(4)+P(2)*P(5)+P(3)*P(6))/PM
      DEV=DSQRT(P(1)**2+P(2)**2+P(3)**2-RIN1(3)**2)
      IF(DEV.GT.(R+DELTA))   GO  TO  26
  104 RC1(3)=(XC(NU,K)*P(4)+YC(NU,K)*P(5)+ZC(NU,K)*P(6))/PM
      DR1=RC1(3)-RIN1(3)
      IF(DR1.LT.0.)   GO  TO  14
      IF(DABS(DR1).LT.0.0001)   GO  TO  14
      DEL2=(XC(NU,K)-P(1))**2+(YC(NU,K)-P(2))**2+(ZC(NU,K)-P(3))**2-
     -DR1**2
      IF(DEL2-DELTA2)   11,11,14
   11 Z1=RC1(3)
      KA=K
      IF(NC.EQ.0)  GO  TO  113
      J=1
   12 IF(Z1.GT.ZKC(J))  GO  TO  13
      Z=ZKC(J)
      ZKC(J)=Z1
      Z1=Z
      KC=NUKC(J)
      NUKC(J)=KA
      KA=KC
   13 IF(J.EQ.NC)   GO  TO  113
      J=J+1
            GO  TO  12
  113 IF(NC.EQ.100)   GO  TO  114
      NC=NC+1
      ZKC(NC)=Z1
      NUKC(NC)=KA
                   GO  TO  14
  114 K2=1
   14 IF(K-NA)  15,26,26
   15 K=K+1
      GO  TO  104
   16 IF(NA-1) 107,107,108
  107 K2=0
           GO TO 26
  108 IF(IP(4).EQ.1)  GO  TO  107
   17 J=1
   18 IF(K1-1)  20,19,20
   19 J=2
   20 R=DSQRT((XC(NU,K1)-XC(NU,J))**2+(YC(NU,K1)-YC(NU,J))**2+
     *(ZC(NU,K1)-ZC(NU,J))**2)
      K2=J
   21 IF(J-K1)  22,24,22
   22 RJ=DSQRT((XC(NU,K1)-XC(NU,J))**2+(YC(NU,K1)-YC(NU,J))**2+
     *(ZC(NU,K1)-ZC(NU,J))**2)
      IF(RJ-R)  23,24,24
   23 R=RJ
      K2=J
   24 IF(J-NA)  25,26,26
   25 J=J+1
              GO  TO  21
   26 RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE  ABEL(PIN,V,U,P1,P2,CT,FI,CM1,CM2)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL *8 PIN,V,U,P1,P2,CT,FI,CM1,CM2
      DIMENSION  PIN(9),V(3),B(3),P1(3),P2(3),PINL(3),PINS(3)
      B(1)=-V(1)
      B(2)=-V(2)
      B(3)=-V(3)
      PINL(1)=PIN(4)
      PINL(2)=PIN(5)
      PINL(3)=PIN(6)
      CALL  CINEMA(PINL,B,PINS,CTS,STS,CFS,SFS,TS,PIN(9))
      ST=DSQRT(1.-CT*CT)
      CF=DCOS(FI)
      SF=DSIN(FI)
      E1=(U*U+CM1**2-CM2**2)/(2.*U)
      T1=DABS(E1-CM1)
      P1M=DSQRT(T1*(T1+2.*CM1))
      P2(1)=P1M*ST*CF
      P2(2)=P1M*ST*SF
      P2(3)=P1M*CT
      CALL  ROTOR(PINS,V,P2,P1)
      P2(1)=-P1(1)
      P2(2)=-P1(2)
      P2(3)=-P1(3)
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE  TINVU(P1,P2,U,V,TIN1)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL *8 P1,P2,U,V,TIN1
      DIMENSION  V(3),P1(9),P2(9)
      E1=P1(8)+P1(9)
      E2=P2(8)+P2(9)
      V(1)=(P1(4)+P2(4))/(E1+E2)
      V(2)=(P1(5)+P2(5))/(E1+E2)
      V(3)=(P1(6)+P2(6))/(E1+E2)
      U=(E1+E2)*DSQRT(1.-V(1)**2-V(2)**2-V(3)**2)
      TIN1=(U**2-(P1(9)+P2(9))**2)/(2.*P2(9))
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE ABSORP(PARTIN,IPATIN,PARTNE,PAR1,NE,IE,MV,NP,V,U)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 PARTIN,PARTNE,PAR1,V,U
      REAL*8 MASN
C     BLOCK OF CALCULATION OF PARTICLE CHARACTERISTICS IN ABSORPTION.
      COMMON/MEMAGT/PMEMO,IMEMO
      DIMENSION PARTIN(9),IPATIN(5),PARTNE(9),PAR1(9),
     *PAF(9),V(3),PIST(3),PNST(3),PMEMO(9,5999),IMEMO(5,5999)
      MASN=0.940
      IF(IPATIN(5).NE.0.AND.IPATIN(4).EQ.1)  GO  TO  13
      PAF(4)=PARTNE(4)+PAR1(4)
      PAF(5)=PARTNE(5)+PAR1(5)
      PAF(6)=PARTNE(6)+PAR1(6)
      PAF(9)=2.*PAR1(9)
      PAF(8)=DSQRT(PAF(4)**2+PAF(5)**2+PAF(6)**2+PAF(9)**2)-PAF(9)
      CALL  TINVU(PARTIN,PAF,U,V,TIN1)
      IF (IPATIN(2)) 14,13,14
   13 CTST = 1.-2.*RNDMD(-1)
      FIST = 6.283185*RNDMD(-1)
                         GO TO 16
   14 IF (TIN1-0.455) 15,15,13
   15 CTST = COSTA(19,TIN1)
      FIST = 6.283185*RNDMD(-1)
                          GO TO 16
   16 CALL ABEL(PARTIN,V,U,PIST,PNST,CTST,FIST,MASN,MASN)
      IF(MV-5997)17,17,18
   18 NP = 0
c-Sob      write(16,19)
c-Sob   19 FORMAT (25X,'MEMAGT IS EXCEEDED IN CASCAD')
      RETURN
   17 CONTINUE
      PMEMO(1,MV+3) = 0.
      PMEMO(2,MV+3) = 0.
      PMEMO(3,MV+3) = 0.
      PMEMO(4,MV+3) = PIST(1)
      PMEMO(5,MV+3) = PIST(2)
      PMEMO(6,MV+3) = PIST(3)
      PMEMO(7,MV+3) = 0.
      PMEMO(9,MV+3) = 0.940
      IMEMO(1,MV+3) = IE
      IMEMO(2,MV+3) = 0
      IMEMO(3,MV+3) = 0
      IMEMO(4,MV+3) = 1
      IMEMO(5,MV+3) = 0
      PMEMO(1,MV+1) = 0.
      PMEMO(2,MV+1) = 0.
      PMEMO(3,MV+1) = 0.
      PMEMO(4,MV+1) = -PIST(1)
      PMEMO(5,MV+1) = -PIST(2)
      PMEMO(6,MV+1) = -PIST(3)
      PMEMO(7,MV+1) = 0.
      PMEMO(9,MV+1) = 0.940
      IMEMO(1,MV+1) = NE
      IMEMO(2,MV+1) = 0
      IMEMO(3,MV+1) = 0
      IMEMO(4,MV+1) = 1
      IMEMO(5,MV+1) = 0
      NP = 2
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE PARTN(NU,I,P,IP)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 P
      COMMON/CENTER/XC(2,240),YC(2,240),ZC(2,240),IZ(2,240)
     */HCASC/ANUCL1,ANUCL2,ZNUCL1,ZNUCL2,T0,EPS1,EPS2,
     *VPI,A1,A2,C1,C2,D1,D2,R0N1,R0N2,TF01,TF02,RM1,RM2
      COMMON/RESULT/AN1,AN2,ZN1,ZN2,ENEXT1,ENEXT2,PNUCL1(3),PNUCL2(3),
     *AMNUC1(3),AMNUC2(3)
      DIMENSION P(9),IP(5)
      P(1)=XC(NU,I)
      P(2)=YC(NU,I)
      P(3)=ZC(NU,I)
      R = DSQRT(P(1)**2+P(2)**2+P(3)**2)
      R=RPOTEN(R)
      IF(NU-1) 10,10,11
   10 IF(ANUCL1-10.) 120,120,110
  110 TFR=TF01*(((1.+DEXP(-A1/C1))/(1.+DEXP((R-A1)/C1)))**0.6666667 )
      GO TO 12
  120 TFR=TF01*DEXP(-(2./3.)*(R/A1)**2)
                               GO TO 12
   11 IF(ANUCL2-10.) 121,121,111
  111 TFR = TF02*(((1.+DEXP(-A2/C2))/(1.+DEXP((R-A2)/C2)))**.6666667)
      GO TO 12
  121 TFR=TF02*DEXP(-(2./3.)*(R/A2)**2)
                               GO TO 12
   12 CONTINUE
c                    !!! 20.06.1995
      IF(NU.EQ.1) TFR=TFR*(AN1/ANUCL1)**(2./3.)
      IF(NU.EQ.2) TFR=TFR*(AN2/ANUCL2)**(2./3.)
c
      TN = TFR*(RNDMD(-1)**(2./3.))
      PN=DSQRT(TN*(TN+1.88))
      CT=1.-2.*RNDMD(-1)
      ST=DSQRT(1.-CT*CT)
      FI=6.283185*RNDMD(-1)
      CF=DCOS(FI)
      SF=DSIN(FI)
      P(4)=PN*ST*CF
      P(5)=PN*ST*SF
      P(6)=PN*CT
      P(7)=0.
      P(8)=TN
      P(9)=0.94
      IP(2)=0
      IP(3)=0
      IP(4)=1
      IP(1)=IZ(NU,I)
      IP(5)=0
      IF(NU-1) 13,13,14
   13 P(9)=P(9)-POTEN(P,IP,A1,C1,D1,TF01,VPI,EPS1)+P(8)
                                            GO TO 15
   14 P(9)=P(9)-POTEN(P,IP,A2,C2,D2,TF02,VPI,EPS2)+P(8)
                                            GO TO 15
   15 CONTINUE
      RETURN
             END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE DIRECT (V,TIN1,MQ,MV,NP,PARTIN,KP,ITH)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 V,TIN1,PARTIN
C     DETERMINING OF DIRECTION OF SECONDARY PARTICLES MOTION.
      DIMENSION PMEMO(9,5999),IMEMO(5,5999),V(3),PAKV(3),PARTIN(9),
     1PLST(3),PL(3),PAKST(3),PIN(3)
     *,B(3),PIL(3)
      COMMON/MEMAGT/PMEMO,IMEMO
      ND = 0
      KP = 0
      IF (MQ-1) 12,12,10
   10 IF (RNDMD(-1) -0.5) 15,11,11
   11 IF (RNDMD(-1) -0.5) 14,16,16
   12 IF (RNDMD(-1) -0.5) 14,13,13
   13 IF (RNDMD(-1) -0.5) 15,16,16
   14 M1 = 2
      M2 = 3
                  GO TO 17
   15 M1 = 1
      M2 = 3
                  GO TO 17
   16 M1 = 1
      M2 = 2
                  GO TO 17
   17 LAMBDA=1
      PAKV(1)=0.
      PAKV(2)=0.
      PAKV(3)=0.
      M1TEMP = MV+M1
      M2TEMP = MV+M2
   18 LTEMP = MV+LAMBDA
      IF (LAMBDA-M1) 19,21,19
   19 IF (LAMBDA-M2) 20,21,20
   20 JA = JTYPA(ITH,MQ,LAMBDA)
      CTL = COSTA(JA,TIN1)
      FL = 6.283185*RNDMD(-1)
      STL = DSQRT (1.-CTL**2)
      TEMP1=DCOS(FL)
      TEMP2=DSIN(FL)
      TEMP3=PMEMO(8,LTEMP)
      PMEMO(1,LTEMP) = 0.
      PMEMO(2,LTEMP) = 0.
      PMEMO(3,LTEMP) = 0.
      PMEMO(4,LTEMP) = TEMP3*STL*TEMP1
      PMEMO(5,LTEMP) = TEMP3*STL*TEMP2
      PMEMO(6,LTEMP) = TEMP3*CTL
      PMEMO(7,LTEMP) = 0.
      PAKV(1) = PAKV(1)+PMEMO(4,LTEMP)
      PAKV(2) = PAKV(2)+PMEMO(5,LTEMP)
      PAKV(3) = PAKV(3)+PMEMO(6,LTEMP)
   21 IF (LAMBDA-NP) 22,23,23
   22 LAMBDA = LAMBDA+1
                    GO TO 18
   23 PAKVM = DSQRT (PAKV(1)**2+PAKV(2)**2+PAKV(3)**2)
      IF (NP-3) 25,24,25
   24 PIL(1)=PARTIN(4)
      PIL(2)=PARTIN(5)
      PIL(3)=PARTIN(6)
      B(1)=-V(1)
      B(2)=-V(2)
      B(3)=-V(3)
      CALL  CINEMA(PIL,B,PIN,CTS,STS,CFS,SFS,TIN,PARTIN(9))
      LAMBDA = 1
               GO TO 27
   25 IF (PMEMO(8,M1TEMP)-PAKVM-PMEMO(8,M2TEMP)) 26,32,32
   26 IF (PMEMO(8,M1TEMP)-DABS(PAKVM-PMEMO(8,M2TEMP))) 32,32,24
   27 LTEMP = MV+LAMBDA
      IF (LAMBDA-M1) 29,34,29
   28 LAMBDA = LAMBDA+1
                    GO TO 27
   29 IF (LAMBDA-M2) 30,34,30
   30 PL(1) = PMEMO(4,LTEMP)
      PL(2) = PMEMO(5,LTEMP)
      PL(3) = PMEMO(6,LTEMP)
      CALL ROTOR (PIN,V,PL,PLST)
      PMEMO(1,LTEMP) = 0.
      PMEMO(2,LTEMP) = 0.
      PMEMO(3,LTEMP) = 0.
      PMEMO(4,LTEMP) = PLST(1)
      PMEMO(5,LTEMP) = PLST(2)
      PMEMO(6,LTEMP) = PLST(3)
      PMEMO(7,LTEMP) = 0.
   34 IF (LAMBDA-NP) 28,31,31
   31 CALL ROTOR (PIN,V,PAKV,PAKST)
      CTM1 = (PMEMO(8,M2TEMP)**2-PMEMO(8,M1TEMP)**2-PAKVM**2)/
     1(2.*PAKVM*PMEMO(8,M1TEMP))
      CTM2 = (PMEMO(8,M1TEMP)**2-PMEMO(8,M2TEMP)**2-PAKVM**2)/
     1(2.*PAKVM*PMEMO(8,M2TEMP))
      FM1 = 6.283185*RNDMD(-1)
      FM2 = 3.141592+FM1
      STM1 = DSQRT (1.-CTM1**2)
      STM2 = DSQRT (1.-CTM2**2)
      CFM1 = DCOS(FM1)
      SFM1 = DSIN(FM1)
      CFM2 = DCOS(FM2)
      SFM2 = DSIN(FM2)
      PL(1) = PMEMO(8,M1TEMP)*STM1*CFM1
      PL(2) = PMEMO(8,M1TEMP)*STM1*SFM1
      PL(3) = PMEMO(8,M1TEMP)*CTM1
      CALL ROTOR (PAKST,V,PL,PLST)
      PMEMO(1,M1TEMP) = 0.
      PMEMO(2,M1TEMP) = 0.
      PMEMO(3,M1TEMP) = 0.
      PMEMO(4,M1TEMP) = PLST(1)
      PMEMO(5,M1TEMP) = PLST(2)
      PMEMO(6,M1TEMP) = PLST(3)
      PMEMO(7,M1TEMP) = 0.
      PL(1) = PMEMO(8,M2TEMP)*STM2*CFM2
      PL(2) = PMEMO(8,M2TEMP)*STM2*SFM2
      PL(3) = PMEMO(8,M2TEMP)*CTM2
      CALL ROTOR (PAKST,V,PL,PLST)
      PMEMO(1,M2TEMP) = 0.
      PMEMO(2,M2TEMP) = 0.
      PMEMO(3,M2TEMP) = 0.
      PMEMO(4,M2TEMP) = PLST(1)
      PMEMO(5,M2TEMP) = PLST(2)
      PMEMO(6,M2TEMP) = PLST(3)
      PMEMO(7,M2TEMP) = 0.
      RETURN
   32 ND = ND+1
      IF (ND-100) 17,33,33
   33 KP = 2
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE PINPN(RN,R0X,R0Y,R0Z,MV,DELTA,NEL)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 RN,R0X,R0Y,R0Z,DELTA
      CHARACTER*8 PNA1
      COMMON/HCASC/ANUCL1,ANUCL2,ZNUCL1,ZNUCL2,T0,EPS1,EPS2,
     *VPI,A1,A2,C1,C2,D1,D2,R0N1,R0N2,TF01,TF02,RM1,RM2
     */CENTER/XC(2,240),YC(2,240),ZC(2,240),IZ(2,240)
      COMMON/ACTIV/MPA(240),MYP(5999),MYT(5999),MYY(5999)
      COMMON/RESULT/AN1,AN2,ZN1,ZN2,ENEXT1,ENEXT2,PNUCL1(3),PNUCL2(3),
     *AMNUC1(3),AMNUC2(3)
       COMMON/MEMAGT/PMEMO(9,5999),IMEMO(5,5999)
C
      COMMON/XBMAX/XBMAX,IFIB0
      COMMON /STIN/ STIN,AMIN
      COMMON/CSLID/CLIDER(5999)
      COMMON /IDPME/ IDPME(5999)
      DIMENSION PARTIN(9),IPATIN(5)
      T=(2.*RN)**2
      IF(ANUCL1-1.1) 98,98,99
   98 RM1=0.
   99 IF(ANUCL1.EQ.2.)  RM1=2.158
      R12=RM1+RM2+DELTA+1./(5.06*DSQRT(T0*(T0+1.88)))
      CALL  RXYZ(R12,R0X,R0Y,R0Z)
      IF(ANUCL1-1.) 100,100,101
  100 PMEMO(1,1)=R0X
      PMEMO(2,1)=R0Y
      PMEMO(3,1)=R0Z
      IMEMO(1,1)=ZNUCL1
      IMEMO(2,1)=0
      IMEMO(3,1)=IDINTG(STIN)
      IMEMO(4,1)=IDINTG(ANUCL1)
      IF(ANUCL1.GT.0.5)  IMEMO(4,1)=1
      IMEMO(5,1)=0
      PMEMO(9,1)=AMIN
      PMEMO(4,1)=0.
      PMEMO(5,1)=0.
      PMEMO(6,1)=DSQRT(T0*(T0+2.*PMEMO(9,1)))
      PMEMO(7,1)=0.
      PMEMO(8,1)=T0
      MYP(1)=0
      MYT(1)=1
      MYY(1)=1
      CLIDER(1)=1.
      DO  1 K=1,9
      PARTIN(K)=PMEMO(K,1)
      IF(K.LE.5)  IPATIN(K)=IMEMO(K,1)
    1 CONTINUE
      CALL PANUN(PARTIN,IPATIN,IK1)
      CALL IDPANU(ID1,IK1,PNA1)
      IDPME(1)=ID1
      MV=1
      AN1=0.
      ZN1=0.
                   GO  TO  102
  101 IF(ANUCL1.GT.2.1)   GO  TO  104
      G=1.+T0/0.940
      CT1=1.-2.*RNDMD(-1)
      FI1=6.283185*RNDMD(-1)
      ST1=DSQRT(1.-CT1*CT1)
      Z1=RM1*CT1/G
      Y1=RM1*ST1*DSIN(FI1)
      X1=RM1*ST1*DCOS(FI1)
      PMEMO(1,1)=R0X+X1
      PMEMO(2,1)=R0Y+Y1
      PMEMO(3,1)=R0Z+Z1
      PMEMO(1,2)=R0X-X1
      PMEMO(2,2)=R0Y-Y1
      PMEMO(3,2)=R0Z-Z1
  103 AL=10.*RNDMD(-1)
      BR=RNDMD(-1)
      FBR=4.*AL*AL/(1.+AL*AL)**2
      IF(BR.GT.FBR)   GO  TO  103
      EPM=0.940*0.00223
      PM1=DSQRT(EPM)*AL
      CT1=1.-RNDMD(-1)
      ST1=DSQRT(1.-CT1*CT1)
      FI1=6.283185*RNDMD(-1)
      PX1=PM1*ST1*DCOS(FI1)
      PY1=PM1*ST1*DSIN(FI1)
      PZ1=PM1*CT1
      PMEMO(4,1)=PX1
      PMEMO(4,2)=-PX1
      PMEMO(5,1)=PY1
      PMEMO(5,2)=-PY1
      PMEMO(6,1)=G*PZ1+DSQRT((PM1*PM1+0.94*0.94)*(G*G-1.))
      PMEMO(6,2)=PMEMO(6,1)-2.*G*PZ1
      PMEMO(7,1)=0.
      PMEMO(7,2)=0.
      PMEMO(8,1)=G*DSQRT(PM1*PM1+0.94*0.94)+PZ1*DSQRT(G*G-1.)-0.940
      PMEMO(8,2)=PMEMO(8,1)-2.*PZ1*DSQRT(G*G-1.)
      PMEMO(9,1)=0.940
      PMEMO(9,2)=0.940
      IMEMO(1,1)=0
      IDPME(1)=1220
      IMEMO(2,1)=0
      IMEMO(3,1)=0
      IMEMO(4,1)=1
      IMEMO(1,2)=1
      IDPME(2)=1120
      IMEMO(2,2)=0
      IMEMO(3,2)=0
      IMEMO(4,2)=1
      IMEMO(5,1)=0
      IMEMO(5,2)=0
      MYP(1)=0
      MYP(2)=0
      MYT(1)=1
      MYY(1)=1
      MYT(2)=1
      MYY(2)=1
      CLIDER(1)=1.
      CLIDER(2)=1.
      MV=2
      AN1=0
      ZN1=0
                   GO  TO  102
  104 CONTINUE
      IF(NEL.NE.0)  CALL  ROTNUC(1,ANUCL1)
      IF(NEL.NE.0)  GO TO 102
      IM=ANUCL1
      NZ1=ZNUCL1
      DO 21 I=1,IM
   11 B1=RNDMD(-1)
      RI=RM1*(B1**(1./3.))
      IF(ANUCL1-10.) 12,12,13
   12 FB=DEXP(-(RI**2)/(A1**2))
                         GO TO 14
   13 FB=(1.+DEXP(-A1/C1))/(1.+DEXP((RI-A1)/C1))
   14 DRND=RNDMD(-1)
      IF(DRND-FB) 15,15,11
   15 CT=1.-2.*RNDMD(-1)
      FI=6.283185*RNDMD(-1)
      ST=DSQRT(1.-CT**2)
      RS=RI*ST
      XC(1,I)=RS*DCOS(FI)
      YC(1,I)=RS*DSIN(FI)
      ZC(1,I)=RI*CT
      IF(I-NZ1) 16,16,17
   16 IZ(1,I)=1
              GO TO 18
   17 IZ(1,I)=0
   18 IF(I-1) 221,221,19
   19 KM=I-1
      DO 20 K=1,KM
      IF((XC(1,I)-XC(1,K))**2+(YC(1,I)-YC(1,K))**2+(ZC(1,I)-ZC(1,K))**2
     *-T) 11,20,20
   20 CONTINUE
  221 CONTINUE
   21 CONTINUE
  102 IF(ANUCL2.GT.1.1)   GO  TO  114
      PMEMO(1,1)=0.
      PMEMO(2,1)=0.
      PMEMO(3,1)=1.E-20
      PMEMO(4,1)=0.
      PMEMO(5,1)=0.
      PMEMO(6,1)=1.E-20
      PMEMO(7,1)=0.
      PMEMO(8,1)=0.
      PMEMO(9,1)=0.940
      IF(ANUCL2.LT.0.5)  PMEMO(9,1)=0.140
      IMEMO(1,1)=ZNUCL2
      IDPME(2)=1120
      IF(ZNUCL2.LE.0.5) IDPME(2)=1220
      IMEMO(2,1)=0
      IMEMO(3,1)=0
      IMEMO(4,1)=ANUCL2
      IMEMO(5,1)=0
      CLIDER(1)=1.
      AN2=0.
      ZN2=0.
      MV=1
      MYP(1)=1
      MYT(1)=0
      MYY(1)=1
      RETURN
  114 CONTINUE
      IF(NEL.NE.0)  GO TO 105
      IM=ANUCL2
      NZ2=ZNUCL2
      DO 33 I=1,IM
   22 B1=RNDMD(-1)
      RI=RM2*(B1**(1./3.))
      IF(ANUCL2-10.) 23,23,24
   23 FB=DEXP(-(RI**2)/(A2**2))
                         GO TO 25
   24 FB=(1.+DEXP(-A2/C2))/(1.+DEXP((RI-A2)/C2))
   25 DRND=RNDMD(-1)
      IF(DRND-FB) 26,26,22
   26 CT=1.-2.*RNDMD(-1)
      FI=6.283185*RNDMD(-1)
      ST=DSQRT(1.-CT**2)
      RS=RI*ST
      XC(2,I)=RS*DCOS(FI)
      YC(2,I)=RS*DSIN(FI)
      ZC(2,I)=RI*CT
      IF (I-NZ2) 27,27,28
   27 IZ(2,I)=1
              GO TO 29
   28 IZ(2,I)=0
   29 IF(I-1) 33,33,30
   30 KM=I-1
      DO 32 K=1,KM
      IF((XC(2,I)-XC(2,K))**2+(YC(2,I)-YC(2,K))**2+(ZC(2,I)-ZC(2,K))**2
     *-T) 22,32,32
   32 CONTINUE
   33 CONTINUE
      RETURN
  105 CALL  ROTNUC(2,ANUCL2)
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE ROTNUC(N,AN)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 AN
      COMMON/CENTER/XC(2,240),YC(2,240),ZC(2,240),IZ(2,240)
      DIMENSION A(3),B(3),R(3),RR(3)
      IF(AN.LT.2.1)   RETURN
      CT=1.-2.*RNDMD(-1)
      ST=DSQRT(1.-CT*CT)
      FI=6.283185*RNDMD(-1)
      A(1)=ST*DCOS(FI)
      A(2)=ST*DSIN(FI)
      A(3)=CT
      B(1)=0.
      B(2)=0.
      B(3)=1.
      IA=AN+0.1
      DO 10 I=1,IA
      R(1)=XC(N,I)
      R(2)=YC(N,I)
      R(3)=ZC(N,I)
      CALL ROTOR(A,B,R,RR)
      XC(N,I)=RR(1)
      YC(N,I)=RR(2)
      ZC(N,I)=RR(3)
   10 CONTINUE
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE   TYPNEW(PARTIN,IPATIN,PARTNE,IPATNE,V,U,TIN1,
     *             SABS,MV,NP,NABS,PAR1,IPA1,N3,NU,N2,NA1,NA2)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 PARTIN,PARTNE,V,U,TIN1,SABS,PAR1
C     REAL*4 ATYP(4)
C     REAL*8 DTYP(11),TDIAG,PNA1,PNA2
      CHARACTER*4 ATYP(4)
      CHARACTER*8 DTYP(11),TDIAG,PNA1,PNA2,PNAK
      LOGICAL YESELA
      COMMON/YESELA/ YESELA
      COMMON/NCASCA/NCAS,NCPRI
      COMMON/INTTYP/ ITYP
      COMMON/ITHEA/ITHEA(11)
      COMMON/LOWMIS/LOWMIS
      COMMON/IACT/ IACT
      COMMON /MEMAGT/ PME(9,5999),IME(5,5999)
      COMMON /HELE/ IHELE
      COMMON /IDN12/ ID1,ID2
      COMMON /IDN120/ ID10,ID20
      COMMON /IDPME/ IDPME(5999)
      COMMON/COMELX/ SEL
      COMMON/COMCRO/ STO
      COMMON /DATA2/ PUD,PS1,SIGMA,CX2
      COMMON /SORI/ SORI(5999),SSOR
      COMMON /ISLERR/ ISLERR   ! HADES1
      DIMENSION PARTIN(9),IPATIN(5),PARTNE(9),IPATNE(5),V(3),
     *PAR1(9),IPA1(5),PA(9),IPA(5)
      DATA ATYP/'ABS:','ELE:','BIN:','HEI:'/
      DATA DTYP/'DIFTRI','PLANAR','UNCYLI','ELAST  ','ANNIH','DIFSMA',
     *          'CHAINS','BINAR ','??????','REGTRI ','DOUBDI'/
C
C     BLOCK OF DETERMINING OF INTERACTION TYPE AND CALCULATION
C     OF SECONDARY PARTICLES CHARACTERISTICS.
C
      AM1=PARTIN(9)
      SSOR=U
   11 CONTINUE
      ID10=ID1
      ID20=ID2
      LOWMIS=0
      ITYP=0
      NABS=0
      YESELA=.TRUE.
      CALL PANUID(ID1,IK1,PNA1)
      CALL PANUID(ID2,IK2,PNA2)
c-Sob      IF(NCAS.GE.NCPRI) THEN
c-Sob        write(16,599) PNA1,ID1,PNA2,ID2,TIN1
c-Sob        write( *,599) PNA1,ID1,PNA2,ID2,TIN1
c-Sob  599 FORMAT(1X,A4,'(',I5,')','+',A4,'(',I5,')','  TIN1=',F7.3)
c-Sob      ENDIF
c-Sob      IF(TIN1.LE.0) THEN
c-Sob        write(16,*) 'PI=', (PARTIN(KK),KK=4,9)
c-Sob        write(16,*) 'PN=', (PARTNE(KK),KK=4,9)
c-Sob      ENDIF
  112 CONTINUE
      IF(IK2.LT.37.OR.IK2.GT.38)           GO  TO  18
      IF(IK1.EQ.1.OR.IK1.EQ.2.OR.IK1.EQ.7) GO  TO  99
      IF(IK1.LT.37.OR.IK1.GT.38)           GO  TO  18
                                 GO  TO  13
   99 CONTINUE
      CALL SLQEK (L,MS,MQ,KSI,ME,IPATIN,IPATNE)
      STO=CROSEG(L,MS,MQ,KSI,0,TIN1,AM1,IPATIN(5))
      BETABS=SABS/(STO+SABS)
      DRND=RNDMD(-1)
      IF(DRND-BETABS) 101,101,13
  101 CONTINUE
      DL=0.
      IF(NU.EQ.1)  CALL  CENUM1(NA1,PARTIN,IPATIN,DL,NC,N2,N3,N2,1)
      IF(NU.EQ.2)  CALL  CENUM1(NA2,PARTIN,IPATIN,DL,NC,N2,N3,N2,2)
      IF(N3.EQ.0)  GO  TO  13
      CALL  PARTN(NU,N3,PAR1,IPA1)
      PAR1(9) = PARTNE(9)
      IE2=IPATNE(1)
      IE3=IPA1(1)
      IE1=IPATIN(1)
  104 IF(IE1) 108,103,105
  103 NE1=IE2
      NE2=IE3
                  GO TO 14
  105 IF(IE2+IE3-1) 107,106,13
  106 NE1=1
      NE2=1
                  GO TO 14
  107 NE1=1
      NE2=0
                  GO TO 14
  108 IF(IE2+IE3-1) 13,109,110
  109 NE1=0
      NE2=0
                  GO TO 14
  110 NE1=0
      NE2=1
   14 CONTINUE
      CALL ABSORP(PARTIN,IPATIN,PARTNE,PAR1,NE1,NE2,MV,NP,V,U)
      NABS = 1
      ITYP=1
      GO  TO  27
   13 CONTINUE
C    ****** TIN1=0.6  IS INTRODUCED *******
***      IF(IACT.GE.2.AND.TIN1.GT.0.600)    GO  TO  18
      IF(IACT.GE.2.AND.TIN1.GT.1.000)    GO  TO  18
      CALL SLQEK (L,MS,MQ,KSI,ME,IPATIN,IPATNE)
      STO=CROSEG(L,MS,MQ,KSI,0,TIN1,AM1,IPATIN(5))
      SEL=CROSEG(L,MS,MQ,KSI,1,TIN1,AM1,IPATIN(5))
      SEX=CROSEG(L,MS,MQ,KSI,2,TIN1,AM1,IPATIN(5))
      BETAEL=(SEL+SEX)/STO
      DRND=RNDMD(-1)
      IF (DRND-BETAEL) 16,16,15
   16 CONTINUE
      CALL ELEX(V,U,TIN1,PARTIN,IPATIN,IPATNE,MV,NP,L,MS,
     *MQ,KSI,ME)
      ITYP=2
      GO  TO  27
   15 CONTINUE
      IF(TIN1.LE.5.0)  GO  TO  20
      YESELA=.FALSE.
      GO  TO  18
   20 CONTINUE
      CALL BINEL(PARTIN,IPATIN,IPATNE,L,MS,MQ,KSI,ME,V,U,TIN1,
     *MV,NP,NIN)
      ITYP=3
      IF (NIN.NE.0)   GO  TO  11
      GO  TO   27
   18 CONTINUE
C
      CALL PINETA(ID1,ID2,V,U,PARTIN,IPATIN,PARTNE,IPATNE,
     &            MV,NP,IETA)
      CALL PINSEX(ID1,ID2,V,U,PARTIN,IPATIN,PARTNE,IPATNE,
     &            MV,NP,ISEX)
      IF(ISEX.NE.0.OR.IETA.NE.0) THEN
      ITYP=2
      NIN=0
      ELSE
C                            12.03.94
      PS1=0.750
C
      CALL HEINEN(PARTIN,IPATIN,PARTNE,IPATNE,MV,NP,NIN)
        if(ISLERR.eq.1)return   ! HADES1
c        if(ISLERR.eq.1)then
c          CONTINUE!WRITE(25,*)'     HEINEN'
c          return
c        end if
      ITYP=4
      DO  19  I=1,11
   19   IF(ITHEA(I).NE.0) TDIAG=DTYP(I)
      IF(NIN.NE.0)  GO  TO  11
      ENDIF
   27 CONTINUE
      IF(ITYP.EQ.4)  GO  TO  3
      DO  2  MM=1,NP
      M=MM
      IF(NP.EQ.2.AND.MM.EQ.2) M=3
      DO  1 K=1,9
      PA(K)=PME(K,MV+M)
      IF(K.LE.5)  IPA(K)=IME(K,MV+M)
    1 CONTINUE
      CALL PANUN(PA,IPA,IKK)
      CALL IDPANU(IDK,IKK,PNAK)
      IDPME(MV+M)=IDK
    2 CONTINUE
    3 CONTINUE
      IF(ITYP.EQ.0) GO  TO  603
      IF(NCAS.LT.NCPRI.AND.LOWMIS.EQ.0)  RETURN
c-Sob      IF(ITYP.NE.4)
c-Sob     *write(16,601) ATYP(ITYP),      PARTIN,IPATIN,PARTNE,IPATNE
c-Sob      IF(ITYP.EQ.4.AND.IHELE.EQ.1)
c-Sob     *write(16,611) ATYP(ITYP),TDIAG,PARTIN,IPATIN,PARTNE,IPATNE
c-Sob      IF(ITYP.EQ.4.AND.IHELE.EQ.2)
c-Sob     *write(16,612) ATYP(ITYP),TDIAG,PARTIN,IPATIN,PARTNE,IPATNE
c-Sob      DO  600  K=1,NP
c-Sob      M=MV+K
c-Sob      IF(NP.EQ.2.AND.K.EQ.2.AND.(IABS(IPATIN(4))+IABS(IPATNE(4))).GT.0)
c-Sob     *M=M+1
c-Sob      write(16,602) M,(PME(I,M),I=4,6),PME(9,M),(IME(J,M),J=1,5)
c-Sob  600 CONTINUE
      RETURN
c-Sob  601 FORMAT(1X,A4,9(1X,F10.3),4I2,I10/
c-Sob     *          5X,9(1X,F10.3),4I2,I10)
c-Sob  611 FORMAT(1X,A4,1X,A8,'(LE)'/
c-Sob     *          5X,9(1X,F10.3),4I2,I10/
c-Sob     *          5X,9(1X,F10.3),4I2,I10)
c-Sob  612 FORMAT(1X,A4,1X,A8,'(HE)'/
c-Sob     *          5X,9(1X,F10.3),4I2,I10/
c-Sob     *          5X,9(1X,F10.3),4I2,I10)
c-Sob  602 FORMAT(1X,I3,4(1X,F10.3),4I2,I10)
  603 CONTINUE
c-Sob      write(16,604)
c-Sob  604 FORMAT(2X,'TYPNEW: ITYP=0 !!!!!'/)
      RETURN
       END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE  PINSEX(ID1,ID2,V,U,PIN,IIN,PN,IPN,MV,NP,ISEX)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 V(3),U,PIN(9),PN(9),P1(3),P2(3),MPI,MN
      LOGICAL YESELA
      COMMON/YESELA/ YESELA
      COMMON/STREXC/STREXC
      COMMON/MEMAGT/PME(9,5999),IME(5,5999)
      DIMENSION  IIN(5),IPN(5)
      DATA  MPI/0.140/,MN/0.940/
      ISEX=0
      STREXC=0.
      IE1=IIN(1)
      IS1=IIN(3)
      IB1=IIN(4)
      IE2=IPN(1)
      IS2=IPN(3)
      IB2=IPN(4)
      IF((IB1+IB2).ne.1.or.IS1.ne.0.or.IS2.ne.0) RETURN
      IF((IE1+IE2).gt.2.or.(IE1+IE2).lt.0)      RETURN
*      IF(PIN(9).gt.1.d0.or.PN(9).gt.1.d0)        RETURN
c      IF(IB1.EQ.0.AND.IS1.EQ.0.AND.(ID2.EQ.1120.OR.ID2.EQ.1220))
c     *   GO  TO  10
c      RETURN
c   10 CONTINUE
c      IF(ID1.EQ.120.OR.ID1.EQ.-120.OR.ID1.EQ.110) THEN
c        TPIN=(U**2-(PIN(9)+PN(9))**2)/(2.*PN(9))
c      ELSE
      TPIN=(U**2-(MPI+MN)**2)/2./MN
c      ENDIF
      IF(TPIN.LE.0.759)           RETURN
      SSP=SST(0,1,IE1,IE2,4,TPIN,0)
      SSM=SST(0,1,IE1,IE2,5,TPIN,0)
      SS0=SST(0,1,IE1,IE2,6,TPIN,0)
      SL =SST(0,1,IE1,IE2,3,TPIN,0)
      SSEX=SSP+SSM+SS0+SL
      IF(SSEX.LE.1.D-10)            RETURN
      STREXC=SSEX
      PX1=PIN(4)
      PY1=PIN(5)
      PZ1=PIN(6)
      AM1=PIN(9)
      PX2=PN(4)
      PY2=PN(5)
      PZ2=PN(6)
      AM2=PN(9)
      CALL  CROSEC(1,ID1,ID2,PX1,PY1,PZ1,AM1,PX2,PY2,PZ2,AM2,SITO,0)
      SINO=SITO
      IF(YESELA)         GO TO 1
      CALL  CROSEC(0,ID1,ID2,PX1,PY1,PZ1,AM1,PX2,PY2,PZ2,AM2,SIEL,0)
      CALL  CROSEC(2,ID1,ID2,PX1,PY1,PZ1,AM1,PX2,PY2,PZ2,AM2,SIEX,0)
      SINO=SITO-SIEL-SIEX
    1 IF(SINO.LE.0.)                RETURN
C
c      WRITE(16,100) TPIN,SL,SS0,SSM,SSP,SSEX,SINO
c  100 FORMAT(2X,'TPIN,SL,SS0,SSM,SSP,SSEX,SINO='/1X,7(1PE10.3,1X))
C
      RNS=RNDMD(-1)
c      WRITE(16,'(5x,''U,ID1,ID2,RNS='',5x,F8.3,1x,2I6,1x,F8.3)')
c     &  U,ID1,ID2,RNS
      IF(RNS.GT.(SSEX/SINO))  RETURN
      RN1=RNDMD(-1)
      IF(RN1.LE.(SL/SSEX))    THEN
C     K L
      AMY=1.116
      IEY= 0
      ELSE IF(RN1.LE.((SL+SS0)/SSEX)) THEN
C     K S0
      AMY=1.189
      IEY= 0
      ELSE IF(RN1.LE.((SL+SS0+SSM)/SSEX)) THEN
C     K S-
      AMY=1.189
      IEY=-1
      ELSE
C     K S+
      AMY=1.189
      IEY=1
      ENDIF
      AMK=0.494
      IF(U.LE.(AMK+AMY+0.001))       RETURN
      if(PN(9).gt.1..or.PIN(9).gt.1.)   AMK=0.896
      IF(U.LE.(AMK+AMY+0.001))       AMK=0.494
      IEK=IE1+IE2-IEY
c-Sob      IF(IEK.LT.0)  THEN
c-Sob        WRITE(16,*) '  PINSEX: K- ???'
c-Sob        WRITE( *,*) '  PINSEX: K- ???'
c-Sob      ENDIF

c      WRITE(16,'(1x,''U,AMK,AMY,ID1,ID2='',3(F8.3,1x),2I6))')
c    & U,AMK,AMY,ID1,ID2
      TMAX=3.5344*TPIN*(TPIN+MPI)/(1.88*(TPIN+MPI)+MPI*MPI+0.8836)
      BET=BHN(TPIN,MPI)
      EXM=BET*TMAX
      IF(EXM.GT.30.)    EXM=30.
      CT=1.+(2.*DLOG(1.+RNDMD(-1)*(DEXP(-EXM)-1.)))/EXM
      FI=6.283185*RNDMD(-1)
      CALL  ABEL(PIN,V,U,P1,P2,CT,FI,AMK,AMY)
      PME(1,MV+3)=0.
      PME(2,MV+3)=0.
      PME(3,MV+3)=0.
      PME(4,MV+3)=P1(1)
      PME(5,MV+3)=P1(2)
      PME(6,MV+3)=P1(3)
      PME(7,MV+3)=0.
      PME(9,MV+3)=AMK
      IME(1,MV+3)=IEK
      IME(2,MV+3)=0
      IME(3,MV+3)=1
      IME(4,MV+3)=0
      IME(5,MV+3)=0
      if(AMK.gt.0.50)  IME(5,MV+3)=IDINTG(1000. *TAUN(12))
*      if(AMK.gt.0.50)  IME(5,MV+3)=IDINTG(1000. *TAUN(14)) !for K*0
      PME(1,MV+1)=0.
      PME(2,MV+1)=0.
      PME(3,MV+1)=0.
      PME(4,MV+1)=-P1(1)
      PME(5,MV+1)=-P1(2)
      PME(6,MV+1)=-P1(3)
      PME(7,MV+1)=0.
      PME(9,MV+1)=AMY
      IME(1,MV+1)=IEY
      IME(2,MV+1)=0
      IME(3,MV+1)=-1
      IME(4,MV+1)=1
      IME(5,MV+1)=0
      ISEX=1
      NP = 2
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE  PINETA(ID1,ID2,V,U,PIN,IIN,PN,IPN,MV,NP,IETA)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 V(3),U,PIN(9),PN(9),P1(3),P2(3),MPI,MN,META
      LOGICAL YESELA
      COMMON/YESELA/ YESELA
      COMMON/SPINEC/SPINEC
      COMMON/MEMAGT/PME(9,5999),IME(5,5999)
      DIMENSION  IIN(5),IPN(5)
      DATA  MPI/0.140/,MN/0.940/,META/0.549/
      IETA=0
      SPINEC=0.
      IF((ID1.EQ.120.OR.ID1.EQ.-120.OR.ID1.EQ.110).AND.
     &   (ID2.EQ.1120.OR.ID2.EQ.1220)) go  to  10
      RETURN
   10 CONTINUE
      IE1=IIN(1)
      IE2=IPN(1)
      IF((IE1+IE2).EQ.2.OR.(IE1+IE2).EQ.-1)     RETURN
      TPIN=(U**2-(PIN(9)+PN(9))**2)/(2.*PN(9))
      SPINEC=SPINET(TPIN)
      IF(SPINEC.LE.1.D-10)                      RETURN
      PX1=PIN(4)
      PY1=PIN(5)
      PZ1=PIN(6)
      AM1=PIN(9)
      PX2=PN(4)
      PY2=PN(5)
      PZ2=PN(6)
      AM2=PN(9)
      CALL  CROSEC(1,ID1,ID2,PX1,PY1,PZ1,AM1,PX2,PY2,PZ2,AM2,SITO,0)
      SINO=SITO
      IF(YESELA)         GO TO 1
      CALL  CROSEC(0,ID1,ID2,PX1,PY1,PZ1,AM1,PX2,PY2,PZ2,AM2,SIEL,0)
      CALL  CROSEC(2,ID1,ID2,PX1,PY1,PZ1,AM1,PX2,PY2,PZ2,AM2,SIEX,0)
      SINO=SITO-SIEL-SIEX
    1 IF(SINO.LE.0.)                RETURN
C
c      WRITE(16,100) TPIN,SPINEC,SINO
c  100 FORMAT(2X,'TPIN,SPINEC,SINO='/1X,3(1PE10.3,1X))
C
      RNS=RNDMD(-1)
      IF(RNS.GT.(SPINEC/SINO))      RETURN
      TMAX=3.5344*TPIN*(TPIN+MPI)/(1.88*(TPIN+MPI)+MPI*MPI+0.8836)
      BET=BHN(TPIN,MPI)
      EXM=BET*TMAX
      IF(EXM.GT.30.)    EXM=30.
      CT=1.+(2.*DLOG(1.+RNDMD(-1)*(DEXP(-EXM)-1.)))/EXM
      FI=6.283185*RNDMD(-1)
      CALL  ABEL(PIN,V,U,P1,P2,CT,FI,META,MN)
      PME(1,MV+3)=0.
      PME(2,MV+3)=0.
      PME(3,MV+3)=0.
      PME(4,MV+3)=P1(1)
      PME(5,MV+3)=P1(2)
      PME(6,MV+3)=P1(3)
      PME(7,MV+3)=0.
      PME(9,MV+3)=META
      IME(1,MV+3)=0
      IME(2,MV+3)=0
      IME(3,MV+3)=0
      IME(4,MV+3)=0
      IME(5,MV+3)=IDINTG(1000. *TAUN(8))
      PME(1,MV+1)=0.
      PME(2,MV+1)=0.
      PME(3,MV+1)=0.
      PME(4,MV+1)=-P1(1)
      PME(5,MV+1)=-P1(2)
      PME(6,MV+1)=-P1(3)
      PME(7,MV+1)=0.
      PME(9,MV+1)=MN
      IME(1,MV+1)=IE1+IE2
      IME(2,MV+1)=0
      IME(3,MV+1)=0
      IME(4,MV+1)=1
      IME(5,MV+1)=0
      IETA=1
      NP = 2
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE ELEX(V,U,TIN1,PARTIN,IPATIN,IPATNE,
     *MV,NP,L,MS,MQ,KSI,ME)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 V,U,TIN1,PARTIN
C     BLOCK OF CALCULATION OF PARTICLE CHARACTERISTICS IN ELASTIC AND
C     CHARGE EXCHANGE SCATTERING.
      DIMENSION PMEMO(9,5999),IMEMO(5,5999),V(3),PARTIN(9),IPATIN(5),
     1IPATNE(5),PIST(3),PNST(3)
      COMMON/MEMAGT/PMEMO,IMEMO
      COMMON/ISOB2/ISOB2
      IS=IPATIN(3)
      NS=IPATNE(3)
      IF (IPATIN(2)) 10,11,10
   10 CMI = 0.14
      CMN = 0.94
                        GO TO 12
   11 CMI = PARTIN(9)
      CMN = 0.94
                        GO TO 12
   12 SEX=CROSEG(L,MS,MQ,KSI,2,TIN1,PARTIN(9),IPATIN(5))
      SEL=CROSEG(L,MS,MQ,KSI,1,TIN1,PARTIN(9),IPATIN(5))
      BETAEX=SEX/(SEX+SEL)
      DRND=RNDMD(-1)
      IF(DRND-BETAEX)14,13,13
   13 IE = IPATIN(1)
      NE = IPATNE(1)
      IF(ISOB2.EQ.0)  GO  TO  113
      ISO=0
      IF(IPATIN(4).EQ.0.AND.IPATIN(3).EQ.0.AND.IPATIN(5).EQ.0)
     *CALL  PINDEL(TIN1,U,MV,IE,NE,NP,ISO)
      IF(ISO.EQ.1)  RETURN
  113 CTSTI = COSEL(L,MQ,KSI,TIN1,PARTIN(9))
                                     GO TO 18
   14 IF(IPATIN(5).NE.0.AND.IPATIN(4).EQ.1)    GO  TO  117
      IF(IPATIN(3).EQ.0)  GO  TO  114
      IF(IPATIN(3).NE.0)
     *CALL STREX(IPATIN(3),IPATIN(1),IPATNE(1),CMI,CMN,IS,NS,IE,NE)
      GO  TO  17
  114 IF (IPATIN(1)) 15,16,15
   15 IE = 0
      NE = ME
                   GO TO 17
   16 NE = 1-IPATNE(1)
      IE = ME-NE
                   GO TO 17
   17 CTSTI = COSEX (L,TIN1,PARTIN(9))
                   GO TO 18
  117 NE=1-IPATNE(1)
      IE=ME-NE
      GO  TO  113
   18 FISTI = 0.
      CALL ABEL(PARTIN,V,U,PIST,PNST,CTSTI,FISTI,CMI,CMN)
      IF(MV-5997) 19,19,20
   20 NP = 0
c-Sob      write(16,21)
c-Sob   21 FORMAT (45X,29H MEMAGT IS EXCEEDED IN CASCAD)
      RETURN
   19 CONTINUE
      PMEMO(1,MV+3)=0.
      PMEMO(2,MV+3)=0.
      PMEMO(3,MV+3)=0.
      PMEMO(4,MV+3) = PIST(1)
      PMEMO(5,MV+3) = PIST(2)
      PMEMO(6,MV+3) = PIST(3)
      PMEMO(7,MV+3)=0.
      PMEMO(9,MV+3)=CMI
      IMEMO(1,MV+3)=IE
      IMEMO(2,MV+3)=0
      IMEMO(3,MV+3) = IS
      IMEMO(4,MV+3) = IPATIN(4)
      IMEMO(5,MV+3) = IPATIN(5)
      PMEMO(1,MV+1)=0.
      PMEMO(2,MV+1)=0.
      PMEMO(3,MV+1)=0.
      PMEMO(4,MV+1) = PNST(1)
      PMEMO(5,MV+1) = PNST(2)
      PMEMO(6,MV+1) = PNST(3)
      PMEMO(7,MV+1)=0.
      PMEMO(9,MV+1)=CMN
      IMEMO(1,MV+1)=NE
      IMEMO(2,MV+1)=0
      IMEMO(3,MV+1) = NS
      IMEMO(4,MV+1) = 1
      IMEMO(5,MV+1) = IPATNE(5)
      NP = 2
      RETURN
       END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE STREX(IS1,IE1,NE1,CMI,CMN,IS,NS,IE,NE)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 CMI,CMN
      IF(IS1.LT.0)  GO  TO  10
      IE=1-IE1
    9 IS=IS1
      NS=0
      CMI=0.492
      CMN=0.940
      GO  TO  23
   10 IF((IE1+NE1).EQ.0)  GO  TO  17
      RND4=4.*RNDMD(-1)
      NBR=IDINTG(RND4)+1
      GO  TO  (11,12,13,14,14),NBR
   11 IE=1
      GO  TO  15
   12 IE=0
      GO  TO  16
   13 IE=1
      GO  TO  16
   14 IE=-1
      GO  TO   9
   15 IS=0
      NS=-1
      CMI=0.140
      CMN=1.116
      GO  TO  23
   16 IS=0
      NS=-1
      CMI=0.140
      CMN=1.189
      GO  TO  23
   17 RND5=5.*RNDMD(-1)
      NBR=IDINTG(RND5)+1
      GO  TO  (18,19,20,21,22,22),NBR
   18 IE=0
      GO  TO  15
   19 IE=0
      GO  TO  16
   20 IE=-1
      GO  TO  16
   21 IE=1
      GO  TO  16
   22 IF(IE1.EQ. 0)  IE=-1
      IF(IE1.EQ.-1)  IE= 0
      GO  TO  9
   23 NE=IE1+NE1-IE
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE PINDEL(T1,DM,MV,IE,NE,NP,ISO)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 T1,DM
      COMMON/MEMAGT/PME(9,5999),IME(5,5999)
      ISO=0
      IF(T1.GT.0.5)                  RETURN
      IF(DM.LE.1.082)                RETURN
      CALL  WMD(DM,TAU0,FMD)
      DRND=RNDMD(-1)
      IF(DRND.GT.FMD) RETURN
      PME(1,MV+1)=0.
      PME(2,MV+1)=0.
      PME(3,MV+1)=0.
      PME(4,MV+1)=0.
      PME(5,MV+1)=0.
      PME(6,MV+1)=0.
      PME(7,MV+1)=0.
      PME(9,MV+1)=DM
      IME(1,MV+1)=IE+NE
      IME(2,MV+1)=0
      IME(3,MV+1)=0
      IME(4,MV+1)=1
      IME(5,MV+1)=IDINTG(1000.*TAU0)
      IF(IME(5,MV+1).EQ.0)  IME(5,MV+1)=1
      NP=1
      ISO=1
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE  TMATUR(NP,V,M,P1,P2,MV,PS,TL,TMAT)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 TMAT,V,P1,P2,PS,TL
      COMMON/MMATUR/MMES,MBAR
      COMMON/INTTYP/ ITYP
      COMMON/ITHEA/ITHEA(11)
       COMMON/MEMAGT/PMEMO(9,5999),IMEMO(5,5999)
      COMMON/NUCSP/VPR(3),VTA(3),RADP(3),RADT(3),VEV(3),VRE(3),GEV,GRE
     *,VEP(3),VET(3),GEP,GET
      COMMON/IACT/ IACT
      COMMON/CSLID/CLIDER(5999)
      COMMON/CMALI/CMALI
      COMMON/ACTIM/TINT
      DIMENSION V(3),P1(9),P2(9),PS(3),RL(3)
      REAL*8  MMAT,MMES,MBAR,MLID
C
C     ITYP=1(ABSORP),2(ELEX),3(BINEL),4(HEINEN)
C
      TMAT=0.
      IF(NP.LE.2.OR.ITYP.LE.2.or.
     &(ITYP.EQ.4.AND.ITHEA(4).EQ.1))  THEN
c   formation time TMAT=0 for  2-body channel
c                              absorption
c                              elastic scattering
c                    high energy collision with.......?
      CLIDER(M)=1.
      ELSE
      V2=V(1)**2+V(2)**2+V(3)**2
      U=(P1(8)+P1(9)+P2(8)+P2(9))*DSQRT(1.-V2)
c      if(U.lt.3.5)                     GO  TO  20
      X=DABS(PS(3))/(U/2.)
C
C  XXX
      IF(IACT.LE.2)                    GO  TO  10
      RTAU=TL-TINT
      TMAT=TL
      IF(CLIDER(M).LT.CMALI)           GO  TO  20
C     IF(IMEMO(4,M).EQ.0)              GO  TO  20
      TMAT=0.
      CLIDER(M)=1.
                               GO  TO  20
C  XXX
   10 CONTINUE
      MLID=MBAR
      IF(IMEMO(4,M).EQ.0)              MMAT=MMES
      IF(IMEMO(4,M).NE.0)              MMAT=MBAR
C  FOR ANTI-BARYON
      IF(IMEMO(4,M).LT.0)              MMAT=MMES
C
      IF(M.EQ.(MV+2))                  MMAT=MLID
      TAU0=1./(5.06*MMAT)
      GL=1.+PMEMO(8,M)/PMEMO(9,M)
      FX=FXTMAT(X)
      TAU0L=GL*TAU0*FX
      BRND=RNDMD(-1)
C     AF=-DLOG(BRND)
      AF=1.
      TMAT=TAU0L*AF
  20  CONTINUE
      if(TMAT.LT.0.) then
c-Sob        write(16,'(2X,''TMATUR:ITYP='',I2,2X,''M='',I5,2X,
c-Sob     &  ''IM(5)='',I15,2X,''T='',F8.3,2X,''MASS='',F8.3,2X,
c-Sob     &  ''X='',E13.6,2X,''TMAT='',E13.6) ')
c-Sob     &  ITYP,M,IMEMO(5,M),(PMEMO(K,M),K=8,9),X,TMAT
      TMAT=0.
      endif
      IF(ITYP.EQ.3)           CLIDER(M)=1.
      ENDIF
      CALL  HTFORL(M,TMAT)
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE HTFORL(M,TL)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 TL
       COMMON/MEMAGT/PMEMO(9,5999),IMEMO(5,5999)
      COMMON/CSLID/CLIDER(5999)
      COMMON /HFOL/ HFOM(20,5),HFOB(20,5),HFNM(3,5),HFNB(3,5),
     *HTMI,HTM,HDT
      CLI=CLIDER(M)
      ILI=0
      IF(CLI.LT.0.10)                  ILI=1
      IF(CLI.GT.0.10.AND.CLI.LT.0.40)  ILI=2
      IF(CLI.GT.0.40.AND.CLI.LT.0.60)  ILI=3
      IF(CLI.GT.0.60.AND.CLI.LT.0.80)  ILI=4
      IF(CLI.GT.0.80)                  ILI=5
      IF(ILI.EQ.0)                     RETURN
      IB=IMEMO(4,M)
      IF(TL.LT.0.)                     GO  TO  10
      T=HTMI
      TMI=10.**HTMI
      IF(TL.GT.TMI)   T=DLOG10(TL)
      W=1.
      IF(IB.EQ.0)
     *CALL DHIST3(T,HTMI,HTM,HDT,HFOM,20,5,W,ILI)
      IF(IB.NE.0)
     *CALL DHIST3(T,HTMI,HTM,HDT,HFOB,20,5,W,ILI)
      IF(IB.EQ.0) HFNM(2,ILI)=HFNM(2,ILI)+TL
      IF(IB.NE.0) HFNB(2,ILI)=HFNB(2,ILI)+TL
      IF(IB.EQ.0) HFNM(3,ILI)=HFNM(3,ILI)+1.
      IF(IB.NE.0) HFNB(3,ILI)=HFNB(3,ILI)+1.
      RETURN
   10 CONTINUE
      IF(IB.EQ.0) HFNM(1,ILI)=HFNM(1,ILI)+1.
      IF(IB.NE.0) HFNB(1,ILI)=HFNB(1,ILI)+1.
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE PRTFL1
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /HFOL/ HFOM(20,5),HFOB(20,5),HFNM(3,5),HFNB(3,5),
     *HTMI,HTM,HDT
      COMMON /TLIMIT/ TLIMIT
      DO 10 J=1,5
      HFNM(1,J)=0.
      HFNB(1,J)=0.
      HFNM(2,J)=0.
      HFNB(2,J)=0.
      HFNM(3,J)=0.
      HFNB(3,J)=0.
      DO 10 I=1,20
      HFOM(I,J)=0.
   10 HFOB(I,J)=0.
      HTMI=-1.
      HTM=DLOG10(TLIMIT)+HTMI
      HDT=(HTM-HTMI)/15.
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE PRTFL2
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /HFOL/ HFOM(20,5),HFOB(20,5),HFNM(3,5),HFNB(3,5),
     *HTMI,HTM,HDT
      DIMENSION ANM(5),ANB(5),HM(5),HB(5)
c-Sob      write(16,100)
      DO  11  K=1,5
      ANM(K)=HFOM(17,K)
      IF(ANM(K).NE.0.)  HFOM(19,K)=HFOM(19,K)/ANM(K)
      IF(HFNM(3,K).NE.0.)  HFNM(2,K)=HFNM(2,K)/HFNM(3,K)
   11 CONTINUE
      DO 13 I=1,15
      DO 12 K=1,5
      HM(K)=0.
      IF(ANM(K).NE.0.)  HM(K)=HFOM(I,K)/ANM(K)
   12 CONTINUE
      TLO=HDT*FLOAT(I)
      T=10.**TLO
c-Sob      write(16,101) TLO,T,HM
   13 CONTINUE
      DO 14 I=16,20
c-Sob      write(16,102)       (HFOM(I,K),K=1,5)
   14 CONTINUE
c-Sob      write(16,103)       (HFNM(1,K),K=1,5)
c-Sob      write(16,104)       (HFNM(2,K),K=1,5)
  100 FORMAT(/5X,'DISTRIBUTION OVER LAB. FORMATION TIME(MESONS)'/
     *4X,'LOG T',3X,'T(FM*C)',6X,'CL=0.0',8X,'CL=1/3',8X,
     *'CL=1/2',8X,'CL=2/3',8X,'CL=1.0'/)
  101 FORMAT(2X,2(F7.3,3X),5(E12.4,2X))
  102 FORMAT(22X,5(E12.4,2X))
  103 FORMAT(17X,'T<0',2X,5(E12.4,2X))
  104 FORMAT(17X,'<T>',2X,5(E12.4,2X))
c-Sob      write(16,200)
      DO  21  K=1,5
      ANB(K)=HFOB(17,K)
      IF(ANB(K).NE.0.)  HFOB(19,K)=HFOB(19,K)/ANB(K)
      IF(HFNB(3,K).NE.0.)  HFNB(2,K)=HFNB(2,K)/HFNB(3,K)
   21 CONTINUE
      DO 23 I=1,15
      DO 22 K=1,5
      HB(K)=0.
      IF(ANB(K).NE.0.)  HB(K)=HFOB(I,K)/ANB(K)
   22 CONTINUE
      TLO=HDT*FLOAT(I)
      T=10.**TLO
c-Sob      write(16,201) TLO,T,HB
   23 CONTINUE
c-Sob      DO 24 I=16,20
c-Sob      write(16,202)       (HFOB(I,K),K=1,5)
c-Sob   24 CONTINUE
c-Sob      write(16,203)       (HFNB(1,K),K=1,5)
c-Sob      write(16,204)       (HFNB(2,K),K=1,5)
  200 FORMAT(5X,'DISTRIBUTION OVER LAB. FORMATION TIME(BARYONS)'/
     *4X,'LOG T',3X,'T(FM*C)',6X,'CL=0.0',8X,'CL=1/3',8X,
     *'CL=1/2',8X,'CL=2/3',8X,'CL=1.0'/)
  201 FORMAT(2X,2(F7.3,3X),5(E12.4,2X))
  202 FORMAT(22X,5(E12.4,2X))
  203 FORMAT(17X,'T<0',2X,5(E12.4,2X))
  204 FORMAT(17X,'<T>',2X,5(E12.4,2X))
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE HTFORP(M)
      IMPLICIT REAL*8 (A-H,O-Z)
       COMMON/MEMAGT/PME(9,5999),IME(5,5999)
      COMMON/CSLID/CLIDER(5999)
      COMMON /HFOP/ HFOM(20,5),HFOB(20,5),HFNM(3,5),HFNB(3,5),
     *HTMI,HTM,HDT
      CLI=CLIDER(M)
      ILI=0
      IF(CLI.LT.0.10)                  ILI=1
      IF(CLI.GT.0.10.AND.CLI.LT.0.40)  ILI=2
      IF(CLI.GT.0.40.AND.CLI.LT.0.60)  ILI=3
      IF(CLI.GT.0.60.AND.CLI.LT.0.80)  ILI=4
      IF(CLI.GT.0.80)                  ILI=5
      IF(ILI.EQ.0)                     RETURN
      TP2=PME(7,M)**2-PME(1,M)**2-PME(2,M)**2-PME(3,M)**2
      TP=DSQRT(DABS(TP2))
      IF(TP.LT.0.)                     GO  TO  10
      T=HTMI
      TMI=10.**HTMI
      IF(TP.GT.TMI)   T=DLOG10(TP)
      IB=IME(4,M)
      W=1.
      IF(IB.EQ.0)
     *CALL DHIST3(T,HTMI,HTM,HDT,HFOM,20,5,W,ILI)
      IF(IB.NE.0)
     *CALL DHIST3(T,HTMI,HTM,HDT,HFOB,20,5,W,ILI)
      IF(IB.EQ.0) HFNM(2,ILI)=HFNM(2,ILI)+TP
      IF(IB.NE.0) HFNB(2,ILI)=HFNB(2,ILI)+TP
      IF(IB.EQ.0) HFNM(3,ILI)=HFNM(3,ILI)+1.
      IF(IB.NE.0) HFNB(3,ILI)=HFNB(3,ILI)+1.
      RETURN
   10 CONTINUE
      IF(IB.EQ.0) HFNM(1,ILI)=HFNM(1,ILI)+1.
      IF(IB.NE.0) HFNB(1,ILI)=HFNB(1,ILI)+1.
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE PRTFP1
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /HFOP/ HFOM(20,5),HFOB(20,5),HFNM(3,5),HFNB(3,5),
     *HTMI,HTM,HDT
      COMMON /TLIMIT/ TLIMIT
      DO 10 J=1,5
      HFNM(1,J)=0.
      HFNB(1,J)=0.
      HFNM(2,J)=0.
      HFNB(2,J)=0.
      HFNM(3,J)=0.
      HFNB(3,J)=0.
      DO 10 I=1,20
      HFOM(I,J)=0.
   10 HFOB(I,J)=0.
      HTMI=-1.
      HTM=DLOG10(TLIMIT)+HTMI
      HDT=(HTM-HTMI)/15.
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE PRTFP2
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /HFOP/ HFOM(20,5),HFOB(20,5),HFNM(3,5),HFNB(3,5),
     *HTMI,HTM,HDT
      DIMENSION ANM(5),ANB(5),HM(5),HB(5)
c-Sob      write(16,100)
      DO  11  K=1,5
      ANM(K)=HFOM(17,K)
      IF(ANM(K).NE.0.)  HFOM(19,K)=HFOM(19,K)/ANM(K)
      IF(HFNM(3,K).NE.0.)  HFNM(2,K)=HFNM(2,K)/HFNM(3,K)
   11 CONTINUE
      DO 13 I=1,15
      DO 12 K=1,5
      HM(K)=0.
      IF(ANM(K).NE.0.)  HM(K)=HFOM(I,K)/ANM(K)
   12 CONTINUE
      TLO=HDT*FLOAT(I)
      T=10.**TLO
c-Sob      write(16,101) TLO,T,HM
   13 CONTINUE
c-Sob      DO 14 I=16,20
c-Sob      write(16,102)       (HFOM(I,K),K=1,5)
c-Sob   14 CONTINUE
c-Sob      write(16,103)       (HFNM(1,K),K=1,5)
c-Sob      write(16,104)       (HFNM(2,K),K=1,5)
  100 FORMAT(/5X,'DISTRIBUTION OVER PROPER FORMATION TIME(MESONS)'/
     *4X,'LOG T',3X,'T(FM*C)',6X,'CL=0.0',8X,'CL=1/3',8X,
     *'CL=1/2',8X,'CL=2/3',8X,'CL=1.0'/)
  101 FORMAT(2X,2(F7.3,3X),5(E12.4,2X))
  102 FORMAT(22X,5(E12.4,2X))
  103 FORMAT(17X,'T<0',2X,5(E12.4,2X))
  104 FORMAT(17X,'<T>',2X,5(E12.4,2X))
c-Sob      write(16,200)
      DO  21  K=1,5
      ANB(K)=HFOB(17,K)
      IF(ANB(K).NE.0.)  HFOB(19,K)=HFOB(19,K)/ANB(K)
      IF(HFNB(3,K).NE.0.)  HFNB(2,K)=HFNB(2,K)/HFNB(3,K)
   21 CONTINUE
      DO 23 I=1,15
      DO 22 K=1,5
      HB(K)=0.
      IF(ANB(K).NE.0.)  HB(K)=HFOB(I,K)/ANB(K)
   22 CONTINUE
      TLO=HDT*FLOAT(I)
      T=10.**TLO
c-Sob      write(16,201) TLO,T,HB
   23 CONTINUE
c-Sob      DO 24 I=16,20
c-Sob      write(16,202)       (HFOB(I,K),K=1,5)
c-Sob   24 CONTINUE
c-Sob      write(16,203)       (HFNB(1,K),K=1,5)
c-Sob      write(16,204)       (HFNB(2,K),K=1,5)
  200 FORMAT(5X,'DISTRIBUTION OVER PROPER FORMATION TIME(BARYONS)'/
     *4X,'LOG T',3X,'T(FM*C)',6X,'CL=0.0',8X,'CL=1/3',8X,
     *'CL=1/2',8X,'CL=2/3',8X,'CL=1.0'/)
  201 FORMAT(2X,2(F7.3,3X),5(E12.4,2X))
  202 FORMAT(22X,5(E12.4,2X))
  203 FORMAT(17X,'T<0',2X,5(E12.4,2X))
  204 FORMAT(17X,'<T>',2X,5(E12.4,2X))
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE TEVOLV(U,TL,IND)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 U,TL
      COMMON /HFT/ HFT(20,4),HTMI,HTM,HDT
      T=DLOG10(TL)
      W=1.
      CALL DHIST3(T,HTMI,HTM,HDT,HFT,20,4,W,1)
      CALL DHIST3(T,HTMI,HTM,HDT,HFT,20,4,U,2)
      IF(IND.EQ.1) GO TO 10
      CALL DHIST3(T,HTMI,HTM,HDT,HFT,20,4,W,3)
      CALL DHIST3(T,HTMI,HTM,HDT,HFT,20,4,U,4)
  10  CONTINUE
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE PRFMT1(NCAS)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /HFT/ HFT(20,4),HTMI,HTM,HDT
      COMMON /TLIMIT/ TLIMIT
      COMMON /N00/ N00
      DO 10 I=1,20
      DO 10 J=1,4
   10 HFT(I,J)=0.
      HTMI=0.
      HTM=DLOG10(TLIMIT)
      HDT=(HTM-HTMI)/15.
      N00=NCAS
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE PRFMT2(NCAS)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /HFT/ HFT(20,4),HTMI,HTM,HDT
      COMMON /TLIMIT/ TLIMIT
      COMMON /N00/ N00
c-Sob      write(16,110) NCAS,N00
c-Sob  110 FORMAT(10X,'NET NUMBER OF EVENTS NCAS-N00 =',I5,'-',I5/)
      IF((NCAS-N00).EQ.0)        RETURN
      AN0=1./FLOAT(NCAS-N00)
      IF(HFT(17,1).LE.0.)  RETURN
      AN1=1./HFT(17,1)
c-Sob      write(16,100)
      DO 20 I=1,15
      TLO=HDT*FLOAT(I)
      T=10.**TLO
      X1=0.
      Y1=0.
      X2=0.
      Y2=0.
      IF(HFT(I,1).LE.0.)   GO  TO  25
      Y1      =HFT(I,2)/HFT(I,1)
      X1      =HFT(I,1)*AN1
      IF(HFT(17,3).LE.0.)  GO TO 25
      IF(HFT(I,3).LE.0.)   GO  TO  25
      Y2      =HFT(I,4)/HFT(I,3)
      X2      =HFT(I,1)/HFT(17,3)
  25  CONTINUE
c-Sob      write(16,101) TLO,T,X1,Y1,X2,Y2
   20 CONTINUE
c-Sob      DO 21 I=16,20
c-Sob   21 write(16,102)       (HFT(I,K),K=1,4)
  100 FORMAT(5X,'DISTRIBUTION OVER INTERACTION TIME'/
     *4X,'LOG T',3X,'T(FM*C)',8X,'NTOT',10X,'UTOT',10X,
     *'NHCC',10X,'UHCC'/)
  101 FORMAT(2X,2(F7.3,3X),4(E12.4,2X))
  102 FORMAT(22X,4(E12.4,2X))
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE KINEMR(V,M,RL,TL)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 V,RL,TL
      COMMON/MEMAGT/PME(9,5999),IME(5,5999)
      DIMENSION V(3),RL(3)
C PROPER TIME
      TP2=PME(7,M)**2-PME(1,M)**2-PME(2,M)**2-PME(3,M)**2
      TP=DSQRT(DABS(TP2))
C TIME IN OBSERVER SYSTEM
      G=(PME(8,M)+PME(9,M))/PME(9,M)
      TL=G*TP
      RL(1)=0.
      RL(2)=0.
      RL(3)=0.
      CALL  HTFORP(M)
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      FUNCTION FXTMAT(X)
      REAL*8 FXTMAT,X
C     FXTMAT=X*(1.-X)
      FXTMAT=1.
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE VMNSP (PARTIN,IPATIN,U,MV,NP,ITH,MQ,TIN1,LP)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 PARTIN,U,TIN1
C     BLOCK OF CALCULATION OF SECONDARY PARTICLES NUMBER AND
C     DETERMINING ABSOLUTE VALUES OF MOMENTUMS IN INELASTIC INTERACTION.
      DIMENSION PMEMO(9,5999),IMEMO(5,5999),PARTIN(9),IPATIN(5)
      COMMON/MEMAGT/PMEMO,IMEMO
      KH = 0
      LP = 0
      AM3=PARTIN(9)
      IF(IPATIN(5).NE.0.AND.IPATIN(4).EQ.1)  AM3=0.940
   10 U1 = U
      LAMBDA = 1
   11 LTEMP = MV+LAMBDA
      IF(LTEMP-5999) 37,37,38
   38 NP = 0
c-Sob      write(16,39)
c-Sob   39 FORMAT (25X,'MEMAGT IS EXCEEDED IN CASCAN')
      RETURN
   37 CONTINUE
      IF (LAMBDA-1) 12,12,13
   12 PMEMO(9,MV+1)=0.94
      IMEMO(2,MV+1)=0
      IMEMO(3,MV+1)=0
      IMEMO(4,MV+1) = 1
      IMEMO(5,MV+1) = 0
                    GO TO 16
   13 IF (LAMBDA-3) 15,14,15
   14 PMEMO(9,MV+3) = AM3
      IMEMO(2,MV+3) = 0
      IMEMO(3,MV+3) = IPATIN(3)
      IMEMO(4,MV+3) = IPATIN(4)
      IMEMO(5,MV+3) = 0
                          GO TO 16
   15 PMEMO(9,LTEMP) = 0.14
      IMEMO(2,LTEMP) = 0
      IMEMO(3,LTEMP) = 0
      IMEMO(4,LTEMP) = 0
      IMEMO(5,LTEMP) = 0
                                    GO TO 16
   16 JB = JTYPB(ITH,MQ,LAMBDA)
      PMEMO(8,LTEMP) = PMOM(JB,TIN1)
      EL = DSQRT (PMEMO(8,LTEMP)**2+PMEMO(9,LTEMP)**2)
      DELTU = U1-EL
      IF (LAMBDA-2) 23,17,23
   17 IF (DELTU-AM3)       35,35,18
   18 IF (ITH) 19,22,19
   19 PMEMO(8,MV+3) = DSQRT (DELTU**2-AM3**2)
      PMEMO(9,MV+3) = AM3
      IMEMO(2,MV+3) = 0
      IMEMO(3,MV+3) = IPATIN(3)
      IMEMO(4,MV+3) = IPATIN(4)
      IMEMO(5,MV+3) = 0
      IF (PMEMO(8,MV+1)-PMEMO(8,MV+2)-PMEMO(8,MV+3)) 20,20,35
   20 IF (PMEMO(8,MV+1)-DABS(PMEMO(8,MV+2)-PMEMO(8,MV+3))) 35,35,21
   21 NP = 3
      RETURN
   22 U1 = DELTU
      LAMBDA = LAMBDA+1
                               GO TO 11
   23 IF (DELTU-0.14) 24,24,22
   24 IF (LAMBDA-1) 35,35,25
   25 IF (LAMBDA-3) 26,35,26
   26 EL=DELTU+EL
      PMEMO(8,LTEMP) = DSQRT (EL**2-PMEMO(9,LTEMP)**2)
      NP = LAMBDA
      I = 1
   28 ITEMP = MV+I
      C = PMEMO(8,ITEMP)
   29 IF (NP-I) 34,34,30
   30 IF (C-PMEMO(8,ITEMP+1)) 32,31,31
   31 I=I+1
      ITEMP=ITEMP+1
            GO  TO  29
   32 I = I+1
            GO TO 28
   34 PMAX = C
      SIGMA = 0.
      DO 33 I=1,NP
      ITEMP = MV+I
      SIGMA = SIGMA+PMEMO(8,ITEMP)
   33 CONTINUE
      IF (2.*PMAX-SIGMA) 27,35,35
   27 RETURN
   35 KH = KH+1
      IF (KH-100) 10,36,36
   36 LP = 2
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE CHINEL (IPATIN,L,MS,MQ,KSI,NP,MV,TIN1,ME,IPATNE,AMIN)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 TIN1,AMIN
C     DETERMINING OF SECONDARY PARTICLES CHARGES IN INELASTIC
C     SCATTERING.
      DIMENSION IPATIN(5),IPATNE(5),PMEMO(9,5999),IMEMO(5,5999)
      COMMON/MEMAGT/PMEMO,IMEMO
      IF(IPATIN(5).NE.0.AND.IPATIN(4).EQ.1)  GO  TO  41
      IF (NP-3) 21,10,21
   10 SPI0 = CROSEG (L,MS,MQ,KSI,4,TIN1,AMIN,IPATIN(5))
      STH = CROSEG (L,MS,MQ,KSI,7,TIN1,AMIN,IPATIN(5))
      BPI0 = SPI0/STH
      BPIEX = (SPI0+CROSEG(L,MS,MQ,KSI,5,TIN1,AMIN,IPATIN(5)))/STH
      TEMP1 = RNDMD(-1)
      IF (TEMP1-BPI0) 19,11,11
   11 IF (TEMP1-BPIEX) 20,12,12
   12 IMEMO(1,MV+1)=IPATNE(1)-(IPATIN(4)-1)*IPATIN(1)
      IMEMO(1,MV+3)=(IPATIN(4)-1)*IPATIN(1)**2-IPATIN(4)*IPATIN(1)+1
      GO TO 18
   18 IMEMO(1,MV+2)=ME-IMEMO(1,MV+1)-IMEMO(1,MV+3)
                                         RETURN
   19 IMEMO(1,MV+1) = IPATNE(1)
                          IMEMO(1,MV+2) = 0
      IMEMO(1,MV+3) = IPATIN(1)
                          RETURN
   20 IMEMO(1,MV+1) = 1-IPATNE(1)
                            IMEMO(1,MV+3) = IPATIN(1)
      IMEMO(1,MV+2)=ME-IMEMO(1,MV+1)-IMEMO(1,MV+3)
                                         RETURN
   21 IF (RNDMD(-1) -0.5) 22,22,23
   22 IMEMO(1,MV+1) = 1
                    GO TO 24
   23 IMEMO(1,MV+1) = 0
                    GO TO 24
   24 IF (MQ-1) 28,28, 25
   25 IF (RNDMD(-1) -0.5) 26,26,27
   26 IMEMO(1,MV+3) = 1
                    GO TO 28
   27 IMEMO(1,MV+3) = 0
                    GO TO 28
   28 LAMBDA = 2
   29 IF (MQ-1) 31,31,30
   30 IF (LAMBDA-3) 31,36,31
   31 TEMP2 = RNDMD(-1)
      MTEMP = MV+LAMBDA
      IF (TEMP2-1./3.) 33,32,32
   32 IF (TEMP2-2./3.) 34,35,35
   33 IMEMO(1,MTEMP) = 1
                     GO TO 36
   34 IMEMO(1,MTEMP) = 0
                     GO TO 36
   35 IMEMO(1,MTEMP) = -1
                     GO TO 36
   36 IF (LAMBDA-NP) 37,38,38
   37 LAMBDA = LAMBDA+1
                    GO TO 29
   38 SIGQ = 0.
      DO 39 I=1,NP
      ITEMP = MV+I
      SIGQ = SIGQ+IMEMO(1,ITEMP)
   39 CONTINUE
      IF (ME-SIGQ) 21,40,21
   40 RETURN
   41 IE1=1
      IF(RNDMD(-1).GT.0.5)  IE1=0
      BR=RNDMD(-1)
      IF(BR.LE.0.33333333)                      IE2=-1
      IF(BR.GT.0.33333333.AND.BR.LE.0.66666667) IE2= 0
      IF(BR.GT.0.66666667)                      IE2= 1
      IE3=ME-IE1-IE2
      IF(IE3.LT.0.OR.IE3.GT.1)  GO  TO  41
      IMEMO(1,MV+1)=IE1
      IMEMO(1,MV+2)=IE2
      IMEMO(1,MV+3)=IE3
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      FUNCTION  POTBAR(NU,X,Y,Z)
*  calculated as U=U0*rho(r)/rho(0) with U0=200 MeV according to the paper
*  by Ye.S.Golubeva, A.S.Iljinov at al. Nucl.Phys. A483 (1988) 539.
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 POTBAR,X,Y,Z
      COMMON/HCASC/ANUCL1,ANUCL2,ZNUCL1,ZNUCL2,T0,EPS1,EPS2,VPI,A1,A2,
     *C1,C2,D1,D2,R0N1,R0N2,TF01,TF02,RM1,RM2
      DATA PBAR0 /0.200/
      POTBAR=0.
      if(NU.eq.1) then
      A=A1
      C=C1
      D=D1
      AN=ANUCL1
      RMAX=RM1
      else
      A=A2
      C=C2
      D=D2
      AN=ANUCL2
      RMAX=RM2
      endif
      R=dsqrt(X**2+Y**2+Z**2)
      if(R.gt.1.5*RMAX)  RETURN
      if(AN.le.10.) then
      RORO0=dexp(-(R/A)**2)
      else
      RORO0=(1.+dexp(-A/C)) / (1.+dexp((R-A)/C))
      endif
      POTBAR=PBAR0*RORO0
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      FUNCTION  POTEN(P,IP,A,C,D,TF0,VPION,EPS)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 POTEN,P,A,C,D,TF0,VPION,EPS
      DIMENSION   IP(5),P(9)
      COMMON/HCASC/ANUCL1,ANUCL2,ZNUCL1,ZNUCL2,T0,EPS1,EPS2,VPI,A1,A2,
     *C1,C2,D1,D2,R0N1,R0N2,TF01,TF02,RM1,RM2
      COMMON/RESULT/AN1,AN2,ZN1,ZN2,ENEXT1,ENEXT2,PNUCL1(3),PNUCL2(3),
     *AMNUC1(3),AMNUC2(3)
      IF(IP(5).NE.0)  GO  TO  10
      IF(IP(3))  10,11,10
   10 POTEN=0.
                 RETURN
   11 IF(IP(4))   13,12,13
   12 POTEN = VPION
                 RETURN
   13 IF(DABS(A-A1)-.001) 100,101,101
  100 AN=ANUCL1
      RMAX=RM1
      NU=1
                 GO TO 105
  101 AN=ANUCL2
      RMAX=RM2
      NU=2
  105 CONTINUE
      R=DSQRT(P(1)**2+P(2)**2+P(3)**2)/RMAX
      R=RPOTEN(R)
      IF(R-1.5)  14,10,10
   14 R=R*RMAX
      IF(AN-10.) 106,106,107
  106 TF=TF0*DEXP(-(2./3.)*(R/A)**2)
                         GO TO 108
  107 TF=TF0*(((1.+DEXP(-A/C))/(1.+DEXP((R-A)/C)))**0.6666667)
  108 CONTINUE
c     !!! 20.06.1995
      IF(NU.EQ.1) TF=TF*(AN1/ANUCL1)**(2./3.)
      IF(NU.EQ.2) TF=TF*(AN2/ANUCL2)**(2./3.)
c
      POTEN=TF+EPS
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE RECUL(NU,PX,PY,PZ,X,Y,Z)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 PX,PY,PZ,X,Y,Z
      COMMON/RESULT/AN1,AN2,ZN1,ZN2,ENEXT1,ENEXT2,PNUCL1(3),PNUCL2(3),
     *AMNUC1(3),AMNUC2(3)
4     IF(NU-1) 10,10,11
   10 PNUCL1(1)=PNUCL1(1)+PX
      PNUCL1(2)=PNUCL1(2)+PY
      PNUCL1(3)=PNUCL1(3)+PZ
      AMNUC1(1)=AMNUC1(1)+Z*PY-Y*PZ
      AMNUC1(2)=AMNUC1(2)+X*PZ-Z*PX
      AMNUC1(3)=AMNUC1(3)+Y*PX-X*PY
                              GO TO 12
   11 PNUCL2(1)=PNUCL2(1)+PX
      PNUCL2(2)=PNUCL2(2)+PY
      PNUCL2(3)=PNUCL2(3)+PZ
      AMNUC2(1)=AMNUC2(1)+Z*PY-Y*PZ
      AMNUC2(2)=AMNUC2(2)+X*PZ-Z*PX
      AMNUC2(3)=AMNUC2(3)+Y*PX-X*PY
   12 CONTINUE
      RETURN
             END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      FUNCTION  RPOTEN(R)
      REAL*8 RPOTEN,R
      RPOTEN=R
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE  COLDEV(R0X,R0Y,R0Z,T0,ANUCL1,ANUCL2,ZNUCL1,ZNUCL2,
     *P01,R01,IPRO)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 R0X,R0Y,R0Z,T0,ANUCL1,ANUCL2,ZNUCL1,ZNUCL2,P01,R01
      DIMENSION  P01(3),R01(3)
      IF(ANUCL1-0.5)  7,7,8
    7 P0=DSQRT(T0*(T0+0.28))
      E0=T0+0.14
      VZ=P0/(E0+ANUCL2*0.940)
      G=1./DSQRT(1.-VZ*VZ)
      E0S=G*(E0-P0*VZ)
      P0S=DSQRT(E0S**2-0.0196)
      ANEFF=0.148936
                   GO TO 9
    8 P0=ANUCL1*DSQRT(T0*(T0+1.88))
      E0=ANUCL1*(T0+0.940)
      VZ=P0/(E0+ANUCL2*0.940)
      G=1./DSQRT(1.-VZ*VZ)
      E0S=G*(E0-P0*VZ)
      P0S=DSQRT(E0S**2-(ANUCL1*0.940)**2)
      ANEFF=ANUCL1
    9 B0=0.940*(ANEFF *ANUCL2/(ANEFF +ANUCL2))*ZNUCL1*ZNUCL2*
     *0.001440/P0S**2
      B=DSQRT(R0X**2+R0Y**2)
      R12=DSQRT(B*B+R0Z**2)
      G1=DSQRT(B*B+B0*B0)
      SB=B0/G1
      CB=B/G1
      CF=R0X/B
      SF=R0Y/B
      PX=P0S*SB*CF
      PY=P0S*SB*SF
      PZ=G*(P0S*CB+B0*VZ)
      B1=B0+G1
      IF(B1-R12)  10,10,11
   11 IPRO=1
             RETURN
   10 IPRO=0
      R0X=B1*CF
      R0Y=B1*SF
      R0Z=-DSQRT(R12**2-B1**2)
      PM=DSQRT(PX**2+PY**2+PZ**2)
      BS=(PX*R0X+PY*R0Y+PZ*R0Z)/PM
      AS=DSQRT(R12**2-BS**2)
      P01(1)=(R0Z-BS*PZ/PM)/AS
      P01(2)=(PX*R0Y-PY*R0X)/(AS*PM)
      P01(3)=PZ/PM
      R01(1)=(R0X-BS*PX/PM)/AS
      R01(2)=(PY*R0Z-PZ*R0Y)/(AS*PM)
      R01(3)=PX/PM
      RETURN
             END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE GHELP(R0N,ANUCL,A,C,D,TF0,RM)        ! 13.02.96
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 R0N,ANUCL,A,C,D,TF0,RM
      DIMENSION   W(8),FIKS(8)
      DATA W/0.1012285363,0.2223810345,0.3137066459,0.3626837834,
     *0.3626837834,0.3137066459,0.2223810345,0.1012285363/,
     * FIKS/0.9602898565,0.7966664774,0.5255324099,0.1834346425,
     *-0.1834346425,-0.5255324099,-0.7966664774,-0.9602898565/
      PI=3.141593
      IF(ANUCL.le.10.) then
      A=R0N
      RM=A*DSQRT(-DLOG(D))
      RO0=ANUCL/((DSQRT(PI)*A)**3)
      ELSE
      S=0.
      A=R0N*ANUCL**0.333333333
      RM=A+C*DLOG((1.-D)/D)
      DO  K=1,8
        SK=(((FIKS(K)+1.)**2)/(DEXP(RM*(FIKS(K)+1.)/(2.*C))+
     &    DEXP(A/C)))*W(K)
        S=S+SK
      ENDDO
      S=(RM**3)*S/8.
      RO0=ANUCL/(12.566370*(DEXP(A/C)+1.)*S)
      ENDIF
      TF0=0.1985*(RO0/2.)**0.666666667
      RH=A
      IF(ANUCL.LE.10.)  RH=A*DSQRT(DLOG(2.D0))
c     write(16,'(10X,''ANUCL = '',F6.0,''  R(1/2)=  '',F6.3,    ! 13.02.96
c     &''  RO0(NUCL)='',F6.4/)') ANUCL,RH,RO0
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      FUNCTION COSEX (L,T,CM)
c   edited by KKG, COSTA is replaced by COSTAN, May 2008 
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 COSEX,T,CM
      data hlf, one, two /0.51d0, 1.d0, 2.d0/
      data emn/0.939d0/ 
C     BLOCK OF COSINUS CALCULATION FOR CHARGE EXCHANGE SCATTERING.
c  kkg 10/28/03        include gamma +N --> pi+- + N1
      if(L.ne.0)              then
c  gamma + N => pi+-  + N1
        if(T.le.0.51)     then
          COSEX = COSTAN(24,T)
c         COSEX = cosgaml(24)
        elseif(T.le.1.0)  then 
          COSEX = COSTAN(25,T)
c         COSEX = cosgaml(25)
        else
          COSEX = COSTAN(26,T)
c         COSEX = cosgaml(26)
        endif
        return
      endif 
   11 IF (T-0.08) 15,15,16
   15 COSEX = COSTAN(12,T)
                            RETURN
   16 IF (T-0.3) 17,17,18
   17 COSEX = COSTAN(13,T)
                            RETURN
   18 IF (T-1.0) 19,19,20
   19 COSEX = COSTAN(10,T)
                            RETURN
   20 IF (T-2.4) 21,21,22
   21 COSEX = COSTAN(11,T)
                            RETURN
   22 continue
c  KKG 02/16/04 (exact tmax):
          tem = two*emn
          tmax = tem*T*tem*(T + two*CM)/(tem*T + (emn + CM)**2) 
          tb1 = 7.5d0*tmax
          r1=RNDMD(-1) 
          COSEX = one + (two*log(one + r1*
     &                  (exp(-tb1) - one)))/tb1
        RETURN
         END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      FUNCTION BHN(T,CM)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 BHN,T,CM
      IF(CM-0.9)  10,10,11
   10 BHN=11.
      IF(T.LT.10.)  BHN=7.5
      RETURN
   11 BHN=8.3+0.56*DLOG(1.88*(T+1.88))
      RETURN
            END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      FUNCTION COSEL (L,MQ,KSI,T,CM)
      IMPLICIT REAL*8 (A-H,O-Z)
c   edited by KKG, COSTA is replaced by COSTAN, May 2008 
      REAL*8 COSEL,T,CM
C     BLOCK OF COSINUS CALCULATION FOR ELASTIC SCATTERING.
c  kkg 10/28/03        include gamma +N --> pi0 + N
      data hlf, one, two /0.51d0, 1.d0, 2.d0/
c
   13 IF(MQ-2) 24,14,78
c  N-He4 scattering
   78 COSEL=COSAL(T)
                       RETURN
c  N-N scattering
   14 GO  TO  (17,15,17,15,15,17),KSI
c   n + p:
   15 IF (T-0.97) 16,16,19
   16 COSEL = COSTAN(3,T)
                           RETURN
c   n + n or p + p:
   17 IF (T-0.46) 18,18,19
   18 COSEL=1.-2.*RNDMD(-1)
                             RETURN
c   All types:
   19 IF (T-2.8) 20,20,21
   20 COSEL = (1.+COSTAN(1,T))/2.
                                   RETURN
   21 IF (T-10.) 22,22,23
   22 COSEL = (3.+COSTAN(2,T))/4.
                                   RETURN
c  KKG 02/16/04:
   23 tmax = two*T*CM
      tm8 = 8.7d0*tmax
      r1 = RNDMD(-1) 
      COSEL = one + (two*log(one + r1*(exp(-tm8) - one)))/tm8
      RETURN
   24 if(L.ne.0)  go  to  42
c  pi-N scattering
      IF (KSI-2) 25,34,33
c  pi+ p or pi- n scattering:
   25 IF (T-0.08) 26,26,27
   26 COSEL = COSTAN(4,T)
                           RETURN
   27 IF (T-0.3) 28,28,29
   28 COSEL = COSTAN(5,T)
                           RETURN
   29 IF (T-1.) 30,30,31
   30 COSEL = COSTAN(6,T)
                           RETURN
   31 IF (T-2.4) 32,32,23
   32 COSEL = COSTAN(7,T)
                           RETURN
c  pi0 p or pi0 n scattering:
   33 IF(RNDMD(-1)-0.5) 25,25,34
c  pi+ n or pi- p scattering:
   34 IF (T-0.08) 35,35,36
   35 COSEL = COSTAN(8,T)
                           RETURN
   36 IF (T-0.3) 37,37,38
   37 COSEL = COSTAN(9,T)
                           RETURN
   38 IF (T-1.0) 39,39,40
   39 COSEL = COSTAN(10,T)
                            RETURN
   40 IF (T-2.4) 41,41,23
   41 COSEL = COSTAN(11,T)
                            RETURN
   42 continue
c  gamma + N => pi0 + N
      if(T.le.0.45)  then
c       COSEL = COSTAN (22,T)
c       COSEL = cosgaml(22)
      else
        COSEL = COSTAN (23,T)
c       COSEL = cosgaml(23)
      endif
      RETURN
         END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c      double precision function cduarteq(tin1,ie1,ie2)
c      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
cc     angular distribution simulation for n+p and p+p
cc     for energies < 2GeV using Duarte's approximations:
cc     n+p : ds/dom ~ exp(B*t)+a*exp(B*u)+c*exp(alc*u)
cc     P+p:  ds/dom ~ exp(B*t)+a*exp(B1*t)
cc     with Mandelstam's variables t,u
cc     H.Duarte http://www.fjfi.cvut.cz/con_adtt99/papers/Mo-o-c17.pdf
c      real*8 mn
c      data mn/0.939/
c      data  pi /3.1415926536d0/
cc
c      elab=tin1*1000.
c      ps2=tin1*mn/2.
c      if(ie1.eq.ie2)   then      ! p +p or n + n
c        if(elab.le.300.)      then
c          B = 0.
c          a = 0.2
c          B1= 0.
c        elseif(elab.le.670.)  then
c          B = 9.87E-8*(elab-300.)**3
c          a = 0.2
c          B1= 0.
c        elseif(elab.le.1100.)  then
c          B = 4.56E-3*(elab-670.)+4.76
c          a = 97.02E+3*EXP(-2.0E-2*elab)+0.053
c          B1= 9.72E-8*EXP(-5.0E-3*(elab-670.))*(elab-670.)**3  
c        else
c          B = 7.4/(1.+3.0E+5/(elab-300.)**2.23)
c          a = 0.28*EXP(-1.5E-3*elab)
c          B1= 1.94*EXP(-7.0E-4*elab)   
c        endif 
c        if(B.le.0.0001)  then
c          CN1 = 4.*pi 
c        else
c          CN1 = 2.*pi*  (1.-EXP(-2.*ps2*B ))/(ps2*B) 
c        endif
c        if(B1.le.0.0001) then
c          CN2 = 4.*pi*a
c        else
c          CN2=  2.*pi*a*(1.-EXP(-2.*ps2*B1))/(ps2*B1)
c        endif
c        CN = CN1+CN2
c        C1 = CN1/CN
c        C2 = CN2/CN 
c        r1=RNDMD(-1)
c        tm=-2.0*ps2  
c        if(r1.le.C1)   then
c          r2=RNDMD(-1)
c          if(B.le.0.0001)  then
c            t= tm*r2
c          else
c            t=LOG(1.0-r2*(1.0-EXP(tm*B)))/B
c          endif
c        else            
c          r2=RNDMD(-1)
c          if(B1.le.0.0001)  then
c            t= tm*r2
c          else
c            t=LOG(1.0-r2*(1.0-EXP(tm*B1)))/B1
c          endif
c        endif 
c        cts=1.0-t/tm
c        if(RNDMD(-1).le.0.5) then
c          cduarteq = cts   
c        else
c          cduarteq =-cts
c        endif 
cc
c      elseif((ie1.eq.0.and.ie2.eq.1).            ! n + p
c     &    or.(ie1.eq.1.and.ie2.eq.0))   then     ! p + n
c        B = 25.0*EXP(-6.0E-3*elab)+2.0E-3*elab+2.8
c        a = (3.0E-3*elab+0.1)*EXP(1.55-4.9E-3*elab)+1.0E-4*elab
c        c = 6.0E-5*elab*elab*EXP(-6.0E-3*elab)+15.0E-5*elab
c        if(elab.lt.100.)  alc=330.-elab 
c        if(elab.ge.100.)  alc=80.+15.0E+3/elab
c        CN1 = pi/(ps2*B)*(1.-EXP(-4.*ps2*B))  
c        CN2 = CN1*a
c        CN3 = pi/(ps2*alc)*c*(1.-EXP(-4.*ps2*alc))
c        CN = CN1+CN2+CN3
c        C1=CN1/CN
c        C2=CN2/CN
c        C3=CN3/CN 
c        tm=-4.0*ps2
c        um=tm  
c        r1 = RNDMD(-1)
c        if(r1.le.C1)         then
c          r2=RNDMD(-1)
c          if(B.le.0.0001)     then
c            t= tm*r2
c          else
c            t=LOG(1.0-r2*(1.0-EXP(tm*B)))/B
c          endif
c          cts = 1.0-2.0*t/tm 
c        elseif(r1.le.(C1+C2)) then   
c          r2=RNDMD(-1)
c          if(B.le.0.0001)     then
c            u= um*r2
c          else
c            u=LOG(1.0-r2*(1.0-EXP(um*B)))/B
c          endif
c          cts =-1.0+2.0*u/um
c        else  
c          r2=RNDMD(-1)
c          if(alc.le.0.0001)     then
c            u= um*r2
c          else
c            u=LOG(1.0-r2*(1.0-EXP(um*alc)))/alc
c          endif
c          cts =-1.0+2.0*u/um
c        endif
c        if(ie1.eq.1.and.ie2.eq.0) then
cc  parametrization is for n + p case        
c          cduarteq=-cts
c        else
c          cduarteq= cts
c        endif  
c
c      else
c        write(*,*) ' cduarteq: ie1 and ie2 are wrong, ie1,ie2=',ie1,ie2
c        stop
c      endif
c      if(ABS(cduarteq).gt.1.0)         then
c        cduarteq=sign(1.d0,cduarteq)
c      elseif(ABS(cduarteq).lt.1.0d-10) then
c        cduarteq=0.0
c      endif
c      return
c      end   
c   
c*********************************************************************
c
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      double precision function costa1(j,t,r1)
      implicit real*8 (a-h, o-z), integer (i-n)
c     written by KKG  13.05.08
cc
C  CALCULATION OF COSINE BY MEANS OF THE TABLE COEFFICIENTS
      COMMON /COEFA/ ankj(4,4,29)
      dimension tk(4), r1n(5)
      data zro, one, two /0.d0, 1.d0, 2.d0/

c ======================================================================

      s1 = zro
      r1n(1) = one
      r1n(2) = r1
      r1n(3) = r1*r1
      r1n(4) = r1*r1n(3)
      r1n(5) = r1n(3)*r1n(3)
      tk(1) = one
      tk(2) = t
      tk(3) = t*t
      tk(4) = t*t*t
      s2 = zro
      do 10 n = 1,4
      do 10 k = 1,4
        term = ankj(n,k,j)*tk(k)
        s1 = s1 + term*r1n(n)
        s2 = s2 + term
   10 continue
      cta = two*sqrt(r1)*(s1 + (one - s2)*r1n(5)) - one
      temp1 = abs(cta)
      if (temp1.le.one) then
        costa1 = cta
      else
        costa1 = sign(one, cta)
      endif
      return
      end
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      double precision function costan (j0, t0)
c     written by KKG  13.05.08

c ======================================================================
c
c     Cosine calculation for elastic and charge-exchange reactions.
c     cmjn "n" for NEW; extrapolates old approximations from finite
c     angles different from 0 and 180 deg.  The old approximations
c     near 0 and 180 exhibit unphysical "pole-like" behavior.
c
c ======================================================================

      implicit real*8 (a-h, o-z), integer (i-n)

c ======================================================================

      common /COEFA/   ankj(4,4,29)
      dimension  tk(4), r1n(5)


      data zro, one, two /0.d0, 1.d0, 2.d0/

c ======================================================================

      t = t0
      j = j0 
      rl1 = 0.2d0
      rl2 = 0.4d0
      rr1 = 0.9d0
      rr2 = 0.95d0
      r1=RNDMD(-1)
      if (r1.lt.rl1 .or. r1.gt.rr2) then
        if (r1.lt.rl1) then
          x1 = zro  
          y1 = -one
          x2 = rl1
          y2 = costa1 (j, t, x2)           
          x3 = rl2
          y3 = costa1 (j, t, x3) 
        else          
          x1 = one
          y1 = one
          x2 = rr2
          y2 = costa1 (j, t, x2)           
          x3 = rr1
          y3 = costa1 (j, t, x3) 
        endif
        d  = (x2 - x3)*x1*x1 + (x3 - x1)*x2*x2 + (x1 - x2)*x3*x3
        da = (x2 - x3)*y1    + (x3 - x1)*y2    + (x1 - x2)*y3
        db = (y2 - y3)*x1*x1 + (y3 - y1)*x2*x2 + (y1 - y2)*x3*x3
        dc = (x2*y3 - x3*y2)*x1*x1 + (x3*y1 - x1*y3)*x2*x2 +
     &       (x1*y2 - x2*y1)*x3*x3
        a = da/d
        b = db/d
        c = dc/d
        cta = a*r1*r1 + b*r1 +c
        go to 20
      else 
c  Old version of angular distributions;
c  Used for r1 >= rl1 and r1 <= rr2
        r1n(1) = one
        r1n(2) = r1
        r1n(3) = r1*r1
        r1n(4) = r1*r1*r1
        r1n(5) = r1n(3)*r1n(3)
        tk(1) = one
        tk(2) = t
        tk(3) = t*t
        tk(4) = t*t*t
        s1 = zro
        s2 = zro
        do n = 1,4
          do k = 1,4
          term = ankj(n,k,j)*tk(k)
          s1 = s1 + term*r1n(n)
          s2 = s2 + term
          end do
        end do
        cta = two*sqrt(r1)*(s1 + (one - s2)*r1n(5)) - one
      endif
   20 temp1 = abs(cta)
      if (temp1.le.one) then
        costan = cta
      else
        costan = sign(one, cta)
      endif
      return
c ======================================================================
      end
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE DHIST1(X,A,B,H,RX,N,W)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 X,A,B,H,RX,W
C     BLOCK OF BILDING OF HISTOGRAMMS.
      DIMENSION RX(N)
      RX(N) = RX(N)+X*W
      IF (X-A) 1,2,2
    1 RX(N-4) = RX(N-4)+W
                      RETURN
    2 IF (X-B) 4,3,3
    3 RX(N-2) = RX(N-2)+W
                      RETURN
    4 L=(X-A)/H
              NL=N-5
                     IF(L.GT.NL) GO TO 5
      RX(L+1)=RX(L+1)+W
      RX(N-1)=RX(N-1)+X*W
                      RX(N-3)=RX(N-3)+W
                                    RETURN
    5 CONTINUE
c-Sob      write(16,6)  X,A,B,H,W
c-Sob    6 FORMAT(25X,'MISTAKE IN DIMENSION OF DHIST '/25X,5(F10.3,5X))
      RETURN
      END
      SUBROUTINE DHIST2(X,A,B,H,RX,N,W,J)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 X,A,B,H,RX,W
C       BLOCK OF BUILDING 2-DIMENSIONAL HISTOGRAMS
      DIMENSION  RX(N,2)
      RX(N,J) = RX(N,J)+X
      IF(X-A) 1,2,2
    1 RX(N-4,J) = RX(N-4,J)+W
                          RETURN
    2 IF(X-B) 4,3,3
    3 RX(N-2,J) = RX(N-2,J)+W
                          RETURN
    4 L = (X-A)/H
                  RX(L+1,J) = RX(L+1,J)+W
      RX(N-1,J) = RX(N-1,J)+X
                          RX(N-3,J) = RX(N-3,J)+W
                                                RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE DHIST3(X,A,B,H,RX,N,M,W,J)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 X,A,B,H,RX,W
C       BLOCK OF BUILDING M-DIMENSIONAL HISTOGRAMS
      DIMENSION  RX(N,M)
      RX(N,J) = RX(N,J)+X*W
      IF(X-A) 1,2,2
    1 RX(N-4,J) = RX(N-4,J)+W
                          RETURN
    2 IF(X-B) 4,3,3
    3 RX(N-2,J) = RX(N-2,J)+W
                          RETURN
    4 L=(X-A)/H
              NL=N-5
                     IF(L.GT.NL) GO TO 5
      RX(L+1,J)=RX(L+1,J)+W
      RX(N-1,J)=RX(N-1,J)+X*W
                          RX(N-3,J) = RX(N-3,J)+W
                                                RETURN
    5 CONTINUE
c-Sob      write(16,6)  X,A,B,H,J
c-Sob    6 FORMAT(25X,'MISTAKE IN DIMENSION OF DHIST3'/25X,5(F10.3,5X),I3)
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE DHIST4(X,A,B,H,RX,N,M,L,W,I,J)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 X,A,B,H,RX,W
      DIMENSION RX(N,M,L)
      RX(I,J,L)=RX(I,J,L)+X*W
      IF(X-A) 1,2,2
    1 RX(I,J,L-4)=RX(I,J,L-4)+W
                          RETURN
    2 IF(X-B) 4,3,3
    3 RX(I,J,L-2)=RX(I,J,L-2)+W
                          RETURN
    4 L1=(X-A)/H
               NL=L-5
                      IF(L1.GT.NL) GO TO 5
      RX(I,J,L1+1) = RX(I,J,L1+1)+W
      RX(I,J,L-1)=RX(I,J,L-1)+X*W
                          RX(I,J,L-3)=RX(I,J,L-3)+W
      RETURN
    5 CONTINUE
c-Sob      write(16,6)  X,A,B,H,W
c-Sob    6 FORMAT(25X,'MISTAKE IN DIMENSION OF DHIST4'/25X,5(F10.3,5X))
      RETURN
             END
C     * * * * * * * * * * * * * * * * *
      BLOCK DATA SIGAR
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /TABELE/SIGMA,ARGUS/TYPECS/ICST,NSICST
      DIMENSION SIGMA(30,37) , ARGUS(30,6) , ICST(37) , NSICST(37)
     *,A1(30,3),A2(30,3),A3(30,3),A4(30,3),
     * A5(30,3),A6(30,3),A7(30,3),A8(30,3),A9(30,2),
     * B1(30,2),B2(30,2),B3(30,2),B4(30,2),EX(30,3),
     *AR1(30,3),AR2(30,2),AR3(30,1)
      EQUIVALENCE (SIGMA(1, 1),A1(1,1)),(SIGMA(1, 4),A2(1,1)),
     *            (SIGMA(1, 7),A3(1,1)),(SIGMA(1,10),A4(1,1)),
     *            (SIGMA(1,13),A5(1,1)),(SIGMA(1,16),A6(1,1)),
     *            (SIGMA(1,19),A7(1,1)),(SIGMA(1,22),A8(1,1)),
     *            (SIGMA(1,25),A9(1,1)),(SIGMA(1,27),B1(1,1)),
     *            (SIGMA(1,29),B2(1,1)),(SIGMA(1,31),B3(1,1)),
     *            (SIGMA(1,33),B4(1,1)),(SIGMA(1,35),EX(1,1)),
     *            (ARGUS(1, 1),AR1(1,1)),(ARGUS(1, 4),AR2(1,1)),
     *            (ARGUS(1, 6),AR3(1,1))
      DATA  A1/
     *  17613.0 ,  330.000 ,  154.000 ,   96.000 ,   70.000 ,
     *  51.000 ,   38.400 ,   30.000 ,   23.600 ,   22.400 ,
     *  22.500 ,   22.700 ,   23.450 ,   24.300 ,   29.000 ,
     *  41.500 ,   47.800 ,   47.600 ,   47.500 ,   46.700 ,
     *  46.000 ,   45.000 ,   42.500 ,   41.200 ,   40.800 ,
     *  40.000 ,   39.800 ,   39.000 ,   39.000 ,   38.500 ,
     *  17613.0 ,  330.000 ,  154.000 ,   96.000 ,   70.000 ,
     *  51.000 ,   38.400 ,   30.000 ,   23.600 ,   22.400 ,
     *  22.500 ,   22.700 ,   23.000 ,   23.500 ,   24.500 ,
     *  25.600 ,   26.000 ,   25.000 ,   23.000 ,   21.500 ,
     *  19.500 ,   17.500 ,   15.000 ,   13.700 ,   12.000 ,
     *  11.000 ,    9.800 ,    8.800 ,    8.500 ,    6.500 ,
     * 20357.0 ,  950.000 ,  480.000 ,  300.000 ,  200.000 ,
     * 156.000 ,  109.000 ,   73.700 ,   49.500 ,   41.000 ,
     *  36.500 ,  34.01067,   32.6228,   31.5561,   30.032 ,
     *  30.761 ,   33.500 ,   35.000 ,   36.500 ,   38.500 ,
     *  39.500 ,   40.500 ,   41.300 ,   42.000 ,   42.500 ,
     *  42.000 ,   41.300 ,   39.750 ,   39.500 ,   39.000 /
      DATA  A2/
     * 20357.0 ,  950.000 ,  480.000 ,  300.000 ,  200.000 ,
     * 156.000 ,  109.000 ,   73.700 ,   49.500 ,   41.000 ,
     *  36.500 ,   34.000,    32.5000,   31.200 ,   28.500 ,
     *  24.000 ,   20.000 ,   19.000 ,   17.300 ,   16.000 ,
     *  15.200 ,   14.000 ,   12.000 ,   11.000 ,   10.900 ,
     *  10.600 ,   10.000 ,    8.000 ,    6.500 ,    6.200 ,
     *   6.000 ,    6.000 ,    6.500 ,    7.000 ,    8.500 ,
     *  10.500 ,   16.000 ,   25.300 ,   57.500 ,   66.000 ,
     *  63.000 ,   50.000 ,   40.500 ,   26.000 ,   29.000 ,
     *  32.100 ,   45.600 ,   39.000 ,   45.000 ,   57.500 ,
     *  57.200 ,   44.000 ,   36.000 ,   35.700 ,   35.000 ,
     *  34.200 ,   32.000 ,   30.200 ,   26.200 ,   25.500 ,
     *   2.000 ,    2.000 ,    2.250 ,    2.300 ,    2.400 ,
     *   2.900 ,    4.00 ,     7.700 ,   20.000 ,   22.000 ,
     *  23.900 ,   17.500 ,   13.970 ,   11.000 ,   11.200 ,
     *  14.000 ,   20.000 ,   15.500 ,   18.800 ,   26.000 ,
     *  25.800 ,   18.300 ,   12.700 ,   11.000 ,   10.300 ,
     * 9.700 ,    6.700 ,    5.900 ,    4.700 ,    4.100 /
      DATA  A3/
     *   4.000 ,    4.000 ,    4.250 ,    4.700 ,    6.100 ,
     *   7.600 ,   12.000 ,   17.600 ,   37.500 ,   44.000 ,
     *  39.100 ,   32.500 ,   26.400 ,   12.500 ,   10.900 ,
     *   9.100 ,    7.800 ,    4.700 ,    6.500 ,    9.000 ,
     *   7.000 ,    3.500 ,    2.500 ,    2.100 ,    1.900 ,
     * 1.400 ,    0.260 ,    0.200 ,    0.040 ,    0.020 ,
     *   1.900 ,    2.300 ,    3.500 ,    5.500 ,    9.000 ,
     *  14.000 ,   28.000 ,   55.000 ,  160.000 ,  195.000 ,
     * 180.000 ,  145.000 ,  113.000 ,   42.000 ,   25.100 ,
     *   20.300 ,   15.000 ,   16.500 ,   22.500 ,   23.200 ,
     *  26.000 ,   27.200 ,   37.800 ,   41.000 ,   39.000 ,
     *  31.500 ,   28.900 ,   27.500 ,   24.800 ,   23.500 ,
     *   1.900 ,    2.300 ,    3.500 ,    5.500 ,    9.000 ,
     *  14.000 ,   28.000 ,   55.000 ,  160.000 ,  195.000 ,
     * 180.000 ,  145.000 ,  113.000 ,   41.680 ,   24.000 ,
     *  18.600 ,   11.250 ,    8.300 ,   10.200 ,   10.800 ,
     *  11.800 ,   13.000 ,   17.300 ,   19.000 ,   16.000 ,
     *  12.800 ,   7.500 ,    6.400 ,    4.800 ,    4.000 /
      DATA  A4/
     *  20.000 ,   17.000 ,   17.000 ,   19.000 ,   22.000 ,
     *  24.500 ,   32.500 ,   40.000 ,   46.500 ,   42.000 ,
     *  32.000 ,   22.000 ,   15.500 ,    4.400 ,    2.000 ,
     *   1.380 ,    0.760 ,    0.550 ,    0.380 ,    0.330 ,
     *   0.270 ,    0.220 ,    0.140 ,    0.    ,    0.    ,
     *   0.    ,    0.    ,    0.    ,    0.    ,    0.    ,
     *   0.    ,    0.    ,    0.    ,    0.050 ,    0.100 ,
     *    .300 ,     .600 ,    1.200 ,    2.400 ,    3.200 ,
     *   3.400 ,    3.600 ,    3.700 ,    3.800 ,    3.900 ,
     *   3.900 ,    3.900 ,    4.000 ,    4.000 ,    4.000 ,
     *   4.000 ,    4.000 ,    3.900 ,    3.800 ,    3.500 ,
     *   3.200 ,    2.900 ,    2.700 ,    2.400 ,    2.100 ,
     *   0.    ,    0.    ,    0.    ,    0.400 ,    0.700 ,
     *   2.000 ,    3.900 ,    5.200 ,    9.000 ,   12.300 ,
     *  14.400 ,   16.100 ,   16.600 ,   17.000 ,   17.300 ,
     *  17.500 ,   17.600 ,   17.700 ,   17.800 ,   17.900 ,
     *  17.800 ,   17.500 ,   16.800 ,   16.000 ,   13.700 ,
     *  12.600 ,  11.600 ,   10.800 ,    9.200 ,    8.000 /
      DATA  A5/
     *   0.    ,    0.    ,  0.009632,   0.09939,   0.2666 ,
     *   0.5429,    0.9929,    1.704 ,    2.707 ,    3.851 ,
     *   4.887 ,    5.689 ,    6.261 ,    6.653 ,    6.915 ,
     *   7.086 ,    7.197 ,    7.324 ,    7.412 ,    7.494 ,
     *   7.576 ,    7.755 ,    7.900 ,    7.700 ,    7.000 ,
     *   6.600 ,    6.200 ,    5.800 ,    4.600 ,    3.500 ,
     *   0.    ,    0.    ,   0.00052,   0.01166,   0.04471,
     *   0.1182,    0.2695,    0.5573,    0.9979,    1.455 ,
     *    1.782 ,    1.982 ,    2.107 ,    2.191 ,    2.255 ,
     *    2.310 ,    2.366 ,    2.504 ,    2.698 ,    2.937 ,
     *    3.199 ,    3.759 ,    4.000 ,    3.900 ,    3.500 ,
     *    3.200 ,    2.700 ,    2.100 ,    1.000 ,    0.600 ,
     *   0.    ,    0.    ,    0.150 ,    0.300 ,    0.500 ,
     *   0.850 ,    1.300 ,    1.900 ,    3.000 ,    4.500 ,
     *   6.900 ,    9.300 ,   10.000 ,   10.200 ,   10.000 ,
     *   9.200 ,    8.900 ,   11.000 ,   12.300 ,   11.700 ,
     *   9.900 ,    7.600 ,    5.300 ,    5.200 ,    4.200 ,
     *   3.800 ,    3.400 ,    3.000 ,    2.300 ,    2.100 /
      DATA  A6/
     *   0.    ,    0.    ,    0.050 ,    0.120 ,    0.180 ,
     *   0.250 ,    0.400 ,    0.550 ,    0.750 ,    1.000 ,
     *   1.300 ,    1.600 ,    1.900 ,    2.000 ,    2.300 ,
     *   2.500 ,    2.600 ,    2.900 ,    3.500 ,    4.000 ,
     *   3.800 ,    3.100 ,    2.500 ,    2.250 ,    2.400 ,
     *   2.500 ,    2.300 ,    2.000 ,    1.500 ,    1.300 ,
     *   0.    ,    0.    ,    0.300 ,    0.800 ,    1.300 ,
     *   1.700 ,    1.900 ,    2.000 ,    2.100 ,    2.200 ,
     *   2.400 ,    2.700 ,    3.100 ,    3.200 ,    3.000 ,
     *   2.500 ,    2.100 ,    1.600 ,    1.400 ,    1.300 ,
     *   1.100 ,    0.700 ,    0.500 ,    0.450 ,    0.300 ,
     *   0.270 ,    0.250 ,    0.240 ,    0.210 ,    0.200 ,
     *   0.    ,    0.010 ,    0.100 ,    0.300 ,    0.600 ,
     *   1.100 ,    2.100 ,    3.500 ,    4.800 ,    4.700 ,
     *   4.200 ,    4.650 ,    5.400 ,    6.500 ,    6.900 ,
     *   6.500 ,    5.400 ,    4.100 ,    4.700 ,    5.000 ,
     *   4.800 ,    4.400 ,    4.000 ,    3.700 ,    3.100 ,
     *   2.800 ,    2.600 ,    2.500 ,    2.200 ,    1.900 /
      DATA  A7/
     *   0.    ,    0.120 ,    0.600 ,    1.400 ,    2.800 ,
     *   4.100 ,    5.000 ,    5.800 ,    7.000 ,    6.200 ,
     *   6.400 ,    8.000 ,   10.000 ,   11.300 ,   11.400 ,
     *   9.800 ,    8.000 ,    7.400 ,    8.000 ,    8.000 ,
     *   7.700 ,    6.800 ,    5.900 ,    5.200 ,    4.200 ,
     *   2.900 ,    2.700 ,    2.500 ,    2.200 ,    1.900 ,
     *   310.00,   270.00,   240.00,   213.00,   172.00,
     *   156.00,   132.00,   116.00,   108.00,   102.00,
     *    98.00,    96.00,    96.00,    98.00,   102.00,
     *   109.00,   116.00,   122.00,   128.00,   132.00,
     *   136.00,   138.00,   140.00,   140.00,   140.00,
     *   141.00,   140.00,   138.00,   134.00,   132.00,
     *   204.00,   165.00,   138.00,   117.00,    87.00,
     *    75.00,    55.00,    40.00,    30.00,    24.00,
     *    22.00,    20.00,    20.00,    20.00,    21.00,
     *    23.00,    26.00,    28.00,    31.00,    33.00,
     *    34.00,    36.00,    36.00,    36.00,    36.00,
     *    36.00,     36.00,     36.00,     34.00,     32.00/
      DATA  A8/
     *   266.00,   228.00,   200.00,   178.00,   148.00,
     *   136.00,   115.00,   102.00,    94.00,    86.00,
     *    79.00,    76.00,    76.00,    78.00,    81.00,
     *    85.00,    91.00,   100.00,   107.00,   112.00,
     *   115.00,   118.00,   118.00,   118.00,   118.00,
     *   118.00,   116.00,   112.00,   110.00,   106.00,
     *   133.00,   105.00,    84.00,    68.00,    46.00,
     *    38.00,    28.00,    22.00,   19.00,    15.00,
     *    11.00,   10.00,     9.50,    10.00,   11.00,
     *    12.00,    15.00,    17.00,    19.00,    21.00,
     *    23.00,    24.00,    24.00,    24.00,    23.00,
     *    22.00,    22.00,    21.00,   20.00,    20.00,
     *     0.00,     0.00,     0.00,     0.20,     0.80,
     *     1.60,     3.20,     6.20,    33.00,    40.30,
     *    41.60,    42.00,    42.00,    42.00,    41.90,
     *    41.50,    40.70,    40.00,    39.30,    38.40,
     *    37.60,    36.60,    35.40,    34.20,    32.80,
     *    31.20,     29.40,     27.70,     26.20,     16.00/
      DATA  A9/
     *000000.00,     0.00,     0.00,     0.30,     1.20,
     *     2.20,     3.40,     4.90,     9.00,    11.20,
     *    13.40,    15.60,    18.10,    20.30,    22.30,
     *    25.40,    27.10,    27.90,    28.30,    28.60,
     *    28.80,    29.00,    29.10,    29.20,    29.30,
     *    29.30,    29.30,    29.30,    29.30,    29.30,
     *     0.00,     1.00,    30.00,    41.00,    43.20,
     *    43.60,    43.40,    43.00,    41.00,    39.00,
     *    32.00,    23.50,    16.00,    12.00,     9.40,
     *     5.60,     3.80,     2.90,     2.20,    1.80,
     *     1.60,     1.30,     1.00,     0.80,     0.70,
     *     0.50,      0.40,      0.40,      0.30,      0.00/
      DATA  B1/
     * 10.0, 11.0, 11.5, 11.9, 12.0, 12.1, 12.2, 12.2, 12.8, 14.5,
     * 16.3, 18.9, 18.0, 17.6, 17.2, 17.1, 17.0, 17.0, 17.1, 17.2,
     * 17.5, 17.8, 18.0, 19.0, 19.4, 20.0, 20.7, 21.4, 22.0, 24.0,
     * 10.0, 11.0, 11.5, 11.9, 12.0, 12.1, 12.2, 12.1, 12.0, 11.9,
     * 11.6, 10.5,  9.0,  6.5,  5.0,  4.1,  3.5,  3.2,  3.2,  3.1,
     *  2.9,  2.5,  2.4,  2.4,  2.4,  2.4,  2.4,  2.4,  2.4,  2.4/
      DATA  B2/
     *  6.0,  7.0,  7.5, 10.0, 12.0, 13.0, 13.0, 13.5, 15.0, 16.5,
     * 18.5, 20.5, 19.0, 18.5, 17.8, 17.3, 17.2, 17.2, 17.2, 17.5,
     * 17.7, 18.0, 18.2, 18.5, 19.0, 19.7, 20.0, 20.3, 20.6, 21.0,
     *  3.0,  4.0,  4.5,  6.0,  8.0,  5.5,  5.5,  5.7,  5.7,  5.7,
     *  5.7,  5.5,  5.2,  5.0,  4.8,  4.5,  4.2,  3.8,  3.5,  3.1,
     *  2.3,  2.5,  2.4,  2.4,  2.4,  2.4,  2.4,  2.4,  2.4,  2.4/
      DATA  B3/
     * 15.0, 15.0, 15.0, 15.0, 15.0, 15.0, 20.0, 26.0, 29.0, 31.0,
     * 31.0, 29.0, 28.0, 23.0, 22.0, 21.0, 20.7, 20.5, 20.3, 20.2,
     * 20.0, 20.0, 20.0, 20.0, 20.1, 20.5, 21.0, 22.0, 23.0, 24.0,
     *  5.0,  5.0,  5.0,  5.0,  5.0,  5.0,  7.0,  9.0, 15.0, 18.0,
     * 17.0, 10.0,  7.5,  5.0,  3.4,  3.1,  3.0,  2.9,  2.8,  2.7,
     * 2.65, 2.60, 2.50, 2.51, 2.52, 2.55, 2.56, 2.57, 2.58, 2.58/
      DATA  B4/
     *450.0,323.0,125.0, 82.0, 65.0, 38.0,30.75, 29.0, 32.5, 34.0,
     * 50.0, 32.0, 33.0, 30.0, 28.0, 27.0, 25.0, 24.0, 23.0, 22.0,
     * 21.0, 20.5, 20.0, 20.0, 20.0, 20.5, 21.0, 21.5, 22.0, 22.5,
     *150.0, 98.0, 60.0, 42.0, 33.0, 23.0, 18.0, 16.0, 17.0, 20.0,
     * 22.0, 15.0,  9.0,  8.0,  5.5,  4.5,  4.0,  3.6,  3.2,  2.9,
     *  2.8,  2.6,  2.5, 2.51, 2.52, 2.55, 2.56, 2.57, 2.58, 2.58/
      DATA  EX/
     *  3.0,  3.0,  3.0,  4.0,  4.0,  7.5,  7.5,  7.8,  7.8,  7.8,
     *  7.8,  4.5,  3.0,  1.8,  0.8, 0.45, 0.30, 0.15,0.075,0.035,
     *  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     * 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0,  9.0,  8.0,  7.0,
     *  6.5,  4.0, 2.75, 1.50, 0.65, 0.30,  0.0,  0.0,  0.0,  0.0,
     *  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     *300.0,225.0, 65.0, 40.0, 32.0, 15.0,12.75, 12.0, 15.5, 14.0,
     * 9.75,  5.0,  4.0,  2.0, 1.25,  0.5,  0.0,  0.0,  0.0,  0.0,
     *  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0/
      DATA  ICST/
     *  210 ,   211 ,   220 ,   221 ,   120 ,   121 ,   122 ,
     *  110 ,   111 ,   123 ,   214 ,   215 ,   224 ,   225 ,
     *  114 ,   115 ,   126 ,   124 ,   125 ,
     *00510 , 00511 , 00410 , 00411 , 00514 , 00515 , 00516,
     * 1110 ,  1111 ,  1120 ,  1121 , -1110 , -1111 , -1120,
     *-1121 ,  1122 , -1112 , -1122 /
      DATA  NSICST/
     *  112 ,   113 ,   116 ,   117 ,   127 ,   130 ,   131 ,
     *  132 ,   133 ,   134 ,   135 ,   136 ,   137 ,   212 ,
     *  213 ,   216 ,   217 ,   222 ,   223 ,   226 ,   227 ,
     *  230 ,   231 ,   232 ,   233 ,
     *  240 ,   241 ,   242 ,   243 ,
     *  250 ,   251 ,   252 ,   253 ,
     *  260 ,   261 ,   262 ,   263 /
      DATA  AR1/
     *      0.0 ,    0.010 ,    0.020 ,    0.030 ,    0.040 ,
     *    .050 ,     .070 ,     .100 ,     .150 ,     .200 ,
     *    .250 ,     .300 ,     .350 ,     .400 ,     .500 ,
     *    .650 ,     .850 ,     .950 ,    1.100 ,    1.300 ,
     *   1.500 ,    2.000 ,    3.000 ,    4.000 ,    5.000 ,
     *   7.000 ,   10.000 ,   16.000 ,   22.000 ,   30.000 ,
     *    .000 ,     .010 ,     .020 ,     .030 ,     .040 ,
     *    .050 ,     .075 ,     .100 ,     .150 ,     .175 ,
     *    .200 ,     .225 ,     .250 ,     .350 ,     .450 ,
     *    .500 ,     .600 ,     .700 ,     .800 ,     .850 ,
     *    .900 ,    1.000 ,    1.200 ,    1.300 ,    1.400 ,
     *   1.600 ,    3.000 ,    4.000 ,   10.000 ,   20.000 ,
     *    .200 ,     .250 ,     .300 ,     .350 ,     .400 ,
     *    .450 ,     .500 ,     .550 ,     .600 ,     .650 ,
     *   0.700 ,    0.750 ,    0.800 ,    0.850 ,    0.900 ,
     *    .950 ,    1.000 ,    1.100 ,    1.200 ,    1.300 ,
     *   1.400 ,    1.600 ,    1.800 ,    2.000 ,    2.400 ,
     *   2.600 ,    2.800 ,    3.000 ,    3.500 ,    4.000 /
      DATA  AR2/
     *   0.05  ,   0.06  ,   0.07  ,   0.08  ,   0.10  ,
     *   0.11  ,   0.13  ,   0.15  ,   0.17  ,   0.20  ,
     *   0.25  ,   0.30  ,   0.35  ,   0.40  ,   0.45  ,
     *   0.50  ,   0.55  ,   0.60  ,   0.65  ,   0.70  ,
     *   0.75  ,   0.80  ,   0.85  ,   0.90  ,    1.00  ,
     *   2.00  ,   3.00  ,   5.00  ,   10.0  ,   30.0  ,
     *   0.0185,   0.0200,   0.0224,   0.0250,   0.0282,
     *   0.0316,   0.0355,   0.0400,   0.0500,   0.0560,
     *   0.0630,   0.0710,   0.0800,   0.0900,   0.1000,
     *   0.1260,   0.1570,   0.2000,   0.2500,   0.3160,
     *   0.4000,   0.5000,   0.6300,   0.8000,   1.0000,
     *   1.2600,   1.5700,   2.0000,   2.5000,   30.000/
      DATA  AR3/
     *   0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
     *   1.0, 1.25, 1.5, 2.0, 3.0, 4.0, 5.0, 7.0, 10.0,15.0,
     *  20.0,30.0, 50.0, 100.0, 150.0, 200.0, 300.0, 400.0,
     * 500.0, 1000.0/
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      BLOCK DATA COEFAB
      IMPLICIT REAL*8 (A-H,O-Z)
C     SUBROUTINE WHICH PUT THE COEFFICIENTS ANKJ,BNKJ AND CKJ
C     IN THE MAIN PROGRAM
c  kkg  10/28/03
      COMMON /COEFA/ ANKJ  /COEFBC/ BNKJ,CKJ
      DIMENSION ANKJ(4,4,29),BNKJ(4,4,8),CKJ(3,8)
     *,A1(4,4,4),A2(4,4,4),A3(4,4,4),A4(4,4,4),
     * A5(4,4,4),A6(4,4,1),A7(4,4,7),A8(4,4,1),
     * B1(4,4,4),B2(4,4,4)
      EQUIVALENCE (ANKJ(1,1, 1),A1(1,1,1)),
     *            (ANKJ(1,1, 5),A2(1,1,1)),
     *            (ANKJ(1,1, 9),A3(1,1,1)),
     *            (ANKJ(1,1,13),A4(1,1,1)),
     *            (ANKJ(1,1,17),A5(1,1,1)),
     *            (ANKJ(1,1,21),A6(1,1,1)),
     *            (ANKJ(1,1,22),A7(1,1,1)),      
     *            (ANKJ(1,1,29),A8(1,1,1)),      
     *            (BNKJ(1,1, 1),B1(1,1,1)),
     *            (BNKJ(1,1, 5),B2(1,1,1))
      DATA  A1/
c  ankj Angular distribution coefficients:
c  j = 1;  N + N elastic scattering; Tlab <= 2.8 GeV:
c          (n + n & p + p isotropic below 0.46 GeV.)
     * 2.7404E 00 , -9.6998E 00 ,  1.0400E 01 ,  2.3882E 00 ,
     *-7.5137E 00 ,  4.4096E 01 , -7.4379E 01 ,  4.6038E 01 ,
     * 7.5479E 00 , -3.9274E 01 ,  6.4835E 01 , -4.1609E 01 ,
     *-1.8369E 00 ,  8.6911E 00 , -1.3060E 01 ,  7.1880E 00 ,
c  j = 2;  N + N elastic scattering; 2.8 < Tlab <= 10. GeV:
     *-3.0853E 01 ,  1.0624E 02 , -1.2939E 02 ,  5.4339E 01 ,
     * 1.9465E 01 , -6.8102E 01 ,  9.6358E 01 , -5.6827E 01 ,
     *-3.4831E 00 ,  1.2341E 01 , -1.8592E 01 ,  1.2024E 01 ,
     * 1.8941E-01 , -6.7880E-01 ,  1.0665E 00 , -7.2910E-01 ,
c  j = 3;  n + p elastic scattering; Tlab <= 0.97 GeV:
     * 1.0258E-01 , -1.0542E 00 ,  1.1389E 01 , -1.6638E 01 ,
     *-4.9607E-01 ,  1.1800E 01 , -9.0857E 01 ,  1.6476E 02 ,
     * 1.5437E 00 , -3.3769E 01 ,  2.5192E 02 , -4.5071E 02 ,
     *-1.2021E 00 ,  2.5336E 01 , -1.8658E 02 ,  3.3254E 02 ,
c  j = 4; pi+ p or pi- n elastic scattering; Tlab <= 0.080 GeV:
     * 1.5789E-01 ,  2.9671E 00 , -5.5251E 00 ,  6.8925E 00 ,
     *-7.0218E 00 , -2.0534E 02 ,  5.6951E 02 , -8.9858E 02 ,
     * 1.3496E 02 ,  4.8722E 03 , -1.4674E 04 ,  2.3924E 04 ,
     *-8.2116E 02 , -3.2586E 04 ,  1.0098E 05 , -1.6553E 05 /
      DATA  A2/
c  j = 5; pi+ p or pi- n elastic scattering; 0.08 < Tlab <= 0.3 GeV:
     * 3.1531E-01 , -7.4981E 00 ,  4.3295E 01 , -7.6360E 01 ,
     *-6.5373E 00 ,  1.9307E 02 , -1.0181E 03 ,  1.7426E 03 ,
     * 4.6864E 01 , -1.3030E 03 ,  6.7291E 03 , -1.1075E 04 ,
     *-9.5192E 01 ,  2.6373E 03 , -1.2857E 04 ,  2.0294E 04 ,
c  j = 6; pi+ p or pi- n elastic scattering; 0.30 < Tlab <= 1.0 GeV:
     *-1.7953E 01 ,  1.0972E 02 , -2.3954E 02 ,  2.2826E 02 ,
     * 9.1968E 01 , -5.1963E 02 ,  1.1266E 03 , -1.0740E 03 ,
     *-1.3270E 02 ,  7.4112E 02 , -1.6000E 03 ,  1.5249E 03 ,
     * 5.8598E 01 , -3.1874E 02 ,  6.7751E 02 , -6.4011E 02 ,
c  j = 7; pi+ p or pi- n elastic scattering; 1.0 < Tlab <= 2.4 GeV:
     * 4.2169E-01 ,  1.4705E 02 , -6.5335E 02 ,  9.1507E 02 ,
     *-3.5198E 00 , -2.6019E 02 ,  1.2250E 03 , -1.7481E 03 ,
     * 3.6373E 00 ,  1.5592E 02 , -7.5201E 02 ,  1.0796E 03 ,
     *-7.8041E-01 , -3.0563E 01 ,  1.4795E 02 , -2.1250E 02 ,
c  j = 8; pi+ n or pi- p elastic scattering; Tlab <= 0.080 GeV:
     *-3.8288E-01 ,  3.7587E 00 , -6.5144E 00 ,  6.7740E 00 ,
     * 1.0381E 02 , -2.7282E 02 ,  4.7759E 02 , -5.1222E 02 ,
     *-1.7882E 03 ,  4.3052E 03 , -7.9314E 03 ,  9.3471E 03 ,
     * 7.1475E 03 , -3.3395E 03 , -4.1392E 03 , -4.4364E 03 /
      DATA  A3/
c  j = 9; pi- p or pi+ n elastic scattering; 0.08 < Tlab <= 0.3 GeV:
     * 2.4991E-01 ,  3.2028E 01 , -1.1882E 02 ,  1.5099E 02 ,
     *-2.6994E 00 , -4.6045E 02 ,  1.8959E 03 , -2.5190E 03 ,
     * 1.6268E 01 ,  2.1384E 03 , -9.1262E 03 ,  1.2431E 04 ,
     *-2.9654E 01 , -3.1823E 03 ,  1.3944E 04 , -1.9342E 04 ,
c  j = 10; pi- p or pi+ n elastic or CX scattering;
c          0.30 < Tlab <= 1.0 GeV:
     * 3.9025E 00 , -9.1126E 01 ,  3.2373E 02 , -4.0048E 02 ,
     *-2.0619E 01 ,  4.9170E 02 , -1.7155E 03 ,  2.1143E 03 ,
     * 3.3004E 01 , -7.6684E 02 ,  2.7003E 03 , -3.3525E 03 ,
     *-1.6367E 01 ,  3.7394E 02 , -1.3202E 03 ,  1.6423E 03 ,
c  j = 11; pi- p or pi+ n elastic or CX scattering; 
c          1.0 < Tlab <= 2.4 GeV:
     * 1.9402E 01 , -2.2446E 02 ,  7.4733E 02 , -9.3570E 02 ,
     *-4.4180E 01 ,  4.7194E 02 , -1.4856E 03 ,  1.8055E 03 ,
     * 3.1567E 01 , -3.0176E 02 ,  9.0763E 02 , -1.0773E 03 ,
     *-6.8648E 00 ,  6.0476E 01 , -1.7520E 02 ,  2.0381E 02 ,
c  j = 12; pi- + p --> pi0 n of pi+ + n --> pi0 + p Charge exchange
c          scattering; Tlab <= 0.08 GeV:
     * 1.4988E-01 ,  2.8753E 00 , -5.3078E 00 ,  6.2233E 00 ,
     *-5.9558E 00 , -1.6203E 02 ,  4.3079E 02 , -6.2548E 02 ,
     * 1.2875E 02 ,  3.1402E 03 , -7.9189E 03 ,  1.0983E 04 ,
     *-8.5161E 02 , -1.8780E 04 ,  4.4607E 04 , -5.8790E 04 /
      DATA  A4/
c  j = 13; pi- + p --> pi0 n of pi+ + n --> pi0 + p Charge exchange
c          scattering; 0.08 < Tlab <= 0.30 GeV:
     * 5.3689E-01 , -1.3216E 01 ,  8.1011E 01 , -1.4285E 02 ,
     *-1.0550E 01 ,  2.9629E 02 , -1.6957E 03 ,  2.8935E 03 ,
     * 6.9621E 01 , -1.9245E 03 ,  1.0620E 04 , -1.7468E 04 ,
     *-1.3865E 02 ,  3.9281E 03 , -2.0293E 04 ,  3.2058E 04 ,
c  j = 14;  N + N --> N + N + pi; nucleon distributions:
     * 8.5591E-02 ,  5.0390E 00 , -1.3782E 01 ,  1.4661E 01 ,
     * 5.4284E-02 , -9.2324E 00 ,  3.6397E 01 , -4.2962E 01 ,
     *-5.1111E-02 ,  4.6003E 00 , -2.0534E 01 ,  2.7731E 01 ,
     * 7.4514E-03 , -6.2529E-01 ,  2.9159E 00 , -4.1101E 00 ,
c  j = 15;  N + N --> N + N + pi; pion distributions:
     * 7.1622E-02 ,  3.0960E 00 , -1.1125E 01 ,  1.8130E 01 ,
     * 9.2581E-02 , -3.2186E 00 ,  2.0273E 01 , -3.3245E 01 ,
     *-5.1531E-02 ,  8.9886E-01 , -7.5084E 00 ,  1.3188E 01 ,
     * 5.8258E-03 , -1.7288E-03 ,  7.0224E-01 , -1.4854E 00 ,
c  j = 16;  N + N --> N + N + n*pi, n > 1; nucleon distributions:
     * 8.2300E-02 ,  1.5854E-01 ,  3.7716E 00 , -4.0562E 00 ,
     * 1.0802E-02 , -3.3688E-01 ,  1.1727E 00 , -6.7476E-01 ,
     *-2.1798E-03 ,  5.2166E-02 , -2.5816E-01 ,  3.2048E-01 ,
     * 6.5764E-05 , -1.4711E-03 ,  7.8209E-03 , -1.0580E-02 /
      DATA  A5/
c  j = 17;  N + N --> N + N + n*pi, n > 1; pion distributions:
     * 1.1138E-01 ,  6.0396E-01 ,  3.0174E 00 , -4.4190E 00 ,
     *-1.7709E-02 ,  2.3015E-01 , -1.8187E 00 ,  3.4518E 00 ,
     * 2.0977E-03 , -2.5458E-02 ,  2.1626E-01 , -4.0692E-01 ,
     *-5.4799E-05 ,  5.9111E-04 , -5.5552E-03 ,  1.0647E-02 ,
c  j = 18;  pi + N --> pi + N + pi; nucleon distributions:
     * 1.7288E-01 ,  7.1080E 00 , -1.7961E 01 ,  1.6403E 01 ,
     *-1.4504E-01 , -1.3032E 01 ,  4.1781E 01 , -4.0799E 01 ,
     * 4.5390E-02 ,  8.3515E 00 , -3.0260E 01 ,  3.2882E 01 ,
     *-4.7961E-03 , -1.4095E 00 ,  5.3505E 00 , -6.0946E 00 ,
c  j = 19;  pi + N --> pi + N + pi; pion distributions:
     * 3.7596E-02 ,  1.4331E 00 , -3.1350E 00 ,  6.4864E 00 ,
     * 2.3827E-01 ,  1.8253E 00 ,  1.7648E 00 , -1.6735E 01 ,
     *-1.5410E-01 , -1.5201E 00 , -1.5692E 00 ,  1.7185E 01 ,
     * 2.5037E-02 ,  3.0588E-01 ,  3.2520E-01 , -3.5277E 00 ,
c  j = 20;  pi + N --> pi + N + n*pi, n > 1; nucleon distributions:
     * 1.2489E-01 ,  1.3573E 00 ,  8.2338E-01 , -1.4595E 00 ,
     *-5.1577E-02 , -3.5778E-01 , -1.1690E 00 ,  1.8078E 00 ,
     * 7.4864E-03 ,  3.2888E-02 ,  2.3744E-01 , -3.9802E-01 ,
     *-2.9880E-04 , -7.5117E-04 , -1.1402E-02 ,  1.9505E-02 /
      DATA  A6/
c  j = 21;  pi + N --> pi + N + n*pi, n > 1; pion distributions:
     * 1.8470E-01 ,  1.9269E 00 , -3.2979E 00 ,  3.6843E 00 ,
     *-7.3932E-02 ,  2.7213E-01 ,  1.0600E 00 , -2.3354E 00 ,
     * 1.8907E-02 , -5.6473E-02 , -1.6487E-01 ,  3.8426E-01 ,
     *-9.2984E-04 ,  2.5506E-03 ,  7.3052E-03 , -1.7220E-02 /
      DATA  A7/
c  j = 22;  gamma + N --> pi0 + N, E-g <= 0.45 GeV. 
     & 4.0693D-01 , -4.1404D 00 ,  1.4044D 01 , -1.7265D 01 ,
     &-3.6799D 00 ,  5.9610D 01 , -1.6269D 02 ,  1.8873D 02 ,
     & 1.4556D 01 , -1.7550D 02 ,  4.5839D 02 , -5.3390D 02 ,
     &-1.2621D 01 ,  1.4964D 02 , -3.8118D 02 ,  4.5141D 02 ,
c  j = 23;  gamma + N --> pi0 + N, E-g > 0.45 GeV. 
     &-4.7554D-01 ,  2.2641D 00 , -1.2528D 01 ,  2.4647D 01 ,
     & 5.1620D 00 , -9.9236D 00 ,  5.5623D 01 , -1.0462D 02 ,
     &-8.1117D 00 ,  1.9315D 01 , -8.4255D 01 ,  1.3908D 02 ,
     & 3.5187D 00 , -9.1783D 00 ,  3.4950D 01 , -5.1243D 01 ,
c  j = 24; gamma + p --> n + pi+; E-g <= 0.51 GeV:
     & 4.8173D-01 ,  5.7726D 00 , -1.3745D 01 ,  2.7125D 01 ,
     &-4.4804D 00 , -3.8582D 01 ,  1.1159D 02 , -2.4305D 02 ,
     & 1.6306D 01 ,  1.1046D 02 , -3.3045D 02 ,  7.2270D 02 ,
     &-1.5968D 01 , -8.0140D 01 ,  2.4616D 02 , -6.0753D 02 ,
c  j = 25; gamma + p --> n + pi+; 0.51 < E-g <= 1.0 GeV:
     &-5.1646D 00 , -6.0776D 00 ,  7.8989D 01 , -1.0705D 02 ,
     & 2.1871D 01 ,  5.6915D 01 , -4.0159D 02 ,  5.1215D 02 ,
     &-2.7993D 01 , -9.4670D 01 ,  5.6928D 02 , -6.9621D 02 ,
     & 1.1587D 01 ,  4.5998D 01 , -2.4566D 02 ,  2.8452D 02 ,
c  j = 26; gamma + p --> n + pi+; 1.0 < E-g <= 10 GeV:
     &-5.3067D 01 ,  5.7612D 02 , -1.5438D 03 ,  1.6455D 05 ,
     & 1.4750D 02 , -1.6638D 03 ,  4.5923D 03 , -4.9949D 03 ,
     &-1.3436D 02 ,  1.5780D 03 , -4.4463D 03 ,  4.9022D 03 ,
     & 4.0253D 01 , -4.8860D 02 ,  1.4001D 03 , -1.5606D 03 ,
c  j = 27; gamma + N --> delta + pi; pion distribution; T < 1.0
     &-1.0306D 00 ,  3.2849D 01 , -7.5052D 01 ,  6.0255D 01 ,
     & 7.9586D 00 , -1.2572D 02 ,  2.5604D 02 , -1.6547D 02 ,
     &-1.4797D 01 ,  1.6590D 02 , -2.7991D 02 ,  1.1333D 02 ,
     & 8.2309D 00 , -6.7871D 01 ,  8.5762D 01 ,  5.9727D 00 ,
c  j = 28; gamma + N --> delta + pi; pion distribution; T > 1.0
     &-2.3722D 02 ,  9.6890D 02 , -1.6219D 03 ,  1.3637D 03 ,
     & 6.5800D 02 , -2.6941D 03 ,  4.5480D 03 , -3.8460D 03 ,
     &-6.0653D 02 ,  2.4983D 03 , -4.2498D 03 ,  3.6136D 03 ,
     & 1.8604D 02 , -7.6933D 02 ,  1.3166D 03 , -1.1242D 03 /
      DATA  A8/
c  j = 29; coefficients for K absorption ang. dist., Tlab <= 0.455
     & 6.5288D-01 ,  3.8977D-01 ,  8.4078D-01 ,  1.8893D-01 ,
     &-4.3964D 00 ,  3.4309D 01 , -7.3692D 01 ,  8.4308D 01 ,
     & 1.4889D 01 , -1.4380D 02 ,  3.1227D 02 , -3.5014D 02 ,
     &-1.5658D 01 ,  1.7160D 02 , -3.7212D 02 ,  4.1299D 02 /
      DATA  B1/
     * 5.0278E-01 ,  3.1442E 00 , -7.8172E 00 ,  8.1667E 00 ,
     * 9.3482E-01 , -1.0590E 01 ,  2.9227E 01 , -3.4550E 01 ,
     *-9.6685E-02 ,  4.7335E 00 , -1.4298E 01 ,  1.7685E 01 ,
     *-2.5041E-02 , -6.2478E-01 ,  2.0282E 00 , -2.5895E 00 ,
     * 1.1965E 00 , -8.2889E-01 ,  1.0426E 00 , -1.9090E 00 ,
     * 2.8703E-01 , -4.9065E 00 ,  1.6264E 01 , -1.9904E 01 ,
     *-2.4492E-01 ,  2.9191E 00 , -9.5776E 00 ,  1.1938E 01 ,
     * 3.7297E-02 , -4.2200E-01 ,  1.3883E 00 , -1.7476E 00 ,
     * 1.3508E 00 , -4.3139E 00 ,  1.2291E 01 , -1.5288E 01 ,
     *-2.0086E-01 ,  1.3641E 00 , -3.4030E 00 ,  3.8559E 00 ,
     * 1.2583E-02 , -8.3492E-02 ,  1.8600E-01 , -2.0043E-01 ,
     *-2.3628E-04 ,  1.3514E-03 , -2.4324E-03 ,  2.1906E-03 ,
     * 1.2419E 00 , -4.3633E 00 ,  1.3743E 01 , -1.8592E 01 ,
     *-2.4404E-01 ,  1.3158E 00 , -3.5691E 00 ,  4.3867E 00 ,
     * 1.5693E-02 , -8.2579E-02 ,  2.1427E-01 , -2.5846E-01 ,
     *-2.9386E-04 ,  1.4060E-03 , -3.3835E-03 ,  3.8664E-03 /
      DATA  B2/
     * 6.3054E-01 , -3.7333E 00 ,  1.3464E 01 , -1.8594E 01 ,
     *2.1801E 00  ,  1.5163E 00 , -1.6380E 01 ,  2.7944E 01 ,
     *-1.2886E 00 , -2.4570E 00 ,  1.5129E 01 , -2.3295E 01 ,
     * 2.0915E-01 ,  5.2279E-01 , -2.8687E 00 ,  4.2688E 00 ,
     * 9.3363E-01 , -1.8181E 00 ,  5.5157E 00 , -8.5216E 00 ,
     * 1.7811E 00 , -8.2927E 00 ,  2.0607E 01 , -2.0827E 01 ,
     *-1.5264E 00 ,  6.8433E 00 , -1.6067E 01 ,  1.6845E 01 ,
     * 2.7128E-01 , -1.1944E 00 ,  2.7495E 00 , -2.9045E 00 ,
     * 1.9439E 00 , -4.6268E 00 ,  9.7879E 00 , -9.6074E 00 ,
     *-3.4640E-01 ,  1.1093E 00 , -1.9313E 00 ,  1.7064E 00 ,
     * 2.7054E-02 , -1.1638E-01 ,  2.6969E-01 , -3.1853E-01 ,
     *-6.6092E-04 ,  5.0728E-03 , -1.4995E-02 ,  1.9605E-02 ,
     * 1.8693E 00 , -5.5678E 00 ,  1.4795E 01 , -1.6903E 01 ,
     *-4.9965E-01 ,  1.7874E 00 , -4.1330E 00 ,  3.8393E 00 ,
     * 4.6194E-02 , -1.8536E-01 ,  4.5315E-01 , -4.6273E-01 ,
     *-1.3341E-03 ,  5.7710E-03 , -1.4554E-02 ,  1.5554E-02 /
      DATA  CKJ/
     * 1.4509E-01 ,  4.6520E-01 , -3.3005E-02 ,  1.5376E-01 ,
     * 2.7436E-01 , -1.4604E-02 ,  6.2959E-01 ,  1.7866E-01 ,
     *-2.6216E-03 ,  8.3810E-01 ,  8.6137E-03 ,  3.2946E-03 ,
     * 9.2852E-02 ,  5.3886E-01 , -5.4493E-02 ,  1.3032E-01 ,
     * 4.0709E-01 , -2.8782E-02 ,  1.4909E-01 ,  3.8502E-01 ,
     *-1.2775E-02 ,  1.8024E-01 ,  3.3022E-01 , -9.4491E-03 /
      END
C     * * * * * * * * * * * * * * * * *
      FUNCTION SIGMAG (L,MS,MQ,KSI,IKS,T)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 SIGMAG,T
C     SUBPROGRAM OF CHOOSING CROSS SECTION TYPE AND
C     CALCULATION CROSS SECTION VALUE FOR GIVEN ENERGY.
      COMMON /TYPECS/ICST,NSICST
      DIMENSION ICST(37),NSICST(37)
      ICS = 10000*L+1000*IABS(MS)+100*MQ+10*KSI+IKS
      IF(MS.LT.0)  ICS=-ICS
      JS = 1
   17 IF (ICS-ICST(JS)) 11,10,11
   10 SIGMAG = QINTG(T,JS)
      RETURN
   11 IF(JS-37)12,13,13
   12 JS = JS + 1
      GO TO 17
   13 NSJS = 1
   16 IF (ICS-NSICST(NSJS)) 15,14,15
   15 IF(NSJS-21)100,39,39
  100 NSJS=NSJS+1
      GO TO 16
   14 KNS = NSJS
      GO TO (18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,
     137,38),KNS
   18 SIGMAG = 0.
      RETURN
   19 SIGMAG = QINTG(T,10)
      RETURN
   20 SIGMAG = 0.
      RETURN
   21 SIGMAG = QINTG(T,15)+QINTG(T,16)
      RETURN
   22 SIGMAG = QINTG(T,18)+QINTG(T,19)+QINTG(T,17)
      RETURN
   23 SIGMAG = (QINTG(T,8)+QINTG(T,5))/2.
      RETURN
   24 SIGMAG = (QINTG(T,9)+QINTG(T,6)-QINTG(T,7))/2.
      RETURN
   25 SIGMAG = QINTG(T,7)
      RETURN
   26 SIGMAG = QINTG(T,10)*0.5
      RETURN
   27 SIGMAG = (QINTG(T,15)+QINTG(T,18))/2.
      RETURN
   28 SIGMAG = (QINTG(T,16)+QINTG(T,19))/2.
      RETURN
   29 SIGMAG = (QINTG(T,17))/2.
      RETURN
   30 SIGMAG = (QINTG(T,15)+QINTG(T,16)+QINTG(T,18)+
     1QINTG(T,19)+QINTG(T,17))/2.
      RETURN
   31 SIGMAG = 0.
      RETURN
   32 SIGMAG = 0.
      RETURN
   33 SIGMAG = 0.
      RETURN
   34 SIGMAG = QINTG(T,11)+QINTG(T,12)
      RETURN
   35 SIGMAG = 0.
      RETURN
   36 SIGMAG = 0.
      RETURN
   37 SIGMAG = QINTG(T,14)
      RETURN
   38 SIGMAG = QINTG(T,13)+2.*QINTG(T,14)
      RETURN
   39 DO  40  ISD=1,16
      KSD=ISD
      IF(ICS.EQ.NSICST(21+ISD))  GO  TO  41
   40 CONTINUE
      GO  TO  35
   41 CONTINUE
      W2=1.88*T+(1.232+0.940)**2
      TN=((DSQRT(W2)+0.94-1.232)**2)/1.88-1.88
      PNN2=W2/4.-0.940**2
      PDN2=(W2+0.940**2-1.232**2)**2/(4.*W2)-0.940**2
      R=0.5*PNN2/PDN2
      HEX=0.5
C  * * * * * * *  PRELIMINARY VERSION  * * * * * * * * * *
      GO  TO  (230,231,232,233,
     *         240,241,242,243,
     *         250,251,252,253,
     *         260,261,262,263),KSD
  230 SIGMAG=QINTG(TN,3)
      RETURN
  240 GO  TO  230
  250 SIGMAG=QINTG(TN,1)
      RETURN
  260 GO  TO  250
  231 SIGMAG=QINTG(TN,4)*(1.-HEX)
      RETURN
  241 GO  TO  231
  251 SIGMAG=QINTG(TN,2)*(1.-HEX)
      RETURN
  261 GO  TO  251
  232 SIGMAG=QINTG(TN,4)*HEX
      RETURN
  242 GO  TO  232
  252 SIGMAG=QINTG(TN,2)*HEX
      RETURN
  262 GO  TO  263
  233 SIGMAG=QINTG(TN,12)/2.*R
      RETURN
  243 SIGMAG=(QINTG(TN,12)/2.+QINTG(TN,11))*R
      RETURN
  253 SIGMAG=(QINTG(TN,14)/2.+QINTG(TN,13)/2.)*R
      RETURN
  263 SIGMAG=0.
      RETURN
       END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      FUNCTION  CROSEG(L,MS,MQ,KSI,IKS,T,AM,IR)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 CROSEG,T,AM
      CROSEG=0.
      X=T
      W=DSQRT(1.88*T+(1.88+AM)**2)
      IF(IR.EQ.0.AND.MS.EQ.0)      GO  TO  10
      IF(IR.NE.0.AND.MS.EQ.0)      GO  TO  11
      IF(IR.EQ.0.AND.MS.NE.0)      GO  TO  12
      IF(IR.NE.0.AND.MS.NE.0)      GO  TO  14
      RETURN
   10 CROSEG=SIGMAG(L,MS,MQ,KSI,IKS,X)
      GO  TO  15
   11 IF(MQ.EQ.2)                  GO  TO  10
      IF(IKS.GE.3)                 RETURN
      T0=((W+0.140-AM)**2-(0.940+0.140)**2)/1.88
      X=T0
      GO  TO  10
   12 IF(MQ.NE.1)                  GO  TO  13
      P0=DSQRT(T*(T+2.*AM))
      X=P0
      GO  TO  10
   13 IF(IKS.GT.1)                 RETURN
      TN=((W+0.940-AM)**2-(0.940+0.940)**2)/1.88
      TP=((W+0.140-AM)**2-(0.940+0.140)**2)/1.88
      TK=((W+0.492-AM)**2-(0.940+0.492)**2)/1.88
      PK=DSQRT(TK*(TK+2.*0.492))
      SPP  =QINTG(TN, 1+IKS)
      SPN  =QINTG(TN, 3+IKS)
      SPIPP=QINTG(TP, 8+IKS)
      SPIMP=QINTG(TP, 5+IKS)
      SKMP =QINTG(PK,33+IKS)
      SKMN =QINTG(PK,31+IKS)
      SLP=SPN+SKMP-SPIMP
      SLN=SPP+SKMN-SPIPP
      IF(MS.EQ.-1.AND.KSI.EQ.1)  CROSEG=   SLP+SPP-SPN
      IF(MS.EQ.-1.AND.KSI.EQ.2)  CROSEG=   SLP+SPN-SPP
      IF(MS.EQ.-1.AND.KSI.EQ.3)  CROSEG=   SLP
      IF(MS.EQ.-1.AND.KSI.EQ.4)  CROSEG=   SLN+SPN-SPP
      IF(MS.EQ.-1.AND.KSI.EQ.5)  CROSEG=   SLN+SPP-SPN
      IF(MS.EQ.-1.AND.KSI.EQ.6)  CROSEG=   SLN
      IF(MS.EQ.-2.AND.KSI.EQ.1)  CROSEG=2.*SLP-SPP
      IF(MS.EQ.-2.AND.KSI.EQ.2)  CROSEG=2.*SLP-SPN
      IF(MS.EQ.-2.AND.KSI.EQ.3)  CROSEG=2.*SLN-SPN
      IF(MS.EQ.-2.AND.KSI.EQ.4)  CROSEG=2.*SLN-SPP
      IF(MS.EQ.-3.AND.KSI.EQ.1)  CROSEG=3.*SLP-SPP-SPN
      IF(MS.EQ.-3.AND.KSI.EQ.2)  CROSEG=3.*SLN-SPP-SPN
      GO  TO  15
   14 IF(MQ.EQ.2)           GO  TO  13
      T0=((W+0.492-AM)**2-(0.940+0.492)**2)/1.88
      P0=DSQRT(T0*(T0+2.*0.492))
      X=P0
      GO  TO  10
   15 IF(CROSEG.LT.0.)  CROSEG=0.
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE ROTOR (AR,BR,PSTAR,PR)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 AR,BR,PSTAR,PR
C    BLOCK OF ROTATION.
      DIMENSION AR(3),BR(3),PSTAR(3),PR(3),AN(3)
      SP = 0.
      DO 31 IR=1,3
      SP = SP+AR(IR)*BR(IR)
   31 CONTINUE
      AMOD = DSQRT (AR(1)**2+AR(2)**2+AR(3)**2)
      IF(AMOD.LT.1.E-30)  GO  TO  10
      ALPHA1 = SP/AMOD
      BMOD2 = BR(1)**2+BR(2)**2+BR(3)**2
      IF(BMOD2.LT.(ALPHA1**2))  GO  TO  10
      ALPHA2 = DSQRT (BMOD2-ALPHA1**2)
      IF(ALPHA2.LT.1.E-30)  GO  TO  10
      AN(1) = AR(2)*BR(3)-AR(3)*BR(2)
      AN(2) = AR(3)*BR(1)-AR(1)*BR(3)
      AN(3) = AR(1)*BR(2)-AR(2)*BR(1)
      PR(1)=PSTAR(1)*BR(1)/ALPHA2+(PSTAR(3)-ALPHA1*PSTAR(1)/ALPHA2)
     1*AR(1)/AMOD+(PSTAR(2)*AN(1))/(ALPHA2*AMOD)
      PR(2)=PSTAR(1)*BR(2)/ALPHA2+(PSTAR(3)-ALPHA1*PSTAR(1)/ALPHA2)
     1*AR(2)/AMOD+(PSTAR(2)*AN(2))/(ALPHA2*AMOD)
      PR(3)=PSTAR(1)*BR(3)/ALPHA2+(PSTAR(3)-ALPHA1*PSTAR(1)/ALPHA2)
     1*AR(3)/AMOD+(PSTAR(2)*AN(3))/(ALPHA2*AMOD)
      RETURN
   10 DO  11  K=1,3
   11 PR(K)=PSTAR(K)
      RETURN
       END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE BINEL(PARTIN,IPATIN,IPATNE,L,MS,MQ,KSI,ME,V,U,
     *TIN1,MV,NP,NIN)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 PARTIN,V,U,TIN1
C     BLOCK OF INELASTIC SCATTERING.
      COMMON /I3ACT/ SIG3,SIGIN,IND,ITH
      COMMON/NCASCA/NCAS,NCPRI
      DIMENSION PARTIN(9),IPATIN(5),IPATNE(5),V(3)
      NIN = 0
      IK = 0
      AM3=PARTIN(9)
      IF(IPATIN(5).EQ.0)  GO  TO  11
      IF(IPATIN(4).EQ.1)  AM3=0.940
      GO  TO  14
   11 IF(TIN1-4.)12,24,24
   24 BETATH=0.
              GO TO 13
   12 SIG3=  CROSEG(L,MS,MQ,KSI,7,TIN1,PARTIN(9),IPATIN(5))
      SIGIN= CROSEG(L,MS,MQ,KSI,0,TIN1,PARTIN(9),IPATIN(5))-
     -       CROSEG(L,MS,MQ,KSI,1,TIN1,PARTIN(9),IPATIN(5))-
     -       CROSEG(L,MS,MQ,KSI,2,TIN1,PARTIN(9),IPATIN(5))
      BETATH=SIG3/SIGIN
   13 DRND=RNDMD(-1)
      IF(DRND-BETATH)14,15,15
   14 ITH=1
      TH=1.
                  GO TO 16
   15 ITH=0
      TH=2.
   16 IF(U- AM3     -0.14*TH-0.96)19,19,17
   17 CONTINUE
c-Sob      IF(NCAS.GE.NCPRI) write( *,101)
      CALL VMNSP (PARTIN,IPATIN,U,MV,NP,ITH,MQ,TIN1,LP)
      IF (NP) 26,26,27
   26 RETURN
   27 CONTINUE
      IF (LP) 19,22,19
   22 CONTINUE
c-Sob      IF(NCAS.GE.NCPRI) write( *,102)
      CALL DIRECT (V,TIN1,MQ,MV,NP,PARTIN,KP,ITH)
      IF (KP) 23,18,23
   18 CONTINUE
c-Sob      IF(NCAS.GE.NCPRI) write( *,103)
      CALL CHINEL (IPATIN,L,MS,MQ,KSI,NP,MV,TIN1,ME,IPATNE,PARTIN(9))
c-Sob      IF(IPATIN(5).EQ.0.AND.ITH.EQ.1.AND.NCAS.GE.NCPRI) write( *,104)
      IF(IPATIN(5).EQ.0.AND.ITH.EQ.1)
     *CALL  DISOB(MV,U,IND,NP,MQ)
      RETURN
   23 IK = IK+1
      IF(IK.LT.50)   GO  TO  17
   19 NIN=2
  101 FORMAT(' ====> VMNSP')
  102 FORMAT(' ====> DIRECT')
  103 FORMAT(' ====> CHINEL')
  104 FORMAT(' ====> DISOB')
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      FUNCTION  COSAL(T)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 COSAL,T
      IF(T-0.147) 1,1,2
    1 A1=25.2*T**(-0.843)
			      GO TO 3
    2 A1=130.*T**0.0145
    3 A2=11.3*T**0.432
      IF(T-0.055) 4,4,5
    4 A3=0.22*T**(-1.35)
			     GO TO 6
    5 A3=0.000043*T**(-4.32)
    6 A4=130.*T**1.33
      TEMP1=1.-DEXP(-A2*3.141592)
      TEMP3=2.-TEMP1
      TEMP2=1.-DEXP(-A4*3.141592)
      TEMP4=2.-TEMP2
      W=A3*TEMP4/((1.+A4**2)*(A1*TEMP3/(1.+A2**2)+A3*TEMP4/(1.+A4**2)))
      DRND=RNDMD(-1)
      IF(DRND-W) 7,7,9
    7 TAU=-DLOG(1.-RNDMD(-1)*TEMP2)/A4
      DRND=RNDMD(-1)
      IF(DRND-DSIN(TAU))8,8,7
    8 TETA=3.141592-TAU
			    GO TO 10
    9 TETA=-DLOG(1.-RNDMD(-1)*TEMP1)/A2
      DRND=RNDMD(-1)
      IF(DRND-DSIN(TETA)) 10,10,9
   10 COSAL= DCOS(TETA)
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      FUNCTION COSTA(J,T)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 COSTA,T
C     BLOCK OF COSINUS CALCULATION.
      COMMON/COEFA/ANKJ
      DIMENSION ANKJ(4,4,29),ANK(4,4)
      IF(T.LT.0.001)  GO  TO  14
      DO 10 K=1,4
      DO 10 N=1,4
      ANK(N,K) = ANKJ(N,K,J)
   10 CONTINUE
      S1 = 0.
      R1 = RNDMD(-1)
      S2 = 0.
      DO 11 N=1,4
      DO 11 K=1,4
      S1 = S1+ANK(N,K)*(T**(K-1))*(R1**(N-1))
   11 CONTINUE
      DO 12 N=1,4
      DO 12 K=1,4
      S2 = S2+ANK(N,K)*T**(K-1)
   12 CONTINUE
      CTA = 2.*DSQRT(R1)*(S1+(1.-S2)*R1**4)-1.
      TEMP1 = DABS(CTA)
      IF (TEMP1-1.) 13,13,14
   13 COSTA = CTA
                RETURN
   14 COSTA=1.-2.*RNDMD(-1)
      RETURN
       END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      FUNCTION QINTG (X,LQ)
      IMPLICIT REAL*8 (A-H,O-Z)
C      SUBPROGRAM OF INTERPOLATION OF CROSS SECTIONS .
      REAL*8 QINTG,X
      COMMON /TABELE/SIGMA,ARGUS
      DIMENSION SIGMA(30,37),ARGUS(30,6)
      IF(LQ-4) 19,19,20
   19 I=1
          GO TO 27
   20 IF(LQ-10) 21,21,22
   21 I=2
          GO TO 27
   22 IF(LQ-19) 23,23,24
   23 I=3
          GO TO 27
   24 IF(LQ-23) 25,25,26
   25 I=4
          GO TO 27
   26 I=5
   27 CONTINUE
      LPHA = 1
   18 IF(X-ARGUS(LPHA,I)) 10,11,12
   10 IF(I-4) 100,101,101
  101 IF(I-5) 102,100,100
  102 IF(LPHA-1) 103,103,100
  103 PHI1=300.
      PSI1=0.005
      PHI2=SIGMA(1,LQ)
      PSI2=ARGUS(1,I)
      PHI3=SIGMA(2,LQ)
      PSI3=ARGUS(2,I)
      GO TO 15
  100 IF(LPHA-1) 13,13,14
   13 QINTG = 0.
               RETURN
   14 IF (LPHA - 29) 17,16,16
   16 PHI1=SIGMA(28,LQ)
      PSI1=ARGUS(28,I)
      PHI2=SIGMA(29,LQ)
      PSI2=ARGUS(29,I)
      PHI3=SIGMA(30,LQ)
      PSI3=ARGUS(30,I)
      GO TO 15
   17 PHI1 = SIGMA (LPHA-1,LQ)
      PSI1 = ARGUS (LPHA-1,I)
      PHI2 = SIGMA (LPHA,LQ)
      PSI2 = ARGUS (LPHA,I)
      PHI3 = SIGMA (LPHA+1,LQ)
      PSI3 = ARGUS (LPHA+1,I)
      GO TO 15
   11 QINTG = SIGMA (LPHA,LQ)
      RETURN
   12 LPHA = LPHA + 1
      IF(LPHA.GT.30)  GO  TO  112
      GO TO 18
   15 DELTA=(PSI2-PSI3)*PSI1**2+(PSI3-PSI1)*PSI2**2+(PSI1-PSI2)*PSI3**2
      DELTAA=PHI1*(PSI2-PSI3)+PHI2*(PSI3-PSI1)+PHI3*(PSI1-PSI2)
      DELTAB=(PHI2-PHI3)*PSI1**2+(PHI3-PHI1)*PSI2**2+(PHI1-PHI2)*PSI3**2
     0DELTAC=(PSI2*PHI3-PSI3*PHI2)*PSI1**2+(PSI3*PHI1-PSI1*PHI3)*PSI2**2
     1+(PSI1*PHI2-PSI2*PHI1)*PSI3**2
      A = DELTAA / DELTA
      B = DELTAB / DELTA
      C = DELTAC / DELTA
      QINTG = A*X**2 + B*X + C
      RETURN
  112 QINTG= SIGMA(30,LQ)/(1.+0.0102*DLOG(1.13*ARGUS(30,I))**2)*
     *                    (1.+0.0102*DLOG(1.13*X          )**2)
      RETURN
       END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE SLQEK(L,MS,MQ,KSI,ME,IIN,IPN)
C     BLOCK OF FORMING CROSS SECTIONS TYPE.
      DIMENSION IIN(5),IPN(5)
      MEIN=IIN(1)
      MS = IIN(3)+IPN(3)
      L =  IIN(2)+IPN(2)
      MQ = IIN(4)+IPN(4)
      ME = IIN(1)+IPN(1)
      IF(IIN(5).NE.0.AND.IIN(4).EQ.1.AND.IIN(3).EQ.0)  GO  TO  30
    9 CONTINUE
      IF (MS) 10,11,10
   10 IF(IIN(4).NE.0)    GO  TO  110
      IF(MS.EQ. 1.AND.IIN(1).EQ. 1.AND.IPN(1).EQ.1)    KSI=1
      IF(MS.EQ. 1.AND.IIN(1).EQ. 0.AND.IPN(1).EQ.1)    KSI=2
      IF(MS.EQ.-1.AND.IIN(1).EQ.-1.AND.IPN(1).EQ.1)    KSI=2
      IF(MS.EQ.-1.AND.IIN(1).EQ. 0.AND.IPN(1).EQ.1)    KSI=1
      IF(MS.EQ. 1.AND.IIN(1).EQ. 1.AND.IPN(1).EQ.0)    KSI=2
      IF(MS.EQ. 1.AND.IIN(1).EQ. 0.AND.IPN(1).EQ.0)    KSI=1
      IF(MS.EQ.-1.AND.IIN(1).EQ.-1.AND.IPN(1).EQ.0)    KSI=1
      IF(MS.EQ.-1.AND.IIN(1).EQ. 0.AND.IPN(1).EQ.0)    KSI=2
      RETURN
  110 IF(MS.EQ.-1.AND.IIN(1).EQ. 1.AND.IPN(1).EQ.1)    KSI=1
      IF(MS.EQ.-1.AND.IIN(1).EQ. 1.AND.IPN(1).EQ.0)    KSI=4
      IF(MS.EQ.-1.AND.IIN(1).EQ.-1.AND.IPN(1).EQ.1)    KSI=2
      IF(MS.EQ.-1.AND.IIN(1).EQ.-1.AND.IPN(1).EQ.0)    KSI=5
      IF(MS.EQ.-1.AND.IIN(1).EQ. 0.AND.IPN(1).EQ.1)    KSI=3
      IF(MS.EQ.-1.AND.IIN(1).EQ. 0.AND.IPN(1).EQ.0)    KSI=6
      IF(MS.EQ.-2.AND.IIN(1).EQ.-1.AND.IPN(1).EQ.1)    KSI=1
      IF(MS.EQ.-2.AND.IIN(1).EQ.-1.AND.IPN(1).EQ.0)    KSI=3
      IF(MS.EQ.-2.AND.IIN(1).EQ. 0.AND.IPN(1).EQ.1)    KSI=2
      IF(MS.EQ.-2.AND.IIN(1).EQ. 0.AND.IPN(1).EQ.0)    KSI=4
      IF(MS.EQ.-3.AND.IIN(1).EQ.-1.AND.IPN(1).EQ.1)    KSI=1
      IF(MS.EQ.-3.AND.IIN(1).EQ.-1.AND.IPN(1).EQ.0)    KSI=2
      RETURN
   11 IF (L) 13,13,12
   12 KSI = 1
      RETURN
   13 IF(MQ-1) 15,15,100
  100 IF(MQ-2) 14,14,16
   14 IF (ME-1)16,17,16
   16 KSI = 1
      RETURN
   17 KSI = 2
      RETURN
   15 IF (ME-2) 19,18,19
   18 KSI = 1
      RETURN
   19 IF (ME+1) 21,20,21
   20 KSI = 1
      RETURN
   21 IF (ME) 22,23,22
   22 IF (MEIN-1) 27,26,27
   26 KSI = 2
      RETURN
   27 KSI = 3
      RETURN
   23 IF (MEIN+1) 25,24,25
   24 KSI = 2
      RETURN
   25 KSI = 3
      RETURN
   30 KSI=5+IIN(1)*(IPN(1)-1)+IPN(1)*(IIN(1)-1)
      RETURN
       END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE CINEMA(PSTAR,V,P,CT,ST,CFI,SFI,T,CM)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 PSTAR,V,P,CT,ST,CFI,SFI,T,CM
C     KINEMATIC BLOCK.
      DIMENSION PSTAR(3),V(3),P(3)
      SPV = PSTAR(1)*V(1)+PSTAR(2)*V(2)+PSTAR(3)*V(3)
      G=1./DSQRT(1.d0-V(1)**2-V(2)**2-V(3)**2)
      ESTAR=DSQRT(PSTAR(1)**2+PSTAR(2)**2+PSTAR(3)**2+CM**2)
      DO  K=1,3
      P(K)=PSTAR(K)+G*V(K)*(SPV*G/(G+1.)+ESTAR)
      ENDDO
      PM = DSQRT(P(1)**2+P(2)**2+P(3)**2)
      IF(PM.LT.1.E-30)   THEN
      CT=1.d0
      CFI=1.d0
      SFI=0.d0
      ST=0.d0
      ELSE
      CT = P(3)/PM
      TEMP4 = 1.d0-CT**2
      if(TEMP4.lt.1.E-30) then
        CT=1.d0
        CFI=1.d0
        SFI=0.d0
        ST=0.d0
      else
        ST = DSQRT(TEMP4)
        TEMP3 = PM*ST
        CFI=P(1)/TEMP3
        SFI=P(2)/TEMP3
      endif
      ENDIF
      T=DSQRT(PM**2+CM**2)-CM
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      FUNCTION PMOM (J,T)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 PMOM,T
C     BLOCK OF CALCULATION OF SECONDARY PARTICLES MOMENTUM.
      COMMON/COEFBC/BNKJ,CKJ
      DIMENSION BNKJ(4,4,8),CKJ(3,8),BNK(4,4)
      DO 10 K=1,4
      DO 10 N=1,4
      BNK(N,K) = BNKJ(N,K,J)
   10 CONTINUE
      S1 = 0.
      R1 = RNDMD(-1)
      S2 = 0.
      S3 = 0.
      DO 11 N=1,4
      DO 11 K=1,4
      S1 = S1+BNK(N,K)*(T**(K-1))*(R1**(N-1))
   11 CONTINUE
      DO 12 N=1,4
      DO 12 K=1,4
      S2 = S2+BNK(N,K)*T**(K-1)
   12 CONTINUE
      DO 13 K=1,3
      S3 = S3+CKJ(K,J)*T**(K-1)
   13 CONTINUE
      PMOM = S3*DSQRT(R1)*(S1+(1.-S2)*R1**4)
      RETURN
       END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      FUNCTION JTYPA (ITH,MQ,LAMB)
C     DETERMINING OF TYPE OF COEFFICIENTS A(N,K).
      IF (ITH) 10,18,10
   10 IF (MQ-1) 11,11,14
   11 IF (LAMB-1) 13,13,12
   12 JTYPA = 19
               RETURN
   13 JTYPA = 18
               RETURN
   14 IF (LAMB-1) 16,16,15
   15 IF (LAMB-3) 17,16,17
   16 JTYPA = 14
               RETURN
   17 JTYPA = 15
               RETURN
   18 IF (MQ-1) 19,19,22
   19 IF (LAMB-1) 20,20,21
   20 JTYPA = 20
               RETURN
   21 JTYPA = 21
               RETURN
   22 IF (LAMB-1) 24,24,23
   23 IF (LAMB-3) 25,24,25
   24 JTYPA = 16
               RETURN
   25 JTYPA = 17
               RETURN
       END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      FUNCTION JTYPB (ITH,MQ,LAMB)
C     DETERMINING OF TYPE OF COEFFICIENTS B(N,K).
      IF (ITH) 10,18,10
   10 IF (MQ-1) 11,11,14
   11 IF (LAMB-1) 13,13,12
   12 JTYPB = 6
              RETURN
   13 JTYPB = 5
              RETURN
   14 IF (LAMB-1) 16,16,15
   15 IF (LAMB-3) 17,16,17
   16 JTYPB = 1
              RETURN
   17 JTYPB = 2
              RETURN
   18 IF (MQ-1) 19,19,22
   19 IF (LAMB-1) 20,20,21
   20 JTYPB = 7
              RETURN
   21 JTYPB = 8
              RETURN
   22 IF (LAMB-1) 24,24,23
   23 IF (LAMB-3) 25,24,25
   24 JTYPB = 3
              RETURN
   25 JTYPB = 4
              RETURN
       END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE COTRAN
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/DTINT/DTAU
      COMMON/NUCSP/VPR(3),VTA(3),RADP(3),RADT(3),VEV(3),VRE(3),GEV,GRE
     *,VEP(3),VET(3),GEP,GET
      COMMON/RESULT/AN1,AN2,ZN1,ZN2,ENEXT1,ENEXT2,PNUCL1(3),PNUCL2(3),
     *AMNUC1(3),AMNUC2(3)
      DIMENSION VCM(3),P(3),P1(3),P2(3),R(3),B(3),RN(3)
      IF((AN1+AN2).EQ.0.)   RETURN
      AMR=AN1*AN2/(AN1+AN2)*0.940
      IF(AMR.LE.0.)  RETURN
      IF(DTAU.LT.0.) RETURN
      AMP=AN1*.940
      AMT=AN2*.940
      GP=1./DSQRT(1.-VPR(1)**2-VPR(2)**2-VPR(3)**2)
      GT=1./DSQRT(1.-VTA(1)**2-VTA(2)**2-VTA(3)**2)
      G1=GP*AN1/(GP*AN1+GT*AN2)
      G2=GT*AN2/(GP*AN1+GT*AN2)
      DO 10  K=1,3
      RN(K)=RADP(K)-RADT(K)
      R(K)=RN(K)-(VPR(K)-VTA(K))*DTAU
      VCM(K)=G1*VPR(K)+G2*VTA(K)
      B(K)=-VCM(K)
   10 P1(K)=GP*AMP*VPR(K)
      CALL CINEMA(P1,B,P,CT,ST,CF,SF,T,AMP)
      P0=DSQRT(P(1)**2+P(2)**2+P(3)**2)
      R0=DSQRT(R(1)**2+R(2)**2+R(3)**2)
      R1=DSQRT(RN(1)**2+RN(2)**2+RN(3)**2)
      PRN=P(1)*RN(1)+P(2)*RN(2)+P(3)*RN(3)
      PL2=(PRN/R1)**2-2.*AMR*(POTNU(R1)-POTNU(R0))
      SPL=1.
      IF(PRN.LT.0.)   SPL=-1.
      IF(PL2.LT.0.)   PL2=0.
      TEMP=(SPL*DSQRT(PL2)-PRN/R1)/R1
      DO  11  K=1,3
   11 P(K)=P(K)+TEMP*RN(K)
      CALL  CINEMA(P,VCM,P1,CT1,ST1,CF1,SF1,T1,AMP)
      P(1)=-P(1)
      P(2)=-P(2)
      P(3)=-P(3)
      CALL  CINEMA(P,VCM,P2,CT2,ST2,CF2,SF2,T2,AMT)
      EP=AMP+T1
      ET=AMT+T2
      DO  12  K=1,3
      VPR(K)=P1(K)/EP
   12 VTA(K)=P2(K)/ET
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      FUNCTION  FORCEN (R1,R2)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 FORCEN,R1,R2
      IF(DABS(R1-R2).EQ.0.)  GO  TO 10
      FORCEN=-(POTNU(R2)-POTNU(R1))/(R2-R1)
      RETURN
   10 FORCEN=0.
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      FUNCTION  POTNU(R)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 POTNU,R
      COMMON/HCASC/ANUCL1,ANUCL2,ZNUCL1,ZNUCL2,T0,EPS1,EPS2,
     *VPI,A1,A2,C1,C2,D1,D2,R0N1,R0N2,TF01,TF02,RM1,RM2
      COMMON/RESULT/AN1,AN2,ZN1,ZN2,ENEXT1,ENEXT2,PNUCL1(3),PNUCL2(3),
     *AMNUC1(3),AMNUC2(3)
      R1=0.
      IF(AN1.GT.0.1)  S1=R0N1*AN1**(1./3.)
      IF(AN1.GT.0.1)  R1=S1*(1.+(C1*3.141592/S1)**2)**(1./3.)
      R2=0.
      IF(AN2.GT.0.1)  S2=R0N2*AN2**(1./3.)
      IF(AN2.GT.0.1)   R2=S2*(1.+(C2*3.141592/S2)**2)**(1./3.)
      IF((R1+R2).LE.0.)  GO  TO  14
      RP=R1
      IF(R1.GT.R2)  RP=R2
      RT=R2
      IF(R1.GE.R2)  RT=R1
      IF(RP.LE.0.0)  GO  TO  14
      X=R/(RP+RT)
      B=0.00144*ZN1*ZN2/(RP+RT)
      A=RT/RP
      BM=1.+1./A
      C=3./5./A/A
      D=(A+2.+1./A)/4.
      X0=(A-1.)/(A+1.)
      IF(X-X0)  10,10,11
   10 POTNU=B*(3.-C-(BM*X)**2)*BM/2.
      RETURN
   11 IF(X-1.) 12,12,13
   12 POTNU=B*(1.-3.*D*D*((1.-X)**4)*(1.-2./15.*D*(1.-X)*(5.+X)))/X
      RETURN
   13 POTNU=B/X
      RETURN
   14 POTNU=0.
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE  KINEMA(PS,V,P,CT,ST,CF,SF,T,CM)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 PS,V,P,CT,ST,CF,SF,T,CM
C   KINEMATIC  BLOCK
      DIMENSION  PS(3),V(3),P(3)
      PSV=PS(1)*V(1)+PS(2)*V(2)+PS(3)*V(3)
      ES=DSQRT(PS(1)**2+PS(2)**2+PS(3)**2+CM**2)
      G=1./DSQRT(1.-V(1)**2-V(2)**2-V(3)**2)
      DO  10  K=1,3
   10 P(K)=PS(K)+G*V(K)*(PSV*G/(G+1.)+ES)
      E=G*(ES+PSV)
      T=E-CM
      PM=DSQRT(P(1)**2+P(2)**2+P(3)**2)
      CT=P(3)/PM
      ST2=1.-CT*CT
      IF(ST2.LE.0.)   GO  TO  11
      ST=DSQRT(ST2)
      CF=P(1)/PM/ST
      SF=P(2)/PM/ST
      RETURN
   11 ST=0.
      CF=1.
      SF=0.
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
       SUBROUTINE  COLEC(NCAS,INTEL,MV,NZON,LUP,ICOL)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*4 UP1,UP2
      REAL*8 MMES,MBAR
      COMMON/HCASC/ANUCL1,ANUCL2,ZNUCL1,ZNUCL2,T0,EPS1,EPS2,VPI,
     *A1,A2,C1,C2,D1,D2,R0N1,R0N2,TF01,TF02,RM1,RM2
      COMMON/RESULT/AN1,AN2,ZN1,ZN2,EN1,EN2,PN1(3),PN2(3),AM1(3),AM2(3)
      COMMON/BIMP/BIMP,BIMX,BIMY/KYPT/KPT,KYP,KYT,KPT1,KYP1,KYT1,KPZ,KTZ
      COMMON/ACTIM/TINT
      COMMON/PIDABS/PIDABS
      COMMON/RCOR/ RCOR
      COMMON/MMATUR/MMES,MBAR
      COMMON /STIN/ STIN,AMIN
      COMMON/INDINT/ INDINT
      COMMON/INTCC/ INTCC
      COMMON /IACT/ IACT/INDDEC/INDDEC/ISOB3/ISOB3/ISOB2/ISOB2
      COMMON/KSYST/ KSYST
      COMMON/EXCIT/TEX1,TEX2,HEX1,HEX2,PEX1,PEX2
      COMMON/UPAC/UP1(66000),IPER/UPAC2/UP2(70000),LU2
      COMMON /NUCSP/VPR(3),VTA(3),RADP(3),RADT(3),VEV(3),VRE(3),GEV,GRE,
     *VEP(3),VET(3),GEP,GET
      COMMON/INDZAP/INDZAP
      COMMON/BGLAB/BLAB,GLAB,KOBR
      INDZAP=0
      IF(ICOL.EQ.1)   GO  TO  604
C
      IF((LU2+50+MV*11+1).GE.70000) GO  TO  604
C
  600 UP2(LU2)=DFLOAT(MV)
      UP2(LU2+1)=AN1
      UP2(LU2+2)=AN2
      UP2(LU2+3)=ZN1
      UP2(LU2+4)=ZN2
      UP2(LU2+5)=EN1
      UP2(LU2+6)=EN2
      DO  601  KL=1,3
      UP2(LU2+6+KL)=PN1(KL)
      UP2(LU2+9+KL)=PN2(KL)
      UP2(LU2+12+KL)=AM1(KL)
  601 UP2(LU2+15+KL)=AM2(KL)
      IF(KOBR.NE.1)  GO  TO  9
      IF(AN1.LT.0.1)  GO  TO  8
      EA1=DSQRT(PN1(1)**2+PN1(2)**2+PN1(3)**2+(0.94*AN1)**2)
      UP2(LU2+9)=-GLAB*(PN1(3)-BLAB*EA1)
    8 IF(AN2.LT.0.1)  GO  TO  9
      EA2=DSQRT(PN2(1)**2+PN2(2)**2+PN2(3)**2+(0.940*AN2)**2)
      UP2(LU2+12)=-GLAB*(PN2(3)-BLAB*EA2)
    9 CONTINUE
      UP2(LU2+19)=TEX1
      UP2(LU2+20)=TEX2
      UP2(LU2+21)=HEX1
      UP2(LU2+22)=HEX2
      UP2(LU2+23)=PEX1
      UP2(LU2+24)=PEX2
      UP2(LU2+25)=FLOAT(KPT)
      UP2(LU2+26)=FLOAT(KYP)
      UP2(LU2+27)=FLOAT(KYT)
      UP2(LU2+28)=FLOAT(KPT1)
      UP2(LU2+29)=FLOAT(KYP1)
      UP2(LU2+30)=FLOAT(KYT1)
      UP2(LU2+31)=FLOAT(KPZ)
      UP2(LU2+32)=FLOAT(KTZ)
      UP2(LU2+33)=BIMP
      UP2(LU2+34)=BIMX
      UP2(LU2+35)=BIMY
      UP2(LU2+36)=TINT
      DO  10  LB=1,3
      UP2(LU2+36+LB)=RADP(LB)
   10 UP2(LU2+39+LB)=RADT(LB)
      UP2(LU2+43)=PIDABS
      UP2(LU2+44)=RCOR
      UP2(LU2+45)=MMES
      UP2(LU2+46)=MBAR
      UP2(LU2+47)=STIN
      UP2(LU2+48)=AMIN
      UP2(LU2+49)=DFLOAT(10*INDINT+INTCC)
      UP2(LU2+50)=DFLOAT(KSYST*10000+IACT*1000+INDDEC*100+ISOB2*10+
     *ISOB3)
      MVK=11*MV
      IF(MV.EQ.0)   GO  TO  603
      DO  602  KL=1,MVK
      UP2(LU2+50+KL)=UP1(KL)
  602 CONTINUE
  603 LU2=LU2+50+MVK+1
      GO  TO  606
  604 IF(LU2.EQ.1)   GO  TO  606
      NCAS1=NCAS-(1-ICOL)
      CALL  WTAP10(NCAS1,INTEL,NZON)
      INDZAP=1
      IF(ICOL.EQ.0)  GO  TO  600
  606 CONTINUE
      LUP=LU2
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE RXYZ(R12,R0X,R0Y,R0Z)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 R12,R0X,R0Y,R0Z
      COMMON /XBMAX/ XBMAX,IFIB0
      COMMON /BIMP/ B00,BX,BY
      COMMON /HCASC/ANUCL1,ANUCL2,ZNUCL1,ZNUCL2,T0,EPS1,EPS2,VPI,A1,A2,
     *C1,C2,D1,D2,R0N1,R0N2,TF01,TF02,RM1,RM2
      DIMENSION VPR(3),VTA(3)
      IF(IFIB0.EQ.0)  GO TO 10
      B1=RNDMD(-1)
      B2=6.283185*RNDMD(-1)
      B3=R12*DSQRT(B1)*XBMAX
      R0X=B3*DCOS(B2)
      R0Y=B3*DSIN(B2)
      GO  TO  11
   10 R0X=R12*XBMAX
      R0Y=0.
      B3=R0X
   11 R0Z=-DSQRT(R12**2-R0X**2-R0Y**2)
      B00=B3
      BX=R0X
      BY=R0Y
      CALL  VINIT(VPR,VTA,ANUCL1,ANUCL2,T0)
      GPR=1./DSQRT(1.-VPR(3)**2)
      GTA=1./DSQRT(1.-VTA(3)**2)
      R0Z=-RM1/GPR-RM2/GTA
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE VINIT(VPR,VTA,A1,A2,T0)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 VPR,VTA,A1,A2,T0
      DIMENSION VPR(3),VTA(3)
      COMMON /KSYST/KSYST
      VPR(1)=0.
      VPR(2)=0.
       VTA(1)=0.
       VTA(2)=0.
      IF(A1.LT.2.1)  GO TO 10
      GO  TO  (10,11,12), KSYST
   10 VPR(3)=DSQRT(T0*(T0+1.88))/(T0+0.94)
       VTA(3)=1.E-6
      RETURN
   11 VPR(3)=DSQRT(T0*(T0+1.88))/(T0+1.88)
      VTA(3)=-VPR(3)
      RETURN
   12 VPR(3)= DSQRT(T0*(T0+1.88))/(T0+0.94*(1.+A1/A2))
      VTA(3)=-DSQRT(T0*(T0+1.88))/(T0+0.94*(1.+A2/A1))
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE SHOCKF(T1,T2,NA1,NA2,MV,IWR)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 T1,T2
      REAL*4     UP2
      COMMON /HCASC/ANUCL1,ANUCL2,ZNUCL1,ZNUCL2,T0,EPS1,EPS2,VPI,A1,A2,
     *C1,C2,D1,D2,R0N1,R0N2,TF01,TF02,RM1,RM2
      COMMON /MEMAGT/  PME(9,5999),  IME(5,5999)
     *       /CENTER/XC(2,240),YC(2,240),ZC(2,240),IZ(2,240)
     *       /KYPT/KPT,KYP,KYT,KPT1,KYP1,KYT1,KPZ,KTZ
      COMMON /RESULT/AN1,AN2,ZN1,ZN2,ENEXT1,ENEXT2,PN1(3),PN2(3),AM1(3),
     *AM2(3)
      COMMON /NUCSP/VPR(3),VTA(3),RADP(3),RADT(3),VEV(3),VRE(3),GEV,GRE,
     *VEP(3),VET(3),GEP,GET
      COMMON /BIMP/ BIMP,BIMX,BIMY
      COMMON /EXCIT/ TEX1,TEX2,HEX1,HEX2,PEX1,PEX2
      COMMON /TIMSH/ TIMF(15),NZON,KTW
      COMMON /UPAC2/ UP2(70000),LUP
      COMMON /TPROD/TPROD(5999)
      COMMON/INDZAP/INDZAP
      COMMON/IDPME/IDPME(5999)
      COMMON/PORIG/IORI(3,5999)
      COMMON/SORI/SORI(5999),SSOR
c
      INDZAP=0
      IF(IWR.EQ.1)  RETURN
      DO 20 KT=1,15
      IF (TIMF(KT).LT.T1/3.)   GO TO 20
      TF=T1
      IF (T2.LT.T1)  GO TO 10
      TF=T2
      IF (TIMF(KT).GT.T2/3.)   GO TO 20
   10 CONTINUE
      IR1=50+11*MV+(NA1+NA2)*4
      KTD=KT
      if((LUP+IR1+1).GT.70000)  then
c-Sob        write(16,'(2X,''MASSIVE UP2 IS EXCEEDED KT ='',I5)') KTD
      RETURN
      endif
      UP2(LUP)=DFLOAT(MV)
      UP2(LUP+1)=AN1
      UP2(LUP+2)=AN2
      UP2(LUP+3)=ZN1
      UP2(LUP+4)=ZN2
      UP2(LUP+5)=ENEXT1
      UP2(LUP+6)=ENEXT2
      GPR=1./DSQRT(1.-VPR(1)**2-VPR(2)**2-VPR(3)**2)
      GTA=1./DSQRT(1.-VTA(1)**2-VTA(2)**2-VTA(3)**2)
      DO KL=1,3
      UP2(LUP+6+KL)=GPR*0.940*AN1*VPR(KL)
      UP2(LUP+9+KL)=GTA*0.940*AN2*VTA(KL)
      UP2(LUP+12+KL)=AM1(KL)
      UP2(LUP+15+KL)=AM2(KL)
      ENDDO
      UP2(LUP+19)=TEX1
      UP2(LUP+20)=TEX2
      UP2(LUP+21)=HEX1
      UP2(LUP+22)=HEX2
      UP2(LUP+23)=PEX1
      UP2(LUP+24)=PEX2
      UP2(LUP+25)=DFLOAT(KPT)
      UP2(LUP+26)=DFLOAT(KYP)
      UP2(LUP+27)=DFLOAT(KYT)
      UP2(LUP+28)=DFLOAT(KPT1)
      UP2(LUP+29)=DFLOAT(KYP1)
      UP2(LUP+30)=DFLOAT(KYT1)
      UP2(LUP+31)=DFLOAT(KPZ)
      UP2(LUP+32)=DFLOAT(KTZ)
      UP2(LUP+33)=BIMP
      UP2(LUP+34)=BIMX
      UP2(LUP+35)=BIMY
      UP2(LUP+36)=TIMF(KT)
      DO KL=1,3
      UP2(LUP+36+KL)=RADP(KL)
      UP2(LUP+39+KL)=RADT(KL)
      ENDDO
      KTW=KT
      UP2(LUP+43)=FLOAT(KTW)
      DO KL=44,50
      UP2(LUP+KL)=0.
      ENDDO
      IF(MV.EQ.0)  GO TO 15
      DR=3.*TIMF(KT)-TF
      DO M=1,MV
      KM=LUP+50+(M-1)*11
      EK=PME(8,M)+PME(9,M)
      UP2(KM+ 1)=DFLOAT(IDPME(M))
      UP2(KM+ 2)=DFLOAT(IORI(1,M))
        ICODE=(IORI(3,M)+4)*10**4 + IABS(IORI(2,M))
        ICODE=ISIGN(ICODE,IORI(2,M))
      UP2(KM+ 3)=DFLOAT(ICODE)
      UP2(KM+ 4)=PME(9,M)
      UP2(KM+ 5)=PME(1,M)+DR*PME(4,M)/EK
      UP2(KM+ 6)=PME(2,M)+DR*PME(5,M)/EK
      UP2(KM+ 7)=PME(3,M)+DR*PME(6,M)/EK
      UP2(KM+ 8)=PME(4,M)
      UP2(KM+ 9)=PME(5,M)
      UP2(KM+10)=PME(6,M)
c        UP2(KM+11)=TPROD(M)   ! production time
c        UP2(KM+11)=SORI(M)    ! collision energy
      UP2(KM+11)=TPROD(M)+PME(7,M) ! with muturity
      ENDDO
   15 IF(NA1.EQ.0)  GO TO 17
      GPR=1./DSQRT(1.-VPR(1)**2-VPR(2)**2-VPR(3)**2)
      DRPX=RADP(1)+VPR(1)*(3.*TIMF(KT)-TF)
      DRPY=RADP(2)+VPR(2)*(3.*TIMF(KT)-TF)
      DRPZ=RADP(3)+VPR(3)*(3.*TIMF(KT)-TF)
      DO I=1,NA1
      VPK=(XC(1,I)*VPR(1)+YC(1,I)*VPR(2)+
     +                      ZC(1,I)*VPR(3))*GPR/(GPR+1.)
      KM=LUP+50+MV*11+(I-1)*4
      UP2(KM+1)=DFLOAT(4+IZ(1,I))
      UP2(KM+2)=XC(1,I)-VPR(1)*VPK+DRPX
      UP2(KM+3)=YC(1,I)-VPR(2)*VPK+DRPY
      UP2(KM+4)=ZC(1,I)-VPR(3)*VPK+DRPZ
      ENDDO
   17 IF(NA2.EQ.0)  GO TO 19
      GTA=1./DSQRT(1.-VTA(1)**2-VTA(2)**2-VTA(3)**2)
      DRTX=RADT(1)+VTA(1)*(3.*TIMF(KT)-TF)
      DRTY=RADT(2)+VTA(2)*(3.*TIMF(KT)-TF)
      DRTZ=RADT(3)+VTA(3)*(3.*TIMF(KT)-TF)
      DO I=1,NA2
      VTK=(XC(2,I)*VTA(1)+YC(2,I)*VTA(2)+
     &                      ZC(2,I)*VTA(3))*GTA/(GTA+1.)
      KM=LUP+50+MV*11+NA1*4+(I-1)*4
      UP2(KM+1)=DFLOAT(4+IZ(2,I))
      UP2(KM+2)=XC(2,I)-VTA(1)*VTK+DRTX
      UP2(KM+3)=YC(2,I)-VTA(2)*VTK+DRTY
      UP2(KM+4)=ZC(2,I)-VTA(3)*VTK+DRTZ
      ENDDO
   19 LUP=LUP+50+MV*11+(NA1+NA2)*4+1
      CALL  WTAP10(NA1,NA2,NZON)
      IF(KTW.EQ.15)  INDZAP=1
   20 CONTINUE
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE  DISOB(MV,U,IND,NP,MQ)
c new version 11-16-95 09:59am
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 U
      REAL*8  MD,M2
      COMMON /PRIMP/ PP(3)
      COMMON /TAUD3/ TAU0
      COMMON /ISOB3/ ISOB3
      COMMON /MEMAGT/ PME(9,5999),IME(5,5999)
      IND=0
      IF(ISOB3.NE.1)  RETURN
      IATT=0
      RND=RNDMD(-1)
      IF(MQ.EQ.2) THEN
      N1=MV+1
      N2=MV+3
      ELSE
      N1=MV+2
      N2=MV+3
      ENDIF
    1 CONTINUE
      IATT=IATT+1
      IF(IATT.GT.2)        RETURN
      IF((RND.LE.0.5.AND.IATT.EQ.1).OR.IATT.EQ.2)  THEN
      DO  K=1,9
       TEMP=PME(K,N1)
       PME(K,N1)=PME(K,N2)
       PME(K,N2)=TEMP
      ENDDO
      DO  K=1,5
       ITEMP=IME(K,N1)
       IME(K,N1)=IME(K,N2)
       IME(K,N2)=ITEMP
      ENDDO
      ENDIF
      L2=MV+3-MQ
      LD=MV+MQ
      IF((IME(4,LD)+IME(4,MV+3)).NE.1)  GO  TO  1
      P2=DSQRT(PME(4,L2)**2+PME(5,L2)**2+PME(6,L2)**2)
      PP(1)=PME(4,MV+2)
      PP(2)=PME(5,MV+2)
      PP(3)=PME(6,MV+2)
      M2=PME(9,L2)
      E2=DSQRT(P2**2+M2**2)
      MD=DSQRT(U*(U-2.*E2)+M2**2)
      IF(MD.LT.1.082)                   GO  TO  1
      CALL  WMD(MD,T0,FMD)
      DRND=RNDMD(-1)
      IF(DRND.GT.FMD)                   GO  TO  1
      PME(4,LD)  =-PME(4,L2)
      PME(5,LD)  =-PME(5,L2)
      PME(6,LD)  =-PME(6,L2)
      PME(9,LD)  =MD
      IME(1,LD)  = IME(1,LD)  +IME(1,MV+3)
      IME(2,LD)  =0
      IME(3,LD)  =0
      IME(4,LD)  =1
      IME(5,LD)  =IDINTG(1000.*T0)
      IF(IME(5,LD).LT.1)  IME(5,LD)=1
      NP=NP-1
      IND=MV+2*MQ-1
      DO  11 K=1,9
      PME(K,MV+3)=PME(K,MV+2)
      IF(K.LE.5) IME(K,MV+3)=IME(K,MV+2)
   11 CONTINUE
      RETURN
      END  
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE  WMD(MD,TAU0,FMD)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 MD,TAU0,FMD,M0
      DATA   M0/1.232D0/
      EN=(MD**2+0.940**2-0.140**2)/(2.*MD)
      Q=DSQRT(EN**2-0.940**2)
      R=Q/0.140
      G=0.47*Q*R*R/(1.+0.6*R*R)
      DRND=RNDMD(-1)
      IF(DRND.LT.1.D-10)  DRND=1.D-10
      TAU0=-1./(5.06*G)*DLOG(DRND)
      IF(TAU0.LE.0.001)  TAU0=0.0011
      A=G*G/4.
      FMD=A/((MD-M0)**2+A)
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE DECREN(LR,NUR,MV,NP)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*8 PNAR,PNAJ
      COMMON/PARTCL/PPTCL(9,499),NPTCL,IORIG(499),IDENT(499)
     *,IDCAY(499)
      COMMON /MEMAGT/PME(9,5999),IME(5,5999)
      COMMON/NCASCA/NCAS,NCPRI
      COMMON/IDPME/ IDPME(5999)
      COMMON/TLIMIT/TLIMIT
      COMMON/ACTIM/TINT
      COMMON /INDRES/ INDRES
      COMMON/PORIG/IORI(3,5999)
      COMMON/IDpart/IDp1,IDp2
      DIMENSION PIN(9),IIN(5)
C
      NDECAY=1
      IF(TINT.GT.TLIMIT)  NDECAY=5
C
      NP=0
      DO   K=1,9
      PIN(K)=PME(K,LR)
      IF(K.LE.5)  IIN(K)=IME(K,LR)
      ENDDO
      IDR=IDPME(LR)
      IDp1=IORI(1,LR)
      IDp2=IORI(2,LR)
      CALL  PANUID(IDR,NUR,PNAR)
c-Sob      IF(NCAS.GE.NCPRI) THEN
c-Sob        write(16,18) PNAR,(PME(I,LR),I=4,9),(IME(K,LR),K=1,5)
c-Sob        write( *,18) PNAR,(PME(I,LR),I=4,9),(IME(K,LR),K=1,5)
c-Sob   18 FORMAT(1X,'DECAY:',A5,6(1X,F7.3),4I2,I10)
c-Sob      ENDIF
      NPTCL=1
      PPTCL(1,1)=PIN(4)
      PPTCL(2,1)=PIN(5)
      PPTCL(3,1)=PIN(6)
      PPTCL(4,1)=PIN(8)+PIN(9)
      PPTCL(5,1)=PIN(9)
      IDENT(1)=IDR
      IORIG(1)=0
      IDCAY(1)=0
C
      do  ND=1,NDECAY
      IF(NPTCL.LE.0)   GO  TO  116
      NDEC=0
      NPTCL1=NPTCL
      DO 114  I=1,NPTCL1
      IDI=IDENT(I)
C       IF(IDI.EQ.110.OR.IDI.EQ.220.OR.IDI.EQ.230.OR.IDI.EQ.-230)
      IF(IDI.EQ.110.OR.IDI.EQ.230.OR.IDI.EQ.-230)
     *                   GO  TO  114
      IF(IDI.EQ.1230.OR.IDI.EQ.-1230)
     *                   GO  TO  114
c       IF(INDRES.EQ.1.AND.(IDI.EQ.220.or.IDI.EQ.221.or.IDI.EQ.331))
c     *  goto 114
      IF(INDRES.EQ.1.AND.(IDI.EQ.220.or.IDI.EQ.221.or.IDI.EQ.331))
     &  then
        if(IDR.ne.IDI)     go to 114
        if(Tint.gt.Tlimit) go to 114
      endif
      CALL DECAY(I,IPOINT)
      IF(IPOINT.LT.0)  GO  TO  114
      NDEC=NDEC+1
      DO  J=1,9
        PPTCL(J,I)=PPTCL(J,NPTCL)
      ENDDO
      IDENT(I)=IDENT(NPTCL)
      IORIG(I)=IORIG(NPTCL)
      IDCAY(I)=IDCAY(NPTCL)
      NPTCL=NPTCL-1
  114   CONTINUE
      IF(NDEC.EQ.0)    GO  TO  116
      enddo
  116 CONTINUE
      IF(NPTCL.LE.1)  THEN
c-Sob        write(16,21) NPTCL,IDI,IPOINT
c-Sob        write( *,21) NPTCL,IDI,IPOINT
c-Sob   21 FORMAT(5X,'DECREN : NPTCL=',I2,2X,'IDI=',I6,2X,'IPOINT=',I3)
      RETURN
      ENDIF
      IF((MV+NPTCL).GT.5999)  THEN
c-Sob        write(16,'(5X,''DECREN : MV+NPTCL>5999'')')
      RETURN
      ENDIF
      DO  J=1,NPTCL
      M=MV+J
      IDJ=IDENT(J)
      CALL  PANUID(IDJ,JP,PNAJ)
      IME(1,M)=IDINT(1.001*CHARGE(IDJ))
      IME(2,M)=0
C  !!!
      IF(JP.GE.57.AND.JP.LE.65) IME(2,M)=IDJ
      IME(3,M)=IS(IDJ)
      IME(4,M)=IB(IDJ)
      IME(5,M)=0
      PME(4,M)=PPTCL(1,J)
      PME(5,M)=PPTCL(2,J)
      PME(6,M)=PPTCL(3,J)
      PME(8,M)=PPTCL(4,J)-PPTCL(5,J)
      PME(9,M)=PPTCL(5,J)
      IDPME(M)=IDJ
c-Sob        IF(NCAS.GE.NCPRI)
c-Sob     *  write(16,'(1X,''====>:'',A5,6(1X,F11.4),4I3,I10)')
c-Sob     &  PNAJ,(PME(I,M),I=4,9),(IME(K,M),K=1,5)
      ENDDO
      NP=NPTCL
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE DECDVM(L,MV)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8  MN,MPI,MD
      COMMON /MEMAGT/ PME(9,5999),IME(5,5999)
      DIMENSION VD(3),PSN(3),PLN(3),IEN(8),IEP(8),CBR(4)
      DATA CBR/1.0,0.6666667,0.6666667,1.0/,
     *IEN/0,0,0,1,1,0,1,1/,
     *IEP/-1,-1,0,-1,0,1,1,1/,MN/0.940/,MPI/0.140/
      MD=PME(9,L)
      ESN=(MD**2+MN**2-MPI**2)/2./MD
      PS=DSQRT(ESN**2-MN**2)
      CTS=1.-2.*RNDMD(-1)
      FIS=6.283185*RNDMD(-1)
      STS=DSQRT(1.-CTS**2)
      PSN(1)=PS*STS*DCOS(FIS)
      PSN(2)=PS*STS*DSIN(FIS)
      PSN(3)=PS*CTS
      ED=PME(8,L)+MD
      VD(1)=PME(4,L)/ED
      VD(2)=PME(5,L)/ED
      VD(3)=PME(6,L)/ED
      CALL  KINEMA(PSN,VD,PLN,CTL,STL,CFL,SFL,TL,MN)
      LN=MV+1
      PME(4,LN)=PLN(1)
      PME(5,LN)=PLN(2)
      PME(6,LN)=PLN(3)
      PME(8,LN)=TL
      PME(9,LN)=MN
      LP=MV+2
      PME(4,LP)=PME(4,L)-PLN(1)
      PME(5,LP)=PME(5,L)-PLN(2)
      PME(6,LP)=PME(6,L)-PLN(3)
      PME(8,LP)=PME(8,L)+MD-TL-MN-MPI
      PME(9,LP)=MPI
      ICD=IME(1,L)+2
      IBR=2*(ICD-1)+1
      DRND=RNDMD(-1)
      IF(DRND.GT.CBR(ICD))  IBR=IBR+1
      IME(1,LN)=IEN(IBR)
      IME(2,LN)=0
      IME(3,LN)=0
      IME(4,LN)=1
      IME(5,LN)=0
      IME(1,LP)=IEP(IBR)
      IME(2,LP)=0
      IME(3,LP)=0
      IME(4,LP)=0
      IME(5,LP)=0
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE RTAP10(NCAS,INTEL,NZON)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*4 UP1,UP2,S1,S2,S3,S4,S5
      COMMON/HCASC/ANUCL1,ANUCL2,ZNUCL1,ZNUCL2,T0,EPS1,EPS2,VPI,A1,A2,
     *C1,C2,D1,D2,R0N1,R0N2,TF01,TF02,RM1,RM2
      COMMON/UPAC/UP1(66000),IPER
      COMMON/UPAC2/UP2(70000),LU2
      COMMON/INDSTR/JPSTR,NREP,JG
      COMMON /IW10/ IW10
      NUNIT=2000
      NZON=0
      LU2=1
      REWIND 10
      IF(JG.eq.0)   JG=500
      IF(JG.GT.0) then
      DO  KG=1,JG
c        write( *,*) 'KG=', KG
      READ(10,END=98) S1,S2,S3,S4,S5,IS6,IS7,IS8,IS9,IS10
C  !!!
      if(IW10.EQ.1)  then
C  !!!
        NUN=IS10/NUNIT
        IF(NUN.ge.1) THEN
          do  LUN=1,NUN
            K1=(LUN-1)*NUNIT+1
            K2=LUN*NUNIT
            READ(10) (UP2(KL),KL=K1,K2)
          enddo
        ENDIF
        K3=NUN*NUNIT+1
c          write( *,*) 'K3,IS10,NUN=', K3,IS10,NUN
        IF(K3.LE.IS10) READ(10) (UP2(KL),KL=K3,IS10)
      endif
      NZON=NZON+1
      ENDDO
      GO TO 99
   98 CONTINUE
c        write( *,*) 'BACKSPACE' (not needed for UNIX )
c      BACKSPACE 10
C  ! THE NEXT LINE IS FOR EC-COMP. and GANIL IBM: END OF RECORD !
c      BACKSPACE 10
 99   A10=S1
      A20=S2
      Z10=S3
      Z20=S4
      T00=S5
      NCAS=IS6
      INTEL=IS7
      NZON=IS8
      IS100=IS10
c-Sob      write(16,'(1X,''FROM10: '',4(F5.0,1X),F7.3,1X,4I8)')
c-Sob     & A10,A20,Z10,Z20,T00,NCAS,INTEL,NZON,IS100
c-Sob      write( *,'(1X,''FROM10: '',4(F5.0,1X),F7.3,1X,4I8)')
c-Sob     & A10,A20,Z10,Z20,T00,NCAS,INTEL,NZON,IS100
      JG=NZON
      IF(DABS(ANUCL1-A10).le.0.1.and.
     &  DABS(ANUCL2-A20).le.0.1.and.
     &  DABS(ZNUCL1-Z10).le.0.1.and.
     &  DABS(ZNUCL2-Z20).le.0.1.and.
     &  DABS(T0-T00).le. 0.0001)     then
        IBRAN=IS9
c          CALL  RDMIN(IBRAN)    ! 13.02.96
      else
        NCAS=0
        INTEL=0
        JG=-1
      endif
      ELSE
      NCAS=0
      INTEL=0
      JG=-1
      ENDIF
      NZONS=NZON+1
c-Sob      write(16,'(/10X,''START NZON='',I5,'' NCAS='',I7,'' INTEL='',I7)')
c-Sob     & NZONS,NCAS,INTEL
c-Sob      write( *,'(/10X,''START NZON='',I5,'' NCAS='',I7,'' INTEL='',I7)')
c-Sob     & NZONS,NCAS,INTEL
      RETURN
C---------------------------------------------
      ENTRY WTAP10(NCAS,INTEL,NZON)
C
      LU2=LU2-1
      NZON=NZON+1
c      CALL  RDMOUT(IBRAN)      ! 13.02.96
      S1=ANUCL1
      S2=ANUCL2
      S3=ZNUCL1
      S4=ZNUCL2
      S5=T0
      WRITE(10) S1,S2,S3,S4,S5,NCAS,INTEL,NZON,IBRAN,LU2
C  !!!
      IF(IW10.EQ.1) THEN
C  !!!
      NUN=LU2/NUNIT
      IF(NUN.ge.1)  then
        DO  LUN=1,NUN
          K1=(LUN-1)*NUNIT+1
          K2=LUN*NUNIT
          WRITE(10) (UP2(KL),KL=K1,K2)
        ENDDO
      endif
      K3=NUN*NUNIT+1
c        write( *,*) 'WRITE: K3,LU2,NUN=', K3,LU2,NUN
      IF(K3.LE.LU2) WRITE(10) (UP2(KL),KL=K3,LU2)
      ENDIF
C  !!!!  THE NEXT TWO LINES ARE FOR EC-COMPUTER ONLY
c      ENDFILE 10        !!for FC
c      BACKSPACE 10      !!for EC
c-Sob      write(16,'(1X,4(F5.0,1X),F8.3,1X,3(I7,1X),I7)')
c-Sob     & ANUCL1,ANUCL2,ZNUCL1,ZNUCL2,T0,NCAS,INTEL,NZON,LU2
c-Sob      write( *,'(1X,4(F5.0,1X),F8.3,1X,3(I7,1X),I7)')
c-Sob     & ANUCL1,ANUCL2,ZNUCL1,ZNUCL2,T0,NCAS,INTEL,NZON,LU2
      LU2=1
      RETURN
      END

