
************************************************************************
*                                                                      *
*  The following subroutines are written by K.K.Gudima and V.D.Toneev  *
*  Joint Institute for Nuclear Research, Dubna,                        *
*  E-mail: gudima@acad.moldnet.md                                      *
*          gudima@cc.acad.md                                           *
*  E-mail: toneev@thsun1.jinr.ru                                       *
*                                                                      *
************************************************************************

C   * * * * * * * * *  11-16-94 03:27PM  MV <==> 5999 * * * * * *
      SUBROUTINE HEINEN(PIN,IIN,PN,IPN,MV,NP,NIN) 
      IMPLICIT REAL*8 (A-H,O-Z) 
      REAL *8 PIN,PN 
      CHARACTER*8 DTYP(11),TDIAG 
      COMMON/MEMAGT/ PME(9,5999),IME(5,5999) 
      DIMENSION  PIN(9),IIN(5),PN(9),IPN(5) 
      COMMON /IACT/ IACT 
      COMMON /NCASCA/ NCAS,NCPRI 
      COMMON /LOWMIS/ LOWMIS 
      COMMON /CSLID/ CLIDER(5999) 
      COMMON/COMENB/ ENBOU 
      COMMON/IORHEI/ IORHEI 
      COMMON/ITHEA/ITHEA(11) 
      COMMON /HELE/ IHELE 
      COMMON /IDPME/ IDPME(5999) 
      COMMON /ISLERR/ ISLERR   ! HADES1
      DIMENSION PSU(3) 
      DATA DTYP/'DIFTRI','PLANAR','UNCYLI','ELAST  ','ANNIH','DIFSMA', 
     *          'CHAINS','BINAR ','??????','REGTRI ','DOUBDI'/ 
c-Sob      IF(NCAS.GE.NCPRI) 
c-Sob     *write(16,601) PIN,IIN,PN,IPN 
c-Sob  601 FORMAT(1X,'HEINEN:',9(1X,E10.3),4I2,I10/ 
c-Sob     *       1X,'      +',9(1X,E10.3),4I2,I10) 
      ES0=DSQRT((PIN(8)+PIN(9)+PN(8)+PN(9))**2-(PIN(4)+PN(4))**2- 
     -(PIN(5)+PN(5))**2-(PIN(6)+PN(6))**2) 
      CS0=FLOAT(IIN(1)+IPN(1)) 
      SS0=FLOAT(IIN(3)+IPN(3)) 
      BS0=FLOAT(IIN(4)+IPN(4)) 
      ESU=0. 
      CSU=0. 
      SSU=0. 
      BSU=0. 
      LOWMIS=0 
      DO  10  K=1,3 
   10 PSU(K)=0. 
      CALL ACTNAM(PIN,PN,MV,NP,NIN) 
        if(ISLERR.eq.1)return   ! HADES1
c        if(ISLERR.eq.1)then
c          CONTINUE!WRITE(25,*)'    ACTNAM'
c          return
c        end if
      CONTINUE 
      IF(NIN.EQ.1)  RETURN 
      DO  14   K=1,NP 
      M=MV+K 
      CSU=CSU+FLOAT(IME(1,M)) 
      SSU=SSU+FLOAT(IME(3,M)) 
      BSU=BSU+FLOAT(IME(4,M)) 
      PSU(1)=PSU(1)+PME(4,M) 
      PSU(2)=PSU(2)+PME(5,M) 
      PSU(3)=PSU(3)+PME(6,M) 
      ETM=DSQRT(PME(4,M)**2+PME(5,M)**2+PME(6,M)**2+PME(9,M)**2) 
      ESU=ESU+ETM 
   14 CONTINUE 
      ESU=DSQRT(ESU**2-PSU(1)**2-PSU(2)**2-PSU(3)**2) 
      IF(DABS(ESU-ES0).GT.0.010.OR.DABS(CSU-CS0).GT.0.1. 
     *OR.DABS(SSU-SS0).GT.0.1.OR.DABS(BSU-BS0).GT.0.1) 
     *LOWMIS=1 
      DO  19  I=1,11 
   19 IF(ITHEA(I).NE.0) TDIAG=DTYP(I) 
c-Sob      IF(LOWMIS.EQ.1) 
c-Sob     *write(16,1000) ES0,CS0,SS0,BS0,ESU,CSU,SSU,BSU 
c-Sob 1000 FORMAT(1X,'ES0,CS0,SS0,BS0,ESU,CSU,SSU,BSU=', 
c-Sob     *F9.3,3(2X,F3.0)/33X,F9.3,3(2X,F3.0)) 
      IF(DABS(PSU(1)).GT.0.001.OR.DABS(PSU(2)).GT.0.001.OR. 
     *DABS(PSU(3)).GT.0.001) 
     *LOWMIS=1 
c-Sob      IF(LOWMIS.EQ.1) 
c-Sob     *write(16,1001) PSU,TDIAG 
c-Sob 1001 FORMAT(1X,'PSU=',3(F9.3,2X),5X,A8) 
c-SobC     write(16,*) 'NP,IORHEI=', NP,IORHEI 
      IF(IORHEI.NE.1)                                  GO  TO  112 
      IF(NP.EQ.1.OR.(IABS(IIN(4))+IABS(IPN(4))).LE.0)  GO  TO  112 
C 
      NBA=0 
      MB=MV+1 
      DO  102  K=1,NP 
      M=MV+K 
      IF(IME(4,M).EQ.0)    GO  TO  102 
      NBA=NBA+1 
      DO  100  L=1,9 
      TEMP=PME(L,MB) 
      PME(L,MB)=PME(L,M) 
  100 PME(L,M)=TEMP 
      DO  101  L=1,5 
      ITEMP=IME(L,MB) 
      IME(L,MB)=IME(L,M) 
  101 IME(L,M)=ITEMP 
      TEMC=CLIDER(MB) 
      CLIDER(MB)=CLIDER(M) 
      CLIDER(M)=TEMC 
      IDTE=IDPME(MB) 
      IDPME(MB)=IDPME(M) 
      IDPME(M)=IDTE 
C      write(16,*) 'M,NBA,MB,IDPME(MB) =', M,NBA,MB,IDPME(MB) 
      MB=MV+3 
  102 CONTINUE 
      IF(NP.NE.2)          GO  TO  105 
      IF(NBA.EQ.2)         GO  TO  107 
      DO 103 L=1,9 
  103 PME(L,MV+3)=PME(L,MV+2) 
      DO 104 L=1,5 
  104 IME(L,MV+3)=IME(L,MV+2) 
      CLIDER(MV+3)=CLIDER(MV+2) 
      IDPME(MV+3)=IDPME(MV+2) 
                        GO  TO  112 
  105 IF(NBA.EQ.2)      GO  TO  107 
      MLI=MV+1 
      YL=YMEM(MLI) 
      DO  106  N=2,NP 
      YN=YMEM(MV+N) 
      IF(YL.GT.YN)      GO  TO  106 
      YL=YN 
      MLI=MV+N 
  106 CONTINUE 
                        GO  TO  109 
  107 Y3=YMEM(MV+3) 
      Y1=YMEM(MV+1) 
      IF(Y3.GT.Y1)      GO  TO  112 
      MLI=MV+1 
  109 DO  110  L=1,9 
      TEMP=PME(L,MV+3) 
      PME(L,MV+3)=PME(L,MLI) 
  110 PME(L,MLI)=TEMP 
      DO  111  L=1,5 
      ITEMP=IME(L,MV+3) 
      IME(L,MV+3)=IME(L,MLI) 
  111 IME(L,MLI)=ITEMP 
      TEMC=CLIDER(MV+3) 
      CLIDER(MV+3)=CLIDER(MLI) 
      CLIDER(MLI)=TEMC 
      IDTE=IDPME(MV+3) 
      IDPME(MV+3)=IDPME(MLI) 
      IDPME(MLI)=IDTE 
  112 CONTINUE 
      IF(NCAS.LT.NCPRI)  RETURN 
c-Sob      write(16,599) NCAS,TDIAG,NP 
c-Sob  599 FORMAT(1X,'NCAS=',I6,' RESULTS FROM ',A8,' NP=',I5) 
      DO  600  L=1,NP 
      M=MV+L 
      IF(NP.EQ.2.AND.L.EQ.2.AND.(IABS(IIN(4))+IABS(IPN(4))).GT.0) 
     *M=M+1 
c-Sob      write(16,602) M,(PME(I,M),I=4,6),PME(9,M),(IME(J,M),J=1,5), 
c-Sob     * CLIDER(M),IDPME(M) 
c-Sob  602 FORMAT(1X,I3,4(1X,E11.4),4I2,I10,1X,F5.3,I5) 
  600 CONTINUE 
      RETURN 
      END 
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
      FUNCTION YMEM(M) 
      IMPLICIT REAL*8 (A-H,O-Z) 
      REAL *8 YMEM 
      COMMON /MEMAGT/ PME(9,5999),IME(5,5999) 
      E=DSQRT(PME(4,M)**2+PME(5,M)**2+PME(6,M)**2+PME(9,M)**2) 
      YMEM=0.5*DLOG((E+PME(6,M))/(E-PME(6,M))) 
      RETURN 
      END 
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
      SUBROUTINE  LOWPEC(PIN,IIN,PN,IPN) 
      IMPLICIT REAL*8 (A-H,O-Z) 
      REAL *8 PIN,PN 
      CHARACTER*8 DTYP(11),TDIAG 
      COMMON /MEMAGT/ PME(9,5999),IME(5,5999) 
      COMMON/ITHEA/ITHEA(11) 
      COMMON /IDPME/ IDPME(5999) 
      DATA DTYP/'DIFTRI','PLANAR','UNCYLI','ELAST  ','ANNIH','DIFSMA', 
     *          'CHAINS','BINAR ','??????','REGTRI ','DOUBDI'/ 
      DIMENSION PSU(3),PIN(9),PN(9),IIN(5),IPN(5) 
      T0=PIN(8) 
      AM1=PIN(9) 
      AM2=PN(9) 
      P0=DSQRT(T0*(T0+2.*AM1)) 
      S0=DSQRT((T0+AM1+AM2)**2-P0**2) 
      IC0=IIN(1)+IPN(1) 
      IS0=IIN(3)+IPN(3) 
      IB0=IIN(4)+IPN(4) 
      NA=10 
      MV=0 
      DO  100 N=1,NA 
      PSU(1)=0. 
      PSU(2)=0. 
      PSU(3)=0. 
      ESU=0. 
      ICSU=0 
      ISSU=0 
      IBSU=0 
      CALL ACTNAM(PIN,PN,MV,NP,NIN) 
c-Sob      write(16,999) 
c-Sob      write( *,999) 
c-Sob  999 FORMAT( 4X,'PX',10X,'PY',10X,'PZ',11X,'E', 
c-Sob     *       9X,'Q',2X,'S',2X,'B',5X,'M',5X,'ID'/) 
      DO  12  L=1,NP 
      M=MV+L 
      PSU(1)=PSU(1)+PME(4,M) 
      PSU(2)=PSU(2)+PME(5,M) 
      PSU(3)=PSU(3)+PME(6,M) 
      EL=DSQRT(PME(4,M)**2+PME(5,M)**2+PME(6,M)**2+PME(9,M)**2) 
      ESU=ESU+EL 
      ICSU=ICSU+IME(1,M) 
      ISSU=ISSU+IME(3,M) 
      IBSU=IBSU+IME(4,M) 
c-Sob      write(16,1000) PME(4,M),PME(5,M),PME(6,M),EL,IME(1,M),IME(3,M), 
c-Sob     *IME(4,M),PME(9,M),IDPME(M) 
c-Sob      write( *,1000) PME(4,M),PME(5,M),PME(6,M),EL,IME(1,M),IME(3,M), 
c-Sob     *IME(4,M),PME(9,M),IDPME(M) 
   12 CONTINUE 
      DO  19  I=1,11 
   19 IF(ITHEA(I).NE.0) TDIAG=DTYP(I) 
c-Sob      write(16,1001) PSU(1),PSU(2),PSU(3),ESU,ICSU,ISSU,IBSU 
c-Sob      write( *,1001) PSU(1),PSU(2),PSU(3),ESU,ICSU,ISSU,IBSU 
c-Sob      write(16,1002) P0,S0,IC0,IS0,IB0,TDIAG 
c-Sob      write( *,1002) P0,S0,IC0,IS0,IB0,TDIAG 
 1000 FORMAT( 1X,4(E11.4,1X),3I3,1X,E11.4,I6) 
 1001 FORMAT(1X,43('-')/1X,4(E11.4,1X),3I5) 
 1002 FORMAT(1X,43('-')/1X,24X,2(E11.4,1X),3I5,1X,A8/1X,78('*')) 
  100 CONTINUE 
      RETURN 
      END 
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
      BLOCK DATA  BPNAME 
C                ----- 
      CHARACTER*8  PNAME 
      COMMON /PNAME/PNAME(65),IDENTP(65) 
      DATA PNAME/ 
     *  'PI+  ','PI-  ','K+   ','K-   ','K0   ','AK0  ','PI0  ', 
     *  'ETA  ','ETAP ','RHO+ ','RHO- ','K*+  ','K*-  ','K*0  ', 
     *  'AK*0 ','RHO0 ','OMEG ','PHI  ', 
     *  'AP   ','AN   ','AS+  ','AS-  ','AS0  ','AXI- ','AXI0 ', 
     *  'AL   ','ADL++','ADL+ ','ADL- ','ADL0 ','AS*+ ','AS*- ', 
     *  'AS*0 ','AXI*-','AXI*0','AOM- ', 
     *  'P    ','N    ','S+   ','S-   ','S0   ','XI-  ','XI0  ', 
     *  'L    ','DL++ ','DL+  ','DL-  ','DL0  ','S*+  ','S*-  ', 
     *  'S*0  ','XI*- ','XI*0 ','OM-  ','KS   ','KL   ', 
     *  'GM   ','E-   ','E+   ','MU+  ','MU-  ','NUE  ','NUM  ', 
     *  'ANUE ','ANUM '/ 
C 
      DATA IDENTP/ 
     *     120 ,  -120 ,   130 ,  -130 ,   230 ,  -230 ,   110 , 
     *     220 ,   330 ,   121 ,  -121 ,   131 ,  -131 ,   231 , 
     *    -231 ,   111 ,   221 ,   331 , 
     *   -1120 , -1220 , -1130 , -2230 , -1230 , -2330 , -1330 , 
     *   -2130 , -1111 , -1121 , -2221 , -1221 , -1131 , -2231 , 
     *   -1231 , -2331 , -1331 , -3331 , 
     *    1120 ,  1220 ,  1130 ,  2230 ,  1230 ,  2330 ,  1330 , 
     *    2130 ,  1111 ,  1121 ,  2221 ,  1221 ,  1131 ,  2231 , 
     *    1231 ,  2331 ,  1331 ,  3331 ,20 ,   -20 , 
     *      10 ,    12 ,   -12 ,   -14 ,    14 ,    11 ,    13 , 
     *     -11 ,   -13 / 
C     ----------------------------------------------------------------- 
      END 
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
      SUBROUTINE INITAM 
      REAL*8 ENBOU 
      COMMON/NODCAY/NODCAY,NOETA,NOPI0,NONUNU,NOEVOL,NOHADR,NOKA0 
      COMMON/ITAPES/ITDKY,ITEVT,ITCOM,ITLIS 
      COMMON/PRINTS/IPRINT 
      COMMON/NOTRE/ NOTRE 
      COMMON/COMENB/ ENBOU 
      COMMON/COMFR/  ICMS 
      COMMON/YESELA/ YESELA 
      COMMON/KEYHH/ KEYHH 
      COMMON/KEYPLA/ KEYPLA 
      COMMON/COMQSE/ QSEE,QVSEE 
      COMMON/CINSID/ INSIDE 
      COMMON/CVALON/ IVALON 
      COMMON/COMWTI/ WTIME 
      COMMON/IORHEI/ IORHEI 
      COMMON/ISOB3/ ISOB3 
      LOGICAL WTIME 
      LOGICAL NODCAY,NOETA,NOPI0,NONUNU,NOEVOL,NOHADR,NOKA0 
      LOGICAL IPRINT 
      LOGICAL LPRNT 
      LOGICAL NOTRE 
      LOGICAL YESELA 
      LOGICAL KEYHH 
      LOGICAL KEYPLA 
      LOGICAL QSEE,QVSEE 
      ITCOM=5 
      ITLIS=6
      ITEVT=17 
      ITDKY=18 
C 
      NODCAY=.FALSE. 
      NOTRE=.FALSE. 
      LPRNT=.FALSE. 
c      LPRNT=.TRUE. 
      IPRINT=.FALSE. 
      YESELA=.FALSE. 
      KEYHH=.FALSE. 
      KEYPLA=.FALSE. 
      QSEE=.FALSE. 
      QVSEE=.FALSE. 
C     WTIME=.TRUE. 
      WTIME=.FALSE. 
      ENBOU= 4.4 
      INSIDE=0 
      IVALON=0 
      IORHEI=1 
      ISOB3=1 
C 
C       IF LAB. FRAME 
      ICMS=0 
C       IF  C.M. FRAME 
      ICMS=1 
      CALL SETCON 
      CALL SETDKY(LPRNT) 
C     WRITE(7,9999) 
 9999 FORMAT(2X,'END INITAM') 
      RETURN 
      END 
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
      SUBROUTINE  ACTNAM(PIN,PN,MV,NP,NIN) 
      IMPLICIT REAL*8 (A-H,O-Z) 
      REAL *8 PIN,PN,MA 
      CHARACTER*8 PNA1,PNA2,PNAJ 
      LOGICAL GH1H2,YESELA 
      COMMON /INDDEC/INDDEC 
      COMMON /NCASCA/ NCAS,NCPRI 
      COMMON/PARTCL/PPTCL(9,499),NPTCL,IORIG(499),IDENT(499) 
     *,IDCAY(499) 
      COMMON /ITHEA/ ITHEA(11) 
      COMMON /H1H2/ GH1H2(11) 
      COMMON /YESELA/ YESELA 
      COMMON /COMLID/ PLIDER(499) 
      COMMON /IDN12/ ID1,ID2 
      COMMON /IDN120/ ID10,ID20 
      COMMON /ISLERR/ ISLERR   ! HADES1
      DIMENSION  PIN(9),PN(9),RS(4),PS(3) 
C 
   1  CONTINUE 
      ID1=ID10 
      ID2=ID20 
      NPTCL=0 
      DO  10  I=1,11 
      ITHEA(I)=0 
   10 GH1H2(I)=.TRUE. 
      GH1H2(4)=YESELA 
C 
      PX1=PIN(4) 
      PY1=PIN(5) 
      PZ1=PIN(6) 
      PX2= PN(4) 
      PY2= PN(5) 
      PZ2= PN(6) 
      AM1=PIN(9) 
      AM2= PN(9) 
      ECM=DSQRT((PIN(8)+PIN(9)+PN(8)+PN(9))**2-(PIN(4)+PN(4))**2- 
     -(PIN(5)+PN(5))**2-(PIN(6)+PN(6))**2) 
      CALL  PANUID(ID1,IK1,PNA1) 
      CALL  PANUID(ID2,IK2,PNA2) 
      CALL  CROSEC(0,ID1,ID2,PX1,PY1,PZ1,AM1,PX2,PY2,PZ2,AM2,SIEL,0) 
      CALL  CROSEC(1,ID1,ID2,PX1,PY1,PZ1,AM1,PX2,PY2,PZ2,AM2,SITO,0) 
      CALL  CROSEC(2,ID1,ID2,PX1,PY1,PZ1,AM1,PX2,PY2,PZ2,AM2,SIEX,0) 
      CALL  CROSEC(3,ID1,ID2,PX1,PY1,PZ1,AM1,PX2,PY2,PZ2,AM2,SIAN,0) 
c-Sob      IF(SITO.LE.0..OR.SIEL.LE.0.) 
c-Sob     *write(16,601) PNA1,PX1,PY1,PZ1,AM1,PNA2,PX2,PY2,PZ2,AM2, 
c-Sob     *             SITO,SIEL,SIEX,SIAN 
c-Sob 601  FORMAT(1X,'ACTNAM1: (',A8,')',4(F8.4,1X)/ 
c-Sob     *8X,'+(',A8,')',4(F8.4,1X)/ 
c-Sob     *8X,'SITO,SIEL,SIEX,SIAN=',4(E11.4,1X)) 
      IF(SITO.LE.0.)  SITO=1. 
      IF(SIEL.LE.0.)  SIEL=SITO 
C 
c-Sob      IF(NCAS.GE.NCPRI) 
c-Sob     *write(16,600) PNA1,PX1,PY1,PZ1,AM1,PNA2,PX2,PY2,PZ2,AM2, 
c-Sob     *             SITO,SIEL,SIEX,SIAN 
c-Sob      IF(NCAS.GE.NCPRI) 
c-Sob     *write( *,600) PNA1,PX1,PY1,PZ1,AM1,PNA2,PX2,PY2,PZ2,AM2, 
c-Sob     *             SITO,SIEL,SIEX,SIAN 
c-Sob 600  FORMAT(1X,'ACTNAM: (',A8,')',4(F8.4,1X)/ 
c-Sob     *8X,'+(',A8,')',4(F8.4,1X)/ 
c-Sob     *8X,'SITO,SIEL,SIEX,SIAN=',4(E11.4,1X)) 
      CALL COLLHH(ID1,AM1,PX1,PY1,PZ1,ID2,AM2,PX2,PY2,PZ2, 
     *SITO,SIAN,SIEL) 
        if(ISLERR.eq.1)return   ! HADES1
c        if(ISLERR.eq.1)then
c          CONTINUE!WRITE(25,*)'   COLLHH'
c          return
c        end if
      NP=0 
      NIN=0 
C 
      IF(INDDEC.EQ.0)  GO  TO 116 
C 
      IF(NPTCL.LE.0)   RETURN 
      DO 115  ND=1,5 
      IF(NPTCL.LE.0)   GO  TO  1 
      NDEC=0 
      NPTCL1=NPTCL 
      DO 114  I=1,NPTCL1 
      IDI=IDENT(I) 
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
      PLIDER(I)=PLIDER(NPTCL) 
      NPTCL=NPTCL-1 
  114 CONTINUE 
      IF(NDEC.EQ.0)     GO  TO  116 
  115 CONTINUE 
  116 IF(NPTCL.LE.0)  GO  TO 1 
      DO  13  J=1,NPTCL 
      MA=PPTCL(5,J) 
      PS(1)=PPTCL(1,J) 
      PS(2)=PPTCL(2,J) 
      PS(3)=PPTCL(3,J) 
      RS(1)=PPTCL(6,J) 
      RS(2)=PPTCL(7,J) 
      RS(3)=PPTCL(8,J) 
      RS(4)=PPTCL(9,J) 
      IDJ=IDENT(J) 
      CALL  PANUID(IDJ,JP,PNAJ) 
      CALL  RASPAN(IDJ,JP,MA,PS,RS,NP,MV) 
   13 CONTINUE 
      DO  14  I=1,11 
      IF(GH1H2(I))   ITHEA(I)=1 
   14 CONTINUE 
      RETURN 
      END 
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
      SUBROUTINE  RASPAN(IDJ,JP,MA,P,R,NP,MV) 
      IMPLICIT REAL*8 (A-H,O-Z) 
      REAL *8 P,MA,R 
      COMMON/IACT/ IACT 
      COMMON /MEMAGT/ PME(9,5999),IME(5,5999) 
      COMMON /COMCHA/ ICHAM(18),ICHAB(18) 
      COMMON /COMSTR/ ISTR(36) 
      COMMON /COMLID/ PLIDER(499) 
      COMMON /CSLID/ CLIDER(5999) 
      COMMON /IDPME/ IDPME(5999) 
      DIMENSION  P(3),R(4) 
      IF(JP.LT.1.OR.JP.GT.65)         GO  TO  12 
      M=MV+NP+1 
      IF(M.GT.5999)                   GO  TO  14 
      CHARx=CHARGE(IDJ)*1.001d0 
      IME(1,M)=IDINT(CHARx) 
c      IME(1,M)=IDNINT(CHARx) 
c      write(16,*) 'from RASPAN: IDJ,CHARx,IME(1,M)=', IDJ,CHARx,IME(1,M) 
      IME(2,M)=0 
      IF(JP.GE.57.AND.JP.LE.65) IME(2,M)=IDJ 
      IME(3,M)=IS(IDJ) 
      IME(4,M)=IB(IDJ) 
      IF(JP.LE. 7)               IME(5,M)=0 
      IF(JP.GE. 8.AND.JP.LE.18)  IME(5,M)=IDINTG(1000.  *TAUN(JP)) 
      IF(JP.GE.19.AND.JP.LE.26)  IME(5,M)=0 
      IF(JP.GE.27.AND.JP.LE.35)  IME(5,M)=IDINTG(1000.  *TAUN(JP)) 
      IF(JP.GE.36.AND.JP.LE.44)  IME(5,M)=0 
      IF(JP.GE.45.AND.JP.LE.53)  IME(5,M)=IDINTG(1000.  *TAUN(JP)) 
      IF(JP.GE.54)               IME(5,M)=0 
   99 CONTINUE 
      IF(IACT.GT.2)  GO  TO  98 
      PME(1,M)=0. 
      PME(2,M)=0. 
      PME(3,M)=0. 
      PME(7,M)=0. 
      GO  TO  97 
   98 PME(1,M)=R(1) 
      PME(2,M)=R(2) 
      PME(3,M)=R(3) 
      PME(7,M)=R(4) 
   97 PME(4,M)=P(1) 
      PME(5,M)=P(2) 
      PME(6,M)=P(3) 
      PME(9,M)=MA 
      NP=NP+1 
      CLIDER(M)=PLIDER(NP) 
      IDPME(M)=IDJ 
      RETURN 
   12 CONTINUE
c-Sob      write(16,13) JP 
c-Sob   13 FORMAT(5X,'RASPAN: JP= ',I5) 
      RETURN 
   14 CONTINUE
c-Sob      write(16,15) 
c-Sob   15 FORMAT(5X,'RASPAN: MV+NP+1>5999') 
      RETURN 
      END 
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
      SUBROUTINE   PANUID(ID,N,PN) 
      CHARACTER*8 PN,PNAME 
      COMMON /PNAME/PNAME(65),IDENTP(65) 
      N=0 
      DO  10  I=1,65 
      IF(IDENTP(I).EQ.ID)  N=I 
      IF(IDENTP(I).EQ.ID)  PN=PNAME(I) 
   10 CONTINUE 
c-Sob      IF(N.EQ.0)   write(16,991) ID 
c-Sob      IF(N.EQ.0)   write( *,991) ID 
c-Sob  991 FORMAT(2X,'PANUID:  ID=',I5) 
      RETURN 
      ENTRY IDPANU(ID,N,PN) 
c-Sob      IF(N.LT.1.OR.N.GT.65)  write(16,992) N 
c-Sob  992 FORMAT(2X,'IDPANU:   N=',I5) 
      ID=IDENTP(N) 
      PN=PNAME(N) 
      RETURN 
      ENTRY LABEID(ID,N,PN) 
      N=0 
      DO  11  I=1,65 
      IF(PNAME(I).EQ.PN)  N=I 
      IF(PNAME(I).EQ.PN)  ID=IDENTP(I) 
   11 CONTINUE 
c-Sob      IF(N.EQ.0)   write(16,993) PN 
c-Sob  993 FORMAT(2X,'LABEID:  PN=',A8) 
      RETURN 
      END 
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
      SUBROUTINE PANUN(P,IP,N) 
      IMPLICIT REAL*8 (A-H,O-Z) 
      REAL *8 P,M 
      CHARACTER*8 PN 
      DIMENSION P(9),IP(5) 
      N=0 
      IC=IP(1) 
      IL=IP(2) 
      IS=IP(3) 
      IB=IP(4) 
      IR=IP(5) 
      M =P(9) 
      IF(IL.NE.0)  CALL  PANUID(IL,N,PN) 
      IF(IP(2).NE.0)  RETURN 
      IF(IC.EQ. 1.AND.IS.EQ. 0.AND.IB.EQ.0.AND.IR.EQ.0)     N=1 
      IF(IC.EQ.-1.AND.IS.EQ. 0.AND.IB.EQ.0.AND.IR.EQ.0)     N=2 
      IF(IC.EQ. 1.AND.IS.EQ. 1.AND.IB.EQ.0.AND.IR.EQ.0)     N=3 
      IF(IC.EQ.-1.AND.IS.EQ.-1.AND.IB.EQ.0.AND.IR.EQ.0)     N=4 
      IF(IC.EQ. 0.AND.IS.EQ. 1.AND.IB.EQ.0.AND.IR.EQ.0)     N=5 
      IF(IC.EQ. 0.AND.IS.EQ.-1.AND.IB.EQ.0.AND.IR.EQ.0)     N=6 
      IF(IC.EQ. 0.AND.IS.EQ. 0.AND.IB.EQ.0.AND.IR.EQ.0.AND. 
     *DABS(M-.140).LT.0.01)                                 N=7 
      IF(IC.EQ. 0.AND.IS.EQ. 0.AND.IB.EQ.0.AND.IR.NE.0.AND. 
     *DABS(M-.549).LT.0.01)                                 N=8 
      IF(IC.EQ. 0.AND.IS.EQ. 0.AND.IB.EQ.0.AND.IR.NE.0.AND. 
     *DABS(M-.958).LT.0.01)                                 N=9 
      IF(IC.EQ. 1.AND.IS.EQ. 0.AND.IB.EQ.0.AND.IR.NE.0)     N=10 
      IF(IC.EQ.-1.AND.IS.EQ. 0.AND.IB.EQ.0.AND.IR.NE.0)     N=11 
      IF(IC.EQ. 1.AND.IS.EQ. 1.AND.IB.EQ.0.AND.IR.NE.0)     N=12 
      IF(IC.EQ.-1.AND.IS.EQ.-1.AND.IB.EQ.0.AND.IR.NE.0)     N=13 
      IF(IC.EQ. 0.AND.IS.EQ. 1.AND.IB.EQ.0.AND.IR.NE.0)     N=14 
      IF(IC.EQ. 0.AND.IS.EQ.-1.AND.IB.EQ.0.AND.IR.NE.0)     N=15 
      IF(IC.EQ. 0.AND.IS.EQ. 0.AND.IB.EQ.0.AND.IR.NE.0.AND. 
     *DABS(M-.770).LT.0.01)                                 N=16 
      IF(IC.EQ. 0.AND.IS.EQ. 0.AND.IB.EQ.0.AND.IR.NE.0.AND. 
     *DABS(M-.783).LT.0.01)                                 N=17 
      IF(IC.EQ. 0.AND.IS.EQ. 0.AND.IB.EQ.0.AND.IR.NE.0.AND. 
     *DABS(M-1.020).LT.0.01)                                N=18 
      IF(IC.EQ.-1.AND.IS.EQ. 0.AND.IB.EQ.-1.AND.IR.EQ.0)    N=19 
      IF(IC.EQ. 0.AND.IS.EQ. 0.AND.IB.EQ.-1.AND.IR.EQ.0)    N=20 
      IF(IC.EQ.-1.AND.IS.EQ. 1.AND.IB.EQ.-1.AND.IR.EQ.0)    N=21 
      IF(IC.EQ. 1.AND.IS.EQ. 1.AND.IB.EQ.-1.AND.IR.EQ.0)    N=22 
      IF(IC.EQ. 0.AND.IS.EQ. 1.AND.IB.EQ.-1.AND.IR.EQ.0)    N=23 
      IF(IC.EQ. 1.AND.IS.EQ. 2.AND.IB.EQ.-1.AND.IR.EQ.0)    N=24 
      IF(IC.EQ. 0.AND.IS.EQ. 2.AND.IB.EQ.-1.AND.IR.EQ.0)    N=25 
      IF(IC.EQ. 0.AND.IS.EQ. 1.AND.IB.EQ.-1.AND.IR.EQ.0)    N=26 
      IF(IC.EQ.-2.AND.IS.EQ. 0.AND.IB.EQ.-1.AND.IR.NE.0)    N=27 
      IF(IC.EQ.-1.AND.IS.EQ. 0.AND.IB.EQ.-1.AND.IR.NE.0)    N=28 
      IF(IC.EQ. 1.AND.IS.EQ. 0.AND.IB.EQ.-1.AND.IR.NE.0)    N=29 
      IF(IC.EQ. 0.AND.IS.EQ. 0.AND.IB.EQ.-1.AND.IR.NE.0)    N=30 
      IF(IC.EQ.-1.AND.IS.EQ. 1.AND.IB.EQ.-1.AND.IR.NE.0)    N=31 
      IF(IC.EQ. 1.AND.IS.EQ. 1.AND.IB.EQ.-1.AND.IR.NE.0)    N=32 
      IF(IC.EQ. 0.AND.IS.EQ. 1.AND.IB.EQ.-1.AND.IR.NE.0)    N=33 
      IF(IC.EQ. 1.AND.IS.EQ. 2.AND.IB.EQ.-1.AND.IR.NE.0)    N=34 
      IF(IC.EQ. 0.AND.IS.EQ. 2.AND.IB.EQ.-1.AND.IR.NE.0)    N=35 
      IF(IC.EQ. 1.AND.IS.EQ. 3.AND.IB.EQ.-1.AND.IR.EQ.0)    N=36 
      IF(IC.EQ. 1.AND.IS.EQ. 0.AND.IB.EQ.1.AND.IR.EQ.0)     N=37 
      IF(IC.EQ. 0.AND.IS.EQ. 0.AND.IB.EQ.1.AND.IR.EQ.0)     N=38 
      IF(IC.EQ. 1.AND.IS.EQ.-1.AND.IB.EQ.1.AND.IR.EQ.0)     N=39 
      IF(IC.EQ.-1.AND.IS.EQ.-1.AND.IB.EQ.1.AND.IR.EQ.0)     N=40 
      IF(IC.EQ. 0.AND.IS.EQ.-1.AND.IB.EQ.1.AND.IR.EQ.0)     N=41 
      IF(IC.EQ.-1.AND.IS.EQ.-2.AND.IB.EQ.1.AND.IR.EQ.0)     N=42 
      IF(IC.EQ. 0.AND.IS.EQ.-2.AND.IB.EQ.1.AND.IR.EQ.0)     N=43 
      IF(IC.EQ. 0.AND.IS.EQ.-1.AND.IB.EQ.1.AND.IR.EQ.0)     N=44 
      IF(IC.EQ. 2.AND.IS.EQ. 0.AND.IB.EQ.1.AND.IR.NE.0)     N=45 
      IF(IC.EQ. 1.AND.IS.EQ. 0.AND.IB.EQ.1.AND.IR.NE.0)     N=46 
      IF(IC.EQ.-1.AND.IS.EQ. 0.AND.IB.EQ.1.AND.IR.NE.0)     N=47 
      IF(IC.EQ. 0.AND.IS.EQ. 0.AND.IB.EQ.1.AND.IR.NE.0)     N=48 
      IF(IC.EQ. 1.AND.IS.EQ.-1.AND.IB.EQ.1.AND.IR.NE.0)     N=49 
      IF(IC.EQ.-1.AND.IS.EQ.-1.AND.IB.EQ.1.AND.IR.NE.0)     N=50 
      IF(IC.EQ. 0.AND.IS.EQ.-1.AND.IB.EQ.1.AND.IR.NE.0)     N=51 
      IF(IC.EQ.-1.AND.IS.EQ.-2.AND.IB.EQ.1.AND.IR.NE.0)     N=52 
      IF(IC.EQ. 0.AND.IS.EQ.-2.AND.IB.EQ.1.AND.IR.NE.0)     N=53 
      IF(IC.EQ.-1.AND.IS.EQ.-3.AND.IB.EQ.1.AND.IR.EQ.0)     N=54 
C 
C  VARIABLE MASS NONSRTANGE RESONANSES ===>RO0 
c-Sob      IF(IC.EQ. 0.AND.IS.EQ. 0.AND.IB.EQ.0.AND.IR.NE.0.AND.N.EQ.0) 
c-Sob     *write(16,1001) M,IC 
c-Sob 1001 FORMAT(1X, 'PANUN: GETTING RHO WITH MASS=',F7.3,', CHARGE=',I3) 
      IF(IC.EQ. 0.AND.IS.EQ. 0.AND.IB.EQ.0.AND.IR.NE.0.AND.N.EQ.0) 
     *                                                      N=16 
C 
C  VARIABLE MASS NONSTRANGE PARTICLE    ===>PI0 
c-Sob      IF(IC.EQ. 0.AND.IS.EQ. 0.AND.IB.EQ.0.AND.IR.EQ.0.AND.N.EQ.0) 
c-Sob     *write(16,1002) M,IC 
c-Sob 1002 FORMAT(1X, 'PANUN: GETTING PI0 WITH MASS=',F7.3,', CHARGE=',I3) 
      IF(IC.EQ. 0.AND.IS.EQ. 0.AND.IB.EQ.0.AND.IR.EQ.0.AND.N.EQ.0) 
     *                                                      N= 7 
C 
c-Sob      IF(N.EQ.0)   write(16,1000) IP,M 
c-Sob 1000 FORMAT(2X,'PANUN: IP=',5I5,2X,'M=',E13.6) 
      RETURN 
      END 
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
      FUNCTION TAUN(IK) 
      IMPLICIT REAL*8 (A-H,O-Z) 
      REAL*8 TAUN 
      COMMON /DATA1/CMIX(6,2),PMASM(18),PMASB(18),PGAMM(18), 
     *PGAMB(18),MESO(9,2) 
      IF(IK.LE.18)   GO  TO  2 
      IF(IK.GE.37)   GO  TO  1 
      IF(PGAMB(IK-18).EQ.0.)  GO  TO  3 
      TAU0=0.197/PGAMB(IK-18) 
      GO  TO  4 
    1 IF(PGAMB(IK-36).EQ.0.)  GO  TO  3 
      TAU0=0.197/PGAMB(IK-36) 
      GO  TO  4 
    2 IF(PGAMM(IK).EQ.0.)     GO  TO  3 
      TAU0=0.197/PGAMM(IK) 
      GO  TO  4 
    3 TAUN=10.0E+9 
      RETURN 
    4 DRND=RNDMD(-1) 
      TAUN=-TAU0*DLOG(DRND) 
      IF(TAUN.LE.0.001)     TAUN=0.0011 
      IF(TAUN.GT.1.0D+3)    TAUN=1.0D+3 
      RETURN 
      END 
