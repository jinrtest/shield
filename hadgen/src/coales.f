
************************************************************************
*                                                                      *
*  The following subroutines are written by K.K.Gudima and V.D.Toneev  *
*  Joint Institute for Nuclear Research, Dubna,                        *
*  E-mail: gudima@acad.moldnet.md                                      *
*          gudima@cc.acad.md                                           *
*  E-mail: toneev@thsun1.jinr.ru                                       *
*                                                                      *
************************************************************************

      SUBROUTINE  COALES(MV)         ! sobol, 15.02.96
      IMPLICIT REAL*8 (A-H,O-Z)      ! sobol, 15.02.96
      COMMON /MEMAGT/ PME(9,5999),IME(5,5999)
      COMMON /IDPME/ IDPME(5999)
      COMMON /PORIG/ IORI(3,5999)
      DIMENSION ID1(200),ID2(200),IT1(50),IT2(50),IT3(50),
     *IA1(50),IA2(50),IA3(50),IA4(50),V(3),PJS(3),PJL(3)
      DIMENSION  P1L(3),P1S(3),P2L(3),P2S(3),P3L(3),P3S(3)
      DIMENSION  IND(4)
C-----------------------------
      P0D=0.090     ! Recommended value P0D=0.090
      P0T=0.108     ! Recommended value P0T=0.108
      P0A=0.115     ! Recommended value P0A=0.115
C-----------------------------
      IF(MV.LE.1)   RETURN
      ND=0
           I=1
   10 IF(IME(4,I).EQ.0.OR.IME(2,I).NE.0.OR.IME(3,I).NE.0) GO TO 12
      J=I+1
   11 IF(J.GT.MV)   GO  TO  12
      IF(IME(4,J).EQ.0.OR.IME(2,J).NE.0.OR.IME(3,J).NE.0) GO TO 13
      I1=I
      P1=SQRT(PME(8,I1)*(PME(8,I1)+1.88))
                                  E1=PME(8,I1)+0.940
      PJ=SQRT(PME(8,J)*(PME(8,J)+1.88))
                                EJ=PME(8,J)+0.940
      PJL(1)=PME(4,J)
      PJL(2)=PME(5,J)
      PJL(3)=PME(6,J)
      ES = EJ+E1
      P1L(1)=PME(4,I1)
      P1L(2)=PME(5,I1)
      P1L(3)=PME(6,I1)
      V(1)=-(P1L(1)+PJL(1))/ES
      V(2)=-(P1L(2)+PJL(2))/ES
      V(3)=-(P1L(3)+PJL(3))/ES
      CALL  CINEMC(PJL,V,PJS,CTJ,STJ,CFJ,SFJ,TJS,0.940D00)
c      print *,'D param=',SQRT(TJS*(TJS+1.88))
      IF(SQRT(TJS*(TJS+1.88)).GT.P0D)   GO  TO  13
      IF(ND.LT.200)   GO  TO  50
      write(16,49)  ND
   49 FORMAT(/20X,'ND=',I5)
      GO TO 14
   50 CONTINUE
      ND=ND+1
            ID1(ND)=I
                      ID2(ND)=J
      IME(2,I)=ND
                IME(2,J)=ND
   12 IF(I.GT.(MV-1))   GO  TO  14
      I=I+1
            GO  TO  10
   13 J=J+1
            GO  TO  11
   14 NA=0
      IF(ND.LT.2)  GO TO 53
      ND1=ND-1
      DO 52 ID=1,ND1
      I1=ID1(ID)
               I2=ID2(ID)
                        IF(I1.EQ.0)  GO TO 52
      IB1=ID+1
      DO 51 IB=IB1,ND
      I3=ID1(IB)
      J=ID2(IB)
      IF(I3.EQ.0)   GO  TO  51
      IQ=IME(1,I1)+IME(1,I2)+IME(1,I3)+IME(1,J)
      IF(IQ.NE.2)  GO TO 51
      P1=SQRT(PME(8,I1)*(PME(8,I1)+1.88))
      E1=PME(8,I1)+0.940
      P2=SQRT(PME(8,I2)*(PME(8,I2)+1.88))
      E2=PME(8,I2)+0.940
      P3=SQRT(PME(8,I3)*(PME(8,I3)+1.88))
      E3=PME(8,I3)+0.940
      PJ=SQRT(PME(8,J)*(PME(8,J)+1.88))
      EJ=PME(8,J)+0.940
      P1L(1)=PME(4,I1)
      P1L(2)=PME(5,I1)
      P1L(3)=PME(6,I1)
      P2L(1)=PME(4,I2)
      P2L(2)=PME(5,I2)
      P2L(3)=PME(6,I2)
      P3L(1)=PME(4,I3)
      P3L(2)=PME(5,I3)
      P3L(3)=PME(6,I3)
      PJL(1)=PME(4,J)
      PJL(2)=PME(5,J)
      PJL(3)=PME(6,J)
      ES=EJ+E1+E2+E3
      V(1)=-(P1L(1)+P2L(1)+P3L(1)+PJL(1))/ES
      V(2)=-(P1L(2)+P2L(2)+P3L(2)+PJL(2))/ES
      V(3)=-(P1L(3)+P2L(3)+P3L(3)+PJL(3))/ES
      CALL  CINEMC(P1L,V,P1S,CT1,ST1,CF1,SF1,T1S,0.940D00)
      CALL  CINEMC(P2L,V,P2S,CT2,ST2,CF2,SF2,T2S,0.940D00)
      CALL  CINEMC(P3L,V,P3S,CT3,ST3,CF3,SF3,T3S,0.940D00)
      CALL  CINEMC(PJL,V,PJS,CTJ,STJ,CFJ,SFJ,TJS,0.940D00)
      IF(SQRT(T1S*(T1S+1.88)).GT.P0A)   GO  TO  51
      IF(SQRT(T2S*(T2S+1.88)).GT.P0A)   GO  TO  51
      IF(SQRT(T3S*(T3S+1.88)).GT.P0A)   GO  TO  51
      IF(SQRT(TJS*(TJS+1.88)).GT.P0A)   GO  TO  51
      IF(NA.LT.50)  GO TO 48
      write(16,'(20X,''NA='',I5)')  NA
      GO TO 51
   48 CONTINUE
      NA=NA+1
      IA1(NA)=I1
      IA2(NA)=I2
      IA3(NA)=I3
      IA4(NA)=J
      IME(2,I1)=NA
      IME(2,I2)=NA
      IME(2,I3)=NA
      IME(2,J)=NA
      ID1(ID)=0
      ID2(ID)=0
      ID1(IB)=0
      ID2(IB)=0
C
      GO  TO  52
   51 CONTINUE
   52 CONTINUE
   53 NT=0
      IF(ND.EQ.0)  RETURN
      DO  16  ID=1,ND
      I1=ID1(ID)
      I2=ID2(ID)
      IF(I1.EQ.0)  GO TO16
      DO  15   J=1,MV
      IF(IME(4,J).EQ.0.OR.IME(2,J).NE.0.OR.IME(3,J).NE.0) GO TO 15
      IQ=IME(1,I1)+IME(1,I2)+IME(1,J)
      IF(IQ.EQ.0.OR.IQ.EQ.3)   GO  TO   15
      P1=SQRT(PME(8,I1)*(PME(8,I1)+1.88))
      E1=PME(8,I1)+0.940
      P2=SQRT(PME(8,I2)*(PME(8,I2)+1.88))
      E2=PME(8,I2)+0.940
      PJ=SQRT(PME(8,J)*(PME(8,J)+1.88))
      EJ=PME(8,J)+0.940
      P1L(1)=PME(4,I1)
      P1L(2)=PME(5,I1)
      P1L(3)=PME(6,I1)
      P2L(1)=PME(4,I2)
      P2L(2)=PME(5,I2)
      P2L(3)=PME(6,I2)
      PJL(1)=PME(4,J)
      PJL(2)=PME(5,J)
      PJL(3)=PME(6,J)
      ES = EJ+E1+E2
      V(1)=-(P1L(1)+P2L(1)+PJL(1))/ES
      V(2)=-(P1L(2)+P2L(2)+PJL(2))/ES
      V(3)=-(P1L(3)+P2L(3)+PJL(3))/ES
      CALL  CINEMC(P1L,V,P1S,CT1,ST1,CF1,SF1,T1S,0.940D00)
      CALL  CINEMC(P2L,V,P2S,CT2,ST2,CF2,SF2,T2S,0.940D00)
      CALL  CINEMC(PJL,V,PJS,CTJ,STJ,CFJ,SFJ,TJS,0.940D00)
      IF(SQRT(T1S*(T1S+1.88)).GT.P0T)   GO  TO  15
      IF(SQRT(T2S*(T2S+1.88)).GT.P0T)   GO  TO  15
      IF(SQRT(TJS*(TJS+1.88)).GT.P0T)   GO  TO  15
      IF(NT.LT.50)  GO  TO  55
      write(16,'(20X,''NT='',I5)')  NT
      GO  TO  15
   55 CONTINUE
      NT=NT+1
      IT1(NT)=I1
      IT2(NT)=I2
      IT3(NT)=J
      IME(2,I1)=NT
      IME(2,I2)=NT
      IME(2,J)=NT
      ID1(ID)=0
      ID2(ID)=0
      GO  TO  16
   15 CONTINUE
   16 CONTINUE
   17 IF(NT.EQ.0)   GO  TO  119
      DO  19  IT=1,NT
      DO  18   J=1,MV
      IF(IME(4,J).EQ.0.OR.IME(2,J).NE.0.OR.IME(3,J).NE.0) GO TO 18
      I1=IT1(IT)
      I2=IT2(IT)
      I3=IT3(IT)
      IQ=IME(1,I1)+IME(1,I2)+IME(1,I3)+IME(1,J)
      IF(IQ.NE.2)   GO  TO  18
      P1=SQRT(PME(8,I1)*(PME(8,I1)+1.88))
      E1=PME(8,I1)+0.940
      P2=SQRT(PME(8,I2)*(PME(8,I2)+1.88))
      E2=PME(8,I2)+0.940
      P3=SQRT(PME(8,I3)*(PME(8,I3)+1.88))
      E3=PME(8,I3)+0.940
      PJ=SQRT(PME(8,J)*(PME(8,J)+1.88))
      EJ=PME(8,J)+0.940
      P1L(1)=PME(4,I1)
      P1L(2)=PME(5,I1)
      P1L(3)=PME(6,I1)
      P2L(1)=PME(4,I2)
      P2L(2)=PME(5,I2)
      P2L(3)=PME(6,I2)
      P3L(1)=PME(4,I3)
      P3L(2)=PME(5,I3)
      P3L(3)=PME(6,I3)
      PJL(1)=PME(4,J)
      PJL(2)=PME(5,J)
      PJL(3)=PME(6,J)
      ES = EJ+E1+E2+E3
      V(1)=-(P1L(1)+P2L(1)+P3L(1)+PJL(1))/ES
      V(2)=-(P1L(2)+P2L(2)+P3L(2)+PJL(2))/ES
      V(3)=-(P1L(3)+P2L(3)+P3L(3)+PJL(3))/ES
      CALL  CINEMC(P1L,V,P1S,CT1,ST1,CF1,SF1,T1S,0.940D00)
      CALL  CINEMC(P2L,V,P2S,CT2,ST2,CF2,SF2,T2S,0.940D00)
      CALL  CINEMC(P3L,V,P3S,CT3,ST3,CF3,SF3,T3S,0.940D00)
      CALL  CINEMC(PJL,V,PJS,CTJ,STJ,CFJ,SFJ,TJS,0.940D00)
      IF(SQRT(T1S*(T1S+1.88)).GT.P0A)   GO  TO  18
      IF(SQRT(T2S*(T2S+1.88)).GT.P0A)   GO  TO  18
      IF(SQRT(T3S*(T3S+1.88)).GT.P0A)   GO  TO  18
      IF(SQRT(TJS*(TJS+1.88)).GT.P0A)   GO  TO  18
      IF(NA.LT.50)  GO TO 57
      write(16,'(20X,''NA*='',I5)')  NA
      GO TO 18
   57 CONTINUE
      NA=NA+1
      IA1(NA)=I1
      IA2(NA)=I2
      IA3(NA)=I3
      IA4(NA)=J
      IME(2,I1)=NA
      IME(2,I2)=NA
      IME(2,I3)=NA
      IME(2,J)=NA
      IT1(IT)=0
      IT2(IT)=0
      IT3(IT)=0
      GO  TO  19
   18 CONTINUE
   19 CONTINUE
  119 IF(NA.EQ.0)   GO  TO  21
C--------------------------------------
      DO  20  IA=1,NA
      I1=IA1(IA)
      I2=IA2(IA)
      I3=IA3(IA)
      I4=IA4(IA)
      IND(1)=I1
      IND(2)=I2
      IND(3)=I3
      IND(4)=I4
      CALL  CODIR(IND,4)
      IME(1,I1)=2
      IME(2,I1)=0
      IME(4,I1)=4
      IDPME(I1)=0
      IORI(1,I1)=0
      IORI(2,I1)=0
      IORI(3,I1)=4
   20 CONTINUE
   21 IF(NT.EQ.0)   GO  TO  23
C--------------------------------------
      DO  22  IT=1,NT
      IF(IT1(IT).EQ.0)   GO  TO  22
      I1=IT1(IT)
      I2=IT2(IT)
      I3=IT3(IT)
      IND(1)=I1
      IND(2)=I2
      IND(3)=I3
      CALL  CODIR(IND,3)
      IME(1,I1)=IME(1,I1)+IME(1,I2)+IME(1,I3)
      IME(2,I1)=0
      IME(4,I1)=3
      IDPME(I1)=0
      IORI(1,I1)=0
      IORI(2,I1)=0
      IORI(3,I1)=4
   22 CONTINUE
C--------------------------------------
   23 DO  25  ID=1,ND
      IF(ID1(ID).EQ.0)   GO  TO  25
      I1=ID1(ID)
      I2=ID2(ID)
      IF((IME(1,I1)+IME(1,I2)).NE.1)   GO  TO  24
      IND(1)=I1
      IND(2)=I2
      CALL  CODIR(IND,2)
      IME(1,I1)=1
      IME(2,I1)=0
      IME(4,I1)=2
      IDPME(I1)=0
      IORI(1,I1)=0
      IORI(2,I1)=0
      IORI(3,I1)=4
C--------------------------------------
      GO  TO  25
   24 IME(2,I1)=0
      IME(2,I2)=0
   25 CONTINUE
      I=1
   26 IF(IME(2,I).EQ.0)   GO  TO  29
      DO  27  K=1,9
   27 PME(K,I)=PME(K,MV)
      DO  28  K=1,5
      IF(K.LE.3)  IORI(K,I)=IORI(K,MV)
   28 IME(K,I)=IME(K,MV)
      IDPME(I)=IDPME(MV)
      MV=MV-1
      IF(I.GT.MV)   GO  TO  30
      GO  TO  26
   29 IF(I.GE.MV)   GO  TO  30
      I=I+1
            GO  TO  26
   30 CONTINUE
      RETURN
      END



      SUBROUTINE  CODIR(IND,NN)
      IMPLICIT REAL*8 (A-H,O-Z)     ! sobol, 15.02.96
      COMMON /MEMAGT/ PME(9,5999),IME(5,5999)
      DIMENSION  PS(3),IND(4)
      PS(1)=0.
      PS(2)=0.
      PS(3)=0.
      DO  K=1,NN
      I=IND(K)
      PI=SQRT(PME(8,I)*(PME(8,I)+1.88))
      PS(1)=PS(1)+PME(4,I)
      PS(2)=PS(2)+PME(5,I)
      PS(3)=PS(3)+PME(6,I)
      ENDDO
      PSM=SQRT(PS(1)**2+PS(2)**2+PS(3)**2)
      I1=IND(1)
      PMM=FLOAT(NN)*0.940
      PME(8,I1)=SQRT(PSM*PSM+PMM*PMM)-PMM
      PME(9,I1)=PMM
      PME(4,I1)=PS(1)
      PME(5,I1)=PS(2)
      PME(6,I1)=PS(3)
      PME(7,I1)=0.
      RETURN
      END



      SUBROUTINE CINEMC(PSTAR,V,P,CT,ST,CFI,SFI,T,CM)
C     KINEMATIC BLOCK.
      IMPLICIT REAL*8 (A-H,O-Z)
C      REAL *4 PSTAR,V,P,CT,ST,CFI,SFI,T,CM
      DIMENSION PSTAR(3),V(3),P(3)
C
      SPV = PSTAR(1)*V(1)+PSTAR(2)*V(2)+PSTAR(3)*V(3)
      V2=V(1)**2+V(2)**2+V(3)**2
      G=1.
      V2S=DSQRT(DABS(1.-V2))
      IF(V2S.GT.0.D0)  G=1./DSQRT(DABS(1.-V2))
      ES2=PSTAR(1)**2+PSTAR(2)**2+PSTAR(3)**2+CM**2
      ES=DSQRT(ES2)
C
      do   K=1,3
        P(K)=PSTAR(K)+G*V(K)*(SPV*G/(G+1.)+ES)
      enddo
C
      E=G*(ES+SPV)
      T=E-CM
      IF(T.LE.0.)  T=1.E-6
      PM2=P(1)**2+P(2)**2+P(3)**2
      PM=DSQRT(PM2)
      CT=P(3)/PM
      ST2=1.-CT*CT
C
      IF(ST2.GT.0.)  then
        ST=DSQRT(ST2)
        CFI=P(1)/PM/ST
        SFI=P(2)/PM/ST
      ELSE
        ST=0.
        CFI=1.
        SFI=0.
      ENDIF
C
      RETURN
       END
