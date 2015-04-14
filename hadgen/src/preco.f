
************************************************************************
*                                                                      *
*  The following subroutines are written by K.K.Gudima and V.D.Toneev  *
*  Joint Institute for Nuclear Research, Dubna,                        *
*  E-mail: gudima@acad.moldnet.md                                      *
*          gudima@cc.acad.md                                           *
*  E-mail: toneev@thsun1.jinr.ru                                       *
*                                                                      *
************************************************************************

      SUBROUTINE PRECO(ENEXT,ATWGHT,CHARGE,PNX,PNY,PNZ,
     *                 AM,AMF,RADNCL,KSMEMO,N0,P0,H0,PZ0,K1)
c Input and Output - in GeV, internal calculations - in MeV (Sobol,25.7.95)
c Input:
c       ENEXT,ATWGHT,CHARGE - input
c       PNX,PNY,PNZ - input and output
c       AM,AMF,RADNCL - input only
c       KSMEMO - input and output            ! former KSTART !
c       N0 - input
c       P0 - input
c       H0 - input
c       PZ0 - input
c       K1 - output (don't used in HADGEN)
c   PRECOM - input,  this parametr turn on the precompound stage
C**********************************************************************
C  THE MAIN CONTROL SUBROUTINE                                        *
C**********************************************************************
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      COMMON /BK1001/ T1Y(130)
      COMMON /BK1002/ T2XY(200)
      COMMON /BK1005/ AJ(6)
      COMMON /BK1006/ ZJ(6)
      COMMON /BK1008/ DLM(6)
      COMMON /BK1014/ GAM(6)
      COMMON /BK1016/ CC(6)
      COMMON /BK1017/ VK(6)
      COMMON /BK1011/ VJ(6)
      COMMON /BK1015/ RJ(7)
        COMMON /BK1003/ U,A,Z     ! Output U,A,Z.
      COMMON /BK1009/ AFJ(7)
      COMMON /BK1010/ ZFJ(6)
        COMMON /LEVERS/ PARLEV(40)  ! Parameters, see Subr. INLEVR
C
      COMMON /PREANG/ ANGL(4)       ! former /BLANGL/!
C
      COMMON /MEMAGT/ PMEMO(9,5999),IMEMO(5,5999)
      REAL*8 PMEMO
        DIMENSION GJ(7)
        DIMENSION PS(3),PL(3),V(3)
      COMMON /BR0J/ R0J(6),BN
      COMMON /BLBJ/ BJ(7)
      COMMON /BLALJ/ ALJ(6)
      COMMON /BLR0/ R0
      COMMON /BLEXN/ EXN
      COMMON /BLAC/ AC
C -----------------------------------------------------------------
        R0=RADNCL
      U=1000.*ENEXT
      A=ATWGHT
      Z=CHARGE
C
cC Debug printout 4.8.95 =================================
c      IF(A.LE.0. .OR. Z.LE.0. .OR. U.LE.0.)THEN
c      CONTINUE!WRITE(25,237)U,A,Z,PNX,PNY,PNZ
c  237 FORMAT(' Invalid PRECO input: U=',F7.2,' A=',F6.1,' Z=',F6.1,
c     *       '  P(vect) - ',3F7.2)
c      END IF
cC =============end debug==================================
C
C Check of PRECO Input --------
      IF(A.LT.0. .OR. Z.LT.0. .OR. U.LT.0.)THEN
        CONTINUE!WRITE(25,238)U,A,Z,PNX,PNY,PNZ
  238     FORMAT(' Invalid PRECO input: U=',F7.2,' A=',F6.1,' Z=',F6.1,
     *           '  P(vect) - ',3F7.2)
        IF(A.LT.0.)A=0.
        IF(Z.LT.0.)Z=0.
        IF(U.LT.0.)U=0.
      END IF
      IF(A.LT.0.1)RETURN   ! Full desintegration will be diagnosed in DEEX
C                          ! (A=1,2,... and Z=0 is OK)
C ------- End of check --------
C
        E=SQRT(PNX**2+PNY**2+PNZ**2+(.94*A)**2)
        V(1)=PNX/E
        V(2)=PNY/E
        V(3)=PNZ/E
      K1=KSMEMO
      PZ=PZ0                          ! charge of excited particles!
C
    3 N=N0
      P=P0
      H=H0
      IF(U.LT.1.)GO TO 7              ! 1 MeV? !
      IF(KSMEMO.GT.5999)GO TO 99
      KST0=KSMEMO
C ------------------------- Main loop ------------------
      DO 48 K=KST0,5999
      IF(A-4.)    101,101,100
  100 IF(Z-2.) 101,101,8
C ----Lightest nucleus------
  101 IF(A-1.1)  102,103,103
C ----Nucleon---------------
  102 PMEMO(4,K)=PNX
      PMEMO(5,K)=PNY
      PMEMO(6,K)=PNZ
      PMEMO(7,K)=0.
      PMEMO(8,K)=SQRT(PNX**2+PNY**2+PNZ**2+0.94**2)-0.94
      PMEMO(9,K)=0.94
      IF(Z-0.1)  111,112,112
  111 IMEMO(1,K)=0
                   GO  TO  113
  112 IMEMO(1,K)=1
  113 IMEMO(2,K)=0
      IMEMO(3,K)=0
      IMEMO(4,K)=1
      KSMEMO=KSMEMO+1
C
C_____Sobolevsky, 16.12.96:_________________________________________
C                       As we have written the nucleon to MEMORY,   !
C                       we should delete this nucleon from residual !
C                       nucleus COMMON /BK1003/ U,A,Z.              !
      U=0.  ! 16.12.96                                              !
      A=0.  ! 16.12.96                                              !
      Z=0.  ! 16.12.96                                              !
        PNX=0.   ! 23.12.96                                       !
        PNY=0.   ! 23.12.96                                       !
        PNZ=0.   ! 23.12.96                                       !
C___________________________________________________________________!
C
      RETURN           ! Memorizing of nucleon and return !
C
  103 IF(A.EQ.4. .AND. Z.EQ.2. .AND. U.LT.20.)RETURN   ! to compound deex.
      IF(A.EQ.3. .AND. Z.GT.0. .AND. U.LT.6.0)RETURN   ! to compound deex.
C
C Simplest statistical break-up of highly excited lightest nuclei (Toneev)
        EP1=U/A
        U=U-EP1
      IF(Z-1.)104,105,105
  104 EP2=0.
            GO TO 106
  105 EP2=1.
  106 Z=Z-EP2
      EP3=940.
      A=A-1.
      GO TO 107
C
C ------------------- Nucleus A>4, Z>2 -------------------
    8 DL=PREDEL(A,Z)                     ! former DELTAM(A,Z) !
      CALL VHELP
C
      DO 9 I=1,6
        AFJ(I)=A-AJ(I)
        ZFJ(I)=Z-ZJ(I)
        VJ(I)=PRECOL(I,RADNCL,AM)          ! former COLOMB !
        BJ(I)=PREDEL(AFJ(I),ZFJ(I))-(DL-DLM(I))
    9 RJ(I)=(U-BJ(I))/(1.+AJ(I)/A)-VJ(I)
C
      R0J(1)=0.76+2.2/AFJ(1)**0.33333
      BN=(2.12/AFJ(1)**0.666667-0.05)/R0J(1)
C
      DO 1000 IR=2,6
 1000 R0J(IR)=1.+CC(IR)
C
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C Parametr PRECOM (0<PRECOM<1) runs the precompound stage:      *
C          PRECOM=0 - precompound is turned off                 *
C          PRECOM=1 - precompound is completely turned on       *
C                                                               *
      PRECOM=PARLEV(22)                     ! see Subr. INLEVR
      NSP=NINT(SQRT((1.19*AM*A*U+0.5)*PRECOM))
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
   10 IF(N.lt.NSP)then
        go to 11
      else
        RETURN       ! all output variables should be defined !
      end if
C
C  **************** PRE-EQUILIBRIUM EMISSION *************************
   11 CONTINUE
      IF(P.GE.1.) GO TO 700
      P=P+1.
      H=H+1.
      N=N+2
  700 CONTINUE
      EXN=P+H
      ALJ(1)=P*(EXN-1.)
      ALJ(2)=ALJ(1)
      ALJ(3)=ALJ(2)*(P-1.)*(EXN-2.)/2.
      ALJ(4)=ALJ(3)*(P-2.)*(EXN-3.)/6.
      ALJ(5)=ALJ(4)
      ALJ(6)=ALJ(5)*(P-3.)*(EXN-4.)/12.
      AC=.595*AM
      CALL GBETA(A)
      DO 16 L=1,6
      IF(P-AJ(L)+.01) 12,12,13
   12 GJ(L)=0.
               GO  TO  16
   13 IF(P+H-AJ(L)-.01)  12,12,50
   50 IF(PZ-ZJ(L)+.01)  12,12,14
   14 IF(RJ(L)) 12,12,15
   15 GJ(L)=GAMAGU(L)
   16 CONTINUE
      G=0.
      DO 17 I=1,6
   17 G=G+GJ(I)
   18 IF(N) 19,19,21
   19 CONTINUE!WRITE(25,123)N
  123 FORMAT(10X,'PRECO: NUMBER OF EXCITONS=',I3)
      RETURN
   21 AC=.595*AM
      CALL TRANS(P,H,AC,C1,C2,C3)
      C=C1+C2+C3
      B1=RRAN(IX)
      B2=C+G
      IF(B1-G/B2)  29,29,22
   22 C11=C1/C
      C21=(C2+C1)/C
      B3=RRAN(IX)
      IF(B3-C11) 23,23,25
   23 N = N+2
      P = P+1.
      H = H+1.
      BZ =RRAN(IX)
      BZ1 = Z/A
      IF(BZ-BZ1)  24,24,10
   24 PZ = PZ+1.
               GO TO 10
   25 IF(B3-C21) 26,26,10
   26 N = N-2
      P = P-1.
      H = H-1.
      IF(PZ) 10,10,27
   27 BZ =RRAN(IX)
      BZ1 = Z/A
      IF(BZ-BZ1)  28,28,10
   28 PZ=PZ-1.
             GO TO 10
   29 DO 30 J=2,6
   30 GJ(J)=GJ(J-1)+GJ(J)
      B=RRAN(IX)*G
      DO 32 J=1,6
      IF(B-GJ(J))  31,31,32
   31 LM=J
           GO TO 33
   32 CONTINUE
   33 EP1=TKINM(LM,P,H)
      EP2=ZJ(LM)
      EP3=940.*AJ(LM)
      K1=K1+1
            GO TO 44
C  *******************************************************************
   44 P=P-AJ(LM)
      N1=AJ(LM)
      N=N-N1
      PZ=PZ-EP2
        U=(RJ(LM)-EP1+VJ(LM))*(1.+AJ(LM)/A)
        A=AFJ(LM)
        Z=ZFJ(LM)
  107 CALL PREISA                                 ! former ISANGL !
      PM=SQRT(EP1*(EP1+2.*EP3))
        PS(1)=PM*ANGL(4)*ANGL(3)
        PS(2)=PM*ANGL(4)*ANGL(2)
        PS(3)=PM*ANGL(1)
      IF(V(1)**2+V(2)**2+V(3)**2-0.00000001)  108,108,109
  108 PMEMO(4,K)=PS(1)/1000.
      PMEMO(5,K)=PS(2)/1000.
      PMEMO(6,K)=PS(3)/1000.
      PMEMO(7,K)=0.
      PMEMO(8,K)=EP1/1000.
      PMEMO(9,K)=EP3/1000.
        IMEMO(1,K)=EP2
        IMEMO(2,K)=0
        IMEMO(3,K)=0
        IMEMO(4,K)=EP3/939.
        PL(1)=PS(1)
        PL(2)=PS(2)
        PL(3)=PS(3)
      GO TO 47
C
  109 CONTINUE
        CALL PRECIN(PS,V,PL,COT,ST,CF,SF,TL,EP3)  ! former CINEMA !
      PMEMO(4,K)=PL(1)/1000.
      PMEMO(5,K)=PL(2)/1000.
      PMEMO(6,K)=PL(3)/1000.
      PMEMO(7,K)=0.
      PMEMO(8,K)=TL/1000.
      PMEMO(9,K)=EP3/1000.
        IMEMO(1,K)=EP2
        IMEMO(2,K)=0
        IMEMO(3,K)=0
        IMEMO(4,K)=EP3/939.
C
   47 PNX=PNX-PL(1)/1000.
      PNY=PNY-PL(2)/1000.
      PNZ=PNZ-PL(3)/1000.
        E=SQRT(PNX**2+PNY**2+PNZ**2+(.94*A)**2)
        V(1)=PNX/E
        V(2)=PNY/E
        V(3)=PNZ/E
      KSMEMO=KSMEMO+1
C
cC Debug printout =========================   
c      CONTINUE!WRITE(25,377) K
c  377 format(' Main Loop ',I3)
c      CONTINUE!WRITE(25,378) U,A,Z
c  378 format(' U,A,Z:  ',F7.2,2F6.1)
c      CONTINUE!WRITE(25,379) PMEMO(8,K),PMEMO(9,K),
c     *              IMEMO(1,K),IMEMO(3,K),IMEMO(4,K)
c  379 format(' T,M; e S B:  ',F7.5,F7.3,4X,3I3)
cC =================end debug===============
   48 CONTINUE
C ------------------- End of main loop ----------------
C
   99 CONTINUE!WRITE(25,49)U,A,Z
   49 FORMAT(10X,'MEMAGT IS EXCEEDED IN PRECO,',
     *       '  U=',F10.5,'  A=',F5.1,'  Z=',F4.1)
    7 RETURN
      END
C
C  ********* last correction 11.16.1994 ****************************
C
      SUBROUTINE GBETA(A)
      COMMON /GAMBET/ GB(6)
      GB(1)=1.
      GB(2)=1.
      GB(3)=16./A
      GB(4)=243./A**2
      GB(5)=GB(4)
      GB(6)=4096./A**3
      RETURN
      END


      FUNCTION PREDEL(X,Y)           ! former DELTAM(A,Z) !
C**********************************************************************
C  CALCULATION OF MASS DEFECT                                         *
C**********************************************************************
      COMMON /BK1001/ T1Y(130) 
      COMMON /BK1002/ T2XY(200)
C
c      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
cC X is A - mass number
cC Y is Z - charge number
c      if(X.le.Y .OR. X.eq.0.)then
c        CONTINUE!WRITE(25,1000)NUMTRE,X,Y
c 1000   format(' PREDEL diagnostic: NUMTRE=',I7,
c     *         '   A=',F5.1,'   Z=',F5.1)
c      end if
C
      ES=(25.8357-44.2355*((X-2.*Y)/X)**2)*
     *(((1.-.62025/X**.66667)**2)*(X**.66667))
      EC=.779*((Y*(Y-1.))/X**.33333)*
     *((1.-1.5849/X**.66667)+(1.2273/X+(1.5772/X**1.33333)))
      EALFA=(-.4323*((Y**1.33333)/(X**.33333)))*
     *((1.+.49597/X)-(.57811/X**.33333+(.14518/X**.66667)))
      EDOB=((8.367*X+31.4506*((X-2.*Y)**2/X))-(.783*Y+17.0354*X))
      I=int(X)
      J=int(Y)
      L=I-J
C
C This bypass do not influence on the results (Sobol, 25.04 2003)
      if(L.le.0)then
        T1=0.
        T2=0.
      else
        T1=T1Y(J)
        T2=T2XY(L)
      end if
C
      PREDEL=((ES+EC)+(EALFA+EDOB))+(T1+T2)
      RETURN
      END



      SUBROUTINE PREISA               ! former ISANGL !
C**********************************************************************
C  CHOOSE ISOTROPIC DISTRIBUTED ANGLE FOR PARTICLE EMITTED            *
C**********************************************************************
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      COMMON /PREANG/ ANGL(4)               ! former /BLANGL/ !
      ANGL(1)=1.-2.*RRAN(IX)
      ANGL(4)=SQRT(1.-ANGL(1)**2)
      F=2.*3.1415926536*RRAN(IX)
      ANGL(2)=SIN(F)
      ANGL(3)=COS(F)
      RETURN
      END
      


      SUBROUTINE VHELP
C**********************************************************************
C  AUXILIARY BLOCK FOR NUCLEAR DATA EXTRACTION                        *
C**********************************************************************
      COMMON /BK1016/ CC(6)
      COMMON /BK1017/ VK(6)
      COMMON /BK1003/ U,A,Z
      COMMON /STMASS/ Z1(5),A1(5),C1(5),A2(5),C2(5)
      CC(1)=0.
      VK(1)=0.
      CC(2)=SUBEV(Z,Z1,C1,5)
      VK(2)=SUBEV(Z,Z1,A1,5)
      CC(6)=SUBEV(Z,Z1,C2,5)
      VK(6)=SUBEV(Z,Z1,A2,5)
      CC(3)=CC(2)/2.
      VK(3)=VK(2)+.06
      CC(4)=CC(2)/3.
      VK(4)=VK(2)+.12
      CC(5)=CC(6)*4./3.
      VK(5)=VK(6)-.06
      RETURN
      END



      FUNCTION PRECOL (L,RADNCL,AM)         ! former COLOMB !
C**********************************************************************
C  CALCULATION OF COULOMB ENERGY                                      *
C***********************************************************************
      COMMON /BK1006/ ZJ(6)
      COMMON /BK1010/ ZFJ(6)
      COMMON /BK1009/ AFJ(7)
      COMMON /BK1005/ AJ(6)
      COMMON /BK1017/ VK(6)
      COMMON /BK1003/ U,A,Z
C
      IF(L.le.1) then
      PRECOL=0.
      RETURN
      END IF
C
      TEMP1=VK(L)*(1.44/RADNCL)
      PRECOL=TEMP1*((ZJ(L)*ZFJ(L))/(AJ(L)**.33333+AFJ(L)**.33333))
      PRECOL=PRECOL*(1.-U/(81.*A*AM))
      IF(PRECOL.lt.0.)PRECOL=0.
      RETURN
      END



      FUNCTION SUBEV(U,E,F,N)
C**********************************************************************
C  QUADRATIC INTERPOLATION                                            *
C**********************************************************************
      DIMENSION E(N),F(N)
        IF(U-E(1))1,1,2
    1 X1=E(1)
      X2=E(2)
      X3=E(3)
      Y1=F(1)
      Y2=F(2)
      Y3=F(3)
      GO TO 9
    2     IF(U-E(N-1))6,3,3
    3     IF(U-E(N))5,5,4
    4 SUBEV = F(N)
      RETURN
C
    5 X1=E(N-2)
      X2=E(N-1)
      X3=E(N)
      Y1=F(N-2)
      Y2=F(N-1)
      Y3=F(N)
      GO TO 9
C
    6 DO 8 J=2,N
        IF(U-E(J))7,8,8
    7 X1=E(J-1)
      X2=E(J)
      X3=E(J+1)
      Y1=F(J-1)
      Y2=F(J)
      Y3=F(J+1)
      GO TO 9
    8 CONTINUE
C
    9 SUBEV=Y1*(((U-X2)*(U-X3))/((X1-X2)*(X1-X3)))+
     *      Y2*(((U-X1)*(U-X3))/((X2-X1)*(X2-X3)))+
     *      Y3*(((U-X1)*(U-X2))/((X3-X1)*(X3-X2)))
      RETURN
      END



      SUBROUTINE TRANS(P,H,AM,C1,C2,C3)
C**********************************************************************
C  CALCULATION OF TRANSITION RATES                                    *
C**********************************************************************
      COMMON /BK1003/ U,A,Z
      TF=45.
      AK=1.
      EST=1.6*TF+U/(P+H)
      B1=SQRT(2.*EST/940.)
      B2=B1**2
      SPP=10.63/B2-29.93/B1+42.9
      SPN=34.10/B2-82.20/B1+82.2
      SF = (SPP+SPN)/2.
      B3=TF/EST
        IF(B3.le.0.5)THEN
            T=1.-7.*B3/5.
        ELSE
            T=1.-7.*B3/5.+0.4*B3*((2.-1./B3)**2.5)
        END IF
      SVV=0.00332*SF*T*SQRT(EST)/(AK*((1.2+1./(4.7*B1))**3))
      C1=SVV
      T1=P+H
      GE=AM*A*U
      C2=(C1*P*H*(T1+1.)*(T1-2.))/GE**2
      C3=C1*(T1+1.)*(P*(P+1.)+4.*P*H+H*(H-1.))/(T1*GE)
        IF(C2.lt.0.)C2=0.
C
      RETURN
      END



      SUBROUTINE GU8(A,B,J,F,Y)
      DIMENSION W(8),FIKS(8)
      DATA W/
     * 0.1012285363,  0.2223810345,  0.3137066459,  0.3626837834,
     * 0.3626837834,  0.3137066459,  0.2223810345,  0.1012285363/
      DATA FIKS/
     * 0.9602898565,  0.7966664774,  0.5255324099,  0.1834346425,
     *-0.1834346425, -0.5255324099, -0.7966664774, -0.9602898565/
      EXTERNAL F     ! KarSjh
C
      Y=0.
      DO K=1,8
        E=.5*((B-A)*FIKS(K)+(B+A))
        Y=Y+W(K)*F(J,E)*(B-A)*.5
      END DO
C
      RETURN
      END



      FUNCTION GAMAGU(J)
      COMMON /BLBJ/ BJ(7)
      COMMON /BK1011/ VJ(6)
      COMMON /BK1003/ U,A,Z
      EXTERNAL WBE
      A1=VJ(J)
      IF(J.EQ.1)A1=0.
      B1=U-BJ(J)
      CALL GU8(A1,B1,J,WBE,Y)
      GAMAGU=Y
      RETURN
      END



      FUNCTION TKINM(J,P,H)
C* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C    KINETIC ENERGY FOR PARTICLES IN PRE-EQUILIBRIUM DECAY   *
C* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      COMMON /BK1015/ RJ(7)
      COMMON /BK1005/ AJ(6)
      COMMON /BK1003/ U,A,Z
      COMMON /BK1011/ VJ(6)
      COMMON /BK1009/ AFJ(7)
      COMMON /BLBJ/ BJ(7)
        IF(J-1)1,1,2
    1 DJ=(2.12/AFJ(J)**0.66667-0.05)/(0.76+2.2/AFJ(J)**0.333333)
      GO TO 3
    2 DJ=-VJ(J)
    3 T=P+H-AJ(J)-1.
      R2=RJ(J)
      R1=R2+VJ(J)
        IF(J-2)12,12,13
   12     IF(T+0.1)4,4,5
    4 TKINM=R1
      RETURN
C
    5     IF(T-0.1)10,10,7
   10     IF(J-1)11,11,6
   11 TKINM=-DJ+SQRT(DJ**2+(RRAN(IX)*(R2**2+2.*DJ*R2)))
      RETURN
C
    6 B1=RRAN(IX)
      TKINM=SQRT(B1)*R2+VJ(J)
      RETURN
C
    7 E1=(R1-DJ*T)/(T+1.)
    8 B1=RRAN(IX)
      E=VJ(J)+B1*R2
      T1=(E+DJ)/(E1+DJ)
      T2=(R1-E)/(R1-E1)
      T3=T1*(T2**T)
      B2=RRAN(IX)
        IF(B2-T3) 9,9,8
    9 TKINM=E
      RETURN
C
   13     IF(T+0.1) 4,4,14
   14     IF(T-0.1) 15,15,16
   15 E1=R1
   17 B1=RRAN(IX)
      E=VJ(J)+B1*R2
      PJ1=AJ(J)-1.5
      ABJ=ABS(BJ(J))
        IF(BJ(J).LE.0. .AND. ABJ.GT.VJ(J))E=ABJ+B1*U
      T3=(((E+BJ(J))/(E1+BJ(J)))**PJ1)*((E+DJ)/(E1+DJ))
      B2=RRAN(IX)
        IF(B2-T3) 18,18,17
   18 TKINM=E
      RETURN
C
   16 PJ1=AJ(J)-1.5
      ES=U*(AJ(J)-0.5)+(U-R2)*(P+H-2.5)
      E1=(ES+SQRT(ES*ES-(U-R2)*(AJ(J)-1.5)*(P+H-1.5)*4.*U))/
     *   ((P+H-1.5)*2.)-U+R1
   19 B1=RRAN(IX)
      E=VJ(J)+B1*R2
      ABJ=ABS(BJ(J))
        IF(BJ(J).LE.0. .AND. ABJ.GT.VJ(J))E=ABJ+B1*U
        IF(((R1-E)/(R1-E1)).LT.(10.**(-50./T)))GO TO 19
      T3=(((E+BJ(J))/(E1+BJ(J)))**PJ1)*((E+DJ)/(E1+DJ))*
     *   (((R1-E)/(R1-E1))**T)
      B2=RRAN(IX)
        IF(B2-T3) 20,20,19
   20 TKINM=E
      RETURN
C
      END



      SUBROUTINE PRECIN(PSTAR,V,P,CT,ST,CFI,SFI,T,CM)  ! former CINEMA !
C KINEMATIC BLOCK.
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*4 PSTAR,V,P,CT,ST,CFI,SFI,T,CM
      DIMENSION PSTAR(3),V(3),P(3)
C
      SPV=PSTAR(1)*V(1)+PSTAR(2)*V(2)+PSTAR(3)*V(3)
      V2=V(1)**2+V(2)**2+V(3)**2
      G=1.
      V2S=DSQRT(DABS(1.-V2))
        IF(V2S.GT.0.D0)G=1./DSQRT(DABS(1.-V2))
      ES2=PSTAR(1)**2+PSTAR(2)**2+PSTAR(3)**2+CM**2
      ES=DSQRT(ES2)
        do K=1,3
            P(K)=PSTAR(K)+G*V(K)*(SPV*G/(G+1.)+ES)
        end do
      E=G*(ES+SPV)
      T=E-CM
        IF(T.LE.0.)T=1.E-6
      PM2=P(1)**2+P(2)**2+P(3)**2
      PM=DSQRT(PM2)
      CT=P(3)/PM
      ST2=1.-CT*CT
        IF(ST2.GT.0.)THEN
            ST=DSQRT(ST2)
            CFI=P(1)/PM/ST
            SFI=P(2)/PM/ST
        ELSE
            ST=0.
            CFI=1.
            SFI=0.
        END IF
      RETURN
      END
      


      FUNCTION WBE(J,E)
      COMMON /BK1009/ AFJ(7)
      COMMON /BR0J/ R0J(6),BN
      COMMON /BK1005/ AJ(6)
      COMMON /GAMBET/ GB(6)
      COMMON /BK1011/ VJ(6)
      COMMON /BLBJ/ BJ(7)
      COMMON /BLALJ/ ALJ(6)
      COMMON /BLR0/ R0
      COMMON /BK1003/ U,A,Z
      COMMON /BLEXN/ EXN
      COMMON /BLAC/ AC
C
*      REAL *8 ACOR
*      COMMON /IR/IR
*      EXTERNAL CORREC
*      CALL ERRSAV(208,ACOR)
*      CALL ERRSET(208,256,-1,1,CORREC)
*      IR=0
C
      IF(J-2)12,11,10
C======================== KarSjh =================================
   10 CONTINUE
c      WBE=GB(J)*R0J(J)*0.104/(R0*AFJ(J)**0.333333*SQRT(AJ(J)*U))*
c     *ALJ(J)*((E-VJ(J))/U)*((E+BJ(J))/U)**(AJ(J)-1.5)*
c     *1.
      WBE1=GB(J)*R0J(J)*0.104
      WBE2=(R0*AFJ(J)**0.333333*SQRT(AJ(J)*U))
      WBE3=((E-VJ(J))/U)
      WBE4=((E+BJ(J))/U)
        if(WBE4.lt.0.)WBE4=0.   ! WBE4 can rarely be < 0. 02.05.01
      WBE5=(AJ(J)-1.5)
      WBE=WBE1/WBE2*ALJ(J)*WBE3*WBE4**WBE5
C====================== end KarSjh ================================
      T1=1.-(E+BJ(J))/U
      T2=30./EXN
        IF(T1.LT.(10.**(-T2)))GO TO 14
      WBE=WBE*T1**(30./T2-AJ(J)-1.)
        GO TO 15
   11 C1=E-VJ(J)
        GO TO 13
   12 C1=E+BN
   13 WBE=0.000234*R0*R0*AFJ(J)**.66667*R0J(J)*
     *ALJ(J)/(AC*U*AFJ(J))*C1
      T1=1.-(E+BJ(J))/U
      T2=30./EXN
        IF(T1.LT.(10.**(-T2)))GO TO 14
      WBE=WBE*T1**(30./T2-2.)
        GO TO 15
   14 WBE=0.
   15 CONTINUE
*      IF(IR.EQ.0)RETURN
*      write(16,16) J,U,E,VJ(J),BJ(J),AJ(J),AFJ(J),EXN
*   16 FORMAT(2X,I2,7(2X,E13.6))
*      CALL  ERRSTR(208,ACOR)
      RETURN
      END



      BLOCK DATA PREDAT
C Data for Precompound emission
      COMMON /BK1016/ CC(6) 
      COMMON /BK1017/ VK(6)
      COMMON /BK1005/ AJ(6)  
      COMMON /BK1006/ ZJ(6) 
      COMMON /BK1008/ DLM(6)
      COMMON /BK1014/ GAM(6)
      COMMON /STMASS/ Z1(5),A1(5),C1(5),A2(5),C2(5)
      COMMON /BK1001/ T1Y(130)
      COMMON /BK1002/ T2XY(200)
C
      DATA AJ  /1., 1., 2., 3., 3., 4./
      DATA ZJ  /0., 1., 1., 1., 2., 2./
      DATA DLM /8.368, 7.569, 13.835, 15.835, 15.817, 3.607/
      DATA CC  /0., 0.2, 0.1, 0.07, 0.13, 0.1/
      DATA VK  /0., 0.7, 0.77, 0.8, 0.8, 0.83/
      DATA GAM /1., 1., 3., 3., 3., 2./
C
      DATA Z1 /10., 20., 30., 50., 70./
      DATA A1 /0.42, 0.58, 0.68, 0.77, 0.80/
      DATA C1 /0.50, 0.28, 0.20, 0.15, 0.10/
      DATA A2 /0.68, 0.82, 0.91, 0.97, 0.98/
      DATA C2 /0.10, 0.10, 0.10, 0.08, 0.06/
C
      DATA   T1Y/
     *20.8,15.8,21.,16.8,19.8,16.5,18.8,16.5,18.5,17.2,18.26,15.05,
     *16.01,12.04,13.27,11.09,12.17,10.26,11.04,8.41,9.79,7.36,8.15,
     *5.63,5.88,3.17,3.32,.82,1.83,.97,2.33,1.27,2.92,1.61,2.91,1.35,
     *2.4,.89,1.74,.36,.95,-.65,-.04,-1.73,-.96,-2.87,-2.05,-4.05,-3.4,
     *-5.72,-3.75,-4.13,-2.42,-2.85,-1.01,-1.33,.54,-.02,1.74,.75,2.24,
     *1.,1.98,.79,1.54,.39,1.08,0.,.78,-.35,.58,-.55,.59,-.61,.59,-.35,
     *.32,-.96,-.52,-2.08,-2.46,-3.64,-1.55,-.96,.97,.88,2.37,1.75,2.72,
     *1.9,2.55,1.46,1.93,.86,1.17,.08,.39,-.76,-.39,-1.51,-1.17,-2.36,
     *-1.95,-3.06,-2.62,-3.55,-2.95,-3.75,-3.07,-3.79,-3.06,-3.77,-3.05,
     *-3.78,-3.12,-3.9,-3.35,-4.24,-3.86,-4.92,-5.06,-6.77,-7.41,-9.18,
     *-10.16,-11.12,-9.76,-9.25,-7.96,-7.65/
      DATA  T2XY/
     *-8.4,-12.9,-8.,-11.9,-9.2,-12.5,-10.8,-13.6,-11.2,-12.2,-12.81,
     *-15.4,-13.07,-15.8,-13.81,-14.98,-12.63,-13.76,-11.37,-12.38,
     *-9.23,-9.65,-7.64,-9.17,-8.05,-9.72,-8.87,-10.76,-8.64,-8.89,-6.6,
     *-7.13,-4.77,-5.33,-3.06,-3.79,-1.72,-2.79,-.93,-2.19,-.52,-1.9,
     *-.45,-2.2,-1.22,-3.07,-2.42,-4.37,-3.94,-6.08,-4.49,-4.50,-3.14,
     *-2.93,-1.04,-1.36,.69,.21,2.11,1.33,3.29,2.46,4.3,3.32,4.79,3.62,
     *4.97,3.64,4.63,3.07,4.06,2.49,3.3,1.46,2.06,.51,.74,-1.18,-1.26,
     *-3.54,-3.97,-5.26,-4.18,-3.71,-2.1,-1.7,-.08,-.18,.94,.27,1.13,
     *.08,.91,-.31,.49,-.78,.08,-1.15,-.23,-1.41,-.42,-1.55,-.55,-1.66,
     *-.66,-1.73,-.75,-1.74,-.78,-1.69,-.78,-1.6,-.75,-1.46,-.67,-1.26,
     *-.51,-1.04,-.53,-1.84,-2.42,-4.52,-4.76,-6.33,-6.76,-7.81,-5.8,
     *-5.37,-3.63,-3.35,-1.75,-1.88,-.61,-.9,.09,-.32,.55,-.13,.7,-.06,
     *.49,-.2,.4,-.22,.36,-.09,.58,.12,.75,.15,.7,.17,1.11,.89,1.85,
     *1.62,2.54,2.29,3.2,2.91,3.84,3.53,4.48,4.15,5.12,4.78,5.75,5.39,
     *6.31,5.91,6.87,6.33,7.13,6.61,7.3,6.31,6.27,4.83,4.49,2.85,2.32,
     *.58,-.11,-.98,.81,1.77,3.37,4.13,5.6,6.15,7.29,7.35,7.95,7.67,
     *8.16,7.83,8.31,8.01,8.53,8.27/
C
      END
