! NOT USED (A.Timofeev)
      SUBROUTINE INMACR
C INPUT OF MEDIA DATA
C
      COMMON /INFEL/ ZNUC(110),ATWEI(110),DENS(110),RIEV(110),SYMB(110)
      CHARACTER*2 SYMB
      COMMON /MACRIN/ ELMIX(8,24,48),RHOMED(48),NELEMD(2,48),NUMMED  ! CORNEL
C       NUMMED - Number of media (1-48)
C         NELEMD(1,J) - Type of the medium No.J (1-4)
C         NELEMD(2,J) - Number of chemical elements in the medium No.J (1-24)
C         RHOMED(J) - Density of the medium No.J (g/cm^3)
C           For the element No.K in the medium No.J:
C         ELMIX(1,K,J) - NUCLID (ID of the chemical element, 1-110)
C         ELMIX(2,K,J) - Concentration (10^27/cm^3)
C         ELMIX(3,K,J) - Partial density (g/cm^3)
C         ELMIX(4,K,J) - Atomic number Z
C         ELMIX(5,K,J) - Atomic weight A
C         ELMIX(6,K,J) - Density of pure material (g/cm^3)
C         ELMIX(7,K,J) - Potential of ionization (eV)
C         ELMIX(8,K,J) - Not used
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
C
      AVOG=6.022169E-04
C
      READ(22,1)NUMMED
    1 FORMAT(7X,I2)
          IF(NUMMED.GE.1 .AND. NUMMED.LE.48)GO TO 3  ! CORNEL
      CONTINUE!WRITE(25,2)NUMMED
    2 FORMAT(' INMACR DIAGNOSTIC: INVALID NUMMED=',I2)
          STOP
C
    3 DO 100 JJ=1,NUMMED ! ---------- Loop on all media ---------------
C
      READ(22,4)MEDTYP,NELEM,RHO
    4 FORMAT(/7X,I1,6X,I2,5X,F7.4)  ! CORNEL (Format changed: 6X,I2)
        IF(MEDTYP.GE.1 .AND. MEDTYP.LE.4)GO TO 6
      CONTINUE!WRITE(25,5)MEDTYP,JJ
    5 FORMAT(' INMACR DIAGNOSTIC: INVALID MEDTYP=',I2,'  FOR MEDIUM',I2)
        STOP
C
    6   IF(NELEM.GE.0 .AND. NELEM.LE.24)GO TO 8  ! CORNEL
      CONTINUE!WRITE(25,7)NELEM,JJ
    7 FORMAT(' INMACR DIAGNOSTIC: INVALID NELEM=',I2,'  FOR MEDIUM ',I2)
        STOP
C
    8 RHOMED(JJ)=RHO  ! Can be zero for MEDTYP=1 but then will be defined
        IF(MEDTYP.EQ.1)NELEM=1
      NELEMD(1,JJ)=MEDTYP
      NELEMD(2,JJ)=NELEM
C
        DO KK=1,NELEM
          READ(22,51)NUCLID,(ELMIX(K,KK,JJ),K=2,7)
   51     FORMAT(9X,I3,6X,E11.5,6X,F7.4,3X,F3.0,3X,F8.4,5X,F7.4,3X,F6.1)
          ELMIX(1,KK,JJ)=FLOAT(NUCLID)
                                      ELMIX(4,KK,JJ)=ZNUC(NUCLID)
                                      ELMIX(5,KK,JJ)=ATWEI(NUCLID)
             IF(ELMIX(6,KK,JJ).EQ.0.) ELMIX(6,KK,JJ)=DENS(NUCLID)
             IF(ELMIX(7,KK,JJ).EQ.0.) ELMIX(7,KK,JJ)=RIEV(NUCLID)
        END DO
C
      GO TO(11,21,31,41),MEDTYP
C
   11 RHO=ELMIX(6,1,JJ)
      RHOMED(JJ)=RHO
      ELMIX(2,1,JJ)=AVOG*(RHO/ELMIX(5,1,JJ))
      ELMIX(3,1,JJ)=RHO
      GO TO 100
C
   21 RMOL=0.
      DO 22 I=1,NELEM
   22 RMOL=RMOL+ELMIX(2,I,JJ)*ELMIX(5,I,JJ)
          DO 23 I=1,NELEM
          CONC=(ELMIX(2,I,JJ)*RHO*AVOG)/RMOL
          RHOT=(ELMIX(2,I,JJ)*ELMIX(5,I,JJ)*RHO)/RMOL
          ELMIX(2,I,JJ)=CONC
          ELMIX(3,I,JJ)=RHOT
   23     CONTINUE
      GO TO 100
C
   31 AVER=0.
      DO 32 I=1,NELEM
   32 AVER=AVER+ELMIX(2,I,JJ)*ELMIX(5,I,JJ)
          DO 33 I=1,NELEM
          CONC=(ELMIX(2,I,JJ)*RHO*AVOG)/AVER
          RHOT=(ELMIX(2,I,JJ)*ELMIX(5,I,JJ)*RHO)/AVER
          ELMIX(2,I,JJ)=CONC
          ELMIX(3,I,JJ)=RHOT
   33     CONTINUE
      GO TO 100
C
   41 GO TO 100
C
  100 CONTINUE  ! ---------- End of loop on all media ---------------
C
C
      CONTINUE!WRITE(26,201)NUMMED
      CONTINUE!WRITE(27,201)NUMMED
      CONTINUE!WRITE(24,201)NUMMED
      CONTINUE!WRITE(46,201)NUMMED   ! userset
  201 FORMAT(/10X,'MEDIA COMPOSITION.'/
     *        10X,'NUMBER OF MEDIA is ',I2)
C
      DO 300 JJ=1,NUMMED
C
      CONTINUE!WRITE(26,210)JJ,NELEMD(1,JJ),NELEMD(2,JJ),RHOMED(JJ)
      CONTINUE!WRITE(27,210)JJ,NELEMD(1,JJ),NELEMD(2,JJ),RHOMED(JJ)
      CONTINUE!WRITE(24,210)JJ,NELEMD(1,JJ),NELEMD(2,JJ),RHOMED(JJ)
      CONTINUE!WRITE(46,210)JJ,NELEMD(1,JJ),NELEMD(2,JJ),RHOMED(JJ)   ! userset
  210 FORMAT(10X,'MEDIUM No.',I2,'. Medium TYPE is ',I1/
     *24X,'NUMBER OF ELEMENTS in this medium is ',I1/
     *24X,'DENSITY of this medium is ',F7.4,' g/cm**3')
C
      CONTINUE!WRITE(26,211)
      CONTINUE!WRITE(27,211)
      CONTINUE!WRITE(24,211)
      CONTINUE!WRITE(46,211)   ! userset
  211 FORMAT('    NCL  Z',9X,'A     RHO  I(eV)  Concentration',
     *'  Partial Density'/21X,'g/cm**3',10X,'10**27/cm**3',6X,
     *'g/cm**3')
C
      NELEM=NELEMD(2,JJ)
      DO 212 KK=1,NELEM
      NUCLID=IFIX(ELMIX(1,KK,JJ))
      CONTINUE!WRITE(26,213)  NUCLID,        ELMIX(4,KK,JJ),ELMIX(5,KK,JJ),
!     *ELMIX(6,KK,JJ),ELMIX(7,KK,JJ),ELMIX(2,KK,JJ),ELMIX(3,KK,JJ)
      CONTINUE!WRITE(27,213)  NUCLID,        ELMIX(4,KK,JJ),ELMIX(5,KK,JJ),
!     *ELMIX(6,KK,JJ),ELMIX(7,KK,JJ),ELMIX(2,KK,JJ),ELMIX(3,KK,JJ)
      CONTINUE!WRITE(24,213)  NUCLID,        ELMIX(4,KK,JJ),ELMIX(5,KK,JJ),
!     *ELMIX(6,KK,JJ),ELMIX(7,KK,JJ),ELMIX(2,KK,JJ),ELMIX(3,KK,JJ)
      CONTINUE!WRITE(46,213)  NUCLID,        ELMIX(4,KK,JJ),ELMIX(5,KK,JJ),   ! userset
!     *ELMIX(6,KK,JJ),ELMIX(7,KK,JJ),ELMIX(2,KK,JJ),ELMIX(3,KK,JJ)    ! userset
  212 CONTINUE
  213 FORMAT(I7,F4.0,F9.4,F8.4,F7.1,4X,E11.5,6X,F7.4)
  300 CONTINUE
C
      RETURN
      END



      BLOCK DATA SLAC02
      COMMON /INFEL/ ZNUC(110),ATWEI(110),DENS(110),RIEV(110),SYMB(110)
      CHARACTER*2 SYMB
      COMMON /SIGMAC/ ENET(60),SIGMAC(60,15,48),NUPAR0     ! CORNEL   !HION
      COMMON /RANMAC/ TR(100),RANMAC(100,13,48),OPTMAC(100,13,48),NUPARZ  ! CORNEL !HION
C
      DATA ZNUC/
     *  1.,  2.,  3.,  4.,  5.,  6.,  7.,  8.,  9., 10., 11., 12.,
     * 13., 14., 15., 16., 17., 18., 19., 20., 21., 22., 23., 24.,
     * 25., 26., 27., 28., 29., 30., 31., 32., 33., 34., 35., 36.,
     * 37., 38., 39., 40., 41., 42., 43., 44., 45., 46., 47., 48.,
     * 49., 50., 51., 52., 53., 54., 55., 56., 57., 58., 59., 60.,
     * 61., 62., 63., 64., 65., 66., 67., 68., 69., 70., 71., 72.,
     * 73., 74., 75., 76., 77., 78., 79., 80., 81., 82., 83., 84.,
     * 85., 86., 87., 88., 89., 90., 91., 92., 93., 94., 95., 96.,
     * 97., 98., 99.,100.,  1.,  1., 92.,  2.,  3.,  5., 94., 3*1./
C
      DATA ATWEI/
     *   1.00797,   4.00260,   6.93900,   9.01220,  10.81100,  12.01115,
     *  14.00670,  15.99940,  18.99840,  20.18300,  22.98980,  24.31200,
     *  26.98150,  28.08800,  30.97380,  32.06400,  35.45300,  39.94800,
     *  39.10200,  40.08000,  44.95600,  47.90000,  50.94200,  51.99800,
     *  54.93800,  55.84700,  58.93320,  58.71000,  63.54000,  65.37000,
     *  69.72000,  72.59000,  74.92160,  78.96000,  79.80800,  83.80000,
     *  85.47000,  87.62000,  88.90500,  91.22000,  92.90600,  95.94000,
     *  99.00000, 101.07000, 102.90500, 106.40000, 107.87000, 112.40000,
     * 114.82000, 118.69000, 121.75000, 127.60000, 126.90440, 131.30000,
     * 132.90500, 137.34000, 138.91000, 140.12000, 140.90700, 144.24001,
     * 147.00000, 150.35001, 151.98000, 157.25000, 158.92400, 162.50000,
     * 164.92999, 167.25999, 168.93401, 173.03999, 174.97000, 178.49001,
     * 180.94800, 183.85001, 186.20000, 190.20000, 192.20000, 195.08000,
     * 196.98700, 200.59000, 204.37000, 207.19000, 208.98000, 210.00000,
     * 210.00000, 222.00000, 223.00000, 226.00000, 227.00000, 232.03600,
     * 231.00000, 238.03000, 237.00000, 242.00000, 243.00000, 247.00000,
     * 247.00000, 248.00000, 254.00000, 253.00000,   2.00000,   3.00000,
     * 235.00000,   3.00000,   6.00000,  10.81100, 240.00000, 3*1.00000/
C
      DATA DENS/
     *   0.0808,    0.1900,    0.5340,    1.8500,    2.5000,    2.2600,
     *   1.1400,    1.5680,    1.5000,    1.0000,    0.9712,    1.7400,
     *   2.7020,    2.4000,    1.8200,    2.0700,    2.2000,    1.6500,
     *   0.8600,    1.5500,    3.0200,    4.5400,    5.8700,    7.1400,
     *   7.3000,    7.8600,    8.7100,    8.9000,    8.9333,    7.1400,
     *   5.9100,    5.3600,    5.7300,    4.8000,    4.2000,    3.4000,
     *   1.5300,    2.6000,    4.4700,    6.4000,    8.5700,    9.0100,
     *  11.5000,   12.2000,   12.5000,   12.0000,   10.5000,    8.6500,
     *   7.3000,    7.3100,    6.6840,    6.2400,    4.9300,    2.7000,
     *   1.8730,    3.5000,    6.1500,    6.9000,    6.7690,    7.0070,
     *   1.0000,    7.5400,    5.1700,    7.8700,    8.2500,    8.5600,
     *   8.8000,    9.0600,    9.3200,    6.9600,    9.8500,   11.4000,
     *  16.6000,   19.3000,   20.5300,   22.4800,   22.4200,   21.4500,
     *  19.3000,   14.1900,   11.8500,   11.3400,    9.7800,    9.3000,
     *   1.0000,    4.0000,    1.0000,    5.0000,    1.0000,   11.0000,
     *  15.3700,   18.9000,   20.5000,   19.7370,   11.7000,    7.0000,
     *   1.0000,    1.0000,    1.0000,    1.0000,    0.0808,    0.0808,
     *  18.9000,    0.1900,    0.5340,    2.5000,   19.7370,  3*1.0000/
C
      DATA RIEV/
     *  19.2,  41.8,  40.0,  63.7,  76.0,  78.0,  82.0,  95.0, 115.0,
     * 137.0, 149.0, 156.0, 166.0, 173.0, 173.0, 180.0, 174.0, 188.0,
     * 190.0, 191.0, 216.0, 233.0, 245.0, 257.0, 272.0, 286.0, 297.0,
     * 311.0, 322.0, 330.0, 334.0, 350.0, 347.0, 348.0, 357.0, 352.0,
     * 363.0, 366.0, 379.0, 393.0, 417.0, 424.0, 428.0, 441.0, 449.0,
     * 470.0, 470.0, 469.0, 488.0, 488.0, 487.0, 485.0, 491.0, 482.0,
     * 488.0, 491.0, 501.0, 523.0, 535.0, 546.0, 560.0, 574.0, 580.0,
     * 591.0, 614.0, 628.0, 650.0, 658.0, 674.0, 684.0, 694.0, 705.0,
     * 718.0, 727.0, 736.0, 746.0, 757.0, 790.0, 790.0, 800.0, 810.0,
     * 823.0, 823.0, 830.0, 825.0, 794.0, 827.0, 826.0, 841.0, 847.0,
     * 878.0, 890.0, 902.0, 921.0, 934.0, 939.0, 952.0, 966.0, 980.0,
     * 994.0,  19.2,  19.2, 890.0,  41.8,  40.0,  76.0, 921.0, 3*10./
C
      DATA SYMB/
     * ' H', 'He', 'Li', 'Be', ' B', ' C', ' N', ' O', ' F', 'Ne',
     * 'Na', 'Mg', 'Al', 'Si', ' P', ' S', 'Cl', 'Ar', ' K', 'Ca',
     * 'Sc', 'Ti', ' V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
     * 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', ' Y', 'Zr',
     * 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn',
     * 'Sb', 'Te', ' I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd',
     * 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',
     * 'Lu', 'Hf', 'Ta', ' W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
     * 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th',
     * 'Pa', ' U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm',
     * ' D', ' T', ' U', 'He', 'Li', ' B', 'Pu', 3*'ZZ'/
C ------------------ENERGY MESH-----------------------------------------------
      DATA ENET/
     *   0.,   3.,  10.,  14.,  15.,  16.,  17.,  18.,  20.,  22.,
     *  25.,  27.,  30.,  33.,  35.,  40.,  45.,  50.,  55.,  60.,
     *  65.,  70.,  80.,  90., 100., 110., 120., 130., 140., 150.,
     * 160., 180., 200., 220., 240., 250., 260., 280., 300., 350.,
     * 400., 450., 500., 550., 600., 700., 800., 900.,1000.,1500.,
     *2000.,3000.,5000.,7000.,10000.,20000.,50000.,100000.,500000.,
     *1000000./
C
      DATA TR/
     *3.,4.,5.,6.,7.,8.,9.,10.,12.,14.,16.,18.,20.,22.5,25.,27.5,30.,
     *35.,40.,45.,50.,55.,60.,65.,70.,75.,80.,90.,100.,110.,120.,130.,
     *140.,150.,160.,180.,200.,225.,250.,275.,300.,325.,350.,375.,400.,
     *450.,500.,550.,600.,700.,800.,900.,1000.,1250.,1500.,1750.,2000.,
     *2250.,2500.,2750.,3000.,3500.,4000.,4500.,5000.,6000.,7000.,8000.,
     *9000.,10000.,12500.,15000.,17500.,20000.,22500.,25000.,27500.,
     *30000.,40000.,50000.,60000.,70000.,80000.,90000.,100000.,
     *125000.,150000.,175000.,200000.,225000.,250000.,275000.,
     *300000.,400000.,500000.,600000.,700000.,800000.,900000.,1000000./
C
      DATA SIGMAC /43200*0.0/           ! CORNEL                !HION
      DATA RANMAC /62400*0.0/           ! CORNEL                !HION
      DATA OPTMAC /62400*0.0/           ! CORNEL                !HION
C
      END


      SUBROUTINE BDPRIN
      COMMON /INFEL/ ZNUC(110),ATWEI(110),DENS(110),RIEV(110),SYMB(110)
      CHARACTER*2 SYMB
C
      CONTINUE!WRITE(26,1)
    1 FORMAT(/5X,'DEFAULT ATOMIC NUMBERS, SYMBOLS, WEIGHTS, DENSITIES'/
     *        5X,'AND IONIZATION POTENTIALS FOR ELEMENTS IN SHIELD/RD')
      CONTINUE!WRITE(26,2)
    2 FORMAT(/'    NUCLID',9X,'Z    SYMBOL  ATWEIGHT   DENSITY',
     *        '     I(EV)')
      DO 3 J=1,110
    3 CONTINUE!WRITE(26,4) J,ZNUC(J),SYMB(J),ATWEI(J),DENS(J),RIEV(J)
    4 FORMAT(I10,F10.0,6X,A2,2X,F10.5,F10.4,F10.1)
C
      RETURN
      END



      SUBROUTINE CALMAC
C CALCULATION OF MACROSCOPIC CROSS-SECTIONS FOR JPART=1-24.
C
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      COMMON /MACRIN/ ELMIX(8,24,48),RHOMED(48),NELEMD(2,48),NUMMED  ! CORNEL
      COMMON /SIGMAC/ ENET(60),SIGMAC(60,15,48),NUPAR0     ! CORNEL   !HION
C                                     !  !  !___ medium
C                                     !  !______ particle
C                                     !_________ energy grid
C
      COMMON /RECJPR/ LPART(25),INVERP(15),LPARTZ(25),INVERZ(13)  !HION
C
C------------------- MACROSCOPICAL CROSS SECTIONS (1/CM) ------------------
      DO 1 JPART=1,24             !HION:  24 but not 25!  For HI in CALMHI!
        IPART=LPART(JPART)
        IF(IPART.LE.-1)GO TO 1
            DO MEDIA=1,NUMMED
              NELEM=NELEMD(2,MEDIA)
                  DO J=1,60
                    TMEV=ENET(J)
                      IF(JPART.eq.21)TMEV=TMEV*2.                    !HION8
                      IF(JPART.eq.22)TMEV=TMEV*3.                    !HION8
                      IF(JPART.eq.23)TMEV=TMEV*3.                    !HION8
                      IF(JPART.eq.24)TMEV=TMEV*4.                    !HION8
                        SMAC=0.
                        DO I=1,NELEM
                          NUCLID=IFIX(ELMIX(1,I,MEDIA))
                          CONC  =ELMIX(2,I,MEDIA)
                          CALL MICROD(JPART,TMEV,NUCLID,SIGTOT,SIGIN)
                          SMAC=SMAC+CONC*SIGTOT
                        END DO
                    SIGMAC(J,IPART,MEDIA)=SMAC
                  END DO
            END DO
    1 CONTINUE
C
C--------------------- PRINT MACROSECTIONS ---------------------------
C     PRINTING IS TRUNCATED
C      DO J=1,NUMMED
C        CONTINUE!WRITE(26,10)J
C   10   FORMAT(/5X,'MACROSCOPIC CROSS SECTIONS (1/CM)',
c     *             ' FOR MEDIUM NO.',I2)
C
C        CONTINUE!WRITE(26,11)
C   11   FORMAT(/'  E(MEV) NEUTRN PROTON    PI-    PI+ ANTI N ANTI P',
C     *         5X,'K-     K+',2X)
C          DO K=1,60
C            CONTINUE!WRITE(26,12) IFIX(ENET(K)),(SIGMAC(K,L,J),L=1,8)
C   12       FORMAT(I8,8F7.4)
C          END DO
C
C       ******************* !HION *********************
C        CONTINUE!WRITE(26,13)
C   13   FORMAT(/' E(MEV/A)    K0 ANTIK0  DEUTR TRITON    He3    He4',
C     *        16X)
c          DO K=1,60
C            CONTINUE!WRITE(26,14) IFIX(ENET(K)),(SIGMAC(K,L,J),L=9,14)
C  14       FORMAT(I8,6F7.4)
C          END DO
C       ****************** end HION *******************
C
C      END DO
C
      RETURN
      END



      SUBROUTINE ROPTPR                                              !HION8
C CALCULATION OF RANGE-ENERGY DEPENDENCE AND OPTICAL DEPTH FOR JPART=1-24.
C
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      COMMON /MACRIN/ ELMIX(8,24,48),RHOMED(48),NELEMD(2,48),NUMMED  ! CORNEL
      COMMON /RANMAC/ TR(100),RANMAC(100,13,48),OPTMAC(100,13,48),NUPARZ  ! CORNEL !HION
C                                    !   !  !___ medium
C                                    !   !_____ charged particle
C                                    !_________ energy grid
C
      COMMON /RECJPR/ LPART(25),INVERP(15),LPARTZ(25),INVERZ(13)  !HION
      DIMENSION RTEMP(100),RTEMP1(100)
C
C--------------------- RANGE-ENERGY DEPENDENCE AND OPTICAL DEPTH --------------
            APART=1.                                                 !HION8
            TMIN=2.
      DO 101 JPART=1,24           !HION:  24 but not 25!  For HI in CALMHI!
        IPARTZ=LPARTZ(JPART)
        IF(IPARTZ.LE.-1)GO TO 101
          IF(JPART.eq.21)APART=2.                                    !HION8
          IF(JPART.eq.22)APART=3.                                    !HION8
          IF(JPART.eq.23)APART=3.                                    !HION8
          IF(JPART.eq.24)APART=4.                                    !HION8
            DO MEDIA=1,NUMMED
C
                  DO J=1,100
                      IF(J.EQ.1)THEN
                        T1=TMIN*APART                                !HION8
                        T2=TR(J)*APART                               !HION8
                      ELSE
                        T1=TR(J-1)*APART                             !HION8
                        T2=TR(J)*APART                               !HION8
                      END IF
                    RTEMP(J)=RT1T2(JPART,MEDIA,T1,T2)
                    RTEMP1(J)=OPT1T2(JPART,MEDIA,T1,T2)
                  END DO
C
                  DO J=2,100
                    RTEMP(J)=RTEMP(J-1)+RTEMP(J)
                    RTEMP1(J)=RTEMP1(J-1)+RTEMP1(J)
                  END DO
C
                  DO J=1,100
                    RANMAC(J,IPARTZ,MEDIA)=RTEMP(J)
                    OPTMAC(J,IPARTZ,MEDIA)=RTEMP1(J)
                  END DO
C
            END DO
  101 CONTINUE
C
C-------------- PRINT RANGE-ENERGY DEPENDENCE AND OPTICAL DEPTH ------------
C     PRINT IS TRUNCATED (A.Timofeev)
C      DO J=1,NUMMED
C        CONTINUE!WRITE(26,110)J
C 110   FORMAT(//10X,'RANGE-ENERGY DEPENDENCE (MEV-CM) AND OPTICAL'/
c    *           10X,'DEPTH (UNDIMENSIONED) FOR MEDIUM NO.',I2/)
c       CONTINUE!WRITE(26,120)
c 120   FORMAT(13X,4('-'),' PROTON ',    4('-'),
c    *          4X,3('-'),' PI MINUS ',  3('-'),
c    *          4X,4('-'),' PI PLUS ',   3('-'),
c    *          4X,2('-'),' ANTIPROTON ',2('-'))
c       CONTINUE!WRITE(26,130)
c 130   FORMAT('   E(MEV)',4('      Range OptDepth'))
C
c       DO K=1,100
c         CONTINUE!WRITE(26,140) TR(K),RANMAC(K,1,J),OPTMAC(K,1,J),
c    *                        RANMAC(K,2,J),OPTMAC(K,2,J),
c    *                        RANMAC(K,3,J),OPTMAC(K,3,J),
c    *                        RANMAC(K,4,J),OPTMAC(K,4,J)
c 140     FORMAT(F10.1,2X,2F14.5,2X,2F14.5,2X,2F14.5,2X,2F14.5)
c 140     FORMAT(E9.3,2X,2E9.3,2X,2E9.3,2X,2E9.3,2X,2E9.3)
c       END DO
C
c       CONTINUE!WRITE(26,150)
c 150   FORMAT()
c       CONTINUE!WRITE(26,160)
c 160   FORMAT(13X,4('-'),' K MINUS ', 3('-'),
c    *          4X,4('-'),' K PLUS ',  4('-'),
c    *          4X,3('-'),' MU MINUS ',3('-'),
c    *          4X,4('-'),' MU PLUS ', 3('-'))
c       CONTINUE!WRITE(26,130)
C
c       DO K=1,100
c         CONTINUE!WRITE(26,140) TR(K),RANMAC(K,5,J),OPTMAC(K,5,J),
c    *                        RANMAC(K,6,J),OPTMAC(K,6,J),
c    *                        RANMAC(K,7,J),OPTMAC(K,7,J),
c    *                        RANMAC(K,8,J),OPTMAC(K,8,J)
c       END DO
C
C       ******************* !HION *********************
c       CONTINUE!WRITE(26,150)
c       CONTINUE!WRITE(26,170)
c 170   FORMAT(13X,4('-'),' DEUTRON ',3('-'),
c    *          4X,4('-'),' TRITON ', 4('-'),
c    *          4X,4('-'),' HELIUM3 ',3('-'),
c    *          4X,5('-'),' ALPHA ',  4('-'))
c       CONTINUE!WRITE(26,131)
c 131   FORMAT(' E(MEV/A)',4('      Range OptDepth'))
C
c       DO K=1,100
c         CONTINUE!WRITE(26,140) TR(K),RANMAC(K, 9,J),OPTMAC(K, 9,J),
c    *                        RANMAC(K,10,J),OPTMAC(K,10,J),
c    *                        RANMAC(K,11,J),OPTMAC(K,11,J),
c    *                        RANMAC(K,12,J),OPTMAC(K,12,J)
c       END DO
C       ****************** end HION *******************
C
c     END DO
C
      RETURN
      END



      SUBROUTINE CALMHI(IPRIT)                                       !HION
C CALCULATION OF MACROSCOPIC CROSS-SECTIONS FOR JPART=25 (Heavy Ion).
C
C INPUT: IPRIT=1 - printout into file 26
C        IPRIT=0 - no printout
C
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      COMMON /MACRIN/ ELMIX(8,24,48),RHOMED(48),NELEMD(2,48),NUMMED  ! CORNEL
      COMMON /INFEL/ ZNUC(110),ATWEI(110),DENS(110),RIEV(110),SYMB(110)
      CHARACTER*2 SYMB
        COMMON /HIPROJ/ APROJ,ZPROJ                                 !HION
      COMMON /SIGMAC/ ENET(60),SIGMAC(60,15,48),NUPAR0     ! CORNEL   !HION
C                                     !  !  !___ medium
C                                     !  !______ particle (15 - Ion)
C                                     !_________ energy grid
C
      COMMON /RECJPR/ LPART(25),INVERP(15),LPARTZ(25),INVERZ(13)  !HION
C
C------------------- MACROSCOPICAL CROSS SECTIONS (1/CM) ------------------
      JPART=25                                            !HION
        IPART=LPART(JPART)
            DO MEDIA=1,NUMMED
              NELEM=NELEMD(2,MEDIA)
                  DO J=1,60
                    TMEV=ENET(J)*APROJ                               !HION8
                        SMAC=0.
                        DO I=1,NELEM
                          NUCLID=IFIX(ELMIX(1,I,MEDIA))
                          CONC  =ELMIX(2,I,MEDIA)
                          CALL MICROD(JPART,TMEV,NUCLID,SIGTOT,SIGIN)
                          if(SIGTOT.lt.0.)SIGTOT=0.                  ! GSITR
                          SMAC=SMAC+CONC*SIGTOT
                        END DO
                    SIGMAC(J,IPART,MEDIA)=SMAC
                  END DO
            END DO
C
C--------------------- PRINT MACROSECTIONS ---------------------------
C
              if(IPRIT.le.0)return
      CONTINUE!WRITE(26,10)SYMB(INT(ZPROJ)),INT(ZPROJ),APROJ
   10 FORMAT(/10X,'Macroscopic cross section (1/cm) for Nucleus '/
     *        10X,A2,'(',I3,',',F5.1,') in Medium No.:')
C
C---------------------------- CORNEL -------------------------
C Printout for MEDIA Nos. 1-12
      CONTINUE!WRITE(26,11)
   11 FORMAT(' E(MeV/A)','  Med1   Med2   Med3   Med4   Med5   Med6',
     *                  '   Med7   Med8   Med9  Med10  Med11  Med12')
      DO K=1,60
        CONTINUE!WRITE(26,12)INT(ENET(K)),(SIGMAC(K,15,J),J=1,12)
   12   FORMAT(I8,12F7.4)
      END DO
C
C Printout for MEDIA Nos.13-24, if any
                   IF(NUMMED.gt.12)THEN
      CONTINUE!WRITE(26,21)
   21 FORMAT(/' E(MeV/A)',' Med13  Med14  Med15  Med16  Med17  Med18',
     *                   '  Med19  Med20  Med21  Med22  Med23  Med24')
      DO K=1,60
        CONTINUE!WRITE(26,12)INT(ENET(K)),(SIGMAC(K,15,J),J=13,24)
      END DO
                   END IF
C
C Printout for MEDIA Nos.25-36, if any
                   IF(NUMMED.gt.24)THEN
      CONTINUE!WRITE(26,31)
   31 FORMAT(/' E(MeV/A)',' Med25  Med26  Med27  Med28  Med29  Med30',
     *                   '  Med31  Med32  Med33  Med34  Med35  Med36')
      DO K=1,60
        CONTINUE!WRITE(26,12)INT(ENET(K)),(SIGMAC(K,15,J),J=25,36)
      END DO
                   END IF
C
C Printout for MEDIA Nos.37-48, if any
                   IF(NUMMED.gt.36)THEN
      CONTINUE!WRITE(26,41)
   41 FORMAT(/' E(MeV/A)',' Med37  Med38  Med39  Med40  Med41  Med42',
     *                   '  Med43  Med44  Med45  Med46  Med47  Med48')
      DO K=1,60
        CONTINUE!WRITE(26,12)INT(ENET(K)),(SIGMAC(K,15,J),J=37,48)
      END DO
                   END IF
C---------------------------- end CORNEL ---------------------
C
      RETURN
      END



      SUBROUTINE ROPTHI(IPRIT)                                       !HION8
C CALCULATION OF RANGE-ENERGY DEPENDENCE AND OPTICAL DEPTH FOR JPART=25
C                                                            (Heavy Ion).
C INPUT: IPRIT=1 - printout into file 26
C        IPRIT=0 - no printout
C
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      COMMON /MACRIN/ ELMIX(8,24,48),RHOMED(48),NELEMD(2,48),NUMMED  ! CORNEL
      COMMON /INFEL/ ZNUC(110),ATWEI(110),DENS(110),RIEV(110),SYMB(110)
      CHARACTER*2 SYMB
        COMMON /HIPROJ/ APROJ,ZPROJ                                 !HION
      COMMON /RANMAC/ TR(100),RANMAC(100,13,48),OPTMAC(100,13,48),NUPARZ  ! CORNEL !HION
C                                    !   !  !___ medium
C                                    !   !_____ charged particle (13 - Ion)
C                                    !_________ energy grid
C
      COMMON /RECJPR/ LPART(25),INVERP(15),LPARTZ(25),INVERZ(13)  !HION
      DIMENSION RTEMP(100),RTEMP1(100)
C
C--------------------- RANGE-ENERGY DEPENDENCE AND OPTICAL DEPTH --------------
            TMIN=2.
      JPART=25                                            !HION
        IPARTZ=LPARTZ(JPART)
            DO MEDIA=1,NUMMED
C
                  DO J=1,100
                      IF(J.EQ.1)THEN
                        T1=TMIN*APROJ                                !HION8
                        T2=TR(J)*APROJ                               !HION8
                      ELSE
                        T1=TR(J-1)*APROJ                             !HION8
                        T2=TR(J)*APROJ                               !HION8
                      END IF
                    RTEMP(J)=RT1T2(JPART,MEDIA,T1,T2)
                    RTEMP1(J)=OPT1T2(JPART,MEDIA,T1,T2)
                  END DO
C
                  DO J=2,100
                    RTEMP(J)=RTEMP(J-1)+RTEMP(J)
                    RTEMP1(J)=RTEMP1(J-1)+RTEMP1(J)
                  END DO
C
                  DO J=1,100
                    RANMAC(J,IPARTZ,MEDIA)=RTEMP(J)
                    OPTMAC(J,IPARTZ,MEDIA)=RTEMP1(J)
                  END DO
C
            END DO
C
C-------------- PRINT RANGE-ENERGY DEPENDENCE AND OPTICAL DEPTH ------------
C
              if(IPRIT.le.0)return
      CONTINUE!WRITE(26,111)SYMB(INT(ZPROJ)),INT(ZPROJ),APROJ
  111 FORMAT(/10X,'Range-Energy dependence (MeV-cm) and optical '/
     *        10X,'depth (dimensionless) for Nucleus ',A2,'(',I3,',',
     *        F5.1,') in Medium No.:')
C
C---------------- Printout of MEDIA 1-4 --------------------------
      CONTINUE!WRITE(26,121)
  121 FORMAT(13X,4('-'),' MED  1 ',4('-'),
     *        4X,4('-'),' MED  2 ',4('-'),
     *        4X,4('-'),' MED  3 ',4('-'),
     *        4X,4('-'),' MED  4 ',4('-'))
      CONTINUE!WRITE(26,130)
  130 FORMAT(' E(MEV/A)',4('      Range OptDepth'))
      DO K=1,100
        CONTINUE!WRITE(26,140) TR(K),RANMAC(K,13,1),OPTMAC(K,13,1),
!     *                      RANMAC(K,13,2),OPTMAC(K,13,2),
!     *                      RANMAC(K,13,3),OPTMAC(K,13,3),
!     *                      RANMAC(K,13,4),OPTMAC(K,13,4)
  140   FORMAT(E9.3,2X,2E9.3,2X,2E9.3,2X,2E9.3,2X,2E9.3)
      END DO
      CONTINUE!WRITE(26,150)
  150 FORMAT()
C
      IF(NUMMED.gt. 4)THEN
C       -------------- Printout of MEDIA 5-8 -----------------
        CONTINUE!WRITE(26,122)
  122   FORMAT(13X,4('-'),' MED  5 ',4('-'),
     *          4X,4('-'),' MED  6 ',4('-'),
     *          4X,4('-'),' MED  7 ',4('-'),
     *          4X,4('-'),' MED  8 ',4('-'))
        CONTINUE!WRITE(26,130)
        DO K=1,100
          CONTINUE!WRITE(26,140) TR(K),RANMAC(K,13,5),OPTMAC(K,13,5),
!     *                        RANMAC(K,13,6),OPTMAC(K,13,6),
!     *                        RANMAC(K,13,7),OPTMAC(K,13,7),
!     *                        RANMAC(K,13,8),OPTMAC(K,13,8)
        END DO
        CONTINUE!WRITE(26,150)
      END IF
C
      IF(NUMMED.gt. 8)THEN
C       -------------- Printout of MEDIA 9-12 -----------------
        CONTINUE!WRITE(26,123)
  123   FORMAT(13X,10('-'),' MED  9 ',4('-'),
     *          4X,10('-'),' MED 10 ',4('-'),
     *          4X,10('-'),' MED 11 ',4('-'),
     *          4X,10('-'),' MED 12 ',4('-'))
        CONTINUE!WRITE(26,130)
        DO K=1,100
          CONTINUE!WRITE(26,140) TR(K),RANMAC(K,13, 9),OPTMAC(K,13, 9),
!     *                        RANMAC(K,13,10),OPTMAC(K,13,10),
!     *                        RANMAC(K,13,11),OPTMAC(K,13,11),
!     *                        RANMAC(K,13,12),OPTMAC(K,13,12)
        END DO
        CONTINUE!WRITE(26,150)
      END IF
C
      IF(NUMMED.gt.12)THEN
C       -------------- Printout of MEDIA 13-16 -----------------
        CONTINUE!WRITE(26,124)
  124   FORMAT(13X,4('-'),' MED 13 ',4('-'),
     *          4X,4('-'),' MED 14 ',4('-'),
     *          4X,4('-'),' MED 15 ',4('-'),
     *          4X,4('-'),' MED 16 ',4('-'))
        CONTINUE!WRITE(26,130)
        DO K=1,100
          CONTINUE!WRITE(26,140) TR(K),RANMAC(K,13,13),OPTMAC(K,13,13),
!     *                        RANMAC(K,13,14),OPTMAC(K,13,14),
!     *                        RANMAC(K,13,15),OPTMAC(K,13,15),
!     *                        RANMAC(K,13,16),OPTMAC(K,13,16)
        END DO
        CONTINUE!WRITE(26,150)
      END IF
C
C ------------------------ CORNEL ------------------------------
      IF(NUMMED.gt.16)THEN
C       -------------- Printout of MEDIA 17-20 -----------------
        CONTINUE!WRITE(26,125)
  125   FORMAT(13X,4('-'),' MED 17 ',4('-'),
     *          4X,4('-'),' MED 18 ',4('-'),
     *          4X,4('-'),' MED 19 ',4('-'),
     *          4X,4('-'),' MED 20 ',4('-'))
        CONTINUE!WRITE(26,130)
        DO K=1,100
          CONTINUE!WRITE(26,140) TR(K),RANMAC(K,13,17),OPTMAC(K,13,17),
!     *                        RANMAC(K,13,18),OPTMAC(K,13,18),
!     *                        RANMAC(K,13,19),OPTMAC(K,13,19),
!     *                        RANMAC(K,13,20),OPTMAC(K,13,20)
        END DO
        CONTINUE!WRITE(26,150)
      END IF
C
      IF(NUMMED.gt.20)THEN
C       -------------- Printout of MEDIA 21-24 -----------------
        CONTINUE!WRITE(26,126)
  126   FORMAT(13X,4('-'),' MED 21 ',4('-'),
     *          4X,4('-'),' MED 22 ',4('-'),
     *          4X,4('-'),' MED 23 ',4('-'),
     *          4X,4('-'),' MED 24 ',4('-'))
        CONTINUE!WRITE(26,130)
        DO K=1,100
          CONTINUE!WRITE(26,140) TR(K),RANMAC(K,13,21),OPTMAC(K,13,21),
!     *                        RANMAC(K,13,22),OPTMAC(K,13,22),
!     *                        RANMAC(K,13,23),OPTMAC(K,13,23),
!     *                        RANMAC(K,13,24),OPTMAC(K,13,24)
        END DO
        CONTINUE!WRITE(26,150)
      END IF
C
      IF(NUMMED.gt.24)THEN
C       -------------- Printout of MEDIA 25-28 -----------------
        CONTINUE!WRITE(26,127)
  127   FORMAT(13X,4('-'),' MED 25 ',4('-'),
     *          4X,4('-'),' MED 26 ',4('-'),
     *          4X,4('-'),' MED 27 ',4('-'),
     *          4X,4('-'),' MED 28 ',4('-'))
        CONTINUE!WRITE(26,130)
        DO K=1,100
          CONTINUE!WRITE(26,140) TR(K),RANMAC(K,13,25),OPTMAC(K,13,25),
!     *                        RANMAC(K,13,26),OPTMAC(K,13,26),
!     *                        RANMAC(K,13,27),OPTMAC(K,13,27),
!     *                        RANMAC(K,13,28),OPTMAC(K,13,28)
        END DO
        CONTINUE!WRITE(26,150)
      END IF
C
      IF(NUMMED.gt.28)THEN
C       -------------- Printout of MEDIA 29-32 -----------------
        CONTINUE!WRITE(26,128)
  128   FORMAT(13X,4('-'),' MED 29 ',4('-'),
     *          4X,4('-'),' MED 30 ',4('-'),
     *          4X,4('-'),' MED 31 ',4('-'),
     *          4X,4('-'),' MED 32 ',4('-'))
        CONTINUE!WRITE(26,130)
        DO K=1,100
          CONTINUE!WRITE(26,140) TR(K),RANMAC(K,13,29),OPTMAC(K,13,29),
!     *                        RANMAC(K,13,30),OPTMAC(K,13,30),
!     *                        RANMAC(K,13,31),OPTMAC(K,13,31),
!     *                        RANMAC(K,13,32),OPTMAC(K,13,32)
        END DO
        CONTINUE!WRITE(26,150)
      END IF
C
      IF(NUMMED.gt.32)THEN
C       -------------- Printout of MEDIA 33-36 -----------------
        CONTINUE!WRITE(26,129)
  129   FORMAT(13X,4('-'),' MED 33 ',4('-'),
     *          4X,4('-'),' MED 34 ',4('-'),
     *          4X,4('-'),' MED 35 ',4('-'),
     *          4X,4('-'),' MED 36 ',4('-'))
        CONTINUE!WRITE(26,130)
        DO K=1,100
          CONTINUE!WRITE(26,140) TR(K),RANMAC(K,13,33),OPTMAC(K,13,33),
!     *                        RANMAC(K,13,34),OPTMAC(K,13,34),
!     *                        RANMAC(K,13,35),OPTMAC(K,13,35),
!     *                        RANMAC(K,13,36),OPTMAC(K,13,36)
        END DO
        CONTINUE!WRITE(26,150)
      END IF
C
      IF(NUMMED.gt.36)THEN
C       -------------- Printout of MEDIA 37-40 -----------------
        CONTINUE!WRITE(26,1291)
 1291   FORMAT(13X,4('-'),' MED 37 ',4('-'),
     *          4X,4('-'),' MED 38 ',4('-'),
     *          4X,4('-'),' MED 39 ',4('-'),
     *          4X,4('-'),' MED 40 ',4('-'))
        CONTINUE!WRITE(26,130)
        DO K=1,100
          CONTINUE!WRITE(26,140) TR(K),RANMAC(K,13,37),OPTMAC(K,13,37),
!     *                        RANMAC(K,13,38),OPTMAC(K,13,38),
!     *                        RANMAC(K,13,39),OPTMAC(K,13,39),
!     *                        RANMAC(K,13,40),OPTMAC(K,13,40)
        END DO
        CONTINUE!WRITE(26,150)
      END IF
C
      IF(NUMMED.gt.40)THEN
C       -------------- Printout of MEDIA 41-44 -----------------
        CONTINUE!WRITE(26,1292)
 1292   FORMAT(13X,4('-'),' MED 41 ',4('-'),
     *          4X,4('-'),' MED 42 ',4('-'),
     *          4X,4('-'),' MED 43 ',4('-'),
     *          4X,4('-'),' MED 44 ',4('-'))
        CONTINUE!WRITE(26,130)
        DO K=1,100
          CONTINUE!WRITE(26,140) TR(K),RANMAC(K,13,41),OPTMAC(K,13,41),
!     *                        RANMAC(K,13,42),OPTMAC(K,13,42),
!     *                        RANMAC(K,13,43),OPTMAC(K,13,43),
!     *                        RANMAC(K,13,44),OPTMAC(K,13,44)
        END DO
        CONTINUE!WRITE(26,150)
      END IF
C
      IF(NUMMED.gt.44)THEN
C       -------------- Printout of MEDIA 45-48 -----------------
        CONTINUE!WRITE(26,1293)
 1293   FORMAT(13X,4('-'),' MED 45 ',4('-'),
     *          4X,4('-'),' MED 46 ',4('-'),
     *          4X,4('-'),' MED 47 ',4('-'),
     *          4X,4('-'),' MED 48 ',4('-'))
        CONTINUE!WRITE(26,130)
        DO K=1,100
          CONTINUE!WRITE(26,140) TR(K),RANMAC(K,13,45),OPTMAC(K,13,45),
!     *                        RANMAC(K,13,46),OPTMAC(K,13,46),
!     *                        RANMAC(K,13,47),OPTMAC(K,13,47),
!     *                        RANMAC(K,13,48),OPTMAC(K,13,48)
        END DO
        CONTINUE!WRITE(26,150)
      END IF
C ------------------------ end CORNEL --------------------------
C
      RETURN
      END



      FUNCTION DEDX(JPART,TMEV,Z0,A0,RHO0,RI)
C CALCULATION OF DE/DX (MEV*CM**2/G) ON BETHE-BLOCH FORMULAE
C (SEE PARTICLE PROPERTIES) FOR PARTICLE TYPE  JPART  WITH ENERGY  TMEV
C IN PURE MEDIUM (ELEMENT NUMBER IS  NUCLID)
C
      COMMON /RECJPR/ LPART(25),INVERP(15),LPARTZ(25),INVERZ(13)  !HION
      COMMON /HIPROJ/ APROJ,ZPROJ                                 !HION
C
      ZPART2=1.   ! Squared charge of particle                    !HION
          IF(LPARTZ(JPART).LE.-1)GO TO 50
      IF(JPART.EQ.2 .OR. JPART.EQ.7)  RM0=940.
      IF(JPART.EQ.3 .OR. JPART.EQ.4)  RM0=140.
      IF(JPART.EQ.8 .OR. JPART.EQ.9)  RM0=495.
      IF(JPART.EQ.15 .OR. JPART.EQ.16)RM0=105.
C
C***************** !HION ********************
      IF(JPART.EQ.21)RM0=1880.
      IF(JPART.EQ.22)RM0=2820.
C      IF(JPART.EQ.23)RM0=2820.        ! Correction of 03.09.2000
        IF(JPART.EQ.23)THEN            ! Correction of 03.09.2000
          RM0=2820.
          ZPART2=4.
        END IF
        IF(JPART.EQ.24)THEN
          RM0=3760.
          ZPART2=4.
        END IF
        IF(JPART.EQ.25)THEN
          RM0=940.*APROJ
          ZPART2=ZPROJ**2
        END IF
C***************** end HION *****************
C
      T0=TMEV
C
C     KINEMATIC CHARACTERISTIC
      P0=SQRT(T0*(T0+2.*RM0))
      E0=T0+RM0
      G2=(E0/RM0)**2
      B2=(G2-1.)/G2
      G2B2=G2*B2
      WMAX=P0**2/(RM0*(RM0/1.022+0.511/(2.*RM0)+E0/RM0))
C
C     DENSITY EFFECT CORRECTION DELTA
      X=ALOG10(P0/RM0)
C
      HNP=28.8*SQRT(RHO0*Z0/A0)
      C=-2.*ALOG(RI/HNP)-1.
      CAV=-C
      XAL=CAV/4.606
C
      IF(RI-100.)1,11,11
    1 X1=2.0
      IF(CAV-3.681)2,3,3
    2 X0=0.2
      GO TO 21
    3 X0=0.326*CAV-1.
      GO TO 21
   11 X1=3.0
      IF(CAV-5.215)12,13,13
   12 X0=0.2
      GO TO 21
   13 X0=0.326*CAV-1.5
      GO TO 21
   21 ALIT=4.606*(XAL-X0)/(X1-X0)**3
C
      IF(X-X0)31,31,32
   31 DELTA=0.
      GO TO 40
   32 IF(X-X1)33,33,34
   33 DELTA=4.606*X+C+ALIT*(X1-X)**3
      GO TO 40
   34 DELTA=4.606*X+C
      GO TO 40
   40 CONTINUE
C
      DEDX=0.3070*(Z0/A0)*(ZPART2/B2)*                            !HION
     *      (ALOG((1.022E+06/RI)*G2B2)-B2-DELTA/2.)
      RETURN
C
   50 CONTINUE!WRITE(25,51) JPART
   51 FORMAT(' DEDX DIAGNOSTIC: INVALID JPART=',I2)
      STOP
      END


      FUNCTION DEDXMD(JPART,TMEV,MEDIA)
C CALCULATION OF DE/DX (MEV/CM) ON BETHE-BLOCH FORMULAE FOR PARTICLE
C TYPE  JPART  WITH ENERGY  TMEV  IN COMPLEX MEDIUM NUMBER  MEDIA.
C
      COMMON /MACRIN/ ELMIX(8,24,48),RHOMED(48),NELEMD(2,48),NUMMED  ! CORNEL
C
      NELEM=NELEMD(2,MEDIA)
            DEDXMD=0.
      DO 1 J=1,NELEM
          Z0  =ELMIX(4,J,MEDIA)
          A0  =ELMIX(5,J,MEDIA)
          RHO0=ELMIX(6,J,MEDIA)
          RI  =ELMIX(7,J,MEDIA)
      DEDXMD=DEDXMD+ELMIX(3,J,MEDIA)*DEDX(JPART,TMEV,Z0,A0,RHO0,RI)
    1 CONTINUE
      RETURN
      END


      FUNCTION RT1T2(JPART,MEDIA,T1,T2)
C INTEGRATION FROM T1 TO T2 (T1<T2) OF 1/(DE/DX).
C RT1T2 (CM) IS THE PASS IN COMPLEX MEDIUM NUMBER MEDIA WHEN THE ENERGY
C OF PARTICLE TYPE  JPART  DECREASE FROM T2 TO T1.
C
C      DIMENSION W(8),X(8),X1(8)
C      DATA X/
C     *-0.96028 98565,-0.79666 64774,-0.52553 24099,-0.1834346425,
C     * 0.18343 46425, 0.52553 24099, 0.79666 64774, 0.96028 98565/
C      DATA W/
C     * 0.10122 85363, 0.22238 10345, 0.31370 66459, 0.36268 37834,
C     * 0.36268 37834, 0.31370 66459, 0.22238 10345, 0.10122 85363/
C
      DIMENSION W(16),X(16),X1(16)
      DATA X/
     *-0.98940 09350,-0.94457 50231,-0.86563 12024,-0.75540 44084,
     *-0.61787 62444,-0.45801 67777,-0.28160 35508,-0.09501 25098,
     * 0.09501 25098, 0.28160 35508, 0.45801 67777, 0.61787 62444,
     * 0.75540 44084, 0.86563 12024, 0.94457 50231, 0.98940 09350/
      DATA W/
     * 0.02715 24594, 0.06225 35239, 0.09515 85117, 0.12462 89713,
     * 0.14959 59888, 0.16915 65194, 0.18260 34150, 0.18945 06105,
     * 0.18945 06105, 0.18260 34150, 0.16915 65194, 0.14959 59888,
     * 0.12462 89713, 0.09515 85117, 0.06225 35239, 0.02715 24594/
C
      DO 1 J=1,16
    1 X1(J)=((T2+T1)+(T2-T1)*X(J))*0.5
          SUM=0.
      DO 2 J=1,16
    2 SUM=SUM+W(J)/DEDXMD(JPART,X1(J),MEDIA)
          RT1T2=(T2-T1)*0.5*SUM
      RETURN
      END


      FUNCTION OPT1T2(JPART,MEDIA,T1,T2)
C INTEGRATION FROM T1 TO T2 (T1<T2) OF SECMAC(E)/(DE/DX).
C OPT1T2  IS THE OPTICALDEPTH IN COMPLEX MEDIUM NUMBER MEDIA WHEN THE ENERGY
C OF PARTICLE TYPE  JPART  DECREASE FROM T2 TO T1.
C
      COMMON /RECJPR/ LPART(25),INVERP(15),LPARTZ(25),INVERZ(13)  !HION
C
C      DIMENSION W(8),X(8),X1(8)
C      DATA X/
C     *-0.96028 98565,-0.79666 64774,-0.52553 24099,-0.1834346425,
C     * 0.18343 46425, 0.52553 24099, 0.79666 64774, 0.96028 98565/
C      DATA W/
C     * 0.10122 85363, 0.22238 10345, 0.31370 66459, 0.36268 37834,
C     * 0.36268 37834, 0.31370 66459, 0.22238 10345, 0.10122 85363/
C
      DIMENSION W(16),X(16),X1(16)
      DATA X/
     *-0.98940 09350,-0.94457 50231,-0.86563 12024,-0.75540 44084,
     *-0.61787 62444,-0.45801 67777,-0.28160 35508,-0.09501 25098,
     * 0.09501 25098, 0.28160 35508, 0.45801 67777, 0.61787 62444,
     * 0.75540 44084, 0.86563 12024, 0.94457 50231, 0.98940 09350/
      DATA W/
     * 0.02715 24594, 0.06225 35239, 0.09515 85117, 0.12462 89713,
     * 0.14959 59888, 0.16915 65194, 0.18260 34150, 0.18945 06105,
     * 0.18945 06105, 0.18260 34150, 0.16915 65194, 0.14959 59888,
     * 0.12462 89713, 0.09515 85117, 0.06225 35239, 0.02715 24594/
C
      IF(LPART(JPART).LE.-1)THEN
          OPT1T2=0.
          RETURN
      END IF
C
      DO 1 J=1,16
    1 X1(J)=((T2+T1)+(T2-T1)*X(J))*0.5
          SUM=0.
      DO 2 J=1,16
      CALL MACROS(JPART,X1(J),MEDIA,SECMAC)
      SUM=SUM+W(J)*(SECMAC/DEDXMD(JPART,X1(J),MEDIA))
    2 CONTINUE
          OPT1T2=(T2-T1)*0.5*SUM
      RETURN
      END



      SUBROUTINE MACROS(JPART,TMEV,MEDIA,SECMAC)
C CALCULATION OF MACROSCOPIC CROSS SECTION.
C     INPUT:
C           JPART - PARTICLE TYPE. JPART=1-24. (25?, HION)
C           TMEV  - PARTICLE ENERGY (MEV).
C           MEDIA - NUMBER OF MEDIUM (1-16).
C     OUTPUT:
C           SECMAC - MACROSCOPICAL CROSS SECTION (1/CM).
C
      COMMON /RECJPR/ LPART(25),INVERP(15),LPARTZ(25),INVERZ(13)    !HION
      COMMON /SIGMAC/ ENET(60),SIGMAC(60,15,48),NUPAR0     ! CORNEL   !HION
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
        COMMON /HIPROJ/ APROJ,ZPROJ            ! SHLDRD6           !HION
C
      IPART=LPART(JPART)
      IF(IPART.LE.-1)GO TO 1
C
      IF(JPART.le.20)TMEVA=TMEV                ! SHLDRD6
      IF(JPART.eq.21)TMEVA=TMEV/2.0            ! SHLDRD6
      IF(JPART.eq.22)TMEVA=TMEV/3.0            ! SHLDRD6
      IF(JPART.eq.23)TMEVA=TMEV/3.0            ! SHLDRD6
      IF(JPART.eq.24)TMEVA=TMEV/4.0            ! SHLDRD6
      IF(JPART.eq.25)TMEVA=TMEV/APROJ          ! SHLDRD6
C
      CALL SEARCH(TMEVA,ENET,60,N1,N2)         ! SHLDRD6
      X1=ENET(N1)
      X2=ENET(N2)
      Y1=SIGMAC(N1,IPART,MEDIA)
      Y2=SIGMAC(N2,IPART,MEDIA)
      SECMAC=Y1+((Y2-Y1)/(X2-X1))*(TMEVA-X1)    ! SHLDRD6
      RETURN
C
    1 CONTINUE!WRITE(25,2) JPART
    2 FORMAT(' MACROS DIAGNOSTIC: INVALID JPART=',I2)
      STOP
      END



C      SUBROUTINE SEARCH(E,ARE,N,N1,N2) -- TRUNCATED (A. Timofeev)
C     Already implemented in libHadgen.so
C     BINARY SEARCH IN ARRANGED (IN INCREASING ORDER) ARRAY ARE(N)
C     INPUT:
C          E - CURRENT VARIABLE
C          ARE(N) - VARIABLE MESH
C     OUTPUT:
C          N1,N2  SUCH, THAT  ARE(N1) < E < ARE(N2)
C     WARNING!!!
C               1. ARE(1).LE.E .AND. E.LE.ARE(N)
C               2. N.GT.2
C-------------------------------------------------------------------

      SUBROUTINE SERCH3(E,ARE,N,M,L,IPART,MEDIA,N1,N2)
C BINARY SEARCH IN ARRANGED (IN INCREASING ORDER) ROW OF 3-DIMENSIONAL
C ARRAY ARE(N,M,L)
C     INPUT:
C          E - CURRENT VARIABLE
C          ARE(N,M,L) - VARIABLE MESH
C                M,L - FIXED INDEXES. SEARCH ONLY ON INDEX N.
C     OUTPUT:
C          N1,N2  SUCH, THAT  ARE(N1,M,L) < E < ARE(N2,M,L)
C-------------------------------------------------------------------
C
      DIMENSION ARE(N,M,L)
C
      N1=1
      N2=N
    1 NM=(N1+N2)/2
      IF(E.LE.ARE(NM,IPART,MEDIA))THEN
          N2=NM
      ELSE
          N1=NM
      END IF
      IF(N2.EQ.N1+1)RETURN
      GO TO 1
      END



C      SUBROUTINE INITIS -- TRUNCATED (A. Timofeev)
C     Already implemented in libHadgen.so
C     Input files are also TRUNCATED (data-in-code)


C      FUNCTION AMIXIS(NUCLID)
C Sampling of nucleus Atomic mass number from natural MIXture of ISotopes
C according to isotope abundance.
C   NUCLID - Atomic number (input)
C   AMIXIS - Atomic mass number (output)
C
C     TRUNCATED (A. Timofeev)
C     Already implemented in libHadgen.so

