


      BLOCK DATA GRIDAGT
      COMMON /EGRIDA/ EGRIDA(119)
C Booked fragments  (Zi,Ai)x22, Zi>2
      COMMON /BOOKFR/ BOOKFR(2,22)
C
C Default energy grid for registration secondary particle spectra 
C in h-A and A-A interactions.
      DATA EGRIDA /
     * .000E+00, .100E-06,                               ! new 1E-7
     * .215E-06, .464E-06, .100E-05,
     * .215E-05, .464E-05, .100E-04,
     * .215E-04, .464E-04, .100E-03,
     * .215E-03, .464E-03, .100E-02,
     * .215E-02, .464E-02, .100E-01,
     * .215E-01, .464E-01, .100E+00,
     * .215E+00, .464E+00, .100E+01,                         ! 1 MeV
C
     * .147E+01, .215E+01, .316E+01, .464E+01, .562E+01,
     * .681E+01, .825E+01, .100E+02,                     ! 10 MeV
C
     * .121E+02, .145E+02, .178E+02, .215E+02, .261E+02, ! 14.5 instead of 14.7
     * .316E+02, .383E+02, .464E+02, .511E+02, .562E+02,
     * .619E+02, .681E+02, .750E+02, .825E+02, .909E+02, 
     * .100E+03,                                         ! 100 MeV
C
     * .110E+03, .121E+03, .133E+03, .147E+03, 
     * .162E+03, .178E+03, .196E+03, .215E+03, .237E+03, 
     * .261E+03, .287E+03, .316E+03, .348E+03, .383E+03, 
     * .422E+03, .464E+03, .511E+03, .562E+03, .619E+03, 
     * .681E+03, .750E+03, .825E+03, .909E+03, .100E+04, ! 1 GeV
C
     * .110E+04, .121E+04, .133E+04, .147E+04, 
     * .162E+04, .178E+04, .196E+04, .215E+04, .237E+04, 
     * .261E+04, .287E+04, .316E+04, .348E+04, .383E+04, 
     * .422E+04, .464E+04, .511E+04, .562E+04, .619E+04, 
     * .681E+04, .750E+04, .825E+04, .909E+04, .100E+05, ! 10 GeV
C
     * .121E+05, .147E+05, .178E+05, .215E+05, .261E+05,
     * .316E+05, .383E+05, .464E+05, .511E+05, .562E+05,
     * .619E+05, .681E+05, .750E+05, .825E+05, .909E+05, 
     * .100E+06,                                         ! 100 GeV
C
     * .147E+06, .215E+06, .316E+06, .464E+06, .562E+06,
     * .681E+06, .825E+06, .100E+07/                     ! 1 TeV
C
      DATA BOOKFR /
C (Zi,Ai), Zi>2.  Ai=1000 means all A; Zi=1000 means all Z
C User is invited to modify this list on his discretion. 
C Alternative way is the modification in dialog mode.
c     *  6., 11.,
c     *  7., 14.,
c     *  8., 16.,
c     * 10., 20.,
c     * 13., 27.,
c     * 14., 28.,
c     * 18., 40.,
c     * 20., 40.,
c     * 22., 48.,
c     * 26., 56.,
c     * 29., 63.,
c     * 36., 84.,
c     * 54., 131.,
c     *  6., 1000.,
c     *  7., 1000.,
c     *  8., 1000.,
c     * 13., 1000.,
c     *  1000., 11.,
c     *  1000., 14.,
c     *  1000., 16.,
c     *  1000., 20.,
c     *  1000., 27./
C
cC For reaction C11+O16
c     *  6., 1000.,
c     *  6.,  9.,
c     *  6., 10.,
c     *  6., 11.,
c     *  6., 12.,
c     *  6., 13.,
c     *  6., 14.,
c     *  5., 1000.,
c     *  5.,  9.,
c     *  5., 10.,
c     *  5., 11.,
c     *  5., 12.,
c     *  4., 1000.,
c     *  4.,  7.,
c     *  4.,  8.,
c     *  4.,  9.,
c     *  4., 10.,
c     *  3., 1000.,
c     *  3.,  5.,
c     *  3.,  6.,
c     *  3.,  7.,
c     *  3.,  8./
C
C For reaction C11+O16
     * 79., 1000.,
     * 72., 1000.,
     * 68., 1000.,
     * 60., 1000.,
     * 55., 1000.,
     * 50., 1000.,
     * 45., 1000.,
     * 40., 1000.,
     * 35., 1000.,
     * 30., 1000.,
     * 25., 1000.,
     * 20., 1000.,
     * 15., 1000.,
     * 13., 1000.,
     * 10., 1000.,
     *  9., 1000.,
     *  8., 1000.,
     *  7., 1000.,
     *  6., 1000.,
     *  5., 1000.,
     *  4., 1000.,
     *  3., 1000./
C
      END



      SUBROUTINE HIST3AGT(E,SPEC,NP,NA,W)
C Histogramming of REAL*4 variable E>0 (energy, MeV, MeV/u)
C with statistical weight W according to fixed logarithmic
C grid EGRIDA(119) in the limits 0<E<1E+6
C in 3-dimensional array SPEC(126,30,180)
C INPUT: E - energy, MeV, MeV/u
C        SPEC - name of 3-dimensional array
C          (126 is equal to 118 bins + 8 comulatives)
C        NP - No. of secondary particle/fragment (1-30)
C        NA - No. of angular interval (1-180)
C
      COMMON /EGRIDA/ EGRIDA(119)
      DIMENSION SPEC(126,30,180)
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
C
C Check of input:
      IF(E.le.0.)THEN
        CONTINUE!WRITE(25,1000)E,NUMTRE
 1000   FORMAT(' HIST3AGT DIAGNOSTIC: E<0, E=',E10.3,',  NUMTRE=',I8)
        RETURN
      END IF
C Registration for all energies E>0.
        SPEC(125,NP,NA)=SPEC(125,NP,NA)+W
        SPEC(126,NP,NA)=SPEC(126,NP,NA)+W*E
C Registration for energies E<A.
      IF(E.lt.EGRIDA(1))THEN
        SPEC(119,NP,NA)=SPEC(119,NP,NA)+W
        SPEC(120,NP,NA)=SPEC(120,NP,NA)+W*E
        RETURN
      END IF
C Registration for energies E>=B.
      IF(E.ge.EGRIDA(119))THEN
        SPEC(123,NP,NA)=SPEC(123,NP,NA)+W
        SPEC(124,NP,NA)=SPEC(124,NP,NA)+W*E
        RETURN
      END IF
C Registration for energies in the interval 0<E<1E+6
        SPEC(121,NP,NA)=SPEC(121,NP,NA)+W
        SPEC(122,NP,NA)=SPEC(122,NP,NA)+W*E
        CALL SEARCH(E,EGRIDA,119,N1,N2)
        SPEC(N1,NP,NA)=SPEC(N1,NP,NA)+W
C
      RETURN
      END



      SUBROUTINE HIST180(JPART,ENER,WGHT,NANGL)
C Histogramming of the particles in each 1 degree interval
C INPUT: JPART - Type of particle
C        ENER  - Energy MeV/A
C        WGHT - Statistical weight
C        NANGL   - No. of angular interval (step 1 degree)
C
C      COMMON /HIPROJ/ APROJ,ZPROJ
      COMMON /HIPROJAGT/ APROJ,ZPROJ                       ! HISTAGT
      COMMON /BOOKFR/ BOOKFR(2,22)  ! (Zi,Ai)x22
      COMMON /SPECAGT/ SPECAGT(126,30,180)
C
C Particles n, p, pi-, pi+, d, t, He3, He4, are registered always.
      IF(JPART.eq.1)THEN                                   ! Neutron
        CALL HIST3AGT(ENER,SPECAGT, 1,NANGL,WGHT)
      END IF
C
      IF(JPART.eq.2)THEN                                   ! proton
        CALL HIST3AGT(ENER,SPECAGT, 2,NANGL,WGHT)
      END IF
C
      IF(JPART.eq.3)THEN                                   ! pi-
        CALL HIST3AGT(ENER,SPECAGT, 3,NANGL,WGHT)
      END IF
C
      IF(JPART.eq.4)THEN                                   ! pi+
        CALL HIST3AGT(ENER,SPECAGT, 4,NANGL,WGHT)
      END IF
C
      IF(JPART.eq.21)THEN                                      ! H2
        CALL HIST3AGT(ENER,SPECAGT, 5,NANGL,WGHT)
      END IF
      IF(JPART.eq.25 .AND. APROJ.eq.2. .AND. ZPROJ.eq.1.)THEN  ! H2
        CALL HIST3AGT(ENER,SPECAGT, 5,NANGL,WGHT)
      END IF
C
      IF(JPART.eq.22)THEN                                      ! H3
        CALL HIST3AGT(ENER,SPECAGT, 6,NANGL,WGHT)
      END IF
      IF(JPART.eq.25 .AND. APROJ.eq.3. .AND. ZPROJ.eq.1.)THEN  ! H3
        CALL HIST3AGT(ENER,SPECAGT, 6,NANGL,WGHT)
      END IF
C
      IF(JPART.eq.23)THEN                                      ! He3
        CALL HIST3AGT(ENER,SPECAGT, 7,NANGL,WGHT)
      END IF
      IF(JPART.eq.25 .AND. APROJ.eq.3. .AND. ZPROJ.eq.2.)THEN  ! He3
        CALL HIST3AGT(ENER,SPECAGT, 7,NANGL,WGHT)
      END IF
C
      IF(JPART.eq.24)THEN                                      ! He4
        CALL HIST3AGT(ENER,SPECAGT, 8,NANGL,WGHT)
      END IF
      IF(JPART.eq.25 .AND. APROJ.eq.4. .AND. ZPROJ.eq.2.)THEN  ! He4
        CALL HIST3AGT(ENER,SPECAGT, 8,NANGL,WGHT)
      END IF
C
C
C The following 22 (ten) fragments are registered according to booking of user
C in the array /BOOKFR/ BOOKFR(2,22): BOOKFR(1,i)=Zi, BOOKFR(2,i)=Ai. 
      DO J=1,22
        IF(JPART.eq.25)THEN                  ! Ion
          IF(BOOKFR(2,J) .eq. 1000.)THEN
C              Fixed Z and any A
            IF(ZPROJ.eq.BOOKFR(1,J))
     *      CALL HIST3AGT(ENER,SPECAGT,8+J,NANGL,WGHT)
          ELSE
            IF(BOOKFR(1,J) .eq. 1000.)THEN
C             Fixed A and any Z
              IF(APROJ.eq.BOOKFR(2,J))
     *        CALL HIST3AGT(ENER,SPECAGT,8+J,NANGL,WGHT)
            ELSE
C             Fixed (Z,A)
              IF(APROJ.eq.BOOKFR(2,J) .AND. ZPROJ.eq.BOOKFR(1,J))
     *        CALL HIST3AGT(ENER,SPECAGT,8+J,NANGL,WGHT)
            END IF
          END IF
        END IF
      END DO
C
      RETURN
      END



      SUBROUTINE STARAGT(NEVENT,IPR)                  ! HISTAGT
C STARAGT is modified program STAR(NEVENT) from the file GENAGT.F
C for scoring of inelastic interaction star. See label HISTAGT.
C
C INPUT: NEVENT - No. of the event
C        IPR - printout of star if(IPR.eq.1)
C
C       SECONDARY PARTICLES CHARACTERISTICS IN ARRAY SPT
C           SPT(1,J) - COS OF POLAR ANGLE OF SECONDARY PARTICLE
C           SPT(2,J), SPT(3,J) - SIN, COS OF AZIMUTAL ANGLE
C           SPT(4,J) - KINETIC ENERGY OF SECONDARY PARTICLE (MEV)
C           SPT(5,J) - SECONDARY PARTICLE TYPE
C           SPT(6,J) - SECONDARY PARTICLE WEIGHT
C       NUCLEUS CHARACTERISTICS
C           SNU(1,1) - NUCLEUS-TARGET ATOMIC WEIGHT
C           SNU(2,1) - NUCLEUS-TARGET ATOMIC NUMBER
C       Nuclei after deexitation
C           SNU(1,J) - A
C           SNU(2,J) - Z
C           SNU(3,J) - Excitation energy (MeV)
C           SNU(4,J) - KINTER
C           SNU(5,J) ]
C           SNU(6,J) ]   -  Nucleus momentum (MeV/C)
C           SNU(7,J) ]
C-------------------------------------------------------------------------
      COMMON /SECHAR/ SPT(6,5000),SNU(10,101),LS6,LS100,LS10,LS11  ! GSITR
      COMMON /HIPROJAGT/ APROJ,ZPROJ                       ! HISTAGT
C
      DIMENSION PARTIC(24)
      CHARACTER*7 PARTIC
      DATA PARTIC /'NEUTRON',' PROTON','  PI-  ','  PI+  ','  PI0  ',
     *             ' ANTI N',' ANTI P','  K-   ','  K+   ','  K0   ',
     *             'ANTI K0',' GAMMA ','ELECTRN','POSITRN','  MU-  ',
     *             '  MU+  ','  NU E ','ANTINUE',' NU MU ','ANTINUM',
     *             'DEUTRON','TRITIUM','HELIUM3',' ALPHA '/
C
      IF(IPR.eq.1)CONTINUE!WRITE(36,10)NEVENT                        ! HISTAGT 
   10 FORMAT(/ 20X,'COMPLETED ===STARAGT=== No.',I3)        ! HISTAGT   
      IF(IPR.eq.1)CONTINUE!WRITE(36,20)                              ! HISTAGT
   20 FORMAT(4X,'  Px(MeV/c)  Py(MeV/c)  Pz(MeV/c)    Ekin(MeV)', 
     *          ' PARTICLE')
C
C Variables for check of conservation laws
      ESUM=0.  ! Summary kinetic energy
      RMSUM=0. ! Summary mass of MESONS only!
      QSUM=0.  ! Summary charge Q
      BSUM=0.  ! Summary barion number
      SSUM=0.  ! Summary strangeness
      RLSUM=0. ! Summary lepton number
      PXSUM=0. ! Summary Px of all particles
      PYSUM=0. ! Summary Py of all particles
      PZSUM=0. ! Summary Pz of all particles
C
C ------------ Loop on all particles in SPT ---------------
      DO 100 J=1,LS100
        IPART=IFIX(SPT(5,J))          ! type of particle
      IF(IPART.EQ.0)GO TO 200         ! go out from the loop
        CT=SPT(1,J)
          NANGL=INT(ACOS(CT)*(180.0/3.1415926536))+1        ! HISTAGT
        ST=SQRT(1.0-CT*CT)
        SF=SPT(2,J)
        CF=SPT(3,J)
            E=SPT(4,J)
            W=SPT(6,J)
            W=1.                                            ! HISTAGT
C
C         Conservation laws
        IF(IPART.eq. 1)THEN    ! Neutron
            ESUM=ESUM+E
            RMSUM=RMSUM+0.
            QSUM=QSUM+0.
            BSUM=BSUM+1.
            SSUM=SSUM+0.
            RLSUM=RLSUM+0.
              PMOD=SQRT(E*(E+2.0*940.))
              PX=PMOD*ST*CF
              PY=PMOD*ST*SF
              PZ=PMOD*CT
                PXSUM=PXSUM+PX
                PYSUM=PYSUM+PY
                PZSUM=PZSUM+PZ
            CALL HIST180(IPART,E,W,NANGL)                   ! HISTAGT
        END IF
        IF(IPART.eq. 2)THEN    ! Proton
            ESUM=ESUM+E
            RMSUM=RMSUM+0.
            QSUM=QSUM+1.
            BSUM=BSUM+1.
            SSUM=SSUM+0.
            RLSUM=RLSUM+0.
              PMOD=SQRT(E*(E+2.0*940.))
              PX=PMOD*ST*CF
              PY=PMOD*ST*SF
              PZ=PMOD*CT
                PXSUM=PXSUM+PX
                PYSUM=PYSUM+PY
                PZSUM=PZSUM+PZ
            CALL HIST180(IPART,E,W,NANGL)                   ! HISTAGT
        END IF
        IF(IPART.eq. 3)THEN    ! PI-
            ESUM=ESUM+E
            RMSUM=RMSUM+140.
            QSUM=QSUM-1.
            BSUM=BSUM+0.
            SSUM=SSUM+0.
            RLSUM=RLSUM+0.
              PMOD=SQRT(E*(E+2.0*140.))
              PX=PMOD*ST*CF
              PY=PMOD*ST*SF
              PZ=PMOD*CT
                PXSUM=PXSUM+PX
                PYSUM=PYSUM+PY
                PZSUM=PZSUM+PZ
            CALL HIST180(IPART,E,W,NANGL)                   ! HISTAGT
        END IF
        IF(IPART.eq. 4)THEN    ! PI+
            ESUM=ESUM+E
            RMSUM=RMSUM+140.
            QSUM=QSUM+1.
            BSUM=BSUM+0.
            SSUM=SSUM+0.
            RLSUM=RLSUM+0.
              PMOD=SQRT(E*(E+2.0*140.))
              PX=PMOD*ST*CF
              PY=PMOD*ST*SF
              PZ=PMOD*CT
                PXSUM=PXSUM+PX
                PYSUM=PYSUM+PY
                PZSUM=PZSUM+PZ
            CALL HIST180(IPART,E,W,NANGL)                   ! HISTAGT
        END IF
        IF(IPART.eq. 5)THEN    ! PI0
            ESUM=ESUM+E
            RMSUM=RMSUM+140.
            QSUM=QSUM+0.
            BSUM=BSUM+0.
            SSUM=SSUM+0.
            RLSUM=RLSUM+0.
              PMOD=SQRT(E*(E+2.0*140.))
              PX=PMOD*ST*CF
              PY=PMOD*ST*SF
              PZ=PMOD*CT
                PXSUM=PXSUM+PX
                PYSUM=PYSUM+PY
                PZSUM=PZSUM+PZ
        END IF
        IF(IPART.eq. 6)THEN    ! AntiN
            ESUM=ESUM+E
            RMSUM=RMSUM+0.
            QSUM=QSUM+0.
            BSUM=BSUM-1.
            SSUM=SSUM+0.
            RLSUM=RLSUM+0.
              PMOD=SQRT(E*(E+2.0*940.))
              PX=PMOD*ST*CF
              PY=PMOD*ST*SF
              PZ=PMOD*CT
                PXSUM=PXSUM+PX
                PYSUM=PYSUM+PY
                PZSUM=PZSUM+PZ
        END IF
        IF(IPART.eq. 7)THEN    ! AntiP
            ESUM=ESUM+E
            RMSUM=RMSUM+0.
            QSUM=QSUM-1.
            BSUM=BSUM-1.
            SSUM=SSUM+0.
            RLSUM=RLSUM+0.
              PMOD=SQRT(E*(E+2.0*940.))
              PX=PMOD*ST*CF
              PY=PMOD*ST*SF
              PZ=PMOD*CT
                PXSUM=PXSUM+PX
                PYSUM=PYSUM+PY
                PZSUM=PZSUM+PZ
        END IF
        IF(IPART.eq. 8)THEN    ! K-
            ESUM=ESUM+E
            RMSUM=RMSUM+495.
            QSUM=QSUM-1.
            BSUM=BSUM+0.
            SSUM=SSUM-1.
            RLSUM=RLSUM+0.
              PMOD=SQRT(E*(E+2.0*495.))
              PX=PMOD*ST*CF
              PY=PMOD*ST*SF
              PZ=PMOD*CT
                PXSUM=PXSUM+PX
                PYSUM=PYSUM+PY
                PZSUM=PZSUM+PZ
        END IF
        IF(IPART.eq. 9)THEN    ! K+
            ESUM=ESUM+E
            RMSUM=RMSUM+495.
            QSUM=QSUM+1.
            BSUM=BSUM+0.
            SSUM=SSUM+1.
            RLSUM=RLSUM+0.
              PMOD=SQRT(E*(E+2.0*495.))
              PX=PMOD*ST*CF
              PY=PMOD*ST*SF
              PZ=PMOD*CT
                PXSUM=PXSUM+PX
                PYSUM=PYSUM+PY
                PZSUM=PZSUM+PZ
        END IF
        IF(IPART.eq.10)THEN    ! K0
            ESUM=ESUM+E
            RMSUM=RMSUM+495.
            QSUM=QSUM+0.
            BSUM=BSUM+0.
            SSUM=SSUM+1.
            RLSUM=RLSUM+0.
              PMOD=SQRT(E*(E+2.0*495.))
              PX=PMOD*ST*CF
              PY=PMOD*ST*SF
              PZ=PMOD*CT
                PXSUM=PXSUM+PX
                PYSUM=PYSUM+PY
                PZSUM=PZSUM+PZ
        END IF
        IF(IPART.eq.11)THEN    ! AntiK0
            ESUM=ESUM+E
            RMSUM=RMSUM+0.
            QSUM=QSUM+0.
            BSUM=BSUM+0.
            SSUM=SSUM-1.
            RLSUM=RLSUM+0.
              PMOD=SQRT(E*(E+2.0*495.))
              PX=PMOD*ST*CF
              PY=PMOD*ST*SF
              PZ=PMOD*CT
                PXSUM=PXSUM+PX
                PYSUM=PYSUM+PY
                PZSUM=PZSUM+PZ
        END IF
        IF(IPART.eq.12)THEN    ! Gamma
            ESUM=ESUM+E
            RMSUM=RMSUM+0.
            QSUM=QSUM+0.
            BSUM=BSUM+0.
            SSUM=SSUM+0.
            RLSUM=RLSUM+0.
              PMOD=SQRT(E*(E+2.0*0.0))
              PX=PMOD*ST*CF
              PY=PMOD*ST*SF
              PZ=PMOD*CT
                PXSUM=PXSUM+PX
                PYSUM=PYSUM+PY
                PZSUM=PZSUM+PZ
        END IF
        IF(IPART.eq.21)THEN    ! Deuteron
            ESUM=ESUM+E
            RMSUM=RMSUM+0.
            QSUM=QSUM+1.
            BSUM=BSUM+2.
            SSUM=SSUM+0.
            RLSUM=RLSUM+0.
              PMOD=SQRT(E*(E+2.0*1880.))
              PX=PMOD*ST*CF
              PY=PMOD*ST*SF
              PZ=PMOD*CT
                PXSUM=PXSUM+PX
                PYSUM=PYSUM+PY
                PZSUM=PZSUM+PZ
            CALL HIST180(IPART,E/2.,W,NANGL)                ! HISTAGT
        END IF
        IF(IPART.eq.22)THEN    ! Triton
            ESUM=ESUM+E
            RMSUM=RMSUM+0.
            QSUM=QSUM+1.
            BSUM=BSUM+3.
            SSUM=SSUM+0.
            RLSUM=RLSUM+0.
              PMOD=SQRT(E*(E+2.0*2820.))
              PX=PMOD*ST*CF
              PY=PMOD*ST*SF
              PZ=PMOD*CT
                PXSUM=PXSUM+PX
                PYSUM=PYSUM+PY
                PZSUM=PZSUM+PZ
            CALL HIST180(IPART,E/3.,W,NANGL)                ! HISTAGT
        END IF
        IF(IPART.eq.23)THEN    ! He3
            ESUM=ESUM+E
            RMSUM=RMSUM+0.
            QSUM=QSUM+2.
            BSUM=BSUM+3.
            SSUM=SSUM+0.
            RLSUM=RLSUM+0.
              PMOD=SQRT(E*(E+2.0*2820.))
              PX=PMOD*ST*CF
              PY=PMOD*ST*SF
              PZ=PMOD*CT
                PXSUM=PXSUM+PX
                PYSUM=PYSUM+PY
                PZSUM=PZSUM+PZ
            CALL HIST180(IPART,E/3.,W,NANGL)                ! HISTAGT
        END IF
        IF(IPART.eq.24)THEN    ! Alpha
            ESUM=ESUM+E
            RMSUM=RMSUM+0.
            QSUM=QSUM+2.
            BSUM=BSUM+4.
            SSUM=SSUM+0.
            RLSUM=RLSUM+0.
              PMOD=SQRT(E*(E+2.0*3760.))
              PX=PMOD*ST*CF
              PY=PMOD*ST*SF
              PZ=PMOD*CT
                PXSUM=PXSUM+PX
                PYSUM=PYSUM+PY
                PZSUM=PZSUM+PZ
            CALL HIST180(IPART,E/4.,W,NANGL)                ! HISTAGT
        END IF
C
        IF(IPR.eq.1)CONTINUE!WRITE(36,30)J,PX,PY,PZ,E,PARTIC(IPART)  ! HISTAGT
   30     FORMAT(I3,1X,3E11.4,E13.4,2X,A7)            
  100 CONTINUE
C -------------------- end of SPT loop -----------------
C
  200 CONTINUE
      IF(IPR.eq.1)CONTINUE!WRITE(36,39)                              ! HISTAGT
   39 FORMAT(66('-'))
      IF(IPR.eq.1)CONTINUE!WRITE(36,40)PXSUM,PYSUM,PZSUM,ESUM,RMSUM, ! HISTAGT
!     *            INT(QSUM),INT(BSUM),INT(SSUM),INT(RLSUM)
   40 FORMAT('Sum:',3E11.4,E13.4,'   Q   B   S   L'/
     *      30X,'MMESON=',E13.4,4I4)
C
C RESIDUAL NUCLEI CHARACTERISTICS
      IF(IPR.eq.1)CONTINUE!WRITE(36,50)                              ! HISTAGT
   50 FORMAT('     Residual nuclei:',34X,'A*   Z*    U*(MeV)')
C
      IF(IPR.eq.1)CONTINUE!WRITE(36,61)1,SNU(1,1),SNU(2,1)           ! HISTAGT
   61 FORMAT(I3,1X,32X,'Nucleus-target  ',2F5.1)
C
      DO 80 J=2,LS11                       ! Nuclei-products
        KINTER=IFIX(SNU(4,J))
        IF(KINTER.EQ.0)RETURN
        IF(KINTER.EQ.6)THEN
          IF(IPR.eq.1)CONTINUE!WRITE(36,62)J,SNU(1,J),SNU(2,J)       ! HISTAGT
   62               FORMAT(I3,1X,27X,'Full desintegration  ',2F5.1)
          GOTO 80      
        END IF
        EKIN=SQRT((SNU(5,J)**2+SNU(6,J)**2+SNU(7,J)**2)+
     *              (940.0*SNU(1,J))**2)-940.0*SNU(1,J)
C
C===================== HISTAGT ============================
        APROJ=SNU(1,J)      ! A of residual nucleus
        ZPROJ=SNU(2,J)      ! Z of residual nucleus
        EKINA=EKIN/APROJ    ! Enrgy of residual nucleus, MeV/A
        PXN=SNU(5,J)        ! Px of residual nucleus
        PYN=SNU(6,J)        ! Py of residual nucleus
        PZN=SNU(7,J)        ! Pz of residual nucleus
        CALL DECARS(PXN,PYN,PZN,PMOD,COST,SINF,COSF)
        NANGL=INT(ACOS(COST)*(180.0/3.1415926536))+1
        CALL HIST180(25,EKINA,W,NANGL)
C===================== end HISTAGT ========================
C
        IF(IPR.eq.1)CONTINUE!WRITE(36,70)J,SNU(5,J),SNU(6,J),SNU(7,J),EKIN, ! HISTAGT
!     *                            SNU(1,J),SNU(2,J),SNU(3,J)
   70     FORMAT(I3,1X,3E11.4,E13.4,2X,2F5.1,E11.4)
   80 CONTINUE
C
      RETURN
      END



      SUBROUTINE BYEVENT(NEVENT)                  ! HADES
C Printout of stars event by event into the file FOR038
C INPUT: NEVENT - No. of the event
C
C SECONDARY PARTICLES CHARACTERISTICS IN ARRAY SPT
C      SPT(1,J) - COS OF POLAR ANGLE OF SECONDARY PARTICLE
C      SPT(2,J), SPT(3,J) - SIN, COS OF AZIMUTAL ANGLE
C      SPT(4,J) - KINETIC ENERGY OF SECONDARY PARTICLE (MEV)
C      SPT(5,J) - SECONDARY PARTICLE TYPE
C      SPT(6,J) - SECONDARY PARTICLE WEIGHT
C NUCLEUS CHARACTERISTICS
C      SNU(1,1) - NUCLEUS-TARGET ATOMIC WEIGHT
C      SNU(2,1) - NUCLEUS-TARGET ATOMIC NUMBER
C Nuclei after deexitation
C      SNU(1,J) - A
C      SNU(2,J) - Z
C      SNU(3,J) - Excitation energy (MeV)
C      SNU(4,J) - KINTER
C      SNU(5,J) ]
C      SNU(6,J) ]   -  Nucleus momentum (MeV/C)
C      SNU(7,J) ]
C
C JPART is the unique index of hadron/nucleus 
C JPART:  Particle:
C     1   NEUTRON 
C     2   PROTON  
C     3   PI-     
C     4   PI+     
C     5   PI0     
C     6   ANTI N  
C     7   ANTI P  
C     8   K-      
C     9   K+      
C    10   K0      
C    11   K0 TILDA
C    21   DEUTERON 
C    22   TRITON   
C    23   He3      
C    24   He4
C    25   HEAVY ION
c---------------------------------------------------------------
      COMMON /SECHAR/ SPT(6,5000),SNU(10,101),LS6,LS100,LS10,LS11  ! GSITR
      COMMON /RIMPACT/ RIMPACT
      REAL*8 RIMPACT
      COMMON /BIMP/ B00,BX,BY
      REAL*8 B00,BX,BY
C
C##################Counter of all particles and fragments
      JLAST=0
      DO J=1,LS100
        JPART=INT(SPT(5,J))          ! type of particle
        IF(JPART.EQ.0)GOTO 1         ! go out from the loop
        JLAST=JLAST+1 ! counter of hadrons and of d,t,He3,He4
      END DO
    1 CONTINUE
      JFRAG=0
      DO 2 J=2,LS11                       ! Nuclei-products
        KINTER=IFIX(SNU(4,J))
        IF(KINTER.EQ.0)GOTO 3
        IF(KINTER.EQ.6)GOTO 2
        JFRAG=JFRAG+1
    2 END DO
    3 CONTINUE
C###########End of Counter of all particles and fragments
C
      WRITE(38,10)NEVENT,B00,JLAST+JFRAG
   10 FORMAT(/' Event',I8,'   Impact Parameter',F12.3,
     *        '   No. of secondaries',I5)
      WRITE(38,20)
   20 FORMAT(' No. JPR   A   Z',
     *       '     COST          SINF          COSF         MeV/A')
C
C ------------ Loop on all particles in SPT ---------------
      JLAST=0
      DO J=1,LS100
        JPART=INT(SPT(5,J))          ! type of particle
          IF(JPART.EQ.0)GOTO 200     ! go out from the loop
        JLAST=JLAST+1 ! counter of hadrons and of d,t,He3,He4
        CT=SPT(1,J)   ! Cos(Theta)
        SF=SPT(2,J)   ! Sin(Phi)
        CF=SPT(3,J)   ! Cos(Phi)
          AA=1.       ! For hadrons JPART=1-11 mass number AA=1
          ZZ=0.       ! For hadrons JPART=1-11 charge ZZ=0 (not essential)
        ET=SPT(4,J)   ! Total kinetic energy (MeV)
        EA=ET/AA      ! Kinetic energy per nucleon (MeV/A)
C
        IF(JPART.eq.21)THEN    ! Deuteron
          AA=2.
          ZZ=1.
          EA=ET/AA
        END IF
        IF(JPART.eq.22)THEN    ! Tritium
          AA=3.
          ZZ=1.
          EA=ET/AA
        END IF
        IF(JPART.eq.23)THEN    ! He3
          AA=3.
          ZZ=2.
          EA=ET/AA
        END IF
        IF(JPART.eq.24)THEN    ! He4
          AA=4.
          ZZ=2.
          EA=ET/AA
        END IF
C
        WRITE(38,30)J,JPART,INT(AA),INT(ZZ),CT,SF,CF,EA
   30   FORMAT(4I4,4E14.7)            
      END DO
C -------------------- end of SPT loop -----------------
C
  200 CONTINUE
      JPART=25
      JFRAG=0
C RESIDUAL NUCLEI CHARACTERISTICS
      DO 80 J=2,LS11                       ! Nuclei-products
        KINTER=IFIX(SNU(4,J))
        IF(KINTER.EQ.0)RETURN
        IF(KINTER.EQ.6)GOTO 80
        JFRAG=JFRAG+1
        ET=SQRT((SNU(5,J)**2+SNU(6,J)**2+SNU(7,J)**2)+
     *              (940.0*SNU(1,J))**2)-940.0*SNU(1,J)
C
        AA=SNU(1,J)      ! A of residual nucleus
        ZZ=SNU(2,J)      ! Z of residual nucleus
        EA=ET/AA         ! Enrgy of residual nucleus, MeV/A
        PXN=SNU(5,J)        ! Px of residual nucleus
        PYN=SNU(6,J)        ! Py of residual nucleus
        PZN=SNU(7,J)        ! Pz of residual nucleus
        CALL DECARS(PXN,PYN,PZN,PMOD,CT,SF,CF)
        WRITE(38,30)JFRAG+JLAST,JPART,INT(AA),INT(ZZ),CT,SF,CF,EA
   80 END DO
C
      RETURN
      END



      SUBROUTINE SEARCH(E,ARE,N,N1,N2)
C BINARY SEARCH IN ARRANGED (IN INCREASING ORDER) ARRAY ARE(N)
C     INPUT:
C          E - CURRENT VARIABLE
C          ARE(N) - VARIABLE MESH
C     OUTPUT:
C          N1,N2  SUCH, THAT  ARE(N1) < E < ARE(N2)
C     WARNING!!!
C               1. ARE(1).LE.E .AND. E.LE.ARE(N)
C               2. N.GT.2
C-------------------------------------------------------------------
C
      DIMENSION ARE(N)
C
      N1=1
      N2=N
    1 NM=(N1+N2)/2
      IF(E.LE.ARE(NM))THEN
          N2=NM
      ELSE
          N1=NM
      END IF
      IF(N2.EQ.N1+1)RETURN
      GO TO 1
      END



      SUBROUTINE MAKEGRID
C This program provides actions for arrangement of output 
C to binary file FOR037.
C 1.Setting up of the grid for registration of enegy spectra
C   of secondary particles and fragments.
C 2.Correction of the default list of secondery fragments.
C
      COMMON /EGRIDA/ EGRIDA(119)  ! Default energy grid
      COMMON /BOOKFR/ BOOKFR(2,22) ! Default (Zi,Ai) of fragments (> He4)
      CHARACTER*1 YN               ! yes, no
C
C
C SETTING UP OF THE ENERGY GRID
C MODIFIED BY A. Timofeev 18.07.2011
C Just one month before wedding X_x
C always use default settings (answer "y")
      GOTO 777
C END OF CHANGE
      write(*,10)
   10 format(/' Do you want use the default energy grid ? (y,n)')            
   11 read(*,12)YN
   12 format(A1)
      if(YN.ne.'Y' .AND. YN.ne.'y' .AND. YN.ne.'N' .AND. YN.ne.'n')
     *  GOTO 11
      if(YN.eq.'Y' .OR. YN.eq.'y')then
        CONTINUE!WRITE(34,13)
   13   FORMAT(/' MAKEGRID message: Default grid is used')
        GOTO 100  ! to list of fragments
      end if
C
      write(*,20)
   20 format(/' What type of the energy greed you need ? (1,2,3)'/
     *        '   1. Uniform'/
     *        '   2. Log10'/
     *        '   3. Geom progression')
   21 read(*,*)ITG
      if(ITG.lt.1 .OR. ITG.gt.3)GOTO 21
C
      write(*,30)
   30 format(/' How many bins Nbin you need ? (2<=Nbin<=118)')
   31 read(*,*)Nbin
      if(Nbin.lt.2 .OR. Nbin.gt.118)GOTO 31
C
      Nnod=Nbin+1  ! number of nodes
C
      IF(ITG.eq.1)CALL MAKEUNI(Nbin)
      IF(ITG.eq.2)CALL MAKELOG(Nbin)
      IF(ITG.eq.3)CALL MAKEGEP(Nbin)
C
      DO J=Nnod+1,119
        EGRIDA(J)=EGRIDA(Nnod)  ! To complete the grid up to 119 nodes
      END DO
C
C
C CORRECTION OF THE LIST OF SECONDARY FRAGMENT
  100 write(*,110)
  110 format(/' Do you want use the default list of fragments ? (y,n)')
  111 read(*,112)YN
  112 format(A1)
      if(YN.ne.'Y' .AND. YN.ne.'y' .AND. YN.ne.'N' .AND. YN.ne.'n')
     *  GOTO 111
      if(YN.eq.'Y' .OR. YN.eq.'y')then
        CONTINUE!WRITE(34,113)
  113   FORMAT(/' MAKEGRID message: Default list of fragments is used')
        GOTO 777
      end if
C
  119 write(*,120)
  120 format(/' Input (Zi,Ai). Ai=1000 for any A. Zi=1000 for any Z.'/)
      DO J=1,22
 1121   write(*,121)J+8,INT(BOOKFR(1,J)),INT(BOOKFR(2,J))
  121   format(' No.',I2,'  (Zi,Ai)=',2I5,'   is OK ? (y,n)')
        read(*,122)YN
  122   format(A1)
      if(YN.ne.'Y' .AND. YN.ne.'y' .AND. YN.ne.'N' .AND. YN.ne.'n')
     *  GOTO 1121
C
        if(YN.eq.'N' .OR. YN.eq.'n')then
          write(*,*)' Zi,Ai ?'
          read(*,*)BOOKFR(1,J),BOOKFR(2,J)
        end if
      END DO
C
      CONTINUE!WRITE(34,123)
      write( *,123)
  123 format(/' MAKEGRID message: List of fragments after correction:')
      DO J=1,22
        CONTINUE!WRITE(34,124)J+8,INT(BOOKFR(1,J)),INT(BOOKFR(2,J))
        write( *,124)J+8,INT(BOOKFR(1,J)),INT(BOOKFR(2,J))
  124   format(' No.',I2,'  (Zi,Ai)=',2I5)
      END DO
C
 1125 write(*,125)
  125 format(/' Is this list OK ? (y,n)')
      read(*,126)YN
  126 format(A1)
      if(YN.ne.'Y' .AND. YN.ne.'y' .AND. YN.ne.'N' .AND. YN.ne.'n')
     *  GOTO 1125
      if(YN.eq.'N' .OR. YN.eq.'n')GOTO 119
C
C MEMORIZING OF ENERGY GRID AND FRAGMENT LIST
  777 CONTINUE!WRITE(37) EGRIDA,BOOKFR
C
      RETURN
      END



      SUBROUTINE MAKEUNI(Nbin)
C Making uniform grid, Nbin from Emin to EmaX
      COMMON /EGRIDA/ EGRIDA(119)
C
      write(*,*)' Emin>=0, Emax for uniform grid ?'
   10 read(*,*)Emin, Emax
        if(Emin.lt.0. .OR. Emin.ge.Emax)GOTO 10
C
      STEP=(Emax-Emin)/FLOAT(Nbin)
C
      DO J=1,Nbin+1
        EGRIDA(J)=Emin+STEP*FLOAT(J-1)
      END DO
C
      CONTINUE!WRITE(34,20)Emin,Emax,Nbin
   20 FORMAT(/' MAKEGRID message: Uniform greed is defined:'/
     *        ' Emin=',E10.4,',  Emax=',E10.4,',  Nbin=',I3)
      RETURN
      END



      SUBROUTINE MAKELOG(Nbin)
C Making Log10 grid, Nbin from Emin to EmaX
      COMMON /EGRIDA/ EGRIDA(119)
C
      write(*,*)' Emin>0, Emax for Log10 grid ?'
   10 read(*,*)Emin, Emax
        if(Emin.le.0. .OR. Emin.ge.Emax)GOTO 10
C
      RELmin=LOG10(Emin)
      RELmax=LOG10(Emax)
      STEPLOG=(RELmax-RELmin)/FLOAT(Nbin)
C
      DO J=1,Nbin+1
        EGRIDA(J)=10.0**(RELmin+STEPLOG*FLOAT(J-1))
      END DO
C
      CONTINUE!WRITE(34,20)Emin,Emax,Nbin
   20 FORMAT(/' MAKEGRID message: Log10 greed is defined:'/
     *        ' Emin=',E10.4,',  Emax=',E10.4,',  Nbin=',I3)
      RETURN
      END



      SUBROUTINE MAKEGEP(Nbin)
C Making Geom progression grid, Nbin from Emin to EmaX
      COMMON /EGRIDA/ EGRIDA(119)
      REAL*8 GPC,GPCN
C
      write(*,*)' Emin>=0, Emax for Geom progression grid ?'
   10 read(*,*)Emin, Emax
        if(Emin.lt.0. .OR. Emin.ge.Emax)GOTO 10
C
      GPC=1.0D0+DBLE(1.0D0/DFLOAT(Nbin))
      GPCN=GPC**Nbin
      DELTA=(Emax-Emin)/((1.-SNGL(GPCN))/(1.0-SNGL(GPC)))
C
      DO J=1,Nbin+1
        IF(J.eq.1)THEN
          EGRIDA(J)=Emin
        ELSE
          EGRIDA(J)=EGRIDA(J-1)+DELTA*SNGL(GPC**(J-2))
        END IF
      END DO      
C      
      CONTINUE!WRITE(34,20)Emin,Emax,Nbin
   20 FORMAT(/' MAKEGRID message: GeomProgressio greed is defined:'/
     *        ' Emin=',E10.4,',  Emax=',E10.4,',  Nbin=',I3)
      RETURN
      END
