
***********************************************************************
*                                                                     *
*        The following subroutines are written by N.M.Sobolevsky,     *
*        Institute for Nuclear Research RAS, Moscow,                  *
*        E-mail: sobolevs@AL20.inr.troitsk.ru                         *
*                                                                     *
***********************************************************************

      SUBROUTINE PRDBUG
C DEBUGGING PRINT OF TREE,SOURN,SOURG AND TRAN
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      COMMON /TREE/ TREE(15,400000),ITREE(7,400000),LTR13,LTR7,LTR300 ! GSITR
C
      COMMON /SOURN/ SOURN(13,100000),ISOURN(7,100000),LSN13,LSN7,LSN100   ! GSITR
      COMMON /SOURG/ SOURG(13,20000),ISOURG(7,20000),LSG13,LSG7,LSG100  ! GSITR
      COMMON /SEMIN/ SEMIN(13,100),ISEMIN(7,100),L1313,L137,L13100
      COMMON /SEPLS/ SEPLS(13,100),ISEPLS(7,100),L1413,L147,L14100
      COMMON /SNUEL/ SNUEL(13,100),ISNUEL(7,100),L1713,L177,L17100
      COMMON /SANUE/ SANUE(13,100),ISANUE(7,100),L1813,L187,L18100
      COMMON /SNUMU/ SNUMU(13,100),ISNUMU(7,100),L1913,L197,L19100
      COMMON /SANUM/ SANUM(13,100),ISANUM(7,100),L2013,L207,L20100
C
      COMMON /MANALL/ NTREE,LINE1,LINE2,NBRMAX,NPRMAX,NNEUTR,NGAMMA,
     *                MELSC,MINTER,MCUT,MOUT,MABSOR,MDECAY,NSTRAN,
     *                NEMIN,NEPLS,NNUEL,NANUE,NNUMU,NANUM
C
      COMMON /TRAN/ TRAN(5,600000),ITRAN(6,600000),LSR5,LSR6,LSR500  ! GSITR
C
      WRITE(24,100)
  100 FORMAT(10X,'TREE')
      WRITE(24,121)
  121 format('   ( 1)   ( 2)   ( 3)   ( 4)   ( 5)   ( 6)   ( 7)   ( 8)',
     *       '   ( 9)   (10)   (11)   (12)   (13)   (14)   (15)',
     *       ' JPR NTR NBR "4" KGN KPS KST')
      WRITE(24,122)
  122 format(6X,'X',6X,'Y',6X,'Z',
     *       '  Cos/A  Sin/Z  Cos/0',
     *       '    T/?   W/Px   ?/Py   ?/Pz',
     *       '     11     12     13     14     15',
     *       ' JPR NTR NBR "4" KGN KPS KST')
      DO 10 K=1,LINE2
c   10 WRITE(24,1)(TREE(J,K),J=1,8),(ITREE(J,K),J=1,7)
   10 WRITE(24,1)(TREE(J,K),J=1,15),(ITREE(J,K),J=1,7)
c    1 FORMAT(3F6.2,3F7.3,F7.1,F6.3,7I4)
    1 FORMAT(3F7.2,3F7.3,9F7.1,7I4)
C
      WRITE(24,101)
  101 FORMAT(10X,'SOURN')
      DO 11 K=1,NNEUTR
   11 WRITE(24,1)(SOURN(J,K),J=1,15),(ISOURN(J,K),J=1,7)
      return
C
      WRITE(24,102)
  102 FORMAT(10X,'SOURG')
      DO 12 K=1,NGAMMA
   12 WRITE(24,1)(SOURG(J,K),J=1,8),(ISOURG(J,K),J=1,7)
C
      WRITE(24,103)
  103 FORMAT(10X,'SEMIN')
      DO 13 K=1,NEMIN
   13 WRITE(24,1)(SEMIN(J,K),J=1,8),(ISEMIN(J,K),J=1,7)
C
      WRITE(24,104)
  104 FORMAT(10X,'SEPLS')
      DO 14 K=1,NEPLS
   14 WRITE(24,1)(SEPLS(J,K),J=1,8),(ISEPLS(J,K),J=1,7)
C
      WRITE(24,105)
  105 FORMAT(10X,'SNUEL')
      DO 15 K=1,NNUEL
   15 WRITE(24,1)(SNUEL(J,K),J=1,8),(ISNUEL(J,K),J=1,7)
C
      WRITE(24,106)
  106 FORMAT(10X,'SANUE')
      DO 16 K=1,NANUE
   16 WRITE(24,1)(SANUE(J,K),J=1,8),(ISANUE(J,K),J=1,7)
C
      WRITE(24,107)
  107 FORMAT(10X,'SNUMU')
      DO 17 K=1,NNUMU
   17 WRITE(24,1)(SNUMU(J,K),J=1,8),(ISNUMU(J,K),J=1,7)
C
      WRITE(24,108)
  108 FORMAT(10X,'SANUM')
      DO 18 K=1,NANUM
   18 WRITE(24,1)(SANUM(J,K),J=1,8),(ISANUM(J,K),J=1,7)
C
C
      WRITE(24,203)
  203 FORMAT(10X,'TRAN')
      DO 213 K=1,NSTRAN
  213 WRITE(24,2)(TRAN(J,K),J=1,5),(ITRAN(J,K),J=1,6)
    2 FORMAT(4F7.2,F8.1,5X,6I5)
C
      RETURN
      END



      SUBROUTINE PRINTR
C PRINTING OF SHIELD/CG EXTRANUCLEAR CESCADE TREE.
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      COMMON /TREE/ TREE(15,400000),ITREE(7,400000),LTR13,LTR7,LTR300 ! GSITR
C
      COMMON /SOURN/ SOURN(13,100000),ISOURN(7,100000),LSN13,LSN7,LSN100   ! GSITR
      COMMON /SOURG/ SOURG(13,20000),ISOURG(7,20000),LSG13,LSG7,LSG100  ! GSITR
      COMMON /SEMIN/ SEMIN(13,100),ISEMIN(7,100),L1313,L137,L13100
      COMMON /SEPLS/ SEPLS(13,100),ISEPLS(7,100),L1413,L147,L14100
      COMMON /SNUEL/ SNUEL(13,100),ISNUEL(7,100),L1713,L177,L17100
      COMMON /SANUE/ SANUE(13,100),ISANUE(7,100),L1813,L187,L18100
      COMMON /SNUMU/ SNUMU(13,100),ISNUMU(7,100),L1913,L197,L19100
      COMMON /SANUM/ SANUM(13,100),ISANUM(7,100),L2013,L207,L20100
C
      COMMON /TREVS/ TREVS(15,400),ITREVS(7,400),LVS13,LVS7,LVS100   !hion6
C
C Number of branches in tree
      NBRMAX=ITREE(4,1)
C Number of lines in tree (excluding first and last lines)
      NLINE =ITREE(3,1)
C Numbers of particles in sources
	  NNEUTR=ISOURN(LSN7,1)
	  NGAMMA=ISOURG(LSG7,1)
	  NEMIN =ISEMIN(L137,1)
	  NEPLS =ISEPLS(L147,1)
	  NNUEL =ISNUEL(L177,1)
	  NANUE =ISANUE(L187,1)
	  NNUMU =ISNUMU(L197,1)
	  NANUM =ISANUM(L207,1)
C
C Tree is empty. No cascade tree was generated. Incident particle
C was removed to source array (20.01.97).
      IF(NBRMAX.eq.0)THEN
	  WRITE(24,111)NUMTRE,IXINIT
  111     FORMAT(/10X,'SHIELD/CG EXTRANUCLEAR CASCADE TREE NO.',I5,
     *           ',    (IXINIT=',I12,')'/)
	  WRITE(24,112)
  112     FORMAT(5X,'NO CASCADE TREE WAS GENERATED. INCIDENT PARTICLE',
     *              ' WAS REMOVED TO SOURCE ARRAY.')
	  RETURN
      END IF
C
C PRINTING OF TREE TITLE
      CALL PRSTAR(0,0,0,0,0)
	  KST=ITREE(7,3)
      IF(KST.EQ.1 .OR. KST.EQ.2 .OR. KST.EQ.5 .OR. KST.EQ.6)THEN
	  WRITE(24,4)
    4     FORMAT()
      ELSE
	  RETURN
      END IF
C
C------------------ CYCLE ON ALL BRANCHES -----------------------------
      JBRST=2
      DO 100 JBR=1,NBRMAX
C
C         SEARCH OF BRANCH, CONTAINING ELASTIC OR INELASTIC INTERACTION
C         and decay. JBRST AND JBRFN - FIRST AND LAST LINES OF ARRAY
C         TREE, CONTAINING THIS BRANCH NO.JBR
		  KSTATE=ITREE(7,JBRST+1)
	      IF(KSTATE.LE.2)THEN
C             Scattering or inelastic interaction
C                 Counting of residual nuclei number NUC
		      NUC=0
		  DO JN=JBRST+3,NLINE+1
		      IF(ITREE(1,JN).EQ.0)THEN
			  NUC=NUC+1
		      ELSE
			  GO TO 10
		      END IF
		  END DO
   10             CONTINUE
C Debug check of residual nuclei number NUC
      IF(NUC.EQ.0)THEN
	  CONTINUE!WRITE(25,11)JBRST,JBRFN
   11     FORMAT(' PRINTR DIAGNOSTIC: THERE IS ONLY NUCLEUS-TARGET',
     *           ' IN BRANCH,'/20X,'JBRST=',I4,',   JBRFN=',I4)
	  CLOSE(25)
	  STOP
      END IF
C End debug check NUC ---------------------
		  JBRFN=JBRST+2+NUC
		  GO TO 20
	      ELSE
		  CONTINUE
	      END IF
C
	  IF(KSTATE.EQ.5 .OR. KSTATE.EQ.6)THEN
C             Decay on fly or at rest
	      JBRFN=JBRST+1
	      GO TO 200
	  ELSE
C             Fly out or absorption
	      JBRFN=JBRST+1
	      GO TO 99
	  END IF
C
C ELASTIC OR INELASTIC INTERACTION PRESENT BRANCH NO.JBR ?
   20     IF(KSTATE.EQ.2)GO TO 50
C
C                       ELASTIC SCATTERING
C SEARCH OF ELASTIC SCATTERED PARTICLE (before range, but allocation -
C                                       after range)
	  DO 30 K=JBRFN+1,NLINE+1
	  IF(ITREE(4,K).NE.JBR)GO TO 30
              DO 31 J=1,LTR13                                        !hion6
   31         TREVS(J,1)=TREE(J,K+1)
              DO 32 J=1,LTR7                                         !hion6
   32         ITREVS(J,1)=ITREE(J,K+1)
	      GO TO 80
   30     CONTINUE
C
C                       INELASTIC INTERACTION
   50     NPCAS=0
	  NNEUT=0
	  NPTOT=0
C
C SEARCH OF CASCADE PARTICLES (after range), CREATED IN INELASTIC
C INTERACTION NO.JBR
	  DO 60 K=JBRFN+1,NLINE+1
	  IF(ITREE(4,K).EQ.JBR .AND. ITREE(7,K).NE.-1)THEN
	      NPCAS=NPCAS+1
	      NPTOT=NPTOT+1
		  IF(NPTOT.GT.LVS100)THEN
		      CONTINUE!WRITE(25,1000)
 1000                 FORMAT(' PRINTR DIAGNOSTIC: ARRAY TREVS EXCEEDED')
		      GO TO 80
		  ELSE
		      CONTINUE
		  END IF
              DO 61 J=1,LTR13                                        !hion6
   61         TREVS(J,NPTOT)=TREE(J,K)
C                      Restoring of energy before range
		       TREVS(7,NPTOT)=TREE(7,K-1)
              DO 62 J=1,LTR7                                         !hion6
   62         ITREVS(J,NPTOT)=ITREE(J,K)
	      GO TO 60
	  ELSE
	      GO TO 60
	  END IF
   60     CONTINUE
C
C SEARCH OF LOW ENERGY NEUTRONS, CREATED IN INELASTIC INTERACTION NO.JBR
	  DO 70 K=2,NNEUTR+1
	  IF(ISOURN(3,K).NE.JBR)GO TO 70
	      NNEUT=NNEUT+1
	      NPTOT=NPTOT+1
		  IF(NPTOT.GT.LVS100)THEN
		      CONTINUE!WRITE(25,1000)
		      GO TO 80
		  ELSE
		      CONTINUE
		  END IF
              DO 71 J=1,LSN13                                        !hion6
   71         TREVS(J,NPTOT)=SOURN(J,K)
              DO 72 J=1,LSN7                                         !hion6
   72         ITREVS(J,NPTOT)=ISOURN(J,K)
   70     CONTINUE
		  GO TO 80
C
C                     Decay on fly and at rest
  200     NPCAS=0
	  NPSOU=0
	  NPTOT=0
C Search of decay products on whole tree and all sources (except 
C neutrons), created in branch No.JBR
	  DO 210 K=JBRFN+1,NLINE+1
	  IF(ITREE(4,K).EQ.JBR .AND. ITREE(7,K).NE.-1)THEN
	      NPCAS=NPCAS+1
	      NPTOT=NPTOT+1
              DO 211 J=1,LTR13                                       !hion6
  211         TREVS(J,NPTOT)=TREE(J,K)
		       TREVS(7,NPTOT)=TREE(7,K-1)
              DO 212 J=1,LTR7                                        !hion6
  212         ITREVS(J,NPTOT)=ITREE(J,K)
	      GO TO 210
	  ELSE
	      GO TO 210
	  END IF
  210     CONTINUE
C
C Search on gamma source
	  DO 220 K=2,NGAMMA+1
	  IF(ISOURG(3,K).NE.JBR)GO TO 220
	      NPSOU=NPSOU+1
	      NPTOT=NPTOT+1
              DO 221 J=1,LSG13                                       !hion6
  221         TREVS(J,NPTOT)=SOURG(J,K)
              DO 222 J=1,LSG7                                        !hion6
  222         ITREVS(J,NPTOT)=ISOURG(J,K)
  220     CONTINUE
C
C Search on electron source
	  DO 230 K=2,NEMIN+1
	  IF(ISEMIN(3,K).NE.JBR)GO TO 230
	      NPSOU=NPSOU+1
	      NPTOT=NPTOT+1
              DO 231 J=1,L1313                                       !hion6
  231         TREVS(J,NPTOT)=SEMIN(J,K)
              DO 232 J=1,L137                                        !hion6
  232         ITREVS(J,NPTOT)=ISEMIN(J,K)
  230     CONTINUE
C
C Search on positron source
	  DO 240 K=2,NEPLS+1
	  IF(ISEPLS(3,K).NE.JBR)GO TO 240
	      NPSOU=NPSOU+1
	      NPTOT=NPTOT+1
              DO 241 J=1,L1413                                       !hion6
  241         TREVS(J,NPTOT)=SEPLS(J,K)
              DO 242 J=1,L147                                        !hion6
  242         ITREVS(J,NPTOT)=ISEPLS(J,K)
  240     CONTINUE
C
C Search on electronic neitrino source
	  DO 250 K=2,NNUEL+1
	  IF(ISNUEL(3,K).NE.JBR)GO TO 250
	      NPSOU=NPSOU+1
	      NPTOT=NPTOT+1
              DO 251 J=1,L1713                                       !hion6
  251         TREVS(J,NPTOT)=SNUEL(J,K)
              DO 252 J=1,L177                                        !hion6
  252         ITREVS(J,NPTOT)=ISNUEL(J,K)
  250     CONTINUE
C
C Search on electronic antineitrino source
	  DO 260 K=2,NANUE+1
	  IF(ISANUE(3,K).NE.JBR)GO TO 260
	      NPSOU=NPSOU+1
	      NPTOT=NPTOT+1
              DO 261 J=1,L1813                                       !hion6
  261         TREVS(J,NPTOT)=SANUE(J,K)
              DO 262 J=1,L187                                        !hion6
  262         ITREVS(J,NPTOT)=ISANUE(J,K)
  260     CONTINUE
C
C Search on muonic neitrtno source
	  DO 270 K=2,NNUMU+1
	  IF(ISNUMU(3,K).NE.JBR)GO TO 270
	      NPSOU=NPSOU+1
	      NPTOT=NPTOT+1
              DO 271 J=1,L1913                                       !hion6
  271         TREVS(J,NPTOT)=SNUMU(J,K)
              DO 272 J=1,L197                                        !hion6
  272         ITREVS(J,NPTOT)=ISNUMU(J,K)
  270     CONTINUE
C
C Search on muonic antineitrino source
	  DO 280 K=2,NANUM+1
	  IF(ISANUM(3,K).NE.JBR)GO TO 280
	      NPSOU=NPSOU+1
	      NPTOT=NPTOT+1
              DO 281 J=1,L2013                                       !hion6
  281         TREVS(J,NPTOT)=SANUM(J,K)
              DO 282 J=1,L207                                        !hion6
  282         ITREVS(J,NPTOT)=ISANUM(J,K)
  280     CONTINUE
		  GO TO 80
C
C
   80     CALL PRSTAR(KSTATE,JBRST,JBRFN,NPCAS,NPTOT)
   99     JBRST=JBRFN+1
C
  100 CONTINUE
C------------------- END OF CYCLE ON JBR --------------------------------
      RETURN
      END
 
 
!  I cannot suppress writing to 24'th file, that's why I had comment whole subroutine
!  If i comment 'write' at label 52, i recive segmentation violation 
      SUBROUTINE PRSTAR(KSTATE,JBRST,JBRFN,NPCAS,NPTOT)
!C PRINTING OF TREE TITLE, ELASTIC AND INELASTIC INTERACTION STARS.
!C INPUT:
!C       KSTATE=0 - PRINT TITLE
!C              1 - PRINT SCATTERING STAR
!C              2 - PRINT INELASTIC STAR
!C              5,6 - print decay star
!C       JBRST - NO. OF LINE OF TREE, WHERE LOCATED THE FIRST LINE OF STAR.
!C       NPCAS - NUMBER OF PARTICLES IN INELASTIC STAR and decay star
!C                EXCLUDING particles from sources
!C       NPTOT - TOTAL NUMBER OF PARTICLES IN INELASTIC and decay STAR.
!C OUTPUT: IN FILE FOR024.DAT
!C--------------------------------------------------------------------------
!      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
!      COMMON /TREE/ TREE(15,400000),ITREE(7,400000),LTR13,LTR7,LTR300 ! GSITR
!      COMMON /TREVS/ TREVS(15,400),ITREVS(7,400),LVS13,LVS7,LVS100   !hion6
!C
!      DIMENSION PARTIC(25),STATUS(8)                                 !HION6
!      CHARACTER*7 PARTIC
!      DATA PARTIC /'NEUTRON',' PROTON','  PI-  ','  PI+  ','  PI0  ',
!     *             ' ANTI N',' ANTI P','  K-   ','  K+   ','  K0   ',
!     *             'ANTI K0',' GAMMA ','ELECTRN','POSITRN','  MU-  ',
!     *             '  MU+  ','  NU E ','ANTINUE',' NU MU ','ANTINUM',
!     *             'DEUTRON','TRITIUM','HELIUM3',' ALPHA ','HEV*ION'/  !HION6
!      CHARACTER*11 STATUS
!      DATA STATUS /' SCATTERING','INTERACTION','  FLY OUT  ',
!     *             ' ABSORPTION',' DECAY(FLY)','DECAY(REST)',
!     *             '  CUT OFF  ','  SOURCE   '/
!      CHARACTER*7 TITLPR,TITK0L,TITK0S
!      DATA TITK0L/'  K0L  '/ TITK0S/'  K0S  '/
!      
!	  IF(KSTATE.EQ.5 .OR. KSTATE.EQ.6)THEN
!	      KST1=4
!	  ELSE
!	      KST1=KSTATE+1
!	  END IF
!      GO TO (10,20,30,200),KST1
!C
!C------------------- PRINTING OF TREE TITLE -----------------------------
!C
!   10 WRITE(24,11)NUMTRE,IXINIT
!   11 FORMAT(/10X,'SHIELD/CG EXTRANUCLEAR CASCADE TREE NO.',I5,
!     *       ',    (IXINIT=',I12,')'/)
!C
!      JPART=ITREE(1,2)
!	  XI=TREE(1,2)
!	  YI=TREE(2,2)
!	  ZI=TREE(3,2)
!	      COST=TREE(4,2)
!	      SINF=TREE(5,2)
!	      COSF=TREE(6,2)
!		  TI=TREE(7,2)
!                      IAPR=INT(TREE(14,2))                           !HION6
!                      IZPR=INT(TREE(15,2))                           !HION6
!C
!      IF(JPART.EQ.10 .OR. JPART.EQ.11)THEN
!	  K0L=IFIX(TREE(10,2))
!	  K0S=IFIX(TREE(11,2))
!	  IF(K0L.EQ.1)TITLPR=TITK0L
!	  IF(K0S.EQ.1)TITLPR=TITK0S
!      ELSE
!	  TITLPR=PARTIC(JPART)
!      END IF
!C
!      if(JPART.eq.25)then                                            !HION6
!        WRITE(24,110)TITLPR,IAPR,IZPR,TI,XI,YI,ZI,COST,SINF,COSF     !HION6
!  110   FORMAT(' PRIMARY PARTICLE IS ',A7,'(',I3,',',I2,')',         !HION6
!     *  ', T=',F8.1,' MeV/A,',                                       !HION6
!     *  ' ENTRY POINT:',3(F7.2,','),' DIRECTION:',2(F7.4,','),F7.4)  !HION6
!      else                                                           !HION6
!        WRITE(24,12)TITLPR,TI,XI,YI,ZI,COST,SINF,COSF
!   12   FORMAT(' PRIMARY PARTICLE IS ',A7,', T=',F8.1,' MEV,',
!     *  ' ENTRY POINT:',3(F7.2,','),' DIRECTION:',2(F7.4,','),F7.4)
!      end if                                                         !HION6
!C
!      ISTATE=ITREE(7,3)
!      IF(ISTATE.LE.2 .OR. ISTATE.EQ.5 .OR. ISTATE.EQ.6)RETURN
!	  XF=TREE(1,3)
!	  YF=TREE(2,3)
!	  ZF=TREE(3,3)
!	      TF=TREE(7,3)
!C
!      WRITE(24,13)STATUS(ISTATE),XF,YF,ZF,TF
!   13 FORMAT(5X,A11,' OF PRIMARY PARTICLE IN POINT:',2(F7.2,','),F7.2,
!     *       '  AT ENERGY ',F9.2,' MEV. END OF TREE.')
!      RETURN
!C
!C---------------------- PRINTING OF ELASTIC STAR --------------------------
!   20 JPART=ITREE(1,JBRST)
!      NBRAN=ITREE(3,JBRST)
!	  XI=TREE(1,JBRST)
!	  YI=TREE(2,JBRST)
!	  ZI=TREE(3,JBRST)
!          TI=TREE(7,JBRST)
!            IAPR=INT(TREE(14,JBRST))                                 !HION6
!            IZPR=INT(TREE(15,JBRST))                                 !HION6
!	      XF=TREE(1,JBRST+1)
!	      YF=TREE(2,JBRST+1)
!	      ZF=TREE(3,JBRST+1)
!	      TF=TREE(7,JBRST+1)
!		  Z=TREE(5,JBRST+2)
!		  A=TREE(4,JBRST+2)
!		      XZ=TREVS(1,1)
!		      YZ=TREVS(2,1)
!		      ZZ=TREVS(3,1)
!		      TZ=TREVS(7,1)
!			  ISTATE=ITREVS(7,1)
!			  NBRNEW=ITREVS(3,1)
!C
!      if(JPART.eq.25)then                                            !HION6
!        WRITE(24,120)PARTIC(JPART),IAPR,IZPR,NBRAN                   !HION6
!  120   FORMAT(' ELASTIC SCATTERING OF ',A7,'(',I3,',',I2,')',       !HION6
!     *         ', MEV==>MeV/A, (NBRAN=',I4,'):')                     !HION6
!      else                                                           !HION6
!        WRITE(24,21)PARTIC(JPART),NBRAN
!   21   FORMAT(' ELASTIC SCATTERING OF ',A7,', (NBRAN=',I4,'):')
!      end if                                                         !HION6
!C
!      WRITE(24,22)PARTIC(JPART),XI,YI,ZI,TI
!   22 FORMAT(5X,A7,' BEFORE RANGE IN POINT:',2(F7.2,','),F7.2,
!     *       ' HAVE ENERGY ',F10.2,' MEV.')
!C
!      WRITE(24,23)PARTIC(JPART),XF,YF,ZF,TF,Z,A
!   23 FORMAT(5X,A7,' IN  SCATTERING  POINT:',2(F7.2,','),F7.2,
!     *' HAVE ENERGY ',F8.2,' MEV AND SCATTAR ON NUCLEUS Z=',
!     *F4.0,' A=',F6.2)
!C
!      WRITE(24,24)PARTIC(JPART),XZ,YZ,ZZ,TZ,STATUS(ISTATE),NBRNEW
!   24 FORMAT(' SCATTERED ',A7,'  ACHIEVE  POINT:',2(F7.2,','),F7.2,
!     *'   AT ENERGY ',F8.2,' MEV AND ',A11,', (NBRAN=',I4,').'/)
!      RETURN
!C
!C-------------------- PRINTING OF INELASTIC STAR ----------------------
!   30 JPART=ITREE(1,JBRST)
!      NBRAN=ITREE(3,JBRST)
!	  XI=TREE(1,JBRST)
!	  YI=TREE(2,JBRST)
!	  ZI=TREE(3,JBRST)
!	  TI=TREE(7,JBRST)
!            IAPR=INT(TREE(14,JBRST))                                 !HION6
!            IZPR=INT(TREE(15,JBRST))                                 !HION6
!	      XF=TREE(1,JBRST+1)
!	      YF=TREE(2,JBRST+1)
!	      ZF=TREE(3,JBRST+1)
!	      TF=TREE(7,JBRST+1)
!		  Z=TREE(5,JBRST+2)
!		  A=TREE(4,JBRST+2)
!C
!      if(JPART.eq.25)then                                            !HION6
!        WRITE(24,130)PARTIC(JPART),IAPR,IZPR,NBRAN                   !HION6
!  130   FORMAT(' INELASTIC INTERACTION OF ',A7,'(',I3,',',I2,')',    !HION6
!     *         ', MEV==>MeV/A, (NBRAN=',I4,'):')                     !HION6
!      else                                                           !HION6
!        WRITE(24,31)PARTIC(JPART),NBRAN
!   31   FORMAT(' INELASTIC INTERACTION OF ',A7,', (NBRAN=',I4,'):')
!      end if                                                         !HION6
!C
!      WRITE(24,32)PARTIC(JPART),XI,YI,ZI,TI,Z,A
!   32 FORMAT(5X,A7,' BEFORE RANGE IN POINT:',2(F7.2,','),F7.2,
!     *' HAVE ENERGY ',F8.2,' MEV. NUCLEUS-TARGET IS   Z=',F4.0,
!     *' A=',F4.0)
!C
!      WRITE(24,33)PARTIC(JPART),XF,YF,ZF,TF
!   33 FORMAT(5X,A7,' IN  INTERACTION POINT:',2(F7.2,','),F7.2,
!     *' HAVE ENERGY ',F8.2,' MEV. NUCLEI-PRODUCTS ARE:   Z*',
!     *3X,'  A*     U*')
!C
!C Printout of residual nuclei
!      IF(ITREE(7,JBRST+3).EQ.6)THEN
!	  WRITE(24,331)
!  331     FORMAT(107X,'FULL DESINTEGRATION')
!      ELSE
!	  DO JN=JBRST+3,JBRFN
!	      WRITE(24,332) TREE(5,JN),TREE(4,JN),TREE(6,JN)
!  332         FORMAT(107X,F4.0,3X,F4.0,F7.2)
!	  END DO
!      END IF
!C
!C Printout of secondary particles
!      WRITE(24,34)
!   34 FORMAT(' SECONDARIES: NO.  PARTICLE    ENERGY',10X,'DIRECTION',
!     *16X,'ACHIEVE POINT',11X,'STATUS    NBRAN')
!C
!      DO 40 I=1,NPCAS
!	  XZ=TREVS(1,I)
!	  YZ=TREVS(2,I)
!	  ZZ=TREVS(3,I)
!	  TZ=TREVS(7,I)
!	      COST=TREVS(4,I)
!	      SINF=TREVS(5,I)
!	      COSF=TREVS(6,I)
!		  ISTATE=ITREVS(7,I)
!		  IBRAN= ITREVS(3,I)
!		  IPART= ITREVS(1,I)
!                      IAPR=TREVS(14,I)                               !HION6
!                      IZPR=TREVS(15,I)                               !HION6
!C
!      IF(IPART.EQ.5)THEN                                             !HION6
!        WRITE(24,41)I,PARTIC(IPART),TZ,COST,SINF,COSF,
!     *                STATUS(ISTATE),IBRAN
!   41   FORMAT(13X,I4,3X,A7,F10.2,'  (',2(F7.4,','),F7.4,')',29X,A11,I7)
!      END IF                                                         !HION6
!C
!      if(IPART.eq.25)then                                            !HION6
!        WRITE(24,140)I,IAPR,IZPR,TZ,COST,SINF,COSF,XZ,YZ,ZZ,         !HION6
!     *                 STATUS(ISTATE),IBRAN                          !HION6
!  140   FORMAT(13X,I4,'HI(',I3,',',I2,')',F10.2,'  (',2(F7.4,','),F7.4, !HION6
!     *  ')  (',2(F7.2,','),F7.2,')  ',A11,I7)                        !HION6
!      end if                                                         !HION6
!C
!      if(IPART.ne.5 .AND. IPART.ne.25)then                           !HION6
!        WRITE(24,42)I,PARTIC(IPART),TZ,COST,SINF,COSF,XZ,YZ,ZZ,
!     *                 STATUS(ISTATE),IBRAN
!   42   FORMAT(13X,I4,3X,A7,F10.2,'  (',2(F7.4,','),F7.4,
!     *  ')  (',2(F7.2,','),F7.2,')  ',A11,I7)
!      end if                                                         !HION6
!C
!   40 CONTINUE
!C
!      DO 50 I=NPCAS+1,NPTOT
!	  TF=TREVS(7,I)
!	      COST=TREVS(4,I)
!	      SINF=TREVS(5,I)
!	      COSF=TREVS(6,I)
!      WRITE(24,51)I,PARTIC(1),TF,COST,SINF,COSF,STATUS(8)
!   51 FORMAT(13X,I4,3X,A7,F10.2,'  (',2(F7.4,','),F7.4,')',29X,A11)
!   50 CONTINUE
!C
!      WRITE(24,52)
!   52 FORMAT()
!      RETURN
!C
!C-------------------- PRINTING OF DECAY STAR ----------------------
!  200 JPART=ITREE(1,JBRST)
!      NBRAN=ITREE(3,JBRST)
!	  XI=TREE(1,JBRST)
!	  YI=TREE(2,JBRST)
!	  ZI=TREE(3,JBRST)
!	  TI=TREE(7,JBRST)
!	      XF=TREE(1,JBRST+1)
!	      YF=TREE(2,JBRST+1)
!	      ZF=TREE(3,JBRST+1)
!	      TF=TREE(7,JBRST+1)
!C
!      WRITE(24,231)PARTIC(JPART),NBRAN
!  231 FORMAT(' DECAY OF ',A7,', (NBRAN=',I4,'):')
!C
!      WRITE(24,232)PARTIC(JPART),XI,YI,ZI,TI
!  232 FORMAT(5X,A7,' BEFORE RANGE IN POINT:',2(F7.2,','),F7.2,
!     *' HAVE ENERGY ',F8.2,' MEV.')
!C
!      WRITE(24,233)PARTIC(JPART),XF,YF,ZF,TF
!  233 FORMAT(5X,A7,' IN DECAY POINT:',7X,2(F7.2,','),F7.2,
!     *' HAVE ENERGY ',F8.2,' MEV.')
!C
!      WRITE(24,234)
!  234 FORMAT(' PRODUCTS:    NO.  PARTICLE    ENERGY',10X,'DIRECTION',
!     *16X,'ACHIEVE POINT',11X,'STATUS    NBRAN')
!C
!      DO 240 I=1,NPCAS
!	  XZ=TREVS(1,I)
!	  YZ=TREVS(2,I)
!	  ZZ=TREVS(3,I)
!	  TZ=TREVS(7,I)
!	      COST=TREVS(4,I)
!	      SINF=TREVS(5,I)
!	      COSF=TREVS(6,I)
!		  ISTATE=ITREVS(7,I)
!		  IBRAN= ITREVS(3,I)
!		  IPART= ITREVS(1,I)
!      WRITE(24,242)I,PARTIC(IPART),TZ,COST,SINF,COSF,XZ,YZ,ZZ,
!     *             STATUS(ISTATE),IBRAN
!  242 FORMAT(13X,I4,3X,A7,F10.2,'  (',2(F7.4,','),F7.4,
!     *')  (',2(F7.2,','),F7.2,')  ',A11,I7)
!  240 CONTINUE
!C
!      DO 250 I=NPCAS+1,NPTOT
!	  TF=TREVS(7,I)
!	      COST=TREVS(4,I)
!	      SINF=TREVS(5,I)
!	      COSF=TREVS(6,I)
!		  IPART=ITREVS(1,I)
!      WRITE(24,251)I,PARTIC(IPART),TF,COST,SINF,COSF,STATUS(8)
!  251 FORMAT(13X,I4,3X,A7,F10.2,'  (',2(F7.4,','),F7.4,')',29X,A11)
!  250 CONTINUE
!
!C
!      WRITE(24,252)
!  252 FORMAT()
      RETURN
C
      END
 
 
 
      SUBROUTINE OUTTOT(LOWENT,TMAX)      
C REGISTRATION OF TOTAL CHARACTERISTICS OF SHIELD/CG EXTRANUCLEAR CASCADE.
C INPUT:
C      LOWENT  =  DON'T USED
C      TMAX - MAXIMAL ENERGY (MEV) FOR HISTOGRAMMING OF OUTGOING PARTICLES.
C OUTPUT:
C       SPNHIG - SPECTRUM OF OUTGOING NEUTRONS WITH T > 14.5 MEV. ! BNAB28
C       SPROT  - SPECTRUM OF OUTGOING PROTONS.
C       SPIM   - SPECTRUM OF OUTGOING PI-.
C       SPIP   - SPECTRUM OF OUTGOING PI+.
C       SANTIN - SPECTRUM OF OUTGOING ANTINEUTRONS.
C       SANTIP - SPECTRUM OF OUTGOING ANTIPROTONS.
C       SKMIN  - SPECTRUM OF OUTGOING K-.
C       SKPLS  - SPECTRUM OF OUTGOING K+.
C       SK0L   - SPECTRUM OF OUTGOING K0L.
C       SK0S   - SPECTRUM OF OUTGOING K0S.
C       SMUM   - SPECTRUM OF OUTGOING MU-.
C       SMUP   - SPECTRUM OF OUTGOING MU+.
C
C       SOURCE - SPECTRUM OF NEUTRON SOURCE INSIDE TARGET.
C       SOGAM  - SPECTRUM OF GAMMA SOURCE INSIDE TARGET.
C       SOEL   - SPECTRUM OF ELECTRON SOURCE INSIDE TARGET.
C       SOPOS  - SPECTRUM OF POSITRON SOURCE INSIDE TARGET.
C       SONUEL - SPECTRUM OF ELECTRONIC NEUTRINO SOURCE INSIDE TARGET.
C       SOANUE - SPECTRUM OF ELECTRONIC ANTINEUTRINO SOURCE INSIDE TARGET.
C       SONUMU - SPECTRUM OF MUONIC NEUTRINO SOURCE INSIDE TARGET.
C       SOANUM - SPECTRUM OF MUONIC ANTINEUTRINO SOURCE INSIDE TARGET.
C
C       TARGIN(I,J), J - NUMBER OF INPUT ZONE
C             (1,J) - NUMBER OF INPUT ZONE
C             (2,J) - NUMBER OF MEDIUM OF INPUT ZONE J
C             (3,J) - IONIZATION LOSS (MEV) IN INPUT ZONE J
C             (4,J) - Recoil nuclei energy (MeV) in input zone J
C             (5,J) - RESIDUAL EXCITATION ENERGY IN INPUT ZONE J
C             (6,J) - NO. OF PI0 IN INPUT ZONE J
C             (7,J) - SUMMAR KINETIC ENERGY OF PI0 IN INPUT ZONE J
C             (8,J) - NO. OF PI+, ABSORBED IN INPUT ZONE J
C             (9,J) - NO. OF K+, ABSORBED IN INPUT ZONE J
C            (10,J) - NO. OF MU-, ABSORBED IN INPUT ZONE J
C            (11,J) - NO. OF MU+, ABSORBED IN INPUT ZONE J
C            (12,J) - NO. OF PROTONS, ABSORBED IN INPUT ZONE J
C            (13,J) - NO. OF DEUTRONS, ABSORBED IN INPUT ZONE J
C            (14,J) - NO. OF TRITONS, ABSORBED IN INPUT ZONE J
C            (15,J) - NO. OF HELIUM3, ABSORBED IN INPUT ZONE J
C            (16,J) - NO. OF ALPHAS, ABSORBED IN INPUT ZONE J
C            (17,J) - No. of elastic interactions in input zone J
C            (18,J) - No. of inelastic interactions in input zone J
C                            (including fissions)
C            (19,J) - No. of fissions in input zone J
C--------------------------------------------------------------------------
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      COMMON /TREE/ TREE(15,400000),ITREE(7,400000),LTR13,LTR7,LTR300 ! GSITR
C
      COMMON /SOURN/ SOURN(13,100000),ISOURN(7,100000),LSN13,LSN7,LSN100   ! GSITR
      COMMON /SOURG/ SOURG(13,20000),ISOURG(7,20000),LSG13,LSG7,LSG100  ! GSITR
      COMMON /SEMIN/ SEMIN(13,100),ISEMIN(7,100),L1313,L137,L13100
      COMMON /SEPLS/ SEPLS(13,100),ISEPLS(7,100),L1413,L147,L14100
      COMMON /SNUEL/ SNUEL(13,100),ISNUEL(7,100),L1713,L177,L17100
      COMMON /SANUE/ SANUE(13,100),ISANUE(7,100),L1813,L187,L18100
      COMMON /SNUMU/ SNUMU(13,100),ISNUMU(7,100),L1913,L197,L19100
      COMMON /SANUM/ SANUM(13,100),ISANUM(7,100),L2013,L207,L20100
C
      COMMON /TRAN/ TRAN(5,600000),ITRAN(6,600000),LSR5,LSR6,LSR500  ! GSITR
      COMMON /TOTOUT/ SPNHIG(28),SPROT(28),SPIM(28),SPIP(28),
     *                SANTIN(28),SANTIP(28),SKMIN(28),SKPLS(28),
     *                SK0L(28),SK0S(28),SMUM(28),SMUP(28),
     *                SOURCE(-1:28),SOGAM(28),SOEL(28),SOPOS(28), ! BNAB28
     *                SONUEL(28),SOANUE(28),SONUMU(28),SOANUM(28),
     *                TARGIN(19,1000),LIZ50
      COMMON /AZDIST/ AZ(25,100),IOUT,IIN
      COMMON /THRTAR/ ETHR,SETHR,WTHR,NTHR    ! Throughout target 14.04.97
C
      DIMENSION BNABEN(-1:27)                                     ! BNAB28
      DATA BNABEN /14.5, 14.0,                                    ! BNAB28
     *10.5,  6.5, 4.0, 2.5, 1.4, 0.8, 0.4, 0.2, 0.1,
     *4.65E-02, 2.15E-02, 1.00E-02, 4.65E-03, 2.15E-03, 1.00E-03,
     *4.65E-04, 2.15E-04, 1.00E-04, 4.65E-05, 2.15E-05, 1.00E-05,
     *4.65E-06, 2.15E-06, 1.00E-06, 4.65E-07, 2.15E-07, 0./
C
C-----------------------------------------------------------------------
CC
C      CONTINUE!WRITE(25,1067)                                            !HIONDEB
C 1067 FORMAT(10X,'TREE')                                        !HIONDEB
C      DO K=1,100                                                !HIONDEB
C        CONTINUE!WRITE(25,1068)K,(TREE(J,K),J=1,8),(ITREE(J,K),J=1,7)    !HIONDEB
C 1068   FORMAT(I2,1X,3F6.2,3F7.3,F7.1,F6.3,7I4)                 !HIONDEB
C      END DO                                                    !HIONDEB
CC
C Number of branches in tree
      NBRMAX=ITREE(4,1)
C Number of lines in tree (excluding first and last lines)
      NLINE =ITREE(3,1)
C Numbers of particles in sources
	  NNEUTR=ISOURN(LSN7,1)
	  NGAMMA=ISOURG(LSG7,1)
	  NEMIN =ISEMIN(L137,1)
	  NEPLS =ISEPLS(L147,1)
	  NNUEL =ISNUEL(L177,1)
	  NANUE =ISANUE(L187,1)
	  NNUMU =ISNUMU(L197,1)
	  NANUM =ISANUM(L207,1)
C Registration of primary particles, which passed through      ! 14.04.97
C the configuration without elastic or inelastic interaction.
      KSTATE=ITREE(7,3)
      IF(KSTATE.eq.3)THEN
          NTHR=NTHR+1                       ! Sum N     !HION7:TREE(7,3)*A?
          WTHR=WTHR+TREE(8,3)               ! Sum W
          ETHR=TREE(7,3)*TREE(8,3)          !     E*W
          SETHR=SETHR+TREE(7,3)*TREE(8,3)   ! Sum E*W
      END IF
C
C HISTOGRAMMING OF LOW ENERGY NEUTRON SOURCE.
      DO 20 J=2,NNEUTR+1
      TPAR  =SOURN(7,J)
      WEIPAR=SOURN(8,J)
         DO 21 I=0,27                                             ! BNAB28
	     IF(TPAR.LT.BNABEN(I))GO TO 21
	     N1=I-1
	     GO TO 22
   21     CONTINUE
   22         SOURCE(N1)=SOURCE(N1)+WEIPAR
	      SOURCE(27)=SOURCE(27)+WEIPAR
	      SOURCE(28)=SOURCE(28)+TPAR*WEIPAR
   20 CONTINUE
C
C Step STEPS for sources energy distributions.
      STEPS=0.05*TMAX
C Energy distributions of sources.
C---------------------------------- Gamma-quant source
      DO 200 J=2,NGAMMA+1
	 TPAR=  SOURG(7,J)
	 WEIPAR=SOURG(8,J)
  200 CALL HIST(TPAR,0.,TMAX,STEPS,SOGAM,28,WEIPAR)
C---------------------------------- Electron source
      DO 210 J=2,NEMIN+1
	 TPAR=  SEMIN(7,J)
	 WEIPAR=SEMIN(8,J)
  210 CALL HIST(TPAR,0.,TMAX,STEPS,SOEL,28,WEIPAR)
C---------------------------------- Positron source
      DO 220 J=2,NEPLS+1
	 TPAR=  SEPLS(7,J)
	 WEIPAR=SEPLS(8,J)
  220 CALL HIST(TPAR,0.,TMAX,STEPS,SOPOS,28,WEIPAR)
C---------------------------------- Electronic neitrino source
      DO 230 J=2,NNUEL+1
	 TPAR=  SNUEL(7,J)
	 WEIPAR=SNUEL(8,J)
  230 CALL HIST(TPAR,0.,TMAX,STEPS,SONUEL,28,WEIPAR)
C---------------------------------- Electronic antineitrino source
      DO 240 J=2,NANUE+1
	 TPAR=  SANUE(7,J)
	 WEIPAR=SANUE(8,J)
  240 CALL HIST(TPAR,0.,TMAX,STEPS,SOANUE,28,WEIPAR)
C---------------------------------- Muonic neitrino source
      DO 250 J=2,NNUMU+1
	 TPAR=  SNUMU(7,J)
	 WEIPAR=SNUMU(8,J)
  250 CALL HIST(TPAR,0.,TMAX,STEPS,SONUMU,28,WEIPAR)
C---------------------------------- Muonic antineitrino source
      DO 260 J=2,NANUM+1
	 TPAR=  SANUM(7,J)
	 WEIPAR=SANUM(8,J)
  260 CALL HIST(TPAR,0.,TMAX,STEPS,SOANUM,28,WEIPAR)
C
C Spectra of flying out particles.
	  STEP=0.05*TMAX
      DO 30 J=2,NLINE+1
      JPART =ITREE(1,J)
      KSTATE=ITREE(7,J)
      IF(JPART.NE.0 .AND. KSTATE.EQ.3)THEN
            IF(JPART.LT.1 .OR. JPART.GT.25)STOP 'OUTTOT1'            !HION6
	  TPAR=  TREE(7,J)
	  WEIPAR=TREE(8,J)
	  GO TO(301,302,303,304, 30,306,307,308,309,
     *          310,311,312,313,314,315,316,                         !HION6
     *           30, 30, 30, 30,321,322,323,324,325),JPART           !HION6
  301     CALL HIST(TPAR,0.,TMAX,STEP,SPNHIG,28,WEIPAR)    ! Neutron
	  GO TO 30
  302     CALL HIST(TPAR,0.,TMAX,STEP,SPROT ,28,WEIPAR)    ! Proton
	  GO TO 30
  303     CALL HIST(TPAR,0.,TMAX,STEP,SPIM  ,28,WEIPAR)    ! PI-
	  GO TO 30
  304     CALL HIST(TPAR,0.,TMAX,STEP,SPIP  ,28,WEIPAR)    ! PI+
	  GO TO 30
  306     CALL HIST(TPAR,0.,TMAX,STEP,SANTIN,28,WEIPAR)    ! Antineutron
	  GO TO 30
  307     CALL HIST(TPAR,0.,TMAX,STEP,SANTIP,28,WEIPAR)    ! Antiproton
	  GO TO 30
  308     CALL HIST(TPAR,0.,TMAX,STEP,SKMIN ,28,WEIPAR)    ! K-
	  GO TO 30
  309     CALL HIST(TPAR,0.,TMAX,STEP,SKPLS ,28,WEIPAR)    ! K+
	  GO TO 30
  310         K0L=IFIX(TREE(10,J))
	      K0S=IFIX(TREE(11,J))  
		  IF((K0L.EQ.0 .AND. K0S.EQ.0) .OR. 
     *               (K0L.EQ.1 .AND. K0S.EQ.1))STOP 'OUTTOT2'
              IF(K0L.EQ.1)CALL HIST(TPAR,0.,TMAX,STEP,SK0L,28,WEIPAR)  ! K0L
              IF(K0S.EQ.1)CALL HIST(TPAR,0.,TMAX,STEP,SK0S,28,WEIPAR)  ! K0S
	  GO TO 30
  311     GO TO 310
  312     CONTINUE   ! Reserv: Gamma
          GO TO 30 
  313     CONTINUE   ! Reserv: Electron
          GO TO 30 
  314     CONTINUE   ! Reserv: Positron
          GO TO 30 
  315     CALL HIST(TPAR,0.,TMAX,STEP,SMUM  ,28,WEIPAR)   ! MU-
	  GO TO 30
  316     CALL HIST(TPAR,0.,TMAX,STEP,SMUP  ,28,WEIPAR)   ! MU+
	  GO TO 30
  321     CONTINUE   ! Reserv: Deuteron
          GO TO 30 
  322     CONTINUE   ! Reserv: Tritium
          GO TO 30 
  323     CONTINUE   ! Reserv: He3
          GO TO 30 
  324     CONTINUE   ! Reserv: Alpha
          GO TO 30 
  325     CONTINUE   ! Reserv: Heavy Ion
          GO TO 30 
      END IF
   30 CONTINUE
C
C----------------------- INTRINSIC CHARACTERISTICS ---------------------
C
	  JI=2
	  LI=2
                          DO 100 JBR=1,NBRMAX  !******Main loop
C
C THE FIRST - JI,LI AND THE LAST - JF,LF LINES OF FRAGMENTS OF ARRAYS
C TREE AND TRAN, CONTAINING THE BRANCH NO.NBR
      JPART =ITREE(1,JI)    ! Type of particle which initiates branch No.JBR
      APART=1.0         ! multiplier for Energy: Etot=E/A * APART  !HION7
      IF(JPART.eq.21)APART=2.0                                     !HION7
      IF(JPART.eq.22)APART=3.0                                     !HION7
      IF(JPART.eq.23)APART=3.0                                     !HION7
      IF(JPART.eq.24)APART=4.0                                     !HION7
      IF(JPART.eq.25)APART=TREE(14,JI)                             !HION7
      KSTATE=ITREE(7,JI+1)
C
      IF(KSTATE.LE.2)THEN
C       Counting of residual nuclei number NUC 
        NUC=0
	  DO JN=JI+3,NLINE+1
          IF(ITREE(1,JN).EQ.0)THEN
            NUC=NUC+1
          ELSE
            GO TO 10
          END IF
	  END DO
   10   CONTINUE
C Debug check of residual nuclei number NUC
        IF(NUC.EQ.0)THEN
        CONTINUE!WRITE(25,11)JI,JF
   11   FORMAT(' OUTTOT DIAGNOSTIC: THERE IS ONLY NUCLEUS-TARGET',
     *         ' IN BRANCH,'/20X,'JI=',I4,',   JF=',I4)
	  CLOSE(25)
	  STOP
        END IF
C End debug check NUC ---------------------
        JF=JI+2+NUC
      ELSE
        JF=JI+1
      END IF
C
      LF=LI+ITRAN(5,LI)
C
C IF BRANCH NO.JBR CONTAIN INTERACTION - REGISTRATE THE KINETIC AND
C EXCITATION ENERGY OF RESIDUAL NUCLEUS.
      IF(KSTATE.le.2)THEN
	  DO JN=JI+3,JF        
	    PNX=TREE( 8,JN)
	    PNY=TREE( 9,JN)
	    PNZ=TREE(10,JN)
            RMAS  =940.*TREE(4,JN)
          WEINUC=TREE(7,JN)
C
            IF(RMAS.LE.0.)THEN
C             FULL DESINTEGRATION OF TARGET NUCLEUS
C             OR TARGET NUCLEUS IS HIDROGEN (12.04.90)
              TNUC=0.
              ESTAR=0.
            ELSE                   !HION6
C-Nonrel.      TNUC=(PNX**2+PNY**2+PNZ**2)/(RMAS*2.0)
              TNUC=SQRT((PNX**2+PNY**2+PNZ**2)+RMAS**2)-RMAS
              ESTAR=TREE(6,JN)
            END IF
C
          INZON=ITRAN(6,LF)
	    TARGIN(4,INZON)=TARGIN(4,INZON)+TNUC *WEINUC
	    TARGIN(5,INZON)=TARGIN(5,INZON)+ESTAR*WEINUC
C
            IF(KSTATE.EQ.2)THEN
C             (A,Z)-distribution
              ARES=TREE(4,JN)
              ZRES=TREE(5,JN)
              CALL RESNHS(ARES,ZRES,WEINUC,AZ,IOUT,IIN,INZON)  ! DCHAIN
            END IF
	  END DO
C
C       Number of inelastic interactions (including fissions)
C       and number of fissions
        WEINUC=TREE(7,JF)
        IF(KSTATE.EQ.2)THEN
          TARGIN(18,INZON)=TARGIN(18,INZON)+WEINUC
          KINTER=ITREE(7,JF)
            IF(KINTER.EQ.3)THEN
              TARGIN(19,INZON)=TARGIN(19,INZON)+WEINUC
            END IF
	  END IF
C         Number of elastic interactions
	  IF(KSTATE.EQ.1)THEN
	      TARGIN(17,INZON)=TARGIN(17,INZON)+WEINUC
	  END IF
      END IF
C
C REGISTRATION OF PI0 DECAY IN FLY.
      IF(JPART.EQ.5)THEN
        TPAR  =TREE(7,JI)
        WEIPAR=TREE(8,JI)
        INZON=ITRAN(6,LF)
        TARGIN(6,INZON)=TARGIN(6,INZON)+WEIPAR
        TARGIN(7,INZON)=TARGIN(7,INZON)+TPAR*WEIPAR
        GO TO 99
      END IF
C
C REGISTRATION OF IONIZATION LOSS OF PROTONS, PI-, PI+, antiprotons,K-, 
C K+, MU-, MU+, DEUTERONS, TRITONS, He3, ALPHAS, HEAVY IONS          !HION6
C AND ACTS OF ABSORPTION OF PI+, K+, MU-, MU+, PROTONS, DEUTERONS,   !HION6
C TRITONS, He3, ALPHAS, and HEAVY IONS.                              !HION6
C
      IF(JPART.EQ. 2 .OR. JPART.EQ. 3 .OR. JPART.EQ.4 .OR.
     *   JPART.EQ. 7 .OR. JPART.EQ. 8 .OR. JPART.EQ.9 .OR.
     *   JPART.EQ.15 .OR. JPART.EQ.16 .OR.                       !HION7
     *   JPART.eq.21 .OR. JPART.eq.22 .OR. JPART.eq.23 .OR.      !HION7
     *                        JPART.eq.24 .OR. JPART.eq.25)THEN  !HION7
C
            WEIPAR=TREE(8,JI)
            DO I=LI+1,LF
              INZON=ITRAN(6,I)
              QZON=TRAN(5,I-1)-TRAN(5,I)
              TARGIN(3,INZON)=TARGIN(3,INZON)+QZON*WEIPAR*APART      !HION7
            END DO
C
            IF(KSTATE.EQ.4 .OR. KSTATE.EQ.7)THEN
              INZON=ITRAN(6,LF)                                      !HION7?
              TARGIN(3,INZON)=TARGIN(3,INZON)+TRAN(5,LF)*WEIPAR*APART !HION7
              IF(JPART.EQ.2)TARGIN(12,INZON)=TARGIN(12,INZON)+WEIPAR
              IF(JPART.EQ.4)TARGIN(8,INZON)=TARGIN(8,INZON)+WEIPAR
              IF(JPART.EQ.9)TARGIN(9,INZON)=TARGIN(9,INZON)+WEIPAR
              IF(JPART.EQ.15)TARGIN(10,INZON)=TARGIN(10,INZON)+WEIPAR
              IF(JPART.EQ.16)TARGIN(11,INZON)=TARGIN(11,INZON)+WEIPAR
              IF(JPART.EQ.21)TARGIN(13,INZON)=TARGIN(13,INZON)+WEIPAR !HION7
C              IF(JPART.EQ.22)TARGIN(14,INZON)=TARGIN(14,INZON)+WEIPAR !HION7
C--------------------------------DCHAIN----------------------
C Scoring of tritium to NMTC_YIELD also
              IF(JPART.EQ.22)THEN
                TARGIN(14,INZON)=TARGIN(14,INZON)+WEIPAR !HION7
                AHION=3. ! tritium
                ZHION=1. ! tritium
                CALL RESNHS(AHION,ZHION,WEIPAR,AZ,IOUT,IIN,INZON)  ! DCHAIN
              END IF
C----------------------------end DCHAIN----------------------
              IF(JPART.EQ.23)TARGIN(15,INZON)=TARGIN(15,INZON)+WEIPAR !HION7
              IF(JPART.EQ.24)TARGIN(16,INZON)=TARGIN(16,INZON)+WEIPAR !HION7
C-reserv      IF(JPART.EQ.25)TARGIN(**,INZON)=TARGIN(**,INZON)+WEIPAR !HION7
C======================= KarSjh ==========================================
C             To add the absorbed HI to (A,Z)-distribution over target.
              IF(JPART.EQ.25)THEN
                AHION=TREE(14,JF)  ! A of the absorbed heavy ion
                ZHION=TREE(15,JF)  ! Z of the absorbed heavy ion
c                CONTINUE!WRITE(25,1234)JPART,AHION,ZHION
c 1234           format(' OUTPUT###JPART,AHION,ZHION',I3,2F6.1)
                CALL RESNHS(AHION,ZHION,WEIPAR,AZ,IOUT,IIN,INZON)  ! DCHAIN
              END IF
C===================== end KarSjh =========================================
            END IF
      END IF
C
C
   99 JI=JF+1
      LI=LF+1
  100                     CONTINUE  !******end of Main loop
C
      RETURN
      END



      SUBROUTINE HIST(X,A,B,H,RX,N,W)                         
      DIMENSION RX(N)
      LE=(B-A)/H+8
      IF(LE-N)11,11,12
   12 CONTINUE!WRITE(25,13)LE,N   
   13 FORMAT(5X,'ERROR IN DIMENSION IN HIST   L=',I3,'N=',I3)
      RETURN
C
   11 CONTINUE   
      RX(N)=RX(N)+X*W
      RX(N-1)=RX(N-1)+W
C  
      IF(X-A)1,2,2 
    1 RX(N-6)=RX(N-6)+X*W
      RX(N-7)=RX(N-7)+W
      RETURN
C
    2 IF(X-B)4,3,3
    3 RX(N-2)=RX(N-2)+X*W
      RX(N-3)=RX(N-3)+W
      RETURN
C
    4 L=(X-A)/H  
      RX(L+1)=RX(L+1)+W
      RX(N-4)=RX(N-4)+X*W
      RX(N-5)=RX(N-5)+W
      RETURN
      END
 
 
 
      SUBROUTINE CLETOT
C CLEANING OF ARRAYS FOR TOTAL OUTPUT.
      COMMON /TOTOUT/ SPNHIG(28),SPROT(28),SPIM(28),SPIP(28),
     *                SANTIN(28),SANTIP(28),SKMIN(28),SKPLS(28),
     *                SK0L(28),SK0S(28),SMUM(28),SMUP(28),
     *                SOURCE(-1:28),SOGAM(28),SOEL(28),SOPOS(28), ! BNAB28
     *                SONUEL(28),SOANUE(28),SONUMU(28),SOANUM(28),
     *                TARGIN(19,1000),LIZ50
      COMMON /AZDIST/ AZ(25,100),IOUT,IIN
      COMMON /AZSPD/ AZSPD(280,100,0:1000),LZS50        ! DCHAIN
      COMMON /ISOSET/ WISOSET(2000),FERTAZ(2000,50),  ! Fertman ! Truba06
     *                ISOSET(2,2000),TDEC(2000),NIS200,NISTRU   ! userset
      CHARACTER*10 TDEC   ! userset
      COMMON /THRTAR/ ETHR,SETHR,WTHR,NTHR    ! Throughout target 14.04.97
C
      COMMON /TLEPA/ TLEPA(122,20,1000),NG73,LTLEPA  ! TLEST
      REAL*8 TLEPA                                  ! TLEST
C
C
        ETHR=0.   ! 14.04.97
        SETHR=0.  ! 14.04.97
        WTHR=0.   ! 14.04.97
        NTHR=0    ! 14.04.97
      DO 1 J=1,28
	  SPNHIG(J)=0.
	  SPROT(J) =0.
	  SPIM(J)  =0.
	  SPIP(J)  =0.
	  SANTIN(J)=0.
	  SANTIP(J)=0.
	  SKMIN(J) =0. 
	  SKPLS(J) =0.
	  SK0L(J)  =0.
	  SK0S(J)  =0.
	  SMUM(J)  =0.
	  SMUP(J)  =0.
          SOURCE(J)=0.
                        SOURCE(-1)=0.        ! BNAB28
                        SOURCE(0)=0.         ! BNAB28
	  SOGAM(J) =0.
	  SOEL(J)  =0.
	  SOPOS(J) =0.
	  SONUEL(J)=0.
	  SOANUE(J)=0.
	  SONUMU(J)=0.
	  SOANUM(J)=0.
    1 CONTINUE
C======= TLEST ===============
      DO K=1,LTLEPA
        DO J=1,20
          DO I=1,NG73
             TLEPA(I,J,K)=0.0D00
          END DO
        END DO
      END DO
C======= end TLEST ===========
      DO 2 J=1,LIZ50 
	  DO 2 I=1,19
    2     TARGIN(I,J)=0.
C
	  IOUT=0
	  IIN=0
      DO 3 I=1,100
	  DO 3 J=1,25
    3     AZ(J,I)=0.
C
      DO K=0,LZS50
        DO I=1,100
          DO J=1,280
            AZSPD(J,I,K)=0.        ! DCHAIN
          END DO
        END DO
      END DO
C
      DO J=1,NIS200
        WISOSET(J)=0.       ! userset
        DO K=1,50
          FERTAZ(J,K)=0.    ! Fertman
        END DO
      END DO
C Preparation action fo (A,Z)-histogramming
      CALL RESNCL
C
C Cleaning of arrays for AGT to SHIELD conversion protocol
      CALL CLECON
C
      RETURN
      END
 
 
 
      SUBROUTINE PRNTOT(TMAX)                                        ! GEMCA
C PRINTING OF TOTAL SHIELD/CG OUTPUT.
C---------------------------COMMONS GEOMETRY                         ! GEMCA
	COMMON /GDATA0/ NUMBOD,NUMZON                                ! GEMCA
C---------------------------COMMONS OF SHIELD/CG
      COMMON /MANALL/ NTREE,LINE1,LINE2,NBRMAX,NPRMAX,NNEUTR,NGAMMA,
     *                MELSC,MINTER,MCUT,MOUT,MABSOR,MDECAY,NSTRAN,
     *                NEMIN,NEPLS,NNUEL,NANUE,NNUMU,NANUM
      COMMON /TOTOUT/ SPNHIG(28),SPROT(28),SPIM(28),SPIP(28),
     *                SANTIN(28),SANTIP(28),SKMIN(28),SKPLS(28),
     *                SK0L(28),SK0S(28),SMUM(28),SMUP(28),
     *                SOURCE(-1:28),SOGAM(28),SOEL(28),SOPOS(28), ! BNAB28
     *                SONUEL(28),SOANUE(28),SONUMU(28),SOANUM(28),
     *                TARGIN(19,1000),LIZ50       ! 21.02.97
      COMMON /AZDIST/ AZ(25,100),IOUT,IIN
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      COMMON /THRTAR/ ETHR,SETHR,WTHR,NTHR    ! Throughout target 14.04.97
C
      DO 1 J=1,NUMZON                                                ! GEMCA
      TARGIN(1,J)=FLOAT(J)
      TARGIN(2,J)=GZMED(J)                                           ! GEMCA
    1 CONTINUE
C
      WRITE(24,10)NTREE
   10 FORMAT(/10X,'TOTAL OUTPUT OF SHIELD/CG CODE. STATISTICS IS',
     *       I8,' EVENTS.'/)
C-------------------------------------------------------------------
      WRITE(24,2)NTHR,ETHR,WTHR,SETHR                          ! 14.04.97
    2 FORMAT(10X,'PRIMARY PARTICLES WITHOUT INTERACTION:'/
     *       12X,'Sum N=',I9/
     *       12X,'E*W  =',E20.13/
     *       14X,'Sum W  =',E20.13/
     *       14X,'Sum W*E=',E20.13/)
C
      WRITE(24,11)TMAX
   11 FORMAT(10X,'SPECTRA OF FLYING OUT PARTICLES (TMIN=0, TMAX=',
     *       F10.1,' MEV, 20 BINS).')
      WRITE(24,12)
   12 FORMAT(' NO.  NEUTRONS   PROTONS',7X,'PI-',7X,'PI+',5X,'ANTI-',5X,       
     *      'ANTI-',8X,'K-',8X,'K+',7X,'K0L',7X,'K0S',7X,'MU-',7X,'MU+')
      WRITE(24,13)
   13 FORMAT(' BIN    E>14.5',32X,'NEUTRONS   PROTONS')           ! BNAB28
C
      DO 100 J=1,20
  100 WRITE(24,14)J,SPNHIG(J),SPROT(J),SPIM(J),SPIP(J),
     *              SANTIN(J),SANTIP(J),SKMIN(J),SKPLS(J),
     *              SK0L(J),SK0S(J),SMUM(J),SMUP(J)
   14 FORMAT(I4,12F10.1)
C
      DO 110 J=21,28
      IF(J.EQ.21 .OR. J.EQ.23 .OR. J.EQ.25 .OR. J.EQ.27)THEN
	  ASSIGN 15 TO LABEL
      ELSE
	  ASSIGN 16 TO LABEL
      END IF
  110 WRITE(24,LABEL)SPNHIG(J),SPROT(J),SPIM(J),SPIP(J),
     *               SANTIN(J),SANTIP(J),SKMIN(J),SKPLS(J),
     *               SK0L(J),SK0S(J),SMUM(J),SMUP(J)
   15 FORMAT(4X,12F10.1)
   16 FORMAT(4X,12E10.3) 
C-------------------------------------------------------------
      WRITE(24,20)TMAX
   20 FORMAT(/5X,'SOUROCES OF PARTICLES INSIDE TARGET (FOR NEUTRONS-',
     *       'BNAB GRID, FOR OTHER PARTICLES: TMIN=0, TMAX=',F10.1,
     *        ' MEV, 20 BINS).') 
      WRITE(24,21)
   21 FORMAT('  NO.',5X,'NEUTRONS',13X,'NO.',5X,'GAMMA ELECTRONS',
     *       ' POSITRONS  NEITRINO  ANTINEI-  NEITRINO  ANTINEI-')
      WRITE(24,22)
   22 FORMAT(' GRUP   E<14.5 MEV',13X,'BIN',4X,'QUANTS',28X,       ! BNAB28
     *       'EL  TRINO EL',8X,'MU  TRINO MU')
C
      DO 199 J=-1,0                                               ! BNAB28
  199 WRITE(24,223)J,SOURCE(J)                                    ! BNAB28
  223 FORMAT(I5,F13.1)                                            ! BNAB28
C
      DO 200 J=1,20
  200 WRITE(24,23)J,SOURCE(J),J,SOGAM(J),SOEL(J),SOPOS(J),
     *              SONUEL(J),SOANUE(J),SONUMU(J),SOANUM(J)
   23 FORMAT(I5,F13.1,14X,I2,7F10.1)
C
      DO 210 J=21,26
      IF(J.EQ.21 .OR. J.EQ.23 .OR. J.EQ.25)THEN
	  ASSIGN 24 TO LABEL
      ELSE
	  ASSIGN 25 TO LABEL
      END IF
  210 WRITE(24,LABEL)J,SOURCE(J),SOGAM(J),SOEL(J),SOPOS(J),
     *                 SONUEL(J),SOANUE(J),SONUMU(J),SOANUM(J)
   24 FORMAT(I5,F13.1,16X,7F10.1)
   25 FORMAT(I5,F13.1,16X,7E10.3)
C
      WRITE(24,26)SOURCE(27),SOGAM(27),SOEL(27),SOPOS(27),
     *                 SONUEL(27),SOANUE(27),SONUMU(27),SOANUM(27)
   26 FORMAT(5X,F13.1,16X,7F10.1)
      WRITE(24,27)SOURCE(28),SOGAM(28),SOEL(28),SOPOS(28),
     *                 SONUEL(28),SOANUE(28),SONUMU(28),SOANUM(28)
   27 FORMAT(5X,E13.6,16X,7E10.3)
C----------------------------------------------------------------
      WRITE(24,30)
   30 FORMAT(/10X,'DISTRIBUTIONS INSIDE TARGET.')
      WRITE(24,31)
   31 FORMAT(' INP. NO.  IONIZATION     RECOIL     RESIDUAL',
     *       '     NO.  TSUM PI0',26X,
     *       'ABSORBED/CUT OFF PARTICLES',28X,'INTERACTIONS')
      WRITE(24,32)
   32 FORMAT(' ZONE MED  LOSS (MEV)  NUCLEI (MEV) EXCITATION',
     *       '    PI0   (MEV)',6X,
     *       '  PI+     K+    MU-    MU+     PROT',
     *       '    DEUT    TRIT     HE3    ALPH',
     *       '     ELSC    INEL    FISS')
C
C Cleaning of variables for summing over all zones.
	      TG3=0.
	      TG4=0.
	      TG5=0.
	  J6 =0
	      TG7=0.
	  J8 =0
	  J9 =0
	  J10=0
	  J11=0
	  J12=0
	  J13=0
	  J14=0
	  J15=0
	  J16=0
	  J17=0
	  J18=0
	  J19=0
C 
      DO 300 J=1,NUMZON                                              ! GEMCA
	  I1 =IFIX(TARGIN( 1,J))
	  I2 =IFIX(TARGIN( 2,J))
	  I6 =IFIX(TARGIN( 6,J))
	  I8 =IFIX(TARGIN( 8,J))
	  I9 =IFIX(TARGIN( 9,J))
	  I10=IFIX(TARGIN(10,J))
	  I11=IFIX(TARGIN(11,J))
	  I12=IFIX(TARGIN(12,J))
	  I13=IFIX(TARGIN(13,J))
	  I14=IFIX(TARGIN(14,J))
	  I15=IFIX(TARGIN(15,J))
	  I16=IFIX(TARGIN(16,J))
	  I17=IFIX(TARGIN(17,J))
	  I18=IFIX(TARGIN(18,J))
	  I19=IFIX(TARGIN(19,J))
C Summation over all zones
	      J6 = J6+I6
	      J8 = J8+I8
	      J9 = J9+I9
	      J10=J10+I10
	      J11=J11+I11
	      J12=J12+I12
	      J13=J13+I13
	      J14=J14+I14
	      J15=J15+I15
	      J16=J16+I16
	      J17=J17+I17
	      J18=J18+I18
	      J19=J19+I19
		  TG3=TG3+TARGIN(3,J)
		  TG4=TG4+TARGIN(4,J)
		  TG5=TG5+TARGIN(5,J)
		  TG7=TG7+TARGIN(7,J)
      WRITE(24,33)I1,I2,TARGIN(3,J),TARGIN(4,J),TARGIN(5,J),I6,
     *TARGIN(7,J),I8,I9,I10,I11,I12,I13,I14,I15,I16,I17,I18,I19
   33 FORMAT(I4,I4,2X,E12.6,1X,E12.6,1X,E10.4,I7,1X,E10.4,1X,4I7,I9,4I8,
     *       1X,3I8)
  300 CONTINUE
      WRITE(24,34)TG3,TG4,TG5,J6,TG7,J8,J9,J10,
     *            J11,J12,J13,J14,J15,J16,J17,J18,J19
   34 FORMAT('  SUMM:   ',
     *                E12.6,1X,E12.6,1X,E10.4,I7,1X,E10.4,1X,4I7,I9,4I8,
     *       1X,3I8)
C
C Printout of (A,Z)-distribution
      CALL RESNPR(AZ,IOUT,IIN)
      CALL PRNSPD                     ! DCHAIN
      CALL PRNSET                     ! userset
C
C Printout of AGT to SHIELD conversion protocol
      CALL PRICON
C
	  LASTIX=IX
      WRITE(24,40)IXFIRS,LASTIX
   40 FORMAT(/' RANDOM GENERATOR STATUS: IXFIRS=',I12/26X,'LASTIX=',I12)
C
      RETURN
      END
 
 
 
      SUBROUTINE RESNHS(A,Z,STATW,AZ,IOUT,IIN,INZON)
C Histogramming of residual nuclei (A,Z)-distribution
C INPUT: A,Z - A and Z of residual nucleus
C        STATW - Statistical weight of residual nucleus
C
C OUTPUT: AZ(25,100) - (A,Z)-distribution
C         IOUT, IIN - number of nuclei outside/inside (A,Z)-distribution
C         AZSPD(280,100,0:1000) - (A,Z,zone)-distribution as input for DCHAIN
C         /ISOSET/ Isotope scorinng over a whole target on user request:
C                  ISOSET(1,J) - Zj of nucleus product requested
C                  ISOSET(2,J) - Aj of nucleus product requested
C                  WISOSET(J) - Sum of weights of the nucleus (Zj,Aj)
C                  FERTAZ(J,Nzon) - Sum of weights in zone Nzon<50   ! Fertman
C-----------------------------------------------------------------------
      COMMON /RESN01/ ASTAB(100),IAA(25)
      COMMON /AZSPD/ AZSPD(280,100,0:1000),LZS50        ! DCHAIN
      COMMON /ISOSET/ WISOSET(2000),FERTAZ(2000,50),  ! Fertman ! Truba06
     *                ISOSET(2,2000),TDEC(2000),NIS200,NISTRU   ! userset
      CHARACTER*10 TDEC   ! userset
      DIMENSION AZ(25,100)
C
C HISTOGRAMMING OF RESIDUAL NUCLEI (A,Z)-DISTRIBUTION
	      IA=NINT(A)
	      IZ=NINT(Z)
        IF(IA.le.0 .OR. IZ.le.0)RETURN  ! Full desintegr. Checked 07.06.2004
	  IF(IZ.GT.100 .OR. IA.gt.280)THEN   ! DCHAIN
            CONTINUE!WRITE(25,1000) IZ,IA
 1000       FORMAT(' RESNHS DIAGNOSTIC: Z=',I3,', A=',I3,
     *             ' OUTSIDE (A,Z)-DISTRIBUTION')
	      RETURN
	  END IF
C
		  JA=IA-(IFIX(ASTAB(IZ))-12)+1
C
	  IF(JA.LT.1 .OR. JA.GT.25)THEN
	      IOUT=IOUT+1
	      RETURN
	  ELSE
	      IIN=IIN+1
	  END IF
C
		  AZ(JA,IZ)=AZ(JA,IZ)+STATW
            AZSPD(IA,IZ,0)=AZSPD(IA,IZ,0)+STATW        ! DCHAIN
C
            IF(INZON.gt.LZS50)THEN
              CONTINUE!WRITE(25,1001)INZON
 1001         FORMAT(' RESNHS DIAGNOSTIC: INZON=',I5)
              RETURN
            END IF
            AZSPD(IA,IZ,INZON)=AZSPD(IA,IZ,INZON)+STATW        ! DCHAIN
C
C//////////////////// userset //////////////////////
      DO JIS=1,NISTRU
        IF(IZ.eq.ISOSET(1,JIS) .AND. IA.eq.ISOSET(2,JIS))THEN
          WISOSET(JIS)=WISOSET(JIS)+STATW
          IF(INZON.ge.1 .AND. INZON.le.50)THEN          ! Fertman
            FERTAZ(JIS,INZON)=FERTAZ(JIS,INZON)+STATW   ! Fertman
          END IF                                        ! Fertman
c          GOTO 1          ! Fertman
        END IF
      END DO
c    1 CONTINUE            ! Fertman
C//////////////////// end userset //////////////////
      RETURN
      END



      SUBROUTINE RESNCL
C Auxiliary operations for residual nuclei histogramming 
      COMMON /RESN01/ ASTAB(100),IAA(25)
      COMMON /INFEL/ ZNUC(110),ATWEI(110),DENS(110),RIEV(110),SYMB(110)
      CHARACTER*2 SYMB
C
C MESH FOR RESIDUAL NUCLEUS (A,Z)-DISTRIBUTION
C     BETHA-STABILITY CURVE A0(Z):
      DO 10 J=1,100
	  Z=ZNUC(J)
	      ALPHA=Z*(1.25E-07*Z**2+0.99)
	      BETHA=Z*SQRT(2.475E-07*Z**2+0.9801)
		  U=(ALPHA+BETHA)**0.333333333
		  IF(ALPHA.GT.BETHA)THEN
		      V=(ALPHA-BETHA)**0.333333333
		  ELSE
		      V=0.
		  END IF
	  ASTAB(J)=ANINT((U+V+0.005*Z)**3)
   10 CONTINUE
C
	  IAA(13)=0
      DO 11 J=1,12
	  IAA(13-J)=-J
	  IAA(13+J)=+J
   11 CONTINUE
C
      RETURN
      END



      SUBROUTINE RESNPR(AZ,IOUT,IIN)
C Printing of residual nuclei (A,Z)-distributions
C Input: AZ - Array of residual nuclei (A,Z) distribution
C        IOUT,IIN - Number of nuclei outside/inside (A,Z) distribution
C
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      COMMON /RESN01/ ASTAB(100),IAA(25)
      COMMON /INFEL/ ZNUC(110),ATWEI(110),DENS(110),RIEV(110),SYMB(110)
      CHARACTER*2 SYMB
      DIMENSION AZ(25,100),IZA(25,100)
C
      WRITE(24,20)
   20 FORMAT(/1X,'RESIDUAL NUCLEI (A,Z)-DISTRIBUTION. A0(Z) IS THE',
     *       ' BETHA STABILITY CURVE.')
      WRITE(24,21)
   21 FORMAT(21X,'DISPOSITION OF  A  RELATIVE  A0(Z):') 
C
C Sending of real array AZ into integer array IZA
      DO 10 I=1,100
      DO 10 J=1,25
	  IZA(J,I)=NINT(AZ(J,I))
   10 CONTINUE
C
	  IMIN0=13
	  IMAX0=13
      DO 500 J=1,100
	  DO 501 I=1,25
	      IF(IZA(I,J).EQ.0)THEN
		  GO TO 501
	      ELSE
		  IMIN=I
		  GO TO 510
	      END IF
  501     CONTINUE
	GO TO 500
  510     DO 511 I=1,25
	      IF(IZA(26-I,J).EQ.0)THEN
		  GO TO 511
	      ELSE
		  IMAX=26-I
		  GO TO 520
	      END IF
  511     CONTINUE
  520   CONTINUE
	  IMIN0=MIN(IMIN0,IMIN)
	  IMAX0=MAX(IMAX0,IMAX)
  500 CONTINUE
C
      WRITE(24,22) (IAA(J),J=IMIN0,IMAX0)
   22 FORMAT(20X,25I6)
C
      DO 530 J=1,100
	  DO 531 I=1,25
	      IF(IZA(I,J).EQ.0)THEN
		  GO TO 531
	      ELSE
		  GO TO 532
	      END IF
  531     CONTINUE
	GO TO 530
  532     WRITE(24,23) SYMB(J),IFIX(ZNUC(J)),IFIX(ASTAB(J)),
     *                 (IZA(I,J),I=IMIN0,IMAX0)
   23     FORMAT(1X,A2,' (Z=',I3,',A0=',I3,'): ',25I6)
  530 CONTINUE
C
      WRITE(24,24) IIN
   24 FORMAT(' NUMBER OF NUCLEI  INSIDE DISTRIBUTION IS ',I10)
      WRITE(24,25) IOUT
   25 FORMAT(' NUMBER OF NUCLEI OUTSIDE DISTRIBUTION IS ',I10)
C
      RETURN
      END



      SUBROUTINE PRNSPD                              ! DCHAIN
C Printout of residual nuclei (A,Z)-distributions 
C as input file NMTC_YIELD for the DCHAIN-SP code
C
	COMMON /GDATA0/ NUMBOD,NUMZON                                ! GEMCA
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      COMMON /AZSPD/ AZSPD(280,100,0:1000),LZS50        ! DCHAIN
      DIMENSION AZ(280,100)
      COMMON /INFEL/ ZNUC(110),ATWEI(110),DENS(110),RIEV(110),SYMB(110)
      CHARACTER*2 SYMB
C
      WRITE(44,1)NUMTRE
    1 FORMAT(' regionwise nuclear yield (per projectile) Nstat=',I8)
      WRITE(44,2)
    2 FORMAT(' -----------------------------------------')
C
C Renormalization per one source projectile
      DO I=1,100
        DO J=1,280
          AZ(J,I)=AZSPD(J,I,0)
          AZ(J,I)=AZ(J,I)/REAL(NUMTRE)
        END DO
      END DO
C
      DO I=1,100   ! Loop on Z
        DO J=1,280
          IF(AZ(J,I).gt.0.)GOTO 10
        END DO
        JMIN=0
        JMAX=0
        GOTO 777   ! There is no yield with Z=I
C
   10   CONTINUE     ! There is yield with Z=I
C
        DO J=1,280
          IF(AZ(J,I).gt.0.)THEN
            JMIN=J                ! Amin for Z=I
            GOTO 20
          END IF
        END DO
C
   20 CONTINUE
        DO J=280,1,-1
          IF(AZ(J,I).gt.0.)THEN
            JMAX=J                ! Amax for Z=I
            GOTO 30
          END IF
        END DO
C
   30 CONTINUE
C Analysing of JMIN and JMAX: JMAX-JMIN+1 must be .LE.36
        IF(JMAX-JMIN+1 .gt. 36)THEN
          CONTINUE!WRITE(25,40)JMIN,JMAX,NUMTRE
   40     FORMAT(' PRNSPD Diagnostic: JMIN,JMAX,NUMTRE:',3I6)
          JMAX=JMIN+35   ! Restriction JMAX by hand
        END IF
C Cleaning of limits for printout by groups by 12 isotopes
        JMIN1=0
        JMAX1=0
        JMIN2=0
        JMAX2=0
        JMIN3=0
        JMAX3=0
C Limits for printout of one, two or three groups by 12 isotopes
        IF(JMAX-JMIN+1 .gt. 24)THEN
          JMIN1=JMIN
          JMAX1=JMIN1+11
          JMIN2=JMAX1+1
          JMAX2=JMIN2+11
          JMIN3=JMAX2+1
          JMAX3=JMAX
        ELSEIF(JMAX-JMIN+1 .gt. 12)THEN
          JMIN1=JMIN
          JMAX1=JMIN1+11
          JMIN2=JMAX1+1
          JMAX2=JMAX
        ELSEIF(JMAX-JMIN+1 .le. 12)THEN
          JMIN1=JMIN
          JMAX1=JMAX
	  END IF
C
C Printout of the first group
        WRITE(44,50)I,SYMB(I)
   50   FORMAT(/I5,'-',A2,' isotope production')
        WRITE(44,60)(J,J=JMIN1,JMAX1)
   60   FORMAT(' reg.',12I10)
                  DO INZON=1,NUMZON-1
        WRITE(44,70)INZON,(AZSPD(J,I,INZON)/REAL(NUMTRE),J=JMIN1,JMAX1)
   70   FORMAT(I5,12E10.4)
                  END DO
C
C Printout of the second group
        IF(JMIN2.gt.0)THEN
          WRITE(44,50)I,SYMB(I)
          WRITE(44,60)(J,J=JMIN2,JMAX2)
                  DO INZON=1,NUMZON-1
       WRITE(44,70)INZON,(AZSPD(J,I,INZON)/REAL(NUMTRE),J=JMIN2,JMAX2)
                  END DO
        END IF
C
C Printout of the third group
        IF(JMIN3.gt.0)THEN
          WRITE(44,50)I,SYMB(I)
          WRITE(44,60)(J,J=JMIN3,JMAX3)
                  DO INZON=1,NUMZON-1
        WRITE(44,70)INZON,(AZSPD(J,I,INZON)/REAL(NUMTRE),J=JMIN3,JMAX3)
                  END DO
        END IF
C
  777 END DO   ! end of Loop on Z
C
      RETURN
      END



      SUBROUTINE PRNSET                     ! userset         ! Fertman
      COMMON /ISOSET/ WISOSET(2000),FERTAZ(2000,50),  ! Fertman ! Truba06
     *                ISOSET(2,2000),TDEC(2000),NIS200,NISTRU   ! userset
      CHARACTER*10 TDEC   ! userset
      DIMENSION SUMZON(50),IZN(50)
C
C Auxiliary arrays
      DO J=1,50
        IZN(J)=J      ! Numbers of zons
        SUMZON(J)=0.  ! Total number of nuclei-products in zone No.J
        DO K=1,NISTRU
          SUMZON(J)=SUMZON(J)+FERTAZ(K,J)
        END DO
      END DO
C
C Header of the Table
      CONTINUE!WRITE(46,10)(IZN(J),J=1,50)
   10 FORMAT(/' No.',2X,'( Z,  A):','  T1/2   ',' All_zons',
     *        9('    #',I1),41('   #',I2)/1X,332('-'))
C
C Printout of the Table
      SUMISOT=0.
      DO J=1,NISTRU
        CONTINUE!WRITE(46,20)J,ISOSET(1,J),ISOSET(2,J),TDEC(J),WISOSET(J),
!     *                (INT(FERTAZ(J,K)),K=1,50)
   20   FORMAT(I4,2X,'(',I2,',',I3,')',A10,F9.1,50I6)
      SUMISOT=SUMISOT+WISOSET(J)
      END DO
C
C Printout of sums of all nuclei-products in each zone
      CONTINUE!WRITE(46,30)SUMISOT,(INT(SUMZON(J)),J=1,50)
   30 FORMAT(1X,332('-')/18X,'TOTAL:',F9.1,50I6)
C
      RETURN
      END



      SUBROUTINE CLECON
C Cleaning of arrays for AGT-conversion protocol
      COMMON /CONVEN/ ECONV(65),NCONV(65),JCONV(65),UEXCIT,ENKIN
      DO J=1,65
	  ECONV(J)=0.
	  NCONV(J)=0
	  JCONV(J)=0
      END DO
	  UEXCIT=0.
	  ENKIN=0.
      RETURN
      END



      SUBROUTINE PRICON
C Printout of AGT-conversion protocol
      CHARACTER*8  PNAME
      COMMON /PNAME/ PNAME(65),IDENTP(65)
      COMMON /CONVEN/ ECONV(65),NCONV(65),JCONV(65),UEXCIT,ENKIN
C
      DIMENSION PARTIC(25)                                           !HION6
      CHARACTER*7 PARTIC
      DATA PARTIC /'NEUTRON',' PROTON','  PI-  ','  PI+  ','  PI0  ',
     *             ' ANTI N',' ANTI P','  K-   ','  K+   ','  K0   ',
     *             'ANTI K0',' GAMMA ','ELECTRN','POSITRN','  MU-  ',
     *             '  MU+  ','  NU E ','ANTINUE',' NU MU ','ANTINUM',
     *             'DEUTRON','TRITIUM','HELIUM3',' ALPHA ','HEV*ION'/ !HION6
C
      DIMENSION IDPSHL(20) ! AGT-codes of particles, accepted in SHIELD
C                   IDENTP[i]  !   i   Particle    JPART
      DATA IDPSHL/             !
     *                  120,   !  [1]   PI+           4
     *                 -120,   !  [2]   PI-           3
     *                  130,   !  [3]   K+            9
     *                 -130,   !  [4]   K-            8
     *                  230,   !  [5]   K0           10
     *                 -230,   !  [6]   AK0          11
     *                  110,   !  [7]   PI0           5
     *                -1120,   ! [19]   AP            7
     *                -1220,   ! [20]   AN            6
     *                 1120,   ! [37]   P             2
     *                 1220,   ! [38]   N             1
     *                   10,   ! [57]   GM           12
     *                   12,   ! [58]   E+           14
     *                  -12,   ! [59]   E-           13
     *                  -14,   ! [60]   MU+          16
     *                   14,   ! [61]   MU-          15
     *                   11,   ! [62]   NUE          17
     *                   13,   ! [63]   NUM          19
     *                  -11,   ! [64]   ANUE         18
     *                  -13/   ! [65]   ANUM         20
C
      WRITE(24,1)
    1 FORMAT(/5X,'----- Protocol of AGT to SHIELD particles',
     *          ' conversion -----')
      WRITE(24,2)
    2 FORMAT(5X,'AGT',23X,'SHIELD     SUM Ekin  No.of Part.')
C
C Counters of total energies and particle numbers
      TSHLD=0.      ! Summary energy of SHIELD particles
      NSHLD=0       ! Total number of SHIELD particles
	TCON=0.     ! Summary energy of CONVERTED particles
	NCON=0      ! Total number of CONVERTED particles
	  TIGN=0.   ! Summary energy of IGNORED particles
	  NIGN=0    ! Total number of IGNORED particles
C      
C ====================== Main loop ============================     
      DO 100 J=1,65
	  IDP=IDENTP(J)   ! AGT particle code
	  JPR= JCONV(J)   ! SHIELD particle code
	     DO K=1,20
		 IF(IDP.EQ.IDPSHL(K))GO TO 20  ! IDP is accepted in SHIELD
	     END DO
	     GO TO 30                       ! IDP isn't accepted in SHIELD
C
   20     CONTINUE
C---------SHIELD particles   
	  IF(NCONV(J).NE.0)THEN    ! Output only if the particles exist
	    TSHLD=TSHLD+ECONV(J)
	    NSHLD=NSHLD+NCONV(J)
	    WRITE(24,21)PARTIC(JCONV(J)),ECONV(J),NCONV(J)
   21       FORMAT(30X,A7,': ',E11.5,I12)
	  END IF
	  GO TO 100
C
   30     CONTINUE
	  IF(NCONV(J).NE.0)THEN    ! Output only if the particles exist
	    IF(JPR.NE.0)THEN  
C-------------CONVERTED particles   
	      TCON=TCON+ECONV(J)
	      NCON=NCON+NCONV(J)
	      WRITE(24,31)PNAME(J),PARTIC(JCONV(J)),ECONV(J),NCONV(J)
   31         FORMAT(5X,A8,' is converted to ',A7,': ',E11.5,I12)
	      GO TO 100
	    ELSE
C-------------IGNORED particles   
	      TIGN=TIGN+ECONV(J)
	      NIGN=NIGN+NCONV(J)
	      WRITE(24,32)PNAME(J),ECONV(J),NCONV(J)
   32         FORMAT(5X,A8,' is ignored',13X,': ',E11.5,I12)
	      GO TO 100
	    END IF
	  END IF
C
  100 CONTINUE
C =======================================================
C
      WRITE(24,40)
   40 FORMAT(5X,'--------------------------------:')
      WRITE(24,41)TSHLD,NSHLD
   41 FORMAT(21X,'SHIELD particles: ',E11.5,I12)
      WRITE(24,42)TCON,NCON
   42 FORMAT(18X,'CONVERTED particles: ',E11.5,I12)
      WRITE(24,43)TIGN,NIGN
   43 FORMAT(20X,'IGNORED particles: ',E11.5,I12)
      WRITE(24,44)UEXCIT
   44 FORMAT(10X,'EXCITATION energy of nuclei: ',E11.5)
      WRITE(24,45)ENKIN
   45 FORMAT(13X,'KINETIC energy of nuclei: ',E11.5)
C
      RETURN
      END



      SUBROUTINE PUTFIL(LUNIT1,LUNIT2)
C Program for including of some text file FILE1 (LREC=76) 
C at unit LUNIT1 into other text File FILE2 at unit LUNIT2
C--------------------------------------------------------
      CHARACTER*76 LINE
C
      WRITE(LUNIT2,20)
   20 FORMAT(/10X,'CG-GEOMETRY INPUT FILE.')
C
C Input of the next line
    1 READ(LUNIT1,10,END=100)LINE
   10 FORMAT(A76)
      WRITE(LUNIT2,10)LINE
      GO TO 1
C
  100 RETURN
      END
