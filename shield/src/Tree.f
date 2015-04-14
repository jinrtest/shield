
***********************************************************************
*                                                                     *
*        The following subroutines are written by N.M.Sobolevsky,     *
*        Institute for Nuclear Research RAS, Moscow,                  *
*        E-mail: sobolevs@AL20.inr.troitsk.ru                         *
*                                                                     *
*        Rewritten by A. Timofeev in 2012                             *
***********************************************************************

      SUBROUTINE DECTRE
C     Decoding cascade tree
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

      CALL TREERES          ! A.Timofeev - external C-call to reset the tree data

C Tree is empty. No cascade tree was generated. Incident particle
C was removed to source array (20.01.97).
      IF(NBRMAX.eq.0)THEN
          
    	  RETURN
      END IF
C
C Primary particle
      CALL TREEDEC(0,0,0,0,0)
	  KST=ITREE(7,3)
      IF(KST.EQ.1 .OR. KST.EQ.2 .OR. KST.EQ.5 .OR. KST.EQ.6)THEN
C	  CONTINUE!WRITE(24,4)
C    4     FORMAT()
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
	      GO TO 80  ! Modified by A.Timofeev
                    ! Fly-outs and absorption are also registered
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
 1000                 FORMAT(' DECTRE DIAGNOSTIC: ARRAY TREVS EXCEEDED')
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
   80     CALL TREEDEC(KSTATE,JBRST,JBRFN,NPCAS,NPTOT)  ! A.Timofeev
C  Calling external c-func
C   80     CALL PRSTAR(KSTATE,JBRST,JBRFN,NPCAS,NPTOT)  ! A.Timofeev
   99     JBRST=JBRFN+1
C
  100 CONTINUE
C------------------- END OF CYCLE ON JBR --------------------------------
      RETURN
      END
