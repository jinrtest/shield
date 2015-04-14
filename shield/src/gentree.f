
***********************************************************************
*                                                                     *
*        The following subroutines are written by N.M.Sobolevsky,     *
*        Institute for Nuclear Research RAS, Moscow,                  *
*        E-mail: sobolevs@AL20.inr.troitsk.ru                         *
*                                                                     *
***********************************************************************

      SUBROUTINE GENTRE(JSTAT)                                       ! GEMCA
C
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      COMMON /TREE/ TREE(15,400000),ITREE(7,400000),LTR13,LTR7,LTR300 ! GSITR
      COMMON /MANALL/ NTREE,LINE1,LINE2,NBRMAX,NPRMAX,NNEUTR,NGAMMA,
     *                MELSC,MINTER,MCUT,MOUT,MABSOR,MDECAY,NSTRAN,
     *                NEMIN,NEPLS,NNUEL,NANUE,NNUMU,NANUM
      COMMON /ISLERR/ ISLERR   ! HADES1
C
C NUMBER OF TREE
      NTREE=JSTAT
      CYCLECOUNTER=0
C
C INITIAL STATUS OF RANDOM NUMBER GENERATOR FOR TREE NO.NTREE
    1 IXINIT=IX
      NUMTRE=NTREE
      
      LOOPCOUNTER=LOOPCOUNTER+1
      IF(LOOPCOUNTER.GT.100000)THEN
        PRINT *, 'GENTRE: EXIT FROM GENTRE BY LOOP COUNTER!'
        JSTAT=JSTAT-1
        CALL CLOSE(0)
        RETURN
      END IF
C
C CLEANING OF TREE AND SOURCES
      CALL CLEAR
C
C ESTABLISHMENT OF CONTROL VARIABLES AND COUNTERS
      CALL ESTABL
C
C SAMPLING OF PRIMARY PARTICLE AND SENDING IT INTO TREE
      CALL INSPAR
C
C TAKING OF THE PARTICLE, WHICH INITIATE THE NEW BRANCH, FROM TREE
  100 CALL TAKPAR(ICOMPL)
	  IF(ICOMPL.EQ.1)THEN
	      CALL CLOSE(0)
	      RETURN
	  ELSE
	      CONTINUE
	  END IF
C
C SAMPLING OF RANGE
      CALL RANGE(INOVER)                                             ! GEMCA
	  IF(INOVER.EQ.1)GO TO 1
      CALL SENS(KOVER)
	  IF(KOVER.EQ.1)GO TO 1
C
C DETERMINATION OF THE PARTICLE STATUS IN THE END OF RANGE
      KSTATE=ITREE(7,LINE2+1)
      IF(KSTATE.EQ.7)GO TO 29
      IF(KSTATE.EQ.3)GO TO 30
      IF(KSTATE.EQ.4)GO TO 31
      IF(KSTATE.EQ.5 .OR. KSTATE.EQ.6)GO TO 40
      IF(KSTATE.LE.2)GO TO 50
C
      CONTINUE!WRITE(25,1001)KSTATE,IXINIT,NUMTRE
 1001 FORMAT(' GENTRE DIAGNOSTIC: UNRECOGNIZED PARTICLE STATUS AFTER',
     *' RANGE.'/20X,'KSTATE=',I3,', IXINIT=',I12,', NUMTRE=',I6)
      CALL CLOSAC
      RETURN
C
C PARTICLE ENERGY BELOW CUT OFF
   29 MCUT=MCUT+1
      GO TO 32
C
C PARTICLE GO OUT FROM TARGET
   30 MOUT=MOUT+1
      GO TO 32
C
C PARTICLE ABSORPTION IN TARGET BECAUSE OF IONISATION LOSS
   31 MABSOR=MABSOR+1
      GO TO 32
C                 BRANCH COMPLETED, NEW LINE2, CONTROL OF OVERFILLING OF TREE
   32 LINE2=LINE2+2
	  IF(LINE2.GE.LTR300-2)THEN
	      CONTINUE!WRITE(25,1002)IXINIT,NUMTRE
 1002         FORMAT(' GENTRE DIAGNOSTIC: ARRAY TREE EXCEEDED, IXINIT=',
     *                I12,', NUMTRE=',I6)
	      GO TO 1
	  ELSE
	      GO TO 100
	  END IF
C
C DECAY OF PARTICLE
   40 MDECAY=MDECAY+1
      CALL DECAY1
      GO TO 80
C
C ELASTIC OR INELASTIC INTERACTION, PREPARATORY ACTIONS
   50 CALL PREPAR
C After PREPAR variable KSTATE can has only three values: 1,2 or 6
      KSTATE=ITREE(7,LINE2+1)
	  IF(KSTATE.EQ.6)GO TO 40
      GO TO(60,70),KSTATE
C
C ELASTIC SCATTERING
   60 MELSC=MELSC+1
      CALL ELSCAT
C      CALL KARSJH(25)   ! Debug printout
      GO TO 80
C
C INELASTIC INTERACTION
   70 MINTER=MINTER+1
      CALL SREACT(INOVERT)                                  ! KTH
      if(ISLERR.eq.1)then   ! HADES1
        CONTINUE!WRITE(25,1004)NUMTRE
 1004   format(' GENTREE DIAGNOSTIC: abnormal return from SREACT, ',
     *           '(HADES1) NUMTRE=',I9)
        goto 1
      end if
C      CALL KARSJH(25)   ! Debug printout
C==================== KTH ======================
      IF(INOVERT.eq.1)THEN
        CONTINUE!WRITE(25,1003)NUMTRE
 1003   FORMAT(' GENTREE DIAGNOSTIC: abnormal return from SREACT, ',
     *           'NUMTRE=',I9)
        GO TO 1
      END IF
C================= end KTH ====================
      GO TO 80
C
C SENDING OF THE INTERACTION PRODUCTS INTO TREE AND SOURCES, BRANCH COMPLETED
   80 CALL SEND(IOVER)
	  IF(IOVER.EQ.1)THEN
	      GO TO 1
	  ELSE
	      GO TO 100
	  END IF
C
      END


      SUBROUTINE CLEAR
C CLEANING OF ARRAYS BEFORE GENERATION OF THE NEW TREE
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
      COMMON /WRPR/ WP(15),IWP(7)                                   !HION1
C
C---------------------------TREE
      DO 10 J=1,LTR300
	  DO 11 I=1,LTR13
   11     TREE(I,J)=0.
	  DO 12 I=1,LTR7
   12     ITREE(I,J)=0
   10 CONTINUE
C---------------------------SOURCES
      DO 20 J=1,LSN100
	  DO 21 I=1,LSN13
   21     SOURN(I,J)=0.
	  DO 22 I=1,LSN7
   22     ISOURN(I,J)=0
   20 CONTINUE
C
      DO 30 J=1,LSG100
	  DO 31 I=1,LSG13
   31     SOURG(I,J)=0.
	  DO 32 I=1,LSG7
   32     ISOURG(I,J)=0
   30 CONTINUE
C
      DO 130 J=1,L13100
	  DO 131 I=1,L1313
  131     SEMIN(I,J)=0.
	  DO 132 I=1,L137
  132     ISEMIN(I,J)=0
  130 CONTINUE
C
      DO 230 J=1,L14100
	  DO 231 I=1,L1413
  231     SEPLS(I,J)=0.
	  DO 232 I=1,L147
  232     ISEPLS(I,J)=0
  230 CONTINUE
C
      DO 330 J=1,L17100
	  DO 331 I=1,L1713
  331     SNUEL(I,J)=0.
	  DO 332 I=1,L177
  332     ISNUEL(I,J)=0
  330 CONTINUE
C
      DO 430 J=1,L18100
	  DO 431 I=1,L1813
  431     SANUE(I,J)=0.
	  DO 432 I=1,L187
  432     ISANUE(I,J)=0
  430 CONTINUE
C
      DO 530 J=1,L19100
	  DO 531 I=1,L1913
  531     SNUMU(I,J)=0.
	  DO 532 I=1,L197
  532     ISNUMU(I,J)=0
  530 CONTINUE
C
      DO 630 J=1,L20100
	  DO 631 I=1,L2013
  631     SANUM(I,J)=0.
	  DO 632 I=1,L207
  632     ISANUM(I,J)=0
  630 CONTINUE
C------------------------------TRAN
      DO 50 J=1,LSR500
	  DO 51 I=1,LSR5
   51     TRAN(I,J)=0.
	  DO 52 I=1,LSR6
   52     ITRAN(I,J)=0
   50 CONTINUE
C-------------------------------WP
      DO 41 J=1,LTR13
   41 WP(J)=0.
      DO 42 J=1,LTR7
   42 IWP(J)=0
C
      RETURN
      END


      SUBROUTINE ESTABL
C ESTABLISHMENT OF CONTROL VARIABLES AND COUNTERS
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      COMMON /MANALL/ NTREE,LINE1,LINE2,NBRMAX,NPRMAX,NNEUTR,NGAMMA,
     *                MELSC,MINTER,MCUT,MOUT,MABSOR,MDECAY,NSTRAN,
     *                NEMIN,NEPLS,NNUEL,NANUE,NNUMU,NANUM
C
      LINE1=2
      LINE2=3
	 NBRMAX=0
	 NPRMAX=0
	    NNEUTR=2
	    NGAMMA=2
	    NEMIN=2
	    NEPLS=2
	    NNUEL=2
	    NANUE=2
	    NNUMU=2
	    NANUM=2
	       MELSC=0
	       MINTER=0
	       MCUT=0
	       MOUT=0
	       MABSOR=0
	       MDECAY=0
		  NSTRAN=2
C
      RETURN
      END



      SUBROUTINE TAKPAR(ICOMPL)
C TAKING OF THE PARTICLE, WHICH INITIATE THE NEW BRANCH, FROM TREE
C     OUTPUT VARIABLE  ICOMPL=+1  IF THE TREE IS COMPLETED, ELSE ICOMPL=0
C
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      COMMON /TREE/ TREE(15,400000),ITREE(7,400000),LTR13,LTR7,LTR300 ! GSITR
      COMMON /MANALL/ NTREE,LINE1,LINE2,NBRMAX,NPRMAX,NNEUTR,NGAMMA,
     *                MELSC,MINTER,MCUT,MOUT,MABSOR,MDECAY,NSTRAN,
     *                NEMIN,NEPLS,NNUEL,NANUE,NNUMU,NANUM
          COMMON /HIPROJ/ APROJ,ZPROJ   ! For AGT                    !HION4
C
C FIND THE PARTICLE WITH KPASS=-1
  100 DO 10 J=LINE1,LINE2-1
      IF(ITREE(1,J).NE.0 .AND. ITREE(6,J).EQ.-1)THEN
	  LINE=J
	  GO TO 20
      ELSE
	  GO TO 10
      END IF
   10 CONTINUE
C
C THERE IS NO PARTICLE WITH KPASS=-1 IN TREE. THE TREE IS COMPLETED.
      ICOMPL=+1
      RETURN
C
C COPY THE PARTICLE WITH KPASS=-1 FROM LINE TO LINE2
   20 DO 21 J=1,LTR13
   21 TREE(J,LINE2)=TREE(J,LINE)
      DO 22 J=1,LTR7
   22 ITREE(J,LINE2)=ITREE(J,LINE)
          APROJ=TREE(14,LINE2)                                       !HION4
          ZPROJ=TREE(15,LINE2)                                       !HION4
C
C REPLACE KPASS=-1 ON KPASS=0 IN LINE2
      ITREE(6,LINE2)=0
C
C NUMBER OF BRANCH FROM WHICH THE INITIATING PARTICLE WAS TAKEN.
      ITREE(4,LINE2)=ITREE(3,LINE)
C
C GIVING OF THE NEW VALUES OF CONTROL VARIABLES AND UNICAL NUMBERS
      LINE1=LINE+1
	  NBRMAX=NBRMAX+1
	     ITREE(3,LINE2)=NBRMAX
C
C DISPLACEMENT OF TREE ON ONE LINE BACKWARDS AND RUB OUT PARTICLE IN LINE
C NO. LINE
      DO 23 J=LINE1,LINE2+1
	  DO 24 I=1,LTR13
   24     TREE(I,J-1)=TREE(I,J)
	  DO 25 I=1,LTR7
   25     ITREE(I,J-1)=ITREE(I,J)
   23 CONTINUE
	      LINE1=LINE1-1
	      LINE2=LINE2-1

C
      ICOMPL=0
      RETURN
      END



      SUBROUTINE CLOSE(JABNOR)
C CLOSE ARRAYS WHEN THE TREE GENERATION IS FINISHED
C     INPUT VARIABLE JABNOR=0 IF NORMAL CLOSING OCCURES.
C     IN CASE ACCIDENTAL CLOSING JABNOR=1.
C
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      COMMON /TREE/ TREE(15,400000),ITREE(7,400000),LTR13,LTR7,LTR300 ! GSITR
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
      COMMON /TRAN/ TRAN(5,600000),ITRAN(6,600000),LSR5,LSR6,LSR500  ! GSITR
      DIMENSION VP(15),IVP(7),UP(5),IUP(6)                            !HION1
C
      DO 1 J=1,LTR13
    1 VP(J)=0.
      DO 2 J=1,LTR7
    2 IVP(J)=0
      DO 3 J=1,LSR5
    3 UP(J)=0.
      DO 4 J=1,LSR6
    4 IUP(J)=0
C
C COUNTERS OF THE DIFFERENT EVENTS IN THE TREE
      VP(1)=FLOAT(MELSC)
      VP(2)=FLOAT(MINTER)
      VP(3)=FLOAT(MCUT)
      VP(4)=FLOAT(MOUT)
      VP(5)=FLOAT(MABSOR)
      VP(6)=FLOAT(MDECAY)
C
C THE CONTROL VARIABLES
      IVP(2)=NTREE
      IVP(3)=LINE2-2
      IVP(4)=NBRMAX
      IVP(5)=NPRMAX
      IVP(6)=NNEUTR-2
      IVP(7)=NGAMMA-2
	  IUP(2)=NTREE
	  IUP(3)=NSTRAN-2
	  IUP(4)=NBRMAX
C
C RECORD OF THE COUNTERS AND CONTROL VARIABLES INTO THE FIRST LINES OF ARRAYS
      DO 11 J=1,LTR13
   11 TREE(J,1)=VP(J)
      DO 12 J=1,LTR7
   12 ITREE(J,1)=IVP(J)
C
      DO 21 J=1,LSN13
   21 SOURN(J,1)=VP(J)
      DO 22 J=1,LSN7
   22 ISOURN(J,1)=IVP(J)
	  ISOURN(LSN7,1)=NNEUTR-2
C
      DO 31 J=1,LSG13
   31 SOURG(J,1)=VP(J)
      DO 32 J=1,LSG7
   32 ISOURG(J,1)=IVP(J)
	  ISOURG(LSG7,1)=NGAMMA-2
C
      DO 131 J=1,L1313
  131 SEMIN(J,1)=VP(J)
      DO 132 J=1,L137
  132 ISEMIN(J,1)=IVP(J)
	  ISEMIN(L137,1)=NEMIN-2
C
      DO 231 J=1,L1413
  231 SEPLS(J,1)=VP(J)
      DO 232 J=1,L147
  232 ISEPLS(J,1)=IVP(J)
	  ISEPLS(L147,1)=NEPLS-2
C
      DO 331 J=1,L1713
  331 SNUEL(J,1)=VP(J)
      DO 332 J=1,L177
  332 ISNUEL(J,1)=IVP(J)
	  ISNUEL(L177,1)=NNUEL-2
C
      DO 431 J=1,L1813
  431 SANUE(J,1)=VP(J)
      DO 432 J=1,L187
  432 ISANUE(J,1)=IVP(J)
	  ISANUE(L187,1)=NANUE-2
C
      DO 531 J=1,L1913
  531 SNUMU(J,1)=VP(J)
      DO 532 J=1,L197
  532 ISNUMU(J,1)=IVP(J)
	  ISNUMU(L197,1)=NNUMU-2
C
      DO 631 J=1,L2013
  631 SANUM(J,1)=VP(J)
      DO 632 J=1,L207
  632 ISANUM(J,1)=IVP(J)
	  ISANUM(L207,1)=NANUM-2
C
      DO 41 J=1,LSR5
   41 TRAN(J,1)=UP(J)
      DO 42 J=1,LSR6
   42 ITRAN(J,1)=IUP(J)
C
C THERE IS QUESTION ABOUT DIMENSIONS LTR13,LTR7 AND LSN13,LSG13,LSN7,LSG7
C
C RECORD OF THE "END OF TREE" IN THE LAST LINES OF ARRAYS
      IF(LINE2.GT.LTR300)THEN
	  ITREE(2,LTR300)=-NTREE
	  ITREE(3,LTR300)=JABNOR
      ELSE
	  ITREE(2,LINE2)=-NTREE
	  ITREE(3,LINE2)=JABNOR
      END IF
C
      IF(NNEUTR.GT.LSN100)THEN
	  ISOURN(2,LSN100)=-NTREE
	  ISOURN(3,LSN100)=JABNOR
      ELSE
	  ISOURN(2,NNEUTR)=-NTREE
	  ISOURN(3,NNEUTR)=JABNOR
      END IF
C
      IF(NGAMMA.GT.LSG100)THEN
	  ISOURG(2,LSG100)=-NTREE
	  ISOURG(3,LSG100)=JABNOR
      ELSE
	  ISOURG(2,NGAMMA)=-NTREE
	  ISOURG(3,NGAMMA)=JABNOR
      END IF
C
      IF(NEMIN.GT.L13100)THEN
	  ISEMIN(2,L13100)=-NTREE
	  ISEMIN(3,L13100)=JABNOR
      ELSE
	  ISEMIN(2,NEMIN)=-NTREE
	  ISEMIN(3,NEMIN)=JABNOR
      END IF
C
      IF(NEPLS.GT.L14100)THEN
	  ISEPLS(2,L14100)=-NTREE
	  ISEPLS(3,L14100)=JABNOR
      ELSE
	  ISEPLS(2,NEPLS)=-NTREE
	  ISEPLS(3,NEPLS)=JABNOR
      END IF
C
      IF(NNUEL.GT.L17100)THEN
	  ISNUEL(2,L17100)=-NTREE
	  ISNUEL(3,L17100)=JABNOR
      ELSE
	  ISNUEL(2,NNUEL)=-NTREE
	  ISNUEL(3,NNUEL)=JABNOR
      END IF
C
      IF(NANUE.GT.L18100)THEN
	  ISANUE(2,L18100)=-NTREE
	  ISANUE(3,L18100)=JABNOR
      ELSE
	  ISANUE(2,NANUE)=-NTREE
	  ISANUE(3,NANUE)=JABNOR
      END IF
C
      IF(NNUMU.GT.L19100)THEN
	  ISNUMU(2,L19100)=-NTREE
	  ISNUMU(3,L19100)=JABNOR
      ELSE
	  ISNUMU(2,NNUMU)=-NTREE
	  ISNUMU(3,NNUMU)=JABNOR
      END IF
C
      IF(NANUM.GT.L20100)THEN
	  ISANUM(2,L20100)=-NTREE
	  ISANUM(3,L20100)=JABNOR
      ELSE
	  ISANUM(2,NANUM)=-NTREE
	  ISANUM(3,NANUM)=JABNOR
      END IF
C
      IF(NSTRAN.GT.LSR500)THEN
	  ITRAN(2,LSR500)=-NTREE
	  ITRAN(3,LSR500)=JABNOR
      ELSE
	  ITRAN(2,NSTRAN)=-NTREE
	  ITRAN(3,NSTRAN)=JABNOR
      END IF
C
      RETURN
      END


      SUBROUTINE CLOSAC
C ACCIDENTAL CLOSING OF THE ARRAYS TREE,SOURN,SOURG,TRAN
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      COMMON /MANALL/ NTREE,LINE1,LINE2,NBRMAX,NPRMAX,NNEUTR,NGAMMA,
     *                MELSC,MINTER,MCUT,MOUT,MABSOR,MDECAY,NSTRAN,
     *                NEMIN,NEPLS,NNUEL,NANUE,NNUMU,NANUM
      COMMON /MANLOC/ KPRINU,KSECNU,KCASC,KNEUTR,KGAMMA,
     *                KEMIN,KEPLS,KNUEL,KANUE,KNUMU,KANUM
      COMMON /STERAN/ STEPR(9,500),LST9,LST50
C
      CONTINUE!WRITE(25,10)NTREE,NBRMAX,NSTRAN
   10 FORMAT(' OVERFILLING OF THE TREE OR SOURCES OR TRAN'/
     *       ' NTREE=',I6,5X,'NBRMAX=',I6,5X,'NSTRAN=',I6)
      CONTINUE!WRITE(25,20)LINE2,KPRINU,KSECNU,KCASC
   20 FORMAT(10X,'LINE2=',I5,5X,'KPRINU=',I1,5X,'KSECNU=',I2,5X,
     *           'KCASC=',I4)
      CONTINUE!WRITE(25,30)NNEUTR,KNEUTR
   30 FORMAT(10X,'NNEUTR=',I4,5X,'KNEUTR=',I4)
      CONTINUE!WRITE(25,40)NGAMMA,KGAMMA
   40 FORMAT(10X,'NGAMMA=',I4,5X,'KGAMMA=',I4)
      CONTINUE!WRITE(25,110)NEMIN,KEMIN
  110 FORMAT(10X,'NEMIN=',I4,5X,'KEMIN=',I4)
      CONTINUE!WRITE(25,120)NEPLS,KEPLS
  120 FORMAT(10X,'NEPLS=',I4,5X,'KEPLS=',I4)
      CONTINUE!WRITE(25,130)NNUEL,KNUEL
  130 FORMAT(10X,'NNUEL=',I4,5X,'KNUEL=',I4)
      CONTINUE!WRITE(25,140)NANUE,KANUE
  140 FORMAT(10X,'NANUE=',I4,5X,'KANUE=',I4)
      CONTINUE!WRITE(25,150)NNUMU,KNUMU
  150 FORMAT(10X,'NNUMU=',I4,5X,'KNUMU=',I4)
      CONTINUE!WRITE(25,160)NANUM,KANUM
  160 FORMAT(10X,'NANUM=',I4,5X,'KANUM=',I4)
	  ISTEP=IFIX(STEPR(2,1))
      CONTINUE!WRITE(25,50)ISTEP
   50 FORMAT(10X,'ISTEP=',I4)
C
      CALL CLOSE(1)
      RETURN
      END


      SUBROUTINE SEND(IOVER)
C SENDING OF THE SECOND NUCLEI AND PARTICLES FROM ARRAYS SNU,SPT
C INTO ARRAYS TREE, SOURN, SOURG
C     OUTPUT VARIABLE  IOVER=+1 IF ARRAYS  TREE, SOURN, SOURG ARE OVERFILLED,
C     ELSE IOVER=0
C
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      COMMON /TREE/ TREE(15,400000),ITREE(7,400000),LTR13,LTR7,LTR300 ! GSITR
c
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
      COMMON /MANLOC/ KPRINU,KSECNU,KCASC,KNEUTR,KGAMMA,
     *                KEMIN,KEPLS,KNUEL,KANUE,KNUMU,KANUM
      COMMON /SECHAR/ SPT(6,5000),SNU(10,101),LS6,LS100,LS10,LS11  ! GSITR
C================================= KarSjh4 ==================================
C Flag INGEN cancels the appropriation of the regular NUMBER OF GENERATION to
C products of ELASTIC nuclear interactions only. Products of other interac-
C tions (inelastic, decays) get the generation number by 1 more then number of
C the initiating projectile.
      COMMON /INGEN/ INGEN                        ! KarSjh4
C=============================== end KarSjh4 ================================
C
C CLEANING OF THE LOCAL CONTROL VARIABLES
	  KSECNU=0
	  KCASC=0
	  KNEUTR=0
	  KGAMMA=0
	      KEMIN=0
	      KEPLS=0
	      KNUEL=0
	      KANUE=0
	      KNUMU=0
	      KANUM=0
C
C DETERMINATION OF THE INTERACTION TYPE
      KSTATE=ITREE(7,LINE2+1)
      INGEN=KSTATE                                 ! KarSjh4
      IF(KSTATE.EQ.5 .OR. KSTATE.EQ.6)GO TO 10
      IF(KSTATE.EQ.1 .OR. KSTATE.EQ.2)GO TO 20
      CONTINUE!WRITE(25,1003)KSTATE
 1003 FORMAT(' UNRECOGNIZED PARTICLE STATUS IN SUBR. SEND,  KSTATE=',I3)
      GO TO 55
C
C DECAY
   10 KPRINU=0
      GO TO 50
C
C ELASTIC SCATTERING OR INELASTIC INTERACTION
   20 KPRINU=1
      GO TO 40
C
C SENDING OF THE NUCLEUS-TARGET AND NUCLEI-PRODUCTS FROM  SNU  INTO  TREE
C AND DETERMINATION OF THE NUCLEI-PRODUCTS NUMBER   KSECNU
   40 CALL SENDNU(1,LINE2+1+KPRINU)
c No. of line containing the first residual nucleus                !HION1
      nucfir=LINE2+1+KPRINU+1                                      !HION1
C
      DO 41 J=2,LS11
      KINTER=IFIX(SNU(4,J))
	  IF(KINTER.GT.0)THEN
	      KSECNU=KSECNU+1
              KK=LINE2+1+KPRINU+KSECNU
		  IF(KK.GE.LTR300)THEN
		      LINE2=KK
		      NNEUTR=NNEUTR+KNEUTR
		      NGAMMA=NGAMMA+KGAMMA
		      NEMIN=NEMIN+KEMIN
		      NEPLS=NEPLS+KEPLS
		      NNUEL=NNUEL+KNUEL
		      NANUE=NANUE+KANUE
		      NNUMU=NNUMU+KNUMU
		      NANUM=NANUM+KANUM
		      GO TO 55
		  ELSE
		      CONTINUE
		  END IF
	      CALL SENDNU(J,KK)
c             No. of line containing the last residual nucleus     !HION1
              nuclas=LINE2+1+KPRINU+KSECNU                         !HION1
	      GO TO 41
	  ELSE
	      GO TO 50
	  END IF
   41 CONTINUE
C
C SENDING OF THE SECOND PARTICLES INTO  TREE,SOURN,SOURG
C DETERMINATION OF THE PARTICLE NUMBERS   KCASC,KNEUTR,KGAMMA
   50 CONTINUE
c  No. of line containing the first cascade particle               !HION1
C      if(KPRINU.eq.1)nprfir=LINE2+1+KPRINU+KSECNU+1                !HION1
C================= KarSjh ==================     
	if(KPRINU.eq.1)then
        nprfir=LINE2+1+KPRINU+KSECNU+1
        nprlas=nprfir-1 ! for the case when secondary particles are absent.
      end if			  ! There was no error, but gives wrong debug printout.
C================ end KarSjh ===============
      DO 51 J=1,LS100
      JPART=IFIX(SPT(5,J))
	  IF(JPART.GT.0)THEN
C                       THERE IS THE NEXT PARTICLE IN ARRAY SPT
              IF(JPART.EQ.1 .AND. SPT(4,J).LT.14.5)THEN           ! BNAB28
C
C             -------------------------------------LOW ENERGY NEUTRON---------
		  KNEUTR=KNEUTR+1
		  KK=NNEUTR+KNEUTR-1
		      IF(KK.GE.LSN100)THEN
			  LINE2=LINE2+1+KPRINU+KSECNU+KCASC+1
			  NNEUTR=KK
			  NGAMMA=NGAMMA+KGAMMA
			  NEMIN=NEMIN+KEMIN
			  NEPLS=NEPLS+KEPLS
			  NNUEL=NNUEL+KNUEL
			  NANUE=NANUE+KANUE
			  NNUMU=NNUMU+KNUMU
			  NANUM=NANUM+KANUM
			  GO TO 55
		      ELSE
			  CONTINUE
		      END IF
		  CALL SENDPR(J,KK,2)
		  GO TO 51
	      ELSE
		  CONTINUE
	      END IF
C
	      IF(JPART.EQ.12)THEN
C
C             --------------GAMMA-QUANT---------------------------------------
		  KGAMMA=KGAMMA+1
		  KK=NGAMMA+KGAMMA-1
		      IF(KK.GE.LSG100)THEN
			  LINE2=LINE2+1+KPRINU+KSECNU+KCASC+1
			  NNEUTR=NNEUTR+KNEUTR
			  NGAMMA=KK
			  NEMIN=NEMIN+KEMIN
			  NEPLS=NEPLS+KEPLS
			  NNUEL=NNUEL+KNUEL
			  NANUE=NANUE+KANUE
			  NNUMU=NNUMU+KNUMU
			  NANUM=NANUM+KANUM
			  GO TO 55
		      ELSE
			  CONTINUE
		      END IF
		  CALL SENDPR(J,KK,3)
		  GO TO 51
	      ELSE
		  CONTINUE
	      END IF
C
	      IF(JPART.EQ.13)THEN
C
C             --------------ELECTRON---------------------------------
		  KEMIN=KEMIN+1
		  KK=NEMIN+KEMIN-1
		      IF(KK.GE.L13100)THEN
			  LINE2=LINE2+1+KPRINU+KSECNU+KCASC+1
			  NNEUTR=NNEUTR+KNEUTR
			  NGAMMA=NGAMMA+KGAMMA
			  NEMIN=KK
			  NEPLS=NEPLS+KEPLS
			  NNUEL=NNUEL+KNUEL
			  NANUE=NANUE+KANUE
			  NNUMU=NNUMU+KNUMU
			  NANUM=NANUM+KANUM
			  GO TO 55
		      ELSE
			  CONTINUE
		      END IF
		  CALL SENDPR(J,KK,4)
		  GO TO 51
	      ELSE
		  CONTINUE
	      END IF
C
	      IF(JPART.EQ.14)THEN
C
C             --------------POSITRON---------------------------------
		  KEPLS=KEPLS+1
		  KK=NEPLS+KEPLS-1
		      IF(KK.GE.L14100)THEN
			  LINE2=LINE2+1+KPRINU+KSECNU+KCASC+1
			  NNEUTR=NNEUTR+KNEUTR
			  NGAMMA=NGAMMA+KGAMMA
			  NEMIN=NEMIN+KEMIN
			  NEPLS=KK
			  NNUEL=NNUEL+KNUEL
			  NANUE=NANUE+KANUE
			  NNUMU=NNUMU+KNUMU
			  NANUM=NANUM+KANUM
			  GO TO 55
		      ELSE
			  CONTINUE
		      END IF
		  CALL SENDPR(J,KK,5)
		  GO TO 51
	      ELSE
		  CONTINUE
	      END IF
C
	      IF(JPART.EQ.17)THEN
C
C             --------------ELECTRONIC NEITRINO---------------------
		  KNUEL=KNUEL+1
		  KK=NNUEL+KNUEL-1
		      IF(KK.GE.L17100)THEN
			  LINE2=LINE2+1+KPRINU+KSECNU+KCASC+1
			  NNEUTR=NNEUTR+KNEUTR
			  NGAMMA=NGAMMA+KGAMMA
			  NEMIN=NEMIN+KEMIN
			  NEPLS=NEPLS+KEPLS
			  NNUEL=KK
			  NANUE=NANUE+KANUE
			  NNUMU=NNUMU+KNUMU
			  NANUM=NANUM+KANUM
			  GO TO 55
		      ELSE
			  CONTINUE
		      END IF
		  CALL SENDPR(J,KK,6)
		  GO TO 51
	      ELSE
		  CONTINUE
	      END IF
C
	      IF(JPART.EQ.18)THEN
C
C             --------------ELECTRONIC ANTINEITRINO---------------------
		  KANUE=KANUE+1
		  KK=NANUE+KANUE-1
		      IF(KK.GE.L18100)THEN
			  LINE2=LINE2+1+KPRINU+KSECNU+KCASC+1
			  NNEUTR=NNEUTR+KNEUTR
			  NGAMMA=NGAMMA+KGAMMA
			  NEMIN=NEMIN+KEMIN
			  NEPLS=NEPLS+KEPLS
			  NNUEL=NNUEL+KNUEL
			  NANUE=KK
			  NNUMU=NNUMU+KNUMU
			  NANUM=NANUM+KANUM
			  GO TO 55
		      ELSE
			  CONTINUE
		      END IF
		  CALL SENDPR(J,KK,7)
		  GO TO 51
	      ELSE
		  CONTINUE
	      END IF
C
	      IF(JPART.EQ.19)THEN
C
C             --------------MUONIC NEITRINO--------------------------
		  KNUMU=KNUMU+1
		  KK=NNUMU+KNUMU-1
		      IF(KK.GE.L19100)THEN
			  LINE2=LINE2+1+KPRINU+KSECNU+KCASC+1
			  NNEUTR=NNEUTR+KNEUTR
			  NGAMMA=NGAMMA+KGAMMA
			  NEMIN=NEMIN+KEMIN
			  NEPLS=NEPLS+KEPLS
			  NNUEL=NNUEL+KNUEL
			  NANUE=NANUE+KANUE
			  NNUMU=KK
			  NANUM=NANUM+KANUM
			  GO TO 55
		      ELSE
			  CONTINUE
		      END IF
		  CALL SENDPR(J,KK,8)
		  GO TO 51
	      ELSE
		  CONTINUE
	      END IF
C
	      IF(JPART.EQ.20)THEN
C
C             --------------MUONIC ANTINEITRINO---------------------
		  KANUM=KANUM+1
		  KK=NANUM+KANUM-1
		      IF(KK.GE.L20100)THEN
			  LINE2=LINE2+1+KPRINU+KSECNU+KCASC+1
			  NNEUTR=NNEUTR+KNEUTR
			  NGAMMA=NGAMMA+KGAMMA
			  NEMIN=NEMIN+KEMIN
			  NEPLS=NEPLS+KEPLS
			  NNUEL=NNUEL+KNUEL
			  NANUE=NANUE+KANUE
			  NNUMU=NNUMU+KNUMU
			  NANUM=KK
			  GO TO 55
		      ELSE
			  CONTINUE
		      END IF
		  CALL SENDPR(J,KK,9)
		  GO TO 51
	      ELSE
		  CONTINUE
	      END IF
C
C-----------------CASCADE PARTICLE--------------------------------------------
		  KCASC=KCASC+1
                  KK=LINE2+1+KPRINU+KSECNU+KCASC
		      IF(KK.GE.LTR300)THEN
			  LINE2=KK
			  NNEUTR=NNEUTR+KNEUTR
			  NGAMMA=NGAMMA+KGAMMA
			  NEMIN=NEMIN+KEMIN
			  NEPLS=NEPLS+KEPLS
			  NNUEL=NNUEL+KNUEL
			  NANUE=NANUE+KANUE
			  NNUMU=NNUMU+KNUMU
			  NANUM=NANUM+KANUM
			  GO TO 55
		      ELSE
			  CONTINUE
		      END IF
		  CALL SENDPR(J,KK,1)
c                 No. of line containing the last cascade particle !HION1
                  if(KPRINU.eq.1)nprlas=LINE2+1+KPRINU+KSECNU+KCASC !HION1
		  GO TO 51
	  ELSE
C             PARTICLES IN ARRAY SPT ARE EXHAUSTED
	      GO TO 52
	  END IF
   51 CONTINUE
C
   52 CONTINUE                                                    !HION!
c Transformation of fast residual nuclei from NUCLEUS to PARTICLE !HION1
      if(KPRINU.eq.1)call nucpar(nucfir,nuclas,nprfir,nprlas,linadd) !HION1
c
C ESTIBLISHMENT OF THE NEW VALUE OF CONTROL VARIABLES  LINE2,NNEUTR,
C NGAMMA,NEMIN,NEPLS,NNUEL,NANUE,NNUMU,NANUM
      if(KPRINU.eq.1)then                                         !HION1
        LINE2=LINE2+1+KPRINU+KSECNU+KCASC+1+LINADD                !HION1
      else                                                        !HION1
        LINE2=LINE2+1+KPRINU+KSECNU+KCASC+1                       !HION1
      end if                                                      !HION1
      NNEUTR=NNEUTR+KNEUTR
      NGAMMA=NGAMMA+KGAMMA
      NEMIN=NEMIN+KEMIN
      NEPLS=NEPLS+KEPLS
      NNUEL=NNUEL+KNUEL
      NANUE=NANUE+KANUE
      NNUMU=NNUMU+KNUMU
      NANUM=NANUM+KANUM
C
      IF( LINE2.GE.LTR300-2 .OR. NNEUTR.GE.LSN100 .OR.
     *   NGAMMA.GE.LSG100   .OR. NEMIN. GE.L13100 .OR.
     *   NEPLS .GE.L14100   .OR. NNUEL. GE.L17100 .OR.
     *   NANUE .GE.L18100   .OR. NNUMU. GE.L19100 .OR.
     *   NANUM .GE.L20100)THEN
	  CONTINUE!WRITE(25,1001)IXINIT,NUMTRE,LINE2,NNEUTR,NGAMMA,
!     *                  NEMIN,NEPLS,NNUEL,NANUE,NNUMU,NANUM
 1001     FORMAT(' SEND DIAGNOSTIC: ARRAYS OVERFILLING, IXINIT=',I12,   ! GSITR
     *         ', NUMTRE=',I6/18X,'LINE2=',I6,', NNEUTR=',I6,
     *         ', NGAMMA=',I5,', NEMIN=',I4/18X,'NEPLS=',I4,', NNUEL=',
     *         I4,', NANUE=',I4,', NNUMU=',I4/18X,'NANUM=',I4)
	  IOVER=+1
      ELSE
	  IOVER=0
      END IF
      RETURN
C
   55 IOVER=+1
      CONTINUE!WRITE(25,1001)IXINIT,NUMTRE,LINE2,NNEUTR,NGAMMA,
!     *              NEMIN,NEPLS,NNUEL,NANUE,NNUMU,NANUM
      RETURN
C
      END



      SUBROUTINE KARSJH(KOUT)
C Debug printout of results of each elastic or inelastic interaction
C into file No. KOUT.
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      COMMON /INREAC/ COST,SINF,COSF,TINT,WINT,JPART,NUCLID,KSTATE
      COMMON /HIPROJ/ APROJ,ZPROJ
      COMMON /SECHAR/ SPT(6,5000),SNU(10,101),LS6,LS100,LS10,LS11  ! GSITR
C
      write(KOUT,10)NUMTRE,KSTATE
   10 format(/' KARSJH: NUMTRE=',I3,'    Interaction KSTATE=',I2)
      write(KOUT,20)JPART,TINT,APROJ,ZPROJ
   20 format(' Projectile: JPART=',I2,'  TMeV=',F8.3,' (Ap,Zp)=',2F4.0)
      write(KOUT,30)NUCLID
   30 format(' Target: NUCLID=',I3)
C
      IF(SPT(4,1).le.0.)THEN
        write(KOUT,41)
   41   format(' Secondary particles: NONE')
      ELSE
        write(KOUT,40)
   40   format(' Secondary particles:')
      END IF
      DO J=1,LS100
        IF(SPT(4,J).gt.0.)THEN
          write(KOUT,50)J,SPT(4,J),SPT(5,J),SPT(6,J)
   50     format(I3,'  TMeV=',F8.3,'  JPART=',F4.1,'  W=',F4.2)
        END IF	        
	END DO
C
      write(KOUT,60)
   60 format(' Nucleus-target and Residual nuclei:')
      DO J=1,LS11
        IF(SNU(4,J).ne.0.)THEN
          IF(SNU(1,J).ne.0.)THEN
            RM=SNU(1,J)*940.
            PX=SNU(5,J)
            PY=SNU(6,J)
            PZ=SNU(7,J)
            TNU=SQRT(PX**2+PY**2+PZ**2+RM**2)-RM
          ELSE
            TNU=0.
          END IF
          write(KOUT,70)J,SNU(1,J),SNU(2,J),SNU(3,J),SNU(4,J),TNU
   70     format(I3,'  A=',F4.0,'  Z=',F3.0,'  U*=',F7.2,
     *              '  KINTER=',F3.0,'  T(MeV)=',F8.3)
        END IF	        
	END DO
C
      RETURN
      END



      SUBROUTINE NUCPAR(NUCFIR,NUCLAS,NPRFIR,NPRLAS,LINADD)          !HION1
C Transformation of fast residual nuclei from NUCLEUS to PARTICLE
C INPUT:
C   Nuclei-products are placed from      line NUCFIR to line NUCLAS.
C   Secondary particles are placed from  line NPRFIR to line NPRLAS.
C OUTPUT:
C   All nuclei-products with kinetic energy > 3 MeV/A will be transformed
C   to "secondary particles" in the /TREE/ with JPART=25.
C
C   LINADD=1, if ALL nuclei-products were transformed to "particles",
C             so there is no line containing nucleus-product. In this
C             case the ADDITIONAL LINE should be inserted into /TREE/,
C             containing ZERO nucleus-product (as at FULL DESINTEGRATION).
C   LINADD=0 otherwise.
C ----------------------------------------------------------------------
C
C ----------------- COMMONs of GEMCA --------------------           ! GEMCA
      COMMON/GDATAP/ X,Y,Z,CX,CY,CZ,XN,YN,ZN,CXN,CYN,CZN,           ! GEMCA
     *               NZONO,MEDOLD,NBPO,NSO,                         ! GEMCA
     *               NZONC,MEDCUR,NBPC,NSCI,NSCO,                   ! GEMCA
     *               PINSFL,IOLEG                                   ! GEMCA
	REAL*8 X,Y,Z,CX,CY,CZ,XN,YN,ZN,CXN,CYN,CZN                  ! GEMCA
	INTEGER PINSFL                                              ! GEMCA
C
      COMMON /TREE/ TREE(15,400000),ITREE(7,400000),LTR13,LTR7,LTR300 ! GSITR
      DIMENSION TR0(15),ITR0(7)
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX  ! KarSjh (for debug only)
C
      LINADD=0                ! Initial state
      KSECNU=NUCLAS-NUCFIR+1  ! Number of nuclei-products
C
*&******Debug printout of the next branch before transformation*****
*&      CONTINUE!WRITE(25,301)NUMTRE
*&  301 format(/' NUMTRE=',I3)
*&      CONTINUE!WRITE(25,1)NUCFIR,NUCLAS,NPRFIR,NPRLAS,KSECNU
*&    1 format(' BRANCH: NUCFIR=',I3,'  NUCLAS=',I3,'  NPRFIR=',I3,
*&     *       '  NPRLAS=',I3,'     KSECNU=',I2/
*&     *        ' Before transformation')
*&        K=NUCFIR-3              ! Initiating particle before range
*&        L=NUCFIR-2              ! Initiating particle after range
*&        CONTINUE!WRITE(25,3)K,(ITREE(I,K),I=1,7),TREE(7,K)
*&        CONTINUE!WRITE(25,3)L,(ITREE(I,L),I=1,7),TREE(7,L)
*&C
*&      do J=NUCFIR-1,NUCLAS      ! Nuclei: target and products
*&C       ....... Kinetic energy of nucleus E/A .......
*&        A=TREE(4,J)
*&        IF(A.gt.0.)THEN
*&          AMAS=A*940.
*&          PX=TREE( 8,J)
*&          PY=TREE( 9,J)
*&          PZ=TREE(10,J)
*&          ETOT=SQRT(PX*PX+PY*PY+PZ*PZ+AMAS**2)-AMAS     ! Etot (MeV)
*&          EA=ETOT/A                                     ! E/A  (MeV)
*&        ELSE
*&          EA=0.
*&        END IF
*&C       .............................................
*&        CONTINUE!WRITE(25,2)J,(ITREE(I,J),I=1,7),TREE(4,J),TREE(5,J),EA
*&    2   format(I6,7I4,5X,2F5.1,F12.4)
*&      end do
*&C
*&      do J=NPRFIR,NPRLAS        ! Secondary CASCADE particles (NO n<14.5!)
*&        CONTINUE!WRITE(25,3)J,(ITREE(I,J),I=1,7),TREE(7,J)
*&    3   format(I6,7I4,15X,F9.1)
*&      end do
*&****************** end of debug printout *************************
C
C Check if this is nucleus or not?
      DO J=NUCFIR,NUCLAS
        IF(ITREE(1,J).ne.0 .or. ITREE(6,J).ne.1)THEN
          CONTINUE!WRITE(25,1000)
 1000     FORMAT(' NUCPAR DIAGNOSTIC: This is not a nucleus.')
          CONTINUE!WRITE(25,1001)(ITREE(I,J),I=1,7)
 1001     FORMAT(6X,7I4)
          CONTINUE!WRITE(25,1002)(TREE(I,J),I=1,7)
 1002     FORMAT(7F10.2)
c          STOP
        END IF
      END DO
C
C Loop over all nuclei-products to find and transform the fast nucleus.
      ITR=0                           ! Counter of tramsformations
      DO 10 J=NUCFIR,NUCLAS
        A=TREE(4,J)
          IF(A .lt. 1.0)GO TO 10      ! Full desintegration
        AMAS=A*940.
        ZION=TREE(5,J)
        W= TREE(7,J)
          PX=TREE( 8,J)
          PY=TREE( 9,J)
          PZ=TREE(10,J)
        ETOT=SQRT(PX*PX+PY*PY+PZ*PZ+AMAS**2)-AMAS     ! Etot (MeV)
        EA=ETOT/A                                     ! E/A  (MeV)
          IF(EA .le. 1.0)GO TO 10     ! Residual nucleus is "at rest" ! KARSJH
C
C       Fast nucleus should be transformed to "secondary particle".
        ITR=ITR+1                     ! Number of transformations
        ITREE(1,J)=25
        ITREE(6,J)=-1
        ITREE(7,J)=-1
C
        CALL DECARS(PX,PY,PZ,PMOD,COST,SINF,COSF)
        TREE( 4,J)=COST
        TREE( 5,J)=SINF
        TREE( 6,J)=COSF
        TREE( 7,J)=EA      ! E/A !!!
        TREE( 8,J)=W       ! Weight of particle = Weight of nucleus ?
C-------------------- ZONE DEFINITION --------- ! GEMCA
        X=DBLE(TREE(1,J))                       ! GEMCA
        Y=DBLE(TREE(2,J))                       ! GEMCA
        Z=DBLE(TREE(3,J))                       ! GEMCA
C                                               ! GEMCA
        COSTGM=TREE(4,J)                        ! GEMCA
        SINTGM=SQRT(1.0-COSTGM**2.)             ! GEMCA
        SINFGM=TREE(5,J)                        ! GEMCA
        COSFGM=TREE(6,J)                        ! GEMCA
C                                               ! GEMCA
        CX=DBLE(SINTGM*COSFGM)                  ! GEMCA
        CY=DBLE(SINTGM*SINFGM)                  ! GEMCA
        CZ=DBLE(COSTGM)                         ! GEMCA
C                                               ! GEMCA
        CALL GNEXTZ(NEXTZ)                      ! GEMCA
        TREE(9,J)=FLOAT(NEXTZ)                  ! GEMCA
C---------------------------------------------- ! GEMCA
        TREE(10,J)=0.
        TREE(11,J)=0.
        TREE(12,J)=0.
        TREE(13,J)=0.
        TREE(14,J)=A
        TREE(15,J)=ZION
   10 CONTINUE
      IF(ITR.eq.0)RETURN  ! No transformation was made
C
C
*&*******Debug printout of the next branch after transformation*****
*&      CONTINUE!WRITE(25,101)ITR
*&  101 format(' After transformation, ITR=',I2)
*&C
*&      do J=NUCFIR-1,NUCLAS      ! Nuclei: target and products
*&        IF(ITREE(1,J).eq.0)THEN
*&C         ....... Kinetic energy of nucleus E/A .......
*&          A=TREE(4,J)
*&          IF(A.gt.0.)THEN
*&            AMAS=A*940.
*&            PX=TREE( 8,J)
*&            PY=TREE( 9,J)
*&            PZ=TREE(10,J)
*&            ETOT=SQRT(PX*PX+PY*PY+PZ*PZ+AMAS**2)-AMAS     ! Etot (MeV)
*&            EA=ETOT/A                                     ! E/A  (MeV)
*&          ELSE
*&            EA=0.
*&          END IF
*&C         .............................................
*&          CONTINUE!WRITE(25,102)J,(ITREE(I,J),I=1,7),TREE(4,J),TREE(5,J),EA
*&  102     format(I6,7I4,5X,2F5.1,F12.4)
*&        ELSE
*&          CONTINUE!WRITE(25,104)J,(ITREE(I,J),I=1,7),TREE(7,J),
*&     *                    TREE(14,J),TREE(15,J)
*&  104     format(I6,7I4,15X,F9.1,2F5.1)
*&        END IF
*&      end do
*&C
*&      do J=NPRFIR,NPRLAS        ! Secondary CASCADE particles (NO n<14.5!)
*&        CONTINUE!WRITE(25,103)J,(ITREE(I,J),I=1,7),TREE(7,J)
*&  103   format(I6,7I4,15X,F9.1)
*&      end do
*&********************** end of debug printout **********************
C
C If ALL nuclei-products was transformed into "particles",
C let us insert ZERO nucleus-product:
      IF(KSECNU.eq.ITR)THEN
        LINADD=1
C       ............. ZERO nucleus ...............
        DO I=1,LTR7
          ITR0(I)=ITREE(I,NUCFIR-1) ! Copy nucleus-target
        END DO
        ITR0(6)=1               ! KBEAFT=1
        ITR0(7)=6               ! KINTER=6 "Full desintegration"
C
        DO I=1,LTR13
          TR0(I)=TREE(I,NUCFIR-1)   ! Copy nucleus-target
        END DO
        TR0(4)=0.               ! A*=0
        TR0(5)=0.               ! Z*=0
        TR0(6)=0.               ! U*=0
        TR0(7)=TREE(8,NUCFIR-2) ! Weight of incident particle
C       ...........................................
C
C       Determine the FIRST and the LAST LINES of fragment of /TREE/
C       which should be shifted: LINF and LINL correspondingly.
        LINF=NUCFIR       ! First line
        IF(NPRFIR.le.NPRLAS)THEN    ! Are there any cascade particles?
          LINL=NPRLAS     ! Last line: There is at least ONE casc. part.
        ELSE
          LINL=NUCLAS     ! Last line: There is NO cascade particles
        END IF
C
C       Shifting of fragment (LINF,LINL) to (LINF+1,LINL+1):
        DO J=LINL,LINF,-1
          DO I=1,LTR7
            ITREE(I,J+1)=ITREE(I,J)
          END DO
          DO I=1,LTR13
            TREE(I,J+1)=TREE(I,J)
          END DO
        END DO
C
C       Insert ZERO nucleus into line NUCFIR=LINF:
        DO I=1,LTR7
          ITREE(I,NUCFIR)=ITR0(I)
        END DO
        DO I=1,LTR13
          TREE(I,NUCFIR)=TR0(I)
        END DO
C
        GO TO 30
      END IF
C
C
C If some (NOT ALL!) nuclei-products was transformed into "particles"
C but some ones wasn't, let us make permutation in /TREE/ to collect
C together the nuclei-products beginning from LINE No. NUCFIR (and to
C collect together the nuclei-"particles" following nuclei-products).
C
      DO 20 J=NUCFIR,NUCLAS
        IF(ITREE(1,J).eq.0)THEN
          GO TO 20    ! No permutation required
        ELSE
          DO 21 K=J+1,NUCLAS
            IF(ITREE(1,K).eq.0)THEN  ! To find the nearest nuc.-prod.
C             Permutation of lines No. J and K --------------
              DO I=1,LTR7           !
                ITR0(I)=ITREE(I,K)  !
              END DO                ! Save line No.K
              DO I=1,LTR13          !
                TR0(I)=TREE(I,K)    !
              END DO                !
C
              DO I=1,LTR7              !
                ITREE(I,K)=ITREE(I,J)  !
              END DO                   ! Recording of line No.J
              DO I=1,LTR13             ! into line No.K
                TREE(I,K)=TREE(I,J)    !
              END DO                   !
C
              DO I=1,LTR7           !
                ITREE(I,J)=ITR0(I)  !
              END DO                ! Recording of line No.K
              DO I=1,LTR13          ! into line No.J
                TREE(I,J)=TR0(I)    !
              END DO                !
C             End of permutation ----------------------------
              GO TO 20  ! To the next J
            ELSE
              GO TO 21
            END IF
   21     CONTINUE
        END IF
   20 CONTINUE
C
C
   30 CONTINUE
C
*&*******Debug printout of the next branch in final view*****
*&      CONTINUE!WRITE(25,201)
*&  201 format(' Final view')
*&C
*&      do J=NUCFIR-1,NUCLAS      ! Nuclei: target and products
*&        IF(ITREE(1,J).eq.0)THEN
*&C         ....... Kinetic energy of nucleus E/A .......
*&          A=TREE(4,J)
*&          IF(A.gt.0.)THEN
*&            AMAS=A*940.
*&            PX=TREE( 8,J)
*&            PY=TREE( 9,J)
*&            PZ=TREE(10,J)
*&            ETOT=SQRT(PX*PX+PY*PY+PZ*PZ+AMAS**2)-AMAS     ! Etot (MeV)
*&            EA=ETOT/A                                     ! E/A  (MeV)
*&          ELSE
*&            EA=0.
*&          END IF
*&C         .............................................
*&          CONTINUE!WRITE(25,202)J,(ITREE(I,J),I=1,7),TREE(4,J),TREE(5,J),EA
*&  202     format(I6,7I4,5X,2F5.1,F12.4)
*&        ELSE
*&          CONTINUE!WRITE(25,204)J,(ITREE(I,J),I=1,7),TREE(7,J),
*&     *                    TREE(14,J),TREE(15,J)
*&  204     format(I6,7I4,15X,F9.1,2F5.1)
*&        END IF
*&      end do
*&C
*&      do J=NPRFIR,NPRLAS+1      ! Secondary CASCADE particles (NO n<14.5!)
*&        if(ITREE(1,J).gt.0)then
*&          if(ITREE(1,J).eq.25)then
*&            CONTINUE!WRITE(25,204)J,(ITREE(I,J),I=1,7),TREE(7,J),
*&     *                      TREE(14,J),TREE(15,J)
*&          else
*&            CONTINUE!WRITE(25,203)J,(ITREE(I,J),I=1,7),TREE(7,J)
*&  203       format(I6,7I4,15X,F9.1)
*&          end if
*&        end if
*&      end do
*&********************** end of debug printout **********************
C
      RETURN
      END


      SUBROUTINE SENDNU(J,K)
C
C SENDING OF LINE NUMBER   J   OF ARRAY SNU(10,J)
C INTO LINE NUMBER   K   OF ARRAYS TREE(15,K),ITREE(7,K)             !HION1
C
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      COMMON /TREE/ TREE(15,400000),ITREE(7,400000),LTR13,LTR7,LTR300 ! GSITR
      COMMON /SECHAR/ SPT(6,5000),SNU(10,101),LS6,LS100,LS10,LS11  ! GSITR
      COMMON /MANALL/ NTREE,LINE1,LINE2,NBRMAX,NPRMAX,NNEUTR,NGAMMA,
     *                MELSC,MINTER,MCUT,MOUT,MABSOR,MDECAY,NSTRAN,
     *                NEMIN,NEPLS,NNUEL,NANUE,NNUMU,NANUM
      COMMON /INGEN/ INGEN                        ! KarSjh4
C
C COORDINATES X,Y,Z
      TREE(1,K)=TREE(1,LINE2+1)
      TREE(2,K)=TREE(2,LINE2+1)
      TREE(3,K)=TREE(3,LINE2+1)
C VALUES A,Z,U
      TREE(4,K)=SNU(1,J)
      TREE(5,K)=SNU(2,J)
      TREE(6,K)=SNU(3,J)
C WEIGHT OF NUCLEUS
      IF(J.EQ.1)THEN
	  TREE(7,K)=1.
      ELSE
	  TREE(7,K)=TREE(8,LINE2+1)
      END IF
C MOMENTUM AND ANGULAR MOMENTUM OF NUCLEUS
      DO 10 I=8,13                                             !HION1
   10 TREE(I,K)=SNU(I-3,J) 
C
C JPART      
      ITREE(1,K)=0
C NTREE
      ITREE(2,K)=ITREE(2,LINE2+1)
C NBRAN
      ITREE(3,K)=ITREE(3,LINE2+1)
C KGENUC
C=========== KarSjh4 =================
      IF(INGEN.eq.1)THEN
        ITREE(5,K)=ITREE(5,LINE2+1)
      ELSE
        ITREE(5,K)=ITREE(5,LINE2+1)+1
      END IF
C========= end KarSjh4 ===============
C KBEAFT
      IF(J.EQ.1)THEN      
	  ITREE(6,K)=0
      ELSE
	  ITREE(6,K)=1
      END IF 
C KINTER
      ITREE(7,K)=IFIX(SNU(4,J))
C
      RETURN
      END
 
 
      SUBROUTINE SENDPR(J,K,LLL)
C
C SENDING OF LINE NUMBER   J   OF ARRAY SPT(6,J)
C INTO LINE NUMBER   K   OF ARRAYS:     TREE(15,K),ITREE(7,K)      !HION1
C                                   OR SOURN(13,K),ISOURN(7,K)
C                                   OR SOURG(13,K),ISOURG(7,K)
C                                   OR SEMIN(13,K),ISEMIN(7,K)
C                                   OR SEPLS(13,K),ISEPLS(7,K)
C                                   OR SNUEL(13,K),ISNUEL(7,K)
C                                   OR SANUE(13,K),ISANUE(7,K)
C                                   OR SNUMU(13,K),ISNUMU(7,K)
C                                   OR SANUM(13,K),ISANUM(7,K)
C
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
      COMMON /SECHAR/ SPT(6,5000),SNU(10,101),LS6,LS100,LS10,LS11  ! GSITR
      COMMON /MANALL/ NTREE,LINE1,LINE2,NBRMAX,NPRMAX,NNEUTR,NGAMMA,
     *                MELSC,MINTER,MCUT,MOUT,MABSOR,MDECAY,NSTRAN,
     *                NEMIN,NEPLS,NNUEL,NANUE,NNUMU,NANUM
      COMMON /INGEN/ INGEN                        ! KarSjh4
C
      GO TO (100,200,300,400,500,600,700,800,900),LLL 
C
C SENDING INTO ARRAYS TREE,ITREE---------------------------
  100 CONTINUE
C COORDINATES X,Y,Z
      TREE(1,K)=TREE(1,LINE2+1)
      TREE(2,K)=TREE(2,LINE2+1)
      TREE(3,K)=TREE(3,LINE2+1)
C DIRECTION OF FLY
      TREE(4,K)=SPT(1,J)
      TREE(5,K)=SPT(2,J)
      TREE(6,K)=SPT(3,J)
C KINETIC ENERGY OF RARTICLE
      TREE(7,K)=SPT(4,J)
C WEIGHT OF PARTICLE
      TREE(8,K)=TREE(8,LINE2+1)*SPT(6,J)      
C GEOMETRY ZONE                                                      ! GEMCA
      TREE(9,K)=TREE(9,LINE2+1)                                      ! GEMCA
C
C JPART      
      ITREE(1,K)=IFIX(SPT(5,J))
C NTREE
      ITREE(2,K)=ITREE(2,LINE2+1)
C NBRAN
      ITREE(3,K)=ITREE(3,LINE2+1)
C KGENER
C=========== KarSjh4 =================
      IF(INGEN.eq.1)THEN
        ITREE(5,K)=ITREE(5,LINE2+1)
      ELSE
        ITREE(5,K)=ITREE(5,LINE2+1)+1
      END IF
C========= end KarSjh4 ===============
C KPASS
      ITREE(6,K)=-1
C KINTER
      ITREE(7,K)=-1
C ********************* !HION7 **********************
C Transformation MeV --> MeV/A for D, T, He3, He4, HI
      IF(ITREE(1,K).eq.21)TREE(7,K)=TREE(7,K)/2.0
      IF(ITREE(1,K).eq.22)TREE(7,K)=TREE(7,K)/3.0
      IF(ITREE(1,K).eq.23)TREE(7,K)=TREE(7,K)/3.0
      IF(ITREE(1,K).eq.24)TREE(7,K)=TREE(7,K)/4.0
      IF(ITREE(1,K).eq.25)TREE(7,K)=TREE(7,K)/TREE(14,K)
C ******************* end of !HION7 ******************
C
      RETURN
C
C SENDING INTO ARRAYS SOURN,ISOURN-----------------------------
  200 CONTINUE
C COORDINATES X,Y,Z
      SOURN(1,K)=TREE(1,LINE2+1)
      SOURN(2,K)=TREE(2,LINE2+1)
      SOURN(3,K)=TREE(3,LINE2+1)
C DIRECTION OF FLY
      SOURN(4,K)=SPT(1,J)
      SOURN(5,K)=SPT(2,J)
      SOURN(6,K)=SPT(3,J)
C KINETIC ENERGY OF RARTICLE
      SOURN(7,K)=SPT(4,J)
C WEIGHT OF PARTICLE
      SOURN(8,K)=TREE(8,LINE2+1)*SPT(6,J)      
C GEOMETRY ZONE                                                      ! GEMCA
      SOURN(9,K)=TREE(9,LINE2+1)                                     ! GEMCA
C
C JPART      
      ISOURN(1,K)=IFIX(SPT(5,J))
C NTREE
      ISOURN(2,K)=ITREE(2,LINE2+1)
C NBRAN
      ISOURN(3,K)=ITREE(3,LINE2+1)
C KGENER
C=========== KarSjh4 =================
      IF(INGEN.eq.1)THEN
        ISOURN(5,K)=ITREE(5,LINE2+1)
      ELSE
        ISOURN(5,K)=ITREE(5,LINE2+1)+1
      END IF
C========= end KarSjh4 ===============
C KPASS
      ISOURN(6,K)=-1
C KINTER
      ISOURN(7,K)=-1
C
      RETURN
C
C SENDING INTO ARRAYS SOURG,ISOURG------------------------------
  300 CONTINUE
C COORDINATES X,Y,Z
      SOURG(1,K)=TREE(1,LINE2+1)
      SOURG(2,K)=TREE(2,LINE2+1)
      SOURG(3,K)=TREE(3,LINE2+1)
C DIRECTION OF FLY
      SOURG(4,K)=SPT(1,J)
      SOURG(5,K)=SPT(2,J)
      SOURG(6,K)=SPT(3,J)
C KINETIC ENERGY OF RARTICLE
      SOURG(7,K)=SPT(4,J)
C WEIGHT OF PARTICLE
      SOURG(8,K)=TREE(8,LINE2+1)*SPT(6,J)      
C GEOMETRY ZONE                                                      ! GEMCA
      SOURG(9,K)=TREE(9,LINE2+1)                                     ! GEMCA
C
C JPART      
      ISOURG(1,K)=IFIX(SPT(5,J))
C NTREE
      ISOURG(2,K)=ITREE(2,LINE2+1)
C NBRAN
      ISOURG(3,K)=ITREE(3,LINE2+1)
C KGENER
C=========== KarSjh4 =================
      IF(INGEN.eq.1)THEN
        ISOURG(5,K)=ITREE(5,LINE2+1)
      ELSE
        ISOURG(5,K)=ITREE(5,LINE2+1)+1
      END IF
C========= end KarSjh4 ===============
C KPASS
      ISOURG(6,K)=-1
C KINTER
      ISOURG(7,K)=-1
C
      RETURN
C
C SENDING INTO ARRAYS SEMIN,ISEMIN---------------------------
  400 CONTINUE
C COORDINATES X,Y,Z
      SEMIN(1,K)=TREE(1,LINE2+1)
      SEMIN(2,K)=TREE(2,LINE2+1)
      SEMIN(3,K)=TREE(3,LINE2+1)
C DIRECTION OF FLY
      SEMIN(4,K)=SPT(1,J)
      SEMIN(5,K)=SPT(2,J)
      SEMIN(6,K)=SPT(3,J)
C KINETIC ENERGY OF RARTICLE
      SEMIN(7,K)=SPT(4,J)
C WEIGHT OF PARTICLE
      SEMIN(8,K)=TREE(8,LINE2+1)*SPT(6,J)      
C GEOMETRY ZONE                                                      ! GEMCA
      SEMIN(9,K)=TREE(9,LINE2+1)                                     ! GEMCA
C
C JPART      
      ISEMIN(1,K)=IFIX(SPT(5,J))
C NTREE
      ISEMIN(2,K)=ITREE(2,LINE2+1)
C NBRAN
      ISEMIN(3,K)=ITREE(3,LINE2+1)
C KGENER
C=========== KarSjh4 =================
      IF(INGEN.eq.1)THEN
        ISEMIN(5,K)=ITREE(5,LINE2+1)
      ELSE
        ISEMIN(5,K)=ITREE(5,LINE2+1)+1
      END IF
C========= end KarSjh4 ===============
C KPASS
      ISEMIN(6,K)=-1
C KINTER
      ISEMIN(7,K)=-1
C
      RETURN
C
C SENDING INTO ARRAYS SEPLS,ISEPLS------------------------------
  500 CONTINUE
C COORDINATES X,Y,Z
      SEPLS(1,K)=TREE(1,LINE2+1)
      SEPLS(2,K)=TREE(2,LINE2+1)
      SEPLS(3,K)=TREE(3,LINE2+1)
C DIRECTION OF FLY
      SEPLS(4,K)=SPT(1,J)
      SEPLS(5,K)=SPT(2,J)
      SEPLS(6,K)=SPT(3,J)
C KINETIC ENERGY OF RARTICLE
      SEPLS(7,K)=SPT(4,J)
C WEIGHT OF PARTICLE
      SEPLS(8,K)=TREE(8,LINE2+1)*SPT(6,J)      
C GEOMETRY ZONE                                                      ! GEMCA
      SEPLS(9,K)=TREE(9,LINE2+1)                                     ! GEMCA
C
C JPART      
      ISEPLS(1,K)=IFIX(SPT(5,J))
C NTREE
      ISEPLS(2,K)=ITREE(2,LINE2+1)
C NBRAN
      ISEPLS(3,K)=ITREE(3,LINE2+1)
C KGENER
C=========== KarSjh4 =================
      IF(INGEN.eq.1)THEN
        ISEPLS(5,K)=ITREE(5,LINE2+1)
      ELSE
        ISEPLS(5,K)=ITREE(5,LINE2+1)+1
      END IF
C========= end KarSjh4 ===============
C KPASS
      ISEPLS(6,K)=-1
C KINTER
      ISEPLS(7,K)=-1
C
      RETURN
C
C SENDING INTO ARRAYS SNUEL,ISNUEL--------------------------------
  600 CONTINUE
C COORDINATES X,Y,Z
      SNUEL(1,K)=TREE(1,LINE2+1)
      SNUEL(2,K)=TREE(2,LINE2+1)
      SNUEL(3,K)=TREE(3,LINE2+1)
C DIRECTION OF FLY
      SNUEL(4,K)=SPT(1,J)
      SNUEL(5,K)=SPT(2,J)
      SNUEL(6,K)=SPT(3,J)
C KINETIC ENERGY OF RARTICLE
      SNUEL(7,K)=SPT(4,J)
C WEIGHT OF PARTICLE
      SNUEL(8,K)=TREE(8,LINE2+1)*SPT(6,J)      
C GEOMETRY ZONE                                                      ! GEMCA
      SNUEL(9,K)=TREE(9,LINE2+1)                                     ! GEMCA
C
C JPART      
      ISNUEL(1,K)=IFIX(SPT(5,J))
C NTREE
      ISNUEL(2,K)=ITREE(2,LINE2+1)
C NBRAN
      ISNUEL(3,K)=ITREE(3,LINE2+1)
C KGENER
C=========== KarSjh4 =================
      IF(INGEN.eq.1)THEN
        ISNUEL(5,K)=ITREE(5,LINE2+1)
      ELSE
        ISNUEL(5,K)=ITREE(5,LINE2+1)+1
      END IF
C========= end KarSjh4 ===============
C KPASS
      ISNUEL(6,K)=-1
C KINTER
      ISNUEL(7,K)=-1
C
      RETURN
C
C SENDING INTO ARRAYS SANUE,ISANUE----------------------------
  700 CONTINUE
C COORDINATES X,Y,Z
      SANUE(1,K)=TREE(1,LINE2+1)
      SANUE(2,K)=TREE(2,LINE2+1)
      SANUE(3,K)=TREE(3,LINE2+1)
C DIRECTION OF FLY
      SANUE(4,K)=SPT(1,J)
      SANUE(5,K)=SPT(2,J)
      SANUE(6,K)=SPT(3,J)
C KINETIC ENERGY OF RARTICLE
      SANUE(7,K)=SPT(4,J)
C WEIGHT OF PARTICLE
      SANUE(8,K)=TREE(8,LINE2+1)*SPT(6,J)      
C GEOMETRY ZONE                                                      ! GEMCA
      SANUE(9,K)=TREE(9,LINE2+1)                                     ! GEMCA
C
C JPART      
      ISANUE(1,K)=IFIX(SPT(5,J))
C NTREE
      ISANUE(2,K)=ITREE(2,LINE2+1)
C NBRAN
      ISANUE(3,K)=ITREE(3,LINE2+1)
C KGENER
C=========== KarSjh4 =================
      IF(INGEN.eq.1)THEN
        ISANUE(5,K)=ITREE(5,LINE2+1)
      ELSE
        ISANUE(5,K)=ITREE(5,LINE2+1)+1
      END IF
C========= end KarSjh4 ===============
C KPASS
      ISANUE(6,K)=-1
C KINTER
      ISANUE(7,K)=-1
C
      RETURN
C
C SENDING INTO ARRAYS SNUMU,ISNUMU-------------------------------
  800 CONTINUE
C COORDINATES X,Y,Z
      SNUMU(1,K)=TREE(1,LINE2+1)
      SNUMU(2,K)=TREE(2,LINE2+1)
      SNUMU(3,K)=TREE(3,LINE2+1)
C DIRECTION OF FLY
      SNUMU(4,K)=SPT(1,J)
      SNUMU(5,K)=SPT(2,J)
      SNUMU(6,K)=SPT(3,J)
C KINETIC ENERGY OF RARTICLE
      SNUMU(7,K)=SPT(4,J)
C WEIGHT OF PARTICLE
      SNUMU(8,K)=TREE(8,LINE2+1)*SPT(6,J)      
C GEOMETRY ZONE                                                      ! GEMCA
      SNUMU(9,K)=TREE(9,LINE2+1)                                     ! GEMCA
C
C JPART      
      ISNUMU(1,K)=IFIX(SPT(5,J))
C NTREE
      ISNUMU(2,K)=ITREE(2,LINE2+1)
C NBRAN
      ISNUMU(3,K)=ITREE(3,LINE2+1)
C KGENER
C=========== KarSjh4 =================
      IF(INGEN.eq.1)THEN
        ISNUMU(5,K)=ITREE(5,LINE2+1)
      ELSE
        ISNUMU(5,K)=ITREE(5,LINE2+1)+1
      END IF
C========= end KarSjh4 ===============
C KPASS
      ISNUMU(6,K)=-1
C KINTER
      ISNUMU(7,K)=-1
C
      RETURN
C
C SENDING INTO ARRAYS SANUM,ISANUM---------------------------------
  900 CONTINUE
C COORDINATES X,Y,Z
      SANUM(1,K)=TREE(1,LINE2+1)
      SANUM(2,K)=TREE(2,LINE2+1)
      SANUM(3,K)=TREE(3,LINE2+1)
C DIRECTION OF FLY
      SANUM(4,K)=SPT(1,J)
      SANUM(5,K)=SPT(2,J)
      SANUM(6,K)=SPT(3,J)
C KINETIC ENERGY OF RARTICLE
      SANUM(7,K)=SPT(4,J)
C WEIGHT OF PARTICLE
      SANUM(8,K)=TREE(8,LINE2+1)*SPT(6,J)      
C GEOMETRY ZONE                                                      ! GEMCA
      SANUM(9,K)=TREE(9,LINE2+1)                                     ! GEMCA
C
C JPART      
      ISANUM(1,K)=IFIX(SPT(5,J))
C NTREE
      ISANUM(2,K)=ITREE(2,LINE2+1)
C NBRAN
      ISANUM(3,K)=ITREE(3,LINE2+1)
C KGENER
C=========== KarSjh4 =================
      IF(INGEN.eq.1)THEN
        ISANUM(5,K)=ITREE(5,LINE2+1)
      ELSE
        ISANUM(5,K)=ITREE(5,LINE2+1)+1
      END IF
C========= end KarSjh4 ===============
C KPASS
      ISANUM(6,K)=-1
C KINTER
      ISANUM(7,K)=-1
C
      RETURN
C
      END
 
 
 
      SUBROUTINE SENS(KOVER)
C SENDING OF THE STEPS OF RANGES FROM ARRAY  STEPR INTO ARRAYS  TRAN,ITRAN
C     OUTPUT VARIABLE  KOVER=+1 IF ARRAYS  TRAN,ITRAN ARE OVERFILLED 
C     ELSE  KOVER=0
C
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      COMMON /MANALL/ NTREE,LINE1,LINE2,NBRMAX,NPRMAX,NNEUTR,NGAMMA,
     *                MELSC,MINTER,MCUT,MOUT,MABSOR,MDECAY,NSTRAN,
     *                NEMIN,NEPLS,NNUEL,NANUE,NNUMU,NANUM
      COMMON /TRAN/ TRAN(5,600000),ITRAN(6,600000),LSR5,LSR6,LSR500  ! GSITR
      COMMON /STERAN/ STEPR(9,500),LST9,LST50
      COMMON /TREE/ TREE(15,400000),ITREE(7,400000),LTR13,LTR7,LTR300 ! GSITR
C
      ISTEP=IFIX(STEPR(2,1))
      IF(NSTRAN+ISTEP .GT. LSR500-1)THEN
	  CONTINUE!WRITE(25,1001)IXINIT,NUMTRE
 1001     FORMAT(' SENS DIAGNOSTIC: ARRAY TRAN EXCEEDED, IXINIT=',I12,
     *           ', NUMTRE=',I6)
	  KOVER=1 
	  RETURN
      ELSE
	  KOVER=0
      END IF
C
      DO 10 J=1,ISTEP+1
	    KK=NSTRAN+J-1
      TRAN(1,KK)=STEPR(4,J)
      TRAN(2,KK)=STEPR(6,J)
      TRAN(3,KK)=STEPR(7,J)
      TRAN(4,KK)=STEPR(8,J)
      TRAN(5,KK)=STEPR(9,J)
C
      ITRAN(1,KK)=ITREE(1,LINE2)
      ITRAN(2,KK)=ITREE(2,LINE2)
      ITRAN(3,KK)=ITREE(3,LINE2)
      ITRAN(4,KK)=IFIX(STEPR(1,J))
      ITRAN(5,KK)=IFIX(STEPR(2,J))
      ITRAN(6,KK)=IFIX(STEPR(3,J))
   10 CONTINUE
C
	  NSTRAN=NSTRAN+ISTEP+1
      RETURN
      END



      SUBROUTINE PREPAR
C SUMPLING OF THE NUCLEUS-TARGET NUMBER  NUCLID  AND TYPE OF INTERACTION
C (ELASTIC OR INELASTIC)  KSTATE FOR PARTICLE TYPE  JPART  IN THE END-POINTE
C OF IT RANGE
C INPUT:
C       LINE NO.LINE2+1 OF ARRAYS  TREE,ITREE
C       LINE NO.NSTRAN-1 OF ARRAYS  TRAN,ITRAN
C OUTPUT:
C        COST,SINF,COSF - PARTICLE-PROJECTILE DIRECTION
C        TINT - PARTICLE-PROJECTILE ENERGY (MEV)
C        WINT - PARTICLE-PROJECTILE WEIGHT
C        JPART - PARTICLE-PROJECTILE TYPE
C        NUCLID - NO. OF NUCLEUS-TARGET
C        KSTATE - INTERACTION TYPE: 1 - ELASTIC SCATTERING
C                                   2 - INELASTIC INTERACTION
C                                   6 - PI- and K- decay at absor-
C                                       btion by H1
C---------------------------------------------------------------
      COMMON /RANDOM/ IX,IXINIT,NUMTRE,IXFIRS,LASTIX
      COMMON /TREE/ TREE(15,400000),ITREE(7,400000),LTR13,LTR7,LTR300 ! GSITR
      COMMON /MANALL/ NTREE,LINE1,LINE2,NBRMAX,NPRMAX,NNEUTR,NGAMMA,
     *                MELSC,MINTER,MCUT,MOUT,MABSOR,MDECAY,NSTRAN,
     *                NEMIN,NEPLS,NNUEL,NANUE,NNUMU,NANUM
      COMMON /TRAN/ TRAN(5,600000),ITRAN(6,600000),LSR5,LSR6,LSR500  ! GSITR
      COMMON /MACRIN/ ELMIX(8,24,48),RHOMED(48),NELEMD(2,48),NUMMED  ! CORNEL
      COMMON /INREAC/ COST,SINF,COSF,TINT,WINT,JPART,NUCLID,KSTATE
	  DIMENSION TEMP(24)  ! CORNEL
C
      COST=TREE(4,LINE2+1)
      SINF=TREE(5,LINE2+1)
      COSF=TREE(6,LINE2+1)
          TINT=TREE(7,LINE2+1)
	  WINT=TREE(8,LINE2+1)
		  JPART=ITREE(1,LINE2+1)
		      MEDIA=ITRAN(5,NSTRAN-1)
			  KSTATE=ITREE(7,LINE2+1)
C ******************* !HION10 *******************
      APART=1.0
      IF(JPART.eq.21)APART=2.0
      IF(JPART.eq.22)APART=3.0
      IF(JPART.eq.23)APART=3.0
      IF(JPART.eq.24)APART=4.0
      IF(JPART.eq.25)APART=TREE(14,LINE2+1)
C ************** end of !HION10 *****************
C
      IF((JPART.EQ.3 .OR. JPART.EQ.7 .OR. JPART.EQ.8) .AND. 
     *   KSTATE.EQ.2)GO TO 20
C
C SAMPLING OF NUCLEUS-TARGET NUMBER NUCLID
      NELEM=NELEMD(2,MEDIA)
	  IF(NELEM.EQ.1)THEN
	      NUCLID=IFIX(ELMIX(1,1,MEDIA))
	      GO TO 11
	  ELSE   
	      CONTINUE
	  END IF
C
      DO 1 I=1,NELEM
      NUCLID=IFIX(ELMIX(1,I,MEDIA))
      CONC  =ELMIX(2,I,MEDIA)
      CALL MICROD(JPART,TINT*APART,NUCLID,SIGTOT,SIGIN)
      TEMP(I)=CONC*SIGTOT
    1 CONTINUE
C
      DO 2 I=2,NELEM
    2 TEMP(I)=TEMP(I-1)+TEMP(I)
C
      BETHA=RRAN(IX)*TEMP(NELEM)
      DO 3 I=1,NELEM
	  IF(BETHA.LE.TEMP(I))THEN
	      NUCLID=IFIX(ELMIX(1,I,MEDIA))
	      GO TO 11
	  ELSE
	      GO TO 3
	  END IF
    3 CONTINUE
C
C SAMPLING OF INTERACTION TYPE   KSTATE
   11 IF(NUCLID.EQ.1)THEN
	  KSTATE=2
	  ITREE(7,LINE2+1)=KSTATE
	  RETURN
      ELSE
	  GO TO 12
      END IF
C
   12 CALL MICROD(JPART,TINT*APART,NUCLID,SIGTOT,SIGIN)
	  IF(RRAN(IX).LE.SIGIN/SIGTOT)THEN
	      KSTATE=2
	  ELSE
	      KSTATE=1
	  END IF
      ITREE(7,LINE2+1)=KSTATE
      RETURN
C
C PI-, K- and antiproton ABSORPTION AT REST
   20 NELEM=NELEMD(2,MEDIA)
	  IF(NELEM.EQ.1)THEN
	      NUCLID=IFIX(ELMIX(1,1,MEDIA))
	      TINT=0.001
	      GO TO 25
	  ELSE   
	      CONTINUE
	  END IF
C
C FERMI-TELLER Z-RULE: interaction probability is proportional to Zi*Ni
      DO 21 I=1,NELEM
   21 TEMP(I)=ELMIX(2,I,MEDIA)*ELMIX(4,I,MEDIA)
C
      DO 22 I=2,NELEM
   22 TEMP(I)=TEMP(I-1)+TEMP(I)
C
      BETHA=RRAN(IX)*TEMP(NELEM)
      DO 23 I=1,NELEM
	  IF(BETHA.LE.TEMP(I))THEN
	      NUCLID=IFIX(ELMIX(1,I,MEDIA))
	      TINT=0.001
	      GO TO 25
	  ELSE
	      GO TO 23
	  END IF
   23 CONTINUE
C
   25 IF((JPART.EQ.3 .OR. JPART.EQ.8) .AND. NUCLID.EQ.1)THEN
	  KSTATE=6
	  ITREE(7,LINE2+1)=KSTATE
	  RETURN
      ELSE
	  RETURN
      END IF
C
      END
