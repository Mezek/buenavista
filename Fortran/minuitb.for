************************************************************************
**                                                                    **
**                           F.James DD/CERN                          **
**                                                                    **
**             Function Minimisation and Error Analysis.              **
**                CERN program library MINUIT (D506).                 **
**                                                                    **
**                  Source Version 1.04, 9-May-1983.                  **
**                     Adapted for IBM PC, 1988.                      **
**                                                                    **
************************************************************************

     	
      SUBROUTINE MINUIT
      CALL MINTSD
      RETURN
      END
		


      SUBROUTINE MINTSD


**
**       This is the main program.  It initializes some COMMON constants
**       then verifies that  FCN  gives the same value when called twice
**       with the same arguments and passes control to COMAND.
**

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)


      CHARACTER*60 TITLE
      CHARACTER*10 CWORD
      CHARACTER*10 NAM

      COMMON
     . /MINERR/ ERP  (100),ERN  (100)
     . /DERIVA/ G    (100),G2   (100)
     . /LIMITS/ ALIM (100),BLIM (100),LCODE (100),LCORSP(100) ,LIMSET
     . /PAREXT/ U    (100),NAM  (100),WERR  (100),MAXEXT     ,NU
     . /PARINT/ X (50)   ,XT   (50),DIRIN (50),MAXINT     ,NPAR
     . /FIX   / XS(50)   ,XTS  (50),WTS   (50),IPFIX(50)  ,NPFIX
     . /SIMVEC/ P (50,51),PSTAR(50),PSTST (50),PBAR (50)  ,PRHO(50)
     . /VARIAN/ V (50,50)
     . /VARIAT/ VT(50,50)
     . /CASC  / JH,JL,Y(51)

      COMMON
     . /UNIT  / ISYSRD,ISYSWR  ,ISYSPU
     . /TITLE / TITLE ,ISW  (7),NBLOCK
     . /CARD  / CWORD ,WORD7(7)
     . /CONVER/ EPSI  ,APSI,VTEST ,NSTEPQ,NFCN,NFCNMX
     . /MINIMA/ AMIN  ,UP  ,NEWMIN,ITAUR ,SIGMA


      MAXINT = 50
      MAXEXT = 100

  100 NFCN = 1
      CALL MIDATA(0)
      CALL INTOEX(X)
      WRITE (ISYSWR,1001)
      CALL FCN(NPAR,G,AMIN,U,1)
      CALL FCN(NPAR,G,AMIN,U,4)
      CALL MPRINT(1,AMIN)
      CALL FCN(NPAR,G,F,U,4)
      IF (F .NE. AMIN) GO TO 200
      NFCN = 3
      CALL COMAND
      IF (CWORD .NE. 'END RETURN') GOTO 100
      RETURN
  200 WRITE (ISYSWR,1002) AMIN,F
      STOP
 1001 FORMAT(' First entry to FCN.')
 1002 FORMAT(' *** FCN is time-dependent:'E20.12' (1)'E20.12' (2).')
      END


      SUBROUTINE MINTIO(IREAD,IWRITE,IPUNCH)
**
**       User callable routine to allow redefinition of I/O streams must
**       be  called  before call to MINUIT to get no I/O  on the default
**       streams  but  may  be  called at any time to switch the current
**       stream(s) to new ones.
**

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)


      CHARACTER*60 TITLE
      CHARACTER*10 CWORD
      CHARACTER*10 NAM

      COMMON
     . /MINERR/ ERP  (100),ERN  (100)
     . /DERIVA/ G    (100),G2   (100)
     . /LIMITS/ ALIM (100),BLIM (100),LCODE (100),LCORSP(100) ,LIMSET
     . /PAREXT/ U    (100),NAM  (100),WERR  (100),MAXEXT     ,NU
     . /PARINT/ X (50)   ,XT   (50),DIRIN (50),MAXINT     ,NPAR
     . /FIX   / XS(50)   ,XTS  (50),WTS   (50),IPFIX(50)  ,NPFIX
     . /SIMVEC/ P (50,51),PSTAR(50),PSTST (50),PBAR (50)  ,PRHO(50)
     . /VARIAN/ V (50,50)
     . /VARIAT/ VT(50,50)
     . /CASC  / JH,JL,Y(51)

      COMMON
     . /UNIT  / ISYSRD,ISYSWR  ,ISYSPU
     . /TITLE / TITLE ,ISW  (7),NBLOCK
     . /CARD  / CWORD ,WORD7(7)
     . /CONVER/ EPSI  ,APSI,VTEST ,NSTEPQ,NFCN,NFCNMX
     . /MINIMA/ AMIN  ,UP  ,NEWMIN,ITAUR ,SIGMA

      ISYSRD = IREAD
      ISYSWR = IWRITE
      ISYSPU = IPUNCH
      RETURN
      END


      SUBROUTINE MIDATA(MODE)
**
**       Reads the data cards (title card and parameter cards) and  sets
**       up  the starting parameter lists. Control then passes to COMAND
**       for reading the command cards.
**

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)


      CHARACTER*60 TITLE
      CHARACTER*10 CWORD
      CHARACTER*10 NAM

      COMMON
     . /MINERR/ ERP  (100),ERN  (100)
     . /DERIVA/ G    (100),G2   (100)
     . /LIMITS/ ALIM (100),BLIM (100),LCODE (100),LCORSP(100) ,LIMSET
     . /PAREXT/ U    (100),NAM  (100),WERR  (100),MAXEXT     ,NU
     . /PARINT/ X (50)   ,XT   (50),DIRIN (50),MAXINT     ,NPAR
     . /FIX   / XS(50)   ,XTS  (50),WTS   (50),IPFIX(50)  ,NPFIX
     . /SIMVEC/ P (50,51),PSTAR(50),PSTST (50),PBAR (50)  ,PRHO(50)
     . /VARIAN/ V (50,50)
     . /VARIAT/ VT(50,50)
     . /CASC  / JH,JL,Y(51)

      COMMON
     . /UNIT  / ISYSRD,ISYSWR  ,ISYSPU
     . /TITLE / TITLE ,ISW  (7),NBLOCK
     . /CARD  / CWORD ,WORD7(7)
     . /CONVER/ EPSI  ,APSI,VTEST ,NSTEPQ,NFCN,NFCNMX
     . /MINIMA/ AMIN  ,UP  ,NEWMIN,ITAUR ,SIGMA

      CHARACTER*1  LETTER
      INTEGER*4  DATE,TIME
      CHARACTER*10 NAMK,NAMZ
      DATA NAMZ   /'          '/
      DATA MNINIT /0/
C
C --- MODE = 0 (new), = 1 (old)
C
      IF (MODE .EQ. 1) GOTO 105
C
C --- Initialize new data block
C
      IF (MNINIT .EQ. 0) NBLOCK=0
      MNINIT = 1
      NBLOCK = 1 + NBLOCK
      NPFIX  = 0
      NINT   = 0
      NU     = 0
      NPAR   = 0
      IFATAL = 0
      DO 50 I= 1,7
   50 ISW(I) = 0
      ISW(5) = 1

      VERSN = 3.0D0
      SIGMA = 0.0D0
      UP    = 1.0D0

      DO 100 I= 1, MAXEXT

      U  (I) = 0.0D0
      ERP(I) = 0.0D0
      ERN(I) = 0.0D0

      NAM   (I) = NAMZ
      LCODE (I) = 0
  100 LCORSP(I) = 0
C
C --- ISW(6) = 0 for file input, = 1 for keyboard input
C
      IUNIT = ISYSRD
      IF (ISYSRD .NE. 0) THEN
        READ (ISYSRD,1000) TITLE
      ELSE
        ISW(6) = 1
        TITLE  = 'Interactive MINUIT'
        WRITE (0,1015)
        READ  (0,1017) LETTER
        CALL   UPCHAR (LETTER)
        IF (LETTER .EQ. 'Y') THEN
          IUNIT = ISYSPU
          WRITE (0,1016)
          READ  (0,1017) LETTER
          CALL   UPCHAR (LETTER)
          IF (LETTER .EQ. 'Y') REWIND IUNIT
          READ (IUNIT,1000) TITLE
        ENDIF
      ENDIF
      CALL DATIMH(DATE,TIME)
      WRITE (ISYSWR,1004) VERSN,MAXINT,MAXEXT,NBLOCK
      WRITE (ISYSWR,1005)
      WRITE (ISYSWR,1018) TITLE,DATE,TIME
      WRITE (ISYSWR,1005)
C
C --- Read parameter cards
C
  105 ISW(2) = 0
      DO 200 I= 1, 200
  110 IF ((ISW(6) .EQ. 1) .AND. (IUNIT .EQ. ISYSRD)) WRITE (0,1014)
      READ (IUNIT,1001,ERR=110) XK,NAMK,UK,WK,A,B

      K = XK + 0.1D0

      NU = MAX0(NU,K)
      IF (K .LE. 0     ) GOTO 250
      IF (K .LE. MAXEXT) GOTO 115
      IFATAL = IFATAL + 1
      WRITE (ISYSWR,1009) K
      WRITE (ISYSWR,1002) K,NAMK,UK,WK,A,B
      GOTO 200
  115 CONTINUE
      IF(NAM(K).EQ.NAMZ) GOTO 120
C
C --- Previously defined parameter is being redefined
C
      WRITE (ISYSWR,1007)

      IF(WERR(K) .GT. 0.0D0) NINT=NINT-1

  120 CONTINUE
      NAM (K) = NAMK
      U   (K) = UK
      WERR(K) = WK

      IF (WK .GT. 0.0D0) GOTO 122

C
C --- Fixed parameter
C
      WRITE (ISYSWR,1002) K,NAMK,UK
      LCODE(K) = 0
      GOTO 160
C
C --- Variable parameter
C
  122 WRITE (ISYSWR,1002) K,NAMK,UK,WK,A,B
      NINT = NINT + 1
      IF (A) 140,130,140
  130 IF (B) 140,135,140
  135 LCODE(K) = 1
      GOTO 160
  140 IF (B-A) 145,142,150
  142 IFATAL = IFATAL + 1
      WRITE (ISYSWR,1010)
      GOTO 150
  145 SAV = B
      B   = A
      A   = SAV
      WRITE (ISYSWR,1003)
  150 ALIM (K) = A
      BLIM (K) = B
      LCODE(K) = 4
      IF ((B-U(K))*(U(K)-A)) 153,155,160
  153 IFATAL = IFATAL + 1
      WRITE (ISYSWR,1011)
      GOTO 160
  155 WRITE (ISYSWR,1006)
  160 CONTINUE
  200 CONTINUE
      IFATAL = IFATAL + 1
      WRITE (ISYSWR,1012)
C
C --- End parameter cards, STOP if fatal error
C
  250 WRITE (ISYSWR,1005)
      IF (NINT .LE. MAXINT) GOTO 255
      WRITE (ISYSWR,1008)
      IFATAL = IFATAL + 1
  255 IF (IFATAL .LE. 0) GOTO 280
      WRITE (ISYSWR,1013) IFATAL
      STOP
C
C --- Calculate step sizes DIRIN
C
  280 NPAR = 0
      DO 300 K= 1, NU
      IF (LCODE(K) .LE. 0) GOTO 300
      NPAR      = NPAR + 1
      LCORSP(K) = NPAR
      SAV       = U(K)
      X (NPAR)  = PINTF(SAV,K)
      XT(NPAR)  = X(NPAR)
      SAV2      = SAV + WERR(K)
      VPLU      = PINTF(SAV2,K) - X(NPAR)
      SAV2      = SAV - WERR(K)
      VMINU     = PINTF(SAV2,K) - X(NPAR)

      DIRIN(NPAR) = 0.5D0*(DABS(VPLU)+DABS(VMINU))

  300 CONTINUE
      IUNIT = ISYSRD
      RETURN
 1000 FORMAT(A60)
 1001 FORMAT(F10.0,A10,4F10.0)
 1002 FORMAT(I10,2X,A10,2X,2G12.6,2X,2G12.6)
 1003 FORMAT(' *** Above limits have been reversed.')
 1004 FORMAT(''
     . /28X'**********************'
     . /28X'* MINUIT Version'F4.1' *'
     . /28X'* Dimensions'I4'/'I3 ' *'
     . /28X'* Data block No.' i4 ' *')
 1005 FORMAT(1X,78('*'))
 1006 FORMAT(' *** Above parameter is at limit.')
 1007 FORMAT(' *** Previous parameter value ignored.')
 1008 FORMAT(' *** Too many variable parameters.')
 1009 FORMAT(' *** Invalid parameter number'I11)
 1010 FORMAT(' *** Upper and lower limits are equal.')
 1011 FORMAT(' *** Parameter outside limits.')
 1012 FORMAT(' *** More than 200 parameter cards.')
 1013 FORMAT(' *** 'I5' errors on parameter cards. Abort.')
 1014 FORMAT(' |   No   ||  Name  ||  Value |',
     .        '| Error  ||Low lim.|| Up lim.|')
 1015 FORMAT(' Are parameters on file PUNCH [Y/N]'/)
 1016 FORMAT(' Shall  I  rewind  file PUNCH [Y/N]'/)
 1017 FORMAT(A1)
 1018 FORMAT(1X,A60,2(1X,I8))
      END


      SUBROUTINE COMAND
**
**       Reads  the  command  cards and takes appropriate action, either
**       directly by skipping to the corresponding code in COMAND, or by
**       setting up a call to a subroutine.
**

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)


      CHARACTER*60 TITLE
      CHARACTER*10 CWORD
      CHARACTER*10 NAM

      COMMON
     . /MINERR/ ERP  (100),ERN  (100)
     . /DERIVA/ G    (100),G2   (100)
     . /LIMITS/ ALIM (100),BLIM (100),LCODE (100),LCORSP(100) ,LIMSET
     . /PAREXT/ U    (100),NAM  (100),WERR  (100),MAXEXT     ,NU
     . /PARINT/ X (50)   ,XT   (50),DIRIN (50),MAXINT     ,NPAR
     . /FIX   / XS(50)   ,XTS  (50),WTS   (50),IPFIX(50)  ,NPFIX
     . /SIMVEC/ P (50,51),PSTAR(50),PSTST (50),PBAR (50)  ,PRHO(50)
     . /VARIAN/ V (50,50)
     . /VARIAT/ VT(50,50)
     . /CASC  / JH,JL,Y(51)

      COMMON
     . /UNIT  / ISYSRD,ISYSWR  ,ISYSPU
     . /TITLE / TITLE ,ISW  (7),NBLOCK
     . /CARD  / CWORD ,WORD7(7)
     . /CONVER/ EPSI  ,APSI,VTEST ,NSTEPQ,NFCN,NFCNMX
     . /MINIMA/ AMIN  ,UP  ,NEWMIN,ITAUR ,SIGMA

      CHARACTER*10 CNAME(30)
      CHARACTER*4  WD,GOOD,BAD,NONE

      DIMENSION GF( 100)

      DATA GOOD     /'GOOD'/
      DATA BAD      /'BAD '/
      DATA NONE     /'NONE'/
      DATA CNAME( 1)/'MINIMIZE  '/
      DATA CNAME( 2)/'SEEK      '/
      DATA CNAME( 3)/'SIMPLEX   '/
      DATA CNAME( 4)/'MIGRAD    '/
      DATA CNAME( 5)/'MINOS     '/
      DATA CNAME( 6)/'PUNCH     '/
      DATA CNAME( 7)/'PRINTOUT  '/
      DATA CNAME( 8)/'END       '/
      DATA CNAME( 9)/'FIX       '/
      DATA CNAME(10)/'RESTORE   '/
      DATA CNAME(11)/'EXIT      '/
      DATA CNAME(12)/'GRADIENT  '/
      DATA CNAME(13)/'CALL FCN  '/
      DATA CNAME(14)/'MATOUT    '/
      DATA CNAME(15)/'COVARIANCE'/
      DATA CNAME(16)/'IMPROVE   '/
      DATA CNAME(17)/'ERROR DEF '/
      DATA CNAME(18)/'END RETURN'/
      DATA CNAME(19)/'HESSE     '/
      DATA CNAME(20)/'UNIT      '/
      DATA CNAME(21)/'STANDARD  '/
      DATA CNAME(22)/'RELEASE   '/
      DATA CNAME(23)/'HELP      '/
      DATA CNAME(24)/'SET PARAM '/
      DATA CNAME(25)/'PARAMETER '/
      DATA CNAME(26)/'LIMITS    '/
      DATA CNAME(27)/'          '/
      DATA CNAME(28)/'          '/
      DATA CNAME(29)/'          '/
      DATA CNAME(30)/'          '/
      DATA NNAME/30/
      IKARD = 0
   50 IKARD = IKARD + 1
      IF (ISW(6) .EQ. 1) WRITE (0,5017)
      READ (ISYSRD,5000,ERR=85) CWORD,WORD7
      DO 55 I=1, 10
   55 CALL UPCHAR(CWORD(I:I))
      K = 1
      DO 60 I= 1, 7
      IF (WORD7(I)) 65,60,65
   65 K = I + 1
   60 CONTINUE
      IF (CWORD .EQ. CNAME(5)) K = 2
      WRITE (ISYSWR,5003)
      WRITE (ISYSWR,5001) IKARD,CWORD,(WORD7(I),I=1,K-1)
      WRITE (ISYSWR,5003)

      NFCNMX = WORD7(1)+0.5D0


      IF (NFCNMX .LE. 0) NFCNMX = 1000

      EPSI = WORD7(2)

      IF (EPSI .LE. 0.D0) EPSI = 0.1D0*UP

      NEWMIN = 0
      ITAUR  = 0
      ISW(1) = 0
      DO 80 I= 1, NNAME
      IF (CWORD .EQ. CNAME(I)) GOTO 90
   80 CONTINUE
   85 WRITE (ISYSWR,5006)
      GOTO 50
C             1    2    3    4    5    6    7    8    9   10
   90 GOTO ( 100, 200, 300, 400, 500, 600, 700,1100, 900,1000,
     1      1100,1200,1300,1400,1500,1600,1700,1100,1402,2000,
     2      2100,2200,2300,2400,2500,2600,2700,2800,2900,3000),I
C
C --- New minimum
C
   95 WRITE (ISYSWR,5011)
      ITAUR = 0

      EPSI = 0.1D0*UP

C
C --- MINIMIZE
C
  100 NF = NFCN
      CALL SIMPLX
      IF (ISW(1) .GE. 1) GOTO 50
      NFCNMX = NFCNMX+NF-NFCN

  150 VTEST = 0.04D0

      GOTO 460
C
C --- SEEK
C
  200 CALL SEEK
      GOTO 50
C
C --- SIMPLEX
C
  300 CALL SIMPLX
      GOTO 50
C
C --- MIGRAD
C
  400 VTEST = WORD7(3)

      IF (VTEST .LE. 0.0D0) VTEST = 0.01D0

  460 NF   = NFCN
      APSI = EPSI
      CALL MIGRAD
      IF (ISW(2) .GT. 2) GOTO 50
      IF (ISW(1) .EQ. 1) GOTO 50
      NFCNMX = NFCNMX+NF-NFCN
      NF     = NFCN
      CALL SIMPLX
      IF (ISW(1) .EQ. 1) GOTO 50
      NFCNMX = NFCNMX+NF-NFCN
      CALL MIGRAD
      GOTO 50
C
C --- MINOS
C
  500 IF (ISW(2) .LT. 1) GOTO 550

      EPSI  = 0.1D0 * UP
      VTEST = 0.1D0

      CALL MINOS
      IF (NEWMIN .LT. 1) GOTO 50
      GOTO 95
  550 WRITE (ISYSWR,5007)
      GOTO 50
C
C --- PUNCH
C
  600 CALL MPUNCH
      GOTO 50
C
C --- PRINTOUT
C

  700 ISW(5) = WORD7(1)+0.5D0

      GOTO 50
C
C --- FIX
C

  900 IT = WORD7(1)+0.5D0

      CALL FIXPAR(IT,0,ILAX)
      IF (ISW(2).GT.1 .AND. ILAX.NE.0) CALL MPRINT(1,AMIN)
      GOTO 50
C
C --- RESTORE
C

 1000 IT = WORD7(1)+0.5D0

      CALL RESTOR(IT)
      GOTO 50
C
C --- END, EXIT, END RETURN
C

 1100 IT = WORD7(1)+0.5D0

      IF (ISW(4) .EQ. 1 .OR. IT .GT. 0) GOTO 1150
      IFLAG = 3
      WRITE (ISYSWR,5009)
      CALL FCN(NPAR,G,F,U,IFLAG)
      NFCN = NFCN+1
 1150 IF (I .EQ. 11) STOP
      RETURN
C
C --- GRADIENT
C
 1200 ISW(3) = 1

      IF (WORD7(1) .GT. 0.0D0) GOTO 50

      DO 1230 I= 1, NU

 1230 GF(I) = 0.0D0

      CALL INTOEX(X)
      CALL FCN(NPAR,GF,AMIN,U,2)
      NFCN = NFCN + 1
      CALL DERIVE(GF,G2)
      ISW(3) = 0
      CALL DERIVE(G, G2)
      WRITE (ISYSWR,5013)
      ISW(3) = 1
      DO 1250 I= 1, NU
      LC = LCORSP(I)
      IF (LC .EQ. 0) GOTO 1250
      WD = GOOD

      IF (DABS(GF(LC)-G(LC)) .GT. DABS(G2(LC))) WD = BAD
      IF (GF(LC) .EQ. 0.D0) WD = NONE

      IF (WD .NE. GOOD) ISW(3) = 0
      WRITE(ISYSWR,5014) I,NAM(I),GF(LC),G(LC),G2(LC),WD
 1250 CONTINUE
      IF (ISW(3) .EQ. 0) WRITE (ISYSWR,5015)
      GOTO 50
C
C --- CALL FCN
C
 1300 IFLAG = WORD7(1)
      IF (IFLAG .EQ. 3) ISW(4) = 1
      CALL FCN(NPAR,G,F,U,IFLAG)
      NFCN = NFCN + 1
      IF(IFLAG.LE.5) GOTO 50
      CALL EXTOIN(X)
      CALL FCN(NPAR,G,AMIN,U,4)
      NFCN=NFCN+1
      IF (ISW(2) .LE. 1) GOTO 1350
      ISW(2) = 1
      WRITE (ISYSWR,5010)
 1350 CALL MPRINT(1,AMIN)
      GOTO 50
C
C --- MATOUT
C
 1400 IF(ISW(2).GE.2) GOTO 1405
C
C --- HESSE
C
 1402 CALL HESSE
      CALL MPRINT(1,AMIN)

 1405 CALL MATOUT(0.0D0,1)

      GOTO 50
C
C --- COVARIANCE
C

 1500 NRAPE = WORD7(1)+0.5D0

      IF (NRAPE .NE. NPAR) GOTO 1550
      READ (ISYSRD,5002) ((V(I,J),I=1,NRAPE),J=1,NRAPE)
      ISW(2) = 3

      CALL MATOUT(0.0D0,1)

      CALL MPRINT(1,AMIN)
      GOTO 50
 1550 WRITE (ISYSWR,5008)
      NRAP2 = NRAPE**2
      READ (ISYSRD,5002) (G(1),I=1,NRAP2)
      GOTO 50
C
C --- IMPROVE
C
 1600 CONTINUE
      IF (ISW(2) .LT. 2) GOTO 550
      CALL IMPROV
      IF (NEWMIN .EQ. 1) GOTO 150
      GOTO 50
C
C --- ERROR DEF
C
 1700 CONTINUE
      UP = WORD7(1)

      IF (UP .LE. 0.D0) UP = 1.0D0

      IF (ISW(2) .GE. 1) CALL MPRINT(1,AMIN)
      GOTO 50
C
C --- UNIT
C
 2000 CONTINUE

      ISYSRD = WORD7(1)+0.5D0
      ISYSWR = WORD7(2)+0.5D0
      ISYSPU = WORD7(3)+0.5D0

      IF (ISYSRD .EQ. 0) THEN
             ISW(6) = 1
        ELSE
             ISW(6) = 0
        ENDIF
      GOTO 50
C
C --- STANDARD
C
 2100 CALL STAND
      GOTO 50
C
C --- RELEASE
C
 2200 CONTINUE
      DO 2220 IRL=1,7
      KRL = WORD7(IRL)
      IF (KRL .EQ. 0) GOTO 50
      KRL = -IABS(KRL)
      CALL RESTOR(KRL)
 2220 CONTINUE
      GOTO 50
C
C --- HELP
C
 2300 CONTINUE
      WRITE (ISYSWR,5016) (CNAME(I),I=1,30)
      GOTO 50
C
C --- SET PARAM
C
 2400 CONTINUE
      IF (NFCNMX .GT. NU) GOTO 2450
      U(NFCNMX) = WORD7(2)
      CALL EXTOIN(X)
      CALL FCN(NPAR,G,AMIN,U,4)
      NFCN=NFCN + 1
      WRITE (ISYSWR,5010)
      ISW(2) = 1
      CALL MPRINT(1,AMIN)
      GOTO 50
 2450 WRITE (ISYSWR,5012)
      GOTO 50
C
C --- PARAMETER
C
 2500 CONTINUE
      CALL MIDATA(1)
      WRITE (ISYSWR,5010)
      CALL INTOEX(X)
      CALL FCN(NPAR,G,AMIN,U,4)
      CALL MPRINT(1,AMIN)
      GOTO 50
C
C --- LIMITS
C
 2600 CONTINUE
      DO 2650 I2= 1, NU
      IF (LCODE(I2) .EQ. 0       ) GOTO 2650
      IF (NFCNMX    .GT. MAXEXT  ) GOTO 2640
      IF (NFCNMX    .NE. I2      ) GOTO 2650
 2640 IF (WORD7(2)  .NE. WORD7(3)) GOTO 2645
      IF (LCODE(I2) .NE. 1       ) WRITE (ISYSWR,5018) I2
      LCODE(I2) = 1
      GOTO 2650

 2645 ALIM(I2) = DMIN1(WORD7(2),WORD7(3))
      BLIM(I2) = DMAX1(WORD7(2),WORD7(3))

      WRITE (ISYSWR,5019) I2,ALIM(I2),BLIM(I2)
      LCODE(I2) = 4
 2650 CONTINUE
      CALL EXTOIN(X)
      CALL MPRINT(1,AMIN)
      GOTO 50
C
C --- Not yet used
C
 2700 CONTINUE
 2800 CONTINUE
 2900 CONTINUE
 3000 CONTINUE
      GOTO 50
 5000 FORMAT(A10,7F10.0)
 5001 FORMAT(I4,': ',A10,7G9.2)
 5002 FORMAT(7F10.0)
 5003 FORMAT(1X,78('*'))
 5006 FORMAT(' *** Unrecognized command. Ignored.')
 5007 FORMAT(' *** Covariance matrix not exist.')
 5008 FORMAT(' *** Invalid size of covariance matrix.')
 5009 FORMAT(' Call to FCN with IFLAG = 3.')
 5010 FORMAT(/15X'New start point assumed. Covariance matrix lost.')
 5011 FORMAT(/15X'New minimum found! Go back to minimization step!')
 5012 FORMAT(' *** Illegal parameter number requested. Ignored.')
 5013 FORMAT(/
     .20X'Check of gradient calculation in FCN.'/20X,37('-')//
     .11X'Parameter      G(in FCN)   G(MINUIT)     Error   Agreement')
 5014 FORMAT(6X,I5,2X,A10,3E12.4,4X,A4)
 5015 FORMAT(/
     .13X'MINUIT does not accept derivative calculations by FCN.')
 5016 FORMAT(/32X'Command list.'//6(10X,5(A10,2X)/))
 5017 FORMAT(' Command  ||  Arg.1 ||  Arg.2 ||  Arg.3 |'
     .        '|  Arg.4 ||  Arg.5 ||  Arg.6 ||  Arg.7')
 5018 FORMAT(20X'Limits removed from parameter'I4)
 5019 FORMAT(11X'Parameter'I4' limits set to '2G15.5)
      END
      SUBROUTINE SEEK
**
**       Performs a rough minimization by Monte Carlo search.  Each time
**       a  new  minimum  is  found,  the  search area is shifted  to be
**       centered at the best value.  Random points are chosen uniformly
**       over a hypercube determined by current step sizes.
**

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)


      CHARACTER*60 TITLE
      CHARACTER*10 CWORD
      CHARACTER*10 NAM

      COMMON
     . /MINERR/ ERP  (100),ERN  (100)
     . /DERIVA/ G    (100),G2   (100)
     . /LIMITS/ ALIM (100),BLIM (100),LCODE (100),LCORSP(100) ,LIMSET
     . /PAREXT/ U    (100),NAM  (100),WERR  (100),MAXEXT     ,NU
     . /PARINT/ X (50)   ,XT   (50),DIRIN (50),MAXINT     ,NPAR
     . /FIX   / XS(50)   ,XTS  (50),WTS   (50),IPFIX(50)  ,NPFIX
     . /SIMVEC/ P (50,51),PSTAR(50),PSTST (50),PBAR (50)  ,PRHO(50)
     . /VARIAN/ V (50,50)
     . /VARIAT/ VT(50,50)
     . /CASC  / JH,JL,Y(51)

      COMMON
     . /UNIT  / ISYSRD,ISYSWR  ,ISYSPU
     . /TITLE / TITLE ,ISW  (7),NBLOCK
     . /CARD  / CWORD ,WORD7(7)
     . /CONVER/ EPSI  ,APSI,VTEST ,NSTEPQ,NFCN,NFCNMX
     . /MINIMA/ AMIN  ,UP  ,NEWMIN,ITAUR ,SIGMA




      DIMENSION AMID(50),N(50)

      WRITE (ISYSWR,1000)

      NUMBER = WORD7(1) + 0.5D0

      IF (NUMBER .LE. 0) NUMBER = 10 * NPAR
      IFLAG = 4
C
C --- Initial values
C
      DO 100 J= 1, NU
      NI = LCORSP(J)
      IF (NI .LE. 0) GOTO 100
      N   (NI) =   J
      AMID(NI) = U(J)
  100 CONTINUE
      NCYCL = 0
C
C --- Monte Carlo search over entire variable parameter space
C
      DO 650 INUM= 1, NUMBER
      DO 200 I2  = 1, NPAR
      I3 = N(I2)
C
C --- Random points in uniform distribution
C

  180 XPLS = 2.0D0* (RNDM(-1) - 0.5D0)

      U(I3) = AMID(I2) + XPLS*WERR(I3)
      IF (LCODE(I3).LE.1                              ) GOTO 200
      IF (U    (I3).GT.BLIM(I3) .OR. U(I3).LT.ALIM(I3)) GOTO 180
  200 CONTINUE
      CALL FCN(NPAR,G,F,U,IFLAG)
      NFCN = NFCN + 1
      IF (F .GE. AMIN) GOTO 650
      AMIN = F
      DO 500 I= 1, NPAR
      J       = N(I)
  500 AMID(I) = U(J)
      NCYCL   = NCYCL + 1
      IF (ISW(5) .LT. 2) GOTO 650
      CALL EXTOIN(X)
      IF (ISW(5).GE.3 .OR. MOD(NCYCL,10).EQ.1) CALL MPRINT(0,AMIN)
  650 CONTINUE
C
C --- Search finished. Set U to best values.
C
      DO 800 I= 1, NPAR
      NI    = N   (I)
  800 U(NI) = AMID(I)
      CALL EXTOIN(X)
      WRITE (ISYSWR,1005)
      CALL MPRINT(1,AMIN)
      RETURN
 1000 FORMAT(/20X'Start SEEK (Monte Carlo minimum search).')
 1005 FORMAT(/25X'Best value found in SEEK is:')
      END


      SUBROUTINE SIMPLX
**
**       Performs a minimization using the simplex method of Nelder  and
**       Mead (Comp. J., 7, 308, (1965)).
**

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)


      CHARACTER*60 TITLE
      CHARACTER*10 CWORD
      CHARACTER*10 NAM

      COMMON
     . /MINERR/ ERP  (100),ERN  (100)
     . /DERIVA/ G    (100),G2   (100)
     . /LIMITS/ ALIM (100),BLIM (100),LCODE (100),LCORSP(100) ,LIMSET
     . /PAREXT/ U    (100),NAM  (100),WERR  (100),MAXEXT     ,NU
     . /PARINT/ X (50)   ,XT   (50),DIRIN (50),MAXINT     ,NPAR
     . /FIX   / XS(50)   ,XTS  (50),WTS   (50),IPFIX(50)  ,NPFIX
     . /SIMVEC/ P (50,51),PSTAR(50),PSTST (50),PBAR (50)  ,PRHO(50)
     . /VARIAN/ V (50,50)
     . /VARIAT/ VT(50,50)
     . /CASC  / JH,JL,Y(51)

      COMMON
     . /UNIT  / ISYSRD,ISYSWR  ,ISYSPU
     . /TITLE / TITLE ,ISW  (7),NBLOCK
     . /CARD  / CWORD ,WORD7(7)
     . /CONVER/ EPSI  ,APSI,VTEST ,NSTEPQ,NFCN,NFCNMX
     . /MINIMA/ AMIN  ,UP  ,NEWMIN,ITAUR ,SIGMA


      DATA ALPHA,BETA,GAMMA,RHOMIN,RHOMAX/1.D0,.5D0,2.D0,4.D0,8.D0/

      IF (NPAR .LE. 0) RETURN
      NPFN   = NFCN
      NPARP1 = NPAR+1

      RHO1 = 1.0D0 + ALPHA

      RHO2 = RHO1 + ALPHA*GAMMA

      WG = 1.0D0 / NPAR

      IFLAG = 4
      WRITE(ISYSWR,100) EPSI
      DO 2 I= 1, NPAR

      IF (ISW(2) .GE. 1) DIRIN(I) = DSQRT(V(I,I)*UP)
      IF (DABS(DIRIN(I)) .LT. 1.0D-10*DABS(X(I))) DIRIN(I)=1.0D-8*X(I)

      IF(ITAUR.LT. 1) V(I,I) = DIRIN(I)**2/UP
    2 CONTINUE
      IF (ITAUR .LT. 1) ISW(2) = 1
C
C --- Choose the initial simplex using single-parameter searches
C
    1 CONTINUE
      YNPP1     = AMIN
      JL        = NPARP1
      Y(NPARP1) = AMIN
      ABSMIN    = AMIN
      DO 10 I= 1, NPAR
      AMING   = AMIN
      PBAR(I) = X(I)
      BESTX   = X(I)
      KG      = 0
      NS      = 0
      NF      = 0
    4 X(I)    = BESTX + DIRIN(I)
      CALL INTOEX(X)
      CALL FCN(NPAR,G,F,U,4)
      NFCN = NFCN + 1
      IF (F .LE. AMING) GOTO 6
C
C --- Failure
C
      IF (KG .EQ. 1) GOTO 8
      KG = -1
      NF = NF + 1

      DIRIN(I) = DIRIN(I)*(-0.4D0)

      IF (NF .LT. 3) GOTO 4
      NS = 6
C
C --- Success
C
    6 BESTX = X(I)

      DIRIN(I) = DIRIN(I) * 3.0D0

      AMING = F
      KG    = 1
      NS    = NS + 1
      IF (NS .LT. 6) GOTO 4
C
C --- Local minimum found in I-th direction
C
    8 Y(I) = AMING
      IF (AMING .LT. ABSMIN)  JL = I
      IF (AMING .LT. ABSMIN)  ABSMIN = AMING
      X(I) = BESTX
      DO 9 K= 1, NPAR
    9 P(K,I) = X(K)
   10 CONTINUE
      JH   = NPARP1
      AMIN = Y(JL)
      CALL RAZZIA(YNPP1,PBAR)
      DO 20 I= 1, NPAR
   20 X(I) = P(I,JL)
      CALL INTOEX(X)
      IF (ISW(5) .GE. 1) CALL MPRINT(0,AMIN)

      SIGMA = SIGMA * 10.D0

      SIG2  = SIGMA
      IGNAL = 0
      NCYCL = 0
C
C --- Start main loop
C
   50 CONTINUE
      IF (IGNAL.GE.10                      ) GOTO 1
      IF (SIG2 .LT.EPSI .AND. SIGMA.LT.EPSI) GOTO 76
      SIG2 = SIGMA
      IF ((NFCN-NPFN) .GT. NFCNMX) GOTO 78
C
C --- Calculate new point (*) by reflection
C
      DO 60 I= 1, NPAR

      PB = 0.D0

      DO 59 J= 1, NPARP1
   59 PB      = PB + WG*P(I,J )
      PBAR(I) = PB - WG*P(I,JH)

   60 PSTAR(I)=(1.D0+ALPHA)*PBAR(I)-ALPHA*P(I,JH)

      CALL INTOEX(PSTAR)
      CALL FCN(NPAR,G,YSTAR,U,4)
      NFCN=NFCN+1
      IF(YSTAR.GE.AMIN) GOTO 70
C
C --- Point (*) better than JL, calculate new point (**)
C
      DO 61 I=1,NPAR

   61 PSTST(I)=GAMMA*PSTAR(I)+(1.D0-GAMMA)*PBAR(I)

      CALL INTOEX(PSTST)
      CALL FCN(NPAR,G,YSTST,U,4)
      NFCN=NFCN+1
C
C --- Try a parabola through PH, PSTAR, PSTST.  MIN = PRHO.
C
      Y1 = (YSTAR-Y(JH)) * RHO2
      Y2 = (YSTST-Y(JH)) * RHO1

      RHO = 0.5D0* (RHO2*Y1 -RHO1*Y2) / (Y1 -Y2)

      IF (RHO .LT. RHOMIN) GOTO 66
      IF (RHO .GT. RHOMAX) RHO = RHOMAX
      DO 64 I= 1, NPAR

   64 PRHO(I) = RHO*PSTAR(I) + (1.D0-RHO)*P(I,JH)

      CALL INTOEX(PRHO)
      CALL FCN(NPAR,G,YRHO,U,4)
      NFCN = NFCN + 1
      IF (YRHO  .LT. Y(JL) .AND. YRHO .LT. YSTST) GOTO 65
      IF (YSTST .LT. Y(JL)                      ) GOTO 67
      IF (YRHO  .GT. Y(JL)                      ) GOTO 66
C
C --- Accept minimum point of parabola, PRHO
C
   65 CALL RAZZIA (YRHO,PRHO)
      IGNAL = MAX0(IGNAL-2,0)
      GOTO 68
   66 IF (YSTST .LT. Y(JL)) GOTO 67
      IGNAL = MAX0(IGNAL-1,0)
      CALL RAZZIA(YSTAR,PSTAR)
      GOTO 68
   67 IGNAL = MAX0(IGNAL-2,0)
  675 CALL RAZZIA(YSTST,PSTST)
   68 NCYCL=NCYCL+1
      IF (ISW(5) .LT. 2) GOTO 50
      IF (ISW(5) .GE. 3 .OR. MOD(NCYCL, 10) .EQ. 0) CALL MPRINT(0,AMIN)
      GOTO 50
C
C --- Point (*) is not as good as JL
C
   70 IF (YSTAR .GE. Y(JH))  GO TO 73
      JHOLD = JH
      CALL RAZZIA(YSTAR,PSTAR)
      IF (JHOLD .NE. JH)  GO TO 50
C
C --- Calculate new point (**)
C
   73 DO 74 I=1,NPAR

   74 PSTST(I)=BETA*P(I,JH)+(1.D0-BETA)*PBAR(I)

      CALL INTOEX (PSTST)
      CALL FCN(NPAR,G,YSTST,U,4)
      NFCN=NFCN+1
      IF(YSTST.GT.Y(JH)) GO TO 1
C
C --- Point (**) is better than JH
C
      IF (YSTST .LT. AMIN) GOTO 675
      IGNAL = IGNAL + 1
      CALL RAZZIA(YSTST,PSTST)
      GOTO 50
C
C --- End main loop
C
   76 WRITE(ISYSWR,120)
      GOTO 80
   78 WRITE(ISYSWR,130)
      ISW(1) = 1
   80 DO 82 I=1,NPAR

      PB = 0.D0

      DO 81 J=1,NPARP1
   81 PB      = PB + WG*P(I,J)
   82 PBAR(I) = PB - WG*P(I,JH)
      CALL INTOEX(PBAR)
      CALL FCN(NPAR,G,YPBAR,U,IFLAG)
      NFCN = NFCN+1
      IF (YPBAR .LT. AMIN) CALL RAZZIA(YPBAR,PBAR)
      CALL INTOEX(X)
      IF (NFCNMX+NPFN-NFCN .LT. 3*NPAR) GOTO 90

      IF (SIGMA .GT. 2.0D0*EPSI) GOTO 1

   90 CALL MPRINT(1-ITAUR, AMIN)
      RETURN
  100 FORMAT(/26X'Start SIMPLEX minimization.'//4X'Convergence'
     .' criterion: Estimated Distance to Minimum (EDM) <'E10.2)
  120 FORMAT(' SIMPLEX minimization has converged.')
  130 FORMAT(' SIMPLEX terminates without convergence.')
      END


      SUBROUTINE MIGRAD
**
**       Performs  a local  function  minimization  using basically  the
**       method   of  Davidon-Fletcher-Powell  as  modified  by Fletcher
**       (Comp. J., 13, 317 (1970)).
**

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)


      CHARACTER*60 TITLE
      CHARACTER*10 CWORD
      CHARACTER*10 NAM

      COMMON
     . /MINERR/ ERP  (100),ERN  (100)
     . /DERIVA/ G    (100),G2   (100)
     . /LIMITS/ ALIM (100),BLIM (100),LCODE (100),LCORSP(100) ,LIMSET
     . /PAREXT/ U    (100),NAM  (100),WERR  (100),MAXEXT     ,NU
     . /PARINT/ X (50)   ,XT   (50),DIRIN (50),MAXINT     ,NPAR
     . /FIX   / XS(50)   ,XTS  (50),WTS   (50),IPFIX(50)  ,NPFIX
     . /SIMVEC/ P (50,51),PSTAR(50),PSTST (50),PBAR (50)  ,PRHO(50)
     . /VARIAN/ V (50,50)
     . /VARIAT/ VT(50,50)
     . /CASC  / JH,JL,Y(51)

      COMMON
     . /UNIT  / ISYSRD,ISYSWR  ,ISYSPU
     . /TITLE / TITLE ,ISW  (7),NBLOCK
     . /CARD  / CWORD ,WORD7(7)
     . /CONVER/ EPSI  ,APSI,VTEST ,NSTEPQ,NFCN,NFCNMX
     . /MINIMA/ AMIN  ,UP  ,NEWMIN,ITAUR ,SIGMA


      DIMENSION GS( 100),R(50),XXS(50),FLNU(50),VG(50),VII(50)


      DATA SLAMIN,SLAMAX,TLAMIN,TLAMAX/0.2D0,3.0D0,0.05D0,6.0D0/

      IF (NPAR .LE. 0)  RETURN
      ISWTR = ISW(5) - ITAUR
      NPFN  = NFCN
      PARN  = NPAR

      RHO2   = 10.D0  * APSI
      ROSTOP = 1.0D-5 * APSI
      TRACE  = 1.D0

      IFLAG=4
      IF (ISW(3) .EQ. 1) IFLAG = 2
      FS = AMIN
      IF (ITAUR .LT. 1) WRITE (ISYSWR,470) ROSTOP,APSI,VTEST
      GOTO 2
    1 WRITE (ISYSWR,520)
C
C --- Step sizes DIRIN
C
    2 NPARD = NPAR
      DO 3 I= 1, NPAR

      D = 0.02D0 * DABS(DIRIN(I))
      IF (ISW(2).GE.1) D = 0.02D0 * DSQRT(DABS(V(I,I))*UP)
      IF (D .LT. 1.0D-8*DABS(X(I))) D = 1.0D-8 * X(I)

    3 DIRIN(I) = D
C
C --- Starting gradient
C
      NTRY  = 0
    4 NEGG2 = 0
      DO 10 ID= 1, NPARD
      I    = ID + NPAR - NPARD
      D    = DIRIN(I)
      XTF  = X    (I)
      X(I) = XTF + D
      CALL INTOEX(X)
      CALL FCN(NPAR,G,FS1,U,4)
      NFCN = NFCN + 1
      X(I) = XTF  - D
      CALL INTOEX(X)
      CALL FCN(NPAR,G,FS2,U,4)
      NFCN = NFCN + 1
      X(I) = XTF

      GS(I) = (FS1 - FS2             ) / (2.0D0*D)
      G2(I) = (FS1 + FS2 - 2.0D0*AMIN) / D**2
      IF (G2(I) .GT. 1.0D-30) GOTO 10

C
C --- Search if G2 <= 0
C
      WRITE (ISYSWR,520)
      NEGG2 = NEGG2 + 1
      NTRY  = NTRY  + 1
      IF (NTRY .GT. 4) GOTO 230

      D = 50.D0*DABS(DIRIN(I))

      XBEG = XTF

      IF (GS(I) .LT. 0.D0) DIRIN(I) = -DIRIN(I)

      KG = 0
      NF = 0
      NS = 0
    5 X(I) = XTF + D
      CALL INTOEX(X)
      CALL FCN(NPAR,G,F,U,4)
      NFCN = NFCN + 1
      IF (F .LE. AMIN) GOTO 6
C
C --- Failure
C
      IF (KG .EQ. 1) GOTO 8
      KG = -1
      NF = NF + 1

      D = -0.4D0*D

      IF (NF .LT. 10) GOTO 5

      D = 1000.D0*D

      GO TO 7
C
C --- Success
C
    6 XTF = X(I)

      D = 3.0D0*D

      AMIN = F
      KG   = 1
      NS   = NS + 1
      IF (NS   .LT. 10) GOTO 5
      IF (AMIN .LT. FS) GOTO 8

      D = 0.001D0*D

    7 XTF = XBEG

      G2(I) = 1.0D0

      NEGG2 = NEGG2 - 1
    8 X(I)  = XTF

      DIRIN(I) = 0.1D0*D

      FS = AMIN
   10 CONTINUE
      IF (NEGG2 .GE. 1) GOTO 4
      NTRY  = 0
      MATGD = 1
C
C --- Diagonal matrix
C
      IF (ISW(2) .GT. 1) GOTO 15
   11 NTRY  = 1
      MATGD = 0
      DO 13 I= 1, NPAR
      DO 12 J= 1, NPAR

   12 V(I,J) = 0.D0
   13 V(I,I) = 2.0D0/G2(I)

C
C --- Get SIGMA and set up loop
C

   15 SIGMA = 0.D0

      DO 18 I= 1, NPAR

      IF (V(I,I) .LE. 0.D0) GOTO 11
      RI = 0.D0

      DO 17 J= 1, NPAR
      XXS(I) = X(I)
   17 RI     = RI + V(I,J)*GS(J)

   18 SIGMA = SIGMA + GS(I)*RI*0.5D0
      IF (SIGMA .GE. 0.D0) GOTO 20

      WRITE (ISYSWR,520)
      IF (NTRY.EQ.0) GOTO 11
      ISW(2) = 0
      GOTO 230
   20 ISW(2) = 1
      ITER   = 0
      CALL INTOEX(X)
      IF (ISWTR .GE. 1) CALL MPRINT(0,AMIN)

      IF (ISWTR .GE. 2) CALL MATOUT(0.0D0,1)

C
C --- Start main loop
C
   24 CONTINUE

      GDEL = 0.D0

      DO 30 I=1,NPAR

      RI = 0.D0

      DO 25 J=1,NPAR
   25 RI = RI + V(I,J)*GS(J)

      DIRIN(I) = -0.5D0*RI

      GDEL = GDEL + DIRIN(I)*GS(I)
C
C --- Linear search along - VG
C
   30 X(I) =XXS(I) + DIRIN(I)
      CALL INTOEX(X)
      CALL FCN (NPAR,G,F,U,4)
      NFCN=NFCN+1
C
C --- Quadr. interp. using slope GDEL
C

      DENOM = 2.0D0*(F-AMIN-GDEL)
      IF (DENOM .LE. 0.D0) GOTO 35

      SLAM = -GDEL/DENOM
      IF (SLAM .GT. SLAMAX) GOTO 35
      IF (SLAM .LT. SLAMIN) SLAM=SLAMIN
      GOTO 40
   35 SLAM = SLAMAX

   40 IF (DABS(SLAM-1.0D0) .LT. 0.1D0) GOTO 70

      DO 45 I= 1, NPAR
   45 X(I) =XXS(I) + SLAM*DIRIN(I)
      CALL INTOEX(X)
      CALL FCN(NPAR,G,F2,U,4)
      NFCN = NFCN + 1
C
C --- Quadr. interp. using 3 points
C
      AA = FS/SLAM

      BB = F /      (1.0D0-SLAM)
      CC = F2/(SLAM*(SLAM-1.0D0))
      DENOM = 2.0D0*(AA+BB+CC)
      IF (DENOM .LE. 0.D0) GOTO 48
      TLAM = (AA*(SLAM+1.0D0)+BB*SLAM+CC)/DENOM

      IF (TLAM .GT. TLAMAX) GOTO 48
      IF (TLAM .LT. TLAMIN) TLAM=TLAMIN
      GOTO 50
   48 TLAM = TLAMAX
   50 CONTINUE
      DO 51 I= 1, NPAR
   51 X(I) = XXS(I)+TLAM*DIRIN(I)
      CALL INTOEX(X)
      CALL FCN(NPAR,G,F3,U,4)
      NFCN = NFCN + 1
      IF (F .GE.AMIN .AND. F2.GE.AMIN .AND. F3.GE.AMIN) GOTO 200
      IF (F .LT.F2   .AND. F .LT. F3                  ) GOTO 61
      IF (F2.LT.F3                                    ) GOTO 58
   55 F    = F3
      SLAM = TLAM
      GOTO 65
   58 F    = F2
      GOTO 65

   61 SLAM = 1.0D0

   65 DO 67 I= 1, NPAR
      DIRIN(I) = DIRIN(I)*SLAM
   67 X    (I) = DIRIN(I)+XXS(I)
   70 AMIN   = F
      ISW(2) = 2
      IF (SIGMA     +FS-AMIN .LT. ROSTOP) GOTO 170
      IF (SIGMA+RHO2+FS-AMIN .GT. APSI  ) GOTO 75
      IF (TRACE .LT. VTEST) GOTO 170
   75 CONTINUE
      IF (NFCN-NPFN .GE. NFCNMX) GOTO 190
      ITER = ITER + 1
      IF (ISWTR.GE. 3 .OR.(ISWTR.EQ. 2 .AND. MOD(ITER,10) .EQ.1))
     .CALL MPRINT(0,AMIN)
C
C --- Get gradient and SIGMA
C
      IF (ISW(3) .NE. 1) GOTO 80
      CALL INTOEX(X)
      CALL FCN(NPAR,G,AMIN,U,IFLAG)
      NFCN = NFCN + 1
   80 CALL DERIVE(G,G2)
      RHO2 = SIGMA

      SIGMA  = 0.D0
      GVG    = 0.D0
      DELGAM = 0.D0

      DO 100 I= 1, NPAR

      RI  = 0.D0
      VGI = 0.D0

      DO 90 J= 1, NPAR
      VGI = VGI + V(I,J)*(G(J)-GS(J))
   90 RI  = RI  + V(I,J)* G(J)

      R (I) = RI *0.5D0
      VG(I) = VGI*0.5D0

      GAMI   = G(I) - GS(I)
      GVG    = GVG + GAMI*VG(I)
      DELGAM = DELGAM + DIRIN(I)*GAMI
  100 SIGMA  = SIGMA + G(I)*R(I)

      IF (SIGMA  .LT. 0.D0) GOTO 1
      IF (GVG    .LE. 0.D0) GOTO 105
      IF (DELGAM .LE. 0.D0) GOTO 105

      GOTO 107

  105 IF (SIGMA .LT. 0.1D0*ROSTOP) GOTO 170

      GOTO 1
  107 CONTINUE
C
C --- Update covariance matrix
C

      TRACE=0.D0

      DO 120 I= 1, NPAR
      VII(I) = V(I,I)
      DO  120  J=1,NPAR
      D = DIRIN(I)*DIRIN(J)/DELGAM - VG(I)*VG(J)/GVG

  120 V(I,J) = V(I,J) + 2.0D0*D

      IF (DELGAM .LE. GVG) GOTO 135
      DO 125 I= 1, NPAR
  125 FLNU(I) = DIRIN(I)/DELGAM - VG(I)/GVG
      DO 130 I= 1, NPAR
      DO 130 J= 1, NPAR

  130 V(I,J) = V(I,J) + 2.0D0*GVG*FLNU(I)*FLNU(J)

  135 CONTINUE
      DO 140 I= 1, NPAR
      XXS(I) = X(I)
      GS(I)  = G(I)
  140 TRACE  = TRACE + ((V(I,I)-VII(I))/(V(I,I)+VII(I)))**2

      TRACE = DSQRT(TRACE/PARN)

      IF (ISWTR .GE. 4)  CALL MATOUT(TRACE,0)
      FS = F
      GOTO 24
C
C --- End main loop
C
  170 WRITE(ISYSWR,500)
      ISW(2) = 3
      IF(ISWTR .GE. 0) CALL MPRINT(1-ITAUR,AMIN)
      ISWTR = ISWTR - 3*ITAUR
      IF (ISWTR .GT. 0) CALL MATOUT(TRACE,1)
      IF (ITAUR .GT. 0) GOTO 435
      IF (MATGD .GT. 0) GOTO 435
      NPARGD = NPAR*(NPAR+5)/2
      IF (NFCN-NPFN .GE. NPARGD) GOTO 435
      WRITE (ISYSWR,180)
      CALL HESSE
      CALL MPRINT(1,AMIN)

      CALL MATOUT(0.0D0,1)

      IF (ISW(2) .GE. 2) ISW(2) = 3
      GOTO 435
  190 ISW(1) = 1
      GOTO 230
  200 WRITE (ISYSWR,650)
      DO 210 I= 1, NPAR
  210 X(I) = XXS(I)
      CALL INTOEX(X)
      ISW(2) = 1
      IF (SIGMA .LT. ROSTOP) GOTO 170
      IF (MATGD .GT. 0     ) GOTO 2
  230 WRITE (ISYSWR,510)
      CALL INTOEX(X)
      CALL MPRINT(1-ITAUR, AMIN)
      ISWTR = ISW(5) - ITAUR*3
      IF (ISWTR  .LT. 1) GOTO 435
      IF (ISW(2) .LE. 1) GOTO 435
      CALL MATOUT(TRACE,1)
  435 RETURN
  180 FORMAT(' Covariance matrix inaccurate. MINUIT will recalculate.')
  470 FORMAT(/26X'Start MIGRAD minimization.'//4X'Convergence'
     .' criteria: Estimated Distance to Minimum (EDM)  <'E10.2/
     .    22X'or  Estimated Distance to Minimum (EDM)  <'E10.2/
     .    22X'and fractional change in variance matrix <'E10.2)
  500 FORMAT(' MIGRAD minimization has converged.')
  510 FORMAT(' MIGRAD terminated without convergence.')
  520 FORMAT(' Covariance matrix is not positive-definite.')
  650 FORMAT(' MIGRAD fails to find improvement.')
      END


      SUBROUTINE IMPROV
**
**       Attempts to improve on a good local minimum by finding a better
**       ONE.  THE  QUADRATIC PART  OF FCN IS REMOVED BY CALFCN AND THIS
**       transformed function is minimized using the simplex method from
**       several random starting points (Goldstein and Price, Math.Comp.
**       25, 569 (1971))
**

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)


      CHARACTER*60 TITLE
      CHARACTER*10 CWORD
      CHARACTER*10 NAM

      COMMON
     . /MINERR/ ERP  (100),ERN  (100)
     . /DERIVA/ G    (100),G2   (100)
     . /LIMITS/ ALIM (100),BLIM (100),LCODE (100),LCORSP(100) ,LIMSET
     . /PAREXT/ U    (100),NAM  (100),WERR  (100),MAXEXT     ,NU
     . /PARINT/ X (50)   ,XT   (50),DIRIN (50),MAXINT     ,NPAR
     . /FIX   / XS(50)   ,XTS  (50),WTS   (50),IPFIX(50)  ,NPFIX
     . /SIMVEC/ P (50,51),PSTAR(50),PSTST (50),PBAR (50)  ,PRHO(50)
     . /VARIAN/ V (50,50)
     . /VARIAT/ VT(50,50)
     . /CASC  / JH,JL,Y(51)

      COMMON
     . /UNIT  / ISYSRD,ISYSWR  ,ISYSPU
     . /TITLE / TITLE ,ISW  (7),NBLOCK
     . /CARD  / CWORD ,WORD7(7)
     . /CONVER/ EPSI  ,APSI,VTEST ,NSTEPQ,NFCN,NFCNMX
     . /MINIMA/ AMIN  ,UP  ,NEWMIN,ITAUR ,SIGMA


      DATA ALPHA,BETA,GAMMA/1.0D0,0.5D0,2.0D0/

      IF (NPAR .LE. 0) RETURN
      ITAUR = 1

      EPSI = 0.1D0*UP

      NPFN  = NFCN
      NLOOP = WORD7(2)
      IF (NLOOP .LE. 0) NLOOP = NPAR + 4
      NPARP1 = NPAR + 1

      WG = 1.0D0 / NPAR

      IFLAG  = 4
      SIGSAV = SIGMA
      APSI   = AMIN
      DO 2 I = 1, NPAR
      XT(I) = X(I)

      G2(I) = DSQRT(V(I,I)*UP)

      DO 2 J = 1, NPAR
    2 VT(I,J) = V(I,J)
      CALL VERMIN(VT,MAXINT,MAXINT,NPAR,IFAIL)
      IF (IFAIL .GE. 1) GOTO 280
      LOOP = 0
    3 DO 4 I= 1, NPAR

      DIRIN(I) = 2.0D0*G2(I)
    4 X(I) = XT(I) + 2.0D0*DIRIN(I)*(RNDM(-1)-0.5D0)

      LOOP = LOOP + 1

      REG = 2.0D0

      WRITE (ISYSWR,1040) LOOP
    5 AMIN = CALFCN(X)
C
C --- Set up random simplex
C
      JL        = NPARP1
      JH        = NPARP1
      Y(NPARP1) = AMIN
      AMAX      = AMIN
      DO 15 I= 1, NPAR
      XI   = X(I)

      X(I) = XI - DIRIN(I)*(RNDM(-1)-0.5D0)

      Y(I) = CALFCN(X)
      IF (Y(I) .GE. AMIN) GOTO 7
      AMIN = Y(I)
      JL   = I
    7 IF (Y(I) .LE. AMAX) GOTO 8
      AMAX = Y(I)
      JH   = I
    8 CONTINUE
      DO 10 J= 1, NPAR
   10 P(J,I) = X(J)
      P(I,NPARP1) = XI
   15 X(I)        = XI
      SIGMA = AMIN
      SIG2  = SIGMA
C
C --- Start main loop
C
   50 CONTINUE

      IF (AMIN .LT. 0.D0) GOTO 95

      IF (ISW(2) .LT. 2) GOTO 280

      EP = 0.1D0*AMIN

      IF (SIG2.LT.EP .AND. SIGMA.LT.EP) GOTO 100
      SIG2 = SIGMA
      IF ((NFCN-NPFN) .GT. NFCNMX) GOTO 300
C
C --- Calculate new point (*) by reflection
C
      DO 60 I= 1, NPAR

      PB = 0.D0

      DO 59 J= 1, NPARP1
   59 PB      = PB + WG * P(I,J)
      PBAR(I) = PB - WG * P(I,JH)

   60 PSTAR(I)=(1.D0+ALPHA)*PBAR(I)-ALPHA*P(I,JH)

      YSTAR = CALFCN(PSTAR)
      IF(YSTAR.GE.AMIN) GOTO 70
C
C --- Point (*) better than JL, calculate new point (**)
C
      DO 61 I=1,NPAR

   61 PSTST(I)=GAMMA*PSTAR(I)+(1.D0-GAMMA)*PBAR(I)

      YSTST = CALFCN(PSTST)
   66 IF (YSTST .LT. Y(JL)) GOTO 67
      CALL RAZZIA(YSTAR,PSTAR)
      GOTO 50
   67 CALL RAZZIA(YSTST,PSTST)
      GOTO 50
C
C --- Point (*) is not as good as JL
C
   70 IF (YSTAR .GE. Y(JH)) GOTO 73
      JHOLD = JH
      CALL RAZZIA(YSTAR,PSTAR)
      IF (JHOLD .NE. JH) GOTO 50
C
C --- Calculate new point (**)
C
   73 DO 74 I=1,NPAR

   74 PSTST(I)=BETA*P(I,JH)+(1.D0-BETA)*PBAR(I)

      YSTST = CALFCN(PSTST)
      IF(YSTST.GT.Y(JH)) GOTO 5
C
C --- Point (**) is better than JH
C
      IF (YSTST .LT. AMIN) GOTO 67
      CALL RAZZIA(YSTST,PSTST)
      GOTO 50
C
C --- End main loop
C
   95 WRITE (ISYSWR,1000)

      REG = 0.1D0

C
C --- Ask if point is new
C
  100 CALL INTOEX(X)
      CALL FCN(NPAR,G,AMIN,U,4)
      NFCN = NFCN + 1
      DO 120 I= 1, NPAR
      DIRIN(I) = REG*G2(I)

      IF (DABS(X(I)-XT(I)) .GT. DIRIN(I)) GOTO 150

  120 CONTINUE
      GOTO 230
  150 NFCNMX = NFCNMX + NPFN - NFCN
      NPFN   = NFCN
      CALL SIMPLX
      IF (AMIN .GE. APSI) GOTO 325
      DO 220 I= 1, NPAR

      DIRIN(I) = 0.1D0*G2(I)
      IF (DABS(X(I)-XT(I)) .GT. DIRIN(I)) GOTO 250

  220 CONTINUE
  230 IF (AMIN .LT. APSI) GOTO 350
      GOTO 325
C
C --- Truly new minimum
C
  250 NEWMIN = 1
      ISW(2) = 0
      ITAUR  = 0
      NFCNMX = NFCNMX + NPFN - NFCN
      WRITE (ISYSWR,1030)
      RETURN
C
C --- Return to previous region
C
  280 WRITE (ISYSWR,1020)
      ISW(2) = 0
      GOTO 325
  300 ISW(1) = 1
  325 DO 330 I= 1, NPAR

      DIRIN(I) = 0.01D0*G2(I)

  330 X(I)  = XT(I)
      AMIN  = APSI
      SIGMA = SIGSAV
  350 CALL INTOEX(X)
      WRITE (ISYSWR,1010)
      IF (ISW(2) .LT. 2) GOTO 380
      IF (LOOP .LT. NLOOP .AND. ISW(1) .LT. 1) GOTO 3
      ISW(2) = 3
  380 CALL MPRINT(1,AMIN)
      RETURN
 1000 FORMAT(' An improvement on the previous minimum has been found.')
 1010 FORMAT(' IMPROVE has returned to region of original minimum.')
 1020 FORMAT(' Covariance matrix was not positive-definite.')
 1030 FORMAT(/20X'IMPROVE has found a truly new minimum.'
     .       /20X,38('*')/)
 1040 FORMAT(/19X'Start attempt No.'I2' to find new minimum.'
     .       /19X,40('-')/)
      END


      SUBROUTINE MINOS
**
**       Performs  a MINOS  error analysis on those parameters for which
**       it  is  requested  on the MINOS command card. The parameter  in
**       question  is  varied,  and  the  minimum  of  the function with
**       respect  to the  other  parameters is followed until it crosses
**       the value FMIN+UP.
**

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)


      CHARACTER*60 TITLE
      CHARACTER*10 CWORD
      CHARACTER*10 NAM

      COMMON
     . /MINERR/ ERP  (100),ERN  (100)
     . /DERIVA/ G    (100),G2   (100)
     . /LIMITS/ ALIM (100),BLIM (100),LCODE (100),LCORSP(100) ,LIMSET
     . /PAREXT/ U    (100),NAM  (100),WERR  (100),MAXEXT     ,NU
     . /PARINT/ X (50)   ,XT   (50),DIRIN (50),MAXINT     ,NPAR
     . /FIX   / XS(50)   ,XTS  (50),WTS   (50),IPFIX(50)  ,NPFIX
     . /SIMVEC/ P (50,51),PSTAR(50),PSTST (50),PBAR (50)  ,PRHO(50)
     . /VARIAN/ V (50,50)
     . /VARIAT/ VT(50,50)
     . /CASC  / JH,JL,Y(51)

      COMMON
     . /UNIT  / ISYSRD,ISYSWR  ,ISYSPU
     . /TITLE / TITLE ,ISW  (7),NBLOCK
     . /CARD  / CWORD ,WORD7(7)
     . /CONVER/ EPSI  ,APSI,VTEST ,NSTEPQ,NFCN,NFCNMX
     . /MINIMA/ AMIN  ,UP  ,NEWMIN,ITAUR ,SIGMA


      DIMENSION LMI(100),KIRSCH(5),XDEV(50),W( 100)

      DATA APOS /'POSI'/
      DATA ANEG /'NEGA'/
      IF (NPAR .LE. 0) GOTO 700
C
C --- Unpack parameter requests
C
      MARC = 0
      DO 5 I=1,30
    5 LMI(I) = 0
      KNT    = 0
      DO 20 I= 2, 7

      LIME = WORD7(I) + 0.5D0

      IF (LIME .EQ. 0) GOTO 20
      DO 10 K= 1, 5
      K2         = 6 - K
      LEMON      = MOD(LIME,100)
      KIRSCH(K2) = LEMON
   10 LIME       = LIME/100
      DO  15 J= 1, 5
      LEMON = KIRSCH(J)
      IF (LEMON         .GT. MAXEXT) GOTO 14
      IF (LEMON         .EQ. 0     ) GOTO 15
      IF (LCORSP(LEMON) .EQ. 0     ) GOTO 14
      KNT      = KNT + 1
      LMI(KNT) = LEMON
      GOTO  15
   14 MARC = 1
   15 CONTINUE
   20 CONTINUE
      IF (KNT .GT. 0) GOTO 40
      DO  30 I= 1, MAXEXT
      IF (LCORSP(I) .LT. 1) GOTO 30
      KNT      = KNT + 1
      LMI(KNT) = I
      IF (KNT .EQ. 9) GOTO 40
   30 CONTINUE
   40 CONTINUE
      IF (MARC .EQ. 1) WRITE (ISYSWR,811)
      WRITE (ISYSWR,810) (LMI(IZ3),IZ3=1,KNT)
C
C --- Save and prepare start values
C
      SIGSAV = SIGMA
      TOLER  = EPSI

      APSI = EPSI*0.5D0

      ITAUR = 1
      ABEST = AMIN
      AIM   = AMIN + UP
      NSAVE = NFCNMX
      MPAR  = NPAR
      DO 130 I= 1, MPAR
      XT(I) = X(I)
      DO 125 J= 1, MPAR
  125 VT(I,J) = V(I,J)
  130 CONTINUE
      DO 135 I= 1, NU

      ERP(I) = 0.D0
      ERN(I) = 0.D0

  135 W(I) = WERR(I)
      KNT  = 0
C
C --- Start main loop
C
  150 KNT    = KNT + 1
      ISW(1) = 0
      NLIMIT = NFCN + NSAVE
      IF (KNT      .GT. 30) GOTO 590
      IF (LMI(KNT) .LT. 1 ) GOTO 590
      ILAX = LMI(KNT)

      ERP(ILAX) = 0.D0

      IT   = LCORSP(ILAX)
      XTIT = XT(IT)
      CALL INTOEX(XT)
      UT   = U(ILAX)
      NSPT = 0
      IF (LCODE(ILAX) .GT. 1) GOTO 160

      ALIM(ILAX) = UT-100.D0*W(ILAX)
      BLIM(ILAX) = UT+100.D0*W(ILAX)

  160 CONTINUE

      XUNIT = DSQRT(UP/VT(IT,IT))

      MARC = 0
      DO 162 I= 1, MPAR
      IF (I .EQ. IT) GOTO 162
      MARC       = MARC + 1
      XDEV(MARC) = XUNIT*VT(IT,I)
  162 CONTINUE
      CALL FIXPAR(IT,1,ILAX)

      SIG = 1.0D0

      ASIG  = APOS
      DULIM = BLIM(ILAX) - UT
      IF(ISW(2).LT.1) GOTO 460
C
C --- SIG = SIGN of error being calculated
C
  165 WRITE (ISYSWR,806) ASIG,ILAX,NAM(ILAX)
      ITER   = 0
      LIMSET = 0
      DU1    = SIG*W(ILAX)

      IF (DABS(DU1) .LE. DULIM) GOTO 180

      LIMSET = 1
      DU1    = SIG * DULIM

      IF (DULIM .LT. 1.0D-3*W(ILAX)) GOTO 440

  180 U(ILAX) = UT + DU1
      IF (NPAR .EQ. 0) GOTO 205
      FAC = DU1/W(ILAX)
      DO 185 I= 1, NPAR
  185 X(I) = XT(I) + FAC*XDEV(I)
  200 CALL INTOEX (X)
  205 WRITE (ISYSWR,801) ILAX,UT,DU1,U(ILAX)
      CALL FCN(NPAR,G,AMIN,U,4)
      NFCN   = NFCN + 1
      NFCNMX = NLIMIT - NFCN
      CALL MIGRAD
      IF (AMIN   .LT. ABEST) GOTO 650
      IF (ISW(1) .GE. 1    ) GOTO 450
      IF (ISW(2) .GE. 2    ) GOTO 240
      NFCNMX = NLIMIT - NFCN
      CALL SIMPLX
      IF (AMIN   .LT. ABEST) GOTO 650
      IF (ISW(1) .GE. 1    ) GOTO 450
      NFCNMX = NLIMIT - NFCN
      CALL MIGRAD
      IF (AMIN   .LT. ABEST) GOTO 650
      IF (ISW(1) .GE. 1    ) GOTO 450
      IF (ISW(2) .LT. 2    ) GOTO 460
  240 CREM = AMIN - ABEST
      NSPT = NSPT + 1

      IF (CREM .LE. 0.0D0) GOTO 650
      SQUC = DSQRT(UP/CREM)
      IF (DABS(AMIN-AIM) .LT. TOLER) GOTO 400

C
C --- Another iteration necessary
C
      ITER = ITER + 1
      IF (ITER .GT. 6) GOTO 430
      IF (ITER .EQ. 1) GOTO 270
C
C --- Check previous iteration to avoid oscillating
C
      SQUC2 = SQUC

      IF ((SQUC2-1.0D0)*(SQUC1-1.0D0) .GT. 0.0D0) GOTO 270
      SQUC   = 0.65D0*SQUC2+0.35D0
      SQUC11 = 1.0 D0/SQUC1
      IF ((SQUC11-SQUC)*(SQUC -1.0D0) .GT. 0.0D0) GOTO 270

      WRITE (ISYSWR,260)

      SQUC = 0.5D0*SQUC11 + 0.5D0

  270 CONTINUE
      SQUC1 = SQUC
      DU1   = DU1 * SQUC

      IF (DABS(DU1) .LE. DULIM) GOTO 280

      IF (LIMSET .EQ. 1)  GO TO 440
      LIMSET  = 1
      DU1     = SIG*DULIM
  280 U(ILAX) = UT + DU1
      DO 290 I= 1, NPAR
  290 X(I) = XT(I) + SQUC*(X(I)-XT(I))
      GOTO 200
C
C --- Error successfully calculated
C
  400 EROS = DU1 * SQUC
      WRITE (ISYSWR,808) ASIG,ILAX,NAM(ILAX),EROS
  410 WRITE (ISYSWR,812)

      IF (SIG .GT. 0.D0) GOTO 420

      ERN(ILAX) = EROS
      GOTO 500
  420 ERP(ILAX) = EROS

      SIG = -1.0D0

      ASIG  = ANEG
      DULIM = UT - ALIM(ILAX)
      GOTO 165
C
C --- Failure returns
C
  430 WRITE (ISYSWR,809)

      EROS = 0.D0

      GOTO 410
  440 WRITE (ISYSWR,807) ASIG,ILAX,NAM(ILAX),DULIM

      EROS = 0.D0

      GOTO 410
  450 WRITE (ISYSWR,802) NSAVE
      GOTO 500
  460 WRITE (ISYSWR,805)
C
C --- Parameter finished. Reset V.
C
  500 CONTINUE
      CALL RESTOR(1)
      ISW(2) = 3
      DO 560 I= 1, MPAR
      DO 550 J= 1, MPAR
  550 V(I,J) = VT(I,J)
  560 CONTINUE
      IF (NSPT .GE. 6) WRITE(ISYSWR,813) ILAX,NAM(ILAX),NSPT
      GOTO 150
C
C --- Printout final values
C
  590 DO 595 I= 1, MPAR
  595 X(I) = XT(I)
      CALL INTOEX(XT)
      SIGMA = SIGSAV
      AMIN  = ABEST
      CALL MPRINT(2,AMIN)

      CALL MATOUT(0.0D0,1)

      GOTO 700
C
C --- New minimum found
C
  650 NEWMIN = 1
      ISW(2) = 0
      CALL RESTOR(1)
      CALL EXTOIN(X)
      DO 670 I= 1, NPAR

  670 DIRIN(I) = DSQRT(VT(I,I)*UP)

  700 RETURN
  260 FORMAT(' MINOS is having trouble with this parameter ...')
  801 FORMAT(9X'Parameter'I4' set to 'G11.3' + 'G10.3' = 'G12.3)
  802 FORMAT(' *** Call limit. Calculation requires > 'I7' FCN calls.')
  805 FORMAT(' *** MINOS ERROR not calculated for this parameter.')
  806 FORMAT(/
     .1X,78('-')/
     .6X'Determination of 'A4'TIVE MINOS error for parameter'I4,2X,A10/
     .1X,78('-'))
  807 FORMAT(/3X
     .'The 'A4'TIVE MINOS error of parameter'I4,2X,A10,' exceeds 'G12.4)
  808 FORMAT(/5X
     .'The 'A4'TIVE MINOS error of parameter'I4,2X,A10,' is '     G12.4)
  809 FORMAT(' *** Too many iterations.')
  810 FORMAT(' MINOS errors requested for parameters:'9I4,6(/39X,9I4))
  811 FORMAT(' *** Invalid MINOS command card.')
  812 FORMAT(1X,78('*'))
  813 FORMAT(/7X'MINOS finds non-parabolic behavior for parameter'I4,
     .2X,A10,'.'/17X'Error calculation required'I2' minimizations.'/)
      END


      SUBROUTINE HESSE
**
**       Calculates the full second-derivative matrix of FCN  by taking
**       finite differences.    Includes   some   safeguards    against
**       non-positive-definite  matrices,  and  it may set off-diagonal
**       elements to zero in attempt to force positiveness.
**

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)


      CHARACTER*60 TITLE
      CHARACTER*10 CWORD
      CHARACTER*10 NAM

      COMMON
     . /MINERR/ ERP  (100),ERN  (100)
     . /DERIVA/ G    (100),G2   (100)
     . /LIMITS/ ALIM (100),BLIM (100),LCODE (100),LCORSP(100) ,LIMSET
     . /PAREXT/ U    (100),NAM  (100),WERR  (100),MAXEXT     ,NU
     . /PARINT/ X (50)   ,XT   (50),DIRIN (50),MAXINT     ,NPAR
     . /FIX   / XS(50)   ,XTS  (50),WTS   (50),IPFIX(50)  ,NPFIX
     . /SIMVEC/ P (50,51),PSTAR(50),PSTST (50),PBAR (50)  ,PRHO(50)
     . /VARIAN/ V (50,50)
     . /VARIAT/ VT(50,50)
     . /CASC  / JH,JL,Y(51)

      COMMON
     . /UNIT  / ISYSRD,ISYSWR  ,ISYSPU
     . /TITLE / TITLE ,ISW  (7),NBLOCK
     . /CARD  / CWORD ,WORD7(7)
     . /CONVER/ EPSI  ,APSI,VTEST ,NSTEPQ,NFCN,NFCNMX
     . /MINIMA/ AMIN  ,UP  ,NEWMIN,ITAUR ,SIGMA


      DIMENSION YY(50),GY( 100)


      DATA DFWANT,DFZERO,DFMIN,DFMAX/0.01D0,0.00000001D0,0.001D0,0.1D0/

      WRITE (ISYSWR, 500)
      IFLAG = 4
      NPFN  = NFCN
      NPARD = NPAR
C
C --- Diagonal elements
C
      MDIAG = 0
      DO 100 ID= 1, NPARD
      I = ID + NPAR - NPARD

      D = 0.02D0 * DABS(DIRIN(I))
      IF (ISW(2) .GE. 1)  D = 0.02D0 * DSQRT(DABS(V(I,I))*UP)
      DIRMIN = DFZERO*DABS(X(I))

      IF (D .LT. DIRMIN) D = DIRMIN
      DO 20 J= 1, NPAR

   20 V(I,J) = 0.D0

      ICYC     = 0
   40 DIRIN(I) = D
      XTF      = X(I)
      X(I)     = XTF + D
      CALL INTOEX(X)
      CALL FCN(NPAR,GY,FS1,U,IFLAG)
      NFCN = NFCN + 1
      X(I) = XTF  - D
      CALL INTOEX(X)
      CALL FCN(NPAR,GY,FS2,U,IFLAG)
      NFCN = NFCN + 1
      X(I) = XTF
C
C --- Check if step sizes appropriate
C
      ICYC = ICYC + 1
      IF (ICYC .GE. 4) GOTO 55

      DF = DMAX1(DABS(FS1-AMIN),DABS(FS2-AMIN))/UP

      IF (DF .GT. DFMIN ) GOTO 45
      IF (DF .GT. DFZERO) GOTO 50

      D = D*1000.D0

      GOTO 40
   45 IF (DF .LT. DFMAX) GOTO 55

   50 CHAN = DSQRT(DFWANT/DF)
      IF (CHAN .LT. 0.001D0) CHAN = 0.001D0

      D = D*CHAN
      GOTO 40
   55 CONTINUE
C
C --- Get first and second derivative
C

      G (I) = (FS1-FS2               )/(2.0D0*D)
      G2(I) = (FS1 + FS2 - 2.0D0*AMIN)/D**2

      YY(I) = FS1

      IF (DABS(G(I)) + DABS(G2(I)).GT. 1.0D-30) GOTO 80

C
C --- Fix a parameter if  G = G2 = 0.0
C
      IF (ITAUR .GE. 1) GOTO 85
      ISW(2) = 0
      CALL FIXPAR(I,1,IFIX)
      WRITE (ISYSWR,460) IFIX,NAM(IFIX),G(I),G2(I)
      IF (NPAR .EQ. 0) MDIAG = 1
      GOTO 100

   80 IF (G2(I) .GT. 1.0D-30) GOTO 90

   85 MDIAG = 1
      WRITE (ISYSWR,510) I
   90 V(I,I) = G2(I)
  100 CONTINUE
      CALL INTOEX(X)
      IF (MDIAG .EQ. 1) GOTO 390
      ISW(2) = 1
C
C --- Off-diagonal elements
C
      IF (NPAR .EQ. 1) GOTO 214
      NPARM1 = NPAR - 1
      DO 200 I= 1, NPARM1
      IP1 = I + 1
      DO 180 J= IP1, NPAR
      IF (NFCNMX-NFCN+NPFN .LT. NPAR) GOTO 210
      XTI  = X(I)
      XTJ  = X(J)
      X(I) = XTI + DIRIN(I)
      X(J) = XTJ + DIRIN(J)
      CALL INTOEX(X)
      CALL FCN(NPAR,GY,FS1,U,IFLAG)
      NFCN = NFCN + 1
      X(I) = XTI
      X(J) = XTJ
      ELEM = (FS1+AMIN-YY(I)-YY(J)) / (DIRIN(I)*DIRIN(J))
      IF (ELEM**2 .LT. G2(I)*G2(J)) GOTO 170

      ELEM = 0.D0

      WRITE (ISYSWR, 470) I,J
  170 V(I,J) = ELEM
      V(J,I) = ELEM
  180 CONTINUE
  200 CONTINUE
      GOTO 214
  210 J = J - 1
      WRITE (ISYSWR, 490) I,J
  214 CALL INTOEX(X)
      CALL VERMIN(V,MAXINT,MAXINT,NPAR,IFAIL)
      IF (IFAIL .LT. 1) GOTO 222
      WRITE (ISYSWR,520)
C
C --- Diagonal matrix only
C
  216 WRITE (ISYSWR,540)
      ISW(2) = 1
      DO 220 I= 1, NPAR
      DO 218 J= 1, NPAR

  218 V(I,J) = 0.D0
  220 V(I,I) = 1.0D0/G2(I)

      MDIAG = 0
      GOTO 223
  222 WRITE (ISYSWR,480)
      ISW(2) = 2
C
C --- Calculate EDM
C
  223 DO 225 I= 1, NPAR
      DO 225 J= 1, NPAR

  225 V(I,J) = 2.0D0* V(I,J)
      SIGMA  = 0.D0

      DO 250 I= 1, NPAR

      IF (V(I,I) .GT. 0.D0) GOTO 228

      WRITE (ISYSWR,510) I
      MDIAG = 1

  228 R = 0.D0

      DO 240 J= 1, NPAR
      IF (I .EQ. J) GOTO 230

      IF (V(I,J)**2 .LT. DABS(V(I,I)*V(J,J))) GOTO 230

      WRITE (ISYSWR, 470) I,J

      V(I,J) = 0.D0
      V(J,I) = 0.D0

  230 CONTINUE
  240 R = R + V(I,J) * G(J)

  250 SIGMA = SIGMA + 0.5D0*R*G(I)

      IF (MDIAG .EQ. 1) GOTO 390

      IF (SIGMA .GT. 0.D0) GOTO 400

      WRITE (ISYSWR,530)
      GOTO 216
  390 ISW(2) = 0
  400 RETURN
  460 FORMAT(/' Parameter'I3,2X,A10' has been fixed.'7X
     .        'First  derivative is 'E11.3/1X,40('*')7X
     .        'Second derivative is 'E11.3/)
  470 FORMAT(' Covariance matrix not positive-definite.'
     .       ' Faulty element in position'2I4)
  480 FORMAT(21X' Second derivative matrix inverted.')
  490 FORMAT(' Call limit in HESSE. Off-diagonal elements'
     .       ' calculated only up to'2I4)
  500 FORMAT(/20X' Start second derivative calculation.')
  510 FORMAT(' Diagonal element 'I5' is zero or negative.')
  520 FORMAT(' Matrix inversion fails.')
  530 FORMAT(' Matrix not positive-definite.')
  540 FORMAT(' Only diagonal matrix produced.')
      END
      SUBROUTINE FIXPAR(I2,KODE,ILAX)
**
**       Removes  parameter  I2  from the  internal (variable) parameter
**       list,  and  arranges the rest of the list to fill the hole.  If
**       KODE=0, I2 is an external number, otherwise internal.  ILAX  is
**       returned as the external number of the parameter.
**

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)


      CHARACTER*60 TITLE
      CHARACTER*10 CWORD
      CHARACTER*10 NAM

      COMMON
     . /MINERR/ ERP  (100),ERN  (100)
     . /DERIVA/ G    (100),G2   (100)
     . /LIMITS/ ALIM (100),BLIM (100),LCODE (100),LCORSP(100) ,LIMSET
     . /PAREXT/ U    (100),NAM  (100),WERR  (100),MAXEXT     ,NU
     . /PARINT/ X (50)   ,XT   (50),DIRIN (50),MAXINT     ,NPAR
     . /FIX   / XS(50)   ,XTS  (50),WTS   (50),IPFIX(50)  ,NPFIX
     . /SIMVEC/ P (50,51),PSTAR(50),PSTST (50),PBAR (50)  ,PRHO(50)
     . /VARIAN/ V (50,50)
     . /VARIAT/ VT(50,50)
     . /CASC  / JH,JL,Y(51)

      COMMON
     . /UNIT  / ISYSRD,ISYSWR  ,ISYSPU
     . /TITLE / TITLE ,ISW  (7),NBLOCK
     . /CARD  / CWORD ,WORD7(7)
     . /CONVER/ EPSI  ,APSI,VTEST ,NSTEPQ,NFCN,NFCNMX
     . /MINIMA/ AMIN  ,UP  ,NEWMIN,ITAUR ,SIGMA


      DIMENSION V1(2),YY(50)

      EQUIVALENCE (V(1,1),V1(1))

      DATA EPSMAC/1.0D-5/

      IF (KODE) 250,50,150
C
C --- Ext. parameter specified
C
   50 I = I2
      IF (I .GT. NU) GOTO 70
      IF (I .LT.  1) GOTO 70
   60 IF (LCORSP(I)) 70,70,80
C
C --- Error return. Parameter already fixed.
C
   70 ILAX = 0
      WRITE (ISYSWR,500) I
      GOTO 300
   80 LC           = LCORSP(I)
      IT           = LC
      LCORSP(I)    = 0
      ILAX         = I
      NPAR         = NPAR  - 1
      NPFIX        = NPFIX + 1
      IPFIX(NPFIX) = I
      XS   (NPFIX) = X (LC)
      XTS  (NPFIX) = XT(LC)

      EPS = DABS(DIRIN(LC))*10.D0
      IF (ISW(2) .GE. 1) EPS = EPS + DSQRT(DABS(V(LC,LC))*UP)
      IF (EPS .LT. EPSMAC*DABS(X(LC))) EPS=EPSMAC*X(LC)
      WTS(NPFIX) = EPS*0.1D0

      DO 100 IK= I, NU
      IF (LCORSP(IK)) 100,100,85
   85 LC         = LCORSP(IK) - 1
      LCORSP(IK) = LC
      X     (LC) = X    (LC+1)
      XT    (LC) = XT   (LC+1)
      DIRIN (LC) = DIRIN(LC+1)
  100 CONTINUE
      IF (ISW(2) .GT. 1) GOTO 250
      ISW(2) = 0
      GO TO 300
C
C --- Int. parameter specified
C
  150 CONTINUE
      DO 200 IQ= 1, NU
      IF (LCORSP(IQ) .NE. I2) GOTO 200
      I = IQ
      GOTO 60
  200 CONTINUE
      GOTO 70
C
C --- Remove one row and one column from variance matrix
C
  250 KON = 0
      IF (NPAR .LE. 0) GOTO 300
      KON2 = 0
      MPAR = NPAR + 1
      DO 260 I= 1, MPAR
  260 YY(I)=V(I,IT)
      DO 294 I= 1, MPAR
      IF (I.EQ.IT) GOTO 294
      KON2 = KON2 + 1
      DO 292 J= 1, MPAR
      IF (J .EQ. IT) GOTO 292
      KON     = KON + 1
      V1(KON) = V(J,I)-YY(J)*YY(I)/YY(IT)
  292 CONTINUE
      KON = MAXINT*KON2
  294 CONTINUE
C
C --- Check for well-behaved final matrix
C
      DO 295 I= 1, NPAR

      IF (V(I,I) .LE. 0.D0) GOTO 296

      DO 295 J= 1, NPAR
      IF (I .EQ. J) GOTO 295

      IF (V(I,J)**2 .GE. V(I,I)*V(J,J)) V(I,J) = 0.D0

  295 CONTINUE
      GOTO 300
  296 ISW(2) = 0
      WRITE (ISYSWR,501)
  300 RETURN
  500 FORMAT(' *** Error in FIXPAR. Parameter'I4' was not variable.')
  501 FORMAT(4X'Covariance matrix was ill-conditioned and has been '
     .         'destroyed by FIXPAR.')
      END


      SUBROUTINE RESTOR(K)
**
**       Restores  a fixed  parameter to variable status by inserting it
**       into the internal parameter list at the appropriate place.
**

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)


      CHARACTER*60 TITLE
      CHARACTER*10 CWORD
      CHARACTER*10 NAM

      COMMON
     . /MINERR/ ERP  (100),ERN  (100)
     . /DERIVA/ G    (100),G2   (100)
     . /LIMITS/ ALIM (100),BLIM (100),LCODE (100),LCORSP(100) ,LIMSET
     . /PAREXT/ U    (100),NAM  (100),WERR  (100),MAXEXT     ,NU
     . /PARINT/ X (50)   ,XT   (50),DIRIN (50),MAXINT     ,NPAR
     . /FIX   / XS(50)   ,XTS  (50),WTS   (50),IPFIX(50)  ,NPFIX
     . /SIMVEC/ P (50,51),PSTAR(50),PSTST (50),PBAR (50)  ,PRHO(50)
     . /VARIAN/ V (50,50)
     . /VARIAT/ VT(50,50)
     . /CASC  / JH,JL,Y(51)

      COMMON
     . /UNIT  / ISYSRD,ISYSWR  ,ISYSPU
     . /TITLE / TITLE ,ISW  (7),NBLOCK
     . /CARD  / CWORD ,WORD7(7)
     . /CONVER/ EPSI  ,APSI,VTEST ,NSTEPQ,NFCN,NFCNMX
     . /MINIMA/ AMIN  ,UP  ,NEWMIN,ITAUR ,SIGMA

C
C --- K  =  0 means restore all parameters
C --- K  =  1 means restore the last parameter fixed
C --- K  = -I means restore external parameter I (if possible)
C --- IQ = fix-location where internal parameters were stored
C --- IR = external number of parameter being restored
C --- IS = internal number of parameter being restored
C
      IF (K     .GT. 1) WRITE (ISYSWR,510)
      IF (NPFIX .LT. 1) WRITE (ISYSWR,500)
      IF (K.EQ.1.OR. K.EQ.0) GOTO 40
C
C --- Release parameter with specified external number
C
      KA = IABS(K)
      IF (LCORSP(KA) .EQ. 0) GOTO 15
      WRITE (ISYSWR,540)
      RETURN
   15 IF (NPFIX .LT. 1) GOTO 21
      DO 20 IK= 1, NPFIX
      IF (IPFIX(IK) .EQ. KA) GOTO 24
   20 CONTINUE
   21 WRITE (ISYSWR,530)
      RETURN
   24 IF (IK .EQ. NPFIX) GOTO 40
C
C --- Move specified parameter to end of list
C
      IPSAV  = IPFIX(IK)
      XSSAV  = XS   (IK)
      XTSSAV = XTS  (IK)
      WTSSAV = WTS  (IK)
      IKP1   = IK + 1
      DO 30 I= IKP1,NPFIX
      IPFIX(I-1) = IPFIX(I)
      XS   (I-1) = XS   (I)
      XTS  (I-1) = XTS  (I)
   30 WTS  (I-1) = WTS  (I)
      IPFIX(NPFIX) = IPSAV
      XS(NPFIX)    = XSSAV
      XTS(NPFIX)   = XTSSAV
      WTS(NPFIX)   = WTSSAV
C
C --- Set step sizes to errors if covar. matrix will be lost
C
   40 IF (ITAUR .GE. 1 .OR. ISW(2) .LT. 1) GOTO 50
      DO 45 I= 1, NPAR

   45 DIRIN(I) = DSQRT(DABS(V(I,I))*UP)

C
C --- Restore last parameter in list
C
   50 CONTINUE
      IF (NPFIX .LT. 1) GOTO 300
      IR = IPFIX(NPFIX)
      IS = 0
      DO 100 IJ= IR, NU
      IK = NU + IR - IJ
      IF (LCORSP(IK)) 100,100,85
   85 LC = LCORSP(IK) + 1
      IS = LC - 1
      LCORSP(IK) = LC
      X     (LC) = X    (LC-1)
      XT    (LC) = XT   (LC-1)
      DIRIN (LC) = DIRIN(LC-1)
  100 CONTINUE
      NPAR = NPAR + 1
      IF (IS .EQ. 0) IS = NPAR
      LCORSP(IR) = IS
      IQ         = NPFIX
      X    (IS)  = XS(IQ)
      XT   (IS)  = XTS(IQ)
      DIRIN(IS)  = WTS(IQ)
      NPFIX      = NPFIX - 1
      ISW(2)     = 0
      IF (ITAUR .LT. 1) WRITE(ISYSWR,520) IR,NAM(IR)
      IF (K.EQ.0) GOTO 50
  300 RETURN
  500 FORMAT(' *** No more fixed parameters.')
  510 FORMAT(' *** Argument greater than one.')
  520 FORMAT(16X'Parameter'I4,2X,A10' restored to variable.')
  530 FORMAT(' *** Parameter specified has never been variable.')
  540 FORMAT(' *** Parameter specified is already variable.')
      END


      SUBROUTINE DERIVE(GG,GG2)
**
**       Calculates  the first derivatives of FCN (GG), either by finite
**       differences  or  by transforming the user-supplied  derivatives
**       to internal coordinates, according to whether ISW(3) is zero or
**       one. If ISW(3) = 0, an error estimate GG2 is available.
**

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)


      CHARACTER*60 TITLE
      CHARACTER*10 CWORD
      CHARACTER*10 NAM

      COMMON
     . /MINERR/ ERP  (100),ERN  (100)
     . /DERIVA/ G    (100),G2   (100)
     . /LIMITS/ ALIM (100),BLIM (100),LCODE (100),LCORSP(100) ,LIMSET
     . /PAREXT/ U    (100),NAM  (100),WERR  (100),MAXEXT     ,NU
     . /PARINT/ X (50)   ,XT   (50),DIRIN (50),MAXINT     ,NPAR
     . /FIX   / XS(50)   ,XTS  (50),WTS   (50),IPFIX(50)  ,NPFIX
     . /SIMVEC/ P (50,51),PSTAR(50),PSTST (50),PBAR (50)  ,PRHO(50)
     . /VARIAN/ V (50,50)
     . /VARIAT/ VT(50,50)
     . /CASC  / JH,JL,Y(51)

      COMMON
     . /UNIT  / ISYSRD,ISYSWR  ,ISYSPU
     . /TITLE / TITLE ,ISW  (7),NBLOCK
     . /CARD  / CWORD ,WORD7(7)
     . /CONVER/ EPSI  ,APSI,VTEST ,NSTEPQ,NFCN,NFCNMX
     . /MINIMA/ AMIN  ,UP  ,NEWMIN,ITAUR ,SIGMA


      DIMENSION GG(100),GG2(50),GY( 100)


      DATA EPSMAC/1.0D-5/

      IF (ISW(3) .EQ. 1) GOTO 100
      IFLAG = 4
      DO 46 I=1,NPAR

      EPS = 0.1D0 * DABS(DIRIN(I))
      IF (ISW(2) .GE. 1) EPS = EPS + 0.005D0*DSQRT(V(I,I)*UP)
      IF (EPS .LT. EPSMAC*DABS(X(I))) EPS = EPSMAC*X(I)

      XTF  = X(I)
      X(I) = XTF + EPS
      CALL INTOEX(X)
      CALL FCN(NPAR,GY,FS1,U,IFLAG)
      NFCN = NFCN + 1
      X(I) = XTF - EPS
      CALL INTOEX(X)
      CALL FCN(NPAR,GY,FS2,U,IFLAG)
      NFCN = NFCN + 1
C
C --- First derivative and error on first derivative
C

      GG (I) = (FS1-FS2           )/(2.0D0*EPS)
      GG2(I) = (FS1+FS2-2.0D0*AMIN)/(2.0D0*EPS)

      X(I) = XTF
   46 CONTINUE
      CALL INTOEX(X)
      GOTO 200
C
C --- Derivatives calc by FCN
C
  100 DO 150 I= 1, NU
      LC = LCORSP(I)
      IF (LC       .LT. 1) GOTO 150
      IF (LCODE(I) .GT. 1) GOTO 120
      GG(LC) = GG(I)
      GOTO 150

  120 DD = (BLIM(I)-ALIM(I))*0.5D0*DCOS(X(LC))

      GG(LC) = GG(I)*DD
  150 CONTINUE
  200 RETURN
      END


      SUBROUTINE MPUNCH
**
**       Punches  current parameter  values and step sizes onto cards in
**       format  which  can  be reread  by MINUIT  for restarting.   The
**       covariance matrix is also punched if it exists.
**

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)


      CHARACTER*60 TITLE
      CHARACTER*10 CWORD
      CHARACTER*10 NAM

      COMMON
     . /MINERR/ ERP  (100),ERN  (100)
     . /DERIVA/ G    (100),G2   (100)
     . /LIMITS/ ALIM (100),BLIM (100),LCODE (100),LCORSP(100) ,LIMSET
     . /PAREXT/ U    (100),NAM  (100),WERR  (100),MAXEXT     ,NU
     . /PARINT/ X (50)   ,XT   (50),DIRIN (50),MAXINT     ,NPAR
     . /FIX   / XS(50)   ,XTS  (50),WTS   (50),IPFIX(50)  ,NPFIX
     . /SIMVEC/ P (50,51),PSTAR(50),PSTST (50),PBAR (50)  ,PRHO(50)
     . /VARIAN/ V (50,50)
     . /VARIAT/ VT(50,50)
     . /CASC  / JH,JL,Y(51)

      COMMON
     . /UNIT  / ISYSRD,ISYSWR  ,ISYSPU
     . /TITLE / TITLE ,ISW  (7),NBLOCK
     . /CARD  / CWORD ,WORD7(7)
     . /CONVER/ EPSI  ,APSI,VTEST ,NSTEPQ,NFCN,NFCNMX
     . /MINIMA/ AMIN  ,UP  ,NEWMIN,ITAUR ,SIGMA

      DIMENSION VC(7)
      CHARACTER*10 NAMZ
      WRITE (ISYSPU,1000) TITLE
      DO 200 I= 1, NU
      IF (NAM(I) .EQ. '          ') GOTO 200
C
C --- Avoid punching minus zero
C

   20 IF (WERR(I) .EQ. 0.0D0) WERR(I) = 0.0D0

      IF (LCODE(I) .GT. 1)  GO TO 100
C
C --- Parameter without limits
C
      WRITE (ISYSPU,1001) I,NAM(I),U(I),WERR(I)
      GOTO 200
C
C --- Parameter with limits
C
  100 CONTINUE
      WRITE (ISYSPU,1001) I,NAM(I),U(I),WERR(I),ALIM(I),BLIM(I)
  200 CONTINUE
      WRITE (ISYSPU,1002)
      IF (ISW(2) .LT. 3) GOTO 300
      WRITE (ISYSPU,1003) NPAR
      K  = 0
      KC = 0
      DO 250 I= 1, NPAR
      DO 250 J= 1, NPAR
      K     = K + 1
      VC(K) = V(I,J)
      IF (K .NE. 7) GOTO 250
      K  = 0
      KC = KC + 1
      WRITE(ISYSPU,1004) VC,TI,NBLOCK,KC
  250 CONTINUE
      IF (K .EQ. 0) GOTO 300
      KP1 = K + 1
      DO 260 I= KP1, 7

  260 VC(I) = 0.D0

      KC = KC + 1
      WRITE(ISYSPU,1004) VC,TI,NBLOCK,KC
  300 RETURN
 1000 FORMAT(A60)
 1001 FORMAT(I10,A10,4E10.4)
 1002 FORMAT(' ')
 1003 FORMAT('COVARIANCE',I10)
 1004 FORMAT(7E10.3,F6.3,2I2)
      END


      SUBROUTINE MATOUT(TRACE,KODE)
**
**       Prints  the  covariance  matrix V.  Calculates  and  prints the
**       individual correlation coefficients and global correlations.
**

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)


      CHARACTER*60 TITLE
      CHARACTER*10 CWORD
      CHARACTER*10 NAM

      COMMON
     . /MINERR/ ERP  (100),ERN  (100)
     . /DERIVA/ G    (100),G2   (100)
     . /LIMITS/ ALIM (100),BLIM (100),LCODE (100),LCORSP(100) ,LIMSET
     . /PAREXT/ U    (100),NAM  (100),WERR  (100),MAXEXT     ,NU
     . /PARINT/ X (50)   ,XT   (50),DIRIN (50),MAXINT     ,NPAR
     . /FIX   / XS(50)   ,XTS  (50),WTS   (50),IPFIX(50)  ,NPFIX
     . /SIMVEC/ P (50,51),PSTAR(50),PSTST (50),PBAR (50)  ,PRHO(50)
     . /VARIAN/ V (50,50)
     . /VARIAT/ VT(50,50)
     . /CASC  / JH,JL,Y(51)

      COMMON
     . /UNIT  / ISYSRD,ISYSWR  ,ISYSPU
     . /TITLE / TITLE ,ISW  (7),NBLOCK
     . /CARD  / CWORD ,WORD7(7)
     . /CONVER/ EPSI  ,APSI,VTEST ,NSTEPQ,NFCN,NFCNMX
     . /MINIMA/ AMIN  ,UP  ,NEWMIN,ITAUR ,SIGMA


      DIMENSION VLINE(50)

      IF (ISW(2) .LT. 2) RETURN
      WRITE (ISYSWR,600)

      IF (TRACE .NE. 0.0D0) WRITE (ISYSWR,610) TRACE

      IF (NPAR  .EQ. 0) GOTO 250
      ISWTR = ISW(5) - ITAUR
      IF (ISWTR .LT. 2) GOTO 120
C
C --- Internal covariance matrix
C
      DO 100 I= 1, NPAR
  100 WRITE (ISYSWR,620) (V(I,J),J=1,I)
  120 CONTINUE
C
C --- Correlation coeffs.
C
      IF (KODE .LT. 1) GOTO 500
      IF (NPAR .LE. 1) GOTO 500
      WRITE (ISYSWR, 650)
      NPARM = MIN0(NPAR-1, 18)
      WRITE (ISYSWR,690) (ID,ID=1,NPARM)
      DO 200 I= 2, NPAR
      IM = I-1
      DO 170 J= 1, IM

  170 VLINE(J) = V(I,J)/DSQRT(DABS(V(I,I)*V(J,J)))

  200 WRITE (ISYSWR,660) I,(VLINE(IZ),IZ=1,IM)
  250 CONTINUE
C
C --- Global correlation coeffs.
C
      DO 300 I= 1, NPAR
      DO 300 J= 1, NPAR
  300 P(I,J) = V(I,J)
      CALL VERMIN(P,MAXINT,MAXINT,NPAR,IERR)
      IF(IERR .GT. 0) RETURN
      WRITE (ISYSWR,670)
      DO 400 I= 1, NU
      L = LCORSP(I)
      IF (L .EQ. 0) GOTO 400

      GCC = 1.0D0- 1.0D0/(V(L,L)*P(L,L))

      WRITE(ISYSWR,680) I,NAM(I),GCC
  400 CONTINUE
  500 RETURN
  600 FORMAT(1X,78('.')//26X'Internal covariance matrix.'/)
  610 FORMAT(21X'Last fractional change was 'F10.6/)
  620 FORMAT(7E11.3)
  650 FORMAT(27X'Correlation coefficients.'/)
  660 FORMAT(I4,12F6.3,4(/4X,12F6.3))
  670 FORMAT(30X'Global correlation.'//26X'Parameter      Coefficient')
  680 FORMAT(20X,I5,2X,A10,F13.5)
  690 FORMAT(' Int.'I4,11I6,4(/3X,12I6))
      END


      SUBROUTINE MPRINT(IKODE,FVAL)
**
**       Prints  the  values  of the parameters at the time of the call.
**       Also prints other relevant information such as function  value,
**       Estimated  Distance  to  Minimum, parameter errors, step sizes.
**       According to the value of IKODE,the printout  is  long  format,
**       short format, or MINOS format (0,1,2)
**

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)


      CHARACTER*60 TITLE
      CHARACTER*10 CWORD
      CHARACTER*10 NAM

      COMMON
     . /MINERR/ ERP  (100),ERN  (100)
     . /DERIVA/ G    (100),G2   (100)
     . /LIMITS/ ALIM (100),BLIM (100),LCODE (100),LCORSP(100) ,LIMSET
     . /PAREXT/ U    (100),NAM  (100),WERR  (100),MAXEXT     ,NU
     . /PARINT/ X (50)   ,XT   (50),DIRIN (50),MAXINT     ,NPAR
     . /FIX   / XS(50)   ,XTS  (50),WTS   (50),IPFIX(50)  ,NPFIX
     . /SIMVEC/ P (50,51),PSTAR(50),PSTST (50),PBAR (50)  ,PRHO(50)
     . /VARIAN/ V (50,50)
     . /VARIAT/ VT(50,50)
     . /CASC  / JH,JL,Y(51)

      COMMON
     . /UNIT  / ISYSRD,ISYSWR  ,ISYSPU
     . /TITLE / TITLE ,ISW  (7),NBLOCK
     . /CARD  / CWORD ,WORD7(7)
     . /CONVER/ EPSI  ,APSI,VTEST ,NSTEPQ,NFCN,NFCNMX
     . /MINIMA/ AMIN  ,UP  ,NEWMIN,ITAUR ,SIGMA

      INTEGER*4  DATE,TIME
      CHARACTER*10 NAMZ
C
C --- Print headings
C
      CALL DATIMH(DATE,TIME)
      IF (IKODE .NE. 2) THEN
        WRITE (ISYSWR,1000) FVAL,NFCN,TIME,SIGMA
      ELSE
        WRITE (ISYSWR,1006) TITLE
        WRITE (ISYSWR,1001) FVAL,NFCN,TIME,SIGMA
      ENDIF
C
C --- Loop over parameters
C
      DO 200 I= 1,NU
      IF (NAM(I) .EQ. '          ') GOTO 200
      L = LCORSP(I)
      IF (L .EQ. 0) GOTO 100
C
C --- Variable parameter. Calculate external error if V exists.
C
      IF (ISW(2) .LT. 1) GOTO 20

      DX = DSQRT(DABS(V(L,L)*UP))

      IF (LCODE(I) .LE. 1) GOTO 10
      AL = ALIM(I)
      BA = BLIM(I)-AL

      DU1 = AL+0.5D0*(DSIN(X(L)+DX)+1.0D0)*BA-U(I)
      DU2 = AL+0.5D0*(DSIN(X(L)-DX)+1.0D0)*BA-U(I)
      IF (DX .GT. 1.0D0) DU1 = BA
      DX = 0.5D0*(DABS(DU1)+DABS(DU2))

   10 WERR(I) = DX
   20 X1 = X    (L)
      X2 = DIRIN(L)
      IF (IKODE .LT. 2) GOTO 30
      X1 = ERP(I)
      X2 = ERN(I)
   30 WRITE (ISYSWR,1002) L,I,NAM(I),U(I),WERR(I),X1,X2
      IF (LCODE(I) .LE. 1) GOTO 200

      IF (DABS(DCOS(X(L))) .LT. 0.001D0) WRITE (ISYSWR,1004)

      GOTO 200
C
C --- Fixed parameter. Print only if IKODE > 0
C
  100 IF (IKODE .EQ. 0) GOTO 200
      WRITE (ISYSWR,1003) I,NAM(I),U(I)
  200 CONTINUE
      IF (IKODE.GE.1 .AND. ISW(2).GE.1) WRITE (ISYSWR,1005) UP
      RETURN
 1000 FORMAT(/20X'FCN Value    Calls     Time        EDM',
     .       /14X,  E16.7     ,  I7 , 2X, I10   ,   E11.2,/
     .       /' Int.Ext. Parameter        Value        '
     .        ' Error        Int.Value   Int.Step Size')
 1001 FORMAT(/20X'FCN Value    Calls     Time        EDM',
     .       /14X,  E16.7     ,  I7 , 2X, I10   ,   E11.2,/
     .       /39X'Parabolic    ...... MINOS errors ......'
     .       /' Int.Ext. Parameter        Value        '
     .        ' Error        Positive      Negative')
 1002 FORMAT(I4,I4,3X,A10,4E14.5)
 1003 FORMAT(4X,I4,3X,A10, E14.5)
 1004 FORMAT(' *** Above parameter is at limit.')
 1005 FORMAT(/14X'Errors correspond to function change of 'F10.4/)
 1006 FORMAT(//21X'Results of full MINOS error analysis.'
     .       /1X,78('*')/10X,A60/1X,78('*')/)
      END


      SUBROUTINE RAZZIA(YNEW,PNEW)
**
**       Called only by SIMPLX and IMPROV to add a new point  and remove
**       an  old one  from  the current simplex,  and get the  Estimated
**       Distance to Minimum.
**

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)


      CHARACTER*60 TITLE
      CHARACTER*10 CWORD
      CHARACTER*10 NAM

      COMMON
     . /MINERR/ ERP  (100),ERN  (100)
     . /DERIVA/ G    (100),G2   (100)
     . /LIMITS/ ALIM (100),BLIM (100),LCODE (100),LCORSP(100) ,LIMSET
     . /PAREXT/ U    (100),NAM  (100),WERR  (100),MAXEXT     ,NU
     . /PARINT/ X (50)   ,XT   (50),DIRIN (50),MAXINT     ,NPAR
     . /FIX   / XS(50)   ,XTS  (50),WTS   (50),IPFIX(50)  ,NPFIX
     . /SIMVEC/ P (50,51),PSTAR(50),PSTST (50),PBAR (50)  ,PRHO(50)
     . /VARIAN/ V (50,50)
     . /VARIAT/ VT(50,50)
     . /CASC  / JH,JL,Y(51)

      COMMON
     . /UNIT  / ISYSRD,ISYSWR  ,ISYSPU
     . /TITLE / TITLE ,ISW  (7),NBLOCK
     . /CARD  / CWORD ,WORD7(7)
     . /CONVER/ EPSI  ,APSI,VTEST ,NSTEPQ,NFCN,NFCNMX
     . /MINIMA/ AMIN  ,UP  ,NEWMIN,ITAUR ,SIGMA

      DIMENSION PNEW(50)
      DO 10 I=1,NPAR
 10   P(I,JH)=PNEW(I)
      Y(  JH)=YNEW
      IF(YNEW .GE. AMIN) GOTO 18
      DO 15 I=1,NPAR
 15   X(I)=PNEW(I)
      CALL INTOEX(X)
      AMIN = YNEW
      JL   = JH
 18   CONTINUE
      JH     = 1
      NPARP1 = 1+NPAR
20    DO 25 J=2, NPARP1
      IF (Y(J) .GT. Y(JH)) JH = J
25    CONTINUE
      SIGMA = Y(JH) - Y(JL)

      IF (SIGMA .LE. 0.D0) GOTO 45
      US = 1.0D0/SIGMA

      DO 35 I= 1, NPAR
      PBIG = P(I,1)
      PLIT = PBIG
      DO 30 J= 2, NPARP1
      IF (P(I,J) .GT. PBIG) PBIG = P(I,J)
      IF (P(I,J) .LT. PLIT) PLIT = P(I,J)
   30 CONTINUE
      DIRIN(I) = PBIG - PLIT

      IF (ITAUR .LT. 1 ) V(I,I) = 0.5D0*(V(I,I) +US*DIRIN(I)**2)

   35 CONTINUE
   40 RETURN
   45 WRITE (ISYSWR,1000) NPAR
      GOTO 40
 1000 FORMAT(' *** Function not depends on any of the'I4' variable'
     .       ' parameters. Check FCN.')
      END


      SUBROUTINE VERMIN(A,L,M,N,IFAIL)
**
**       Inverts a symmetric matrix.  Matrix is first scaled to have all
**       ones  on  the diagonal  (equivalent to change of units)  but no
**       pivoting is done since matrix is positive-definite.
**

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)


      CHARACTER*60 TITLE
      CHARACTER*10 CWORD
      CHARACTER*10 NAM

      COMMON
     . /MINERR/ ERP  (100),ERN  (100)
     . /DERIVA/ G    (100),G2   (100)
     . /LIMITS/ ALIM (100),BLIM (100),LCODE (100),LCORSP(100) ,LIMSET
     . /PAREXT/ U    (100),NAM  (100),WERR  (100),MAXEXT     ,NU
     . /PARINT/ X (50)   ,XT   (50),DIRIN (50),MAXINT     ,NPAR
     . /FIX   / XS(50)   ,XTS  (50),WTS   (50),IPFIX(50)  ,NPFIX
     . /SIMVEC/ P (50,51),PSTAR(50),PSTST (50),PBAR (50)  ,PRHO(50)
     . /VARIAN/ V (50,50)
     . /VARIAT/ VT(50,50)
     . /CASC  / JH,JL,Y(51)

      COMMON
     . /UNIT  / ISYSRD,ISYSWR  ,ISYSPU
     . /TITLE / TITLE ,ISW  (7),NBLOCK
     . /CARD  / CWORD ,WORD7(7)
     . /CONVER/ EPSI  ,APSI,VTEST ,NSTEPQ,NFCN,NFCNMX
     . /MINIMA/ AMIN  ,UP  ,NEWMIN,ITAUR ,SIGMA


      DIMENSION A(L,M),PP(50),Q(50),S(50)

      IFAIL = 0
      IF (N .LT. 1     ) GOTO 100
      IF (N .GT. MAXINT) GOTO 100
C
C --- Scale matrix by SQRT of diag. elements
C
      DO 8  I=1,N
      SI = A(I,I)
      IF (SI) 100,100,8

    8 S(I) = 1.0D0/DSQRT(SI)

      DO 20 I= 1, N
      DO 20 J= 1, N
   20 A(I,J) = A(I,J)*S(I)*S(J)
C
C --- Start main loop
C
      DO 65 I=1,N
      K = I
C
C --- Preparation for elimination STEP1
C

      Q  (K) = 1.D0/A(K,K)
      PP (K) = 1.0D0
      A(K,K) = 0.0D0

      KP1 = K+1
      KM1 = K-1
      IF (KM1) 100,50,40
   40 DO 49 J=1,KM1
      PP(J) = A(J,K)
      Q (J) = A(J,K)*Q(K)

   49 A(J,K) = 0.D0

   50 IF (K-N) 51,60,100
   51 DO 59 J = KP1,N
      PP(J) = A(K,J)
      Q (J) =-A(K,J)*Q(K)

   59 A(K,J) = 0.0D0

C
C --- Elimination proper
C
   60 DO 65 J=1,N
      DO 65 K=J,N
   65 A(J,K) = A(J,K)+PP(J)*Q(K)
C
C --- Elements of left diagonal and unscaling
C
      DO 70 J= 1, N
      DO 70 K= 1, J
      A(K,J) = A(K,J)*S(K)*S(J)
   70 A(J,K) = A(K,J)
      RETURN
C
C --- Failure return
C
  100 IFAIL = 1
      RETURN
      END


      SUBROUTINE INTOEX(PINT)
**
**       Transforms   from   internal  coordinates  (PINT)  to  external
**       parameters (U).  The minimizing routines which work in internal
**       coordinates call this routine before calling FCN.
**

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)


      CHARACTER*60 TITLE
      CHARACTER*10 CWORD
      CHARACTER*10 NAM

      COMMON
     . /MINERR/ ERP  (100),ERN  (100)
     . /DERIVA/ G    (100),G2   (100)
     . /LIMITS/ ALIM (100),BLIM (100),LCODE (100),LCORSP(100) ,LIMSET
     . /PAREXT/ U    (100),NAM  (100),WERR  (100),MAXEXT     ,NU
     . /PARINT/ X (50)   ,XT   (50),DIRIN (50),MAXINT     ,NPAR
     . /FIX   / XS(50)   ,XTS  (50),WTS   (50),IPFIX(50)  ,NPFIX
     . /SIMVEC/ P (50,51),PSTAR(50),PSTST (50),PBAR (50)  ,PRHO(50)
     . /VARIAN/ V (50,50)
     . /VARIAT/ VT(50,50)
     . /CASC  / JH,JL,Y(51)

      COMMON
     . /UNIT  / ISYSRD,ISYSWR  ,ISYSPU
     . /TITLE / TITLE ,ISW  (7),NBLOCK
     . /CARD  / CWORD ,WORD7(7)
     . /CONVER/ EPSI  ,APSI,VTEST ,NSTEPQ,NFCN,NFCNMX
     . /MINIMA/ AMIN  ,UP  ,NEWMIN,ITAUR ,SIGMA

      DIMENSION PINT(2)
      DO 100 I= 1, NU
      J = LCORSP(I)
      IF (J) 100,100,50
   50 CONTINUE
      IF (LCODE(I) .EQ. 1) GOTO 80
      AL = ALIM(I)

      U(I) = AL+0.5D0*(DSIN(PINT(J))+1.0D0)*(BLIM(I)-AL)

      GOTO 100
   80 U(I) = PINT(J)
  100 CONTINUE
      RETURN
      END


      SUBROUTINE EXTOIN(PINT)
**
**       Transforms the external parameter values  X  to internal values
**       in the dense array PINT. Function PINTF is used.
**

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)


      CHARACTER*60 TITLE
      CHARACTER*10 CWORD
      CHARACTER*10 NAM

      COMMON
     . /MINERR/ ERP  (100),ERN  (100)
     . /DERIVA/ G    (100),G2   (100)
     . /LIMITS/ ALIM (100),BLIM (100),LCODE (100),LCORSP(100) ,LIMSET
     . /PAREXT/ U    (100),NAM  (100),WERR  (100),MAXEXT     ,NU
     . /PARINT/ X (50)   ,XT   (50),DIRIN (50),MAXINT     ,NPAR
     . /FIX   / XS(50)   ,XTS  (50),WTS   (50),IPFIX(50)  ,NPFIX
     . /SIMVEC/ P (50,51),PSTAR(50),PSTST (50),PBAR (50)  ,PRHO(50)
     . /VARIAN/ V (50,50)
     . /VARIAT/ VT(50,50)
     . /CASC  / JH,JL,Y(51)

      COMMON
     . /UNIT  / ISYSRD,ISYSWR  ,ISYSPU
     . /TITLE / TITLE ,ISW  (7),NBLOCK
     . /CARD  / CWORD ,WORD7(7)
     . /CONVER/ EPSI  ,APSI,VTEST ,NSTEPQ,NFCN,NFCNMX
     . /MINIMA/ AMIN  ,UP  ,NEWMIN,ITAUR ,SIGMA

      DIMENSION PINT(2)
      LIMSET=0
      DO 100  I= 1, NU
      J = LCORSP(I)
      IF  (J) 100,100,50
   50 PINT(J) = PINTF(U(I),I)
  100 CONTINUE
      RETURN
      END


      SUBROUTINE STAND
**
**       Optional  user-supplied  subroutine  is  called  whenever   the
**       command "STANDARD" appears.
**
      RETURN
      END



      DOUBLE PRECISION FUNCTION PINTF(PEXTI,I)

**
**       Calculates  the internal parameter value PINTF corresponding to
**       the external value PEXTI for parameter I.
**

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)


      CHARACTER*60 TITLE
      CHARACTER*10 CWORD
      CHARACTER*10 NAM

      COMMON
     . /MINERR/ ERP  (100),ERN  (100)
     . /DERIVA/ G    (100),G2   (100)
     . /LIMITS/ ALIM (100),BLIM (100),LCODE (100),LCORSP(100) ,LIMSET
     . /PAREXT/ U    (100),NAM  (100),WERR  (100),MAXEXT     ,NU
     . /PARINT/ X (50)   ,XT   (50),DIRIN (50),MAXINT     ,NPAR
     . /FIX   / XS(50)   ,XTS  (50),WTS   (50),IPFIX(50)  ,NPFIX
     . /SIMVEC/ P (50,51),PSTAR(50),PSTST (50),PBAR (50)  ,PRHO(50)
     . /VARIAN/ V (50,50)
     . /VARIAT/ VT(50,50)
     . /CASC  / JH,JL,Y(51)

      COMMON
     . /UNIT  / ISYSRD,ISYSWR  ,ISYSPU
     . /TITLE / TITLE ,ISW  (7),NBLOCK
     . /CARD  / CWORD ,WORD7(7)
     . /CONVER/ EPSI  ,APSI,VTEST ,NSTEPQ,NFCN,NFCNMX
     . /MINIMA/ AMIN  ,UP  ,NEWMIN,ITAUR ,SIGMA


      DATA BIG,SMALL /1.570796326795D0,-1.570796326795D0/

      IGO = LCODE(I)
      GOTO (100,200,300,400),IGO
C
C --- IGO = 1 means no limits
C
  100 PINTF = PEXTI
      GOTO 800
  200 CONTINUE
  300 CONTINUE
C
C --- IGO = 4 means there are two limits
C
  400 ALIMI = ALIM(I)
      BLIMI = BLIM(I)
      IF (PEXTI-ALIMI) 440,500,460
  440 A     = SMALL
  450 PINTF = A

      PEXTI = ALIMI+0.5D0*(BLIMI-ALIMI)*(DSIN(A)+1.0D0)

      LIMSET=1
      WRITE (ISYSWR,1000) I
      GOTO 800
  460 IF (BLIMI-PEXTI) 470,520,480
  470 A = BIG
      GOTO 450

  480 YY    = 2.0D0*(PEXTI-ALIMI)/(BLIMI-ALIMI)-1.0D0
      PINTF = DATAN(YY/DSQRT(1.0D0-YY**2))

      GOTO 800
  500 PINTF = SMALL
      GOTO 800
  520 PINTF = BIG
  800 RETURN
 1000 FORMAT(' *** Variable'I4' corrected by PINTF.')
      END



      DOUBLE PRECISION FUNCTION CALFCN(PVEC)

**
**       Called  only  from  IMPROV.   Transforms  the function  FCN  by
**       dividing  out  the  quadratic part  in  order  to find  further
**       minima.  Calculates (F-FMIN)/(X-XMIN)*V*(X-XMIN).
**

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)


      CHARACTER*60 TITLE
      CHARACTER*10 CWORD
      CHARACTER*10 NAM

      COMMON
     . /MINERR/ ERP  (100),ERN  (100)
     . /DERIVA/ G    (100),G2   (100)
     . /LIMITS/ ALIM (100),BLIM (100),LCODE (100),LCORSP(100) ,LIMSET
     . /PAREXT/ U    (100),NAM  (100),WERR  (100),MAXEXT     ,NU
     . /PARINT/ X (50)   ,XT   (50),DIRIN (50),MAXINT     ,NPAR
     . /FIX   / XS(50)   ,XTS  (50),WTS   (50),IPFIX(50)  ,NPFIX
     . /SIMVEC/ P (50,51),PSTAR(50),PSTST (50),PBAR (50)  ,PRHO(50)
     . /VARIAN/ V (50,50)
     . /VARIAT/ VT(50,50)
     . /CASC  / JH,JL,Y(51)

      COMMON
     . /UNIT  / ISYSRD,ISYSWR  ,ISYSPU
     . /TITLE / TITLE ,ISW  (7),NBLOCK
     . /CARD  / CWORD ,WORD7(7)
     . /CONVER/ EPSI  ,APSI,VTEST ,NSTEPQ,NFCN,NFCNMX
     . /MINIMA/ AMIN  ,UP  ,NEWMIN,ITAUR ,SIGMA

      DIMENSION PVEC(50)
      CALL INTOEX(PVEC)
      CALL FCN(NPAR,G,F,U,4)
      NFCN = NFCN + 1
      DO 200 I= 1, NPAR

      G(I) = 0.D0

      DO 200 J= 1, NPAR
  200 G(I) = G(I) + VT(I,J) * (XT(J)-PVEC(J))

      DENOM = 0.D0

      DO 210 I= 1, NPAR
  210 DENOM = DENOM + G(I) * (XT(I)-PVEC(I))

      IF (DENOM .LE. 0.D0) ISW(2) = 0
      IF (DENOM .LE. 0.D0) DENOM  = 1.0D0

      CALFCN = (F-APSI) / DENOM
      RETURN
      END
      SUBROUTINE UPCHAR(LETTER)
**
**       Translate character to upper case.
**
      CHARACTER*1 LETTER
      I = ICHAR  (LETTER)
      IF (I .LE.  96) RETURN
      IF (I .GE. 128) RETURN
      LETTER = CHAR(I-32)
      RETURN
      END

      DOUBLE PRECISION FUNCTION RNDM(DUMMY)
      INTEGER DUMMY
      RNDM=0.D0
      WRITE(*,*) 'WARNING!! ATTEMPT TO USE MONTE-CARLO METHOD!!!!'
      RETURN
      END

      SUBROUTINE DATIMH (DATE,TIME)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      INTEGER*4 DATE,TIME
      INTEGER*2 HR,MIN,SEC,YR,MON,DAY,DUMMY

      CALL GETTIM(HR,MIN,SEC,DUMMY)
      CALL GETDAT(YR,MON,DAY)

      DATE = DAY*10000+MON*100+YR
      TIME = HR*10000+MIN*100+SEC

      RETURN
      END
