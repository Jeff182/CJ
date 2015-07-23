      PROGRAM ALTPFIT
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*10 IPNAME(60),ITITLE(30),IQUIT
     2,LABEL,LABEL2
      Character*30 output
      character*80 lcomment
      logical set
      DIMENSION DATA(6,4500),V(4),NFLAG(30),TX(1500),TQ2(1500)
      DIMENSION XQ(1500),YDAT(1500),ERD(1500),Q2TH(1500),YTH(1500)
     2,YOUT(1500),CHI(1500),THOUT(1500),ITOUT(1500),itype(30),errn(30)
     3,ITGT(30),ICORR(30),ITMC(30) 
      DIMENSION Q2CUT(30),XPOS(350),INORM(30) 
     2,RSPOS(15),XLO(30),XHI(30),W2MIN(30) 
     3,XMLO(15) 
      COMMON/GFUNC/CALC(8,40,50) 
      COMMON/MINIMB/CURPAR(60),UNCRT(60),ERMIN(40),PWATE(60),IFREP(40), 
     * NPARAM,NVAR,IREJ 
      COMMON/MINIMC/IPNAME,PLB(60),PUB(60),TOLRNC,IPAR(60),NSTEPS,
     * IPRNT,IERR,IBND,CHISQR 
      COMMON/Q2STUF/ Q02,Q2MAX
      COMMON/INPUT/INS,NVL,NMAX,DELTA,S0,IORD 
      COMMON/GRPTHY/FLAVOR
      common/cuts/q2cut,xlo,xhi,w2min
      common/expts/ititle,errn,iserr
      COMMON/CONSTANTS/PI,PI2
      COMMON/DYTEST/XSECLO,GLUCOR
      COMMON/FLAGS/ITYPE,NFLAG,INORM,ITGT,ICORR,ITMC
      common/qcdpar/al,nf,norder,set
      common/normalization/nfit
      common/cdf/covinv(33,33),cdfchisqr,ncdf,ncdfl,ncov
      common/outputfile/output
      DATA IQUIT/'END'/
C     
C  FITTING PROGRAM FOR DEEP INELASTIC AND DRELL YAN DATA
C  BASED AROUND THE ALTARELLI-PARISI INTEGRATION SUBROUTINE 
C  MODIFIED TO INCLUDE DRELL YAN DATA 1/13/83.
C 
C  VAX VERSION MODIFIED APRIL 1985 FROM ORIGINAL CYBER SOURCE
C
C  Converted to Linux, Summer 2002. Numerous modifications and upgrades 
c  incorporated.
c
c  Modified to read CTEQ data files
c  Modified to use CTEQ6 parametrization form
c
C  **************************************************************
C  MAXIMUM OF 60 PARAMETERS ALLOWED. OF THESE, A MAXIMUM OF 
C  40 CAN BE VARIED (THIS CAN BE EASILY INCREASED). THE FIRST 
C  35 PARAMETERS ARE RESERVED FOR THE PARTON DISTRIBUTIONS. THE 
C  OTHER 25 ARE FOR NORMALIZATIONS AND OTHER PARAMETERS WHICH DO 
C  NOT AFFECT THE PARTON DISTRIBUTIONS. VARIOUS MODS TO INCREASE 
C  THE EFFICIENCY RELY ON THIS SPLIT.
C  **************************************************************
C
C
C  *****************************************************************
C
      OPEN(UNIT=7,FILE='input.dat',STATUS='unknown')
      read(7,*) output
      open(unit=6,file=output,status='unknown')
C
C 
C  SET UP VARIOUS CONSTANTS AND INITIALIZE INTEGRATION SUBROUTINES
C 
      PI=4.*DATAN(1.D0)
      PI2=PI**2
      NV=6
      CALL WATE4
      CALL WATE8
      CALL WATE16
      CALL WATE32 
      CALL WATE8L 
c      DELTA=0.10
      DELTA=.04
      IERR=0
      IDRV=2
      NSTEPS=15 
      TOLRNC=0.5
      IPRNT=1 
      IPLOT=0
c
c initialize parameters for the jet routines
c
      call setqcd

C 
C  READ IN FIT PARAMETERS 
C
      I=0 
   32 I=I+1 
      READ(7,*) IPNAME(I),CURPAR(I),PWATE(I)
      IPAR(I)=I 
      IF(IPNAME(I).EQ.IQUIT) GO TO 31 
      GO TO 32
   31 NPARAM=I-1
C 
C  READ IN TYPES OF DATA TO BE FIT AND THE CUTS 
C 
      I=0 
   33 I=I+1 
      READ(7,*) ITYPE(I),NFLAG(I),ITITLE(I),INORM(I),ITGT(I),
     $Q2CUT(I),XLO(I),XHI(I),W2MIN(I),icorr(i),itmc(i) 
C
C  TERMINATE LIST WITH ITYPE = 999
C
      IF(ITYPE(I).EQ.999) GO TO 34 
      if(itype(i).le.100)then
         WRITE(6,38) ITITLE(I),Q2CUT(I),XLO(I),XHI(I),W2MIN(I) 
 38      FORMAT(1H ,A10,' DIS data with Q2 greater than',F6.2,
     2   ' and x between ',F5.3,' and ',F5.3,' and W2 >= ',F6.2) 
      If(nflag(i).eq.2.or.nflag(i).eq.9)then
         if(itgt(i).eq.0)then
            write(6,47)
         else if(itgt(i).eq.1)then
            write(6,48)
         endif
      else if(nflag(i).eq.5.or.nflag(i).eq.8.or.nflag(i).eq.10
     2.or.nflag(i).eq.11)then
         if(itgt(i).eq.0)then
            write(6,47)
         else if(itgt(i).eq.1)then
            write(6,35)
         else if(itgt(i).eq.2)then
            write(6,49)
         else if(itgt(i).eq.3)then
            write(6,52)              
         endif 
      endif
 47   format(1h ,'No nuclear correction applied')
 48   format(1h ,'Gomez et al deuteron correction applied')
 49   format(1h ,'EmcNmc and Gomez nuclear corrections applied')
 52   format(1h ,'Kulagin-Petti nuclear correction applied')
 35   format(1h ,'EmcNmc nuclear corrections applied')
         if(itmc(i).eq.1)write(6,53)
 53   format(1h ,'Target Mass Corrections applied')
      else if(itype(i).le.200)then
         write(6,39) ititle(i),q2cut(i)
 39      format(1h ,a10,'vector boson data with Q2 greater than', 
     2   f7.2)
      else if(itype(i).le.300)then
         if(icorr(i).eq.0)then
            write(6,37) ititle(i)
 37         format(1h ,a10,'Single jet inclusive data')
         else
            write(6,36) ititle(i)
 36         format(1h ,a10,
     2'Single jet inclusive data with correlated errors')
         endif
      endif
      goto 33
 34   NFIT=I-1
C 
C  READ IN FLAGS FOR THE FITS 
C  INS=1....NON-SINGLET ONLY
C     =0....INCLUDES SOME SINGLET TERMS 
C  NVL= NUMBER OF NON-SINGLET ARRAYS TO BE EVOLVED (6 in this version)
C  ARRAYS ARE UV,DV, u+, d+, s+, c+,SINGLET,GLUON 
C 
C  IORD=1....NEXT TO LEADING ORDER TERMS INCLUDED 
C      =0....LEADING ORDER TERMS ONLY 
C  ISERR=0...NO SYSTEMATIC ERRORS INCLUDED
C       =1...SYSTEMATIC ERRORS ADDED IN QUADRATURE, WHERE AVAILABLE
C
      READ(7,*) INS,NVL,IORD,Q02,ISERR
      norder=iord+1
C  
C  THE UNIVERSAL PARTON DISTRIBUTION IS USED HERE
C
c      IORDP=IORD
C  IORDP PASSED THROUGH PATRICK'S COMMON BLOCK
      IF(ISERR.NE.0) WRITE(6,1004)
 1004 FORMAT(/' SYSTEMATIC ERRORS ADDED IN QUADRATURE')
c
      Read(7,51) lcomment
      write(6,51) lcomment
 51   format(A80)

c  Read in parameters for the PDF error band gneration
c
      Read(7,*) QSQ,Iband
c
c  Iband = 1  Tables produced with PDFs and errors
c  Iband = 2  Results are normalized to the respective PDFs (bands 
c             given relative to 1.0)
c
C 
C  NXPOS IS THE NUMBER OF DIFFERENT X VALUES APPEARING IN THE 
C  DIS DATA SET. THE VALUES ARE STORED IN THE ARRAY XPOS. 
C  THIS LIST DETERMINES WHICH X VALUES WILL BE PRINTED IN THE 
C  OUTPUT TABLES. 
c
c  Currently set up for data sets without correlated errors
c
C 
C  ITYPE= 3  NFLAG = 1  is  bcdms hydrogen 
C         4        = 2  is  bcdms deuterium
c         5        = 1  is  slac hydrogen
c         6        = 2  is  slac deuterium
c         9        = 1  is  h1 f2 (94 data) h1f2
c        10        = 1  is  zeus f2 zeusf2
c        11        = 1  is  zeus f2 (96/97 data) zeus01c
c        12        = 1  is  h1 f2 (94-97 data) H1F2.99c
c        13        = 1  is  h1 f2 (86-97 e+ data - low x) H1_181c
c        14        = 1  is  h1 f2 (98-99 data) H1_187c
c        51        = 1  is  Nmc F2 p
c        52        = 2  is  Nmc F2 d
c        53        = 9  is  NmcF2r  F2d/F2p (x and Q2)
c        54        = 9  is  NmcF2rX F2d/F2p (x and <Q2>)
c        65        = 5  is  ccfr f2
c        66        = 8  is  ccfr xf3
c        70        =10  is  nutev neutrino iron cross section
c        71        =11  is  nutev antineutrino iron cross section
c
c       106        = 1  is  e605 dimuon on cu target
c       108        = 1  is  e866 dimuon xf pp 
c       110        = 1  is  e866 dimuon xf pd
c       ---        = 2  is  e866 dimuon pp or pn
c       125        = 3  is  NA51 (pp-pn)/(pp+pn) dimuon ratio
c       133        = 4  is  e866f pd/2pp dimuon asymmetry 
c       126        = 5  is  CDF W lepton asymmetry (1996)
c       127        = 5  is  CDF W lepton asymmetry (1998)
c       128        = 5  is  CDF W lepton asymmetry (2005)
c       129        = 5  is  D0  W lepton asymmetry (2008)
c
c       201        = 1  is  CDF inclusive single jet
c       202        = 1  is  D0 inclusive single jet 
c
c       301        = 1  is D0 gamma + jet, rapidity interval 1
c       302        = 2  is D0 gamma + jet, rapidity interval 2
c       303        = 3  is D0 gamma + jet, rapidity interval 3
c       304        = 4  is D0 gamma + jet, rapidity interval 4
C     
      call fildat(nfit,npts,data,q2max)
      call rdxpos(xpos,nxpos,rspos,nrspos,npts,data,ndis,ndy,njet,
     2ngamjet)
c
      IF(IORD) 42,42,43 
   42 WRITE(6,44)
   44 FORMAT(///,' THIS IS A LEADING ORDER FIT')
      GO TO 45
   43 WRITE(6,46)
   46 FORMAT(///,' THIS IS A NEXT TO LEADING ORDER FIT')
   45 CONTINUE
C 
C  THE ROUTINE TEST PERFORMS VARIOUS TESTS ON THE ALTARELLI-
C  PARISI EQUATION INTEGRATOR.
C 
C      CALL TEST(CURPAR)
C
C  CLOSE DATA FILES
C
      CLOSE(7)
      WRITE(6,71) Q02
   71 FORMAT(/,' Q02=',F8.2)
c
c  prepare for fitting
c
      FIT=0.
      DO 66 IF=1,NPARAM 
   66 FIT=FIT+PWATE(IF)
      if(fit.eq.0.d0)then
         V(1)=DATA(3,1)
         V(2)=DATA(4,1)
         v(3)=data(5,1)
         V(4)=DATA(6,1)
         DUM=THEORY(1,1,V,CURPAR)
         WRITE(6,74)
 74      FORMAT(/,' INPUT PARAMETERS') 
         DO 72 L=1,NPARAM
 72      WRITE(6,73) IPNAME(L),CURPAR(L)
 73      FORMAT(1H ,A10,5X,E15.4)
      else
         CALL MINIM(NV,NPTS,DATA,IDRV) 
      endif
C 
C  PRINT OUTPUT NOW.
C  JT=1....PROTON TARGET F2 (MUON OR ELECTRON)
C    =2....ISOSCALAR TARGET F2 (MUON OR ELECTRON)
C    =3....F2 NEUTRINO PROTON
C    =4....F2 NEUTRINO NEUTRON
C    =5....ISOSCALAR TARGET F2 NEUTRINO
C    =6....XF3 NEUTRINO PROTON
C    =7....XF3 NEUTRINO NEUTRON
C    =8....ISOSCALAR TARGET XF3 NEUTRINO 
C    =9....F2 NEUTRON/PROTON RATIO 
C    =10...NEUTRINO IRON CROSS SECTION
C    =11...ANTINEUTRINO IRON CROSS SECTION
C
      CHITOT=0. 
      NEXP=0
      CHITEMP=0.
      schi=0.
      resid=0.
      oldit=0
      DO 80 K=1,NDIS      
      IT=itype(DATA(6,K))
      if(it.ne.oldit)then
         if(oldit.ne.0)then
            rat=chitemp/nexp
            write(6,27) chitemp,nexp,rat
            write(6,28) schi
            write(6,29) resid
            chitot=chitot+chitemp
         endif
         oldit=it
         write(6,83)
         chitemp=0.
         schi=0.
         resid=0.
         nexp=0
      endif
      nexp=nexp+1
      V(1)=DATA(3,K)
      V(2)=DATA(4,K)
      v(3)=data(5,k)
      V(4)=data(6,k)
      TH=THEORY(1,2,V,CURPAR)
      CHISQ=((DATA(1,k)-TH)/DATA(2,k))**2
      CHITEMP=CHITEMP+CHISQ
      IS=-1
      IF(TH.GE.DATA(1,K)) IS=1
      CHISQ=IS*CHISQ
      schi=schi+chisq
      resid=resid+(data(1,k)-th)*th/data(2,k)**2
      write(6,84) data(3,k),data(4,k),th,data(1,k),data(2,k),
     2chisq,ititle(data(6,k))
 80   continue
 83      FORMAT(///,5X,' X',8X,' Q2',5X,'THEORY',6X,'DATA',5X,
     2   'ERROR',2X,'CHI SQUARE',3X,'ITYPE')
 84      FORMAT(f10.6,5F10.4,5X,a10)
 27   FORMAT(///,' CHI SQUARE=',F8.2,' FOR ',I4,' POINTS',
     2'   CHISQ/NPTS=',f8.3)
 28    format(' Signed chisquare total =',f8.2)
 29    format(' Modified residual =',f8.2)
c
c  write out last chi square entry for dis section
c
      if(ndis.ne.0)then
         rat=chitemp/nexp
         write(6,27) chitemp,nexp,rat
         write(6,28) schi
         write(6,29) resid
         chitot=chitot+chitemp
      endif

c
c  iout=1...m^3 dsig/dm/dy
c  iout=2...m^3 dsig/dm/dxf
c  iout=3...(pp-pn)/(pp+pn)
c  iout=4...pd/2pp
c  iout=5...W asymmetry
c
      ND=NDIS+1
      NDU=ndis+ndy
      DO 201 KT=1,5
      DO 203 KRS=1,NRSPOS
      NEXP=0
      CHITEMP=0.
      DO 202 K=ND,NDU 
      IEXP=DATA(6,K)
      IT=ITYPE(IEXP)
      nfl=nflag(iexp)
      IF(nfl.NE.KT) GOTO 202
      IOUT=KT
      TST=RSPOS(KRS)**2 
      IF(abs(DATA(3,K)-TST).gt..1) GO TO 202
      NEXP=NEXP+1
      XQ(NEXP)=DATA(4,K)/sqrt(data(3,k))
      YDAT(NEXP)=DATA(1,K)
      ERD(NEXP)=DATA(2,K) 
      ITOUT(NEXP)=DATA(6,K) 
      YOUT(NEXP)=DATA(5,K)
      V(1)=DATA(3,K)
      V(2)=DATA(4,K)
      V(3)=DATA(5,K)
      V(4)=DATA(6,K)
      THOUT(NEXP)=THEORY(1,2,V,CURPAR)
      CHI(NEXP)=((YDAT(NEXP)-THOUT(NEXP))/ERD(NEXP))**2 
      CHITEMP=CHITEMP+CHI(NEXP) 
      IS=-1 
      IF(THOUT(NEXP).GE.YDAT(NEXP)) IS=1
      CHI(NEXP)=IS*CHI(NEXP)
  202 CONTINUE
      IF(NEXP) 222,222,223
  223 CONTINUE
      WRITE(6,224) RSPOS(KRS)
  224 FORMAT(///,10X,' SQRT(S)=',F8.2)
      GO TO (231,232,231,232,231),IOUT
 231  CONTINUE
      index=v(4)
      if(itype(index).eq.108.or.itype(index).eq.110)then
         write(6,244)
      else
         WRITE(6,241)
      endif 
 241  FORMAT(///,' SQRT(TAU)',5X,'Y',6X,'THEORY',6X,'DATA', 
     25X,'ERROR',2X,'CHI SQUARE',3X,'ITYPE')
      WRITE(6,242) (XQ(N),YOUT(N),THOUT(N),YDAT(N),ERD(N),CHI(N),
     2ITITLE(ITOUT(N)),N=1,NEXP)
 242  FORMAT(6F10.4,5X,a10)
      GO TO 234 
 232  CONTINUE
      WRITE(6,244) 
 244  FORMAT(///,' SQRT(TAU)',5X,'XF',5X,'THEORY',6X,'DATA',5X,'ERROR',
     22X,'CHI SQUARE',3X,'ITYPE') 
      WRITE(6,245) (YOUT(N),THOUT(N),YDAT(N),ERD(N),CHI(N),
     2ITITLE(ITOUT(N)),N=1,NEXP)
 245  FORMAT(6F10.4,5X,a10)
  234 CONTINUE
      rat=chitemp/nexp
      WRITE(6,27) CHITEMP,NEXP,rat 
      CHITOT=CHITOT+CHITEMP 
  222 CONTINUE
  203 CONTINUE
  201 CONTINUE
  204 CONTINUE
c
c print out jet output, if any
c
      ND=NDU+1
      NDU=ndis+ndy+njet
      NEXP=0
      CHITEMP=0.
      oldit=0.
      DO 302 K=ND,NDU 
      IEXP=DATA(6,K)
      IT=ITYPE(IEXP)
      if(oldit.ne.it)then
         if(oldit.ne.0)then
            rat=chitemp/nexp
            write(6,27) chitemp,nexp,rat
            if(oldit.eq.201.and.ncov.ne.0)then
               write(6,304)cdfchisqr
               chitot=chitot+cdfchisqr
            else
               chitot=chitot+chitemp
            endif
         endif
         oldit=it
         write(6,303)
 303     format(///,3x,' ymin',5x,'ymax ',5x,'  ET ',6x,'Theory',8x,
     $'Data',7x,'Error',3x,'Chi Square',3x,'Itype')       
         chitemp=0.
         nexp=0 
      endif
      NEXP=NEXP+1
      V(1)=DATA(3,K)
      V(2)=DATA(4,K)
      V(3)=DATA(5,K)
      V(4)=DATA(6,K)
      TH=THEORY(1,2,V,CURPAR)
      CHISQ=((Data(1,k)-TH)/Data(2,k))**2 
      CHITEMP=CHITEMP+CHISQ
      IS=-1 
      IF(TH.GE.Data(1,k)) IS=1
      CHISQ=IS*CHISQ
      WRITE(6,342) data(3,k),data(4,k),data(5,k),TH,data(1,k),
     2data(2,k),CHISQ,ITITLE(data(6,k))
 342  FORMAT(3F10.4,3e12.4,f10.4,5X,a10)
  302 CONTINUE
      if(oldit.ne.0.)then
         rat=chitemp/nexp
         write(6,27) chitemp,nexp,rat
         if(oldit.eq.201.and.ncov.ne.0)then
            write(6,304)cdfchisqr
 304     format(1h ,'cdf chi square with correlations = ',f10.4)
            chitot=chitot+cdfchisqr
         else
            chitot=chitot+chitemp
         endif
      endif
c
c print out gamma + jet output, if any
c
      ND=NDU+1
      NDU=ndis+ndy+njet+ngamjet
      NEXP=0
      CHITEMP=0.
      oldit=0.
      DO 402 K=ND,NDU 
      IEXP=DATA(6,K)
      IT=ITYPE(IEXP)
      if(oldit.ne.it)then
         if(oldit.ne.0)then
            rat=chitemp/nexp
            write(6,27) chitemp,nexp,rat
            chitot=chitot+chitemp
         endif
         oldit=it
         write(6,403)
 403     format(///,3x,'  rs ',5x,' pt  ',6x,'Theory',8x,
     $'Data',7x,'Error',3x,'Chi Square',3x,'Itype')       
         chitemp=0.
         nexp=0 
      endif
      NEXP=NEXP+1
      V(1)=DATA(3,K)
      V(2)=DATA(4,K)
      V(3)=DATA(5,K)
      V(4)=DATA(6,K)
      TH=THEORY(1,2,V,CURPAR)
      CHISQ=((Data(1,k)-TH)/Data(2,k))**2 
      CHITEMP=CHITEMP+CHISQ
      IS=-1 
      IF(TH.GE.Data(1,k)) IS=1
      CHISQ=IS*CHISQ
      WRITE(6,442) data(3,k),data(4,k),TH,data(1,k),
     2data(2,k),CHISQ,ITITLE(data(6,k))
 442  FORMAT(2F10.4,3e12.4,f10.4,5X,a10)
  402 CONTINUE
      if(oldit.ne.0.)then
         rat=chitemp/nexp
         write(6,27) chitemp,nexp,rat
         chitot=chitot+chitemp
      endif
c
c  write out the normalization chi square contributions
c
      write(6,314)
 314  format(//,1x,'data set',5x,'normalization',8x,'error',
     2     5x,'chi square')
      do 311 i=1,nfit
      if(errn(i).eq.0.d0.or.inorm(i).eq.60)goto311
      chitemp=((1.d0-curpar(inorm(i)))/errn(i))**2
      write(6,312) ititle(i),curpar(inorm(i)),errn(i),chitemp
 312  format(A10,3f15.4)
      chitot=chitot+chitemp
 311  continue
      WRITE(6,2001) CHITOT,NPTS
   10 CONTINUE
 2003 FORMAT('   LAMBDA =',F10.3) 
 2001 FORMAT('  TOTAL CHI-SQUARE =',F10.3,' FOR ',I4,' POINTS') 
 1000 FORMAT(I10,5F10.0,10X,A10)
 1001 FORMAT('  ITYPE=',I2,'  X=',F5.3,'  Q2=',F7.2,'  DATA=',
     $ F7.4,'  ERROR=',F7.4,'  LABEL=',A10) 
 1002 FORMAT('  THERE ARE',I5,' DATA POINTS ENTERED') 
      Q2=10.
c      IF(INS.EQ.1) GOTO 5000
      DO 5000 J=1,2 
      WRITE(6,5003) Q2 
      WRITE(6,5004)
      ALS0=DLOG(Q02/CURPAR(1)**2) 
      ALQ=DLOG(Q2/CURPAR(1)**2) 
      ALS=DLOG(ALQ/ALS0)
      DO 5001 K=1,19
      X=.05*K 
      call pdf(0,x,als,u,d,ub,db,sb,cb,bb,glue)
      dou=d/u
      WRITE(6,5002) X,GLUE,U,D,SB,UB,DB,CB,BB,dou
 5002 FORMAT(F6.2,9E12.4) 
 5001 CONTINUE
      Q2=64. 
 5000 CONTINUE
 5003 FORMAT(//,' Q2=',F6.2)
 5004 FORMAT(/,2X,' X',7X,'XG',10X,'XU',10X,'XD',10X,'XS',
     210X,'XUB',9X,'XDB',9X,'XC ',9X,'XB',8x,'d/u') 
      x=.1
c      IF(INS.EQ.1) GOTO 5000
      DO 6000 J=1,2 
      WRITE(6,6003) x 
      WRITE(6,6004)
      DO 6001 K=1,21
      Q2=15.+(k-1)*1.
      ALS0=DLOG(Q02/CURPAR(1)**2) 
      ALQ=DLOG(Q2/CURPAR(1)**2) 
      ALS=DLOG(ALQ/ALS0)
      call pdf(0,x,als,u,d,ub,db,sb,cb,bb,glue)
      dou=d/u
      WRITE(6,6002) Q2,GLUE,U,D,SB,UB,DB,CB,BB,dou
 6002 FORMAT(F6.2,9E12.4) 
 6001 CONTINUE
      x=.3 
 6000 CONTINUE
 6003 FORMAT(//,' x=',F6.2)
 6004 FORMAT(/,2X,' Q2',7X,'XG',10X,'XU',10X,'XD',10X,'XS',
     210X,'XUB',9X,'XDB',9X,'XC ',9X,'XB',8x,'d/u') 
      call pdfband(qsq,iband)
C
C  DRELL-YAN TEST SECTION
C  PRINTS ONLY IF E-605 DATA ARE IN THE NLO FIT
C
c      iend=0
c      do 6005 j=1,nfit
c      if(itype(j).eq.106.and.iord.eq.1)then
c          v(4)=j
c          iend=1
c      endif
c 6005 continue
c      if(iend.eq.0)call exit
c      RS=38.8
c      Y=0.
c      WRITE(6,6001) RS,Y
c 6001 FORMAT(/' DRELL-YAN TEST AT SQRT(S)=',F6.2,' GEV AND Y=',F5.2)
c      WRITE(6,6004)
c 6004 FORMAT(/' SQRT(TAU)',4X,'ANSWER',11x,'Born',11X,'RATIO',10X,
c     2'QUARK',10X,'GLUON')
c      DO 6002 I=2,9
c      RTAU=.05*I
c      XM=RTAU*RS
c      V(1)=RS**2
c      V(2)=XM
c      V(3)=Y
c      index=v(4)
c      in=inorm(index)
c      ANS=THEORY(1,2,V,CURPAR)
c      fac=rtau**3
c      xseclo=xseclo*fac/curpar(in)
c      glucor=glucor*fac/curpar(in)
c      AKFAC=ANS/XSECLO
c      HO=ANS-XSECLO
c      QCOR=HO-GLUCOR
c      AKQ=QCOR/XSECLO
c      AKG=GLUCOR/XSECLO
c      WRITE(6,6003) RTAU,ANS,xseclo,AKFAC,AKQ,AKG
c 6002 CONTINUE
c 6003 FORMAT(F10.2,5E15.4)
      STOP
      END 
      subroutine pdfband(qsq,iband)
      implicit real*8 (a-h,o-z)
      character*30, output, file1
      character*4, tmp
      dimension derv(-5:2,50,35), err(-5:2,50), pdf0(-5:2,50), xgrid(50)
      common/gint/gf(8,40,1800)
      COMMON/MINIMB/CURPAR(60),UNCRT(60),ERMIN(40),PWATE(60),IFREP(40),
     * NPARAM,NVAR,IREJ 
      COMMON/Q2STUF/ Q02,Q2MAX
      COMMON/INPUT/INS,NVL,NMAX,DELTA,S0,IORD
      DATA NX,XGRID/50,.0001,.0002,.0004,.0006,.0008,.001,.002,.004,
     2.008,.016,.032,.064,.1,.13,.16,.19,.22,.25,.28,.31,.34,.37,.40,
     3.43,.46,.49,.52,.55,.58,.61,.64,.67,.70,.73,.76,.79,.82,.84,.86,
     4.88,.90,.91,.92,.93,.94,.95,.96,.97,.98,.99/
      common/outputfile/output      
      common/errors/cov(40,40)
      ALS0=DLOG(Q02/CURPAR(1)**2) 
      ALQ=DLOG(QSQ/CURPAR(1)**2) 
      ALS=DLOG(ALQ/ALS0)
      call trmstr(output,len1)
      tmp='.pdf'
      call trmstr(tmp,len2)
      len1=len1-4
      file1=output(1:len1)//tmp(1:len2)
      len=len1+len2
      open(unit=3,file=file1(1:len),status='unknown')
      do 100 j=1,nx
      call pdf(0,xgrid(j),als,u0,d0,ub0,db0,sb0,cb0,bb0,glue0)
      pdf0(-5,j)=bb0
      pdf0(-4,j)=cb0
      pdf0(-3,j)=sb0
      pdf0(-2,j)=db0
      pdf0(-1,j)=ub0
      pdf0(-0,j)=glue0
      pdf0(1,j)=u0
      pdf0(2,j)=d0
      do 101 k=1,nvar
      call pdf(k,xgrid(j),als,u,d,ub,db,sb,cb,bb,glue)
      kk=ifrep(k)
      derv(-5,j,k)=(bb-bb0)/pwate(kk)
      derv(-4,j,k)=(cb-cb0)/pwate(kk)
      derv(-3,j,k)=(sb-sb0)/pwate(kk)
      derv(-2,j,k)=(db-db0)/pwate(kk)
      derv(-1,j,k)=(ub-ub0)/pwate(kk)
      derv(0,j,k)=(glue-glue0)/pwate(kk)
      derv(1,j,k)=(u-u0)/pwate(kk)
      derv(2,j,k)=(d-d0)/pwate(kk)
 101  continue
      do 102 i=-5,2
      err(i,j)=0.d0
      do 103 l=1,nvar
      do 103 m=1,nvar
      err(i,j)=err(i,j)+derv(i,j,l)*cov(l,m)*derv(i,j,m)
 103  continue
      err(i,j)=sqrt(err(i,j))
 102  continue
 100  continue
      do 200 i=-5,2
      write(3,210) i,qsq
 210  format(/ ,'parton flavor =',i5,' at Q**2 = ',f10.2,' GeV**2')
      write(2,211)
 211  format(/ ,5x,' x',5x,'pdf',5x,'error')
      do 201 j=1,nx
      if (iband.eq.1)then
         pdfout=pdf0(i,j)
         errout=err(i,j)
      else if(iband.eq.2)then
         pdfout=1.d0
         errout=err(i,j)/pdf0(i,j)
      endif
      write(3,212) xgrid(j),pdfout,errout
 212  format(f10.4,2e12.4)
 201  continue
 200  continue
      call exit
      end

c                     The function below "EMC(x, iver)" gives the
c	Fe/D correction for the "x" bin. Note that the
c	.07< x < 1 is what you must use this function for.
	
c                   For x=.045, the Fe/D is 0.95. It is the
c	"shadowing" region and the factor was obtained for
c	Ca/D, and correcting Ca to Fe using a prescription
c	by Strikman.

	function emcslac(f,x,q2)

c two versions of SLAC fits for the F2(Fe)/F2(D2) ratio

	implicit double precision (a-z)

c version 2 is the fit to SLAC E-139 and E-140 and STEIN

	data d0, d1, d2, d3, d4, d5, d6, d7, d8
     $  /4.58558707D-01,
     $   1.62185596D+01,
     $  -1.79392859D+02,
     $   1.04313998D+03,
     $  -3.53408342D+03,
     $   7.17002801D+03,
     $  -8.56431003D+03,
     $   5.54039709D+03,
     $  -1.49318167D+03/

        x2 = x*x
        x4 = x2*x2
        x8 = x4*x4

        if(x.ge..07) then
           emc1 = d0 + d1*x + d2*x2 + d3*x*x2 + d4*x4
     $        +  d5*x4*x + d6*x4*x2 + d7*x4*x2*x + d8*x8
        else
c
c the following form is a linear interpolation between the x = .07
c value of the above polynomial and the Strikman nuclear shadowing
c calculation - emc1 = .95 @ x = .045
c
           emc1 = .95 + (x - .045)*1.9544

        endif
        emcslac = f / emc1

	return
	end

      FUNCTION F2DTOA(F,X,Q2)
 

      Implicit Double Precision (A-H, O-Z)
 

C *** GIVEN F2 MEASURED ON FE OR C AT A GIVEN X IT WILL GIVE THE
C *** EQUIVALENT DEUTERIUM VALUE.  IT WILL ALSO CALCULATE THE POSSIBLE
C *** Q**2 DEPEndENCE OF THE EFFECT FOR IRON IF THE VALUE OF Q GIVEN TO
C *** FUNCTION IS GREATER THAN SQRT(5).
 

C     DATA Q20 / 5.0 /
 

      X2 = X * X
      X4 = X2 * X2
 

      DF2 = 1.18 - 4.02*X + 24.35*X2 - 61.5*X2*X + 72.61*X4 -32.35*X4*X
      F2D = DF2 * F
 

C      IF(Q2 .LE. Q20) GO TO 20
C      DQ2 = Q2 / 5.0
C      ALDQ2 = LOG(DQ2)
C      F2D = F2D + (0.077 - 0.244*X) * ALDQ2
C  20  Continue
 

      F2DTOA = F2D
      Return
C                        ****************************
      End

      function emcnmc(f,x,q2)

c Concerning the "heavy target correction" I used a parametrisation of
c our measurement of Ca/D2 and the Slac result for Fe/D2.
c As far as I can see it,this should be also valid for the iron data
c of CCFR (from our preliminay data on Sn/C we conclude that the
c amount of shadowing seems to saturate at values of A around 40).
c I include this parametrisation into this mail.


      implicit double precision(a-z)

      emc = 1.118 - 0.4199*x - 0.3597*exp(-22.88*x) +
     &        1.872*(x**11.27)

      emcnmc = f /  emc

      return
      end

      function gomez(f,x)
      implicit double precision (a-h,o-z)
      data p1,p2,p3,p4,p5,p6,p7/1.0164d0,-.047762d0,-.13354d0,
     2.35303d0,.22719d0,-1.2906d0,5.6075d0/
c
c  Fit to the data for F2D/F2N where N=(p+n)/2
c  as extractred by Gomez et al PRD 49, 4348 (1994)
c
      THEORY=P1+P2*X+P3*x**2+P4*x**3+P5*x**4
     2+P6*x**5
      deufac=p7*(1./x-1.)
      if(deufac.ge.20.)deufac=20.
      deucor=1.-exp(-deufac)
      theory=theory/deucor
      gomez=f/theory
      RETURN
      END



      function emce665(f,x,q2)

      implicit double precision (a-z)

c Data from E-665 (Xenon) and NMC (Ca scaled) used for large x fit
c and E-665 for very small x.  The small x fit went to x = 10**-5

      x2 = x*x
      x4 = x2*x2

      if(x.ge..01) then
         rat = 0.73 + 7.0*x - 52.9*x2 + 173.8*x2*x - 263.9*x4 +
     $        150.*x4*x
      else
         rat = 0.71 + 6.8*x - 31.1*x2
      endif

      emce665 = f / rat

      return

      end
      subroutine fildat(nfit,npts,data,q2max)
      implicit real*8 (a-h,o-z)
      character*80 title(4)
      character*10 ititle(30),fln,tmp
      character*30 file1
      dimension data(6,2500), ifit(30),q2cut(30),xlo(30),xhi(30),
     $w2min(30),errn(30),dum(6)
      common/cuts/q2cut,xlo,xhi,w2min
      common/flags/itype(30),nflag(30),inorm(30),itgt(30),icorr(30)
     2,itmc(30)
      common/expts/ititle,errn,iserr
      common/Onejetset/Jset 
      common/DeltaYmax/DelYmax
      common/cdf/covinv(33,33),cdfchisqr,ncdf,ncdfl,ncov
      data fln/'../data/'/
      Jset=1
      DelYmax=7.0d0
      ncdf=0
      ncdfl=0
      ncov=0
      npts=0
      q2max=0.d0
      do 100 j=1,nfit
      call trmstr(fln,len1)
      tmp=ititle(j)
      call trmstr(tmp,len2)
      file1=fln(1:len1)//tmp(1:len2)
      len=len1+len2
      open(unit=8,file=file1(1:len),status='old')
c      print*,'opening file=',file1(1:len),j,itype(j)
      do 99 l=1,3
 99   read(8,*) title(l)
c      read(8,*,err=98) errn(j),ncorr
c
c  temporarily disabled ncorr because of format problems in some data files
c
      read(8,*,err=98) errn(j)
c      goto 98
 97   ncorr=0.
 98   continue
      read(8,*) title(4)
 200  continue
      read(8,*,end=96) (dum(jk),jk=1,6),ak0,ak1
      if(iserr.eq.0)then
         err=dum(5)
      else if(icorr(j).eq.0)then
         err=sqrt(dum(5)**2+dum(6)**2)
      else if(icorr(j).ne.0)then
         err=dum(5)
      endif
      if(itype(j).lt.70)then
         x=dum(1)
         q2=dum(2)
         dat=dum(3)
         if(itgt(j).eq.1)then
            if(nflag(j).eq.5.or.nflag(j).eq.8.)then
               dat=emcnmc(dat,x,q2)
               err=emcnmc(err,x,q2)
            else if(nflag(j).eq.2.or.nflag(j).eq.9)then
               dat=gomez(dat,x)
               err=gomez(err,x)
            endif            
         else if(itgt(j).eq.2)then
            dat=emcnmc(dat,x,q2)
            err=emcnmc(err,x,q2)
            dat=gomez(dat,x)
            err=gomez(err,x)
         endif
         if(x.lt.xlo(j))goto 200
         if(x.gt.xhi(j))goto 200
         if(q2.lt.q2cut(j))goto 200
         w2=q2*(1.d0/x-1.d0)+.88d0
         if(w2.lt.w2min(j))goto 200
         npts=npts+1
         data(1,npts)=dat
         data(2,npts)=err
         data(3,npts)=x
         data(4,npts)=q2
         data(5,npts)=w2
         data(6,npts)=j
      else if(itype(j).eq.70.or.itype(j).eq.71)then
c
c  nutev iron (anti)neutrino cross section data
c
c         e=dum(1)
c         if(e.eq.85..or.e.eq.170.)then
c            continue
c         else
c            goto 200
c         endif
c         x=dum(2)
c         y=dum(3)         
c
c  modified to read ordered nutev data files
c
         x=dum(1)
         q2=dum(2)
         dat=dum(4)
         y=dum(3)
c         q2=2.*.938*e*x*y
         if(itgt(j).eq.1)then
            dat=emcnmc(dat,x,q2)
            err=emcnmc(err,x,q2)
         else if(itgt(j).eq.2)then
            dat=emcnmc(dat,x,q2)
            err=emcnmc(err,x,q2)
            dat=gomez(dat,x)
            err=gomez(err,x)
         endif
         if(x.lt.xlo(j))goto 200
         if(x.gt.xhi(j))goto 200
         if(q2.lt.q2cut(j))goto 200
         w2=q2*(1.d0/x-1.d0)+.88d0
         if(w2.lt.w2min(j))goto 200
         npts=npts+1
         data(1,npts)=dat
         data(2,npts)=err
         data(3,npts)=x
         data(4,npts)=q2
         data(5,npts)=y
         data(6,npts)=j
      else if(itype(j).gt.100.and.itype(j).le.200)then
         rs=dum(1)
         if(itype(j).eq.133)then
            x2=dum(2)
            xf=dum(3)
            x1=xf+x2
            xm=sqrt(x1*x2)*rs
            y=xf
         else
            rtau=dum(2)
            xm=rs*rtau
            y=dum(3)
         endif
c         if(itype(j).eq.106)y=dum(3)
         dat=dum(4)
         if(xm**2.lt.q2cut(j))goto 200
         npts=npts+1
         data(1,npts)=dat
         data(2,npts)=err
         data(3,npts)=rs**2
         data(4,npts)=xm
         data(5,npts)=y
         data(6,npts)=j
         q2=xm**2
      else if(itype(j).gt.200.and.itype(j).le.300)then
         y1min=dum(1)
         y1max=dum(2)
c
c  adjust y1min, y1max  for CDF data set since first entry is sqrt(s)
c
         if(itype(j).eq.201)then
            y1min=0.1d0
            y1max=0.7d0
            if(ncdf.eq.0)ncdfl=npts+1
            ncdf=ncdf+1
         endif
         pt=dum(3)
         dat=dum(4)
         npts=npts+1
         data(1,npts)=dat
         data(2,npts)=err
         data(3,npts)=y1min
         data(4,npts)=y1max
         data(5,npts)=pt
         data(6,npts)=j
         q2=pt**2
      else if(itype(j).gt.300.and.itype(j).lt.400)then
         rs=dum(1)
         pt=dum(2)
         dat=dum(4)
         npts=npts+1
         data(1,npts)=dat
         data(2,npts)=err
         data(3,npts)=rs
         data(4,npts)=pt
         data(5,npts)=itype(j)-300
         data(6,npts)=j
         q2=4.*pt**2
      endif
      if(q2.gt.q2max)q2max=q2
      goto 200
 96   continue
      close(8)
 100  continue
      njet=0
      do 101 j=1,nfit
      if(itype(j).gt.200.and.itype(j).lt.300)njet=1
 101  continue
c
c  if there is at least one jet data set, read in the jet kfactor table
c
      if(njet.ne.0)then
         call trmstr(fln,len1)
         file1=fln(1:len1)//'kf1jet.msb'
         call trmstr(file1,len1)
         call read1jet(file1(1:len1))
      endif      
c
c  read CDF covariance matrix if needed
c
      if(ncdf.ne.0)then
         do j=1,nfit
         if(itype(j).eq.201.and.icorr(j).eq.1)then
           ncov=1
           ConvUnit=1D00
           IU=NextUn()
           Call TrmStr(fln,  Len)
           file1= fln(1:Len)//'cdf_01_07_covmtx.dat'
           OPEN(IU,file=file1,status='old',err=999)
           READ(IU,*) n
           If(n.ne.ncdf) 
     2     stop '# of jet data does not match Error Matrix.'
           DO k = 1,n*n
             READ (IU,*) i,l,c,e
             CovInv(i,l) = e * ConvUnit
           ENDDO
           CLOSE(IU)
           print*,'CDF covariance matrix read in'
         endif
         enddo
      endif
      return
 999  Print*,'Cannot open file:',file1
      stop
C                      *******************************       

      end
      subroutine rdxpos(xpos,nxpos,rspos,nrspos,npts,data,ndis,ndy,njet,
     2ngamjet)
      implicit real*8 (a-h,o-z)
      dimension data(6,2500),xpos(350),rspos(15)
      common/flags/itype(30),nflag(30),inorm(30),itgt(30),icorr(30)
     2,itmc(30)
      nxpos=0
      nrspos=0
      ndis=0
      ndy=0
      njet=0
      ngamjet=0
      do 100 j=1,npts
      jt=data(6,j)
      it=itype(jt)
      if(it.le.100)then
         x=data(3,j)
         ndis=ndis+1
         if(nxpos.eq.0)then
            nxpos=1
            xpos(nxpos)=x
         else
            do 90 l=1,nxpos
            if(x.eq.xpos(l))goto 91
 90         continue
            nxpos=nxpos+1
            if(nxpos.gt.npts)return
            xpos(nxpos)=x
 91         continue
         endif
      else if(it.gt.100.and.it.le.200)then
         ndy=ndy+1
         rs=sqrt(data(3,j))
         if(nrspos.eq.0)then
            nrspos=1
            rspos(nrspos)=rs
         else
            do 80 l=1,nrspos
            if(rs.eq.rspos(l))goto 81
 80         continue
            nrspos=nrspos+1
            rspos(nrspos)=rs
 81         continue
         endif
      else if(it.gt.200.and.it.lt.300)then
         njet=njet+1
      else if(it.gt.300.and.it.lt.400)then
         ngamjet=ngamjet+1
      endif
 100  continue
c      print*,nxpos,xpos
      return
      end
