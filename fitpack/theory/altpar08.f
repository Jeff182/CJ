************************************************************************
*     ALTPAR08                                                         *
*                                                                      * 
*     PDF parametrization and evolution                                *
*                                                                      *
*     HISTORY:                                                         *
*                                                                      *
*     altpar08   : based on altpar02                                   *
*     (12 Nov 08)  modifeis the d-quark parametrization to             *
*                  allow d/u->finite as x->1                           *
*                  NOTE: requirs a new d-quark parameter in position   *
*                  #32, called 'a6dv'.                                 *
*                                                                      *
************************************************************************


      SUBROUTINE INTQCD(XPAR) 
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION FVL(30),FVL1(30),FVLSTO(30),PRD(30),PRD1(30),
     $ COR(30),COR1(30),FUDGE(30),PRDM(30)
      DIMENSION FS(30),FS1(30),FSSTO(30),PRDS(30),PRDS1(30),CORS(30), 
     $ CORS1(30),FUDGES(30),PRDMS(30)
      DIMENSION FG(30),FG1(30),FGSTO(30),PRDG(30),PRDG1(30),CORG(30), 
     $ CORG1(30),FUDGEG(30),PRDMG(30)
      DIMENSION XPAR(60)
      COMMON/CONSTANTS/PI,PI2
      COMMON/ALQ2/T 
      COMMON/GRPTHY/FLAVOR
      COMMON/PARAM/PARA(60)
      COMMON/GAUS16/XI(16),WI(16),NTERMS,XX(17)
      COMMON/INPUT/INS,NVL,NMAX,DELTA,S0,IORD
      COMMON/GFUNC/CALC(8,40,30)
      COMMON/GRID/NX,XGRID(30)
      common/threshold/sb
      DATA NX,XGRID/30,.0001,.0002,.0004,.0006,.0008,.001,.002,.004,
     2.008,.016,.032,.064,.1,.15,.2,.25,.3,.35,.4,.45,.5,.55,.6,.65,
     3.7,.75,.80,.85,.9,.95/
C  GFUNC IS THE ARRAY WHICH HOLDS THE EVOLVED DISTRIBUTIONS.
C  ALLOWANCE HAS BEEN MADE FOR EIGHT PARTON DISTRIBUTIONS.
C  THE SECOND DIMENSION IS THE NUMBER OF STEPS IN S
C  S=LN(LN(Q2/LAMBDA2)/LN(Q02/LAMBDA2)) 
C  THE THIRD DIMENSION IS THE NUMBER OF X POINTS POINTS IN THE GRID 
      DO 61 I=1,8
      DO 61 J=1,40
      DO 61 K=1,30
   61 CALC(I,J,K)=0.D0
C  ALLOWANCE HAS BEEN MADE FOR UP TO 60 PARAMETERS
      DO 48 J=1,60
   48 PARA(J)=XPAR(J)
      PI=4.*DATAN(1.D0)
      PI2=PI**2
c
c  set b threshold
c
      sb=dlog(dlog(4.5**2/para(1)**2)/s0)
      ith=0
c
c  renormalize uv, dv, g to satisfy sum rules
c
c  check nmax
c
      if(nmax.gt.40)then
         print*,'nmax=',nmax
         nmax=40
         print*,'nmax reset to 40'
      endif
      para(2)=1.
      para(8)=1.
      para(25)=1.
      call pdfnrm
C  CALCULATE THE SINGLET (FS) AND GLUON (FG)
C  DISTRIBUTIONS FIRST UNLESS INS IS GREATER
C  THAN ZERO
      S=0.0
C**  CALCULATE FS AND FG AT S=0.0
      IF(INS) 30,30,31
   30 CONTINUE
      DO 1 I=1,NX
      X=XGRID(I)
      FS(I)=FCNFS(X)
      FG(I)=FCNFG(X)
      CALC(7,1,I)=FS(I)
      CALC(8,1,I)=FG(I)
    1 CONTINUE
c      print2001,(calc(7,1,i),i=1,nx)
c 2001 format(7e11.2)
      tst=calc(7,1,5)
C**  CALCULATE D(FS,FG)/DS AT S=0.0
      T=S0
      DO 2 I=1,NX
      X=XGRID(I)
      CALL FINTGS(X,FS,FG,FS1(I))
      CALL FINTGG(X,FG,FS,FG1(I))
    2 CONTINUE
c      print2001,fs1
C**  CALCULATE FS AND FG AT S=0.0-DELTA 
      T=S0*DEXP(-DELTA)
      DLT=DELTA
      CALL STARTS(FS,FS1,FSSTO,FG,FG1,FGSTO,DLT)
C**  CALCULATE PREDICTOR AT S=0.0+DELTA 
      DO 3 I=1,NX
      PRDS(I)=FSSTO(I)+2.*DELTA*FS1(I)
      PRDG(I)=FGSTO(I)+2.*DELTA*FG1(I)
    3 CONTINUE
C**  CALCULATE D(PRD)/DS AT S=0.0+DELTA 
      T=S0*DEXP(DELTA)
      DO 4 I=1,NX
      X=XGRID(I)
      CALL FINTGS(X,PRDS,PRDG,PRDS1(I)) 
      CALL FINTGG(X,PRDG,PRDS,PRDG1(I)) 
C**  CALCULATE CORRECTOR AT S=0.0+DELTA 
      CORS(I)=FS(I)+0.5*DELTA*(FS1(I)+PRDS1(I))
      CORG(I)=FG(I)+0.5*DELTA*(FG1(I)+PRDG1(I))
      FUDGES(I)=PRDS(I)-CORS(I)
      FUDGEG(I)=PRDG(I)-CORG(I)
    4 CONTINUE
      DO 6 I=1,NX
      X=XGRID(I)
      CALL FINTGS(X,CORS,CORG,CORS1(I)) 
      CALL FINTGG(X,CORG,CORS,CORG1(I)) 
    6 CONTINUE
C**  RESHUFFLE THE DECK
      DO 7 I=1,NX
      FSSTO(I)=FS(I)
      FGSTO(I)=FG(I)
      FS(I)=CORS(I) 
      FG(I)=CORG(I) 
      FS1(I)=CORS1(I)
      FG1(I)=CORG1(I)
      CALC(7,2,I)=FS(I)
      CALC(8,2,I)=FG(I)
    7 CONTINUE
      DO 10 NTIMES=2,NMAX-1
      NT1=NTIMES+1
      DO 11 I=1,NX
      PRDS(I)=FSSTO(I)+2.*DELTA*FS1(I)
      PRDG(I)=FGSTO(I)+2.*DELTA*FG1(I)
      PRDMS(I)=PRDS(I)-0.8*FUDGES(I)
      PRDMG(I)=PRDG(I)-0.8*FUDGEG(I)
   11 CONTINUE
      T=T*DEXP(DELTA)
      DO 12 I=1,NX
      X=XGRID(I)
      FSSTO(I)=FS(I)
      FGSTO(I)=FG(I)
      CALL FINTGS(X,PRDMS,PRDMG,PRDS1(I))
      CALL FINTGG(X,PRDMG,PRDMS,PRDG1(I))
      CORS(I)=FS(I)+0.5*DELTA*(FS1(I)+PRDS1(I))
      CORG(I)=FG(I)+0.5*DELTA*(FG1(I)+PRDG1(I))
      FUDGES(I)=PRDS(I)-CORS(I)
      FUDGEG(I)=PRDG(I)-CORG(I)
      FS(I)=CORS(I)+0.2*FUDGES(I)
      FG(I)=CORG(I)+0.2*FUDGEG(I)
      CALC(7,NT1,I)=FS(I)
      CALC(8,NT1,I)=FG(I)
   12 CONTINUE
      if(tst.ne.calc(7,1,5))then
         print*,ntimes, nmax
         stop
      endif
      DO 14 I=1,NX
      X=XGRID(I)
      CALL FINTGS(X,FS,FG,FS1(I))
      CALL FINTGG(X,FG,FS,FG1(I))
   14 CONTINUE      
      S=NTIMES*DELTA
   10 CONTINUE
C**  START THE NONSINGLET CALCULATION
C**  CALCULATE FVL AT S=0.
   31 CONTINUE
      DO 100 IVL=1,NVL
C  IPM IS A PHASE THAT APPEARS IN THE
C  NEXT-TO-LEADING ORDER CALCULATION
C  IPM=-1 FOR THE U AND D VALENCE DISTRIBUTIONS
C  IPM=+1 FOR THE 'PLUS' TYPE DIsTRIBUTIONS DEFINED AS
C         (Q+QBAR)-SINGLET/FLAVORS
      IPM=-1
      IF(IVL.GT.2) IPM=1
      ith=0
      S=0.
      DO 41 I=1,NX
      X=XGRID(I)
      FVL(I)=FCNVL(IVL,X)
      CALC(IVL,1,I)=FVL(I)
   41 CONTINUE
C**  CALCULATE D(FVL)/DS AT S=0.
      T=S0
      DO 42 I=1,NX
      X=XGRID(I)
      CALL FINTGV(X,IPM,FVL,FVL1(I))
   42 CONTINUE
C**  CALCULATE FVL AT S=0.-DELTA
      T=S0*DEXP(-DELTA)
      DLT=DELTA
      CALL STARTV(FVLSTO,DLT,FVL,FVL1,IPM)
C**  CALCULATE PREDICTOR AT S=0.+DELTA
      DO 43 I=1,NX
      PRD(I)=FVLSTO(I)+2.*DELTA*FVL1(I) 
   43 CONTINUE
C**  CALCULATE D(PRD)/DS AT S=0.+DELTA
      T=S0*DEXP(DELTA)
      DO 44 I=1,NX
      X=XGRID(I)
      CALL FINTGV(X,IPM,PRD,PRD1(I))
      COR(I)=FVL(I)+0.5*DELTA*(FVL1(I)+PRD1(I))
      FUDGE(I)=PRD(I)-COR(I)
   44 CONTINUE
      DO 45 I=1,NX
      X=XGRID(I)
      CALL FINTGV(X,IPM,COR,COR1(I))
   45 CONTINUE
      DO 46 I=1,NX
      FVLSTO(I)=FVL(I)
      FVL(I)=COR(I) 
      FVL1(I)=COR1(I)
      CALC(IVL,2,I)=FVL(I)
   46 CONTINUE
C**  INITIALIZATION COMPLETE...ITERATE UNTIL S=SMAX
      DO 47 NTIMES=2,NMAX-1
      NT1=NTIMES+1
      DO 51 I=1,NX
      X=XGRID(I)
      PRD(I)=FVLSTO(I)+2.*DELTA*FVL1(I) 
      PRDM(I)=PRD(I)-0.8*FUDGE(I)
   51 CONTINUE
      T=T*DEXP(DELTA)
      DO 52 I=1,NX
      X=XGRID(I)
      FVLSTO(I)=FVL(I)
      CALL FINTGV(X,IPM,PRDM,PRD1(I))
      COR(I)=FVL(I)+0.5*DELTA*(FVL1(I)+PRD1(I))
      FUDGE(I)=PRD(I)-COR(I)
      FVL(I)=COR(I)+0.2*FUDGE(I)
   52 CONTINUE
c
c  adjust for crossing b threshold
c
      S=NTIMES*DELTA
      if(s.ge.sb.and.ith.eq.0.and.ivl.gt.2)then
         ith=1
         do 71 i=1,nx
         fvl(i)=fvl(i)+1./20.*calc(7,nt1,i)  
         fvlsto(i)=fvlsto(i)+1./20.*calc(7,nt1-1,i)
 71   continue
      endif
      DO 53 I=1,NX
      X=XGRID(I)
      CALL FINTGV(X,IPM,FVL,FVL1(I))
      CALC(IVL,NT1,I)=FVL(I)
   53 CONTINUE
c      S=NTIMES*DELTA
   47 CONTINUE
  100 CONTINUE
      DO 49 J=1,60
   49 XPAR(J)=PARA(J)
      RETURN
      END 



      subroutine pdfnrm
      implicit real*8 (a-h,o-z)
c      common/gaus16/xi(16),wi(16),nterms,xx(17)
      common/gaus32/xi(32),wi(32),nterms,xx(33)
      common/param/xp(60)
c
c  quantum number sum rule for uv, dv normalization
c
      auv=0.
      adv=0.
      do 100 l=1,nterms
      z=.5*(1.+xi(l))
      x=z**3
c
c  integral of uv from x= 0-1 without overall coefficient
c
      auv=auv+.5*wi(l)*fuv(x)/x/xp(2)*3.*z**2
c
c  integral of dv from x= 0-1 without overall coefficient
c
      adv=adv+.5*wi(l)*fdv(x)/x/xp(8)*3.*z**2
 100  continue
c
c  renormalize uv, dv to satisfy sum rules
c
      xp(2)=2./auv
      xp(8)=1./adv
c
c  momentum sum rule for gluon normalization
c
      puv=0.
      pdv=0.
      as=0.
      ag=0.
      do 200 l=1,nterms
      z=.5*(1.+xi(l))
      x=z**3
c
c  integral of x(singlet)
c
      as=as+.5*wi(l)*fcnfs(x)*3.*z**2
c
c  integrals of xuv and xdv
c
      puv=puv+.5*wi(l)*fuv(x)*3.*z**2
      pdv=pdv+.5*wi(l)*fdv(x)*3.*z**2
c
c  integral of xgluon without overall coefficient
c
      ag=ag+.5*wi(l)*fcnfg(x)/xp(25)*3.*z**2
 200  continue
c
c  renormalize gluon distribution
c
      xp(25)=(1.-as)/ag
      psea=as-puv-pdv
      pglu=ag*xp(25)
c      write(6,201)puv,pdv,psea,as,pglu
 201  format(/,'<x>uv=',f8.4,2x,'<x>dv=',f8.4,2x,'<x>sea=',f8.4,
     2     2x,'<x>singlet=',f8.4,2x,'<x>glue=',f8.4)
      return
      end



      FUNCTION FCNVL(IVL,X)
      IMPLICIT REAL*8 (A-H,O-Z) 
C  CALCULATES THE INPUT QUARK DISTRIBUTIONS AT Q0**2
C  THE SINGLET AND GLUON DISTRIBUTIONS ARE HANDLED IN
C  SEPARATE ROUTINES (FCNFS AND FCNFG)
      COMMON/PARAM/XP(60)
      COMMON/GRID/NX,XGRID(30)
      COMMON/GRPTHY/FLAVOR
C  IVL = 1  XUV
C  IVL = 2  XDV
c  IVL = 3  XU+
C  IVL = 4  XD+
C  IVL = 5  XS+
C  IVL = 6  XC+
C
C  Revised 6/3/02 for five flavors
c  Parametrization changed to CTEQ6 style
c  Revised 2/27/03 to implement improved treatment of b threshold
c
      if(ivl.eq.1)then
c
c  xuv
c
         fcnvl=fuv(x)
         return
      else if(ivl.eq.2)then
c
c  xdv
c
         fcnvl=fdv(x)
         return
      else if(ivl.eq.3)then
c
c  x(uv+2.*ubar-.25*singlet)
c
         ubpdb=fubpdb(x)
         dboub=fdboub(x)
         fcnvl=fuv(x)+2.*ubpdb/(1.+dboub)-.25*fcnfs(x)
         return
      else if(ivl.eq.4)then
c
c  x(dv+2.*dbar-.25*singlet)
c
         ubpdb=fubpdb(x)
         dboub=fdboub(x)
         fcnvl=fdv(x)+2.*ubpdb/(1.+1./dboub)-.25*fcnfs(x)
         return
      else if(ivl.eq.5)then
c
c  x(s + sbar) = kappa*(ubar + dbar)
c  x(s+sbar-.25*singlet)
         fcnvl=xp(31)*fubpdb(x)-.25*fcnfs(x)
         return
      else if(ivl.eq.6)then
c
c  xcplus = -1/4*singlet
c
         fcnvl=-.25*fcnfs(x)
         return
      endif
      end



      function fuv(x)
*     Up valence parametrization at Q0
      implicit real*8 (a-h,o-z)
      common/param/xp(60)
      fuv=xp(2)*x**xp(3)*(1.-x)**xp(4)*exp(xp(5)*x)*
     $(1.+exp(xp(6))*x)**xp(7)
      return
      end



      function fdv(x)
*     Down valence -  parametrization at Q0
      implicit real*8 (a-h,o-z)
      common/param/xp(60)
*     d-valence parametrization - generalizes the CTEQ6.1 
*     parametrization by allowing d/u->finite as x->1.
*     Fix xp(32)=0.0 to get back CTEQ6.1.
*     (A.Accardi 11 Nov 2008)
      fdv = xp(8)*x**xp(9) * ((1.-x)**xp(10)+(xp(32)**2)*(1.-x)**xp(4))
     $     *exp(xp(11)*x)*(1.+exp(xp(12))*x)**xp(13)

c --------------------------------------------------
* (12 Nov 08) alternative d-valence modification
* for d/u->finite as x->1. Doesn't seem to work properly.
*      dum = dexp(xp(32)*(1-x))
*      fdv = xp(8)*x**xp(9) 
*     &     * (1.-x)**(dum*xp(10) + (1-dum)*xp(32))
*     $     * exp(xp(11)*x) * (1.+exp(xp(12))*x)**xp(13)

c --------------------------------------------------
* (12 Nov 08) standard CTEQ6.1 d-valence parametrization 
* superseded by the large-x parametrization 08.11 that allows 
* d/u->finite as x->1.  
*
*      fdv=xp(8)*x**xp(9)*(1.-x)**xp(10)*exp(xp(11)*x)*
*     $(1.+exp(xp(12))*x)**xp(13)

c --------------------------------------------------
C MOD for NUTEV Analysis removed 10/19/05
C Fix dv to have same shape as uv
C
c      fdv=xp(8)*x**xp(3)*(1.-x)**xp(4)*exp(xp(5)*x)*
c     $(1.+exp(xp(7))*x)**xp(7)
C

      return
      end



      function fubpdb(x)
*     ubar+dbar -  parametrization at Q0
      implicit real*8 (a-h,o-z)
      common/param/xp(60)
      fubpdb=xp(14)*x**xp(15)*(1.-x)**xp(16)*exp(xp(17)*x)*
     $(1.+exp(xp(18))*x)**xp(19)
      return
      end



      function fdboub(x)
*     dbar/ubar -  parametrization at Q0
      implicit real*8 (a-h,o-z)
      common/param/xp(60)
      fdboub=xp(20)*x**(xp(21))*(1.-x)**xp(22)
     $+(1.+xp(23)*x)*(1.-x)**xp(24)
      return
      end



      function fcnfs(x)
      implicit real*8 (a-h,o-z)
      common/param/xp(60)
c
c  singlet at Q0
c
      fcnfs=fuv(x)+fdv(x)+fubpdb(x)*(2.+xp(31))
      return
      end



      function fcnfg(x)
      implicit real*8 (a-h,o-z)
      common/param/xp(60)
c
c  gluon at Q0
c
      fcnfg=xp(25)*x**xp(26)*(1.-x)**xp(27)*exp(xp(28)*x)*
     $(1.+exp(xp(29))*x)**xp(30)
      return
      end



      SUBROUTINE FINTGV(X,IPM,FVL,FINTG)
      IMPLICIT REAL*8 (A-H,O-Z)
C  INTEGRATES VALENCE TERM(S) - SEE NOTES
      DIMENSION FVL(30)
      COMMON/LAGUER/XL(8),WL(8),NTERML
      COMMON/GAUS16/XI(16),WI(16),NTERMS,XX(17)
      COMMON/GRPTHY/ FLAVOR
      COMMON/INPUT/INS,NVL,NMAX,DELTA,S0,IORD
      COMMON/ALQ2/T 
      COMMON/CONSTANTS/PI,PI2
      common/param/xc(60)
c      AL=ALPHA(T)/(2.*PI)
      xmu2=xc(1)**2*dexp(t)
      al=alpha_s(iord+1,xmu2,xc(1),neff)/(2.*pi)
      flavor=neff
      FX=GETFV(X,FVL)
      FINTG=T*AL*FX*(2.+8./3.*DLOG(1.-X))
      IF(IORD) 3,3,4
    4 C1=16./9.
      C2=2.
      C3=2./3.*FLAVOR
      C4=1.202056903D0
C      C5=-.7512853D0
      C5=-5./8.*C4
      A1=C2*2.*(67./9.-PI2/3.)-C3*20./9.
      C=C1*(3/8.-PI2/2.+C4-8.*C5)+C2*(17./12.+11.*PI2/9.
     2-C4+8.*C5)-C3*(1./6.+2.*PI2/9.)
      FINTG=FINTG+T*AL**2*FX*(C+A1*DLOG(1.-X))
    3 CONTINUE
      DO 1 I=1,NTERMS
      Z=0.5*(1.-X)*XI(I)+0.5*(1.+X)
      FZ=GETFV(X/Z,FVL)
      FINTG=FINTG+0.5*(1.-X)*T*AL*4./3.*((1.+Z*Z)*FZ-2.*FX)/(1.-Z)*WI(I)
      IF (IORD) 1,1,2
    2 ALZ=DLOG(Z)
      AL1=DLOG(1.-Z)
      Z2=1.+Z*Z
      ZM=1.-Z
      ZP=1.+Z
      AZ=C1*(-2.*Z2*ALZ*AL1-3.*ALZ)+C2*Z2*(ALZ*(ALZ+11./3.)+67./9.-PI2/
     23.)-C3*2./3.*Z2*(ALZ+5./3.)
      PA=4.*ZM+2.*ZP*ALZ+2.*Z2/ZP*S2(Z)
      BZ=C1*(-ALZ*(ALZ*ZP/2.+2.*Z)-5.*ZM+IPM*PA)+C2*(2.*ZP*ALZ
     2+40./3.*ZM-IPM*PA)-C3*4./3.*ZM
      FINTG=FINTG+0.5*(1.-X)*T*AL**2*((AZ*FZ-A1*FX)/ZM+BZ*FZ)*WI(I)
    1 CONTINUE
      RETURN
      END 



      FUNCTION S2(X)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION FN(51)
      DATA FN/1.644934076,1.588625448,1.545799712,1.507899041,
     X1.473125860,1.440633797,1.409928300,1.380685041,1.352675161,
     X1.325728728,1.299714723,1.274529160,1.250087584,1.226320101,
     X1.203167961,1.180581124,1.158516487,1.136936560,1.115808451,
     X1.095103088,1.074794600,1.054859830,1.035277934,1.016030062,
     X0.997099088,0.978469393,0.960126675,0.942057798,0.924250654,
     X0.906694053,0.889377624,0.872291733,0.855427404,0.838776261,
     X0.822330471,0.806082689,0.790026024,0.774153992,0.758460483,
     X0.742939737,0.727586308,0.712395042,0.697361058,0.682479725,
     X0.667746644,0.653157631,0.638708705,0.624396071,0.610216108,
     X0.596165361,0.582240526/
      COMMON/CONSTANTS/PI,PI2
C
C  THESE ARE THE VALUES OF F(X)=LI2(1-X) FOR X BETWEEN O AND .50 IN STEPS
C  OF O.O1 TAKEN FROM ABRAMOWITZ AND STEGUN, PG. 1005, TABLE 27.7.
C
C  S2(X)=INTEGRAL OF LN((1-Z)/Z)/Z FROM X/(1+X) TO 1/(1+X)
C
C  REWRITE S2(X) IN TERMS OF F(X/(1+X))
C
C  USE A LINEAR INTERPOLATION TO OBTAIN F(Z), Z=X/(1.+X)
C  IN THE REGION OF X BELOW 0.8
C
C  NEAR X=1 SMALL ERRORS IN THE INTERPOLATION ARE AMPLIFIED 
C  SINCE THE ANSWER IS A SMALL DIFFERENCE BETWEEN MUCH LARGER 
C  NUMBERS. USE A TAYLOR SERIES FOR F(X) EXPANDED ABOUT THE 
C  NEAREST X VALUE.
C
      Z=X/(1.+X)
      N=100.*Z+1
      Z1=(N-1)/100.
      Z2=N/100.
      IF(X.LT.0.8) THEN
         F=FN(N)*(Z-Z2)/(Z1-Z2)+FN(N+1)*(Z-Z1)/(Z2-Z1)
      ELSE
         DELT1=DABS(Z-Z1)
         DELT2=DABS(Z-Z2)
         ZT=Z1
         IF(DELT2.LT.DELT1) THEN
            ZT=Z2
            N=N+1
         ENDIF
         OMZT=1.-ZT
         F0=FN(N)
         F1=DLOG(ZT)/OMZT
         F2=(1./ZT+F1)/OMZT
         F3=(1./ZT**2+2.*F2)/OMZT
         F=F0+(Z-ZT)*F1+.5*(Z-ZT)**2*F2+(Z-ZT)**3/6.*F3
      ENDIF
      S2=PI**2/6.+.5*DLOG(X)**2-DLOG(1.+X)**2-2.*F
      RETURN
      END



      SUBROUTINE FINTGG(X,FG,FS,FINTG)
      IMPLICIT REAL*8 (A-H,O-Z)
C  INTEGRATES GLUON TERM - SEE NOTES
      DIMENSION FG(30),FS(30) 
      COMMON/GAUS16/XI(16),WI(16),NTERMS,XX(17)
      COMMON/GRPTHY/ FLAVOR
      COMMON/INPUT/INS,NVL,NMAX,DELTA,S0,IORD
      COMMON/ALQ2/T 
      COMMON/CONSTANTS/PI,PI2
      common/param/xc(60)
c      AL=ALPHA(T)/2./PI
      xmu2=xc(1)**2*dexp(t)
      al=alpha_s(iord+1,xmu2,xc(1),neff)/(2.*pi)
      flavor=neff
      FX=GETFV(X,FG)
      FINTG=T*AL*FX*((33.-2.*FLAVOR)/6.+6.*DLOG(1.-X))
      CF=4./3.
      IF(IORD) 2,2,3
    3 CONTINUE
      TR=FLAVOR/2.
      CA=3.
      Z3=1.202056903D0
      AL1=DLOG(1.-X)
      PGG1=-CF*TR-4./3.*CA*TR*(1.+5./3.*AL1)+CA*CA/3.*
     2(8.+9.*Z3+(67/3.-PI2)*AL1)
      FINTG=FINTG+T*AL*AL*FX*PGG1
    2 CONTINUE
      DO 1 I=1,NTERMS
      Z=0.5*(1.-X)*XI(I)+0.5*(1.+X)
      FGZ=GETFV(X/Z,FG)
      FSZ=GETFV(X/Z,FS)
      FINTG=FINTG+.5*(1.-X)*T*AL*WI(I)*(CF*FSZ*(Z*Z-2.*Z+2.)/Z+
     26.*(FGZ*Z-FX)/(1.-Z)+6.*FGZ*(1.-Z+Z*Z-Z**3)/Z)
      IF(IORD) 4,4,5
    5 CONTINUE
      AL1=DLOG(1.-Z)
      ALZ=DLOG(Z)
      PGG2=CF*TR*(4./3.*(1./Z-12.+6.*Z+5.*Z*Z)-2.*(3.+5.*Z)*ALZ
     2-2.*(1.+Z)*ALZ*ALZ)+4./3.*CA*TR*(-1./6.*(23./Z-29.+19.*Z-23.*Z*Z)
     3-(1.+Z)*ALZ)+CA*CA*(-1./18.*(25.+109.*Z)-1./3.*(25.-11.*Z+44.*Z*Z)
     4*ALZ+4.*(1.+Z)*ALZ*ALZ-PI2/3.*(1./Z-2.+Z-Z*Z)+(1./Z-2.+Z-Z*Z+
     51./(1.-Z))*ALZ*(ALZ-4.*AL1)-2.*(1./Z+2.+Z+Z*Z-1./(1.+Z))*S2(Z))
      PGG3=(-20./9.*CA*TR+CA*CA/3.*(67./3.-PI2))/(1.-Z)
      PGQ=4./3.*CF*TR*(2./3.*(5./Z-5.+4.*Z)+(2./Z-2.+Z)*AL1)
     2+.5*CF*CF*(5.+7.*Z-(4.+7.*Z)*ALZ+2.*(6./Z-6.+5*Z)*AL1 
     3+(2.-Z)*ALZ**2+2.*(2./Z-2+Z)*AL1**2)+CA*CF*(-1./9.*(9./Z+19.
     4+37.*Z+44.*Z*Z)+1./3.*(36.+15.*Z+8.*Z*Z)*ALZ-.5*(2./Z+6.
     5+3.*Z)*ALZ**2-1./3.*(22./Z-22.+17.*Z)*AL1+(2./Z-2.+Z)*(
     62.*ALZ*AL1+PI2/6.-AL1**2)+(2./Z+2.+Z)*S2(Z))
      PGQ=-PGQ
      FINTG=FINTG+.5*(1.-X)*T*AL*AL*WI(I)*(FSZ*PGQ+FGZ*PGG2+PGG3*(
     2FGZ-FX))
    4 CONTINUE
    1 CONTINUE
      RETURN
      END 



      SUBROUTINE FINTGS(X,FS,FG,FINTG)
      IMPLICIT REAL*8 (A-H,O-Z)
C  INTEGRATES SINGLET TERM - SEE NOTES
      DIMENSION FS(30),FG(30) 
      COMMON/GAUS16/XI(16),WI(16),NTERMS,XX(17)
      COMMON/GRPTHY/ FLAVOR
      COMMON/INPUT/INS,NVL,NMAX,DELTA,S0,IORD
      COMMON/ALQ2/T 
      COMMON/CONSTANTS/PI,PI2
      common/param/xc(60)
c      AL=ALPHA(T)/2./PI
      xmu2=xc(1)**2*dexp(t)
      al=alpha_s(iord+1,xmu2,xc(1),neff)/(2.*pi)
      flavor=neff
      FX=GETFV(X,FS)
      FINTG=FX*T*AL*(2.+8./3.*DLOG(1.-X))
      CF=4./3.
      IF(IORD) 2,2,3
    3 CONTINUE
      TR=FLAVOR/2.
      CA=3.
      Z3=1.202056903D0
      AL1=DLOG(1.-X)
      PFF1=-2./9.*CF*TR*(3./4.+PI2+10.*AL1)+CF*CF*(3.-
     24.*PI2+48.*Z3)/8.+CA*CF*(17./3.+44./9.*PI2-24.*Z3+
     3(536./9.-8.*PI2/3.)*AL1)/8.
      FINTG=FINTG+T*AL*AL*FX*PFF1
    2 CONTINUE
      DO 1 I=1,NTERMS
      Z=0.5*(1.-X)*XI(I)+0.5*(1.+X)
      FSZ=GETFV(X/Z,FS)
      FGZ=GETFV(X/Z,FG)
      FINTG=FINTG+.5*(1.-X)*T*AL*WI(I)*(CF*((1.+Z*Z)*FSZ-2.*FX)/(1.-Z)
     2+FLAVOR*((1.-Z)**2+Z*Z)*FGZ)
      IF(IORD) 4,4,5
    5 CONTINUE
      AL1=DLOG(1.-Z)
      ALZ=DLOG(Z)
      PFF2=2./9.*CF*TR*(20./Z-19.+65.*Z-56.*Z*Z+6.*(2.+8.*Z+4.*
     2Z*Z-1./(1.-Z))*ALZ-9.*(1.+Z)*ALZ**2)+CF*CF*(-1.+Z+(2.-3./
     3(1.-Z))*ALZ-.5*(1.+Z)*ALZ**2-2.*(1.+Z*Z)/(1.-Z)*ALZ*AL1+
     42.*(1.+Z*Z)/(1.+Z)*S2(Z))+CA*CF/8.*(4./9.*(17.-151*Z)
     5+4.*(1.+Z*Z)/(1.-Z)*ALZ*(11./3.+ALZ)+4.*(1.+Z)*PI2/3.-8.*
     6(1.+Z*Z)/(1.+Z)*S2(Z))
      PFF3=(-20./9.*CF*TR+(67./9.-PI2/3.)*CA*CF)/(1.-Z)
      PQG=-CF*TR*(-14.+29.*Z-20.*Z*Z-(3.-4.*Z+8.*Z*Z)*ALZ-
     2(1.-2.*Z+4.*Z*Z)*ALZ**2-8.*Z*(1.-Z)*AL1+2.*(1.-2.*Z+2.*Z*Z)*
     3(2.*ALZ*AL1+PI2/3.-AL1**2))-CA*TR*(-2./9.*(20./Z-18.+225.*Z-
     4218.*Z*Z)-2./3.*(3.+24.*Z+44.*Z*Z)*ALZ+(3.+6.*Z+2.*Z*Z)*
     5ALZ**2+8.*Z*(1.-Z)*AL1+(1.-2.*Z+2.*Z*Z)*(2*AL1**2-PI2/3.)
     6-2.*(1.+2.*Z+2.*Z*Z)*S2(Z))
      FINTG=FINTG+.5*(1.-X)*T*AL*AL*WI(I)*(FSZ*PFF2+PFF3*(FSZ-FX)
     2+PQG*FGZ)
    4 CONTINUE
    1 CONTINUE
      RETURN
      END 



      FUNCTION GETFV(X,FVL)
      IMPLICIT REAL*8 (A-H,O-Z)
C  LOGARITHMIC INTERPOLATOR - WATCH OUT FOR NEGATIVE
C  FUNCTIONS AND/OR X VALUES OUTSIDE THE RANGE 0 TO 1.
C  NOTE: DIMENSION OF FVL IS OVERWRITTEN BY VALUE USED
C  IN MAIN ROUTINE. 
      DIMENSION FVL(30),xx(3),fx(3)
      COMMON/GRID/NX,XGRID(30)
      IREP=-1
      DO 1 I=1,NX 
      IF(X.LT.XGRID(I)) GO TO 2
    1 CONTINUE
    2 I=I-1
      IF(I.EQ.0) THEN
         I=I+1
         IREP=IREP+1
      ELSE IF(I.GT.28) THEN
         I=28
      ENDIF
      xx(1)=xgrid(i)
      xx(2)=xgrid(i+1)
      xx(3)=xgrid(i+2)
      fx(1)=fvl(i)
      fx(2)=fvl(i+1)
      fx(3)=fvl(i+2)
      call polint(xx,fx,3,x,ans,dy)
      getfv=ans
      RETURN
      END 



      SUBROUTINE STARTV(FVLSTO,DELTA,FVL,FVL1,IPM)
      IMPLICIT REAL*8 (A-H,O-Z)
C  CALCULATES STARTING VALUES OF VALENCE DISTRIBUTIONS
C  AS NEEDED FOR THE PREDICTOR-CORRECTOR ALGORITHM
      DIMENSION FVLSTO(30),FVL(30),FVL1(30),TVLSTO(30)
      COMMON/GRID/NX,XGRID(30)
      DO 1 I=1,NX
    1 TVLSTO(I)=FVL(I)-DELTA*FVL1(I)
      DO 2 I=1,NX
      X=XGRID(I)
      CALL FINTGV(X,IPM,TVLSTO,FVLST1)
      FVLSTO(I)=FVL(I)-0.5*DELTA*(FVL1(I)+FVLST1) 
    2 CONTINUE
      RETURN
      END 



      SUBROUTINE STARTS(FS,FS1,FSSTO,FG,FG1,FGSTO,DELTA)
      IMPLICIT REAL*8 (A-H,O-Z)
C  CALULATES SINGLET AND GLUON STARTING VALUES AS NEEDED
C  FOR THE PREDICTOR-CORRECTOR ALGORITHM
      DIMENSION FS(30),FS1(30),FSSTO(30),TSSTO(30),
     $          FG(30),FG1(30),FGSTO(30),TGSTO(30)
      COMMON/GRID/NX,XGRID(30)
      DO 1 I=1,NX
      TSSTO(I)=FS(I)-DELTA*FS1(I)
      TGSTO(I)=FG(I)-DELTA*FG1(I)
    1 CONTINUE
      DO 2 I=1,NX
      X=XGRID(I)
      CALL FINTGS(X,TSSTO,TGSTO,FSSTO1)
      CALL FINTGG(X,TGSTO,TSSTO,FGSTO1) 
      FSSTO(I)=FS(I)-0.5*DELTA*(FS1(I)+FSSTO1)
      FGSTO(I)=FG(I)-0.5*DELTA*(FG1(I)+FGSTO1)
    2 CONTINUE
      RETURN
      END 



      SUBROUTINE WATE96
  !*******************************************************************
  !*****              *****
  !***** THE X(I) AND W(I) ARE THE DIRECT OUTPUT FROM A PROGRAM  *****
  !***** USING NAG ROUTINE D01BCF TO CALCULATE THE        *****
  !***** GAUSS-LEGENDRE WEIGHTS FOR 96 POINT INTEGRATION.        *****
  !***** THEY AGREE TO TYPICALLY 14 DECIMAL PLACES WITH THE      *****
  !***** TABLE IN ABRAMOWITZ & STEGUN, PAGE 919.         *****
  !*****              *****
  !***** ---->   PETER HARRIMAN, APRIL 3RD 1990.         *****
  !*****              *****
  !*******************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(48),W(48)
      COMMON/GAUS96/XI(96),WI(96),nterms,XX(97)
      NTERMS=96
      X( 1)=   0.01627674484960183561
      X( 2)=   0.04881298513604856015
      X( 3)=   0.08129749546442434360
      X( 4)=   0.11369585011066471632
      X( 5)=   0.14597371465489567682
      X( 6)=   0.17809688236761733390
      X( 7)=   0.21003131046056591064
      X( 8)=   0.24174315616383866556
      X( 9)=   0.27319881259104774468
      X(10)=   0.30436494435449495954
      X(11)=   0.33520852289262397655
      X(12)=   0.36569686147231213885
      X(13)=   0.39579764982890709712
      X(14)=   0.42547898840729897474
      X(15)=   0.45470942216774136446
      X(16)=   0.48345797392059470382
      X(17)=   0.51169417715466604391
      X(18)=   0.53938810832435567233
      X(19)=   0.56651041856139533470
      X(20)=   0.59303236477757022282
      X(21)=   0.61892584012546672523
      X(22)=   0.64416340378496526886
      X(23)=   0.66871831004391424358
      X(24)=   0.69256453664216964528
      X(25)=   0.71567681234896561582
      X(26)=   0.73803064374439816819
      X(27)=   0.75960234117664555964
      X(28)=   0.78036904386743123629
      X(29)=   0.80030874413913884180
      X(30)=   0.81940031073792957139
      X(31)=   0.83762351122818502758
      X(32)=   0.85495903343459936363
      X(33)=   0.87138850590929436968
      X(34)=   0.88689451740241818933
      X(35)=   0.90146063531585023110
      X(36)=   0.91507142312089592706
      X(37)=   0.92771245672230655266
      X(38)=   0.93937033975275308073
      X(39)=   0.95003271778443564022
      X(40)=   0.95968829144874048809
      X(41)=   0.96832682846326217918
      X(42)=   0.97593917458513455843
      X(43)=   0.98251726356301274934
      X(44)=   0.98805412632962202890
      X(45)=   0.99254390032376081654
      X(46)=   0.99598184298720747465
      X(47)=   0.99836437586317963722
      X(48)=   0.99968950388322870559
      W( 1)=   0.03255061449236316962
      W( 2)=   0.03251611871386883307
      W( 3)=   0.03244716371406427668
      W( 4)=   0.03234382256857594104
      W( 5)=   0.03220620479403026124
      W( 6)=   0.03203445623199267876
      W( 7)=   0.03182875889441101874
      W( 8)=   0.03158933077072719007
      W( 9)=   0.03131642559686137819
      W(10)=   0.03101033258631386231
      W(11)=   0.03067137612366917839
      W(12)=   0.03029991542082762553
      W(13)=   0.02989634413632842385
      W(14)=   0.02946108995816795100
      W(15)=   0.02899461415055528410
      W(16)=   0.02849741106508543861
      W(17)=   0.02797000761684838950
      W(18)=   0.02741296272602931385
      W(19)=   0.02682686672559184485
      W(20)=   0.02621234073567250055
      W(21)=   0.02557003600534944960
      W(22)=   0.02490063322248370695
      W(23)=   0.02420484179236479915
      W(24)=   0.02348339908592633665
      W(25)=   0.02273706965832950717
      W(26)=   0.02196664443874448477
      W(27)=   0.02117293989219144572
      W(28)=   0.02035679715433347898
      W(29)=   0.01951908114014518992
      W(30)=   0.01866067962741165898
      W(31)=   0.01778250231604547316
      W(32)=   0.01688547986424539715
      W(33)=   0.01597056290256253144
      W(34)=   0.01503872102699521608
      W(35)=   0.01409094177231515264
      W(36)=   0.01312822956696188190
      W(37)=   0.01215160467108866759
      W(38)=   0.01116210209983888144
      W(39)=   0.01016077053500880978
      W(40)=   0.00914867123078384552
      W(41)=   0.00812687692569928101
      W(42)=   0.00709647079115442616
      W(43)=   0.00605854550423662775
      W(44)=   0.00501420274292825661
      W(45)=   0.00396455433844564804
      W(46)=   0.00291073181793626202
      W(47)=   0.00185396078894924657
      W(48)=   0.00079679206555731759
      DO I=1,48
         XI(I)=-X(49-I)
         WI(I)=W(49-I)
         XI(I+48)=X(I)
         WI(I+48)=W(I)
      END DO
      DO I=1,96
         XX(I)=0.5*(XI(I)+1.)
      END DO
      XX(97)=1.0
      EXPON=1.0
      DO I=1,96
         YI=2.*(0.5*(1.+XI(I)))**EXPON-1.
         WI(I)=WI(I)/(1.+XI(I))*(1.+YI)*EXPON
         XI(I)=YI
         XX(I)=0.5*(1.+YI)
      END DO
      RETURN
      END SUBROUTINE WATE96




      SUBROUTINE WATE32
      IMPLICIT REAL*8 (A-H,O-Z)
C  32 POINT GAUSSIAN QUADRATURE ROUTINE 
      DIMENSION X(16),W(16)
      COMMON/GAUS32/XI(32),WI(32),NTERMS,XX(33)
      NTERMS=32
      X(1)=0.048307665687738316235
      X(2)=0.144471961582796493485
      X(3)=0.239287362252137074545
      X(4)=0.331868602282127649780
      X(5)=0.421351276130635345364
      X(6)=0.506899908932229390024
      X(7)=0.587715757240762329041
      X(8)=0.663044266930215200975
      X(9)=0.732182118740289680387
      X(10)=0.794483795967942406963
      X(11)=0.849367613732569970134
      X(12)=0.896321155766052123965
      X(13)=0.934906075937739689171
      X(14)=0.964762255587506430774
      X(15)=0.985611511545268335400
      X(16)=0.997263861849481563545
      W(1)=0.096540088514727800567
      W(2)=0.095638720079274859419
      W(3)=0.093844399080804565639
      W(4)=0.091173878695763884713
      W(5)=0.087652093004403811143
      W(6)=0.083311924226946755222
      W(7)=0.078193895787070306472
      W(8)=0.072345794108848506225
      W(9)=0.065822222776361846838
      W(10)=0.058684093478535547145
      W(11)=0.050998059262376176196
      W(12)=0.042835898022226680657
      W(13)=0.034273862913021433103
      W(14)=0.025392065309262059456
      W(15)=0.016274394730905670605
      W(16)=0.007018610009470096600
      DO 1 I=1,16
      XI(I)=-X(17-I)
      WI(I)=W(17-I) 
      XI(I+16)=X(I) 
      WI(I+16)=W(I) 
    1 CONTINUE
      DO 2 I=1,32
    2 XX(I)=0.5*(XI(I)+1.)
      XX(33)=1.0
      RETURN
      END 



      SUBROUTINE WATE16
      IMPLICIT REAL*8 (A-H,O-Z)
C  16 POINT GAUSSIAN QUADRATURE ROUTINE 
      DIMENSION X(8),W(8)
      COMMON/GAUS16/XI(16),WI(16),NTERMS,XX(17)
      NTERMS=16
      X(1)=0.095012509837637440185D0
      X(2)=0.281603550779258913230D0
      X(3)=0.458016777657227386342D0
      X(4)=0.617876244402643748447D0
      X(5)=0.755404408355003033895D0
      X(6)=0.865631202387831743880D0
      X(7)=0.944575023073232576078D0
      X(8)=0.989400934991649932596D0
      W(1)=0.189450610455068496285D0
      W(2)=0.182603415044923588867D0
      W(3)=0.169156519395002538189D0
      W(4)=0.149595988816576732081D0
      W(5)=0.124628971255533872052D0
      W(6)=0.095158511682492784810D0
      W(7)=0.062253523938647892863D0
      W(8)=0.027152459411754094852D0
      DO 1 I=1,8
      XI(I)=-X(9-I)
      WI(I)=W(9-I) 
      XI(I+8)=X(I) 
      WI(I+8)=W(I) 
    1 CONTINUE
      DO 2 I=1,16
    2 XX(I)=0.5*(XI(I)+1.)
      XX(17)=1.0
      RETURN
      END 



      SUBROUTINE WATE8
      IMPLICIT REAL*8 (A-H,O-Z)
C  8 POINT GAUSSIAN QUADRATURE ROUTINE
      COMMON/GAUSS8/XI(8),WI(8),NTERMS,XX(9)
      NTERMS=8
      XI(4)=-0.183434642495650D0
      XI(3)=-0.525532409916329D0
      XI(2)=-0.796666477413627D0
      XI(1)=-0.960289856497536D0
      WI(4)=0.362683783378362D0
      WI(3)=0.313706645877887D0
      WI(2)=0.222381034453374D0
      WI(1)=0.101228536290376D0
      DO 1 I=5,8
      XI(I)=-XI(9-I)
    1 WI(I)=WI(9-I) 
      DO 2 I=1,8
    2 XX(I)=0.5*(XI(I)+1.D0)
      XX(9)=1.
      RETURN
      END 



      FUNCTION XMNT32(F,N)
      IMPLICIT REAL*8 (A-H,O-Z)
C  CALULATES NTH MOMENT OF F USING 32 POINT GAUSSIAN
C  WARNING!! THIS IS NOT ENTIRELY ADEQUATE FOR N=2
C  AT HIGH Q**2 WHERE A LOW-X SPIKE DEVELOPS. A VARIABLE
C  TRANSFORMATION IS NEEDED IF ONE WANTS TO TEST THE
C  PROGRAM AT LOW X AND HIGH Q**2.
      DIMENSION F(1)
      COMMON/GAUS32/XI(32),WI(32),NTERMS,XX(33)
      XMNT32=0.0D0
      M=N-2
C
C  MODIFIED FOR N=1 TO OBTAIN GREATER ACCURACY
C
      IF(N.GT.1)THEN
      DO 1 I=1,NTERMS
      X=XX(I)
    1 XMNT32=XMNT32+0.5*X**M*F(I)*WI(I)
      ELSE
      DO 2 I=1,NTERMS
      X=XX(I)**2
      XMNT32=XMNT32+0.5*XX(I)**M*F(I)*WI(I)*2.
    2 CONTINUE
      ENDIF
      RETURN
      END 



      SUBROUTINE ADIMEN
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/ASFMOM/XMSTRS(25),XMSTRG(25),GPSI(25),GPSIA(25),GAPSI(25),
     $ GAA(25),GPLUS(25),GMINUS(25),ALPHAN(25),BETAN(25),EPSILN(25)
      COMMON/GRPTHY/ FLAVOR
 2000 FORMAT(10X,'CONSTANTS OF ASYMPTOTIC FREEDOM FOR',F3.0,' FLAVORS',
     $/10X,'NOTATION FOLLOWS A. J. BURAS, NUCL. PHYS. B125,125(1977)')
C    $/10X,'NOTATION FOLLOWS GROSS AND WILCZEK')
C    $/10X,'NOTATION FOLLOWS FLORATOS, ROSS AND SACHRAJDA') 
      WRITE(6,2001) 
 2001 FORMAT(/3X,'N',5X,'G(PSI,PSI)',2X,'G(PSI,A)',4X,'G(A,PSI)',4X,
     $'G(A,A)',6X,'G(+)',8X,'G(-)',8X,'ALPHA',7X,'BETA',8X,'EPSILON'/)
      B=(33.-2.*FLAVOR)/3.
      S=0.0
      DO 1 N=2,25
      S=S+1./N
      GPSI(N)=8./(6.*B)*(1.-2./N/(N+1.)+4.*S)
      IF(IPOLZN) 60,60,70
   60 CONTINUE
      GPSIA(N)=FLAVOR/(2.*B)*(8./(N+2.)+16./N/(N+1.)/(N+2.))
C     GPSIA(N)=-0.5*GPSIA(N)
      GAPSI(N)=8./(6.*B)*(1./(N+1.)+2./N/(N-1.))
C     GAPSI(N)=-2.*GAPSI(N)
      GAA(N)=3./B*(1./3.-4./N/(N-1.)-4./(N+1.)/(N+2.)+4.*S)+2.*FLAVOR/3.
     $/B
      GO TO 80
   70 CONTINUE
      GPSIA(N)=FLAVOR/(2.*B)* 8.*(N-1.)/N/(N+1.)
      GAPSI(N)=8./(6.*B)*(N+2.)/N/(N+1.)
      GAA(N)=3./B*(1./3.-8./N/(N+1.)+4.*S)+2.*FLAVOR/3./B
   80 CONTINUE
C        1         2         3         4         5         6         7 *
      GPLUS(N)=0.5*(GPSI(N)+GAA(N)+DSQRT((GPSI(N)-GAA(N))**2+4.*GPSIA(N)
     $ *GAPSI(N)))
      GMINUS(N)=0.5*(GPSI(N)+GAA(N)-DSQRT((GPSI(N)-GAA(N))**2+
     $4.*GPSIA(N)*GAPSI(N)))
      ALPHAN(N)=GPSIA(N)*GAPSI(N)/(GPSIA(N)*GAPSI(N)+(GMINUS(N)-GPSI(N))
     $**2)
      BETAN(N)=0.5*ALPHAN(N)*(GPSI(N)-GMINUS(N))/GAPSI(N)
      EPSILN(N)=ALPHAN(N)*(1.-ALPHAN(N))/BETAN(N) 
      WRITE(6,1000) N,GPSI(N),GPSIA(N),GAPSI(N),GAA(N),GPLUS(N),
     $GMINUS(N),ALPHAN(N),BETAN(N),EPSILN(N)
 1000 FORMAT(1X,I3,3X,9(F10.7,2X))
    1 CONTINUE
      RETURN
      END 



      SUBROUTINE WATE8L
      IMPLICIT REAL*8 (A-H,O-Z)
C  8 POINT LAGUERRE INTGRATOR (FOR EXPONENTIALS OR OTHER
C  STEEPLY FALLING FUNCTIONS
      COMMON/LAGUER/XL(8),WL(8),NTERML
      XL(1)=0.1702796
      XL(2)=0.9037018
      XL(3)=2.2510866
      XL(4)=4.2667002
      XL(5)=7.0459054
      XL(6)=10.758516
      XL(7)=15.740679
      XL(8)=22.863132
      WL(1)=0.4377234
      WL(2)=1.0338693
      WL(3)=1.6697098
      WL(4)=2.3769247
      WL(5)=3.2085409
      WL(6)=4.2685755
      WL(7)=5.8180834
      WL(8)=8.9062262
      NTERML=8
      RETURN
      END 
