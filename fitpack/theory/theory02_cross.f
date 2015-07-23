      FUNCTION THEORY(MODE,NPT,V,XC)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*4 PAR(60),V4(4),ANSWER
      real*8 nuke_cteq
      DIMENSION V(4),XC(60),GF(8,40,1080)
     2,VP(4),FTEMP(2),NFLAG(30),INORM(30),itype(30),ITGT(30),icorr(30)
     3,itmc(30),iht(30)
      COMMON/GINT/GF
      COMMON/INPUT/INS,NVL,NMAX,DELTA,S0,IORD
      COMMON/GFUNC/CALC(8,40,30)
      COMMON/Q2STUF/ Q02,Q2MAX
      COMMON/GAUSS4/XI(4),WI(4),NTERMS,XX(5)
      COMMON/GRPTHY/FLAVOR
      COMMON/FLAGS/ITYPE,NFLAG,INORM,ITGT,icorr,itmc,iht
      COMMON/MINIMB/CURPAR(60),UNCRT(60),ERMIN(40),PWATE(60),IFREP(40)
      COMMON/CONSTANTS/PI,PI2
      common/target/az,an
      common/e866xf/ixfx
      common/jetpdfs/ndrv,nfl,als

*     Higher-twist parameter indexes
      integer nht(4)
      common/ht/nht


C  THIS ROUTINE CONTROLS THE CALCULATION OF THE VARIOUS
C  OBSERVABLES WHICH ARE TO BE FITTED. IT MAKES USE OF A NUMBER OF
C  ASSOCIATED ROUTINES AND CAN EASILY BE EXPANDED FOR THE CALCULATION 
C  OF ADDITIONAL QUANTITIES.
c
c  [A.Accardi] Updated 5/19/08 to include TMC in collinear factorization 
c  and approximate Georgi-Politzer TMC in a unified way
c
c  [A.Accardi] 9/9/08 inclusion of Higher-Twist corrections
c
c  Approximate Georgi-Politzer target mass corrections added 5/5/08
c
c  modified 12/05 to calculate nutev cross section data
c  corrected DIS bug for calculating f_2^d/f_2^p
c  corrected bug in gluon coefficient for DIS with nf=5
c  generalized 'theory2' and 'fetch' to make the routines more streamlined 
c  and versatile
C
C  THEORY2 CALCULATES THE DIS OBSERVABLES
C  DYANTH CALCULATES THE DRELL-YAN AND J/PSI OBSERVABLES
C
C  DECIDE WHETHER OR NOT TO CALL INTQCD 
C
C  MODE=1....NO DERIVATIVES NEEDED
C      =0....CALCULATE DERIVATIVES
C
C  NPT=1....FILL GF 
C     >1....PROCEED 
C
C  MODIFIED 11/87 SO THAT EVOLUTION IS NOT PERFORMED WHEN 
C  NORMALIZATION PARAMETERS ARE VARIED...SAVES TIME.
C
*  NOTES:
*
*  (09 Sep 08) HT only for F2
*              GP only for F2, F3
*              CF, NV only for F1, F2, F3
*

      IF(NPT-1)1,1,2
    1 IF(MODE-1) 3,4,3
    4 NDRV=0
      NOLD=1
      GO TO 7
    3 CONTINUE
      NDRV=NDRV+1
      IF(IFREP(NDRV)-35)7,7,5 
    7 CONTINUE
      S0=DLOG(Q02/XC(1)**2)
      SMAX=DLOG(Q2MAX/XC(1)**2)
      SMAX=DLOG(SMAX/S0)
      NMAX=SMAX/DELTA+2
C
C  FILL ARRAY GF
C
      CALL INTQCD(XC)
      DO 6 J=1,40
      DO 6 K=1,30
      KK=K+NDRV*30
      DO 20 IL=1,8
      GF(IL,J,KK)=CALC(IL,J,K)
   20 CONTINUE
    6 CONTINUE
      GO TO 5
    2 IF(MODE-1) 8,9,8
    9 NDRV=0
      GO TO 5
    8 IF(NOLD-NPT) 10,11,10
   11 NDRV=NDRV+1
      GO TO 5
   10 NOLD=NPT
      NDRV=1
    5 CONTINUE
      index=V(4)
      it=itype(index)
      nfl=nflag(index)
      itg=itgt(index)
      itm=itmc(index)
      ihtfl=iht(index)
*    ... old (pre-080519) vp assignments
      vp(1)=v(1)
      vp(2)=v(2)
      vp(3)=v(3)
      vp(4)=v(4)

      IF(IT.LT.100)then

*       *** [AA 080519] inclusion of DIS TMC in col.fact. Modified the 
*       ... meaning of VP(3), instead of VP(3)=V(3)=W^2 as previously
*       ... done. Not a problem yet for "THEORY2" subroutine because
*       ... W^2 is not used in the str.fn. computations. 
*       ... NOTE: this part was moved after the IF(IT.LT.100) 
*       ... because it is specific to DIS. 
         x=v(1)
         q2=v(2)
         amu=.8836/q2
         r=sqrt(1.+4.*x**2*amu)
         xtmc=2.*x/(1.+r)
         if(itm.eq.0)then
*          *** No TMC: 
*          ... vp(1)=xB  vp(3)=xB
            vp(1) = v(1)
            vp(3) = v(1)
            f2cor=1.
            f3cor=1.
         else if(itm.eq.1)then
*          *** approximate GP TMC:
*          ... VP(1) = xi_Nachtmann, and VP(3) = xi_Nachtmann
            vp(1)=xtmc
            vp(3)=xtmc
*          ... f1cor not yet implemented !!!!
            f1cor = -1d0
            f2cor=(x/xtmc)**2/r**3*(1.+6.*amu*x*xtmc/r*(1.-xtmc)**2)
            f3cor=x/xtmc/r**2*(1.-amu*x*xtmc/r*(1.-xtmc)*dlog(xtmc))
         else if (itm.eq.2) then
*          *** TMC in col.fact.
*          ... VP(1) = xi_Nachtmann, and VP(3) = xB
            vp(1)=xtmc
            vp(3)=x
            f1cor=1.
            f2cor=x/(xtmc*r**2)
            f3cor=x/(xtmc*r)
         else 
*          *** naive TMC in col.fact.
*          ... VP(1) = xi_Nachtmann, and VP(3) = xi_Nachtmann
            vp(1)=xtmc
            vp(3)=xtmc
            f1cor=1.
            f2cor=x/(xtmc*r**2)
            f3cor=x/(xtmc*r)
         endif

*       ... Higher-twist term
*       ... (09 Sep 08) for the moment only for F2
*       ... ihtfl  =  0      no HT
*       ... ihtfl  =  1-10   F2(TMC) = F2(TMC) + ht_f2
*       ... ihtfl  = 11-20   F2(TMC) = F2(TMC) * (1+ht_f2)
*                            with same ht_f2 as ihtfl-10
         if (ihtfl.le.10) then
            xht = x
         else 
            xht = xtmc
         end if

         if (ihtfl.eq.0) then
            ht_f2 = 0d0
         else if (ihtfl.eq.1) then  ! a * xB^b * (1+c*xB+d*xB^2) / Q^2
            ht_f2 = curpar(nht(1)) 
     &           * xht**curpar(nht(2)) 
     &           * (1d0 + curpar(nht(3))*xht) 
     &           * (1-xht)**curpar(nht(4))
     &           / Q2
         end if
         
*       ... Computes structure functions
         IF(nfl.eq.9)then
C
C  N/P RATIO CALCULATED HERE for charged lepton DIS
C
C  first get f_2 (n+p)/2
c
            az=.5
            an=.5
            CALL THEORY2(1,2,2,NDRV,VP,XC,F2NP)
            f2np=f2np*f2cor * (1d0 + ht_f2)
C
C  CALCULATE F2 PROTON HERE
C
            az=1.
            an=0.
            CALL THEORY2(1,1,2,NDRV,VP,XC,F2P)
            f2p=f2p*f2cor * (1d0 + ht_f2)
*               print*, '*',x,q2,ht_f2,index,ihtfl,(curpar(nht(i)),i=1,3)
c
c  modified to F2D/F2P  10/21/05
c
            theory=f2np/f2p
            RETURN
         endif
C
C  DIS STRUCTURE FUNCTIONS CALCULATED HERE
C
c  the logic here must be updated as new data sets are added 
c
c  note: s=sbar, c=cbar, and b=bar are assumed for this version
c
c        higher-twist only for F2 (09 sep 2008)
c
         if(nfl.eq.1)then
            call theory2(1,1,2,ndrv,vp,xc,ans)
            ans=ans*f2cor * (1d0 + ht_f2)
*               print*, x,q2,ht_f2, (curpar(nht(i)),i=1,4)
         else if(nfl.eq.2)then
            az=.5
            an=.5
            call theory2(1,2,2,ndrv,vp,xc,ans)
            ans=ans*f2cor * (1d0 + ht_f2)
         else if(nfl.eq.5)then
            az=.5
            an=.5
            call theory2(2,2,2,ndrv,vp,xc,ans1)
            call theory2(3,2,2,ndrv,vp,xc,ans2)
            ans1=ans1*f2cor * (1d0 + ht_f2)
            ans2=ans2*f2cor * (1d0 + ht_f2)
            if(itg.eq.3)then
*             ... nuclear corrections
               ans1=ans1*nuke_cteq(v(1),v(2),2,0,2,3,2)
               ans2=ans2*nuke_cteq(v(1),v(2),2,0,-2,3,2)
            endif
            ans=(ans1+ans2)/2.
         else if(nfl.eq.8)then
*          ... NOTE: no HT for F3 yet
            az=.5
            an=.5
            call theory2(2,2,3,ndrv,vp,xc,ans1)
            call theory2(3,2,3,ndrv,vp,xc,ans2)
            ans1=ans1*f3cor
            ans2=ans2*f3cor
            if(itg.eq.3)then
               ans1=ans1*nuke_cteq(v(1),v(2),3,0,2,3,2)
               ans2=ans2*nuke_cteq(v(1),v(2),3,0,-2,3,2)
            endif
            ans=(ans1+ans2)/2.

c
c  no TMC nor HT for nfl=10 or 11 yet
c
         else if(nfl.eq.10)then
            x=v(1)
            q2=v(2)
            y=v(3)
            if(itg.eq.3)then
               az=0.5
               an=0.5
               call theory2(2,2,1,ndrv,v,xc,ans1)
               call theory2(2,2,2,ndrv,v,xc,ans2)
               call theory2(2,2,3,ndrv,v,xc,ans3)
               ans1=ans1*nuke_cteq(x,q2,1,0,2,3,2)
               ans2=ans2*nuke_cteq(x,q2,2,0,2,3,2)
               ans3=ans3*nuke_cteq(x,q2,3,0,2,3,2)
            else
               az=0.465
               an=0.535
               call theory2(2,2,1,ndrv,v,xc,ans1)
               call theory2(2,2,2,ndrv,v,xc,ans2)
               call theory2(2,2,3,ndrv,v,xc,ans3)
            endif
            ans=(1.-y-(.938*x*y)**2/q2)*ans2+y**2*x*ans1
     2          +y*(1-y/2.)*ans3
            ans=ans*1.5816/(1.+q2/80.22**2)**2
         else if(nfl.eq.11)then
            x=v(1)
            q2=v(2)
            y=v(3)
            if(itg.eq.3)then
               az=0.5
               an=0.5
               call theory2(3,2,1,ndrv,v,xc,ans1)
               call theory2(3,2,2,ndrv,v,xc,ans2)
               call theory2(3,2,3,ndrv,v,xc,ans3nubar)
               call theory2(2,2,3,ndrv,v,xc,ans3nu)
               ans1=ans1*nuke_cteq(x,q2,1,0,-2,3,2)
               ans2=ans2*nuke_cteq(x,q2,2,0,-2,3,2)
               ans3=(ans3nu+ans3nubar)*nuke_cteq(x,q2,3,0,-2,3,2)
     2                         -ans3nu*nuke_cteq(x,q2,3,0,2,3,2)
            else
               az=0.465
               an=0.535
               call theory2(3,2,1,ndrv,v,xc,ans1)
               call theory2(3,2,2,ndrv,v,xc,ans2)
               call theory2(3,2,3,ndrv,v,xc,ans3)
            endif
            ans=(1.-y-(.938*x*y)**2/q2)*ans2+y**2*x*ans1
     2          -y*(1-y/2.)*ans3
            ans=ans*1.5816/(1.+q2/80.22**2)**2
         endif

      else if(it.gt.100.and.it.le.200)then
C  AS CURRENTLY SET UP ITYPE < 100 CORRESPONDS TO DEEP
C  INELASTIC SCATTERING WHILE DRELL-YAN AND W PRODUCTION HAVE ITYPE = 100
C  OR GREATER. THIS CONVENTION MUST BE REMEMBERED IF ADDITIONS OR
C  ALTERATIONS ARE MADE.
         az=1.
         an=0. 
         ixfx=0
         IF(nfl.eq.1.or.nfl.eq.2)THEN
            if(it.eq.106)then
               az=0.456
               an=0.544
            else if(it.eq.108)then
               az=1.
               an=0.
               ixfx=1
            else if(it.eq.110)then
               az=0.5
               an=0.5 
               ixfx=1
            endif
            CALL DYANTH(NDRV,V,XC,THEORY)
            IF(IT.EQ.106.or.it.eq.108.or.it.eq.110) THEN
               THEORY=THEORY*(V(2)/SQRT(V(1)))**3
               in=inorm(index)
               theory=theory/xc(in)
            ENDIF
         ELSE IF(IT.EQ.125)THEN         
            CALL DYANTH(NDRV,V,XC,THEORYP)
            az=0.
            an=1.
            CALL DYANTH(NDRV,V,XC,THEORYN)
            THEORY=(THEORYP-THEORYN)/(THEORYP+THEORYN)
         ELSE IF(IT.GE.126.and.IT.LE.129)THEN
            S=V(1)
            Y=V(3)
            RS=SQRT(S)
            ALS=DLOG(DLOG(80.**2/XC(1)**2)/S0)
            CALL WASYM(RS,Y,ALS,NDRV,THEORY,it)
         ENDIF
         return

      else if(it.gt.200.and.it.le.300)then
* JET observables
         ymin=v(1)
         ymax=v(2)
         pt=v(3)
c
c  hardwire these for now - may pass from input file later on
c
         rs=1800.
         smucoef=.5
c
         q2=(smucoef*pt)**2
         s0=dlog(q02/xc(1)**2)
         als=dlog(dlog(q2/xc(1)**2)/s0)
         call jetdet(rs,pt,smucoef,ans,ymin,ymax,iord)         
      else if(it.gt.300.and.it.lt.400)then
c
c  gamma + jet
c
         rs=v(1)
         pt=v(2)
         nregion=v(3)
c
c  hardwire scale for now
c
         scale=2.0
         q2=(scale*pt)**2
         s0=dlog(q02/xc(1)**2)
         als=dlog(dlog(q2/xc(1)**2)/s0)
         s=rs**2
         call d0gamjet(s,pt,nregion,ndrv,als,ans,scale)
         if(ans.eq.0.d0)then
            print*,pt,nregion
         endif
      endif   
      IN=INORM(index)
      THEORY=ANS/XC(IN)
      RETURN
      END 



      SUBROUTINE THEORY2(ibeam,itarg,istruc,NDRV,V,XC,THEORY)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION V(4),XC(60),FTEMP(2)
      COMMON/GAUS16/XI(16),WI(16),NTERMS,XX(17)
c      COMMON/GAUS32/XI(32),WI(32),NTERMS,XX(33)
      COMMON/GRPTHY/FLAVOR
      COMMON/INPUT/INS,NVL,NMAX,DELTA,S0,IORD
      COMMON/Q2STUF/Q02,Q2MAX 
      common/flags/itype(30),nflag(30),inorm(30),itgt(30),icorr(30)
     2     ,itmc(30),iht(30)
      COMMON/CONSTANTS/PI,PI2
C  CALCULATES DEEP INELASTIC OBSERVABLES
C  FETCH IS USED TO UNFOLD THE PROPER QUARK DISTRIBUTIONS
C  FROM THOSE WHICH ARE EVOLVED BY INTQCD.
C  FOR THE HIGHER ORDER CALCULATION THE CONVOLUTION WITH
C  THE COEFFICIENT FUNCTIONS IS HANDLED HERE. THE CONVENTIONS
C  FOR THE GLUON COEFFICIENT FUNCTION ARE EXPLAINED BELOW.
c
c  [A.Accardi] Updated 5/19/08 to include TMC in collinear factorization 
c  and approximate Georgi-Politzer TMC in a unified way.
c
c  updated 12/05 to make the package more streamlined and flexible and 
c  to allow the calculation of the nutev cross section data
c

*    *** input vars, and derived quantities
*    ... NOTE: x here can be x=xB or x=xi_Nachtmann depending on 
*    ... the calling routine
      X=V(1)
      Q2=V(2)
      S0=DLOG(Q02/XC(1)**2)
      S=DLOG(DLOG(Q2/XC(1)**2)/S0)

*    *** tree level str.fn. 
      CALL FETCH(ibeam,itarg,istruc,NDRV,X,S,FTEMP)
      THEORY=FTEMP(1)
      IF(IORD.EQ.0)then
         if(istruc.eq.1)theory=theory/2./x
         return
      endif

*    *** 1 loop contribution
c
c  flags for F1 or xF3
c
      if1=0
      IF3=0
      if(istruc.eq.1)then
        if1=1
      else if(istruc.eq.3)then
        if3=1
      endif
      CF=4./3.
      AL1=DLOG(1.-X)
      al=alpha_s(iord+1,q2,xc(1),neff)/(4.*pi)
      flavor=neff
      FAC=FLAVOR
c
c  factor for sum over quark squared charges
c  flavor starts at 4 at Q0 and switches to 5 at the b threshold
c  flavor = 3 or 6 are never encountered in the current kinematic range
c
      ymax = v(1)/v(3)
      ymin = v(1)
      IF(ibeam.eq.1)fac=10./9.+(flavor-4.)/9. 
      FX=THEORY
      THEORY=FX+FX*AL*CF*(-9.-2.*PI2/3.+AL1*(-3.+2.*AL1))
      DO 23 I=1,NTERMS
         Y=0.5*(ymax-ymin)*XI(I)+0.5*(ymax+ymin)
         XY=X/Y
         AL1=DLOG(1.-Y)
         CALL FETCH(ibeam,itarg,istruc,NDRV,XY,S,FTEMP)
         C22=CF*(6.+4.*Y-2.*(1.+Y*Y)/(1.-Y)*DLOG(Y)-2.*(1.+Y)*AL1
     2        -IF3*2.*(1.+Y)-if1*4.*y)
         C23=CF*(-3.+4.*AL1)/(1.-Y)
         CG2=2.*FAC*(-1.+8.*Y*(1.-Y)+(1.-2.*Y+2.*Y*Y)*DLOG(1./Y-1.)
     2        -if1*4.*Y*(1.-Y))
C     THE ABOVE GLUON COEFFICIENT FUNCTION CORRESPONDS TO THE
C     CONVENTIONS OF FLORATOS,HERROD,WADA,ETC. THE FOLLOWING
C     EXPRESSION CORRESPONDS TO THE CONVENTION OF ALTARELLI,
C     ELLIS,AND MARTINELLI.
C     CG2=2.*FAC*(6.*Y*(1.-Y)+(1.-2.*Y+2.*Y*Y)*DLOG(1./Y-1.))
         THEORY = THEORY + .5*(ymax-ymin)*WI(I)*AL*(C22*FTEMP(1)
     &        +C23*(FTEMP(1)-FX))
         if(istruc.lt.3)THEORY=THEORY
     &        +.5*(ymax-ymin)*WI(I)*AL*CG2*FTEMP(2)
   23 CONTINUE
      if(istruc.eq.1)theory=theory/2./x
      RETURN
      END 


      SUBROUTINE FETCH(ibeam,itarg,istruc,NDRV,X,S,FTEMP)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION FTEMP(2)
      COMMON/INPUT/INS,NVL,NMAX,DELTA,S0,IORD
      COMMON/GRPTHY/FLAVOR
      common/target/az,an
C
C  UNFOLDS QUARK AND GLUON DISTRIBUTIONS FOR USE IN THEORY2.
C
      call pdf(ndrv,x,s,u,d,ub,db,sb,cb,bb,glue)
c
c  s=sbar, c=cbar, b=bbar
c
      sq=sb
      cq=cb
      bq=bb
c
c  correct for neutron/proton ratio unless it is a proton target
c
      if(itarg.eq.2)then
         tu=u
         tub=ub
         td=d
         tdb=db
         u=az*tu+an*td
         ub=az*tub+an*tdb
         d=az*td+an*tu
         db=az*tdb+an*tub
      endif
c
c  change sign of antiquarks if for xf3
c
      iqbar=1
      if(istruc.eq.3) iqbar=-1
c
c  set up q and qbar coefficients
c
      if(ibeam.eq.1)then
         fu=4./9.
         fub=4./9.
         fd=1./9.
         fdb=1./9.
      else if(ibeam.eq.2)then
         fu=0.
         fub=2.*iqbar
         fd=2.
         fdb=0.
      else if(ibeam.eq.3)then
         fu=2.
         fub=0.
         fd=0.
         fdb=2.*iqbar
      endif
      ftemp(1)=fu*(u+cq)+fub*(ub+cb)+fd*(d+sq)+fdb*(db+sb)
c
c  add initial b quark contribution if for charged lepton DIS
c
      if(ibeam.eq.1) ftemp(1)=ftemp(1)+fd*bq+fdb*bb
      ftemp(2)=glue
      return
      end


      SUBROUTINE DYANTH(NDRV,V,XC,THEORY)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION V(4),XC(60)
      COMMON/Q2STUF/Q02,Q2MAX 
      COMMON/GAUS16/XI(16),WI(16),NTERMS,XX(17)
      COMMON/GRPTHY/FLAVOR
      COMMON/INPUT/INS,NVL,NMAX,DELTA,S0,IORD
      common/flags/itype(30),nflag(30),inorm(30),itgt(30),icorr(30)
     2     ,itmc(30),iht(30)
      index=V(4)
      IFL=nflag(index)
      S=V(1)
      Q2=V(2)**2
      IF(IFL.EQ.5) Q2=Q2/2.
      ALS0=DLOG(Q02/XC(1)**2) 
      ALS=DLOG(DLOG(Q2/XC(1)**2)/ALS0)
      Y=V(3)
      CALL HODY(IFL,NDRV,S,Q2,Y,ALS,THEORY)
      RETURN
      END 


      SUBROUTINE HODY(IFL,NDRV,S,QS,Y,ALS,THEORY)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION Q10(-5:5),Q1(-5:5,32),Q20(-5:5),Q2(-5:5,32),QT(-5:5)
     2,HQ1Q(32),HQQ2(32),HQG2(32),HG1Q(32)
      COMMON/GAUS32/XI(32),WI(32),NTERMS,XX(33)
      COMMON/CONSTANTS/PI,PI2
      COMMON/PARAM/XC(60)
      COMMON/DYTEST/XSECLO,GLUCOR
      COMMON/INPUT/INS,NVL,NMAX,DELTA,S0,IORD
      common/e866xf/ixfx
C
C  HIGHER ORDER DRELL-YAN SUBROUTINE ADDED APRIL 29, 1991
C  CALCULATES S**1.5 DSIGMA/DQ/DY
C
c  modifications added to allow calculation of s**1.5 dsigma/dq dxf
c  for E866  3/5/03
c
      TAU=QS/S
      if(ixfx.eq.1)then
         xf=y
         rt=sqrt(xf**2+4.*tau)
         x10=(xf+rt)/2.
         x20=(-xf+rt)/2.
         xffac=1./(x10+x20)
      else
         EX=EXP(Y)
         X10=SQRT(TAU)*EX
         X20=X10/EX**2
         xffac=1.
      endif
      al=alpha_s(iord+1,qs,xc(1),neff)
      flavor=neff
      ACF=AL*4./3./(2.*PI)
      AT=3.*ACF/4.
      FAC=8.*PI/137.**2/9./SQRT(TAU)*389.E03
      CALL HOQUARK(1,IFL,X10,ALS,NDRV,Q10)
      CALL HOQUARK(2,IFL,X20,ALS,NDRV,Q20)
      HQQ=(4.*Q10(1)*Q20(-1)+Q10(2)*Q20(-2)+Q10(3)*Q20(3)
     2    +4.*Q10(-1)*Q20(1)+Q10(-2)*Q20(2)+Q10(-3)*Q20(3)
     3    +8.*Q10(4)*Q20(4)+2.*q10(5)*q20(5))/9.
      THEORY=FAC*HQQ*(1.+
     2       ACF*(-8.+PI2+DLOG(X10*X20/(1.-X10)/(1.-X20))**2))
      theory=theory*xffac
      XSECLO=FAC*HQQ*xffac
      tmp0=(qs/s)**1.5*xseclo
      if(iord.eq.0)then
         theory=xseclo
         return
      endif
      QQ1=0.
      QG1=0.
      DO 100 I=1,NTERMS
      X2=.5*(1.-X20)*XI(I)+.5*(1.+X20)
      Z=X20/X2
      CALL HOQUARK(2,IFL,X2,ALS,NDRV,QT)
      DO 101 IQ=-5,5
  101 Q2(IQ,I)=QT(IQ)
      HQQ2(I)=(4.*Q10(1)*Q2(-1,I)+Q10(2)*Q2(-2,I)+Q10(3)*Q2(-3,I)
     2    +4.*Q10(-1)*Q2(1,I)+Q10(-2)*Q2(2,I)+Q10(-3)*Q2(3,I)
     3    +8.*Q10(4)*Q2(4,I)+2.*q10(5)*q2(5,i))/9.
C        1         2         3         4         5         6         7 *
      if(ixfx.eq.1)then
         TEMP=((1.+Z*Z)*DLOG((x10+x20)*(1.-X10)/X10/x20/(X2+X10))
     2        *HQQ2(I)-2.*DLOG((1.-X10)/X10/X20)
     3         *HQQ)/(X2-X20)+2.*DLOG(X2-X20)/(X2-X20)*(Z*HQQ2(I)-HQQ)
     4         +HQQ2(I)*(X2-X20)/X2**2*(1.+DLOG(X2-X20))
      else
         TEMP=((1.+Z*Z)*DLOG(2.*(1.-X10)/X10/(X2+X20))*HQQ2(I)
     2        -2.*DLOG((1.-X10)/X10/X20)
     3         *HQQ)/(X2-X20)+2.*DLOG(X2-X20)/(X2-X20)*(Z*HQQ2(I)-HQQ)
     4         +HQQ2(I)*(X2-X20)/X2**2*(1.+DLOG(X2-X20))
      endif
      QQ1=QQ1+.5*(1.-X20)*WI(I)*FAC*ACF*TEMP*xffac
      HQG2(I)=(4.*Q10(1)+Q10(2)+Q10(3)+4.*Q10(-1)+Q10(-2)+Q10(-3)
     2        +8.*Q10(4)+2.*q10(5))*Q2(0,I)/9.
      if(ixfx.eq.1)then
         TEMP=((X20**2+(X2-X20)**2)/(2.*X2**3)*
     2         DLOG((X2-X20)*(1.-X10)*(x10+x20)/X10/x20/(X2+X10))
     3         +X20*(X2-X20)/X2**3)*HQG2(I)
      else
         TEMP=((X20**2+(X2-X20)**2)/(2.*X2**3)*
     2         DLOG(2.*(X2-X20)*(1.-X10)/X10/(X2+X20))
     3         +X20*(X2-X20)/X2**3)*HQG2(I)
      endif
      QG1=QG1+.5*(1.-X20)*WI(I)*FAC*AT*TEMP*xffac
  100 CONTINUE
      QQ2=0.
      QG2=0.
      DO 200 I=1,NTERMS
      X1=.5*(1.-X10)*XI(I)+.5*(1.+X10)
      Z=X10/X1
      CALL HOQUARK(1,IFL,X1,ALS,NDRV,QT)
      DO 201 IQ=-5,5
  201 Q1(IQ,I)=QT(IQ)
      HQ1Q(I)=(4.*Q1(1,I)*Q20(-1)+Q1(2,I)*Q20(-2)+Q1(3,I)*Q20(-3)
     2    +4.*Q1(-1,I)*Q20(1)+Q1(-2,I)*Q20(2)+Q1(-3,I)*Q20(3)
     3    +8.*Q1(4,I)*Q20(4)+2.*q1(5,i)*q20(5))/9.
      if(ixfx.eq.1)then
         TEMP=((1.+Z*Z)*DLOG((x10+x20)*(1.-X20)/x10/X20/(X1+X20))
     2        *HQ1Q(I)-2.*DLOG((1.-X20)/X10/X20)
     3        *HQQ)/(X1-X10)+2.*DLOG(X1-X10)/(X1-X10)*(Z*HQ1Q(I)-HQQ)
     4        +HQ1Q(I)*(X1-X10)/X1**2*(1.+DLOG(X1-X10))
      else
         TEMP=((1.+Z*Z)*DLOG(2.*(1.-X20)/X20/(X1+X10))*HQ1Q(I)
     2         -2.*DLOG((1.-X20)/X10/X20)
     3         *HQQ)/(X1-X10)+2.*DLOG(X1-X10)/(X1-X10)*(Z*HQ1Q(I)-HQQ)
     4         +HQ1Q(I)*(X1-X10)/X1**2*(1.+DLOG(X1-X10))
      endif
      QQ2=QQ2+.5*(1.-X10)*WI(I)*FAC*ACF*TEMP*xffac
      HG1Q(I)=Q1(0,I)*(4.*Q20(1)+Q20(2)+Q20(3)+4.*Q20(-1)+Q20(-2)+
     2                 Q20(-3)+8.*Q20(4)+2.*q20(5))/9.
      if(ixfx.eq.1)then
         TEMP=((X10**2+(X1-X10)**2)/(2.*X1**3)*
     2         DLOG((x10+x20)*(X1-X10)*(1.-X20)/x10/X20/(X1+X20))
     3         +X10*(X1-X10)/X1**3)*HG1Q(I)
      else
         TEMP=((X10**2+(X1-X10)**2)/(2.*X1**3)*
     2         DLOG(2.*(X1-X10)*(1.-X20)/X20/(X1+X10))
     3         +X10*(X1-X10)/X1**3)*HG1Q(I)
      endif
      QG2=QG2+.5*(1.-X10)*WI(I)*FAC*AT*TEMP*xffac
  200 CONTINUE
C        1         2         3         4         5         6         7 *
      QQ3=0.
      QG3=0.
      DO 300 I=1,NTERMS
      X2=.5*(1.-X20)*XI(I)+.5*(1.+X20)
      TEMPQQ=0.
      TEMPQG=0.
      TEMPGQ=0.
      DO 400 J=1,NTERMS
      X1=.5*(1.-X10)*XI(J)+.5*(1.+X10)
      HQ1Q2=(4.*Q1(1,J)*Q2(-1,I)+Q1(2,J)*Q2(-2,I)+Q1(3,J)*Q2(-3,I)
     2      +4.*Q1(-1,J)*Q2(1,I)+Q1(-2,J)*Q2(2,I)+Q1(-3,J)*Q2(3,I)
     3      +8.*Q1(4,J)*Q2(4,I)+2.*q1(5,j)*q2(5,i))/9.
      HQ1G2=(4.*Q1(1,J)+Q1(2,J)+Q1(3,J)+4.*Q1(-1,J)+Q1(-2,J)+Q1(-3,J)
     2       +8.*Q1(4,J)+2.*q1(5,j))*Q2(0,I)/9.
      HG1Q2=(4.*Q2(1,I)+Q2(2,I)+Q2(3,I)+4.*Q2(-1,I)+Q2(-2,I)+Q2(-3,I)
     2      +8.*Q2(4,I)+2.*q2(5,i))*Q1(0,J)/9.
      TEMP=HQ1Q2*HA(X1,X2,X10,X20)+(GA(X1,X2,X10,X20)*HQ1Q2
     2    +GA(X10,X20,X10,X20)*HQQ-GA(X1,X20,X10,X20)*HQ1Q(J)
     3    -GA(X10,X2,X10,X20)*HQQ2(I))/(X1-X10)/(X2-X20)
      TEMPQQ=TEMPQQ+.5*(1.-X10)*WI(J)*TEMP*FAC*ACF*2.
      TEMP=HQ1G2*HC(X1,X2,X10,X20)+(GC(X1,X2,X10,X20)*HQ1G2
     2+   -GC(X10,X2,X10,X20)*HQG2(I))/(X1-X10)
      TEMPQG=TEMPQG+.5*(1.-X10)*WI(J)*TEMP*FAC*AT
      TEMP=HG1Q2*HC(X2,X1,X20,X10)+(GC(X2,X1,X20,X10)*HG1Q2
     2+   -GC(X20,X1,X20,X10)*HG1Q(J))/(X2-X20)
      TEMPGQ=TEMPGQ+.5*(1.-X10)*WI(J)*TEMP*FAC*AT
  400 CONTINUE
      QQ3=QQ3+.5*(1.-X20)*WI(I)*TEMPQQ
      QG3=QG3+.5*(1.-X20)*WI(I)*(TEMPQG+TEMPGQ)
  300 CONTINUE
      THEORY=THEORY+QQ1+QG1+QQ2+QG2+QQ3+QG3
      GLUCOR=QG1+QG2+QG3
      RETURN
      END


      FUNCTION GA(X1,X2,X10,X20)
      implicit real*8 (a-h,o-z)
      common/e866xf/ixfx
      TAU=X10*X20
      if(ixfx.eq.1)then
         temp=(x1+x2)*(tau**2+(x1*x2)**2)
         ga=temp/2./(x1*x2)**2/(x1+x20)/(x2+x10)
      else
         TEMP=(X1*X2+TAU)*(TAU**2+(X1*X2)**2)
         GA=TEMP/(X1*X2)**2/(X1+X10)/(X2+X20)
      endif
      RETURN
      END
      FUNCTION HA(X1,X2,X10,X20)
      implicit real*8 (a-h,o-z)
      common/e866xf/ixfx
      TAU=X10*X20
      if(ixfx.eq.1)then
         ha=-1./x1/x2/(x1+x2)
      else
         HA=-2.*TAU*(X1*X2+TAU)/(X1*X2)/(X10*X2+X20*X1)**2
      endif
      RETURN
      END


      FUNCTION GC(X1,X2,X10,X20)
      implicit real*8 (a-h,o-z)
      common/e866xf/ixfx
c
c  this term corresponds to my notation with 1/(x1-x10)_R
c  so 1<-->2 compared to KLMP here and in HC
c 
      TAU=X10*X20
      if(ixfx.eq.1)then
         temp=tau**2+(x1*x2-tau)**2
         gc=temp/2./x1**2/x2**3/(x1+x20)
      else
         TEMP=X10*(X1*X2+TAU)*(TAU**2+(TAU-X1*X2)**2)
         GC=TEMP/X1**2/X2**3/(X1+X10)/(X10*X2+X20*X1)
      endif
      RETURN
      END
      FUNCTION HC(X1,X2,X10,X20)
      implicit real*8 (a-h,o-z)
      common/e866xf/ixfx
      TAU=X10*X20
      if(ixfx.eq.1)then
         temp=x2*(x1+x20)*(x1-x10)+2.*tau*(x1+x2)
         hc=temp/2./(x1*x2)**2/(x1+x2)**2
      else
         TEMP=X20*X2*X1**2+TAU*(X10*X2+2.*X1*X20)
         HC=TEMP*TAU*(X1*X2+TAU)/(X1*X2)**2/(X10*X2+X20*X1)**3
      endif
      RETURN
      END


      SUBROUTINE HOQUARK(INDEX,IFL,X,ALS,NDRV,QUARK)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION QUARK(-5:5)
      common/grpthy/flavor
      common/target/az,an
C
c  unfolds quark distributions
c  
      call pdf(ndrv,x,als,u,d,ub,db,sb,cb,bb,glue)
      QUARK(0)=GLUE/X
      quark(5)=bb/x
      QUARK(4)=cb/X
      QUARK(3)=sb/X
      quark(-5)=quark(5)
      quark(-4)=quark(4)
      QUARK(-3)=QUARK(3)
      IF(INDEX.EQ.1) THEN
         QUARK(-1)=ub/X
         QUARK(-2)=db/X
         QUARK(1)=u/x
         QUARK(2)=d/x
      ELSE IF(INDEX.EQ.2) THEN
         TUB2=ub/x
         TDB2=db/x
         TU2=u/x
         TD2=d/x
         QUARK(-1)=AZ*TUB2+AN*TDB2
         QUARK(-2)=AZ*TDB2+AN*TUB2
         QUARK(1)=AZ*TU2+AN*TD2
         QUARK(2)=AZ*TD2+AN*TU2
      ENDIF
      RETURN
      END




      SUBROUTINE PDF(NDRV,X,S,U,D,UB,DB,SB,CB,BB,GLUE)
      IMPLICIT REAL*8 (A-H,O-Z)
c      common/input/ins,nvl,nmax,delta,s0,iord
      common/threshold/sbbar
      CALL GINTERP(1,NDRV,X,S,UV)
      CALL GINTERP(2,NDRV,X,S,DV)
      CALL GINTERP(3,NDRV,X,S,UPLUS)
      CALL GINTERP(4,NDRV,X,S,DPLUS)
      CALL GINTERP(5,NDRV,X,S,SPLUS)
      call ginterp(6,ndrv,x,s,cplus)
      CALL GINTERP(7,NDRV,X,S,SING)
      CALL GINTERP(8,NDRV,X,S,GLUE)
c      q2=alam**2*dexp(s0*dexp(s))
c      al=alpha_s(iord+1,q2,alam,neff)
c      flavor=neff
      flavor=5.
      if(s.lt.sbbar)flavor=4.
      if(flavor.eq.4.)then
         bb=0.
      else if(flavor.eq.5.)then
         bplus=-(uplus+dplus+splus+cplus)
         bb=.5*(bplus+sing/flavor)
      endif
      cb=.5*(cplus+sing/flavor)
      sb=.5*(splus+sing/flavor)
      db=.5*(dplus-dv+sing/flavor)
      ub=.5*(uplus-uv+sing/flavor)
      u=uv+ub
      d=dv+db
C      write(65,*) x,u,d,sb,ub,db,sb,glue,cb
      RETURN
      END



      SUBROUTINE GINTERP(I,NDRV,X,S,ANS)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION F1(30),F2(30)
      COMMON/GAUS16/XI(16),WI(16),NTERMS,XX(17)
      COMMON/INPUT/INS,NVL,NMAX,DELTA,S0,IORD
      COMMON/MINIMB/CURPAR(60),UNCRT(60),ERMIN(40),PWATE(60),IFREP(40)
      COMMON/GINT/GF(8,40,1080)
      common/threshold/sb
C
C  THIS ROUTINE INTERPOLATES AS NEEDED IN THE ARRAY GF TO OBTAIN THE
C  EVOLVED DISTRIBUTIONS AT THE REQUIRED Q2 AND X VALUES. NOTE THAT ONLY
C  40 PARAMETERS CAN BE VARIED IN A GIVEN FIT. OF THESE, AT MOST 35 CAN 
C  BE PARAMETERS RELATED TO THE PARTON DISTRIBUTIONS. THIS IS GOVERNED BY THE
C  THIRD DIMENSION OF GF (1085 ABOVE). THIS IS 30+30*(NO. OF VARIED PARTON 
C  PARAMETERS).
C
      IS=S/DELTA+1
      isign=1
c
c  Adjust so that interpolation is not done across the b threshold
c
      if(s.lt.sb.and.(s+delta).gt.sb)isign=-1
      IS1=IS+isign
      if(s.gt.sb.and.(delta*(is-1)).lt.sb)then
         is=is+1
         is1=is+1
      endif   
      NDRV2=NDRV
      IF(IFREP(NDRV).GT.35) NDRV2=0
      DO 1 L=1,30
      KL=L+30*NDRV2
      F1(L)=GF(I,IS,KL)
      F2(L)=GF(I,IS1,KL)
    1 CONTINUE
      A1=GETFV(X,F1)
      A2=GETFV(X,F2)
C      A1=DLOG(A1)
C      A2=DLOG(A2)
      S1=(IS-1)*DELTA
      S2=S1+DELTA*isign
      ANS=A1*(S-S2)/(S1-S2)+A2*(S-S1)/(S2-S1)
C      ANS=DEXP(ANS)
      RETURN
      END 
      SUBROUTINE TEST(XC)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XC(60),FV(32),XMSTV(10),CN31(10),DN(10),DN1(10)
     2,V(4)
      COMMON/Q2STUF/Q02,Q2MAX 
      COMMON/INPUT/INS,NVL,NMAX,DELTA,S0,IORD
      COMMON/GAUS32/XI(32),WI(32),NTERMS,XX(33)
      COMMON/GRPTHY/FLAVOR
      COMMON/NORMAL/INORM(30)
      DATA CN31/0.,4.282,6.047,7.209,8.094,8.820,9.440,9.983,
     210.468,10.907/
      DATA DN/0.,.427,.667,.837,.971,1.080,1.174,1.255,1.327,1.392/
      COMMON/CONSTANTS/PI,PI2
C
C  THIS ROUTINE TESTS THE ACCURACY OF THE EVOLUTION ROUTINE BY
C  COMPARING THE EVOLVED MOMENTS WITH TH EMOMENTS OF THE EVOLVED
C  DISTRIBUTIONS. AT THIS TIME IT IS IMPLEMENTED ONLY FOR THE
C  NONSINGLET DISTRIBUTION X(UV+DV).
C  FOR TEST PURPOSES SOME SINGLET AND GLUON IDSTRIBUTIONS ARE
C  PRINTED OUT AS WELL.
C
C      IF(INS) 200,200,100
  100 CONTINUE
      WRITE(6,1001)
 1001 FORMAT(///,' NON-SINGLET MOMENTS TEST') 
      WRITE(6,1002)
 1002 FORMAT(5X,' N',5X,'EVOLVED MOMENT',1X,'INTQCD MOMENT',
     24X,'PERCENT ERROR') 
      B0=11.-2.*FLAVOR/3.
      B1=102.-38.*FLAVOR/3.
      DO 10 N=1,10
      DN1(N)=CN31(N)*B0/B1
   10 CONTINUE
c      T=DLOG(Q02/XC(1)**2)
c      A0=ALPHA(T)
      a0=alpha_s(iord+1,q02,xc(1),neff)
      flavor=neff
      ALB0=1.+B1/B0*A0/4./PI
C      DO 1 I=1,NTERMS
C      V(1)=XX(I)
C      V(2)=Q02
C      V(4)=7
C      FV(I)=THEORY(1,I,V,XC)
C    1 CONTINUE
C      DO 2 N=2,10
C      XMSTV(N)=XMNT32(FV,N)
C    2 CONTINUE
      C=3./(BETA(XC(2),XC(3)+1.)+XC(4)*BETA(XC(2)+1.,XC(3)+1.)
     2+XC(5)*BETA(XC(2)+2.,XC(3)+1.))
      DO 2 N=1,10
      XMSTV(N)=C*(BETA(XC(2)+N-1.,XC(3)+1.)+XC(4)*BETA(XC(2)+N,XC(3)+1.)
     2+XC(5)*BETA(XC(2)+N+1.,XC(3)+1.))/XC(INORM(7))
    2 CONTINUE
c      T=DLOG(Q2MAX/XC(1)**2)
c      AL=ALPHA(T)
      al=alpha_s(iord+1,q2,xc(1),neff)
      flavor=neff
      ALB=1.+B1*AL/B0/4./PI
      DO 4 N=1,10
      IF(N.LE.2) THEN
      DO 3 I=1,NTERMS
      V(1)=XX(I)
C
C  MODIFIED N=1 MOMENT INTEGRATION
C
      IF(N.EQ.1) V(1)=XX(I)**2
      V(2)=Q2MAX
      V(4)=7
      FV(I)=THEORY(1,I,V,XC)
    3 CONTINUE
      ENDIF
      XMV=XMNT32(FV,N)
      S1=SS(N,1)
      S2=SS(N,2)
      CN=4./3.*(2.*S1**2-2.*S2+3.*S1-2.*S1/N/(N+1.)
     2+1./N+2./(N+1.)+2./N/N-9.)
      IF(IORD) 6,6,7
    6 CONTINUE
      ALB=1.
      ALB0=1.
      CN=0.
    7 CONTINUE
      TV=XMSTV(N)*(AL/A0)**DN(N)*(ALB/ALB0)**(DN1(N)-DN(N)) 
C      TV=TV*(1.+CN*AL/4./PI)/(1.+CN*A0/4./PI)
      TV=TV*(1.+CN*AL/4./PI)
      ERV=(TV-XMV)/TV*100.
    4 WRITE(6,5) N,TV,XMV,ERV
    5 FORMAT(5X,I2,5X,3E15.4) 
      RETURN
C  200 CONTINUE
C      WRITE(6,201)
C  201 FORMAT(///,' SINGLET AND GLUON TEST OUTPUT')
C      ITIME=1
C      Q2=Q02
C      DO 211 K=1,2
C      S=DLOG(DLOG(Q2/XC(1)**2)/DLOG(Q02/XC(1)**2))
C      WRITE(6,202) Q2
C  202 FORMAT(///,' Q2=',F8.2) 
C      WRITE(6,203)
C  203 FORMAT(/,5X,' X',8X,'SINGLET',9X,'GLUON',12X,'F2',
C     212X,'XF3',12X,'F2EN') 
C      DO 210 J=1,19 
C      X=.05*J
C      V(1)=X
C      V(2)=Q2
C      V(4)=6
C      F2=THEORY(1,ITIME,V,XC) 
C      V(4)=7
C      ITIME=2
C      XF3=THEORY(1,ITIME,V,XC)
C     V(4)=3
C     F2EN=THEORY(1,ITIME,V,XC)
C      F2EN=0.
C      CALL GINTERP(7,0,X,S,SING)
C      CALL GINTERP(8,0,X,S,GLUE)
C      WRITE(6,204) X,SING,GLUE,F2,XF3,F2EN 
C  204 FORMAT(F10.4,5E15.4)
C  210 CONTINUE
C      Q2=Q2MAX
C  211 CONTINUE
C      RETURN
      END 


      FUNCTION SS(N,I)
      implicit real*8 (a-h,o-z)
      SS=0.
      DO 1 J=1,N
      SS=SS+1./J**I 
    1 CONTINUE
      RETURN
      END 


      SUBROUTINE PLOTLL(X1,X2,Y1,Y2,X,Y,ER,XTH,TH,NEXP,NTH) 
      implicit real*8 (a-h,o-z)
      REAL LN(100)
      DIMENSION X(1),Y(1),ER(1),XTH(1),TH(1)
      DATA AB,AS,AX,AE,AP,AI/1H ,1H*,1HX,1H-,1H+,1HI/
C
C  THIS IS A RELATIVELY CRUDE PLOTTER WHICH WORKS ON A
C  LINE PRINTER. IT CAN BE MODERATELY USEFUL BUT IT
C  DOES RESULT IN A LOT OF OUTPUT IF YOU FIT A LOT OF
C  DIFFERENT X POINTS.
C
      B=100./(DLOG(X2)-DLOG(X1))
      A=1.-B*DLOG(X1)
      D=49./(DLOG(Y1)-DLOG(Y2))
      C=1.-D*DLOG(Y2)
      WRITE(6,9)
    9 FORMAT(1H1)
      DO 1 I=1,100
    1 LN(I)=AS
      DO 2 I=10,90,10
    2 LN(I)=AP
      WRITE(6,7) (LN(I),I=1,100) 
    7 FORMAT(1H ,15X,100A1)
      DO 100 I=1,50 
      LN(1)=AS
      LN(100)=AS
      DO 5 J=2,99
    5 LN(J)=AB
      YY=EXP((I-C)/D)
      DO 99 K=1,NEXP
      IX=A+B*DLOG(X(K))
      IYU=C+D*DLOG(Y(K)+ER(K))
      IY=C+D*DLOG(Y(K))
      IF(ER(K).GE.Y(K)) ER(K)=Y(K)-.001 
      IYL=C+D*DLOG(Y(K)-ER(K))
      IF(I.EQ.IYU.OR.I.EQ.IYL) LN(IX)=AE
      IF(I.EQ.IY) LN(IX)=AX
   99 CONTINUE
      DO 98 L=1,NTH 
      IX=A+B*DLOG(XTH(L))
      ITH=C+D*DLOG(TH(L))
      IF(I.EQ.ITH) LN(IX)=AS
   98 CONTINUE
  100 WRITE(6,3) YY,(LN(K),K=1,100)
    3 FORMAT(1H ,F10.4,5X,100A1)
      DO 6 I=1,100
      LN(I)=AS
      IP=I/10*10
      IF(I.EQ.IP) LN(I)=AP
    6 CONTINUE
      LN(100)=AS
      WRITE(6,7) (LN(K),K=1,100) 
      WRITE(6,8) X1,X2 
    8 FORMAT(1H ,10X,F10.3,90X,F10.3)
      RETURN
      END 



      SUBROUTINE WATE4
      implicit real*8 (a-h,o-z)
      COMMON/GAUSS4/XI(4),WI(4),NTERMS,XX(5)
      NTERMS=4
C
C  4 POINT GAUSSIAN QUADRATURE
C
      XI(1)=-.8611363116
      XI(2)=-.3399810436
      XI(3)= .3399810436
      XI(4)= .8611363116
      WI(1)=.3478548451
      WI(2)=.6521451549
      WI(3)=.6521451549
      WI(4)=.3478548451
      DO 1 J=1,4
    1 XX(J)=.5*(XI(J)+1.)
      XX(5)=1.
      RETURN
      END 


      subroutine wasym(rs,y,ALS,NDRV,asym,it)
      IMPLICIT REAL*8 (A-H,O-Z)
      dimension xa(24),argp(24),argm(24),qa(-5:5),qb(-5:5)
C      common/qs/q2
C      common/cteqpdf/ist
c
c  calculates the w lepton asymmetry. See Barger and Phillips 
c  "Collider Physics" pgs 254-256. Note: theta-hat on pgs 255, 256 
c  is - (theta-hat) on page 254!  
c
c  J.F. Owens    June 30, 1994
c
c  rs = sqrt(s)  y=lepton rapidity
c
      Ak1ResCt4m(x)= 1.07949- .0041764*x+ .017066*x**2+ .016434*x**3
     >    + .00074073*x**4 - .0015393*x**5 - .00040555*x**6
      CDF1800(x)= 1.0904 + 0.02835 * x + 0.023048 * x**2 
     >- 0.0075403 * x**3 - 0.0022601 * x**4 + 0.0020092 * x**5
      CDF1960(x) = 1.0795 + 0.014884 * x + 0.025591 * x**2 
     >+ 0.00015113 * x**3 - 0.0029232 * x**4 + 0.00082595 * x**5
      D01960(x) = 1.1209 + 0.022313 * x + 0.018531 * x**2 
     >- 0.001236 * x**3 - 0.0019652 * x**4 + 0.00041654 * x**5

c  cteq2m flag
      ist=1
c  ckm factors (using just the Cabibbo angle)
      c2tc=.95
      s2tc=.05
c  w mass
      xmw=80.
      q2=xmw**2
      s=rs**2
c
c  CDF uses pt(muon)>25 x=xmw/(2.*ptmin)
c  D0 uses 20

      x=xmw/50.
      if(it.eq.129) x=xmw/40.
      w=x+sqrt(x**2-1.)
      xi=Dlog(w)
c
c  In lowest order, the pt cut translates into simple 
c  limits on the xa integral
c
      xamin=xmw/rs*dexp(y-xi)
      xamax=xmw/rs*dexp(y+xi)
c
c  6-pt gaussian quadrature routine -- supplied below
c
      call gq11(xamin,xamax,4,xa,argp,ansp)
      do 100 j=1,24
      xb=xmw**2/(s*xa(j))
c
c  my parton distribution calls
c  replace with yours
c
C      call dist(xa(j),qa)
C      call dist(xb,qb)
      CALL HOQUARK(1,16,XA(J),ALS,NDRV,QA)
      CALL HOQUARK(1,16,XB,ALS,NDRV,QB)
      yh=y-.5*Dlog(xa(j)/xb)
      st=1./cosh(yh)
      ct=tanh(yh)
      facp=(1.-ct)**2*st**2
      facm=(1.+ct)**2*st**2
      argp(j)=c2tc*((qa(1)*qb(2)+qa(4)*qb(3))*facp
     2+(qa(-2)*qb(-1)+qa(-3)*qb(4))*facm)
     3+s2tc*((qa(1)*qb(3)+qa(4)*qb(2))*facp
     4+(qa(-3)*qb(-1)+qa(-2)*qb(4))*facm)
      argm(j)=c2tc*((qa(2)*qb(1)+qa(3)*qb(4))*facm
     2+(qa(-1)*qb(-2)+qa(4)*qb(-3))*facp)
     3+s2tc*((qa(3)*qb(1)+qa(2)*qb(4))*facm
     4+(qa(-1)*qb(-3)+qa(4)*qb(-2))*facp)
      argp(j)=argp(j)/xa(j)
      argm(j)=argm(j)/xa(j)
  100 continue
      call gq11(xamin,xamax,0,xa,argp,ansp)
      call gq11(xamin,xamax,0,xa,argm,ansm)
      if(it.eq.126)then
         ansp=ansp*Ak1ResCt4m(y)
         ansm=ansm*Ak1ResCt4m(-y)
      else if(it.eq.127)then
         ansp=ansp*CDF1800(y)
         ansm=ansm*CDF1800(-y)
      else if(it.eq.128)then
         ansp=ansp*CDF1960(y)
         ansm=ansm*CDF1960(-y)
      else if(it.eq.129)then
         ansp=ansp*D01960(y)
         ansm=ansm*D01960(-y)
      endif
      asym=(ansp-ansm)/(ansp+ansm)
      return
      end
      SUBROUTINE GQ11(XMIN,XMAX,N,X,Y,W)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(24),Y(24),V(6),R(6) 
      DATA V/0.03376524,0.16939531,0.38069041,0.61930959,0.83060469,
     10.96623476/,R/0.08566225,0.18038079,0.23395697,0.23395697,
     20.18038079,0.08566225/
      save d,nsave
      IF(N) 250,250,150 
  150 K=0 
      NSAVE=N 
      D=(XMAX-XMIN)/N 
      XL=XMIN-D 
      DO 200 I=1,N
      XL=XL+D 
      DO 200 J=1,6
      K=K+1 
  200 X(K)=XL+D*V(J)
      RETURN
  250 W=0.
      K=0 
      DO 300 I=1,NSAVE
      DO 300 J=1,6
      K=K+1 
  300 W=W+Y(K)*R(J) 
      W=W*D 
      RETURN
      END 
      FUNCTION BETA(X1,X2)
      IMPLICIT REAL*8 (A-H,O-Z)
C  EULER BETA FUNCTION SUBROUTINE
      CALL GAMMA(X1,G1,IER)
      CALL GAMMA(X2,G2,IER)
      X3=X1+X2
      CALL GAMMA(X3,G3,IER)
      BETA=G1*G2/G3 
      RETURN
      END 
