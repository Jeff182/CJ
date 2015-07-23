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
c  Aug2008: modified to include gamma+jet cross sections (Jeff)
c           modified to include TMC in DIS (Alberto)
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
      PROGRAM ALTPFIT
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*10 IPNAME(60),ITITLE(30),datname(30),IQUIT
     2,LABEL,LABEL2
      character lcomment*80, car*9, input*50, infile*50, string*100
      integer len,jtgt,ioff,IUout
      double precision sgntot,restot

      logical set
      dimension datchis(30),datchis_w(30)
     &     ,datresid(30),datresid_w(30)
     &     ,datsignd(30),datsignd_w(30)
     &     ,idatpts(30)
     &     ,xdum(100)
      DIMENSION DATA(6,4500),xerr(4500),V(4),NFLAG(30),TX(1500)
     &     ,TQ2(1500)
      DIMENSION XQ(1500),YDAT(1500),ERD(1500),Q2TH(1500),YTH(1500)
     2,YOUT(1500),CHI(1500),THOUT(1500),ITOUT(1500),itype(30),errn(30)
     3,ITGT(30),ICORR(30),ITMC(30),iht(30)
      DIMENSION Q2CUT(30),XPOS(350),INORM(30) 
     2,RSPOS(15),XLO(30),XHI(30),W2MIN(30) 
     3,XMLO(15) 

*    *** common blocks 

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
      COMMON/FLAGS/ITYPE,NFLAG,INORM,ITGT,ICORR,ITMC,iht
      common/qcdpar/al,nf,norder,set
      common/normalization/nfit
      common/cdf/covinv(33,33),cdfchisqr,ncdf,ncdfl,ncov

      Character*50 output
      common/outputfile/output
      DATA IQUIT/'END'/

*     Higher-twist parameters
      integer nht(4)
      common/ht/nht

*     vars for pdf errors
      character*80, file1
      character*4, tmp
      dimension derv(-5:2,50,35), err(-5:2,50), pdf0(-5:2,50)
     &     , xgrid(50),q2grid(30)
      DATA NX,XGRID/30,.0001,.0002,.0004,.0006,.0008,.001,.002,.004,
     2.008,.016,.032,.064,.1,.15,.2,.25,.3,.35,.4,.45,.5,.55,.6,.65,
     3.7,.75,.80,.85,.9,.95,20*0./
*     finer xgrid, corresponding to the xgrid used in the dglap evolution
*      DATA NX,XXGRID/50,.0001,.0002,.0004,.0006,.0008,.001,.002,.004,
*     2.008,.016,.032,.064,.1,.13,.16,.19,.22,.25,.28,.31,.34,.37,.40,
*     3.43,.46,.49,.52,.55,.58,.61,.64,.67,.70,.73,.76,.79,.82,.84,.86,
*     4.88,.90,.91,.92,.93,.94,.95,.96,.97,.98,.99/
      DATA NQ2,Q2GRID/29,
     &     1.7,2.0,2.5,3.2,4.0,5.0,6.3,7.9,10
     &     ,13,17,20,25,32,40,50,63,79,100
     &     ,130,170,200,250,320,400,500,630,790,1000,0 /
      common/gint/gf(8,40,1800)
      common/errors/cov(40,40)

*     Q2 values for PDF output
      double precision q2pdf(10),xpdf(10)
      integer nq2pdf,nxpdf
      data  q2pdf/ 1.7,2.,10.,25.,64.,100.,1000.,3*0 /
      data nq2pdf/ 7 /
      data   xpdf/ 0.1, 0.3, 0.5, 0.7, 0.85,5*0 /
      data  nxpdf/ 5 /
C     

C  *****************************************************************
C

*    *** reads input file name from std input and opens it

      write(6,11)
 11   format(/' Enter input filename for fit >')
      read(5,*) infile
      call trmstr(infile,len1)
      if (infile(len1-3:len1).eq.'.dat') len1=len1-4

      input = infile(1:len1)//'.dat'
      open(unit=10,file=input(1:len1+4),status='old')

*    *** opens output files

      write(6,*) 'Output in '''//infile(1:len1)//'.###'''

      tmp='.out'
      len2=4
      file1=infile(1:len1)//tmp(1:len2)
      len=len1+len2
      IUout=4
      open(IUout,file=file1(1:len),status='new')

      if (infile(1:4).ne.'comp') then
         tmp='.pdf'
         call trmstr(tmp,len2)
         file1=infile(1:len1)//tmp(1:len2)
         len=len1+len2
         open(unit=3,file=file1(1:len),status='new')
      end if

C 
C  READ IN FIT PARAMETERS 
C

      read(10,*) string

      I=0 
   32 I=I+1 
      READ(10,*) IPNAME(I),CURPAR(I),PWATE(I)
      IPAR(I)=I 
      IF(IPNAME(I).EQ.IQUIT) GO TO 31 
*    ... assign higher-twist parameter numbers
      if(ipname(i).eq.'ht1') then
         nht(1)=i
      else if(ipname(i).eq.'ht2') then
         nht(2)=i
      else if(ipname(i).eq.'ht3') then
         nht(3)=i
      else if(ipname(i).eq.'ht4') then
         nht(4)=i
      end if
      GO TO 32
   31 NPARAM=I-1

C 
C  READ IN TYPES OF DATA TO BE FIT AND THE CUTS 
C 
      I=0 
   33 I=I+1 
      READ(10,*) ITYPE(I),NFLAG(I),ITITLE(I),INORM(I),ITGT(I),
     $Q2CUT(I),XLO(I),XHI(I),W2MIN(I),icorr(i),itmc(i),iht(i) 
C
C  TERMINATE LIST WITHs ITYPE = 999
C
      IF(ITYPE(I).EQ.999) GO TO 34 
      if(itype(i).le.100)then
         WRITE(IUout,38) ITITLE(I),Q2CUT(I),XLO(I),XHI(I),W2MIN(I) 
 38      FORMAT(1H ,A10,' DIS data with Q2 > ',F6.2,
     2        ' and ',F5.3,' < xB < ',F5.3,' and W2 >= ',F6.2) 
         if(nflag(i).eq.2.or.nflag(i).eq.9)then
            call split_nuke(itgt(i),ibj,ioff,jtgt)  
                             ! selects Bj limit, off-shell and smearing model
            if(jtgt.eq.0)then
               write(IUout,47)
            else if(jtgt.eq.1)then
               write(IUout,48)
            else if(jtgt.eq.2)then
               write(IUout,50)
            else if(jtgt.eq.3)then
               write(IUout,60)
            endif
            if (ibj.ge.1) then
               write(IUout,70)
            end if
            if (ioff.eq.1) then
               write(IUout,61) 
            else if (ioff.eq.2) then
               write(IUout,62) 
            else if (ioff.eq.3) then
               write(IUout,63) 
            else if (.not.ioff.eq.0) then
               write(*,*) 'ERROR(altpfit): ioff out of range =',ioff
     &              ,itgt(i),i
               stop
            end if
         else if(nflag(i).eq.5.or.nflag(i).eq.8.or.nflag(i).eq.10
     2           .or.nflag(i).eq.11)then
            if(itgt(i).eq.0)then
               write(IUout,47)
            else if(itgt(i).eq.-1)then
               write(IUout,35)
            else if(itgt(i).eq.-2)then
               write(IUout,49)
            else if(itgt(i).eq.-3)then
               write(IUout,52)              
            endif 
         endif
 47      format(1h ,'  No nuclear correction applied')
 48      format(1h ,'  Density Model nuclear correction applied')
 50      format(1h ,'  Nucleon smearing correction - KMK')
 60      format(1h ,'  Nucleon smearing correction - LC Frank/Strik')
 61      format(1h ,'  Off-shell corr. by Kulagin-Petti')
 62      format(1h ,'  Off-shell corrections by W.Melnitchouk')
 63      format(1h ,'  Bjorken limit for on-shell Deuteron corrections')
 49      format(1h ,'  EmcNmc and Gomez nuclear corrections applied')
 52      format(1h ,'  Kulagin-Petti nuclear correction applied')
 35      format(1h ,'  EmcNmc nuclear corrections applied')
 70      format(1h ,'  Bjorken limit for nuclear corrections imposed')
 
        if (itmc(i).eq.1) then 
            write(IUout,53) 
         else if (itmc(i).eq.2) then 
            write(IUout,54)
         else if (itmc(i).eq.3) then 
            write(IUout,55)
         end if
 53      format(1h ,'  GP Target Mass Corrections applied')
 54      format(1h ,'  CF Target Mass Corrections applied')
 55      format(1h ,'  NV Target Mass Corrections applied')
         if (iht(i).eq.1) then 
            write(IUout,56) 
         else if (iht(i).eq.11) then 
            write(IUout,58)
         else if (iht(i).ne.0) then
            write(IUout,100), 'iht out of range:', iht
            stop
         end if
 56      format(1h ,'  HT corrections, C=C(xB)')
 58      format(1h ,'  HT corrections, C=C(xi)')
      else if(itype(i).le.200)then
         write(IUout,39) ititle(i),q2cut(i)
 39      format(1h ,a10,' vector boson data with Q2 >=', 
     2        f7.2)
      else if(itype(i).le.300)then
         if(icorr(i).eq.0)then
            write(IUout,37) ititle(i)
 37         format(1h ,a10,' Single jet inclusive data')
         else
            write(IUout,36) ititle(i)
 36         format(1h ,a10,
     2           ' Single jet inclusive data with correlated errors')
         endif
      else if(itype(i).le.400)then
         WRITE(IUout,40) ITITLE(I)
 40      FORMAT(1H ,A10,' Photon data -- cuts to be written yet ...')
      endif
 100  format('ERROR (altpfit08):',X,A,I3)
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
      READ(10,*) INS,NVL,IORD,Q02,ISERR
      norder=iord+1
C  
C  THE UNIVERSAL PARTON DISTRIBUTION IS USED HERE
C
c      IORDP=IORD
C  IORDP PASSED THROUGH PATRICK'S COMMON BLOCK
      IF(ISERR.NE.0) WRITE(IUout,1004)
 1004 FORMAT(/' SYSTEMATIC ERRORS ADDED IN QUADRATURE')
c
      read(10,51) lcomment
      write(IUout,*)
      write(IUout,51) lcomment
 51   format(A80)

      close(10)


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
      CALL WATE96
      CALL WATE8L 
c      DELTA=0.10
      DELTA=.04
      IERR=0
      IDRV=2
      NSTEPS=50 
      TOLRNC=0.5
      IPRNT=1 
      IPLOT=0
c
c initialize parameters for the jet routines
c
      call setqcd


*     Read experimental data files
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
c         7        = 1  is  JLab E-00-106 hydrogen
c         8        = 2  is  JLab E-00-106 deuterium
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
      call fildat(nfit,npts,data,xerr,q2max)
      
      call rdxpos(xpos,nxpos,rspos,nrspos,npts,data,ndis,ndy,njet,
     2ngamjet)
c
      IF(IORD) 42,42,43 
 42   WRITE(IUout,44)
 44   FORMAT(///,' THIS IS A LEADING ORDER FIT')
      GO TO 45
 43   WRITE(IUout,46)
 46   FORMAT(///,' THIS IS A NEXT TO LEADING ORDER FIT')
 45   CONTINUE
C 
C  THE ROUTINE TEST PERFORMS VARIOUS TESTS ON THE ALTARELLI-
C  PARISI EQUATION INTEGRATOR.
C 
C      CALL TEST(CURPAR)
C
C  CLOSE DATA FILES
C
      WRITE(IUout,71) Q02
 71   FORMAT(/,' Q02=',F8.2)

c     
c     prepare for fitting
c     

      FIT=0.
      DO 66 IF=1,NPARAM 
         FIT=FIT+PWATE(IF)
 66   continue
      if(fit.eq.0.d0)then
         V(1)=DATA(3,1)
         V(2)=DATA(4,1)
         v(3)=data(5,1)
         V(4)=DATA(6,1)
         DUM=THEORY(1,1,V,CURPAR)
         WRITE(IUout,74)
 74      FORMAT(/,' INPUT PARAMETERS') 
         DO 72 L=1,NPARAM
            WRITE(IUout,73) IPNAME(L),CURPAR(L)
 72      continue
 73      FORMAT(1H ,A10,5X,E15.4)
      else
         CALL MINIM(NV,NPTS,DATA,IDRV) 
      endif

C 
C  PRINT OUTPUT NOW.
C

      ndat = 0

*    *** DIS output

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
      CHITOT_w=0. 
      NEXP=0
      CHITEMP=0.
      CHITEMP_w=0.
      schi=0.
      schi_w=0.
      resid=0.
      resid_w=0.
      oldit=0
      ND = 1
      NDU = NPTS
      DO K=1,NDU+1      
         IEXP=DATA(6,K)
         IT=ITYPE(IEXP)
         nfl=nflag(iexp)
         if(it.ne.oldit)then
            if(oldit.ne.0)then
               rat=chitemp/nexp
               rat_w=chitemp_w/nexp
               datname(ndat) = ititle(data(6,K-1))
               idatpts(ndat) = nexp
               datchis(ndat) = chitemp
               datchis_w(ndat) = chitemp_w
               datsignd(ndat) = schi
               datsignd_w(ndat) = schi_w
               datresid(ndat) = resid
               datresid_w(ndat) = resid_w
               chitot=chitot+datchis(ndat)
               chitot_w=chitot_w+datchis_w(ndat)
               if (chitemp.eq.chitemp_w) then
                  write(IUout,26) chitemp,nexp,rat
               else
                  write(IUout,27) chitemp,chitemp_w,nexp,rat,rat_w
               end if
               write(IUout,28) schi
               write(IUout,29) resid
               if(oldit.eq.201.and.ncov.ne.0)then
*                ... cdf chi^2 for jets
                  datchis(ndat) = cdfchi2(data,NV,NPTS,curpar,xdum)
                  write(IUout,304)datchis(ndat)
 304              format(1h ,'cdf chi^2 with correlations = ',f10.4)
                  datchis_w(ndat) = datchis(ndat)
                  datsignd(ndat) = 0d0
                  datsignd_w(ndat) = 0d0
                  datresid(ndat) = 0d0
                  datresid_w(ndat) = 0d0
               end if
            endif
            ndat = ndat + 1
            oldit=it
            chitemp=0.
            chitemp_w=0.
            schi=0.
            schi_w=0.
            resid=0.
            resid_w=0.
            nexp=0
            if (k.le.NDIS) then
*             ... DIS headers
               write(IUout,83)
 83            FORMAT(///,5x,'X',6x,'Q2',6x,'THEORY',9x,'DATA'
     &              ,9x,'ERROR',8x,'chi^2',5x,'_w',4x,'ITYPE')
            else if (k.le.NDIS+NDY) then
*             ... DY headers
               if (nfl.eq.1.or.nfl.eq.3.or.nfl.eq.5) then
                  index=v(4)
                  if(IT.eq.108.or.IT.eq.110)then
                     write(IUout,244)
                  else
                     WRITE(IUout,241)
                  endif 
               else 
                  WRITE(IUout,244) 
               end if
 241           FORMAT(///,' SQRT(TAU)',5X,'Y',6X,'THEORY',6X
     2              ,'DATA',5X,'ERROR',5X,'chi2    _w',4X,'ITYPE')
 244           FORMAT(///,' SQRT(TAU)',5X,'XF',5X,'THEORY',6X
     2              ,'DATA',5X,'ERROR',5X,'chi2    _w',4X,'ITYPE') 
            else if (k.le.NDIS+NDY+NJET) then
*             ... JET headers
               write(IUout,303)
 303           format(///,3x,' ymin',5x,'ymax ',5x,'  ET ',6x,'Theory'
     $              ,8x,'Data',7x,'Error',5x,'chi2    _w',5x,'Itype')       
            else if (k.le.NDIS+NDY+NJET+ngamjet) then
*             ... GAMMA+JET headers
               write(IUout,403)
 403           format(///,3x,'  rs ',5x,' pt  ',6x,'Theory',8x,
     $              'Data',7x,'Error',5x,'chi2    _w',5x,'Itype')       
            endif
         end if
         if (k.le.NDU) then
            nexp=nexp+1
            V(1)=DATA(3,K)
            V(2)=DATA(4,K)
            v(3)=data(5,k)
            V(4)=data(6,k)
            TH=THEORY(1,2,V,CURPAR)
            chisq=chi2(th,DATA(1,k),xerr(k))
            chisq_w=chi2(th,DATA(1,k),DATA(2,k))
            CHITEMP=CHITEMP+CHISQ
            CHITEMP_w=CHITEMP_w+CHISQ_w
            IS=-1
            IF(TH.GE.DATA(1,K)) IS=1
            CHISQ=IS*CHISQ
            CHISQ_w=IS*CHISQ_w
            schi=schi+chisq
            schi_w=schi_w+chisq_w
            resid=resid+residual(th,DATA(1,k),xerr(k))
            resid_w=resid_w+residual(th,DATA(1,k),DATA(2,k))
*          ... data output
            if (k.le.NDIS) then
*             ... writes DIS data
               if (chisq_w.eq.chisq) then
                  write(IUout,84) data(3,k),data(4,k),th,data(1,k)
     2                 ,xerr(k),chisq,ititle(data(6,k))
 84               FORMAT(f8.4,F9.3,3E14.5,F7.2,10X,a10)
               else
                  write(IUout,85) data(3,k),data(4,k),th,data(1,k)
     2                 ,xerr(k),chisq,chisq_w,ititle(data(6,k))
 85               FORMAT(f8.4,F9.3,3E14.5,2F7.2,3X,a10)
               end if  
            else if (k.le.NDIS+NDY) then
*             ... writes DY data
               xqq = data(4,k)/sqrt(data(3,k))
               ydatt=data(1,k)
               erdd=data(2,k)
               itoutt=data(6,k)
               youtt = data(5,k)
               if (nfl.eq.1.or.nfl.eq.3.or.nfl.eq.5) then
*                ... data and chi^2 for nflag = 1,2,5
                  if (chisq_w.eq.chisq) then
                     WRITE(IUout,242) xqq,youtt,TH
     2                    ,YDATT,xerr(k),CHISQ,ititle(data(6,k))
 242                 FORMAT(5F10.4,F7.2,10X,a10)
                  else
                     WRITE(IUout,243) xqq,youtt,TH
     2                    ,YDATT,xerr(k),erdd,CHISQ,ititle(data(6,k))
 243                 FORMAT(5F10.4,2F7.2,3X,a10)
                  end if
               else 
*                ... data and chi^2 nflag = 2,4
                  if (chisq_w.eq.chisq) then
                     WRITE(IUout,245) YOUTT,TH,YDATT
     2                    ,xerr(k),CHISQ,ITITLE(data(6,k))
 245                 FORMAT(6F10.4,5X,a10)
                  else
                     WRITE(IUout,245) YOUTT,TH,YDATT
     2                    ,xerr(k),erdd,CHISQ,ITITLE(data(6,k))
                  end if
               end if
            else if (k.le.NDIS+NDY+NJET) then 
*             ... writes JET data
              if (chisq_w.eq.chisq) then
                  WRITE(IUout,342) data(3,k),data(4,k),data(5,k)
     2                 ,TH,data(1,k),xerr(k),CHISQ,ITITLE(data(6,k))
 342              FORMAT(3F10.4,3e12.4,f7.2,10X,a10)
               else
                  WRITE(IUout,343) data(3,k),data(4,k),data(5,k)
     2                 ,TH,data(1,k),xerr(k),data(2,k)
     2                 ,CHISQ,ITITLE(data(6,k))
 343              FORMAT(3F10.4,3e12.4,2f7.2,3X,a10)
               end if
            else if (k.le.NDIS+NDY+NJET+NGAMJET) then
*             ... writes GAMMA+JET data
               if (chisq_w.eq.chisq) then
                  WRITE(IUout,442) data(3,k),data(4,k),TH,data(1,k),
     2                 xerr(k),CHISQ,ITITLE(data(6,k))
 442              FORMAT(2F10.4,3e12.4,f7.2,10X,a10)
               else
                  WRITE(IUout,443) data(3,k),data(4,k),TH,data(1,k),
     2                 xerr(k),CHISQ,ITITLE(data(6,k))
 443              FORMAT(2F10.4,3e12.4,2f7.2,3X,a10)
               end if
            end if
         else
            ndat = ndat-1
            oldit = 0
         end if
      end do
 26   FORMAT(///,' CHI SQUARE=',F8.2,' FOR ',I4,' POINTS',
     2     '   CHISQ/NPTS=',f8.2)
 27   FORMAT(///,' CHI SQUARE=',F8.2,' /',F8.2,' FOR ',I4,' POINTS',
     2     '   CHISQ/NPTS=',f8.2,' /',f8.2)
 28   format(' Signed chisquare total =',f8.2)
 29   format(' Modified residual =',f8.2)


c
c     write out summary chi square: total & for each data set 
c     

      write(IUout,315)
 315  format(//,'Summary normalization chi2 contributions')
      write(IUout,314)
 314  format(/,1x,'data set',5x,'normalization',8x,'error',
     2     9x,'chisq')
      xnormtot=0d0
      do 311 i=1,nfit
         if(errn(i).eq.0.d0.or.inorm(i).eq.60)goto311
         chitemp=chi2(curpar(inorm(i)),1d0,errn(i))
         write(IUout,312) ititle(i),curpar(inorm(i)),errn(i),chitemp
 312     format(A10,3f15.4)
         xnormtot=xnormtot+chitemp
 311  continue

*    ... data chi^s
      npts     = 0
      chitot   = 0d0
      chitot_w = 0d0
      sgntot   = 0d0
      restot   = 0d0
      write(IUout,316)
 316  format(//,'Summary chi2 per data set')
      write(IUout,317)
 317  format(/,' data set      chisq    _w     signd    _w'
     &     ,'     resid     _w   npts   chisq/pts  _w')
      do i=1,ndat
         if( datchis(i).eq.datchis_w(i)
     &        .and.datsignd(i).eq.datsignd_w(i)
     &        .and.datresid(i).eq.datresid_w(i) ) then
*            print*, '* data = ',datname(i)
            write(IUout,318) datname(i),datchis(i),datsignd(i)
     &        ,datresid(i),idatpts(i),datchis(i)/idatpts(i)
         else
            write(IUout,319) datname(i)
     &           ,datchis(i),datchis_w(i)
     &           ,datsignd(i),datsignd_w(i)
     &           ,datresid(i),datresid_w(i)
     &           ,idatpts(i)
     &           ,datchis(i)/idatpts(i),datchis_w(i)/idatpts(i)
         end if
         chitot = chitot + datchis(i)
         chitot_w = chitot_w + datchis_w(i)
         sgntot = sgntot + datsignd(i)
         restot = restot + datresid(i)
         npts = npts + idatpts(i)
      end do
      write(IUout,320)
      if ( chitot.eq.chitot_w
     &     .and.sgntot.eq.sgntot_w
     &     .and.restot.eq.restot_w ) then
         WRITE(IUout,318) 'TOTAL     ',CHITOT,sgntot,restot
     &        ,NPTS,chitot/npts
      else
         WRITE(IUout,319) 'TOTAL     ',CHITOT,CHITOT_w
     &        ,sgntot,sgntot_w,restot,restot_w
     &        ,NPTS,chitot/npts,chitot_w/npts
      end if
      chitotnorm = CHITOT+xnormtot
      chitotnorm_w = CHITOT_w+xnormtot
      WRITE(IUout,321) 'TOTAL+norm',CHITOTnorm,CHITOTnorm_w
     &     ,NPTS,CHITOTnorm/npts,CHITOTnorm_w/npts

 318  format(A10,2X,3(F8.1,8X),I6,2F8.2)
 319  format(A10,2X,6F8.1,I6,2F8.2)
 320  format('--------------------------------------------------------'
     &     ,'----------------------------')
 321  format(A10,2X,2F8.1,32X,I6,2F8.2)



*    *** writes PDF and errors to ##.pdf file 
      if (infile(1:4).ne.'comp') then

         write(3,*) 'PDFs and errorbands'
         write(3,*) lcomment

*       ... PDF at given Q2 vs. x
         DO 5000 jq2=1,nq2pdf
            Q2 = Q2pdf(jq2)
            WRITE(3,5103) Q2 
 5103       FORMAT(//,' Q2=',F10.2)
            WRITE(3,5104)
 5104       FORMAT(/,2X,' X',10X,'xd',22X,'xu',22X,'xg',22X,'xub',
     2           21X,'xdb',21X,'xs',22X,'xc',22X,'xb')
            ALS0=DLOG(Q02/CURPAR(1)**2) 
            ALQ=DLOG(Q2/CURPAR(1)**2) 
            ALS=DLOG(ALQ/ALS0)
            do 5001 j=1,nx
               call pdf(0,xgrid(j),als,u0,d0,ub0,db0,sb0,cb0,bb0,glue0)
               pdf0(-5,j)=bb0
               pdf0(-4,j)=cb0
               pdf0(-3,j)=sb0
               pdf0(-2,j)=db0
               pdf0(-1,j)=ub0
               pdf0(-0,j)=glue0
               pdf0(1,j)=u0
               pdf0(2,j)=d0
               do 5002 k=1,nvar
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
 5002          continue
               do 5003 i=-5,2
                  err(i,j)=0.d0
                  do 5004 l=1,nvar
                     do 5004 m=1,nvar
                        err(i,j)=err(i,j)
     &                       +derv(i,j,l)*cov(l,m)*derv(i,j,m)
 5004                continue
                     err(i,j)=sqrt(err(i,j))
 5003             continue
 5001          continue
*          ... write to file
               DO 5021 j = 1, nx
                  WRITE(3,5031) xgrid(j)
     &                 ,(pdf0(ipt,j),err(ipt,j),ipt=2,-5,-1)
 5031             FORMAT(E10.3,16E12.4) 
 5021          CONTINUE
*          ...on to the next Q2 value...
 5000       continue
            
*    ... PDF at given x vs. Q2
            WRITE(3,*)
            WRITE(3,*) '-----------------------------------------------'
            DO jx=1,nxpdf
               x = xpdf(jx)
               WRITE(3,5203) x 
 5203          FORMAT(//,'  x=',E15.3)
               WRITE(3,5204)
 5204          FORMAT(/,2X,'Q2',10X,'xdb',21X,'xub',21X,'xg',22X,'xu',
     2              22X,'xd',22X,'xs',22X,'xc',22X,'xb')
               do j=1,nq2
                  Q2 = q2grid(j)
                  ALS0=DLOG(Q02/CURPAR(1)**2) 
                  ALQ=DLOG(Q2/CURPAR(1)**2) 
                  ALS=DLOG(ALQ/ALS0)
                  call pdf(0,x,als,u0,d0,ub0,db0,sb0,cb0,bb0,glue0)
                  pdf0(-5,j)=bb0
                  pdf0(-4,j)=cb0
                  pdf0(-3,j)=sb0
                  pdf0(-2,j)=db0
                  pdf0(-1,j)=ub0
                  pdf0(-0,j)=glue0
                  pdf0(1,j)=u0
                  pdf0(2,j)=d0
                  do k=1,nvar 
                     call pdf(k,x,als,u,d,ub,db,sb,cb,bb,glue)
                     kk=ifrep(k)
                     derv(-5,j,k)=(bb-bb0)/pwate(kk)
                     derv(-4,j,k)=(cb-cb0)/pwate(kk)
                     derv(-3,j,k)=(sb-sb0)/pwate(kk)
                     derv(-2,j,k)=(db-db0)/pwate(kk)
                     derv(-1,j,k)=(ub-ub0)/pwate(kk)
                     derv(0,j,k)=(glue-glue0)/pwate(kk)
                     derv(1,j,k)=(u-u0)/pwate(kk)
                     derv(2,j,k)=(d-d0)/pwate(kk)
                  end do        ! loop over k
                  do i=-5,2
                     err(i,j)=0.d0
                     do l=1,nvar
                        do m=1,nvar
                           err(i,j)=err(i,j)
     &                          +derv(i,j,l)*cov(l,m)*derv(i,j,m)
                        end  do
                     end do     ! loops over l,m
                     err(i,j)=sqrt(err(i,j))
                  end do        ! loop over i
               end do           ! loop over j
*     ... write to file
               DO  j = 1, nq2
                  WRITE(3,5031) q2grid(j)
     &                 ,(pdf0(ipt,j),err(ipt,j),ipt=2,-5,-1)
               end do
            end do              ! on to the next xB value...
         
         end if                 !END printing PDFs to file

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
c      WRITE(IUout,6001) RS,Y
c 6001 FORMAT(/' DRELL-YAN TEST AT SQRT(S)=',F6.2,' GEV AND Y=',F5.2)
c      WRITE(IUout,6004)
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
c      WRITE(IUout,6003) RTAU,ANS,xseclo,AKFAC,AKQ,AKG
c 6002 CONTINUE
c 6003 FORMAT(F10.2,5E15.4)

      call exit 
      STOP
      END 


c      subroutine pdfband(qsq,iband)
c      implicit real*8 (a-h,o-z)
c      character*30, output, file1
c      character*4, tmp
c      dimension derv(-5:2,30,35), err(-5:2,30), pdf0(-5:2,30), xgrid(30)
c      common/gint/gf(8,40,1080)
c      COMMON/MINIMB/CURPAR(60),UNCRT(60),ERMIN(40),PWATE(60),IFREP(40),
c     * NPARAM,NVAR,IREJ 
c      COMMON/Q2STUF/ Q02,Q2MAX
c      COMMON/INPUT/INS,NVL,NMAX,DELTA,S0,IORD
c      common/outputfile/output      
c      common/errors/cov(40,40)
c      DATA NX,XGRID/30,.0001,.0002,.0004,.0006,.0008,.001,.002,.004,
c     2.008,.016,.032,.064,.1,.15,.2,.25,.3,.35,.4,.45,.5,.55,.6,.65,
c     3.7,.75,.80,.85,.9,.95/
c      ALS0=DLOG(Q02/CURPAR(1)**2) 
c      ALQ=DLOG(QSQ/CURPAR(1)**2) 
c      ALS=DLOG(ALQ/ALS0)
c      call trmstr(output,len1)
c      tmp='.pdf'
c      call trmstr(tmp,len2)
c      len1=len1-4
c      file1=output(1:len1)//tmp(1:len2)
c      len=len1+len2
c      open(unit=3,file=file1(1:len),status='unknown')
c      do 5000 j=1,nx
c         call pdf(0,xgrid(j),als,u0,d0,ub0,db0,sb0,cb0,bb0,glue0)
c         pdf0(-5,j)=bb0
c         pdf0(-4,j)=cb0
c         pdf0(-3,j)=sb0
c         pdf0(-2,j)=db0
c         pdf0(-1,j)=ub0
c         pdf0(-0,j)=glue0
c         pdf0(1,j)=u0
c         pdf0(2,j)=d0
c         do 5001 k=1,nvar
c            call pdf(k,xgrid(j),als,u,d,ub,db,sb,cb,bb,glue)
c            kk=ifrep(k)
c            derv(-5,j,k)=(bb-bb0)/pwate(kk)
c            derv(-4,j,k)=(cb-cb0)/pwate(kk)
c            derv(-3,j,k)=(sb-sb0)/pwate(kk)
c            derv(-2,j,k)=(db-db0)/pwate(kk)
c            derv(-1,j,k)=(ub-ub0)/pwate(kk)
c            derv(0,j,k)=(glue-glue0)/pwate(kk)
c            derv(1,j,k)=(u-u0)/pwate(kk)
c            derv(2,j,k)=(d-d0)/pwate(kk)
c 5001    continue
c         do 5002 i=-5,2
c            err(i,j)=0.d0
c            do 5003 l=1,nvar
c            do 5003 m=1,nvar
c               err(i,j)=err(i,j)+derv(i,j,l)*cov(l,m)*derv(i,j,m)
c 5003       continue
c            err(i,j)=sqrt(err(i,j))
c 5002    continue
c 5000 continue
c      do 5004 i=-5,2
c         write(3,5010) i,qsq
c 5010    format(/ ,'parton flavor =',i5,' at Q**2 = ',f10.2,' GeV**2')
c         write(2,5011)
c 5011    format(/ ,5x,' x',5x,'pdf',5x,'error')
c         do 5005 j=1,nx
c            if (iband.eq.1)then
c               pdfout=pdf0(i,j)
c               errout=err(i,j)
c            else if(iband.eq.2)then
c               pdfout=1.d0
c               errout=err(i,j)/pdf0(i,j)
c            endif
c            write(3,5012) xgrid(j),pdfout,errout
c 5012       format(f10.4,2e12.4)
c 5005    continue
c 5004 continue
c      call exit
c      end




      subroutine fildat(nfit,npts,data,xerr,q2max)
      implicit real*8 (a-h,o-z)
      character*80 title(4)
      character*10 ititle(30),tmp
      character*200 file1,fln
      dimension data(6,2500),xerr(2500),ifit(30),q2cut(30),xlo(30)
     $     ,xhi(30),w2min(30),errn(30),dum(6),v(4)
      common/cuts/q2cut,xlo,xhi,w2min
      common/flags/itype(30),nflag(30),inorm(30),itgt(30),icorr(30)
     2     ,itmc(30),iht(30)
      common/expts/ititle,errn,iserr
      common/Onejetset/Jset 
      common/DeltaYmax/DelYmax
      common/cdf/covinv(33,33),cdfchisqr,ncdf,ncdfl,ncov

C     Try setting fln to the environment path of the data directory.
C     Assumes that the path for the data directory is set 
C     [Use one of the following commands in your terminal:
C        setenv cteqx_dat <path>
C        export cteqx_dat='<path>'  ]

      call GETENV('cteqx_dat',fln)
*      write(6,*),'data directory path = ',fln
      call trmstr(fln,len1)
      fln = fln(1:len1)//'/'
      len1 = len1+1

C     Read in the chosen data sets

      Jset=1
      DelYmax=7.0d0
      ncdf=0
      ncdfl=0
      ncov=0
      npts=0
      q2max=0.d0
      do 100 j=1,nfit
         tmp=ititle(j)
         call trmstr(tmp,len2)
         file1=fln(1:len1)//tmp(1:len2)
         len=len1+len2
         open(unit=8,file=file1(1:len),status='old')
***         print*,'* opening file=',file1(1:len),j,itype(j)
*       ... skip title lines
         do 99 l=1,3
            read(8,*) title(l)
 99      continue
*       ... read normalization error and number of crrelated errors 
c     read(8,*,err=98) errn(j),ncorr
c     
c     temporarily disabled ncorr because of format problems in some data files
c     
         read(8,*,err=98) errn(j)
***         print*, '* Norm. error =', errn(j)
c      goto 98
 97      ncorr=0.
 98      continue
*       ... reads variable names
         read(8,*) title(4)
 200     continue
*       ... reads data points until end of file
***         read(8,*,end=96) dum(1),dum(2),dum(3),dum(4)
***         print*, '* dum =', (dum(jk),jk=1,4)

         read(8,*,end=96) (dum(jk),jk=1,6),ak0,ak1

         if(iserr.eq.0)then
            err=dum(5)
         else if(icorr(j).eq.0)then
            err=sqrt(dum(5)**2+dum(6)**2)
         else if(icorr(j).ne.0)then
            err=dum(5)
         endif
         if(itype(j).lt.70)then
*          ... elm DIS data
            x=dum(1)
            q2=dum(2)
            dat=dum(3)
c$$$*          ... phenomenological Fe correction for itgt=-1,-2
c$$$            if(itgt(j).eq.-1)then
c$$$               dat=emcnmc(dat,x,q2)
c$$$               err=emcnmc(err,x,q2)
c$$$            else if(itgt(j).eq.-2)then
c$$$               dat=emcnmc(dat,x,q2)
c$$$               err=emcnmc(err,x,q2)
c$$$               dat=gomez(dat,x)
c$$$               err=gomez(err,x)
c$$$            endif
            if(x.lt.xlo(j))goto 200
            if(x.gt.xhi(j))goto 200
            if(q2.lt.q2cut(j))goto 200
            w2=q2*(1.d0/x-1.d0)+.88d0
            if(w2.lt.w2min(j))goto 200
            npts=npts+1
*          ... if reading JLab e00106 F2(D), divides F2 and error by 2
            if (itype(j).eq.8) then
               dat = dat/2
               err=err/2
            end if
            data(1,npts)=dat
            data(2,npts)=err
            data(3,npts)=x
            data(4,npts)=q2
            data(5,npts)=w2
            data(6,npts)=j
         else if(itype(j).eq.70.or.itype(j).eq.71)then
*          ...nutev iron (anti)neutrino cross section data
c
c            e=dum(1)
c            if(e.eq.85..or.e.eq.170.)then
c               continue
c            else
c               goto 200
c            endif
c            x=dum(2)
c            y=dum(3)         
c
c  modified to read ordered nutev data files
c
            x=dum(1)
            q2=dum(2)
            dat=dum(4)
            y=dum(3)
c            q2=2.*.938*e*x*y
*          ... phenomenological Fe correction for itgt=-1,-2
c$$$            if(itgt(j).eq.1)then
c$$$               dat=emcnmc(dat,x,q2)
c$$$               err=emcnmc(err,x,q2)
c$$$            else if(itgt(j).eq.2)then
c$$$               dat=emcnmc(dat,x,q2)
c$$$               err=emcnmc(err,x,q2)
c$$$               dat=gomez(dat,x)
c$$$               err=gomez(err,x)
c$$$            endif
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
*          ... DY data
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
c            if(itype(j).eq.106)y=dum(3)
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
*          ... JET data
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
C       ... photon data
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
 96      continue
         close(8)
 100  continue
      njet=0
      do 101 j=1,nfit
         if(itype(j).gt.200.and.itype(j).lt.300)njet=1
 101  continue

*     fake data point, used for output purposes
      data(1,npts+1)=0d0
      data(2,npts+1)=0d0
      data(3,npts+1)=0d0
      data(4,npts+1)=0d0
      data(5,npts+1)=0d0
      data(6,npts+1)=nfit+1

*    *** weights the errors according to a chi^2 weighting function
*    ... and saves the original data 
      do ipt = 1, npts
         do i=1,4
            v(i) = data(i+2,ipt)
         end do
         w = chi2wt(itype(data(6,ipt)),v)
         xerr(ipt) = data(2,ipt)
         data(2,ipt) = data(2,ipt)/dsqrt(w)
      end do
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
     2                 stop '# of jet data does not match Error Matrix.'
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
 999     Print*,'Cannot open file:',file1
         stop
C     *******************************       

      end



      subroutine rdxpos(xpos,nxpos,rspos,nrspos,npts,data,ndis,ndy,njet,
     2ngamjet)

      implicit real*8 (a-h,o-z)
      dimension data(6,2500),xpos(350),rspos(15)
      common/flags/itype(30),nflag(30),inorm(30),itgt(30),icorr(30)
     2     ,itmc(30),iht(30)
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
 90            continue
               nxpos=nxpos+1
               if(nxpos.gt.npts)return
               xpos(nxpos)=x
 91            continue
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
 80            continue
               nrspos=nrspos+1
               rspos(nrspos)=rs
 81            continue
            endif
      else if(it.gt.200.and.it.lt.300)then
         njet=njet+1
      else if(it.gt.300.and.it.lt.400)then
         ngamjet=ngamjet+1
         endif
 100  continue
c       print*,nxpos,xpos
      return
      end



