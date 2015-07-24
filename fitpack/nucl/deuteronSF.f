C **********************************************************************
      SUBROUTINE DEUTERONSF (ibeam,istruc,inuke,itm,iht,ndrv,v,xc,Fd)
C
C  Compute 0.5* deuteron F2 structure function (F2d) per nucleon 
C  at x (=v(1)) and Q^2 (=v(2)) by smearing nucleon structure function 
C  (array F2Narr dimension nx defined over x array xarr) with nucleon 
C  momentum distribution.
C
C    ibeam             = (i) 1=e 2=nu 3=nubar
C    istruc            = (i)  0=FL  1=F1  2=F2  3=x*F3
C                            -1=nPDFs as smeared PDFs; results returned in 
C                               common/nPDF/UA,DA,UBA,DBA,SBA,CBA,BBA,GLUEA
C                               and Fd set to -9999999d0
C    inuke             = (i) nuclear corrections: ABCDE
C                            A: TMC in smearing functions: 
C                                   gamma^2=1+4*xN^2*Lambda^2/Q^2
C                               0: Lambda^2=M
C                               1: Lambda^2=0  ("Bjorken limit", no TMC)
C                            B: Shadowing correction
C                               0 = no
C                               1 = Melnitchouk-Thomas 
C                            C: Off-shell correction
C                               0=none
C                               1=KP phenomenological fit
C                               2=MST 
C                               3=KP fit parametrized
C                               4=mKP lo
C                               5=mKP central
C                               6=mKP hi
C                            D: Deuteron wave-function
C                               0=Paris
C                               1=AV18
C                               2=CDBonn
C                               3=WJC1
C                               4=WJC2
C                            E: Deuteron correction model
C                               0=isospin average
C                               1=density model
C                               2= --- not supported anymore ----
C                               3=nuclear smearing WBA
C                               4=nuclear smearing WBAREL (WBA relativized)
C                               5=nuclear smearing AQV  (M/p0 norm)
C                               6=nuclear smearing AQVc (const. norm)
C    itm,iht,ndrv,v,xc = see 'strfn' routine
C    Fd                = deuteron structure function per nucleon
C
C  This version compues F2(xN) as needed iside the convolution
C  integral, instead of precomputing it and interpolating an array.
C  (for the latter, see deuteronSF_pre.f)
C
*     EXAMPLES:
*
*     + free, on-shell,              inuke = 0         (00000 ==  0)
*     + free, MST off-sh., Bj.limit  inuke = 10200       
*     + km on-shell Paris            inuke = 2         (00002 ==  2)
*     + wba AV18, mKP centr, shadow. inuke = 1513      (01513 ==  1513)
*     + dmc on-shell                 inuke = 1         (00001 ==  1)
*
C **********************************************************************
      IMPLICIT NONE
      INTEGER	ibeam,inuke,istruc,itm,iht,ndrv,idum1,idum2,idum3,idum4
      double precision v(4),xc(100),Fd
      INTEGER nx,ix, ismear
      PARAMETER (nx=100)
      REAL*8	x,Q2,varr(4),xarr(nx),FNarr(nx),dsf,dmc

*    *** Functions
      double precision gomez

*    *** common blocks

*     Proton, neutron fraction    
      double precision az,an
      common/target/az,an
      
      az=0.5D0                  ! p + n isoscalar, i.e., (p+n)/2
      an=0.5D0

      call split_nuke(inuke,idum1,idum2,idum3,idum4,ismear) ! selects smearing model
	
      if (ismear.le.1) then
*    *** p+n isoscalar 
         CALL strfn (ibeam,istruc,itm,iht,ndrv,v,xc,Fd)
         if (ismear.eq.1) then
*       *** Density model corrections (for backward compatibility)
            dmc=gomez(1d0,v(1))
            Fd = Fd/dmc
         end if
      else if ((ismear.ge.2).and.(ismear.le.9)) then
*       *** Nuclear smearing model
C       ...Compute nucleon str.fn.to be (interpolated and) smeared
*         varr(2) = v(2)	! store Q^2 value in new array
*	  DO ix=1,nx
*	     x = DFLOAT(ix)/100.D0 - 1d-4
*	     xarr(ix) = x	! x array for interpolation 
*	     varr(1) = x	! value of x at which SF to be calculated here
*	     CALL strfn (ibeam,istruc,itm,iht,ndrv,varr,xc,FNarr(ix))
*	  ENDDO
C        ...Compute deuteron str.fn. by smearing nucleon SF (F2 only for now)
         !if (istruc.ne.2) then
         !   write(*,*) 'ERROR (deuteronsf): F1 and F3 not implemented'
         !   stop
         !end if
         CALL SMEARF2 (ibeam,istruc,inuke,itm,iht,ndrv,v,xc,Fd)
      else
         print*, 'ERROR(DEUTERONSF): inuke out of range =',inuke
      end if
	
      RETURN
      END


C **********************************************************************
      SUBROUTINE SMEARF2 (ibeam,istruc,inuke,itm,iht,ndrv,v,xc,FFd)
C
C  Smear nucleon structure function (array FNarr dimension nx defined
C  over x array xarr) with nucleon light-cone momentum distribution.
C
C  xarr is free-nucleon x in [0,1]
C  Convolution set up for "deuteron" xD in [0,1], rather than
C    "nucleon" x in [0,2]
C    => need factor 2 conversion xN = 2 xD [change later].
C
C  If istruc=-1, smears the PDFs at gamma=1, results are returned
C  in common/nPDF/
C
C  NOTE: az,an must have been previously set in common/target/az,an
C
C **********************************************************************
      IMPLICIT NONE

      INTEGER	ibeam,inuke,istruc,itm,iht,ndrv
      double precision v(4),xc(100),FFd,weight
      
      INTEGER i,ibj,ioff,ismear,iwfn,ishad,ilam
      REAL*8  xD,x,yD,yDmin,yDmax,fy_diag,fy_off,FFN,F2N,err,Q2,gamma
      REAL*8  PHI_INT2D,PHI_INT3D
      REAL*8  OFF_MST,OFF_KP,OFF_KPMDL,DEL_SHAD,shad,offsh
      double precision S0,S,U,D,UB,DB,SB,CB,BB,GLUE

*    *** common blocks

*     Gaussian integration
*      double precision xi(16),WI(16),xx(17)
*      double precision xi(32),WI(32),xx(33)
      double precision xi(96),WI(96),xx(97)
      integer nterms
*      COMMON/GAUS16/XI,WI,NTERMS,XX
*      COMMON/GAUS32/XI,WI,NTERMS,XX
      COMMON/GAUS96/XI,WI,NTERMS,XX

*     Nuclear PDFS
      double precision UA,DA,UBA,DBA,SBA,CBA,BBA,GLUEA
      common/nPDF/UA,DA,UBA,DBA,SBA,CBA,BBA,GLUEA

*     Q^2 min and max
      double precision Q02,Q2MAX
      COMMON/Q2STUF/ Q02,Q2MAX

C    *** Constants

C    ... hc [MeV.fm]           - conversion factor
C    ... mN [MeV] mNGeV [GeV]  - nucleon mass 
C    ... MD [MeV]              - deuterium mass, including binding energy
      REAL*8  pi,hc,mN,mNGeV,MD 
      parameter(pi=3.141592653589793d0,hc=197.327D0,mN=938.91897D0
     $     ,mNGeV=mN/1d3,mD=2*mN-2.224575D0 )

c...Smearing model 

      ! selects Bj limit, off-shell and smearing model
      call split_nuke(inuke,ilam,ishad,ioff,iwfn,ismear)  
      
C...Value of x,gamma at which convolution to be made
      x = v(1)	
      xD = x/2                  ! deuteron Bjorken variable
      Q2 = v(2)

      if (ilam.ge.1) then
*    ... Bjorken limit
         gamma = 1d0
      else
*    ... finite Q^2
         gamma = DSQRT(1.D0 + 4*x**2 * mNGeV**2/Q2)
      end if
	
C...Convolution approximation (x used in convolution here is xD in [0,1])
      yDmax = 1d0 
      yDmin = xD	
      FFD = 0d0
      if(istruc.eq.-1) then
         UA = 0
         DA = 0
         UBA = 0
         DBA = 0
         SBA = 0
         CBA = 0
         BBA = 0
         GLUEA  = 0
      end if

      DO I=1,NTERMS
         yD=0.5*(yDmax-yDmin)*XI(I)+0.5*(yDmax+yDmin) ! y = yD in [x,1]
         v(1)=xD/yD		! computes FN at xN=xD/yD
         if ((ioff.eq.0).or.(ioff.ge.2)) then 
            ! 2D interpolate smearing fn. in y and gamma 
            if (istruc.eq.0) then 
               !!!! FL = f0 * FLN + f1 * F2N
               fy_diag = PHI_INT2D (yD,gamma,0,ismear,iwfn) 
               CALL strfn (ibeam,0,itm,iht,ndrv,v,xc,FFN)
               fy_off = PHI_INT2D (yD,gamma,1,ismear,iwfn) 
               CALL strfn (ibeam,2,itm,iht,ndrv,v,xc,F2N)
            else if (istruc.eq.1) then 
               !!!! (xF1) = (f0 * (xN*F1N) + (1/4)*f1 * F2N
               fy_diag = PHI_INT2D (yD,gamma,0,ismear,iwfn) 
               CALL strfn (ibeam,1,itm,iht,ndrv,v,xc,FFN)
               FFN = v(1)*FFN
               fy_off = PHI_INT2D (yD,gamma,1,ismear,iwfn) / 4d0
               CALL strfn (ibeam,2,itm,iht,ndrv,v,xc,F2N)
            else if (istruc.eq.2) then 
               !!! F2 = f2 * F2N
               fy_diag = PHI_INT2D (yD,gamma,2,ismear,iwfn)
               CALL strfn (ibeam,2,itm,iht,ndrv,v,xc,FFN)
               fy_off  = 0d0
               F2N = 0d0
            else if (istruc.eq.3) then 
               !!!! xF3 = f0 * xF3N
               print*, 'ERROR(SMEARF2): F3 convolutin formula yet to be'
     &              //'checked analytically'
               stop
               !fy_diag = PHI_INT2D (yD,gamma,0,ismear,iwfn)
               !CALL strfn (ibeam,3,itm,iht,ndrv,v,xc,FFN)
               !fy_off  = 0d0
               !F2N = 0d0
            else if (istruc.eq.-1) then
               !!! Nuclear PDFs as convolution of nucleon PDFs
               !!! Beware: it makes physical sense only at very large Q^2
               fy_diag = PHI_INT2D (yD,1d0,0,ismear,iwfn)  ! gamma=1
               !!! the following 3 lines are for compatibility with fitpack...
               S0=DLOG(Q02/XC(1)**2)
               S=DLOG(DLOG(Q2/XC(1)**2)/S0)
               call fsupdf(0,x,s,u,d,ub,db,sb,cb,bb,glue)
               !!! if working with stand-alone package, could use this 1-line:
               !!! call PDFsa(x,Q2,U,D,UB,DB,SB,CB,BB,GLUE) 
            else
               fy_diag = 0d0
               FFN = 0d0
               fy_off = 0d0
               F2N = 0d0
            end if
         else if (ioff.eq.1) then            ! 3D interpolate smearing fn
            print*, 'WARNING!'
            print*, 'KP off-shell smearing has become OBSOLETE!!! '
            print*, '(Needs to be rewritten as 2D interpolation '
            print*, 'using the p2-m2 averaged smearing functions)'
            print*, 'STOPPING HERE!'
            stop
*               fy = PHI_INT3D (yD,gamma,xD/yD,ismear,ioff)  
*                                                ! in y, gamma, xN=xB/y=xD/yD
         else 
            write(*,*) 'ERROR (smearf2): ioff out of range = ',ioff
            stop
         end if

         if (v(1).lt.0.995d0) then
            weight = .5*(yDmax-yDmin)*WI(I)
            if (istruc.ne.-1) then ! Structure function convolution
               FFd = FFd + weight * (fy_diag*FFN + fy_off*F2N)
            else                   ! PDF convolution
               UA    = UA    + weight * fy_diag*U
               DA    = DA    + weight * fy_diag*D
               UBA   = UBA   + weight * fy_diag*UB
               DBA   = DBA   + weight * fy_diag*DB
               SBA   = SBA   + weight * fy_diag*SB
               CBA   = CBA   + weight * fy_diag*CB
               BBA   = BBA   + weight * fy_diag*BB
               GLUEA = GLUEA + weight * fy_diag*GLUE
            end if
         end if
         
         !print*, '* fy_diag _off =',fy_diag, fy_off
         !print*, '* FFN, F2N =', FFN,F2N
         !print*, '* FFd =',FFd
         !stop

      ENDDO
      if (istruc.eq.1) then       ! convolution done for x*F1_A, returns F1_A
         FFd = FFd/x
      else if (istruc.eq.-1) then ! nPDFs are output in common/nPDF/  
         FFd = -9999999d0
      end if


*    *** Parametrized off-shell corrections (F2 only!!)
      offsh=0d0
      if (ioff.eq.2) then
         ! Off-shell corrections by Melnitchouk-Schreiber-Thomas
         ! (parametrization by W.Melnitchouk) 
         offsh = OFF_MST(x)
      else if (ioff.eq.3) then
         ! approximate Kulagin-Petti NPA(2007) 
         ! (parametrization by W.Melnitchouk) 
         offsh = OFF_KP(x)
      else if (ioff.eq.4) then
         ! modified Kulagin-Petti (by W.Melnitchouck) 
         ! lowest correction  
         offsh = OFF_KPMDL(x,-1)
      else if (ioff.eq.5) then
         ! modified Kulagin-Petti (by W.Melnitchouck) 
         ! central correction  
         offsh = OFF_KPMDL(x,0)
      else if (ioff.eq.6) then
         ! modified Kulagin-Petti (by W.Melnitchouck) 
         ! highest correction  
         offsh = OFF_KPMDL(x,1)
      end if
      FFd = FFd/(1.D0-offsh)

*    *** Parametrized deuterium shadowing
      shad = 0d0
      if (ishad.eq.1) then
         ! Melnitchouk and Thomas, PRD 47, 3783 (1993)
         shad = DEL_SHAD(x,Q2)         
      end if
      FFd = FFd/(1.d0+shad)


      v(1) = x
      
      RETURN
      END


************************************************
*     divides nuclear correction flag into its components
      subroutine split_nuke(inuke,ilam,ishad,ioff,iwfn,ismear) 
*     Author: A.Accardi
*     date: 18 Feb 2010
*
*     INPUT:
*
*       inuke   = (i) nuclear correction: ABCDE
C                     see 'deuteronsf' routine
*
*     OUTPUT:
*
*       ilam    = (i) A
*       ishad   = (i) B
*       ioff    = (i) C
*       iwfn    = (i) D
*       ismear  = (i) E
*

      implicit none

      integer inuke,ilam,ishad,ioff,iwfn,ismear

      ilam   = inuke/10000               ! treatment of TMC scale Lambda 
      ishad  = (inuke-ilam*10000)/1000   ! shadowing corrections
      ioff   = (inuke-ilam*10000         ! off-shell corrections
     &               -ishad*1000)/100     
      iwfn   = (inuke-ilam*10000         ! nuclear wave function
     &               -ishad*1000
     &               -ioff*100)/10     
      ismear = (inuke-ilam*10000         ! convolution formula
     &               -ishad*1000         ! or deuteron correction model 
     &               -ioff*100
     &               -iwfn*10)     

      return
      end



************************************************
*     Density model correction gomez(f) = N/D
      function gomez(f,x)
      implicit double precision (a-h,o-z)
      data p1,p2,p3,p4,p5,p6,p7/1.0164d0,-.047762d0,-.13354d0,
     2.35303d0,.22719d0,-1.2906d0,5.6075d0/
c
c  Fit to the data for F2D/F2N where N=(p+n)/2
c  as extractred by Gomez et al PRD 49, 4348 (1994)
c  [Original reference by Frankfurt and Strikman]
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


C ***********************************************************************
      FUNCTION OFF_MST (x)
C
C  Nucleon off-shell correction in MST model.
C
C  Defined such that F2d = F2d(conv) + del^off F2d
C       with OFF_MST = del^off F2d / F2d
C                    = off-shell correction w.r.t. F2d
C
C  Parameters defined for x = nucleon scaling variable in [0,2].
C
C  Ref: Melnitchouk, Schreiber, Thomas, PRD 49, 1183 (1994);
C       PLB 335, 11 (1994)
C ***********************************************************************
      IMPLICIT NONE
      REAL*8  OFF_MST,x
      REAL*8  aoff(0:5)
      DATA    aoff /-0.014D0, 3.D0, 20.D0, 1.067D0, 1.5D0, 18.D0/

      OFF_MST = 0.D0
      IF (x.LE.0.16D0) RETURN
      
      OFF_MST = aoff(0) * (1.D0 + aoff(1)*x**aoff(2))
     &     * (1.D0 - (aoff(3) - x**aoff(4))**aoff(5))
      
      RETURN
      END



C ***********************************************************************
      FUNCTION OFF_KP (x)
C
C  Analytic parametrization of nucleon off-shell correction from
C    Kulagin-Petti fit
C
C  Defined such that F2d = F2d(conv) + del^off F2d
C       with OFF_KO = del^off F2d / F2d
C                   = off-shell correction w.r.t. F2d
C
C  Parameters defined for x = nucleon scaling variable in [0,2]. 
C
C  Ref: Kulagin, Petti, NPA 765, 126 (2006).
C
C  July 2010.
C ***********************************************************************
      IMPLICIT NONE
      REAL*8  OFF_KP,x
      REAL*8  aoff(0:4),boff(0:13),coff(0:5)  
      DATA    aoff /-0.010675D0, 0.27315D0, -0.96047D0, 1.2396D0,
     &     -0.71897D0/
      DATA    boff /-0.11185D-1, 0.29485D0, -1.3024D0, 4.1337D0,
     &     -17.047D0, 66.378D0, -184.54D0, 293.49D0,
     &     -64.581D0, -785.17D0, 1779.8D0, -1908.6D0,
     &     1066.9D0, -250.54D0/
      DATA    coff /3.9913D0, 4.3786D0, 2.4799D0, 2.5043D0,
     &     -3.9996D0, 0.018932D0/
      
      OFF_KP = 0.D0
      
C...13th order polynomial fit valid for x < 0.91
c       OFF_KP = boff(0) + boff(1)*x + boff(2)*x**2 + boff(3)*x**3
c     &        + boff(4)*x**4 + boff(5)*x**5 + boff(6)*x**6
c     &        + boff(7)*x**7 + boff(8)*x**8 + boff(9)*x**9
c     &        + boff(10)*x**10 + boff(11)*x**11 + boff(12)*x**12
c     &        + boff(13)*x**13
c
c       IF (x.GT.0.85D0) RETURN
C...4th order polynomial fit valid for x < 0.86
c       OFF_KP = aoff(0) + aoff(1)*x + aoff(2)*x**2 + aoff(3)*x**3
c     &        + aoff(4)*x**4
        
C...6-parameter fit  (courtesy of Simona Malace)
      OFF_KP = coff(0) + coff(1)*x + coff(2)*x**coff(3)
     &     + coff(4)*DEXP(x) + coff(5)/DLOG(x)

      RETURN
      END


C ***********************************************************************
      FUNCTION OFF_KPMDL (x,imode)
C       
C  Analytic parametrization of nucleon off-shell correction from
C    Kulagin-Petti model, modified for realistic deuteron parameters.
C  
C  Defined such that F2d = F2d(conv) + del^off F2d
C       with OFF_KPMDL = del^off F2d / F2d  
C                      = off-shell correction w.r.t. F2d
C
C     imode = -1,0,1 --> lowest, central, highest correction
C
C  Parameters defined for x = nucleon scaling variable in [0,2].
C
C  Ref: Kulagin, Petti, NPA 765, 126 (2006).
C
C  August 2010
C ***********************************************************************
      IMPLICIT NONE
      REAL*8  OFF_KPMDL,x
      integer imode
      REAL*8  poff(10),pmin(10),pmax(10)
      DATA    poff /0.00917324D0, -4.87488D0, -0.129925D0,
     &     24.1748D0, 0.249259D0, -4.24231D0,
     &     -0.0353244D0, -59.7096D0, -0.214936D0, 46.3766D0/
      DATA    pmin /0.00606295D0, -3.97072D0, -0.0551534D0,
     &     10.5675D0, -0.245408D0, -5.09068D0,
     &     0.699754D0, -18.0466D0, -0.426636D0, 24.6164D0/
      DATA    pmax /0.00964982D0, -6.84144D0, -0.0985649D0,
     &     43.0597D0, -0.0181488D0, -23.9927D0,
     &     0.616283D0, -75.3024D0, -0.711494D0, 65.4527D0/
      
      OFF_KPMDL = 0.D0

      if (imode.eq.-1) then 
         ! Lower limit (least negative)
         OFF_KPMDL = ( pmin(1) + pmin(3)*x + pmin(5)*x**2
     &        + pmin(7)*x**3 + pmin(9)*x**4 )
     &        * ( 1.D0    + pmin(2)*x + pmin(4)*x**2
     &        + pmin(6)*x**3 + pmin(8)*x**4 + pmin(10)*x**5 )
      else if (imode.eq.0) then
         ! Central values (fits courtesy of Simona Malace)
         OFF_KPMDL = ( poff(1) + poff(3)*x + poff(5)*x**2
     &        + poff(7)*x**3 + poff(9)*x**4 )
     &        * ( 1.D0    + poff(2)*x + poff(4)*x**2
     &        + poff(6)*x**3 + poff(8)*x**4 + poff(10)*x**5 )
      else if (imode.eq.1) then
         ! Upper limit (most negative)
         OFF_KPMDL = ( pmax(1) + pmax(3)*x + pmax(5)*x**2 
     &        + pmax(7)*x**3 + pmax(9)*x**4 )
     &        * ( 1.D0    + pmax(2)*x + pmax(4)*x**2
     &        + pmax(6)*x**3 + pmax(8)*x**4 + pmax(10)*x**5 )
      end if

      RETURN
      END


C ***********************************************************************
        FUNCTION DEL_SHAD (x,Q2)
C
C  Nuclear shadowing correction to F2d, including VMD, Pomeron and
C    meson (antishadowing) exchange.
C
C  Defined such that F2d = F2d(noshad) - \delta(shad) F2d
C
C  Parameters defined for x = nucleon(!) scaling variable in [0,2].
C          
C  Ref: Melnitchouk and Thomas, PRD 47, 3783 (1993)
C ***********************************************************************
        IMPLICIT NONE
        REAL*8  DEL_SHAD,x,Q2
        REAL*8  Q02,DEL_V,DEL_P,DEL_M
        REAL*8  A(0:4),B(0:4),C(0:4)
        DATA    A /-0.038D0, -0.04D0, 1.8D0, -3.4D0, 0.9D0/
        DATA    B /-0.003D0, -0.13D0, 5.D0,  -2.2D0, 0.4D0/
        DATA    C / 0.002D0,  0.03D0, 6.D0,   0.D0,  0.D0 /

        DEL_SHAD = 0.D0
        IF (x.GT.0.3D0) RETURN          ! param. valid for x <~ 0.3
   
C...VMD (Q2 in GeV^2, Q02 = 0.7 GeV^2)
        Q02 = 0.7D0
        DEL_V = Q2/(1 + Q2/Q02)**2
     &        * A(0) * x**A(1) * (1.D0-x)**A(2) * (1.D0 + A(3)*x**A(4))

C...Pomeron-exchange
        DEL_P = B(0) * x**B(1) * (1.D0-x)**B(2) * (1.D0 + B(3)*x**B(4))
      
C...Meson-exchange (antishadowing)
        DEL_M = C(0) * x**C(1) * (1.D0-x)**C(2) * (1.D0 + C(3)*x**C(4))
      
        DEL_SHAD = DEL_V + DEL_P + DEL_M
      
        RETURN
        END



**********************************************

      subroutine smearfile(outfile,outlen,ismear,iwfn,ioff)

      implicit none

      character outfile*150

      integer len,oldlen,outlen,ismear,iwfn,ioff


      len = 4
      outfile = 'phi.'
      oldlen = len+1

c$$$      if (ismear.eq.2) then
c$$$         len = oldlen + 1
c$$$         outfile(oldlen:len) = 'km'
c$$$         oldlen = len + 1
c$$$      else 
      if (ismear.eq.3) then
         len = oldlen + 2
         outfile(oldlen:len) = 'wba'
         oldlen = len + 1
      else if (ismear.eq.4) then
         len = oldlen + 5
         outfile(oldlen:len) = 'wbarel'
         oldlen = len + 1
      else if (ismear.eq.5) then
         len = oldlen + 5
         outfile(oldlen:len) = 'aqv'
         oldlen = len + 1
      else if (ismear.eq.6) then
         len = oldlen + 5
         outfile(oldlen:len) = 'aqvc'
         oldlen = len + 1
      else if ((ismear.ne.0).and.(ismear.ne.1)) then
         len = oldlen + 7
         outfile(oldlen:len) = '_smear??'
         oldlen = len + 1
      end if

      if (iwfn.eq.0) then
         len = oldlen + 5
         outfile(oldlen:len) = '_paris'
         oldlen = len + 1
      else if (iwfn.eq.1) then
         len = oldlen + 4
         outfile(oldlen:len) = '_AV18'
         oldlen = len + 1
      else if (iwfn.eq.2) then
         len = oldlen + 6
         outfile(oldlen:len) = '_CDBONN'
         oldlen = len + 1
      else if (iwfn.eq.3) then
         len = oldlen + 4
         outfile(oldlen:len) = '_WJC1'
         oldlen = len + 1
      else if (iwfn.eq.4) then
         len = oldlen + 4
         outfile(oldlen:len) = '_WJC2'
         oldlen = len + 1
      else 
         len = oldlen + 5
         outfile(oldlen:len) = '_wfn??'
         oldlen = len + 1
      end if

      if (ioff.eq.1) then
         len = oldlen + 4
         outfile(oldlen:len) = '_KPoff'
         oldlen = len + 1
      else if (ioff.ne.0) then
         len = oldlen + 7
         outfile(oldlen:len) = '_offsh??'
         oldlen = len + 1
      end if

      outlen = len

      return
      end





