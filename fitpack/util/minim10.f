      SUBROUTINE MINIM(NX,NNPTS,DATA,IDRV) 
      implicit real*8(a-h,o-z)
      CHARACTER*10 IPNAME(100)
      character*4 fnot,froz,blnk,tol,err
      REAL*8 MODALF 
      DIMENSION ALPHA(50,50),MODALF(50,50),BETA(50),COV(50,50) 
      DIMENSION OLDPAR(100),DATA(1,1)
      COMMON/ERRORS/COV
      COMMON/AMAT/ALPHA,BETA
      COMMON/MINIMC/IPNAME,PLB(100),PUB(100),TOLRNC,IPAR(100),NSTEPS,
     1IPRNT,IERR,IBND,CHISQR 
      COMMON/MINIMB/CURPAR(100),UNCRT(100),ERMIN(50),PWATE(100),
     1IFREP(50),NPARAM,NVAR,IREJ 
      DATA FROZ/'FIX '/,BLNK/'    '/,FNOT/' NOT'/ 
C
C  MINIMIZATION ROUTINE USING THE MARQUARDT ALGORITHM
C  FROM BEVINGTON, PAGE 235. SEE PROGRAM CURFIT.
C  VERSION CREATED BY D.W. DUKE CIRCA 1978.
C  COMMENTS AND SOME MODS ADDED BY J.F. OWENS
C
C INPUT AND OUTPUT UNIT NUMBERS
      MI=5
      MO=4
      FLAMDA=0.001 
C COUNT NUMBER OF PARAMETERS VARIED
      NVAR=0 
      DO 340 I=1,NPARAM 
      IF(PWATE(I).EQ.0.0) GO TO 340 
      NVAR=NVAR+1 
      IFREP(NVAR)=I 
C IFREP IS THE SEQUENCE OF PARAMETER NUMBERS FOR THOSE BEING VARIED
  340 CONTINUE 
C CALCULATE CHI SQUARE
      CALL DATSCN(1,CHISQR,NX,NNPTS,DATA)
c      write(6,*), ierr,iprnt
      IF(IERR.EQ.1) GO TO 1 
      IF(IPRNT.EQ.0) GO TO 1
      WRITE(MO,1111)(IP,IPNAME(IP),CURPAR(IP),PWATE(IP),IP=1,NPARAM)
 1111 FORMAT(5X,I5,2X,A10,'  VALUE=',G14.5,'  WEIGHT=',G14.5)
      WRITE(MO,1305) TOLRNC,NSTEPS 
      WRITE(MO,1300) CHISQR 
      WRITE(MO,1310)
    1 DO 100 NOSTEP=1,NSTEPS
      IF(FLAMDA.GT.0.001) FLAMDA=0.001 
      IF(FLAMDA.LT.1.E-10) FLAMDA=1.E-10
C SAVE OLD PARAMETER VALUES
      OLDCHI=CHISQR 
      DO 5 I=1,NPARAM 
    5 OLDPAR(I)=CURPAR(I) 
C EVALUATE ALPHA AND BETA
      CALL DATSCN(IDRV,DUM,NX,NNPTS,DATA) 
      NTRIES=0 
   10 DO 20 I=1,NVAR
      DO 15 J=1,NVAR
c      WRITE(6,2001) I,J
c 2001 FORMAT(2I10)
   15 MODALF(I,J)=ALPHA(I,J)/ SQRT(ALPHA(I,I)*ALPHA(J,J)) 
   20 MODALF(I,I)=1.+FLAMDA
C INVERT MODALF 
      CALL MINV(MODALF,NVAR,DET)
      NTRIES=NTRIES+1 
      NISTEP=NOSTEP
C UPDATE PARAMETERS 
      DO 25 I=1,NVAR
      K=IFREP(I) 
      DO 26 J=1,NVAR      
      CURPAR(K)=CURPAR(K)+MODALF(I,J)*BETA(J)/ SQRT(ALPHA(I,I)* 
     1ALPHA(J,J))
 26   continue
 25   continue
   32 CALL DATSCN(1,CHISQR,NX,NNPTS,DATA)
      IF(CHISQR.LT.OLDCHI) GO TO 40 
      IF(NTRIES.GT.6) GO TO 110 
      FLAMDA=10.*FLAMDA 
      DO 35 I=1,NPARAM 
   35 CURPAR(I)=OLDPAR(I) 
      GO TO 10 
C CHECK FOR CONVERGENCE
   40 IF((OLDCHI-CHISQR).LT.TOLRNC) GO TO 110 
C CHECK NUMBER OF STEPS
      IF(NOSTEP.EQ.NSTEPS) GO TO 110 
      IF(IPRNT.EQ.0) GO TO 100
      ERR=BLNK 
      IF(PWATE(1).EQ.0.0) ERR=FROZ
      IF(IPRNT.NE.0) WRITE(MO,1320)NISTEP,CHISQR
      DO 60 I=1,NPARAM 
      ERR=BLNK 
      IF(PWATE(I).EQ.0.0) ERR=FROZ 
      IF(IPRNT.NE.0) WRITE(MO,1330) IPAR(I),IPNAME(I),CURPAR(I),ERR 
   60 CONTINUE 
  100 FLAMDA=FLAMDA/10.0 
  110 IF(IERR.EQ.1) RETURN
C EVALUATE ERRORS AND PRINT OUT RESULTS OF FIT 
      CALL DATSCN(IDRV,DUM,NX,NNPTS,DATA) 
      DO 120 I=1,NVAR 
      DO 120 J=1,I 
      MODALF(I,J)=ALPHA(I,J)/ SQRT(ALPHA(I,I)*ALPHA(J,J)) 
  120 MODALF(J,I)=MODALF(I,J) 
      CALL MINV(MODALF,NVAR,DET) 
      DO 130 I=1,NVAR 
      AMOD=MODALF(I,I)/ALPHA(I,I) 
      IMOD=1 
      IF(AMOD.LT.0.0) IMOD=-1 
  130 UNCRT(IFREP(I))=IMOD*SQRT(IMOD*AMOD) 
      CHIPNT=CHISQR/(NNPTS-IREJ-NVAR) 
      IF(IPRNT.EQ.0) GO TO 170
      WRITE(MO,1320)NISTEP,CHISQR 
      WRITE(MO,1335)(IPAR(I),IPNAME(I),CURPAR(I),UNCRT(I),I=1,NPARAM) 
      WRITE(MO,1340) NISTEP 
      TOL=BLNK 
      IF(dabs(OLDCHI-CHISQR).GT.TOLRNC) TOL=FNOT 
      WRITE(MO,1350) TOL 
      IF(IREJ.NE.0) WRITE(MO,1360) IREJ 
  160 WRITE(MO,1345)CHIPNT 
      DO 180 I=1,NVAR 
      DO 180 J=1,NVAR 
  180 MODALF(I,J)=MODALF(I,J)/SQRT(ALPHA(I,I)*ALPHA(J,J)) 
      DO 182 I=1,NVAR 
      DO 181 J=1,NVAR 
      COV(I,J)=MODALF(I,J)
  181 ALPHA(I,J)=MODALF(I,J)/SQRT(MODALF(I,I)*MODALF(J,J))
  182 WRITE(MO,1380)(ALPHA(I,J),J=1,NVAR) 
      WRITE(MO,1370)
  170 CONTINUE
      RETURN
 1300 FORMAT(28H0THE INITIAL CHI SQUARED IS ,1PE11.4) 
 1305 FORMAT(1H0,17HTHE TOLERANCE IS ,F9.4/32H THE MAXIMUM NUMBER OF STEPS IS ,I3) 
 1310 FORMAT(33H0  STEP  CHI SQUARED    PARAMETER,6X,5HVALUE,9X,5HERROR) 
 1320 FORMAT(1H ,60(1H-)/3H   ,I3,3X,1PE11.4) 
 1330 FORMAT(22X,I2,2X,A8,2X,1PE11.4,7X,A4) 
 1335 FORMAT(22X,I2,2X,A8,2X,1PE11.4,3X,1PE11.4) 
 1340 FORMAT(35H0 MINIMIZATION TERMINATED WITH STEP,I3) 
 1345 FORMAT(39H0 CHI SQUARED PER DEGREE OF FREEDOM IS ,1PE10.3) 
 1350 FORMAT(22H  TOLERANCE CRITERION ,A4,10H SATISFIED)
 1360 FORMAT(12H0 THERE WERE,I3,40H DATA POINTS WITH ZERO OR NEGATIVE ERROR)
 1370 FORMAT(/1H0,16X,32H***   ***   ***  ***   ***   ***) 
 1380 FORMAT(1X,20F7.3) 
      END 


      SUBROUTINE DATSCN(MODE,CHISQR,NX,NNPTS,DATA)
      implicit real*8 (a-h,o-z)
      character*10 ititle(35)
      dimension errn(35),iserr(35),itype(35),nflag(35),inorm(35),
     2itgt(35),icorr(35)
      DIMENSION ALPHA(50,50),BETA(50) 
      COMMON/MINIMB/CURPAR(100),UNCRT(100),ERMIN(50),PWATE(100),
     1IFREP(50),NPARAM,NVAR,IREJ 
      COMMON/AMAT/ALPHA,BETA
      common/expts/ititle,errn,iserr
      common/flags/itype,nflag,inorm,itgt,icorr
      common/normalization/nfit 
      common/cdf/covinv(33,33),cdfchisqr,ncdf,ncdfl,ncov
      DIMENSION DATA(NX,NNPTS),V(10),DRV(100),cdfxsec(100),
     2cdfdrv(100,100)
      save cdfxsec
      NN=NX-2 
C CALCULATE CHI SQUARE
      IF(MODE.NE.1) GO TO 100 
      CHISQR=0.0 
      DO 10 NPT=1,NNPTS 
C FILL ARRAY V FOR DATA POINT LABELLED BY NPT
         DO 5 I=1,NN 
            V(I)=DATA(I+2,NPT)
 5       continue
      IF(DATA(2,NPT).LE.0.d0) GO TO 10
      index=v(4)      
      if(itype(index).eq.201.and.icorr(index).ne.0)goto 10 
      chisqr = chisqr 
     &     + chi2(THEORY(1,NPT,V,CURPAR),DATA(1,NPT),DATA(2,NPT))
   10 CONTINUE 
c
c  Add CDF chi square if correlations are used
c
      chisqr0 =chisqr
      if(ncov.ne.0)then
         chisqr=chisqr0+cdfchi2(data,NX,NNPTS,curpar,cdfxsec)
      endif
c
c  add normalization chi square penalty
c
      do 11 i=1,nfit
      if(errn(i).eq.0.d0.or.inorm(i).eq.100)goto 11
      chisqr=chisqr+chi2(curpar(inorm(i)),1d0,errn(i))
 11   continue
      RETURN 

  100 IREJ=0 
C CALCULATE ALPHA AND BETA
      DO 20 I=1,NVAR
      BETA(I)=0.0 
      DRV(I)=0.0 
      DO 20 J=1,NVAR
   20 ALPHA(I,J)=0.0
      DO 50 NPT=1,NNPTS 
C IGNORE DATA POINTS WITH NEGATIVE OR ZERO ERRORS
      IF(DATA(2,NPT).LE.0.0) GO TO 45 
      DO 41 II=1,NN 
   41 V(II)=DATA(II+2,NPT) 
      FVAL=THEORY(1,NPT,V,CURPAR) 
C IF MODE EQUALS 3 USER SUPPLIES DERIVATIVES IN ARRAY DRV
C OTHERWISE, PROGRAM CALCULATES DERIVATIVES USING FINITE DIFFERENCE
      IF(MODE.EQ.3) GO TO 125 
      DO 30 I=1,NVAR
      K=IFREP(I)
      SV=CURPAR(K) 
      CURPAR(K)=SV+PWATE(K)
C CALCULATE DERIVATIVE
      DRV(K)=(THEORY(MODE,NPT,V,CURPAR)-FVAL)/PWATE(K)
      if(itype(index).eq.201.or.itype(index).eq.202)then
      endif
      index=v(4)
      if(itype(index).eq.201.and.icorr(index).ne.0)then
         nc=npt-ncdfl+1
         cdfdrv(k,nc)=drv(k)
      endif
   30 CURPAR(K)=SV 
c
c  The parameter is reset to its original value, but if an evolution was 
c  performed, the valence and gluon normalizations are not reset. these are 
c  currently parameters 2, 8, and 25. We can put a call to pdfnrm here 
c  to reset them to the correct values. This only affects the output 
c  tables, since the values are overwritten each time an evolution is 
c  performed.
c
      call pdfnrm
c
c 
  125 continue
c
c  choose method to calculate alpha and beta
c
      if(icorr(index).eq.0)then
         DO 40 I=1,NVAR
         K=IFREP(I)
         BETA(I)=BETA(I)+(DATA(1,NPT)-FVAL)/DATA(2,NPT)*DRV(K)/
     2   DATA(2,NPT) 
         DO 40 J=1,I 
         L=IFREP(J) 
         ALPHA(I,J)=ALPHA(I,J)+DRV(K)*DRV(L)/DATA(2,NPT)**2
   40    ALPHA(J,I)=ALPHA(I,J)
      endif
      GO TO 50 
   45 IREJ=IREJ+1 
   50 CONTINUE 
c
c  calculate alpha and beta for the cdf points if required
c
      if(ncov.ne.0)then
      do ik=1,nvar
         k=ifrep(ik)
         do i=1,ncdf
            do ip=1,ncdf
               beta(ik)=beta(ik)+(data(1,ncdfl+i-1)-cdfxsec(i))
     2         *covinv(i,ip)*cdfdrv(k,ip)
            if(ik.eq.9)then
            endif
            enddo
         enddo
         do il=1,ik
            l=ifrep(il)
            do i=1,ncdf
               do ip=1,ncdf
                  alpha(ik,il)=alpha(ik,il)+
     2            cdfdrv(k,i)*covinv(i,ip)*cdfdrv(l,ip)
               enddo
            enddo
         enddo
         alpha(il,ik)=alpha(ik,il)
      enddo
      endif
c
c  calculate alpha and beta for the normalization parameter terms
c
      do 12 i=1,nvar
      do 13 k=1,nfit
      if(errn(k).eq.0.d0.or.inorm(k).eq.100)goto 13
      if(ifrep(i).eq.inorm(k))then
         beta(i)=beta(i)+(1.d0-curpar(ifrep(i)))/errn(k)**2
         alpha(i,i)=alpha(i,i)+1.d0/errn(k)**2
      endif
 13   continue
 12   continue
      RETURN 
      END 



      SUBROUTINE MINV(ARRAY,NORDER,DET) 
      implicit real*8 (a-h,o-z)
      DIMENSION ARRAY(50,50),IK(50),JK(50)
C
C  MATRIX INVERSION ROUTINE FROM BEVINGTON, PAGE 302.
C
      DET=1. 
      DO 100 K=1,NORDER 
      AMAX=0.D0 
   21 DO 30 I=K,NORDER 
      DO 30 J=K,NORDER 
c      IF( ABS(AMAX).GT. ABS(ARRAY(I,J))) GO TO 30 
      IF(DABS(AMAX).GT.DABS(ARRAY(I,J))) GO TO 30 
   24 AMAX=ARRAY(I,J)
      IK(K)=I 
      JK(K)=J 
   30 CONTINUE 
      IF(AMAX.NE.0.D0) GO TO 41 
   32 DET=0.
      GO TO 140 
   41 I=IK(K) 
      IF(I-K) 21,51,43 
   43 DO 50 J=1,NORDER 
      SAVE=ARRAY(K,J) 
      ARRAY(K,J)=ARRAY(I,J) 
   50 ARRAY(I,J)=-SAVE 
   51 J=JK(K) 
      IF(J-K) 21,61,53 
   53 DO 60 I=1,NORDER
      SAVE=ARRAY(I,K) 
      ARRAY(I,K)=ARRAY(I,J) 
   60 ARRAY(I,J)=-SAVE 
   61 DO 70 I=1,NORDER 
      IF(I.EQ.K) GO TO 70 
   63 ARRAY(I,K)=-ARRAY(I,K)/AMAX 
   70 CONTINUE 
      DO 80 I=1,NORDER 
      DO 80 J=1,NORDER 
      IF(I.EQ.K.OR.J.EQ.K) GO TO 80 
   75 ARRAY(I,J)=ARRAY(I,J)+ARRAY(I,K)*ARRAY(K,J) 
   80 CONTINUE 
      DO 90 J=1,NORDER 
      IF(J.EQ.K) GO TO 90 
   83 ARRAY(K,J)=ARRAY(K,J)/AMAX 
   90 CONTINUE 
      ARRAY(K,K)=1./AMAX 
  100 DET=DET*AMAX 
      DO 130 L=1,NORDER 
      K=NORDER-L+1
      J=IK(K) 
      IF(J.LE.K) GO TO 111 
  105 DO 110 I=1,NORDER 
      SAVE=ARRAY(I,K) 
      ARRAY(I,K)=-ARRAY(I,J) 
  110 ARRAY(I,J)=SAVE 
  111 I=JK(K) 
      IF(I.LE.K) GO TO 130 
  113 DO 120 J=1,NORDER 
      SAVE=ARRAY(K,J) 
      ARRAY(K,J)=-ARRAY(I,J) 
  120 ARRAY(I,J)=SAVE 
  130 CONTINUE 
  140 RETURN 
      END 
