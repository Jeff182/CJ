      SUBROUTINE DIST(X,pdf)
      implicit real*8 (a-h,o-z)
      DIMENSION pdf(-5:5)
      COMMON/QS/Q2 
      COMMON/GFUNC/calc(11,60,60)
c      common/input/xmc,xmb,alambda5,ncb,delta,s0,sb
      common/input/ILOOP,IORD,NMAX,IVL,XMC,XMB,NCB
      character*50 filename
      common/inputfile/filename
      common/threshold/sb
      common/stepsize/delta
      DATA ITIME/0/
      save s0,alambda5
      IF(ITIME.EQ.0) THEN
         OPEN(UNIT=4,FILE=filename,STATUS='UNKNOWN')
         read(4,11) xmc,xmb,alambda5,ncb
 11      format(3f10.4,i5)
         READ(4,10) calc
   10    FORMAT(8E15.9)
         CLOSE(4)
         s0=dlog(xmc**2/alambda5**2)
         sb=dlog(dlog(xmb**2/alambda5**2)/s0)
         delta=sb/ncb
         ITIME=1
      ENDIF
      sbar=dlog(dlog(q2/alambda5**2)/s0)
      CALL GINTERPsa(1,X,Sbar,glue) 
      CALL GINTERPsa(2,X,Sbar,sing) 
      CALL GINTERPsa(3,X,Sbar,up) 
      CALL GINTERPsa(4,X,Sbar,dp) 
      CALL GINTERPsa(5,X,Sbar,sp) 
      CALL GINTERPsa(6,X,Sbar,cp)
      call ginterpsa(7,x,sbar,singm)
      call ginterpsa(8,x,sbar,um)
      call ginterpsa(9,x,sbar,dm)
      nf=5
      if(sbar.lt.sb) nf=4
      uv=um+singm/nf
      dv=dm+singm/nf
      upub=up+sing/nf
      dpdb=dp+sing/nf
      ub=.5*(upub-uv)
      db=.5*(dpdb-dv)       
      s=.5*(sp+sing/nf)
      c=.5*(cp+sing/nf)
      if(nf.eq.5)then
         bp=-(up+dp+sp+cp)
         b=.5*(bp+sing/nf)
      else
         b=0.
      endif
c
c  added 1/8/10 as a test
c
      if(x.lt.0.80d0)then
         pdf(1)=uv+ub
         pdf(2)=dv+db
      else
         pdf(1)=uv
         pdf(2)=dv
      endif
      pdf(0)=glue
      pdf(-1)=ub
      pdf(-2)=db
      pdf(-3)=s
      pdf(-4)=c
      pdf(-5)=b
      do 100 j=-5,-1
      if(pdf(j).lt.0.d0) pdf(j)=0.d0
 100  continue
      pdf(3)=pdf(-3)
      pdf(4)=pdf(-4)
      pdf(5)=pdf(-5)
      RETURN
      END 


      SUBROUTINE GINTERPsa(I,X,S,ANS)
C     THIS ROUTINE INTERPOLATES AS NEEDED IN THE ARRAY GF TO OBTAIN THE
C     EVOLVED DISTRIBUTIONS AT THE REQUIRED Q2 AND X VALUES
C     This is a stand-alone version of the ginterp routine used 
C     in the fitting package
      implicit real*8 (a-h,o-z)
      DIMENSION F1(60),F2(60) 
      COMMON/GFUNC/GF(11,60,60) 
c      common/input/xmc,xmb,alambda5,ncb,delta,s0,sb
      common/input/iloop,iord,nmax,ivl,xmc,xmb,ncb
      common/threshold/sb
      common/stepsize/delta
      if(s.lt.sb)then         
         is=s/delta+1
         s1=(is-1)*delta
         s2=s1+delta
      else
         is=(s-sb)/delta+14
         s1=(is-14)*delta+sb
         s2=s1+delta
      endif
      IS1=IS+1
      DO 1 L=1,60
      F1(L)=GF(I,IS,L) 
      F2(L)=GF(I,IS1,L)
    1 CONTINUE
      A1=GETFV(X,F1)
      A2=GETFV(X,F2)
      ANS=A1*(S-S2)/(S1-S2)+A2*(S-S1)/(S2-S1) 
      RETURN
      END 
      FUNCTION GETFV(X,FVL)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION FVL(60),xx(4),fx(4),xgrid(60)
      DATA NX,XGRID/60,.00001,.000015,.00002,.00004,.00006,.00008,
     2.0001,.00015,.0002,.0004,.0006,.0008,.001,.0015,.002,.004,.006,
     3.008,.01,.015,.02,.04,.06,.08,.1,.125,.15,.175,.2,.225,.25,.275,
     4.3,.325,.35,.375,.4,.425,.45,.475,.50,.525,.55,.575,.60,.625,
     5.65,.675,.70,.725,.75,.775,.80,.825,.85,.875,.90,.925,.95,.975/
      DO 1 I=1,NX 
      IF(X.LT.XGRID(I)) GO TO 2
    1 CONTINUE
    2 I=I-1
      IF(I.EQ.0) I=1
      if(i.gt.nx-2)i=nx-2
c      IF(I.GT.NX-3) I=NX-3
      R12=FVL(I+1)/FVL(I)
      R13=FVL(I+2)/FVL(I)
      IF(R12.GT.0.D0.AND.R13.GT.0.D0)THEN
         J=I+1
         K=J+1
         AXI=DLOG(XGRID(I))
         BXI=DLOG(1.-XGRID(I))
         AXJ=DLOG(XGRID(J))
         BXJ=DLOG(1.-XGRID(J))
         AXK=DLOG(XGRID(K))
         BXK=DLOG(1.-XGRID(K))
         FI=DLOG(DABS(FVL(I)))
         FJ=DLOG(DABS(FVL(J)))
         FK=DLOG(DABS(FVL(K)))
         DET=AXI*(BXJ-BXK)+AXJ*(BXK-BXI)+AXK*(BXI-BXJ)
         ALOGA=(FI*(AXJ*BXK-AXK*BXJ)+FJ*(AXK*BXI-AXI*BXK)+
     $   FK*(AXI*BXJ-AXJ*BXI))/DET
         ALPHA=(FI*(BXJ-BXK)+FJ*(BXK-BXI)+FK*(BXI-BXJ))/DET
         BETA=(FI*(AXK-AXJ)+FJ*(AXI-AXK)+FK*(AXJ-AXI))/DET
         if(dabs(alpha).gt.50.d0.or.dabs(beta).gt.50.d0)goto 3
         isign=fvl(i)/dabs(fvl(i))
         ANS=isign*DEXP(ALOGA)*X**ALPHA*(1.-X)**BETA
      ELSE
 3       xx(1)=xgrid(i)
         xx(2)=xgrid(i+1)
         xx(3)=xgrid(i+2)
         fx(1)=fvl(i)
         fx(2)=fvl(i+1)
         fx(3)=fvl(i+2)
         if(i.eq.nx-2)then
            xx(4)=1.
            fx(4)=0.
         else
            xx(4)=xgrid(i+3)
            fx(4)=fvl(i+3)
         endif
         CALL POLINT(XX,FX,4,X,ANS,DY)
      ENDIF
      getfv=ans
      RETURN
      END 

