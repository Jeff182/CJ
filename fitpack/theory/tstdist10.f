      PROGRAM TSTDIST
      implicit real*8 (a-h,o-z)
      character*50 filename
      DIMENSION pdf(-5:5),xgrid(30)
      COMMON/QS/Q2
      common/inputfile/filename
      DATA NX,XGRID/30,.0001,.0002,.0004,.0006,.0008,.001,.002,.004,
     2.008,.016,.032,.064,.1,.15,.2,.25,.3,.35,.4,.45,.5,.55,.6,.65,
     3.7,.75,.80,.85,.9,.95/
*
C
C  SET UP FLAGS, I/O FILES, ETC.
C
      filename='calc_tstfit10_nofit2'
      Q02=1.69
      DO 100 J=1,3
      Q2=Q02
      IF(J.EQ.2) Q2=100.
      if(j.eq.3) q2=1000.
      PRINT 90,Q2
   90 FORMAT(/,' Q2=',F7.2,' GEV**2')
      PRINT 91
   91 FORMAT(/,7X,' x',9X,' xd',9X,' xu',10X,' xg',9X,'xub',
     2     9X,'xdb',9x,' xs',9x,' xc',9x,' xb')
      DO 80 K=1,nx
      x=xgrid(k)
      CALL DIST(x,pdf)
      print 92, x,(pdf(l),l=2,-5,-1)
   92 FORMAT(9E12.4)
   80 CONTINUE
  100 CONTINUE
      CALL EXIT
      END
