      PROGRAM FILLGF
      implicit real*8 (a-h,o-z)
      character*50 parafile, output
      character*6 tmp
      character*10 ipname(100)
      DIMENSION PAR(100), pwate(100) 
      dimension sp(60),sm(60),cp(60),cm(60),bp(60),bm(60)
c      COMMON/INPUT/ILOOP,IORD,NMAX,IVL,xmc,xmb,ncb,ndeg 
      COMMON/INPUT/ILOOP,IORD,NMAX,IVL,xmc,xmb,ncb
      COMMON/GRPTHY/FLAVOR
      COMMON/GFUNC/CALC(11,60,60)
      common/grid/nx,xgrid(60) 
      common/param/par
      common/threshold/sb
      common/stepsize/delta
      data xmc,xmb,ncb/1.3d0,4.5d0,12/


*    ***  OPEN I/O FILES
      write(6,11)
 11   format(/' Enter filename for PDF parameter file')
      read(5,*) parafile
      print*, '* input = ',parafile
      tmp='calc_'
      call trmstr(tmp,len2)
      call trmstr(parafile,len1)
      if(parafile(len1-3:len1).eq.'.par')len1=len1-4
      output=tmp(1:len2)//parafile(1:len1)
      open(unit=2,file=output,status='unknown')
      open(unit=3,file=parafile,status='old')

*    *** initializations
*    ... initializes PAR array
      DO 20 J=1,100
 20      PAR(J)=0.
*    ... Gaussian integrations needed by 'intqcd' subroutine    
      CALL WATE16
      CALL WATE32
*    ... CALC matrix
      print*,'nx =',nx
      DO 10 I=1,11 
      DO 10 J=1,60
      DO 10 K=1,nx
 10      CALC(I,J,K)=0.
*    *** reads PDF params from file
      do 31 j=1,100
         read(3,*,end=32)ipname(j),par(j),pwate(j)
 31   continue
 32   jmax=j-1
      do 40 j=1,jmax
         print*,ipname(j),par(j),pwate(j)
 40   continue

*    *** DGLAP evolution
*    ... sets parameters for evolution subroutine
c      NDEG=2
      FLAVOR=4. 
      ILOOP=2
      IORD=1                    ! LO: IORD=0    NLO: IORD=1
      S0=DLOG(xmc**2/PAR(1)**2) 
      SMAX=DLOG(DLOG(462400./PAR(1)**2)/S0) 
      SB=DLOG(DLOG(XMB**2/PAR(1)**2)/S0)
      DELTA=SB/NCB
      NMAX=SMAX/DELTA+3 
      IF(NMAX.GT.60) NMAX=60
      print*,s0,sb,delta,ncb,smax,nmax
*    ... calls PDF evolution routine
      CALL INTQCD(PAR)

*    *** writes output to file
      write(2,101) xmc,xmb,par(1),ncb
  101 format(3f10.4,i5)
      WRITE(2,100) CALC
  100 FORMAT(8E15.9)

      CALL EXIT 
      END 
