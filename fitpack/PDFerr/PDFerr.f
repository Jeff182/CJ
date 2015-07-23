****************************************************************
*     Calculation of the CJ PDF error sets
*
*     produce 2*NFREE .par files and .gf files corresponidng
*     to the "eigenPDFs" needed for PDF error calculations
*     (Diagonalizes the error matrix, and shifts the parameters according
*     to the eigenvectors normalized to the eigenvalue).
*
*     Programmers: aa = A.Accardi (accardi@jlab.org)
*
*     HISTORY:
*
*     PDFerr      _the first one_
*     (28 Jul 11)
*
****************************************************************
      program PDFerrors

      implicit none

      character*50 file,parfile,errfile,gffile,tblfile,pdffile
      character*2 cpdf 
 
      integer len,len1,inuke,itmc,iht,ipar,i,j,ind_i,ind_j,npdf,k
     &     ,ipdf, itest
      integer*4 today(3),now(3)

      integer np
      parameter (np=50)
      double precision Amat(np,np),Bvec(np),epsilon(np,np),eval(np)
     &     ,E(np),norm,norm2,par0(100),indx(np),d,delta,tmp,fac 

      double precision del(np,np),det,myeps(np)
       
*     alpha matrix and parameter errors
      integer pmax
      parameter(pmax=50) ! must be = dimension of alpha in 'minim'
      double precision alfa(pmax,pmax),alfainv(pmax,pmax)
     &     ,central(pmax),eps(pmax),cov(pmax,pmax)
      integer pos(pmax),nfree
      common/PDFerr/alfa,cov,central,eps,pos,nfree

*     parameters for QCD evolution
      double precision par(100)
      common/param/par

*     Other PDF parameters 
      character*10 pname(100)  ! needed if writing a par file with 'writepar'
      double precision uncrt(100)
      common/oparam/pname,uncrt

*     Total number of parameters
      integer npar              
      common/npar/npar

*     comment lines in .par file
      character*100 line(3)
      common/parfil/line


      if (iargc().eq.0) then    ! HELP
         print*, 'USAGE: PDFerr <file> [<ipdf>]'
         print*, 'where file = .par file from fitpack (no suffix)'
         print*, '      ipdf = 1 --> writes .pdf files'
         stop
      end if

      call getarg(1,file)       ! gets file name (no suffix, please)
      if(iargc().eq.1) then
         ipdf = 1               ! default: 1 = writes sample .pdf files
      else
         call getarg(2,cpdf)
         read(cpdf,'(i10)') ipdf
      end if


      call trmstr(file,len)
      
      parfile = file(1:len)//'.par'

      ! reads the central value pars, evolves and writes .par and .gf file
      write(*,*) 'Working on '//file
      call readpar(parfile,par,inuke,itmc,iht)

      print*, '*1'

c jfo alfa is the Hessian matrix from the fit

      call tredn2(alfa,nfree,pmax,eval,E) ! NOTE: alfa is destroyed 
      print*, '*2'
      call tqli(eval,E,nfree,pmax,alfa)
c jfo
c jfo copy the eigenvectors into epsilon to be compatible with the 
c jfo previous code
c jfo
      do k=1,nfree
      do i=1,nfree
         epsilon(k,i)=alfa(k,i)
      enddo
      enddo  
      print*, '*3'
      ! on exit, eval contains eigenvalues, and the Kth column of epsilon 
      ! the Kth eigenvector, i.e., epsilon(n,K) = K-eigenvec(n)

c jfo Eigenvectors are normalized to unity
c jfo Change normalization to be equivalent to choosing the scale factors 
c jfo s_k=1/dsqrt(eval(k)) as introduced in the Tung Hessian paper. 
c jfo 
      do k=1,nfree                ! cycles thorugh the eigenvectors
         norm = 1./dsqrt(eval(k))
         do i=1,nfree
            epsilon(i,k) = norm*epsilon(i,k)
         end do
      end do                    ! Now epsilon(*,k) contains the normalized 
                                ! k-th eigenvector
c jfo      print*,(eval(i),i=1,nfree)
c jfo      do k=1,nfree
c jfo         print*,k,(epsilon(i,k),i=1,nfree)
c jfo      enddo

*    *** PDFs central value

      call QCDev      
      ! Output .par, .gf, .tbl, .pdf files for central value
      parfile = file(1:len)//'_00.par'
      call trmstr(line(1),len1)
      tmp = alfa(1,1)

      alfa(1,1) = 0d0           ! to suppress output of alfa
      call writepar(parfile,line,pname,par,uncrt,alfa,nfree
     &     ,inuke,itmc,iht)
      alfa(1,1) = tmp

      gffile  = file(1:len)//'_00.gf'
      call writegf(gffile)

      tblfile  = file(1:len)//'_00.tbl'
      call writetbl(tblfile)

      if (ipdf.eq.1) then
         pdffile  = file(1:len)//'_00.pdf'
         call writepdf(pdffile)
      end if


*    *** Calculates the PDF error sets
      itest=1
      if(itest.eq.2)stop
      npdf = 0
      do i = 1, npar
         par0(i) = par(i)
      end do
      do ipar = 1, nfree

         print*
         print*, 'Working on parameter #',ipar,' of',nfree

         ! creates and writes new par file
         do i = 1, npar
            par(i) = par0(i)
         end do
         do i=1,-1,-2           ! positive and negative shifts
            npdf = npdf + 1
            ind_j = 0
            do j=1,nfree        ! shifted parameter set
c jfo               par(pos(j)) = par0(pos(j)) + i*epsilon(i,j)
c jfo
c jfo For a chi square tolerance T, one should shift the parameters 
c jfo by an amount sqrt(T)/2 times the error
c jfo T=50 is used here
c jfo
c               fac=3.5355
               fac=1.d0
c
c  fac=1 corresponds to delta chi = 1 provided that one divides the
c  error on an observable by 2
c  delta sigma = 1/2 sqrt(sum(sigma(a+)-sigma(a-)**2)
c  where the sum runs over the varied parameters denoted by 'a'
c  and a+ (a-) are the parameters in the odd(even) numbered PDF error sets
c
 
               par(pos(j)) = par0(pos(j)) + i*fac*epsilon(j,ipar)
            end do

            call QCDev          ! performs QCD evolution 
                                ! (and normalizes valence and glue)
            
            write(cpdf,'(I2.2)'),npdf
            errfile = file(1:len)//'_'//cpdf//'.par'

            line(1) = line(1)(1:len1)//' - pdf set #'//cpdf
            call idate(today)   ! today(1)=day, (2)=month, (3)=year
            call itime(now)     ! now(1)=hour, (2)=minute, (3)=second
            write(line(3),1000) today(2), today(1), today(3), now
 1000       format ( '# Created on ', i2.2, '/', i2.2, '/', i4.4, ' at '
     &           ,i2.2, ':', i2.2, ':', i2.2 )
            tmp = alfa(1,1)
            alfa(1,1) = 0d0     ! to suppress output of alfa
            call writepar(errfile,line,pname,par,uncrt,alfa,nfree
     &           ,inuke,itmc,iht)
            alfa(1,1) = tmp

            ! CJ grids
            gffile  = file(1:len)//'_'//cpdf//'.gf'
            call writegf(gffile)

            ! CTEQ6-series grids
            tblfile  = file(1:len)//'_'//cpdf//'.tbl'
            call writetbl(tblfile)

            ! sample .pdf files for quality control
            if (ipdf.eq.1) then
               pdffile = file(1:len)//'_'//cpdf//'.pdf'
               call writepdf(pdffile)
            end if

         end do


      end do

      stop
      end
