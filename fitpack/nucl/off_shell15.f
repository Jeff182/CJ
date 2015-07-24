C  This File contains code necessary for the flavor modified off shell
C  corrections.  These off shell corrections follow the same general
C  model as the modified Kulagin Petti corrections but modified even
C  further to allow different flavors of quarks.
C
C  Created 7/10/13 by L.T.Brady [lucas_brady@hmc.edu]
C  Currently supports valence quarks, sea quarks (u and d), and gluons
C  Implemented in such a way that it will propogate correctly through
C  any possible observable in DIS13.f
C
C  Usage:
C        Access this option in the input file by using an inuke of the
C        form:   inuke = FED5BA   (i.e. include a 5 in the C space)
C        
C        Before accessing anything else, call set_offshell with the
C        correct parameters and .true. in the last spot
C
C IMPORTANT: whenever you are done with calculating offshell corrections
C            call set_offshell again with .false. in the last flag
C            this turns off the off shell corrections and reverts
C            DIS10 to its original state.  Things will be bad otherwise



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine get_qder (x,iflavor,qder)
C
C     A specific combination of things including the quark derivative
C     this combination is needed for the offshell correction model
C     specifically it calculates
C     (1/q)*\pder{q}{x}*h(x)
C     iflavor is an integer flag indicating which flavor to use:
C           1     =     valence
C           2     =     sea
C           3     =     gluons
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Implicit None
      double precision x,qder, delt
      double precision h
      double precision fdv, fuv, fubpdb, fcnfg
      integer iflavor
      double precision a0,a1,a2,a3,a4
      double precision b,c,a0d,a1d,a2d,a3d,a4d
      double precision xp(100)
      common/param/xp
      call get_h(x,iflavor,h)
      delt = 1.0d-4
c
c  Don't allow x+delt to be greater than 1.d0
c
      if(1.0-x.le.delt)delt=-delt
      if (iflavor.eq.1) then        ! valence
c            a0=xp(2)
c            a1=xp(3)    ! up part
c            a2=xp(4)
c            a3=xp(5)
c            a4=xp(6)
c            b = xp(32)  ! down part
c            c = xp(33)
c            a0d=xp(8)
c            a1d=xp(9)
c            a2d=xp(10)
c            a3d=xp(11)
c            a4d=xp(12)
c            qder=(h*(a0d*(1 - x)**a2d*x**a1d*
c     -            (2 + a3d*Sqrt(x)*(1 - x + 2*a2d*x) + 
c     -            2*a1d*(-1 + x)*
c     -            (1 + a3d*Sqrt(x) + a4d*x) + 
c     -            2*x*(-1 + a2d + a2d*a4d*x)) + 
c     -            a0*(1 - x)**a2*x**a1*
c     -            (2 + 2*x*(-1 + a2 + a2*a4*x) + 
c     -            a3*(Sqrt(x) + (-1 + 2*a2)*x**1.5 - 
c     -            a0d*b*(-1 + 2*c)*x**(0.5 + c) + 
c     -            a0d*b*(-1 + 2*a2 + 2*c)*x**(1.5 + c))
c     -            + 2*a0d*b*x**c*
c     -            (1 + c*(-1 + x)*(1 + a4*x) + 
c     -            x*(-1 + a2 + a2*a4*x)) + 
c     -            2*a1*(-1 + x)*
c     -            ((1 + a4*x)*(1 + a0d*b*x**c) + 
c     -            a3*(Sqrt(x) + a0d*b*x**(0.5 + c))))))/
c     -            (2.*(-1 + x)*x*
c     -            (a0d*(1 - x)**a2d*x**a1d*
c     -            (1 + a3d*Sqrt(x) + a4d*x) + 
c     -            a0*(1 - x)**a2*x**a1*(1 + a3*Sqrt(x) + a4*x)*
c     -            (1 + a0d*b*x**c)))
c         delt=.001*(1.d0-x)
         qder=(fdv(x+delt)+fuv(x+delt))/(fdv(x)+fuv(x))-1.
         qder=qder/delt
c
c  fsupdf returns xq, xg, etc the -1./x corrects the derivative to 
c  1/q dq/dx etc
c
         qder=h*(qder-1./x)
      else if (iflavor.eq.2) then   ! sea
c            a1=xp(15)
c            a2=xp(16)
c            a3=xp(17)
c            a4=xp(18)
c            qder = (h*(2 + a3*Sqrt(x)*(1 - x + 2*a2*x) + 
c     -            2*a1*(-1 + x)*(1 + a3*Sqrt(x) + a4*x) + 
c     -            2*x*(-1 + a2 + a2*a4*x)))/
c     -            (2.*(-1 + x)*x*(1 + a3*Sqrt(x) + a4*x))
c         delt=.001*(1.d0-x)
         qder=fubpdb(x+delt)/fubpdb(x)-1.
         qder=qder/delt
         qder=h*(qder-1./x)
      else if (iflavor.eq.3) then   ! gluons
c            a1=xp(26)
c            a2=xp(27)
c            a3=xp(28)
c            a4=xp(29)
c            qder = (h*(2 + a3*Sqrt(x)*(1 - x + 2*a2*x) + 
c     -            2*a1*(-1 + x)*(1 + a3*Sqrt(x) + a4*x) + 
c     -            2*x*(-1 + a2 + a2*a4*x)))/
c     -            (2.*(-1 + x)*x*(1 + a3*Sqrt(x) + a4*x))
c         delt=.001*(1.d0-x)
         qder=fcnfg(x+delt)/fcnfg(x)-1.
         qder=qder/delt
         qder=h*(qder-1./x)
      else
            write(*,*) "Error: Offshell model: iflavor out of range"
      endif
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine get_h (x,iflavor,h)
C
C     Computes h(x) which is a set of constants needed for the fmKP offshell
C     corrections.  This function relies on sbar, the average spectator mass squared
C     outputs the result as h
C     iflavor is an integer flag:
C           1     =     valence
C           2     =     sea
C           3     =     gluons
C     Note, h(x) here is equivalent to x(1-x)h(x) from Kulagin-Petti
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      double precision x, h
      integer iflavor
      double precision s, sv, ssq, sg, lambda
      double precision Mp
      
      Mp = 0.938272D0
      
      call get_offshell_s(sv,ssq,sg,lambda)
      if (iflavor.eq.1) then
            s=sv
      else if (iflavor.eq.2) then
            s=ssq
      else if (iflavor.eq.3) then
            s=sg
      else
            write(*,*) "Error: Offshell model: iflavor out of range"
      endif
      
      h=x*(1-x)*((lambda*s)/Mp**2 + (1 - lambda)*(1 - x))/
     &  (-(s/Mp**2) + (1 - x)**2)
     
      return
      end
      

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine get_delta (ndrv,x,dfv,dfsq,dfg,dfsb,dfcb,dfbb)
C
C     determines the offshell correction to the valence quark
C     input parameters are a specific x, s, c, and lambda
C     outputs in dfv.  This was made in Mathematica
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Implicit None
      double precision x,dfv,dfsq,dfg,dfsb,dfcb,dfbb
      double precision sv, ssq, sg, lambda, old_lambda
      double precision c_v, c_sq, c_g
      double precision logis
      double precision cv(0:35), csq(0:35), cg(0:35), xlambda(0:35)
C     Determine whether we need this offshell correction      
      logical onoff
      integer ndrv
      data old_lambda/-1.0d9/
      common/C_array/cv, csq, cg, xlambda
      save old_lambda
      save c_v,c_sq,c_g

c
c modified get_delta routines --jfo 3/24/15
c allows for variable off_shell strength to be fitted
c combines all 'get_delta' type calls into one routine
c
c      call get_offshell_on(onoff)
      
c      call get_offshell_c(c_v,c_sq,c_g)
           
c      if (onoff) then 
c         call set_offshell_par
c         call get_offshell_s(sv,ssq,sg,lambda)
c         if(lambda.ne.old_lambda)then
c            old_lambda=lambda
c            call init_offshell
c         endif

          call get_qder(x,1,dfv)
          dfv= dfv + cv(ndrv)
          call get_logistic(x, logis)
          dfv = dfv*logis
          call get_qder(x,2,dfsq)
          dfsq= dfsq + csq(ndrv)
          call get_logistic(x, logis)
          dfsq = dfsq*logis
          call get_qder(x,3,dfg)
          dfg= dfg + cg(ndrv)
          call get_logistic(x, logis)
          dfg = dfg*logis
          dfsb=dfsq
          dfcb=0.
          dfbb=0.
c      else
c          dfv = 1.D0
c          dfsq=1.d0
c          dfg=1.d0
c          dfsb=1.d0
c          dfcb=1.d0
c          dfbb=1.d0
c      endif
      return
      end


cCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c      subroutine get_delta_v (x,dfv)
cC
cC     determines the offshell correction to the valence quark
cC     input parameters are a specific x, s, c, and lambda
cC     outputs in dfv.  This was made in Mathematica
cCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c      Implicit None
c      double precision x,dfv
c      double precision sv, ssq, sg, lambda
c      double precision c_v, c_sq, c_g
c      double precision logis
cC     Determine whether we need this offshell correction      
c      logical onoff
c      call get_offshell_on(onoff)
c      
c      call get_offshell_c(c_v,c_sq,c_g)
c      call get_offshell_s(sv,ssq,sg,lambda)
c      
c      if (onoff) then 
c          call get_qder(x,1,dfv)
c          dfv= dfv + c_v
c          call get_logistic(x, logis)
c          dfv = dfv*logis
c      else
c          dfv = 1.D0
c      endif
c      
c      return
c      end



cCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c      subroutine get_delta_sq (x,dfsq)
cC
cC     determines the offshell correction to the sea quarks
cC     input parameters are a specific x, s, c, and lambda
cC     outputs in dfsq.  This was made in Mathematica
cCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c      Implicit None
c      double precision x,dfsq
c      double precision sv, ssq, sg, lambda
c      double precision c_v, c_sq, c_g
c      double precision logis
cC     Determine whether we need this offshell correction      
c      logical onoff
c      call get_offshell_on(onoff)
c      
c      call get_offshell_c(c_v,c_sq,c_g)
c      call get_offshell_s(sv,ssq,sg,lambda)
c      
c      if (onoff) then 
c          call get_qder(x,2,dfsq)
c          dfsq= dfsq + c_sq
c          call get_logistic(x, logis)
c          dfsq = dfsq*logis
c      else
c          dfsq = 1.D0
c      endif
c
c     return
c     end



cCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c     subroutine get_delta_g (x,dfg)
cC
cC     determines the offshell correction to the gluons
cC     input parameters are a specific x, s, c, and lambda
cC     outputs in dfg.  This was made in Mathematica
cCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c      Implicit None
c      double precision x,dfg
c      double precision sv, ssq, sg, lambda
c      double precision c_v, c_sq, c_g
c      double precision logis
C     Determine whether we need this offshell correction      
c      logical onoff
c      call get_offshell_on(onoff)
c      
c      call get_offshell_c(c_v,c_sq,c_g)
c      call get_offshell_s(sv,ssq,sg,lambda)
c      
c      if (onoff) then 
c          call get_qder(x,3,dfg)
c          dfg= dfg + c_g
c          call get_logistic(x, logis)
c          dfg = dfg*logis
c      else
c          dfg = 1.D0
c      endif
c
c      return
c      end
      
      
cCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c      subroutine get_delta_cb (x,dfcb)
cC
cC     determines the offshell correction to the charm quarks
cC     input parameters are a specific x, s, c, and lambda
cC     outputs in dfg.
cC     Not yet implemented fully - assumes offshell charm contribution is zero
cCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c      Implicit None
c      double precision x,dfcb
c      logical onoff
c      call get_offshell_on(onoff)
c      
c      if (onoff) then 
c          dfcb = 0.d0
c      else
c          dfcb = 1.D0
c      endif
c
c      return
c      end
      
      
cCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c      subroutine get_delta_sb (x,dfsb)
cC
cC     determines the offshell correction to the strange quarks
cC     input parameters are a specific x, s, c, and lambda
cC     outputs in dfsb.
cCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c      Implicit None
c      double precision x,dfsb
c      logical onoff
c      call get_offshell_on(onoff)
c      
c      if (onoff) then 
c          call get_delta_sq(x,dfsb)
c      else
c          dfsb = 1.D0
c      endif
c
c      return
c      end
      
      
cCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c      subroutine get_delta_bb (x,dfbb)
cC
cC     determines the offshell correction to the bottom quarks
cC     input parameters are a specific x, s, c, and lambda
cC     outputs in dfg.
cC     Not yet implemented fully - assumes offshell bottom contribution is zero
cCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c      Implicit None
c      double precision x,dfbb
c      logical onoff
c      call get_offshell_on(onoff)
c      
c      if (onoff) then 
c          dfbb = 0.d0
c      else
c          dfbb = 1.D0
c      endif
c
c      return
c      end
      
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine get_offshell_on (onoff)
C     returns the current value of offshell_on
C     onoff - Logical - whether or not the offshell corrections
C                          are on or off.  If true, then DIS10 will 
C                          compute the off shell corrections instead of
C                          the normal DIS results
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      Logical onoff
      logical offshell_on
C     If the offshell common block has not bee initialized offshell_on=false
      data offshell_on /.false./
      common /offshell/ offshell_on
      
      onoff = offshell_on
      
      return
      end
      
      
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine get_offshell_c (cv, csq, cg)
C     returns the current value of offshell normalization constants
C     All the outputs are double precision
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      double precision cv, csq, cg
      double precision c_v, c_sq, c_g
      common /C_block/ c_v, c_sq, c_g
      
      cv=c_v
      csq=c_sq
      cg=c_g
      
      return
      end
      
    
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine get_offshell_s (s_v, s_sq, s_g,l)
C     returns the current value of offshell s and lambda constants
C     All the outputs are double precision
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      double precision s_v,s_sq,s_g,l
      double precision sv, ssq, sg, lambda
      common /S_block/ sv, ssq, sg, lambda
      
      s_v=sv
      s_sq=ssq
      s_g=sg
      l=lambda
      
      return
      end
    

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine mk_cv (cv)
C     
C     makes the normalization constant for the valence offshell correction
C     This requires an integral which is currently computed using 
C     Simpson's rule.  Returns the normalization constant as cv
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c      implicit none
c      double precision x, cv, fdv, fuv, qder, dx, coeff, N
c      integer ix,nx
      implicit real*8 (a-h,o-z)
      common/gaus32/xi(32),wi(32),nterms,xx(33)
c      nx=25 ! No. of points for Simpson's rule integration
c      dx=1.D0/nx
      
      cv=0.D0
      den=0.D0
      
c      do ix=1,nx-1 ! Simpson's rule integration
c            if (ix.Eq.1.or.ix.Eq.nx-1) then 
c                  coeff=1.
c            else if (mod(ix,2).Eq.0) then
c                  coeff=4.
c            else
c                  coeff=2.
c            endif
            
c            x=ix*dx
c            call get_qder(x,1,qder)
c            qder = qder*(1./x)*(fdv(x)+fuv(x))
c            cv = cv+coeff*qder
c            N = N + coeff*(1./x)*(fdv(x)+fuv(x))
c      enddo
      do l=1,nterms
         z=0.5*(1.+xi(l))
         x=z**3
         call get_qder(x,1,qder)
         tmp=(fuv(x)+fdv(x))/x*3.*z**2         
         qder=qder*tmp
         cv=cv+0.5*wi(l)*qder
         den=den+0.5*wi(l)*tmp
       enddo      
c      cv=-cv/N
      cv=-cv/den      
      return
      end   
      
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine mk_cg_csq (cv, cg, csq)
C     
C     makes the normalization constant for the gluon and sea offshell
C     corrections.  This requires an integral which is currently 
C     computed using Simpson's rule.  Returns the normalization 
C     constants as cg and csq.  Takes in cv as an input
C     At the moment, this assumes that cg=csq
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c      implicit none
c      double precision x, cv, cg, csq, qder, dx, coeff
c      double precision fdv, fuv, fubpdb, fcnfg, FCNVL ! down valence, 
c      up valence, sea, gluon
c      double precision Delta ! defined as csq/cg.  Assumed to be 1 for now
c      double precision num, denom
c      integer ix,nx

c      nx=25 ! No. of points for Simpson's rule integration
c      dx=1.D0/nx
c      Delta = 1.D0

c      num=0D0
c      denom=0D0
      
c      do ix=1,nx-1 ! Simpson's rule integration
c            if (ix.Eq.1.or.ix.Eq.nx-1) then 
c                  coeff=1.
c            else if (mod(ix,2).Eq.0) then
c                  coeff=4.
c            else
c                  coeff=2.
c            endif
c            x=ix*dx
            
c            call get_qder(x,1,qder) ! Valence part
c            num = num + (cv+qder)*(fdv(x)+fuv(x))
c            call get_qder(x,2,qder) ! Sea part
c            num = num + 2*qder*fubpdb(x)
c            call get_qder(x,3,qder) ! Gluon part
c            num = num + qder*fcnfg(x)
c            call get_qder(x,2,qder) ! Strange part
c            num = num + qder*fcnvl(5,x)
            
c            denom = denom + 2*Delta*(fubpdb(x)+fcnvl(5,x)) + fcnfg(x)
c      enddo
      
c      cg  = -num/denom
c      csq =  Delta*cg
      implicit real*8 (a-h,o-z)
      common/gaus32/xi(32),wi(32),nterms,xx(33)
      cg=0.D0
      csq=0.D0
      anum=0.D0
      denom=0.D0
      Delta=1.D0
      do l=1,nterms
         z=0.5*(1.+xi(l))
         x=z**3
         fac=0.5*wi(l)*3.*z**2
         call get_qder(x,1,qder)
         c1=fdv(x)+fuv(x)
         anum=anum+fac*(cv+qder)*c1
         call get_qder(x,2,qder)
         c2=fubpdb(x)
         anum=anum+fac*2.*qder*c2
         call get_qder(x,3,qder)
         c3=fcnfg(x)
         anum=anum+fac*qder*c3
         call get_qder(x,2,qder)
         c4=(fcnvl(5,x)+.25*fcnfs(x))
         anum=anum+fac*qder*c4
         denom=denom+fac*Delta*(2.*c2+c3+c4)
      enddo
      cg=-anum/denom
      csq=Delta*cg  
      return
      end   
    
    
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine get_logistic (x, ans)
C     
C     returns a logistic curve.  This curve is designed to keep large x
C     behavior unchanged but suppress small x behavior to zero.  At the
C     moment, the function turns at x=0.15.  It is 92% by x=0.2 and is
C     7% by x=0.1.  If we want to change this function in the future,
C     the parameters are tunable.
C     A = Lower Asymptote, currently 0
C     m = mid-point of the downturn, currently x=0.15
C     B = strength of the downturn, currently 50
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC    
      implicit none
      double precision x, ans
      double precision A, m, B
      
      A = 0D0
      m = 0.15D0
      B = 50D0
      
      ans = A + (1D0-A)/(1D0+exp(-B*(x-m)))

      return
      end



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC    
      subroutine findoff(jmax,par,ipname,pwate)
*     finds which parameters are for off-shell corrections (if any)
*     and stores their index in noff. If no parameter nht(i)=0.

      implicit none
      integer jmax
      character*10 ipname(100)
      double precision par(100),pwate(100)

      integer i, j
      character*2 num

*     Off-shell model parameter indexes
      integer maxpar
      parameter(maxpar=14)
      integer noff(maxpar) ! max of 14 off-shell parameters 
      common/offpar/noff

      do i=1,maxpar
         noff(i)=0
         write(num,'(I0)') i
         do j=1,jmax
            if(ipname(j).eq.'off'//trim(num)) then
               noff(i)=j
             end if
         end do
      end do

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC    
      subroutine get_off(offpar)
*     returns the current value of the offshell parameters.
*     If a parameter is not present in the original input file
*     it returns -1D99 for that parameter. The calling routine shoudl take
*     care to test for missing parameters by checking against this value 

      implicit none

      double precision offpar(14) ! Note: size should be = maxpar below

      integer i

*     Off-shell model parameter indexes
      integer maxpar
      parameter(maxpar=14)
      integer noff(maxpar) ! max of 14 off-shell parameters 
      common/offpar/noff
*     Current fit parameters
      double precision xc(100)
      common/param/xc

      do i=1,maxpar
         if (noff(i).gt.0) then
            offpar(i) = xc(noff(i))
         else
            offpar(i) = -1D99
         end if
      end do

      return
      end


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine set_offshell_on (onoff)
C
C     This function sets on and off the off-shell corrections to the PDFs
C     as needed by the user. 
c
C        onoff - Logical - sets whether or not the offshell corrections
C                          are on or off.  If true, then DIS10 will 
C                          compute the off shell corrections instead of
C                          the normal DIS results
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC NOTE CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     onoff should be set to false using this routine as soon as you are
C     done calculating the offshell corrections.  It WILL interfere in
C     the rest of the code otherwise.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      Implicit None
      logical onoff

      logical offshell_on
      common /offshell/ offshell_on

      offshell_on = onoff

      ! Sets the needed fmKP coefficients 
      ! (Needed for the free-floating option, to make sure the current dRN
      ! is used - negligible overhead otherwise)
      call set_offshell_par 
         
      return
      end


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine set_offshell_par
C
C     This function sets all the possible constants we will need
C     It specifically sets the values of sbar for each parton,
C     lambda, and the normalization constants
C
C     It retrieves the nuclear model previously stored in a common block
C     and uses the following switches. (NOTE: it assumes the inuke code
C     was stored using 'set_nuke' in the calling routines.)
C
C        wfn - Integer - Denotes the wavefunction being used
C                    0 = Paris
C                    1 = AV18
C                    3 = CDBonn   
C                    4 = WJC-1  
C                    5 = WJC-2   
C
C        ioff - offshell model  4 = mKP, 5 = fmKP, 6 = "negative" fmKP
C                               7 = free floating fmKP (iDRN ignored)
C 
C        idRN - Integer - Denotes the size of the Nuclear Radius decrease
C
C                      mkP          fmKP      neg fmKP
C
C                    1 = 0.3%     1 = 0.1%    1 = -0.1% 
C                    2 = 0.6%     2 = 0.2%    2 = -0.2% 
C                    3 = 0.9%     3 = 0.3%    3 = -0.3% 
C                    4 = 1.2%     4 = 0.4%    4 = -0.4% 
C                    5 = 1.5%     5 = 0.5%    5 = -0.5% 
C                    6 = 1.8%     6 = 0.6%    6 = -0.6% 
C                    7 = 2.1%     7 = 0.7%    7 = -0.7% 
C                    8 = 2.4%     8 = 0.8%    8 = -0.8% 
C                    9 = 2.7%     9 = 0.9%    9 = -0.9% 
C                    0 = 3.0%     0 = 0%      0 =  0% 
c
      Implicit None
      integer inuke,ilam,ishad,idRN,ioff,iwfn,wfn,itgt
      double precision dp2_array(5)
      double precision dRN_array_mKP(10),dRN_array_fmKP(10)
      double precision offpar(14)
      double precision c_uv_array(5,10),c_dv_array(5,10)
      double precision c_g_array(5,10),c_s_array(5,10),c_v_array(5,10)
      double precision sv, ssq, sg, lambda
      common /S_block/ sv, ssq, sg, lambda
      data dp2_array/-3.60D-2, -4.26D-2, -6.21D-2, -4.86D-2, -4.29D-2/
      data dRN_array_mKP/0.3D-2, 0.6D-2, 0.9D-2, 1.2D-2, 1.5D-2, 1.8D-2,
     &                   2.1D-2, 2.4D-2, 2.7D-2, 3.0D-2/
      data dRN_array_fmKP/0.1D-2, 0.2D-2, 0.3D-2, 0.4D-2, 0.5D-2, 0.6D-2,
     &                    0.7D-2, 0.8D-2, 0.9D-2, 0D-2/

      double precision old_lambda
      data old_lambda/-2000d0/
      save old_lambda
      
      ! Retrieves and unfolds the nuclear model code
      call get_nuke(inuke)
      call split_nuke(inuke,ilam,ishad,idRN,ioff,wfn,itgt)

      ! if fmKP off-shell model requested, sets its parameters
      if (ioff.ge.5.and.ioff.le.7) then

         if (wfn.Eq.2) iwfn=1
         if (wfn.Eq.1) iwfn=2
         if (wfn.Eq.3) iwfn=3
         if (wfn.Eq.4) iwfn=4
         if (wfn.Eq.0) iwfn=5
         
         if (idRN.eq.0) then
            idRN = 10
         endif
         
         ssq = 5.5              ! 2.907 commented GRV, using GJR
         sg  = 8.0              ! 1.0
         sv =  2.2              ! 2.4
         if (ioff.eq.4) then    ! mKP off-shell parameter grid 
            lambda = -2*dRN_array_mKP(idRN)/dp2_array(iwfn)
         else if (ioff.eq.5) then ! finer grid for fmKP model
            lambda = -2*dRN_array_fmKP(idRN)/dp2_array(iwfn)
         else if (ioff.eq.6) then ! "negative" fmKP model (plus sign below)
            lambda = 2*dRN_array_fmKP(idRN)/dp2_array(iwfn)
         else if (ioff.eq.7) then ! fmKP with running offshell param
            call get_off(offpar)
            lambda = -2*offpar(1)/dp2_array(iwfn)
         else
            print*, 'ERROR(set_off): ioff out of range: ',ioff
            stop
         end if

         if (lambda.ne.old_lambda) then
            old_lambda=lambda
c            print*, '* -----------  New lambda =',lambda,offpar(1)
            if (isnan(offpar(1))) stop
         end if

      end if

      return 
      end




CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine init_offshell

C     Initializes the fmKP off-shell model. 
C     (Does no harm if called for other models)

      implicit none

      double precision c_v, c_sq, c_g
      common /C_block/ c_v, c_sq, c_g

      integer count
      data count/0/
      save count


      ! Switches off offshell effects as a precaution.
      !   (The user is always required to switch them on only when needed, 
      !   and off as done with off-shell PDF calculations as to not mess up
      !   possible subsequent on-shell PDF evaluations.)
      call set_offshell_on(.false.) 

      ! Sets the needed fmKP coefficients
      call set_offshell_par 

      ! calculuates needed PDF integrals, store in common block
      call mk_cv(c_v)
      call mk_cg_csq(c_v, c_g, c_sq)  

      count = count+1
C      print*, '* init_offshell call numebr', count, c_v,c_g,c_sq
      return
      end
