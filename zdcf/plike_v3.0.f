ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                             c
c                  Calculating the Peak Likelihood of the ZDCF                c
c                                                                             c
c        A special version for processing the ZDCF output file *.dcf          c
c                                                                             c
c                                Version 3.0                                  c
c                                                                             c
c    (Bayesian prior assumption of uniform distribution in z-space)           c
c                                                                             c
c                                                                             c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                             c
c SAMPLE RUN (user input marked by "<-USER". See explanations below)          c
c ----------                                                                  c
c                                                                             c
c PLIKE V3.0 begins.                                                          c
c Enter dcf file name:                                                        c
c test.dcf       <-USER                                                       c
c Enter lower bound on peak location:                                         c
c -200           <-USER                                                       c
c Enter upper bound on peak location:                                         c
c 200            <-USER                                                       c
c Calculating ML in range  -200.  to  200.                                    c
c                                                                             c
c ZDCF peak at  +0.133     r_max =  +0.881     , ZDCF C.O.M at   -4.31        c
c                                                                             c
c  #      lag         r         -dr        +dr    likelihood                  c
c --- ---------- ---------- ---------- ---------- ----------                  c
c  1  -76.8      0.355      0.314      0.275      2.701E-03                   c
c  2  -19.8      0.780      0.121      9.680E-02  0.111                       c
c  3  0.133      0.881      8.032E-02  6.053E-02  0.458                       c
c  4   5.62      0.821      0.124      9.322E-02  0.205                       c
c  5   24.1      0.392      0.306      0.265      3.669E-03                   c
c  6   87.5     -0.162      0.319      0.339      2.365E-05                   c
c  7   200.     -0.435      0.251      0.295      7.652E-07                   c
c                                                                             c
c ML Peak at  +0.133     , L at peak =  +0.458                                c
c 1 sigma ML interval =  +0.133       +10.0       -23.7                       c
c                     = (  -23.6     ,   +10.2    )                           c
c                                                                             c
c Program ended.                                                              c
c                                                                             c
c EXPLANATION                                                                 c
c -----------                                                                 c
c                                                                             c
c o The program uses the *.dcf output file of the ZDCF.                       c
c                                                                             c
c o The user has to give lower and upper bounds on the peak location          c
c                                                                             c
c o The output consists of the likelihood function and the fiducial interval, c
c   given both as t_peak + dt_up - dt_low, and as the interval (t_min, t_max) c
c                                                                             c
c FORTRAN ISSUES                                                              c
c --------------                                                              c
c                                                                             c
c o Compile with static variables flag:                                       c
c      Intel Fortran: ifort -save -o plike plike_v3.0.f                       c
c      g77 Fortran: g77 -fno-automatic  -o plike plike_v3.0.f                 c
c      g95 Fortran: g95 -fstatic  -o plike plike_v3.0.f                       c
c      Other compilers: -static (?)                                           c
c                                                                             c
c o The program can be VERY slow - don't panic!                               c
c                                                                             c
c o This program incorporates the functions QROMO, POLINT and MIDPNT          c
c   from Numerical Recipes / Press, Flannery, Teukolsky & Vetterling.         c
c                                                                             c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      integer i,np,maxp
      parameter (maxp = 200)
      integer n(maxp)
      real t(maxp),dtm(maxp),dtp(maxp),r(maxp),drm(maxp),drp(maxp)
      real liklhd(maxp),mxlik,tmxlik
      real tlow,thigh,dtlikl,dtliku,tcom
      character*72 dcf
      integer oerr,cerr
      real com
c     
      print *,'PLIKE V3.0 begins.'
c     
      print *,'Enter dcf file name:'
      read '(a72)',dcf
      open (unit=11,file=dcf,status='OLD',iostat=oerr,err=901)
c     
      print *,'Enter lower bound on peak location:'
      read *,tlow
      print *,'Enter upper bound on peak location:'
      read *,thigh
      print *,'Calculating ML in range ',tlow,' to ', thigh
      np = 0
1     continue
         if (np .ge. maxp) then
            print *,'Only 1st ',maxp,' points were read!'
            goto 2
         endif 
         read (11,*,end=2) 
     >      t(np+1),dtm(np+1),dtp(np+1),
     >      r(np+1),drm(np+1),drp(np+1),n(np+1)
         if (t(np+1) .lt. tlow .or. t(np+1) .gt. thigh) goto 1
         np = np+1
         goto 1
2     continue
      close (unit=11,status='KEEP',iostat=cerr,err=902)
c     
      print *
      tcom = com(t,r,np,i)
      print '(''ZDCF peak at '',sp,1p,g11.3,'' r_max = '',g11.3,
     >        '' , ZDCF C.O.M at '',g11.3)',t(i),r(i),tcom
      call plike3(t,dtm,dtp,r,drm,drp,n,np,tlow,thigh,
     >            liklhd,tmxlik,mxlik,dtlikl,dtliku,2)
c     
      print *
      print *,
     >' #      lag         r         -dr        +dr    likelihood '
      print *,
     >'--- ---------- ---------- ---------- ---------- ---------- '
      print '(i3,'' '',1p,g10.3,1x,g10.3,1x,g10.3,1x,g10.3,1x,g10.3)',
     >   (i,t(i),r(i),drm(i),drp(i),liklhd(i),i=1,np)
c     
      print *
      print '(''ML Peak at '',sp,1p,g11.3,'' , L at peak = '',g11.3)',
     >   tmxlik,mxlik
      print '(''1 sigma ML interval = '',sp,1p,g11.3,1x,g11.3,1x,
     >        g11.3,/,19x,'' = ('',g11.3,'' , '',g11.3,'')'')',
     >   tmxlik,dtliku,-dtlikl,tmxlik-dtlikl,tmxlik+dtliku
      print *
      print *,'Program ended.'
      stop
 901  print *,'ERROR ON OPENING FILE - ERRCOD=',oerr
      stop
 902  print *,'ERROR ON CLOSING FILE - ERRCOD=',cerr
      stop
      end
c//////////////////////////////////////////////////////////////////////////////
c     
c com calculates the center-of-mass (centroid) around the ZDCF peak.
c The sum is taken over all the points lying between the two nearest 
c the maximum with r < r_max/2.
c     
      real function com(t,r,n,imax)
      implicit none
      integer n,imax,nerr
      real t(n),r(n)
      real sumr,sumtr,rhalf
      integer i,i1,i2
c Locating the maximum
      imax = 1
      do i = 1,n
         if (r(imax) .lt. r(i)) imax = i
      enddo
c     
      rhalf = r(imax)/2.
c Locating the peak
      i1 = imax-1
      do while (i1.ge.1 .and. r(i1).ge.rhalf)
         i1 = i1-1
      enddo
      i1 = i1+1
c Warning 
      if (i1.eq.1 .and. r(1).ge.rhalf) then
         print *,
     >      'COM: Warning - Peak extends beyond range to the left...'
      endif
c     
      i2 = imax+1
      do while (i2.le.n .and. r(i2).ge.rhalf)
         i2 = i2+1
      enddo
      i2 = i2-1
c Warning
      if (i2.eq.n .and. r(n).ge.rhalf) then
         print *,
     >      'COM: Warning - Peak extends beyond range to the right...'
         nerr = nerr+1
      endif 
c Calculating the center of mass
      sumtr = 0.0
      sumr = 0.0
      do i = i1,i2
         sumtr = sumtr+t(i)*r(i)
         sumr = sumr+r(i)
      enddo
      if (sumr .le. 0.0) then
         print *,'COM: Problems. n = ',n,' imax = ',imax,
     >           ' rmax = ',r(imax),
     >           ' i1 = ',i1,' i2 = ',i2,' sumr = ',sumr
         print '(1x,f10.3,1x,f10.3)',(t(i),r(i),i=1,n)
         stop      
      else
         com = sumtr/sumr
      endif
      return
      end
c//////////////////////////////////////////////////////////////////////////////
cold      include "/home/tal/ZDCF/src/plike_func3.f"
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                             c
c                  Calculating the Peak Likelihood of the ZDCF                c
c                                                                             c
c        A special version for processing the ZDCF output file *.dcf          c
c                                                                             c
c                                Version 3.0                                  c
c                                                                             c
c    (Bayesian prior assumption of uniform distribution in z-space)           c
c                                                                             c
c                                                                             c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                             c
c SAMPLE RUN (user input marked by "<-USER". See explanations below)          c
c ----------                                                                  c
c                                                                             c
c PLIKE V3.0 begins.                                                          c
c Enter dcf file name:                                                        c
c test.dcf      <-USER                                                        c
c Enter lower bound on peak location:                                         c
c -200           <-USER                                                       c
c Enter upper bound on peak location:                                         c
c 200            <-USER                                                       c
c Working on   7 points                                                       c
c                                                                             c
c   i       tau        z          s_z       likelihood                        c
c -----------------------------------------------------                       c
c    1  -76.77000    0.39168    0.35056     3.2769E-03                        c
c    2  -19.79000    1.07518    0.28485     0.1270                            c
c    3    0.13330    1.42651    0.32573     0.5923                            c
c    4    5.61800    1.20773    0.34551     0.2729                            c
c    5   24.10000    0.43595    0.35034     4.4965E-03                        c
c    6   87.48000   -0.17234    0.35143     2.2088E-05                        c
c    7  199.89999   -0.49122    0.35009     4.8453E-07                        c
c                                                                             c
c Peak at 0.133     , ML(peak) = 0.592                                        c
c 1 sigma ML interval =  0.133    + 13.4    -  23.9                           c
c                     = ( -23.8     ,   13.5    )                             c
c                                                                             c
c Program ended.                                                              c
c                                                                             c
c EXPLANATION                                                                 c
c -----------                                                                 c
c                                                                             c
c o The program uses the *.dcf output file of the ZDCF.                       c
c                                                                             c
c o The user has to give lower and upper bounds on the peak location          c
c                                                                             c
c o The output consists of the likelihood function and the fiducial interval, c
c   given both as t_peak + dt_up - dt_low, and as the interval (t_min, t_max) c
c                                                                             c
c FORTRAN ISSUES                                                              c
c --------------                                                              c
c                                                                             c
c o The program is VERY slow - don't panic!                                   c
c                                                                             c
c o This program incorporates the functions QROMO, POLINT and MIDPNT          c
c   from Numerical Recipes / Press, Flannery, Teukolsky & Vetterling.         c
c                                                                             c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine plike3(t,dtm,dtp,r,drm,drp,n,np0,tlow,thigh,
c                      i i   i   i i  i  i i i   i    i    
     >                 liklhd,tmxlik,mxlik,dtlikl,dtliku,zrflag)
c                      o      o      o     o      o      i     
      implicit none
      integer np0,zrflag
      real tlow,thigh,tmxlik,mxlik,dtlikl,dtliku
      real t(np0),dtm(np0),dtp(np0),r(np0),drm(np0),drp(np0),liklhd(np0)
      integer n(np0)
c     
      integer i,maxp,np,ilow
      parameter (maxp = 200)
      integer pnt(maxp)
      real z(maxp),pcum(maxp)
      real p,rpeak
      integer imxlik,ipeak
      real eps,epslik,tiny,maxsig
      parameter (eps = 1e-6, epslik = 1e-10, tiny = 1e-30, maxsig = 5)
      real likint,likin0,likin1
      external likint,midpt1
c     
      np = 0
      ilow = 0
      rpeak = -1.0
c     Input given in z-space (zrflag = 1): vector r holds z
      if (zrflag .eq. 1) then
         do i = 1,np0
            if (t(i) .ge. tlow .and. t(i) .le. thigh) then
               np = np+1
               if (np .gt. maxp) then
                  print *,'PLIKE: Work area too small ! ',
     >                 'Increase maxp to >= ',maxp
                  stop
               endif
               z(np) = r(i)
               liklhd(np) = 0.0
               if (r(i) .gt. rpeak) then
                  ipeak = np
                  rpeak = r(i)
               endif
            else if (t(i) .lt. tlow) then 
               ilow = i
            endif
         enddo
c     Input given in r-space (rzflag <> 1)
      else
         do i = 1,np0
            if (t(i) .ge. tlow .and. t(i) .le. thigh) then
               np = np+1
               if (np .gt. maxp) then
                  print *,'PLIKE: Work area too small ! ',
     >                 'Increase maxp to >= ',maxp
                  stop
               endif
               z(np)  = log((1.+r(np))/(1.-r(np)))/2
               liklhd(np) = 0.0
               if (r(i) .gt. rpeak) then
                  ipeak = np
                  rpeak = r(i)
               endif
            else if (t(i) .lt. tlow) then 
               ilow = i
            endif
         enddo
      endif
c Creating a vector of pointers to cells in t, dtm, dtp which relate to z
      do i = 1,np
         pnt(i) = ilow+i
      enddo
      ilow = ilow+1
c     
c Calculating the likelihood function
c     
      p = likin0(z,n,np)
      do i = ipeak-1,1,-1
         p = likin1(i)
         print '(''.'',$)'
c Note small modifications of original qromo
         call qromo1(likint,-1d0,1d0,liklhd(i),midpt1)
cdbg( 
c         print '(1x,i3,1x,0p,f10.5,1x,f10.5,5x,i3,7x,1p,g10.4)',
c     >         i,t(i),z(i),n(i),liklhd(i)
cdbg) 
         if (liklhd(i) .lt. epslik) goto 301
      enddo
 301  continue 
      do i = ipeak,np
         p = likin1(i)
         print '(''.'',$)'
c Note small modifications of original qromo
         call qromo1(likint,-1d0,1d0,liklhd(i),midpt1)
cdbg( 
c         print '(1x,i3,1x,0p,f10.5,1x,f10.5,5x,i3,7x,1p,g10.4)',
c     >         i,t(i),z(i),n(i),liklhd(i)
cdbg) 
         if (liklhd(i) .lt. epslik) goto 302
      enddo
 302  continue 
      print '(/)'
c     
c Finding the maximal likelihood value and error interval
c     
      call conf1(liklhd,t(ilow),dtm(ilow),dtp(ilow),np,
     >           pcum,imxlik,tmxlik,mxlik,dtlikl,dtliku) 
c     
cdbg      print *
cdbg      print '(''Peak at '',1p,g9.3,'' , ML(peak) = '',g9.3)',
cdbg     >   tmxlik,mxlik
cdbg      print '(''1 sigma ML interval = '',1p,g10.3,''+'',g9.3,''- '',g9.3,
cdbg     >        /,19x,'' = ('',g10.3,'' , '',g10.3,'')'')',
cdbg     >   tmxlik,t1sigu-tmxlik,tmxlik-t1sigl,t1sigl,t1sigu
c     
      return 
      end
c//////////////////////////////////////////////////////////////////////////////
c     
c The integrand of the likelyhood function
c     
      real function likint(rho)
      implicit none
      integer maxp
      parameter (maxp = 200)
      real*8 rho
      real z0(*),z(maxp),cn,zbar,sz,dzdr
      integer n0(*),n(maxp)
      integer i0,i,np0,np
c cpi = 1/sqrt(2*Pi)
      integer j
      real cpi,tiny
      parameter (cpi = 0.398942280402, tiny = 1e-10)
      real clike,exp,log
      real likin0,likin1
c     
cdbg( 
      if (abs(rho) .ge. 1.0) then
	print *,'LIKINT: |rho| >= 1'
	stop
      endif
cdbg) 
      zbar = log((1.+rho)/(1.-rho))/2.+rho/2./(n(i)-1)*
     >         (1.+(5.+rho**2)/4./(n(i)-1.)+
     >             (11.+2.*rho**2+3.*rho**4)/8./(n(i)-1)**2) 
      sz = sqrt(1./(n(i)-1)*(1.+(4-rho**2)/2./(n(i)-1)+
     >                       (22.-6*rho**2-3*rho**4)/6./(n(i)-1)**2))
      dzdr = 1./(1.+rho)/(1.-rho)
      likint = log(dzdr)+log(cpi)-log(sz)-((z(i)-zbar)/sz)**2/2
      do j = 1,i-1
         cn = clike(rho,z(j),n(j))
         if (cn .le. tiny) then
            likint = 0.0
            return
         endif
         likint = likint+log(cn)
      enddo
      do j = i+1,np
         cn = clike(rho,z(j),n(j))
         if (cn .le. tiny) then
            likint = 0.0
            return
         endif
         likint = likint+log(cn)
      enddo
      likint = exp(likint)
      return
c     
      entry likin0(z0,n0,np0)
      if (np .gt. maxp) then
         print *,'LIKINT: maxp too small. np = ',np
         stop
      endif
      np = np0
      do j = 1,np
         z(j) = z0(j)
         n(j) = n0(j)
      enddo
      return
c     
      entry likin1(i0)
      i = i0
      return
      end
c///////////////////////////////////////////////////////////////////////
c Integrating over the likelihood function from rho=-1 to rho=rho0 for a
c given value of the estimator r evaluated on n pairs.
      real function clike(rho0,z,n)
      implicit none
      real*8 rho0
      real z
      integer n
c     
      real phiz0,phiz1
      external midpt2,phiz1
c     
      clike = phiz0(n,z)
cdbg( 
      if (abs(rho0) .gt. 1.0) then
	print *,'CLIKE: |rho| > 1'
	stop
      endif
cdbg) 
	
c Note small modifications of original qromo
      call qromo2(phiz1,-1d0,rho0,clike,midpt2)
      return 
      end
c///////////////////////////////////////////////////////////////////////
c The integrand of the cumulative likelihood function
c exp(-(z(r)-zbar(rho))/sz(rho))**2/2)/sqrt(2Pi)/sz(rho)
      real function phiz1(rho)
      implicit none
      integer n0,n
      real z0,z
      real*8 rho
      real dzdr,zbar,sz
c sqrpi2 = 1/sqrt(2*Pi)
      real sqrpi2
      parameter (sqrpi2 = .39894228040143267794)
      real sqrt,exp,phiz0
c     
cdbg( 
      if (abs(rho) .ge. 1.0) then
	print *,'PHIZ1: |rho| >= 1'
	stop
      endif
cdbg) 
      zbar = log((1.+rho)/(1.-rho))/2.+rho/2./(n-1)*
     >         (1.+(5.+rho**2)/4./(n-1.)+
     >             (11.+2.*rho**2+3.*rho**4)/8./(n-1)**2) 
      sz = sqrt(1./(n-1)*(1.+(4-rho**2)/2./(n-1)+
     >                       (22.-6*rho**2-3*rho**4)/6./(n-1)**2))
      dzdr = 1./(1.+rho)/(1.-rho)
      phiz1 = dzdr*sqrpi2/sz*exp(-((z-zbar)/sz)**2/2)
cdbg      phiz1 = log(dzdr*sqrpi2)-log(sz)-((z-zbar)/sz)**2/2
cdbg      phiz1 = exp(phiz1)
      return 
c     
      entry phiz0(n0,z0)
      n = n0
      z = z0
      phiz0 = 0.0
      return
      end
c//////////////////////////////////////////////////////////////////////////////
c Finding the +/- 1 sigma intervals on the interpolated likelihood function
c which (after normalization) is the "Fiducial Probability Function" (Kendall 
c & Stuart, Vol 2, 3rd ed. p. 141) 
      subroutine conf1(liklhd,t,dtm,dtp,np,pcum,imxlik,tmxlik,mxlik,
c                      i      i i   i   i  o    o      o      o
     >                 dtlikl,dtliku)
c                      o      o
      implicit none
      integer np,imxlik
      real liklhd(np),t(np),dtm(np),dtp(np),pcum(np)
      real tmxlik,mxlik,dtlikl,dtliku
      integer i
      real cn,s,p,t1sigl,t1sigu
c     
      mxlik = -1.0
      do i = 1,np
         if (mxlik .lt. liklhd(i)) then
            tmxlik = t(i)
            imxlik = i
            mxlik = liklhd(i)
         endif
      enddo
      cn = 0.0
      pcum(1) = 0.0
      do i = 2,np
         cn = cn+(liklhd(i)+liklhd(i-1))/2.0*
     >           (t(i)-t(i-1))
         pcum(i) = cn
      enddo
      do i = 1,np
         pcum(i) = pcum(i)/cn
      enddo
c     
      p = pcum(imxlik)+0.3413*(pcum(np)-pcum(imxlik))*2
      do i = imxlik,np
         if (pcum(i) .ge. p) then
            s = (liklhd(i)-liklhd(i-1))/(t(i)-t(i-1))
            t1sigu = t(i-1)+
     >             (-liklhd(i-1)+
     >              sqrt(liklhd(i-1)**2+2*(p-pcum(i-1))*s*cn))/s
            goto 100
         endif
      enddo
cdbg  
      print *,'t1sigu = end point...'
      t1sigu = t(np)
 100  continue
c     
      p = pcum(imxlik)-0.3413*(pcum(imxlik)-pcum(1))*2
      do i = imxlik,1,-1
         if (pcum(i) .le. p) then
            s = (liklhd(i+1)-liklhd(i))/(t(i+1)-t(i))
            t1sigl = t(i)+
     >             (-liklhd(i)+
     >              sqrt(liklhd(i)**2+2*(p-pcum(i))*s*cn))/s
            goto 200
         endif
      enddo
cdbg  
      print *,'t1sigl = end point...'
      t1sigl = t(1)
 200  continue
c     
c???      dtliku = amax0(t1sigu-tmxlik,dtp(imxlik))
c???      dtlikl = amax0(tmxlik-t1sigl,dtm(imxlik))
      dtliku = amax1(t1sigu-tmxlik,dtp(imxlik))
      dtlikl = amax1(tmxlik-t1sigl,dtm(imxlik))
cdbg  
cdbg      print *,'CONF1: t1sigl = ',t1sigl,' t1sigu = ',t1sigu
cdbg  
      return 
      end
c///////////////////////////////////////////////////////////////////////
      SUBROUTINE QROMO1(FUNC,A,B,SS,CHOOSE)
      save all
      PARAMETER (EPS=1.E-6,JMAX=14,JMAXP=JMAX+1,KM=4,K=KM+1)
      DIMENSION S(JMAXP),H(JMAXP)
cTal( 
      real*8 a,b
cTal)
      H(1)=1.
      DO 11 J=1,JMAX
        CALL CHOOSE(FUNC,A,B,S(J),J)
        IF (J.GE.K) THEN
          CALL POLIN1(H(J-KM),S(J-KM),K,0.0,SS,DSS)
cTal          IF (ABS(DSS).LT.EPS*ABS(SS)) RETURN
          IF (ABS(DSS).LE.EPS*ABS(SS)) RETURN
        ENDIF
        S(J+1)=S(J)
        H(J+1)=H(J)/9.
11    CONTINUE
      print *,'QROMO1: Too many steps...'
      END
c///////////////////////////////////////////////////////////////////////
      SUBROUTINE MIDPT1(FUNC,A,B,S,N)
cTal( 
      save all
      real*8 a,b,tnm,del,ddel,x,sum
cTal) 
      IF (N.EQ.1) THEN
        S=(B-A)*FUNC(0.5*(A+B))
        IT=1
      ELSE
        TNM=IT
        IF(TNM.EQ.0.)PAUSE 'MIDPT1: attempting to divide by den = 0'
        DEL=(B-A)/(3.*TNM)
        DDEL=DEL+DEL
        X=A+0.5*DEL
        SUM=0.
        DO 11 J=1,IT
          SUM=SUM+FUNC(X)
          X=X+DDEL
          SUM=SUM+FUNC(X)
          X=X+DEL
11      CONTINUE
        S=(S+(B-A)*SUM/TNM)/3.
        IT=3*IT
      ENDIF
      RETURN
      END
c//////////////////////////////////////////////////////////////////////////////
      SUBROUTINE POLIN1(XA,YA,N,X,Y,DY)
      save all
      PARAMETER (NMAX=10) 
      DIMENSION XA(N),YA(N),C(NMAX),D(NMAX)
      NS=1
      DIF=ABS(X-XA(1))
      DO 11 I=1,N 
        DIFT=ABS(X-XA(I))
        IF (DIFT.LT.DIF) THEN
          NS=I
          DIF=DIFT
        ENDIF
        C(I)=YA(I)
        D(I)=YA(I)
11    CONTINUE
      Y=YA(NS)
      NS=NS-1
      DO 13 M=1,N-1
        DO 12 I=1,N-M
          HO=XA(I)-X
          HP=XA(I+M)-X
          W=C(I+1)-D(I)
          DEN=HO-HP
          IF(DEN.EQ.0.)PAUSE 'POLIN1: attempting to divide by den = 0'
          DEN=W/DEN
          D(I)=HP*DEN
          C(I)=HO*DEN
12      CONTINUE
        IF (2*NS.LT.N-M)THEN
          DY=C(NS+1)
        ELSE
          DY=D(NS)
          NS=NS-1
        ENDIF
        Y=Y+DY
13    CONTINUE
      RETURN
      END
c///////////////////////////////////////////////////////////////////////
      SUBROUTINE QROMO2(FUNC,A,B,SS,CHOOSE)
      save all
cTal      PARAMETER (EPS=1.E-6,JMAX=14,JMAXP=JMAX+1,KM=4,K=KM+1)
      PARAMETER (EPS=1.E-6,JMAX=15,JMAXP=JMAX+1,KM=4,K=KM+1)
      DIMENSION S(JMAXP),H(JMAXP)
cTal( 
      real*8 a,b
cTal) 
      H(1)=1.
      DO 11 J=1,JMAX
        CALL CHOOSE(FUNC,A,B,S(J),J)
        IF (J.GE.K) THEN
          CALL POLIN2(H(J-KM),S(J-KM),K,0.0,SS,DSS)
cTal          IF (ABS(DSS).LT.EPS*ABS(SS)) RETURN
          IF (ABS(DSS).LE.EPS*ABS(SS)) RETURN
        ENDIF
        S(J+1)=S(J)
        H(J+1)=H(J)/9.
11    CONTINUE
      print *,'QROMO2: Too many steps...'
      END
c///////////////////////////////////////////////////////////////////////
      SUBROUTINE MIDPT2(FUNC,A,B,S,N)
cTal( 
      save all
      real*8 a,b,tnm,del,ddel,x,sum
cTal) 
      IF (N.EQ.1) THEN
        S=(B-A)*FUNC(0.5*(A+B))
        IT=1
      ELSE
        TNM=IT
        IF(TNM.EQ.0.)PAUSE 'MIDPT2: attempting to divide by den = 0'
        DEL=(B-A)/(3.*TNM)
        DDEL=DEL+DEL
        X=A+0.5*DEL
        SUM=0.
        DO 11 J=1,IT
          SUM=SUM+FUNC(X)
          X=X+DDEL
          SUM=SUM+FUNC(X)
          X=X+DEL
11      CONTINUE
        S=(S+(B-A)*SUM/TNM)/3.
        IT=3*IT
      ENDIF
      RETURN
      END
c//////////////////////////////////////////////////////////////////////////////
      SUBROUTINE POLIN2(XA,YA,N,X,Y,DY)
      save all
      PARAMETER (NMAX=10) 
      DIMENSION XA(N),YA(N),C(NMAX),D(NMAX)
      NS=1
      DIF=ABS(X-XA(1))
      DO 11 I=1,N 
        DIFT=ABS(X-XA(I))
        IF (DIFT.LT.DIF) THEN
          NS=I
          DIF=DIFT
        ENDIF
        C(I)=YA(I)
        D(I)=YA(I)
11    CONTINUE
      Y=YA(NS)
      NS=NS-1
      DO 13 M=1,N-1
        DO 12 I=1,N-M
          HO=XA(I)-X
          HP=XA(I+M)-X
          W=C(I+1)-D(I)
          DEN=HO-HP
          IF(DEN.EQ.0.)PAUSE 'POLIN2: attempting to divide by den = 0'
          DEN=W/DEN
          D(I)=HP*DEN
          C(I)=HO*DEN
12      CONTINUE
        IF (2*NS.LT.N-M)THEN
          DY=C(NS+1)
        ELSE
          DY=D(NS)
          NS=NS-1
        ENDIF
        Y=Y+DY
13    CONTINUE
      RETURN
      END
c//////////////////////////////////////////////////////////////////////////////
