C***********************************************************************
C    Module:  sgutil.f
C 
C    Copyright (C) 2020 Mark Drela, Harold Youngren
C 
C    This program is free software; you can redistribute it and/or modify
C    it under the terms of the GNU General Public License as published by
C    the Free Software Foundation; either version 2 of the License, or
C    (at your option) any later version.
C
C    This program is distributed in the hope that it will be useful,
C    but WITHOUT ANY WARRANTY; without even the implied warranty of
C    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C    GNU General Public License for more details.
C
C    You should have received a copy of the GNU General Public License
C    along with this program; if not, write to the Free Software
C    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
C***********************************************************************

      subroutine akima ( x, y, n, xx, yy, slp )
c---------------------------------------------------------------------
c     General-purpose monovariate interpolation routine
c     using a locally-fitted cubic spline. One point is
c     interpolated from the input data arrays.
c
c inputs:
c     x(.)  array of abscissas in ascending or descending value order
c     y(.)  array of corresponding ordinates. y(x) must be single-valued.
C     n     size of the x,y arrays
c     xx    x point at which an interpolated y value is desired.
C
C Outputs:
c     yy    interpolated y ordinate
c     slp   interpolated slope dy/dx
C
C...DISCUSSION This spline method produces an interpolated curve that
C              is relatively free of the oscillatory behavior normally
C              associated with cubic splines. The curve produced will
C              be continuous in 'Y' and 'DY/DX' but further derivatives
C              may be discontinuous. The interpolated slope should be
C              treated with some caution, however, due to the tendency
C              of this spline curve to concentrate changes of curvature
C              at the input points (interval ends). If only 2 points are
C              input the curve will be linear, if 3 points are given
C              the curve will be parabolic, more than 3 points will
C              produce a cubic. Extrapolation beyond the bounds of the
C              input data is done with a quadratic through the last 3
C              points at that data boundary.
C
C              This routine is intended as a replacement for B.Wainfan's
C              AKIMAD. It is shorter, twice as fast and it works on
C              input data in ascending or descending order. The calling
C              sequence is more natural and is easier to use for the
C              bulk of applications.
C
C              The coding is compatible with IBM or VAX Fortran 77 
C              and IBM Fortran IV G or H with optimization.
C
C...ORIGIN     Harold Youngren  CALAC Dept. 72-71   3/81
C
C...REFERENCE  Akima, Hiroshi, "A New Method of Interpolation and
C              Smooth Curve Fitting Based on Local Procedures",
C              Journal of the Association for Computing Machines,
C              Vol. 17, No. 4, Oct 1970, pages 589-602
C
C              Wainfan, B. S., "The Akima Subroutines:  Nonlinear
C              Interpolation by Local Polynomial Fit", LR 29244,
C              Oct 30, 1979.
C---------------------------------------------------------------------
      real x(n), y(n), d(5), t(2)
c 
c
c...check for a degenerate case ( x(1)=x(n) ).    
c 
      if (x(1) .ne. x(n))  go to 10
         yy  = y(1) 
         slp = 0. 
         go to 70
c
c
c...find which interval contains the point by binary search.
c...the binary search loop is terminated when the search residual
c...(nstep) is zero. the index 'i' will point to the input point
c...lower than or equal to the desired point for the 'x' values
c...in ascending order or to the input point greater than or equal
c...to the desired point for descending order 'x' values.       
c
   10 xordr = 1.0
      if (x(1) .gt. x(n))  xordr = -1.0
c
      ibot = 1
      itop = n
      xxo  = xx * xordr
c
   20 nstep = ( itop - ibot ) / 2
      i     = ibot + nstep
      xo    = x(i) * xordr
      if ( xxo .ge. xo )  ibot = i
      if ( xxo .lt. xo )  itop = i
      if ( nstep .ne. 0 )  go to 20
c
c
c...calculate the straight line slopes between adjacent input points.
c...d(3) is the slope on the interval of interpolation. if the other
c...slopes d(1), d(2), d(4) or d(5) are not defined they will be
c...created by quadratic extrapolation (only at start and end of data).
c
      do 30  j = 1, 5
         k  = i + (j-2)
         if ( ((k-1) .ge. 1) .and. (k .le. n) )
     &      d(j) = ( y(k) - y(k-1) ) / ( x(k) - x(k-1) )
   30 continue
c
c...synthesize upper and lower slopes if required. check for
c...single line segment input (n=2).
c
      if (n .eq. 2)  d(2) = d(3)
c
      if ((i+2) .gt. n)  d(4) = 2. * d(3) - d(2)
      if ((i+3) .gt. n)  d(5) = 2. * d(4) - d(3)
      if ((i-1) .lt. 1)  d(2) = 2. * d(3) - d(4)
      if ((i-2) .lt. 1)  d(1) = 2. * d(2) - d(3)
c
c
c...calculate the slopes (t(1),t(2)) at the lower and upper
c...points bounding the interval of interpolation. if the point is
c...at an intersection of straight line segments the slope is
c...defined by the average of the adjacent segment slopes.
c
      do 50 j = 1, 2
         a = abs( d(j+3)  -  d(j+2) )
         b = abs( d(j+1)  -  d(j)   )
         if ((a + b) .ne. 0.)  go to 40
            a = 1.
            b = 1.
   40    t(j) = ( a*d(j+1) + b*d(j+2) ) / ( a + b )
   50 continue
c
c
c...check if desired point is on upper point of interval. this
c...reduces error at the transition points between intervals.
c
      if (xx .ne. x(i+1))  go to 60
         yy  = y(i+1)
         slp = t(2)
         go to 70
c
c...calculate the cubic coefficients.
c
   60 xint  =  x(i+1) - x(i)
      xdif  =  xx     - x(i)
      p0    =  y(i)
      p1    =  t(1)
      p2    =  ( 3.*d(3) - 2.*t(1) - t(2) ) / xint
      p3    =  ( t(1)  + t(2)  -  2.*d(3) ) / (xint*xint)
c
c...calculate the y-value and the slope.
c
      yy  =  p0 + xdif*( p1 + xdif*( p2 + xdif*p3 ) )
      slp =  p1 + xdif*( 2.*p2 + xdif*( 3.*p3 ) )
c
   70 return
      end ! akima



      function trp1 (n,x,y,xtrp)
c
c...purpose  to linearly interpolate a value from an 
c            array of data
c
c...input    n         number of points in array 
c            x,y       arrays of abscissae and ordinates
c            xtrp      coordinate at which ytrp is desired 
c
c...output   trp1      interpolated value
c
c...comments   
c
      real x(n), y(n)
c
      if (n.lt.1)  then
        trp1 = 0.
        return
      endif
      if (n.lt.2)  then
        trp1 = y(1)
        return
      endif
c
c...find the interval containing the point
      i = 1
   10 if (x(i+1).gt.xtrp .or. i+1.eq.n)  go to 20 
        i = i + 1
        go to 10
c
   20 trp1 = y(i) + (y(i+1)-y(i))*(xtrp-x(i))/(x(i+1)-x(i))
      return
      end ! trp1


      subroutine nrmliz (n,x)
c
c...purpose  to normalize an array of data
c
c...input    n         number of points in array 
c            x         array of data in ascending or descending order
c
c...output   x         normalized array (0<=x<=1.)
c
c...comments   
c
      real x(n)
c
      if (n.le.1)  return
c
      dx = x(n) - x(1)
      if (dx .eq. 0.0) dx = 1.0
c
      x1 = x(1)
      do 10 i = 1, n
        x(i) = (x(i)-x1) / dx
   10 continue
      return
      end ! nrmliz


      subroutine sspacer (n,pspace,x)
c-------------------------------------------------------------
c     Generates spanwise-spacing arrays.
C     Divides 0..1 interval into N segements with variable
c     spacing specified by the pspace parameter
c
c Inputs:  n       Number of points in array
c          pspace  Spacing parameter
c                  = 0  : equal spacing
c                  = 1  : cosine spacing.
c                  = 2  : sine spacing (concentrating points near 0)
c                  = 3  : equal spacing.
C
c                 negative values of pspace produce spacing
c                 which is reversed (affects only sine spacing).
c                 intermediate values of pspace will produce
c                 a spacing which is a linear combination
c                 of the corresponding integer values.
c     
c Output:  x(.)   Normalized spacing array 0 .. 1
c                 x(1) = 0.0
c                 x(n) = 1.0
C-------------------------------------------------------------
      real  x(n)

      pi = 4.0*atan(1.0)
c
      pabs = abs (pspace)
      nabs = ifix (pabs) + 1
c
      go to (10,20,30,30), nabs
c
   10   pequ = 1.0 - pabs
        pcos = pabs
        psin = 0.
      go to 50
c
   20   pequ = 0.
        pcos = 2.0 - pabs
        psin = pabs - 1.0
        go to 50
c
   30   pequ = pabs - 2.0
        pcos = 0.
        psin = 3.0 - pabs
c
   50 continue
      do 100  k = 1, n
         frac = float(k-1)/float(n-1)
         theta =  frac * pi
         if (pspace .ge. 0.0) then
          x(k) = pequ * frac
     &         + pcos * ( 1.0 - cos (theta)     ) / 2.0
     &         + psin * ( 1.0 - cos (theta/2.0) )
         else
           x(k) = pequ * frac
     &          + pcos * ( 1.0 - cos (theta) ) / 2.0
     &          + psin * sin (theta/2.0)
         endif
  100 continue
c
      return
      end ! sspacer



      subroutine cspacer(n,cspace,claf, xpt,xvr,xsr,xcp)
c---------------------------------------------------------------
c     Divides 0..1 interval into 4*n + 1 subintervals with
c     blended uniform, sine, or cosine spacing.
c
c Inputs:  n       Number of intervals in arrays
c          cspace  Spacing parameter
c                  = 0.0  : equal spacing
c                  = 1.0  : cosine spacing.
c                  = 2.0  : sine spacing (concentrating points near 0)
c                  = 3.0  : equal spacing.
c
c                  Negative values of pspace produce spacing
c                  which is reversed (affects only sine spacing).
c                  Intermediate values of pspace will produce
c                  a spacing which is a linear combination
c                  of the corresponding integer values.
c
c Output:  xpt(.)  points at 0,1 sub-interval points
c          xvr(.)  points at 1/4 sub-interval points
c          xsr(.)  points at 1/2 sub-interval points
c          xcp(.)  points at 3/4 sub-interval points
c
c          will always return xpt(1) = 0.0, xpt(n+1) = 1.0
c
c     For the case n=2, cspace=0.0 (uniform), the returned points are:
c          xpt(1) = 0.0
c          xvr(1) = 0.125
c          xsr(1) = 0.25
c          xcp(1) = 0.375
c          xpt(2) = 0.5
c          xvr(2) = 0.625
c          xsr(2) = 0.75
c          xcp(2) = 0.875
c          xpt(3) = 1.0
c
C---------------------------------------------------------------
      real xpt(*), xvr(*), xsr(*), xcp(*)
c
      pi = 4.0*atan(1.0)
c
c---- set blending weights
      acsp = abs(cspace)
      ncsp = ifix(acsp)
      if    (ncsp.eq.0) then
       f0 = 1.0 - acsp
       f1 = acsp
       f2 = 0.
      elseif(ncsp.eq.1) then
       f0 = 0.
       f1 = 2.0 - acsp
       f2 = acsp - 1.0
      else
       f0 = acsp - 2.0
       f1 = 0.
       f2 = 3.0 - acsp
      endif
c
c---- cosine chordwise spacing
      dth1 =     pi/float(4*n + 2)
      dth2 = 0.5*pi/float(4*n + 1)
      dxc0 =    1.0/float(4*n)
c
      do i = 1, n
c------ uniform
        xc0 = int(4*i - 4) * dxc0
        xpt0 = xc0        
        xvr0 = xc0 +     dxc0
        xsr0 = xc0 + 2.0*dxc0
        xcp0 = xc0 +     dxc0 + 2.0*dxc0*claf
c
c------ cosine
        th1 = int(4*i - 3) * dth1
        xpt1 = 0.5*(1.0 - cos(th1         ))
        xvr1 = 0.5*(1.0 - cos(th1+    dth1))
        xsr1 = 0.5*(1.0 - cos(th1+2.0*dth1))
        xcp1 = 0.5*(1.0 - cos(th1+    dth1+2.0*dth1*claf))
c
        if(cspace .gt. 0.0) then
c------- sine
         th2 = int(4*i - 3) * dth2
         xpt2 = 1.0 - cos(th2         )
         xvr2 = 1.0 - cos(th2+    dth2)
         xsr2 = 1.0 - cos(th2+2.0*dth2)
         xcp2 = 1.0 - cos(th2+    dth2+2.0*dth2*claf)
        else
c------- -sine
         th2 = int(4*i - 4) * dth2
         xpt2 = sin(th2         )
         xvr2 = sin(th2+    dth2)
         xsr2 = sin(th2+2.0*dth2)
         xcp2 = sin(th2+    dth2+2.0*dth2*claf)
        endif
c
c------ blend 'em
        xpt(i) = f0*xpt0 + f1*xpt1 + f2*xpt2
        xvr(i) = f0*xvr0 + f1*xvr1 + f2*xvr2
        xsr(i) = f0*xsr0 + f1*xsr1 + f2*xsr2
        xcp(i) = f0*xcp0 + f1*xcp1 + f2*xcp2
c
      enddo
      xpt(1) = 0.0
      xpt(n+1) = 1.0
c
      return
      end ! cspacer

