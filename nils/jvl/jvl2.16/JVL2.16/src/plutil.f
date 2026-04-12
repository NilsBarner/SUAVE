C***********************************************************************
C    Module:  plutil.f
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
 
      subroutine xyplot(n,x,y,xoff,xsf,yoff,ysf,ilin,sh,isym)
      real x(n), y(n)

      if(isym.le.0) call xyline(n,x,y,xoff,xsf,yoff,ysf,ilin)
      if(isym.ne.0) call xysymb(n,x,y,xoff,xsf,yoff,ysf,sh,iabs(isym))

      return
      end ! xyplot


      subroutine xaxis2(x1,y1,xaxt,dxann,fann,dann,iflag,cst,ndig)
c---- xaxis2 differs from libplt xaxis by having a flag to suppress 
c     both end annotations rather than the zero annotation.
c------------------------------------------------------------------
c     x1,y1  starting point of x axis
c     xaxt   length of x axis 
c     dxann  distance between annotations
c     fann   first annotation value
c     dann   delta annotation value
c     iflag  flag to suppress end annotations 
c            = 0    all annotations
c            = 1    suppress first annotation
c            = 2    suppress last  annotation
c            = 3    suppress first and last annotations
c     cst    character size   ( - = annotation above axis)
c     ndig   number of digits to right of decimal point
c            = -1   no decimal point
c            = -2   number of digits determined internally
c------------------------------------------------------------------
      xax = abs(xaxt)
      if(xax.le.0) return
      cs = abs(cst)

c---- determine # of digits to use for annotations
      if(ndig.le.-2) then
        nd = 1 - max( 0 , int(log10(dann)) )
        if(dann*10**nd - aint(dann*10**nd+0.01) .gt. 0.01) nd = nd + 1
        if(dann*10**nd - aint(dann*10**nd+0.01) .gt. 0.01) nd = nd + 1
      else
        nd = ndig
      endif

c---- x-axis
      call plot(x1,y1,3)
      call plot(x1+xax,y1,2)
      nann = 1 + ifix(xax/dxann + 0.1)

c---- annotate x-axis
      do 10 nt=1, nann
        xt = x1 + dxann*float(nt-1)
c---- skip annotation for first or last annotation position as given by iflg
        if(mod(iflag,2).eq.1 .and. nt.eq.1)    go to 10
        if(iflag.gt.1        .and. nt.eq.nann) go to 10

        call plot(xt,y1-0.2*cs,3)
        call plot(xt,y1+0.2*cs,2)
        rn = fann + dann*float(nt-1)
        grn = 0.
        if(rn.ne.0.0) grn = alog10(abs(rn)+0.5/10.0**nd)
        grn = max(grn,0.0)
        nabc = int(grn) + 2 + nd
        width = 0.95*cs*float(nabc)
        if(rn.lt.0.0) width = width + cs
        xnum = xt - 0.5*width
        ynum = y1 - 2.1*cs
        if(cst.lt.0.0) ynum = y1 + 0.9*cs

        call plnumb(xnum,ynum,cs,rn,0.0,nd)
   10 continue

      return
      end ! xaxis2
 

      subroutine viewinit(azim, elev, tilt, robinv)
c----------------------------------------------------------------------------
c     Sets up projection for points in 3-D cartesian space
c     onto a 2-D plane from the viewpoint of an observer
c     at a specified location.  this can be used to "view" 
c     a 3-D object described by a set of points on a planar
c     2-D graphics screen. 
c        The viewing  plane, which has its own x,y coordinate
c     system, always faces the observer but can be turned
c     around the viewing axis, thus simulating the observer
c     tilting his head while looking at the object.  this tilt
c     is specified by giving a vector which "points up" relative
c     to the observer.
c        The distance of the observer from the object is specified
c     explicitly.  this does not affect much the size of the viewed
c     object, since the viewing plane contains the 3-D space origin
c     and hence is at or near the object.  it does however affect the
c     apparent distortion of the object due to perspective.  this
c     is very useful to convey the 3-dimensionality of the object.
c     if the observer is very very far away, there is no distortion
c     (as in a mechanical drawing).
c
c     azim      azimuth   angle of observer in degrees   (input)
c     elev      elevation angle of observer in degrees   (input)
c     tilt      tilt      angle of observer in degress   (input)
c     robinv    1/(distance to observer)                 (input)
c
c     xyzob(.)  cartesian vector pointing towards observer  (internal output)
c               (magnitude irrelevant)
c     xyzup(.)  cartesian vector which points "up" from the 
c               observer's viewpoint (magnitude irrelevant) (internal output)
c----------------------------------------------------------------------------
      common /viewdata/ rinv,
     &       xihat, yihat, zihat,
     &       xjhat, yjhat, zjhat,
     &       xkhat, ykhat, zkhat

      real xyzob(3), xyzup(3)

      data dtr / 0.01745329 /

      xyzob(1) = -cos(dtr*azim) * cos(dtr*elev)
      xyzob(2) =  sin(dtr*azim) * cos(dtr*elev)
      xyzob(3) =                  sin(dtr*elev)

      xyzup(1) =  0.
      xyzup(2) = -sin(dtr*tilt)
      xyzup(3) =  cos(dtr*tilt)

c---- unit view vector out of viewing plane (towards observer)
      xkhat = xyzob(1)/sqrt(xyzob(1)**2 + xyzob(2)**2 + xyzob(3)**2)
      ykhat = xyzob(2)/sqrt(xyzob(1)**2 + xyzob(2)**2 + xyzob(3)**2)
      zkhat = xyzob(3)/sqrt(xyzob(1)**2 + xyzob(2)**2 + xyzob(3)**2)

c---- vector along plane's local x coordinate: (up vector)x(view vector)
      xip = xyzup(2)*zkhat - xyzup(3)*ykhat
      yip = xyzup(3)*xkhat - xyzup(1)*zkhat
      zip = xyzup(1)*ykhat - xyzup(2)*xkhat

c---- normalize plane's x coordinate vector
      xihat = xip/sqrt(xip**2 + yip**2 + zip**2)
      yihat = yip/sqrt(xip**2 + yip**2 + zip**2)
      zihat = zip/sqrt(xip**2 + yip**2 + zip**2)

c---- unit vector along plane's y coord.: (view vector)x(x unit vector)
      xjhat = ykhat*zihat - zkhat*yihat
      yjhat = zkhat*xihat - xkhat*zihat
      zjhat = xkhat*yihat - ykhat*xihat

      rinv = robinv

      return
      end ! viewinit


      subroutine viewproj(xyz,n,xyzproj)
      real xyz(3,n), xyzproj(3,n)
c-------------------------------------------------------------------------
c     Projects one or more points in 3-D cartesian space
c     onto a 2-D plane from the viewpoint of an observer
c     at a specified location using the projection transformation
c     setup in viewinit.
c
c     This can be used to "view" a 3-D object described by a set 
c     of points on a planar 2-D graphics screen.  this routine also 
c     returns the out-of-plane vector for depth comparison.
c        The viewing  plane, which has its own x,y coordinate
c     system, always faces the observer but can be turned
c     around the viewing axis, thus simulating the observer
c     tilting his head while looking at the object.  this tilt
c     is specified by giving a vector which "points up" relative
c     to the observer.
c        The distance of the observer from the object is specified
c     explicitly.  this does not affect much the size of the viewed
c     object, since the viewing plane contains the 3-D space origin
c     and hence is at or near the object.  it does however affect the
c     apparent distortion of the object due to perspective.  this
c     is very useful to convey the 3-dimensionality of the object.
c     if the observer is very very far away, there is no distortion
c     (as in a mechanical drawing).
c
c     xyz(.)       cartesian point coordinates                  (input)
c     n            number of points                             (input)
c     xyzproj(.)   projected point coordinates on viewing plane (output)
c-------------------------------------------------------------------------
      common /viewdata/ rinv,
     &       xihat, yihat, zihat,
     &       xjhat, yjhat, zjhat,
     &       xkhat, ykhat, zkhat

c---- go over all points
      do i=1, n
        rdoti = xyz(1,i)*xihat + xyz(2,i)*yihat + xyz(3,i)*zihat
        rdotj = xyz(1,i)*xjhat + xyz(2,i)*yjhat + xyz(3,i)*zjhat
        rdotk = xyz(1,i)*xkhat + xyz(2,i)*ykhat + xyz(3,i)*zkhat

c------ viewing-axis component of vector
        rkx = rdotk*xkhat
        rky = rdotk*ykhat
        rkz = rdotk*zkhat

c------ projected vector scaling factor due to perspective
        vscal = 1.0 / sqrt( (xkhat-rinv*rkx)**2
     &                    + (ykhat-rinv*rky)**2
     &                    + (zkhat-rinv*rkz)**2 )

c------ dot vector into plane coordinate system unit vectors, and scale
        xyzproj(1,i) = vscal * rdoti
        xyzproj(2,i) = vscal * rdotj
        xyzproj(3,i) = vscal * rdotk
      end do

      return
      end ! viewproj


      subroutine oplset(idev,idevh,ipslu,
     &                  plsize,plotar,
     &                  xmarg,ymarg,xpage,ypage,
     &                  csize,scrnfr,lcurs,lcrev)
      logical lcurs,lcrev
c-----------------------------------------------------------
c     Allows user modification of various plot parameters.
c-----------------------------------------------------------
      character*1 var
      character*4 comand
      character*128 comarg
      character*10 chcurs, chland
      dimension iinput(20)
      dimension rinput(20)
      logical error, lgraph, lcolor, lfiles

 1000 format(a)

 1    continue
      if(lcurs) then
       chcurs = 'cursor    '
      else
       chcurs = 'keyboard  '
      endif

      if(scrnfr.gt.0.0) then
       chland = 'Landscape '
      else
       chland = 'Portrait  '
      endif

      asf = abs(scrnfr)

      lgraph = idev  .ge. 1
      lcolor = idevh .eq. 4
      lfiles = ipslu .lt. 0

      write(*,2000) lgraph, lcolor, lfiles,
     &              plotar, plsize,
     &              xpage,ypage, xmarg,ymarg, 
     &              csize, asf,
     &              chland, chcurs
 2000 format(' ...............................................'
     &     //'  G raphics-enable flag        ', l2,
     &      /'  C olor postscript output?    ', l2,
     &      /'  I ndividual ps file output?  ', l2,
     &      /'  A spect ratio of plot object ', f8.4
     &      /'  S ize of plot object         ', f6.2,'"'
     &      /'  P age dimensions             ', f6.2,' x',f6.2,'"'
     &      /'  M argins from page edges     ', f6.2,'",',f6.2,'"'
     &      /'  F ont size (relative)        ', f8.4
     &      /'  W indow/screen size fraction ', f8.4
     &      /'  O rientation of plot:        ', a 
     &      /'  B lowup input method:        ', a )
ccc     &      /'  R everse-video output?       ', l2 )
c
c   A B C D E F G H I J K L M N O P Q R S T U V W X Y Z
c   x x x     x x           x   x x     x       x     

 5    call askc('      Option, Value   (or <return>) ^',comand,comarg)

      do i=1, 20
        iinput(i) = 0.0
        rinput(i) = 0.0
      enddo
      ninput = 0
      call getint(comarg,iinput,ninput,error)
      ninput = 0
      call getflt(comarg,rinput,ninput,error)

      var = comand(1:1)
      if (var.eq.'0' .or. var.eq.' ') then
       return

      elseif (index('Gg',var).ne.0) then
        if(idev.eq.0) then
         idev = 1
        else
         idev = 0
        endif

      elseif (index('Cc',var).ne.0) then
        if(idevh.eq.2) then
         idevh = 4
        else
         idevh = 2
        endif

      elseif (index('Ii',var).ne.0) then
        if(ipslu.eq.0) then
         ipslu = -1
        else
         ipslu = 0
        endif

      elseif (index('Ss',var).ne.0) then
        if(ninput.ge.1) then
          plsize = rinput(1)
        else
          call askr('Enter size (in)^',plsize)
        endif

      elseif (index('Aa',var).ne.0) then
        if(ninput.ge.1) then
          plotar = rinput(1)
        else
          call askr('Enter aspect ratio^',plotar)
        endif

      elseif (index('Pp',var).ne.0) then
        if(ninput.ge.2) then
          xpage = rinput(1)
          ypage = rinput(2)
        elseif(ninput.ge.1) then
          xpage = rinput(1)
          call askr('Enter page y dimension (in)^',ypage)
        else
          call askr('Enter page x dimension (in)^',xpage)
          call askr('Enter page y dimension (in)^',ypage)
        endif

      elseif (index('Mm',var).ne.0) then
        if(ninput.ge.2) then
          xmarg = rinput(1)
          ymarg = rinput(2)
        elseif(ninput.ge.1) then
          xmarg = rinput(1)
          call askr('Enter page y margin (in)^',ymarg)
        else
          call askr('Enter page x margin (in)^',xmarg)
          call askr('Enter page y margin (in)^',ymarg)
        endif

      elseif (index('Ff',var).ne.0) then
        if(ninput.ge.1) then
          csize = rinput(1)
        else
          call askr('Enter character font size / plot size^',csize)
        endif

      elseif (index('Ww',var).ne.0) then
        if(ninput.ge.1) then
          asf = rinput(1)
        else
          call askr('Enter window/screen size fraction^',asf)
        endif
        scrnfr = sign( scrnfr , asf )

      elseif (index('Bb',var).ne.0) then
        lcurs = .not. lcurs

      elseif (index('Oo',var).ne.0) then
        scrnfr = -scrnfr
        write(*,*)
        write(*,*) 'Swapping x,y page dimensions'
        xtmp = xpage
        ytmp = ypage
        xpage = ytmp
        ypage = xtmp

      elseif (index('Rr',var).ne.0) then
        lcrev = .not. lcrev

      else
        write(*,*) '*** item not recognized ***'
      endif
      go to 1

      end ! oplset





      subroutine axisadj2(xmin,xmax,xspan,deltax,ntics)
c...make scaled axes with engineering increments between tics
c
c   input:    xmin, xmax - input range for which scaled axis is desired
c
c   output:   xmin, xmax - adjusted range for scaled axis
c             xspan      - adjusted span of scaled axis
c             deltax     - increment to be used for scaled axis
c             ntics      - number of axis tics (each deltax long)
c
      real    xmin,xmax,xspan,deltax,xinc,xinctbl(4)
      integer ntics,i
      data    xinctbl / 0.1, 0.2, 0.5, 1.0 /
c
      xspan1 = xmax-xmin
      if (xspan1.eq.0.) xspan1 = 1.
c
      xpon = ifix(log10(xspan1))
      xspan = xspan1 / 10.**xpon
c
      do i = 1, 4
        xinc = xinctbl(i)
        ntics = 1 + ifix(xspan/xinc)
        if (ntics.le.12) go to 1
      end do
c
   1  deltax = xinc*10.**xpon
      xmin = deltax*  ifloor(xmin/deltax)
      xmax = deltax*iceiling(xmax/deltax)
      xspan = xmax - xmin
      return
      end
