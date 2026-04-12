C***********************************************************************
C    Module:  pltops.f
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

      subroutine pltini(idev1)
      include 'jvlplt.inc'

c---- terminate old plot if any
      if(lplot) call plend
c
c---- initialize new plot
      call plopen(scrnfrac,ipslu,idev1)
      lplot = .true.

      call drawtobuffer

c---- set x-window size in inches (might have been resized by user)
      call getwinsize(xwind,ywind)
c
      if(lcrev) then
       call bgfill
      endif
c
c---- draw plot page outline offset by margins
      call newpen(5)
      if(xmarg .gt. 0.0) then
        call plotabs(      xmarg,      ymarg,3)
        call plotabs(      xmarg,ypage-ymarg,2)
        call plotabs(xpage-xmarg,      ymarg,3)
        call plotabs(xpage-xmarg,ypage-ymarg,2)
      endif
      if(ymarg .gt. 0.0) then
        call plotabs(      xmarg,      ymarg,3)
        call plotabs(xpage-xmarg,      ymarg,2)
        call plotabs(      xmarg,ypage-ymarg,3)
        call plotabs(xpage-xmarg,ypage-ymarg,2)
      endif
      call newpen(1)
c
      call plotabs(xmarg,ymarg,-3)
      call newclipabs( xmarg, xpage-xmarg, ymarg, ypage-ymarg )

      return
      end ! pltini



      subroutine pltseg(vec,alf,n)
c...purpose:    plot portions of vectors in vec given by normalized
c               segment start and end points in alf
c
      include 'jvlplt.inc'
      real vec(3,2),alf(2,n)
c
      xmod(xtmp) = sf * (xtmp - xoff)
      ymod(ytmp) = sf * (ytmp - yoff)
c
      x1 = vec(1,1)
      y1 = vec(2,1)
c
      dx = vec(1,2) - vec(1,1)
      dy = vec(2,2) - vec(2,1)
c
      i = 1
      xa = x1 + alf(1,i)*dx
      ya = y1 + alf(1,i)*dy
      call plot(xmod(xa),ymod(ya),3)
c
      xa = x1 + alf(2,i)*dx
      ya = y1 + alf(2,i)*dy
      call plot(xmod(xa),ymod(ya),2)
c
      do 10 i=2, n
c
        if(abs(alf(1,i)-alf(2,i-1)) .gt. 1.0e-3) then
         xa = x1 + alf(1,i)*dx
         ya = y1 + alf(1,i)*dy
         call plot(xmod(xa),ymod(ya),3)
        endif
c
        xa = x1 + alf(2,i)*dx
        ya = y1 + alf(2,i)*dy
        call plot(xmod(xa),ymod(ya),2)
c
   10 continue
c
      return
      end ! pltseg



      subroutine pltpoly(pts,npts)
c...purpose:   plot polygon given by npts vertices in pts array of xyz points
c
      include 'jvlplt.inc'
      real pts(3,*)
c
      xmod(xtmp) = sf * (xtmp - xoff)
      ymod(ytmp) = sf * (ytmp - yoff)
c
      call plot(xmod(pts(1,1)),ymod(pts(2,1)),3)
      do k=2, npts
        call plot(xmod(pts(1,k)),ymod(pts(2,k)),2)
      end do
      call plot(xmod(pts(1,1)),ymod(pts(2,1)),2)
      return
      end ! pltpoly


      subroutine pltint(pt,n,cs,loffset)
c...purpose:   plot integer n at location pt 
c              with character size cs
c
      include 'jvlplt.inc'
      real pt(3)
      logical loffset
c
      xmod(xtmp) = sf * (xtmp - xoff)
      ymod(ytmp) = sf * (ytmp - yoff)
c
      fltn = float(n)
      x = xmod(pt(1))
      y = ymod(pt(2))
      if(loffset) then
       x = x + 0.80*cs
       y = y - 0.50*cs
      else
       nchar = 1 + int(log10(abs(fltn))+0.01)
       x = x - 0.5*cs*float(nchar)
       y = y - 0.5*cs
       if(fltn.lt.0.0) x = x - cs
      endif
      call plnumb(x,y,cs,fltn,0.,-1)
c
      return
      end ! pltint


      subroutine pltflt(pt,r,cs,loffset,ndig)
c...purpose:   plot float r at location pt 
c              with character size cs
c
      include 'jvlplt.inc'
      real pt(3)
      logical loffset
c
      xmod(xtmp) = sf * (xtmp - xoff)
      ymod(ytmp) = sf * (ytmp - yoff)
c
      absr = abs(r)
c
c---- determine # of digits to use for diplay
      if(ndig.le.-2) then
       nd = 1 - max( 0 , int(log10(absr)) )
       if(r*10**nd - aint(r*10**nd+0.01) .gt. 0.01) nd = nd + 1
       if(r*10**nd - aint(r*10**nd+0.01) .gt. 0.01) nd = nd + 1
      else
       nd = ndig
      endif
c
      x = xmod(pt(1))
      y = ymod(pt(2))
      if(loffset) then
       x = x + 0.80*cs
       y = y - 0.50*cs
      else
       nchar = 1 + int(log10(absr)) + 1 + nd
       x = x - 0.5*cs*float(nchar)
       y = y - 0.5*cs
       if(r.lt.0.0) x = x - cs
      endif
      call plnumb(x,y,cs,r,0.0,nd)
c
      return
      end ! pltflt


      subroutine pltarrow(pt1,pt2)
c...purpose:   plot an arrow vector from pt1 to pt2
c
      include 'jvlplt.inc'
      real pt1(3), pt2(3)
c
      xmod(xtmp) = sf * (xtmp - xoff)
      ymod(ytmp) = sf * (ytmp - yoff)
c
      xave = pt1(1)
      yave = pt1(2)
      xhed = pt2(1)
      yhed = pt2(2)
c
      dx = xhed - xave
      dy = yhed - yave
c
      call plot(xmod(xave),ymod(yave),3)
      call plot(xmod(xhed),ymod(yhed),2)
c
      x1 = xave + 0.8*dx + 0.02*dy
      y1 = yave + 0.8*dy - 0.02*dx
      x2 = xave + 0.8*dx - 0.02*dy
      y2 = yave + 0.8*dy + 0.02*dx
      call plot(xmod(x1  ),ymod(y1  ),2)
      call plot(xmod(x2  ),ymod(y2  ),2)
      call plot(xmod(xhed),ymod(yhed),2)
c
      return
      end ! pltarrow




      subroutine rotatept(pt0,pt1,dir,cosr,sinr)
c...purpose:  rotate point pt1 about vector in direction dir through
c             point pt0 by angle with cosine cosr and sine sinr.
c             rotated point is returned in pt1.
c
c
      real pt0(3), pt1(3), dir(3), dpt(3), ep(3), eq(3)
c
      do l = 1, 3
        dpt(l) = pt1(l) - pt0(l)
      end do
c-----ep = normal-vector component perpendicular to hinge line
      endot = dot(dpt,dir)
      ep(1) = dpt(1) - endot*dir(1)
      ep(2) = dpt(2) - endot*dir(2)
      ep(3) = dpt(3) - endot*dir(3)
c-----eq = unit vector perpendicular to both ep and dir
      call cross(dir,ep,eq)
c-----rotated vector consists of sin,cos parts from ep and eq,
c-    with hinge-parallel component endot restored
      do l = 1, 3
        dpt(l) = ep(l)*cosr + eq(l)*sinr + endot*dir(l)
        pt1(l) = pt0(l) + dpt(l)
      end do
c
      return
      end ! rotatept


      subroutine offini       
c---- set initial scaling and offset parameters   
c
      include 'jvlplt.inc'
c
      xrange = xmax-xmin
      yrange = ymax-ymin
      if(xrange.eq.0.) xrange = 1.0
      if(yrange.eq.0.) yrange = 1.0
      sf = 0.95*min( 1.0/xrange , plotar/yrange )
      xoff = xmin - 0.5*(1.0    - sf*xrange)/sf
      yoff = ymin - 0.5*(plotar - sf*yrange)/sf
c
      return        
      end ! offini  


      subroutine offget       
c---- get offsets for zoom from user interaction (mouse)
c
      include 'jvlplt.inc'
      character*1 ckey
c
      sh = 2.0
c
      write(*,*)
      write(*,*) 'mark off corners of blowup area'
      write(*,*) '(2 spaces default to current area)'       
      call getcursorxy(xx1,yy1,ckey)
      call plsymb(xx1,yy1,sh,3,0.0,0)
      call plflush
      call getcursorxy(xx2,yy2,ckey)
      call plsymb(xx2,yy2,sh,3,0.0,0)
      call plflush
      if(abs(xx1-xx2).lt.0.01 .and. abs(yy1-yy2).lt.0.01) return      
c
      xoff = min(xx1,xx2)/sf + xoff
      yoff = min(yy1,yy2)/sf + yoff
      xdif = abs(xx2 - xx1)/sf
      ydif = abs(yy2 - yy1)/sf   
      if(xdif.eq.0.0) xdif = 1.0e-5
      if(ydif.eq.0.0) ydif = 1.0e-5
      sf = min( 1.0/xdif , plotar/ydif )     
c
c
c---- re-center the blowup
c      xdif = max(xdif,ydif/0.75)
      ydif = max(0.75*xdif,ydif)
      xoff = xoff - 0.5*(1.0    - sf*xdif)/sf
      yoff = yoff - 0.5*(plotar - sf*ydif)/sf
c
      return        
      end ! offget  


      subroutine bgfill
      include 'jvlplt.inc'
      real xbox(5), ybox(5)
      data xbox / 0.0 , 11.0 , 11.0 , 0.0 , 0.0 /
      data ybox / 0.0 ,  0.0 ,  8.5 , 8.5 , 0.0 /
c
      call newcolorname('black')
      if(scrnfrac .gt. 0.0) then
       call polyline(xbox,ybox,5,1)
      else
       call polyline(ybox,xbox,5,1)
      endif
      call newcolorname('white')
c
      return
      end ! bgfill
