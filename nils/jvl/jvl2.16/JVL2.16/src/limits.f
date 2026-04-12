C***********************************************************************
C    Module:  limits.f
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

      subroutine glims(xyzmin,xyzmax, lproj)
c--------------------------------------------------------------
c     Finds geometry limits
c     If lproj=t, does projection before finding limits
c--------------------------------------------------------------
      include 'jvl.inc'
      include 'jvlplt.inc'
      real xyzmin(3), xyzmax(3)
      logical lproj
c
      real pts(3,4)
c
      xyzmin(1) =  1.0e20
      xyzmin(2) =  1.0e20
      xyzmin(3) =  1.0e20
      xyzmax(1) = -1.0e20
      xyzmax(2) = -1.0e20
      xyzmax(3) = -1.0e20
c
      do n = 1, nsurf
       if(lpltsurf(n)) then
c
        j1 = jfrst(n)
        jn = jfrst(n) + nj(n)-1
        jdel = 1  !  max( (jn-j1)/5 , 1)
c
        do j = j1, jn, jdel
          pts(1,1) = rle1(1,j)
          pts(2,1) = rle1(2,j)
          pts(3,1) = rle1(3,j)
          pts(1,2) = rle1(1,j) + chord1(j)
          pts(2,2) = rle1(2,j)
          pts(3,2) = rle1(3,j)

          pts(1,3) = rle2(1,j)
          pts(2,3) = rle2(2,j)
          pts(3,3) = rle2(3,j)
          pts(1,4) = rle2(1,j) + chord1(j)
          pts(2,4) = rle2(2,j)
          pts(3,4) = rle2(3,j)

          if(lproj) call viewproj(pts,4,pts)

          do k = 1, 3
            xyzmin(k) = min( xyzmin(k), 
     &                       pts(k,1), pts(k,2), pts(k,3), pts(k,4) )
            xyzmax(k) = max( xyzmax(k),
     &                       pts(k,1), pts(k,2), pts(k,3), pts(k,4) )
          enddo
        enddo
       endif
      end do
c
c
      do n = 1, nbody
       if(lpltbody(n)) then
c
        l1 = lfrst(n)
        ln = llast(n)
c
        pts(1,1) = rl(1,l1)
        pts(2,1) = rl(2,l1)
        pts(3,1) = rl(3,l1)
        pts(1,2) = rl(1,ln)
        pts(2,2) = rl(2,ln)
        pts(3,2) = rl(3,ln)
        if(lproj) call viewproj(pts,2,pts)
c
        do k = 1, 3
          xyzmin(k) = min( xyzmin(k), pts(k,1), pts(k,2) )
          xyzmax(k) = max( xyzmax(k), pts(k,1), pts(k,2) )
        enddo
c
       endif
      end do
c
      return
      end ! glims



      subroutine grlims(xyzmin,xyzmax, lproj, tt,xyzr,dxyz)
c--------------------------------------------------------------
c     Finds rotated,translated geometry limits
c     If lproj=t, does plot projection before finding limits
c--------------------------------------------------------------
      include 'jvl.inc'
      include 'jvlplt.inc'
      real xyzmin(3), xyzmax(3)
      logical lproj
      real tt(3,3), xyzr(3), dxyz(3)
c
c
      real pts(3,2)
c
      xyzmin(1) =  1.0e20
      xyzmin(2) =  1.0e20
      xyzmin(3) =  1.0e20
      xyzmax(1) = -1.0e20
      xyzmax(2) = -1.0e20
      xyzmax(3) = -1.0e20
c
      do n = 1, nsurf
       if(lpltsurf(n)) then
c
        j1 = jfrst(n)
        jn = jfrst(n) + nj(n)-1
c
        pts(1,1) = rle1(1,j1)
        pts(2,1) = rle1(2,j1)
        pts(3,1) = rle1(3,j1)
        pts(1,2) = rle1(1,j1) + chord1(j1)
        pts(2,2) = rle1(2,j1)
        pts(3,2) = rle1(3,j1)
        call tetran(pts(1,1),tt,xyzr,dxyz)
        call tetran(pts(1,2),tt,xyzr,dxyz)
        if(lproj) call viewproj(pts,2,pts)
c
        do k = 1, 3
          xyzmin(k) = min( xyzmin(k), pts(k,1), pts(k,2) )
          xyzmax(k) = max( xyzmax(k), pts(k,1), pts(k,2) )
        enddo
c
        pts(1,1) = rle2(1,jn)
        pts(2,1) = rle2(2,jn)
        pts(3,1) = rle2(3,jn)
        pts(1,2) = rle2(1,jn) + chord2(jn)
        pts(2,2) = rle2(2,jn)
        pts(3,2) = rle2(3,jn)
        call tetran(pts(1,1),tt,xyzr,dxyz)
        call tetran(pts(1,2),tt,xyzr,dxyz)
        if(lproj) call viewproj(pts,2,pts)
c
        do k = 1, 3
          xyzmin(k) = min( xyzmin(k), pts(k,1), pts(k,2) )
          xyzmax(k) = max( xyzmax(k), pts(k,1), pts(k,2) )
        enddo
c
       endif
      end do
c
c
      do n = 1, nbody
       if(lpltbody(n)) then
c
        l1 = lfrst(n)
        ln = llast(n)
c
        pts(1,1) = rl(1,l1)
        pts(2,1) = rl(2,l1)
        pts(3,1) = rl(3,l1)
        pts(1,2) = rl(1,ln)
        pts(2,2) = rl(2,ln)
        pts(3,2) = rl(3,ln)
        call tetran(pts(1,1),tt,xyzr,dxyz)
        call tetran(pts(1,2),tt,xyzr,dxyz)
        if(lproj) call viewproj(pts,2,pts)
c
        do k = 1, 3
          xyzmin(k) = min( xyzmin(k), pts(k,1), pts(k,2) )
          xyzmax(k) = max( xyzmax(k), pts(k,1), pts(k,2) )
        enddo
c
       endif
      end do
c
      return
      end ! grlims



      subroutine axlims
c---------------------------------------------------------
c     Sets min,max,delta annotations for geometry axes
c---------------------------------------------------------
      include 'jvl.inc'
      include 'jvlplt.inc'
c
      do k = 1, 3
        axmin(k) = min( gmin(k) , 0.0 )
        axmax(k) = max( gmax(k) , 0.0 )
      enddo
c
      axdmax = max( axmax(1)-axmin(1),
     &              axmax(2)-axmin(2),
     &              axmax(3)-axmin(3) )
      do k = 1, 3
        if(axmax(k)-axmin(k) .lt. 0.25*axdmax) then
         axmin(k) = min( axmin(k) , -0.125*axdmax )
         axmax(k) = max( axmax(k) ,  0.125*axdmax )
        endif
        call axisadj(axmin(k),axmax(k),axspan(k),axdel(k),naxann(k))
      enddo
c
      return
      end ! axlims



      subroutine hidinit(lreset)
c-------------------------------------------------------------
c     Sets up surface triangle data for hidden line routine
c
c     If lreset=t resets ntri=0 to start accumulating 
c     a new set of triangles
c-------------------------------------------------------------
      include 'jvl.inc'
      include 'jvlplt.inc'
      logical lreset
c
c
      real pts(3,4)
c
      if(lreset) ntri = 0
c
      do n = 1, nsurf
        if(lpltsurf(n)) then
         do j = jfrst(n), jlast(n)
           id = j
           pts(1,1) = rle1(1,j)
           pts(2,1) = rle1(2,j)
           pts(3,1) = rle1(3,j)
           pts(1,2) = rle1(1,j) + chord1(j)
           pts(2,2) = rle1(2,j)
           pts(3,2) = rle1(3,j)
c
           pts(1,3) = rle2(1,j)
           pts(2,3) = rle2(2,j)
           pts(3,3) = rle2(3,j)
           pts(1,4) = rle2(1,j) + chord2(j)
           pts(2,4) = rle2(2,j)
           pts(3,4) = rle2(3,j)
c
           call viewproj(pts,4,pts)
           call triinit(id,1,1,pts, ntri,tri)
         end do
        endif
      end do
c
c
c      do n = 1, nbody
c       if(lpltbody(n)) then
cc
c        l1 = lfrst(n)
c        ln = llast(n)
cc
c        pts(1,1) = rl(1,l1)
c        pts(2,1) = rl(2,l1)
c        pts(3,1) = rl(3,l1)
c        pts(1,2) = rl(1,ln)
c        pts(2,2) = rl(2,ln)
c        pts(3,2) = rl(3,ln)
cc
c        call viewproj(pts,2,projpts)
cc
c       endif
c      end do
c
      return
      end ! hidinit


      subroutine hidinite(lreset, ang,pos,xyzr)
c-------------------------------------------------------------
c     Sets up surface triangle data for hidden line routine
c
c     If lreset=t resets ntri=0 to start accumulating 
c     a new set of triangles
c-------------------------------------------------------------
      include 'jvl.inc'
      include 'jvlplt.inc'
      logical lreset
      real ang(3), pos(3), xyzr(3)
c
      real pts(3,4)
      real tt(3,3), tt_ang(3,3,3)
c
      if(lreset) ntri = 0
c
      call rotens3(ang,tt,tt_ang)
c
      do n = 1, nsurf
       if(lpltsurf(n)) then
c
        j1 = jfrst(n)
        jn = jfrst(n) + nj(n)-1
        jinc = 1
c
        do j = j1, jn, jinc
          id = j
          pts(1,1) = rle1(1,j)
          pts(2,1) = rle1(2,j)
          pts(3,1) = rle1(3,j)
          pts(1,2) = rle1(1,j) + chord1(j)
          pts(2,2) = rle1(2,j)
          pts(3,2) = rle1(3,j)
c
          pts(1,3) = rle2(1,j)
          pts(2,3) = rle2(2,j)
          pts(3,3) = rle2(3,j)
          pts(1,4) = rle2(1,j) + chord2(j)
          pts(2,4) = rle2(2,j)
          pts(3,4) = rle2(3,j)
c
          call tetran(pts(1,1),tt,xyzr,pos)
          call tetran(pts(1,2),tt,xyzr,pos)
          call tetran(pts(1,3),tt,xyzr,pos)
          call tetran(pts(1,4),tt,xyzr,pos)
c
          call viewproj(pts,4,pts)
          call triinit(id,1,1,pts, ntri,tri)
c
        end do
c
       endif
      end do
c
c
c      do n = 1, nbody
c       if(lpltbody(n)) then
cc
c        l1 = lfrst(n)
c        ln = llast(n)
cc
c        pts(1,1) = rl(1,l1)
c        pts(2,1) = rl(2,l1)
c        pts(3,1) = rl(3,l1)
c        pts(1,2) = rl(1,ln)
c        pts(2,2) = rl(2,ln)
c        pts(3,2) = rl(3,ln)
cc
c        call viewproj(pts,2,projpts)
cc
c       endif
c      end do
c
      return
      end ! hidinit
