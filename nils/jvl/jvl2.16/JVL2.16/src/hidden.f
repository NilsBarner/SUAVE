C***********************************************************************
C    Module:  hidden.f
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

      subroutine triinit(id,nrows,ncols,pts, ntri,tri)
c--------------------------------------------------------------------------
c     Sets up triangle data for hidden line routine from 
c     quadrilateral input polygons.
c
c          this routine takes a grid of quadrilaterals and defines
c          two triangles per quad and puts them into the triangle 
c          arrays.  the input grid of quads is in the form of nrowsxncols
c          of quads, defined by nrows+1 x ncols+1 points in the pts array. 
c          the id index corresponds to the first quad in the lattice and 
c          is incremented as triinit steps through the quad array to define
c          triangles.  thus the output value of id is id + nrowsxncols - 1.
c          this allows the user to use the id option in the hidden line 
c          routines to avoid visibility tests for input features that  
c          lie on triangles with the same id index.
c--------------------------------------------------------------------------
      real pts(3,*)
      real xyzlim(3,2), tri(16,*)
c
      it = ntri
      idx = id - 1
c
c--- for a rectangular grid of quadrilaterals defined by nrows, ncols
c    polygons, find the vertices of each quadrilateral and divide it into
c    two triangles.
c
      do 110 j = 1, ncols
        do 1110 k = 1, nrows
c          
c--- increment id# to get id for next element
          idx  = idx + 1
c
          ip1 = (nrows+1)*(j-1) + k
          ip2 = ip1 + 1
          ip3 = ip1 + nrows+1
          ip4 = ip3 + 1
c
c          1  = x1
c          2  = y1
c          3  = z1
c          4  = x2
c          5  = y2
c          6  = z2
c          7  = x3
c          8  = y3
c          9  = z3
c         10  = min(x)
c         11  = min(y)
c         12  = min(z)
c         13  = max(x)
c         14  = max(y)
c         15  = max(z)
c         16  = id reference number
c
c--- first triangle
          it = it+1
c 
          tri(1,it) = pts(1,ip1)
          tri(2,it) = pts(2,ip1)
          tri(3,it) = pts(3,ip1)
c         
          tri(4,it) = pts(1,ip2)
          tri(5,it) = pts(2,ip2)
          tri(6,it) = pts(3,ip2)
c         
          tri(7,it) = pts(1,ip3)
          tri(8,it) = pts(2,ip3)
          tri(9,it) = pts(3,ip3)
c         
          tri(10,it) = min( tri(1,it) , tri(4,it) , tri(7,it) )
          tri(11,it) = min( tri(2,it) , tri(5,it) , tri(8,it) )
          tri(12,it) = min( tri(3,it) , tri(6,it) , tri(9,it) )
          tri(13,it) = max( tri(1,it) , tri(4,it) , tri(7,it) )
          tri(14,it) = max( tri(2,it) , tri(5,it) , tri(8,it) )
          tri(15,it) = max( tri(3,it) , tri(6,it) , tri(9,it) )
c
          tri(16,it) = float(idx)
c
c--- second triangle          
          it = it+1
c         
          tri(1,it) = pts(1,ip3)
          tri(2,it) = pts(2,ip3)
          tri(3,it) = pts(3,ip3)
c         
          tri(4,it) = pts(1,ip2)
          tri(5,it) = pts(2,ip2)
          tri(6,it) = pts(3,ip2)
c         
          tri(7,it) = pts(1,ip4)
          tri(8,it) = pts(2,ip4)
          tri(9,it) = pts(3,ip4)
c         
          tri(10,it) = min( tri(1,it) , tri(4,it) , tri(7,it) )
          tri(11,it) = min( tri(2,it) , tri(5,it) , tri(8,it) )
          tri(12,it) = min( tri(3,it) , tri(6,it) , tri(9,it) )
          tri(13,it) = max( tri(1,it) , tri(4,it) , tri(7,it) )
          tri(14,it) = max( tri(2,it) , tri(5,it) , tri(8,it) )
          tri(15,it) = max( tri(3,it) , tri(6,it) , tri(9,it) )
c
          tri(16,it) = float(idx)
c         
 1110   continue
  110 continue
c
      ntri = it
      id = idx 
c
      return
      end ! triinit



      subroutine hidlin(xyzlin,id,ngrp,ilstgrp,xyzgrp,
     &                   ntri,xyztri,nseg,alfseg)
c
c      this is a hidden line routine.  given a set of triangles, and a
c      single line as input, it computes the line segments (if any)
c      which are visible.
c
c      xyzlin = array of data about input line (unchanged by routine)
c          1,1  = x1
c          2,1  = y1
c          3,1  = z1
c          1,2  = x2
c          2,2  = y2
c          3,2  = z2
c      id = id index for input line 
c           may be used to designate triangle or collection of triangles 
c           that contain this line and should be skipped for visibility 
c           tests.  if id=0 no triangles will be skipped.
c
c      ngrp = number of triangle groups 
c             (if ngrp<=0 no group minmax tests will be done)
c      ilstgrp(n) = index of last triangle in group n
c      xyzgrp(.n) = x,y group limits
c             1   = min(x)
c             2   = min(y)
c             3   = min(z)
c             4   = max(x)
c             5   = max(y)
c             6   = max(z)
c
c      ntri = number of triangles (used if ngrp=0)
c      xyztri = array of triangle data
c          1  = x1
c          2  = y1
c          3  = z1
c          4  = x2
c          5  = y2
c          6  = z2
c          7  = x3
c          8  = y3
c          9  = z3
c         10  = min(x)
c         11  = min(y)
c         12  = min(z)
c         13  = max(x)
c         14  = max(y)
c         15  = max(z)
c         16  = id field, id associated with this triangle 
c
c      nseg   = (input)  dimension of alfseg = max. no. of segments
c               (output) number of visible line segments
c      alfseg = output values describing the visible line segments as
c               fractions of the original line
c          1  = beginning of segment
c          2  =       end of segment
c
       real xyzlin(3,2), xyztri(16,*), alfseg(2,nseg)
       real xyzgrp(6,ngrp)
       integer ilstgrp(ngrp)
c
       data eps/1.0e-7/
c
       nmxseg = nseg
       nseg = 0
c
c----- initialize a few things
       x1 = xyzlin(1,1)
       y1 = xyzlin(2,1)
       z1 = xyzlin(3,1)
       x2 = xyzlin(1,2)
       y2 = xyzlin(2,2)
       z2 = xyzlin(3,2)
c
       xmin = min(x1,x2)
       ymin = min(y1,y2)
       zmin = min(z1,z2)
       xmax = max(x1,x2)
       ymax = max(y1,y2)
       zmax = max(z1,z2)
c---check for degenerate line
       if(xmax.eq.xmin .and. ymax.eq.ymin .and. zmax.eq.zmin) return
c
       nseg = 1
       alfseg(1,1) = 0.0
       alfseg(2,1) = 1.0
c
c----- main processing loop
       klast = 0
       nngrp = max(1,ngrp)
c
       do 100 igrp=1, nngrp
c
       kt1 = klast + 1
       klast = ntri
c
       if(ngrp.ge.1) then   
c--- check for group with no triangles (last index <= 0)
       if(ilstgrp(igrp).le.0) go to 100
        klast = ilstgrp(igrp)       
c..... preliminary purge of whole group
        if( (xyzgrp(1,igrp).ge.xmax .or. 
     &       xyzgrp(4,igrp).le.xmin) .or.
     &      (xyzgrp(2,igrp).ge.ymax .or. 
     &       xyzgrp(5,igrp).le.ymin) ) go to 100
c     &      (xyzgrp(6,igrp).le.zmin) ) go to 100
       endif
c
       do 1 kt=kt1, klast
c
c....... second purge of most cells
         if( (xyztri(10,kt).ge.xmax .or. 
     &        xyztri(13,kt).le.xmin)   .or.
     &       (xyztri(11,kt).ge.ymax .or. 
     &        xyztri(14,kt).le.ymin)   .or.
     &       (xyztri(15,kt).le.zmin) ) go to 1
c
c....... if id of this line matches a triangle id skip visibility tests
         if(id.ne.0) then
          idtri = ifix(xyztri(16,kt))
          if(id.eq.idtri) go to 1
         endif
c
c---check for line as a side of masking triangle
         nsame=0
         do is=1,2
           do l=1, 3
             if( abs(xyzlin(1,is)-xyztri(1+3*(l-1),kt)).lt.eps .and.
     &           abs(xyzlin(2,is)-xyztri(2+3*(l-1),kt)).lt.eps .and.
     &           abs(xyzlin(3,is)-xyztri(3+3*(l-1),kt)).lt.eps ) then
               nsame=nsame+1
             endif
           end do
         end do
         if(nsame.ge.2) then
          go to 1
         endif
c
c--- find where cell nodes are relative to line
         n1 = 0
           if( (xyztri(1,kt)-x1)*(y2-y1) .gt. 
     &         (xyztri(2,kt)-y1)*(x2-x1) )    n1 = n1+1
           if( (xyztri(4,kt)-x1)*(y2-y1) .gt. 
     &         (xyztri(5,kt)-y1)*(x2-x1) )    n1 = n1+2
           if( (xyztri(7,kt)-x1)*(y2-y1) .gt. 
     &         (xyztri(8,kt)-y1)*(x2-x1) )    n1 = 3-n1
c
c--- skip tests if triangle points all on one side of line
         if(n1.eq.0) go to 1
c
c--- set up temporary variables in canonical orientation where the line
c    passes through sides 1-2 and 1-3 of triangle (point 1 is common point)
         n2 = mod(n1,3) + 1           
         n3 = mod(n2,3) + 1           
c
         xt1 = xyztri(3*n1-2,kt)
         yt1 = xyztri(3*n1-1,kt)
         zt1 = xyztri(3*n1  ,kt)
         xt2 = xyztri(3*n2-2,kt)
         yt2 = xyztri(3*n2-1,kt)
         zt2 = xyztri(3*n2  ,kt)
         xt3 = xyztri(3*n3-2,kt)
         yt3 = xyztri(3*n3-1,kt)
         zt3 = xyztri(3*n3  ,kt)
c
c--- compute xi, eta, zeta values for line endpoints in triangle coords
         det = ( (xt2-xt1)*(yt3-yt1) - (xt3-xt1)*(yt2-yt1) )
         if(det.eq.0.0) go to 1
c
         detinv = 1.0 / det
         xi1 = detinv * ( (x1 -xt1)*(yt3-yt1) - (xt3-xt1)*(y1 -yt1) )
         et1 = detinv * ( (xt2-xt1)*(y1 -yt1) - (x1 -xt1)*(yt2-yt1) )
         xi2 = detinv * ( (x2 -xt1)*(yt3-yt1) - (xt3-xt1)*(y2 -yt1) )
         et2 = detinv * ( (xt2-xt1)*(y2 -yt1) - (x2 -xt1)*(yt2-yt1) )
c
c--- check for line parallel to triangle sides (skip hide tests)
         if( (abs(xi1).lt.eps .and. abs(xi2).lt.eps) .or.
     &       (abs(et1).lt.eps .and. abs(et2).lt.eps) ) go to 1
c
c--- get z coordinates relative to triangle plane
         if(abs(xi1-1.0).lt.eps) xi1 = 1.0
         if(abs(et1-1.0).lt.eps) et1 = 1.0
         ze1 = (z1-zt1) - xi1*(zt2-zt1) - et1*(zt3-zt1)
c
         if(abs(xi2-1.0).lt.eps) xi2 = 1.0
         if(abs(et2-1.0).lt.eps) et2 = 1.0
         ze2 = (z2-zt1) - xi2*(zt2-zt1) - et2*(zt3-zt1)
c
c--- skip this triangle test if line is outside xi=0,eta=0 or xi+eta-1 edges
         if(max(xi1,xi2).le.0.0 .or.
     &      max(et1,et2).le.0.0 .or.
     &      min(xi1+et1,xi2+et2) .gt. 1.0) go to 1
c
c--- check for line completely above triangle plane
         if(min(ze1,ze2).ge.0.) go to 1
c
         alf1 = 0.0
         alf2 = 1.0
c
c--- intersect the line with sides 1-2 and 1-3 
c    this defines the alf (relative coordinate on line 0->1)
         denbet = (x2-x1)*(yt2-yt1) - (y2-y1)*(xt2-xt1)
         dengam = (x2-x1)*(yt3-yt1) - (y2-y1)*(xt3-xt1)
         if( abs(denbet).lt.eps .or.
     &       abs(dengam).lt.eps ) go to 1

         bet = ( (0.5*(xt1+xt2)-x1)*(yt2-yt1) -
     &           (0.5*(yt1+yt2)-y1)*(xt2-xt1) ) / denbet
c
         gam = ( (0.5*(xt1+xt3)-x1)*(yt3-yt1) -
     &           (0.5*(yt1+yt3)-y1)*(xt3-xt1) ) / dengam
c
         alf1 = max(alf1,min(bet,gam))
         alf2 = min(alf2,max(bet,gam))
c
c--- if line endpoint z values straddle the triangle plane, also 
c    check the line-plane intersection point
         if(ze1.gt.0.0 .and. ze2.lt.0.0) then
           alf1 = max(alf1, ze1/(ze1-ze2))
          elseif(ze1.lt.0.0 .and. ze2.gt.0.0) then
           alf2 = min(alf2, ze1/(ze1-ze2))
         endif
c
c
c--- skip to next triangle if no part is hidden
         if(alf1.ge.alf2) go to 1
         if(alf1 .lt.     eps) alf1 = 0.0
         if(alf2 .gt. 1.0-eps) alf2 = 1.0
c
c....... combine new hidden part with previous visible bits
         nseg2 = nseg
c
         do ns = nseg2, 1, -1
           if(alf1.ge.alfseg(2,ns) .or. alf2.le.alfseg(1,ns)) then
           else if(alf1.gt.alfseg(1,ns) .and. alf2.lt.alfseg(2,ns)) then
            if(nseg.eq.nmxseg) then
             write(*,*) '** HIDLIN: Too many segments'
             go to 100
            endif
            nseg = nseg+1
            alfseg(1,nseg) = alf2
            alfseg(2,nseg) = alfseg(2,ns)
            alfseg(2,ns) = alf1
           else if(alf1.le.alfseg(1,ns) .and. alf2.ge.alfseg(2,ns)) then
            alfseg(1,ns) = alfseg(1,nseg)
            alfseg(2,ns) = alfseg(2,nseg)
            nseg = nseg-1
           else if(alf1.le.alfseg(1,ns) .and. alf2.lt.alfseg(2,ns)) then
            alfseg(1,ns) = alf2
           else if(alf1.gt.alfseg(1,ns) .and. alf2.ge.alfseg(2,ns)) then
            alfseg(2,ns) = alf1
           endif
         end do
c
         if(nseg.eq.0) return
 1     continue
c
 100   continue
       return
       end ! hidlin


       subroutine hidpnt(xyzpt,id,ngrp,ilstgrp,xyzgrp,ntri,xyztri,lvis)
c
c      this routine determines the visibility of a single point.
c      given a set of triangles it computes visibility of an input point
c
c      xyzpt = test point coordinate array (x,y,z)
c      id = id index for input point 
c           may be used to designate triangle or collection of triangles 
c           that contain this point and should be skipped for visibility 
c           tests.  if id=0 no triangles will be skipped.
c
c      ngrp = number of triangle groups 
c             (if ngrp<=0 no group minmax tests will be done)
c      ilstgrp(n) = index of last triangle in group n
c      xyzgrp(.n) = x,y group limits
c             1   = min(x)
c             2   = min(y)
c             3   = min(z)
c             4   = max(x)
c             5   = max(y)
c             6   = max(z)
c
c      ntri = number of triangles (used if ngrp=0)
c      xyztri = array of triangle data
c          1  = x1
c          2  = y1
c          3  = z1
c          4  = x2
c          5  = y2
c          6  = z2
c          7  = x3
c          8  = y3
c          9  = z3
c         10  = min(x)
c         11  = min(y)
c         12  = min(z)
c         13  = max(x)
c         14  = max(y)
c         15  = max(z)
c         16  = id field, id associated with this triangle 
c
c      lvis = logical flag  lvis=.true.  if the point is visible
c
       real xyzpt(3), xyztri(16,*), xyzgrp(6,ngrp)
       integer ilstgrp(ngrp)
       logical lvis
c
       data eps/1.0e-7/
c
c----- initialize a few things
       x1 = xyzpt(1)
       y1 = xyzpt(2)
       z1 = xyzpt(3)
       lvis=.true.
c
c----- main processing loop
       klast = 0
       nngrp = max(1,ngrp)
c
       do 100 igrp=1, nngrp
c
       kt1 = klast + 1
       klast = ntri
c
       if(ngrp.ge.1) then
c--- check for group with no triangles (last index <= 0)
        if(ilstgrp(igrp).le.0) go to 100
        klast = ilstgrp(igrp)       
c--- preliminary purge of whole group
        if( (xyzgrp(1,igrp).ge.x1 .or. 
     &       xyzgrp(4,igrp).le.x1) .or.
     &      (xyzgrp(2,igrp).ge.y1 .or. 
     &       xyzgrp(5,igrp).le.y1) ) go to 100
c     &      (xyzgrp(6,igrp).le.z1) ) go to 100
       endif
c
       do 1 kt=kt1, klast
c
c--- second purge of most cells
         if( (xyztri(10,kt).ge.x1 .or. 
     &        xyztri(13,kt).le.x1)   .or.
     &       (xyztri(11,kt).ge.y1 .or. 
     &        xyztri(14,kt).le.y1)   .or.
     &       (xyztri(15,kt).le.z1) ) go to 1
c
c--- if id of this point matches a triangle id skip visibility tests
         if(id.ne.0) then
          idtri = ifix(xyztri(16,kt))
          if(id.eq.idtri) go to 1
         endif
c
         nsame=0
         do l=1, 3
           if( abs(x1-xyztri(1+3*(l-1),kt)).lt.eps .and.
     &         abs(y1-xyztri(2+3*(l-1),kt)).lt.eps .and.
     &         abs(z1-xyztri(3+3*(l-1),kt)).lt.eps )  nsame=nsame+1
         end do
         if(nsame.ge.1) go to 1
c
         n1 = 1
c--- set up temporary variables in canonical orientation where the line
c    passes through sides 1-2 and 1-3 of triangle (point 1 is common point)
         n2 = mod(n1,3) + 1           
         n3 = mod(n2,3) + 1           
c
         xt1 = xyztri(3*n1-2,kt)
         yt1 = xyztri(3*n1-1,kt)
         zt1 = xyztri(3*n1  ,kt)
         xt2 = xyztri(3*n2-2,kt)
         yt2 = xyztri(3*n2-1,kt)
         zt2 = xyztri(3*n2  ,kt)
         xt3 = xyztri(3*n3-2,kt)
         yt3 = xyztri(3*n3-1,kt)
         zt3 = xyztri(3*n3  ,kt)
c
c--- compute xi, eta, zeta values for point in triangle coords
         det = ( (xt2-xt1)*(yt3-yt1) - (xt3-xt1)*(yt2-yt1) )
         if(det.eq.0.0) go to 1
c
         detinv = 1.0 / det
         xi1 = detinv * ( (x1 -xt1)*(yt3-yt1) - (xt3-xt1)*(y1 -yt1) )
         et1 = detinv * ( (xt2-xt1)*(y1 -yt1) - (x1 -xt1)*(yt2-yt1) )
c
         if(abs(xi1-1.0).lt.eps) xi1 = 1.0
         if(abs(et1-1.0).lt.eps) et1 = 1.0
         ze1 = (z1-zt1) - xi1*(zt2-zt1) - et1*(zt3-zt1)
c
c---point is visible if it is outside triangle
         if(xi1.le.0. .or. et1.le.0. .or. (xi1+et1).gt.1.0) go to 1
c---or above it...
         if(ze1.gt.0.) go to 1
c
         lvis=.false.
         return
 1     continue
c
 100   continue
       return
       end ! hidpnt


