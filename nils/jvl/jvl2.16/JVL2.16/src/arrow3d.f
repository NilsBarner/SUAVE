
      subroutine arwset(r,sn,slen,hlen,rhead,nhead,rlines,nlines)
      real r(3), sn(3), rlines(3,2,*)
c------------------------------------------------------------------------
c     Creates a 3d wireframe arrow with a conical arrowhead.
c     the last rlines dimension must be at least 2*nhead+1
c
c     Input:
c        r(.)   location of arrow tail
c        sn(.)  unit vector along arrow
c        slen   length of arrowshaft 
c        hlen   length of arrowhead  (total length is slen+hlen)
c        rhead  arrowhead radius
c        nhead  number of "pie slices" forming conical arrowhead
c
c     Output:
c        rlines(.1s)  first  point of line segment
c        rlines(.2s)  second point of line segment
c        nlines       number of line segments  ( = 2*nhead+1 )
c
c        The last segment (s=nlines) is the arrowshaft
c------------------------------------------------------------------------
      real a(3), b(3)
      data pi / 3.1415926 /
c
      kmin = 1
      if(abs(sn(2)) .lt. abs(sn(kmin))) kmin = 2
      if(abs(sn(3)) .lt. abs(sn(kmin))) kmin = 3
c
      b(1) = 0.
      b(2) = 0.
      b(3) = 0.
      b(kmin) = 1.0
c
      call cross(b,sn,a)
      amag = sqrt(a(1)**2 + a(2)**2 + a(3)**2)
      if(amag .gt. 0.0) then
       a(1) = rhead*a(1)/amag
       a(2) = rhead*a(2)/amag
       a(3) = rhead*a(3)/amag
      endif
c
      call cross(a,sn,b)
c
c---- set arrowhead perimeter lin segments
      il = 1
      rlines(1,1,il) = r(1) + a(1)
      rlines(2,1,il) = r(2) + a(2)
      rlines(3,1,il) = r(3) + a(3)
c
      nt = nhead
      dt = 2.0*pi / float(nt)
      do it = 1, nt-1
        t = dt * float(it)
        cost = cos(t)
        sint = sin(t)
c
        il = it
        rlines(1,2,il) = r(1) + a(1)*cost + b(1)*sint
        rlines(2,2,il) = r(2) + a(2)*cost + b(2)*sint
        rlines(3,2,il) = r(3) + a(3)*cost + b(3)*sint
c
        rlines(1,1,il+1) = rlines(1,2,il)
        rlines(2,1,il+1) = rlines(2,2,il)
        rlines(3,1,il+1) = rlines(3,2,il)
      enddo
c
      il = nt
      rlines(1,2,il) = rlines(1,1,1)
      rlines(2,2,il) = rlines(2,1,1)
      rlines(3,2,il) = rlines(3,1,1)
      nlines = nt
c
c---- set arrowhead radial segments
      do it = 1, nt
        il = nt + it
c
c------ from point on perimeter ...
        rlines(1,1,il) = rlines(1,1,it)
        rlines(2,1,il) = rlines(2,1,it)
        rlines(3,1,il) = rlines(3,1,it)
c
c------  ... to arrowhead point
        rlines(1,2,il) = r(1) + sn(1)*hlen
        rlines(2,2,il) = r(2) + sn(2)*hlen
        rlines(3,2,il) = r(3) + sn(3)*hlen
      enddo
c
      nlines = 2*nt
c
c---------------------------------------------------
c---- add arrowshaft
c
c---- move entire arrowhead by shaft distance
      do il = 1, nlines
        do k = 1, 3
          rlines(k,1,il) = rlines(k,1,il) + slen*sn(k)
          rlines(k,2,il) = rlines(k,2,il) + slen*sn(k)
        enddo
      enddo
c
c---- add shaft segment
      nlines = nlines + 1
      il = nlines
c
      rlines(1,1,il) = r(1)
      rlines(2,1,il) = r(2)
      rlines(3,1,il) = r(3)
c
      rlines(1,2,il) = r(1) + slen*sn(1)
      rlines(2,2,il) = r(2) + slen*sn(2)
      rlines(3,2,il) = r(3) + slen*sn(3)
c
      return
      end ! arwset



      subroutine arwplt(xo,yo,afac,rlines,sn,rhead,nhead)
c---------------------------------------------------------
c     plots the 3d arrow created by arwset, presumably
c     after it was projected into screen x,y coordinates.
c     the z coordinates are used here for hidden-line 
c     logic to give the arrow a "solid" appearance.
c---------------------------------------------------------
      real rlines(3,2,*), sn(3)
c
      real dlhead(3)
c
c---- end of arrowshaft
      nl = 2*nhead + 1
      zshaft = rlines(3,2,nl)
c
c---- arrowshaft length
      sshaft = sqrt(  (rlines(1,2,nl)-rlines(1,1,nl))**2
     &              + (rlines(2,2,nl)-rlines(2,1,nl))**2
     &              + (rlines(3,2,nl)-rlines(3,1,nl))**2 )
c
c---- normalize direction vector just in case it was deformed by perspective
      snsq = sn(1)**2 + sn(2)**2 + sn(3)**2
      if(snsq .gt. 0.0) then
       sni = 1.0/sqrt(snsq)
      else
       sni = 1.0
      endif
c
c---- plot arrowhead perimiter
      do ih = 1, nhead
        il = ih
c------ midpoint of perimter segment
        zmid = (rlines(3,2,il) + rlines(3,1,il))*0.5
c
c------ plot segment only if arrow points down (z component into screen),
c        or arrowshaft end's z location is below segment's z midpoint
        if(sn(3) .lt. 0.0 .or. zshaft .lt. zmid) then
         x1 = xo + afac*rlines(1,1,il)
         y1 = yo + afac*rlines(2,1,il)
         x2 = xo + afac*rlines(1,2,il)
         y2 = yo + afac*rlines(2,2,il)
         call plot(x1,y1,3)
         call plot(x2,y2,2)
        endif
      enddo
c
c---- plot arrowhead radial segments
      do ih = 1, nhead
        il = nhead + ih
c
c------ set unit vector along radial segment
        dlhead(1) = rlines(1,2,il) - rlines(1,1,il)
        dlhead(2) = rlines(2,2,il) - rlines(2,1,il)
        dlhead(3) = rlines(3,2,il) - rlines(3,1,il)
        slhead = sqrt(dlhead(1)**2 + dlhead(2)**2 + dlhead(3)**2)
        if(slhead .gt. 0.0) then
         dlhead(1) = dlhead(1)/slhead
         dlhead(2) = dlhead(2)/slhead
         dlhead(3) = dlhead(3)/slhead
        endif
c
c------ plot segment only if it points more into screen than arrowshaft
        if(dlhead(3) .lt. sn(3)*sni) then
         x1 = xo + afac*rlines(1,1,il)
         y1 = yo + afac*rlines(2,1,il)
         x2 = xo + afac*rlines(1,2,il)
         y2 = yo + afac*rlines(2,2,il)
         call plot(x1,y1,3)
         call plot(x2,y2,2)
        endif
      enddo
c
c---- don't bother plotting zero-length arrowshaft
      if(sshaft .eq. 0.0) return
c
c---- find visible fraction of arrowshaft
      sn12 = sqrt(sn(1)**2 + sn(2)**2)
      if(sn12 .eq. 0.0) then
       tant = 0.0
      else
       tant = sn(3) / sn12
      endif
      sseen = max( sshaft - rhead*max( tant , 0.0 ) , 0.0 )
c
      if(sseen .eq. 0.0) return
c
c---- plot visible fraction of arrowshaft
      sfrac = sseen/sshaft
c
      il = nl
      x1 = xo + afac* rlines(1,1,il)
      y1 = yo + afac* rlines(2,1,il)
      x2 = x1 + afac*(rlines(1,2,il)-rlines(1,1,il))*sfrac
      y2 = y1 + afac*(rlines(2,2,il)-rlines(2,1,il))*sfrac
      call plot(x1,y1,3)
      call plot(x2,y2,2)
c
      return
      end ! arwplt

