C***********************************************************************
C    Module:  airutil.f
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

      subroutine readbl(fname,ibx,nbx,xb,yb,iib,nbl,
     &                  name,xinl,xout,ybot,ytop)
      character*(*) fname
      real xb(ibx,nbx) ,yb(ibx,nbx)
      integer iib(nbx)
      character*(*) name
c----------------------------------------------
c     Reads in blade.xxx dataset
c----------------------------------------------
      character*80 line
      logical error

      real ainput(20)

      character*19 nchars
      data nchars / '0123456789-+.,edED ' /

c---- default top/bottom flow area ratio for "old" blade.xxx files
      data yrat / 1.3 /

c---- first assume that there will be a read error
      nbl = 0

      nf = index(fname,' ') + 1

      lu = 3
      open(lu,file=fname,status='old',err=98)

      read(lu,1000) line

c---- if first line has any non-numeric character, go treat it as the name
      do k=1, 80
        if(index(nchars,line(k:k)) .eq. 0) go to 20
      enddo

c---- plain unlabeled file: rewind, and just read in x,y coordinates
      name = ' '
      ninput = 0

      rewind(lu)
cc      write(*,*)
cc      write(*,*) 'Reading plain coordinate file'
      go to 40

c---- first line interpreted as label string
   20 read(line,1000) name

c---- read and decode second line --- grid domain limits will be set later
      read(lu,1000) line
      ninput = 4
      call getflt(line,ainput,ninput,error)

      if(error) go to 99

      if(ninput.lt.4) then
c------ no domain parameters: re-read name string and then read x,y coordinates
        rewind(lu)
        read(lu,1000) line
      endif

cc      write(*,1010) name
cc 1010 format(/1x,'reading in coordinate file for: ',a/)

 40   continue

c---- read in airfoil coordinates
      do 55 n=1, nbx+1
        ib = 1
 50     continue
          read(lu,*,end=56,err=99) xbt, ybt
          if(xbt.eq.999.0) then
           iib(n) = ib-1
           go to 55
          endif
          if(n.gt.nbx) then
           write(*,*) '*** readbl: too many elements. increase nbx to',n
           stop
          endif
          iblim = min( ib , ibx )
          xb(iblim,n) = xbt
          yb(iblim,n) = ybt

          ib = ib + 1
          go to 50
   55 continue
      n = nbx

   56 continue
      if(ib.eq.1) then
c----- coordinate file has "999.0 999.0" at the end ...
       nbl = n-1
      else
c----- coordinate file has no ending line (single element ises file)
       nbl = n
       iib(n) = ib-1
      endif

      close(lu)


      do 80 n = 1, nbl
        if(iib(n).gt.ibx) then
         write(*,*)
     &    '*** READBL: Too many airfoil points. Increase ibx to', iib(n)
         stop
        endif

c------ calculate airfoil element area
        area = 0.
        do 802 ib=1, iib(n)-1
          rx = xb(ib+1,n) + xb(ib,n)
          ry = yb(ib+1,n) + yb(ib,n)
          dx = xb(ib+1,n) - xb(ib,n)
          dy = yb(ib+1,n) - yb(ib,n)
          da = 0.25*(rx*dy - ry*dx)
          area = area + da
 802    continue

        if(area.lt.0.0) then
c------- if area is negative (clockwise order), reverse coordinate order
         do 804 ib=1, iib(n)/2
           iback = iib(n) - ib + 1
           xtmp = xb(ib,n)
           ytmp = yb(ib,n)
           xb(ib,n) = xb(iback,n)
           yb(ib,n) = yb(iback,n)
           xb(iback,n) = xtmp
           yb(iback,n) = ytmp
 804     continue
        endif

 80   continue

      if     (ninput.lt.4) then
c------ plain or labeled airfoil file -- no domain parameters specified
        xinl = 0.
        xout = 0.
        ybot = 0.
        ytop = 0.
      elseif (ninput.eq.4) then
c------ for "new" blade.xxx file, grid size is input directly
        xinl = ainput(1)
        xout = ainput(2)
        ybot = ainput(3)
        ytop = ainput(4)
      else
c------ "old" blade.xxx file, grid size is implied from airfoil limits
        chinl = ainput(3)
        chout = ainput(4)
        chwid = ainput(5)
c
        xmin = xb(1,1) + 1.0
        xmax = xb(1,1) - 1.0
        do 84 n=1, nbl
          do 842 ib=1, iib(n)
            if(xb(ib,n).le.xmin) then
             xmin = xb(ib,n)
             ymin = yb(ib,n)
            endif
            if(xb(ib,n).ge.xmax) then
             xmax = xb(ib,n)
             ymax = yb(ib,n)
            endif
 842      continue
 84     continue

        xinl = xmin - chinl
        xout = xmax + chout
        ybot = ymin - chwid *  1.0/(1.0 + yrat)
        ytop = ymax + chwid * yrat/(1.0 + yrat)
      endif


      return

   98 continue
      write(*,1050) fname(1:nf)
      nbl = 0
      return

   99 continue
      write(*,1100) fname(1:nf)
      return
c...............................................................
 1000 format(a)
 1050 format(/' file open error:  ', a)
 1100 format(/' file read error:  ', a)
      end ! readbl



      subroutine getcam(x,y,n,xc,yc,tc,nc,lnorm)
c--------------------------------------------------------
c     Takes airfoil x,y surface points and returns 
c     the camber defined in xc,yc at nc points
c--------------------------------------------------------
      real x(*), y(*), xc(*), yc(*), tc(*)
      logical lnorm

      parameter(nsiz=300)
      real xp(nsiz), yp(nsiz), s(nsiz)

      pi = 4.0*atan(1.0)
      if(n.gt.nsiz) then
       write(*,*) '*** getcam: array overflow. increase nsiz to', n
       stop
      endif

c---- spline coordinates
      call scalc(x,y,s,n)
      call segspl(x,xp,s,n)
      call segspl(y,yp,s,n)

c---- find arc length position of leading edge
      call lefind(sle,x,xp,y,yp,s,n)

c---- normalize airfoil and its spline data to unit chord
      if(lnorm) call normit(sle,x,xp,y,yp,s,n)

      xle = seval(sle,x,xp,s,n)
      yle = seval(sle,y,yp,s,n)
      xte = 0.5*(x(1)+x(n))
      yte = 0.5*(y(1)+y(n))

c---- number of output points defaults to 30
      if(nc.le.0) nc = 30

      su = sle - 0.01
      sl = sle + 0.01
      fnc1 = float(nc-1)
      xc(1) = xle
      yc(1) = yle
      tc(1) = 0.0
      do i = 2, nc
        xout = xle + (xte-xle)*0.5*(1.0 - cos(pi*float(i-1)/fnc1))
        call sinvrt(su,xout,x,xp,s,n)
        yu = seval(su,y,yp,s,n)
        call sinvrt(sl,xout,x,xp,s,n)
        yl = seval(sl,y,yp,s,n)
        xc(i) = xout
        yc(i) = 0.5*(yu+yl)
        tc(i) = yu-yl
      end do

      return
      end ! getcam


      subroutine lefind(sle,x,xp,y,yp,s,n)
c------------------------------------------------
c     finds the spline parameter value sle
c     at the leftmost point of the airfoil
c     (i.e. the leading edge)
c------------------------------------------------
      real x(n),y(n),xp(n),yp(n),s(n)

c---- initial guess for sle = leftmost point
      do i = 2, n
        if(x(i).gt.x(i-1)) go to 6
      end do
    6 sle = s(i-1)

c---- newton solution for the exact sle value (at least to within machine zero)
      sref = s(n) - s(1)
      do iter = 1, 20
        res  = deval(sle,x,xp,s,n)
        resp = d2val(sle,x,xp,s,n)
        dsle = -res/resp
        sle = sle + dsle
        if(abs(dsle)/sref .lt. 1.0e-5) return
      end do
      write(*,*) '** lefind: leading edge not found.  continuing...'
      sle = s(i-1)
      return
      end ! lefind
 
 

 
      subroutine normit(sle,x,xp,y,yp,s,n)
c--------------------------------------------------------
c     Normalizes airfoil coordinates and their spline
c     parameter to unit chord. the x coordinates are
c     also offset so that the leading edge is at x = 0.
c--------------------------------------------------------
      real x(n),y(n),xp(n),yp(n),s(n)

c---- leading edge and trailing edge x coordinates
      xle = seval(sle,x,xp,s,n)
      xte = 0.5*(x(1)+x(n))

c---- normalizing factor
      dnorm = 1.0/(xte-xle)
      do i = 1, n
        x(i) = (x(i)-xle)*dnorm
        y(i) = y(i)*dnorm
        s(i) = s(i)*dnorm
      end do
      sle = sle*dnorm

      return
      end ! normit
