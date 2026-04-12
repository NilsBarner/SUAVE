C***********************************************************************
C    Module:  matrix.f
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

      subroutine ludcmp(nsiz,n,a)
c     *******************************************************
c     *                                                     *
c     *   factors a full nxn matrix into an lu form.        *
c     *   subr. baksub can back-substitute it with some rhs.*
c     *   assumes matrix is non-singular...                 *
c     *    ...if it isn't, a divide by zero will result.    *
c     *                                                     *
c     *   a is the matrix...                                *
c     *     ...replaced with its lu factors.                *
c     *                                                     *
c     *                              mark drela  1988       *
c     *******************************************************
      implicit none

      integer nsiz,n
      real a(nsiz,nsiz)

      integer i,j,k
      real sum, dum
c
      do j=1, n
        do i=1, j-1
          sum = a(i,j)
          do k=1, i-1
            sum = sum - a(i,k)*a(k,j)
          enddo
          a(i,j) = sum
        enddo
c
        do i=j, n
          sum = a(i,j)
          do k=1, j-1
            sum = sum - a(i,k)*a(k,j)
          enddo
          a(i,j) = sum
        enddo
c
        dum = 1.0/a(j,j)
        do i=j+1, n
          a(i,j) = a(i,j)*dum
        enddo
      enddo
c
      return
      end ! ludcmp



      subroutine baksub(nsiz,n,a,b)
      implicit none

      integer nsiz,n
      real a(nsiz,nsiz), b(nsiz)
c
      integer i,j,k,ii
      real sum, dum

      ii = 0
      do i=1, n
        sum = b(i)
        if(ii.ne.0) then
         do j=ii, i-1
           sum = sum - a(i,j)*b(j)
         enddo
        else if(sum.ne.0.0) then
         ii = i
        endif
        b(i) = sum
      enddo
c
      do i=n, 1, -1
        sum = b(i)
        if(i.lt.n) then
         do j=i+1, n
           sum = sum - a(i,j)*b(j)
         enddo
        endif
        b(i) = sum/a(i,i)
      enddo
c
      return
      end ! baksub



      subroutine ludcmpi(nsiz,n,a,indx,work)
c     *******************************************************
c     *   factors a full nxn matrix into an lu form.        *
c     *   subr. baksub can back-substitute it with some rhs.*
c     *   assumes matrix is non-singular...                 *
c     *    ...if it isn't, a divide by zero will result.    *
c     *                                                     *
c     *   a is the matrix...                                *
c     *     ...replaced with its lu factors.                *
c     *                                                     *
c     *   stolen from numerical recipes.                    *
c     *******************************************************
c
      real a(nsiz,nsiz), work(nsiz)
      integer indx(nsiz)
c
      do 12 i=1, n
        aamax = 0.
        do 11 j=1, n
          aamax = max( abs(a(i,j)) , aamax )
   11   continue
        work(i) = 1.0/aamax
   12 continue
c
      do 19 j=1, n
        do 14 i=1, j-1
          sum = a(i,j)
          do 13 k=1, i-1
            sum = sum - a(i,k)*a(k,j)
   13     continue
          a(i,j) = sum
   14   continue
c
        aamax = 0.
        do 16 i=j, n
          sum = a(i,j)
          do 15 k=1, j-1
            sum = sum - a(i,k)*a(k,j)
   15     continue
          a(i,j) = sum
c
          dum = work(i)*abs(sum)
          if(dum.ge.aamax) then
           imax = i
           aamax = dum
          endif
   16   continue
c
        if(j.ne.imax) then
         do 17 k=1, n
           dum = a(imax,k)
           a(imax,k) = a(j,k)
           a(j,k) = dum
   17    continue
         work(imax) = work(j)
        endif
c
        indx(j) = imax
        if(j.ne.n) then
         dum = 1.0/a(j,j)
         do 18 i=j+1, n
           a(i,j) = a(i,j)*dum
   18    continue
        endif
c
   19 continue
c
      return
      end ! ludcmpi


      subroutine baksubi(nsiz,n,a,indx,b)
      real a(nsiz,nsiz), b(nsiz)
      integer indx(nsiz)
c
      ii = 0
      do 12 i=1, n
        ll = indx(i)
        sum = b(ll)
        b(ll) = b(i)
        if(ii.ne.0) then
         do 11 j=ii, i-1
           sum = sum - a(i,j)*b(j)
   11    continue
        else if(sum.ne.0.0) then
         ii = i
        endif
        b(i) = sum
   12 continue
c
      do 14 i=n, 1, -1
        sum = b(i)
        if(i.lt.n) then
         do 13 j=i+1, n
           sum = sum - a(i,j)*b(j)
   13    continue
        endif
        b(i) = sum/a(i,i)
   14 continue
c
      return
      end ! baksubi


