C***********************************************************************
C    Module:  jutil.f
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

      subroutine m3inv(a,ainv)
      real a(3,3), ainv(3,3)
c----------------------------------------------
c     Calculates inverse of 3x3 matrix with 
c     special treatment...
c
c       If a diagonal element is exactly zero, 
c       then it is interpreted as infinity.
c----------------------------------------------
      real t(3,6)
c
c---- set augmented matrix
      do k=1, 3
        t(k,1) = a(k,1)
        t(k,2) = a(k,2)
        t(k,3) = a(k,3)
        t(k,4) = 0.
        t(k,5) = 0.
        t(k,6) = 0.
      enddo
      t(1,4) = 1.0
      t(2,5) = 1.0
      t(3,6) = 1.0
c
      do n=1, 3
        pivot = t(n,n)
c
c------ normalize pivot row...
        if(pivot .eq. 0.0) then
c-------- interpret pivot as infinity
          do l=n+1, 6
            t(n,l) = 0.0
          enddo
        else
c-------- interpret pivot as-is
          do l=n+1, 6
            t(n,l) = t(n,l)/pivot
          enddo
        endif
c
c------ eliminate pivot column
        do k=n+1, 3
          tel = t(k,n)
          do l=n+1, 6
            t(k,l) = t(k,l) - tel*t(n,l)
          enddo
        enddo
      enddo
c
c---- back-substitute
      do n=3, 2, -1
        do k=n-1, 1, -1
          tel = t(k,n)
          do l=4, 6
            t(k,l) = t(k,l) - tel*t(n,l)
          enddo
        enddo
      enddo
c
      do k=1, 3
        ainv(k,1) = t(k,4)
        ainv(k,2) = t(k,5)
        ainv(k,3) = t(k,6)
      enddo
c
      return
      end ! m3inv


      subroutine m6inv(a,ainv)
      real a(6,6), ainv(6,6)
c----------------------------------------------
c     Calculates inverse of 6x6 matrix with 
c     special treatment...
c
c       If a diagonal element is exactly zero, 
c       then it is interpreted as infinity.
c----------------------------------------------
      real t(6,12)
c
c---- set augmented matrix
      do k = 1, 6
        do l = 1, 6
          t(k,l) = a(k,l)
          t(k,l+6) = 0.
        enddo
        t(k,k+6) = 1.0
      enddo

      do n = 1, 6
        pivot = t(n,n)
c
c------ normalize pivot row...
        if(pivot .eq. 0.0) then
c-------- interpret pivot as infinity
          do l = n+1, 12
            t(n,l) = 0.0
          enddo
        else
c-------- interpret pivot as-is
          do l = n+1, 12
            t(n,l) = t(n,l)/pivot
          enddo
        endif
c
c------ eliminate pivot column
        do k = n+1, 6
          tel = t(k,n)
          do l = n+1, 12
            t(k,l) = t(k,l) - tel*t(n,l)
          enddo
        enddo
      enddo
c
c---- back-substitute
      do n = 6, 2, -1
        do k = n-1, 1, -1
          tel = t(k,n)
          do l = 7, 12
            t(k,l) = t(k,l) - tel*t(n,l)
          enddo
        enddo
      enddo
c
      do k = 1, 6
        do l = 1, 6
          ainv(k,l) = t(k,l+6)
        enddo
      enddo
c
      return
      end ! m6inv



      subroutine rateki3(a, r,r_a)
      real a(3)
      real r(3,3), r_a(3,3,3)
c-----------------------------------------------------
c                                    t
c     Calculates inverse of tensor  T K
c
c                            t
c     d{phi theta psi}  =  [T K] d{p q r}
c
c-----------------------------------------------------
      c1 = cos(a(1))
      c2 = cos(a(2))
c
      s1 = sin(a(1))
c
      t2 = tan(a(2))
c
      r(1,1)     = -1.0
      r(2,1)     = 0.
      r(3,1)     = 0.
c
      r(1,2)     =  s1*t2
      r(2,2)     =  c1
      r(3,2)     =  s1/c2
c
      r(1,3)     = -c1*t2
      r(2,3)     =  s1
      r(3,3)     = -c1/c2
c
      r_a(1,1,1) = 0.
      r_a(2,1,1) = 0.
      r_a(3,1,1) = 0.
c
      r_a(1,2,1) =  c1*t2
      r_a(2,2,1) = -s1
      r_a(3,2,1) =  c1/c2
c
      r_a(1,3,1) =  s1*t2
      r_a(2,3,1) =  c1
      r_a(3,3,1) =  s1/c2
c
      r_a(1,1,2) = 0.
      r_a(2,1,2) = 0.
      r_a(3,1,2) = 0.
c
      r_a(1,2,2) =  s1/c2**2
      r_a(2,2,2) = 0.
      r_a(3,2,2) =  s1*t2/c2
c
      r_a(1,3,2) = -c1/c2**2
      r_a(2,3,2) = 0.
      r_a(3,3,2) = -c1*t2/c2
c
      r_a(1,1,3) = 0.
      r_a(2,1,3) = 0.
      r_a(3,1,3) = 0.
c
      r_a(1,2,3) = 0.
      r_a(2,2,3) = 0.
      r_a(3,2,3) = 0.
c
      r_a(1,3,3) = 0.
      r_a(2,3,3) = 0.
      r_a(3,3,3) = 0.
c
      return
      end ! rateki3



      subroutine rotens3(a,t,t_a)
      real a(3), t(3,3), t_a(3,3,3)
c------------------------------------------
c     Calculates net rotation tensor T(..) 
c     from the three rotation angles a(.).
c
c     The order and sense of rotation from
c     x axes to x' axes is
c
c      x'  =  a3 [ -a2 [ a1 [ x ] ] ]
c      x'  =  T x
c
c     Inputs:
c       a(.)   phi, theta, psi (radians)
c
c     Outputs:
c       t(..)      T matrix
c       t_a(..1)   dT/dphi
c       t_a(..2)   dT/dtheta
c       t_a(..3)   dT/dpsi
c------------------------------------------
c
      c1 = cos(a(1))   !  phi
      c2 = cos(a(2))   !  theta
      c3 = cos(a(3))   !  psi
c
      s1 = sin(a(1))
      s2 = sin(a(2))
      s3 = sin(a(3))
c
      t(1,1) =     c2*c3
      t(2,1) =    -c2*s3
      t(3,1) =    -s2
c
      t(1,2) = -s1*s2*c3 + c1*s3
      t(2,2) =  s1*s2*s3 + c1*c3
      t(3,2) = -s1*c2
c
      t(1,3) =  c1*s2*c3 + s1*s3
      t(2,3) = -c1*s2*s3 + s1*c3
      t(3,3) =  c1*c2
c
c
      t_a(1,1,1) = 0.
      t_a(2,1,1) = 0.
      t_a(3,1,1) = 0.
c
      t_a(1,2,1) = -c1*s2*c3 - s1*s3
      t_a(2,2,1) =  c1*s2*s3 - s1*c3
      t_a(3,2,1) = -c1*c2
c
      t_a(1,3,1) = -s1*s2*c3 + c1*s3
      t_a(2,3,1) =  s1*s2*s3 + c1*c3
      t_a(3,3,1) = -s1*c2
c
c
      t_a(1,1,2) =    -s2*c3
      t_a(2,1,2) =     s2*s3
      t_a(3,1,2) =    -c2
c
      t_a(1,2,2) = -s1*c2*c3
      t_a(2,2,2) =  s1*c2*s3
      t_a(3,2,2) =  s1*s2
c
      t_a(1,3,2) =  c1*c2*c3
      t_a(2,3,2) = -c1*c2*s3
      t_a(3,3,2) = -c1*s2
c
c
      t_a(1,1,3) =    -c2*s3
      t_a(2,1,3) =    -c2*c3
      t_a(3,1,3) = 0.
c
      t_a(1,2,3) =  s1*s2*s3 + c1*c3
      t_a(2,2,3) =  s1*s2*c3 - c1*s3
      t_a(3,2,3) = 0.
c
      t_a(1,3,3) = -c1*s2*s3 + s1*c3
      t_a(2,3,3) = -c1*s2*c3 - s1*s3
      t_a(3,3,3) = 0.
c
c
      return
      end ! rotens3



