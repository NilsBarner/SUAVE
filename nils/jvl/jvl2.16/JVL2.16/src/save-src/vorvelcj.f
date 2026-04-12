
      subroutine vorvelc(x,y,z,lbound,x1,y1,z1,x2,y2,z2,beta,
     &                   u,v,w, rcore)
c----------------------------------------------------------
c     Same as vorvel, with finite core radius
c----------------------------------------------------------
      logical lbound

      real a(3), b(3), axb(3)

c---- 1 / (4 pi)
      data pi4inv / 0.079577471545947667884441881686257184 /

      a(1) = (x - x1)/beta
      a(2) =  y - y1
      a(3) =  z - z1

      b(1) = (x - x2)/beta
      b(2) =  y - y2
      b(3) =  z - z2

      asq = a(1)**2 + a(2)**2 + a(3)**2
      bsq = b(1)**2 + b(2)**2 + b(3)**2

      amag = sqrt(asq)
      bmag = sqrt(bsq)

      u = 0.
      v = 0.
      w = 0.

c---- contribution from the transverse bound leg
      if (lbound  .and.  amag*bmag .ne. 0.0) then
        axb(1) = a(2)*b(3) - a(3)*b(2)
        axb(2) = a(3)*b(1) - a(1)*b(3)
        axb(3) = a(1)*b(2) - a(2)*b(1)
        axbsq = axb(1)**2 + axb(2)**2 + axb(3)**2

       if (axbsq .ne. 0.0) then
         adb = a(1)*b(1) + a(2)*b(2) + a(3)*b(3)
         alsq = asq + bsq - 2.0*adb

ccc      rsq = axbsq / alsq

         ab = amag*bmag
c        t = (amag+bmag)*(1.0 - adb/ab) / (axbsq + alsq*rcore**2)
         t = (  (bsq-adb)/sqrt(bsq+rcore**2)
     &        + (asq-adb)/sqrt(asq+rcore**2) ) / (axbsq + alsq*rcore**2)

         u = axb(1)*t
         v = axb(2)*t
         w = axb(3)*t
       endif
      endif

c---- trailing leg attached to a
      if (amag .ne. 0.0) then
        axisq = a(3)**2 + a(2)**2

        adi = a(1)
        rsq = axisq

        t = (1.0 + adi/amag) / (rsq + rcore**2)

        v = v + a(3)*t
        w = w - a(2)*t
      endif

c---- trailing leg attached to b
      if (bmag .ne. 0.0) then
        bxisq = b(3)**2 + b(2)**2

        bdi = b(1)
        rsq = bxisq

        t =  - (1.0 + bdi/bmag) / (rsq + rcore**2)

        v = v + b(3)*t
        w = w - b(2)*t
      endif

      u = u*pi4inv / beta
      v = v*pi4inv 
      w = w*pi4inv 

      return
      end ! vorvelc

