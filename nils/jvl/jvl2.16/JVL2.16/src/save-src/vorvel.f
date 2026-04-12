
      subroutine vorvel(x,y,z,lbound,x1,y1,z1,x2,y2,z2,beta,
     &                   u,v,w )
c----------------------------------------------------------
c     Horseshoe vortex velocity influence, w/o core radius
c----------------------------------------------------------
      logical lbound
c
      real a(3), b(3), axb(3)
c
c---- 1 / (4 pi)
      data pi4inv / 0.079577471545947667884441881686257184 /
c
      a(1) = (x1 - x)/beta
      a(2) =  y1 - y
      a(3) =  z1 - z
c
      b(1) = (x2 - x)/beta
      b(2) =  y2 - y
      b(3) =  z2 - z
c
      asq = a(1)**2 + a(2)**2 + a(3)**2
      bsq = b(1)**2 + b(2)**2 + b(3)**2
c
      amag = sqrt(asq)
      bmag = sqrt(bsq)
c
      u = 0.
      v = 0.
      w = 0.
c
c---- contribution from the transverse bound leg
      if (lbound .and.  amag*bmag .ne. 0.0) then
        axb(1) = a(2)*b(3) - a(3)*b(2)
        axb(2) = a(3)*b(1) - a(1)*b(3)
        axb(3) = a(1)*b(2) - a(2)*b(1)
c
        adb = a(1)*b(1) + a(2)*b(2) + a(3)*b(3)
c
        den = amag*bmag + adb
c
        if(den .ne. 0.0) then
         t = (1.0/amag + 1.0/bmag) / den
c
         u = axb(1)*t
         v = axb(2)*t
         w = axb(3)*t
        endif
      endif
c
c---- trailing leg attached to a
      if (amag .ne. 0.0) then
        axisq = a(3)**2 + a(2)**2
c
        adi = a(1)
        rsq = axisq
c
        t = - (1.0 - adi/amag) / rsq
c
        v = v + a(3)*t
        w = w - a(2)*t
      endif
c
c---- trailing leg attached to b
      if (bmag .ne. 0.0) then
        bxisq = b(3)**2 + b(2)**2
c
        bdi = b(1)
        rsq = bxisq
c
        t =   (1.0 - bdi/bmag) / rsq
c
        v = v + b(3)*t
        w = w - b(2)*t
      endif
c
      u = u*pi4inv / beta
      v = v*pi4inv 
      w = w*pi4inv 
c
      return
      end ! vorvel

