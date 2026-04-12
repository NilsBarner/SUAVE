
      subroutine vorvelc(x,y,z,lbound,x1,y1,z1,x2,y2,z2,beta,
     &                   u,v,w, rcore)
c----------------------------------------------------------
c     Same as vorvel, with finite core radius
C     Original Scully (AKA Burnham-Hallock) core model 
C       Vtan = Gam/2*pi . r/(r^2 +rcore^2)
C      
C     Uses Leishman's R^4 variant of Scully (AKA Burnham-Hallock) core model 
C       Vtan = Gam/2*pi . r/sqrt(r^4 +rcore^4)
C----------------------------------------------------------
      logical lbound
c
c
      real a(3), b(3), axb(3)
c
c---- 1 / (4 pi)
      data pi4inv / 0.079577471545947667884441881686257184 /

c---- prandtl-glauert coordinates 
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
      rcore2 = rcore**2
      rcore4 = rcore2**2
c
      u = 0.
      v = 0.
      w = 0.
c
c---- contribution from the transverse bound leg
      if (lbound  .and.  amag*bmag .ne. 0.0) then
        axb(1) = a(2)*b(3) - a(3)*b(2)
        axb(2) = a(3)*b(1) - a(1)*b(3)
        axb(3) = a(1)*b(2) - a(2)*b(1)
        axbsq = axb(1)**2 + axb(2)**2 + axb(3)**2
c
        if(axbsq .ne. 0.0) then
         adb = a(1)*b(1) + a(2)*b(2) + a(3)*b(3)
         alsq = asq + bsq - 2.0*adb
ccc      rsq = axbsq / alsq
c     
         abmag = amag*bmag
c---- scully core model      
cccc        t = (amag+bmag)*(1.0 - adb/abmag) / (axbsq + alsq*rcore2)
cc        t = (  (bsq-adb)/sqrt(bsq+rcore2)
cc     &       + (asq-adb)/sqrt(asq+rcore2) ) / (axbsq + alsq*rcore2)
c---- leishman core model
         t = (  (bsq-adb)/sqrt(sqrt(bsq**2+rcore4))
     &        + (asq-adb)/sqrt(sqrt(asq**2+rcore4)) )
     &        / sqrt(axbsq**2 + alsq**2*rcore4)
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
        adx = a(1)
        rsq = axisq
c
c---- scully core model      
cc        t = - (1.0 - adx/amag) / (rsq + rcore2)
c---- leishman core model
        t = - (1.0 - adx/amag) / sqrt(rsq**2 + rcore4)
c
        v = v + a(3)*t
        w = w - a(2)*t
      endif
c
c---- trailing leg attached to b
      if (bmag .ne. 0.0) then
        bxisq = b(3)**2 + b(2)**2
        bdx = b(1)
        rsq = bxisq
c
c---- scully core model      
cc        t =   (1.0 - bdx/bmag) / (rsq + rcore2)
c---- leishman core modeld
        t =   (1.0 - bdx/bmag) / sqrt(rsq**2 + rcore4)
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
      end ! vorvelc

