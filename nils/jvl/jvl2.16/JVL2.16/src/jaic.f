C***********************************************************************
C    Module:  jaic.f
C 
C    Copyright (C) 2021 Mark Drela, Harold Youngren
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


      subroutine vvor(betm,iysym,ysym,zsym,vrcorec,vrcorew,
     &     nv,rv1,rv2,ncompv,chordv,izimv,
     &     nc,rc ,    ncompc,lvtest,
     &     vc_gam,ncdim)
c--------------------------------------------------------------------
c     Calculates the velocity influence matrix for a collection 
c     of horseshoe vortices and field points
c     
c Inputs:
c     betm      sqrt(1-mach*mach)
c     iysym     plane of symmetry xz 
c                = 0 no symmetry
c                = 1 regular symmetry
c                =-1 free-surface symmetry
c     ysym      y coordinate of symmetry plane
c     izimv(.)  second plane of symmetry xy for each vortex
c                = 0 no second plane
c                = 1 regular symmetry
c                =-1 free-surface symmetry
c     zsym      z coordinate of symmetry plane
c
c     vrcorec   vortex-line core radius factor (*semichord)
c     vrcorew   vortex-line core radius factor (*vortex width)
c
c     nv        number of vortices
c     rv1(3,v)  coordinates of endpoint #1 of the vortices
c     rv2(3,v)  coordinates of endpoint #2 of the vortices
c     ncompv(v) index of surface containing h.v.
c     chordv(v) chord of strip containing h.v.
c
c     nc        number of field points
c     rc(3,c)   coordinates of the field points
c     ncompc(c) index of surface containing field point
c     lvtest    t if core-radius test is to be applied
c
c     ncdim     declared size of vc_gam matrix
c     
c Outputs:
c     vc_gam(3..)   velocity/gamma influence matrix
c     
c--------------------------------------------------------------------
      real rv1(3,nv),
     &     rv2(3,nv),
     &     chordv(nv)
      integer izimv(nv)
      real rc(3,nc),
     &     vc_gam(3,ncdim,*)
      integer ncompv(nv), ncompc(nc)
      logical lvtest

      logical lbound

      fysym = float(iysym)
c     fzsym = float(izsym)

c...  nested pair of loops to calculate the normalwash influence matrix
c     the outer loop runs over all the control points
c     the inner loop runs over all the vortex elements 

      do 200 i = 1, nc
c------- field point location
         x = rc(1,i)
         y = rc(2,i)
         z = rc(3,i)

         u = 0.0
         v = 0.
         w = 0.

         zisgn = z - zsym

         do 100 j = 1, nv
           izsym = izimv(j)
           fzsym = float(izsym)
           if(izsym .ne. 0) then
c---------- no influence if rc and rv are on opposite sides of z=zsym           
            zjsgn = 0.5*(rv1(3,j)+rv2(3,j)) - zsym
            if(zisgn*zjsgn .le. 0.0) go to 100
           endif
           
c--------- set vortex core
           dsyz = sqrt(  (rv2(2,j)-rv1(2,j))**2
     &                 + (rv2(3,j)-rv1(3,j))**2 )
           if(ncompc(i) .eq. ncompv(j)) then
c--------- Original same component vortex core scaling
ccc            rcore = 0.001*chordv(j)
c--------- NEW same component vortex core scaling
            rcore = 0.0001*dsyz
c---------BUG? setting rcore to 0.001*dsyz fails with errors in CLa - WHY?
c         occurs for high AR1000 wing geometry
ccc            rcore = 0.001*dsyz
c---------BUG? setting rcore to 0.0 fails with NaN's - WHY?
ccc            rcore = 0.0
           else
c--------- Original vortex core scaling
ccc            rcore = max( vrcore*chordv(j) , 2.0*vrcore*dsyz )
c--------- NEW independent core scaling on max of *chord and *strip width
            rcore = max( vrcorec*chordv(j) , vrcorew*dsyz )
           endif

           ui = 0.0
           vi = 0.0
           wi = 0.0

           yoff = 2.0* ysym
           zoff = 2.0* zsym
ccc   zoff = 2.0*(zsym + alfa*0.5*(rv1(1,j)+rv2(1,j)) )

c--------- influence of the real vortex
           lbound = .not.(lvtest  .and. i.eq.j)
           call vorvelc(x,y,z,lbound,
     &           rv1(1,j),rv1(2,j),rv1(3,j),
     &           rv2(1,j),rv2(2,j),rv2(3,j),
     &           betm,u,v,w,rcore)

           if(iysym.ne.0) then
c----------- influence of the y-image vortex
             lbound = .true.
c...  for sym/asym matrices check for vortex midpoints of image vortices
             if(iysym.eq.1) then
               xave =        0.5*(rv1(1,j)+rv2(1,j))
               yave = yoff - 0.5*(rv1(2,j)+rv2(2,j))
               zave =        0.5*(rv1(3,j)+rv2(3,j))
               if(x.eq.xave .and. 
     &            y.eq.yave .and.
     &            z.eq.zave       ) lbound = .false.
ccc   if(.not.lbound) write(*,*) 'pos self vortex i,j ',i,j
             endif
             call vorvelc(x,y,z,lbound,
     &             rv2(1,j),yoff-rv2(2,j),rv2(3,j),
     &             rv1(1,j),yoff-rv1(2,j),rv1(3,j),
     &             betm,ui,vi,wi,rcore)
             ui = ui*fysym
             vi = vi*fysym
             wi = wi*fysym 
           endif

           if(izsym.ne.0) then
c----------- influence of the z-image vortex
             lbound = .true.
             call vorvelc(x,y,z,lbound,
     &             rv2(1,j),rv2(2,j),zoff-rv2(3,j),
     &             rv1(1,j),rv1(2,j),zoff-rv1(3,j),
     &             betm,uii,vii,wii,rcore)
             u = u + uii*fzsym
             v = v + vii*fzsym
             w = w + wii*fzsym
c----------- influence of the y,z-image vortex
             if(iysym.ne.0) then
               lbound = .true.
               call vorvelc(x,y,z,lbound,
     &              rv1(1,j),yoff-rv1(2,j),zoff-rv1(3,j),
     &              rv2(1,j),yoff-rv2(2,j),zoff-rv2(3,j),
     &              betm,uii,vii,wii,rcore)

               ui = ui + uii*fysym*fzsym
               vi = vi + vii*fysym*fzsym
               wi = wi + wii*fysym*fzsym
             endif
           endif

           us = u + ui
           vs = v + vi
           ws = w + wi

           vc_gam(1,i,j) = us
           vc_gam(2,i,j) = vs
           vc_gam(3,i,j) = ws

 100     continue ! next vortex j
 200  continue ! next field point i

      return
      end ! vvor



      subroutine vsrd(betm,iysym,ysym,zsym,srcore,
     &                nbody,lfrst,llast,nldim,
     &                rl1,rl2,ncompl,drbl,darl,radl,iziml,
     &                nu,src_u,dbl_u,
     &                nc,rc,ncompc,
     &                wc_u,ncdim)
c--------------------------------------------------------------------
c     Calculates the velocity influence matrix for a collection 
c     of source+doublet lines, at specified field points
c     
c Input
c -----
c       betm     sqrt(1-mach*mach)
c       iysym    plane of symmetry xz 
c                 = 0 no symmetry
c                 = 1 regular symmetry
c                 =-1 free-surface symmetry
c       ysym     y coordinate of symmetry plane
c       iziml(.) second plane of symmetry xy for each body segment
c                 = 0 no second plane
c                 = 1 regular symmetry
c                 =-1 free-surface symmetry
c       zsym     z coordinate of symmetry plane
c
c       srcore   source-line core radius / body radius ratio
c
c       nbody      number of bodies
c       lfrst(b)   index of first segment in body b
c       llast(b)   index of last  segment in body b
c       nldim      size of src_u, dbl_u matrices
c       rl1(3,b)   1st source-line node
c       rl2(3,b)   2nd source-line node
c       ncompl(b)  index of component containing source line
c       drbl(b)    body radius change over segment
c       darl(b)    area change over segment
c       radl(b)    body radius at line segment
c
c       nu         number of apparent-freestream components
c       src_u(u)   source  strength per unit freestream component
c       dbl_u(3,u) doublet strength per unit freestream component
c
c       nc         number of field points
c       rc(3,c)    field point node where velocity is evaluated
c       ncompc(c)  index of component containing field point
c
c       ncdim      size of wc matrix
c          
c Output
c ------
c       wc_u(3,c,u)  velocity for unit freestream, rotation components
c
c--------------------------------------------------------------------
      integer lfrst(*), llast(*)
      integer ncompl(*), ncompc(*)
      real rl1(3,*), rl2(3,*), drbl(*), darl(*), radl(*)
      integer iziml(*)
      real src_u(nldim,*), dbl_u(3,nldim,*),
     &     rc(3,ncdim),
     &     wc_u(3,ncdim,*)

      real vsrc(3), vdbl(3,3)
      data pi / 3.14159265358979323846264338327950280 /

      fysym = float(iysym)
c     fzsym = float(izsym)

      yoff = 2.0* ysym
      zoff = 2.0* zsym
ccc   zoff = 2.0*(zsym + alfa*0.5*(rv1(1,j)+rv2(1,j)) )

      do i = 1, nc
        do iu = 1, nu
          wc_u(1,i,iu) = 0.
          wc_u(2,i,iu) = 0.
          wc_u(3,i,iu) = 0.
        enddo
      enddo

      do 10 ibody = 1, nbody
c------ go over body segments, and l=llast(ibody)+1  wake segment
        do 105 l = lfrst(ibody), llast(ibody)+1
          izsym = iziml(l)
          fzsym = float(izsym)

          zlsgn = 0.5*(rl1(3,l)+rl2(3,l)) - zsym

          do 1005 i = 1, nc
            if(ncompc(i) .eq. ncompl(l)) then
c----------- sources and doublets of a body do not contribute to velocity
c-            on that body (the doublet near-fields are too singular)
             go to 1005
            endif

            zisgn = rc(3,i) - zsym
            if(izsym .ne. 0) then
c----------- no influence if rc and rv are on opposite sides of z=zsym           
             if(zisgn*zlsgn .le. 0.0) go to 1005
            endif

            
c---------- original core scaling with radius
ccc            rcore = srcore*radl(l)
c---------- NEW scale core on length of segment (HHY 3/23)
            rlavg = sqrt( (rl2(1,l2)-rl1(1,l1))**2 +
     &                    (rl2(2,l2)-rl1(2,l1))**2 +
     &                    (rl2(3,l2)-rl1(3,l1))**2 )
            rcore = srcore*rlavg
c------------------------------------------------------------
c---------- influence of real segment
            call srdvelc(rc(1,i) ,rc(2,i) ,rc(3,i) ,
     &                   rl1(1,l),rl1(2,l),rl1(3,l),
     &                   rl2(1,l),rl2(2,l),rl1(3,l),
     &                   betm,rcore,
     &                   vsrc,vdbl  )
            do iu = 1, nu
              do k = 1, 3
                wc_u(k,i,iu) = wc_u(k,i,iu)
     &               + vsrc(k)*src_u(l,iu)
     &               + vdbl(k,1)*dbl_u(1,l,iu)
     &               + vdbl(k,2)*dbl_u(2,l,iu)
     &               + vdbl(k,3)*dbl_u(3,l,iu)
              enddo
            enddo

c------------------------------------------------------------
            if (iysym.ne.0) then
c----------- influence of y-image
             call srdvelc(rc(1,i) ,     rc(2,i) ,rc(3,i) ,
     &                    rl1(1,l),yoff-rl1(2,l),rl1(3,l),
     &                    rl2(1,l),yoff-rl2(2,l),rl2(3,l),
     &                    betm,rcore,
     &                    vsrc,vdbl  )
             do iu = 1, nu
               do k = 1, 3
                 wc_u(k,i,iu) = wc_u(k,i,iu)
     &                  + ( vsrc(k)*src_u(l,iu)
     &                    + vdbl(k,1)*dbl_u(1,l,iu)
     &                    - vdbl(k,2)*dbl_u(2,l,iu)
     &                    + vdbl(k,3)*dbl_u(3,l,iu) )*fysym
               enddo
             enddo
            endif

c------------------------------------------------------------
            if (izsym.ne.0) then
c----------- influence of z-image
             call srdvelc(rc(1,i) ,rc(2,i) ,     rc(3,i) ,
     &                    rl1(1,l),rl1(2,l),zoff-rl1(3,l),
     &                    rl2(1,l),rl2(2,l),zoff-rl2(3,l),
     &                    betm,rcore,
     &                    vsrc,vdbl  )
             do iu = 1, nu
               do k = 1, 3
                 wc_u(k,i,iu) = wc_u(k,i,iu)
     &                  + ( vsrc(k)*src_u(l,iu)
     &                    + vdbl(k,1)*dbl_u(1,l,iu)
     &                    + vdbl(k,2)*dbl_u(2,l,iu)
     &                    - vdbl(k,3)*dbl_u(3,l,iu) )*fzsym
               enddo
             enddo

c--------------------------------------
             if (iysym.ne.0) then
c------------ influence of y,z-image
              call srdvelc(rc(1,i) ,     rc(2,i) ,     rc(3,i) ,
     &                     rl1(1,l),yoff-rl1(2,l),zoff-rl1(3,l),
     &                     rl2(1,l),yoff-rl2(2,l),zoff-rl2(3,l),
     &                     betm,rcore,
     &                     vsrc,vdbl  )
              do iu = 1, nu
                do k = 1, 3
                  wc_u(k,i,iu) = wc_u(k,i,iu)
     &                   + ( vsrc(k)*src_u(l,iu)
     &                     + vdbl(k,1)*dbl_u(1,l,iu)
     &                     - vdbl(k,2)*dbl_u(2,l,iu)
     &                     - vdbl(k,3)*dbl_u(3,l,iu) )*fysym*fzsym
                enddo
              enddo
             endif
            endif

 1005     continue ! next field point i
 105    continue ! next body segment l
 10   continue ! next ibody

      return
      end ! vsrd



      subroutine srdu(betm,xyzref,
     &               nbody,lfrst,llast,nldim,
     &               ibcent, iysym,
     &               rl1,rl2,rl,drbl,darl,radl,
     &               src_u,dbl_u )
c----------------------------------------------------------
c     Sets source+doublet line segment AICs w.r.t.
c     freestream and rotation components
c
c Input
c -----
c       betm     sqrt(1-mach*mach)
c
c       nbody      number of body segments
c       lfrst(b)   index of first segment in body b
c       llast(b)   index of last  segment in body b
c       nldim      size of src_u, dbl_u matrices
c       rl1(3,b)   source-line 1st node
c       rl2(3,b)   source-line 2nd node
c       rl(3,b)    source-line control point
c       drbl(b)    body radius change over segment
c       darl(b)    body area change over segment
c       radl(b)    body radius at control point
c
c Output
c ------
c       src_u(u)   source  strength per unit freestream component
c       dbl_u(3,u) doublet strength per unit freestream component
c----------------------------------------------------------
      real xyzref(3)
      integer lfrst(*), llast(*), ibcent(*)
      real rl1(3,*), rl2(3,*), rl(3,*), drbl(*), darl(*), radl(*)
      real src_u(nldim,*), dbl_u(3,nldim,*)

      real drl(3), esl(3), un(3)
      real wrot(3), urel(3), rlrot(3)

      data pi / 3.14159265358979323846264338327950280 /

      do 10 ibody = 1, nbody
        if(iysym .eq. 1 .and.
     &     ibcent(ibody) .eq. 1) then
c------- body y-image will be added on, so use only half the area
         sdfac = 0.5
        else
c------- no y-image, so use entire area
         sdfac = 1.0
        endif

        do 105 l = lfrst(ibody), llast(ibody)
          drl(1) = (rl2(1,l) - rl1(1,l))/betm
          drl(2) =  rl2(2,l) - rl1(2,l)
          drl(3) =  rl2(3,l) - rl1(3,l)
          drlmag = sqrt(drl(1)**2 + drl(2)**2 + drl(3)**2)
          if(drlmag.eq.0.0) then
           drlmi = 0.0
          else
           drlmi = 1.0/drlmag
          endif

c-------- unit vector along line segment
          esl(1) = drl(1) * drlmi
          esl(2) = drl(2) * drlmi
          esl(3) = drl(3) * drlmi

          drdl =         drbl(l)*drlmi
          adel = sdfac * darl(l)
          area = sdfac * pi*radl(l)**2

          rlrot(1) = rl(1,l) - xyzref(1)
          rlrot(2) = rl(2,l) - xyzref(2)
          rlrot(3) = rl(3,l) - xyzref(3)

c-------- go over freestream velocity and rotation components
          do iu = 1, 6
            urel(1) = 0.
            urel(2) = 0.
            urel(3) = 0.
            wrot(1) = 0.
            wrot(2) = 0.
            wrot(3) = 0.

            if  (iu.le.3) then
             urel(iu) = 1.0
            else
             wrot(iu-3) = 1.0
             call cross(rlrot,wrot,urel)
            endif
            urel(1) = urel(1)/betm

c---------- u.es
            us = urel(1)*esl(1) + urel(2)*esl(2) + urel(3)*esl(3)

c---------- velocity projected on normal plane = u - (u.es) es
            un(1) = urel(1) - us*esl(1)
            un(2) = urel(2) - us*esl(2)
            un(3) = urel(3) - us*esl(3)

c---------- source-line and doublet-line strengths of this segment
            asrc = adel
            adbl = area / (1.0 + 16.0*drdl**2)
            src_u(l,iu)   = asrc*us*drlmi
            dbl_u(1,l,iu) = adbl*un(1)*2.0
            dbl_u(2,l,iu) = adbl*un(2)*2.0
            dbl_u(3,l,iu) = adbl*un(3)*2.0
          enddo
 105    continue ! next body segment l

c------ base area terms
        l = llast(ibody)

c------ set doublet strength at end of body
        if(lfrst(ibody) .eq. llast(ibody)) then
c------- only a single body segment... extrapolate from nose,midpoint to end
         fo = 2.0
         do k = 1, 3
ccc        dbl(k,l+1) = fo*dbl(k,l)
           do iu = 1, 6
             dbl_u(k,l+1,iu) = fo*dbl_u(k,l,iu)
           enddo
         enddo
        else
c------- extrapolate from l,l-1 points to end
         fo = 1.5
         fm = -0.5
         do k = 1, 3
ccc        dbl(k,l+1) = fo*dbl(k,l) + fm*dbl(k,l-1)
           do iu = 1, 6
             dbl_u(k,l+1,iu) = fo*dbl_u(k,l,iu) + fm*dbl_u(k,l-1,iu)
           enddo
         enddo
        endif

 10   continue ! next ibody

      return
      end ! srdu


 
      subroutine cross (u,v,w)
      real u(3), v(3), w(3)
      w(1) = u(2)*v(3) - u(3)*v(2)
      w(2) = u(3)*v(1) - u(1)*v(3)
      w(3) = u(1)*v(2) - u(2)*v(1)
      return
      end      


      function dot (u,v)
      real u(3),v(3)
      dot = u(1)*v(1) + u(2)*v(2) + u(3)*v(3)
      return
      end





      subroutine vorvel(x,y,z,lbound,x1,y1,z1,x2,y2,z2,beta,
     &                   u,v,w )
c----------------------------------------------------------
c     Computes velocity at x,y,z of a horseshoe vortex
c     with corners at r1, r2, with trailing legs along x.
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
      if (lbound .and.  amag*bmag .ne. 0.0) then
        axb(1) = a(2)*b(3) - a(3)*b(2)
        axb(2) = a(3)*b(1) - a(1)*b(3)
        axb(3) = a(1)*b(2) - a(2)*b(1)

        adb = a(1)*b(1) + a(2)*b(2) + a(3)*b(3)

        den = amag*bmag + adb
        if(den.ne.0.0) then
         t = (1.0/amag + 1.0/bmag) / den
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

        t =  (1.0 + adi/amag) / rsq

        v = v + a(3)*t
        w = w - a(2)*t
      endif

c---- trailing leg attached to b
      if (bmag .ne. 0.0) then
        bxisq = b(3)**2 + b(2)**2

        bdi = b(1)
        rsq = bxisq

        t = -(1.0 + bdi/bmag) / rsq

        v = v + b(3)*t
        w = w - b(2)*t
      endif

      u = u*pi4inv / beta
      v = v*pi4inv 
      w = w*pi4inv 

      return
      end ! vorvel


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


      subroutine srdvelc(x,y,z, x1,y1,z1, x2,y2,z2,
     &                   beta,rcore,
     &                   uvws,uvwd  )
c-------------------------------------------------------------------
c     Same as srdvel, but with finite core radius
c-------------------------------------------------------------------
      real uvws(3), uvwd(3,3)

      real r1(3), r2(3)
      real rxr(3)

c---- 1 / (4 pi)
      data pi4inv / 0.079577471545947667884441881686257184 /

      r1(1) = (x-x1)/beta
      r1(2) =  y-y1
      r1(3) =  z-z1

      r2(1) = (x-x2)/beta
      r2(2) =  y-y2
      r2(3) =  z-z2

      ds = sqrt( (r2(1)-r1(1))**2
     &         + (r2(2)-r1(2))**2
     &         + (r2(3)-r1(3))**2 )

      rcsq = rcore**2

      r1sq = r1(1)**2 + r1(2)**2 + r1(3)**2
      r2sq = r2(1)**2 + r2(2)**2 + r2(3)**2

      r1eps = sqrt(r1sq + rcsq)
      r2eps = sqrt(r2sq + rcsq)

      if(r1eps .eq. 0.0) then
       r1i = 0.
      else
       r1i = 1.0/r1eps
      endif

      if(r2eps .eq. 0.0) then
       r2i = 0.
      else
       r2i = 1.0/r2eps
      endif

      rdr = r1(1)*r2(1) + r1(2)*r2(2) + r1(3)*r2(3)
      rxr(1) = r1(2)*r2(3) - r1(3)*r2(2)
      rxr(2) = r1(3)*r2(1) - r1(1)*r2(3)
      rxr(3) = r1(1)*r2(2) - r1(2)*r2(1)

      rxrsq = rxr(1)**2 + rxr(2)**2 + rxr(3)**2

      all = r1sq + r2sq - 2.0*rdr

      den = rxrsq + rcsq*all
      if(den .eq. 0.0) then
       deni = 0.
      else
       deni = 1.0/den
      endif

      ai1 = (r2eps - (rdr+rcsq)*r1i)*deni
      ai2 = (r1eps - (rdr+rcsq)*r2i)*deni

c---- set velocity components for unit source and doublet
      do k = 1, 3
        uvws(k) = r1(k)*ai1 + r2(k)*ai2

        rr1 =  (r1(k)+r2(k))    *r1i
     &        - r1(k)*(rdr+rcsq)*r1i**3
     &        - r2(k)           *r2i

        rr2 =  (r1(k)+r2(k))    *r2i
     &        - r2(k)*(rdr+rcsq)*r2i**3
     &        - r1(k)           *r1i

        rrt = 2.0*r1(k)*(r2sq  - rdr)
     &      + 2.0*r2(k)*(r1sq  - rdr)  

        aj1 = -(rr1 + ai1*rrt)*deni

        aj2 = -(rr2 + ai2*rrt)*deni

        do j = 1, 3
          uvwd(k,j) = aj1*r1(j)
     &              + aj2*r2(j)
        enddo

        uvwd(k,k) = uvwd(k,k) + ai1 + ai2
      enddo

      uvws(1) = ds*uvws(1)*pi4inv/beta
      uvws(2) = ds*uvws(2)*pi4inv
      uvws(3) = ds*uvws(3)*pi4inv
      do l = 1, 3
        uvwd(1,l) = ds*uvwd(1,l)*pi4inv/beta
        uvwd(2,l) = ds*uvwd(2,l)*pi4inv
        uvwd(3,l) = ds*uvwd(3,l)*pi4inv
      enddo

      return
      end ! srdvelc
