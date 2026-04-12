C***********************************************************************
C    Module:  jbdforc.f
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

      subroutine bdforc
c-----------------------------------------------------------------------
c     Calculate the forces on the bodies
c
c  Inputs:   
c     alfa     angle of attack (for stability-axis definition)
c     vbar(.)  freestream velocity components
c     wbar(.)  roll,pitch,yaw  rates
c     mach     mach number
c     nbody    number of bodies
c     lfrst(.) first node in each body
c     llast(.) last  node in each body
c     rl(..)   line segment coordinates
c     radl(.)  body radii
c          
c  Outputs:  
c     dcpb(..) body element loadings
c     cfbi(..)  body force/q  (not normalized by Sref)
c     cmbi(..)  body moment/q (not normalized by Sref*cref)
c     cfbi_u(..)  d cfbi / d (Vbar,Wbar)
c     cfbi_d(..)  d cfbi / d control
c       .
c       .     
c-----------------------------------------------------------------------
      include 'jvl.inc'

      real vinf, vinf_u(6)
      real rrot(3)
      real wxr(3), wxr_w(3,3)

      real vtot(3),
     &     vtot_u(3,numax),
     &     vtot_d(3,ndmax),
     &     vtot_j(3,njmax),
     &     vtot_g(3,ndmax)

      real r(3), dr(3), esl(3)
      real vs,
     &     vs_u(numax),
     &     vs_d(ndmax),
     &     vs_j(njmax),
     &     vs_g(ngmax)
      real vn(3),
     &     vn_u(3,numax),
     &     vn_d(3,ndmax),
     &     vn_j(3,njmax),
     &     vn_g(3,ngmax)

      real fsym(3), msym(3)

      integer icrs(3), jcrs(3)
      data icrs / 2, 3, 1 / , jcrs / 3, 1, 2 /

      vinf = sqrt(vbar(1)**2 + vbar(2)**2 + vbar(3)**2)
      vinf_u(1) = vbar(1)/vinf
      vinf_u(2) = vbar(2)/vinf
      vinf_u(3) = vbar(3)/vinf
      vinf_u(4) = 0.
      vinf_u(5) = 0.
      vinf_u(6) = 0.

      betm = sqrt(1.0 - mach**2)

c---- forces and moments for each body
      do 100 ib = 1, nbody
        dmax = 0.0
c------ clear for accumulation over body elements
        do k = 1, 3
          fbi(k,ib) = 0.
          mbi(k,ib) = 0.
          fbv(k,ib) = 0.
          mbv(k,ib) = 0.
          do n = 1, numax
            fbi_u(k,ib,n) = 0.
            mbi_u(k,ib,n) = 0.
            fbv_u(k,ib,n) = 0.
            mbv_u(k,ib,n) = 0.
          enddo
          do n = 1, ncontrol
            fbi_d(k,ib,n) = 0.
            mbi_d(k,ib,n) = 0.
            fbv_d(k,ib,n) = 0.
            mbv_d(k,ib,n) = 0.
          enddo
          do n = 1, nvarjet
            fbi_j(k,ib,n) = 0.
            mbi_j(k,ib,n) = 0.
            fbv_j(k,ib,n) = 0.
            mbv_j(k,ib,n) = 0.
          enddo
          do n = 1, ndesign
            fbi_g(k,ib,n) = 0.
            mbi_g(k,ib,n) = 0.
            fbv_g(k,ib,n) = 0.
            mbv_g(k,ib,n) = 0.
          enddo
        enddo ! next k

        if(iysym .eq. 1 .and.
     &     ibcent(ib) .eq. 1) then
c------- image will be added on, so set only half the strength
         sdfac = 0.5
        else
         sdfac = 1.0
        endif
ccc        write(14,*) ib,sdfac,lfrst(ib),llast(ib)
        
c------ go over body segments
        do 110 l = lfrst(ib), llast(ib)
          dr(1) = rl2(1,l) - rl1(1,l)
          dr(2) = rl2(2,l) - rl1(2,l)
          dr(3) = rl2(3,l) - rl1(3,l)
          drmag = sqrt(dr(1)**2 + dr(2)**2 + dr(3)**2)
          if(drmag.eq.0.0) then
           drmi = 0.0
          else
           drmi = 1.0/drmag
          endif

          drdl =         drbl(l)*drmi
          adel = sdfac * darl(l)
          area = sdfac * pi*radl(l)**2

          dia = 2.0*radl(l)
          dmax = max(dmax,dia)
          if(dia.le.0.0) then
           dinv = 0.0
          else
           dinv = 1.0/dia
          endif

c-------- unit vector along line segment
          esl(1) = dr(1) * drmi
          esl(2) = dr(2) * drmi
          esl(3) = dr(3) * drmi

c-------- rotation arm relative to xyzref
          rrot(1) = rl(1,l) - xyzref(1)
          rrot(2) = rl(2,l) - xyzref(2)
          rrot(3) = rl(3,l) - xyzref(3)

c-------- Wbar x r
          do k = 1, 3
            ic = icrs(k)
            jc = jcrs(k)
            wxr(k) = wbar(ic)*rrot(jc)
     &             - wbar(jc)*rrot(ic)
            wxr_w(k,ic) =  rrot(jc)
            wxr_w(k,jc) = -rrot(ic)
            wxr_w(k,k) = 0.
          enddo

c-------- total relative velocity at r
          do k = 1, 3
            vtot(k) = vbar(k) - wxr(k) + vl(k,l) + wl(k,l)
            do n = 1, numax
              vtot_u(k,n) = wl_u(k,l,n)
              do i = 1, nvor
                vtot_u(k,n) = vtot_u(k,n) + vl_gam(k,l,i)*gam_u(i,n)
              enddo
            enddo

            do n = 1, ncontrol
              vtot_d(k,n) = 0.
              do i = 1, nvor
                vtot_d(k,n) = vtot_d(k,n) + vl_gam(k,l,i)*gam_d(i,n)
              enddo
            enddo

            do n = 1, nvarjet
              vtot_j(k,n) = 0.
              do i = 1, nvor
                vtot_j(k,n) = vtot_j(k,n) + vl_gam(k,l,i)*gam_j(i,n)
              enddo
            enddo

            do n = 1, ndesign
              vtot_g(k,n) = 0.
              do i = 1, nvor
                vtot_g(k,n) = vtot_g(k,n) + vl_gam(k,l,i)*gam_g(i,n)
              enddo
            enddo

c           vtot(k)     = vbar(k) - wxr(k)     +  vl(k,l) + wl(k,l)
            vtot_u(k,k) = vtot_u(k,k) + 1.0
            vtot_u(k,4) = vtot_u(k,4) - wxr_w(k,1)
            vtot_u(k,5) = vtot_u(k,5) - wxr_w(k,2)
            vtot_u(k,6) = vtot_u(k,6) - wxr_w(k,3)
          enddo ! next k
ccc          write(15,98) ib,l,vtot,wxr,(vl(k,l),k=1,3),(wl(k,l),k=1,3)
 98       format(i4,i4,12g12.4)

c-------- V.es
          vs = vtot(1)*esl(1) + vtot(2)*esl(2) + vtot(3)*esl(3)
          do n = 1, numax
            vs_u(n) = vtot_u(1,n)*esl(1)
     &              + vtot_u(2,n)*esl(2)
     &              + vtot_u(3,n)*esl(3)
          enddo
          do n = 1, ncontrol
            vs_d(n) = vtot_d(1,n)*esl(1)
     &              + vtot_d(2,n)*esl(2)
     &              + vtot_d(3,n)*esl(3)
          enddo
          do n = 1, nvarjet
            vs_j(n) = vtot_j(1,n)*esl(1)
     &              + vtot_j(2,n)*esl(2)
     &              + vtot_j(3,n)*esl(3)
          enddo
          do n = 1, ndesign
            vs_g(n) = vtot_g(1,n)*esl(1)
     &              + vtot_g(2,n)*esl(2)
     &              + vtot_g(3,n)*esl(3)
          enddo

          do k = 1, 3
c---------- velocity projected on normal plane = V - (V.es) es
            vn(k) = vtot(k) - vs*esl(k)
            do n = 1, numax
              vn_u(k,n) = vtot_u(k,n) - vs_u(n)*esl(k)
            enddo
            do n = 1, ncontrol
              vn_d(k,n) = vtot_d(k,n) - vs_d(n)*esl(k)
            enddo
            do n = 1, nvarjet
              vn_j(k,n) = vtot_j(k,n) - vs_j(n)*esl(k)
            enddo
            do n = 1, ndesign
              vn_g(k,n) = vtot_g(k,n) - vs_g(n)*esl(k)
            enddo
          enddo

c-------- segment source and doublet strengths from dA/dl and area
          asrc = adel
          adbl = area 
c-----segment area correction (??) based on radius slope dR/dl
c          adbl = adbl / (1.0 + 16.0*drdl**2)
          src(l)   = asrc*vs*drmi
          dbl(1,l) = adbl*vn(1)*2.0
          dbl(2,l) = adbl*vn(2)*2.0
          dbl(3,l) = adbl*vn(3)*2.0

ccc          write(14,99) ib,l,vs,vn,vtot,dbl(1,l),dbl(2,l),dbl(3,l)
 99       format(i4,i4,10g12.4)
c          if(l.eq.lfrst(ib)) then
c           vint0 = 0.
c           vint1 = 0.
c          endif
c          vint0 = vint0 + drmag*area
c          vint1 = vint1 + drmag*adbl
c          if(l.eq.llast(ib)) then
c            write(*,*) vint0, vint1
c          endif

          do n = 1, numax
            src_u(l,n)   = asrc*vs_u(n)*drmi
            dbl_u(1,l,n) = adbl*vn_u(1,n)*2.0
            dbl_u(2,l,n) = adbl*vn_u(2,n)*2.0
            dbl_u(3,l,n) = adbl*vn_u(3,n)*2.0
          enddo
          do n = 1, ncontrol
            src_d(l,n)   = asrc*vs_d(n)*drmi
            dbl_d(1,l,n) = adbl*vn_d(1,n)*2.0
            dbl_d(2,l,n) = adbl*vn_d(2,n)*2.0
            dbl_d(3,l,n) = adbl*vn_d(3,n)*2.0
          enddo
          do n = 1, nvarjet
            src_j(l,n)   = asrc*vs_j(n)*drmi
            dbl_j(1,l,n) = adbl*vn_j(1,n)*2.0
            dbl_j(2,l,n) = adbl*vn_j(2,n)*2.0
            dbl_j(3,l,n) = adbl*vn_j(3,n)*2.0
          enddo
          do n = 1, ndesign
            src_g(l,n)   = asrc*vs_g(n)*drmi
            dbl_g(1,l,n) = adbl*vn_g(1,n)*2.0
            dbl_g(2,l,n) = adbl*vn_g(2,n)*2.0
            dbl_g(3,l,n) = adbl*vn_g(3,n)*2.0
          enddo

          rho = rhol(l)

          do k = 1, 3
            ic = icrs(k)
            jc = jcrs(k)
            mbi(k,ib) = mbi(k,ib)
     &                 - 0.5*rho*vinf*( dr(ic)*dbl(jc,l)
     &                                - dr(jc)*dbl(ic,l) )
            do n = 1, numax
              mbi_u(k,ib,n) = mbi_u(k,ib,n)
     &                 - 0.5*rho*vinf*( dr(ic)*dbl_u(jc,l,n)
     &                                - dr(jc)*dbl_u(ic,l,n))
     &                 - 0.5*rho*vinf_u(n)*( dr(ic)*dbl(jc,l)
     &                                     - dr(jc)*dbl(ic,l) )
            enddo
            do n = 1, ncontrol
              mbi_d(k,ib,n) = mbi_d(k,ib,n)
     &                 - 0.5*rho*vinf*( dr(ic)*dbl_d(jc,l,n)
     &                                - dr(jc)*dbl_d(ic,l,n))
            enddo
            do n = 1, nvarjet
              mbi_j(k,ib,n) = mbi_j(k,ib,n)
     &                 - 0.5*rho*vinf*( dr(ic)*dbl_j(jc,l,n)
     &                                - dr(jc)*dbl_j(ic,l,n))
            enddo
            do n = 1, ndesign
              mbi_g(k,ib,n) = mbi_g(k,ib,n)
     &                 - 0.5*rho*vinf*( dr(ic)*dbl_g(jc,l,n)
     &                                - dr(jc)*dbl_g(ic,l,n))
            enddo
          enddo
 110    continue ! next ilseg

c------ base area terms...if needed
        l = llast(ib)
        if(rad2(l).gt.0.0001*dmax) then

c------ set doublet strength at end of body
         if(lfrst(ib) .eq. llast(ib)) then
c------- only a single body segment... extrapolate from nose,midpoint to end
          fo = 2.0
          do k = 1, 3
           dbl(k,l+1) = fo*dbl(k,l)
ccc           write(21,*) 'wake1 ib,l,k,dbl(k,l+1),dbl(k,l) ',
ccc     &                 ib,l,k,dbl(k,l),dbl(k,l+1)
           do n = 1, numax
             dbl_u(k,l+1,n) = fo*dbl_u(k,l,n)
           enddo
           do n = 1, ncontrol
             dbl_d(k,l+1,n) = fo*dbl_d(k,l,n)
           enddo
           do n = 1, nvarjet
             dbl_j(k,l+1,n) = fo*dbl_j(k,l,n)
           enddo
           do n = 1, ndesign
             dbl_g(k,l+1,n) = fo*dbl_g(k,l,n)
           enddo
          enddo

         else
c------- extrapolate from l,l-1 points to end
          fo = 1.5
          fm = -0.5
          do k = 1, 3
           dbl(k,l+1) = fo*dbl(k,l) + fm*dbl(k,l-1)
ccc           write(21,*) 'wake2 ib,l,k,dbl(k,l-1),dbl(k,l),dbl(k,l+1) ',
ccc     &                 ib,l,k,dbl(k,l-1),dbl(k,l),dbl(k,l+1)
           do n = 1, numax
             dbl_u(k,l+1,n) = fo*dbl_u(k,l,n) + fm*dbl_u(k,l-1,n)
           enddo
           do n = 1, ncontrol
             dbl_d(k,l+1,n) = fo*dbl_d(k,l,n) + fm*dbl_d(k,l-1,n)
           enddo
           do n = 1, nvarjet
             dbl_j(k,l+1,n) = fo*dbl_j(k,l,n) + fm*dbl_j(k,l-1,n)
           enddo
           do n = 1, ndesign
             dbl_g(k,l+1,n) = fo*dbl_g(k,l,n) + fm*dbl_g(k,l-1,n)
           enddo
          enddo
         endif

        endif

c------ moment arm
        r(1) = rl2(1,l) - xyzmom(1)
        r(2) = rl2(2,l) - xyzmom(2)
        r(3) = rl2(3,l) - xyzmom(3)

        rho = rhol(l)

        do k = 1, 3
          ic = icrs(k)
          jc = jcrs(k)
          fbi(k,ib) = 0.5*rho*vinf*dbl(k,l+1)
          mbi(k,ib) = mbi(k,ib)
     &              + 0.5*rho*vinf*( r(ic)*dbl(jc,l+1)
     &                             - r(jc)*dbl(ic,l+1) )
          do n = 1, numax
            fbi_u(k,ib,n) =
     &      fbi_u(k,ib,n) + 0.5*rho*vinf_u(n)*dbl(k,l+1)
     &                    + 0.5*rho*vinf     *dbl_u(k,l+1,n)
            mbi_u(k,ib,n) =
     &      mbi_u(k,ib,n) + 0.5*rho*vinf_u(n)*( r(ic)*dbl(jc,l+1)
     &                                        - r(jc)*dbl(ic,l+1)     )
     &                    + 0.5*rho*vinf     *( r(ic)*dbl_u(jc,l+1,n)
     &                                        - r(jc)*dbl_u(ic,l+1,n) )
          enddo
          do n = 1, ncontrol
            fbi_d(k,ib,n) =
     &      fbi_d(k,ib,n) + 0.5*rho*vinf*dbl_d(k,l+1,n)
            mbi_d(k,ib,n) =
     &      mbi_d(k,ib,n) + 0.5*rho*vinf*( r(ic)*dbl_d(jc,l+1,n)
     &                                   - r(jc)*dbl_d(ic,l+1,n) )
          enddo
          do n = 1, nvarjet
            fbi_j(k,ib,n) =
     &      fbi_j(k,ib,n) + 0.5*rho*vinf*dbl_j(k,l+1,n)
            mbi_j(k,ib,n) =
     &      mbi_j(k,ib,n) + 0.5*rho*vinf*( r(ic)*dbl_j(jc,l+1,n)
     &                                   - r(jc)*dbl_j(ic,l+1,n) )
          enddo
          do n = 1, ndesign
            fbi_g(k,ib,n) =
     &      fbi_g(k,ib,n) + 0.5*rho*vinf*dbl_g(k,l+1,n)
            mbi_g(k,ib,n) =
     &      mbi_g(k,ib,n) + 0.5*rho*vinf*( r(ic)*dbl_g(jc,l+1,n)
     &                                   - r(jc)*dbl_g(ic,l+1,n) )
          enddo
        enddo

c------ set average delta(cp) on area of each body segment (for plotting)
        do l = lfrst(ib), llast(ib)
          dr(1) = rl2(1,l) - rl1(1,l)
          dr(2) = rl2(2,l) - rl1(2,l)
          dr(3) = rl2(3,l) - rl1(3,l)
          drmag = sqrt(dr(1)**2 + dr(2)**2 + dr(3)**2)
          if(drmag.eq.0.0) then
           drmi = 0.0
          else
           drmi = 1.0/drmag
          endif

          dia = 2.0*radl(l)
          if(dia.le.0.0) then
           dinv = 0.0
          else
           dinv = 1.0/dia
          endif

          do k = 1, 3
            if    (l .eq. lfrst(ib)) then
              ddbl = 4.0*dbl(k,l)*drmi          ! assumes body LE zero area
            elseif(l .eq. llast(ib)) then
              if(rad2(l).lt.0.0001*dmax) then   ! closed body at TE?
               ddbl = -4.0*dbl(k,l)*drmi         ! assumes body TE zero area
              else  
               ddbl = (dbl(k,l) - dbl(k,l-1))*drmi
              endif
            else
              ddbl = 0.5*(dbl(k,l+1) - dbl(k,l-1))*drmi
            endif
            dcpb(k,l) = 2.0*vinf*ddbl*dinv
          enddo
        enddo
 100  continue ! next ib

      if(lbforce) then
c---- add on body forces to total forces

      do 300 ib = 1, nbody
        if(.not.lbload(ib)) go to 300

        do k = 1, 3
          fti(k) = fti(k) + fbi(k,ib)
          mti(k) = mti(k) + mbi(k,ib)
          ftv(k) = ftv(k) + fbv(k,ib)
          mtv(k) = mtv(k) + mbv(k,ib)
          do n = 1, numax
            fti_u(k,n) = fti_u(k,n) + fbi_u(k,ib,n)
            mti_u(k,n) = mti_u(k,n) + mbi_u(k,ib,n)
            ftv_u(k,n) = ftv_u(k,n) + fbv_u(k,ib,n)
            mtv_u(k,n) = mtv_u(k,n) + mbv_u(k,ib,n)
          enddo
          do n = 1, ncontrol
            fti_d(k,n) = fti_d(k,n) + fbi_d(k,ib,n)
            mti_d(k,n) = mti_d(k,n) + mbi_d(k,ib,n)
            ftv_d(k,n) = ftv_d(k,n) + fbv_d(k,ib,n)
            mtv_d(k,n) = mtv_d(k,n) + mbv_d(k,ib,n)
          enddo
          do n = 1, nvarjet
            fti_j(k,n) = fti_j(k,n) + fbi_j(k,ib,n)
            mti_j(k,n) = mti_j(k,n) + mbi_j(k,ib,n)
            ftv_j(k,n) = ftv_j(k,n) + fbv_j(k,ib,n)
            mtv_j(k,n) = mtv_j(k,n) + mbv_j(k,ib,n)
          enddo
          do n = 1, ndesign
            fti_g(k,n) = fti_g(k,n) + fbi_g(k,ib,n)
            mti_g(k,n) = mti_g(k,n) + mbi_g(k,ib,n)
            ftv_g(k,n) = ftv_g(k,n) + fbv_g(k,ib,n)
            mtv_g(k,n) = mtv_g(k,n) + mbv_g(k,ib,n)
          enddo
        enddo
 300  continue ! next ib

      endif

      return
      end ! bdforc

