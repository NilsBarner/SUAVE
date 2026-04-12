C***********************************************************************
C    Module:  jaero.f
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

      subroutine aero
c----------------------------------------------------------------------
c     Computes aero forces and moments for current state variables
c     also computes force and moment sensitivities to state variables
c
c   Inputs:  gam(.)
c            vbar(.)
c            wbar(.)
c            delcon(.)
c            deljet(.)
c
c   Outputs: All forces and moments in /forc_r/ common block
c----------------------------------------------------------------------
      include 'jvl.inc'

      real fsym(3), msym(3)
      real r(3), vmag_u(3)

      real sa_u(3), ca_u(3)

      integer icrs(3), jcrs(3)
      data icrs / 2, 3, 1 / , jcrs / 3, 1, 2 /

c---------------------------------------------------------
c---- clear total forces and moments for accumulation
      do k = 1, 3
        fti(k) = 0.
        mti(k) = 0.
        ftj(k) = 0.
        mtj(k) = 0.
        ftp(k) = 0.
        mtp(k) = 0.
        ftv(k) = 0.
        mtv(k) = 0.
        do n = 1, numax
          fti_u(k,n) = 0.
          mti_u(k,n) = 0.
          ftj_u(k,n) = 0.
          mtj_u(k,n) = 0.
          ftp_u(k,n) = 0.
          mtp_u(k,n) = 0.
          ftv_u(k,n) = 0.
          mtv_u(k,n) = 0.
        enddo
        do n = 1, ncontrol
          fti_d(k,n) = 0.
          mti_d(k,n) = 0.
          ftj_d(k,n) = 0.
          mtj_d(k,n) = 0.
          ftp_d(k,n) = 0.
          mtp_d(k,n) = 0.
          ftv_d(k,n) = 0.
          mtv_d(k,n) = 0.
        enddo
        do n = 1, nvarjet
          fti_j(k,n) = 0.
          mti_j(k,n) = 0.
          ftj_j(k,n) = 0.
          mtj_j(k,n) = 0.
          ftp_j(k,n) = 0.
          mtp_j(k,n) = 0.
          ftv_j(k,n) = 0.
          mtv_j(k,n) = 0
        enddo
        do n = 1, ndesign
          fti_g(k,n) = 0.
          mti_g(k,n) = 0.
          ftj_g(k,n) = 0.
          mtj_g(k,n) = 0.
          ftp_g(k,n) = 0.
          mtp_g(k,n) = 0.
          ftv_g(k,n) = 0.
          mtv_g(k,n) = 0.
        enddo
      enddo

c---------------------------------------------------------
c---- clear total jet quatities for accumulation
      dqt = 0.
      djt = 0.
      det = 0.
      do n = 1, numax
        dqt_u(n) = 0.
        djt_u(n) = 0.
        det_u(n) = 0.
      enddo
      do n = 1, nvarjet
        dqt_j(n) = 0.
        djt_j(n) = 0.
        det_j(n) = 0.
      enddo

c---------------------------------------------------------
c---- calculate all near field force components

c---- inviscid(vortex), jet, viscous forces for elements, strips, surfaces,
c-    accumulate to total
c      print *, 'calling sfforc'
      call sfforc

c---- inviscid and viscous forces for segments, bodies,
c-    accumulate to total
c      print *, 'calling bdforc'
      call bdforc
      
cc---- jet forces  (hybrid near/far field)
c      call sjforc

c---------------------------------------------------------
c---- double the x,z forces, zero y force for a y symmetric case
      if(iysym.eq.1) then
       fsym(1) = 2.0
       fsym(2) = 0.
       fsym(3) = 2.0
       msym(1) = 0.
       msym(2) = 2.0
       msym(3) = 0.

       do k = 1, 3
         fti(k) = fsym(k)*fti(k)
         mti(k) = msym(k)*mti(k)
         ftj(k) = fsym(k)*ftj(k)
         mtj(k) = msym(k)*mtj(k)
         ftp(k) = fsym(k)*ftp(k)
         mtp(k) = msym(k)*mtp(k)
         ftv(k) = fsym(k)*ftv(k)
         mtv(k) = msym(k)*mtv(k)
         do n = 1, numax
           fti_u(k,n) = fsym(k) * fti_u(k,n)
           mti_u(k,n) = msym(k) * mti_u(k,n)
           ftj_u(k,n) = fsym(k) * ftj_u(k,n)
           mtj_u(k,n) = msym(k) * mtj_u(k,n)
           ftp_u(k,n) = fsym(k) * ftp_u(k,n)
           mtp_u(k,n) = msym(k) * mtp_u(k,n)
           ftv_u(k,n) = fsym(k) * ftv_u(k,n)
           mtv_u(k,n) = msym(k) * mtv_u(k,n)
         enddo
         do n = 1, ncontrol
           fti_d(k,n) = fsym(k) * fti_d(k,n)
           mti_d(k,n) = msym(k) * mti_d(k,n)
           ftj_d(k,n) = fsym(k) * ftj_d(k,n)
           mtj_d(k,n) = msym(k) * mtj_d(k,n)
           ftp_d(k,n) = fsym(k) * ftp_d(k,n)
           mtp_d(k,n) = msym(k) * mtp_d(k,n)
           ftv_d(k,n) = fsym(k) * ftv_d(k,n)
           mtv_d(k,n) = msym(k) * mtv_d(k,n)
         enddo
         do n = 1, nvarjet
           fti_j(k,n) = fsym(k) * fti_j(k,n)
           mti_j(k,n) = msym(k) * mti_j(k,n)
           ftj_j(k,n) = fsym(k) * ftj_j(k,n)
           mtj_j(k,n) = msym(k) * mtj_j(k,n)
           ftp_j(k,n) = fsym(k) * ftp_j(k,n)
           mtp_j(k,n) = msym(k) * mtp_j(k,n)
           ftv_j(k,n) = fsym(k) * ftv_j(k,n)
           mtv_j(k,n) = msym(k) * mtv_j(k,n)
         enddo
         do n = 1, ndesign
           fti_g(k,n) = fsym(k) * fti_g(k,n)
           mti_g(k,n) = msym(k) * mti_g(k,n)
           ftj_g(k,n) = fsym(k) * ftj_g(k,n)
           mtj_g(k,n) = msym(k) * mtj_g(k,n)
           ftp_g(k,n) = fsym(k) * ftp_g(k,n)
           mtp_g(k,n) = msym(k) * mtp_g(k,n)
           ftv_g(k,n) = fsym(k) * ftv_g(k,n)
           mtv_g(k,n) = msym(k) * mtv_g(k,n)
         enddo
       enddo
      endif

c---------------------------------------------------------
c---- double the jet mass, momentum and energy totals
      if(iysym.eq.1) then
        dqt = 2.0 * dqt
        djt = 2.0 * djt
        det = 2.0 * det
        do n = 1, numax
          dqt_u(n) = 2.0 * dqt_u(n)
          djt_u(n) = 2.0 * djt_u(n)
          det_u(n) = 2.0 * det_u(n)
        enddo
        do n = 1, nvarjet
          dqt_j(n) = 2.0 * dqt_j(n)
          djt_j(n) = 2.0 * djt_j(n)
          det_j(n) = 2.0 * det_j(n)
        enddo
      endif
      
c---------------------------------------------------------
c---- add baseline profile drag force along vbar, applied at xyzref
      vsq = vbar(1)**2 + vbar(2)**2 + vbar(3)**2
      vmag = sqrt(vsq)
      vmag_u(1) = vbar(1)/vmag
      vmag_u(2) = vbar(2)/vmag
      vmag_u(3) = vbar(3)/vmag

      r(1) = xyzref(1)
      r(2) = xyzref(2)
      r(3) = xyzref(3)

      cda = cdref*sref

      do k = 1, 3
        ftv(k)     = ftv(k)     + cda*0.5*vbar(k)*vmag
        ftv_u(k,k) = ftv_u(k,k) + cda*0.5        *vmag
        ftv_u(k,1) = ftv_u(k,1) + cda*0.5*vbar(k)*vmag_u(1)
        ftv_u(k,2) = ftv_u(k,2) + cda*0.5*vbar(k)*vmag_u(2)
        ftv_u(k,3) = ftv_u(k,3) + cda*0.5*vbar(k)*vmag_u(3)
C---- Bug 03272024 HHY
C---- if baseline profile drag acts at xyzref the moments are 0.0
c        ic = icrs(k)
c        jc = jcrs(k)
c        rxv = r(ic)*vbar(jc) - r(jc)*vbar(ic)
c        mtv(k)      = mtv(k)      + cda*0.5*rxv  *vmag
c        mtv_u(k,ic) = mtv_u(k,ic) - cda*0.5*r(jc)*vmag
c        mtv_u(k,jc) = mtv_u(k,jc) + cda*0.5*r(ic)*vmag
c        mtv_u(k,1)  = mtv_u(k,1)  + cda*0.5*rxv  *vmag_u(1)
c        mtv_u(k,2)  = mtv_u(k,2)  + cda*0.5*rxv  *vmag_u(2)
c        mtv_u(k,3)  = mtv_u(k,3)  + cda*0.5*rxv  *vmag_u(3)
      enddo

c---------------------------------------------------------
c---- total forces and moments
      do k = 1, 3
        ft(k) = fti(k) + ftj(k)
     &        + ftp(k) + ftv(k)

ccccc
        qs  = 0.5*sref
c        print *,"ft(",k,") tot ",ft(k)/qs
c        print *,"ft(",k,") i j ",fti(k)/qs,ftj(k)/qs
c        print *,"ft(",k,") p v ",ftp(k)/qs,ftv(k)/qs

        mt(k) = mti(k) + mtj(k)
     &        + mtp(k) + mtv(k) 
        do n = 1, numax
          ft_u(k,n) = fti_u(k,n) + ftj_u(k,n)
     &              + ftp_u(k,n) + ftv_u(k,n)
          mt_u(k,n) = mti_u(k,n) + mtj_u(k,n)
     &              + mtp_u(k,n) + mtv_u(k,n)
        enddo                                                               
        do n = 1, ncontrol
          ft_d(k,n) = fti_d(k,n) + ftj_d(k,n)
     &              + ftp_d(k,n) + ftv_d(k,n)
          mt_d(k,n) = mti_d(k,n) + mtj_d(k,n)
     &              + mtp_d(k,n) + mtv_d(k,n)
        enddo                                                               
        do n = 1, nvarjet
          ft_j(k,n) = fti_j(k,n) + ftj_j(k,n)
     &              + ftp_j(k,n) + ftv_j(k,n)
          mt_j(k,n) = mti_j(k,n) + mtj_j(k,n)
     &              + mtp_j(k,n) + mtv_j(k,n)
        enddo                                                               
        do n = 1, ndesign
          ft_g(k,n) = fti_g(k,n) + ftj_g(k,n)
     &              + ftp_g(k,n) + ftv_g(k,n)
          mt_g(k,n) = mti_g(k,n) + mtj_g(k,n)
     &              + mtp_g(k,n) + mtv_g(k,n)
        enddo
      enddo
c
c---------------------------------------------------------
c---- set total forces in stability axes
      v13sq = vbar(1)**2 + vbar(3)**2
      v13 = sqrt(v13sq)
      sa = vbar(3)/v13
      ca = vbar(1)/v13

      sa_u(1) = -sa*vbar(1)/v13sq
      sa_u(2) = 0.
      sa_u(3) = -sa*vbar(3)/v13sq + 1.0/v13

      ca_u(1) = -ca*vbar(1)/v13sq + 1.0/v13
      ca_u(2) = 0.
      ca_u(3) = -ca*vbar(3)/v13sq

      dtot = ca*ft(1) + sa*ft(3)
      ytot =    ft(2)
      ltot = ca*ft(3) - sa*ft(1)

      do n = 1, numax
        dtot_u(n) = ca*ft_u(1,n) + sa*ft_u(3,n)
        ytot_u(n) =    ft_u(2,n)
        ltot_u(n) = ca*ft_u(3,n) - sa*ft_u(1,n)
      enddo
      do n = 1, 3
        dtot_u(n) = dtot_u(n) + ca_u(n)*ft(1) + sa_u(n)*ft(3)
        ltot_u(n) = ltot_u(n) + ca_u(n)*ft(3) - sa_u(n)*ft(1)
      enddo

      do n = 1, ncontrol
        dtot_d(n) = ca*ft_d(1,n) + sa*ft_d(3,n)
        ytot_d(n) =    ft_d(2,n)
        ltot_d(n) = ca*ft_d(3,n) - sa*ft_d(1,n)
      enddo
      do n = 1, nvarjet
        dtot_j(n) = ca*ft_j(1,n) + sa*ft_j(3,n)
        ytot_j(n) =    ft_j(2,n)
        ltot_j(n) = ca*ft_j(3,n) - sa*ft_j(1,n)
      enddo
      do n = 1, ndesign
        dtot_g(n) = ca*ft_g(1,n) + sa*ft_g(3,n)
        ytot_g(n) =    ft_g(2,n)
        ltot_g(n) = ca*ft_g(3,n) - sa*ft_g(1,n)
      enddo

c================================================================
c---- calculate total farfield (trefftz-plane) forces
      call tpforc

c---- add baseline profile drag force (this doesn't get ysym doubling)
      dffv = dffv + 0.5*cdref*sref

c---------------------------------------------------------
c---- set total trefftz-plane forces
      lff = lffi + lffj + lffb
      yff = yffi + yffj + yffb
      dff = dffi + dffj + dffb + dffv

      do n = 1, numax
        lff_u(n) = lffi_u(n) + lffj_u(n) + lffb_u(n) 
        yff_u(n) = yffi_u(n) + yffj_u(n) + yffb_u(n) 
        dff_u(n) = dffi_u(n) + dffj_u(n) + dffb_u(n) + dffv_u(n)
      enddo
      do n = 1, ncontrol
        lff_d(n) = lffi_d(n) + lffj_d(n) + lffb_d(n) 
        yff_d(n) = yffi_d(n) + yffj_d(n) + yffb_d(n) 
        dff_d(n) = dffi_d(n) + dffj_d(n) + dffb_d(n) + dffv_d(n)
      enddo
      do n = 1, nvarjet
        lff_j(n) = lffi_j(n) + lffj_j(n) + lffb_j(n) 
        yff_j(n) = yffi_j(n) + yffj_j(n) + yffb_j(n) 
        dff_j(n) = dffi_j(n) + dffj_j(n) + dffb_j(n) + dffv_j(n)
      enddo
      do n = 1, ndesign
        lff_g(n) = lffi_g(n) + lffj_g(n) + lffb_g(n) 
        yff_g(n) = yffi_g(n) + yffj_g(n) + yffb_g(n) 
        dff_g(n) = dffi_g(n) + dffj_g(n) + dffb_g(n) + dffv_g(n)
      enddo

      return
      end ! aero


      subroutine ab2vbar
c--------------------------------------------------------------------------
c     Calculates freestream vector components and sensitivities
c
c   Inputs:   
c      alfa       angle of attack (for stability-axis definition)
c      beta       sideslip angle (positive wind on right cheek facing fwd)
c
c   Outputs:  
c      vbar(3)    velocity components of free stream
c      vbar_a(3)  dvbar()/dalfa
c      vbar_b(3)  dvbar()/dbeta
c--------------------------------------------------------------------------
      include 'jvl.inc'

      sina = sin(alfa)
      cosa = cos(alfa)
      sinb = sin(beta)
      cosb = cos(beta)

      vbar(1) =  cosa*cosb
      vbar(2) =      -sinb
      vbar(3) =  sina*cosb

      vbar_a(1) = -sina*cosb
      vbar_a(2) =  0.
      vbar_a(3) =  cosa*cosb

      vbar_b(1) = -cosa*sinb
      vbar_b(2) =      -cosb
      vbar_b(3) = -sina*sinb

      return
      end ! ab2vbar


      subroutine vbar2ab
c--------------------------------------------------------------------------
c     Calculates alfa,beta from freestream vector components
c
c   Inputs:   
c      vbar(3)    velocity components of free stream
c
c   Outputs:  
c      alfa       angle of attack
c      beta       sideslip angle
c      alfa_u(3)  dalfa/dvbar()
c      beta_u(3)  dbeta/dvbar()
c--------------------------------------------------------------------------
      include 'jvl.inc'

      alfa = atan2( vbar(3) , vbar(1) )
      alfa_u(1) = -vbar(3) / (vbar(1)**2 + vbar(3)**2)
      alfa_u(2) = 0.
      alfa_u(3) =  vbar(1) / (vbar(1)**2 + vbar(3)**2)

      v13 = sqrt(vbar(1)**2 + vbar(3)**2)
      vsq = vbar(1)**2 + vbar(2)**2 + vbar(3)**2
      beta = atan2( -vbar(2) , v13 )
      beta_u(1) = vbar(1)*vbar(2)/(v13*vsq)
      beta_u(2) = -v13/vsq
      beta_u(3) = vbar(3)*vbar(2)/(v13*vsq)

      return
      end ! vbar2ab
