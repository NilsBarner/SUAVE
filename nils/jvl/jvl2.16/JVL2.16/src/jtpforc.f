C***********************************************************************
C    Module:  jtpforc.f
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

      subroutine tpforc
c------------------------------------------------------------------------
c     Calculate the far-field forces on the configuration using
c            a Trefftz Plane method.
c
c  Inputs:
c     geometry data from labelled commons
c     cnc     strip span loading 
c
c  Outputs:
c     lff    total far-field lift
c     yff    total far-field side force
c     dff    total far-field drag
c     spanef  span efficiency
c     dwwake  far-field wake downwash at center of strip 
c
c  Comments:
c     The far-field drag is calculated using the Trefftz
c     Plane (kinetic energy integral in the far wake).
c     the span-loading cnc is all that is required, plus
c     geometry data defining the wake position.
c     since the wakes are just the horseshoe legs extending 
c     downstream from the bound legs, only the y and z 
c     coordinates are used. the normalwash in the cross-plane
c     is evaluated over the real and "image" sides.
c------------------------------------------------------------------------
      include 'jvl.inc'

      real ny, nz
      real vy_u(numax), vz_u(numax), 
     &     vy_d(ndmax), vz_d(ndmax), 
     &     vy_j(njmax), vz_j(njmax),
     &     vy_g(ngmax), vz_g(ngmax)
      real p(3,3), p_m(3,3), p_a(3,3), p_b(3,3)
      real rt1(3,nsmax),
     &     rt2(3,nsmax),
     &     rtc(3,nsmax)
      real gams(nsmax), 
     &     gams_u(nsmax,numax),
     &     gams_d(nsmax,ndmax),
     &     gams_j(nsmax,njmax),
     &     gams_g(nsmax,ngmax)

      real vinf, vinf_u(numax)
      real dvjf, dvjf_u(numax)
      real dvjet, dvjet_j(njmax)
      real vjet, vjet_u(numax), vjet_j(njmax)
      real dqp, dqp_u(numax),  dqp_j(njmax),
     &     djp, djp_u(numax),  djp_j(njmax),
     &     dep, dep_u(numax),  dep_j(njmax)

      real jvec(3),jvecsq, jvecai,
     &     jvec_u(3,numax),
     &     jvec_d(3,ndmax),
     &     jvec_j(3,njmax),
     &     jvec_g(3,ngmax)
      real jhat(3),
     &     jhat_u(3,numax),
     &     jhat_d(3,ndmax),
     &     jhat_j(3,njmax),
     &     jhat_g(3,ngmax)

      real sina_u(numax), cosa_u(numax)

      hpi = 1.0 / (2.0*pi)

      vinf = sqrt(vbar(1)**2 + vbar(2)**2 + vbar(3)**2)
      vinf_u(1) = vbar(1)/vinf
      vinf_u(2) = vbar(2)/vinf
      vinf_u(3) = vbar(3)/vinf
      vinf_u(4) = 0.
      vinf_u(5) = 0.
      vinf_u(6) = 0.

      do n = 1,numax
        sina_u(n) = 0.
        cosa_u(n) = 0.
      enddo

      call vbar2ab
      sina = sin(alfa)
      cosa = cos(alfa)
      do n = 1, 3
        sina_u(n) =  cosa*alfa_u(n)
        cosa_u(n) = -sina*alfa_u(n)
      enddo

c---- set prandtl-glauert transformation matrix
      alfat = 0.
      betat = 0.
      call pgmat(amach,alfat,betat,p,p_m,p_a,p_b)

      yoff = 2.0*ysym
      zoff = 2.0*zsym

      lffi = 0.
      yffi = 0.
      dffi = 0.
      lffj = 0.
      yffj = 0.
      dffj = 0.
      lffb = 0.
      yffb = 0.
      dffb = 0.
      dffv = 0.
      do n = 1, numax
        lffi_u(n) = 0.
        yffi_u(n) = 0.
        dffi_u(n) = 0.
        lffj_u(n) = 0.
        yffj_u(n) = 0.
        dffj_u(n) = 0.
        lffb_u(n) = 0.
        yffb_u(n) = 0.
        dffb_u(n) = 0.
        dffv_u(n) = 0.
      enddo
      do n = 1, ncontrol
        lffi_d(n) = 0.
        yffi_d(n) = 0.
        dffi_d(n) = 0.
        lffj_d(n) = 0.
        yffj_d(n) = 0.
        dffj_d(n) = 0.
        lffb_d(n) = 0.
        yffb_d(n) = 0.
        dffb_d(n) = 0.
        dffv_d(n) = 0.
      enddo
      do n = 1, nvarjet
        lffi_j(n) = 0.
        yffi_j(n) = 0.
        dffi_j(n) = 0.
        lffj_j(n) = 0.
        yffj_j(n) = 0.
        dffj_j(n) = 0.
        lffb_j(n) = 0.
        yffb_j(n) = 0.
        dffb_j(n) = 0.
        dffv_j(n) = 0.
      enddo
      do n = 1, ndesign
        lffi_g(n) = 0.
        yffi_g(n) = 0.
        dffi_g(n) = 0.
        lffj_g(n) = 0.
        yffj_g(n) = 0.
        dffj_g(n) = 0.
        lffb_g(n) = 0.
        yffb_g(n) = 0.
        dffb_g(n) = 0.
        dffv_g(n) = 0.
      enddo

      dqff = 0.
      djff = 0.
      deff = 0.
      do n = 1, nvarjet
        dqff_j(n) = 0.
        djff_j(n) = 0.
        deff_j(n) = 0.
      enddo

      do js = 1, nstrip
        gams(js) = 0.
        do n = 1, numax
          gams_u(js,n) = 0.
        enddo
        do n = 1, ncontrol
          gams_d(js,n) = 0.
        enddo
        do n = 1, nvarjet
          gams_j(js,n) = 0.
        enddo
        do n = 1, ndesign
          gams_g(js,n) = 0.
        enddo

        do i = ifrsts(js), ilasts(js)
          gams(js) = gams(js) + gam(i)
          do n = 1, numax
            gams_u(js,n) = gams_u(js,n) + gam_u(i,n)
          enddo
          do n = 1, ncontrol
            gams_d(js,n) = gams_d(js,n) + gam_d(i,n)
          enddo
          do n = 1, nvarjet
            gams_j(js,n) = gams_j(js,n) + gam_j(i,n)
          enddo
          do n = 1, ndesign
            gams_g(js,n) = gams_g(js,n) + gam_g(i,n)
          enddo
        enddo

        do i = ifrstu(js), ilastu(js)
          gams(js) = gams(js) + gam(i)
          do n = 1, numax
            gams_u(js,n) = gams_u(js,n) + gam_u(i,n)
          enddo
          do n = 1, ncontrol
            gams_d(js,n) = gams_d(js,n) + gam_d(i,n)
          enddo
          do n = 1, nvarjet
            gams_j(js,n) = gams_j(js,n) + gam_j(i,n)
          enddo
          do n = 1, ndesign
            gams_g(js,n) = gams_g(js,n) + gam_g(i,n)
          enddo
        enddo

        do i = ifrstw(js), ilastw(js)
          gams(js) = gams(js) + gam(i)
          do n = 1, numax
            gams_u(js,n) = gams_u(js,n) + gam_u(i,n)
          enddo
          do n = 1, ncontrol
            gams_d(js,n) = gams_d(js,n) + gam_d(i,n)
          enddo
          do n = 1, nvarjet
            gams_j(js,n) = gams_j(js,n) + gam_j(i,n)
          enddo
          do n = 1, ndesign
            gams_g(js,n) = gams_g(js,n) + gam_g(i,n)
          enddo
        enddo
      enddo

c---- x,y,z in wind axes (y,z are then in trefftz plane)
      do js = 1, nstrip
c------ TE h.v. in strip defines trefftz plane trace
        ic = ilasts(js)
        do k = 1, 3
          rt1(k,js) = p(k,1)*rv1(1,ic)+p(k,2)*rv1(2,ic)+p(k,3)*rv1(3,ic)
          rt2(k,js) = p(k,1)*rv2(1,ic)+p(k,2)*rv2(2,ic)+p(k,3)*rv2(3,ic)
          rtc(k,js) = p(k,1)*rc (1,ic)+p(k,2)*rc (2,ic)+p(k,3)*rc (3,ic)
        enddo
      enddo

c---- normal velocity across each strip at the projected c.p. location
      do 40 js = 1, nstrip
        dxt = rt2(1,js) - rt1(1,js)
        dyt = rt2(2,js) - rt1(2,js)
        dzt = rt2(3,js) - rt1(3,js)
        dst = sqrt(dyt*dyt + dzt*dzt)

        ny = -dzt / dst
        nz =  dyt / dst
        ycntr = rtc(2,js)
        zcntr = rtc(3,js)

        zisgn = zcntr - zsym

        vy = 0.
        vz = 0.
        do n = 1, numax
          vy_u(n) = 0.
          vz_u(n) = 0.
        enddo
        do n = 1, ncontrol
          vy_d(n) = 0.
          vz_d(n) = 0.
        enddo
        do n = 1, nvarjet
          vy_j(n) = 0.
          vz_j(n) = 0.
        enddo
        do n = 1, ndesign
          vy_g(n) = 0.
          vz_g(n) = 0.
        enddo

c------ sum velocity contributions from wake vortices
        do 30 jv = 1, nstrip
          izsym = izims(jv)

          if(izsym .ne. 0) then
c--------- no influence if rtc and rt are on opposite sides of z=zsym           
           zjsgn = 0.5*(rt1(3,jv)+rt2(3,jv)) - zsym
           if(zisgn*zjsgn .le. 0.0) go to 30
          endif

          dsyz = sqrt(  (rt2(2,jv)-rt1(2,jv))**2
     &                + (rt2(3,jv)-rt1(3,jv))**2 )
          if(lscomp(isurfs(js)) .eq. lscomp(isurfs(jv))) then
ccc        rcore = 0.00001*dsyz
           rcore = 0.
          else
           rcore = max( vrcorec*chord(jv) , vrcorew*dsyz )
          endif

          rcore = 0.

          dy1 = ycntr - rt1(2,jv)
          dy2 = ycntr - rt2(2,jv)
          dz1 = zcntr - rt1(3,jv)
          dz2 = zcntr - rt2(3,jv)
          rsq1 = dy1*dy1 + dz1*dz1 + rcore**2
          rsq2 = dy2*dy2 + dz2*dz2 + rcore**2
          vy = vy + hpi*gams(jv)*( dz1/rsq1 - dz2/rsq2)
          vz = vz + hpi*gams(jv)*(-dy1/rsq1 + dy2/rsq2)
          do n = 1, numax
            vy_u(n) = vy_u(n) + hpi*gams_u(jv,n)*( dz1/rsq1 - dz2/rsq2)
            vz_u(n) = vz_u(n) + hpi*gams_u(jv,n)*(-dy1/rsq1 + dy2/rsq2)
          enddo
          do n = 1, ncontrol
            vy_d(n) = vy_d(n) + hpi*gams_d(jv,n)*( dz1/rsq1 - dz2/rsq2)
            vz_d(n) = vz_d(n) + hpi*gams_d(jv,n)*(-dy1/rsq1 + dy2/rsq2)
          enddo
          do n = 1, nvarjet
            vy_j(n) = vy_j(n) + hpi*gams_j(jv,n)*( dz1/rsq1 - dz2/rsq2)
            vz_j(n) = vz_j(n) + hpi*gams_j(jv,n)*(-dy1/rsq1 + dy2/rsq2)
          enddo
          do n = 1, ndesign
            vy_g(n) = vy_g(n) + hpi*gams_g(jv,n)*( dz1/rsq1 - dz2/rsq2)
            vz_g(n) = vz_g(n) + hpi*gams_g(jv,n)*(-dy1/rsq1 + dy2/rsq2)
          enddo

          if(izsym .ne. 0) then
            fzsym = float(izsym)
            dy1 = ycntr -       rt1(2,jv)
            dy2 = ycntr -       rt2(2,jv)
            dz1 = zcntr - (zoff-rt1(3,jv))
            dz2 = zcntr - (zoff-rt2(3,jv))
ccc         dz1 = zcntr - (zoff-rt1(3,jv)+alfa*rt1(1,jv))
ccc         dz2 = zcntr - (zoff-rt2(3,jv)+alfa*rt2(1,jv))
            rsq1 = dy1*dy1 + dz1*dz1
            rsq2 = dy2*dy2 + dz2*dz2
c---------- BUG, must recalculate rri1,rri2  (added  21 Apr 22)
            if(rsq1 .eq. 0.0) then
             rri1 = 0.
            else
             rri1 = 1.0/rsq1
            endif
            if(rsq2 .eq. 0.0) then
             rri2 = 0.
            else
             rri2 = 1.0/rsq2
            endif
c------------------------------------
            vy = vy - hpi*gams(jv)*( dz1/rsq1 - dz2/rsq2)*fzsym
            vz = vz - hpi*gams(jv)*(-dy1/rsq1 + dy2/rsq2)*fzsym
            do n = 1, numax
              vy_u(n) = vy_u(n)
     &          - hpi*gams_u(jv,n)*( dz1/rsq1 - dz2/rsq2)*fzsym
              vz_u(n) = vz_u(n)
     &          - hpi*gams_u(jv,n)*(-dy1/rsq1 + dy2/rsq2)*fzsym
            enddo
            do n = 1, ncontrol
              vy_d(n) = vy_d(n)
     &          - hpi*gams_d(jv,n)*( dz1/rsq1 - dz2/rsq2)*fzsym
              vz_d(n) = vz_d(n)
     &          - hpi*gams_d(jv,n)*(-dy1/rsq1 + dy2/rsq2)*fzsym
            enddo
            do n = 1, nvarjet
              vy_j(n) = vy_j(n)
     &          - hpi*gams_j(jv,n)*( dz1/rsq1 - dz2/rsq2)*fzsym
              vz_j(n) = vz_j(n)
     &          - hpi*gams_j(jv,n)*(-dy1/rsq1 + dy2/rsq2)*fzsym
            enddo
            do n = 1, ndesign
              vy_g(n) = vy_g(n)
     &          - hpi*gams_g(jv,n)*( dz1/rsq1 - dz2/rsq2)*fzsym
              vz_g(n) = vz_g(n)
     &          - hpi*gams_g(jv,n)*(-dy1/rsq1 + dy2/rsq2)*fzsym
            enddo
          endif 

          if(iysym .ne. 0) then
            fysym = float(iysym)
            dy1 = ycntr - (yoff-rt1(2,jv))
            dy2 = ycntr - (yoff-rt2(2,jv))
            dz1 = zcntr -       rt1(3,jv)
            dz2 = zcntr -       rt2(3,jv)
            rsq1 = dy1*dy1 + dz1*dz1
            rsq2 = dy2*dy2 + dz2*dz2
            vy = vy - hpi*gams(jv)*( dz1/rsq1 - dz2/rsq2)*fysym
            vz = vz - hpi*gams(jv)*(-dy1/rsq1 + dy2/rsq2)*fysym
            do n = 1, numax
              vy_u(n) = vy_u(n)
     &          - hpi*gams_u(jv,n)*( dz1/rsq1 - dz2/rsq2)*fysym
              vz_u(n) = vz_u(n)
     &          - hpi*gams_u(jv,n)*(-dy1/rsq1 + dy2/rsq2)*fysym
            enddo
            do n = 1, ncontrol
              vy_d(n) = vy_d(n)
     &          - hpi*gams_d(jv,n)*( dz1/rsq1 - dz2/rsq2)*fysym
              vz_d(n) = vz_d(n)
     &          - hpi*gams_d(jv,n)*(-dy1/rsq1 + dy2/rsq2)*fysym
            enddo
            do n = 1, nvarjet
              vy_j(n) = vy_j(n)
     &          - hpi*gams_j(jv,n)*( dz1/rsq1 - dz2/rsq2)*fysym
              vz_j(n) = vz_j(n)
     &          - hpi*gams_j(jv,n)*(-dy1/rsq1 + dy2/rsq2)*fysym
            enddo
            do n = 1, ndesign
              vy_g(n) = vy_g(n)
     &          - hpi*gams_g(jv,n)*( dz1/rsq1 - dz2/rsq2)*fysym
              vz_g(n) = vz_g(n)
     &          - hpi*gams_g(jv,n)*(-dy1/rsq1 + dy2/rsq2)*fysym
            enddo

            if(iysym .ne. 0 .and. izsym .ne. 0) then
              fysym = float(iysym)
              fzsym = float(izsym)
              dy1 = ycntr - (yoff-rt1(2,jv))
              dy2 = ycntr - (yoff-rt2(2,jv))
              dz1 = zcntr - (zoff-rt1(3,jv))
              dz2 = zcntr - (zoff-rt2(3,jv))
ccc           dz1 = zcntr - (zoff-rt1(3,jv)+alfa*rt1(1,jv))
ccc           dz2 = zcntr - (zoff-rt2(3,jv)+alfa*rt2(1,jv))
              rsq1 = dy1*dy1 + dz1*dz1
              rsq2 = dy2*dy2 + dz2*dz2
              vy = vy + hpi*gams(jv)*( dz1/rsq1 - dz2/rsq2)*fysym*fzsym
              vz = vz + hpi*gams(jv)*(-dy1/rsq1 + dy2/rsq2)*fysym*fzsym
              do n = 1, numax
                vy_u(n) = vy_u(n)
     &            - hpi*gams_u(jv,n)*( dz1/rsq1 - dz2/rsq2)*fysym*fzsym
                vz_u(n) = vz_u(n)
     &            - hpi*gams_u(jv,n)*(-dy1/rsq1 + dy2/rsq2)*fysym*fzsym
              enddo
              do n = 1, ncontrol
                vy_d(n) = vy_d(n)
     &            - hpi*gams_d(jv,n)*( dz1/rsq1 - dz2/rsq2)*fysym*fzsym
                vz_d(n) = vz_d(n)
     &            - hpi*gams_d(jv,n)*(-dy1/rsq1 + dy2/rsq2)*fysym*fzsym
              enddo
              do n = 1, nvarjet
                vy_j(n) = vy_j(n)
     &            - hpi*gams_j(jv,n)*( dz1/rsq1 - dz2/rsq2)*fysym*fzsym
                vz_j(n) = vz_j(n)
     &            - hpi*gams_j(jv,n)*(-dy1/rsq1 + dy2/rsq2)*fysym*fzsym
              enddo
              do n = 1, ndesign
                vy_g(n) = vy_g(n)
     &            - hpi*gams_g(jv,n)*( dz1/rsq1 - dz2/rsq2)*fysym*fzsym
                vz_g(n) = vz_g(n)
     &            - hpi*gams_g(jv,n)*(-dy1/rsq1 + dy2/rsq2)*fysym*fzsym
              enddo
            endif
          endif
   30   continue ! next strip jv


        isurf = isurfs(js)
        if(lfload(isurf)) then
c------- add on forces from this strip 
         dwwake(js) = -(ny*vy + nz*vz)
         rho = rhos(js)
         lffi = lffi +     rho*gams(js)*          dyt*vinf
         yffi = yffi -     rho*gams(js)* dzt*vinf         
         dffi = dffi + 0.5*rho*gams(js)*(dzt*vy - dyt*vz)
         do n = 1, numax
           lffi_u(n) = lffi_u(n) + rho*gams_u(js,n)*dyt*vinf
     &                           + rho*gams(js)    *dyt*vinf_u(n)
           yffi_u(n) = yffi_u(n) - rho*gams_u(js,n)*dzt*vinf
     &                           - rho*gams(js)    *dzt*vinf_u(n)
           dffi_u(n) = dffi_u(n)
     &          + 0.5*(  rho*gams_u(js,n)*(dzt*vy      - dyt*vz     )
     &                 + rho*gams(js)    *(dzt*vy_u(n) - dyt*vz_u(n)) )
         enddo
         do n = 1, ncontrol
           lffi_d(n) = lffi_d(n) + rho*gams_d(js,n)*dyt*vinf
           yffi_d(n) = yffi_d(n) - rho*gams_d(js,n)*dzt*vinf
           dffi_d(n) = dffi_d(n)
     &          + 0.5*(  rho*gams_d(js,n)*(dzt*vy      - dyt*vz     )
     &                 + rho*gams(js)    *(dzt*vy_d(n) - dyt*vz_d(n)) )
         enddo
         do n = 1, nvarjet
           lffi_j(n) = lffi_j(n) + rho*gams_j(js,n)*dyt*vinf
           yffi_j(n) = yffi_j(n) - rho*gams_j(js,n)*dzt*vinf
           dffi_j(n) = dffi_j(n)
     &          + 0.5*(  rho*gams_j(js,n)*(dzt*vy      - dyt*vz     )
     &                 + rho*gams(js)    *(dzt*vy_j(n) - dyt*vz_j(n)) )
         enddo
         do n = 1, ndesign
           lffi_g(n) = lffi_g(n) + rho*gams_g(js,n)*dyt*vinf
           yffi_g(n) = yffi_g(n) - rho*gams_g(js,n)*dzt*vinf
           dffi_g(n) = dffi_g(n)
     &          + 0.5*(  rho*gams_g(js,n)*(dzt*vy      - dyt*vz     )
     &                 + rho*gams(js)    *(dzt*vy_g(n) - dyt*vz_g(n)) )
         enddo

         sstrip = chord(js)*wstrip(js)
         cdv = cdv_strp(js)
         dffv = dffv + 0.5*rho*cdv*sstrip
        endif

c------ skip jet sheet calcs if there's no jet sheet
        if(ifrstw(js) .eq. 0) go to 40

c------ jet velocity (if any) for this strip
        dvjf = vinf**dvjexp
        do n = 1, numax
          dvjf_u(n) = dvjexp*(dvjf/vinf)*vinf_u(n)
        enddo

        dvjet = 0.
        do n = 1, nvarjet
          dvjet = dvjet + deljet(n)*gjstrp(js,n)
          dvjet_j(n) =              gjstrp(js,n)
        enddo

        vjet = vinf + dvjet*dvjf

        do n = 1, numax
          vjet_u(n) = vinf_u(n) + dvjet*dvjf_u(n)
        enddo
        do n = 1, nvarjet
          vjet_j(n) = dvjet_j(n)*dvjf
        enddo

        vjmin = 0.1
        if(vjet .lt. vjmin) then
         write(*,*) '?  Vjet <', vjmin
         vjet = vjmin
        endif

c------ jet height, including jet contraction from propulsor-disk height
        rhjet = 1.0
        fh = fhstrp(js)
        hjet = hdstrp(js)*(1.0 + 0.5*fh*(vinf/vjet - 1.0))
        hjet_vinf =   hdstrp(js)*0.5*fh      /vjet
        hjet_vjet =  -hdstrp(js)*0.5*fh* vinf/vjet**2

c------ mass, momentum, energy flow excesses (per span)
        rho = rhos(js)   ! rhoinf/rhoref
        dqp      = (rhjet*vjet    - vinf   )*hjet     * rho
        djp      = (rhjet*vjet**2 - vinf**2)*hjet     * rho
        dep      = (rhjet*vjet**3 - vinf**3)*hjet*0.5 * rho

        dqp_vjet =  rhjet                   *hjet      * rho
     &           + (rhjet*vjet    - vinf   )*hjet_vjet * rho
        djp_vjet =  rhjet*vjet*2.0          *hjet      * rho
     &           + (rhjet*vjet**2 - vinf**2)*hjet_vjet * rho
        dep_vjet =  rhjet*vjet**2 * 3.0     *hjet      * 0.5 * rho
     &           + (rhjet*vjet**3 - vinf**3)*hjet_vjet * 0.5 * rho

        dqp_vinf = (              - 1.0    )*hjet      * rho
     &           + (rhjet*vjet    - vinf   )*hjet_vinf * rho
        djp_vinf = (              - vinf*2.)*hjet      * rho
     &           + (rhjet*vjet**2 - vinf**2)*hjet_vinf * rho
        dep_vinf = (         - vinf**2 * 3.)*hjet      * 0.5 * rho
     &           + (rhjet*vjet**3 - vinf**3)*hjet_vinf * 0.5 * rho
        do n = 1, numax
          dqp_u(n) = dqp_vjet*vjet_u(n) + dqp_vinf*vinf_u(n)
          djp_u(n) = djp_vjet*vjet_u(n) + djp_vinf*vinf_u(n)
          dep_u(n) = dep_vjet*vjet_u(n) + dep_vinf*vinf_u(n)
        enddo
        do n = 1, nvarjet
          dqp_j(n) = dqp_vjet*vjet_j(n)
          djp_j(n) = djp_vjet*vjet_j(n)
          dep_j(n) = dep_vjet*vjet_j(n)
        enddo

c------ jet vector in trefftz plane
        jvec(1) = p(1,1)
        jvec(2) = p(1,2) + vy/vinf
        jvec(3) = p(1,3) + vz/vinf

        do n = 1, numax
          jvec_u(1,n) = 0.
          jvec_u(2,n) = vy_u(n)/vinf - (vy/vinf**2)*vinf_u(n)
          jvec_u(3,n) = vz_u(n)/vinf - (vz/vinf**2)*vinf_u(n)
        enddo
        do n = 1, ncontrol
          jvec_d(1,n) = 0.
          jvec_d(2,n) = vy_d(n)/vinf
          jvec_d(3,n) = vz_d(n)/vinf
        enddo
        do n = 1, nvarjet
          jvec_j(1,n) = 0.
          jvec_j(2,n) = vy_j(n)/vinf
          jvec_j(3,n) = vz_j(n)/vinf
        enddo
        do n = 1, ndesign
          jvec_g(1,n) = 0.
          jvec_g(2,n) = vy_g(n)/vinf
          jvec_g(3,n) = vz_g(n)/vinf
        enddo

        jvecsq = jvec(1)**2 + jvec(2)**2 + jvec(3)**2
        if(jvecsq .ne. 0.0) then
         jvecai = 1.0 / sqrt(jvecsq)
        else
         jvecai = 1.0
        endif

        do k = 1, 3
          jhat(k) = jvec(k)*jvecai
          do n = 1, numax
            jhat_u(k,n) = jvec_u(k,n)*jvecai
     &                  - jvec(k)*( jvec(1)*jvec_u(1,n)
     &                            + jvec(2)*jvec_u(2,n)
     &                            + jvec(3)*jvec_u(3,n) )*jvecai**3
          enddo
          do n = 1, ncontrol
            jhat_d(k,n) = jvec_d(k,n)*jvecai
     &                  - jvec(k)*( jvec(1)*jvec_d(1,n)
     &                            + jvec(2)*jvec_d(2,n)
     &                            + jvec(3)*jvec_d(3,n) )*jvecai**3
          enddo
          do n = 1, nvarjet
            jhat_j(k,n) = jvec_j(k,n)*jvecai
     &                  - jvec(k)*( jvec(1)*jvec_j(1,n)
     &                            + jvec(2)*jvec_j(2,n)
     &                            + jvec(3)*jvec_j(3,n) )*jvecai**3
          enddo
          do n = 1, ndesign
            jhat_g(k,n) = jvec_g(k,n)*jvecai
     &                  - jvec(k)*( jvec(1)*jvec_g(1,n)
     &                            + jvec(2)*jvec_g(2,n)
     &                            + jvec(3)*jvec_g(3,n) )*jvecai**3
          enddo
        enddo


        lffj = lffj - (djp*jhat(3)           )*dst
        yffj = yffj - (djp*jhat(2)           )*dst
        dffj = dffj - (djp*jhat(1) - dqp*vinf)*dst

        dqff  = dqff  +  dqp*dst      
        djff  = djff  +  djp*dst
        deff  = deff  +  dep*dst

        do n = 1, numax
          lffj_u(n) = lffj_u(n) - ( djp     *jhat_u(3,n)
     &                            + djp_u(n)*jhat(3)    )*dst
          yffj_u(n) = yffj_u(n) - ( djp     *jhat_u(2,n)
     &                            + djp_u(n)*jhat(2)    )*dst
          dffj_u(n) = dffj_u(n) - ( djp     *jhat_u(1,n)
     &                            + djp_u(n)*jhat(1)
     &                            - dqp     *vinf_u(n)
     &                            - dqp_u(n)*vinf       )*dst
          dqff_u(n) = dqff_u(n) + dqp_u(n)*dst    
          djff_u(n) = djff_u(n) + djp_u(n)*dst
          deff_u(n) = deff_u(n) + dep_u(n)*dst    
        enddo
        do n = 1, ncontrol
          lffj_d(n) = lffj_d(n) - (djp*jhat_d(3,n) )*dst
          yffj_d(n) = yffj_d(n) - (djp*jhat_d(2,n) )*dst
          dffj_d(n) = dffj_d(n) - (djp*jhat_d(1,n) )*dst
        enddo
        do n = 1, nvarjet
          lffj_j(n) = lffj_j(n) - ( djp     *jhat_j(3,n)
     &                            + djp_j(n)*jhat(3)    )*dst
          yffj_j(n) = yffj_j(n) - ( djp     *jhat_j(2,n)
     &                            + djp_j(n)*jhat(2)    )*dst
          dffj_j(n) = dffj_j(n) - ( djp     *jhat_j(1,n)
     &                            + djp_j(n)*jhat(1)
     &                            - dqp_j(n)*vinf       )*dst

          dqff_j(n) = dqff_j(n) + dqp_j(n)*dst    
          djff_j(n) = djff_j(n) + djp_j(n)*dst
          deff_j(n) = deff_j(n) + dep_j(n)*dst    
        enddo
        do n = 1, ndesign
          lffj_g(n) = lffj_g(n) - (djp*jhat_g(3,n) )*dst
          yffj_g(n) = yffj_g(n) - (djp*jhat_g(2,n) )*dst
          dffj_g(n) = dffj_g(n) - (djp*jhat_g(1,n) )*dst
        enddo

   40 continue ! strip js

c---- use near-field body lift,sideforce to obtain TP forces
      do 60 ib = 1, nbody
        l = llast(ib)
        if(iysym .eq. 0) then
         area = pi*radl(l)**2
        else
         area = 0.5*pi*radl(l)**2
        endif

        if(area .lt. sref*1.0e-8) then
          ainv = 0.0
        else
          ainv = 1.0/area
        endif

        dlff = cosa*fbi(3,ib) - sina*fbi(1,ib)
        dyff =      fbi(2,ib)

        ddff = (dlff**2 + dyff**2) * ainv / (4.0*vinf**2)
        ddff_dlff =    2.0*dlff    * ainv / (4.0*vinf**2)
        ddff_dyff =    2.0*dyff    * ainv / (4.0*vinf**2)
        ddff_vinf = -2.0*ddff/vinf

        lffb = lffb + dlff
        yffb = yffb + dyff
        dffb = dffb + ddff

        do n = 1, numax
          dlff_u = cosa     *fbi_u(3,ib,n) - sina     *fbi_u(1,ib,n)
     &         + ( cosa_u(n)*fbi(3,ib)     - sina_u(n)*fbi(1,ib)    )
          dyff_u =      fbi_u(2,ib,n)
          ddff_u = ddff_dlff*dlff_u
     &           + ddff_dyff*dyff_u
     &           + ddff_vinf*vinf_u(n)
          lffb_u(n) = lffb_u(n) + dlff_u
          yffb_u(n) = yffb_u(n) + dyff_u
          dffb_u(n) = dffb_u(n) + ddff_u
        enddo

        do n = 1, ncontrol
          dlff_d = cosa*fbi_d(3,ib,n) - sina*fbi_d(1,ib,n)
          dyff_d =      fbi_d(2,ib,n)
          ddff_d = ddff_dlff*dlff_d
     &           + ddff_dyff*dyff_d
          lffb_d(n) = lffb_d(n) + dlff_d
          yffb_d(n) = yffb_d(n) + dyff_d
          dffb_d(n) = dffb_d(n) + ddff_d
        enddo
        do n = 1, nvarjet
          dlff_j = cosa*fbi_j(3,ib,n) - sina*fbi_j(1,ib,n)
          dyff_j =      fbi_j(2,ib,n)
          ddff_j = ddff_dlff*dlff_j
     &           + ddff_dyff*dyff_j
          lffb_j(n) = lffb_j(n) + dlff_j
          yffb_j(n) = yffb_j(n) + dyff_j
          dffb_j(n) = dffb_j(n) + ddff_j
        enddo
        do n = 1, ndesign
          dlff_g = cosa*fbi_g(3,ib,n) - sina*fbi_g(1,ib,n)
          dyff_g =      fbi_g(2,ib,n)
          ddff_g = ddff_dlff*dlff_g
     &           + ddff_dyff*dyff_g
          lffb_g(n) = lffb_g(n) + dlff_g
          yffb_g(n) = yffb_g(n) + dyff_g
          dffb_g(n) = dffb_g(n) + ddff_g
        enddo
 60   continue

c---- double the x,z forces, zero y force for a y symmetric case
      if(iysym.eq.1) then
       lffi = 2.0 * lffi
       yffi = 0.0
       dffi = 2.0 * dffi
       lffj = 2.0 * lffj
       yffj = 0.0
       dffj = 2.0 * dffj
       lffb = 2.0 * lffb
       yffb = 0.0
       dffb = 2.0 * dffb
       dffv = 2.0 * dffv
       do n = 1, numax
         lffi_u(n) = 2.0 * lffi_u(n)
         yffi_u(n) = 0.0
         dffi_u(n) = 2.0 * dffi_u(n)
         lffj_u(n) = 2.0 * lffj_u(n)
         yffj_u(n) = 0.0
         dffj_u(n) = 2.0 * dffj_u(n)
         lffb_u(n) = 2.0 * lffb_u(n)
         yffb_u(n) = 0.0
         dffb_u(n) = 2.0 * dffb_u(n)
         dffv_u(n) = 2.0 * dffv_u(n)
       enddo
       do n = 1, ncontrol
         lffi_d(n) = 2.0 * lffi_d(n)
         yffi_d(n) = 0.0
         dffi_d(n) = 2.0 * dffi_d(n)
         lffj_d(n) = 2.0 * lffj_d(n)
         yffj_d(n) = 0.0
         dffj_d(n) = 2.0 * dffj_d(n)
         lffb_d(n) = 2.0 * lffb_d(n)
         yffb_d(n) = 0.0
         dffb_d(n) = 2.0 * dffb_d(n)
         dffv_d(n) = 2.0 * dffv_d(n)
       enddo
       do n = 1, nvarjet
         lffi_j(n) = 2.0 * lffi_j(n)
         yffi_j(n) = 0.0
         dffi_j(n) = 2.0 * dffi_j(n)
         lffj_j(n) = 2.0 * lffj_j(n)
         yffj_j(n) = 0.0
         dffj_j(n) = 2.0 * dffj_j(n)
         lffb_j(n) = 2.0 * lffb_j(n)
         yffb_j(n) = 0.0
         dffb_j(n) = 2.0 * dffb_j(n)
         dffv_j(n) = 2.0 * dffv_j(n)
       enddo
       do n = 1, ndesign
         lffi_g(n) = 2.0 * lffi_g(n)
         yffi_g(n) = 0.0
         dffi_g(n) = 2.0 * dffi_g(n)
         lffj_g(n) = 2.0 * lffj_g(n)
         yffj_g(n) = 0.0
         dffj_g(n) = 2.0 * dffj_g(n)
         lffb_g(n) = 2.0 * lffb_g(n)
         yffb_g(n) = 0.0
         dffb_g(n) = 2.0 * dffb_g(n)
         dffv_g(n) = 2.0 * dffv_g(n)
       enddo

       dqff = 2.0 * dqff
       djff = 2.0 * djff
       deff = 2.0 * deff
       do n = 1, nvarjet
         dqff_j(n) = 2.0 * dqff_j(n)
         djff_j(n) = 2.0 * djff_j(n)
         deff_j(n) = 2.0 * deff_j(n)
       enddo
      endif

c---- bsq = aspect ratio * sref
      bsq = bref**2

c---- span efficiency
      if(dffi .eq. 0.0) then
       spanef = 0.
       spanef_a = 0.
       do n = 1, numax
         spanef_u(n) = 0.
       enddo
       do n = 1, ncontrol
         spanef_d(n) = 0.
       enddo
       do n = 1, nvarjet
         spanef_j(n) = 0.
       enddo
       do n = 1, ndesign
         spanef_g(n) = 0.
       enddo

      else
       dact = dffi + dffj + djff - vinf*dqff

       spanef     = 2.0*((lffi+lffj)**2 + (yffi+yffj)**2)
     &            / ((pi*bsq + 4.0*djff) * dact * vinf**2)
       spanef_lff = 4.0*(lffi+lffj)
     &            / ((pi*bsq + 4.0*djff) * dact * vinf**2)
       spanef_yff = 4.0*(yffi+yffj)
     &            / ((pi*bsq + 4.0*djff) * dact * vinf**2)
       spanef_dff = -spanef/dact
       spanef_dj  = -spanef/dact - spanef/(pi*bsq + 4.0*djff) * 4.0
       spanef_dq  =  spanef/dact * vinf
       spanef_vi  =  spanef/dact * dqff
     &             - spanef/vinf * 2.0

       spanev     = 2.0*((lffi+lffj)**2 + (yffi+yffj)**2)
     &                              / (pi*bsq * dact * vinf**2)
       spanev_lff = 4.0*(lffi+lffj) / (pi*bsq * dact * vinf**2)
       spanev_yff = 4.0*(yffi+yffj) / (pi*bsq * dact * vinf**2)
       spanev_dff = -spanev/dact
       spanev_dj  = -spanev/dact
       spanev_dq  =  spanev/dact * vinf
       spanev_vi  =  spanev/dact * dqff
     &             - spanev/vinf * 2.0

       spanef_a = 0.
       spanev_a = 0.
       do n = 1, numax
         spanef_u(n) = spanef_lff*(lffi_u(n)+lffj_u(n))
     &               + spanef_yff*(yffi_u(n)+yffj_u(n))
     &               + spanef_dff*(dffi_u(n)+dffj_u(n))
     &               + spanef_vi*vinf_u(n)
         spanev_u(n) = spanev_lff*(lffi_u(n)+lffj_u(n))
     &               + spanev_yff*(yffi_u(n)+yffj_u(n))
     &               + spanev_dff*(dffi_u(n)+dffj_u(n))
     &               + spanev_vi*vinf_u(n)
       enddo
       do n = 1, ncontrol
         spanef_d(n) = spanef_lff*(lffi_d(n)+lffj_d(n))
     &               + spanef_yff*(yffi_d(n)+yffj_d(n))
     &               + spanef_dff*(dffi_d(n)+dffj_d(n))
         spanev_d(n) = spanev_lff*(lffi_d(n)+lffj_d(n))
     &               + spanev_yff*(yffi_d(n)+yffj_d(n))
     &               + spanev_dff*(dffi_d(n)+dffj_d(n))
       enddo
       do n = 1, nvarjet
         spanef_j(n) = spanef_lff*(lffi_j(n)+lffj_j(n))
     &               + spanef_yff*(yffi_j(n)+yffj_j(n))
     &               + spanef_dff*(dffi_j(n)+dffj_j(n))
     &               + spanef_dj*djff_j(n)
     &               + spanef_dq*dqff_j(n)
         spanev_j(n) = spanev_lff*(lffi_j(n)+lffj_j(n))
     &               + spanev_yff*(yffi_j(n)+yffj_j(n))
     &               + spanev_dff*(dffi_j(n)+dffj_j(n))
     &               + spanev_dj*djff_j(n)
     &               + spanev_dq*dqff_j(n)
       enddo
       do n = 1, ndesign
         spanef_g(n) = spanef_lff*(lffi_g(n)+lffj_g(n))
     &               + spanef_yff*(yffi_g(n)+yffj_g(n))
     &               + spanef_dff*(dffi_g(n)+dffj_g(n))
         spanev_g(n) = spanev_lff*(lffi_g(n)+lffj_g(n))
     &               + spanev_yff*(yffi_g(n)+yffj_g(n))
     &               + spanev_dff*(dffi_g(n)+dffj_g(n))
       enddo
      endif

      return
      end ! tpforc



      subroutine pgmat(mach,alfa,beta,p,p_m,p_a,p_b)
c-------------------------------------------------------
c     calculates prandtl-glauert transformation matrix.
c      
c      xi      [       ] x
c              [       ]  
c      eta  =  [   p   ] y
c              [       ]
c      zeta    [       ] z
c
c-------------------------------------------------------

      real mach, alfa, beta
      real p(3,3), p_m(3,3), p_a(3,3), p_b(3,3)

      binv = 1.0 / sqrt(1.0 - mach**2)
      bi_m = mach * binv**3

      sina = sin(alfa)
      cosa = cos(alfa)

      sinb = sin(beta)
      cosb = cos(beta)


      p(1,1) =  cosa*cosb*binv
      p(1,2) =      -sinb*binv
      p(1,3) =  sina*cosb*binv

      p(2,1) =  cosa*sinb
      p(2,2) =       cosb
      p(2,3) =  sina*sinb

      p(3,1) = -sina
      p(3,2) = 0.
      p(3,3) =  cosa


      p_m(1,1) =  cosa*cosb*bi_m
      p_m(1,2) =      -sinb*bi_m
      p_m(1,3) =  sina*cosb*bi_m

      p_m(2,1) = 0.
      p_m(2,2) = 0.
      p_m(2,3) = 0.

      p_m(3,1) = 0.
      p_m(3,2) = 0.
      p_m(3,3) = 0.


      p_a(1,1) = -sina*cosb*binv
      p_a(1,2) = 0.
      p_a(1,3) =  cosa*cosb*binv

      p_a(2,1) = -sina*sinb
      p_a(2,2) = 0.
      p_a(2,3) =  cosa*sinb

      p_a(3,1) = -cosa
      p_a(3,2) = 0.
      p_a(3,3) = -sina


      p_b(1,1) = -cosa*sinb*binv
      p_b(1,2) =      -cosb*binv
      p_b(1,3) = -sina*sinb*binv

      p_b(2,1) =  cosa*cosb
      p_b(2,2) =      -sinb
      p_b(2,3) =  sina*cosb

      p_b(3,1) = 0.
      p_b(3,2) = 0.
      p_b(3,3) = 0.

      return
      end ! pgmat
