C***********************************************************************
C    Module:  jsfforc.f
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

     
      subroutine sfforc
c-------------------------------------------------------------------------
c     Calculates near-field forces and moments on the vortices and
c     body segments, and on the overall strips, surfaces, and bodies.
c
c  Inputs:
c    gam(.)     h.v. strengths
c    gam_u(..)  dgam/d(vbar,wbar)
c    gam_d(..)  dgam/d(delcon)
c    gam_j(..)  dgam/d(deljet)
c    gam_g(..)  dgam/d(deldes)
c    vbar(.)    freestream velocity components
c    wbar(.)    roll,pitch,yaw rates
c          
c  Outputs:  
c    dcpv     vortex sheet  loadings  rho V gamma / q
c    dcpj     jet-curvature loadings  DJ` kappa / q
c    ltot     total L
c    ytot     total Y
c    dtot     total D
c    ft(.)    total F
c    mt(.)    total M
c-------------------------------------------------------------------------
      include 'jvl.inc'
c
      real rc4(3), rrot(3)
      real wxr(3), wxr_w(3,3)
      real vtot(3),
     &     vtot_u(3,numax), 
     &     vtot_d(3,ndmax), 
     &     vtot_j(3,njmax), 
     &     vtot_g(3,ngmax)
      real vtota,
     &     vtota_u(numax)
      real dl(3), r(3)
      real rh(3), f(3), m(3)
      real vxl(3),
     &     vxl_u(3,numax),
     &     vxl_d(3,ndmax),
     &     vxl_j(3,njmax),
     &     vxl_g(3,ngmax)
      real fi(3),
     &     fi_u(3,numax),
     &     fi_d(3,ndmax),
     &     fi_j(3,njmax),
     &     fi_g(3,ngmax)
      real fj(3),
     &     fj_u(3,numax),
     &     fj_d(3,ndmax),
     &     fj_j(3,njmax),
     &     fj_g(3,ngmax)
      real spn(3)
      real udrag(3),
     &     udrag_u(3,numax)
      real ulift(3), 
     &     ulift_u(3,numax)
      real dxs(3),
     &     dxs_u(3,numax)
      real dxsa,
     &     dxsa_u(numax)

      real enctot(3)

      real vinf, vinf_u(numax)
      real dvjf, dvjf_u(numax)
      real dvjet, dvjet_j(njmax)
      real vjet, vjet_u(numax), vjet_j(njmax)
      real dqp, dqp_u(numax),  dqp_j(njmax),
     &     djp, djp_u(numax),  djp_j(njmax),
     &     dep, dep_u(numax),  dep_j(njmax)
      real fp(3),
     &     fp_u(3,numax), fp_j(3,njmax)

      real dvcn, dvcn_u(numax), dvcn_d(ndmax), 
     &           dvcn_j(njmax), dvcn_g(ngmax)
      real dwcn, dwcn_u(numax), dwcn_d(ndmax),
     &           dwcn_j(njmax), dwcn_g(ngmax)
      real delf, delf_d(ndmax), delf_g(ngmax)
      real ddte, ddte_d(ndmax), ddte_g(ndmax)
      real gamj,
     &     gamj_u(numax),
     &     gamj_d(ndmax),
     &     gamj_j(njmax),
     &     gamj_g(ngmax)

      real clv,
     &     clv_u(numax),
     &     clv_d(ndmax),
     &     clv_j(njmax),
     &     clv_g(ngmax)
c
c---- indices for forming cross-products
      integer icrs(3), jcrs(3)
      data icrs / 2, 3, 1 / , jcrs / 3, 1, 2 /

c---- clear element vortex and jet DCp loadings
      do i = 1, nvor
        dcpv(i) = 0.
        dcpj(i) = 0.

        do n = 1, numax
          dcpv_u(i,n) = 0.
          dcpj_u(i,n) = 0.
        enddo

        do n = 1, ncontrol
          dcpv_d(i,n) = 0.
          dcpj_d(i,n) = 0.
        enddo

        do n = 1, nvarjet
          dcpv_j(i,n) = 0.
          dcpj_j(i,n) = 0.
        enddo

        do n = 1, ndesign
          dcpv_g(i,n) = 0.
          dcpj_g(i,n) = 0.
        enddo
      enddo

      vinf = sqrt(vbar(1)**2 + vbar(2)**2 + vbar(3)**2)
      vinf_u(1) = vbar(1)/vinf
      vinf_u(2) = vbar(2)/vinf
      vinf_u(3) = vbar(3)/vinf
      vinf_u(4) = 0.
      vinf_u(5) = 0.
      vinf_u(6) = 0.

c=====================================================================
c---- strip forces
      do 100 js = 1, nstrip
        sstrip = chord(js)*wstrip(js)

        cr = chord(js)
        sr = sstrip

c------ define local strip lift and drag directions
c------ the "spanwise" vector is cross product of strip normal with x chordline 
        spn(1) =  0.0
        spn(2) =  ensz(js)
        spn(3) = -ensy(js)

        do k = 1, 3
          udrag(k) = vbar(k)
          do n = 1, numax
            udrag_u(k,n) = 0.
          enddo
          udrag_u(k,k) = 1.0
        enddo          

c------ strip lift direction is vector product of "stream" and spanwise vector
        do k = 1, 3
          ic = icrs(k)
          jc = jcrs(k)
          dxs(k) = udrag(ic)*spn(jc) - udrag(jc)*spn(ic)
          do n = 1, numax
            dxs_u(k,n) = udrag_u(ic,n)*spn(jc) - udrag_u(jc,n)*spn(ic)
          enddo
        enddo

        dxsa = sqrt(dxs(1)**2 + dxs(2)**2 + dxs(3)**2)
        if(dxsa .eq. 0.0) then
         ulift(1) = 0.
         ulift(2) = 0.
         ulift(3) = 1.0
         do n = 1, numax
           ulift_u(1,n) = 0.
           ulift_u(2,n) = 0.
           ulift_u(3,n) = 0.
         enddo

        else
         do n = 1, numax
           dxsa_u(n) = ( dxs(1)*dxs_u(1,n)
     &                 + dxs(2)*dxs_u(2,n)
     &                 + dxs(3)*dxs_u(3,n) ) / dxsa
         enddo
         do k = 1, 3
           ulift(k) = dxs(k)/dxsa
           do n = 1, numax
             ulift_u(k,n) = (dxs_u(k,n) - ulift(k)*dxsa_u(n))/dxsa
           enddo
         enddo

        endif

c------ use the strip 1/4 chord location for strip moments
        rc4(1) = rle(1,js) + 0.25*cr
        rc4(2) = rle(2,js)
        rc4(3) = rle(3,js)

c------ clear strip forces and moments for accumulation
        do k = 1, 3
          fsi(k,js) = 0.
          msi(k,js) = 0.
          fsj(k,js) = 0.
          msj(k,js) = 0.
          fsp(k,js) = 0.
          msp(k,js) = 0.
          fsv(k,js) = 0.
          msv(k,js) = 0.
          do n = 1, numax
            fsi_u(k,js,n) = 0.
            msi_u(k,js,n) = 0.
            fsj_u(k,js,n) = 0.
            msj_u(k,js,n) = 0.
            fsp_u(k,js,n) = 0.
            msp_u(k,js,n) = 0.
            fsv_u(k,js,n) = 0.
            msv_u(k,js,n) = 0.
          enddo
          do n = 1, ncontrol
            fsi_d(k,js,n) = 0.
            msi_d(k,js,n) = 0.
            fsj_d(k,js,n) = 0.
            msj_d(k,js,n) = 0.
            fsp_d(k,js,n) = 0.
            msp_d(k,js,n) = 0.
            fsv_d(k,js,n) = 0.
            msv_d(k,js,n) = 0.
          enddo
          do n = 1, nvarjet
            fsi_j(k,js,n) = 0.
            msi_j(k,js,n) = 0.
            fsj_j(k,js,n) = 0.
            msj_j(k,js,n) = 0.
            fsp_j(k,js,n) = 0.
            msp_j(k,js,n) = 0.
            fsv_j(k,js,n) = 0.
            msv_j(k,js,n) = 0.
          enddo
          do n = 1, ndesign
            fsi_g(k,js,n) = 0.
            msi_g(k,js,n) = 0.
            fsj_g(k,js,n) = 0.
            msj_g(k,js,n) = 0.
            fsp_g(k,js,n) = 0.
            msp_g(k,js,n) = 0.
            fsv_g(k,js,n) = 0.
            msv_g(k,js,n) = 0.
          enddo
        enddo

        do l = 1, ncontrol
          chinge(l) = 0.
          do n = 1, numax
            chinge_u(l,n) = 0.
          enddo
          do n = 1, ncontrol
            chinge_d(l,n) = 0.
          enddo
          do n = 1, nvarjet
            chinge_j(l,n) = 0.
          enddo
          do n = 1, ndesign
            chinge_g(l,n) = 0.
          enddo
        enddo

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
         write(*,*) '****  Vjet < ', vjmin
         write(*,*) '****  Vjet reset to ', vjmin
         vjet = vjmin
        endif

c------ jet height, including jet contraction from propulsor-disk height
        rhjet = 1.0      ! rhojet/rhoinf
        fh = fhstrp(js)
        hjet = hdstrp(js)*(1.0 + 0.5*fh*(vinf/vjet - 1.0))
        hjet_vinf =   hdstrp(js)*0.5*fh      /vjet
        hjet_vjet =  -hdstrp(js)*0.5*fh* vinf/vjet**2
c---HHY
cc        write(*,*) 'vjet0, hjet0, dq0 ', vjet,hdstrp(js),vjet*hdstrp(js)
cc        write(*,*) 'vjet1, hjet,  dq1 ', vjet,hjet,vjet*hjet
        
c------ mass, momentum, energy flow excesses (per span)
        rho = rhos(js)   ! rhoinf/rhoref
        dqp      = (rhjet*vjet    - vinf   )*hjet      * rho
        djp      = (rhjet*vjet**2 - vinf**2)*hjet      * rho
        dep      = (rhjet*vjet**3 - vinf**3)*hjet*0.5  * rho

        dqp_vjet =  rhjet                   *hjet      * rho
     &           + (rhjet*vjet    - vinf   )*hjet_vjet * rho
        djp_vjet =  rhjet*vjet*2.0          *hjet      * rho
     &           + (rhjet*vjet**2 - vinf**2)*hjet_vjet * rho
        dep_vjet =  rhjet*vjet**2 * 3.0     *hjet      * 0.5 * rho
     &           + (rhjet*vjet**3 - vinf**3)*hjet_vjet * 0.5 * rho

        dqp_vinf = (              - 1.0    )*hjet      * rho
     &           + (rhjet*vjet    - vinf   )*hjet_vinf * rho
        djp_vinf = (         - vinf    * 2.)*hjet      * rho
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

c        djp = 5.0*deljet(1)
c        djp_j(1) = 5.0

        do 20 isuw = 1, 3

        if    (isuw .eq. 1) then
         ifrst = ifrsts(js)
         ilast = ilasts(js)
        elseif(isuw .eq. 2) then
         ifrst = ifrstu(js)
         ilast = ilastu(js)
        else
         ifrst = ifrstw(js)
         ilast = ilastw(js)
        endif
       print *, '***** js isuw ifrst ilast ',js,isuw,ifrst,ilast
      

c------ go over all elements (surface and jets) in this strip
        do 25 i = ifrst, ilast
c-------- index "i-1" of element directly upstream (-999 for first element)
          im = ijetm(i)

c-------- element area (for computing dcpv,dcpi)
          selem = dxv(i)*wstrip(js)

c-------- set up effective vortex strength from jet curvature
c
c        gamj(u,d,J,g) =  DJp * (dvc.n + dwc.n + ddte)
c
c   vc = vc(u,d,J,g)
c   wc = wc(u,d,J,g)
c   dvc = vc(i) - vc(i-1)
c   dwc = wc(i) - wc(i-1)
c   n = n(d,g)
c   DJp(J)  = (Vjet^2 - 1)*hdisk*(1 + 1/Vjet)/2
c   Vjet(J) = 1  +  sum_n [ Jn * gain(n) ]
c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          if(i .eq. ifrstw(js)) then
c--------- first point in wake jet... set TE jet deviation angle

c--------- TE flap deflection
           delf = 0.
           do n = 1, ncontrol
             delf = delf + delcon(n)*enc_d(1,im,n)
             delf_d(n) =             enc_d(1,im,n)
           enddo
           do n = 1, ndesign
             delf = delf + deldes(n)*enc_g(1,im,n)
             delf_g(n) =             enc_g(1,im,n)
           enddo

c--------- TE jet deviation angle
           ddte = dj0strp(js)
     &          + dj1strp(js)*delf
     &          + dj3strp(js)*delf**3
           do n = 1, ncontrol
             ddte_d(n) = 
     &            dj1strp(js)*delf_d(n)
     &          + dj3strp(js)*delf_d(n) * 3.0*delf**2
           enddo          
           do n = 1, ndesign
             ddte_g(n) = 
     &            dj1strp(js)*delf_g(n)
     &          + dj3strp(js)*delf_g(n) * 3.0*delf**2
           enddo          
          
          else
c--------- all other points... no jet deviation angle
           ddte = 0.
           do n = 1, ncontrol
             ddte_d(n) = 0.
           enddo          
           do n = 1, ndesign
             ddte_g(n) = 0.
           enddo          

          endif

c-------- deflected normal vector n(d,g)
          enctot(1) = enc(1,i)
          enctot(2) = enc(2,i)
          enctot(3) = enc(3,i)
          do n = 1, ncontrol
            enctot(1) = enctot(1) + enc_d(1,i,n)*delcon(n)
            enctot(2) = enctot(2) + enc_d(2,i,n)*delcon(n)
            enctot(3) = enctot(3) + enc_d(3,i,n)*delcon(n)
          enddo
          do n = 1, ndesign
            enctot(1) = enctot(1) + enc_g(1,i,n)*deldes(n)
            enctot(2) = enctot(2) + enc_g(2,i,n)*deldes(n)
            enctot(3) = enctot(3) + enc_g(3,i,n)*deldes(n)
          enddo

c-------- flow curvature as difference of perturbation-velocity directions
          if(im .lt. 0) then
c--------- no upstream point
           dvpn = 0.

           dvcn = vc(1,i)*enc(1,i)
     &          + vc(2,i)*enc(2,i)
     &          + vc(3,i)*enc(3,i)
           dwcn = wc(1,i)*enc(1,i)
     &          + wc(2,i)*enc(2,i)
     &          + wc(3,i)*enc(3,i)
           do n = 1, numax
             dvcn_u(n) = vc_u(1,i,n)*enc(1,i)
     &                 + vc_u(2,i,n)*enc(2,i)
     &                 + vc_u(3,i,n)*enc(3,i)
             dwcn_u(n) = wc_u(1,i,n)*enc(1,i)
     &                 + wc_u(2,i,n)*enc(2,i)
     &                 + wc_u(3,i,n)*enc(3,i)
           enddo
           do n = 1, ndmax
             dvcn_d(n) = vc_d(1,i,n)*enc(1,i)
     &                 + vc_d(2,i,n)*enc(2,i)
     &                 + vc_d(3,i,n)*enc(3,i)
             dwcn_d(n) = 0.
           enddo
           do n = 1, nvarjet
             dvcn_j(n) = vc_j(1,i,n)*enc(1,i)
     &                 + vc_j(2,i,n)*enc(2,i)
     &                 + vc_j(3,i,n)*enc(3,i)
             dwcn_j(n) = 0.
           enddo
           do n = 1, ndesign
             dvcn_g(n) = vc_g(1,i,n)*enc(1,i)
     &                 + vc_g(2,i,n)*enc(2,i)
     &                 + vc_g(3,i,n)*enc(3,i)
             dwcn_g(n) = 0.
           enddo

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          else
c--------- all other points except first one
           dvpn = 0.

           dvcn = (vc(1,i)-vc(1,im))*enctot(1)
     &          + (vc(2,i)-vc(2,im))*enctot(2)
     &          + (vc(3,i)-vc(3,im))*enctot(3)
           dwcn = (wc(1,i)-wc(1,im))*enctot(1)
     &          + (wc(2,i)-wc(2,im))*enctot(2)
     &          + (wc(3,i)-wc(3,im))*enctot(3)
           do n = 1, numax
             dvcn_u(n) = (vc_u(1,i,n)-vc_u(1,im,n))*enctot(1)
     &                 + (vc_u(2,i,n)-vc_u(2,im,n))*enctot(2)
     &                 + (vc_u(3,i,n)-vc_u(3,im,n))*enctot(3)
             dwcn_u(n) = (wc_u(1,i,n)-wc_u(1,im,n))*enctot(1)
     &                 + (wc_u(2,i,n)-wc_u(2,im,n))*enctot(2)
     &                 + (wc_u(3,i,n)-wc_u(3,im,n))*enctot(3)
           enddo
           do n = 1, ncontrol
             dvcn_d(n) = (vc_d(1,i,n)-vc_d(1,im,n))*enctot(1)
     &                 + (vc_d(2,i,n)-vc_d(2,im,n))*enctot(2)
     &                 + (vc_d(3,i,n)-vc_d(3,im,n))*enctot(3)
     &                 + (vc(1,i)-vc(1,im))*enc_d(1,i,n)
     &                 + (vc(2,i)-vc(2,im))*enc_d(2,i,n)
     &                 + (vc(3,i)-vc(3,im))*enc_d(3,i,n)
             dwcn_d(n) = (wc(1,i)-wc(1,im))*enc_d(1,i,n)
     &                 + (wc(2,i)-wc(2,im))*enc_d(2,i,n)
     &                 + (wc(3,i)-wc(3,im))*enc_d(3,i,n)
           enddo
           do n = 1, nvarjet
             dvcn_j(n) = (vc_j(1,i,n)-vc_j(1,im,n))*enctot(1)
     &                 + (vc_j(2,i,n)-vc_j(2,im,n))*enctot(2)
     &                 + (vc_j(3,i,n)-vc_j(3,im,n))*enctot(3)
             dwcn_j(n) = 0.
           enddo
           do n = 1, ndesign
             dvcn_g(n) = (vc_g(1,i,n)-vc_g(1,im,n))*enctot(1)
     &                 + (vc_g(2,i,n)-vc_g(2,im,n))*enctot(2)
     &                 + (vc_g(3,i,n)-vc_g(3,im,n))*enctot(3)
     &                 + (vc(1,i)-vc(1,im))*enc_g(1,i,n)
     &                 + (vc(2,i)-vc(2,im))*enc_g(2,i,n)
     &                 + (vc(3,i)-vc(3,im))*enc_g(3,i,n)
             dwcn_g(n) = (wc(1,i)-wc(1,im))*enc_g(1,i,n)
     &                 + (wc(2,i)-wc(2,im))*enc_g(2,i,n)
     &                 + (wc(3,i)-wc(3,im))*enc_g(3,i,n)
           enddo
          endif

c-------- apparent vortex strength from jet curvature
          gamj = djp*((dvcn + dwcn)/vinf + dtr*ddte)/rho

          do n = 1, numax
            gamj_u(n) =
     &       djp_u(n)*((dvcn      + dwcn     )/vinf + dtr*ddte)   /rho
     &     + djp     * (dvcn_u(n) + dwcn_u(n))/vinf               /rho
     &     - djp     * (dvcn      + dwcn     )/vinf**2 * vinf_u(n)/rho
          enddo

          do n = 1, ncontrol
            gamj_d(n) = 
     &           djp*((dvcn_d(n) + dwcn_d(n))/vinf + dtr*ddte_d(n))/rho
          enddo

          do n = 1, nvarjet
            gamj_j(n) =
     &           djp_j(n)*((dvcn      + dwcn     )/vinf + dtr*ddte)/rho
     &         + djp     * (dvcn_j(n) + dwcn_j(n))/vinf            /rho
          enddo

          do n = 1, ndesign
            gamj_g(n) =
     &           djp*((dvcn_g(n) + dwcn_g(n))/vinf + dtr*ddte_g(n))/rho
          enddo

c-------- local moment arm vector to vortex midpoint
          r(1) = rv(1,i) - xyzref(1)
          r(2) = rv(2,i) - xyzref(2)
          r(3) = rv(3,i) - xyzref(3)

c-------- vector from rotation reference point
          rrot(1) = rv(1,i) - xyzref(1)
          rrot(2) = rv(2,i) - xyzref(2)
          rrot(3) = rv(3,i) - xyzref(3)

c-------- set total effective velocity = induced + freestream + rotation
          do k = 1, 3
            ic = icrs(k)
            jc = jcrs(k)
            wxr(k) = wbar(ic)*rrot(jc)
     &             - wbar(jc)*rrot(ic)
            wxr_w(k,ic) =  rrot(jc)
            wxr_w(k,jc) = -rrot(ic)
            wxr_w(k,k) = 0.
          enddo

c-------- total body-relative velocity
          do k = 1, 3
            vtot(k) = vv(k,i) + vbar(k) - wxr(k)
c            vtot(k) =           vbar(k) - wxr(k)   ! @@@
            do n = 1, numax
              vtot_u(k,n) = vv_u(k,i,n)
c              vtot_u(k,n) = 0.           ! @@@
            enddo
            vtot_u(k,k) = vtot_u(k,k) + 1.0
            vtot_u(k,4) = vtot_u(k,4) - wxr_w(k,1)
            vtot_u(k,5) = vtot_u(k,5) - wxr_w(k,2)
            vtot_u(k,6) = vtot_u(k,6) - wxr_w(k,3)

            do n = 1, ncontrol
              vtot_d(k,n) = vv_d(k,i,n)
c             vtot_d(k,n) = 0.           ! @@@
            enddo
            do n = 1, nvarjet
              vtot_j(k,n) = vv_j(k,i,n)
c             vtot_j(k,n) = 0.           ! @@@
            enddo
            do n = 1, ndesign
              vtot_g(k,n) = vv_g(k,i,n)
c             vtot_g(k,n) = 0.           ! @@@
            enddo
          enddo

c-------- force on vortex segment is rho (vtot x gamma)
          dl(1) = rv2(1,i) - rv1(1,i)
          dl(2) = rv2(2,i) - rv1(2,i)
          dl(3) = rv2(3,i) - rv1(3,i)
c          call cross(vtot, dl, vxl)
c          do n = 1, numax
c            call cross(vtot_u(1,n), dl, vxl_u(1,n))
c          enddo
c          do n = 1, ncontrol
c            call cross(vtot_d(1,n), dl, vxl_d(1,n))
c          enddo
c          do n = 1, nvarjet
c            call cross(vtot_j(1,n), dl, vxl_j(1,n))
c          enddo
c          do n = 1, ndesign
c            call cross(vtot_g(1,n), dl, vxl_g(1,n))
c          enddo

          do k = 1, 3
            ic = icrs(k)
            jc = jcrs(k)
            vxl(k) = vtot(ic)*dl(jc)
     &             - vtot(jc)*dl(ic)
            do n = 1, numax
              vxl_u(k,n) = vtot_u(ic,n)*dl(jc)
     &                   - vtot_u(jc,n)*dl(ic)
            enddo
            do n = 1, ncontrol
              vxl_d(k,n) = vtot_d(ic,n)*dl(jc)
     &                   - vtot_d(jc,n)*dl(ic)
            enddo
            do n = 1, nvarjet
              vxl_j(k,n) = vtot_j(ic,n)*dl(jc)
     &                   - vtot_j(jc,n)*dl(ic)
            enddo
            do n = 1, ndesign
              vxl_g(k,n) = vtot_g(ic,n)*dl(jc)
     &                   - vtot_g(jc,n)*dl(ic)
            enddo
          enddo

          rho = rhov(i)
          
c-------- jet force on jets over surface     !!HHY 07052023
          if(i .ge. ifrsts(js) .and.
     &       i .le. ilasts(js)       ) then
           do k = 1, 3
            fj(k) = -rho*gamj*vxl(k)
            do n = 1, numax
              fj_u(k,n) = -rho*(gamj_u(n)*vxl(k) + gamj*vxl_u(k,n))
            enddo
            do n = 1, ncontrol
              fj_d(k,n) = -rho*(gamj_d(n)*vxl(k) + gamj*vxl_d(k,n))
            enddo
            do n = 1, nvarjet
              fj_j(k,n) = -rho*(gamj_j(n)*vxl(k) + gamj*vxl_j(k,n))
            enddo
            do n = 1, ndesign
              fj_g(k,n) = -rho*(gamj_g(n)*vxl(k) + gamj*vxl_g(k,n))
            enddo
           enddo
          else 
           do k = 1, 3
            fj(k) = 0.0
            do n = 1, numax
              fj_u(k,n) = 0.0
            enddo
            do n = 1, ncontrol
              fj_d(k,n) = 0.0
            enddo
            do n = 1, nvarjet
              fj_j(k,n) = 0.0
            enddo
            do n = 1, ndesign
              fj_g(k,n) = 0.0
            enddo
           enddo
          endif

c-------- force on lifting surface itself    
          if(i .ge. ifrsts(js) .and.
     &       i .le. ilasts(js)       ) then
           do k = 1, 3
             fi(k) = rho*gam(i)*vxl(k)
             do n = 1, numax
               fi_u(k,n) = rho*(gam_u(i,n)*vxl(k) + gam(i)*vxl_u(k,n))
             enddo
             do n = 1, ncontrol
               fi_d(k,n) = rho*(gam_d(i,n)*vxl(k) + gam(i)*vxl_d(k,n))
             enddo
             do n = 1, nvarjet
               fi_j(k,n) = rho*(gam_j(i,n)*vxl(k) + gam(i)*vxl_j(k,n))
             enddo
             do n = 1, ndesign
               fi_g(k,n) = rho*(gam_g(i,n)*vxl(k) + gam(i)*vxl_g(k,n))
             enddo
           enddo

           fn = dot(env(1,i),fi)
           dcpv(i) = 2.0*fn / (rho*selem)

           do n = 1, numax
             fn_u = dot(env(1,i),fi_u(1,n))
             dcpv_u(i,n) = 2.0*fn_u / (rho*selem)
           enddo

           do n = 1, ncontrol
             fn_d = dot(env(1,i),fi_d(1,n)) + dot(env_d(1,i,n),fi)
             dcpv_d(i,n) = 2.0*fn_d / (rho*selem)
           enddo

           do n = 1, nvarjet
             fn_j = dot(env(1,i),fi_j(1,n))
             dcpv_j(i,n) = 2.0*fn_j / (rho*selem)
           enddo

           do n = 1, ndesign
             fn_g = dot(env(1,i),fi_g(1,n)) + dot(env_g(1,i,n),fi)
             dcpv_g(i,n) = 2.0*fn_g / (rho*selem)
           enddo

c--------- accumulate strip forces and moments
           do k = 1, 3
             ic = icrs(k)
             jc = jcrs(k)
             fsi(k,js) = fsi(k,js) + fi(k)
             fsj(k,js) = fsj(k,js) + fj(k)
             msi(k,js) = msi(k,js) + ( r(ic)*fi(jc)
     &                                -r(jc)*fi(ic))
             msj(k,js) = msj(k,js) + ( r(ic)*fj(jc)
     &                                -r(jc)*fj(ic))
             do n = 1, numax
               fsi_u(k,js,n) = fsi_u(k,js,n) +  fi_u(k,n)
               fsj_u(k,js,n) = fsj_u(k,js,n) +  fj_u(k,n)
               msi_u(k,js,n) = msi_u(k,js,n) + (r(ic)*fi_u(jc,n)
     &                                         -r(jc)*fi_u(ic,n))
               msj_u(k,js,n) = msj_u(k,js,n) + (r(ic)*fj_u(jc,n)
     &                                         -r(jc)*fj_u(ic,n))
             enddo
             do n = 1, ncontrol
               fsi_d(k,js,n) = fsi_d(k,js,n) +  fi_d(k,n)
               fsj_d(k,js,n) = fsj_d(k,js,n) +  fj_d(k,n)
               msi_d(k,js,n) = msi_d(k,js,n) + (r(ic)*fi_d(jc,n)
     &                                         -r(jc)*fi_d(ic,n))
               msj_d(k,js,n) = msj_d(k,js,n) + (r(ic)*fj_d(jc,n)
     &                                         -r(jc)*fj_d(ic,n))
             enddo
             do n = 1, nvarjet
               fsi_j(k,js,n) = fsi_j(k,js,n) +  fi_j(k,n)
               fsj_j(k,js,n) = fsj_j(k,js,n) +  fj_j(k,n)
               msi_j(k,js,n) = msi_j(k,js,n) + (r(ic)*fi_j(jc,n)
     &                                         -r(jc)*fi_j(ic,n))
               msj_j(k,js,n) = msj_j(k,js,n) + (r(ic)*fj_j(jc,n)
     &                                         -r(jc)*fj_j(ic,n))
             enddo
             do n = 1, ndesign
               fsi_g(k,js,n) = fsi_g(k,js,n) +  fi_g(k,n)
               fsj_g(k,js,n) = fsj_g(k,js,n) +  fj_g(k,n)
               msi_g(k,js,n) = msi_g(k,js,n) + (r(ic)*fi_g(jc,n)
     &                                         -r(jc)*fi_g(ic,n))
               msj_g(k,js,n) = msj_g(k,js,n) + (r(ic)*fj_g(jc,n)
     &                                         -r(jc)*fj_g(ic,n))
             enddo
           enddo


c--------- hinge moments
           do l = 1, ncontrol
            dfac = dcontrol(i,l) / (sref*cref)

            rh(1) = rv(1,i) - phinge(1,js,l)
            rh(2) = rv(2,i) - phinge(2,js,l)
            rh(3) = rv(3,i) - phinge(3,js,l)

            f(1) = fi(1) + fj(1)
            f(2) = fi(2) + fj(2)
            f(3) = fi(3) + fj(3)
            call cross(rh,f,m)
            chinge(l) = chinge(l) + dot(m,vhinge(1,js,l))*dfac

            do n = 1, numax
              f(1) = fi_u(1,n) + fj_u(1,n)
              f(2) = fi_u(2,n) + fj_u(2,n)
              f(3) = fi_u(3,n) + fj_u(3,n)
              call cross(rh,f,m)
              chinge_u(l,n) = chinge_u(l,n) +dot(m,vhinge(1,js,l))*dfac
            enddo
            do n = 1, ncontrol
              f(1) = fi_d(1,n) + fj_d(1,n)
              f(2) = fi_d(2,n) + fj_d(2,n)
              f(3) = fi_d(3,n) + fj_d(3,n)
              call cross(rh,f,m)
              chinge_d(l,n) = chinge_d(l,n) +dot(m,vhinge(1,js,l))*dfac
            enddo
            do n = 1, nvarjet
              f(1) = fi_j(1,n) + fj_j(1,n)
              f(2) = fi_j(2,n) + fj_j(2,n)
              f(3) = fi_j(3,n) + fj_j(3,n)
              call cross(rh,f,m)
              chinge_j(l,n) = chinge_j(l,n) +dot(m,vhinge(1,js,l))*dfac
            enddo
            do n = 1, ndesign
              f(1) = fi_g(1,n) + fj_g(1,n)
              f(2) = fi_g(2,n) + fj_g(2,n)
              f(3) = fi_g(3,n) + fj_g(3,n)
              call cross(rh,f,m)
              chinge_g(l,n) = chinge_g(l,n) +dot(m,vhinge(1,js,l))*dfac
            enddo
           enddo

          endif

c-------- element DCp
          if(i .ge. ifrsts(js) .and.
     &       i .le. ilasts(js)      ) then
c--------- surface
           dcpsgn = 1.0
          else
c--------- jet
           dcpsgn = -1.0
          endif

          fn = dot(env(1,i),fj)
          dcpj(i) = dcpsgn*2.0*fn / (rho*selem)

          do n = 1, numax
            fn_u = dot(env(1,i),fj_u(1,n))
            dcpj_u(i,n) = dcpsgn*2.0*fn_u / (rho*selem)
          enddo

          do n = 1, ncontrol
            fn_d = dot(env(1,i),fj_d(1,n)) + dot(env_d(1,i,n),fj)
            dcpj_d(i,n) = dcpsgn*2.0*fn_d / (rho*selem)
          enddo

          do n = 1, nvarjet
            fn_j = dot(env(1,i),fj_j(1,n))
            dcpj_j(i,n) = dcpsgn*2.0*fn_j / (rho*selem)
          enddo

          do n = 1, ndesign
            fn_g = dot(env(1,i),fj_g(1,n)) + dot(env_g(1,i,n),fj)
            dcpj_g(i,n) = dcpsgn*2.0*fn_g / (rho*selem)
          enddo

c------------------------------------------------
 25     continue ! next i
 20     continue    ! next isuw
        
c------------------------------------------------------------------------
c------ if no trailing leg forces, skip it
        if(.not.ltrforce) go to 80

         xte1 = rle1(1,js) + chord1(js)
         xte2 = rle2(1,js) + chord2(js)

c------ sum forces in the strip as generated by velocity (freestream + rotation)
c-        the parts of trailing legs which lie on the surface
        do 72 i = ifrsts(js), ilasts(js)
          do 71 ileg = 1, 2

          if(ileg.eq.1) then
c--------- moment vector
           r(1) = 0.5*(rv1(1,i) + xte1) - xyzref(1)
           r(2) =      rv1(2,i)         - xyzref(2)
           r(3) =      rv1(3,i)         - xyzref(3)

c--------- vector from rotation axes
           rrot(1) = 0.5*(rv1(1,i) + xte1) - xyzref(1)
           rrot(2) =      rv1(2,i)         - xyzref(2)
           rrot(3) =      rv1(3,i)         - xyzref(3)

c--------- part of trailing leg lying on surface
           dl(1) = rv1(1,i) - xte1
           dl(2) = 0.
           dl(3) = 0.

          else
c--------- moment vector
           r(1) = 0.5*(rv2(1,i) + xte2) - xyzref(1)
           r(2) =      rv2(2,i)         - xyzref(2)
           r(3) =      rv2(3,i)         - xyzref(3)

c--------- vector from rotation axes
           rrot(1) = 0.5*(rv2(1,i) + xte2) - xyzref(1)
           rrot(2) =      rv2(2,i)         - xyzref(2)
           rrot(3) =      rv2(3,i)         - xyzref(3)

c--------- part of trailing leg lying on surface
           dl(1) = xte2 - rv2(1,i)
           dl(2) = 0.
           dl(3) = 0.
          endif

c-------- set total effective velocity = freestream + rotation
          do k = 1, 3
            ic = icrs(k)
            jc = jcrs(k)
            wxr(k) = wbar(ic)*rrot(jc)
     &             - wbar(jc)*rrot(ic)
            wxr_w(k,ic) =  rrot(jc)
            wxr_w(k,jc) = -rrot(ic)
            wxr_w(k,k) = 0.
          enddo

          do k = 1, 3
            vtot(k) = vbar(k) - wxr(k)
            do n = 1, 6
              vtot_u(k,n) = 0.
            enddo
            vtot_u(k,k) = 1.0
            vtot_u(k,4) = -wxr_w(k,1)
            vtot_u(k,5) = -wxr_w(k,2)
            vtot_u(k,6) = -wxr_w(k,3)
          enddo

          rho = rhov(i)

c-------- force-area on vortex segment is rho (vtot x gamma)
          call cross (vtot, dl, vxl)

          do n = 1, numax
            call cross(vtot_u(1,n), dl, vxl_u(1,n))
          enddo

          do k = 1, 3
            fi(k) = rho*gam(i)*vxl(k)
            do n = 1, numax
              fi_u(k,n) = rho*(gam_u(i,n)*vxl(k) + gam(i)*vxl_u(k,n))
            enddo
            do n = 1, ncontrol
              fi_d(k,n) = rho*gam_d(i,n)*vxl(k)
            enddo
            do n = 1, nvarjet
              fi_j(k,n) = rho*gam_j(i,n)*vxl(k)
            enddo
            do n = 1, ndesign
              fi_g(k,n) = rho*gam_g(i,n)*vxl(k)
            enddo
          enddo

cc-------- delta cp (loading across lifting surface) due to vortex 
c          fn = dot(env(1,i),fi)
c          dcpi(i) = dcpi(i) + 2.0*fn / (rho*dxv(i)*wstrip(js))
cc
c          do n = 1, numax
c            fn_u = dot(env(1,i),fi_u(1,n))
c            dcpi_u(i,n) = dcpi_u(i,n) + 2.0*fn_u / (rho*dxv(i)*wstrip(js))
c          enddo
cc
c          do n = 1, ncontrol
c            fn_d = dot(env(1,i),fi_d(1,n)) + dot(env_d(1,i,n),fi)
c            dcpi_d(i,n) = dcpi_d(i,n) + 2.0*fn_d / (rho*dxv(i)*wstrip(js))
c          enddo
cc
c          do n = 1, nvarjet
c            fn_j = dot(env(1,i),fi_j(1,n)) + dot(env_j(1,i,n),fi)
c            dcpi_j(i,n) = dcpi_j(i,n) + 2.0*fn_j / (rho*dxv(i)*wstrip(js))
c          enddo

c          do n = 1, ndesign
c            fn_g = dot(env(1,i),fi_g(1,n)) + dot(env_g(1,i,n),fi)
c            dcpi_g(i,n) = dcpi_g(i,n) + 2.0*fn_g / (rho*dxv(i)*wstrip(js))
c          enddo
c------------------------------------------------------------------------
c-------- accumulate strip forces and moments
          do k = 1, 3
            ic = icrs(k)
            jc = jcrs(k)

            fsi(k,js) = fsi(k,js) + fi(k)
            msi(k,js) = msi(k,js) + ( r(ic)*fi(jc)
     &                               -r(jc)*fi(ic))
            do n = 1, numax
              fsi_u(k,js,n) = fsi_u(k,js,n) +  fi_u(k,n)
              msi_u(k,js,n) = msi_u(k,js,n) + (r(ic)*fi_u(jc,n)
     &                                        -r(jc)*fi_u(ic,n))
            enddo
            do n = 1, ncontrol
              fsi_d(k,js,n) = fsi_d(k,js,n) +  fi_d(k,n)
              msi_d(k,js,n) = msi_d(k,js,n) + (r(ic)*fi_d(jc,n)
     &                                        -r(jc)*fi_d(ic,n))
            enddo
            do n = 1, nvarjet
              fsi_j(k,js,n) = fsi_j(k,js,n) +  fi_j(k,n)
              msi_j(k,js,n) = msi_j(k,js,n) + (r(ic)*fi_j(jc,n)
     &                                        -r(jc)*fi_j(ic,n))
            enddo
            do n = 1, ndesign
              fsi_g(k,js,n) = fsi_g(k,js,n) +  fi_g(k,n)
              msi_g(k,js,n) = msi_g(k,js,n) + (r(ic)*fi_g(jc,n)
     &                                        -r(jc)*fi_g(ic,n))
            enddo
          enddo

cc-------- hinge moments
c          do l=1, ncontrol
c            rh(1) = rv(1,i) - phinge(1,js,l)
c            rh(2) = rv(2,i) - phinge(2,js,l)
c            rh(3) = rv(3,i) - phinge(3,js,l)
cc
c            dfac = dcontrol(i,l) / (sref * cref)
cc
c            f(1) = fi(1)
c            f(2) = fi(2)
c            f(3) = fi(3)
c            call cross(rh,f,m)
c            chinge(l) = chinge(l) + dot(m,vhinge(1,js,l))*dfac
cc
c            do n = 1, numax
c              f(1) = fi_u(1,n)
c              f(2) = fi_u(2,n)
c              f(3) = fi_u(3,n)
c              call cross(rh,f,m)
c              chinge_u(l,n) = chinge_u(l,n) + dot(m,vhinge(1,js,l))*dfac
c            enddo
c            do n = 1, ncontrol
c              f(1) = fi_d(1,n)
c              f(2) = fi_d(2,n)
c              f(3) = fi_d(3,n)
c              call cross(rh,f,m)
c              chinge_d(l,n) = chinge_d(l,n) + dot(m,vhinge(1,js,l))*dfac
c            enddo
c            do n = 1, nvarjet
c              f(1) = fi_j(1,n)
c              f(2) = fi_j(2,n)
c              f(3) = fi_j(3,n)
c              call cross(rh,f,m)
c              chinge_j(l,n) = chinge_j(l,n) + dot(m,vhinge(1,js,l))*dfac
c            enddo
c            do n = 1, ndesign
c              f(1) = fi_g(1,n)
c              f(2) = fi_g(2,n)
c              f(3) = fi_g(3,n)
c              call cross(rh,f,m)
c              chinge_g(l,n) = chinge_g(l,n) + dot(m,vhinge(1,js,l))*dfac
c            enddo
c          enddo

   71   continue
   72   continue
 80     continue

c
c*******************************************************************
c--- drag terms due to viscous effects
c    drag forces are assumed to be characterized by velocity at the c/4 
c    point and are assumed to act thru the same point. cd is defined by 
c    user-specified cd(cl) polar.  drag comes from function lookup on 
c    section polar drag using local lift coefficient.  

        if(lvisc.and.lviscstrp(js)) then
c------- viscous force is assumed to be applied at c/4 point
         r(1) = rc4(1) - xyzref(1)
         r(2) = rc4(2) - xyzref(2)
         r(3) = rc4(3) - xyzref(3)

         rrot(1) = rc4(1) - xyzref(1)
         rrot(2) = rc4(2) - xyzref(2)
         rrot(3) = rc4(3) - xyzref(3)

c------- set total effective onset velocity = freestream + rotation
         do k = 1, 3
           ic = icrs(k)
           jc = jcrs(k)
           wxr(k) = wbar(ic)*rrot(jc)
     &            - wbar(jc)*rrot(ic)
           wxr_w(k,ic) =  rrot(jc)
           wxr_w(k,jc) = -rrot(ic)
           wxr_w(k,k) = 0.
         enddo

         do k = 1, 3
           vtot(k) = vbar(k) - wxr(k)
           vtot_u(k,1) = 0.
           vtot_u(k,2) = 0.
           vtot_u(k,3) = 0.
           vtot_u(k,k) = 1.0
           vtot_u(k,4) = -wxr_w(k,1)
           vtot_u(k,5) = -wxr_w(k,2)
           vtot_u(k,6) = -wxr_w(k,3)
         enddo

         vtota = sqrt(vtot(1)**2 + vtot(2)**2 + vtot(3)**2)
         do n = 1, numax
           vtota_u(n) = ( vtot(1)*vtot_u(1,n)
     &                  + vtot(2)*vtot_u(2,n)
     &                  + vtot(3)*vtot_u(3,n) ) / vtota
         enddo

c------- generate cd from stored function using strip cl as parameter

c---- for clv to use only inviscid lift (corresponding to CLair)
         rho = rhos(js)
         clv = ( ulift(1)*fsi(1,js)
     &         + ulift(2)*fsi(2,js)
     &         + ulift(3)*fsi(3,js) )*2.0 / (rho*sstrip)

         do n = 1, numax
           clv_u(n) = ( ulift(1)*fsi_u(1,js,n)
     &                + ulift(2)*fsi_u(2,js,n)
     &                + ulift(3)*fsi_u(3,js,n)
     &                + ulift_u(1,n)*fsi(1,js)
     &                + ulift_u(2,n)*fsi(2,js)
     &                + ulift_u(3,n)*fsi(3,js) )*2.0 / (rho*sstrip)
         enddo
         do n = 1, ncontrol
           clv_d(n) = ( ulift(1)*fsi_d(1,js,n)
     &                + ulift(2)*fsi_d(2,js,n)
     &                + ulift(3)*fsi_d(3,js,n) )*2.0 / (rho*sstrip)
         enddo
         do n = 1, nvarjet
           clv_j(n) = ( ulift(1)*fsi_j(1,js,n)
     &                + ulift(2)*fsi_j(2,js,n)
     &                + ulift(3)*fsi_j(3,js,n) )*2.0 / (rho*sstrip)
         enddo
         do n = 1, ndesign
           clv_g(n) = ( ulift(1)*fsi_g(1,js,n)
     &                + ulift(2)*fsi_g(2,js,n)
     &                + ulift(3)*fsi_g(3,js,n) )*2.0 / (rho*sstrip)
         enddo
        
c===================================================================
c---- for clv to use total lift (inviscid+jet) add jet contributions  
         if(lclvjet) then
          clv = clv +
     &          ( ulift(1)*fsj(1,js)
     &          + ulift(2)*fsj(2,js)
     &          + ulift(3)*fsj(3,js) )*2.0 / (rho*sstrip)
          do n = 1, numax
             clv_u(n) = clv_u(n) +
     &                 ( ulift(1)*fsj_u(1,js,n)
     &                 + ulift(2)*fsj_u(2,js,n)
     &                 + ulift(3)*fsj_u(3,js,n)
     &                 + ulift_u(1,n)*fsj(1,js)
     &                 + ulift_u(2,n)*fsj(2,js)
     &                 + ulift_u(3,n)*fsj(3,js) )*2.0 / (rho*sstrip)
          enddo
          do n = 1, ncontrol
            clv_d(n) = clv_d(n) +
     &                 ( ulift(1)*fsj_d(1,js,n)
     &                 + ulift(2)*fsj_d(2,js,n)
     &                 + ulift(3)*fsj_d(3,js,n) )*2.0 / (rho*sstrip)
          enddo
          do n = 1, nvarjet
             clv_j(n) = clv_j(n) +
     &                 ( ulift(1)*fsj_j(1,js,n)
     &                 + ulift(2)*fsj_j(2,js,n)
     &                 + ulift(3)*fsj_j(3,js,n) )*2.0 / (rho*sstrip)
          enddo
          do n = 1, ndesign
             clv_g(n) = clv_g(n) +
     &                 ( ulift(1)*fsj_g(1,js,n)
     &                 + ulift(2)*fsj_g(2,js,n)
     &                 + ulift(3)*fsj_g(3,js,n) )*2.0 / (rho*sstrip)
          enddo
         endif
c===================================================================
c     
         call cdcl(clv,cdv,cdv_clv,cdcls(1,js))

         cdv_strp(js) = cdv

         cda     = cdv    *sstrip
         cda_clv = cdv_clv*sstrip

c------- viscous strip force
         do k = 1, 3
           fsv(k,js) = 0.5*rho*vtot(k)*vtota * cda
           do n = 1, numax
             fsv_u(k,js,n) = 0.5*rho*( vtot_u(k,n)*vtota
     &                               + vtot(k)    *vtota_u(n) )*cda
     &                     + 0.5*rho*vtot(k)*vtota * cda_clv*clv_u(n)
           enddo
           do n = 1, ncontrol
             fsv_d(k,js,n) =  0.5*rho*vtot(k)*vtota * cda_clv*clv_d(n)
           enddo
           do n = 1, nvarjet
             fsv_j(k,js,n) =  0.5*rho*vtot(k)*vtota * cda_clv*clv_j(n)
           enddo
           do n = 1, ndesign
             fsv_g(k,js,n) =  0.5*vtot(k)*vtota * cda_clv*clv_g(n)
           enddo
         enddo

c------- viscous strip moment
         do k = 1, 3
           ic = icrs(k)
           jc = jcrs(k)
           msv(k,js) = r(ic)*fsv(jc,js)
     &               - r(jc)*fsv(ic,js)
           do n = 1, numax
             msv_u(k,js,n) = r(ic)*fsv_u(jc,js,n)
     &                     - r(jc)*fsv_u(ic,js,n)
           enddo
           do n = 1, ncontrol
             msv_d(k,js,n) = r(ic)*fsv_d(jc,js,n)
     &                     - r(jc)*fsv_d(ic,js,n)
           enddo
           do n = 1, nvarjet
             msv_j(k,js,n) = r(ic)*fsv_j(jc,js,n)
     &                     - r(jc)*fsv_j(ic,js,n)
           enddo
           do n = 1, ndesign
             msv_g(k,js,n) = r(ic)*fsv_g(jc,js,n)
     &                     - r(jc)*fsv_g(ic,js,n)
           enddo
         enddo

        else
         cdv_strp(js) =  0.

        endif        

c------ thrust force on propulsor strip
        do k = 1, 3
          fp(k) = taxstrp(k,js)*(djp-dqp)*wstrip(js)
          do n = 1, numax
            fp_u(k,n) = taxstrp(k,js)*(djp_u(n)-dqp_u(n))*wstrip(js)
          enddo
          do n = 1, nvarjet
            fp_j(k,n) = taxstrp(k,js)*(djp_j(n)-dqp_j(n))*wstrip(js)
          enddo
        enddo

c------ thrust moment arm
        dy = rle2(2,js) - rle1(2,js)
        dz = rle2(3,js) - rle1(3,js)
        ds = sqrt(dy**2 + dz**2)
        eny = -dz/ds
        enz =  dy/ds
        r(1) = rle(1,js) + dxdstrp(js)     - xyzref(1)
        r(2) = rle(2,js) + dndstrp(js)*eny - xyzref(2)
        r(3) = rle(3,js) + dndstrp(js)*enz - xyzref(3)

        do k = 1, 3
          ic = icrs(k)
          jc = jcrs(k)
          fsp(k,js) = fp(k)
          msp(k,js) = r(ic)*fp(jc)
     &              - r(jc)*fp(ic)
          do n = 1, numax
            fsp_u(k,js,n) = fp_u(k,n)
            msp_u(k,js,n) = r(ic)*fp_u(jc,n)
     &                    - r(jc)*fp_u(ic,n)
          enddo
          do n = 1, nvarjet
            fsp_j(k,js,n) = fp_j(k,n)
            msp_j(k,js,n) = r(ic)*fp_j(jc,n)
     &                    - r(jc)*fp_j(ic,n)
          enddo
        enddo

c------ store strip mass, momentum, energy defects
        dqs(js) = dqp*wstrip(js)
        djs(js) = djp*wstrip(js)
        des(js) = dep*wstrip(js)
        do n = 1, numax
          dqs_u(js,n) = dqp_u(n)*wstrip(js)
          djs_u(js,n) = djp_u(n)*wstrip(js)
          des_u(js,n) = dep_u(n)*wstrip(js)
        enddo
        do n = 1, nvarjet
          dqs_j(js,n) = dqp_j(n)*wstrip(js)
          djs_j(js,n) = djp_j(n)*wstrip(js)
          des_j(js,n) = dep_j(n)*wstrip(js)
        enddo
  100 continue

c---- accumulate surface forces from strip forces
      do 200 is = 1, nsurf
c------ clear for accumulation
        do k = 1, 3
          fni(k,is) = 0.
          mni(k,is) = 0.
          fnj(k,is) = 0.
          mnj(k,is) = 0.
          fnp(k,is) = 0.
          mnp(k,is) = 0.
          fnv(k,is) = 0.
          mnv(k,is) = 0.
          do n = 1, numax
            fni_u(k,is,n) = 0.
            mni_u(k,is,n) = 0.
            fnj_u(k,is,n) = 0.
            mnj_u(k,is,n) = 0.
            fnp_u(k,is,n) = 0.
            mnp_u(k,is,n) = 0.
            fnv_u(k,is,n) = 0.
            mnv_u(k,is,n) = 0.
          enddo
          do n = 1, ncontrol
            fni_d(k,is,n) = 0.
            mni_d(k,is,n) = 0.
            fnj_d(k,is,n) = 0.
            mnj_d(k,is,n) = 0.
            fnp_d(k,is,n) = 0.
            mnp_d(k,is,n) = 0.
            fnv_d(k,is,n) = 0.
            mnv_d(k,is,n) = 0.
          enddo
          do n = 1, nvarjet
            fni_j(k,is,n) = 0.
            mni_j(k,is,n) = 0.
            fnj_j(k,is,n) = 0.
            mnj_j(k,is,n) = 0.
            fnp_j(k,is,n) = 0.
            mnp_j(k,is,n) = 0.
            fnv_j(k,is,n) = 0.
            mnv_j(k,is,n) = 0.
          enddo
          do n = 1, ndesign
            fni_g(k,is,n) = 0.
            mni_g(k,is,n) = 0.
            fnj_g(k,is,n) = 0.
            mnj_g(k,is,n) = 0.
            fnp_g(k,is,n) = 0.
            mnp_g(k,is,n) = 0.
            fnv_g(k,is,n) = 0.
            mnv_g(k,is,n) = 0.
          enddo
        enddo
        dqn(is) = 0.
        djn(is) = 0.
        den(is) = 0.
        do n = 1, numax
          dqn_u(is,n) = 0.
          djn_u(is,n) = 0.
          den_u(is,n) = 0.
        enddo
        do n = 1, nvarjet
          dqn_j(is,n) = 0.
          djn_j(is,n) = 0.
          den_j(is,n) = 0.
        enddo

c------ add up all the strip contributions for this surface 
        do 220 js = jfrst(is), jlast(is)
          do k = 1, 3
            fni(k,is) = fni(k,is) + fsi(k,js)
            mni(k,is) = mni(k,is) + msi(k,js)
            fnj(k,is) = fnj(k,is) + fsj(k,js)
            mnj(k,is) = mnj(k,is) + msj(k,js)
            fnp(k,is) = fnp(k,is) + fsp(k,js)
            mnp(k,is) = mnp(k,is) + msp(k,js)
            fnv(k,is) = fnv(k,is) + fsv(k,js)
            mnv(k,is) = mnv(k,is) + msv(k,js)
            do n = 1, numax
              fni_u(k,is,n) = fni_u(k,is,n) + fsi_u(k,js,n)
              mni_u(k,is,n) = mni_u(k,is,n) + msi_u(k,js,n)
              fnj_u(k,is,n) = fnj_u(k,is,n) + fsj_u(k,js,n)
              mnj_u(k,is,n) = mnj_u(k,is,n) + msj_u(k,js,n)
              fnp_u(k,is,n) = fnp_u(k,is,n) + fsp_u(k,js,n)
              mnp_u(k,is,n) = mnp_u(k,is,n) + msp_u(k,js,n)
              fnv_u(k,is,n) = fnv_u(k,is,n) + fsv_u(k,js,n)
              mnv_u(k,is,n) = mnv_u(k,is,n) + msv_u(k,js,n)
            enddo
            do n = 1, ncontrol
              fni_d(k,is,n) = fni_d(k,is,n) + fsi_d(k,js,n)
              mni_d(k,is,n) = mni_d(k,is,n) + msi_d(k,js,n)
              fnj_d(k,is,n) = fnj_d(k,is,n) + fsj_d(k,js,n)
              mnj_d(k,is,n) = mnj_d(k,is,n) + msj_d(k,js,n)
              fnp_d(k,is,n) = fnp_d(k,is,n) + fsp_d(k,js,n)
              mnp_d(k,is,n) = mnp_d(k,is,n) + msp_d(k,js,n)
              fnv_d(k,is,n) = fnv_d(k,is,n) + fsv_d(k,js,n)
              mnv_d(k,is,n) = mnv_d(k,is,n) + msv_d(k,js,n)
            enddo
            do n = 1, nvarjet
              fni_j(k,is,n) = fni_j(k,is,n) + fsi_j(k,js,n)
              mni_j(k,is,n) = mni_j(k,is,n) + msi_j(k,js,n)
              fnj_j(k,is,n) = fnj_j(k,is,n) + fsj_j(k,js,n)
              mnj_j(k,is,n) = mnj_j(k,is,n) + msj_j(k,js,n)
              fnp_j(k,is,n) = fnp_j(k,is,n) + fsp_j(k,js,n)
              mnp_j(k,is,n) = mnp_j(k,is,n) + msp_j(k,js,n)
              fnv_j(k,is,n) = fnv_j(k,is,n) + fsv_j(k,js,n)
              mnv_j(k,is,n) = mnv_j(k,is,n) + msv_j(k,js,n)
            enddo
            do n = 1, ndesign
              fni_g(k,is,n) = fni_g(k,is,n) + fsi_g(k,js,n)
              mni_g(k,is,n) = mni_g(k,is,n) + msi_g(k,js,n)
              fnj_g(k,is,n) = fnj_g(k,is,n) + fsj_g(k,js,n)
              mnj_g(k,is,n) = mnj_g(k,is,n) + msj_g(k,js,n)
              fnp_g(k,is,n) = fnp_g(k,is,n) + fsp_g(k,js,n)
              mnp_g(k,is,n) = mnp_g(k,is,n) + msp_g(k,js,n)
              fnv_g(k,is,n) = fnv_g(k,is,n) + fsv_g(k,js,n)
              mnv_g(k,is,n) = mnv_g(k,is,n) + msv_g(k,js,n)
            enddo
          enddo

          dqn(is) = dqn(is) + dqs(js) 
          djn(is) = djn(is) + djs(js) 
          den(is) = den(is) + des(js) 
          do n = 1, numax
            dqn_u(is,n) = dqn_u(is,n) + dqs_u(js,n)
            djn_u(is,n) = djn_u(is,n) + djs_u(js,n)
            den_u(is,n) = den_u(is,n) + des_u(js,n)
          enddo
          do n = 1, nvarjet
            dqn_j(is,n) = dqn_j(is,n) + dqs_j(js,n)
            djn_j(is,n) = djn_j(is,n) + djs_j(js,n)
            den_j(is,n) = den_j(is,n) + des_j(js,n)
          enddo
  220   continue
  200 continue

c---- add on surface forces to total forces
      do 300 is = 1, nsurf
        if(.not.lfload(is)) go to 300

        do k = 1, 3
          fti(k) = fti(k) + fni(k,is)
          mti(k) = mti(k) + mni(k,is)
          ftj(k) = ftj(k) + fnj(k,is)
          mtj(k) = mtj(k) + mnj(k,is)
          ftp(k) = ftp(k) + fnp(k,is)
          mtp(k) = mtp(k) + mnp(k,is)
          ftv(k) = ftv(k) + fnv(k,is)
          mtv(k) = mtv(k) + mnv(k,is)
          do n = 1, numax
            fti_u(k,n) = fti_u(k,n) + fni_u(k,is,n)
            mti_u(k,n) = mti_u(k,n) + mni_u(k,is,n)
            ftj_u(k,n) = ftj_u(k,n) + fnj_u(k,is,n)
            mtj_u(k,n) = mtj_u(k,n) + mnj_u(k,is,n)
            ftp_u(k,n) = ftp_u(k,n) + fnp_u(k,is,n)
            mtp_u(k,n) = mtp_u(k,n) + mnp_u(k,is,n)
            ftv_u(k,n) = ftv_u(k,n) + fnv_u(k,is,n)
            mtv_u(k,n) = mtv_u(k,n) + mnv_u(k,is,n)
          enddo
          do n = 1, ncontrol
            fti_d(k,n) = fti_d(k,n) + fni_d(k,is,n)
            mti_d(k,n) = mti_d(k,n) + mni_d(k,is,n)
            ftj_d(k,n) = ftj_d(k,n) + fnj_d(k,is,n)
            mtj_d(k,n) = mtj_d(k,n) + mnj_d(k,is,n)
            ftp_d(k,n) = ftp_d(k,n) + fnp_d(k,is,n)
            mtp_d(k,n) = mtp_d(k,n) + mnp_d(k,is,n)
            ftv_d(k,n) = ftv_d(k,n) + fnv_d(k,is,n)
            mtv_d(k,n) = mtv_d(k,n) + mnv_d(k,is,n)
          enddo
          do n = 1, nvarjet
            fti_j(k,n) = fti_j(k,n) + fni_j(k,is,n)
            mti_j(k,n) = mti_j(k,n) + mni_j(k,is,n)
            ftj_j(k,n) = ftj_j(k,n) + fnj_j(k,is,n)
            mtj_j(k,n) = mtj_j(k,n) + mnj_j(k,is,n)
            ftp_j(k,n) = ftp_j(k,n) + fnp_j(k,is,n)
            mtp_j(k,n) = mtp_j(k,n) + mnp_j(k,is,n)
            ftv_j(k,n) = ftv_j(k,n) + fnv_j(k,is,n)
            mtv_j(k,n) = mtv_j(k,n) + mnv_j(k,is,n)
          enddo
          do n = 1, ndesign
            fti_g(k,n) = fti_g(k,n) + fni_g(k,is,n)
            mti_g(k,n) = mti_g(k,n) + mni_g(k,is,n)
            ftj_g(k,n) = ftj_g(k,n) + fnj_g(k,is,n)
            mtj_g(k,n) = mtj_g(k,n) + mnj_g(k,is,n)
            ftp_g(k,n) = ftp_g(k,n) + fnp_g(k,is,n)
            mtp_g(k,n) = mtp_g(k,n) + mnp_g(k,is,n)
            ftv_g(k,n) = ftv_g(k,n) + fnv_g(k,is,n)
            mtv_g(k,n) = mtv_g(k,n) + mnv_g(k,is,n)
          enddo
        enddo

        dqt = dqt + dqn(is)
        djt = djt + djn(is)
        det = det + den(is)
        do n = 1, numax
          dqt_u(n) = dqt_u(n) + dqn_u(is,n)
          djt_u(n) = djt_u(n) + djn_u(is,n)
          det_u(n) = det_u(n) + den_u(is,n)
        enddo
        do n = 1, nvarjet
          dqt_j(n) = dqt_j(n) + dqn_j(is,n)
          djt_j(n) = djt_j(n) + djn_j(is,n)
          det_j(n) = det_j(n) + den_j(is,n)
        enddo

 300  continue ! next is

      return
      end ! sfforc

