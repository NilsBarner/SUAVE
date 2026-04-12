C***********************************************************************
C    Module:  jsetup.f
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

      subroutine setup
c--------------------------------------------------------------------------
C     Sets up the flow-tangency and jet-curvature equation 
c     residuals and Jacobians.
c
c     Inputs:   Global data in commons, for configuration and current state
c
c     Outputs:  Residual and Jacobian matrices
c       resn(i)
c       aicn(i,j)
c       resn_u(i.)
c       resn_d(i.)
c       resn_j(i.)
c       resn_g(i.)
c--------------------------------------------------------------------------
      use jvl_inc
      include 'jvl.inc'
      real enctot(3)
      real rrot(3), wxr(3), wxr_w(3,3)

      real vinf, vinf_u(numax)
      real dvjf, dvjf_u(numax)
      real dvjet, dvjet_j(njmax)
      real vjet_u(numax), vjet_j(njmax)
      real djp_u(numax), djp_j(njmax)
      real dvcn, dvcn_gam(nvmax)
      real dwcn, dwcn_u(numax)
      real delf, delf_d(ndmax), delf_g(ngmax)
      real ddte, ddte_d(ndmax), ddte_g(ndmax)

c---- indices for forming cross-products
      integer icrs(3), jcrs(3)
      data icrs / 2, 3, 1 / , jcrs / 3, 1, 2 /

c---- set c.p. velocities vc,wc for current gam(j), vbar,wbar
      do i = 1, nvor
        vc(1,i) = 0.
        vc(2,i) = 0.
        vc(3,i) = 0.
        do j = 1, nvor
          vc(1,i) = vc(1,i) + vc_gam(1,i,j)*gam(j)
          vc(2,i) = vc(2,i) + vc_gam(2,i,j)*gam(j)
          vc(3,i) = vc(3,i) + vc_gam(3,i,j)*gam(j)
        enddo

        do k = 1, 3
          wc(k,i) = wc_u(k,i,1)*vbar(1)
     &            + wc_u(k,i,2)*vbar(2)
     &            + wc_u(k,i,3)*vbar(3)
     &            + wc_u(k,i,4)*wbar(1)
     &            + wc_u(k,i,5)*wbar(2)
     &            + wc_u(k,i,6)*wbar(3)
        enddo
      enddo

c---- clear gamma residuals and jacobians
      do i = 1, nvor
        resn(i) = 0.
        do j = 1, nvor
          aicn(i,j) = 0.
        enddo
        do n = 1, numax
          resn_u(i,n) = 0.
        enddo
        do n = 1, ncontrol
          resn_d(i,n) = 0.
        enddo
        do n = 1, nvarjet
          resn_j(i,n) = 0.
        enddo
        do n = 1, ndesign
          resn_g(i,n) = 0.
        enddo
      enddo

      vinf = sqrt(vbar(1)**2 + vbar(2)**2 + vbar(3)**2)
      vinf_u(1) = vbar(1)/vinf
      vinf_u(2) = vbar(2)/vinf
      vinf_u(3) = vbar(3)/vinf
      vinf_u(4) = 0.
      vinf_u(5) = 0.
      vinf_u(6) = 0.

c---- go over all elements
      do 100 js = 1, nstrip

c------ go over surface elements in strip
        do 110 i = ifrsts(js), ilasts(js)
c-------- set up flow tangency equations

c          v.n + V.n - Vn = 0
c
c    v(i)  =  sum_j vc_gam(i,j) * gam(j)    ;  h.v. velocity
c    V(i)  =  sum_k vc_u(i,k) * u(k)        ;  freestream + body velocity
c          =  Vbar - Wbar x r(i)  +  w(i) 
c    Vn    =  Vstrip fstall(cl)             ;  viscous "leakage" velocity
c    cl = 2*gam(i) / (Vstrip c)             ;  strip lift coefficient
c
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

ccc       resn(i) = [vc(i) + wc(i) + Vbar - Wbar x r(i)] . n(i)
ccc
ccc         vc( gam(j) )    = vc_gam(j)*gam(j)
ccc         wc( Vbar,Wbar ) = wc_u(1:3)*Vbar + wc_w(4:6)*Wbar
ccc         n( D )          = n0  +  n_D*D

          resn(i) = resn(i) + (vc(1,i) + wc(1,i))*enctot(1)
     &                      + (vc(2,i) + wc(2,i))*enctot(2)
     &                      + (vc(3,i) + wc(3,i))*enctot(3)
          do j = 1, nvor
c---------- dresn(i)/dgam(j)
            aicn(i,j) = vc_gam(1,i,j)*enctot(1)
     &                + vc_gam(2,i,j)*enctot(2)
     &                + vc_gam(3,i,j)*enctot(3)
          enddo
          do n = 1, numax
            resn_u(i,n) = resn_u(i,n) 
     &                  + wc_u(1,i,n)*enctot(1)
     &                  + wc_u(2,i,n)*enctot(2)
     &                  + wc_u(3,i,n)*enctot(3)
          enddo
          do n = 1, ncontrol
            resn_d(i,n) = resn_d(i,n)
     &                  + (vc(1,i) + wc(1,i))*enc_d(1,i,n)
     &                  + (vc(2,i) + wc(2,i))*enc_d(2,i,n)
     &                  + (vc(3,i) + wc(3,i))*enc_d(3,i,n)
          enddo
          do n = 1, ndesign
            resn_g(i,n) = resn_g(i,n)
     &                  + (vc(1,i) + wc(1,i))*enc_g(1,i,n)
     &                  + (vc(2,i) + wc(2,i))*enc_g(2,i,n)
     &                  + (vc(3,i) + wc(3,i))*enc_g(3,i,n)
          enddo

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - 
          if(lvalbe(i)) then
c--------- body-rotation velocity at c.p.,  Wbar x (rc-rref)
           rrot(1) = rc(1,i) - xyzref(1)
           rrot(2) = rc(2,i) - xyzref(2)
           rrot(3) = rc(3,i) - xyzref(3)
           do k = 1, 3
             ic = icrs(k)
             jc = jcrs(k)
             wxr(k) = wbar(ic)*rrot(jc)
     &              - wbar(jc)*rrot(ic)
             wxr_w(k,ic) =  rrot(jc)
             wxr_w(k,jc) = -rrot(ic)
             wxr_w(k,k) = 0.
           enddo

c--------- add onset-flow Vbar,Wbar contributions
           resn(i) = resn(i) + (vbar(1) - wxr(1))*enctot(1)
     &                       + (vbar(2) - wxr(2))*enctot(2)
     &                       + (vbar(3) - wxr(3))*enctot(3)
           resn_u(i,1) = resn_u(i,1) + enctot(1)
           resn_u(i,2) = resn_u(i,2) + enctot(2)
           resn_u(i,3) = resn_u(i,3) + enctot(3)
           resn_u(i,4) = resn_u(i,4) - wxr_w(1,1)*enctot(1)
     &                               - wxr_w(2,1)*enctot(2)
     &                               - wxr_w(3,1)*enctot(3)
           resn_u(i,5) = resn_u(i,5) - wxr_w(1,2)*enctot(1)
     &                               - wxr_w(2,2)*enctot(2)
     &                               - wxr_w(3,2)*enctot(3)
           resn_u(i,6) = resn_u(i,6) - wxr_w(1,3)*enctot(1)
     &                               - wxr_w(2,3)*enctot(2)
     &                               - wxr_w(3,3)*enctot(3)
           do n = 1, ncontrol
             resn_d(i,n) = resn_d(i,n)
     &            + (vbar(1) - wxr(1))*enc_d(1,i,n)
     &            + (vbar(2) - wxr(2))*enc_d(2,i,n)
     &            + (vbar(3) - wxr(3))*enc_d(3,i,n)
           enddo

           do n = 1, ndesign
             resn_g(i,n) = resn_g(i,n)
     &            + (vbar(1) - wxr(1))*enc_g(1,i,n)
     &            + (vbar(2) - wxr(2))*enc_g(2,i,n)
     &            + (vbar(3) - wxr(3))*enc_g(3,i,n)
           enddo

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - 
          else
c--------- this c.p. does not see alpha, beta, rotation ...
c-         ... only freestream along x
           resn(i)     = resn(i)     + 1.0*enctot(1)
           resn_u(i,1) = resn_u(i,1) +     enctot(1)

           do n = 1, ncontrol
             resn_d(i,n) = resn_d(i,n) + 1.0*enc_d(1,i,n)
           enddo

           do n = 1, ndesign
             resn_g(i,n) = resn_g(i,n) + 1.0*enc_g(1,i,n)
           enddo

          endif
 110    continue ! next strip element

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

c------ freestream-normalized jet momentum excess per span
        djp      = (rhjet*vjet**2 - vinf**2)*hjet
        djp_vjet =  rhjet*vjet*2.           *hjet 
     &           + (rhjet*vjet**2 - vinf**2)*hjet_vjet
        djp_vinf = (              - vinf*2.)*hjet
     &           + (rhjet*vjet**2 - vinf**2)*hjet_vinf
        do n = 1, numax
          djp_u(n) = djp_vjet*vjet_u(n)
     &             + djp_vinf*vinf_u(n)
        enddo
        do n = 1, nvarjet
          djp_j(n) = djp_vjet*vjet_j(n)
        enddo

c        djp = 5.0*deljet(1)
c        djp_j(1) = 5.0

c------ go over upstream jet elements, if any
        do 120 i = ifrstu(js), ilastu(js)
c-------- index "i-1" of element directly upstream (-999 for first element)
          im = ijetm(i)

c-------- set up jet curvature/loading equation
c
c        resn(i)  =  DJp * ((dvc.n + dwc.n)/Vinf + djet) - gam  =  0
c
c   dvc = vc(i) - vc(i-1) = sum_j [ vc_gam(i,j) - vc_gam(i-1,j) ]*gam(j)
c   dwc = wc(i) - wc(i-1) = sum_n [ wc_u(i,n)   - wc_u(i-1,n)   ]*(Vbar,Wbar)
c   DJp(J)  = (Vjet^2 - Vinf^2)*hdisk*(1 + Vinf/Vjet)/2
c   Vjet(J) = Vinf  +  sum_n [ Jn * gain(n) ]
c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c-------- flow curvature as difference of perturbation-velocity directions
          if(i .eq. ifrstu(js)) then
c--------- upstream element is on propulsor
           dvcn = vc(1,i)*enc(1,i)
     &          + vc(2,i)*enc(2,i)
     &          + vc(3,i)*enc(3,i)
           do j = 1, nvor
             dvcn_gam(j) = vc_gam(1,i,j)*enc(1,i)
     &                   + vc_gam(2,i,j)*enc(2,i)
     &                   + vc_gam(3,i,j)*enc(3,i)
           enddo

           dwcn = wc(1,i)*enc(1,i)
     &          + wc(2,i)*enc(2,i)
     &          + wc(3,i)*enc(3,i)
           do n = 1, numax
             dwcn_u(n) = wc_u(1,i,n)*enc(1,i)
     &                 + wc_u(2,i,n)*enc(2,i)
     &                 + wc_u(3,i,n)*enc(3,i)
           enddo

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          else
c--------- on upstream jet
           dvcn = (vc(1,i)-vc(1,im))*enc(1,i)
     &          + (vc(2,i)-vc(2,im))*enc(2,i)
     &          + (vc(3,i)-vc(3,im))*enc(3,i)
           do j = 1, nvor
             dvcn_gam(j) = (vc_gam(1,i,j)-vc_gam(1,im,j))*enc(1,i)
     &                   + (vc_gam(2,i,j)-vc_gam(2,im,j))*enc(2,i)
     &                   + (vc_gam(3,i,j)-vc_gam(3,im,j))*enc(3,i)
           enddo

           dwcn = (wc(1,i)-wc(1,im))*enc(1,i)
     &          + (wc(2,i)-wc(2,im))*enc(2,i)
     &          + (wc(3,i)-wc(3,im))*enc(3,i)
           do n = 1, numax
             dwcn_u(n) = (wc_u(1,i,n)-wc_u(1,im,n))*enc(1,i)
     &                 + (wc_u(2,i,n)-wc_u(2,im,n))*enc(2,i)
     &                 + (wc_u(3,i,n)-wc_u(3,im,n))*enc(3,i)
           enddo
          endif

          resn(i) = djp*(dvcn + dwcn)/vinf - gam(i)
          do j = 1, nvor
            aicn(i,j) = djp*dvcn_gam(j)/vinf
          enddo
          aicn(i,i) = aicn(i,i) - 1.0

          do n = 1, numax
            resn_u(i,n) = djp_u(n)*(dvcn + dwcn     )/vinf
     &                  + djp     *        dwcn_u(n) /vinf
     &                  - djp     *(dvcn + dwcn     )/vinf**2*vinf_u(n)
          enddo

          do n = 1, ncontrol
            resn_d(i,n) = 0.
          enddo

          do n = 1, nvarjet
            resn_j(i,n) = djp_j(n)*(dvcn + dwcn)/vinf
          enddo

          do n = 1, ndesign
            resn_g(i,n) = 0.
          enddo
 120    continue ! next jet element

c------ go over wake jet elements, if any
        do 130 i = ifrstw(js), ilastw(js)

c-------- index "i-1" of element directly upstream
          im = ijetm(i)

c-------- set up jet curvature/loading equation
c
c        resn(i)  =  DJp * (dvc.n + dwc.n + djet)/vinf - gam  =  0
c
c   dvc = vc(i) - vc(i-1) = sum_j [ vc_gam(i,j) - vc_gam(i-1,j) ]*gam(j)
c   dwc = wc(i) - wc(i-1) = sum_n [ wc_u(i,n)   - wc_u(i-1,n)   ]*(Vbar,Wbar)
c   DJp(J)  = (Vjet^2 - Vinf^2)*hdisk*(1 + Vinf/Vjet)/2
c   Vjet(J) = Vinf  +  sum_n [ Jn * gain(n) ]
c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          if(i .eq. ifrstw(js)) then
c--------- upstream element is on surface TE -- set jet deviation angle term

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

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          else
c--------- upstream element is not on surface --- no TE jet deviation term here
           ddte = 0.
           do n = 1, ncontrol
             ddte_d(n) = 0.
           enddo
           do n = 1, ndesign
             ddte_g(n) = 0.
           enddo

          endif
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         
c-------- flow curvature as difference of perturbation-velocity directions
          dvcn = (vc(1,i)-vc(1,im))*enc(1,i)
     &         + (vc(2,i)-vc(2,im))*enc(2,i)
     &         + (vc(3,i)-vc(3,im))*enc(3,i)
          do j = 1, nvor
            dvcn_gam(j) = (vc_gam(1,i,j)-vc_gam(1,im,j))*enc(1,i)
     &                  + (vc_gam(2,i,j)-vc_gam(2,im,j))*enc(2,i)
     &                  + (vc_gam(3,i,j)-vc_gam(3,im,j))*enc(3,i)
          enddo

          dwcn = (wc(1,i)-wc(1,im))*enc(1,i)
     &         + (wc(2,i)-wc(2,im))*enc(2,i)
     &         + (wc(3,i)-wc(3,im))*enc(3,i)
          do n = 1, numax
            dwcn_u(n) = (wc_u(1,i,n)-wc_u(1,im,n))*enc(1,i)
     &                + (wc_u(2,i,n)-wc_u(2,im,n))*enc(2,i)
     &                + (wc_u(3,i,n)-wc_u(3,im,n))*enc(3,i)
          enddo

          resn(i) = djp*((dvcn + dwcn)/vinf + dtr*ddte) - gam(i)
          do j = 1, nvor
            aicn(i,j) = djp*dvcn_gam(j)/vinf
          enddo
          aicn(i,i) = aicn(i,i) - 1.0

          do n = 1, numax
            resn_u(i,n) = djp_u(n)*((dvcn + dwcn     )/vinf + dtr*ddte)
     &                  + djp     *         dwcn_u(n) /vinf
     &                  - djp     * (dvcn + dwcn     )/vinf**2*vinf_u(n)
          enddo

          do n = 1, ncontrol
            resn_d(i,n) = djp*dtr*ddte_d(n)
          enddo

          do n = 1, nvarjet
            resn_j(i,n) = djp_j(n)*((dvcn + dwcn)/vinf + dtr*ddte)
          enddo

          do n = 1, ndesign
            resn_g(i,n) = djp*dtr*ddte_g(n)
          enddo
 130    continue ! next jet element

 100  continue ! next strip

C---- process each surface which does not shed a wake
      do 200 isurf = 1, nsurf
        if(lfwake(isurf)) go to 200

c------ go over strips on this surface
        do 210 js = jfrst(isurf), jlast(isurf)
c-------- clear system row for te control point
          i = ilasts(js)

          resn(i) = 0.
          do n = 1, numax
            resn_u(i,n) = 0.
          enddo
          do n = 1, ncontrol
            resn_d(i,n) = 0.
          enddo
          do n = 1, nvarjet
            resn_j(i,n) = 0.
          enddo
          do n = 1, ndesign
            resn_g(i,n) = 0.
          enddo
          do j = 1, nvor
            aicn(i,j) = 0.
          enddo

c-------- set  sum_strip(gamma) = 0  for this strip
          do j = ifrsts(js), ilasts(js)
            resn(i) = resn(i) + gam(j)
            aicn(i,j) = 1.0
          enddo
 210    continue ! next strip
 200  continue ! next surface


cc...holdover from hpv hydro project for forces near free surface
cc...eliminates excluded vortices from eqns which are below z=zsym 
c     call mungea

      laic = .true.

      return
      end ! setup


      subroutine mungea
c----------------------------------------------------------------
c     Removes hidden vortex equations in aic matrix
c          
c  Input:  aicn(..)  AIC matrix
c
c  Output: aicn(..)  AIC matrix with affected rows cleared, 
c                      with 1 put on diagonal
c----------------------------------------------------------------
      use jvl_inc
      include 'jvl.inc'

      do js = 1, nstrip
        if (lstripoff(js)) then
c------- process all elements in strip, including jet if any
         do i = ifrsts(js), ilasts(js)
           resn(i) = gam(i)
           do j = 1, nvor
             aicn(i,j) = 0.0
           enddo
           aicn(i,i) = 1.0

           do n = 1, numax
             resn_u(i,n) = 0.
           enddo
           do n = 1, ncontrol
             resn_d(i,n) = 0.
           enddo
           do n = 1, nvarjet
             resn_j(i,n) = 0.
           enddo
           do n = 1, ndesign
             resn_g(i,n) = 0.
           enddo
         enddo

         do i = ifrstu(js), ilastu(js)
           resn(i) = gam(i)
           do j = 1, nvor
             aicn(i,j) = 0.0
           enddo
           aicn(i,i) = 1.0

           do n = 1, numax
             resn_u(i,n) = 0.
           enddo
           do n = 1, ncontrol
             resn_d(i,n) = 0.
           enddo
           do n = 1, nvarjet
             resn_j(i,n) = 0.
           enddo
           do n = 1, ndesign
             resn_g(i,n) = 0.
           enddo
         enddo

         do i = ifrstw(js), ilastw(js)
           resn(i) = gam(i)
           do j = 1, nvor
             aicn(i,j) = 0.0
           enddo
           aicn(i,i) = 1.0

           do n = 1, numax
             resn_u(i,n) = 0.
           enddo
           do n = 1, ncontrol
             resn_d(i,n) = 0.
           enddo
           do n = 1, nvarjet
             resn_j(i,n) = 0.
           enddo
           do n = 1, ndesign
             resn_g(i,n) = 0.
           enddo
         enddo

        endif
      enddo

      return
      end ! mungea

