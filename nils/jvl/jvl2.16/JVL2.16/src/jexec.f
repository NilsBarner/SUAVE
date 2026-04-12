C***********************************************************************
C    Module:  jexec.f
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

      subroutine exec(niter1,info,ir)
c---------------------------------------------------
c     Solves JVL system by Newton iteration
c---------------------------------------------------
      use jvl_inc
      include 'jvl.inc'
      real vsys(ivmax,ivmax), vres(ivmax), ddc(ndmax), ddj(njmax)
      real vwork(ivmax)
      integer ivsys(ivmax)
      integer indx(nvmax)
 
      real vbarm(3), dvbar(3), dwbar(3)
      real cax, cax_u(3)
      real sax, sax_u(3)

c---- convergence tolerance
c     data eps / 1.0e-5  /
      data eps / 1.0e-7  /
c     data eps / 1.0e-8  /
c     data eps / 1.0e-11 /

c---- max angle change limit (radians)
      data dmax / 0.1 /

      niter = niter1

      if(lnasa_sa) then
c----- nasa std. stability axes, x fwd, z down
       dirax = -1.0
      else
c----- geometric stability axes, x aft, z up
       dirax =  1.0
      endif

c---- set working variables from varval(.ir) array
      call varini(ir)

      xyzref(1) = parval(ipxcg,ir)
      xyzref(2) = parval(ipycg,ir)
      xyzref(3) = parval(ipzcg,ir)

      xyzmom(1) = xyzref(1)
      xyzmom(2) = xyzref(2)
      xyzmom(3) = xyzref(3)

      cdref = parval(ipcd0,ir)

      mach = parval(ipmach,ir)

      if(mach.ne.amach) then
c----- new mach number invalidates aic matrices and solution
       laic = .false.
       lsol = .false.
      endif

      dvjexp = parval(ipdvjx,ir)

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if(niter.gt.0) then
c----- directly set operating variables if they are known
c-      (should speed up convergence)

       if(icon(ivvelz,ir).eq.icalfa) then
        alfa = conval(icalfa,ir)*dtr
c------ set corresponding vbar(alfa,beta)
        call ab2vbar
       endif

       if(icon(ivvely,ir).eq.icbeta) then
        beta = conval(icbeta,ir)*dtr
c------ set corresponding vbar(alfa,beta)
        call ab2vbar
       endif

       if(icon(ivomgx,ir).eq.icomgx .and.
     &    icon(ivomgz,ir).eq.icomgz      ) then
        pax = dirax*conval(icomgx,ir)*2.0/bref
        rax = dirax*conval(icomgz,ir)*2.0/bref
        if(lsa_rates) then
c------- rates specified in NASA stability-axes, transform to body axes
         cax = cos(alfa)
         sax = sin(alfa)
        else
c------- rates specified in body-axes, no transformation
         cax = 1.0
         sax = 0.0
        endif
        wbar(1) = pax*cax - rax*sax
        wbar(3) = rax*cax + pax*sax
       endif

       if(icon(ivomgx,ir).eq.icomgx .and. .not.lsa_rates) then
        wbar(1) = dirax*conval(icomgx,ir)*2.0/bref
       endif

       if(icon(ivomgy,ir).eq.icomgy) then
        wbar(2) = conval(icomgy,ir)*2./cref
       endif

       if(icon(ivomgz,ir).eq.icomgz .and. .not.lsa_rates) then
        wbar(3) = dirax*conval(icomgz,ir)*2.0/bref
       endif

       do n = 1, ncontrol
         ic = ictot + n
         iv = ivtot + n
         if(icon(iv,ir).eq.ic) delcon(n) = conval(ic,ir)
       enddo

       do n = 1, nvarjet
         ic = ictot + ncontrol + n
         iv = ivtot + ncontrol + n
         if(icon(iv,ir).eq.ic) deljet(n) = conval(ic,ir)
       enddo
      endif

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      if(.not.laic) then
c----- mach number defining prandtl-glauert factor for AIC matrices
       amach = mach
       betm = sqrt(1.0 - amach**2)

       write(*,*) ' Building control-point velocity AIC matrix...'
       call vvor(betm,iysym,ysym,zsym,vrcorec,vrcorew,
     &           nvor,rv1,rv2,lscompv,chordv, izimv,
     &           nvor,rc ,    lscompv,.false.,
     &           vc_gam,nvmax)

       write(*,*) ' Building bound-vortex  velocity AIC matrix...'
       call vvor(betm,iysym,ysym,zsym,vrcorec,vrcorew,
     &           nvor,rv1,rv2,lscompv,chordv, izimv,
     &           nvor,rv ,    lscompv,.true.,
     &           vv_gam,nvmax)

c---------------------------------
       if(nlbody .gt. 0) then
        write(*,*) ' Building body-point velocity AIC matrix...'
        call vvor(betm,iysym,ysym,zsym,vrcorec,vrcorew,
     &            nvor  ,rv1,rv2,lscompv,chordv ,izimv,
     &            nlbody,rl ,    lbcompl,.false.,
     &            vl_gam,nlmax)

        write(*,*) ' Building source+doublet strength AIC matrix...'
        call srdu(betm,xyzref,
     &            nbody,lfrst,llast,nlmax,
     &            ibcent, iysym,
     &            rl1,rl2,rl,drbl,darl,radl,
     &            src_u,dbl_u)

        write(*,*) ' Building source+doublet velocity AIC matrices...'
        nu = 6
        call vsrd(betm,iysym,ysym,zsym,srcore,
     &            nbody,lfrst,llast,nlmax,
     &            rl1,rl2,lbcompl,drbl,darl,radl,iziml,
     &            nu,src_u,dbl_u,
     &            nvor,rc,lscompv,
     &            wc_u,nvmax)
        call vsrd(betm,iysym,ysym,zsym,srcore,
     &            nbody,lfrst,llast,nlmax,
     &            rl1,rl2,lbcompl,drbl,darl,radl,iziml,
     &            nu,src_u,dbl_u,
     &            nvor,rv,lscompv,
     &            wv_u,nvmax)
        call vsrd(betm,iysym,ysym,zsym,srcore,
     &            nbody,lfrst,llast,nlmax,
     &            rl1,rl2,lbcompl,drbl,darl,radl,iziml,
     &            nu,src_u,dbl_u,
     &            nlbody,rl,lbcompl,
     &            wl_u,nlmax)
       endif
c---------------------------------

       laic = .true.
      endif

c---- newton iteration loop
      do 500 iter = 1, niter
c------ set alfa, beta for current Vbar
        call vbar2ab

c------ set up newton system
        call setup
c       write(*,*) ' Factoring gamma AIC matrix...'

ccc     call sgetrf(nvor,nvmax,aicn,nvmax,indx,infos)
        call ludcmp(nvmax,nvor,aicn,indx,work)

c       write(*,*) ' Back-sub global variable influence vectors...'

c------ back-substitute surface-bc residuals and jacobians to get 
c-       gamma Newton changes and its global-variable sensitivities
        do i = 1, nvor
          dgam(i) = -resn(i)
          do n = 1, numax
            gam_u(i,n) = -resn_u(i,n)
          enddo
          do n = 1, ncontrol
            gam_d(i,n) = -resn_d(i,n)
          enddo
          do n = 1, nvarjet
            gam_j(i,n) = -resn_j(i,n)
          enddo
          do n = 1, ndesign
            gam_g(i,n) = -resn_g(i,n)
          enddo
        enddo

ccc     call sgetrs('N',nvor,1,aicn,nvmax,indx,dgam,nvmax,infos)
        call baksub(nvmax,nvor,aicn,indx,dgam)

        do n = 1, numax
ccc       call sgetrs('N',nvor,1,aicn,nvmax,indx,gam_u(1,n),nvmax,infos)
          call baksub(nvmax,nvor,aicn,indx,gam_u(1,n))
        enddo

        do n = 1, ncontrol
ccc       call sgetrs('N',nvor,1,aicn,nvmax,indx,gam_d(1,n),nvmax,infos)
          call baksub(nvmax,nvor,aicn,indx,gam_d(1,n))
        enddo

        do n = 1, nvarjet
ccc       call sgetrs('N',nvor,1,aicn,nvmax,indx,gam_j(1,n),nvmax,infos)
          call baksub(nvmax,nvor,aicn,indx,gam_j(1,n))
        enddo

        do n = 1, ndesign
ccc       call sgetrs('N',nvor,1,aicn,nvmax,indx,gam_g(1,n),nvmax,infos)
          call baksub(nvmax,nvor,aicn,indx,gam_g(1,n))
        enddo

c        write(*,*) ' Computing gamma, forces...'

c------ Partial update of h.v. circulations, 
c-      not accounting for global variable changes (determined next)
c-      This improves convergence of nonlinear cases with jet effects
        do i = 1, nvor
          gam(i) = gam(i) + dgam(i)
        enddo

c------ set h.v., c.p., body velocities vc,vv,vl for current gam(j),
c-      and their u,d,j,g sensitivities
        call vvsum
        call vcsum
        call vlsum
        
c------ set h.v. and c.p. velocities wv,wc,wl for current Vbar,Wbar
        call wvsum
        call wcsum
        call wlsum
        
c------ compute forces and moments, 
c-       and their global-variable u,d,j,g sensitivities
c      print *, 'calling aero'
        call aero
        
        if(lsa_rates) then
c------- rates specified in NASA stability-axes, transform to body axes
c-       sax = sin(alfa)
c-       cax = cos(alfa)
         v13sq = vbar(1)**2 + vbar(3)**2
         v13 = sqrt(v13sq)
         sax = vbar(3)/v13
         cax = vbar(1)/v13

         sax_u(1) = -sax*vbar(1)/v13sq
         sax_u(2) = 0.
         sax_u(3) = -sax*vbar(3)/v13sq + 1.0/v13

         cax_u(1) = -cax*vbar(1)/v13sq + 1.0/v13
         cax_u(2) = 0.
         cax_u(3) = -cax*vbar(3)/v13sq

        else
c------- rates specified in body-axes, no transformation
         cax = 1.0
         sax = 0.0
         cax_u(1) = 0.
         cax_u(2) = 0.
         cax_u(3) = 0.
         sax_u(1) = 0.
         sax_u(2) = 0.
         sax_u(3) = 0.

        endif

c------ clear small global-variable Newton matrix
        do k = 1, nvtot
          vres(k) = 0.
          do l = 1, nvtot
            vsys(k,l) = 0.
          enddo
        enddo

c        write(*,*) ' Setting up global variable system...'

c------ set up Newton system for all global variables
c-       by imposing the specified constraint for each variable

c------------------------------------
c------ |Vbar| = 1 is always enforced, in effect mainly drives vbar(1)
        iv = 1
        vres(iv) = vbar(1)**2 + vbar(2)**2 + vbar(3)**2  - 1.0
        vsys(iv,ivvelx) = 2.0*vbar(1)
        vsys(iv,ivvely) = 2.0*vbar(2)
        vsys(iv,ivvelz) = 2.0*vbar(3)

        do 100 iv = 2, nvtot
c-------- index of constraint for this variable
          ic = icon(iv,ir)

c------------------------------------
          if    (ic.eq.icalfa) then
c--------- specified alpha
           vres(iv) = alfa - conval(ic,ir)*dtr
           vsys(iv,ivvelx) = alfa_u(1)
           vsys(iv,ivvely) = alfa_u(2)
           vsys(iv,ivvelz) = alfa_u(3)

c------------------------------------
          elseif(ic.eq.icbeta) then
c--------- specified beta
           vres(iv) = beta - conval(ic,ir)*dtr
           vsys(iv,ivvelx) = beta_u(1)
           vsys(iv,ivvely) = beta_u(2)
           vsys(iv,ivvelz) = beta_u(3)

c------------------------------------
          elseif(ic.eq.icomgx) then
c--------- specified roll rate
           vres(iv) = wbar(1)*cax + wbar(3)*sax
     &              - dirax*conval(ic,ir)*2.0/bref
           vsys(iv,ivomgx) = cax
           vsys(iv,ivomgz) = sax
           vsys(iv,ivvelx) = wbar(1)*cax_u(1) + wbar(3)*sax_u(1)
           vsys(iv,ivvely) = wbar(1)*cax_u(2) + wbar(3)*sax_u(2)
           vsys(iv,ivvelz) = wbar(1)*cax_u(3) + wbar(3)*sax_u(3)

c------------------------------------
          elseif(ic.eq.icomgy) then
c--------- specified pitch rate
           vres(iv) = wbar(2)  - conval(ic,ir)*2.0/cref
           vsys(iv,ivomgy) = 1.0

c------------------------------------
          elseif(ic.eq.icomgz) then
c--------- specified yaw rate
           vres(iv) = wbar(3)*cax - wbar(1)*sax
     &              - dirax*conval(ic,ir)*2.0/bref
           vsys(iv,ivomgx) = -sax
           vsys(iv,ivomgz) =  cax
           vsys(iv,ivvelx) = wbar(3)*cax_u(1) - wbar(1)*sax_u(1)
           vsys(iv,ivvely) = wbar(3)*cax_u(2) - wbar(1)*sax_u(2)
           vsys(iv,ivvelz) = wbar(3)*cax_u(3) - wbar(1)*sax_u(3)

c------------------------------------
          elseif(ic.eq.iccl  ) then
c--------- specified CL
           vres(iv) = ltot - conval(ic,ir)*0.5*sref
           vsys(iv,ivvelx) = ltot_u(1)
           vsys(iv,ivvely) = ltot_u(2)
           vsys(iv,ivvelz) = ltot_u(3)
           vsys(iv,ivomgx) = ltot_u(4)
           vsys(iv,ivomgy) = ltot_u(5)
           vsys(iv,ivomgz) = ltot_u(6)

           do n = 1, ncontrol
             nv = ivtot + n
             vsys(iv,nv) = ltot_d(n)
           enddo

           do n = 1, nvarjet
             nv = ivtot + ncontrol + n
             vsys(iv,nv) = ltot_j(n)
           enddo

c           if(iter.le.2) then
c------------ fix alpha
c             do l = 1, nvtot
c                vsys(iv,l) = 0.
c              enddo
c              vres(iv) = 0.
c              vsys(iv,ivvelx) = alfa_u(1)
c              vsys(iv,ivvely) = alfa_u(2)
c              vsys(iv,ivvelz) = alfa_u(3)
c           endif
c------------------------------------
          elseif(ic.eq.iccy  ) then
c--------- specified CY
           vres(iv) = ytot - conval(ic,ir)*0.5*sref
           vsys(iv,ivvelx) = ytot_u(1)
           vsys(iv,ivvely) = ytot_u(2)
           vsys(iv,ivvelz) = ytot_u(3)
           vsys(iv,ivomgx) = ytot_u(4)
           vsys(iv,ivomgy) = ytot_u(5)
           vsys(iv,ivomgz) = ytot_u(6)

           do n = 1, ncontrol
             nv = ivtot + n
             vsys(iv,nv) = ytot_d(n)
           enddo

           do n = 1, nvarjet
             nv = ivtot + ncontrol + n
             vsys(iv,nv) = ytot_j(n)
           enddo

c------------------------------------
          elseif(ic.eq.icmomx) then
c--------- specified roll moment
           vres(iv) = mt(1)*cax + mt(3)*sax
     &              - dirax*conval(ic,ir)*0.5*sref*bref
           vsys(iv,ivvelx) = mt_u(1,1)*cax  + mt_u(3,1)*sax
     &                     + mt(1)*cax_u(1) + mt(3)*sax_u(1)
           vsys(iv,ivvely) = mt_u(1,2)*cax  + mt_u(3,2)*sax
     &                     + mt(1)*cax_u(2) + mt(3)*sax_u(2)
           vsys(iv,ivvelz) = mt_u(1,3)*cax  + mt_u(3,3)*sax
     &                     + mt(1)*cax_u(3) + mt(3)*sax_u(3)
           vsys(iv,ivomgx) = mt_u(1,4)*cax + mt_u(3,4)*sax
           vsys(iv,ivomgy) = mt_u(1,5)*cax + mt_u(3,5)*sax
           vsys(iv,ivomgz) = mt_u(1,6)*cax + mt_u(3,6)*sax

           do n = 1, ncontrol
             nv = ivtot + n
             vsys(iv,nv) = mt_d(1,n)*cax + mt_d(3,n)*sax
           enddo

           do n = 1, nvarjet
             nv = ivtot + ncontrol + n
             vsys(iv,nv) = mt_j(1,n)*cax + mt_j(3,n)*sax
           enddo

c------------------------------------
          elseif(ic.eq.icmomy) then
c--------- specified pitching moment
           vres(iv) = mt(2)
     &              - conval(ic,ir)*0.5*sref*cref
           vsys(iv,ivvelx) = mt_u(2,1)
           vsys(iv,ivvely) = mt_u(2,2)
           vsys(iv,ivvelz) = mt_u(2,3)
           vsys(iv,ivomgx) = mt_u(2,4)
           vsys(iv,ivomgy) = mt_u(2,5)
           vsys(iv,ivomgz) = mt_u(2,6)

           do n = 1, ncontrol
             nv = ivtot + n
             vsys(iv,nv) = mt_d(2,n)
           enddo

           do n = 1, nvarjet
             nv = ivtot + ncontrol + n
             vsys(iv,nv) = mt_j(2,n)
           enddo

c------------------------------------
          elseif(ic.eq.icmomz) then
c--------- specified yaw moment
           vres(iv) = mt(3)*cax - mt(1)*sax
     &              - dirax*conval(ic,ir)*0.5*sref*bref
           vsys(iv,ivvelx) = mt_u(3,1)*cax  - mt_u(1,1)*sax
     &                     + mt(3)*cax_u(1) - mt(1)*sax_u(1)
           vsys(iv,ivvely) = mt_u(3,2)*cax  - mt_u(1,2)*sax
     &                     + mt(3)*cax_u(2) - mt(1)*sax_u(2)
           vsys(iv,ivvelz) = mt_u(3,3)*cax  - mt_u(1,3)*sax
     &                     + mt(3)*cax_u(3) - mt(1)*sax_u(3)
           vsys(iv,ivomgx) = mt_u(3,4)*cax - mt_u(1,4)*sax
           vsys(iv,ivomgy) = mt_u(3,5)*cax - mt_u(1,5)*sax
           vsys(iv,ivomgz) = mt_u(3,6)*cax - mt_u(1,6)*sax

           do n = 1, ncontrol
             nv = ivtot + n
             vsys(iv,nv) = mt_d(3,n)*cax - mt_d(1,n)*sax
           enddo

           do n = 1, nvarjet
             nv = ivtot + ncontrol + n
             vsys(iv,nv) = mt_j(3,n)*cax - mt_j(1,n)*sax
           enddo

c------------------------------------
          elseif(ic.eq.icjetj) then
c--------- specified Delta CJ
           vres(iv) = djt - conval(ic,ir)*0.5*sref
           vsys(iv,ivvelx) = djt_u(1)
           vsys(iv,ivvely) = djt_u(2)
           vsys(iv,ivvelz) = djt_u(3)
           vsys(iv,ivomgx) = djt_u(4)
           vsys(iv,ivomgy) = djt_u(5)
           vsys(iv,ivomgz) = djt_u(6)
           do n = 1, nvarjet
             nv = ivtot + ncontrol + n
             vsys(iv,nv) = djt_j(n)
           enddo

c------------------------------------
          elseif(ic.eq.icjett) then
c--------- specified CT
           n = 1
           vres(iv) = djt - dqt - conval(ic,ir)*0.5*sref
           vsys(iv,ivvelx) = djt_u(1) - dqt_u(1)
           vsys(iv,ivvely) = djt_u(2) - dqt_u(2)
           vsys(iv,ivvelz) = djt_u(3) - dqt_u(3)
           vsys(iv,ivomgx) = djt_u(4) - dqt_u(4)
           vsys(iv,ivomgy) = djt_u(5) - dqt_u(5)
           vsys(iv,ivomgz) = djt_u(6) - dqt_u(6)
           do n = 1, nvarjet
             nv = ivtot + ncontrol + n
             vsys(iv,nv) = djt_j(n) - dqt_j(n)
           enddo

c------------------------------------
          elseif(ic.eq.icjetp) then
c--------- specified CP
           vres(iv) = det - dqt - conval(ic,ir)*sref
           vsys(iv,ivvelx) = det_u(1) - dqt_u(1)
           vsys(iv,ivvely) = det_u(2) - dqt_u(2)
           vsys(iv,ivvelz) = det_u(3) - dqt_u(3)
           vsys(iv,ivomgx) = det_u(4) - dqt_u(4)
           vsys(iv,ivomgy) = det_u(5) - dqt_u(5)
           vsys(iv,ivomgz) = det_u(6) - dqt_u(6)
           do n = 1, nvarjet
             nv = ivtot + ncontrol + n
             vsys(iv,nv) = det_j(n) - dqt_j(n)
           enddo

c------------------------------------
          else
c--------- direct control variable constraints
           do n = 1, ncontrol
             iccon = ictot + n
             ivcon = ivtot + n
             if(ic.eq.iccon) then
              vres(iv) = delcon(n) - conval(iccon,ir)
              vsys(iv,ivcon) = 1.0
              go to 100
             endif
           enddo

c--------- direct jet variable constraints
           do n = 1, nvarjet
             iccon = ictot + ncontrol + n
             ivcon = ivtot + ncontrol + n
             if(ic.eq.iccon) then
              vres(iv) = deljet(n) - conval(iccon,ir)
              vsys(iv,ivcon) = 1.0
              go to 100
             endif
           enddo

           write(*,*) '? illegal constraint index: ', ic
          endif

 100    continue

c      do k = 1, nvtot
c        write(*,3344) (vsys(k,l), l=1, nvtot), vres(k)
c 3344   format(1x,20f9.5)
c      enddo

c------ LU-factor, and back-substitute residual rhs
ccc     call sgetrf(nvtot,nvtot,vsys,ivmax,ivsys,infos)
c      print *, 'calling ludcmp'
        call ludcmp(ivmax,nvtot,vsys,ivsys,vwork)

ccc     call sgetrs('N',nvtot,1,vsys,ivmax,ivsys,vres,ivmax,infos)
c      print *, 'calling baksub'
        call baksub(ivmax,nvtot,vsys,ivsys,vres)

c------ set Newton deltas
        dvbar(1) = -vres(ivvelx)
        dvbar(2) = -vres(ivvely)
        dvbar(3) = -vres(ivvelz)
        dwbar(1) = -vres(ivomgx)
        dwbar(2) = -vres(ivomgy)
        dwbar(3) = -vres(ivomgz)
        do n = 1, ncontrol
          iv = ivtot + n
          ddc(n) = -vres(iv)
        enddo
        do n = 1, nvarjet
          iv = ivtot + ncontrol + n
          ddj(n) = -vres(iv)
        enddo

        dal = alfa_u(1)*dvbar(1)
     &      + alfa_u(2)*dvbar(2)
     &      + alfa_u(3)*dvbar(3)
        dbe = beta_u(1)*dvbar(1)
     &      + beta_u(2)*dvbar(2)
     &      + beta_u(3)*dvbar(3)

c------ changes to gamma forced by global variable changes
        do i = 1, nvor
          dgam(i) = gam_u(i,1)*dvbar(1)
     &            + gam_u(i,2)*dvbar(2)
     &            + gam_u(i,3)*dvbar(3)
     &            + gam_u(i,4)*dwbar(1)
     &            + gam_u(i,5)*dwbar(2)
     &            + gam_u(i,6)*dwbar(3)
          do n = 1, ncontrol
            dgam(i) = dgam(i) + gam_d(i,n)*ddc(n)
          enddo
          do n = 1, nvarjet
            dgam(i) = dgam(i) + gam_j(i,n)*ddj(n)
          enddo
        enddo

c------ find max change to an overall strip circulation, for convergence test
        dgsmax = 0.
        do j = 1, nstrip
          dgams = 0.
          do i = ifrsts(j), ilasts(j)
            dgams = dgams + dgam(i)
          enddo
          do i = ifrstu(j), ilastu(j)
            dgams = dgams + dgam(i)
          enddo
          do i = ifrstw(j), ilastw(j)
            dgams = dgams + dgam(i)
          enddo
          if(abs(dgams) .gt. abs(dgsmax)) then
            dgsmax = dgams
          endif
        enddo

c------ limits on angles, rates, control deflections
        dmaxa = dmax
        dmaxr = 5.0*dmax
        dmaxd = 2.0*dmax

        rlx = 1.0

        if(rlx*abs(dal) .gt. dmaxa) then
         rlx = dmaxa/abs(dal)
        endif
        if(rlx*abs(dbe) .gt. dmaxa) then
         rlx = dmaxa/abs(dbe)
        endif
        if(rlx*abs(dwbar(1)) .gt. dmaxr) then
         rlx = dmaxr/abs(dwbar(1))
        endif
        if(rlx*abs(dwbar(2)) .gt. dmaxr) then
         rlx = dmaxr/abs(dwbar(2))
        endif
        if(rlx*abs(dwbar(3)) .gt. dmaxr) then
         rlx = dmaxr/abs(dwbar(3))
        endif
        do n = 1, ncontrol
          if(rlx*abs(ddc(n)*gconmax(n)) .gt. dmaxd) then
           rlx = dmaxd/abs(ddc(n)*gconmax(n))
          endif
        enddo

        if(info .ge. 1) then
c------- display Newton deltas
         if(iter.eq.1) then
          write(*,*)
          write(*,1902)
     &            'iter','  rlx ',
     &            ' d(Gams/c) ',
     &            ' d(alpha)  ',
     &            ' d(beta)   ',
     &            ' d(pb/2V)  ',
     &            ' d(qc/2V)  ',
     &            ' d(rb/2V)  ',
     &            (' d(', dname(k),')', k=1, ncontrol),
     &            (' d(', jname(k),')', k=1, nvarjet)
 1902     format(1x,a4,a6, 6a11,1x, 30(a3,a8,a1) )
         endif
         write(*,1903) iter, 
     &                 rlx,
     &                 dgsmax/cref,
     &                 dal/dtr,
     &                 dbe/dtr, 
     &                 dwbar(1)*bref/2.0, 
     &                 dwbar(2)*cref/2.0, 
     &                 dwbar(3)*bref/2.0,
     &                (ddc(k), k=1, ncontrol),
     &                (ddj(k), k=1, nvarjet)
 1903    format(1x,i3, f6.3, 6(e10.2,1x), 1x, 50(e10.2,2x))

c         write(*,1905) iter, 
c     &                 dgsmax/cref,
c     &                 dal/dtr,
c     &                 dbe/dtr, 
c     &                 dwbar(1)*bref/2.0, 
c     &                 dwbar(2)*cref/2.0, 
c     &                 dwbar(3)*bref/2.0
c         write(*,1910) (ddc(k), k=1, ncontrol)
c         write(*,1920) (ddj(k), k=1, nvarjet)
c 1905    format(1x,i3,'dGam =', e11.3,'  da =', e11.3,'  db =', e11.3,
c     &          1x,3x,'  dp =', e11.3,'  dq =', e11.3,'  dr =', e11.3)
c 1910    format(1x,3x,'ddel =', 50e11.3 )
c 1920    format(1x,3x,'djet =', 50e11.3 )
        endif

        if(iter .gt. 1+nitmax/2) then
c------ if changes are still too big after half the max # of iterations, 
c-        configuration is probably untrimmable
         if(abs(alfa+dal) .gt. dmaxa) then
          write(*,*) 'Cannot trim.  Alpha too large.  a =',
     &                (alfa+dal)/dtr
          return
         endif

         if(abs(beta+dbe) .gt. dmaxa) then
          write(*,*) 'Cannot trim.  Beta too large.  b =',
     &                (beta+dbe)/dtr
          return
         endif

         if(abs(wbar(1)+dwbar(1))*bref .gt. dmaxr) then
          write(*,*) 'Cannot trim.  Roll rate too large.  pb/2V =', 
     &                (wbar(1)+dwbar(1))*bref*0.5
          return
         endif

         if(abs(wbar(2)+dwbar(2))*cref .gt. dmaxr) then
          write(*,*) 'Cannot trim.  Pitch rate too large.  qc/2V =',
     &                (wbar(2)+dwbar(2))*cref*0.5
          return
         endif

         if(abs(wbar(3)+dwbar(3))*bref .gt. dmaxr) then
          write(*,*) 'Cannot trim.  Yaw rate too large.  rb/2V =',
     &                (wbar(3)+dwbar(3))*bref*0.5
          return
         endif
        endif

c------ update
        do i = 1, nvor
          gam(i) = gam(i) + rlx*dgam(i)
        enddo
        vbar(1) = vbar(1) + rlx*dvbar(1)
        vbar(2) = vbar(2) + rlx*dvbar(2)
        vbar(3) = vbar(3) + rlx*dvbar(3)

c        dalfa = alfa_u(1)*dvbar(1)
c     &        + alfa_u(2)*dvbar(2)
c     &        + alfa_u(3)*dvbar(3)
c        dbeta = beta_u(1)*dvbar(1)
c     &        + beta_u(2)*dvbar(2)
c     &        + beta_u(3)*dvbar(3)
c        alfa = alfa + rlx*dalfa
c        beta = beta + rlx*dbeta
c        vbar(1) = cos(alfa)*cos(beta)
c        vbar(2) = -sin(beta)
c        vbar(3) = sin(alfa)*cos(beta)

        wbar(1) = wbar(1) + rlx*dwbar(1)
        wbar(2) = wbar(2) + rlx*dwbar(2)
        wbar(3) = wbar(3) + rlx*dwbar(3)
        do n = 1, ncontrol
          delcon(n) = delcon(n) + rlx*ddc(n)
        enddo
        do n = 1, nvarjet
          deljet(n) = deljet(n) + rlx*ddj(n)
        enddo

        do ib = 1, nbody
          do l = lfrst(ib), llast(ib)-1
            src(l) = src_u(l,1)*vbar(1)
     &             + src_u(l,2)*vbar(2)
     &             + src_u(l,3)*vbar(3)
     &             + src_u(l,4)*wbar(1)
     &             + src_u(l,5)*wbar(2)
     &             + src_u(l,6)*wbar(3)
            do k = 1, 3
              dbl(k,l) = dbl_u(k,l,1)*vbar(1)
     &                 + dbl_u(k,l,2)*vbar(2)
     &                 + dbl_u(k,l,3)*vbar(3)
     &                 + dbl_u(k,l,4)*wbar(1)
     &                 + dbl_u(k,l,5)*wbar(2)
     &                 + dbl_u(k,l,6)*wbar(3)
            enddo
          enddo
        enddo

c------- optional rescaling to force unit normalized freestream
c        vmag = sqrt(vbar(1)**2 + vbar(2)**2 + vbar(3)**2)
c        write(*,*) '|V| =', vmag
c        vbar(1) = vbar(1)/vmag
c        vbar(2) = vbar(2)/vmag
c        vbar(3) = vbar(3)/vmag

c------ store updated variables in run-case array
        call varset(ir)

c------ convergence check
        delmax = max( abs(dal), 
     &                abs(dbe),
     &                abs(dwbar(1)*bref/2.0),
     &                abs(dwbar(2)*cref/2.0),
     &                abs(dwbar(3)*bref/2.0),
     &                abs(dgsmax/cref) )
c        do n = 1, ncontrol
c          delmax = max( delmax , abs(ddc(n)) )
c        enddo
        do n = 1, nvarjet
          delmax = max( delmax , abs(ddj(n)) )
        enddo

        if(delmax.lt.eps) then
         lsol = .true.
c------- mark trim case as being converged
         itrim(ir) = iabs(itrim(ir))
         lsolr(ir) = .true.
         go to 510
        endif

 500  continue
      if(niter.gt.0) then
       write(*,*) 'Trim convergence failed'
       lsol = .false.
       return
      endif

c------------------------------------------------------
c---- pick up here on convergence
 510  continue

c---- compute final velocities, forces, moments
      call vvsum
      call vcsum
      call vlsum

      call wvsum
      call wcsum
      call wlsum

      call aero

c---- store values in run-case array
      parval(ipalfa,ir) = alfa/dtr
      parval(ipbeta,ir) = beta/dtr
      parval(iprotx,ir) = wbar(1)*0.5*bref
      parval(iproty,ir) = wbar(2)*0.5*cref
      parval(iprotz,ir) = wbar(3)*0.5*bref
      parval(ipcl  ,ir) = ltot/(0.5*sref)

      return
      end ! exec

      subroutine varini(ir)
c------------------------------------------------------------
c     Set working variables from arrays for run ir,
c------------------------------------------------------------
      include 'jvl.inc'

      vbar(1) = varval(ivvelx,ir)
      vbar(2) = varval(ivvely,ir)
      vbar(3) = varval(ivvelz,ir)
      wbar(1) = varval(ivomgx,ir)
      wbar(2) = varval(ivomgy,ir)
      wbar(3) = varval(ivomgz,ir)

c---- set alfa,beta from vbar(.)
      call vbar2ab

      do n = 1, ncontrol
        iv = n + 6
        delcon(n) = varval(iv,ir)
      enddo

      do n = 1, nvarjet
        iv = n + 6 + ncontrol
        deljet(n) = varval(iv,ir)
      enddo

      do n = 1, ndesign
c       iv = n + 6 + ncontrol + nvarjet
c       deldes(n) = varval(iv,ir)
        deldes(n) = 0.
      enddo

      lsolr(ir) = .false.
      lsol = .false.

      return
      end ! varini


      subroutine varset(ir)
c------------------------------------------------------------
c     Set stored variable arrays for run ir,
c     from current working variables
c------------------------------------------------------------
      include 'jvl.inc'

      varval(ivvelx,ir) = vbar(1)
      varval(ivvely,ir) = vbar(2)
      varval(ivvelz,ir) = vbar(3)
      varval(ivomgx,ir) = wbar(1)
      varval(ivomgy,ir) = wbar(2)
      varval(ivomgz,ir) = wbar(3)

      do n = 1, ncontrol
        iv = n + 6
        varval(iv,ir) = delcon(n)
      enddo

      do n = 1, nvarjet
        iv = n + 6 + ncontrol
        varval(iv,ir) = deljet(n)
      enddo

c      do n = 1, ndesign
c       iv = n + 6 + ncontrol + nvarjet
c       varval(iv,ir) = deldes(n)
c      enddo

      return
      end ! varset

