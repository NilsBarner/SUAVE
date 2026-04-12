C***********************************************************************
C    Module:  jmake.f
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

      subroutine makesurf(isurf, ibx,nsec, 
     &       nvc1,cspace, nvs1,sspace, nvcu1,cspaceu, nvcw1,cspacew,
     &       xyzscal,xyztran,addinc,
     &       xyzles,chords,aincs,
     &       hdisks,fhjets,djets0,djets1,djets3,dxdisks,dndisks,tdisks,
     &       sspaces,nspans,
     &       xasec,sasec,tasec,nasec,
     &       cdclsec,claf,
     &       iconx, 
     &       icontd,nscon,gaind,xhinged,vhinged,refld,
     &       ijettd,nsjet,gainj,reflj,
     &       idestd,nsdes,gaing,
     &       rho,izimage,
     &       ldupl,ydupl)
c--------------------------------------------------------------
c     Sets up all stuff for surface isurf, 
c     using info from configuration input file.
c--------------------------------------------------------------
      include 'jvl.inc'

      real xyzscal(3), xyztran(3), addinc,
     &     xyzles(3,*),chords(*),aincs(*),
     &     hdisks(*),fhjets(*),
     &     djets0(*),djets1(*),djets3(*),
     &     dxdisks(*),dndisks(*),tdisks(3,*),
     &     sspaces(*)
      integer nspans(*), nasec(*)
      real xasec(ibx,*), sasec(ibx,*), tasec(ibx,*)
      real cdclsec(6,*), claf(*)
      integer icontd(iconx,*),nscon(*),
     &        ijettd(iconx,*),nsjet(*),
     &        idestd(iconx,*),nsdes(*)
      real gaind(iconx,*),xhinged(iconx,*),vhinged(3,iconx,*),
     &     refld(iconx,*)
      real gainj(iconx,*),
     &     reflj(iconx,*)
      real gaing(iconx,*)
      logical ldupl
      real ydupl

c--------------------------------------------------------------
c---- local arrays
      real xyzlel(3), xyzler(3)

      parameter (kcmax=50,
     &           ksmax=500)
      real xpt(kcmax), xcp(kcmax), xvr(kcmax), xsr(kcmax),
     &     ypt(ksmax), ycp(ksmax)
      real yzlen(ksmax)
      integer iptloc(ksmax)

      parameter (kpmax=4*kcmax+2*ksmax)
      real fspace(kpmax)

      real chsinl_g(ngmax),chcosl_g(ngmax),
     &     chsinr_g(ngmax),chcosr_g(ngmax)
      integer isconl(ndmax), isconr(ndmax),
     &        isjetl(njmax), isjetr(njmax)
      real xled(ndmax), xted(ndmax), gainda(ndmax)

      real tdiskl(3),
     &     tdiskr(3)

      if(nsec.lt.2) then
       write(*,*) '*** Need at least 2 SECTIONS per SURFACE.'
       stop
      endif

      nvc = nvc1
      nvs = nvs1
      nvcu = nvcu1
      nvcw = nvcw1

      if(nvc.gt.kcmax) then
       write(*,*) '* MAKESURF: Array overflow.  Increase kcmax to', nvc
       nvc = kcmax
      endif

      if(nvs.gt.ksmax) then
       write(*,*) '* MAKESURF: Array overflow.  Increase ksmax to', nvs
       nvs = ksmax
      endif

      if(nvcu.gt.kcmax) then
       write(*,*) '* MAKESURF: Array overflow.  Increase kcmax to', nvcu
       nvcu = kcmax
      endif

      if(nvcw.gt.kcmax) then
       write(*,*) '* MAKESURF: Array overflow.  Increase kcmax to', nvcw
       nvcw = kcmax
      endif

c--- image flag set to indicate section definition direction
c    idups= 1  defines edge 1 located at surface root edge 
c    idups=-1  defines edge 2 located at surface root edge (reflected surfaces)
      idups(isurf) = 1

      jfrst(isurf) = nstrip + 1

      nk(isurf) = nvc

c-----------------------------------------------------------------
c---- section arc lengths of wing trace in y-z plane
      yzlen(1) = 0.
      do isec = 2, nsec
        dy = xyzles(2,isec) - xyzles(2,isec-1)
        dz = xyzles(3,isec) - xyzles(3,isec-1)
        yzlen(isec) = yzlen(isec-1) + sqrt(dy*dy + dz*dz)
      enddo

c
      if(nvs1.eq.0) then
c----- set spanwise spacing using spacing parameters for each section interval
       nvs = 0
       do isec = 1, nsec-1
         nvs = nvs + nspans(isec)
       enddo
       if(nvs.gt.ksmax) then
        write(*,*) '*** MAKESURF: Array overflow. Increase ksmax to',nvs
        stop
       endif

       nvs = 0
       ypt(1) = yzlen(1)
       iptloc(1) = 1

       do isec = 1, nsec-1
         dyzlen = yzlen(isec+1) - yzlen(isec)

         nvint = nspans(isec)

c------- set spanwise spacing array
         nspace = 2*nvint + 1
         if(nspace.gt.kpmax) then
          write(*,*) '*** MAKESURF: Array overflow. Increase kpmax to', 
     &                 nspace
          stop
         endif
         call sspacer(nspace,sspaces(isec),fspace)

         do n = 1, nvint
           ivs = nvs + n
           ycp(ivs)   = ypt(nvs+1) + dyzlen*fspace(2*n)
           ypt(ivs+1) = ypt(nvs+1) + dyzlen*fspace(2*n+1)
         enddo
         iptloc(isec+1) = nvs + nvint + 1

         nvs = nvs + nvint
       enddo

      else
c----- set spanwise spacing using overall parameters nvs, sspace

       nspace = 2*nvs + 1
       if(nspace.gt.kpmax) then
        write(*,*) '*** MAKESURF: Array overflow. Increase kpmax to', 
     &              nspace
        stop
       endif
       call sspacer(nspace,sspace,fspace)

       ypt(1) = yzlen(1)
       do ivs = 1, nvs
         ycp(ivs)   = yzlen(1) + (yzlen(nsec)-yzlen(1))*fspace(2*ivs)
         ypt(ivs+1) = yzlen(1) + (yzlen(nsec)-yzlen(1))*fspace(2*ivs+1)
       enddo

       npt = nvs + 1

c----- find node nearest each section
       do isec = 2, nsec-1
         yptloc = 1.0e9
         iptloc(isec) = 1
         do ipt = 1, npt
           yptdel = abs(yzlen(isec) - ypt(ipt))
           if(yptdel .lt. yptloc) then
            yptloc = yptdel
            iptloc(isec) = ipt
           endif
         enddo
       enddo
       iptloc(1)    = 1
       iptloc(nsec) = npt

c----- fudge glauert angles to make nodes match up exactly with interior sections
       do isec = 2, nsec-1
         ipt1 = iptloc(isec-1)
         ipt2 = iptloc(isec  )
         if(ipt1.eq.ipt2) then
          call bstrip(stitle(isurf),nst)
          write(*,7000) isec, stitle(isurf)(1:nst)
          stop
         endif

         ypt1 = ypt(ipt1)
         yscale = (yzlen(isec)-yzlen(isec-1)) / (ypt(ipt2)-ypt(ipt1))
         do ipt = ipt1, ipt2-1
           ypt(ipt) = yzlen(isec-1) + yscale*(ypt(ipt)-ypt1)
         enddo
         do ivs = ipt1, ipt2-1
           ycp(ivs) = yzlen(isec-1) + yscale*(ycp(ivs)-ypt1)
         enddo

         ipt1 = iptloc(isec  )
         ipt2 = iptloc(isec+1)
         if(ipt1.eq.ipt2) then
          call bstrip(stitle(isurf),nst)
          write(*,7000) isec, stitle(isurf)(1:nst)
          stop
         endif

         ypt1 = ypt(ipt1)
         yscale = (ypt(ipt2)-yzlen(isec)) / (ypt(ipt2)-ypt(ipt1))
         do ipt = ipt1, ipt2-1
           ypt(ipt) = yzlen(isec) + yscale*(ypt(ipt)-ypt1)
         enddo
         do ivs = ipt1, ipt2-1
           ycp(ivs) = yzlen(isec) + yscale*(ycp(ivs)-ypt1)
         enddo

 7000    format(
     &   /' *** Cannot adjust spanwise spacing at SECTION', i3, 
     &    ', on SURFACE ', a
     &   /' *** Insufficient number of spanwise vortices to work with')
       enddo

      endif

c====================================================
c---- define strips between input sections

      nj(isurf) = 0

      if(ncontrol.gt.ndmax) then
       write(*,*) '*** Too many control variables.  Increase ndmax to',
     &            ncontrol
       stop
      endif

      if(nvarjet.gt.njmax) then
       write(*,*) '*** Too many jet variables.  Increase njmax to',
     &            nvarjet
       stop
      endif

      if(ndesign.gt.ngmax) then
       write(*,*) '*** Too many design variables.  Increase ngmax to',
     &            ndesign
       stop
      endif

c---- go over section intervals
      do 200 isec = 1, nsec-1
        xyzlel(1) = xyzscal(1)*xyzles(1,isec)    + xyztran(1)
        xyzlel(2) = xyzscal(2)*xyzles(2,isec)    + xyztran(2)
        xyzlel(3) = xyzscal(3)*xyzles(3,isec)    + xyztran(3)
        xyzler(1) = xyzscal(1)*xyzles(1,isec+1)  + xyztran(1)
        xyzler(2) = xyzscal(2)*xyzles(2,isec+1)  + xyztran(2)
        xyzler(3) = xyzscal(3)*xyzles(3,isec+1)  + xyztran(3)

        width = sqrt(  (xyzler(2)-xyzlel(2))**2
     &               + (xyzler(3)-xyzlel(3))**2 )

        chordl = xyzscal(1)*chords(isec)
        chordr = xyzscal(1)*chords(isec+1)

        clafl = claf(isec)
        clafr = claf(isec+1)

        aincl = aincs(isec)   + addinc - 4.0*dtr*(clafl-1.0)
        aincr = aincs(isec+1) + addinc - 4.0*dtr*(clafr-1.0)

        chsinl = chordl*sin(aincl)
        chsinr = chordr*sin(aincr)
        chcosl = chordl*cos(aincl)
        chcosr = chordr*cos(aincr)

        hdiskl = hdisks(isec)
        hdiskr = hdisks(isec+1)

        fhjetl = fhjets(isec)
        fhjetr = fhjets(isec+1)

        djet0l = djets0(isec)
        djet0r = djets0(isec+1)

        djet1l = djets1(isec)
        djet1r = djets1(isec+1)

        djet3l = djets3(isec)
        djet3r = djets3(isec+1)

        dxdiskl = dxdisks(isec)
        dxdiskr = dxdisks(isec+1)

        dndiskl = dndisks(isec)
        dndiskr = dndisks(isec+1)

        do k = 1, 3
          tdiskl(k) = tdisks(k,isec)
          tdiskr(k) = tdisks(k,isec+1)
        enddo

c------ set control-declaration lines for each control variable
        do n = 1, ncontrol
          isconl(n) = 0
          isconr(n) = 0
          do iscon = 1, nscon(isec)
            if(icontd(iscon,isec)  .eq.n) isconl(n) = iscon
          enddo
          do iscon = 1, nscon(isec+1)
            if(icontd(iscon,isec+1).eq.n) isconr(n) = iscon
          enddo
        enddo

c------ set jet control-declaration lines for each jet control variable
        do n = 1, nvarjet
          isjetl(n) = 0
          isjetr(n) = 0
c          write(*,*)
c          write(*,*) 'isec', isec
          do isjet = 1, nsjet(isec)
            if(ijettd(isjet,isec)  .eq.n) isjetl(n) = isjet
c            write(*,*) 'isjet ijettd', isjet, ijettd(isjet,isec), n,
c     &                  isjetl(n)
          enddo
          do isjet = 1, nsjet(isec+1)
            if(ijettd(isjet,isec+1).eq.n) isjetr(n) = isjet
          enddo
        enddo

c------ set design-variable sensitivities of chsin and chcos
        do n = 1, ndesign
          chsinl_g(n) = 0.
          chsinr_g(n) = 0.
          chcosl_g(n) = 0.
          chcosr_g(n) = 0.

          do isdes = 1, nsdes(isec)
            if(idestd(isdes,isec).eq.n) then
             chsinl_g(n) =  chcosl * gaing(isdes,isec)*dtr
             chcosl_g(n) = -chsinl * gaing(isdes,isec)*dtr
            endif
          enddo

          do isdes = 1, nsdes(isec+1)
            if(idestd(isdes,isec+1).eq.n) then
             chsinr_g(n) =  chcosr * gaing(isdes,isec+1)*dtr
             chcosr_g(n) = -chsinr * gaing(isdes,isec+1)*dtr
            endif
          enddo
        enddo

c------ go over chord strips
        iptl = iptloc(isec)
        iptr = iptloc(isec+1)
        nspan = iptr - iptl
        do 150 ispan = 1, nspan
c-------- define left and right edges of vortex strip
c-          note that incidence angle is set by atan of chord projections,
c-          not by linear interpolation of ainc
          ipt1 = iptl + ispan - 1
          ipt2 = iptl + ispan
          ivs  = iptl + ispan - 1
          f1 = (ypt(ipt1)-ypt(iptl))/(ypt(iptr)-ypt(iptl))
          f2 = (ypt(ipt2)-ypt(iptl))/(ypt(iptr)-ypt(iptl))
          fc = (ycp(ivs) -ypt(iptl))/(ypt(iptr)-ypt(iptl))

c-------- store strip in global data arrays
          nstrip = nstrip + 1
          nj(isurf) = nj(isurf) + 1

          j = nstrip

          rle1(1,j) = (1.0-f1)*xyzlel(1) + f1*xyzler(1)
          rle1(2,j) = (1.0-f1)*xyzlel(2) + f1*xyzler(2)
          rle1(3,j) = (1.0-f1)*xyzlel(3) + f1*xyzler(3)
          chord1(j) = (1.0-f1)*chordl    + f1*chordr

          rle2(1,j) = (1.0-f2)*xyzlel(1) + f2*xyzler(1)
          rle2(2,j) = (1.0-f2)*xyzlel(2) + f2*xyzler(2)
          rle2(3,j) = (1.0-f2)*xyzlel(3) + f2*xyzler(3)
          chord2(j) = (1.0-f2)*chordl    + f2*chordr

          rle(1,j)  = (1.0-fc)*xyzlel(1) + fc*xyzler(1)
          rle(2,j)  = (1.0-fc)*xyzlel(2) + fc*xyzler(2)
          rle(3,j)  = (1.0-fc)*xyzlel(3) + fc*xyzler(3)
          chord(j)  = (1.0-fc)*chordl    + fc*chordr

          wstrip(j) = abs(f2-f1)*width
          tanle(j)  = (xyzler(1)-xyzlel(1))/width
          tante(j)  = (xyzler(1)+chordr - xyzlel(1)-chordl)/width

          chsin = (1.0-fc)*chsinl + fc*chsinr
          chcos = (1.0-fc)*chcosl + fc*chcosr
          ainc(j) = atan2(chsin,chcos)

          do n = 1, ndesign
            chsin_g = (1.0-fc)*chsinl_g(n) + fc*chsinr_g(n)
            chcos_g = (1.0-fc)*chcosl_g(n) + fc*chcosr_g(n)
            ainc_g(j,n) = (chcos*chsin_g - chsin*chcos_g)
     &                       / (chsin**2 + chcos**2)
          enddo

c-------- strip has same density factor as its surface
          rhos(j) = rho

c-------- strip z-image flag
          izims(j) = izimage
          
c-------- surface strip
          istype(j) = 0

c-------- assume strip will not have a jet sheet
          ifrstu(j) = 0
          ilastu(j) = -1
          ifrstw(j) = 0
          ilastw(j) = -1

c-------- interpolate jet sheet parameters (if any)
          hdstrp(j)  = (1.0-fc)*hdiskl + fc*hdiskr
          fhstrp(j)  = (1.0-fc)*fhjetl + fc*fhjetr
          dj0strp(j) = (1.0-fc)*djet0l + fc*djet0r
          dj1strp(j) = (1.0-fc)*djet1l + fc*djet1r
          dj3strp(j) = (1.0-fc)*djet3l + fc*djet3r

          dxdstrp(j)  = (1.0-fc)*dxdiskl + fc*dxdiskr
          dxdstrp1(j) = (1.0-f1)*dxdiskl + f1*dxdiskr
          dxdstrp2(j) = (1.0-f2)*dxdiskl + f2*dxdiskr

          dndstrp(j)  = (1.0-fc)*dndiskl + fc*dndiskr
          dndstrp1(j) = (1.0-f1)*dndiskl + f1*dndiskr
          dndstrp2(j) = (1.0-f2)*dndiskl + f2*dndiskr

          taxstrp(1,j) = (1.0-fc)*tdiskl(1) + fc*tdiskr(1)
          taxstrp(2,j) = (1.0-fc)*tdiskl(2) + fc*tdiskr(2)
          taxstrp(3,j) = (1.0-fc)*tdiskl(3) + fc*tdiskr(3)
          taxdim = sqrt( taxstrp(1,j)**2
     &                 + taxstrp(2,j)**2
     &                 + taxstrp(3,j)**2 )
          taxstrp(1,j) = taxstrp(1,j)/taxdim
          taxstrp(2,j) = taxstrp(2,j)/taxdim
          taxstrp(3,j) = taxstrp(3,j)/taxdim

          do n = 1, ncontrol
            icl = isconl(n)
            icr = isconr(n)

            if(icl.eq.0 .or. icr.eq.0) then
c----------- no control effect here
             gainda(n) = 0.
             xled(n) = 0.
             xted(n) = 0.

             vhinge(1,j,n) = 0.
             vhinge(2,j,n) = 0.
             vhinge(3,j,n) = 0.

             reflsd(j,n) = 0.

             phinge(1,j,n) = 0.
             phinge(2,j,n) = 0.
             phinge(3,j,n) = 0.

            else
c----------- control variable # n is active here
             gainda(n) = gaind(icl,isec  )*(1.0-fc)
     &                 + gaind(icr,isec+1)*     fc
             gconmax(n) = max( dtr*abs(gainda(n)) , gconmax(n) )

             xhd = chordl*xhinged(icl,isec  )*(1.0-fc)
     &           + chordr*xhinged(icr,isec+1)*     fc
             if(xhd.ge.0.0) then
c------------ te control surface, with hinge at xhd
              xled(n) = xhd
              xted(n) = chord(j)
             else
c------------ le control surface, with hinge at -xhd
              xled(n) =  0.0
              xted(n) = -xhd
             endif

             vhx = vhinged(1,icl,isec)*xyzscal(1)
             vhy = vhinged(2,icl,isec)*xyzscal(2)
             vhz = vhinged(3,icl,isec)*xyzscal(3)
             vsq = vhx**2 + vhy**2 + vhz**2
             if(vsq.eq.0.0) then
c------------ default: set hinge vector along hingeline
              vhx = xyzles(1,isec+1) + abs(chordr*xhinged(icr,isec+1))
     &            - xyzles(1,isec  ) - abs(chordl*xhinged(icl,isec  ))
              vhy = xyzles(2,isec+1)
     &            - xyzles(2,isec  )
              vhz = xyzles(3,isec+1)
     &            - xyzles(3,isec  )
              vhx = vhx*xyzscal(1)
              vhy = vhy*xyzscal(2)
              vhz = vhz*xyzscal(3)
              vsq = vhx**2 + vhy**2 + vhz**2
             endif

             vmod = sqrt(vsq)
             vhinge(1,j,n) = vhx/vmod
             vhinge(2,j,n) = vhy/vmod
             vhinge(3,j,n) = vhz/vmod

             reflsd(j,n) = refld(icl,isec)

             if(xhd .ge. 0.0) then
              phinge(1,j,n) = rle(1,j) + xhd
              phinge(2,j,n) = rle(2,j)
              phinge(3,j,n) = rle(3,j)
             else
              phinge(1,j,n) = rle(1,j) - xhd
              phinge(2,j,n) = rle(2,j)
              phinge(3,j,n) = rle(3,j)
             endif

            endif
          enddo

          do n = 1, nvarjet
            icl = isjetl(n)
            icr = isjetr(n)

c           write(*,*) n, icl, icr

            if(icl.eq.0 .or. icr.eq.0) then
c----------- jet control variable n has no effect here
             reflsj(j,n) = 0.
             gjstrp(j,n) = 0.

            else
c----------- jet control variable n is active here
             reflsj(j,n) = reflj(icl,isec)
             gjstrp(j,n) = gainj(icl,isec  )*(1.0-fc)
     &                   + gainj(icr,isec+1)*     fc
             gjetmax(n) = max( abs(gjstrp(j,n)) , gjetmax(n) )

             ilastu(j) = -2
             ilastw(j) = -2
            endif
          enddo

c-------- interpolate cd-cl polar from input sections to strip
          do l = 1, 6
            cdcls(l,j) = (1.0-fc)*cdclsec(l,isec) 
     &                 +      fc *cdclsec(l,isec+1)
          end do

c-------- if the min drag is zero, flag the strip as no-viscous data
          lviscstrp(j) = (cdcls(4,j) .ne. 0.0)

          ifrsts(j) = nvor + 1

          isurfs(j) = isurf

          nsl = nasec(isec  )
          nsr = nasec(isec+1)

          chordc = chord(j)

          clafc = (1.0-fc)*(chordl/chordc)*clafl
     &          +      fc *(chordr/chordc)*clafr

c-------- set chordwise spacing fraction arrays
          call cspacer(nvc,cspace,clafc, xpt,xvr,xsr,xcp)

c-------- go over vortices in this strip
          do 1505 ivc = 1, nvc
            nvor = nvor + 1
            i = nvor

            rv1(1,i) = rle1(1,j) + xvr(ivc)*chord1(j)
            rv1(2,i) = rle1(2,j)
            rv1(3,i) = rle1(3,j)

            rv2(1,i) = rle2(1,j) + xvr(ivc)*chord2(j)
            rv2(2,i) = rle2(2,j)
            rv2(3,i) = rle2(3,j)

            rv(1,i) = rle(1,j) + xvr(ivc)*chordc
            rv(2,i) = rle(2,j)
            rv(3,i) = rle(3,j)

            rc(1,i) = rle(1,j) + xcp(ivc)*chordc
            rc(2,i) = rle(2,j)
            rc(3,i) = rle(3,j)

            rs(1,i) = rle(1,j) + xsr(ivc)*chordc
            rs(2,i) = rle(2,j)
            rs(3,i) = rle(3,j)

c            if(j.eq.1) then
c            write(*,*) ivc, xvr(ivc)/(1.0-fac*xvr(ivc)), 
c     &                      xcp(ivc)/(1.0-fac*xcp(ivc))
c            endif


            call akima(xasec(1,isec  ),sasec(1,isec  ),nsl,
     &                 xcp(ivc),slopel, dsdx)
            call akima(xasec(1,isec+1),sasec(1,isec+1),nsr,
     &                 xcp(ivc),sloper, dsdx)
            slopec(i) =  (1.-fc)*(chordl/chordc)*slopel 
     &                +     fc *(chordr/chordc)*sloper

            call akima(xasec(1,isec  ),sasec(1,isec  ),nsl,
     &                 xvr(ivc),slopel, dsdx)
            call akima(xasec(1,isec+1),sasec(1,isec+1),nsr,
     &                 xvr(ivc),sloper, dsdx)
            slopev(i) =  (1.-fc)*(chordl/chordc)*slopel 
     &                +     fc *(chordr/chordc)*sloper

            dxoc = xpt(ivc+1) - xpt(ivc)
            dxv(i) = dxoc*chordc
            chordv(i) = chordc

c---------- element has same density factor and z-image flag as its strip
            rhov(i) = rhos(j)
            izimv(i) = izims(j)

c---------- element inherits alpha,beta flag from surface
            lvalbe(i) = lfalbe(isurf)

c---------- pointers indicating strip and surface containing this element
            jstripv(i) = j
            isurfv(i) = isurf

c---------- component index of this element is same as that of the surface
            lscompv(i) = lscomp(isurf)

            do n = 1, ncontrol
c------------ scale control gain by factor 0..1, 
c-            (fraction of element on control surface)
              fracle = (xled(n)/chordc-xpt(ivc)) / dxoc
              fracte = (xted(n)/chordc-xpt(ivc)) / dxoc

              fracle = min( 1.0 , max( 0.0 , fracle ) )
              fracte = min( 1.0 , max( 0.0 , fracte ) )

              dcontrol(i,n) = gainda(n)*(fracte-fracle)
            enddo

            ijetm(i) = i-1
 1505     continue ! next vortex in strip

c-------- last-element index in this strip on surface
          ilasts(j) = nvor

c-------- assume LE element has no upstream element
          i = ifrsts(j)
          ijetm(i) = -999
 150    continue ! next strip j
 200  continue ! next section interval isec...isec+1

c---- last-strip index in this surface
      jlast(isurf) = j

c---- add any required upstream and wake jet strips
      do j = jfrst(isurf), jlast(isurf)
        ldstrp(j) = .false.

        if(ilastu(j) .eq. -2) then
c------- x/c interval at LE, to be matched by jet strip spacing
         dxle1 = (xvr(1) - xpt(1))*chord1(j)
         dxle2 = (xvr(1) - xpt(1))*chord2(j)
         dxle  = (xvr(1) - xpt(1))*chord(j)
         call makejetu(j,dxle1,dxle2,dxle,nvcu,cspaceu)
         ldstrp(j) = .true.
        endif

        i = ifrsts(j)
        if(ilastu(j)-ifrstu(j) .ge. 0) then
c------- there is an upstream jet
         ijetm(i) = ilastu(j)
        else
c------- there is no upstream jet
         ijetm(i) = -999
        endif

        if(ilastw(j) .eq. -2) then
c------- x/c interval at TE, to be matched by jet strip spacing
         dxte1 = (xpt(nvc+1) - xcp(nvc))*chord1(j)
         dxte2 = (xpt(nvc+1) - xcp(nvc))*chord2(j)
         dxte  = (xpt(nvc+1) - xcp(nvc))*chord(j)
         call makejetw(j,dxte1,dxte2,dxte,nvcw,cspacew)
         ldstrp(j) = .true.
        endif
      enddo

c---- find wetted surface area (one side)
      sum  = 0.0
      wtot = 0.0
      do j = jfrst(isurf), jlast(isurf)
        astrp = wstrip(j)*chord(j)
        sum  = sum + astrp
        wtot = wtot + wstrip(j)
      enddo
      ssurf(isurf) = sum

      if(wtot .eq. 0.0) then
       cavesurf(isurf) = 0.0
      else
       cavesurf(isurf) = sum/wtot
      endif

      rhon(isurf) = rho
      
      if(ldupl) then
       call sdupl(isurf,ydupl,'YDUP')
      endif

      return
      end ! makesurf



      subroutine makebody(ibody, ibx,
     &       nvb1, bspace,
     &       xyzscal,xyztran,
     &       xbod,ybod,abod,nbod,
     &       rho,izimage,
     &       ldupl,ydupl)

c--------------------------------------------------------------
c     Sets up all stuff for body ibody,
c     using info from configuration input file.
c--------------------------------------------------------------
      include 'jvl.inc'

      real xyzscal(3), xyztran(3)
      real xbod(nbod), ybod(nbod), abod(nbod)
      logical ldupl
      real ydupl

      parameter (klmax=201)
      real xpt(klmax), xcp(klmax)
      real fspace(2*klmax)

c      if(nsec.lt.2) then
c       write(*,*) '*** Need at least 2 SECTIONS per body.'
c       stop
c      endif

      nvb = nvb1

      if(nvb.gt.klmax) then
       write(*,*) '* MAKEBODY: Array overflow.  Increase klmax to', nvb
       nvb = klmax
      endif


      lfrst(ibody) = nlbody + 1 

      if(nlbody+nvb+1 .gt. nlmax) then
       write(*,*) '*** MAKEBODY: Array overflow. Increase nlmax to',
     &             nlbody+nvb+1
       stop
      endif

      ascal = xyzscal(2)*xyzscal(3)

c---- assume body is on centerline
      ibcent(ibody) = 1

      cbody(ibody) = xbod(nbod) - xbod(1)
      sbody(ibody) = 0.
      vbody(ibody) = 0.

c-----------------------------------------------------------------
c---- lengthwise spacing fraction arrays
      nspace = 2*nvb + 1
      if(nspace .gt. 2*klmax) then
       write(*,*) '*** MAKEBODY: Array overflow. Increase klmax to', 
     &             nspace/2
       stop
      endif
      call sspacer(nspace,bspace,fspace)

      do ivb = 1, nvb
        xpt(ivb) = fspace(2*ivb-1)
        xcp(ivb) = fspace(2*ivb  )
      enddo
      xpt(1) = 0.0
      xpt(nvb+1) = 1.0

c      do j = 1, nbod
c        write(10,*) xbod(j), ybod(j), abod(j)
c      enddo

c---- body nodes, radii, area changes
      do ivb = 1, nvb
        nlbody = nlbody + 1
        l = nlbody

        lbcompl(l) = lbcomp(ibody)

        xvb1 = xbod(1) + (xbod(nbod)-xbod(1))*xpt(ivb)
        xvb2 = xbod(1) + (xbod(nbod)-xbod(1))*xpt(ivb+1)
        xvb  = xbod(1) + (xbod(nbod)-xbod(1))*xcp(ivb)

        call akima(xbod,ybod,nbod,xvb1,yvb1,dydx)
        call akima(xbod,abod,nbod,xvb1,avb1,dadx)

        call akima(xbod,ybod,nbod,xvb2,yvb2,dydx)
        call akima(xbod,abod,nbod,xvb2,avb2,dadx)

        call akima(xbod,ybod,nbod,xvb,yvb,dydx)
        call akima(xbod,abod,nbod,xvb,avb,dadx)

c        write(*,*) ivb, nvb, avb1, avb, avb2
c        write(11,*) xvb1, yvb1, avb1
c        write(12,*) xvb2, yvb2, avb2
c        write(13,*) xvb, yvb, avb

        rl1(1,l) = xyztran(1) + xyzscal(1)*xvb1
        rl1(2,l) = xyztran(2)
        rl1(3,l) = xyztran(3) + xyzscal(3)*yvb1
        rad1(l) = sqrt(ascal*avb1/pi)

        rl2(1,l) = xyztran(1) + xyzscal(1)*xvb2
        rl2(2,l) = xyztran(2)
        rl2(3,l) = xyztran(3) + xyzscal(3)*yvb2
        rad2(l) = sqrt(ascal*avb2/pi)

        rl(1,l)  = xyztran(1) + xyzscal(1)*xvb
        rl(2,l)  = xyztran(2)
        rl(3,l)  = xyztran(3) + xyzscal(3)*yvb

        avb = max( avb , 0.0 )
        darl(l) = ascal*(avb2-avb1)
        radl(l) = sqrt(ascal*avb/pi)
        drbl(l) = sqrt(ascal*avb2/pi) - sqrt(ascal*avb1/pi)

c        write(19,*) xvb, avb, radl(l), darl(l), drbl(l) 
c        write(20,*) xvb1, avb1, xvb2, avb2

        delx = (xvb2 - xvb1)*xyzscal(1)
        abody(ibody) = abody(ibody) + delx*2.0*radl(l)
        sbody(ibody) = sbody(ibody) + 
     &    pi*(rad1(l)+rad2(l)) * sqrt(delx**2 + (rad1(l)-rad2(l))**2)
        vbody(ibody) = vbody(ibody) + delx*ascal*avb

c------ segment has same density factor and z-image flag as its body
        rhol(l) = rho
        iziml(l) = izimage

        if(abs(rl(2,l)) .gt. 0.0001*bref) then
c------- body is not on centerline
         ibcent(ibody) = 0
        endif
      enddo
      llast(ibody) = l

c---- create body-wake doublet line segment
      nlbody = nlbody + 1
      l = nlbody

      lbcompl(l) = lbcomp(ibody)

      bwakel = max( 100.0*(xbod(nbod)-xbod(1)) , 100.0*bref )
      avb = max( abod(nbod) , 0.0 )

      rl1(1,l) = rl2(1,l-1)
      rl1(2,l) = rl2(2,l-1)
      rl1(3,l) = rl2(3,l-1)
      rad1(l) = rad2(l-1)

      rl2(1,l) = rl2(1,l-1) + bwakel
      rl2(2,l) = rl2(2,l-1)
      rl2(3,l) = rl2(3,l-1)
      rad2(l) = rad2(l-1)

      darl(l) = 0.
      radl(l) = sqrt(ascal*avb/pi)
      drbl(l) = 0.

      rhol(l) = rho
      iziml(l) = izimage
    
      rhob(ibody) = rho
      izimb(ibody) = izimage

      if(ldupl) then
       call bdupl(ibody,ydupl,'YDUP')
      endif

      return
      end ! makebody


      subroutine sdupl(nn,ydupl,msg)
c-----------------------------------------------------------------
c     Creates new image of surface nn, reflected about y=ydupl.
c     name of new surface has "( msg )" appended to it.
c-----------------------------------------------------------------
      include 'jvl.inc'
      character*(*) msg

      nni = nsurf + 1
      nni = min( nni , nsmax )
c      if(nni.gt.nfmax) then
c       write(*,*) 'SDUPL: Surface array overflow. Increase nfmax.'
c       stop
c      endif

      klen = len(stitle(nn))
      do k = klen, 1, -1
        if(stitle(nn)(k:k) .ne. ' ') go to 6
      enddo
 6    stitle(nni) = stitle(nn)(1:k) // ' (' // msg // ')'
      write(*,*) ' '
      write(*,*) '  Building duplicate image-surface: ',stitle(nni)

c---- duplicate surface is assumed to be the same logical surface
      lscomp(nni) = lscomp(nn)

c---- same various logical flags
      lfwake(nni) = lfwake(nn)
      lfalbe(nni) = lfalbe(nn)
      lfload(nni) = lfload(nn)

      jfrst(nni) = nstrip + 1
      nj(nni) = nj(nn)
      nk(nni) = nk(nn)

      nvc = nk(nni)
      nvs = nj(nni)

      ssurf(nni)    = ssurf(nn)
      cavesurf(nni) = cavesurf(nn)

c--- note hinge axis is flipped to reverse the y component of the hinge
c    vector.   this means that deflections need to be reversed for image
c    surfaces.

c--- image flag reversed (set to -idups) for imaged surfaces
      idups(nni) = -idups(nn)

      yoff = 2.0*ydupl

c--- create image strips, to maintain the same sense of positive gamma
c    these have the 1 and 2 strip edges reversed (i.e. root is edge 2, 
c    not edge 1 as for a strip with idups=1
      do 100 jj = jfrst(nn), jlast(nn)
        nstrip = nstrip + 1
        if(nstrip.gt.nsmax) then
          write(*,*) 'SDUPL: Strip array overflow. Increase nsmax.'
          stop
        endif

        jji = nstrip

        rle1(1,jji) =  rle2(1,jj)
        rle1(2,jji) = -rle2(2,jj) + yoff
        rle1(3,jji) =  rle2(3,jj)
        chord1(jji) =  chord2(jj)

        rle2(1,jji) =  rle1(1,jj)
        rle2(2,jji) = -rle1(2,jj) + yoff
        rle2(3,jji) =  rle1(3,jj)
        chord2(jji) =  chord1(jj)

        rle(1,jji) =  rle(1,jj)
        rle(2,jji) = -rle(2,jj) + yoff
        rle(3,jji) =  rle(3,jj)

        chord(jji)  =  chord(jj)
        wstrip(jji) =  wstrip(jj)
        tanle(jji)  = -tanle(jj)
        ainc (jji)  =  ainc(jj)

        ldstrp(jji) = ldstrp(jj)
        hdstrp(jji) = hdstrp(jj)
        fhstrp(jji) = fhstrp(jj)

        dj0strp(jji) = dj0strp(jj)
        dj1strp(jji) = dj1strp(jj)
        dj3strp(jji) = dj3strp(jj)

        dxdstrp(jji)  = dxdstrp(jj)
        dxdstrp1(jji) = dxdstrp1(jj)
        dxdstrp2(jji) = dxdstrp2(jj)

        dndstrp(jji)  = dndstrp(jj)
        dndstrp1(jji) = dndstrp1(jj)
        dndstrp2(jji) = dndstrp2(jj)

        taxstrp(1,jji) =  taxstrp(1,jj)
        taxstrp(2,jji) = -taxstrp(2,jj)
        taxstrp(3,jji) =  taxstrp(3,jj)

        rhos(jji) = rhos(jj)
        izims(jji) = izims(jj)
        
        isurfs(nstrip) = nni

        do n = 1, ndesign
          ainc_g(jji,n) = ainc_g(jj,n)
        enddo

        do n = 1, nvarjet
ccc       rsgn = sign( 1.0 , reflsj(jj,n) )
          rsgn = reflsj(jj,n)
          gjstrp(jji,n) = gjstrp(jj,n)*rsgn
        enddo

        do n = 1, ncontrol
          reflsd(jji,n) = reflsd(jj,n)

          vhinge(1,jji,n) =  vhinge(1,jj,n)
          vhinge(2,jji,n) = -vhinge(2,jj,n)
          vhinge(3,jji,n) =  vhinge(3,jj,n)

          phinge(1,jji,n) =  phinge(1,jj,n)
          phinge(2,jji,n) = -phinge(2,jj,n) + yoff
          phinge(3,jji,n) =  phinge(3,jj,n)
        enddo

c------ the defined section for image strip is flagged with (-)
        ifrsts(jji) = nvor + 1
        do l = 1, 6
          cdcls(l,jji) = cdcls(l,jj) 
        end do
        lviscstrp(jji) = lviscstrp(jj)

        do 80 ii = ifrsts(jj), ilasts(jj)
          nvor = nvor + 1
          if(nvor.gt.nvmax) then
            write(*,*) 'SDUPL: Vortex array overflow. Increase nvmax.'
            stop
          endif

          iii = nvor

          rv1(1,iii)  =  rv2(1,ii)
          rv1(2,iii)  = -rv2(2,ii) + yoff
          rv1(3,iii)  =  rv2(3,ii)

          rv2(1,iii)  =  rv1(1,ii)
          rv2(2,iii)  = -rv1(2,ii) + yoff
          rv2(3,iii)  =  rv1(3,ii)

          rv(1,iii)   =  rv(1,ii)
          rv(2,iii)   = -rv(2,ii) + yoff
          rv(3,iii)   =  rv(3,ii)

          rc(1,iii)   =  rc(1,ii)
          rc(2,iii)   = -rc(2,ii) + yoff
          rc(3,iii)   =  rc(3,ii)

          slopec(iii) = slopec(ii)
          slopev(iii) = slopev(ii)

          dxv(iii)    = dxv(ii)
          chordv(iii) = chordv(ii)

          rhov(iii) = rhov(ii)
          izimv(iii) = izimv(ii)

          isurfv(iii) = nni
          jstripv(iii) = nstrip
          lscompv(iii) = lscomp(nni)

          lvalbe(iii) = lvalbe(ii)

          do n = 1, ncontrol
ccc         rsgn = sign( 1.0 , reflsd(jj,n) )
            rsgn = reflsd(jj,n)
            dcontrol(iii,n) = -dcontrol(ii,n)*rsgn
          enddo

          ijetm(iii) = iii-1
   80   continue ! next vortex element
        ilasts(jji) = nvor
  100 continue ! next spanwise strip

      jlast(nni) = nstrip

c---- also duplicate upstream jet strips, if any
      do 200 ivs = 1, nvs
        jji = jfrst(nni) + ivs-1
        jj  = jfrst(nn)  + ivs-1

        if(ifrstu(jj) .gt. 0) then
c------- create image upstream jet strip elements
         ifrstu(jji) = nvor + 1

         do 210 ii = ifrstu(jj), ilastu(jj)
           nvor = nvor + 1
           if(nvor.gt.nvmax) then
             write(*,*) 'SDUPL: Vortex array overflow. Increase nvmax.'
             stop
           endif

           iii = nvor
         
           rv1(1,iii) =  rv2(1,ii)
           rv1(2,iii) = -rv2(2,ii) + yoff
           rv1(3,iii) =  rv2(3,ii)

           rv2(1,iii) =  rv1(1,ii)
           rv2(2,iii) = -rv1(2,ii) + yoff
           rv2(3,iii) =  rv1(3,ii)

           rv(1,iii) =  rv(1,ii)
           rv(2,iii) = -rv(2,ii) + yoff
           rv(3,iii) =  rv(3,ii)

           rc(1,iii) =  rc(1,ii)
           rc(2,iii) = -rc(2,ii) + yoff
           rc(3,iii) =  rc(3,ii)

           rs(1,iii) =  rs(1,ii)
           rs(2,iii) = -rs(2,ii) + yoff
           rs(3,iii) =  rs(3,ii)

           slopec(iii) = slopec(ii)
           slopev(iii) = slopev(ii)

           dxv(iii)    = dxv(ii)
           chordv(iii) = chordv(ii)

           rhov(iii) = rhov(ii)
           izimv(iii) = izimv(ii)

           isurfv(iii) = nni
           lscompv(iii) = lscomp(nni)
           jstripv(iii) = nstrip
 210     continue ! next jet-strip element
         ilastu(jji) = nvor

         rhos(jji) = rhos(jj)
         izims(jji) = izims(jj)

c------- index of element just upstream of each jet element
         iii = ifrstu(jji)
         ijetm(iii) = -999
         do iii = ifrstu(jji)+1, ilastu(jji)
           ijetm(iii) = iii-1
         enddo

         iii = ifrsts(jji)
         ijetm(iii) = ilastu(jji)

        else
c------- this strip has no upstream jet to duplicate
         ifrstu(jji) = 0
         ilastu(jji) = -1

         iii = ifrsts(jji)
         ijetm(iii) = -999

        endif

 200  continue ! next strip

c---- also duplicate wake jet strips, if any
      do 300 ivs = 1, nvs
        jji = jfrst(nni) + ivs-1
        jj  = jfrst(nn)  + ivs-1

        if(ifrstw(jj) .gt. 0) then
c------- create image jet strip elements
         ifrstw(jji) = nvor + 1

         do 310 ii = ifrstw(jj), ilastw(jj)
           nvor = nvor + 1
           if(nvor.gt.nvmax) then
             write(*,*) 'SDUPL: Vortex array overflow. Increase nvmax.'
             stop
           endif

           iii = nvor
         
           rv1(1,iii) =  rv2(1,ii)
           rv1(2,iii) = -rv2(2,ii) + yoff
           rv1(3,iii) =  rv2(3,ii)

           rv2(1,iii) =  rv1(1,ii)
           rv2(2,iii) = -rv1(2,ii) + yoff
           rv2(3,iii) =  rv1(3,ii)

           rv(1,iii) =  rv(1,ii)
           rv(2,iii) = -rv(2,ii) + yoff
           rv(3,iii) =  rv(3,ii)

           rc(1,iii) =  rc(1,ii)
           rc(2,iii) = -rc(2,ii) + yoff
           rc(3,iii) =  rc(3,ii)

           rs(1,iii) =  rs(1,ii)
           rs(2,iii) = -rs(2,ii) + yoff
           rs(3,iii) =  rs(3,ii)

           slopec(iii) = slopec(ii)
           slopev(iii) = slopev(ii)

           dxv(iii)    = dxv(ii)
           chordv(iii) = chordv(ii)

           rhov(iii) = rhov(ii)
           izimv(iii) = izimv(ii)

           isurfv(iii) = nni
           lscompv(iii) = lscomp(nni)
           jstripv(iii) = nstrip
 310     continue ! next jet-strip element
         ilastw(jji) = nvor

         rhos(jji) = rhos(jj)
         izims(jji) = izims(jj)

c------- index of element just upstream of each jet element
         iii = ifrstw(jji)
         ijetm(iii) = ilasts(jji)
         do iii = ifrstw(jji)+1, ilastw(jji)
           ijetm(iii) = iii-1
         enddo

        else
c------- this strip has no jet to duplicate
         ifrstw(jji) = 0
         ilastw(jji) = -1

        endif

 300  continue ! next strip

c---- image surface has same density factor      
      rhon(nni) = rhon(nn)

      nsurf = nsurf + 1

c---- return index of duplicate surface
      nn = nni

      return
      end ! sdupl



      subroutine bdupl(nn,ydupl,msg)
c-----------------------------------
c     Adds image of body nn,
c     reflected about y=ydupl.
c-----------------------------------
      include 'jvl.inc'
      character*(*) msg

      nni = nbody + 1
c      nni = min( nni , nbmax )
      if(nni.gt.nbmax) then
       write(*,*) 'BDUPL: Body array overflow. Increase nbmax.'
       stop
      endif

      klen = len(btitle(nn))
      do k = klen, 1, -1
        if(btitle(nn)(k:k) .ne. ' ') go to 6
      enddo
 6    btitle(nni) = btitle(nn)(1:k) // ' (' // msg // ')'
      write(*,*) ' '
      write(*,*) '  Building duplicate image-body: ',btitle(nni)

c---- same various logical flags
      lbwake(nni) = lbwake(nn)
      lbalbe(nni) = lbalbe(nn)
      lbload(nni) = lbload(nn)

      lfrst(nni) = nlbody + 1

      nvb = llast(nn) - lfrst(nn) + 1

      if(nlbody+nvb+1.gt.nlmax) then
       write(*,*) '*** MAKEBODY: Array overflow. Increase nlmax to',
     &             nlbody+nvb+1
       stop
      endif

      cbody(nni) = cbody(nn)
      abody(nni) = abody(nn)
      sbody(nni) = sbody(nn)
      vbody(nni) = vbody(nn)

      yoff = 2.0*ydupl

c---- set body nodes and radii
      do ivb = 1, nvb
        nlbody = nlbody + 1

        lli = lfrst(nni) + ivb-1
        ll  = lfrst(nn)  + ivb-1

        rl1(1,lli) =  rl1(1,ll)
        rl1(2,lli) = -rl1(2,ll) + yoff
        rl1(3,lli) =  rl1(3,ll)

        rl2(1,lli) =  rl2(1,ll)
        rl2(2,lli) = -rl2(2,ll) + yoff
        rl2(3,lli) =  rl2(3,ll)

        rl(1,lli) =  rl(1,ll)
        rl(2,lli) = -rl(2,ll) + yoff
        rl(3,lli) =  rl(3,ll)

        darl(lli) =  darl(ll)
        radl(lli) =  radl(ll)
        drbl(lli) =  drbl(ll)

        rhol(lli) = rhol(ll)

        lbcompl(lli) = -nni
      enddo
      llast(nni) = lli

c---- wake doublet segment
      nlbody = nlbody + 1

      lli = nlbody
      ll  = llast(nn) + 1

      rl1(1,lli) =  rl1(1,ll)
      rl1(2,lli) = -rl1(2,ll) + yoff
      rl1(3,lli) =  rl1(3,ll)
      rad1(lli)  =  rad1(ll)

      rl2(1,lli) =  rl2(1,ll)
      rl2(2,lli) = -rl2(2,ll) + yoff
      rl2(3,lli) =  rl2(3,ll)
      rad2(lli)  =  rad2(ll)

      darl(lli) =  darl(ll)
      radl(lli) =  radl(ll)
      drbl(lli) =  drbl(ll)

      rhol(lli) = rhol(ll)
      iziml(lli) = iziml(ll)
      lbcompl(lli) = -nni
      
c---- image body has same density factor      
      rhob(nni) = rhob(nn)
      izimb(nni) = izimb(nn)

      nbody = nbody + 1
      
      return
      end ! bdupl


      subroutine makejetu(istrip,dxle1,dxle2,dxle,nvcu,cspaceu)
c----------------------------------------------------------------
c     Sets up jet elements ahead of surface strip istrip
c----------------------------------------------------------------
      include 'jvl.inc'

      parameter (kcmax=50)
      real xptj(kcmax), xcpj(kcmax), xvrj(kcmax), xsrj(kcmax)

      real rpr1(3), rpr2(3), rpr(3)

      data dxdf / 0.001 /

      if(nvcu.gt.kcmax) then
       write(*,*) '* MAKEJETU: Array overflow.  Increase kcmax to', nvcu
       nvc = kcmax
      endif

      if(nvor+nvcu .gt. nvmax) then
       write(*,*) '* MAKEJETU: Array overflow.  Increase nvmax to',
     &   nvor+nvcu
       stop
      endif

c====================================================
c---- make sure all prop dx offsets are sufficiently negative
      dxdmin = -dxdf*chord(istrip)
      if(dxdstrp1(istrip) .gt. dxdmin .and.
     &   dxdstrp2(istrip) .gt. dxdmin .and.
     &   dxdstrp(istrip)  .gt. dxdmin      ) then
c----- don't add any upstream jet elements to this strip
       ifrstu(istrip) = 0
       ifrstu(istrip) = -1
       return
      endif

c====================================================
      isurf = isurfs(istrip)

      ifrstu(istrip) = nvor + 1
      ilastu(istrip) = nvor + nvcu

      i1 = ifrstu(istrip)
      in = ilastu(istrip)

c---- index of element just upstream of each jet element
      ijetm(i1) = -999
      do i = i1+1, in
        ijetm(i) = i-1
      enddo

      rpr1(1) = rle1(1,istrip) + dxdstrp1(istrip)
      rpr1(2) = rle1(2,istrip)
      rpr1(3) = rle1(3,istrip)

      rpr2(1) = rle2(1,istrip) + dxdstrp2(istrip)
      rpr2(2) = rle2(2,istrip)
      rpr2(3) = rle2(3,istrip)

      rpr(1) = rle(1,istrip) + dxdstrp(istrip)
      rpr(2) = rle(2,istrip)
      rpr(3) = rle(3,istrip)

      clafc = 1.0
      call cspacer(nvcu,cspaceu,clafc, xptj,xvrj,xsrj,xcpj)

      ulen1 = -dxdstrp1(istrip)
      ulen2 = -dxdstrp2(istrip)
      ulen  = -dxdstrp(istrip)
      dxc = xptj(nvcu+1) - xcpj(nvcu)

c---- x-stretching lengths for jet strip
      al1 = ulen1*(ulen1-dxle1)/dxle1 * dxc/(1.0-dxc)
      al2 = ulen2*(ulen2-dxle2)/dxle2 * dxc/(1.0-dxc)
      al  = ulen *(ulen -dxle )/dxle  * dxc/(1.0-dxc)

c---- go over upstream jet vortex elements ahead of LE
      do ivc = 1, nvcu
        nvor = nvor + 1
        i = nvor

        xvr = xvrj(ivc)
        xcp = xcpj(ivc)
        xsr = xsrj(ivc)

c------ new spacing formulation
        xvr = xsrj(ivc)
        xcp = xptj(ivc+1)

        rv1(1,i) = rpr1(1) + ulen1*al1*xvr/(ulen1*(1.0-xvr) + al1*xvr)
        rv1(2,i) = rpr1(2)
        rv1(3,i) = rpr1(3)

        rv2(1,i) = rpr2(1) + ulen2*al2*xvr/(ulen2*(1.0-xvr) + al2*xvr)
        rv2(2,i) = rpr2(2)
        rv2(3,i) = rpr2(3)

        rv(1,i) = rpr(1) + ulen*al*xvr/(ulen*(1.0-xvr) + al*xvr)
        rv(2,i) = rpr(2)
        rv(3,i) = rpr(3)

        rc(1,i) = rpr(1) + ulen*al*xcp/(ulen*(1.0-xcp) + al*xcp)
        rc(2,i) = rpr(2)
        rc(3,i) = rpr(3)

        rs(1,i) = rpr(1) + ulen*al*xsr/(ulen*(1.0-xsr) + al*xsr)
        rs(2,i) = rpr(2)
        rs(3,i) = rpr(3)

        slopec(i) = 0.0
        slopev(i) = 0.0

c------ streamwise length of element
        dxv(i) = ulen*al
     &    * (  xptj(ivc+1)/(al*(1.0-xptj(ivc+1)) + ulen*xptj(ivc+1))
     &       - xptj(ivc  )/(al*(1.0-xptj(ivc  )) + ulen*xptj(ivc  )) )

        chordv(i) = chord(istrip)
        isurfv(i) = isurf
        jstripv(i) = istrip
        lscompv(i) = lscomp(isurf)

        rhov(i) = rhos(istrip)
        izimv(i) = izims(istrip)
      enddo

      return
      end ! makejetu


      subroutine makejetw(istrip,dxte1,dxte2,dxte,nvcw,cspacew)
c----------------------------------------------------------------
c     Sets up jet elements trailing from surface strip istrip
c----------------------------------------------------------------
      include 'jvl.inc'

      parameter (kcmax=50)
      real xptj(kcmax), xcpj(kcmax), xvrj(kcmax), xsrj(kcmax)

      real rte1(3), rte2(3), rte(3)

      if(nvcw.gt.kcmax) then
       write(*,*) '* MAKEJETW: Array overflow.  Increase kcmax to', nvcw
       nvc = kcmax
      endif

      if(nvor+nvcw .gt. nvmax) then
       write(*,*) '* MAKEJETW: Array overflow.  Increase nvmax to',
     &   nvor+nvcw
       stop
      endif

c====================================================
      isurf = isurfs(istrip)

      ifrstw(istrip) = nvor + 1
      ilastw(istrip) = nvor + nvcw

      i1 = ifrstw(istrip)
      in = ilastw(istrip)

c---- index of element just upstream of each jet element
      ijetm(i1) = ilasts(istrip)
      do i = i1+1, in
        ijetm(i) = i-1
      enddo

      rte1(1) = rle1(1,istrip) + chord1(istrip)
      rte1(2) = rle1(2,istrip)
      rte1(3) = rle1(3,istrip)

      rte2(1) = rle2(1,istrip) + chord2(istrip)
      rte2(2) = rle2(2,istrip)
      rte2(3) = rle2(3,istrip)

      rte(1) = rle(1,istrip) + chord(istrip)
      rte(2) = rle(2,istrip)
      rte(3) = rle(3,istrip)

      clafc = 1.0
      call cspacer(nvcw,cspacew,clafc, xptj,xvrj,xsrj,xcpj)

      wlen1 = 2.0*bref
      wlen2 = 2.0*bref
      wlen  = 2.0*bref
      dxc = xvrj(1) - xptj(1)

c---- x-stretching lengths for jet strip
      al1 = wlen1*dxte1/(wlen1-dxte1) * (1.0-dxc)/dxc
      al2 = wlen2*dxte2/(wlen2-dxte2) * (1.0-dxc)/dxc
      al  = wlen *dxte /(wlen -dxte)  * (1.0-dxc)/dxc

c---- go over jet vortex elements trailing from surface strip
      do ivc = 1, nvcw
        nvor = nvor + 1
        i = nvor

        xvr = xvrj(ivc)
        xcp = xcpj(ivc)
        xsr = xsrj(ivc)

c------ new spacing formulation
        xvr = xsrj(ivc)
        xcp = xptj(ivc+1)

        rv1(1,i) = rte1(1) + wlen1*al1*xvr/(wlen1*(1.0-xvr) + al1*xvr)
        rv1(2,i) = rte1(2)
        rv1(3,i) = rte1(3)

        rv2(1,i) = rte2(1) + wlen2*al2*xvr/(wlen2*(1.0-xvr) + al2*xvr)
        rv2(2,i) = rte2(2)
        rv2(3,i) = rte2(3)

        rv(1,i) = rte(1) + wlen*al*xvr/(wlen*(1.0-xvr) + al*xvr)
        rv(2,i) = rte(2)
        rv(3,i) = rte(3)

        rc(1,i) = rte(1) + wlen*al*xcp/(wlen*(1.0-xcp) + al*xcp)
        rc(2,i) = rte(2)
        rc(3,i) = rte(3)

        rs(1,i) = rte(1) + wlen*al*xsr/(wlen*(1.0-xsr) + al*xsr)
        rs(2,i) = rte(2)
        rs(3,i) = rte(3)

        slopec(i) = 0.0
        slopev(i) = 0.0

c------ streamwise length of element
        dxv(i) = wlen*al
     &    * (  xptj(ivc+1)/(wlen*(1.0-xptj(ivc+1)) + al*xptj(ivc+1))
     &       - xptj(ivc  )/(wlen*(1.0-xptj(ivc  )) + al*xptj(ivc  )) )

        chordv(i) = chord(istrip)
        isurfv(i) = isurf
        jstripv(i) = istrip
        lscompv(i) = lscomp(isurf)

        rhov(i) = rhos(istrip)
        izimv(i) = izims(istrip)
      enddo

      return
      end ! makejetw


      subroutine encalc
c-------------------------------------------------------------------
c     Calculate normal vectors for strips, h.v.'s, and c.p.'s
c     incorporates surface deflections.
c
c   Inputs:
c            nvor      number of vortices
c            x1        coordinates of endpoint #1 of the vortices
c            x2        coordinates of endpoint #2 of the vortices
c            slopev    slope at bound vortices
c            slopec    slope at control points
c            nstrip    number of strips
c            ifrst    index of first element in strip
c            ainc      angle of incidence of strip
c            ldes      include design-variable deflections if true
c
c   Outputs: enc(.)        normal vector at control point
c            env(.)        normal vector at bound vortices
c            ensy, ensz    strip normal vector (ensx=0)
c            lstripoff     non-used strip (t) (below z=zsym)
c-------------------------------------------------------------------
      include 'jvl.inc'

      real ep(3), eq(3), es(3), eb(3), ec(3), ecxb(3)
      real ec_g(3,ndmax), ecxb_g(3)

c---- normal vectors at control points and bound vortex midpoints
      do 100 j = 1, nstrip
c------ normal vector for the strip (normal to x axis)
        i = ifrsts(j)
        dxle =  rv2(1,i)-rv1(1,i)
        dyle =  rv2(2,i)-rv1(2,i)
        dzle =  rv2(3,i)-rv1(3,i)
c       axle = (rv2(1,i)+rv1(1,i))*0.5
c       ayle = (rv2(2,i)+rv1(2,i))*0.5
c       azle = (rv2(3,i)+rv1(3,i))*0.5
        axle = rv(1,i)
        ayle = rv(2,i)
        azle = rv(3,i)

        i = ilasts(j)
        dxte =  rv2(1,i)-rv1(1,i)
        dyte =  rv2(2,i)-rv1(2,i)
        dzte =  rv2(3,i)-rv1(3,i)
c       axte = (rv2(1,i)+rv1(1,i))*0.5
c       ayte = (rv2(2,i)+rv1(2,i))*0.5
c       azte = (rv2(3,i)+rv1(3,i))*0.5
        axte = rv(1,i)
        ayte = rv(2,i)
        azte = rv(3,i)

        dxt = (1.0-saxfr)*dxle + saxfr*dxte
        dyt = (1.0-saxfr)*dyle + saxfr*dyte
        dzt = (1.0-saxfr)*dzle + saxfr*dzte

        ess(1,j) =  dxt/sqrt(dxt*dxt + dyt*dyt + dzt*dzt)
        ess(2,j) =  dyt/sqrt(dxt*dxt + dyt*dyt + dzt*dzt)
        ess(3,j) =  dzt/sqrt(dxt*dxt + dyt*dyt + dzt*dzt)

        ensy(j) = -dzt/sqrt(dyt*dyt + dzt*dzt)
        ensz(j) =  dyt/sqrt(dyt*dyt + dzt*dzt)

        xyzrefs(1,j) = (1.0-saxfr)*axle + saxfr*axte
        xyzrefs(2,j) = (1.0-saxfr)*ayle + saxfr*ayte
        xyzrefs(3,j) = (1.0-saxfr)*azle + saxfr*azte

        es(1) = 0.
        es(2) = ensy(j)
        es(3) = ensz(j)

        lstripoff(j) = .false.

c------ go over surface elements in this strip
        do 10 i = ifrsts(j), ilasts(j)
          do n = 1, ncontrol
            env_d(1,i,n) = 0.
            env_d(2,i,n) = 0.
            env_d(3,i,n) = 0.
            enc_d(1,i,n) = 0.
            enc_d(2,i,n) = 0.
            enc_d(3,i,n) = 0.
          enddo

          do n = 1, ndesign
            env_g(1,i,n) = 0.
            env_g(2,i,n) = 0.
            env_g(3,i,n) = 0.
            enc_g(1,i,n) = 0.
            enc_g(2,i,n) = 0.
            enc_g(3,i,n) = 0.
          enddo

c...define unit vector along bound leg
          dxb = rv2(1,i)-rv1(1,i)
          dyb = rv2(2,i)-rv1(2,i)
          dzb = rv2(3,i)-rv1(3,i)
          emag = sqrt(dxb**2 + dyb**2 + dzb**2)
          eb(1) = dxb/emag
          eb(2) = dyb/emag
          eb(3) = dzb/emag

c...define direction of normal vector at control point 
c   the yz projection of the normal vector matches the camber slope
c   + section local incidence in the yz defining plane for the section
          ang = ainc(j) - atan(slopec(i))
cc          if(ldes) then
c--------- add design-variable contribution to angle
c           do n = 1, ndesign
c             ang = ang + ainc_g(j,n)*deldes(n)
c           enddo
cc          endif

          sinc = sin(ang)
          cosc = cos(ang)
          ec(1) =  cosc
          ec(2) = -sinc*es(2)
          ec(3) = -sinc*es(3)
          do n = 1, ndesign
            ec_g(1,n) = -sinc      *ainc_g(j,n)
            ec_g(2,n) = -cosc*es(2)*ainc_g(j,n)
            ec_g(3,n) = -cosc*es(3)*ainc_g(j,n)
          enddo

c...normal vector is perpendicular to camberline vector and to the bound leg
          call cross(ec,eb,ecxb)
          emag = sqrt(ecxb(1)**2 + ecxb(2)**2 + ecxb(3)**2)
          if(emag.ne.0.0) then
            enc(1,i) = ecxb(1)/emag
            enc(2,i) = ecxb(2)/emag
            enc(3,i) = ecxb(3)/emag
            do n = 1, ndesign
              call cross(ec_g(1,n),eb,ecxb_g)
              emag_g = enc(1,i)*ecxb_g(1)
     &               + enc(2,i)*ecxb_g(2)
     &               + enc(3,i)*ecxb_g(3)
              enc_g(1,i,n) = (ecxb_g(1) - enc(1,i)*emag_g)/emag
              enc_g(2,i,n) = (ecxb_g(2) - enc(2,i)*emag_g)/emag
              enc_g(3,i,n) = (ecxb_g(3) - enc(3,i)*emag_g)/emag
            enddo
          else
            enc(1,i) = es(1)
            enc(2,i) = es(2)
            enc(3,i) = es(3)
          endif

c...define direction of normal vector at vortex mid-point. 
c   the yz projection of the normal vector matches the camber slope
c   + section local incidence in the yz defining plane for the section
          ang = ainc(j) - atan(slopev(i)) 
cc          if(ldes) then
c--------- add design-variable contribution to angle
c           do n = 1, ndesign
c             ang = ang + ainc_g(j,n)*deldes(n)
c           enddo
cc          endif

          sinc = sin(ang)
          cosc = cos(ang)
          ec(1) =  cosc
          ec(2) = -sinc*es(2)
          ec(3) = -sinc*es(3)
          do n = 1, ndesign
            ec_g(1,n) = -sinc      *ainc_g(j,n)
            ec_g(2,n) = -cosc*es(2)*ainc_g(j,n)
            ec_g(3,n) = -cosc*es(3)*ainc_g(j,n)
          enddo

c-------- normal vector is perpendicular to camberline vector 
c-         and to the bound leg
          call cross(ec,eb,ecxb)
          emag = sqrt(ecxb(1)**2 + ecxb(2)**2 + ecxb(3)**2)
          if(emag.ne.0.0) then
            env(1,i) = ecxb(1)/emag
            env(2,i) = ecxb(2)/emag
            env(3,i) = ecxb(3)/emag
            do n = 1, ndesign
              call cross(ec_g(1,n),eb,ecxb_g)
              emag_g = enc(1,i)*ecxb_g(1)
     &               + enc(2,i)*ecxb_g(2)
     &               + enc(3,i)*ecxb_g(3)
              env_g(1,i,n) = (ecxb_g(1) - env(1,i)*emag_g)/emag
              env_g(2,i,n) = (ecxb_g(2) - env(2,i)*emag_g)/emag
              env_g(3,i,n) = (ecxb_g(3) - env(3,i)*emag_g)/emag
            enddo
          else
            env(1,i) = es(1)
            env(2,i) = es(2)
            env(3,i) = es(3)
          endif

c=======================================================
c-------- rotate normal vectors for control surface
          do 15 n = 1, ncontrol

c---------- skip everything if this element is unaffected by control variable n
            if(dcontrol(i,n).eq.0.0) go to 15

            ang     = dtr*dcontrol(i,n)*delcon(n)
            ang_ddc = dtr*dcontrol(i,n)

            cosd = cos(ang)
            sind = sin(ang)

c---------- ep = normal-vector component perpendicular to hinge line
            endot = dot(enc(1,i),vhinge(1,j,n))
            ep(1) = enc(1,i) - endot*vhinge(1,j,n)
            ep(2) = enc(2,i) - endot*vhinge(2,j,n)
            ep(3) = enc(3,i) - endot*vhinge(3,j,n)
c---------- eq = unit vector perpendicular to both ep and hinge line
            call cross(vhinge(1,j,n),ep,eq)

c---------- rotated vector would consist of sin,cos parts from ep and eq,
c-          with hinge-parallel component endot restored 
cc          enc(1,i) = ep(1)*cosd + eq(1)*sind + endot*vhinge(1,j,n)
cc          enc(2,i) = ep(2)*cosd + eq(2)*sind + endot*vhinge(2,j,n)
cc          enc(3,i) = ep(3)*cosd + eq(3)*sind + endot*vhinge(3,j,n)

c---------- linearize about zero deflection (cosd=1, sind=0)
            enc_d(1,i,n) = enc_d(1,i,n) + eq(1)*ang_ddc
            enc_d(2,i,n) = enc_d(2,i,n) + eq(2)*ang_ddc
            enc_d(3,i,n) = enc_d(3,i,n) + eq(3)*ang_ddc

c
c---------- repeat for env vector

c---------- ep = normal-vector component perpendicular to hinge line
            endot = dot(env(1,i),vhinge(1,j,n))
            ep(1) = env(1,i) - endot*vhinge(1,j,n)
            ep(2) = env(2,i) - endot*vhinge(2,j,n)
            ep(3) = env(3,i) - endot*vhinge(3,j,n)
c---------- eq = unit vector perpendicular to both ep and hinge line
            call cross(vhinge(1,j,n),ep,eq)

c---------- rotated vector would consist of sin,cos parts from ep and eq,
c-          with hinge-parallel component endot restored 
cc          env(1,i) = ep(1)*cosd + eq(1)*sind + endot*vhinge(1,j,n)
cc          env(2,i) = ep(2)*cosd + eq(2)*sind + endot*vhinge(2,j,n)
cc          env(3,i) = ep(3)*cosd + eq(3)*sind + endot*vhinge(3,j,n)

c---------- linearize about zero deflection (cosd=1, sind=0)
            env_d(1,i,n) = env_d(1,i,n) + eq(1)*ang_ddc
            env_d(2,i,n) = env_d(2,i,n) + eq(2)*ang_ddc
            env_d(3,i,n) = env_d(3,i,n) + eq(3)*ang_ddc
 15       continue ! next control variable n
 10     continue ! next surface element i

c------ go over upstream jet elements in this strip, if any
        do 20 i = ifrstu(j), ilastu(j)
          do n = 1, ncontrol
            env_d(1,i,n) = 0.
            env_d(2,i,n) = 0.
            env_d(3,i,n) = 0.
            enc_d(1,i,n) = 0.
            enc_d(2,i,n) = 0.
            enc_d(3,i,n) = 0.
          enddo

          do n = 1, ndesign
            env_g(1,i,n) = 0.
            env_g(2,i,n) = 0.
            env_g(3,i,n) = 0.
            enc_g(1,i,n) = 0.
            enc_g(2,i,n) = 0.
            enc_g(3,i,n) = 0.
          enddo

c...define unit vector along bound leg
          dxb = rv2(1,i)-rv1(1,i)
          dyb = rv2(2,i)-rv1(2,i)
          dzb = rv2(3,i)-rv1(3,i)
          emag = sqrt(dxb**2 + dyb**2 + dzb**2)
          eb(1) = dxb/emag
          eb(2) = dyb/emag
          eb(3) = dzb/emag

          ec(1) = 1.0
          ec(2) = 0.
          ec(3) = 0.

          call cross(ec,eb,ecxb)
          emag = sqrt(ecxb(1)**2 + ecxb(2)**2 + ecxb(3)**2)
          if(emag.ne.0.0) then
            enc(1,i) = ecxb(1)/emag
            enc(2,i) = ecxb(2)/emag
            enc(3,i) = ecxb(3)/emag
          else
            enc(1,i) = es(1)
            enc(2,i) = es(2)
            enc(3,i) = es(3)
          endif
c
c...define direction of normal vector at vortex mid-point. 
c   the yz projection of the normal vector matches the camber slope (0)
c   + section local incidence in the yz defining plane for the section
          ec(1) = 1.0
          ec(2) = 0.
          ec(3) = 0.

c-------- normal vector is perpendicular to camberline vector 
c-         and to the bound leg
          call cross(ec,eb,ecxb)
          emag = sqrt(ecxb(1)**2 + ecxb(2)**2 + ecxb(3)**2)
          if(emag.ne.0.0) then
            env(1,i) = ecxb(1)/emag
            env(2,i) = ecxb(2)/emag
            env(3,i) = ecxb(3)/emag
          else
            env(1,i) = es(1)
            env(2,i) = es(2)
            env(3,i) = es(3)
          endif
 20     continue ! next upstream jet element i

c------ go over jet elements in this strip, if any
        do 30 i = ifrstw(j), ilastw(j)
          do n = 1, ncontrol
            env_d(1,i,n) = 0.
            env_d(2,i,n) = 0.
            env_d(3,i,n) = 0.
            enc_d(1,i,n) = 0.
            enc_d(2,i,n) = 0.
            enc_d(3,i,n) = 0.
          enddo

          do n = 1, ndesign
            env_g(1,i,n) = 0.
            env_g(2,i,n) = 0.
            env_g(3,i,n) = 0.
            enc_g(1,i,n) = 0.
            enc_g(2,i,n) = 0.
            enc_g(3,i,n) = 0.
          enddo

c...define unit vector along bound leg
          dxb = rv2(1,i)-rv1(1,i)
          dyb = rv2(2,i)-rv1(2,i)
          dzb = rv2(3,i)-rv1(3,i)
          emag = sqrt(dxb**2 + dyb**2 + dzb**2)
          eb(1) = dxb/emag
          eb(2) = dyb/emag
          eb(3) = dzb/emag

          ec(1) = 1.0
          ec(2) = 0.
          ec(3) = 0.

          call cross(ec,eb,ecxb)
          emag = sqrt(ecxb(1)**2 + ecxb(2)**2 + ecxb(3)**2)
          if(emag.ne.0.0) then
            enc(1,i) = ecxb(1)/emag
            enc(2,i) = ecxb(2)/emag
            enc(3,i) = ecxb(3)/emag
          else
            enc(1,i) = es(1)
            enc(2,i) = es(2)
            enc(3,i) = es(3)
          endif

c
c...define direction of normal vector at vortex mid-point. 
c   the yz projection of the normal vector matches the camber slope (0)
c   + section local incidence in the yz defining plane for the section
          ec(1) = 1.0
          ec(2) = 0.
          ec(3) = 0.

c-------- normal vector is perpendicular to camberline vector 
c-         and to the bound leg
          call cross(ec,eb,ecxb)
          emag = sqrt(ecxb(1)**2 + ecxb(2)**2 + ecxb(3)**2)
          if(emag.ne.0.0) then
            env(1,i) = ecxb(1)/emag
            env(2,i) = ecxb(2)/emag
            env(3,i) = ecxb(3)/emag
          else
            env(1,i) = es(1)
            env(2,i) = es(2)
            env(3,i) = es(3)
          endif
 30     continue ! next jet element i

 100  continue ! next strip j

      lenc = .true.

      return
      end ! encalc
