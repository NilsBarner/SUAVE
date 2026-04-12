
C***********************************************************************
C    Module:  joutput.f
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

      subroutine outtot(lun)
c---------------------------------------------------------------
c     Writes formatted overall forces and moments
c     for current solution to logical unit lun
c---------------------------------------------------------------
      include 'jvl.inc'
      character*50 satype
      real msax(3), rotsax(3)
      real lti, ltj, ltv
      real drm(3), mtm(3)

 1000 format (a)

      if (lun.eq.0) return

      qs  = 0.5*sref
      qsb = 0.5*sref*bref
      qsc = 0.5*sref*cref

      drm(1) = xyzref(1) - xyzmom(1) 
      drm(2) = xyzref(2) - xyzmom(2) 
      drm(3) = xyzref(3) - xyzmom(3) 
      mtm(1) = mt(1) + drm(2)*ft(3) - drm(3)*ft(2)
      mtm(2) = mt(2) + drm(3)*ft(1) - drm(1)*ft(3)
      mtm(3) = mt(3) + drm(1)*ft(2) - drm(2)*ft(1)

      call getsa(lnasa_sa,satype,dir)

      ca = cos(alfa)
      sa = sin(alfa)

c---- set rates in stability axes
      rotsax(1) = wbar(1)*ca + wbar(3)*sa
      rotsax(2) = wbar(2)                
      rotsax(3) = wbar(3)*ca - wbar(1)*sa

c---- set moments in stability axes
      msax(1) = mtm(1)*ca + mtm(3)*sa
      msax(2) = mtm(2) 
      msax(3) = mtm(3)*ca - mtm(1)*sa  

c---- inviscid near-field forces
      dti = ca*fti(1) + sa*fti(3)
      yti =    fti(2)            
      lti = ca*fti(3) - sa*fti(1)

c---- jet near-field forces
      dtj = ca*ftj(1) + sa*ftj(3)
      ytj =    ftj(2)            
      ltj = ca*ftj(3) - sa*ftj(1)

c---- viscous near-field forces
      dtv = ca*ftv(1) + sa*ftv(3)
      ytv =    ftv(2)            
      ltv = ca*ftv(3) - sa*ftv(1)

c---- dump it
      write(lun,200)
      write(lun,201) 
      write(lun,202) title(1:60),nsurf,nstrip,nvor
      if(iysym.gt.0) write (*,2034) ysym
      if(iysym.lt.0) write (*,2035) ysym
      if(izsym.gt.0) write (*,2036) zsym
      if(izsym.lt.0) write (*,2037) zsym

      write(lun,204) sref,cref,bref, 
     &               xyzref(1), xyzref(2), xyzref(3),
     &               xyzmom(1), xyzmom(2), xyzmom(3)

      write(lun,205) satype

      write(lun,218) rtitle(irun)
      write(lun,220) 
     &  alfa/dtr, dir*wbar(1)*bref/2.0, dir*rotsax(1)*bref/2.0,
     &  beta/dtr,     wbar(2)*cref/2.0,
     &  amach   , dir*wbar(3)*bref/2.0, dir*rotsax(3)*bref/2.0

      write (lun,226)
     & dir*ft(1)/qs, dir*mtm(1)/qsb, dir*msax(1)/qsb,
     &     ft(2)/qs,     mtm(2)/qsc,
     & dir*ft(3)/qs, dir*mtm(3)/qsb, dir*msax(3)/qsb,
     & ltot/qs,  lff/qs,
     & ytot/qs,  yff/qs,
     & dtot/qs,  dff/qs, spanef,
     &  lti/qs, lffi/qs,
     &  ltj/qs, lffj/qs,
     &  dti/qs, dffi/qs, spanev,
     &  dtj/qs, dffj/qs,
     &  dtv/qs, dffv/qs

      if(nvarjet .gt. 0) then
       dcq = dqt/sref
       dcj = djt/qs
       dce = 0.5*det/qs
       ct = dcj - 2.0*dcq
       cp = dce - dcq
       write(lun,222) dcq, dcj, dce, ct, cp
      endif

 226  format (
     & /2x,'CXtot =',f10.5,5x,'Cltot =',f10.5,5x,'Cl''tot =',f10.5
     & /2x,'CYtot =',f10.5,5x,'Cmtot =',f10.5 
     & /2x,'CZtot =',f10.5,5x,'Cntot =',f10.5,5x,'Cn''tot =',f10.5
     &//2x,'CLtot =',f10.5,5x,'CLff  =',f10.5,4x,
     & /2x,'CYtot =',f10.5,5x,'CYff  =',f10.5,4x,
     & /2x,'CDtot =',f10.5,5x,'CDff  =',f10.5,4x,'    e =',f10.5,
     &//2x,'CLcir =',f10.5,5x,'CLffc =',f10.5,4x,
     & /2x,'CLjet =',f10.5,5x,'CLffj =',f10.5,4x,
     &//2x,'CDind =',f10.5,5x,'CDffi =',f10.5,4x,'e_vec =',f10.5,
     & /2x,'CDjet =',f10.5,5x,'CDffj =',f10.5,4x,
     & /2x,'CDvis =',f10.5,5x,'CDffv =',f10.5,4x )

      write(lun,*)
      do k = 1, ncontrol
        if(k.le.9) then
         write(lun,231) dname(k),'d',k, delcon(k)
        else
         write(lun,232) dname(k),'d',k, delcon(k)
        endif
      enddo

      if(ndesign.gt.0) write(lun,*)
      do k = 1, ndesign
        if(k.le.9) then
         write(lun,231) gname(k),'g',k, deldes(k)
        else
         write(lun,232) gname(k),'g',k, deldes(k)
        endif
      enddo

      if(nvarjet.gt.0) write(lun,*)
      do k = 1, nvarjet
        if(k.le.9) then
         write(lun,231) jname(k),'j',k, deljet(k)
        else
         write(lun,232) jname(k),'j',k, deljet(k)
        endif
      enddo

      write(lun,200)
 200  format(1x,
     &'---------------------------------------------------------------')
 201  format(' Jet Vortex Lattice Output -- Total Forces')
 202  format(/' Configuration: ',a
     &       /5x,'# Surfaces =',i4
     &       /5x,'# Strips   =',i4
     &       /5x,'# Vortices =',i4)

 2034 format(' Y symmetry: Wall plane   at Ysym =',g10.4)
 2035 format(' Y symmetry: Free surface at Ysym =',g10.4)
 2036 format(' Z symmetry: Ground plane at Zsym =',g10.4)
 2037 format(' Z symmetry: Free surface at Zsym =',g10.4)

 204  format(/2x, 'Sref =',g13.5,3x,'Cref =',g13.5,3x,'Bref =',g13.5
     &       /2x, 'Xref =',g13.5,3x,'Yref =',g13.5,3x,'Zref =',g13.5
     &       /2x, 'Xmom =',g13.5,3x,'Ymom =',g13.5,3x,'Zmom =',g13.5 )

 205  format(/1x, a)
 218  format(/' Run case: ', a)
 220  format(
     & /2x,'alpha =',f10.5,5x,'pb/2V =',f10.5,5x,'p''b/2V =',f10.5
     & /2x,'beta  =',f10.5,5x,'qc/2V =',f10.5
     & /2x,'Mach  =',f10.3,5x,'rb/2V =',f10.5,5x,'r''b/2V =',f10.5)
 222  format(
     & /1x,'DCQtot =',g12.5,3x,'DCJtot =',g12.5,3x,'DCEtot =',g12.5
     & /1x,'        ',12x  ,3x,' CTtot =',g12.5,3x,' CPtot =',g12.5)

 231  format(1x,a,1x,a,i1,'  =',g12.5)
 232  format(1x,a,1x,a,i2, ' =',g12.5)

 240  format (/)

      return
      end ! outtot


      subroutine outsurf(lun)
c---------------------------------------------------------------
c     Writes formatted surface forces and moments
c     for current solution to logical unit lun
c---------------------------------------------------------------
      include 'jvl.inc'
      character*50 satype
      real r(3)
      real fn(3), mn(3)
      real ln
      real drm(3)

      integer icrs(3), jcrs(3)
      data icrs / 2, 3, 1 / , jcrs / 3, 1, 2 /

      call getsa(lnasa_sa,satype,dir)

      sa = sin(alfa)
      ca = cos(alfa)

      qs  = 0.5*sref
      qsb = 0.5*sref*bref
      qsc = 0.5*sref*cref

      drm(1) = xyzref(1) - xyzmom(1) 
      drm(2) = xyzref(2) - xyzmom(2) 
      drm(3) = xyzref(3) - xyzmom(3) 

c========================================================================
c---- force components from each surface
      write(lun,200)
 200  format(1x,
     &'---------------------------------------------------------------')
      write (lun,205) satype,
     &                sref,cref,bref, 
     &                xyzmom(1), xyzmom(2), xyzmom(3)
 205  format ( ' Surface Forces (referred to Sref,Cref,Bref',
     &                                   ' about Xmom,Ymom,Zmom)',
     &        /1x,a //
     &      5x,'Sref =',g13.5,3x,'Cref =',g13.5,3x,'Bref =',g13.5/
     &      5x,'Xmom =',g13.5,3x,'Ymom =',g13.5,3x,'Zmom =',g13.5 )

      write (lun,208)
 208  format(/' n',4x,'Area',9x,'CL',8x,'CD',8x,'Cm',
     &                       8x,'CY',8x,'Cn',8x,'Cl',8x,'CDi',7x,'CDv')
      do is = 1, nsurf
        do k = 1, 3
          fn(k) = fni(k,is) + fnv(k,is) + fnj(k,is)
          mn(k) = mni(k,is) + mnv(k,is) + mnj(k,is)
        enddo

        mn(1) = mn(1) + drm(2)*fn(3) - drm(3)*fn(2)
        mn(2) = mn(2) + drm(3)*fn(1) - drm(1)*fn(3)
        mn(3) = mn(3) + drm(1)*fn(2) - drm(2)*fn(1)

        dn = ca*fn(1) + sa*fn(3)
        yn =    fn(2)
        ln = ca*fn(3) - sa*fn(1)

        dni = ca*fni(1,is) + sa*fni(3,is)
        dnv = ca*fnv(1,is) + sa*fnv(3,is)

        call bstrip(stitle(is),nt)
        write (lun,211) is,ssurf(is),
     &      ln    / qs,
     &      dn    / qs,
     &      mn(2) / qsc,
     &      yn    / qs,
     &  dir*mn(3) / qsb,
     &  dir*mn(1) / qsb,
     &      dni   / qs,
     &      dnv   / qs,
     &      stitle(is)(1:nt)
      end do

cc      write(lun,212)
  211 format (i2,1x,g12.5,8f10.5,3x,a)
  212 format (/)

c========================================================================
c---- surface forces normalized by local reference quantities
      write (lun,220) 
  220 format (/' Surface Forces (referred to Ssurf, Cave ',
     &         'about root LE on hinge axis)'//
     &        2x,' n',2x,'Ssurf',8x,'Cave',
     &        8x,'cl_surf',4x,'cd_surf',4x,'cdv_surf',2x,'cmLE_surf')
  221 format (2x,i2,g12.5,g12.5,4f11.6,3x,a)

      do is = 1, nsurf
c------ set total surface force and moment
        do k = 1, 3
          fn(k) = fni(k,is) + fnv(k,is) + fnj(k,is)
          mn(k) = mni(k,is) + mnv(k,is) + mnj(k,is)
        enddo

c------ reference point for surface le (hinge) moments
c        defined by surface hinge vector direction thru first strip le point
        if(idups(is).ge.0) then
         r(1) = xyzref(1) - rle1(1,jfrst(is))
         r(2) = xyzref(2) - rle1(2,jfrst(is))
         r(3) = xyzref(3) - rle1(3,jfrst(is))
        else
         r(1) = xyzref(1) - rle2(1,jfrst(is))
         r(2) = xyzref(2) - rle2(2,jfrst(is))
         r(3) = xyzref(3) - rle2(3,jfrst(is))
        endif

c------ set moment about point r(.)
        mn(1) = mn(1) + r(2)*fn(3) - r(3)*fn(2)
        mn(2) = mn(2) + r(3)*fn(1) - r(1)*fn(3)
        mn(3) = mn(3) + r(1)*fn(2) - r(2)*fn(1)

c------ rotate forces into stability axes
        dn = ca*fn(1) + sa*fn(3)
        yn =    fn(2)
        ln = ca*fn(3) - sa*fn(1)

        dni = ca*fni(1,is) + sa*fni(3,is)
        dnv = ca*fnv(1,is) + sa*fnv(3,is)

c------ surface hinge moment defined about hinge vector 
        mle = dot(mn,ess(1,is))
        
c       rho = 1.0
        rho = rhon(is)
        call bstrip(stitle(is),nt)
        write (lun,221) is,
     &            ssurf(is),cavesurf(is),
     &            ln  / (0.5*rho*ssurf(is)),
     &            dn  / (0.5*rho*ssurf(is)),
     &            dnv / (0.5*rho*ssurf(is)),
     &            mle / (0.5*rho*ssurf(is)*cavesurf(is)),
     &            stitle(is)(1:nt)
      end do
      write(lun,200)

      return
      end ! outsurf


      subroutine outbody(lun)
c---------------------------------------------------------------
c     Writes formatted body forces and moments
c     for current solution to logical unit lun
c---------------------------------------------------------------
      include 'jvl.inc'
      character*50 satype
      real r(3)
      real fb(3), mb(3)
      real lb
      real drm(3)

      integer icrs(3), jcrs(3)
      data icrs / 2, 3, 1 / , jcrs / 3, 1, 2 /

      call getsa(lnasa_sa,satype,dir)

      sa = sin(alfa)
      ca = cos(alfa)

      qs  = 0.5*sref
      qsb = 0.5*sref*bref
      qsc = 0.5*sref*cref

      drm(1) = xyzref(1) - xyzmom(1) 
      drm(2) = xyzref(2) - xyzmom(2) 
      drm(3) = xyzref(3) - xyzmom(3) 

c========================================================================
c---- force components from each body
      write(lun,200)
 200  format(1x,
     &'---------------------------------------------------------------')
      write (lun,205) satype,
     &                sref,cref,bref, 
     &                xyzmom(1), xyzmom(2), xyzmom(3)
 205  format ( ' Surface Forces (referred to Sref,Cref,Bref',
     &                                   ' about Xmom,Ymom,Zmom)',
     &        /1x,a //
     &      5x,'Sref =',g11.5,3x,'Cref =',g11.5,3x,'Bref =',g11.5/
     &      5x,'Xmom =',g11.5,3x,'Ymom =',g11.5,3x,'Zmom =',g11.5 )

      write (lun,208)
 208  format(/' n',4x,'Area',9x,'CL',8x,'CD',8x,'Cm',
     &                       8x,'CY',8x,'Cn',8x,'Cl',8x,'CDi',7x,'CDv')

c---- body forces
      do ib = 1, nbody
        do k = 1, 3
          fb(k) = fbi(k,ib) + fbv(k,ib)
          mb(k) = mbi(k,ib) + mbv(k,ib)
        enddo
        mb(1) = mb(1) + drm(2)*fb(3) - drm(3)*fb(2)
        mb(2) = mb(2) + drm(3)*fb(1) - drm(1)*fb(3)
        mb(3) = mb(3) + drm(1)*fb(2) - drm(2)*fb(1)

        db = ca*fb(1) + sa*fb(3)
        yb =    fb(2)
        lb = ca*fb(3) - sa*fb(1)

        dbi = ca*fbi(1,ib) + sa*fbi(3,ib)
        dbv = ca*fbv(1,ib) + sa*fbv(3,ib)

        call bstrip(btitle(ib),nt)
        write (lun,211) ib,abody(ib),
     &      lb    / qs,
     &      db    / qs,
     &      mb(2) / qsc,
     &      yb    / qs,
     &  dir*mb(3) / qsb,
     &  dir*mb(1) / qsb,
     &      dbi   / qs,
     &      dbv   / qs,
     &      btitle(ib)(1:nt)
      end do

cc      write(lun,212)
  211 format (i2,1x,g12.5,8f10.6,3x,a)
  212 format (/)

c========================================================================
c---- body forces normalized by local reference quantities
      write (lun,220) 
  220 format (/' Body Forces (referred to Abody, Cbody ',
     &         'about nose point)'//
     &        2x,' n',3x,'Abody',7x,'Cbody',7x,'Asurf',7x,'Vbody',
     &     8x,'cl_body',4x,'cd_body',4x,'cdv_body',
     &     2x,'cmLE_nose',2x,'cnLE_nose')
 221  format (2x,i2,4g12.5,5f11.6,3x,a)

      do ib = 1, nbody
c------ set total surface force and moment
        do k = 1, 3
          fb(k) = fbi(k,ib) + fbv(k,ib)
          mb(k) = mbi(k,ib) + mbv(k,ib)
        enddo

c------ reference point for body
        l = lfrst(ib)
        r(1) = xyzref(1) - rl1(1,l)
        r(2) = xyzref(2) - rl1(2,l)
        r(3) = xyzref(3) - rl1(3,l)

        mb(1) = mb(1) + r(2)*fb(3) - r(3)*fb(2)
        mb(2) = mb(2) + r(3)*fb(1) - r(1)*fb(3)
        mb(3) = mb(3) + r(1)*fb(2) - r(2)*fb(1)

c------ rotate forces into stability axes
        db = ca*fb(1) + sa*fb(3)
        yb =    fb(2)
        lb = ca*fb(3) - sa*fb(1)

        dbi = ca*fbi(1,ib) + sa*fbi(3,ib)
        dbv = ca*fbv(1,ib) + sa*fbv(3,ib)

c       rho = 1.0
        rho = rhob(ib)
        call bstrip(btitle(ib),nt)
        write (lun,221) ib,
     &            abody(ib),cbody(ib),sbody(ib),vbody(ib),
     &            lb  / (0.5*rho*abody(ib)),
     &            db  / (0.5*rho*abody(ib)),
     &            dbv / (0.5*rho*abody(ib)),
     &            mb(2) / (0.5*rho*abody(ib)*cbody(ib)),
     &           -mb(3) / (0.5*rho*abody(ib)*cbody(ib)),
     &            btitle(ib)(1:nt)
      end do
      write(lun,200)

      return
      end ! outbody


      subroutine outstrp(lun)
c---------------------------------------------------------------
c     Writes formatted strip forces and moments
c     for current solution to logical unit lun
c---------------------------------------------------------------
      include 'jvl.inc'
      character*50 satype
      real fn(3), mn(3),
     &     fs(3), ms(3),
     &     r(3)
      real ln
      real drm(3)

      integer icrs(3), jcrs(3)
      data icrs / 2, 3, 1 / , jcrs / 3, 1, 2 /

      call getsa(lnasa_sa,satype,dir)

      if (lun.eq.0) return

      ca = cos(alfa)
      sa = sin(alfa)

      qs  = 0.5*sref
      qsb = 0.5*sref*bref
      qsc = 0.5*sref*cref

      drm(1) = xyzref(1) - xyzmom(1) 
      drm(2) = xyzref(2) - xyzmom(2) 
      drm(3) = xyzref(3) - xyzmom(3) 

c...print out the results -> forces by surface and strip
      write(lun,200)
      write(lun,210) 
      write(lun,211) satype
  200 format(1x,
     &'---------------------------------------------------------------')
  210 format (' Surface and Strip forces by surface')
  211 format (/'  Forces referred to Sref, Cref, Bref ',
     &        'about Xmom, Ymom, Zmom'/
     &         ' ',a)

      do is = 1, nsurf
        do k = 1, 3
          fn(k) = fni(k,is) + fnj(k,is) + fnv(k,is)
          mn(k) = mni(k,is) + mnj(k,is) + mnv(k,is)
        enddo
        mn(1) = mn(1) + drm(2)*fn(3) - drm(3)*fn(2)
        mn(2) = mn(2) + drm(3)*fn(1) - drm(1)*fn(3)
        mn(3) = mn(3) + drm(1)*fn(2) - drm(2)*fn(1)

        dn = ca*fn(1) + sa*fn(3)
        yn =    fn(2)
        ln = ca*fn(3) - sa*fn(1)

        dni = ca*fni(1,is) + sa*fni(3,is)
        dnv = ca*fnv(1,is) + sa*fnv(3,is)

        ns = nj(is)
        nv = nk(is)
        j1 = jfrst(is)

        write (lun,212) is,stitle(is),nv,ns,j1,ssurf(is),cavesurf(is)
        write (lun,213)
     &    ln/qs, dir*mn(1)/qsb,
     &    yn/qs,     mn(2)/qsc,
     &    dn/qs, dir*mn(3)/qsb,
     &    dni/qs,
     &    dnv/qs
 212    format (
     &       /2x,'Surface #',i2,5x,a/
     &        5x,'# chordwise =',i3,3x,'# spanwise =',i3,
     &        5x,'First strip =',i3/
     &        5x,'Surface area =',g12.5,5x,'  ave. chord =',g12.5)
 213    format (
     &         5x,'CLsurf  =',f11.6,5x,'Clsurf  =',f11.6,
     &        /5x,'CYsurf  =',f11.6,5x,'Cmsurf  =',f11.6,
     &        /5x,'CDsurf  =',f11.6,5x,'Cnsurf  =',f11.6, 
     &        /5x,'CDisurf =',f11.6,5x,'CDvsurf =',f11.6)

        rho = rhon(is)
        write (lun,214) ln/(0.5*rho*ssurf(is)), dn/(0.5*rho*ssurf(is))
        write (lun,216) 
  214   format (/'  Forces referred to Ssurf, Cave ',
     &         'about hinge axis thru LE'/
     &         5x,'cl_surf =',f11.6,5x,'cd_surf =',f11.6 )
  216   format (/' Strip Forces referred to strip area, chord'/
     &          2x,'  j ',3x,'xLE',9x,'yLE',9x,'zLE',
     &          9x,'Chord',7x,'Area',
     &          8x,'c cl',9x,'ai',10x,'cl',10x,'clj',10x,'cd',
     &          10x,'cdv',7x,'cm_c/4',5x,'c.p.x/c',8x,'Dcj')
  217   format (2x,i4,16g12.4)
        
        do j = jfrst(is), jlast(is)
          call fstrip(j,
     &                fs,ms, 
     &                caxial,cnorml, 
     &                cl_strp,cd_strp,
     &                clj_strp,cdj_strp,
     &                clt_strp,cla_strp,
     &                cmc4_strp,cmle_strp,
     &                ccl_strp )
          rho = rhos(j)
          astrp = wstrip(j)*chord(j)
          dcj_strp = djs(j)/(0.5*rho*astrp)
          xcp = 999.
          if(cl_strp.ne.0.0) then
            xcp = 0.25 - cmc4_strp/cl_strp
            xcp = min( 99.0, max( -99.0 , xcp ) )
          else
            xcp = 0.
          endif
          write (lun,217)
     &         j,(rle(l,j),l=1,3),
     &         chord(j),astrp,
     &         ccl_strp,dwwake(j),cla_strp,clj_strp,cd_strp,
     &         cdv_strp(j),cmc4_strp,xcp,dcj_strp
        end do
      end do
      write(lun,200)

      return
      end ! outstrp


      subroutine outstrpb(lun)
c---------------------------------------------------------------
c     Writes formatted strip forces and moments in body axes
c     for current solution to logical unit lun
c---------------------------------------------------------------
      include 'jvl.inc'
      character*50 satype
      real fn(3), mn(3),
     &     fs(3), ms(3),
     &     r(3)
      real ln
      real drm(3)

      integer icrs(3), jcrs(3)
      data icrs / 2, 3, 1 / , jcrs / 3, 1, 2 /

      call getsa(lnasa_sa,satype,dir)

      if (lun.eq.0) return

c...Body axes forces and moments

      qs  = 0.5*sref
      qsb = 0.5*sref*bref
      qsc = 0.5*sref*cref

      drm(1) = xyzref(1) - xyzmom(1) 
      drm(2) = xyzref(2) - xyzmom(2) 
      drm(3) = xyzref(3) - xyzmom(3) 

c...print out the results -> forces by surface and strip
      write(lun,200)
      write(lun,210) 
      write(lun,211) satype
  200 format(1x,
     &'---------------------------------------------------------------')
  210 format (' Surface and Strip forces in Body axes by surface')
  211 format (/'  Forces referred to Sref, Cref, Bref ',
     &        'about Xmom, Ymom, Zmom'/
     &         ' ',a)

      do is = 1, nsurf
        do k = 1, 3
          fn(k) = fni(k,is) + fnj(k,is) + fnv(k,is)
          mn(k) = mni(k,is) + mnj(k,is) + mnv(k,is)
        enddo
        mn(1) = mn(1) + drm(2)*fn(3) - drm(3)*fn(2)
        mn(2) = mn(2) + drm(3)*fn(1) - drm(1)*fn(3)
        mn(3) = mn(3) + drm(1)*fn(2) - drm(2)*fn(1)

        ns = nj(is)
        nv = nk(is)
        j1 = jfrst(is)

        write (lun,212) is,stitle(is),nv,ns,j1,ssurf(is),cavesurf(is)
        write (lun,213)
     &    dir*fn(1)/qs, dir*mn(1)/qsb,
     &        fn(2)/qs,     mn(2)/qsc,
     &    dir*fn(3)/qs, dir*mn(3)/qsb

 212  format (/2x,'Surface #',i2,5x,a/
     &        5x,'# chordwise =',i3,3x,'# spanwise =',i3,
     &        5x,'First strip =',i3/
     &        5x,'Surface area =',g12.5,5x,'  ave. chord =',g12.5)
  212 format (/' Body Axes Forces referred to Sref, Cref, Bref ',
     &        'about Xref, Yref, Zref'/' ',A)
  213 format ( 5X,'CFXsurf  =',G12.5,5X,'CMXsurf  =',G12.5,
     &        /5X,'CFYsurf  =',G12.5,5X,'CMYsurf  =',G12.5,
     &        /5X,'CFZsurf  =',G12.5,5X,'CMZsurf  =',G12.5)

        write (lun,216) 
  214   format (/'  Forces referred to Srf, Cave ',
     &         'about hinge axis thru LE'/
     &         5x,'cl_surf =',f11.6,5x,'cd_surf =',f11.6 )
  216   format (/' Body axes strip Forces referred to strip area, chord'
     &         ' about strip c/4'/
     &          2x,'  j ',5X,'Xc/4',7X,'Yc/4',7X,'Zc/4',
     &          9x,'Chord',7x,'Area',7x,'Dihedral',
     &          6x,'CFXstrp',6x,'CFYstrp',6x,'CFZstrp',
     &          6x,'CMXstrp',6x,'CMYstrp',6x,'CMZstrp')
  217   format (2x,i4,12g12.4)
        
        do j = jfrst(is), jlast(is)
          call fstrip(j,
     &                fs,ms, 
     &                caxial,cnorml, 
     &                cl_strp,cd_strp,
     &                clj_strp,cdj_strp,
     &                clt_strp,cla_strp,
     &                cmc4_strp,cmle_strp,
     &                ccl_strp )
          rho = rhos(j)

          xco4 = rle(1,j)+0.25*chord(j)
          dihed = -atan2(ensy(j),ensz(j))/dtr
          
          astrp = wstrip(j)*chord(j)
          dcj_strp = djs(j)/(0.5*rho*astrp)
          write (lun,217)
     &         j,(rle(l,j),l=1,3),
     &         chord(j),astrp,dihed,
     &         dir*fs(1),fs(2),dir*fs(3),
     &         dir*ms(1),ms(2),dir*ms(3))
        end do
      end do
      write(lun,200)

      return
      end ! outstrpb


      subroutine outelem(lun)
c---------------------------------------------------------------
c     Writes formatted vortex element forces and moments
c     for current solution to logical unit lun
c---------------------------------------------------------------
      include 'jvl.inc'
      character*50 satype
      real fn(3), mn(3),
     &     fs(3), ms(3),
     &     r(3)
      real ln
      real drm(3)

      integer icrs(3), jcrs(3)
      data icrs / 2, 3, 1 / , jcrs / 3, 1, 2 /

      call getsa(lnasa_sa,satype,dir)

      ca = cos(alfa)
      sa = sin(alfa)

      qs  = 0.5*sref
      qsb = 0.5*sref*bref
      qsc = 0.5*sref*cref

      drm(1) = xyzref(1) - xyzmom(1) 
      drm(2) = xyzref(2) - xyzmom(2) 
      drm(3) = xyzref(3) - xyzmom(3) 

c...forces on each strip and element (long output, and slow to printout)
      write(lun,200)
      write(lun,202)
      write(lun,205) satype

  200 format(1x,
     &'---------------------------------------------------------------')
  202 format (' Vortex strengths (by surface, by strip)')
  205 format (/'  Forces referred to Sref, Cref, Bref ',
     &        'about Xmom, Ymom, Zmom'
     &        /'  ',a)

      do is = 1, nsurf
        do k = 1, 3
          fn(k) = fni(k,is) + fnj(k,is) + fnv(k,is)
          mn(k) = mni(k,is) + mnj(k,is) + mnv(k,is)
        enddo
        mn(1) = mn(1) + drm(2)*fn(3) - drm(3)*fn(2)
        mn(2) = mn(2) + drm(3)*fn(1) - drm(1)*fn(3)
        mn(3) = mn(3) + drm(1)*fn(2) - drm(2)*fn(1)

        dn = ca*fn(1) + sa*fn(3)
        yn =    fn(2)
        ln = ca*fn(3) - sa*fn(1)

        dni = ca*fni(1,is) + sa*fni(3,is)
        dnv = ca*fnv(1,is) + sa*fnv(3,is)

c
        ns = nj(is)
        nv = nk(is)
        j1 = jfrst(is)

        write(lun,212) is,stitle(is),nv,ns,j1,
     &                  ssurf(is),cavesurf(is)
  212   format(/1x,78('*')/2x,'Surface #',i2,5x,a
     &        /5x,'# chordwise  =',i3,3x,'# spanwise =',i3,
     &         3x,'First strip  =',i3
     &        /5x,'Surface area =',g12.6,5x,'  ave. chord =',g12.6)

        write(lun,213)
     &    ln/qs, dir*mn(1)/qsb,
     &    yn/qs,     mn(2)/qsc,
     &    dn/qs, dir*mn(3)/qsb,
     &    dni/qs,
     &    dnv/qs
  213   format( 5x,'CLsurf  =',g11.5,4x,'Clsurf  =',g11.5,
     &         /5x,'CYsurf  =',g11.5,4x,'Cmsurf  =',g11.5,
     &         /5x,'CDsurf  =',g11.5,4x,'Cnsurf  =',g11.5, 
     &         /5x,'CDisurf =',g11.5,4x,'CDvsurf =',g11.5)

        rho = rhon(is)
        write(lun,214) ln/(0.5*rho*ssurf(is)), dn/(0.5*rho*ssurf(is))
  214   format(/'  Forces referred to Ssurf, Cave ',
     &          'about hinge axis thru LE'
     &        /5x,'cl_surf =',g11.5,4x,'cd_surf =',g11.5
     &        /1x,78('*'))

  232   format(/1x,'Strip #',i3,5x,'# chordwise =',i3,
     &          3x,'First vortex =',i5
     &         /4x,'xLE =',g12.4,4x,'ave. chord   =',g12.4,
     &          3x,'Incidence  =',g12.4,' deg'
     &         /4x,'yLE =',g12.4,4x,'strip width  =',g12.4,
     &          3x,'strip area =',f12.6/
     &          4x,'zLE =',g12.4,4x,'strip dihed. =',g12.4)
  233   format(/4x,'cl  =',g12.4,3x,'   cd  =',g12.4,3x,'  cdv =',g12.4,
     &         /4x,'cn  =',g12.4,3x,'   ca  =',g12.4,
     &         4x,'c cl =',g12.4,3x,'wake dnwsh =',g12.4,
     &         /4x,'cmLE=',g12.4,3x,'cm_c/4 =',g12.4,
     &        //4x,'i',
     &          4x,'x   ',8x,'y   ',8x,'z   ',8x,'dx  ',
     &          8x,'nx  ',8x,'A   ',8x,'dCp ',
     &          8x,'u   ',8x,'v   ',8x,'w   ')
  234   format(1x,i4,13g12.4)

        do jj = 1, ns
          j = j1 + jj-1
          i1 = ifrsts(j)
          rho = rhos(j)
          astrp = wstrip(j)*chord(j)
          dihed = -atan2(ensy(j),ensz(j))/dtr
          write (lun,232) j,nv,i1,
     &                    rle(1,j),chord(j),ainc(j)/dtr,
     &                    rle(2,j),wstrip(j),astrp,
     &                    rle(3,j),dihed
          call fstrip(j,
     &                fs,ms, 
     &                caxial,cnorml, 
     &                cl_strp,cd_strp,
     &                clj_strp,cdj_strp,
     &                clt_strp,cla_strp,
     &                cmc4_strp,cmle_strp,
     &                ccl_strp )
          write (lun,233) cl_strp, cd_strp, cdv_strp(j),
     &                    cnorml , caxial, 
     &                    ccli(j),  dwwake(j),
     &                    cmle_strp, cmc4_strp
          do ii = 1, nv
            i = i1 + (ii-1)
            xm = 0.5*(rv1(1,i)+rv2(1,i))
            ym = 0.5*(rv1(2,i)+rv2(2,i))
            zm = 0.5*(rv1(3,i)+rv2(3,i))
            slopen = -enc(1,i)/enc(3,i)
            anv = dxv(i)*wstrip(j)
            write (lun,234) i,xm,ym,zm,dxv(i),enc(1,i),anv,dcpv(i),
     &            vv(1,i),vv(2,i),vv(3,i)
c            write (lun,234) i,xm,ym,zm,dxv(i),slopec(i),slopen,dcp(i)
          end do
        end do
      end do
      write(lun,200)

      return
      end ! outelem

      subroutine outblin(lun)
c---------------------------------------------------------------
c     Writes formatted body line element forces and moments
c     for current solution to logical unit lun
c---------------------------------------------------------------
      include 'jvl.inc'
      character*50 satype
      real fb(3), mb(3),
     &     r(3)
      real lb
      real drm(3)

      integer icrs(3), jcrs(3)
      data icrs / 2, 3, 1 / , jcrs / 3, 1, 2 /

      call getsa(lnasa_sa,satype,dir)

      ca = cos(alfa)
      sa = sin(alfa)

      qs  = 0.5*sref
      qsb = 0.5*sref*bref
      qsc = 0.5*sref*cref

      drm(1) = xyzref(1) - xyzmom(1) 
      drm(2) = xyzref(2) - xyzmom(2) 
      drm(3) = xyzref(3) - xyzmom(3) 

c...forces on each strip and element (long output, and slow to printout)
      write(lun,200)
      write(lun,202)
      write(lun,205) satype

  200 format(1x,
     &'---------------------------------------------------------------')
  202 format (' Body line element strengths')
  205 format (/'  Forces referred to Sref, Cref, Bref ',
     &        'about Xmom, Ymom, Zmom'/
     &         '  ',a)

      do ib = 1, nbody
        do k = 1, 3
          fb(k) = fbi(k,ib) + fbv(k,ib)
          mb(k) = mbi(k,ib) + mbv(k,ib)
        enddo
        mb(1) = mb(1) + drm(2)*fb(3) - drm(3)*fb(2)
        mb(2) = mb(2) + drm(3)*fb(1) - drm(1)*fb(3)
        mb(3) = mb(3) + drm(1)*fb(2) - drm(2)*fb(1)

        db = ca*fb(1) + sa*fb(3)
        yb =    fb(2)
        lb = ca*fb(3) - sa*fb(1)

        dbi = ca*fbi(1,ib) + sa*fbi(3,ib)
        dbv = ca*fbv(1,ib) + sa*fbv(3,ib)

        write (lun,212) ib,btitle(ib), abody(ib),cbody(ib),sbody(ib)
  212   format (/1x,78('*')/2x,'Body #',i2,5x,a/
     &       5x,'Planform area =',g12.6,5x,'  Length =',g12.6,
     &       5x,'  Surface area =',g12.6)

        write (lun,213)
     &    lb/qs, dir*mb(1)/qsb,
     &    yb/qs,     mb(2)/qsc,
     &    db/qs, dir*mb(3)/qsb,
     &    dbi/qs,
     &    dbv/qs
  213   format ( 5x,'CLbody  =',g12.4,4x,'Clbody  =',g12.4,
     &          /5x,'CYbody  =',g12.4,4x,'Cmbody  =',g12.4,
     &          /5x,'CDbody  =',g12.4,4x,'Cnbody  =',g12.4, 
     &          /5x,'CDibody =',g12.4,4x,'CDvbody =',g12.4)

        write (lun,220)
 220    format(/1x,'   l ',
     &      '   x           y           z     ',
     &   '     DCpx        DCpy        DCpz   ',
     &   '     nux         nuy         nuz    ',
     &   '    sigma ')
c         12345678901 12345678901 12345678901 
        do l = lfrst(ib), llast(ib)
          write (lun,234) l,(rl(k,l),k=1,3),
     &                    (dcpb(k,l),k=1,3),
     &                     (dbl(k,l),k=1,3), src(l)
  234 format (1x,i4,12g12.4)
        enddo
      enddo ! next body ib
      write(lun,200)

      return
      end ! outblin


      subroutine outhinge(lun)
c---------------------------------------------------------------
c     Writes formatted control hinge forces and moments
c     for current solution to logical unit lun
c---------------------------------------------------------------
      include 'jvl.inc'

      write(lun,200)
 200  format(1x,
     &'---------------------------------------------------------------')

      write (lun,210) sref,cref
  210 format (
     & ' Control hinge moments' /
     & ' (referred to    Sref =',g12.5, 3x,'Cref =',g12.5,')' )

      write (lun,212) 
 212  format(/' Control          Chinge'
     &       /' ---------------- -----------')

      do n = 1, ncontrol
        write (lun,220) dname(n),chinge(n)
  220   format (1x,a16,g12.5)
      end do

      write(lun,200)

      return
      end ! outhinge


      subroutine outcnc
c---------------------------------------------------------------
c     Writes formatted strip cl*c forces
c     for current solution to logical unit lun
c---------------------------------------------------------------
      include 'jvl.inc'
      character*1 ans
      character*80 fnam
      real fs(3), ms(3)

      save fnam
      data fnam /' '/

 1000 format (a)

      call bstrip(fnam,nfn)
   10 write(*,2080) fnam(1:nfn)
 2080 format('Enter forces output file: ', a)
      read (*,1000) fnam

      lu = 18
      if(fnam.ne.' ') then
        open(lu,file=fnam,status='unknown',err=10)
        rewind(lu)
       else
c-------- just a <return> was entered...
        return
      endif

      write(*,2090) 
 2090 format('output file simple(y,ccl,cl) or full(x,y,z,ccl,cl,c)? ',$)
      read (*,1000) ans
      call touper(ans)

c...print out the results -> strip loadings
c     write (lu,210) 
      do j=1, nstrip
        call fstrip(j,
     &                fs,ms, 
     &                caxial,cnorml, 
     &                cl_strp,cd_strp,
     &                clj_strp,cdj_strp,
     &                clt_strp,cla_strp,
     &                cmc4_strp,cmle_strp,
     &                ccl_strp )
        i = ifrsts(j)
        xm = 0.5*(rv1(1,i)+rv2(1,i))
        ym = 0.5*(rv1(2,i)+rv2(2,i))
        zm = 0.5*(rv1(3,i)+rv2(3,i))
        cclm = ccl_strp
        clm  = cl_strp
        chm  = chord(j)
        dym  = wstrip(j)
        asm  = dym*chm
        if(ans.eq.'s') then
          write (lu,213) ym,cclm,clm,chm
         else
          write (lu,212) xm,ym,zm,cclm,clm,chm,dym,asm
        endif
      end do
      close(lu)

  210 format (//' *** strip loadings')
  212 format (1x,3(g10.4,1x),2(g10.4,1x),2(g10.4,1x),g10.4)
  213 format (1x,g10.4,1x,3(g10.4,1x))

      return
      end ! outcnc


      subroutine outdcp(lun,jstrip)
c---------------------------------------------------------------
c     Writes formatted strip dCp loadings 
c     for current solution to logical unit lun
c---------------------------------------------------------------
      include 'jvl.inc'

      qs  = 0.5*sref
      qsb = 0.5*sref*bref
      qsc = 0.5*sref*cref

      j = jstrip

      xle = rle(1,j)
      xte = rle(1,j) + chord(j)

      write(lun,231) alfa/dtr, dqt/sref, djt/qs, 0.5*det/qs,
     &               ltot/qs,mt(2)/qsc

      write(lun,232) j, xle, chord(j)

      astrp = wstrip(j)*chord(j)
      dihed = -atan2(ensy(j),ensz(j))/dtr
      do i = ifrsts(j), ilasts(j)
        xm = 0.5*(rv1(1,i)+rv2(1,i))
        ym = 0.5*(rv1(2,i)+rv2(2,i))
        zm = 0.5*(rv1(3,i)+rv2(3,i))
        slopen = -enc(1,i)/enc(3,i)
        xoc = (xm - xle) / chord(j)
        write (lun,234) xoc, dcpv(i)
      end do

      return

 231  format('# a =',F7.3,
     &        3x,'DCQ =',g10.4,
     &        3x,'DCJ =',g10.4,
     &        3x,'DCE =',g10.4,
     &        3x,'CL =',g10.4,
     &        3x,'CM =',g10.4)
 232  format ('#  Strip',i3,
     &        4x,'XLE =',g10.4,
     &        4x,'chord =',g10.4)
 234  format (2x,5g10.4)

      end ! outdcp


      subroutine dermatm(lu)
c------------------------------------------------------------
c     Calculates stability derivative matrix for current
c     alfa,beta,p,q,r, and outputs is to logical unit lu
c------------------------------------------------------------
      include 'jvl.inc'
      character*50 satype

      real clt_u(numax),
     &     cyt_u(numax),
     &     cdt_u(numax)
      real clt_d(ndmax),
     &     cyt_d(ndmax),
     &     cdt_d(ndmax)
      real clt_j(njmax),
     &     cyt_j(njmax),
     &     cdt_j(njmax)
      real clt_g(ngmax),
     &     cyt_g(ngmax),
     &     cdt_g(ngmax)

      real drm(3)
      real cmt_u(3,numax)
      real cmt_d(3,ndmax)
      real cmt_j(3,njmax)
      real cmt_g(3,ngmax)

      real cdff_d(ndmax)
      real cdff_j(njmax)
      real cdff_g(ngmax)

      real wbar_a(3), wbar_rx(3), wbar_rz(3)

      call getsa(lnasa_sa,satype,dir)

c---- set freestream velocity components from alpha, beta
      call ab2vbar

c---- calculate forces and sensitivities
      call aero

      qs  = 0.5*sref
      qsb = 0.5*sref*bref
      qsc = 0.5*sref*cref

      drm(1) = xyzref(1) - xyzmom(1)
      drm(2) = xyzref(2) - xyzmom(2)
      drm(3) = xyzref(3) - xyzmom(3)

      do n = 1, numax
        clt_u(n) = ltot_u(n)/qs
        cyt_u(n) = ytot_u(n)/qs
        cdt_u(n) = dtot_u(n)/qs
        cmt_u(1,n) = ( mt_u(1,n) + drm(2)*ft_u(3,n)
     &                           - drm(3)*ft_u(2,n) )/qsb
        cmt_u(2,n) = ( mt_u(2,n) + drm(3)*ft_u(1,n)
     &                           - drm(1)*ft_u(3,n) )/qsc
        cmt_u(3,n) = ( mt_u(3,n) + drm(1)*ft_u(2,n)
     &                           - drm(2)*ft_u(1,n) )/qsb
      enddo

      do n = 1, ncontrol
        clt_d(n) = ltot_d(n)/qs
        cyt_d(n) = ytot_d(n)/qs
        cdt_d(n) = dtot_d(n)/qs
        cmt_d(1,n) = ( mt_d(1,n) + drm(2)*ft_d(3,n)
     &                           - drm(3)*ft_d(2,n) )/qsb
        cmt_d(2,n) = ( mt_d(2,n) + drm(3)*ft_d(1,n)
     &                           - drm(1)*ft_d(3,n) )/qsc
        cmt_d(3,n) = ( mt_d(3,n) + drm(1)*ft_d(2,n)
     &                           - drm(2)*ft_d(1,n) )/qsb
        cdff_d(n) = dff_d(n)/qs
      enddo

      do n = 1, nvarjet
        clt_j(n) = ltot_j(n)/qs
        cyt_j(n) = ytot_j(n)/qs
        cdt_j(n) = dtot_j(n)/qs
        cmt_j(1,n) = ( mt_j(1,n) + drm(2)*ft_j(3,n)
     &                           - drm(3)*ft_j(2,n) )/qsb
        cmt_j(2,n) = ( mt_j(2,n) + drm(3)*ft_j(1,n)
     &                           - drm(1)*ft_j(3,n) )/qsc
        cmt_j(3,n) = ( mt_j(3,n) + drm(1)*ft_j(2,n)
     &                           - drm(2)*ft_j(1,n) )/qsb
        cdff_j(n) = dff_j(n)/qs
      enddo

      do n = 1, ndesign
        clt_g(n) = ltot_g(n)/qs
        cyt_g(n) = ytot_g(n)/qs
        cdt_g(n) = dtot_g(n)/qs
        cmt_g(1,n) = ( mt_g(1,n) + drm(2)*ft_g(3,n)
     &                           - drm(3)*ft_g(2,n) )/qsb
        cmt_g(2,n) = ( mt_g(2,n) + drm(3)*ft_g(1,n)
     &                           - drm(1)*ft_g(3,n) )/qsc
        cmt_g(3,n) = ( mt_g(3,n) + drm(1)*ft_g(2,n)
     &                           - drm(2)*ft_g(1,n) )/qsb
        cdff_g(n) = dff_g(n)/qs
      enddo

c---- set stability-axes rates (rx,ry,rz) in terms of body-axes rates
      ca = cos(alfa)
      sa = sin(alfa)
      rx = (wbar(1)*ca + wbar(3)*sa) * dir
      ry =  wbar(2)
      rz = (wbar(3)*ca - wbar(1)*sa) * dir

c---- now vice-versa, and set sensitivities (which is what's really needed)
cc    wbar(1)    =  rx*ca - rz*sa
cc    wbar(2)    =  ry
cc    wbar(3)    =  rz*ca + rx*sa

      wbar_rx(1) = ca     * dir
      wbar_rx(2) = 0.
      wbar_rx(3) =     sa * dir

      wbar_rz(1) =    -sa * dir
      wbar_rz(2) = 0.
      wbar_rz(3) = ca     * dir

      wbar_a(1)  = -rx*sa - rz*ca   !!! = -wbar(3)
      wbar_a(2)  =  0.
      wbar_a(3)  = -rz*sa + rx*ca   !!! =  wbar(1)

c---- set force derivatives in stability axes
      cl_al = clt_u(1)*vbar_a(1) + clt_u(4)*wbar_a(1)
     &      + clt_u(2)*vbar_a(2) + clt_u(5)*wbar_a(2)
     &      + clt_u(3)*vbar_a(3) + clt_u(6)*wbar_a(3)
      cl_be = clt_u(1)*vbar_b(1)
     &      + clt_u(2)*vbar_b(2)
     &      + clt_u(3)*vbar_b(3)
      cl_rx = clt_u(4)*wbar_rx(1) + clt_u(6)*wbar_rx(3)
      cl_ry = clt_u(5)
      cl_rz = clt_u(6)*wbar_rz(3) + clt_u(4)*wbar_rz(1)

      cy_al = cyt_u(1)*vbar_a(1) + cyt_u(4)*wbar_a(1)
     &      + cyt_u(2)*vbar_a(2) + cyt_u(5)*wbar_a(2)
     &      + cyt_u(3)*vbar_a(3) + cyt_u(6)*wbar_a(3)
      cy_be = cyt_u(1)*vbar_b(1)
     &      + cyt_u(2)*vbar_b(2)
     &      + cyt_u(3)*vbar_b(3)
      cy_rx = cyt_u(4)*wbar_rx(1) + cyt_u(6)*wbar_rx(3)
      cy_ry = cyt_u(5)
      cy_rz = cyt_u(6)*wbar_rz(3) + cyt_u(4)*wbar_rz(1)

      cd_al = cdt_u(1)*vbar_a(1) + cdt_u(4)*wbar_a(1)
     &      + cdt_u(2)*vbar_a(2) + cdt_u(5)*wbar_a(2)
     &      + cdt_u(3)*vbar_a(3) + cdt_u(6)*wbar_a(3)
      cd_be = cdt_u(1)*vbar_b(1)
     &      + cdt_u(2)*vbar_b(2)
     &      + cdt_u(3)*vbar_b(3)
      cd_rx = cdt_u(4)*wbar_rx(1) + cdt_u(6)*wbar_rx(3)
      cd_ry = cdt_u(5)
      cd_rz = cdt_u(6)*wbar_rz(3) + cdt_u(4)*wbar_rz(1)

      cr_al = cmt_u(1,1)*vbar_a(1) + cmt_u(1,4)*wbar_a(1)
     &      + cmt_u(1,2)*vbar_a(2) + cmt_u(1,5)*wbar_a(2)
     &      + cmt_u(1,3)*vbar_a(3) + cmt_u(1,6)*wbar_a(3)
      cr_be = cmt_u(1,1)*vbar_b(1)
     &      + cmt_u(1,2)*vbar_b(2)
     &      + cmt_u(1,3)*vbar_b(3)
      cr_rx = cmt_u(1,4)*wbar_rx(1) + cmt_u(1,6)*wbar_rx(3)
      cr_ry = cmt_u(1,5)
      cr_rz = cmt_u(1,6)*wbar_rz(3) + cmt_u(1,4)*wbar_rz(1)

      cm_al = cmt_u(2,1)*vbar_a(1) + cmt_u(2,4)*wbar_a(1)
     &      + cmt_u(2,2)*vbar_a(2) + cmt_u(2,5)*wbar_a(2)
     &      + cmt_u(2,3)*vbar_a(3) + cmt_u(2,6)*wbar_a(3)
      cm_be = cmt_u(2,1)*vbar_b(1)
     &      + cmt_u(2,2)*vbar_b(2)
     &      + cmt_u(2,3)*vbar_b(3)
      cm_rx = cmt_u(2,4)*wbar_rx(1) + cmt_u(2,6)*wbar_rx(3)
      cm_ry = cmt_u(2,5)
      cm_rz = cmt_u(2,6)*wbar_rz(3) + cmt_u(2,4)*wbar_rz(1)

      cn_al = cmt_u(3,1)*vbar_a(1) + cmt_u(3,4)*wbar_a(1)
     &      + cmt_u(3,2)*vbar_a(2) + cmt_u(3,5)*wbar_a(2)
     &      + cmt_u(3,3)*vbar_a(3) + cmt_u(3,6)*wbar_a(3)
      cn_be = cmt_u(3,1)*vbar_b(1)
     &      + cmt_u(3,2)*vbar_b(2)
     &      + cmt_u(3,3)*vbar_b(3)
      cn_rx = cmt_u(3,4)*wbar_rx(1) + cmt_u(3,6)*wbar_rx(3)
      cn_ry = cmt_u(3,5)
      cn_rz = cmt_u(3,6)*wbar_rz(3) + cmt_u(3,4)*wbar_rz(1)

      call outtot(lu)

      write(lu,7004)
 7004 format(/' derivatives...')

      write(lu,7006)
 7006 format(14x, 4x,'           alpha',
     &            4x,'            beta'
     &      /14x, 4x,'----------------',
     &            4x,'----------------')

      write(lu,7010) cl_al, cl_be
 7010 format(' z force     |','    CLa =',f11.6,'    CLb =',f11.6)

      write(lu,7020) cy_al, cy_be
 7020 format(' y force     |','    CYa =',f11.6,'    CYb =',f11.6)

      write(lu,7010) cd_al, cd_be
 7030 format(' x force     |','    CDa =',f11.6,'    CDb =',f11.6)

      write(lu,7040) dir*cr_al, dir*cr_be
 7040 format(' roll  x mom.|','    Cla =',f11.6,'    Clb =',f11.6)

      write(lu,7050) cm_al, cm_be
 7050 format(' pitch y mom.|','    Cma =',f11.6,'    Cmb =',f11.6)

      write(lu,7060) dir*cn_al, dir*cn_be
 7060 format(' yaw   z mom.|','    Cna =',f11.6,'    Cnb =',f11.6)

      write(lu,7106)
 7106 format(/14x, 4x,'    roll rate  p',
     &             4x,'   pitch rate  q',
     &             4x,'     yaw rate  r'
     &       /14x, 4x,'----------------',
     &             4x,'----------------',
     &             4x,'----------------' )

      write(lu,7110) cl_rx*2.0/bref, 
     &               cl_ry*2.0/cref, 
     &               cl_rz*2.0/bref
 7110 format(' z force     |','    CLp =',f11.6,
     &                        '    CLq =',f11.6,
     &                        '    CLr =',f11.6 )

      write(lu,7120) cy_rx*2.0/bref,
     &               cy_ry*2.0/cref,
     &               cy_rz*2.0/bref
 7120 format(' y force     |','    CYp =',f11.6,
     &                        '    CYq =',f11.6,
     &                        '    CYr =',f11.6 )

      write(lu,7140) dir*cr_rx*2.0/bref, 
     &               dir*cr_ry*2.0/cref,
     &               dir*cr_rz*2.0/bref
 7140 format(' roll  x mom.|','    Clp =',f11.6,
     &                        '    Clq =',f11.6,
     &                        '    Clr =',f11.6 )

      write(lu,7150) cm_rx*2.0/bref,
     &               cm_ry*2.0/cref,
     &               cm_rz*2.0/bref
 7150 format(' pitch y mom.|','    Cmp =',f11.6,
     &                        '    Cmq =',f11.6,
     &                        '    Cmr =',f11.6 )

      write(lu,7160) dir*cn_rx*2.0/bref,
     &               dir*cn_ry*2.0/cref,
     &               dir*cn_rz*2.0/bref
 7160 format(' yaw   z mom.|','    Cnp =',f11.6,
     &                        '    Cnq =',f11.6,
     &                        '    Cnr =',f11.6 )


      if(ncontrol.gt.0) then

      write(lu,8106) (dname(k), k, k=1, ncontrol)
 8106 format(/14x,20(4x,a12, ' d',i1,' '))
      write(lu,8107) (' ',k=1, ncontrol)
 8107 format( 14x,20(3x,a,'----------------'))

      write(lu,8110) (' ',k,clt_d(k), k=1, ncontrol)
 8110 format(' z force     |',20(a,'  CLd',i1,' =',f11.6))

      write(lu,8120) (' ',k,cyt_d(k), k=1, ncontrol)
 8120 format(' y force     |',20(a,'  CYd',i1,' =',f11.6))

      write(lu,8140) (' ',k,dir*cmt_d(1,k), k=1, ncontrol)
 8140 format(' roll  x mom.|',20(a,'  Cld',i1,' =',f11.6))

      write(lu,8150) (' ',k,    cmt_d(2,k), k=1, ncontrol)
 8150 format(' pitch y mom.|',20(a,'  Cmd',i1,' =',f11.6))

      write(lu,8160) (' ',k,dir*cmt_d(3,k), k=1, ncontrol)
 8160 format(' yaw   z mom.|',20(a,'  Cnd',i1,' =',f11.6))

      write(lu,8170) (' ',k,    cdff_d(k), k=1, ncontrol)
 8170 format(' Trefftz drag|',20(a,'CDffd',i1,' =',f11.6))

      write(lu,8180) (' ',k,  spanef_d(k), k=1, ncontrol)
 8180 format(' span eff.   |',20(a,'   ed',i1,' =',f11.6))

      write(lu,*)
      write(lu,*)

      endif

      if(nvarjet.gt.0) then

      write(lu,8206) (jname(k), k, k=1, nvarjet)
 8206 format(/14x,20(4x,a12, ' j',i1,' '))
      write(lu,8207) (' ',k=1, nvarjet)
 8207 format( 14x,20(3x,a,'----------------'))

      write(lu,8210) (' ',k,clt_j(k), k=1, nvarjet)
 8210 format(' z force     |',20(a,'  CLj',i1,' =',f11.6))

      write(lu,8220) (' ',k,cyt_j(k), k=1, nvarjet)
 8220 format(' y force     |',20(a,'  CYj',i1,' =',f11.6))

      write(lu,8230) (' ',k,dir*cmt_j(1,k), k=1, nvarjet)
 8230 format(' roll  x mom.|',20(a,'  Clj',i1,' =',f11.6))

      write(lu,8240) (' ',k,    cmt_j(2,k), k=1, nvarjet)
 8240 format(' pitch y mom.|',20(a,'  Cmj',i1,' =',f11.6))

      write(lu,8250) (' ',k,dir*cmt_j(3,k), k=1, nvarjet)
 8250 format(' yaw   z mom.|',20(a,'  Cnj',i1,' =',f11.6))

      write(lu,8260) (' ',k,    cdff_j(k), k=1, nvarjet)
 8260 format(' trefftz drag|',20(a,'CDffj',i1,' =',f11.6))

      write(lu,8270) (' ',k,  spanef_j(k), k=1, nvarjet)
 8270 format(' span eff.   |',20(a,'   ej',i1,' =',f11.6))

      write(lu,*)
      write(lu,*)

      endif

      if(ndesign.gt.0) then

      write(lu,8306) (gname(k), k, k=1, ndesign)
 8306 format(/14x,20(4x,a12, ' g',i1,' '))
      write(lu,8307) (' ',k=1, ndesign)
 8307 format( 14x,20(3x,a,'----------------'))

      write(lu,8310) (' ',k,clt_g(k), k=1, ndesign)
 8310 format(' z force     |',20(a,'  CLg',i1,' =',f11.6))

      write(lu,8320) (' ',k,cyt_g(k), k=1, ndesign)
 8320 format(' y force     |',20(a,'  CYg',i1,' =',f11.6))

      write(lu,8330) (' ',k,dir*cmt_g(1,k), k=1, ndesign)
 8330 format(' roll  x mom.|',20(a,'  Clg',i1,' =',f11.6))

      write(lu,8340) (' ',k,    cmt_g(2,k), k=1, ndesign)
 8340 format(' pitch y mom.|',20(a,'  Cmg',i1,' =',f11.6))

      write(lu,8350) (' ',k,dir*cmt_g(3,k), k=1, ndesign)
 8350 format(' yaw   z mom.|',20(a,'  Cng',i1,' =',f11.6))

      write(lu,8360) (' ',k,    cdff_g(k), k=1, ndesign)
 8360 format(' trefftz drag|',20(a,'CDffg',i1,' =',f11.6))

      write(lu,8370) (' ',k,  spanef_g(k), k=1, ndesign)
 8370 format(' span eff.   |',20(a,'   eg',i1,' =',f11.6))

      write(lu,*)
      write(lu,*)

      endif

      if(cl_al .ne. 0.0) then
       xnp = xyzmom(1) - cref*cm_al/cl_al
       write(lu,8401) xnp
 8401  format(/' Neutral point  Xnp =', g14.6)
      endif

      if(abs(cr_rz*cn_be*(2.0/bref)) .gt. 0.0001) then
       bb = cr_be*cn_rz / (cr_rz*cn_be)
       write(lu,8402) bb 
 8402  format(/' Clb Cnr / Clr Cnb  =', f11.6,
     &    '    (  > 1 if spirally stable )')
      endif


      return
      end ! dermatm



      subroutine dermats(lu)
c------------------------------------------------------------
c     Calculates stability derivative matrix for current
c     alfa, and outputs is to logical unit lu
c     Calculation uses simplified form from e.g. Etkin
c------------------------------------------------------------
      include 'jvl.inc'
      character*50 satype
      real wbar_rx(3), wbar_rz(3), wbar_a(3)
      real cmsax(3),
     &     cmsax_a(3),
     &     cmsax_u(3,numax),
     &     cmsax_d(3,ndmax),
     &     cmsax_j(3,njmax),
     &     cmsax_g(3,ngmax)

      real clt_u(numax),
     &     cyt_u(numax),
     &     cdt_u(numax)
      real clt_d(ndmax),
     &     cyt_d(ndmax),
     &     cdt_d(ndmax)
      real clt_j(njmax),
     &     cyt_j(njmax),
     &     cdt_j(njmax)
      real clt_g(ngmax),
     &     cyt_g(ngmax),
     &     cdt_g(ngmax)

      real drm(3)
      real mtm(3)
      real mtm_u(3,numax)
      real mtm_d(3,ndmax)
      real mtm_j(3,njmax)
      real mtm_g(3,ngmax)

      real cmt(3)
      real cmt_u(3,numax)
      real cmt_d(3,ndmax)
      real cmt_j(3,njmax)
      real cmt_g(3,ngmax)

      real cdff_d(ndmax)
      real cdff_j(njmax)
      real cdff_g(ngmax)

      call getsa(lnasa_sa,satype,dir)

c---- set freestream velocity components from alpha, beta
      call ab2vbar

c---- calculate forces and sensitivities
      call aero


      qs  = 0.5*sref
      qsb = 0.5*sref*bref
      qsc = 0.5*sref*cref

      drm(1) = xyzref(1) - xyzmom(1)
      drm(2) = xyzref(2) - xyzmom(2)
      drm(3) = xyzref(3) - xyzmom(3)
      mtm(1) = mt(1) + drm(2)*ft(3) - drm(3)*ft(2)
      mtm(2) = mt(2) + drm(3)*ft(1) - drm(1)*ft(3)
      mtm(3) = mt(3) + drm(1)*ft(2) - drm(2)*ft(1)
      do n = 1, numax
        mtm_u(1,n) = mt_u(1,n) + drm(2)*ft_u(3,n) - drm(3)*ft_u(2,n)
        mtm_u(2,n) = mt_u(2,n) + drm(3)*ft_u(1,n) - drm(1)*ft_u(3,n)
        mtm_u(3,n) = mt_u(3,n) + drm(1)*ft_u(2,n) - drm(2)*ft_u(1,n)
      enddo

      do n = 1, ncontrol
        mtm_d(1,n) = mt_d(1,n) + drm(2)*ft_d(3,n) - drm(3)*ft_d(2,n)
        mtm_d(2,n) = mt_d(2,n) + drm(3)*ft_d(1,n) - drm(1)*ft_d(3,n)
        mtm_d(3,n) = mt_d(3,n) + drm(1)*ft_d(2,n) - drm(2)*ft_d(1,n)
      enddo

      do n = 1, nvarjet
        mtm_j(1,n) = mt_j(1,n) + drm(2)*ft_j(3,n) - drm(3)*ft_j(2,n)
        mtm_j(2,n) = mt_j(2,n) + drm(3)*ft_j(1,n) - drm(1)*ft_j(3,n)
        mtm_j(3,n) = mt_j(3,n) + drm(1)*ft_j(2,n) - drm(2)*ft_j(1,n)
      enddo

      do n = 1, ndesign
        mtm_g(1,n) = mt_g(1,n) + drm(2)*ft_g(3,n) - drm(3)*ft_g(2,n)
        mtm_g(2,n) = mt_g(2,n) + drm(3)*ft_g(1,n) - drm(1)*ft_g(3,n)
        mtm_g(3,n) = mt_g(3,n) + drm(1)*ft_g(2,n) - drm(2)*ft_g(1,n)
      enddo


      cmt(1) = mtm(1)/qsb
      cmt(2) = mtm(2)/qsc
      cmt(3) = mtm(3)/qsb
      do n = 1, numax
        clt_u(n) = ltot_u(n)/qs
        cyt_u(n) = ytot_u(n)/qs
        cdt_u(n) = dtot_u(n)/qs
        cmt_u(1,n) = mtm_u(1,n)/qsb
        cmt_u(2,n) = mtm_u(2,n)/qsc
        cmt_u(3,n) = mtm_u(3,n)/qsb
      enddo

      do n = 1, ncontrol
        clt_d(n) = ltot_d(n)/qs
        cyt_d(n) = ytot_d(n)/qs
        cdt_d(n) = dtot_d(n)/qs
        cmt_d(1,n) = mtm_d(1,n)/qsb
        cmt_d(2,n) = mtm_d(2,n)/qsc
        cmt_d(3,n) = mtm_d(3,n)/qsb
        cdff_d(n) = dff_d(n)/qs
      enddo

      do n = 1, nvarjet
        clt_j(n) = ltot_j(n)/qs
        cyt_j(n) = ytot_j(n)/qs
        cdt_j(n) = dtot_j(n)/qs
        cmt_j(1,n) = mtm_j(1,n)/qsb
        cmt_j(2,n) = mtm_j(2,n)/qsc
        cmt_j(3,n) = mtm_j(3,n)/qsb
        cdff_j(n) = dff_j(n)/qs
      enddo

      do n = 1, ndesign
        clt_g(n) = ltot_g(n)/qs
        cyt_g(n) = ytot_g(n)/qs
        cdt_g(n) = dtot_g(n)/qs
        cmt_g(1,n) = mtm_g(1,n)/qsb
        cmt_g(2,n) = mtm_g(2,n)/qsc
        cmt_g(3,n) = mtm_g(3,n)/qsb
        cdff_g(n) = dff_g(n)/qs
      enddo

c---- set stability-axes rates (rx,ry,rz) in terms of body-axes rates
      ca = cos(alfa)
      sa = sin(alfa)

      rx = (wbar(1)*ca + wbar(3)*sa) * dir
      ry =  wbar(2)
      rz = (wbar(3)*ca - wbar(1)*sa) * dir

c---- now vice-versa, and set sensitivities (which is what's really needed)
cc    wbar(1)    = (rx*ca - rz*sa) * dir
cc    wbar(2)    =  ry
cc    wbar(3)    = (rz*ca + rx*sa) * dir

      wbar_rx(1) = ca     * dir
      wbar_rx(2) = 0.
      wbar_rx(3) =     sa * dir

      wbar_rz(1) =    -sa * dir
      wbar_rz(2) = 0.
      wbar_rz(3) = ca     * dir

      wbar_a(1)  = -rx*sa - rz*ca   !!! = -wbar(3)
      wbar_a(2)  =  0.
      wbar_a(3)  = -rz*sa + rx*ca   !!! =  wbar(1)


      cmsax(1)   =  cmt(1)*ca + cmt(3)*sa
      cmsax(2)   =  cmt(2)
      cmsax(3)   =  cmt(3)*ca - cmt(1)*sa
      cmsax_a(1) = -cmt(1)*sa + cmt(3)*ca
      cmsax_a(2) =  0.      
      cmsax_a(3) = -cmt(3)*sa - cmt(1)*ca

      do n = 1, numax
        cmsax_u(1,n) = cmt_u(1,n)*ca + cmt_u(3,n)*sa
        cmsax_u(2,n) = cmt_u(2,n)                  
        cmsax_u(3,n) = cmt_u(3,n)*ca - cmt_u(1,n)*sa
      enddo
      do n = 1, ncontrol
        cmsax_d(1,n) = cmt_d(1,n)*ca + cmt_d(3,n)*sa
        cmsax_d(2,n) = cmt_d(2,n)                  
        cmsax_d(3,n) = cmt_d(3,n)*ca - cmt_d(1,n)*sa
      enddo
      do n = 1, nvarjet
        cmsax_j(1,n) = cmt_j(1,n)*ca + cmt_j(3,n)*sa
        cmsax_j(2,n) = cmt_j(2,n)                  
        cmsax_j(3,n) = cmt_j(3,n)*ca - cmt_j(1,n)*sa
      enddo
      do n = 1, ndesign
        cmsax_g(1,n) = cmt_g(1,n)*ca + cmt_g(3,n)*sa
        cmsax_g(2,n) = cmt_g(2,n)                  
        cmsax_g(3,n) = cmt_g(3,n)*ca - cmt_g(1,n)*sa
      enddo


c---- set force derivatives in stability axes
      cl_al = clt_u(1)*vbar_a(1) + clt_u(4)*wbar_a(1)
     &      + clt_u(2)*vbar_a(2) + clt_u(5)*wbar_a(2)
     &      + clt_u(3)*vbar_a(3) + clt_u(6)*wbar_a(3)
      cl_be = clt_u(1)*vbar_b(1)
     &      + clt_u(2)*vbar_b(2)
     &      + clt_u(3)*vbar_b(3)
      cl_rx = clt_u(4)*wbar_rx(1) + clt_u(6)*wbar_rx(3)
      cl_ry = clt_u(5)
      cl_rz = clt_u(6)*wbar_rz(3) + clt_u(4)*wbar_rz(1)

      cy_al = cyt_u(1)*vbar_a(1) + cyt_u(4)*wbar_a(1)
     &      + cyt_u(2)*vbar_a(2) + cyt_u(5)*wbar_a(2)
     &      + cyt_u(3)*vbar_a(3) + cyt_u(6)*wbar_a(3)
      cy_be = cyt_u(1)*vbar_b(1)
     &      + cyt_u(2)*vbar_b(2)
     &      + cyt_u(3)*vbar_b(3)
      cy_rx = cyt_u(4)*wbar_rx(1) + cyt_u(6)*wbar_rx(3)
      cy_ry = cyt_u(5)
      cy_rz = cyt_u(6)*wbar_rz(3) + cyt_u(4)*wbar_rz(1)

      cd_al = cdt_u(1)*vbar_a(1) + cdt_u(4)*wbar_a(1)
     &      + cdt_u(2)*vbar_a(2) + cdt_u(5)*wbar_a(2)
     &      + cdt_u(3)*vbar_a(3) + cdt_u(6)*wbar_a(3)
      cd_be = cdt_u(1)*vbar_b(1)
     &      + cdt_u(2)*vbar_b(2)
     &      + cdt_u(3)*vbar_b(3)
      cd_rx = cdt_u(4)*wbar_rx(1) + cdt_u(6)*wbar_rx(3)
      cd_ry = cdt_u(5)
      cd_rz = cdt_u(6)*wbar_rz(3) + cdt_u(4)*wbar_rz(1)

      cr_al = cmsax_u(1,1)*vbar_a(1) + cmsax_u(1,4)*wbar_a(1)
     &      + cmsax_u(1,2)*vbar_a(2) + cmsax_u(1,5)*wbar_a(2)
     &      + cmsax_u(1,3)*vbar_a(3) + cmsax_u(1,6)*wbar_a(3)
     &      + cmsax_a(1)
      cr_be = cmsax_u(1,1)*vbar_b(1)
     &      + cmsax_u(1,2)*vbar_b(2)
     &      + cmsax_u(1,3)*vbar_b(3)
      cr_rx = cmsax_u(1,4)*wbar_rx(1) + cmsax_u(1,6)*wbar_rx(3)
      cr_ry = cmsax_u(1,5)
      cr_rz = cmsax_u(1,6)*wbar_rz(3) + cmsax_u(1,4)*wbar_rz(1)

      cm_al = cmsax_u(2,1)*vbar_a(1) + cmsax_u(2,4)*wbar_a(1)
     &      + cmsax_u(2,2)*vbar_a(2) + cmsax_u(2,5)*wbar_a(2)
     &      + cmsax_u(2,3)*vbar_a(3) + cmsax_u(2,6)*wbar_a(3)
     &      + cmsax_a(2)
      cm_be = cmsax_u(2,1)*vbar_b(1)
     &      + cmsax_u(2,2)*vbar_b(2)
     &      + cmsax_u(2,3)*vbar_b(3)
      cm_rx = cmsax_u(2,4)*wbar_rx(1) + cmsax_u(2,6)*wbar_rx(3)
      cm_ry = cmsax_u(2,5)
      cm_rz = cmsax_u(2,6)*wbar_rz(3) + cmsax_u(2,4)*wbar_rz(1)

      cn_al = cmsax_u(3,1)*vbar_a(1) + cmsax_u(3,4)*wbar_a(1)
     &      + cmsax_u(3,2)*vbar_a(2) + cmsax_u(3,5)*wbar_a(2)
     &      + cmsax_u(3,3)*vbar_a(3) + cmsax_u(3,6)*wbar_a(3)
     &      + cmsax_a(3)
      cn_be = cmsax_u(3,1)*vbar_b(1)
     &      + cmsax_u(3,2)*vbar_b(2)
     &      + cmsax_u(3,3)*vbar_b(3)
      cn_rx = cmsax_u(3,4)*wbar_rx(1) + cmsax_u(3,6)*wbar_rx(3)
      cn_ry = cmsax_u(3,5)
      cn_rz = cmsax_u(3,6)*wbar_rz(3) + cmsax_u(3,4)*wbar_rz(1)

      call outtot(lu)

      write(lu,7004)
 7004 format(/' Stability-axis derivatives...')

      write(lu,7006)
 7006 format(/14x, 4x,'           alpha',
     &             4x,'            beta'
     &       /14x, 4x,'----------------',
     &             4x,'----------------')

      write(lu,7010) cl_al, cl_be
 7010 format(' z'' force CL |' ,'    CLa =',f11.6,'    CLb =',f11.6)

      write(lu,7020) cy_al, cy_be
 7020 format(' y  force CY |'  ,'    CYa =',f11.6,'    CYb =',f11.6)

      write(lu,7030) cd_al, cd_be
 7030 format(' x'' force CD |' ,'    CDa =',f11.6,'    CDb =',f11.6)

      write(lu,7040) dir*cr_al, dir*cr_be
 7040 format(' x'' mom.  Cl''|','    Cla =',f11.6,'    Clb =',f11.6)

      write(lu,7050)     cm_al,     cm_be
 7050 format(' y  mom.  Cm |'  ,'    Cma =',f11.6,'    Cmb =',f11.6)

      write(lu,7060) dir*cn_al, dir*cn_be
 7060 format(' z'' mom.  Cn''|','    Cna =',f11.6,'    Cnb =',f11.6)

      write(lu,7106)
 7106 format(/14x, 4x,'   roll rate  p''',
     &             4x,'  pitch rate  q''',
     &             4x,'    yaw rate  r'''
     &       /14x, 4x,'----------------',
     &             4x,'----------------',
     &             4x,'----------------' )

      write(lu,7110) cl_rx*(2.0/bref), 
     &               cl_ry*(2.0/cref), 
     &               cl_rz*(2.0/bref)
 7110 format(' z'' force CL |','    CLp =',f11.6,
     &                         '    CLq =',f11.6,
     &                         '    CLr =',f11.6 )

      write(lu,7120) cy_rx*(2.0/bref),
     &               cy_ry*(2.0/cref),
     &               cy_rz*(2.0/bref)
 7120 format(' y  force CY |','    CYp =',f11.6,
     &                        '    CYq =',f11.6,
     &                        '    CYr =',f11.6 )

      write(lu,7130) cd_rx*(2.0/bref), 
     &               cd_ry*(2.0/cref), 
     &               cd_rz*(2.0/bref)
 7130 format(' x'' force CD |','    CDp =',f11.6,
     &                         '    CDq =',f11.6,
     &                         '    CDr =',f11.6 )


      write(lu,7140) dir*cr_rx*(2.0/bref), 
     &               dir*cr_ry*(2.0/cref), 
     &               dir*cr_rz*(2.0/bref)
 7140 format(' x'' mom.  Cl''|','    Clp =',f11.6,
     &                        '    Clq =',f11.6,
     &                        '    Clr =',f11.6 )

      write(lu,7150)     cm_rx*(2.0/bref), 
     &                   cm_ry*(2.0/cref), 
     &                   cm_rz*(2.0/bref)
 7150 format(' y  mom.  Cm |','    Cmp =',f11.6,
     &                        '    Cmq =',f11.6,
     &                        '    Cmr =',f11.6 )

      write(lu,7160) dir*cn_rx*(2.0/bref),
     &               dir*cn_ry*(2.0/cref), 
     &               dir*cn_rz*(2.0/bref)
 7160 format(' z'' mom.  Cn''|','    Cnp =',f11.6,
     &                        '    Cnq =',f11.6,
     &                        '    Cnr =',f11.6 )

      if(ncontrol.gt.0) then

      write(lu,8106) (dname(n), n, n=1, ncontrol)
 8106 format(/14x,20(4x,a12, ' d',i1,' '))
      write(lu,8107) (' ',n=1, ncontrol)
 8107 format( 14x,20(3x,a,'----------------'))

      write(lu,8110) (' ',n,clt_d(n), n=1, ncontrol)
 8110 format(' z'' force CL |' ,20(a,'  CLd',i1,' =',f11.6))

      write(lu,8120) (' ',n,cyt_d(n), n=1, ncontrol)
 8120 format(' y  force CY |'  ,20(a,'  CYd',i1,' =',f11.6))

      write(lu,8130) (' ',n,cdt_d(n), n=1, ncontrol)
 8130 format(' x'' force CD |' ,20(a,'  CDd',i1,' =',f11.6))

      write(lu,8140) (' ',n,dir*cmt_d(1,n), n=1, ncontrol)
 8140 format(' x'' mom.  Cl''|',20(a,'  Cld',i1,' =',f11.6))

      write(lu,8150) (' ',n,    cmt_d(2,n), n=1, ncontrol)
 8150 format(' y  mom.  Cm |'  ,20(a,'  Cmd',i1,' =',f11.6))

      write(lu,8160) (' ',n,dir*cmt_d(3,n), n=1, ncontrol)
 8160 format(' z'' mom.  Cn''|',20(a,'  Cnd',i1,' =',f11.6))

      write(lu,8170) (' ',n,    cdff_d(n), n=1, ncontrol)
 8170 format(' Trefftz drag|'  ,20(a,'CDffd',i1,' =',f11.6))

      write(lu,8180) (' ',n,  spanef_d(n), n=1, ncontrol)
 8180 format(' span eff.   |'  ,20(a,'   ed',i1,' =',f11.6))

      write(lu,*)
      write(lu,*)

      endif


      if(nvarjet.gt.0) then

      write(lu,8206) (jname(n), n, n=1, nvarjet)
 8206 format(/14x,20(4x,a12, ' j',i1,' '))
      write(lu,8207) (' ',n=1, nvarjet)
 8207 format( 14x,20(3x,a,'----------------'))

      write(lu,8210) (' ',n,clt_j(n), n=1, nvarjet)
 8210 format(' z'' force CL |'  ,20(a,'  CLj',i1,' =',f11.6))

      write(lu,8220) (' ',n,cyt_j(n), n=1, nvarjet)
 8220 format(' y  force CY |'   ,20(a,'  CYj',i1,' =',f11.6))

      write(lu,8230) (' ',n,cdt_j(n), n=1, nvarjet)
 8230 format(' x'' force CD |'   ,20(a,'  CDj',i1,' =',f11.6))

      write(lu,8240) (' ',n,dir*cmt_j(1,n), n=1, nvarjet)
 8240 format(' x'' mom.  Cl''|' ,20(a,'  Clj',i1,' =',f11.6))

      write(lu,8250) (' ',n,    cmt_j(2,n), n=1, nvarjet)
 8250 format(' y  mom.  Cm |'    ,20(a,'  Cmj',i1,' =',f11.6))

      write(lu,8260) (' ',n,dir*cmt_j(3,n), n=1, nvarjet)
 8260 format(' z'' mom.  Cn''|',20(a,'  Cnj',i1,' =',f11.6))

      write(lu,8270) (' ',n,    cdff_j(n), n=1, nvarjet)
 8270 format(' Trefftz drag|',20(a,'CDffj',i1,' =',f11.6))

      write(lu,8280) (' ',n,  spanef_j(n), n=1, nvarjet)
 8280 format(' span eff.   |',20(a,'   ej',i1,' =',f11.6))

      write(lu,*)
      write(lu,*)

      endif


      if(ndesign.gt.0) then

      write(lu,8306) (gname(n), n, n=1, ndesign)
 8306 format(/14x,20(4x,a12, ' g',i1,' '))
      write(lu,8307) (' ',n=1, ndesign)
 8307 format( 14x,20(3x,a,'----------------'))

      write(lu,8310) (' ',n,clt_g(n), n=1, ndesign)
 8310 format(' z'' force CL |'  ,20(a,'  CLg',i1,' =',f11.6))

      write(lu,8320) (' ',n,cyt_g(n), n=1, ndesign)
 8320 format(' y  force CY |'   ,20(a,'  CYg',i1,' =',f11.6))

      write(lu,8330) (' ',n,clt_g(n), n=1, ndesign)
 8330 format(' x'' force CD |'  ,20(a,'  CDg',i1,' =',f11.6))

      write(lu,8340) (' ',n,dir*cmt_g(1,n), n=1, ndesign)
 8340 format(' x'' mom.  Cl''|' ,20(a,'  Clg',i1,' =',f11.6))

      write(lu,8350) (' ',n,    cmt_g(2,n), n=1, ndesign)
 8350 format(' y  mom.  Cm |'    ,20(a,'  Cmg',i1,' =',f11.6))

      write(lu,8360) (' ',n,dir*cmt_g(3,n), n=1, ndesign)
 8360 format(' z'' mom.  Cn''|',20(a,'  Cng',i1,' =',f11.6))

      write(lu,8370) (' ',n,    cdff_g(n), n=1, ndesign)
 8370 format(' Trefftz drag|',20(a,'CDffg',i1,' =',f11.6))

      write(lu,8380) (' ',n,  spanef_g(n), n=1, ndesign)
 8380 format(' span eff.   |',20(a,'   eg',i1,' =',f11.6))

      write(lu,*)
      write(lu,*)

      endif

      if(cl_al .ne. 0.0) then
       xnp = xyzmom(1) - cref*cm_al/cl_al
       write(lu,8401) xnp
 8401  format(/' Neutral point  Xnp =', g14.6)
      endif

      if(abs(cr_rz*cn_be*(2.0/bref)) .gt. 0.0001) then
       bb = cr_be*cn_rz / (cr_rz*cn_be)
       write(lu,8402) bb 
 8402  format(/' Clb Cnr / Clr Cnb  =', f11.6,
     &    '    (  > 1 if spirally stable )')
      endif

      return
      end ! dermats



      subroutine dermatb(lu)
c---------------------------------------------------------
c     Outputs stability derivative matrix
c     for current alfa, beta, p,q,r, in body axes.
c     outputs it to logical unit lu
c---------------------------------------------------------
      include 'jvl.inc'
      character*50 satype
      real wbar_rx(3), wbar_rz(3), wbar_a(3)
      real drm(3)

      call getsa(lnasa_sa,satype,dir)

      drm(1) = xyzref(1) - xyzmom(1)
      drm(2) = xyzref(2) - xyzmom(2)
      drm(3) = xyzref(3) - xyzmom(3)

c---- set freestream velocity components from alpha, beta
      call ab2vbar

c---- calculate forces and sensitivities
      call aero

c---- coefficients and derivatives in geometry axes
      call gcoeffs(irun,drm)

      call outtot(lu)

      write(lu,7004)
 7004 format(/' geometry-axis derivatives...')

      write(lu,7006)
 7006 format(/14x, 4x,'  axial   vel. u',
     &             4x,' sideslip vel. v',
     &             4x,'  normal  vel. w'
     &       /14x, 4x,'----------------',
     &             4x,'----------------',
     &             4x,'----------------' )

      write(lu,7010)     cf_v(1,1),
     &               dir*cf_v(1,2),
     &                   cf_v(1,3)
 7010 format(' x force CX  |','    CXu =',f11.6,
     &                        '    CXv =',f11.6,
     &                        '    CXw =',f11.6 )

      write(lu,7020) dir*cf_v(2,1),
     &                   cf_v(2,2),
     &               dir*cf_v(2,3)
 7020 format(' y force CY  |','    CYu =',f11.6,
     &                        '    CYv =',f11.6,
     &                        '    CYw =',f11.6 )

      write(lu,7030)     cf_v(3,1),
     &               dir*cf_v(3,2),
     &                   cf_v(3,3)
 7030 format(' z force CZ  |','    CZu =',f11.6,
     &                        '    CZv =',f11.6,
     &                        '    CZw =',f11.6 )

      write(lu,7040)     cm_v(1,1),
     &               dir*cm_v(1,2),
     &                   cm_v(1,3)
 7040 format(' x mom.  Cl  |','    Clu =',f11.6,
     &                        '    Clv =',f11.6,
     &                        '    Clw =',f11.6 )

      write(lu,7050) dir*cm_v(2,1),
     &                   cm_v(2,2),
     &               dir*cm_v(2,3)
 7050 format(' y mom.  Cm  |','    Cmu =',f11.6,
     &                        '    Cmv =',f11.6,
     &                        '    Cmw =',f11.6 )

      write(lu,7060)     cm_v(3,1),
     &               dir*cm_v(3,2),
     &                   cm_v(3,3)
 7060 format(' z mom.  Cn  |','    Cnu =',f11.6,
     &                        '    Cnv =',f11.6,
     &                        '    Cnw =',f11.6 )

      write(lu,7106)
 7106 format(/14x, 4x,'    roll rate  p',
     &             4x,'   pitch rate  q',
     &             4x,'     yaw rate  r'
     &       /14x, 4x,'----------------',
     &             4x,'----------------',
     &             4x,'----------------' )

      write(lu,7110)     cf_w(1,1),
     &               dir*cf_w(1,2), 
     &                   cf_w(1,3)
 7110 format(' x force CX  |','    CXp =',f11.6,
     &                        '    CXq =',f11.6,
     &                        '    CXr =',f11.6 )

      write(lu,7120) dir*cf_w(2,1),
     &                   cf_w(2,2),
     &               dir*cf_w(2,3)
 7120 format(' y force CY  |','    CYp =',f11.6,
     &                        '    CYq =',f11.6,
     &                        '    CYr =',f11.6 )

      write(lu,7130)     cf_w(3,1),
     &               dir*cf_w(3,2),
     &                   cf_w(3,3)
 7130 format(' z force CZ  |','    CZp =',f11.6,
     &                        '    CZq =',f11.6,
     &                        '    CZr =',f11.6 )

      write(lu,7140)     cm_w(1,1),
     &               dir*cm_w(1,2),
     &                   cm_w(1,3)
 7140 format(' x mom.  Cl  |','    Clp =',f11.6,
     &                        '    Clq =',f11.6,
     &                        '    Clr =',f11.6 )

      write(lu,7150) dir*cm_w(2,1),
     &                   cm_w(2,2),
     &               dir*cm_w(2,3)
 7150 format(' y mom.  Cm  |','    Cmp =',f11.6,
     &                        '    Cmq =',f11.6,
     &                        '    Cmr =',f11.6 )

      write(lu,7160)     cm_w(3,1),
     &               dir*cm_w(3,2),
     &                   cm_w(3,3)
 7160 format(' z mom.  Cn  |','    Cnp =',f11.6,
     &                        '    Cnq =',f11.6,
     &                        '    Cnr =',f11.6 )

      if(ncontrol.gt.0) then

      write(lu,8106) (dname(n), n, n=1, ncontrol)
 8106 format(/14x,20(4x,a12, ' d',i1,' '))
      write(lu,8107) (' ',n=1, ncontrol)
 8107 format( 14x,20(3x,a,'----------------'))

      write(lu,8110) (' ',n,dir*cf_d(1,n), n=1, ncontrol)
 8110 format(' x force CX  |',20(a,'  CXd',i1,' =',f11.6))

      write(lu,8120) (' ',n,    cf_d(2,n), n=1, ncontrol)
 8120 format(' y force CY  |',20(a,'  CYd',i1,' =',f11.6))

      write(lu,8130) (' ',n,dir*cf_d(3,n), n=1, ncontrol)
 8130 format(' z force CZ  |',20(a,'  CZd',i1,' =',f11.6))

      write(lu,8140) (' ',n,dir*cm_d(1,n), n=1, ncontrol)
 8140 format(' x mom.  Cl  |',20(a,'  Cld',i1,' =',f11.6))

      write(lu,8150) (' ',n,    cm_d(2,n), n=1, ncontrol)
 8150 format(' y mom.  Cm  |',20(a,'  Cmd',i1,' =',f11.6))

      write(lu,8160) (' ',n,dir*cm_d(3,n), n=1, ncontrol)
 8160 format(' z mom.  Cn  |',20(a,'  Cnd',i1,' =',f11.6))

      write(lu,*)
      write(lu,*)

      endif


      if(nvarjet.gt.0) then

      write(lu,8206) (jname(n), n, n=1, nvarjet)
 8206 format(/14x,20(4x,a12, ' j',i1,' '))
      write(lu,8207) (' ',n=1, nvarjet)
 8207 format( 14x,20(3x,a,'----------------'))

      write(lu,8210) (' ',n,dir*cf_j(1,n), n=1, nvarjet)
 8210 format(' x force CX  |',20(a,'  CXj',i1,' =',f11.6))

      write(lu,8220) (' ',n,    cf_j(2,n), n=1, nvarjet)
 8220 format(' y force CY  |',20(a,'  CYj',i1,' =',f11.6))

      write(lu,8230) (' ',n,dir*cf_j(3,n), n=1, nvarjet)
 8230 format(' z force CZ  |',20(a,'  CZj',i1,' =',f11.6))

      write(lu,8240) (' ',n,dir*cm_j(1,n), n=1, nvarjet)
 8240 format(' x mom.  Cl  |',20(a,'  Clj',i1,' =',f11.6))

      write(lu,8250) (' ',n,    cm_j(2,n), n=1, nvarjet)
 8250 format(' y mom.  Cm  |',20(a,'  Cmj',i1,' =',f11.6))

      write(lu,8260) (' ',n,dir*cm_j(3,n), n=1, nvarjet)
 8260 format(' z mom.  Cn  |',20(a,'  Cnj',i1,' =',f11.6))

      write(lu,*)
      write(lu,*)

      endif


      if(ndesign.gt.0) then

      write(lu,8306) (gname(n), n, n=1, ndesign)
 8306 format(/14x,20(4x,a12, ' g',i1,' '))
      write(lu,8307) (' ',n=1, ndesign)
 8307 format( 14x,20(3x,a,'----------------'))

      write(lu,8310) (' ',n,dir*cf_g(1,n), n=1, ndesign)
 8310 format(' x force CX  |',20(a,'  CXg',i1,' =',f11.6))

      write(lu,8320) (' ',n,    cf_g(2,n), n=1, ndesign)
 8320 format(' y force CY  |',20(a,'  CYg',i1,' =',f11.6))

      write(lu,8330) (' ',n,dir*cf_g(3,n), n=1, ndesign)
 8330 format(' z force CZ  |',20(a,'  CZg',i1,' =',f11.6))

      write(lu,8340) (' ',n,dir*cm_g(1,n), n=1, ndesign)
 8340 format(' x mom.  Cl  |',20(a,'  Clg',i1,' =',f11.6))

      write(lu,8350) (' ',n,    cm_g(2,n), n=1, ndesign)
 8350 format(' y mom.  Cm  |',20(a,'  Cmg',i1,' =',f11.6))

      write(lu,8360) (' ',n,dir*cm_g(3,n), n=1, ndesign)
 8360 format(' z mom.  Cn  |',20(a,'  Cng',i1,' =',f11.6))

      write(lu,*)
      write(lu,*)

      endif

      return
      end ! dermatb



      subroutine dumpit(lu,nf,vbar,alpha, beta,
     &                  omegax, omegay, omegaz,
     &                  cfx, cfy, cfz, cmx, cmy, cmz)
c-------------------------------------------------------------
c     Writes flow condition header to logical unit lu
c-------------------------------------------------------------
      real vbar(nf),alpha(nf), beta(nf),
     &     omegax(nf), omegay(nf), omegaz(nf),
     &     cfx(nf), cfy(nf), cfz(nf), cmx(nf), cmy(nf), cmz(nf)

      do if=1, nf
        write(lu,2050) if, vbar(if)
        write(lu,2060) alpha(if), beta(if),
     &                 omegax(if), omegay(if), omegaz(if)
        write(lu,2070) cfx(if), cfy(if), cfz(if),
     &                 cmx(if), cmy(if), cmz(if)
      end do

 2050 format(/1x,'Flow condition', i3, '    Vinf =', f8.3)
 2060 format(/1x,6x,'alpha' ,7x,'beta',
     &           5x,'Omegax',5x,'Omegay',5x,'Omegaz' /1x,5f11.6)
 2070 format(/1x,8x,'CFx',8x,'CFy',8x,'CFz',
     &           8x,'CMx',8x,'CMy',8x,'CMz' /1x,6f11.6 / )

      return
      end ! dumpit


      subroutine getsa(lsa,satype,dir)
c-----------------------------------------------------------
c     Sets label string and multiplier dir corresponding to
c     standard flight-stability axes or geometry axes
c-----------------------------------------------------------
      logical lsa
      character*(*) satype

      if(lsa) then
       satype = 'Standard axis orientation,  X fwd, Z down'
       dir = -1.0
      else
       satype = 'Geometric axis orientation,  X aft, Z up  '
       dir =  1.0
      endif

      return
      end ! getsa



      subroutine fstrip(jstrip,
     &                  fs, ms,
     &                  caxial,cnorml, 
     &                  cl_strp,cd_strp,
     &                  clj_strp,cdj_strp,
     &                  clt_strp,cla_strp,
     &                  cmc4_strp,cmle_strp,
     &                  ccl_strp)
c-------------------------------------------------------------------
c     Returns various locally-normalized strip forces and moments
c
c   Input: 
c      jstrip     strip index
c
c   Output:
c      fs(3)      body axis strip forces
c      ms(3)      body axis strip moments about strip's c/4 point
c      caxial     force coeff in x direction
c      cnorml     force coeff normal to strip, along to ensy(j),ensz(j)
c      cl_strp    cl referenced to Vref,c (includes jet force)
c      cd_strp    cd referenced to Vref,c (includes jet force)
c      clj_strp   cl_jet referenced Vref,c (part of cl)
c      cdj_strp   cd_jet referenced Vref,c (part of cd)
c      clt_strp   cl referenced to span-normal local Vperp,c
c      cla_strp   cl referenced to local V,c
c      cmc4_strp  cm about strip's c/4 point, about spanwise axis
c      cmle_strp  cm about strip's le  point, about spanwise axis
c      ccl_strp   c cl
c-------------------------------------------------------------------
      include 'jvl.inc'

      real r(3), rc4(3), rrot(3)
      real vrot(3), vtot(3), vperp(3)
      real fs(3), ms(3)

      real spn(3)
      real udrag(3), udrag_u(3,numax),
     &     ulift(3), ulift_u(3,numax)
      real dxs(3), dxs_u(3,numax),
     &     dxsa, dxsa_u(numax)

      real drle(3)

      integer icrs(3), jcrs(3)
      data icrs / 2, 3, 1 / , jcrs / 3, 1, 2 /

      j = jstrip

      rho = rhos(j)
      cstrip = chord(j)
      sstrip = chord(j)*wstrip(j)

      rc4(1)  = rle(1,j) + 0.25*chord(j)
      rc4(2)  = rle(2,j)
      rc4(3)  = rle(3,j)

c---- define local strip lift and drag directions
c-      the "spanwise" vector is cross product of strip normal with x chordline 
      spn(1) =  0.0
      spn(2) =  ensz(j)
      spn(3) = -ensy(j)

      do k = 1, 3
        udrag(k) = vbar(k)
        do n = 1, numax
          udrag_u(k,n) = 0.
        enddo
        udrag_u(k,k) = 1.0
      enddo          

c---- strip lift direction is vector product of "stream" and spanwise vector
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
     &               + dxs(2)*dxs_u(2,n)
     &               + dxs(3)*dxs_u(3,n) ) / dxsa
       enddo
       do k = 1, 3
         ulift(k) = dxs(k)/dxsa
         do n = 1, numax
           ulift_u(k,n) = (dxs_u(k,n) - ulift(k)*dxsa_u(n))/dxsa
         enddo
       enddo

      endif

c---- translate moments from xyzref to c/4 point, to create local strip moments
      r(1) = xyzmom(1) - rc4(1)
      r(2) = xyzmom(2) - rc4(2)
      r(3) = xyzmom(3) - rc4(3)

      do k = 1, 3
        fs(k) = fsi(k,j) + fsv(k,j) + fsj(k,j)
        ms(k) = msi(k,j) + msv(k,j) + msj(k,j)
      enddo
      do k = 1, 3
        ic = icrs(k)
        jc = jcrs(k)
        ms(k) = ms(k) + r(ic)*fs(jc) - r(jc)*fs(ic)
      enddo

      caxial = fs(1) / (0.5*rho*sstrip)
      cnorml = (ensy(j)*fs(2)
     &        + ensz(j)*fs(3)) / (0.5*rho*sstrip)

      cmc4_strp = ( ensz(j)*ms(2)
     &            - ensy(j)*ms(3) ) / (0.5*rho*sstrip*cstrip)

c---- vector at chord reference point from rotation axes
      rrot(1) = xyzrefs(1,j) - xyzref(1)
      rrot(2) = xyzrefs(2,j) - xyzref(2)
      rrot(3) = xyzrefs(3,j) - xyzref(3)

c---- set total effective velocity = freestream + rotation
      call cross(rrot,wbar,vrot)
      vtot(1) = vbar(1) + vrot(1)
      vtot(2) = vbar(2) + vrot(2)
      vtot(3) = vbar(3) + vrot(3)

      vsq = vtot(1)**2 + vtot(2)**2 + vtot(3)**2
      if(vsq .eq. 0.0) then
       vsqi = 1.0
      else
       vsqi = 1.0 / vsq
      endif

c---- spanwise and perpendicular velocity components
      vspan = vtot(1)*ess(1,j) + vtot(2)*ess(2,j) + vtot(3)*ess(3,j)
      vperp(1) = vtot(1) - ess(1,j)*vspan
      vperp(2) = vtot(2) - ess(2,j)*vspan
      vperp(3) = vtot(3) - ess(3,j)*vspan

      vpsq = vperp(1)**2 + vperp(2)**2 + vperp(3)**2
      if(vpsq .eq. 0.0) then
       vpsqi = 1.0
      else
       vpsqi = 1.0 / vpsq
      endif

      cl_strp = ( ulift(1)*fs(1)
     &          + ulift(2)*fs(2)
     &          + ulift(3)*fs(3) ) / (0.5*rho*sstrip)
      cd_strp = ( udrag(1)*fs(1)
     &          + udrag(2)*fs(2)
     &          + udrag(3)*fs(3) ) / (0.5*rho*sstrip)

      clt_strp = cl_strp * vpsqi
      cla_strp = cl_strp * vsqi

c---- moment about strip le midpoint in direction of le segment
      r(1) = rc4(1) - rle(1,j)
      r(2) = rc4(2) - rle(2,j)
      r(3) = rc4(3) - rle(3,j)
      drle(1) = rle2(1,j) - rle1(1,j)
      drle(2) = rle2(2,j) - rle1(2,j)
      drle(3) = rle2(3,j) - rle1(3,j)

      if(idups(isurfs(j)).lt.0) then 
       drle(1) = -drle(1)
       drle(2) = -drle(2)
       drle(3) = -drle(3)
      endif
      dmag = sqrt(drle(1)**2+drle(2)**2+drle(3)**2)
      if(dmag.eq.0.0) then
       cmle_strp = 0.0
      else
       cmle_strp =
     &   ( drle(1)/dmag*(ms(1) + (r(2)*fs(3) - r(3)*fs(2)))
     &   + drle(2)/dmag*(ms(2) + (r(3)*fs(1) - r(1)*fs(3)))
     &   + drle(3)/dmag*(ms(3) + (r(1)*fs(2) - r(2)*fs(1))) )
     &   / (0.5*rho*sstrip*cstrip)
      endif

      ccl_strp = ( ulift(1)*(fsi(1,j)+fsj(1,j))
     &           + ulift(2)*(fsi(2,j)+fsj(2,j))
     &           + ulift(3)*(fsi(3,j)+fsj(3,j)) ) / (0.5*rho*wstrip(j))
c     ccl_strp = ( ensy(j)*fsi(2,j)
c    &           + ensz(j)*fsi(3,j) ) / (0.5*rho*wstrip(j))


c---- strip on surface has this strip trailing from it
      clj_strp = ( ulift(1)*fsj(1,j)
     &           + ulift(2)*fsj(2,j)
     &           + ulift(3)*fsj(3,j) ) / (0.5*rho*sstrip)
      cdj_strp = ( udrag(1)*fsj(1,j)
     &           + udrag(2)*fsj(2,j)
     &           + udrag(3)*fsj(3,j) ) / (0.5*rho*sstrip)

      return
      end ! fstrip

