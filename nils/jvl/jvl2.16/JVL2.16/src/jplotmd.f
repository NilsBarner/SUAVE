C***********************************************************************
C    Module:  jplotmd.f
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

      subroutine plotmd(azim, elev, tilt, rinv, keig, ir)
c---------------------------------------------------------------------
c     Plots eigenmode map and mode animations
c---------------------------------------------------------------------
      include 'jvl.inc'
      include 'jvlplt.inc'

      logical lmove, leview, lcpan
      logical linitview, lfirst
      save    linitview, lcpan
      logical error
      character*4 opt
      character*1 chkey

      real ang0(3), pos0(3)
      real ange(3), pose(3), vele(3), rote(3)
      real ang(3), pos(3), angp(3),dang(3)
      real tt(3,3), tt_ang(3,3,3),
     &     rt(3,3), rt_ang(3,3,3)

      real evr(jemax)

      real rinp(10)

      real*8 tcpu, tcpu0, tcpu1

c---- viewpoint changes (deg), zoom/unzoom, perspective scale factors
      data dazim, delev, zufac,  dpan , prat 
     &    / 5.0 , 5.0  , 1.5  ,  0.075, 1.1 /

c---- phase step and scale factor step for interactive phase plots
      data dphase, scalef / 5.0 , 1.25 /

      call getwinsize(xwind,ywind)
      ccs = 0.8*csize
      xplt0 = xabs2usr(pmarg)
      yplt0 = yabs2usr(ywind-ymarg-pmarg) - 1.2*ccs

c---- line y-spacing factor
      ysp = 2.0

c---- initialization for plot program variables
      if(lpltnew .or. (.not.lplot)) then 
        lpltnew = .false.
        linitview = .false.
      endif

      lfirst = .true.

c---- default is plot at zero baseline euler angles
      leview = .true.

c---- default is camera pans with aircraft
      lcpan = .true.

c---- no movie mode
      lmove = .false.

c---- find geometry limits
      call glims(gmin,gmax,.false.)

c---- set stuff for this mode
      vee = parval(ipvee,ir)
      refl = bref*unitl
      refv = vee

      evmin = (refv/refl) * 1.0e-5
      call evnorm(evec(1,keig,ir),esf,refl,refv)

      sigma = real(eval(keig,ir))
      omega = imag(eval(keig,ir))
      evmag = sqrt(sigma**2 + omega**2)

      frq = omega / (2.0*pi)
      if(sigma .eq. 0.0) then
       dampr = 0.
      else
       dampr = -sigma / evmag
      endif

c---- initial position
      pos0(1) = 0.
      pos0(2) = 0.
      pos0(3) = 0.

      ang0(1) = parval(ipphi,ir)*dtr
      ang0(2) = parval(ipthe,ir)*dtr
      ang0(3) = parval(ippsi,ir)*dtr

      xyzref(1) = parval(ipxcg,ir)
      xyzref(2) = parval(ipycg,ir)
      xyzref(3) = parval(ipzcg,ir)

c---- sat baseline velocities and rotation rates
      alfa = parval(ipalfa,ir)*dtr
      beta = parval(ipbeta,ir)*dtr
      call ab2vbar
      wbar(1) = parval(iprotx,ir)*2.0/bref
      wbar(2) = parval(iproty,ir)*2.0/cref
      wbar(3) = parval(iprotz,ir)*2.0/bref

c---- time step sign
      tsign = 1.0

c---- initial time, phase, and eigenvector scale
      timed = 0.

      pos(1) = pos0(1)
      pos(2) = pos0(2)
      pos(3) = pos0(3)
      ang(1) = ang0(1)
      ang(2) = ang0(2)
      ang(3) = ang0(3)

c---- set time step based on max euler angle change
      wrmax = max( abs(wbar(1)) , abs(wbar(2)) , abs(wbar(3)) )

c***************************************************
    1 continue
      write(*,2041) lcpan, lmwait
 2041 format(/'  ------------------------------'
     &       /'   L eft           R ight       '
     &       /'   U p             D own        '
     &       /'   C lear                       '
     &      //'   Z oom           N ormal size '
     &       /'   I ngress        O utgress    '
     &       /'   H ardcopy       A nnotate    '
     &      //'   P anning camera toggle: ', l1
     &       /'   W ait-for-clock toggle: ', l1
     &      //'  < > 0  mode play              '
     &       /'  - + 1  mode scale             '
     &       /'  S      mode sign change       '
     &       /'  { }    slower-mo / faster-mo  '
     &      //' Type in plot window:  command,  or  <space> to exit')

c   A B C D E F G H I J K L M N O P Q R S T U V W X Y Z
C   x   x x       x x     x   x x x   x x   x         x

c***************************************************
c---- setup view transformation 
    4 call viewinit(azim, elev, tilt, rinv)

c---- scale and project x,y,z orientation arrows uaxarw(....)
c-    also project x,y,z unit vectors uaxdir(..)
      arwlen = 0.2*bref

      rheadp = rhead*arwlen
      do ixyz = 1, 3
        do k = 1, 3
          do il = 1, nlinax
            uaxarwp(k,1,il,ixyz) = uaxarw(k,1,il,ixyz)*arwlen
            uaxarwp(k,2,il,ixyz) = uaxarw(k,2,il,ixyz)*arwlen
          enddo
        enddo
      enddo
      nlin = 3 * 2*nlinax
      naxd = 3   
      call viewproj(uaxarwp,nlin,uaxarwp)
      call viewproj(uaxdir ,naxd,uaxdirp)

      do iaxd = 1, naxd
        dirmag = sqrt( uaxdirp(1,iaxd)**2
     &               + uaxdirp(2,iaxd)**2
     &               + uaxdirp(3,iaxd)**2 )
        if(dirmag .gt. 0.0) then
         uaxdirp(1,iaxd) = uaxdirp(1,iaxd)/dirmag
         uaxdirp(2,iaxd) = uaxdirp(2,iaxd)/dirmag
         uaxdirp(3,iaxd) = uaxdirp(3,iaxd)/dirmag
        endif
      enddo

c***************************************************
c---- compute perturbed position at current or new time
    6 continue

      if(lmove) then
c----- advance or retreat in time
       timed = timed + tsign*dtimed*tmofac

       if(lcpan) then
c------ leave aircraft position unchanged (since camera is following it)

       else
c------ integrate velocities and angular rates to next movie frame time

c------ time interval in avl units
        dtime = dtimed*vee/unitl * tmofac

c------ set number of time steps for suitable accuracy (limit max angle change)
        damax = wrmax*dtime
        ntime = max( 1 , int(damax/0.02) )

c------ integrate over time interval  t = time..time+dtime
        do itime = 1, ntime
          dt = tsign * dtime/float(ntime)

c-------- predictor step, slopes evaluated at  t
          call rateki3(ang,rt,rt_ang)
          do k = 1, 3
            dang(k) = dt*( rt(k,1)*wbar(1)
     &                   + rt(k,2)*wbar(2)
     &                   + rt(k,3)*wbar(3) )
            angp(k) = ang(k) + dang(k)
          enddo

c-------- corrector step, slopes evaluated at  t + dt
          call rateki3(angp,rt,rt_ang)
          do k = 1, 3
            dangp = dt*( rt(k,1)*wbar(1)
     &                 + rt(k,2)*wbar(2)
     &                 + rt(k,3)*wbar(3) )
c---------- midpoint angles at  t + dt/2
            angp(k) = ang(k) + 0.25*(dang(k) + dangp)
          enddo

c-------- use midpoint-angle matrices
          call rateki3(angp,rt,rt_ang)
          call rotens3(angp,tt,tt_ang)

c-------- final integration step, using midpoint slopes
          do k = 1, 3
            pos(k) = pos(k) - dt*( tt(k,1)*vbar(1)
     &                           + tt(k,2)*vbar(2)
     &                           + tt(k,3)*vbar(3) )
            ang(k) = ang(k) + dt*( rt(k,1)*wbar(1)
     &                           + rt(k,2)*wbar(2)
     &                           + rt(k,3)*wbar(3) )
          enddo
        enddo
       endif
      endif

      ephase = omega*timed / dtr
      efac = esf*eigenf
      call evreal(evec(1,keig,ir),eval(keig,ir), efac,timed, evr)

c---- scale from standard (si or english) to avl units

      evr(jex) = evr(jex)/unitl
      evr(jey) = evr(jey)/unitl
      evr(jez) = evr(jez)/unitl

      evr(jeu) = evr(jeu)/vee
      evr(jev) = evr(jev)/vee
      evr(jew) = evr(jew)/vee

      evr(jep) = evr(jep)*unitl/vee
      evr(jeq) = evr(jeq)*unitl/vee
      evr(jer) = evr(jer)*unitl/vee

ccc   evr(jeph) = evr(jeph)
ccc   evr(jeth) = evr(jeth)
ccc   evr(jeps) = evr(jeps)

c      write(*,*)
c      write(*,*) 'xyz', evr(jex),evr(jey),evr(jez)
c      write(*,*) 'uvw', evr(jeu),evr(jev),evr(jew)
c      write(*,*) 'pqr', evr(jep),evr(jeq),evr(jer)
c      write(*,*) 'reh', evr(jeph),evr(jeth),evr(jeps)

c
c---- set final position, including eigenvector perturbation
      pose(1) = pos(1) + evr(jex)
      pose(2) = pos(2) + evr(jey)
      pose(3) = pos(3) + evr(jez)

c---- set final angles, including eigenvector perturbation
      ange(1) = ang(1) + evr(jeph)
      ange(2) = ang(2) + evr(jeth)
      ange(3) = ang(3) + evr(jeps)

      vele(1) = vbar(1) - evr(jeu)
      vele(2) = vbar(2) - evr(jev)
      vele(3) = vbar(3) - evr(jew)

      rote(1) = wbar(1) + evr(jep)
      rote(2) = wbar(2) + evr(jeq)
      rote(3) = wbar(3) + evr(jer)

      vmage = sqrt(vele(1)**2 + vele(2)**2 + vele(3)**2)
      v13   = sqrt(vele(1)**2              + vele(3)**2)
      alfae = atan2(  vele(3) , vele(1) )
      betae = atan2( -vele(2) , v13     )

c      call rateki3(ang,rt,rt_ang)
c      do k = 1, 3
c        dang(k) =    rt(k,1)*wbar(1)
c     &             + rt(k,2)*wbar(2)
c     &             + rt(k,3)*wbar(3)
c      enddo
c      write(1,1234) ( pos(k)/12.0, k=1, 3),
c     &              ( ang(k)/dtr , k=1, 3),
c     &              ( dang(k)*bref/dtr    , k=1, 3)
c      write(2,1234) (pose(k)/12.0, k=1, 3), (ange(k)/dtr, k=1, 3)
c 1234 format(1x,9f10.3)
c
c---- set limits for baseline position
      call rotens3(ang,tt,tt_ang)
      call grlims(vminp,vmaxp,.true. ,tt ,xyzref,pos )
      if(lfirst) then
c----- first frame... set plot limits directly
       xmin = vminp(1)
       ymin = vminp(2)
       xmax = vmaxp(1)
       ymax = vmaxp(2)
      else
c----- clip to new limits
       xmin = min(xmin,vminp(1))
       ymin = min(ymin,vminp(2))
       xmax = max(xmax,vmaxp(1))
       ymax = max(ymax,vmaxp(2))
      endif
      lfirst = .false.

c
c---- also enforce limits for perturbed position
      call rotens3(ange,tt,tt_ang)
      call grlims(vminp,vmaxp,.true. ,tt ,xyzref,pose)
      xmin = min(xmin,vminp(1))
      ymin = min(ymin,vminp(2))
      xmax = max(xmax,vmaxp(1))
      ymax = max(ymax,vmaxp(2))

      call offini

c***************************************************
c---- plot the baseline and displaced geometry
    8 continue
      call pvlini(title,azim,elev,tilt,version,lsvmov)

      xplt = xplt0
      yplt = yplt0
      call plchar(xplt,yplt,ccs,'Run  ',0.0,5)
      call plnumb(999.,yplt,ccs,float(ir),0.0,-1)

      yplt = yplt - ysp*ccs
      call plchar(xplt,yplt,ccs,'Mode ',0.0,5)
      call plnumb(999.,yplt,ccs,float(keig),0.0,-1)

      yplt = yplt - ysp*ccs
      call plchar(xplt,yplt,ccs,'f = ',0.0,4)
      call plnumb(999.,yplt,ccs,frq,0.0,4)
      call plchar(999.,yplt,ccs,' cycles/' ,0.0,8)
      call plchar(999.,yplt,ccs,uncht(1:nut),0.0,nut)

      yplt = yplt - ysp*ccs
      call plmath(xplt,yplt,ccs,'z = ',0.0,4)
      call plnumb(999.,yplt,ccs,dampr,0.0,6)

      yplt = yplt - ysp*ccs
      call plchar(xplt,yplt,ccs,'t = ',0.0,4)
      call plnumb(999.,yplt,ccs,timed ,0.0,3)
      call plchar(999.,yplt,ccs,uncht(1:nut),0.0,nut)

      yplt = yplt - ysp*ccs
      call plmath(xplt,yplt,ccs,'w = ',0.0, 4)
      call plchar(xplt,yplt,ccs,' t  ',0.0, 4)
      call plnumb(999.,yplt,ccs,ephase,0.0,-1)
      call plmath(999.,yplt,ccs,   '"',0.0, 1)

      yplt = yplt - 1.0*ccs

      yplt = yplt - ysp*ccs
      call plchar(xplt,yplt,ccs,'x = ',0.0, 4)
      call plnumb(999.,yplt,ccs,pose(1)*unitl,0.0,1)
      call plchar(999.,yplt,ccs,unchl,0.0,nul)

      yplt = yplt - ysp*ccs
      call plchar(xplt,yplt,ccs,'y = ',0.0, 4)
      call plnumb(999.,yplt,ccs,pose(2)*unitl,0.0,1)
      call plchar(999.,yplt,ccs,unchl,0.0,nul)

      yplt = yplt - ysp*ccs
      call plchar(xplt,yplt,ccs,'z = ',0.0, 4)
      call plnumb(999.,yplt,ccs,pose(3)*unitl,0.0,1)
      call plchar(999.,yplt,ccs,unchl,0.0,nul)

      yplt = yplt - 1.0*ccs

      yplt = yplt - ysp*ccs
      call plmath(xplt,yplt,ccs,'f = ',0.0, 4)
      call plnumb(999.,yplt,ccs,ange(1)/dtr,0.0,2)
      call plmath(999.,yplt,ccs,'"',0.0,1)

      yplt = yplt - ysp*ccs
      call plmath(xplt,yplt,ccs,'q = ',0.0, 4)
      call plnumb(999.,yplt,ccs,ange(2)/dtr,0.0,2)
      call plmath(999.,yplt,ccs,'"',0.0,1)

      yplt = yplt - ysp*ccs
      call plmath(xplt,yplt,ccs,'y = ',0.0, 4)
      call plnumb(999.,yplt,ccs,ange(3)/dtr,0.0,2)
      call plmath(999.,yplt,ccs,'"',0.0,1)

      yplt = yplt - 1.0*ccs

      yplt = yplt - ysp*ccs
      call plchar(xplt,yplt,ccs,'V = ',0.0, 4)
      call plnumb(999.,yplt,ccs,vmage*vee  ,0.0,2)
      call plchar(999.,yplt,ccs,unchv,0.0,nuv)

      yplt = yplt - ysp*ccs
      call plmath(xplt,yplt,ccs,'a = ',0.0, 4)
      call plnumb(999.,yplt,ccs,alfae/dtr,0.0,3)
      call plmath(999.,yplt,ccs,'"',0.0,1)

      yplt = yplt - ysp*ccs
      call plmath(xplt,yplt,ccs,'b = ',0.0, 4)
      call plnumb(999.,yplt,ccs,betae/dtr,0.0,3)
      call plmath(999.,yplt,ccs,'"',0.0,1)

c---- setup hidden line data
      call hidinite(.true., ange,pose,xyzref)

c---- plot baseline position
      call getcolor(icol0)
      call plotmode(ang ,pos ,xyzref,0)

c---- plot mode-displaced position
      call newcolorname('red')
      call plotmode(ange,pose,xyzref,1)

      call newcolor(icol0)
      call plflush

      if(lmove .and. lmwait) then
c----- waste time if in movie mode
       call seconds(tcpu1)
 9     continue
       call seconds(tcpu)
       dtcpu = tcpu - tcpu0
ccc       write(*,*) tcpu0, tcpu, dtcpu, dtimed
       if(dtcpu .lt. dtimed) go to 9

c       dtcpu1 = tcpu1 - tcpu0
c       write(*,*) 'Plot time frac:', dtcpu1/dtimed
      endif

ccc      call drawtoscreen
c---------------------------------------------
c---- set user command key
 10   continue

      chkey = ' '
      call getcursorxy(xx,yy,chkey)

      lmove = .false.
      if(chkey .eq. ' ') then
        return

c---- change observer's azimuth and elevation angles
      elseif(index('Ll',chkey).ne.0) then
        azim = azim + dazim
        if(azim .gt.  180.01) azim = azim - 360.0
        go to 4

      elseif(index('Rr',chkey).ne.0) then
        azim = azim - dazim
        if(azim .lt. -180.01) azim = azim + 360.0
        go to 4

      elseif(index('Uu',chkey).ne.0) then
        elev = elev + delev
        if(elev .gt.  180.01) elev = elev - 360.0
        go to 4

      elseif(index('Dd',chkey).ne.0) then
        elev = elev - delev
        if(elev .lt. -180.01) elev = elev + 360.0
        go to 4

      elseif(index('Cc',chkey).ne.0) then
        elev = 0.
        azim = 0.
        lfirst = .true.
        go to 4

      elseif(index('Zz',chkey).ne.0) then
c------ put center at current cursor location, increase magnification
        xoff = xx/sf + xoff
        yoff = yy/sf + yoff
        sf = sf*zufac
        xdif =    1.0/sf
        ydif = plotar/sf
        xoff = xoff - 0.5*xdif
        yoff = yoff - 0.5*ydif
        go to 8

      elseif(index('nn',chkey).ne.0) then
c------ force resetting of scale,offset factors
        lfirst = .true.
        go to 6

      elseif(index('Ii',chkey).ne.0) then
c------ increase perspective distortion
        if(rinv.eq.0.0) then
         rinv = 0.02/bref
        else
         rinv = rinv*prat
        endif
        write(*,8010) 1.0/rinv
        go to 4

      elseif(index('Oo',chkey).ne.0) then
c------ decrease perspective distortion
        if(rinv .lt. 0.02/bref) then
         rinv = 0.0
         write(*,*) '  Observer distance = infinity  (no perspective)'
        else
         rinv = rinv/prat
         write(*,8010) 1.0/rinv
        endif
        go to 4

      elseif(index('Aa',chkey).ne.0) then
c------ annotate plot
        call newpen(2)
        call annot(csize)
        go to 1

      elseif(index('Hh',chkey).ne.0) then
c------ write postscript version of current screen
        call replot(idevh)
        go to 1

c      elseif(index('Ee',chkey).ne.0) then
cc------ toggle viewing at true/zero euler angles
c        leview = .not.leview
c        if(leview) then
c         write(*,*) 'View at actual Euler angles'
c        else
c         write(*,*) 'View at zero Euler angles'
c        endif
c        go to 6

      elseif(index('Pp',chkey).ne.0) then
c------ toggle camera panning
        lcpan = .not.lcpan
        if(lcpan) then
         write(*,*) 'Camera will pan with aircraft'
        else
         write(*,*) 'Camera will not pan with aircraft'
        endif
        ephase = 0.
        timed = 0.
        pos(1) = pos0(1)
        pos(2) = pos0(2)
        pos(3) = pos0(3)
        ang(1) = ang0(1)
        ang(2) = ang0(2)
        ang(3) = ang0(3)
        eigenf = 1.
        lfirst = .true.
        go to 6

      elseif(index('Ww',chkey).ne.0) then
c------ toggle clock waiting
        lmwait = .not.lmwait
        if(lmwait) then
         write(*,*) 'Frame plotting will wait for clock'
        else
         write(*,*) 'Frame plotting will not wait for clock'
        endif
        go to 6

      elseif(index('0' ,chkey).ne.0) then
c------ set zero phase
ccc        ephase = 0.0
        timed = 0.
        pos(1) = pos0(1)
        pos(2) = pos0(2)
        pos(3) = pos0(3)
        ang(1) = ang0(1)
        ang(2) = ang0(2)
        ang(3) = ang0(3)
        go to 6

      elseif(index('1' ,chkey).ne.0) then
c------ set unity scale factor
        eigenf = 1.0
        go to 6

      elseif(index('<,',chkey).ne.0) then
c------ decrease phase by dphase amount
ccc        ephase = ephase - dphase
        lmove = .true.
        tsign = -1.0
        call seconds(tcpu0)
        go to 6

      elseif(index('>.',chkey).ne.0) then
c------ increase phase by dphase amount
ccc        ephase = ephase + dphase
        lmove = .true.
        tsign = +1.0
        call seconds(tcpu0)
        go to 6

      elseif(index('-_',chkey).ne.0) then
c------ decrease response vector coefficient
        eigenf = eigenf / scalef
        go to 6

      elseif(index('+=',chkey).ne.0) then
c------ increase response vector coefficient
        eigenf = eigenf * scalef
        go to 6

      elseif(index('{[',chkey).ne.0) then
c------ slow down time advance
        tmofac = tmofac / scalef
        if(abs(tmofac-1.0) .lt. 0.01) then
         tmofac = 1.0
         write(*,8020) tmofac, ' (real time)'
        else
         write(*,8020) tmofac
        endif
        go to 6

      elseif(index('}]',chkey).ne.0) then
c------ speed up time advance
        tmofac = tmofac * scalef
        if(abs(tmofac-1.0) .lt. 0.01) then
         tmofac = 1.0
         write(*,8020) tmofac, ' (real time)'
        else
         write(*,8020) tmofac
        endif
        go to 6

      elseif(index('Ss',chkey).ne.0) then
c------ change mode sign
        eigenf = -eigenf
        go to 6

      else
        write(*,*) '* Key command not recognized'
        go to 1

      endif
 8010 format('  Observer distance =', f11.3,1x,a)
 8020 format('  Time advance factor =', f11.3,1x,a)

      end ! plotmd




      subroutine evreal(evec,eval,efac,timed,evr)
      complex evec(*), eval
      real evr(*)
      include 'jindex.inc'

      complex evfac

c
c---- set complex mode scale factor
      evfac = cexp( eval*timed ) * efac

c---- real part of phsyically-scaled scaled eigenmode
      nj = jetot
      do j = 1, nj
        evr(j) = real( evfac * evec(j) )
      enddo

      end ! evreal



      subroutine evnorm(ev,esf,refl,refv)
      complex ev(*)

      include 'jindex.inc'

      data efrac / 0.5 /

      refw = refv/refl
      refa = 0.5

      rmax = 0.

c---- test x,y,z changes
      rmax = max( rmax , cabs( ev(jex)/refl ) )
      rmax = max( rmax , cabs( ev(jey)/refl ) )
      rmax = max( rmax , cabs( ev(jez)/refl ) )

c---- test u,v,w changes
      rmax = max( rmax , cabs( ev(jeu)/refv ) )
      rmax = max( rmax , cabs( ev(jev)/refv ) )
      rmax = max( rmax , cabs( ev(jew)/refv ) )

c---- test p,q,r changes
      rmax = max( rmax , cabs( ev(jep)/refw ) )
      rmax = max( rmax , cabs( ev(jeq)/refw ) )
      rmax = max( rmax , cabs( ev(jer)/refw ) )

c---- test euler angle changes
      rmax = max( rmax , cabs( ev(jeph)/refa ) )
      rmax = max( rmax , cabs( ev(jeth)/refa ) )
      rmax = max( rmax , cabs( ev(jeps)/refa ) )

      if(rmax .eq. 0.0) then
        esf = 1.0
      else
        esf = efrac/rmax
      endif

      return
      end ! evnorm



      subroutine plotmode(ang,pos,xyzr,ipat)
      include 'jvl.inc'
      include 'jvlplt.inc'

      real ang(3), pos(3), xyzr(3)

      real tt(3,3), tt_ang(3,3,3)

      real pts_lines(3,2,nvmax), id_lines(nvmax),
     &     pts_lproj(3,2,nvmax)

      data id_lines / nvmax * 0.0 /

c---- number of line segments for drawing a body section "hoop"
      data nhoop / 18 /

      include 'masks.inc'

      if(ipat.eq.0) then
       lmout = lmask2
       lmchd = lmask3
       ipout = 3
       ipchd = 1
      else
       lmout = lmask0
       lmchd = lmask1
       ipout = 4
       ipchd = 1
      endif

      call rotens3(ang,tt,tt_ang)

c---- go over surfaces
      do 100 n = 1, nsurf
c------ do only surfaces which are to be plotted
        if(.not.lpltsurf(n)) go to 100

        j1 = jfrst(n)
        jn = jfrst(n) + nj(n)-1
        jinc = 1
        if(idups(n).lt.0) then
         j1 = jfrst(n) + nj(n)-1
         jn = jfrst(n)
         jinc = -1
        endif

c------ chordlines plotted
        if(lchordplt) then
         ip = 0
         do j = j1, jn, jinc
           ip = ip+1
           pts_lines(1,1,ip) = rle1(1,j)
           pts_lines(2,1,ip) = rle1(2,j)
           pts_lines(3,1,ip) = rle1(3,j)
           pts_lines(1,2,ip) = rle1(1,j) + chord1(j)
           pts_lines(2,2,ip) = rle1(2,j)
           pts_lines(3,2,ip) = rle1(3,j)
           id_lines(ip) = j
           call tetran(pts_lines(1,1,ip),tt,xyzr,pos)
           call tetran(pts_lines(1,2,ip),tt,xyzr,pos)
         end do
         ip = ip+1
         pts_lines(1,1,ip) = rle2(1,jn)
         pts_lines(2,1,ip) = rle2(2,jn)
         pts_lines(3,1,ip) = rle2(3,jn)
         pts_lines(1,2,ip) = rle2(1,jn) + chord2(jn)
         pts_lines(2,2,ip) = rle2(2,jn)
         pts_lines(3,2,ip) = rle2(3,jn)
         id_lines(ip) = jn
         call tetran(pts_lines(1,1,ip),tt,xyzr,pos)
         call tetran(pts_lines(1,2,ip),tt,xyzr,pos)

         nlines = ip
         nproj = 2*nlines
         call viewproj(pts_lines,nproj,pts_lproj)

         call newpat(lmchd)
         call newpen(ipchd)
         call plotlines(nlines,pts_lproj,id_lines)
c         do ip = 1, nlines
c           call plot(pts_lines(1,1,ip),pts_lines(2,1,ip),3)
c           call plot(pts_lines(1,2,ip),pts_lines(2,2,ip),2)
c         enddo
        endif

c--- surface boundary lines plotted
        ip = 0
        ip = ip+1
        pts_lines(1,1,ip) = rle1(1,j1)
        pts_lines(2,1,ip) = rle1(2,j1)
        pts_lines(3,1,ip) = rle1(3,j1)
        pts_lines(1,2,ip) = rle1(1,j1) + chord1(j1)
        pts_lines(2,2,ip) = rle1(2,j1)
        pts_lines(3,2,ip) = rle1(3,j1)
        id_lines(ip) = j1
        call tetran(pts_lines(1,1,ip),tt,xyzr,pos)
        call tetran(pts_lines(1,2,ip),tt,xyzr,pos)

        ip = ip+1
        pts_lines(1,1,ip) = rle2(1,jn)
        pts_lines(2,1,ip) = rle2(2,jn)
        pts_lines(3,1,ip) = rle2(3,jn)
        pts_lines(1,2,ip) = rle2(1,jn) + chord2(jn)
        pts_lines(2,2,ip) = rle2(2,jn)
        pts_lines(3,2,ip) = rle2(3,jn)
        id_lines(ip) = jn
        call tetran(pts_lines(1,1,ip),tt,xyzr,pos)
        call tetran(pts_lines(1,2,ip),tt,xyzr,pos)

        do j = j1, jn, jinc
          ip = ip+1
          pts_lines(1,1,ip) = rle1(1,j)
          pts_lines(2,1,ip) = rle1(2,j)
          pts_lines(3,1,ip) = rle1(3,j)
          pts_lines(1,2,ip) = rle2(1,j)
          pts_lines(2,2,ip) = rle2(2,j)
          pts_lines(3,2,ip) = rle2(3,j)
          id_lines(ip) = j
          call tetran(pts_lines(1,1,ip),tt,xyzr,pos)
          call tetran(pts_lines(1,2,ip),tt,xyzr,pos)

          ip = ip+1
          pts_lines(1,1,ip) = rle1(1,j) + chord1(j)
          pts_lines(2,1,ip) = rle1(2,j)
          pts_lines(3,1,ip) = rle1(3,j)
          pts_lines(1,2,ip) = rle2(1,j) + chord2(j)
          pts_lines(2,2,ip) = rle2(2,j)
          pts_lines(3,2,ip) = rle2(3,j)
          id_lines(ip) = j
          call tetran(pts_lines(1,1,ip),tt,xyzr,pos)
          call tetran(pts_lines(1,2,ip),tt,xyzr,pos)
        end do

        nlines = ip
        nproj = 2*nlines
        call viewproj(pts_lines,nproj,pts_lproj)
        call newpat(lmout)
        call newpen(ipout)
        call plotlines(nlines,pts_lproj,id_lines)
c        do ip = 1, nlines
c          call plot(pts_lines(1,1,ip),pts_lines(2,1,ip),3)
c          call plot(pts_lines(1,2,ip),pts_lines(2,2,ip),2)
c        enddo
 100  continue

c---- go over bodies
      do 300 n = 1, nbody
c------ do only bodies which are to be plotted
        if(.not.lpltbody(n)) go to 300

        l1 = lfrst(n)
        ln = llast(n)
        linc = 1

c------ source lines plotted
        ip = 0
        do l = l1, ln-linc, linc
          ip = ip+1
          pts_lines(1,1,ip) = rl(1,l)
          pts_lines(2,1,ip) = rl(2,l)
          pts_lines(3,1,ip) = rl(3,l)
          pts_lines(1,2,ip) = rl(1,l+linc)
          pts_lines(2,2,ip) = rl(2,l+linc)
          pts_lines(3,2,ip) = rl(3,l+linc)
          call tetran(pts_lines(1,1,ip),tt,xyzr,pos)
          call tetran(pts_lines(1,2,ip),tt,xyzr,pos)
        end do

        nlines = ip
        nproj = 2*nlines
        call viewproj(pts_lines,nproj,pts_lproj)

        call newpat(lmchd)
        call plotlines(nlines,pts_lproj,id_lines)
c        do ip = 1, nlines
c          call plot(pts_lines(1,1,ip),pts_lines(2,1,ip),3)
c          call plot(pts_lines(1,2,ip),pts_lines(2,2,ip),2)
c        enddo

c------ body hoops plotted
ccc        if(lchordplt) then
         ip = 0
         do l = l1, ln, linc
c--------- set "horizontal" and "vertical" body section unit vectors
           if(linc.gt.0) then
            lm = max(l-1,l1)
            lp = min(l+1,ln)
           else
            lm = max(l-1,ln)
            lp = min(l+1,l1)
           endif
           dxl = rl(1,lp) - rl(1,lm)
           dyl = rl(2,lp) - rl(2,lm)
           dzl = rl(3,lp) - rl(3,lm)
           dsl = sqrt(dxl**2 + dyl**2 + dzl**2)
           if(dsl.ne.0.0) then
            dxl = dxl/dsl
            dyl = dyl/dsl
            dzl = dzl/dsl
           endif
c     
           uhx = -dyl
           uhy =  dxl
           uhz = 0.

           uvx = dyl*uhz - dzl*uhy
           uvy = dzl*uhx - dxl*uhz
           uvz = dxl*uhy - dyl*uhx

           do k = 1, nhoop
             thk1 = 2.0*pi*float(k-1)/float(nhoop)
             thk2 = 2.0*pi*float(k  )/float(nhoop)

             hx1 = uhx*cos(thk1)
             hy1 = uhy*cos(thk1)
             hz1 = uhz*cos(thk1)

             vx1 = uvx*sin(thk1)
             vy1 = uvy*sin(thk1)
             vz1 = uvz*sin(thk1)

             hx2 = uhx*cos(thk2)
             hy2 = uhy*cos(thk2)
             hz2 = uhz*cos(thk2)

             vx2 = uvx*sin(thk2)
             vy2 = uvy*sin(thk2)
             vz2 = uvz*sin(thk2)

             ip = ip+1
             pts_lines(1,1,ip) = rl(1,l) + hx1*radl(l) + vx1*radl(l)
             pts_lines(2,1,ip) = rl(2,l) + hy1*radl(l) + vy1*radl(l)
             pts_lines(3,1,ip) = rl(3,l) + hz1*radl(l) + vz1*radl(l)
             pts_lines(1,2,ip) = rl(1,l) + hx2*radl(l) + vx2*radl(l)
             pts_lines(2,2,ip) = rl(2,l) + hy2*radl(l) + vy2*radl(l)
             pts_lines(3,2,ip) = rl(3,l) + hz2*radl(l) + vz2*radl(l)
             call tetran(pts_lines(1,1,ip),tt,xyzr,pos)
             call tetran(pts_lines(1,2,ip),tt,xyzr,pos)
           enddo
         end do

         nlines = ip
         nproj = 2*nlines
         call viewproj(pts_lines,nproj,pts_lproj)

         call newpat(lmout)
         call newpen(ipout)
         call plotlines(nlines,pts_lproj,id_lines)
c         do ip = 1, nlines
c           call plot(pts_lines(1,1,ip),pts_lines(2,1,ip),3)
c           call plot(pts_lines(1,2,ip),pts_lines(2,2,ip),2)
c         enddo
ccc        endif
 300  continue

      call newpat(lmask0)
      return
      end ! plotmode


      subroutine plemap(xorg0,yorg0, 
     &                  irun1,irun2,
     &                  keview,irview,
     &                  lproot,lprnum,overlay,
     &                  tmin,tmax,tdel,tfac,
     &                  fmin,fmax,fdel,ffac )
      include 'jvl.inc'
      include 'jvlplt.inc'

      logical lproot(jemax,*), lprnum(jemax,*), overlay

      logical lconvt(nrmax), lline, lbase
      logical lpgrid

      real parv(ipmax,nrmax)

      cs  = 0.55*csize
      csl = 0.70*csize
      csn = 0.75*csize
      csp = 0.55*csize

      dxlab = 2.0*csp
      dylab = 2.3*csp

c---- data root overlay symbol size and type
      css = 0.55*cs
      isymb = 3

      lpgrid = .true.

c---- initialize plot window and move origin away from lower left corner
      call plopen(scrnfrac,ipslu,idev)

      if(lcrev) then
       call bgfill
      endif

      call plot(1.0,0.75,-3)
      call newfactor(plsize)

      call getcolor(icol0)

c---- aspect ratio, lower-left plot location
      xorg = xorg0
      yorg = yorg0

c---- draw plot
      call plroot(xorg,yorg, 
     &            irun1, irun2, ircolor, 
     &            jemax,neigen,eval,lproot,lprnum,
     &            uncht(1:nut), plotar,cs, lpgrid, 
     &            tmin,tmax,tdel,tfac,fmin,fmax,fdel,ffac)

cc---- plot "o" over shift location
c      call newpen(1)
c      xc1 = xorg + (real(zshift)-tmin)*tfac
c      yc1 = yorg + (imag(zshift)-fmin)*ffac
c      call plsymb(xc1,yc1,cs,1,0.0,0)

      if(irview.ge.1 .and. irview.le.nrun) then
       ir = irview
       if(keview.gt.0 .and. keview.le.neigen(ir)) then
c------ plot "+" over currently-examined root, if any
        call newpen(1)
        xc1 = xorg + (real(eval(keview,ir))-tmin)*tfac
        yc1 = yorg + (imag(eval(keview,ir))-fmin)*ffac
        call plsymb(xc1,yc1,2.0*cs,3,0.0,0)
       endif
      endif

      if(overlay) then
c----- data exists... overlay it
       call getcolor(icol0)

       call newpen(2)
       xorg = xorg0
       yorg = yorg0

       do ir = 1, nrmax
         do keig = 1, neigendat(ir)
           call newcolor(ircolor(ir))
           xplt = xorg + (real(evaldat(keig,ir))-tmin)*tfac
           yplt = yorg + (imag(evaldat(keig,ir))-fmin)*ffac
           call plsymb(xplt,yplt,css,isymb,0.0,0)
         enddo
       enddo
       call newpen(1)
       call newcolor(icol0)
      endif

      call plflush

c---- put up parameter list
      do ir = irun1, irun2
        do ip = 1, ipmax
          parv(ip,ir) = parval(ip,ir)
          if(abs(parv(ip,ir)) .lt. 1.0e-8) parv(ip,ir) = 0.
        enddo
      enddo

      xplt = xorg
      yplt = yorg + (fmax-fmin)*ffac

      call pltpar(xplt,yplt, irun1, irun2,
     &            parv, lppar,
     &            ircolor, csp, dxlab, dylab )

      yplt = yplt + 3.0*csp
      call bstrip(title,nti)
      call newpen(3)
      call plchar(xorg,yplt,csn,title,0.0,nti)

      call plflush
      return
      end ! plemap



      subroutine tetran(r,tt,rref,dr)
      real r(3), tt(3,3), rref(3), dr(3)

      real rb(3)

      rb(1) = r(1) - rref(1)
      rb(2) = r(2) - rref(2)
      rb(3) = r(3) - rref(3)

      do k = 1, 3
        r(k) = tt(k,1)*rb(1)
     &       + tt(k,2)*rb(2)
     &       + tt(k,3)*rb(3) + rref(k) + dr(k)
      enddo

      return
      end ! tetran

