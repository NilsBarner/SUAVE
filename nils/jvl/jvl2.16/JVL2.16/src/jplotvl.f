C***********************************************************************
C    Module:  jplotvl.f
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

      subroutine plotvl(azim, elev, tilt, rinv)
c------------------------------------------------------------------------
c     Plots geometry and loading results for vortex lattice geometry plots:
c
c       surfaces     (straight tapered panels of chordwise strips of vortices)
c       strips       (chordwise strips of vortices)
c       vortex legs  (bound vortex legs)
c       control pts  (vortex element control points)
c       camber slope (camber of each strip - used for establishing bc's)
c       hinge lines  (surface deflection axis and deflected surface outline)
c       strip loading (chordwise plot of vortex loading on each strip)
c------------------------------------------------------------------------
      include 'jvl.inc'
      include 'jvlplt.inc'

      logical lkeys, lcolhc, error
      logical linitview, ldebug
      save    linitview, ldebug
      character*80 line
      character*4 opt
      character*1 ckey

      real rinp(10), arg(5)

c---- viewpoint changes (deg), zoom/unzoom, perspective scale factors
      data dazim, delev, zufac, dpan , prat 
     &     / 5.0 , 5.0 , 1.25 , 0.075, 1.1 /

    1 format(a)

c***************************************************
c---- initialization for plot program variables
      if(lpltnew .or. (.not.lplot)) then 
        lpltnew = .false.
        linitview = .false.
        ldebug = .false.
      endif

      lkeys = .false.

c---- find geometry and axis limits
      call glims(gmin,gmax,.false.)
      call axlims

c***************************************************
c---- setup view transformation 
    4 continue
      call viewinit(azim, elev, tilt, rinv)

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

c---- setup hidden line data
      call hidinit(.true.)

c---- find plot limits
      call glims(gminp,gmaxp,.true.)
      xmin = gminp(1)
      ymin = gminp(2)
      xmax = gmaxp(1)
      ymax = gmaxp(2)

c---- set plot offsets and scale
      call offini

      linitview = .false.

c
c***************************************************
c---- start geometry plot
    6 continue
      call pvlini(title,azim,elev,tilt,version,.false.)

c---- debug: check the triangle array
      if(ldebug) then
       call getcolor(icol0)
       call newcolorname('red')
       do i = 1, ntri
        call pltpoly(tri(1,i),3)
       end do
       call newcolor(icol0)
       call plflush
cc      pause
      endif

c---- plot the selected geometry items
      call plotgeom
      call plflush

 7    continue

      if(lkeys) then
c------ process keystroke-mode commands, skipping normal menu loop farther below
        ckey = ' '
        call getcursorxy(xx,yy,ckey)

        if(ckey.eq.' ') then
c------- exit to normal menu loop
         lkeys = .false.

        elseif(index('Ll',ckey).ne.0) then
         azim = azim + dazim
         if(azim .gt.  180.01) azim = azim - 360.0
         go to 4

        elseif(index('Rr',ckey).ne.0) then
         azim = azim - dazim
         if(azim .lt. -180.01) azim = azim + 360.0
         go to 4

        elseif(index('Uu',ckey).ne.0) then
         elev = elev + delev
         if(elev .gt.  90.01) then
          elev = 90.0
          write(*,*) 'Elevation angle is limited to +/- 90 deg'
         endif
         go to 4

        elseif(index('Dd',ckey).ne.0) then
         elev = elev - delev
         if(elev .lt. -90.01) then
          elev = -90.0
          write(*,*) 'Elevation angle is limited to +/- 90 deg'
         endif
         go to 4

c        elseif(index('Tt',ckey).ne.0) then
c         tilt = tilt + 5.0
c         if(tilt .gt.  180.01) tilt = tilt - 360.0
c         go to 4
cc
c        elseif(index('Ss',ckey).ne.0) then
c         tilt = tilt - 5.0
c         if(tilt .lt. -180.01) tilt = tilt + 360.0
c         go to 4

        elseif(index('Cc',ckey).ne.0) then
         azim = 0.
         elev = 0.
         tilt = 0.
         go to 4

        elseif(index('Zz',ckey).ne.0) then
         xoff = xx/sf + xoff
         yoff = yy/sf + yoff
         sf = sf*zufac
         xdif =    1.0/sf
         ydif = plotar/sf
         xoff = xoff - 0.5*xdif
         yoff = yoff - 0.5*ydif
         go to 6

        elseif(index('Ee',ckey).ne.0) then
         xoff = xx/sf + xoff
         yoff = yy/sf + yoff
         sf = sf/zufac
         xdif =    1.0/sf
         ydif = plotar/sf
         xoff = xoff - 0.5*xdif
         yoff = yoff - 0.5*ydif
         go to 6

        elseif(index('Pp',ckey).ne.0) then
         dx = xx - 0.5
         dy = yy - 0.5*plotar

         xoff = xoff + dpan*dx/sf
         yoff = yoff + dpan*dy/sf
         go to 6

        elseif(index('Nn',ckey).ne.0) then
         go to 4

        elseif(index('Ii',ckey).ne.0) then
         if(rinv.eq.0.0) then
          rinv = 0.02/bref
         else
          rinv = rinv*prat
         endif
         write(*,8010) unitl/rinv, unchl(1:nul)
         go to 4

        elseif(index('Oo',ckey).ne.0) then
         if(rinv .lt. 0.02/bref) then
          rinv = 0.0
          write(*,*) '  Observer distance = infinity  (no perspective)'
         else
          rinv = rinv/prat
          write(*,8010) unitl/rinv, unchl(1:nul)
         endif
         go to 4

        elseif(index('Aa',ckey).ne.0) then
c------- annotate
         call annot(csize)
         write(*,2040)
         go to 7

        elseif(index('Hh',ckey).ne.0) then
c------- make hardcopy
         if(lplot) call plend
         lplot = .false.
         call replot(idevh)
         go to 6

        else
         write(*,*)
         write(*,*) '"', ckey,'"  key not recognized'
         write(*,2040)

        endif

        go to 4
      endif

 8010 format('   Observer distance =', g11.3,1x,a)

c***************************************************
c---- top of normal user interaction loop

    8 write(*,2010)
      write(*,2020) lchordplt , lcambplt ,
     &              lcntlplt  , lwakeplt,
     &              lboundplt , lnrmlplt,
     &              lloadplt  , lwaklplt,
     &              laxesplt  , lvelcplt,
     &              ldiskplt  , lwsegplt,
     &              lhidden
      write(*,2030) 
      read(*,1) opt 

 2010 format(
     &  /' ========================================='
     &  /'  K eystroke mode      V iewpoint        '
     &  /'  A nnotate plot       O ptions          '
     &  /'  H ardcopy plot       S elect surfaces  '
     &  /'  Z oom                U nzoom           ')
 2020 format(
     &  /'  CH ordline  ',l2,7x,'CA mber         ',l2,
     &  /'  CN tlpoint  ',l2,7x,'TR ailing legs  ',l2,
     &  /'  BO ound leg ',l2,7x,'NO rmal vector  ',l2,
     &  /'  LO ading    ',l2,7x,'WA ke loading   ',l2,
     &  /'  AX es       ',l2,7x,'VC ntl velocity ',l2,
     &  /'  JE t disk   ',l2,7x,'VV ort velocity ',l2,
     &  /'  HI de lines ',l2                           )

 2030 format(/' Geometry plot command: ', $)

      call touper(opt) 
      if(opt.eq.'  ' .or. opt.eq.'Q ') then
       call plend
       return
      endif

      if(opt.eq.'K ') then
c----- enter keystroke mode
       lkeys = .true.
       write(*,2040)
 2040  format(
     &   /'  ------------------------------------------------'
     &   /'  Type keys in graphics window ...'
     &  //'    L eft               R ight        (azimuth  ) '
     &   /'    U p                 D own         (elevation) '
     &   /'    C lear'
     &  //'    P an from cursor    Z oom on cursor           '
     &   /'    I ngress            E xpand on cursor         '
     &   /'    O utgress           N ormal size              '
     &  //'    H ardcopy           A nnotate plot            '
     &  //'  ... <space> to exit  '
     &   /'  ------------------------------------------------')
       go to 6
      endif

c===================================================
      if(opt.eq.'V ') then
 20     write(*,2050) azim, elev
 2050   format(
     &   /'  Enter viewpoint azimuth,elevation angles:',2f7.0)
        rinp(1) = azim
        rinp(2) = elev
        ninp = 2
        call readr(ninp,rinp,error)
        if(error) go to 20

        azim = rinp(1)
        elev = rinp(2)
        if(azim .gt.  180.01) azim = azim - 360.0
        if(azim .lt. -180.01) azim = azim + 360.0
        if(elev .gt.   90.01) elev =  90.0
        if(elev .lt.  -90.01) elev = -90.0
        go to 4
c===================================================
      else if(opt.eq.'CH') then
c------ toggle chordline plotting
        lchordplt = .not.lchordplt

c===================================================
      else if(opt.eq.'CA') then
c------ toggle camber
        lcambplt = .not.lcambplt

c===================================================
      else if(opt.eq.'BO') then
c------ toggle plotting of bound legs
        lboundplt = .not.lboundplt

c===================================================
      else if(opt.eq.'CN') then
c------ toggle control point plotting
        lcntlplt = .not.lcntlplt

c===================================================
      else if(opt.eq.'TR') then
c------ toggle wake plotting
        lwakeplt = .not.lwakeplt

c===================================================
      else if(opt.eq.'LO') then
c------ toggle plotting of surface loading
        lloadplt = .not.lloadplt

c===================================================
      else if(opt.eq.'WA') then
c------ toggle plotting of wake loading
        lwaklplt = .not.lwaklplt

c===================================================
      else if(opt.eq.'NO') then
c------ toggle plotting of surface normals
        lnrmlplt = .not.lnrmlplt

c===================================================
      else if(opt.eq.'VC') then
c------ toggle plotting of surface c.p. velocities
        lvelcplt = .not.lvelcplt

c===================================================
      else if(opt.eq.'VV') then
c------ toggle plotting of surface vortex  velocities
        lwsegplt = .not.lwsegplt

c===================================================
      else if(opt.eq.'JE') then
c------ toggle plotting of actuator disks
        ldiskplt = .not.ldiskplt

c===================================================
      else if(opt.eq.'AX') then
c------ toggle plotting of x,y,z axes
        laxesplt = .not.laxesplt

c===================================================
       else if(opt.eq.'HI') then
c------ toggle hidden line
        lhidden = .not.lhidden

c===================================================
      else if(opt.eq.'RE') then
c------ toggle plotting of x,y,z reference location
        lrrefplt = .not.lrrefplt

c===================================================
cc    else if(opt.eq.'HI') then
ccc---- toggle plotting of hinge axes 
cc      lhingeplt = .not.lhingeplt

c===================================================
      else if(opt.eq.'DE') then
c------ toggle debug
        ldebug = .not.ldebug

c===================================================
      elseif(opt.eq.'A ') then
c----- annotate
       if(lplot) then
        call annot(csize)
       else
        write(*,*) 'No active plot'
       endif

c===================================================
      elseif(opt.eq.'H ') then
c----- make hardcopy
       if(lplot) call plend
       lplot = .false.
       call replot(idevh)
       go to 8

c===================================================
      elseif(opt.eq.'Z ') then
c----- zoom in on the plot
       if (lplot) then
        call offget   
       else
        write(*,*) 'No active plot'
       endif

c===================================================
      elseif(opt.eq.'U ') then
c----- reset to default offset and scaling factors
       call offini   

c===================================================
      elseif(opt.eq.'S ') then
c----- set surfaces to plot
       call selsurf(azim, elev, tilt, rinv)
       go to 4

c===================================================
      elseif(opt.eq.'O ') then
c----- display options 
 50    continue
       write(*,*) ' '
       rob = 999999.
       if(rinv.ne.0.) rob = 1.0/rinv
       lcolhc = idevh.eq.4
       write(*,5500) rob, plsize, cpfac, enfac, wsfac, 
     &               dazim, delev,
     &               lcolhc,
     &               label_surf, 
     &               label_strp, 
     &               label_vrtx, 
     &               label_body
 5500  format(/' ================================================'
     &        /'   D ist. to observer      currently ',f12.2
     &        /'   S ize of plot           currently ',f10.2
     &        /'   L oading  scale factor  currently ',f10.4
     &        /'   N ormal   scale factor  currently ',f10.4
     &        /'   V elocity scale factor  currently ',f10.4
     &        /'   R otation angle steps   currently ',2f10.4
     &        /'   C olor hardcopy         currently ',l2
     &        /'   LN label surfaces       currently ',l2
     &        /'   LS label strips         currently ',l2
     &        /'   LV label vortices       currently ',l2
     &        /'   LB label bodies         currently ',l2
     &       //' Select item to change (or <return>):  ', $)
       read(*,1,err=50) line
       call bstrip(line,nlin)
       opt = line(1:2)
       call touper(opt)

       line(1:2) = '  '
       call bstrip(line,nlin)
       narg = 2
       call getflt(line,arg,narg,error)

c----------------------------------------------------
       if(opt.eq.'  ') then
         if(linitview) then
           go to 4
          else
           go to 6
         endif

c----------------------------------------------------
       else if (opt.eq.'S ') then 
c------ change size
        if(narg.gt.0) then
         plsize = max( arg(1) , 0.0 )
        else
   52    write(*,1052)
 1052    format(' Enter new plot size:  ',$)
         read(*,*,err=52,end=50) plsize
        endif

c----------------------------------------------------
       else if(opt.eq.'D ') then
c------ set distance to observer (controls perspective)
        if(narg.gt.0) then
         rob = max( arg(1) , 0.0 )
        else
   53    write(*,1053)
 1053    format(' Enter new distance to observer (0 for infinity):  ',$)
         read (*,*,err=53,end=50) rob
        endif

        if(rob.eq.0.0) then
         rinv = 0.
        else
         rinv = 1.0/rob
        endif
        linitview = .true.

c----------------------------------------------------
       else if(opt.eq.'L ') then
c------ loading scale factor for dcp values (as fraction of cref)
        if(narg.gt.0) then
         cpfac = arg(1)
        else
   54    write(*,1054)
 1054    format(' Enter new loading scale factor:  ',$)
         read (*,*,err=54,end=50) cpfac
        endif

c----------------------------------------------------
       else if(opt.eq.'N ') then
c------ normal vector scale factor (as fraction of cref)
        if(narg.gt.0) then
         enfac = arg(1)
        else
 55      write(*,1055)
 1055    format(' Enter new normal scale factor:  ',$)
         read (*,*,err=55,end=50) enfac
        endif

c----------------------------------------------------
       else if(opt.eq.'V ') then
c------ velocity scale factor
        if(narg.gt.0) then
         wsfac = arg(1)
        else
 56      write(*,1056)
 1056    format(' Enter new velocity scale factor:  ',$)
         read (*,*,err=56,end=50) wsfac
        endif

c----------------------------------------------------
       else if(opt.eq.'R ') then
c------ rotation angle steps
        if    (narg.ge.2) then
         dazim = arg(1)
         delev = arg(2)
        elseif(narg.ge.1) then
         dazim = arg(1)
         delev = arg(1)
        else
 57      write(*,1057)
 1057    format(' Enter new rotation angles (deg):  ',$)
         read (*,*,err=57,end=50) dazim, delev
        endif

c----------------------------------------------------
       else if(opt.eq.'C ') then
c------ toggle color hardcopy
        lcolhc = .not.lcolhc
        if(lcolhc) then
         idevh = 4
        else
         idevh = 2
        endif

c----------------------------------------------------
       else if(opt.eq.'LN') then
c------ labeling for surfaces
        label_surf = .not.label_surf
c----------------------------------------------------
       else if(opt.eq.'LS') then
c------ labeling for strips
        label_strp = .not.label_strp
c----------------------------------------------------
       else if(opt.eq.'LV') then
c------ labeling for vortex elements
        label_vrtx = .not.label_vrtx
c----------------------------------------------------
       else if(opt.eq.'LB') then
c------ labeling for bodies
        label_body = .not.label_body
       endif
       go to 50

      else
       write(*,*) ' * Unrecognized command'

      endif

      go to 6

      end ! plotvl


      subroutine selsurf(azim, elev, tilt, rinv)
c--- selects surfaces for plotting

      include 'jvl.inc'
      include 'jvlplt.inc'
      character*40 opt, opt2
      character*1 ckey
      logical error

    1 format(a)

c---- select surfaces to plot
   60 continue
      call viewinit(azim, elev, tilt, rinv)
      call hidinit(.true.)
      call glims(gmin,gmax,.false.)
      call axlims
cc    call offini
      call pvlini(title,azim,elev,tilt,version,.false.)
      call plotgeom
      call plflush

      write(*,*) ' '
      write(*,*) '==================================================='
      write(*,*) '#surf  comp.  Nchord Nspan   plot?        name'
      do n = 1, nsurf
        write(*,61) n, lscomp(n), nk(n), nj(n), lpltsurf(n), stitle(n)
      end do
   61 format(1x,i3,4x,i3, 4x,i3,4x,i3,4x,l3,6x,a)


 62   write(*,*) '---------------------------------------------------'
      write(*,1062) imarksurf
 1062 format(/'   #    Select surface #  or  #:#   selects range',
     &       /'  -#  DeSelect surface #  or -#:# deselects range,'
     &       /'   A    Select all surfaces',
     &       /'   N  DeSelect all surfaces',
     &       /'   M    Mark surface, currently', i3 )

 65   write(*,1065)
 1065 format(/' Enter selection: (P prints surface list):  ',$)

      read(*,1,err=65) opt
cc      call lc2uc(opt)
      kc = index(opt,':')
      km = index(opt,',')

      ckey = opt(1:1)

      if(ckey.eq.' ') then
        go to 100

c--- read colon-separated number specs - i.e. 1:3 or -5:8
       elseif(kc.gt.1) then
        read(opt(1:kc-1),*,err=62) n1
        read(opt(kc+1:len(opt)),*,err=62) n2
        l1 = abs(n1)
        l2 = abs(n2)
        l1 = min(l1,nsurf)
        l2 = min(l2,nsurf)
        if(l1.gt.0 .and. l2.ge.l1) then
          do n = l1, l2
            lpltsurf(n) = (n1.gt.0)
          end do
        endif

c--- read comma-separated number string - i.e. 1,2,3,-5,-8
       elseif(km.gt.1) then
        opt2 = opt
 67     read(opt2(1:km-1),*,err=62) n
        nabs = abs(n)
        if(nabs.le.0 .or. nabs.gt.nsurf) go to 65
        lpltsurf(nabs) = (n.gt.0) 
        opt2 = opt2(km+1:len(opt2))
        km = index(opt2,',')
        if(km.gt.1) go to 67

        read(opt2,*,err=62) n
        nabs = abs(n)
        if(nabs.le.0 .or. nabs.gt.nsurf) go to 65
        lpltsurf(nabs) = (n.gt.0) 

      elseif(index('Pp',ckey).ne.0) then
        go to 60

      elseif(index('Aa',ckey).ne.0) then
        do n=1,nsurf
          lpltsurf(n) = .true.
        end do

      elseif(index('Nn',ckey).ne.0) then
        do n=1,nsurf
          lpltsurf(n) = .false.
        end do

      elseif(index('Mm',ckey).ne.0) then
        ninp = 1
        call getint(opt(2:40),imarksurf,ninp,error)

        if(ninp.eq.0 .or. error) then
         write(*,1068)
 1068    format(/' Enter surface # to mark:  ', $)
         read (*,*,err=65,end=65) imarksurf
        endif

        if(imarksurf.lt.0)     imarksurf = 0
        if(imarksurf.gt.nsurf) imarksurf = 0

      else
        read(opt,*,err=65) n
        nabs = abs(n)
        if(nabs.le.0 .or. nabs.gt.nsurf) go to 65
        lpltsurf(nabs) = (n.gt.0) 
      endif
      go to 60
cc
 100  return
      end ! selsurf




      subroutine pvlini(title,azim,elev,tilt,version,ldevm)
c...purpose:    initializes plot frame for plotting and 
c               puts up axis orientation and labels

      include 'jvlplt.inc'

      character*80 title
      logical ldevm

      character*1 chxyz(3)
      data chxyz / 'x', 'y', 'z' /

      if (lplot) call plend

      scs = 0.7*plsize*csize
      dlen = 0.045*plsize

      if(ldevm) then
       call pltini(idevm)
      else
       call pltini(idev)
      endif

c
c      if(.not.lplot) call colormap(ncolors)
      lplot = .true. 
      xplt = pmarg
      yplt = pmarg
      call newpen(2)
      call plchar(xplt,yplt,0.6*scs,'JVL ',0.0,4)
      call plnumb(999.,999.,0.6*scs,version,0.0,2)

      call getcolor(icol0)

c---- small axis arrows
      call newpen(3)
c     call newcolorname('green')
      call newcolorrgb(0,190,0)
      xo = pmarg + dlen + 1.5*scs
      yo = pmarg + dlen + 8.5*scs

      afac = dlen/arwlen
      do ixyz = 1, 3
        call arwplt(xo,yo,afac,
     &              uaxarwp(1,1,1,ixyz),
     &              uaxdirp(1,ixyz),rheadp,nhead)
      enddo

c---- xyz labels
      ltip = nhead+1
      do ixyz = 1, 3
        xplt = xo + uaxarwp(1,2,ltip,ixyz)*afac
     &            + uaxdirp(1       ,ixyz)*1.3*scs - 0.5*scs
        yplt = yo + uaxarwp(2,2,ltip,ixyz)*afac
     &            + uaxdirp(2       ,ixyz)*1.3*scs - 0.5*scs
        call plchar(xplt,yplt,0.9*scs,chxyz(ixyz),0.0,1)
      enddo
      call newcolor(icol0)

c---- azimuth and elevation angles
      call newpen(3)
      xplt = pmarg
      yplt = pmarg + 3.0*scs

      call plchar(xplt             ,yplt,0.75*scs,'elev',0.0, 4)
      call plchar(xplt+4.5*0.75*scs,yplt,0.75*scs,'='   ,0.0, 1)
      call plnumb(xplt+6.0*0.75*scs,yplt,0.75*scs,elev  ,0.0,-1)
      call plmath(999.             ,999.,0.75*scs,'"'   ,0.0, 1)

      yplt = yplt + 1.8*scs
      call plchar(xplt             ,yplt,0.75*scs,'azim',0.0, 4)
      call plchar(xplt+4.5*0.75*scs,yplt,0.75*scs,'='   ,0.0, 1)
      call plnumb(xplt+6.0*0.75*scs,yplt,0.75*scs,azim  ,0.0,-1)
      call plmath(999.             ,999.,0.75*scs,'"'   ,0.0, 1)

c---- configuration title
      call newpen(4)
      xplt = pmarg + 10.0*scs
      yplt = pmarg
      call plchar(xplt,yplt,1.1*scs,title,0.0,80)

c---- move plot origin to miss axis block
      call plot(pmarg+2.0*dlen+2.0*scs,pmarg+5.0*scs,-3)   
      call newfactor(plsize)
      call newpen(2)

c      call plot(0.,0.,3)   
c      call plot(1.,0.,2)   
c      call plot(1.,ar,2)   
c      call plot(0.,ar,2)   
c      call plot(0.,0.,2)   

      return        
      end ! pvlini  



      subroutine plotgeom
c--------------------------------------------------------------------------
c     Plots selected geometry data:
c      surfaces     (panels of chordwise strips of vortices)
c      bodies       (rings along centerline)
c      strips       (chordwise strips of vortices)
c      vortex legs  (bound vortex legs)
c      control pts  (vortex element control points)
c      normal vector (normal vector at control points)
c      camber slope (camber of each strip - used for establishing bc's)
c      hinge lines  (surface deflection axis and deflected surface outline)
c      strip loading (chordwise plot of vortex loading on each strip)
c      actuator disks

c     User can optionally display an identifying index for each item 
c      (surface, strip, element)
c--------------------------------------------------------------------------
      include 'jvl.inc'
      include 'jvlplt.inc'

      logical lvis

      real entot(3)
      real rrot(3), vrot(3), vctot(3), vvtot(3)
      real plab(3), plab2(3), svec1(3), svec2(3)
      real ptsorg(3)
      real pts_lines(3,2,nvmax), id_lines(nvmax),
     &     pts_lproj(3,2,nvmax)
      real pts_point(3,nlmax),
     &     pts_pproj(3,nlmax),
     &     rnum_point(nlmax)
      logical loff(nlmax)
      integer idum(1)
      real dum(1,1)

      real rcp0(3,nvmax), 
     &     rcps(3,nvmax), 
     &     rcpj(3,nvmax)

      real en(3)

      data id_lines / nvmax * 0.0 /

c---- number of line segments for drawing a body section "hoop"
      data nhoop / 18 /

      include 'masks.inc'

      eps    = 1.0e-5
      chsize = 0.01

      if(laxesplt) then
c----- plot x,y,z axes
       dtick = 0.08 * min( axdel(1), axdel(2), axdel(3) )
       chn = 1.0*dtick*sf

       ip = 0
       in = 0

c----- go over x,y,z axes
       do 50 iax = 1, 3
c------- go over axis tick marks
         do iann = 1, naxann(iax)
c--------- previous and current tick mark location (x, y, or z value)
           xyzm = axmin(iax) + axdel(iax)*float(iann-2)
           xyzo = axmin(iax) + axdel(iax)*float(iann-1)

c--------- go over x,y,z directions
           do 505 idir = 1, 3
             if(iann.eq.1 .and. idir.eq.iax) go to 505

             if(idir.eq.iax) then
c------------ axis segment, along iax direction
              ip = ip+1
              do k = 1, 3
                pts_lines(k,1,ip) = 0.
                pts_lines(k,2,ip) = 0.
              enddo
              pts_lines(iax,1,ip) = xyzm
              pts_lines(iax,2,ip) = xyzo

             else
c------------ tick marks, along the other two directions
              ip = ip+1
              do k = 1, 3
                pts_lines(k,1,ip) = 0.
                pts_lines(k,2,ip) = 0.
              enddo
              pts_lines(iax ,1,ip) = xyzo
              pts_lines(iax ,2,ip) = xyzo
              pts_lines(idir,1,ip) = -0.5*dtick
              pts_lines(idir,2,ip) =  0.5*dtick
             endif
 505       continue

c--------- save point on axis for annotation number placement
           if(xyzo .ne. 0.0) then
            in = in+1
            rnum_point(in) = xyzo
            if    (iax.eq.1) then
             pts_point(1,in) = xyzo
             pts_point(2,in) = 0.
             pts_point(3,in) = -0.5*dtick - 1.3*chn/sf
             loff(in) = .false.
            elseif(iax.eq.2) then
             pts_point(1,in) = 0.
             pts_point(2,in) = xyzo
             pts_point(3,in) = -0.5*dtick - 1.3*chn/sf
             loff(in) = .false.
            else
             pts_point(1,in) = 0.
             pts_point(2,in) = 0.
             pts_point(3,in) = xyzo
             loff(in) = .true.
            endif
ccc            id_lines(in) = iax
            id_lines(in) = 0
           endif

         enddo
 50    continue

c-----------------------------------
c----- plot axes
       call getpen(ipen0)
       call getcolor(icol0)

       call newpen(1)
       call newcolorrgb(0,128,0)

       nlines = ip
       nproj = 2*nlines
       call viewproj(pts_lines,nproj,pts_lproj)
       call plotlines(nlines,pts_lproj,id_lines)

       npoint = in
       call viewproj(pts_point,npoint,pts_pproj)

       do in = 1, npoint
         lvis = .true.
         if(lhidden) then
          id = 0
          ngrp = 0
          idum(1) = 0
          dum(1,1) = 0.
          call hidpnt(pts_pproj(1,in),id,ngrp,idum,dum,ntri,tri,lvis)
         endif

         rnum = rnum_point(in)
         if(lvis) call pltflt(pts_pproj(1,in),rnum,chn,loff(in),-2)
       enddo

c----- plot moment-reference point
       call newpen(3)
       call newcolorname('red')
       ip = 0
       do idir = 1, 3
         ip = ip+1
         pts_lines(1,1,ip) = xyzref(1)
         pts_lines(2,1,ip) = xyzref(2)
         pts_lines(3,1,ip) = xyzref(3)
         pts_lines(idir,1,ip) = pts_lines(idir,1,ip) - 0.7*dtick
         id_lines(ip) = 0

         pts_lines(1,2,ip) = xyzref(1)
         pts_lines(2,2,ip) = xyzref(2)
         pts_lines(3,2,ip) = xyzref(3)
         pts_lines(idir,2,ip) = pts_lines(idir,2,ip) + 0.7*dtick
         id_lines(ip) = 0
       enddo

       nlines = ip
       nproj = 2*nlines
       call viewproj(pts_lines,nproj,pts_lproj)
       call plotlines(nlines,pts_lproj,id_lines)

       call newpen(ipen0)
       call newcolor(icol0)
      endif


c---- go over surfaces
      do 100 n = 1, nsurf
c------ skip this surface if it is not to be plotted
        if(.not.lpltsurf(n)) go to 100

c------ emphasize panel 'imarksurf' with color
        if(imarksurf.eq.n) then
         call newcolorname('red')
        else
         call newcolor(icol0)
        endif

        if(idups(n).ge.0) then
         j1 = jfrst(n)
         jn = jlast(n)
         jinc = 1
        else
         j1 = jlast(n)
         jn = jfrst(n)
         jinc = -1
        endif

c------------------------------------------------------------------
c------ bound vortex legs plotted
        if(lboundplt) then
         ip = 0
c         write(*,*) 'surf j1 jn', n, j1, jn, jinc
         do j = j1, jn, jinc
           do iv = ifrsts(j), ilasts(j)
             ip = ip+1
             xave = 0.5*(rv1(1,iv) + rv2(1,iv))
             yave = 0.5*(rv1(2,iv) + rv2(2,iv))
             zave = 0.5*(rv1(3,iv) + rv2(3,iv))

c--- optionally regress vortex edges slightly for visibility
             regr = 0.0
             pts_lines(1,1,ip) = (1.0-regr)*(rv1(1,iv) + regr*xave)
             pts_lines(2,1,ip) = (1.0-regr)*(rv1(2,iv) + regr*yave)
             pts_lines(3,1,ip) = (1.0-regr)*(rv1(3,iv) + regr*zave)
             pts_lines(1,2,ip) = (1.0-regr)*(rv2(1,iv) + regr*xave)
             pts_lines(2,2,ip) = (1.0-regr)*(rv2(2,iv) + regr*yave)
             pts_lines(3,2,ip) = (1.0-regr)*(rv2(3,iv) + regr*zave)
             id_lines(ip) = j
           end do

           do iv = ifrstu(j), ilastu(j)
             ip = ip+1
             xave = 0.5*(rv1(1,iv) + rv2(1,iv))
             yave = 0.5*(rv1(2,iv) + rv2(2,iv))
             zave = 0.5*(rv1(3,iv) + rv2(3,iv))

c--- optionally regress vortex edges slightly for visibility
             regr = 0.0
             pts_lines(1,1,ip) = (1.0-regr)*(rv1(1,iv) + regr*xave)
             pts_lines(2,1,ip) = (1.0-regr)*(rv1(2,iv) + regr*yave)
             pts_lines(3,1,ip) = (1.0-regr)*(rv1(3,iv) + regr*zave)
             pts_lines(1,2,ip) = (1.0-regr)*(rv2(1,iv) + regr*xave)
             pts_lines(2,2,ip) = (1.0-regr)*(rv2(2,iv) + regr*yave)
             pts_lines(3,2,ip) = (1.0-regr)*(rv2(3,iv) + regr*zave)
             id_lines(ip) = j
           end do

           do iv = ifrstw(j), ilastw(j)
             ip = ip+1
             xave = 0.5*(rv1(1,iv) + rv2(1,iv))
             yave = 0.5*(rv1(2,iv) + rv2(2,iv))
             zave = 0.5*(rv1(3,iv) + rv2(3,iv))

c--- optionally regress vortex edges slightly for visibility
             regr = 0.0
             pts_lines(1,1,ip) = (1.0-regr)*(rv1(1,iv) + regr*xave)
             pts_lines(2,1,ip) = (1.0-regr)*(rv1(2,iv) + regr*yave)
             pts_lines(3,1,ip) = (1.0-regr)*(rv1(3,iv) + regr*zave)
             pts_lines(1,2,ip) = (1.0-regr)*(rv2(1,iv) + regr*xave)
             pts_lines(2,2,ip) = (1.0-regr)*(rv2(2,iv) + regr*yave)
             pts_lines(3,2,ip) = (1.0-regr)*(rv2(3,iv) + regr*zave)
             id_lines(ip) = j
           end do
         end do

         nlines = ip
         nproj = 2*nlines
         call viewproj(pts_lines,nproj,pts_lproj)
c 
         call getpen(ipen)
         call getcolor(icol)

c         if(n.le.nsurf) then
          call newpen(4)
          call newcolorname('magenta')
c         else
c          call newpen(2)
c          call newcolorrgb(128,0,128)
c         endif
         call plotlines(nlines,pts_lproj,id_lines)

         call newpen(ipen)
         call newcolor(icol)
        endif

c------------------------------------------------------------------
c------ control points plotted
        if(lcntlplt) then
         scl = 0.01*min(cref,0.5*bref)

         ip = 0
         do j = j1, jn, jinc
c--------- define vector on strip plane for control point "x" mark
           svec1(1) =  1.0    *scl
           svec1(2) = -ensz(j)*scl
           svec1(3) =  ensy(j)*scl
           svec2(1) = -svec1(1)
           svec2(2) =  svec1(2)
           svec2(3) =  svec1(3)
           do iv = ifrsts(j), ilasts(j)
             ip = ip+1
             pts_lines(1,1,ip) = rc(1,iv) - svec1(1)  
             pts_lines(2,1,ip) = rc(2,iv) - svec1(2)  
             pts_lines(3,1,ip) = rc(3,iv) - svec1(3)  
             pts_lines(1,2,ip) = rc(1,iv) + svec1(1)  
             pts_lines(2,2,ip) = rc(2,iv) + svec1(2)  
             pts_lines(3,2,ip) = rc(3,iv) + svec1(3)  
             id_lines(ip) = j
             ip = ip+1
             pts_lines(1,1,ip) = rc(1,iv) - svec2(1)  
             pts_lines(2,1,ip) = rc(2,iv) - svec2(2)  
             pts_lines(3,1,ip) = rc(3,iv) - svec2(3)  
             pts_lines(1,2,ip) = rc(1,iv) + svec2(1)  
             pts_lines(2,2,ip) = rc(2,iv) + svec2(2)  
             pts_lines(3,2,ip) = rc(3,iv) + svec2(3)  
             id_lines(ip) = j
           end do
           do iv = ifrstu(j), ilastu(j)
             ip = ip+1
             pts_lines(1,1,ip) = rc(1,iv) - svec1(1)  
             pts_lines(2,1,ip) = rc(2,iv) - svec1(2)  
             pts_lines(3,1,ip) = rc(3,iv) - svec1(3)  
             pts_lines(1,2,ip) = rc(1,iv) + svec1(1)  
             pts_lines(2,2,ip) = rc(2,iv) + svec1(2)  
             pts_lines(3,2,ip) = rc(3,iv) + svec1(3)  
             id_lines(ip) = j
             ip = ip+1
             pts_lines(1,1,ip) = rc(1,iv) - svec2(1)  
             pts_lines(2,1,ip) = rc(2,iv) - svec2(2)  
             pts_lines(3,1,ip) = rc(3,iv) - svec2(3)  
             pts_lines(1,2,ip) = rc(1,iv) + svec2(1)  
             pts_lines(2,2,ip) = rc(2,iv) + svec2(2)  
             pts_lines(3,2,ip) = rc(3,iv) + svec2(3)  
             id_lines(ip) = j
           end do
           do iv = ifrstw(j), ilastw(j)
             ip = ip+1
             pts_lines(1,1,ip) = rc(1,iv) - svec1(1)  
             pts_lines(2,1,ip) = rc(2,iv) - svec1(2)  
             pts_lines(3,1,ip) = rc(3,iv) - svec1(3)  
             pts_lines(1,2,ip) = rc(1,iv) + svec1(1)  
             pts_lines(2,2,ip) = rc(2,iv) + svec1(2)  
             pts_lines(3,2,ip) = rc(3,iv) + svec1(3)  
             id_lines(ip) = j
             ip = ip+1
             pts_lines(1,1,ip) = rc(1,iv) - svec2(1)  
             pts_lines(2,1,ip) = rc(2,iv) - svec2(2)  
             pts_lines(3,1,ip) = rc(3,iv) - svec2(3)  
             pts_lines(1,2,ip) = rc(1,iv) + svec2(1)  
             pts_lines(2,2,ip) = rc(2,iv) + svec2(2)  
             pts_lines(3,2,ip) = rc(3,iv) + svec2(3)  
             id_lines(ip) = j
           end do
         end do

         nlines = ip
         nproj = 2*nlines
         call viewproj(pts_lines,nproj,pts_lproj)
c 
         call getpen(ipen)
         call getcolor(icol)

         call newpen(4)
         call newcolorname('yellow')
         call plotlines(nlines,pts_lproj,id_lines)

         call newpen(ipen)
         call newcolor(icol)


         ip = 0
         do j = j1, jn, jinc
c--------- define vector on strip plane for bound vortex midpoint "x" mark
           svec1(1) =  1.0    *scl
           svec1(2) = -ensz(j)*scl
           svec1(3) =  ensy(j)*scl
           svec2(1) = -svec1(1)
           svec2(2) =  svec1(2)
           svec2(3) =  svec1(3)
           do iv = ifrsts(j), ilasts(j)
             ip = ip+1
             pts_lines(1,1,ip) = rv(1,iv) - svec1(1)  
             pts_lines(2,1,ip) = rv(2,iv) - svec1(2)  
             pts_lines(3,1,ip) = rv(3,iv) - svec1(3)  
             pts_lines(1,2,ip) = rv(1,iv) + svec1(1)  
             pts_lines(2,2,ip) = rv(2,iv) + svec1(2)  
             pts_lines(3,2,ip) = rv(3,iv) + svec1(3)  
             id_lines(ip) = j
             ip = ip+1
             pts_lines(1,1,ip) = rv(1,iv) - svec2(1)  
             pts_lines(2,1,ip) = rv(2,iv) - svec2(2)  
             pts_lines(3,1,ip) = rv(3,iv) - svec2(3)  
             pts_lines(1,2,ip) = rv(1,iv) + svec2(1)  
             pts_lines(2,2,ip) = rv(2,iv) + svec2(2)  
             pts_lines(3,2,ip) = rv(3,iv) + svec2(3)  
             id_lines(ip) = j
           end do

           do iv = ifrstu(j), ilastu(j)
             ip = ip+1
             pts_lines(1,1,ip) = rv(1,iv) - svec1(1)  
             pts_lines(2,1,ip) = rv(2,iv) - svec1(2)  
             pts_lines(3,1,ip) = rv(3,iv) - svec1(3)  
             pts_lines(1,2,ip) = rv(1,iv) + svec1(1)  
             pts_lines(2,2,ip) = rv(2,iv) + svec1(2)  
             pts_lines(3,2,ip) = rv(3,iv) + svec1(3)  
             id_lines(ip) = j
             ip = ip+1
             pts_lines(1,1,ip) = rv(1,iv) - svec2(1)  
             pts_lines(2,1,ip) = rv(2,iv) - svec2(2)  
             pts_lines(3,1,ip) = rv(3,iv) - svec2(3)  
             pts_lines(1,2,ip) = rv(1,iv) + svec2(1)  
             pts_lines(2,2,ip) = rv(2,iv) + svec2(2)  
             pts_lines(3,2,ip) = rv(3,iv) + svec2(3)  
             id_lines(ip) = j
           end do

           do iv = ifrstw(j), ilastw(j)
             ip = ip+1
             pts_lines(1,1,ip) = rv(1,iv) - svec1(1)  
             pts_lines(2,1,ip) = rv(2,iv) - svec1(2)  
             pts_lines(3,1,ip) = rv(3,iv) - svec1(3)  
             pts_lines(1,2,ip) = rv(1,iv) + svec1(1)  
             pts_lines(2,2,ip) = rv(2,iv) + svec1(2)  
             pts_lines(3,2,ip) = rv(3,iv) + svec1(3)  
             id_lines(ip) = j
             ip = ip+1
             pts_lines(1,1,ip) = rv(1,iv) - svec2(1)  
             pts_lines(2,1,ip) = rv(2,iv) - svec2(2)  
             pts_lines(3,1,ip) = rv(3,iv) - svec2(3)  
             pts_lines(1,2,ip) = rv(1,iv) + svec2(1)  
             pts_lines(2,2,ip) = rv(2,iv) + svec2(2)  
             pts_lines(3,2,ip) = rv(3,iv) + svec2(3)  
             id_lines(ip) = j
           end do
         end do

         nlines = ip
         nproj = 2*nlines
         call viewproj(pts_lines,nproj,pts_lproj)

         call getpen(ipen)
         call getcolor(icol)

         call newcolorname('orange')
         call plotlines(nlines,pts_lproj,id_lines)

         call newpen(ipen)
         call newcolor(icol)
        endif


c------------------------------------------------------------------
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

c--------- plot trailing legs back 2.0*bref if enabled
           if(lwakeplt .and. lfwake(n)) then
             ip = ip+1
             pts_lines(1,1,ip) = pts_lines(1,2,ip-1)
             pts_lines(2,1,ip) = pts_lines(2,2,ip-1)
             pts_lines(3,1,ip) = pts_lines(3,2,ip-1)
             pts_lines(1,2,ip) = pts_lines(1,1,ip) + 2.0*bref
             pts_lines(2,2,ip) = pts_lines(2,1,ip)
             pts_lines(3,2,ip) = pts_lines(3,1,ip)
             id_lines(ip) = 0
           endif
         end do
         ip = ip+1
         pts_lines(1,1,ip) = rle2(1,jn)
         pts_lines(2,1,ip) = rle2(2,jn)
         pts_lines(3,1,ip) = rle2(3,jn)
         pts_lines(1,2,ip) = rle2(1,jn) + chord2(jn)
         pts_lines(2,2,ip) = rle2(2,jn)
         pts_lines(3,2,ip) = rle2(3,jn)
         id_lines(ip) = jn

c------- plot trailing legs back 2.0*bref if enabled
         if(lwakeplt .and. lfwake(n)) then
           ip = ip+1
           pts_lines(1,1,ip) = pts_lines(1,2,ip-1)
           pts_lines(2,1,ip) = pts_lines(2,2,ip-1)
           pts_lines(3,1,ip) = pts_lines(3,2,ip-1)
           pts_lines(1,2,ip) = pts_lines(1,1,ip) + 2.0*bref
           pts_lines(2,2,ip) = pts_lines(2,1,ip)
           pts_lines(3,2,ip) = pts_lines(3,1,ip)
           id_lines(ip) = 0
         endif

         nlines = ip
         nproj = 2*nlines
         call viewproj(pts_lines,nproj,pts_lproj)

c------- plot chordlines dashed 
         call getpen(ipen)

         call newpen(2)
         call newpat(lmask2)
         call plotlines(nlines,pts_lproj,id_lines)
         call newpat(lmask0)
         call newpen(ipen)
        endif

c------------------------------------------------------------------
c------ actuator disks plotted
        if(ldiskplt) then
c------- vertical chord lines
         ip = 0
         do j = j1, jn, jinc
           if(ldstrp(j)) then
            dy = rle2(2,j) - rle1(2,j)
            dz = rle2(3,j) - rle1(3,j)
            ds = sqrt(dy**2 + dz**2)
            esy =  dy/ds
            esz =  dz/ds
            eny = -dz/ds
            enz =  dy/ds
            dhy = 0.5*hdstrp(j) * eny
            dhz = 0.5*hdstrp(j) * enz

            ip = ip+1
            pts_lines(1,1,ip) = rle1(1,j) + dxdstrp1(j)
            pts_lines(2,1,ip) = rle1(2,j) + dndstrp1(j)*eny + dhy
            pts_lines(3,1,ip) = rle1(3,j) + dndstrp1(j)*enz + dhz
            pts_lines(1,2,ip) = rle1(1,j) + dxdstrp1(j)
            pts_lines(2,2,ip) = rle1(2,j) + dndstrp1(j)*eny - dhy
            pts_lines(3,2,ip) = rle1(3,j) + dndstrp1(j)*enz - dhz
            id_lines(ip) = j

            ip = ip+1
            pts_lines(1,1,ip) = rle2(1,j) + dxdstrp2(j)
            pts_lines(2,1,ip) = rle2(2,j) + dndstrp2(j)*eny + dhy
            pts_lines(3,1,ip) = rle2(3,j) + dndstrp2(j)*enz + dhz
            pts_lines(1,2,ip) = rle2(1,j) + dxdstrp2(j)
            pts_lines(2,2,ip) = rle2(2,j) + dndstrp2(j)*eny - dhy
            pts_lines(3,2,ip) = rle2(3,j) + dndstrp2(j)*enz - dhz
            id_lines(ip) = j

            ip = ip+1
            pts_lines(1,1,ip) = rle1(1,j) + dxdstrp1(j)
            pts_lines(2,1,ip) = rle1(2,j) + dndstrp1(j)*eny + dhy
            pts_lines(3,1,ip) = rle1(3,j) + dndstrp1(j)*enz + dhz
            pts_lines(1,2,ip) = rle2(1,j) + dxdstrp2(j)
            pts_lines(2,2,ip) = rle2(2,j) + dndstrp2(j)*eny + dhy
            pts_lines(3,2,ip) = rle2(3,j) + dndstrp2(j)*enz + dhz
            id_lines(ip) = j

            ip = ip+1
            pts_lines(1,1,ip) = rle1(1,j) + dxdstrp1(j)
            pts_lines(2,1,ip) = rle1(2,j) + dndstrp1(j)*eny - dhy
            pts_lines(3,1,ip) = rle1(3,j) + dndstrp1(j)*enz - dhz
            pts_lines(1,2,ip) = rle2(1,j) + dxdstrp2(j)
            pts_lines(2,2,ip) = rle2(2,j) + dndstrp2(j)*eny - dhy
            pts_lines(3,2,ip) = rle2(3,j) + dndstrp2(j)*enz - dhz
            id_lines(ip) = j

            ip = ip+1
            pts_lines(1,1,ip) = rle1(1,j) + dxdstrp1(j)
            pts_lines(2,1,ip) = rle1(2,j) + dndstrp1(j)*eny
            pts_lines(3,1,ip) = rle1(3,j) + dndstrp1(j)*enz
            pts_lines(1,2,ip) = rle2(1,j) + dxdstrp2(j)
            pts_lines(2,2,ip) = rle2(2,j) + dndstrp2(j)*eny
            pts_lines(3,2,ip) = rle2(3,j) + dndstrp2(j)*enz
            id_lines(ip) = j

            ip = ip+1
            pts_lines(1,1,ip) = rle1(1,j) + dxdstrp1(j)
            pts_lines(2,1,ip) = rle1(2,j) + dndstrp1(j)*eny + dhy
            pts_lines(3,1,ip) = rle1(3,j) + dndstrp1(j)*enz + dhz
            pts_lines(1,2,ip) = rle2(1,j) + dxdstrp2(j)
            pts_lines(2,2,ip) = rle2(2,j) + dndstrp2(j)*eny + dhy
            pts_lines(3,2,ip) = rle2(3,j) + dndstrp2(j)*enz + dhz
            id_lines(ip) = j

           endif
         enddo

         nlines = ip
         if(nlines .ge. 0) then
          nproj = 2*nlines
          call viewproj(pts_lines,nproj,pts_lproj)

          call getpen(ipen)
          call newpen(2)
          call newpat(lmask1)
          call plotlines(nlines,pts_lproj,id_lines)
          call newpat(lmask0)
          call newpen(ipen)
         endif

        endif

c------------------------------------------------------------------
c------ strip camberline plotted (approx., integrated from camber slopes)
        if(lcambplt) then
         call getpen(ipen)
         call getcolor(icol)

         do ipass = 1, 2
           ip = 0
           do j = j1, jn, jinc
             xx0 = rle(1,j)
             yy0 = rle(2,j)
             zz0 = rle(3,j) + eps
             do iv = ifrsts(j), ilasts(j)
               dent = 0.
               if(ipass.eq.1) then
                do k = 1, ncontrol
                  dent = dent
     &               + 0.5*(enc_d(1,iv,k) + env_d(1,iv,k))*delcon(k)
                enddo
               endif
               entot(1) = 0.5*(enc(1,iv) + env(1,iv)) + dent

               tc = -dxv(iv)*entot(1)
               xx1 = xx0 + dxv(iv)
               yy1 = yy0 + tc*ensy(j)
               zz1 = zz0 + tc*ensz(j)

               ip = ip+1
               pts_lines(1,1,ip) = xx0
               pts_lines(2,1,ip) = yy0
               pts_lines(3,1,ip) = zz0
               pts_lines(1,2,ip) = xx1
               pts_lines(2,2,ip) = yy1
               pts_lines(3,2,ip) = zz1
               id_lines(ip) = 0

               xx0 = xx1
               yy0 = yy1
               zz0 = zz1
             enddo
           enddo

           nlines = ip
           nproj = 2*nlines
           call viewproj(pts_lines,nproj,pts_lproj)
c 
           if(ipass.eq.1) then
            call newpen(4)
            call newcolorname('orange')
           else
            call newpen(4)
            call newcolorname('cyan')
           endif
           call plotlines(nlines,pts_lproj,id_lines)
         enddo

         call newpen(ipen)
         call newcolor(icol)
        endif

c------------------------------------------------------------------
c------ strip loading plotted
        if(lloadplt) then

         cpscl = cpfac*cref
         do j = j1, jn, jinc
           do isuw = 1, 3
           if    (isuw .eq. 1) then
            ifrst = ifrsts(j)
            ilast = ilasts(j)
           elseif(isuw .eq. 2) then
            ifrst = ifrstu(j)
            ilast = ilastu(j)
           else
            ifrst = ifrstw(j)
            ilast = ilastw(j)
           endif

           do iv = ifrst, ilast
             en(1) = 0.
             en(2) = ensy(j)
             en(3) = ensz(j)
             do k = 1, 3
               rcp0(k,iv) = rv(k,iv)
               rcpj(k,iv) = rv(k,iv) + cpscl* dcpj(iv)          *en(k)
               rcps(k,iv) = rv(k,iv) + cpscl*(dcpj(iv)+dcpv(iv))*en(k)
             enddo
           enddo ! next iv
           enddo ! next isuw
         enddo ! next j

         call getpen(ipen)
         call getcolor(icol)

c------- plot surface load vectors
         call newpen(5)
         call newcolorname('green')
         ip = 0
         do j = j1, jn, jinc
           do iv = ifrsts(j), ilasts(j)
             ip = ip+1
             pts_lines(1,1,ip) = rcp0(1,iv)
             pts_lines(2,1,ip) = rcp0(2,iv)
             pts_lines(3,1,ip) = rcp0(3,iv)
             pts_lines(1,2,ip) = rcps(1,iv)
             pts_lines(2,2,ip) = rcps(2,iv)
             pts_lines(3,2,ip) = rcps(3,iv)
             id_lines(ip) = 0
           end do
         end do
         nlines = ip
         nproj = 2*nlines
         call viewproj(pts_lines,nproj,pts_lproj)
         call plotlines(nlines,pts_lproj,id_lines)

         if(lwaklplt) then
c-------- plot jet load vectors
          call newpen(5)
          call newcolorname('cyan')
          ip = 0
          do j = j1, jn, jinc
            do iv = ifrstu(j), ilastu(j)
              ip = ip+1
              pts_lines(1,1,ip) = rcp0(1,iv)
              pts_lines(2,1,ip) = rcp0(2,iv)
              pts_lines(3,1,ip) = rcp0(3,iv)
              pts_lines(1,2,ip) = rcpj(1,iv)
              pts_lines(2,2,ip) = rcpj(2,iv)
              pts_lines(3,2,ip) = rcpj(3,iv)
              id_lines(ip) = 0
            end do
            do iv = ifrstw(j), ilastw(j)
              ip = ip+1
              pts_lines(1,1,ip) = rcp0(1,iv)
              pts_lines(2,1,ip) = rcp0(2,iv)
              pts_lines(3,1,ip) = rcp0(3,iv)
              pts_lines(1,2,ip) = rcpj(1,iv)
              pts_lines(2,2,ip) = rcpj(2,iv)
              pts_lines(3,2,ip) = rcpj(3,iv)
              id_lines(ip) = 0
            end do
          end do
          nlines = ip
          nproj = 2*nlines
          call viewproj(pts_lines,nproj,pts_lproj)
          call plotlines(nlines,pts_lproj,id_lines)
         endif

c------- plot jet-curvature surface load vectors
         call newpen(5)
         call newcolorname('orange')
         ip = 0
         do j = j1, jn, jinc
           do iv = ifrsts(j), ilasts(j)
             ip = ip+1
             pts_lines(1,1,ip) = rcp0(1,iv)
             pts_lines(2,1,ip) = rcp0(2,iv)
             pts_lines(3,1,ip) = rcp0(3,iv)
             pts_lines(1,2,ip) = rcpj(1,iv)
             pts_lines(2,2,ip) = rcpj(2,iv)
             pts_lines(3,2,ip) = rcpj(3,iv)
             id_lines(ip) = 0
           end do
         end do
         nlines = ip
         nproj = 2*nlines
         call viewproj(pts_lines,nproj,pts_lproj)
         call plotlines(nlines,pts_lproj,id_lines)

c------- plot the loading lines 
         call newpen(4)
         call newcolorname('red')

         ip = 0
         do j = j1, jn, jinc
           do iv = ifrsts(j), ilasts(j)-1
             ip = ip+1
             pts_lines(1,1,ip) = rcps(1,iv)
             pts_lines(2,1,ip) = rcps(2,iv)
             pts_lines(3,1,ip) = rcps(3,iv)
             pts_lines(1,2,ip) = rcps(1,iv+1)
             pts_lines(2,2,ip) = rcps(2,iv+1)
             pts_lines(3,2,ip) = rcps(3,iv+1)
             id_lines(ip) = 0
           end do

           if(lwaklplt) then
            do iv = ifrstu(j), ilastu(j)-1
              ip = ip+1
              pts_lines(1,1,ip) = rcpj(1,iv)
              pts_lines(2,1,ip) = rcpj(2,iv)
              pts_lines(3,1,ip) = rcpj(3,iv)
              pts_lines(1,2,ip) = rcpj(1,iv+1)
              pts_lines(2,2,ip) = rcpj(2,iv+1)
              pts_lines(3,2,ip) = rcpj(3,iv+1)
              id_lines(ip) = 0
            enddo

            do iv = ifrstw(j), ilastw(j)-1
              ip = ip+1
              pts_lines(1,1,ip) = rcpj(1,iv)
              pts_lines(2,1,ip) = rcpj(2,iv)
              pts_lines(3,1,ip) = rcpj(3,iv)
              pts_lines(1,2,ip) = rcpj(1,iv+1)
              pts_lines(2,2,ip) = rcpj(2,iv+1)
              pts_lines(3,2,ip) = rcpj(3,iv+1)
              id_lines(ip) = 0
            enddo
           endif
         end do
         nlines = ip
         nproj = 2*nlines
         call viewproj(pts_lines,nproj,pts_lproj)
         call plotlines(nlines,pts_lproj,id_lines)

         call newpen(ipen)
         call newcolor(icol)
        endif

c------------------------------------------------------------------
c------ control-point total velocities
        if(lvelcplt) then
         call getcolor(icol0)
         call newcolorrgb(128,128,0)

         wref = 1.0
         wscl = wsfac/wref

         ip = 0
         do j = j1, jn, jinc
           do ic = ifrsts(j), ilasts(j)
             rrot(1) = rc(1,ic) - xyzref(1)
             rrot(2) = rc(2,ic) - xyzref(2)
             rrot(3) = rc(3,ic) - xyzref(3)
             call cross(rrot,wbar,vrot)

             do k = 1, 3
               vctot(k) = vbar(k)
     &                  + vrot(k)
     &                  + wc_u(k,ic,1)*vbar(1)
     &                  + wc_u(k,ic,2)*vbar(2)
     &                  + wc_u(k,ic,3)*vbar(3)
     &                  + wc_u(k,ic,4)*wbar(1)
     &                  + wc_u(k,ic,5)*wbar(2)
     &                  + wc_u(k,ic,6)*wbar(3)
     &                  + vc(k,ic)
             enddo

             ip = ip+1
             pts_lines(1,1,ip) = rc(1,ic)
             pts_lines(2,1,ip) = rc(2,ic)
             pts_lines(3,1,ip) = rc(3,ic)
             pts_lines(1,2,ip) = rc(1,ic) + wscl*vctot(1)
             pts_lines(2,2,ip) = rc(2,ic) + wscl*vctot(2)
             pts_lines(3,2,ip) = rc(3,ic) + wscl*vctot(3)
             id_lines(ip) = 0
           end do

           do ic = ifrstu(j), ilastu(j)
             rrot(1) = rc(1,ic) - xyzref(1)
             rrot(2) = rc(2,ic) - xyzref(2)
             rrot(3) = rc(3,ic) - xyzref(3)
             call cross(rrot,wbar,vrot)

             do k = 1, 3
               vctot(k) = vbar(k)
     &                  + vrot(k)
     &                  + wc_u(k,ic,1)*vbar(1)
     &                  + wc_u(k,ic,2)*vbar(2)
     &                  + wc_u(k,ic,3)*vbar(3)
     &                  + wc_u(k,ic,4)*wbar(1)
     &                  + wc_u(k,ic,5)*wbar(2)
     &                  + wc_u(k,ic,6)*wbar(3)
     &                  + vc(k,ic)
             enddo

             ip = ip+1
             pts_lines(1,1,ip) = rc(1,ic)
             pts_lines(2,1,ip) = rc(2,ic)
             pts_lines(3,1,ip) = rc(3,ic)
             pts_lines(1,2,ip) = rc(1,ic) + wscl*vctot(1)
             pts_lines(2,2,ip) = rc(2,ic) + wscl*vctot(2)
             pts_lines(3,2,ip) = rc(3,ic) + wscl*vctot(3)
             id_lines(ip) = 0
           end do

           do ic = ifrstw(j), ilastw(j)
             rrot(1) = rc(1,ic) - xyzref(1)
             rrot(2) = rc(2,ic) - xyzref(2)
             rrot(3) = rc(3,ic) - xyzref(3)
             call cross(rrot,wbar,vrot)

             do k = 1, 3
               vctot(k) = vbar(k)
     &                  + vrot(k)
     &                  + wc_u(k,ic,1)*vbar(1)
     &                  + wc_u(k,ic,2)*vbar(2)
     &                  + wc_u(k,ic,3)*vbar(3)
     &                  + wc_u(k,ic,4)*wbar(1)
     &                  + wc_u(k,ic,5)*wbar(2)
     &                  + wc_u(k,ic,6)*wbar(3)
     &                  + vc(k,ic)
             enddo

             ip = ip+1
             pts_lines(1,1,ip) = rc(1,ic)
             pts_lines(2,1,ip) = rc(2,ic)
             pts_lines(3,1,ip) = rc(3,ic)
             pts_lines(1,2,ip) = rc(1,ic) + wscl*vctot(1)
             pts_lines(2,2,ip) = rc(2,ic) + wscl*vctot(2)
             pts_lines(3,2,ip) = rc(3,ic) + wscl*vctot(3)
             id_lines(ip) = 0
           end do
         end do

         nlines = ip
         nproj = 2*nlines
         call viewproj(pts_lines,nproj,pts_lproj)
         call plotlines(nlines,pts_lproj,id_lines)
         call newcolor(icol0)

        endif

c------------------------------------------------------------------
c------ bound vortex total velocities
        if(lwsegplt) then
         call getpen(ipen)
         call getcolor(icol)

         call newpen(4)
         call newcolorname('yellow')

         wref = 1.0
         wscl = wsfac/wref

         ip = 0
         do j = j1, jn, jinc
           do iv = ifrsts(j), ilasts(j)
             rrot(1) = rv(1,iv) - xyzref(1)
             rrot(2) = rv(2,iv) - xyzref(2)
             rrot(3) = rv(3,iv) - xyzref(3)
             call cross(rrot,wbar,vrot)

             do k = 1, 3
               vvtot(k) = vbar(k)
     &                  + vrot(k)
     &                  + wv_u(k,iv,1)*vbar(1)
     &                  + wv_u(k,iv,2)*vbar(2)
     &                  + wv_u(k,iv,3)*vbar(3)
     &                  + wv_u(k,iv,4)*wbar(1)
     &                  + wv_u(k,iv,5)*wbar(2)
     &                  + wv_u(k,iv,6)*wbar(3)
     &                  + vv(k,iv)
             enddo

             ip = ip+1
             pts_lines(1,1,ip) = rv(1,iv)
             pts_lines(2,1,ip) = rv(2,iv)
             pts_lines(3,1,ip) = rv(3,iv)
             pts_lines(1,2,ip) = rv(1,iv) + wscl*vvtot(1)
             pts_lines(2,2,ip) = rv(2,iv) + wscl*vvtot(2)
             pts_lines(3,2,ip) = rv(3,iv) + wscl*vvtot(3)
             id_lines(ip) = 0
           end do

           do iv = ifrstu(j), ilastu(j)
             rrot(1) = rv(1,iv) - xyzref(1)
             rrot(2) = rv(2,iv) - xyzref(2)
             rrot(3) = rv(3,iv) - xyzref(3)
             call cross(rrot,wbar,vrot)

             do k = 1, 3
               vvtot(k) = vbar(k)
     &                  + vrot(k)
     &                  + wv_u(k,iv,1)*vbar(1)
     &                  + wv_u(k,iv,2)*vbar(2)
     &                  + wv_u(k,iv,3)*vbar(3)
     &                  + wv_u(k,iv,4)*wbar(1)
     &                  + wv_u(k,iv,5)*wbar(2)
     &                  + wv_u(k,iv,6)*wbar(3)
     &                  + vv(k,iv)
             enddo

             ip = ip+1
             pts_lines(1,1,ip) = rv(1,iv)
             pts_lines(2,1,ip) = rv(2,iv)
             pts_lines(3,1,ip) = rv(3,iv)
             pts_lines(1,2,ip) = rv(1,iv) + wscl*vvtot(1)
             pts_lines(2,2,ip) = rv(2,iv) + wscl*vvtot(2)
             pts_lines(3,2,ip) = rv(3,iv) + wscl*vvtot(3)
             id_lines(ip) = 0
           end do

           do iv = ifrstw(j), ilastw(j)
             rrot(1) = rv(1,iv) - xyzref(1)
             rrot(2) = rv(2,iv) - xyzref(2)
             rrot(3) = rv(3,iv) - xyzref(3)
             call cross(rrot,wbar,vrot)

             do k = 1, 3
               vvtot(k) = vbar(k)
     &                  + vrot(k)
     &                  + wv_u(k,iv,1)*vbar(1)
     &                  + wv_u(k,iv,2)*vbar(2)
     &                  + wv_u(k,iv,3)*vbar(3)
     &                  + wv_u(k,iv,4)*wbar(1)
     &                  + wv_u(k,iv,5)*wbar(2)
     &                  + wv_u(k,iv,6)*wbar(3)
     &                  + vv(k,iv)
             enddo

             ip = ip+1
             pts_lines(1,1,ip) = rv(1,iv)
             pts_lines(2,1,ip) = rv(2,iv)
             pts_lines(3,1,ip) = rv(3,iv)
             pts_lines(1,2,ip) = rv(1,iv) + wscl*vvtot(1)
             pts_lines(2,2,ip) = rv(2,iv) + wscl*vvtot(2)
             pts_lines(3,2,ip) = rv(3,iv) + wscl*vvtot(3)
             id_lines(ip) = 0
           end do
         end do

         nlines = ip
         nproj = 2*nlines
         call viewproj(pts_lines,nproj,pts_lproj)
         call plotlines(nlines,pts_lproj,id_lines)

         call newpen(ipen)
         call newcolor(icol)
        endif

c------------------------------------------------------------------

c------ surface outline
        ip = 0
        ip = ip+1
        pts_lines(1,1,ip) = rle1(1,j1)
        pts_lines(2,1,ip) = rle1(2,j1)
        pts_lines(3,1,ip) = rle1(3,j1)
        pts_lines(1,2,ip) = rle1(1,j1) + chord1(j1)
        pts_lines(2,2,ip) = rle1(2,j1)
        pts_lines(3,2,ip) = rle1(3,j1)
        id_lines(ip) = j1

c------ plot trailing legs back 2.0*bref if enabled
        if(lwakeplt .and. lfwake(n)) then
          ip = ip+1
          pts_lines(1,1,ip) = pts_lines(1,2,ip-1)
          pts_lines(2,1,ip) = pts_lines(2,2,ip-1)
          pts_lines(3,1,ip) = pts_lines(3,2,ip-1)
          pts_lines(1,2,ip) = pts_lines(1,1,ip) + 2.0*bref
          pts_lines(2,2,ip) = pts_lines(2,1,ip)
          pts_lines(3,2,ip) = pts_lines(3,1,ip)
          id_lines(ip) = 0
        endif

        ip = ip+1
        pts_lines(1,1,ip) = rle2(1,jn)
        pts_lines(2,1,ip) = rle2(2,jn)
        pts_lines(3,1,ip) = rle2(3,jn)
        pts_lines(1,2,ip) = rle2(1,jn) + chord2(jn)
        pts_lines(2,2,ip) = rle2(2,jn)
        pts_lines(3,2,ip) = rle2(3,jn)
        id_lines(ip) = jn

c------ plot trailing legs back 2.0*bref if enabled
        if(lwakeplt .and. lfwake(n)) then
          ip = ip+1
          pts_lines(1,1,ip) = pts_lines(1,2,ip-1)
          pts_lines(2,1,ip) = pts_lines(2,2,ip-1)
          pts_lines(3,1,ip) = pts_lines(3,2,ip-1)
          pts_lines(1,2,ip) = pts_lines(1,1,ip) + 2.0*bref
          pts_lines(2,2,ip) = pts_lines(2,1,ip)
          pts_lines(3,2,ip) = pts_lines(3,1,ip)
          id_lines(ip) = 0
        endif

        do j = j1, jn, jinc
          ip = ip+1
          pts_lines(1,1,ip) = rle1(1,j)
          pts_lines(2,1,ip) = rle1(2,j)
          pts_lines(3,1,ip) = rle1(3,j)
          pts_lines(1,2,ip) = rle2(1,j)
          pts_lines(2,2,ip) = rle2(2,j)
          pts_lines(3,2,ip) = rle2(3,j)
          id_lines(ip) = j
          ip = ip+1
          pts_lines(1,1,ip) = rle1(1,j) + chord1(j)
          pts_lines(2,1,ip) = rle1(2,j)
          pts_lines(3,1,ip) = rle1(3,j)
          pts_lines(1,2,ip) = rle2(1,j) + chord2(j)
          pts_lines(2,2,ip) = rle2(2,j)
          pts_lines(3,2,ip) = rle2(3,j)
          id_lines(ip) = j
        end do

        nlines = ip
        nproj = 2*nlines
        call viewproj(pts_lines,nproj,pts_lproj)
        call plotlines(nlines,pts_lproj,id_lines)

c---------------------------------------------------
c------ now do any labels (integer indices for surfaces,strips or vortices)

c------ surface labels
        if(label_surf) then
          jmid = jfrst(n) + ifix(0.5*float(nj(n))) - 1
          plab(1) = rle(1,jmid) + 0.6*chord(jmid)
          plab(2) = rle(2,jmid)
          plab(3) = rle(3,jmid)
          call viewproj(plab,1,pts_lproj)
          id = jmid
          lvis=.true.
          if(lhidden) then
            ngrp = 0
            idum(1) = 0
            dum(1,1) = 0.
            call hidpnt(pts_lproj,id,ngrp,idum,dum,ntri,tri,lvis)
          endif
          ilab = n
          if (lvis) call pltint(pts_lproj,ilab,1.5*chsize,.false.)
        endif

c------ strip labels
        if(label_strp) then
          ngrp = 0
          idum(1) = 0
          dum(1,1) = 0.
          call getpen(ipen)
          call getcolor(icol)

          call newpen(2)
          call newcolorname('cyan')
          do j = j1, jn, jinc
            plab(1) = rle(1,j) + 0.4*chord(j)
            plab(2) = rle(2,j)
            plab(3) = rle(3,j)
            call viewproj(plab,1,pts_lproj)
            id = j
            lvis=.true.
            if(lhidden) then
              call hidpnt(pts_lproj,id,ngrp,idum,dum,ntri,tri,lvis)
            endif
            ilab = j
            if(lvis) call pltint(pts_lproj,ilab,1.2*chsize,.false.)
          end do

          call newpen(ipen)
          call newcolor(icol)
        endif

c------ vortex labels
        if(label_vrtx) then
          ngrp = 0
          idum(1) = 0
          dum(1,1) = 0.
          call getcolor(icol)
          call getpen(ipen)

          call newpen(2)
          call newcolorname('magenta')
          do j = j1, jn, jinc
            do iv = ifrsts(j), ilasts(j)
              ip = ip+1
              xave = rv(1,iv)
              yave = rv(2,iv)
              zave = rv(3,iv)
              plab(1) = 0.5*(xave+rc(1,iv))
              plab(2) = 0.5*(yave+rc(2,iv))
              plab(3) = 0.5*(zave+rc(3,iv))
              call viewproj(plab,1,pts_lproj)
              id = j
              lvis=.true.
              if(lhidden) then
               call hidpnt(pts_lproj,id,ngrp,idum,dum,ntri,tri,lvis)
              endif
              ilab = iv
              if(lvis) call pltint(pts_lproj,ilab,chsize,.false.)
            end do

            do iv = ifrstu(j), ilastu(j)
              ip = ip+1
              xave = rv(1,iv)
              yave = rv(2,iv)
              zave = rv(3,iv)
              plab(1) = 0.5*(xave+rc(1,iv))
              plab(2) = 0.5*(yave+rc(2,iv))
              plab(3) = 0.5*(zave+rc(3,iv))
              call viewproj(plab,1,pts_lproj)
              id = j
              lvis=.true.
              if(lhidden) then
               call hidpnt(pts_lproj,id,ngrp,idum,dum,ntri,tri,lvis)
              endif
              ilab = iv
              if(lvis) call pltint(pts_lproj,ilab,chsize,.false.)
            end do

            do iv = ifrstw(j), ilastw(j)
              ip = ip+1
              xave = rv(1,iv)
              yave = rv(2,iv)
              zave = rv(3,iv)
              plab(1) = 0.5*(xave+rc(1,iv))
              plab(2) = 0.5*(yave+rc(2,iv))
              plab(3) = 0.5*(zave+rc(3,iv))
              call viewproj(plab,1,pts_lproj)
              id = j
              lvis=.true.
              if(lhidden) then
               call hidpnt(pts_lproj,id,ngrp,idum,dum,ntri,tri,lvis)
              endif
              ilab = iv
              if(lvis) call pltint(pts_lproj,ilab,chsize,.false.)
            end do
          end do
          call newpen(ipen)
          call newcolor(icol)
        endif

 100  continue

c======================================================================
      if(lnrmlplt) then
      do ipass = 1, 2

      call getcolor(icol)
      if(ipass.eq.1) then
c----- with control deflections
       call newcolorname('red')
      else
c----- without control deflections
       call newcolorname('blue')
      endif

c---- go over surfaces again
      do 200 n = 1, nsurf
c------ do only surfaces which are to be plotted
        if(.not.lpltsurf(n)) go to 200

        j1 = jfrst(n)
        jn = jfrst(n) + nj(n)-1
        jinc = 1
        if(idups(n).lt.0) then
         j1 = jfrst(n) + nj(n)-1
         jn = jfrst(n)
         jinc = -1
        endif

c------ normal vectors plotted
        scl = enfac*cref

        ip = 0
        do j = j1, jn, jinc
          do iv = ifrsts(j), ilasts(j)
            entot(1) = enc(1,iv)
            entot(2) = enc(2,iv)
            entot(3) = enc(3,iv)
            if(ipass.eq.1) then
             do k = 1, ncontrol
               entot(1) = entot(1) + enc_d(1,iv,k)*delcon(k)
               entot(2) = entot(2) + enc_d(2,iv,k)*delcon(k)
               entot(3) = entot(3) + enc_d(3,iv,k)*delcon(k)
             enddo
            endif

c---------- define vector from control point in direction of normal vector 
            svec1(1) = entot(1) * scl
            svec1(2) = entot(2) * scl
            svec1(3) = entot(3) * scl
            svec2(1) = 0.02*scl
            svec2(2) = 0.0
            svec2(3) = 0.0

c---------- vector out from control point
            ip = ip+1
            pts_lines(1,1,ip) = rc(1,iv)
            pts_lines(2,1,ip) = rc(2,iv)
            pts_lines(3,1,ip) = rc(3,iv)
            pts_lines(1,2,ip) = rc(1,iv) + svec1(1)  
            pts_lines(2,2,ip) = rc(2,iv) + svec1(2)  
            pts_lines(3,2,ip) = rc(3,iv) + svec1(3)  
            id_lines(ip) = 0

c---------- arrow head for normal vector
            ip = ip+1
            pts_lines(1,1,ip) = rc(1,iv) + svec1(1)  
            pts_lines(2,1,ip) = rc(2,iv) + svec1(2)  
            pts_lines(3,1,ip) = rc(3,iv) + svec1(3)  
            pts_lines(1,2,ip) = rc(1,iv) + 0.8*svec1(1) + svec2(1) 
            pts_lines(2,2,ip) = rc(2,iv) + 0.8*svec1(2) + svec2(2) 
            pts_lines(3,2,ip) = rc(3,iv) + 0.8*svec1(3) + svec2(3)  
            id_lines(ip) = 0
            ip = ip+1
            pts_lines(1,1,ip) = rc(1,iv) + svec1(1)  
            pts_lines(2,1,ip) = rc(2,iv) + svec1(2)  
            pts_lines(3,1,ip) = rc(3,iv) + svec1(3)  
            pts_lines(1,2,ip) = rc(1,iv) + 0.8*svec1(1) - svec2(1) 
            pts_lines(2,2,ip) = rc(2,iv) + 0.8*svec1(2) - svec2(2) 
            pts_lines(3,2,ip) = rc(3,iv) + 0.8*svec1(3) - svec2(3)  
            id_lines(ip) = 0
          end do

          do iv = ifrstu(j), ilastu(j)
            entot(1) = enc(1,iv)
            entot(2) = enc(2,iv)
            entot(3) = enc(3,iv)
            if(ipass.eq.1) then
             do k = 1, ncontrol
               entot(1) = entot(1) + enc_d(1,iv,k)*delcon(k)
               entot(2) = entot(2) + enc_d(2,iv,k)*delcon(k)
               entot(3) = entot(3) + enc_d(3,iv,k)*delcon(k)
             enddo
            endif

c---------- define vector from control point in direction of normal vector 
            svec1(1) = entot(1) * scl
            svec1(2) = entot(2) * scl
            svec1(3) = entot(3) * scl
            svec2(1) = 0.02*scl
            svec2(2) = 0.0
            svec2(3) = 0.0

c---------- vector out from control point
            ip = ip+1
            pts_lines(1,1,ip) = rc(1,iv)
            pts_lines(2,1,ip) = rc(2,iv)
            pts_lines(3,1,ip) = rc(3,iv)
            pts_lines(1,2,ip) = rc(1,iv) + svec1(1)  
            pts_lines(2,2,ip) = rc(2,iv) + svec1(2)  
            pts_lines(3,2,ip) = rc(3,iv) + svec1(3)  
            id_lines(ip) = 0

c---------- arrow head for normal vector
            ip = ip+1
            pts_lines(1,1,ip) = rc(1,iv) + svec1(1)  
            pts_lines(2,1,ip) = rc(2,iv) + svec1(2)  
            pts_lines(3,1,ip) = rc(3,iv) + svec1(3)  
            pts_lines(1,2,ip) = rc(1,iv) + 0.8*svec1(1) + svec2(1) 
            pts_lines(2,2,ip) = rc(2,iv) + 0.8*svec1(2) + svec2(2) 
            pts_lines(3,2,ip) = rc(3,iv) + 0.8*svec1(3) + svec2(3)  
            id_lines(ip) = 0
            ip = ip+1
            pts_lines(1,1,ip) = rc(1,iv) + svec1(1)  
            pts_lines(2,1,ip) = rc(2,iv) + svec1(2)  
            pts_lines(3,1,ip) = rc(3,iv) + svec1(3)  
            pts_lines(1,2,ip) = rc(1,iv) + 0.8*svec1(1) - svec2(1) 
            pts_lines(2,2,ip) = rc(2,iv) + 0.8*svec1(2) - svec2(2) 
            pts_lines(3,2,ip) = rc(3,iv) + 0.8*svec1(3) - svec2(3)  
            id_lines(ip) = 0
          end do

          do iv = ifrstw(j), ilastw(j)
            entot(1) = enc(1,iv)
            entot(2) = enc(2,iv)
            entot(3) = enc(3,iv)
            if(ipass.eq.1) then
             do k = 1, ncontrol
               entot(1) = entot(1) + enc_d(1,iv,k)*delcon(k)
               entot(2) = entot(2) + enc_d(2,iv,k)*delcon(k)
               entot(3) = entot(3) + enc_d(3,iv,k)*delcon(k)
             enddo
            endif

c---------- define vector from control point in direction of normal vector 
            svec1(1) = entot(1) * scl
            svec1(2) = entot(2) * scl
            svec1(3) = entot(3) * scl
            svec2(1) = 0.02*scl
            svec2(2) = 0.0
            svec2(3) = 0.0

c---------- vector out from control point
            ip = ip+1
            pts_lines(1,1,ip) = rc(1,iv)
            pts_lines(2,1,ip) = rc(2,iv)
            pts_lines(3,1,ip) = rc(3,iv)
            pts_lines(1,2,ip) = rc(1,iv) + svec1(1)  
            pts_lines(2,2,ip) = rc(2,iv) + svec1(2)  
            pts_lines(3,2,ip) = rc(3,iv) + svec1(3)  
            id_lines(ip) = 0

c---------- arrow head for normal vector
            ip = ip+1
            pts_lines(1,1,ip) = rc(1,iv) + svec1(1)  
            pts_lines(2,1,ip) = rc(2,iv) + svec1(2)  
            pts_lines(3,1,ip) = rc(3,iv) + svec1(3)  
            pts_lines(1,2,ip) = rc(1,iv) + 0.8*svec1(1) + svec2(1) 
            pts_lines(2,2,ip) = rc(2,iv) + 0.8*svec1(2) + svec2(2) 
            pts_lines(3,2,ip) = rc(3,iv) + 0.8*svec1(3) + svec2(3)  
            id_lines(ip) = 0
            ip = ip+1
            pts_lines(1,1,ip) = rc(1,iv) + svec1(1)  
            pts_lines(2,1,ip) = rc(2,iv) + svec1(2)  
            pts_lines(3,1,ip) = rc(3,iv) + svec1(3)  
            pts_lines(1,2,ip) = rc(1,iv) + 0.8*svec1(1) - svec2(1) 
            pts_lines(2,2,ip) = rc(2,iv) + 0.8*svec1(2) - svec2(2) 
            pts_lines(3,2,ip) = rc(3,iv) + 0.8*svec1(3) - svec2(3)  
            id_lines(ip) = 0
          end do
        end do
        nlines = ip
        nproj = 2*nlines
        call viewproj(pts_lines,nproj,pts_lproj)

        call plotlines(nlines,pts_lproj,id_lines)
 200  continue

      call newcolor(icol)
      enddo
      endif

c---- go over bodies
      do 300 ib = 1, nbody
c------ do only bodies which are to be plotted
        if(.not.lpltbody(ib)) go to 300

        l1 = lfrst(ib)
        ln = llast(ib)

c------ source lines plotted
        ip = 0
        do l = l1, ln
          ip = ip+1
          pts_lines(1,1,ip) = rl1(1,l)
          pts_lines(2,1,ip) = rl1(2,l)
          pts_lines(3,1,ip) = rl1(3,l)
          pts_lines(1,2,ip) = rl2(1,l)
          pts_lines(2,2,ip) = rl2(2,l)
          pts_lines(3,2,ip) = rl2(3,l)
          id_lines(ip) = 0
        end do

        nlines = ip
        nproj = 2*nlines
        call viewproj(pts_lines,nproj,pts_lproj)

        call getpen(ipen)
        call getcolor(icol)
        call newpen(4)
        call newcolorname('magenta')
        call plotlines(nlines,pts_lproj,id_lines)
        call newcolor(icol)
        call newpen(ipen)


c------- plot body hoops
ccc        if(lchordplt) then
         if(iysym .eq. 1 .and.
     &      ibcent(ib) .eq. 1) then
c------- y-symmetry with body at y=0: plot only right side of body
          nks = nhoop/2
         else
c------- plot entire body
          nks = nhoop
         endif
         ip = 0
         do l = l1, ln
           dxl = rl2(1,l) - rl1(1,l)
           dyl = rl2(2,l) - rl1(2,l)
           dzl = rl2(3,l) - rl1(3,l)
           dsl = sqrt(dxl**2 + dyl**2 + dzl**2)
           if(dsl.ne.0.0) then
            dxl = dxl/dsl
            dyl = dyl/dsl
            dzl = dzl/dsl
           endif
c     
c--------- "vertical" body section unit vector
           uvx = -dzl
           uvy = 0.
           uvz =  dxl

c--------- "horizontal" body section unit vector
           uhx = uvy*dzl - uvz*dyl
           uhy = uvz*dxl - uvx*dzl
           uhz = uvx*dyl - uvy*dxl

           do k = 1, nks
             thk1 = 2.0*pi*float(k-1)/float(nhoop)
             thk2 = 2.0*pi*float(k  )/float(nhoop)

             rx1 = radl(l)*(uvx*cos(thk1) + uhx*sin(thk1))
             ry1 = radl(l)*(uvy*cos(thk1) + uhy*sin(thk1))
             rz1 = radl(l)*(uvz*cos(thk1) + uhz*sin(thk1))

             rx2 = radl(l)*(uvx*cos(thk2) + uhx*sin(thk2))
             ry2 = radl(l)*(uvy*cos(thk2) + uhy*sin(thk2))
             rz2 = radl(l)*(uvz*cos(thk2) + uhz*sin(thk2))

             ip = ip+1
             pts_lines(1,1,ip) = rl(1,l) + rx1
             pts_lines(2,1,ip) = rl(2,l) + ry1
             pts_lines(3,1,ip) = rl(3,l) + rz1
             pts_lines(1,2,ip) = rl(1,l) + rx2
             pts_lines(2,2,ip) = rl(2,l) + ry2
             pts_lines(3,2,ip) = rl(3,l) + rz2
             id_lines(ip) = 0
           enddo
         enddo

         nlines = ip
         nproj = 2*nlines

         call getpen(ipen)
         call newpen(2)
         call viewproj(pts_lines,nproj,pts_lproj)
         call plotlines(nlines,pts_lproj,id_lines)
         call newpen(ipen)
ccc        endif

c------ source line loading plotted
        if(lloadplt) then
         call getpen(ipen)
         call getcolor(icol)

         call newpen(5)
         call newcolorname('green')

c--- plot the load vectors (up from surface)
         ip = 0
         cpscl = cpfac*cref
         do l = l1, ln
           ip = ip+1
           xave = rl(1,l)
           yave = rl(2,l)
           zave = rl(3,l)
           xload = xave + dcpb(1,l) * cpscl
           yload = yave + dcpb(2,l) * cpscl
           zload = zave + dcpb(3,l) * cpscl

           pts_lines(1,1,ip) = xave
           pts_lines(2,1,ip) = yave
           pts_lines(3,1,ip) = zave
           pts_lines(1,2,ip) = xload
           pts_lines(2,2,ip) = yload
           pts_lines(3,2,ip) = zload
           id_lines(ip) = 0
         end do
         nlines = ip
         nproj = 2*nlines
         call viewproj(pts_lines,nproj,pts_lproj)
         call plotlines(nlines,pts_lproj,id_lines)

         call newpen(ipen)
         call newcolor(icol)
        endif


c*************************************************************************
c---now do any labels (integer indices for surfaces,strips or vortices)

c--- body labels
        if(label_body) then
          jmid = jfrst(n) + ifix(0.5*float(nj(n))) - 1
          plab(1) = rle(1,jmid) + 0.6*chord(jmid)
          plab(2) = rle(2,jmid)
          plab(3) = rle(3,jmid)
          call viewproj(plab,1,pts_lproj)
          id = jmid
          lvis=.true.
          if(lhidden) then
            ngrp = 0
            idum(1) = 0
            dum(1,1) = 0.
            call hidpnt(pts_lproj,id,ngrp,idum,dum,ntri,tri,lvis)
          endif
          ilab = n
          if (lvis) call pltint(pts_lproj,ilab,1.5*chsize,.false.)
        endif

 300  continue

      return
      end ! plotgeom


      subroutine plotlines(nlines,pts_lproj,id_lines)
c...purpose:    plot line segments with hidden line 
c               routine.
      include 'jvl.inc'
      include 'jvlplt.inc'

      real pts_lproj(3,2,*), id_lines(*)
c     
      parameter (nvx=32)
      real alf(2,nvx)

      integer idum(1)
      real dum(1,1)

      alf(1,1) = 0.0
      alf(2,1) = 1.0
      nalf = 1

      do l=1, nlines
c-----  check for hidden line treatment
        if(lhidden) then
          nalf = nvx
          ngrp = 0
          id = id_lines(l)
          idum(1) = 0
          dum(1,1) = 0.0
          call hidlin(pts_lproj(1,1,l),id,ngrp,idum,dum,
     &                ntri,tri,nalf,alf)
        endif
c-----  if any of the line segment is visible, plot the visible pieces
        if(nalf.gt.0) call pltseg(pts_lproj(1,1,l),alf,nalf)
      end do

      return
      end ! plotlines
