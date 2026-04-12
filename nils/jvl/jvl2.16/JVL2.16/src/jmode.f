C***********************************************************************
C    Module:  jmode.f
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

      subroutine mode
c---------------------------------------------
c     Flight dynamics analysis driver
c---------------------------------------------
      include 'jvl.inc'
      include 'jvlplt.inc'
      logical error, lok, saved, overlay, lowrit

      character*1 chkey
      character*2 opt, chplot
      character*4 comand, itemc
      character*80 fnnew, fnsys
      character*80 comarg, prompt, crun

      logical lproot(jemax,nrmax),
     &        lprnum(jemax,nrmax)
      save lproot,lprnum

      real*8 asys(jemax,jemax),
     &       bsys(jemax,ndmax+njmax),rsys(jemax)
      complex zcroot

      real rinput(20), rinp(20)
      integer iinput(20), iinp(20)

      irun0 = irun

 1000 format (a)

c---- no viewed eigenmode yet
      keview = 0
      irview = 0

c---- force computation of plot limits
      tmin = 0.
      tmax = 0.
      fmin = 0.
      fmax = 0.

      lplot = .false.
      chplot = '  '

c      eyeptz = 180.
c      eyepty = 90.
c      eyeptx = 0.
c      robinv = 0.

      overlay = .false.

c=================================================================
c---- start of user interaction loop
 800  continue
      call cfrac(irune,nrun,crun,npr)

      write(*,1050) 
 1050 format(
     &  /'  Run-case parameters for eigenmode analyses ... ')

      lu = 6
      call runlst(lu,irune)

      write(*,1052) dvjexp
 1052 format(
     &   ' =========================================================='
     & //' "#" select run case for eigenmode analysis (0 = all)'
     & //'  M odify parameters'
     & //'  N ew eigenmode calculation'
     & //'  J et velocity-excess exponent  (currently =', f8.4,')'
     & //'  P lot root locus'
     &  /'  B lowup window'
     &  /'  R eset to normal size'
     &  /' eX amine selected eigenmode'
     & //'  A nnotate current plot'
     &  /'  H ardcopy current plot'
     &  /'  T ime-integration parameters'
     & //'  S ystem matrix output'
     &  /'  W rite eigenvalues to file'
     &  /'  D ata file overlay toggle'
     & //'  Z oom'
     &  /'  U nzoom')

c   A B C D E F G H I J K L M N O P Q R S T U V W X Y Z
c   x x   x       x   x     x x   x   x x x x   x x   x

 810  continue
      if(irune .eq. 0) then
       prompt = ' .MODE (all cases)^'
      else
       prompt = ' .MODE (case ' // crun(1:npr) // ')^'
      endif
      call askc(prompt,comand,comarg)

 815  continue
c------------------------------------------------------
      if    (comand.eq.'    ') then
       if(lplot) call plend
       lplot = .false.
       call clrzoom
       irun = irun0
       return

      elseif(comand.eq.'?   ') then
       go to 800

      endif

c-------------------------------------------------------------------
c---- see if command is an integer
      ninput = 1
      call getint(comand,iinput,ninput,error)
      if(.not.error .and. ninput.ge.1 
     &  .and. comand(1:1).ne.'T'
     &  .and. comand(1:1).ne.'F' ) then
c----- command is an integer... new case index?
       if(iinput(1).lt.0 .or. iinput(1).gt.nrun) then
        write(*,*)
        write(*,*) '* Selected new run case is not defined'
        go to 800
       else
c------ valid new run case selected... go back to menu
        irune = iinput(1)
        go to 800
       endif
      endif

c-------------------------------------------------------------------
c---- extract command line numeric arguments
      do i=1, 20
        iinput(i) = 0
        rinput(i) = 0.0
      enddo
      ninput = 20
      call getint(comarg,iinput,ninput,error)
      ninput = 20
      call getflt(comarg,rinput,ninput,error)

c---- set up for multiple or single run cases
      if(irune.eq.0) then
       irun1 = 1
       irun2 = nrun
      else
       irun1 = irune
       irun2 = irune
      endif

c-----------------------------------------------------------------------
c      if    (comand.eq.'R   ') then
cc----- set moment reference point to cg location
c       dxcref = (xyzmass(1)/unitl - xyzref(1)) / cref
c       write(*,*)
c       write(*,1210) '  Previous XYZref  = ', (xyzref(k), k=1,3),'lunit'
c       xyzref(1) = xyzmass(1)/unitl
c       xyzref(2) = xyzmass(2)/unitl
c       xyzref(3) = xyzmass(3)/unitl
c       write(*,1210) '  New (cg) XYZref  = ', (xyzref(k), k=1,3),'lunit'
cc
c       write(*,*)
c       write(*,1210) '  delta(Xref)/Cref = ', dxcref
c 1210  format(1x, a, 3g10.4,3x,a)
cc
c-----------------------------------------------------------------------
      if    (comand.eq.'M   ') then
       call parmod(irune)

c-------------------------------------------------------------------
      elseif(comand .eq. 'N   ' .or.
     &       comand .eq. 'C   '      ) then
c------ execute eigenmode calculation
        do 100 ir = irun1, irun2
          call runchk(ir,lok)
          if(.not.lok) then
           write(*,*) '** Skipping ill-posed run case', ir
           go to 100
          endif

          niter = 10
          info = 0
          call exec(niter,info,ir)

          if(comand .eq. 'N   ') then
           call sysmat(ir,asys,bsys,rsys,nsys)
          else
           call appmat(ir,asys,bsys,rsys,nsys)
          endif

cc          lu = 6
cc          call syssho(lu,asys,bsys,rsys,nsys)

          info = 1
          etol = 1.0e-5
          call eigsol(info,ir,etol,asys,nsys)

c
          do keig = 1, neigen(ir)
            lproot(keig,ir) = .true.
            lprnum(keig,ir) = irun1 .ne. irun2
          enddo

          xorg0 = 0.0
          yorg0 = 0.0
          tmin = 0.
          tmax = 0.
          fmin = 0.
          fmax = 0.
          call plemap(xorg0,yorg0, 
     &                irun1,ir,
     &                keview,irview,
     &                lproot,lprnum,overlay,
     &                tmin,tmax,tdel,tfac,
     &                fmin,fmax,fdel,ffac )

          write(*,*)
          call eiglst(6,ir)

 100    continue

        chplot = 'p '

c-------------------------------------------------------------------
      elseif(comand .eq. 'J   ') then
        if(ninput .lt. 1) then
         rinput(1) = dvjexp
         write(*,3010) rinput(1)
 3010    format(/' Enter exponent of DVjet variation with Vinf:',
     &           f8.4)
         call readr(1,rinput(1),error)
         if(error) then
          go to 810
         endif
        endif

        dvjexp = rinput(1)

c        write(*,3020) dvjexp
c 3020   format(/1x,' DVjet ~ (Vinf/Vref)^', f8.4)
c        go to 810

        go to 800

c-------------------------------------------------------------------
      elseif(comand .eq. 'P   ') then
        do ir = irun1, irun2
          if(neigen(ir) .gt. 0) go to 51
        enddo
        write(*,*) 'First compute eigenmodes with N command'
        go to 810

 51     xorg0 = 0.0
        yorg0 = 0.0
        call plemap(xorg0,yorg0, 
     &              irun1,irun2,
     &              keview,irview,
     &              lproot,lprnum,overlay,
     &              tmin,tmax,tdel,tfac,
     &              fmin,fmax,fdel,ffac )
        chplot = comand(1:2)

c-------------------------------------------------------------------
      elseif(comand .eq. 'X   ') then
       do ir = irun1, irun2
         if(neigen(ir) .gt. 0) go to 61
       enddo
       write(*,*) 'first compute eigenmodes with n command'
       go to 810

 61    continue
       if(index('p',chplot(1:1)) .eq. 0) then
c------ make root map plot first
        xorg0 = 0.0
        yorg0 = 0.0
        tmin = 0.
        tmax = 0.
        fmin = 0.
        fmax = 0.
        call plemap(xorg0,yorg0, 
     &              irun1,irun2,
     &              keview,irview,
     &              lproot,lprnum,overlay,
     &              tmin,tmax,tdel,tfac,
     &              fmin,fmax,fdel,ffac )
       endif

       write(*,*)
       write(*,*) 'Click on root of eigenmode to be viewed...'

       call getcursorxy(xc1,yc1,chkey)
       xe1 = (xc1-xorg0)/tfac + tmin
       ye1 = (yc1-yorg0)/ffac + fmin

       if(xe1 .lt. tmin .or. xe1 .gt. tmax  .or.
     &    ye1 .lt. fmin .or. ye1 .gt. fmax       ) then
c------ click was outside of grid... manual input
 65     rinput(1) = 0.
        rinput(2) = 0.
        write(*,1165) rinput(1), rinput(2)
 1165   format(/' Enter  sigma,omega  of root to be viewed:', 2f12.6)
        call readr(2,rinput,error)
        if(error) go to 65

        zcroot = cmplx( rinput(1) , rinput(2) )

       else
c------ set root from its symbol-plot coordinates
        zcroot = cmplx( xe1 , ye1 )

       endif

c----- search current roots for one nearest to specified location
       ir = 0
       keig = 0
       dist = 1.0e32
       do jr = irun1, irun2
         do k = 1, neigen(jr)
           dist1 = abs( zcroot - eval(k,jr) )
           if(dist1.le.dist) then
            ir = jr
            keig = k
            dist = dist1
           endif
         enddo
       enddo

       if(ir.eq.0 .or. keig.eq.0) then
        write(*,*) 'Nearest root not identified'
        go to 810
       endif

       call plotmd(azimob, elevob, tiltob, robinv, keig, ir)

       keview = keig
       irview = ir

       chplot = comand(1:2)

c-------------------------------------------------------------------
c---- blowup window
      elseif(comand.eq.'B   ') then
       if(index('P',chplot(1:1)) .eq. 0) then
        write(*,*) 'No plot to blow up'
        go to 810
       endif

       dcross = 2.0
       write(*,*) 'Mark off corners of blowup area'
       call getcursorxy(xc1,yc1,chkey)
       call plsymb(xc1,yc1,dcross,3,0.0,0)
       call plflush
       call getcursorxy(xc2,yc2,chkey)
       call plsymb(xc2,yc2,dcross,3,0.0,0)
       call plflush

       xe1 = tmin+(xc1-xorg0)/tfac
       ye1 = fmin+(yc1-yorg0)/ffac
       xe2 = tmin+(xc2-xorg0)/tfac
       ye2 = fmin+(yc2-yorg0)/ffac

       tmin = min( xe1 , xe2 )
       tmax = max( xe1 , xe2 )
       fmin = min( ye1 , ye2 )
       fmax = max( ye1 , ye2 )

       call axisadj(tmin,tmax, tspan, tdel, nann)
       call axisadj(fmin,fmax, fspan, fdel, nann)

       comand = 'p   '
       go to 815

c-------------------------------------------------------------------
c---- normal size
      elseif(comand.eq.'R   ') then
       tmin = 0.
       tmax = 0.
       fmin = 0.
       fmax = 0.

       comand = 'P   '
       go to 815

c-------------------------------------------------------------------
c---- annotate
      elseif(comand.eq.'A   ') then
       if(lplot) then
        call annot(csize)
       else
        write(*,*) 'No active plot'
       endif

c-------------------------------------------------------------------
c---- hardcopy
      elseif(comand.eq.'H   ') then
       if(lplot) call plend
       lplot = .false.
       call replot(idevh)

c-------------------------------------------------------------------
c---- change time integration parameters
      elseif(comand.eq.'T   ') then
 70     write(*,2070) dtimed, dtmovie, tmovie
 2070   format(/' ==========================='
     &         /'  I ntegration  delta(t): ', g11.5,
     &         /'  O utput-frame delta(t): ', g11.5,
     &         /'  M aximum movie time   : ', g11.5 )

 72     call askc(' Select item,value^',itemc,comarg)
 2072   format(' Enter new ',a,': ', $)

        if(itemc.eq.'    ') then
          go to 810
        endif

        ninp = 1
        call getflt(comarg,rinp,ninp,error)

        if    (index('II',itemc(1:1)).ne.0) then
          if(ninp.eq.0) then
 74        write(*,2072) 'Integration delta(t)'
           read (*,*,err=74) dtimed
          else
           dtimed = rinp(1)
          endif

        elseif(index('OO',itemc(1:1)).ne.0) then
          if(ninp.eq.0) then
 75        write(*,2072) 'Output-frame delta(t)'
           read (*,*,err=75) dtmovie
          else
           dtmovie = rinp(1)
          endif

        elseif(index('MM',itemc(1:1)).ne.0) then
          if(ninp.eq.0) then
 76        write(*,2072) 'Maximum movie time'
           read (*,*,err=76) tmovie
          else
           tmovie = rinp(1)
          endif

        else
          write(*,*) 'Item not recognized'
          go to 70

        endif
        go to 70

c-------------------------------------------------------------------
c---- write system matrices
      elseif(comand.eq.'S   ' .or.
     &       comand.eq.'SC  '      ) then
        do 80 ir = irun1, irun2
          write(*,2080) ir
 2080     format(/' Run case', i3,'...')

          call runchk(ir,lok)
          if(.not.lok) then
           write(*,*) '** Skipping ill-posed run case', ir
           go to 80
          endif

          alfa    = parval(ipalfa,ir)*dtr
          beta    = parval(ipbeta,ir)*dtr
          wbar(1) = parval(iprotx,ir)*2.0/bref
          wbar(2) = parval(iproty,ir)*2.0/cref
          wbar(3) = parval(iprotz,ir)*2.0/bref

          xyzref(1) = parval(ipxcg,ir)
          xyzref(2) = parval(ipycg,ir)
          xyzref(3) = parval(ipzcg,ir)

          cdref     = parval(ipcd0,ir)

          hrotor(1) = parval(iphx,ir)
          hrotor(2) = parval(iphy,ir)
          hrotor(3) = parval(iphz,ir)

          niter = 10
          info = 0
          call exec(niter,info,ir)

          if(comand.eq.'S   ') then
           call sysmat(ir,asys,bsys,rsys,nsys)
          else
           call appmat(ir,asys,bsys,rsys,nsys)
          endif

          lu = 6
          call syssho(lu,asys,bsys,rsys,nsys)

          write(*,2082)
 2082     format(/' Enter output filename (or <return>): ', $)
          read(*,1000) fnsys
          call bstrip(fnsys,nfs)
          if(nfs.eq.0) go to 80

          open(lusys,file=fnsys,status='old',err=81)
          if(lowrit(fnsys)) then
           rewind(lusys)
           go to 82
          else
           close(lusys)
           go to 80
          endif

 81       open(lusys,file=fnsys,status='new',err=85)
 82       call syssho(lusys,asys,bsys,rsys,nsys)
          close(lusys)
          go to 800

 85       continue
          write(*,*) '* File OPEN error'

 80     continue

c-------------------------------------------------------------------
c---- write eigenvalues
      elseif(comand.eq.'W   ') then
       if(fevdef(1:1).eq.' ') then
c------ set default filename
        kdot = index(fildef,'.')
        if(kdot.eq.0) then
         call slen(fildef,nfil)
         fevdef = fildef(1:nfil) // '.eig'
        else
         fevdef = fildef(1:kdot) // 'eig'
        endif
       endif

       call slen(fevdef,nfe)
       nfe = max( nfe , 1 )
       write(*,2040) fevdef(1:nfe)
 2040  format(' Enter eigenvalue save filename: ', a)
       read (*,1000) fnnew
       if(fnnew.ne.' ') fevdef = fnnew
       call eigout(fevdef, irun1,irun2, saved)

c-------------------------------------------------------------------
c---- toggle overlay flag
      elseif(comand.eq.'D   ') then
       overlay = .not.overlay

       if(overlay) then
        if(fevdef(1:1).eq.' ') then
c------- set default filename
         kdot = index(fildef,'.')
         if(kdot.eq.0) then
          call slen(fildef,nfil)
          fevdef = fildef(1:nfil) // '.eig'
         else
          fevdef = fildef(1:kdot) // 'eig'
         endif
        endif

        call slen(fevdef,nfe)
        nfe = max( nfe , 1 )
        write(*,2050) fevdef(1:nfe)
 2050   format(' Enter eigenvalue data filename: ', a)
        read (*,1000) fnnew
        if(fnnew.ne.' ') fevdef = fnnew

c------ read the data file
        call eiginp(fevdef, error)
        if(error) then
         overlay = .false.
         go to 800
        endif

        overlay = .true.
        write(*,*) 'Overlay plotting enabled'

        if(index('P',chplot(1:1)).ne.0) then
         comand = chplot
         go to 815
        endif

       else
        write(*,*) 'Overlay plotting disabled'

       endif

c-------------------------------------------------------------------
c---- zoom in on plot
      elseif(comand.eq.'Z   ') then
       if(index('P',chplot(1:1)) .eq. 0) then
        write(*,*) 'No plot to zoom'
        go to 810
       endif

       call usetzoom(.true.,.true.)
       call replot(idev)

c-------------------------------------------------------------------
c---- reset zoom on plot
      elseif(comand.eq.'U   ') then
       call clrzoom
       call replot(idev)

c-------------------------------------------------------------------
      else
        write(*,*)
        write(*,*) '* option not recognized'

      endif
      go to 800

      end ! mode


      subroutine eiglst(lu,ir)
      include 'jvl.inc'

      write(lu,3100) ir, rtitle(ir)

      do j = 1, neigen(ir)
      write(lu,3200) j, real(eval(j,ir)), imag(eval(j,ir))
      write(lu,3300)
     & 'u  ', real(evec(jeu,j,ir)), imag(evec(jeu,j,ir)),
     & 'v  ', real(evec(jev,j,ir)), imag(evec(jev,j,ir)),
     & 'x  ', real(evec(jex,j,ir)), imag(evec(jex,j,ir))
      write(lu,3300)
     & 'w  ', real(evec(jew,j,ir)), imag(evec(jew,j,ir)),
     & 'p  ', real(evec(jep,j,ir)), imag(evec(jep,j,ir)),
     & 'y  ', real(evec(jey,j,ir)), imag(evec(jey,j,ir))
      write(lu,3300)
     & 'q  ', real(evec(jeq,j,ir)), imag(evec(jeq,j,ir)),
     & 'r  ', real(evec(jer,j,ir)), imag(evec(jer,j,ir)),
     & 'z  ', real(evec(jez,j,ir)), imag(evec(jez,j,ir))
      write(lu,3300)
     & 'the', real(evec(jeth,j,ir)), imag(evec(jeth,j,ir)),
     & 'phi', real(evec(jeph,j,ir)), imag(evec(jeph,j,ir)),
     & 'psi', real(evec(jeps,j,ir)), imag(evec(jeps,j,ir))
      enddo

 3100 format(/1x,'Run case',i3,':  ',a)
 3200 format(/1x,' Mode',i2,':', 2g14.6)
 3300 format( 1x, a,':', 2f11.4, 6x, a,':', 2f11.4, 6x, a,':', 2g12.4)

      return
      end ! eiglst



      subroutine eigout(filnam, irun1,irun2, saved)
      include 'jvl.inc'

      character*(*) filnam
      logical saved

      logical lowrit

      lu = 2
      open(lu,file=filnam,status='old',err=17)
      if(lowrit(filnam)) then
       go to 18
      else
       close(lu)
       write(*,*) 'Eigenvalues not saved.'
       saved = .false.
       return
      endif

 17   open(lu,file=filnam,status='unknown',err=80)
 18   continue
      rewind(lu)

c---- write header
      write(lu,1100) title
 1100 format('# ', a
     &      /'# '
     &      /'#   Run case     Eigenvalue')

 2300 format(1x, i7, 2g18.8, f12.4 )

      do ir = irun1, irun2
        do keig = 1, neigen(ir)
          write(lu,2300) ir,
     &                 real(eval(keig,ir)),
     &                 imag(eval(keig,ir))
        enddo
      enddo

      close(lu)
      saved = .true.
      return

 80   continue
      saved = .false.
      return
      end ! eigout




      subroutine eiginp(filnam, error)
      include 'jvl.inc'

      character*(*) filnam
      logical error

      character*1 dummy

c---- assume that read will fail
      error = .true.

 1000 format(a)

      do ir = 1, nrmax
        neigendat(ir) = 0
      enddo

      lu = 2
      open(lu,file=filnam,status='old',err=85)
      rewind(lu)

c---- read header
      read(lu,1000,err=90) dummy
      read(lu,1000,err=90) dummy
      read(lu,1000,err=90) dummy

c---- find index of forcing parameter of this file
      do iline = 1, 123456
 15     read(lu,*,end=20,err=90) ir, sigma, omega
        if(ir.lt.1 .or. ir.gt.nrmax) go to 15

        keig = neigendat(ir) + 1
        if(keig .le. jemax) then
         evaldat(keig,ir) = cmplx(sigma,omega)
         neigendat(ir) = keig
        endif
      enddo

 20   continue
      close(lu)
      error = .false.
      return

 85   continue
      write(*,*) 'file open error'
      return

 90   continue
      write(*,*) 'file read error'
      close(lu)
      return
      end ! eiginp





      subroutine runchk(jrun,ok)
c------------------------------------------------------------------
c     returns ok=t if run case jrun has no redundant constraints.
c------------------------------------------------------------------
      include 'jvl.inc'
      logical ok

      ok = .true.
      do iv = 2, nvtot
        do jv = 2, nvtot
          if(iv.ne.jv .and. icon(iv,jrun).eq.icon(jv,jrun)) then
           ok = .false.
          endif
        enddo
      enddo

      return
      end ! runchk



      subroutine sysmat(ir,asys,bsys,rsys,nsys)
c------------------------------------------------------------------
c     Computes the 12 dx/dt components of the state vector x(t)
c     from the current geometry-axis coefficients and derivatives, 
c     and from the dimensional parameters for run case ir.
c
c     Also computes the instantaneous A and B matrices.
c
c     All quantities are in JVL's geometry axes.
c     To convert to standard body axes of flight dynamics, 
c     negate all x and z components, i.e.  u,w, p,r, Rx,Rz
c
c     Inputs (from the common blocks in jvl.inc):
c       cf(.)     Cx,Cy,Cz
c       cm(.)     Cl,Cm,Cn
c
c       cf_v(..)  Cxu,Cyu,Czu, Cxv,Cyv ...   Czw
c       cf_w(..)  Cxp,Cyp,Czp, Cxq,Cyq ...   Czr
c       cf_d(..)  Cxd1,Cyd1 ...  Czdn
c       cf_j(..)  Cxj1,Cyj1 ...  Czjn
c
c       cm_v(..)  Clu,Cmu ...  Cnw
c       cm_w(..)  Clp,Cmp ...  Cnr
c
c       vbar(.)   Vinf /Vref = -U/Vref = (-u/Vref,-v/Vref,-w/Vref)
c       wbar(.)   Omega/Vref           = ( p/Vref, q/Vref, r/Vref)
c       parval(.ir)  dimensional parameters
c
c     Ouputs:
c       rsys(.)   dx/dt components
c       asys(..)  A matrix
c       bsys(..)  B matrix
c
c------------------------------------------------------------------
      include 'jvl.inc'
      real*8 asys(jemax,jemax),
     &       bsys(jemax,ndmax+njmax),rsys(jemax)

      real v(3), w(3)

      real ang(3), tt(3,3), tt_ang(3,3,3), rt(3,3), rt_ang(3,3,3)

      real lrefd(3)
      real vd(3), wd(3)
      real qsd, qsd_vd(3)

      real cf_vd(3,3),
     &     cf_wd(3,3)
      real cm_vd(3,3),
     &     cm_wd(3,3)

      real fd(3), fd_vd(3,3), fd_wd(3,3), fd_d(3,ndmax), fd_j(3,njmax)
      real md(3), md_vd(3,3), md_wd(3,3), md_d(3,ndmax), md_j(3,njmax)

      real riner(3,3),
     &     mamat(3,3),
     &     rimat(3,3)
      real mmat(6,6),
     &     minv(6,6)
      real p(3), p_vd(3,3),
     &     h(3), h_wd(3,3)
      real wxp(3), wxp_vd(3,3), wxp_wd(3,3),
     &     wxh(3), wxh_vd(3,3), wxh_wd(3,3)
      real mif(3), mif_vd(3,3), mif_wd(3,3),
     &             mif_d(3,ndmax), mif_j(3,njmax),
     &     rim(3), rim_vd(3,3), rim_wd(3,3), 
     &             rim_d(3,ndmax), rim_j(3,njmax)
      real prf(3), prf_vd(3,3), prf_wd(3,3),
     &     prm(3), prm_vd(3,3), prm_wd(3,3)
      real mg(3) , mg_ang(3,3),
     &     mgf(3), mgf_ang(3,3),
     &     mgm(3), mgm_ang(3,3)
      real drm(3)

      integer icrs(3), jcrs(3)
      data icrs / 2, 3, 1 / , jcrs / 3, 1, 2 /

c---- dimensional parameters for run case ir
      vrefd = parval(ipvee,ir)
      gee = parval(ipgee,ir)
      rho = parval(iprho,ir)
      phi = parval(ipphi,ir)
      the = parval(ipthe,ir)
      psi = parval(ippsi,ir)
      xcg = parval(ipxcg,ir)
      ycg = parval(ipycg,ir)
      zcg = parval(ipzcg,ir)
      rmass = parval(ipmass,ir)
      riner(1,1) = parval(ipixx,ir)
      riner(2,2) = parval(ipiyy,ir)
      riner(3,3) = parval(ipizz,ir)
      riner(1,2) = parval(ipixy,ir)
      riner(2,3) = parval(ipiyz,ir)
      riner(3,1) = parval(ipizx,ir)
      riner(2,1) = parval(ipixy,ir)
      riner(3,2) = parval(ipiyz,ir)
      riner(1,3) = parval(ipizx,ir)

      hrotor(1) = parval(iphx,ir)
      hrotor(2) = parval(iphy,ir)
      hrotor(3) = parval(iphz,ir)

c---- added stability derivatives from aero lags
      dcl_adot = parval(ipclad,ir)
      dcd_adot = parval(ipcdad,ir)
      dcm_adot = parval(ipcmad,ir)

c---- unit-freestream vector Vbar from alfa,beta
      alfa = parval(ipalfa,ir)*dtr
      beta = parval(ipbeta,ir)*dtr
      call ab2vbar
      sina = sin(alfa)
      cosa = cos(alfa)

c---- Euler angles and T,K matrices
      ang(1) = phi*dtr
      ang(2) = the*dtr
      ang(3) = psi*dtr
      call rotens3(ang,tt,tt_ang)
      call rateki3(ang,rt,rt_ang)

c-----------------------------------------------------------------
c---- normalized flight dynamics velocities, rotation rates,
c-      from JVL variables
c-    v(1:3) = (  u/Vref ,  v/Vref ,  w/Vref  )          
c-    w(1:3) = ( pb/2Vref, qc/2Vref, rb/2Vref )          
      v(1) = -vbar(1)
      v(2) = -vbar(2)
      v(3) = -vbar(3)
      w(1) =  wbar(1)*0.5*bref
      w(2) =  wbar(2)*0.5*cref
      w(3) =  wbar(3)*0.5*bref

c---- force and moment coefficients, 
c-     and stability and control derivatives in geometry xyz axes
      drm(1) = 0.
      drm(2) = 0.
      drm(3) = 0.
      call gcoeffs(ir,drm)

c---- dimensional reference quantities
      srefd = sref*unitl**2
      brefd = bref*unitl
      crefd = cref*unitl

      xyzref(1) = xcg
      xyzref(2) = ycg
      xyzref(3) = zcg

      lrefd(1) = brefd
      lrefd(2) = crefd
      lrefd(3) = brefd

c---- dimensional u,v,w velocities, p,q,r rotation rates
      do k = 1, 3
        vd(k) = v(k)*vrefd
        wd(k) = w(k)*2.0*vrefd/lrefd(k)
      enddo

c---- stability derivatives wrt dimensional velocities, rotation rates
      do k = 1, 3
        do n = 1, 3
          cf_vd(k,n) = cf_v(k,n)/vrefd
          cm_vd(k,n) = cm_v(k,n)/vrefd

          cf_wd(k,n) = cf_w(k,n)*0.5*lrefd(n)/vrefd
          cm_wd(k,n) = cm_w(k,n)*0.5*lrefd(n)/vrefd
        enddo
      enddo

c---- dimensional dynamic pressure*Sref
      qsd  = 0.5*rho*vrefd**2 * srefd
      qsd_vd(1) = 0.
      qsd_vd(2) = 0.
      qsd_vd(3) = 0.

c---- mass and inertia tensors (real + apparent)
      do k = 1, 3
        mamat(k,1) = amass(k,1)*rho
        mamat(k,2) = amass(k,2)*rho
        mamat(k,3) = amass(k,3)*rho
        mamat(k,k) = mamat(k,k) + rmass

        rimat(k,1) = riner(k,1) + ainer(k,1)*rho
        rimat(k,2) = riner(k,2) + ainer(k,2)*rho
        rimat(k,3) = riner(k,3) + ainer(k,3)*rho
      enddo

c---- overall mass+inertia matrix (without aero rate terms)
      do k = 1, 3
        do l = 1, 6
          mmat(k,l) = 0.
          mmat(k+3,l) = 0.
        enddo

        mmat(k,1) = mamat(k,1)
        mmat(k,2) = mamat(k,2)
        mmat(k,3) = mamat(k,3)

        mmat(k+3,4) = rimat(k,1)
        mmat(k+3,5) = rimat(k,2)
        mmat(k+3,6) = rimat(k,3)
      enddo

c---- additional apparent inertia from aero rates
      dcx_wdot = -(dcd_adot*cosa - dcl_adot*sina)
      dcz_wdot = -(dcd_adot*sina + dcl_adot*cosa)
      dcm_wdot = -dcm_adot
      mmat(1,3) = mmat(1,3) - dcx_wdot * qsd*0.5*crefd/vrefd**2
      mmat(3,3) = mmat(3,3) - dcz_wdot * qsd*0.5*crefd/vrefd**2
      mmat(5,3) = mmat(5,3) - dcm_wdot * qsd*0.5*crefd/vrefd**2 * crefd

c---- invert total mass+inertia tensor
      call m6inv(mmat,minv)

c---- linear and angular momentum vectors in JVL axes
      do k = 1, 3
c                        -     -   
c------ linear momentum  P = m V 
        p(k) = rmass*vd(k)
        p_vd(k,1) = 0.
        p_vd(k,2) = 0.
        p_vd(k,3) = 0.
        p_vd(k,k) = rmass
c                         -   = -   -
c------ angular momentum  H = I W + h
        h(k) = riner(k,1)*wd(1)
     &       + riner(k,2)*wd(2)
     &       + riner(k,3)*wd(3)  +  hrotor(k)
        h_wd(k,1) = riner(k,1)
        h_wd(k,2) = riner(k,2)
        h_wd(k,3) = riner(k,3)

c------ gravity force
        mg(k)       = -rmass*gee*tt(3,k)
        mg_ang(k,1) = -rmass*gee*tt_ang(3,k,1)
        mg_ang(k,2) = -rmass*gee*tt_ang(3,k,2)
        mg_ang(k,3) = -rmass*gee*tt_ang(3,k,3)
      enddo

c---- momentum rates from rotation, in JVL axes
      do k = 1, 3
        ic = icrs(k)
        jc = jcrs(k)

c------ linear momentum rate  W x P
        wxp(k)      = wd(ic)*p(jc)      - wd(jc)*p(ic)
        wxp_vd(k,1) = wd(ic)*p_vd(jc,1) - wd(jc)*p_vd(ic,1)
        wxp_vd(k,2) = wd(ic)*p_vd(jc,2) - wd(jc)*p_vd(ic,2)
        wxp_vd(k,3) = wd(ic)*p_vd(jc,3) - wd(jc)*p_vd(ic,3)
        wxp_wd(k,ic) =       p(jc)
        wxp_wd(k,jc) =                  -        p(ic)
        wxp_wd(k,k)  = 0.

c------ angular momentum rate  W x H
        wxh(k)      = wd(ic)*h(jc)      - wd(jc)*h(ic)
        wxh_vd(k,1) = 0.
        wxh_vd(k,2) = 0.
        wxh_vd(k,3) = 0.
        wxh_wd(k,1) = wd(ic)*h_wd(jc,1) - wd(jc)*h_wd(ic,1)
        wxh_wd(k,2) = wd(ic)*h_wd(jc,2) - wd(jc)*h_wd(ic,2)
        wxh_wd(k,3) = wd(ic)*h_wd(jc,3) - wd(jc)*h_wd(ic,3)
        wxh_wd(k,ic) = wxh_wd(k,ic) + h(jc)
        wxh_wd(k,jc) = wxh_wd(k,jc) - h(ic)
      enddo

c---- F, M
      do k = 1, 3
        fd(k) = cf(k)*qsd
        md(k) = cm(k)*qsd*lrefd(k)
        do n = 1, 3
          fd_vd(k,n) =  cf_vd(k,n)*qsd + cf(k)*qsd_vd(n)
          md_vd(k,n) = (cm_vd(k,n)*qsd + cm(k)*qsd_vd(n))*lrefd(k)
        enddo
        do n = 1, 3
          fd_wd(k,n) = cf_wd(k,n)*qsd
          md_wd(k,n) = cm_wd(k,n)*qsd*lrefd(k)
        enddo
        do n = 1, ncontrol
          fd_d(k,n) = cf_d(k,n)*qsd
          md_d(k,n) = cm_d(k,n)*qsd*lrefd(k)
        enddo
        do n = 1, nvarjet
          fd_j(k,n) = cf_j(k,n)*qsd
          md_j(k,n) = cm_j(k,n)*qsd*lrefd(k)
        enddo
      enddo

      do k = 1, 3
c------ m^-1 F
        mif(k) = minv(k,1)*fd(1)
     &         + minv(k,2)*fd(2)
     &         + minv(k,3)*fd(3)
     &         + minv(k,4)*md(1)
     &         + minv(k,5)*md(2)
     &         + minv(k,6)*md(3)

c------ i^-1 M
        rim(k) = minv(k+3,1)*fd(1)
     &         + minv(k+3,2)*fd(2)
     &         + minv(k+3,3)*fd(3)
     &         + minv(k+3,4)*md(1)
     &         + minv(k+3,5)*md(2)
     &         + minv(k+3,6)*md(3)

c------ m^-1 (WxP)
        prf(k) = minv(k,1)*wxp(1)
     &         + minv(k,2)*wxp(2)
     &         + minv(k,3)*wxp(3)
     &         + minv(k,4)*wxh(1)
     &         + minv(k,5)*wxh(2)
     &         + minv(k,6)*wxh(3)

c------ i^-1 (WxH)
        prm(k) = minv(k+3,1)*wxp(1)
     &         + minv(k+3,2)*wxp(2)
     &         + minv(k+3,3)*wxp(3)
     &         + minv(k+3,4)*wxh(1)
     &         + minv(k+3,5)*wxh(2)
     &         + minv(k+3,6)*wxh(3)

c------ m^-1 (mg)
        mgf(k) = minv(k,1)*mg(1)
     &         + minv(k,2)*mg(2)
     &         + minv(k,3)*mg(3)

c------ i^-1 (mg)
        mgm(k) = minv(k+3,1)*mg(1)
     &         + minv(k+3,2)*mg(2)
     &         + minv(k+3,3)*mg(3)

        do l = 1, 3
          mgf_ang(k,l) = minv(k,1)*mg_ang(1,l)
     &                 + minv(k,2)*mg_ang(2,l)
     &                 + minv(k,3)*mg_ang(3,l)
          mgm_ang(k,l) = minv(k+3,1)*mg_ang(1,l)
     &                 + minv(k+3,2)*mg_ang(2,l)
     &                 + minv(k+3,3)*mg_ang(3,l)
        enddo
        do n = 1, 3
          mif_vd(k,n) = 
     &           minv(k,1)*fd_vd(1,n)
     &         + minv(k,2)*fd_vd(2,n)
     &         + minv(k,3)*fd_vd(3,n)
     &         + minv(k,4)*md_vd(1,n)
     &         + minv(k,5)*md_vd(2,n)
     &         + minv(k,6)*md_vd(3,n)
          rim_vd(k,n) = 
     &           minv(k+3,1)*fd_vd(1,n)
     &         + minv(k+3,2)*fd_vd(2,n)
     &         + minv(k+3,3)*fd_vd(3,n)
     &         + minv(k+3,4)*md_vd(1,n)
     &         + minv(k+3,5)*md_vd(2,n)
     &         + minv(k+3,6)*md_vd(3,n)
          prf_vd(k,n) =
     &           minv(k,1)*wxp_vd(1,n)
     &         + minv(k,2)*wxp_vd(2,n)
     &         + minv(k,3)*wxp_vd(3,n)
     &         + minv(k,4)*wxh_vd(1,n)
     &         + minv(k,5)*wxh_vd(2,n)
     &         + minv(k,6)*wxh_vd(3,n)
          prm_vd(k,n) =
     &           minv(k+3,1)*wxp_vd(1,n)
     &         + minv(k+3,2)*wxp_vd(2,n)
     &         + minv(k+3,3)*wxp_vd(3,n)
     &         + minv(k+3,4)*wxh_vd(1,n)
     &         + minv(k+3,5)*wxh_vd(2,n)
     &         + minv(k+3,6)*wxh_vd(3,n)
        enddo

        do n = 1, 3
          mif_wd(k,n) = 
     &           minv(k,1)*fd_wd(1,n)
     &         + minv(k,2)*fd_wd(2,n)
     &         + minv(k,3)*fd_wd(3,n)
     &         + minv(k,4)*md_wd(1,n)
     &         + minv(k,5)*md_wd(2,n)
     &         + minv(k,6)*md_wd(3,n)
          rim_wd(k,n) = 
     &           minv(k+3,1)*fd_wd(1,n)
     &         + minv(k+3,2)*fd_wd(2,n)
     &         + minv(k+3,3)*fd_wd(3,n)
     &         + minv(k+3,4)*md_wd(1,n)
     &         + minv(k+3,5)*md_wd(2,n)
     &         + minv(k+3,6)*md_wd(3,n)
          prf_wd(k,n) =
     &           minv(k,1)*wxp_wd(1,n)
     &         + minv(k,2)*wxp_wd(2,n)
     &         + minv(k,3)*wxp_wd(3,n)
     &         + minv(k,4)*wxh_wd(1,n)
     &         + minv(k,5)*wxh_wd(2,n)
     &         + minv(k,6)*wxh_wd(3,n)
          prm_wd(k,n) =
     &           minv(k+3,1)*wxp_wd(1,n)
     &         + minv(k+3,2)*wxp_wd(2,n)
     &         + minv(k+3,3)*wxp_wd(3,n)
     &         + minv(k+3,4)*wxh_wd(1,n)
     &         + minv(k+3,5)*wxh_wd(2,n)
     &         + minv(k+3,6)*wxh_wd(3,n)
        enddo

        do n = 1, ncontrol
          mif_d(k,n) = 
     &           minv(k,1)*fd_d(1,n)
     &         + minv(k,2)*fd_d(2,n)
     &         + minv(k,3)*fd_d(3,n)
     &         + minv(k,4)*md_d(1,n)
     &         + minv(k,5)*md_d(2,n)
     &         + minv(k,6)*md_d(3,n)
          rim_d(k,n) = 
     &           minv(k+3,1)*fd_d(1,n)
     &         + minv(k+3,2)*fd_d(2,n)
     &         + minv(k+3,3)*fd_d(3,n)
     &         + minv(k+3,4)*md_d(1,n)
     &         + minv(k+3,5)*md_d(2,n)
     &         + minv(k+3,6)*md_d(3,n)
        enddo
        do n = 1, nvarjet
          mif_j(k,n) = 
     &           minv(k,1)*fd_j(1,n)
     &         + minv(k,2)*fd_j(2,n)
     &         + minv(k,3)*fd_j(3,n)
     &         + minv(k,4)*md_j(1,n)
     &         + minv(k,5)*md_j(2,n)
     &         + minv(k,6)*md_j(3,n)
          rim_j(k,n) = 
     &           minv(k+3,1)*fd_j(1,n)
     &         + minv(k+3,2)*fd_j(2,n)
     &         + minv(k+3,3)*fd_j(3,n)
     &         + minv(k+3,4)*md_j(1,n)
     &         + minv(k+3,5)*md_j(2,n)
     &         + minv(k+3,6)*md_j(3,n)
        enddo
      enddo

      nsys = 12
      do ieq = 1, nsys
        do je = 1, nsys
          asys(ieq,je) = 0.
        enddo
        do n = 1, ncontrol+nvarjet
          bsys(ieq,n) = 0.
        enddo
      enddo
      
c---- du/dt
      ieq = jeu
      k = 1
      rsys(ieq)      = mif(k)      - prf(k)     + mgf(k)
      asys(ieq,jeu)  = mif_vd(k,1) - prf_vd(k,1)
      asys(ieq,jev)  = mif_vd(k,2) - prf_vd(k,2)
      asys(ieq,jew)  = mif_vd(k,3) - prf_vd(k,3)
      asys(ieq,jep)  = mif_wd(k,1) - prf_wd(k,1)
      asys(ieq,jeq)  = mif_wd(k,2) - prf_wd(k,2)
      asys(ieq,jer)  = mif_wd(k,3) - prf_wd(k,3)
      asys(ieq,jeph) =                          + mgf_ang(k,1)
      asys(ieq,jeth) =                          + mgf_ang(k,2)
      asys(ieq,jeps) =                          + mgf_ang(k,3)
      do n = 1, ncontrol
        bsys(ieq,n)  = mif_d(k,n)
      enddo
      do n = 1, nvarjet
        nvj = ncontrol+n
        bsys(ieq,nvj)  = mif_j(k,n)
      enddo

c---- dv/dt
      ieq = jev
      k = 2
      rsys(ieq)      = mif(k)      - prf(k)     + mgf(k)
      asys(ieq,jeu)  = mif_vd(k,1) - prf_vd(k,1)
      asys(ieq,jev)  = mif_vd(k,2) - prf_vd(k,2)
      asys(ieq,jew)  = mif_vd(k,3) - prf_vd(k,3)
      asys(ieq,jep)  = mif_wd(k,1) - prf_wd(k,1)
      asys(ieq,jeq)  = mif_wd(k,2) - prf_wd(k,2)
      asys(ieq,jer)  = mif_wd(k,3) - prf_wd(k,3)
      asys(ieq,jeph) =                          + mgf_ang(k,1)
      asys(ieq,jeth) =                          + mgf_ang(k,2)
      asys(ieq,jeps) =                          + mgf_ang(k,3)
      do n = 1, ncontrol
        bsys(ieq,n)  = mif_d(k,n)
      enddo
      do n = 1, nvarjet
        nvj = ncontrol+n
        bsys(ieq,nvj)  = mif_j(k,n)
      enddo

c---- dw/dt
      ieq = jew
      k = 3
      rsys(ieq)      = mif(k)      - prf(k)     + mgf(k)
      asys(ieq,jeu)  = mif_vd(k,1) - prf_vd(k,1)
      asys(ieq,jev)  = mif_vd(k,2) - prf_vd(k,2)
      asys(ieq,jew)  = mif_vd(k,3) - prf_vd(k,3)
      asys(ieq,jep)  = mif_wd(k,1) - prf_wd(k,1)
      asys(ieq,jeq)  = mif_wd(k,2) - prf_wd(k,2)
      asys(ieq,jer)  = mif_wd(k,3) - prf_wd(k,3)
      asys(ieq,jeph) =                          + mgf_ang(k,1)
      asys(ieq,jeth) =                          + mgf_ang(k,2)
      asys(ieq,jeps) =                          + mgf_ang(k,3)
      do n = 1, ncontrol
        bsys(ieq,n)  = mif_d(k,n)
      enddo
      do n = 1, nvarjet
        nvj = ncontrol+n
        bsys(ieq,nvj)  = mif_j(k,n)
      enddo

c---- dp/dt
      ieq = jep
      k = 1
      rsys(ieq)      = rim(k)      - prm(k)     + mgm(k)
      asys(ieq,jeu)  = rim_vd(k,1) - prm_vd(k,1)
      asys(ieq,jev)  = rim_vd(k,2) - prm_vd(k,2)
      asys(ieq,jew)  = rim_vd(k,3) - prm_vd(k,3)
      asys(ieq,jep)  = rim_wd(k,1) - prm_wd(k,1)
      asys(ieq,jeq)  = rim_wd(k,2) - prm_wd(k,2)
      asys(ieq,jer)  = rim_wd(k,3) - prm_wd(k,3)
      asys(ieq,jeph) =                          + mgm_ang(k,1)
      asys(ieq,jeth) =                          + mgm_ang(k,2)
      asys(ieq,jeps) =                          + mgm_ang(k,3)
      do n = 1, ncontrol
        bsys(ieq,n)  = rim_d(k,n)
      enddo
      do n = 1, nvarjet
        nvj = ncontrol+n
        bsys(ieq,nvj) = rim_j(k,n)
      enddo

c---- dq/dt
      ieq = jeq
      k = 2
      rsys(ieq)      = rim(k)      - prm(k)     + mgm(k)
      asys(ieq,jeu)  = rim_vd(k,1) - prm_vd(k,1)
      asys(ieq,jev)  = rim_vd(k,2) - prm_vd(k,2)
      asys(ieq,jew)  = rim_vd(k,3) - prm_vd(k,3)
      asys(ieq,jep)  = rim_wd(k,1) - prm_wd(k,1)
      asys(ieq,jeq)  = rim_wd(k,2) - prm_wd(k,2)
      asys(ieq,jer)  = rim_wd(k,3) - prm_wd(k,3)
      asys(ieq,jeph) =                          + mgm_ang(k,1)
      asys(ieq,jeth) =                          + mgm_ang(k,2)
      asys(ieq,jeps) =                          + mgm_ang(k,3)
      do n = 1, ncontrol
        bsys(ieq,n)  = rim_d(k,n)
      enddo
      do n = 1, nvarjet
        nvj = ncontrol+n
        bsys(ieq,nvj) = rim_j(k,n)
      enddo

c---- dr/dt
      ieq = jer
      k = 3
      rsys(ieq)      = rim(k)      - prm(k)     + mgm(k)
      asys(ieq,jeu)  = rim_vd(k,1) - prm_vd(k,1)
      asys(ieq,jev)  = rim_vd(k,2) - prm_vd(k,2)
      asys(ieq,jew)  = rim_vd(k,3) - prm_vd(k,3)
      asys(ieq,jep)  = rim_wd(k,1) - prm_wd(k,1)
      asys(ieq,jeq)  = rim_wd(k,2) - prm_wd(k,2)
      asys(ieq,jer)  = rim_wd(k,3) - prm_wd(k,3)
      asys(ieq,jeph) =                          + mgm_ang(k,1)
      asys(ieq,jeth) =                          + mgm_ang(k,2)
      asys(ieq,jeps) =                          + mgm_ang(k,3)
      do n = 1, ncontrol
        bsys(ieq,n)  = rim_d(k,n)
      enddo
      do n = 1, nvarjet
        nvj = ncontrol+n
        bsys(ieq,nvj) = rim_j(k,n)
      enddo

c---- dphi/dt
      ieq = jeph
      k = 1
      rsys(ieq)      = rt(k,1)*wd(1)
     &               + rt(k,2)*wd(2)
     &               + rt(k,3)*wd(3)
      asys(ieq,jep)  = rt(k,1)
      asys(ieq,jeq)  = rt(k,2)
      asys(ieq,jer)  = rt(k,3)
      asys(ieq,jeph) = rt_ang(k,1,1)*wd(1)
     &               + rt_ang(k,2,1)*wd(2)
     &               + rt_ang(k,3,1)*wd(3)
      asys(ieq,jeth) = rt_ang(k,1,2)*wd(1)
     &               + rt_ang(k,2,2)*wd(2)
     &               + rt_ang(k,3,2)*wd(3)
      asys(ieq,jeps) = rt_ang(k,1,3)*wd(1)
     &               + rt_ang(k,2,3)*wd(2)
     &               + rt_ang(k,3,3)*wd(3)

c---- dtheta/dt
      ieq = jeth
      k = 2
      rsys(ieq)      = rt(k,1)*wd(1)
     &               + rt(k,2)*wd(2)
     &               + rt(k,3)*wd(3)
      asys(ieq,jep)  = rt(k,1)
      asys(ieq,jeq)  = rt(k,2)
      asys(ieq,jer)  = rt(k,3)
      asys(ieq,jeph) = rt_ang(k,1,1)*wd(1)
     &               + rt_ang(k,2,1)*wd(2)
     &               + rt_ang(k,3,1)*wd(3)
      asys(ieq,jeth) = rt_ang(k,1,2)*wd(1)
     &               + rt_ang(k,2,2)*wd(2)
     &               + rt_ang(k,3,2)*wd(3)
      asys(ieq,jeps) = rt_ang(k,1,3)*wd(1)
     &               + rt_ang(k,2,3)*wd(2)
     &               + rt_ang(k,3,3)*wd(3)

c---- dpsi/dt
      ieq = jeps
      k = 3
      rsys(ieq)      = rt(k,1)*wd(1)
     &               + rt(k,2)*wd(2)
     &               + rt(k,3)*wd(3)
      asys(ieq,jep)  = rt(k,1)
      asys(ieq,jeq)  = rt(k,2)
      asys(ieq,jer)  = rt(k,3)
      asys(ieq,jeph) = rt_ang(k,1,1)*wd(1)
     &               + rt_ang(k,2,1)*wd(2)
     &               + rt_ang(k,3,1)*wd(3)
      asys(ieq,jeth) = rt_ang(k,1,2)*wd(1)
     &               + rt_ang(k,2,2)*wd(2)
     &               + rt_ang(k,3,2)*wd(3)
      asys(ieq,jeps) = rt_ang(k,1,3)*wd(1)
     &               + rt_ang(k,2,3)*wd(2)
     &               + rt_ang(k,3,3)*wd(3)

c---- dRx/dt
      ieq = jex
      k = 1
      rsys(ieq)      = tt(k,1)*vd(1)
     &               + tt(k,2)*vd(2)
     &               + tt(k,3)*vd(3)
      asys(ieq,jeu)  = tt(k,1)
      asys(ieq,jev)  = tt(k,2)
      asys(ieq,jew)  = tt(k,3)
      asys(ieq,jeph) = tt_ang(k,1,1)*vd(1)
     &               + tt_ang(k,2,1)*vd(2)
     &               + tt_ang(k,3,1)*vd(3)
      asys(ieq,jeth) = tt_ang(k,1,2)*vd(1)
     &               + tt_ang(k,2,2)*vd(2)
     &               + tt_ang(k,3,2)*vd(3)
      asys(ieq,jeps) = tt_ang(k,1,3)*vd(1)
     &               + tt_ang(k,2,3)*vd(2)
     &               + tt_ang(k,3,3)*vd(3)

c---- dRy/dt
      ieq = jey
      k = 2
      rsys(ieq)      = tt(k,1)*vd(1)
     &               + tt(k,2)*vd(2)
     &               + tt(k,3)*vd(3)
      asys(ieq,jeu)  = tt(k,1)
      asys(ieq,jev)  = tt(k,2)
      asys(ieq,jew)  = tt(k,3)
      asys(ieq,jeph) = tt_ang(k,1,1)*vd(1)
     &               + tt_ang(k,2,1)*vd(2)
     &               + tt_ang(k,3,1)*vd(3)
      asys(ieq,jeth) = tt_ang(k,1,2)*vd(1)
     &               + tt_ang(k,2,2)*vd(2)
     &               + tt_ang(k,3,2)*vd(3)
      asys(ieq,jeps) = tt_ang(k,1,3)*vd(1)
     &               + tt_ang(k,2,3)*vd(2)
     &               + tt_ang(k,3,3)*vd(3)

c---- dRz/dt
      ieq = jez
      k = 3
      rsys(ieq)      = tt(k,1)*vd(1)
     &               + tt(k,2)*vd(2)
     &               + tt(k,3)*vd(3)
      asys(ieq,jeu)  = tt(k,1)
      asys(ieq,jev)  = tt(k,2)
      asys(ieq,jew)  = tt(k,3)
      asys(ieq,jeph) = tt_ang(k,1,1)*vd(1)
     &               + tt_ang(k,2,1)*vd(2)
     &               + tt_ang(k,3,1)*vd(3)
      asys(ieq,jeth) = tt_ang(k,1,2)*vd(1)
     &               + tt_ang(k,2,2)*vd(2)
     &               + tt_ang(k,3,2)*vd(3)
      asys(ieq,jeps) = tt_ang(k,1,3)*vd(1)
     &               + tt_ang(k,2,3)*vd(2)
     &               + tt_ang(k,3,3)*vd(3)

      return
      end ! sysmat


      subroutine appmat(ir,asys,bsys,rsys,nsys)
c------------------------------------------------------------------
c     Computes system matrices for run case ir.
c     Current forces and derivatives are assumed to be correct.
c------------------------------------------------------------------
      include 'jvl.inc'
      real*8 asys(jemax,jemax),
     &       bsys(jemax,ndmax+njmax),rsys(jemax)

      logical lterr

      real riner(3,3),
     &     mamat(3,3),
     &     mainv(3,3),
     &     rimat(3,3),
     &     riinv(3,3),
     &     p(3), p_u(3,numax),
     &     h(3), h_u(3,numax),
     &     wxp(3), wxp_u(3,numax),
     &     wxh(3), wxh_u(3,numax),
     &     mif(3), mif_u(3,numax), mif_d(3,ndmax),
     &     rim(3), rim_u(3,numax), rim_d(3,ndmax),
     &     prf(3), prf_u(3,numax),
     &     prm(3), prm_u(3,numax)
      real cft(3), cft_u(3,numax), cft_d(3,ndmax)
      real cmt(3), cmt_u(3,numax), cmt_d(3,ndmax)

      real ang(3), tt(3,3), tt_ang(3,3,3), rt(3,3), rt_ang(3,3,3)

      integer icrs(3), jcrs(3)
      data icrs / 2, 3, 1 / , jcrs / 3, 1, 2 /

c---- standard force and moment coefficients
      cft(1) = ft(1) * 2.0/ sref
      cft(2) = ft(2) * 2.0/ sref
      cft(3) = ft(3) * 2.0/ sref
      cmt(1) = mt(1) * 2.0/(sref*bref)
      cmt(2) = mt(2) * 2.0/(sref*cref)
      cmt(3) = mt(3) * 2.0/(sref*bref)
      do n = 1, numax
        cft_u(1,n) = ft_u(1,n) * 2.0/ sref
        cft_u(2,n) = ft_u(2,n) * 2.0/ sref
        cft_u(3,n) = ft_u(3,n) * 2.0/ sref
        cmt_u(1,n) = mt_u(1,n) * 2.0/(sref*bref)
        cmt_u(2,n) = mt_u(2,n) * 2.0/(sref*cref)
        cmt_u(3,n) = mt_u(3,n) * 2.0/(sref*bref)
      enddo
      do n = 1, ncontrol
        cft_d(1,n) = ft_d(1,n) * 2.0/ sref
        cft_d(2,n) = ft_d(2,n) * 2.0/ sref
        cft_d(3,n) = ft_d(3,n) * 2.0/ sref
        cmt_d(1,n) = mt_d(1,n) * 2.0/(sref*bref)
        cmt_d(2,n) = mt_d(2,n) * 2.0/(sref*cref)
        cmt_d(3,n) = mt_d(3,n) * 2.0/(sref*bref)
      enddo


      vrefd = parval(ipvee,ir)
      gee = parval(ipgee,ir)
      rho = parval(iprho,ir)
      phi = parval(ipphi,ir)
      the = parval(ipthe,ir)
      psi = parval(ippsi,ir)
      xcg = parval(ipxcg,ir)
      ycg = parval(ipycg,ir)
      zcg = parval(ipzcg,ir)
      rmass = parval(ipmass,ir)
      riner(1,1) = parval(ipixx,ir)
      riner(2,2) = parval(ipiyy,ir)
      riner(3,3) = parval(ipizz,ir)
      riner(1,2) = parval(ipixy,ir)
      riner(2,3) = parval(ipiyz,ir)
      riner(3,1) = parval(ipizx,ir)
      riner(2,1) = parval(ipixy,ir)
      riner(3,2) = parval(ipiyz,ir)
      riner(1,3) = parval(ipizx,ir)

      hrotor(1) = parval(iphx,ir)
      hrotor(2) = parval(iphy,ir)
      hrotor(3) = parval(iphz,ir)

      dcl_u = parval(ipclu,ir)
      dcm_u = parval(ipcmu,ir)
      dcl_a = parval(ipcla,ir)
      dcm_a = parval(ipcma,ir)

      lterr = .false.
      if(vrefd .le. 0.0) then
       write(*,*)
       write(*,*) '** Zero velocity.  Specify with run file or M menu'
       lterr = .true.
      endif
      if(rmass .le. 0.0) then
       write(*,*)
       write(*,*) '** Zero mass.  Specify with mass file or M menu'
       lterr = .true.
      endif
      if(riner(1,1) .le. 0.0) then
       write(*,*)
       write(*,*) '** Zero Ixx.  Specify with mass file or M menu'
       lterr = .true.
      endif
      if(riner(2,2) .le. 0.0) then
       write(*,*)
       write(*,*) '** Zero Iyy.  Specify with mass file or M menu'
       lterr = .true.
      endif
      if(riner(3,3) .le. 0.0) then
       write(*,*)
       write(*,*) '** Zero Izz.  Specify with mass file or M menu'
       lterr = .true.
      endif

      if(lterr) then
       write(*,*)
       write(*,*) 'Eigenmodes not computed for run case', ir
       return
      endif

      srefd = sref*unitl**2
      brefd = bref*unitl
      crefd = cref*unitl
c      srefd = unitl**2
c      brefd = unitl
c      crefd = unitl

      xyzref(1) = xcg
      xyzref(2) = ycg
      xyzref(3) = zcg

      qsd  = 0.5*rho*vrefd**2 * srefd

      qsdc = qsd*crefd
      qsdb = qsd*brefd

      rinxx = riner(1,1)
      rinyy = riner(2,2)
      rinzz = riner(3,3)

c
      rotd = vrefd/unitl

c
      ang(1) = phi*dtr
      ang(2) = the*dtr
      ang(3) = psi*dtr
      call rotens3(ang,tt,tt_ang)
      call rateki3(ang,rt,rt_ang)

      nsys = 12
      do ieq = 1, nsys
        do je = 1, nsys
          asys(ieq,je) = 0.
        enddo
        do n = 1, ncontrol+nvarjet
          bsys(ieq,n) = 0.
        enddo
      enddo

c---- x-acceleration
      ieq = jeu
      asys(ieq,jeu) = -cft_u(1,1)*qsd /rmass / vrefd
      asys(ieq,jew) = -cft_u(1,3)*qsd /rmass / vrefd
      asys(ieq,jeq) =  cft_u(1,5)*qsd /rmass / rotd  +  vbar(3)*vrefd
      asys(ieq,jeth)=  gee
      do n = 1, ncontrol
        bsys(ieq,n) =  cft_d(1,n)*qsd /rmass
      enddo

c---- z-acceleration
      ieq = jew
      asys(ieq,jeu) = -cft_u(3,1)*qsd /rmass / vrefd
      asys(ieq,jew) = -cft_u(3,3)*qsd /rmass / vrefd
      asys(ieq,jeq) =  cft_u(3,5)*qsd /rmass / rotd  -  vbar(1)*vrefd
      do n = 1, ncontrol
        bsys(ieq,n) =  cft_d(3,n)*qsd /rmass
      enddo

c---- y-ang.accel.
      ieq = jeq
      asys(ieq,jeu) = -cmt_u(2,1)*qsdc/rinyy / vrefd
      asys(ieq,jew) = -cmt_u(2,3)*qsdc/rinyy / vrefd
      asys(ieq,jeq) =  cmt_u(2,5)*qsdc/rinyy / rotd
      do n = 1, ncontrol
        bsys(ieq,n) =  cmt_d(2,n)*qsdc/rinyy
      enddo

c---- theta rate
      ieq = jeth
      asys(ieq,jeq) =  1.0

c

c---- y-acceleration
      ieq = jev
      asys(ieq,jev) = -cft_u(2,2)*qsd /rmass / vrefd
      asys(ieq,jep) =  cft_u(2,4)*qsd /rmass / rotd  -  vbar(3)*vrefd
      asys(ieq,jer) =  cft_u(2,6)*qsd /rmass / rotd  +  vbar(1)*vrefd
      asys(ieq,jeph)=  gee
      do n = 1, ncontrol
        bsys(ieq,n) =  cft_d(2,n)*qsd /rmass
      enddo

c---- x-ang.accel.
      ieq = jep
      asys(ieq,jev) = -cmt_u(1,2)*qsdb/rinxx / vrefd
      asys(ieq,jep) =  cmt_u(1,4)*qsdb/rinxx / rotd
      asys(ieq,jer) =  cmt_u(1,6)*qsdb/rinxx / rotd
      do n = 1, ncontrol
        bsys(ieq,n) =  cmt_d(1,n)*qsdb/rinxx
      enddo

c---- z-ang.accel.
      ieq = jer
      asys(ieq,jev) = -cmt_u(3,2)*qsdb/rinzz / vrefd
      asys(ieq,jep) =  cmt_u(3,4)*qsdb/rinzz / rotd
      asys(ieq,jer) =  cmt_u(3,6)*qsdb/rinzz / rotd
      do n = 1, ncontrol
        bsys(ieq,n) =  cmt_d(3,n)*qsdb/rinzz
      enddo

c---- phi rate
      ieq = jeph
      asys(ieq,jep) =  -1.0



c---- psi rate
      ieq = jeps
      k = 3
      rsys(ieq)      = (  rt(k,1)*wbar(1)
     &                  + rt(k,2)*wbar(2)
     &                  + rt(k,3)*wbar(3) )*rotd
      asys(ieq,jep)  =    rt(k,1)
      asys(ieq,jeq)  =    rt(k,2)
      asys(ieq,jer)  =    rt(k,3)
      asys(ieq,jeph) = (  rt_ang(k,1,1)*wbar(1)
     &                  + rt_ang(k,2,1)*wbar(2)
     &                  + rt_ang(k,3,1)*wbar(3) )*rotd
      asys(ieq,jeth) = (  rt_ang(k,1,2)*wbar(1)
     &                  + rt_ang(k,2,2)*wbar(2)
     &                  + rt_ang(k,3,2)*wbar(3) )*rotd
      asys(ieq,jeps) = (  rt_ang(k,1,3)*wbar(1)
     &                  + rt_ang(k,2,3)*wbar(2)
     &                  + rt_ang(k,3,3)*wbar(3) )*rotd

c---- x-velocity
      ieq = jex
      k = 1
      rsys(ieq)     = -(  tt(k,1)*vbar(1)
     &                  + tt(k,2)*vbar(2)
     &                  + tt(k,3)*vbar(3) )*vrefd
      asys(ieq,jeu) =     tt(k,1)
      asys(ieq,jev) =     tt(k,2)
      asys(ieq,jew) =     tt(k,3)
      asys(ieq,jeph)= -(  tt_ang(k,1,1)*vbar(1)
     &                  + tt_ang(k,2,1)*vbar(2)
     &                  + tt_ang(k,3,1)*vbar(3) )*vrefd
      asys(ieq,jeth)= -(  tt_ang(k,1,2)*vbar(1)
     &                  + tt_ang(k,2,2)*vbar(2)
     &                  + tt_ang(k,3,2)*vbar(3) )*vrefd
      asys(ieq,jeps)= -(  tt_ang(k,1,3)*vbar(1)
     &                  + tt_ang(k,2,3)*vbar(2)
     &                  + tt_ang(k,3,3)*vbar(3) )*vrefd

c---- y-velocity
      ieq = jey
      k = 2
      rsys(ieq)     = -(  tt(k,1)*vbar(1)
     &                  + tt(k,2)*vbar(2)
     &                  + tt(k,3)*vbar(3) )*vrefd
      asys(ieq,jeu) =     tt(k,1)
      asys(ieq,jev) =     tt(k,2)
      asys(ieq,jew) =     tt(k,3)
      asys(ieq,jeph)= -(  tt_ang(k,1,1)*vbar(1)
     &                  + tt_ang(k,2,1)*vbar(2)
     &                  + tt_ang(k,3,1)*vbar(3) )*vrefd
      asys(ieq,jeth)= -(  tt_ang(k,1,2)*vbar(1)
     &                  + tt_ang(k,2,2)*vbar(2)
     &                  + tt_ang(k,3,2)*vbar(3) )*vrefd
      asys(ieq,jeps)= -(  tt_ang(k,1,3)*vbar(1)
     &                  + tt_ang(k,2,3)*vbar(2)
     &                  + tt_ang(k,3,3)*vbar(3) )*vrefd

c---- z-velocity
      ieq = jez
      k = 3
      rsys(ieq)     = -(  tt(k,1)*vbar(1)
     &                  + tt(k,2)*vbar(2)
     &                  + tt(k,3)*vbar(3) )*vrefd
      asys(ieq,jeu) =     tt(k,1)
      asys(ieq,jev) =     tt(k,2)
      asys(ieq,jew) =     tt(k,3)
      asys(ieq,jeph)= -(  tt_ang(k,1,1)*vbar(1)
     &                  + tt_ang(k,2,1)*vbar(2)
     &                  + tt_ang(k,3,1)*vbar(3) )*vrefd
      asys(ieq,jeth)= -(  tt_ang(k,1,2)*vbar(1)
     &                  + tt_ang(k,2,2)*vbar(2)
     &                  + tt_ang(k,3,2)*vbar(3) )*vrefd
      asys(ieq,jeps)= -(  tt_ang(k,1,3)*vbar(1)
     &                  + tt_ang(k,2,3)*vbar(2)
     &                  + tt_ang(k,3,3)*vbar(3) )*vrefd

      return
      end ! appmat


      subroutine syssho(lu,asys,bsys,rsys,nsys)
c------------------------------------------------------------------
c     Computes eigenvalues and eigenvectors for run case ir.
c     Current forces and derivatives are assumed to be correct.
c------------------------------------------------------------------
      include 'jvl.inc'
      real*8 asys(jemax,jemax),
     &       bsys(jemax,ndmax+njmax),rsys(jemax)

      write(lu,*)
      write(lu,1100)
     &'     u         w         q        the   ',
     &'     v         p         r        phi   ',
     &'     x         y         z        psi   ',
     & (dname(n), n=1, ncontrol),
     & (jname(n), n=1, nvarjet)
c      1234567890123456789012345678901234567890
 1100 format(1x,a,a,a,1x,'|',2x,50a12)

      do i = 1, nsys
        write(lu,1200) (asys(i,j), j=1, nsys),
     &                 (bsys(i,n), n=1, ncontrol+nvarjet)
 1200   format(1x,12f10.4,3x,50g12.4)
      enddo

      return
      end ! syssho


      subroutine eigsol(info,ir,etol,asys,nsys)
c------------------------------------------------------------------
c     Computes eigenvalues and eigenvectors for run case ir.
c     Current forces and derivatives are assumed to be correct.
c------------------------------------------------------------------
      include 'jvl.inc'
      real*8 asys(jemax,jemax)

      real*8 wr(jemax),wi(jemax),wvec(jemax,jemax)
      real*8 rwork(jemax)
      integer iwork(jemax)
      logical lterr

c---- call eispack's real/general eigenvalue routine
c     icalc = 0   ! get eigenvalues only
      icalc = 1   ! get eigenvalues and eigenvectors
      call rg(jemax,nsys,asys,wr,wi,icalc,wvec,iwork,rwork,ierr)

c-------------------------------------------------------
c---- store goodies

      vrefd = parval(ipvee,ir)
      brefd = bref*unitl
      etolsq = (etol*vrefd/brefd)**2

      keig = 0

      do 100 j=1, nsys
c------ don't store eigenvalue smaller than tolerance
        emagsq = wr(j)**2 + wi(j)**2
        if(emagsq .lt. etolsq) go to 100

        keig = keig + 1
        eval(keig,ir) = cmplx( wr(j), wi(j) )

        if    (wi(j) .eq. 0.0) then
c------- real eigenvalue... just store real eigenvector
         do i = 1, nsys
           evr = wvec(i,j)
           evi = 0.
           evec(i,keig,ir) = cmplx( evr , evi )
         enddo

        elseif(wi(j) .gt. 0.0) then
c------- positive imaginary part of eigenvalue... store complex eigenvector
         jp1 = min( j+1 , nsys )
         do i = 1, nsys
           evr = wvec(i,j)
           evi = wvec(i,jp1)
           evec(i,keig,ir) = cmplx( evr , evi )
         enddo

        elseif(wi(j) .lt. 0.0) then
c------- negative imaginary part of eigenvalue... store complex eigenvector
         jm = max( j-1 , 1 )
         do i = 1, nsys
           evr =  wvec(i,jm)
           evi = -wvec(i,j)
           evec(i,keig,ir) = cmplx( evr , evi )
         enddo
        endif
 100  continue

      neigen(ir) = keig

      return
      end ! eigsol



      subroutine runlst(lu,ire)
      include 'jvl.inc'
      character*1 rsym

      write(lu,1050)
 1050 format(1x,' ')
      write(lu,1100) '  Run  ',
     &    parnam(ipalfa),
     &    parnam(ipbeta), 
     &    parnam(ipcl  ), 
     &    parnam(ipcd0 ), 
     &    parnam(ipphi ), 
     &    parnam(ipvee ), 
     &    parnam(iprho ), 
     &    parnam(iprad ),
     &    parnam(ipfac ),
     &    parnam(ipxcg ),
     &    parnam(ipzcg ),
     &    parnam(ipmass),
     &    parnam(ipixx ),
     &    parnam(ipiyy ),
     &    parnam(ipizz ),
     &    parnam(iphx  ),
     &    parnam(iphy  ),
     &    parnam(iphz  ),
     &    parnam(ipdvjx)

      write(lu,1100) '       ',
     &    parunch(ipalfa),
     &    parunch(ipbeta), 
     &    parunch(ipcl  ), 
     &    parunch(ipcd0 ), 
     &    parunch(ipphi ), 
     &    parunch(ipvee ), 
     &    parunch(iprho ), 
     &    parunch(iprad ),
     &    parunch(ipfac ),
     &    parunch(ipxcg ),
     &    parunch(ipzcg ),
     &    parunch(ipmass),
     &    parunch(ipixx ),
     &    parunch(ipiyy ),
     &    parunch(ipizz ),
     &    parunch(iphx  ),
     &    parunch(iphy  ),
     &    parunch(iphz  ),
     &    parunch(ipdvjx)

 1100 format(1x,a, 25(1x,a9))

      do ir = 1, nrun
        if(ire.eq.ir .or. ire.eq.0) then
         rsym = '>'
        else
         rsym = ' '
        endif
        write(lu,1300) rsym, ir, 
     &    parval(ipalfa,ir),
     &    parval(ipbeta,ir), 
     &    parval(ipcl  ,ir), 
     &    parval(ipcd0 ,ir), 
     &    parval(ipphi ,ir), 
     &    parval(ipvee ,ir), 
     &    parval(iprho ,ir), 
     &    parval(iprad ,ir),
     &    parval(ipfac ,ir),
     &    parval(ipxcg ,ir),
     &    parval(ipzcg ,ir),
     &    parval(ipmass,ir),
     &    parval(ipixx ,ir),
     &    parval(ipiyy ,ir),
     &    parval(ipizz ,ir),
     &    parval(iphx  ,ir),
     &    parval(iphy  ,ir),
     &    parval(iphz  ,ir),
     &    parval(ipdvjx,ir)
 1300   format(1x,a, i3, 2x, 25g10.3)
      enddo

      return
      end ! runlst


