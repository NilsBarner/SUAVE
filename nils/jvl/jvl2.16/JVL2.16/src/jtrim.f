C***********************************************************************
C    Module:  jtrim.f
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

      subroutine trmset(comand,comarg,lmatch,ir)
      include 'jvl.inc'
      character*(*) comand, comarg
      logical lmatch

      character*2 cnum

      character*4 com
      character*80 carg, prompt, rtnew
      logical error, repkey

      integer iinp(10)
      real rinp(10)

c---- for trim command, first number in arg string should be part of command
      if(comand(1:1) .ne. 'C') then
       lmatch = .false.
       return
      endif

      cnum = comarg(1:2)

c
      if    (cnum.eq.'1 ') then
c----- set level or banked  horizontal flight constraints 
       ktrim = 1
       lmatch = .true.

      elseif(cnum.eq.'2 ') then
c----- set steady pitch rate (looping) flight constraints
       ktrim = 2
       lmatch = .true.

      else
       lmatch = .false.
       return

      endif

c---- start here with new run case
 100  continue

      if(ktrim.eq.1 .or.
     &   ktrim.eq.2      ) then
       if(parval(ipcl,ir) .eq. 0.0) then
c------ current case has cl=0 ... set it using cl constraint, if it's present
        do iv = 2, nvtot
          if(icon(iv,ir) .eq. iccl) then
           write(*,*)
           write(*,*)'       Setting trim CL from current CL constraint'
           parval(ipcl,ir) = conval(iccl,ir)
           go to 101
          endif
        enddo
       endif
 101   continue
      endif

      write(*,*)

c---- tag this run case with trim type (not converged yet)
      itrim(ir) = -ktrim

      if(parval(iprho,ir) .le. 0.0) parval(iprho,ir)  = rho0
      if(parval(ipgee,ir) .le. 0.0) parval(ipgee,ir)  = gee0
      if(parval(ipmass,ir).le. 0.0) parval(ipmass,ir) = rmass0

      ir1 = ir
      ir2 = ir

c----------------------------------------------------------------
c---- jump back here to calculate stuff depending on case
 5    continue

      crefd = cref*unitl
      brefd = bref*unitl
      srefd = sref*unitl**2

c- - - - - - - - - - - - - - - - - - - - - - - - 
 7500 format('     .. setting new ', a,' for run case', i3)

      if    (ktrim.eq.1) then
c------ level or banked horizontal flight:    recalculate v or cl, radius

        do 7 jr = ir1, ir2
          phi = parval(ipphi,jr)
          the = parval(ipthe,jr)
          cl  = parval(ipcl ,jr)
          cd0 = parval(ipcd0,jr)
          vee = parval(ipvee,jr)
          rad = parval(iprad,jr)
          rho = parval(iprho,jr)
          gee = parval(ipgee,jr)
          fac = parval(ipfac,jr)
          rmass = parval(ipmass,jr)

          sinp = sin(phi*dtr)
          cosp = cos(phi*dtr)

c-------- set velocity ?
          if(vee .le. 0.0 .and.
     &       cl  .gt. 0.0       ) then
           vee = sqrt( 2.0*rmass*gee / (rho*srefd*cl*cosp) )
           parval(ipvee,jr) = vee
           write(*,7500) 'velocity', jr
          endif

c-------- set cl ?
          if(cl  .le. 0.0 .and.
     &       vee .gt. 0.0       ) then
           cl  = 2.0*rmass*gee / (rho*srefd*vee**2*cosp)
           parval(ipcl,jr) = cl
           write(*,7500) 'CL', jr
          endif

          if(sinp .eq. 0.0) then
           rad = 0.
          else
           rad = vee**2 * cosp / (gee * sinp)
          endif
          parval(iprad,jr) = rad
          write(*,7500) 'turn radius', jr

          fac = 1.0/cosp
          parval(ipfac,jr) = fac
          write(*,7500) 'load factor', jr

          the = 0.
          parval(ipthe,jr) = the
cc        write(*,7500) 'elevation', jr

c-------- set up cl and rotation rate constraints
          if(rad .gt. 0.0) then
           whx = 0.
           why = sinp * crefd/(2.0*rad)
           whz = cosp * brefd/(2.0*rad)
          else
           whx = 0.
           why = 0.
           whz = 0.
          endif

          conval(iccl  ,jr) = cl 
          conval(icomgx,jr) = whx
          conval(icomgy,jr) = why
          conval(icomgz,jr) = whz

          icon(ivvelz,jr) = iccl
          icon(ivomgx,jr) = icomgx
          icon(ivomgy,jr) = icomgy
          icon(ivomgz,jr) = icomgz
 7      continue

c- - - - - - - - - - - - - - - - - - - - - - - - 
      elseif(ktrim.eq.2) then
c------ steady pitch rate (looping flight

        do 8 jr = ir1, ir2
          phi = parval(ipphi,jr)
          the = parval(ipthe,jr)
          cl  = parval(ipcl ,jr)
          cd0 = parval(ipcd0,jr)
          vee = parval(ipvee,jr)
          rad = parval(iprad,jr)
          rho = parval(iprho,jr)
          gee = parval(ipgee,jr)
          fac = parval(ipfac,jr)
          rmass = parval(ipmass,jr)

          sinp = sin(phi*dtr)
          cosp = cos(phi*dtr)

c-------- set radius ?
          if(rad .eq. 0.0 .and.
     &       cl  .gt. 0.0      ) then
           rad = rmass / (0.5*rho * srefd * cl)
           parval(iprad,jr) = rad
           write(*,7500) 'turn radius', jr
          endif

c-------- set cl ?
          if(rad .gt. 0.0 .and.
     &       cl  .eq. 0.0      ) then
           cl = rmass / (0.5*rho * srefd * rad)
           parval(ipcl,jr) = cl
           write(*,7500) 'cl', jr
          endif

c-------- set load factor ?
          if(fac .eq. 0.0 .and.
     &       cl  .gt. 0.0 .and. 
     &       vee .gt. 0.0 .and. 
     &       gee .gt. 0.0      ) then
           fac = 0.5*rho*vee**2 * srefd * cl / (rmass*gee)
           parval(ipfac,jr) = fac
           write(*,7500) 'load factor', jr
          endif

c-------- set velocity ?
          if(fac .gt. 0.0 .and.
     &       cl  .gt. 0.0 .and. 
     &       vee .eq. 0.0 .and. 
     &       gee .gt. 0.0      ) then
           vee = sqrt( fac * rmass*gee / (0.5*rho*srefd*cl) )
           parval(ipvee,jr) = vee
           write(*,7500) 'velocity', jr
          endif

          the = 0.
          parval(ipthe,jr) = the
cc        write(*,7500) 'elevation', jr

c-------- set up cl and rotation rate constraints
          if(rad .gt. 0.0) then
           whx = 0.
           why = crefd/(2.0*rad)
           whz = 0.
          else
           whx = 0.
           why = 0.
           whz = 0.
          endif

          conval(iccl  ,jr) = cl 
          conval(icomgx,jr) = whx
          conval(icomgy,jr) = why
          conval(icomgz,jr) = whz
c 
          icon(ivvelz,jr) = iccl
          icon(ivomgx,jr) = icomgx
          icon(ivomgy,jr) = icomgy
          icon(ivomgz,jr) = icomgz
 8      continue

      endif

      phi = parval(ipphi,ir)
      the = parval(ipthe,ir)
      cl  = parval(ipcl ,ir)
      cd0 = parval(ipcd0,ir)
      vee = parval(ipvee,ir)
      rad = parval(iprad,ir)
      rho = parval(iprho,ir)
      gee = parval(ipgee,ir)
      fac = parval(ipfac,ir)
      xcg = parval(ipxcg,ir)
      ycg = parval(ipycg,ir)
      zcg = parval(ipzcg,ir)
      rmass = parval(ipmass,ir)

      sinp = sin(phi*dtr)
      cosp = cos(phi*dtr)

 1000 format(a)
 2000 format(/'     Setup of trimmed run case ',a,':  ', a)
 2005 format( '     ', a
     &       /'     =================================================')
 2105 format(6x, a, g10.4, 2x, a)
 2110 format(6x, a, a )

c--------------------------------------------------------------------------
c----- jump back here just for menu
 10    continue
       call cfrac(ir,nrun,prompt,npr)
       write(*,2000) prompt(1:npr),rtitle(ir)

       if    (ktrim.eq.1) then
        write(*,2005) '(level or banked horizontal flight)'
        write(*,2105) 'B  bank angle = ', phi  , 'deg'
        write(*,2105) 'C  CL         = ', cl   , ' '
        write(*,2105) 'V  velocity   = ', vee  , unchv(1:nuv)
        write(*,2105) 'M  mass       = ', rmass, unchm(1:num)
        write(*,2105) 'D  air dens.  = ', rho  , unchd(1:nud)
        write(*,2105) 'G  grav.acc.  = ', gee  , uncha(1:nua)
        write(*,2105) '   turn rad.  = ', rad  , unchl(1:nul)
        write(*,2105) '   load fac.  = ', fac  , ' '
        write(*,2105) 'X  X_cg       = ', xcg  , 'Lunit'
        write(*,2105) 'Y  Y_cg       = ', ycg  , 'Lunit'
        write(*,2105) 'Z  Z_cg       = ', zcg  , 'Lunit'

       elseif(ktrim.eq.2) then
        write(*,2005) '(steady pitch rate - looping flight)'
        write(*,2105) 'C  CL        = ', cl   , ' '
        write(*,2105) 'V  velocity  = ', vee  , unchv(1:nuv) 
        write(*,2105) 'M  mass      = ', rmass, unchm(1:num)
        write(*,2105) 'D  air dens. = ', rho  , unchd(1:nud) 
        write(*,2105) 'G  grav.acc. = ', gee  , uncha(1:nua) 
        write(*,2105) 'R  turn rad. = ', rad  , unchl(1:nul) 
        write(*,2105) 'L  load fac. = ', fac  , ' '
        write(*,2105) 'X  X_cg      = ', xcg  , 'Lunit'
        write(*,2105) 'Y  Y_cg      = ', ycg  , 'Lunit'
        write(*,2105) 'Z  Z_cg      = ', zcg  , 'Lunit'

       endif

       call askc('     Enter parameter, value  (or  # - + N )',com,carg)

       if(com.eq.'    ') then
c------ just a return entered... go back
 2200   format(1x,'Setting ',a,
     &            ' constraint to specified ',a,' =',f12.5)
        write(*,*)
        write(*,2200) 'alpha     ', 'CL   ', cl
        write(*,2200) 'roll  rate', 'pb/2V', whx
        write(*,2200) 'pitch rate', 'qc/2V', why
        write(*,2200) 'yaw   rate', 'rb/2V', whz
        return
       endif

c------------------------------------------------------
c---- check for run case commands
      if(com.eq.'+   ') then
c----- add new case after current one

       if(nrun.eq.nrmax) then
        write(*,*)
        write(*,*) '* Run case array limit nrmax reached'
       else
        nrun = nrun + 1

        do jr = nrun, ir+1, -1
          call rcopy(jr,jr-1)
        enddo
        write(*,*) 'Initializing new run case from current one'
c        
        ir = ir + 1
       endif

       go to 100

      elseif(com.eq.'-   ') then
c----- delete current case

       if(nrun.le.1) then
        write(*,*)
        write(*,*) '* Cannot delete one remaining run case'
       else
        do jr = ir, nrun-1
          call rcopy(jr,jr+1)
        enddo
        nrun = nrun - 1
        irun = max( 1 , min( irun , nrun ) )
       endif

       go to 100

      endif

c------------------------------------------------------
c---- see if command is an integer (new run case selection)
      ninp = 1
      call getint(com,iinp,ninp,error)
      if(.not.error .and. ninp.ge.1 
     &  .and. com(1:1).ne.'T'
     &  .and. com(1:1).ne.'F' ) then
c----- command is an integer... new case index?
       if(iinp(1).lt.1 .or. iinp(1).gt.nrun) then
        write(*,*)
        write(*,*) '* Selected new run case is not defined'
        go to 10
       else
c------ valid new run case selected... go back to top of menu
        ir = iinp(1)
        go to 100
       endif
      endif


c------------------------------------------------------
c---- extract argument, if any
      ninp = 1
      call getflt(carg,rinp,ninp,error)

      if    (com(1:1) .eq. com(2:2)) then
       repkey = .true.
       ir1 = 1
       ir2 = nrun
       com(2:2) = ' '
      elseif(com(1:2) .eq. com(3:4)) then
       repkey = .true.
       ir1 = 1
       ir2 = nrun
       com(3:4) = '  '
      else
       repkey = .false.
       ir1 = ir
       ir2 = ir
      endif

c
c==============================================================================
c---- now decode regular parameter value commands

c------------------------------------
      if(com(1:1).eq.'B') then
        if    (ktrim.eq.1) then
 11      continue
         if(ninp.ge.1) then
          phi = rinp(1)
         else
          call askr('      Enter bank angle^',phi)
         endif

         if(phi.le.-90.0 .or. phi.ge.90.0) then
          write(*,*) '    * Must have  -90 < bank < +90'
          ninp = 0
          go to 11
         endif

         do jr = ir1, ir2
           parval(ipphi,jr) = phi

c--------- recalculate v, r, and n
           parval(ipvee,jr) = 0.
         enddo

        elseif(ktrim.eq.2) then
         write(*,*) 'bank angle not used for this trim case'
         go to 10

        endif
c-------------------------------------
      elseif(com(1:1).eq.'C' .or.
     &       com(1:2).eq.'CL'      ) then
 21     continue
        if(ninp.ge.1) then
         cl = rinp(1)
        else
         call askr('      Enter CL^',cl)
        endif

        if(cl .le. 0.0) then
         write(*,*) '    * Must have  CL > 0'
         ninp = 0
         go to 21
        endif

        do jr = ir1, ir2
          parval(ipcl,jr) = cl

          if    (ktrim.eq.1) then
c--------- go recalculate v and r
           parval(ipvee,jr) = 0.
           parval(iprad,jr) = 0.

          elseif(ktrim.eq.2) then
c--------- go recalculate r and n
           parval(iprad,jr) = 0.
           parval(ipfac,jr) = 0.

          endif
        enddo

c-------------------------------------
      elseif(com(1:1).eq.'V') then
 31     continue
        if(ninp.ge.1) then
         vee = rinp(1)
        else
         call askr('      Enter velocity^',vee)
        endif

        if(vee.le.0.0) then
         write(*,*) '    * Must have  V > 0'
         ninp = 0
         go to 31
        endif

        do jr = ir1, ir2
          parval(ipvee,jr) = vee

          if    (ktrim.eq.1) then
c--------- go recalculate cl
           parval(ipcl,jr) = 0.

          elseif(ktrim.eq.2) then
c--------- go recalculate n
           parval(ipfac,jr) = 0.

          endif
        enddo

c-------------------------------------
      elseif(com(1:1).eq.'M') then
 41     continue
        if(ninp.ge.1) then
         rmass = rinp(1)
        else
         call askr('      Enter mass^',rmass)
        endif

        if(rmass.le.0.0) then
         write(*,*) '    * Must have  m > 0'
         ninp = 0
         go to 41
        endif

        do jr = ir1, ir2
          parval(ipmass,jr) = rmass

          if    (ktrim.eq.1) then
c--------- go recalculate v
           parval(ipvee,jr) = 0.

          elseif(ktrim.eq.2) then
c--------- go recalculate r and n
           parval(iprad,jr) = 0.
           parval(ipfac,jr) = 0.

          endif
        enddo

c-------------------------------------
      elseif(com(1:1).eq.'D') then
 51     continue
        if(ninp.ge.1) then
         rho = rinp(1)
        else
         call askr('      Enter air density^',rho)
        endif

        if(rho.le.0.0) then
         write(*,*) '    * Must have  rho > 0'
         ninp = 0
         go to 51
        endif

        do jr = ir1, ir2
          parval(iprho,jr) = rho

          if    (ktrim.eq.1) then
c--------- go recalculate v
           parval(ipvee,jr) = 0.

          elseif(ktrim.eq.2) then
c--------- go recalculate r and n
           parval(iprad,jr) = 0.
           parval(ipfac,jr) = 0.

          endif
        enddo

c-------------------------------------
      elseif(com(1:1).eq.'G') then
 61     continue
        if(ninp.ge.1) then
         gee = rinp(1)
        else
         call askr('      Enter gravity^',gee)
        endif

        if(gee.le.0.0) then
         write(*,*) '    * Must have  g > 0'
         ninp = 0
         go to 61
        endif

        do jr = ir1, ir2
          parval(ipgee,jr) = gee

          if    (ktrim.eq.1) then
c--------- go recalculate v
           parval(ipvee,jr) = 0.

          elseif(ktrim.eq.2) then
c--------- go recalculate n
           parval(ipfac,jr) = 0.

          endif
        enddo

c-------------------------------------
      elseif(com(1:1).eq.'R') then
        if    (ktrim.eq.1) then
         write(*,*) 'Turn radius not specifiable for this trim case'
         go to 10

        elseif(ktrim.eq.2) then
 71      continue
         if(ninp.ge.1) then
          rad = rinp(1)
         else
          call askr('      Enter turn radius^',rad)
         endif

         if(rad.le.0.0) then
          write(*,*) '    * Must have  R > 0'
          ninp = 0
          go to 71
         endif

         do jr = ir1, ir2
           parval(iprad,jr) = rad

c--------- go recalculate cl, n
           parval(ipcl ,jr) = 0.
           parval(ipfac,jr) = 0.
         enddo

        endif

c-------------------------------------
      elseif(com(1:1).eq.'L') then
        if    (ktrim.eq.1) then
         write(*,*) 'load factor not specifiable for this trim case'
         go to 10

        elseif(ktrim.eq.2) then
 81      continue
         if(ninp.ge.1) then
          fac = rinp(1)
         else
          call askr('      Enter load factor^',fac)
         endif

         if(fac.le.0.0) then
          write(*,*) '    * Must have  N > 0'
          ninp = 0
          go to 81
         endif

         do jr = ir1, ir2
           parval(ipfac,jr) = fac

c--------- go recalculate v
           parval(ipvee,jr) = 0.
         enddo

        endif

c-------------------------------------
      elseif(com(1:1).eq.'X') then
        if(ninp.ge.1) then
         xcg = rinp(1)
        else
         call askr('      Enter X_cg location^',xcg)
        endif
        do jr = ir1, ir2
          parval(ipxcg,jr) = xcg
          write(*,7500) 'x_cg', jr
        enddo

c-------------------------------------
      elseif(com(1:1).eq.'Y') then
        if(ninp.ge.1) then
         ycg = rinp(1)
        else
         call askr('      Enter Y_cg location^',ycg)
        endif
        do jr = ir1, ir2
          parval(ipycg,jr) = ycg
          write(*,7500) 'y_cg', jr
        enddo

c-------------------------------------
      elseif(com(1:1).eq.'Z') then
        if(ninp.ge.1) then
         zcg = rinp(1)
        else
         call askr('      Enter Z_cg location^',zcg)
        endif
        do jr = ir1, ir2
          parval(ipzcg,jr) = zcg
          write(*,7500) 'z_cg', jr
        enddo

c-------------------------------------
      elseif(com(1:1).eq.'P') then
        if(ninp.ge.1) then
         cd0 = rinp(1)
        else
         call askr('      Enter profile CDo^',cd0)
        endif
        do jr = ir1, ir2
          parval(ipcd0,jr) = cd0
          write(*,7500) 'cdo', jr
        enddo

c------------------------------------------------------
      elseif(com .eq. 'N') then
c------ change name of run case
        if(carg.ne.' ') then
         rtitle(ir) = carg
        else
         write(*,830) rtitle(ir)
 830     format(/' Enter run case name:  ', a)
         read(*,1000) rtnew
         if(rtnew.ne.' ') rtitle(ir) = rtnew
        endif

c-------------------------------------
      else
        write(*,*) '     * Unrecognized parameter'
        go to 10

      endif

      go to 5

      end ! trmset
