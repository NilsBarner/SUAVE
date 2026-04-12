C***********************************************************************
C    Module:  joper.f
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

      subroutine oper
c---------------------------------------
c     Main driver routine for JVL
c---------------------------------------
      include 'jvl.inc'
      include 'jvlplt.inc'
      logical error, lcerr, lcerri, lwrit, lmatch

      character*4 comand, itemc
      character*80 fnout, fnder, fnnew, fndcp
      character*80 line, fnvb, comarg, crun, prompt, rtnew

      logical lowrit

      real    rinput(20), rinp(20)
      integer iinput(20), iinp(20)

      if(.not.lgeo) then
       write(*,*)
       write(*,*) '* Configuration not defined'
       return
      endif

      fnvb = ' '

 1000 format (a)

      lplot = .false.
      lwrit = .false.

      fnout = ' '

c=================================================================
c---- start of user interaction loop
 800  continue

      lcerr = .false.

      call cfrac(irun,nrun,crun,npr)

      write(*,1050) crun(1:npr), rtitle(irun)
      call conlst(irun)
      write(*,1052)

 1050 format(
     &  /' Operation of run case ',a,':  ', a
     &  /' ==========================================================')

 1052 format(
     &  /'  C1  set level or banked  horizontal flight constraints '
     &  /'  C2  set steady pitch rate (looping) flight constraints '
     &  /'  M odify parameters                                     '
     & //' "#" select  run case         L ist defined run cases   '
     &  /'  +  add new run case         S ave run cases to file   '
     &  /'  -  delete  run case         F etch run cases from file'
     &  /'  N ame current run case      W rite forces to file     '
     & //' eX ecute run case            I nitialize variables     '
     &  /'  G eometry plot              T refftz plane plot       '
     &  /'  O ptions                                              '
     & //'  FT  total   forces          FS  surface strip forces '
     &                                          '(FSB bodyaxes)'
     &  /'  FN  surface forces          FE  surface element forces  '
     &  /'  FB  body    forces          FL  body line segment forces'
     & //'  ST  stability derivatives   RE  reference quantities    '
     &  /'  SB  body-axis derivatives   DP  dCp(x) output           '
     & //'  VM  strip shear,moment      DE  design changes          '
     &  /'  HM  hinge moments                                       ')
c 
c   A B C D E F G H I J K L M N O P Q R S T U V W X Y Z
c   x x x x   x x   x     x x x x x   x x x   x x x x


 810  continue
      prompt = ' .OPER (case ' // crun(1:npr) // ')^'
      call askc(prompt,comand,comarg)

c------------------------------------------------------
      if    (comand .eq. '    ') then
       if(lplot) call plend
       lplot = .false.
       call clrzoom
       return

      elseif(comand .eq. '?   ') then
       go to 800

      endif

c------------------------------------------------------
c---- check for run case commands
      if(comand .eq. '+   ') then
c----- add new case after current one

       if(nrun.eq.nrmax) then
        write(*,*)
        write(*,*) '* Run case array limit nrmax reached'
       else
        nrun = nrun + 1

        do jr = nrun, irun+1, -1
          call rcopy(jr,jr-1)
        enddo
        write(*,*) 'Initializing new run case from current one'
c        
        irun = irun + 1
       endif

       go to 800

      elseif(comand .eq. '-   ') then
c----- delete current case

       if(nrun.le.1) then
        write(*,*)
        write(*,*) '* Cannot delete one remaining run case'
       else
        do jr = irun, nrun-1
          call rcopy(jr,jr+1)
        enddo
        nrun = nrun - 1
        irun = max( 1 , min( irun , nrun ) )
       endif

       go to 800

      endif

c------------------------------------------------------
c---- see if command is an integer
      ninput = 1
      call getint(comand,iinput,ninput,error)
      if(.not.error .and. ninput.ge.1 
     &  .and. comand(1:1).ne.'T'
     &  .and. comand(1:1).ne.'F' ) then
c----- command is an integer... new case index?
       irun = max( 1 , min( nrun , iinput(1) ) )
       go to 800
      endif

c------------------------------------------------------
c---- extract command line numeric arguments
      do i=1, 20
        iinput(i) = 0
        rinput(i) = 0.0
      enddo
      ninput = 20
      call getint(comarg,iinput,ninput,error)
      ninput = 20
      call getflt(comarg,rinput,ninput,error)

 14   continue
c------------------------------------------------------
c---- check for parameter toggle/set command
      call conset(comand,comarg,lmatch,irun)

      if(lmatch) then
c----- match found... go back to oper menu
       go to 800
      endif

c------------------------------------------------------
c---- check for trim set command
      call trmset(comand,comarg,lmatch,irun)

      if(lmatch) then
c----- match found... go back to oper menu
       go to 800
      endif

c------------------------------------------------------
c---- pick up here to try decoding for remaining commands

      if(comand .eq. 'X   ') then
c------ execute calculation
        if(lcerr) then
         write(*,*) '** Flow solution is not possible.'
         write(*,*) '** Cannot impose a constraint more than once.'
         go to 800
        endif

        info = 1
        call exec(nitmax,info,irun)

c------ check residual linearizations
c        call linchkr

c------ check gam,vc,vv,ft,mt,dtot,ltot(u,d,j,g) derivatives
c        call linchkf

ccc     if(.not.lsol) go to 810

        if(lptot)   call outtot(6)
        if(lpsurf)  call outsurf(6)
        if(lpbody)  call outbody(6)
        if(lpstrp)  call outstrp(6)
        if(lpelem)  call outelem(6)
        if(lpblin)  call outblin(6)
        if(lphinge) call outhinge(6)

c------------------------------------------------------
      elseif(comand .eq. 'XX  ') then
c------ execute calculation for all run cases
        do 24 ir = 1, nrun
c-------- check for well-posedness
          lcerri = .false.
          do iv = 2, nvtot
            do jv = 2, nvtot
              if(iv.ne.jv .and. icon(iv,ir).eq.icon(jv,ir)) then
               lcerri = .true.
              endif
            enddo
          enddo
          if(lcerri) then
           write(*,*) '** Run case', ir,' ...'
           write(*,*) '** Flow solution is not possible.'
           write(*,*) '** Cannot impose a constraint more than once.'
           go to 24          
          endif

          info = 1
          call exec(nitmax,info,ir)

ccc       if(.not.lsol) go to 24

          if(lptot)  call outtot(6)
          if(lpsurf) call outsurf(6)
          if(lpbody) call outbody(6)
          if(lpstrp) call outstrp(6)
          if(lpelem) call outelem(6)
          if(lpblin)  call outblin(6)
          if(lphinge) call outhinge(6)
 24     continue

c------------------------------------------------------
      elseif(comand .eq. 'M   ') then
        call parmod(irun)

c------------------------------------------------------
      elseif(comand .eq. 'N   ') then
c------ change name of run case
        if(comarg.ne.' ') then
         rtitle(irun) = comarg
        else
         write(*,830) rtitle(irun)
 830     format(/' Enter run case name:  ', a)
         read(*,1000) rtnew
         if(rtnew.ne.' ') rtitle(irun) = rtnew
        endif

c------------------------------------------------------
      elseif(comand .eq. 'DE  ') then
c------ design changes
        if(ndesign.eq.0) then
         write(*,*) '* No design parameters are declared'
         go to 810
        endif

 30     continue
        write(*,1036)
        do n = 1, ndesign
          write(*,1037) n, gname(n), deldes(n)
        enddo
 1036   format(/' ================================================'
     &         /' Current design parameter changes:' /
     &         /'    k   Parameter      change')
ccc              x1234xxx123456789012345612345678901234
 1037   format(1x, i4,3x, a, g14.5) 

        write(*,*)
 35     write(*,*) 'Enter  k, design changes (<return> if done) ...'
 37     read (*,1000) line
        call bstrip(line,nlin)
        if(line(1:1).eq.'?') go to 30
        if(line(1:1).eq.' ') go to 800

        ninp = 40
        call getflt(line,rinp,ninp,error)
        if(error) then
         write(*,*) '* bad input'
         go to 35
        endif

        do i = 1, ninp, 2
          n = int(rinp(i))
          if(n.lt.1 .or. n.gt.ndesign) then
           write(*,*) 'Index k out of bounds. Input ignored.'
          endif
          deldes(n) = rinp(i+1)
          lsol = .false.
          go to 37
        enddo
        go to 30

c------------------------------------------------------
      elseif(comand .eq. 'I   ') then
c------ clear operating parameters
        call runini(irun)
        lsol = .false.

c------------------------------------------------------
      elseif(comand .eq. 'II  ') then
c------ clear operating parameters
        do ir = 1, nrun
          call runini(ir)
        enddo
        lsol = .false.

c------------------------------------------------------
      elseif(comand .eq. 'G   ') then
c------ plot geometry
        call plotvl(azimob, elevob, tiltob, robinv)
        lplot = .true.

c------------------------------------------------------
      elseif(comand .eq. 'T   ') then
c------ plot spanloadings in Trefftz plane
c        if(lsol) then
          call plottp
          lplot = .true.
c        else
c          write(*,*) '* Execute flow calculation first!'
c        endif

c------------------------------------------------------
      elseif(comand .eq. 'FE  ') then
c------ print vortex element forces
c        if(lsol) then
c          call outelem(6)
c        else
c          write(*,*) '* Execute flow calculation first!'
c        endif
        IF(LSOL) THEN
         CALL GETFILE(LU,COMARG)
C
         IF(LU.LE.-1) THEN
          WRITE(*,*) '* Filename error *'
         ELSEIF(LU.EQ.0) THEN
          WRITE(*,*) '* Data not written'
         ELSE
          CALL OUTELEM(LU)
          IF(LU.NE.5 .AND. LU.NE.6) CLOSE(LU)
         ENDIF
C
        ELSE
         WRITE(*,*) '* Execute flow calculation first!'
        ENDIF

c------------------------------------------------------
      elseif(comand .eq. 'FL  ') then
c------ print body line-element forces
c        if(lsol) then
c          call outblin(6)
c        else
c          write(*,*) '* Execute flow calculation first!'
c        endif
        IF(LSOL) THEN
         CALL GETFILE(LU,COMARG)
C
         IF(LU.LE.-1) THEN
          WRITE(*,*) '* Filename error *'
         ELSEIF(LU.EQ.0) THEN
          WRITE(*,*) '* Data not written'
         ELSE
          CALL OUTBLIN(LU)
          IF(LU.NE.5 .AND. LU.NE.6) CLOSE(LU)
         ENDIF
C
        ELSE
         WRITE(*,*) '* Execute flow calculation first!'
        ENDIF

c------------------------------------------------------
      elseif(comand .eq. 'FN  ') then
c------ print surface forces
c        if(lsol) then
c          call outsurf(6)
c        else
c          write(*,*) '* Execute flow calculation first!'
c        endif
        IF(LSOL) THEN
         CALL GETFILE(LU,COMARG)
C
         IF(LU.LE.-1) THEN
          WRITE(*,*) '* Filename error *'
         ELSEIF(LU.EQ.0) THEN
          WRITE(*,*) '* Data not written'
         ELSE
          CALL OUTSURF(LU)
          IF(LU.NE.5 .AND. LU.NE.6) CLOSE(LU)
         ENDIF
C
        ELSE
         WRITE(*,*) '* Execute flow calculation first!'
        ENDIF

c------------------------------------------------------
      elseif(comand .eq. 'FB  ') then
c------ print body forces
c        if(lsol) then
c          call outbody(6)
c        else
c          write(*,*) '* Execute flow calculation first!'
c        endif
        IF(LSOL) THEN
         CALL GETFILE(LU,COMARG)
C
         IF(LU.LE.-1) THEN
          WRITE(*,*) '* Filename error *'
         ELSEIF(LU.EQ.0) THEN
          WRITE(*,*) '* Data not written'
         ELSE
          CALL OUTBODY(LU)
          IF(LU.NE.5 .AND. LU.NE.6) CLOSE(LU)
         ENDIF
C
        ELSE
         WRITE(*,*) '* Execute flow calculation first!'
        ENDIF

c------------------------------------------------------
      elseif(comand .eq. 'FSB ') then
c------ print body axes strip forces
c        if(lsol) then
c          call outstrpb(6)
c        else
c          write(*,*) '* Execute flow calculation first!'
c        endif
        IF(LSOL) THEN
         CALL GETFILE(LU,COMARG)
C
         IF(LU.LE.-1) THEN
          WRITE(*,*) '* Filename error *'
         ELSEIF(LU.EQ.0) THEN
          WRITE(*,*) '* Data not written'
         ELSE
          CALL OUTSTRPB(LU)
           IF(LU.NE.5 .AND. LU.NE.6) CLOSE(LU)
         ENDIF
C
        ELSE
         WRITE(*,*) '* Execute flow calculation first!'
        ENDIF

c------------------------------------------------------
      elseif(comand .eq. 'FS  ') then
c------ print strip forces
c        if(lsol) then
c          call outstrp(6)
c        else
c          write(*,*) '* Execute flow calculation first!'
c        endif
        IF(LSOL) THEN
         CALL GETFILE(LU,COMARG)
C
         IF(LU.LE.-1) THEN
          WRITE(*,*) '* Filename error *'
         ELSEIF(LU.EQ.0) THEN
          WRITE(*,*) '* Data not written'
         ELSE
          CALL OUTSTRP(LU)
          IF(LU.NE.5 .AND. LU.NE.6) CLOSE(LU)
         ENDIF
C
        ELSE
         WRITE(*,*) '* Execute flow calculation first!'
        ENDIF

c------------------------------------------------------
      elseif(comand .eq. 'FT  ') then
c------ print total forces
c        if(lsol) then
c          call outtot(6)
c        else
c          write(*,*) '* Execute flow calculation first!'
c        endif
        IF(LSOL) THEN
         CALL GETFILE(LU,COMARG)
C
         IF(LU.LE.-1) THEN
          WRITE(*,*) '* Filename error *'
         ELSEIF(LU.EQ.0) THEN
          WRITE(*,*) '* Data not written'
         ELSE
          CALL OUTTOT(LU)
           IF(LU.NE.5 .AND. LU.NE.6) CLOSE(LU)
         ENDIF
C
        ELSE
         WRITE(*,*) '* Execute flow calculation first!'
        ENDIF

c------------------------------------------------------
      elseif(comand .eq. 'HM  ') then
c------ print hinge moments
c        if(lsol) then
c          call outhinge(6)
c        else
c          write(*,*) '* Execute flow calculation first!'
c        endif
        IF(LSOL) THEN
         CALL GETFILE(LU,COMARG)
C
         IF(LU.LE.-1) THEN
          WRITE(*,*) '* Filename error *'
         ELSEIF(LU.EQ.0) THEN
          WRITE(*,*) '* Data not written'
         ELSE
          CALL OUTHINGE(LU)
          IF(LU.NE.5 .AND. LU.NE.6) CLOSE(LU)
         ENDIF
C
        ELSE
         WRITE(*,*) '* Execute flow calculation first!'
        ENDIF

c------------------------------------------------------
      elseif(comand .eq. 'VM  ') then
c------ calculate and print shear and bending on surfaces
c        if(lsol) then
          call getvm(fnvb)
c        else
c          write(*,*) '* Execute flow calculation first!'
c        endif

c------------------------------------------------------
      elseif(comand .eq. 'DP  ') then
c------ print vortex element forces
c        if(lsol) then
          if(ninput .ge. 1) then
           jstrip = iinput(1)
          else
 70        call aski('Enter strip index^',jstrip)
           if(jstrip .le. 0) go to 800
           if(jstrip .gt. nstrip) then
            write(*,*) 'Number of strips = ', nstrip
            go to 70
           endif
          endif

          lu = 6
          call outdcp(lu,jstrip)

          call asks('Enter outout filename, or <enter>^',fndcp)
          if(fndcp .eq. ' ') go to 800

          lu = 70
          open(lu,file=fndcp,status='unknown',err=800)
          call outdcp(lu,jstrip)
          close(lu)

c        else
c          write(*,*) '* Execute flow calculation first!'
c        endif

c------------------------------------------------------
      elseif(comand .eq. 'CN  ') then
c------ print a spanloading file
c        if(lsol) then
          call outcnc
c        else
c          write(*,*) '* Execute flow calculation first!'
c        endif

c------------------------------------------------------
      elseif(comand .eq. 'L   ') then
c------ list run cases
        lu = 6
        call runsav(lu)

c------------------------------------------------------
      elseif(comand .eq. 'S   ') then
c------ save run case file
        call bstrip(frndef,nfr)
        write(*,2040) frndef(1:nfr)
 2040   format(' Enter run case filename: ', a)
        read (*,1000) fnnew

        if(fnnew.ne.' ') frndef = fnnew
        open(lurun,file=frndef,status='old',err=42)

        if(lowrit(frndef)) then
         rewind(lurun)
         go to 45
        else
         write(*,*) 'Run cases not saved'
         go to 810
        endif

 42     open(lurun,file=frndef,status='new')

 45     call runsav(lurun)
        close(lurun)

c------------------------------------------------------
      elseif(comand .eq. 'F   ') then
c------ fetch run case file
        call bstrip(frndef,nfr)
        write(*,2050) frndef(1:nfr)
 2050   format(' Enter run case filename: ', a)
        read (*,1000) fnnew
        if(fnnew.ne.' ') frndef = fnnew

        call runget(lurun,frndef,error)
        if(error) then
         go to 810
        else
         write(*,2055) (ir, rtitle(ir), ir=1, nrun)
 2055    format(' Run cases read in ...',
     &          100(/1x,i4,': ',a))
        endif

        irun = min( irun, nrun )

c------------------------------------------------------
      elseif(comand .eq. 'ST  ') then
c------ create stability derivatives
c        if(lsol) then
          call dermats(6)

          write(*,2060)
 2060     format(/' Enter output filename (or <return>): ', $)
          read(*,1000) fnder
          call bstrip(fnder,nfd)
          if(nfd.eq.0) go to 800

          open(lustd,file=fnder,status='old',err=51)
          if(lowrit(fnder)) then
           rewind(lustd)
           go to 52
          else
           close(lustd)
           go to 800
          endif

 51       open(lustd,file=fnder,status='new',err=55)
 52       call dermats(lustd)
          close(lustd)
          go to 800

 55       continue
          write(*,*) '* File OPEN error'

c         else
c          write(*,*) '* Execute flow calculation before derivatives'
c        endif

c------------------------------------------------------
      elseif(comand .eq. 'SM  ') then
c------ create stability derivatives
c        if(lsol) then
          call dermatm(6)

          write(*,2060)
          read(*,1000) fnder
          call bstrip(fnder,nfd)
          if(nfd.eq.0) go to 800

          open(lustd,file=fnder,status='old',err=61)
          if(lowrit(fnder)) then
           rewind(lustd)
           go to 62
          else
           close(lustd)
           go to 800
          endif

 61       open(lustd,file=fnder,status='new',err=65)
 62       call dermatm(lustd)
          close(lustd)
          go to 800

 65       continue
          write(*,*) '* File OPEN error'

c         else
c          write(*,*) '* Execute flow calculation before derivatives'
c        endif

c------------------------------------------------------
      elseif(comand .eq. 'SB  ') then
c------ create stability derivatives
c        if(lsol) then
          call dermatb(6)

          write(*,2060)
          read(*,1000) fnder
          call bstrip(fnder,nfd)
          if(nfd.eq.0) go to 800

          open(lustd,file=fnder,status='old',err=71)
          if(lowrit(fnder)) then
           rewind(lustd)
           go to 72
          else
           close(lustd)
           go to 800
          endif

 71       open(lustd,file=fnder,status='new',err=75)
 72       call dermatb(lustd)
          close(lustd)
          go to 800

 75       continue
          write(*,*) '* File OPEN error'

c         else
c          write(*,*) '* Execute flow calculation before derivatives'
c        endif

c------------------------------------------------------
      elseif(comand .eq. 'W   ') then
c------ write force  data to a file
        call bstrip(fnout,nfn)
        write(*,1080) fnout(1:nfn)
 1080   format(' Enter forces output file: ', a)
        read (*,1000) fnnew

        if(fnnew.ne.' ') then
c-------- new filename was entered...
c-------- if previous file is open, close it
          if(lwrit) close(luout)
          fnout = fnnew
c-------- open new file and write header
          open(luout,file=fnout,status='unknown')
          lwrit = .true.

        else
c-------- just a <return> was entered...
          if(.not.lwrit) then
            write(*,*) 'No action taken.'
            go to 800
          endif

        endif

        if(lptot)   call outtot(luout)
        if(lpsurf)  call outsurf(luout)
        if(lpbody)  call outbody(luout)
        if(lpstrp)  call outstrp(luout)
        if(lpelem)  call outelem(luout)
        if(lpblin)  call outblin(luout)
ccc     if(lpderiv) call dermat(luout)

c------------------------------------------------------
      elseif(comand .eq. 'RE  ') then
c------ change reference data 
 89     write(*,2090) sref,cref,bref, (xyzmom(k), k=1, 3)
 2090   format(/' ==========================='
     &         /'  S ref: ', g11.5,
     &         /'  C ref: ', g11.5,
     &         /'  B ref: ', g11.5,
     &         /'  X mom: ', g11.5,
     &         /'  Y mom: ', g11.5,
     &         /'  Z mom: ', g11.5 )

 90     call askc(' Select item,value^',itemc,comarg)
 2100   format(' Enter new ',a,': ', $)

        if(itemc.eq.'    ') then
          go to 800
        endif

        ninp = 1
        call getflt(comarg,rinp,ninp,error)

        if    (index('Ss',itemc(1:1)).ne.0) then
          if(ninp.eq.0) then
 91        write(*,2100) 'Reference area Sref'
           read (*,*,err=91) sref
          else
           sref = rinp(1)
          endif
          lsol = .false.
          do ir = 1, nrun
            lsolr(ir) = .false.
          enddo

        elseif(index('Cc',itemc(1:1)).ne.0) then
          if(ninp.eq.0) then
 92        write(*,2100) 'Reference chord Cref'
           read (*,*,err=92) cref
          else
           cref = rinp(1)
          endif
          lsol = .false.
          do ir = 1, nrun
            lsolr(ir) = .false.
          enddo

        elseif(index('Bb',itemc(1:1)).ne.0) then
          if(ninp.eq.0) then
 93        write(*,2100) 'Reference span Bref'
           read (*,*,err=93) bref
          else
           bref = rinp(1)
          endif
          lsol = .false.
          do ir = 1, nrun
            lsolr(ir) = .false.
          enddo

        elseif(index('Xx',itemc(1:1)).ne.0) then
          if(ninp.eq.0) then
 94        write(*,2100) 'Moment reference Xmom'
           read (*,*,err=94) xyzmom(1)
          else
           xyzmom(1) = rinp(1)
          endif

        elseif(index('Yy',itemc(1:1)).ne.0) then
          if(ninp.eq.0) then
 95        write(*,2100) 'Moment reference Ymom'
           read (*,*,err=95) xyzmom(2)
          else
           xyzmom(2) = rinp(1)
          endif

        elseif(index('Zz',itemc(1:1)).ne.0) then
          if(ninp.eq.0) then
 96        write(*,2100) 'Moment reference Zmom'
           read (*,*,err=96) xyzmom(3)
          else
           xyzmom(3) = rinp(1)
          endif

        else
          write(*,*) 'item not recognized'
          go to 89

        endif
        go to 89

c------------------------------------------------------
      elseif(comand .eq. 'DJ  ') then
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

        write(*,3020) dvjexp
 3020   format(/1x,' DVjet  ~  (Vinf/Vref)^', f8.4)
        go to 810

c------------------------------------------------------
      elseif(comand .eq. 'O   ') then
        call optget

c------------------------------------------------------
      else
        write(*,*)
        write(*,*) '* Option not recognized'

      endif
      go to 800

      end ! oper



      subroutine conlst(ir)
c-----------------------------------------------------------------
c     Lists the constraint associated with each global variable, 
c     for run case ir
c-----------------------------------------------------------------
      include 'jvl.inc'

      character*4 chss

      write(*,1010)

      do iv = 2, nvtot
        ic = icon(iv,ir)
        chss = '  '
        do jv = 2, nvtot
          if(iv.ne.jv .and. icon(iv,ir).eq.icon(jv,ir)) then
           chss = '**'
          endif
        enddo
        write(*,1020) varkey(iv), connam(ic), conval(ic,ir), chss
      enddo

      write(*,1030)
      return
c      
 1010 format(
     &  /'  variable          constraint              '
     &  /'  ------------      ------------------------')
 1020 format(
     &   '  ',a,'  ->  ', a, '=', g12.4, 1x, a)
 1030 format(
     &   '  ------------      ------------------------')
      end ! conlst



      subroutine conset(comand,comarg,lmatch,ir)
      include 'jvl.inc'
      character*(*) comand, comarg
      logical lmatch

      character*80 prompt
      character*4 arrow
      real    rinp(20)
      integer iinp(20)
      logical error

c---- for control or jet variable, 
c-     first number in arg string should be part of command
      if(comand(1:2) .eq. 'D ' .or.
     &   comand(1:2) .eq. 'J '      ) then
       comand(2:3) = comarg(1:2)
       comarg(1:2) = '  '
       call bstrip(comarg,narg)
      endif

c---- length of non-blank part of command, if any
      kclen = index(comand,' ') - 1
      if(kclen.le.0) kclen = len(comand)

c---- test command against variable keys, using only non-blank part of command
      do iv = 2, nvtot
        kvlen = index(varkey(iv),' ') - 1
        if(kclen .eq. kvlen .and.
     &     comand(1:kclen) .eq. varkey(iv)(1:kvlen)) go to 16
      enddo

c---- no variable key matched... go test for regular commands
      lmatch = .false.
      return

c------------------------------------------------------
c---- found a variable-key match!
 16   continue
      lmatch = .true.
      call touper(comarg) 

      write(*,*) comarg


c---- see if constraint was already specified as second command argument
      kclen = index(comarg,' ') - 1
      if(kclen.le.0) kclen = len(comarg)
      kclen = min( kclen , len(conkey(1)) )
      do ic = 2, nctot
        if(comarg(1:kclen) .eq. conkey(ic)(1:kclen)) go to 18
      enddo

c---- constraint not given... get it from constraint-selection menu
      write(*,1081)
      do ic = 2, nctot
        if(ic.eq.icon(iv,ir)) then
         arrow = '->  '
        else
         arrow = '    '
        endif
        write(*,1082) arrow, conkey(ic), connam(ic), conval(ic,ir)
      enddo
 1081 format(/'       Constraint            Value     '
     &       /'      - - - - - - - - - - - - - - - - -')
 1082 format( '   ', a, a, 2x, a, '=', g12.4)

      prompt= '      Select new  constraint,value  for '
     &        // varnam(iv) // '^'
      call askc(prompt,comand,comarg)
      if(comand .eq. ' ') return

      if(comand(1:2) .eq. 'D ' .or.
     &   comand(1:2) .eq. 'J '      ) then
       comand(2:3) = comarg(1:2)
       comarg(1:2) = '  '
       call bstrip(comarg,narg)
      endif

c---- try to parse command again
      comarg = comand(1:3) // ' ' // comarg
      go to 16

c----------------------------------
c---- pick up here to set new constraint
 18   continue

c---- set new constraint index for selected variable iv
      icon(iv,ir) = ic

c---- see if constraint value was already specified in command argument
      ninp = 1
      call getflt(comarg(kclen+1:80),rinp,ninp,error)
      if(error) ninp = 0

      if(ninp.ge.1) then
c----- yep...  set constraint value to command argument
       conval(ic,ir) = rinp(1)
      else
c----- nope... get constraint value from user (current value is the default)
 19    write(*,1090) connam(ic), conval(ic,ir)
 1090  format(/' enter specified ', a,':', g12.4)
       call readr(1,conval(ic,ir),error)
       if(error) go to 19
      endif

c---- go back to oper menu
      return
      end ! conset




      subroutine shores
      use jvl_inc
      include 'jvl.inc'
      real vinf, rrot(3), vrot(3)

      vinf = sqrt(vbar(1)**2 + vbar(2)**2 + vbar(3)**2)

c----- go over solid surfaces
       do i = 1, nvor
        im = ijetm(i)

        if(im .eq. 0) then
c------- solid-surface element... flow tangency equation
         rrot(1) = rc(1,i) - xyzref(1)
         rrot(2) = rc(2,i) - xyzref(2)
         rrot(3) = rc(3,i) - xyzref(3)
         call cross(rrot,wbar,vrot)

         res = (vbar(1) + vrot(1) + vc(1,i))*enc(1,i)
     &       + (vbar(2) + vrot(2) + vc(2,i))*enc(2,i)
     &       + (vbar(3) + vrot(3) + vc(3,i))*enc(3,i)

c         write(*,*) i, res

       else
c------- jet element... jet curvature/loading equation
         js = jstripv(i)

         dvjf = vinf**dvjexp
         dvjet = 0.
         do n = 1, nvarjet
           dvjet = dvjet + deljet(n)*gjstrp(js,n)
         enddo

         vjet = vinf + dvjet*dvjf

         vjmin = 0.1
         if(vjet .lt. vjmin) then
          write(*,*) '?  vjet <', vjmin
          vjet = vjmin
         endif

         rhjet = 1.0
         fh = fhstrp(js)
         hjet = hdstrp(js)*(1.0 + 0.5*fh*(vinf/vjet - 1.0))
         djp = (rhjet*vjet**2 - vinf**2)*hjet

         js = jstripv(im)
         if(istype(js) .eq. 0) then
c-------- upstream element is on solid te
          delf = 0.
          do n = 1, ncontrol
            delf = delf + delcon(n)*enc_d(1,im,n)
          enddo
          ddte = dj0strp(js)
     &         + dj1strp(js)*delf
     &         + dj3strp(js)*delf**3

         else
          delf = 0.
          ddte = 0.

         endif

         dvcn = 0.
         do j = 1, nvor
           dvcn_gam = (vc_gam(1,i,j)-vc_gam(1,im,j))*enc(1,i)
     &              + (vc_gam(2,i,j)-vc_gam(2,im,j))*enc(2,i)
     &              + (vc_gam(3,i,j)-vc_gam(3,im,j))*enc(3,i)
           dvcn = dvcn + dvcn_gam*gam(j)
         enddo

         dwcn = 0.
         do k = 1, 3
           n = k
           dwcn_u = (wc_u(1,i,n)-wc_u(1,i,n))*enc(1,i)
     &            + (wc_u(2,i,n)-wc_u(2,i,n))*enc(2,i)
     &            + (wc_u(3,i,n)-wc_u(3,i,n))*enc(3,i)
           dwcn = dwcn + dwcn_u*vbar(k)
         enddo
         do k = 1, 3
           n = k+3
           dwcn_u = (wc_u(1,i,n)-wc_u(1,i,n))*enc(1,i)
     &            + (wc_u(2,i,n)-wc_u(2,i,n))*enc(2,i)
     &            + (wc_u(3,i,n)-wc_u(3,i,n))*enc(3,i)
           dwcn = dwcn + dwcn_u*wbar(k)
         enddo

         res = djp*((dvcn + dwcn)/vinf + ddte) - gam(i)

c         write(*,*) -i, res, djp

        endif

      enddo

      return
      end ! shores


      subroutine optget
c-------------------------------------------------
c     Allows toggling and setting of various 
c     printing and plotting stuff.
c-------------------------------------------------
      include 'jvl.inc'
      character*4 itemc
      character*80 comarg
      character*50 satype, rottype
      logical error, err1

      real    rinput(20)
      integer iinput(20)
      logical linput(20)

 1000 format(a)

      call getsa(lnasa_sa,satype,dir)

 100  continue
      if(lsa_rates) then
        rottype =
     &         'Rates,moments about Stability axes, x along Vinf'
      else
        rottype =
     &         'Rates,moments about Body axes, x along geometric x axis'
      endif
  
      write(*,1110) lptot,lpsurf,lpstrp,lpelem,lpbody,lpblin,
     &              lphinge,lpderiv,
     &              ltrforce,lvisc,lbforce,lclvjet,
     &                satype,rottype,izsym,zsym,
     &                saxfr,
     &                vrcorec,vrcorew,srcore,
     &                nitmax
 1110   format(
     &  /'   ======================================'
     &  /'    P rint default output for...'
     &  /'        total     :  ',l2,
     &  /'        surfaces  :  ',l2,
     &  /'        strips    :  ',l2,
     &  /'        elements  :  ',l2,
     &  /'        bodies    :  ',l2,
     &  /'        body pts. :  ',l2,
     & //'    H inge mom. output:  ',l2,
     &  /'    D erivative output:  ',l2,
     & //'    T rail.leg forces:  ',l2,
     &  /'    V iscous forces  :  ',l2,
     &  /'    B ody forces     :  ',l2,
     &  /'    CLV jet forces   :  ',l2,
     & //'    A xis orient. :  ', a, 
     &  /'    R ate,mom axes:  ', a,
     &  /'    Z  symmetry   :  ',i2,' @ z =',f10.4
     &  /'    S pan axis x/c:  ',f10.4
     &  /'    CC vortex core radius/chord:',f10.4,
     &  /'    CW vortex core radius/stripwidth:',f10.4,
     &  /'    CS source core radius/segmentlength:',f10.4,
     & //'    I teration limit:',i5 )
c
c   A B C D E F G H I J K L M N O P Q R S T U V W X Y Z
c   x x x x       x x             x   x x x   x       x

      call askc(' ..Select item to change^',itemc,comarg)

c------------------------------------------------------
      if    (itemc.eq.'    ') then
        return

c---------------------------------
      elseif(itemc.eq.'A   ') then
        lnasa_sa = .not.lnasa_sa
        call getsa(lnasa_sa,satype,dir)

c---------------------------------
      elseif(itemc.eq.'R   ') then
        lsa_rates = .not.lsa_rates

c---------------------------------
      elseif(itemc.eq.'V   ') then
        lvisc = .not.lvisc
        if(lvisc) then
          write(*,*) 'Forces will include profile drag'
         else
          write(*,*) 'Forces will not include profile drag'
        endif
        do ir = 1, nrun
          lsolr(ir) = .false.
        enddo

c---------------------------------
      elseif(itemc.eq.'B   ') then
        lbforce = .not.lbforce
        if(lbforce) then
          write(*,*) 'Forces will include body forces'
         else
          write(*,*) 'Forces will not include body forces'
        endif
        do ir = 1, nrun
          lsolr(ir) = .false.
        enddo

c---------------------------------
      elseif(itemc.eq.'CLV ') then
        lclvjet = .not.lclvjet
        if(lclvjet) then
          write(*,*) 'Viscous CD(CL) includes jet forces'
         else
          write(*,*) 'Viscous CD(CL) will not include jet forces'
        endif
        do ir = 1, nrun
          lsolr(ir) = .false.
        enddo

c---------------------------------
      elseif(itemc.eq.'T   ') then
        ltrforce = .not.ltrforce
        if(ltrforce) then
          write(*,*) 'Forces on trailing legs will be included'
         else
          write(*,*) 'Forces on trailing legs will not be included'
        endif
        do ir = 1, nrun
          lsolr(ir) = .false.
        enddo

c---------------------------------
      elseif(itemc.eq.'P   ') then
 128   if(comarg(1:1).ne.' ') then
         ninp = 6
         call getlgv(comarg,linput,ninp,error)
         if(error) go to 130

         if(ninp.ge.1) lptot   = linput(1)
         if(ninp.ge.2) lpsurf  = linput(2)
         if(ninp.ge.3) lpstrp  = linput(3)
         if(ninp.ge.4) lpelem  = linput(4)
         if(ninp.ge.5) lpbody  = linput(5)
         if(ninp.ge.6) lpblin  = linput(6)

         go to 100
        endif

 130    write(*,2100)
     &  'Enter print flags T/F (total,surf,strip,elem,body,bodypt.)'
 2100   format(1x,a,': ', $)

        read(*,1000) comarg
        if(comarg.eq.' ') then
         go to 100
        else
         go to 128
        endif

c---------------------------------
      elseif(itemc.eq.'H   ') then
        lphinge = .not.lphinge

c---------------------------------
      elseif(itemc.eq.'D   ') then
        lpderiv = .not.lpderiv

c---------------------------------
      elseif(itemc.eq.'Z   ') then
       write(*,*) ' '
       write(*,*) 'Currently:'
       if(izsym.eq.0) then
         write (*,1015)
        elseif(izsym.gt.0) then
         write (*,1016) zsym
        elseif(izsym.lt.0) then
         write (*,1017) zsym
       endif
       write(*,*) 'Enter symmetry flag: -1 Free surface'
       write(*,*) '                      0 no Z symmetry'
       write(*,*) '                      1 Ground plane'
       zsymin = 0.0
       read(*,*,err=100) izsymin
       if(izsymin.ne.0.0) then
         write(*,2100) 'Enter Z for symmetry plane'
         read(*,*,err=100) zsymin
       endif
       izsym = izsymin
       zsym  = zsymin
       
c---- reset all z symmetry flags to new symmetry setting
       do n = 1, nsurf
         izimn(n) = izsym         ! z-image flag for surfaces
         do js = jfrst(n), jlast(n)
           izims(js) = izsym      ! z-image flag for strips
           do i = ifrsts(js), ilasts(js)
             izimv(i) = izsym     ! z-image flag for vortices    
           enddo
         enddo
       enddo
       do n = 1, nbody
         izimb(n) = izsym         ! z-image flag for bodies
         do l = lfrst(n), llast(n)
           iziml(l) = izsym       ! z-image flag for body segments
         enddo
       enddo
              
       laic = .false.
       lsol = .false.
       do ir = 1, nrun
         lsolr(ir) = .false.
       enddo

 1015  format(' Z symmetry: No symmetry assumed')
 1016  format(' Z symmetry: Ground plane at Zsym =',f10.4)
 1017  format(' Z symmetry: Free surface at Zsym =',f10.4)

c---------------------------------
      elseif(itemc.eq.'S   ') then
        ninp = 1
        call getflt(comarg,rinput,ninp,error)

        if(error .or. ninp.le.0) then
         rinput(1) = saxfr
         write(*,1030) rinput(1)
 1030    format(/' Enter x/c location of spanwise ref. axis:', f10.4)
         call readr(1,rinput(1),err1)
         if(err1) go to 100
        endif

        saxfr = max( 0.0 , min(1.0,rinput(1)) )
        call encalc
        call aero

c---------------------------------
      elseif(itemc.eq.'CC  ') then
        ninp = 1
        call getflt(comarg,rinput,ninp,error)

        if(error .or. ninp.lt.1) then
         rinput(1) = vrcorec
         write(*,1041) rinput(1)
 1041    format(/' Enter vortex core radius/stripchord ratio:', f10.4)
         call readr(1,rinput(1),err1)
         if(err1) go to 100
        endif

        vrcorec = max( 0.0 , rinput(1) )
        call encalc
        laic = .false.
        lsol = .false.
        do ir = 1, nrun
          lsolr(ir) = .false.
        enddo

c---------------------------------
      elseif(itemc.eq.'CW  ') then
        ninp = 1
        call getflt(comarg,rinput,ninp,error)

        if(error .or. ninp.lt.1) then
         rinput(1) = vrcorew
         write(*,1042) rinput(1)
 1042    format(/' Enter vortex core radius/stripwidth ratio:', f10.4)
         call readr(1,rinput(1),err1)
         if(err1) go to 100
        endif
 
        vrcorew = max( 0.0 , rinput(1) )
        call encalc
        laic = .false.
        lsol = .false.
        do ir = 1, nrun
          lsolr(ir) = .false.
        enddo

c---------------------------------
      elseif(itemc.eq.'CS  ') then
        ninp = 1
        call getflt(comarg,rinput,ninp,error)

        if(error .or. ninp.lt.1) then
         rinput(1) = srcore
         write(*,1043) rinput(1)
 1043    format(/' Enter body core/segment length ratio:', f10.4)
         call readr(1,rinput(1),err1)
         if(err1) go to 100
        endif

        srcore = max( 0.0 , rinput(1) )
        call encalc
        laic = .false.
        lsol = .false.
        do ir = 1, nrun
          lsolr(ir) = .false.
        enddo

c---------------------------------
      elseif(itemc.eq.'I   ') then
        ninp = 1
        call getint(comarg,iinput,ninp,error)

        if(error .or. ninp.le.0) then
         iinput(1) = nitmax
         write(*,1050) iinput(1)
 1050    format(/' Enter max number of iterations:', i5)
         call readi(1,iinput,error)
         if(error) go to 100
        endif

        nitmax = iinput(1)
c---------------------------------
      else
        write(*,*) 'Item not recognized'
        go to 100
      endif
      go to 100

      end ! optget



      subroutine cfrac(irun,nrun,cpr,npr)
c------------------------------------------------------------------
c     Generates a run-case number string cpr for the OPER prompt.
c     Examples of the returned string and string length:
c       for irun=2, nrun=7 :  cpr  =  '2/7',   ncpr = 3
c       for irun=3, nrun=12:  cpr  =  '3/12',  ncpr = 4
c------------------------------------------------------------------
      character*(*) cpr

      izero = ichar('0')
      iten = irun/10
      ione = irun - 10*(irun/10)
      if(iten.le.0) then
       cpr = char(izero+ione) // '/'
      else
       cpr = char(izero+iten) // char(izero+ione) // '/'
      endif

      npr = index(cpr,'/')
      iten = nrun/10
      ione = nrun - 10*(nrun/10)
      if(iten.le.0) then
       cpr = cpr(1:npr)
     &    // char(izero+ione) // '^'
      else
       cpr = cpr(1:npr)
     &    // char(izero+iten) // char(izero+ione) // '^'
      endif

      npr = index(cpr,'^') - 1

      return
      end ! cfrac



      subroutine rcopy(irset,ir)
      include 'jvl.inc'

      do iv = 1, ivmax
        varval(iv,irset) = varval(iv,ir)
      enddo
      do iv = 2, nvtot
        icon(iv,irset) = icon(iv,ir)
      enddo
      do ic = 2, nctot
        conval(ic,irset) = conval(ic,ir)
      enddo
      do ip = 1, nptot
        parval(ip,irset) = parval(ip,ir)
      enddo

      rtitle(irset) = rtitle(ir)
      itrim(irset) = itrim(ir)
      neigen(irset) = neigen(ir)

      do ke = 1, jemax
        eval(ke,irset) = eval(ke,ir)
        do je = 1, jemax
          evec(ke,je,irset) = evec(ke,je,ir)
        enddo
      enddo

      lsolr(irset) = lsolr(ir)

      return
      end ! rcopy



      subroutine parmod(ire)
c----------------------------------------------------
c     Routine for interactive modification 
c     of flight-dynamics system parameters
c----------------------------------------------------
      include 'jvl.inc'

      character*4 com
      character*80 carg, prompt, rtnew
      logical error, repkey

      integer iinp(10)
      real rinp(10)

      if(ire.eq.0) then
       ir = 1
      else
       ir = ire
      endif

c
 1000 format(a)

 2000 format(/'     Parameters of run case ',a,':  ', a)
 2005 format( '     ', a
     &       /'     =================================================')
 2105 format(6x, a, g10.4, 2x, a)
 2110 format(6x, a, a )

c--------------------------------------------------------------------------
c----- jump back here just for menu
 10    continue
       call cfrac(ir,nrun,prompt,npr)
       write(*,2000) prompt(1:npr),rtitle(ir)

       phi   = parval(ipphi ,ir)
       the   = parval(ipthe ,ir)
       mach  = parval(ipmach,ir)
       vee   = parval(ipvee ,ir)
       rho   = parval(iprho ,ir)
       gee   = parval(ipgee ,ir)
       rmass = parval(ipmass,ir)
       rixx  = parval(ipixx ,ir)
       riyy  = parval(ipiyy ,ir)
       rizz  = parval(ipizz ,ir)
       xcg   = parval(ipxcg ,ir)
       ycg   = parval(ipycg ,ir)
       zcg   = parval(ipzcg ,ir)
       cd0   = parval(ipcd0 ,ir)

       hx    = parval(iphx  ,ir)
       hy    = parval(iphy  ,ir)
       hz    = parval(iphz  ,ir)

       dvjexp = parval(ipdvjx,ir)

       dcl_a  = parval(ipcla ,ir)
       dcl_u  = parval(ipclu ,ir)
       dcl_ad = parval(ipclad,ir)

       dcd_a  = parval(ipcda ,ir)
       dcd_u  = parval(ipcdu ,ir)
       dcd_ad = parval(ipcdad,ir)

       dcm_a  = parval(ipcma ,ir)
       dcm_u  = parval(ipcmu ,ir)
       dcm_ad = parval(ipcmad,ir)

       write(*,2105) 'B  bank      = ', phi  , 'deg'
       write(*,2105) 'E  elevation = ', the  , 'deg'
       write(*,2105) 'MN Mach no.  = ', mach , ' '
       write(*,2105) 'V  velocity  = ', vee  , unchv(1:nuv)
       write(*,2105) 'D  air dens. = ', rho  , unchd(1:nud)
       write(*,2105) 'G  grav.acc. = ', gee  , uncha(1:nua)
       write(*,2105) 'M  mass      = ', rmass, unchm(1:num)
       write(*,2105) 'IX Ixx       = ', rixx , unchi(1:nui)
       write(*,2105) 'IY Iyy       = ', riyy , unchi(1:nui)
       write(*,2105) 'IZ Izz       = ', rizz , unchi(1:nui)
       write(*,2105) 'X  X_cg      = ', xcg  , 'lunit'
       write(*,2105) 'Y  Y_cg      = ', ycg  , 'lunit'
       write(*,2105) 'Z  Z_cg      = ', zcg  , 'lunit'
       write(*,2105) 'HX hx        = ', hx   , unchh(1:nuh)
       write(*,2105) 'HY hy        = ', hy   , unchh(1:nuh)
       write(*,2105) 'HZ hz        = ', hz   , unchh(1:nuh)
       write(*,2105) 'CD CDo       = ', cd0  , ' '
       write(*,2105) 'DJ DVj exp.  = ', dvjexp, ' '
       write(*,2105) 'LA dCL_a     = ', dcl_a, ' '
       write(*,2105) 'LU dCL_u     = ', dcl_u, ' '
       write(*,2105) 'LT dCL_adot  = ', dcl_ad, ' '
       write(*,2105) 'DA dCD_a     = ', dcd_a, ' '
       write(*,2105) 'DU dCD_u     = ', dcd_u, ' '
       write(*,2105) 'DT dCD_adot  = ', dcd_ad, ' '
       write(*,2105) 'MA dCM_a     = ', dcm_a, ' '
       write(*,2105) 'MU dCM_u     = ', dcm_u, ' '
       write(*,2105) 'MT dCM_adot  = ', dcm_ad, ' '

       call askc('     Enter parameter, value  (or  # - + N )',com,carg)

       if(com.eq.'    ') then
c------ just a return entered... go back
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

       go to 10

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

       go to 10

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
        write(*,*) '* selected new run case is not defined'
        go to 10
       else
c------ valid new run case selected... go back to top of menu
        ir = iinp(1)
        go to 10
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

c==============================================================================
c---- now decode regular parameter value commands

c------------------------------------
      if(com(1:2) .eq. 'B ') then
 11     continue
        if(ninp.ge.1) then
         phi = rinp(1)
        else
         call askr('      Enter bank angle^',phi)
        endif

        if(phi.le.-90.0.or.phi.ge.90.0) then
         write(*,*) '    * Must have  -90 < bank < +90'
         ninp = 0
         go to 11
        endif

        do jr = ir1, ir2
          parval(ipphi,jr) = phi
        enddo

c-------------------------------------
      elseif(com(1:2) .eq. 'E ') then
 21     continue
        if(ninp.ge.1) then
         the = rinp(1)
        else
         call askr('      Enter elevation angle^',the)
        endif

        if(the.le.-90.0.or.the.ge.90.0) then
         write(*,*) '    * Must have  -90 < elevation < +90'
         ninp = 0
         go to 21
        endif

        do jr = ir1, ir2
          parval(ipthe,jr) = the
        enddo

c-------------------------------------
      elseif(com(1:2) .eq. 'MN') then
 24     continue
        if(ninp.ge.1) then
         mach = rinp(1)
        else
         call askr('      Enter Mach Number^',mach)
        endif

        if(mach .lt. 0.0 .or. mach .gt. 0.999) then
         write(*,*) '    * Must have 0 < Mach < 0.999'
         go to 24
        endif

        do jr = ir1, ir2
          parval(ipmach,jr) = mach
          lsolr(jr) = .false.
        enddo

c-------------------------------------
      elseif(com(1:2) .eq. 'V ') then
 31     continue
        if(ninp.ge.1) then
         vee = rinp(1)
        else
         call askr('      Enter Velocity^',vee)
        endif

        if(vee.le.0.0) then
         write(*,*) '    * Must Have  V > 0'
         ninp = 0
         go to 31
        endif

        do jr = ir1, ir2
          parval(ipvee,jr) = vee
        enddo

c-------------------------------------
      elseif(com(1:2) .eq. 'M ') then
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
        enddo

c-------------------------------------
      elseif(com(1:2) .eq. 'IX') then
 42     continue
        if(ninp.ge.1) then
         rixx = rinp(1)
        else
         call askr('      Enter Ixx^',rixx)
        endif

        if(rixx.le.0.0) then
         write(*,*) '    * Must have  Ixx > 0'
         ninp = 0
         go to 42
        endif

        do jr = ir1, ir2
          parval(ipixx,jr) = rixx
        enddo

c-------------------------------------
      elseif(com(1:2) .eq. 'IY') then
 43     continue
        if(ninp.ge.1) then
         riyy = rinp(1)
        else
         call askr('      Enter Iyy^',riyy)
        endif

        if(riyy.le.0.0) then
         write(*,*) '    * Must have  Iyy > 0'
         ninp = 0
         go to 43
        endif

        do jr = ir1, ir2
          parval(ipiyy,jr) = riyy
        enddo

c-------------------------------------
      elseif(com(1:2) .eq. 'IZ') then
 44     continue
        if(ninp.ge.1) then
         rizz = rinp(1)
        else
         call askr('      Enter Izz^',rizz)
        endif

        if(rizz.le.0.0) then
         write(*,*) '    * Must have  Izz > 0'
         ninp = 0
         go to 44
        endif

        do jr = ir1, ir2
          parval(ipizz,jr) = rizz
        enddo

c-------------------------------------
      elseif(com(1:2) .eq. 'D ') then
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
        enddo

c-------------------------------------
      elseif(com(1:2) .eq. 'G ') then
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
        enddo

c-------------------------------------
      elseif(com(1:2) .eq. 'X ') then
        if(ninp.ge.1) then
         xcg = rinp(1)
        else
         call askr('      Enter X_cg location^',xcg)
        endif

        do jr = ir1, ir2
          parval(ipxcg,jr) = xcg
          lsolr(jr) = .false.
        enddo

c-------------------------------------
      elseif(com(1:2) .eq. 'Y ') then
        if(ninp.ge.1) then
         ycg = rinp(1)
        else
         call askr('      Enter Y_cg location^',ycg)
        endif

        do jr = ir1, ir2
          parval(ipycg,jr) = ycg
          lsolr(jr) = .false.
        enddo

c-------------------------------------
      elseif(com(1:2) .eq. 'Z ') then
        if(ninp.ge.1) then
         zcg = rinp(1)
        else
         call askr('      Enter Z_cg location^',zcg)
        endif

        do jr = ir1, ir2
          parval(ipzcg,jr) = zcg
          lsolr(jr) = .false.
        enddo

c-------------------------------------
      elseif(com(1:2) .eq. 'HX') then
        if(ninp.ge.1) then
         hx = rinp(1)
        else
         call askr('      Enter onboard hx component^',hx)
        endif

        do jr = ir1, ir2
          parval(iphx,jr) = hx
        enddo

c-------------------------------------
      elseif(com(1:2) .eq. 'HY') then
        if(ninp.ge.1) then
         hy = rinp(1)
        else
         call askr('      Enter onboard hy component^',hy)
        endif

        do jr = ir1, ir2
          parval(iphy,jr) = hy
        enddo

c-------------------------------------
      elseif(com(1:2) .eq. 'HZ') then
        if(ninp.ge.1) then
         hz = rinp(1)
        else
         call askr('      Enter onboard hz component^',hz)
        endif

        do jr = ir1, ir2
          parval(iphz,jr) = hz
        enddo

c-------------------------------------
      elseif(com(1:2) .eq. 'CD') then
        if(ninp.ge.1) then
         cd0 = rinp(1)
        else
         call askr('      Enter profile CDo^',cd0)
        endif

        do jr = ir1, ir2
          parval(ipcd0,jr) = cd0
          lsolr(jr) = .false.
        enddo

c-------------------------------------
      elseif(com(1:2) .eq. 'DJ') then
        if(ninp.ge.1) then
         dvjexp = rinp(1)
        else
         call askr(
     &  '      Enter exponent of DVjet variation with Vinf:',dvjexp)
        endif

        do jr = ir1, ir2
         parval(ipdvjx,jr) = dvjexp
         lsolr(jr) = .false.
        enddo

c-------------------------------------
      elseif(com(1:2) .eq. 'LA') then
        if(ninp.ge.1) then
         dcl_a = rinp(1)
        else
         call askr('      Enter added CL_a^',dcl_a)
        endif

        do jr = ir1, ir2
          parval(ipcla,jr) = dcl_a
        enddo

c-------------------------------------
      elseif(com(1:2) .eq. 'LU') then
        if(ninp.ge.1) then
         dcl_u = rinp(1)
        else
         call askr('      Enter added CL_u^',dcl_u)
        endif

        do jr = ir1, ir2
          parval(ipclu,jr) = dcl_u
        enddo

c-------------------------------------
      elseif(com(1:2) .eq. 'LT') then
        if(ninp.ge.1) then
         dcl_ad = rinp(1)
        else
         call askr('      Enter added CL_adot^',dcl_ad)
        endif

        do jr = ir1, ir2
          parval(ipclad,jr) = dcl_ad
        enddo

c-------------------------------------
      elseif(com(1:2) .eq. 'DA') then
        if(ninp.ge.1) then
         dcd_a = rinp(1)
        else
         call askr('      Enter added CD_a^',dcd_a)
        endif

        do jr = ir1, ir2
          parval(ipcda,jr) = dcd_a
        enddo

c-------------------------------------
      elseif(com(1:2) .eq. 'DU') then
        if(ninp.ge.1) then
         dcd_u = rinp(1)
        else
         call askr('      Enter added CD_u^',dcd_u)
        endif

        do jr = ir1, ir2
          parval(ipcdu,jr) = dcd_u
        enddo

c-------------------------------------
      elseif(com(1:2) .eq. 'DT') then
        if(ninp.ge.1) then
         dcd_at = rinp(1)
        else
         call askr('      Enter added CD_adot^',dcd_at)
        endif

        do jr = ir1, ir2
          parval(ipcdad,jr) = dcd_at
        enddo

c-------------------------------------
      elseif(com(1:2) .eq. 'MA') then
        if(ninp.ge.1) then
         dcm_a = rinp(1)
        else
         call askr('      Enter added CM_a^',dcm_a)
        endif

        do jr = ir1, ir2
          parval(ipcma,jr) = dcm_a
        enddo

c-------------------------------------
      elseif(com(1:2) .eq. 'MU') then
        if(ninp.ge.1) then
         dcm_u = rinp(1)
        else
         call askr('      Enter added CM_u^',dcm_u)
        endif

        do jr = ir1, ir2
          parval(ipcmu,jr) = dcm_u
        enddo

c-------------------------------------
      elseif(com(1:2) .eq. 'MT') then
        if(ninp.ge.1) then
         dcm_ad = rinp(1)
        else
         call askr('      Enter added CM_adot^',dcm_ad)
        endif

        do jr = ir1, ir2
          parval(ipcmad,jr) = dcm_ad
        enddo

c------------------------------------------------------
      elseif(com(1:1) .eq. 'N') then
c------ change name of run case
        if(carg.ne.' ') then
         rtitle(ir) = carg
        else
         write(*,830) rtitle(ir)
 830     format(/' Enter run case name:  ', a)
         read(*,1000) rtnew
         if(rtnew.ne.' ') rtitle(ir) = rtnew
        endif

      else
        write(*,*) '     * Unrecognized parameter'

      endif

      if(repkey) write(*,*) '      Value set for all run cases'
      go to 10

      end ! parmod


      subroutine gcoeffs(ir,drm)
c-------------------------------------------------------------------
c     Shifts all moments by offset arm drm.
c
c     Computes force and moment coefficients, and their 
c     stability and control derivatives for case ir.
c     These are used for display and output in OPER,
c     and for eigenmode calculation in MODE.
c
c Inputs: 
c     ft(.)   total force vector / (rhoref Vref^2)
c     ft_u(.,1:3)  d ft(.) / d Vinf(1:3)
c     ft_u(.,4:6)  d ft(.) / d Wrot(1:3)
c     ft_d(.,.)    d ft(.) / d control(.)
c     ft_j(.,.)    d ft(.) / d jetcontrol(.)
c
c     mt(.)   total moment vector/ (rhoref Vref^2)
c     mt_u(.,1:3)  d mt(.) / d Vinf(1:3)
c     mt_u(.,4:6)  d mt(.) / d Wrot(1:3)
c     mt_d(.,.)    d mt(.) / d control(.)
c     mt_j(.,.)    d mt(.) / d jetcontrol(.)
c
c Outputs: 
c     cf(.)      Cx,Cy,Cz
c     cf_v(.,.)  Cxu,Cyu,Czu, Cxv,Cyv,...  Czw
c     cf_w(.,.)  Cxp,Cyp,Czp, Cyq,Cyq,...  Czr
c     cf_d(.,.)  Cxd1,Cyd1,Czd1, Cxd2,Cyd2, ...
c     cf_j(.,.)  Cxj1,Cyj1,Czj1, Cxj2,Cyj2, ...

c     cm(.)      Cl,Cm,Cn
c     cm_v(.,.)  Clu,Cmu,Cnu, Clv,Cmv,...  Cnw
c     cm_w(.,.)  Clp,Cmp,Cnp, Cmq,Cmq,...  Cnr
c     cm_d(.,.)  Cld1,Cmd1,Cnd1, Cld2,Cmd2, ...
c     cm_j(.,.)  Clj1,Cmj1,Cnj1, Clj2,Cmj2, ...
c
c
c     All these outputs can be overwritten here by external data.
c     But note that all forces,moments,velocities,rotation rates 
c     here are in JVL's geometry axes, which have x,z reversed
c     from the standard body axes which have x forward and z down.
c     So the following quantities here have reversed signs from 
c     their conventional body-axis definitions:
c
c         Cx,Cz, Cl,Cn, u,v, p,r
c
c     This therefore requires negating some derivatives.
c     For example,
c      Cmu,Cmw  will be negated
c      Cxu,Cxw  will be negated twice (i.e. unchanged)
c      all Cx,Cz,Cl,Cn control derivatives  will be negated
c
c-------------------------------------------------------------------
      include 'jvl.inc'
      real drm(3)

      real lref(3), dr(3)
      
      integer icrs(3), jcrs(3)
      data icrs / 2, 3, 1 / , jcrs / 3, 1, 2 /

c---- added-on derivative contributions
c-     (the *_adot derivatives are added on in sysmat)
      dcl_u    = parval(ipclu,ir)
      dcl_a    = parval(ipcla,ir)
c      dcl_adot = parval(ipclad,ir)

      dcd_u    = parval(ipcdu,ir)
      dcd_a    = parval(ipcda,ir)
c      dcd_adot = parval(ipcdad,ir)

      dcm_u    = parval(ipcmu,ir)
      dcm_a    = parval(ipcma,ir)
c      dcm_adot = parval(ipcmad,ir)

c---- trim alpha, beta
      alfa = parval(ipalfa,ir)*dtr
      beta = parval(ipbeta,ir)*dtr
      sina = sin(alfa)
      cosa = cos(alfa)

      lref(1) = bref
      lref(2) = cref
      lref(3) = bref

c---- force and moment coefficients in JVL geometry axes
      do k = 1, 3
c------ Cx,Cy,Cz
        cf(k) = ft(k) * 2.0/sref

c------ Cl,Cm,Cn
        cm(k) = mt(k) * 2.0/(sref*lref(k))

c------ stability and control derivatives, in JVL geometry axes
c-       u(1:3) is Vinf(1:3), which is -u,-v,-w
c-       u(4:6) is Wrot(1:3), which is  p, q, r
        do n = 1, 3
c-------- Cxu,Cxv ... Cnw
          cf_v(k,n) = -ft_u(k,n) * 2.0/ sref
          cm_v(k,n) = -mt_u(k,n) * 2.0/(sref*lref(k))

c-------- Cxp,Cxq ... Cnr
          cf_w(k,n) =  ft_u(k,n+3) * 2.0/ sref          * 2.0/lref(n)
          cm_w(k,n) =  mt_u(k,n+3) * 2.0/(sref*lref(k)) * 2.0/lref(n)
        enddo

c------ Cxd ... Cnd
        do n = 1, ncontrol
          cf_d(k,n) = ft_d(k,n) * 2.0/ sref
          cm_d(k,n) = mt_d(k,n) * 2.0/(sref*lref(k))
        enddo

c------ Cxj ... Cnj
        do n = 1, nvarjet
          cf_j(k,n) = ft_j(k,n) * 2.0/ sref
          cm_j(k,n) = mt_j(k,n) * 2.0/(sref*lref(k))
        enddo

c------ Cxg ... Cng
        do n = 1, ndesign
          cf_g(k,n) = ft_g(k,n) * 2.0/ sref
          cm_g(k,n) = mt_g(k,n) * 2.0/(sref*lref(k))
        enddo
      enddo

      dr(1) = drm(1)/lref(1)
      dr(2) = drm(2)/lref(2)
      dr(3) = drm(3)/lref(3)
      do k = 1, 3
        ic = icrs(k)
        jc = jcrs(k)

        cm(k) = cm(k) + dr(ic)*cf(jc) - dr(jc)*cf(ic)

        do n = 1, 3
          cm_v(k,n) = cm_v(k,n) + dr(ic)*cf_v(jc,n) - dr(jc)*cf_v(ic,n)
          cm_w(k,n) = cm_w(k,n) + dr(ic)*cf_w(jc,n) - dr(jc)*cf_w(ic,n)
        enddo

        do n = 1, ncontrol
          cm_d(k,n) = cm_d(k,n) + dr(ic)*cf_d(jc,n) - dr(jc)*cf_d(ic,n)
        enddo

        do n = 1, nvarjet
          cm_j(k,n) = cm_j(k,n) + dr(ic)*cf_j(jc,n) - dr(jc)*cf_j(ic,n)
        enddo

        do n = 1, ndesign
          cm_g(k,n) = cm_g(k,n) + dr(ic)*cf_g(jc,n) - dr(jc)*cf_g(ic,n)
        enddo
      enddo

c---------------------------------------------------------------------------
c---- extra imposed derivatives, rotated into geometry axes by trim alpha_0

c---- negate signs since u here is +back, and standard CDu has u +forward
      dcx_u = -(dcd_u*cosa - dcl_u*sina)
      dcz_u = -(dcd_u*sina + dcl_u*cosa)
      dcm_u =  -dcm_u

c---- negate signs since w here is +up, and standard CLa has a = w/V +down
      dcx_w = -(dcd_a*cosa - dcl_a*sina)
      dcz_w = -(dcd_a*sina + dcl_a*cosa)
      dcm_w =  -dcm_a

c---- add on extra imposed derivatives 
      cf_v(1,1) = cf_v(1,1) + dcx_u
      cf_v(1,3) = cf_v(1,3) + dcx_w

      cf_v(3,1) = cf_v(3,1) + dcz_u
      cf_v(3,3) = cf_v(3,3) + dcz_w

      cm_v(2,1) = cm_v(2,1) + dcm_u + dr(3)*dcx_u - dr(1)*dcz_u
      cm_v(2,3) = cm_v(2,3) + dcm_w + dr(3)*dcx_w - dr(1)*dcz_w

      return
      end ! gcoeffs


      SUBROUTINE GETFILE(LU,FNAME)
      CHARACTER(*) FNAME
C
      CHARACTER*1 ANS, DUMMY
C
 1000 FORMAT(A)
C
      IF(FNAME.EQ.' ') THEN
       CALL ASKS('Enter filename, or <return> for screen output^',FNAME)
      ENDIF
C
      IF(FNAME.EQ.' ') THEN
       LU = 6
       RETURN
C
      ELSE
       LU = 11
       OPEN(LU,FILE=FNAME,STATUS='OLD',FORM='FORMATTED',ERR=44)
       WRITE(*,*) 
       WRITE(*,*) 'File exists.  Append/Overwrite/Cancel  (A/O/C)?  C'
       READ(*,1000) ANS
       IF (INDEX('Aa',ANS).NE.0) THEN
C---- reopen file with append status HHY 4/17/18
         CLOSE(LU)
         OPEN(LU,FILE=FNAME,STATUS='OLD',FORM='FORMATTED',
     &        ACCESS='APPEND',ERR=44)
C---- old append code (re-reads to EOF)
cc 40     CONTINUE
cc        READ(LU,1000,END=42) DUMMY
cc        GO TO 40
cc 42     CONTINUE
       ELSEIF(INDEX('Oo',ANS).NE.0) THEN
         REWIND(LU) 
       ELSE
         CLOSE(LU)
         LU = 0
       ENDIF
       RETURN
C
 44    OPEN(LU,FILE=FNAME,STATUS='UNKNOWN',FORM='FORMATTED',ERR=48)
       REWIND(LU)
       RETURN
C
 48    CONTINUE
       LU = -1
       RETURN
      ENDIF
      END ! GETFILE
