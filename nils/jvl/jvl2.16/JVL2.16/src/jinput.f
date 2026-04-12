C***********************************************************************
C    Module:  jinput.f
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

      subroutine input(lun,fname1,ferr)
c---------------------------------------------------------
c     Reads and processes a jvl configuration input file
c---------------------------------------------------------
      include 'jvl.inc'

      character*(*) fname1
      logical ferr

      character*4  keywd
      character*80 cname, aname
      character*128 line
      logical ldupl, lhinge

      real clx(3), cdx(3)
      real xyzscal(3), xyztran(3)

      parameter (nwrk=nsmax, ibx=300)
      real xyzles(3,nwrk),chords(nwrk),aincs(nwrk),
     &     sspaces(nwrk),
     &     brady(nwrk), bradz(nwrk)

      real clafsrf, claf(nwrk)

      integer nspans(nwrk)

      integer nasrf, nasec(nwrk)
      real xasrf(ibx), xasec(ibx,nwrk), 
     &     sasrf(ibx), sasec(ibx,nwrk), 
     &     tasrf(ibx), tasec(ibx,nwrk)

      real hdiskf, hdisks(nwrk),
     &     fhjetf, fhjets(nwrk),
     &     djetf0, djets0(nwrk),
     &     djetf1, djets1(nwrk),
     &     djetf3, djets3(nwrk),
     &     dxdiskf,dxdisks(nwrk),
     &     dndiskf,dndisks(nwrk),
     &     tdiskf(3),tdisks(3,nwrk)

      real cdclsrf(6), cdclsec(6,nwrk)

      real xb(ibx), yb(ibx)
      real xin(ibx), yin(ibx), tin(ibx)
      real xbod(ibx), ybod(ibx), tbod(ibx), abod(ibx)

c---- max number of control, jet, or design variable lines per section
      parameter (iconx = 20)

      integer icontd(iconx,nwrk), nscon(nwrk),
     &        ijettd(iconx,nwrk), nsjet(nwrk),
     &        idestd(iconx,nwrk), nsdes(nwrk)

      real xhinged(iconx,nwrk), 
     &     vhinged(3,iconx,nwrk), 
     &     refld(iconx,nwrk),
     &     reflg(iconx,nwrk),
     &     reflj(iconx,nwrk),
     &     gaind(iconx,nwrk),
     &     gaing(iconx,nwrk),
     &     gainj(iconx,nwrk)

      real    rinput(10)
      integer iinput(10)
      logical error

      character*128 fname

c----------------------------------------------------
      ferr = .false.

      call bstrip(fname1,nfn1)

      fname = fname1
      nfn = nfn1
ccc      print *,'Reading file: ',fname
      open(unit=lun,file=fname,status='old',err=3)
      go to 6

 3    continue
      write(*,*) 
      write(*,*) '** OPEN error on file: ', fname(1:nfn)
      fname = fname1(1:nfn1) // '.jvl'
      nfn = nfn1+4
      write(*,*) '   Trying alternative: ', fname(1:nfn)
      open(unit=lun,file=fname,status='old',err=4)
      go to 6

 4    continue
      write(*,*) 
      write(*,*) '** OPEN error on file: ', fname(1:nfn)
      fname = fname1(1:nfn1) // '.avl'
      nfn = nfn1+4
      write(*,*) '   Trying alternative: ', fname(1:nfn)
      open(unit=lun,file=fname,status='old',err=5)
      go to 6

 5    continue
      write(*,*) 
      write(*,*) '** OPEN error on file: ', fname(1:nfn)
      ferr = .true.
      return

c----------------------------------------------------
 6    continue
      call bstrip(fname,nfn)
      write(*,*) 
      write(*,*) 'Reading file: ', fname(1:nfn), '  ...'

      dcl_a0 = 0.
      dcl_u0 = 0.
      dcl_ad0 = 0.

      dcd_a0 = 0.
      dcd_u0 = 0.
      dcd_ad0 = 0.

      dcm_a0 = 0.
      dcm_u0 = 0.
      dcm_ad0 = 0.

c---- initialize all entity counters
      nsec = 0

      nsurf = 0
      nstrip = 0
      nvor = 0

      nbody = 0
      nlbody = 0

      ncontrol = 0
      nvarjet = 0
      ndesign = 0

c---- initialize counters and active-entity indicators
      isurf = 0
      ibody = 0

c---- initialize input-file line counter
      iline = 0

c------------------------------------------------------------------------------
c---- start reading file

c---------------------------------------------------
      call rdline(lun,line,nline,iline)
      title = line(1:nline)

      write(*,1001) title(1:60)
 1001 format(/' Configuration: ', a)

c---------------------------------------------------
      call rdline(lun,line,nline,iline)
      read(line,*,err=990) mach0
      mach = mach0
      amach = 0.0

c---------------------------------------------------
      call rdline(lun,line,nline,iline)
c-----read symmetry inputs
 108  read(line,*,err=990) iysym, izsym, zsym
 109  continue
      
c---- y-symmetry plane hard-wired at y=0
      ysym = 0.

      if(iysym.gt.0) iysym =  1
      if(iysym.lt.0) iysym = -1
      if(izsym.gt.0) izsym =  1
      if(izsym.lt.0) izsym = -1

c---------------------------------------------------
      call rdline(lun,line,nline,iline)
      read(line,*,err=990) sref,cref,bref

      if(sref .le. 0.) sref = 1.
      if(cref .le. 0.) cref = 1.
      if(bref .le. 0.) bref = 1.

c---------------------------------------------------
      call rdline(lun,line,nline,iline)
      ninput = 3
      call getflt(line,rinput,ninput,error)
      if(error .or. ninput.lt.3) go to 990

      xyzref0(1) = rinput(1)
      xyzref0(2) = rinput(2)
      xyzref0(3) = rinput(3)

c---------------------------------------------------
c---- try to read cd data which may or may not be present
      call rdline(lun,line,nline,iline)
      read(line,*,err=8) cdref0

c---- drag data was read ok... just keep going
      go to 10

 8    continue
c---- read error occurred (drag data wasn't present)...
c      ... interpret the line as keyword
      cdref0 = 0.
      go to 11

c==============================================================================
c---- start of keyword-interpretation loop
 10   continue
      call rdline(lun,line,nline,iline)

 11   continue
      keywd = line(1:4)
      call touper(keywd)

c===========================================================================
      if    (keywd.eq.'CORE') then
       call rdline(lun,line,nline,iline)
       ninput = 3
       call getflt(line,rinput,ninput,error)
       if(error .or. ninput.lt.3) go to 990

       vrcorec = rinput(1)
       vrcorew = rinput(2)
       srcore  = rinput(3)
       
      elseif (keywd.eq.'EOF ') then
c------ end of file... clean up loose ends

        if(isurf.ne.0) then
c------- "old" surface is still active, so build it before finishing
         call makesurf(isurf, ibx,nsec, 
     &       nvc,cspace, nvs,sspace, nvcu,cspaceu, nvcw,cspacew,
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
         isurf = 0
        endif

        if(ibody.ne.0) then
c------- "old" body is still active, so build it before finishing

c------ check for body duplicated on symmetry plane
         if(ibody.ne.0 .and. ldupl .and.
     &      ydupl.eq.0.0 .and. xyztran(2).eq.0.0) then
          write(*,*)    '** Cannot duplicate body on symmetry plane'
          go to 990
         endif

         call makebody(ibody, ibx,
     &       nvb, bspace,
     &       xyzscal,xyztran,
     &       xbod,ybod,abod,nbod,
     &       rho,izimage,
     &       ldupl,ydupl)
         ibody = 0
        endif

c------ go finish up
        go to 900

c===========================================================================
      elseif(keywd.eq.'SURF') then
c------ new surface is about to start

        if(isurf.ne.0) then
c------- "old" surface is still active, so build it before starting new one
         call makesurf(isurf, ibx,nsec, 
     &       nvc,cspace, nvs,sspace, nvcu,cspaceu, nvcw,cspacew,
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
         isurf = 0
        endif

        if(ibody.ne.0) then
c------- "old" body is still active, so build it before finishing

c------ check for body duplicated on symmetry plane
          if(ibody.ne.0 .and. ldupl .and.
     &       ydupl.eq.0.0 .and. xyztran(2).eq.0.0) then
           write(*,*)    '** Cannot duplicate body on symmetry plane'
           go to 990
          endif

          call makebody(ibody, ibx,
     &       nvb, bspace,
     &       xyzscal,xyztran,
     &       xbod,ybod,abod,nbod,
     &       rho,izimage,
     &       ldupl,ydupl)
         ibody = 0
        endif

c------ new surface  (isurf.ne.0 denotes surface accumulation is active)
        nsurf = nsurf + 1
        isurf = min( nsurf , nsmax )

c------ default logical surface index is just the surface number
        lscomp(isurf) = isurf

c------ clear indices for accumulation
        nsec = 0
        isec = 0

c------ set surface defaults
        ydupl  = 0.0
        ldupl  = .false.
        lhinge = .false.

c------ assume this will be a conventional loaded surface
        lfwake(isurf) = .true.
        lfalbe(isurf) = .true.
        lfload(isurf) = .true.

c------ assume unity density factor
        rho = 1.0

c------ default z-image flag same as global flag
        izimage = izsym

c------ unity scales and zero offsets
        xyzscal(1) = 1.0
        xyzscal(2) = 1.0
        xyzscal(3) = 1.0
        xyztran(1) = 0.
        xyztran(2) = 0.
        xyztran(3) = 0.
        addinc = 0.

c------ dcl/da = 2 pi
        clafsrf = 1.0

c------ flat camberlines
        nasrf    = 2
        xasrf(1) = 0.
        xasrf(2) = 1.0
        sasrf(1) = 0.
        sasrf(2) = 0.
        tasrf(1) = 0.
        tasrf(2) = 0.

c------ no profile drag
        do l = 1, 6
          cdclsrf(l) = 0.
        enddo
ccc        lvisc  = .false.

c------ default jet parameters
        hdiskf = 0.
        fhjetf = 0.
        djetf0 = 0.
        djetf1 = 0.
        djetf3 = 0.
        dxdiskf = 0.
        dndiskf = 0.
        tdiskf(1) = -1.0
        tdiskf(2) =  0.
        tdiskf(3) =  0.
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c------ read surface name
        call rdline(lun,line,nline,iline)
        stitle(isurf) = line(1:nline)
        write(*,*)
        write(*,*) '  Building surface: ', stitle(isurf)

c------ read surface vortex-spacing parameters
        call rdline(lun,line,nline,iline)
        ninput = 8
        call getflt(line,rinput,ninput,error)
        if(error .or. ninput.lt.2) go to 990

        nvc = int( rinput(1) + 0.001 )
        cspace = rinput(2)

        if(ninput.ge.4) then
         nvs = int( rinput(3) + 0.001 )
         sspace = rinput(4)
        else
         nvs = 0
         sspace = 0.0
        endif

        if(ninput.lt.5) then
         nvcu = nvc/2
         nvcu = max( 1 , nvcu )
        else
         nvcu = int( rinput(5) + 0.5 )
        endif

        if(ninput.lt.6) then
c        cspaceu = cspace
         cspaceu = 0.0
        else
         cspaceu = rinput(6)
        endif

        if(ninput.lt.7) then
         nvcw = nvc
        else
         nvcw = int( rinput(7) + 0.5 )
         nvcw = max( 1 , nvcw )
        endif

        if(ninput.lt.8) then
c        cspacew = cspace
         cspacew = 0.0
        else
         cspacew = rinput(8)
        endif

c===========================================================================
      elseif(keywd.eq.'BODY') then
c------ new body is about to start

        if(isurf.ne.0) then
c------- "old" surface is still active, so build it before starting new one
         call makesurf(isurf, ibx,nsec, 
     &       nvc,cspace, nvs,sspace, nvcu,cspaceu, nvcw,cspacew,
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
         isurf = 0
        endif

        if(ibody.ne.0) then
c------- "old" body is still active, so build it before finishing

c------ check for body duplicated on symmetry plane
         if(ibody.ne.0 .and. ldupl .and.
     &      ydupl.eq.0.0 .and. xyztran(2).eq.0.0) then
          write(*,*)    '** Cannot duplicate body on symmetry plane'
          go to 990
         endif

         call makebody(ibody, ibx,
     &       nvb, bspace,
     &       xyzscal,xyztran,
     &       xbod,ybod,abod,nbod,
     &       rho,izimage,
     &       ldupl,ydupl)
         ibody = 0
        endif

c------ new body  (ibody.ne.0 denotes body accumulation is active)
        nbody = nbody + 1
        if(nbody.gt.nbmax) then
         write(*,*) 'Body array overflow. Increase nbmax.'
         stop
        endif
        ibody = nbody

c------ default logical body index is just minus the body number
        lbcomp(ibody) = -ibody

        nsec = 0
        isec = 0
        nbod = 0

        nin = 0

c------ assume this will be a conventional loaded body
        lbwake(ibody) = .true.
        lbalbe(ibody) = .true.
        lbload(ibody) = .true.

        ydupl  = 0.0
        ldupl  = .false.
        lhinge = .false.
ccc     ljet   = .false.

c------ assume unity density factor
        rho = 1.0

c------ default image flag same as global flag
        izimage = izsym

c------ unity scales and zero offsets
        xyzscal(1) = 1.0
        xyzscal(2) = 1.0
        xyzscal(3) = 1.0
        xyztran(1) = 0.
        xyztran(2) = 0.
        xyztran(3) = 0.

        call rdline(lun,line,nline,iline)
        btitle(ibody) = line(1:nline)
        write(*,*)
        write(*,*) '  Building body: ', btitle(ibody)

        call rdline(lun,line,nline,iline)
        read(line,*,err=990) nvb, bspace

c===========================================================================
      elseif(keywd.eq.'YDUP') then
c------ this surface is to be duplicated with an image surface
        if    (isurf.ne.0) then
cc         write(*,*) '  + duplicate surface ',stitle(isurf)
        elseif(ibody.ne.0) then
cc         write(*,*) '  + duplicate body ',btitle(ibody)
        else
         write(*,9000) '** Misplaced line', iline, line(1:nline)
         write(*,*)    '** No active surface or body to duplicate'
         go to 10
        endif
        
        call rdline(lun,line,nline,iline)
        read(line,*,err=990) ydupl

        ldupl = .true.

        if(iysym.ne.0 .and. ydupl.eq.0.0) then
         write(*,*) 'ERROR: Redundant y-symmetry specifications...'
         write(*,*) '       IYSYM /= 0'
         write(*,*) '       YDUPLICATE  0.0'
         write(*,*) 'Can use one or the other, but not both.'
         stop
        endif

c===========================================================================
      elseif (keywd.eq.'INDE' .or. keywd.eq.'COMP') then
c------ set surface index
        if(isurf.eq.0) then
         write(*,9000) '** Misplaced line', iline, line(1:nline)
         write(*,*)    '** No active surface for index'
         go to 10
        endif

        call rdline(lun,line,nline,iline)
        read(line,*,err=990) lscomp(isurf)

c===========================================================================
      elseif (keywd.eq.'SCAL') then
c------ read scaling factors
        if(isurf.eq.0 .and. ibody.eq.0) then
         write(*,9000) '** Misplaced line', iline, line(1:nline)
         write(*,*)    '** No active surface or body for scaling'
         go to 10
        endif

        call rdline(lun,line,nline,iline)
        read(line,*,err=990) xyzscal(1), xyzscal(2), xyzscal(3)

c===========================================================================
      elseif (keywd.eq.'TRAN') then
c------ read translation vector
        if(isurf.eq.0 .and. ibody.eq.0) then
         write(*,9000) '** Misplaced line', iline, line(1:nline)
         write(*,*)    '** No active surface or body for translation'
         go to 10
        endif

        call rdline(lun,line,nline,iline)
        read(line,*,err=990) xyztran(1), xyztran(2), xyztran(3)

c===========================================================================
      elseif (keywd.eq.'ANGL' .or. keywd.eq.'AINC') then
c------ read surface angle change
        if(isurf.eq.0 .and. ibody.eq.0) then
         write(*,9000) '** Misplaced line', iline, line(1:nline)
         write(*,*)    '** No active surface or body for rotation'
         go to 10
        endif

        call rdline(lun,line,nline,iline)
        read(line,*,err=990) addinc
        addinc = addinc*dtr

c===========================================================================
      elseif (keywd.eq.'DENS') then
c------ read density factor
        if(isurf.eq.0 .and. ibody.eq.0) then
         write(*,9000) '** Misplaced line', iline, line(1:nline)
         write(*,*)    '** No active surface or body for density factor'
         go to 10
        endif

        call rdline(lun,line,nline,iline)
        read(line,*,err=990) rho

c===========================================================================
      elseif (keywd.eq.'ZIMA') then
c------ read image flag
        if(isurf.eq.0 .and. ibody.eq.0) then
         write(*,9000) '** Misplaced line', iline, line(1:nline)
         write(*,*)    '** No active surface or body for z-image flag'
         go to 10
        endif

        call rdline(lun,line,nline,iline)
        read(line,*,err=990) izimage
        if(izimage .gt. 0) izimage =  1
        if(izimage .lt. 0) izimage = -1

c===========================================================================
      elseif (keywd.eq.'NOWA') then
c------ disable wake shedding for this surface
        if(isurf.eq.0) then
         write(*,9000)'** Misplaced line', iline, line(1:nline)
         write(*,*)   '** No active surface for wake-shedding flag'
         go to 10
        endif

        lfwake(isurf) = .false.

c===========================================================================
      elseif (keywd.eq.'NOAL') then
c------ disable freestream angles for this surface
        if(isurf.eq.0) then
         write(*,9000)'** Misplaced line', iline, line(1:nline)
         write(*,*)   '** No active surface for freestream-angles flag'
         go to 10
        endif

        lfalbe(isurf) = .false.

c===========================================================================
      elseif (keywd.eq.'NOLO') then
c------ disable total-load contributions for this surface
        if(isurf.eq.0) then
         write(*,9000)'** Misplaced line', iline, line(1:nline)
         write(*,*)   '** No active surface for load-disable flag'
         go to 10
        endif

        lfload(isurf) = .false.

c===========================================================================
c      elseif (keywd.eq.'BSEC') then
cc------ read body section
c        if(ibody.eq.0) then
c         write(*,9000) '** Misplaced line', iline, line(1:nline)
c         write(*,*)    '** No active body for this section'
c         go to 10
c        endif
c
c        call rdline(lun,line,nline,iline)
c
cc------ store section data for current body
c        nsec = nsec + 1
c        isec = min( nsec , nwrk )
c
c        ninput = 5
c        call getflt(line,rinput,ninput,error)
c        if(error .or. ninput.lt.4) go to 990
cc
c
c        xyzles(1,isec) = rinput(1)
c        xyzles(2,isec) = rinput(2)
c        xyzles(3,isec) = rinput(3)
c        brady(isec) = rinput(4)
c
c        if(ninput.ge.5) then
c         bradz(isec) = rinput(5)
c        else
c         bradz(isec) = brady(isec)
c        endif
c
cc        nbod = isec
cc        i = nbod
c
cc        call getcam(xb,yb,nb,xbod,ybod,tbod,nbod,.false.)
cc        do i = 1, nbod
cc          abod(i) = 0.25*pi*tbod(i)**2
cc        enddo

c===========================================================================
      elseif (keywd.eq.'SECT') then
c------ read surface section
        if(isurf.eq.0) then
         write(*,9000) '** Misplaced line', iline, line(1:nline)
         write(*,*)    '** No active surface for this section'
         go to 10
        endif

        call rdline(lun,line,nline,iline)

c------ store section data for current surface
        nsec = nsec + 1
        isec = min( nsec , nwrk )

        ninput = 7
        call getflt(line,rinput,ninput,error)
        if(error .or. ninput.lt.5) go to 990

        xyzles(1,isec) = rinput(1)
        xyzles(2,isec) = rinput(2)
        xyzles(3,isec) = rinput(3)
        chords(isec) = rinput(4)
        aincs(isec)  = rinput(5)*dtr

        if(ninput.ge.7) then
         nspans(isec) = int( rinput(6) + 0.001 )
         sspaces(isec) = rinput(7)
        else
         nspans(isec) = 0
         sspaces(isec) = 0.
        endif

c------ default section parameters are those for entire surface...

c   ... no control variables
        nscon(isec) = 0

c   ... no jet variables
        nsjet(isec) = 0

c   ... no design variables
        nsdes(isec) = 0

c   ... camberline
        nasec(isec) = nasrf
        do ia = 1, nasrf
          xasec(ia,isec) = xasrf(ia) 
          sasec(ia,isec) = sasrf(ia) 
          tasec(ia,isec) = tasrf(ia) 
        enddo

c   ... dcl/da factor
        claf(isec) = clafsrf

c   ... polar data
        do l = 1, 6
          cdclsec(l,isec) = cdclsrf(l)
        end do

c   ... jet parameters
        hdisks(isec) = hdiskf
        fhjets(isec) = fhjetf
        djets0(isec) = djetf0
        djets1(isec) = djetf1
        djets3(isec) = djetf3
        dxdisks(isec) = dxdiskf
        dndisks(isec) = dndiskf
        tdisks(1,isec) = tdiskf(1)
        tdisks(2,isec) = tdiskf(2)
        tdisks(3,isec) = tdiskf(3)
c===========================================================================
      elseif (keywd.eq.'NACA') then 
c------ input naca camberline
        if(isurf.eq.0) then
         write(*,9000) '** Misplaced line', iline, line(1:nline)
         write(*,*)    '** No active surface for this airfoil'
         go to 10
        endif

        call rdline(lun,line,nline,iline)

        ib = index(line,' ')
        read(line(1:ib-1),*,err=990) ides
        if(line(ib:nline).ne.' ') then
         read(line(ib:nline),*,err=990) xfmin, xfmax
ccc           write(*,*) '   Using data in normalized range ',xfmin,xfmax
        else
         xfmin = 0.0
         xfmax = 1.0
        endif

        icam = ides/1000
        ipos = (ides-1000*icam)/100
        ithk =  ides-1000*icam-100*ipos
        c = float(icam) / 100.0
        p = float(ipos) / 10.0
        t = float(ithk) / 100.0


        napts = min( 50 , ibx )
        do i = 1, napts
          frac = float(i-1)/float(napts-1)
          xf = xfmin*(1.0-frac) + xfmax*frac
          if(xf.lt.p) then
           slp = c/p**2 * 2.0*(p - xf)
          else
           slp = c/(1.0-p)**2 * 2.0*(p - xf)
          endif
          thk = (0.29690*sqrt(xf)
     &          - 0.12600*xf
     &          - 0.35160*xf**2
     &          + 0.28430*xf**3
     &          - 0.10150*xf**4) * t / 0.10

          if(isec.eq.0) then
c--------- no SECTION active yet... store camberline for entire surface
           xasrf(i) = xf
           sasrf(i) = slp
           tasrf(i) = thk
          else
c--------- store for this SECTION only
           xasec(i,isec) = xf
           sasec(i,isec) = slp
           tasec(i,isec) = thk
          endif
        enddo

        if(isec.eq.0) then
         nasrf = napts
         call nrmliz(nasrf,xasrf(1))
        else
         nasec(isec) = napts
         call nrmliz(nasec(isec),xasec(1,isec))
        endif

c===========================================================================
      else if (keywd.eq.'AIRF') then 
c------ input y(x) for an airfoil, get camber then slopes via spline
        if(isurf.eq.0) then
         write(*,9000) '** Misplaced line', iline, line(1:nline)
         write(*,*)    '** No active surface for this airfoil'
         go to 10
        endif

        call rdline(lun,line,nline,iline)

        ib = index(line,' ')
        read(line(ib:nline),*,err=990) xfmin, xfmax

        do i = 1, 999999
          ib = min(i,ibx)

          call rdline(lun,line,nline,iline)
          ninput = 2
          call getflt(line,rinput,ninput,error)
          if(error .or. ninput.lt.2) then
           nb = ib-1
           go to 40
          else
           xb(ib) = rinput(1)
           yb(ib) = rinput(2)
          endif
        enddo

 40     continue
        if(i.gt.ibx) then
         write(*,*) 
     &    '*** INPUT: Airfoil array overflow.  Increase ibx to', i
         stop
        endif

c------ set camber and thickness, normalized to unit chord
        nin = min( 50 , ibx )
        call getcam(xb,yb,nb,xin,yin,tin,nin,.true.)

        if(isec.eq.0) then
c------- no section active yet... store camberline for the entire surface
         nasrf = nin
         do i = 1, nin
           xf = xfmin + (xfmax-xfmin)*float(i-1)/float(nasrf-1)
           xasrf(i) = xin(1) + xf*(xin(nin)-xin(1))
           call akima(xin,yin,nin,xasrf(i),dummy,sasrf(i))
           call akima(xin,tin,nin,xasrf(i),tasrf(i),dummy)
         end do
         call nrmliz(nasrf,xasrf(1))

        else
c------- store camberline for this section
         nasec(isec) = nin
         do i = 1, nin
           xf = xfmin + (xfmax-xfmin)*float(i-1)/float(nasec(isec)-1)
           xasec(i,isec) = xin(1) + xf*(xin(nin)-xin(1))
           call akima(xin,yin,nin,xasec(i,isec),dummy,sasec(i,isec))
           call akima(xin,tin,nin,xasec(i,isec),tasec(i,isec),dummy)
         end do
         call nrmliz(nasec(isec),xasec(1,isec))

        endif

c------ go to top of keyword-reading loop, with last-read line
        go to 11

c===========================================================================
      elseif (keywd.eq.'AFIL') then 
c------ input y(x) camberline from an airfoil coordinate file
        if(isurf.eq.0) then
         write(*,9000) '** Misplaced line', iline, line(1:nline)
         write(*,*)    '** No active surface for this airfoil'
         go to 10
        endif

        call rdline(lun,line,nline,iline)

        ib = index(line,' ')
        cname = line(1:ib)

        if(line(ib:nline).ne.' ') then
         read(line(ib:nline),*,err=990) xfmin, xfmax
ccc         write(*,*) '     Using data in normalized range ',xfmin,xfmax
        else
         xfmin = 0.
         xfmax = 1.
        endif

        call bstrip(cname,ncn)
        write(*,*) '    Reading airfoil from file: ',cname(1:ncn)
        nblds = 1
        call readbl(cname,ibx,nblds,xb,yb,nb,nbl,
     &               aname,xinl,xout,ybot,ytop)

        if(nbl.eq.0) then
         write(*,*) '**   Airfoil file not found  : ',cname(1:ncn)
         write(*,*) '**   Using default camberline'

        else
c------- camber and thickness
         nin = min( 50 , ibx )
         call getcam(xb,yb,nb,xin,yin,tin,nin,.true.)

         if(isec.eq.0) then
c-------- no SECTION active yet... set camberline for entire surface
          nasrf = nin
          do i = 1, nin
            xf = xfmin + (xfmax-xfmin)*float(i-1)/float(nasrf-1)
            xasrf(i) = xin(1) + xf*(xin(nin)-xin(1))
            call akima(xin,yin,nin,xasrf(i),dummy,sasrf(i))
            call akima(xin,tin,nin,xasrf(i),tasrf(i),dummy)
          end do
          call nrmliz (nasrf,xasrf(1))

         else
c-------- store camberline for this SECTION
          nasec(isec) = nin
          do i = 1, nin
            xf = xfmin + (xfmax-xfmin)*float(i-1)/float(nasec(isec)-1)
            xasec(i,isec) = xin(1) + xf*(xin(nin)-xin(1))
            call akima(xin,yin,nin,xasec(i,isec),dummy,sasec(i,isec))
            call akima(xin,tin,nin,xasec(i,isec),tasec(i,isec),dummy)
          end do
          call nrmliz (nasec(isec),xasec(1,isec))

         endif

        endif

c===========================================================================
      elseif (keywd.eq.'BFIL') then 
c------ input body r(x) from an airfoil coordinate file
        if(ibody.eq.0) then
         write(*,9000) '** Misplaced line', iline, line(1:nline)
         write(*,*)    '** No active body for this shape'
         go to 10
        endif

        call rdline(lun,line,nline,iline)

        ib = index(line,' ')
        cname = line(1:ib)
        if(line(ib:nline).ne.' ') then
         read(line(ib:nline),*,err=990) xfmin, xfmax
ccc         write(*,*) '     using data in normalized range ',xfmin,xfmax
        else
         xfmin = 0.
         xfmax = 1.
        endif

        call bstrip(cname,ncn)
        write(*,*) '    Reading body shape from file: ',cname(1:ncn)
        nblds = 1
        call readbl(cname,ibx,nblds,xb,yb,nb,nbl,
     &               aname,xinl,xout,ybot,ytop)

c------ set thread line y, and thickness t ( = 2r)
        nbod = min( 50 , ibx )
        call getcam(xb,yb,nb,xbod,ybod,tbod,nbod,.false.)
        do i = 1, nbod
          abod(i) = 0.25*pi*tbod(i)**2
        enddo

c===========================================================================
      elseif (keywd.eq.'CDCL') then 
c------ input approximate cd(cl) polar defining data
        IF(ISURF.EQ.0 .AND. ISEC.EQ.0) THEN
         write(*,9000) '** Misplaced line', iline, line(1:nline)
         write(*,*)    '** No active surface for this polar'
         go to 10
        endif

C------ if defining surface before sections store polar for surface
        IF(ISURF.NE.0 .AND. ISEC.EQ.0) THEN
         WRITE(*,*)    '** Polar data for surface (all sections)'

        CALL RDLINE(LUN,LINE,NLINE,ILINE)
        READ(LINE,*,ERR=990) CLX(1),CDX(1),CLX(2),CDX(2),CLX(3),CDX(3)
C
        LMAX = 1
        LMIN = 1
        DO L = 2, 3
          IF(CLX(L).GT.CLX(LMAX)) LMAX = L
          IF(CLX(L).LT.CLX(LMIN)) LMIN = L
        END DO
C
C------ Trick: sum must be 6 so we can get the "other" index
         LMID = 6 - (LMIN+LMAX)
         CDCLSRF(1) = CLX(LMIN)
         CDCLSRF(2) = CDX(LMIN)
         CDCLSRF(3) = CLX(LMID)
         CDCLSRF(4) = CDX(LMID)
         CDCLSRF(5) = CLX(LMAX)
         CDCLSRF(6) = CDX(LMAX)
         WRITE(*,1700) CLX(LMIN),CDX(LMIN),
     &                 CLX(LMID),CDX(LMID),
     &                 CLX(LMAX),CDX(LMAX)
 1700    FORMAT('    Reading CD(CL) data for surface',
     &         /'     CLneg    = ',F8.3,'  CD@CLneg = ',F10.5,
     &         /'     CL@CDmin = ',F8.3,'  CDmin    = ',F10.5,
     &         /'     CLpos    = ',F8.3,'  CD@CLpos = ',F10.5)
         LVISC = .TRUE.
C
C------ define polar for the current active section
        ELSEIF(ISURF.NE.0 .AND. ISEC.NE.0) THEN
C
         CALL RDLINE(LUN,LINE,NLINE,ILINE)
         READ(LINE,*,ERR=990) CLX(1),CDX(1),CLX(2),CDX(2),CLX(3),CDX(3)
C
         LMAX = 1
         LMIN = 1
         DO L = 2, 3
           IF(CLX(L).GT.CLX(LMAX)) LMAX = L
           IF(CLX(L).LT.CLX(LMIN)) LMIN = L
         END DO
C
         IF(ISEC.GT.1) THEN
          IF(CDCLSEC(4,ISEC-1).LE.0.0) THEN
           WRITE(*,*) '* AINPUT: previous section defined with no polar'
          ENDIF
         ENDIF
C
C------ Trick: sum must be 6 so we can get the "other" index
         LMID = 6 - (LMIN+LMAX)
         CDCLSEC(1,ISEC) = CLX(LMIN)
         CDCLSEC(2,ISEC) = CDX(LMIN)
         CDCLSEC(3,ISEC) = CLX(LMID)
         CDCLSEC(4,ISEC) = CDX(LMID)
         CDCLSEC(5,ISEC) = CLX(LMAX)
         CDCLSEC(6,ISEC) = CDX(LMAX)
         WRITE(*,1701) CLX(LMIN),CDX(LMIN),
     &                 CLX(LMID),CDX(LMID),
     &                 CLX(LMAX),CDX(LMAX)
 1701    FORMAT('    Reading CD(CL) data for section',
     &         /'     CLneg    = ',F8.3,'  CD@CLneg = ',F10.5,
     &         /'     CL@CDmin = ',F8.3,'  CDmin    = ',F10.5,
     &         /'     CLpos    = ',F8.3,'  CD@CLpos = ',F10.5)

         LVISC = .TRUE.
       ENDIF

c===========================================================================
      elseif (keywd.eq.'CLAF') then 
c------ input dcl/da scaling factor
        if(isurf.eq.0) then
         write(*,9000) '** Misplaced line', iline, line(1:nline)
         write(*,*)    '** No active surface for dcl/da factor'
         go to 10
        endif

        call rdline(lun,line,nline,iline)
        read(line,*,err=990) claf1

        if(claf1 .le. 0.0 .or. 
     &     claf1 .ge. 2.0      ) then
         write(*,*) '** dcl/da factor must be in the range 0..2 --',
     &              ' Setting factor to 1.0'
         claf1 = 1.0
        endif

        if(isec.eq.0) then
         clafsrf = claf1
        else
         claf(isec) = claf1
        endif

c===========================================================================
      elseif (keywd.eq.'JETP') then 
c------ input jet parameters
        if(isurf.eq.0) then
         write(*,9000) '** Misplaced line', iline, line(1:nline)
         write(*,*)    '** No active surface for jet parameters'
         go to 10
        endif

        call rdline(lun,line,nline,iline)

        ninput = 10
        call getflt(line,rinput,ninput,error)
        if(error .or. ninput.lt.1) then
         write(*,*) '*** Bad jet parameter data line:  ', line
         stop
        endif

        if(isec.eq.0) then
c------- no SECTION defined yet... store jet parameters for entire surface
         if(ninput .ge. 1) hdiskf = rinput(1)
         if(ninput .ge. 2) fhjetf = rinput(2)
         if(ninput .ge. 3) djetf0 = rinput(3)
         if(ninput .ge. 4) djetf1 = rinput(4)
         if(ninput .ge. 5) djetf3 = rinput(5)
         if(ninput .ge. 6) dxdiskf = rinput(6)
         if(ninput .ge. 7) dndiskf = rinput(7)
         if(ninput .ge. 8 ) tdiskf(1) = rinput(8)
         if(ninput .ge. 9 ) tdiskf(2) = rinput(9)
         if(ninput .ge. 10) tdiskf(3) = rinput(10)
        else
c------- store jet parameters for this SECTION only
         if(ninput .ge. 1) hdisks(isec) = rinput(1)
         if(ninput .ge. 2) fhjets(isec) = rinput(2)
         if(ninput .ge. 3) djets0(isec) = rinput(3)
         if(ninput .ge. 4) djets1(isec) = rinput(4)
         if(ninput .ge. 5) djets3(isec) = rinput(5)
         if(ninput .ge. 6) dxdisks(isec) = rinput(6)
         if(ninput .ge. 7) dndisks(isec) = rinput(7)
         if(ninput .ge. 8 ) tdisks(1,isec) = rinput(8)
         if(ninput .ge. 9 ) tdisks(2,isec) = rinput(9)
         if(ninput .ge. 10) tdisks(3,isec) = rinput(10)
         if(hdisks(isec) .lt. 0.0) then
          write(*,*) '* hdisk < 0 specified'
         endif
        endif

c===========================================================================
      elseif (keywd.eq.'CONT') then
c------ link section to control variables
        if(isurf.eq.0 .or. isec.eq.0) then
         write(*,9000) '** Misplaced line', iline, line(1:nline)
         write(*,*)    '** No active section for this control'
         go to 10
        endif

        call rdline(lun,line,nline,iline)

c------ increment control-declaration counter for this section
        nscon(isec) = nscon(isec) + 1
        iscon = min( nscon(isec) , iconx )

c------ extract control name
        nname = index(line,' ') - 1
        if(nname.le.0) then
         write(*,*) '** Bad control declaration line:  ', line
         stop
        endif

c------ see if this control variable has already been declared
        do n = 1, ncontrol
          ndname = index(dname(n),' ') - 1      ! added 7 Dec 2021
          if(nname .eq. ndname .and.            ! added 7 Dec 2021
     &      line(1:nname) .eq. dname(n)(1:nname)) then   ! modified 7 Dec 2021
ccc      if(line(1:nname) .eq. dname(n)(1:nname)) then
           icontrol = n
           go to 62
          endif
        enddo

c------ new control variable... assign slot for it
        ncontrol = ncontrol + 1
        icontrol = min( ncontrol , ndmax )
        dname(icontrol) = line(1:nname)
        gconmax(icontrol) = 0.

 62     continue
        icontd(iscon,isec) = icontrol

c------ read numbers after control variable name
        ninput = 6
        call getflt(line(nname+1:120),rinput,ninput,error)
        if(error) then
         write(*,*) '*** Bad control data line:  ', line
         stop
        endif

        if(ninput.lt.1) then
         gaind(iscon,isec) = 1.0
        else
         gaind(iscon,isec) = rinput(1)
        endif

        if(ninput.lt.2) then
         xhinged(iscon,isec) = 0.0
        else
         xhinged(iscon,isec) = rinput(2)
        endif

        if(ninput.lt.5) then
         vhinged(1,iscon,isec) = 0.0
         vhinged(2,iscon,isec) = 0.0
         vhinged(3,iscon,isec) = 0.0
        else
         vhinged(1,iscon,isec) = rinput(3)
         vhinged(2,iscon,isec) = rinput(4)
         vhinged(3,iscon,isec) = rinput(5)
        endif

        if(ninput.lt.6) then
         refld(iscon,isec) = 1.0
        else
         refld(iscon,isec) = rinput(6)
        endif

c===========================================================================
      elseif (keywd.eq.'JETC') then
c------ link section to jet control variables
        if(isurf.eq.0 .or. isec.eq.0) then
         write(*,9000) '** Misplaced line', iline, line(1:nline)
         write(*,*)    '** No active section for this jet control'
         go to 10
        endif

        call rdline(lun,line,nline,iline)

c------ increment jet control-declaration counter for this section
        nsjet(isec) = nsjet(isec) + 1
        isjet = min( nsjet(isec) , iconx )

c------ extract jet control name
        nname = index(line,' ') - 1
        if(nname.le.0) then
         write(*,*) '** Bad jet control declaration line:  ', line
         stop
        endif

c------ see if this jet control variable has already been declared
        do n = 1, nvarjet
          if(line(1:nname) .eq. jname(n)(1:nname)) then
           ivarjet = n
           go to 72
          endif
        enddo

c------ new jet control variable... assign slot for it
        nvarjet = nvarjet + 1
        ivarjet = min( nvarjet , njmax )
        jname(ivarjet) = line(1:nname)
        gjetmax(ivarjet) = 0.

 72     continue
c------ jet variable index for this JETC line 
        ijettd(isjet,isec) = ivarjet

c------ read numbers after jet control variable name
        ninput = 2
        call getflt(line(nname+1:120),rinput,ninput,error)
        if(error) then
         write(*,*) '*** Bad control data line:  ', line
         stop
        endif

        if(ninput.lt.1) then
         gainj(isjet,isec) = 1.0
        else
         gainj(isjet,isec) = rinput(1)
        endif

        if(ninput.lt.2) then
         reflj(isjet,isec) = 1.0
        else
         reflj(isjet,isec) = rinput(2)
        endif

c===========================================================================
      elseif (keywd.eq.'DESI') then 
c------ link section to design variable and weight

        if(isurf.eq.0 .or. isec.eq.0) then
         write(*,9000) '** Misplaced line', iline, line(1:nline)
         write(*,*)    '** No active section for this design var.'
         go to 10
        endif

        call rdline(lun,line,nline,iline)

c------ increment design-declaration counter for this section
        nsdes(isec) = nsdes(isec) + 1
        isdes = min( nsdes(isec) , iconx )

c------ extract design name
        nname = index(line,' ') - 1
        if(nname.le.0) then
         write(*,9000) '   *** Bad design declaration line', 
     &                  iline, line(1:nline)
         stop
        endif

c------ see if this control variable has already been declared
        do k = 1, ndesign
          if(line(1:nname) .eq. gname(k)(1:nname)) then
           idesign = k
           go to 82
          endif
        enddo

        ndesign = ndesign + 1
        idesign = min( ndesign , ngmax )
        gname(idesign) = line(1:nname)

 82     continue
        idestd(isdes,isec) = idesign

c------ read numbers after control variable name
        ninput = 1
        call getflt(line(nname+1:120),rinput,ninput,error)
        if(error) go to 990

        if(ninput.lt.1) then
         gaing(isdes,isec) = 1.0
        else
         gaing(isdes,isec) = rinput(1)
        endif

c===========================================================================
      else
c------ line not recognized or unassignable ... keep reading file
        write(*,8000) iline, line(1:nline)
 8000   format('  * Line',i5,' ignored: ', a)

      endif
      go to 10

c===========================================================================
c---- normal end-of-file exit point
 900  continue
      close(unit=lun)

      if(nsurf .gt. nsmax) then
       write(*,*) 'INPUT: Array overflow. Increase nsmax to ', nsurf
       stop
      endif

      if(nbody .gt. nbmax) then
       write(*,*) 'INPUT: Array overflow. Increase nbmax to ', nbody
       stop
      endif

c===================================================================
      njets = 0
      nvoru = 0
      nvorw = 0
      do js = 1, nstrip
        if(ifrstw(js) .gt. 0) then
         njets = njets + 1
         nvorw = nvorw + (ilastu(js) - ifrstu(js) + 1)
     &                 + (ilastw(js) - ifrstw(js) + 1)
        endif
      enddo

      nvorj = nvoru + nvorw
      write (*,2018) mach,nbody,nsurf,nstrip,njets,nvor,nvorj
      write (*,2019) ncontrol,nvarjet,ndesign

      if(iysym.gt.0) write (*,2024) ysym
      if(iysym.lt.0) write (*,2025) ysym
      if(izsym.gt.0) write (*,2026) zsym
      if(izsym.lt.0) write (*,2027) zsym

 2018 format (/' Mach =',f10.4,'  (default)'
     &       //1x,i4,' Bodies'
     &        /1x,i4,' Surfaces'
     &        /1x,i4,' Strips,   ', i4,' of these have Jet Strips'
     &        /1x,i4,' Vortices, ', i4,' of these are on Jet Strips')
 2019 format (/1x,i4,' Control variables'
     &        /1x,i4,' Jet variables'
     &        /1x,i4,' Design parameters')

 2024 format (/' Y symmetry: Wall plane   at Ysym =',f10.4)
 2025 format (/' Y symmetry: Free surface at Ysym =',f10.4)
 2026 format (/' Z symmetry: Ground plane at Zsym =',f10.4)
 2027 format (/' Z symmetry: Free surface at Zsym =',f10.4)

      lgeo = .true.
      return

c*********************************************************************
 990  continue
      write(*,9000) '** READ error on line', iline, line(1:nline)
      ferr = .true.
      return

 9000 format(/ 1x,a,i5,' ...' / 1x,a)
      end ! input


      subroutine rdline(lun,line,nline,iline)
c-----------------------------------------------------------------------
c     Reads next non-comment line from logical unit lun
c     Strips off leading blanks
c     Ignores everything after and including "!"
c
c     line   returns the line
c     nline  returns the number of characters in non-blank portion
c
c     If e.o.f. is reached, line returns 'EOF'
c     If read error occurs, line returns 'ERR'
c-----------------------------------------------------------------------
      character*(*) line

   20 continue
      read(lun,1000,end=80,err=90) line
 1000 format(a)
      iline = iline + 1

c---- skip comment line
      if(index('!#',line(1:1)) .ne. 0) go to 20

c---- skip blank line
      if(line.eq.' ') go to 20

c---- strip off leading blanks and do normal return after significant line
      call bstrip(line,nline)
      kexl = index(line(1:nline),'!')
      if(kexl.gt.1) nline = kexl-1
      return

   80 line = 'EOF '
      return

   90 line = 'ERR '
      return
      end ! rdline


      subroutine touper(input)
      character*(*) input

      character*26 lcase, ucase
      data lcase / 'abcdefghijklmnopqrstuvwxyz' /
      data ucase / 'ABCDEFGHIJKLMNOPQRSTUVWXYZ' /

      n = len(input)

      do i=1, n
        k = index( lcase , input(i:i) )
        if(k.gt.0) input(i:i) = ucase(k:k)
      enddo

      return
      end ! touper

