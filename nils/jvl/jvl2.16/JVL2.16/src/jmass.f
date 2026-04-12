C***********************************************************************
C    Module:  jmass.f
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
 
      subroutine masini
      include 'jvl.inc'
c---------------------------------------------------
c     initializes default mass, inertia tensor
c---------------------------------------------------
      rmass0 = 1.0

      do k = 1, 3
        riner0(k,1) = 0.
        riner0(k,2) = 0.
        riner0(k,3) = 0.
        riner0(k,k) = 1.0

        amass(k,1) = 0.
        amass(k,2) = 0.
        amass(k,3) = 0.

        ainer(k,1) = 0.
        ainer(k,2) = 0.
        ainer(k,3) = 0.

        xyzmass0(k) = 0.
        hrotor0(k) = 0.
      enddo

      lmass = .false.

      return
      end ! masini



      subroutine masget(lu,fname,error)
c--------------------------------------------------
c     reads mass distributions file,
c     computes default mass, inertia tensor
c--------------------------------------------------
      include 'jvl.inc'
      character*(*) fname

      real mass, mi
      real ixx , iyy , izz , ixy , ixz , iyz
      real ixxi, iyyi, izzi, ixyi, ixzi, iyzi

      character*80 case
      character*256 line, lineu

      character*32 unchgee, unchrho, unch

      logical error
      parameter (nidim=13)
      real rinp(nidim), fac(nidim), add(nidim)

 1000 format(a)

      open(lu,file=fname,status='old',err=90)

      sum_m   = 0.

      sum_mx  = 0.
      sum_my  = 0.
      sum_mz  = 0.

      sum_mxx = 0.
      sum_myy = 0.
      sum_mzz = 0.
      sum_mxy = 0.
      sum_mxz = 0.
      sum_myz = 0.
      sum_mzz = 0.

      sum_ixx = 0.
      sum_iyy = 0.
      sum_izz = 0.
      sum_ixy = 0.
      sum_iyz = 0.
      sum_ixz = 0.

      sum_hx  = 0.
      sum_hy  = 0.
      sum_hz  = 0.

      xcg = 0.
      ycg = 0.
      zcg = 0.

c---- default multipliers and adders
      do k = 1, nidim
        fac(k) = 1.0
        add(k) = 0.0
      enddo

      unitl = 1.
      unitm = 1.
      unitt = 1.
      unchl = 'Lunit'
      unchm = 'Munit'
      uncht = 'Tunit'
      nul = 5
      num = 5
      nut = 5

      gee0 = 1.
      rho0 = 1.
      unchgee = 'Lunit/Tunit^2'
      unchrho = 'Munit/Lunit^3'
      nugee = 13
      nurho = 13

      iline = 0

c---- search all lines of file for data values
c=============================================
 10   continue
        read(lu,1000,end=50) line
        iline = iline + 1

        if(index('#!',line(1:1)).ne.0) go to 10

        if(index('*',line(1:1)).ne.0) then
c------- read multiplier line
         kinp = nidim
         call getflt(line(2:80),rinp,kinp,error)
         if(error) go to 40

         do k = 1, kinp
           fac(k) = rinp(k)
         enddo
         go to 10
        endif

        if(index('+',line(1:1)).ne.0) then
c------- read adder line
         kinp = nidim
         call getflt(line(2:80),rinp,kinp,error)
         if(error) go to 40

         do k = 1, kinp
           add(k) = rinp(k)
         enddo
         go to 10
        endif

c
c------ don't look for special parameters if there's no "=" character
        keq = index(line,'=')
        if(keq.le.1) go to 20

c------------------------------------------
c------ look for parameter keywords in front of "=" character

        k = index(line(1:keq-1),'Lunit')
        if(k .ne. 0) then
         lineu = line(keq+1:256)
         call bstrip(lineu,nln)
         kb = index(lineu,' ')
         read(lineu(1:kb),*,err=40) unitl
         unchl = lineu(kb+1:256)
         call bstrip(unchl,nul)
         nul = max(nul,1)
         go to 10
        endif

        k = index(line(1:keq-1),'Munit')
        if(k .ne. 0) then
         lineu = line(keq+1:256)
         call bstrip(lineu,nln)
         kb = index(lineu,' ')
         read(lineu(1:kb),*,err=40) unitm
         unchm = lineu(kb+1:256)
         call bstrip(unchm,num)
         num = max(num,1)
         go to 10
        endif

        k = index(line(1:keq-1),'Tunit')
        if(k .ne. 0) then
         lineu = line(keq+1:256)
         call bstrip(lineu,nln)
         kb = index(lineu,' ')
         read(lineu(1:kb),*,err=40) unitt
         uncht = lineu(kb+1:256)
         call bstrip(uncht,nut)
         nut = max(nut,1)
         go to 10
        endif

c
        k = index(line(1:keq-1),'g')
        if(k .ne. 0) then
         lineu = line(keq+1:256)
         call bstrip(lineu,nln)
         kb = index(lineu,' ')
         read(lineu(1:kb),*,err=40) gee0
         unchgee = lineu(kb+1:256)
         call bstrip(unchgee,nugee)
         nugee = max(nugee,1)
         go to 10
        endif

        k = index(line(1:keq-1),'rho')
        if(k .ne. 0) then
         lineu = line(keq+1:256)
         call bstrip(lineu,nln)
         kb = index(lineu,' ')
         read(lineu(1:kb),*,err=40) rho0
         unchrho = lineu(kb+1:256)
         call bstrip(unchrho,nurho)
         nurho = max(nurho,1)
         go to 10
        endif

c------------------------------
 20     continue

        do i=1, nidim
          rinp(i) = 0.
        enddo

        kinp = nidim
        call getflt(line,rinp,kinp,error)
        if(error) go to 40

        mi = fac(1)*rinp(1) + add(1)
        xi = fac(2)*rinp(2) + add(2)
        yi = fac(3)*rinp(3) + add(3)
        zi = fac(4)*rinp(4) + add(4)
        ixxi = fac(5)*rinp(5) + add(5)
        iyyi = fac(6)*rinp(6) + add(6)
        izzi = fac(7)*rinp(7) + add(7)
        ixyi = fac(8)*rinp(8) + add(8)
        ixzi = fac(9)*rinp(9) + add(9)
        iyzi = fac(10)*rinp(10) + add(10)
        hx = fac(11)*rinp(11) + add(11)
        hy = fac(12)*rinp(12) + add(12)
        hz = fac(13)*rinp(13) + add(13)

        sum_m   = sum_m   + mi
        sum_mx  = sum_mx  + mi*xi
        sum_my  = sum_my  + mi*yi
        sum_mz  = sum_mz  + mi*zi
        sum_mxx = sum_mxx + mi*xi*xi
        sum_myy = sum_myy + mi*yi*yi
        sum_mzz = sum_mzz + mi*zi*zi
        sum_mxy = sum_mxy + mi*xi*yi
        sum_mxz = sum_mxz + mi*xi*zi
        sum_myz = sum_myz + mi*yi*zi

        sum_ixx = sum_ixx + ixxi
        sum_iyy = sum_iyy + iyyi
        sum_izz = sum_izz + izzi
        sum_ixy = sum_ixy + ixyi
        sum_ixz = sum_ixz + ixzi
        sum_iyz = sum_iyz + iyzi

        sum_hx = sum_hx + hx
        sum_hy = sum_hy + hy
        sum_hz = sum_hz + hz

      go to 10

c-----------------------------------
 40   write(*,1500) iline, line(1:80)
 1500 format(' * Bad data line',i5,'  ignored:  ',a)
      go to 10

c=============================================

 50   continue
      close(lu)

      mass = sum_m

      xcg = sum_mx/sum_m
      ycg = sum_my/sum_m
      zcg = sum_mz/sum_m

      ixx = sum_ixx + (sum_myy + sum_mzz) - sum_m*(ycg**2 + zcg**2)
      iyy = sum_iyy + (sum_mzz + sum_mxx) - sum_m*(zcg**2 + xcg**2)
      izz = sum_izz + (sum_mxx + sum_myy) - sum_m*(xcg**2 + ycg**2)
      ixy = sum_ixy +  sum_mxy            - sum_m* xcg*ycg
      ixz = sum_ixz +  sum_mxz            - sum_m* xcg*zcg
      iyz = sum_iyz +  sum_myz            - sum_m* ycg*zcg

      rmass0 = mass * unitm

      riner0(1,1) =  ixx * unitm*unitl**2
      riner0(1,2) = -ixy * unitm*unitl**2
      riner0(1,3) = -ixz * unitm*unitl**2

      riner0(2,1) = -ixy * unitm*unitl**2
      riner0(2,2) =  iyy * unitm*unitl**2
      riner0(2,3) = -iyz * unitm*unitl**2

      riner0(3,1) = -ixz * unitm*unitl**2
      riner0(3,2) = -iyz * unitm*unitl**2
      riner0(3,3) =  izz * unitm*unitl**2

      xyzmass0(1) = xcg * unitl
      xyzmass0(2) = ycg * unitl
      xyzmass0(3) = zcg * unitl

      hrotor0(1) = sum_hx * unitm*unitl**2/unitt
      hrotor0(2) = sum_hy * unitm*unitl**2/unitt
      hrotor0(3) = sum_hz * unitm*unitl**2/unitt

c----------------------------------------
c---- set derived unit values and names
      call unitset

c---- set parameter unit names
      call paruset

      error = .false.
      lmass = .true.
      return

 90   continue
      call bstrip(fname,nfn)
      write(*,9000) fname(1:nfn)
 9000 format(/' Mass file  ',a,'  OPEN error')
      error = .true.
      lmass = .false.
      return
      end ! masget



      subroutine masput(ir1,ir2)
      include 'jvl.inc'
c-------------------------------------------
c     stores default mass, inertias 
c     in run case parameter array
c-------------------------------------------

      do ir = ir1, ir2
        parval(ipmass,ir) = rmass0

        parval(ipixx,ir) = riner0(1,1)
        parval(ipiyy,ir) = riner0(2,2)
        parval(ipizz,ir) = riner0(3,3)
        parval(ipixy,ir) = riner0(1,2)
        parval(ipiyz,ir) = riner0(2,3)
        parval(ipizx,ir) = riner0(3,1)

        parval(ipgee,ir) = gee0
        parval(iprho,ir) = rho0

        parval(ipxcg,ir) = xyzmass0(1)/unitl
        parval(ipycg,ir) = xyzmass0(2)/unitl
        parval(ipzcg,ir) = xyzmass0(3)/unitl

        parval(iphx,ir) = hrotor0(1)
        parval(iphy,ir) = hrotor0(2)
        parval(iphz,ir) = hrotor0(3)
      enddo

      return
      end ! masput



      subroutine massho(lu)
      include 'jvl.inc'

      write(lu,*)
      if(unchm(1:num).ne.'Munit')
     &write(lu,1250) 'Mass        = ', rmass0/unitm, 'Munit'
      write(lu,1250) 'Mass        = ', rmass0, unchm(1:num)

      write(lu,*)
      write(lu,1260) 'Ref. x,y,z  = ',(xyzref0(k)       , k=1,3),'Lunit'
      if(unchl(1:nul).ne.'Lunit')
     &write(lu,1260) 'C.G. x,y,z  = ',(xyzmass0(k)/unitl, k=1,3),'Lunit'
      write(lu,1260) 'C.G. x,y,z  = ',(xyzmass0(k), k=1,3), unchl(1:nul)

      write(lu,*)
      write(lu,1271) 'Ixx -Ixy -Ixz   | ',
     &            riner0(1,1),riner0(1,2),riner0(1,3),'|'
      write(lu,1272) '     Iyy -Iyz = | ',
     &                        riner0(2,2),riner0(2,3),'|', unchi(1:nui)
      write(lu,1273) '          Izz   | ',
     &                                    riner0(3,3),'|'

      write(lu,*)
      write(lu,1260) 'Rotor hx,hy,hz = ',(hrotor0(k),k=1,3),unchh(1:nuh)

 1250 format(1x, a,  g12.4,2x,a)
 1260 format(1x, a, 3g12.4,2x,a)
 1271 format(1x, a,      3g12.4, 2x, a, 2x, a)
 1272 format(1x, a, 12x, 2g12.4, 2x, a, 2x, a)
 1273 format(1x, a, 24x,  g12.4, 2x, a, 2x, a)

      return
      end ! massho


      subroutine appsho(lu,rho)
      include 'jvl.inc'

      write(lu,*) 'Apparent mass, inertia'
      write(lu,*)
      write(lu,1271) 'mxx  mxy  mxz   | ',
     &  amass(1,1)*rho,amass(1,2)*rho,amass(1,3)*rho,'|'
      write(lu,1272) '     myy  myz = | ',
     &                 amass(2,2)*rho,amass(2,3)*rho,'|', unchm(1:num)
      write(lu,1273) '          mzz   | ',
     &                                amass(3,3)*rho,'|'
      write(lu,*)
      write(lu,1271) 'Ixx -Ixy -Ixz   | ',
     &  ainer(1,1)*rho,ainer(1,2)*rho,ainer(1,3)*rho,'|'
      write(lu,1272) '     Iyy -Iyz = | ',
     &                 ainer(2,2)*rho,ainer(2,3)*rho,'|', unchi(1:nui)
      write(lu,1273) '          Izz   | ',
     &                                ainer(3,3)*rho,'|'

 1271 format(1x, a,      3g12.4, 2x, a, 2x, a)
 1272 format(1x, a, 12x, 2g12.4, 2x, a, 2x, a)
 1273 format(1x, a, 24x,  g12.4, 2x, a, 2x, a)

      return
      end ! appsho



      subroutine appget
c----------------------------------------------------------------
c     calculates apparent mass and inertia for unit air density
c----------------------------------------------------------------
      include 'jvl.inc'
      real uc(3), us(3), un(3), rm(3), rxun(3)

      do k = 1, 3
        do l = 1, 3
          amass(k,l) = 0.
          ainer(k,l) = 0.
        enddo
      enddo

      do 100 j = 1, nstrip
        cr = chord(j)
        sr = chord(j)*wstrip(j)

c------ strip normal unit vector un
        un(1) = 0.0
        un(2) = ensy(j)
        un(3) = ensz(j)

c------ spanwise unit vector us (along midchord)
        us(1) = rle2(1,j) - rle1(1,j) + 0.5*(chord2(j) - chord1(j))
        us(2) = rle2(2,j) - rle1(2,j)
        us(3) = rle2(3,j) - rle1(3,j)
        umag = sqrt(us(1)*us(1) + us(2)*us(2) + us(3)*us(3))
        if(umag.gt.0.0) then
         us(1) = us(1)/umag
         us(2) = us(2)/umag
         us(3) = us(3)/umag
        endif

c------ perpendicular-chord unit vector uc
        call cross(us,un,uc)

c------ use the strip 1/2 chord location for moment arm
        rm(1) = rle(1,j) + 0.5*cr
        rm(2) = rle(2,j)
        rm(3) = rle(3,j)

        call cross(rm,un,rxun)

c------ perpendicular chord
        cperp = cr * (us(2)*un(3) - us(3)*un(2))

c------ apparent mass and spanwise-axis apparent inertia of strip
        appm = sr * 0.25*pi*cperp          
        appi = sr * 0.25*pi*cperp**3 / 64.0

c------ add apparent mass contribution
        do k = 1, 3
          do l = 1, 3
            amass(k,l) = amass(k,l) + appm*  un(k)*  un(l) * unitl**3
            ainer(k,l) = ainer(k,l) + appm*rxun(k)*rxun(l) * unitl**5
     &                              + appi*  us(k)*  us(l) * unitl**5
          enddo
        enddo
  100 continue

      return
      end ! appget



      subroutine unitset
      include 'jvl.inc'

c----------------------------------------
c---- set force unit
      unitf = unitm*unitl/unitt**2
      if(unchm(1:1).ne.' ' .and.
     &   unchl(1:1).ne.' ' .and.
     &   uncht(1:1).ne.' '      ) then
       unchf = unchm(1:num)//'-'//unchl(1:nul)//'/'//uncht(1:nut)//'^2'
       call bstrip(unchf,nuf)
       if(unchf(1:nuf).eq.'slug-ft/s^2') unchf = 'lb'
       if(unchf(1:nuf).eq.'kg-m/s^2'   ) unchf = 'N'
       if(unchf(1:nuf).eq.'g-cm/s^2'   ) unchf = 'dyn'
       call bstrip(unchf,nuf)
      endif

c---- set area unit
      units = unitl**2
      if(unchl(1:1).ne.' ') then
       unchs = unchl(1:nul)//'^2'
       call bstrip(unchs,nus)
      endif
c 
c---- set velocity unit
      unitv = unitl/unitt
      if(unchl(1:1).ne.' ' .and.
     &   uncht(1:1).ne.' '      ) then
       unchv = unchl(1:nul)//'/'//uncht(1:nut)
       call bstrip(unchv,nuv)
      endif
c 
c---- set acceleration unit name
      unita = unitl/unitt**2
      if(unchl(1:1).ne.' ' .and.
     &   uncht(1:1).ne.' '      ) then
       uncha = unchl(1:nul)//'/'//uncht(1:nut)//'^2'
       call bstrip(uncha,nua)
      endif
c 
c---- set inertia unit name
      uniti = unitm*unitl**2
      if(unchm(1:1).ne.' ' .and.
     &   unchl(1:1).ne.' '      ) then
       unchi = unchm(1:num)//'-'//unchl(1:nul)//'^2'
       call bstrip(unchi,nui)
      endif

c---- set density unit name
      unitd = unitm/unitl**3
      if(unchm(1:1).ne.' ' .and.
     &   unchl(1:1).ne.' '      ) then
       unchd = unchm(1:num)//'/'//unchl(1:nul)//'^3'
       call bstrip(unchd,nud)
      endif

c---- set angular momentum unit name
      unith = unitm*unitl**2/unitt
      if(unchm(1:1).ne.' ' .and.
     &   unchl(1:1).ne.' ' .and.
     &   uncht(1:1).ne.' '      ) then
       unchh = unchm(1:num)//'-'//unchl(1:nul)//'^2'//'/'//uncht(1:nut)
       call bstrip(unchh,nuh)
      endif

      return
      end ! unitset

