C***********************************************************************
C    Module:  getvm.f
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

      subroutine getvm(fnamvb)
c---------------------------------------------------------------------
c     Prints strip shear and bending quantities, ie. v, bm
c            integrates spanload to get shear and bending moment
c     Note:  only works for single surface at at time (ie. v,bm on each panel)
c---------------------------------------------------------------------
      include 'jvl.inc'
      logical lfile
      real fs(3), ms(3)
      real v(nsmax), bm(nsmax), ystrp(nsmax)
      character*72 fnamvb

    1 format(a)
    2 format(a,i3,a)

c---- open output file for loading info, or default to previously used file
      write(*,*) '  '
      call bstrip(fnamvb,nfn)
      write(*,3) fnamvb(1:nfn)
   3  format('Enter shear/bending output file (return for screen): ', a)
      read (*,1) fnamvb
      call bstrip(fnamvb,nfn)

      lfile = fnamvb.ne.' '

      if(lfile) then
       lu = 19
       open(lu,file=fnamvb,status='unknown')
      else
       lu = 6
      endif

      write(lu,10) title(1:60),mach,alfa/dtr,
     &             ltot/(0.5*sref),beta/dtr,sref,bref
 10   format(/' Shear/q and Bending Moment/q vs y'
     &       /'  Configuration: ',a
     &       /'  Mach  = ',f8.3,
     &       /'  alpha = ',f8.3,'    CLtot = ',f8.3,
     &       /'  beta  = ',f8.3,
     &      //'  Sref  = ',f11.5
     &       /'  Bref  = ',f11.5)

c---- process the surfaces one by one, calculating shear and bending on each, 
c      with moments refered to y=0 (centerline)

      do n = 1, nsurf
        ymin =  1.0e10
        ymax = -1.0e10
        do j = jfrst(n), jlast(n)
          iv = ifrsts(j)
          ymin = min( ymin , rv1(2,iv) , rv2(2,iv) )
          ymax = max( ymax , rv1(2,iv) , rv2(2,iv) )
        enddo

c------ integrate spanload from first strip to last strip defined for 
c        this surface to get shear and bending moment
        ccllst = 0.0
        bmlst  = 0.0
        wlst   = 0.0
        vlst   = 0.0

        jf = jfrst(n)
        jl = jlast(n)
        jinc = 1

c------ integrate from first to last strip in surface
        do j = jl, jf, -jinc
          call fstrip(j,
     &              fs,ms,
     &              caxial,cnorml, 
     &              cl_strp,cd_strp,
     &              clj_strp,cdj_strp,
     &              cltstrp,clastrp,
     &              cmc4_strp,cmle_strp,
     &              ccl_strp )

          jj = jinc*(j - jf + jinc)

          dy = 0.5*(wstrip(j)+wlst)
          ystrp(jj) = rle(2,j)
          v(jj)     = vlst  + 0.5*(ccl_strp+ccllst) * dy 
          bm(jj)    = bmlst + 0.5*(v(jj)+vlst)    * dy

          vlst   = v(jj) 
          bmlst  = bm(jj)
          ccllst = ccl_strp
          wlst   = wstrip(j)
        enddo

c------ inboard edge y,vz,mx
        vroot  = vlst  +      ccllst      * 0.5*dy 
        bmroot = bmlst + 0.5*(vroot+vlst) * 0.5*dy
        vtip   = 0.0
        bmtip  = 0.0
        if(idups(n).ge.0) then
         yroot = rle1(2,jf)
         ytip  = rle2(2,jl)
        else
         yroot = rle2(2,jf)
         ytip  = rle1(2,jl)
        endif

        dir = 1.0
        if(ymin+ymax.lt.0.0) dir = -1.0

        write(lu,*) '  '
        write(lu,2) ' Surface: ',n,'  ',stitle(n)
        write(lu,*) '    2Ymin/Bref = ',2.0*ymin/bref
        write(lu,*) '    2Ymax/Bref = ',2.0*ymax/bref
        write(lu,*) '  2y/Bref      Vz/(q*Sref)      Mx/(q*Bref*Sref)'

        write(lu,4)   2.*yroot/bref,vroot/sref,dir*bmroot/sref/bref
        do j = 1, nj(n)
          write(lu,4) 2.*ystrp(j)/bref,v(j)/sref,dir*bm(j)/sref/bref
        enddo
        write(lu,4)   2.*ytip/bref,vtip/sref,dir*bmtip/sref/bref
   4    format(f10.4,g14.6,3x,g14.6)

c       n = 0
c       call grphin(n,ystrp,v)
c       call grphin(nj(n),ystrp,v)
c       call grphin(nj(n),ystrp,bm)
c       call grphpl
      enddo

      if(lfile) close(lu)

      return
      end ! getvm
