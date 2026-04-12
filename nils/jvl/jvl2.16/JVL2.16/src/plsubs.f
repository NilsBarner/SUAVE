C***********************************************************************
C    Module:  plsubs.f
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

      subroutine plroot(xorg,yorg, 
     &                  irun1, irun2, ircolor,
     &                  jemax,neigen, eval, lproot, lprnum,
     &                  uncht, par,cs, lpgrid,
     &                  xmin,xmax,xdel,xfac,ymin,ymax,ydel,yfac)
      integer ircolor(*)
      integer neigen(*)
      complex eval(jemax,*)
      character*(*) uncht
      logical lproot(jemax,*), lprnum(jemax,*), lpgrid
c-----------------------------------------------------------------
c     Makes root-locus plot
c     If xmin=xmax, then xmin,xmax,xdel are set here.
c     If ymin=ymax, then ymin,ymax,ydel are set here.
c-----------------------------------------------------------------
      data tpi / 6.2831853071795864769 /
      include 'masks.inc'

c---- character sizes
      csl = 1.5*cs
      csu = 1.2*cs
      csa = 1.0*cs
      css = 0.5*cs
      csn = 0.8*cs

c---- root symbol index
      isymb = 1

      call getcolor(icol0)

      if(xmin.eq.xmax) then
       xmin =  1.0e23
       xmax = -1.0e23
       do ir = irun1, irun2
         do keig = 1, neigen(ir)
           if(lproot(keig,ir)) then
            xmin = min(xmin,real(eval(keig,ir)))
            xmax = max(xmax,real(eval(keig,ir)))
           endif
         enddo
       enddo
       call axisadj(xmin,xmax, xtot, xdel, nann)
      endif

      if(ymin.eq.ymax) then
       ymin =  1.0e23
       ymax = -1.0e23
       do ir = irun1, irun2
         do keig = 1, neigen(ir)
           if(lproot(keig,ir)) then
            ymin = min(ymin,imag(eval(keig,ir)))
            ymax = max(ymax,imag(eval(keig,ir)))
           endif
         enddo
       enddo

c----- plot only upper half of plane
       ymin = 0.0
       call axisadj(ymin,ymax, ytot, ydel, nann)
      endif

c---- make sure aspect ratio of plot isn't absurd
      if(xmax-xmin .lt. 2.0*ydel) then
       xmin = xmax - 2.0*ydel
       call axisadj(xmin,xmax, xtot, xdel, nann)
      endif
      if(ymax-ymin .lt. 2.0*xdel) then
c       ymin = 0.5*(ymin+ymax) - xdel
       ymax = ymin + 2.0*xdel
       call axisadj(ymin,ymax, ytot, ydel, nann)
      endif

c---- annotations for hz axis
      ytmin = ymin/tpi
      ytmax = ymax/tpi
      call axisadj2(ytmin,ytmax, yttot, ytdel, nann)
      if(ytdel.gt.0.49*(ytmax-ytmin)) then
       ytdel = 0.5*ytdel
       nann = 2*nann
      endif
      if(ytmin.lt.ymin/tpi) ytmin = ytmin + ytdel
      if(ytmax.gt.ymax/tpi) ytmax = ytmax - ytdel

      xfac = 1.0/(xmax-xmin)
      yfac = par/(ymax-ymin)

      sfac = min( xfac , yfac )
      xfac = sfac
      yfac = sfac

      xlen = (xmax-xmin)*xfac
      ylen = (ymax-ymin)*yfac

      ytfac = yfac*tpi

      call newpen(1)
      if(lpgrid) then
       dxg = 0.5*xdel*xfac
       dyg = 0.5*ydel*yfac

       nxg = int( 2.0*(xmax-xmin)/xdel + 0.5 )
       nyg = int( 2.0*(ymax-ymin)/ydel + 0.5 )
       call plgrid(xorg,yorg, nxg,dxg, nyg,dyg, lmask2 )

      else
       call plot(xorg     ,yorg+ylen,3)
       call plot(xorg+xlen,yorg+ylen,2)
       call plot(xorg+xlen,yorg     ,3)
       call plot(xorg+xlen,yorg+ylen,2)

      endif

      call newpen(2)
      call xaxis(xorg     ,yorg,xlen, xdel*xfac,  xmin, xdel, csa,-2)
      call yaxis(xorg     ,yorg,ylen, ydel*yfac,  ymin, ydel, csa,-2)
      call yaxis(xorg+xlen,yorg-ymin*yfac+ytmin*ytfac,
     &          (ytmax-ytmin)*ytfac,ytdel*ytfac, ytmin,ytdel,-csa,-2)

c---- plot x,y axes as heavy dotted lines
      call newpen(2)
      call newpat(lmask1)
      call plot(xorg     ,yorg-ymin*yfac,3)
      call plot(xorg+xlen,yorg-ymin*yfac,2)
      call plot(xorg-xmin*xfac,yorg     ,3)
      call plot(xorg-xmin*xfac,yorg+ylen,2)
      call newpat(lmask0)

c---- plot axis parameter name annotations
      call newpen(3)
      xl = xorg + xlen - 0.5*csl - 1.5*xdel*xfac
      yl = yorg        - 2.3*csl
      call plmath(xl  ,yl, csl,'s'    ,0.0, 1)
      xl = xorg + xlen - 1.5*csu - 0.5*xdel*xfac
      call plchar(xl  ,yl, csu,'1/'   ,0.0, 2)
      call plchar(999.,yl, csu, uncht ,0.0,-1)

      xl = xorg        - 3.0*csl
      yl = yorg + ylen - 0.5*csl - 0.5*ydel*yfac
      call plmath(xl  ,yl, csl,'w',0.0, 1)
      xl = xorg        - 4.2*csu
      yl = yorg + ylen - 0.5*csl - 1.5*ydel*yfac
      call plchar(xl  ,yl, csu,'1/',0.0,2)
      call plchar(999.,yl, csu, uncht,0.0,-1)

      xl = xorg + xlen      + 2.0*csl
      yl = yorg - ymin*yfac - 0.5*csl + (ytmax - 0.5*ytdel)*ytfac
      call plmath(xl  ,yl, csl,'w  p',0.0, 4)
      call plchar(xl  ,yl, csl,' /2 ',0.0, 4)
      xl = xorg + xlen      + 1.2*csu
      yl = yorg - ymin*yfac - 0.5*csl + (ytmax - 1.5*ytdel)*ytfac
      call plchar(xl  ,yl, csu,'cycles/',0.0,7)
      call plchar(999.,yl, csu, uncht,0.0,-1)

c---- plot root symbols
      xoff = xmin - xorg/xfac
      yoff = ymin - yorg/yfac
      do ir = irun1, irun2
        do keig = 1, neigen(ir)
          if(lproot(keig,ir)) then
           xev = real(eval(keig,ir))
           yev = imag(eval(keig,ir))
           if(yev .ge. -1.0e-4) then
            call newpen(5)
            call newcolor(ircolor(ir))
            xplt = (xev - xoff)*xfac
            yplt = (yev - yoff)*yfac
            call xysymb(1,xplt,yplt,0.0,1.0,0.0,1.0,css,isymb)
            if(lprnum(keig,ir)) then
             call newpen(2)
             xnum = xplt + 0.65*css
             ynum = yplt + 0.65*css
             call plnumb(xnum,ynum,csn,float(ir),0.0,-1)
            endif
           endif
          endif
        enddo
      enddo
      call newcolor(icol0)

      call plflush

      return
      end ! plroot



      subroutine pltpar(xplt,yplt, irun1, irun2,
     &                  parval,lppar,
     &                  ircolor, csiz, delx, dely )
c-----------------------------------------------
c     plots operating conditions in table form
c-----------------------------------------------
      include 'jindex.inc'

      real parval(iptot,*)
      logical lppar(iptot)
      integer ircolor(*)

c
      integer npdig(iptot), npwid(iptot)
      real xpar(iptot), facpar(iptot)

      data dtr / 0.0174532925 /

c
      cs  =     csiz
      csl = 1.3*csiz

c---- set number of significant digits for each variable
      do ip = 1, iptot
        npdig(ip) = 4
      enddo

c---- set total number of digits for each variable (including decimal point)
      do ip = 1, iptot
        npwid(ip) = 3 + npdig(ip)
      enddo

c---- set scale factors of displayed value
      do ip = 1, iptot
        facpar(ip) = 1.0
      enddo
cc      facpar(ipalfa) = 1.0/dtr
cc      facpar(ipbeta) = 1.0/dtr
cc      facpar(ipthe ) = 1.0/dtr
cc      facpar(ipphi ) = 1.0/dtr

c---- starting x-position
      xpos = xplt + 2.5*cs + delx

c---- set starting x position for each variable
      do ip = 1, iptot
        if(lppar(ip)) then
         xpar(ip) = xpos
         xpos = xpos + cs*npwid(ip) + delx
        else
         xpar(ip) = 0.
        endif
      enddo

c
      call getcolor(icol0)

      y = yplt
      do ir = irun2, irun1, -1
        y = y + dely

        call newcolor(ircolor(ir))

        call newpen(2)
        x = xplt
        call plnumb(x,y,cs,float(ir),0.0,-1)
        call plchar(999.,y,cs,':',0.0,1)

c------ plot numerical parameter values
        call newpen(2)
        do ip = 1, iptot
          if(lppar(ip)) then
           x = xpar(ip)
           rnum = parval(ip,ir)*facpar(ip)
           call plnums(x,y,cs,rnum,0.0,npdig(ip))
          endif
        enddo
      enddo

      call newcolor(icol0)
      call newpen(3)

      y = y + dely + 0.3*cs

      do ip = 1, iptot
        if(lppar(ip)) then
         x = xpar(ip) + 0.5*cs*float(npwid(ip)-3) - 0.5*csl
         call ppname(x,y,csl, ip)
        endif
      enddo

c---- save upper y limit for passing back
      yplt = y + csl

      return
      end



      subroutine ppname(x,y,csiz,ip)
      include 'jindex.inc'
c------------------------------------------------------------------
c     plots symbol name of parameter specified by ip index
c
c      x,y   position of lower-left corner of symbol
c      csiz  size of symbol
c      ip    index of parameter
c------------------------------------------------------------------
      csp =      csiz
      csc = 0.70*csiz

      call getlastxy(xplt,yplt)
      if(x .ne. 999.0) xplt = x
      if(y .ne. 999.0) yplt = y

c------ plot scalar parameter name

      if    (ip.eq.ipalfa) then
       call plmath(xplt,yplt,csp,'a"' ,0.0,2)

      elseif(ip.eq.ipbeta) then
       call plmath(xplt,yplt,csp,'b"' ,0.0,2)

      elseif(ip.eq.iprotx) then
       call plchar(xplt,yplt,csp,'p'  ,0.0,1)

      elseif(ip.eq.iproty) then
       call plchar(xplt,yplt,csp,'q'  ,0.0,1)

      elseif(ip.eq.iprotz) then
       call plchar(xplt,yplt,csp,'r'  ,0.0,1)

      elseif(ip.eq.ipcl  ) then
       call plchar(xplt        ,yplt        ,csp,'C' ,0.0,1)
       call plchar(xplt+0.9*csp,yplt-0.4*csp,csc, 'L',0.0,1)

      elseif(ip.eq.ipcd0 ) then
       call plchar(xplt        ,yplt        ,csp,'C'  ,0.0,1)
       call plchar(xplt+0.9*csp,yplt-0.4*csp,csc, 'Do',0.0,2)

      elseif(ip.eq.ipthe ) then
ccc       call plmath(xplt,yplt,csp,'Q"' ,0.0,2)
       call plchar(xplt,yplt,csp,'elev' ,0.0,4)

      elseif(ip.eq.ipphi ) then
ccc       call plmath(xplt,yplt,csp,'F"' ,0.0,2)
       call plchar(xplt,yplt,csp,'bank' ,0.0,4)

      elseif(ip.eq.ipvee ) then
       call plchar(xplt,yplt,csp,'V'  ,0.0,1)

      elseif(ip.eq.iprho ) then
       call plmath(xplt,yplt,csp,'r'  ,0.0,1)

      elseif(ip.eq.ipgee ) then
       call plchar(xplt,yplt,csp,'g'  ,0.0,1)

      elseif(ip.eq.iprad ) then
       call plchar(xplt        ,yplt        ,csp,'R'   ,0.0,1)
       call plchar(xplt+0.9*csp,yplt-0.4*csp,csc,'turn',0.0,4)


      elseif(ip.eq.ipfac ) then
       call plchar(xplt,yplt,csp,'N'  ,0.0,1)

      elseif(ip.eq.ipxcg ) then
       call plchar(xplt        ,yplt        ,csp,'X'  ,0.0,1)
       call plchar(xplt+0.9*csp,yplt-0.4*csp,csc, 'cg',0.0,2)

      elseif(ip.eq.ipycg ) then
       call plchar(xplt        ,yplt        ,csp,'Y'  ,0.0,1)
       call plchar(xplt+0.9*csp,yplt-0.4*csp,csc, 'cg',0.0,2)

      elseif(ip.eq.ipzcg ) then
       call plchar(xplt        ,yplt        ,csp,'Z'  ,0.0,1)
       call plchar(xplt+0.9*csp,yplt-0.4*csp,csc, 'cg',0.0,2)

      elseif(ip.eq.ipmass) then
       call plchar(xplt        ,yplt        ,csp,'mass',0.0,4)

      elseif(ip.eq.ipixx ) then
       call plchar(xplt        ,yplt        ,csp,'I'  ,0.0,1)
       call plchar(xplt+0.9*csp,yplt-0.4*csp,csc, 'xx',0.0,2)

      elseif(ip.eq.ipiyy ) then
       call plchar(xplt        ,yplt        ,csp,'I'  ,0.0,1)
       call plchar(xplt+0.9*csp,yplt-0.4*csp,csc, 'yy',0.0,2)

      elseif(ip.eq.ipizz ) then
       call plchar(xplt        ,yplt        ,csp,'I'  ,0.0,1)
       call plchar(xplt+0.9*csp,yplt-0.4*csp,csc, 'zz',0.0,2)

      elseif(ip.eq.ipixy ) then
       call plchar(xplt        ,yplt        ,csp,'I'  ,0.0,1)
       call plchar(xplt+0.9*csp,yplt-0.4*csp,csc, 'xy',0.0,2)

      elseif(ip.eq.ipiyz ) then
       call plchar(xplt        ,yplt        ,csp,'I'  ,0.0,1)
       call plchar(xplt+0.9*csp,yplt-0.4*csp,csc, 'yz',0.0,2)

      elseif(ip.eq.ipizx ) then
       call plchar(xplt        ,yplt        ,csp,'I'  ,0.0,1)
       call plchar(xplt+0.9*csp,yplt-0.4*csp,csc, 'zx',0.0,2)

      elseif(ip.eq.ipcla ) then
       call plmath(xplt        ,yplt        ,csp,'D'  ,0.0,1)
       call plchar(xplt+1.0*csp,yplt        ,csp,'C'  ,0.0,1)
       call plchar(xplt+1.9*csp,yplt-0.4*csp,csc, 'L' ,0.0,1)
       call plmath(xplt+2.7*csp,yplt-0.6*csp,csc,  'a',0.0,1)

      elseif(ip.eq.ipclu ) then
       call plmath(xplt        ,yplt        ,csp,'D'  ,0.0,1)
       call plchar(xplt+1.0*csp,yplt        ,csp,'C'  ,0.0,1)
       call plchar(xplt+1.9*csp,yplt-0.4*csp,csc, 'L' ,0.0,1)
       call plchar(xplt+2.7*csp,yplt-0.6*csp,csc,  'u',0.0,1)

      elseif(ip.eq.ipcma ) then
       call plmath(xplt        ,yplt        ,csp,'D'  ,0.0,1)
       call plchar(xplt+1.0*csp,yplt        ,csp,'C'  ,0.0,1)
       call plchar(xplt+1.9*csp,yplt-0.4*csp,csc, 'M' ,0.0,1)
       call plmath(xplt+2.7*csp,yplt-0.6*csp,csc,  'a',0.0,1)

      elseif(ip.eq.ipcmu ) then
       call plmath(xplt        ,yplt        ,csp,'D'  ,0.0,1)
       call plchar(xplt+1.0*csp,yplt        ,csp,'C'  ,0.0,1)
       call plchar(xplt+1.9*csp,yplt-0.4*csp,csc, 'M' ,0.0,1)
       call plchar(xplt+2.7*csp,yplt-0.6*csp,csc,  'u',0.0,1)

      endif

c---- make sure subscript didn't screw up vertical pen location
      call plchar(999.0,yplt,0.001*csp,' ',0.0,1)

      return
      end ! ppname


      subroutine plnums(x,y,cs,rnum,ang,ndig)
c----------------------------------------------------------
c     plots real number with ndig significant digits
c----------------------------------------------------------

      if(rnum.eq.0.0) then
       xplt = x + cs*float(ndig/2-1)
       call plnumb(xplt,y,cs,rnum,ang,1)
       return
      endif

      arnum = abs(rnum)
      rlog = log10(arnum)

c---- set number of digits to left and right of decimal point
      nl = int( rlog + 101.0 ) - 100
      nr = ndig - nl

      if(nl.gt.ndig .or. nr.gt.ndig) then
c----- plot with scientific notation
       nexp = 1 - nl
       fexp = 10.0**nexp
       nd = max( ndig - 1 , 1 )
       call plnumb(x   ,y,cs,rnum*fexp,ang,nd)
       call plchar(999.,y,cs,'e',ang,1)
       call plnumb(999.,y,cs,float(-nexp),ang,-1)
      else
       nd = max( nr , 0 )
       call plnumb(x,y,cs,rnum,ang,nd)
      endif

      return
      end
