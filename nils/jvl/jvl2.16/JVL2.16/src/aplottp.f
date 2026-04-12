C***********************************************************************
C    Module:  aplottp.f
C 
C    Copyright (C) 2002 Mark Drela, Harold Youngren
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

      subroutine plottp
c-------------------------------------------------
c     plots circulation, cl, downwash versus y
c-------------------------------------------------
      include 'avl.inc'
      include 'avlplt.inc'
c
      character opt*2, line*80, yaxtype*1, satype*16, plttype*16
      real xtmp(2), ytmp(2)
      real sle(nsmax), xplt, yplt
      logical lcolhc
      character*24 col1, col2, col3, col4
c
      real rinp(3)
      logical error
c
      data  yaxtype / 'y' /
      data  plttype / 'yzcnc' /
c
      include 'masks.inc'
c
c--- plfac is scale factor adjustment for this whole plot 
c--- wfac  is scale factor for w axis (as fraction of c axis)
c--- vtitlefrac is fraction of vertical axis to use for big label block
      data plfac, wfac, vtitlefrac / 0.85, 0.5, 0.15 /
c
      smod(ss) = ssf*(ss - soff)
      xmod(xx) = xsf*(xx - xoff)
      ymod(yy) = ysf*(yy - yoff)
c
      if(.not.lsol) then
        write(*,*) '*** no flow solution...'
        return
      endif
c
c---- initialize plot stuff
      sizeold = size
      csold = cs
      size  = 9.0
      cs    = 0.012
c---- number of control-variable lines in excess of 6
ccc      nconlin = (ncontrol+1)/2
      nconlin = max( ncontrol-6 , 0 )
c
c---- line y-spacing factor
      ysp = 2.0
c
c--- define colors used for plot data
      col1 = 'green'   ! cl c/cref
      col2 = 'blue'    ! alpha_i
      col3 = 'red'     ! cl_perp
ccc   col4 = 'yellow'  ! cl
      col4 = 'orange'  ! cl
c
c
      call getsa(lnasa_sa,satype,dir)
c
      ca = cos(alfa)
      sa = sin(alfa)
      crsax = crtot*ca + cntot*sa
      cmsax = cmtot              
      cnsax = cntot*ca - crtot*sa  
c
c---- find the min/max of y,z,cl,cnc,downwash
 50   ymin =  rv1(2,1)
      ymax =  rv1(2,1)
      zmin =  rv1(3,1)
      zmax =  rv1(3,1)
      cymin = ymin
      cymax = ymax
      czmin = zmin
      czmax = zmax
      fmin =  cnc(1)
      fmax =  cnc(1)
      cmin = 0.0
      cmax = 0.0
      wmax = 0.0
      wmin = 0.0
      do i=1, nstrip
        iv = ijfrst(i)
        ymin = min( ymin, rv1(2,iv), rv2(2,iv), 0.0 )
        ymax = max( ymax, rv1(2,iv), rv2(2,iv), 0.0 )
        zmin = min( zmin, rv1(3,iv), rv2(3,iv), 0.0 )
        zmax = max( zmax, rv1(3,iv), rv2(3,iv), 0.0 )
        fmin = min( fmin, cnc(i), 0.0 )
        fmax = max( fmax, cnc(i), 0.0 )
c---- use either min-max or l8 norms for axis length determination
c     l8 prevents large local spikes from giving inappropriate axes
c--- l8 norm
        wmin = wmin + min( -dwwake(i)  , 0.0 ) ** 8
        wmax = wmax + max( -dwwake(i)  , 0.0 ) ** 8
c--- regular min/max
        cmin = min( cmin, cltstrp(i)   , 0.0 )
        cmax = max( cmax, cltstrp(i)   , 0.0 )
c        wmin = min( wmin, -dwwake(i)  , 0.0 )
c        wmax = max( wmax, -dwwake(i)  , 0.0 )
      end do
c
c--- process l8 sums to norms
      if(wmin .ne. 0.0) wmin = -(abs(wmin)/nstrip) ** 0.125
      if(wmax .ne. 0.0) wmax =  (abs(wmax)/nstrip) ** 0.125
c
c---- normalized spanload cnc/cref and cl will have the same axis, so...
      cmin = min(fmin/cref,cmin)
      cmax = max(fmax/cref,cmax)
c
      if(cmax-cmin .lt. 1.0e-5) then
       cmin = 0.0
       cmax = 0.1
      endif
c
      if(wmax-wmin .lt. 1.0e-5) then
       wmin = 0.0
       wmax = 0.1
      endif
c
c---- determine "nice" upper,lower bounds on y, z, cnc/cref, cl, 
c-    and corresponding "nice" number of annotations
      call axisadj(ymin,ymax,yspan,dely,nyann)
      call axisadj(zmin,zmax,zspan,delz,nzann)
      call axisadj(wmin,wmax,wspan,delw,nwann)
      call axisadj(cmin,cmax,cspan,delc,ncann)
c---- ratio of w and cl scale factors
      woc = wfac*(cmax-cmin)/(wmax-wmin)
c
c
c
c---------------------------------------------------
c---- plot yz trefftz plane geometry and loading
 100  if(plttype.ne.'trace') go to 150 
c
c---- start with geometry sizing scale
      ysf  = 1.0    / (ymax-ymin)
      zsf  = plotar / (zmax-zmin)
      xoff = ymin
      yoff = zmin
      sf = min(ysf,zsf)
      xsf = sf
      ysf = sf
c      write(*,*) 'ymin,ymax ',ymin,ymax,xsf
c      write(*,*) 'zmin,zmax ',zmin,zmax,ysf
c      vecscl = 0.3*plotar/sf
c      vecscl = 0.5*plotar/sf
      vecscl = plotar/sf
c      fmag = fmax-fmin
c      if(fmag.ne.0) vecscl = vecscl/fmag
      write(*,*) 'fmin,fmax ',fmin,fmax,fmag
c      write(*,*) 'vecscl ',vecscl
c
      ymin =  rv1(2,1)
      ymax =  rv1(2,1)
      zmin =  rv1(3,1)
      zmax =  rv1(3,1)
      cymin = ymin
      cymax = ymax
      czmin = zmin
      czmax = zmax
      cncyzav = 0.0
      cncyzmx = 0.0
      do n=1, nsurf
        j1 = jfrst(n)
        do j = 1, nj(n)
          jj = j1+j-1
          y0 = rle(2,jj)
          z0 = rle(3,jj)
          y1 = y0 + vecscl*cnc(jj)*ensy(jj) 
          z1 = z0 + vecscl*cnc(jj)*ensz(jj) 
          ymin  = min( ymin, y0 )
          ymax  = max( ymax, y0 )
          zmin  = min( zmin, z0 )
          zmax  = max( zmax, z0 )
          cymin = min( cymin, y1 )
          cymax = max( cymax, y1 )
          czmin = min( czmin, z1 )
          czmax = max( czmax, z1 )
          cncyzav = cncyzav + sqrt((y1-y0)**2 + (z1-z0)**2)
          cncyzmx = max(cncyzmx,sqrt((y1-y0)**2 + (z1-z0)**2))
        end do
      end do
      cncyzav = cncyzav/nstrip
      write(*,*) 'ymin,ymax ',ymin,ymax
      write(*,*) 'zmin,zmax ',zmin,zmax
      write(*,*) 'cymin,cymax ',cymin,cymax
      write(*,*) 'czmin,czmax ',czmin,czmax
      write(*,*) 'cncyzav,cncyzmx ',cncyzav,cncyzmx
c
c---- rescale cnc vectors and geometry to fit
      dcymin = cymin - ymin
      dcymax = cymax - ymax
      dczmin = czmin - zmin
      dczmax = czmax - zmax
cc      call axisadj(cymin,cymax,cyspan,delcy,ncyann)
cc      call axisadj(czmin,czmax,czspan,delcz,nczann)
      yspan = ymax - ymin
      zspan = zmax - zmin
      cyspan = cymax - cymin
      czspan = czmax - czmin
      write(*,*) 'yspan ',yspan
      write(*,*) 'zspan ',zspan
      write(*,*) 'cymin,cymax,cyspan ',cymin,cymax,cyspan
      write(*,*) 'czmin,czmax,czspan ',czmin,czmax,czspan

c---- scale cnc vector envelope to 30% beyond geometry y,z
cc      cscl = 1.3/(max(cyspan/yspan,czspan/zspan))
cc      cscl = 0.3*max(yspan,zspan)/cncyzav 
cc      cscl = 0.5*max(yspan,zspan)/cncyzmx 
      cscl = 0.3*max(yspan,zspan)/cncyzmx 
c      cscl = 0.1*max(yspan,zspan)/cncyzmx 
      write(*,*) 'cscl ',cscl
c      ymin = ymin - max(0.0,cscl*(cyspan-yspan))
c      ymax = ymax + max(0.0,cscl*(cyspan-yspan))
c      zmin = zmin - max(0.0,cscl*(czspan-zspan))
c      zmax = zmax + max(0.0,cscl*(czspan-zspan))
c
      ymin = ymin - max(0.0,cscl*(cymin-ymin))
      ymax = ymax + max(0.0,cscl*(cymax-ymax))
      zmin = zmin - max(0.0,cscl*(czmin-zmin))
      zmax = zmax + max(0.0,cscl*(czmax-zmax))
      write(*,*) 'ymin,ymax ',ymin,ymax
      write(*,*) 'zmin,zmax ',zmin,zmax
c---- reset "nice" y,z plot limits
      call axisadj(ymin,ymax,yspan,dely,nyann)
      call axisadj(zmin,zmax,zspan,delz,nzann)
      write(*,*) 'adj ymin,ymax ',ymin,ymax
      write(*,*) 'adj zmin,zmax ',zmin,zmax
c
      xsf  = 1.0    / (ymax-ymin)
      ysf  = plotar / (zmax-zmin)
      xoff = ymin
      yoff = zmin
      sf = min(xsf,ysf)
      xsf = sf
      ysf = sf
      vecscl = cscl*vecscl
c      write(*,*) 'cscl, vecscl ',cscl,vecscl
c
c---- start the y,z,cnc plot
      call pltini(idev)
      call newfactor(plfac*size)
      call getcolor(icol0)
c
c----offset the x-axis, leaving room on left for labels
      call plot(12.0*cs,8.0*cs,-3)
c
c---put up the xaxis (suppressing the end annotations)
      call newpen(3)
      xlab = xmod(ymax - 1.5*dely)
      ylab = ymod(zmax - 1.5*delz)
      iflg = 1
      call xaxis2(0.0,0.0,xsf*(ymax-ymin),xsf*dely,
     &            ymin,dely,iflg,0.9*cs,-2)
      call plchar(xlab-0.7*cs,-3.5*cs,1.2*cs,'y',0.,1)
c
      call yaxis(0.0,0.0,ysf*(zmax-zmin),ysf*delz,
     &           zmin,delz, 0.9*cs,-2)
      call plchar(-6.5*cs,ylab-0.7*cs,1.2*cs,'z',0.,1)
c
c---put up a reference grid on y,z 
      nx = ifix(0.1 + (ymax-ymin)/dely)
      ny = ifix(0.1 + (zmax-zmin)/delz)
      call plgrid(0.0,0.0,nx,xsf*dely,ny,ysf*delz,lmask2)

c---plot the y,z trace
      call newpen(4)
      call newcolorname(col1)
      call getcolor(icol1)
      call newcolorname(col4)
      call getcolor(icol4)
c
      do n = 1, nsurf
        j1 = jfrst(n)
cc        write(*,*) 'n,j,nj ',n,j,nj(n)
c
c---- plot y,z geometry of far-field trace 
        call newcolor(icol4)
        do j = 1, nj(n)
          iv = ijfrst(j1+j-1)
          if(imags(n).eq.1) then
           if(j.eq.1) then
            xplt = xmod(rv1(2,iv))
            yplt = ymod(rv1(3,iv))
            call plot(xplt,yplt,3)
           else
            xplt = xmod(rv2(2,iv))
            yplt = ymod(rv2(3,iv))
            call plot(xplt,yplt,2)
           endif
          else
           if(j.eq.1) then
            xplt = xmod(rv2(2,iv))
            yplt = ymod(rv2(3,iv))
            call plot(xplt,yplt,3)
           else
            xplt = xmod(rv1(2,iv))
            yplt = ymod(rv1(3,iv))
            call plot(xplt,yplt,2)
           endif
          endif
        end do
c---- plot cnc loading as arrows
        call newcolor(icol1)
        do j = 1, nj(n)
          jj = j1+j-1
          x0 = rle(2,jj)
          y0 = rle(3,jj)
          x1 = x0 + vecscl*cnc(jj)*ensy(jj) 
          y1 = y0 + vecscl*cnc(jj)*ensz(jj) 
          xplt = xmod(x0)
          yplt = ymod(y0)
          call plot(xplt,yplt,3)
          xplt = xmod(x1)
          yplt = ymod(y1)
          call plot(xplt,yplt,2)
          dx = x1 - x0
          dy = y1 - y0
          x2 = x0 + 0.8*dx + 0.02*dy
          y2 = y0 + 0.8*dy - 0.02*dx
          x3 = x0 + 0.8*dx - 0.02*dy
          y3 = y0 + 0.8*dy + 0.02*dx
          call plot(xmod(x2),ymod(y2),2)
          call plot(xmod(x3),ymod(y3),2)
          call plot(xmod(x1),ymod(y1),2)
        end do
c
        if(llabsurf) then
          jlab = j1 + 0.50*nj(n)
          xlab = xmod(rle(2,jlab))         
          ylab = ymod(rle(3,jlab))         
          call plnumb(xlab+0.5*cs,ylab+0.5*cs,0.8*cs,float(n),0.0,-1)
        endif
      enddo
c
      call plflush
      call newcolor(icol0)
c
c--- forces block at top of y,z,cnc plot
      xlab = 0.0
      ylab = ysf*(zmax-zmin) + ysp*cs*(8.0+float(nconlin))
      go to 60
c
c---------------------------------------------------
c---- plot cnc/cref, cl, dw vs y or z
 150  continue
      if(yaxtype.eq.'y') then
        smin = ymin
        smax = ymax
        dels  = dely
        do j = 1, nstrip
          sle(j) = rle(2,j)
        end do
       else
        smin = zmin
        smax = zmax
        dels  = delz
        do j = 1, nstrip
          sle(j) = rle(3,j)
        end do
      endif
c
 155  continue
      ssf  = 1.0 / (smax-smin)
      soff = smin
      if(ldwashplt ) then
       csf = plotar / (max(cmax,woc*wmax)-min(cmin,woc*wmin))
c---- modify csf to allow portion of vertical space for a big label block
       csf = (1.0 - vtitlefrac)*csf
       wsf = woc * csf
      else
       csf = plotar / (cmax-cmin)
c---- modify csf to allow portion of vertical space for a big label block
       csf = (1.0 - vtitlefrac)*csf
       wsf = 0.0
      endif
c
      call pltini(idev)
      call newfactor(plfac*size)
      call getcolor(icol0)
c
c----offset the x-axis, leaving room on left for labels
      fmin = min(wsf*wmin,csf*cmin)
      call plot(14.0*cs,6.0*cs-fmin,-3)
c
c---put up the xaxis (suppressing the end annotations)
      call newpen(3)
      xlab = smod(smax - 1.5*dels)
      if(ldwashplt ) then
       iflg = 3
      else
       iflg = 1
      endif
      call xaxis2(0.0,0.0,ssf*(smax-smin),ssf*dels,
     &            smin,dels,iflg,0.9*cs,-2)
      call plchar(xlab-0.7*cs,-3.5*cs,1.2*cs,yaxtype,0.,1)
c
c--- cl_perp, cl, cl*c/cref axis
      call newpen(2)
      xlab = smod(smin)
c
      if(lclperplt) then
       ylab = csf*(cmax - 0.5*delc)
       call newcolorname(col3)
       call plchar(xlab-6.5*cs,ylab-0.5*cs,1.2*cs,'c'  ,0., 1)
       call plmath(xlab-5.5*cs,ylab-0.9*cs,0.8*cs, 'v' ,0., 1)
       call plchar(xlab-4.2*cs,ylab-0.4*cs,0.7*cs,'t'  ,180.0, 1)
      endif
c
      ylab = csf*(cmax - 1.5*delc)
      call newcolorname(col4)
      call plchar(xlab-6.5*cs,ylab-0.5*cs,1.2*cs,'c'  ,0., 1)
      call plmath(xlab-5.5*cs,ylab-0.9*cs,0.8*cs, 'v' ,0., 1)
c
      ylab = csf*(cmax - 2.5*delc)
      call newcolorname(col1)
      call plchar(xlab-8.5*cs,ylab-0.5*cs,1.2*cs,'c c/c',0.,5)
      call plchar(999.,999.,0.65*cs,'ref',0.,3)
      call plmath(xlab-7.5*cs,ylab-0.9*cs,0.8*cs, 'v'    ,0.,1)
      call newcolor(icol0)
      call yaxis(0.0,csf*cmin,csf*(cmax-cmin),csf*delc,
     &           cmin,delc, 0.9*cs,-2)
c
c---put up a reference grid on the y/z and cl_perp, cl, cl*c/cref axes
      ny = ifix(0.1 + (smax-smin)/dels)
      nc = ifix(0.1 + (cmax-cmin)/delc)
      call plgrid(0.0,csf*cmin,ny,ssf*dels,nc,csf*delc,lmask2)
c
c---downwash axis
      if(ldwashplt ) then
        call newpen(3)
        call newcolorname(col2)
        xlab = smod(smax)
        ylab = wsf*(wmin + 0.5*delw)
        call plmath(xlab+6.5*cs,ylab-0.4*cs,1.2*cs,'a',0., 1)
        call plchar(xlab+7.6*cs,ylab-0.8*cs,0.9*cs,'i',0., 1)
        call newpen(4)
        call yaxis(smod(smax),wsf*wmin,wsf*(wmax-wmin),wsf*delw,
     &             wmin,delw,-0.9*cs,-2)
      endif
c
c---plot the cl,clc/cref aand downwash curves
      call newpen(4)
      do n = 1, nsurf
        j  = jfrst(n)
        jlabclc = j + 0.33*nj(n)
        jlabdw  = j + 0.50*nj(n)
        jlabcl  = j + 0.66*nj(n)
c
        call newcolorname(col4)
        call xyline(nj(n),sle(j),clastrp(j),smin,ssf,0.0,csf,2)
        if(llabsurf) then
          xlab = smod(sle(jlabcl))         
          ylab = csf*clastrp(jlabcl)         
          call plnumb(xlab+0.5*cs,ylab+0.5*cs,0.8*cs,float(n),0.0,-1)
        endif
c
        if(lclperplt) then
         call newcolorname(col3)
         call xyline(nj(n),sle(j),cltstrp(j),smin,ssf,0.0,csf,3)
         if(llabsurf) then
           xlab = smod(sle(jlabcl))         
           ylab = csf*cltstrp(jlabcl)         
           call plnumb(xlab+0.5*cs,ylab+0.5*cs,0.8*cs,float(n),0.0,-1)
         endif
        endif
c
        call newcolorname(col1)
        call xyline(nj(n),sle(j),cnc(j),smin,ssf,0.0,csf/cref,1)
        if(llabsurf) then
          xlab = smod(sle(jlabclc))         
          ylab = csf/cref*cnc(jlabclc)         
          call plnumb(xlab+0.5*cs,ylab+0.5*cs,0.8*cs,float(n),0.0,-1)
        endif
c
        call newcolorname(col2)
        if(ldwashplt ) then
          call xyline(nj(n),sle(j),dwwake(j),smin,ssf,0.0,-wsf,4)
          if(llabsurf) then
            xlab = smod(sle(jlabdw))         
            ylab = -wsf*dwwake(jlabdw)         
            call plnumb(xlab+0.5*cs,ylab+0.5*cs,0.8*cs,float(n),0.0,-1)
          endif
        endif
      end do
c
c---- put up curve legends
      call newpen(4)
c
      xtmp(1) = 0.75*smax 
      xtmp(2) = 0.85*smax
      ylab = max(wsf*wmax,csf*cmax) + 12.55*cs
c
      if(lclperplt) then
       call newcolorname(col3)
       ytmp(1) = ylab/csf
       ytmp(2) = ytmp(1)
       call xyline(2,xtmp,ytmp,smin,ssf,0.0,csf,3)
       call plchar(smod(xtmp(2))+1.5*cs,csf*ytmp(2)-0.5*cs,
     &            1.2*cs,'c'  ,0., 1)
       call plmath(smod(xtmp(2))+2.5*cs,csf*ytmp(2)-0.9*cs,
     &            0.8*cs, 'v' ,0., 1)
       call plchar(smod(xtmp(2))+3.8*cs,csf*ytmp(2)-0.4*cs,
     &            0.7*cs,'t'  ,180.0, 1)
      endif
c
      call newcolorname(col4)
      ytmp(1) = ylab/csf - ysp*cs/csf
      ytmp(2) = ytmp(1)
      call xyline(2,xtmp,ytmp,smin,ssf,0.0,csf,2)
      call plchar(smod(xtmp(2))+1.5*cs,csf*ytmp(2)-0.5*cs,
     &            1.2*cs,'c'  ,0., 1)
      call plmath(smod(xtmp(2))+2.5*cs,csf*ytmp(2)-0.9*cs,
     &            0.8*cs, 'v' ,0., 1)
c
      call newcolorname(col1)
      ytmp(1) = ylab/csf - 4.4*cs/csf
      ytmp(2) = ytmp(1)
      call xyline(2,xtmp,ytmp,smin,ssf,0.0,csf,1)
      call plchar(smod(xtmp(2))+1.5*cs,csf*ytmp(2)-0.5*cs,
     &            1.2*cs,'c c/c',0.,5)
      call plchar(999.,999.,0.6*cs,'ref',0.,3)
      call plmath(smod(xtmp(2))+2.6*cs,csf*ytmp(2)-0.9*cs,
     &            0.8*cs, 'v'    ,0.,1)
c
      ytmp(1) = ylab/csf - 6.6*cs/csf
      ytmp(2) = ytmp(1)
      if(ldwashplt) then
        call newcolorname(col2)
        call xyline(2,xtmp,ytmp,smin,ssf,0.0,csf,4)
        call plmath(smod(xtmp(2))+1.5*cs,csf*ytmp(2)-0.7*cs,
     &              1.2*cs,'a',0., 1)
        call plchar(smod(xtmp(2))+2.6*cs,csf*ytmp(2)-1.1*cs,
     &              0.9*cs,'i',0., 1)
      endif
c
      call newpen(2)
      call newcolor(icol0)
c
c---- trefftz plane plot labels
      xlab2 = smod(smax) - 13.0*0.8*cs
      ylab2 = csf*cmax + 0.9*cs
      call plchar(xlab2,ylab2,0.8*cs,'trefftz plane',0.0,-1)
c
      xlab2 = smod(smax) - 8.0*0.8*cs
      ylab2 = csf*cmax + 2.5*cs
      call plchar(xlab2,ylab2,0.8*cs,'avl ',0.0,4)
      call plnumb(999.,999.,0.8*cs,version,0.0,2)
c
c      xlab = smod(smin) + 45.5*cs
c      ylab = csf*cmax + 1.0*cs
c      call plchar(xlab,ylab,0.8*cs,satype,0.0,16)
c
      xlab = smod(smin)
      ylab = max(wsf*wmax,csf*cmax) + ysp*cs*(8.0+float(nconlin))
c
c     
c---- label block with flow condition and forces located above plot
c
c---- case title, run title
 60   call newpen(3)
      call plchar(xlab,ylab,1.2*cs,title,0.0,len(title))
      ylab = ylab - ysp*cs
      if(index(rtitle(irun),'unnamed') .eq. 0) then
       call plchar(xlab,ylab,1.1*cs,rtitle(irun),0.0,len(rtitle(irun)))
      endif
c
      call newpen(2)
      xl1 = xlab
      xl2 = xlab + 15.0*cs
      xl3 = xlab + 29.0*cs
      xl4 = xlab + 45.0*cs
c
      ylab = ylab - 0.6*cs
c
c--- flow condition and forces
      rx = wrot(1)*bref/2.0
      ry = wrot(2)*cref/2.0
      rz = wrot(3)*bref/2.0
c
c--- row 1
      ylab = ylab - ysp*cs
      ylab1 = ylab
      call plmath(xl1       ,ylab,1.1*cs,'a'   ,0.0,1)
      call plchar(xl1       ,ylab,cs,'  = ',0.0,4)
      call plnumb(xl1+4.0*cs,ylab,cs, alfa/dtr,0.0,4)
c
      call plchar(xl2       ,ylab,cs,'  cl = ',0.0,7)
      call plnumb(xl2+7.0*cs,ylab,cs, cltot  ,0.0,4)
c
      call plchar(xl3       ,ylab,cs,'  cl = ',0.0,7)
      call plmath(xl3       ,ylab,cs,'    `'  ,0.0,5)
cc    call plnumb(xl3+7.0*cs,ylab,cs, dir*crtot ,0.0,4)
      call plnumb(xl3+7.0*cs,ylab,cs, dir*crsax ,0.0,4)
c
c--- row 2
      ylab = ylab - ysp*cs
      call plmath(xl1        ,ylab,1.1*cs,'b'   ,0.0,1)
      call plchar(xl1        ,ylab,cs,'  = ',0.0,4)
      call plnumb(xl1+4.0*cs,ylab,cs, beta/dtr,0.0,4)
c
      call plchar(xl2       ,ylab,cs,'  cy = ',0.0,7)
      call plnumb(xl2+7.0*cs,ylab,cs, cytot   ,0.0,4)
c
      call plchar(xl3       ,ylab,cs,'  cm = ',0.0,7)
      call plnumb(xl3+7.0*cs,ylab,cs, cmtot   ,0.0,4)
c
c--- row 3
      ylab = ylab - ysp*cs
      call plchar(xl1       ,ylab,cs,'m = ',0.0,4)
      call plnumb(xl1+4.0*cs,ylab,cs,amach ,0.0,3)
c
      call plchar(xl2       ,ylab,cs,'  cd = ',0.0,7)
      call plnumb(xl2+7.0*cs,ylab,cs, cdtot   ,0.0,5)
c
      call plchar(xl3       ,ylab,cs,'  cn = ',0.0,7)
      call plmath(xl3       ,ylab,cs,'    `'  ,0.0,5)
cc    call plnumb(xl3+7.0*cs,ylab,cs, dir*cntot, 0.0,4)
      call plnumb(xl3+7.0*cs,ylab,cs, dir*cnsax, 0.0,4)
c
c--- row 4
      ylabi = ylab
      ylabi = ylabi - ysp*cs
      call plchar(xl1       ,ylabi,cs,'pb/2v = ',0.0,8)
      call plnumb(xl1+8.0*cs,ylabi,cs, dir*rx   ,0.0,4)
c
      call plchar(xl2       ,ylabi,cs,'  cd = ',0.0,7)
      call plchar(xl2+3.8*cs,ylabi-0.4*cs,
     &                         0.8*cs,'i'      ,0.0,1)
      call plnumb(xl2+7.0*cs,ylabi,cs, cdff    ,0.0,5)
c
      call plchar(xl3       ,ylabi,cs,'   e = ',0.0,7)
      call plnumb(xl3+7.0*cs,ylabi,cs, spanef  ,0.0,4)
c
c--- row 5
      ylabi = ylabi - ysp*cs
      call plchar(xl1       ,ylabi,cs,'qc/2v = ',0.0,8)
      call plnumb(xl1+8.0*cs,ylabi,cs, ry   ,0.0,4)
c
      call plchar(xl2       ,ylabi,cs,'  cd = ',0.0,7)
      call plchar(xl2+3.8*cs,ylabi-0.4*cs,
     &                         0.8*cs,'p'      ,0.0,1)
      call plnumb(xl2+7.0*cs,ylabi,cs, cdvtot  ,0.0,5)
c
c---  row 6
      ylabi = ylabi - ysp*cs
      call plchar(xl1       ,ylabi,cs,'rb/2v = ',0.0,8)
      call plnumb(xl1+8.0*cs,ylabi,cs, dir*rz  ,0.0,4)
c
c      ylab = ylab - 0.3*cs
      ylab = ylab1 + ysp*cs

      numd = 0
      do n = 1, ncontrol
        call strip(dname(n),numdk)
        numd = max(numd,numdk)
      enddo
      do n = 1, ncontrol
        xlab = xl4
        ylab = ylab - ysp*cs
        call plchar(xlab,ylab,cs,dname(n) ,0.0,numd)
        call plchar(999.,ylab,cs,' = '    ,0.0,3)
        call plnumb(999.,ylab,cs,delcon(n),0.0,4)
      enddo
c
      call plflush
c
c*********************************************************
c
   15 lcolhc = idevh.eq.4
      write(*,1010) lclperplt, ldwashplt, llabsurf, lcolhc
   16 write(*,1030)
c
 1010 format(/' ======================================================'
     &       /'   t plot yz trace'
     &       /'   y plot data vs y'
     &       /'   z plot data vs z'
     &       /'   p erpendicular cl plot toggle (currently ',l2,')'
     &       /'   d ownwash angle   plot toggle (currently ',l2,')'
     &      //'   l imits for plot'
     &       /'   r eset plot limits'
     &      //'   n umber surfaces toggle (currently ',l2,')'
     &       /'   c olor hardcopy  toggle (currently ',l2,')'
     &       /'   a nnotate plot'
     &       /'   h ardcopy current plot'
     &      //'   zm zoom'
     &       /'   u nzoom'
     &       /'   s ize change'/)
 1030 format( ' trefftz plot command: ',$)
c
      read(*,1000) opt
      call touper(opt) 
 1000 format(a)
c
      if(opt.eq.' ') then
        call clrzoom
        call plend
        size = sizeold
        cs = csold
        write(*,*) ' '
        return
c
c---- reset plot limits
      elseif(opt.eq.'R') then
        go to 50
c
c---- set plot limits
       else if(opt.eq.'L') then
        write(*,*) 'enter new plot limits (<return> for no change)'
        rinp(1) = smin
        rinp(2) = smax
        rinp(3) = dels
        write(*,32) smin,smax,dels
 32     format('    xmin,xmax,xdel: ',3g12.6)
        call readr(3,rinp,error)
        if(error) go to 15
        smin = rinp(1)
        smax = rinp(2)
        dels = rinp(3)
c
        rinp(1) = cmin
        rinp(2) = cmax
        rinp(3) = delc
        write(*,34) cmin,cmax,delc
 34     format('    ymin,ymax,ydel: ',3g12.6)
        call readr(3,rinp,error)
        if(error) go to 15
        cmin = rinp(1)
        cmax = rinp(2)
        delc = rinp(3)
        go to 155
c
c---- plot trace
      elseif(opt.eq.'T') then
        plttype = 'trace'
        go to 100
c
c---- use y as abscissa for plot
      elseif(opt.eq.'Y') then
        plttype = 'yzcnc'
        yaxtype = 'y'
        go to 100
c
c---- use z as abscissa for plot
      elseif(opt.eq.'Z') then
        plttype = 'yzcnc'
        yaxtype = 'z'
        go to 100
c
c---- zoom in on plot
      elseif(opt.eq.'ZM') then
        call usetzoom(.false.,.true.)
        call replot(idev)
        go to 15
c
c---- reset zoom on plot
      elseif(opt.eq.'U') then
        call clrzoom
        call replot(idev)
c
c---- set plot size
      elseif(opt.eq.'S') then
   10   write(*,*)
        write(*,*) 'currently size = ',size*plfac
   12   write(*,5050)
 5050   format(' enter new value:  ',$)
        read (*,*,err=12,end=15) splf
        plfac = splf/size
        call clrzoom
        go to 100
c
c---- number loadings by surface index
      elseif(opt.eq.'N') then
        llabsurf = .not.llabsurf
        go to 100
c
c---- annotate plot
      elseif(opt.eq.'A') then
        write(*,*)
        write(*,*) '================================='
        call annot(cs)
        go to 15
c
c---- hardcopy color toggle
      elseif(opt.eq.'C') then
        lcolhc = .not. lcolhc
        if(lcolhc) then
         idevh = 4
        else
         idevh = 2
        endif
        go to 15
c
c---- hardcopy plot
      elseif(opt.eq.'H') then
        call replot(idevh)
c
c---- toggle display of cl_perp on plot
      elseif(opt.eq.'P') then
        lclperplt = .not.lclperplt
        call clrzoom
        go to 100
c
c---- toggle display of downwash on plot
      elseif(opt.eq.'D') then
        ldwashplt = .not.ldwashplt
        call clrzoom
        go to 100
c
      endif
c
      go to 16
      end ! plottp

