C***********************************************************************
C    Module:  jplottp.f
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

      subroutine plottp
c-------------------------------------------------
c     plots circulation, cl, downwash versus y
c-------------------------------------------------
      include 'jvl.inc'
      include 'jvlplt.inc'

      character opt*2, line*80, yaxtype*1, satype*16
      real xtmp(2), ytmp(2)
      real sle(nsmax)
      logical lcolhc
      character*24 col1, col2, col3, col4, col5, col6

      real cl_tot(nsmax), cl_air(nsmax), cl_perp(nsmax), cclc(nsmax)
      real dcj(nsmax)
      real fs(3),ms(3)

      real rinp(3)
      logical error

      data  yaxtype / 'Y' /

      include 'masks.inc'

c--- plfac is scale factor adjustment for this whole plot 
c--- wfac  is scale factor for w axis (as fraction of c axis)
c--- vtitlefrac is fraction of vertical axis to use for big label block
      data plfac, wfac, vtitlefrac / 0.85, 0.5, 0.15 /

      smod(ss) = ssf*(ss - soff)

c      if(.not.lsol) then
c       write(*,*) '*** No flow solution...'
c       return
c      endif

c---- initialize plot stuff
      cs = 0.85*csize

c--- define colors used for plot data
      col1 = 'green'   ! cl c/cref
      col2 = 'blue'    ! alpha_i
      col3 = 'red'     ! cl
      col4 = 'orange'  ! cla
      col5 = 'magenta' ! cl_perp (optional)
      col6 = 'cyan'    ! Dcj

c
      call getsa(lnasa_sa,satype,dir)

c---- calculate strip quantities to be plotted
      do js = 1, nstrip
        call fstrip(js,
     &              fs,ms,
     &              caxial,cnorml, 
     &              cl ,cd,
     &              clj,cdj,
     &              clperp,clabs,
     &              cmc4strp,cmlestrp,
     &              ccl )
        cl_tot(js) = cl
        cl_air(js) = cl - clj
        cl_perp(js) = clperp
        cclc(js) = ccl/cref
        dcj(js) = djs(js)/(0.5*wstrip(js)*chord(js))
      enddo

c---- find the min/max of y,cl,ccl/cref,downwash
 50   continue
      ymin =  rv1(2,1)
      ymax =  rv1(2,1)
      zmin =  rv1(3,1)
      zmax =  rv1(3,1)
      fmin =  cclc(1)
      fmax =  cclc(1)
      cmin = 0.0
      cmax = 0.0
      wmax = 0.0
      wmin = 0.0
      do js = 1, nstrip
        iv = ifrsts(js)
        ymin = min( ymin, rv1(2,iv), rv2(2,iv), 0.0 )
        ymax = max( ymax, rv1(2,iv), rv2(2,iv), 0.0 )
        zmin = min( zmin, rv1(3,iv), rv2(3,iv), 0.0 )
        zmax = max( zmax, rv1(3,iv), rv2(3,iv), 0.0 )
        fmin = min( fmin, cclc(js), 0.0 )
        fmax = max( fmax, cclc(js), 0.0 )
c---- use either min-max or l8 norms for axis length determination
c     l8 prevents large local spikes from giving inappropriate axes
c--- l8 norm
        wmin = wmin + min( -dwwake(js)  , 0.0 ) ** 8
        wmax = wmax + max( -dwwake(js)  , 0.0 ) ** 8
c--- regular min/max
        cmin = min( cmin, cl_perp(js), cl_tot(js) , 0.0 )
        cmax = max( cmax, cl_perp(js), cl_tot(js) , 0.0 )
c        wmin = min( wmin, -dwwake(js)  , 0.0 )
c        wmax = max( wmax, -dwwake(js)  , 0.0 )
      end do

c--- process l8 sums to norms
      if(wmin .ne. 0.0) wmin = -(abs(wmin)/nstrip) ** 0.125
      if(wmax .ne. 0.0) wmax =  (abs(wmax)/nstrip) ** 0.125

c---- normalized spanload ccl/cref and cl will have the same axis, so...
      cmin = min(fmin,cmin)
      cmax = max(fmax,cmax)

      if(cmax-cmin .lt. 1.0e-5) then
       cmin = cmin - 0.0
       cmax = cmax + 0.1
      endif

      if(wmax-wmin .lt. 1.0e-5) then
       wmin = wmin - 0.0
       wmax = wmax + 0.1
      endif

      if(ymax-ymin .lt. 0.001*bref) then
       ymin = ymin - 0.1*bref
       ymax = ymax + 0.1*bref
      endif
      
      if(zmax-zmin .lt. 0.001*bref) then
       zmin = zmin - 0.05*bref
       zmax = zmax + 0.05*bref
      endif
      
c---- determine "nice" upper,lower bounds on y, z, ccl/cref, cl, 
c-    and corresponding "nice" number of annotations
      call axisadj(ymin,ymax,yspan,dely,nyann)
      call axisadj(zmin,zmax,zspan,delz,nzann)
      call axisadj(wmin,wmax,wspan,delw,nwann)
      call axisadj(cmin,cmax,cspan,delc,ncann)

c---- ratio of w and c scale factors
      woc = wfac*(cmax-cmin)/(wmax-wmin)

c---- get down to the plot

 100  continue
      if(yaxtype.eq.'Y') then
        smin = ymin
        smax = ymax
        dels  = dely
        do js = 1, nstrip
          sle(js) = rle(2,js)
        end do
       else
        smin = zmin
        smax = zmax
        dels  = delz
        do js = 1, nstrip
          sle(js) = rle(3,js)
        end do
      endif

 105  continue
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

      call pltini(idev)

      call newfactor(plfac*plsize)
      call getcolor(icol0)

c----offset the x-axis, leaving room on left for labels
      fmin = min(wsf*wmin,csf*cmin)
      call plot(14.0*cs,6.0*cs-fmin,-3)

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

c---- cl_perp, cl, cl*c/cref axis
      call newpen(3)
      xlab = smod(smin)

      if(lclperplt) then
       ylab = csf*(cmax - 0.5*delc)
       call newcolorname(col5)
       call plchar(xlab-6.5*cs,ylab-0.5*cs,1.2*cs,'c'  ,0., 1)
       call plmath(xlab-5.5*cs,ylab-0.9*cs,0.8*cs, 'V' ,0., 1)
       call plchar(xlab-4.2*cs,ylab-0.4*cs,0.7*cs,'T'  ,180.0, 1)
      endif

      ylab = csf*(cmax - 1.5*delc)
      call newcolorname(col1)
      call plchar(xlab-8.5*cs,ylab-0.5*cs,1.2*cs,'c c/c',0.,5)
      call plchar(999.,999.,0.65*cs,'ref',0.,3)
      call plmath(xlab-7.5*cs,ylab-0.9*cs,0.8*cs, 'V'    ,0.,1)
      call newcolor(icol0)
      call yaxis(0.0,csf*cmin,csf*(cmax-cmin),csf*delc,
     &           cmin,delc, 0.9*cs,-2)

      ylab = csf*(cmax - 2.5*delc)
      call newcolorname(col3)
      call plchar(xlab-6.5*cs,ylab-0.5*cs,1.2*cs,'c'  ,0., 1)
      call plmath(xlab-5.5*cs,ylab-0.9*cs,0.8*cs, 'V' ,0., 1)

      ylab = csf*(cmax - 3.5*delc)
      call newcolorname(col4)
      call plchar(xlab-6.5*cs,ylab-0.5*cs,1.2*cs,'c'  ,0., 1)
      call plmath(xlab-5.5*cs,ylab-0.9*cs,0.8*cs, 'V' ,0., 1)
      call plchar(xlab-4.5*cs,ylab-1.1*cs,0.8*cs,'air',0., 3)

      if(nvarjet .gt. 0) then
       ylab = csf*(cmax - 4.5*delc)
       call newcolorname(col6)
       call plmath(xlab-7.5*cs,ylab-0.5*cs,1.0*cs,'D'  ,0., 1)
       call plchar(xlab-6.5*cs,ylab-0.5*cs,1.2*cs,'c'  ,0., 1)
       call plchar(xlab-5.5*cs,ylab-0.9*cs,0.8*cs, 'J' ,0., 1)
      endif

c---put up a reference grid on the y/z and cl_perp, cl, cl*c/cref axes
      call newcolor(icol0)
      ny = ifix(0.1 + (smax-smin)/dels)
      nc = ifix(0.1 + (cmax-cmin)/delc)
      call plgrid(0.0,csf*cmin,ny,ssf*dels,nc,csf*delc,lmask2)

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

c---plot the cl,clc/cref and downwash curves
      call newpen(4)
      do n = 1, nsurf
        jlabclc = jfrst(n) + 0.33*(jlast(n)-jfrst(n)+1)
        jlabdw  = jfrst(n) + 0.50*(jlast(n)-jfrst(n)+1)
        jlabcl  = jfrst(n) + 0.66*(jlast(n)-jfrst(n)+1)
        jlabdcj = jfrst(n) + 0.22*(jlast(n)-jfrst(n)+1)

c        do jj = jfrst(n), jlast(n)
c          clastrpt(jj) = clastrp(jj)
c          cltstrpt(jj) = cltstrp(jj)
c        enddo

c------ plot cl
        j = jfrst(n)
        call newcolorname(col3)
        call xyline(nj(n),sle(j),cl_tot(j),smin,ssf,0.0,csf,2)
        if(llabsurf) then
          xlab = smod(sle(jlabcl))         
          ylab = csf*cl_tot(jlabcl)         
          call plnumb(xlab+0.5*cs,ylab+0.5*cs,0.8*cs,float(n),0.0,-1)
        endif

c------ plot cl_air
        call newcolorname(col4)
        call xyline(nj(n),sle(j),cl_air(j),smin,ssf,0.0,csf,2)
        if(llabsurf) then
          xlab = smod(sle(jlabcl))         
          ylab = csf*cl_air(jlabcl)         
          call plnumb(xlab+0.5*cs,ylab+0.5*cs,0.8*cs,float(n),0.0,-1)
        endif

c------ plot cl_perp
        if(lclperplt) then
         call newcolorname(col5)
         call xyline(nj(n),sle(j),cl_perp(j),smin,ssf,0.0,csf,2)
         if(llabsurf) then
           xlab = smod(sle(jlabcl))         
           ylab = csf*cl_perp(jlabcl)         
           call plnumb(xlab+0.5*cs,ylab+0.5*cs,0.8*cs,float(n),0.0,-1)
         endif
        endif

        call newcolorname(col1)
        call xyline(nj(n),sle(j),cclc(j),smin,ssf,0.0,csf,1)
        if(llabsurf) then
          xlab = smod(sle(jlabclc))         
          ylab = csf*cclc(jlabclc)         
          call plnumb(xlab+0.5*cs,ylab+0.5*cs,0.8*cs,float(n),0.0,-1)
        endif

        if(nvarjet .gt. 0) then
         call newcolorname(col6)
         call xyline(nj(n),sle(j),dcj(j),smin,ssf,0.0,csf,2)
         if(llabsurf) then
           xlab = smod(sle(jlabdcj))         
           ylab = csf*dcj(jlabdcj)         
           call plnumb(xlab+0.5*cs,ylab+0.5*cs,0.8*cs,float(n),0.0,-1)
         endif
        endif

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

c---- put up curve legends
c      call newpen(4)
c
c      xtmp(1) = 0.75*smax 
c      xtmp(2) = 0.85*smax
c      ylab = max(wsf*wmax,csf*cmax) + 14.25*cs
cc
c      if(lclperplt) then
c       call newcolorname(col3)
c       ytmp(1) = ylab/csf + 2.2*cs/csf
c       ytmp(2) = ytmp(1)
c       call xyline(2,xtmp,ytmp,smin,ssf,0.0,csf,2)
c       call plchar(smod(xtmp(2))+1.5*cs,csf*ytmp(2)-0.5*cs,
c     &            1.2*cs,'c'  ,0., 1)
c       call plmath(smod(xtmp(2))+2.5*cs,csf*ytmp(2)-0.9*cs,
c     &            0.8*cs, 'V' ,0., 1)
c       call plchar(smod(xtmp(2))+3.8*cs,csf*ytmp(2)-0.4*cs,
c     &            0.7*cs,'T'  ,180.0, 1)
c      endif
cc
c      call newcolorname(col5)
c      ytmp(1) = ylab/csf
c      ytmp(2) = ytmp(1)
c      call xyline(2,xtmp,ytmp,smin,ssf,0.0,csf,2)
cc
c      call plchar(smod(xtmp(2))+1.5*cs,csf*ytmp(2)-0.5*cs,
c     &           1.2*cs,'c'  ,0., 1)
c      call plmath(smod(xtmp(2))+2.5*cs,csf*ytmp(2)-0.9*cs,
c     &           0.8*cs, 'V' ,0., 1)
cc
c      dxl = 2.2*cs
c      call plchar(smod(xtmp(2))+1.5*cs+dxl,csf*ytmp(2)-0.5*cs,
c     &           1.2*cs,'+c'  ,0., 2)
c      dxl = dxl +1.2*cs
c      call plmath(smod(xtmp(2))+2.5*cs+dxl,csf*ytmp(2)-0.9*cs,
c     &           0.8*cs, 'V' ,0., 1)
c      call plchar(smod(xtmp(2))+3.2*cs+dxl,csf*ytmp(2)-1.1*cs,
c     &           0.7*cs,'j'  ,0., 1)
cc
c      call newcolorname(col4)
c      ytmp(1) = ylab/csf - 2.2*cs/csf
c      ytmp(2) = ytmp(1)
c      call xyline(2,xtmp,ytmp,smin,ssf,0.0,csf,2)
c      call plchar(smod(xtmp(2))+1.5*cs,csf*ytmp(2)-0.5*cs,
c     &            1.2*cs,'c'  ,0., 1)
c      call plmath(smod(xtmp(2))+2.5*cs,csf*ytmp(2)-0.9*cs,
c     &            0.8*cs, 'V' ,0., 1)
cc
c      call newcolorname(col1)
c      ytmp(1) = ylab/csf - 4.4*cs/csf
c      ytmp(2) = ytmp(1)
c      call xyline(2,xtmp,ytmp,smin,ssf,0.0,csf,1)
c      call plchar(smod(xtmp(2))+1.5*cs,csf*ytmp(2)-0.5*cs,
c     &            1.2*cs,'c c/c',0.,5)
c      call plchar(999.,999.,0.6*cs,'ref',0.,3)
c      call plmath(smod(xtmp(2))+2.6*cs,csf*ytmp(2)-0.9*cs,
c     &            0.8*cs, 'V'    ,0.,1)
cc
c      ytmp(1) = ylab/csf - 6.6*cs/csf
c      ytmp(2) = ytmp(1)
c      if(ldwashplt) then
c        call newcolorname(col2)
c        call xyline(2,xtmp,ytmp,smin,ssf,0.0,csf,4)
c        call plmath(smod(xtmp(2))+1.5*cs,csf*ytmp(2)-0.7*cs,
c     &              1.2*cs,'a',0., 1)
c        call plchar(smod(xtmp(2))+2.6*cs,csf*ytmp(2)-1.1*cs,
c     &              0.9*cs,'i',0., 1)
c      endif
cc
      call newpen(2)
      call newcolor(icol0)

c
c---- label block with flow condition and forces located above plot
      xlab2 = smod(smax) - 13.0*0.8*cs
      ylab2 = csf*cmax + 0.9*cs
      call plchar(xlab2,ylab2,0.8*cs,'Trefftz Plane',0.0,-1)

      xlab2 = smod(smax) - 8.0*0.8*cs
      ylab2 = csf*cmax + 2.5*cs
      call plchar(xlab2,ylab2,0.8*cs,'JVL ',0.0,4)
      call plnumb(999.,999.,0.8*cs,version,0.0,2)

c      xlab = smod(smin) + 45.5*cs
c      ylab = csf*cmax + 1.0*cs
c      call plchar(xlab,ylab,0.8*cs,satype,0.0,16)

c---- number of control-variable lines
ccc      nconlin = (ncontrol+1)/2
      nconlin = max( ncontrol-4 , 0 )

c---- case title
      xlab = smod(smin)
      ylab = max(wsf*wmax,csf*cmax) + 20.2*cs + 2.2*cs*float(nconlin)

      call newpen(3)
      call plchar(xlab,ylab,1.2*cs,title,0.0,len(title))
      ylab = ylab - 2.2*cs
      if(index(rtitle(irun),'unnamed') .eq. 0) then
       call plchar(xlab,ylab,1.1*cs,rtitle(irun),0.0,len(rtitle(irun)))
      endif

      call newpen(2)
      xl1 = xlab
      xl2 = xlab + 13.5*cs
      xl3 = xlab + 29.0*cs
      xl4 = xlab + 43.5*cs

      ylab = ylab - 0.6*cs
ccc   ylab = ylab - 2.2*cs

c---- flow condition and forces
      rx = dir*wbar(1)*bref/2.0
      ry =     wbar(2)*cref/2.0
      rz = dir*wbar(3)*bref/2.0

      cl = ltot/(0.5*sref)
      cy = ytot/(0.5*sref)
      cd = dtot/(0.5*sref)

      crt = dir*mt(1)/(0.5*sref*bref)
      cmt =     mt(2)/(0.5*sref*cref)
      cnt = dir*mt(3)/(0.5*sref*bref)

      cdi = dffi/(0.5*sref)
      cdv = dffv/(0.5*sref)

      dcjt = djt/(0.5*sref)
      dcqt = dqt/ sref
      ct = (djt - dqt)/(0.5*sref)
      cp = (det - dqt)/ sref

      call plmath(xl1       ,ylab,1.1*cs,'a'   ,0.0,1)
      call plchar(xl1       ,ylab,cs,'  = ',0.0,4)
      call plnumb(xl1+4.0*cs,ylab,cs, alfa/dtr,0.0,4)

      call plchar(xl2       ,ylab,cs,'pb/2V = ',0.0,8)
      call plnumb(xl2+8.0*cs,ylab,cs, rx       ,0.0,4)

      call plchar(xl3       ,ylab,cs,'  CL = ',0.0,7)
      call plnumb(xl3+7.0*cs,ylab,cs, cl      ,0.0,4)

      call plchar(xl4       ,ylab,cs,'  Cl = ',0.0,7)
      call plnumb(xl4+7.0*cs,ylab,cs, crt     ,0.0,4)

c
      ylab = ylab - 2.2*cs
      call plmath(xl1        ,ylab,1.1*cs,'b'   ,0.0,1)
      call plchar(xl1        ,ylab,cs,'  = ',0.0,4)
      call plnumb(xl1+4.0*cs,ylab,cs, beta/dtr,0.0,4)

      call plchar(xl2       ,ylab,cs,'qc/2V = ',0.0,8)
      call plnumb(xl2+8.0*cs,ylab,cs, ry       ,0.0,4)

      call plchar(xl3       ,ylab,cs,'  CY = ',0.0,7)
      call plnumb(xl3+7.0*cs,ylab,cs, cy      ,0.0,4)

      call plchar(xl4       ,ylab,cs,'  Cm = ',0.0,7)
      call plnumb(xl4+7.0*cs,ylab,cs, cmt     ,0.0,4)

c
      ylab = ylab - 2.2*cs
      call plchar(xl1       ,ylab,cs,'M = ',0.0,4)
      call plnumb(xl1+4.0*cs,ylab,cs,amach ,0.0,3)

      call plchar(xl2       ,ylab,cs,'rb/2V = ',0.0,8)
      call plnumb(xl2+8.0*cs,ylab,cs, rz       ,0.0,4)

      call plchar(xl3       ,ylab,cs,'  CD = ',0.0,7)
      call plnumb(xl3+7.0*cs,ylab,cs, cd      ,0.0,5)

      call plchar(xl4       ,ylab,cs,'  Cn = ',0.0,7)
      call plnumb(xl4+7.0*cs,ylab,cs, cnt    , 0.0,4)

c
      ylabi = ylab
      ylabi = ylabi - 2.2*cs
      call plchar(xl3       ,ylabi,cs,'  CD = ',0.0,7)
      call plchar(xl3+3.8*cs,ylabi-0.4*cs,
     &                         0.8*cs,'i'      ,0.0,1)
      call plnumb(xl3+7.0*cs,ylabi,cs, cdi     ,0.0,5)

      call plchar(xl4       ,ylabi,cs,'   e = ',0.0,7)
      call plnumb(xl4+7.0*cs,ylabi,cs, spanef  ,0.0,4)

      ylabi = ylabi - 2.2*cs
      call plchar(xl3       ,ylabi,cs,'  CD = ',0.0,7)
      call plchar(xl3+3.8*cs,ylabi-0.4*cs,
     &                         0.8*cs,'p'      ,0.0,1)
      call plnumb(xl3+7.0*cs,ylabi,cs, cdv     ,0.0,5)

      if(nvarjet .gt. 0) then
       ylabi = ylabi - 2.2*cs
       call plmath(xl3       ,ylabi,cs,' D'     ,0.0,2)
       call plchar(xl3       ,ylabi,cs,'  CJ = ',0.0,7)
       call plnumb(xl3+7.0*cs,ylabi,cs, dcjt    ,0.0,4)

       call plchar(xl4       ,ylabi,cs,'  CT = ',0.0,7)
       call plnumb(xl4+7.0*cs,ylabi,cs, ct      ,0.0,4)

c       call plchar(xl4       ,ylabi,cs,' e   = ',0.0,7)
c       call plchar(xl4+1.9*cs,ylabi-0.4*cs,0.7*cs,'vec',0.0,3)
c       call plnumb(xl4+7.0*cs,ylabi,cs, spanev  ,0.0,4)

c
       ylabi = ylabi - 2.2*cs
       call plmath(xl3       ,ylabi,cs,' D'     ,0.0,2)
       call plchar(xl3       ,ylabi,cs,'  CQ = ',0.0,7)
       call plnumb(xl3+7.0*cs,ylabi,cs, dcqt    ,0.0,4)

       call plchar(xl4       ,ylabi,cs,'  CP = ',0.0,7)
       call plnumb(xl4+7.0*cs,ylabi,cs, cp      ,0.0,4)
      endif

c
      ylab = ylab - 0.3*cs

      numd = 0
      do n = 1, ncontrol
        call bstrip(dname(n),numdk)
        numd = max(numd,numdk)
      enddo
      do n = 1, ncontrol
        xlab = xl1
        ylab = ylab - 2.2*cs
        call plchar(xlab,ylab,cs,dname(n) ,0.0,numd)
        call plchar(999.,ylab,cs,' = '    ,0.0,3)
        call plnumb(999.,ylab,cs,delcon(n),0.0,4)
      enddo

      ylab = ylab - 0.3*cs

      numj = 0
      do n = 1, nvarjet
        call bstrip(jname(n),numjk)
        numj = max(numj,numjk)
      enddo
      do n = 1, nvarjet
        xlab = xl1
        ylab = ylab - 2.2*cs
        call plchar(xlab,ylab,cs,jname(n) ,0.0,numj)
        call plchar(999.,ylab,cs,' = '    ,0.0,3)
        call plnumb(999.,ylab,cs,deljet(n),0.0,4)
      enddo

      call plflush

c*********************************************************

   15 lcolhc = idevh.eq.4
      write(*,1010) lclperplt, ldwashplt, llabsurf, lcolhc
   16 write(*,1030)

 1010 format(/' ======================================================'
     &       /'   Y plot data vs Y'
     &       /'   Z plot data vs Z'
     &       /'   P erpendicular cl plot toggle (currently ',l2,')'
     &       /'   D ownwash angle   plot toggle (currently ',l2,')'
     &      //'   L imits for plot'
     &       /'   R eset plot limits'
     &      //'   N umber surfaces toggle (currently ',l2,')'
     &       /'   C olor hardcopy  toggle (currently ',l2,')'
     &       /'   A nnotate plot'
     &       /'   H ardcopy current plot'
     &      //'   ZM zoom'
     &       /'   U nzoom'
     &       /'   S ize change'/)
 1030 format( ' Trefftz Plot command: ',$)

      read(*,1000) opt
      call touper(opt) 
 1000 format(a)

      if(opt.eq.' ') then
        call clrzoom
        call plend
        write(*,*) ' '
        return

c---- reset plot limits
      elseif(opt.eq.'R') then
        go to 50

c---- set plot limits
       else if(opt.eq.'L') then
        write(*,*) 'Enter new plot limits (<return> for no change)'
        rinp(1) = smin
        rinp(2) = smax
        rinp(3) = dels
        write(*,32) smin,smax,dels
 32     format('    Xmin,Xmax,Xdel: ',3g12.6)
        call readr(3,rinp,error)
        if(error) go to 15
        smin = rinp(1)
        smax = rinp(2)
        dels = rinp(3)

        rinp(1) = cmin
        rinp(2) = cmax
        rinp(3) = delc
        write(*,34) cmin,cmax,delc
 34     format('    Ymin,Ymax,Ydel: ',3g12.6)
        call readr(3,rinp,error)
        if(error) go to 15
        cmin = rinp(1)
        cmax = rinp(2)
        delc = rinp(3)
        go to 105

c---- use y as abscissa for plot
      elseif(opt.eq.'Y') then
        yaxtype = 'Y'
        go to 100

c---- use z as abscissa for plot
      elseif(opt.eq.'Z') then
        yaxtype = 'Z'
        go to 100

c---- zoom in on plot
      elseif(opt.eq.'ZM') then
        call usetzoom(.false.,.true.)
        call replot(idev)
        go to 15

c---- reset zoom on plot
      elseif(opt.eq.'U') then
        call clrzoom
        call replot(idev)

c---- set plot size
      elseif(opt.eq.'S') then
   10   write(*,*)
        write(*,*) 'Currently plsize = ',plsize*plfac
   12   write(*,5050)
 5050   format(' Enter new value:  ',$)
        read (*,*,err=12,end=15) splf
        plfac = splf/plsize
        call clrzoom
        go to 100

c---- number loadings by surface index
      elseif(opt.eq.'n') then
        llabsurf = .not.llabsurf
        go to 100

c---- annotate plot
      elseif(opt.eq.'A') then
        write(*,*)
        write(*,*) '================================='
        call annot(cs)
        go to 15

c---- hardcopy color toggle
      elseif(opt.eq.'C') then
        lcolhc = .not. lcolhc
        if(lcolhc) then
         idevh = 4
        else
         idevh = 2
        endif
        go to 15

c---- hardcopy plot
      elseif(opt.eq.'H') then
        call replot(idevh)

c---- toggle display of cl_perp on plot
      elseif(opt.eq.'P') then
        lclperplt = .not.lclperplt
        call clrzoom
        go to 100

c---- toggle display of downwash on plot
      elseif(opt.eq.'D') then
        ldwashplt = .not.ldwashplt
        call clrzoom
        go to 100

      endif

      go to 16
      end ! plottp

