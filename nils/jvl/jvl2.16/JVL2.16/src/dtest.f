
      program dtest
C---------------------------------------------------------------------------
C     Checks JVL stability and control derivatives by finite differencing.
C     Usage: 
C             % coef xxx yyy
C
C     where "xxx", "yyy" are JVL stability derivative output filenames
C     with one parameter changed slightly.
C
C---------------------------------------------------------------------------
      implicit real(a-h,o-z)
      logical found, ljet

      parameter (ndx=9,
     &           njx=9,
     &           ngx=9)

      logical ld(ndx),
     &        lj(njx),
     &        lg(ngx)

      real dvar(ndx), dvar$(ndx), dd(ndx),
     &     jvar(njx), jvar$(njx), dj(njx),
     &     gvar(ngx), gvar$(ngx), dg(ngx)

      real CLd(ndx),   CLd$(ndx),
     &     CYd(ndx),   CYd$(ndx),
     &     Crd(ndx),   Crd$(ndx),
     &     Cmd(ndx),   Cmd$(ndx),
     &     Cnd(ndx),   Cnd$(ndx),
     &     CDffd(ndx), CDffd$(ndx),
     &     ed(ndx),    ed$(ndx),
     &    DCQd(ndx),  DCQd$(ndx),
     &    DCJd(ndx),  DCJd$(ndx),
     &    DCEd(ndx),  DCEd$(ndx),
     &     CTd(ndx),   CTd$(ndx),
     &     CPd(ndx),   CPd$(ndx)    

      real CLj(njx),   CLj$(njx),
     &     CYj(njx),   CYj$(njx),
     &     Crj(njx),   Crj$(njx),
     &     Cmj(njx),   Cmj$(njx),
     &     Cnj(njx),   Cnj$(njx),
     &     CDffj(njx), CDffj$(njx),
     &     ej(njx),    ej$(njx),
     &    DCQj(njx),  DCQj$(njx),
     &    DCJj(njx),  DCJj$(njx),
     &    DCEj(njx),  DCEj$(njx),
     &     CTj(njx),   CTj$(njx),
     &     CPj(njx),   CPj$(njx)    

      real CLg(ngx),   CLg$(ngx),
     &     CYg(ngx),   CYg$(ngx),
     &     Crg(ngx),   Crg$(ngx),
     &     Cmg(ngx),   Cmg$(ngx),
     &     Cng(ngx),   Cng$(ngx),
     &     CDffg(ngx), CDffg$(ngx),
     &     eg(ngx),    eg$(ngx),
     &    DCQg(ngx),  DCQg$(ngx),
     &    DCJg(ngx),  DCJg$(ngx),
     &    DCEg(ngx),  DCEg$(ngx),
     &     CTg(ngx),   CTg$(ngx),
     &     CPg(ngx),   CPg$(ngx)    

      character*80 argp1, argp2
      character*128 line

      character*5 vkey
      character*10 dkey, cnum
c
      pi = 4.0*atan(1.0)
c
      do id = 1, ndx
        ld(id) = .false.
        DCQd(id) = 0.
        DCJd(id) = 0.
        DCEd(id) = 0.
        CTd(id) = 0.
        CPd(id) = 0.
        DCQd$(id) = 0.
        DCJd$(id) = 0.
        DCEd$(id) = 0.
        CTd$(id) = 0.
        CPd$(id) = 0.
      enddo

      do ij = 1, njx
        lj(ij) = .false.
        DCQj(ij) = 0.
        DCJj(ij) = 0.
        DCEj(ij) = 0.
        CTj(ij) = 0.
        CPj(ij) = 0.
        DCQj$(ij) = 0.
        DCJj$(ij) = 0.
        DCEj$(ij) = 0.
        CTj$(ij) = 0.
        CPj$(ij) = 0.
      enddo

      do ig = 1, ngx
        lg(ig) = .false.
        DCQg(ig) = 0.
        DCJg(ig) = 0.
        DCEg(ig) = 0.
        CTg(ig) = 0.
        CPg(ig) = 0.
        DCQg$(ig) = 0.
        DCJg$(ig) = 0.
        DCEg$(ig) = 0.
        CTg$(ig) = 0.
        CPg$(ig) = 0.
      enddo

      ljet = .false.

c--------------------------------------
c
 1000 format(a)
c
c---- get Unix command arguments
      call getarg(1,argp1)
      call getarg(2,argp2)
c
c-----------------------------------------------------------------
      open(1,file=argp1,status='old')
      iline = 0
c
c==================================================================
c---- top of file-line reading loop
 10   continue
      read(1,1000,end=11) line
      iline = iline + 1
      if(index('#!',line(1:1)).ne.0) go to 10
C
      call getval(line,'alpha =',adeg,found)
      call getval(line,'beta  =',bdeg,found)
      a = adeg*pi/180.0
      b = bdeg*pi/180.0
C
      call getval(line,'p''b/2V =',p,found)
      call getval(line,  'qc/2V =',q,found)
      call getval(line,'r''b/2V =',r,found)

      call getval(line,'Cl''tot =',Cr,found)
      call getval(line,  'Cmtot =',Cm,found)
      call getval(line,'Cn''tot =',Cn,found)

      call getval(line,'CYtot =',CY,found)
      call getval(line,'CLtot =',CL,found)
      call getval(line,'CDtot =',CD,found)
      call getval(line,'CDff  =',CDff,found)
      call getval(line,'    e =', e,found)

      call getval(line,'DCQtot =',DCQ,found)
      call getval(line,'DCJtot =',DCJ,found)
      call getval(line,'DCEtot =',DCE,found)
      call getval(line,'CTtot =',CT,found)
      call getval(line,'CPtot =',CP,found)
      if(found) ljet = .true.

      izero = ichar('0')

      do id = 1, ndx
        if(id .le. 9) then
         vkey = 'd' // char(izero+id) // '  ='
        else
         vkey = 'd' // char(izero+1) // char(izero+id-10) // ' ='
        endif
        call getval(line,vkey,dvar(id),found)
        if(found) ld(id) = .true.
      enddo

c      write(*,*) ld
c      write(*,*) dvar

      do ij = 1, njx
        if(ij .le. 9) then
         vkey = 'j' // char(izero+ij) // '  ='
        else
         vkey = 'j' // char(izero+1) // char(izero+ij-10) // ' ='
        endif
        call getval(line,vkey,jvar(ij),found)
        if(found) lj(ij) = .true.
      enddo

      do ig = 1, ngx
        if(ig .le. 9) then
         vkey = 'g' // char(izero+ig) // '  ='
        else
         vkey = 'g' // char(izero+1) // char(izero+ig-10) // ' ='
        endif
        call getval(line,vkey,gvar(ig),found)
        if(found) lg(ig) = .true.
      enddo

      call getval(line,'CLa =',CLa,found)
      call getval(line,'CYa =',CYa,found)
      call getval(line,'Cla =',Cra,found)
      call getval(line,'Cma =',Cma,found)
      call getval(line,'Cna =',Cna,found)
                
      call getval(line,'CLb =',CLb,found)
      call getval(line,'CYb =',CYb,found)
      call getval(line,'Clb =',Crb,found)
      call getval(line,'Cmb =',Cmb,found)
      call getval(line,'Cnb =',Cnb,found)
                
      call getval(line,'CLp =',CLp,found)
      call getval(line,'CYp =',CYp,found)
      call getval(line,'Clp =',Crp,found)
      call getval(line,'Cmp =',Cmp,found)
      call getval(line,'Cnp =',Cnp,found)
                
      call getval(line,'CLq =',CLq,found)
      call getval(line,'CYq =',CYq,found)
      call getval(line,'Clq =',Crq,found)
      call getval(line,'Cmq =',Cmq,found)
      call getval(line,'Cnq =',Cnq,found)
                
      call getval(line,'CLr =',CLr,found)
      call getval(line,'CYr =',CYr,found)
      call getval(line,'Clr =',Crr,found)
      call getval(line,'Cmr =',Cmr,found)
      call getval(line,'Cnr =',Cnr,found)

      izero = ichar('0')

      do id = 1, ndx
        cnum = char(izero+id) // ' =       '
        if(ld(id)) then
         dkey = 'CLd' // cnum
         call getval(line,dkey,CLd(id),found)
         dkey = 'CYd' // cnum
         call getval(line,dkey,CYd(id),found) 
         dkey = 'Cld' // cnum
         call getval(line,dkey,Crd(id),found)
         dkey = 'Cmd' // cnum
         call getval(line,dkey,Cmd(id),found)
         dkey = 'Cnd' // cnum
         call getval(line,dkey,Cnd(id),found)
         dkey = 'CDffd' // cnum
         call getval(line,dkey,CDffd(id),found)
         dkey = 'ed' // cnum
         call getval(line,dkey,ed(id),found)
         dkey = 'DCQtotd' // cnum
         call getval(line,dkey,DCQd(id),found)
         dkey = 'DCJtotd' // cnum
         call getval(line,dkey,DCJd(id),found)
         dkey = 'DCEtotd' // cnum
         call getval(line,dkey,DCEd(id),found)
         dkey = 'CTtotd' // cnum
         call getval(line,dkey,CTd(id),found)
         dkey = 'CPtotd' // cnum
         call getval(line,dkey,CPd(id),found)
        endif
      enddo

      do ij = 1, njx
        cnum = char(izero+ij) // ' =       '
        if(lj(ij)) then
         dkey = 'CLj' // cnum
         call getval(line,dkey,CLj(ij),found)
         dkey = 'CYj' // cnum
         call getval(line,dkey,CYj(ij),found) 
         dkey = 'Clj' // cnum
         call getval(line,dkey,Crj(ij),found)
         dkey = 'Cmj' // cnum
         call getval(line,dkey,Cmj(ij),found)
         dkey = 'Cnj' // cnum
         call getval(line,dkey,Cnj(ij),found)
         dkey = 'CDffj' // cnum
         call getval(line,dkey,CDffj(ij),found)
         dkey = 'ej' // cnum
         call getval(line,dkey,ej(ij),found)
         dkey = 'DCQtotj' // cnum
         call getval(line,dkey,DCQj(ij),found)
         dkey = 'DCJtotj' // cnum
         call getval(line,dkey,DCJj(ij),found)
         dkey = 'DCEtotj' // cnum
         call getval(line,dkey,DCEj(ij),found)
         dkey = 'CTtotj' // cnum
         call getval(line,dkey,CTj(ij),found)
         dkey = 'CPtotj' // cnum
         call getval(line,dkey,CPj(ij),found)
        endif
      enddo

      do ig = 1, ngx
        cnum = char(izero+ig) // ' =       '
        if(lg(ig)) then
         dkey = 'CLg' // cnum
         call getval(line,dkey,CLg(ig),found)
         dkey = 'CYg' // cnum
         call getval(line,dkey,CYg(ig),found) 
         dkey = 'Clg' // cnum
         call getval(line,dkey,Crg(ig),found)
         dkey = 'Cmg' // cnum
         call getval(line,dkey,Cmg(ig),found)
         dkey = 'Cng' // cnum
         call getval(line,dkey,Cng(ig),found)
         dkey = 'CDffg' // cnum
         call getval(line,dkey,CDffg(ig),found)
         dkey = 'eg' // cnum
         call getval(line,dkey,eg(ig),found)
         dkey = 'DCQtotg' // cnum
         call getval(line,dkey,DCQg(ig),found)
         dkey = 'DCJtotg' // cnum
         call getval(line,dkey,DCJg(ig),found)
         dkey = 'DCEtotg' // cnum
         call getval(line,dkey,DCEg(ig),found)
         dkey = 'CTtotg' // cnum
         call getval(line,dkey,CTg(ig),found)
         dkey = 'CPtotg' // cnum
         call getval(line,dkey,CPg(ig),found)
        endif
      enddo

      go to 10
C
C=============================================
C
 11   continue
      close(1)
c
c-----------------------------------------------------------------
      open(2,file=argp2,status='old')
      iline = 0
c
c==================================================================
c---- top of file-line reading loop
 20   continue
      read(2,1000,end=21) line
      iline = iline + 1
      if(index('#!',line(1:1)).ne.0) go to 20
C
      call getval(line,'alpha =',adeg,found)
      call getval(line,'beta  =',bdeg,found)
      a$ = adeg*pi/180.0
      b$ = bdeg*pi/180.0
C
      call getval(line,'p''b/2V =',p$,found)
      call getval(line,  'qc/2V =',q$,found)
      call getval(line,'r''b/2V =',r$,found)

      call getval(line,'Cl''tot =',Cr$,found)
      call getval(line,  'Cmtot =',Cm$,found)
      call getval(line,'Cn''tot =',Cn$,found)

      call getval(line,'CYtot =',CY$,found)
      call getval(line,'CLtot =',CL$,found)
      call getval(line,'CDtot =',CD$,found)
      call getval(line,'CDff  =',CDff$,found)
      call getval(line,'    e =', e$,found)

      do id = 1, ndx
        if(id .le. 9) then
         vkey = 'd' // char(izero+id) // '  ='
        else
         vkey = 'd' // char(izero+1) // char(izero+id-10) // ' ='
        endif
        call getval(line,vkey,dvar$(id),found)
        if(found) ld(id) = .true.
      enddo

      do ij = 1, njx
        if(ij .le. 9) then
         vkey = 'j' // char(izero+ij) // '  ='
        else
         vkey = 'j' // char(izero+1) // char(izero+ij-10) // ' ='
        endif
        call getval(line,vkey,jvar$(ij),found)
        if(found) lj(ij) = .true.
      enddo

      do ig = 1, ngx
        if(ig .le. 9) then
         vkey = 'g' // char(izero+ig) // '  ='
        else
         vkey = 'g' // char(izero+1) // char(izero+ig-10) // ' ='
        endif
        call getval(line,vkey,gvar$(ig),found)
        if(found) lg(ig) = .true.
      enddo

      call getval(line,'CLa =',CLa$,found)
      call getval(line,'CYa =',CYa$,found)
      call getval(line,'Cla =',Cra$,found)
      call getval(line,'Cma =',Cma$,found)
      call getval(line,'Cna =',Cna$,found)
                
      call getval(line,'CLb =',CLb$,found)
      call getval(line,'CYb =',CYb$,found)
      call getval(line,'Clb =',Crb$,found)
      call getval(line,'Cmb =',Cmb$,found)
      call getval(line,'Cnb =',Cnb$,found)
                
      call getval(line,'CLp =',CLp$,found)
      call getval(line,'CYp =',CYp$,found)
      call getval(line,'Clp =',Crp$,found)
      call getval(line,'Cmp =',Cmp$,found)
      call getval(line,'Cnp =',Cnp$,found)
                
      call getval(line,'CLq =',CLq$,found)
      call getval(line,'CYq =',CYq$,found)
      call getval(line,'Clq =',Crq$,found)
      call getval(line,'Cmq =',Cmq$,found)
      call getval(line,'Cnq =',Cnq$,found)
                
      call getval(line,'CLr =',CLr$,found)
      call getval(line,'CYr =',CYr$,found)
      call getval(line,'Clr =',Crr$,found)
      call getval(line,'Cmr =',Cmr$,found)
      call getval(line,'Cnr =',Cnr$,found)

      call getval(line,'DCQtot =',DCQ$,found)
      call getval(line,'DCJtot =',DCJ$,found)
      call getval(line,'DCEtot =',DCE$,found)
      call getval(line,'CTtot =',CT$,found)
      call getval(line,'CPtot =',CP$,found)

      izero = ichar('0')

      do id = 1, ndx
        cnum = char(izero+id) // ' =       '
        if(ld(id)) then
         dkey = 'CLd' // cnum
         call getval(line,dkey,CLd$(id),found)
         dkey = 'CYd' // cnum
         call getval(line,dkey,CYd$(id),found) 
         dkey = 'Cld' // cnum
         call getval(line,dkey,Crd$(id),found)
         dkey = 'Cmd' // cnum
         call getval(line,dkey,Cmd$(id),found)
         dkey = 'Cnd' // cnum
         call getval(line,dkey,Cnd$(id),found)
         dkey = 'CDffd' // cnum
         call getval(line,dkey,CDffd$(id),found)
         dkey = 'ed' // cnum
         call getval(line,dkey,ed$(id),found)
         dkey = 'DCQtotd' // cnum
         call getval(line,dkey,DCQd$(id),found)
         dkey = 'DCJtotd' // cnum
         call getval(line,dkey,DCJd$(id),found)
         dkey = 'DCEtotd' // cnum
         call getval(line,dkey,DCEd$(id),found)
         dkey = 'CTtotd' // cnum
         call getval(line,dkey,CTd$(id),found)
         dkey = 'CPtotd' // cnum
         call getval(line,dkey,CPd$(id),found)
        endif
      enddo

      do ij = 1, njx
        cnum = char(izero+ij) // ' =       '
        if(lj(ij)) then
         dkey = 'CLj' // cnum
         call getval(line,dkey,CLj$(ij),found)
         dkey = 'CYj' // cnum
         call getval(line,dkey,CYj$(ij),found) 
         dkey = 'Clj' // cnum
         call getval(line,dkey,Crj$(ij),found)
         dkey = 'Cmj' // cnum
         call getval(line,dkey,Cmj$(ij),found)
         dkey = 'Cnj' // cnum
         call getval(line,dkey,Cnj$(ij),found)
         dkey = 'CDffj' // cnum
         call getval(line,dkey,CDffj$(ij),found)
         dkey = 'ej' // cnum
         call getval(line,dkey,ej$(ij),found)
         dkey = 'DCQtotj' // cnum
         call getval(line,dkey,DCQj$(ij),found)
         dkey = 'DCJtotj' // cnum
         call getval(line,dkey,DCJj$(ij),found)
         dkey = 'DCEtotj' // cnum
         call getval(line,dkey,DCEj$(ij),found)
         dkey = 'CTtotj' // cnum
         call getval(line,dkey,CTj$(ij),found)
         dkey = 'CPtotj' // cnum
         call getval(line,dkey,CPj$(ij),found)
        endif
      enddo

      do ig = 1, ngx
        cnum = char(izero+ig) // ' =       '
        if(lg(ig)) then
         dkey = 'CLg' // cnum
         call getval(line,dkey,CLg$(ig),found)
         dkey = 'CYg' // cnum
         call getval(line,dkey,CYg$(ig),found) 
         dkey = 'Clg' // cnum
         call getval(line,dkey,Crg$(ig),found)
         dkey = 'Cmg' // cnum
         call getval(line,dkey,Cmg$(ig),found)
         dkey = 'Cng' // cnum
         call getval(line,dkey,Cng$(ig),found)
         dkey = 'CDffg' // cnum
         call getval(line,dkey,CDffg$(ig),found)
         dkey = 'eg' // cnum
         call getval(line,dkey,eg$(ig),found)
         dkey = 'DCQtotg' // cnum
         call getval(line,dkey,DCQg$(ig),found)
         dkey = 'DCJtotg' // cnum
         call getval(line,dkey,DCJg$(ig),found)
         dkey = 'DCEtotg' // cnum
         call getval(line,dkey,DCEg$(ig),found)
         dkey = 'CTtotg' // cnum
         call getval(line,dkey,CTg$(ig),found)
         dkey = 'CPtotg' // cnum
         call getval(line,dkey,CPg$(ig),found)
        endif
      enddo

      go to 20
C
C=============================================
c
 21   continue
      close(2)
C
c---- independent variable changes
      da = a$ - a
      db = b$ - b
      dp = p$ - p
      dq = q$ - q
      dr = r$ - r

      do id = 1, ndx
        if(ld(id)) then
         dd(id) = dvar$(id) - dvar(id)
        else
         dd(id) = 0.
        endif
      enddo

      do ij = 1, njx
        if(lj(ij)) then
         dj(ij) = jvar$(ij) - jvar(ij)
        else
         dj(ij) = 0.
        endif
      enddo

      do ig = 1, ngx
        if(lg(ig)) then
         dg(ig) = gvar$(ig) - gvar(ig)
        else
         dg(ig) = 0.
        endif
      enddo

      dsum = da + db + dp + dq + dr
      do id = 1, ndx
        if(ld(id)) dsum = dsum + dd(id)
      enddo
      do ij = 1, njx
        if(lj(ij)) dsum = dsum + dj(ij)
      enddo
      do ig = 1, ngx
        if(lg(ig)) dsum = dsum + dg(ig)
      enddo
c
c---- accumulate predicted output changes via stability derivatives
      dCL = 0.
      dCY = 0.
      dCr = 0.
      dCm = 0.
      dCn = 0.
      dCDff = 0.
      de = 0.
      dDCQ = 0.
      dDCJ = 0.
      dDCE = 0.
      dCT = 0.
      dCP = 0.

      dCL = dCL + 0.5*(CLa$+CLa)*da
      dCY = dCY + 0.5*(CYa$+CYa)*da
      dCr = dCr + 0.5*(Cra$+Cra)*da
      dCm = dCm + 0.5*(Cma$+Cma)*da
      dCn = dCn + 0.5*(Cna$+Cna)*da
                
      dCL = dCL + 0.5*(CLb$+CLb)*db
      dCY = dCY + 0.5*(CYb$+CYb)*db
      dCr = dCr + 0.5*(Crb$+Crb)*db
      dCm = dCm + 0.5*(Cmb$+Cmb)*db
      dCn = dCn + 0.5*(Cnb$+Cnb)*db
                
      dCL = dCL + 0.5*(CLp$+CLp)*dp
      dCY = dCY + 0.5*(CYp$+CYp)*dp
      dCr = dCr + 0.5*(Crp$+Crp)*dp
      dCm = dCm + 0.5*(Cmp$+Cmp)*dp
      dCn = dCn + 0.5*(Cnp$+Cnp)*dp
                
      dCL = dCL + 0.5*(CLq$+CLq)*dq
      dCY = dCY + 0.5*(CYq$+CYq)*dq
      dCr = dCr + 0.5*(Crq$+Crq)*dq
      dCm = dCm + 0.5*(Cmq$+Cmq)*dq
      dCn = dCn + 0.5*(Cnq$+Cnq)*dq
                
      dCL = dCL + 0.5*(CLr$+CLr)*dr
      dCY = dCY + 0.5*(CYr$+CYr)*dr
      dCr = dCr + 0.5*(Crr$+Crr)*dr
      dCm = dCm + 0.5*(Cmr$+Cmr)*dr
      dCn = dCn + 0.5*(Cnr$+Cnr)*dr

      do id = 1, ndx
        if(ld(id)) then
         dCL = dCL + 0.5*(CLd$(id)+CLd(id))*dd(id)
         dCY = dCY + 0.5*(CYd$(id)+CYd(id))*dd(id)
         dCr = dCr + 0.5*(Crd$(id)+Crd(id))*dd(id)
         dCm = dCm + 0.5*(Cmd$(id)+Cmd(id))*dd(id)
         dCn = dCn + 0.5*(Cnd$(id)+Cnd(id))*dd(id)
         dCDff = dCDff + 0.5*(CDffd$(id)+CDffd(id))*dd(id)
         de = de + 0.5*(ed$(id)+ed(id))*dd(id)
         dDCQ = dDCQ + 0.5*(DCQd$(id)+DCQd(id))*dd(id)
         dDCJ = dDCJ + 0.5*(DCJd$(id)+DCJd(id))*dd(id)
         dDCE = dDCE + 0.5*(DCEd$(id)+DCEd(id))*dd(id)
         dCT = dCT + 0.5*(CTd$(id)+CTd(id))*dd(id)
         dCP = dCP + 0.5*(CPd$(id)+CPd(id))*dd(id)
        endif
      enddo

      do ij = 1, njx
        if(lj(ij)) then
         dCL = dCL + 0.5*(CLj$(ij)+CLj(ij))*dj(ij)
         dCY = dCY + 0.5*(CYj$(ij)+CYj(ij))*dj(ij)
         dCr = dCr + 0.5*(Crj$(ij)+Crj(ij))*dj(ij)
         dCm = dCm + 0.5*(Cmj$(ij)+Cmj(ij))*dj(ij)
         dCn = dCn + 0.5*(Cnj$(ij)+Cnj(ij))*dj(ij)
         dCDff = dCDff + 0.5*(CDffj$(ij)+CDffj(ij))*dj(ij)
         de = de + 0.5*(ej$(ij)+ej(ij))*dj(ij)
         dDCQ = dDCQ + 0.5*(DCQj$(ij)+DCQj(ij))*dj(ij)
         dDCJ = dDCJ + 0.5*(DCJj$(ij)+DCJj(ij))*dj(ij)
         dDCE = dDCE + 0.5*(DCEj$(ij)+DCEj(ij))*dj(ij)
         dCT = dCT + 0.5*(CTj$(ij)+CTj(ij))*dj(ij)
         dCP = dCP + 0.5*(CPj$(ij)+CPj(ij))*dj(ij)
        endif
      enddo

      do ig = 1, ngx
        if(lg(ig)) then
         dCL = dCL + 0.5*(CLg$(ig)+CLg(ig))*dg(ig)
         dCY = dCY + 0.5*(CYg$(ig)+CYg(ig))*dg(ig)
         dCr = dCr + 0.5*(Crg$(ig)+Crg(ig))*dg(ig)
         dCm = dCm + 0.5*(Cmg$(ig)+Cmg(ig))*dg(ig)
         dCn = dCn + 0.5*(Cng$(ig)+Cng(ig))*dg(ig)
         dCDff = dCDff + 0.5*(CDffg$(ig)+CDffg(ig))*dg(ig)
         de = de + 0.5*(eg$(ig)+eg(ig))*dg(ig)
         dDCQ = dDCQ + 0.5*(DCQg$(ig)+DCQg(ig))*dg(ig)
         dDCJ = dDCJ + 0.5*(DCJg$(ig)+DCJg(ig))*dg(ig)
         dDCE = dDCE + 0.5*(DCEg$(ig)+DCEg(ig))*dg(ig)
         dCT = dCT + 0.5*(CTg$(ig)+CTg(ig))*dg(ig)
         dCP = dCP + 0.5*(CPg$(ig)+CPg(ig))*dg(ig)
        endif
      enddo

c---- display state where derivatives are tested
      write(*,*)
      write(*,2100) 'CL', CL
      write(*,2100) '  ', CL$
      write(*,*)
      write(*,2100) 'CY', CY
      write(*,2100) '  ', CY$
      write(*,*)
      write(*,2100) 'Cl', Cr
      write(*,2100) '  ', Cr$
      write(*,*)
      write(*,2100) 'Cm', Cm
      write(*,2100) '  ', Cm$
      write(*,*)
      write(*,2100) 'Cn', Cn
      write(*,2100) '  ', Cn$


c---- display independent variable changes
      write(*,*)
      write(*,2100) 'da',da * 180.0/pi
      write(*,2100) 'db',db * 180.0/pi
      write(*,2100) 'dp',dp 
      write(*,2100) 'dq',dq 
      write(*,2100) 'dr',dr 

 2100 format(1x,a,1x,' =',f11.7)
 2200 format(1x,a,i1,' =',f11.7)

      do id = 1, ndx
        if(ld(id)) then
         write(*,2200) 'dd',id, dd(id)
        endif
      enddo

      do ij = 1, njx
        if(lj(ij)) then
         write(*,2200) 'dj',ij, dj(ij)
        endif
      enddo

      do ig = 1, ngx
        if(lg(ig)) then
         write(*,2200) 'dg',ig, dg(ig)
        endif
      enddo

 2500 format(1x,9(4x,a,4x))
 2510 format(1x,9f12.8)

c---- compare changes from derivatives with actual changes

      write(*,*)
      write(*,2500)   ' dCL',        ' dCY'
      write(*,2510) (   dCL)/dsum, (   dCY)/dsum
      write(*,2510) (CL$-CL)/dsum, (CY$-CY)/dsum

      write(*,*)
      write(*,2500)   ' dCl',        ' dCm',        ' dCn'
      write(*,2510) (   dCr)/dsum, (   dCm)/dsum, (   dCn)/dsum
      write(*,2510) (Cr$-Cr)/dsum, (Cm$-Cm)/dsum, (Cn$-Cn)/dsum

      write(*,*)
      write(*,2500)      'dCDff',        'de '
      write(*,2510) (     dCDff)/dsum, (  de)/dsum
      write(*,2510) (CDff$-CDff)/dsum, (e$-e)/dsum

      if(ljet) then

      write(*,*)
      write(*,2500)     'dDCQ',     'dDCJ',     'dDCE'
      write(*,2510) (    dDCQ)/dsum, (    dDCJ)/dsum, (    dDCE)/dsum
      write(*,2510) (DCQ$-DCQ)/dsum, (DCJ$-DCJ)/dsum, (DCE$-DCE)/dsum

      write(*,*)
      write(*,2500 )  ' dCT',   ' dCP'
      write(*,2510) (   dCT)/dsum, (   dCP)/dsum
      write(*,2510) (CT$-CT)/dsum, (CP$-CP)/dsum

      endif

      stop
      end


      subroutine getval(line,key,val,found)
      character*(*) line, key
      logical found
c
      found = .false.

      llen = len(line)

      klen = index(key,'=')
      if(klen.eq.0) klen = len(key)

      k = index(line,key(1:klen))
      if(k.gt.0) then
c
        k1 = k + klen
        k2 =     llen
        if(k2.lt.k1) return
c
        read(line(k1:k2),*,err=10) val
        found = .true.
      endif
      return
c
 10   write(*,*) 'READ error:', line(k1:k2)
      return
      end

