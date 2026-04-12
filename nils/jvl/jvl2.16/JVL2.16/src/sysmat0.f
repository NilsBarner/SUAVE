
      subroutine sysmat0(ir,asys,bsys,rsys,nsys)
c------------------------------------------------------------------
c     Computes linearized-system matrices A,B for run case ir.
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
     &     riinv(3,3)
      real p(3), p_u(3,numax),
     &     h(3), h_u(3,numax),
     &     wxp(3), wxp_u(3,numax),
     &     wxh(3), wxh_u(3,numax),
     &     mif(3), mif_u(3,numax), mif_d(3,ndmax),
     &     rim(3), rim_u(3,numax), rim_d(3,ndmax),
     &     prf(3), prf_u(3,numax),
     &     prm(3), prm_u(3,numax),
     &     mg(3) , mg_ang(3,3),
     &     mmg(3), mmg_ang(3,3)

      real cft(3), cft_u(3,numax), cft_d(3,ndmax), cft_j(3,njmax)
      real cmt(3), cmt_u(3,numax), cmt_d(3,ndmax), cmt_j(3,njmax)

      real ang(3), tt(3,3), tt_ang(3,3,3), rt(3,3), rt_ang(3,3,3)

      real cosa_u(3), sina_u(3),
     &     cosb_u(3), sinb_u(3)

      integer icrs(3), jcrs(3)
      data icrs / 2, 3, 1 / , jcrs / 3, 1, 2 /


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

      dcl_u    = parval(ipclu,ir)
      dcl_a    = parval(ipcla,ir)
      dcl_adot = parval(ipclad,ir)

      dcd_u    = parval(ipcdu,ir)
      dcd_a    = parval(ipcda,ir)
      dcd_adot = parval(ipcdad,ir)

      dcm_u    = parval(ipcmu,ir)
      dcm_a    = parval(ipcma,ir)
      dcm_adot = parval(ipcmad,ir)

      alfa = parval(ipalfa,ir)*dtr
      beta = parval(ipbeta,ir)*dtr
      call ab2vbar
      call vbar2ab

      sina = sin(alfa)
      cosa = cos(alfa)
      sinb = sin(beta)
      cosb = cos(beta)
      do k = 1, 3
        sina_u(k) =  cosa*alfa_u(k)
        cosa_u(k) = -sina*alfa_u(k)
        sinb_u(k) =  cosb*beta_u(k)
        cosb_u(k) = -sinb*beta_u(k)
      enddo

      ang(1) = phi*dtr
      ang(2) = the*dtr
      ang(3) = psi*dtr
      call rotens3(ang,tt,tt_ang)
      call rateki3(ang,rt,rt_ang)


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
       write(*,*) '** Zero Iyy.  specify with mass file or M menu'
       lterr = .true.
      endif
      if(riner(3,3) .le. 0.0) then
       write(*,*)
       write(*,*) '** Zero Izz.  specify with mass file or M menu'
       lterr = .true.
      endif

      if(lterr) then
       write(*,*)
       write(*,*) 'Eigenmodes not computed for run case', ir
       return
      endif

c---- force and moment coefficients in JVL geometry axes
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

c---- factor for dimensionalizing JVL rotation-rate variables
      rotd = vrefd/unitl

c---- total mass and inertia tensors (real + apparent)
      do k = 1, 3
        mamat(k,1) = amass(k,1)*rho
        mamat(k,2) = amass(k,2)*rho
        mamat(k,3) = amass(k,3)*rho
        mamat(k,k) = mamat(k,k) + rmass

        rimat(k,1) = riner(k,1) + ainer(k,1)*rho
        rimat(k,2) = riner(k,2) + ainer(k,2)*rho
        rimat(k,3) = riner(k,3) + ainer(k,3)*rho
      enddo

c---- invert total tensors
      call m3inv(mamat,mainv)
      call m3inv(rimat,riinv)

c---- linear and angular momentum vectors in JVL axes
      do k = 1, 3
c                        -     -   
c------ linear momentum  P = m V 
        p(k) = -rmass*vbar(k) * vrefd
        p_u(k,1) = 0.
        p_u(k,2) = 0.
        p_u(k,3) = 0.
        p_u(k,4) = 0.
        p_u(k,5) = 0.
        p_u(k,6) = 0.
        p_u(k,k) = -rmass*vrefd
c                         -   = -   -
c------ angular momentum  H = I W + h
        h(k) = (  riner(k,1)*wbar(1)
     &          + riner(k,2)*wbar(2)
     &          + riner(k,3)*wbar(3) )*rotd  +  hrotor(k)
        h_u(k,1) = 0.
        h_u(k,2) = 0.
        h_u(k,3) = 0.
        h_u(k,4) = riner(k,1)*rotd
        h_u(k,5) = riner(k,2)*rotd
        h_u(k,6) = riner(k,3)*rotd

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
        wxp(k) = (wbar(ic)*p(jc) - wbar(jc)*p(ic))*rotd
        wxh(k) = (wbar(ic)*h(jc) - wbar(jc)*h(ic))*rotd

c------ linear momentum rate  W x P
        wxp_u(k,1) = (wbar(ic)*p_u(jc,1) - wbar(jc)*p_u(ic,1))*rotd
        wxp_u(k,2) = (wbar(ic)*p_u(jc,2) - wbar(jc)*p_u(ic,2))*rotd
        wxp_u(k,3) = (wbar(ic)*p_u(jc,3) - wbar(jc)*p_u(ic,3))*rotd
        wxp_u(k,4) = 0.
        wxp_u(k,5) = 0.
        wxp_u(k,6) = 0.
        wxp_u(k,ic+3) = wxp_u(k,ic+3) + p(jc)*rotd
        wxp_u(k,jc+3) = wxp_u(k,jc+3) - p(ic)*rotd

c------ angular momentum rate  W x H
        wxh_u(k,1) = 0.
        wxh_u(k,2) = 0.
        wxh_u(k,3) = 0.
        wxh_u(k,4) = (wbar(ic)*h_u(jc,4) - wbar(jc)*h_u(ic,4))*rotd
        wxh_u(k,5) = (wbar(ic)*h_u(jc,5) - wbar(jc)*h_u(ic,5))*rotd
        wxh_u(k,6) = (wbar(ic)*h_u(jc,6) - wbar(jc)*h_u(ic,6))*rotd
        wxh_u(k,ic+3) = wxh_u(k,ic+3) + h(jc)*rotd
        wxh_u(k,jc+3) = wxh_u(k,jc+3) - h(ic)*rotd
      enddo

c          =    -   =    -   =     - -    =     - -   =      -
c---- set  m^-1 F,  i^-1 M,  m^-1 (WxP),  i^-1 (WxH), m^-1 (mg)
      do k = 1, 3
        mif(k) = mainv(k,1)*cft(1)*qsd
     &         + mainv(k,2)*cft(2)*qsd
     &         + mainv(k,3)*cft(3)*qsd
        rim(k) = riinv(k,1)*cmt(1)*qsd*brefd
     &         + riinv(k,2)*cmt(2)*qsd*crefd
     &         + riinv(k,3)*cmt(3)*qsd*brefd
        prf(k) = mainv(k,1)*wxp(1)
     &         + mainv(k,2)*wxp(2)
     &         + mainv(k,3)*wxp(3)
        prm(k) = riinv(k,1)*wxh(1)
     &         + riinv(k,2)*wxh(2)
     &         + riinv(k,3)*wxh(3)
        mmg(k) = mainv(k,1)*mg(1)
     &         + mainv(k,2)*mg(2)
     &         + mainv(k,3)*mg(3)
        do l = 1, 3
          mmg_ang(k,l) = mainv(k,1)*mg_ang(1,l)
     &                 + mainv(k,2)*mg_ang(2,l)
     &                 + mainv(k,3)*mg_ang(3,l)
        enddo
        do iu = 1, 6
          mif_u(k,iu) = 
     &           mainv(k,1)*cft_u(1,iu)*qsd
     &         + mainv(k,2)*cft_u(2,iu)*qsd
     &         + mainv(k,3)*cft_u(3,iu)*qsd
          rim_u(k,iu) = 
     &           riinv(k,1)*cmt_u(1,iu)*qsd*brefd
     &         + riinv(k,2)*cmt_u(2,iu)*qsd*crefd
     &         + riinv(k,3)*cmt_u(3,iu)*qsd*brefd
          prf_u(k,iu) =
     &           mainv(k,1)*wxp_u(1,iu)
     &         + mainv(k,2)*wxp_u(2,iu)
     &         + mainv(k,3)*wxp_u(3,iu)
          prm_u(k,iu) =
     &           riinv(k,1)*wxh_u(1,iu)
     &         + riinv(k,2)*wxh_u(2,iu)
     &         + riinv(k,3)*wxh_u(3,iu)
        enddo

        do n = 1, ncontrol
          mif_d(k,n) = 
     &           mainv(k,1)*cft_d(1,n)*qsd
     &         + mainv(k,2)*cft_d(2,n)*qsd
     &         + mainv(k,3)*cft_d(3,n)*qsd
          rim_d(k,n) = 
     &           riinv(k,1)*cmt_d(1,n)*qsd*brefd
     &         + riinv(k,2)*cmt_d(2,n)*qsd*crefd
     &         + riinv(k,3)*cmt_d(3,n)*qsd*brefd
        enddo

c------ add additional derivatives, from viscous or mach effects, or whatever
        iu = 1
        mif_u(k,iu) = mif_u(k,iu)  -  mainv(k,1)*dcl_u*qsd
        rim_u(k,iu) = rim_u(k,iu)  -  riinv(k,2)*dcm_u*qsd*crefd

        iu = 3
        mif_u(k,iu) = mif_u(k,iu)  +  mainv(k,3)*dcl_a*qsd
        rim_u(k,iu) = rim_u(k,iu)  +  riinv(k,2)*dcm_a*qsd*crefd
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
      
c---- x-acceleration
      ieq = jeu
      k = 1
      rsys(ieq)     =   mif(k)     - prf(k)      + mmg(k)
      asys(ieq,jeu) = -(mif_u(k,1) - prf_u(k,1)) / vrefd
      asys(ieq,jev) = -(mif_u(k,2) - prf_u(k,2)) / vrefd
      asys(ieq,jew) = -(mif_u(k,3) - prf_u(k,3)) / vrefd
      asys(ieq,jep) =  (mif_u(k,4) - prf_u(k,4)) / rotd
      asys(ieq,jeq) =  (mif_u(k,5) - prf_u(k,5)) / rotd
      asys(ieq,jer) =  (mif_u(k,6) - prf_u(k,6)) / rotd
      asys(ieq,jeph)=                            + mmg_ang(k,1)
      asys(ieq,jeth)=                            + mmg_ang(k,2)
      asys(ieq,jeps)=                            + mmg_ang(k,3)
      do n = 1, ncontrol
        bsys(ieq,n) =   mif_d(k,n)
      enddo

c---- y-acceleration
      ieq = jev
      k = 2
      rsys(ieq)     =   mif(k)     - prf(k)      + mmg(k)
      asys(ieq,jeu) = -(mif_u(k,1) - prf_u(k,1)) / vrefd
      asys(ieq,jev) = -(mif_u(k,2) - prf_u(k,2)) / vrefd
      asys(ieq,jew) = -(mif_u(k,3) - prf_u(k,3)) / vrefd
      asys(ieq,jep) =  (mif_u(k,4) - prf_u(k,4)) / rotd
      asys(ieq,jeq) =  (mif_u(k,5) - prf_u(k,5)) / rotd
      asys(ieq,jer) =  (mif_u(k,6) - prf_u(k,6)) / rotd
      asys(ieq,jeph)=                            + mmg_ang(k,1)
      asys(ieq,jeth)=                            + mmg_ang(k,2)
      asys(ieq,jeps)=                            + mmg_ang(k,3)
      do n = 1, ncontrol
        bsys(ieq,n) =   mif_d(k,n)
      enddo

c---- z-acceleration
      ieq = jew
      k = 3
      rsys(ieq)     =   mif(k)     - prf(k)      + mmg(k)
      asys(ieq,jeu) = -(mif_u(k,1) - prf_u(k,1)) / vrefd
      asys(ieq,jev) = -(mif_u(k,2) - prf_u(k,2)) / vrefd
      asys(ieq,jew) = -(mif_u(k,3) - prf_u(k,3)) / vrefd
      asys(ieq,jep) =  (mif_u(k,4) - prf_u(k,4)) / rotd
      asys(ieq,jeq) =  (mif_u(k,5) - prf_u(k,5)) / rotd
      asys(ieq,jer) =  (mif_u(k,6) - prf_u(k,6)) / rotd
      asys(ieq,jeph)=                            + mmg_ang(k,1)
      asys(ieq,jeth)=                            + mmg_ang(k,2)
      asys(ieq,jeps)=                            + mmg_ang(k,3)
      do n = 1, ncontrol
        bsys(ieq,n) =   mif_d(k,n)
      enddo

c---- x-ang.accel.
      ieq = jep
      k = 1
      rsys(ieq)     =   rim(k)     - prm(k)
      asys(ieq,jeu) = -(rim_u(k,1) - prm_u(k,1)) / vrefd
      asys(ieq,jev) = -(rim_u(k,2) - prm_u(k,2)) / vrefd
      asys(ieq,jew) = -(rim_u(k,3) - prm_u(k,3)) / vrefd
      asys(ieq,jep) =  (rim_u(k,4) - prm_u(k,4)) / rotd
      asys(ieq,jeq) =  (rim_u(k,5) - prm_u(k,5)) / rotd
      asys(ieq,jer) =  (rim_u(k,6) - prm_u(k,6)) / rotd
      do n = 1, ncontrol
        bsys(ieq,n) =   rim_d(k,n)
      enddo

c---- y-ang.accel.
      ieq = jeq
      k = 2
      rsys(ieq)     =   rim(k)     - prm(k)
      asys(ieq,jeu) = -(rim_u(k,1) - prm_u(k,1)) / vrefd
      asys(ieq,jev) = -(rim_u(k,2) - prm_u(k,2)) / vrefd
      asys(ieq,jew) = -(rim_u(k,3) - prm_u(k,3)) / vrefd
      asys(ieq,jep) =  (rim_u(k,4) - prm_u(k,4)) / rotd
      asys(ieq,jeq) =  (rim_u(k,5) - prm_u(k,5)) / rotd
      asys(ieq,jer) =  (rim_u(k,6) - prm_u(k,6)) / rotd
      do n = 1, ncontrol
        bsys(ieq,n) =   rim_d(k,n)
      enddo

c---- z-ang.accel.
      ieq = jer
      k = 3
      rsys(ieq)     =   rim(k)     - prm(k)
      asys(ieq,jeu) = -(rim_u(k,1) - prm_u(k,1)) / vrefd
      asys(ieq,jev) = -(rim_u(k,2) - prm_u(k,2)) / vrefd
      asys(ieq,jew) = -(rim_u(k,3) - prm_u(k,3)) / vrefd
      asys(ieq,jep) =  (rim_u(k,4) - prm_u(k,4)) / rotd
      asys(ieq,jeq) =  (rim_u(k,5) - prm_u(k,5)) / rotd
      asys(ieq,jer) =  (rim_u(k,6) - prm_u(k,6)) / rotd
      do n = 1, ncontrol
        bsys(ieq,n) =   rim_d(k,n)
      enddo

c---- phi rate
      ieq = jeph
      k = 1
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

c---- theta rate
      ieq = jeth
      k = 2
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


c      write(*,*) 'h  ', h
c      write(*,*) 'wxh', wxh
c      write(*,*) 'i-1 m  ', rim
c      write(*,*) 'i-1 wxh', prm

c      write(*,*) nzmax, nsys

      return
      end ! sysmat0

