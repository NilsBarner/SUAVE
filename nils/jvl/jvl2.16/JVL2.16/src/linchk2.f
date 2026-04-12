
      subroutine linchka
c-------------------------------------------------
c     Checks alfa(vbar),beta(vbar) derivatives
c     Checks Vbar(alfa,beta) derivatives
c-------------------------------------------------
      include 'jvl.inc'

      real vbar$(3), vbar_a$(3), vbar_b$(3)
      real alfa$, alfa_u$(3)
      real beta$, beta_u$(3)

      eps = 1.0e-5

c------------------------------------------
      v1 = 1.8
      v2 = 0.6
      v3 = -0.25
c      va = sqrt(v1**2 + v2**2 + v3**2)
      va = 1.0
      vbar(1) = v1/va
      vbar(2) = v2/va
      vbar(3) = v3/va
      call vbar2ab

      alfa$ = alfa
      beta$ = beta
      do k=1, 3
        vbar$(k) = vbar(k)
        alfa_u$(k) = alfa_u(k)
        beta_u$(k) = beta_u(k)
      enddo

      do k = 1, 3

      vbar(k) = vbar$(k) + eps
      call vbar2ab
      write(*,*)
      write(*,1100) 'da/dV', k, (alfa_u(k)+alfa_u$(k))*0.5,
     &              '      '  , (alfa-alfa$)/eps
      write(*,*)
      write(*,1100) 'db/dV', k, (beta_u(k)+beta_u$(k))*0.5,
     &              '      '  , (beta-beta$)/eps
 1100 format(1x,a,i1,f15.8,
     &      /1x,a,   f15.8 )
      vbar(k) = vbar$(k)

      enddo

c------------------------------------------

      alfa = 0.4
      beta = -0.25
      call ab2vbar
      alfa$ = alfa
      beta$ = beta
      do k = 1, 3
        vbar$(k) = vbar(k)
        vbar_a$(k) = vbar_a(k)
        vbar_b$(k) = vbar_b(k)
      enddo

      alfa = alfa$ + eps
      call ab2vbar
      write(*,*) 
      write(*,1200) 'a', ((vbar_a(k)+vbar_a$(k))*0.5, k=1, 3),
     &              ' ', ((vbar(k)-vbar$(k))/eps, k=1, 3)
      alfa = alfa$

      beta = beta$ + eps
      call ab2vbar
      write(*,*) 
      write(*,1200) 'b', ((vbar_b(k)+vbar_b$(k))*0.5, k=1, 3),
     &              ' ', ((vbar(k)-vbar$(k))/eps, k=1, 3)
      beta = beta$

 1200 format(1x,'dV/d',a, 3(f15.8,2x),
     &      /1x,'    ',a, 3(f15.8,2x) )

      stop
      end ! linchka



      subroutine linchkr
c---------------------------------------------------------------------
c     Checks R(g,u,d,j,g) derivatives
c
c     Assumes baseline quantities are currently defined
c---------------------------------------------------------------------
      use jvl_inc
      include 'jvl.inc'

      real resn$(nvmax), 
     &     aicn$(nvmax,nvmax),
     &     resn_u$(nvmax,numax),
     &     resn_d$(nvmax,ndmax),
     &     resn_j$(nvmax,njmax),
     &     resn_g$(nvmax,ngmax)

      real gam$(nvmax)

      real vbar$(3), vbar_a$(3), vbar_b$(3)
      real wbar$(3)
 
      real delcon$(ndmax)
      real deljet$(njmax)
      real deldes$(ngmax)

      character*80 line
      integer iinp(10), ir(9)
      logical error

      eps = 1.0e-6
      niter = 1

      alfa =  8.0*dtr
      beta = 12.0*dtr
      call ab2vbar
      wbar(1) = 0.8
      wbar(2) = -0.7
      wbar(3) =  1.2
      do n = 1, ncontrol
        delcon(n) = 5.0
      enddo
      do n = 1, nvarjet
        deljet(n) = 0.8
      enddo
      do i = 1, nvor
        gam(i) = 0.3 + 0.2*rc(1,i)
      enddo
      call setup

c------------------------------------------
c---- save baseline solution
      do i = 1, nvor
        gam$(i) = gam(i)
        resn$(i) = resn(i)
        do j = 1, nvor
          aicn$(i,j) = aicn(i,j)
        enddo
        do n = 1, numax
          resn_u$(i,n) = resn_u(i,n)
        enddo
        do n = 1, ncontrol
          resn_d$(i,n) = resn_d(i,n)
        enddo
        do n = 1, nvarjet
          resn_j$(i,n) = resn_j(i,n)
        enddo
        do n = 1, ndesign
          resn_g$(i,n) = resn_g(i,n)
        enddo
      enddo

      do k=1, 3
        vbar$(k) = vbar(k)
        wbar$(k) = wbar(k)
      enddo
      do n = 1, ncontrol
        delcon$(n) = delcon(n)
      enddo
      do n = 1, nvarjet
        deljet$(n) = deljet(n)
      enddo
      do n = 1, ndesign
        deldes$(n) = deldes(n)
      enddo

 1100 format(/1x,a,i1,f15.8)
 1101 format( 1x,a,1x,f15.8)
 1300 format(/1x,a,i1,3f15.8)
 1301 format( 1x,a,1x,3f15.8)
 1500 format(/1x,a,i1,5f15.8)
 1501 format( 1x,a,1x,5f15.8)
 1900 format(/1x,a,i1,9f15.8)
 1901 format( 1x,a,1x,9f15.8)

      js = 2
      write(*,*) 'surf', ifrsts(js),  ilasts(js)
      write(*,*) 'ujet', ifrstu(js),  ilastu(js)
      write(*,*) 'wjet', ifrstw(js),  ilastw(js)

      write(*,*) 'Enter j, i1, i2, i3 ... '
      read(*,'(a)') line

      ninp = 10
      call getint(line,iinp,ninp,error)
      if(error) stop

      j = iinp(1)

      ni = ninp-1
      do ii = 1, ni
        ir(ii) = iinp(ii+1)
      enddo

cc---- gam(j) will be examined
c      j = 9
c      j = 69
c      write(*,*) 'rc', rc(1,j), rc(2,j), rc(3,j), enc_d(1,j,1)

c---- examine resn(i1:i2)  (no more than 9 allowed by 1900,1901 formats)
c     i1 = i-2
c     i2 = i+2


c---- ping gam(j)
      gam(j) = gam$(j) + eps
      call setup
      write(*,*) '_______________________'
      write(*,1900) 
     &  'dR/dgam',0,((aicn(ir(i),j)+aicn$(ir(i),j))*0.5,i=1,ni)
      write(*,1901) 
     &  '       ',  ((resn(ir(i))  -resn$(ir(i))  )/eps,i=1,ni)
      gam(j) = gam$(j)
    
c---- ping Vbar(.)
      write(*,*) '_______________________'
      do l = 1, 3
        vbar(l) = vbar$(l) + eps
        call setup
        write(*,1900)
     &    'dR/dV',l,((resn_u(ir(i),l)+resn_u$(ir(i),l))*0.5,i=1,ni)
        write(*,1901)
     &    '     ',  ((resn(ir(i))    -resn$(ir(i))    )/eps,i=1,ni)
        vbar(l) = vbar$(l) 
      enddo

c---- ping Wbar(.)
      write(*,*) '_______________________'
      do l = 1, 3
        wbar(l) = wbar$(l) + eps
        call setup
        write(*,1900) 
     &    'dR/dW',l,
     &         ((resn_u(ir(i),l+3)+resn_u$(ir(i),l+3))*0.5,i=1,ni)
        write(*,1901)
     &    '     ',  ((resn(ir(i))    -resn$(ir(i))    )/eps,i=1,ni)
        wbar(l) = wbar$(l) 
      enddo

c---- ping Delcon(.)
      write(*,*) '_______________________'
      do n = 1, ncontrol
        delcon(n) = delcon$(n) + eps
        call setup
        write(*,1900)
     &     'dR/dD',n,((resn_d(ir(i),n)+resn_d$(ir(i),n))*0.5,i=1,ni)
        write(*,1901)
     &     '     ',  ((resn(ir(i))    -resn$(ir(i))    )/eps,i=1,ni)
        delcon(n) = delcon$(n)
      enddo


c---- ping Deljet(.)
      write(*,*) '_______________________'
      do n = 1, nvarjet
        deljet(n) = deljet$(n) + eps
        call setup
        write(*,1900)
     &    'dR/dJ',n,((resn_j(ir(i),n)+resn_j$(ir(i),n))*0.5,i=1,ni)
        write(*,1901)
     &    '     ',  ((resn(ir(i))    -resn$(ir(i))    )/eps,i=1,ni)
        deljet(n) = deljet$(n)
      enddo

      stop
      end ! linchkr



      subroutine linchkf
c---------------------------------------------------------------------
c     Checks gam,vc,vv,ft,mt,dtot,ltot(u,d,j,g) derivatives
c
c     Assumes baseline quantities are currently defined
c---------------------------------------------------------------------
      use jvl_inc
      include 'jvl.inc'

      real vbar$(3), vbar_a$(3), vbar_b$(3)
      real wbar$(3)
 
      real delcon$(ndmax)
      real deljet$(njmax)
      real deldes$(ngmax)

      real
     & gam$(nvmax),
     & gam_u$(nvmax,numax),
     & gam_d$(nvmax,ndmax),
     & gam_g$(nvmax,ngmax),
     & gam_j$(nvmax,njmax)

      real
     & vc$(3,nvmax),
     & vc_u$(3,nvmax,numax),
     & vc_d$(3,nvmax,ndmax),
     & vc_g$(3,nvmax,ngmax),
     & vc_j$(3,nvmax,njmax)

      real
     & vv$(3,nvmax),
     & vv_u$(3,nvmax,numax),
     & vv_d$(3,nvmax,ndmax),
     & vv_g$(3,nvmax,ngmax),
     & vv_j$(3,nvmax,njmax)

      real
     & wc$(3,nvmax),
     & wc_u$(3,nvmax,numax),
     & wc_d$(3,nvmax,ndmax),
     & wc_g$(3,nvmax,ngmax),
     & wc_j$(3,nvmax,njmax)

      real
     & wv$(3,nvmax),
     & wv_u$(3,nvmax,numax),
     & wv_d$(3,nvmax,ndmax),
     & wv_g$(3,nvmax,ngmax),
     & wv_j$(3,nvmax,njmax)

      real
     & vl$(3,nlmax),
     & vl_u$(3,nlmax,numax),
     & vl_d$(3,nlmax,ndmax),
     & vl_g$(3,nlmax,ngmax),
     & vl_j$(3,nlmax,njmax)

      real dtot$,
     &     dtot_u$(numax),
     &     dtot_d$(ndmax),
     &     dtot_j$(njmax),
     &     dtot_g$(ngmax)

      real ytot$,
     &     ytot_u$(numax),
     &     ytot_d$(ndmax),
     &     ytot_j$(njmax),
     &     ytot_g$(ngmax)

      real ltot$,
     &     ltot_u$(numax),
     &     ltot_d$(ndmax),
     &     ltot_j$(njmax),
     &     ltot_g$(ngmax)

      real dff$,
     &     dff_u$(numax),
     &     dff_d$(ndmax),
     &     dff_j$(njmax),
     &     dff_g$(ngmax)

      real yff$,
     &     yff_u$(numax),
     &     yff_d$(ndmax),
     &     yff_j$(njmax),
     &     yff_g$(ngmax)

      real lff$,
     &     lff_u$(numax),
     &     lff_d$(ndmax),
     &     lff_j$(njmax),
     &     lff_g$(ngmax)


      real ft$(3),
     &     ft_u$(3,numax),
     &     ft_d$(3,ndmax),
     &     ft_j$(3,njmax),
     &     ft_g$(3,ngmax)
      real fti$(3),
     &     fti_u$(3,numax),
     &     fti_d$(3,ndmax),
     &     fti_j$(3,njmax),
     &     fti_g$(3,ngmax)
      real ftj$(3),
     &     ftj_u$(3,numax),
     &     ftj_d$(3,ndmax),
     &     ftj_j$(3,njmax),
     &     ftj_g$(3,ngmax)
      real ftp$(3),
     &     ftp_u$(3,numax),
     &     ftp_d$(3,ndmax),
     &     ftp_j$(3,njmax),
     &     ftp_g$(3,ngmax)


      real mt$(3),
     &     mt_u$(3,numax),
     &     mt_d$(3,ndmax),
     &     mt_j$(3,njmax),
     &     mt_g$(3,ngmax)
      real mti$(3),
     &     mti_u$(3,numax),
     &     mti_d$(3,ndmax),
     &     mti_j$(3,njmax),
     &     mti_g$(3,ngmax)
      real mtj$(3),
     &     mtj_u$(3,numax),
     &     mtj_d$(3,ndmax),
     &     mtj_j$(3,njmax),
     &     mtj_g$(3,ngmax)
      real mtp$(3),
     &     mtp_u$(3,numax),
     &     mtp_d$(3,ndmax),
     &     mtp_j$(3,njmax),
     &     mtp_g$(3,ngmax)

      integer indx(nvmax)

      niter = 1
 
      eps = 2.0e-6
      deps = eps
c      deps = 0.


      alfa =  8.0*dtr
      beta = 12.0*dtr
c      alfa = 0.   !###
c      beta = 0.   !###
      call ab2vbar
      call vbar2ab

      wbar(1) = 0.2
      wbar(2) = -0.1
      wbar(3) =  0.5
c      wbar(1) = 0.   !###
c      wbar(2) = 0.   !###
c      wbar(3) = 0.   !###

      do n = 1, ncontrol
        delcon(n) = 5.0
c        delcon(n) = 0.   !###
      enddo
      do n = 1, nvarjet
        deljet(n) = 0.8
c        deljet(n) = 0.    !###
      enddo

      do iter = 1, niter+2
        call setup
        call dgetrf(nvor,nvmax,aicn,nvmax,indx,info)
ccc        call ludcmp(nvmax,nvor,aicn,indx,work)
        do i = 1, nvor
          dgam(i) = -resn(i)
        enddo
        call dgetrs('N',nvor,1,aicn,nvmax,indx,dgam,nvmax,info)
ccc        call baksub(nvmax,nvor,aicn,indx,dgam)
        do i = 1, nvor
          gam(i) = gam(i) + dgam(i)
        enddo

        dgmax = 0.
        do i = 1, nvor
          dgmax = max( abs(dgam(i)) , dgmax )
        enddo
        write(*,*) 'it', iter, dgmax
      enddo ! next iter


c      do i = 1, nvor
c      write(*,*) i, gam(i)
c      enddo
c      pause

      call setup
      call dgetrf(nvor,nvmax,aicn,nvmax,indx,info)
ccc      call ludcmp(nvmax,nvor,aicn,indx,work)
      do i = 1, nvor
        do n = 1, numax
          gam_u(i,n) = -resn_u(i,n)
        enddo
        do n = 1, ncontrol
          gam_d(i,n) = -resn_d(i,n)
        enddo
        do n = 1, nvarjet
          gam_j(i,n) = -resn_j(i,n)
        enddo
        do n = 1, ndesign
          gam_g(i,n) = -resn_g(i,n)
        enddo
      enddo

      do n = 1, numax
        call dgetrs('N',nvor,1,aicn,nvmax,indx,gam_u(1,n),nvmax,info)
ccc        call baksub(nvmax,nvor,aicn,indx,gam_u(1,n))
      enddo
      do n = 1, ncontrol
        call dgetrs('N',nvor,1,aicn,nvmax,indx,gam_d(1,n),nvmax,info)
ccc        call baksub(nvmax,nvor,aicn,indx,gam_d(1,n))
      enddo
      do n = 1, nvarjet
        call dgetrs('N',nvor,1,aicn,nvmax,indx,gam_j(1,n),nvmax,info)
ccc        call baksub(nvmax,nvor,aicn,indx,gam_j(1,n))
      enddo
      do n = 1, ndesign
        call dgetrs('N',nvor,1,aicn,nvmax,indx,gam_g(1,n),nvmax,info)
ccc        call baksub(nvmax,nvor,aicn,indx,gam_g(1,n))
      enddo

      call vcsum
      call vvsum
      call vlsum
      call wcsum
      call wvsum
      call aero

c      write(*,*) 'dyl', dtot, ytot, ltot
c      write(*,*) 'ft ', ft
c      write(*,*) 'mt ', mt
c      write(*,*) 'fi ', fti
c      write(*,*) 'mi ', mti

c------------------------------------------
c---- save baseline solution
      do k=1, 3
        vbar$(k) = vbar(k)
        wbar$(k) = wbar(k)
      enddo
      do n = 1, ncontrol
        delcon$(n) = delcon(n)
      enddo
      do n = 1, nvarjet
        deljet$(n) = deljet(n)
      enddo
      do n = 1, ndesign
        deldes$(n) = deldes(n)
      enddo

      do i = 1, nvor
        gam$(i) = gam(i)
        do n = 1, numax
          gam_u$(i,n) = gam_u(i,n)
        enddo
        do n = 1, ncontrol
          gam_d$(i,n) = gam_d(i,n)
        enddo
        do n = 1, nvarjet
          gam_j$(i,n) = gam_j(i,n)
        enddo
        do n = 1, ndesign
          gam_g$(i,n) = gam_g(i,n)
        enddo

        do k = 1, 3
          vc$(k,i) = vc(k,i)
          vv$(k,i) = vv(k,i)
          do n = 1, numax
            vc_u$(k,i,n) = vc_u(k,i,n)
            vv_u$(k,i,n) = vv_u(k,i,n)
          enddo
          do n = 1, ncontrol
            vc_d$(k,i,n) = vc_d(k,i,n)
            vv_d$(k,i,n) = vv_d(k,i,n)
          enddo
          do n = 1, nvarjet
            vc_j$(k,i,n) = vc_j(k,i,n)
            vv_j$(k,i,n) = vv_j(k,i,n)
          enddo
          do n = 1, ndesign
            vc_g$(k,i,n) = vc_g(k,i,n)
            vv_g$(k,i,n) = vv_g(k,i,n)
          enddo

          wc$(k,i) = wc(k,i)
          wv$(k,i) = wv(k,i)
          do n = 1, numax
            wc_u$(k,i,n) = wc_u(k,i,n)
            wv_u$(k,i,n) = wv_u(k,i,n)
          enddo
        enddo
      enddo

      do i = 1, nlbody
        do k = 1, 3
          vl$(k,i) = vl(k,i)
          do n = 1, numax
            vl_u$(k,i,n) = vl_u(k,i,n)
          enddo
          do n = 1, ncontrol
            vl_d$(k,i,n) = vl_d(k,i,n)
          enddo
          do n = 1, nvarjet
            vl_j$(k,i,n) = vl_j(k,i,n)
          enddo
          do n = 1, ndesign
            vl_g$(k,i,n) = vl_g(k,i,n)
          enddo
        enddo
      enddo

      do k = 1, 3
        ft$(k) = ft(k)
        fti$(k) = fti(k)
        ftj$(k) = ftj(k)
        ftp$(k) = ftp(k)
        mt$(k) = mt(k)
        mti$(k) = mti(k)
        mtj$(k) = mtj(k)
        mtp$(k) = mtp(k)
        do n = 1, numax
          ft_u$(k,n) = ft_u(k,n)
          fti_u$(k,n) = fti_u(k,n)
          ftj_u$(k,n) = ftj_u(k,n)
          ftp_u$(k,n) = ftp_u(k,n)
          mt_u$(k,n) = mt_u(k,n)
          mti_u$(k,n) = mti_u(k,n)
          mtj_u$(k,n) = mtj_u(k,n)
          mtp_u$(k,n) = mtp_u(k,n)
        enddo
        do n = 1, ncontrol
          ft_d$(k,n) = ft_d(k,n)
          fti_d$(k,n) = fti_d(k,n)
          ftj_d$(k,n) = ftj_d(k,n)
          ftp_d$(k,n) = ftp_d(k,n)
          mt_d$(k,n) = mt_d(k,n)
          mti_d$(k,n) = mti_d(k,n)
          mtj_d$(k,n) = mtj_d(k,n)
          mtp_d$(k,n) = mtp_d(k,n)
        enddo
        do n = 1, nvarjet
          ft_j$(k,n) = ft_j(k,n)
          fti_j$(k,n) = fti_j(k,n)
          ftj_j$(k,n) = ftj_j(k,n)
          ftp_j$(k,n) = ftp_j(k,n)
          mt_j$(k,n) = mt_j(k,n)
          mti_j$(k,n) = mti_j(k,n)
          mtj_j$(k,n) = mtj_j(k,n)
          mtp_j$(k,n) = mtp_j(k,n)
        enddo
        do n = 1, ndesign
          ft_g$(k,n) = ft_g(k,n)
          fti_g$(k,n) = fti_g(k,n)
          ftj_g$(k,n) = ftj_g(k,n)
          ftp_g$(k,n) = ftp_g(k,n)
          mt_g$(k,n) = mt_g(k,n)
          mti_g$(k,n) = mti_g(k,n)
          mtj_g$(k,n) = mtj_g(k,n)
          mtp_g$(k,n) = mtp_g(k,n)
        enddo

      enddo

      ltot$ = ltot
      ytot$ = ytot
      dtot$ = dtot
      do n = 1, numax
        ltot_u$(n) = ltot_u(n)
        ytot_u$(n) = ytot_u(n)
        dtot_u$(n) = dtot_u(n)
      enddo
      do n = 1, ncontrol
        ltot_d$(n) = ltot_d(n)
        ytot_d$(n) = ytot_d(n)
        dtot_d$(n) = dtot_d(n)
      enddo
      do n = 1, nvarjet
        ltot_j$(n) = ltot_j(n)
        ytot_j$(n) = ytot_j(n)
        dtot_j$(n) = dtot_j(n)
      enddo
      do n = 1, ndesign
        ltot_g$(n) = ltot_g(n)
        ytot_g$(n) = ytot_g(n)
        dtot_g$(n) = dtot_g(n)
      enddo

      lff$ = lff
      yff$ = yff
      dff$ = dff
      do n = 1, numax
        lff_u$(n) = lff_u(n)
        yff_u$(n) = yff_u(n)
        dff_u$(n) = dff_u(n)
      enddo
      do n = 1, ncontrol
        lff_d$(n) = lff_d(n)
        yff_d$(n) = yff_d(n)
        dff_d$(n) = dff_d(n)
      enddo
      do n = 1, nvarjet
        lff_j$(n) = lff_j(n)
        yff_j$(n) = yff_j(n)
        dff_j$(n) = dff_j(n)
      enddo
      do n = 1, ndesign
        lff_g$(n) = lff_g(n)
        yff_g$(n) = yff_g(n)
        dff_g$(n) = dff_g(n)
      enddo

 1100 format(/1x,a,i1,f15.8)
 1101 format( 1x,a,1x,f15.8)
 1300 format(/1x,a,i1,3f15.8)
 1301 format( 1x,a,1x,3f15.8)
 1500 format(/1x,a,i1,5f15.8)
 1501 format( 1x,a,1x,5f15.8)
 1900 format(/1x,a,i1,9f15.8)
 1901 format( 1x,a,1x,9f15.8)

c---- vc(j), vv(j) will be examined
c      j = 50
      j = 12

c      write(*,*) 'rc', rc(1,j), rc(2,j), rc(3,j), enc_d(1,j,1)

c---- examine gam(j1:j2)  (no more than 5 allowed by 1500,1501 formats)
      j1 = j-2
      j2 = j+2

c---- ping Vbar(.)
      do l = 1, 3
        vbar(l) = vbar$(l) + deps
        call vbar2ab

        do iter = 1, niter
          call setup
          call dgetrf(nvor,nvmax,aicn,nvmax,indx,info)
ccc          call ludcmp(nvmax,nvor,aicn,indx,work)
          do i = 1, nvor
            dgam(i) = -resn(i)
            do n = 1, numax
              gam_u(i,n) = -resn_u(i,n)
            enddo
          enddo
          call dgetrs('N',nvor,1,aicn,nvmax,indx,dgam,nvmax,info)
ccc          call baksub(nvmax,nvor,aicn,indx,dgam)
          do n = 1, numax
           call dgetrs('N',nvor,1,aicn,nvmax,indx,gam_u(1,n),nvmax,info)
ccc            call baksub(nvmax,nvor,aicn,indx,gam_u(1,n))
          enddo
          do i = 1, nvor
            gam(i) = gam(i) + dgam(i)
          enddo

          dgmax = 0.
          do i = 1, nvor
            dgmax = max( abs(dgam(i)) , dgmax )
          enddo
c          write(*,*) iter, dgmax
        enddo ! next iter

        call vvsum
        call vcsum
        call vlsum
        call wvsum
        call wcsum
        call aero

        write(*,*) '_______________________'
        write(*,1900) 'dgam/dV',l,((gam_u(i,l)+gam_u$(i,l))*.5,i=j1,j2)
        write(*,1901) '       ',  ((gam(i)-gam$(i))/eps,i=j1,j2)
        write(*,1300) 'dvc /dV',l,((vc_u(k,j,l)+vc_u$(k,j,l))*0.5,k=1,3)
        write(*,1301) '       ',  ((vc(k,j)-vc$(k,j))/eps, k=1,3 )
        write(*,1300) 'dvv /dV',l,((vv_u(k,j,l)+vv_u$(k,j,l))*0.5,k=1,3)
        write(*,1301) '       ',  ((vv(k,j)-vv$(k,j))/eps, k=1,3 )
        write(*,1300) 'dft /dV',l,((ft_u(k,l)+ft_u$(k,l))*0.5,k=1,3)
        write(*,1301) '       ',  ((ft(k)-ft$(k))/eps, k=1,3 )
        write(*,1100) 'ddt /dV',l,(dtot_u(l)+dtot_u$(l))*0.5
        write(*,1101) '       ',  (dtot-dtot$)/eps
        write(*,1100) 'dlt /dV',l,(ltot_u(l)+ltot_u$(l))*0.5
        write(*,1101) '       ',  (ltot-ltot$)/eps
        write(*,1100) 'ddff/dV',l,(dff_u(l)+dff_u$(l))*0.5
        write(*,1101) '       ',  (dff-dff$)/eps
        write(*,1100) 'dlff/dV',l,(lff_u(l)+lff_u$(l))*0.5
        write(*,1101) '       ',  (lff-lff$)/eps
        write(*,1300) 'dmt /dV',l,((mt_u(k,l)+mt_u$(k,l))*0.5,k=1,3)
        write(*,1301) '       ',  ((mt(k)-mt$(k))/eps, k=1,3 )
        write(*,1300) 'dfi /dV',l,((fti_u(k,l)+fti_u$(k,l))*0.5,k=1,3)
        write(*,1301) '       ',  ((fti(k)-fti$(k))/eps, k=1,3 )
        write(*,1300) 'dfj /dV',l,((ftj_u(k,l)+ftj_u$(k,l))*0.5,k=1,3)
        write(*,1301) '       ',  ((ftj(k)-ftj$(k))/eps, k=1,3 )
        write(*,1300) 'dmi /dV',l,((mti_u(k,l)+mti_u$(k,l))*0.5,k=1,3)
        write(*,1301) '       ',  ((mti(k)-mti$(k))/eps, k=1,3 )
        write(*,1300) 'dmj /dV',l,((mtj_u(k,l)+mtj_u$(k,l))*0.5,k=1,3)
        write(*,1301) '       ',  ((mtj(k)-mtj$(k))/eps, k=1,3 )

        vbar(l) = vbar$(l)
        call vbar2ab
      enddo


c---- ping Wbar(.)
      do l = 1, 3
        wbar(l) = wbar$(l) + deps

        do iter = 1, niter
          call vbar2ab
          call setup
          call dgetrf(nvor,nvmax,aicn,nvmax,indx,info)
ccc          call ludcmp(nvmax,nvor,aicn,indx,work)
          do i = 1, nvor
            dgam(i) = -resn(i)
            do n = 1, numax
              gam_u(i,n) = -resn_u(i,n)
            enddo
          enddo
          call dgetrs('N',nvor,1,aicn,nvmax,indx,dgam,nvmax,info)
ccc          call baksub(nvmax,nvor,aicn,indx,dgam)
          do n = 1, numax
           call dgetrs('N',nvor,1,aicn,nvmax,indx,gam_u(1,n),nvmax,info)
ccc           call baksub(nvmax,nvor,aicn,indx,gam_u(1,n))
          enddo
          do i = 1, nvor
            gam(i) = gam(i) + dgam(i)
          enddo

          dgmax = 0.
          do i = 1, nvor
            dgmax = max( abs(dgam(i)) , dgmax )
          enddo
          write(*,*) iter, dgmax
        enddo ! next iter

        call vvsum
        call vcsum
        call vlsum
        call wvsum
        call wcsum
        call aero
        
        write(*,*) '_______________________'
        lw = l+3
        write(*,1900)'dgam/dW',l,((gam_u(i,lw)+gam_u$(i,lw))*.5,i=j1,j2)
        write(*,1901)'       ',  ((gam(i)-gam$(i))/eps,i=j1,j2)
        write(*,1300)'dvc /dW',l,((vc_u(k,j,lw)+vc_u$(k,j,lw))*.5,k=1,3)
        write(*,1301)'       ',  ((vc(k,j)-vc$(k,j))/eps, k=1,3 )
        write(*,1300)'dvv /dW',l,((vv_u(k,j,lw)+vv_u$(k,j,lw))*.5,k=1,3)
        write(*,1301)'       ',  ((vv(k,j)-vv$(k,j))/eps, k=1,3 )
        write(*,1300)'dft /dW',l,((ft_u(k,lw)+ft_u$(k,lw))*0.5,k=1,3)
        write(*,1301)'       ',  ((ft(k)-ft$(k))/eps, k=1,3 )
        write(*,1100)'ddt /dW',l,(dtot_u(lw)+dtot_u$(lw))*0.5
        write(*,1101)'       ',  (dtot-dtot$)/eps
        write(*,1100)'dlt /dW',l,(ltot_u(lw)+ltot_u$(lw))*0.5
        write(*,1101)'       ',  (ltot-ltot$)/eps
        write(*,1100)'ddff/dW',l,(dff_u(lw)+dff_u$(lw))*0.5
        write(*,1101)'       ',  (dff-dff$)/eps
        write(*,1100)'dlff/dW',l,(lff_u(lw)+lff_u$(lw))*0.5
        write(*,1101)'       ',  (lff-lff$)/eps
        write(*,1300)'dmt /dW',l,((mt_u(k,lw)+mt_u$(k,lw))*0.5,k=1,3)
        write(*,1301)'       ',  ((mt(k)-mt$(k))/eps, k=1,3 )
        write(*,1300)'dfi /dW',l,((fti_u(k,lw)+fti_u$(k,lw))*.5,k=1,3)
        write(*,1301)'       ',  ((fti(k)-fti$(k))/eps, k=1,3 )
        write(*,1300)'dfj /dW',l,((ftj_u(k,lw)+ftj_u$(k,lw))*.5,k=1,3)
        write(*,1301)'       ',  ((ftj(k)-ftj$(k))/eps, k=1,3 )
        write(*,1300)'dmi /dW',l,((mti_u(k,lw)+mti_u$(k,lw))*.5,k=1,3)
        write(*,1301)'       ',  ((mti(k)-mti$(k))/eps, k=1,3 )
        write(*,1300)'dmj /dW',l,((mtj_u(k,lw)+mtj_u$(k,lw))*.5,k=1,3)
        write(*,1301)'       ',  ((mtj(k)-mtj$(k))/eps, k=1,3 )

        wbar(l) = wbar$(l)
      enddo

c---- ping delcon(.)
      do l = 1, ncontrol
        delcon(l) = delcon$(l) + deps

        do iter = 1, niter
          call setup
          call dgetrf(nvor,nvmax,aicn,nvmax,indx,info)
ccc          call ludcmp(nvmax,nvor,aicn,indx,work)
          do i = 1, nvor
            dgam(i) = -resn(i)
            do n = 1, ncontrol
              gam_d(i,n) = -resn_d(i,n)
            enddo
          enddo
          call dgetrs('N',nvor,1,aicn,nvmax,indx,dgam,nvmax,info)
ccc          call baksub(nvmax,nvor,aicn,indx,dgam)
          do n = 1, ncontrol
           call dgetrs('N',nvor,1,aicn,nvmax,indx,gam_d(1,n),nvmax,info)
ccc           call baksub(nvmax,nvor,aicn,indx,gam_d(1,n))
          enddo
          do i = 1, nvor
            gam(i) = gam(i) + dgam(i)
          enddo

          dgmax = 0.
          do i = 1, nvor
            dgmax = max( abs(dgam(i)) , dgmax )
          enddo
          write(*,*) iter, dgmax
        enddo ! next iter

        call vvsum
        call vcsum
        call vlsum
        call wvsum
        call wcsum
        call aero

        write(*,*) '_______________________'
        write(*,1900) 'dgam/dD',l,((gam_d(i,l)+gam_d$(i,l))*0.5,i=j1,j2)
        write(*,1901) '       ',  ((gam(i)-gam$(i))/eps, i=j1,j2)
        write(*,1300) 'dvc /dD',l,((vc_d(k,j,l)+vc_d$(k,j,l))*0.5,k=1,3)
        write(*,1301) '       ',  ((vc(k,j)-vc$(k,j))/eps, k=1,3)
        write(*,1300) 'dvv /dD',l,((vv_d(k,j,l)+vv_d$(k,j,l))*0.5,k=1,3)
        write(*,1301) '       ',  ((vv(k,j)-vv$(k,j))/eps, k=1,3)
        write(*,1300) 'dft /dD',l,((ft_d(k,l)+ft_d$(k,l))*0.5, k=1,3)
        write(*,1301) '       ',  ((ft(k)-ft$(k))/eps, k=1,3)
        write(*,1100) 'ddt /dD',l,(dtot_d(l)+dtot_d$(l))*0.5
        write(*,1101) '       ',  (dtot-dtot$)/eps
        write(*,1100) 'dlt /dD',l,(ltot_d(l)+ltot_d$(l))*0.5
        write(*,1101) '       ',  (ltot-ltot$)/eps
        write(*,1100) 'ddff/dD',l,(dff_d(l)+dff_d$(l))*0.5
        write(*,1101) '       ',  (dff-dff$)/eps
        write(*,1100) 'dlff/dD',l,(lff_d(l)+lff_d$(l))*0.5
        write(*,1101) '       ',  (lff-lff$)/eps
        write(*,1300) 'dmt /dD',l,((mt_d(k,l)+mt_d$(k,l))*0.5, k=1,3)
        write(*,1301) '       ',  ((mt(k)-mt$(k))/eps, k=1,3)
        write(*,1300) 'dfi /dD',l,((fti_d(k,l)+fti_d$(k,l))*0.5,k=1,3)
        write(*,1301) '       ',  ((fti(k)-fti$(k))/eps, k=1,3)
        write(*,1300) 'dfj /dD',l,((ftj_d(k,l)+ftj_d$(k,l))*0.5,k=1,3)
        write(*,1301) '       ',  ((ftj(k)-ftj$(k))/eps, k=1,3)
        write(*,1300) 'dfp /dD',l,((ftp_d(k,l)+ftp_d$(k,l))*0.5,k=1,3)
        write(*,1301) '       ',  ((ftp(k)-ftp$(k))/eps, k=1,3)

        delcon(l) = delcon$(l)
      enddo

c---- ping deljet(.)
      do l = 1, nvarjet
        deljet(l) = deljet$(l) + deps

        do iter = 1, niter
          call setup
          call dgetrf(nvor,nvmax,aicn,nvmax,indx,info)
ccc          call ludcmp(nvmax,nvor,aicn,indx,work)
          do i = 1, nvor
            dgam(i) = -resn(i)
            do n = 1, nvarjet
              gam_j(i,n) = -resn_j(i,n)
            enddo
          enddo
          call dgetrs('N',nvor,1,aicn,nvmax,indx,dgam,nvmax,info)
ccc          call baksub(nvmax,nvor,aicn,indx,dgam)
          do n = 1, nvarjet
           call dgetrs('N',nvor,1,aicn,nvmax,indx,gam_j(1,n),nvmax,info)
ccc           call baksub(nvmax,nvor,aicn,indx,gam_j(1,n))
          enddo
          do i = 1, nvor
            gam(i) = gam(i) + dgam(i)
          enddo

          dgmax = 0.
          do i = 1, nvor
            dgmax = max( abs(dgam(i)) , dgmax )
          enddo
          write(*,*) iter, dgmax
        enddo ! next iter

        call vvsum
        call vcsum
        call vlsum
        call wvsum
        call wcsum
        call aero

        write(*,*) '_______________________'
        write(*,1900) 'dgam/dJ',l,((gam_j(i,l)+gam_j$(i,l))*0.5,i=j1,j2)
        write(*,1901) '       ',  ((gam(i)-gam$(i))/eps, i=j1,j2)
        write(*,1300) 'dvc /dJ',l,((vc_j(k,j,l)+vc_j$(k,j,l))*0.5,k=1,3)
        write(*,1301) '       ',  ((vc(k,j)-vc$(k,j))/eps, k=1,3 )
        write(*,1300) 'dvv /dJ',l,((vv_j(k,j,l)+vv_j$(k,j,l))*0.5,k=1,3)
        write(*,1301) '       ',  ((vv(k,j)-vv$(k,j))/eps, k=1,3 )
        write(*,1300) 'dft /dJ',l,((ft_j(k,l)+ft_j$(k,l))*0.5,k=1,3)
        write(*,1301) '       ',  ((ft(k)-ft$(k))/eps, k=1,3 )
        write(*,1100) 'ddt /dJ',l,(dtot_j(l)+dtot_j$(l))*0.5
        write(*,1101) '       ',  (dtot-dtot$)/eps
        write(*,1100) 'dlt /dJ',l,(ltot_j(l)+ltot_j$(l))*0.5
        write(*,1101) '       ',  (ltot-ltot$)/eps
        write(*,1100) 'ddff/dJ',l,(dff_j(l)+dff_j$(l))*0.5
        write(*,1101) '       ',  (dff-dff$)/eps
        write(*,1100) 'dlff/dJ',l,(lff_j(l)+lff_j$(l))*0.5
        write(*,1101) '       ',  (lff-lff$)/eps
        write(*,1300) 'dmt /dJ',l,((mt_j(k,l)+mt_j$(k,l))*0.5,k=1,3)
        write(*,1301) '       ',  ((mt(k)-mt$(k))/eps, k=1,3 )
        write(*,1300) 'dfi /dJ',l,((fti_j(k,l)+fti_j$(k,l))*0.5,k=1,3)
        write(*,1301) '       ',  ((fti(k)-fti$(k))/eps, k=1,3)
        write(*,1300) 'dfj /dJ',l,((ftj_j(k,l)+ftj_j$(k,l))*0.5,k=1,3)
        write(*,1301) '       ',  ((ftj(k)-ftj$(k))/eps, k=1,3)
        write(*,1300) 'dfp /dJ',l,((ftp_j(k,l)+ftp_j$(k,l))*0.5,k=1,3)
        write(*,1301) '       ',  ((ftp(k)-ftp$(k))/eps, k=1,3)

        deljet(l) = deljet$(l)
      enddo

c---- ping deldes(.)
      do l = 1, ndesign
        deldes(l) = deldes$(l) + deps

        do iter = 1, niter
          call setup
          call dgetrf(nvor,nvmax,aicn,nvmax,indx,info)
ccc          call ludcmp(nvmax,nvor,aicn,indx,work)
          do i = 1, nvor
            dgam(i) = -resn(i)
            do n = 1, ndesign
              gam_g(i,n) = -resn_g(i,n)
            enddo
          enddo
          do n = 1, ndesign
           call dgetrs('N',nvor,1,aicn,nvmax,indx,gam_g(1,n),nvmax,info)
ccc           call baksub(nvmax,nvor,aicn,indx,gam_g(1,n))
          enddo
          do i = 1, nvor
            gam(i) = gam(i) + dgam(i)
          enddo

          dgmax = 0.
          do i = 1, nvor
            dgmax = max( abs(dgam(i)) , dgmax )
          enddo
          write(*,*) iter, dgmax
        enddo ! next iter

        call vvsum
        call vcsum
        call vlsum
        call wvsum
        call wcsum
        call aero

        write(*,*) '_______________________'
        write(*,1900) 'dgam/dG',l,((gam_g(i,l)+gam_g$(i,l))*0.5,i=j1,j2)
        write(*,1901) '       ',  ((gam(i)-gam$(i))/eps, i=j1,j2)
        write(*,1300) 'dvc /dG',l,((vc_g(k,j,l)+vc_g$(k,j,l))*0.5,k=1,3)
        write(*,1301) '       ',  ((vc(k,j)-vc$(k,j))/eps, k=1,3 )
        write(*,1300) 'dvv /dG',l,((vv_g(k,j,l)+vv_g$(k,j,l))*0.5,k=1,3)
        write(*,1301) '       ',  ((vv(k,j)-vv$(k,j))/eps, k=1,3 )
        write(*,1300) 'dft /dG',l,((ft_g(k,l)+ft_g$(k,l))*0.5,k=1,3)
        write(*,1301) '       ',  ((ft(k)-ft$(k))/eps, k=1,3 )
        write(*,1100) 'ddt /dG',l,(dtot_g(l)+dtot_g$(l))*0.5
        write(*,1101) '       ',  (dtot-dtot$)/eps
        write(*,1100) 'dlt /dG',l,(ltot_g(l)+ltot_g$(l))*0.5
        write(*,1101) '       ',  (ltot-ltot$)/eps
        write(*,1100) 'ddff/dG',l,(dff_g(l)+dff_g$(l))*0.5
        write(*,1101) '       ',  (dff-dff$)/eps
        write(*,1100) 'dlff/dG',l,(lff_g(l)+lff_g$(l))*0.5
        write(*,1101) '       ',  (lff-lff$)/eps
        write(*,1300) 'dmt /dG',l,((mt_g(k,l)+mt_g$(k,l))*0.5,k=1,3)
        write(*,1301) '       ',  ((mt(k)-mt$(k))/eps, k=1,3 )

        deldes(l) = deldes$(l)
      enddo

      stop
      end ! linchkf

