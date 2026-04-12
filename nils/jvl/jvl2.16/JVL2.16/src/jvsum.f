
      subroutine vcsum
c--------------------------------------------------
c     Sums AIC components to get vc at all c.p.'s
c     and their total sensitivities
c
c     v(u,d,j,g) = v_gam*gam(u,d,j,g)
c     v_u        = v_gam*gam_u
c     v_d        = v_gam*gam_d
c     v_j        = v_gam*gam_j
c     v_g        = v_gam*gam_g
c--------------------------------------------------
      use jvl_inc
      include 'jvl.inc'

      do i = 1, nvor
        vc(1,i) = 0.
        vc(2,i) = 0.
        vc(3,i) = 0.
        do n = 1, numax
          vc_u(1,i,n) = 0.
          vc_u(2,i,n) = 0.
          vc_u(3,i,n) = 0.
        enddo
        do n = 1, ncontrol
          vc_d(1,i,n) = 0.
          vc_d(2,i,n) = 0.
          vc_d(3,i,n) = 0.
        enddo
        do n = 1, nvarjet
          vc_j(1,i,n) = 0.
          vc_j(2,i,n) = 0.
          vc_j(3,i,n) = 0.
        enddo
        do n = 1, ndesign
          vc_g(1,i,n) = 0.
          vc_g(2,i,n) = 0.
          vc_g(3,i,n) = 0.
        enddo

        do j = 1, nvor
          vc(1,i) = vc(1,i) + vc_gam(1,i,j)*gam(j)
          vc(2,i) = vc(2,i) + vc_gam(2,i,j)*gam(j)
          vc(3,i) = vc(3,i) + vc_gam(3,i,j)*gam(j)

          do n = 1, numax
            vc_u(1,i,n) = vc_u(1,i,n) + vc_gam(1,i,j)*gam_u(j,n)
            vc_u(2,i,n) = vc_u(2,i,n) + vc_gam(2,i,j)*gam_u(j,n)
            vc_u(3,i,n) = vc_u(3,i,n) + vc_gam(3,i,j)*gam_u(j,n)
          enddo
          do n = 1, ncontrol
            vc_d(1,i,n) = vc_d(1,i,n) + vc_gam(1,i,j)*gam_d(j,n)
            vc_d(2,i,n) = vc_d(2,i,n) + vc_gam(2,i,j)*gam_d(j,n)
            vc_d(3,i,n) = vc_d(3,i,n) + vc_gam(3,i,j)*gam_d(j,n)
          enddo
          do n = 1, nvarjet
            vc_j(1,i,n) = vc_j(1,i,n) + vc_gam(1,i,j)*gam_j(j,n)
            vc_j(2,i,n) = vc_j(2,i,n) + vc_gam(2,i,j)*gam_j(j,n)
            vc_j(3,i,n) = vc_j(3,i,n) + vc_gam(3,i,j)*gam_j(j,n)
          enddo
          do n = 1, ndesign
            vc_g(1,i,n) = vc_g(1,i,n) + vc_gam(1,i,j)*gam_g(j,n)
            vc_g(2,i,n) = vc_g(2,i,n) + vc_gam(2,i,j)*gam_g(j,n)
            vc_g(3,i,n) = vc_g(3,i,n) + vc_gam(3,i,j)*gam_g(j,n)
          enddo
        enddo
      enddo

      return
      end ! vcsum



      subroutine vvsum
c--------------------------------------------------
c     Sums AIC components to get vv at all h.v.`s
c     and their total sensitivities
c
c     v(u,d,j,g) = v_gam*gam(u,d,j,g)
c     v_u        = v_gam*gam_u
c     v_d        = v_gam*gam_d
c     v_j        = v_gam*gam_j
c     v_g        = v_gam*gam_g
c--------------------------------------------------
      use jvl_inc
      include 'jvl.inc'

      do i = 1, nvor
        vv(1,i) = 0.
        vv(2,i) = 0.
        vv(3,i) = 0.
        do n = 1, numax
          vv_u(1,i,n) = 0.
          vv_u(2,i,n) = 0.
          vv_u(3,i,n) = 0.
        enddo
        do n = 1, ncontrol
          vv_d(1,i,n) = 0.
          vv_d(2,i,n) = 0.
          vv_d(3,i,n) = 0.
        enddo
        do n = 1, nvarjet
          vv_j(1,i,n) = 0.
          vv_j(2,i,n) = 0.
          vv_j(3,i,n) = 0.
        enddo
        do n = 1, ndesign
          vv_g(1,i,n) = 0.
          vv_g(2,i,n) = 0.
          vv_g(3,i,n) = 0.
        enddo

        do j = 1, nvor
          vv(1,i) = vv(1,i) + vv_gam(1,i,j)*gam(j)
          vv(2,i) = vv(2,i) + vv_gam(2,i,j)*gam(j)
          vv(3,i) = vv(3,i) + vv_gam(3,i,j)*gam(j)

          do n = 1, numax
            vv_u(1,i,n) = vv_u(1,i,n) + vv_gam(1,i,j)*gam_u(j,n)
            vv_u(2,i,n) = vv_u(2,i,n) + vv_gam(2,i,j)*gam_u(j,n)
            vv_u(3,i,n) = vv_u(3,i,n) + vv_gam(3,i,j)*gam_u(j,n)
          enddo
          do n = 1, ncontrol
            vv_d(1,i,n) = vv_d(1,i,n) + vv_gam(1,i,j)*gam_d(j,n)
            vv_d(2,i,n) = vv_d(2,i,n) + vv_gam(2,i,j)*gam_d(j,n)
            vv_d(3,i,n) = vv_d(3,i,n) + vv_gam(3,i,j)*gam_d(j,n)
          enddo
          do n = 1, nvarjet
            vv_j(1,i,n) = vv_j(1,i,n) + vv_gam(1,i,j)*gam_j(j,n)
            vv_j(2,i,n) = vv_j(2,i,n) + vv_gam(2,i,j)*gam_j(j,n)
            vv_j(3,i,n) = vv_j(3,i,n) + vv_gam(3,i,j)*gam_j(j,n)
          enddo
          do n = 1, ndesign
            vv_g(1,i,n) = vv_g(1,i,n) + vv_gam(1,i,j)*gam_g(j,n)
            vv_g(2,i,n) = vv_g(2,i,n) + vv_gam(2,i,j)*gam_g(j,n)
            vv_g(3,i,n) = vv_g(3,i,n) + vv_gam(3,i,j)*gam_g(j,n)
          enddo
        enddo
      enddo

      return
      end ! vvsum



      subroutine vlsum
c------------------------------------------------------------------
c     Sums AIC components to get vl at all body segment points rl,
c     and their total sensitivities
c
c     v(u,d,j,g) = v_gam*gam(u,d,j,g)
c     v_u        = v_gam*gam_u
c     v_d        = v_gam*gam_d
c     v_j        = v_gam*gam_j
c     v_g        = v_gam*gam_g
c------------------------------------------------------------------
      include 'jvl.inc'

      do i = 1, nlbody
        vl(1,i) = 0.
        vl(2,i) = 0.
        vl(3,i) = 0.
        do n = 1, numax
          vl_u(1,i,n) = 0.
          vl_u(2,i,n) = 0.
          vl_u(3,i,n) = 0.
        enddo
        do n = 1, ncontrol
          vl_d(1,i,n) = 0.
          vl_d(2,i,n) = 0.
          vl_d(3,i,n) = 0.
        enddo
        do n = 1, nvarjet
          vl_j(1,i,n) = 0.
          vl_j(2,i,n) = 0.
          vl_j(3,i,n) = 0.
        enddo
        do n = 1, ndesign
          vl_g(1,i,n) = 0.
          vl_g(2,i,n) = 0.
          vl_g(3,i,n) = 0.
        enddo

        do j = 1, nvor
          vl(1,i) = vl(1,i) + vl_gam(1,i,j)*gam(j)
          vl(2,i) = vl(2,i) + vl_gam(2,i,j)*gam(j)
          vl(3,i) = vl(3,i) + vl_gam(3,i,j)*gam(j)

          do n = 1, numax
            vl_u(1,i,n) = vl_u(1,i,n) + vl_gam(1,i,j)*gam_u(j,n)
            vl_u(2,i,n) = vl_u(2,i,n) + vl_gam(2,i,j)*gam_u(j,n)
            vl_u(3,i,n) = vl_u(3,i,n) + vl_gam(3,i,j)*gam_u(j,n)
          enddo
          do n = 1, ncontrol
            vl_d(1,i,n) = vl_d(1,i,n) + vl_gam(1,i,j)*gam_d(j,n)
            vl_d(2,i,n) = vl_d(2,i,n) + vl_gam(2,i,j)*gam_d(j,n)
            vl_d(3,i,n) = vl_d(3,i,n) + vl_gam(3,i,j)*gam_d(j,n)
          enddo
          do n = 1, nvarjet
            vl_j(1,i,n) = vl_j(1,i,n) + vl_gam(1,i,j)*gam_j(j,n)
            vl_j(2,i,n) = vl_j(2,i,n) + vl_gam(2,i,j)*gam_j(j,n)
            vl_j(3,i,n) = vl_j(3,i,n) + vl_gam(3,i,j)*gam_j(j,n)
          enddo
          do n = 1, ndesign
            vl_g(1,i,n) = vl_g(1,i,n) + vl_gam(1,i,j)*gam_g(j,n)
            vl_g(2,i,n) = vl_g(2,i,n) + vl_gam(2,i,j)*gam_g(j,n)
            vl_g(3,i,n) = vl_g(3,i,n) + vl_gam(3,i,j)*gam_g(j,n)
          enddo
        enddo
      enddo

      return
      end ! vlsum



      subroutine wcsum
c--------------------------------------------------
c     Sums AIC components to get wc at all c.p.`s
c
c     w(u) = w_u*u
c--------------------------------------------------
      include 'jvl.inc'

      do i = 1, nvor
        do k = 1, 3
          wc(k,i) = wc_u(k,i,1)*vbar(1)
     &            + wc_u(k,i,2)*vbar(2)
     &            + wc_u(k,i,3)*vbar(3)
     &            + wc_u(k,i,4)*wbar(1)
     &            + wc_u(k,i,5)*wbar(2)
     &            + wc_u(k,i,6)*wbar(3)
        enddo
      enddo

      return
      end ! wcsum



      subroutine wvsum
c--------------------------------------------------
c     Sums AIC components to get wv at all h.v.`s
c
c     w(u) = w_u*u
c--------------------------------------------------
      include 'jvl.inc'

      do i = 1, nvor
        do k = 1, 3
          wv(k,i) = wv_u(k,i,1)*vbar(1)
     &            + wv_u(k,i,2)*vbar(2)
     &            + wv_u(k,i,3)*vbar(3)
     &            + wv_u(k,i,4)*wbar(1)
     &            + wv_u(k,i,5)*wbar(2)
     &            + wv_u(k,i,6)*wbar(3)
        enddo
      enddo

      return
      end ! wvsum


      subroutine wlsum
c-------------------------------------------------------------------
c     Sums AIC components to get wl at all body segment points rl.
c
c     w(u) = w_u*u
c-------------------------------------------------------------------
      include 'jvl.inc'

      do i = 1, nlbody
        do k = 1, 3
          wl(k,i) = wl_u(k,i,1)*vbar(1)
     &            + wl_u(k,i,2)*vbar(2)
     &            + wl_u(k,i,3)*vbar(3)
     &            + wl_u(k,i,4)*wbar(1)
     &            + wl_u(k,i,5)*wbar(2)
     &            + wl_u(k,i,6)*wbar(3)
        enddo
      enddo

      return
      end ! wlsum
