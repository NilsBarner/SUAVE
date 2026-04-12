!=============================================================================80
! Allocates jvl AIC variables that may be too big for COMMONS
!=============================================================================80
subroutine jvlheap_init()

  use jvl_inc

! Allocate AIC variable storage

  if (.not. allocated(aicn)) then
    allocate(aicn(NVX,NVX))
    allocate(vc_gam(3,NVX,NVX))
    allocate(vv_gam(3,NVX,NVX))
  endif
end subroutine jvlheap_init

!=============================================================================80
! Free jvl AIC variable storage
!=============================================================================80
subroutine jvlheap_clean()

  use jvl_inc

! Deallocate heap storage for AIC's 

  if (allocated(aicn)) then
    deallocate(aicn)
    deallocate(vc_gam)
    deallocate(vv_gam)
  endif

end subroutine jvlheap_clean

