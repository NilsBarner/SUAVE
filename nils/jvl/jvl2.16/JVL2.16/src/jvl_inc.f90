! Declarations for heap memory storage for JVL
! Harold Youngren 2/25/2023

module jvl_inc

! Heap array include file for large jvl AIC matrices

  ! All non-constant variables are declared as threadprivate for OpenMP
  INTEGER, PARAMETER :: NVX=5000    ! must match nvmax in jdimen.inc)

  
  REAL*8, DIMENSION(:,:), ALLOCATABLE :: aicn
!!  !$omp threadprivate(aicn)

  REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: vc_gam
!!  !$omp threadprivate(vc_gam)

  REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: vv_gam
!!  !$omp threadprivate(vv_gam)


end module jvl_inc
