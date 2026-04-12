!-------------------------------------------------------------------
!   Parameter and global variables for Jet Vortex Lattice code
!
!   Primary array limits
!    nvmax   number of horseshoe vortices

!    nsmax   number of chord strips
!    nfmax   number of surfaces

!    nlmax   number of source/doublet line nodes
!    nbmax   number of bodies

!    numax   number of aircraft motion variables (u,omega)
!    ndmax   number of surface-deflection control variables
!    njmax   number of jet control variables
!    ngmax   number of design variables

!    nrmax   number of stored run cases
!    ntmax   number of stored time levels

      integer nvmax
      parameter (nvmax=5000)  ! must match NVX defined in jvl_inc.f90

      integer nsmax,nfmax,nlmax,nbmax
      integer numax,ndmax,njmax,ngmax 
      integer nrmax,ntmax

      parameter (nsmax=500, nfmax=80)
      parameter (nlmax=500, nbmax=20)
      parameter (numax=6,   ndmax=30)
      parameter (njmax=6,   ngmax=20)
      parameter (nrmax=25,  ntmax=500)
