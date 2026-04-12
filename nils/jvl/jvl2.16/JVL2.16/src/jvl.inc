
      implicit real (a-h,m,o-z)
      include 'jindex.inc'
c-------------------------------------------------------------------
c   Array limits
      include 'jdimen.inc'
c-------------------------------------------------------------------
c   Derived array limits
c
c    itmax  number of triangles for hidden-line routines
c    ivmax  number of global variables 
c    ipmax  number of parameters
c    icmax  number of constraints
c    jemax  number of eigenmode components

      parameter (itmax=2*nvmax,
     &           ivmax=ivtot+ndmax+njmax,
     &           icmax=ictot+ndmax+njmax,
     &           ipmax=iptot,
     &           jemax=jetot )

c-------------------------------------------------------------------
c---- unit values, names, and namelengths
      character*32 
     &  unchl,unchm,uncht,unchf,unchs,unchv,uncha,unchi,unchd,unchh
      common /un_r/
     &  unitl,unitm,unitt,unitf,units,unitv,unita,uniti,unitd,unith
      common /un_c/
     &  unchl,unchm,uncht,unchf,unchs,unchv,uncha,unchi,unchd,unchh
      common /un_i/
     &  nul  ,num  ,nut  ,nuf  ,nus  ,nuv  ,nua  ,nui  ,nud  ,nuh

      character*80 fildef, frndef, fmsdef, fprdef, fevdef
      character*80 title
      character*40 stitle, btitle, rtitle
      character*16 dname, jname, gname
      character*12 varnam
      character*16 connam
      character*12 varkey
      character*3  conkey
      character*10 parnam
      character*32 parunch
      common /case_c/
     & fildef,         ! default configuration save file
     & frndef,         ! default run case save file
     & fmsdef,         ! default mass distribution file
     & fprdef,         ! default dimensional parameter file
     & fevdef,         ! default eigenvalue save file
     & title,          ! configuration title
     & stitle(nfmax),  ! surface title
     & btitle(nbmax),  ! body title
     & rtitle(nrmax),  ! run case title
     & dname(ndmax),   ! control variable name
     & jname(njmax),   ! jet control variable name
     & gname(ngmax),   ! design  variable name
     & varnam(ivmax),  ! variable   name
     & connam(icmax),  ! constraint name
     & varkey(ivmax),  ! variable   selection key
     & conkey(icmax),  ! constraint selection key
     & parnam(ipmax),  ! run case parameter name
     & parunch(ipmax)  ! run case parameter unit name

      common /case_i/
     & luinp,             ! logical unit for configuration file
     & lurun,             ! logical unit for run case file
     & luout,             ! logical unit for output dump file
     & lustd,             ! logical unit for stability deriv. dump file
     & lusys,             ! logical unit for dynamic system matrix dump file
     & iysym,             ! y image flag (0 = no image, 1 = image, -1 = anti-image)
     & izsym,             ! z image flag (0 = no image, 1 = image, -1 = anti-image)
     & nvor,              ! number of horseshoe vortices
     & nstrip,            ! number of chordwise strips
     & nsurf,             ! number of surfaces
     & nlbody,            ! number of body source+doublet line segments
     & nbody,             ! number of bodies
     & icon(ivmax,nrmax), ! index of constraint for each variable
     & nvtot,             ! total number of available variables
     & nctot,             ! total number of available constraints
     & nptot,             ! total number of available parameters
     & ncontrol,          ! number of control variables
     & nvarjet,           ! number of jet variables
     & ndesign,           ! number of design variables
     & nitmax,            ! max number of newton iterations
     & irun, nrun,        ! current run case, number of run cases stored
     & irune,             ! target run case for eigenmode calculations
     & irunt,             ! target run case for time march initial state
     & itrim(nrmax),      ! trim type used for run case (if any)
     & neigen(nrmax),     ! number of valid eigenmodes available for run case
     & neigendat(nrmax)   ! number reference data eigenvalues

      logical lgeo,lenc,
     &        laic,lsol,
     &        lvisc,lmass,
     &        lsolr,
     &        lcondef, ljetdef, ldesdef,
     &        lptot,lpsurf,lpstrp,lpelem,lpbody,lpblin,
     &        lphinge,lpderiv,
     &        lbforce,lclvjet,
     &        lnasa_sa, lsa_rates,
     &        lmwait
      logical lppar
      common /case_l/
     & lgeo,   ! t if geometry exists
     & lenc,   ! t if normal vectors exist
     & laic,   ! t if valid velocity AIC matrices exist
     & lsol,   ! t if valid solution exists
     & lvisc,  ! t if viscous profile drag terms are to be added
     & lmass,  ! t if mass data has been read in
     & lsolr(nrmax),    ! T if run case is converged
     & lcondef(ndmax),  ! t if control variable has been declared
     & ljetdef(njmax),  ! t if jet     variable has been declared
     & ldesdef(ngmax),  ! t if design  variable has been declared
     & lptot,  ! t if total   forces are to be printed
     & lpsurf, ! t if surface forces are to be printed
     & lpstrp, ! t if strip   forces are to be printed
     & lpelem, ! t if element forces are to be printed
     & lpbody, ! t if body    forces are to be printed
     & lpblin,  ! t if body    forces are to be printed
     & lphinge, ! t if hinge moments are to be printed
     & lpderiv, ! t if stability and control derivs are to be printed
     & lbforce, ! t if body forces are to be included in total forces
     & lclvjet, ! t if viscous CL includes jet forces
     & lnasa_sa,    ! t if nasa-std stability axes are to be used
     & lsa_rates,   ! t if stability-axis rates are to be used
     & lmwait,      ! t if mode display is to wait for real time
     & lppar(ipmax)          ! t if parameter value is to be plotted

      common /case_r/
     & version,    ! jvl version number
     & dtr, pi,    ! pi/180 ,  3.1415926...
     & ysym,      ! y location of y-symmetry plane (set to zero in jinput)
     & zsym,      ! z location of z-symmetry plane
     & alfa,      ! angle of attack  = atan2( vbarz / vbarx )
     & alfa_u(3), ! d(alfa)/d(vbar)
     & beta,      ! sideslip angle   = atan2( -vbary , sqrt(vbarx^2 + vinz^2) )
     & beta_u(3), ! d(beta)/d(vbar)
     & vbar(3),   ! freestream velocies in body axes
     & vbar_a(3), ! d(vbar)/d(alpha)
     & vbar_b(3), ! d(vbar)/d(beta)
     & wbar(3),             ! rotation rates in body axes
     & dtimed,              ! time step for eigenmode movie integration
     & varval(ivmax,nrmax),  ! variable  values for run cases
     & parval(ipmax,nrmax),  ! parameter values for run cases
     & conval(icmax,nrmax),  ! imposed constraint value
     & delcon(ndmax),        ! control variable value
     & deljet(njmax),        ! jet control variable value
     & deldes(ngmax),        ! design  variable value
     & gconmax(ndmax),       ! max gain for control variable
     & gjetmax(njmax),       ! max gain for jet variable
     & sref, cref, bref,     ! area, chord, span  reference values
     & xyzref(3),            ! x,y,z location for moments, rotation
     & xyzref0(3),           ! x,y,z location for moments, rotation (default)
     & xyzmom(3),            ! x,y,z location for alt. moment reference
     & xyzmom0(3),           ! x,y,z location for alt. moment reference (default)
     & hrotor(3),            ! onboard angular momentum
     & hrotor0(3),           ! onboard angular momentum (default)
     & mach,                 ! freestream mach number
     & mach0,                ! freestream mach number (default)
     & cdref,                ! baseline profile cd
     & cdref0,               ! baseline profile cd (default)
     & vrcorec, vrcorew,     ! vortex core radius factors (*chord, *strip width)
     & srcore,               ! source core radius factor  (*body segment length)
     & dvjexp,  ! exponent for speed dependence of DVjet ~ (Vinf/Vref)^dvjexp
     & dvjexp0  ! default dvjexp

c----- NOTE: All forces and moments normalized with  rhoref*Vref^2  =  2 qref

      real ltot, ltot_u, ltot_d, ltot_j, ltot_g
      real lff , lff_u , lff_d , lff_j , lff_g
      real lffi, lffi_u, lffi_d, lffi_j, lffi_g
      real lffj, lffj_u, lffj_d, lffj_j, lffj_g
      real lffb, lffb_u, lffb_d, lffb_j, lffb_g

      common /forc_r/
     & lff, yff, dff,  ! trefftz-plane  total L,Y,D / rho Vref^2
     & lff_u(numax),yff_u(numax),dff_u(numax), ! deriv wrt wbar
     & lff_d(ndmax),yff_d(ndmax),dff_d(ndmax), ! deriv wrt control
     & lff_j(njmax),yff_j(njmax),dff_j(njmax), ! deriv wrt jet control
     & lff_g(ngmax),yff_g(ngmax),dff_g(ngmax), ! deriv wrt design
     & lffi, yffi, dffi,  ! trefftz-plane  inviscid (vortex) L,Y,D / rho Vref^2
     & lffi_u(numax),yffi_u(numax),dffi_u(numax), ! deriv wrt vbar,wbar
     & lffi_d(ndmax),yffi_d(ndmax),dffi_d(ndmax), ! deriv wrt control
     & lffi_j(njmax),yffi_j(njmax),dffi_j(njmax), ! deriv wrt jet control
     & lffi_g(ngmax),yffi_g(ngmax),dffi_g(ngmax), ! deriv wrt design
     & lffj, yffj, dffj,  ! trefftz-plane  jet L,Y,D / rho Vref^2
     & lffj_u(numax),yffj_u(numax),dffj_u(numax), ! deriv wrt vbar,wbar
     & lffj_d(ndmax),yffj_d(ndmax),dffj_d(ndmax), ! deriv wrt control
     & lffj_j(njmax),yffj_j(njmax),dffj_j(njmax), ! deriv wrt jet control
     & lffj_g(ngmax),yffj_g(ngmax),dffj_g(ngmax), ! deriv wrt design
     & lffb, yffb, dffb,  ! trefftz-plane  body L,Y,D / rho Vref^2
     & lffb_u(numax),yffb_u(numax),dffb_u(numax), ! deriv wrt vbar,wbar
     & lffb_d(ndmax),yffb_d(ndmax),dffb_d(ndmax), ! deriv wrt control
     & lffb_j(njmax),yffb_j(njmax),dffb_j(njmax), ! deriv wrt jet control
     & lffb_g(ngmax),yffb_g(ngmax),dffb_g(ngmax), ! deriv wrt design
     & dffv,                ! trefftz-plane viscous D
     & dffv_u(numax),
     & dffv_d(ndmax),
     & dffv_j(njmax),
     & dffv_g(ngmax),
     & dqff,           ! trefftz-plane jet Dm / rho Vref
     & djff,           ! trefftz-plane jet DJ / rho Vref^2
     & deff,           ! trefftz-plane jet DE / rho Vref^3
     & dqff_u(numax),  ! deriv wrt freestream u,v,w
     & djff_u(numax),
     & deff_u(numax),
     & dqff_j(njmax),  ! deriv wrt jet control
     & djff_j(njmax),
     & deff_j(njmax),
     & spanef,           ! span efficiency 
     & spanef_u(numax),  ! d(spanef)/d(vbar,wbar)
     & spanef_d(ndmax),  ! d(spanef)/d(control)
     & spanef_j(njmax),  ! d(spanef)/d(jet)
     & spanef_g(ngmax),  ! d(spanef)/d(design)
     & spanev,           ! span efficiency of vectored jet
     & spanev_u(numax),  ! d(spanev)/d(vbar,wbar)
     & spanev_d(ndmax),  ! d(spanev)/d(control)
     & spanev_j(njmax),  ! d(spanev)/d(jet control)
     & spanev_g(ngmax),  ! d(spanev)/d(design)
     & dtot, ytot, ltot,          ! total D,Y,L
     & dtot_u(numax),ytot_u(numax),ltot_u(numax),  !sens wrt vbar,wbar
     & dtot_d(ndmax),ytot_d(ndmax),ltot_d(ndmax),  !sens wrt control
     & dtot_j(njmax),ytot_j(njmax),ltot_j(njmax),  !sens wrt jet control
     & dtot_g(ngmax),ytot_g(ngmax),ltot_g(ngmax),  !sens wrt design
     & ft(3),   ! total Fxyz / rho Vref^2
     & mt(3),   ! total Mxyz / rho Vref^2
     & ft_u(3,numax), !sens wrt vbar,wbar
     & mt_u(3,numax), !sens wrt vbar,wbar
     & ft_d(3,ndmax), !sens wrt control
     & mt_d(3,ndmax), !sens wrt control
     & ft_j(3,njmax), !sens wrt jet control
     & mt_j(3,njmax), !sens wrt jet control
     & ft_g(3,ngmax), !sens wrt design
     & mt_g(3,ngmax), !sens wrt design
     & fti(3),   ! inviscid (vortex) Fxyz / rho Vref^2
     & mti(3),   ! inviscid (vortex) Mxyz / rho Vref^2
     & fti_u(3,numax), !sens wrt vbar,wbar
     & mti_u(3,numax), !sens wrt vbar,wbar
     & fti_d(3,ndmax), !sens wrt control
     & mti_d(3,ndmax), !sens wrt control
     & fti_j(3,njmax), !sens wrt jet control
     & mti_j(3,njmax), !sens wrt jet control
     & fti_g(3,ngmax), !sens wrt design
     & mti_g(3,ngmax), !sens wrt design
     & ftj(3),   ! jet Fxyz / rho Vref^2
     & mtj(3),   ! jet Mxyz / rho Vref^2
     & ftj_u(3,numax), !sens wrt vbar,wbar
     & mtj_u(3,numax), !sens wrt vbar,wbar
     & ftj_d(3,ndmax), !sens wrt control
     & mtj_d(3,ndmax), !sens wrt control
     & ftj_j(3,njmax), !sens wrt jet control
     & mtj_j(3,njmax), !sens wrt jet control
     & ftj_g(3,ngmax), !sens wrt design
     & mtj_g(3,ngmax), !sens wrt design
     & ftp(3),   ! propulsor Fxyz / rho Vref^2
     & mtp(3),   ! propulsor Mxyz / rho Vref^2
     & ftp_u(3,numax), !sens wrt vbar,wbar
     & mtp_u(3,numax), !sens wrt vbar,wbar
     & ftp_d(3,ndmax), !sens wrt control
     & mtp_d(3,ndmax), !sens wrt control
     & ftp_j(3,njmax), !sens wrt jet control
     & mtp_j(3,njmax), !sens wrt jet control
     & ftp_g(3,ngmax), !sens wrt design
     & mtp_g(3,ngmax), !sens wrt design
     & ftv(3),   ! viscous Fxyz / rho Vref^2		
     & mtv(3),   ! viscous Mxyz / rho Vref^2
     & ftv_u(3,numax), !sens wrt vbar,wbar
     & mtv_u(3,numax), !sens wrt vbar,wbar
     & ftv_d(3,ndmax), !sens wrt control
     & mtv_d(3,ndmax), !sens wrt control
     & ftv_j(3,njmax), !sens wrt jet control
     & mtv_j(3,njmax), !sens wrt jet control
     & ftv_g(3,ngmax), !sens wrt design
     & mtv_g(3,ngmax), !sens wrt design
     & ttot,  ! total jet thrust (isolated)
     & dqt,   ! total jet mass flow excess
     & djt,   ! total jet momentum excess
     & det,   ! total jet K.E. excess
     & ttot_j(njmax),
     & dqt_u(numax),
     & djt_u(numax),
     & det_u(numax),
     & dqt_j(njmax),
     & djt_j(njmax),
     & det_j(njmax),
     & chinge(ndmax), 
     & chinge_u(ndmax,numax),
     & chinge_d(ndmax,ndmax),
     & chinge_g(ndmax,ngmax),
     & chinge_j(ndmax,njmax),
     & cf(3),      ! force coefficients Cx,Cy,Cz in geometry axes
     & cf_v(3,3),  ! stability derivatives (Cx,Cy,Cz)_(u,v,w) * Vref
     & cf_w(3,3),  ! stability derivatives (Cx,Cy,Cz)_(p,q,r) * 2Vref/lref
     & cf_d(3,ndmax),  ! control derivatives
     & cf_j(3,njmax),  ! jet control derivatives
     & cf_g(3,ngmax),  ! design derivatives
     & cm(3),      ! moment coefficients Cl,Cm,Cn in geometry axes
     & cm_v(3,3),
     & cm_w(3,3),
     & cm_d(3,ndmax),
     & cm_j(3,njmax),
     & cm_g(3,ngmax),
     & dcl_a0,  ! default additional stability derivatives
     & dcl_u0,
     & dcl_ad0,
     & dcd_a0,
     & dcd_u0,
     & dcd_ad0,
     & dcm_a0,
     & dcm_u0,
     & dcm_ad0

      complex eval, evec, evaldat
      common /case_z/ 
     &  eval(jemax,nrmax),       ! mode eigenvalue
     &  evec(jemax,jemax,nrmax), ! mode eigenvector
     &  evaldat(jemax,nrmax)     ! mode eigenvalue reference data

      common /time_i/
     &  itlev,           ! current time level
     &  ntlev,           ! number of stored time levels
     &  ntsteps          ! default number of time steps to run

      common /time_r/
     &  deltat,                 ! integration time step
     &  tder(3),                ! time-derivative differencing coefficient
     &  tlev(ntmax),            ! time values
     &  tpars(kptot,ntmax),     ! scalar parameter time trace
     &  tparv(3,kpvtot,ntmax),  ! vector parameter time trace
     &  tpard(ndmax,ntmax)      ! control parameter time traces

      common /mass_r/
     &   rho0, gee0,     ! density, gravity   | from .mass file
     &   xyzmass0(3),    ! mass centroid      | from .mass file
     &   rmass0     ,    ! real mass          | from .mass file
     &   riner0(3,3),    ! real inertia       | from .mass file
     &   amass(3,3),     ! apparent mass/rho     | from geometry
     &   ainer(3,3)      ! apparent inertia/rho  | from geometry

      logical lfwake, lfalbe, lfload
      common /surf_l/  
     & lfwake(nfmax),   ! T if surface sheds a wake
     & lfalbe(nfmax),   ! T if surface is to see freestream alpha,beta
     & lfload(nfmax)    ! T if surface contributes to overall loads

      common /surf_i/  
     & nj(nfmax),       ! number of elements along span  in surface
     & nk(nfmax),       ! number of elements along chord in surface
     & jfrst(nfmax),    ! index of first strip in surface
     & jlast(nfmax),    ! index of last  strip in surface
     & idups(nfmax),    ! indicates whether surface is a yduplicated one
     & lscomp(nfmax),   ! logical surface component index
     & izimn(nfmax)     ! z-image flag for surface

      common /surf_r/
     & ssurf(nfmax),     ! surface area
     & cavesurf(nfmax),  ! surface average chord
     & rhon(nfmax),      ! surface fluid density factor
     & fni(3,nfmax),  ! surface inviscid (vortex) Fxyz / rho Vref^2
     & mni(3,nfmax),  ! surface inviscid (vortex) Mxyz / rho Vref^2
     & fni_u(3,nfmax,numax),
     & mni_u(3,nfmax,numax),
     & fni_d(3,nfmax,ndmax),
     & mni_d(3,nfmax,ndmax),
     & fni_g(3,nfmax,ngmax),
     & mni_g(3,nfmax,ngmax),
     & fni_j(3,nfmax,njmax),
     & mni_j(3,nfmax,njmax),
     & fnj(3,nfmax),  ! surface jet Fxyz / rho Vref^2
     & mnj(3,nfmax),  ! surface jet Mxyz / rho Vref^2
     & fnj_u(3,nfmax,numax),
     & mnj_u(3,nfmax,numax),
     & fnj_d(3,nfmax,ndmax),
     & mnj_d(3,nfmax,ndmax),
     & fnj_g(3,nfmax,ngmax),
     & mnj_g(3,nfmax,ngmax),
     & fnj_j(3,nfmax,njmax),
     & mnj_j(3,nfmax,njmax),
     & fnp(3,nfmax),  ! surface propulsor Fxyz / rho Vref^2
     & mnp(3,nfmax),  ! surface propulsor Mxyz / rho Vref^2
     & fnp_u(3,nfmax,numax),
     & mnp_u(3,nfmax,numax),
     & fnp_d(3,nfmax,ndmax),
     & mnp_d(3,nfmax,ndmax),
     & fnp_g(3,nfmax,ngmax),
     & mnp_g(3,nfmax,ngmax),
     & fnp_j(3,nfmax,njmax),
     & mnp_j(3,nfmax,njmax),
     & fnv(3,nfmax),  ! surface viscous Fxyz / rho Vref^2
     & mnv(3,nfmax),  ! surface viscous Mxyz / rho Vref^2
     & fnv_u(3,nfmax,numax),
     & mnv_u(3,nfmax,numax),
     & fnv_d(3,nfmax,ndmax),
     & mnv_d(3,nfmax,ndmax),
     & fnv_g(3,nfmax,ngmax),
     & mnv_g(3,nfmax,ngmax),
     & fnv_j(3,nfmax,njmax),
     & mnv_j(3,nfmax,njmax),
     & dqn(nfmax),               ! surface Dm / rho Vref
     & djn(nfmax),               ! surface DJ / rho Vref^2
     & den(nfmax),               ! surface DE / rho Vref^3
     & dqn_u(nfmax,numax),
     & djn_u(nfmax,numax), 
     & den_u(nfmax,numax),
     & dqn_j(nfmax,njmax),
     & djn_j(nfmax,njmax), 
     & den_j(nfmax,njmax)

      common /strp_i/
     & isurfs(nsmax),  ! index of surface which contains strip
     & ifrstu(nsmax),  ! index of first jet element upstream of strip LE
     & ilastu(nsmax),  ! index of last  jet element upstream of strip LE
     & ifrsts(nsmax),  ! index of first surface element in strip 
     & ilasts(nsmax),  ! index of last  surface element in strip 
     & ifrstw(nsmax),  ! index of first jet element downstream of strip TE
     & ilastw(nsmax),  ! index of last  jet element downstream of strip TE
     & istype(nsmax),  ! strip-type index  (0=solid, 1=jet)
     & izims(nsmax)    ! z-image flag for strip

      logical lstripoff,lviscstrp,lj1sect,lj2sect,ldstrp
      common /strp_l/
     & lstripoff(nsmax),  ! t if strip is "turned off" (outside of fluid)
     & lviscstrp(nsmax),  ! t is strip has viscous drag data
     & lj1sect(nsmax),    ! t if station 1 is a section
     & lj2sect(nsmax),    ! t if station 2 is a section
     & ldstrp(nsmax)      ! t if strip has an actuator disk

      common /strp_r/
     & rhos(nsmax),       ! strip fluid density factor
     & rle(3,nsmax),  chord(nsmax),    ! strip c.p. line LE point, chord
     & rle1(3,nsmax), chord1(nsmax),   ! strip left  end LE point, chord
     & rle2(3,nsmax), chord2(nsmax),   ! strip right end LE point, chord
     & wstrip(nsmax),                  ! strip y-z width
     & tanle(nsmax),  tante(nsmax),    ! strip le,te sweep slopes
     & vhinge(3,nsmax,ndmax),  ! hinge vector for control rot. of normal
     & phinge(3,nsmax,ndmax),  ! point on hingeline for hinge moment calculation
     & reflsd(nsmax,ndmax),    ! sign applied to hinge vec. of refl. strip
     & reflsj(nsmax,ndmax),    ! sign applied to jet control of refl. strip
     & cdcls(6,nsmax),               ! strip viscous polar
     & saxfr,                       ! x/c of spanwise axis for vperp def
     & ess(3,nsmax),                ! spanwise unit vector for vperp def
     & ensy(nsmax), ensz(nsmax),    ! strip normal vector in trefftz-plane
     & xyzrefs(3,nsmax),
     & cdv_strp(nsmax),          ! strip's profile cd
     & ainc(nsmax),              ! strip's incidence twist angle
     & ainc_g(nsmax,ngmax),      ! dainc/dg
     & ccli(nsmax), dwwake(nsmax),
     & fsi(3,nsmax),            ! strip inviscid (vortex) Fxyz / rho Vref^2
     & msi(3,nsmax),            ! strip inviscid (vortex) Mxyz / rho Vref^2
     & fsi_u(3,nsmax,numax), 
     & msi_u(3,nsmax,numax), 
     & fsi_d(3,nsmax,ndmax), 
     & msi_d(3,nsmax,ndmax), 
     & fsi_j(3,nsmax,njmax), 
     & msi_j(3,nsmax,njmax), 
     & fsi_g(3,nsmax,ngmax), 
     & msi_g(3,nsmax,ngmax), 
     & fsj(3,nsmax),            ! strip jet Fxyz / rho Vref^2
     & msj(3,nsmax),            ! strip jet Mxyz / rho Vref^2
     & fsj_u(3,nsmax,numax),
     & msj_u(3,nsmax,numax),
     & fsj_d(3,nsmax,ndmax),
     & msj_d(3,nsmax,ndmax),
     & fsj_j(3,nsmax,njmax),
     & msj_j(3,nsmax,njmax),
     & fsj_g(3,nsmax,ngmax),
     & msj_g(3,nsmax,ngmax),
     & fsp(3,nsmax),            ! strip propulsor Fxyz / rho Vref^2
     & msp(3,nsmax),            ! strip propulsor Mxyz / rho Vref^2
     & fsp_u(3,nsmax,numax), 
     & msp_u(3,nsmax,numax), 
     & fsp_d(3,nsmax,ndmax), 
     & msp_d(3,nsmax,ndmax), 
     & fsp_j(3,nsmax,njmax), 
     & msp_j(3,nsmax,njmax), 
     & fsp_g(3,nsmax,ngmax), 
     & msp_g(3,nsmax,ngmax), 
     & fsv(3,nsmax),            ! strip viscous Fxyz / rho Vref^2
     & msv(3,nsmax),            ! strip viscous Mxyz / rho Vref^2
     & fsv_u(3,nsmax,numax),
     & msv_u(3,nsmax,numax),
     & fsv_d(3,nsmax,ndmax),
     & msv_d(3,nsmax,ndmax),
     & fsv_j(3,nsmax,njmax),
     & msv_j(3,nsmax,njmax),
     & fsv_g(3,nsmax,ngmax),
     & msv_g(3,nsmax,ngmax),
     & dqs(nsmax),        ! strip jet mass excess     = Dm' / rho Vref
     & djs(nsmax),        ! strip jet momentum excess = DJ' / rho Vref^2
     & des(nsmax),        ! strip jet K.E. excess     = DE' / rho Vref^3
     & dqs_u(nsmax,numax),
     & djs_u(nsmax,numax),
     & des_u(nsmax,numax),
     & dqs_j(nsmax,njmax),
     & djs_j(nsmax,njmax),
     & des_j(nsmax,njmax),
     & cncj(nsmax),
     & hdstrp(nsmax),        ! strip propulsor-disk height
     & fhstrp(nsmax),        ! strip jet/disk height-relation factor
     & dj0strp(nsmax),       ! strip jet constant deviation parameter
     & dj1strp(nsmax),       ! strip jet linear   deviation parameter
     & dj3strp(nsmax),       ! strip jet cubic    deviation parameter
     & dxdstrp(nsmax),       ! strip jet x offset from LE
     & dxdstrp1(nsmax),      ! strip jet x offset from LE on left  end
     & dxdstrp2(nsmax),      ! strip jet x offset from LE on right end
     & dndstrp(nsmax),       ! strip jet n offset from LE
     & dndstrp1(nsmax),      ! strip jet n offset from LE on left  end
     & dndstrp2(nsmax),      ! strip jet n offset from LE on right end
     & taxstrp(3,nsmax),     ! thrust axis unit vector, towards + thrust
     & gjstrp(nsmax,njmax)   ! strip dvjet/vbar / jet variable  (gain)

      logical lvalbe
      common /vrtx_l/
     & lvalbe(nvmax)     ! T if c.p. is to see freestream alpha,beta, rotation

      common /vrtx_i/
     & lscompv(nvmax),  ! logical index of surface which contains vortex element
     & isurfv(nvmax),   ! index of surface which contains vortex element
     & jstripv(nvmax),  ! index of strip   which contains vortex element
     & ijetm(nvmax),    ! index of jet element just upstream (0 for surface)
     & izimv(nvmax)     ! z-image flag for vortex element

      common /vrtx_r/
     & rhov(nvmax),     ! h.v. fluid density factor
     & rv1(3,nvmax),    ! h.v. bound-leg left  points
     & rv2(3,nvmax),    ! h.v. bound-leg right points
     & rv(3,nvmax),     ! h.v. bound-leg center points
     & rc(3,nvmax),     ! h.v. control points
     & rs(3,nvmax),     ! h.v. source points
     & dxv(nvmax),               ! chord of element
     & chordv(nvmax),            ! chord of element-containing strip
     & slopev(nvmax),            ! camber slopes at h.v. bound leg
     & slopec(nvmax),            ! camber slopes at c.p.
     & dcontrol(nvmax,ndmax),  ! d(normal-vector angle)/dcontrol
     & enc(3,nvmax),           ! c.p. normal vector
     & env(3,nvmax),           ! h.v. normal vector
     & enc_d(3,nvmax,ndmax),   ! sensitivities
     & enc_g(3,nvmax,ngmax),
     & env_d(3,nvmax,ndmax),
     & env_g(3,nvmax,ngmax),
     & dcpv(nvmax),             ! delta(cp) on vortex element
     & dcpv_u(nvmax,numax),
     & dcpv_d(nvmax,ndmax),
     & dcpv_j(nvmax,njmax),
     & dcpv_g(nvmax,ngmax),
     & dcpj(nvmax),       ! additional delta(cp) from surface jet curvature
     & dcpj_u(nvmax,numax),
     & dcpj_d(nvmax,ndmax),
     & dcpj_j(nvmax,njmax),
     & dcpj_g(nvmax,ngmax),
     & gam(nvmax),             ! circulation of h.v.
     & dgam(nvmax),            ! gam Newton change
     & gam_u(nvmax,numax),
     & gam_d(nvmax,ndmax),
     & gam_g(nvmax,ngmax),
     & gam_j(nvmax,njmax),
     & wc(3,nvmax),          ! body velocity at c.p. locations rc
     & wv(3,nvmax),          ! body velocity at h.v. locations rv
     & wl(3,nlmax),          ! body velocity at h.v. locations rl
     & wc_u(3,nvmax,numax),
     & wv_u(3,nvmax,numax),
     & wl_u(3,nlmax,numax)

      logical lbwake, lbalbe, lbload
      common /body_l/  
     & lbwake(nbmax),   ! T if body sheds a wake
     & lbalbe(nbmax),   ! T if body is to see freestream alpha,beta
     & lbload(nbmax)    ! T if body contributes to overall loads

      common /body_i/  
     & lfrst(nbmax),    ! index of first line node in body
     & llast(nbmax),    ! index of last  line node in body
     & ibcent(nbmax),   ! = 1 if body is on centerline, = 0 otherwise
     & lbcomp(nbmax),   ! logical index of body
     & lbcompl(nlmax),  ! logical index of body which contains source element
     & izimb(nbmax),    ! z-image flag for body
     & iziml(nlmax)     ! z-image flag for body segment

      common /body_r/
     & rhob(nbmax),     ! body fluid density factor
     & cbody(nbmax),    ! body "chord"  (length)
     & abody(nbmax),    ! body planform area
     & sbody(nbmax),    ! body surface area
     & vbody(nbmax),    ! body volume
     & rhol(nlmax),     ! source line fluid density factor
     & rl1(3,nlmax),    ! source line segment upstream endpoints
     & rl2(3,nlmax),    ! source line segment downstream endpoints
     & rad1(nlmax),     ! source line segment upstream radii
     & rad2(nlmax),     ! source line segment downstream radii
     & rl(3,nlmax),     ! source line segment midpoints
     & darl(nlmax),     ! source line segment area changes
     & drbl(nlmax),     ! source line segment radius changes
     & radl(nlmax),     ! source line segment radii
     & src(nlmax),     ! source  strength of source+doublet line element
     & dbl(3,nlmax),   ! doublet strength of source+doublet line element
     & src_u(nlmax,numax),
     & src_d(nlmax,ndmax),
     & src_j(nlmax,njmax),
     & src_g(nlmax,ngmax),
     & dbl_u(3,nlmax,numax),
     & dbl_d(3,nlmax,ndmax),
     & dbl_j(3,nlmax,njmax),
     & dbl_g(3,nlmax,ngmax),
     & dcpb(3,nlmax),
     & fbi(3,nbmax),  ! total body inviscid Fxyz / rho Vref^2			
     & mbi(3,nbmax),  ! total body inviscid Mxyz / rho Vref^2
     & fbi_u(3,nbmax,numax),
     & mbi_u(3,nbmax,numax),
     & fbi_d(3,nbmax,ndmax),
     & mbi_d(3,nbmax,ndmax),
     & fbi_j(3,nbmax,njmax),
     & mbi_j(3,nbmax,njmax),
     & fbi_g(3,nbmax,ngmax),
     & mbi_g(3,nbmax,ngmax),
     & fbv(3,nbmax),  ! total body viscous Fxyz / rho Vref^2			
     & mbv(3,nbmax),  ! total body viscous Mxyz / rho Vref^2
     & fbv_u(3,nbmax,numax),
     & mbv_u(3,nbmax,numax),
     & fbv_d(3,nbmax,ndmax),
     & mbv_d(3,nbmax,ndmax),
     & fbv_j(3,nbmax,njmax),
     & mbv_j(3,nbmax,njmax),
     & fbv_g(3,nbmax,ngmax),
     & mbv_g(3,nbmax,ngmax)

      common /solv_i/ 
     & iapiv(nvmax)           ! pivot indices for pivoting LU solver

      common /solv_r/  
     & amach,                 ! mach number at which aic matrices were computed
     & work(nvmax),           ! work array for factoring aicn matrix
     & resn(nvmax),           ! control-point residual
ccc     & aicn(nvmax,nvmax),     ! d resn(i)/d gam(j) -> now dynamically allocated
     & resn_u(nvmax,numax),   ! d resn(i)/d vbar,wbar
     & resn_d(nvmax,ndmax),   ! d resn(i)/d delcon(l)
     & resn_j(nvmax,njmax),   ! d resn(i)/d deljet(l)
     & resn_g(nvmax,ngmax),   ! d resn(i)/d desvar(l)
     & vc(3,nvmax),           ! vortex velocities at c.p. locations rc
ccc     & vc_gam(3,nvmax,nvmax), ! d vc/d gam    influence matrix  -> now dynamically allocated
     & vc_u(3,nvmax,numax),   ! d vc/d vbar,wbar
     & vc_d(3,nvmax,ndmax),
     & vc_j(3,nvmax,njmax),
     & vc_g(3,nvmax,ngmax),
     & vv(3,nvmax),           ! vortex velocities at h.v. locations rv
ccc     & vv_gam(3,nvmax,nvmax), ! d vv/d gam    influence matrix -> now dynamically allocated
     & vv_u(3,nvmax,numax),   ! d vv/d vbar,wbar
     & vv_d(3,nvmax,ndmax),   ! d vv/d delcon
     & vv_j(3,nvmax,njmax),   ! d vv/d deljet
     & vv_g(3,nvmax,ngmax),   ! d vv/d deldes
     & vl(3,nlmax),           ! vortex velocities at body locations rl
     & vl_gam(3,nlmax,nvmax), ! d vb/d gam    influence matrix
     & vl_u(3,nlmax,numax),   ! d vb/d vbar,wbar
     & vl_d(3,nlmax,ndmax),   ! d vb/d delcon
     & vl_j(3,nlmax,njmax),   ! d vb/d deljet
     & vl_g(3,nlmax,ngmax)    ! d vb/d deldes

c---- plot-selection flags
      logical
     & lwakeplt,   lloadplt,   lwaklplt,  lhingeplt,
     & lboundplt,  lchordplt,  lcambplt,  lcntlplt, 
     & lnrmlplt,   laxesplt,   lvelcplt,  lwsegplt,
     & lrrefplt,   ldiskplt,
     & lclperplt,  ldwashplt,  llabsurf,   ltrforce,
     & label_surf, label_strp, label_vrtx, label_body,
     & lpltnew,    
     & lpltsurf,
     & lpltbody,
     & lsvmov
      common /plot_l/
     & lwakeplt,   lloadplt,   lwaklplt,  lhingeplt,
     & lboundplt,  lchordplt,  lcambplt,  lcntlplt, 
     & lnrmlplt,   laxesplt,   lvelcplt,  lwsegplt,
     & lrrefplt,   ldiskplt,
     & lclperplt,  ldwashplt,  llabsurf,   ltrforce,
     & label_surf, label_strp, label_vrtx, label_body,
     & lpltnew,    
     & lpltsurf(nfmax),
     & lpltbody(nbmax),
     & lsvmov

      common /plot_i/ ntri, imarksurf, naxann(3), ircolor(nrmax)

      common /plot_r/ tri(16,itmax),
     &                axmin(3),axmax(3),axdel(3),axspan(3),
     &                gmin(3),gmax(3),gminp(3),gmaxp(3),
     &                vmin(3),vmax(3),vminp(3),vmaxp(3),
     &                dtmovie,tmovie, slomof,
     &                ephase, eigenf, tmofac
