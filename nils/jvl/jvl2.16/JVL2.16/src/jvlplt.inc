
c---- global variables for JVL plotting routines

      logical lplot, lcurs, lhidden, lcrev
      common /gplt_l/ lplot, lcurs, lhidden, lcrev

      common /gplt_i/
     &   idev, idevh, idevm, ipslu, ncolors

      common /gplt_r/
     &   scrnfrac, plsize, plotar,
     &   xmin,  xmax,  ymin, ymax,
     &   xoff,  yoff,  sf,  csize,
     &   xmarg,ymarg,
     &   xpage,ypage,
     &   xwind,ywind,
     &   pmarg,cpfac,hnfac,enfac,wsfac,
     &   azimob,elevob,tiltob,robinv

      parameter (nhaxis = 12)  ! number of segments forming x,y,z arrowheads
      parameter (nhx = 2*nhaxis+1)
      common /axes_i/
     &    nhead ,
     &    nlinax

      common /axes_r/
     &    arwlen,
     &    rhead,
     &    rheadp,
     &    uaxdir(3,3),
     &    uaxdirp(3,3),
     &    uaxarw(3,2,nhx,3),
     &    uaxarwp(3,2,nhx,3)
