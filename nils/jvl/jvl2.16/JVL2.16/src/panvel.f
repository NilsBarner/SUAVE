
      SUBROUTINE PANVEL(rp1,rp2,rp3,rp4,nc,rc,BETA,
     &                  SPOT,DPOT,SVEL,DVEL)
C
C
C...PURPOSE    TO CALCULATE THE POTENTIAL AND VELOCITY INDUCED BY 
C              SOURCE AND DOUBLET QUADRILATERAL ELEMENTS OF 
C              UNIFORM STRENGTH AT A SET OF FIELD POINTS 
C              (THE DOUBLET AXIS POINTS OPPOSITE TO THE DIRECTION
C              OF THE OUTWARD NORMAL).
C
C...INPUT      rp1,rp2,  =   COORDINATES OF THE CORNER POINTS
C               rp3,rp4      OF THE INFLUENCING ELEMENT
C                              (INPUT IN THE SENSE DEFINED BY THE
C                              OUTWARD NORMAL VIA RIGHT-HAND RULE)
C              nc      =   NUMBER OF FIELD POINTS
C              rc      =   COORDINATES OF SET OF FIELD POINTS
C              BETA        =   SQRT( 1 - MINF**2 )
C
C...OUTPUT     SPOT, DPOT  =   SOURCE AND DOUBLET POTENTIAL INDUCED
C                              BY THE ELEMENT AT EACH FIELD POINT
C              SVEL, DVEL  =   SOURCE AND DOUBLET VELOCITY INDUCED
C                              BY THE ELEMENT AT EACH FIELD POINT
C
C...SUBROUTINES
C   CALLED     CROSS,DOT.
C
C...DISCUSSION THE POTENTIAL AND VELOCITY IN COMPRESSIBLE FLOW ARE 
C              COMPUTED USING AN EQUIVALENT INCOMPRESSIBLE CALCULATION
C              IN WHICH THE CORNER POINTS AND FIELD POINTS ARE SCALED
C              VIA A PRANDTL-GLAUERT TRANSFORMATION.
C              THE LINES JOINING THE MIDPOINTS OF THE OPPOSITE SIDES
C              OF THE ELEMENTS DEFINE A PLANE, AND THEIR POINT OF
C              INTERSECTION DEFINES THE CENTROID OF THE ELEMENT.
C              THE ELEMENT IS MADE PLANAR BY PROJECTION INTO THE
C              PLANE OF THE MIDPOINTS.
C              IF THE FIELD POINT IS SUFFICIENTLY FAR FROM THE
C              ELEMENT CENTROID, THE POTENTIAL AND VELOCITY ARE 
C              COMPUTED AS IF THE ELEMENT WAS A POINT SINGULARITY
C              LOCATED AT ITS CENTROID.
C              THE ELEMENT AND FIELD POINT ARE ROTATED SO THAT
C              THE ELEMENT NORMAL POINTS ALONG THE +Z-AXIS AND
C              TRANSLATED SO THAT THE ELEMENT LIES IN THE X-Y PLANE.
C              THE PLANAR FORMULAS OF HESS AND SMITH ARE THEN USED.
C
C...FIX 1/86  HHY  CORRECTED NUMERICAL BUG (FOR R=0, NEARLY)
C...FIX 4/95  HHY  CORRECTED NUMERICAL BUG FOR COSINE SPACING NEAR WAKES
C                  CALCULATIONS INTERNAL TO POTVEL DONE DOUBLE PRECISION
C...FIX 5/13  MD   Corrected numerical precision problem when distant 
C                  field point is nearly along edge vector. Now switches
C                  asymptotic formula when threshold < EDGTOL2
C
      IMPLICIT REAL (A-H,O-Z)

c---- passed variables
      INTEGER  nc
      REAL rp1(3), rp2(3), rp3(3), rp4(3), rc(3,nc),
     &     BETA, SPOT(nc), DPOT(nc),
     &     SVEL(3,nc), DVEL(3,nc)
C
C---- local arrays
      REAL P0(3), d1(3), d2(3), R0(3), dxd(3), T(3,3),
     &     PG(3,4), XPT(4), YPT(4), D(4), C(4), S(4), 
     &     SV(3), DV(3)
      REAL DIST(4)
C
c---- pi  and  1/(4 pi)
      DATA PI / 3.14159265358979323846264338327950 /
      DATA PI4I / 7.957747154594766788444188168625718E-02 /

      DATA  FACTOR, EDGTOL, EDGTOL2 / 5.0, 1.0E-03, 5.0E-5 /
C
C
C...INPUT CORNER POINTS OF NONPLANAR QUADRILATERAL USING
C...PRANDTL-GLAUERT SUBSONIC COMPRESSIBILITY TRANSFORMATION.
C
      PG(1,1) = rp1(1)
      PG(1,4) = rp2(1)
      PG(1,3) = rp3(1)
      PG(1,2) = rp4(1)
C
      DO K = 2, 3
         PG(K,1) = rp1(K) * BETA
         PG(K,4) = rp2(K) * BETA
         PG(K,3) = rp3(K) * BETA
         PG(K,2) = rp4(K) * BETA
      END DO
C
C
C...COMPUTE MEAN POINT AND MIDPOINT VECTORS.
C...REFERENCE THE CORNER POINTS TO THE MEAN POINT (CENTER).
C
      DO K = 1, 3
         P0(K) = 0.25*( PG(K,1) + PG(K,2) + PG(K,3) + PG(K,4) )
         d1(K) = 0.25*( PG(K,1) + PG(K,2) - PG(K,3) - PG(K,4) )
         d2(K) = 0.25*( PG(K,1) - PG(K,2) - PG(K,3) + PG(K,4) )
         DO I = 1, 4
            PG(K,I) = PG(K,I) - P0(K)
         END DO
      END DO
C
C
C...COMPUTE NORMAL AT MEAN POINT. 
C
cc    CALL DCROSS ( d1,d2,an )
      dxd(1) = d1(2)*d2(3) - d1(3)*d2(2)
      dxd(2) = d1(3)*d2(1) - d1(1)*d2(3)
      dxd(3) = d1(1)*d2(2) - d1(2)*d2(1)

cc    dxdmod = SQRT( DPDOT(dxd,dxd) )
      dxdmod = SQRT( dxd(1)**2 + dxd(2)**2 + dxd(3)**2 )
C
      DO ic=1,nc
        SVEL(1,ic) = 0.
        SVEL(2,ic) = 0.
        SVEL(3,ic) = 0.
      END DO
C
      if(dxdmod .eq. 0.0) then
c----- degenerate panel
       return
      endif
C
C...COMPUTE ELEMENT AREA, REFERENCE LENGTH AND FAR-rc DISTANCE.
C
      area = 4.0 * dxdmod
cc    d1mod = sqrt ( dpdot(d1,d1) )
cc    d2mod = sqrt ( dpdot(d2,d2) )
      d1mod = sqrt( d1(1)**2 + d1(2)**2 + d1(3)**2 )
      d2mod = sqrt( d2(1)**2 + d2(2)**2 + d2(3)**2 )
      reflen =  2.0 * max (d1mod,d2mod)
      ffdist = factor * reflen
c
c
c---- d1,d2,dxd = unit vectors of panel axis system
cc    CALL DCROSS ( dxd,d1,d2 )
      d2(1) = dxd(2)*d1(3) - dxd(3)*d1(2)
      d2(2) = dxd(3)*d1(1) - dxd(1)*d1(3)
      d2(3) = dxd(1)*d1(2) - dxd(2)*d1(1)

cc    d2mod = sqrt ( dpdot(d2,d2) )
      d2mod = sqrt ( d2(1)**2 + d2(2)**2 + d2(3)**2 )
c
      do k = 1, 3
         dxd(k) = dxd(k) / dxdmod
         d1(k)  = d1(k)  / d1mod
         d2(k)  = d2(k)  / d2mod
      end do
C
C
C...COMPUTE area JACOBIAN OF COMPRESSIBLE COORDINATE SYSTEM.
c
      ajcbn = beta / sqrt( (dxd(1)/beta)**2 + dxd(2)**2 + dxd(3)**2 )
      pi4ij = pi4i / ajcbn
c
c
c---- transformation matrix from global to panel coordinates 
      do k = 1, 3
         t(k,1) = d1(k)
         t(k,2) = d2(k)
         t(k,3) = dxd(k)
      enddo
C
C---- panel corners in panel axes
      do i = 1, 4
         xpt(i) = t(1,1)*pg(1,i) + t(2,1)*pg(2,i) + t(3,1)*pg(3,i)
         ypt(i) = t(1,2)*pg(1,i) + t(2,2)*pg(2,i) + t(3,2)*pg(3,i)
      end do
c
c
c
c---- parameters d,c,s in formulas of Hess and Smith
c-    which are independent of field point.
      do i = 1, 4
        ip = mod(i,4) + 1
        dx = xpt(ip) - xpt(i)
        dy = ypt(ip) - ypt(i)
c
        d(i) = sqrt ( dx*dx + dy*dy )
        if( d(i) .ne. 0.0 )  then
         c(i) = dx / d(i)
         s(i) = dy / d(i)
        endif
      enddo
c
c***********************************************************************
c***********************************************************************
c
c...for each rc POINT:
C
      DO 400 ic = 1,nc
C
      DIST(1) = SQRT( (rp1(1)-rc(1,ic))**2
     &              + (rp1(2)-rc(2,ic))**2
     &              + (rp1(3)-rc(3,ic))**2 )
      DIST(4) = SQRT( (rp2(1)-rc(1,ic))**2
     &              + (rp2(2)-rc(2,ic))**2
     &              + (rp2(3)-rc(3,ic))**2 )
      DIST(3) = SQRT( (rp3(1)-rc(1,ic))**2
     &              + (rp3(2)-rc(2,ic))**2
     &              + (rp3(3)-rc(3,ic))**2 )
      DIST(2) = SQRT( (rp4(1)-rc(1,ic))**2
     &              + (rp4(2)-rc(2,ic))**2
     &              + (rp4(3)-rc(3,ic))**2 )
C
C...INPUT rc POINT USING PRANDTL-GLAUERT
C...SUBSONIC COMPRESSIBILITY TRANSFORMATION.
C
      R0(1) = rc(1,ic)         - P0(1)
      R0(2) = rc(2,ic) * BETA  - P0(2)
      R0(3) = rc(3,ic) * BETA  - P0(3)
C
C
C...COMPUTE DISTANCE FROM CENTROID TO rc POINT.
C
C...IF THE rc POINT IS SUFFICIENTLY DISTANT, USE THE POINT
C...SOURCE AND DOUBLET APPROXIMATION FOR FAR rc.
C
cc    R0MOD = SQRT ( DPDOT(R0,R0) )
      R0MOD = SQRT( R0(1)**2 + R0(2)**2 + R0(3)**2 )
C
      IF ( R0MOD .LE. FFDIST )   GO TO 220
C
cc       RNMODR = DPDOT(dxd,R0) / R0MOD**2
         RNMODR = (dxd(1)*R0(1)+dxd(2)*R0(2)+dxd(3)*R0(3)) / R0MOD**2
C
C...FAR-rc POTENTIALS.
C
         IF (ISING .EQ. 1)  GO TO 200
C
            SPOT(ic) =        - area * PI4IJ / R0MOD
            DPOT(ic) = RNMODR * area * PI4I  / R0MOD
C
C
C...FAR-rc VELOCITIES.
C
  200    IF (ISING .EQ. 0)  GO TO 400
C
            RMOD3I = 1. / R0MOD**3
C
            DO K = 1, 3
               SV(K) =                     R0(K)  * area * RMOD3I
               DV(K) = (dxd(K) - 3.*RNMODR*R0(K)) * area * RMOD3I
            END DO
C
            SVEL(1,ic) = SV(1)        *  PI4IJ
            SVEL(2,ic) = SV(2) * BETA *  PI4IJ
            SVEL(3,ic) = SV(3) * BETA *  PI4IJ
C
            DVEL(1,ic) = DV(1)        *  PI4I
            DVEL(2,ic) = DV(2) * BETA *  PI4I
            DVEL(3,ic) = DV(3) * BETA *  PI4I
            GO TO 400
C
C
C
C...COMPUTE TRANSFORMED COORDINATES OF rc POINT IN LOCAL SYSTEM.
C
  220 X = T(1,1)*R0(1) + T(2,1)*R0(2) + T(3,1)*R0(3)
      Y = T(1,2)*R0(1) + T(2,2)*R0(2) + T(3,2)*R0(3)
      Z = T(1,3)*R0(1) + T(2,3)*R0(2) + T(3,3)*R0(3)
C
C
C...APPLY FORMULAS OF HESS AND SMITH.
C
      ZZ   = Z*Z
      ZMAG = ABS(Z)
C
      SGNZ = 0.
      IF (Z .NE. 0.)  SGNZ = SIGN (1.0, Z)
C
      SUMJ   = 0.
      SUMT   = 0.
      SUMRQ  = 0.
C
      SVX = 0.
      SVY = 0.
      DVX = 0.
      DVY = 0.
      DVZ = 0.
C
C
C
C...PROCESS INFLUENCES OF EACH SIDE.
C
      DO 300 I = 1, 4
C
         IF ( D(I) .EQ. 0. )  GO TO 300
C
         IP1 = MOD(I,4) + 1
C
         DX1 = XPT(I)   - X
         DY1 = YPT(I)   - Y
         DX2 = XPT(IP1) - X
         DY2 = YPT(IP1) - Y
C
C
C...COMPUTE R.
C
         R  = -DX1 * S(I)  +  DY1 * C(I)
C
C...COMPUTE DD1,DD2.
C
         DD1 = DX1**2 + DY1**2
         DD2 = DX2**2 + DY2**2
C
C...ENSURE R IS 0. FOR POINT ON A CORNER
C
         IF ((DD1 .EQ. 0.) .OR.
     &       (DD2 .EQ. 0.) )   R = 0.
         RR = R*R
         RRZZ = RR + ZZ
C
C...COMPUTE RR1, RR2, R1, R2, R1R2
C
         RR1 = DD1 + ZZ
         RR2 = DD2 + ZZ
         R1 = SQRT(RR1)
         R2 = SQRT(RR2)
C
         R1R2 = R1*R2
C
C
C...IF THE rc POINT IS ON A CORNER POINT, DO NOT CALCULATE 
C...INFLUENCES FOR THIS SIDE.
C
         IF (R1R2 .EQ. 0.)  GO TO 300
C
C
C...COMPUTE S1,S2.
C
         S1 = DX1 * C(I)  +  DY1 * S(I)
         S2 = DX2 * C(I)  +  DY2 * S(I)
C
         S1S2 = S1*S2
         RS12 = R1*S2 - R2*S1
C
C------- numerical error "bug" fix,  MD  27 May 13

         SSUM = ABS(S1) + ABS(S2)
         RRATIO = MIN(SQRT(RRZZ),DIST(I),DIST(IP1)) / SSUM
         IF(S1S2 .GT. 0.0) THEN
          IF(RRATIO .LT. EDGTOL) THEN
C--------- use asymptotic formula when distant field point is nearly along edge vector 
           RS12 = 0.5*RRZZ/(R1*R2)
          ENDIF

         ELSE
          IF(RRATIO .LT. EDGTOL2) THEN
           R = 0.
           Z = 0.
           RR = 0.
           ZZ = 0.
           RRZZ = 0.
          ENDIF
         ENDIF
C

C
C...COMPUTE Q.
C...CHECK FOR rc POINTS VERY CLOSE TO THE PANEL EDGE.
C...USE A NEAR-rc EXPANSION FOR THESE CASES.
C...THE SINGULARITY IS REMOVED (Q=0) FOR rc POINTS ON EDGE.
C
C
         Q1 = R1 + R2 + D(I)
         Q2 = R1 + R2 - D(I)
         DRATIO = ABS(Q2) / Q1 
C
         Q = 0.
         IF (DRATIO .LE. EDGTOL)  GO TO 250
            Q = LOG ( 1.0 / DRATIO )
            GO TO 260
C
  250    IF (RRZZ .NE. 0.)  Q = LOG ( 4.0*R1R2 / RRZZ )
C
C
C
  260    IF (ISING .EQ. 0)  GO TO 290
C
C...COMPUTE SOURCE VELOCITY TERMS FROM PANEL EDGE.
C
            SVX = SVX - S(I) * Q
            SVY = SVY + C(I) * Q
C
C
C...COMPUTE EDGE TERMS FROM LINE VORTEX ON DOUBLET PANEL EDGE.
C...THE SINGULARITY IS REMOVED FOR rc POINTS ON EDGE.
C
            EDGLV = 0.
            IF (RRZZ .NE. 0.0)  EDGLV = -RS12 / (R1R2*RRZZ)
C
            DVX = DVX + S(I) * (Z * EDGLV)
            DVY = DVY - C(I) * (Z * EDGLV)
            DVZ = DVZ        -  R * EDGLV

C
C
C...THE REMAINDER OF THESE TERMS WILL BE SET TO ZERO FOR A rc POINT 
C...IN LINE WITH THE EDGE.
C
  290    IF ( R .EQ. 0.0 )  GO TO 300
C
         SGNR = SIGN (1.0, R)
C
C
C...COMPUTE R*Q TERMS.
C
         SUMRQ = SUMRQ + R*Q
C
C
C...COMPUTE J.
C
         SINA = R*ZMAG*RS12
         COSA = R1R2*RR + ZZ*S1S2 
         EJ = ATAN2 ( SINA, COSA )
C
         SUMJ = SUMJ + EJ
C
C
C...COMPUTE THETA.
C
         ARG = ( RR + S1S2 ) / SQRT ( DD1*DD2 )
         IF ( ARG .GT.  1.0 )   ARG =  1.0
         IF ( ARG .LT. -1.0 )   ARG = -1.0
C
         THETA = SGNR * ACOS (ARG)
         SUMT  = SUMT + THETA
C
C
  300 CONTINUE
C
C
C...SLIGHTLY NON-ZERO VALUES OF THETA CAUSED BY NUMERICAL ERROR
C...MAY SERIOUSLY AFFECT PRECISION IF J IS ALSO SMALL.
C
      IF ( ABS(SUMT) .LT. 0.01 )  SUMT = 0.0
C
      SUMJT = SUMJ - SUMT
C
C
      IF (ISING .EQ. 1)  GO TO 320
C
C...COMPUTE SOURCE AND DOUBLET POTENTIALS.
C
         SP =  SUMRQ  +  ZMAG * SUMJT
         DP =            SGNZ * SUMJT
C
         SPOT(ic) = -SP * PI4IJ
         DPOT(ic) = -DP * PI4I
C
C
C
  320 IF (ISING .EQ. 0)  GO TO 400
C
C...COMPUTE SOURCE AND DOUBLET VELOCITIES.
C...TRANSFORM VELOCITIES BACK INTO GLOBAL REFERENCE FRAME.
C
         SVZ = - SGNZ * SUMJT
C
         DO K = 1, 3
            SV(K) = T(K,1)*SVX + T(K,2)*SVY + T(K,3)*SVZ
            DV(K) = T(K,1)*DVX + T(K,2)*DVY + T(K,3)*DVZ
         END DO
C
         SVEL(1,ic) =   SV(1)        *  PI4IJ
         SVEL(2,ic) =   SV(2) * BETA *  PI4IJ
         SVEL(3,ic) =   SV(3) * BETA *  PI4IJ
         DVEL(1,ic) =  -DV(1)        *  PI4I
         DVEL(2,ic) =  -DV(2) * BETA *  PI4I 
         DVEL(3,ic) =  -DV(3) * BETA *  PI4I
C
  400 CONTINUE
C
      RETURN
      END ! POTVEL



      SUBROUTINE DCROSS (U,V,W)
C
C...PURPOSE    TO CALCULATE THE VECTOR PRODUCT OF TWO VECTORS.
C
C...INPUT      U,V
C
C...OUTPUT     W = U X V
C
C...SUBROUTINES
C   CALLED     NONE.
C
      REAL  U(3), V(3), W(3)
C
C
      W(1) =  U(2)*V(3) - U(3)*V(2)
      W(2) = -U(1)*V(3) + U(3)*V(1)
      W(3) =  U(1)*V(2) - U(2)*V(1)
C
      RETURN
      END


      REAL FUNCTION DPDOT (U,V)
C
C...PURPOSE    TO CALCULATE THE SCALAR PRODUCT OF TWO VECTORS.
C
C...INPUT      U,V
C
C...OUTPUT     DOT = U O V
C
C...SUBROUTINES
C   CALLED     NONE.
C
      REAL  U(3), V(3)
C
      DPDOT = U(1)*V(1) + U(2)*V(2) + U(3)*V(3)
C
      RETURN
      END

