
      SUBROUTINE LUDCMP(NSIZ,N,A)
C     *******************************************************
C     *   Factors a full NxN matrix into an LU form.        *
C     *   Uses LINPACK routines for linear algebra          *
C     *******************************************************
      REAL A(NSIZ,NSIZ)

      PARAMETER (NVMAX = 4000)
      INTEGER INDX(NVMAX)
      common / lap_i / INDX

      IF(N .GT. NVMAX) STOP 'LUDCMP: Local work array overflow'

      CALL SGETRF(N,NSIZ,A,NSIZ,INDX,INFO)

      RETURN
      END ! LUDCMP


      SUBROUTINE BAKSUB(NSIZ,N,A,B)
      DIMENSION A(NSIZ,NSIZ), B(NSIZ)
C     *******************************************************
C     *   BAKSUB does back-substitution with RHS using      *
C     *   stored LU decomposition.                          *
C     *   Uses LINPACK routines for linear algebra          *
C     *******************************************************

      PARAMETER (nvmax = 4000)
      INTEGER INDX(NVMAX)
      common / lap_i / INDX

      MRHS = 1
      CALL SGETRS('N',N,MRHS,A,NSIZ,INDX,B,NSIZ,INFO)

      RETURN
      END ! BAKSUB



      SUBROUTINE LUDCMPI(NSIZ,N,A,INDX,WORK)
C     *******************************************************
C     *   Factors a full NxN matrix into an LU form.        *
C     *   Uses LINPACK routines for linear algebra          *
C     *******************************************************
C
      REAL A(NSIZ,NSIZ), WORK(NSIZ)
      INTEGER INDX(NSIZ)
C
      CALL SGETRF(N,N,A,NSIZ,INDX,INFO)
C
      RETURN
      END ! LUDCMPI


      SUBROUTINE BAKSUBI(NSIZ,N,A,INDX,B)
      DIMENSION A(NSIZ,NSIZ), B(NSIZ), INDX(NSIZ)
C     *******************************************************
C     *   BAKSUB does back-substitution with RHS using      *
C     *   stored LU decomposition.                          *
C     *   Uses LINPACK routines for linear algebra          *
C     *******************************************************
C
      MRHS = 1
      CALL SGETRS('N',N,MRHS,A,NSIZ,INDX,B,NSIZ,INFO)
C
      RETURN
      END ! BAKSUBI



