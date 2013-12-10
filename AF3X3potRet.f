c***********************************************************************
      SUBROUTINE AF3X3potRet(RDIST,DELTAE,C3val,C3valP,C6val,
     1  C8val,De,ULR,DEIGM1,DEIGM2,DEIGM3,DEIGM5,DEIGR,DEIGDe)
c***********************************************************************
c** Subroutine to generate the lowest eigenvalue of the 3x3 long-range
c  Li2 interaction matrix of Eq.(25) of J.Mol.Spectr. 268, 199 (2011), 
c  and its derivatives w.r.t. the potential function parameters.
c==> Input: distance  r= RDIST, DELTAE, C3, C6eff, C8 & De.
c==> Output:  ULR= -lambda(min) and iuts partial derivatives DEIGxx
c** Version from Nik in June 2011 for Real root case; cosmetically 
c    modified by RJL in Aug 2011
c-----------------------------------------------------------------------
      REAL*8  A(3,3),DM1(3,3),DM2(3,3),DM3(3,3),DM5(3,3),DR(3,3),
     1              DDe(3,3),Q(3,3)
      REAL*8  DEIGM1(1,1),DEIGM2(1,1),DEIGM3(1,1),DEIGM5(1,1),
     1        DEIGR(1,1),DEIGDe(1,1), EIGVEC(3,1), RESID(3,1), W(3) 
      REAL*8  RDIST,RDIST2,RDIST3,DELTAE,C3val,C3valP,C6val,C8val,De,
     1   ULR,RET,RETSig,RETPi,Modulus,M1,M2,M3,M5,Z
      INTEGER          I,J,L,K
      M1= C3val
      M2= C3valP
      M3= C6val
      M5= C8val
      RET= 9.36423830d-4*RDIST
      RETSig= DCOS(RET) + (RET)*DSIN(RET)
      RETPi= RETSig - RET**2 *DCOS(RET)
      RDIST2= RDIST**2
      RDIST3= RDIST*RDIST2
c      WRITE(25,*) 'Variables = "r", "U(r)","U(r)-U(r)^2/(4De)" ' 
c      WRITE(25,*) 'zone T = "U(r)"'
c  Initialize interaction matrix to 0.d0
      DO  I= 1,3
          A(I,I)=0.0D0
          ENDDO
ccccc Prepare interation matrix  A 
      A(1,1)= -(M1*RETSig+ M3/(RDIST3)+M5/(RDIST3*RDIST2))/(3.d0*RDIST3)
      A(1,2)= -(DSQRT(2.D0))*A(1,1)
      A(2,1)= A(1,2)
      A(1,3)= 2*M2*RETPi/(DSQRT(6.D0)*RDIST3)
      A(3,1)= A(1,3)
      A(2,2)= 2*A(1,1) + DELTAE
      A(2,3)= A(1,3)/DSQRT(2.d0)
      A(3,2)= A(2,3)
      A(3,3)= DELTAE
cccccc Prepare radial derivative of interaction matrix (? is it needed ?)
      DR(1,1)= (3.d0*M1*RETSig + 6.d0*M3/RDIST3 
     1                  + 8.D0*M5/(RDIST3*RDIST2))/(3.d0*RDIST3*RDIST)
      DR(1,2)= -DSQRT(2.d0)*DR(1,1)
      DR(2,1)= DR(1,2)
      DR(2,2)= 2.d0*DR(1,1)
      DR(1,3)= -3.d0*A(1,3)/RDIST
      DR(3,1)= DR(1,3)
      DR(2,3)= -3.d0*A(2,3)/RDIST
      DR(3,2)= DR(2,3)
      DR(3,3)= 0.d0 
cccccc Partial derivative of interaction matric  H  w.r.t.  C3
      DM1(1,1)= -(RETSig + M1/(2.d0*De*RDIST3))/(3.d0*RDIST3)
      DM1(1,2)= -DSQRT(2.d0)*DM1(1,1)
      DM1(2,1)= DM1(1,2)
      DM1(2,2)= 2.d0*DM1(1,1)
      DM1(1,3)= 0.d0
      DM1(3,1)= 0.d0
      DM1(2,3)= 0.d0
      DM1(3,2)= 0.d0
      DM1(3,3)= 0.d0
cccccc Partial derivative of interaction matric  H  w.r.t.  C3
      DM2(1,1)= 0.d0
      DM2(1,2)= 0.d0
      DM2(2,1)= 0.d0
      DM2(2,2)= 0.d0
      DM2(1,3)= RETPi/(DSQRT(6.d0)*RDIST3) 
      DM2(3,1)= DM2(1,3)
      DM2(2,3)= DM2(1,3)/DSQRT(2.d0)
      DM2(3,2)= DM2(2,3)
      DM2(3,3)= 0.d0
cccccc Partial derivative of interaction matric  H  w.r.t.  C6
      DM3(1,1)= -1.d0/(3.d0*RDIST3**2)
      DM3(1,2)= -SQRT(2.d0)*DM3(1,1)
      DM3(1,3)= 0.D0
      DM3(2,1)= DM3(1,2)
      DM3(2,2)= 2.d0*DM3(1,1)
      DM3(2,3)= 0.D0
      DM3(3,1)= DM3(1,3)
      DM3(3,2)= DM3(2,3)
      DM3(3,3)= 0.D0
cccccc Partial derivative of interaction matric  H  w.r.t.  C8
      DM5(1,1)= DM3(1,1)/(RDIST2)
      DM5(1,2)= DM3(1,2)/(RDIST2)
      DM5(1,3)= 0.D0
      DM5(2,1)= DM3(1,2)
      DM5(2,2)= DM3(2,2)/(RDIST2)
      DM5(2,3)= 0.D0
      DM5(3,1)= DM5(1,3)
      DM5(3,2)= DM5(2,3)
      DM5(3,3)= 0.D0
cccccc Partial derivative of interaction matric  H  w.r.t.  De
      DDe(1,1)= M1**2/(12.D0*(RDIST3*De)**2)
      DDe(1,2)= -SQRT(2.D0)*DDe(1,1)
      DDe(1,3)= 0.D0
      DDe(2,1)= DDe(1,2)
      DDe(2,2)= 2.D0*DDe(1,1)
      DDe(2,3)= 0.d0
      DDe(3,1)= DDe(1,3)
      DDe(3,2)= DDe(2,3)
      DDe(3,3)= 0.D0
cccccc Call subroutine to prepare and invert interaction matrix  H
      CALL DSYEVJ3(A,Q,W)
      L=1
ccc Nor - identify the lowest eigenvalue of  H  and label it  L
      DO J=2,3
          IF (W(J) .LT. W(L)) THEN
              L=J
              ENDIF
          ENDDO  
      ULR= -W(L)
      DO I=1,3      
          EIGVEC(I,1) = Q(I,L)
          ENDDO  
      DEIGM1= -MATMUL(TRANSPOSE(EIGVEC),MATMUL(DM1,EIGVEC))
      DEIGM2= -MATMUL(TRANSPOSE(EIGVEC),MATMUL(DM2,EIGVEC))
      DEIGM3= -MATMUL(TRANSPOSE(EIGVEC),MATMUL(DM3,EIGVEC))
      DEIGM5= -MATMUL(TRANSPOSE(EIGVEC),MATMUL(DM5,EIGVEC))           
      DEIGR = -MATMUL(TRANSPOSE(EIGVEC),MATMUL(DR,EIGVEC))
      DEIGDe= -MATMUL(TRANSPOSE(EIGVEC),MATMUL(DDe,EIGVEC))
c     WRITE(25,600) RDIST ,ULR 
c 600 FORMAT(2D16.7)
c     WRITE(26,601) RDIST , DEIGM1, DEIGR ,DEIGDe
c 601 FORMAT(4D16.7)  
      RETURN
      CONTAINS
c=======================================================================
      SUBROUTINE DSYEVJ3(A, Q, W)
c ----------------------------------------------------------------------------
c** Subroutine to setup and diagonalize the matrix  A  and return 
c   eigenvalues W and eigenvector matric  Q
      INTEGER N, I, X, Y, R
      PARAMETER (N=3)
      REAL*8 A(3,3), Q(3,3), W(3)
      REAL*8 SD, SO, S, C, T, G, H, Z, THETA, THRESH
c     Initialize Q to the identitity matrix
c     --- This loop can be omitted if only the eigenvalues are desired ---
      DO  X = 1, N
          Q(X,X) = 1.0D0
          DO  Y = 1, X-1
              Q(X, Y) = 0.0D0
              Q(Y, X) = 0.0D0
              ENDDO
          ENDDO
c Initialize W to diag(A)
      DO  X = 1, N
          W(X) = A(X, X)
          ENDDO
c Calculate SQR(tr(A))  
      SD = 0.0D0
      DO  X = 1, N
          SD = SD + ABS(W(X))
          ENDDO
      SD = SD**2
c Main iteration loop
      DO 40 I = 1, 50
c Test for convergence
          SO = 0.0D0
          DO  X = 1, N
              DO  Y = X+1, N
                  SO = SO + ABS(A(X, Y))
                  ENDDO
              ENDDO
          IF(SO .EQ. 0.0D0)  RETURN
          IF(I .LT. 4) THEN
              THRESH = 0.2D0 * SO / N**2
            ELSE
              THRESH = 0.0D0
            END IF
c Do sweep
          DO 60 X = 1, N
              DO 61 Y = X+1, N
                  G = 100.0D0 * ( ABS(A(X, Y)) )
                  IF ( I .GT. 4 .AND. ABS(W(X)) + G .EQ. ABS(W(X))
     $                          .AND. ABS(W(Y)) + G .EQ. ABS(W(Y))) THEN
                      A(X, Y) = 0.0D0
                    ELSE IF (ABS(A(X, Y)) .GT. THRESH) THEN
c Calculate Jacobi transformation
                      H = W(Y) - W(X)
                      IF ( ABS(H) + G .EQ. ABS(H) ) THEN
                          T = A(X, Y) / H
                        ELSE
                          THETA = 0.5D0 * H / A(X, Y)
                          IF (THETA .LT. 0.0D0) THEN
                              T= -1.0D0/(SQRT(1.0D0 + THETA**2)-THETA)
                            ELSE
                              T= 1.0D0/(SQRT(1.0D0 + THETA**2) + THETA)
                            END IF
                        END IF
                      C = 1.0D0 / SQRT( 1.0D0 + T**2 )
                      S = T * C
                      Z = T * A(X, Y)
c Apply Jacobi transformation
                      A(X, Y) = 0.0D0
                      W(X)    = W(X) - Z
                      W(Y)    = W(Y) + Z
                      DO  R = 1, X-1
                          T       = A(R, X)
                          A(R, X) = C * T - S * A(R, Y)
                          A(R, Y) = S * T + C * A(R, Y)
                          ENDDO
                      DO  R = X+1, Y-1
                          T       = A(X, R)
                          A(X, R) = C * T - S * A(R, Y)
                          A(R, Y) = S * T + C * A(R, Y)
                          ENDDO
                      DO  R = Y+1, N
                          T       = A(X, R)
                          A(X, R) = C * T - S * A(Y, R)
                          A(Y, R) = S * T + C * A(Y, R)
                          ENDDO
c Update eigenvectors
c --- This loop can be omitted if only the eigenvalues are desired ---
                      DO  R = 1, N
                          T       = Q(R, X)
                          Q(R, X) = C * T - S * Q(R, Y)
                          Q(R, Y) = S * T + C * Q(R, Y)
                          ENDDO
                    END IF
   61             CONTINUE
   60         CONTINUE
   40     CONTINUE
c      PRINT *, "DSYEVJ3: No convergence."
      END SUBROUTINE DSYEVJ3
      END SUBROUTINE AF3X3potRet
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
