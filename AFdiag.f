c***********************************************************************
      SUBROUTINE AFdiag(RDIST,VLIM,NCMM,NCMMax,MMLR,Cm,rhoAB,IVSR,
     1                                       IDSTT,ULR,dULRdCm,dULRdR)
c***********************************************************************
c**   Aubert-Frecon Potential Model for u_{LR}(r)
c***********************************************************************
c** Subroutine to generate, at the onee distance RDIST, an eigenvalue 
c  of the 2x2 or 3x3 long-range interaction matrix described by Eqs.1
c and 10, resp., of J.Mol.Spec.188, 182 (1998) (Aubert-Frecon et al)
c** and its derivatives w.r.t. the potential function parameters.
c***********************************************************************
c==> Input: distance r = RDIST,VLIM,NCMM,m=MMLR & Cm.
c==> Output: ULR= -lambda(min or max) and its partial derivatives DEIGxx
c-----------------------------------------------------------------------
c** Original Version from Nike Dattani in June 2011 for 3x3 case
c** Generalized to incorporate 2x2 case and remove retardation terms
c                      by Kai Slaughter: July 2014
c-----------------------------------------------------------------------
      INTEGER NCMMax
c-----------------------------------------------------------------------
      REAL*8 RDIST,VLIM,Cm(NCMMax),ULR,dULRdCm(NCMMax),dULRdR,R2,R3,R5,
     1       R6,R8,R9,T1,T0,T2,T0P,T0P23,Dm(NCMMax),Dmp(NCMMax),
     2       Dmpp(NCMMax),rhoAB,A(3,3),DR(3,3),Q(3,3),DMx(NCMMax,3,3),
     3       DMtemp(3,3),DEIGMx(NCMMax,1,1),DEIGMtemp(1,1),DEIGR(1,1),
     4       EIGVEC(3,1),RESID(3,1),W(3),RPOW(NCMMax),DELTAE,Modulus,Z
      INTEGER H,I,J,K,L,M,X,NCMM,MMLR(NCMMax),IVSR,IDSTT,MMtemp
c-----------------------------------------------------------------------
      DELTAE=Cm(1)
c-----------------------------------------------------------------------
      IF(rhoAB.GT.0.d0) THEN
          CALL dampF(RDIST,rhoAB,NCMM,MMLR,IVSR,IDSTT,Dm,Dmp,Dmpp)
      ELSE
          DO  m=2,NCMM
              Dm(m)= 1.d0
              Dmp(m)=0.d0
              Dmpp(m)=0.d0
              ENDDO
          ENDIF
c-----------------------------------------------------------------------
      IF(MMLR(1).GE.-1) THEN
c***********************************************************************
c************* Aubert Frecon 2x2 case   NCMM= 7  and  ...
c***              Cm(1) = DELTAE
c***              Cm(2) = C3Sig
c***              Cm(3) = C3Pi
c***              Cm(4) = C6Sig
c***              Cm(5) = C6Pi
c***              Cm(6) = C8Sig
c***              Cm(7) = C8Pi
c***********************************************************************
          R2= 1.d0/RDIST**2
          R3= R2/RDIST
          R5= R2*R3
          R6= R3*R3
          R8= R6*R2
c-----------------------------------------------------------------------
          T1= R3*(Dm(2)*(Cm(2)-Cm(3)) + R3*Dm(4)*(Cm(4)-Cm(5)) + 
     1        R5*Dm(6)*(Cm(6)-Cm(7)))/3.d0
          T0= DSQRT((T1 - Cm(1))**2 + 8.d0*T1**2)
          ULR= 0.5d0*(-Cm(1) + R3*(Dm(2)*(Cm(2)+Cm(3)) + 
     1         R3*Dm(4)*(Cm(4)+Cm(5)) + R5*Dm(6)*(Cm(6)+Cm(7))) + T0)
c-----------------------------------------------------------------------
c...  adjustment for the b-state
          IF(MMLR(2).EQ.-1) ULR=ULR-T0
c...  now get derivatives
          T0P= 0.5d0*(9.d0*T1 - Cm(1))/T0
          T0P23= 0.5d0 + T0P/3.d0
c...  another adjustment for the b-state
          IF(MMLR(2).EQ.-1) T0P23=T0P23-2.d0*T0P/3.d0
          dULRdCm(1)= 0.d0
          dULRdCm(2)= R3*(T0P23)
          dULRdCm(3)= R3*(1.d0-T0P23)
          dULRdCm(4)= R6*(T0P23)
          dULRdCm(5)= R6*(1.d0 - T0P23)
          dULRdCm(6)= R8*T0P23
          dULRdCm(7)= R8*(1.d0-T0P23)
          T2        =-T0P*R3*((Dm(2)*(Cm(2)-Cm(3))+R3*(Dm(4)*2.d0*(Cm(4)
     1                -Cm(5))+R2*Dm(6)*8.d0/3.d0*(Cm(6)-Cm(7))))/RDIST
     2                +(Dmp(2)*(Cm(2)-Cm(3))+R3*Dmp(4)*(Cm(4)-Cm(5))+
     3                R2*R3*Dmp(6)*(Cm(6)-Cm(7)))/3.d0)
          dULRdR    = -R3*((1.5d0*Dm(2)*(Cm(2)+Cm(3)) + R3*(Dm(4)*3.d0*
     1                (Cm(4)+Cm(5))+4.d0*Dm(6)*R2*(Cm(6)+Cm(7))))/RDIST
     2                + 0.5d0*(Dmp(2)*(Cm(2)+Cm(3)) + Dmp(4)*R3*(Cm(4)+
     3                Cm(5)) + Dmp(6)*R3*R2*(Cm(6)+Cm(7)))) + T2
c... and a final adjustment for the b-state
          IF(MMLR(2).EQ.-1) dULRdR=dULRdR-2.d0*T2
c-----------------------------------------------------------------------
      ELSE

c***********************************************************************
c********* Aubert Frecon 3x3 case   NCMM= 10  and ...
c*********        Cm(1) = DELTAE
c*********        Cm(2) = C3Sig
c*********        Cm(3) = C3Pi1
c*********        Cm(4) = C3Pi3
c*********        Cm(5) = C6Sig
c*********        Cm(6) = C6Pi1
c*********        Cm(7) = C6Pi3
c*********        Cm(8) = C8Sig
c*********        Cm(9) = C8Pi1
c*********        Cm(10)= C8Pi3
c***********************************************************************      
c...      Initialize interaction matrix to 0.d0
          DO  I= 1,3
              DO J= 1,3
                  A(I,J)=0.0D0
                  DR(I,J)=0.d0
                  DO  K= 1,NCMMax
                      DMx(K,I,J)=0.d0
                      ENDDO
                  ENDDO
              ENDDO
c...      Prepare interaction matrix  A
          DO  I= 2,NCMM,3
              RPOW(I)= RDIST**MMLR(I)
           A(1,1)= A(1,1) - Dm(I)*(Cm(I)+Cm(I+1)+Cm(I+2))/(3.d0*RPOW(I))
           A(1,2)= A(1,2) - Dm(I)*(Cm(I+2)+Cm(I+1)-2.d0*Cm(I))/(RPOW(I))
           A(1,3)= A(1,3) - Dm(I)*(Cm(I+2)-Cm(I+1))/(RPOW(I))
           A(2,2)= A(2,2) - Dm(I)*(Cm(I+2)+Cm(I+1)+4.d0*Cm(I))
     1                             /(6.d0*RPOW(I))
           A(3,3)= A(3,3) - Dm(I)*(Cm(I+2)+Cm(I+1))/(2.d0*RPOW(I))
           ENDDO
          A(1,1) = A(1,1) + VLIM
          A(1,2) = A(1,2)/(3.d0*DSQRT(2.d0))
          A(2,1) = A(1,2)
          A(2,2) = A(2,2) + VLIM + DELTAE
          A(2,3) = A(1,3)/(2.d0*DSQRT(3.d0))
          A(1,3) = A(1,3)/(DSQRT(6.d0))
          A(3,1) = A(1,3)
          A(3,2) = A(2,3)
          A(3,3) = A(3,3) + VLIM + DELTAE

c...      Prepare radial derivative of interaction matrix (? is it needed ?)
          DO  I= 2,NCMM,3
              DR(1,1)= DR(1,1) + Dm(I)*MMLR(I)*(Cm(I)+Cm(I+1)+Cm(I+2))
     1                             /(3.d0*RPOW(I)*RDIST)
     2                    -Dmp(I)*(Cm(I)+Cm(I+1)+Cm(I+2))/(3.d0*RPOW(I))
              DR(1,2)= DR(1,2) + Dm(I)*MMLR(I)*(Cm(I+2)+Cm(I+1)-2.d0*
     1                           Cm(I))/(RPOW(I)*RDIST)
     2                    -Dmp(I)*(Cm(I+2)+Cm(I+1)-2.d0*Cm(I))/(RPOW(I))
              DR(1,3)= DR(1,3) + Dm(I)*MMLR(I)*(Cm(I+2)-Cm(I+1))
     1                            /(RPOW(I)*RDIST)
     2                        -Dmp(I)*(Cm(I+2)-Cm(I+1))/(RPOW(I))
              DR(2,2)= DR(2,2) + Dm(I)*MMLR(I)*(Cm(I+2)+Cm(I+1)+
     1                           4.d0*Cm(I))/(6.d0*RPOW(I)*RDIST)
     2                        -Dmp(I)*(Cm(I+2)+Cm(I+1)+4.d0*Cm(I))
     3                            /(6.d0*RPOW(I))
              DR(3,3)= DR(3,3) + Dm(I)*MMLR(I)*(Cm(I+2)+Cm(I+1))
     1                            /(2.d0*RPOW(I)*RDIST)
     2                        -Dmp(I)*(Cm(I+2)+Cm(I+1))/(2.d0*RPOW(I)) 
              ENDDO
          DR(1,2) = DR(1,2)/(3.d0*DSQRT(2.d0))
          DR(2,1) = DR(1,2)
          DR(2,3) = DR(1,3)/(2.d0*DSQRT(3.d0))
          DR(1,3) = DR(1,3)/(DSQRT(6.d0))
          DR(3,1) = DR(1,3)
          DR(3,2) = DR(2,3)

c...      Partial derivatives of interaction matrix A  w.r.t.  Cx
          DO  I= 2,NCMM,3 
              DMx(I,1,1)= -Dm(I)/(3.d0*RPOW(I))
              DMx(I+1,1,1)= DMx(I,1,1) 
              DMx(I+2,1,1)= DMx(I,1,1)
              DMx(I,1,2)= 2.d0*Dm(I)/(3.d0*DSQRT(2.d0)*RPOW(I))
              DMx(I+1,1,2)= -DMx(I,1,2)/2.d0
              DMx(I+2,1,2)= DMx(I+1,1,2)
              DMx(I,2,1)= DMx(I,1,2)
              DMx(I+1,2,1)= DMx(I+1,1,2)
              DMx(I+2,2,1)= DMx(I+2,1,2)
              DMx(I,1,3)= 0.d0
              DMx(I,3,1)= 0.d0
              DMx(I+1,1,3)= Dm(I)/(DSQRT(6.d0)*RPOW(I))
              DMx(I+1,3,1)= DMx(I+1,1,3)
              DMx(I+2,1,3)= -DMx(I+1,1,3)
              DMx(I+2,3,1)= DMx(I+2,1,3)
              DMx(I,2,2)= 2.d0*Dm(I)/(3.d0*RPOW(I))
              DMx(I+1,2,2)= DMx(I,2,2)/4.d0
              DMx(I+2,2,2)= DMx(I+1,2,2)
              DMx(I,2,3)= 0.d0
              DMx(I,3,2)= 0.d0
              DMx(I+1,2,3)= Dm(I)/(2.d0*DSQRT(3.d0)*RPOW(I))
              DMx(I+1,3,2)= DMx(I+1,2,3)
              DMx(I+2,2,3)= -DMx(I+1,2,3)
              DMx(I+2,3,2)= DMx(I+2,2,3)
              DMx(I,3,3)= 0.d0
              DMx(I+1,3,3)= Dm(I)/(2.d0*RPOW(I))
              DMx(I+2,3,3)= DMx(I+1,3,3)
              ENDDO
c...      Call subroutine to prepare and invert interaction matrix  A
          CALL DSYEVJ3(A,Q,W)
          L=1
c...      Now - identify the lowest eigenvalue of  A  and label it  L
          DO J=2,3
              IF (W(J) .LT. W(L)) THEN
                  L=J
                  ENDIF
              ENDDO
c...      Identifiy the highest eigenvalue of A and label it H
          H=1 
          DO J=2,3
              IF(W(J).GT.W(H)) THEN
                  H=J
                  ENDIF
              ENDDO
c...      Identify the middle eigenvalue of A and label it M
          M=1 
          DO J=2,3
              IF((J.NE.L).AND.(J.NE.H)) M= J
              ENDDO
c...      Select which eigenvalue to use based on user input
          IF(MMLR(1).EQ.-2) THEN 
              X = L
          ELSEIF(MMLR(1).EQ.-3) THEN
              X = M
          ELSE         
              X = H
              ENDIF

c...      determine ULR and eigenvectors
          ULR= -W(X)
          IF((MMLR(1).EQ.-3).OR.(MMLR(1).EQ.-4)) ULR = ULR + DELTAE
          DO I=1,3      
              EIGVEC(I,1) = Q(I,X)
              ENDDO 
ccc    loop over values of m to determine partial derivatives for each C value
          DO I=2,NCMM
             DMtemp(1:3,1:3) = DMx(I,1:3,1:3) 
             DEIGMtemp= -MATMUL(TRANSPOSE(EIGVEC),MATMUL(DMtemp,EIGVEC))
             dULRdCm(I)= DEIGMtemp(1,1)
             ENDDO
          DEIGR = -MATMUL(TRANSPOSE(EIGVEC),MATMUL(DR,EIGVEC))
          dULRdR= DEIGR(1,1)
c------------------------------------------------------------------------
          ENDIF
c------------------------------------------------------------------------
      RETURN

      CONTAINS
c=======================================================================
      SUBROUTINE DSYEVJ3(A, Q, W)
c ----------------------------------------------------------------------
c** Subroutine to setup and diagonalize the matrix  A  and return 
c   eigenvalues W and eigenvector matrix  Q
      INTEGER N, I, X, Y, R
      PARAMETER (N=3)
      REAL*8 A(3,3), Q(3,3), W(3)
      REAL*8 SD, SO, S, C, T, G, H, Z, THETA, THRESH
c     Initialize Q to the identitity matrix
c --- This loop can be omitted if only the eigenvalues are desired ---
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
      WRITE(6,'("DSYEVJ3: No convergence.")')
      END SUBROUTINE DSYEVJ3
      END SUBROUTINE AFdiag
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
