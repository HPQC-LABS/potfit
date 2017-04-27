c***********************************************************************
      SUBROUTINE DVIRDP(IDAT,ISTATE,ZMU,BVIR,dBVIRdP)
c=======================================================================
c  This subroutine calculates the second virial coefficient BVIR and its
c  partial derivatives with respect to the various parameters. It is 
c  used when virial data has been input into the program.  It performs a
c  classical calcultion with two quantum corrections.   Update: 15/05/16
c=======================================================================
      INCLUDE 'arrsizes.h'
      INCLUDE 'BLKPOT.h'
      INCLUDE 'BLKCOUNT.h'
      INCLUDE 'BLKDVDP.h'
      INCLUDE 'BLKISOT.h'
      INCLUDE 'BLKDATA.h'
c* Define types for local variables
      INTEGER check,j,counter,nParams,ISTATE,ISTART,IDAT,i,m,k,kk,jj,n,
     1  VIRCNT(NSTATEMX) 
      REAL*8  XG(8),WG(8),x(4),w(4),BTemp,BTempInv,
     1 CONST(3),jump,a,b,YVAL(8),RDIST(8),Class,Q1corr,Q2corr,
     2 VDIST(8),EXP_TERM,int_fact,Vsq,VPsq,Rsq,RINV,XTEMP,Bclass,Bq1,
     3 Bq2,INTEGRALS(NPARMX+3),error,BVIR,dVdR(8),d2VdR(8),
     4 XTEMP1,dBcdP(NPARMX),dBq1dP(NPARMX),dBVIRdP(NPARMX),ZMU,dBVIR
c
       REAL*8, PARAMETER :: k_boltz=1.3806488D-23, NA=6.02214129D23,
     1  h = 6.62606957D-34, c = 2.99792458D10, k_cm = k_boltz/(h*c),
     2  h_cm = 1.d0/c, pi = 3.14159265358979323846, amu=1.660538921D-27,
     3  Cu= 16.857629206D0
c* Common block data for quadrature weights and points
      data x/0.960289856497536d0, 0.796666477413627d0,
     1       0.525532409916329d0, 0.183434642495650d0/,
     2     w/0.101228536290376d0, 0.222381034453374d0,
     3       0.313706645877887d0, 0.362683783378362d0/
c***********************************************************************
      ERROR= 0.001d0
      IF(VIRCNT(ISTATE).EQ.0) THEN
          VIRCNT(ISTATE)= VIRCNT(ISTATE) + 1
          WRITE(32,100)
          ENDIF
      ISTART = POTPARI(ISTATE)- 1
c.. runs a loop to set all quadrature weights and points
      DO j= 1,4
          XG(j)= -x(j)
          WG(j)= w(j)
          XG(9-j)= x(j)
          WG(9-j)= w(j)
          ENDDO
c.. initializes the array dBVIRdP to zero.
      DO J= 1,HPARMX+3
          dBVIRdP(J)= 0.d0
          ENDDO
c.. takes the current temperature from the virial data and changes it to E cm
      BTemp = k_cm*Temp(IDAT)
c.. inverts BTemp
      BTempInv = 1.d0/BTemp
c.. sets the constants for the intergration of each correction term
      Const(1)= -2.d0*pi*0.6022140857d0           !! updated Avogadro No.
      Const(2)= pi*0.6022140857d0*Cu*BTempInv**3/(6.d0*ZMU)
      Const(3)= -(pi*0.6022140857d0/6.d0)*(Cu*BTempInv**2/ZMU)**2
c.. sets the number of paramters used for the integration and initializes
c.  the integrals to zero
      nParams = 3 + NCMM(ISTATE) + Nbeta(ISTATE)
      check = 1
c.. this first loop is here to repeat the calculations until the values
c** Outermost loop: repeatedly bisect overall [-1,+1] interval NBISMX
c  times until convergence to within absolute error  ERROR  is achieved
      n= 12
      DO kk= 1,n
          DO J= 1,nParams+3
              INTEGRALS(J)= 0.d0
              ENDDO
          IF(check.EQ.-1) EXIT
          counter= 2**kk
c.. futher subdivides the interval for gaussian integration formula
          jump= 2.d0/DBLE(counter)
          b= -1.d0
c. the first loop is over each subinterval within (-1,1)
          DO m= 1,counter
              a= b
              b= a+ jump
              IF(m.EQ.counter) b= 1.d0
              DO i= 1,8
c.. sets the y value for the gaussian formula
                  YVAL(i)= 0.5d0*((b - a)*XG(i) + (b + a))
c... q=1 mapping  r <-> y= (r/re - 1)/(r/re + 1) 
                  RDIST(i)= Re(ISTATE)*DSQRT((1.0d0 + YVAL(i))/
     1                                                (1.d0 - YVAL(i)))
                  VDIST(i)= 0.d0
                  ENDDO
c.. calls VGENP to find the neccesary PEC value and radial derivatives
              CALL VGENP(ISTATE,RDIST,VDIST,dVdR,d2VdR,IDAT)
c.. the next loop over each Gaussian point within each subinterval
              DO  i= 1,8
c.. some of the terms that will be used in later calculations of the
c.  virial coefficients are constructed here
                  IF(RDIST(I).LT.0.8d0*RE(ISTATE)) THEN
c** Check for and correct for potential turnover problems
                      IF((VDIST(I).LT.0.d0).OR.(DVDR(I).GT.0.d0)
c... or for potential becoming singular at very short dist (r < 0.1 Ang)
     1             .OR.(RDIST(I).LE.0.1d0).OR.(D2VDR(I).LT.0.d0)) THEN
                          VDIST(I)= 1.d6
                          dVdR(I)= -1.d6
                          d2VdR(I)= 1.d6
                          ENDIF
                      ENDIF
                  EXP_TERM= DEXP(-VDIST(i)*BTempINV)
                  RINV= 1.d0/RDIST(i)
                  int_fact= (Re(ISTATE)/(1.d0 - YVAL(i)))**2
     1                                     *RINV*(b - a)*0.5d0
                  Vsq= VDIST(i)**2
                  VPsq= dVdR(i)**2
                  Rsq= RDIST(i)**2
                  XTEMP= d2VdR(i)**2/10.d0 + VPsq*RINV**2/5.d0
     1      + dVdR(i)**3*BTempInv*RINV/9.d0 - (VPsq*BTempInv)**2/72.d0
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c*     Finally the relevant partial derivatives of the virial 
c*     coefficients w.r.t. potential param are calculated
                  DO J = 1,nParams
                      jj = ISTART + J
c* the derivative of the classical expression
                      dBcdP(J)= -BTempInv*EXP_TERM
     1                                   *DVtot(jj,i)*Rsq*int_fact
c* and of the first quantum correction 
                      XTEMP1 = -VPsq*BTempInv*DVtot(jj,i)
     1                                     + 2*dVdR(i)*dVpdP(jj,i)
                      dBq1dP(J)= EXP_TERM*XTEMP1*Rsq*int_fact
c..  As the final step these terms all added together in a weighted sum
                      INTEGRALS(J) = INTEGRALS(J) + Const(1)
     1                  *dBcdP(J)*WG(I) + Const(2)*dBq1dP(J)*WG(I)
                      ENDDO
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c. now the integrands are evaluated at each particular Gaussian point
c. and summed together with the proper weighting
                  Bclass= (EXP_TERM - 1.d0)*Rsq*int_fact
                  INTEGRALS(nParams+1)= INTEGRALS(NParams+1) 
     1                                                 + Bclass*WG(i)
                  Bq1= EXP_TERM*VPsq*Rsq*int_fact
                  INTEGRALS(nParams+2)= INTEGRALS(nParams+2) 
     1                                                    + Bq1*WG(i)
                  Bq2= EXP_TERM*XTEMP*Rsq*int_fact
                  INTEGRALS(nParams+3)= INTEGRALS(nParams+3)
     1                                                    + Bq2*WG(i)
                  ENDDO     !! end of 'i' loop over 8 gaussian points 
              ENDDO    !! end of 'm' loop over 'counter' radial subintervals
          IF(PSEL(ISTATE).LE.3) THEN
                  DO j= 1,nParams      ! partial derivatives for fit cases
                  dBVIRdP(j)= INTEGRALS(j)
                  ENDDO
              ENDIF
          dBVIR=BVIR
          Class= Const(1)*INTEGRALS(nParams+1)
          Q1corr= Const(2)*INTEGRALS(nParams+2)
          Q2corr= Const(3)*Integrals(nParams+3)
          BVIR= Class + Q1corr+ Q2corr
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c. Check to see if BVIR has converged
          IF(counter.GE.2) THEN
              error=0.000001
              IF(DABS(BVIR-dBVIR).LE.(abs(error*BVIR))) THEN
                  check= -1 
                  ENDIF
              ENDIF 
          ENDDO               !! end of 'kk' loop over # sets of subdivisions
      WRITE(32,102) Temp(IDAT),Class,Q1corr,Q2corr,BVIR,counter
  100 FORMAT('Temperature Classical   FirstQ   SecondQ   Total  counter
     1', /24('==='))
  102 FORMAT(2X,F7.2,2X,F8.3,2X,F8.3,2X,F8.3,2X,F8.3,I6)
      RETURN
      END 
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

