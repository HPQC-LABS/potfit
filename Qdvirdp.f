c***********************************************************************
      SUBROUTINE QDVIRDP(IDAT,ISTATE,ZMU,BVIR,dBVIRdP)
c=======================================================================
c  This subroutine calculates the second virial coefficient BVIR and its
c  partial derivatives with respect to the various parameters. It is 
c  used when virial data has been input into the program.  It performs a
c  quantum phase shift calculation with quantum corrections
c .................INCOMPLETE.................................
c=======================================================================
      INCLUDE 'arrsizes.h'
      INCLUDE 'BLKPOT.h'
      INCLUDE 'BLKCOUNT.h'
      INCLUDE 'BLKDVDP.h'
      INCLUDE 'BLKISOT.h'
      INCLUDE 'BLKDATA.h'
c* Define types for local variables
       INTEGER check,j,counter,nParams,ISTATE,ISTART,
     1 IDAT,i,m,k,kk,jj,n
       REAL*8  XG(8),WG(8),x(4),w(4),BTemp,BTempInv,
     1 CONST(3),jump,a,b,YVAL(8),RVAL(8),
     2 VDIST(8),EXP_TERM,int_fact,Vsq,VPsq,Rsq,RINV,XTEMP,Bclass,Bq1,
     3 Bq2,INTEGRALS(NPARMX+3),error,BVIR,BETADISTa,dVdR(8),d2VdR(8),
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
       integer NLGP
      NLGP= 15
       REAL*8 xL25(15),wL15(15)
15-point Laguerre Integration points and weights
      data/xL15/0.093307812017d0,0.4926917403002d0,1.215595412071d0,
     1  1.269940526204d0,3.667622721751d00, 5.425336627414d0,
     2  7.565916226613d0, 10.120228568019d0, 13.130282482176d0,
     3  16.654407708330d0, 20.776478899449d0, 25.623894226729d0,
     4  31.407519169754d0, 38.530683306486d0, 48.026085572686d0/
      data/wL15/2.18234885940d-01, 3.42210177923d-01, 2.63027577942d-02
     1 1.26425818106d-01, 4.02068649210d-02, 8.56387780362d-03,
     2 1.21243614721d-03, 1.11674392344d-04, 6.45992676202d-06,
     3 2.22631690710d-07, 4.22743038498d-09, 3.92189726704d-11,
     4 1.45651526407d-13, 1.48302705111d-16, 1.60059490621d-20/
c***********************************************************************
      ERROR= 0.001d0
c.. runs a loop to set all quadrature weights and points
      DO j= 1,4
          XG(j)= -x(j)
          WG(j)= w(j)
          XG(9-j)= x(j)
          WG(9-j)= w(j)
          ENDDO
c** zero contbn. to Q. Virial and transp.
c.. initializes the array dBVIRdP to zero.
      DO J= 1,HPARMX
          dBVIRdP(J)= 0.d0
          ENDDO
c** copy 2-D potential array VPOT  for state ISTATE into 1D array VV and create
XM2=1/RD**2 array

c.. takes the current temperature from the virial data and changes it to E cm
      BTemp = k_cm*Temp(IDAT)
c.. invert
      DO IE= 1,15
      E= xL15 (IE)*BTEMP
c** Loop over partial waves
     Do  L=1,1000
     CALL PHSHWF( .... )


      BTempInv = 1.d0/BTemp
c.. sets the constants for the intergration of each correction term

c.. sets the number of paramters used for the integration and initializes
c.  the integrals to zero
      nParams = 3 + NCMM(ISTATE) + Nbeta(ISTATE)
      check = 1

c.. this first loop is here to repeat the calculations until the values
c** Outermost loop: repeatedly bisect overall [-1,+1] interval NBISMX
c  times until convergence to within absolute error  ERROR  is achieved
      n= 12
      DO  kk= 1,n
          DO J= 1,nParams+3
              INTEGRALS(J)= 0.d0
              ENDDO
          IF(check.EQ.-1) EXIT
          counter= 2.d0**kk

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
c... mapping  r<->y=(r/re - 0.9999)/(r/re + 1.0001) 
c   where 'real' range is  0.01*Re  to infty - with  'p=2'
                  RVAL(i)= Re(ISTATE)*DSQRT((1.0001d0 
     1                           + 0.9999d0*YVAL(i))/(1.d0 - YVAL(i)))
                  VDIST(i)= 0.d0
                  ENDDO
c.. calls VGENp to find the neccesary derivatives
              CALL VGENp(ISTATE,RVAL,VDIST,BETADISTa,IDAT,dVdR,d2VdR)
c.. the next loop is over each Gaussian point within each subinterval
              DO i= 1,8
c.. some of the terms that will be used in later calculations of the
c.  virial coefficients are constructed here
                  EXP_TERM= DEXP(-VDIST(i)*BTempINV)
                  RINV= 1.d0/RVAL(i)
                  int_fact= (Re(ISTATE)/(1.d0 - YVAL(i)))**2
     1                                          *RINV*(b - a)*0.5d0
                  Vsq= VDIST(i)**2
                  VPsq= dVdR(i)**2
                  Rsq= RVAL(i)**2
                  XTEMP= d2VdR(i)**2/10.d0 + VPsq*RINV**2/5.d0
     1      + dVdR(i)**3*BTempInv*RINV/9.d0 - (VPsq*BTempInv)**2/72.d0
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c*     Finally the relevant partial derivatives of the virial 
c*     coefficients are calculated
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
                  ENDDO
              ENDDO
              
          DO j= 1,nParams
                  dBVIRdP(j)= INTEGRALS(j)
                  ENDDO
          dBVIR=BVIR
          BVIR= Const(1)*INTEGRALS(nParams+1) + 
     2    Const(2)*INTEGRALS(nParams+2) +
     3    Const(3)*Integrals(nParams+3)

c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c. Check to see if BVIR has converged
          IF(counter.GE.2) THEN
              error=0.000001
              IF(DABS(BVIR-dBVIR).LE.(abs(error*BVIR))) THEN
                  check= -1 
                  ENDIF
              ENDIF 
          ENDDO
      END 
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
