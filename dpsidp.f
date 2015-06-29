c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE dPSIdp(ISTATE,IISTP,EO,NBEG,NEND,NDIMR,BvWN,V,WF0,
     1                                              dEdPk,dBdPk,dVdPk)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c  Subroutine for calculating the partial derivative of a vibrational
c  eigenfunction w.r.t. one of the fitting parmeters determining its
c  PEC.  Solves the linear inhomogeneous differential equation for
c  WFk= d\PSI/dp_k using the approach that was formulated by J.M. Hutson 
c  [J.Phys.B14, 851 (1982)] for calculating centrifugal distortion 
c  constants of a diatomic molecule.  We use the improved algorithm of
c  J. Tellinghuisen [J.Mol.Spectrosc. 122, 455 (1987)].  
c**  Use this partial derivative to calculate  dBdPk  the partial 
c    derivative of  Bv  for that level w.r.t parmerer p_k
c
c** On entry:   EO    is the eigenvalue (in units [cm-1])
c              WF0     is the zero'th order wavefunction array
c            NBEG & NEND   are the end points of the range of  WF0
c             dEdPk  is the eigenvalue partial derivative: dE/dp_k
c              dVdPk  is the array of  dV(r)/dp_k values
c               RH    is the integration stepsize (in units [Ang])
c               V(i)  is the effective potential (including centrifugal
c                     term if calculation performed at  J > 0) in 
c                     'internal' units, including the factor  RH**2/BvWN
c** return:   dBdPk  partial deriv of Bv w.r.t. param Pk ???
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c               COPYRIGHT 2013  by  Robert J. Le Roy              +
c   Dept. of Chemistry, Univ. of Waterloo, Waterloo, Ontario, Canada   +
c    This software may not be sold or any other commercial use made    +
c      of it without the express written permission of the author.     +
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c  Based of subdoutine CDJOEL for calculating centrifugal distortion 
c  constants using approach of J. Tellinghuisen [J.Mol.Spec.122,455(1987]
c                        Version of 04/18/2013  ?? did it ever work ??
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Dimension:  potential arrays  and  vib. level arrays.
      INTEGER I,ISTATE,IISTP,M,M1,M2,NBEG,NEND,NDIMR
      REAL*8 V(NEND),WF0(NEND),dVdPk(NEND),P(NDIMR),WF1(NDIMR)
      REAL*8 BvWN,Bv,dBdPk,DV,EO,E,RM2,RHSQ,DR,ZTW,AR,G1,G2,G3,P0,P1,
     1  P2,P3,PI,PIF,PRS,PRT,V1,V2,V3,Y1,Y2,Y3,OV,OV01,OV11,PER01,PER11,
     2  AMB,AMB1,AMB2,DEDPK,DDEDPK,gCENT
c-----------------------------------------------------------------------
      INCLUDE 'arrsizes.h'
      INCLUDE 'BLKPOT.h'
      INCLUDE 'BLKISOT.h'
      INCLUDE 'BLKBOBRF.h'
c-----------------------------------------------------------------------
c
      IF(NEND.GT.NDIMR) THEN
          WRITE(6,602) NEND,NDIMR
          I= 0
          RETURN
          ENDIF
      ZTW= 1.D0/12.d0
      DR= RH(ISTATE)
      RHSQ = DR*DR
      DV = RHSQ/12.D0
      E= EO*RHSQ/BvWN
      OV01 = 0.D0
      OV11 = 0.D0
      PER01 = 0.D0
      PER11 = 0.D0
c** First, replace the 'dimensionless dEdPk value' to replace the
c  'dimensionless' Bv value (R2IN) used in the CDC calculations 
      DDEDPK= DEDPK*RHSQ/BvWN
   10 P1= 0.D0
      P2= 0.D0
c??? A possible alternate starting condition - identical to zero-order
c     P1= WF0(NEND)
c     P2= WF0(NEND-1)
c----------------------------------------------------
      P(NEND) = P1
      P(NEND-1) = P2
      V1 = V(NEND) - E
      V2 = V(NEND-1) - E
      G1= (dVdPk(NEND) - DDEDPK)*WF0(NEND)
      Y1 = P1*(1.D0 - ZTW*V1) - DV*G1
      G2 = (dVdPk(NEND-1) - DDEDPK)*WF0(NEND-1)
      Y2 = P2*(1.D0 - ZTW*V2) - DV*G2
      M= NEND-1
c** Now - integrate inward from outer end of range
      DO  I = NBEG+2,NEND
          M = M-1
          Y3 = Y2 + Y2 - Y1 + RHSQ*G2 + V2*P2
          G3 = (dVdPk(M) - DDEDPK)*WF0(M)
          V3 = V(M) - E
          P3 = (Y3 + DV*G3)/(1.D0 - ZTW*V3)
          IF(V3.LT.0.D0)  GO TO 32
          P(M) = P3
          Y1 = Y2
          Y2 = Y3
          V2 = V3
          P2 = P3
          G2 = G3
          ENDDO
      GO TO 90
c** Escaped loop at outer turning point:  initialize outward integration
   32 PRS = P3
      PRT = P(M+1)
      P1 = 0.D0
      P2 = 1.D0
c
c     P1 = WF0(NBEG)
c     P2 = WF0(NBEG+1)
c
      P(NBEG) = P1
      P(NBEG+1) = P2
      V1 = V(NBEG) - E
      V2 = V(NBEG+1) - E
      G1= (dVdPk(NBEG) - DDEDPK)*WF0(NBEG)
      Y1 = P1*(1.D0 - ZTW*V1) - DV*G1
      G2 = (dVdPk(NBEG+1) - DDEDPK)*WF0(NBEG+1)
      Y2 = P2*(1.D0 - ZTW*V2) - DV*G2
      AR = 0.D0
      M1 = M+1
c** Now ... integrate outward from inner end of range
      DO  I = NBEG+2,M1
          Y3 = Y2 + Y2 - Y1 + RHSQ*G2 + V2*P2
          P0 = WF0(I)
          G3 = (dVdPk(I) - DDEDPK)*P0
          V3 = V(I) - E
          P3 = (Y3 + DV*G3)/(1.D0 - ZTW*V3)
          P(I) = P3
          Y1 = Y2
          Y2 = Y3
          V2 = V3
          P2 = P3
          G2 = G3
          AR = AR + P0*P3
          ENDDO
c** Average for 2 adjacent mesh points to get Joel's "(a-b)"
      AMB2 = (P3-PRT)/P0
      AMB1 = (P(M)-PRS)/WF0(M)
      AMB = (AMB1+AMB2)*0.5D0
      M2 = M+2
c** Find the rest of the overlap with zero-th order solution ...
      DO  I = M2,NEND
          P0 = WF0(I)
          PI = P(I) + AMB*P0
          P(I) = PI
          AR = AR + PI*P0
          ENDDO
      OV = AR*DR
      Bv= 0.d0
      dBdPk= 0.d0
      DO  I = NBEG,NEND
          P0 = WF0(I)
          RM2= 1.d0/RD(I,ISTATE)**2
c ... Now ... project out contribution of zero'th-order part of solution
          PI = P(I) - OV*P0
          PIF = PI*dVdPk(I)
c** Now - accumulate integrals for Bv and dBv/dPk & orthogonality test
          WF1(I) = PI
          gCENT= (1.0d0 + ZMTA(IISTP,ISTATE)*TAR(I,ISTATE)
     1                         + ZMTB(IISTP,ISTATE)*TBR(I,ISTATE))*RM2
          Bv= Bv + gCENT*P0**2
          dBdPk= dBdPk + WF1(I)*gCENT*WF0(I)
          OV01 = OV01 + PI*P0
          OV11 = OV11 + PI*PI
          PER01 = PER01 + PIF*P0
          PER11 = PER11 + PI*PIF
          ENDDO
      Bv = Bv*DR*BvWN
      dBdPk= 2.d0*dBdPk*DR*BvWN/RHSQ
      IF(DABS(OV01).GT.1.D-9) WRITE(6,604) OV01
      RETURN
   90 WRITE(6,601) EO
      RETURN
  601 FORMAT(' *** ERROR in CDJOEL *** for input energy  E =',f12.4,
     1   '  never reach outer turning point')
  602 FORMAT(/' *** Dimensioning PROBLEM in CDJOEL ***   NEND=',i6,
     1  ' > NDIMR=',i6)
  604 FORMAT(' ** CAUTION ** dPSIdp  orthogonality test gives OV01:',
     1                                                         1Pd9.1)
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

