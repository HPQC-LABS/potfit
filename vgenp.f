c***********************************************************************
      SUBROUTINE VGENP(ISTATE,RDIST,VDIST,BETADIST,IDAT,dVdR,d2VdR)
c***********************************************************************
c** This subroutine will generate function values and derivatives
c   of Morse/Long-Range potentials as required for semiclassical 
c   calculation (with quantum corrections) of virial coefficients and
c   their analytical derivatives in direct hamiltonian fitting
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++  COPYRIGHT 2009-2012  by  R.J. Le Roy and Aleksander Cholewinski ++
c   Dept. of Chemistry, Univ. of Waterloo, Waterloo, Ontario, Canada   +
c    This software may not be sold or any other commercial use made    +
c      of it without the express written permission of the authors.    +
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c           ----- Version of 18 November 2012 -----
c        (after removal of RREFad & RREFad & RREFw)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** On entry:
c   ISTATE  is the electronic state being considered in this CALL.
c   RDIST:  at the 8 input RDIST(i) distances, calculate potl & derivs
c     * return potential function at those points as VDIST, and the
c       first and second radial derivartives as  dVdR  &  d2VdR
c     * return potential function exponent coefficientw as BETADIS
c???  * skip partial derivative calculation if  IDAT.le.0
c     * If RDIST.le.0  calculate partial derivatives at distances
c             given by array RD(i,ISTATE) & return them in array DVtot
c** On entry via common blocks:
c  APSE(s).le.0  to use {p,q}-type exponent polynomial of order Nbeta(s)
c     if APSE(s) > 0 \beta(r) is Pashov spline defined by Nbeta(s) points
c* Nbeta(s) is order of the beta(r) exponent polynomial or # spline points
c    MMLR(j,s)  are long-range inverse-powers for an MLR or DELR potential
c    nPB(s)  the basic value of power p for the beta(r)  exponent function
c    nQB(s)  the power p for the power series expansion variable in beta(r)
c    pAD(s) & qAD(s) the values of power p for adiabatic u(r) BOB functions
c    nNA(s) & qNA(s) the values of power p for centrifugal q(r) BOB functions
c    Qqw(s)  the power defining the radial variable y_{Pqw}(r) in the 
c           Lambda-doubling radial strength function  f_{\Lambda}(r)
c    DE     is the Dissociation Energy for each state.
c    RE     is the Equilibrium Distance for each state.
c    BETA    is the array of potential (exponent) expansion parameters 
c    NDATPT is the number of meshpoints used for the array.
c-----------------------------------------------------------------------
c** On exit via common blocks:
c    R      is the distance array
c    VPOT   is the potential that is generated.
c    BETAFX  is used to contain the  beta(r)  function.
c** Internal partial derivative arrays ...
c    DUADRe & DUBDRe are p.derivs of adiabatic fx. w.r.t.  Re
c    DVDQA & DVDQB are p.derivs of non-adiabatic fx. wrt q_A(i) & q_B(i)
c    DTADRe & DTBDRe are p.derivs of non-adiabatic fx. w.r.t.  Re
c    dVdL & dLDDRe are p.derivatives of f_\lambda(r) w.r.t. beta_i & Re
c    DBDB & DBDRe are p.derives of beta(r) w.r.t. \beta_i & Re, respectively
c
c** Temp:
c    BTEMP     is used to represent the sum used for dV/dRe.
c              is used in GPEF for De calculations.
c    BINF      is used to represent the  beta(\infty)  value.
c    YP        is used to represent (R^p-Re^p)/(R^p+Re^p) or R-Re.
c    XTEMP     is used to represent (uLR/uLR_e)* exp{-BINF*RTEMP}
c    PBTEMP    is used to calculate dV/dBi.
c    PETEMP    is used to calculate dV/dBi.
c    AZERO     is used for the trial exponential calculations.
c    AONE      is used for the trial exponential calculations.
c    ATWO      is used for the trial exponential calculations.
c    AZTEMP    is used in the MMO trial exponential calculations.
c              is used in the GPEF = (a+b)/k
c    AOTEMP    is used in the GPEF = [a(k+1)-b(k-1)]/k
c    ATTEMP    is used in the GPEF = [a^2(k+1)-b^2(k-1)]/k
c    ARTEMP    is used in the GPEF = [a^3(k+1)-b^3(k-1)]/k
c    FSW       is used to represent the MLJ switching function.
c=======================================================================
      INCLUDE 'arrsizes.h'
      INCLUDE 'BLKISOT.h'
      INCLUDE 'BLKDATA.h'
      INCLUDE 'BLKPOT.h'
      INCLUDE 'BLKDVDP.h'
      INCLUDE 'BLKBOB.h'
      INCLUDE 'BLKBOBRF.h'
      INCLUDE 'BLKCOUNT.h'
c-----------------------------------------------------------------------
c** Define local variables ...
      INTEGER I,J,I1,ISTATE,IPV,IPVSTART,ISTART,ISTOP,LAMB2,m,npow,
     1  IDAT, NBAND, IISTP,MMLR1D(NCMMax),KDER
      REAL*8 BTEMP,BINF,RVAL(8),RTEMP,RM2,XTEMP,PBTEMP,PETEMP,RET,
     1 FSW,Xtemp2,Btemp2,BMtemp,BMtemp2,RMF,PBtemp2,C3VAL,C3bar,C6bar,
     2 C6adj,C9adj,YP,YQ,YPA,YPB,YQA,YQB,YPE,YPM,YPMA,YPMB,YPP,YQP,YQPA,
     3 YQPB,REp,Req,RDp,RDq,DYPDRE,DYQDRE,VAL,DVAL,HReP,HReQ,SL,SLB,
     4 AREF,AREFp,AREFq, RE3,RE6,RE8,T0,T0P,T1,ULRe,Scalc,dLULRedCm(9),
     5 dLULRedRe,dLULRedDe,dULRdDe,dULRdCm(9),RD3,RD6,RD8,DVDD,RDIST(8),
     6 VDIST(8),BETADIST,BFCT,JFCT,JFCTLD,RETSig,RETPi,RETp,RETm,
     7 REpADA,REpADB,REqADA,REqADB,D2VAL,dYPdR,
     8 dYPEdR,dYQdR,d2YPdR,d2YQdR,d2YPEdR,RINV,dDULRdR,d2DULRdR,dULRdR,
     9 d2ULRdR,DXTEMP,D2XTEMP,dVdR(8),d2VdR(8),dLULRdR,YPPP,dBdR,d2BdR,
     x dULRdRCm(9),dXdP(HPARMX),dXpdP(HPARMX),dLULRdCm(9),
     y DYPEDRE,dVALdRe,dYBdRe,dBpdRe,DYPpDRE,DYPEpdRE,DYQpDRE,dYBpdRe,
     z xBETA(NbetaMX),rKL(NbetaMX,NbetaMX)
c***********************************************************************
c** Common block for partial derivatives at the one distance RDIST
      REAL*8 dVdPk(HPARMX)
      COMMON /dVdPkBLK/dVdPk
c***********************************************************************
c** Temporary variables for MLR and DELR potentials
      INTEGER MMLRP,IDATLAST
      REAL*8 ULR,dAAdRe,dBBdRe,dVdBtemp,CmVALL,tDm,tDmp,tDmpp,
     1  Dm(NCMMAX),Dmp(NCMMAX),Dmpp(NCMMAX)
ccc   DATA IDATLAST/999999999/
ccc   SAVE IDATLAST,REP,AREF,AREFp,AREFq
c***********************************************************************
c** Initializing variables on first entry for each cycle.
ccc   IF(IDAT.LE.IDATLAST) THEN
ccc       IDATLAST= IDAT
ccc put much of the initialization stuff here to be done once per cycle

      REP= RE(ISTATE)**nPB(ISTATE)
      IF(RREF(ISTATE).LE.0) AREF= RE(ISTATE)
      IF(RREF(ISTATE).GT.0) AREF= RREF(ISTATE)
      AREFp= AREF**nPB(ISTATE)
      AREFq= AREF**nQB(ISTATE)
c** Normally data point starts from 1
      ISTART= 1
      ISTOP= 8
c** When calculating only one potential point
          VDIST= 0.0d0
          BETADIST= 0.d0
      PBTEMP= 0.0d0
      PETEMP= 0.0d0
      DO  I= ISTART,ISTOP
          BETAFX(I,ISTATE)= 0.0d0
          UAR(I,ISTATE)= 0.d0
          UBR(I,ISTATE)= 0.d0
          TAR(I,ISTATE)= 0.d0
          TBR(I,ISTATE)= 0.d0
          UAR(I,ISTATE)= 0.d0
          WRAD(I,ISTATE)= 0.d0
          ENDDO
c** Initialize parameter counter for this state ...
      IPVSTART= POTPARI(ISTATE) - 1
c** First - define values & derivatives of uLR at Re fro MLR potential
      ULRe= 0.d0
      T1= 0.d0
      IF(rhoAB(ISTATE).GT.0.d0) THEN
c ... save uLR powers in a 1D array
          DO  m= 1, NCMM(ISTATE)
              MMLR1D(m)= MMLR(m,ISTATE)
              ENDDO
          KDER= 1 
          CALL dampF(RE(ISTATE),rhoAB(ISTATE),NCMM(ISTATE),
     1      MMLR1D,IDF(ISTATE),IDSTT(ISTATE),KDER,Dm,Dmp,Dmpp)
          ENDIF
      DO  m= 1,NCMM(ISTATE)
          IF(rhoAB(ISTATE).LE.0.d0) THEN
              dLULRedCm(m)= 1.d0/RE(ISTATE)**MMLR(m,ISTATE)
            ELSE
              dLULRedCm(m)= Dm(m)/RE(ISTATE)**MMLR(m,ISTATE)
            ENDIF
          T0= CmVAL(m,ISTATE)*dLULRedCm(m)
          ULRe= ULRe + T0
          T1= T1 + MMLR(m,ISTATE)*T0
          ENDDO
      dLULRedRe= -T1/(ULRe*RE(ISTATE))
      DO  m= 1,NCMM(ISTATE)
          dLULRedCm(m)= dLULRedCm(m)/ULRe
          IF(rhoAB(ISTATE).GT.0) THEN
              dLULRedRe= dLULRedRe + dLULRedCm(m)*Dmp(m)/Dm(m)
              ENDIF
          ENDDO
      BINF= DLOG(2.0d0*DE(ISTATE)/ULRe)
      betaINF(ISTATE)= BINF
      DO  I= ISTART,ISTOP
          RVAL(I)= RDIST(I)
          RINV= 1.d0/RVAL(I)
          RDp= RVAL(I)**nPB(ISTATE)
          RDq= RVAL(I)**nQB(ISTATE)
          YPE= (RDp-REP)/(RDp+REP)
          YP= (RDp-AREFp)/(RDp+AREFp)
          YQ= (RDq-AREFq)/(RDq+AREFq)
          YPM= 1.d0 - YP
          DYPDRE= -0.5d0*nPB(ISTATE)*(1.d0 - YP**2)/RE(ISTATE)
          DYQDRE= -0.5d0*nQB(ISTATE)*(1.d0 - YQ**2)/RE(ISTATE)
          DYPEDRE= -0.5d0*nPB(ISTATE)*(1.d0 - YPE**2)/RE(ISTATE)
          DYPDR= -DYPDRE*RE(ISTATE)*RINV
          DYPEDR= 0.5d0*nPB(ISTATE)*RINV*(1.d0 - YPE**2)
          DYQDR= -DYQDRE*RE(ISTATE)*RINV
          D2YPDR= -DYPDR*RINV*(1.d0 + nPB(ISTATE)*YP)
          D2YPEDR= -DYPEDR*RINV*(1.d0 + nPB(ISTATE)*YPE)
          D2YQDR= -DYQDR*RINV*(1.d0 + nQB(ISTATE)*YQ)
          DYPpDRE= -nPB(ISTATE)*YP*RINV*DYPDRE
          DYPEpDRE= -nPB(ISTATE)*YPE*RINV*DYPEDRE
          DYQpDRE= -nQB(ISTATE)*YQ*RINV*DYQDRE
          D2VAL= 0.d0
          YPP= 1.d0
          DVAL= 0.d0
          DBDB(0,I,ISTATE)= 1.0d0
          VAL= BETA(0,ISTATE) + YQ*BETA(1,ISTATE)
          DVAL= BETA(1,ISTATE)
          npow= Nbeta(ISTATE)
c-------------------------------------------------------------------
          DO  J= 2,npow
c... now calculate power series part of the Morse-like exponent,along
c    with its radial derivatives
              D2VAL= D2VAL + BETA(J,ISTATE)* DBLE(J)
     1                                            *DBLE(J - 1) *YPP
              YPP= YPP*YQ
              DVAL= DVAL + BETA(J,ISTATE)* DBLE(J)* YPP
              YPPP= YPP* YQ
              VAL= VAL + BETA(J,ISTATE)*YPPP
              DBDB(J,I,ISTATE)= YPM*YPPP
              ENDDO
          YPP= YPPP
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c*** DBDB & DBDRe= dBeta/dRe  used in uncertainty calculation in fununc.f
          DBDRe(I,ISTATE)= -YP*dLULRedRe
          dVALdRe= DBDRe(I,ISTATE) + (BINF - VAL)*DYPDRE
     1                                   + (1.d0 - YP)*DVAL*DYQDRE
          IF(RREF(ISTATE).LE.0.d0) DBDRe(I,ISTATE)= dVALdRe
c-----------------------------------------------------------------------
c... now the power series and its radial derivatives are used in the
c    construction of the derivatives  with respect to the parameters
          dBpdRe= DYPpDRE*(BINF - VAL) - DYPDR*dLULRedRe
     1   + (-DYPDR*DYQDRE + (1.d0 - YP)*DYQpDRE - DYPDRE*DYQDR)*DVAL
     2                              + (1.d0 - YP)*DYQDR*DYQDRE*D2VAL
          D2VAL= (BINF - VAL)*D2YPDR - 2.d0*DYPDR*DYQDR*DVAL
     1                   + (1.d0- YP)*(D2YQDR*DVAL + DYQDR**2*D2VAL)
          DVAL= (BINF - VAL)*DYPDR + (1.d0- YP)*DYQDR*DVAL
          VAL= YP*BINF + (1.d0- YP)*VAL
          dBdR= dYPEdR*VAL + YPE*DVAL
          d2BdR= d2YPEdR*VAL + 2.d0*dYPEdR*DVAL + YPE*D2VAL
          dYBdRe= DYPEDRE*VAL + YPE*dVALdRe
          dYBpdRe= VAL*DYPEpDRE + DYPEDRE*DVAL + DYPEDR*dVALdRe
     1                                                  + YPE*dBpdRe
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          BETAFX(I,ISTATE)= VAL
          XTEMP= DEXP(-VAL*YPE)
c** Now begin by generating  uLR(r)
          ULR= 0.d0
c-------------------------------------------------------------------
          dULRdR= 0.d0
          d2ULRdR= 0.d0
          dULRdRCm= 0.d0
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          IF(rhoAB(ISTATE).GT.0.d0) THEN
              KDER= 2
              CALL dampF(RVAL(I),rhoAB(ISTATE),NCMM(ISTATE),MMLR1D,
     1                     IDF(ISTATE),IDSTT(ISTATE),KDER,Dm,Dmp,Dmpp)
              ENDIF
          DO  m= 1,NCMM(ISTATE)
              IF(rhoAB(ISTATE).LE.0.d0) THEN
c-----------------------------------------------------------------------
                  dULRdCm(m)= 1.d0*RINV**MMLR(m,ISTATE)
                  dULRdRCm(m)= -dULRdCm(m)*RINV*DBLE(MMLR(m,ISTATE))
                  dDULRdR= 0.d0
                  d2DULRdR= 0.d0
                ELSE
                  dULRdCm(m)= Dm(m)*RINV**MMLR(m,ISTATE)
                  dULRdRCm(m)= -dULRdCm(m)*RINV *DBLE(MMLR(m,ISTATE)) 
     2                                   + Dmp(m)*RINV**MMLR(m,ISTATE)
                  dDULRdR= Dmp(m)*RINV**MMLR(m,ISTATE)
                  d2DULRdR= Dmpp(m)*RINV**MMLR(m,ISTATE)
                ENDIF
              ULR= ULR + CmVAL(m,ISTATE)*dULRdCm(m)
              dULRdR= dULRdR + CmVAL(m,ISTATE)*(dDULRdR
     1                         - dULRdCm(m)*RINV*DBLE(MMLR(m,ISTATE)))
              d2ULRdR= d2ULRdR + CmVAL(m,ISTATE)*(d2DULRdR
     1   - 2.d0*dDULRdR*RINV*DBLE(MMLR(m,ISTATE)) + dULRdCm(m)*RINV**2
     2               *DBLE(MMLR(m,ISTATE))*DBLE((MMLR(m,ISTATE) + 1)))
              ENDDO
          dLULRdR= dULRdR/ULR
          DO m= 1,NCMM(ISTATE)
              dLULRdCm(m)= dULRdCm(m)/ULR
              ENDDO
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          XTEMP= XTEMP*ULR/ULRe
c... note ... reference energy for each state is asymptote ...
          DVDD= XTEMP*(XTEMP - 2.D0)  
c---      VPOT(I,ISTATE)= DE(ISTATE)*DVDD + VLIM(ISTATE)
c---      VDIST(I)= VPOT(I,ISTATE)
          VDIST(I)= DE(ISTATE)*DVDD + VLIM(ISTATE)
          BETADIST= VAL
          IF(IDAT.LE.0) GO TO 999
c-              ENDDO
          YPP= 2.d0*DE(ISTATE)*(1.0d0-XTEMP)*XTEMP
          IPV= IPVSTART+2
c... derivatives w.r.t R
          DXTEMP= XTEMP*(dLULRdR - dBdR)
          D2XTEMP= XTEMP*(dBdR**2 - d2BdR + (d2ULRdR 
     1                                           - 2*dBdR*dULRdR)/ULR)
          dVdR(I)= 2.d0*DE(ISTATE)*DXTEMP*(XTEMP - 1.d0)
          d2VdR(I)= 2.d0*DE(ISTATE)*(DXTEMP**2 + D2XTEMP
     1                                                *(XTEMP - 1.d0))
c *** This is just to write the derivatives for testing
c             IF(RDIST.LT.0) WRITE (40,640) (RVAL,VVAL,dVdR,d2VdR,
c    1                                                           YVAL)
c 640         FORMAT(G12.5, G18.10, G18.10, G18.10, G14.7)
c ... derivative w.r.t. Cm's
          DO  m= 1, NCMM(ISTATE)
              IPV= IPV+ 1
              dXdP(IPV)= XTEMP*(dLULRdCm(m) + (YPE*YP - 1.d0)
     1                                                  *dLULRedCm(m))
              dXpdP(IPV)= DEXP(-VAL*YPE)/ULRe*(dULRdRCm(m)
     1                    - dBdR*dULRdCm(m)) + (DXTEMP*(YPE*YP - 1.d0)
     2                   + XTEMP*(dYPEdR*YP + YPE*dYPdR))*dLULRedCm(m)
              dVpdP(IPV,I)= 2.d0*DE(ISTATE)*(dXdP(IPV)*DXTEMP
     1                                    + (XTEMP - 1.d0)*dXpdP(IPV))
              DVtot(IPV,I)= -YPP*(dLULRedCm(m)*(YP*YPE- 1.d0)
     1                                               + dULRdCm(m)/ULR)
              ENDDO
c... derivative w.r.t. Re  
          dXdP(IPVSTART+2)= -XTEMP*(dYBdRe + dLULRedRe)
          dXpdP(IPVSTART+2)= -DXTEMP*(dYBdRe + dLULRedRe)
     1                                                 - XTEMP*dYBpdRe
          dVpdP(IPVSTART+2,I)= 2.d0*DE(ISTATE)*(dXdP(IPVSTART+2)
     1                     *DXTEMP + (XTEMP - 1.d0)*dXpdP(IPVSTART+2))
          DVtot(IPVSTART+2,I)= YPP*(dYBdRe + dLULRedRe)
c... derivative w.r.t. De
          dXdP(IPVSTART+1)= -XTEMP*YPE*YP
          dXpdP(IPVSTART+1)= -(XTEMP*(YPE*DYPDR + DYPEDR*YP)
     1                                            + YPE*YP*DXTEMP)
          DVDD= DVDD + YPP*YP*YPE/DE(ISTATE)
          YPP= YPP*YPE*(1.d0 - YP)
          dVpdP(IPVSTART+1,I)= 2.d0*(dXdP(IPVSTART+1)*DXTEMP
     1                             + (XTEMP - 1.d0)*dXpdP(IPVSTART+1)) 
     2                                    + 2.d0*(XTEMP - 1.d0)*DXTEMP
          DVtot(IPVSTART+1,I)= DVDD
c... finally ... derivatives w.r.t. exponent expansion coefficients
          DO  J= 0,npow
              IPV= IPV+1
              dXdP(IPV)= XTEMP*YPE*(1.d0 - YP)*YQ**J
              dXpdP(IPV)= (XTEMP*((1.d0 - YP)*DYPEDR - DYPDR*YQ)
     1                  + YPE*(1.d0 - YP)*DXTEMP)*YQ**J + XTEMP*J*(YPE
     2                                       *(1.d0 - YP))*YQ**(J - 1)
              dVpdP(IPV,I)= 2.d0*DE(ISTATE)*(dXdP(IPV)*DXTEMP
     1                                    + (XTEMP - 1.d0)*dXpdP(IPV))
              DVtot(IPV,I)= YPP
              YPP= YPP*YQ
              ENDDO
          ENDDO
ccccc Print for testing
      rewind(10)
      write(10,610) (RD(i,ISTATE),vpot(i,istate),BETAFX(i,istate),
     1                                        i= 1, NDATPT(ISTATE),50)
  610 FORMAT(/(f10.4,f15.5,f12.6))
ccccc End of Print for testing

      IF((NUA(ISTATE).GE.0).OR.(NUB(ISTATE).GT.0)) THEN
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Treat any 'adiabatic' BOB radial potential functions here ...
c        u_A(r) = yp*uA_\infty + [1 - yp]\sum_{i=0,NUA} {uA_i yq^i} 
c   where the  u_\infty  values stored/fitted as  UA(NUA(ISTATE))
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          REp= RE(ISTATE)**pAD(ISTATE)
          REq= RE(ISTATE)**qAD(ISTATE)
          HReP= 0.5d0*pAD(ISTATE)/RE(ISTATE)
          HReQ= 0.5d0*qAD(ISTATE)/RE(ISTATE)
          REpADA= RE(ISTATE)**pAD(ISTATE)
          REqADA= RE(ISTATE)**qAD(ISTATE)
          REpADB= RE(ISTATE)**pAD(ISTATE)
          REqADB= RE(ISTATE)**qAD(ISTATE)
          IF((BOBCN(ISTATE).GE.1).AND.(pAD(ISTATE).EQ.0)) THEN
              HReP= 2.d0*HReP
              HReQ= 2.d0*HReQ
              ENDIF
c ... reset parameter counter ...
          IPVSTART= IPV
          DO  I= ISTART,ISTOP
              RVAL(I)= RD(I,ISTATE)
              IF(RDIST(I).GT.0.d0) RVAL(I)= RDIST(I)
              RDp= RVAL(I)**pAD(ISTATE)
              RDq= RVAL(I)**qAD(ISTATE)
              YPA= (RDp - REpADA)/(RDp + REpADA)
              YQA= (RDq - REqADA)/(RDq + REqADA)
              YPB= (RDp - REpADB)/(RDp + REpADB)
              YQB= (RDq - REqADB)/(RDq + REqADB)
              YPMA= 1.d0 - YPA
              YPMB= 1.d0 - YPB
              IF(BOBCN(ISTATE).GE.1) THEN
c** If  BOBCN > 0  &  p= 1,  assume use of Ogilvie-Tipping vble.
                  IF(pAD(ISTATE).EQ.1) THEN
                      YPA= 2.d0*YPA
                      YPB= 2.d0*YPB
                      ENDIF
                  ENDIF
              IF(NUA(ISTATE).GE.0) THEN
c ... Now ... derivatives of UA w.r.t. expansion coefficients
                  VAL= UA(0,ISTATE)
                  DVAL= 0.d0
                  IPV= IPVSTART + 1
                  DVtot(IPV,I)= YPMA
                  YQPA= 1.d0
                  IF(NUA(ISTATE).GE.2) THEN
                      DO  J= 1,NUA(ISTATE)-1
                          DVAL= DVAL+ DBLE(J)*YQPA*UA(J,ISTATE)
                          YQPA= YQPA*YQA
                          VAL= VAL+ UA(J,ISTATE)*YQPA
                          IPV= IPV+ 1
                          DVtot(IPV,I)= YPMA*YQPA
                          ENDDO
                      ENDIF
                  IPV= IPV + 1
                  DVtot(IPV,I)= YPA
                  UAR(I,ISTATE)= VAL*YPMA + YPA*UA(NUA(ISTATE),ISTATE)
                  DUADRe(I,ISTATE)= 0.d0
c ... and derivative of UA w.r.t. Re ... 
                   DUADRe(I,ISTATE)= -HReQ*(1.d0 - YQA**2)*YPMA*DVAL
     1            + HReP*(1.d0 - YPA**2)*(VAL- UA(NUA(ISTATE),ISTATE))
                  ENDIF
              IF(NUB(ISTATE).GE.0) THEN
c ... Now ... derivatives of UB w.r.t. expansion coefficients
                  VAL= UB(0,ISTATE)
                  DVAL= 0.d0
                  IF(NUA(ISTATE).LT.0) THEN
                      IPV= IPVSTART + 1
                    ELSE
                      IPV= IPV + 1
                    ENDIF
                  DVtot(IPV,I)= YPMB
                  YQPB= 1.d0
                  IF(NUB(ISTATE).GE.2) THEN
                      DO  J= 1,NUB(ISTATE)-1
                          DVAL= DVAL+ DBLE(J)*YQPB*UB(J,ISTATE)
                          YQPB= YQPB*YQB
                          VAL= VAL+ UB(J,ISTATE)*YQPB
                          IPV= IPV + 1
                          DVtot(IPV,I)= YPMB*YQPB
                          ENDDO
                      ENDIF
                  IPV= IPV + 1
                  DVtot(IPV,I)= YPB
                  UBR(I,ISTATE)= VAL*YPMB + YPB*UB(NUB(ISTATE),ISTATE)
                  DUBDRe(I,ISTATE)= 0.d0
c ... and derivative of UB w.r.t. Re ... 
                  DUBDRe(I,ISTATE)= -HReQ*(1.d0 - YQB**2)*YPMB*DVAL
     1            + HReP*(1.d0 - YPB**2)*(VAL- UB(NUB(ISTATE),ISTATE))
                  ENDIF
              ENDDO
          ENDIF
c++++ END of treatment of adiabatic potential BOB function++++++++++++++

      IF((NTA(ISTATE).GE.0).OR.(NTB(ISTATE).GE.0)) THEN
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Treat any 'non-adiabatic' centrifugal BOB functions here ...
c       q_A(r) = yp*qA_\infty + [1 - yp]\sum_{i=0,NTA} {qA_i yq^i}
c   where the  q_\infty  values stored/fitted as  TA(NTA(ISTATE))
c  Incorporate the 1/r^2 factor into the partial derivatives (but not in
c     the g(r) functions themselves, since pre-SCHRQ takes care of that).
c    Need to add  M_A^{(1)}/M_A^{(\alpha)}  factor later too
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          REp= RE(ISTATE)**pNA(ISTATE)
          REq= RE(ISTATE)**qNA(ISTATE)
          HReP= 0.5d0*pNA(ISTATE)/RE(ISTATE)
          HReQ= 0.5d0*qNA(ISTATE)/RE(ISTATE)
          IF((BOBCN(ISTATE).GE.1).AND.(pNA(ISTATE).EQ.0)) THEN
              HReP= 2.d0*HReP
              HReQ= 2.d0*HReQ
              ENDIF
          IPVSTART= IPV
          DO  I= ISTART,ISTOP
              RVAL(I)= RD(I,ISTATE)
              IF(RDIST(I).GT.0.d0) RVAL(I)= RDIST(I)
              RM2= 1/RVAL(I)**2
              RDp= RVAL(I)**pNA(ISTATE)
              RDq= RVAL(I)**qNA(ISTATE)
              YP= (RDp - REp)/(RDp + REp)
              YQ= (RDq - REq)/(RDq + REq)
              YPM= 1.d0 - YP
              IF(BOBCN(ISTATE).GE.1) THEN
                  YPM= 1.d0
                  YP= 2.d0*YP
                  ENDIF
              IF(NTA(ISTATE).GE.0) THEN 
c ... Now ... derivatives of TA w,r,t, expansion coefficients
                  VAL= TA(0,ISTATE)
                  DVAL= 0.d0
                  IPV= IPVSTART + 1
                  DVtot(IPV,I)= YPM*RM2
                  YQP= 1.d0
                  IF(NTA(ISTATE).GE.2) THEN
                      DO  J= 1,NTA(ISTATE)-1
                          DVAL= DVAL+ DBLE(J)*YQP*TA(J,ISTATE)
                          YQP= YQP*YQ
                          VAL= VAL+ TA(J,ISTATE)*YQP
                          IPV= IPV + 1
                          DVtot(IPV,I)= YPM*YQP*RM2
                          ENDDO
                      ENDIF
                  IPV= IPV + 1
                  DVtot(IPV,I)= YP*RM2
                  TAR(I,ISTATE)= VAL*YPM + YP*TA(NTA(ISTATE),ISTATE)
c ... and derivative of TA w.r.t. Re ...
                  DTADRe(I,ISTATE)= (-HReQ*(1.d0 - YQ**2)*YPM*DVAL
     1        + HReP*(1.d0 - YP**2)*(VAL- TA(NTA(ISTATE),ISTATE)))*RM2
                  ENDIF
              IF(NTB(ISTATE).GE.0) THEN
c ... Now ... derivatives of TB w,r,t, expansion coefficients
                  VAL= TB(0,ISTATE)
                  DVAL= 0.d0
                  IF(NTA(ISTATE).LT.0) THEN
                      IPV= IPVSTART + 1
                    ELSE
                      IPV= IPV + 1
                    ENDIF
                  DVtot(IPV,I)= YPM*RM2
                  YQP= 1.d0
                  IF(NTB(ISTATE).GE.2) THEN
                      DO  J= 1,NTB(ISTATE)-1
                          DVAL= DVAL+ DBLE(J)*YQP*TB(J,ISTATE)
                          YQP= YQP*YQ 
                          VAL= VAL+ TB(J,ISTATE)*YQP
                          IPV= IPV + 1
                          DVtot(IPV,I)= YPM*YQP*RM2
                          ENDDO
                      ENDIF
                  IPV= IPV + 1
                  DVtot(IPV,I)= YP*RM2
                  TBR(I,ISTATE)= VAL*YPM + YP*TB(NTB(ISTATE),ISTATE)
c ... and derivative of TB w.r.t. Re ...
                  DTBDRe(I,ISTATE)= (-HReQ*(1.d0 - YQ**2)*YPM*DVAL
     1        + HReP*(1.d0 - YP**2)*(VAL- TB(NTB(ISTATE),ISTATE)))*RM2
                  ENDIF
              ENDDO
          ENDIF
c++++ END of treatment of non-adiabatic centrifugal BOB function++++++++

      IF(NwCFT(ISTATE).GE.0) THEN
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Treat any Lambda- or 2\Sigma-doubling radial strength functions here
c    representing it as   f(r)= Sum{ w_i * yp^i}
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          LAMB2= 2*IOMEG(ISTATE)
          REP= RE(ISTATE)**Pqw(ISTATE)
          HReP= 0.5d0*Pqw(ISTATE)/RE(ISTATE)
          IPVSTART= IPV
          DO  I= ISTART,ISTOP
              RVAL(I)= RD(I,ISTATE)
              IF(RDIST(I).GT.0.d0) RVAL(I)= RDIST(I)
              RMF= 1.d0/RVAL(I)**2
              IF(IOMEG(ISTATE).GT.0)  RMF= RMF**LAMB2
              RDp= RVAL(I)**Pqw(ISTATE)
              YP= (RDp - REP)/(RDp + REP)
              DVAL= 0.d0
              YPP= RMF
              VAL= wCFT(0,ISTATE)*YPP
              IPV= IPVSTART + 1
              DVtot(IPV,I)= YPP
              IF(NwCFT(ISTATE).GE.1) THEN
                  DO  J= 1,NwCFT(ISTATE)
                      DVAL= DVAL+ DBLE(J)*YPP*wCFT(J,ISTATE)
                      YPP= YPP*YP
                      IPV= IPV + 1
                      DVtot(IPV,I)= YPP
                      VAL= VAL+ wCFT(J,ISTATE)*YPP
                      ENDDO
                  ENDIF
              wRAD(I,ISTATE)= VAL
              dLDDRe(I,NSTATEMX)= -HReP*(1.d0 - YP**2)*DVAL
              ENDDO
          ENDIF
c++++ END of treatment of Lambda/2-sigma centrifugal BOB function+++++++

c++++ Test for inner wall inflection above the asymptote, and if it ++++
c++++ occurs, replace inward potential with linear approximation +++++++
cc    I1= (RE(ISTATE)-RD(1,ISTATE))/(RD(2,ISTATE)-RD(1,ISTATE))
cc    IF(I1.GT.3) THEN
cc        SL= 0.d0
cc        DO  I= I1-2, 1, -1
cc            SLB= SL
cc            SL= VPOT(I,ISTATE) - VPOT(I+1,ISTATE)
cc            IF((SL.LE.SLB).AND.(VPOT(I,ISTATE).GT.VLIM(ISTATE))) THEN
cc                DO  J= I,1,-1
cc                    VPOT(J,ISTATE)= VPOT(I,ISTATE) + (I-J)*SL
cc                    ENDDO
cc                WRITE(6,606) SLABL(ISTATE),RD(I,ISTATE),VPOT(I,ISTATE)
cc                GOTO 66
cc                ENDIF
cc            ENDDO
cc        ENDIF
cc 66 CONTINUE
cc606 FORMAT(9('===')/'!!!! Extrapolate to correct ',A2,' inner-wall inf
cc   1lection at   R=',f6.4,'   V=',f8.0/9('==='))
c++++++++++++End of Inner Wall Test/Correction code+++++++++++++++++++++

c======================================================================
c** At the one distance RDIST calculate total effective potential VDIST
c  including (!!) centrifugal and Lambda/2Sigma doubling terms,
c  and its partial derivatives w.r.t. Hamiltonian parameters dVdPk.
c** This case only for simulation & fitting of tunneling width data.
c
      DO I= 1,8
      IF((RDIST(I).GT.0).AND.(IDAT.GT.0)) THEN 
          NBAND= IB(IDAT)
          IISTP= ISTP(NBAND)
cccccccc
c         WRITE (40,644) IISTP,RDIST,RVAL,VDIST,I,NDATPT(ISTATE)
c 644 FORMAT ('IISTP =',I3,' RDIST =',G16.8,' RVAL =',G16.8,
c    &          ' VDIST =',G16.8,' I =',I6,' NDATPT =',I6)
cccccccc
          BFCT= 16.857629205d0/(ZMASS(3,IISTP)*RDIST(I)**2)
          JFCT= DBLE(JPP(IDAT)*(JPP(IDAT)+1))
          IF(IOMEG(ISTATE).GT.0) JFCT= JFCT - IOMEG(ISTATE)**2
          IF(IOMEG(ISTATE).EQ.-2) JFCT= JFCT + 2.D0
          JFCT= JFCT*BFCT
c ... First get total effective potential, including BOB terms
          VDIST(I)= VDIST(I) + JFCT
          IF(NUA(ISTATE).GE.0) VDIST= VDIST 
     1                          + ZMUA(IISTP,ISTATE)*UAR(ISTOP,ISTATE)
          IF(NUB(ISTATE).GE.0) VDIST= VDIST 
     1                          + ZMUB(IISTP,ISTATE)*UBR(ISTOP,ISTATE)
          IF(NTA(ISTATE).GE.0) VDIST= VDIST 
     1                     + JFCT*ZMTA(IISTP,ISTATE)*TAR(ISTOP,ISTATE)
          IF(NTB(ISTATE).GE.0) VDIST= VDIST 
     1                     + JFCT*ZMTB(IISTP,ISTATE)*TBR(ISTOP,ISTATE)
          JFCTLD= 0.d0
          IF(IOMEG(ISTATE).NE.0) THEN
              IF(IOMEG(ISTATE).GT.0) THEN
c ... for Lambda doubling case ...
                  JFCTLD= (EFPP(IDAT)-EFREF(ISTATE)) 
     1         *(DBLE(JPP(IDAT)*(JPP(IDAT)+1))*BFCT**2)**IOMEG(ISTATE)
                  ENDIF
              IF(IOMEG(ISTATE).EQ.-1) THEN
c ... for doublet Sigma doubling case ...
                  IF(EFPP(IDAT).GT.0) JFCTLD= 0.5d0*JPP(IDAT)*BFCT
                  IF(EFPP(IDAT).EQ.0) JFCTLD= 0.d0
                  IF(EFPP(IDAT).LT.0) JFCTLD= -0.5d0*(JPP(IDAT)+1)*BFCT
                  ENDIF
              VDIST(I)= VDIST(I) + JFCTLD* WRAD(ISTOP,ISTATE)
             ENDIF
cccccccc
c       WRITE (40,648) JPP(IDAT),EFPP(IDAT),RDIST,VDIST
c 648   FORMAT ('J =',I3,' efPARITY =',I3,' RDIST =',G16.8,' VDIST =',
c    1                                                           G16.8/)
cccccccc
              DO  IPV= 1,TOTPOTPAR
                  dVdPk(IPV)= 0.d0
                  ENDDO
c** Now ... generate requisite partial derivatives.
              DO IPV= POTPARI(ISTATE),POTPARF(ISTATE)
                  dVdPk(IPV)= DVtot(IPV,ISTOP)
                  ENDDO
              IF(NUA(ISTATE).GE.0) THEN
                  DO IPV= UAPARI(ISTATE),UAPARF(ISTATE)
                      dVdPk(IPV)= ZMUA(IISTP,ISTATE)*DVtot(IPV,ISTOP)
                      ENDDO
                  ENDIF
              IF(NUB(ISTATE).GE.0) THEN
                  DO IPV= UBPARI(ISTATE),UBPARF(ISTATE)
                      dVdPk(IPV)= ZMUB(IISTP,ISTATE)*DVtot(IPV,ISTOP)
                      ENDDO
                  ENDIF
              IF(NTA(ISTATE).GE.0) THEN
                  DO IPV= TAPARI(ISTATE),TAPARF(ISTATE)
                     dVdPk(IPV)=JFCT*ZMTA(IISTP,ISTATE)*DVtot(IPV,ISTOP)
                     ENDDO
                  ENDIF
              IF(NTB(ISTATE).GE.0) THEN
                  DO IPV= TBPARI(ISTATE),TBPARF(ISTATE)
                     dVdPk(IPV)=JFCT*ZMTB(IISTP,ISTATE)*DVtot(IPV,ISTOP)
                     ENDDO
                  ENDIF
              IF(NwCFT(ISTATE).GE.0) THEN
                  DO IPV= LDPARI(ISTATE),LDPARF(ISTATE)
                      dVdPk(IPV)= JFCTLD*DVtot(IPV,ISTOP)
                      ENDDO
                  ENDIF
          ENDIF
          ENDDO
c*****7********************** BLOCK END ******************************72
  999 RETURN
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
