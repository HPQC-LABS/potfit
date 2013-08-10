c***********************************************************************
      SUBROUTINE VGEN(ISTATE,RDIST,VDIST,BETADIST,IDAT)
c***********************************************************************
c** This subroutine will generate one of six possible families of
c   analytical molecular potentials for the direct hamiltonian fitting
c   program, as specified by parameter PSEL
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++  COPYRIGHT 2006-2011  by  R.J. Le Roy, Jenning Seto, Yiye Huang  ++
c  and N. Dattani, Department of Chemistry, University of of Waterloo, +
c                        Waterloo, Ontario, Canada                     +
c    This software may not be sold or any other commercial use made    +
c      of it without the express written permission of the authors.    +
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c               ----- Version of  09 February 2013 -----
c           (after removal or RREFna, RREFad & RREFw variables!)
c  (and aftter replacing NLpow by Nbeta & NSpow by APSE 
c                 not finished for DELR - see L. 624 !!!
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** On entry:
c   ISTATE  is the electronic state being considered in this CALL.
c   RDIST: If RDIST > 0, calculate potl & derivs only @ that onee distance
c   -----    * return potential function at that point as VDIST and
c              potential function exponent as BETADIST
c            * skip partial derivative calculation if  IDAT.le.0
c        * If RDIST.le.0  calculate partial derivatives at distances
c             given by array RD(i,ISTATE) & return them in array DVtot
c** On entry via common blocks: 
c* PSEL specifies how data for state ISTATE are to be represented!
c     PSEL = -2 : Represent {v,J,p} levels of state ISTATE by term values
c     PSEL = -1 : Represent {v,J} levels with band constants
c     PSEL = 0 : Use a fixed potential function defined in READPOT.
c     PSEL = 1 : Use an Expanded Morse Oscillator(p) potential.
c     PSEL = 2 : Use a Morse/Lennard-Jones(p) potential.
c     PSEL = 3 : Use a Double-Exponential Long-Range Potential.
c     PSEL = 4 : Use a [Tiemann Polynomial]
c     PSEL = 5 : Use a Tang-Toennies potential
c     PSEL = 6 : Use a Surkus "Generalized Potential Energy Function".
c  APSE(s).le.0  to use {p,q}-type exponent polynomial of order Nbeta(s)
c     if APSE(s) > 0 \beta(r) is Pashov spline defined by Nbeta(s) points
c* Nbeta(s) is order of the beta(r) exponent polynomial or # spline points
c*  MMLR(j,s)  are long-range inverse-powers for an MLR or DELR potential
c*  nPB(s)  the basic value of power p for the beta(r)  exponent function
c*  nQB(s)  the power q for the power series expansion variable in beta(r)
c*  pAD(s) & qAD(s) the values of power p for adiabatic u(r) BOB functions
c*  nNA(s) & qNA(s) the values of power p for centrifugal q(r) BOB functions
c*  Pqw(s)  the power of r defining the y_{Pqw}(r) expansion variable in 
c*       the f_{Lambda}(r) strength function
c*  DE     is the Dissociation Energy for each state.
c*  RE     is the Equilibrium Distance for each state.
c*  BETA    is the array of potential (exponent) expansion parameters 
c*  NDATPT is the number of meshpoints used for the array.
c-----------------------------------------------------------------------
c** On exit via common blocks:
c->  R      is the distance array
c->  VPOT   is the potential that is generated.
c->  BETAFX  is used to contain the  beta(r)  function.
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
c    XTEMP     is used to represent (uLR/uLR_e)* exp{-beta*RTEMP}
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
c** Common block for partial derivatives of potential at the one distance RDIST
      REAL*8 dVdPk(HPARMX)
      COMMON /dVdPkBLK/dVdPk
c=======================================================================
c** Define local variables ...
      INTEGER I,J,I1,ISTATE,IPV,IPVSTART,ISTART,ISTOP,LAMB2,m,npow,
     1  IDAT, NBAND, IISTP, KDER, MMLR1D(NCMMax)
      REAL*8 BTEMP,BINF,RVAL,RTEMP,RM2,XTEMP,PBTEMP,PETEMP,RET,Xtemp2,
     1 Btemp2,BMtemp,BMtemp2,RMF,PBtemp2,C3VAL,C3bar,C6bar,C6adj,C6Pi,
     2 C9adj,C8VAL,C8Pi,YP,YQ,YPA,YPB,YQA,YQB,YPE,YPM,YPMA,YPMB,YPP,YQP,
     3 REP,RDp,RDq,DYPDRE,DYQDRE,VAL,DVAL,HReP,HReQ,SL,
     4 SLB,AREF,AREFp,AREFq, RE3,RE6,RE8,RE9,T0,T0P,T0P23,T1,ULRe,Scalc,
     5 dLULRedCm(9),dLULRedRe,dLULRedDe,dULRdDe,dULRdCm(9),RD3,RD6,RD8,
     6 RD9,DVDD,RDIST,VDIST,BETADIST,BFCT,JFCT,JFCTLD,RETSig,RETPi,RETp,
     7 RETm,REadAp,REadBp,REadAq,REadBq,REnaAp,REnaBp,REnaAq,REnaBq,
     8 REwp,dC6dDe,dC9dC3,dC9dC6,dC9dDe,BT,Rinn,Rout,A1,A2,A3,B1,B2,B3,
     9 B4,B5,xBETA(NbetaMX),rKL(NbetaMX,NbetaMX),C1LIM,BETA0,BETAN,TM,
     x VATT,DTT,ATT,BTT,REQ,
     x DEIGM1(1,1),DEIGM3(1,1),DEIGM5(1,1),DEIGR(1,1),DEIGRe(1,1),
     y DEIGDe(1,1)

c... temporary values for checking 3x3
ccc    real*8  ulr1,ulr2,ulr3,ulr4

c***********************************************************************
c** Temporary variables for MLR and DELR potentials
      INTEGER MMLRP
      REAL*8 ULR,dAAdRe,dBBdRe,dVdBtemp,CmVALL, Dm(NCMMax),Dmp(NCMMax),
     1  Dmpp(NCMMax)
c***********************************************************************
c** Initializing variables.
      REP= RE(ISTATE)**nPB(ISTATE)
      IF(RREF(ISTATE).LE.0) AREF= RE(ISTATE)
      IF(RREF(ISTATE).GT.0) AREF= RREF(ISTATE)
      AREFp= AREF**nPB(ISTATE)
      AREFq= AREF**nQB(ISTATE)
c** Normally data point starts from 1
      ISTART= 1
      ISTOP= NDATPT(ISTATE)
c** When calculating only one potential point
      IF(RDIST.GT.0) THEN
          ISTART= NPNTMX
          ISTOP= NPNTMX
          VDIST= 0.0d0
          BETADIST= 0.d0
          ENDIF
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
      IF(PSEL(ISTATE).EQ.0) THEN
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** For forward data simulation using a known input pointwise potential 
c?????? Prepare distance array ... ?? {Should it then call PPREPOT ?}
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
ccc       DO  I= ISTART,ISTOP
ccc           RVAL= RD(I,ISTATE)
ccc           VDIST= VPOT(I,ISTATE)
ccc           ENDDO
          ENDIF
******7***** End Forward Calculation ***********************************

      IF(PSEL(ISTATE).EQ.1) THEN
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** For the Expanded Morse Oscillator:  exponent polynomial order /Nbeta
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** First ... calculate the Generalized Morse Oscillator exponent
          DO  I= ISTART,ISTOP
              RVAL= RD(I,ISTATE)
              IF(RDIST.GT.0.d0) RVAL= RDIST
              RDp= RVAL**nPB(ISTATE)
              YP= (RDp - AREFp)/(RDp + AREFp)
              VAL= BETA(0,ISTATE) 
              DVAL= 0.d0
              DBDB(0,I,ISTATE)= 1.0d0
              YPP= 1.d0
              DO  J= 1, Nbeta(ISTATE)
                  DVAL= DVAL + BETA(J,ISTATE)* DBLE(J)* YPP
                  YPP= YPP*YP
                  VAL= VAL + BETA(J,ISTATE)*YPP 
                  DBDB(J,I,ISTATE)= YPP
                  ENDDO
c*** DBDB & DBDRe= dBeta/dRe  used in uncertainty calculation in fununc.f
              DBDRe(I,ISTATE)= 0.d0
              BETAFX(I,ISTATE)= VAL
              XTEMP= DEXP(-VAL*(RVAL-RE(ISTATE)))
c** First calculate the partial derivative w.r.t. DE
              IPV= IPVSTART+ 1
              DVtot(IPV,I)= XTEMP*(XTEMP- 2.d0)
              DVDD= DVtot(IPV,I)
c** Now calculate the actual potential
              VPOT(I,ISTATE)= DE(ISTATE)*DVDD + VLIM(ISTATE)
              IF(RDIST.GT.0.d0) THEN
                  VDIST= VPOT(I,ISTATE)
                  BETADIST= VAL
                  IF(IDAT.LE.0) GO TO 999
                  ENDIF
c... now generate remaining partial derivatives
              YPP= 2.0d0 * DE(ISTATE) * XTEMP * (1.d0 - XTEMP) 
              IF(RREF(ISTATE).LE.0.d0) THEN
                  DBDRe(I,ISTATE)= -0.5d0*nPB(ISTATE)*(1.d0-YP**2)
     1                                                *DVAL/RE(ISTATE)
                  VAL= VAL - (RVAL- RE(ISTATE))*DBDRe(I,ISTATE)
                  ENDIF
              IPV= IPV+1
              DVtot(IPV,I)= -YPP*VAL
              YPP= YPP*(RVAL - RE(ISTATE))
              DO  J= 0, Nbeta(ISTATE)
                  IPV= IPV+1
                  DVtot(IPV,I)= YPP
                  YPP= YPP*YP
                  ENDDO
              ENDDO
          ENDIF
ccccc Print for testing
        rewind(10)
        write(10,610) (RD(i,ISTATE),vpot(i,istate),BETAFX(i,istate),
     1                                        i= 1, NDATPT(ISTATE),50)
ccccc
c******** End preparation of Expanded Morse Potential Function *********

      IF(PSEL(ISTATE).EQ.2) THEN
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** For the  {Morse/Long-Range}_p  potential.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** First - define value & derivatives of uLR at Re
          IF((NCMM(ISTATE).GE.4).AND.(MMLR(2,ISTATE).LE.0)) THEN
c** First - define  C6adj & C9adj for the 2x2 and 3x3 Aubert-Frecon cases
              C3VAL= CmVAL(1,ISTATE)
              C3bar= C3VAL/DE(ISTATE)
              C6bar= CmVAL(3,ISTATE)/DE(ISTATE)
              C6adj= CmVAL(3,ISTATE) + 0.25D0*C3VAL*C3bar
              C6Pi= CmVAL(4,ISTATE)
              C8val= CmVAL(5,ISTATE)
              C8Pi= CmVAL(6,ISTATE)
              C9adj= 0.5d0*C3bar*C6adj
              dC9dC6= 0.5d0*C3bar
              dC9dC3= 0.5d0*C6bar + 1.5d0*dC9dC6**2
              dC9dDe= - dC9dC6*(C6bar + C3bar*dC9dC6)
              dC6dDe= - dC9dC6**2
              RE3= 1.d0/RE(ISTATE)**3
              RE6= RE3*RE3
              RE8= RE6/RE(ISTATE)**2
              RE9= RE6*RE3
c** Include QED retardation function in C3 terms for Li2(A) !!
c  NOTE ... the numerical factor is  2\pi/\lambda  for this case
c  Ignore effect of retardation on the partial derivatives & C6adj & C9adj.
              RET= 9.36423830d-04*RE(ISTATE)
              RETSig= DCOS(RET) + (RET)*DSIN(RET)
              RETPi= RETSig - RET**2 *DCOS(RET)
              RETp= RETSig + 0.5d0*RETPi
              RETm= RETSig - 0.5d0*RETPi
              IF((MMLR(2,ISTATE).EQ.0).OR.(MMLR(2,ISTATE).EQ.-2)) THEN
c ... for Aubert-Frecon 2x2 treatment of {C3,C6,C8} for Li2 A- or b-state
                  T1= RE3*(C3VAL*RETm + RE3*(C6adj - C6Pi))/3.d0
                  T1= T1+ (C8VAL - C8Pi)*RE8/3.d0
                  T0= DSQRT((T1- CmVAL(2,ISTATE))**2 + 8.d0*T1**2)
                  ULRe= 0.5d0*(-CmVAL(2,ISTATE) + RE3*(C3VAL*RETp
     1                   + RE3*(C6adj + C6Pi)) + T0) + RE9*C9adj 
                  ULRe= ULRe+ 0.5d0*(C8VAL + C8Pi)*RE8
c.... ajdustment for the b-state
                  IF(MMLR(2,ISTATE).EQ.-2) ULRe= ULRe - T0
c ... now collect derivative terms ...
                  T0P= 0.5d0*(9.d0*T1- CmVAL(2,ISTATE))/T0
                  T0P23= 0.5d0 + T0P/3.d0
                  dLULRedCm(1)= RE3*((0.5d0*RETp + T0P*RETm/3.d0)
     1                         + RE3*(dC9dC6*T0P23 + RE3*dC9dC3))/ULRe
                  dLULRedCm(3)= RE6*(T0P23 + RE3*dC9dC6)/ULRe
                  dLULRedCm(4)= RE6*(1.d0 - T0P23)/ULRe
                  dLULRedDe= RE6*(dC6dDe*T0P23 + RE3*dC9dDe)/ULRe
                  dLULRedRe= - RE3*(C3VAL*(1.5d0*RETp + T0P*RETm)
     1    + RE3*(C6adj*6.d0*T0P23 + RE3*9.d0*C9adj))/(RE(ISTATE)*ULRe)
                  dLULRedCm(5)= RE8*T0P23/ULRe
                  dLULRedCm(6)= RE8*(1.d0 - T0P23)/ULRe
                  dLULRedRe= dLULRedRe - RE9*4.d0*(C8VAL + C8Pi
     1                            + 2.d0*T0P*(C8VAL - C8Pi)/3.d0)/ULRe
                  ENDIF
              IF(MMLR(2,ISTATE).EQ.-1) THEN
c ... extension for Li2(c) {3,0,6,6,8,8} 3x3 Aubert-Frecon case ...
                  CALL AF3X3potRet(RE(ISTATE),CmVAL(2,ISTATE),
     1               CmVAL(1,ISTATE),C6adj,CmVAL(5,ISTATE),DE(ISTATE),
     2               ULRe,DEIGM1,DEIGM3,DEIGM5,DEIGR,DEIGDe)
                  ULRe= ULRe + C9adj*RE9
                  dLULRedCm(5)= DEIGM5(1,1)/ULRe
                  dLULRedCm(3)= (DEIGM3(1,1) + dC9dC6*RE9)/ULRe
                  dLULRedCm(1)= (DEIGM1(1,1) + dC9dC3*RE9)/ULRe
                  dLULRedDe= (DEIGDe(1,1) + dC9dDe*RE9)/ULRe
                  dLULRedRe=(DEIGR(1,1)- 9.d0*C9adj*RE9/RE(ISTATE))/ULRe
                  ENDIF
            ELSE
c*** for 'ordinary' NCMM-term MLR uLR(r) ...  with damping [if rhoAB > 0]
              ULRe= 0.d0
              T1= 0.d0
              IF(rhoAB(ISTATE).GT.0.d0) THEN
c ... save uLR powers in a 1D array
                  DO  m= 1, NCMM(ISTATE)
                      MMLR1D(m)= MMLR(m,ISTATE)
                      ENDDO
                  KDER= 1
                  CALL dampF(RE(ISTATE),rhoAB(ISTATE),NCMM(ISTATE),
     1              MMLR1D,IDF(ISTATE),IDSTT(ISTATE),KDER,Dm,Dmp,Dmpp)
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
            ENDIF
          BINF= DLOG(2.0d0*DE(ISTATE)/ULRe)
          betaINF(ISTATE)= BINF
          IF(APSE(ISTATE).GT.0) THEN
c*** For Pashov-natural-spline exponent coefficient ...
              DO  I= 1,Nbeta(ISTATE)
                  xBETA(I)= ypBETA(I,ISTATE)
                  ENDDO
              BETA(Nbeta(ISTATE),ISTATE)= BINF
              CALL Lkoef(Nbeta(ISTATE),xBETA,rKL,NbetaMX)
              ENDIF
          KDER= 0
c-----------------------------------------------------------------------
          DO  I= ISTART,ISTOP
c** Now - generate potential while looping over radial array
              RVAL= RD(I,ISTATE)
              IF(RDIST.GT.0.d0) RVAL= RDIST
              RDp= RVAL**nPB(ISTATE)
              RDq= RVAL**nQB(ISTATE)
              YPE= (RDp-REP)/(RDp+REP)
              YP= (RDp-AREFp)/(RDp+AREFp)
              YQ= (RDq-AREFq)/(RDq+AREFq)
              YPM= 1.d0 - YP
              DYPDRE= -0.5d0*nPB(ISTATE)*(1.d0 - YP**2)/RE(ISTATE)
              DYQDRE= -0.5d0*nQB(ISTATE)*(1.d0 - YQ**2)/RE(ISTATE)
              YPP= 1.d0
              DVAL= 0.d0
              DBDB(0,I,ISTATE)= 1.0d0
              npow= Nbeta(ISTATE)
              IF(APSE(ISTATE).LE.0) THEN
c** For 'conventional' Huang power-series exponent function ...
                  VAL= BETA(0,ISTATE)
                  DO  J= 1, Nbeta(ISTATE)
c... now calculate power series part of the MLR exponent
                      DVAL= DVAL + BETA(J,ISTATE)* DBLE(J)* YPP
                      YPP= YPP*YQ
                      VAL= VAL + BETA(J,ISTATE)*YPP
                      DBDB(J,I,ISTATE)= YPM*YPP
                      ENDDO
c*** DBDB & DBDRe= dBeta/dRe  used in uncertainty calculation in fununc.f
                  DBDRe(I,ISTATE)= -YP*dLULRedRe
                  IF(RREF(ISTATE).LE.0.d0) DBDRe(I,ISTATE)= 
     1                            DBDRe(I,ISTATE)+ (BINF - VAL)*DYPDRE 
     2                                         + (1.d0-YP)*DVAL*DYQDRE
                  VAL= YP*BINF + (1.d0- YP)*VAL
                ELSE
c... now calculate Pashov-spline exponent coefficient & its derivatives
                  VAL= 0.d0
                  DO  J= 1, Nbeta(ISTATE)
                      DBDB(J,I,ISTATE)= 
     1                               Scalc(YP,J,npow,xBETA,rKL,NbetaMX)
                      VAL= VAL+ DBDB(J,I,ISTATE)*BETA(J,ISTATE)
                      ENDDO
                  DBDRe(I,ISTATE)= -DBDB(npow,I,ISTATE)*dLULRedRe
                ENDIF
              BETAFX(I,ISTATE)= VAL
              XTEMP= DEXP(-VAL*YPE)
c** Now begin by generating  uLR(r)
              IF((NCMM(ISTATE).GE.4).AND.(MMLR(2,ISTATE).LE.0)) THEN
c ... generate ULR for Aubert-Frecon Li2(A) {3,0,6,6} type case ...
                  RD3= 1.d0/RVAL**3
                  RD6= RD3*RD3
                  RD9= RD6*RD3
c** Include QED retardation function in definition of C3 ... as above!!
                  RET= 9.36423830d-4*RVAL
                  RETSig= DCOS(RET) + RET*DSIN(RET)
                  RETPi= RETSig - RET**2 *DCOS(RET)
                  RETp= RETSig + 0.5d0*RETPi
                  RETm= RETSig - 0.5d0*RETPi
                  IF((MMLR(2,ISTATE).EQ.0).OR.(MMLR(2,ISTATE).EQ.-2))
     1                                                            THEN
c ... extension for Li2(A) 2x2 {3,0,6,6,8,8} case ...
                      T1= (C3VAL*RETm + (C6adj - C6Pi)*RD3)*RD3/3.d0
                      RD8= RD6/RVAL**2
                      T1= T1+ (C8VAL - C8Pi)*RD8/3.d0
                      T0= DSQRT((T1- CmVAL(2,ISTATE))**2 + 8.d0*T1**2)
                      ULR= 0.5d0*(-CmVAL(2,ISTATE) + RD3*(C3VAL*RETp
     1                   + RD3*(C6adj + C6Pi))) + 0.5d0*T0 + C9adj*RD9
                      ULR= ULR+ 0.5d0*(C8VAL + C8Pi)*RD8
c.....  adjustment for the b-state .........
                      IF(MMLR(2,ISTATE).EQ.-2) ULR= ULR- T0
                      T0P= 0.5d0*(9.d0*T1 - CmVAL(2,ISTATE))/T0
                      T0P23= 0.5d0 + T0P/3.d0
                      dULRdCm(1)= RD3*(0.5d0*RETp + T0P*RETm/3.d0
     1                              + RD3*(dC9dC6*T0P23 + RD3*dC9dC3))
                      dULRdCm(2)= 0.d0
                      dULRdCm(3)= RD6*(T0P23 + RD3*dC9dC6)
                      dULRdCm(4)= RD6*(1.d0 - T0P23)
                      dULRdDe= RD6*(dC6dDe*T0P23 + RD3*dC9dDe)
                      IF(NCMM(ISTATE).GT.4) THEN
                          dULRdCm(5)= RD8*T0P23
                          dULRdCm(6)= RD8*(1.d0 - T0P23)
                          ENDIF
                      ENDIF
                  IF(MMLR(2,ISTATE).EQ.-1) THEN
c ... generate ULR for Aubert-Frecon Li2(c) {3,-1,6,6} type case ...
                      CALL AF3X3potRet(RVAL,CmVAL(2,ISTATE),
     1               CmVAL(1,ISTATE),C6adj,CmVAL(5,ISTATE),DE(ISTATE),
     2                          ULR,DEIGM1,DEIGM3,DEIGM5,DEIGR,DEIGDe)
                      ULR= ULR + C9adj*RD9
                      dULRdCm(5) = DEIGM5(1,1)
                      dULRdCm(3) = DEIGM3(1,1) + dC9dC6*RD9
                      dULRdCm(1) = DEIGM1(1,1) + dC9dC3*RD9
                      dULRdDe    = DEIGDe(1,1) + dC9dDe*RD9
                      ENDIF
c----- End of special Aubert-Frecon Li2 cases ------------------------
                ELSE
c ... for the case of a 'normal' MLR or MLJ function
                  ULR= 0.d0
                  IF(rhoAB(ISTATE).GT.0.d0) 
     1            CALL dampF(RVAL,rhoAB(ISTATE),NCMM(ISTATE),
     2             MMLR1D,IDF(ISTATE),IDSTT(ISTATE),KDER,Dm,Dmp,Dmpp)
                  DO  m= 1,NCMM(ISTATE)
                      IF(rhoAB(ISTATE).LE.0.d0) THEN
                          dULRdCm(m)= 1.d0/RVAL**MMLR(m,ISTATE)
                        ELSE
                          dULRdCm(m)= Dm(m)/RVAL**MMLR(m,ISTATE)
                        ENDIF
                      ULR= ULR + CmVAL(m,ISTATE)*dULRdCm(m)
                      ENDDO
                ENDIF
              XTEMP= XTEMP*ULR/ULRe
c... note ... reference energy for each state is asymptote ...
              DVDD= XTEMP*(XTEMP - 2.D0)  
              VPOT(I,ISTATE)= DE(ISTATE)*DVDD + VLIM(ISTATE)
              IF(RDIST.GT.0.d0) THEN
                  VDIST= VPOT(I,ISTATE)
                  BETADIST= VAL
                  IF(IDAT.LE.0) GO TO 999
                  ENDIF
              YPP= 2.d0*DE(ISTATE)*(1.0d0-XTEMP)*XTEMP
              IPV= IPVSTART+2
              IF((NCMM(ISTATE).GE.4).AND.(MMLR(2,ISTATE).EQ.0)) THEN
c ... derivatives w.r.t. long-range parameters for Aubert-Frecon  uLR
                  IPV= IPV+1
                  DVtot(IPV,I)= YPP*((1.d0 - YP*YPE)*dLULRedCm(1)
     1                                               - dULRdCm(1)/ULR)
                  IPV= IPV+1
c... note that derivative w.r.t. atomic level splitting FIXED
                  DVtot(IPV,I)= 0.d0
                  IPV= IPV+1
                  DVtot(IPV,I)= YPP*((1.d0 - YP*YPE)*dLULRedCm(3)
     1                                               - dULRdCm(3)/ULR)
                  IPV= IPV+1
                  DVtot(IPV,I)= YPP*((1.d0 - YP*YPE)*dLULRedCm(4)
     1                                               - dULRdCm(4)/ULR)
c... and then derivatives w.r.t. Aubert-Frecon 2x2 C8 values ...
                  IPV= IPV+1
                  DVtot(IPV,I)= YPP*((1.d0 - YP*YPE)*dLULRedCm(5)
     1                                               - dULRdCm(5)/ULR)
                  IPV= IPV+1
                  DVtot(IPV,I)= YPP*((1.d0 - YP*YPE)*dLULRedCm(6)
     1                                               - dULRdCm(6)/ULR)
                ELSEIF((NCMM(ISTATE).GE.4).AND.(MMLR(2,ISTATE).EQ.-1))
     1                                                            THEN
c ... derivatives w.r.t. long-range prameters for 3x3 A-F case uLR
                  IPV= IPV+1
                  DVtot(IPV,I)= YPP*((1.d0 - YP*YPE)*dLULRedCm(1)
     1                                             - dULRdCm(1)/ULR)

c                   DVtot(IPV,I)= 2*DE(ISTATE)*((1-
c     1                ULR*DEXP(BETADIST*YP)/ULRe)**2)*(-
c     2                (dULRdCm(1))*DEXP(BETADIST*YP)/ULRe + 
c     3               ULR*DEXP(BETADIST*YP)*dLULRedCm(1)/ULRe**2) 

                  IPV= IPV+1
c ... note that derivative w.r.t. atomic level splitting FIXED = 0.0`
                  DVtot(IPV,I)= 0.d0
                  IPV= IPV+1
                  DVtot(IPV,I)= YPP*((1.d0 - YP*YPE)*dLULRedCm(3)
     1                                               - dULRdCm(3)/ULR)
                  IPV= IPV+1
                  DVtot(IPV,I)= YPP*((1.d0 - YP*YPE)*dLULRedCm(4)
     1                                               - dULRdCm(4)/ULR)
c ... and then derivatives w.r.t. Aubert-Frecon 3x3 C8 values ...
                  IPV= IPV+1
                  DVtot(IPV,I)= YPP*((1.d0 - YP*YPE)*dLULRedCm(5)
     1                                               - dULRdCm(5)/ULR)
                  IPV= IPV+1
                  DVtot(IPV,I)= YPP*((1.d0 - YP*YPE)*dLULRedCm(6)
     1                                               - dULRdCm(6)/ULR)

                ELSE
c ... derivative w.r.t. Cm's for ordinary MLR/MLJ case ...
                  DO  m= 1, NCMM(ISTATE)
                      IPV= IPV+ 1
                      DVtot(IPV,I)= -YPP*(dLULRedCm(m)*(YP*YPE- 1.d0)
     1                                               + dULRdCm(m)/ULR)
                      ENDDO
                ENDIF
c... derivative w.r.t. Re  
              DVtot(IPVSTART+2,I)= YPP*(YPE*DBDRe(I,ISTATE)
     1                                       + VAL*DYPDRE + dLULRedRe)
              IF(APSE(ISTATE).LE.0) THEN
c... derivative w.r.t. De  for 'conventional' power-series exponent
                  DVDD= DVDD + YPP*YP*YPE/DE(ISTATE)
                  IF((NCMM(ISTATE).GE.4).AND.(MMLR(2,ISTATE).EQ.0))
c... derivative w.r.t. De  for Aubert-Frecon 2x2 exponent
     1        DVDD= DVDD+ YPP*((1.d0- YP*YPE)*dLULRedDe - dULRdDe/ULR)
c... final value of derivative w.r.t. De [ignoring beta(0)]
                  DVtot(IPVSTART+1,I)= DVDD
                  YPP= YPP*YPE*(1.d0 - YP)

c???????  RJL ... check this out!  ????
                  IF((IDSTT(ISTATE).GT.1).AND.(IDF(ISTATE).EQ.-1))
     1                                            YPP= YPP*(1.d0 + YP)
c... finally ... derivatives w.r.t. exponent expansion coefficients
                  DO  J= 0, Nbeta(ISTATE)
                      IPV= IPV+1
                      DVtot(IPV,I)= YPP
                      YPP= YPP*YQ
                      ENDDO
                ELSE
c... for Pashov-spline exponent cases...
                  YPP= YPP*YPE
                  DO  J= 1,Nbeta(ISTATE)
                      IPV= IPV+ 1
                      DVtot(IPV,I)= DBDB(J,I,ISTATE)*YPP
                      ENDDO
                  DVtot(IPVSTART+1,I)= DVDD
     1                             + YPP*DBDB(npow,I,ISTATE)/DE(ISTATE)
                ENDIF
              ENDDO
ccccc Print for testing
      rewind(10)
      write(10,610) (RD(i,ISTATE),vpot(i,istate),BETAFX(i,istate),
     1                              i= 1, NDATPT(ISTATE),OSEL(ISTATE))
  610 FORMAT(/(f10.4,f15.5,f12.6))
ccccc End of Print for testing
          ENDIF
c********* End for Morse/Lennard-Jones(p) potential function ***********

      IF(PSEL(ISTATE).EQ.3) THEN
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** For the Double-Exponential Long-Range Oscillator.
          AA(ISTATE)= 0.0d0
          BB(ISTATE)= 0.0d0
          dAAdRe= 0.0d0
          dBBdRe= 0.0d0
          dVdBtemp= 0.0d0
          KDER=2
c ... save uLR powers in a 1D array
          DO  m= 1, NCMM(ISTATE)
              MMLR1D(m)= MMLR(m,ISTATE)
              ENDDO
          CALL dampF(RE(ISTATE),rhoAB(ISTATE),NCMM(ISTATE),
     1              MMLR1D,IDF(ISTATE),IDSTT(ISTATE),KDER,Dm,Dmp,Dmpp)
c... first, get  AA & BB and their derivatives!
          DO  m= 1,NCMM(ISTATE)
              MMLRP= MMLR(m,ISTATE)
              CmVALL= CmVAL(m,ISTATE)
              AA(ISTATE)= AA(ISTATE) - CmVALL/RE(ISTATE)**MMLRP
     2     *(Dm(m) + (Dmp(m) - MMLRP*Dm(m)/RE(ISTATE))/BETA(0,ISTATE))
              BB(ISTATE)= BB(ISTATE) - CmVALL*Dm(m)/RE(ISTATE)**MMLRP
              dAAdRe= dAAdRe - CmVALL/RE(ISTATE)**MMLRP * (Dmp(m) 
     1        - MMLRP*Dm(m)/RE(ISTATE) - (Dmpp(m) - MMLRP*(2.d0*Dmp(m)
     2             - (MMLRP+1)/RE(ISTATE))/RE(ISTATE))/BETA(0,ISTATE))
              dBBdRe= dBBdRe - CmVALL/RE(ISTATE)**MMLRP 
     1                              *(Dmp(m) - MMLRP*Dm(m)/RE(ISTATE))
c ...  dVdBtemp  is   d{AA}/d{beta(0)}
              dVdBtemp= dVdBtemp - CmVALL*(Dm(m)*MMLRP/RE(ISTATE)
     1                 - Dmp(m))/(RE(ISTATE)**MMLRP*BETA(0,ISTATE)**2)
              ENDDO
          dBBdRe= dBBdRe + dAAdRe
          AA(ISTATE)= AA(ISTATE) + DE(ISTATE) 
          BB(ISTATE)= AA(ISTATE) + BB(ISTATE) + DE(ISTATE)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          DO  I= ISTART,ISTOP
c** First, calculate the generalized Morse exponent
              RVAL= RD(I,ISTATE)
              IF(RDIST.GT.0.d0) RVAL= RDIST
              RDp= RVAL**NPB(ISTATE)
              YP= (RDp-REP)/(RDp+REP)
              npow= Nbeta(ISTATE)
              BETAFX(I,ISTATE)= BETA(0,ISTATE)
              DO  J= 1, Nbeta(ISTATE)
                  BETAFX(I,ISTATE)= BETAFX(I,ISTATE) +
     1                                         BETA(J,ISTATE)*YP**J
                  ENDDO
c** Calculate some temporary variables.
              XTEMP= DEXP(-BETAFX(I,ISTATE)*(RVAL-RE(ISTATE)))
              Xtemp2= XTEMP*XTEMP
              BMtemp= BETAFX(I,ISTATE)
              BMtemp2= BETAFX(I,ISTATE)
c
c** Now to calculate the actual potential and partial derivatives:
c
              ULR= 0.0d0
              CALL dampF(RVAL,rhoAB(ISTATE),NCMM(ISTATE),
     1              MMLR1D,IDF(ISTATE),IDSTT(ISTATE),KDER,Dm,Dmp,Dmpp)
                  DO  m= 1,NCMM(ISTATE)
                  MMLRP= MMLR(m,ISTATE)
                  CmVALL= CmVAL(m,ISTATE)
                  ULR= ULR + CmVALL*Dm(m)/RVAL**MMLRP
                  ENDDO
              VPOT(I,ISTATE)= (AA(ISTATE)*XTEMP - BB(ISTATE))*XTEMP 
     1                                            - ULR + VLIM(ISTATE)
              IF(RDIST.GT.0.d0) THEN
                  VDIST= VPOT(I,ISTATE)
                  BETADIST= BETAFX(I,ISTATE)
                  IF(IDAT.LE.0) GO TO 999
                  ENDIF
c
c** Now to calculate the partial derivatives
c
              IPV= IPVSTART + 1
c ... first, derivative of the potential w.r.t. De
              DVtot(IPV,I)= (XTEMP - 2.0d0)*XTEMP
              DVDD= (XTEMP - 2.0d0)*XTEMP 
              BTEMP= 0.0d0
              Btemp2= 0.0d0
              DO  J= 1,Nbeta(ISTATE)
                  Btemp2= Btemp2 + BETA(J,ISTATE)*DBLE(J)*YP**(J-1)
                  ENDDO
              Btemp2= Btemp2*(RVAL-RE(ISTATE))
              BTEMP= BTEMP*(RVAL-RE(ISTATE))
c
c** Now to calculate the remaining partial derivatives.
c
              PETEMP= 2.0d0*nPB(ISTATE)*RE(ISTATE)**(nPB(ISTATE)-1)
     &                                          * RVAL**nPB(ISTATE) /
     &                 (RVAL**nPB(ISTATE)+RE(ISTATE)**nPB(ISTATE))**2
              PBTEMP= BB(ISTATE) * XTEMP * (RVAL-RE(ISTATE))
              PBtemp2= 2.0d0 * AA(ISTATE) * Xtemp2* (RVAL-RE(ISTATE))
              IPV= IPV + 1
c ... then, derivative of the potential w.r.t. Re
              DVtot(IPV,I)= 2.0d0*AA(ISTATE)*Xtemp2*
     1                       (BMtemp2+PETEMP* Btemp2) + Xtemp2*dAAdRe
     2        - BB(ISTATE)*XTEMP*(BMtemp+PETEMP*BTEMP) - XTEMP*dBBdRe
c ... finally, derivatives of the potential w.r.t. the \beta_i
              IPV= IPV+ 1
              DVtot(IPV,I)= -PBtemp2 + PBTEMP + dVdBtemp*(Xtemp2-XTEMP)
              DO  J= 1,APSE(ISTATE)
                  IPV= IPV+ 1
                  DVtot(IPV,I)= (-PBtemp2 + PBTEMP) * YP**J
                  ENDDO
c?????? PROBLEM here !!!!!!!!!!!!!!
c?????? why is deriv different for  R > Re vs. R < Re ???  09/-02/13 ??/
              DO  J= APSE(ISTATE)+1, Nbeta(ISTATE)
                  IPV= IPV+ 1
                  IF (RVAL .LT. RE(ISTATE)) THEN
                      DVtot(IPV,I)= 0.0d0
                    ELSE
                      DVtot(IPV,I)= PBTEMP * YP**J
                    ENDIF
                  ENDDO
c *** DBDRe and DBDB is used in uncertainty calculation, see fununc.f
              DBDRe(I,ISTATE)= -PETEMP*BTEMP/(RVAL-RE(ISTATE))
              DBDB(0,I,ISTATE)= 1.0d0
              DO  J= 1, Nbeta(ISTATE)
                  DBDB(J,I,ISTATE)= DBDB(J-1,I,ISTATE)*YP
                  ENDDO
              ENDDO
          ENDIF
******7* End Double-Exponential Long-Range Potential Function ********

      IF(PSEL(ISTATE).EQ.4) THEN
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** For the Tiemann Polynomial Potential Energy Function.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          BT= BETA(Nbeta(ISTATE)+1, ISTATE)
          Rinn= BETA(Nbeta(ISTATE)+2, ISTATE)
          Rout= BETA(Nbeta(ISTATE)+3, ISTATE)
c** Determine analytic function attaching smoothly to inner wall of 
c  polynomial expansion at  R= Rinn < Rm
          YP= (Rinn - RE(ISTATE))/(Rinn+ BT*RE(ISTATE))
          YPP= 1.d0
          A1= BETA(0,ISTATE)
          A2= 0.d0
          DO  J= 1, Nbeta(ISTATE)
              A2= A2+ J*YPP*BETA(J,ISTATE)
              YPP= YPP*YP
              A1= A1+ YPP*BETA(J,ISTATE)
              ENDDO
          A2= A2*(RE(ISTATE)+ BT*RE(ISTATE))/(Rinn + BT*RE(ISTATE))**2
c* If inward extrapolation is exponential:   A1*exp(-A2*(R-Rinn))
          A2= -A2/A1
c** With long-range tail an NCMM-term inverse-power sum, add 1 exponential
c   term to ensure continuity and smoothness at  Rout
          YP= (Rout - RE(ISTATE))/(Rout+ BT*RE(ISTATE))
          YPP= 1.d0
          B1= BETA(0,ISTATE) + VLIM(ISTATE) - DE(ISTATE)
          B2= 0.d0
          DO  J= 1, Nbeta(ISTATE)
              B2= B2+ J*YPP*BETA(J,ISTATE)
              YPP= YPP*YP
              B1= B1+ YPP*BETA(J,ISTATE)
              ENDDO
          B2= B2*(RE(ISTATE)+ BT*RE(ISTATE))/(Rout + BT*RE(ISTATE))**2
          B3= VLIM(ISTATE)
          B4= 0.d0
          DO  J= 1, NCMM(ISTATE)
              B5= CmVAL(J,ISTATE)/Rout**MMLR(J,ISTATE)
              B3= B3- B5
              B4= B4+ MMLR(J,ISTATE)*B5
              ENDDO
          B3= B3- B1
          B4= (B2- B4/Rout)/B3
c ... now generate potential as a Tiemann-type expansion
c??
c               WRITE(6,698) IPVSTART
c698        FORMAT('IPV before Tiemann poly. loop:'I5)
c??
          DO  I= ISTART, ISTOP
              RVAL= RD(I,ISTATE)
              IF(RDIST.GT.0) RVAL= RDIST
              IF(RVAL.LE.Rinn) THEN
c ... for exponential inward extrapolation ...
                  VPOT(I,ISTATE)= A1*DEXP(-A2*(RVAL- Rinn))
     1                                     + VLIM(ISTATE) - De(ISTATE)
                ELSEIF(RVAL.LE.Rout) THEN
                  YP= (RVAL - RE(ISTATE))/(RVAL + BT*RE(ISTATE))
                  A3= BETA(0,ISTATE) + VLIM(ISTATE) - DE(ISTATE)
                  YPP= 1.d0
                  DO  J= 1,Nbeta(ISTATE)
                      YPP= YPP*YP
                      A3= A3+ BETA(J,ISTATE)*YPP
                      ENDDO
                  VPOT(I,ISTATE)= A3 
                ELSEIF(RVAL.GT.Rout) THEN
                  A3= VLIM(ISTATE)
                  DO  J= 1, NCMM(ISTATE)
                      A3= A3- CmVAL(J,ISTATE)/RVAL**MMLR(J,ISTATE)
                      ENDDO
                  VPOT(I,ISTATE)= A3 - B3*DEXP(-B4*(RVAL- Rout))
cc     1                                 + VLIM(ISTATE) - De(ISTATE)
                ENDIF
              IF(RDIST.GT.0) THEN
                  VDIST= VPOT(I,ISTATE)
                  BETADIST= 0.d0 
                  IF(IDAT.LE.0) GO TO 999
                  ENDIF
              IPV= IPVSTART+ 1
              DVtot(IPV,I)= -1.d0
              IPV= IPV+1
              DVtot(IPV,I)= 0.d0
              DO  J=1, NCMM(ISTATE)
                  IPV= IPV+1
                  DVtot(IPV,I)= 0.d0
                  ENDDO 
              IPV= IPV+ 1
              DVtot(IPV,I)= 1.d0
              DO  J= 1, Nbeta(ISTATE)
                 IPV= IPV+ 1
                 DVtot(IPV,I)= DVtot(IPV-1,I)*YP
                 ENDDO
              ENDDO
c???
c               WRITE(6,699) IPV
c699        FORMAT('IPV after Tiemann poly. loop:'I5)
c???
      rewind(10)
      write(10,612) (RD(i,ISTATE),vpot(i,istate),
     1                                        i= 1, NDATPT(ISTATE),50)
          ENDIF
c*********** End Tiemann Potential Energy Function *****************

      IF(PSEL(ISTATE).EQ.5) THEN
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** For the TANG-Toennies-type potential with simple exponential term 
c   and damped (s=+1) repulsion terms:  A= A(b,r_e,{C_m})
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c ... first ... save uLR powers in a 1D array
            DO  m= 1, NCMM(ISTATE)
                MMLR1D(m)= MMLR(m,ISTATE)
                ENDDO
c** Now, calculate A and De  using the input values of b and r_e 
            REQ= RE(ISTATE)
            BTT= BETA(0,ISTATE)
            rhoAB(ISTATE)= BTT/3.13d0
            KDER= 1
            CALL dampF(REQ,rhoAB(ISTATE),NCMM(ISTATE),
     1              MMLR1D,IDF(ISTATE),IDSTT(ISTATE),KDER,Dm,Dmp,Dmpp)
            VATT= 0.d0
            DTT= 0.d0
            DO  m= 1,NCMM(ISTATE)
                TM= CMval(m,ISTATE)/REQ**MMLR1D(m)
                VATT= VATT + Dm(m)*TM
                DTT= DTT+TM*(MMLR1D(m)*Dm(MMLR1D(m))/REQ-Dmp(MMLR1D(m)))
                ENDDO
c*** If getting BTT and AII from D_e abd r_e, ...
c            BTT= DTT/(VATT - DE(ISTATE))
c            ATT= (VATT - DE(ISTATE))*DEXP(BTT*RE(ISTATE))
c** if getting ATT and D_e from R_e and b, ...
            ATT= -DTT*DEXP(-BTT*REQ)/BTT
            DE(ISTATE)= VATT - ATT*DEXP(-BTT*REQ)
c
            WRITE(6,600) SLABL(ISTATE),DE(ISTATE),ATT
  600 FORMAT('  Tang-Toennies potential for state-',A2,'    has  D_e=',
     1  f11.6,'  and A_{TT}=',F13.2)
            DO  I= ISTART, ISTOP
                RVAL= RD(I,ISTATE)
                CALL dampF(RVAL,rhoAB(ISTATE),NCMM(ISTATE),
     1              MMLR1D,IDF(ISTATE),IDSTT(ISTATE),KDER,Dm,Dmp,Dmpp)
                VATT= 0.d0
                DO M= 1,NCMM(ISTATE)
                    VATT= VATT+ CmVAL(m,ISTATE)*Dm(m)/RVAL**MMLR1D(m)
                    ENDDO
                VPOT(I,ISTATE)= ATT*EXP(-BTT*RVAL)- VATT
              ENDDO
            ENDIF
      IF(PSEL(ISTATE).EQ.6) THEN
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** For the Surkus-type Generalized Potential Energy Function.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** First, we calculate the implied Dissociation Energy (if it exists)
        IF(AGPEF(ISTATE).NE.0.0d0) THEN
            YPP= 1.d0/AGPEF(ISTATE)**2
            VAL= YPP
            DO  I= 1, Nbeta(ISTATE)
                YPP= YPP/AGPEF(ISTATE)
                VAL= VAL + BETA(I,ISTATE)*YPP
                ENDDO
            DE(ISTATE)= VAL*BETA(0,ISTATE)
            ENDIF
        DO  I= ISTART, ISTOP
            RVAL= RD(I,ISTATE)
            IF(RDIST.GT.0.d0) RVAL= RDIST
            RDp= RVAL**nPB(ISTATE)
            YP= (RDp-REP)/(AGPEF(ISTATE)*RDp + BGPEF(ISTATE)*REP)
c** Now to calculate the actual potential
            YPP= 1.d0
            VAL= 1.d0
            DVAL= 2.d0
            DO  J= 1, Nbeta(ISTATE)
                YPP= YPP*YP
                VAL= VAL + BETA(J,ISTATE)*YPP
                DVAL= DVAL+ (J+2)*BETA(J,ISTATE)*YPP
                ENDDO
            VPOT(I,ISTATE)= VAL*BETA(0,ISTATE)*YP**2 + VLIM(ISTATE) 
            IF(RDIST.GT.0) THEN
                VDIST= VPOT(I,ISTATE)
                BETADIST= 0.d0
                IF(IDAT.LE.0) GO TO 999
                ENDIF
            DVAL= DVAL*BETA(0,ISTATE)*YP
c** Now to calculate the partial derivatives
            DVDD= 0.d0
            IPV= IPVSTART + 1
c ... derivative of the potential w.r.t. Re
            DVtot(IPV,I)= -DVAL*REP*RDp*(AGPEF(ISTATE)+BGPEF(ISTATE))
     1                  *(nPB(ISTATE)/RE(ISTATE))/(AGPEF(ISTATE)*RDp +
     2                                           BGPEF(ISTATE)*REP)**2
c ... and derivatives w.r.t. the beta_i expansion coefficients ...
            IPV= IPV+ 1
            DVtot(IPV,I)= VAL*YP**2
            IPV= IPV+ 1
            DVtot(IPV,I)= BETA(0,ISTATE)*YP**3
            DO  J= 2, Nbeta(ISTATE)
                IPV= IPV+ 1
                DVtot(IPV,I)= DVtot(IPV-1,I)*YP
                ENDDO
            ENDDO
        IF(RDIST.LE.0.d0) VLIM(ISTATE)= VPOT(NDATPT(ISTATE),ISTATE)
c????
        rewind(10)
        write(10,612) (RD(i,ISTATE),vpot(i,istate),
     1                                        i= 1, NDATPT(ISTATE),50)
  612 FORMAT(/(f10.4,f15.5))
c????
        ENDIF
c*********** End Generalized Potential Energy Function *****************

      IF((NUA(ISTATE).GE.0).OR.(NUB(ISTATE).GT.0)) THEN
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Treat any 'adiabatic' BOB radial potential functions here ...
c        u_A(r) = yp*uA_\infty + [1 - yp]\sum_{i=0,NUA} {uA_i yq^i} 
c   where the  u_\infty  values stored/fitted as  UA(NUA(ISTATE))
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          HReP= 0.5d0*pAD(ISTATE)/RE(ISTATE)
          HReQ= 0.5d0*qAD(ISTATE)/RE(ISTATE)
          REadAp= RE(ISTATE)**pAD(ISTATE)
          REadAq= RE(ISTATE)**qAD(ISTATE)
          REadBp= RE(ISTATE)**pAD(ISTATE)
          REadBq= RE(ISTATE)**qAD(ISTATE)
          IF((BOBCN(ISTATE).GE.1).AND.(pAD(ISTATE).EQ.0)) THEN
              HReP= 2.d0*HReP
              HReQ= 2.d0*HReQ
              ENDIF
c ... reset parameter counter ...
          IPVSTART= IPV
          DO  I= ISTART,ISTOP
              RVAL= RD(I,ISTATE)
              IF(RDIST.GT.0.d0) RVAL= RDIST
              RDp= RVAL**pAD(ISTATE)
              RDq= RVAL**qAD(ISTATE)
              YPA= (RDp - REadAp)/(RDp + REadAp)
              YQA= (RDq - REadAq)/(RDq + REadAq)
              YPB= (RDp - REadBp)/(RDp + REadBp)
              YQB= (RDq - REadBq)/(RDq + REadBq)
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
                  YQP= 1.d0
                  IF(NUA(ISTATE).GE.2) THEN
                      DO  J= 1,NUA(ISTATE)-1
                          DVAL= DVAL+ DBLE(J)*YQP*UA(J,ISTATE)
                          YQP= YQP*YQA
                          VAL= VAL+ UA(J,ISTATE)*YQP
                          IPV= IPV+ 1
                          DVtot(IPV,I)= YPMA*YQP
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
                  YQP= 1.d0
                  IF(NUB(ISTATE).GE.2) THEN
                      DO  J= 1,NUB(ISTATE)-1
                          DVAL= DVAL+ DBLE(J)*YQP*UB(J,ISTATE)
                          YQP= YQP*YQB
                          VAL= VAL+ UB(J,ISTATE)*YQP
                          IPV= IPV + 1
                          DVtot(IPV,I)= YPMB*YQP
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
          HReP= 0.5d0*pNA(ISTATE)/RE(ISTATE)
          HReQ= 0.5d0*qNA(ISTATE)/RE(ISTATE)
          REnaAp= Re(ISTATE)**pNA(ISTATE)
          REnaAq= Re(ISTATE)**qNA(ISTATE)
          REnaBp= Re(ISTATE)**pNA(ISTATE)
          REnaBq= Re(ISTATE)**qNA(ISTATE)
          IF((BOBCN(ISTATE).GE.1).AND.(pNA(ISTATE).EQ.0)) THEN
              HReP= 2.d0*HReP
              HReQ= 2.d0*HReQ
              ENDIF
          IPVSTART= IPV
          DO  I= ISTART,ISTOP
              RVAL= RD(I,ISTATE)
              IF(RDIST.GT.0.d0) RVAL= RDIST
              RM2= 1/RVAL**2
              RDp= RVAL**pNA(ISTATE)
              RDq= RVAL**qNA(ISTATE)
              YPA= (RDp - REnaAp)/(RDp + REnaAp)
              YQA= (RDq - REnaAq)/(RDq + REnaAq)
              YPB= (RDp - REnaBp)/(RDp + REnaBp)
              YQB= (RDq - REnaBq)/(RDq + REnaBq)
              YPMA= 1.d0 - YPA
              YPMB= 1.d0 - YPB
              IF(BOBCN(ISTATE).GE.1) THEN
c** If  BOBCN > 0  &  p= 1,  assume use of Ogilvie-Tipping vble.
                  YPMA= 1.d0
                  YPA= 2.d0*YPA
                  ENDIF
              IF(NTA(ISTATE).GE.0) THEN 
c ... Now ... derivatives of TA w,r,t, expansion coefficients
                  VAL= TA(0,ISTATE)
                  DVAL= 0.d0
                  IPV= IPVSTART + 1
                  DVtot(IPV,I)= YPMA*RM2
                  YQP= 1.d0
                  IF(NTA(ISTATE).GE.2) THEN
                      DO  J= 1,NTA(ISTATE)-1
                          DVAL= DVAL+ DBLE(J)*YQP*TA(J,ISTATE)
                          YQP= YQP*YQA
                          VAL= VAL+ TA(J,ISTATE)*YQP
                          IPV= IPV + 1
                          DVtot(IPV,I)= YPMA*YQP*RM2
                          ENDDO
                      ENDIF
                  IPV= IPV + 1
                  DVtot(IPV,I)= YPA*RM2
                  TAR(I,ISTATE)= VAL*YPMA + YPA*TA(NTA(ISTATE),ISTATE)
c ... and derivative of TA w.r.t. Re ... 
                  DTADRe(I,ISTATE)= (-HReQ*(1.d0 - YQA**2)*YPMA*DVAL
     1       + HReP*(1.d0 - YPA**2)*(VAL- TA(NTA(ISTATE),ISTATE)))*RM2
                  ENDIF
              IF(NTB(ISTATE).GE.0) THEN
c ... Now ... derivatives of TB w.r.t. expansion coefficients
                  VAL= TB(0,ISTATE)
                  DVAL= 0.d0
                  IF(NTA(ISTATE).LT.0) THEN
                      IPV= IPVSTART + 1
                    ELSE
                      IPV= IPV + 1
                    ENDIF
                  DVtot(IPV,I)= YPMB*RM2
                  YQP= 1.d0
                  IF(NTB(ISTATE).GE.2) THEN
                      DO  J= 1,NTB(ISTATE)-1
                          DVAL= DVAL+ DBLE(J)*YQP*TB(J,ISTATE)
                          YQP= YQP*YQB 
                          VAL= VAL+ TB(J,ISTATE)*YQP
                          IPV= IPV + 1
                          DVtot(IPV,I)= YPMB*YQP*RM2
                          ENDDO
                      ENDIF
                  IPV= IPV + 1
                  DVtot(IPV,I)= YPB*RM2
                  TBR(I,ISTATE)= VAL*YPMB + YPB*TB(NTB(ISTATE),ISTATE)
c ... and derivative of TA w.r.t. Re ...
                  DTBDRe(I,ISTATE)= (-HReQ*(1.d0 - YQB**2)*YPMB*DVAL
     1       + HReP*(1.d0 - YPB**2)*(VAL- TB(NTB(ISTATE),ISTATE)))*RM2
                  ENDIF
              ENDDO
          ENDIF
c++++ END of treatment of non-adiabatic centrifugal BOB function++++++++

      IF(NwCFT(ISTATE).GE.0) THEN
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Treat any Lambda- or 2\Sigma-doubling radial strength functions here
c    representing it as   f(r)= Sum{ w_i * y_{Pqw}^i}
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          LAMB2= 2*IOMEG(ISTATE)
          HReP= 0.5d0*Pqw(ISTATE)/RE(ISTATE)
          REwp= RE(ISTATE)**Pqw(ISTATE)
          IPVSTART= IPV
          DO  I= ISTART,ISTOP
              RVAL= RD(I,ISTATE)
              IF(RDIST.GT.0.d0) RVAL= RDIST
              RMF= 1.d0/RVAL**2
              IF(IOMEG(ISTATE).GT.0)  RMF= RMF**LAMB2
              RDp= RVAL**Pqw(ISTATE)
              YP= (RDp - REwp)/(RDp + REwp)
              DVAL= 0.d0
              YQP= RMF
              VAL= wCFT(0,ISTATE)*YQP
              IPV= IPVSTART + 1
              DVtot(IPV,I)= YQP
              IF(NwCFT(ISTATE).GE.1) THEN
                  DO  J= 1,NwCFT(ISTATE)
                      DVAL= DVAL+ DBLE(J)*YQP*wCFT(J,ISTATE)
                      YQP= YQP*YP
                      IPV= IPV + 1
                      DVtot(IPV,I)= YQP
                      VAL= VAL+ wCFT(J,ISTATE)*YQP
                      ENDDO
                  ENDIF
              wRAD(I,ISTATE)= VAL 
              dLDDRe(I,NSTATEMX)= -HReP*(1.d0 - YP**2)*DVAL
              ENDDO
          ENDIF
c++++ END of treatment of Lambda/2-Sigma centrifugal BOB function+++++++

c++++ Test for inner wall inflection , and if it occurs, replace inward
c++++ potential with linear approximation +++++++
      I1= (RE(ISTATE)-RD(1,ISTATE))/(RD(2,ISTATE)-RD(1,ISTATE))
      IF(I1.GT.3) THEN
          SL= 0.d0
          DO  I= I1-2, 1, -1
              SLB= SL
              SL= VPOT(I,ISTATE) - VPOT(I+1,ISTATE)
              IF(SL.LE.SLB) THEN
                  DO  J= I,1,-1
                      VPOT(J,ISTATE)= VPOT(I,ISTATE) + (I-J)*SL
                      ENDDO
                  WRITE(6,606) SLABL(ISTATE),RD(I,ISTATE),VPOT(I,ISTATE)
                  GOTO 66
                  ENDIF
              ENDDO
          ENDIF
   66 CONTINUE
  606 FORMAT(9('===')/'!*!* Extrapolate to correct ',A2,' inner-wall inf
     1lection at   R=',f6.4,'   V=',f8.0/9('==='))
c++++++++++++End of Inner Wall Test/Correction code+++++++++++++++++++++
c======================================================================
c** At the one distance RDIST calculate total effective potential VDIST
c  including (!!) centrifugal and Lambda/2Sigma doubling terms,
c  and its partial derivatives w.r.t. Hamiltonian parameters dVdPk.
c** This case only for simulation & fitting of tunneling width data.
c
      IF((RDIST.GT.0).AND.(IDAT.GT.0)) THEN 
          NBAND= IB(IDAT)
          IISTP= ISTP(NBAND)
cccccccc
c         WRITE (40,644) IISTP,RDIST,RVAL,VDIST,I,NDATPT(ISTATE)
c 644 FORMAT ('IISTP =',I3,' RDIST =',G16.8,' RVAL =',G16.8,
c    &          ' VDIST =',G16.8,' I =',I6,' NDATPT =',I6)
cccccccc
          BFCT= 16.857629205d0/(ZMASS(3,IISTP)*RDIST**2)
          JFCT= DBLE(JPP(IDAT)*(JPP(IDAT)+1))
          IF(IOMEG(ISTATE).GT.0) JFCT= JFCT - IOMEG(ISTATE)**2
          IF(IOMEG(ISTATE).EQ.-2) JFCT= JFCT + 2.D0
          JFCT= JFCT*BFCT
c ... First get total effective potential, including BOB terms
          VDIST= VDIST + JFCT
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
              VDIST= VDIST + JFCTLD* WRAD(ISTOP,ISTATE)
             ENDIF
cccccccc
c       WRITE (40,648) JPP(IDAT),EFPP(IDAT),RDIST,VDIST
c 648   FORMAT ('J =',I3,' efPARITY =',I3,' RDIST =',G16.8,' VDIST =',
c    1                                                           G16.8/)
cccccccc
          IF(PSEL(ISTATE).GT.0) THEN
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
          ENDIF
c*****7********************** BLOCK END ******************************72
  999 RETURN
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
