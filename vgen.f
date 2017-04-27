c***********************************************************************
      SUBROUTINE VGEN(ISTATE,RDIST,VDIST,BETADIST,IDAT)
c***********************************************************************
c** This subroutine will generate one of seven possible families of
c   analytical molecular potentials for the program dPotFit,as specified
c    by parameter PSEL
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++  COPYRIGHT 2006-2016  by  R.J. Le Roy, Jenning Seto, Yiye Huang  ++
c  and N. Dattani, Department of Chemistry, University of of Waterloo, +
c                        Waterloo, Ontario, Canada                     +
c    This software may not be sold or any other commercial use made    +
c      of it without the express written permission of the authors.    +
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c               ----- Version of  10 April 2017 - with RREFq/p -----
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** On entry:
c   ISTATE  is the electronic state being considered in this CALL.
c   RDIST: If RDIST > 0, calculate potl & derivs only @ that onee distance
c   -----    * return potential function at that point as VDIST and
c              potential function exponent as BETADIST
c        * skip partial derivative calculation if  IDAT.LT.0 {pointless??}
c        * If RDIST.le.0  calculate PEC & BOB fx. & partial derivatives 
c             at distances given by array RD(i,ISTATE) & return them
c        ->  RDIST > 0 & IDAT > 0 for tunneling width & derivative calc'n
c        ->  RDIST > 0 & IDAT.le.0 for PEC only w. IDAT<0 to omit derivs
c** On entry via common blocks: 
c* PSEL specifies how data for state ISTATE are to be represented!
c     PSEL = -2 : Represent {v,J,p} levels of state ISTATE by term values
c     PSEL = -1 : Represent {v,J} levels with band constants
c     PSEL =  0 : Use a fixed potential function defined in READPOT.
c     PSEL =  1 : Use an Expanded Morse Oscillator(p) potential.
c     PSEL =  2 : Use a Morse/Lennard-Jones(p) potential.
c     PSEL =  3 : Use a Double-Exponential Long-Range Potential.
c     PSEL =  4 : Use a Surkus "Generalized Potential Energy Function".
c     PSEL =  5 : Use a [Tiemann Polynomial]
c     PSEL =  6 : Use a Generalized Tang-Toennies-type potential
c     PSEL =  7 : Use an Aziz'ian HFD-ABC or HFD-D potential
c* Nbeta(s) is order of the beta(r) {exponent} polynomial or # spline points
c  APSE(s).le.0  to use PE-MLR {p,q}-type exponent polynomial of order Nbeta(s)
c   if APSE(s) > 0 \beta(r) is SE-MLR spline defined by Nbeta(s) points
c          \beta_i=PARM(i) at distances defined by  y_q^{Rref} 
c*  MMLR(j,s)  are long-range inverse-powers for an MLR or DELR potential
c*  pPOT(s)  the basic value of power p for the beta(r)  exponent function
c*  qPOT(s)  the power q for the power series expansion variable in beta(r)
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
c** Common block for partial derivatives of potential at the one 
c   distance RDIST and HPP derivatives for uncertainties
      REAL*8 dVdPk(HPARMX),dDe(0:NbetaMX),dDedRe
      COMMON /dVdPkBLK/dVdPk,dDe,dDedRe
c=======================================================================
c** Define local variables ...
      INTEGER I,J,I1,ISTATE,IPV,IPVSTART,ISTART,ISTOP,LAMB2,m,m1,npow,
     1  POWmax,IDAT,IISTP,MMp2, NIFL,MCMM,MMLR1D(NCMMax)
      REAL*8 BTEMP,BINF,RVAL,RTEMP,RM2,XTEMP,PBTEMP,PETEMP,Btemp2,RMF,
     1 PBtemp2, bohr, Cm1D(NCMMax),CmEFF1D(NCMMax),
     2 C3VAL,C3bar,C3Pi,C6bar,C6adj,C6Pi,C8Pi,C9adj,C11adj,C8VAL,YP,YQ,
     3 YPA,YPB,YQA,YQB,YPE,YPM,YPMA,YPMB,YPP,YQP,REp,RDp,RDq,dYPdRRp,
     4 DYPDRE,VAL,BVAL,DVAL,HReP,HReQ,yqRe,dyqRedRe,betaRe,DbetaRe,yPOW,
     4 dAAdb0,dbetaFX,ULRe,dULRe,d2ULRe,SL,SLB,SLBB,AREFpp,AREFqq,T0,
     5 T1,T2,Scalc,dLULRedRe,dLULRedCm(NCMMax),dLULRedDe,dULRdDe,rhoINT,
     6 dULRdCm(NCMMax),dULRepdCm(NCMMax),dULRedCm(NCMMax),DVDD,RDIST,
     7 VDIST,BETADIST,X,BFCT,JFCT,JFCTLD,REadAp,REadBp,REadAq,REadBq,
     8 REnaAp,REnaBp,REnaAq,REnaBq,REwp,dC6dDe,dC9dC3,dC9dC6,dC9dDe,BT,
     9 Rinn,Rout,A0,A1,A2,A3,xBETA(NbetaMX),rKL(NbetaMX,NbetaMX),C1LIM,
     o B5,BETA0,BETAN,TM,VATT,dATTdRe,dATTdb,ATT,REq,VMIN,dYQdRRq,
     a REQQ,XRI,dXRI,fRO,XRIpw,XRO,dXRO,XROpw,ROmp2,dXROdRe,d2XROdRe,
     b DXRIdRe,d2XRIdRe,dCmp2dRe,EXPBI,BIrat,CMMp2,RMMp2,dAIdRe,dBIdRe,
     c VX,dVX,dVdRe,dDeROdRe,dDeRIdRe,dULRdR,dCmASUM, dCmBSUM,AREFq,
     d AREFp,  dAI(0:NbetaMX),dBI(0:NbetaMX),dCmp2(0:NbetaMX),
     e DEIGM1(1,1),DEIGM3(1,1),DEIGM5(1,1),DEIGR(1,1),DEIGRe(1,1),
     f DEIGDe(1,1),DEIGMx(NCMMax,1,1)
c***********************************************************************
c** Temporary variables for MLR and DELR potentials
      REAL*8 ULR,dAAdRe,dBBdRe,dVdBtemp,CmVALL, Dm(NCMMax),Dmp(NCMMax),
     1  Dmpp(NCMMax)
c***********************************************************************
cc    SAVE MCMM,Cm1D,MMLR1D      ??? not needed - regenerate every call
c** Initializing variables.
      DATA bohr/0.52917721092d0/     !! 2010 physical constants d:mohr12
      REq= RE(ISTATE)**qPOT(ISTATE)
      REp= RE(ISTATE)**pPOT(ISTATE)
      AREFq= RREFq(ISTATE)
      AREFp= RREFp(ISTATE)
      IF(AREFq.LE.0.d0) THEN
          AREFq= RE(ISTATE)
          IFXrefq(ISTATE)= 1
          ENDIF
      IF(AREFp.LE.0.d0) THEN
          AREFp= AREFq
          IFXrefp(ISTATE)= 1
          ENDIF
      AREFqq= AREFq**qPOT(ISTATE)
      AREFpp= AREFp**pPOT(ISTATE)
c** Normally data point starts from 1
      ISTART= 1
      ISTOP= NDATPT(ISTATE)
c** When calculating only one potential point
      IF(RDIST.GT.0.d0) THEN
          ISTART= NPNTMX
          ISTOP= NPNTMX
          VDIST= 0.0d0
          BETADIST= 0.d0
          ENDIF
      PETEMP= 0.0d0
      DO  I= ISTART,ISTOP
          BETAFX(I,ISTATE)= 0.0d0
          UAR(I,ISTATE)= 0.d0
          UBR(I,ISTATE)= 0.d0
          TAR(I,ISTATE)= 0.d0
          TBR(I,ISTATE)= 0.d0
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
c** First ... calculate the Extended Morse Oscillator exponent
          DO  I= ISTART,ISTOP
              RVAL= RD(I,ISTATE)
              IF(RDIST.GT.0.d0) RVAL= RDIST
              RDQ= RVAL**qPOT(ISTATE)
              YQ= (RDQ - AREFqq)/(RDQ + AREFqq)
              VAL= BETA(0,ISTATE) 
              DVAL= 0.d0
              DBDB(0,I,ISTATE)= 1.0d0
              YQP= 1.d0
              DO  J= 1, Nbeta(ISTATE)
                  DVAL= DVAL + BETA(J,ISTATE)* DBLE(J)* YQP
                  YQP= YQP*YQ
                  VAL= VAL + BETA(J,ISTATE)*YQP 
                  DBDB(J,I,ISTATE)= YQP
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
c... branch to skip derivatives and inclusion of centrifugal & BOB terms
                  IF(IDAT.LE.-1) GOTO 999  
                  ENDIF
c... now generate remaining partial derivatives
              YPP= 2.0d0*DE(ISTATE)*XTEMP*(1.d0 - XTEMP)
              dYQdRRq= 0.50d0*DVAL*(YQ**2-1.d0)*qPOT(ISTATE)
     1                                                  /RREFq(ISTATE)
              IF(RREFq(ISTATE).LE.0.d0) THEN
                  DBDRe(I,ISTATE)= dYQdRRq   
                  VAL= VAL - (RVAL- RE(ISTATE))*DBDRe(I,ISTATE) 
                  ENDIF
              IPV= IPV+1
              DVtot(IPV,I)= -YPP*VAL       !! derivative w.r.t Re 
              IPV= IPV+1            !! derivative w.r.t. RREFq
              DVtot(IPV,I) = YPP*(RVAL- RE(ISTATE))*dYQdRRq       
              YQP= YPP*(RVAL - RE(ISTATE)) 
              DO  J= 0, Nbeta(ISTATE)
                  IPV= IPV+1
                  DVtot(IPV,I)= YQP        !! derivative w.r.t. \beta_i 
                  YQP= YQP*YQ 
                  ENDDO
ccccc Print for testing
              IF(MOD(I,OSEL(ISTATE)).EQ.0) 
     1  WRITE(10,610) RD(i,ISTATE),vpot(i,istate),BETAFX(i,istate),ULR
              ENDDO
          rewind(10) 
c          write(10,610) (RD(i,ISTATE),vpot(i,istate),BETAFX(i,istate),
c                              i= 1, NDATPT(ISTATE),OSEL(ISTATE))
ccccc
          ENDIF
c........ End preparation of Expanded Morse Potential Function .........

      IF(PSEL(ISTATE).EQ.2) THEN
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** For the  {Morse/Long-Range}  potential.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** For 'normal' inverse-power sum MLR case, with or without damping,
c   set up and prepare 'Dattani-corrected' effective Cm values 
          CmEFF1d= 0.d0
          DO m= 1, NCMM(ISTATE)
              CmEFF(m,ISTATE)= CmVAL(m,ISTATE)
              Cm1D(m)= CmVAL(m,ISTATE)
              CmEFF1D(m)= Cm1D(m)
              MMLR1D(m)= MMLR(m,ISTATE)  !! powers in 1D array for dampF
              ENDDO
          MCMM= NCMM(ISTATE)
cc   ! Assuming {adj} corrections remembered from a previous call .....
cc        IF((RDIST.GT.0.d0).AND.((IDAT.GT.1).OR.(IDAT.LT.0))) GO TO 50
c*** ! ELSE - make (& write) 'Dattani' Cm{adj} correctons here *********
              IF((RDIST.LE.0.d0).OR.(IDAT.EQ.0)) THEN 
                  CALL quadCORR(NCMM(ISTATE),MCMM,NCMMAX,
     1                                 MMLR1D,DE(ISTATE),Cm1D,CmEFF1D)
c... avoids unnecess. quadCORR calls when PEF parameters not changing
                  IF(MCMM.GT.NCMM(ISTATE)) THEN
                      DO m= NCMM(ISTATE),MCMM
                          CmEFF(m,ISTATE)= CmEFF1D(m)
                          ENDDO
                      ENDIF
                  ENDIF
          flush(6)
   50     IF(MMLR(1,ISTATE).LE.0) THEN
c------------------------------------------------------------------------
c** Define value & derivatives of uLR at Re ... first for A-F cases
c... Cm       1    2    3     4     5    6     7     8    9     10     
c... 2x2 {DELTAE, C3s, C3p,  C6s,  C6p, C8s,  C8p}
c    Aubert-Frecon 2x2 treatment of {C3,C6,C8} for Li2 A- or b-state
c... 3x3 {DELTAE, C3s, C3p1, C3p3, C6s, C6p1, C6p3, C8s, C8p1, C8p3}
c    Aubert-Frecon 3x3 treatment of {C3,C6,C8} for Li2 1^3\Pi_g or B-state
c------------------------------------------------------------------------
              CALL AFdiag(RE(ISTATE),NCMM(ISTATE),NCMMax,MMLR1D,CmEFF1D,
     1                       rhoAB(ISTATE),IVSR(ISTATE),IDSTT(ISTATE),
     2                                       ULRe,dLULRedCm,dLULRedRe)
              DO m= 1,NCMM(ISTATE)
                  dLULRedCm(m)= dLULRedCm(m)/ULRe
                  ENDDO
              dLULRedRe= dLULRedRe/ULRe
            ELSE
c** and then for 'normal' inverse-power sum uLR fx.
              CALL dampF(RE(ISTATE),rhoAB(ISTATE),MCMM,NCMMAX,
     1                   MMLR1D,IVSR(ISTATE),IDSTT(ISTATE),Dm,Dmp,Dmpp)
              ULRe= 0.d0
              T1= 0.d0
              DO  m= 1,MCMM
                  IF(rhoAB(ISTATE).LE.0.d0) THEN
                      dLULRedCm(m)= 1.d0/RE(ISTATE)**MMLR(m,ISTATE)
                    ELSE
                      dLULRedCm(m)= Dm(m)/RE(ISTATE)**MMLR(m,ISTATE)
                    ENDIF
                  T0= CmEFF(m,ISTATE)*dLULRedCm(m)
                  ULRe= ULRe + T0
                  T1= T1 + MMLR(m,ISTATE)*T0
                  ENDDO
              dLULRedRe= -T1/(ULRe*RE(ISTATE))
              DO  m= 1,MCMM
                  dLULRedCm(m)= dLULRedCm(m)/ULRe
                  IF(rhoAB(ISTATE).GT.0) dLULRedRe= dLULRedRe + 
     1                                       dLULRedCm(m)*Dmp(m)/Dm(m)
                  ENDDO
            ENDIF
          BINF= DLOG(2.0d0*DE(ISTATE)/ULRe)
          betaINF(ISTATE)= BINF
          IF(APSE(ISTATE).GT.0) THEN
c*** For Pashov-natural-spline exponent coefficient ...
              DO  I= 1,Nbeta(ISTATE)
                  xBETA(I)= yqBETA(I,ISTATE)
                  ENDDO
              BETA(Nbeta(ISTATE),ISTATE)= BINF
              CALL Lkoef(Nbeta(ISTATE),xBETA,rKL,NbetaMX)
              ENDIF
c-----------------------------------------------------------------------
          rewind(10)
          rewind(23)
          DO  I= ISTART,ISTOP
c** Now - generate potential while looping over radial array
              RVAL= RD(I,ISTATE)
              IF(RDIST.GT.0.d0) RVAL= RDIST
              RDp= RVAL**pPOT(ISTATE)
              RDq= RVAL**qPOT(ISTATE)
              YPE= (RDp-REp)/(RDp+REp)
              YP= (RDp-AREFpp)/(RDp+AREFpp)
              YQ= (RDq-AREFqq)/(RDq+AREFqq)
              YPM= 1.d0 - YP
              DYPDRE= 0.5d0*(YP**2 - 1.d0)*pPOT(ISTATE)/RE(ISTATE)
              dYQdRRq= 0.5d0*(YQ**2 - 1.d0)*qPOT(ISTATE)/AREFq
              dYPdRRp= 0.5d0*(YP**2 - 1.d0)*pPOT(ISTATE)/AREFp
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
                  IF(RREFq(ISTATE).LE.0.d0) DBDRe(I,ISTATE)= 
     1                            DBDRe(I,ISTATE)+ (BINF - VAL)*dYPdRRp 
     2                                         + (1.d0-YP)*DVAL*dYQdRRq
                  BVAL= YP*BINF + (1.d0- YP)*VAL
                ELSE
c... now calculate Pashov-spline exponent coefficient & its derivatives
                  BVAL= 0.d0
                  DO  J= 1, Nbeta(ISTATE)
                      DBDB(J,I,ISTATE)= 
     1                               Scalc(YQ,J,npow,xBETA,rKL,NbetaMX)
                      BVAL= BVAL+ DBDB(J,I,ISTATE)*BETA(J,ISTATE)
                      ENDDO
                  DBDRe(I,ISTATE)= -DBDB(npow,I,ISTATE)*dLULRedRe
                ENDIF
              BETAFX(I,ISTATE)= BVAL
              XTEMP= DEXP(-BVAL*YPE)
c** Now begin by generating  uLR(r)
c----- Special Aubert-Frecon cases ------------------------------------
              IF(MMLR(1,ISTATE).LE.0) THEN
c ... generate ULR for Aubert-Frecon type case ...
                  CALL AFdiag(RVAL,NCMM(ISTATE),NCMMax,MMLR1D,CmEFF1D,
     1                                     rhoAB(ISTATE),IVSR(ISTATE),
     2                               IDSTT(ISTATE),ULR,dULRdCm,dULRdR)
c----- End of special Aubert-Frecon Li2 cases ------------------------
                ELSE
c ... for the case of a 'normal' inverse-power sum u_{LR}(r) function
                  ULR= 0.d0
                  CALL dampF(RVAL,rhoAB(ISTATE),MCMM,NCMMAX,
     2                  MMLR1D,IVSR(ISTATE),IDSTT(ISTATE),Dm,Dmp,Dmpp)
                  DO  m= 1,MCMM      
                      IF(rhoAB(ISTATE).LE.0.d0) THEN
                          dULRdCm(m)= 1.d0/RVAL**MMLR(m,ISTATE)
                        ELSE
                          dULRdCm(m)= Dm(m)/RVAL**MMLR(m,ISTATE)
                        ENDIF
                      ULR= ULR + CmEFF(m,ISTATE)*dULRdCm(m)
                      ENDDO
                ENDIF
              XTEMP= XTEMP*ULR/ULRe
c... note ... reference energy for each state is its asymptote ...
              DVDD= XTEMP*(XTEMP - 2.D0)  
              VPOT(I,ISTATE)= DE(ISTATE)*DVDD + VLIM(ISTATE)
              IF(RDIST.GT.0.d0) THEN
                  VDIST= VPOT(I,ISTATE)
                  BETADIST= BVAL
c... branch to skip derivatives and inclusion of centrifugal & BOB terms
                  IF(IDAT.LE.-1) GOTO 999  
                  ENDIF
              YPP= 2.d0*DE(ISTATE)*(1.0d0-XTEMP)*XTEMP   !! == DER
             IPV= IPVSTART+4          !! save count for De, Re, RREFq & RREFp
c... First, generate derivatives w.r.t long-range parameters ... 
              m1= 1
              IF(MMLR(1,ISTATE).LE.0) THEN
c... for Aubert-Frecon diagonalization uLR .....
                  IPV=IPV+1
                  DVtot(IPV,I)= 0.d0      !! derivative w.r.t. splitting=0.0
                  m1= 2
                  ENDIF
c... now derivative w.r.t. Cm's
              DO  m= m1, NCMM(ISTATE)    
                  IPV= IPV+ 1
                  DVtot(IPV,I)= YPP*(dLULRedCm(m)*(1.d0 - YP*YPE)
     1                                               - dULRdCm(m)/ULR)
                  ENDDO
              IF(MCMM.GT.NCMM(ISTATE)) THEN
c ... ideally should ajdust dV/dC3 for C6{adj} and C9{adj} ... but ...
                  ENDIF
c... derivative w.r.t. Re
              DVtot(IPVSTART+2,I)= YPP*(YPE*DBDRe(I,ISTATE)
     1                                       + BVAL*DYPDRE + dLULRedRe)
              IF(APSE(ISTATE).LE.0) THEN
c... derivative w.r.t. De  for 'conventional' power-series exponent
                  DVtot(IPVSTART+1,I)= DVDD + YPP*YPE*YP/DE(ISTATE)
c... derivative w.r.t. RREFp  for 'conventional' power-series exponent
                  DVtot(IPVSTART+4,I)= YPP*YPE*(BINF - VAL)*dYPdRRp
c... derivative w.r.t. RREFq  for 'conventional' power-series exponent
                  DVtot(IPVSTART+3,I)= YPP*YPE*(1.d0-YP)*DVAL*dYQdRRq
                  IF(RREFp(ISTATE).LE.0.d0) THEN
                      DVtot(IPVSTART+3,I)= DVtot(IPVSTART+3,I) 
     1                                           + DVtot(IPVSTART+4,I)
                      DVtot(IPVSTART+4,I)= 0.d0
                      ENDIF
c... if appropriate, add  RREFq {& RREFp}  contribn to  dV/dr_e
                  IF(RREFq(ISTATE).LE.0.d0) THEN
                      DVtot(IPVSTART+2,I) = DVtot(IPVSTART+2,I) 
     1                                           + DVtot(IPVSTART+3,I)
                      DVtot(IPVSTART+3,I)= 0.d0
                      ENDIF
c... derivatives w.r.t. \beta_i for PE-MLR cases...
                  YPP= YPP*YPE*(1 - YP)
                  DO  J= 0, Nbeta(ISTATE)
                      IPV= IPV+1
                      DVtot(IPV,I)= YPP
                      YPP= YPP*YQ
                      ENDDO
                ELSE
c... derivatives w.r.t. \beta_i for Pashov-spline exponent cases...
                  YPP= YPP*YPE
                  DO  J= 1,Nbeta(ISTATE)
                      IPV= IPV+ 1
                      DVtot(IPV,I)= DBDB(J,I,ISTATE)*YPP
                      ENDDO
c... derivative w.r.t. De  for Pashov-Spline expoinent 
                  DVtot(IPVSTART+1,I)= DVDD             !! OK (I think)
     1                   + YPP*DBDB(Nbeta(ISTATE),I,ISTATE)/DE(ISTATE)
                ENDIF
              IF(MOD(I,OSEL(ISTATE)).EQ.0) 
     1  WRITE(23,610) RD(i,ISTATE),vpot(i,istate),BETAFX(i,istate),ULR
              ENDDO                !! end of loop over MLR radial array
ccccc Print for testing
c         write(10,610) (RD(i,ISTATE),vpot(i,istate),BETAFX(i,istate),
c    1                              i= ISTART,ISTOP,OSEL(ISTATE))
  610 FORMAT(f10.4,f15.5,3f20.6)
ccccc End of Print for testing
          ENDIF
          flush(6)
c......... End for Morse/Lennard-Jones(p) potential function ...........
      IF(PSEL(ISTATE).EQ.3) THEN
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** For the Double-Exponential Long-Range (DELR) potential
          AA(ISTATE)= 0.0d0
          BB(ISTATE)= 0.0d0
          dAAdRe= 0.0d0
          dBBdRe= 0.0d0
          dVdBtemp= 0.0d0
c ... First, save uLR powers & coefficients in 1D arrays
          DO  m= 1, NCMM(ISTATE)
              MMLR1D(m)= MMLR(m,ISTATE)
              Cm1D(m)= CmVAL(m,ISTATE)
              ENDDO
c... then get  AA & BB and their derivatives!
          yqRe= (REq - AREFqq)/(REq + AREFqq)       !! next - dyq/dr @ r_e
          dyqRedRe = 2.d0*qPOT(ISTATE)*REq*AREFqq/(Re(ISTATE)*
     1                                               (REq + AREFqq)**2)
          IF(RREFq(ISTATE).LE.0.d0) dyqRedRe= 0.d0
          betaRe= beta(0,ISTATE)
          DbetaRe= 0.d0                !! this is d{beta}/d{y}  at r= Re
          yPOW= 1.d0
          npow= Nbeta(ISTATE)
          POWmax= npow
          IF(npow.GE.1) THEN
              DO j= 1,npow
                  DbetaRe= DbetaRe + j*beta(J,ISTATE)*yPOW
                  yPOW= yPOW*yqRe
                  betaRe= betaRe + yPOW*beta(J,ISTATE)
                  ENDDO
              ENDIF
          CALL dampF(RE(ISTATE),rhoAB(ISTATE),NCMM(ISTATE),NCMMAX,
     1              MMLR1D,IVSR(ISTATE),IDSTT(ISTATE),Dm,Dmp,Dmpp)
          ULRe= 0.d0
          dULRe= 0.d0
          d2ULRe= 0.d0
          IF(MMLR(1,ISTATE).LE.0) THEN
c ... for Aubert-Frecon diagonalization for u_{LR}(r)
              CALL AFdiag(RE(ISTATE),NCMM(ISTATE),NCMMax,MMLR1D,Cm1D,
     1                                     rhoAB(ISTATE),IVSR(ISTATE),
     2                              IDSTT(ISTATE),ULRe,dULRedCm,dULRe)
              dLULRedRe=dULRe/ULRe
            ELSE
c ... for ordinary inverse-power sum u_{LR}(r) 
              DO  m= 1,NCMM(ISTATE)
                  T0= CmVAL(m,ISTATE)/RE(ISTATE)**MMLR(m,ISTATE)
                  IF(rhoAB(ISTATE).GT.0.d0) THEN
                      ULRe= ULRe+ T0*DM(m)
                      dULRe= dULRe+ T0*(Dmp(m) - 
     1                                     Dm(m)*MMLR1D(m)/RE(ISTATE))
                      d2ULRe= d2ULRe + T0*(Dmpp(m) - 
     1                 2.d0*MMLR1D(m)*Dmp(m)/Re(ISTATE) + 
     2                 MMLR1D(m)*(MMLR1D(m)+1.d0)*Dm(m)/Re(ISTATE)**2)
                    ELSE
                      ULRe= ULRe+ T0
                      dULRe= dULRe - T0*MMLR1D(m)/RE(ISTATE)
                      d2ULRe= d2ULRe + T0*MMLR1D(m)*(MMLR1D(m)+1.d0)
     1                                                  /Re(ISTATE)**2
                    ENDIF
                  ENDDO
            ENDIF
c
          AA(ISTATE)= DE(ISTATE) - ULRe - dULRe/betaRe
          BB(ISTATE)= AA(ISTATE) + DE(ISTATE) - ULRe
          dAAdb0 = dULRe/betaRe**2          !! this is d{AA}/d{beta(0)}
          dAAdRe= -dULRe - d2ULRe/dbetaRe + dAAdb0*DbetaRe*dyqRedRe
          dBBdRe= dAAdRe - dULRe
c===== end of calcn. for properties at Re performed, for 1'st point ====
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Now to calculate the actual potential and partial derivatives:
          DO  I= ISTART,ISTOP
c** Start by generating the exponent and its derivative w.r.t. yq
              RVAL= RD(I,ISTATE)
              IF(RDIST.GT.0.d0) RVAL= RDIST    !! for calc at onee point
              RDQ= RVAL**qPOT(ISTATE)
              YQ= (RDQ-AREFqq)/(RDQ+AREFqq)
              YPOW= 1.d0
              npow= NBETA(ISTATE)
              betaFX(I,ISTATE)= beta(0,ISTATE)
              dbetaFX= 0.d0               !! this is  d{beta(r)}/dy_p(r)
              DO  J= 1,npow 
                  dbetaFX= dbetaFX + DBLE(J)*beta(J,ISTATE)*YPOW
                  YPOW= YPOW*YQ
                  betaFX(I,ISTATE)= betaFX(I,ISTATE) +
     1                                         beta(J,ISTATE)*YPOW
                  ENDDO
c** Calculate some temporary variables.
              XTEMP= DEXP(-betaFX(I,ISTATE)*(RVAL-RE(ISTATE)))
c** Now to calculate uLR and the actual potential function value
              ULR= 0.0d0
c*** For Aubert-Frecon alkali dimer nS + nP diagonalization u_{LR}(r)
              IF(MMLR(1,ISTATE).LE.0) THEN
                  CALL AFdiag(RVAL,NCMM(ISTATE),NCMMax,MMLR1D,Cm1D,
     1                                     rhoAB(ISTATE),IVSR(ISTATE),
     2                               IDSTT(ISTATE),ULR,dULRdCm,dULRdR)
                ELSE
c... now for 'regular' inverse-power sum u_{LR}(r)
                  CALL dampF(RVAL,rhoAB(ISTATE),NCMM(ISTATE),NCMMAX,
     1                  MMLR1D,IVSR(ISTATE),IDSTT(ISTATE),Dm,Dmp,Dmpp)
                  DO  m= 1,NCMM(ISTATE)
                      ULR= ULR + CmVAL(m,ISTATE)*Dm(m)/
     1                                            RVAL**MMLR(m,ISTATE)
                      ENDDO
                ENDIF
c... END of u_{LR}(r) calculation ........
cc            REWIND(30)
cc            WRITE(30,*) RVAL,ULR     !! test printout for error check
              VPOT(I,ISTATE)= (AA(ISTATE)*XTEMP - BB(ISTATE))*XTEMP 
     1                                            - ULR + VLIM(ISTATE)
              IF((I.EQ.ISTART).AND.(I.EQ.ISTOP)) THEN    
                  VDIST= VPOT(I,ISTATE)   !! if getting V(r) at onee point
                  betaDIST= betaFX(I,ISTATE)
                  ENDIF                
c... branch to skip derivatives and inclusion of centrifugal & BOB terms
              IF(IDAT.LE.-1) GOTO 999  
c-----------------------------------------------------------------------
c** Now, calculate the partial derivatives ...
c-----------------------------------------------------------------------
              IPV= IPVSTART + 1
c ... first, derivative of the potential w.r.t. De
              DVtot(IPV,I)= (XTEMP - 2.0d0)*XTEMP 
c** Now to calculate the derivative of the potential w.r.t. Re
              Btemp= (2.0d0*AA(ISTATE)*XTEMP - BB(ISTATE))*XTEMP
              IPV= IPV + 1
              DVtot(IPV,I)= betaFX(I,ISTATE)*Btemp 
     1                                 + XTEMP*(dAAdRe*XTEMP - dBBdRe)
              Btemp= (RVAL- RE(ISTATE))*Btemp
              IF(RREFq(ISTATE).LE.0.d0) 
     1             DVtot(IPV,I)= DVtot(IPV,I) - Btemp*DbetaRe*dyqRedRe
c... ** when calculating Cm derivatives, dULRe/dCm has been excluded ** 
              DO m= 1,NCMM(ISTATE)
                  dULRepdCm(m)=0.d0
                  ENDDO
c...
              IF((NCMM(ISTATE).GE.4).AND.(MMLR(1,ISTATE).LE.0)) THEN
c... derivatives w.r.t long-range parameters for Aubert-Frecon uLR
                  IPV=IPV+1
                  DVtot(IPV,I)= 0.d0
                  DO m=2,NCMM(ISTATE)
                      IPV=IPV+1
                      DVtot(IPV,I)= XTEMP*(2*dULRedCm(m)+
     1                              dULRepdCm(m)/betaRe - XTEMP*
     2                              (dULRedCm(m)+dULRepdCm(m)/betaRe))
     3                              - dULRdCm(m)
                      ENDDO
              ELSE
c ... derivative w.r.t. Cm's for ordinary MLR/MLJ case ...
                  DO  m= 1, NCMM(ISTATE)
                      IPV= IPV+ 1
                      DVtot(IPV,I)= XTEMP*(2*dULRedCm(m)+
     1                              dULRepdCm(m)/betaRe - XTEMP*
     2                              (dULRedCm(m)+dULRepdCm(m)/betaRe))
     3                              - dULRdCm(m)
                      ENDDO
                ENDIF
c
c ... finally, derivatives of the potential w.r.t. the \beta_i
              Btemp2= (Xtemp - 1.d0)*Xtemp*dAAdb0
              IPV= IPV+ 1
              DVtot(IPV,I)= - Btemp + Btemp2
              DO  J= 1,npow
                  IPV= IPV+ 1
                  Btemp= Btemp*YQ
                  Btemp2= Btemp2*yqRe
                  DVtot(IPV,I)= - Btemp + Btemp2
                  ENDDO
              IF(npow.LT.POWmax) THEN
                  DO J= npow+1, POWmax
                      IPV= IPV+1
                      DVtot(IPV,I)= 0.0d0
                      ENDDO
                  ENDIF
c *** DBDRe and DBDB is used in uncertainty calculation, see fununc.f
c
c??? QUESTION ,,, IS the parameter count correct here ?????
c
              DBDRe(I,ISTATE)= 0.d0
              IF(RREFq(ISTATE).LE.0) DBDRe(I,ISTATE)= 1.d0
              DBDB(0,I,ISTATE)= 1.0d0
              DO  J= 1, npow
                  DBDB(J,I,ISTATE)= DBDB(J-1,I,ISTATE)*YQ
                  ENDDO
              IF(npow.LT.POWmax) THEN
                  DO J= npow+1,POWmax
                      DBDB(J,I,ISTATE)= 0.0d0
                      ENDDO
                  ENDIF
              ENDDO
ccccc Print for testing
      rewind(10)
      write(10,610) (RD(i,ISTATE),vpot(i,istate),betaFX(i,istate),
     1                              i= 1, NDATPT(ISTATE),OSEL(ISTATE))
ccccc End of Print for testing
          ENDIF
c....... End Double-Exponential Long-Range Potential Function ........
      IF(PSEL(ISTATE).EQ.4) THEN
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
            RDp= RVAL**pPOT(ISTATE)
            YP= (RDp-REp)/(AGPEF(ISTATE)*RDp + BGPEF(ISTATE)*REp)
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
c... branch to skip derivatives and inclusion of centrifugal & BOB terms
                IF(IDAT.LE.-1) GOTO 999  
                ENDIF
            DVAL= DVAL*BETA(0,ISTATE)*YP
c** Now to calculate the partial derivatives
            DVDD= 0.d0
            IPV= IPVSTART + 1
c ... derivative of the potential w.r.t. Re
            DVtot(IPV,I)= -DVAL*REp*RDp*(AGPEF(ISTATE)+BGPEF(ISTATE))
     1                  *(pPOT(ISTATE)/RE(ISTATE))/(AGPEF(ISTATE)*RDp +
     2                                           BGPEF(ISTATE)*REp)**2
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
        IF(RDIST.LE.0.d0) THEN
            rewind(10)
            write(10,612) (RD(i,ISTATE),vpot(i,istate),
     1                              i= 1, NDATPT(ISTATE),OSEL(ISTATE))
            ENDIF
  612 FORMAT(/(f10.4,f15.5))
c????
        ENDIF
c.......................... End Surkus GPEF Potential Function .........
      IF(PSEL(ISTATE).EQ.5) THEN
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** For the Tiemann 'HPP' Polynomial Potential Energy Function.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          BT= BETA(Nbeta(ISTATE)+1, ISTATE)
          Rinn= BETA(Nbeta(ISTATE)+2, ISTATE)
          Rout= BETA(Nbeta(ISTATE)+3, ISTATE)
c** With long-range tail an NCMM-term inverse-power sum, adjust De and 
c  add 1 more inverse-power term  CMMp2/r**{m_{last}+2}} to ensure continuity
c  and smoothness at  Rout
          XRO= (Rout - RE(ISTATE))/(Rout+ BT*RE(ISTATE))
          YPP= 1.d0
          VX= 0.d0
          dVX= 0.d0
          DO  J= 1, Nbeta(ISTATE)
              dVX= dVX+ J*YPP*BETA(J,ISTATE)
              YPP= YPP*XRO
              VX= VX+ YPP*BETA(J,ISTATE)
              ENDDO
          dXRO= (RE(ISTATE)+ BT*RE(ISTATE))/(Rout + BT*RE(ISTATE))**2
          dXRI= (RE(ISTATE)+ BT*RE(ISTATE))/(Rinn + BT*RE(ISTATE))**2
c***  dXRO= dX(r)/dr @ r=R_{out}   &  dXRORe= dX(r)/dr_e @ r=R_{out}
          dXROdRe= -dXRO*Rout/RE(ISTATE)
          dXRIdRe= -dXRI*Rinn/RE(ISTATE)
          d2XROdRe = (1.d0 + BT)*(Rout - BT*RE(ISTATE))/
     1                                       (Rout + BT*RE(ISTATE))**3
          d2XRIdRe = (1.d0 + BT)*(Rinn - BT*RE(ISTATE))/
     1                                       (Rinn + BT*RE(ISTATE))**3
          dVX= dVX*dXRO
c  VX={polynomial part V_X @ Rout} and dVX is its derivative w.r.t. r
          uLR= 0.d0
          CMMp2= 0.d0
          DO  J= 1, NCMM(ISTATE)
              B5= CmVAL(J,ISTATE)/Rout**MMLR(J,ISTATE)
              uLR= uLR+ B5
              CMMp2= CMMp2+ MMLR(J,ISTATE)*B5
              ENDDO
          MMp2= MMLR(NCMM(ISTATE),ISTATE)+2    
          fRO= Rout**(MMp2+1)/MMp2            !! factor for derivatives
          CMMp2= (dVX- CMMp2/Rout)*fRO
c??? zero out C5(A) for Mg2 to match KNoeckel ??  !! as per Marcel          
ccc       IF(ISTATE.GT.1) CMMp2= 0.d0
          DE(ISTATE)= uLR + VX + CMMp2/Rout**MMp2
c** CMMp2= C_{m_{last}+2}:  now get the updated value of  DE(ISTATE)
c** now ... Determine analytic function attaching smoothly to inner wall
c   of polynomial expansion at  R= Rinn < Rm
          XRI= (Rinn - RE(ISTATE))/(Rinn+ BT*RE(ISTATE))
          YPP= 1.d0
          B5= VLIM(ISTATE) - DE(ISTATE)
          A1= 0.d0
          A2= 0.d0
          DO  J= 1, Nbeta(ISTATE)
              A2= A2+ J*YPP*BETA(J,ISTATE)
              YPP= YPP*XRI
              A1= A1+ YPP*BETA(J,ISTATE)
              ENDDO
          A2= A2*dXRI                     !! dXRI= dX(r)/dr @ r=R_{inn}
          A2= -A2/A1
c** Extrapolate inwardly with the exponential: B5 + A1*exp(-A2*(R-Rinn))
c... but first collect some common factors for the derivatives
          dCmp2dRe= 0.d0
          dDeROdRe= 0.d0
          dDeRIdRe= 0.d0
          XROpw= 1.d0/XRO**2
          XRIpw= 1.d0/(A1*XRI**2)
          ROmp2= 1.d0/Rout**MMp2
          BIrat= A2/A1
          DO J=1, Nbeta(ISTATE)
c... first ... outer boundary factors & derivatives w.r.t. beta_i             
              dCmp2dRe= dCmp2dRe + J*(J-1)*BETA(J,ISTATE)*XROpw
              XROpw= XROpw*XRO                 !! power now   (J-1)
              dDeROdRe= dDeROdRe + J*BETA(J,ISTATE)*XROpw
              dCmp2(J)= J*XROpw*dXRO*fRO
              dDe(J)= XROpw*XRO + ROmp2*dCmp2(J)   !! uses power  J
c... then ... inner boundary factors & derivatives w.r.t. beta_i 
              dBIdRe= dBIdRe + J*(J-1)*BETA(J,ISTATE)*XRIpw
              XRIpw= XRIpw*XRI                 !! power now   (J-1)
              dDeRIdRe= dDeRIdRe + J*BETA(J,ISTATE)*XRIpw
              dAI(J)= XRIpw*XRI - dDe(J)       !! add term with power J
              dBI(J)= J*XRIpw*dXRI - BIrat*dAI(J)  
              ENDDO
          dCmp2dRe= (dDeROdRe*d2XROdRe + dCmp2dRe*dXRO*dXROdRe)*fRO
          dDedRe=  dDeROdRe*dXROdRe + ROmp2*dCmp2dRe
          dAIdRe= -dDeDRe + dDeRIdRe*A1*dXRI
          dBIdRe= - dDeRIdRe*d2XRIdRe- dBIdRe*dXRI*dXRIdRe- BIrat*dAIdRe
          DO  I= ISTART, ISTOP
c*** Now ... loop to generate the potential ...
              RVAL= RD(I,ISTATE)
              IF(RDIST.GT.0) RVAL= RDIST
              YP= (RVAL - RE(ISTATE))/(RVAL + BT*RE(ISTATE))
              IF(RVAL.LE.Rinn) THEN
c ... for exponential inward extrapolation2 ...
                  EXPBI= DEXP(-A2*(RVAL- Rinn))
                  VPOT(I,ISTATE)= B5 + A1*EXPBI
                  IPV= IPVSTART+ 1            !! count for Re, NOT for De
                  DO  J=1, NCMM(ISTATE)
                      IPV= IPV+1      !! count for derivatives w.r.t. Cm's
                      DVtot(IPV,I)= 0.d0
                      ENDDO 
                  DO J=1, Nbeta(ISTATE)
                      IPV= IPV+ 1            !! counter for \beta_i
                      DVtot(IPV,I)= - dDe(J) + EXPBI*(dAI(J) 
     1                                       - A1*(RVAL- Rinn)*dBI(J))
                      ENDDO 
                  DVTOT(IPVSTART+1,I)= -dDe(J) + EXPBI*(dAIdRe 
     1                                       - A1*dBIdRe*(RVAL- Rinn))
                ELSEIF(RVAL.LE.Rout) THEN
c ... for 'middle' well region X-polynomial power series ...                
                  IPV= IPVSTART+ 1            !! count for Re, NOT for De
                  DO  J=1, NCMM(ISTATE)
                      IPV= IPV+1              !! for derivatives w.r.t. Cm's
                      DVtot(IPV,I)= 0.d0
                      ENDDO 
                  VX= VLIM(ISTATE) - DE(ISTATE)
                  dVdRe= 0.d0
                  YPOW= 1.d0 
cc                IF(DABS(YP).GT.0.d0) YPOW= 1.d0/YP     !! if start @ J=0
                  DO J=1, Nbeta(ISTATE)
                      IPV= IPV+ 1            !! counter for \beta_i
                      dVdRe= dVdRe + J*BETA(J,ISTATE)*YPOW  !! 
                      YPOW= YPOW* YP         !! brings power up to J
                      DVTOT(IPV,I)= -dDe(J) + YPOW
                      VX= VX + BETA(J,ISTATE)*YPOW
                      ENDDO
                  VPOT(I,ISTATE)= VX
                  DVTOT(IPVSTART+1,I)= - dDedRe - dVdRe*RVAL*(BT+1.d0)
     1                                      /(RVAL + BT*RE(ISTATE))**2
                ELSEIF(RVAL.GT.Rout) THEN
c ... for Van der Waals tail region with added inverse-power term
                  IPV= IPVSTART+ 1            !! count for Re, NOT for De
                  A3= VLIM(ISTATE)
                  DO  J= 1, NCMM(ISTATE)
                      IPV= IPV+1              !! for derivatives w.r.t. Cm's
                      DVtot(IPV,I)= 0.d0
                      A3= A3- CmVAL(J,ISTATE)/RVAL**MMLR(J,ISTATE)
                      ENDDO
                  RMMp2= 1.d0/RVAL**MMp2
                  VPOT(I,ISTATE)= A3 - CMMp2*RMMp2
                  DO J=0, Nbeta(ISTATE)
                      IPV= IPV+ 1            !! counter for \beta_i
                      DVTOT(IPV,I)= -dCmp2(J)*RMMp2
                      ENDDO
                  DVTOT(IPVSTART+1,I)= -dCmp2dRe*RMMp2 
                ENDIF
              ENDDO
** end of loop over distance array
          IF(RDIST.GT.0) THEN
              VDIST= VPOT(I,ISTATE)
              BETADIST= 0.d0 
c... branch to skip derivatives and inclusion of centrifugal & BOB terms
              IF(IDAT.LE.-1) GOTO 999  
              ENDIF
          IF(RDIST.LE.0.d0) THEN
              rewind(10)
              write(10,612) (RD(i,ISTATE),vpot(i,istate),
     1                              i= 1, NDATPT(ISTATE),OSEL(ISTATE))
              ENDIF
          ENDIF
c........... End Tiemann Potential Energy Function .................
      IF(PSEL(ISTATE).EQ.6) THEN
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** For the generalized TANG-Toennies-type potential with 4-term exponent
c  & 5-term pre-exponential factor minus and damped (s=+1) repulsion terms
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c ... first ... save uLR powers in a 1D array
          DO  m= 1, NCMM(ISTATE)
              MMLR1D(m)= MMLR(m,ISTATE)
              ENDDO
         rhoINT= rhoAB(ISTATE)/3.13d0     !! remove btt(IVSR(ISTATE)/2) 
c** Now, calculate A and De  using the input values of b and r_e 
          REQQ= RE(ISTATE)
          RVAL= RDIST
          DO I= ISTART, ISTOP
              IF(RDIST.LE.0.d0) RVAL= RD(I,ISTATE) 
              T0= RVAL*(BETA(1,ISTATE) + RVAL*(BETA(2,ISTATE)))
     1                   + (BETA(3,ISTATE) + BETA(4,ISTATE)/RVAL)/RVAL
              ATT= (BETA(5,ISTATE) + RVAL*(BETA(6,ISTATE) + RVAL*
     1                       (BETA(8,ISTATE) + RVAL*BETA(9,ISTATE)))) 
     2                                           + BETA(7,ISTATE)/RVAL
              CALL dampF(RVAL,rhoINT,NCMM(ISTATE),NCMMAX,
     1                  MMLR1D,IVSR(ISTATE),IDSTT(ISTATE),Dm,Dmp,Dmpp)
              VATT= 0.d0
              DO M= 1,NCMM(ISTATE)
                  VATT= VATT+ CmVAL(m,ISTATE)*Dm(m)/RVAL**MMLR1D(m)
                  ENDDO
              VPOT(I,ISTATE)= ATT*EXP(-T0)- VATT
c!! Special insert for Shen-Tang Be2  PRA 88, 011517 (2013)
c             VPOT(I,ISTATE)= VPOT(I,ISTATE) - 9.486575760D+05
c    1               *DEXP(-RVAL*(1.113237666d0+ RVAL*0.2764004206d0))
c             write(32,832) RVAL, Vatt, VPOT(I,ISTATE)
c 832 Format(F8.3, 2F10.3)
              IF(VPOT(I,ISTATE).LT.VMIN) THEN
                  VMIN= VPOT(I,ISTATE)
                  REQQ= RVAL  
                  ENDIF
              ENDDO
          IF(RDIST.GT.0.d0) THEN
              VDIST= VPOT(ISTOP,ISTATE)
              BETADIST= 0.d0
              ENDIF    
          IF(RDIST.LE.0.d0) THEN
              IF(ISTOP.GT.ISTART) WRITE(6,602) VMIN,REQQ
  602 FORMAT('    Extended TT potential has   VMIN=',f9.4,'   at  RMIN='
     1                                                          f8.5)
c...... Print for testing            
              rewind(10)
              write(10,612) (RD(i,ISTATE),vpot(i,istate),
     1                              i= 1, NDATPT(ISTATE),OSEL(ISTATE))
              ENDIF
          ENDIF
c........... End Tang-Toennies Potential Energy Function .................
      IF(PSEL(ISTATE).EQ.7) THEN
          IF(Nbeta(ISTATE).EQ.5) THEN    
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** For the Aziz'ian HFD-ABC potential
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
              A1= BETA(1,ISTATE)
              A2= BETA(2,ISTATE)
              A3= BETA(3,ISTATE)
              X= RDIST
              DO  I= ISTART, ISTOP
                  IF(RDIST.LE.0.d0) X= RD(I,ISTATE)
                  VATT= 0.d0
                  DO M= 1,NCMM(ISTATE)
                      VATT= VATT+ CmVAL(m,ISTATE)/X**MMLR(m,ISTATE)
                      ENDDO
                  IF(X.LT.A2) VATT= VATT*DEXP(-A1*(A2/X -1.d0)**A3)
                  VPOT(I,ISTATE)= AA(ISTATE)*
     1                                    (X/RE(ISTATE))**BETA(5,ISTATE)
     1                   *EXP(-X*(BB(ISTATE) + X*BETA(4,ISTATE))) - VATT
                  ENDDO
              IF(RDIST.GT.0.d0) THEN
                  VDIST= VPOT(ISTOP,ISTATE)
                  BETADIST= 0.d0
                  ENDIF    
            ELSEIF(Nbeta(ISTATE).EQ.2) THEN
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** For the Aziz'ian HFD-D potential
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c ... first ... save uLR powers in a 1D array
              DO  m= 1, NCMM(ISTATE)
                  MMLR1D(m)= MMLR(m,ISTATE)
                  ENDDO
              X= RDIST
              DO  I= ISTART, ISTOP
                  IF(RDIST.LE.0.d0) X= RD(I,ISTATE)
                  CALL dampF(X,rhoAB(ISTATE),NCMM(ISTATE),NCMMAX,
     1                  MMLR1D,IVSR(ISTATE),IDSTT(ISTATE),Dm,Dmp,Dmpp)
                  VATT= 0.d0
                  DO M= 1,NCMM(ISTATE)
                      VATT= VATT+ CmVAL(m,ISTATE)*Dm(m)/X**MMLR1D(m)
                      ENDDO
                  VATT= VATT*(1 - (rhoAB(ISTATE)*X/bohr)**1.68d0
     1                            *DEXP(-0.78d0*rhoAB(ISTATE)*X/bohr))
                  VPOT(I,ISTATE)= AA(ISTATE)*
     1                                  (X/RE(ISTATE))**BETA(2,ISTATE)
     1                 *EXP(-X*(BB(ISTATE) + X*BETA(1,ISTATE))) - VATT
                  ENDDO
              IF(RDIST.GT.0.d0) THEN
                  VDIST= VPOT(ISTOP,ISTATE)
                  BETADIST= 0.d0
                  ENDIF    
            ENDIF
c...... Print for testing            
          IF(RDIST.LE.0.d0) THEN  
              rewind(10)
              write(10,612) (RD(i,ISTATE),vpot(i,istate),
     1                              i= 1, NDATPT(ISTATE),OSEL(ISTATE))
              ENDIF
          ENDIF
c........... End Aziz'ian HFD-ABC & D Potential Energy Function ........
  700 IF((IDAT.LE.0).AND.(RDIST.GT.0)) GOTO 999
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
                          IPV= IPV+ 1            !! VAL is power series sum
                          DVtot(IPV,I)= YPMA*YQP
                          ENDDO
                      ENDIF
                  IPV= IPV + 1
                  DVtot(IPV,I)= YPA
                  IF(LRad(ISTATE).LE.0) THEN
                      UAR(I,ISTATE)= VAL*YPMA + 
     1                                      YPA*UA(NUA(ISTATE),ISTATE)
                    ELSE                !! Add in the \delta{Cm} terms
                      dCmASUM = UA(NUA(ISTATE),ISTATE)
                      DO m= 1,NCMM(ISTATE)
                          dCmASUM = dCmASUM +
     1                           dCmA(m,ISTATE)/(RVAL**MMLR(m,ISTATE))
                          ENDDO
                      UAR(I,ISTATE)= VAL*YPMA + YPA*dCmASUM
                    ENDIF  
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
                  IF(LRad(ISTATE).EQ.0) THEN  !! NO  \delta{Cm} terms
                      UBR(I,ISTATE)= VAL*YPMB + 
     1                                      YPB*UB(NUB(ISTATE),ISTATE)
                    ELSE                !! Add up the \delta{Cm} terms
                      dCmBSUM = UA(NUB(ISTATE),ISTATE)
                      DO m= 1,NCMM(ISTATE)
                          dCmBSUM = dCmBSUM +
     1                           dCmB(m,ISTATE)/(RVAL**MMLR(m,ISTATE))
                          ENDDO
                      UBR(I,ISTATE)= VAL*YPMB + YPB*dCmBSUM
                    ENDIF  
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
c ... Now ... derivatives of R_{na}(A) w,r,t, expansion coefficients
                  VAL= TA(0,ISTATE)
                  DVAL= 0.d0
                  IPV= IPVSTART + 1
                  DVtot(IPV,I)= YPMA*RM2   !! deriv. w.r.t. t_0
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
                  DVtot(IPV,I)= YPA*RM2     !! deriv w.r.r. t_{\inf}
                  TAR(I,ISTATE)= VAL*YPMA + YPA*TA(NTA(ISTATE),ISTATE)
c ... and derivative of R_{na}(A) w.r.t. Re ... 
                  DTADRe(I,ISTATE)= (-HReQ*(1.d0 - YQA**2)*YPMA*DVAL
     1       + HReP*(1.d0 - YPA**2)*(VAL- TA(NTA(ISTATE),ISTATE)))*RM2
c!!! temorary test printing !!!!!n!!!!!!
cc                write(14,699) RVAL,TAR(I,ISTATE),(DVtot(J,I),
cc   1                                       J=IPVSTART+1,IPV)
cc699 FORMAT(f9.4,1P,10D15.7)
c!!! temorary test printing !!!!!!!!!!!
c!!!!!!!!!!!!!! incomplete -how is IPVSTART initialized for NUA, NUB, NTA, NTB
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
c!!! temorary test printing !!!!!!!!!!!
c!!               write(15,699) RVAL,TAR(I,ISTATE),(DVtot(J,I),
c!!  1                                       J=IPVSTART+1,IPV)
c!!! temorary test printing !!!!!!!!!!!
              ENDDO
          ENDIF
c.... END of treatment of non-adiabatic centrifugal BOB function........
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
c.... END of treatment of Lambda/2-Sigma centrifugal BOB function.......
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++++ Test for inner wall inflection , and if it occurs, replace inward
c++++ potential with linear approximation +++++++
      IF(PSEL(ISTATE).NE.5) REQQ= RE(ISTATE)
c!!!! temporary fix to handle Sheng/Tang Be2 case       
      I1= (REQQ -RD(1,ISTATE))/(RD(2,ISTATE)-RD(1,ISTATE))
      IF((I1.GT.3).AND.(RDIST.LE.0)) THEN !! skip check on 1-point CALLs
          NIFL=0    !! NIFL is No. (+) to (-) curv. inflection points @ R < Re
          SLB= +1.d0
          SL= 0.d0
          DO  I= I1-2, 1, -1
              SLBB= SLB
              SLB= SL
              SL= VPOT(I,ISTATE) - VPOT(I+1,ISTATE)
              IF((SL.LE.SLB).AND.(SLB.GE.SLBB)) THEN
                  NIFL= NIFL+ 1
                  WRITE(6,606) SLABL(ISTATE),RD(I,ISTATE),VPOT(I,ISTATE)
                  IF(NIFL.LE.MAXMIN(ISTATE))  THEN  !? prob if inner well deeper
                      IF(VPOT(I,ISTATE).GE.VLIM(ISTATE)) THEN
                          DO  J= I,1,-1             !! Only for wall above VLIM
                              VPOT(J,ISTATE)= VPOT(I,ISTATE) + (I-J)*SL
                              ENDDO
                          WRITE(6,608)
                          GOTO 66
                          ENDIF
                      ENDIF
                  ENDIF
              ENDDO
          ENDIF
   66 CONTINUE
  606 FORMAT(12('===')/'!*!* Find State ',A3,' inner-wall inflection at 
     1   R=', f6.4,'   V=',f11.1 '..... !*!*')
  608 FORMAT(5x,'... and ... extrapolate repulsive wall inward from the
     1re as a LINEAR function'/12('==='))
c++++++++++++End of Inner Wall Test/Correction code+++++++++++++++++++++
c======================================================================
c** For simulation & fitting of tunneling width data ......
c** At the one distance RDIST calculate total effective potential VDIST
c     including (!!) centrifugal and Lambda/2Sigma doubling terms, and 
c     get their partial derivatives w.r.t. Hamiltonian parameters dVdPk.
c
      IF((RDIST.GT.0).AND.(IDAT.GT.0).AND.(IDAT.LT.NDATAMX)) THEN 
          IISTP= ISTP(IB(IDAT))
cccccccc
c         WRITE (40,644) IISTP,RDIST,RVAL,VDIST,I,NDATPT(ISTATE)
c 644 FORMAT ('IISTP =',I3,' RDIST =',G16.8,' RVAL =',G16.8,
c    &          ' VDIST =',G16.8,' I =',I6,' NDATPT =',I6)
cccccccc
          BFCT= 16.857629206d0/(ZMASS(3,IISTP)*RDIST**2)
          JFCT= DBLE(JPP(IDAT)*(JPP(IDAT)+1))
          IF(IOMEG(ISTATE).GT.0) JFCT= JFCT - IOMEG(ISTATE)**2
          IF(IOMEG(ISTATE).EQ.-2) JFCT= JFCT + 2.D0 !! ??? for Li2(A,c)??
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

