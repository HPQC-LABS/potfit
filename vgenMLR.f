c***********************************************************************
      SUBROUTINE VGENmlr(RDIST,VDIST,dVdR,d2VdR2)
c***********************************************************************
c** This subroutine will generate and return function values and radial
c    derivatives of Morse/Long-Range potentials as required for 
c   semiclassical calculation (with quantum corrections) of virial 
c   coefficients and their (analytical) derivatives w.r.t. the potential
c   parameters for use in direct potential fits
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++  COPYRIGHT 2015  by  R.J. Le Roy, Dept. of Chemistry, 
c    Univ. of Waterloo, Waterloo, Ontario, Canada   +
c    This software may not be sold or any other commercial use made    +
c      of it without the express written permission of the authors.    +
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c           ----- Version of 29 March 2015 --------
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** On entry:
c   RDIST:  at the 1 input distances RDIST, calculate potl & derivs
c  * return potential function at those points as VDIST, and the
c     first and second radial derivartives as  dVdR  &  d2VdR2
c  * If RDIST.le.0  read in MLR potential parameters and print description
c    of model
c-----------------------------------------------------------------------
c** On exit via common blocks:
c    dVtot is array or p. derivatives w.r.t. the  2+NCMM+Nbeta+1  param
c    DBDB & DBDRe are p.derives of beta(r) w.r.t. \beta_i & Re, respectively
c
c** Temp:
c    BTEMP     is used to represent the sum used for dV/dRe.
c              is used in GPEF for De calculations.
c    betaINF      is used to represent the  beta(\infty)  value.
c    YP        is used to represent (R^p-Re^p)/(R^p+Re^p) or R-Re.
c    XTEMP     is used to represent (uLR/uLR_e)* exp{-betaINF*RTEMP}
c    AZERO     is used for the trial exponential calculations.
c    AONE      is used for the trial exponential calculations.
c    ATWO      is used for the trial exponential calculations.
c    AZTEMP    is used in the MMO trial exponential calculations.
c              is used in the GPEF = (a+b)/k
c    AOTEMP    is used in the GPEF = [a(k+1)-b(k-1)]/k
c    ATTEMP    is used in the GPEF = [a^2(k+1)-b^2(k-1)]/k
c    ARTEMP    is used in the GPEF = [a^3(k+1)-b^3(k-1)]/k
c    FSW       is used to represent the MLJ switching function.
c======================================================================
        REAL*8 RDIST,VDIST,dVdR,d2VdR2,DSCM,REQ,Rref,rhoAB,ULR,ULRE,
     1  betaINF,PVSR,VAL,DVAL,D2VAL,DVALdRe,AREF,REP,AREFp,AREFq,
     2  BETADIST,RVAL,RINV,RDp,RDq,YPE,YP,YQ,YPM,DYPDRE,DYQDRE,DYPEDRE,
     3  DBDRe,DBDR,dULRdR,d2ULRdr,DYPDR,DYPEDR,DYQDR,D2YPDR,D2YPEDR,
     4  D2YQDR,DYPpDRe,DYPEpDRe,DYPEdRe,DYQpDRe,YPP,YPPP,XTEMP,DXTEMP,
     5  D2XTEMP,dULRdR,d2ULRdR,dLULRdR,dDULRdr,dLULRedRe,T0,T1,dBpdRe,
     6  d2BdR,dYBdRe,dYBpdRe,d2ULRdR,d2DULRdR,DVDD,
     6  CMM(20),CmEFF(20),DM(30),DMP(20),DMPP(20),BETA(0:20),
     7  dULRdCm(20),dLULRdCm(20),dLULRedCm(20),dULRdRCm(20),bTT(-2:2),
     8  cDS(-4:4),bDS(-4:4)
      INTEGER Nbeta,PPAR,QPAR,NCMM,MCMM,IVSR,LVSR,IDSTT,IPVSTART,I,J,m,
     1  IPV,MMLR(20)
c** Common block to return derivatives w.r.t. potential parameters
      REAL*8 dVdPk(30),DVtot(30),DBDB(0:30),dXdP(30),dXpdP(30),dVpdP(30)
      COMMON /VDER/dVtot,DBDB,DBDRe
c-----------------------------------------------------------------------
      DATA bTT/2.10d0,2.44d0,2.78d0,3.13d0,3.47d0/
      DATA bDS/2.50d0,2.90d0,3.30d0,3.69d0,3.95d0,0.d0,4.53d0,0.d0,
     1          4.99d0/
      DATA cDS/0.468d0,0.446d0,0.423d0,0.405d0,0.390d0,0.d0,0.360d0,
     1            0.d0,0.340d0/
      SAVE 

      IF(RDIST.LE.0.d0) THEN
c** if(RDIST.LE.0) READ MLR parameters & write description of potential
c=======================================================================
c** Define MLR potential [as per J.Phys.Chem. A117, 13373(2013).
c=======================================================================
          READ(5,*) Nbeta, PPAR, QPAR
          READ(5,*) DSCM, REQ, Rref
c** For MLR, DELR, HFD, Tang-Toennies or Tiemann-polynomial potentials .....
          READ(5,*) NCMM, rhoAB, IVSR, IDSTT
          DO  I=1,NCMM
              READ(5,*) MMLR(I), CMM(I)
              ENDDO
c-----------------------------------------------------------------------
          DO  I=0,Nbeta
              READ(5,*) BETA(I)
              ENDDO
c-----------------------------------------------------------------------
c** for first need to define ULRE an print potential description
          ULRe= 0.d0
c!! Initialize CmEFF on VERY first call from MAIN prog.
          DO m= 1, NCMM
              CmEFF(m)= CMM(m)
              ENDDO
          MCMM= NCMM
c=======================================================================
cc   ! As appropriate - make (& write) 'Dattani' Cm{adj} correctons here
          IF((MMLR(1).EQ.6).AND.(NCMM.GE.4)) THEN
c... First, consider C6/C12adj(C14adj) for MMLR(m)={6,8,10,(11),12,14} case
              IF(MMLR(4).EQ.12) THEN      ! explicitly MMLR(4)=12
                  CmEFF(4)= CMM(4) + 0.25D0*CMM(1)**2/DSCM
                  WRITE(6,710) MMLR(4),MMLR(4),CmEFF(4)
                  IF(NCMM.GE.5) THEN    ! assuming MMLR(5)=14
                      CmEFF(5)= CMM(5) + 0.5d0*CMM(1)*CMM(2)/DSCM
                      WRITE(6,710) MMLR(5),MMLR(5),CmEFF(5)
                      ENDIF
                ELSE      !! Assuming explicitly MMLR(2)=11 & MMLR(5)=12
                  CmEFF(5)= CMM(5) + 0.25D0*CMM(1)**2/DSCM
                  WRITE(6,710) MMLR(5),MMLR(5),CmEFF(5)
                  IF(NCMM.GE.6) THEN             ! implicitly MMLR(6)=14
                      CmEFF(6)= CMM(6)+0.5D0*CMM(1)*CMM(2)/DSCM
                      WRITE(6,710) MMLR(6),MMLR(6),CmEFF(6)
                      ENDIF
                ENDIF
              ENDIF
          IF((MMLR(1).EQ.5).AND.(NCMM.GE.4)) THEN
c... Then, consider C5/C10adj + C12adj for MMLR(m)={5,6,8,10,12,14} cases
              CmEFF(4)= CMM(4) + 0.25D0*CMM(1)**2/DSCM
              WRITE(6,710) MMLR(4),MMLR(4),CmEFF(4)
              IF(NCMM.GE.5) THEN                   ! introduce C12^{adj}
                  CmEFF(5)= CMM(5) + 0.25D0*CMM(2)**2/DSCM
                  WRITE(6,710) MMLR(5),MMLR(5),CmEFF(5)
                  IF(NCMM.GE.6) THEN               ! introduce C14^{adj}
                      CmEFF(6)= CMM(6) + 0.5D0*CMM(2)*CMM(3)/DSCM
                      WRITE(6,710) MMLR(6),MMLR(6),CmEFF(6)
                      ENDIF
                  ENDIF
              ENDIF
          IF((MMLR(1).EQ.4).AND.(NCMM.GE.3)) THEN
c... Then, consider C4/C8adj + C12adj for MMLR(m)={4,6,7,8,10,12,14} cases
c... first, allowing for a C7:   C8 is m=4, C10 is m=5 ... etc.
              IF((MMLR(3).EQ.7).AND.(NCMM.GE.4)) THEN
                  CmEFF(4)= CMM(4) + 0.25D0*CMM(1)**2/DSCM
                  WRITE(6,712) MMLR(4),MMLR(4),CmEFF(4)
                  IF(NCMM.GE.5) THEN             ! implicitly MMLR(5)=10
                      CmEFF(5)= CMM(5) + 0.5D0*CMM(1)*CMM(2)/DSCM
                      WRITE(6,710)MMLR(5),MMLR(5), CmEFF(5)
                      IF(NCMM.GE.6) THEN         ! implicitly MMLR(6)=12
                          CmEFF(6)= CMM(6) + 0.5D0*CMM(1)*CmEFF(4)/DSCM
     1                                         + 0.25D0*CMM(3)**2/DSCM
                          WRITE(6,710) MMLR(6),MMLR(6),CmEFF(6)
                          IF(NCMM.GE.7) THEN     ! implicitly MMLR(7)=14
                              CmEFF(7)= CMM(7) + 0.25d0*CMM(3)**2/DSCM
     2           + 0.5D0*CMM(2)*CMM(4)/DSCM + 0.5D0*CMM(1)*CMM(5)/DSCM
                              WRITE(6,710) MMLR(7),MMLR(7),CmEFF(7)
                              ENDIF
                          ENDIF
                      ENDIF
                ELSE
c... Now, consider C4/C8adj + C12adj for MMLR(m)={4,6,8,10,12,14} cases
                  CmEFF(3)= CMM(3) + 0.25D0*CMM(1)**2/DSCM
                  WRITE(6,712) MMLR(3),MMLR(3),CmEFF(3)
                  IF(NCMM.GE.4) THEN           ! implicitly MMLR(4)=10
                      CmEFF(4)= CMM(4) + 0.5D0*CMM(1)*CMM(2)/DSCM
                      WRITE(6,710) MMLR(4),MMLR(4),CmEFF(4)
                      IF(NCMM.GE.5) THEN       ! implicitly MMLR(5)=12
                          CmEFF(5)= CMM(5) + 0.5D0*CMM(1)*CMM(3)/DSCM
     1                                         + 0.25D0*CMM(2)**2/DSCM
                          WRITE(6,710) MMLR(5),MMLR(5),CmEFF(5)
                          IF(NCMM.GE.6) THEN     ! implicitly MMLR(6)=14
                              CmEFF(6)= CMM(6)+0.5D0*CMM(2)*CMM(3)/DSCM
     1                                      + 0.5D0*CMM(1)*CMM(4)/DSCM
                              WRITE(6,710) MMLR(6),MMLR(6),CmEFF(6)
                              ENDIF
                          ENDIF
                      ENDIF
                ENDIF
              ENDIF
          IF((MMLR(1).EQ.3).AND.(NCMM.GE.2)) THEN
c... Then, consider C3/C6adj & C9adj for MMLR(m)={3,6,8,(9),10,12,14} cases
              CmEFF(2)= CMM(2) + 0.25D0*CMM(1)**2/DSCM
              WRITE(6,712) MMLR(2),MMLR(2),CmEFF(2)
              IF(NCMM.GE.3) THEN            ! introduce C9adj & MMLR=9
                  MCMM= NCMM+1              ! placing m=9 as last power
                  MMLR(MCMM)= 9
                  CmEFF(MCMM)= 0.5d0*CMM(1)*CmEFF(2)/DSCM
                  WRITE(6,714) MMLR(MCMM),CmEFF(MCMM)
                  IF(NCMM.GE.5) THEN           ! implicitly MMLR(5)=12
                      CmEFF(5)= CMM(5) + 0.5D0*CMM(1)*CmEFF(MCMM)/DSCM
     1                                       + 0.25D0*CmEFF(2)**2/DSCM
                      WRITE(6,710) MMLR(5),MMLR(5),CmEFF(5)
                      IF(NCMM.GE.6) THEN        ! implicitly MMLR(6)=14
                          CmEFF(6)= CMM(6)+0.5D0*CmEFF(2)*CMM(3)/DSCM
                          WRITE(6,710) MMLR(6),MMLR(6),CmEFF(6)
                          ENDIF
                      ENDIF
                  ENDIF
              ENDIF
c** End of  CmEFF= Cm + CmADJ  setup ===================================
c** Now - initialize for Damping at r= REQ

          CALL dampF(REQ,rhoAB,MCMM,MMLR,IVSR,IDSTT,DM,DMP,DMPP)
          DO  m= 1,MCMM
              dLULRedCm(m)= 1.d0/REQ**MMLR(m)
              IF(rhoAB.GT.0.d0) dLULRedCm(m)= Dm(m)*dLULRedCm(m)
              T0= CmEFF(m)*dLULRedCm(m)
              ULRe= ULRe + T0
              T1= T1 + MMLR(m)*T0
              ENDDO
          dLULRedRe= -T1/(ULRe*REQ)
          DO  m= 1,MCMM
              dLULRedCm(m)= dLULRedCm(m)/ULRe
              IF(rhoAB.GT.0) dLULRedRe= dLULRedRe + 
     1                                       dLULRedCm(m)*Dmp(m)/Dm(m)
              ENDDO
          betaINF= DLOG(2.d0*DSCM/ULRe)
c*** print for MLR form
          WRITE(6,602) PPAR,QPAR,DSCM,REQ
c... for Huang form: \beta(yp)= Binf*yp + [1-yp]*{power series in yq}
          WRITE(6,607) PPAR,PPAR,QPAR,Nbeta,Nbeta+1,(BETA(J),J= 0,Nbeta)
          IF(Rref.GT.0) THEN
              WRITE(6,613) Rref
            ELSE
              WRITE(6,615) REQ
              Rref= REQ
            ENDIF  
          IF(rhoAB.GT.0.d0) THEN    !! describe type of Damping Fx
              PVSR= 0.5d0*IVSR
              IF(IDSTT.GT.0) THEN     !! first option: Douketis-Scoles
                  WRITE(6,664) rhoAB,PVSR,bDS(IVSR),cDS(IVSR),PVSR
                ELSE                  !! second option: Tang-Toennies
                  LVSR= IVSR/2
                  WRITE(6,666) rhoAB,LVSR,bTT(LVSR)
                ENDIF
            ELSE                     !!  ELSE ... for no damping ...
              WRITE(6,668)
            ENDIF
          WRITE(6,617) betaINF,MMLR(1),CmEFF(1),MMLR(1)
          IF(NCMM.GT.1) THEN
              DO  I= 2,NCMM
                  IF(MMLR(I).LE.9) WRITE(6,619) MMLR(I),CMM(I)
     1                                                    ,MMLR(I)
                  IF(MMLR(I).GT.9) WRITE(6,621) MMLR(I),CMM(I)
     1                                                    ,MMLR(I)
                  ENDDO
              ENDIF
c====================================================================
  602 FORMAT(/' MLR(p=',I1,', q=',I1,') Potential with:   De='
     1 ,F10.4,'[cm-1]    Re=',F12.8,'[A]')
  607 FORMAT('   with exponent coefficient   beta(r)= beta{INF}*y',I1,
     1  ' + [1-y',i1,']*Sum{beta_i*y',i1,'^i}'/6x,'exponent coefft. powe
     2r series order',I3/6x,'and',i3,' coefficients:',1PD16.8,2D16.8:/
     3  (10x,4D16.8:))
  613 FORMAT(9x,'with radial variables  y_p and/or y_q  defined w.r.t.',
     1 '  Rref=',F10.7)
  615 FORMAT(6x,'radial variables  y_p and/or y_q  defined w.r.t.',
     1 '  Rref= Re=' F10.7)
  617 FORMAT('      betaINF=',f16.12,'  & uLR defined by  C',i1,' =',
     1  1PD13.6,'[cm-1 Ang','^',0P,I1,']')
  619 FORMAT(50x,'C',I1,' =',1PD13.6,'[cm-1 Ang','^',0P,I1,']')
  621 FORMAT(50x,'C',I2,'=',1PD13.6,'[cm-1 Ang','^',0P,I2,']')
  664 FORMAT(4x,'uLR inverse-power terms incorporate DS-type damping wit
     1h   rhoAB=',f10.7/8x,'defined to give very short-range  Dm(r)*Cm/r
     2^m  behaviour   r^{',SS,f4.1,'}'/8x,'Dm(r)= [1 - exp(-',f5.2, 
     3 '(rhoAB*r)/m -',f6.3,'(rhoAB*r)^2/sqrt{m})]^{m',SP,F4.1,'}')
  666 FORMAT(4x,'uLR inverse-power terms incorporate TT-type damping wit
     1h   rhoAB=',f10.7/8x,'defined to give very short-range  Dm(r)*Cm/r
     2^m  behaviour   r^{',I2,'}'/8x,'Dm(r)= [1 - exp(-bTT*r)*SUM{(bTT*r
     3)^k/k!}]   where   bTT=',f6.3,'*rhoAB')
  668 FORMAT(4x,'uLR inverse-power terms incorporate NO damping function
     1s')
  710 Format("  'Dattani adjustment' for MLR  C",I2,'  yields   internal
     1  C',  I2,'{adj}=',1PD14.7)
  712 Format("  'Dattani adjustment' for MLR  C",I1, '   yields internal
     1   C',  I1,'{adj}=',1PD15.8)
  714 Format("  'Dattani adjustment' for MLR(m1=3) yields internal   C",
     1  I1,'{adj}=',1PD15.8)
c** End of initialization & printout
          REP= REQ**PPAR
          IF(Rref.LE.0) AREF= REQ
          IF(Rref.GT.0) AREF= Rref
          AREFp= AREF**PPAR
          AREFq= AREF**QPAR
          ENDIF
c=======================================================================
      VDIST= 0.0d0
      BETADIST= 0.d0
c** Initialize parameter counter for this state ...
      IPVSTART= 0
      RVAL= RDIST
      RINV= 1.d0/RVAL
      RDp= RVAL**PPAR
      RDq= RVAL**QPAR
      YPE= (RDp-REP)/(RDp+REP)
      YP= (RDp-AREFp)/(RDp+AREFp)
      YQ= (RDq-AREFq)/(RDq+AREFq)
      YPM= 1.d0 - YP
      DYPDRE= -0.5d0*PPAR*(1.d0 - YP**2)/REQ
      DYQDRE= -0.5d0*QPAR*(1.d0 - YQ**2)/REQ
      DYPEDRE= -0.5d0*PPAR*(1.d0 - YPE**2)/REQ
      DYPDR= -DYPDRE*REQ*RINV
      DYPEDR= 0.5d0*PPAR*RINV*(1.d0 - YPE**2)
      DYQDR= -DYQDRE*REQ*RINV
      D2YPDR= -DYPDR*RINV*(1.d0 + PPAR*YP)
      D2YPEDR= -DYPEDR*RINV*(1.d0 + PPAR*YPE)
      D2YQDR= -DYQDR*RINV*(1.d0 + QPAR*YQ)
      DYPpDRE= -PPAR*YP*RINV*DYPDRE
      DYPEpDRE= -PPAR*YPE*RINV*DYPEDRE
      DYQpDRE= -QPAR*YQ*RINV*DYQDRE
      D2VAL= 0.d0
      YPP= 1.d0
      DVAL= 0.d0
      DBDB(0)= 1.0d0
      VAL= BETA(0) + YQ*BETA(1)
      DVAL= BETA(1)
      Nbeta= Nbeta
c-------------------------------------------------------------------
      DO  J= 2,Nbeta
c... now calculate power series part of the Morse-like exponent,along
c    with its radial derivatives
          D2VAL= D2VAL + BETA(J)* DBLE(J)*DBLE(J - 1) *YPP
          YPP= YPP*YQ
          DVAL= DVAL + BETA(J)* DBLE(J)* YPP
          YPPP= YPP* YQ
          VAL= VAL + BETA(J)*YPPP
          DBDB(J)= YPM*YPPP
          ENDDO
      YPP= YPPP
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c*** DBDB & DBDRe= dBeta/dRe  used in uncertainty calculation in fununc.f
      DBDRe= -YP*dLULRedRe
      dVALdRe= DBDRe + (betaINF - VAL)*DYPDRE + (1.d0 - YP)*DVAL*DYQDRE
      IF(Rref.LE.0.d0) DBDRe= dVALdRe
c-----------------------------------------------------------------------
c... now the power series and its radial derivatives are used in the
c    construction of the derivatives  with respect to the parameters
      dBpdRe= DYPpDRE*(betaINF - VAL) - DYPDR*dLULRedRe
     1   + (-DYPDR*DYQDRE + (1.d0 - YP)*DYQpDRE - DYPDRE*DYQDR)*DVAL
     2                              + (1.d0 - YP)*DYQDR*DYQDRE*D2VAL
      D2VAL= (betaINF - VAL)*D2YPDR - 2.d0*DYPDR*DYQDR*DVAL
     1                   + (1.d0- YP)*(D2YQDR*DVAL + DYQDR**2*D2VAL)
      DVAL= (betaINF - VAL)*DYPDR + (1.d0- YP)*DYQDR*DVAL
      VAL= YP*betaINF + (1.d0- YP)*VAL
      dBdR= dYPEdR*VAL + YPE*DVAL
      d2BdR= d2YPEdR*VAL + 2.d0*dYPEdR*DVAL + YPE*D2VAL
      dYBdRe= DYPEDRE*VAL + YPE*dVALdRe
      dYBpdRe= VAL*DYPEpDRE + DYPEDRE*DVAL + DYPEDR*dVALdRe+ YPE*dBpdRe
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      XTEMP= DEXP(-VAL*YPE)
c** Now begin by generating  uLR(r)
      ULR= 0.d0
c-------------------------------------------------------------------
      dULRdR= 0.d0
      d2ULRdR= 0.d0
      dULRdRCm= 0.d0
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IF(rhoAB.GT.0.d0) THEN
          CALL dampF(RVAL,rhoAB,NCMM,MMLR,IVSR,IDSTT,Dm,Dmp,Dmpp)
          ENDIF
      DO  m= 1,NCMM
          IF(rhoAB.LE.0.d0) THEN
c-----------------------------------------------------------------------
              dULRdCm(m)= 1.d0*RINV**MMLR(m)
              dULRdRCm(m)= -dULRdCm(m)*RINV*DBLE(MMLR(m))
              dDULRdR= 0.d0
              d2DULRdR= 0.d0
            ELSE
              dULRdCm(m)= Dm(m)*RINV**MMLR(m)
              dULRdRCm(m)= -dULRdCm(m)*RINV *DBLE(MMLR(m)) 
     2                                   + Dmp(m)*RINV**MMLR(m)
              dDULRdR= Dmp(m)*RINV**MMLR(m)
              d2DULRdR= Dmpp(m)*RINV**MMLR(m)
            ENDIF
          ULR= ULR + CmEFF(m)*dULRdCm(m)
          dULRdR= dULRdR + CmEFF(m)*(dDULRdR
     1                         - dULRdCm(m)*RINV*DBLE(MMLR(m)))
          d2ULRdR= d2ULRdR + CmEFF(m)*(d2DULRdR
     1   - 2.d0*dDULRdR*RINV*DBLE(MMLR(m)) + dULRdCm(m)*RINV**2
     2               *DBLE(MMLR(m))*DBLE((MMLR(m) + 1)))
          ENDDO
      dLULRdR= dULRdR/ULR
      DO m= 1,NCMM
          dLULRdCm(m)= dULRdCm(m)/ULR
          ENDDO
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      XTEMP= XTEMP*ULR/ULRe
c... note ... reference energy for each state is asymptote ...
      DVDD= XTEMP*(XTEMP - 2.D0)  
      VDIST= DSCM*DVDD 
      BETADIST= VAL
      YPP= 2.d0*DSCM*(1.0d0-XTEMP)*XTEMP
      IPV= IPVSTART+2
c... derivatives w.r.t R
      DXTEMP= XTEMP*(dLULRdR - dBdR)
      D2XTEMP= XTEMP*(dBdR**2 - d2BdR + (d2ULRdR - 2*dBdR*dULRdR)/ULR)
      dVdR= 2.d0*DSCM*DXTEMP*(XTEMP - 1.d0)
      d2VdR2= 2.d0*DSCM*(DXTEMP**2 + D2XTEMP*(XTEMP - 1.d0))
c... derivative w.r.t. De
      IPV= IPVSTART+1
      dXdP(IPV)= -XTEMP*YPE*YP
      dXpdP(IPV)= -(XTEMP*(YPE*DYPDR + DYPEDR*YP)+ YPE*YP*DXTEMP)
      DVDD= DVDD + YPP*YP*YPE/DSCM
      YPP= YPP*YPE*(1.d0 - YP)
      dVpdP(IPV)= 2.d0*(dXdP(IPV)*DXTEMP + (XTEMP - 1.d0)*dXpdP(IPV)) 
     1                                    + 2.d0*(XTEMP - 1.d0)*DXTEMP
      DVtot(IPV)= DVDD
c... derivative w.r.t. Re
      IPV= IPV+1
      dXdP(IPV)= -XTEMP*(dYBdRe + dLULRedRe)
      dXpdP(IPV)= -DXTEMP*(dYBdRe + dLULRedRe) - XTEMP*dYBpdRe
      dVpdP(IPV)= 2.d0*DSCM*(dXdP(IPV)*DXTEMP 
     1                                    + (XTEMP - 1.d0)*dXpdP(IPV))
      DVtot(IPV)= YPP*(dYBdRe + dLULRedRe)
c ... derivative w.r.t. Cm's
      DO  m= 1, NCMM
          IPV= IPV+1
          dXdP(IPV)= XTEMP*(dLULRdCm(m) + (YPE*YP - 1.d0)*dLULRedCm(m))
          dXpdP(IPV)= DEXP(-VAL*YPE)/ULRe*(dULRdRCm(m)- dBdR*dULRdCm(m))
     1                    + (DXTEMP*(YPE*YP - 1.d0) + XTEMP*(dYPEdR*YP
     2                                      + YPE*dYPdR))*dLULRedCm(m)
          dVpdP(IPV)= 2.d0*DSCM*(dXdP(IPV)*DXTEMP
     1                                    + (XTEMP - 1.d0)*dXpdP(IPV))
          ENDDO
c... finally ... derivatives w.r.t. exponent expansion coefficients

      DO  J= 0,Nbeta
          IPV= IPV+1
          dXdP(IPV)= XTEMP*YPE*(1.d0 - YP)*YQ**J
          dXpdP(IPV)= (XTEMP*((1.d0 - YP)*DYPEDR - DYPDR*YQ)
     1                  + YPE*(1.d0 - YP)*DXTEMP)*YQ**J + XTEMP*J*(YPE
     2                                       *(1.d0 - YP))*YQ**(J - 1)
          dVpdP(IPV)= 2.d0*DSCM*(dXdP(IPV)*DXTEMP
     1                                    + (XTEMP - 1.d0)*dXpdP(IPV))
          DVtot(IPV)= YPP
          YPP= YPP*YQ
          ENDDO
ccccc Print for testing
cc    rewind(10)
cc    write(10,610) (RD(i),vpot(i,istate),BETAFX(i,istate),
cc   1                                        i= 1, NDATPT,50)
cc610 FORMAT(/(f10.4,f15.5,f12.6))
ccccc End of Print for testing
c======================================================================

c**************************** BLOCK END ******************************72
  999 RETURN
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE dampF(r,rhoAB,NCMM,MMLR,IVSR,IDSTT,DM,DMP,DMPP)
c** Subroutine to generate values 'Dm' and its first `Dmp' and second
c   'Dmpp' derivatives w.r.t. r of the chosen form of the damping
c    function, for  m= 1 to MMAX.
c---------------------- RJL Version of 05 August 2014 ------------------
c-----------------------------------------------------------------------
c                 Upon Input
c* r - the radial distance in Angsroms (!) 
c* RHOab  'universal' scaling coefficient used for systems other than H_2
c       RHOab= 2*(RHOa*RHOb)/(RHOa+RHOb) where RHOa = (I_p^A/I_p^H)^0.66
c              where I_p^A is the ionization potential of atom A
c              and I_p^H is the ionization potential of atomic hydrogen
c* NCMM  the number of inverse-power terms to be considered
c* MMLR  are the powers of the NCMM inverse-power terms
c* IVSR defines damping s.th.  Dm(r)/r^m --> r^{IVSR/2} as r --> 0
c* IDSTT specifies damping function type:  > 0  use Douketis et al. form 
c                               if  IDSTT .LE. 0  use Tang-Toennies form
c-----------------------------------------------------------------------
c                 Upon Output
c  DM(m) - The value of the damping function for the long range term 
c          C_MMLR(m)/r^MMLR(m)    {m= 1, NCMM}
c  DMP(m): first derivative of the damping function  DM(m) w.r.t. r
c  DMPP(m): second derivative of the damping function  DM(m) w.r.t. r
c-----------------------------------------------------------------------
      INTEGER NCMM,NCMMax,MMLR(NCMM),IVSR,IDSTT,IVSRF,FIRST, Lsr,m,
     1  MM,MMAX,MMTEMP
      REAL*8 r,rhoAB,bTT(-2:2),cDS(-4:4),bDS(-4:4),aTT,br,XP,YP,
     1  TK, DM(NCMM),DMP(NCMM),DMPP(NCMM),SM(-3:25),
     2  bpm(20,-4:0), cpm(20,-4:0),ZK
c------------------------------------------------------------------------
c  The following values for the numerical factors used in both TT and DS
c  were  normalized to the Hydrogen data presented
c  by Kreek and Meath in J.Chem.Phys. 50, 2289 (1969).
c  The ratio has been chosen such that  b= FACTOR*(I_p^X / I_p^H)^{2/3}
c  for the homoatomic diatomic species X_2, where I_p^A is the ionization
c------------------------------------------------------------------------
       DATA bTT/2.10d0,2.44d0,2.78d0,3.13d0,3.47d0/
       DATA bDS/2.50d0,2.90d0,3.30d0,3.69d0,3.95d0,0.d0,4.53d0,0.d0,
     1          4.99d0/
       DATA cDS/0.468d0,0.446d0,0.423d0,0.405d0,0.390d0,0.d0,0.360d0,
     1            0.d0,0.340d0/
       DATA FIRST/ 1/
       SAVE FIRST, bpm, cpm
c------------------------------------------------------------------------
       MMTEMP = MMLR(1)              ! temporary substn.: corrected@end
       IF(MMLR(1).LE.0) MMLR(1) = 1  ! enforce stability for A-F cases
      IF(RHOab.LE.0) THEN
          DO  m=1,NCMMax
              DM(m)=1.d0
              DMP(m)= 0.d0
              DMPP(m)= 0.d0
              ENDDO
          WRITE(6,602) RHOab
          RETURN
          ENDIF
      IF(IDSTT.LE.0) THEN
c===========================================
c** For Tang-Toennies type damping functions
c===========================================
          Lsr= IVSR/2
          IF((IVSR.LT.-4).OR.(IVSR.GT.4).OR.((2*LSR).NE.IVSR)) THEN
                WRITE(6,600) 'TT',IVSR
                STOP
                ENDIF
          MMAX= MMLR(NCMM) + Lsr - 1
          aTT= RHOab*bTT(Lsr)
          br= aTT*r
          XP= DEXP(-br)
          SM(-3)= 0.d0
          SM(-2)= 0.d0
          SM(-1)= 0.d0
          SM(0)=  1.d0
          TK= 1.d0
          IF(br.GT.0.5d0) THEN
              DO  m= 1,MMAX
                  TK= TK*br/DFLOAT(m)
                  SM(m)= SM(m-1)+ TK
                  ENDDO
              DO m= 1, NCMM
                  MM= MMLR(m) - 1 + Lsr
                  DM(m)= 1.d0 - XP*SM(MM)
                  DMP(m)= aTT*XP*(SM(MM) - SM(MM-1))
                  DMPP(m)= -aTT*aTT*XP*(SM(MM) 
     1                                     - 2.d0*SM(MM-1) + SM(MM-2))
                  ENDDO
c-----------------------------------------------------------------------
c  The above section handles the calculation of the value of the damping
c  function for most values of r.  However, at very small r that algorithm
c  becomes unstable due to numerical noise.  To avoid this, if the 
c  argument is very small it is re-evaluated as a finite sum ...
c-----------------------------------------------------------------------
            ELSE
              MMAX= MMAX+5
              DO  m= 1, MMAX
c... NOTE that here SM(m) is the m'th term  (b*r)^m/m!  [not a sum]
                  SM(m)= SM(m-1)*br/DFLOAT(m)
                  ENDDO
              DO  m= 1, NCMM
                  MM= MMLR(m) + Lsr
                  DM(m)= XP*(SM(MM)+ SM(MM+1)+ SM(MM+2)+ SM(MM+3) 
     1                                                     + SM(MM+4))
                  DMP(m)= aTT*XP*SM(m-1)
                  DMPP(m)= aTT*aTT*XP*(SM(m-2)-SM(m-1))
                  ENDDO
            ENDIF
          ENDIF
c
      IF(IDSTT.GT.0) THEN
c=======================================================================
c** For Douketis-Scoles-Marchetti-Zen-Thakkar type damping function ...
c=======================================================================
          IF((IVSR.LT.-4).OR.(IVSR.GT.4).OR.(IVSR.EQ.1).OR.(IVSR.EQ.3)) 
     1    THEN
              WRITE(6,600) 'DS',IVSR
              STOP
              ENDIF
          IF(FIRST.EQ.1) THEN
              DO m= 1, 20
                  DO  IVSRF= -4,0
                      bpm(m,IVSRF)= bDS(IVSRF)/DFLOAT(m)
                      cpm(m,IVSRF)= cDS(IVSRF)/DSQRT(DFLOAT(m))
                      ENDDO
                  ENDDO
              FIRST= 0 
              ENDIF
          br= rhoAB*r
          DO m= 1, NCMM
              MM= MMLR(m)
              XP= DEXP(-(bpm(MM,IVSR) + cpm(MM,IVSR)*br)*br)
              YP= 1.d0 - XP
              ZK= MM + 0.5d0*IVSR
              DM(m)= YP**ZK
              TK= (bpm(MM,IVSR) + 2.d0*cpm(MM,IVSR)*br)*rhoAB
              DMP(m) = ZK*XP*TK*DM(m)/YP
c ... calculate second derivative [for DELR case] {check this!}
              DMPP(m)= (ZK-1.d0)*DMP(m)*(XP*TK)/YP
     1               - DMP(m)*TK + DMP(m)*2.d0*cpm(MM,IVSR)*rhoAB**2/TK
              ENDDO   
          ENDIF  
      MMLR(1) = MMTEMP
      RETURN
  600 FORMAT(/,' *** ERROR ***  For  ',A2,'-damping functions not yet de
     1fined for   IVSR=',i3)
  602 FORMAT( /,' ***ERROR ***  should not call dampF when rhoAB=',
     1  F7.4)
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
