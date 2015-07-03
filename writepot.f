c***********************************************************************
      SUBROUTINE WRITEPOT(NPASS,SLABL,NAME,DECM,PV,PU,PS,CM,VMAXIN, 
     1                    BANDNAME)
c***********************************************************************
c** Subroutine to print out complete description of the potential fx.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c               ----- Version of  29 June 2015 -----
c          (after allowing p,q,r^p_ref,r^q_ref to be free)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** On entry:
c    NPASS   is the number of times WRITEPOT was called.
c    NAME    is the name of the molecule.
c    PU      are the parameter uncertainties.
c    PS      are the parameter sensitivities.
c----------------------------------------------------
c    NSTATES is number of states being considered (in COMMON BLKPARAM)
c=======================================================================
      INCLUDE 'arrsizes.h'
      INCLUDE 'BLKPOT.h'
      INCLUDE 'BLKPARAM.h'
      INCLUDE 'BLKISOT.h'
      INCLUDE 'BLKBOB.h'
      INCLUDE 'BLKCOUNT.h'
c-----------------------------------------------------------------------
c** Common block for partial derivatives of potential at the one distance RDIST
c   and HPP derivatives for uncertainties
      REAL*8 dVdPk(HPARMX),dDe(0:NbetaMX),dDedRe
      COMMON /dVdPkBLK/dVdPk,dDe,dDedRe
c=======================================================================
      INTEGER MMAX, VTST, IISTP
      PARAMETER (MMAX= 20)
      CHARACTER*2  NAME(2),LAB4(-4:0)
      CHARACTER*3  SLABL(-5:NSTATEMX)
      CHARACTER*5 DASH
      CHARACTER*6  BCNAM(8),QCNAM(8),LAB2(7),LAB3(10)
      CHARACTER*7  NAMEDBLE(2)
      CHARACTER*30 BANDNAME(NPARMX)
      INTEGER NPASS, NALL,ISTATE,IPV,I,I4,ISOT,J,MMN,m,m1,LSR,MMp2,
     1  JSTATE,NUApr,NUBpr,NTApr,NTBpr,IPVRe(NSTATEMX),MMLR1D(NCMMAX),
     2  VMAXIN(NSTATEMX),IPVRrefP(NSTATEMX),IPVRrefQ(NSTATEMX)
      REAL*8 DECM(NSTATEMX),BTEMP,UAT,UBT,SAT,BINF,RE3,RE6,RE8,T0,T1,
     1 DX,DX1,ULRe,C3VAL,C6adj,C9adj,RET,RETSig,RETPi,RETp,RETm,Tm,
     2 VATTRe,dVATTRe,PVSR,ATT,YP,YPP,BT,Rinn,Rout,A1,A2,A3,B5,VX,dVX,
     3 uLR,CMMp2,uCMMp2,uDe,XRO,dXRO,dXROdRe,d2XROdRe,fRO,XROpw,ROmp2,
     4 dCmp2dRe,dDeROdRe,yqRe,betaRe,yPOW,XRI,AREFP,AREFQ,
     5 dCmp2(0:NbetaMX),dULRdCm(NCMMax),DM(MMAX),DMP(MMAX),DMPP(MMAX),
     6 bTT(-2:2),cDS(-4:4),bDS(-4:4)
      DATA QCNAM/' QB(v=','QD(v=',' QH(v=',' QL(v=',' QM(v=',
     1    ' QN(v=',' QO(v=',' '/
      DATA BCNAM/' Tv(v=',' Bv(v=','-Dv(v=',' Hv(v=',' Lv(v=',' Mv(v=',
     1    ' Nv(v=',' Ov(v='/
c** Damping function factors from Table 1 of Mol.Phys. 109, 435 (2011)
       DATA bTT/2.10d0,2.44d0,2.78d0,3.13d0,3.47d0/
       DATA bDS/2.50d0,2.90d0,3.3d0,3.69d0,3.95d0,0.d0,4.53d0,0.d0,
     1          4.99d0/
cc     DATA cDS/0.468d0,0.446d0,0.423d0,0.405d0,0.390d0,0.d0,0.360d0,
       DATA cDS/0.468d0,0.446d0,0.423d0,0.400d0,0.390d0,0.d0,0.360d0,
     1            0.d0,0.340d0/
      DATA DASH/'-----'/ 
      SAVE bTT, bDS, cDS
c c** NLLSSRR variables used in output c
      REAL*8 PV(NPARMX), PU(NPARMX), PS(NPARMX),CM(NPARMX,NPARMX),
     1                                                      PT(HPARMX)
      DATA NAMEDBLE/'wLambda',' wSigma'/
c** labels for matrix cases of MLR: LAB2 for 2x2, LAB3 for 3x3, LAB4 for names
      DATA LAB2/'DELTAE',' C3Sig','  C3Pi',' C6Sig','  C6Pi',
     1                   ' C8Sig','  C8Pi'/
      DATA LAB3/'DELTAE',' C3Sig',' C3Pi1',' C3Pi3',' C6Sig',' C6Pi1',
     1          ' C6Pi3',' C8Sig',' C8Pi1',' C8Pi3'/
      DATA LAB4/' x',' B',' c',' b',' A'/
c-----------------------------------------------------------------------
c** Writing out state specific information.
c-----------------------------------------------------------------------
      IPV= 0 
      DO  90 ISTATE=1,NSTATES
          VATTRe= 0.d0
          dVATTRe= 0.d0
          WRITE(6,600)
          IF(PSEL(ISTATE).LT.0) THEN
c** Write .20 file heading for term-value or band-constant states
              IF(NPASS.GT.1) THEN
                  WRITE(20,*)
                  WRITE(20,700) SLABL(ISTATE), IOMEG(ISTATE), 
     1                   VMIN(ISTATE,1), VMAX(ISTATE,1),JTRUNC(ISTATE), 
     2                                             EFSEL(ISTATE),ISTATE
                  WRITE(20,701) PSEL(ISTATE),VLIM(ISTATE),
     1                       MAXMIN(ISTATE),BOBCN(ISTATE),OSEL(ISTATE)
                  ENDIF
              ENDIF 
          IF(PSEL(ISTATE).EQ.0) THEN
              WRITE(6,603) SLABL(ISTATE)
              WRITE(6,636) 'VLIM',VLIM(ISTATE)
              GOTO 90
              ENDIF
          IF(PSEL(ISTATE).EQ.-2) THEN
              WRITE(6,601) SLABL(ISTATE)
              GO TO 90
              ENDIF
          IF(PSEL(ISTATE).EQ.-1) THEN
c** If fitting to band constants for this state ....
              IF(NPASS.GT.1) WRITE(6,684) SLABL(ISTATE)
              IF(NPASS.EQ.1) THEN       !! in first initialising call
                      WRITE(6,6062) SLABL(ISTATE),(I,I=1,NISTP)
                      WRITE(6,6072) (DASH,I=1,NISTP)

                  DO  I= VMIN(ISTATE,1),VMAX(ISTATE,1)
                      DO IISTP= 1,NISTP       !! Check bounds on NBC & NQC
                          IF(NBC(I,IISTP,ISTATE).GT.NBCMX)
     1                                      NBC(I,IISTP,ISTATE)= NBCMX
                          IF(NQC(I,IISTP,ISTATE).GT.NBCMX)
     1                                      NQC(I,IISTP,ISTATE)= NBCMX
                          ENDDO
                      WRITE(6,6082) I,(NBC(I,IISTP,ISTATE),
     1                                                  IISTP= 1,NISTP)
                      IF(IOMEG(ISTATE).GT.0) WRITE(6,6092)
     1                            (NQC(I,IISTP,ISTATE),IISTP= 1,NISTP)
                      ENDDO
                  ENDIF
              IF(NPASS.GT.1) THEN   !! in final call after fit is done
                  WRITE(6,690)
                  DO  ISOT= 1, NISTP
                      DO  I= VMIN(ISTATE,ISOT),VMAX(ISTATE,ISOT)
                          IF(NBC(I,ISOT,ISTATE).GT.0) THEN
                              DO  J= 1,NBC(I,ISOT,ISTATE)
                                  IPV= IPV+1
                                  WRITE(6,686) BCNAM(J),I,ISOT,
     1                                         PV(IPV),PU(IPV),PS(IPV)
                                  ENDDO
                              IF(NQC(I,ISOT,ISTATE).GT.0) THEN
                                  DO  J= 1,NQC(I,ISOT,ISTATE)
                                      IPV= IPV+1   !! Lambda Count
                                      WRITE(6,686) QCNAM(J),I,ISOT,
     1                                         PV(IPV),PU(IPV),PS(IPV)
                                      ENDDO
                                  ENDIF
                              ENDIF
                          ENDDO
                      WRITE(6,688)
                      ENDDO
                  DO I= VMIN(ISTATE,1),VMAX(ISTATE,1)
                      WRITE(20,687) I,(NBC(I,ISOT,ISTATE),ISOT= 1,
     1                                                         NISTP)
                      IF(IOMEG(ISTATE).GT.0) 
     1            WRITE(20,6872) (NQC(I,ISOT,ISTATE),ISOT= 1,NISTP)
                      ENDDO
                  ENDIF
              GOTO 90
              ENDIF
          AREFP= RREFP(ISTATE)
          AREFQ= RREFQ(ISTATE)
          IF(AREFP.LE.0.d0) AREFP= RE(ISTATE)
          IF(AREFQ.LE.0.d0) AREFQ= RE(ISTATE)
          IF(PSEL(ISTATE).EQ.1) THEN
c** Header printout for EMO potential
              WRITE(6,602) SLABL(ISTATE),Nbeta(ISTATE),nPB(ISTATE),
     1                                       nPB(ISTATE),Nbeta(ISTATE)
              IF(RREFP(ISTATE).LE.0.d0) WRITE(6,552) (nQB(ISTATE),i=1,5)
              IF(RREFP(ISTATE).GT.0.d0) WRITE(6,555) AREFP,AREFP
              IF(RREFQ(ISTATE).LE.0.d0) WRITE(6,552) (nQB(ISTATE),i=1,5)
              IF(RREFQ(ISTATE).GT.0.d0) WRITE(6,555) AREFQ,AREFQ
              ENDIF

          IF(PSEL(ISTATE).EQ.2) THEN
c** Header printout for MLR potential
              BINF= betaINF(ISTATE)
              WRITE(6,604) SLABL(ISTATE),nPB(ISTATE),nQB(ISTATE) 
              IF(NSR(ISTATE).GE.0) THEN
                  WRITE(6,605) nPB(ISTATE),nPB(ISTATE),nQB(ISTATE),
     1                                       Nbeta(ISTATE)
                ELSE
                  WRITE(6,680) Nbeta(ISTATE)
                  BETA(Nbeta(ISTATE),ISTATE)= BINF
                ENDIF
              IF(RREFP(ISTATE).LE.0.d0) WRITE(6,552) (nQB(ISTATE),i=1,5)
              IF(RREFP(ISTATE).GT.0.d0) WRITE(6,555) AREFP,AREFP
              IF(RREFQ(ISTATE).LE.0.d0) WRITE(6,552) (nQB(ISTATE),i=1,5)
              IF(RREFQ(ISTATE).GT.0.d0) WRITE(6,555) AREFQ,AREFQ
              ENDIF
c
          IF(PSEL(ISTATE).EQ.3) THEN
c** Header printout for DELR potential form ...
              WRITE(6,612) SLABL(ISTATE),Nbeta(ISTATE),nQB(ISTATE),
     1                           nQB(ISTATE),NSR(ISTATE),Nbeta(ISTATE)
              IF(RREFP(ISTATE).LE.0.d0) WRITE(6,552) (nQB(ISTATE),i=1,5)
              IF(RREFP(ISTATE).GT.0.d0) WRITE(6,555) AREFP,AREFP
              IF(RREFQ(ISTATE).LE.0.d0) WRITE(6,552) (nQB(ISTATE),i=1,5)
              IF(RREFQ(ISTATE).GT.0.d0) WRITE(6,555) AREFQ,AREFQ
              ENDIF
c
          IF(PSEL(ISTATE).EQ.4) THEN
c** Header printout for Tiemann HPP potential ...
c... First, need to define long-and short-range connections .... 
c ... Begin by getting V(r) and V'(r) of polynomial VX at R_i and R_o
              WRITE(6,623) SLABL(ISTATE), BETA(Nbeta(ISTATE)+1,ISTATE),
     1        BETA(Nbeta(ISTATE)+2,ISTATE),BETA(Nbeta(ISTATE)+3,ISTATE),
     2                      (nPB(ISTATE)),BETA(Nbeta(ISTATE)+1, ISTATE)
              BT= BETA(Nbeta(ISTATE)+1, ISTATE)
              Rinn= BETA(Nbeta(ISTATE)+2, ISTATE)
              Rout= BETA(Nbeta(ISTATE)+3, ISTATE)
c** With long-range tail an NCMM-term inverse-power sum, define De and 
c  add 1 more inverse-power term  CMMp2/r**{m_{last}+2}} to ensure 
c  continuity and smoothness at  Rout
              XRO= (Rout - RE(ISTATE))/(Rout+ BT*RE(ISTATE))
              YPP= 1.d0
              VX= 0.d0
              dVX= 0.d0
              DO  J= 1, Nbeta(ISTATE)
                  dVX= dVX+ J*YPP*BETA(J,ISTATE)
                  YPP= YPP*XRO
                  VX= VX+ YPP*BETA(J,ISTATE)
                  ENDDO
              dXRO=(RE(ISTATE)+ BT*RE(ISTATE))/(Rout + BT*RE(ISTATE))**2
c***  dXRO= dX(r)/dr @ r=R_{out}   &  dXRORe= dX(r)/dr_e @ r=R_{out}
              dXROdRe= -dXRO*Rout/RE(ISTATE)
              d2XROdRe = (1.d0 + BT)*(Rout - BT*RE(ISTATE))/
     1                                       (Rout + BT*RE(ISTATE))**3
              dVX= dVX*dXRO
c  VX={polynomial part V_X @ Rout} and dVX is its radial derivative
              uLR= 0.d0                       !!   VLIM(ISTATE)
              CMMp2= 0.d0
              DO  J= 1, NCMM(ISTATE)
                  B5= CmVAL(J,ISTATE)/Rout**MMLR(J,ISTATE)
                  uLR= uLR + B5
                  CMMp2= CMMp2 + MMLR(J,ISTATE)*B5
                  ENDDO
              MMp2= MMLR(NCMM(ISTATE),ISTATE)+2
              fRO= Rout**(MMp2+1)/MMp2            !! factor for derivatives
              CMMp2= (dVX - CMMp2/Rout)*Rout**(MMp2+1)/MMp2
c!!! zero our C5(A) for Mg2 to try to match Knoeckel
cc            IF(ISTATE.EQ.2) CMMp2= 0.d0
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
              A2= A2*(RE(ISTATE)+ BT*RE(ISTATE))/(Rinn+BT*RE(ISTATE))**2
              A2= -A2/A1
c** Extrapolate inwardly with the exponential: B5+ A1*exp(-A2*(R-Rinn))
c... Now ... printout for HPP inward and outward extrapolations              
              WRITE(6,676) Rinn,B5,A1,A2,Rinn,Rout,DE(ISTATE),CMMp2,MMp2
  676 FORMAT(5x,'Extrapolate smoothly inward from   Rinn=',f6.3/10x,
     1 'as  ',F14.7'  + ',1PD14.7,'*exp[-',d14.7,'(r -',0PF6.3,')]'/5x,
     2  'Extrapolate smoothly  outward from   Rout=',F6.2,/15x,
     2  'by  setting   De=',  F14.7,' and adding ',1PD15.7,'/r**',I2)
              IF(NPASS.GT.1) THEN
                  dCmp2dRe= 0.d0
                  dDeROdRe= 0.d0
                  XROpw= 1.d0/XRO**2
                  ROmp2= 1.d0/Rout**MMp2
                  DO J=1, Nbeta(ISTATE)
c... first get outer boundary factors & derivatives w.r.t.  beta_i             
                      dCmp2dRe= dCmp2dRe + J*(J-1)*BETA(J,ISTATE)*XROpw
                      XROpw= XROpw*XRO         !! power now   (J-1)
                      dDeROdRe= dDeROdRe + J*BETA(J,ISTATE)*XROpw
                      dCmp2(J)= J*XROpw*dXRO*fRO
                      dDe(J)= XROpw*XRO + ROmp2*dCmp2(J)   !!  uses power !J
                      ENDDO
                  dCmp2dRe=(dDeROdRe*d2XROdRe+dCmp2dRe*dXRO*dXROdRe)*fRO
                  dDedRe=  dDeROdRe*dXROdRe + ROmp2*dCmp2dRe
cc                  if(npass.gt.1) write(6,500) DE(ISTATE),ddedre,
cc   2                                    (j,dDe(j),j=0,nbeta(istate))
c** NOW ... uncertainties in De
                  J= POTPARI(ISTATE)
                  PT(J)= dDedRe*PU(J)
                  DO I= 0, Nbeta(ISTATE)
                      J= J+1
                      PT(J)= PU(J)*dDe(I)
                      ENDDO
                  CALL MMCALC(POTPARI(ISTATE),POTPARF(ISTATE),PT,CM,uDe)
c** NOW ... uncertainties in CMMp2
                  J= POTPARI(ISTATE)
                  PT(J)= dCmp2dRe*PU(J)
                  DO I= 1, Nbeta(ISTATE)
                      J= J+1
                      PT(J)= PU(J)*dCmp2(I)
                      ENDDO
                  CALL MMCALC(POTPARI(ISTATE),POTPARF(ISTATE),PT,CM,
     1                                                         uCMMp2)
                  ENDIF
              ENDIF
          IF(PSEL(ISTATE).EQ.5) THEN
c** Header printout for Tang-Toennies type potential ...  
c** if getting ATT= AA and D_e from R_e and bTT=BB=\beta(0), 
c  ... first get Vatt and Vatt'
              IF(Nbeta(ISTATE).EQ.0) THEN
                  DO m=1,NCMM(ISTATE)
                      MMLR1D(m)= MMLR(m,ISTATE)
                      ENDDO
c... define TT's 'b' coefficient as BB(ISTATE)                      
                  BB(ISTATE)= beta(0,ISTATE)
                  rhoAB(ISTATE)= BB(ISTATE)/3.13d0
                  CALL DAMPF(RE(ISTATE),rhoAB(ISTATE),NCMM(ISTATE),
     1                  MMLR1D,IVSR(ISTATE),IDSTT(ISTATE),DM,DMP,DMPP)
                  DO m=1,NCMM(ISTATE)
                      Tm= CmVAL(m,ISTATE)/Re(ISTATE)**MMLR1D(m)
                      VATTRe= VATTRe + DM(m)*Tm
                      dVATTRe= dVATTRe + Tm*(DMP(m) 
     1                                  - Dm(m)*MMLR1D(m)/RE(ISTATE))
                      ENDDO
c ... define TT's pre-exponential factor as AA      
                  AA(ISTATE)= -dVATTRe*DEXP(+BB(ISTATE)*RE(ISTATE))/
     1                                                    BB(ISTATE)
                  DE(ISTATE)= VATTRE - AA(ISTATE)*DEXP(-BB(ISTATE)*
     1                                                     RE(ISTATE))
     1                                                       
                  WRITE(6,626) SLABL(ISTATE),RE(ISTATE),BB(ISTATE),
     1    DE(ISTATE),AA(ISTATE),rhoAB(ISTATE),rhoAB(ISTATE)*3.13D0
                  ENDIF
  626 FORMAT(/' State ',A3,' represented by a Tang-Toennies type potenti
     1al'/'  Input values of r_e=',F9.6,'  and  b=',F9.6/10x,'are used t
     1o determine   D_e=',F11.4,'   and A_{TT}=',1PD14.7/10x,'where   rh
     3oAB=',0P,F9.6,'   yields  bTT(eff)=',f10.6)
              ENDIF
          IF(PSEL(ISTATE).EQ.6) THEN
c** For Aziz'ian HFD-type potential: print header and derive leading 
c   exponent coefficient \beta_1  and pre-exponential factor A for use
c   in subroutines 'vgen' & 'vgenp'
              A1= BETA(0,ISTATE)
              A2= BETA(1,ISTATE)
              A3= BETA(2,ISTATE)
              DX= 1.d0
              DX1= 0.d0
              IF(A2.GT.1.d0) THEN
                  DX= DEXP(-A1*(A2- 1.d0)**A3)
                  DX1= A1*A2*A3*DX*(A2- 1.d0)**(A3- 1.d0)
                  ENDIF
              T0= 0.d0
              DO m= 1,NCMM(ISTATE)
                  T0= T0+ CMVAL(m,ISTATE)
                  ENDDO
               T0= T0*DX - 1.d0
               IF(T0.LE.0.d0) THEN
                   WRITE(6,624) T0,(MMLR(m,ISTATE),CmVAL(m,ISTATE),
     1                                             m= 1, NCMM(ISTATE))
                   STOP
                   ENDIF
              T1= 0.d0
              DO  m= 1, NCMM(ISTATE)
                  T1= T1+ (MMLR(m,ISTATE)*DX - DX1)*CMVAL(m,ISTATE)
                  ENDDO
              T1= BETA(4,ISTATE) - 2.d0*BETA(3,ISTATE) + T1/T0
              T0= T0*DEXP(T1+BETA(3,ISTATE))
              BB(ISTATE)= T1
              AA(ISTATE)= T0*DE(ISTATE)
              WRITE(6,628) BETA(4,ISTATE),T1,BETA(3,ISTATE),AA(ISTATE),
     1             (MMLR(m,ISTATE),CmVAL(m,ISTATE),m= 1, NCMM(ISTATE))
              WRITE(6,629) DE(ISTATE),RE(ISTATE),A1,A2,A3
              ENDIF
  624 FORMAT(/' *** ERROR in generating HFD potential *** generate   ALF
     1A=',G15.7,'  from reduced  Cm  coefficients:'/(3x,3(' C',I2,'=',
     2 1PD15.7:)) )
 628  FORMAT(/' Potential is Generalized HFD with exponent factors   gam
     1ma=',f9.6/'   alpha*=',f12.8,'   beta*=',f9.6,' A=',1PD16.9,
     2 " & reduced Cm's:"/(3x,3('   C',I2,' =',D15.7:)) )
 629  FORMAT('   De=',f10.4,'[cm-1]   Re=',f9.6,'[Angst.]   and'/
     1  '     Damping function  D(r)= exp[ -',0P,f6.4,'*(',f7.4,
     2  '/X -1.0)**',f5.2,']')
c-----------------------------------------------------------------------
c** Common uLR(r) printout for the MLR, DELR,  HPP and TT potentials
          IF((PSEL(ISTATE).GE.2).AND.(PSEL(ISTATE).LE.5)) THEN
c... first, specify choice of damping fx. {if damping included}
              IF(rhoAB(ISTATE).GT.0.d0) THEN
                  IF(IDSTT(ISTATE).GT.0) THEN
                      PVSR= DFLOAT(IVSR(ISTATE))*0.5d0
                      WRITE(6,607) rhoAB(ISTATE),PVSR,bDS(IVSR(ISTATE)),
     1                                           cDS(IVSR(ISTATE)),PVSR
                    ELSE
                      LSR= IVSR(ISTATE)/2
                      WRITE(6,617) rhoAB(ISTATE), LSR, bTT(LSR)
                    ENDIF
                ELSE
                  WRITE(6,664)
cc                WRITE(6,608) MMLR(1,ISTATE),CmVAL(1,ISTATE),
cc   1                                                  MMLR(1,ISTATE)
                ENDIF
c** List (inverse) power and coefficients of terms contributing to uLR(r)
c... First ... header (& lead coefft.) for all A-F diagonalization cases
              m1= 1
              I4= 99
              IF(MMLR(1,ISTATE).LE.0) THEN
                  I4= MMLR(1,ISTATE)
                  IF(I4.GE.-1) WRITE(6,606) LAB4(I4),CmVAL(1,ISTATE)
                  IF(I4.LE.-2) WRITE(6,6066) LAB4(I4),CmVAL(1,ISTATE)
                  m1=2
                  ENDIF
              DO  m= m1,NCMM(ISTATE)
              IF((I4.LE.0).AND.(m.GT.1)) THEN           !! A-F cases
c... now,... Cm's for A-F 2x2 cases
                      IF(I4.GE.-1) WRITE(6,708)LAB2(m),CmVAL(m,ISTATE),
     1                                                  MMLR(m,ISTATE)
c... now,... Cm's for A-F 3x3 cases
                      IF(I4.LE.-2) WRITE(6,708) LAB3(m),CmVAL(m,ISTATE),
     1                                                  MMLR(m,ISTATE)
                    ELSE
c... Finally, print Cm's for simple {damped} inverse-power sum cases
                      IF(MMLR(m,ISTATE).LE.9)
     1      WRITE(6,608) MMLR(m,ISTATE),CmVAL(m,ISTATE),MMLR(m,ISTATE)
                      IF(MMLR(m,ISTATE).GT.9)
     1      WRITE(6,609) MMLR(m,ISTATE),CmVAL(m,ISTATE),MMLR(m,ISTATE)
                    ENDIF
                  ENDDO
              IF(PSEL(ISTATE).EQ.2) WRITE(6,682) BINF
              ENDIF
c-----------------------------------------------------------------------

          IF(PSEL(ISTATE).EQ.7) THEN
c** Header printout for Surkus GPEF potential form ...
              WRITE(6,610) SLABL(ISTATE),(nPB(ISTATE),i=1,3),
     1             AGPEF(ISTATE),nPB(ISTATE),BGPEF(ISTATE),nPB(ISTATE)
              ENDIF
c============== End of potential form header printout ==================
          IF((IOMEG(ISTATE).NE.0).AND.(PSEL(ISTATE).GE.0)) 
     1         WRITE(6,683) IOMEG(ISTATE), IOMEG(ISTATE)*IOMEG(ISTATE)
          IF((NUA(ISTATE).GE.0).OR.(NUB(ISTATE).GE.0)) THEN
c** Print description of 'adiabatic' BOB functional forms ...
              IF(BOBCN(ISTATE).GT.0) WRITE(6,556) qAD(ISTATE)
              IF(BOBCN(ISTATE).LE.0)WRITE(6,557) pAD(ISTATE),qAD(ISTATE)
              IF(NUA(ISTATE).GE.0) THEN
                  IF(BOBCN(ISTATE).GT.0) WRITE(6,564) '\tilde{S}(',
     1                               NAME(1),qAD(ISTATE),NUA(ISTATE)-1
                  IF(BOBCN(ISTATE).LE.0) WRITE(6,558) '\tilde{S}(',
     1       NAME(1),pAD(ISTATE),pAD(ISTATE),qAD(ISTATE),NUA(ISTATE)-1
                  WRITE(6,554) NAME(1),(qAD(ISTATE),i= 1,5)
                  ENDIF
              IF(NUB(ISTATE).GE.0) THEN
                  IF(BOBCN(ISTATE).GT.0) WRITE(6,564) '\tilde{S}(',
     1                               NAME(1),qAD(ISTATE),NUB(ISTATE)-1
                  IF(BOBCN(ISTATE).LE.0) WRITE(6,558) '\tilde{S}(',
     1       NAME(1),pAD(ISTATE),pAD(ISTATE),qAD(ISTATE),NUB(ISTATE)-1
                  WRITE(6,554) NAME(2),(qAD(ISTATE),i= 1,5)
                  ENDIF
              ENDIF
          IF((NTA(ISTATE).GE.0).OR.(NTB(ISTATE).GE.0)) THEN
c** Print description of centrifugal BOB functional forms ...
              IF(BOBCN(ISTATE).GT.0) WRITE(6,560) qNA(ISTATE)
              IF(BOBCN(ISTATE).LE.0)WRITE(6,559) pNA(ISTATE),qNA(ISTATE)
              IF(NTA(ISTATE).GE.0) THEN
                  IF(BOBCN(ISTATE).GT.0) WRITE(6,564) '\tilde{R}(',
     1                               NAME(1),qNA(ISTATE),NTA(ISTATE)-1
                  IF(BOBCN(ISTATE).LE.0) WRITE(6,558) '\tilde{R}(',
     1       NAME(1),pNA(ISTATE),pNA(ISTATE),qNA(ISTATE),NTA(ISTATE)-1
                  WRITE(6,554) NAME(1),(qNA(ISTATE),i=1,5)
                      ENDIF
              IF(NTB(ISTATE).GE.0) THEN
                  IF(BOBCN(ISTATE).GT.0) WRITE(6,564) '\tilde{R}(',
     1                               NAME(2),qNA(ISTATE),NTB(ISTATE)-1
                  IF(BOBCN(ISTATE).LE.0) WRITE(6,558) '\tilde{R}(',
     1       NAME(1),pNA(ISTATE),pNA(ISTATE),qNA(ISTATE),NTB(ISTATE)-1
                  WRITE(6,554) NAME(2),(qNA(ISTATE),i=1,5)
                  ENDIF
              ENDIF
          IF((NwCFT(ISTATE).GE.0).AND.(PSEL(ISTATE).GT.0).OR.
     1                                      (PSEL(ISTATE).EQ.-1)) THEN
c** Print description of Lambda/2-Sigma doubling functional forms ...
              IF(IOMEG(ISTATE).GT.0) THEN
                  WRITE(6,618) 'Lambda',Pqw(ISTATE),NwCFT(ISTATE),
     1                                            (Pqw(ISTATE),i= 1,5)
                  IF(efREF(ISTATE).EQ.-1) WRITE(6,692) SLABL(ISTATE)
                  IF(efREF(ISTATE).EQ.0) WRITE(6,694) SLABL(ISTATE)
                  IF(efREF(ISTATE).EQ.1) WRITE(6,696) SLABL(ISTATE)
                  ENDIF
              IF(IOMEG(ISTATE).EQ.-1) THEN
                  WRITE(6,618) ' Gamma',Pqw(ISTATE),NwCFT(ISTATE),
     1                                            (Pqw(ISTATE),i= 1,5)
                  ENDIF
               ENDIF
          IF(IOMEG(ISTATE).LE.-2) WRITE(6,619) -IOMEG(ISTATE)
c c** Write out headings for parameter list
          IF(PSEL(ISTATE).GT.0) THEN
              IF(NPASS.EQ.1) WRITE(6,614)
              IF(NPASS.EQ.2) WRITE(6,615)
              ENDIF
c-----------------------------------------------------------------------
c** Writing out heading for the .20 file
c-----------------------------------------------------------------------
          IF(NPASS.GT.1) THEN
              WRITE(20,*)
              IF(VMAXIN(ISTATE).GE.0) THEN
                  WRITE(20,700) SLABL(ISTATE), IOMEG(ISTATE), 
     1                      VMIN(ISTATE,1), VMAX(ISTATE,1), 
     2                      JTRUNC(ISTATE), EFSEL(ISTATE),ISTATE
                ELSE
                  WRITE(20,700) SLABL(ISTATE), IOMEG(ISTATE),  
     1                      VMIN(ISTATE,1), VMAXIN(ISTATE),      
     2                      JTRUNC(ISTATE), EFSEL(ISTATE),ISTATE
                  WRITE(20,705) (VMAX(ISTATE,I), I=1,NISTP)
                ENDIF 
              WRITE(20,701) PSEL(ISTATE),VLIM(ISTATE),MAXMIN(ISTATE),
     1                                      BOBCN(ISTATE),OSEL(ISTATE)
              WRITE(20,702) RMIN(ISTATE), RMAX(ISTATE), RH(ISTATE)
              IF(PSEL(ISTATE).GE.2) THEN 
                  WRITE(20,*)
                  WRITE(20,703) NCMM(ISTATE), 
     1                        rhoAB(ISTATE),IVSR(ISTATE),IDSTT(ISTATE)
                  IF(NCMM(ISTATE).GT.0) THEN
                      IF(MMLR(1,ISTATE).GT.0) THEN
                          DO I= 1,NCMM(ISTATE)
                              WRITE(20,704) MMLR(I,ISTATE), 
     1                          CmVAL(I,ISTATE),IFXCM(I,ISTATE), I,I,I
                              ENDDO
                        ELSE
                          IF((MMLR(1,ISTATE).EQ.0).OR.
     1                                    (MMLR(1,ISTATE).EQ.-1)) THEN
                              DO I= 1,NCMM(ISTATE)
                                  WRITE(20,706) MMLR(I,ISTATE),
     1                         CmVAL(I,ISTATE),IFXCM(I,ISTATE),LAB2(I)
                                  ENDDO
                            ELSE
                              DO I= 1,NCMM(ISTATE)
                                  WRITE(20,706) MMLR(I,ISTATE),
     1                         CmVAL(I,ISTATE),IFXCM(I,ISTATE),LAB3(I)
                                  ENDDO
                           ENDIF
                        ENDIF
                      ENDIF
                  ENDIF
              ENDIF
          IF(PSEL(ISTATE).EQ.-1) GOTO 90 
c-----------------------------------------------------------------------
c** Writing out the absolute energy information.
c-----------------------------------------------------------------------
          WRITE(6,636) 'VLIM',VLIM(ISTATE)
c-----------------------------------------------------------------------
c** Writing out the Te information.
c-----------------------------------------------------------------------
          IF(ISTATE.GT.1) THEN
              UAT = 0.0d0
              UBT = 0.0d0
              IF(IFXDE(1).LE.0) UAT = PU(1)
              IF(IFXDE(ISTATE).LE.0) UBT = -PU(IPV+1)
cc            UBT = (UAT+UBT)*DSQRT(DECM(ISTATE)+1.0d0)
              UBT = DSQRT(UAT**2 + UBT**2 + 2.d0*DECM(ISTATE)*UAT*UBT)
              UAT = DE(1) - VLIM(1) + VLIM(ISTATE) - DE(ISTATE)
              SAT= DSQRT(PS(1)**2 + PU(IPV+1)**2)
              WRITE(6,620) 'Te',UAT,UBT,SAT
              END IF
c-----------------------------------------------------------------------
c** Writing out the De information.
c-----------------------------------------------------------------------
          IF((PSEL(ISTATE).GE.1).AND.(PSEL(ISTATE).LT.4)) THEN
              IPV= IPV + 1
              IF(IFXDE(ISTATE).LE.0) THEN
                  IF(DABS(DE(ISTATE)).GT.PU(IPV)) THEN
                      WRITE(6,620) 'De',DE(ISTATE),PU(IPV),PS(IPV)
                    ELSE
                      WRITE(6,621) 'De',DE(ISTATE),PU(IPV),PS(IPV)
                    ENDIF
                ELSE
                  WRITE(6,622) 'De',DE(ISTATE)
                ENDIF
              IF(NPASS.GT.1) THEN
                  WRITE(20,670)
                  WRITE(20,670) DE(ISTATE),IFXDE(ISTATE),'De','De'
                  ENDIF
              ENDIF
          IF((PSEL(ISTATE).EQ.4).AND.(NPASS.GT.1))
     1                                WRITE(6,620) 'De',DE(ISTATE),uDe
c-----------------------------------------------------------------------
c** Writing out the Re information.
c-----------------------------------------------------------------------
          IPV= IPV + 1
          IPVRe(ISTATE) = IPV
          IF(IFXRE(ISTATE).LE.0) THEN
              IF (DABS(RE(ISTATE)).GT.PU(IPV)) THEN
                  WRITE(6,620) 'Re',RE(ISTATE),PU(IPV),PS(IPV)
                ELSE
                  WRITE(6,621) 'Re',RE(ISTATE),PU(IPV),PS(IPV)
                ENDIF
            ELSE
              WRITE(6,622) 'Re',RE(ISTATE)
            ENDIF
          IF(NPASS.GT.1) WRITE(20,670)RE(ISTATE),IFXRE(ISTATE),'Re','Re'
c
c-----------------------------------------------------------------------
c** Writing out the RrefP and RrefQ information
c-----------------------------------------------------------------------
          PV= IPV + 1
          IPVRrefP(ISTATE) = IPV
          IF(IFXRREFP(ISTATE).LE.0) THEN
              IF (DABS(RREFP(ISTATE)).GT.PU(IPV)) THEN
                  WRITE(6,620) 'Rp',RREFP(ISTATE),PU(IPV),PS(IPV)
                ELSE
                  WRITE(6,621) 'Rp',RREFP(ISTATE),PU(IPV),PS(IPV)
                ENDIF
            ELSE
              WRITE(6,622) 'Rp',RREFP(ISTATE)
            ENDIF
          IF(NPASS.GT.1) 
     1     WRITE(20,670)RREFP(ISTATE),IFXRREFP(ISTATE),'RRefP','RRefP'
c
          PV= IPV+ 1
          IPVRrefQ(ISTATE) = IPV
          IF(IFXRREFQ(ISTATE).LE.0) THEN
              IF (DABS(RREFQ(ISTATE)).GT.PU(IPV)) THEN
                  WRITE(6,620) 'Rq',RREFQ(ISTATE),PU(IPV),PS(IPV)
                ELSE
                  WRITE(6,621) 'Rq',RREFQ(ISTATE),PU(IPV),PS(IPV)
                ENDIF
              ELSE
c                WRITE(6,222) 'Rq',RREFQ(ISATE) !!don't know why compile
c time error when this is uncommented!
              ENDIF
           IF(NPASS.GT.1) 
     2      WRITE(20,670)RREFQ(ISTATE),IFXRREFQ(ISTATE),'RrefQ','RrefQ'
c
          IF((PSEL(ISTATE).GE.2).AND.(PSEL(ISTATE).LE.5)) THEN
c-----------------------------------------------------------------------
c** For MLR or DELR OR HPP or TT, write out the  Cm  information.
c-----------------------------------------------------------------------
              IF(MMLR(1,ISTATE).LE.0) THEN
c** For Aubert-Frecon treatment of C3(r):C6(r) for alkali dimers
                  DO m=1,NCMM(ISTATE)
                      IPV= IPV+1
                      IF((MMLR(1,ISTATE).EQ.0)
     1                                 .OR.(MMLR(1,ISTATE).EQ.-1))THEN
                           IF(IFXCm(m,ISTATE).LE.0) THEN
                               IF(DABS(CmVAL(m,ISTATE)).GT.PU(IPV))
     1                             WRITE(6,720) LAB2(m),CmVAL(m,ISTATE),
     2                                                   PU(IPV),PS(IPV)
                               IF(DABS(CmVAL(m,ISTATE)).LE.PU(IPV))
     1                             WRITE(6,721) LAB2(m),CmVAL(m,ISTATE),
     2                                                  PU(IPV),PS(IPV)
                           ELSE
                               WRITE(6,722) LAB2(m),CmVAL(m,ISTATE)
                               ENDIF
                        ELSE
                           IF(IFXCm(m,ISTATE).LE.0) THEN
                               IF(DABS(CmVAL(m,ISTATE)).GT.PU(IPV)) THEN
                                   WRITE(6,720) LAB3(m),CmVAL(m,ISTATE),
     1                                                  PU(IPV),PS(IPV)
                                 ELSE
                                   WRITE(6,721) LAB3(m),CmVAL(m,ISTATE),
     1                                                 PU(IPV),PS(IPV)
                                 ENDIF
                             ELSE
                               WRITE(6,722) LAB3(m),CmVAL(m,ISTATE)
                             ENDIF
                        ENDIF
                     ENDDO

                ELSE
c ... For 'regular' MLJ or MLR or DELR or HPP or TT cases ...
                  DO  m= 1,NCMM(ISTATE)
                      IPV= IPV+ 1
                      IF(IFXCm(m,ISTATE).LE.0) THEN
                          IF(DABS(CmVAL(m,ISTATE)).GT.PU(IPV))
     1     WRITE(6,660) MMLR(m,ISTATE),CmVAL(m,ISTATE),PU(IPV),PS(IPV)
                          IF(DABS(CmVAL(m,ISTATE)).LE.PU(IPV))
     1     WRITE(6,661) MMLR(m,ISTATE),CmVAL(m,ISTATE),PU(IPV),PS(IPV)
                        ELSE
                          WRITE(6,662) MMLR(m,ISTATE),CmVAL(m,ISTATE)
                        ENDIF
                      ENDDO
                  IF((PSEL(ISTATE).EQ.4).AND.(NPASS.GT.1)) 
     1                                   WRITE(6,660) MMp2,CMMp2,uCMMp2
                ENDIF
c** Check & do printouts re. Cm values constrained in fits
             IPV= IPVRe(ISTATE)
             DO m= 1, NCMM(ISTATE)
                 IPV= IPV+1
                 IF(IFXCm(m,ISTATE).GT.1) THEN
c... Print re. a fitted Cm value constrained to equal that from another
c   state (with smaller ISTATE).  Input value of IFXCm(m,ISTATE) is IPV
c   parameter-counter value for that earlier Cm value.  c NOTE !!!! Need
c   to fix ISTATE count label  !!!!!!!!!
                     DO JSTATE= ISTATE,1,-1
                         IF((IFXCm(m,ISTATE).LT.IPVRe(JSTATE)).AND.
     1                      (IFXCm(m,ISTATE).GT.IPVRE(JSTATE-1))) THEN
                             CmVAL(m,ISTATE)= CmVAL(m,JSTATE-1)
                             WRITE(6,666) MMLR(m,ISTATE),IPV,
     1                         IFXCm(m,ISTATE),MMLR(m,ISTATE),JSTATE-1,
     2                                                  CmVAL(m,ISTATE)
                             ENDIF
                         ENDDO
                     ENDIF
                 ENDDO
  666 FORMAT('  Constrain C_',I1,' = PV(',i3,')  to equal fitted  PV('
     1    ,I3,') = C_',I1,'(ISTATE=',I2,')'/53x,'=',1Pd14.7)
              ENDIF

c-----------------------------------------------------------------------
c** For DELR, calculate and write out the  A  and  B  coefficients
c-----------------------------------------------------------------------
          IF(PSEL(ISTATE).EQ.3) THEN
              yqRe=(RE(ISTATE)**nQB(ISTATE) - AREFQ**nQB(ISTATE))
     1                  /(RE(ISTATE)**nQB(ISTATE) + AREFQ**nQB(ISTATE))
              betaRe= beta(0,ISTATE)
              yPOW= 1.d0
              DO i= 1, Nbeta(ISTATE)
                  YPOW= YPOW*yqRe
                  betaRe= betaRe+ YPOW*beta(I,ISTATE)
                  ENDDO
              DO m=1,NCMM(ISTATE)
                  MMLR1D(m)= MMLR(m,ISTATE)
                  ENDDO
              CALL DAMPF(RE(ISTATE),rhoAB(ISTATE),NCMM(ISTATE),
     1                  MMLR1D,IVSR(ISTATE),IDSTT(ISTATE),DM,DMP,DMPP)
              IF(MMLR1D(1).LE.0) THEN
                  CALL AFdiag(RE(ISTATE),VLIM(ISTATE),NCMM(ISTATE),
     1             NCMMax,MMLR1D,CmVAL(1:NCMMax,ISTATE),rhoAB(ISTATE),
     2              IVSR(ISTATE),IDSTT(ISTATE),VATTRe,dULRdCm,dVATTRE)
                ELSE
                  DO m=1,NCMM(ISTATE)
                      Tm= CmVAL(m,ISTATE)/Re(ISTATE)**MMLR1D(m)
                      VATTRe= VATTRe + DM(m)*Tm
                      dVATTRe= dVATTRe + Tm*(DMP(m) 
     1                                  - Dm(m)*MMLR1D(m)/RE(ISTATE))
                      ENDDO
                ENDIF
              AA(ISTATE)= DE(ISTATE) - VATTRe - dVATTRE/betaRe
              BB(ISTATE)= AA(ISTATE) + DE(ISTATE) - VATTRe
              WRITE(6,633) 'A(DELR)',AA(ISTATE)
              WRITE(6,633) 'B(DELR)',BB(ISTATE)
              ENDIF
c-----------------------------------------------------------------------
c** Writing out the exponent expansion parameter information.
c-----------------------------------------------------------------------
          BTEMP= 0.d0
          IF(NPASS.GT.1) WRITE(20,671) NSR(ISTATE), Nbeta(ISTATE),
     1            nPB(ISTATE),nQB(ISTATE),RREFP(ISTATE),RREFQ(ISTATE)
          NALL= MAX(NSR(ISTATE),Nbeta(ISTATE)) 
          J=0 
          IF(NSR(ISTATE).LT.0) J=1 
          DO  I=J, NALL
              IPV= IPV + 1
              IF(IFXBETA(I,ISTATE).LE.0) THEN
                  IF(DABS(BETA(I,ISTATE)).GT.PU(IPV)) THEN
                      IF(NSR(ISTATE).GE.0) WRITE(6,640) 'be','ta',I,
     1                                   BETA(I,ISTATE),PU(IPV),PS(IPV)
                      IF(NSR(ISTATE).LT.0) WRITE(6,640) 'be','ta',I,
     1                  BETA(I,ISTATE),PU(IPV),PS(IPV),ypBETA(I,ISTATE)
                    ELSE
                      IF(NSR(ISTATE).GE.0) WRITE(6,641) 'be','ta',I,
     1                                   BETA(I,ISTATE),PU(IPV),PS(IPV)
                      IF(NSR(ISTATE).LT.0) WRITE(6,641) 'be','ta',I,
     1                  BETA(I,ISTATE),PU(IPV),PS(IPV),ypBETA(I,ISTATE)
                    ENDIF
                ELSE
                  IF(NSR(ISTATE).GE.0) WRITE(6,638) 'be','ta',I,
     1                                                   BETA(I,ISTATE)
                  IF(NSR(ISTATE).LT.0) WRITE(6,638) 'be','ta',I,
     1                                  BETA(I,ISTATE),ypBETA(I,ISTATE)
                ENDIF
              IF(NPASS.GT.1) THEN
                  IF(NSR(ISTATE).GE.0) THEN
                      WRITE(20,669) BETA(I,ISTATE),IFXBETA(I,ISTATE),
     1                                               'BETA',I,'BETA',I
                    ELSE
                      WRITE(20,668) ypBETA(I,ISTATE),BETA(I,ISTATE),
     1                                                IFXBETA(I,ISTATE)
                    ENDIF
                  ENDIF
              BTEMP= BTEMP+ BETA(I,ISTATE)
              ENDDO

c???       IF(PSEL(ISTATE).EQ.4) THEN   ! HPP ??why bother doing it here?
c              DO I= NALL+1, NALL+3
c                    IPV= IPV+1
c                    IF(IFXBETA(I,ISTATE).LE.0) THEN
c                        IF(DABS(BETA(I,ISTATE)).GT.PU(IPV)) THEN
c                           IF(NSR(ISTATE).GE.0) 
c     1        WRITE(6,640) 'pa','rm',I,BETA(I,ISTATE),PU(IPV),PS(IPV)
c                           IF(NSR(ISTATE).LT.0) 
c     1        WRITE(6,640) 'pa','rm',I,BETA(I,ISTATE),PU(IPV),PS(IPV),
c     2                                               ypBETA(I,ISTATE)
c                          ELSE
c                           IF(NSR(ISTATE).GE.0) 
c     1        WRITE(6,641) 'pa','rm',I,BETA(I,ISTATE),PU(IPV),PS(IPV)
c                           IF(NSR(ISTATE).LT.0) 
c     1        WRITE(6,641) 'pa','rm',I,BETA(I,ISTATE),PU(IPV),PS(IPV),
c     2                                               ypBETA(I,ISTATE)
c                          ENDIF
c                      ELSE
c                        IF(NSR(ISTATE).GE.0) 
c     1                        WRITE(6,638) 'pa','rm',I,BETA(I,ISTATE)
c                        IF(NSR(ISTATE).LT.0) 
c     1       WRITE(6,638) 'pa','rm',I,BETA(I,ISTATE),ypBETA(I,ISTATE)
c                      ENDIF
c                   ENDDO
c               ENDIF
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                  if(npass.gt.1) then 
cc         write(6,500) DE(ISTATE),ddedre,(j,dDe(j),j=0,nbeta(istate))
cc500 FORMAT(/'    De=',f12.7,'     d{De}/d{Re}=',1p,d12.4,'   and'/
cc   1    3('  d{De}/db(',i2,')=',d12.4):)
cc         write(6,510) CMMp2,dCmp2dRe, (j,dCmp2(j),j=0,nbeta(istate))
cc510 FORMAT(/'    CMMp2=',1P,d14.7,'   d{CMMp2}/d{Re}=',D12.4,'   and'/
cc   1    3('  d{Cmmp2}/db(',i2,')=',d12.4):)
cc                       endif
cccccccccccccccccccccccccccccccccccCCcccccccccccccccccccccccccccccccccccc
c
c** Write out  phi_\infty  constant for the EMO or DELR forms
          IF((PSEL(ISTATE).EQ.1).OR.(PSEL(ISTATE).EQ.3)) THEN 
              BINF= BTEMP
              WRITE(6,648) BINF
              IF(BINF.LT.0.d0) WRITE(6,647)
  647 FORMAT(' *** CAUTION *** negative beta_INf means potential blows u
     1p at large r ***')
              ENDIF
          IF(PSEL(ISTATE).EQ.2) THEN
              WRITE(6,648) BINF
              IF(MMLR(1,ISTATE).GT.0) THEN
                  BTEMP= CmVAL(1,ISTATE)*2.d0*(2.d0*BINF - BTEMP)
     1                                        *RE(ISTATE)**nPB(ISTATE)
                  WRITE(6,652) MMLR(1,ISTATE)+nPB(ISTATE),BTEMP
                ELSE
                  BTEMP= CmVAL(2,ISTATE)*2.d0*(2.d0*BINF - BTEMP)
     1                                        *RE(ISTATE)**nPB(ISTATE)
                  WRITE(6,652) MMLR(2,ISTATE)+nPB(ISTATE),BTEMP
                ENDIF
               ENDIF
c-----------------------------------------------------------------------
c** Writing out the adiabatic BOB radial function for atom A.
c-----------------------------------------------------------------------
          IF(NPASS.GT.1) THEN      !! next 4 - to stablize printout
               NUApr= NUA(ISTATE)-1
               NUBpr= NUB(ISTATE)-1
               NTApr= NTA(ISTATE)-1
               NTBpr= NTB(ISTATE)-1
               IF(NUA(ISTATE).LT.0)  NUApr= -1
               IF(NUB(ISTATE).LT.0)  NUBpr= -1
               IF(NTA(ISTATE).LT.0)  NTApr= -1
               IF(NTB(ISTATE).LT.0)  NTBpr= -1
               WRITE(20,672) NUApr, NUBpr,pAD(ISTATE),qAD(ISTATE),
     1                                                    LRad(ISTATE)
               ENDIF
          IF(NUA(ISTATE).GE.1) THEN
              DO  I= 0,NUA(ISTATE)-1
                  IPV= IPV + 1
                  IF(IFXUA(I,ISTATE).LE.0) THEN
                      IF(DABS(UA(I,ISTATE)).GT.PU(IPV)) THEN
                          WRITE(6,640) ' u',NAME(1),I,UA(I,ISTATE),
     1                                               PU(IPV),PS(IPV)
                        ELSE
                          WRITE(6,641) ' u',NAME(1),I,UA(I,ISTATE),
     1                                               PU(IPV),PS(IPV)
                        ENDIF
                    ELSE
                      WRITE(6,650) ' u',NAME(1),I,UA(I,ISTATE)
                    ENDIF
                  IF(NPASS.GT.1) WRITE(20,667) UA(I,ISTATE),
     1                                   IFXUA(I,ISTATE),'UA',I,'UA',I
                  ENDDO
              IPV= IPV + 1
              IF(NPASS.GT.1) WRITE(20,674) UA(NUA(ISTATE),ISTATE),
     1                       IFXUA(NUA(ISTATE),ISTATE),'uAinf','uAinf'
              IF(IFXUA(NUA(ISTATE),ISTATE).LE.0) THEN
                  WRITE(6,644) ' u',NAME(1),UA(NUA(ISTATE),ISTATE),
     1                                                 PU(IPV),PS(IPV)
                ELSE
                  WRITE(6,646) ' u',NAME(1),UA(NUA(ISTATE),ISTATE)
                ENDIF
              ENDIF    
  667 FORMAT(1Pd20.12,0P,I3,9x,'% ',A2,I2,'  IFX',A2,I2)
  668 FORMAT(1Pd20.12,d20.12,0P,I3)
  669 FORMAT(1Pd20.12,0P,I3,9x,'% ',A4,I2,'  IFX',A4,I2)
  670 FORMAT(1Pd20.12,0P,I3,9x,'% ',A2,'  IFX',A2)
  671 FORMAT(/2I3,I4,I3,1PD11.2,8x,'% NSR Nbeta nPB nQB RREFP RREFQ')
  672 FORMAT(/2I3,I4,I3,I5,14x,'% NUA NUB pAD qAD LRad')
  673 FORMAT(/2I3,I4,I3,19x,'% NTA NTB pNA qNA')
  674 FORMAT(1Pd20.12,0P,I3,9x,'% ',a5,'  IFX',A5)
  675 FORMAT(/3I3,24x,'% NwCFT Pqw efREF')
c-----------------------------------------------------------------------
c** Writing out the adiabatic BOB radial function for atom B.
c-----------------------------------------------------------------------
          IF(NUB(ISTATE).GE.1) THEN
              DO  I= 0, NUB(ISTATE)- 1
                  IPV= IPV + 1
                  IF(IFXUB(I,ISTATE).EQ.0) THEN
                      IF(DABS(UB(I,ISTATE)).GT.PU(IPV)) THEN
                          WRITE(6,640) ' u',NAME(2),I,UB(I,ISTATE),
     1                                                 PU(IPV),PS(IPV)
                        ELSE
                          WRITE(6,641) ' u',NAME(2),I,UB(I,ISTATE),
     1                                             PU(IPV),PS(IPV)
                        ENDIF
                    ELSE
                      WRITE(6,650) ' u',NAME(2),I,UB(I,ISTATE)
                    ENDIF
                  IF(NPASS.GT.1) WRITE(20,667) UB(I,ISTATE),
     1                                   IFXUB(I,ISTATE),'UB',I,'UB',I
                  END DO
              IPV= IPV + 1
              IF(NPASS.GT.1) WRITE(20,674) UB(NUB(ISTATE),ISTATE),
     1                       IFXUB(NUB(ISTATE),ISTATE),'uBinf','UBinf'
              IF(IFXUB(NUB(ISTATE),ISTATE).LE.0) THEN
                  WRITE(6,644) ' u',NAME(2),UB(NUB(ISTATE),ISTATE),
     1                                             PU(IPV),PS(IPV)
                ELSE
                  WRITE(6,646) ' u',NAME(2),UB(NUB(ISTATE),ISTATE)
                ENDIF
              ENDIF
c-----------------------------------------------------------------------
c** Writing out the Rotational Non-Adiabatic information for atom A.
c-----------------------------------------------------------------------
          IF(NPASS.GT.1) WRITE(20,673) NTApr, NTBpr,pNA(ISTATE),
     1                                                     qNA(ISTATE)
          IF(NTA(ISTATE).GE.1) THEN
              DO  I= 0, NTA(ISTATE)-1
                  IPV= IPV + 1
                  IF(IFXTA(I,ISTATE).LE.0) THEN
                      IF(DABS(TA(I,ISTATE)).GT.PU(IPV)) THEN
                          WRITE(6,640) ' t',NAME(1),I,TA(I,ISTATE),
     1                                             PU(IPV),PS(IPV)
                        ELSE
                          WRITE(6,641) ' t',NAME(1),I,TA(I,ISTATE),
     1                                             PU(IPV),PS(IPV)
                        ENDIF
                    ELSE
                      WRITE(6,650) ' t',NAME(1),I,TA(I,ISTATE)
                    ENDIF
                  IF(NPASS.GT.1) WRITE(20,667) TA(I,ISTATE),
     1                                   IFXTA(I,ISTATE),'TA',I,'TA',I
                  END DO
              IPV= IPV + 1
              IF(NPASS.GT.1) WRITE(20,674) TA(NTA(ISTATE),ISTATE),
     1                       IFXTA(NTA(ISTATE),ISTATE),'tAinf','TAinf'
              IF(IFXTA(NTA(ISTATE),ISTATE).LE.0) THEN
                  WRITE(6,644) ' t',NAME(1),TA(NTA(ISTATE),ISTATE),
     1                                             PU(IPV),PS(IPV)
                ELSE
                  WRITE(6,646) ' t',NAME(1),TA(NTA(ISTATE),ISTATE)
                  ENDIF
              ENDIF
c-----------------------------------------------------------------------
c** Writing out the Rotational Non-Adiabatic information for atom B.
c-----------------------------------------------------------------------
          IF(NTB(ISTATE).GE.1) THEN
              DO  I= 0, NTB(ISTATE)-1
                  IPV= IPV + 1
                  IF(IFXTB(I,ISTATE).LE.0) THEN
                      IF(DABS(TB(I,ISTATE)).GT.PU(IPV)) THEN
                          WRITE(6,640) ' t',NAME(2),I,TB(I,ISTATE),
     1                                             PU(IPV),PS(IPV)
                        ELSE
                          WRITE(6,641) ' t',NAME(2),I,TB(I,ISTATE),
     1                                             PU(IPV),PS(IPV)
                        ENDIF
                    ELSE
                      WRITE(6,650) ' t',NAME(2),I,TB(I,ISTATE)
                    ENDIF
                  IF(NPASS.GT.1) WRITE(20,667) TB(I,ISTATE),
     1                                   IFXTB(I,ISTATE),'TB',I,'TB',I
                  END DO
              IPV= IPV + 1
              IF(NPASS.GT.1) WRITE(20,674) TB(NTB(ISTATE),ISTATE),
     1                       IFXTB(NTB(ISTATE),ISTATE),'tBinf','TBinf'
              IF(IFXTB(NTB(ISTATE),ISTATE).LE.0) THEN
                  WRITE(6,644) ' t',NAME(2),TB(NTB(ISTATE),ISTATE),
     1                                             PU(IPV),PS(IPV)
                ELSE
                  WRITE(6,646) ' t',NAME(2),TB(NTB(ISTATE),ISTATE)
                ENDIF
              ENDIF
c-----------------------------------------------------------------------
c** Writing out Lambda-doubling/2-Sigma coefficients
c-----------------------------------------------------------------------
          IF((IOMEG(ISTATE).GT.0).OR.(IOMEG(ISTATE).EQ.-1)) THEN
              IF(NPASS.GT.1) WRITE(20,675) NwCFT(ISTATE),Pqw(ISTATE),
     1                                                   efREF(ISTATE)
              ENDIF
          IF(NwCFT(ISTATE).GE.0) THEN
              J= 1
              IF(IOMEG(ISTATE).EQ.-1) J=2
              DO  I= 0, NwCFT(ISTATE)
                  IPV= IPV + 1
                  IF(IFXwCFT(I,ISTATE).LE.0) THEN
                      IF(DABS(wCFT(I,ISTATE)).GT.PU(IPV)) THEN
                          WRITE(6,642) NAMEDBLE(J),I,wCFT(I,ISTATE),
     1                                                 PU(IPV),PS(IPV)
                        ELSE
                          WRITE(6,643) NAMEDBLE(J),I,wCFT(I,ISTATE),
     1                                                 PU(IPV),PS(IPV)
                        ENDIF
                    ELSE         
                      WRITE(6,651) NAMEDBLE(J),I,wCFT(I,ISTATE)
                    ENDIF
                  IF(NPASS.GT.1) WRITE(20,669) wCFT(I,ISTATE),
     1                            IFXwCFT(I,ISTATE),'wCFT',I,'wCFT,',I
                  END DO
              ENDIF
   90     CONTINUE
      WRITE(6,600)
      RETURN
c-----------------------------------------------------------------------
  552 FORMAT(11x,'using radial expansion variable:   y',I1,' = (R^',I1,
     1   ' - Re^',I1,')/(R^',I1,' + Re^',I1,')')
  555 FORMAT(8x,'with radial variable:   y_{p,q} = (R^q -',F9.6,'^q)/(R^
     1q +',F9.6,'^q)')
  554 FORMAT(8x,'with ',A2,'-atom radial expansion variable:   y',I1,
     1  ' = (R^',I1,' - Re^',I1,')/(R^',I1,' + Re^',I1,')')
  556 FORMAT(' Adiabatic BOB functions are simple power series in  y_'
     1  I1,'(r) scaled by  m_e/M(A):') 
  557 FORMAT(' Adiabatic BOB functions with {p=',i2,', q=',i2,
     1  '} are scaled by  DELTA{M(A)}/M(A):')
  564 FORMAT(5x,A10,A2,';R) = \sum\{u_i * [y',I1,']^i}   for  i= 0 to',
     1 i3)
  558 FORMAT(5x,A10,A2,';R) = u(inf)*y',i1,' + (1 - y',I1, ')*\Sum{u_i *
     1 [y',I1,']^i}   for  i= 0 to',i3)
  559 FORMAT(' Non-Adiabatic centrifugal BOB fx. with {p=',i2,', q=',i2,
     1  '} are scaled by  M(1)/M(A):')
  560 FORMAT(' Non-Adiabatic centrifugal BOB fx are power series in y_',
     1  I2,'(r) scaled by  m_e/M(A):')
  600 FORMAT(1X,39('=='))
  601 FORMAT(' All distinct levels of State ',A3,' fitted as independent
     1 term values')
  603 FORMAT(/' FIXED State ',A3,' potential defined by interpolating ov
     1er input turning points')
  684 FORMAT(/' For state ',A3,' fits represents level with Band Constan
     1ts'/1x,6('=='))
  686 FORMAT(2X,A6,I3,'; IS=',I2,')=',1PD20.12,3X,1PD8.1,6X,1PD8.1,10x)
  687 FORMAT(I4,10I4)
 6872 FORMAT(4x,10I4)
  688 FORMAT('  ')
  690 FORMAT(7x,'Parameter',10x,'Final Value',5x'Uncertainty   Sensitivi
     1ty')
  602 FORMAT(/' State ',A3,' represented by an EMO(N=',I2,' q=',i2,') po
     1tential defined in terms of'/1x,4('=='),  '  exponent coefficient:
     2   beta(R)= Sum{beta_i*y',i1,'^i}  for  i= 0 -',i3)
  604 FORMAT(/' State ',A3,' represented by an MLR(p=',i2,', q=',i2,
     1  ') potential defined in terms of')
  605 FORMAT(1x,4('=='), '  exponent coefficient:  beta(R)= betaINF*y',
     1 i1,' +(1-y',I1,')*Sum{beta_i*y',i1,'^i}'/62x,'for  i= 0 to',I3)
  610 FORMAT(/' For state ',A3,"  use Surkus' Generalized Potential Ener
     1gy Function GPEF with"/1x,6('=='),'  expansion vble:   y_',i1, 
     2 '(r) = (r^',i1,' - re^',i1,')/(',F5.2,'*r^',i1,' +',F5.2,'*re^',
     3 i1,'p)')
 6102 FORMAT(' *** Input ERROR *** band constant specification  v=',I3,
     1  ' .NE.', I3)
  612 FORMAT(/' State ',A3,' represented by a DELR(N=',I2,'  q=',i2,') 
     1potential with'/1x,4('=='),'   exponent coefficient:   beta(r)= Su
     2m{beta_i*y',i1,'^i}'/8x,'with polynomial order   N=',I3,' for  r <
     3 r_e  and   N=',I3,'  for  r.ge.r_e')
  607 FORMAT(4x,'uLR inverse-power terms incorporate DS-type damping wit
     1h   rhoAB=',f10.7/11x,'defined to give very short-range  Dm(r)*Cm/
     2r^m  behaviour   r^{',SS,f4.1,'}'/8x,'Dm(r)= [1 - exp(-',f5.2,
     3 '(rhoAB*r)/m -',f6.3,'(rhoAB*r)^2/sqrt{m})]^{m',SP,F4.1,'}')
 6072 FORMAT(4x,18('-'),11(A5: ))
  617 FORMAT(4x,'uLR inverse-power terms incorporate TT-type damping wit
     1h   rhoAB=',f10.7/8x,'defined to give very short-range  Dm(r)*Cm/r
     2^m  behaviour   r^{',I2,'}'/8x,'Dm(r)= [1 - exp(-bTT*r)*SUM{(bTT*r
     1)^k/k!}]   where   bTT=',f6.3,'*rhoAB') 
  680 FORMAT(1x,4('=='),'  exponent coefficient defined as a Pashov natu
     1ral spline through'/ 10x,i3,' specified points including the fixed
     2 value of  betaINF  at  yp= 1')
  682 FORMAT(20x,'These constants yield:   betaINF=',F14.10 )
  683 FORMAT(/4x,'Since this state has (projected) electronic angular mo
     1mentum  OMEGA=',I2/10x,'eigenvalue calculations use centrifugal po
     2tential  [J*(J+1) -',I2,']/r**2')
  608 FORMAT(48x,'C',I1,'=',1PD14.7,'[cm-1 Ang^',i1,']')
 6082 Format(3x,I4,9x,12I5)
  609 FORMAT(47x,'C',I2,'=',1PD14.7,'[cm-1 Ang^{',i2,'}]')
 6092 Format('    n(q{lambda})  ',12I5)
  606 FORMAT(' Use Aubert-Frecon 2x2 ',A2,'-state uLR(r)  with   Aso=',
     1 F11.6,'[cm-1]')
 6062 FORMAT(/' For state ',A3,' represent level energies by independent
     1 band constant for each'/1x,6('=='),27x,'vibrational levels of eac
     2h isotologue'/3x,'No. band constants'/6x,'v','   isotop:  #',I1:
     3   11('   #',I1:):)
 6066 FORMAT(' Use Aubert-Frecon 3x3 ',A2,'-state uLR(r)  with   Aso=',
     1 F11.6,'[cm-1]')
  614 FORMAT(/'   Parameter    Initial Value    Uncertainty   Sensitivit
     1y')
  615 FORMAT(/'   Parameter      Final Value    Uncertainty   Sensitivit
     1y')
  618 FORMAT(1x,A6,'-doubling splitting strength function is expanded as
     1'/7x,'f(r) =  Sum{w_i * (y',i1,')^i}   for   i= 0 to',i3/
     2 11x,'with radial expansion variable:   y',I1,' = (R^',I1,
     3  ' - Re^',I1,')/(R^',I1,' + Re^',I1,')')
  619 FORMAT(/' Including BOB term makes centrifugal potential strength 
     1factor   [J(J+1) +',I2,']')
  620 FORMAT(6X,A2,3X,1PD21.12,3X,1PD8.1,6X,1PD8.1,10x)
  621 FORMAT(5X,'*',A2,3X,1PD21.12,3X,1PD8.1,6X,1PD8.1)
  622 FORMAT(6X,A2,3X,1PD21.12,7X,'--',12X,'--')
  623 FORMAT(/' State ',A3,' represented by a Tiemann polynomial (b=',
     1   F5.2,', R_in=',F4.2,', R_out=',F4.2,')'/1x,4('=='),5x,'with exp
     2ansion variable:  y_',i1,'(r) = (r - re)/(r ',SP,F5.2,'*re)') 
c    1,3726(1984)]'/5x,'[1 - exp(-3.13*RHO*r)*SUM{(3.13*RHO*r)^k/k!}]',
c    2 '   with  rhoAB=',F9.6)
  630 FORMAT(3X,A5,3X,1PD21.12,3X,1PD8.1,6X,1PD8.1) 
  631 FORMAT(2X,'*',A5,3X,1PD21.12,3X,1PD8.1,6X,1PD8.1) 
  632 FORMAT(3X,A5,3X,1PD21.12,7X,'--',12X,'--')
  633 FORMAT(4x,A7,1PD21.12)
  636 FORMAT(4X,A4,3X,1PD21.12,7X,'--',12X,'--')
  638 FORMAT(3X,A2,A2,'(',I2,')',1PD21.12,7X,'--',12X,'--':
     1  '     at   yp=',0pF14.10)
  640 FORMAT(3X,A2,A2,'(',I2,')',1PD21.12,3X,1PD8.1,6X,1PD8.1:
     1  '   at   yp=',0pF14.10)
  641 FORMAT('  *',2A2,'(',I2,')',1PD21.12,3X,1PD8.1,6X,1PD8.1:
     1  '   at  yp=',0pF14.10)
  642 FORMAT(1X,A6,'(',I2,')',1PD21.12,3X,1PD8.1,6X,1PD8.1)
  643 FORMAT('*',A6,'(',I2,')',1PD21.12,3X,1PD8.1,6X,1PD8.1)
  644 FORMAT(1x,A2,'_inf(',A2,')',1PD21.12,1PD11.1,6X,1PD8.1)
  646 FORMAT(1x,A2,'_inf(',A2,')',1PD21.12,7X,'--',12X,'--')
  648 FORMAT(3X,'beta_INF',1PD21.12,7X,'--',12X,'--') 
cc649 FORMAT(3X,'beta_INF',1PD21.12,7X,'--',12X,'--')
cc   ,5x,'at   yp=  1.000 1000000')
  650 FORMAT(3X,A2,A2,'(',I2,')',1PD21.12,7X,'--',12X,'--')
  651 FORMAT(1X,A6,'(',I2,')',1PD21.12,7X,'--',12X,'--')
  652 FORMAT('   C',I2,'{exp}',1PD21.12,7X,'--',12X,'--')
  654 FORMAT(1X,A6,'_inf',1PD21.12,3X,1PD8.1,6X,1PD8.1)
  656 FORMAT('*',A6,'_inf',1PD21.12,3X,1PD8.1,6X,1PD8.1)
  658 FORMAT(1X,A6,'_inf',1PD21.12,7X,'--',12X,'--') 
  660 FORMAT(5X,'C',I2,3X,1PD21.12,3X,1PD8.1,6X,1PD8.1) 
  661 FORMAT(3X,'* C',I2,3X,1PD21.12,3X,1PD8.1,6X,1PD8.1) 
  662 FORMAT(5X,'C',I2,3X,1PD21.12,7X,'--',12X,'--':)
  664 FORMAT(4x,'uLR inverse-power terms incorporate NO damping function
     1s')
  692 FORMAT(' ', A3,' state  Lambda-doubling split levels referenced to
     1 f-parity levels')
  694 FORMAT(' ', A3,' state  Lambda-doubling split levels referenced to
     1 the e/f level mid-point')
  696 FORMAT(' ', A3,' state  Lambda-doubling split levels referenced to
     1 e-parity levels')
  700 FORMAT(1x,"'",A3,"'",2I3,2I5,I3,6x,' % ('I1')SLABL IOMEG VMIN' 
     1         ' VMAX JTRUNC EFSEL')
  701 FORMAT(I3,F12.4,3I3,7x,' % PSEL VLIM MAXMIN BOBCN OSEL')
  702 FORMAT(2F7.2,F8.4,9x,' % RMIN RMAX RH')
  703 FORMAT(I3,F7.2,2I3,15x,' % NCMM RHOab IDF IDSTT')
  704 FORMAT(I3,1P,D17.8,0P,I3,8x' % MMLR(',I1,'), CmVAL(',I1,
     1  '), IFXCm(',I1,')') 
  705 FORMAT(10I4)
  706 FORMAT(I3,1P,D17.8,0P,I3,8x' % ',A6)
  708 FORMAT(44x,A6,'=',1PD14.7,'[cm-1 Ang^{',i2,'}]')
  709 FORMAT(44x,A6,'=',1PD14.7,'[cm-1]')
  720 FORMAT(2X,A6,3X,1PD21.12,3X,1PD8.1,6X,1PD8.1,10x)
  721 FORMAT(1X,'*',A6,3X,1PD21.12,3X,1PD8.1,6X,1PD8.1)
  722 FORMAT(2X,A6,3X,1PD21.12,7X,'--',12X,'--')
 
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
