c***********************************************************************
      SUBROUTINE WRITEPOT(NPASS,SLABL,NAME,DECM,PV,PU,PS,CM,VMAXIN) 
c***********************************************************************
c** Subroutine to print out complete description of the potential fx.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++       Version of  16 May 2016  {after generatized TT ipgrade}
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
      CHARACTER*3  SLABL(-6:NSTATEMX)
      CHARACTER*5 DASH
      CHARACTER*6  BCNAM(8),QCNAM(8)
      CHARACTER*7  NAMEDBLE(2)
      CHARACTER*10 LAB2(7),LAB3(10)
      INTEGER NPASS, ISTATE,IPV,I,I1,I4,ISOT,J,MMN,m,m1,LSR,MMp2,
     1  JSTATE,NUApr,NUBpr,NTApr,NTBpr,MCMM,IPVRe(NSTATEMX),
     2  MMLR1D(NCMMAX),VMAXIN(NSTATEMX)
      REAL*8 DECM(NSTATEMX),BTEMP,UAT,UBT,SAT,BINF,RE3,RE6,RE8,T0,T1,
     1 DX,DX1,ULRe,C3VAL,C6adj,C9adj,RET,RETSig,RETPi,RETp,RETm,Tm,
     2 VATTRe,dVATTRe,PVSR,ATT,YP,YPP,BT,Rinn,Rout,A1,A2,A3,B5,VX,dVX,
     3 uLR,CMMp2,uCMMp2,uDe,XRO,dXRO,dXROdRe,d2XROdRe,fRO,XROpw,ROmp2,
     4 dCmp2dRe,dDeROdRe,yqRe,betaRe,yPOW,XRI,AREF,ttVMIN,ttRMIN,RR,
     5 bohr,f2,f2p,rhoINT,
     5 dCmp2(0:NbetaMX),dULRdCm(NCMMax),DM(MMAX),DMP(MMAX),DMPP(MMAX),
     6 bTT(-2:2),cDS(-4:4),bDS(-4:4),CmVAL1D(NCMMAX),CmEFF1D(NCMMAX)
      DATA QCNAM/' QB(v=','QD(v=',' QH(v=',' QL(v=',' QM(v=',
     1    ' QN(v=',' QO(v=',' '/
      DATA BCNAM/' Tv(v=',' Bv(v=','-Dv(v=',' Hv(v=',' Lv(v=',' Mv(v=',
     1    ' Nv(v=',' Ov(v='/
c** Damping function factors from Table 1 of Mol.Phys. 109, 435 (2011)
       DATA bTT/2.10d0,2.44d0,2.78d0,3.13d0,3.47d0/
       DATA bDS/2.50d0,2.90d0,3.3d0,3.69d0,3.95d0,0.d0,4.53d0,0.d0,
     1          4.99d0/
       DATA cDS/0.468d0,0.446d0,0.423d0,0.405d0,0.390d0,0.d0,0.360d0,
     1            0.d0,0.340d0/
c...For testing: precise Scolegian values of 'b' and 'c' for s=0 ......
cc    DATA bDS/2.50d0,2.90d0,3.30d0,3.69d0,3.968424883d0,4*0.d0/
cc    DATA cDS/0.468d0,0.446d0,0.423d0,0.405d0,0.3892460703d0,4*0.d0/

      DATA DASH/'-----'/ 
      SAVE bTT, bDS, cDS
c c** NLLSSRR variables used in output c
      REAL*8 PV(NPARMX), PU(NPARMX), PS(NPARMX),CM(NPARMX,NPARMX),
     1                                                      PT(NPARMX)
      DATA NAMEDBLE/'wLambda',' wSigma'/
c** labels for matrix cases of MLR: LAB2 for 2x2, LAB3 for 3x3, LAB4 for names
      DATA LAB2/'  DELTAE  ',' C3(^1Sig)',' C3(^3Pi)' ,' C6(^1Sig)',
     1                       ' C6(^3Pi) ',' C8(^1Sig)',' C8(^3Pi) '/
      DATA LAB3/'DELTAE  ',' C3(^3Sig)',' C3(^1Pi) ',' C3(^3Pi) ',
     1                     ' C6(^3Sig)',' C6(^1Pi) ',' C6(^3Pi) ',
     2                     ' C8(^3Sig)',' C8(^1Pi) ',' C8(^3Pi) '/
      DATA LAB4/' ?',' B',' c',' b',' A'/
      DATA bohr/0.52917721092d0/     !! 2010 physical constants d:mohr12
c-----------------------------------------------------------------------
c** Writing out state specific information.
c-----------------------------------------------------------------------
      IPV= 0 
      DO  90 ISTATE=1,NSTATES
          VATTRe= 0.d0
          dVATTRe= 0.d0
          WRITE(6,600)
          IF(PSEL(ISTATE).EQ.0) THEN
              WRITE(6,603) SLABL(ISTATE)
              WRITE(6,636) 'VLIM',VLIM(ISTATE)
              GOTO 50
              ENDIF
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
                                  IF(DABS(PU(IPV)).GT.DABS(PV(IPV)))THEN
                                      WRITE(6,685) BCNAM(J),I,ISOT,
     1                                         PV(IPV),PU(IPV),PS(IPV)
                                    ELSE
                                      WRITE(6,686) BCNAM(J),I,ISOT,
     1                                         PV(IPV),PU(IPV),PS(IPV)
                                    ENDIF
                                  ENDDO
                              IF(NQC(I,ISOT,ISTATE).GT.0) THEN
                                  DO  J= 1,NQC(I,ISOT,ISTATE)
                                      IPV= IPV+1   !! Lambda Count
                                      IF(DABS(PU(IPV)).GT.DABS(PV(IPV)))
     1                                                            THEN
                                          WRITE(6,685) QCNAM(J),I,ISOT,
     1                                         PV(IPV),PU(IPV),PS(IPV)
                                        ELSE  
                                          WRITE(6,686) QCNAM(J),I,ISOT,
     1                                         PV(IPV),PU(IPV),PS(IPV)
                                        ENDIF  
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
          AREF= RREFq(ISTATE)
          IF(AREF.LE.0.d0) AREF= RE(ISTATE)
          IF(PSEL(ISTATE).EQ.1) THEN
c** Header printout for EMO potential
              WRITE(6,602) SLABL(ISTATE),Nbeta(ISTATE),qPOT(ISTATE),
     1                                       qPOT(ISTATE),Nbeta(ISTATE)
              IF(RREFq(ISTATE).LE.0.d0) WRITE(6,552)(qPOT(ISTATE),i=1,4)
              IF(RREFq(ISTATE).GT.0.d0) WRITE(6,555) (qPOT(ISTATE),AREF,
     1  qPOT(ISTATE),i=1,2)
              ENDIF
          IF(PSEL(ISTATE).EQ.2) THEN
c** Header printout for MLR potential
              BINF= betaINF(ISTATE)
              WRITE(6,604) SLABL(ISTATE),qPOT(ISTATE),pPOT(ISTATE) 
              IF(APSE(ISTATE).LE.0) THEN
                  WRITE(6,605) pPOT(ISTATE),pPOT(ISTATE),qPOT(ISTATE),
     1                                                   Nbeta(ISTATE)
                ELSE
                  WRITE(6,680) Nbeta(ISTATE)
                  BETA(Nbeta(ISTATE),ISTATE)= BINF
                ENDIF
              IF(RREFq(ISTATE).LE.0.d0) WRITE(6,552)(qPOT(ISTATE),i=1,4)
              IF(RREFq(ISTATE).GT.0.d0) WRITE(6,555) (qPOT(ISTATE),AREF,
     1                                             qPOT(ISTATE),i=1,2)
              IF(RREFp(ISTATE).LE.0.d0) WRITE(6,553) (pPOT(ISTATE),
     1                               RREFq(ISTATE),pPOT(ISTATE),i=1,2)
              IF(RREFp(ISTATE).GT.0.d0) WRITE(6,551) (pPOT(ISTATE),
     1                               RREFp(ISTATE),pPOT(ISTATE),i=1,2)
              ENDIF
c
          IF(PSEL(ISTATE).EQ.3) THEN
c** Header printout for DELR potential form ...
              WRITE(6,612) SLABL(ISTATE),Nbeta(ISTATE),qPOT(ISTATE),
     1                                       qPOT(ISTATE),Nbeta(ISTATE)
              IF(RREFq(ISTATE).LE.0.d0) WRITE(6,552)(qPOT(ISTATE),i=1,4)
              IF(RREFq(ISTATE).GT.0.d0) WRITE(6,555) (qPOT(ISTATE),AREF,
     1  qPOT(ISTATE),i=1,2)
              ENDIF
c
          IF(PSEL(ISTATE).EQ.4) THEN
c** Header printout for Surkus GPEF potential form ...
              WRITE(6,610) SLABL(ISTATE),(qPOT(ISTATE),i=1,3),
     1             AGPEF(ISTATE),qPOT(ISTATE),BGPEF(ISTATE),qPOT(ISTATE)
              ENDIF
c
          IF(PSEL(ISTATE).EQ.5) THEN
c** Header printout for Tiemann HPP potential ...
c... First, need to define long-and short-range connections .... 
c ... Begin by getting V(r) and V'(r) of polynomial VX at R_i and R_o
              WRITE(6,623) SLABL(ISTATE), BETA(Nbeta(ISTATE)+1,ISTATE),
     1        BETA(Nbeta(ISTATE)+2,ISTATE),BETA(Nbeta(ISTATE)+3,ISTATE),
     2                      (pPOT(ISTATE)),BETA(Nbeta(ISTATE)+1, ISTATE)
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
              ENDIF
c  ..................................  end of HPP printout .............
          IF(PSEL(ISTATE).EQ.6) THEN
c** Header printout for Generalized Tang-Toennies type potential ...  
c  first locate ACTUAL minimum to compare with input D_e and r_e values 
              I1= (0.9*RE(ISTATE)- RMIN(ISTATE))/RH(ISTATE) 
              RR= RD(I1,ISTATE)
              DO m=1,NCMM(ISTATE)
                  MMLR1D(m)= MMLR(m,ISTATE)
                  ENDDO
              ttVMIN= 0.d0
              ttRMIN= RR
              rhoINT= rhoAB(ISTATE)/bTT(IVSR(ISTATE)/2)
              A1= 0.d0
              A2= 0.d0    
   55         CALL dampF(RR,rhoINT,NCMM(ISTATE),NCMMAX,MMLR1D,
     1                        IVSR(ISTATE),IDSTT(ISTATE),DM,DMP,DMPP)
              A1= A2
              A2= A3
c....calculate the (damped) long range tail
              T0= 0.d0
              DO m= 1, NCMM(ISTATE)
                  T0= T0+ DM(m)*CmVAL(m,ISTATE)/RR**MMLR(m,ISTATE)
                  ENDDO
c.... Now evaluate Generalized TT model
              T1= BETA(1,ISTATE)*RR+ BETA(2,ISTATE)*RR**2+
     1                        BETA(3,ISTATE)/RR + BETA(4,ISTATE)/RR**2
              A3= (BETA(5,ISTATE)+ BETA(6,ISTATE)*RR+ BETA(7,ISTATE)/RR 
     1     + BETA(8,ISTATE)*RR**2+ BETA(9,ISTATE)*RR**3)*DEXP(-T1)- T0
              IF(A3.LE.ttVMIN) THEN
c... search for potential minimum ...
                  ttVMIN= A3
                  ttRMIN= RR
                  RR= RR+ RH(ISTATE)
                  GOTO 55
                  ENDIF
              WRITE(6,626) (BETA(i,ISTATE),i=1,9)
c*** Use quadratic approximation to determine REQ and DSCM
              T0= (A3- 2.d0*A2 + A1)/(2.d0*RH(ISTATE)**2)   !! curvature
              RR= ttRMIN
              ttRMIN= ttRMIN+ 0.5d0*RH(ISTATE)
     1                                   -(A3-A2)/(2.d0*RH(ISTATE)*T0)
              ttVMIN=  T0*(RR- ttRMIN)**2 - A2
              WRITE(6,627) DE(ISTATE), RE(ISTATE), ttVMIN, ttRMIN
              ENDIF
  626 FORMAT(/' Use a Generalized Tang-Tonnies Potential function with e
     1xponent'/' - {{',SP,F15.11,'*r',F15.11,'*r^2',F15.11,'/r',F15.11,
     2  '/r^2}}'/' and pre-exponential factor:'/3x,'{{',SP,1PD15.8,
     3  D16.8,'*r',d16.8,'/r',d16.8,'*r^2'/21x,D16.8,'*r^3}}',S)
  627 FORMAT(10x,'Input    DSCM=',F10.4,'   REQ=',f9.6/ 10x,
     1 'Actual   DSCM=',F10.4,'   REQ=',f9.6)
c======================================================================
          IF(PSEL(ISTATE).EQ.7) THEN
              IF(Nbeta(ISTATE).EQ.5) THEN
c** For Aziz'ian HFD-ABC type potential: print header and derive leading 
c exponent coefficient \beta_1  and pre-exponential factor A for use 
c   in subroutines 'vgen' & 'vgenp'; all in units cm-1 and \AA
                  A1= BETA(1,ISTATE) 
                  A2= BETA(2,ISTATE) 
                  A3= BETA(3,ISTATE)
                  DX= 1.d0
                  DX1= 0.d0
                  IF(A2.GT.RE(ISTATE)) THEN
                      DX= DEXP(-A1*(A2/RE(ISTATE) - 1.d0)**A3) 
                      DX1= A1*A2*A3*DX*(A2/RE(ISTATE)- 1.d0)**(A3- 1.d0)
     1                                                  /RE(ISTATE)**2
                      ENDIF
                  T0= 0.d0
                  T1= 0.d0
                  DO m= 1,NCMM(ISTATE)
                      Tm= CMVAL(m,ISTATE)/RE(ISTATE)**MMLR(m,ISTATE)
                      T0= T0+ Tm
                      T1= T1+ Tm*(DX1 - MMLR(m,ISTATE)*DX/RE(ISTATE))
                      ENDDO
                  T0= T0*DX - DE(ISTATE)
                  IF(T0.LE.0.d0) THEN
                       WRITE(6,624) T0,(MMLR(m,ISTATE),CmVAL(m,ISTATE),
     1                                             m= 1, NCMM(ISTATE))
                       STOP
                       ENDIF
                  BB(ISTATE) = BETA(5,ISTATE)/RE(ISTATE) 
     1                        - 2.d0*BETA(4,ISTATE)*RE(ISTATE) - T1/T0
                  AA(ISTATE)= T0 * DEXP(RE(ISTATE)
     1                      *(BB(ISTATE) + BETA(4,ISTATE)*RE(ISTATE)))
                  WRITE(6,628) 'ABC',BETA(5,ISTATE),BB(ISTATE),
     1                                      BETA(4,ISTATE),AA(ISTATE),
     2             (MMLR(m,ISTATE),CmVAL(m,ISTATE),m= 1, NCMM(ISTATE))
                  WRITE(6,629) DE(ISTATE),RE(ISTATE),A2,A1,A2,A3
                ELSEIF(Nbeta(ISTATE).EQ.2) THEN
c------------------------------------------------------------------------
c** For Aziz'ian HFD-ID type potential: print header and derive leading 
c exponent coefficient \beta_1  and pre-exponential factor A for use 
c   in subroutines 'vgen' & 'vgenp'; all in units cm-1 and \AA
                  f2= 1.d0 - (rhoAB(ISTATE)*RE(ISTATE)/bohr)**1.68d0 
     1                   *EXP(-0.78d0*rhoAB(ISTATE)*(RE(ISTATE)/bohr))
                  f2p= (f2- 1.d0)*(1.68d0/RE(ISTATE) -
     1                                      0.78d0*rhoAB(ISTATE)/bohr)
                  DO m=1,NCMM(ISTATE)
                      MMLR1D(m)= MMLR(m,ISTATE)
                      ENDDO
                  CALL dampF(RE(ISTATE),rhoAB(ISTATE),NCMM(ISTATE),
     1          NCMMAX,MMLR1D,IVSR(ISTATE),IDSTT(ISTATE),DM,DMP,DMPP)
                  T0= 0.d0
                  T1= 0.d0
                  DO m= 1,NCMM(ISTATE)
                      Tm= CMVAL(m,ISTATE)/RE(ISTATE)**MMLR(m,ISTATE)
                      T0= T0+ Dm(m)*Tm
                      T1= T1+ Tm*(f2p*Dm(m) + f2*(Dmp(m) 
     1                                  - Dm(m)*MMLR1D(m)/RE(ISTATE)))
                      ENDDO
                  T0= T0*f2
                  BB(ISTATE) = BETA(2,ISTATE)/RE(ISTATE)
     1        - 2.d0*BETA(1,ISTATE)*RE(ISTATE) - T1/(T0 - DE(ISTATE))
                  AA(ISTATE)= (T0 - DE(ISTATE))*EXP(RE(ISTATE)
     1                        *(BB(ISTATE)+ BETA(1,ISTATE)*RE(ISTATE)))
                  WRITE(6,628) 'ID ',BETA(2,ISTATE),BB(ISTATE),
     1                                      BETA(1,ISTATE),AA(ISTATE),
     2             (MMLR(m,ISTATE),CmVAL(m,ISTATE),m= 1, NCMM(ISTATE))
                  WRITE(6,629) DE(ISTATE),RE(ISTATE)
                ELSEIF((Nbeta(ISTATE).NE.2).AND.(Nbeta(ISTATE).NE.5))
     1                                                            THEN
                  Write (6,625) Nbeta(ISTATE)
                  STOP
                ENDIF
              ENDIF
  624 FORMAT(/' *** ERROR in generating HFD potential *** generate   VAT
     1T=',G15.7,'  from Cm  coefficients:'/(3x,3(' C',I2,'=',1PD15.7:)))
  625 FORMAT(/' *** ERROR *** The number of parameters',I3,'  does not e
     1qual the the number needed for HFD-ABC or HFD-ID')
  628 FORMAT(/' Potential is Generalized HFD-',A3,' with exponent factor
     1s   gamma=',f9.6/'   beta1=',f12.8,'   beta2=',f10.6,5x,'A=',
     2 1PD16.9:"  &  Cm's:"/(3x,3('   C',I2,' =',D15.7:)))
  629 FORMAT('   De=',f10.4,'[cm-1]   Re=',f9.6,'[Angst.]':'   and   for
     1  r <',F9.6/'     Damping function  D(r)= exp[ -',f6.4,'*(',f9.6, 
     2  '/r -1.0)**',f5.2,']')
 6299 FORMAT(15x,'and overall damping function:'/20x,'f2(r)= 1 - rhoAB*r
     1[bohr]^1.68 *exp{0.78*rhoAB*r[bohr]}'/)
c=======================================================================
c** Common uLR(r) printout for the MLR, DELR,  HPP, TT and HDF potentials
          IF((PSEL(ISTATE).GE.2).AND.(PSEL(ISTATE).NE.4)) THEN
c... first, specify choice of damping fx. {if damping included}
              IF(rhoAB(ISTATE).GT.0.d0) THEN
                  IF(IDSTT(ISTATE).GT.0) THEN
                      PVSR= DFLOAT(IVSR(ISTATE))*0.5d0 
                      WRITE(6,607) rhoAB(ISTATE),PVSR,bDS(IVSR(ISTATE)),
     1                                           cDS(IVSR(ISTATE)),PVSR
                    ELSE
                      LSR= IVSR(ISTATE)/2
                      IF(PSEL(ISTATE).NE.6) WRITE(6,617) rhoAB(ISTATE),
     1                                                    LSR, bTT(LSR)
                      IF(PSEL(ISTATE).EQ.6) WRITE(6,617) rhoAB(ISTATE),
     1                                                              LSR
                    ENDIF
                ELSE
                  WRITE(6,664)
cc                WRITE(6,608) MMLR(1,ISTATE),CmVAL(1,ISTATE),MMLR(1,ISTATE)
                ENDIF
c** List (inverse) power and coefficients of terms contributing to uLR(r)
c... First ... header (& lead coefft.) for all A-F diagonalization cases
              m1= 1
              IF(MMLR(1,ISTATE).LE.0) THEN
                  I4= MMLR(1,ISTATE)
                  IF(I4.GE.-1) WRITE(6,606) LAB4(I4),CmVAL(1,ISTATE) 
                  IF(I4.LE.-2) WRITE(6,6066) LAB4(I4),CmVAL(1,ISTATE)
                  m1=2
                  ENDIF
              DO  m= m1,NCMM(ISTATE)
                  IF(m1.GT.1) THEN           !! A-F cases
c... now,... Cm's for A-F 2x2 cases
                      IF(I4.GE.-1) WRITE(6,708) LAB2(m),CmVAL(m,ISTATE),
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
              IF((PSEL(ISTATE).EQ.7).AND.
     1                             (Nbeta(ISTATE).EQ.2)) WRITE(6,6299)
              IF(PSEL(ISTATE).EQ.2) WRITE(6,682) BINF
c** quadratic corrections for MLR
              IF(PSEL(ISTATE).EQ.2) THEN
c... First define 1D arrays for L-R powers & coefficients
                  DO m=1, NCMM(ISTATE)
                      MMLR1D(m)= MMLR(m,ISTATE)
                      CmVAL1D(m)= CmVAL(m,ISTATE)
                      CmEFF1D(m)= CmVAL1D(m)
                      ENDDO
                  CALL quadCORR(NCMM(ISTATE),MCMM,NCMMAX,MMLR1D,
     1                                     DE(ISTATE),CmVAL1D,CmEFF1D)
                  DO m=1, NCMM(ISTATE)
                      MMLR1D(m)= MMLR(m,ISTATE)
                      CmEFF(m,ISTATE)= CmEFF1D(m)
                      ENDDO
                  ENDIF
              ENDIF
c============== End of potential form header printout ==================
          IF((IOMEG(ISTATE).GE.0).AND.(PSEL(ISTATE).GE.0))
     1         WRITE(6,683) IOMEG(ISTATE), IOMEG(ISTATE)*IOMEG(ISTATE)
   50     IF((NUA(ISTATE).GE.0).OR.(NUB(ISTATE).GE.0)) THEN
c** Print description of 'adiabatic' BOB functional forms ...
              IF(BOBCN(ISTATE).GT.0) WRITE(6,556) qAD(ISTATE)
              IF(BOBCN(ISTATE).LE.0)WRITE(6,557) pAD(ISTATE),qAD(ISTATE)
              IF(NUA(ISTATE).GE.0) THEN
                  IF(BOBCN(ISTATE).GT.0) WRITE(6,564) '\tilde{S}(',
     1                               NAME(1),qAD(ISTATE),NUA(ISTATE)-1
                  IF(BOBCN(ISTATE).LE.0) WRITE(6,558) '\tilde{S}(',
     1       NAME(1),pAD(ISTATE),pAD(ISTATE),qAD(ISTATE),NUA(ISTATE)-1
                  WRITE(6,554) NAME(1),(qAD(ISTATE),i= 1,5)
                  IF(LRad(ISTATE).GT.0) THEN
                      WRITE(6,570) 
                      DO  m= 1,NCMM(ISTATE)
                          WRITE(6,571) MMLR(m,ISTATE),dCmA(m,ISTATE),
     1                                                  MMLR(m,ISTATE)
                          ENDDO
                      ENDIF
                  ENDIF
              IF(NUB(ISTATE).GE.0) THEN
                  IF(BOBCN(ISTATE).GT.0) WRITE(6,564) '\tilde{S}(',
     1                               NAME(2),qAD(ISTATE),NUB(ISTATE)-1
                  IF(BOBCN(ISTATE).LE.0) WRITE(6,558) '\tilde{S}(',
     1       NAME(2),pAD(ISTATE),pAD(ISTATE),qAD(ISTATE),NUB(ISTATE)-1
                  WRITE(6,554) NAME(2),(qAD(ISTATE),i= 1,5)
                  IF(LRad(ISTATE).GT.0) THEN
                      WRITE(6,570)
                      DO  m= 1,NCMM(ISTATE)
                          WRITE(6,571) MMLR(m,ISTATE),dCmB(m,ISTATE),
     1                                                  MMLR(m,ISTATE)
                          ENDDO
                      ENDIF
                  ENDIF
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
          IF(IOMEG(ISTATE).EQ.-2) WRITE(6,619) -IOMEG(ISTATE)
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
     2                            JTRUNC(ISTATE), EFSEL(ISTATE),ISTATE
                ELSE
                  WRITE(20,700) SLABL(ISTATE), IOMEG(ISTATE),
     1                      VMIN(ISTATE,1), VMAXIN(ISTATE), 
     2                            JTRUNC(ISTATE), EFSEL(ISTATE),ISTATE
                  WRITE(20,705) (VMAX(ISTATE,I), I=1,NISTP)
                ENDIF
              WRITE(20,701) PSEL(ISTATE),VLIM(ISTATE),MAXMIN(ISTATE),
     1                                      BOBCN(ISTATE),OSEL(ISTATE)
              WRITE(20,702) RMIN(ISTATE), RMAX(ISTATE), RH(ISTATE)
              IF((PSEL(ISTATE).GE.2).AND.(PSEL(ISTATE).LE.5)) THEN
                  WRITE(20,*)
                  WRITE(20,703) NCMM(ISTATE),rhoAB(ISTATE),
     1                                      IVSR(ISTATE),IDSTT(ISTATE)
                  IF(NCMM(ISTATE).GT.0) THEN
                      IF(MMLR(1,ISTATE).GT.0) THEN
                          DO I= 1,NCMM(ISTATE)
                              WRITE(20,704) MMLR(I,ISTATE),
     1                          CmVAL(I,ISTATE),IFXCM(I,ISTATE), I,I,I
                              ENDDO
                        ELSE
                          IF(MMLR(1,ISTATE).GE.-1) THEN
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
              ENDIF
c-----------------------------------------------------------------------
c** Write out PEF information for fitted  potentials.  1'st for De
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
          IF(PSEL(ISTATE).EQ.0) GO TO 60
          IF(PSEL(ISTATE).LE.4) THEN
              IPV= IPV + 1
              IPVRe(ISTATE) = IPV
              IF(IFXRE(ISTATE).LE.0) THEN
                  IF(DABS(RE(ISTATE)).GT.PU(IPV)) THEN
                      WRITE(6,620) 'Re',RE(ISTATE),PU(IPV),PS(IPV)
                    ELSE
                      WRITE(6,621) 'Re',RE(ISTATE),PU(IPV),PS(IPV)
                    ENDIF
                ELSE
                  WRITE(6,622) 'Re',RE(ISTATE)
                ENDIF
              IF(NPASS.GT.1) WRITE(20,670) RE(ISTATE),IFXRE(ISTATE),
     1                                                       'Re','Re'
              ENDIF
c-----------------------------------------------------------------------
c** Writing out the expansion variable information.
c-----------------------------------------------------------------------
          IF(PSEL(ISTATE).LE.4) THEN
              IPV= IPV + 1
              IF(IFXrefq(ISTATE).LE.0) THEN
                  IF(DABS(RREFq(ISTATE)).GT.PU(IPV)) THEN
                      WRITE(6,613) 'q',qPOT(ISTATE),RREFq(ISTATE),
     1                                                 PU(IPV),PS(IPV)
                    ELSE
                      WRITE(6,634) 'q',qPOT(ISTATE),RREFq(ISTATE),
     1                                                 PU(IPV),PS(IPV)
                    ENDIF
                ELSE
                  WRITE(6,635) 'q',qPOT(ISTATE),RREFq(ISTATE)
                ENDIF
              IF(NPASS.GT.1) WRITE(20,678) qPOT(ISTATE),RREFq(ISTATE),
     1                                     IFXrefq(ISTATE),('q',i=1,3)
  678  FORMAT(I3,1Pd17.10,0P,I3,9x,'% ',A1,'POT RREF',a1,' IFXref',A1) 
              ENDIF
c-----------------------------------------------------------------------
c** Writing out the switching function information.
c-----------------------------------------------------------------------
          IF(PSEL(ISTATE).EQ.2) THEN
              IPV= IPV + 1
              IF(IFXrefp(ISTATE).LE.0) THEN
                  IF(DABS(RREFp(ISTATE)).GT.PU(IPV)) THEN
                      WRITE(6,613) 'p',pPOT(ISTATE),RREFp(ISTATE),
     1                                                 PU(IPV),PS(IPV)
                    ELSE
                      WRITE(6,634) 'p',pPOT(ISTATE),RREFp(ISTATE),
     1                                                 PU(IPV),PS(IPV)
                    ENDIF
                ELSE
                  IF(RREFp(ISTATE).GT.0) THEN
                      WRITE(6,635) 'p',pPOT(ISTATE),RREFp(ISTATE)
                    ELSE  
                      WRITE(6,635) 'p',pPOT(ISTATE),RREFq(ISTATE)
                    ENDIF  
                ENDIF
              IF(NPASS.GT.1) WRITE(20,679) pPOT(ISTATE),RREFp(ISTATE),
     1                      IFXrefp(ISTATE),APSE(ISTATE), ('p',i=1,3)
  679  FORMAT(I3,1Pd17.10,0P,I3,I4,5x,'% ',A1,'POT RREF',a1,' IFXref',
     1  A1,' APSE')
              ENDIF
c
          IF((PSEL(ISTATE).GE.2).OR.(PSEL(ISTATE).EQ.3)) THEN
c-----------------------------------------------------------------------
c** For MLR or DELR, write out the  Cm  information.
c-----------------------------------------------------------------------
              IF(MMLR(1,ISTATE).LE.0) THEN
c** For Aubert-Frecon treatment of C3(r):C6(r) for alkali dimers
                  DO m=1,NCMM(ISTATE)
                      IPV= IPV+1
                      IF((MMLR(1,ISTATE).EQ.0)
     1                                .OR.(MMLR(1,ISTATE).EQ.-1)) THEN
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
c!! remove Cm's from parameter count for TT. HFD, ... etc. PEFs
                ELSEIF(PSEL(ISTATE).LE.4) THEN
c ... For 'regular' MLJ or MLR or DELR cases ... count & print Cm's
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
!!           IPV= IPVRe(ISTATE)       ??? section needing Work !!!
!!           DO m= 1, NCMM(ISTATE)
!!               IF(IFXCm(m,ISTATE).GT.1) THEN      !! ???????? huh ??
!!                   IPV= IPV+1   !! don't increase count unless...
c... Print re. a fitted Cm value constrained to equal that from another
c   state (with smaller ISTATE).  Input value of IFXCm(m,ISTATE) is IPV
c   parameter-counter value for that earlier Cm value.  c NOTE !!!! Need
c   to fix ISTATE count label  !!!!!!!!!
!!                   DO JSTATE= ISTATE,1,-1
!!                       IF((IFXCm(m,ISTATE).LT.IPVRe(JSTATE)).AND.
!!   1                      (IFXCm(m,ISTATE).GT.IPVRE(JSTATE-1))) THEN
!!                           CmVAL(m,ISTATE)= CmVAL(m,JSTATE-1) 
!!                           WRITE(6,666) MMLR(m,ISTATE),IPV,
!!   1                         IFXCm(m,ISTATE),MMLR(m,ISTATE),JSTATE-1,
!!   2                                                  CmVAL(m,ISTATE)
!!                           ENDIF
!!                       ENDDO
!!                   ENDIF
!!               ENDDO
!!666 FORMAT('  Constrain C_',I1,' = PV(',i3,')  to equal fitted  PV('
!!   1    ,I3,') = C_',I1,'(ISTATE=',I2,')'/53x,'=',1Pd14.7)
              ENDIF
c-----------------------------------------------------------------------
c** For DELR, calculate and write out the  A  and  B  coefficients
c-----------------------------------------------------------------------
          IF(PSEL(ISTATE).EQ.3) THEN
              yqRe=(RE(ISTATE)**qPOT(ISTATE) - AREF**qPOT(ISTATE))
     1                  /(RE(ISTATE)**qPOT(ISTATE) + AREF**qPOT(ISTATE))
              betaRe= beta(0,ISTATE)
              yPOW= 1.d0
              DO i= 1, Nbeta(ISTATE)
                  YPOW= YPOW*yqRe
                  betaRe= betaRe+ YPOW*beta(I,ISTATE)
                  ENDDO
              DO m=1,NCMM(ISTATE)
                  MMLR1D(m)= MMLR(m,ISTATE) 
                  ENDDO
              CALL DAMPF(RE(ISTATE),rhoAB(ISTATE),NCMM(ISTATE),NCMMAX,
     1                  MMLR1D,IVSR(ISTATE),IDSTT(ISTATE),DM,DMP,DMPP)
              IF(MMLR1D(1).LE.0) THEN
                  CALL AFdiag(RE(ISTATE),NCMM(ISTATE),NCMMax,MMLR1D,
     1              CmEFF1D,rhoAB(ISTATE),IVSR(ISTATE),IDSTT(ISTATE),
     2                                        VATTRe,dULRdCm,dVATTRE)
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
          IF(NPASS.GT.1) WRITE(20,671) Nbeta(ISTATE), APSE(ISTATE)
          J=0
          IF((PSEL(ISTATE).EQ.2).AND.(APSE(ISTATE).GT.0)) J=1
          IF(PSEL(ISTATE).GE.6) J=1
          DO  I=J, Nbeta(ISTATE)
              IPV= IPV + 1
              IF(IFXBETA(I,ISTATE).LE.0) THEN
                  IF(DABS(BETA(I,ISTATE)).GT.PU(IPV)) THEN
                      IF(APSE(ISTATE).LE.0) WRITE(6,640) 'be','ta',I,
     1                                   BETA(I,ISTATE),PU(IPV),PS(IPV)
                      IF(APSE(ISTATE).GT.0) WRITE(6,640) 'be','ta',I,
     1                  BETA(I,ISTATE),PU(IPV),PS(IPV),yqBETA(I,ISTATE)
                    ELSE
                      IF(APSE(ISTATE).LE.0) WRITE(6,641) 'be','ta',I,
     1                                   BETA(I,ISTATE),PU(IPV),PS(IPV)
                      IF(APSE(ISTATE).GT.0) WRITE(6,641) 'be','ta',I,
     1                  BETA(I,ISTATE),PU(IPV),PS(IPV),yqBETA(I,ISTATE)
                    ENDIF
                ELSE
                  IF(APSE(ISTATE).LE.0) WRITE(6,638) 'be','ta',I,
     1                                                   BETA(I,ISTATE)
                  IF(APSE(ISTATE).GT.0) WRITE(6,638) 'be','ta',I,
     1                                  BETA(I,ISTATE),yqBETA(I,ISTATE)
                ENDIF
              IF(NPASS.GT.1) THEN
                  IF(APSE(ISTATE).LE.0) THEN
                      WRITE(20,669) BETA(I,ISTATE),IFXBETA(I,ISTATE),
     1                                               'BETA',I,'BETA',I
                    ELSE
                      WRITE(20,668) yqBETA(I,ISTATE),BETA(I,ISTATE),
     1                                                IFXBETA(I,ISTATE)
                    ENDIF
                  ENDIF
              BTEMP= BTEMP+ BETA(I,ISTATE)
              ENDDO
c???       IF(PSEL(ISTATE).EQ.4) THEN   ! HPP ??why bother doing it here?
c              DO I= Nbeta(ISTATE)+1, Nbeta(ISTATE)+3 
c                  IPV= IPV+1 
c                  IF(IFXBETA(I,ISTATE).LE.0) THEN
c                        IF(DABS(BETA(I,ISTATE)).GT.PU(IPV)) THEN
c                            WRITE(6,640) 'pa','rm',I,BETA(I,ISTATE),
c     1                                                PU(IPV),PS(IPV) 
c                            WRITE(6,641) 'pa','rm',I,BETA(I,ISTATE),
c     1                                                PU(IPV),PS(IPV)
c                          ENDIF c                      ELSE
c                        WRITE(6,638) 'pa','rm',I,BETA(I,ISTATE)
c                      ENDIF
c                   ENDDO
c               ENDIF
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                  if(npass.gt.1) then
cc         write(6,500) DE(ISTATE),ddedre,(j,dDe(j),j=0,nbeta(istate)) 
cc500 FORMAT(/'    De=',f12.7,'     d{De}/d{Re}=',1p,d12.4,'   and'/
cc   1    3('  d{De}/db(',i2,')=',d12.4):)
cc         write(6,510) CMMp2,dCmp2dRe, (j,dCmp2(j),j=0,nbeta(istate))
cc510 FORMAT(/' CMMp2=',1P,d14.7,'   d{CMMp2}/d{Re}=',D12.4,'   and'/
cc   1    3(' d{Cmmp2}/db(',i2,')=',d12.4):)
cc                       endif
ccccccccccccccccccccccccccccccccccCCcccccccccccccccccccccccccccccccccccc
c c** Write out  phi_\infty  constant for the EMO or DELR forms
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
     1                                        *RE(ISTATE)**pPOT(ISTATE)
                  WRITE(6,652) MMLR(1,ISTATE)+pPOT(ISTATE),BTEMP
                ELSE
                  BTEMP= CmVAL(2,ISTATE)*2.d0*(2.d0*BINF - BTEMP)
     1                                        *RE(ISTATE)**pPOT(ISTATE)
                  WRITE(6,652) MMLR(2,ISTATE)+pPOT(ISTATE),BTEMP
                ENDIF
               ENDIF
c-----------------------------------------------------------------------
c** Writing out the adiabatic BOB radial function for atom A.
c-----------------------------------------------------------------------
   60     IF(NPASS.GT.1) THEN      !! next 4 - to stablize printout
               NUApr= NUA(ISTATE)-1 
               NUBpr= NUB(ISTATE)-1 
               NTApr= NTA(ISTATE)-1 
               NTBpr= NTB(ISTATE)-1 
               IF(NUA(ISTATE).LT.0) NUApr= -1 
               IF(NUB(ISTATE).LT.0) NUBpr= -1 
               IF(NTA(ISTATE).LT.0) NTApr= -1 
               IF(NTB(ISTATE).LT.0)  NTBpr= -1
               WRITE(20,672) NUApr, NUBpr,qAD(ISTATE),pAD(ISTATE),
     1                                                    LRad(ISTATE)
               IF(LRad(ISTATE).GT.0) THEN
                   DO m= 1,NCMM(ISTATE)
                       WRITE(20,677) dCmA(m,ISTATE),'A',MMLR(m,ISTATE) 
                       ENDDO
                   DO m= 1,NCMM(ISTATE)
                       WRITE(20,677) dCmB(m,ISTATE),'B',MMLR(m,ISTATE) 
                       ENDDO
                   ENDIF
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
  670 FORMAT(1Pd20.12,0P,I3,9x,'% ',A2,' IFX',A2) 
  671 FORMAT(/2I3,I4,25x,'% Nbeta APSE') 
  672 FORMAT(/2I3,I4,I3,I5,14x,'% NUA NUB qAD pAD LRad') 
  673 FORMAT(/2I3,I4,I3,19x,'% NTA NTB qNA pNA')
  674 FORMAT(1Pd20.12,0P,I3,9x,'% ',a5,'  IFX',A5)
  675 FORMAT(/3I3,24x,'% NwCFT Pqw efREF') 
  677 FORMAT(F15.4,21x,'% dCm',A1,'('I2,')'  )
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
                  ENDDO
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
          IF(NPASS.GT.1) WRITE(20,673) NTApr, NTBpr,qNA(ISTATE),
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
  552 FORMAT(8x,'using radial expansion variable:   y_{q} = (R^',I1,
     1   ' - Re^',I1,')/(R^',I1,' + Re^',I1,')')
  553 FORMAT(6x,'and switching function variable:   y_{p} = (R^',I1,
     1  ' -',F6.3,'^',I1')/(R^',I1' +',F6.3,'^',I1,')'/15x,'while fixing
     2 RREFp= RREFq')
  551 FORMAT(6x,'and switching function variable:   y_{p} = (R^',I1,
     1  ' -',F6.3,'^',I1')/(R^',I1' +',F6.3,'^',I1,')')
  555 FORMAT(5x,'using radial expansion variable:    y_{q} = (R^',I1,
     1  ' -',F6.3,'^',I1')/(R^',I1' +',F6.3,'^',I1,')')
  554 FORMAT(8x,'with ',A2,'-atom radial expansion variable:   y',I1,
     1  ' = (R^',I1,' - Re^',I1,')/(R^',I1,' + Re^',I1,')')
  556 FORMAT(' Adiabatic BOB functions are simple power series in  y_'
     1  I1,'(r) scaled by  m_e/M(A):')
  557 FORMAT(' Adiabatic BOB functions with {p=',i2,', q=',i2,
     1  '} are scaled by  DELTA{M(A)}/M(A):')
  564 FORMAT(4x,A10,A2,';R) = \sum\{u_i * [y',I1,']^i}   for  i= 0 to',
     1 i3)
  558 FORMAT(4x,A10,A2,';R) = u(inf)*y',i1,' + (1 - y',I1, ')*\Sum{u_i *
     1 [y',I1,']^i}   for  i= 0 to',i3)
  559 FORMAT(' Non-Adiabatic centrifugal BOB fx. with {p=',i2,', q=',i2,
     1  '} are scaled by  M(1)/M(A):')
  560 FORMAT(' Non-Adiabatic centrifugal BOB fx are power series in y_',
     1  I1,'(r) scaled by  m_e/M(A)')
  570 FORMAT(8x,'but replace  u(inf)  with  u(inf) + Sum{dCm/r**m}  wher
     1e:')
  571 FORMAT(46x,'dC',I1,'=',1PD16.9,'[cm-1 Ang^',i1,']') 
  600 FORMAT(1X,39('=='))
  601 FORMAT(' All distinct levels of State ',A3,' fitted as independent
     1 term values')
  603 FORMAT(/' FIXED State ',A3,' potential defined by interpolating ov
     1er input turning points')
  684 FORMAT(/' For state ',A3,' fits represents level with Band Constan
     1ts'/1x,6('=='))
  685 FORMAT(' *',A6,I3,'; IS=',I2,')=',1PD20.12,3X,1PD8.1,6X,1PD8.1,
     1                                                            10x)
  686 FORMAT(2X,A6,I3,'; IS=',I2,')=',1PD20.12,3X,1PD8.1,6X,1PD8.1,10x)
  687 FORMAT(I4,10I4)
 6872 FORMAT(4x,10I4)
  688 FORMAT('  ')
  690 FORMAT(7x,'Parameter',10x,'Final Value',5x,'Uncertainty Sensitivit
     1y')
  602 FORMAT(/' State ',A3,' represented by an EMO(Nbeta=',I2,' q=',i2,
     1 ') potential defined in terms of'/1x,4('=='),  '  exponent coeffi
     2cient:   beta(R)= Sum{beta_i*y',i1,'^i}  for  i= 0 -',i3)
  604 FORMAT(/' State ',A3,' represented by an MLR(q=',i2,', p=',i2,
     1  ') potential defined in terms of')
  605 FORMAT(1x,4('=='), '  exponent coefficient:  beta(R)= betaINF*y',
     1 i1,' +(1-y',I1,')*Sum{beta_i*y',i1,'^i}'/62x,'for  i= 0 to',I3)
  610 FORMAT(/' For state ',A3,"  use Surkus' Generalized Potential Ener
     1gy Function GPEF with"/1x,6('=='),'  expansion vble:   y_q(r) = r^
     2',i1,' - re^',i1,')/(',F5.2,'*r^',i1,' +',F5.2,'*re^',i1,'p)')
 6102 FORMAT(' *** Input ERROR *** band constant specification  v=',I3,
     1  ' .NE.', I3)
  612 FORMAT(/' State ',A3,' represented by a DELR(N=',I2,'  q=',i2,')
     1potential with'/1x,4('=='),'   exponent coefficient:   beta(r)=
     Su 2m{beta_i*y',i1,'^i}'/48x,'with polynomial order   Nbeta=',I3)
  607 FORMAT(4x,'uLR inverse-power terms incorporate DS-type damping wit
     1h   rhoAB=',f10.7/11x,'defined to give very short-range  Dm(r)*Cm/
     2r^m  behaviour   r^{',SP,f4.1,'}'/8x,SS,'Dm(r)= [1 - exp(-',f5.3, 
     3 '(rhoAB*r)/m -',f6.4,'(rhoAB*r)^2/sqrt{m})]^{m',SP,F4.1,'}')
 6072 FORMAT(4x,18('-'),11(A5: ))
  617 FORMAT(4x,'uLR inverse-power terms incorporate TT-type damping wit
     1h   rhoAB=',f10.7/8x,'defined to give very short-range  Dm(r)*Cm/r
     2^m  behaviour   r^{',SP,I2,'}'/8x,'Dm(r)= [1 - exp(-bTT*r)*SUM{(bT
     3T*r)^k/k!}]   where   bTT= rhoAB': '*',SS,f6.3)
  680 FORMAT(1x,4('=='),'  exponent coefficient defined as a Pashov natu
     1ral spline through'/ 10x,i3,' specified points including the fixed
     2 value of  betaINF  at  yp= 1')
  682 FORMAT(20x,'These constants yield:   betaINF=',F14.10 ) 
  683 FORMAT(/4x,'Since this state has (projected) electronic angular mo
     1mentum  OMEGA=',I2/10x,'eigenvalue calculations use centrifugal po
     2tential  [J*(J+1) -',I2,']/r**2')
  608 FORMAT(46x,'C',I1,'=',1PD16.9,'[cm-1 Ang^',i1,']')
 6082 Format(3x,I4,9x,12I5)
  609 FORMAT(45x,'C',I2,'=',1PD16.9,'[cm-1 Ang^{',i2,'}]')
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
  619 FORMAT(/' Centrifugal BOB term adds ',I2,'*[\hbar^2/{2\mu r^2}]  t
     1o the PEF')
  620 FORMAT(6X,A2,3X,1PD21.12,3X,1PD8.1,6X,1PD8.1,10x)
  621 FORMAT(5X,'*',A2,3X,1PD21.12,3X,1PD8.1,6X,1PD8.1) 
  622 FORMAT(6X,A2,3X,1PD21.12,7X,'--',12X,'--') 
  613 FORMAT(' RREF(',A1,'=',I2,')',1PD21.12,3X,1PD8.1,6X,1PD8.1,10x)
  634 FORMAT('*RREF(',A1,'=',I2,')',1PD21.12,3X,1PD8.1,6X,1PD8.1,10x)
  635 FORMAT(' RREF(',A1,'=',I2,')',1PD21.12,7X,'--',12X,'--') 
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
cc   ,5x,'at   yp= 1.000 1000000')
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
  701 FORMAT(I3,F12.4,2I3,I4,6x,' % PSEL VLIM MAXMIN BOBCN OSEL') 
  702 FORMAT(2F7.2,F8.4,9x,' % RMIN RMAX RH')
  703 FORMAT(I3,F11.6,2I3,15x,' % NCMM RHOab IDF IDSTT') 
  704 FORMAT(I3,1P,D17.8,0P,I3,8x,' % MMLR(',I1,'), CmVAL(',I1,
     1  '), IFXCm(',I1,')')
  705 FORMAT(10I4)
  706 FORMAT(I3,1P,D17.8,0P,I3,8x,' % ',A10)
  708 FORMAT(39x,A10,'=',1PD14.7,'[cm-1 Ang^{',i2,'}]')
  720 FORMAT(2X,A10,1PD20.12,3X,1PD8.1,6X,1PD8.1,10x)
  721 FORMAT(1X,'*',A10,3X,1PD21.12,3X,1PD8.1,6X,1PD8.1) 
  722 FORMAT(2X,A10,1PD20.12,7X,'--',12X,'--')
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

