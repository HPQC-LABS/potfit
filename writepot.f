c***********************************************************************
      SUBROUTINE WRITEPOT(NPASS,SLABL,NAME,DECM,PV,PU,PS)
c***********************************************************************
c** Subroutine to print out complete description of the potential fx.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++                   Version of  18 November 2012
c               (just after removal of RREFad, RREFna & RREFw)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** On entry:
c    NPASS   is the number of times WRIDATA was called.
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
      INTEGER MMAX
      PARAMETER (MMAX= 20)
      CHARACTER*2  NAME(2),SLABL(-5:NSTATEMX)
      CHARACTER*6  BCNAM(8)
      CHARACTER*7  NAMEDBLE(2)
      INTEGER NPASS, NALL,ISTATE,IPV,I,ISOT,J,MMN,m,LSR,
     1                    JSTATE,IPVRe(NSTATEMX)
      REAL*8 DECM(NSTATEMX),BTEMP,UAT,UBT,BINF,RE3,RE6,RE8,T0,T1,ULRe,
     1 C3VAL,C6adj,C9adj,RET,RETSig,RETPi,RETp,RETm,DM(MMAX),DMP(MMAX),
     2 DMPP(MMAX),bTT(-2:2),cDS(-4:0),bDS(-4:0), PVSR,ATT
      DATA BCNAM/' Tv(v=',' Bv(v=','-Dv(v=',' Hv(v=',' Lv(v=',' Mv(v=',
     1    ' Nv(v=',' Ov(v='/
c** Damping function factors from Table 1 of Mol.Phys. 109, 435 (2011)
       DATA bTT/2.10d0,2.44d0,2.78d0,3.13d0,3.47d0/
       DATA bDS/2.50d0,2.90d0,3.3d0,3.69d0,3.95d0/
       DATA cDS/0.468d0,0.446d0,0.423d0,0.405d0,0.390d0/
       SAVE bTT, bDS, cDS
c
c** NLLSSRR variables used in output
c
      REAL*8 PV(NPARMX), PU(NPARMX), PS(NPARMX)
      DATA NAMEDBLE/'wLambda',' wSigma'/
c-----------------------------------------------------------------------
c** Writing out state specific information.
c-----------------------------------------------------------------------
      IPV= 0
      DO  90 ISTATE=1,NSTATES
          WRITE(6,600)
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
              WRITE(6,684) SLABL(ISTATE)
              IF(NPASS.GT.1) THEN
                  WRITE(6,690)
                  DO  ISOT= 1, NISTP
                      DO  I= VMIN(ISTATE,ISOT),VMAX(ISTATE,ISOT)
                          IF(NBC(I,ISOT,ISTATE).GT.0) THEN
                              DO  J= 1,NBC(I,ISOT,ISTATE)
                                  IPV= IPV+1
                                  WRITE(6,686) BCNAM(J),I,ISOT,PV(IPV),
     1                                                 PU(IPV),PS(IPV)
                                  ENDDO
                              ENDIF
                          ENDDO
                      WRITE(6,688)
                      ENDDO
                  ENDIF
              GOTO 90
              ENDIF
          IF(PSEL(ISTATE).EQ.1) THEN
c** Header printout for EMO potential
              WRITE(6,602) SLABL(ISTATE),Nbeta(ISTATE),nPB(ISTATE),
     1                                       nPB(ISTATE),Nbeta(ISTATE)
              IF(RREF(ISTATE).LE.0.d0) WRITE(6,553) 
              IF(RREF(ISTATE).GT.0.d0) WRITE(6,555) RREF(ISTATE),
     1                                                    RREF(ISTATE)
              ENDIF
          IF(PSEL(ISTATE).EQ.2) THEN
c** Header printout for MLR potential
              BINF= betaINF(ISTATE)
              WRITE(6,604) SLABL(ISTATE),nPB(ISTATE),nQB(ISTATE)
              IF(APSE(ISTATE).LE.0) THEN
                  WRITE(6,605) nPB(ISTATE),nPB(ISTATE),nQB(ISTATE),
     1                                                   Nbeta(ISTATE)
                ELSE
                  WRITE(6,680) Nbeta(ISTATE)
                  BETA(Nbeta(ISTATE),ISTATE)= BINF
                ENDIF
              IF(RREF(ISTATE).LE.0.d0) WRITE(6,553) 
              IF(RREF(ISTATE).GT.0.d0) WRITE(6,555) RREF(ISTATE),
     1                                                    RREF(ISTATE)
              IF(rhoAB(ISTATE).GT.0.d0) THEN
                  IF(IDSTT(ISTATE).GT.0) THEN
                      PVSR= DFLOAT(IDF(ISTATE))*0.5d0
                      WRITE(6,607) rhoAB(ISTATE),PVSR,bDS(IDF(ISTATE)),
     1                                           cDS(IDF(ISTATE)),PVSR
                    ELSE
                      LSR= IDF(ISTATE)/2
                      WRITE(6,617) rhoAB(ISTATE), LSR, bTT(LSR)
                    ENDIF
                ELSE
                  WRITE(6,664) 
                ENDIF
c** for Aubert-Frecon 2x2 or 3x3 diagonalization uLR(r) funcrions
              IF((NCMM(ISTATE).GE.4).AND.(MMLR(2,ISTATE).LE.0)) THEN
                  IF(MMLR(2,ISTATE).EQ.0) WRITE(6,606) 'A-State',
     1CmVAL(2,ISTATE),CmVAL(1,ISTATE),(CmVAL(i,ISTATE),i=3,NCMM(ISTATE))
                  IF(MMLR(2,ISTATE).EQ.-2) WRITE(6,606) 'b-State',
     1CmVAL(2,ISTATE),CmVAL(1,ISTATE),(CmVAL(i,ISTATE),i=3,NCMM(ISTATE))
                  IF(MMLR(2,ISTATE).EQ.-1) WRITE(6,6066)
     1CmVAL(2,ISTATE),CmVAL(1,ISTATE),(CmVAL(i,ISTATE),i=3,NCMM(ISTATE))
                ELSE
c** for 'regular inverse-power sum uLR(r)
                  DO  m= 1,NCMM(ISTATE)
                      IF(MMLR(m,ISTATE).LE.9) 
     1      WRITE(6,608) MMLR(m,ISTATE),CmVAL(m,ISTATE),MMLR(m,ISTATE)
                      IF(MMLR(m,ISTATE).GT.9) 
     1      WRITE(6,609) MMLR(m,ISTATE),CmVAL(m,ISTATE),MMLR(m,ISTATE)
                      ENDDO
                ENDIF
              WRITE(6,682) BINF
              ENDIF
c
          IF(PSEL(ISTATE).EQ.3) THEN
c** Header printout for DELR potential form ...
              WRITE(6,612) SLABL(ISTATE),Nbeta(ISTATE),nPB(ISTATE),
     1                                                   Nbeta(ISTATE)
              IF(RREF(ISTATE).LE.0.d0) WRITE(6,553) 
              IF(RREF(ISTATE).GT.0.d0) WRITE(6,555) RREF(ISTATE),
     2                                                    RREF(ISTATE)
              WRITE(6,624)  rhoAB(ISTATE)
              WRITE(6,628) NCMM(ISTATE),MMLR(1,ISTATE),CmVAL(1,ISTATE)
              IF(NCMM(ISTATE).GT.1) WRITE(6,629) (MMLR(J,ISTATE),
     1                             CmVAL(j,ISTATE),j= 2,NCMM(ISTATE))
              ENDIF
c
          IF(PSEL(ISTATE).EQ.4) THEN
c** Header printout for Tiemann polynomial potential ...
              WRITE(6,623) SLABL(ISTATE), BETA(Nbeta(ISTATE)+1,ISTATE),
     1        BETA(Nbeta(ISTATE)+2,ISTATE),BETA(Nbeta(ISTATE)+3,ISTATE),
     2                      (nPB(ISTATE)),BETA(Nbeta(ISTATE)+1, ISTATE)
              DO  m= 1,NCMM(ISTATE)
                  IF(MMLR(m,ISTATE).LE.9) 
     1      WRITE(6,608) MMLR(m,ISTATE),CmVAL(m,ISTATE),MMLR(m,ISTATE)
                  IF(MMLR(m,ISTATE).GT.9) 
     1      WRITE(6,609) MMLR(m,ISTATE),CmVAL(m,ISTATE),MMLR(m,ISTATE)
                  ENDDO
              ENDIF
c
          IF(PSEL(ISTATE).EQ.5) THEN
c** Header printout for Tang-Toennies type potential ...
c*** If getting BTT and ATT from D_e and r_e, ...
c             BTT= DATT/(VATT - DE(ISTATE))
c             ATT= (VATT - DE(ISTATE))*DEXP(BTT*RE(ISTATE))
c** if getting ATT and D_e from R_e and b, ...
c             BTT= beta(0,ISTATE)
c             ATT= -DTT*DEXP(-BTT*REQ)/BTT
c             DE(ISTATE)= VTT - ATT*DEXP(-BTT*REQ)
c             rhoAB(ISTATE)= BTT/3.13d0
c
              WRITE(6,626)  RE(ISTATE), beta(0,ISTATE)
  626 FORMAT(/' State ',A2,' represented by a Tang-Toennies type potenti
     1al'/'  Input values of r_e=',F9.6,'  and  b=',F9.6,'  Used to dete
     2rmine   D_e   and A_{TT}')
              ENDIF
c
          IF(PSEL(ISTATE).EQ.6) THEN
c** Header printout for Surkus GPEF potential form ...
              WRITE(6,610) SLABL(ISTATE),(nPB(ISTATE),i=1,3),
     1             AGPEF(ISTATE),nPB(ISTATE),BGPEF(ISTATE),nPB(ISTATE)
              ENDIF
c
          IF((NUA(ISTATE).GE.0).OR.(NUB(ISTATE).GE.0)) THEN
c** Print description of 'adiabatic' BOB functional forms ...
              IF(BOBCN(ISTATE).LE.0) THEN
                  WRITE(6,557) 
                  IF(NUA(ISTATE).GE.0) THEN
                     WRITE(6,558) NAME(1),pAD(ISTATE),pAD(ISTATE),
     1                                       qAD(ISTATE),NUA(ISTATE)-1
                      WRITE(6,554) NAME(1),(qAD(ISTATE),i= 1,5)
                      ENDIF
                  IF(NUB(ISTATE).GE.0) THEN
                     WRITE(6,558) NAME(2),pAD(ISTATE),pAD(ISTATE),
     1                                       qAD(ISTATE),NUB(ISTATE)-1
                      WRITE(6,554) NAME(2),(qAD(ISTATE),i= 1,5)
                      ENDIF
                ELSE
                  WRITE(6,564) (qAD(ISTATE),i=1,5)
                ENDIF
              ENDIF
          IF((NTA(ISTATE).GE.0).OR.(NTB(ISTATE).GE.0)) THEN
c** Print description of centrifugal BOB functional forms ...
              IF(BOBCN(ISTATE).LE.0) THEN
                  WRITE(6,559) 
                  IF(NTA(ISTATE).GE.0) THEN
                      WRITE(6,560) NAME(1),pNA(ISTATE),pNA(ISTATE),
     1                                       qNA(ISTATE),NTA(ISTATE)-1
                      WRITE(6,554) NAME(1),(qNA(ISTATE),i= 1,5)
                      ENDIF
                  IF(NTB(ISTATE).GE.0) THEN
                      WRITE(6,560) NAME(2),pNA(ISTATE),pNA(ISTATE),
     1                                       qNA(ISTATE),NTB(ISTATE)-1
                      WRITE(6,554) NAME(2),(qNA(ISTATE),i=1,5)
                      ENDIF
                ELSE
                  WRITE(6,562) (qNA(ISTATE),i=1,5)
                ENDIF
              ENDIF
          IF(NwCFT(ISTATE).GE.0) THEN
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
c
c** Write out headings for parameter list
          IF(NPASS.EQ.1) WRITE(6,614)
          IF(NPASS.EQ.2) WRITE(6,615)
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
              IF(IFXDE(ISTATE).LE.0) UBT = PU(IPV+1)
              UBT = (UAT+UBT)*DSQRT(DECM(ISTATE)+1.0d0)
              UAT = DE(1) - VLIM(1) + VLIM(ISTATE) - DE(ISTATE)
              WRITE(6,620) 'Te',UAT,UBT,0.0d0,'Te'
              END IF
c-----------------------------------------------------------------------
c** Writing out the De information.
c-----------------------------------------------------------------------
          IF((PSEL(ISTATE).GE.1).AND.(PSEL(ISTATE).LE.5)) THEN
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
          IF((PSEL(ISTATE).GE.2).AND.(PSEL(ISTATE).LE.5)) THEN
c-----------------------------------------------------------------------
c** For MLR or DELR, write out the  Cn  information.
c-----------------------------------------------------------------------
              IF((NCMM(ISTATE).GE.4).AND.(MMLR(2,ISTATE).EQ.0)) THEN
c** For Aubert-Frecon treatment of C3(r):C6(r) for alkali dimers
                  IPV= IPV+ 1
                  IF(IFXCm(1,ISTATE).LE.0) THEN
                      IF(DABS(CmVAL(1,ISTATE)).GT.PU(IPV))
     1                  WRITE(6,620)'C3',CmVAL(1,ISTATE),PU(IPV),PS(IPV)
                      IF(DABS(CmVAL(1,ISTATE)).LE.PU(IPV))
     1                  WRITE(6,621)'C3',CmVAL(1,ISTATE),PU(IPV),PS(IPV)
                    ELSE
                      WRITE(6,622) 'C3',CmVAL(1,ISTATE)
                    ENDIF
                  IPV= IPV+2
c ... now write out Aubert-Frecon C6 coefficients for alkali dimers
                  IF(IFXCm(3,ISTATE).LE.0) THEN
                      IF(DABS(CmVAL(3,ISTATE)).GT.PU(IPV)) THEN
                          WRITE(6,630) 'C6sig',CmVAL(3,ISTATE),
     1                                                 PU(IPV),PS(IPV)
                        ELSE
                          WRITE(6,631) 'C6sig',CmVAL(3,ISTATE),
     1                                                 PU(IPV),PS(IPV)
                        ENDIF
                    ELSE
                      WRITE(6,632) 'C6sig',CmVAL(3,ISTATE)
                    ENDIF
                  IPV= IPV+1
                  IF(IFXCm(4,ISTATE).LE.0) THEN
                      IF(DABS(CmVAL(4,ISTATE)).GT.PU(IPV)) THEN
                          WRITE(6,630) 'C6_pi',CmVAL(4,ISTATE),
     1                                                 PU(IPV),PS(IPV)
                        ELSE
                          WRITE(6,631) 'C6_pi',CmVAL(4,ISTATE),
     1                                                 PU(IPV),PS(IPV)
                        ENDIF
                    ELSE
                      WRITE(6,632) 'C6_pi',CmVAL(4,ISTATE)
                    ENDIF
                  IF(NCMM(ISTATE).GT.4) THEN
c ... if needed write out Aubert-Frecon C8 coefficients for alkali dimers
                      IPV= IPV+ 1
                      IF(IFXCm(5,ISTATE).LE.0) THEN
                          IF(DABS(CmVAL(5,ISTATE)).GT.PU(IPV)) THEN
                              WRITE(6,630) 'C8sig',CmVAL(5,ISTATE),
     1                                                 PU(IPV),PS(IPV)
                            ELSE
                              WRITE(6,631) 'C8sig',CmVAL(5,ISTATE),
     1                                                 PU(IPV),PS(IPV)
                            ENDIF
                        ELSE
                          WRITE(6,632) 'C8sig',CmVAL(5,ISTATE)
                        ENDIF
                      IPV= IPV+1
                      IF(IFXCm(6,ISTATE).LE.0) THEN
                          IF(DABS(CmVAL(6,ISTATE)).GT.PU(IPV)) THEN
                              WRITE(6,630) 'C8_pi',CmVAL(6,ISTATE),
     1                                                 PU(IPV),PS(IPV)
                            ELSE
                              WRITE(6,631) 'C8_pi',CmVAL(6,ISTATE),
     1                                                 PU(IPV),PS(IPV)
                            ENDIF
                        ELSE
                          WRITE(6,632) 'C8_pi',CmVAL(6,ISTATE)
                        ENDIF
                      ENDIF
                ELSE
c ... For 'regular' MLJ or MLR or DELR or Tiemann case ...
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
                ENDIF
c** Check & do printouts re. Cm values constrained in fits
             IPV= IPVRe(ISTATE)
             DO m= 1, NCMM(ISTATE)
                 IPV= IPV+1
                 IF(IFXCm(m,ISTATE).GT.1) THEN
c... Print re. a fitted Cm value constrained to equal that from another
c   state (with smaller ISTATE).  Input value of IFXCm(m,ISTATE) is IPV
c   parameter-counter value for that earlier Cm value.
c NOTE !!!! Need to fix ISTATE count label  !!!!!!!!!
                     DO JSTATE= ISTATE,1,-1
                         IF((IFXCm(m,ISTATE).LT.IPVRe(JSTATE)).AND.
     1                      (IFXCm(m,ISTATE).GT.IPVRE(JSTATE-1))) THEN
                             CmVAL(m,ISTATE)= CmVAL(m,JSTATE-1)
                             WRITE(6,666) MMLR(m,ISTATE),IPV,
     1                         IFXCm(m,ISTATE),MMLR(m,ISTATE),JSTATE-1,
     2                         CmVAL(m,ISTATE)
                             ENDIF
                         ENDDO
                     ENDIF
                 ENDDO
  666 FORMAT('  Constrain C_',I1,' = PV(',i3,')  to equal fitted  PV('
     1    ,I3,') = C_',I1,'(ISTATE=',I2,')'/53x,'=',1Pd14.7)
              ENDIF 
c-----------------------------------------------------------------------
c** For DELR write out the  A  and  B  coefficients 
c-----------------------------------------------------------------------
          IF(PSEL(ISTATE).EQ.3) THEN
              WRITE(6,633) 'A(DELR)',AA(ISTATE)
              WRITE(6,633) 'B(DELR)',BB(ISTATE)
              WRITE(6,662) (MMLR(m,ISTATE),CmVAL(m,ISTATE), 
     1                                              m= 1,NCMM(ISTATE))
              ENDIF
c-----------------------------------------------------------------------
c** Writing out the exponent expansion parameter information.
c-----------------------------------------------------------------------
          BTEMP= 0.d0
          IF(NPASS.GT.1) WRITE(20,671) APSE(ISTATE), Nbeta(ISTATE),
     1                          nPB(ISTATE), nQB(ISTATE), RREF(ISTATE)
          NALL= Nbeta(ISTATE)
          J=0
          IF(APSE(ISTATE).GT.0) J=1
          DO  I=J, NALL
              IPV= IPV + 1
              IF(IFXBETA(I,ISTATE).LE.0) THEN
                  IF(DABS(BETA(I,ISTATE)).GT.PU(IPV)) THEN
                      IF(APSE(ISTATE).LE.0) WRITE(6,640) 'be','ta',I,
     1                                   BETA(I,ISTATE),PU(IPV),PS(IPV)
                      IF(APSE(ISTATE).GT.0) WRITE(6,640) 'be','ta',I,
     1                  BETA(I,ISTATE),PU(IPV),PS(IPV),ypBETA(I,ISTATE)
                    ELSE
                      IF(APSE(ISTATE).LE.0) WRITE(6,641) 'be','ta',I,
     1                                   BETA(I,ISTATE),PU(IPV),PS(IPV)
                      IF(APSE(ISTATE).GT.0) WRITE(6,641) 'be','ta',I,
     1                  BETA(I,ISTATE),PU(IPV),PS(IPV),ypBETA(I,ISTATE)
                    ENDIF
                ELSE
                  IF(APSE(ISTATE).LE.0) WRITE(6,638) 'be','ta',I,
     1                                                   BETA(I,ISTATE)
                  IF(APSE(ISTATE).GT.0) WRITE(6,638) 'be','ta',I,
     1                                  BETA(I,ISTATE),ypBETA(I,ISTATE)
                ENDIF
              IF(NPASS.GT.1) THEN
                  IF(APSE(ISTATE).LE.0) THEN
                      WRITE(20,669) BETA(I,ISTATE),IFXBETA(I,ISTATE),
     1                                               'BETA',I,'BETA',I
                    ELSE
                      WRITE(20,668) ypBETA(I,ISTATE),BETA(I,ISTATE),
     1                                                IFXBETA(I,ISTATE)
                    ENDIF
                  ENDIF
              BTEMP= BTEMP+ BETA(I,ISTATE)
              ENDDO

c???       IF(PSEL(ISTATE).EQ.4) THEN                ! ???????
c              DO I= NALL+1, NALL+3
c                    IPV= IPV+1
c                    IF(IFXBETA(I,ISTATE).LE.0) THEN
c                        IF(DABS(BETA(I,ISTATE)).GT.PU(IPV)) THEN
c                           IF(APSE(ISTATE).GE.0) 
c     1        WRITE(6,640) 'pa','rm',I,BETA(I,ISTATE),PU(IPV),PS(IPV)
c                           IF(APSE(ISTATE).LT.0) 
c     1        WRITE(6,640) 'pa','rm',I,BETA(I,ISTATE),PU(IPV),PS(IPV),
c     2                                               ypBETA(I,ISTATE)
c                          ELSE
c                           IF(APSE(ISTATE).GE.0) 
c     1        WRITE(6,641) 'pa','rm',I,BETA(I,ISTATE),PU(IPV),PS(IPV)
c                           IF(APSE(ISTATE).LT.0) 
c     1        WRITE(6,641) 'pa','rm',I,BETA(I,ISTATE),PU(IPV),PS(IPV),
c     2                                               ypBETA(I,ISTATE)
c                          ENDIF
c                      ELSE
c                        IF(APSE(ISTATE).GE.0) 
c     1                        WRITE(6,638) 'pa','rm',I,BETA(I,ISTATE)
c                        IF(APSE(ISTATE).LT.0) 
c     1       WRITE(6,638) 'pa','rm',I,BETA(I,ISTATE),ypBETA(I,ISTATE)
c                      ENDIF
c                   ENDDO
c               ENDIF

          IF(PSEL(ISTATE).EQ.1) BINF= BTEMP
c** Write out  phi_\infty  constant for the MLR/MLJ form
          IF(PSEL(ISTATE).EQ.2) THEN
               IF(APSE(ISTATE).GE.0) WRITE(6,648) BINF
               IF(APSE(ISTATE).LT.0) WRITE(6,649) BINF
               BTEMP= CmVAL(1,ISTATE)*2.d0*(2.d0*BINF - BTEMP)
     1                                        *RE(ISTATE)**nPB(ISTATE)
               WRITE(6,652) MMLR(1,ISTATE)+nPB(ISTATE),BTEMP
               ENDIF
c-----------------------------------------------------------------------
c** Writing out the adiabatic BOB radial function for atom A.
c-----------------------------------------------------------------------
          IF(NPASS.GT.1) WRITE(20,672) NUA(ISTATE)-1, NUB(ISTATE)-1,
     1                 pAD(ISTATE),qAD(ISTATE),'NUA','NUB','pAD','qAD'
          IF(NUA(ISTATE).GE.0) THEN
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
                  END DO
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
  670 FORMAT(1Pd20.12,0PI3,9x,'% ',A2,'  IFX',A2)
  671 FORMAT(/2I3,I4,I3,1PD11.2,8x,'% APSE Nbeta nPB nQB RREF')
  672 FORMAT(/2I3,I4,I3,19x,'% ',A3,1x,A3,2x,A3,1x,A3)
  674 FORMAT(1Pd20.12,0P,I3,9x,'% ',a5,'  IFX',A5)
  675 FORMAT(/3I3,24x,'% NwCFT Pqw efREF')
c-----------------------------------------------------------------------
c** Writing out the adiabatic BOB radial function for atom B.
c-----------------------------------------------------------------------
          IF(NUB(ISTATE).GE.0) THEN
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
          IF(NPASS.GT.1) WRITE(20,672) NTA(ISTATE)-1, NTB(ISTATE)-1,
     1                 pNA(ISTATE),qNA(ISTATE),'NTA','NTB','pNA','qNA'
          IF(NTA(ISTATE).GE.0) THEN
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
          IF(NTB(ISTATE).GE.0) THEN
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
          IF(NwCFT(ISTATE).GE.0) THEN
              IF(NPASS.GT.1) WRITE(20,675) NwCFT(ISTATE),Pqw(ISTATE),
     1                                     efREF(ISTATE)
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
  553 FORMAT(11x,'with radial expansion variable:   yq = (R^q - Re^q)/(R
     1^q + Re^q)')
  554 FORMAT(7x,'with ',A2,'-atom radial expansion variable:   y',I1,
     1  ' = (R^',I1,' - Re^',I1,')/(R^',I1,' + Re^',I1,')')
  555 FORMAT(8x,'with radial variable:   y_{p,q} = (R^q -',F9.6,'^q)/(R^
     1q +',F9.6,'^q)')
  557 FORMAT(' Adiabatic BOB functions are scaled by  DELTA{M(A)}/M(A)',
     1 '  and expanded as:')
  564 FORMAT(' Adiabatic BOB functions are scaled by  m_e/M(A)  and expa
     1nded as'/5x,' simple power series in   y',i1,'= (R^',i1,' - Re^',
     2 i1,')/(R^',i1,' + Re^',i1,')' )
  558 FORMAT(5x,'dV_ad(',A2,';R) = u(inf)*y',i1,' + (1 - y',I1,')*Sum{u_
     1i * [y',I1,']^i}   for  i= 0 to',i3)
  559 FORMAT(' Non-Adiabatic centrifugal BOB fx. are scaled by  M(1)/M(A
     1)  and expanded as:')
  560 FORMAT(5x,'t_nad(',A2,';r) = t(inf)*y',i1,' + (1 - y',I1,')*Sum{t_
     1i * [y',I1,']^i}   for  i= 0 to',i3)
  562 FORMAT(' Non-Adiabatic centrifugal BOB fx. are scaled by  m_e/M(A)
     1  and expanded as'/5x,' simple power series in   y_',i1, 
     2 '(r) = (r^',i1,' - re^',i1,')/(r^',i1,' + re^',i1,')' )
  600 FORMAT(1X,39('=='))
  601 FORMAT(' All distinct levels of State ',A2,' fitted as independent
     1 term values')
  603 FORMAT(/' FIXED State ',A2,' potential defined by interpolating ov
     1er input turning points')
  684 FORMAT(/' For state ',A2,' fits represents level with Band Constan
     1ts'/1x,6('=='))
  686 FORMAT(2X,A6,I3,'; IS=',I2,')=',1PD20.12,3X,1PD8.1,6X,1PD8.1,10x)
  688 FORMAT('  ')
  690 FORMAT(7x,'Parameter',10x,'Final Value',5x'Uncertainty   Sensitivi
     1ty')
  602 FORMAT(/' State ',A2,' represented by an EMO(N=',I2,' q=',i2,') po
     1tential defined in terms of'/1x,4('=='),  '  exponent coefficient:
     2   beta(R)= Sum{beta_i*y',i1,'^i}  for  i= 0 -',i3)
  604 FORMAT(/' State ',A2,' represented by an MLR(p=',i2,', q=',i2,
     1  ') potential defined in terms of')
  605 FORMAT(1x,4('=='), '  exponent coefficient:  beta(R)= betaINF*y',
     1 i1,' +(1-y',I1,')*Sum{beta_i*y',i1,'^i}'/64x,'for  i= 0 -',I3)
  610 FORMAT(/' For state ',A2,"  use Surkus' Generalized Potential Ener
     1gy Function GPEF with"/1x,6('=='),'  expansion vble:   y_',i1, 
     2 '(r) = (r^',i1,' - re^',i1,')/(',F5.2,'*r^',i1,' +',F5.2,'*re^',
     3 i1,'p)')
  612 FORMAT(/' State ',A2,' represented by a DELR(N=',I2,'  q=',i2,') 
     1potential with an exponent'/1x,4('=='),   '  exponent coefficient:
     2   beta(R)= Sum{beta_i*y',i1,'^i}  for  i= 0 -',i3)
  607 FORMAT(4x,'uLR inverse-power terms incorporate DS-type damping wit
     1h   rhoAB=',f10.7/8x,'defined to give very short-range  Dm(r)*Cm/r
     2^m  behaviour   r^{',SS,f4.1,'}'/8x,'Dm(r)= [1 - exp(-',f5.2,
     3 '(rhoAB*r)/m -',f6.3,'(rhoAB*r)^2/sqrt{m})]^{m',SP,F4.1,'}')
  617 FORMAT(4x,'uLR inverse-power terms incorporate TT-type damping wit
     1h   rhoAB=',f10.7/8x,'defined to give very short-range  Dm(r)*Cm/r
     2^m  behaviour   r^{',I2,'}'/8x,'Dm(r)= [1 - exp(-bTT*r)*SUM{(bTT*r
     1)^k/k!}]   where   bTT=',f6.3,'*rhoAB') 
  680 FORMAT(1x,4('=='),'  exponent coefficient defined as a Pashov natu
     1ral spline through'/ 10x,i3,' specified points including the fixed
     2 value of  betaINF  at  yp= 1')
  682 FORMAT(20x,'These constants yield:   betaINF=',F14.10 )
  608 FORMAT(48x,'C',I1,'=',1PD14.7,'[cm-1 Ang^',i1,']')
  609 FORMAT(47x,'C',I2,'=',1PD14.7,'[cm-1 Ang^{',i2,'}]')
  606 FORMAT(5x,'Use Lyon 2x2  ',A7,'  uLR(r)  with   Aso=',F10.6/47x,
     1'C_3(1Sigma)='  ,1PD15.7:/47x,'C_6(Sigma) =',D15.7:/47x,'C_6(Pi)',
     2 '    =',D15.7:/  47x,'C_8(Sigma) =',1PD15.7:/47x,'C_8(Pi)    =',
     3  D15.7)
 6066 FORMAT(' Use Lyon 3x3  uLR(r)  with   Aso=',F10.6,'   C_3(Sigma)='
     1  ,1PD15.7:/47x,'C_6(Sigma)=',D15.7:/47x,'C_6(Pi)   =',D15.7:/
     2  47x,'C_8(Sigma)=',1PD15.7:/47x,'C_8(Pi)   =',D15.7)

  614 FORMAT(/'   Parameter    Initial Value    Uncertainty   Sensitivit
     1y')
  615 FORMAT(/'   Parameter      Final Value    Uncertainty   Sensitivit
     1y')
  618 FORMAT(1x,A6,'-doubling splitting strength function is expanded as
     1'/7x,'f(r) =  Sum{w_i * (y',i1,')^i}   for   i= 0 to',i3/
     2 11x,'with radial expansion variable:   y',I1,' = (R^',I1,
     3  ' - Re^',I1,')/(R^',I1,' + Re^',I1,')')
  619 FORMAT(' Centrifugal potential strength factor is  [J(J+1) +',I2,
     1   ']')
  620 FORMAT(6X,A2,3X,1PD21.12,3X,1PD8.1,6X,1PD8.1,10x)
  621 FORMAT(5X,'*',A2,3X,1PD21.12,3X,1PD8.1,6X,1PD8.1)
  622 FORMAT(6X,A2,3X,1PD21.12,7X,'--',12X,'--')
  623 FORMAT(/' State ',A2,' represented by a Tiemann polynomial (b=',
     1   F5.2,', R_in=',F4.2,', R_out=',F4.2,')'/1x,4('=='),5x,'with exp
     2ansion variable:  y_',i1,'(r) = (r - re)/(r ',F5.2,'*re)') 
  624 FORMAT(3x,'u_{LR} term uses Tang-Toennies damping function [JCP 80
     1,3726(1984)]'/5x,'[1 - exp(-3.13*RHO*r)*SUM{(3.13*RHO*r)^k/k!}]',
     2 '   with  rhoAB=',F9.6)
  628 FORMAT('   and',i2,' inverse-powers terms with coefficients:',
     1  '   C',0P,I2,'=',1Pd14.6)
  629 FORMAT((51x,'C',0P,I2,'=',1Pd14.6:)) 
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
  649 FORMAT(3X,'beta_INF',1PD21.12,7X,'--',12X,'--')
cc   ,5x,'at   yp=  1.000 1000000')
  650 FORMAT(3X,A2,A2,'(',I2,')',1PD21.12,7X,'--',12X,'--')
  651 FORMAT(1X,A6,'(',I2,')',1PD21.12,7X,'--',12X,'--')
  652 FORMAT('   C',I2,'{eff}',1PD21.12,7X,'--',12X,'--')
  654 FORMAT(1X,A6,'_inf',1PD21.12,3X,1PD8.1,6X,1PD8.1)
  656 FORMAT('*',A6,'_inf',1PD21.12,3X,1PD8.1,6X,1PD8.1)
  658 FORMAT(1X,A6,'_inf',1PD21.12,7X,'--',12X,'--') 
  660 FORMAT(5X,'C',I2,3X,1PD21.12,3X,1PD8.1,6X,1PD8.1) 
  661 FORMAT(3X,'* C',I2,3X,1PD21.12,3X,1PD8.1,6X,1PD8.1) 
  662 FORMAT(5X,'C',I2,3X,1PD21.12,7X,'--',12X,'--':)
  664 FORMAT(4x,'uLR inverse-power terms incorporate NO damping function
     1s')
  692 FORMAT(' ', A2,' state  Lambda-doubling split levels referenced to
     1 f-parity levels')
  694 FORMAT(' ', A2,' state  Lambda-doubling split levels referenced to
     1 the e/f level mid-point')
  696 FORMAT(' ', A2,' state  Lambda-doubling split levels referenced to
     1 e-parity levels')
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
