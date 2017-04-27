c***********************************************************************
      SUBROUTINE DYIDPJ(IDAT,NDATA,NPTOT,YOBS,YC,PV,PD)
c***********************************************************************
c** This program assumes that the upper state IEP is at a higher energy
c   than the lower state IEPP.
c** This subroutine returns the calculated value YC of datum IDAT, and
c  its partial derivatives PD(k) w.r.t. the NPTOT parameters PV.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++  COPYRIGHT 2007-16  by  R.J. Le Roy, Jenning Seto and Yiye Huang +++
c   Dept. of Chemistry, Univ. of Waterloo, Waterloo, Ontario, Canada   +
c    This software may not be sold or any other commercial use made    +
c      of it without the express written permission of the authors.    +
c++++++++++++++++++++ (version of 18/02/2016) ++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** On entry:
c    IDAT     is the number of the current observable being considered. 
c    NPTOT    is the total number of parameters in the model
c    PV(i)    is the array of parameters being varied.
c-----------------------------------------------------------------------
c** On exit:
c    YC       is the calculated value of the IDATth observable.
c    PV(i)    is the array of parameters being varied.
c    PD(i)    is the partial derivative array (de/dp).
c=======================================================================
      INCLUDE 'arrsizes.h'
      INCLUDE 'BLKISOT.h'
      INCLUDE 'BLKDATA.h'
      INCLUDE 'BLKPOT.h'
      INCLUDE 'BLKPARAM.h'
      INCLUDE 'BLKBOB.h'
      INCLUDE 'BLKCOUNT.h'
c-----------------------------------------------------------------------
c** Common block for partial derivatives of potential at the one distance RDIST
c   and HPP derivatives for uncertainties
      REAL*8 dVdPk(HPARMX),dDe(0:NbetaMX),dDedRe
      COMMON /dVdPkBLK/dVdPk,dDe,dDedRe
c=======================================================================
      INTEGER IDAT,NPTOT,NBAND,IISTP,ISTATE,I,J,NDATA,fcount
c
c** Define parameters required locally and from NLLSSRR.
      REAL*8 RDIST,VDIST,BETADIST,VLAST,EUP,ELW,YOBS,YC,EO,width,
     1  VMAXX(NSTATEMX),PV(NPARMX),PD(NPARMX),UPPER(NPARMX),
     2  LOWER(NPARMX),DEDPK(HPARMX),BVIR,dBVIRdP(NPARMX)
c
c----------------------------------------------------------------------- c
      INTEGER INNR(0:NVIBMX)
      SAVE VMAXX
c
      IF(IDAT.EQ.1) THEN
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++ At beginning of each fit cycle (datum #1), re-map internal NLLSSRR
c  parameter array PV onto external (physical) variable set, and get
c  updated band constant array ZK for estimating trial eigenvalues.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          fcount= 0
          CALL MAPPAR(NISTP,PV,1)
          REWIND(22)
          REWIND(7)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** If fitting to fluoresence data, convert fluoresence series origin
c   parameters back to external (logical) variable system.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          IF((NFSTOT.GT.0).OR.(NTVALL(0).GT.0)) THEN
              DO J= TOTPOTPAR+1,NPTOT
                  TVALUE(J)= PV(J)
                  ENDDO
              ENDIF
          DO  ISTATE= 1,NSTATES
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Call subroutine to update potential functions AND their partial 
c   derivatives w.r.t. fitting parameter for each state.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
              IF(PSEL(ISTATE).GT.0) CALL VGEN(ISTATE,-1.d0,VDIST,
     1                                                  BETADIST,IDAT)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Now generate band constants for each state and isotopologue for 
c  generating trial eigenvalues in fit calculations.  To take account of
c  the different vibrational ranges for different isotopologues, only 
c  generate values for levels to the input VMAX for isotopologue-1.
              IF(PSEL(ISTATE).GE.0) THEN
                  DO  IISTP= 1,NISTP
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Call subroutine INITDD to calculate band constants for initial
c   trial eigenvalue estimates as well as (? what) the partial derivatives.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                    CALL INITDD(ISTATE,IISTP,VMAX(ISTATE,IISTP),
     1                                              VMAXX(ISTATE),INNR)
                    ENDDO
                  ENDIF
              ENDDO
          ENDIF
c++++++++++++++end of datum-1 set-up ++++++++++end of datum-1 set-up +++
c** Initialize variables for current datum
      EUP= 0.0D0
      ELW= 0.0D0
      DO  I= 1,NPTOT
          PD(I)= 0.0d0
          UPPER(I)= 0.0d0
          LOWER(I)= 0.0d0
          ENDDO
      NBAND= IB(IDAT) 
      IISTP= ISTP(NBAND)
c
c** Now to determine partials with respect to the upper and lower
c
      IF(IEP(NBAND).EQ.0) THEN
c=======================================================================
c** For fluorescence series data ...
c=======================================================================
          I= NFS1-1
          DO J= NFS1,NPTOT
              IF(FSBAND(J-I).EQ.NBAND) THEN
c ... PV(HPARMX+I) is the energy of the I'th fluorescence band.
                  IF(IFXFS(NFS(NBAND)).LE.0) THEN
                      EUP= TVALUE(J)
                      UPPER(J)= 1.0d0
                    ELSEIF(FSSame.GT.0) THEN
c ... if this FS band shares its origin with some earlier FS band ...
                      EUP= TVALUE(IFXFS(NFS(NBAND)))
                      UPPER(IFXFS(NFS(NBAND)))= 1.d0
                    ENDIF
                  ENDIF
              ENDDO
          CALL DEDP(IDAT,IEPP(NBAND),IISTP,ZMASS(3,IISTP),JP(IDAT),
     1  JPP(IDAT),EFPP(IDAT),ELW,VMAXX(IEPP(NBAND)),width,LOWER,fcount)
              IF(ELW.LT.-9.d9) THEN
c... if eigenvalue search failed, remove this datumn from the fit
                  WRITE(6,600) SLABL(IEPP(NBAND)),JP(IDAT),JPP(IDAT),
     1                                                       IDAT,YOBS
                  YC= YOBS
                  RETURN
                  ENDIF
              IF(width.GT.0.d0) THEN
c*** Reduce weight of Qbdd levels assuming unc(Airy) = Fqb*width 
                  UFREQ(IDAT) = DSQRT(YUNC(IDAT)**2 + (Fqb*width)**2)
                  WRITE(22,623) 'LOWER',IDAT,FREQ(IDAT),VP(NBAND),
     1                           JP(IDAT),width,YUNC(IDAT),UFREQ(IDAT)
                  ENDIF
c
        ELSEIF((IEP(NBAND).EQ.-1).AND.(PSEL(IEPP(NBAND)).GE.0)) THEN
c=======================================================================
c*** For PAS data ...
c=======================================================================
          CALL DEDP(IDAT,IEPP(NBAND),IISTP,ZMASS(3,IISTP),JP(IDAT),
     1  JPP(IDAT),EFPP(IDAT),ELW,VMAXX(IEPP(NBAND)),width,LOWER,fcount)
          IF(ELW.LT.-9.d9) THEN
c... if eigenvalue search failed, remove this datum from the fit
              WRITE(6,600) SLABL(IEPP(NBAND)),JP(IDAT),JPP(IDAT),
     1                                                        IDAT,YOBS
              YC= YOBS 
              RETURN 
              ENDIF
          ISTATE= IEPP(NBAND) 
          EUP= VLIM(ISTATE)
c... As appropriate, add isotopic u_\infty adjustment to the binding energy
          IF(NUA(ISTATE).GE.0) EUP= EUP +
     1                        ZMUA(IISTP,ISTATE)*UA(NUA(ISTATE),ISTATE)
          IF(NUB(ISTATE).GE.0) EUP= EUP +
     1                        ZMUB(IISTP,ISTATE)*UB(NUB(ISTATE),ISTATE)
c
        ELSEIF((IEP(NBAND).EQ.-2).AND.(PSEL(IEPP(NBAND)).GE.0)) THEN
c=======================================================================
c*** If datum is width of a tunneling-predissociation quasibound level 
c ... for forward calculation, use widths from SCHRQ in DEDP
c=======================================================================
          CALL DEDP(IDAT,IEPP(NBAND),IISTP,ZMASS(3,IISTP),JP(IDAT),
     1  JPP(IDAT),EFPP(IDAT),EO,VMAXX(IEPP(NBAND)),width,DEDPK,fcount)
          IF(EO.LT.-9.d9) THEN
c... if eigenvalue search failed, remove this datumn from the fit
              WRITE(6,600) SLABL(IEPP(NBAND)),JP(IDAT),JPP(IDAT),
     1                                                       IDAT,YOBS
              YC= YOBS 
              RETURN 
              ENDIF
c ... otherwise ... calculate 'width' and its derivatives in DWDP
          IF(PSEL(IEPP(NBAND)).GT.0) THEN
              CALL DWDP(IDAT,IEPP(NBAND),FREQ(IDAT),ZMASS(3,IISTP),
     1                           JP(IDAT),JPP(IDAT),EO,width,DEDPK,PD)
              YC= width 
              RETURN 
              ENDIF
c
        ELSEIF((IEP(NBAND).EQ.-3).AND.(PSEL(IEPP(NBAND)).GT.0)) THEN
c=======================================================================
c*** If datum is the potential energy function at some specific distance
c=======================================================================
          ISTATE= IEPP(NBAND) 
          RDIST= TEMP(IDAT) 
          VLAST= VPOT(NPNTMX,ISTATE)
          CALL VGEN(IEPP(NBAND),RDIST,VDIST,BETADIST,IDAT) 
            YC= VDIST 
          DO  J= POTPARI(IEPP(NBAND)), POTPARF(IEPP(NBAND))
              PD(J)= dVdPk(J) 
              ENDDO
          VPOT(NPNTMX,ISTATE)= VLAST 
          RETURN
        ELSEIF((IEP(NBAND).EQ.-4).AND.(PSEL(IEPP(NBAND)).GT.0)) THEN
c=======================================================================
c*** For Pressure Virial coefficient data ...
c=======================================================================
          CALL DVIRDP(IDAT,IEPP(NBAND),ZMASS(3,IISTP),BVIR,PD)
          YC= BVIR
          RETURN
        ELSEIF((IEP(NBAND).EQ.-5).AND.(PSEL(IEPP(NBAND)).GT.0)) THEN
c=======================================================================
c*** For Acoustic Virial coefficient data ...
c=======================================================================
          CALL DVACDP(IDAT,IEPP(NBAND),ZMASS(3,IISTP),BVIR,PD)
          YC= BVIR
          RETURN
c 
        ELSEIF(IEP(NBAND).GT.0) THEN
c=======================================================================
c*** For 'normal' microwave, infrared, and electronic data,
c*** determine the partials for the upper and lower levels`
c=======================================================================
          IF(PSEL(IEP(NBAND)).EQ.-2) THEN
c... if upper state being represented by term values ...
              UPPER(TVUP(IDAT))= 1.d0
              EUP= TVALUE(TVUP(IDAT))
            ELSE
c... for normal case of UPPER level being represented by a potential
c=======================================================================
              CALL DEDP(IDAT,IEP(NBAND),IISTP,ZMASS(3,IISTP),VP(NBAND),
     1    JP(IDAT),EFP(IDAT),EUP,VMAXX(IEP(NBAND)),width,UPPER,fcount)
              IF(EUP.LT.-9.d9) THEN
c... if eigenvalue search failed, remove this datumn from the fit
                  WRITE(6,600) SLABL(IEP(NBAND)),VP(NBAND),JP(IDAT),
     1                                                        IDAT,YOBS
                  YC= YOBS
                  RETURN
                  ENDIF
              IF(width.GT.0.d0) THEN
c*** Reduce weight of Qbdd levels assuming unc(Airy) = Fqb*width 
                  UFREQ(IDAT) = DSQRT(YUNC(IDAT)**2 + (Fqb*width)**2)
                  WRITE(22,623) 'UPPER',IDAT,FREQ(IDAT),VP(NBAND),
     1                           JP(IDAT),width,YUNC(IDAT),UFREQ(IDAT)
                  ENDIF
  623 FORMAT('  ',A5', level of Datum(',I5,')=',F9.2,'   v=',I3,'  J=',
     1 I3,'  has width=',1Pd9.2/9x,'so increase datum uncertainty from',
     2   0Pf9.6,' to',f9.6)
            ENDIF
          IF(PSEL(IEPP(NBAND)).EQ.-2) THEN
c... if lower state being represented by term values ...
              LOWER(TVLW(IDAT))= +1.d0
              ELW= TVALUE(TVLW(IDAT))
            ELSE
c... for normal case of LOWER level being represented by a potential
              CALL DEDP(IDAT,IEPP(NBAND),IISTP,ZMASS(3,IISTP),
     1         VPP(NBAND),JPP(IDAT),EFPP(IDAT),ELW,VMAXX(IEPP(NBAND)),
     2                                             width,LOWER,fcount)
              IF(ELW.LT.-9.d9) THEN
c... if eigenvalue search failed, remove this datumn from the fit
                  WRITE(6,600) SLABL(IEPP(NBAND)),VPP(NBAND),JPP(IDAT),
     1                                                        IDAT,YOBS
                  YC= YOBS
                  RETURN
                  ENDIF
              IF(width.GT.0.d0) THEN
c*** Reduce weight of Qbdd levels assuming unc(Airy) = Fqb*width 
                  UFREQ(IDAT) = DSQRT(YUNC(IDAT)**2 + (Fqb*width)**2)
                  WRITE(22,623) 'LOWER',IDAT,FREQ(IDAT),VP(NBAND),
     1                           JP(IDAT),width,YUNC(IDAT),UFREQ(IDAT)
                  ENDIF
            ENDIF
        ENDIF
      DO  I= 1,NPTOT
          PD(I)= UPPER(I) - LOWER(I)
          ENDDO
c----------------------------------------------------------------------
c** Get calculated value for the IDAT'th observable from energy levels
c----------------------------------------------------------------------
      YC= EUP - ELW
cc    if((widthLW.GT.0.d0).OR.(widthUP.GT.0.d0)) THEN
cc        WRITE(22,622) IDAT,FREQ(IDAT),widthUP,widthLW
cc622 FORMAT('  Datum(',I5')=',f10.2,'   widthUP=',1Pd9.2,'   widthLW=',
cc   1  d9.2) 
cc        ENDIF
c----------------------------------------------------------------------
      RETURN
  600 FORMAT(' *** FAIL to find level(',A3,')   v=',I3,'   J=',I3,
     1  '  so ignore  YOBS(',i5,')=',f12.4)
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE INITDD(ISTATE,IISTP,VIMX,VMAXX,INNR)
c***********************************************************************
c** This subroutine updates the trial vibrational energies & rotational
c        constants on each iteration
c** On entry:   ISTATE  is the states being considered.
c               IISTP  is the isotopologue being considered.
c               VIMX   is the upper vibrational bound for this state.
c** On exit:   ZK (in BLKISOT) are the band constants for this state & isotope
c              VMAXX  is  MAX{barrier maximum, VLIM} for this state.
c** Internal:  RM2   is the  (1+ ZMTA*TAR + ZMTB*TBR)/r**2  array for 
c                    this state required by CDJOEL for CDC calculation
c=======================================================================
      INCLUDE 'arrsizes.h'
      INCLUDE 'BLKBOB.h'
      INCLUDE 'BLKBOBRF.h'
      INCLUDE 'BLKDATA.h'
      INCLUDE 'BLKISOT.h'
      INCLUDE 'BLKPOT.h'
c-----------------------------------------------------------------------
c
      INTEGER ISTATE,IISTP,VIMX,AFLAG,VIN,NBEG,NEND,WARN,IWR,LPRWF,I,J,
     1   INNODE,KV,KVtop,NCN
c
      REAL*8 GV(0:NVIBMX),RCNST(NROTMX),RR(NPNTMX),RM2(NPNTMX),
     1  V(NPNTMX),SWF(NPNTMX),FWHM,PMAX,BFCT,BvWN,RHSQ,C3gu,T3,VMAXX
c
      INTEGER INNR(0:NVIBMX)
      REAL*8 SWF2,qDBL,qFCT
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      VIN= VIMX
      AFLAG= 0
      WARN= 0
      IWR= -1           !!! should in general be -1
ccc   IWR=  3           !!! use this when trouble shooting level search
      LPRWF= 0
      INNODE= 1
c
c** Calculate potential scaling factors for ALF and CDJOEL
c
      RHSQ= RH(ISTATE)*RH(ISTATE) 
      BFCT= (ZMASS(3,IISTP)/16.857629206D0)*RHSQ
      BvWN= 16.857629206D0/ZMASS(3,IISTP)
      qFCT= RH(ISTATE)*BvWN**(2*IOMEG(ISTATE))
c
c** Now generate the BFCT-scaled and adiabatically corrected potential
c   for the current isotope for use in SCHRQ, plus the 1/R**2 array for
c   CDC calculations in CDJOEL
c
      DO  I=1,NDATPT(ISTATE)
          RR(I)= RD(I,ISTATE)
          V(I)= BFCT*(VPOT(I,ISTATE)+ ZMUA(IISTP,ISTATE)*UAR(I,ISTATE)
     1                             + ZMUB(IISTP,ISTATE)*UBR(I,ISTATE))
          RM2(I)= 1.d0/RD(I,ISTATE)**2
c** Special BOB correction for A-state Li2 and analogous cases. !!!!!!!!
          IF(IOMEG(ISTATE).EQ.-2) V(I)= V(I) + 2.d0*RHSQ*RM2(I)
          ENDDO 
      IF((NCMM(ISTATE).GE.3).AND.(MMLR(1,ISTATE).LE.0).AND.(AN(1).EQ.3)
     1        .AND.(AN(2).EQ.3).AND.(MN(1,IISTP).NE.MN(2,IISTP))) THEN
c** Add g/u symmetry breakdown correction for special case {6,7}Li2(A) !!
          C3gu= (2.d0/3.d0)*CmVAL(2,ISTATE)
          DO I= 1,NDATPT(ISTATE)
              T3= C3gu/RD(I,ISTATE)**3
              V(I)= V(I) + BFCT*(T3 - DSQRT(T3**2 + 0.03085959756d0))
              ENDDO
          ENDIF
      IF((NTA(ISTATE).GE.0).OR.(NTB(ISTATE).GE.0)) THEN
          DO  I= 1, NDATPT(ISTATE)
              RM2(I)= RM2(I)*(1.d0 + ZMTA(IISTP,ISTATE)*TAR(I,ISTATE) 
     1                             + ZMTB(IISTP,ISTATE)*TBR(I,ISTATE))
              ENDDO
          ENDIF
      NCN= MMLR(1,ISTATE)
      IF(MMLR(1,ISTATE).LE.0) NCN= MMLR(2,ISTATE)
      IF(PSEL(ISTATE).LE.1) NCN= 999
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Call subroutine ALF that will locate the needed vibrational levels 
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      CALL ALF(NDATPT(ISTATE),RH(ISTATE),NCN,RR,V,SWF,VLIM(ISTATE),
     1  MAXMIN(ISTATE),VIN,KVtop,NVIBMX,VMAXX,AFLAG,ZMASS(3,IISTP),
     2  EPS(ISTATE),GV,INNODE,INNR,IWR)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** If a serious error occured during within ALF, then print out a 
c   warning, record the constants that we have and hope the program
c   doesn't call on the constants for which ALF could not calculate.
c
      IF (AFLAG.LT.0) THEN
          WRITE(6,600) ISTATE,IISTP
          IF(AFLAG.EQ.-1) THEN
              WRITE(6,601) KVtop,VIN, AFLAG
              STOP              !! ??  need to stop or just Print Warning
              ENDIF
          IF (AFLAG.EQ.-2) WRITE(6,602)
          IF (AFLAG.EQ.-3) WRITE(6,603)
          IF (AFLAG.EQ.-4) WRITE(6,604)
          IF (AFLAG.EQ.-5) WRITE(6,606)
          IF (AFLAG.EQ.-6) WRITE(6,608)
          IF (AFLAG.EQ.-8) THEN
              WRITE(6,610)
        ELSE
          WRITE(6,612) ISTATE,VIN
        ENDIF
      ENDIF
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Now to calculate the rotational constants for each vibrational level
c   of this state for this isotopologue
      IF(IOMEG(ISTATE).GT.0) THEN
          IF(NwCFT(ISTATE).GT.0) THEN
              WRITE(7,615)  NAME(1),MN(1,IISTP),NAME(2),MN(2,IISTP),
     1                      IOMEG(ISTATE), IOMEG(ISTATE)*IOMEG(ISTATE)
            ELSE
              WRITE(7,617)  NAME(1),MN(1,IISTP),NAME(2),MN(2,IISTP),
     1                      IOMEG(ISTATE), IOMEG(ISTATE)*IOMEG(ISTATE)
            ENDIF
        ELSEIF(IOMEG(ISTATE).LE.-2) THEN
          IF(NwCFT(ISTATE).GT.0) THEN
              WRITE(7,6615)  NAME(1),MN(1,IISTP),NAME(2),MN(2,IISTP),
     1                      IOMEG(ISTATE), -IOMEG(ISTATE)
            ELSE
              WRITE(7,6617)  NAME(1),MN(1,IISTP),NAME(2),MN(2,IISTP),
     1                      IOMEG(ISTATE), -IOMEG(ISTATE)
            ENDIF
        ELSE 
          IF(NwCFT(ISTATE).LT.0) THEN
              WRITE(7,614)  NAME(1),MN(1,IISTP),NAME(2),MN(2,IISTP)
            ELSE
              WRITE(7,618)  NAME(1),MN(1,IISTP),NAME(2),MN(2,IISTP)
            ENDIF
        ENDIF 
      DO  KV= 0,VIN
          AFLAG= 0
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Call subroutine SCHRQ to calculate the wavefunction required by
c   CDJOEL to calculate the rotational constants.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          CALL SCHRQ(KV,AFLAG,GV(KV),FWHM,PMAX,VLIM(ISTATE),V,SWF,BFCT,
     1             EPS(ISTATE),RMIN(ISTATE),RH(ISTATE),NDATPT(ISTATE),
     2                             NBEG,NEND,INNODE,INNR(KV),IWR,LPRWF)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Call subroutine CDJOEL to determine the rotational constants.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          CALL CDJOEL(GV(KV),NBEG,NEND,BvWN,RH(ISTATE),WARN,V,SWF,RM2,
     1                                                          RCNST)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Store molecular constants in array   ZK(v,J,isotope,state)
c
          ZK(KV,0,IISTP,ISTATE)= GV(KV)
          DO  J= 1,NROTMX
              ZK(KV,J,IISTP,ISTATE)= RCNST(J)
              ENDDO
          IF(NwCFT(ISTATE).LT.0) THEN
c*** Write band constants for 'normal' (Lambda= 0) case
              WRITE(7,616) KV,GV(KV),(RCNST(J),J=1,NROTMX)
            ELSE
c** For Lambda- or 2-Sigma doubling, calculate the q{B} parameter and
c   then print it with the other band constants
              qDBL = 0.d0
              DO  I= NBEG,NEND
                  SWF2= SWF(I)**2
                  qDBL= qDBL + SWF2*wRAD(I,ISTATE)
                  ENDDO
              qDBL= qDBL*qFCT
              WRITE(7,616) KV,GV(KV),(RCNST(J),J=1,NROTMX),qDBL
            ENDIF
          ENDDO
      flush(7)
      J=1
c
c** If all is well, then (without further ado) continue with the
c   calculations.
c
      RETURN
c-----------------------------------------------------------------------
  600 FORMAT(/'  *** INITDD ERROR ***',/4X,'For state',I3,
     1'  of isotope',I3,'  a serious error has occured:')
  601 FORMAT(' *** WARNING !! ALF finds highest level   v=',i3,'  is bel
     1ow desired (v=',I3,', J=',I3,')')
  602 FORMAT(4X,'The Schrodinger Solver was unable to use the initial tr
     1ial energy.')
  603 FORMAT(4X,'The Schrodinger Solver was unable to use the calculated
     1 trial energy')
  604 FORMAT(4X,'The Automatic Level Finder could not find the first vib
     1rational level.')
  606 FORMAT(4X,'The next calculated trial energy was too low.')
  608 FORMAT(4X,'The next calculated trial energy was too high.')
  610 FORMAT(4X,'A second minimum exists in this potential')
  612 FORMAT(4X,'Could not find vibrational levels of state',i3,
     1  ' beyond (v=',I3,')')
  614 FORMAT(/' For  ',A2,'(',I3,') - ',A2,'(',I3,')'/1x,11('--')/
     1  '   v       ','E',12x,'Bv',11x,'-Dv',13x,'Hv',13x,'Lv',
     2  12x,'Mv',13x,'Nv',13x,'Ov'/1x,58('=='))
  615 FORMAT(/' For  ',A2,'(',I3,') - ',A2,'(',I3,')'/1x,11('==')/
     1 ' Although   IOMEGA=',I2,', these band constants were obtained fo
     2r  [J(J+1) ',SP,I2,'] = 0'/1x,39('--')/'   v       ','E',12x,'Bv',
     3 11x,'-Dv',13x,'Hv',13x,'Lv',12x,'Mv',13x,'Nv',13x,'Ov'/
     4 1x,58('=='))
 6615 FORMAT(/' For  ',A2,'(',I3,') - ',A2,'(',I3,')'/1x,11('==')/
     1 ' Since   IOMEGA=',I3,', these band constants were obtained for',
     2'  [J(J+1) ',SP,I2,'] = 2'/1x,39('--')/'   v       ','E',12x,'Bv',
     3 11x,'-Dv',13x,'Hv',13x,'Lv',12x,'Mv',13x,'Nv',13x,'Ov'/
     4 1x,58('=='))
  616 FORMAT(I4,f12.4,f14.10,7(1PD15.7))
  617 FORMAT(/' For  ',A2,'(',I3,') - ',A2,'(',I3,')'/1x,11('==')/
     1 ' Although   IOMEGA=',I2,', these band constants were obtained fo
     2r  [J(J+1) ',SP,I2,'] = 0'/1x,39('--')/'   v       ','E',12x,'Bv',
     3  11x,'-Dv',13x,'Hv',13x,'Lv',13x,'Mv',13x,'Nv',13x,'Ov'13x,
     4  'qB(v)'/1x,67('=='))
 6617 FORMAT(/' For  ',A2,'(',I3,') - ',A2,'(',I3,')'/1x,11('==')/
     1 ' Since   IOMEGA=',I3,', these band constants were obtained for',
     2'  [J(J+1) ',SP,I2,'] = 2'/1x,39('--')/'   v       ','E',12x,'Bv',
     3  11x,'-Dv',13x,'Hv',13x,'Lv',13x,'Mv',13x,'Nv',13x,'Ov'13x,
     4  'qB(v)'/1x,67('=='))
  618 FORMAT(/' For  ',A2,'(',I3,') - ',A2,'(',I3,')'/1x,11('--')/
     1  '   v       ','E',12x,'Bv',11x,'-Dv',13x,'Hv',13x,'Lv',
     2  13x,'Mv',13x,'Nv',13x,'Ov'13x,'qB(v)'/1x,67('=='))
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE DEDP(IDAT,ISTATE,IISTP,ZMU,KVLEV,JROT,efPARITY,EO,
     1                                        VMAXX,FWHM,DEDPK,fcount)
c***********************************************************************
c** This subroutine calculates  dE/dp from the expectation values of the
c   partial derivatives of an analytical potential with respect to its
c   parameters {p(k)}:    dE/dp(k) = <vj|dV/dp(k)|vj>  stored in DEDPK.
c
c** On entry:
c    ISTATE  is the molecular state being considered.
c    IISTP   is the isotopologue being considered.
c    ZMU     is the reduced mass of the diatom in atomic units.
c    KVLEV   is the vibrational quantum number.
c    JROT    is the rotational quantum number.
c    EO      is the initial trial energy (DE in cm-1).
c    VMAXX  = is  MAX{barrier maximum} or VLIM for this state
c    ZK (in BLKISOT) is matrix of band constants for all levels of all ISOT
c
c** On exit:
c    EO       is the final calculated energy (DE in cm-1).
c    DEDPK(i) are the values of the partial derivative dE/dP(k).
c
c** Flags: Use only when debugging.
c    INNER specifies wave function matching (& initiation) conditions.
c          = 0 : Match inward & outward solutions at outermost wave
c                function maximum
c          <>0 : Match at inner edge of classically allowed region.
c          < 0 : uses zero slope inner boundary condition.
c          For most normal cases set INNER = 0,  but ......
c            To find "inner-well-dominated" solutions of an asymmetric
c            double minimum potential, set  INNER > 0.
c            To find symmetric eigenfunctions of a symmetric potential,
c            set INNER < 0  & start integration (set RMIN) at potential
c            mid point.
c    IWR   specifies the level of printing inside SCHRQ
c          <> 0 : print error & warning descriptions.
c          >= 1 : also print final eigenvalues & node count.
c          >= 2 : also show end-of-range wave function amplitudes.
c          >= 3 : print also intermediate trial eigenvalues, etc.
c    LPRWF specifies option of printing out generated wavefunction
c          > 0 : print wave function every LPRWF-th  point.
c          < 0 : compactly write to channel-7 every |LPRWF|-th wave
c                function value.
c          The first line identifies the level, gives the position of
c          1-st point and radial mesh, & states No. of  points.
c    fcount  counts the number of failed attempts to get desired level
c=======================================================================
      INCLUDE 'arrsizes.h'
      INCLUDE 'BLKISOT.h'
      INCLUDE 'BLKDATA.h'
      INCLUDE 'BLKPOT.h'
      INCLUDE 'BLKPARAM.h'
      INCLUDE 'BLKDVDP.h'
      INCLUDE 'BLKBOB.h'
      INCLUDE 'BLKBOBRF.h'
      INCLUDE 'BLKCOUNT.h'
c=======================================================================
      INTEGER IDAT, bandN, efPARITY, JRe, fcount,Idble,INNsav(NDATAMX)
c
      INTEGER I,ICOR,ISTATE,IISTP,INNER,J,JROT,JIN,KVLEV,KVtop,KV, 
     1  NBEG,NEND,INNODE, IWR, LPRWF, NCN, AFLAG, INNR(0:NVIBMX)
      REAL*8 ZMU,EO,BFCT,JFCT,JFCTP,JFCTL,FWHM,UMAX, RM2, ETRY,VMAXX,
     1  DGDV2,SWF2,JFCTA,JFCTB,DUARe,DUBRe,DTARe,DTBRe,DLDRe,DEROT,
     2  DEROTB,JFCTDBL,JFCTD,muFCT,C3gu,T3,RHSQ,  GV(0:NVIBMX),
     3  Vtotal(NPNTMX),SWF(NPNTMX), V(NPNTMX), DEDPK(HPARMX),
     4  RR(NPNTMX),EvjSAV(NDATAMX)
c
      COMMON /VBLIK/Vtotal
c
c** Vibrational Band Constants for generating trial energies.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c** Initializing the arrays and variables
      DATA LPRWF/0/,INNODE/1/,Idble/1/
c
c** To calculate the values for <vj|dV/dP(k)|vj>
c   we must first determine the wave equation for the vj state.
      IWR= 0                              !! normal setting
!!    IWR= 1                              !! when testing for problems
!!    IWR= 2                              !! when testing for problems
!!    IWR= 3                              !! when testing for problems
      IF(Idble.EQ.1) THEN
          INNsav= -9
          Idble= Idble+ 1
          ENDIF
      muFCT= 16.857629206d0/ZMU
      RHSQ= RH(ISTATE)*RH(ISTATE)
      BFCT= RHSQ/muFCT
      JFCT= DBLE(JROT*(JROT+1))
      JFCTL= JFCT-IOMEG(ISTATE)**2
      DO  J= 1,HPARMX
          DEDPK(J)= 0.d0
          ENDDO
      IF(IOMEG(ISTATE).GT.0) JFCT= JFCT - 
     1                               DBLE(IOMEG(ISTATE)*IOMEG(ISTATE))
c** If using band constants for this state ....
      IF(PSEL(ISTATE).EQ.-1) THEN
          IF(NBC(KVLEV,IISTP,ISTATE).GT.0) THEN
              I= BCPARI(KVLEV,IISTP,ISTATE)
              DEDPK(I)= 1.d0
              EO= ZBC(KVLEV,0,IISTP,ISTATE)
              JFCTP= JFCTL
              IF(NBC(KVLEV,IISTP,ISTATE).GT.1) THEN
                  DO  J= 2, NBC(KVLEV,IISTP,ISTATE)
                      I= I+ 1
                      DEDPK(I)= JFCTP
                      EO= EO+ ZBC(KVLEV,J-1,IISTP,ISTATE)*JFCTP
                      JFCTP= JFCTP*JFCTL 
                      ENDDO
                  IF(NQC(KVLEV,IISTP,ISTATE).GT.0) THEN
                      JFCTP= JFCT*0.5d0*(efPARITY-efREF(ISTATE))
                      IF(IOMEG(ISTATE).EQ.-1) THEN
                          JFCTL= JFCT
                          IF(efPARITY.GT.0) JFCTP= 0.5d0*JROT
                          IF(efPARITY.EQ.0) JFCTP=  0.d0
                          IF(efPARITY.LT.0) JFCTP= -0.5d0*(JROT+1)
                          ENDIF
                      DO  J= 1, NQC(KVLEV,IISTP,ISTATE)
                          I= I+ 1
                          DEDPK(I)= JFCTP
                          EO= EO+ ZQC(KVLEV,J-1,IISTP,ISTATE)*JFCTP
                          JFCTP= JFCTP*JFCTL 
                          ENDDO
                      ENDIF 
                  ENDIF
              ENDIF
          RETURN
          ENDIF
c** Calculating the trial energy value from the band constants
      ETRY= ZK(KVLEV,0,IISTP,ISTATE)
      IF(JROT.GT.0) THEN
          DEROT= 9.d9
          DO I=1,NROTMX
              DEROTB= DEROT
              JFCTP= JFCT**I
              IF(IOMEG(ISTATE).EQ.-1) JFCTP= JFCTP - 2**I  !! ??? huh ??
              DEROT= ZK(KVLEV,I,IISTP,ISTATE) * JFCTP
c... if centrifugal term bigger than previous one - truncate summation
              ETRY= ETRY + DEROT
              ENDDO
    4     ENDIF
      IF(IOMEG(ISTATE).GT.0) THEN 
c** For Lambda doubling, prepare rotational/mass factors
          JFCTDBL= 0.5d0*(efPARITY-efREF(ISTATE))  
     1               * (DBLE(JROT*(JROT+1)) * muFCT**2)**IOMEG(ISTATE)
          ENDIF
      IF(IOMEG(ISTATE).EQ.-1) THEN
c** For doublet Sigma splitting, prepare rotational/mass factors
          IF(efPARITY.GT.0) JFCTDBL=  0.5d0*JROT*muFCT
          IF(efPARITY.EQ.0) JFCTDBL=  0.d0
          IF(efPARITY.LT.0) JFCTDBL= -0.5d0*(JROT+1)*muFCT
          ENDIF
c
c** Generating potential function including centrifugal, adiabatic BOB,
c   and rotational non-adiabatic BOB, and if appropriate, Lambda 
c   doubling or 2\Sigma doubling radial functions.
cc
c==> 'Model-B' <== for Li2(A) *** skip smartass centrifugal stuff and put
c              centrifugal BOB term 2*B(r) in PEF directly
      JFCT= JFCT*RHSQ
      JFCTD= JFCTDBL*RHSQ/muFCT
      bandN= IB(IDAT)
      C3gu= (2.d0/3.d0)*CmVAL(2,ISTATE)  !! for special case {6,7}Li2(A)
      DO I= 1,NDATPT(ISTATE)
          RM2= 1/RD(I,ISTATE)**2
          V(I)= BFCT*(VPOT(I,ISTATE) + ZMUA(IISTP,ISTATE)*UAR(I,ISTATE)
     1                              + ZMUB(IISTP,ISTATE)*UBR(I,ISTATE))
     2          + JFCT*(1.0d0 + ZMTA(IISTP,ISTATE)*TAR(I,ISTATE)
     3                         + ZMTB(IISTP,ISTATE)*TBR(I,ISTATE))*RM2
c** Special BOB correction for A-state Li2 and analogous cases. !!!!!!!!
          IF(IOMEG(ISTATE).EQ.-2) V(I)= V(I) + 2.d0*RHSQ*RM2
c** Add radial potential for Lambda- or 2-Sigma splitting
c ... note that power of 1/r**2  included in  wRAD  array ...
          IF((IOMEG(ISTATE).GT.0).OR.(IOMEG(ISTATE).EQ.-1)) 
     1                               V(I)= V(I) + JFCTD*wRAD(I,ISTATE)
          Vtotal(I)= V(I)/BFCT         !! in [cm-1] needed for TP search
          IF((MMLR(1,ISTATE).LE.0).AND.(AN(1).EQ.AN(2))
     1                         .AND.(MN(1,IISTP).NE.MN(2,IISTP))) THEN
cc* Add g/u symmetry breakdown correction for special case {6,7}Li2(A) !!
              T3= C3gu/RD(I,ISTATE)**3
              Vtotal(I)= Vtotal(I)+ T3- DSQRT(T3**2+ 3.085959756d-2)
              V(I)= Vtotal(I)*BFCT
              ENDIF
          ENDDO
      ICOR= 0
      INNER= 0
      EO= ETRY
   10 KV= KVLEV
      IF(INNsav(IDAT).GE.0) THEN
          EO= EvjSAV(IDAT)
          INNER= INNsav(IDAT)
          ENDIF
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      CALL SCHRQ(KV,JROT,EO,FWHM,UMAX,VLIM(ISTATE),V,SWF,BFCT,
     1             EPS(ISTATE),RMIN(ISTATE),RH(ISTATE),NDATPT(ISTATE),
     2                               NBEG,NEND,INNODE,INNER,IWR,LPRWF)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IF(KV.NE.KVLEV) THEN
c** If SCHRQ found the wrong level - usually only in double-well cases
          IF(IWR.NE.0) write(6,666) KVLEV,JROT,ETRY,KV,EO
          ICOR= ICOR+1
c???  new expt ?????????????????????????????????????????
          IF(ICOR.EQ.1) THEN
              NCN= 99
              AFLAG= JROT
              DO  I=1,NDATPT(ISTATE)
                  RR(I)= RD(I,ISTATE)    !! 1D distance array for ALF
                  ENDDO
c** On first failure, perform brute force ALF search for missing level
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
              CALL ALF(NDATPT(ISTATE),RH(ISTATE),NCN,RR,V,SWF,
     1           VLIM(ISTATE),MAXMIN(ISTATE),KVLEV,KVtop,NVIBMX,VMAXX,
     2           AFLAG,ZMASS(3,IISTP),EPS(ISTATE),GV,INNODE,INNR,IWR)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
              IF(AFLAG.GE.0) THEN
                  EO= GV(KVLEV)
                  INNER= INNR(KVLEV)
                  GO TO 20
                ELSE                  !! if cannot find desired level
                  EO= -9.9d9
                  IWR= 0                          !! normal setting
                  RETURN
                ENDIF
              ENDIF
c???  new expt ?????????????????????????????????????????
          IF((ICOR.LE.10).AND.(KV.GE.0)) THEN
c... SCECOR uses semiclassical methods to estimate correct energy
            CALL SCECOR(KV,KVLEV,JROT,INNER,ICOR,IWR,EO,RH(ISTATE),BFCT,
     1       NDATPT(ISTATE),MMLR(1,ISTATE),V,VMAXX,VLIM(ISTATE),DGDV2)
c***********************************************************************
              V= KVLEV
              GOTO 10
              ENDIF 
c** If the calculated eigenvalue is still for the wrong vibrational
c   level, then write out a warning and skip the calculation of the
c   partial derivatives (hence setting them to zero).
          fcount= fcount+1
          WRITE(6,610) fcount,KVLEV,JROT,KV
!!        IF(fcount.ge.100) THEN
          IF(fcount.ge.100000) THEN
              WRITE(6,612)
  612 FORMAT(/' *** Excessive SCECOR failures, so stop and figure out wh
     1y *** ' //)
              STOP
              ENDIF        
c.. eigenvalue of -9.9d9 indicates that eigenvalue search failed completely
          EO= -9.9d9
          IWR= 0                          !! normal setting
          RETURN
          ENDIF
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** If the calculated wavefunction is for the right vibrational level,
c   then continue with calculation of the partial derivatives by
c   integration using the modified trapezoidal rule.
c
c** First determine which partials need to be changed, increment so
c   that the (fluorescence term values and) lower states are skipped.
c
   20 JRe= POTPARI(ISTATE)+ 1
      IF(MAXMIN(ISTATE).EQ.2) THEN
!!        WRITE(24,622) IDAT,KVLEV,JROT,EO,INNER
          EvjSAV(IDAT)= EO         !! Save trial energy for next cycle!!
          INNsav(IDAT)= INNER
          ENDIF
!!622 format('  IDAT=',I5,'   KVLEV=',I3,'   J=',I4,'   E=',
!!   1 f10.3,' INNER=',I3)
      IF(PSEL(ISTATE).EQ.6) JRe= POTPARI(ISTATE)
c** Now calculate the rotational factors JFCT, JFCTA & JFCTB
      JFCT= JFCT/BFCT
      JFCTA= JFCT * ZMTA(IISTP,ISTATE)
      JFCTB= JFCT * ZMTB(IISTP,ISTATE)
      DUARe= 0.d0
      DUBRe= 0.d0
      DTARe= 0.d0
      DTBRe= 0.d0
      DLDRe= 0.d0
c** Eigenvalue derivative calc using compact partial derivative array
      IF(PSEL(ISTATE).LE.0) RETURN
      DO  I= NBEG,NEND
          SWF2= SWF(I)**2
          IF((I.EQ.NBEG).OR.(I.EQ.NEND)) SWF2= 0.5d0*SWF2
c ... collect contributions of BOB terms to derivatives w.r.t. Re
          IF(NUA(ISTATE).GE.0) DUARe= DUARe+ SWF2*DUADRe(I,ISTATE)
          IF(NUB(ISTATE).GE.0) DUBRe= DUBRe+ SWF2*DUBDRe(I,ISTATE)
          IF(NTA(ISTATE).GE.0) DTARe= DTARe+ SWF2*DTADRe(I,ISTATE)
          IF(NTB(ISTATE).GE.0) DTBRe= DTBRe+ SWF2*DTBDRe(I,ISTATE)
          IF(NwCFT(ISTATE).GE.0) DLDRe= DLDRe+ SWF2*DLDDRe(I,ISTATE)
          DO  J= POTPARI(ISTATE), POTPARF(ISTATE)
              DEDPK(J)= DEDPK(J) + SWF2*DVtot(J,I)
              ENDDO
          ENDDO
      DEDPK(JRe)= DEDPK(JRe) + DUARe*ZMUA(IISTP,ISTATE)
     1      + DUBRe*ZMUB(IISTP,ISTATE) + JFCTA*DTARe + JFCTB*DTBRe
     2      + JFCTDBL*DLDRE
      DO  J= POTPARI(ISTATE), POTPARF(ISTATE)
          DEDPK(J)= DEDPK(J)*RH(ISTATE)
          ENDDO
      IF(NUA(ISTATE).GE.0) THEN
          DO  I= NBEG,NEND
              SWF2= SWF(I)**2
              IF((I.EQ.NBEG).OR.(I.EQ.NEND)) SWF2= 0.5d0*SWF2
              DO  J= UAPARI(ISTATE), UAPARF(ISTATE)
                  DEDPK(J)= DEDPK(J) + SWF2*DVtot(J,I)
                  ENDDO
              ENDDO
          DO  J= UAPARI(ISTATE), UAPARF(ISTATE)
              DEDPK(J)= DEDPK(J)* ZMUA(IISTP,ISTATE)*RH(ISTATE)
              ENDDO
          ENDIF 
      IF(NUB(ISTATE).GE.0) THEN
          DO  I= NBEG,NEND
              SWF2= SWF(I)**2
              IF((I.EQ.NBEG).OR.(I.EQ.NEND)) SWF2= 0.5d0*SWF2
              DO  J= UBPARI(ISTATE), UBPARF(ISTATE)
                  DEDPK(J)= DEDPK(J) + SWF2*DVtot(J,I)
                  ENDDO
              ENDDO
          DO  J= UBPARI(ISTATE), UBPARF(ISTATE)
              DEDPK(J)= DEDPK(J)* ZMUB(IISTP,ISTATE)* RH(ISTATE)
              ENDDO
          ENDIF 
      IF(NTA(ISTATE).GE.0) THEN
          DO  I= NBEG,NEND
              SWF2= SWF(I)**2
              IF((I.EQ.NBEG).OR.(I.EQ.NEND)) SWF2= 0.5d0*SWF2
              DO  J= TAPARI(ISTATE), TAPARF(ISTATE)
                  DEDPK(J)= DEDPK(J) + SWF2*DVtot(J,I)
                  ENDDO
              ENDDO
          DO  J= TAPARI(ISTATE), TAPARF(ISTATE)
              DEDPK(J)= DEDPK(J)* JFCTA* RH(ISTATE)
              ENDDO
          ENDIF 
      IF(NTB(ISTATE).GE.0) THEN
          DO  I= NBEG,NEND
              SWF2= SWF(I)**2
              IF((I.EQ.NBEG).OR.(I.EQ.NEND)) SWF2= 0.5d0*SWF2
              DO  J= TBPARI(ISTATE), TBPARF(ISTATE)
                  DEDPK(J)= DEDPK(J) + SWF2*DVtot(J,I)
                  ENDDO
              ENDDO
          DO  J= TBPARI(ISTATE), TBPARF(ISTATE)
              DEDPK(J)= DEDPK(J)* JFCTB* RH(ISTATE)
              ENDDO
          ENDIF 
      IF(NwCFT(ISTATE).GE.0) THEN
          DO  I= NBEG,NEND
              SWF2= SWF(I)**2
              IF((I.EQ.NBEG).OR.(I.EQ.NEND)) SWF2= 0.5d0*SWF2
              DO  J= LDPARI(ISTATE), LDPARF(ISTATE)
                  DEDPK(J)= DEDPK(J) + SWF2*DVtot(J,I)
                  ENDDO
              ENDDO
          DO  J= LDPARI(ISTATE), LDPARF(ISTATE)
              DEDPK(J)= DEDPK(J)* JFCTDBL* RH(ISTATE)
              ENDDO
          ENDIF
      RETURN
c-----------------------------------------------------------------------
  610 FORMAT(' *** SCECOR failed',I4,' times,   Currently Seeking   v=',
     1   i3,', J=',i3,';  Found  v=',I3)
  666 FORMAT(' Search for v=',I3,'   J=',I3,'  starting from  E=',
     1  f9.2,' finds  E(v=',I3,')=', f9.2)
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE DWDP(IDAT,ISTATE,EXPT,ZMU,vb,Jr,EO,width,DEDPK,dWdPk)
c=======================================================================
c** This subroutine calculates dW/dP using Pajunen's quadrature method 
c   Ref: JCP. 71(6), 2618, 1979.
c** The values of dW/dP(k) for each parameter P(k) are stored in dWdPk.
c====================RJL's Version of 11 May 2002=======================
c** On entry:
c    IDAT    is the experimental data number
c    ISTATE  is the molecular state being considered.
c    IISTP   is the isotopologue being considered.
c    ZMU     is the reduced mass of the diatom in atomic units.
c    vb      is the vibrational quantum number.
c    Jr      is the rotational quantum number.
c    EO      is the energy for this level (DE in cm-1).
c    width   is the FWHM level width calculated in SCHRQ
c    DEDPK(i) are the values of the partial derivative dE/dP(k)
c    
c** On exit:
c    width   is the width calculated herein (DE in cm-1).
c    dWdPk(i) are the values of the partial derivative dW/dP(k).
c=======================================================================
      INCLUDE 'arrsizes.h'
      INCLUDE 'BLKCOUNT.h'
      INTEGER I,ISTATE,IDAT,NN,NST,vb,Jr,IPV
c ** Phase integral and partial derivative arrays ....
      REAL*8 Iwell(0:HPARMX,-1:1),Ibarr(0:HPARMX,-1:1), DEDPK(HPARMX),
     1       dWdPk(HPARMX)
      REAL*8 ARG,dARG,ZMU,EO,FCTOR,FWHM,Pi,EXPT,width, R1, R2, R3,
     1  KAPPA,KAPPAcl,EMSC,EPSRJ,dKdEPS, XX,TI,COR,dvdG2pi
      DATA pi/3.141592653589793D0/
*----------------------------------------------------------------------*
*** Begin by calling subroutine locateTP to find the three turning 
*   points of the quasi-bound level 
      FCTOR= DSQRT(ZMU/16.857629206d0)
      DO IPV= POTPARI(ISTATE),HPARF(ISTATE)
          dWdPk(IPV)= 0.0d0
          ENDDO 
      CALL locateTP(IDAT,ISTATE,vb,Jr,EO,R1,R2,R3)
c *** if some turning points are not found, ignore this datum
      IF((R1.LT.0.d0) .OR. (R2.LT.0.d0) .OR. (R3.LT.0.d0)) THEN
          width= EXPT
          WRITE(6,600) IDAT
          RETURN
  600 FORMAT(' <<Energy of  WIDTH ',I5,' is above barrier maximum, so om
     &it this datum!>>')
          ENDIF
c** Call subroutine  phaseIntegral  to calculate the phase integrals 
c   using Pajunen's quadrature method
      CALL PhaseIntegral(IDAT,ISTATE,R2,R3,-1,EO,DEDPK,Ibarr)
      CALL PhaseIntegral(IDAT,ISTATE,R1,R2,0,EO,DEDPK,Iwell)
c ... save WIDTH calculated in SCHRQ to test vs. new value
      FWHM= width
*** Calculate the width
      EPSRJ= 2.0d0*FCTOR*Ibarr(0,-1)
      KAPPAcl= DEXP(-EPSRJ)
      KAPPA= DSQRT(1.d0 + KAPPAcl) - 1.d0
c ... alternate calculation to give better precision for small TUN0
      IF(KAPPAcl.LT.1.d-5) KAPPA= KAPPAcl*(0.5d0- KAPPAcl*(0.125d0-
     1                                                0.0625d0*KAPPAcl))
      KAPPA= 4.d0* KAPPA /(KAPPA+ 2.d0)
c** Derivative of complex gamma function argument calculated as
c .....  EPSRJ= -2.* PI* EMSC
c  per eq.(6.1.27) in Abramowitz and Stegun.
      EMSC= - EPSRJ/(2.d0*pi)
      NST= DABS(EMSC)*1.d2
      NST= MAX0(NST,4)
      ARG= -1.963510026021423d0
      dARG= 0.d0
      DO  I= 0,NST
          NN= I
          XX= I + 0.5d0
          TI= 1.d0/((XX/EMSC)**2 + 1.d0)
          dARG= dARG+ XX*TI*TI
          TI= TI/XX
          ARG= ARG+TI
          IF(DABS(TI).LT.1.d-10) GO TO 233
          ENDDO
c ... and use integral approximation for tails of summations ...
  233 COR= 0.5d0*(EMSC/(NN + 1.d0))**2
      ARG= ARG + COR - COR**2
      dARG= dARG + EMSC**2 *(COR - 2.d0* COR**2)
      dvdG2pi= FCTOR * (Iwell(0,0) 
     1           + (DLOG(DABS(EMSC)) - ARG)* Ibarr(0,0)/(2.d0*pi) )
      width= KAPPA/dvdG2pi
c??????
c      WRITE(32,320) vb,Jr,EO,width, width/FWHM- 1.d0, FWHM
c 320  FORMAT ('Level v=',I3,' J=',I3,'  Ep=',G16.8,'    FWHM(new)=',
c    1  1Pd16.8/ 12x,'Relative Diff.=',d9.2,'    FWHM(SCHRQ)=', d16.8)
c??????
      dKdEPS= DSQRT(1.d0 + KAPPAcl)
      dKdEPS= -4.d0*KAPPAcl/(dKdEPS*(dKdEPS+ 1.d0)**2)
      DO  IPV= POTPARI(ISTATE),HPARF(ISTATE)
          dWdPk(IPV)= dKdEPS * FCTOR * Ibarr(IPV,0)/dvdG2pi
     1             - (KAPPA/dvdG2pi**2)* (0.5d0 * FCTOR * Iwell(IPV,1)
     2         + (FCTOR**2 *Ibarr(0,0)/(2.d0*pi* EPSRJ))* Ibarr(IPV,0) 
     3                            * (1.d0 - dARG* 8.d0*pi**2/EPSRJ**2)
     4    - (FCTOR* Ibarr(IPV,1)/(4.d0*pi))* (DLOG(DABS(EMSC)- ARG)) )
          ENDDO
c?????
cc    write(32,321)  (dWdPk(IPV), IPV= POTPARI(ISTATE),HPARF(ISTATE)
cc    write(33,321)  (dEdPk(IPV), IPV= POTPARI(ISTATE),HPARF(ISTATE)
cc    write(33,*)
cc321 format(1Pd11.3,6d11.3/(11x,6d11.3))
c??????
      RETURN
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE locateTP(IDAT,ISTATE,vb,Jr,EO,R1,R2,R3)
c=======================================================================
c  Subroutine locateTP locates the turning points R1, R2 and R3 for    
c  level  v= vb, J= Jr at energy  EO.  The search for each starts with 
c  a scan over the pre-prepared total potential function array  Vtotal, 
c  and then applies iterative local linear interpolation until the 
c  numerical-noise machine precision limit is reached.
c======================================================================= 
c** On entry
c    IDAT    is the experimental data number
c    ISTATE  is the molecular state being considered.
c    vb      is the vibrational quantum number.
c    Jr      is the rotational quantum number.
c    EO      is the energy for this level (DE in cm-1).
c** On exit
c  R1          is the array of inner turning points.                 
c  R2          is the array of second turning points.                
c  R3          is the array of third turning points.                 
c=======================================================================
      INCLUDE 'arrsizes.h'
      INCLUDE 'BLKISOT.h'
      INCLUDE 'BLKDATA.h'
      INCLUDE 'BLKPOT.h'
c-----------------------------------------------------------------------
c** Define types for local variables
      INTEGER  I,IDAT,ISTATE,vb,Jr,index
      REAL*8 VLIMT,R1,R2,R3,EO,bMi,Vmid,EMV,EMVB,EMVBB,RTB,RTBB,
     1   RR(1),RM2(1),VV(1),Vtotal(NPNTMX)
      COMMON /VBLIK/Vtotal
c*** Start search for innermost turning point ...
      index= 1
      EMV= -99.d0
c ... First, scan to find first mesh point past innermost turning point
      DO  I= 1, NDATPT(ISTATE)
          EMVB= EMV
          EMV= EO- Vtotal(I)
          IF(EMV.GT.0.d0) THEN
              INDEX= I
              GOTO 4 
              ENDIF
          ENDDO
      WRITE(6,600) 'R1',vb,Jr,EO
      R1= -1.0d0
      GO TO 999
    4 R1= RD(index,ISTATE)
c ... test to see if by accident, mesh point IS exactly  R1
      IF(DABS(EO-Vtotal(INDEX)).GT.0.d0) THEN
          R1= RD(index,ISTATE)
          RTB= RD(index-1,ISTATE)
  10      RTBB= RTB
          EMVBB= EMVB
          RTB= R1
          EMVB= EMV
c ... interpolate linearly between the two most recent estimates
          R1= RTB - (RTB-RTBB)*EMVB/(EMVB-EMVBB)      
          IF(PSEL(ISTATE).GT.0) THEN
c ... obtain potential function value at new turning point estimate
              CALL VGEN(ISTATE,R1,Vmid,bMi,IDAT)
            ELSE
              RR(1)= R1
              RM2(1)= 1.0d0/R1**2
c ... for fixed pointwise potential, interpolate for potential value
              CALL PREPOTT(0,AN(1),AN(2),MN(1,1),MN(2,1),1,VLIMT,RR,VV)
              Vmid= VV(1) 
            ENDIF
          EMV= EO- VMID
ccccc
ccc       WRITE (30,21) RTB,'R1',R1,RTBB,EO,Vmid
ccc21 FORMAT (' RTB=',F11.8,2x,A2,'=',F11.8,'  RTBB=',F11.8,'  E=',
ccc  1  1PD16.9,' Vmid=',D16.9)
ccccc
c ... test for convergence to machine precision 
          IF(DABS(EMV).LT.DABS(EMVB)) GOTO 10
ccc
cc        WRITE(30,100) vb,Jr,'R1',R1,EO,EMV,EMVB
cc100 FORMAT('For  v=',I3,'  J=',I3,2x,a2,'=',F9.6,'  E =',1PD15.8,
cc   &  '  EmV=',D15.8,'  EmVB=',D15.8)
ccc
          R1= RTB
          ENDIF
c
c****************now Locate second turning point R2 *****************
      EMV= +99.d0
c ... Now, scan to find first mesh point past the second turning point
      DO  I= INDEX, NDATPT(ISTATE)
          EMVB= EMV
          EMV= EO- Vtotal(I)
          IF(EMV.LT.0.d0) THEN
              INDEX= I
              GOTO 14
              ENDIF
          ENDDO
      WRITE(6,600) 'R2',vb,Jr,EO
      R2= -1.0d0
      GO TO 999
   14 R2= RD(index,ISTATE)
c ... test to see if by accident, mesh point IS exactly  R2
      IF(DABS(EO-Vtotal(INDEX)).GT.0.d0) THEN
          R2= RD(index,ISTATE)
          RTB= RD(index-1,ISTATE)
  20      RTBB= RTB
          EMVBB= EMVB
          RTB= R2
          EMVB= EMV
c ... interpolate linearly between the two most recent estimates
          R2= RTB - (RTB-RTBB)*EMVB/(EMVB-EMVBB)
          IF(PSEL(ISTATE).GT.0) THEN
c ... obtain PEC value at new turning point estimate: no deriv or BOB neeeded
              CALL VGEN(ISTATE,R2,Vmid,bMi,IDAT)
            ELSE
              RR(1)= R2
              RM2(1)= 1.0d0/R2**2
c ... for fixed pointwise potential, interpolate for potential value
              CALL PREPOTT(0,AN(1),AN(2),MN(1,1),MN(2,1),1,VLIMT,RR,VV)
              Vmid= VV(1)
            ENDIF
          EMV= EO- VMID
ccc  
ccc       WRITE (30,21) RTB,'R2',R2,RTBB,EO,Vmid
ccc  
c ... test for conergence to machine precision
          IF(DABS(EMV).LT.DABS(EMVB)) GOTO 20
ccc
ccc       WRITE(30,100) vb,Jr,'R2',R2,EO,EMV,EMVB
ccc
          R2= RTB
          ENDIF
c
c****************now Locate third turning point R3 *****************
      EMV= +99.d0
c ... Now, scan to find first mesh point past the third turning point
      DO  I= INDEX, NDATPT(ISTATE)
          EMVB= EMV
          EMV= EO- Vtotal(I)
          IF(EMV.GT.0.d0) THEN
              INDEX= I
              GOTO 24
              ENDIF
          ENDDO
      WRITE(6,600) 'R3',vb,Jr,EO
      R3= -1.0d0
      GO TO 999
   24 R3= RD(index,ISTATE)
c ... test to see if by accident, mesh point IS exactly  R3
      IF(DABS(EO-Vtotal(INDEX)).GT.0.d0) THEN
          R3= RD(index,ISTATE)
          RTB= RD(index-1,ISTATE)
  30      RTBB= RTB
          EMVBB= EMVB
          RTB= R3
          EMVB= EMV
c ... interpolate linearly between the two most recent estimates
          R3= RTB - (RTB-RTBB)*EMVB/(EMVB-EMVBB)
          IF(PSEL(ISTATE).GT.0) THEN
c ... obtain potential function value at new turning point estimate
              CALL VGEN(ISTATE,R3,Vmid,bMi,IDAT)
            ELSE
              RR(1)= R3
              RM2(1)= 1.0d0/R3**2
c ... for fixed pointwise potential, interpolate for potential value
              CALL PREPOTT(0,AN(1),AN(2),MN(1,1),MN(2,1),1,VLIMT,RR,VV)
              Vmid= VV(1)
            ENDIF
          EMV= EO- VMID
ccc  
ccc       WRITE (30,21) RTB,'R3',R3,RTBB,EO,Vmid
ccc  
c ... test for conergence to machine precision
          IF(DABS(EMV).LT.DABS(EMVB)) GOTO 30
ccc
ccc       WRITE(30,100) vb,Jr,'R3',R3,EO,EMV,EMVB
ccc
          R3= RTB
          ENDIF
  999 RETURN
  600 FORMAT(' Turning point  ',A2,'  at   E(v=',I3,' J=',I3,')=',
     1  G14.7,'  lies beyond  RMAX.')
      END 
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE PhaseIntegral(IDAT,ISTATE,Ri,Ro,k1,EO,DEDPK,Jntegral)
c=======================================================================
c  Subroutine PhaseIntegral calculates phase integrals of the integrand
c  f(R)/[E-V(R)]**(k+1/2), for k=-1,0 or 1  on interval between turning 
c  points  Ri and Ro, where f(R) is either 1 or a derivative of the
c  potential w.r.t. some potential parameter. 
c             Ro                                         N
c  Jntegral = Int dr f(r)/|E-V(r)|^{k+1/2} = ½(Ro - Ri) Sum{Wi*F(Zi)}
c             Ri                                        i=1
c  Evaluate integrals for  k=k1 to 1  where for k=k1 f(r)=1
c                                           for k>k1 f(r)=[dE/dp-dV/dp]
c====================RJL's Version of 11 May 2002=======================
c  Tests against proper Gaussian (25.4.38 & 25.4.40) gave agreement of
c    ca. 10^{-10} to 10^{-8) for I_{-1}^{barr}({1})  and  ranging from 
c    10^{-13} to [rare worst case] 10^{-6} fo I_0^{barr}({deriv})
c    typically 1-^{-6} for I_0^{well}({1}).
c* For k=1 integrals, comparing PP "k=1" vs. "k=3" integration gave
c  typical agreement for I_1^{barr}({deriv.}) of ca. 10^{-6}, but
c  occasionally as bad as 10^{-1} and as good as 10^{-9}
c   [?? check cgce. of turning point search!].
c** I_1^{well}({f})  derivatives always poor cgce. w.r.t. k=1 vs. k=3
c               PP method integration tests ... ???
c* Notice deriv.'s w.r.t. q(r) parameters typically orders of magnitude
c     poorer than others (??)
c
c* On entry:  IDAT    is the experimental datum number
c-----------  ISTATE  is the molecular state being considered.
c             Ri      is the inner turning point
c             Ro      is the outer turning point
c             EO      is the energy for this level (DE in cm-1).
c             k       is the powers in the phase integral Jntegral
c             n=0  for  f(r)= 1;  n=1  for  f(r)= [dV(r)/dp_k - dE/dp_k]
c             DEDPK(i) are the values of the partial derivative dE/dp_k.
c* On exit:  Jntegral(j,k)    are the phase integral(s)
c=======================================================================
      INCLUDE 'arrsizes.h'
      INCLUDE 'BLKCOUNT.h'
c** Type statements & common block for data
      INTEGER i,IPV,ISTATE, k,k1,kmx,M,IDAT
c** M  is the number of quadrature mesh points
      PARAMETER (M=21)
c** Common block for partial derivatives of potential at the one distance RDIST
c   and HPP derivatives for uncertainties
      REAL*8 dVdPk(HPARMX),dDe(0:NbetaMX),dDedRe
      COMMON /dVdPkBLK/dVdPk,dDe,dDedRe
c=======================================================================
      REAL*8 Ri,Ro,EO,Fzt,RDIST,VDIST,BETADIST,EMinusV,
     1  Jntegral(0:HPARMX,-1:1),DEDPK(HPARMX)
c
      REAL*8 Zp21(M),Wp21k1(M)
ccc  1   ,Wp21k3(M)
ccc  2   ,Zp11(M),Wp11k1(M),Wp11k3(M),Zp15(M),Wp15k1(M),Wp15k3(M)
ccc
ccc   DATA Zp11/0.9659258262890683, 0.8660254037844386,
ccc  1          0.7071067811865475, 0.5000000000000000,
ccc  2          0.2588190451025208, 
ccc  3          0.0000000000000000,
ccc  4         -0.2588190451025208,-0.5000000000000000,
ccc  5         -0.7071067811865475,-0.8660254037844386,
ccc  6         -0.9659258262890683,
ccc  7          0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0/
ccc
ccc   DATA Wp11k1/22.06633397656787,-37.69911184307752,
ccc  1            35.60471674068432,-37.69911184307752,
ccc  2            36.57672889044161,
ccc  3           -37.69911184307752,
ccc  4            36.57672889044161,-37.69911184307752,
ccc  5            35.60471674068432,-37.69911184307752,
ccc  6            22.06633397656787,
ccc  7          0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0/
ccc
ccc   DATA Wp11k3/228.3214545745810, -720.4719152232592,
ccc  1           1139.3509357018984,-1357.1680263507907,
ccc  2           1447.1946273399755,
ccc  3          -1474.4541520848097,
ccc  4           1447.1946273399755,-1357.1680263507907,
ccc  5           1139.3509357018984, -720.4719152232592,
ccc  6            228.3214545745810,
ccc  7          0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0/
ccc
ccc   DATA Zp15/0.9807852804032304, 0.9238795325112868, 
ccc  1          0.8314696123025452, 0.7071067811865475,
ccc  2          0.5555702330196022, 0.3826834323650898,
ccc  3          0.1950903220161283, 
ccc  4          0.0000000000000000,
ccc  5         -0.1950903220161283,-0.3826834323650898,
ccc  6         -0.5555702330196022,-0.7071067811865475,
ccc  7         -0.8314696123025452,-0.9238795325112868,
ccc  8         -0.9807852804032304,
ccc  9          0.0,0.0,0.0,0.0,0.0,0.0/
ccc
ccc   DATA Wp15k1/29.62981929591175,-50.26548245743668,
ccc  1            47.72092686124880,-50.26548245743664,
ccc  2            49.12943331558201,-50.26548245743658,
ccc  3            49.44900912828539,
ccc  4           -50.26548245743656,
ccc  5            49.44900912828539,-50.26548245743658,
ccc  6            49.12943331558201,-50.26548245743664,
ccc  7            47.72092686124880,-50.26548245743668,
ccc  8            29.62981929591175, 
ccc  9          0.0,0.0,0.0,0.0,0.0,0.0/
ccc
ccc   DATA Wp15k3/1164.639428963841,-3580.432803129281,
ccc  1            5525.620597073791,-6534.512719466765,
ccc  2            7020.542275282852,-7276.911407677031,
ccc  3            7400.700330802901,
ccc  4           -7439.291403700618,
ccc  5            7400.700330802901,-7276.911407677031,
ccc  6            7020.542275282852,-6534.512719466765,
ccc  7            5525.620597073791,-3580.432803129281,
ccc  8            1164.639428963841,
ccc  9          0.0,0.0,0.0,0.0,0.0,0.0/
ccc
      DATA Zp21/0.9898214418809327, 0.9594929736144974,
     1          0.9096319953545184, 0.8412535328311812,
     2          0.7557495743542583, 0.6548607339452851,
     3          0.5406408174555976, 0.4154150130018864,
     4          0.2817325568414297, 0.1423148382732851,
     5          0.0000000000000000,
     6         -0.1423148382732851,-0.2817325568414297,
     7         -0.4154150130018864,-0.5406408174555976,
     8         -0.6548607339452851,-0.7557495743542583, 
     9         -0.8412535328311812,-0.9096319953545184,
     a         -0.9594929736144974,-0.9898214418809327/
c
      DATA Wp21k1/40.91258980361040,-69.11503837897816,
     1            65.80507790523560,-69.11503837898373,
     2            67.78308420797106,-69.11503837899778,
     3            68.30792716563759,-69.11503837900795,
     4            68.49459295516724,-69.11503837901213,
     5            68.54383971474920,
     6           -69.11503837901213, 68.49459295516724,
     7           -69.11503837900795, 68.30792716563759,
     8           -69.11503837899778, 67.78308420797106,
     9           -69.11503837898373, 65.80507790523560,
     a           -69.11503837897816, 40.91258980361040/
ccc
ccc   DATA Wp21k3/6364.91821744174,-19267.83229098791,
ccc  1           29317.26172550868,-34478.42376106038,
ccc  2           37049.70046340271,-38499.30026119479,
ccc  3           39357.74970048209,-39887.54536375314,
ccc  4           40205.64256060925,-40378.03411697539,
ccc  5           40431.72625305376,
ccc  6          -40378.03411697539, 40205.64256060925,
ccc  7          -39887.54536375314, 39357.74970048209,
ccc  8          -38499.30026119479, 37049.70046340271,
ccc  9          -34478.42376106038, 29317.26172550868,
ccc  a          -19267.83229098791,  6364.91821744174/
c
ccc   DATA Pi/3.141592653589793D0/
ccc   SAVE Pi, Zp21, Wp21k1, Wp21l3
******7***************************************************************72
* In Pajunen method, there is a option to choose different points of 
* weight. To do this, just change the variable names for the Zp* & Wp*
* e.g., in order to use only 15 points
*                   Zi(i) = Zp15(i) and Wi(i) = -Wp15k1(i), 
* and to use only 11 points 
*                   Zi(i) = Zp11(i) and Wi(i) = -Wp11k1(i)
******7***************************************************************72  
c*** Zero integral arrays
      kmx= k1
      kmx= 1
      IF(k1 .eq. -1) kmx= 0
      DO  k= k1,kmx
          Jntegral(0,k)= 0.d0
          DO  IPV= POTPARI(ISTATE),HPARF(ISTATE)
              Jntegral(IPV,k)= 0.0d0
              ENDDO
          ENDDO
c
c*** Begin quadrature loop for sums over M mesh points
      DO  i= 1,M
c ... first get potential and derivatives at Pajunen k=1 point
          RDIST= 0.5d0*(Ro+Ri + (Ro-Ri)*Zp21(i))
          CALL VGEN(ISTATE,RDIST,VDIST,BETADIST,IDAT)
          EMinusV= DABS(EO- VDIST)
          DO  k= k1, kmx
c ... loop over posible  k  values
              IF(k .EQ. -1) THEN
c* For k= -1, using Pajunen k=1 quadrature method ......
                  Fzt= DSQRT(EminusV/(1-Zp21(i)**2))* (1-Zp21(i)**2)**2
                ELSEIF(k .EQ. 0) THEN 
c* For k= 0, using Pajunen k=1 quadrature method ......
                  Fzt= DSQRT((1-Zp21(i)**2)/EminusV) * (1-Zp21(i)**2)
c* For k= 1, use Pajunen quadrature method [JCP, 71(6), 2618, 1979]
                ELSEIF(k .EQ. 1) THEN
c* If  k = 1 points are used
                  Fzt= ( DSQRT((1-Zp21(i)**2)/EminusV) )**3
                ENDIF
c*** Accumulate the phase integrals
              Jntegral(0,k)= Jntegral(0,k) - Wp21k1(i)* Fzt
ccc           IF(k .gt. k1) THEN
                  DO IPV= POTPARI(ISTATE),HPARF(ISTATE)
c ... Loop over free parameters accumulating phase integral derivs.
                      Jntegral(IPV,k)= Jntegral(IPV,k)- Wp21k1(i)* Fzt
     1                                       *(dVdPk(IPV)- DEDPK(IPV))
                      ENDDO
ccc               ENDIF
c .... end of loop over  k values
              ENDDO
c ... end of quadrature loop over mesh points ....................
          ENDDO
      DO  k= k1,kmx
          Jntegral(0,k)= 0.25d0*(Ro-Ri) * Jntegral(0,k)
          DO  IPV= POTPARI(ISTATE),HPARF(ISTATE)
c*** In Pajunen's method, a factor of (1/2) is incorporated because 
c    the phase integral is a contour integral
              Jntegral(IPV,k)= 0.25d0*(Ro-Ri) * Jntegral(IPV,k)
              ENDDO
          ENDDO
c????????????????????????
c** Test that 'regular' & Pajunen integrals agree for k= -1 & 0
cc  and that Pajunen  k=1 & 3 integrals agree for k=1 case
cc    DO  IPV= POTPARI(ISTATE),HPARF(ISTATE)
cc        WRITE(32, 320) EO, (Jntegral(IPV,k), k= k1, kmx)
cc        WRITE(32, 321)  (Jntegral(IPV,k), k= k1, kmx)
cc        Do  k= k1, kmx
cc            tst(k)= 0.d0
cc            IF(DABS(Ink(IPV,k)) .gt. 0.d0) 
cc   1                         tst(k)= Jntegral(IPV,k)/Ink(IPV,k)-1.d0
cc            enddo
cc        write(32, 322) (tst(k), k= k1,kmx)
cc        ENDDO
cc320 FORMAT('  At   E=',F12.6,1P3D16.8)
cc321 FORMAT(21x,1P3D16.8)
cc322 FORMAT(21x,3(1PD10.2,6x))
c????????????????????????
      RETURN
      END      
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

