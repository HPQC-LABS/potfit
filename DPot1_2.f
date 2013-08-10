c***********************************************************************
        PROGRAM DPotFit
c***********************************************************************
c** Program "D(iatomic)Pot(ential)Fit" (DPotFit) for performing least-
c   squares fits of diatomic spectral data to molecular potential
c   energy functions for one or multiple electronic states.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++       COPYRIGHT 2006  by R.J. Le Roy, J.Y. Seto and Y. Huang     +++
c   Dept. of Chemistry, Univ. of Waterloo, Waterloo, Ontario, Canada   +
c    This software may not be sold or any other commercial use made    +
c      of it without the express written permission of the authors.    +
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
c++ Uses least-squares subroutine NLLSSRR written by Le Roy & Dulick +++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
c** This program can perform the following types of calculations:
c (i)  From a set of read-in constants, make predictions for any chosen
c    input data set consisting of diatomic singlet-singlet transitions,
c    and calculate deviations [calc.-obs.]
c (ii)  Fit a data set made up of any combination of MW, IR or
c    electronic vibrational bands, and/or fluorescence series, involving
c    one or more electronic states and one or more isotopomers, to
c    parameters defining the observed levels of each state.
c=======================================================================
c** Dimensioning parameters intrinsic to the program are input through
c          'arrsizes.h'
c** Parameters characterizing the problem and governing the fits are
c   read on  channel-5  while the raw data are read on  channel-4 .
c   Principle output goes to  channel-6  while higher channel numbers
c   are used for secondary or more detailed/voluminous output.
c***********************************************************************
      INCLUDE 'arrsizes.h'
      INCLUDE 'BLKISOT.h'
      INCLUDE 'BLKDATA.h'
      INCLUDE 'BLKPOT.h'
      INCLUDE 'BLKBOB.h'
      INCLUDE 'BLKCOUNT.h'
c-----------------------------------------------------------------------
      CHARACTER*40 DATAFILE,MAKEPRED
      CHARACTER*24 WRITFILE,TVNAME(NPARMX)
      CHARACTER*23 FN4,FN6,FN7,FN8,FN10,FN11,FN12,FN13,FN14,FN15,FN16
cc   1                                                           ,FN32
      INTEGER*4 lnblnk
      INTEGER I,J,ISTATE,IISTP,CHARGE,IPV,MKPRED,PRINP,
     1 JTRUNC(NSTATEMX),EFSEL(NSTATEMX),PASok(NSTATEMX),OSEL(NSTATEMX),
     2 NDAT(0:NVIBMX,NISTPMX,NSTATEMX)
      REAL*8 UCUTOFF,ZMASE,DECM(NSTATEMX)
c
      INTEGER NOWIDTHS
      COMMON /WIDTHBLK/NOWIDTHS
c                                                      
c** Parameters required for NLLSSRR.
c
      INTEGER NPTOT,IROUND,ROBUST,LPRINT,SIROUND,NFPAR,IFXPV(NPARMX),
     1  SIFXPV(NPARMX)
      REAL*8 PV(NPARMX),PU(NPARMX),PS(NPARMX),CM(NPARMX,NPARMX),
     1  PUSAV(NPARMX),PSSAV(NPARMX),TSTPS,TSTPU,DSE
c-----------------------------------------------------------------------
c** Common block carrying term values and term-value count labels
      REAL*8 TVALUE(NPARMX)
      INTEGER NSTATES,NTVALL(0:NSTATEMX),NTVI(NSTATEMX),NTVF(NSTATEMX),
     1  VMIN(NSTATEMX),VMAX(NSTATEMX)
      COMMON /NLSBLK/TVALUE,NSTATES,NTVALL,NTVI,NTVF,VMIN,VMAX
c
c** Set type statements for unused MASSES variables.
c
      CHARACTER*2  CATOM
      INTEGER GELGS(2,NISTPMX),GNS(2,NISTPMX)
      REAL*8 ABUND(2,NISTPMX)
c------------------------------------------------------------------------
      REAL*8 RDIST,VDIST,PHIDIST
c
      REAL*8 Vsr(NPNTMX,NSTATEMX),Bsr(NPNTMX,NSTATEMX)
      INTEGER nPointSR(NSTATEMX)
      COMMON /VsrBLK/Vsr,Bsr,nPointSR
c
      REAL*8 Plr(NPNTMX,NSTATEMX),Blr(NPNTMX,NSTATEMX)
      INTEGER nPointLR(NSTATEMX)
      COMMON /PlrBLK/Plr,Blr,nPointLR
c**************************
      REAL*8 RR(NPNTMX),RM2(NPNTMX),VV(NPNTMX),VLIMT
      INTEGER NCNN
      DATA ZMASE /5.4857990945D-04/
      DATA MAKEPRED/'MAKEPRED                                    '/
c=======================================================================
      SLABL(-3)='BV'
      SLABL(-2)='WI'
      SLABL(-1)='PA'
      SLABL(0)='FS'
c=======================================================================
      NFS1= 0
      DO  I=1,NPARMX
          TVALUE(I)= 0.d0
          PV(I)= 0.0d0
          PU(I)= 0.0d0
          PS(I)= 0.0d0
          IFXPV(I)= 1
          ENDDO
      SIROUND= 0
c=======================================================================
c** Start by reading parameters describing the overall nature of the
c   case and placing chosen restrictions on the data set to be used.
c
c  AN(1) & AN(2) are atomic numbers identifying the atoms forming the
c          molecule.
c
c  CHARGE (+/- integer) is the charge on the molecule (=0 for neutral).
c  If(CHARGE.ne.0) use Watson's(JMS 1980) charge-modified reduced mass.
c
c  NISTP   is the number of isotopomers to be simultaneously considered.
c
c  NSTATES  is the number of electronic states associated with the data
c        set to be analysed:  NSTATES = 1  for fits to IR/MW and/or
c        fluorescence data for a single electronic state, while
c        NSTATES > 1  for multi-state fits.
c        Upper states of fluorescence series NOT included in this count.
c
c  DATAFILE  is the (character variable) name of the file containing the
c     experimental data to be used in the fit.  If it is not located in
c     the current directory, the name 'DATAFILE' must include the
c     relative path.  The valiable name may (currently) consist of up to
c     40 characters.  READ ON A SEPARATE LINE!
c
c !! To make predictions using a completely specified set of parameters,
c      the input value of parameter DATAFILE must be  'MAKEPRED'
c
c  WRITFILE  is the (character variable) name of the file to which the
c     output will be written.  Channel-6 outut goes to  WRITFILE.6,
c     channel-7 output to WRITFILE.7, channel-8 to WRITFILE.8, ... etc.
c     If not in the current directory, the name 'WRITFILE' must include the
c     relative path.  The valiable name may (currently) consist of up to
c     40 characters, enclosed in single quotes, with no leading spaces.
c=======================================================================
      READ(5,*)  AN(1), AN(2), CHARGE, NISTP, NSTATES
      READ(5,*)  DATAFILE
      READ(5,*)  WRITFILE
c=======================================================================
c** These statements construct and define the names of output files
c   associated with WRITE's to channels 6-10 used by the program.
      WRITE(FN6,*) WRITFILE(1:lnblnk(WRITFILE)),'.6'
      WRITE(FN7,*) WRITFILE(1:lnblnk(WRITFILE)),'.7'
      WRITE(FN8,*) WRITFILE(1:lnblnk(WRITFILE)),'.8'
      OPEN(UNIT= 6, FILE= FN6)
      OPEN(UNIT= 7, FILE= FN7)
      OPEN(UNIT= 8, FILE= FN8)
      MKPRED= 0
      IF(DATAFILE.EQ.MAKEPRED) THEN
          MKPRED= 1
          ENDIF
c-----------------------------------------------------------------------
c  UCUTOFF   Neglect any input data with uncertainties > UCUTOFF (cm-1)
c
c  NOWIDTHS >  0  causes the program to ignore any tunneling widths in
c                 the data set.
c           <= 0  causes the program to fit to tunneling widths
c           <  0  use simple version of dWdP, ignoring the partial 
c                 derivative of t_vib which involves k = 1 phase integral
c  IROUND  specifies the level of rounding inside NLLSSRR if:
c          > 0 : requires that Sequential Rounding & Refitting be
c                performed, with each parameter being rounded at the
c                IROUND'th sig. digit of its local uncertainty.
c          <=0 : simply stops after full convergence (without rounding).
c
c  ROBUST > 0  (integer) causes "Robust" least-squares weighting (as per
c              Watson [J.Mol.Spectrosc. 219, 326 (2003)] to be used
c         = 0  uses normal data weights  1/[uncertainty(i)]**2
c
c  LPRINT  specifies the level of printing inside NLLSSRR if:
c            = 0 : no print except for failed convergence.
c           <  0 : only converged, unrounded parameters, PU & PS's
c           >= 1 : print converged parameters, PU & PS's
c           >= 2 : also print parameter change each rounding step
c           >= 3 : also indicate nature of convergence
c           >= 4 : also print convergence tests on each cycle
c           >= 5 : also parameters changes & uncertainties, each cycle
c
c  PRINP > 0  causes a summary of the input data to be printed before
c            the fitting starts.  Normally set =0 unless troubleshooting
c=======================================================================
      READ(5,*)  UCUTOFF, NOWIDTHS, IROUND, ROBUST, LPRINT, PRINP
c=======================================================================
c!!!
cc    IF(NOWIDTHS.LE.0) THEN
cc        WRITE(FN32,*) WRITFILE(1:lnblnk(WRITFILE)),'.32'
cc        OPEN(UNIT= 32, FILE= FN32)
cc        ENDIF
c!!!  
      I= 999
      WRITE(6,601) NISTP
      DO  IISTP= 1,NISTP
c
c** Read the mass numbers of the atoms in each of the isotopomers
c
c  MN(i,IISTP)  is the mass number for atom with atomic number AN(i)
c       [NOTE: be sure order of MN values consistent with that of AN's].
c       Choosing it .ne. value for some known isotope if that species
c       causes the average atomic mass to be used.
c=======================================================================
          READ(5,*) MN(1,IISTP), MN(2,IISTP)
c=======================================================================
          I= MIN(I,MN(1,IISTP),MN(2,IISTP))
          CALL MASSES(AN(1),MN(1,IISTP),CATOM,GELGS(1,IISTP),
     1                     GNS(1,IISTP),ZMASS(1,IISTP),ABUND(1,IISTP))
          IF (IISTP.EQ.1) NAME(1)= CATOM
          CALL MASSES(AN(2),MN(2,IISTP),CATOM,GELGS(2,IISTP),
     1                     GNS(2,IISTP),ZMASS(2,IISTP),ABUND(2,IISTP))
          IF (IISTP.EQ.1) NAME(2)= CATOM
          ZMASS(3,IISTP)= (ZMASS(1,IISTP)*ZMASS(2,IISTP))/
     1                    (ZMASS(1,IISTP)+ZMASS(2,IISTP)-CHARGE*ZMASE)
          WRITE(6,602) NAME(1),MN(1,IISTP),NAME(2),MN(2,IISTP),
     1                                          (ZMASS(J,IISTP),J=1,3)
          RSQMU(IISTP)= DSQRT(ZMASS(3,1)/ZMASS(3,IISTP))
          ENDDO
c... end of loop over isotopologues ....................................
      IF(I.EQ.0) WRITE(6,603)
      IF(CHARGE.NE.0) WRITE(6,597) CHARGE
      WRITE(6,599) DATAFILE
cc    WRITE(6,*)
      IF(AN(1).EQ.AN(2)) WRITE(6,604)
  599 FORMAT(/' Use experimental data input file:  ',a30)
  597 FORMAT(1x,67('-')/' Since this is an ion with charge',SP,i3,
     1  ", use Watson's charge-modified reduced mass.")
  601 FORMAT(2X,'Input data for',I3,'  isotopomer(s)'/2X,16('**')/2X,
     1  '    Isotopomer       Mass of atom-1   Mass of atom-2    Reduced
     2 mass'/ 2X,'----------------- ',3('   --------------'))
  602 FORMAT(2X,A2,'(',I3,') - ',A2,'(',I3,')',3(3X,F14.9))
  603 FORMAT('  Note that   (Mass Number) = 0   causes the average atomi
     1c mass to be used.')
  604 FORMAT('  For homonuclear diatomics, BO correction terms are the s
     1ame for both atoms.'/'  Only the first set of correction terms wil
     2l be used, UA(R) and TA(R),'/'    with a mass scaling factor equal
     3 to the sum of the two individual terms.')
c
      DO  ISTATE= 1,NSTATES 
c-----------------------------------------------------------------------
c** Read parameters to characterize state & possibly restrict data used
c  SLABL(s)  is a 2-character alphameric label enclosed in single quotes
c            to identify the electronic state; e.g., 'X0', 'A1', ... etc.
c  IOMEG(s) .GE.0  is electronic angular momentum of singlet state with
c                    projection quantum number  Lambda= IOMEG
c  IOMEG(s) .LT.0  if it indicates a doublet SIGMA electronic state
c                 [other spin multiplets not yet coded]
c  V(MIN/MAX)(s) Neglect data for electronic state vibrational levels
c                outside the range  VMIN  to  VMAX.
c  JTRUNC(s)     data with J > JTRUNC are not included in the fit.
c  EFSEL(s)  allows a user to consider data for:
c          * ONLY the e-parity levels of this state, if EFSEL > 0
c          * ONLY the f-parity levels of this state, if EFSEL < 0
c          * BOTH e- and f-parity levels of thsi state, if EFSEL = 0
c=======================================================================
          READ(5,*) SLABL(ISTATE), IOMEG(ISTATE), VMIN(ISTATE), 
     1                     VMAX(ISTATE), JTRUNC(ISTATE), EFSEL(ISTATE)
c=======================================================================
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          CALL READPOT(ISTATE,SLABL,OSEL)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** These statements construct and define the names of output files
c   associated with WRITE's to channels 6-10 used by the program.
          IF(OSEL(ISTATE).GT.0) THEN
              WRITE(FN10,*) WRITFILE(1:lnblnk(WRITFILE)),'.10'
              WRITE(FN11,*) WRITFILE(1:lnblnk(WRITFILE)),'.11'
              OPEN(UNIT=10, FILE= FN10)
              OPEN(UNIT=11, FILE= FN11)
              IF(NUA(ISTATE).GE.0) THEN
                  WRITE(FN12,*) WRITFILE(1:lnblnk(WRITFILE)),'.12'
                  OPEN(UNIT=12, FILE= FN12)
                  ENDIF
              IF(NUB(ISTATE).GE.0) THEN
                  WRITE(FN13,*) WRITFILE(1:lnblnk(WRITFILE)),'.13'
                  OPEN(UNIT=13, FILE= FN13)
                  ENDIF
              IF(NTA(ISTATE).GE.0) THEN
                  WRITE(FN14,*) WRITFILE(1:lnblnk(WRITFILE)),'.14'
                  OPEN(UNIT=14, FILE= FN14)
                  ENDIF
              IF(NTB(ISTATE).GE.0) THEN
                  WRITE(FN15,*) WRITFILE(1:lnblnk(WRITFILE)),'.15'
                  OPEN(UNIT=15, FILE= FN15)
                  ENDIF
              IF((IOMEG(ISTATE).NE.0).AND.(NwCFT(ISTATE).GE.0)) THEN
                  WRITE(FN16,*) WRITFILE(1:lnblnk(WRITFILE)),'.16'
                  OPEN(UNIT=16, FILE= FN16)
                  ENDIF
              ENDIF
          PASok(ISTATE)= 1
          IF(PSEL(ISTATE).EQ.4) PASok(ISTATE)= 0
          IF(PSEL(ISTATE).EQ.3) CALL VGEN(ISTATE,1.0d0,VDIST,PHIDIST,0)
          ENDDO
c** Now write summary of the initial potential parameters for each state
      CALL WRITEPOT(1,NSTATES,SLABL,NAME,DECM,PU,PS)
c
c** Now ... count potential parameters of various types for each state
c=======================================================================
c** Counters for numbers of potential parameters of different types for 
c   each state
c     COMMON /BLKCOUNT/TOTPOTPAR,POTPARI,POTPARF,UAPARI,UAPARF,
c    1  UBPARI,UBPARF,TAPARI,TAPARF,TBPARI,TBPARF,LDPARI,LDPARF
c=======================================================================
      TOTPOTPAR= 0
      IPV= 0
      DO 90 ISTATE= 1,NSTATES
          IF(PSEL(ISTATE).LE.0) GOTO 90
c... For all potential types, count  Re
          IPV= IPV+ 1
          POTPARI(ISTATE)= IPV
          POTPARF(ISTATE)= IPV
          IFXPV(IPV)= IFXRe(ISTATE)
          UAPARI(ISTATE)= 0
          UAPARF(ISTATE)= 0
          UBPARI(ISTATE)= 0
          UBPARF(ISTATE)= 0
          TAPARI(ISTATE)= 0
          TAPARF(ISTATE)= 0
          TBPARI(ISTATE)= 0
          TBPARF(ISTATE)= 0
          LDPARI(ISTATE)= 0
          LDPARF(ISTATE)= 0
          HPARF(ISTATE)= 0
c... For all cases, except GPEF (where it doesn't appear), count De
          IF(PSEL(ISTATE).NE.4) THEN
              IFXPV(IPV)= IFXDe(ISTATE)
              IPV= IPV+ 1
              POTPARF(ISTATE)= IPV
              IFXPV(IPV)= IFXRe(ISTATE)
              ENDIF
          IF(PSEL(ISTATE).EQ.2) THEN
c... For MLJ/MLR family, count long-range parameters: first count Cn
              DO  J= 1,NCMM(ISTATE)
                  IPV= IPV+ 1
                  POTPARF(ISTATE)= IPV
                  IFXPV(IPV)= IFXCm(J,ISTATE)
c... additional Aubert-Frecon{3,6} parameters included in this loop
                  ENDDO
              ENDIF
c... Now count [exponent] \phi_i expansion coefficients
          DO J= 0,MAX0(NSphi(ISTATE),NLphi(ISTATE))
              IPV= IPV+ 1
              POTPARF(ISTATE)= IPV
              IFXPV(IPV)= IFXPHI(J,ISTATE)
              ENDDO
          IF(NUA(ISTATE).GE.0) THEN
c... Count adiabatic parameters for atom A (if appropriate)
              UAPARI(ISTATE)= IPV + 1
              DO  J= 0,NUA(ISTATE)
                  IPV= IPV+ 1
                  UAPARF(ISTATE)= IPV
                  IFXPV(IPV)= IFXUA(J,ISTATE)
                  ENDDO
              ENDIF
          IF(NUB(ISTATE).GE.0) THEN
c... Count adiabatic parameters for atom B (if appropriate)
              UBPARI(ISTATE)= IPV + 1
              DO  J= 0,NUB(ISTATE)
                  IPV= IPV+ 1
                  UBPARF(ISTATE)= IPV
                  IFXPV(IPV)= IFXUB(J,ISTATE)
                  ENDDO
              ENDIF
          IF(NTA(ISTATE).GE.0) THEN
c... Count centrifugal BOB parameters for atom A (if appropriate)
              TAPARI(ISTATE)= IPV + 1
              DO  J= 0,NTA(ISTATE)
                  IPV= IPV+ 1
                  TAPARF(ISTATE)= IPV
                  IFXPV(IPV)= IFXTA(J,ISTATE)
                  ENDDO
              ENDIF
          IF(NTB(ISTATE).GE.0) THEN
c... Count centrifugal BOB parameters for atom B (if appropriate)
              TBPARI(ISTATE)= IPV + 1
              DO  J= 0,NTB(ISTATE)
                  IPV= IPV+ 1
                  TBPARF(ISTATE)= IPV
                  IFXPV(IPV)= IFXTB(J,ISTATE)
                  ENDDO
              ENDIF
          IF(NwCFT(ISTATE).GE.0) THEN
c... Count Lambda/doublet-sigma doubling parameters (if appropriate)
              LDPARI(ISTATE)= IPV + 1
              DO  J= 0,NwCFT(ISTATE)
                  IPV= IPV+ 1
                  LDPARF(ISTATE)= IPV
                  IFXPV(IPV)= IFXwCFT(J,ISTATE)
                  ENDDO
              ENDIF
          HPARF(ISTATE)= IPV
   90     CONTINUE
      TOTPOTPAR= IPV
      IF(TOTPOTPAR.GT.HPARMX) THEN
          WRITE(6,626) TOTPOTPAR,HPARMX
          STOP
          ENDIF
      NPTOT= TOTPOTPAR
      NFPAR= 0
c** Count total free Hamiltonian fitting parameters
      DO  IPV= 1, TOTPOTPAR
          IF(IFXPV(IPV).LE.0) NFPAR= NFPAR+ 1
          ENDDO
c------------ Finished counting Hamiltonian Parameters------------------
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Call subroutine to input experimental data in specified
c   band-by-band, and do bookkeeping to characterize amounts of data or
c   each type.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IF(MKPRED.LE.0) OPEN(UNIT= 4, STATUS= 'OLD', FILE= DATAFILE)
c** when COMMON blocks check out ... introduce MKPRED option ......
      IF(MKPRED.GT.0) THEN
          WRITE(FN4,*) WRITFILE(1:lnblnk(WRITFILE)),'.4'
          OPEN(UNIT= 4, FILE= FN4)
          IF(UCUTOFF.LT.1.d0) UCUTOFF= 1.d0
          CALL MKPREDICT(NSTATES,NDAT)
          REWIND(4)
          ENDIF
      CALL READATA(NSTATES,PASok,UCUTOFF,JTRUNC,EFSEL,VMIN,VMAX,NDAT,
     &                                                 NOWIDTHS,PRINP)
      NTVALL(0)= 0
      DO  ISTATE= 1,NSTATES
          IF(PSEL(ISTATE).EQ.-2) THEN
c... If this state to be represented by term values, determine the number
c    and add them to the parameter count
              NTVI(ISTATE)= NPTOT+ 1
              CALL TVSORT(ISTATE,NPTOT,VMAX,NTVALL,TVNAME)
              NTVALL(0)= NTVALL(0) + NTVALL(ISTATE)
              IF(NTVALL(ISTATE).GT.0) THEN
                  NTVF(ISTATE)= NPTOT
                  ENDIF
              ENDIF
          ENDDO
c** Add number of fluorescence series origins to total parameter count
c   and set initial values of any fluorescence series origins to zero.
      IF(NFSTOT.GT.0) THEN
          NFS1= NPTOT+ 1
          NPTOT= NPTOT+ NFSTOT
          ENDIF
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Call subroutine that will generate the potential function and its 
c   partial derivatives w.r.t. any free parameters for each state.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      DO  ISTATE=1, NSTATES
          IF(PSEL(ISTATE).EQ.0) THEN
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Call subroutine PREPOT that will generate the potential function for
c  each state for which perform forward calculation from known potential
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
              DO  I= 1, NDATPT(ISTATE)
                  RR(I)= RMIN(ISTATE)+ DBLE(I)* RH(ISTATE)- RH(ISTATE)
                  RM2(I)= 1.0d0 / (RR(I) ** 2)
                  ENDDO
              CALL PREPOT(1,AN(1),AN(2),MN(1,1),MN(2,1),NDATPT(ISTATE),
     1                                           RR,RM2,VLIMT,VV,NCNN)
              DO  I= 1, NDATPT(ISTATE)
                  VPOT(I,ISTATE)= VV(I) - DE(ISTATE) + VLIM(ISTATE)
                  ENDDO
              J= 0.05d0/RH(ISTATE)
              DO  I= 1,NDATPT(ISTATE),J
                  WRITE(17,900) RR(I),VPOT(I,ISTATE)
                  ENDDO
            ENDIF
          ENDDO
c** Set the energy convergence criterion to be 1/100th of the smallest
c   experimental uncertainty. [UCUTOFF reset by READATA to that min. unc.]
      DO  ISTATE=1,NSTATES
          EPS(ISTATE)= DMIN1(UCUTOFF/100.0d0,1.d-06)
          WRITE(6,638) SLABL(ISTATE), EPS(ISTATE)
c** Initialize the dissociation energy ????
          DECM(ISTATE)= 0.0d0
          ENDDO
      IF(IROUND.NE.0) WRITE(6,685) ABS(IROUND)
      IF(IROUND.GT.0) WRITE(6,686)
      IF(IROUND.LT.0) WRITE(6,687)
      IF(ROBUST.GT.0) THEN
          ROBUST= 2
          WRITE(6,596)
        ELSE
          WRITE(6,598)
        ENDIF
  596 FORMAT(/" Fit uses Watson's",' "Robust" data weighting [J.Mol/Spec
     1trosc. 219, 326 (2003)] '/20x,'1/[{unc(i)}^2 + {calc.-obs.}^2/3]')
  598 FORMAT(/' Fit uses standard  1/[uncertainty(i)]**2  data weighting
     1')
  626 FORMAT(/' *** Dimension Error *** [(total No. Hamiltonian parmaete
     1rs)=',i4,'] >  HPARMX=',i4)
  638 FORMAT(' State ',A2,'  Energy Convergence criterion EPS is',
     1  1PD8.1,' cm-1')
  685 FORMAT(/' Apply "Sequential Rounding & Refitting" at digit-',
     1  i1,' of the (local) parameter')
  686 FORMAT(4x,'uncertainty, selecting remaining parameter with largest
     1 relative uncertainty')
  687 FORMAT(4x,'uncertainty, proceeding sequentially from the LAST para
     1meter to the FIRST.')
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Now Generate internal NLLSSRR variables {PV} from the external ones 
      CALL MAPPAR(NSTATES,PV,0)
      SIROUND= IROUND
      IROUND= 0
      IF((NFSTOT.GT.0).OR.(NTVALL(0).GT.0)) THEN
c** If HAVE fluorescence series ... first fix ALL potential parameters
c  and fit to determine the series origins, and only THEN free potential
c  parameters too.  First, save read-in values of  'IF(fix)' parameters
          DO  I= 1,TOTPOTPAR
              SIFXPV(I)= IFXPV(I)
              IFXPV(I)= 1
              ENDDO
          DO  I= TOTPOTPAR+1,NPTOT
              IFXPV(I)= 0
              PV(I)= 0.d0
              ENDDO
c** Call NLLSSRR to get Fluorescence series origins ....
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          CALL NLLSSRR(COUNTOT,NPTOT,NPARMX,IROUND,ROBUST,LPRINT,IFXPV,
     1                FREQ,UFREQ,DFREQ,PV,PU,PS,CM,TSTPS,TSTPU,DSE)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c...  Now reset "IFX(fix)" parameters to read-in values ... & proceed ..
          DO  I= 1,TOTPOTPAR
              IFXPV(I)= SIFXPV(I)
              ENDDO
c** Now, set TVALUE values & reset parameter array for global fit
          DO  I= TOTPOTPAR+1,NPTOT
              TVALUE(I-TOTPOTPAR)= PV(I)
              IFXPV(I)= 0
              ENDDO
          NFPAR= NFPAR+ NFSTOT+ NTVALL(0)
          CALL MAPPAR(NSTATES,PV,0)
          ENDIF
c--- End of section to determine preliminary values of any fluorescence 
c    series origins
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Call subroutine NLLSSRR to calculate converged parameters from trial
c   values and spectroscopic data.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      CALL NLLSSRR(COUNTOT,NPTOT,NPARMX,IROUND,ROBUST,LPRINT,IFXPV,
     1                   FREQ,UFREQ,DFREQ,PV,PU,PS,CM,TSTPS,TSTPU,DSE)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IF(SIROUND.NE.0) THEN
c** If SRR rounding is to be performed, first save global uncertainties
          DO  I= 1, NPTOT
              PUSAV(I)= PU(I)
              PSSAV(I)= PS(I)
              ENDDO
c** Perform group rounding of fluorescence series origins in single step
          IF(NFSTOT.GT.0) THEN
              IROUND= IABS(SIROUND) + 1
              CALL GPROUND(IROUND,NPTOT,NPARMX,TOTPOTPAR+1,NPTOT,
     1                                             LPRINT,IFXPV,PV,PU)
              ENDIF
c ... and then call NLLSSRR again to sequentially round remaining parm.
          IROUND= SIROUND
          CALL NLLSSRR(COUNTOT,NPTOT,NPARMX,IROUND,ROBUST,LPRINT,IFXPV,
     1                   FREQ,UFREQ,DFREQ,PV,PU,PS,CM,TSTPS,TSTPU,DSE)
c ... finally, reset all parameter uncertainties at pre-rounding values
          DO  I= 1, NPTOT
              PU(I)= PUSAV(I)
              PS(I)= PSSAV(I)
              ENDDO
ccc       DSE= DSE*DSQRT(DFLOAT(COUNTOT- (NFPAR- NFSTOT- NTVALL(0))/
          DSE= DSE*DSQRT(DFLOAT(COUNTOT- (NFPAR- NFSTOT))/
     1                                         DFLOAT(COUNTOT- NFPAR))
          ENDIF
c** Writing out the general information of the fit.
c-----------------------------------------------------------------------
      WRITE(6,690)
      WRITE(6,691) NFPAR,COUNTOT,DSE
c-----------------------------------------------------------------------
c** Writing out the fluorescence band results.
c-----------------------------------------------------------------------
      IF(NFSTOT.GT.0) THEN
          WRITE(6,690)
          WRITE(6,692) NFSTOT
          J= NPTOT - NFSTOT
          DO  I= 1,NFSTOT
              WRITE(6,694) VP(FSBAND(I)),VPP(FSBAND(I)),
     1            EFP(IFIRST(FSBAND(I))),ISTP(FSBAND(I)),TVALUE(J+I),
     2                                                 PU(J+I),PS(J+I)
              ENDDO
          ENDIF
      DO  ISTATE= 1, NSTATES
          IF(PSEL(ISTATE).EQ.-2) THEN
c** For states represented by independent term values for each level ...
              WRITE(6,690)
              WRITE(6,696) SLABL(ISTATE),NTVALL(ISTATE)
              WRITE(6,698) (TVNAME(I),PV(I),PU(I),PS(I),I= 
     1                                      NTVI(ISTATE),NTVF(ISTATE))
              ENDIF
          ENDDO
  690 FORMAT(/,1X,34('=='))
  691 FORMAT(' For fit of',I5,'  parameters to',I6,
     1  '  transitions,  DSE=',G15.8)
  692 FORMAT('  The following',I5,' Fluorescence Series Origins were det
     1ermined'/1x,31('--')/"  ( v', J', p'; ISTP)",5x,'T(value)',5x,
     2  'Uncertainty  Sensitivity'/1x,31('--'))
cc694 FORMAT(3X,'(',I3,',',I3,',',SP,I3,SS,';',I2,')',1X,1PD19.10,
  694 FORMAT(2X,'(',I4,',',I3,',',SP,I3,SS,';',I2,')',1X,1PD19.10,
     1  D11.1,D12.1)
  696 FORMAT('  State ',A2,' represented by the',I5,' individual term va
     1lues:'/1x,34('--')/" T(es: v', J', p';IS)  #dat",5x,'T(value)',5x,
     2  'Uncertainty  Sensitivity'/1x,34('--'))
  698 FORMAT(2X,A24,1PD19.10,D11.1,D12.1)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Calculation of the uncertainties for Te for each potential require 
c   elements from the correlation matrix.
c
      DO  ISTATE= 1,NSTATES
          IF((IFXDE(1).LE.0).AND.(IFXDE(ISTATE).LE.0)) THEN
              DECM(ISTATE)= CM(1,POTPARI(ISTATE))
            ELSE
              DECM(ISTATE)= 0.0d0
            ENDIF
          ENDDO
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Call subroutine MAPPAR to convert internal NLLSSRR parameter array
c   back into external (logical) variable system.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      CALL MAPPAR(NSTATES,PV,1)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Call subroutine VGEN to generate the potential function from the
c   final calculated converged parameters.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IF(PSEL(NSTATES).NE.0) THEN
          DO  ISTATE= 1,NSTATES
              IF(PSEL(ISTATE).GT.0) THEN
                  nPointSR(ISTATE)= RMIN(ISTATE)/RH(ISTATE) 
                  IF (OSEL(ISTATE) .NE. 0) THEN
                      IF (RMAX(ISTATE) .GT. 100.0) THEN
                          nPointLR(ISTATE)= 0
                        ELSE
                          nPointLR(ISTATE)= (100.0-RMAX(ISTATE))
     &                                        /(RH(ISTATE)*OSEL(ISTATE))
                        ENDIF
                      ENDIF
                  CALL VGEN(ISTATE,-1.0d0,VDIST,PHIDIST,0)
                  IF(OSEL(ISTATE).NE.0) THEN
                      DO  I= 1,nPointSR(ISTATE),OSEL(ISTATE)
c ... generate potential & exponent values in inner extrapolation region
                          RDIST= RH(1)*DBLE(I-1) 
                          CALL VGEN(ISTATE,RDIST,VDIST,PHIDIST,0)
                          Vsr(I,ISTATE)= VDIST
                          Bsr(I,ISTATE)= PHIDIST
                          ENDDO
                      DO  I= 1,nPointLR(ISTATE)
c ... generate potential & exponent values in outer extrapolation region
                          RDIST= RMAX(ISTATE) + 
     1                                   RH(ISTATE)*DBLE(I*OSEL(ISTATE))
                          CALL VGEN(ISTATE,RDIST,VDIST,PHIDIST,0)
                          Plr(I,ISTATE)= VDIST
                          Blr(I,ISTATE)= PHIDIST
                          ENDDO
                      ENDIF
                  ENDIF
              ENDDO
          ENDIF
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Call subroutine to print out a summary of the converged and fixed
c   values to standard output (channel-6).
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      CALL WRITEPOT(2,NSTATES,SLABL,NAME,DECM,PU,PS)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** If chosen, output file(s) will be created for the export of the
c   generated functions: V, PHIFX, UAR/UBR, or TAR/TBR and their
c   respective uncertainties.
      DO  ISTATE= 1, NSTATES
          IF(OSEL(ISTATE).GT.0) THEN
              IF(PSEL(ISTATE).GT.0) THEN
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Call subroutine to print out the generated functions and their
c   respective uncertainties.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                  CALL FUNUNC(ISTATE,WRITFILE,OSEL,PU,CM)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                ELSE 
                  J= 0.05d0/ RH(NSTATES)
                  DO  I= 1,NDATPT(NSTATES),J
                      WRITE(17,900) RR(I),VPOT(I,NSTATES)
                      ENDDO
                ENDIF
              ENDIF
          ENDDO
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Call subroutine to print out summary of dimensionless standard
c   errors on a band-by-band basis, and (if desired) print [obs.-calc.]
c   values to channel-8.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      CALL DIFFSTATS(NSTATES,ROBUST,MKPRED)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      STOP
  900 FORMAT(5X,G18.8,5X,G18.8)
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE MASSES(IAN,IMN,NAME,GELGS,GNS,MASS,ABUND)
c***********************************************************************
c** For isotope with (input) atomic number IAN and mass number IMN,
c  return (output):  (i) as the right-adjusted 2-character variable NAME
c  the alphabetic symbol for that element,  (ii) the ground state
c  electronic degeneracy GELGS, (iii) the nuclear spin degeneracy GNS,
c  (iv) the atomic mass MASS [amu], and  (v) the natural isotopic
c  abundance ABUND [in percent].   GELGS values based on atomic states
c  in Moore's "Atomic Energy Level" tables, the isotope masses are taken
c  from the 2003 mass table [Audi, Wapstra & Thibault, Nucl.Phys. A729,
c  337-676 (2003)] and other quantities from Tables 6.2 and 6.3 of 
c  "Quantities, Units and Symbols in Physical Chemistry", by Mills et 
c  al. (Blackwell, 2'nd Edition, Oxford, 1993).
c** If the input value of IMN does not equal one of the tabulated values
c  for atomic species IAN, return the abundance-averaged standard atomic
c  weight of that atom and set GNS=-1 and ABUND=-1.
c                          COPYRIGHT 2005
c** By R.J. Le Roy (with assistance from G.T. Kraemer & J.Y. Seto).
c                      Last modified  1 June 2005
c***********************************************************************
      REAL*8 zm(123,0:10),mass,ab(123,10),abund
      INTEGER i,ian,imn,gel(123),nmn(123),mn(123,10),ns2(123,10),
     1        gelgs, gns
      CHARACTER*2 NAME,AT(123)
c
      DATA  at(1),gel(1),nmn(1),(mn(1,i),i=1,3)/' H',2,3,1,2,3/
      DATA  (zm(1,i),i=0,3)/1.00794d0, 1.00782503207d0, 2.0141017778d0,
     1                      3.0160492777d0/
      DATA  (ns2(1,i),i=1,3)/1,2,1/
      DATA  (ab(1,i),i=1,3)/99.985d0,0.015d0,0.d0/
c
      DATA  at(2),gel(2),nmn(2),(mn(2,i),i=1,2)/'He',1,2,3,4/
      DATA  (zm(2,i),i=0,2)/4.002602d0, 3.0160293191d0, 4.00260325415d0/
      DATA  (ns2(2,i),i=1,2)/1,0/
      DATA  (ab(2,i),i=1,2)/0.000137d0,99.999863d0/
c
      DATA  at(3),gel(3),nmn(3),(mn(3,i),i=1,2)/'Li',2,2,6,7/
      DATA  (zm(3,i),i=0,2)/6.941d0, 6.015122795d0, 7.01600455d0/
      DATA  (ns2(3,i),i=1,2)/2,3/
      DATA  (ab(3,i),i=1,2)/7.5d0,92.5d0/
c
      DATA  at(4),gel(4),nmn(4),(mn(4,i),i=1,1)/'Be',1,1,9/
      DATA  (zm(4,i),i=0,1)/9.012182d0, 9.0121822d0/
      DATA  (ns2(4,i),i=1,1)/3/
      DATA  (ab(4,i),i=1,1)/100.d0/
c
      DATA at(5),gel(5),nmn(5),(mn(5,i),i=1,2)/' B',2,2,10,11/
      DATA (zm(5,i),i=0,2)/10.811d0, 10.0129370d0, 11.0093054d0/
      DATA  (ns2(5,i),i=1,2)/6,3/
      DATA  (ab(5,i),i=1,2)/19.9d0,80.1d0/
c
      DATA at(6),gel(6),nmn(6),(mn(6,i),i=1,3)/' C',1,3,12,13,14/
      DATA (zm(6,i),i=0,3)/12.011d0, 12.d0, 13.0033548378d0, 
     1                      14.003241989d0/
      DATA  (ns2(6,i),i=1,3)/0,1,0/
      DATA  (ab(6,i),i=1,3)/98.90d0,1.10d0, 0.d0/
c
      DATA at(7),gel(7),nmn(7),(mn(7,i),i=1,2)/' N',4,2,14,15/
      DATA (zm(7,i),i=0,2)/14.00674d0, 14.0030740048d0, 15.0001088982d0/
      DATA (ns2(7,i),i=1,2)/2,1/
      DATA (ab(7,i),i=1,2)/99.634d0,0.366d0/
c
      DATA at(8),gel(8),nmn(8),(mn(8,i),i=1,3)/' O',5,3,16,17,18/
      DATA (zm(8,i),i=0,3)/15.9994d0, 15.99491461956d0, 16.99913170d0,
     1                      17.9991610d0/
      DATA (ns2(8,i),i=1,3)/0,5,0/
      DATA (ab(8,i),i=1,3)/99.762d0, 0.038d0, 0.200d0/
c
      DATA at(9),gel(9),nmn(9),(mn(9,i),i=1,1)/' F',4,1,19/
      DATA (zm(9,i),i=0,1)/18.9984032d0, 18.99840322d0/
      DATA (ns2(9,i),i=1,1)/1/
      DATA (ab(9,i),i=1,1)/100.d0/
c
      DATA at(10),gel(10),nmn(10),(mn(10,i),i=1,3)/'Ne',1,3,20,21,22/
      DATA (zm(10,i),i=0,3)/20.1797d0, 19.9924401754d0, 20.99384668d0,
     1                       21.991385114d0/
      DATA (ns2(10,i),i=1,3)/0,3,0/
      DATA (ab(10,i),i=1,3)/90.48d0, 0.27d0, 9.25d0/
c
      DATA at(11),gel(11),nmn(11),(mn(11,i),i=1,1)/'Na',2,1,23/
      DATA (zm(11,i),i=0,1)/22.989768d0, 22.9897692809d0/
      DATA (ns2(11,i),i=1,1)/3/
      DATA (ab(11,i),i=1,1)/100.d0/
c
      DATA at(12),gel(12),nmn(12),(mn(12,i),i=1,3)/'Mg',1,3,24,25,26/
      DATA (zm(12,i),i=0,3)/24.3050d0, 23.985041700d0, 24.98583692d0,
     1                       25.982592929d0/
      DATA (ns2(12,i),i=1,3)/0,5,0/
      DATA (ab(12,i),i=1,3)/78.99d0, 10.00d0, 11.01d0/
c
      DATA at(13),gel(13),nmn(13),(mn(13,i),i=1,1)/'Al',2,1,27/
      DATA (zm(13,i),i=0,1)/26.981539d0, 26.98153863d0/
      DATA (ns2(13,i),i=1,1)/5/
      DATA (ab(13,i),i=1,1)/100.d0/
c
      DATA at(14),gel(14),nmn(14),(mn(14,i),i=1,3)/'Si',1,3,28,29,30/
      DATA (zm(14,i),i=0,3)/28.0855d0, 27.9769265325d0, 28.976494700d0,
     1                       29.97377017d0/
      DATA (ns2(14,i),i=1,3)/0,1,0/
      DATA (ab(14,i),i=1,3)/92.23d0, 4.67d0, 3.10d0/
 
      DATA at(15),gel(15),nmn(15),(mn(15,i),i=1,1)/' P',4,1,31/
      DATA (zm(15,i),i=0,1)/30.973762d0, 30.97376163d0/
      DATA (ns2(15,i),i=1,1)/1/
      DATA (ab(15,i),i=1,1)/100.d0/
c
      DATA at(16),gel(16),nmn(16),(mn(16,i),i=1,4)/' S',5,4,32,33,34,36/
      DATA (zm(16,i),i=0,4)/32.066d0, 31.97207100d0, 32.97145876d0,
     1                       33.96786690d0, 35.96708076d0/
      DATA (ns2(16,i),i=1,4)/0,3,0,0/
      DATA (ab(16,i),i=1,4)/95.02d0, 0.75d0, 4.21d0, 0.02d0/
c
      DATA at(17),gel(17),nmn(17),(mn(17,i),i=1,2)/'Cl',4,2,35,37/
      DATA (zm(17,i),i=0,2)/35.4527d0, 34.96885268d0, 36.96590259d0/
      DATA (ns2(17,i),i=1,2)/3,3/
      DATA (ab(17,i),i=1,2)/75.77d0, 24.23d0/
c
      DATA at(18),gel(18),nmn(18),(mn(18,i),i=1,3)/'Ar',1,3,36,38,40/
      DATA (zm(18,i),i=0,3)/39.948d0, 35.967545106d0, 37.9627324d0,
     1                       39.9623831225d0/
      DATA (ns2(18,i),i=1,3)/0,0,0/
      DATA (ab(18,i),i=1,3)/0.337d0, 0.063d0, 99.600d0/
c
      DATA at(19),gel(19),nmn(19),(mn(19,i),i=1,3)/' K',2,3,39,40,41/
      DATA (zm(19,i),i=0,3)/39.0983d0, 38.96370668d0, 39.96399848d0,
     1                       40.96182576d0/
      DATA (ns2(19,i),i=1,3)/3,8,3/
      DATA (ab(19,i),i=1,3)/93.2581d0, 0.0117d0, 6.7302d0/
 
      DATA at(20),gel(20),nmn(20),(mn(20,i),i=1,6)/'Ca',1,6,40,42,43,44,
     1                                              46,48/
      DATA (zm(20,i),i=0,6)/40.078d0, 39.96259098d0, 41.95861801d0,
     1         42.9587666d0, 43.9554818d0, 45.9536926d0, 47.952534d0/
      DATA (ns2(20,i),i=1,6)/0,0,7,0,0,0/
      DATA (ab(20,i),i=1,6)/96.941d0, 0.647d0, 0.135d0, 2.086d0,
     1                      0.004d0, 0.187d0/
c
      DATA at(21),gel(21),nmn(21),(mn(21,i),i=1,1)/'Sc',4,1,45/
      DATA (zm(21,i),i=0,1)/44.955910d0, 44.9559119d0/
      DATA (ns2(21,i),i=1,1)/7/
      DATA (ab(21,i),i=1,1)/100.d0/
c
      DATA at(22),gel(22),nmn(22),(mn(22,i),i=1,5)/'Ti',5,5,46,47,48,49,
     1                                              50/
      DATA (zm(22,i),i=0,5)/47.88d0, 45.9526316d0, 46.9517631d0,
     1         47.9479463d0, 48.9478700d0, 49.9447912d0/
      DATA (ns2(22,i),i=1,5)/0,5,0,7,0/
      DATA (ab(22,i),i=1,5)/8.0d0, 7.3d0, 73.8d0, 5.5d0, 5.4d0/
c
      DATA at(23),gel(23),nmn(23),(mn(23,i),i=1,2)/' V',4,2,50,51/
      DATA (zm(23,i),i=0,2)/50.9415d0, 49.9471585d0, 50.9439595d0/
      DATA (ns2(23,i),i=1,2)/12,7/
      DATA (ab(23,i),i=1,2)/0.250d0, 99.750d0/
c
      DATA at(24),gel(24),nmn(24),(mn(24,i),i=1,4)/'Cr',7,4,50,52,53,54/
      DATA (zm(24,i),i=0,4)/51.9961d0, 49.9460442d0, 51.9405075d0,
     1                       52.9406494d0, 53.9388804d0/
      DATA (ns2(24,i),i=1,4)/0,0,3,0/
      DATA (ab(24,i),i=1,4)/4.345d0, 83.789d0, 9.501d0, 2.365d0/
c
      DATA at(25),gel(25),nmn(25),(mn(25,i),i=1,1)/'Mn',6,1,55/
      DATA (zm(25,i),i=0,1)/54.93805d0, 54.9380451d0/
      DATA (ns2(25,i),i=1,1)/5/
      DATA (ab(25,i),i=1,1)/100.d0/
c
      DATA at(26),gel(26),nmn(26),(mn(26,i),i=1,4)/'Fe',9,4,54,56,57,58/
      DATA (zm(26,i),i=0,4)/55.847d0, 53.9396105d0, 55.9349375d0,
     1                       56.9353940d0, 57.9332756d0/
      DATA (ns2(26,i),i=1,4)/0,0,1,0/
      DATA (ab(26,i),i=1,4)/5.8d0, 91.72d0, 2.2d0, 0.28d0/
c
      DATA at(27),gel(27),nmn(27),(mn(27,i),i=1,1)/'Co',10,1,59/
      DATA (zm(27,i),i=0,1)/58.93320d0, 58.9331950d0/
      DATA (ns2(27,i),i=1,1)/7/
      DATA (ab(27,i),i=1,1)/100.d0/
c
      DATA at(28),gel(28),nmn(28),(mn(28,i),i=1,5)/'Ni',9,5,58,60,61,62,
     1                                              64/
      DATA (zm(28,i),i=0,5)/58.69d0, 57.9353429d0, 59.9307864d0,
     1         60.9310560d0, 61.9283451d0, 63.9279660d0/
      DATA (ns2(28,i),i=1,5)/0,0,3,0,0/
      DATA (ab(28,i),i=1,5)/68.077d0,26.223d0,1.140d0,3.634d0,0.926d0/
c
      DATA at(29),gel(29),nmn(29),(mn(29,i),i=1,2)/'Cu',2,2,63,65/
      DATA (zm(29,i),i=0,2)/63.546d0, 62.9295975d0,64.9277895d0/
      DATA (ns2(29,i),i=1,2)/3,3/
      DATA (ab(29,i),i=1,2)/69.17d0, 30.83d0/
c
      DATA at(30),gel(30),nmn(30),(mn(30,i),i=1,5)/'Zn',1,5,64,66,67,68,
     1                                              70/
      DATA (zm(30,i),i=0,5)/65.40d0, 63.9291422d0, 65.9260334d0,
     1         66.9271273d0, 67.9248442d0, 69.9253193d0/
      DATA (ns2(30,i),i=1,5)/0,0,5,0,0/
      DATA (ab(30,i),i=1,5)/48.6d0, 27.9d0, 4.1d0, 18.8d0, 0.6d0/
c
      DATA at(31),gel(31),nmn(31),(mn(31,i),i=1,2)/'Ga',2,2,69,71/
Coxon   DATA (zm(31,i),i=0,2)/69.723d0, 68.925581d0, 70.9247073d0/
      DATA (zm(31,i),i=0,2)/69.723d0, 68.9255736d0, 70.9247013d0/
      DATA (ns2(31,i),i=1,2)/3,3/
      DATA (ab(31,i),i=1,2)/60.108d0, 39.892d0/
c
      DATA at(32),gel(32),nmn(32),(mn(32,i),i=1,5)/'Ge',1,5,70,72,73,74,
     1                                              76/
      DATA (zm(32,i),i=0,5)/72.61d0, 69.9242474d0, 71.9220758d0,
     1         72.9234589d0, 73.9211778d0, 75.9214026d0/
      DATA (ns2(32,i),i=1,5)/0,0,9,0,0/
      DATA (ab(32,i),i=1,5)/21.23d0, 27.66d0, 7.73d0, 35.94d0, 7.44d0/
c
      DATA at(33),gel(33),nmn(33),(mn(33,i),i=1,1)/'As',4,1,75/
      DATA (zm(33,i),i=0,1)/74.92159d0, 74.9215965d0/
      DATA (ns2(33,i),i=1,1)/3/
      DATA (ab(33,i),i=1,1)/100.d0/
c
      DATA at(34),gel(34),nmn(34),(mn(34,i),i=1,6)/'Se',5,6,74,76,77,78,
     1                                              80,82/
      DATA (zm(34,i),i=0,6)/78.96d0, 73.9224764d0, 75.9192136d0,
     1         76.9199140d0, 77.9173091d0, 79.9165213d0, 81.9166994d0/
      DATA (ns2(34,i),i=1,6)/0,0,1,0,0,0/
      DATA (ab(34,i),i=1,6)/0.89d0, 9.36d0, 7.63d0, 23.78d0, 49.61d0,
     1                      8.73d0/
c
      DATA at(35),gel(35),nmn(35),(mn(35,i),i=1,2)/'Br',4,2,79,81/
      DATA (zm(35,i),i=0,2)/79.904d0, 78.9183371d0, 80.9162906d0/
      DATA (ns2(35,i),i=1,2)/3,3/
      DATA (ab(35,i),i=1,2)/50.69d0, 49.31d0/
c
      DATA at(36),gel(36),nmn(36),(mn(36,i),i=1,6)/'Kr',1,6,78,80,82,83,
     1                                              84,86/
      DATA (zm(36,i),i=0,6)/83.80d0, 77.9203648d0, 79.9163790d0,
     1          81.9134836d0, 82.914136d0, 83.911507d0, 85.91061073d0/
      DATA (ns2(36,i),i=1,6)/0,0,0,9,0,0/
      DATA (ab(36,i),i=1,6)/0.35d0, 2.25d0, 11.6d0, 11.5d0, 57.0d0,
     1                      17.3d0/
c
      DATA at(37),gel(37),nmn(37),(mn(37,i),i=1,2)/'Rb',2,2,85,87/
      DATA (zm(37,i),i=0,2)/85.4678d0, 84.911789738d0, 86.909180527d0/
      DATA (ns2(37,i),i=1,2)/5,3/
      DATA (ab(37,i),i=1,2)/72.165d0, 27.835d0/
c
      DATA at(38),gel(38),nmn(38),(mn(38,i),i=1,4)/'Sr',1,4,84,86,87,88/
      DATA (zm(38,i),i=0,4)/87.62d0, 83.913425d0, 85.9092602d0,
     1                      86.9088771d0, 87.9056121d0/
      DATA (ns2(38,i),i=1,4)/0,0,9,0/
      DATA (ab(38,i),i=1,4)/0.56d0, 9.86d0, 7.00d0, 82.58d0/
c
      DATA at(39),gel(39),nmn(39),(mn(39,i),i=1,1)/' Y',4,1,89/
      DATA (zm(39,i),i=0,1)/88.90585d0, 88.9058483d0/
      DATA (ns2(39,i),i=1,1)/1/
      DATA (ab(39,i),i=1,1)/100.d0/
c
      DATA at(40),gel(40),nmn(40),(mn(40,i),i=1,5)/'Zr',5,5,90,91,92,94,
     1                                              96/
      DATA (zm(40,i),i=0,5)/91.224d0, 89.9047044d0, 90.9056458d0,
     1                      91.9050408d0, 93.9063152d0, 95.9082734d0/
      DATA (ns2(40,i),i=1,5)/0,5,0,0,0/
      DATA (ab(40,i),i=1,5)/51.45d0, 11.22d0, 17.15d0, 17.38d0, 2.80d0/
c
      DATA at(41),gel(41),nmn(41),(mn(41,i),i=1,1)/'Nb',2,1,93/
      DATA (zm(41,i),i=0,1)/92.90638d0, 92.9063781d0/
      DATA (ns2(41,i),i=1,1)/9/
      DATA (ab(41,i),i=1,1)/100.d0/
c
      DATA at(42),gel(42),nmn(42),(mn(42,i),i=1,7)/'Mo',7,7,92,94,95,96,
     1                                              97,98,100/
      DATA (zm(42,i),i=0,7)/95.94d0, 91.906811d0, 93.9050883d0,
     1        94.9058421d0, 95.9046795d0, 96.9060215d0, 97.9054082d0,
     2        99.907477d0/
      DATA (ns2(42,i),i=1,7)/0,0,5,0,5,0,0/
      DATA (ab(42,i),i=1,7)/14.84d0, 9.25d0, 15.92d0, 16.68d0, 9.55d0,
     1                      24.13d0, 9.63d0/
c
      DATA at(43),gel(43),nmn(43),(mn(43,i),i=1,1)/'Tc',6,1,98/
      DATA (zm(43,i),i=0,1)/97.907215d0, 97.907216d0/
      DATA (ns2(43,i),i=1,1)/12/
      DATA (ab(43,i),i=1,1)/100.d0/
c
      DATA at(44),gel(44),nmn(44),(mn(44,i),i=1,7)/'Ru',11,7,96,98,99,
     1                                              100,101,102,104/
      DATA (zm(44,i),i=0,7)/101.07d0, 95.907598d0, 97.905287d0,
     1     98.9059393d0, 99.9042195d0, 100.9055821d0, 101.9043493d0,
     2     103.905433d0/
      DATA (ns2(44,i),i=1,7)/0,0,5,0,5,0,0/
      DATA (ab(44,i),i=1,7)/5.52d0, 1.88d0, 12.7d0, 12.6d0, 17.0d0,
     1                      31.6d0, 18.7d0/
c
      DATA at(45),gel(45),nmn(45),(mn(45,i),i=1,1)/'Rh',10,1,103/
      DATA (zm(45,i),i=0,1)/102.90550d0, 102.905504d0/
      DATA (ns2(45,i),i=1,1)/1/
      DATA (ab(45,i),i=1,1)/100.d0/
c
      DATA at(46),gel(46),nmn(46),(mn(46,i),i=1,6)/'Pd',1,6,102,104,105,
     1                                              106,108,110/
      DATA (zm(46,i),i=0,6)/106.42d0, 101.905609d0, 103.904036d0,
     1       104.905085d0, 105.903486d0, 107.903892d0, 109.905153d0/
      DATA (ns2(46,i),i=1,6)/0,0,5,0,0,0/
      DATA (ab(46,i),i=1,6)/1.02d0, 11.14d0, 22.33d0, 27.33d0, 26.46d0,
     1                      11.72d0/
c
      DATA at(47),gel(47),nmn(47),(mn(47,i),i=1,2)/'Ag',2,2,107,109/
      DATA (zm(47,i),i=0,2)/107.8682d0, 106.905097d0, 108.904752d0/
      DATA (ns2(47,i),i=1,2)/1,1/
      DATA (ab(47,i),i=1,2)/51.839d0, 48.161d0/
c
      DATA at(48),gel(48),nmn(48),(mn(48,i),i=1,8)/'Cd',1,8,106,108,110,
     1                                             111,112,113,114,116/ 
      DATA (zm(48,i),i=0,8)/112.411d0, 105.906459d0, 107.904184d0, 
     1       109.9030021d0, 110.9041781d0, 111.9027578d0, 112.9044017d0,
     2       113.9033585d0, 115.904756d0/
      DATA (ns2(48,i),i=1,8)/0,0,0,1,0,1,0,0/
      DATA (ab(48,i),i=1,8)/1.25d0, 0.89d0, 12.49d0, 12.80d0, 24.13d0,
     1                      12.22d0, 28.73d0, 7.49d0/
c
      DATA at(49),gel(49),nmn(49),(mn(49,i),i=1,2)/'In',2,2,113,115/
      DATA (zm(49,i),i=0,2)/114.818d0, 112.904058d0, 114.903878d0/
      DATA  (ns2(49,i),i=1,2)/9,9/
      DATA (ab(49,i),i=1,2)/4.3d0, 95.7d0/
c
      DATA at(50),gel(50),nmn(50),(mn(50,i),i=1,10)/'Sn',1,10,112,114,
     1                                 115,116,117,118,119,120,122,124/
      DATA (zm(50,i),i=0,10)/118.710d0, 111.904818d0, 113.902779d0,
     1     114.903342d0, 115.901741d0, 116.902952d0, 117.901603d0,
     2     118.903308d0, 119.9021947d0, 121.9034390d0, 123.9052739d0/
      DATA (ns2(50,i),i=1,10)/0,0,1,0,1,0,1,0,0,0/
      DATA (ab(50,i),i=1,10)/0.97d0, 0.65d0, 0.34d0, 14.53d0, 7.68d0,
     1                       24.23d0, 8.59d0, 32.59d0, 4.63d0, 5.79d0/
c
      DATA at(51),gel(51),nmn(51),(mn(51,i),i=1,2)/'Sb',4,2,121,123/
      DATA (zm(51,i),i=0,2)/121.757d0, 120.9038157d0, 122.9042140d0/
      DATA (ns2(51,i),i=1,2)/5,7/
      DATA (ab(51,i),i=1,2)/57.36d0, 42.64d0/
c
      DATA at(52),gel(52),nmn(52),(mn(52,i),i=1,8)/'Te',5,8,120,122,123,
     1                                             124,125,126,128,130/
      DATA (zm(52,i),i=0,8)/127.60d0, 119.904020d0, 121.9030439d0,
     1    122.9042700d0, 123.9028179d0, 124.9044307d0, 125.9033117d0,
     2    127.9044631d0, 129.9062244d0/
      DATA (ns2(52,i),i=1,8)/0,0,1,0,1,0,0,0/
      DATA (ab(52,i),i=1,8)/0.096d0, 2.603d0, 0.908d0, 4.816d0,
     1                      7.139d0, 18.95d0, 31.69d0, 33.80d0/
c
      DATA at(53),gel(53),nmn(53),(mn(53,i),i=1,2)/' I',4,2,127,129/
      DATA (zm(53,i),i=0,2)/126.90447d0, 126.904473d0, 128.904988d0/
      DATA (ns2(53,i),i=1,2)/5,7/
      DATA (ab(53,i),i=1,2)/100.d0,0.d0/
c
      DATA at(54),gel(54),nmn(54),(mn(54,i),i=1,9)/'Xe',1,9,124,126,128,
     1                                          129,130,131,132,134,136/
      DATA (zm(54,i),i=0,9)/131.29d0, 123.9058930d0, 125.904274d0,
     1    127.9035313d0, 128.9047794d0, 129.9035080d0, 130.9050824d0,
     2    131.9041535d0, 133.9053945d0, 135.907219d0/
      DATA (ns2(54,i),i=1,9)/0,0,0,1,0,3,0,0,0/
      DATA (ab(54,i),i=1,9)/0.10d0, 0.09d0, 1.91d0, 26.4d0, 4.1d0,
     1                      21.2d0, 26.9d0, 10.4d0, 8.9d0/
c
      DATA at(55),gel(55),nmn(55),(mn(55,i),i=1,1)/'Cs',2,1,133/
      DATA (zm(55,i),i=0,1)/132.90543d0, 132.905451933d0/
      DATA (ns2(55,i),i=1,1)/7/
      DATA (ab(55,i),i=1,1)/100.d0/
c
      DATA at(56),gel(56),nmn(56),(mn(56,i),i=1,7)/'Ba',1,7,130,132,134,
     1                                             135,136,137,138/
      DATA (zm(56,i),i=0,7)/137.327d0, 129.9063208d0, 131.9050613d0,
     1    133.9045084d0, 134.9056886d0, 135.9045759d0, 136.9058274d0,
     2    137.9052472d0/
      DATA (ns2(56,i),i=1,7)/0,0,0,3,0,3,0/
      DATA (ab(56,i),i=1,7)/0.106d0, 0.101d0, 2.417d0, 6.592d0, 
     1                      7.854d0, 11.23d0, 71.70d0/
c
      DATA at(57),gel(57),nmn(57),(mn(57,i),i=1,2)/'La',4,2,138,139/
      DATA (zm(57,i),i=0,2)/138.9055d0, 137.907112d0, 138.9063533d0/
      DATA (ns2(57,i),i=1,2)/10,7/ 
      DATA (ab(57,i),i=1,2)/0.0902d0, 99.9098d0/
c
      DATA at(58),gel(58),nmn(58),(mn(58,i),i=1,4)/'Ce',9,4,136,138,140,
     1                                             142/
      DATA (zm(58,i),i=0,4)/140.115d0, 135.907172d0, 137.905991d0,
     1    139.9054387d0, 141.909244d0/
      DATA (ns2(58,i),i=1,4)/0,0,0,0/
      DATA (ab(58,i),i=1,4)/0.19d0, 0.25d0, 88.48d0, 11.08d0/
c
      DATA at(59),gel(59),nmn(59),(mn(59,i),i=1,1)/'Pr',10,1,141/
      DATA (zm(59,i),i=0,1)/140.90765d0, 140.9076528d0/
      DATA (ns2(59,i),i=1,1)/5/
      DATA (ab(59,i),i=1,1)/100.d0/
c
      DATA at(60),gel(60),nmn(60),(mn(60,i),i=1,7)/'Nd',9,7,142,143,144,
     1                                             145,146,148,150/
      DATA (zm(60,i),i=0,7)/144.24d0, 141.9077233d0, 142.9098143d0,
     1    143.9100873d0, 144.9125736d0, 145.9131169d0, 147.916893d0,
     2    149.920891d0/
      DATA (ns2(60,i),i=1,7)/0,7,0,7,0,0,0/
      DATA (ab(60,i),i=1,7)/27.13d0, 12.18d0, 23.80d0, 8.30d0, 17.19d0,
     1                       5.76d0, 5.64d0/
c
      DATA at(61),gel(61),nmn(61),(mn(61,i),i=1,1)/'Pm',6,1,145/
      DATA (zm(61,i),i=0,1)/144.912743d0, 144.912749d0/
      DATA (ns2(61,i),i=1,1)/5/
      DATA (ab(61,i),i=1,1)/100.d0/
c
      DATA at(62),gel(62),nmn(62),(mn(62,i),i=1,7)/'Sm',1,7,144,147,148,
     1                                             149,150,152,154/
      DATA (zm(62,i),i=0,7)/150.36d0, 143.911999d0, 146.9148979d0,
     1    147.9148227d0, 148.9171847d0, 149.9172755d0, 151.9197324d0,
     2    153.9222093d0/
      DATA (ns2(62,i),i=1,7)/0,7,0,7,0,0,0/
      DATA (ab(62,i),i=1,7)/3.1d0, 15.0d0, 11.3d0, 13.8d0, 7.4d0,
     1                      26.7d0, 22.7d0/
c
      DATA at(63),gel(63),nmn(63),(mn(63,i),i=1,2)/'Eu',8,2,151,153/
      DATA (zm(63,i),i=0,2)/151.965d0, 150.9198502d0, 152.9212303d0/
      DATA (ns2(63,i),i=1,2)/5,5/
      DATA (ab(63,i),i=1,2)/47.8d0, 52.2d0/
c
      DATA at(64),gel(64),nmn(64),(mn(64,i),i=1,7)/'Gd',5,7,152,154,155,
     1                                              156,157,158,160/
      DATA (zm(64,i),i=0,7)/157.25d0, 151.9197910d0, 153.92086560,
     1    154.9226220d0, 155.9221227d0, 156.9239601d0, 157.9241039d0,
     2    159.9270541d0/
      DATA (ns2(64,i),i=1,7)/0,0,3,0,3,0,0/
      DATA (ab(64,i),i=1,7)/0.20d0, 2.18d0, 14.80d0, 20.47d0, 15.65d0,
     1                      24.84d0, 21.86d0/
c
      DATA at(65),gel(65),nmn(65),(mn(65,i),i=1,1)/'Tb',16,1,159/
      DATA (zm(65,i),i=0,1)/158.92534d0, 158.9253468d0/
      DATA (ns2(65,i),i=1,1)/3/
      DATA (ab(65,i),i=1,1)/100.d0/
c
      DATA at(66),gel(66),nmn(66),(mn(66,i),i=1,7)/'Dy',17,7,156,158,
     1                                           160,161,162,163,164/
      DATA (zm(66,i),i=0,7)/162.50d0, 155.924283d0, 157.924409d0,
     1    159.9251975d0, 160.9269334d0, 161.9267984d0, 162.9287312d0,
     2    163.9291748d0/
      DATA (ns2(66,i),i=1,7)/0,0,0,5,0,5,0/
      DATA (ab(66,i),i=1,7)/0.06d0, 0.10d0, 2.34d0, 18.9d0, 25.5d0,
     1                      24.9d0, 28.2d0/
c
      DATA at(67),gel(67),nmn(67),(mn(67,i),i=1,1)/'Ho',16,1,165/
      DATA (zm(67,i),i=0,1)/164.93032d0, 164.9303221d0/
      DATA (ns2(67,i),i=1,1)/7/
      DATA (ab(67,i),i=1,1)/100.d0/
     
      DATA at(68),gel(68),nmn(68),(mn(68,i),i=1,6)/'Er',13,6,162,164,
     1                                            166,167,168,170/
      DATA (zm(68,i),i=0,6)/167.26d0, 161.928778d0, 163.929200d0,
     1    165.9302931d0, 166.9320482d0, 167.9323702d0, 169.9354643d0/
      DATA (ns2(68,i),i=1,6)/0,0,0,7,0,0/
      DATA (ab(68,i),i=1,6)/0.14d0, 1.61d0, 33.6d0, 22.95d0, 26.8d0,
     1                      14.9d0/
c
      DATA at(69),gel(69),nmn(69),(mn(69,i),i=1,1)/'Tm',8,1,169/  
      DATA (zm(69,i),i=0,1)/168.93421d0, 168.9342133d0/
      DATA (ns2(69,i),i=1,1)/1/
      DATA (ab(69,i),i=1,1)/100.d0/
c
      DATA at(70),gel(70),nmn(70),(mn(70,i),i=1,7)/'Yb',1,7,168,170,171,
     1                                            172,173,174,176/
      DATA (zm(70,i),i=0,7)/173.04d0, 167.933897d0, 169.9347618d0,
     1    170.936323580, 171.9363815d0, 172.9382108d0, 173.9388621d0,
     2    175.9425717d0/
      DATA (ns2(70,i),i=1,7)/0,0,1,0,5,0,0/
      DATA (ab(70,i),i=1,7)/0.13d0, 3.05d0, 14.3d0, 21.9d0, 16.12d0,
     1                      31.8d0, 12.7d0/
c
      DATA at(71),gel(71),nmn(71),(mn(71,i),i=1,2)/'Lu',4,2,175,176/
      DATA (zm(71,i),i=0,2)/174.967d0, 174.9407718d0, 175.9426863d0/
      DATA (ns2(71,i),i=1,2)/7,14/
      DATA (ab(71,i),i=1,2)/97.41d0, 2.59d0/
c
      DATA at(72),gel(72),nmn(72),(mn(72,i),i=1,6)/'Hf',5,6,174,176,177,
     1                                             178,179,180/
      DATA (zm(72,i),i=0,6)/178.49d0, 173.940046d0, 175.9414086d0,
     1    176.9432207d0, 177.9436988d0, 178.9458161d0, 179.9465500d0/
      DATA (ns2(72,i),i=1,6)/0,0,7,0,9,0/
      DATA (ab(72,i),i=1,6)/0.162d0, 5.206d0, 18.606d0, 27.297d0,
     1                      13.629d0, 35.100d0/
c
      DATA at(73),gel(73),nmn(73),(mn(73,i),i=1,2)/'Ta',4,2,180,181/
      DATA (zm(73,i),i=0,2)/180.9479d0, 179.9474648d0, 180.9479958d0/
      DATA (ns2(73,i),i=1,2)/16,7/
      DATA (ab(73,i),i=1,2)/0.012d0, 99.988d0/
c
      DATA at(74),gel(74),nmn(74),(mn(74,i),i=1,5)/' W',1,5,180,182,183,
     1                                             184,186/
      DATA (zm(74,i),i=0,5)/183.84d0, 179.946704d0, 181.9482042d0,
     1    182.9502230d0, 183.9509312d0, 185.9543641d0/
      DATA (ns2(74,i),i=1,5)/0,0,1,0,0/
      DATA (ab(74,i),i=1,5)/0.13d0, 26.3d0, 14.3d0, 30.67d0, 28.6d0/
c
      DATA at(75),gel(75),nmn(75),(mn(75,i),i=1,2)/'Re',6,2,185,187/
      DATA (zm(75,i),i=0,2)/186.207d0, 184.9529550d0, 186.9557531d0/
      DATA (ns2(75,i),i=1,2)/5,5/
      DATA (ab(75,i),i=1,2)/37.40d0, 62.60d0/
c
      DATA at(76),gel(76),nmn(76),(mn(76,i),i=1,7)/'Os',9,7,184,186,187,
     1                                             188,189,190,192/
      DATA (zm(76,i),i=0,7)/190.23d0, 183.9524891d0, 185.9538382d0,
     1    186.9557505d0, 187.9558382d0, 188.9581475d0, 189.9584470d0,
     2    191.9614807d0/
      DATA (ns2(76,i),i=1,7)/0,0,1,0,3,0,0/
      DATA (ab(76,i),i=1,7)/0.02d0, 1.58d0, 1.6d0, 13.3d0, 16.1d0,
     1                      26.4d0, 41.0d0/
c
      DATA at(77),gel(77),nmn(77),(mn(77,i),i=1,2)/'Ir',10,2,191,193/
      DATA (zm(77,i),i=0,2)/192.22d0, 190.9605940d0, 192.9629264d0/
      DATA (ns2(77,i),i=1,2)/3,3/
      DATA (ab(77,i),i=1,2)/37.3d0, 62.7d0/
c
c
      DATA at(78),gel(78),nmn(78),(mn(78,i),i=1,6)/'Pt',7,6,190,192,194,
     1                                            195,196,198/
      DATA (zm(78,i),i=0,6)/195.08d0, 189.959932d0, 191.9610380d0,
     1    193.9626803d0, 194.9647911d0, 195.9649515d0, 197.967893d0/
      DATA (ns2(78,i),i=1,6)/0,0,0,1,0,0/
      DATA (ab(78,i),i=1,6)/0.01d0,0.79d0,32.9d0,33.8d0,25.3d0,7.2d0/
c
      DATA at(79),gel(79),nmn(79),(mn(79,i),i=1,1)/'Au',2,1,197/
      DATA (zm(79,i),i=0,1)/196.96654d0, 196.9665687d0/
      DATA (ns2(79,i),i=1,1)/3/
      DATA (ab(79,i),i=1,1)/100.d0/
c
      DATA at(80),gel(80),nmn(80),(mn(80,i),i=1,7)/'Hg',1,7,196,198,199,
     1                                            200,201,202,204/
      DATA (zm(80,i),i=0,7)/200.59d0, 195.965833d0, 197.9667690d0,
     1    198.9682799d0, 199.9683260d0, 200.9703023d0, 201.9706430d0,
     2    203.9734939d0/
      DATA (ns2(80,i),i=1,7)/0,0,1,0,3,0,0/
      DATA (ab(80,i),i=1,7)/0.15d0, 9.97d0, 16.87d0, 23.10d0, 13.18d0,
     1                      29.86d0, 6.87d0/
c
      DATA at(81),gel(81),nmn(81),(mn(81,i),i=1,2)/'Tl',2,2,203,205/
      DATA (zm(81,i),i=0,2)/204.3833d0, 202.9723442d0, 204.9744275d0/
      DATA (ns2(81,i),i=1,2)/1,1/
      DATA (ab(81,i),i=1,2)/29.524d0, 70.476d0/
c
      DATA at(82),gel(82),nmn(82),(mn(82,i),i=1,4)/'Pb',1,4,204,206,207,
     1                                             208/
      DATA (zm(82,i),i=0,4)/207.2d0, 203.9730436d0, 205.9744653d0,
     1    206.9758969d0, 207.9766521d0/
      DATA (ns2(82,i),i=1,4)/0,0,1,0/
      DATA (ab(82,i),i=1,4)/1.4d0, 24.1d0, 22.1d0, 52.4d0/
c
      DATA at(83),gel(83),nmn(83),(mn(83,i),i=1,1)/'Bi',4,1,209/
      DATA (zm(83,i),i=0,1)/208.98037d0, 208.9803987d0/
      DATA (ns2(83,i),i=1,1)/9/
      DATA (ab(83,i),i=1,1)/100.d0/
c
      DATA at(84),gel(84),nmn(84),(mn(84,i),i=1,1)/'Po',5,1,209/
      DATA (zm(84,i),i=0,1)/208.982404d0, 208.9824304d0/
      DATA (ns2(84,i),i=1,1)/1/
      DATA (ab(84,i),i=1,1)/100.d0/
c
      DATA at(85),gel(85),nmn(85),(mn(85,i),i=1,1)/'At',-1,1,210/
      DATA (zm(85,i),i=0,1)/209.987126d0, 209.987148d0/
      DATA (ns2(85,i),i=1,1)/10/
      DATA (ab(85,i),i=1,1)/100.d0/
c
      DATA at(86),gel(86),nmn(86),(mn(86,i),i=1,1)/'Rn',1,1,222/
      DATA (zm(86,i),i=0,1)/222.017571d0, 222.0175777d0/
      DATA (ns2(86,i),i=1,1)/0/
      DATA (ab(86,i),i=1,1)/100.d0/
c
      DATA at(87),gel(87),nmn(87),(mn(87,i),i=1,1)/'Fr',-1,1,223/
      DATA (zm(87,i),i=0,1)/223.019733d0, 223.0197359d0/
      DATA (ns2(87,i),i=1,1)/3/
      DATA (ab(87,i),i=1,1)/100.d0/
c
      DATA at(88),gel(88),nmn(88),(mn(88,i),i=1,1)/'Ra',1,1,226/
      DATA (zm(88,i),i=0,1)/226.025403d0, 226.0254098d0/
      DATA (ns2(88,i),i=1,1)/0/
      DATA (ab(88,i),i=1,1)/100.d0/
c
      DATA at(89),gel(89),nmn(89),(mn(89,i),i=1,1)/'Ac',4,1,227/
      DATA (zm(89,i),i=0,1)/227.027750d0, 227.0277521d0/
      DATA (ns2(89,i),i=1,1)/3/
      DATA (ab(89,i),i=1,1)/100.d0/
c
      DATA at(90),gel(90),nmn(90),(mn(90,i),i=1,1)/'Th',-1,1,232/
      DATA (zm(90,i),i=0,1)/232.038d0, 232.0380553d0/
      DATA (ns2(90,i),i=1,1)/0/
      DATA (ab(90,i),i=1,1)/100.d0/
c
      DATA at(91),gel(91),nmn(91),(mn(91,i),i=1,1)/'Pa',-1,1,231/
      DATA (zm(91,i),i=0,1)/231.03588d0, 231.0358840d0/
      DATA (ns2(91,i),i=1,1)/3/
      DATA (ab(91,i),i=1,1)/100.d0/
c
      DATA at(92),gel(92),nmn(92),(mn(92,i),i=1,4)/' U',-1,4,233,234,
     1                                             235,238/
      DATA (zm(92,i),i=0,4)/238.0289d0, 233.0396352d0, 234.0409521d0,
     1    235.0439299d0, 238.0507882d0/
      DATA (ns2(92,i),i=1,4)/5,0,7,0/
      DATA (ab(92,i),i=1,4)/0.d0, 0.0055d0, 0.7200d0, 99.2745d0/
c
      DATA at(93),gel(93),nmn(93),(mn(93,i),i=1,1)/'Np',-1,1,237/
      DATA (zm(93,i),i=0,1)/237.0481678d0, 237.0481734d0/
      DATA (ns2(93,i),i=1,1)/5/
      DATA (ab(93,i),i=1,1)/100.d0/
c
      DATA at(94),gel(94),nmn(94),(mn(94,i),i=1,1)/'Pu',-1,1,244/
      DATA (zm(94,i),i=0,1)/244.064199d0, 244.064204d0/
      DATA (ns2(94,i),i=1,1)/0/
      DATA (ab(94,i),i=1,1)/100.d0/
c
      DATA at(95),gel(95),nmn(95),(mn(95,i),i=1,1)/'Am',-1,1,243/
      DATA (zm(95,i),i=0,1)/243.061375d0, 243.0613811d0/
      DATA (ns2(95,i),i=1,1)/5/
      DATA (ab(95,i),i=1,1)/100.d0/
c
      DATA at(96),gel(96),nmn(96),(mn(96,i),i=1,1)/'Cm',-1,1,247/
      DATA (zm(96,i),i=0,1)/247.070347d0, 247.070354d0/
      DATA (ns2(96,i),i=1,1)/9/
      DATA (ab(96,i),i=1,1)/100.d0/
c
      DATA at(97),gel(97),nmn(97),(mn(97,i),i=1,1)/'Bk',-1,1,247/
      DATA (zm(97,i),i=0,1)/247.070300d0, 247.070307d0/
      DATA (ns2(97,i),i=1,1)/3/
      DATA (ab(97,i),i=1,1)/100.d0/
c
      DATA at(98),gel(98),nmn(98),(mn(98,i),i=1,1)/'Cf',-1,1,251/
      DATA (zm(98,i),i=0,1)/251.079580d0, 251.079587d0/
      DATA (ns2(98,i),i=1,1)/1/
      DATA (ab(98,i),i=1,1)/100.d0/
c
      DATA at(99),gel(99),nmn(99),(mn(99,i),i=1,1)/'Es',-1,1,252/
      DATA (zm(99,i),i=0,1)/252.082944d0, 252.082980d0/
      DATA (ns2(99,i),i=1,1)/10/
      DATA (ab(99,i),i=1,1)/100.d0/
c
      DATA at(100),gel(100),nmn(100),(mn(100,i),i=1,1)/'Fm',-1,1,257/
      DATA (zm(100,i),i=0,1)/257.095099d0, 257.095105d0/
      DATA (ns2(100,i),i=1,1)/9/
      DATA (ab(100,i),i=1,1)/100.d0/
c
      DATA at(101),gel(101),nmn(101),(mn(101,i),i=1,1)/'Md',-1,1,258/
      DATA (zm(101,i),i=0,1)/258.09857d0, 258.098431d0/
      DATA (ns2(101,i),i=1,1)/16/
      DATA (ab(101,i),i=1,1)/100.d0/
c
      DATA at(102),gel(102),nmn(102),(mn(102,i),i=1,1)/'No',-1,1,259/
      DATA (zm(102,i),i=0,1)/259.100931d0, 259.101030d0/
      DATA (ns2(102,i),i=1,1)/9/
      DATA (ab(102,i),i=1,1)/100.d0/
c
      DATA at(103),gel(103),nmn(103),(mn(103,i),i=1,1)/'Lr',-1,1,260/
      DATA (zm(103,i),i=0,1)/260.105320d0, 260.105500d0/
      DATA (ns2(103,i),i=1,1)/-1/
      DATA (ab(103,i),i=1,1)/100.d0/
c
      DATA at(104),gel(104),nmn(104),(mn(104,i),i=1,1)/'Rf',-1,1,261/
      DATA (zm(104,i),i=0,1)/261.10869d0, 261.108770d0/
      DATA (ns2(104,i),i=1,1)/-1/
      DATA (ab(104,i),i=1,1)/100.d0/
c
      DATA at(105),gel(105),nmn(105),(mn(105,i),i=1,1)/'Db',-1,1,262/
      DATA (zm(105,i),i=0,1)/262.11376d0, 262.114080d0/
      DATA (ns2(105,i),i=1,1)/-1/
      DATA (ab(105,i),i=1,1)/100.d0/
c
      DATA at(106),gel(106),nmn(106),(mn(106,i),i=1,1)/'Sg',-1,1,263/
      DATA (zm(106,i),i=0,1)/263.11822d0, 263.118320d0/
      DATA (ns2(106,i),i=1,1)/-1/
      DATA (ab(106,i),i=1,1)/100.d0/
c
      DATA at(107),gel(107),nmn(107),(mn(107,i),i=1,1)/'Bh',-1,1,262/
      DATA (zm(107,i),i=0,1)/262.12293d0, 262.122890d0/
      DATA (ns2(107,i),i=1,1)/-1/
      DATA (ab(107,i),i=1,1)/100.d0/
c
      DATA at(108),gel(108),nmn(108),(mn(108,i),i=1,1)/'Hs',-1,1,265/
      DATA (zm(108,i),i=0,1)/265.13016d0, 265.130090d0/
      DATA (ns2(108,i),i=1,1)/-1/
      DATA (ab(108,i),i=1,1)/100.d0/
c
      DATA at(109),gel(109),nmn(109),(mn(109,i),i=1,1)/'Mt',-1,1,266/
      DATA (zm(109,i),i=0,1)/266.13764d0, 266.137300d0/
      DATA (ns2(109,i),i=1,1)/-1/
      DATA (ab(109,i),i=1,1)/100.d0/
c
      IF((IAN.LE.0).OR.(IAN.GT.109)) THEN
          MASS= 0.d0
          NAME= 'XX'
          IMN= 0
          WRITE(6,601) IAN
          RETURN
        ELSE
          NAME= AT(IAN)
        ENDIF
      IF((IAN.EQ.1).AND.(IMN.NE.1)) THEN
c** Special case: insert common name for deuterium or tritium
          IF(IMN.EQ.2) NAME=' D'
          IF(IMN.EQ.3) NAME=' T'
          ENDIF
      GELGS= GEL(IAN)
      MASS= -1.d0
      GNS= -1
	ABUND = -1.d0
      DO  I= 1,NMN(IAN)
          if(i.gt.10)  write(6,606) ian,imn,nmn(ian)
  606  format(3i9)
          IF(IMN.EQ.MN(IAN,I)) THEN
              MASS= ZM(IAN,I)
              GNS= NS2(IAN,I)+1
              ABUND = AB(IAN,I)
              ENDIF
          ENDDO
      IF(MASS.LT.0.d0) THEN
          MASS= ZM(IAN,0)
          IF(IMN.NE.0) WRITE(6,602) AT(IAN),IMN
          IMN= 0
          ENDIF
      RETURN
  601 FORMAT(' *** MASSES Data base does not include Atomic Number=',i4)
  602 FORMAT(' *** MASSES Data base does not include ',A2,'(',i3,
     1 '), so use average atomic mass.')
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE READATA(NSTATES,PASok,UCUTOFF,JTRUNC,EFSEL,VMIN,VMAX,
     1                                            NDAT,NOWIDTHS,PRINP)
c***********************************************************************
c** Subroutine to read, do book-keeping for, and print summary of
c  experimental data used in fits to spectroscopic data for one or more
c  electronic states and one or more isotopomers. 
c             ********* Version of 6 June 2006 *********
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++  COPYRIGHT 1997-2006 by  Robert J. Le Roy & Dominique R.T. Appadoo +
c   Dept. of Chemistry, Univ. of Waterloo, Waterloo, Ontario, Canada   +
c    This software may not be sold or any other commercial use made    +
c      of it without the express written permission of the authors.    +
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** The present program version can treat seven types of experimental
c   experimental data, for up to NISTPMX isotopomers of a given species.
c   The data are read in grouped as "bands", as (fluorescence) series, 
c   as binding energies (from photoassociation spectroscopy), as a set
c   of Bv values for a given electronic state, and [in a potential-fit
c   aanalysis] as tunneling predissociation level widths.  The types are
c   identified by the values of the 'electronic state label' parameters
c   IEP & IEPP.  They are:
c (i)  microwave transitions within a given electronic state;
c (ii)  infrared bands among the vibrational levels a given state;
c (iii) fluorescence series from some initial excited state level into 
c    vibration-rotation levels of a given electronic state
c (iv)  visible (electronic) absorption or emission bands between vib.
c    levels of two electronic state.
c (v)  binding energies - as from photoassociation spectroscopy
c (vi) "experimental" B_v values for vibrational levels of one of the
c    electronic states.
c (vii) Widths of tunneling predissociation quasibound levels (this 
c    option only meaningful for program DSPotFit).  
c-----------------------------------------------------------------------
c** On Entry:
c  NSTATES is the number of electronic states involved in the data set
c    considered (don't count states giving rise to fluorescence series).
c  PASok indicates how photoassociation data to be treated in analysis:
c    If(PASok(ISTATE).GE.1) treat it as proper PA binding energy data.
c    If(PASok(ISTATE).LE.0) treat PAS data as fluorescence series.
c    Set PASok= 0 if potential model has no explicit Dissoc. Energy
c  Data cutoffs:  for levels of electronic state  s , neglect data with:
c     J(s) > JTRUNC(s),  or vibrational levels lying outside the range
c     VMIN(s)  to  VMAX(s),  AND  NEGLECT any data for which the read-
c     in uncertainty is  > UCUTOFF (cm-1).  EFSEL(s) > 0 causes f-parity
c     levels to be neglected, EFSEL(s) < 0 omits e-parity levels
c     while  EFSEL(s) = 0  allows both types of parity to be included.
c  NOWIDTHS > 0  causes the program to ignore any tunneling widths in
c            the data set.
c  PRINP > 0  turns on the printing of a summary description of the data.
c** On Return:
c  UCUTOFF (cm-1)  is the smallest uncertainty in the (accepted) data
c  NDAT(v,i,s)  is the number of transitions associated with 
c    vibrational level-v of isotopomer-i of state-s [for NDEGB < 0 case]
c** This subroutine reads in the experimental data on channel-4
c-----------------------------------------------------------------------
      INCLUDE 'arrsizes.h'
c
      INTEGER I,IBB,NTRANS,COUNT,IBAND,JMAX(NPARMX),JMIN(NPARMX),
     1  VMX(NSTATEMX),ISOT,NBND,ESP,ESPP,ISTATE,ISTATEE,MN1,MN2,PRINP,
     2  FSOMIT,VMAXesp,VMINesp,VMAXespp,VMINespp,JTRUNCesp,JTRUNCespp
      INTEGER NSTATES,NOWIDTHS,JTRUNC(NSTATEMX),EFSEL(NSTATEMX),
     1  VMIN(NSTATEMX),VMAX(NSTATEMX),NDAT(0:NVIBMX,NISTPMX,NSTATEMX),
     2  PASok(NSTATES)
      REAL*8 UCUTOFF,UMIN,TOTUFREQ
      CHARACTER*3 NEF(-1:1)
      CHARACTER*2 LABLP,LABLPP

      INCLUDE 'BLKISOT.h'
      INCLUDE 'BLKDATA.h'
      INCLUDE 'BLKTYPE.h'
c
      DATA NEF/'  f','   ','  e'/
c-----------------------------------------------------------------------
      WRITE(6,603) UCUTOFF 
      DO  ISTATE= 1,NSTATES
          IF(JTRUNC(ISTATE).GE.0) THEN
              WRITE(6,607) SLABL(ISTATE),JTRUNC(ISTATE),VMIN(ISTATE),
     1                                                    VMAX(ISTATE)
            ELSE
              WRITE(6,605) SLABL(ISTATE),-JTRUNC(ISTATE),VMIN(ISTATE),
     1                                                    VMAX(ISTATE)
            ENDIF 
          IF(EFSEL(ISTATE).GT.0) WRITE(6,601) NEF(-1)
          IF(EFSEL(ISTATE).LT.0) WRITE(6,601) NEF(1)
          ENDDO
      UMIN= UCUTOFF
c** Initialize counters for book-keeping on input data
      COUNT= 0
      DO  ISOT= 1,NISTP
          DO  ISTATE= 1,NSTATES
              NTRANSFS(ISOT,ISTATE)= 0
              NTRANSIR(ISOT,ISTATE)= 0
              NTRANSMW(ISOT,ISTATE)= 0
              NBANDFS(ISOT,ISTATE)= 0
              NBANDVIS(ISOT,ISTATE)= 0
              NBANDIR(ISOT,ISTATE)= 0
              NBANDMW(ISOT,ISTATE)= 0
              NBVPP(ISOT,ISTATE)= 0
              NWIDTH(ISOT,ISTATE)= 0
              NEBPAS(ISOT,ISTATE)= 0
              DO  I= 1,NSTATES
                  NTRANSVIS(ISOT,ISTATE,I)= 0
                  NBANDEL(ISOT,ISTATE,I)= 0
                  ENDDO
              ENDDO
          NBANDS(ISOT)= 0
          ENDDO
      DO  ISTATE= 1,NSTATES
          VMX(ISTATE)= 0
          ENDDO
      NFSTOT= 0
      FSOMIT= 0
c========================================================================
c** Begin loop to read in data, band(or series)-by-band(or series).
c  STOP when run out of bands or when encounter a negative vibrational
c  quantum number.
c** Read all data for each isotopomer at one time.
      IBAND= 0
   10 CONTINUE
      IBAND= IBAND+1
      IF(IBAND.GT.NPARMX) THEN
            IF(PRINP.GT.0) WRITE(6,609) IBAND,NPARMX
            IBAND= IBAND-1
            GOTO 20
            ENDIF
c
c For each "band", read in:  (i) upper/lower vibrational quantum numbers
c   VP & VPP,  (ii) a two-character electronic-state alphameric label 
c   {enclosed in single quotes; e.g., 'X0' or 'A1'} for the upper
c   (LABLP) and lower (LABLP) state, and  (iii) integers NM1 & NM2 are
c   the mass numbers [corresponding to input atomic numbers AN(1) & 
c   AN(2)] identifying the particular isotopomer.  Note that LABLP also
c   identifies the type of data in the 'band' or data-group (see below).
c
c** LABLP = LABLPP  and  VP = VPP  for a microwave band
c   LABLP = LABLPP  and  VP.ne.VPP  for an infrared band 
c   LABLP = 'FS'  identifies this data group/band as a fluorescence 
c           series from a single emitting level into vibrational levels
c           of electronic state LABLPP.  In this case: VP is the quantum
c           number v' for the emitting level, while VPP is actually the 
c           rotational quantum number J' for the emitting level and JP
c           [see below] the lower state vibrational quantum number v".
c   LABLP = 'PA'  identifies this data group/band as a set of binding
c           energies [D-E(v,J,p)] for a given state.  Labels as for 'FS'
c   LABLP = 'BV'  identifies this data group/band as a set of Bv values
c           for electronic state LABLPP.  In this case, parameters  VP
c           & VPP are dummy variables, as are EFP, JPP and EFPP [see
c           below],  JP is actually the vibrational quantum number v",
c           FREQ the Bv value & UFREQ its uncertainty
c   LABLP = 'WI'  identifies this data group/band as a set of tunneling 
c           predissociation widths for electronic state LABLPP.  In this
c           case, parameters VP, VPP and EFP are dummy variables, while
c           the predissociating level is identified as: v"=JP, J"=JPP,
c           and parity p"=EFPP.
c   NOTE: !!!!!!!!!!! This option is ignored by DSParFit !!!!!!!!!!!!!!!
c** STOP reading when run out of bands OR when read-in VPP is negative   
c-----------------------------------------------------------------------
      READ(4,*,END=20) VP(IBAND), VPP(IBAND), LABLP, LABLPP, MN1,MN2
c-----------------------------------------------------------------------
      IF(VP(IBAND).LT.0) GO TO 20
      IEP(IBAND)= -99
      IEPP(IBAND)= -99
      DO  I= -3,NSTATES
          IF(LABLP.EQ.SLABL(I)) IEP(IBAND)= I
          IF(LABLPP.EQ.SLABL(I)) IEPP(IBAND)= I
          ENDDO
c** Check that this isotopomer is one of those chosen to be fitted ...
      ISOT= 0
      DO  I= 1,NISTP
          IF((MN1.EQ.MN(1,I)).AND.(MN2.EQ.MN(2,I))) ISOT= I
          ENDDO
      ISTP(IBAND)= ISOT
      TOTUFREQ= 0.D0
      MAXUFREQ(IBAND)= 0
      JMAX(IBAND)= 0
      JMIN(IBAND)= 9999
      COUNT= COUNT+1
      IF(COUNT.GT.NDATAMX) THEN
          WRITE(6,640) COUNT,NDATAMX
          STOP
          ENDIF
      NTRANS= 0
      IFIRST(IBAND)= COUNT
      ESP= IEP(IBAND)
      ESPP= IEPP(IBAND)
      VMAXespp= VMAX(ESPP)
      VMINespp= VMIN(ESPP)
      JTRUNCespp= JTRUNC(ESPP)
      IF(ISOT.GT.1) THEN
          VMAXespp= INT((VMAX(ESPP)+0.5d0)/RSQMU(ISOT)-0.5d0)
          VMINespp= INT((VMIN(ESPP)+0.5d0)/RSQMU(ISOT)-0.5d0)
          JTRUNCespp= INT(JTRUNC(ESPP)/RSQMU(ISOT))
          ENDIF
      VMAXesp= VMAX(ESPP)
      IF(ESP.GT.0) THEN
          VMAXesp= VMAX(ESP)
          VMINesp= VMIN(ESP)
          JTRUNCesp= JTRUNC(ESP)
          IF(ISOT.GT.1) THEN
              VMAXesp= INT((VMAX(ESP)+ 0.5d0)/RSQMU(ISOT) - 0.5d0)
              VMINesp= INT((VMIN(ESP)+ 0.5d0)/RSQMU(ISOT) - 0.5d0)
              JTRUNCesp= INT(JTRUNC(ESP)/RSQMU(ISOT))
              ENDIF
          ENDIF
c** For each of the lines in a given band/series, read upper level
c  rotational quantum number (JP) and e/f parity [EFP= +1 for e, -1 for
c  f, and  0 if e/f splitting unresolved and to  be ignored], and lower
c  level rotational quantum number (JPP) and parity [EFPP, as above],
c  the transition frequency  FREQ, and its uncertainty UFREQ.
c** For PAS or Tunneling Width data,  JP(COUNT)=v", JPP(COUNT)=J", 
c  EFPP(COUNT)=p", FREQ is the observable (a positive No.), while 
c  EFP(COUNT), VP(IBAND) & VPP(IBAND) are dummy variables.
c** For Bv values, JP(COUNT)=v" while JPP(COUNT), EFP(COUNT) and
c   EFPP(COUNT) as well as VP(IBAND) & VPP(IBAND) are dummy variables.
c-----------------------------------------------------------------------
   15 READ(4,*) JP(COUNT), EFP(COUNT), JPP(COUNT), EFPP(COUNT), 
     1                                       FREQ(COUNT), UFREQ(COUNT)
c-----------------------------------------------------------------------
c=======================================================================
c   Sample IR band data of HF for the '.4' file:
c   --------------------------------------------                          
c   1 0  'X0' 'X0'  1 19             % VP VPP LABLP LABLPP MN1 MN2
c   8 1   9 1  266.0131002  0.005    % JP EFP JPP EFPP FREQ UFREQ
c   9 1  10 1  265.8885896  0.003
c  10 1  11 1  265.7716591  0.002
c   .    .      .            .
c   .    .      .            .
c   [end of a band indicated by -ve JP and/or JPP value(s)]
c  -1 1  -1 1  -1.1         -1.1
c=======================================================================
      IF(EFP(COUNT).GT.1) EFP(COUNT)= 1
      IF(EFP(COUNT).LT.-1) EFP(COUNT)= -1
      IF(EFPP(COUNT).GT.1) EFPP(COUNT)= 1
      IF(EFPP(COUNT).LT.-1) EFPP(COUNT)= -1
c** At end of a band, exit from implicit loop
      IF((JPP(COUNT).LT.0).OR.(JP(COUNT).LT.0)) GOTO 18
c** If this band is not for one of the isotopomers chosen to be fitted,
c  omit its data from the fit
      IF(ISOT.EQ.0) GO TO 15
c** If this band involves electronic states other than those chosen to 
c   be treated, omit its data from the fit
      IF((ESP.EQ.-99).OR.(ESPP.EQ.-99)) GO TO 15
c** If a datum uncertainty of zero is accidentally read in, STOP
      IF(DABS(UFREQ(COUNT)).LE.0.d0) THEN
          WRITE(6,600) COUNT,FREQ(COUNT),IBAND
          STOP
          ENDIF
c** Omit data  with uncertainties outside specified limit UCUTOFF
      IF(UFREQ(COUNT).GT.UCUTOFF) GOTO 15
c** Require that datum lies within specified J & v ranges
      IF(ESP.GE.-2) THEN
          IF(((JTRUNCespp.GE.0).AND.(JPP(COUNT).GT.JTRUNCespp)).OR.
     1       ((JTRUNCespp.LT.0).AND.(JPP(COUNT).LT.-JTRUNCespp)))
     2                                                         GOTO 15
          IF((EFPP(COUNT)*EFSEL(ESPP)).LT.0) GOTO 15
          ENDIF
      IF(ESP.GT.0) THEN
          IF(VPP(IBAND).GT.VMAXespp) GOTO 15
          IF(VPP(IBAND).LT.VMINespp) GOTO 15
          IF(VP(IBAND).GT.VMAXesp) GOTO 15
          IF(VP(IBAND).LT.VMINesp) GOTO 15
          IF((JTRUNCesp.GE.0).AND.(JP(COUNT).GT.JTRUNCesp)) GOTO 15
          IF((JTRUNCesp.LT.0).AND.(JP(COUNT).LT.-JTRUNCesp)) GOTO 15
          IF((EFP(COUNT)*EFSEL(ESP)).LT.0) GOTO 15
        ELSE
          IF(JP(COUNT).GT.VMAXespp) GOTO 15
          IF(JP(COUNT).LT.VMINespp) GOTO 15
        ENDIF
c** If NOWIDTHS > 0  omit any tunneling width data from the fit.
      IF((ESP.EQ.-2).AND.(NOWIDTHS.GT.0)) GOTO 15
c
c** End of tests for datum inclusion.  Now count/sort data
c=======================================================================
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c%%% Convert  MHz  to  cm-1
c     freq(count)=freq(count)/2.99792458d+4
c     ufreq(count)=ufreq(count)/2.99792458d+4
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      TVUP(COUNT)= 0
      TVLW(COUNT)= 0
      IF(ESP.GE.-1) UMIN= MIN(UMIN,UFREQ(COUNT))
c** Determine actual v & J range of data & count data for each v
c  JMIN & JMAX needed for printout summary & data-count for testing
c  no. parameters allowed in Band Constant fit.
c??? This segment imperfect & needs re-examination ?????????????
      IF(ESP.GT.0) THEN
          IF(JPP(COUNT).LT.JMIN(IBAND)) JMIN(IBAND)= JPP(COUNT) 
          IF(JPP(COUNT).GT.JMAX(IBAND)) JMAX(IBAND)= JPP(COUNT)
          IF(JP(COUNT).LT.JMIN(IBAND)) JMIN(IBAND)= JP(COUNT) 
          IF(JP(COUNT).GT.JMAX(IBAND)) JMAX(IBAND)= JP(COUNT)
          VMX(ESP)= MAX(VMX(ESP),VP(IBAND))
          VMX(ESPP)= MAX(VMX(ESPP),VPP(IBAND))
c
c** Accumulate count of data associated with each vibrational level ...
          NDAT(VPP(IBAND),ISTP(IBAND),ESPP)=
     1                            NDAT(VPP(IBAND),ISTP(IBAND),ESPP)+ 1
          NDAT(VP(IBAND),ISTP(IBAND),ESP)=
     1                              NDAT(VP(IBAND),ISTP(IBAND),ESP)+ 1
        ELSEIF((ESP.LE.0).OR.(ESP.GE.-2)) THEN
          IF(JP(COUNT).LT.JMIN(IBAND)) JMIN(IBAND)= JP(COUNT) 
          IF(JP(COUNT).GT.JMAX(IBAND)) JMAX(IBAND)= JP(COUNT)
          VMX(ESPP)= MAX(VMX(ESPP),JP(COUNT))
          NDAT(JP(COUNT),ISTP(IBAND),ESPP)=
     1                      NDAT(JP(COUNT),ISTP(IBAND),ESPP)+ 1
        ELSEIF(ESP.EQ.-3) THEN
c... and for Bv data ...
          IF(JPP(COUNT).LT.JMIN(IBAND)) JMIN(IBAND)= JPP(COUNT) 
          IF(JPP(COUNT).GT.JMAX(IBAND)) JMAX(IBAND)= JPP(COUNT)
          NDAT(JPP(COUNT),ISTP(IBAND),ESPP)=
     1                     NDAT(JPP(COUNT),ISTP(IBAND),ESPP)+ 1
        ENDIF
      DFREQ(COUNT)= 0.d0
      IB(COUNT)= IBAND 
      TOTUFREQ= TOTUFREQ+UFREQ(COUNT) 
      IF(UFREQ(COUNT).GT.MAXUFREQ(IBAND)) MAXUFREQ(IBAND)= UFREQ(COUNT)
      COUNT= COUNT+1 
      IF(COUNT.GT.NDATAMX) THEN
          WRITE(6,640) COUNT,NDATAMX
          STOP
          ENDIF
      GOTO 15 
c** End of loop reading data for a given band/series 
c
c** Tidy up at end of reading for a given band
   18 COUNT= COUNT-1
      ILAST(IBAND)= COUNT 
      NTRANS= ILAST(IBAND)-IFIRST(IBAND)+1
      IF(NTRANS.GT.0) THEN
c** Treat PAS data as Fluorescence series unless  PASok > 0
          IF((IEP(IBAND).EQ.-1).AND.(PASok(IEPP(IBAND)).LE.0)) 
     1                                                    IEP(IBAND)=0
          IF((NTRANS.EQ.1).AND.(LABLP.EQ.'FS')) THEN
c** Ignore any fluorescence series consisting of only one datum
              COUNT= COUNT-1
              IBAND= IBAND-1
              FSOMIT= FSOMIT+1
              GOTO 10
              ENDIF
          AVEUFREQ(IBAND)= TOTUFREQ/NTRANS
          NBANDS(ISTP(IBAND))= NBANDS(ISTP(IBAND))+1
        ELSE
          IBAND= IBAND-1
          GOTO 10
        ENDIF
c=======================================================================
c** Accumulate counters for bands/series of different types
      IF(ESP.EQ.0) THEN
c** For Fluorescence Series ... first enumerate the No. of bands & lines
          NFSTOT= NFSTOT+1
          FSBAND(NFSTOT)= IBAND
c** Define counter to label which f.s. is associated with band IBAND 
          NFS(IBAND)= NFSTOT
          NBANDFS(ISOT,ESPP)= NBANDFS(ISOT,ESPP)+1
          NBND= NBANDFS(ISOT,ESPP)
          NTRANSFS(ISOT,ESPP)= NTRANSFS(ISOT,ESPP)+NTRANS
c ... and then set up labels/ranges/properties for each band
          YPR(ISOT,ESPP,1,1,NBND)= VP(IBAND)
          YPR(ISOT,ESPP,1,2,NBND)= VPP(IBAND)
          YPR(ISOT,ESPP,1,3,NBND)= NTRANS
          YPR(ISOT,ESPP,1,4,NBND)= IBAND
          YPR(ISOT,ESPP,1,5,NBND)= JMIN(IBAND)
          YPR(ISOT,ESPP,1,6,NBND)= JMAX(IBAND)
          ENDIF
c
      IF((ESP.GT.0).AND.(ESP.NE.ESPP)) THEN
c** For vibrational band of a normal 2-state electronic transition
c ... count bands and transitions in visible (electronic) spectrum
          NBANDEL(ISOT,ESP,ESPP)= NBANDEL(ISOT,ESP,ESPP)+ 1
          NBANDVIS(ISOT,ESPP)= NBANDVIS(ISOT,ESPP)+ 1
          NBND= NBANDVIS(ISOT,ESPP)
          NTRANSVIS(ISOT,ESP,ESPP)= NTRANSVIS(ISOT,ESP,ESPP)+NTRANS
c ... and then set up labels/ranges/properties for each of them
          YPR(ISOT,ESPP,2,1,NBND)= VPP(IBAND)
          YPR(ISOT,ESPP,2,2,NBND)= VP(IBAND)
          YPR(ISOT,ESPP,2,3,NBND)= NTRANS
          YPR(ISOT,ESPP,2,4,NBND)= IBAND
          YPR(ISOT,ESPP,2,5,NBND)= JMIN(IBAND)
          YPR(ISOT,ESPP,2,6,NBND)= JMAX(IBAND)
          ENDIF 
c
      IF((ESP.EQ.ESPP).AND.(VP(IBAND).NE.VPP(IBAND))) THEN
c** For an Infrared band of electronic state  s=ESPP=ESP
c** First cumulatively count the number of IR bands & transitions
          NBANDIR(ISOT,ESPP)= NBANDIR(ISOT,ESPP)+1
          NBND= NBANDIR(ISOT,ESPP)
          NTRANSIR(ISOT,ESPP)= NTRANSIR(ISOT,ESPP)+NTRANS 
c ... and then set up labels/ranges/properties for each of them
          YPR(ISOT,ESPP,3,1,NBND)= VPP(IBAND)
          YPR(ISOT,ESPP,3,2,NBND)= VP(IBAND)
          YPR(ISOT,ESPP,3,3,NBND)= NTRANS
          YPR(ISOT,ESPP,3,4,NBND)= IBAND
          YPR(ISOT,ESPP,3,5,NBND)= JMIN(IBAND)
          YPR(ISOT,ESPP,3,6,NBND)= JMAX(IBAND)
          ENDIF
c
      IF((ESP.EQ.ESPP).AND.(VP(IBAND).EQ.VPP(IBAND))) THEN
c** For Microwave transitions in electronic state  s=ESPP=ESP
c** First cumulatively count the number of MW bands & transitions
          NBANDMW(ISOT,ESPP)= NBANDMW(ISOT,ESPP)+1
          NBND= NBANDMW(ISOT,ESPP)
          NTRANSMW(ISOT,ESPP)= NTRANSMW(ISOT,ESPP)+NTRANS
c ... and then set up labels/ranges/properties for each of them
          YPR(ISOT,ESPP,4,1,NBND)= VPP(IBAND)
          YPR(ISOT,ESPP,4,2,NBND)= VP(IBAND)
          YPR(ISOT,ESPP,4,3,NBND)= NTRANS
          YPR(ISOT,ESPP,4,4,NBND)= IBAND
          YPR(ISOT,ESPP,4,5,NBND)= JMIN(IBAND)
          YPR(ISOT,ESPP,4,6,NBND)= JMAX(IBAND)
          ENDIF
c
c** NOTE ... in YPR array a last index counts bands of this type for 
c  this isotopomer of this electronic state ... and put all Bv's, 
c  Tunneling Widths or PAS binding energies in one group.
      IF(ESP.EQ.-3) THEN
c** Data are not transition energies, but rather the values of Bv in
c  electronic state s=IEPP  [As in the published IBr(A-X) analysis].
ccc       IF((NBVPP(ISOT,ESPP).GT.0).AND.(NTRANS.GT.0)) THEN
              WRITE(6,612) ESPP,ISOT
ccc           STOP
ccc           ENDIF
          NBVPP(ISOT,ESPP)= NTRANS
          YPR(ISOT,ESPP,5,3,1)= NTRANS
          YPR(ISOT,ESPP,5,4,1)= IBAND
          YPR(ISOT,ESPP,5,5,1)= JMIN(IBAND)
          YPR(ISOT,ESPP,5,6,1)= JMAX(IBAND)
          ENDIF
c
      IF(ESP.EQ.-2) THEN
c** Data are tunneling predissociation linewidths (in cm-1) for levels
c  of electronic state IEPP=ESPP
ccc       IF((NWIDTH(ISOT,ESPP).GT.0).AND.(NTRANS.GT.0)) THEN
              WRITE(6,626) ESPP,ISOT
ccc           STOP
ccc           ENDIF
          NWIDTH(ISOT,ESPP)= NTRANS
          YPR(ISOT,ESPP,6,3,1)= NTRANS
          YPR(ISOT,ESPP,6,4,1)= IBAND
          YPR(ISOT,ESPP,6,5,1)= JMIN(IBAND)
          YPR(ISOT,ESPP,6,6,1)= JMAX(IBAND)
          ENDIF
c
      IF(ESP.EQ.-1) THEN
c** Data are PhotoAssociation Binding Energies (in cm-1) for levels
c  of electronic state IEPP=ESPP
          WRITE(6,636) LABLPP,ISOT
          NEBPAS(ISOT,ESPP)= NTRANS
          YPR(ISOT,ESPP,7,3,1)= NTRANS
          YPR(ISOT,ESPP,7,4,1)= IBAND
          YPR(ISOT,ESPP,7,5,1)= JMIN(IBAND)
          YPR(ISOT,ESPP,7,6,1)= JMAX(IBAND)
          ENDIF
c** Now return to read the next band
      GOTO 10
c========================================================================
c** Now, write a summary of the input data to the output file
   20 COUNTOT= COUNT
      NBANDTOT= 0
      DO  I= 1,NISTP
          NBANDTOT= NBANDTOT+ NBANDS(I)
          ENDDO
      ISOT= 1
      UCUTOFF= UMIN
      IF(FSOMIT.GT.0) WRITE(6,650) FSOMIT
      IF(PRINP.LE.0) RETURN
c** Print a summary of the data, one isotopomer at a time.
   26 WRITE(6,602) NBANDS(ISOT), (NAME(I),MN(I,ISOT),I=1,2)
c
      DO 50 ISTATE= 1,NSTATES
c ... For internal use, may wish to update VMAX(ISTATE) to the actual 
c  highest v in the data set for this state. ** Reactivate as needed.
c      VMAX(ISTATE)= VMX(ISTATE)
c ... and separately list data for each (lower) electronic state in turn
      IF(NTRANSMW(ISOT,ISTATE).GT.0) THEN
c** Book-keeping for Micowave data
          WRITE(6,604) NTRANSMW(ISOT,ISTATE),SLABL(ISTATE),(NAME(I),
     1                          MN(I,ISOT),I=1,2),NBANDMW(ISOT,ISTATE)
          DO  I= 1,NBANDMW(ISOT,ISTATE)
              IBB= YPR(ISOT,ISTATE,4,4,I)
              WRITE(6,606) YPR(ISOT,ISTATE,4,2,I),
     1                     YPR(ISOT,ISTATE,4,1,I),
     2                  YPR(ISOT,ISTATE,4,3,I),YPR(ISOT,ISTATE,4,5,I),
     3                  YPR(ISOT,ISTATE,4,6,I), 
     3                  AVEUFREQ(YPR(ISOT,ISTATE,4,4,I)),
     4                  MAXUFREQ(YPR(ISOT,ISTATE,4,4,I))
              ENDDO
	    ENDIF
c
      IF(NTRANSIR(ISOT,ISTATE).GT.0)THEN
c** Book-keeping for Infrared data
          WRITE(6,608) NTRANSIR(ISOT,ISTATE),SLABL(ISTATE),(NAME(I),
     1                          MN(I,ISOT),I=1,2),NBANDIR(ISOT,ISTATE)
          DO  I= 1,NBANDIR(ISOT,ISTATE)
              IBB= YPR(ISOT,ISTATE,3,4,I)
              WRITE(6,606) YPR(ISOT,ISTATE,3,2,I),
     1                     YPR(ISOT,ISTATE,3,1,I),
     2                  YPR(ISOT,ISTATE,3,3,I),YPR(ISOT,ISTATE,3,5,I),
     3                  YPR(ISOT,ISTATE,3,6,I), 
     4                  AVEUFREQ(YPR(ISOT,ISTATE,3,4,I)),
     5                  MAXUFREQ(YPR(ISOT,ISTATE,3,4,I))
              ENDDO
          ENDIF
c
c** Book-keeping for electronic vibrational band data
      DO  ISTATEE= 1,NSTATES
          IF((ISTATEE.NE.ISTATE).AND.
     1                 (NTRANSVIS(ISOT,ISTATEE,ISTATE).GT.0)) THEN
c ... for ISTATEE{upper}-ISTATE{lower} electronic vibrational bands
              WRITE(6,610) NTRANSVIS(ISOT,ISTATEE,ISTATE),
     1         (NAME(I),MN(I,ISOT),I=1,2),SLABL(ISTATEE),SLABL(ISTATE),
     2                                    NBANDEL(ISOT,ISTATEE,ISTATE)
              DO  I= 1,NBANDVIS(ISOT,ISTATE)
                  IBB= YPR(ISOT,ISTATE,2,4,I)
                  IF(IEP(IBB).EQ.ISTATEE) THEN
                      WRITE(6,606) YPR(ISOT,ISTATE,2,2,I),
     1                            YPR(ISOT,ISTATE,2,1,I),
     2                  YPR(ISOT,ISTATE,2,3,I),YPR(ISOT,ISTATE,2,5,I),
     3                  YPR(ISOT,ISTATE,2,6,I), 
     4                  AVEUFREQ(YPR(ISOT,ISTATE,2,4,I)),
     5                  MAXUFREQ(YPR(ISOT,ISTATE,2,4,I))
                      ENDIF
                  ENDDO
              ENDIF
          ENDDO
      IF(NTRANSFS(ISOT,ISTATE).GT.0)THEN
c** Book-keeping for Fluorescence data
          WRITE(6,614) NTRANSFS(ISOT,ISTATE),SLABL(ISTATE),(NAME(I),
     1                          MN(I,ISOT),I=1,2),NBANDFS(ISOT,ISTATE)
          DO  I= 1,NBANDFS(ISOT,ISTATE)
              IBB= YPR(ISOT,ISTATE,1,4,I)
              WRITE(6,616) YPR(ISOT,ISTATE,1,1,I),
     1                   YPR(ISOT,ISTATE,1,2,I),NEF(EFP(IFIRST(IBB))),
     2                  YPR(ISOT,ISTATE,1,3,I),YPR(ISOT,ISTATE,1,5,I),
     3                  YPR(ISOT,ISTATE,1,6,I), 
     4                  AVEUFREQ(YPR(ISOT,ISTATE,1,4,I)),
     5                  MAXUFREQ(YPR(ISOT,ISTATE,1,4,I))
              ENDDO
          ENDIF
      IF(NBVPP(ISOT,ISTATE).GT.0)THEN
c** Book-keeping for  Bv  data
          WRITE(6,618) NBVPP(ISOT,ISTATE),SLABL(ISTATE),(NAME(I),
     1                                               MN(I,ISOT),I=1,2)
          IBB= YPR(ISOT,ISTATE,5,4,1)
          WRITE(6,620) YPR(ISOT,ISTATE,5,3,1),YPR(ISOT,ISTATE,5,5,1),
     1       YPR(ISOT,ISTATE,5,6,1),AVEUFREQ(YPR(ISOT,ISTATE,5,4,1)),
     2                               MAXUFREQ(YPR(ISOT,ISTATE,5,4,1))
          ENDIF
      IF(NWIDTH(ISOT,ISTATE).GT.0) THEN
c** Book-keeping for  Tunneling Width  data
          WRITE(6,628) NWIDTH(ISOT,ISTATE),SLABL(ISTATE),(NAME(I),
     1                                               MN(I,ISOT),I=1,2)
          IBB= YPR(ISOT,ISTATE,6,4,1)
          WRITE(6,630) YPR(ISOT,ISTATE,6,3,1),
     1                  YPR(ISOT,ISTATE,6,5,1),YPR(ISOT,ISTATE,6,6,1),
     2                               AVEUFREQ(YPR(ISOT,ISTATE,6,4,1)),
     3                                MAXUFREQ(YPR(ISOT,ISTATE,6,4,1))
          ENDIF
      IF(NEBPAS(ISOT,ISTATE).GT.0) THEN
c** Book-keeping for  PAS Binding Energy  data
          WRITE(6,632) NEBPAS(ISOT,ISTATE),SLABL(ISTATE),(NAME(I),
     1                                               MN(I,ISOT),I=1,2)
          IBB= YPR(ISOT,ISTATE,6,4,1)
          WRITE(6,630) YPR(ISOT,ISTATE,7,3,1),
     1                  YPR(ISOT,ISTATE,7,5,1),YPR(ISOT,ISTATE,7,6,1),
     2                               AVEUFREQ(YPR(ISOT,ISTATE,7,4,1)),
     3                                MAXUFREQ(YPR(ISOT,ISTATE,7,4,1))
          ENDIF
   50 CONTINUE
      IF(ISOT.LT.NISTP) THEN
c** If NISTP > 1, return to print data summaries for other isotopomers
          ISOT= ISOT+1
          GO TO 26
          ENDIF 
      WRITE(6,622)
      RETURN
  600 FORMAT(/' *** INPUT ERROR ***  Datum   FREQ(',i5,')=',f12.4,
     1 '  in   IBAND=',i4,'   has zero uncertainty!!!')
  601 FORMAT(23x,'or with',A3,'-parity.')
  603 FORMAT(/' Neglect data with:  Uncertainties > UCUTOFF=',G12.3,
     1  ' (cm-1)')
  605 FORMAT(7x,'and State ',A2,' data with  J < JTRUNC=',I4,
     1  '  or  v  outside range',i3,'  to',i4)
  607 FORMAT(7x,'and State ',A2,' data with  J > JTRUNC=',I4,
     2  '  or  v  outside range',i3,'  to',i4)
  602 FORMAT(/1x,20('===')/'  *** Input data for',i5,' bands/series of '
     1  ,A2,'(',I3,')-',A2,'(',I3,') ***'/1x,20('==='))
  604 FORMAT(1x,28('--')/I5,' State ',A2,1x,A2,'(',I3,')-',A2,'(',I3,
     1 ') MW transitions in',i4,' sets'/1x,28('--')/"   v'  ",
     1 'v"  #data   Jmin   Jmax  Avge.Unc.  Max.Unc.'/1x,25('--'))
  606 FORMAT(I4,I4,3I7,1x,1P2D10.1)
  608 FORMAT(1x,32('--')/I6,' State ',A2,1x,A2,'(',I3,')-',A2,'(',I3,
     1 ') InfraRed transitions in',I4,' bands'/1x,32('--')/
     2 "   v'  ",'v"  #data   Jmin   Jmax  Avge.Unc.  Max.Unc.'/
     3 1x,25('--'))
  609 FORMAT(/' *** ERROR *** Dimension allocated for number of bands ex
     1ceeded:'/' (IBAND=',i4,') > (NBANDMX=',i4,')   so truncate input a
     2nd TRY to continue ...')
  610 FORMAT(/1x,35('==')/I6,1x,A2,'(',I3,')-',A2,'(',i3,')  {State ',
     1  A2,'}--{State ',A2,'} Transitions in',i4,' Bands'/1x,35('--')/
     2 "   v'",'  v"  #data   Jmin   Jmax  Avge.Unc.  Max.Unc.'/
     3 1x,25('--'))
  612 FORMAT(/" NOTE that all read-in Bv's for   ISTATE=",i2,'   ISOT=',
     1  i2/32x,' must be input as a single "band" or data group')
cc612 FORMAT(/" *** STOP INPUT *** and put all read-in Bv's for   ISTATE
cc   1=",i2,'   ISOT=',i2/ 10x,'into one "band" or data group.')
  614 FORMAT(1x,38('==')/I5,' Fluorescence transitions into State ',
     1 A2,2x,A2,'(',I3,')-',A2,'(',I3,')  in',i5,' series'/
     2 1x,38('==')/"   v'  j' p' ",'#data  v"min  v"max  Avge.Unc.  Max.
     3Unc.'/1x,51('-'))
  616 FORMAT(2I4,A3,I6,2I7,1x,1P2D10.1)
  618 FORMAT(1x,65('=')/1x,I3,' State ',A2,1x,A2,'(',I3,')-',A2,'(',I3,
     1 ') Bv values treated as independent data'/1x,24('--')/
     2 '  #values   v(min)  v(max)  Avge.Unc.   Max.Unc.'/
     3 1x,24('--'))
  620 FORMAT(I7,I9,I8,3x,1P2D11.1)
  622 FORMAT(1x,25('===')/1x,25('==='))
  626 FORMAT(/" NOTE that all read-in Tunneling Widths for   ISTATE=",
     1 i2,'   ISOT=',i2/10x,' must be in a single "band" or data group')
cc626 FORMAT(/" *** STOP INPUT *** and put all read-in Tunneling Widths'
cc   1  '  for   ISTATE=",i2,'   ISOT=',i2/ 
cc   2  10x,'into one "band" or data group.')
  628 FORMAT(1x,61('=')/1x,I3,' State ',A2,1x,A2,'(',I3,')-',A2,'(',I3,
     1 ') Tunneling Widths included as data'/
     2 1x,61('-')/'  #values   v(min)  v(max)   Avge.Unc.   Max.Unc.'/
     3 1x,24('--'))
  630 FORMAT(I7,I9,I8,2x,1P2D11.1)
  632 FORMAT(1x,70('=')/I4,' State ',A2,1x,A2,'(',I3,')-',A2,'(',I3,
     1 ') PAS Binding Energies included in data set'/
     2 1x,70('-')/'  #values   v(min)  v(max)   Avge.Unc.   Max.Unc.'/
     3 1x,24('--'))
  636 FORMAT(/' NOTE that all read-in PAS Binding Energies for   ISTATE=
     1 ',a2,'  ISOT=',i2/10x,' must be in a single "band" or data group'
     2 )
  640 FORMAT(/' *** Input Data Count reaches',i6,' which EXCEEDS ARRAY L
     1IMIT of',i6)
  650 FORMAT(/' Data input IGNORES',i4,' fluorescence series consisting'
     1 ,' of only  onee  line!')
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE TVSORT(ISTATE,NPTOT,VMAX,NTVALL,TVNAME)
c***********************************************************************
c** Subroutine to sort through global data file, and for each isotopomer
c  in state ISTATE:  (1) find the number of transitions coupled to each
c  level (v,J,p),  (2) for levels in order (v,J,p), add a free parameter
c  for each level involved in one or more transitions, and  (3) label each
c  transition involving one of these levels by the index/counter of the
c  parameter associated with that term value.
c             ********* Version of 27 August 2004 *********
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** On Entry:
c------------
c  ISTATE is the electronic state being considered.
c  NPTOT  enters as the cumulative count of parameters prior to entry
c  TVUP(i) and TVLW(i) in COMMON equal zero for all data
c
c** On Return:
c-------------
c  NPTOT  is updated to include the number of term values for this state
c  TVUP(i) & TVLW(i): if the upper and/or lower level of transition-i is
c      to be represented by a term value, TVUP and TVLW (respectively)
c      is the associated parameter index; otherwise they = 0.
c  NTVALL  is the number of term value parameters for this state
c  TVNAME(j)  is the alphameric name identifying term value parameter j
c
c** Internally
c-------------
c  NLV(v,J.p) * initially, counts transitions for level {v,J,p} of a 
c                          given isotopologue
c           * later reset it as the parameter index for that term value
c-----------------------------------------------------------------------
      INCLUDE 'arrsizes.h'
      INCLUDE 'BLKISOT.h'
      INCLUDE 'BLKDATA.h'
c
      INTEGER I,J,P,IBAND,ISOT,ISTATE,NPTOT,LOWEST,VMAX(NSTATEMX),
     1 NLV(0:NVIBMX,0:NVIBMX,-1:1),NTVS(NSTATEMX,NISTPMX),
     2 NTVALL(0:NSTATEMX)
      CHARACTER*24 TVNAME(NPARMX)
c=======================================================================
      WRITE(6,600) SLABL(ISTATE) 
      LOWEST= 1
      IF(ISTATE.GT.1) LOWEST= 0
      NTVALL(ISTATE)= 0
      DO  ISOT= 1, NISTP
c** First ... zero transition counter array for this isotopomer
          DO  I= 0, VMAX(ISTATE)
cc            DO  J= 0, NROTMX 
              DO  J= 0, NVIBMX 
                  DO  P= -1,1
                      NLV(I,J,P)= 0
                      ENDDO
                  ENDDO
              ENDDO
          DO  IBAND= 1, NBANDTOT
c** Then ... search for bands involving isotopomer ISOT in this state
              IF(((IEP(IBAND).EQ.ISTATE).OR.(IEPP(IBAND).EQ.ISTATE))
     1          .AND.(ISTP(IBAND).EQ.ISOT).AND.(IEP(IBAND).GE.0)) THEN
                  DO  I= IFIRST(IBAND), ILAST(IBAND)
c ... for each such band, loop over all transitions, and increment NLV 
c     for each {v,J,p} level encountered in a transision
                      IF(IEP(IBAND).EQ.ISTATE) THEN
                          IF(JP(I).GT.NVIBMX) THEN
c ... check for array dimension overruns
                              WRITE(6,602) ISTATE,ISOT,JP(I),NROTMX
                              STOP
                              ENDIF
                          NLV(VP(IBAND),JP(I),EFP(I))= 
     1                                  NLV(VP(IBAND),JP(I),EFP(I))+ 1
                          ENDIF
                      IF(IEPP(ISTATE).EQ.ISTATE) THEN
                          IF(JPP(I).GT.NVIBMX) THEN
                              WRITE(6,604) ISTATE,ISOT,JPP(I),NROTMX
                              STOP
                              ENDIF
                          NLV(VPP(IBAND),JPP(I),EFPP(I))
     1                             = NLV(VPP(IBAND),JPP(I),EFPP(I))+ 1
                          ENDIF
                      ENDDO
                  ENDIF
c** Finished scan over all data set for this isotopologue
              ENDDO
c
c** Now ... count a free parameter for each level in a transition
c** NTV  is the total number of term values for case (ISTATE,ISOT) 
c   NTVS is the no. of them involved in only a single transition
          NTV(ISTATE,ISOT)= 0
          NTVS(ISTATE,ISOT)= 0
          DO  I= 0, VMAX(ISTATE)
ccc           DO  J= 0, NROTMX
              DO  J= 0, NVIBMX
                  DO  P= -1,1
                      IF(NLV(I,J,P).GT.0) THEN
                          IF(LOWEST.EQ.1) THEN
c** If using term values for `lowest' state (defined as the first state
c  considered), its lowest observed level for isotopologue-1 defines the
c   absolute energy zero
                              WRITE(6,606) I,J,P,ISOT,SLABL(ISTATE)
                              LOWEST= 0
                              NLV(I,J,P)= 0
                              GOTO 20
                              ENDIF
                          NPTOT= NPTOT+ 1
                          NTV(ISTATE,ISOT)= NTV(ISTATE,ISOT)+ 1
                          IF(NLV(I,J,P).EQ.1) NTVS(ISTATE,ISOT)=
     1                                           NTVS(ISTATE,ISOT) +1 
                          REWIND(30)
                          WRITE(30,700) SLABL(ISTATE),I,J,P,ISOT,
     1                                                      NLV(I,J,P)
                          REWIND(30)
                          READ(30,*) TVNAME(NPTOT)
c ... reset NLV(v,J,p) as the parameter index for that term value
                          NLV(I,J,P)= NPTOT
                          ENDIF
   20                 CONTINUE
                      ENDDO
                  ENDDO
              ENDDO
c
c** Finally - label each transition with term-value parameter index for
c   (as appropriate) upper & lower level of each transition
          DO  IBAND= 1, NBANDTOT
              IF(((IEP(IBAND).EQ.ISTATE).OR.(IEPP(IBAND).EQ.ISTATE))
     1          .AND.(ISTP(IBAND).EQ.ISOT).AND.(IEP(IBAND).GE.0)) THEN
c ... for each band involving state ISTATE of this isotopologue, label 
c     each transition with the term value parameter index (which is zero
c     if the state is not represented by term values!).
                  DO  I= IFIRST(IBAND), ILAST(IBAND)
                      IF(IEP(IBAND).EQ.ISTATE) 
     1                          TVUP(I)= NLV(VP(IBAND),JP(I),EFP(I))
                      IF(IEPP(IBAND).EQ.ISTATE) 
     1                          TVLW(I)= NLV(VPP(IBAND),JPP(I),EFP(I))
                      ENDDO
                  ENDIF
              ENDDO
          WRITE(6,608) NAME(1),MN(1,ISOT),NAME(2),MN(2,ISOT),
     1                              NTV(ISTATE,ISOT),NTVS(ISTATE,ISOT)
          NTVALL(ISTATE)= NTVALL(ISTATE)+ NTV(ISTATE,ISOT)
          ENDDO
      RETURN
  600 FORMAT(/' For State ',A2,'  fit to individual term values for each
     1  {v,J,p,isot}'/1x,6('******'))
  602 FORMAT(/' *** ARRAY DIMENSION PROBLEM ***  JP(ISTATE)=',i2,
     1  ',ISOT=',I2,')=',i3,'  greater than  NVIBMX=',i4)
  604 FORMAT(/' *** ARRAY DIMENSION PROBLEM ***  JPP(ISTATE)=',i2,
     1  ',ISOT=',I2,')=',i3,'  greater than  NVIBMX=',i4)
  606 FORMAT(/'  Absolute zero of energy is fixed at level {v=',i3,
     1 ', J=',i3,', p=',i2,'}'/1x,12('**'),10x,'of isotopomer ',i2,
     2 ' of  State ',A2)
  608 FORMAT(' For ',A2,'(',i3,')-',A2,'(',I3,')  fit to',i4,
     1 ' T(v,J,p) term values,'/20x,'of which',i4,' are involved in only
     2 one transition')
  700 FORMAT("'",'T(',A2,':',i3,',',i3,',',SP,i2,';',SS,i2,')',I5,"'")
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c**********************************************************************N(
      SUBROUTINE READPOT(ISTATE,SLABL,OSEL)
c**********************************************************************
c** This subroutine will generate one of four possible types of
c   analytical molecular potentials for the direct Hamiltonian fitting
c   program. These four families are:
c     - The Expanded Morse Oscillator (EMO potential
c     - The Morse/Long-Range (MLR) or Morse/Lennard-Jones (MLJ) potential
c     - The Double-Exponential Long-Range (DELR) potential
c     - The Generalized Potential Energy Function (GPEF)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c                   Version of 7 November 2006
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** On entry:
c   ISTATE  is the electronic state being fitted to
c   SLABL   is the two-character label identifying that state
c-----------------------------------------------------------------------
      INCLUDE 'arrsizes.h'
      INCLUDE 'BLKISOT.h'
      INCLUDE 'BLKPOT.h'
      INCLUDE 'BLKBOB.h'
c-----------------------------------------------------------------------
c** Type statements for input or local variables
      INTEGER I, ISTATE, IISTP, IFXUAinf, IFXUBinf, IFXTAinf, IFXTBinf,
     1  m, MMN, npow, OSEL(NSTATEMX)
      CHARACTER*2 SLABL(-3:NSTATEMX)
      CHARACTER*3 MLX
      REAL*8 UAinf,UBinf,TAinf,TBinf,ZMASE
      DATA ZMASE /5.4857990945D-04/
c
c** Set some defaults for parameters not common to all models ...
      IFXDe(ISTATE)= 1
      IFXRe(ISTATE)= 1
      DO  m= 1, NCMMmax
          IFXCm(m,ISTATE)= 1
          ENDDO
c-----------------------------------------------------------------------
c** First choose potential model and select form of BOB representation
c  PSEL(s)  choses the type of analytical potential to be fitted to:
c         = -2 : represent each distinct observed level of this state as
c               an independent term value [an alternative to an 'FS' 
c               treatmentof transitions involving that state]
c         = 0 : Do forward calculation from pointwise input potential
c         = 1 : Use the Expanded Morse Oscillator (p)
c         = 2 : Use the Morse Lennard-Jones Potential.
c         = 3 : Use the Double-Exponential Long-Range Potential.
c         = 4 : Use the Generalized Potential Energy Function.
c  VLIM(s)  is the fixed absolute energy of the potential asymptote
c  BOBCN  is a flag to denote reference & scaling for BOB corrections 
c      = 0  using differences as per RJL [JMS 194,189(1999)]
c      = 1  use 'clamped nuclei' limit, m_e1/MASS scaling.
c  OSEL(s)  controls printout of radial function arrays to Ch. 10-16.
c       OSEL > 0: Export to file every OSEL'th point of final function
c=======================================================================
      READ(5,*) PSEL(ISTATE), VLIM(ISTATE), BOBCN(ISTATE), OSEL(ISTATE)
c=======================================================================
c** For term value fits ... no further reads needed!
      IF(PSEL(ISTATE).EQ. -2) THEN
          WRITE(6,604) SLABL(ISTATE)
          RETURN
          ENDIF
c-----------------------------------------------------------------------
c** Next define form of potential dunction expansion variable  y_p(r)
c  nPB(s)  is the power  p  for  phi(r)  exponent expansion variable
c  nPu(s)  is the power  p  for  u(r) function expansion variable
c  nPt(s)  is the power  p  for  g(r) function expansion variable 
c  RREF(s)  defines the reference distance in the potential exponent 
c    expansion variable:  * for  RREF.le.0 , define parameter  RREF = Re
c      * for  RREF.gt.0 , fix parameter  RREF   at its read-in value
c=======================================================================
      READ(5,*) nPB(ISTATE), nPu(ISTATE), nPt(ISTATE), RREF(ISTATE)
c=======================================================================
c** Now to read in the range and mesh for the numerical integration
c  RMIN/MAX(s)  define the range over which this potential is defined.
c  RH(s)  specifies radial mesh for numerical integration for this state
c=======================================================================
      READ(5,*) RMIN(ISTATE), RMAX(ISTATE), RH(ISTATE)
c=======================================================================
      NDATPT(ISTATE)= (RMAX(ISTATE)-RMIN(ISTATE))/RH(ISTATE)+1.0001d0
      NDATPT(ISTATE)= MIN(NPNTMX-1,NDATPT(ISTATE))
      RMAX(ISTATE)= RMIN(ISTATE) + RH(ISTATE)*DBLE(NDATPT(ISTATE)-1)
      IF(PSEL(ISTATE).EQ.0) WRITE(6,600) SLABL(ISTATE)
      IF(PSEL(ISTATE).EQ.2) THEN
c-----------------------------------------------------------------------
c** For MLJ or MLR, read number of inverse-power terms NCMM in long-range
c   form   V(R) -> VLIM - CmVAL(1)/R**MMLR(1) - CmVAL(2)/R**MMLR(2) - ...
c  If Asw > 0.0  fix  Cn = CmVAL(1)  using the switching fx.
c              fsw(R) = 1/[e^(Asw*(R - Rsw)) + 1]  
c     where  Rsw(s)  defines distance at which  fsw(Rsw) = 0.5
c  If Asw .le.0   phi(R)= yp*phiINF + [1-yp]* Sum{ phi_i * yp^i }
c=======================================================================
c=======================================================================
          READ(5,*) NCMM(ISTATE), Asw(ISTATE),Rsw(ISTATE)
          DO  m= 1,NCMM(ISTATE)
              READ(5,*) MMLR(m,ISTATE), CmVAL(m,ISTATE), IFXCm(m,ISTATE)
              ENDDO
c=======================================================================
c** For alkali dimer A-state: Aubert-Frecon ULR(r) [PRA 55, 3458 (1997)]
c  is special case with MMLR= 3, 0, 6, 6 {and 8, 8 possibly?}
c=======================================================================
          MMN= MMLR(NCMM(ISTATE),ISTATE)- MMLR(1,ISTATE)
          IF((NCMM(ISTATE).GT.1).AND.
     1                      ((nPB(ISTATE).LE.MMN).OR.(MMN.LE.0))) THEN
              WRITE(6,628) nPB(ISTATE),MMN
              NCMM(ISTATE)= 1
              MMN= 0 
              ENDIF
          MLX= 'MLR'
          IF(NCMM(ISTATE).GT.1) MLX= 'MLR'
          ENDIF
      IF(PSEL(ISTATE).EQ.3) THEN
c-----------------------------------------------------------------------
c** For DELR potential, read parameters to define form of `LR' part
c  IDF  define the form of damping function:  IDF=1 - Tang-Toennies fx.
c  IDF=2 use Douketis-Scoles Dm fx.  
c!!!!  IDF=3 also Scoles' global fx. !!!! this option commented out !!!!
c=======================================================================
          READ(5,*) NCMM(ISTATE), RHOd(ISTATE), IDF(ISTATE)
          DO  m= 1, NCMM(ISTATE)
              READ(5,*) MMLR(m,ISTATE), CmVAL(m,ISTATE)
              ENDDO
c=======================================================================
          IF(IDF(ISTATE).GE.3) IDF(ISTATE)= 2
          ENDIF
      IF(PSEL(ISTATE).EQ.4) THEN
c-----------------------------------------------------------------------
c** For GPEF potential, read parameters defining the expansion variable
c                                 p    p      p      p
c                  y(R;k,a,b) = (R - Re )/(a*R + b*Re )
c=======================================================================
          READ(5,*) AGPEF(ISTATE), BGPEF(ISTATE)
c=======================================================================
          RREF(ISTATE)= -1.d0
          ENDIF
      WRITE(6,626) SLABL(ISTATE),RMIN(ISTATE),RMAX(ISTATE),RH(ISTATE)
c
c** Now to read in the trial dissociation energy and equilibrium
c   radial distance for the state.
c  De(s)   is the dissociation energy for each state.
c  Re(s)   is the equilibrium radial distance for each state.
c  IFDe(s) indicates whether the dissociation energy will be:
c           = 1: held fixed at read-in values.
c          <= 0: determined from fits.
c  IFRe(s) indicates whether the equilibrium radial distance will be:
c           = 1: held fixed at read-in values.
c          <= 0: determined from fits.
c=======================================================================
      READ(5,*) DE(ISTATE), IFXDE(ISTATE)
      READ(5,*) RE(ISTATE), IFXRE(ISTATE)
c=======================================================================
c** Read in trial initial trial parameters for exponent \phi(r)
c
c  NSphi(s)  is the order of the phi_i exponent expansion for  R < Re
c  NLphi(s)  is the order of the phi_i exponent expansion for R > Re
c  PHI(b,s)   contains the parameters that define the Born-Oppenheimer
c              potential functions:
c        PSEL .ne. 4 : read-in values are  phi_i  parameters for the
c               exponent expansion of the 'clamped nuclei' potential.
c        PSEL = 4 :    read-in values are leading coefficients in 
c                Surkus' Generalized Potential Energy Function (GPEF).
c  IFXPHI(b,s) indicates whether each exponent potential expansion 
c             coefficient will be:    = 1: held fixed at read-in values.
c                                  .LE. 0: determined from fits.
c=======================================================================
      IF(PSEL(ISTATE).EQ.4) IFXDE(ISTATE)= 1
      READ(5,*) NSphi(ISTATE), NLphi(ISTATE)
      npow= MAX(NSphi(ISTATE),NLphi(ISTATE))
      DO  I= 0,npow
          READ(5,*) PHI(I,ISTATE), IFXPHI(I,ISTATE)
          ENDDO
c=======================================================================
c** Read parameters defining the BOB adiabatic radial functions
c  NUA/NUB(s)  specifies the order of the polynomial in  yp  defining 
c              the adiabatic BOB function for atom A/B
c         if < 0   do not read in any adiabatic BOB function parameters
c  UA/UB(a,s)   are the adiabatic BOB function expansion coefficients
c  IFU(A/B)(a,s) indicates whether each expansion coefficient is to be
c           > 0 :  held fixed at read-in value, or
c        .le. 0 :  varied in the fit
c  UAinf/UBinf  is the limiting asymptotic value of uA(r)/uB(r), as per
c               Theochem paper [internally stored as  UA(NUA+1), etc.]
c  IFXUAinf/IFXUBinf  specifies whether (>0) or not (.le.0)  UAinf/UBinf 
c               is to be held fixed at the read-in value
c=======================================================================
      READ(5,*) NUA(ISTATE)
      IF(NUA(ISTATE).GE.0) THEN
          DO  I= 0, NUA(ISTATE)
              READ(5,*) UA(I,ISTATE), IFXUA(I,ISTATE)
              ENDDO
          READ(5,*) UAinf, IFXUAinf
c=======================================================================
          IF(BOBCN(ISTATE).GT.0) THEN
              UAinf= 0.d0
              IFXUAinf= 1
              ENDIF
          NUA(ISTATE)= NUA(ISTATE)+ 1
          UA(NUA(ISTATE),ISTATE)= UAinf
          IFXUA(NUA(ISTATE),ISTATE)= IFXUAinf
          ENDIF
c=======================================================================
      READ(5,*) NUB(ISTATE)
      IF(NUB(ISTATE).GE.0) THEN
          DO  I= 0, NUB(ISTATE)
              READ(5,*) UB(I,ISTATE), IFXUB(I,ISTATE)
              ENDDO
	    READ(5,*) UBinf, IFXUBinf
c=======================================================================
          IF(BOBCN(ISTATE).GT.0) THEN
              UBinf= 0.d0
              IFXUBinf= 1
              ENDIF
          NUB(ISTATE)= NUB(ISTATE)+ 1
          UB(NUB(ISTATE),ISTATE)= UBinf
          IFXUB(NUB(ISTATE),ISTATE)= IFXUBinf
          ENDIF 
c***********************************************************************
c** Read parameters defining the BOB non-adiabatic centrifugal functions
c** If NISTP= 1 , read only one set of non-adiabatic parameters
c
c  NTA/NTB(s)  specifies the order of the polynomial in  yp  defining
c              the non-adiabatic centrifugal BOB functions for atom A/B
c         if < 0   do not read in any non-adiabatic BOB parameters
c  TA/TB(a,s)   are the non-adiabatic centrifugal BOB expansion coeffts
c  IFXTA/IFXTB(a,s) indicates whether each expansion coefficient is to be
c           > 0 :  held fixed at read-in value, or
c        .le. 0 :  varied in the fit
c  TAinf/TBinf  is the limiting asymptotic value of qA(r)/qB(r), as per 
c               Theochem paper [internally stored as  TA(NTA+1), etc.]
c  IFXTAinf/IFXTBinf  specifies whether (>0) or not (.le.0)  TAinf/TBinf 
c               is to be held fixed at the read-in value
c=======================================================================
      READ(5,*) NTA(ISTATE)
      IF(NTA(ISTATE).GE.0) THEN
          DO  I= 0, NTA(ISTATE)
              READ(5,*) TA(I,ISTATE), IFXTA(I,ISTATE)
              ENDDO
          READ(5,*) TAinf, IFXTAinf
c=======================================================================
          IF(BOBCN(ISTATE).GT.0) THEN
              TAinf= 0.d0
              IFXTAinf= 1
              ENDIF
          NTA(ISTATE)= NTA(ISTATE)+ 1
          TA(NTA(ISTATE),ISTATE)= TAinf
          IFXTA(NTA(ISTATE),ISTATE)= IFXTAinf 
          ENDIF 
c=======================================================================
      READ(5,*) NTB(ISTATE)
      IF(NTB(ISTATE).GE.0) THEN
          DO  I= 0, NTB(ISTATE)
              READ(5,*) TB(I,ISTATE), IFXTB(I,ISTATE)
              ENDDO
          READ(5,*) TBinf, IFXTBinf
c=======================================================================
          IF(BOBCN(ISTATE).GT.0) THEN
              TBinf= 0.d0
              IFXTBinf= 1
              ENDIF
          NTB(ISTATE)= NTB(ISTATE)+ 1
          TB(NTB(ISTATE),ISTATE)= TBinf
          IFXTB(NTB(ISTATE),ISTATE)= IFXTBinf
          ENDIF
c
      NwCFT(ISTATE)= -1
      IF(IOMEG(ISTATE).NE.0) THEN
c-----------------------------------------------------------------------
c** If electronic angular momentum not zero for this state, read Lambda 
c       doubling or doublet Sigma radialfunction parameters.
c*  NwCFT(s)  is order of the polynomial representing the radial fx.
c*  nPw(s)  defined nature of radial expansion variable:  
c         yp = [R^p - Re^p]/[R^p + Re^p]    with  p= nPw
c*  efREF(s) defines reference level for the Lambda doubling splitting
c            = -1  treats f level as the reference
c            =  0  treats the mid-point between e and f as reference
c            =  1  treats e level as the reference
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          READ(5,*) NwCFT(ISTATE), nPw(ISTATE), efREF(ISTATE)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          IF(NwCFT(ISTATE).GE.0) THEN
              DO  I= 0, NwCFT(ISTATE)
c** Read in  nLamba+1  Lambda-doubling correction coefficients
c=======================================================================
                  READ(5,*) wCFT(I,ISTATE), IFXwCFT(I,ISTATE)
                  ENDDO
              NwCFT(ISTATE)= NwCFT(ISTATE)+ 1
              READ(5,*) wCFT(NwCFT(ISTATE),ISTATE),
     1                                   IFXwCFT(NwCFT(ISTATE),ISTATE)
c=======================================================================
              IF(IOMEG(ISTATE).GT.0) THEN
                  IF(EFREF(ISTATE).EQ.-1) WRITE(6,640) SLABL(ISTATE)
                  IF(EFREF(ISTATE).EQ.0) WRITE(6,642) SLABL(ISTATE)
                  IF(EFREF(ISTATE).GE.1) WRITE(6,644) SLABL(ISTATE)
                  ENDIF
              ENDIF
          ENDIF
c
c** Calculate BOB mass scaling factors for the adiabatic (ZMUA, ZMUB) &
c   non-adiabaric centrifugal (ZMTA, ZMTB) BOB functions
      DO  IISTP= 1, NISTP
          IF(BOBCN(ISTATE).GE.1) THEN
c** For Watson/Coxon/Ogilvie-type clamped-nuclei reference species:
              ZMUA(IISTP,ISTATE)= ZMASE/ZMASS(1,IISTP)
              ZMUB(IISTP,ISTATE)= ZMASE/ZMASS(2,IISTP)
              ZMTA(IISTP,ISTATE)= ZMASE/ZMASS(1,IISTP)
              ZMTB(IISTP,ISTATE)= ZMASE/ZMASS(2,IISTP)
            ELSE
c** Using RJL's mass differences for adiabatic corrections (ZMUA, ZMUB):
              ZMUA(IISTP,ISTATE)= 1.0d0 - ZMASS(1,1)/ZMASS(1,IISTP)
              ZMUB(IISTP,ISTATE)= 1.0d0 - ZMASS(2,1)/ZMASS(2,IISTP)
c   and mass ratios for the rotational corrections (ZMTA, ZMTB):
              ZMTA(IISTP,ISTATE)= ZMASS(1,1)/ZMASS(1,IISTP)
              ZMTB(IISTP,ISTATE)= ZMASS(2,1)/ZMASS(2,IISTP)
            END IF
c
c** For homonuclear diatomics, set the first mass scaling term for each
c   set of correction terms to be the sum of the two original mass
c   scaling factors, and set the second mass term to zero.
c
          IF(AN(1).EQ.AN(2)) THEN
              ZMUA(IISTP,ISTATE)= ZMUA(IISTP,ISTATE)+ ZMUB(IISTP,ISTATE)
              ZMUB(IISTP,ISTATE)= 0.0d0
              ZMTA(IISTP,ISTATE)= ZMTA(IISTP,ISTATE)+ ZMTB(IISTP,ISTATE)
              ZMTB(IISTP,ISTATE)= 0.0d0
              END IF
          ENDDO
c-----------------------------------------------------------------------
      RETURN
  600 FORMAT(/'For state ',A2,' perform a forward calculation from a pot
     1ential defined using PREPOT')
  604 FORMAT(/' For state ',A2,' represent level energies by independent
     1 term values')
  626 FORMAT(/' For state   ',A2/4x,'integrate from   RMIN=',f5.2,
     1  '   to   RMAX=',f6.2,'   with mesh   RH=',f8.5)
  628 FORMAT(' ** Since  p=',i2,' .LE. [MMLR(NCMM)-MMLR(1)]=',i2,
     1  '  or   [MMLR(NCMM)-MMLR(1)].le.0,  reset   NCMM= 1')
  640 FORMAT(/' ', A2,' state energies referenced to f-parity levels')
  642 FORMAT(/' ', A2,' state energies referenced to the mid-point betwe
     1en e and f-parity levels')
  644 FORMAT(/' ', A2,' state energies referenced to e-parity levels')
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE WRITEPOT(NPASS,NSTATES,SLABL,NAME,DECM,PU,PS)
c***********************************************************************
c** Subroutine to print out complete description of the potential fx.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++                   Version of  12 January 2007
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** On entry:
c    NPASS   is the number of times WRIDATA was called.
c    NSTATES is the number of states being considered.
c    NAME    is the name of the molecule.
c    PU      are the parameter uncertainties.
c    PS      are the parameter sensitivities.
c=======================================================================
      INCLUDE 'arrsizes.h'
      INCLUDE 'BLKPOT.h'
      INCLUDE 'BLKBOB.h'
      INCLUDE 'BLKCOUNT.h'
c-----------------------------------------------------------------------
      CHARACTER*2  NAME(2),SLABL(-3:NSTATEMX)
      CHARACTER*3  MLX
      CHARACTER*7  NAMEDBLE(2)
      INTEGER NPASS, NSTATES,ISTATE,IPV, I,J, MMN, m, mpow,npow
      REAL*8 DECM(NSTATEMX), BTEMP, UAT, UBT, BINF, RE3,T1,ULRe
c
c** NLLSSRR variables used in output
c
      REAL*8 PU(NPARMX), PS(NPARMX)
      DATA NAMEDBLE/'wLambda',' wSigma'/
c-----------------------------------------------------------------------
c** Writing out state specific information.
c-----------------------------------------------------------------------
      IPV= 0
      DO  90 ISTATE=1,NSTATES
          WRITE(6,600)
          IF(PSEL(ISTATE).EQ.-2) THEN
              WRITE(6,601) SLABL(ISTATE)
              GO TO 90
              ENDIF
          IF(PSEL(ISTATE).EQ.1) WRITE(6,602) SLABL(ISTATE),nPB(ISTATE),
     1                                     NSphi(ISTATE),NLphi(ISTATE)
          IF(PSEL(ISTATE).EQ.2) THEN
              MLX= 'MLJ'
              IF(NCMM(ISTATE).GT.1) MLX= 'MLR'
              IF((NCMM(ISTATE).EQ.4).AND.(MMLR(2,ISTATE).EQ.0)) THEN
c** Fix up long range uLR for Aubert-Frecon {3,ASO,6,6} case ...
                  RE3= RE(ISTATE)**3
                  T1= CmVAL(1,ISTATE)/(9.d0*RE3) 
     1       + CmVAL(3,ISTATE)*(5.d0 + CmVAL(4,ISTATE))/(45.d0*RE3**2)
                  ULRe= 0.5d0*CmVAL(1,ISTATE)/RE3
     1 + CmVAL(3,ISTATE)*(5.d0 + 8.2d0*CmVAL(4,ISTATE))/(18.d0*RE3**2)
     2            + 0.5d0*DSQRT((T1- CmVAL(2,ISTATE))**2 + 8.d0*T1**2)
                ELSE
c** Fix up long range for 'normal' MLR/MLJ potential
                  ULRe= 0.d0
                  DO  m= 1,NCMM(ISTATE)
                      ULRe= ULRe+ CmVAL(m,ISTATE)/
     1                                      RE(ISTATE)**MMLR(m,ISTATE)
                      ENDDO
                  MMN= MMLR(NCMM(ISTATE),ISTATE) - MMLR(1,ISTATE)
                ENDIF
              BINF= DLOG(2.d0*DE(ISTATE)/ULRe)
              WRITE(6,604) SLABL(ISTATE),MLX,MMLR(1,ISTATE),nPB(ISTATE),
     1                NSphi(ISTATE),NLphi(ISTATE),BINF,MMLR(1,ISTATE),
     2                                  CmVAL(1,ISTATE),MMLR(1,ISTATE)
              IF(NCMM(ISTATE).GT.1) THEN
                  DO  m= 2,NCMM(ISTATE)
                      IF(MMLR(m,ISTATE).LE.9) 
     1      WRITE(6,608) MMLR(m,ISTATE),CmVAL(m,ISTATE),MMLR(m,ISTATE)
                      IF(MMLR(m,ISTATE).GT.9) 
     1      WRITE(6,609) MMLR(m,ISTATE),CmVAL(m,ISTATE),MMLR(m,ISTATE)
                      ENDDO
                  ENDIF
              IF((NCMM(ISTATE).GE.4).AND.(MMLR(2,ISTATE).EQ.0))
     1                          WRITE(6,606) (CmVAL(2,ISTATE), m= 2,4)
              IF(ASW(ISTATE).LE.0.d0) WRITE(6,605)
              IF(ASW(ISTATE).GT.0.d0)WRITE(6,607)ASW(ISTATE),RSW(ISTATE)
              ENDIF
          IF(PSEL(ISTATE).EQ.3) WRITE(6,612) SLABL(ISTATE),nPB(ISTATE),
     1                                     NSphi(ISTATE),NLphi(ISTATE)
          IF(PSEL(ISTATE).LE.3) THEN
              IF(RREF(ISTATE).LE.0.d0) WRITE(6,553) nPB(ISTATE),
     1                 nPB(ISTATE),nPB(ISTATE),nPB(ISTATE),nPB(ISTATE)
              IF(RREF(ISTATE).GT.0.d0) WRITE(6,555) nPB(ISTATE),
     1             nPB(ISTATE),RREF(ISTATE),nPB(ISTATE),nPB(ISTATE),
     2                                      RREF(ISTATE),nPB(ISTATE)
              ENDIF
          IF(PSEL(ISTATE).EQ.3) THEN
              IF(IDF(ISTATE).EQ.1) WRITE(6,624) RHOd(ISTATE)
              IF(IDF(ISTATE).EQ.2) WRITE(6,626) RHOd(ISTATE)
              WRITE(6,628) NCMM(ISTATE),MMLR(1,ISTATE),CmVAL(1,ISTATE)
              IF(NCMM(ISTATE).GT.1) WRITE(6,629) (MMLR(J,ISTATE),
     1                             CmVAL(j,ISTATE),j= 2,NCMM(ISTATE))
              ENDIF
          IF(PSEL(ISTATE).EQ.4) THEN
              WRITE(6,610) SLABL(ISTATE),(nPB(ISTATE),i=1,3),
     1             AGPEF(ISTATE),nPB(ISTATE),BGPEF(ISTATE),nPB(ISTATE)
              ENDIF
c
          IF((NUA(ISTATE).GE.0).OR.(NUB(ISTATE).GE.0)) THEN
              IF(BOBCN(ISTATE).LE.0) THEN
                  npow= nPu(ISTATE)
                  mpow= npow
                  IF(PSEL(ISTATE).EQ.2) mpow= MMLR(1,ISTATE)
                  WRITE(6,557) 
                  IF(NUA(ISTATE).GE.0)
     1               WRITE(6,558) NAME(1),mpow,mpow,npow,NUA(ISTATE)-1
                  IF(NUB(ISTATE).GE.0)
     1               WRITE(6,558) NAME(2),mpow,mpow,npow,NUB(ISTATE)-1
                ELSE
                  WRITE(6,564) (nPu(ISTATE),i=1,5)
                ENDIF
              ENDIF
          IF((NTA(ISTATE).GE.0).OR.(NTB(ISTATE).GE.0)) THEN
              IF(BOBCN(ISTATE).LE.0) THEN
                  WRITE(6,559) 
                  IF(NTA(ISTATE).GE.0) WRITE(6,560) NAME(1),nPt(ISTATE),
     1                           nPt(ISTATE),nPt(ISTATE),NTA(ISTATE)-1
                  IF(NTB(ISTATE).GE.0) WRITE(6,560) NAME(2),nPt(ISTATE),
     1                           nPt(ISTATE),nPt(ISTATE),NTB(ISTATE)-1
                ELSE
                  WRITE(6,562) (nPt(ISTATE),i=1,5)
                ENDIF
              ENDIF
          IF(NwCFT(ISTATE).GE.0) THEN
              IF(IOMEG(ISTATE).GT.0) WRITE(6,618) 'Lambda',
     1         (nPW(ISTATE),i=1,3),NwCFT(ISTATE)-1,(nPW(ISTATE),i=1,5)
              IF(IOMEG(ISTATE).LT.0) WRITE(6,618) ' Gamma',
     1         (nPW(ISTATE),i=1,3),NwCFT(ISTATE)-1,(nPW(ISTATE),i=1,5)
              ENDIF
c
c** Write out potential constraints for each state (if any).
c
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
              WRITE(6,620) 'Te',UAT,UBT,0.0d0
              END IF
c-----------------------------------------------------------------------
c** Writing out the De information.
c-----------------------------------------------------------------------
          IF(PSEL(ISTATE).NE.4) THEN
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
              ENDIF
c-----------------------------------------------------------------------
c** Writing out the Re information.
c-----------------------------------------------------------------------
          IPV= IPV + 1
          IF(IFXRE(ISTATE).LE.0) THEN
              IF (DABS(RE(ISTATE)).GT.PU(IPV)) THEN
                  WRITE(6,620) 'Re',RE(ISTATE),PU(IPV),PS(IPV)
                ELSE
                  WRITE(6,621) 'Re',RE(ISTATE),PU(IPV),PS(IPV)
                ENDIF
            ELSE
              WRITE(6,622) 'Re',RE(ISTATE)
            ENDIF
c
          IF(PSEL(ISTATE).EQ.2) THEN
c-----------------------------------------------------------------------
c** For MLR/MLJ, write out the  Cn  information.
c-----------------------------------------------------------------------
              IF((NCMM(ISTATE).EQ.4).AND.(MMLR(2,ISTATE).EQ.0)) THEN
c ... For Aubert-Frecon treatment of C3(r):C6(r) for alkali dimers
                  IPV= IPV+ 1
                  IF(IFXCm(1,ISTATE).LE.0) THEN
                      IF(DABS(CmVAL(1,ISTATE)).GT.PU(IPV))
     1                  WRITE(6,620)'M2',CmVAL(1,ISTATE),PU(IPV),PS(IPV)
                      IF(DABS(CmVAL(1,ISTATE)).LE.PU(IPV))
     1                  WRITE(6,621)'M2',CmVAL(1,ISTATE),PU(IPV),PS(IPV)
                    ELSE
                      WRITE(6,622) 'M2',CmVAL(1,ISTATE)
                    ENDIF
c** For Aubert-Frecon treatment of C3(r):C6(r) for alkali dimers
                  IPV= IPV+1
                  IF(IFXCm(3,ISTATE).LE.0) THEN
                      IF(DABS(CmVAL(3,ISTATE)).GT.PU(IPV)) THEN
                          WRITE(6,630) 'C6_01',CmVAL(3,ISTATE),
     1                                                 PU(IPV),PS(IPV)
                        ELSE
                          WRITE(6,631) 'C6_01',CmVAL(3,ISTATE),
     1                                                 PU(IPV),PS(IPV)
                        ENDIF
                    ELSE
                      WRITE(6,632) 'C6_01',CmVAL(3,ISTATE)
                    ENDIF
                  IPV= IPV+1
                  IF(IFXCm(4,ISTATE).LE.0) THEN
                      IF(DABS(CmVAL(4,ISTATE)).GT.PU(IPV)) THEN
                          WRITE(6,630) 'Q12C6',CmVAL(4,ISTATE),
     1                                                 PU(IPV),PS(IPV)
                        ELSE
                          WRITE(6,631) 'Q12C6',CmVAL(4,ISTATE),
     1                                                 PU(IPV),PS(IPV)
                        ENDIF
                    ELSE
                      WRITE(6,632) 'Q12C6',CmVAL(4,ISTATE)
                    ENDIF
                ELSE
c ... For 'regular' MLJ or MLR case ...
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
          DO  I=0, NLphi(ISTATE)
              IPV= IPV + 1
              IF(IFXPHI(I,ISTATE).LE.0) THEN
                  IF(DABS(PHI(I,ISTATE)).GT.PU(IPV)) THEN
                      WRITE(6,640) ' p','hi',I,PHI(I,ISTATE),PU(IPV),
     1                                                         PS(IPV)
                    ELSE
                      WRITE(6,641) ' p','hi',I,PHI(I,ISTATE),PU(IPV),
     1                                                         PS(IPV)
                    ENDIF
                ELSE
                  WRITE(6,640) ' p','hi',I,PHI(I,ISTATE)
                ENDIF
              BTEMP= BTEMP+ PHI(I,ISTATE)
              ENDDO
          IF(PSEL(ISTATE).EQ.1) BINF= BTEMP
c** Write out  phi_\infty  constant for the MLR/MLJ form
          IF(PSEL(ISTATE).EQ.2) THEN
               WRITE(6,648) BINF
               BTEMP= CmVAL(1,ISTATE)*2.d0*(2.d0*BINF - BTEMP)
     1                                        *RE(ISTATE)**nPB(ISTATE)
               WRITE(6,652) MMLR(1,ISTATE)+nPB(ISTATE),BTEMP
               ENDIF
c-----------------------------------------------------------------------
c** Writing out the adiabatic BOB radial function for atom A.
c-----------------------------------------------------------------------
          IF(NUA(ISTATE).GE.0) THEN
              DO  I= 0,NUA(ISTATE)-1
                  IPV= IPV + 1
                  IF(IFXUA(I,ISTATE).LE.0) THEN
                      IF (DABS(UA(I,ISTATE)).GT.PU(IPV)) THEN
                          WRITE(6,640) ' u',NAME(1),I,UA(I,ISTATE),
     1                                               PU(IPV),PS(IPV)
                        ELSE
                          WRITE(6,641) ' u',NAME(1),I,UA(I,ISTATE),
     1                                               PU(IPV),PS(IPV)
                        ENDIF
                    ELSE
                      WRITE(6,650) ' u',NAME(1),I,UA(I,ISTATE)
                    ENDIF
                  END DO
              IPV= IPV + 1
              IF(IFXUA(NUA(ISTATE),ISTATE).LE.0) THEN
                  WRITE(6,644) ' u',NAME(1),UA(NUA(ISTATE),ISTATE),
     1                                                 PU(IPV),PS(IPV)
                ELSE
                  WRITE(6,646) ' u',NAME(1),UA(NUA(ISTATE),ISTATE)
                ENDIF
              ENDIF    
c-----------------------------------------------------------------------
c** Writing out the adiabatic BOB radial function for atom B.
c-----------------------------------------------------------------------
        IF(NUB(ISTATE).GE.0) THEN
            DO  I= 0, NUB(ISTATE)- 1
                IPV= IPV + 1
                IF(IFXUB(I,ISTATE).EQ.0) THEN
                    IF(DABS(UB(I,ISTATE)).GT.PU(IPV)) THEN
                        WRITE(6,640) ' u',NAME(2),I,UB(I,ISTATE),
     1                                             PU(IPV),PS(IPV)
                      ELSE
                        WRITE(6,641) ' u',NAME(2),I,UB(I,ISTATE),
     1                                             PU(IPV),PS(IPV)
                      ENDIF
                  ELSE
                    WRITE(6,650) ' u',NAME(2),I,UB(I,ISTATE)
                  ENDIF
                END DO
            IPV= IPV + 1
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
                END DO
            IPV= IPV + 1
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
                END DO
            IPV= IPV + 1
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
            J= 1
            IF(IOMEG(ISTATE).LT.0) J=2
            DO  I= 0, NwCFT(ISTATE)-1
                IPV= IPV + 1
                IF(IFXwCFT(I,ISTATE).LE.0) THEN
                    IF(DABS(wCFT(I,ISTATE)).GT.PU(IPV)) THEN
                        WRITE(6,642) NAMEDBLE(J),I,wCFT(I,ISTATE),
     1                                             PU(IPV),PS(IPV)
                      ELSE
                        WRITE(6,643) NAMEDBLE(J),I,wCFT(I,ISTATE),
     1                                             PU(IPV),PS(IPV)
                      ENDIF
                  ELSE         
                    WRITE(6,651) NAMEDBLE(J),I,wCFT(I,ISTATE)
                  ENDIF
                END DO
            IPV= IPV + 1
            IF(IFXwCFT(NwCFT(ISTATE),ISTATE).LE.0) THEN
                IF(DABS(wCFT(NwCFT(ISTATE),ISTATE)).GT.PU(IPV)) THEN
                    WRITE(6,654) NAMEDBLE(J),wCFT(NwCFT(ISTATE),ISTATE),
     1                                             PU(IPV),PS(IPV)
                  ELSE
                    WRITE(6,656) NAMEDBLE(J),wCFT(NwCFT(ISTATE),ISTATE),
     1                                             PU(IPV),PS(IPV)
                  ENDIF
              ELSE
                WRITE(6,658) NAMEDBLE(J),wCFT(NwCFT(ISTATE),ISTATE)
              ENDIF
            ENDIF
   90     CONTINUE
      WRITE(6,600)
      RETURN
c-----------------------------------------------------------------------
  553 FORMAT(4x,'using radial expansion variable:   yp = y',i1,' = (R^',
     1  i1,' - Re^',i1,')/(R^',i1,' + Re^',i1,')')
  555 FORMAT(4x,'using radial expansion variable:   yp = y',i1,
     1  ' = (R^',i1,' -',F5.2,'^',i1,')/(R^',i1,' +',f5.2,'^',i1,')')
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
     1 '(r) = (r^',i1,' - re^',i1,')/(r^',i1,' + re^',i1,')' )
  600 FORMAT(1X,39('=='))
  601 FORMAT(' Distinct levels of State ',A2,' fitted as independent ter
     1m values')
  602 FORMAT(/' State ',A2,' is represented by an Expanded Morse Oscilla
     1tor EMO(p=',i2,') potential'/1x,4('=='),'  with exponent  power se
     2ries order',i3,' for  R < Re','  and',i3,' for  R > Re')
  604 FORMAT(/' State ',A2,' represented by an ',A3,'(n=',i2,';p=',i2,
     1  ') potential defined in terms of an'/1x,4('=='),'  exponent powe
     2r series of order',i3,' for  R < Re','  and',i3,' for  R > Re'/
     3 7x,'with   phiINF=',F12.8,'  defined by:   C',I1,'=',1PD12.5,
     4 '[cm-1 Ang^',i1,']')
  608 FORMAT(49x,'C',I1,'=',1PD12.5,'[cm-1 Ang^',i1,']')
  609 FORMAT(48x,'C',I2,'=',1PD12.5,'[cm-1 Ang^{',i2,'}]')
  605 FORMAT(4x,'& exponent coefficient given by:   phi(R)= phiINF*yp +
     1 (1-yp)*Sum{phi*yp^i}')
  607 FORMAT(4x,'& exponent based on switching function with  Asw=',
     1 F8.6,'  &  Rsw=',F9.6)
  606 FORMAT(7x,'where Aubert-Frecon ULR(r) has   ASO=',f10.6,
     1  '   C6_01=',1PD16.8/57x,'Q12C6=',D16.8)
  610 FORMAT(/' For state ',A2,"  use Surkus' Generalized Potential Ener
     1gy Function GPEF with"/1x,6('=='),'  expansion vble:   y_',i1,
     2 '(r) = (r^',i1,' - re^',i1,')/(',F5.2,'*r^',i1,' +',F5.2,'*re^',
     3 i1,'p)')
  612 FORMAT(/' State ',A2,' represented by a DELR(p=',i2,') potential w
     1ith an exponent'/1x,4('=='),'  coefficient power series of order',
     2  i3,' for  R < Re','  and',i3,' for  R >= Re') 
  614 FORMAT(/'  Parameter:   Initial Value:   Uncertainty:  Sensitivity
     1:')
  615 FORMAT(/'  Parameter:    Final Value:    Uncertainty:  Sensitivity
     1:')
  618 FORMAT(1x,A6,'-doubling splitting strength function is expanded as
     1'/7x,'f(r) =  w(inf)*y',i1,' + (1 - y',i1,')*Sum{w_i * (y',i1,')^i
     2} for  i= 0 to',i3 /5x,'where   y',i1,'= (R^',i1,' - Re^',i1,
     3 ')/(R^',i1,' + Re^',i1,')' )
  620 FORMAT(6X,A2,3X,1PD19.10,3X,1PD8.1,6X,1PD8.1)
  621 FORMAT(5X,'*',A2,3X,1PD19.10,3X,1PD8.1,6X,1PD8.1)
  622 FORMAT(6X,A2,3X,1PD19.10,7X,'--',12X,'--')
  624 FORMAT(3x,'U_{LR} term uses Tang-Toennies damping function [JCP 80
     1,3726(1984)]'/5x,'[1 - exp(-RHOd*r)*SUM{(RHOd*r)^k/k!}]   with   R
     2HOd =',F9.6)
  626 FORMAT(3x,'U_{LR} term uses Douketis dispersion damping function [
     1Mol.P. 52,763(1984)]'/5x,'[1 - exp(-3.97*RHOd*r/m - 0.39*(RHOd*r)^
     22/sqrt{m})]^m  with  RHOd=',F9.6)
  628 FORMAT('   and',i2,' inverse-powers terms with coefficients:',
     1  '   C',0Pi2,'=',1Pd14.6)
  629 FORMAT((51x,'C',0Pi2,'=',1Pd14.6:))
  630 FORMAT(3X,A5,3X,1PD19.10,3X,1PD8.1,6X,1PD8.1)
  631 FORMAT(2X,'*',A5,3X,1PD19.10,3X,1PD8.1,6X,1PD8.1)
  632 FORMAT(3X,A5,3X,1PD19.10,7X,'--',12X,'--')
  633 FORMAT(4x,A7,1PD19.10)
  636 FORMAT(4X,A4,3X,1PD19.10,7X,'--',12X,'--')
  640 FORMAT(3X,A2,A2,'(',I2,')',1PD19.10,3X,1PD8.1,6X,1PD8.1)
  641 FORMAT('  *',2A2,'(',I2,')',1PD19.10,3X,1PD8.1,6X,1PD8.1)
  642 FORMAT(1X,A6,'(',I2,')',1PD19.10,3X,1PD8.1,6X,1PD8.1)
  643 FORMAT('*',A6,'(',I2,')',1PD19.10,3X,1PD8.1,6X,1PD8.1)
  644 FORMAT(1x,A2,'_inf(',A2,')',1PD19.10,1PD11.1,6X,1PD8.1)
  646 FORMAT(1x,A2,'_inf(',A2,')',1PD19.10,7X,'--',12X,'--')
  648 FORMAT(4X,'phi_INF',1PD19.10,7X,'--',12X,'--')
  650 FORMAT(3X,A2,A2,'(',I2,')',1PD19.10,7X,'--',12X,'--')
  651 FORMAT(1X,A6,'(',I2,')',1PD19.10,7X,'--',12X,'--')
  652 FORMAT('   C',I2,'{eff}',1PD19.10,7X,'--',12X,'--')
  654 FORMAT(1X,A6,'_inf',1PD19.10,3X,1PD8.1,6X,1PD8.1)
  656 FORMAT('*',A6,'_inf',1PD19.10,3X,1PD8.1,6X,1PD8.1)
  658 FORMAT(1X,A6,'_inf',1PD19.10,7X,'--',12X,'--')
  660 FORMAT(5X,'C',I2,3X,1PD19.10,3X,1PD8.1,6X,1PD8.1)
  661 FORMAT(3X,'* C',I2,3X,1PD19.10,3X,1PD8.1,6X,1PD8.1)
  662 FORMAT(5X,'C',I2,3X,1PD19.10,7X,'--',12X,'--':)
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE VGEN(ISTATE,RDIST,VDIST,PHIDIST,IDAT)
c***********************************************************************
c** This subroutine will generate one of six possible families of
c   analytical molecular potentials for the direct hamiltonian fitting
c   program. These four families are:
c     - The Extended Morse Oscillator
c     - The Modified Lennard-Jones Oscillator
c     - The Double-Exponential/Long-Range potential
c     - The Generalized Potential Energy Function
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++  COPYRIGHT 2006  by  R.J. Le Roy, Jenning Seto & Yiye Huang     +
c   Dept. of Chemistry, Univ. of Waterloo, Waterloo, Ontario, Canada   +
c    This software may not be sold or any other commercial use made    +
c      of it without the express written permission of the authors.    +
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c              ----- Version of 12 January 2007 -----
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** On entry:
c   ISTATE  is the electronic state being considered in this CALL.
c   RDIST: If RDIST > 0, calculate potl & derivs only @ that onee distance
c   -----    * return potential function at that point as VDIST and
c              potential function exponent as PHIDIST
c            * skip partial derivative calculation if  IDAT.le.0
c        * If RDIST.le.0   use {RMIN,RMAX,RH} to create array of values,
c          calculate partial derivatives & return them in array DVtot
c** On entry via common blocks:
c PSEL   specifies the type of potential used.
c      PSEL = 0 : Do forward calculation from input pointwise potential
c      PSEL = 1 : Use the Expanded Morse Oscillator(p) potential.
c      PSEL = 2 : Use the Morse/Lennard-Jones(p) potential.
c      PSEL = 3 : Use the Double-Exponential Long-Range Potential.
c      PSEL = 4 : Use the Generalized Potential Energy Function.
c  NSphi(s)  is the order of the  phi(r)  exponent expansion for R.le.Re
c  NLphi(s)  is the order of the  phi(r)  exponent expansion for R > Re
c!!!  The code (currently) assumes  NLphi .ge. NSphi !!!!
c    MMLR(j,s)  are long-range inverse-powers for an MLR or DELR potential
c    nPB(s)  the value of power p for  phi(yp)  exponent function
c    nPu(s)  the value of power p for adiabatic u(r) BOB functions
c    nPt(s)  the value of power p for centrifugal q(r) BOB functions
c    nPw(s)  the value of power m for Lambda(r) function
c    DE     is the Dissociation Energy for each state.
c    RE     is the Equilibrium Distance for each state.
c    PHI    is the array of potential (exponent) expansion parameters 
c    NDATPT is the number of meshpoints used for the array.
c-----------------------------------------------------------------------
c** On exit via common blocks:
c    R      is the distance array
c    VPOT   is the potential that is generated.
c    PHIFX  is used to contain the  phi(r)  function.
c** Internal partial derivative arrays ...
c    DUADRe & DUBDRe are p.derivs of adiabatic fx. w.r.t.  Re
c    DVDQA & DVDQB are p.derivs of non-adiabatic fx. wrt q_A(i) & q_B(i)
c    DTADRe & DTBDRe are p.derivs of non-adiabatic fx. w.r.t.  Re
c    dVdL & dLDDRe are p.derivatives of f_\lambda(r) w.r.t. phi_i & Re
c    DBDB & DBDRe are p.derives of phi(r) w.r.t. \phi_i & Re, respectively
c
c** Temp:
c    BTEMP     is used to represent the sum used for dV/dRe.
c              is used in GPEF for De calculations.
c    BINF      is used to represent the  phi(\infty)  value.
c    YP        is used to represent (R-Re)/(R+Re) or R-Re.
c                                         n    -PHIFX*RTEMP
c    XTEMP     is used to represent (Re/R)  * e 
c                                         n    -BINF*RTEMP
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
c    Rsw(ISTATE) is the Rx value used in the MLJ switching function.
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
      INTEGER I,J,I1,ISTATE,IPV,IPVSTART,ISTART,ISTOP,LAMB2,m,MMN,npow,
     1  IDAT, NBAND, IISTP, ICOXON
      REAL*8 BTEMP,BINF,RVAL,RTEMP,RM2,XTEMP,PBTEMP,PETEMP,
     1  FSW,Xtemp2,Btemp2,BMtemp,BMtemp2,RMF,PBtemp2,
     2  YP,YPE,YPM,YPP,YM,YMM,REP,REm,RDp,RDm,DYPDRE,VAL,DVAL,
     3  HRe,HRem,SL,SLB, AREF,AREFp, RE3,T0,T1,ULRe,
     4  dLULRedCm(9),dLULRedRe,dULRdCm(9),RD3, DVDD, RDIST,VDIST,
     5  PHIDIST,BFCT,JFCT,JFCTLD
c***********************************************************************
c** Common block for partial derivatives at the one distance RDIST
      REAL*8 dVdPk(HPARMX)
      COMMON /dVdPkBLK/dVdPk
c***********************************************************************
c** Temporary variables for DELR potentials
      INTEGER MMLRP
      REAL*8 ULR,Dm,Dmp,Dmpp,dAAdRe,dBBdRe,dVdBtemp,CmVALL
c***********************************************************************
c** Initializing variables.
      REP= RE(ISTATE)**nPB(ISTATE)
      IF(RREF(ISTATE).LE.0) AREF= RE(ISTATE)
      IF(RREF(ISTATE).GT.0) AREF= RREF(ISTATE)
      AREFp= AREF**nPB(ISTATE)
c** Normally data point starts from 1
      ISTART= 1
      ISTOP= NDATPT(ISTATE)
c** When calculating only one potential point
      IF(RDIST.GT.0) THEN
          ISTART= NPNTMX
          ISTOP= NPNTMX
          VDIST= 0.0d0
          PHIDIST= 0.d0
          ENDIF
      PBTEMP= 0.0d0
      PETEMP= 0.0d0
      DO  I= ISTART,ISTOP
          PHIFX(I,ISTATE)= 0.0d0
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
c** For forward data simulation using a known input potential 
c?????? Most/all of this section is dummy ?!?!?!?!?!?!?!?!?!/
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          DO  I= ISTART,ISTOP
              RVAL= RMIN(ISTATE) + RH(ISTATE)*DBLE(I-1)
              RD(I,ISTATE)= RVAL
              RTEMP= (RVAL**nPB(ISTATE)-REP)/(RVAL**nPB(ISTATE)+REP)
              PHIFX(I,ISTATE)= PHI(0,ISTATE)
              DO  J= 1,NLphi(ISTATE)
                  PHIFX(I,ISTATE)=PHIFX(I,ISTATE)+PHI(J,ISTATE)*RTEMP**J
                  ENDDO
              XTEMP= DEXP(-PHIFX(I,ISTATE)*(RVAL-RE(ISTATE)))
c** Now to calculate the actual potential
              VPOT(I,ISTATE)= DE(ISTATE)*(1.0d0-XTEMP)**2
     1                   + VLIM(ISTATE) - DE(ISTATE)
              VDIST= VPOT(I,ISTATE)
              ENDDO
          ENDIF
******7***** End Forward Calculation ***********************************

      IF(PSEL(ISTATE).EQ.1) THEN
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** For the Expanded Morse Oscillator - with  NSphi/NLphi
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** First ... calculate the Generalized Morse Oscillator exponent
          DO  I= ISTART,ISTOP
              RVAL= RDIST
              IF(RDIST.LE.0.d0) THEN
                  RVAL= RMIN(ISTATE) + RH(ISTATE)*DBLE(I-1)
                  ENDIF
              RD(I,ISTATE)= RVAL
              RDp= RVAL**nPB(ISTATE)
              YP= (RDp - AREFp)/(RDp + AREFp)
              VAL= PHI(0,ISTATE) 
              DVAL= 0.d0
              DBDB(0,I,ISTATE)= 1.0d0
              YPP= 1.d0
              npow= NSphi(ISTATE)
              IF(RVAL.GT.RE(ISTATE)) npow= NLphi(ISTATE)
              DO  J= 1,npow
                  DVAL= DVAL + PHI(J,ISTATE)* DBLE(J)* YPP
                  YPP= YPP*YP
                  VAL= VAL + PHI(J,ISTATE)*YPP 
                  DBDB(J,I,ISTATE)= YPP
                  ENDDO
c*** DBDB & DBDRe= dBeta/dRe  used in uncertainty calculation in fununc.f
              DBDRe(I,ISTATE)= 0.d0
              PHIFX(I,ISTATE)= VAL
              XTEMP= DEXP(-VAL*(RVAL-RE(ISTATE)))
c** First calculate the partial derivative w.r.t. DE
              IPV= IPVSTART+ 1
              DVtot(IPV,I)= XTEMP*(XTEMP- 2.d0)
              DVDD= DVtot(IPV,I)
c** Now calculate the actual potential
              VPOT(I,ISTATE)= DE(ISTATE)*DVDD + VLIM(ISTATE)
              IF(RDIST.GT.0.d0) THEN
                  VDIST= VPOT(I,ISTATE)
                  PHIDIST= VAL
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
              DO  J= 0,npow
                  IPV= IPV+1
                  DVtot(IPV,I)= YPP
                  YPP= YPP*YP
                  ENDDO
              ENDDO
          ENDIF
ccccc Print for testing
cc      write(10,610) (RD(i,ISTATE),vpot(i,istate),PHIFX(i,istate),
cc   1                                        i= 1, NDATPT(ISTATE),50)
ccccc
c******** End preparation of Expanded Morse Potential Function *********

      IF(PSEL(ISTATE).EQ.2) THEN
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** For the Morse/Lennard-Jones(p)  or  Morse/Long-Range  potential.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c... if necessary ... set flag for Coxon use of 'true' O-T variable
          ICOXON= 0
          IF((nPB(ISTATE).EQ.1).AND.(BOBCN(ISTATE).GE.1)) ICOXON= 1
          MMN= MMLR(NCMM(ISTATE),ISTATE) - MMLR(1,ISTATE)
c** First - define values & derivatives of uLR at Re
          IF((NCMM(ISTATE).EQ.4).AND.(MMLR(2,ISTATE).EQ.0)) THEN
c ... for Aubert-Frecon treatment of C3/C6 for Alkali A-state
              MMN= 9
              RE3= RE(ISTATE)**3
              T1= CmVAL(1,ISTATE)/(9.d0*RE3)
     1       + CmVAL(3,ISTATE)*(5.d0 + CmVAL(4,ISTATE))/(45.d0*RE3**2)
              T0= DSQRT((T1- CmVAL(2,ISTATE))**2 + 8.d0*T1**2)
              ULRe= 0.5d0*(CmVAL(1,ISTATE)/RE3 - CmVAL(2,ISTATE)) +
     1    CmVAL(3,ISTATE)*(5.d0+ 8.2d0*CmVAL(4,ISTATE))/(18.d0*RE3**2)
     2                                                      + 0.5d0*T0
              T0= (9.d0*T1-CmVAL(2,ISTATE))/T0
              dLULRedCm(1)= (9.d0 + T0)/(18.D0*RE3 * ULRe)
              dLULRedCm(3)= (25.d0 + 41.d0*CmVAL(4,ISTATE) + 
     1               T0*(5.d0+ CmVAL(4,ISTATE)))/(90.D0*RE3**2 * ULRe)
              dLULRedCm(4)= CmVAL(3,ISTATE)*(41.d0+ T0)/
     1                                             (90.d0*RE3**2*ULRe)
              dLULRedRe= (-1.5d0*CmVAL(1,ISTATE)
     1     - CmVAL(3,ISTATE)*(5.d0 + 8.2d0*CmVAL(4,ISTATE))/(3.d0*RE3)
     2     - 0.5d0*T0*(CmVAL(1,ISTATE)/3.d0 + CmVAL(3,ISTATE)*(10.d0
     3   + 2.d0*CmVAL(4,ISTATE))/(15.d0*Re3)))/(Re3*RE(ISTATE) * ULRe)
            ELSE
c ... for ordinary NCMM-term MLR uLR(r)
              ULRe= 0.d0
              T1= 0.d0
              DO  m= 1,NCMM(ISTATE)
                  dLULRedCm(m)= 1.d0/RE(ISTATE)**MMLR(m,ISTATE)
                  T0= CmVAL(m,ISTATE)*dLULRedCm(m)
                  ULRe= ULRe + T0
                  T1= T1 + MMLR(m,ISTATE)*T0
                  ENDDO
              DO  m= 1,NCMM(ISTATE)
                  dLULRedCm(m)= dLULRedCm(m)/ULRe
                  ENDDO
              dLULRedRe= -T1/(ULRe*RE(ISTATE))
            ENDIF
          BINF= DLOG(2.0d0*DE(ISTATE)/ULRe)
          IF(ICOXON.GE.1) BINF= 0.5d0*BINF
          DO  I= ISTART,ISTOP
              RVAL= RDIST
              IF(RDIST.LE.0.d0) THEN
                  RVAL= RMIN(ISTATE) + RH(ISTATE)*DBLE(I-1)
                  ENDIF
              RD(I,ISTATE)= RVAL
              RDp= RVAL**nPB(ISTATE)
              YPE= (RDp-REP)/(RDp+REP)
              YP= (RDp-AREFp)/(RDp+AREFp)
              IF(ICOXON.GE.1) THEN
                  YP= 2.d0*YP
                  YPE= 2.d0*YPE
                  ENDIF
              YPM= 1.d0 - YP
              DYPDRE= -0.5d0*nPB(ISTATE)*(1.d0 - YP**2)/RE(ISTATE)
              IF(ICOXON.GE.1) DYPDRE= 2.d0*DYPDRE
              YPP= 1.d0
              VAL= PHI(0,ISTATE)
              DVAL= 0.d0
              DBDB(0,I,ISTATE)= 1.0d0
              npow= NSphi(ISTATE)
              IF(RVAL.GT.RE(ISTATE)) npow= NLphi(ISTATE)
              DO  J= 1,npow
c ... now calculate power series part of the Morse-like exponent
                  DVAL= DVAL + PHI(J,ISTATE)* DBLE(J)* YPP
                  YPP= YPP*YP
                  VAL= VAL + PHI(J,ISTATE)*YPP
                  DBDB(J,I,ISTATE)= YPP*YPM
                  ENDDO
c*** DBDB & DBDRe= dBeta/dRe  used in uncertainty calculation in fununc.f
              IF(Asw(ISTATE).GT.0.0d0) THEN
c-----------------------------------------------------------------------
c** For switching function version of constrained exponent ...
c           fsw(R) = 1 / [e^(Asw'*(R - Rx)) + 1]
c-----------------------------------------------------------------------
                  FSW= 1.0d0/(1.0d0 + 
     1                         DEXP(Asw(ISTATE)*(RVAL-Rsw(ISTATE))))
                  DBDRe(I,ISTATE)= (1.d0-FSW)*DBLE(MMLR(1,ISTATE))
     1                                                     /RE(ISTATE)
                  IF(RREF(ISTATE).LE.0.d0) DBDRe(I,ISTATE)=
     1                               DBDRe(I,ISTATE) + FSW*DVAL*DYPDRE
                  VAL= FSW*(VAL - BINF) + BINF
                ELSEIF(Asw(ISTATE).LE.0.0d0) THEN
c-----------------------------------------------------------------------
c** But preferably - constrain leading Cm value as per Huang thesis ...
c-----------------------------------------------------------------------
                  DBDRe(I,ISTATE)= -YP*dLULRedRe
                  IF(RREF(ISTATE).LE.0.d0) DBDRe(I,ISTATE)= 
     1           DBDRe(I,ISTATE)+ (BINF - VAL + (1.d0-YP)*DVAL)*DYPDRE
                  VAL= YP*BINF + (1.d0- YP)*VAL
                ENDIF
              PHIFX(I,ISTATE)= VAL
              XTEMP= DEXP(-VAL*YPE)
c** Now begin by generating  uLR(r)
              IF((NCMM(ISTATE).EQ.4).AND.(MMLR(2,ISTATE).EQ.0)) THEN
c ... generate ULR for Aubert-Frecon Li2(A) {3,0,6,6} type case ...
                  RD3= RVAL**3
                  T1= CmVAL(1,ISTATE)/(9.d0*RD3)
     1       + CmVAL(3,ISTATE)*(5.d0 + CmVAL(4,ISTATE))/(45.d0*RD3**2)
                  T0= DSQRT((T1- CmVAL(2,ISTATE))**2 + 8.d0*T1**2)
                  ULR= 0.5d0*(CmVAL(1,ISTATE)/RD3 - CmVAL(2,ISTATE)) +
     1    CmVAL(3,ISTATE)*(5.d0+ 8.2d0*CmVAL(4,ISTATE))/(18.d0*RD3**2)
     2                                                      + 0.5d0*T0
                  T0= (9.d0*T1- CmVAL(2,ISTATE))/T0
                  dULRdCm(1)= (9.d0 + T0)/(18.D0*RD3)
                  dULRdCm(2)= 0.d0
                  dULRdCm(3)= (25.d0+ 41.d0*CmVAL(4,ISTATE) 
     1                   + T0*(5.d0 + CmVAL(4,ISTATE)))/(90.D0*RD3**2)
                  dULRdCm(4)= CmVAL(3,ISTATE)*(41.d0+ T0)/(90.d0*RD3**2)
                ELSE
c ... for the case of a 'normal' MLR or MLJ function
                  ULR= 0.d0
                  DO  m= 1,NCMM(ISTATE)
                      dULRdCm(m)= 1.d0/RVAL**MMLR(m,ISTATE)
                      ULR= ULR + CmVAL(m,ISTATE)*dULRdCm(m)
                      ENDDO
                ENDIF
              XTEMP= XTEMP*ULR/ULRe
              DVDD= XTEMP*(XTEMP - 2.D0)  
              VPOT(I,ISTATE)= DE(ISTATE)*DVDD + VLIM(ISTATE)
              IF(RDIST.GT.0.d0) THEN
                  VDIST= VPOT(I,ISTATE)
                  PHIDIST= VAL
                  IF(IDAT.LE.0) GO TO 999
                  ENDIF
              YPP= 2.d0*DE(ISTATE)*(1.0d0-XTEMP)*XTEMP
              IPV= IPVSTART+2
              IF((NCMM(ISTATE).EQ.4).AND.(MMLR(2,ISTATE).EQ.0)) THEN
c ... derivatives w.r.t. long-range parameters for Aubert-Frecon  uLR
                  IPV= IPV+1
                  DVtot(IPV,I)= YPP*((1.d0 - YP*YPE)*dLULRedCm(1)
     1                                               - dULRdCm(1)/ULR)
                  IPV= IPV+1
                  DVtot(IPV,I)= YPP*((1.d0 - YP*YPE)*dLULRedCm(3)
     1                                               - dULRdCm(3)/ULR)
                  IPV= IPV+1
                  DVtot(IPV,I)= YPP*((1.d0 - YP*YPE)*dLULRedCm(4)
     1                                               - dULRdCm(4)/ULR)
                ELSE
c ... derivative w.r.t. Cm's for ordinary MLR/MLJ case ...
                  DO  m= 1, NCMM(ISTATE)
                      IPV= IPV+ 1
                      DVtot(IPV,I)= -YPP*YP*YPE*dLULRedCm(m)
                      ENDDO
                ENDIF
c... derivative w.r.t. Re  
              DVtot(IPVSTART+2,I)= YPP*(YPE*DBDRe(I,ISTATE)
     1                                       + VAL*DYPDRE + dLULRedRe)
              YPP= YPP*YPE
              IF(Asw(ISTATE).LE.0.d0) THEN
c... derivative w.r.t. De  for Huang-type exponent
                  DVDD= DVDD + YPP*YP/DE(ISTATE)
                  YPP= YPP*(1.d0 - YP)
                ELSE
c... derivative w.r.t. De  for switching fx. exponent
                  DVDD= DVDD + YPP*(1.d0- FSW)/DE(ISTATE)
                  YPP= YPP*FSW
                ENDIF
              DVtot(IPVSTART+1,I)= DVDD
c... finally ... derivatives w.r.t. exponent expansion coefficients
              DO  J= 0,npow
                  IPV= IPV+1
                  DVtot(IPV,I)= YPP
                  YPP= YPP*YP
                  ENDDO
              ENDDO
ccccc Print for testing
cc      rewind(10)
cc      write(10,610) (RD(i,ISTATE),vpot(i,istate),PHIFX(i,istate),
cc   1                                        i= 1, NDATPT(ISTATE),50)
cc610 FORMAT(/(f10.4,f15.5,f12.6))
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
c... first, get  AA & BB and their derivatives!
          DO  m= 1,NCMM(ISTATE)
              MMLRP= MMLR(m,ISTATE)
              CmVALL= CmVAL(m,ISTATE)
              CALL dampF(RE(ISTATE),RHOd(ISTATE),MMLRP,Dm,Dmp,
     1                                               Dmpp,IDF(ISTATE))
              AA(ISTATE)= AA(ISTATE)- CmVALL*Dm/RE(ISTATE)**MMLRP
     2                                *(1.0d0 + Dmp/(PHI(0,ISTATE)*Dm) 
     3                             - MMLRP/(PHI(0,ISTATE)*RE(ISTATE)))
              BB(ISTATE)= BB(ISTATE)- CmVALL*Dm/RE(ISTATE)**MMLRP
     2                                *(2.0d0 + Dmp/(PHI(0,ISTATE)*Dm)
     3                             - MMLRP/(PHI(0,ISTATE)*RE(ISTATE)))
              dAAdRe= dAAdRe - CmVALL*((Dmp*RE(ISTATE) - MMLRP*Dm)
     2            /RE(ISTATE)**(MMLRP+1)*(1.0d0+Dmp/(PHI(0,ISTATE)*Dm)
     3             - MMLRP/(PHI(0,ISTATE)*RE(ISTATE))) 
     4             + Dm/RE(ISTATE)**MMLRP*((Dm*Dmpp-Dmp**2)
     5  /(PHI(0,ISTATE)*Dm**2) + MMLRP/(PHI(0,ISTATE)*RE(ISTATE)**2)))
              dBBdRe= dBBdRe - CmVALL*((Dmp*RE(ISTATE) - MMLRP*Dm)
     1           /RE(ISTATE)**(MMLRP+ 1)*(2.0d0+Dmp/(PHI(0,ISTATE)*Dm)
     2             - MMLRP/(PHI(0,ISTATE)*RE(ISTATE)))
     3             + Dm/RE(ISTATE)**MMLRP*((Dm*Dmpp-Dmp**2)
     4  /(PHI(0,ISTATE)*Dm**2) + MMLRP/(PHI(0,ISTATE)*RE(ISTATE)**2)))
              dVdBtemp= dVdBtemp- CmVALL*Dm/
     1  (RE(ISTATE)**MMLRP*PHI(0,ISTATE)**2)*(MMLRP/RE(ISTATE)-Dmp/Dm)
              ENDDO
          AA(ISTATE)= AA(ISTATE) + DE(ISTATE) 
          BB(ISTATE)= BB(ISTATE) + 2.0d0*DE(ISTATE)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          DO  I= ISTART,ISTOP
c** First, calculate the generalized morse exponent
              RVAL= RDIST
              IF(RDIST.LE.0.d0) THEN
                  RVAL= RMIN(ISTATE) + RH(ISTATE)*DBLE(I-1)
                  ENDIF
              RD(I,ISTATE)= RVAL
              RDp= RVAL**NPB(ISTATE)
              YP= (RDp-REP)/(RDp+REP)
              npow= NSphi(ISTATE)
              IF(RVAL.GT.RE(ISTATE)) npow= NLphi(ISTATE)
              PHIFX(I,ISTATE)= PHI(0,ISTATE)
              DO  J= 1,npow
                  PHIFX(I,ISTATE)= PHIFX(I,ISTATE) +
     1                                         PHI(J,ISTATE)*YP**J
                  ENDDO
c** Calculate some temporary variables.
              XTEMP= DEXP(-PHIFX(I,ISTATE)*(RVAL-RE(ISTATE)))
              Xtemp2= XTEMP*XTEMP
              BMtemp= PHIFX(I,ISTATE)
              BMtemp2= PHIFX(I,ISTATE)
c
c** Now to calculate the actual potential and partial derivatives:
c
              ULR= 0.0d0
              DO  m= 1,NCMM(ISTATE)
                  MMLRP= MMLR(m,ISTATE)
                  CmVALL= CmVAL(m,ISTATE)
                  CALL dampF(RVAL,RHOd(ISTATE),MMLRP,Dm,Dmp,Dmpp,
     1                                                    IDF(ISTATE))
                  ULR= ULR + CmVALL*Dm/RVAL**MMLRP
                  ENDDO
              VPOT(I,ISTATE)= (AA(ISTATE)*XTEMP - BB(ISTATE))*XTEMP 
     1                                            - ULR + VLIM(ISTATE)
              IF(RDIST.GT.0.d0) THEN
                  VDIST= VPOT(I,ISTATE)
                  PHIDIST= PHIFX(I,ISTATE)
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
              IF(RVAL .LT. RE(ISTATE)) THEN
                  DO  J= 1,NSphi(ISTATE)
                      Btemp2= Btemp2 + PHI(J,ISTATE)*DBLE(J)
     1                                                   *YP**(J-1)
                      ENDDO
                  BTEMP= Btemp2
                ELSE
                  DO  J= 1,NLphi(ISTATE)
                      Btemp2= Btemp2 + PHI(J,ISTATE)*DBLE(J) 
     1                                                   *YP**(J-1)
                      ENDDO
                  BTEMP= Btemp2
                ENDIF
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
c ... finally, derivatives of the potential w.r.t. the \phi_i
              IPV= IPV+ 1
              DVtot(IPV,I)= -PBtemp2 + PBTEMP + dVdBtemp*(Xtemp2-XTEMP)
              DO  J= 1,NSphi(ISTATE)
                  IPV= IPV+ 1
                  DVtot(IPV,I)= (-PBtemp2 + PBTEMP) * YP**J
                  ENDDO
              DO  J= NSphi(ISTATE)+1, NLphi(ISTATE)
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
              DO  J= 1, NLphi(ISTATE)
                  DBDB(J,I,ISTATE)= DBDB(J-1,I,ISTATE)*YP
                  ENDDO
              ENDDO
          ENDIF
******7* End Double-Exponential Long-Range Potential Function ********

      IF(PSEL(ISTATE).EQ.4) THEN
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** For the Surkus-type Generalized Potential Energy Function.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** First, we calculate the implied Dissociation Energy
        IF(AGPEF(ISTATE).NE.0.0d0) THEN
            YPP= 1.d0/AGPEF(ISTATE)**2
            VAL= YPP
            DO  I= 1, NLphi(ISTATE)
                YPP= YPP/AGPEF(ISTATE)
                VAL= VAL + PHI(I,ISTATE)*YPP
                ENDDO
            DE(ISTATE)= VAL*PHI(0,ISTATE)
            ENDIF
        DO  I= 1, NDATPT(ISTATE) 
            RVAL= RDIST
            IF(RDIST.LE.0.d0) THEN
                RVAL= RMIN(ISTATE) + RH(ISTATE)*DBLE(I-1)
                ENDIF
            RD(I,ISTATE)= RVAL
            RDp= RVAL**nPB(ISTATE)
            YP= (RDp-REP)/(AGPEF(ISTATE)*RDp + BGPEF(ISTATE)*REP)
c** Now to calculate the actual potential
            YPP= 1.d0
            VAL= 1.d0
            DVAL= 2.d0
            DO  J= 1, NLphi(ISTATE)
                YPP= YPP*YP
                VAL= VAL + PHI(J,ISTATE)*YPP
                DVAL= DVAL+ (J+2)*PHI(J,ISTATE)*YPP
                ENDDO
            VPOT(I,ISTATE)= VAL*PHI(0,ISTATE)*YP**2 + VLIM(ISTATE) 
            IF(RDIST.GT.0) THEN
                VDIST= VPOT(I,ISTATE)
                PHIDIST= 0.d0
                IF(IDAT.LE.0) GO TO 999
                ENDIF
            DVAL= DVAL*PHI(0,ISTATE)*YP
c** Now to calculate the partial derivatives
            DVDD= 0.d0
            IPV= IPVSTART + 1
c ... derivative of the potential w.r.t. Re
            DVtot(IPV,I)= -DVAL*REP*RDp*(AGPEF(ISTATE)+BGPEF(ISTATE))
     1                  *(nPB(ISTATE)/RE(ISTATE))/(AGPEF(ISTATE)*RDp +
     2                                           BGPEF(ISTATE)*REP)**2
c ... and derivatives w.r.t. the phi_i expansion coefficients ...
            IPV= IPV+ 1
            DVtot(IPV,I)= VAL*YP**2
            IPV= IPV+ 1
            DVtot(IPV,I)= PHI(0,ISTATE)*YP**3
            DO  J= 2, NLphi(ISTATE)
                IPV= IPV+ 1
                DVtot(IPV,I)= DVtot(IPV-1,I)*YP
                ENDDO
            ENDDO
        IF(RDIST.LE.0.d0) VLIM(ISTATE)= VPOT(NDATPT(ISTATE),ISTATE)
c????
cc      rewind(10)
cc      write(10,612) (RD(i,ISTATE),vpot(i,istate),
cc   1                                        i= 1, NDATPT(ISTATE),50)
cc612 FORMAT(/(f10.4,f15.5))
c????
        ENDIF
c*********** End Generalized Potential Energy Function *****************

      IF((NUA(ISTATE).GE.0).OR.(NUB(ISTATE).GT.0)) THEN
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Treat any 'adiabatic' BOB radial potential functions here ...
c        u_A(r) = ym*uA_\infty + [1 - ym]\sum_{i=0,NUA} {uA_i yp^i} 
c   where the  u_\infty  values stored/fitted as  UA(NUA(ISTATE))
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          REP= RE(ISTATE)**nPu(ISTATE)
          REm= REP
          HRe= 0.5d0*nPu(ISTATE)/RE(ISTATE)
          HRem= HRe
          IF(PSEL(ISTATE).EQ.2) THEN
              REm= RE(ISTATE)**MMLR(1,ISTATE)
              HRem= 0.5d0*MMLR(1,ISTATE)/RE(ISTATE)
              ENDIF
          IF((BOBCN(ISTATE).GE.1).AND.(nPu(ISTATE).EQ.0)) THEN
              HRe= 2.d0*HRe
              HRem= 2.d0*HRem
              ENDIF
c ... reset parameter counter ...
          IPVSTART= IPV
          DO  I= ISTART,ISTOP
              RVAL= RD(I,ISTATE)
              RDp= RVAL**nPu(ISTATE)
              YP= (RDp - REP)/(RDp + REP)
              YPM= 1.d0 - YP
              YM= YP
              IF(PSEL(ISTATE).EQ.2) THEN
                  RDm= RVAL**MMLR(1,ISTATE)
                  YM= (RDm - REm)/(RDm + REm)
                  ENDIF
              IF(BOBCN(ISTATE).GE.1) THEN
                  YM= 0.d0
c** If  BOBCN > 0  &  p= 1,  assume use of Ogilvie-Tipping vble.
                  IF(nPu(ISTATE).EQ.1) YP= 2.d0*YP
                  ENDIF
              YMM= 1.d0 - YM
              IF(NUA(ISTATE).GE.0) THEN
c ... Now ... derivatives of UA w,r,t, expansion coefficients
                  VAL= UA(0,ISTATE)
                  DVAL= 0.d0
                  IPV= IPVSTART + 1
                  DVtot(IPV,I)= YMM
                  YPP= 1.d0
                  IF(NUA(ISTATE).GE.2) THEN
                      DO  J= 1,NUA(ISTATE)-1
                          DVAL= DVAL+ DBLE(J)*YPP*UA(J,ISTATE)
                          YPP= YPP*YP
                          VAL= VAL+ UA(J,ISTATE)*YPP
                          IPV= IPV+ 1
                          DVtot(IPV,I)= YMM*YPP
                          ENDDO
                      ENDIF
                  IPV= IPV + 1
                  DVtot(IPV,I)= YM
                  UAR(I,ISTATE)= VAL*YMM + YM*UA(NUA(ISTATE),ISTATE)
c ... and derivative of UA w.r.t. Re ...
                  DUADRe(I,ISTATE)= -HRe*(1.d0 - YP**2)*YMM*DVAL
                  IF(BOBCN(ISTATE).LE.0) DUADRe(I,ISTATE)= 
     1                          DUADRe(I,ISTATE) + HRem*(1.d0 - YM**2)
     2                                 *(VAL - UA(NUA(ISTATE),ISTATE))
                  ENDIF
              IF(NUB(ISTATE).GE.0) THEN
c ... Now ... derivatives of UB w,r,t, expansion coefficients
                  VAL= UB(0,ISTATE)
                  DVAL= 0.d0
                  IF(NUA(ISTATE).LT.0) THEN
                      IPV= IPVSTART + 1
                    ELSE
                      IPV= IPV + 1
                    ENDIF
                  DVtot(IPV,I)= YMM
                  YPP= 1.d0
                  IF(NUB(ISTATE).GE.2) THEN
                      DO  J= 1,NUB(ISTATE)-1
                          DVAL= DVAL+ DBLE(J)*YPP*UB(J,ISTATE)
                          YPP= YPP*YP
                          VAL= VAL+ UB(J,ISTATE)*YPP
                          IPV= IPV + 1
                          DVtot(IPV,I)= YMM*YPP
                          ENDDO
                      ENDIF
                  IPV= IPV + 1
                  DVtot(IPV,I)= YM
                  UBR(I,ISTATE)= VAL*YMM + YM*UB(NUB(ISTATE),ISTATE)
c ... and derivative of UB w.r.t. Re ...
                  DUBDRe(I,ISTATE)= -HRe*(1.d0 - YP**2)*YMM*DVAL 
                  IF(BOBCN(ISTATE).LE.0) DUBDRe(I,ISTATE)= 
     1                          DUBDRe(I,ISTATE) + HRem*(1.d0 - YM**2)
     2                                 *(VAL - UB(NUB(ISTATE),ISTATE))
                  ENDIF
              ENDDO
          ENDIF
c++++ END of treatment of adiabatic potential BOB function++++++++++++++

      IF((NTA(ISTATE).GE.0).OR.(NTB(ISTATE).GE.0)) THEN
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Treat any 'non-adiabatic' centrifugal BOB functions here ...
c       q_A(r) = yp*qA_\infty + [1 - yp]\sum_{i=0,NTA} {qA_i yp^i}
c   where the  q_\infty  values stored/fitted as  TA(NTA(ISTATE))
c  Incorporate the 1/r^2 factor into the partial derivatives (but not in
c     the g(r) functions themselves, since pre-SCHRQ takes care of that).
c    Need to add  M_A^{(1)}/M_A^{(\alpha)}  factor later too
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          REP= RE(ISTATE)**nPt(ISTATE)
          HRe= 0.5d0*nPt(ISTATE)/RE(ISTATE)
          IF((BOBCN(ISTATE).GE.1).AND.(nPt(ISTATE).EQ.0)) HRe= 2.d0*HRe
          IPVSTART= IPV
          DO  I= ISTART,ISTOP
              RVAL= RD(I,ISTATE)
              RM2= 1/RVAL**2
              RDp= RVAL**nPt(ISTATE)
              YP= (RDp - REP)/(RDp + REP)
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
                  YPP= 1.d0
                  IF(NTA(ISTATE).GE.2) THEN
                      DO  J= 1,NTA(ISTATE)-1
                          DVAL= DVAL+ DBLE(J)*YPP*TA(J,ISTATE)
                          YPP= YPP*YP
                          VAL= VAL+ TA(J,ISTATE)*YPP
                          IPV= IPV + 1
                          DVtot(IPV,I)= YPM*YPP*RM2
                          ENDDO
                      ENDIF
                  IPV= IPV + 1
                  DVtot(IPV,I)= YP*RM2
                  TAR(I,ISTATE)= VAL*YPM 
                  YPP= -HRe*(1.d0 - YP**2)*RM2
                  DTADRe(I,ISTATE)= YPP*YPM*DVAL
                  IF(BOBCN(ISTATE).LE.0) THEN
                      DTADRe(I,ISTATE)= DTADRe(I,ISTATE)
     1                 + YPP*(YPM*DVAL - VAL + TA(NTA(ISTATE),ISTATE))
                      TAR(I,ISTATE)= TAR(I,ISTATE) + 
     1                                       YP*TA(NTA(ISTATE),ISTATE)
                      ENDIF
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
                  YPP= 1.d0
                  IF(NTB(ISTATE).GE.2) THEN
                      DO  J= 1,NTB(ISTATE)-1
                          DVAL= DVAL+ DBLE(J)*YPP*TB(J,ISTATE)
                          YPP= YPP*YP 
                          VAL= VAL+ TB(J,ISTATE)*YPP
                          IPV= IPV + 1
                          DVtot(IPV,I)= YPM*YPP*RM2
                          ENDDO
                      ENDIF
                  IPV= IPV + 1
                  DVtot(IPV,I)= YP*RM2
                  TBR(I,ISTATE)= VAL*YPM 
                  YPP= -HRe*(1.d0 - YP**2)*RM2
                  DTBDRe(I,ISTATE)= YPP*YPM*DVAL
                  IF(BOBCN(ISTATE).LE.0) THEN
                      DTBDRe(I,ISTATE)= DTBDRe(I,ISTATE)
     1                 + YPP*(YPM*DVAL - VAL + TB(NTB(ISTATE),ISTATE))
                      TBR(I,ISTATE)= TBR(I,ISTATE) + 
     1                                       YP*TB(NTB(ISTATE),ISTATE)
                      ENDIF
                  ENDIF
              ENDDO
          ENDIF
c++++ END of treatment of non-adiabatic centrifugal BOB function++++++++

      IF(NwCFT(ISTATE).GE.0) THEN
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Treat any Lambda- or 2\Sigma-doubling radial strength functions here
c    representing it as   f(r)= yp*w_inf + [1-yp]*Sum{ w_i * yp^i}
c    where the  w_inf  values stored/fitted as  wCFT(NwCFT(ISTATE))
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          LAMB2= 2*IOMEG(ISTATE)
          REP= RE(ISTATE)**nPw(ISTATE)
          HRe= 0.5d0*nPw(ISTATE)/RE(ISTATE)
          IPVSTART= IPV
          DO  I= ISTART,ISTOP
              RVAL= RD(I,ISTATE)
              RMF= 1.d0/RVAL**2
              IF(IOMEG(ISTATE).GT.0)  RMF= RMF**LAMB2
              RDp= RVAL**nPw(ISTATE)
              YP= (RDp - REP)/(RDp + REP)
              YPM= 1.d0 - YP
              DVAL= 0.d0
              YPP= YPM*RMF
              VAL= wCFT(0,ISTATE)*YPP
              IPV= IPVSTART + 1
              DVtot(IPV,I)= YPP
              IF(NwCFT(ISTATE).GE.2) THEN
                  DO  J= 1,NwCFT(ISTATE)-1
                      DVAL= DVAL+ DBLE(J)*YPP*wCFT(J,ISTATE)
                      YPP= YPP*YP
                      IPV= IPV + 1
                      DVtot(IPV,I)= YPP
                      VAL= VAL+ wCFT(J,ISTATE)*YPP
                      ENDDO
                  ENDIF
              IPV= IPV + 1
              DVtot(IPV,I)= YP*RMF
              wRAD(I,ISTATE)= VAL + YP*wCFT(NwCFT(ISTATE),ISTATE)*RMF
              dLDDRe(I,NSTATEMX)= -HRe*(1.d0 - YP**2)*(DVAL - VAL +
     1                                 wCFT(NwCFT(ISTATE),ISTATE)*RMF)
              ENDDO
          ENDIF
c++++ END of treatment of Lambda/2-sigma centrifugal BOB function+++++++

c++++ Test for inner wall inflection, and if it appears,  ++++++++++++++
c++++ replace inward potential with linear approximation +++++++++++++++
      I1= (RE(ISTATE)-RMIN(ISTATE))/RH(ISTATE)
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
   66 CONTINUE
  606 FORMAT(9('===')/'!*!* Extrapolate to correct ',A2,' inner-wall inf
     1lection at   R=',f6.4,'   V=',f8.0/9('==='))
c++++++++++++End of Inner Wall Test/Correction code+++++++++++++++++++++
c======================================================================
c** At the one distance RDIST calculate total effective potential VDIST
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
          BFCT= 16.85762920d0/(ZMASS(3,IISTP)*RDIST**2)
          JFCT= DBLE(JPP(IDAT)*(JPP(IDAT)+1))
          IF(IOMEG(ISTATE).GT.0) JFCT= JFCT - IOMEG(ISTATE)**2
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
          IF(IOMEG(ISTATE).NE.0) THEN
              IF(IOMEG(ISTATE).GT.0) THEN
c ... for Lambda doubling case ...
                  JFCTLD= (EFPP(IDAT)-EFREF(ISTATE)) 
     1         *(DBLE(JPP(IDAT)*(JPP(IDAT)+1))*BFCT**2)**IOMEG(ISTATE)
                  ENDIF
              IF(IOMEG(ISTATE).EQ.-2) THEN
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

c***********************************************************************
      SUBROUTINE dampF(R,RHOd,m,Dm,Dmp,Dmpp,IDF)
c** Subroutine to generate values 'Dm' and its first `Dmp' and second
c  'Dmpp' derivatives w.r.t. R of the chosen potential damping function
c   IDF = 1 :  Tang-Toennies damping function [JCP 80,3726(1984)]
c   IDF = 2 :  damped-dispersion damping function [MolP 52,763(1984)]
c   IDF = 3 : extend IDF=2 by multiplying by Scoles' global damping fx.
c  all assuming  R  in units Angst. 
c***********************************************************************
      INTEGER i,m,IDF
      REAL*8 R,RHOd,RHOdR,Dm,Dmp,Dmpp,fd,fdp,fdpp,Coef,XPOW
      REAL*8 sum(m)
c
      RHOdR= RHOd*R
      IF(IDF.EQ.1) THEN
c** Tang-Toennies damping function: numerical factor 3.16 makes these
c  functions equivalent to Scoles functions, and allows the  RHOd  
c  scaling factor to have the same significance in both cases.
          RHOdR= 3.16d0*RHOdR
          XPOW = DEXP(-RHOdR)
          DO  i= 1,m
              sum(i) = 1.0d0
              END DO
          fd = 1.0d0
          DO  i= 1,m
              fd = fd*RHOdR/i
              sum(i) = sum(i) + fd
              ENDDO
          Dm = 1.0d0-XPOW*sum(m)
          Dmp = 3.16d0*RHOd*XPOW*(sum(m) - sum(m-1))
          Dmpp = 3.16d0*RHOd**2*XPOW*(-sum(m) + 2*sum(m-1) - sum(m-2))
          ENDIF
c
      IF((IDF.EQ.2).OR.(IDF.EQ.3)) THEN
c** Douketis et al. [MolP 52, 763 (1984)]  damping function.
          XPOW= DEXP(-3.97d0*RHOdR/m-0.39d0*(RHOdR)**2/DSQRT(m*1.0d0))
ccc       XPOW= DEXP(-3.968424881d0*RHOdR/m
ccc  1                      - 0.389246070d0*(RHOdR)**2/DSQRT(m*1.0d0))
          Dm = (1.0d0 - XPOW)**m
          Coef = 3.97d0*RHOd/m + 0.78d0*RHOd*RHOdR/DSQRT(m*1.0d0)
ccc       Coef = 3.968424881d0*RHOd/m 
ccc  1                      + 0.7784921396d0*RHOd*RHOdR/DSQRT(m*1.0d0)
          Dmp = m*Coef*XPOW*(1.0d0 - XPOW)**(m-1)
          Dmpp = XPOW*(1.0d0 - XPOW)**(m-1)
     1                                  *(0.78d0*SQRT(m*1.0d0)*RHOd**2
ccc  1                          *(0.7784921396d0*SQRT(m*1.0d0)*RHOd**2 
     2              - m*Coef**2*(1.0d0-m*XPOW)/(1.0d0 - XPOW))
          ENDIF
c
      IF(IDF.EQ.3) THEN
c** If desired, also apply Scoles' global scaling/damping factor
          XPOW= DEXP(-1.474d0*RHOdR)
          Coef= (RHOdR)**1.68d0
          fd= Coef*XPOW
          fdp= RHOd*fd*(1.474d0 - 1.68/(RHOdR))
          fdpp= RHOd**2 *fd*(-2.172676d0 + 
     1                       (4.95264d0 - 1.1424d0/(RHOdR))/(RHOdR))
          fd= 1.d0 - fd
          Dmpp= Dmpp*fd + 2.d0*Dmp*fdp + Dm*fdpp
          Dmp= Dmp*fd + Dm*fdp
          Dm= Dm*fd
          END IF
      RETURN
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE MKPREDICT(NSTATES,NDAT)
c***********************************************************************
c** Subroutine to prepare fake input data array which will cause ParFit
c  to make transition energy predictions for electronic or infrared band
c  or microwave transitions.  On entry:
c  NSTATES  is the number of states involved in the data set.  
c    NSTATES= 1 generates infrared or microwave bands for state SLABL(1)
c    NSTATES= 2 generates electronic bands from lower state SLABL(1) 
c               into upper state SLABL(2)
c  VMIN(s) and VMAX(s) are the bounds on the vibrational energy range 
c      for state 's' specified in the main input file.
c** On return:
c  NDAT(v,i,s)  is the number of transitions associated with
c    vibrational level-v of isotopomer-i of state-s [for NDEGB < 0 case]
c** This subroutine reads in band specifications on Channel-5 and writes
c   the transition energy specifications to channel-4
c-----------------------------------------------------------------------
c                         Version of 1 September 2005
c-----------------------------------------------------------------------
      INCLUDE 'arrsizes.h'
      INCLUDE 'BLKISOT.h'
      INCLUDE 'BLKDATA.h'
      INCLUDE 'BLKTYPE.h'
c-----------------------------------------------------------------------
c
      CHARACTER*2 LABLP,LABLPP
      INTEGER I,J,J2,JD,J2DL,J2DU,J2DD,JMAXX,PP,PPP,NTRANS,COUNT,
     1  IBAND,JMAX(NPARMX),JMIN(NPARMX),
     1  VMX(NSTATEMX),ISOT,ESP,ESPP,ISTATE,MN1,MN2
      INTEGER NSTATES,NDAT(0:NVIBMX,NISTPMX,NSTATEMX)
c-----------------------------------------------------------------------
c** Initialize counters for book-keeping on input data
      COUNT= 0
      DO  ISOT= 1,NISTP
          DO  ISTATE= 1,NSTATES
              NTRANSFS(ISOT,ISTATE)= 0
              NTRANSIR(ISOT,ISTATE)= 0
              NTRANSMW(ISOT,ISTATE)= 0
              NBANDFS(ISOT,ISTATE)= 0
              NBANDVIS(ISOT,ISTATE)= 0
              NBANDIR(ISOT,ISTATE)= 0
              NBANDMW(ISOT,ISTATE)= 0
              NBVPP(ISOT,ISTATE)= 0
              NWIDTH(ISOT,ISTATE)= 0
              DO  I= 1,NSTATES
                  NTRANSVIS(ISOT,ISTATE,I)= 0
                  NBANDEL(ISOT,ISTATE,I)= 0
                  ENDDO
              ENDDO
          NBANDS(ISOT)= 0
          ENDDO
      DO  ISTATE= 1,NSTATES
          VMX(ISTATE)= 0
          ENDDO
      NFSTOT= 0
      IBAND= 0
   70 IBAND= IBAND+ 1
      IF(IBAND.GT.NPARMX) THEN
            WRITE(6,609) IBAND,NPARMX
            IBAND= IBAND-1
            GOTO 99
            ENDIF
c** Generate "empty" band data sets to allow ParFit to make predictions
c  for those sets of transitions.  
c** LABLP & LABLPP are the two-character variables identifying the upper
c     and lower electronic states, respectively.  LABLP=LABLPP for IR or
c     MW transitions within a given electronic state
c** VP & VPP are the v' & v" values identifying the band;
c** PP & PPP specify rotational parities (+/- 1) of upper and lower levels
c** MN1 & MN2 identify the isotopomer
c** Generate 'lines' for  J"= 0 to JMAXX subject to selection rule that
c  Delta(J) runs from J2DL to J2DU in steps of J2DD
c-----------------------------------------------------------------------
      READ(5,*,end=99) VP(IBAND),VPP(IBAND),LABLP,LABLPP,MN1,MN2,PP,PPP,
     1                                            JMAXX,J2DL,J2DU,J2DD
c-----------------------------------------------------------------------
      IF(VP(IBAND).LT.0) GO TO 99
c** Set electronic state number for upper & lower levels.  
c* Always set lower state as 1'st state considered in input [SLABL(1)]
c* For NSTATES= 1, upper state is the same one.  For NSTATES= 2 the 
c  upper state is 2'nd one considered [SLABL(2)]
      IEPP(IBAND)= 1
      IEP(IBAND)= NSTATES
      WRITE(4,400) VP(IBAND),VPP(IBAND),LABLP,LABLPP,MN1,MN2
      ISOT= 0
c** Determine the correct isotopomer-number for this band.
      DO  I= 1,NISTP
          IF((MN1.EQ.MN(1,I)).AND.(MN2.EQ.MN(2,I))) ISOT= I
          ENDDO
      ISTP(IBAND)= ISOT
      MAXUFREQ(IBAND)= 0
      JMAX(IBAND)= JMAXX
      JMIN(IBAND)= 0
      NTRANS= 0
      IFIRST(IBAND)= COUNT+ 1
      ESP= IEP(IBAND)
      ESPP= IEPP(IBAND)
c** Now - loop over J to generate all possible transitions ...
      DO  J= 0, JMAXX
          DO  JD= J2DL, J2DU, J2DD
              J2= J+ JD
              IF((J2.GE.0).AND.((J.NE.0).OR.(J2.NE.0))) THEN
                  COUNT= COUNT+1
                  IF(COUNT.GT.NDATAMX) THEN
                      WRITE(6,640) COUNT,NDATAMX
                      STOP
                      ENDIF
                  WRITE(4,402) J2,PP,J,PPP
                  JP(COUNT)= J2
                  EFP(COUNT)= PP
                  JPP(COUNT)= J
                  EFPP(COUNT)= PPP
                  FREQ(COUNT)= 0.d0
                  UFREQ(COUNT)= 0.001d0
                  DFREQ(COUNT)= 0.d0
                  IB(COUNT)= IBAND
c** Accumulate count of data associated with each vibrational level ...
                  NDAT(VPP(IBAND),ISTP(IBAND),ESPP)=
     1                            NDAT(VPP(IBAND),ISTP(IBAND),ESPP)+ 1
                  NDAT(VP(IBAND),ISTP(IBAND),ESP)=
     1                              NDAT(VP(IBAND),ISTP(IBAND),ESP)+ 1
                  ENDIF
              ENDDO
          ENDDO
      WRITE(4,404)
  400 FORMAT(2I4,"   '",A2,"'  '",A2,"'   ",2I4)
  402 FORMAT(I4,I3,I5,I3,'    0.d0     1.0d-3')
  404 FORMAT('  -1 -1   -1 -1    -1.d0    -1.d-3'/)
      VMX(ESP)= MAX(VMX(ESP),VP(IBAND))
      VMX(ESPP)= MAX(VMX(ESPP),VPP(IBAND))
      ILAST(IBAND)= COUNT
      NTRANS= ILAST(IBAND)-IFIRST(IBAND)+1
      GOTO 70
   99 RETURN
  609 FORMAT(/' *** ERROR *** Dimension allocated for number of bands ex
     1ceeded:'/' (IBAND=',i4,') > (NBANDMX=',i4,')   so truncate input a
     2nd TRY to continue ...')
  640 FORMAT(/' *** Input Data Count reaches',i6,' which EXCEEDS ARRAY L
     1IMIT of',i6)
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE DYIDPJ(IDAT,NDATA,NPTOT,JFXP,YOBS,YC,PV,PD,PSS,RMSR)
c***********************************************************************
c** This program assumes that the upper state IEP is at a higher energy
c   than the lower state IEPP.
c** This subroutine returns the calculated value YC of datum IDAT, and
c  its partial derivatives PD(k) w.r.t. the NPTOT parameters PV.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++  COPYRIGHT 2007  by  Robert J. Le Roy, Jenning Seto and Yiye Huang +
c   Dept. of Chemistry, Univ. of Waterloo, Waterloo, Ontario, Canada   +
c    This software may not be sold or any other commercial use made    +
c      of it without the express written permission of the authors.    +
c++++++++++++++++++++ (version of 12/01/2007) ++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** On entry:
c    IDAT     is the number of the current observable being considered. 
c    NPTOT    is the total number of parameters in the model
c    PV(i)    is the array of parameters being varied.
c
c** On exit:
c    YC       is the calculated value of the IDATth observable.
c    PV(i)    is the array of parameters being varied.
c    PD(i)    is the partial derivative array (de/dp).
c=======================================================================
      INCLUDE 'arrsizes.h'
      INCLUDE 'BLKISOT.h'
      INCLUDE 'BLKDATA.h'
      INCLUDE 'BLKPOT.h'
      INCLUDE 'BLKCOUNT.h'
c-----------------------------------------------------------------------
      INTEGER IDAT,NPTOT,NBAND,IISTP,ISTATE,VIMX,I,J,JFXP(NPARMX),NDATA
c
c** Define parameters required locally and from NLLSSRR.
      REAL*8 VDIST,PHIDIST,EUP,ELW,YOBS,YC,EO,RMSR,width,
     1  PV(NPARMX),PD(NPARMX),PSS(NPARMX),UPPER(NPARMX),LOWER(NPARMX),
     2  DEDPK(HPARMX)
c
c----------------------------------------------------------------------- c
c** Common block carrying term values and term-value count labels
      REAL*8 TVALUE(NPARMX)
      INTEGER NSTATES,NTVALL(0:NSTATEMX),NTVI(NSTATEMX),NTVF(NSTATEMX),
     1  VMIN(NSTATEMX),VMAX(NSTATEMX)
      COMMON /NLSBLK/TVALUE,NSTATES,NTVALL,NTVI,NTVF,VMIN,VMAX
c 
c** Vibrational energies and rotational constants.  
c
      REAL*8 ZK(0:NVIBMX,0:NROTMX,NISTPMX,NSTATEMX)
      COMMON /ZKBLK/ZK
c
      IF(IDAT.EQ.1) THEN
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++ At beginning of each fit cycle (datum #1), re-map internal NLLSSRR
c  parameter array PV onto external (physical) variable set, and get
c  updated band constant array ZK for estimating trial eigenvalues.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          CALL MAPPAR(NSTATES,PV,1)
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
c** Call subroutine to update potential functions and their partial 
c   derivatives w.r.t. fitting parameter for each state.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
              IF(PSEL(ISTATE).GT.0) CALL VGEN(ISTATE,-1.d0,VDIST,
     1                                                      PHIDIST,0)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Now generate band constants for each state and isotopomer for 
c  generating trial eigenvalues in fit calculations.  To take account of
c  the different vibrational ranges for different isotopomers, only 
c  generate values for levels to the input VMAX for isotopomer-1.
              IF(PSEL(ISTATE).GE.0) THEN
                  DO  IISTP= 1,NISTP
                      VIMX=VMAX(ISTATE)*DSQRT(ZMASS(3,IISTP)/ZMASS(3,1))
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Call subroutine INITDD to calculate band constants for initial
c   trial eigenvalue estimates as well as the partial derivatives.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                      CALL INITDD(ISTATE,IISTP,VIMX,ZK)
                      ENDDO
                  ENDIF
              ENDDO
          ENDIF
c ++++++++++++++++++++++++++++++++++++++++++end of datum-1 set-up ++++++
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
c   state for fluorescence data.
c
      IF(IEP(NBAND).EQ.-1) THEN
c** First ... for PAS data ...
          CALL DEDP(IDAT,IEPP(NBAND),IISTP,ZMASS(3,IISTP),
     1               JP(IDAT),JPP(IDAT),EFPP(IDAT),ELW,width,ZK,LOWER)
              IF(ELW.LT.-9.d9) THEN
c... if eigenvalue search failed, remove this datumn from the fit
                  WRITE(6,600) SLABL(IEPP(NBAND)),JP(IDAT),JPP(IDAT),
     1                                                        IDAT,YOBS
                  YC= YOBS
                  RETURN
                  ENDIF
          EUP= VLIM(IEPP(NBAND))
        ELSEIF(IEP(NBAND).EQ.0) THEN
c** For fluorescence series data ...
          I= NFS1-1
          DO J= NFS1,NPTOT
              IF(FSBAND(J-I).EQ.NBAND) THEN
c ... PV(HPARMX+I) is the energy of the I'th fluorescence band.
                  EUP= TVALUE(J)
                  UPPER(J)= 1.0d0
                  ENDIF
              ENDDO
          CALL DEDP(IDAT,IEPP(NBAND),IISTP,ZMASS(3,IISTP),
     1               JP(IDAT),JPP(IDAT),EFPP(IDAT),ELW,width,ZK,LOWER)
              IF(ELW.LT.-9.d9) THEN
c... if eigenvalue search failed, remove this datumn from the fit
                  WRITE(6,600) SLABL(IEPP(NBAND)),JP(IDAT),JPP(IDAT),
     1                                                       IDAT,YOBS
                  YC= YOBS
                  RETURN
                  ENDIF
        ELSEIF(IEP(NBAND).EQ.-2) THEN
c** For widths of tunneling predissociation quasibound levels
c ... for forward calculation, use widths from SCHRQ in DEDP
          CALL DEDP(IDAT,IEPP(NBAND),IISTP,ZMASS(3,IISTP),
     1                JP(IDAT),JPP(IDAT),EFPP(IDAT),EO,width,ZK,DEDPK)
              IF(EO.LT.-9.d9) THEN
c... if eigenvalue search failed, remove this datumn from the fit
                  WRITE(6,600) SLABL(IEPP(NBAND)),JP(IDAT),JPP(IDAT),
     1                                                       IDAT,YOBS
                  YC= YOBS
                  RETURN
                  ENDIF
c ... otherwise ... calculate 'width' and its detivatives in DWDP
          IF(PSEL(IEPP(NBAND)).GT.0) THEN
              CALL DWDP(IDAT,IEPP(NBAND),FREQ(IDAT),ZMASS(3,IISTP),
     1                           JP(IDAT),JPP(IDAT),EO,width,DEDPK,PD)
              YC= width 
              RETURN
              ENDIF
c 
        ELSEIF(IEP(NBAND).GT.0) THEN
c** Alternately ... determine the partials for the upper and lower
c   states for microwave, infrared, and electronic data.
          IF(PSEL(IEP(NBAND)).EQ.-2) THEN
c... if upper state being represented by term values ...
              UPPER(TVUP(IDAT))= 1.d0
              EUP= TVALUE(TVUP(IDAT))
            ELSE
c... for normal case of lower state being represented by a potential
              CALL DEDP(IDAT,IEP(NBAND),IISTP,ZMASS(3,IISTP),
     1                VP(NBAND),JP(IDAT),EFP(IDAT),EUP,width,ZK,UPPER)
              IF(EUP.LT.-9.d9) THEN
c... if eigenvalue search failed, remove this datumn from the fit
                  WRITE(6,600) SLABL(IEP(NBAND)),VP(NBAND),JP(IDAT),
     1                                                        IDAT,YOBS
                  YC= YOBS
                  RETURN
                  ENDIF
            ENDIF
          IF(PSEL(IEPP(NBAND)).EQ.-2) THEN
c... if lower state being represented by term values ...
              LOWER(TVLW(IDAT))= -1.d0
              ELW= TVALUE(TVLW(IDAT))
            ELSE
c... for normal case of lower state being represented by a potential
              CALL DEDP(IDAT,IEPP(NBAND),IISTP,ZMASS(3,IISTP),
     1             VPP(NBAND),JPP(IDAT),EFPP(IDAT),ELW,width,ZK,LOWER)
              IF(ELW.LT.-9.d9) THEN
c... if eigenvalue search failed, remove this datumn from the fit
                  WRITE(6,600) SLABL(IEPP(NBAND)),VPP(NBAND),JPP(IDAT),
     1                                                        IDAT,YOBS
                  YC= YOBS
                  RETURN
                  ENDIF
            ENDIF
        ENDIF
      DO  I= 1,NPTOT
          PD(I)= UPPER(I) - LOWER(I)
          ENDDO
c----------------------------------------------------------------------
c** Get calculated value for the IDATth observable from energy levels
c----------------------------------------------------------------------
      YC= EUP - ELW
c----------------------------------------------------------------------
      RETURN
  600 FORMAT(' *** FAIL to find level(',A2,')   v=',I3,'   J=',I3,
     1  '  so ignore  YOBS(',i4,')=',f12.4)
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE INITDD(ISTATE,IISTP,VIMX,ZK)
c***********************************************************************
c** This subroutine updates the trial vibrational energies & rotational
c        constants on each iteration
c** On entry:   ISTATE  is the states being considered.
c               IISTP  is the isotopomer being considered.
c               VIMX   is the upper vibrational bound for each state.
c** On exit:   ZK       are the band constants for this state & isotope
c** Internal:  RM2   is the  (1+ ZMTA*TAR + ZMTB*TBR)/r**2  array for 
c                    this state required by CDJOEL for CDC calculation
c=======================================================================
      INCLUDE 'arrsizes.h'
      INCLUDE 'BLKPOT.h'
      INCLUDE 'BLKISOT.h'
      INCLUDE 'BLKBOB.h'
      INCLUDE 'BLKBOBRF.h'
c-----------------------------------------------------------------------
c
      INTEGER ISTATE, IISTP, VIMX, AFLAG, INNR(0:NVIBMX),
     1  VIN, NBEG, NEND, WARN, IWR, LPRWF, I, J, INNODE
c
      REAL*8 ZK(0:NVIBMX,0:NROTMX,NISTPMX,NSTATEMX),
     1  GV(0:NVIBMX),RCNST(NROTMX),
     2  RM2(NPNTMX),FWHM,PMAX,BFCT,BvWN,V(NPNTMX),SWF(NPNTMX)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IF((ISTATE.EQ.1).AND.(IISTP.EQ.1)) THEN
          REWIND(7)
          REWIND(21)
          ENDIF
      VIN= VIMX
      AFLAG= 0
      WARN= 0
      IWR= -1
      LPRWF= 0
      INNODE= 1
c
c** Calculate potential scaling factors for ALF and CDJOEL
c
      BFCT= (ZMASS(3,IISTP)/16.85762920D0)*RH(ISTATE)*RH(ISTATE)
      BvWN= 16.85762920D0/ZMASS(3,IISTP)
c
c** Now generate the BFCT-scaled and adiabatically corrected potential
c   for the current isotope for use in SCHRQ, plus the 1/R**2 array for
c   CDC calculations in CDJOEL
c
      DO  I=1,NDATPT(ISTATE)
          V(I)= BFCT*(VPOT(I,ISTATE)+ ZMUA(IISTP,ISTATE)*UAR(I,ISTATE)
     1                             + ZMUB(IISTP,ISTATE)*UBR(I,ISTATE))
          RM2(I)= 1.d0/RD(I,ISTATE)**2
          ENDDO
      IF((NTA(ISTATE).GE.0).OR.(NTB(ISTATE).GE.0)) THEN
          DO  I= 1, NDATPT(ISTATE)
              RM2(I)= RM2(I)*(1.d0 + ZMTA(IISTP,ISTATE)*TAR(I,ISTATE) 
     1                             + ZMTB(IISTP,ISTATE)*TBR(I,ISTATE))
              ENDDO
          ENDIF
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Call subroutine ALF that will locate the needed vibrational levels 
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      CALL ALF(NDATPT(ISTATE),RMIN(ISTATE),RH(ISTATE),V,SWF,
     1 VLIM(ISTATE),VIN,AFLAG,ZMASS(3,IISTP),EPS(ISTATE),MMLR(1,ISTATE),
     2                                            GV,INNODE,INNR,IWR)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** If a serious error occured during within ALF, then print out a 
c   warning, record the constants that we have and hope the program
c   doesn't call on the constants for which ALF could not calculate.
c
      IF (AFLAG.LT.0) THEN
          WRITE(6,600) ISTATE,IISTP
          IF (AFLAG.EQ.-1) WRITE(6,601) VIMX, AFLAG
          IF (AFLAG.EQ.-2) WRITE(6,602)
          IF (AFLAG.EQ.-3) WRITE(6,603)
          IF (AFLAG.EQ.-4) WRITE(6,604)
          IF (AFLAG.EQ.-5) WRITE(6,606)
          IF (AFLAG.EQ.-6) WRITE(6,608)
          IF (AFLAG.EQ.-8) THEN
              WRITE(6,610)
        ELSE
          WRITE(6,612) VIN
        ENDIF
      ENDIF
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Now to calculate the rotational constants for each vibrational level
c   of each state and isotopomer.
      WRITE(7,614)
      DO  I= 0,VIN
          AFLAG= 0
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Call subroutine SCHRQ to calculate the wavefunction required by
c   CDJOEL to calculate the rotational constants.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          CALL SCHRQ(I,AFLAG,GV(I),FWHM,PMAX,VLIM(ISTATE),V,SWF,BFCT,
     1             EPS(ISTATE),RMIN(ISTATE),RH(ISTATE),NDATPT(ISTATE),
     2                             NBEG,NEND,INNODE,INNR(I),IWR,LPRWF)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Call subroutine CDJOEL to determine the rotational constants.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          CALL CDJOEL(GV(I),NBEG,NEND,BvWN,RH(ISTATE),WARN,V,SWF,RM2,
     1                                                          RCNST)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Store molecular constants in array   ZK(v,J,isotope,state)
c
          ZK(I,0,IISTP,ISTATE)= GV(I)
          DO  J= 1,NROTMX
              ZK(I,J,IISTP,ISTATE)= RCNST(J)
              ENDDO
          WRITE(7,616) I,GV(I),(RCNST(J),J=1,NROTMX)
          ENDDO
      J=1
c
c** If all is well, then (without further ado) continue with the
c   calculations.
c
      RETURN
c-----------------------------------------------------------------------
  600 FORMAT(/'  *** INITDD ERROR ***',/4X,'For state',I3,
     1'  of isotope',I3,'  a serious error has occured:')
  601 FORMAT(4X,'Highest found level is less than desired (v=',I3,',J='
     1,I3,')')
  602 FORMAT(4X,'The Schrodinger Solver was unable to use the initial tr
     1ial energy.')
  603 FORMAT(4X,'The Schrodinger Solver was unable to use the calculated
     1 trial energy')
  604 FORMAT(4X,'The Automatic Level Finder could not find the first vib
     1rational level.')
  606 FORMAT(4X,'The next calculated trial energy was too low.')
  608 FORMAT(4X,'The next calculated trial energy was too high.')
  610 FORMAT(4X,'A second minimum exists in this potential')
  612 FORMAT(4X,'Could not find vibrational levels beyond (v=',I3,')')
  614 FORMAT(/'   v       ','E',9x,'Bv',10x,'-Dv',12x,'Hv',12x,'Lv',
     1  12x,'Mv',12x,'Nv',12x,'Ov'/1x,58('=='))
  616 FORMAT(I4,f12.4,f12.8,6(1PD14.6))
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE DEDP(IDAT,ISTATE,IISTP,ZMU,KVLEV,JROT,efPARITY,EO,
     1                                                 FWHM,ZK,DEDPK)
c***********************************************************************
c** This subroutine calculates  dE/dp from the expectation values of the
c   partial derivatives of an analytical potential with respect to its
c   parameters {p(k)}:    dE/dp(k) = <vj|dV/dp(k)|vj>  stored in DEDPK.
c
c** On entry:
c    ISTATE  is the molecular state being considered.
c    IISTP   is the isotopomer being considered.
c    ZMU     is the reduced mass of the diatom in atomic units.
c    KVLEV   is the vibrational quantum number.
c    JROT    is the rotational quantum number.
c    EO      is the initial trial energy (DE in cm-1).
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
c          A lead "card" identifies the level, gives the position of
c          1-st point and radial mesh, & states No. of  points.
c=======================================================================
      INCLUDE 'arrsizes.h'
      INCLUDE 'BLKISOT.h'
      INCLUDE 'BLKDATA.h'
      INCLUDE 'BLKPOT.h'
      INCLUDE 'BLKDVDP.h'
      INCLUDE 'BLKBOB.h'
      INCLUDE 'BLKBOBRF.h'
      INCLUDE 'BLKCOUNT.h'
c=======================================================================
      INTEGER IDAT, bandN, efPARITY, FIX, LAMB2, JRe
c
      INTEGER I,ICOR,ISTATE,IISTP,INNER, J,JROT,JIN, K,KVLEV,KV, 
     1  NBEG,NEND,INNODE,INNR(0:NVIBMX), IWR, LPRWF
c
      REAL*8 ZMU, EO, BFCT, JFCT, FWHM, VMAX, RM2, S0,
     1  ETRY,EE,SUM, SWF2,JFCTA,JFCTB, DUARe,DUBRe,DTARe,DTBRe,DLDRe,
     2   GV(0:NVIBMX),Vtotal(NPNTMX),JFCTDBL,JFCTD, muFCT,
     3   SWF(NPNTMX), V(NPNTMX), DEDPK(HPARMX)
c
      COMMON /VBLIK/Vtotal
c
c** Vibrational energies and rotational constants.
c
      REAL*8 ZK(0:NVIBMX,0:NROTMX,NISTPMX,NSTATEMX)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c** Initializing the arrays and variables
c
      DATA IWR/0/,LPRWF/0/,INNODE/1/
c
c** To calculate the values for <vj|dV/dP(k)|vj>
c   we must first determine the wave equation for the vj state.
      muFCT= 16.85762920d0/ZMU
      BFCT= RH(ISTATE)*RH(ISTATE)/muFCT
      JFCT= DBLE(JROT*(JROT+1))
      IF(IOMEG(ISTATE).GT.0) JFCT= JFCT - 
     1                               DBLE(IOMEG(ISTATE)*IOMEG(ISTATE))
      DO  J= 1,HPARMX
          DEDPK(J)= 0.d0
          ENDDO
c** Calculating the trial energy value from the band constants
      ETRY= ZK(KVLEV,0,IISTP,ISTATE)
      IF(JROT.GT.0) THEN
          DO I=1,NROTMX
              ETRY= ETRY + ZK(KVLEV,I,IISTP,ISTATE) * JFCT**I
              ENDDO
          ENDIF
      IF(IOMEG(ISTATE).GT.0) THEN 
c** For Lambda doubling, prepare rotational/mass factors
          LAMB2= 2*IOMEG(ISTATE)
          JFCTDBL= (efPARITY-EFREF(ISTATE))  
     1               * (DBLE(JROT*(JROT+1)) * muFCT**2)**IOMEG(ISTATE)
          ENDIF
      IF(IOMEG(ISTATE).EQ.-2) THEN
c** For doublet Sigma splitting, prepare rotational/mass factors
          LAMB2= 1
          IF(efPARITY.GT.0) JFCTDBL=  0.5d0*DBLE(JROT)*muFCT
          IF(efPARITY.EQ.0) JFCTDBL=  0.d0
          IF(efPARITY.LT.0) JFCTDBL= -0.5d0*DBLE(JROT+1)*muFCT
          ENDIF
c
c** Generating potential function including centrifugal, adiabatic BOB,
c   and rotational non-adiabatic BOB, and if appropriate, Lambda 
c   doubling or 2\Sigma doubling radial functions.
      JFCT= JFCT* RH(ISTATE)* RH(ISTATE)
      JFCTD= JFCTDBL* RH(ISTATE)* RH(ISTATE)/muFCT
      bandN= IB(IDAT)
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c!!!! Special modification for BOB  s0 = -0.8449 [GaH] for atom-2
ccc  Coxon UNrounded value for GaH:  s0(H)  = -0.844923030354E+00 
      FIX= 0
ccc36 S0= -0.844477105993d0* 5.4857990945d-04/ZMASS(2,IISTP)
ccc36 S0= -0.8d0* 5.4857990945d-04/ZMASS(2,IISTP)
   36 S0= 0.d0
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DO I= 1,NDATPT(ISTATE)
          RM2= 1/RD(I,ISTATE)**2
          V(I)= BFCT*(VPOT(I,ISTATE) + ZMUA(IISTP,ISTATE)*UAR(I,ISTATE)
     1                              + ZMUB(IISTP,ISTATE)*UBR(I,ISTATE))
     2          + JFCT* (1.0d0 + ZMTA(IISTP,ISTATE)*TAR(I,ISTATE)
     1                         + ZMTB(IISTP,ISTATE)*TBR(I,ISTATE))*RM2
          IF(IOMEG(ISTATE).NE.0) THEN
c ... note that power of 1/r**2  included in  wRAD  array ...
              V(I)= V(I) + JFCTD * wRAD(I,ISTATE)
              ENDIF
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c!!!! Special modification for BOB  s0 = -0.8449 [GaH] for atom-2
          IF(DABS(S0).GT.0.d0) V(I)= V(I) -
ccc  1     S0*(BFCT*(VPOT(I,ISTATE) - ETRY) + JFCT*RM2)
     1     S0*(1.d0- S0+ S0**2)*(V(I)- BFCT*ETRY)
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          Vtotal(I)= V(I)/BFCT
          ENDDO
      KV= KVLEV
      INNER= 0
      EO= ETRY
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      CALL SCHRQ(KV,JROT,EO,FWHM,VMAX,VLIM(ISTATE),V,SWF,BFCT,
     1             EPS(ISTATE),RMIN(ISTATE),RH(ISTATE),NDATPT(ISTATE),
     2                               NBEG,NEND,INNODE,INNER,IWR,LPRWF)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c!!!! Special Case ... modification for BOB  s0  for atom-2
      IF(DABS(S0).GT.0.d0) THEN
          EO= EO + (EO-ETRY)*S0
          IF(DABS(EO-ETRY).GT.1.d-04) THEN
              ETRY= EO
              FIX= FIX+1
              IF(FIX.LT.20) GOTO 36
cc          else
cc            write(33,333) kv,jrot,eO
cc333 format('   v=',i2,'   J=',i3,'   E=',f10.4)
              ENDIF
          ENDIF
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ICOR= 0
      IF(KV.NE.KVLEV) THEN
c** If the calculated wavefunction is for the wrong vibrational level,
c .... First try to estimate correct energy using WKB energy derivative
   38     ICOR= ICOR+ 1
          EE= EO*BFCT
          DO  I= NEND,1,-1
              K= I
              IF(V(I).LT.EE) GOTO 40
              ENDDO
   40     SUM= 0.d0
          DO  I= K,1,-1
              IF(V(I).GT.EE) GO TO 42
              SUM= SUM+ 1.d0/DSQRT(EE - V(I))
              ENDDO
   42     ETRY= EO - DBLE(KV-KVLEV)*6.28318531d0/(BFCT*SUM)
          EO= ETRY
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          CALL SCHRQ(KV,JROT,EO,FWHM,VMAX,VLIM(ISTATE),V,SWF,BFCT,
     1             EPS(ISTATE),RMIN(ISTATE),RH(ISTATE),NDATPT(ISTATE),
     2                               NBEG,NEND,INNODE,INNER,IWR,LPRWF)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          IF((KV.NE.KVLEV).AND.(ICOR.LE.2)) GOTO 38
          ENDIF
      IF(KV.NE.KVLEV) THEN
c ... if that fails, next try finding the level in the inner well 
          WRITE(6,600) ISTATE, IISTP
          WRITE(6,610) KVLEV, JROT
          IF(KV.GE.0) THEN
              WRITE(6,620) KV, JROT, EO, ETRY
              KV= KVLEV
              INNER= 1
              EO= ETRY
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
              CALL SCHRQ(KV,JROT,EO,FWHM,VMAX,VLIM(ISTATE),V,SWF,BFCT,
     1             EPS(ISTATE),RMIN(ISTATE),RH(ISTATE),NDATPT(ISTATE),
     2                               NBEG,NEND,INNODE,INNER,IWR,LPRWF)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            ELSE
              WRITE(6,621)
            ENDIF
          ENDIF
c ... If that also fails, then try calling ALF to recalculate the
c   vibrational energy from that effective radial potential.
      IF(KV.NE.KVLEV) THEN
          IF(INNER.EQ.1) THEN
              WRITE(6,611) KVLEV, JROT
              IF (KV.GE.0) THEN
                  WRITE(6,620) KV, JROT, EO, ETRY
                ELSE
                  WRITE(6,621)
                ENDIF
              ENDIF
          WRITE(6,612)
          GV(KVLEV)= 0.0d0
          KV= KVLEV
          JIN= JROT
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          CALL ALF(NDATPT(ISTATE),RMIN(ISTATE),RH(ISTATE),V,SWF,
     1    VLIM(ISTATE),KV,JIN,ZMU,EPS(ISTATE),MMLR(1,ISTATE),GV,INNODE,
     2                                                     INNR,IWR)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          KV= KVLEV
          ETRY= GV(KVLEV)
          EO= ETRY
c
c** Call SCHRQ again with the new calculated energy.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          CALL SCHRQ(KV,JROT,EO,FWHM,VMAX,VLIM(ISTATE),V,SWF,BFCT,
     1             EPS(ISTATE),RMIN(ISTATE),RH(ISTATE),NDATPT(ISTATE),
     2                         NBEG,NEND,INNODE,INNR(KVLEV),IWR,LPRWF)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          ENDIF
c
c** If the calculated wavefunction is still for the wrong vibrational
c   level, then write out a warning and skip the calculation of the
c   partial derivatives (hence setting them to zero).
c
      IF(KV.NE.KVLEV) THEN
          IF (INNER.EQ.1) THEN
              WRITE(6,611) KVLEV, JROT
              IF (KV.GE.0) THEN
                  WRITE(6,620) KV, JROT, EO, ETRY
                ELSE
                  WRITE(6,621)
                ENDIF
              ENDIF
c.. an eigenvalue of  -9.9d9 signifies that eigenvalue search failed
          EO= -9.9d9
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** If the calculated wavefunction is for the right vibrational level,
c   then continue with calculation of the partial derivatives by
c   integration using the modified trapezoidal rule.
c
c** First determine which partials need to be changed, increment so
c   that the (fluorescence term values and) lower states are skipped.
c
        ELSE
          JRe= POTPARI(ISTATE)+ 1
          IF(PSEL(ISTATE).EQ.4) JRe= POTPARI(ISTATE)
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
     1          + DUBRe*ZMUB(IISTP,ISTATE) + JFCTA*DTARe + JFCTB*DTBRe
     2          + JFCTDBL*DLDRE
          DO  J= POTPARI(ISTATE), POTPARF(ISTATE)
              DEDPK(J)= DEDPK(J)* RH(ISTATE)
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
                  DEDPK(J)= DEDPK(J)* ZMUA(IISTP,ISTATE)* RH(ISTATE)
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
          ENDIF
      RETURN
c-----------------------------------------------------------------------
  600 FORMAT(/'  *** DEDP ERROR ***',/4X,'For state',I3,
     1'  of isotope',I3,'  a serious error has occured:')
  610 FORMAT(4X,'Unable to find v=',I3,' J=',I3,' level in outer well.')
  611 FORMAT(4X,'Unable to find v=',I3,' J=',I3,' level in inner well.')
  612 FORMAT(4X,'Molecular constants give improper trial energy.')
  620 FORMAT('   SCHRQ found   E(v=',I3,',J=',I3,')=',F10.3,'  with  E(t
     1rial)=',F10.3) 
  621 FORMAT(4X,'SCHRQ was unable to calculate any level.')
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
c    IISTP   is the isotopomer being considered.
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
      FCTOR= DSQRT(ZMU/16.85762920d0)
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
      INTEGER  I,IDAT,ISTATE,vb,Jr,index,NCNN
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
              CALL PREPOT(0,AN(1),AN(2),MN(1,1),MN(2,1),1,RR,RM2,
     1                                                  VLIMT,VV,NCNN)
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
c ... obtain potential function value at new turning point estimate
              CALL VGEN(ISTATE,R2,Vmid,bMi,IDAT)
            ELSE
              RR(1)= R2
              RM2(1)= 1.0d0/R2**2
c ... for fixed pointwise potential, interpolate for potential value
              CALL PREPOT(0,AN(1),AN(2),MN(1,1),MN(2,1),1,RR,RM2,
     1                                                  VLIMT,VV,NCNN)
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
              CALL PREPOT(0,AN(1),AN(2),MN(1,1),MN(2,1),1,RR,RM2,
     1                                                  VLIMT,VV,NCNN)
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
      SUBROUTINE PhaseIntegral(IDAT,ISTATE,Ri,Ro,k1,EO,DEDPK,Int)
c=======================================================================
c  Subroutine PhaseIntegral calculates phase integrals of the integrand
c  f(R)/[E-V(R)]**(k+1/2), for k=-1,0 or 1  on interval between turning 
c  points  Ri and Ro, where f(R) is either 1 or a derivative of the
c  potential w.r.t. some potential parameter. 
c             Ro                                         N
c      Int = Int dr f(r)/|E-V(r)|^{k+1/2} = ½(Ro - Ri) Sum{Wi*F(Zi)}
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
c             k       is the powers in the phase integral Int
c             n=0  for  f(r)= 1;  n=1  for  f(r)= [dV(r)/dp_k - dE/dp_k]
c             DEDPK(i) are the values of the partial derivative dE/dp_k.
c* On exit:  Int(j,k)    are the phase integral(s)
c=======================================================================
      INCLUDE 'arrsizes.h'
      INCLUDE 'BLKCOUNT.h'
c** Type statements & common block for data
      INTEGER i,IPV,ISTATE, k,k1,kmx,M,IDAT
c** M  is the number of quadrature mesh points
      PARAMETER (M= 21)
c
      REAL*8 dVdPk(HPARMX)
      COMMON /dVdPkBLK/dVdPk
c
      REAL*8 Ri,Ro,EO,Fzt,RDIST,VDIST,PHIDIST,EMinusV,
     1  Int(0:HPARMX,-1:1),DEDPK(HPARMX)
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
      IF(k1 .eq. -1) kmx= 0
c???     IF( NTFP(ISTATE) .GT. 0 ) kmx= 1
      kmx= 1
      DO  k= k1,kmx
          Int(0,k)= 0.d0
          DO  IPV= POTPARI(ISTATE),HPARF(ISTATE)
              Int(IPV,k)= 0.0d0
              ENDDO
          ENDDO
c
c*** Begin quadrature loop for sums over M mesh points
      DO  i= 1,M
c ... first get potential and derivatives at Pajunen k=1 point
          RDIST= 0.5d0*(Ro+Ri + (Ro-Ri)*Zp21(i))
          CALL VGEN(ISTATE,RDIST,VDIST,PHIDIST,IDAT)
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
              Int(0,k)= Int(0,k) - Wp21k1(i)* Fzt
ccc           IF(k .gt. k1) THEN
                  DO IPV= POTPARI(ISTATE),HPARF(ISTATE)
c ... Loop over free parameters accumulating phase integral derivs.
                      Int(IPV,k)= Int(IPV,k) - Wp21k1(i)* Fzt
     1                                       *(dVdPk(IPV)- DEDPK(IPV))
                      ENDDO
ccc               ENDIF
c .... end of loop over  k values
              ENDDO
c ... end of quadrature loop over mesh points ....................
          ENDDO
      DO  k= k1,kmx
          Int(0,k)= 0.25d0*(Ro-Ri) * Int(0,k)
          DO  IPV= POTPARI(ISTATE),HPARF(ISTATE)
c*** In Pajunen's method, a factor of (1/2) is incorporated because 
c    the phase integral is a contour integral
              Int(IPV,k)= 0.25d0*(Ro-Ri) * Int(IPV,k)
              ENDDO
          ENDDO
c????????????????????????
c** Test that 'regular' & Pajunen integrals agree for k= -1 & 0
cc  and that Pajunen  k=1 & 3 integrals agree for k=1 case
cc    DO  IPV= POTPARI(ISTATE),HPARF(ISTATE)
cc        WRITE(32, 320) EO, (Int(IPV,k), k= k1, kmx)
cc        WRITE(32, 321)  (Int(IPV,k), k= k1, kmx)
cc        Do  k= k1, kmx
cc            tst(k)= 0.d0
cc            IF(DABS(Ink(IPV,k)) .gt. 0.d0) 
cc   1                              tst(k)= Int(IPV,k)/Ink(IPV,k)-1.d0
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

c**********************************************************************
      SUBROUTINE DIFFSTATS(NSTATES,ROBUST,MKPRED)
c** Subroutine to summarise dimensionless standard errors on a band-by-
c  band basis, and (if desired) print [obs.-calc.] values to channel-8.
c-----------------------------------------------------------------------
c                 Version of  XX February 2006
c-----------------------------------------------------------------------
      INCLUDE 'arrsizes.h'
      INCLUDE 'BLKISOT.h'
      INCLUDE 'BLKDATA.h'
      INCLUDE 'BLKTYPE.h'
c
      INTEGER I,IBB,ISOT,ISTATE,ISTATEE,J,K,NSTATES,MKPRED,ROBUST
      REAL*8 AVE,AVETOT,DIV,RMSR,RMSTOT,SSQTOT
      CHARACTER*3 MARKER,NEF(-1:1)
c
      DATA  NEF/'  f',' ef','  e'/
c========================================================================
      ISOT= 1
      SSQTOT= 0.d0
      IF(MKPRED.GT.0) THEN
          WRITE(6,600)
          REWIND 4
          WRITE(4,601)
          ENDIF
c** Summarize data discrepancies for one isotopomer at a time.
   10 WRITE(6,602) NBANDS(ISOT),(NAME(I),MN(I,ISOT),I= 1,2)
c
c** Loop over bands for each (lower) electronic state, in turm
      DO 90 ISTATE= 1,NSTATES
      IF(NTRANSMW(ISOT,ISTATE).GT.0)THEN
c** Book-keeping for Micowave data
          WRITE(6,604) NTRANSMW(ISOT,ISTATE),SLABL(ISTATE),(NAME(I),
     1                         MN(I,ISOT),I= 1,2),NBANDMW(ISOT,ISTATE)
          WRITE(6,605)
          WRITE(8,604) NTRANSMW(ISOT,ISTATE),
     1  SLABL(ISTATE),(NAME(I),MN(I,ISOT),I= 1,2),NBANDMW(ISOT,ISTATE)
          RMSTOT= 0.d0
          AVETOT= 0.d0
          DO  I= 1,NBANDMW(ISOT,ISTATE)
              IBB= YPR(ISOT,ISTATE,4,4,I)
              IF(MKPRED.LE.0) THEN
                  CALL BNDERR(IFIRST(IBB),ILAST(IBB),ROBUST,AVE,RMSR,
     1                                             SSQTOT,DFREQ,UFREQ)
                  RMSTOT= RMSTOT+ YPR(ISOT,ISTATE,4,3,I)*RMSR**2
                  AVETOT= AVETOT+ YPR(ISOT,ISTATE,4,3,I)*AVE
                  WRITE(6,606)YPR(ISOT,ISTATE,4,2,I),
     1                  YPR(ISOT,ISTATE,4,1,I),YPR(ISOT,ISTATE,4,3,I),
     2                  YPR(ISOT,ISTATE,4,5,I),YPR(ISOT,ISTATE,4,6,I),
     3                  AVEUFREQ(IBB),MAXUFREQ(IBB),AVE,RMSR
                  ENDIF
              WRITE(8,605)
              IF(MKPRED.LE.0) WRITE(8,606) YPR(ISOT,ISTATE,4,2,I),
     1                  YPR(ISOT,ISTATE,4,1,I),YPR(ISOT,ISTATE,4,3,I),
     2                  YPR(ISOT,ISTATE,4,5,I),YPR(ISOT,ISTATE,4,6,I),
     4                  AVEUFREQ(IBB),MAXUFREQ(IBB),AVE,RMSR
              IF(MKPRED.GT.0) THEN
                  WRITE(8,606) YPR(ISOT,ISTATE,4,2,I),
     1                  YPR(ISOT,ISTATE,4,1,I),YPR(ISOT,ISTATE,4,3,I),
     2                  YPR(ISOT,ISTATE,4,5,I),YPR(ISOT,ISTATE,4,6,I)
                  WRITE(4,640) YPR(ISOT,ISTATE,4,2,I),
     1             YPR(ISOT,ISTATE,4,1,I),SLABL(ISTATE),SLABL(ISTATE),
     2                                              (MN(K,ISOT),K=1,2)
                  ENDIF
  640 FORMAT(/2I4,2(2x,"'",A2,"'"),2x,2I4)
              CALL PBNDERR(IBB,MKPRED,NEF)
              ENDDO
          RMSTOT= DSQRT(RMSTOT/NTRANSMW(ISOT,ISTATE))
          AVETOT= AVETOT/NTRANSMW(ISOT,ISTATE)
          IF(MKPRED.LE.0) WRITE(6,630) NTRANSMW(ISOT,ISTATE),AVETOT,
     1                                                          RMSTOT
	    ENDIF
c
      IF(NTRANSIR(ISOT,ISTATE).GT.0)THEN
c** Book-keeping for Infrared data
          WRITE(6,608) NTRANSIR(ISOT,ISTATE),SLABL(ISTATE),(NAME(I),
     1                         MN(I,ISOT),I= 1,2),NBANDIR(ISOT,ISTATE)
          WRITE(6,605)
          WRITE(8,608) NTRANSIR(ISOT,ISTATE),
     1  SLABL(ISTATE),(NAME(I),MN(I,ISOT),I= 1,2),NBANDIR(ISOT,ISTATE)
          RMSTOT= 0.d0
          AVETOT= 0.d0
          DO  I= 1,NBANDIR(ISOT,ISTATE)
              IBB= YPR(ISOT,ISTATE,3,4,I)
              IF(MKPRED.LE.0) THEN
                  CALL BNDERR(IFIRST(IBB),ILAST(IBB),ROBUST,AVE,RMSR,
     1                                             SSQTOT,DFREQ,UFREQ)
                  RMSTOT= RMSTOT+ YPR(ISOT,ISTATE,3,3,I)*RMSR**2
                  AVETOT= AVETOT+ YPR(ISOT,ISTATE,3,3,I)*AVE
                  WRITE(6,606) YPR(ISOT,ISTATE,3,2,I),
     1                  YPR(ISOT,ISTATE,3,1,I),YPR(ISOT,ISTATE,3,3,I),
     2                  YPR(ISOT,ISTATE,3,5,I),YPR(ISOT,ISTATE,3,6,I),
     3                  AVEUFREQ(IBB),MAXUFREQ(IBB),AVE,RMSR
                  ENDIF
              WRITE(8,605)
              IF(MKPRED.LE.0) WRITE(8,606) YPR(ISOT,ISTATE,3,2,I),
     1                  YPR(ISOT,ISTATE,3,1,I),YPR(ISOT,ISTATE,3,3,I),
     2                  YPR(ISOT,ISTATE,3,5,I),YPR(ISOT,ISTATE,3,6,I),
     3                  AVEUFREQ(IBB),MAXUFREQ(IBB),AVE,RMSR
              IF(MKPRED.GT.0) THEN
                  WRITE(8,606) YPR(ISOT,ISTATE,3,2,I),
     1                  YPR(ISOT,ISTATE,3,1,I),YPR(ISOT,ISTATE,3,3,I),
     2                  YPR(ISOT,ISTATE,3,5,I),YPR(ISOT,ISTATE,3,6,I)
                  WRITE(4,640) YPR(ISOT,ISTATE,3,2,I),
     1             YPR(ISOT,ISTATE,3,1,I),SLABL(ISTATE),SLABL(ISTATE),
     2                                              (MN(K,ISOT),K=1,2)
                  ENDIF
              CALL PBNDERR(IBB,MKPRED,NEF)
              ENDDO
          RMSTOT= DSQRT(RMSTOT/NTRANSIR(ISOT,ISTATE))
          AVETOT= AVETOT/NTRANSIR(ISOT,ISTATE)
          IF(MKPRED.LE.0) WRITE(6,630) NTRANSIR(ISOT,ISTATE),AVETOT,
     1                                                          RMSTOT
          ENDIF
c
c** Book-keeping for Electronic vibrational band data
      DO  ISTATEE= 1,NSTATES
          IF((ISTATEE.NE.ISTATE).AND.
     1                 (NTRANSVIS(ISOT,ISTATEE,ISTATE).GT.0)) THEN
c ... for ISTATEE{upper}-ISTATE{lower} electronic vibrational bands
              WRITE(6,610) NTRANSVIS(ISOT,ISTATEE,ISTATE),
     1                      (NAME(I),MN(I,ISOT),I=1,2),SLABL(ISTATEE),
     2                      SLABL(ISTATE),NBANDEL(ISOT,ISTATEE,ISTATE)
              WRITE(6,605)
              WRITE(8,610) NTRANSVIS(ISOT,ISTATEE,ISTATE),
     1                      (NAME(I),MN(I,ISOT),I=1,2),SLABL(ISTATEE),
     2                      SLABL(ISTATE),NBANDEL(ISOT,ISTATEE,ISTATE)
              RMSTOT= 0.d0
              AVETOT= 0.d0
              DO  I= 1,NBANDVIS(ISOT,ISTATE)
                  IBB= YPR(ISOT,ISTATE,2,4,I)
                  IF(IEP(IBB).EQ.ISTATEE) THEN
                      IF(MKPRED.LE.0) THEN
                          CALL BNDERR(IFIRST(IBB),ILAST(IBB),ROBUST,AVE,
     1                                        RMSR,SSQTOT,DFREQ,UFREQ)
                          RMSTOT= RMSTOT+ YPR(ISOT,ISTATE,2,3,I)*RMSR**2
                          AVETOT= AVETOT+ YPR(ISOT,ISTATE,2,3,I)*AVE
                          WRITE(6,606) YPR(ISOT,ISTATE,2,2,I),
     1                  YPR(ISOT,ISTATE,2,1,I),YPR(ISOT,ISTATE,2,3,I),
     2                  YPR(ISOT,ISTATE,2,5,I),YPR(ISOT,ISTATE,2,6,I),
     3                            AVEUFREQ(IBB),MAXUFREQ(IBB),AVE,RMSR
                          ENDIF
                      WRITE(8,605)
                      IF(MKPRED.LE.0) WRITE(8,606) 
     1                  YPR(ISOT,ISTATE,2,2,I),YPR(ISOT,ISTATE,2,1,I),
     2                  YPR(ISOT,ISTATE,2,3,I),YPR(ISOT,ISTATE,2,5,I),
     3             YPR(ISOT,ISTATE,2,6,I),AVEUFREQ(IBB),MAXUFREQ(IBB),
     4                                                        AVE,RMSR
                      IF(MKPRED.GT.0) THEN
                           WRITE(8,606) YPR(ISOT,ISTATE,2,2,I),
     1                  YPR(ISOT,ISTATE,2,1,I),YPR(ISOT,ISTATE,2,3,I),
     2                  YPR(ISOT,ISTATE,2,5,I),YPR(ISOT,ISTATE,2,6,I)
                           WRITE(4,640) YPR(ISOT,ISTATE,2,2,I),
     1             YPR(ISOT,ISTATE,2,1,I),SLABL(ISTATE),SLABL(ISTATE),
     2                                              (MN(K,ISOT),K=1,2)
                  ENDIF
                      CALL PBNDERR(IBB,MKPRED,NEF)
                      ENDIF
                  ENDDO
              RMSTOT= DSQRT(RMSTOT/NTRANSVIS(ISOT,ISTATEE,ISTATE))
              AVETOT= AVETOT/NTRANSVIS(ISOT,ISTATEE,ISTATE)
              IF(MKPRED.LE.0) WRITE(6,630) 
     1                    NTRANSVIS(ISOT,ISTATEE,ISTATE),AVETOT,RMSTOT
              ENDIF
          ENDDO
c
      IF(NTRANSFS(ISOT,ISTATE).GT.0)THEN
c** Book-keeping for Fluorescence data
          WRITE(6,612) NTRANSFS(ISOT,ISTATE),SLABL(ISTATE),(NAME(I),
     1                          MN(I,ISOT),I=1,2),NBANDFS(ISOT,ISTATE)
          WRITE(6,617)
          WRITE(8,612) NTRANSFS(ISOT,ISTATE),SLABL(ISTATE),
     1                 (NAME(I),MN(I,ISOT),I=1,2),NBANDFS(ISOT,ISTATE)
          RMSTOT= 0.d0
          AVETOT= 0.d0
          DO  I= 1,NBANDFS(ISOT,ISTATE)
              IBB= YPR(ISOT,ISTATE,1,4,I)
              CALL BNDERR(IFIRST(IBB),ILAST(IBB),ROBUST,AVE,RMSR,
     1                                             SSQTOT,DFREQ,UFREQ)
              RMSTOT= RMSTOT+ YPR(ISOT,ISTATE,1,3,I)*RMSR**2
              AVETOT= AVETOT+ YPR(ISOT,ISTATE,1,3,I)*AVE
              WRITE(6,614) YPR(ISOT,ISTATE,1,1,I),
     1                   YPR(ISOT,ISTATE,1,2,I),NEF(EFP(IFIRST(IBB))),
     2                  YPR(ISOT,ISTATE,1,3,I),YPR(ISOT,ISTATE,1,5,I),
     3                  YPR(ISOT,ISTATE,1,6,I),
     4                  AVEUFREQ(IBB),MAXUFREQ(IBB),AVE,RMSR
              WRITE(8,617)
              WRITE(8,614) YPR(ISOT,ISTATE,1,1,I),
     1                  YPR(ISOT,ISTATE,1,2,I),NEF(EFP(IFIRST(IBB))),
     2                  YPR(ISOT,ISTATE,1,3,I),YPR(ISOT,ISTATE,1,5,I),
     3                  YPR(ISOT,ISTATE,1,6,I),
     4                  AVEUFREQ(IBB),MAXUFREQ(IBB),AVE,RMSR
              CALL PBNDERR(IBB,MKPRED,NEF)
              ENDDO
              RMSTOT= DSQRT(RMSTOT/NTRANSFS(ISOT,ISTATE))
              AVETOT= AVETOT/NTRANSFS(ISOT,ISTATE)
              WRITE(6,632) NTRANSFS(ISOT,ISTATE),AVETOT,RMSTOT
          ENDIF
c
      IF(NEBPAS(ISOT,ISTATE).GT.0) THEN
c** Book-keeping for  PAS  data
          IBB= YPR(ISOT,ISTATE,7,4,1)
          CALL BNDERR(IFIRST(IBB),ILAST(IBB),ROBUST,AVE,RMSR,SSQTOT,
     1                                                    DFREQ,UFREQ)
          WRITE(6,626) NEBPAS(ISOT,ISTATE),SLABL(ISTATE),(NAME(I),
     1 MN(I,ISOT),I=1,2),YPR(ISOT,ISTATE,7,3,1),YPR(ISOT,ISTATE,7,5,1),
     2      YPR(ISOT,ISTATE,7,6,1),AVEUFREQ(IBB),MAXUFREQ(IBB),AVE,RMSR
          WRITE(8,626) NEBPAS(ISOT,ISTATE),SLABL(ISTATE),(NAME(I),
     1 MN(I,ISOT),I=1,2),YPR(ISOT,ISTATE,7,3,1),YPR(ISOT,ISTATE,7,5,1),
     2      YPR(ISOT,ISTATE,7,6,1),AVEUFREQ(IBB),MAXUFREQ(IBB),AVE,RMSR
          WRITE(8,627)
          DO  I= IFIRST(IBB),ILAST(IBB)
              DIV= DABS(DFREQ(I)/UFREQ(I))
              marker='   '
              IF( (DIV.GE.2.d0).AND.(DIV.LT.5.d0) ) marker='*  '
              IF( (DIV.GE.5.d0).AND.(DIV.LT.10.d0) ) marker='** '
              IF( (DIV.GE.10.d0) ) marker='***'
              WRITE(8,628) JP(I),JPP(I),NEF(EFPP(I)),FREQ(I),
     1                      UFREQ(I),DFREQ(I),DFREQ(I)/UFREQ(I),MARKER
              ENDDO
          WRITE(6,629)
          WRITE(8,629)
          ENDIF
c
      IF(NBVPP(ISOT,ISTATE).GT.0) THEN
c** Book-keeping for  Bv  data
          IBB= YPR(ISOT,ISTATE,5,4,1)
          CALL BNDERR(IFIRST(IBB),ILAST(IBB),ROBUST,AVE,RMSR,SSQTOT,
     1                                                    DFREQ,UFREQ)
          WRITE(6,616) NBVPP(ISOT,ISTATE),SLABL(ISTATE),
     1              (NAME(I),MN(I,ISOT),I=1,2),YPR(ISOT,ISTATE,5,3,1),
     2              YPR(ISOT,ISTATE,5,5,1),YPR(ISOT,ISTATE,5,6,1),
     3              AVEUFREQ(IBB),MAXUFREQ(IBB),AVE,RMSR
          WRITE(8,616) NBVPP(ISOT,ISTATE),SLABL(ISTATE),
     1              (NAME(I),MN(I,ISOT),I=1,2),YPR(ISOT,ISTATE,5,3,1),
     2              YPR(ISOT,ISTATE,5,5,1),YPR(ISOT,ISTATE,5,6,1),
     3              AVEUFREQ(IBB),MAXUFREQ(IBB),AVE,RMSR
          DO  J= IFIRST(IBB),ILAST(IBB)
              WRITE(6,618) JP(J),NEF(EFPP(J)),FREQ(J),UFREQ(J),
     1                                      DFREQ(J),DFREQ(J)/UFREQ(J)
              WRITE(8,618) JP(J),NEF(EFPP(J)),FREQ(J),UFREQ(J),
     1                                      DFREQ(J),DFREQ(J)/UFREQ(J)
              ENDDO
          ENDIF
c
      IF(NWIDTH(ISOT,ISTATE).GT.0) THEN
c** Book-keeping for  Tunneling Width  data
          IBB= YPR(ISOT,ISTATE,6,4,1)    
          CALL BNDERR(IFIRST(IBB),ILAST(IBB),ROBUST,AVE,RMSR,SSQTOT,
     1                                                    DFREQ,UFREQ)
          WRITE(6,620) NWIDTH(ISOT,ISTATE),SLABL(ISTATE),
     1              (NAME(I),MN(I,ISOT),I=1,2),YPR(ISOT,ISTATE,6,3,1),
     2              YPR(ISOT,ISTATE,6,5,1),YPR(ISOT,ISTATE,6,6,1),
     3              AVEUFREQ(IBB),MAXUFREQ(IBB),AVE,RMSR
          WRITE(8,620) NWIDTH(ISOT,ISTATE),SLABL(ISTATE),
     1              (NAME(I),MN(I,ISOT),I=1,2),YPR(ISOT,ISTATE,6,3,1),
     2              YPR(ISOT,ISTATE,6,5,1),YPR(ISOT,ISTATE,6,6,1),
     3              AVEUFREQ(IBB),MAXUFREQ(IBB),AVE,RMSR
          DO  J= IFIRST(IBB),ILAST(IBB)
              WRITE(6,622) JP(J),JPP(J),NEF(EFPP(J)),FREQ(J),UFREQ(J),
     1                                      DFREQ(J),DFREQ(J)/UFREQ(J)
              WRITE(8,622) JP(J),JPP(J),NEF(EFPP(J)),FREQ(J),
     1                             UFREQ(J),DFREQ(J),DFREQ(J)/UFREQ(J)
              ENDDO
          ENDIF
c** End of loop over the various (lower) electronic states
   90 CONTINUE
c=======================================================================
      IF(ISOT.LT.NISTP) THEN
c** If NISTP > 1, return to print data summaries for other isotopomers
          ISOT= ISOT+1
          GO TO 10
          ENDIF 
      RMSR= DSQRT(SSQTOT/COUNTOT)
      WRITE(6,624) COUNTOT,RMSR
      RETURN
  600 FORMAT(/1x,36('**')/'  Write to Channel-8 Predictions From Complet
     1e Set of Input Parameters!'/1x,36('**'))
  601 FORMAT(/1x,25('**')/'  Predictions From Complete Set of Input Para
     1meters!'/1x,25('**'))
  602 FORMAT(/1x,21('===')/'  *** Discrepancies for',I5,' bands/series o
     1f ',A2,'(',I3,')-',A2,'(',I3,') ***'/1x,21('==='))
  604 FORMAT(/1x,21('===')/I5,' State ',A2,1x,A2,'(',I3,')-',A2,'(',I3,
     1 ') MW transitions in',i4,' vib. levels')
  605 FORMAT(1x,16('==='),'== Avge. ========'/"   v' ",
     2  ' v" #data  J"min  J"max  Av.Unc.  Max.Unc.   Err/Unc   DRMSD'/
     1  1x,13('-----'))
  606 FORMAT(2I4,I6,3x,I4,3x,I4,1x,1P2D9.1,0PF11.5,F8.3)
  608 FORMAT(/1x,63('=')/I5,' State ',A2,1x,A2,'(',I3,')-',A2,'(',I3,
     1 ') InfraRed transitions in',I4,' bands')
  610 FORMAT(/1x,35('==')/I6,1x,A2,'(',I3,')-',A2,'(',i3,')  {State ',
     1  A2,'}--{State ',A2,'} Transitions in',i4,' bands')
  612 FORMAT(/1x,75('=')/I5,' Fluorescence transitions into State ',A2,
     1 2x,A2,'(',I3,')-',A2,'(',I3,') in',i5,' series')
  617 FORMAT(1x,52('='),'= Avge. ',15('=')/"   v'  j' p' ",
     2  '#data  v"min  v"max','  AvgeUnc  Max.Unc.   Err/Unc   DRMSD'/
     3  1x,25('---'))
  614 FORMAT(2I4,A3,I6,2I7,1x,1P2D9.1,0PF11.5,F8.3)
  616 FORMAT(/1x,66('=')/1x,I3,' State ',A2,1x,A2,'(',I3,')-',A2,'(',
     1 I3,') Bv values treated as independent data'/1x,20('=='),
     2 '  Avge.  ',17('=')/' #data   v"min  v"max  AvgeUnc  Max.Unc.  Er
     3r/Unc  DRMSD'/1x,55('-')/I5,2I7,2x,1P2D9.1,0PF9.3,F8.3/
     4 1x,30('==')/'    v  p',8x,'Bv',7x,'u(Bv)',4x,
     5 '[calc-obs]  [calc-obs]/unc',/1x,30('--'))
  618 FORMAT(I5,A3,2x,F12.8,1PD9.1,0PF13.8,F12.4)
  620 FORMAT(/1x,73('=')/1x,I3,' State ',A2,1x,A2,'(',I3,')-',A2,'(',
     1 I3,') Tunneling Widths treated as independent data'/1x,20('=='),
     2 '  Avge.  ',24('=')/' #data   v"min  v"max  AvgeUnc  Max.Unc.  Er
     3r/Unc  DRMSD'/1x,55('-')/I5,2I7,2x,1P2D9.1,0PF9.3,F8.3/
     4 1x,59('=')/'   v   J  p     Width',7x,'u(Width)  [calc-obs]  [cal
     5c-obs]/unc'/1x,59('-'))
  622 FORMAT(2I4,A3,1PD14.6,D10.1,D13.2,0PF10.3)
  624 FORMAT(/1x,29('==')/' For overall fit to',i6,' data,  DRMS(deviati
     1ons)=',G13.6/1x,30('==')) 
  626 FORMAT(/1x,29('==')/I5,' PAS Binding Energies for State ',A2,2x,
     1 A2,'(',I3,')-',A2,'(',I3,')'/1x,50('='),' Avge. ',('=')/
     2 ' #data  v_min  v_max   AvgeUnc  Max.Unc.  Err/Unc  DRMSD'/
     3  1x,29('--')/I5,2I7,2x,1P2D9.1,0PF9.3,F8.3)
  627 FORMAT(1x,48('='),'  calc-obs'/'   v   j  p      PAS(Eb)     u(Eb)
     1      calc-obs   /u(FREQ)'/1x,29('--'))
  628 FORMAT(2I4,A3,F14.6,1PD10.1,D13.2,0PF11.4,1X,A3)
  629 FORMAT(1x,29('=='))
  630 FORMAT(1x,7('--'),' For these',i6,' lines, overall:',F11.5,F8.3)
  632 FORMAT(1x,17('-'),' For these',i6,' lines, overall:',F11.5,F8.3)
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE BNDERR(FIRST,LAST,ROBUST,AVEDD,RMSDD,SSQTOT,DFREQ,
     1                                                          UFREQ)
c** Calculate the average (AVEDD) & the root mean square dimensionless 
c  deviation (RSMDD) for the band running from datum # FIRST to LAST.
      INCLUDE 'arrsizes.h'
c
      REAL*8 DFREQ(NDATAMX),UFREQ(NDATAMX),AVEDD,RMSDD,SSQTOT
      INTEGER FIRST,LAST,NDAT,I,ROBUST
c
      AVEDD= 0.d0
      RMSDD= 0.d0
      DO  I= FIRST,LAST
          AVEDD= AVEDD+ DFREQ(I)/UFREQ(I)
          IF(ROBUST.LE.0) RMSDD= RMSDD+ (DFREQ(I)/UFREQ(I))**2
          IF(ROBUST.GT.0) RMSDD= RMSDD+ DFREQ(I)**2/
     1                                (UFREQ(I)**2 + DFREQ(I)**2/3.d0)
          ENDDO
      SSQTOT= SSQTOT+ RMSDD
      NDAT= LAST-FIRST+1
      AVEDD= AVEDD/NDAT
      RMSDD= DSQRT(RMSDD/NDAT)
      RETURN
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE PBNDERR(IBB,MKPRED,NEF)
c** Print to channel-8 a listing of the [obs.-calc.] values for the band
c  running from datum # FIRST to LAST.                           
      INCLUDE 'arrsizes.h'             
      INCLUDE 'BLKDATA.h'
      REAL*8 DIV
      INTEGER IBB,I,MKPRED
      CHARACTER*3 marker, NEF(-1:1)
c----------------------------------------------------------------------- 
      IF(MKPRED.LE.0) WRITE(8,600)
      IF(MKPRED.GT.0) WRITE(8,601)
      DO  I= IFIRST(IBB),ILAST(IBB)
          IF(MKPRED.LE.0) THEN
              DIV= DABS(DFREQ(I)/UFREQ(I))
              marker='   '
              IF( (DIV.GE.2.d0).AND.(DIV.LT.4.d0) ) marker='*  '
              IF( (DIV.GE.4.d0).AND.(DIV.LT.8.d0) ) marker='** '
              IF( (DIV.GE.8.d0) ) marker='***'
              IF(IEP(IBB).GT.0) WRITE(8,602) VP(IBB),JP(I),NEF(EFP(I)),
     1         VPP(IBB),JPP(I),NEF(EFPP(I)),FREQ(I),UFREQ(I),DFREQ(I),
     2                                      DFREQ(I)/UFREQ(I),marker
              IF(IEP(IBB).EQ.0) WRITE(8,602) VP(IBB),VPP(IBB),
     1                  NEF(EFP(I)),JP(I),JPP(I),NEF(EFPP(I)),FREQ(I),
     2                      UFREQ(I),DFREQ(I),DFREQ(I)/UFREQ(I),marker
            ELSE
              WRITE(8,602) VP(IBB),JP(I),NEF(EFP(I)),VPP(IBB),JPP(I),
     1                                           NEF(EFPP(I)),DFREQ(I)
              WRITE(4,608) JP(I),EFP(I),JPP(I),EFPP(I),DFREQ(I),UFREQ(I)
c* Print predictions in alternate (Lyon) format
c             WRITE(11,606)VP(IBB),VPP(IBB),JPP(I),JP(I)-JPP(I),DFREQ(I)
c 606 FORMAT(2I4,I5,I4,f13.4)
            ENDIF
          ENDDO
      WRITE(8,604)
      RETURN                       
  600 FORMAT(1x,59('='),'  calc-obs'/ "   v'  J' p'",
     1  '  v"  J" p"    FREQ(obs)     u(FREQ)    calc-obs  /u(FREQ)'/
     2  1x,69('-'))
  601 FORMAT(1x,36('=')/ "   v'  J' p'",'  v"  J" p"   FREQ(calc)'/
     1  1x,36('-'))
  602 FORMAT(2(2I4,A3),f14.6,2f12.6,f10.4,1x,A3)
  604 FORMAT(1x,69('-'))
  608 FORMAT(I5,I3,I5,I3,F13.4,F9.4)
      END   
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE PREPOT(LNPT,IAN1,IAN2,IMN1,IMN2,NPP,RR,RM2,VLIM,
     1                                                        VV,NCN)
c** Driver subroutine of package to read parameters and/or generate
c  values of a potential V(I) at the NPP input distances RR(I).
c====================== Version of 20 July 2001 ====================
c**** Subroutine Input:
c----------------------
c  LNPT  is an integer specifying the operational mode:
c      *  LNPT > 0  : for a new case for which all potential-defining
c                     parameters are read in & a description printed
c      *  LNPT.le.0 : if potential points are to be generated in exactly
c                     the same manner as on preceding call, but at
c                     different distances RR(I) (no reads or writes)
c  IAN1 & IAN2 are the atomic numbers and IMN1 & IMN2 the mass numbers
c        of atoms #1 & 2, used (if needed) to specify isotope masses for
c        calculating adiabatic and/or non-adiabatic B-O-B correction fx.
c  NPP (integer) is the number of input distances  RR(i) (in Angstroms)
c        at which potential values  VV(i) (in cm-1) are to be generated
c  RR  (real array) is set of NPP distances where potential calculated
c  RM2 (real array) on input is the (centrifugal) array of  1/RR(i)**2
c----------------------
c**** Subroutine Output:
c----------------------
c  VLIM (cm-1) is the absolute energy at the potential asymptote
c  VV (real array) is the set of function values generated (in cm-1)
c  RM2 values returned may (if appropriate) be modified to include B-O-B
c      corrections to the (centrifugal) potential  1/RR(i)**2
c  NCN is an integer power defining the asymptotically-dominant 
c      inverse-power long-range potential tail:  CNN/R**NCN 
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+ Calls GENINT (which calls PLYINTRP, SPLINT & SPLINE) ,  or POTGEN ++
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Set maximum array dimension for the input function values to be
c  interpolated over & extrapolated beyong
      INTEGER NTPMX
      PARAMETER (NTPMX= 1600) 
      INTEGER I,J,IAN1,IAN2,IMN1,IMN2,INPTS,ILR,IR2,JWR,LNPT,LPPOT,LWR,
     1  NCN,NLIN,NPP,NROW,NTP,NUSE,NPRS,NPRF
      REAL*8  RFACT,EFACT,RH,RMIN,VLIM,VSHIFT,VV(NPP),RR(NPP),RM2(NPP),
     1  XI(NTPMX),YI(NTPMX),RWR(20),VWR(20),VWRB(3),D1V(3),D1VB(3),
     2  D2V(3),CNN
c
c** Save variables needed for 'subsequent' LNPT.le.0 calls
      SAVE ILR,IR2,LPPOT,NTP,NUSE
      SAVE CNN,VSHIFT,XI,YI
c
      DATA VWRB/3*0.D0/,D1VB/3*0.D0/
      LPPOT= 0
      NPRS= 1
      NPRF= NPP
c
      IF(LNPT.GT.0) THEN
c** If NTP > 0 :  define potential by interpolation over & extrapolation
c        beyond the NTP read-in turning points using subroutine GENINT.
c   If NTP.le.0 : generate a (fully analytic) potential in POTGEN.
c** If LPPOT > 0 : at every |LPPOT|-th point, print potential and 
c      derivatives-by-differences. ***  If  LPPOT < 0  write potential
c      at every |LPPOT|-th point to channel-8 in a compact format **
c** VLIM (cm-1) is the energy associated with the potential asymptote.
c-----------------------------------------------------------------------
          READ(5,*) NTP, LPPOT, VLIM
c-----------------------------------------------------------------------
          WRITE(6,600) VLIM
          IF(NTP.GT.0) THEN
c** For a pointwise potential (NTP > 0), now read points & parameters
c  controlling how the interpolation/extrapolation is to be done.
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** NTP (read above) is number of turning points (XI,YI) to be read in.
c** If NUSE > 0  interpolate with NUSE-point piecewise polynomials
c    (usually choose NUSE even, say, = 6, 8 or 10). ***  If(NUSE.LE.0)
c    interpolate with cubic spline instead of local polynomials.
c** If IR2 > 0 , interpolate over  YI*XI**2 ; otherwise on  YI  itself
c   This may help if interpolation has trouble on steep repulsive wall.
c** ILR specifies how to extrapolate beyond largest input distance XI(i)
c  If ILR < 0 , fit last 3 points to:  VLIM - A*exp(-b*(R-R0)**2)
c  If ILR = 0 , fit last 3 points to:  VLIM - A*R**p *exp(-b*R)
c  If ILR = 1 : fit last two points to:  VLIM - A/R**B .
c** If(ILR > 1) fit last turning points to:  VLIM - sum{of ILR
c  inverse-power terms beginning with  1/R**NCN}. *** If CNN.ne.0 ,
c  leading coefficient fixed at  CNN ; otherwise get it from points too.
c* Assume read-in CNN value has units:  [(cm-1)(Angstroms)**'NCN'].
c* If ILR = 2 or 3 , successive higher power terms differ by  1/R**2
c* If ILR > 3 : successive higher power terms differ by factor  1/R
c-----------------------------------------------------------------------
              READ(5,*) NUSE, IR2, ILR, NCN, CNN
c-----------------------------------------------------------------------
              IF(NTP.GT.NTPMX) THEN
                  WRITE(6,602) NTP,NTPMX
                  STOP
                  ENDIF
              IF(NUSE.GT.0) WRITE(6,604) NUSE,NTP
              IF(NUSE.LE.0) WRITE(6,606) NTP
              IF(IR2.GT.0) WRITE(6,608)
              IF((ILR.GT.1).AND.(DABS(CNN).GT.0.D0))WRITE(6,610)CNN,NCN
c** Read in turning points to be interpolated over
c** RFACT & EFACT are factors required to convert units of input turning
c       points (XI,YI) to Angstroms & cm-1, respectively (may = 1.d0)
c** Turning points (XI,YI) must be ordered with increasing XI(I)
c** Energy VSHIFT (cm-1) is added to the input potential points to
c   make their absolute energy consistent with VLIM (often VSHIFT=Te).
c-----------------------------------------------------------------------
              READ(5,*) RFACT, EFACT, VSHIFT
              READ(5,*) (XI(I), YI(I), I= 1,NTP)
c-----------------------------------------------------------------------
              WRITE(6,612) RFACT,EFACT
              NROW= (NTP+2)/3
              DO  J= 1,NROW
                  IF(EFACT.LE.10.D0) THEN
                      WRITE(6,614) (XI(I),YI(I),I= J,NTP,NROW)
                    ELSE
                      WRITE(6,616) (XI(I),YI(I),I= J,NTP,NROW)
                    ENDIF
                  ENDDO
              WRITE(6,618) VSHIFT
              DO  I= 1,NTP
                  YI(I)= YI(I)*EFACT+ VSHIFT
                  XI(I)= XI(I)*RFACT
                  ENDDO
              IF(IR2.GT.0) THEN
                  DO  I= 1,NTP
                      YI(I)= YI(I)*XI(I)**2
                      ENDDO
                  ENDIF
              ENDIF
          ENDIF
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IF(NTP.GT.0) THEN
          CALL GENINT(LNPT,NPP,RR,VV,NUSE,IR2,NTP,XI,YI,VLIM,ILR,
     1                                              NCN,CNN,NPRS,NPRF)
        ELSE
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Generate a fully analytic potential in subroutine POTGEN ***********
c* Potentials generated in cm-1 with potential asymptote at energy VLIM
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** IPOTL specifies the type of potential function to be generated.
c** MPAR & NPAR are integers for specifying potential types.
c** NVARB is number of (real*8) potential parameters read in.
c** IBOB specifies whether (if > 0) or not (if .le. 0) atomic mass 
c      dependent Born-Oppenheimer breakdown corrections will be included
c** For all functions considered, well depth and equilibrium distance
c  are read as  DSCM (cm-1)  and REQ (Angstroms), respectively.
c* [Most read-in parameters are dimensionless (scaled by DSCM & REQ).]
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c** If IPOTL=1  generate an L.J.(MPAR,NPAR) potential.
c** If IPOTL=2  generate an MLJ(NPAR) potential [JCP 112, 3949 (2000)]
c      If MPAR.ge.0  exponent parameter is polynomial of order (NVARB-1)
c           in z=(R-Re)/(R+Re), with the NVARB coefficients  PARM(j)
c      If MPAR < 0  exponent polynomial in  z  has order (NVARB-4) with 
c           coefficients PARM(i) (i= 1,NVARB-3), & includes a switching
c           function with exponent coefficient  ALPHA= PARM(NVARB)  and  
c           RSW= PARM(NVARB-1),  defined to yield limiting inverse-power
c           potential coefficient  Cn= PARM(NVARB-2).
c** If IPOTL=3  generate a Morse or Extended Morse Oscillator potential
c       with exponent factor "beta" defined as a power series of order
c       (NVARB-1) in  z=(R-Re)/(R+Re)  with NVARB coefficients PARM(i).
c       Set  NVARB= 1 for conventional "simple" Morse potential.
c*  Special option #1: set  MPAR= -1  to produce Wei Hua's 4-parameter 
c      modified Morse function with  b= PARM(1)  and C= PARM(2).
c*  Special option #2: set  MPAR= -2  to produce Coxon's "Generalized
c      Morse Oscillator" potential with exponent expansion in (R-Re)]
c ...  otherwise, set  MPAR.ge.0
c** If IPOTL=4  use Seto's modification of Surkus' GPEF expansion in
c       z = [R^NPAR - Re^NPAR]/[a*R^NPAR + b*Re^NPAR] where 
c       a=PARM(NVARB-1) & b=PARM(NVARB), which incorporates Dunham, SPF,
c       O-T and other forms: V(z) = c_0 z^2 [1 + c_1 z + c_2 z^2 + ...]
c       where  c_0 [cm-1] is read in as DSCM, and the first (NVARB-2)
c       PARM(i)'s are the  c_i  (i > 0).  [MPAR is dummy parameter here]
c  * For Dunham case:  NPAR=1, PARM(NVARB-1)= 0.0, PARM(NVARB)= 1.0
c  * For SPF case:  NPAR=1, PARM(NVARB-1)= 1.0, PARM(NVARB)= 0.0
c  * For Ogilvie-Tipping:  NPAR=1, PARM(NVARB-1)= 0.5 = PARM(NVARB)
c  * NOTE that for Surkus NPAR < 0 case:  z(NPAR,a,b)= z(|NPAR|,-b,-a)
c      Generate & return the  D_e  value implied by these coefficients.
c** If IPOTL=5  generate generalized HFD(NPAR,6,8,10,12,14) potential.
c       PARM(1-3) are the parameters defining the HFD damping function
c       D(x)=exp[-pparm(1)*(PARM(2)/x - 1)**PARM(3)] {for x < PARM(2)}
c       PARM(4) the quadratic coefficient in the exponent, and
c       PARM(5) is the power of  x=R/Req  multiplying the repulsive term
c              AREP*x**PARM(5) *exp[-beta*x - PARM(4)*x**2] ;
c       PARM(6-11)  are the reduced C_NPAR, C_6, C_8, C_10, C_12 and C14
c       parameters (NPAR < 6), while  AREP and  beta  are defined
c       by having the potential minimum at  x=1.  For NVARB < 11, higher
c       C_m coefficients automatically zero;  necessarily  NVARB.ge.7 .
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c** IBOB > 0, add atomic-mass-dependent Born-Openheimer breakdown 
c  correction functions to rotationless and/or centrifugal potential(s).
c  Both expressed as power series in  z= (R-Re)/(R+Re) starting with the
c  constant term, using the mass shift convention of Le Roy [J.Mol.Spec. 
c  194, 189 (1999)].  Adiabatic B-O-B potential correction fx. defined 
c  by polynomials of order NC1 with (NC1+1) coefficients {CA1(i)} for 
c  atom-1  and order NC2 with (NC2+1) coefficients {CA2(i)} for atom-2, 
c  while centrifugal correction fx. defined polynomial of order NG1 with
c  (NG1+1) coefficients {GA1(i)} for atom-1 and order NG2 with (NG2+1)
c  coefficients {GA2(i)} for atom-2.
c** Input parameters IANi & IMNi are the atomic & mass number of atom-i
c  (i=1,2), while integers RMN1 & RMN2 read here are the mass numbers of
c  the reference isotopes defining the B-O-B correction functions.
c** NC1 & NC2 are orders of polynomials DELTA(V,atom-i) defining 
c  'adiabatic' corrections to the rotationless potential for atoms 1 & 2
c  DELTA(V)= (1-M1ref/M1)*DELTA(V,atom-1) + (1-M2ref/M2)*DELTA(V,atom-2)
c** NG1 & NG2 are orders of polynomials q1(z) & q2(z) defining B-O-B
c   correction to the centrifugal potential:
c      V(centrifugal)= [1 + (M1ref/M1)*q1(z) + (M2ref/M2)*q2(z)]/R**2
c ... to omit a particular correction set associated NCi or NGi .lt.0
c** RX > 0.0  invokes Coxon's (older) expansions in (R-Re) for potential
c     correction and in  [(R-Rx)**j - (Re-Rx)**j] for centrifugal corrn.
c ... OTHERWISE (to use Le Roy B-O-B formalism) set  RX.le.0.d0 !!
c-----------------------------------------------------------------------
c** Read inside subroutine POTGEN
c         IF(LNPT.GT.0) THEN
c             READ(5,*) IPOTL, MPAR, NPAR, NVARB, IBOB, DSCM, REQ
c             IF(NVARB.GT.0) READ(5,*) (PARM(I), I=1,NVARB)
c             IF(IBOB.GT.0) THEN
c                 READ(5,*) RMN1, RMN2, NC1, NC2, NG1, NG2, RX
c                 IF(NC1.GE.0) READ(5,*) (CA1(I), I=0,NC1)
c                 IF(NC2.GE.0) READ(5,*) (CA2(I), I=0,NC2)
c                 IF(NG1.GE.0) READ(5,*) (GA1(I), I=0,NG1)
c                 IF(NG2.GE.0) READ(5,*) (GA2(I), I=0,NG2)
c             ENDIF
c         ENDIF
c-----------------------------------------------------------------------
          NCN= 99
          CALL POTGEN(LNPT,NPP,IAN1,IAN2,IMN1,IMN2,VLIM,RR,RM2,VV,
     1                                                        NCN,CNN)
        ENDIF
      IF(LPPOT.NE.0) THEN
c** If desired, on the first pass (i.e. if LNPT > 0) print the potential
          RH= RR(2)-RR(1)
          INPTS= IABS(LPPOT)
          IF(LPPOT.LT.0) THEN
c** Option to write resulting function compactly to channel-8. 
              RMIN= RR(1)
              NLIN= NPP/INPTS+ 1
              WRITE(8,800) NLIN,VLIM
              WRITE(8,802) (RR(I),VV(I),I= 1,NPP,INPTS)
            ELSE
c** Option to print potential & its 1-st three derivatives, the latter
c  calculated by differences, assuming equally spaced RR(I) values.
              WRITE(6,620)
              NPRS= MAX(1,(NPRS- 15*INPTS))
              NPRF= MIN(NPP,(NPRF+ 15*INPTS))
              NLIN= (NPRF-NPRS+1)/(2*INPTS)+1
              RH= INPTS*RH
              DO  I= 1,NLIN
                  LWR= NPRS+INPTS*(I-1)
                  DO  J= 1,2
                      JWR= LWR+(J-1)*NLIN*INPTS
                      IF(JWR.LE.NPP) THEN
                          RWR(J)= RR(JWR)
                          VWR(J)= VV(JWR)
                          D1V(J)= (VWR(J)-VWRB(J))/RH
                          VWRB(J)= VWR(J)
                          D2V(J)= (D1V(J)-D1VB(J))/RH
                          D1VB(J)= D1V(J)
                        ELSE
                          RWR(J)= 0.d0
                          VWR(J)= 0.d0
                        ENDIF
                      IF(I.LE.2) THEN
                          D2V(J)= 0.d0
                          IF(I.EQ.1) D1V(J)= 0.d0
                          ENDIF
                      ENDDO
                  WRITE(6,622) (RWR(J),VWR(J),D1V(J),D2V(J),J= 1,2)
                  ENDDO
            ENDIF
          ENDIF
      IF(LNPT.GT.0) WRITE(6,624)
      RETURN
  600 FORMAT(' Absolute energy at asymptote:   Y(lim)=',F12.4,'(cm-1)')
  602 FORMAT(/' **** ERROR in dimensioning of arrays required'
     1 ,' by GENINT;   No. input points ',I5,' > NTPMX =',I4)
  604 FORMAT(' Perform',I3,'-point piecewise polynomial interpolation ov
     1er',I5,' input points' )
  606 FORMAT(' Perform cubic spline interpolation over the',I5,' input p
     1oints' )
  608 FORMAT(' Interpolation performed over modified input array:   Y(I)
     1 * R(I)**2')
  610 FORMAT(/' Beyond read-in points extrapolate to limiting asymptotic
     1 behaviour:'/20x,'Y(R)  =  Y(lim) - (',D16.7,')/R**',I2)
  612 FORMAT(/' Scale input points:  (distance)*',1PD16.9,'   &  (energy
     1)*',D16.9/13x,'to get required internal units  [Angstroms & cm-1 f
     2or potentials]'/3('      R(i)         Y(i)  ')/3(3X,11('--')))
  614 FORMAT((3(F13.8,F12.4)))
  616 FORMAT((3(F12.6,F13.8)))
  618 FORMAT(3x,24('---')/' To make above input points consistent with',
     1  '  Y(lim),  add  Y(shift)=',F12.4)
  620 FORMAT(/'  Function and first 2 derivatives by differences'/
     1  2('     R       Y(R)     d1Y/dR1    d2Y/dR2')/2(2X,19('--')))
  622 FORMAT(2(0PF8.3,F11.3,1PD11.3,D10.2))
  624 FORMAT(1x,39('--'))
  800 FORMAT(I7,' function values with asymptotic value:',F14.6)
  802 FORMAT((1X,3(F12.8,F14.6)))
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE GENINT(LNPT,NPP,XX,YY,NUSE,IR2,NTP,XI,YI,VLIM,ILR,
     1                                              NCN,CNN,NPRS,NPRF)
c** GENINT produces a smooth function YY(i) at the NPP input distances
c  XX(i) by performing numerical interpolation over the range of the 
c  NTP input function values YI(j) at the distances XI(j), and using
c  analytic functions to extrapolate beyond their range to with an
c  exponential at short range and a form specified by ILR, NCN & CNN
c** ILR specifies how to extrapolate beyond largest given turning pts
c   If ILR < 0 , fit last 3 points to:  VLIM - A*exp(-b*(R-R0)**2)
c   If ILR = 0 , fit last 3 points to:  VLIM - A*R**p *exp(-b*R)
c   If ILR = 1 : fit last two points to:  VLIM - A/R**B .
c* If(ILR.ge.2) fit last turning points to:  VLIM - sum(of ILR
c  inverse-power terms beginning with  1/R**NCN). *** If CNN.ne.0 ,
c  leading coefficient fixed at  CNN ; otherwise get it from points too.
c* Assume read-in CNN value has units:  ((cm-1)(Angstroms)**'NCN').
c  If ILR = 2 or 3 , successive higher power terms differ by  1/R**2
c  If ILR > 3 : this factor is  1/R .
c=== Calls subroutines PLYINTRP, SPLINT & SPLINE ==================
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      INTEGER  I,J,IFXCN,IDER,IR2,ILR,ISR,LNPT,MBEG,MFIN,MINNER,
     1  NN,NPP,NUSE,NUST,NORD,NCN,NCN2,NCN4,NPRS,NPRF,NTP,
     2  IMX1,NMX,JR2,JMAX,MI(10),MF(10)
      REAL*8  ASR,BSR,CSR,ALR,BLR,CLR,DCSR,ADCSR,PDCSR,VRAT,
     1  DX1,DX2,DX3,EX1,EX2,EX3,CNN,VLIM,X1,X2,X3,Y1,Y2,Y3,
     1  XX(NPP),YY(NPP),XI(NTP),YI(NTP),XJ(20),YJ(20),DUMM(20)
c
      SAVE ASR,BSR,CSR,ISR,ALR,BLR,CLR,IMX1,NMX,JR2,JMAX
c
      NUST= NUSE/2
      IF(NUSE.LE.0) NUST= 2
      IDER= 0
      NPRS= 1
      NPRF= NPP
c** Determine if/where need to begin extrapolation beyond input data
c  XX(MI(J))  is the 1-st mesh point past turning point  XI(J) .
c  XX(MF(J))  is the last mesh point before turning point  XI(NTP+1-J)
      DO 6 J = 1,NUST
          MI(J)= 1
          MF(J)= 0
          DO  I= 1,NPP
              IF(XX(I).LE.XI(J)) MI(J)= I+ 1
              IF(XX(I).GE.XI(NTP+1-J)) GO TO 6
              MF(J)= I
              ENDDO
    6     CONTINUE
      IF(NUST.LT.2) THEN
          MFIN= MI(1)-1
        ELSE
          MFIN= MI(2)-1
        ENDIF
      IF(LNPT.GT.0) THEN
c-----------------------------------------------------------------------
c** For a new case determine analytic functions for extrapolating beyond
c  the range of the input points (if necessary) on this or later calls.
c** Try to fit three innermost turning points to  V(R)=A+B*DEXP(-C*R).
c** If unsatisfactory, extrapolate inward with inverse power function
          IF(IR2.LE.0) THEN
              DO  I= 1,4
                  YJ(I)= YI(I)
                  ENDDO
            ELSE
              DO  I= 1,4
                  YJ(I)= YI(I)/XI(I)**2
                  ENDDO
            ENDIF
          X1= XI(1)
          X2= XI(2)
          X3= XI(3)
          Y1= YJ(1)
          Y2= YJ(2)
          Y3= YJ(3)
          IF((Y1-Y2)*(Y2-Y3).LE.0.d0) THEN
c** If 3 innermost points not monotonic, use A+B/X inward extrapoln.
              ISR= 0
              WRITE(6,600)
            ELSE
c** Use cubic through innermost points to get initial trial exponent
c  from ratio of derivatives,  Y''/Y'
              IDER= 2
              ISR= 4
              CALL PLYINTRP(XI,YJ,ISR,X2,XJ,ISR,IDER)
              CSR= XJ(3)/XJ(2)
              DCSR= DABS(CSR*X2)
              IF(DCSR.GT.1.5D+2) THEN
c** If exponential causes overflows, use inverse power inward extrapoln.
                  ISR= 0
                  WRITE(6,602) CSR
                  GO TO 20
                  ENDIF
c** Prepare parameters for inward exponential extrapolation
              VRAT= (Y3- Y2)/(Y1- Y2)
              DX1= X1- X2
              DX3= X3- X2
              EX2= 1.D0
              ADCSR= 1.d99
c** Now iterate (with actual point) to get exact exponent coefficient 
              DO  J= 1,15
                  PDCSR= ADCSR
                  EX1= DEXP( CSR*DX1)
                  EX3= DEXP( CSR*DX3)
                  DCSR= (VRAT- (EX3- EX2)/(EX1- EX2)) /
     1   ((X3*EX3- X2 - (X1*EX1- X2)*(EX3-EX2)/(EX1- EX2))/(EX1- EX2))
                  ADCSR= ABS(DCSR)
                  IF((ADCSR.GT.PDCSR).AND.(ADCSR.LT.1.d-8)) GO TO 12
                  IF(ADCSR.LT.1.d-12) GO TO 12
                  CSR= CSR+ DCSR 
                  ENDDO
              WRITE(6,604) DCSR
   12         BSR= (Y1-Y2)/(EX1-EX2)
              ASR= Y2-BSR*EX2
              BSR= BSR*DEXP(-CSR*X2)
              WRITE(6,606) X2,ASR,BSR,CSR
            ENDIF
   20     IF(ISR.LE.0) THEN
              IF((X1*X2).LE.0.d0) THEN
c** If 1'st two mesh points of opposite sign, extrapolate linearly
                  ISR= -1
                  ASR= Y2
                  BSR= (Y2- Y1)/(X2- X1)
                  CSR= X2
                  WRITE(6,608) X2,ASR,BSR,CSR
                ELSE
c** For inward extrapolation as inverse power through 1'st two points ..
                  BSR= (Y1-Y2)* X1*X2/(X2- X1)
                  ASR= Y1-BSR/X1
                  CSR= X2
                  WRITE(6,610) X2,ASR,BSR
                ENDIF
              ENDIF
          ENDIF
  600 FORMAT('  ** CAUTION ** Exponential inward extrapolation fails sin
     1ce first 3 points not monotonic, ... so ...')
  602 FORMAT(' *** CAUTION ** inward extrapolation exponent coefficient
     1   C=',D12.4/10x,'could cause overflows, ... so ...')
  604 FORMAT(' *** CAUTION ** after 15 tries inward extrap. exponent coe
     1fft change is',1PD9.1)
  606 FORMAT(' Extrapolate to   X .le.',F7.4,'  with'/'   Y=',F13.3,
     1  SP,1PD15.6,' * exp(',SS,D13.6,'*R)')
  608 FORMAT(' Extrapolate to   X .le.',F8.4,'   with'/'   Y=',F13.3,
     1  SP,1PD16.7,' * [R - (',SS,F8.4,')]')
  610 FORMAT(' Extrapolate to   X .le.',F8.4,'   with   Y=',F12.3,
     1  SP,1PD15.6,')/R**1')
c
      IF(MFIN.GT.0) THEN
c** If needed, calculate function in inner extrapolation region
          IF(ISR.GT.0) THEN
c ... either as an exponential
              DO  I= 1,MFIN
                  EX1= CSR*XX(I)
                  IF(DABS(EX1).GT.1.D+2) EX1= 1.D+2*DSIGN(1.d0,EX1)
                  YY(I)= ASR+BSR*DEXP(EX1)
                  ENDDO
            ELSEIF(ISR.EQ.0) THEN
c ... or if that fails, as an inverse power
              DO  I= 1,MFIN
                  YY(I)= ASR+BSR/XX(I)
                  ENDDO
            ELSEIF(ISR.LT.0) THEN
c ... or if X changes sign, extrapolate inward linearly
              DO  I= 1,MFIN
                  YY(I)= ASR+ BSR*(XX(I)- CSR)
                  ENDDO
            ENDIF
          ENDIF
c** End of inward extrapolation procedure
c-----------------------------------------------------------------------
      MINNER= MFIN
      IF(NUST.GT.2) THEN
c** If(NUSE.gt.5) minimize spurious behaviour by interpolating with
c  order less than NUSE on intervals near inner end of range
          DO  J= 3,NUST
              NORD= 2*(J-1)
              MBEG= MI(J-1)
              MFIN= MI(J)-1
              IF(MFIN.GE.MBEG) THEN
                  DO  I=  MBEG,MFIN
                      CALL PLYINTRP(XI,YI,NTP,XX(I),DUMM,NORD,IDER)
                      YY(I)= DUMM(1)
                      ENDDO
                  ENDIF
              ENDDO
          ENDIF
c** Main interpolation step begins here
c=======================================================================
      MBEG= MI(NUST)
      MFIN= MF(NUST)
      IF(MFIN.GE.MBEG) THEN
          IF(NUSE.LE.0) THEN
c** Either ... use cubic spline for main interpolation step
              CALL SPLINT(LNPT,NTP,XI,YI,MBEG,MFIN,XX,YY)
            ELSE
c ... or use piecewise polynomials for main interpolation step
              DO  I= MBEG,MFIN
                  CALL PLYINTRP(XI,YI,NTP,XX(I),DUMM,NUSE,IDER)
                  YY(I)= DUMM(1)
                  ENDDO
            ENDIF
          ENDIF
      IF(MFIN.LT.NPP) THEN
          IF(NUST.LE.2) THEN
c** If(NUSE.gt.5) minimize spurious behaviour by interpolating with
c  order less than NUSE on intervals near outer end of range
              MBEG= MF(NUST)+1
            ELSE
              NN= NUST-2
              DO  J= 1,NN
                  NORD= 2*(NUST-J)
                  MBEG= MF(NUST-J+1)+1
                  MFIN= MF(NUST-J)
                  IF(MFIN.GE.MBEG) THEN
                      DO  I= MBEG,MFIN
                          CALL PLYINTRP(XI,YI,NTP,XX(I),DUMM,NORD,IDER)
                          YY(I)= DUMM(1)
                          ENDDO
                      END IF
                  ENDDO
            ENDIF
          ENDIF
      MBEG= MFIN+1
      IF((MFIN.GT.MINNER).AND.(IR2.GT.0)) THEN
c** In (IR2.gt.0) option, now remove X**2 from the interpolated function
          DO  I= MINNER+1,MFIN
              YY(I)= YY(I)/XX(I)**2
              ENDDO
          ENDIF
c** Print test of smoothness at join with analytic inward extrapolation
c     IF(LNPT.GT.0) THEN
c         MST= MAX0(MINNER-4,1)
c         NPRS= MST
c         MFN= MST+8
c         IF(MFN.GT.NPP) MFN= NPP
c         IF(MFN.GT.MFIN) MFN= MFIN
c         IF(MINNER.GT.0) WRITE(6,611) X2,((XX(I),YY(I),I= J,MFN,3),
c    1        J= MST,MST+2)
c 611 FORMAT('     Verify smoothness of inner join at   X=',F9.5/
c    1  (3X,3(F10.5,G15.7)))
c         ENDIF
c-----------------------------------------------------------------------
c** To extrapolate potential beyond range of given turning points ...
      IF(LNPT.GT.0) THEN
c** On first entry, calculate things needed for extrapolation constants
          Y1= YI(NTP)
          Y2= YI(NTP-1)
          Y3= YI(NTP-2)
          X1= XI(NTP)
          X2= XI(NTP-1)
          X3= XI(NTP-2)
          IF(IR2.GT.0) THEN
              Y1= Y1/X1**2
              Y2= Y2/X2**2
              Y3= Y3/X3**2
              ENDIF
          ENDIF
c** Check inverse-power tail power ...
      IF(NCN.LE.0) NCN= 6
      IF(ILR.LT.0) THEN
          IF(LNPT.GT.0) THEN
C** For  ILR.lt.0  use  Y = VLIM - ALR * exp[-CLR*(X - BLR)**2]
              EX1= DLOG((VLIM-Y1)/(VLIM-Y2))/(X1-X2)
              EX2= DLOG((VLIM-Y2)/(VLIM-Y3))/(X2-X3)
              BLR= (X1+X2 - (X2+X3)*EX1/EX2)/(2.d0- 2.d0*EX1/EX2)
              CLR= -EX1/(X1+X2-2.d0*BLR)
              ALR= (VLIM-Y1)*DEXP(CLR*(X1-BLR)**2)
              WRITE(6,614) X2,VLIM,ALR,CLR,BLR
              IF(CLR.LT.0.d0) THEN
c ... but replace it by an inverse power of exponent constant negative
                  WRITE(6,612)
                  ILR= 1
                  GO TO 50
                  ENDIF
              ENDIF
          IF(MBEG.LE.NPP) THEN
              DO  I= MBEG,NPP
                  YY(I)= VLIM- ALR*DEXP(-CLR*(XX(I) - BLR)**2)
                  ENDDO
              ENDIF
          GO TO 90
          ENDIF
      IF(ILR.EQ.0) THEN
c** For ILR.le.0  use  Y = VLIM - ALR * X**p * exp(-CLR*X)
          IF(LNPT.GT.0) THEN
              EX1= DLOG((VLIM-Y1)/(VLIM-Y2))/(X1-X2)
              EX2= DLOG((VLIM-Y2)/(VLIM-Y3))/(X2-X3)
              DX1= DLOG(X1/X2)/(X1-X2)
              DX2= DLOG(X2/X3)/(X2-X3)
              BLR= (EX1-EX2)/(DX1-DX2)
              CLR= BLR*DX1- EX1
              ALR= (VLIM-Y1)* DEXP(CLR*X1)/X1**BLR 
              WRITE(6,616) X2,VLIM,ALR,BLR,CLR
              IF(CLR.LT.0.d0) THEN
c ... but replace it by an inverse power of exponent constant negative
                  WRITE(6,612)
                  ILR= 1
                  GO TO 50
                  ENDIF
              ENDIF
          IF(MBEG.LE.NPP) THEN
              DO  I= MBEG,NPP
                  YY(I)= VLIM- ALR*XX(I)**BLR *DEXP(-CLR*XX(I))
                  ENDDO
              ENDIF
          GO TO 90
          ENDIF
   50 IF(ILR.EQ.1) THEN
c** For  ILR=1 ,  use     Y = VLIM + ALR/X**BLR
          IF(LNPT.GT.0) THEN
              BLR= DLOG((VLIM-Y2)/(VLIM-Y1))/DLOG(X1/X2)
              ALR= (Y1- VLIM)*X1**BLR
              NCN= BLR
              WRITE(6,618) X2,VLIM,ALR,BLR,NCN
              ENDIF
          IF(MBEG.LE.NPP) THEN
              DO  I= MBEG,NPP
                  YY(I)= VLIM+ ALR/XX(I)**BLR
                  ENDDO
              ENDIF
          GO TO 90
          ENDIF
c** Set constants for long-range extrapolation
      IFXCN= 0
      IF((CNN.GT.0.d0).OR.(CNN.LT.0.d0)) IFXCN= 1
      NCN2= NCN+2
      IF(ILR.EQ.2) THEN
c** For ILR=2 ,  use   Y = VLIM - CNN/X**NCN - BSR/X**(NCN+2)
c*  If CNN held fixed need ILR > 2  to prevent discontinuity
          IF(LNPT.GT.0) THEN
              IF(IFXCN.LE.0) THEN
                  CNN= ((VLIM-Y1)*X1**NCN2 -
     1                 (VLIM-Y2)*X2**NCN2)/(X1**2-X2**2)
                  ENDIF
              ALR= CNN
              BLR= (VLIM-Y1)*X1**NCN2 - CNN*X1**2
              WRITE(6,620) X2,VLIM,CNN,NCN,BLR,NCN2
              ENDIF
          IF(MBEG.LE.NPP) THEN
              DO  I= MBEG,NPP
                  YY(I)= VLIM-(ALR+BLR/XX(I)**2)/XX(I)**NCN
                  ENDDO
              ENDIF
          GO TO 90
          ENDIF
      IF(ILR.EQ.3) THEN
c** For ILR=3 , use   Y = VLIM - (CN + CN2/X**2 + CN4/X**4)/X**NCN
          IF(LNPT.GT.0) THEN
              NCN4= NCN+4
              IF(IFXCN.GT.0) THEN
                  ALR= CNN
                  BLR= (((VLIM-Y1)*X1**NCN-ALR)*X1**4-((VLIM-Y2)
     1                     *X2**NCN-ALR)*X2**4)/(X1**2-X2**2)
                  CLR= ((VLIM-Y1)*X1**NCN-ALR-BLR/X1**2)*X1**4
                ELSE
                  EX1= X1**2
                  EX2= X2**2
                  EX3= X3**2
                  DX1= (VLIM-Y1)*X1**NCN4
                  DX2= (VLIM-Y2)*X2**NCN4
                  DX3= (VLIM-Y3)*X3**NCN4
                  BLR= (DX1-DX2)/(EX1-EX2)
                  ALR= (BLR-(DX2-DX3)/(EX2-EX3))/(EX1-EX3)
                  BLR= BLR-ALR*(EX1+EX2)
                  CLR= DX1-(ALR*EX1+BLR)*EX1
                ENDIF
              WRITE(6,622) X2,VLIM,ALR,NCN,BLR,NCN2,CLR,NCN4
              ENDIF
          IF(MBEG.LE.NPP) THEN
              DO  I= MBEG,NPP
                  EX2= 1.d0/XX(I)**2
                  YY(I)= VLIM-(ALR+EX2*(BLR+EX2*CLR))/XX(I)**NCN
                  ENDDO
              ENDIF
          GO TO 90
          ENDIF
      IF(ILR.GE.4) THEN
c** For ILR.ge.4,   Y = VLIM-SUM(BB(K)/X**K) , (K=NCN,NMX=NCN+ILR-1)
          IF(LNPT.GT.0) THEN
              IF(NCN.LE.0) NCN= 1
              IMX1= ILR-1
              NMX= NCN+IMX1
              JR2= 0
              IF(IR2.GT.0) JR2= 2
              IDER= 0
              JMAX= ILR
              IF(IFXCN.GT.0) JMAX= IMX1
              WRITE(6,624) X2,ILR,NCN,VLIM
              IF(IFXCN.GT.0) WRITE(6,626) NCN,CNN
              ENDIF
c** Actually extrapolate with polynomial fitted to the last JMAX
c  values of  (VLIM - YI(I))*XI(I)**NMX  , & then convert back to  YY(I).
          IF(MBEG.LE.NPP) THEN
              J= NTP- JMAX
              DO  I= 1,JMAX
                  J= J+1
                  XJ(I)= XI(J)
                  YJ(I)= (VLIM-YI(J)/XI(J)**JR2)*XI(J)**NMX
                  IF(IFXCN.GT.0) YJ(I)= YJ(I)- CNN*XI(J)**IMX1
                  ENDDO
              DO  I= MBEG,NPP
                  CALL PLYINTRP(XJ,YJ,JMAX,XX(I),DUMM,JMAX,IDER)
                  YY(I)= DUMM(1)
                  IF(IFXCN.GT.0) YY(I)= YY(I)+ CNN*XX(I)**IMX1
                  YY(I)= VLIM-YY(I)/XX(I)**NMX
                  ENDDO
              ENDIF
          ENDIF
c** Finished extrapolation section.
   90 CONTINUE
c** Test smoothness at outer join to analytic extrapolation function
c     IF((LNPT.GT.0).AND.(MBEG.LE.NPP)) THEN
c         MST= MBEG-5
c         IF(MST.LT.1) MST= 1
c         MFN= MST+8
c         IF(MFN.GT.NPP) MFN= NPP
c         WRITE(6,627) X2,((XX(I),YY(I),I= J,MFN,3),J= MST,MST+2)
c         NPRF= MFN
c         ENDIF
c 627 FORMAT('     Verify smoothness of outer join at   X=',F9.5/
c    1  (3X,3(F10.5,G15.7)))
      RETURN
  612 FORMAT('  *** BUT *** since exponent has positive coefficient, swi
     1tch form ...')
  614 FORMAT(' Function for  X .GE.',F8.4,'   generated as'/'   Y=',
     1  F12.4,' - (',1PD13.6,') * exp{-',0PF10.6,' * (R -',F9.6,')**2}')
  616 FORMAT(' Function for  X .GE.',F8.4,'   generated as'/'   Y=',
     1 F12.4,' - (',1PD13.6,') * R**',0PF10.6,'  * exp{-(',F11.6,'*R)}')
  618 FORMAT(' Extrapolate to  X .GE.',F8.4,'  using'/'   Y=',
     1  F12.4,SP,1PD15.6,'/X**(',SS,D13.6,')] ,  yielding   NCN=',I3)
  620 FORMAT(' Extrapolate to  X .GE.',F8.4,'  using'/'   Y=',
     1  F12.4,' - [',1PD13.6,'/X**',I1,SP,D14.6,'/X**',SS,I1,']')
  622 FORMAT(' Extrapolate to  X .GE.',F8.4,'  using'/
     1  '   Y=',F12.4,' - [',1PD13.6,'/X**',I1,SP,D14.6,'/X**',
     2  SS,I1,SP,D14.6,'/X**',SS,I2,']')
  624 FORMAT(' Function for  X .GE.',F7.3,'  generated by',I3,
     1 '-point inverse-power interpolation'/'   with leading term  1/R**
     2',I1,'  relative to dissociation limit   YLIM=',F11.3)
  626 FORMAT('   and (dimensionless) leading coefficient fixed as   C',
     1  I1,'=',G15.8)
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE PLYINTRP(XI,YI,NPT,RR,C,NCFT,IDER)
c* From the NPT known mesh points (XI,YI) ,given in order of increasing
c  or decreasing XI(I), select the NCFT points (XJ,YJ) surrounding the 
c  given point RR, and by fitting an (NCFT-1)-th degree polynomial through
c  them, interpolate to find the function CC(1) and its first IDER 
c  derivatives (CC(I+1),I=1,IDER) evaluated at RR.
c* Adapted by  R.J. Le Roy  from algorithm #416,Comm.A.C.M.;  27/02/1988
c=======================================================================
      INTEGER  I,J,K,I1,I2,IFC,IM,IDER,J1,NH,NPT,NCFT
      REAL*8  RR,XX,XI(NPT),YI(NPT),C(NCFT),XJ(20),YJ(20)
c
      IF((NCFT.GT.20).OR.(NCFT.GT.NPT)) GO TO 101
      NH= NCFT/2
c** First locate the known mesh points (XJ,YJ) bracketing RR
      I1= 1
      I2= NCFT
      IF(NCFT.NE.NPT) THEN
          IF(XI(NPT).LE.XI(1)) THEN
              DO  I= 1,NPT
                  IM= I
                  IF(XI(I).LT.RR) GO TO 20
                  ENDDO
            ELSE
              DO  I= 1,NPT
                  IM= I
                  IF(XI(I).GT.RR) GO TO 20
                  ENDDO
            ENDIF
   20     I1= IM-NH
          IF(I1.LE.0) I1= 1
          I2= I1+NCFT-1
          IF(I2.GT.NPT) THEN
              I2= NPT
              I1= I2-NCFT+1
              ENDIF
          ENDIF
      J= 0
      DO  I= I1,I2
          J= J+1
          XJ(J)= XI(I)-RR
          YJ(J)= YI(I)
          ENDDO
c** Now determine polynomial coefficients C(I).
      DO  I= 2,NCFT
          I1= I-1
          K= I1+1
          DO  J= 1,I1
              K= K-1
              YJ(K)= (YJ(K+1)-YJ(K))/(XJ(I)-XJ(K))
              ENDDO
          ENDDO
      C(1)= YJ(1)
      DO  I= 2,NCFT
          XX= XJ(I)
          C(I)= C(I-1)
          IF(I.NE.2) THEN
              I1= I-1
              K= I1+1
              DO  J= 2,I1
                  K= K-1
                  C(K)= -XX*C(K)+C(K-1)
                  ENDDO
              ENDIF
          C(1)= YJ(I)-XX*C(1)
          ENDDO
c** Finally, convert polynomial coefficients to derivatives at RR.
      IFC= 1
      IF(IDER.GE.NCFT) IDER= NCFT-1
      IF(IDER.LE.1) GO TO 99
      DO  I= 2,IDER
          J= I+1
          IFC= IFC*I
          C(J)= C(J)*IFC
          ENDDO
      IF(J.LT.NCFT) THEN
          J1= J+1
          DO  I= J1,NCFT
              C(I)= 0.D+0
              ENDDO
          ENDIF
   99 RETURN
  101 WRITE(6,601) NCFT,NCFT,NPT
      STOP
  601 FORMAT(/' *** Dimensioning ERROR in PLYINTRP :  either   (NCFT=',
     1  I2,' .GT. 20)   or   (NCFT=',I2,' .GT. NPT=',I3,')')
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c**********************************************************************
      SUBROUTINE SPLINT(LNPT,NTP,R1,V1,MBEG,MEND,XX,YY)
c** Subroutine to generate (if LNPT.ge.0) 4*NTP coefficients CSP(J)
c  of a cubic spline passing through the NTP points (R1(J),V1(J))
c  and to then calculate values of the resulting function YY(I) at the
c  entering abscissae values XX(I) for  I=MBEG to MEND.
c** If LNPT < 0 , generate function values at the given XX(I) using
c  the coefficients CSP(J) obtained and SAVEd on a preceding call.
c** Assumes both R1(J) & XX(I) are monotonic increasing.
c+++++ Calls only subroutine SPLINE +++++++++++++++++++++++++++++++++++
c======================================================================
      INTEGER MAXSP
      PARAMETER (MAXSP=6400)
      INTEGER  I,IER,I1ST,IDER,JK,K,KK,LNPT,N2,N3,NIPT,NTP,MBEG,MEND
      REAL*8 EPS,R2,RI,RRR,TTMP,R1(NTP),V1(NTP),CSP(MAXSP),
     1  YY(MEND),XX(MEND)
      SAVE CSP
c
      IF(4*NTP.GT.MAXSP) THEN
          WRITE(6,602) MAXSP,NTP
          STOP
          ENDIF
      EPS= 1.D-6*(R1(2)-R1(1))
      N2= 2*NTP
      N3= 3*NTP
      IF(LNPT.GT.0) THEN
c** On first pass for a given data set, generate spline function
c  coefficients in subroutine SPLINE
c** Start by using a cubic polynomial at each end of the range to get
c  the first derivative at each end for use in defining the spline.
          IDER= 1
          NIPT= 4
          I1ST= NTP-3
          CALL PLYINTRP(R1(I1ST),V1(I1ST),NIPT,R1(NTP),CSP,NIPT,IDER)
          TTMP= CSP(2)
          CALL PLYINTRP(R1,V1,NIPT,R1(1),CSP,NIPT,IDER)
          CSP(1)= CSP(2)
          CSP(2)= TTMP
c** Now call routine to actually generate spline coefficients
          CALL SPLINE(R1,V1,NTP,3,CSP,MAXSP,IER)
          IF(IER .NE. 0) THEN
              WRITE(6,604)
              STOP
              ENDIF
          ENDIF
      IF(MEND.LT.MBEG) GO TO 99
c** Now, use spline to generate function at desired points XX(I)
      DO  I= MBEG,MEND
          RI= XX(I)
          RRR= RI-EPS
          KK= 1
c** For a monotonic increasing distance array XX(I),  this statement 
c  speeds up the search for which set of cubic coefficients to use.
          IF(I.GT.MBEG) THEN
              IF(XX(I).GT.XX(I-1)) KK= JK
              ENDIF
          DO  K= KK,NTP
              JK= K
              IF(R1(K).GE.RRR) GO TO 64
              ENDDO
   64     CONTINUE
          JK= JK-1
          IF(JK.LT.1) JK= 1
          R2= RI-R1(JK)
          YY(I)= CSP(JK)+R2*(CSP(NTP+JK)+R2*(CSP(N2+JK)+R2*CSP(N3+JK)))
          ENDDO
   99 RETURN
  602 FORMAT(' *** ERROR in SPLINT ***  Array dimension  MAXSP=',I4,
     1  ' cannot contain spline coefficients for  NTP=',I4)
  604 FORMAT(' *** ERROR in generating spline coefficients in SPLINE')
      END
c**********************************************************************
      SUBROUTINE SPLINE(X,Y,N,IOPT,C,N4,IER)
c** Subroutine for generating cubic spline coefficients
c  C(J), (J=1,N4=4*N) through the N points X(I), Y(I).
c** C(I+M*N), M=0-3  are the coefficients of order  0-3  of cubic
c  polynomial expanded about X(I) so as to describe the interval:
c             -  X(I) to X(I+1)  , if  X(I)  in increasing order
c             -  X(I-1) to X(I)  , if  X(I)  in decreasing order.
c** IOPT indicates boundary conditions used in creating the  spline .
c*  If (IOPT=0)  second derivatives = zero at both ends of range.
c*  If (IOPT=1)  1st derivative at first point X(1) fixed at C(1),
c                and 2nd derivative at X(N) = zero.
c*  If (IOPT=2)  1st derivative at last point X(N) fixed at C(2),
c                and 2nd derivative at X(1) = zero.
c*  If (IOPT=3)  constrain first derivatives at end points to have
c                (read in) values  C(1)  at  X(1)  &  C(2)  at  X(N)
c** IER is the error flag.  IER=0  on return if routine successful.
c-----------------------------------------------------------------------
      INTEGER I,II,IER,IOH,IOL,IOPT,J,J1,J2,J3,NER,N,N4,JMP
      REAL*8  A,H,R,DY2,DYA,DYB,XB,XC,YA,YB, X(N),Y(N),C(N4)
c
      JMP= 1
      NER= 1000
      IF(N.LE.1) GO TO 250
c** Initialization
      XC= X(1)
      YB= Y(1)
      H= 0.D0
      A= 0.D0
      R= 0.D0
      DYB= 0.D0
      NER= 2000
c
c  IOL=0 - given derivative at firstpoint
c  IOH=0 - given derivative at last point
c
      IOL= IOPT-1
      IOH= IOPT-2
      IF(IOH.EQ.1) THEN
          IOL= 0
          IOH= 0
          ENDIF
      DY2= C(2)
c
c  Form the system of linear equations
c  and eliminate subsequentially
c
      J= 1
      DO  I= 1,N
          J2= N+I
          J3= J2+N
          A= H*(2.D0-A)
          DYA= DYB+H*R
          IF(I.GE.N) THEN
c
c  set derivative dy2 at last point
c
              DYB= DY2
              H= 0.D0
              IF(IOH.EQ.0) GOTO 200
              DYB= DYA
              GOTO 220
              ENDIF
          J= J+JMP
          XB= XC
          XC= X(J)
          H= XC-XB
c
c  II= 0 - increasing abscissae
c  II= 1 - decreasing abscissae
c
          II= 0
          IF(H.LT.0) II= 1
          IF(H.EQ.0) GO TO 250
          YA= YB
          YB= Y(J)
          DYB= (YB-YA)/H
          IF(I.LE.1) THEN
              J1= II
              IF(IOL.NE.0) GO TO 220
              DYA= C(1)
              ENDIF
200       IF(J1.NE.II) GO TO 250
          A= 1.D0/(H+H+A)
220       R= A*(DYB-DYA)
          C(J3)= R
          A= H*A
          C(J2)= A
          C(I)= DYB
          ENDDO
c
c  back substitution of the system of linear equations
c     and computation of the other coefficients
c
      A= 1.D0
      J1= J3+N+II-II*N
      I= N
      DO  IOL= 1,N
          XB= X(J)
          H= XC-XB
          XC= XB
          A= A+H
          YB= R
          R= C(J3)-R*C(J2)
          YA= R+R
          C(J3)= YA+R
          C(J2)= C(I)-H*(YA+YB)
          C(J1)= (YB-R)/A
          C(I)= Y(J)
          A= 0.D0
          J= J-JMP
          I= I-1
          J2= J2-1
          J3= J3-1
          J1= J3+N+II
          ENDDO
      IER= 0
      RETURN
250   IER= NER
      RETURN
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE POTGEN(LNPT,NPP,IAN1,IAN2,IMN1,IMN2,VLIM,XO,RM2,VV,
     1                                                        NCN,CNN)
c** Generate analytic potential  VV(i)  as specified by the choice
c  of parameter IPOTL (see comments in PREPOT (& in main program))
c** All potentials generated in units cm-1 with absolute asymptote at
c  (input) energy VLIM for distance array  X0(i) Angstroms.
c** Return with NCN equal to power of asymptotically dominant inverse
c  power term in long range part of potential
c** Born-Oppenheimer correction functions in IPOTL=3 option may have up
c  to NBOB+1 terms.
c-----------------------------------------------------------------------
      INTEGER NBOB
      PARAMETER (NBOB=20)
      INTEGER  I,J,M,IBOB,IAN1,IAN2,IMN1,IMN2,RMN1,RMN2,IORD,IPOTL,
     1  NC1,NC2,NG1,NG2,NCMAX,NPAR,MPAR,NVARB,NPP,LNPT,NCN,GNS,GEL
      CHARACTER*2 NAME1,NAME2
      REAL*8  A0,A1,A2,A3,ALFA,BETA,BINF,B1,B2,CSAV,
     1  ABUND,CNN,DSCM,DX,DX1,FCT,
     2  FC1,FC2,MASS1,MASS2,RMASS1,RMASS2,RC6,RC8,RC10,RC12,RC14,
     3  RCNPAR,RD,RDIF,REQ,RPOW,RX,SC1,SC2,SG1,SG2,VLIM,
     4  XDF,X1,XM,XN,XM2C,XP1,ZZ,ZP, CA1(0:NBOB),CA2(0:NBOB),
     5  GA1(0:NBOB),GA2(0:NBOB),PARM(20),XO(NPP),VV(NPP),RM2(NPP)
c
      SAVE IBOB,IPOTL,IORD,MPAR,NPAR,NVARB,DSCM,REQ,PARM,CA1,CA2,GA1,
     1     GA2,RX,CSAV
c
      IF(LNPT.GT.0) THEN
c** Parameter definitions listed preceeding CALL in subroutine PREPOT
c-----------------------------------------------------------------------
          READ(5,*) IPOTL, MPAR, NPAR, NVARB, IBOB, DSCM, REQ
          IF(IPOTL.EQ.1) NVARB= 0
          IF((IPOTL.EQ.3).AND.(NPAR.EQ.-1)) NVARB= 2
          IF(NVARB.GT.0)  READ(5,*) (PARM(I), I=1,NVARB)
          IF(IBOB.GT.0) THEN
              READ(5,*) RMN1, RMN2, NC1, NC2, NG1, NG2, RX
c-----------------------------------------------------------------------
              NCMAX= MAX0(NC1,NC2,NG1,NG2)
              IF(NCMAX.LT.0) THEN
                  IBOB= 0
                ELSE
c** If appropriate, read parameters & prepare to add mass-dep. BOB corrn
                  CALL MASSES(IAN1,IMN1,NAME1,GEL,GNS,MASS1,ABUND)
                  CALL MASSES(IAN1,RMN1,NAME1,GEL,GNS,RMASS1,ABUND)
                  CALL MASSES(IAN2,IMN2,NAME2,GEL,GNS,MASS2,ABUND)
                  CALL MASSES(IAN2,RMN2,NAME2,GEL,GNS,RMASS2,ABUND)
                  WRITE(6,628)
c  For simplicity, first zero out all correction function coefficients
                  DO  I=0,NCMAX
                      CA1(I)= 0.d0
                      CA2(I)= 0.d0
                      GA1(I)= 0.d0
                      GA2(I)= 0.d0
                      ENDDO
                  FC1= 0.d0
                  FC2= 0.d0
c=======================================================================
c** Read actual B-O-B polynomial expansion coefficients
c=======================================================================
                  IF(NC1.GE.0) THEN
c-----------------------------------------------------------------------
                      READ(5,*) (CA1(I), I=0,NC1)
c-----------------------------------------------------------------------
                      IF(RX.LE.0.d0) THEN
                          WRITE(6,630) 1,MASS1,NC1,NAME1,RMN1,NAME1,
     1                                    IMN1,NC1+1,(CA1(I),I= 0,NC1)
                          FC1= 1.d0 - RMASS1/MASS1
                        ELSE
                          WRITE(6,632) 1,MASS1,NC1,NAME1,IMN1,NC1+1,
     1                                               (CA1(I),I= 0,NC1)
                          FC1= 1.d0/MASS1
                        ENDIF
                      ENDIF
                  IF(NC2.GE.0) THEN
c-----------------------------------------------------------------------
                      READ(5,*) (CA2(I), I=0,NC2)
c-----------------------------------------------------------------------
                      IF(RX.LE.0.d0) THEN
                          WRITE(6,630) 2,MASS2,NC2,NAME2,RMN2,NAME2,
     1                                    IMN2,NC2+1,(CA2(I),I= 0,NC2)
                          FC2= 1.d0 - RMASS2/MASS2
                        ELSE
                          WRITE(6,632) 2,MASS2,NC2,NAME2,IMN2,NC2+1,
     1                                                (CA2(I),I=0,NC2)
                          FC2= 1.d0/MASS2
                        ENDIF
                      ENDIF
                  IF(NG1.GE.0) THEN
c-----------------------------------------------------------------------
                      READ(5,*) (GA1(I), I=0,NG1)
c-----------------------------------------------------------------------
                      IF(RX.LE.0.d0) THEN
                          WRITE(6,634) 1,MASS1,NG1,NAME1,RMN1,NAME1,
     1                                    IMN1,NG1+1,(GA1(I),I= 0,NG1)
                        ELSE
                          WRITE(6,636) 1,MASS1,NG1,NAME1,IMN1,RX,
     1                                         NG1+1,(GA1(I),I= 0,NG1)
                        ENDIF                     
                      ENDIF
                  IF(NG2.GE.0) THEN
c-----------------------------------------------------------------------
                      READ(5,*) (GA2(I), I=0,NG2)
c-----------------------------------------------------------------------
                      IF(RX.LE.0.d0) THEN
                          WRITE(6,634) 2,MASS2,NG2,NAME2,RMN2,NAME2,
     1                                    IMN2,NG2+1,(GA2(I),I= 0,NG2)
                        ELSE
                          WRITE(6,636) 2,MASS2,NG2,NAME2,IMN2,RX,
     1                                         NG2+1,(GA2(I),I= 0,NG2)
                        ENDIF
                      ENDIF
                  DO  I=0,NCMAX
                      CA1(I)= CA1(I)*FC1
                      CA2(I)= CA2(I)*FC2
                      IF(RX.LE.0.d0) THEN
                          GA1(I)= GA1(I)*(1.d0-FC1)
                          GA2(I)= GA2(I)*(1.d0-FC2)
                        ELSE
                          GA1(I)= GA1(I)*FC1
                          GA2(I)= GA2(I)*FC2
                        ENDIF
                      ENDDO
                ENDIF
              ENDIF
          ENDIF
      IF(IPOTL.EQ.1) THEN 
c=======================================================================
c** Generate a  Lennard-Jones(MPAR,NPAR)  potential here.
c=======================================================================
          XM= MPAR
          XN= NPAR
          XDF= DSCM/(XM-XN)
          IF(LNPT.GE.0) WRITE(6,600) MPAR,NPAR,DSCM,REQ
          NCN= NPAR
          CNN= XM*XDF*REQ**NPAR
          DO  I= 1,NPP
              VV(I)= (XN*(REQ/XO(I))**MPAR - XM*(REQ/XO(I))**NPAR)*XDF
     1                  +VLIM
              ENDDO
          ENDIF
      IF(IPOTL.EQ.2) THEN
c=======================================================================
c** Generate an MLJ potential [as per JCP 112, 3949 (2000)] here ...
c=======================================================================
          NCN= NPAR
          IORD= NVARB-1
          IF(MPAR.LT.0) THEN
c  If appropriate, prepare to calculate switching function
              IORD= IORD- 3
              CNN= PARM(NVARB-2)
              BINF= DLOG(2.d0*DSCM*REQ**NPAR/CNN)
              ALFA= PARM(NVARB- 1)
              RX= PARM(NVARB)
            ELSE
c  Generate limiting Cn value for non-switching case from BINF
              BINF= 0.D0
              DO  I= 1,NVARB
                  BINF= BINF+ PARM(I)
                  ENDDO
              CNN= 2.d0*DSCM*REQ**NPAR *DEXP(-BINF)
            ENDIF
          IF(LNPT.GT.0) THEN
              WRITE(6,602) DSCM,REQ,IORD,(PARM(J),J= 1,IORD+1)
              IF(MPAR.LT.0) WRITE(6,604) NPAR,NPAR,NPAR,PARM(NVARB-2),
     1                                       PARM(NVARB-1),PARM(NVARB)
              ENDIF
          NCN= NPAR
c  Loop over distance array XO(I)
          DO  I= 1,NPP
              ZZ= (XO(I)- REQ)/(XO(I)+ REQ)
              BETA= 0.d0
              DO  J= IORD,0,-1
                  BETA= BETA*ZZ+ PARM(J+1)
                  ENDDO
c  Calculate and apply switching function to MLJ exponent coefficient
              IF(MPAR.LT.0) BETA= BINF+ (BETA- BINF)/
     1                                 (1.d0+ DEXP(ALFA*(XO(I) - RX)))
              VV(I)= DSCM*(1.d0 - (REQ/XO(I))**NPAR *DEXP(-BETA*ZZ))**2
     1                                                    - DSCM+ VLIM
              ENDDO
          ENDIF
      IF(IPOTL.EQ.3) THEN
c=======================================================================
c** Generate a simple Morse, or Extended (EMO) Morse potential, or as
c  special cases, Coxon's GMO or Wei Hua's generalized Morse
c=======================================================================
          BETA= PARM(1)
          NCN= 99
          IF(LNPT.GE.0) THEN
              IF(MPAR.EQ.-1) THEN
c** Option to generate Wei Hua's extended 4-parameter Morse-type potl.
                  CSAV= PARM(2)
                  WRITE(6,605) DSCM,REQ,CSAV,BETA
                ELSE 
                  IF(NVARB.LE.1) WRITE(6,606) DSCM,REQ,BETA
                  IF(NVARB.GT.1) THEN
                      IF(MPAR.GE.0) WRITE(6,608) DSCM,REQ,NVARB-1,
     1                                      NVARB,(PARM(i),i= 1,NVARB)
                      IF(MPAR.EQ.-2) WRITE(6,610) DSCM,REQ,NVARB-1,
     1                                            (PARM(i),i= 1,NVARB)
                      ENDIF
                ENDIF
              ENDIF
c  Loop over distance array XO(I)
          DO  I= 1,NPP
c ... for Wei Hua's extended Morse function ...
              IF(MPAR.EQ.-1) THEN
                  VV(I)= DSCM*((1.d0 - DEXP(-BETA*(XO(I)-REQ)))/(1.d0 
     1                - CSAV*DEXP(-BETA*(XO(I)-REQ))))**2 - DSCM+ VLIM
                ELSE 
                  IF(NVARB.GT.1) THEN
                      ZZ= (XO(I)- REQ)/(XO(I)+ REQ)
c ... for Coxon-Hajigeorgiou "GMO" potential
                      IF(MPAR.EQ.-2) ZZ= (XO(I)- REQ)
                      BETA= 0.d0
                      DO  J= NVARB,1,-1             
                          BETA= BETA*ZZ+ PARM(J)
                          ENDDO
                      ENDIF
                  VV(I)= DSCM*(1.d0 - DEXP(-BETA*(XO(I)-REQ)))**2 
     1                                                    - DSCM+ VLIM
                ENDIF
              ENDDO
          ENDIF
      IF(IPOTL.EQ.4) THEN
c=======================================================================
c** Generate Seto-modified form of Surkus' GPEF function which includes
c  Dunham, SPF and OT forms as special cases.
c=======================================================================
          A0= DSCM
          IORD= NVARB-2
          X1= 1.d0
          FCT= PARM(NVARB-1)
          IF(DABS(FCT).GT.0.d0) THEN
              FCT= 1.d0/PARM(NVARB-1)
              DO  J=1,IORD
                  X1= X1+ PARM(J)*FCT**J
                  ENDDO
              DSCM= DSCM*X1
              ENDIF
          IF(NPAR.EQ.1) THEN
c  Cases with power =1 (including Dunham, SPF & O-T expansions).
              IF(DABS(PARM(NVARB-1)).LE.0.d0) THEN
c ... print for Dunham expansion ...
                  WRITE(6,612) PARM(NVARB),REQ,A0,NVARB-2,
     1                                          (PARM(I),I= 1,NVARB-2)
                  NCN= -99
                  CNN= 0.d0
                  ENDIF
              IF(DABS(PARM(NVARB)).LE.0.d0) THEN
c ... print for Simons-Parr-Finlan expansion ...
                  WRITE(6,614) PARM(NVARB-1),REQ,DSCM,A0,NVARB-2,
     1                                          (PARM(I),I= 1,NVARB-2)
                  NCN= 1
                  ENDIF
              IF(DABS(PARM(NVARB)-PARM(NVARB-1)).LE.0.d0) THEN
c ... print for Ogilvie-Tipping expansion ...
                  WRITE(6,616) PARM(NVARB),REQ,DSCM,A0,NVARB-2,
     1                                          (PARM(I),I= 1,NVARB-2)
                  NCN= 1
                  ENDIF
              ENDIF
          IF((NPAR.NE.1).OR.((DABS(PARM(NVARB)-PARM(NVARB-1)).GT.0.d0)
     1           .AND.(DABS(PARM(NVARB)*PARM(NVARB-1)).GT.0.d0))) THEN
c ... print for general GPEF expansion variable ...
              IF(NPAR.LT.0) THEN
c ... for negative NPAR, convert to equivalent positive NPAR case
                  NPAR= -NPAR
                  A1= PARM(NVARB)
                  PARM(NVARB)= -PARM(NVARB-1)
                  PARM(NVARB-1)= -A1
                  ENDIF
              WRITE(6,618) NPAR,NPAR,PARM(NVARB-1),NPAR,PARM(NVARB),
     1                 NPAR,REQ,DSCM,A0,NVARB-2,(PARM(I),I= 1,NVARB-2)
              NCN= NPAR
              ENDIF
          DO  I= 1, NPP
              ZZ= (XO(I)**NPAR - REQ**NPAR)/(PARM(NVARB-1)*XO(I)**NPAR
     1                                       + PARM(NVARB)*REQ**NPAR)
              A1= 1.d0
              ZP= 1.d0
              DO  J=1, NVARB-2
                  ZP= ZP*ZZ
                  A1= A1+ PARM(J)*ZP
                  ENDDO
              VV(I)= A0*ZZ*ZZ*A1 + VLIM
              ENDDO
          ENDIF
      IF(IPOTL.EQ.5) THEN
c=======================================================================
c** For generalized  H.F.D.(NPAR,6,8,10,12,14)  potential with reduced 
c  form   VBAR = ALFA*x**PARM(5) * exp[-BETR*x - PARM(4)*x**2] - D(x)*
c       [PARM(6)/x**NPAR + PARM(7)/x**6 + PARM(8)/x**8 + PARM(9)/x**10 
c       + PARM(10)/X**12 + PARM(11)/X**14]   where   x=r/R_e ,  
c  VBAR= V/epsilon   and   D(x)= exp[-PARM(1)*(PARM(2)/x - 1)**PARM(3)]
c  for  x < PARM(2)
c=======================================================================
          A1= PARM(1)
          A2= PARM(2)
          A3= PARM(3)
          B2= PARM(4)
          RC8= 0.d0
          RC10= 0.d0
          RC12= 0.d0
          RC14= 0.d0
          RCNPAR= PARM(6)
          NCN= 6
          IF(RCNPAR.GT.0.d0) NCN= NPAR
          RC6= PARM(7)
          IF(NVARB.ge.8)  RC8= PARM(8)
          IF(NVARB.ge.9)  RC10= PARM(9)
          IF(NVARB.ge.10) RC12= PARM(10)
          IF(NVARB.ge.11) RC14= PARM(11)
          DX= 1.d0
          DX1= 0.d0
          IF(A2.GT.1.d0) THEN
              DX= DEXP(-A1*(A2- 1.d0)**A3)
              DX1= A1*A2*A3*DX*(A2- 1.d0)**(A3- 1.d0)
              ENDIF
          ALFA= -1.D0+ (RCNPAR+ RC6+ RC8+ RC10+ RC12+ RC14)*DX
          IF(ALFA.LE.0.d0) THEN
              WRITE(6,622) RCNPAR,RC6,RC8,RC10,RC12,RC14,ALFA
              STOP
              ENDIF
          B1= ((NPAR*RCNPAR+6.D0*RC6+8.D0*RC8+10.D0*RC10+12.d0*RC12+
     1      14.d0*RC14)*DX - (RCNPAR+RC6+RC8+RC10+RC12+RC14)*DX1)/ALFA
     2      + PARM(5) - 2.D0*B2
          ALFA= ALFA*DEXP(B1+B2)
          IF(LNPT.GE.0) WRITE(6,624) NPAR,PARM(5),B1,B2,ALFA*DSCM,
     1                 RCNPAR,RC6,RC8,RC10,RC12,RC14,A1,A2,A3,DSCM,REQ
          DO  I= 1,NPP
              X1= XO(I)/REQ
              XP1= 0.0D0
              IF((B1*X1+ B2*X1**2).LT.170.D0) XP1= DEXP(-X1*(B1+ B2*X1))
              XP1= XP1*X1**PARM(5)
              FC1= 1.D0
              IF(A2.GT.X1) FC1= DEXP(-A1*(A2/X1- 1.d0)**A3)
              XM2C= (REQ/XO(I))**2
              VV(I)= DSCM*(ALFA*XP1- FC1*(((((RC14*XM2C+RC12)*XM2C+RC10)
     1         *XM2C+ RC8)*XM2C+ RC6)*XM2C**3 + RCNPAR/X1**NPAR)) + VLIM
              ENDDO
          ENDIF
      IF(IBOB.GT.0) THEN
c=======================================================================
c** If appropriate, generate Born-Oppenheimer breakdown correction 
c  functions to rotationless and/or centrifugal potential(s).
c [Special "Coxon" option: if  RX > 0.0, expand as per older Coxon work]
c=======================================================================
          IF(RX.GE.0.D0) RDIF= REQ-RX
          DO  I=1,NPP
              IF(RX.LE.0.d0) THEN
                  ZZ= (XO(I)-REQ)/(XO(I)+REQ)
                ELSE
                  ZZ= XO(I)- REQ
                  RD= XO(I)- RX
                ENDIF
              SC1= 0.d0
              SC2= 0.d0
              SG1= 0.d0
              SG2= 0.d0
              RPOW= 1.d0
              DO  J= 0,NCMAX
                  SC1= SC1+ RPOW*CA1(J)
                  SC2= SC2+ RPOW*CA2(J)
                  IF(RX.LE.0.d0) THEN
                      SG1= SG1+ RPOW*GA1(J)
                      SG2= SG2+ RPOW*GA2(J)
                    ELSE 
                      M= J-1
                      SG1= SG1+ (RD**J -RDIF**J)*GA1(J)
                      SG2= SG2+ (RD**J -RDIF**J)*GA2(J)
                    ENDIF
                  RPOW= RPOW*ZZ
                  ENDDO
              RM2(I)= (1.d0+ SG1+ SG2)/XO(i)**2
              VV(I)= VV(I) + SC1 + SC2
              ENDDO
          ENDIF
      RETURN
  600 FORMAT(/' Lennard-Jones(',I2,',',I2,') potential with   De=',
     1  F10.3,'(cm-1)   Re =',F10.6,'(A)')
  602 FORMAT(/' Use an MLJ potential with   De =',F10.3,
     1  '(cm-1)    Re =',F12.8,'(A)'/3x,'with exponent parameter BETA an
     2 order-',i2,' polynomial in  z=(R-Re)/(R+Re)  with'/
     3  '   coefficients:',1PD16.8,3D16.8:/(5D16.8:))
  604 FORMAT(' & exponent switching function yielding limiting C',i1,
     1 '/R^',i1,' with   C_',i1,'=',1PD13.6/10x,'defined by   ALPHA_s=',
     2  0Pf9.6,'   R_s=',f10.6)
  605 FORMAT(/' Potential is a Hua-Wei 4-parameter Morse type function w
     1ith   De =',F11.4/11x,'Re =',F12.9,'   C=',f7.4,'   &   beta=',
     1  F13.10,' [1/Angstroms]')
  606 FORMAT(/' Potential is a simple Morse function with   De =',F11.4,
     1  '    Re =',F12.9/39x,'and   beta =',F13.10,' [1/Angstroms]')
  608 FORMAT(/' Potential is Extended Morse Oscillator with   De=',
     1  F11.4,'    Re=',F12.9/5x,'Exponent factor "beta" is order-',i2,
     2 ' power series in  z=(R-Re)/(R+Re)  with'/5x,I2,' coefficients:',
     3  1x,1PD18.9,2D18.9:/(4X,4D18.9:))
  610 FORMAT(/' Potential is Generalized Morse Oscillator with   De=',
     1 F10.3,'   Re=',F11.8/4x,'Exponent factor "beta" is',i3,' order po
     2wer series in (R-Re) with coefficients:'/4x,1PD18.9,3D18.9:/
     3 (4X,4D18.9:))
  612 FORMAT(/' Potential is a Dunham expansion in  (R-Re)/(',f5.2,
     1  ' * Re)  with   Re=',f12.9/'  a0=',1PD16.9,'   and',i3,
     2  '  a_i coefficients:'/(5D16.8))
  614 FORMAT(/' Potential is an SPF expansion in  (R-Re)/(',F5.2,
     1  '* R)  with   Re=',f12.9/5x,'De=',g18.10,'   b0=',
     2  1PD16.9,'   and',i3,'  b_i  coefficients:'/(5D16.8))
  616 FORMAT(/' Potential is an O-T expansion in  (R-Re)/[',f5.2,
     1  '*(R+Re)]  with   Re=',f12.9/5x,'De=',G18.10,
     2  '   c0=',1PD16.9,'   and',i3,'  c_i coefficients:'/(5D16.8))
  618 FORMAT(/' Potential is a general GPEF expansion in  (R^',i1,
     1  ' - Re^',i1,')/(',SP,F5.2,'*R^',SS,i1,SP,F6.2,'*Re^',SS,i1,')'/
     2  5x,'with   Re=',f12.9,'   De=',g18.10,'   g0=',1PD16.9/
     2  5x,'and',i3,'  g_i coefficients:  ',3D16.8/(5D16.8:))
  622 FORMAT(/' *** ERROR in generating HFD potential *** C',i1,
     1 ', C6, C8, C10, C12, C14  =',6G15.7/10X,'yield    ALFA =',G15.7)
  624 FORMAT(/' Potential is Generalized HFD(',i1,',6,8,10,12,14) with',
     1 '  gamma=',f9.6/'    beta1=',f12.8,'    beta2=',f9.6,'    A=',
     2 1PD16.9/"    reduced {Cn's}:",3D14.6/19x,3D14.6/'    Damping func
     3tion  D(R)= exp[ -',0Pf6.4,'*(',f7.4,'/X -1.0)**',f5.2,']' /
     4 '  & DSCM=',f10.4,'[cm-1]   Re=',f9.6,'[Angst.]')
  628 FORMAT(' ')
  630 FORMAT(' B-O-B correction to rotationless potential for atom-',
     1  I1,'  of mass ',f14.10/5x,'is  [order-',I2,' polynomial in {(R-R
     2e)/(R+Re)}] * [1- MASS(',A2,i3,')/MASS(',A2,I3,')]'/5x,'with',i3,
     3  ' coefficients:',3G18.10:/(8x,4G18.10:))
  632 FORMAT(' B-O-B correction to rotationless potential for atom-',
     1 I1,'  of mass ',f14.10/5x,'is  [order-',I2,' polynomial in (R-Re)
     2]/[MASS(',A2,I3,')]   with',i3,' coefficients:'/(5x,4G18.10:))
  634 FORMAT(' B-O-Breakdown correction to centrifugal term for atom-',
     1 I1,'  of mass ',f14.10/5x,'is  [order-',I2,' polynomial in {(R-Re
     2/(R+Re)}] * [MASS(',A2,I3,')/MASS(',A2,I3,')]'/5x,'with',i3,' coef
     3ficients:',3G18.10:/(8x,4G18.10:))
  636 FORMAT(' B-O-Breakdown correction to centrifugal term for atom-',
     1 I1,'  of mass ',f14.10/5x,'is  [order-',I2,' polynomial in {(R-Rx
     2)**i - (Re-Rx)**i}] / [MASS(',A2,I3,')]'/5x,'with  Rx=',f6.3,
     3 '  &',i3,' coefficients: ',1PD18.9,D18.9:/(5x,4D18.9:))
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE MAPPAR(NSTATES,PV,FORBAC)
c***********************************************************************
c** This subroutine will convert external logical physical parameters 
c  into the generic NLLSSRR parameter array PV or the reverse.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++ COPYRIGHT 1997-2007  by  J.Y. Seto & R.J. Le Roy  (ver. 12/04/2007)+
c   Dept. of Chemistry, Univ. of Waterloo, Waterloo, Ontario, Canada   +
c    This software may not be sold or any other commercial use made    +
c      of it without the express written permission of the authors.    +
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c    NSTATES is the number of states being considered.
c    PV(i)   is the NLLSSRR parameter array.
c    FORBAC  is a flag to determine which way the parameters are mapped
c            FORBAC = 0 : Map internal PV to external variables.
c            FORBAC = 1 : Map external varuables to internal PV.
c=======================================================================
      INCLUDE 'arrsizes.h'
      INCLUDE 'BLKPOT.h'
      INCLUDE 'BLKBOB.h'
c-----------------------------------------------------------------------
      INTEGER NSTATES, m, FORBAC, ISTATE, IPV, I
      REAL*8 PV(NPARMX)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Map external free parameters (De, Re, etc.) onto internal NLLSSRR
c  parameters PV(j)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IF(FORBAC.EQ.0) THEN
          IPV= 0
          DO  ISTATE= 1,NSTATES
              IF(PSEL(ISTATE).NE.-2) THEN
                  IF(PSEL(ISTATE).NE.4) THEN
                      IPV= IPV+ 1
                      PV(IPV)= DE(ISTATE)
                      ENDIF
                  IPV= IPV+ 1
                  PV(IPV)= RE(ISTATE)
                  IF(PSEL(ISTATE).EQ.2) THEN
                      DO  m= 1,NCMM(ISTATE)
                          IPV= IPV+ 1
                          PV(IPV)= CmVAL(m,ISTATE)
                          ENDDO
                      ENDIF
                  DO  I= 0,MAX(NSphi(ISTATE),NLphi(ISTATE))
                      IPV= IPV+ 1
                      PV(IPV)= PHI(I,ISTATE)
                      ENDDO
                  IF(NUA(ISTATE).GE.0) THEN
                      DO  I= 0,NUA(ISTATE)
                          IPV= IPV+ 1
                          PV(IPV) = UA(I,ISTATE)
                          ENDDO
                      ENDIF
                  IF(NUB(ISTATE).GE.0) THEN
                      DO  I= 0,NUB(ISTATE)
                          IPV= IPV+ 1
                          PV(IPV) = UB(I,ISTATE)
                          ENDDO
                      ENDIF
                  IF(NTA(ISTATE).GE.0) THEN
                      DO  I= 0,NTA(ISTATE)
                          IPV= IPV+ 1
                          PV(IPV) = TA(I,ISTATE)
                          ENDDO
                      ENDIF
                  IF(NTB(ISTATE).GE.0) THEN
                      DO  I= 0,NTB(ISTATE)
                          IPV= IPV+ 1
                          PV(IPV) = TB(I,ISTATE)
                          ENDDO
                      ENDIF
                  IF(NwCFT(ISTATE).GE.0) THEN
                      DO  I= 0, NwCFT(ISTATE)
                          IPV= IPV+ 1
                          PV(IPV) = wCFT(I,ISTATE)
                          ENDDO
                      ENDIF
                  ENDIF
              ENDDO
        ELSEIF(FORBAC.EQ.1) THEN
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Convert internal NLLSSRR parameter array back into external
c   (logical) variable system (De, Re, etc.).
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          IPV = 0
          DO  ISTATE=1,NSTATES
              IF(PSEL(ISTATE).NE.-2) THEN
                  IF(PSEL(ISTATE).NE.4) THEN
                      IPV= IPV + 1
                      DE(ISTATE)= PV(IPV)
                      ENDIF
                  IPV= IPV + 1
                  RE(ISTATE) = PV(IPV)
                  IF(PSEL(ISTATE).EQ.2) THEN
                      DO  m= 1,NCMM(ISTATE)
                          IPV= IPV+ 1
                          CmVAL(m,ISTATE)= PV(IPV) 
                          ENDDO
                      ENDIF
                  DO I=0, MAX(NSphi(ISTATE),NLphi(ISTATE))
                      IPV = IPV + 1
                      PHI(I,ISTATE) = PV(IPV)
                      ENDDO
                  IF(NUA(ISTATE).GE.0) THEN
                      DO I= 0,NUA(ISTATE)
                          IPV = IPV + 1
                          UA(I,ISTATE) = PV(IPV)
                          ENDDO
                      ENDIF
                  IF(NUB(ISTATE).GE.0) THEN
                      DO I= 0,NUB(ISTATE)
                          IPV = IPV + 1
                          UB(I,ISTATE) = PV(IPV)
                          ENDDO
                      ENDIF
                  IF(NTA(ISTATE).GE.0) THEN
                      DO I= 0,NTA(ISTATE)
                          IPV = IPV + 1
                          TA(I,ISTATE) = PV(IPV)
                          ENDDO
                      ENDIF
                  IF(NTB(ISTATE).GE.0) THEN
                      DO I= 0,NTB(ISTATE)
                          IPV = IPV + 1
                          TB(I,ISTATE) = PV(IPV)
                          ENDDO
                      ENDIF
                  IF(NwCFT(ISTATE).GE.0) THEN
                      DO I=0,NwCFT(ISTATE)
                          IPV = IPV + 1
                          wCFT(I,ISTATE) = PV(IPV)
                          ENDDO
                      ENDIF
                  ENDIF
              ENDDO
        ENDIF
      RETURN
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE ALF(NDP,RMIN,RH,V,SWF,VLIM,KVMAX,AFLAG,ZMU,EPS,NCN,GV,
     1                                                 INNODE,INNR,IWR) 
c***********************************************************************
c   Version 1.0 dated July 6, 2003. Cosmetic formatting changes 11/04/07
c-----------------------------------------------------------------------
c** The subroutine ALF (Automatic vibrational Level Finder) will
c   automatically generate the eigenvalues from the first vibrational
c   level (v=0) to a user specified level (v=KVMAX) or the highest
c   allowed vibrational level of a given smooth single (or double)
c   minimum potential (V). These energies are stored and returned to the
c   calling program in the molecular constants array GV(v=0-KVMAX).
c** For any errors that cannot be resolved within the subroutine, ALF
c   returns AFLAG with a value that defines which error had occured.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++   COPYRIGHT 1998 - 2003  by  Jenning Seto and Robert J. Le Roy   +++
c   Dept. of Chemistry, Univ. of Waterloo, Waterloo, Ontario, Canada   +
c    This software may not be sold or any other commercial use made    +
c     of it without the express written permission of the authors.     +
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+ Please inform me of any bugs, by phone at: (519)888-4567, ext. 4051 +
c++++++++ by e-mail to: leroy@uwaterloo.ca , or write me at: ++++++++++
c+++ Dept. of Chemistry, Univ. Waterloo, Waterloo, Ontario  N2L 3G1 ++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Based on the automatic level finding routine from LEVEL 6.0 written 
c   by Robert J. Le Roy
c** Uses the Schrodinger solver subroutine SCHRQ.
c
c** On entry:
c    NDP    is the number of datapoints used for the potential.
c    RMIN   is the inner radial distance of the potential (ang).
c    RH     is the meshvalue (ang).
c           NDP, RMIN, and RH define the radial range over which the
c           potential is defined.
c    V(i)   is the scaled input potential (cm-1).
c           The scaling factor BFCT is (2*mu/hbar^2)*RH^2.
c    VLIM   is the potential asymptote (cm-1).
c    KVMAX  is v for the highest vibrational level we wish to find.
c    AFLAG  is rot.quantum J for the (centrifugally distorted) potential
c    ZMU    is the reduced mass of the diatom (amu).
c    EPS    is the energy convergence criterion (cm-1).
c    NCN    is the near dissociation limit radial exponential.
c    INNODE specifies whether wave fx. initiation @ RMIN starts with a
c        note (normal case: INNODE > 0) or zero slope (when INNODE.le.0)
c    IWR    specifies the level of printing inside SCHRQ
c           <> 0 : print error & warning descriptions.
c           >= 1 : also print final eigenvalues & node count.
c           >= 2 : also show end-of-range wave function amplitudes.
c           >= 3 : print also intermediate trial eigenvalues, etc.
c
c** On exit:
c    KVMAX   is vib.quantum number for the highest allowed vibrational 
c            level found (may be less than the input value of KVMAX).
c    AFLAG   returns calculation outcome to calling program.
c            >=  0 : Subroutine found all levels to v=KVMAX(input).
c             = -1 : KVMAX larger than number of allowed levels.
c             = -2 : Initial trial energy is unusable.
c             = -3 : Calculated trial energy is unusable.
c             = -4 : Cannot find first vibrational level.
c             = -5 : Calculated trial energy too low.
c             = -6 : Calculated trial energy too high.
c             = -7 : An impossible situation occured.
c             = -8 : Potential found to have a second minimum.
c    GV(v)   contains the vibrational energy levels found for v=0-KVMAX
c    INNR(v) labels each level as belonging to the inner (INNR = 1) or
c            outer (INNR = 0) well.
c
c** Flags: Modify only when debugging.
c    AWO   specifies the level of printing inside ALF
c          <> 0 : print error & warning descriptions.
c          >  0 : also print intermediate ALF messages.
c    MCO   specifies the level of printing of molecular constants.
c          >  0 : print out vibrational energies to channel-21.
c    INNER specifies wave function matching (& initiation) conditions.
c        .le.0 : Match inward & outward solutions at outermost well t.p.
c          > 0 : Match at innermost well inner turning point
c        For most normal cases set INNER = 0,  but ......
c            To find "inner-well-dominated" solutions of an asymmetric
c            double minimum potential, set  INNER > 0.
c    LPRWF specifies option of printing out generated wavefunction
c          > 0 : print wave function every LPRWF-th  point.
c          < 0 : compactly write to channel-7 every |LPRWF|-th wave
c                function value.
c          A lead "card" identifies the level, gives the position of
c          1-st point and radial mesh, & states No. of  points.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** The dimensioning parameters must be consistant with the sizes of the
c   arrays used in the calling program.
c
c    NVIBMX  is the maximum number of vibrational levels considered.
c            Note: NVIBMX should be larger than KVMAX.
c
      INTEGER NVIBMX
      PARAMETER (NVIBMX= 400)
c
c** NF counts levels found in automatic search option
c
c** OWL holds the vibrational levels that are contained in the outer
c   well.
c** IWL holds the vibrational levels that are contained in the inner
c   well (if present).
c
      INTEGER NDP,KVMAX,NCN,KV,AFLAG,NF,NBEG,NEND,INNR(0:KVMAX),IWR,
     1  I,IZPE,IVDIF,IVCOR,IQT,IEG,LTRY,TRIALE,AWO,MCO,INNODE,INNER,
     2  LPRWF,JROT,NPMIN, NPMAX, NIWL,IWL(0:NVIBMX),NOWL,OWL(0:NVIBMX)
c
      REAL*8 RMIN,RMAX,RH,V(NDP),SWF(NDP),VLIM,EO,ZMU,EPS,LHIE,LLOE,
     1  BZ,BFCT,PW,PWI,GAMA,VMIN,VMAX,RE,PMAX,VDMV,VDL,VDU,DRAVOUT,
     2  DRAVIN,AO,VD,GV(0:KVMAX),RAVG(0:NVIBMX),VPMIN(10),RPMIN(10),
     3  VPMAX(10),RPMAX(10)
c
      DATA AWO/-1/,MCO/0/,LPRWF/0/
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Check that the array dimensions are adequate.
      IF(KVMAX.GT.NVIBMX) THEN
          WRITE(6,610)
          WRITE(6,613) KVMAX, NVIBMX
          STOP
          ENDIF
c
c** Initialize level counters for each well.
      DO  I= 0,KVMAX
          INNR(I)= -1
          IWL(I)= 0
          OWL(I)= 0
          END DO
c
c** Initialize the remaining variables and flags.
      NF= 0
      NIWL= 0
      NOWL= 0
      KV= 0
      INNER= 0
      LTRY= 0
      LHIE= 1.0d99
      CALL INITVAL(IQT,IVDIF,IZPE,IEG,IVCOR)
c
c** Store rotational quantum number.
      JROT= AFLAG
c
c** Numerical factor  16.85762920 (+/- 0.00000011) based on Compton
c  wavelength of proton & proton mass (u) from 2002 physical constants.
      BZ= ZMU/16.85762920d0
      BFCT= BZ*RH*RH
c
c** RMAX is the outer radial distance over which potential is defined. 
      RMAX= RMIN + DBLE(NDP-1)*RH
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Locate the potential minima.
      NPMIN= 0
      DO  I= 2,NDP-1
          IF((V(I).LT.V(I-1)).AND.(V(I).LT.V(I+1))) THEN
              NPMIN= NPMIN + 1
              RPMIN(NPMIN)= RMIN + DBLE(I-1)*RH
              VPMIN(NPMIN)= V(I) / BFCT
              IF(NPMIN.EQ.10) GOTO 100
              ENDIF
          END DO
c** If NO minimum can be found, then print a warning and stop.
  100 IF(NPMIN.EQ.0) THEN
          IF(V(2).LE.V(1)) THEN
              WRITE(6,614) JROT
              KVMAX= -1
              RETURN
              ENDIF
          NPMIN= 1
          VPMIN(NPMIN)= V(1)/BFCT
          RPMIN(NPMIN)= RMIN
          WRITE(6,618) VPMIN(1),RMIN
          ENDIF
c** If more than two minima are found, then print a warning and hope ...
      IF(NPMIN.GT.2) THEN
          WRITE(6,605)
          WRITE(6,615) NPMIN, 'minima'
c         STOP
          ENDIF
      NPMAX= 0
c** Locate the potential maxima (if it exists).
      DO  I= 2,NDP-1
          IF((V(I).GT.V(I-1)).AND.(V(I).GT.V(I+1))) THEN
              NPMAX= NPMAX + 1
              RPMAX(NPMAX)= RMIN + DBLE(I-1)*RH
              VPMAX(NPMAX)= V(I) / BFCT
              IF(NPMAX.EQ.10) GOTO 150
              ENDIF
          END DO
c
  150 IF(NPMAX.EQ.0) THEN
c** If no maxima were found, then set energy maximum to be the value 
c   at the end of the radial range.
          NPMAX= 1
          RPMAX(NPMAX)= RMAX
          VPMAX(NPMAX)= V(NDP) / BFCT
          ENDIF
c
c** If more than three maxima found, then print a warning and hope ...
      IF(NPMAX.GT.3) THEN
          WRITE(6,605)
          WRITE(6,615) NPMAX, 'maxima'
c         STOP
          ENDIF
c
c** If there is no barrier to dissociation, then set the
c   final VPMAX to be the value at the end of the range. [huh ????]
      IF(RPMAX(NPMAX).LT.RPMIN(NPMIN)) THEN
          NPMAX= NPMAX + 1
          RPMAX(NPMAX)= RMAX
          VPMAX(NPMAX)= V(NDP) / BFCT
          ENDIF
c
c** If innermost maximum occurs before innermost minimum, the potential 
c   turns over in short range region and should not be used.  
c   Print a warning and STOP.
      IF(RPMAX(1).LT.RPMIN(1)) THEN
          WRITE(6,610)
          WRITE(6,616) RPMAX(1)
          STOP
          ENDIF
c
c** Now find the absolute potential minimum.
      VMIN= VPMIN(1)
      RE= RPMIN(1)
      DO  I= 2,NPMIN
          IF(VMIN.GT.VPMIN(I)) THEN
              VMIN= VPMIN(I)
              RE= RPMIN(I)
              ENDIF
          END DO 
c
c** Now find the absolute potential maximum.
      VMAX= VPMAX(1)
      DO  I= 2,NPMAX
          IF(VMAX.LT.VPMAX(I)) VMAX= VPMAX(I)
          END DO 
c
c** If the absolute potential maximum is lower than the absolute
c   potential minimum, then print out an error statement and quit.
      IF(VMAX.LE.VMIN) THEN
          WRITE(6,610)
          WRITE(6,617)
          STOP
          ENDIF
c
c** Otherwise, print out potential extrema count, if desired.
      IF(AWO.GT.0) THEN
          WRITE(6,650) NPMIN, VMIN
          WRITE(6,651) NPMAX, VMAX
          ENDIF
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Calculate 2*NCN/(NCN - 2) for use when calculating trial energies. 
      PW= 20.0d0
      IF((NCN.GT.0).AND.(NCN.NE.2)) PW= 2.d0*DBLE(NCN)/(DBLE(NCN)-2.d0)
      IF(VMAX.GT.VLIM) PW= 2.d0
      PWI= 1.d0/PW
c
c** Use Lennard-Jones estimation of zero point energy to determine the
c   initial trial energy.
c                            _____________________________
c            vD + 0.5 = ao \/ZMU * De * Re^2 / 16.85762920
c                 De = A (vD - v)^3 = A (vD + 0.5)^3
c              E(v=0) = VMIN + A [(vD + 0.5)^3 - vD^3]
c** Choose AO to have a value of 0.25.
c
      AO= 0.25d0
      VD= AO*DSQRT(BZ*(VMAX-VMIN))
      VD= VD*RE - 0.5d0
      AO= (VMAX-VMIN)*(1.d0 - (VD/(VD+0.5d0))**3)
      EO= VMIN + AO
c
      IF(MCO.GE.1) THEN
c** If desired, write out SCHRQ control information.
          WRITE(21,2100)
          WRITE(21,2110) RMIN, RMAX, RH, BZ, ZMU
          WRITE(21,2111) EPS
          WRITE(21,2112)
          WRITE(21,2101)
          ENDIF
c=========== Begin Actual Eigenvalue Calculation Loop Here =============
c** Compute eigenvalues ... etc. up to the KVMAX'th vibrational level.
c** When attempts to find the next eigenvalue fails, then perhaps the
c   next level is located in a second (inner) well. If so, then the
c   subroutine will set INNER = 1, and attempt to find that level.
c
   10 IF(AWO.GT.0) THEN
          IF(INNER.EQ.0) WRITE(6,601)
          IF(INNER.EQ.1) WRITE(6,603)
          ENDIF
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Call subroutine SCHRQ to find eigenvalue EO and eigenfunction SWF(I).
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      CALL SCHRQ(KV,JROT,EO,GAMA,PMAX,VLIM,V,SWF,BFCT,EPS,RMIN,RH,NDP,
     1                               NBEG,NEND,INNODE,INNER,IWR,LPRWF)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** If a level higher than NF has already been found and the level found
c   not the correct level, ignore all errors and continue with the 
c   binary search.
c
      IF((KV.GE.0).AND.(KV.LT.NF).AND.((IVCOR.GT.0).OR.(IEG.GT.3))) THEN
          KV= NF + 1
          EO= VMAX + 1.0d0
          IF(IVCOR.GT.0) THEN
              LLOE= (LHIE + LLOE) / 2.0d0
            ELSE
              LLOE= GV(NF-1)
              LHIE= VMAX
              IVCOR= 1
            ENDIF
          TRIALE= 0
          ENDIF
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** The SCHRQ error condition is KV < 0.
c   There are three possible situations to consider:
c     EO > VMAX : Trial energy greater than potential maximum
c     NF = 0 : Looking for the first vibrational level (v = 0)
c     NF > 0 : Looking for the other vibrational levels (v > 0)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** For the case when the next trial energy is higher than the potential
c   maximum, try one last ditch attempt to find the highest bound level
c   (quasi or otherwise) in the potential.
c
      IF((KV.LT.0).AND.(EO.GT.VMAX)) THEN
          IF(IVCOR.LT.20) THEN
              KV= NF + 1
              IF(IVCOR.EQ.0) THEN
                  LLOE= GV(NF-1)
                  LHIE= VMAX
                  IVCOR= 1
                  ENDIF
              ENDIF
          IF(LTRY.LT.1) THEN
              LTRY= 1
              KV= 999
              EO= VMAX - 1.0d-2
c
c** If unsuccessful, then print out a warning and exit.
            ELSE
              IF(AWO.NE.0) THEN
                  WRITE(6,605)
                  WRITE(6,606) NF, EO, VMAX
                  ENDIF
              AFLAG= -1
              GOTO 200
            ENDIF
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** If trying to find the first vibrational level (v=0), then double the
c   zero point energy estimation (AO):   E(v=0) = VMIN + IQT*AO
c
        ELSEIF((KV.LT.0).AND.(NF.EQ.0)) THEN
          IF(IQT.GT.1) THEN
              IF(AWO.NE.0) THEN
                  WRITE(6,610)
                  WRITE(6,611)
                  WRITE(6,620) IQT, EO
                  ENDIF
c
c** If this fails, then try changing the wavefunction matching
c   condition (INNER) to see if a possible second minimum contains the
c   zero-point level.
              IF(INNER.EQ.0) THEN
                  INNER= 1
                  CALL INITVAL(IQT,IVDIF,IZPE,IEG,IVCOR)
c** If both attempts fail, then print out warning message and exit the
c   subroutine.
                ELSE
                  AFLAG= -2
                  GOTO 200
                ENDIF
              ENDIF
          IQT= IQT + 1
          EO= VMIN + DBLE(IQT)*AO
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** If trying to find other vibrational levels (v > 0) then switch to
c   use of differences for estimating spacing.
c
        ELSEIF((KV.LT.0).AND.(NF.GT.0)) THEN
          IF(IVDIF.GT.0) THEN
              IF(AWO.NE.0) THEN
                  WRITE(6,610)
                  WRITE(6,612)
                  WRITE(6,621) NF,IVDIF
                ENDIF
c
c** If differences fails, then try changing the wavefunction matching
c   condition (INNER) to see if a possible second minimum contains the
c   zero point level.
c
              IF(INNER.EQ.0) THEN
                  INNER= 1
                  CALL INITVAL(IQT,IVDIF,IZPE,IEG,IVCOR)
c
c** If both attempts fail, then print out warning message and exit the
c   subroutine.
                ELSE
                  AFLAG= -3
                  GOTO 200
                ENDIF
              ENDIF
          IVDIF= 1
          IF(INNER.EQ.0) THEN
              CALL DTENG(IEG,NF,NOWL,OWL,NVIBMX,VMIN,GV,EO)
            ELSE
              CALL DTENG(IEG,NF,NIWL,IWL,NVIBMX,VMIN,GV,EO)
            ENDIF
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** If first level found isn't v=0, try up to 3 times to 'harmonically'
c   estimate improved trial ground state energy.
c           E(v=0)= E(v=KV) - (E(v=KV) - VMIN)/(1 + KV/2)
c
        ELSEIF((KV.GT.0).AND.(NF.EQ.0)) THEN
          IF(IZPE.GT.3) THEN
              IF(AWO.NE.0) THEN
                  WRITE(6,610)
                  WRITE(6,611)
                  WRITE(6,622) IZPE,GV(0),KV,EO
                  ENDIF
c
c** If differences fails, then try changing the wavefunction matching
c   condition (INNER) to see if a possible second minimum contains the
c   zero point level.
              IF(INNER.EQ.0) THEN
                  INNER= 1
                  CALL INITVAL(IQT,IVDIF,IZPE,IEG,IVCOR)
c
c** If both attempts fail, then print out warning message and exit the
c   subroutine.
                ELSE
                  AFLAG= -4
                  GOTO 200
                ENDIF
              ENDIF
          IZPE= IZPE + 1
          EO= EO - (EO-VMIN)/(1.d0+0.5d0/KV)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** For the next three cases, KV >= 0 and NF > 0.
c** If the calculated vibrational level is less than the next expected
c   level, then the estimated trial energy is too low.
c** Perhaps the difference in energy between vibrational levels v and
c   v-1 is much greater than the energy between levels v-1 and v-2. 
c                 E(v) - E(v-1) >> E(v-1) - E(v-2)
c   In which case (most likely a potential with a shelf), try twice to
c   estimate a higher trial energy.
c
c   E(v)= E(v-1) + (1+IEG/2) * (2*(E(v-1)-E(v-2)) - (E(v-2)-E(v-3)))
c
        ELSEIF(KV.LT.NF) THEN
          IF(IEG.GT.3) THEN
              IF(AWO.NE.0) THEN
                  WRITE(6,610)
                  WRITE(6,612)
                  WRITE(6,623) NF, KV
                  ENDIF
c
c** If this fails, then try changing the wavefunction matching
c   condition (INNER) to see if a possible second minimum contains the
c   zero point level.
              IF(INNER.EQ.0) THEN
                  INNER= 1
                  CALL INITVAL(IQT,IVDIF,IZPE,IEG,IVCOR)
c
c** If both attempts fail, then print out warning message and exit the
c   subroutine.
                ELSE
                  AFLAG= -5
                  GOTO 200
                ENDIF
              ENDIF
          IEG= IEG + 1
c
c** If a second minimum is present, then the next vibrational level may
c   be in the inner well. If so, use the inner well vibrational levels
c   to estimate the next trial energy.
          IF(INNER.EQ.0) THEN
              CALL DTENG(IEG,NF,NOWL,OWL,NVIBMX,VMIN,GV,EO)
            ELSE
              CALL DTENG(IEG,NF,NIWL,IWL,NVIBMX,VMIN,GV,EO)
            ENDIF
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** If the calculated vibrational level is the next expected level, ...
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ELSEIF(KV.EQ.NF) THEN
          GV(NF)= EO
          INNR(NF)= INNER
          LTRY= 0
          TRIALE= 1
          LLOE= EO
          LHIE= VMAX
          CALL INITVAL(IQT,IVDIF,IZPE,IEG,IVCOR)
c-----------------------------------------------------------------------
c
c** To ease confusion when using a potential with two minima, use
c  <R>=RAVG  to keep track of which is which.
c
   20     RAVG(NF)= 0.5d0*(SWF(NBEG)*(RMIN+DBLE(NBEG-1)*RH)
     1                        + SWF(NEND)*(RMIN+DBLE(NEND-1)*RH))
          DO  I= NBEG+1,NEND-1
              RAVG(NF)= RAVG(NF) + SWF(I)*(RMIN+DBLE(I-1)*RH)
              END DO
          RAVG(NF)= RAVG(NF)*RH
c
c** Double check that the calculated level is in fact located in the
c   correct well. This can be done (for v.ne.0) by comparing <R> for
c   the new level and with previous values in each well. If the
c   difference is greater than 1.5 times the difference in the other
c   well, assume the calculated level is probably in the wrong well.
c
          IF(NOWL.GT.0) THEN
              DRAVOUT= DABS(RAVG(NF) - RAVG(OWL(NOWL-1)))
            ELSE
              DRAVOUT= 9999.9d0
            ENDIF
          IF(NIWL.GT.0) THEN
              DRAVIN= DABS(RAVG(NF) - RAVG(IWL(NIWL-1)))
            ELSE
              DRAVIN= 9999.9d0
            ENDIF
          IF(INNER.EQ.0) THEN
              IF((NOWL.GT.0).AND.(DRAVOUT.GT.(1.5d0*DRAVIN))) THEN
                  IF(MCO.GE.1) 
     1             WRITE(21,2113) NF,'Inner',NIWL,GV(NF)-VMIN,RAVG(NF)
                  IWL(NIWL)= NF
                  NIWL= NIWL + 1
                ELSE
                  IF(MCO.GE.1) 
     1             WRITE(21,2113) NF,'Outer',NOWL,GV(NF)-VMIN,RAVG(NF)
                  OWL(NOWL)= NF
                  NOWL= NOWL + 1
                ENDIF
            ELSE
              IF((NIWL.GT.0).AND.(DRAVIN.GT.(1.5d0*DRAVOUT))) THEN
                  IF(MCO.GE.1) 
     1             WRITE(21,2113) NF,'Outer',NOWL,GV(NF)-VMIN,RAVG(NF)
                  OWL(NOWL)= NF
                  NOWL= NOWL + 1
                ELSE
                  IF(MCO.GE.1) 
     1             WRITE(21,2113) NF,'Inner',NIWL,GV(NF)-VMIN,RAVG(NF)
                  IWL(NIWL)= NF
                  NIWL= NIWL + 1
                ENDIF
              INNER= 0
            ENDIF
c
c** Look for the next uncalculated level.
          NF= NF + 1
          IF(NF.LE.KVMAX) THEN
              IF(INNR(NF).GE.0) THEN
c ... if this next level was found earlier, its INNER > initial -1 value
c  so skip up one, but first calculate its Bv on the way.
                  GOTO 20
                  ENDIF
              ENDIF
c-----------------------------------------------------------------------
c** Now estimate trial energy for next higher vibrational energy level
c   by using the Near-Dissociation Theory result that:
c                  (binding energy)**((NCN-2)/(2*NCN))
c   is (at least locally) linear in vibrational quantum number.
c
          IF(NF.EQ.1) THEN
              VDMV= 0.5d0/(((VMAX-VMIN)/(VMAX-GV(0)))**PWI - 1.d0)
            ELSE
              VDMV= 1.d0/(((VMAX-GV(NF-2))/(VMAX-GV(NF-1)))**PWI- 1.d0)
            ENDIF
c
c** If unable to calculate the next trial energy, see if all of the
c   desired levels have been calculated. If not then turn on the warning
c   flag and quit, otherwise print out success message and quit.
c
          IF((VDMV.LT.1.d0).AND.(NCN.GT.2)) THEN
              IF(NF.LE.KVMAX) THEN
                  AFLAG= -1
                  WRITE(6,640) JROT, NF - 1 + VDMV
                ELSEIF(AWO.GT.0) THEN
                  WRITE(6,630) KVMAX
                ENDIF
              GOTO 200
              ENDIF
c
c** Now calculate the next trial energy.
c
          EO= VMAX - (VMAX-GV(NF-1))*(1.d0-1.d0/VDMV)**PW
c
c** However, if the level is above the dissociation limit (for
c   potentials with barriers) then use differences to calculate the
c   next trial energy.
c
          IF(EO.GT.VMAX) CALL DTENG(IEG,NF,NOWL,OWL,NVIBMX,VMIN,GV,EO)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** If the calculated vibrational level is higher then the next expected
c   level, then try once to interpolate harmonically for the missed
c   level:      E(v)= E(v-1) + (E(KV) - E(v-1)) / 2
c
        ELSEIF(KV.GT.NF) THEN
c
c** Record vibrational level (if haven't already) for posterity.
c
          IF((KV.LE.KVMAX).AND.(EO.LT.VMAX)) THEN
              IF(INNR(KV).EQ.-1) THEN
                  GV(KV)= EO
                  INNR(KV)= INNER
                  ENDIF
              ENDIF
          IF(IVCOR.GT.19) THEN
              IF(AWO.NE.0) THEN
                  WRITE(6,610)
                  WRITE(6,612)
                  WRITE(6,624) IVCOR,KV,EO,(NF-1),GV(NF-1)
                  ENDIF
c
c** If interpolation fails, then try changing the wavefunction matching
c   condition (INNER) to see if a possible second minimum contains the
c   missing level. 
              IF(INNER.EQ.0) THEN
                  INNER= 1
                  LLOE= GV(NF-1)
                  LHIE= VMAX
                  CALL INITVAL(IQT,IVDIF,IZPE,IEG,IVCOR)
c
c** If both attempts fail, then print a warning message and exit 
                ELSE
                  AFLAG= -6
                  GOTO 200
                ENDIF
              ENDIF
c
c** Use NDE theory to predict eigenvalue for the missing level.
          IF(IVCOR.EQ.0) THEN
              VDU= (VMAX-EO)**PWI
              VDL= (VMAX-GV(NF-1))**PWI
              EO= VMAX - (VDL + (VDU - VDL) / DBLE(KV - NF + 1))**PW
c         IF((NPMIN.EQ.1).AND.(IVCOR.EQ.0)) THEN
c             VDU= (VPMAX(1)-EO)**PWI
c             VDL= (VPMAX(1)-GV(OWL(NOWL-1)))**PWI
c             EO= VPMAX(1)- (VDL+ (VDU - VDL)/ DBLE(KV - NF + 1))**PW
c           ELSEIF((((INNER.EQ.1).AND.(NIWL.GT.0)).OR.
c    1          ((INNER.EQ.0).AND.(NOWL.EQ.0))).AND.(IVCOR.EQ.0)) THEN
c             VDU= (VPMAX(1)-EO)**PWI
c             VDL= (VPMAX(1)-GV(IWL(NIWL-1)))**PWI
c             EO= VPMAX(1)- (VDL+ (VDU - VDL)/ DBLE(KV - NF + 1))**PW
c           ELSEIF((((INNER.EQ.0).AND.(NOWL.GT.0)).OR.
c     1         ((INNER.EQ.1).AND.(NIWL.EQ.0))).AND.(IVCOR.EQ.0)) THEN
c             VDU= (VPMAX(2)-EO)**PWI
c             VDL= (VPMAX(2)-GV(OWL(NOWL-1)))**PWI
c             EO= VPMAX(2)- (VDL+ (VDU - VDL)/ DBLE(KV - NF + 1))**PW
c
c** Otherwise use a binary search.
c
            ELSE
              IF((TRIALE.EQ.1).AND.(EO.LT.LHIE)) THEN
                  LHIE= EO
                  IVCOR= 0
                ELSEIF(TRIALE.EQ.1) THEN
                  LHIE= (LHIE + LLOE) / 2.0d0
                ENDIF
              TRIALE= 1
              EO= (LHIE + LLOE) / 2.0d0
            ENDIF
          IVCOR= IVCOR + 1
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** If an unknown case occurs (quite impossible but don't quote me on
c   it) then write out an error message and exit.
        ELSE
          IF(AWO.NE.0) THEN
              WRITE(6,610)
              WRITE(6,666) KV,NF
              ENDIF
          AFLAG= -7
          GOTO 200
        ENDIF
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Set KV to the next uncalculated vibrational level to be found unless
c   looking for the highest vibrational level.
c
      IF(KV.NE.999) KV= NF
c
c** If still haven't found all of the vibrational levels then
c   look for the next vibrational level.
c
      IF((KV.LE.KVMAX).OR.(KV.EQ.999)) GOTO 10
c
c** Otherwise, print out a message saying that all is well.
c
      IF((KV.GT.KVMAX).AND.(AWO.GT.0)) WRITE(6,630) KVMAX
c
c** If the potential has levels in a second minimum, then print out a
c   list of those levels to channel-21 if desired.
c
      IF((NIWL.GT.0).AND.(NOWL.GT.0)) THEN
          IF(MCO.GE.1) WRITE(21,2114) NIWL, NOWL
          IF(AWO.NE.0) THEN
              WRITE(6,605)
              WRITE(6,607)
              ENDIF
          AFLAG= -8
          ENDIF
c
  200 IF(AFLAG.LT.0) THEN
c** If unable to find all KVMAX+1 levels requested, then return KVMAX as
c  v for the highest vibrational level actually found, and print out the
c  the energy of that level.
          KVMAX= NF - 1
          IF(AWO.NE.0) WRITE(6,626) KVMAX, GV(KVMAX)
          ENDIF
      IF(MCO.GE.1) THEN
c         WRITE(21,2100)
c         DO  I= 0,NF-1
c             IF(INNR(I).EQ.0)WRITE(21,2113)I,'Outer',I,GV(I)-VMIN,RAVG(I)
c             IF(INNR(I).EQ.1)WRITE(21,2113)I,'Inner',I,GV(I)-VMIN,RAVG(I)
c             END DO
          WRITE(21,2100)
          ENDIF
      RETURN
c-----------------------------------------------------------------------
  601 FORMAT(/' Solve by matching inward and outward solutions at the'/
     1  5x,'outermost wave function maximum, S(max), where  R = RR(M)')
  603 FORMAT(/' Solve by matching inward and outward solutions at the'/
     1  5x,'innermost turning point   R1 = R(M)')
  605 FORMAT(/'  *** ALF WARNING ***')
  606 FORMAT(4x,'Next estimated trial energy  E(v=',I3,') =',G15.8/4X,
     1  'lies above potential maximum  VMAX =',G15.8)
  607 FORMAT(4X,'Potential found to have a second minimum.')
  610 FORMAT(/'  *** ALF ERROR ***')
  611 FORMAT(4X,'Attempt to find zero point level fails!')
  612 FORMAT(4X,'Attempt to find next higher vibrational level fails!')
  613 FORMAT(4X,'Number of vib levels requested=',i4,' exceeds internal 
     1ALF array dimension  NVIBMX=',i4)
  614 FORMAT(/'  *** ALF ERROR ***   Unable to find a potential minimum
     1 for   J=',i4)
  615 FORMAT(4X,'There are',I3,'  potential ',A6,' in this potential. St
     1op searching after 10.')
  616 FORMAT(4X,'The potential turns over in the short range region at R
     1 = ',G15.8)
  617 FORMAT(4X,'VMAX =',G15.8,' found to be less than VMIN =',G15.8)
  618 FORMAT(/'  ALF  finds onee potential minimum of',1PD15.7,
     1  '  at  R(1)=',0Pf9.6)
  620 FORMAT(4X,'Use of energy ',I1,'0% up the potential well (E =',
     1G15.8,')'/4X,' fails to produce a viable vibrational eigenstate.')
  621 FORMAT(4X,'Use of differences to estimate the energy for the next'
     1/4X,' vibrational level (v=',I3,') failed after',I3,'  attempt.')
  622 FORMAT(4X,'After',I3,' tries to harmonically estimate the zero-poi
     1nt energy,'/4X,' initial trial energy',G15.8,'   had yielded   E(v
     2=',I3,') =',G15.8)
  623 FORMAT(4X,'Expecting to find level (v=',I3,') but found level (v='
     1,I3,')')
  624 FORMAT(4X,'After',I3,' tries, failed to interpolate trial energy b
     1etween'/4X,'E(v=',I3,') =',G15.8,'   and   E(v=',I3,') =',G15.8)
  626 FORMAT(4X,'The highest calculated level is  E(v=',I3,') =',G15.8)
  630 FORMAT(/' ALF successfully finds all vibrational levels up to   v=
     1 KVMAX=',I3)
  640 FORMAT(/' ALF finds all  J=',i3,'  vib. levels below  vD=',F7.3,
     1  '  estimated by N-D theory')
  650 FORMAT(/' There were',I3,'  potential minima found with the absolu
     1te minimum'/4X,'VMIN =',G15.8,'  cm-1.')
  651 FORMAT(/' There were',I3,'  potential maxima found with the absolu
     1te maximum'/4X,'VMAX =',G15.8,'  cm-1.')
  666 FORMAT(4X,'Undefined case for automatic search.'/,4X,'Values of KV
     1 =',I3,'  and NF =',I3)
 2100 FORMAT(/1X,39('=='))
 2101 FORMAT(/1X,39('--'))
 2110 FORMAT(/' Limits and increment of integration (in Angstroms):'
     1 /'    RMIN =',F6.3,'    RMAX =',F7.3,'    RH =',F9.6,
     2 //' Generate    BZ =',G19.12,' ((1/cm-1)(1/Angstroms**2))'
     3 /' from ZMU:',F15.11,' (amu)')
 2111 FORMAT(/' Eigenvalue convergence criterion is   EPS =',G11.4,'(cm-
     11)')
 2112 FORMAT(/' Calculating properties of the potential described above.
     1 '/' Use Airy function at 3-rd turning point as outer boundary'
     2 /' condition for quasibound levels.')
 2113 FORMAT(' v=',I3,4X,'v(',A5,')=',I3,4X,'Gv=',F16.9,4X,'Bv=',F16.12)
 2114 FORMAT(/' Found',I4,' level(s) in the inner well and',I4,' level(s
     1) in the outer well.')
      END
c***********************************************************************
      SUBROUTINE INITVAL(IQT,IVDIF,IZPE,IEG,IVCOR)
c***********************************************************************
c** This subroutine reinitializes the condition flags when considering a
c   new case (found next vibrational level or finding level in inner
c   well - INNER = 1).
c
c** On entry and exit:
c    IQT     Case when KV < 0 and NF = 0
c            determines the value used for the initial trial energy.
c    IVDIF   Case when KV < 0 and NF > 0
c            is the flag denoting the use of differences to calculate
c            trial energies.
c    IZPE    Case when KV > 0 and NF = 0
c            is the number of times the zero point energy (v = 0) has
c            been estimated harmonically.
c    IEG     Case when KV < NF and NF > 0
c            are the number of times that a larger trial energy is used
c            to find the next level.
c    IVCOR   Case when KV > NF and NF > 0
c            are the number of times that a smaller trial energy is used
c            to find the next level.
c
      INTEGER IZPE,IVDIF,IVCOR,IQT,IEG
c
      IQT= 1
      IVDIF= 0
      IZPE= 0
      IEG= 0
      IVCOR= 0
      RETURN
      END
c***********************************************************************
      SUBROUTINE DTENG(IEG,NF,NVEL,VEL,NVIBMX,VMIN,GV,EO)
c***********************************************************************
c** This subroutine calculates the next trial energy using differences.
c
c** On entry:
c    IEG     factor by which a larger trial energy should be calculated:
c              NVEL = 2 : Increase correction by increments of 25%
c              NVEL > 2 : Increase correction by increments of 50%
c    NF      is the highest calculated vibrational level.
c    NVEL    is the number of levels found in the potential well.
c    VEL(v)  keeps track of all levels in the potential well.
c    NVIBMX  is the maximum number of vibrational levels (dimension).
c    VMIN    is the absolute value of the potential minimum (cm-1).
c    GV(v)   contains the vibrational energy level spacings
c            and rotational constants for each level (cm-1).
c
c** On exit:
c    EO      is the calculated trial energy.
c
      INTEGER IEG,NF,NVEL,NVIBMX,VEL(0:NVIBMX)
c
      REAL*8 VMIN,GV(0:NVIBMX),EO
c
c** If determining the first (non-zero point energy) level in the well,
c   then use the last determined level in the other well plus a larger
c   than harmonic correction that becomes smaller with each new
c   iteration:      E(v=0)= E(v=NF-1) + (E(v=NF-1)-VMIN)/(NF-1+IEG/4) 
c
      IF(NVEL.EQ.0) THEN
        EO= GV(NF-1) + (GV(NF-1) - VMIN)/(NF - 1 + 0.25d0*DBLE(IEG))
c
c** Try to get v= 1 using smaller-than-harmonic spacing.
c
c              E(v=1)= E(v=0) + (1.0+IEG/4)*(E(v=0)-VMIN)
c
      ELSEIF(NVEL.EQ.1) THEN
        EO= GV(VEL(0)) + (1.d0+DBLE(IEG)*0.25d0)*(GV(VEL(0))-VMIN)
c
c** Try to get v= 2 using a sequentially increasing correction.
c
c             E(v=2)= E(v=1) + (0.8+IEG/4)*(E(v=1)-E(v=0))
c
      ELSEIF(NVEL.EQ.2) THEN
        EO= GV(VEL(1))+ (0.8d0+DBLE(IEG)*0.25d0)*(GV(VEL(1))-GV(VEL(0)))
c
c** Try to get v > 2 using a sequentially increasing correction.
c
c        E(v)= E(v-1) + (1.0+IEG)*(2.0*E(v-1)-3.0*E(v-2)+E(v-3))
c
      ELSE
        EO= GV(VEL(NVEL-1)) + (1.d0+DBLE(IEG))
     1    *(2.d0*GV(VEL(NVEL-1))-3.0d0*GV(VEL(NVEL-2))+GV(VEL(NVEL-3)))
      ENDIF
      RETURN
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c**********************************************************************
      SUBROUTINE CDJOEL(EO,NBEG,NEND,BvWN,RH,WARN,V,WF0,RM2,RCNST)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c  Subroutine solving the linear inhomogeneous differential equations
c  formulated by J.M. Hutson [J.Phys.B14, 851 (1982)] for treating 
c  centrifugal distortion as a perturbation, to determine centrifugal 
c  distortion constants of a diatomic molecule.  Uses the algorithm of
c  J. Tellinghuisen [J.Mol.Spectrosc. 122, 455 (1987)].  The current
c  version calculates Bv, Dv, Hv, Lv, Mv, Nv and Ov and writes them out, 
c  but does not return values to the calling program.
c
c** On entry:   EO    is the eigenvalue (in units [cm-1])
c               NBEG & NEND  the mesh point range over which the input
c wavefunction  WF0  (in units 1/sqrt(Ang))  has non-negligible values
c               BvWn  is the numerical factor (hbar^2/2mu) [cm-1 Ang^2]
c               RH    is the integration stepsize (in units [Ang])
c               WARN  is an integer flag: > 0 print internal warnings,
c               V(i)  is the effective potential (including centrifugal
c                     term if calculation performed at  J > 0) in 
c                     'internal' units, including the factor  RH**2/BvWN
c               RM2(i) is the array  1/(distance**2) in units [1/Ang**2]
c** On exit:    RCNST(i)  is the set of 7 rotational constants: Bv, -Dv,
c                       Hv, Lv, Mv, Nv & Ov
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c                COPYRIGHT 1994  by  Robert J. Le Roy                  +
c   Dept. of Chemistry, Univ. of Waterloo, Waterloo, Ontario, Canada   +
c    This software may not be sold or any other commercial use made    +
c      of it without the express written permission of the author.     +
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c  Authors: R.J. Le Roy & J. Tellinghuisen         Version of 30/09/1999
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Dimension:  potential arrays  and  vib. level arrays.
      INTEGER NDIMR
      PARAMETER (NDIMR=100001)
      INTEGER I,M,IPASS,M1,M2,NBEG,NEND,WARN
      REAL*8 V(NEND),WF0(NEND),RM2(NEND),P(NDIMR),WF1(NDIMR),
     1                                             WF2(NDIMR),RCNST(7)
      REAL*8 BvWN,DV,DVV,HVV,HV2,LVV,LV2,MVV,MV2,NVV,OVV,EO,E,RH,RHSQ,
     1  ZTW,AR,R2IN,G2,G3,P0,P1,P2,P3,PI,PIF,PRS,PRT,V1,V2,V3,Y1,Y2,Y3,
     2  TSTHv,TSTLv,TSTMv,AMB,AMB1,AMB2,
     3  OV,OV01,OV02,OV03,OV11,OV12,OV13,OV22,OV23,OV33,
     4  PER01,PER02,PER03,PER11,PER12,PER13,PER22,PER23,PER33
c
      IF(NEND.GT.NDIMR) THEN
          WRITE(6,602) NEND,NDIMR
          RETURN
          ENDIF
      ZTW= 1.D0/12.d0
      RHSQ = RH*RH
      DV = RHSQ/12.D0
      E= EO*RHSQ/BvWN
      IPASS = 1
      OV01 = 0.D0
      OV02 = 0.D0
      OV03 = 0.D0
      OV11 = 0.D0
      OV22 = 0.D0
      OV12 = 0.D0
      OV33 = 0.D0
      OV23 = 0.D0
      OV13 = 0.D0
      PER01 = 0.D0
      PER02 = 0.D0
      PER03 = 0.D0
      PER11 = 0.D0
      PER12 = 0.D0
      PER13 = 0.D0
      PER22 = 0.D0
      PER23 = 0.D0
      PER33 = 0.D0
c** First, calculate the expectation value of  1/r**2  and hence Bv
      R2IN= 0.5D0*(RM2(NBEG)*WF0(NBEG)**2 + RM2(NEND)*WF0(NEND)**2)
      DO   I= NBEG+1, NEND-1
         R2IN= R2IN+ RM2(I)*WF0(I)**2
         ENDDO
      R2IN = R2IN*RH
      RCNST(1)= R2IN*BvWN
c
c** On First pass  IPASS=1  and calculate first-order wavefx., Dv & Hv
c  On second pass  IPASS=2  and calculate second-order wavefx., Lv & Mv
c  On third pass   IPASS=3  and calculate third-order wavefx., Nv & Ov
c
   10 P1= 0.D0
      P2= 0.D0
c
c     P1= WF0(NEND)
c     P2= WF0(NEND-1)
c
      P(NEND) = P1
      P(NEND-1) = P2
      V1 = V(NEND) - E
      V2 = V(NEND-1) - E
      IF(IPASS.EQ.1) THEN
          Y1 = P1*(1.D0 - ZTW*V1) - DV*(RM2(NEND) - R2IN)*WF0(NEND)
          G2 = (RM2(NEND-1) - R2IN)*WF0(NEND-1)
        ELSEIF(IPASS.EQ.2) THEN
          Y1 = P1*(1.D0 - ZTW*V1) - DV*((RM2(NEND) - R2IN)*WF1(NEND)
     1                                                - DVV*WF0(NEND))
          G2 = (RM2(NEND-1) - R2IN)*WF1(NEND-1) - DVV*WF0(NEND-1)
        ELSEIF(IPASS.EQ.3) THEN
          Y1 = P1*(1.D0 - ZTW*V1) - DV*((RM2(NEND) - R2IN)*WF2(NEND)
     1                                - DVV*WF1(NEND) - HVV*WF0(NEND))
          G2 = (RM2(NEND-1) - R2IN)*WF2(NEND-1) - DVV*WF1(NEND-1)
     1                                               - HVV*WF0(NEND-1)
        ENDIF
      Y2 = P2*(1.D0 - ZTW*V2) - DV*G2
      M= NEND-1
c** Now - integrate inward from outer end of range
      DO  I = NBEG+2,NEND
          M = M-1
          Y3 = Y2 + Y2 - Y1 + RHSQ*G2 + V2*P2
          IF(IPASS.EQ.1) G3 = (RM2(M) - R2IN)*WF0(M)
          IF(IPASS.EQ.2) G3 = (RM2(M) - R2IN)*WF1(M) - DVV*WF0(M)
          IF(IPASS.EQ.3) G3 = (RM2(M) - R2IN)*WF2(M) - DVV*WF1(M) 
     1                                                    - HVV*WF0(M)
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
      P2 = 0.D0
c
c     P1 = WF0(NBEG)
c     P2 = WF0(NBEG+1)
c
      P(NBEG) = P1
      P(NBEG+1) = P2
      V1 = V(NBEG) - E
      V2 = V(NBEG+1) - E
      IF(IPASS.EQ.1) THEN
          Y1 = P1*(1.D0 - ZTW*V1) - DV*(RM2(NBEG) - R2IN)*WF0(NBEG)
          G2 = (RM2(NBEG+1) - R2IN)*WF0(NBEG+1)
        ELSEIF(IPASS.EQ.2) THEN
          Y1 = P1*(1.D0 - ZTW*V1) - DV*((RM2(NBEG) - R2IN)*WF1(NBEG)
     1                                                - DVV*WF0(NEND))
          G2 = (RM2(NBEG+1) - R2IN)*WF1(NBEG+1) - DVV*WF0(NBEG+1)
        ELSEIF(IPASS.EQ.3) THEN
          Y1 = P1*(1.D0 - ZTW*V1) - DV*((RM2(NBEG) - R2IN)*WF2(NBEG)
     1                                - DVV*WF1(NEND) - HVV*WF0(NEND))
          G2 = (RM2(NBEG+1) - R2IN)*WF2(NBEG+1) - DVV*WF1(NBEG+1)
     2                                               - HVV*WF0(NBEG+1)
        ENDIF
      Y2 = P2*(1.D0 - ZTW*V2) - DV*G2
      AR = 0.D0
      M1 = M+1
c** Now ... integrate outward from inner end of range
      DO  I = NBEG+2,M1
          Y3 = Y2 + Y2 - Y1 + RHSQ*G2 + V2*P2
          P0 = WF0(I)
          IF(IPASS.EQ.1) G3 = (RM2(I) - R2IN)*P0
          IF(IPASS.EQ.2) G3 = (RM2(I)-R2IN)*WF1(I) - DVV*P0
          IF(IPASS.EQ.3) G3 = (RM2(I)-R2IN)*WF2(I) - DVV*WF1(I) - HVV*P0
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
      OV = AR*RH
      DO  I = NBEG,NEND
          P0 = WF0(I)
c ... and project out contribution of zero'th-order part of solution
          PI = P(I) - OV*P0
          PIF = PI*RM2(I)
          IF(IPASS.EQ.1) THEN
c** Now - on first pass accumulate integrals for Dv and Hv
              WF1(I) = PI
              OV01 = OV01 + PI*P0
              OV11 = OV11 + PI*PI
              PER01 = PER01 + PIF*P0
              PER11 = PER11 + PI*PIF
            ELSEIF(IPASS.EQ.2) THEN
c ... and on next pass, accumulate integrals for Lv and Mv
              WF2(I) = PI
              P1 = WF1(I)
              OV02 = OV02 + PI*P0
              OV12 = OV12 + PI*P1
              OV22 = OV22 + PI*PI
              PER02 = PER02 + PIF*P0
              PER12 = PER12 + PIF*P1
              PER22 = PER22 + PI*PIF
            ELSEIF(IPASS.EQ.3) THEN
c ... and on next pass, accumulate integrals for Nv and Ov
              P1 = WF1(I)
              P2 = WF2(I)
              OV03 = OV03 + PI*P0
              OV13 = OV13 + PI*P1
              OV23 = OV23 + PI*P2
              OV33 = OV33 + PI*PI
              PER03 = PER03 + PIF*P0
              PER13 = PER13 + PIF*P1
              PER23 = PER23 + PIF*P2
              PER33 = PER33 + PIF*PI
            ENDIF
          ENDDO
      IF(IPASS.EQ.1) THEN
          DVV = RH*PER01
          HVV = RH*(PER11 - R2IN*OV11)
          IPASS = 2
          RCNST(2) = DVV*BvWN
          RCNST(3) = HVV*BvWn
          GO TO 10
        ELSEIF(IPASS.EQ.2) THEN
          HV2 = RH*PER02*BvWN
          LVV = RH*(PER12 - R2IN*OV12 - DVV*OV11)
          MVV = RH*(PER22 - R2IN*OV22 - 2.D0*DVV*OV12 - HVV*OV11)
          IPASS = 3
          RCNST(4) = LVV*BvWN
          RCNST(5) = MVV*BvWN
          GO TO 10
        ELSEIF(IPASS.EQ.3) THEN
          LV2 = RH*PER03*BvWN
          MV2 = RH*(PER13 - R2IN*OV13 - DVV*OV12 - HVV*OV11)*BvWN
          NVV = RH*(PER23 - R2IN*OV23 - DVV*(OV13 + OV22) 
     1                                     - 2.D0*HVV*OV12 - LVV*OV11)
          OVV = RH*(PER33 - R2IN*OV33 - 2.D0*DVV*OV23 
     1             - HVV*(2.D0*OV13+ OV22) - 2.D0*LVV*OV12 - MVV*OV11)
          RCNST(6) = NVV*BvWN
          RCNST(7) = OVV*BvWN
        ENDIF
      IF(WARN.GT.0) THEN
          IF(DMAX1(DABS(OV01),DABS(OV02),DABS(OV01)).GT.1.D-9)
     1                                     WRITE(6,604) OV01,OV02,OV03
          TSTHV= dabs(RCNST(3)/HV2-1.D0)
          TSTLV= dabs(RCNST(4)/LV2-1.D0)
          TSTMV= dabs(RCNST(5)/MV2-1.D0)
          IF(DMAX1(TSTHV,TSTLV,TSTMV).GT.1.d-5)
     1                                  WRITE(6,603) TSTHV,TSTLV,TSTMV
          ENDIF
      RETURN
   90 WRITE(6,601) EO
      RETURN
  601 FORMAT(' *** ERROR in CDJOEL *** for input energy  E =',f12.4,
     1   '  never reach outer turning point')
  602 FORMAT(/' *** Dimensioning PROBLEM in CDJOEL ***   NEND=',i6,
     1  ' > NDIMR=',i6)
  603 FORMAT(' ** CAUTION ** Comparison tests for Hv, Lv & Mv give:',
     1 3(1Pd9.1))
  604 FORMAT(' ** CAUTION ** CDJOEL orthogonality tests OV01,OV02 & OV03
     1:',3(1Pd9.1))
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
c***** R.J. Le Roy  subroutine SCHRQ, last updated  8 April 2007 *******
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c                COPYRIGHT 2007  by  Robert J. Le Roy                  +
c   Dept. of Chemistry, Univ. of Waterloo, Waterloo, Ontario, Canada   +
c    This software may not be sold or any other commercial use made    +
c      of it without the express written permission of the author.     +
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** SCHRQ solves radial Schrodinger equation in dimensionless form
c  d2WF/dR2 = - (E-V(R))*WF(R) ,  where WF(I) is the wave function.
c** Integrate by Numerov method over N mesh points with increment
c  H=RH across range beginning at RMIN .
c** Input trial energy EO, eigenvalue convergence criterion EEPS
c  potential asymptote VLIM, and all returned energies (EO, GAMA & VMAX)
c  have units (cm-1).
c** On entry, the input potential V(I) must include the centrifugal
c  term and the factor:  'BFCT'=2*mu*(2*pi*RH/hPLANCK)**2  (1/cm-1) ,
c  which is also internally incorporated into EO, VLIM & EEPS.
c* Note that these reduced quantities (& the internal eigenvalue E)
c  contain a factor of the squared integration increment  RH**2 .
c  This saves arithmetic work in the innermost loop of the algorithm.
c** For energy in (cm-1), BFCT=ZMU(u)*H(Angst)**2/16.85762920 (1/cm-1)
c** INNODE > 0  specifies that wavefx. initiates at RMIN with a node 
c     (normal default case);  INNODE.le.0  specifies  zero slope  at
c     RMIN (for finding symmetric eigenfunctions of symmetric potential
c     with potential mid-point @ RMIN).
c** INNER specifies wave function matching condition: INNER = 0  makes
c     matching of inward & outward solutions occur at outermost turning
c     point;  INNER > 0 makes matching occur at innermost turning point.
c * Normally use  INNER=0 ,  but to find inner-well levels of double 
c     minimum potential, set  INNER > 0 .
c----------------------------------------------------------------------
      SUBROUTINE SCHRQ(KV,JROT,EO,GAMA,VMAX,VLIM,V,WF,BFCT,EEPS,RMIN,
     1                          RH,N,NBEG,NEND,INNODE,INNER,IWR,LPRWF)
c----------------------------------------------------------------------
c** Output vibrational quantum number KV, eigenvalue EO, normalized
c  wave function WF(I), and range, NBEG .le. I .le. NEND  over
c  which WF(I) is defined. *** Have set  WF(I)=0  outside this range.
c* (NBEG,NEND), defined by requiring  abs(WF(I)) < RATST=1.D-9  outside.
c** If(LPRWF.gt.0) print wavefunction WF(I) every LPRWF-th point.
c* If(LPRWF.lt.0) "punch" (i.e., WRITE(10,XXX)) every |LPRWF|-th point
c  of the wave function on disk starting at R(NBEG) with step size
c  of  IPSIQ=|LPRWF|*RH. 
c** For energies above the potential asymptote VLIM, locate quasibound
c  levels using Airy function boundary condition and return the level
c  width GAMA and barrier height VMAX, as well as EO.
c** ERROR condition on return is  KV < 0 ; usually KV=-1, but return
c  KV=-2 if error appears to arise from too low trial energy.
c** If(IWR.ne.0) print error & warning descriptions
c  If (IWR.gt.0) also print final eigenvalues & node count.
c  If (IWR.ge.2) also show end-of-range wave function amplitudes
c  If (IWR.ge.3) print also intermediate trial eigenvalues, etc.
c** If input KV.ge.998 , tries to find highest bound level, and
c  trial energy should be only slightly less than VLIM.
c** If input KV < -10 , use log-derivative outer boundary condition at
c  mesh point |KV| , based on incoming value of wave function WF(|KV|)
c  and of the wavefunction derivative at that point, SPNEND, which is
c  brought in as WF(|KV|-1).  For a hard wall condition at mesh point
c  |KV|, set WF(|KV|)=0 and WF(|KV|-1)= -1 before entry.
c----------------------------------------------------------------------
c++ "SCHRQ" calls subroutineas "QBOUND" and "WIDTH", and the latter
c++ calls "LEVQAD" .
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      INTEGER  I,IBEGIN,ICOR,IJ,IJK,INNODE,INNER,IPSID,IQTST,IT,
     1         ITER,ITP1,ITP1P,ITP3,IWR,J,JJ,J1,J2,JPSIQ,JQTST,JROT,
     2         KKV,KV,KVIN,LPRWF,M,MS,MSAVE,
     3         N,NBEG,NBEGB,NBEG2,NDN,NEND,NENDCH,NLINES,NPR
      REAL*8  BFCT,DE,DEP,DEPRN,DF,DOLD,DSOC,
     2        E,EEPS,EO,EPS,F,FX,GAMA,GI,GN,H,H2,HT,PROD,PPROD,
     3        RATIN,RATOUT,RATST,RH,RINC,RMIN,RMINN,RR,RSTT,RWR(20),
     4        WF(N),SB,SI,SM,SN,SNEND,SPNEND,SRTGI,SRTGN,SWR(20),
     5        V(N),VLIM,VMAX,VMX,VPR,
     6        WKBTST,XEND,XPR,XPW,DXPW,Y1,Y2,Y3,YIN,YM,YOUT
      DATA RATST/1.D-9/,XPW/20.72d0/
      DATA NDN/15/
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      DXPW= (XPW+ 2.30d0)/NDN
      ICOR= 0
      KVIN= KV
      KV= -1
      RMINN= RMIN-RH
      GAMA= 0.d0
      VMAX= VLIM
      VMX= VMAX*BFCT
      H= RH
      H2= H*H
      HT= 1.d0/12.D+0
      E= EO*BFCT
      EPS= EEPS*BFCT
      DSOC= VLIM*BFCT
      DE= 0.d0
      RATIN= 0.d0
      RATOUT= 0.d0
      IF(IWR.GT.2) THEN
          IF(KVIN.GE.998) then
              WRITE(6,610) EO
            ELSE
              WRITE(6,601) KVIN,JROT,EO,INNER
            ENDIF
          WRITE(6,602)
        ENDIF
      NEND= N
      IF(KVIN.LT.-10) THEN
          NEND= -KVIN
          SNEND= WF(NEND)
          SPNEND= WF(NEND-1)
          ENDIF
      JQTST = 0
c** Start iterative loop; try to converge for up to 15 iterations.
      DO 90 IT= 1,15
          ITER= IT
          IF(INNER.GT.0) GO TO 38
   10     IF(KVIN.LT.-10) THEN
c** If desired, (KVIN < -10) outer boundary set at NEND=|KVIN| and 
c  initialize wavefunction with log-derivative condition based on value
c  WF(NEND) & derivative SPNEND at that mesh point (brought in in CALL)
              GN= V(NEND)-E
              GI= V(NEND-1)-E
              SB= SNEND
              SI= SB*(1.d0+ 0.5d0*GN)- RH*SPNEND
              GO TO 24
              END IF
          IF(E.GE.DSOC) THEN
c** For quasibound levels, initialize wave function in "QBOUND"
              CALL QBOUND(KVIN,JROT,E,EO,VMX,DSOC,V,RMIN,H,GN,GI,
     1                                 SB,SI,N,ITP3,IWR,IQTST,BFCT,IT)
              NEND= ITP3
              VMAX= VMX/BFCT
              IF(IQTST.GT.0) GO TO 24
              IF(IQTST.LT.0) THEN
                  JQTST = JQTST+IQTST
                  IF((JQTST.LE.-2).OR.(VMAX.LT.VLIM)) GO TO 999
c** Try up to once to find level using trial value just below maximum
                  EO = VMAX-0.1D0
                  E = EO*BFCT
                  GO TO 90
                  ENDIF
              GO TO 20
              ENDIF
c** For  E < DSOC  begin inward integration by using JWKB to estimate
c  optimum (minimum) inward starting point which will still give 
c  RATOUT < RATST = exp(-XPW) (ca. 1.d-9) [not needed after 1'st 2 ITER]
          IF(ITER.LE.2) THEN
              NEND= N
c ... first do rough inward search for outermost turning point
              DO  M= N,1,-NDN
                  MS= M
                  GI= V(M)- E
                  IF(GI.LE.0.D0) GO TO 12
                  GN= GI
                  ENDDO
              IF(IWR.NE.0) WRITE(6,611) JROT,EO
              GO TO 999
   12         IF(MS.GE.N) GO TO 998
              FX= GN/(GI-GN)
              SM= 0.5d0*(1.d0+ FX)*DSQRT(GN)
              MS= MS+ 2*NDN
              IF(MS.GE.N) GO TO 20
c ... now integrate exponent till JWKB wave fx. would be negligible
              DO  M= MS,N,NDN
                  NEND= M
                  SM= SM+ DSQRT(V(M)- E)
                  IF(SM.GT.DXPW) GO TO 18
                  ENDDO
   18         IF(NEND.LT.N) NEND= NEND+ NDN
              ENDIF
c** For truly bound state initialize wave function as 1-st order WKB
c   solution increasing inward
   20     GN= V(NEND)- E
          GI= V(NEND-1)- E
          MS= NEND-1
          IF(GI.LT.0.d0) GO TO 998
          SRTGN= DSQRT(GN)
          SRTGI= DSQRT(GI)
          SB= 1.d0
          SI= SB*DSQRT(SRTGN/SRTGI)*DEXP((SRTGN+SRTGI)*0.5d0)
          IF(SB.GT.SI) THEN
c WOOPS - JWKB gives inward DEcreasing solution, so initialize with node
              IF(IWR.NE.0) WRITE(6,618) JROT,EO,SB/SI
              SI= 1.d0
              SB= 0.d0
              ENDIF
   24     M= NEND-1
          Y1= (1.d0-HT*GN)*SB
          Y2= (1.d0-HT*GI)*SI
          WF(NEND)= SB
          WF(NEND-1)= SI
          MS= NEND
          NENDCH= NEND
          IBEGIN= 3
          IF(INNER.GT.0) IBEGIN= ITP1+2
c** Actual inward integration loop starts here
          DO  I= IBEGIN,NEND
              M= M-1
              Y3= Y2+Y2-Y1+GI*SI
              GI= V(M)-E
              SB= SI
              SI= Y3/(1.d0-HT*GI)
              WF(M)= SI
              IF(DABS(SI).GE.1.D+17) THEN
c** Renormalize to prevent overflow of  WF(I)  in classically
c  forbidden region where  (V(I) .gt. E)
                  SI= 1.d0/SI
                  DO  J= M,MS
                      WF(J)= WF(J)*SI
                      ENDDO
                  NENDCH= MS
                  MS= M
                  Y2= Y2*SI
                  Y3= Y3*SI
                  SB= SB*SI
                  SI= 1.d0
                  ENDIF
              Y1= Y2
              Y2= Y3
c** Test for outermost maximum of wave function.
c ... old matching condition - turning point works OK & is simpler.
cc            IF((INNER.EQ.0).AND.(SI.LE.SB)) GO TO 32
c** Test for outer well turning point 
              IF((INNER.EQ.0).AND.(GI.lt.0.d0)) GO TO 32
              ENDDO
          IF(INNER.EQ.0) THEN
c** Error mode ... inward propagation finds no turning point
              KV= -2
              IF(IWR.NE.0) WRITE(6,616) KV,JROT,EO
              GO TO 999
              ENDIF
c** Scale outer part of wave function before proceding
   32     SI= 1.d0/SI
          MSAVE= M
          RR= RMINN+MSAVE*H
          YIN= Y1*SI
          RATOUT= WF(NEND)*SI
          NEND= NENDCH
          DO  J= MSAVE,NEND
              WF(J)= WF(J)*SI
              ENDDO
          IF(INNER.GT.0) GO TO 70
c-------------------------------------------------------------------
c** Set up to prepare for outward integration **********************
   38     NBEG= 1
          IF(INNODE.LE.0) THEN
c** Option to initialize with zero slope at beginning of the range
              SB= 1.d0
              GN= V(1)-E
              Y1= SB*(1.d0-HT*GN)
              Y2= Y1+GN*SB*0.5d0
              GI= V(2)-E
              SI= Y2/(1.d0-HT*GI)
            ELSE
c** Initialize outward integration with a node at beginning of range
   40         GN= V(NBEG)-E
              IF(GN.GT.10.D0) THEN
c** If potential has [V(1)-E] so high that H is (locally) much too
c  large, then shift inner starting point outward.
                  NBEG= NBEG+1
                  IF(NBEG.LT.N) GO TO 40
                  IF(IWR.NE.0) WRITE(6,613)
                  GO TO 999
                  ENDIF
              IF((ITER.LE.1).AND.(IWR.NE.0)) THEN
                  IF(NBEG.GT.1) WRITE(6,609) JROT,EO,NBEG
                  IF(GN.LE.0.d0) WRITE(6,604) JROT,EO,NBEG,V(NBEG)/BFCT
                  ENDIF
c** Initialize outward wave function with a node:  WF(NBEG) = 0.
              SB= 0.d0
              SI= 1.d0
              GI= V(NBEG+1)-E
              Y1= SB*(1.d0- HT*GN)
              Y2= SI*(1.d0- HT*GI)
            ENDIF
c
          WF(NBEG)= SB
          WF(NBEG+1)= SI
          NBEGB= NBEG
          NBEG2= NBEG+2
          IF(INNER.GT.0) MSAVE= N
c** Actual outward integration loops start here
          DO  I= NBEG2,MSAVE
              Y3= Y2+Y2-Y1+GI*SI
              GI= V(I)-E
              SI= Y3/(1.d0- HT*GI)
              WF(I)= SI
              IF(DABS(SI).GE.1.D+17) THEN
c** Renormalize to prevent overflow of  WF(I)  in classically forbidden
c  region where  V(I) .gt. E
                  SI= 1.d0/SI
                  NBEG= NBEGB
                  DO  J= NBEG,I
                      WF(J)= WF(J)*SI
                      ENDDO
                  NBEGB= I
                  Y2= Y2*SI
                  Y3= Y3*SI
                  SI= 1.d0
                  ENDIF
              Y1= Y2
              Y2= Y3
              ITP1= I
c** Exit from this loop at onset of classically allowed region
              IF(GI.LE.0.d0) GO TO 52
              ENDDO
          MS= MSAVE
          IF((INNER.EQ.0).AND.(GN.LE.0.d0)) GO TO 60
          IF(IWR.NE.0) WRITE(6,612) KVIN,JROT,EO,MSAVE
          GO TO 999
   52     ITP1P= ITP1+1
          MS= ITP1
          IF(INNER.GT.0) GO TO 60
          DO  I= ITP1P,MSAVE
              Y3= Y2+Y2-Y1+GI*SI
              GI= V(I)-E
              SI= Y3/(1.d0- HT*GI)
              WF(I)= SI
              IF(DABS(SI).GT.1.D+17) THEN
c** Renormalize to prevent overflow of  WF(I) , as needed.
                  SI= 1.d0/SI
                  NBEG= NBEGB
                  DO  J= NBEG,I
                      WF(J)= WF(J)*SI
                      ENDDO
                  NBEGB= I
                  Y2= Y2*SI
                  Y3= Y3*SI
                  SI= 1.d0
                  ENDIF
              Y1= Y2
              Y2= Y3
              ENDDO
          MS= MSAVE
c** Finished outward integration.  Normalize w.r.t. WF(MSAVE)
   60     SI= 1.d0/SI
          YOUT= Y1*SI
          YM= Y2*SI
          RATIN= WF(NBEG+1)*SI
          DO  I= NBEG,MS
              WF(I)= WF(I)*SI
              ENDDO
          IF(INNER.GT.0) GO TO 10
c----- Finished numerical integration ... now correct trial energy
c** DF*H  is the integral of  (WF(I))**2 dR
   70     DF= 0.d0
          DO  J= NBEG,NEND
              DF= DF+WF(J)**2
              ENDDO
c** Add edge correction to DF assuming wave function dies off as simple
c  exponential past R(NEND);  matters only if WF(NEND) unusually large.
          IF((E.LE.DSOC).AND.(WF(NEND).NE.0)) THEN
              IF((KVIN.GE.-10).AND.(WF(NEND-1)/WF(NEND).GT.1.d0))
     1              DF= DF+ WF(NEND)**2/(2.d0*DLOG(WF(NEND-1)/WF(NEND)))
              ENDIF
          F= (-YOUT-YIN+2.d0*YM+GI)
          DOLD= DE
          IF(DABS(F).LE.1.D+30) THEN
              DE= F/DF
            ELSE
              F= 9.9D+30
              DF= F
              DE= DABS(0.01D+0 *(DSOC-E))
            ENDIF
          IF(IWR.GT.2) THEN
              DEPRN = DE/BFCT
              XEND= RMINN+NEND*H
c** RATIN & RATOUT  are wave fx. amplitude at inner/outer ends of range
c  relative to its value at outermost extremum.
              WRITE(6,603) IT,EO,F,DF,DEPRN,MSAVE,RR,RATIN,RATOUT,
     1                                                  XEND,NBEG,ITP1
              ENDIF
c** Test trial eigenvalue for convergence
          IF(DABS(DE).LE.DABS(EPS)) GO TO 100
          E= E+DE
c** KV.ge.998  Option ... Search for highest bound level.  Adjust new
c  trial energy downward if it would have been above dissociation.
          IF((KVIN.GE.998).AND.(E.GT.VMX)) E= VMX- 2.d0*(VMX-E+DE)
          EO= E/BFCT
          IF((IT.GT.4).AND.(DABS(DE).GE.DABS(DOLD)).AND.
     1                                       ((DOLD*DE).LE.0.d0)) THEN
c** Adjust energy increment if having convergence difficulties.  Not
c  usually needed except for some quasibounds extremely near  VMAX .
              ICOR= ICOR+1
              DEP= DE/BFCT
              IF(IWR.NE.0) WRITE(6,617) IT,DEP
              DE= 0.5d0*DE
              E= E-DE
              EO= E/BFCT
              ENDIF
   90     CONTINUE
c** End of iterative loop which searches for eigenvalue ************
c-------------------------------------------------------------------*
c** Convergence fails, so return in error condition
      E= E-DE
      EO= E/BFCT
      DEPRN= DE/BFCT
      IF(IWR.NE.0) WRITE(6,620) KVIN,JROT,ITER,DEPRN
      GO TO 999
  100 IF(IWR.NE.0) THEN
          IF(IWR.GE.3) WRITE(6,619)
          IF((DABS(RATIN).GT.RATST).AND.(INNODE.GT.0)
     1                  .AND.(RMIN.GT.0.d0)) WRITE(6,614) JROT,EO,RATIN
          IF((E.LT.DSOC).AND.(DABS(RATOUT).GT.RATST)) THEN
              WKBTST=0.5d0*DABS(V(NEND)-V(NEND-1))/DSQRT((V(NEND)-E)**3)
              IF(WKBTST.GT.1.d-3)WRITE(6,615)JROT,EO,RATOUT,RATST,WKBTST
              ENDIF
          ENDIF
      KKV = 0
c** Perform node count on converged solution
      PROD= WF(ITP1)*WF(ITP1-1)
      J1= ITP1+1
      J2= NEND-1
      DO  J= J1, J2
          PPROD= PROD
          PROD= WF(J)*WF(J-1)
          IF((PPROD.LE.0.d0).AND.(PROD.GT.0.d0)) KKV= KKV+1
          ENDDO
      KV = KKV
c** Normalize & find interval (NBEG,NEND) where WF(I) is non-negligible
      SN= 1.d0/DSQRT(H*DF)
      DO  I= NBEG,NEND
          WF(I)= WF(I)*SN
          ENDDO
      IF(ITP1.LE.1) GO TO 122
      J= ITP1P
      DO  I= 1,ITP1
          J= J-1
          IF(DABS(WF(J)).LT.RATST) GO TO 119
          ENDDO
  119 NBEG= J
      IF(NBEG.LE.1) GO TO 122
      J= J-1
      DO  I= 1,J
          WF(I)= 0.d0
          ENDDO
  122 IF(KVIN.GE.-10) THEN
c** For "non-wall" cases, move NEND inward to where wavefunction 
c  "non-negligible"
          J= NEND-1
          DO  I= NBEG,NEND
              IF(DABS(WF(J)).GT.RATST) GO TO 126
              J= J-1
              ENDDO
  126     NEND= J+1
          END IF
      IF(NEND.LT.N) THEN
c** Zero out wavefunction array at distances past NEND
          DO  I= NEND+1,N
              WF(I)= 0.d0
              ENDDO
          ENDIF
      IF(LPRWF.LT.0) THEN
c** If desired, write every |LPRWF|-th point of the wave function 
c  to a file on channel-10, starting at the NBEG-th mesh point.
          JPSIQ= -LPRWF
          NPR= 1+(NEND-NBEG)/JPSIQ
          RINC= RH*JPSIQ
          RSTT= RMINN+NBEG*RH
c** Write every JPSIQ-th point of the wave function for level  v=KV
c  J=JROT , beginning at mesh point NBEG & distance RSTT where
c  the NPR values written separated by mesh step RINC=JPSIQ*RH
          WRITE(10,701) KV,JROT,EO,NPR,RSTT,RINC,NBEG,JPSIQ
          WRITE(10,702) (RMINN+I*RH,WF(I),I=NBEG,NEND,JPSIQ)
          GO TO 140
          ENDIF
c** Print solutions every  LPRWF-th  point, 6 to a line, in columns.
      IF(LPRWF.GT.0) THEN
          NLINES= ((1+(NEND-NBEG)/LPRWF)+3)/4
          IPSID= LPRWF*NLINES
          WRITE(6,605) KV,JROT,EO
          DO  J= 1,NLINES
              JJ= NBEG+(J-1)*LPRWF
              IJK= 0
              DO  IJ= JJ,NEND,IPSID
                  IJK= IJK+1
                  RWR(IJK)= RMINN+IJ*H
                  SWR(IJK)= WF(IJ)
                  ENDDO
              WRITE(6,606) (RWR(I),SWR(I),I= 1,IJK)
              ENDDO
          ENDIF
  140 IF(IWR.EQ.1) WRITE(6,607) KV,JROT,EO
      IF(IWR.GE.2) WRITE(6,607) KV,JROT,EO,ITER,RR,RATIN,RATOUT
c** For quasibound levels, calculate width in subroutine "WIDTH"
      IF((E.GT.DSOC).AND.(KVIN.GT.-10)) CALL WIDTH(KV,JROT,E,EO,DSOC,
     1  V,WF,VMX,RMIN,H,BFCT,IWR,ITP1,ITP3,INNER,N,GAMA)
      RETURN
c** ERROR condition if  E.gt.V(R)  at outer end of integration range.
  998 XPR= RMINN+MS*H
      VPR= V(MS)/BFCT
      IF(IWR.NE.0) WRITE(6,608) EO,MS,VPR,XPR,IT
c** Return in error mode
  999 KV= -1
      RETURN
  601 FORMAT(/' Solve for  v=',I3,'   J=',I3,'   ETRIAL=',1PD15.7,
     1   '  INNER=',i2,'   WF(1st) WF(NEND)' )
  602 FORMAT(' ITER    ETRIAL',8X,'F(E)      DF(E)     D(E)',
     1 5X,'M    R(M)  /WF(M)   /WF(M)  R(NEND) NBEG ITP1'/
     2  1X,96('-'))
  603 FORMAT(I4,1PD15.7,3D10.2,0PI5,F7.3,1P2D9.1,0PF8.2,I4,I5)
  604 FORMAT('   NOTE:  for  J=',I3,'   EO=',F12.4,' .ge. V(',i3,')=',
     1  F12.4)
  605 FORMAT(/' Solution of radial Schr. equation for   E(v=',I3,',J=',
     1  I3,') =',F15.7/2x,4('    R(I)   WF(I)   ')/2X,38('--') )
  606 FORMAT(2X,4(F8.3,F11.7))
  607 FORMAT('E(v=',I3,',J=',I3,')=',F11.4,1x,I3,' Iterations',
     1  '   R(M)=',F6.3,'  WF(NBEG)/WF(M)=',1PD8.1/
     2  57x,'WF(NEND)/WF(M)=',D8.1)
  608 FORMAT(' *** SCHRQ Error:  E=',F9.2,' > V(',I5,')=',F9.2,
     1  '  at  Rmax=',F6.2,'  for  IT=',I2)
  609 FORMAT(' *** For  J=',I3,'   E=',1PD15.7,"  integration can't",
     1 ' start till past mesh'/37x,'point',I5,',  so RMIN smaller than n
     2eeded')
  610 FORMAT(/' Attempt to find the highest bound level starting from',
     1 '   ETRIAL =',1PD9.2)
  611 FORMAT(' *** SCHRQ inward search at   J=',i3,'   E=',f11.2,
     1  ' finds no classical region')
  612 FORMAT(/' *** ERROR *** for   v =',I3,'   J =',I3,'   E =',
     1  F12.4,'  Innermost turning point not found by   M = MSAVE =',I5)
  613 FORMAT(/' *** ERROR in potential array ... V(I) everywhere',
     1 ' too big to integrate with given  increment')
  614 FORMAT(' *** CAUTION *** For  J=',I3,'  E=',G15.8/16x,
     1 'WF(first)/WF(Max)=',D9.2,'  suggests  RMIN  may be too large')
  615 FORMAT(' ** CAUTION ** For  J=',I3,'  E=',1PD13.6,
     1 '  WF(NEND)/WF(Max)=',D8.1,' >',D8.1/4X,'& initialization ',
     2 'quality test ',1PD8.1,' > 1.D-3   so RMAX may be too small')
  616 FORMAT(' ** WARNING *** For  v=',I2,', J=',I3,' at  E=',G14.7,
     1  ':  inward propagation finds no turning point ... Energy too low
     2 or potential too weak' )
  617 FORMAT(' *** SCHRQ has a convergence problem, so for  IT=',I2,
     1 '  cut  DE=',1PD10.2,'  in HALF' )
  618 FORMAT(' *** For  J=',I3,'  E=',F9.2,'  JWKB start gives  SB/SI=',
     1  1PD10.3,'  so use a node.')
  619 FORMAT(1X,96('-'))
  620 FORMAT(' *** CAUTION for  v=',I3,'  J=',I3,"  SCHRQ doesn't conver
     1ge by  ITER=",I2,'  DE=',1PD9.2)
  701 FORMAT(/2x,'Level  v=',I3,'   J=',I3,'   E=',F12.4,' ,  wave funct
     1ion at',I6,' points.'/7x,'R(1-st)=',F12.8,'   mesh=',F12.8,
     2  '   NBEG=',I4,'   |LPRWF|=',I3)
  702 FORMAT((1X,4(f9.4,f10.6)))
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c*******************************************************************
      SUBROUTINE QBOUND(KV,JROT,E,EO,VMX,DSOC,V,RMIN,H,GB,GI,SB,SI,N,
     1  ITP3,IWR,IQTST,BFCT,IT)
c*******************************************************************
c** Subroutine to initialize quasibound level wave function as Airy
c  function at third turning point (if possible). For the theory see 
c  J.Chem.Phys. 54, 5114 (1971),  J.Chem.Phys. 69, 3622-31 (1978) 
c----------------------------------------------------------------------
c** IQTST  is error flag. *** If (IQTST.lt.0) initialization fails
c  so eigenvalue calculation aborts *** (IQTST.gt.0) for successful
c  Airy function initialization. *** (IQTST=0) if Airy function
c  initialization prevented because 3-rd turning point beyond
c  range, so that WKB initialization is used.
c----------------------------------------------------------------------
      INTEGER I,II,IQTST,IT,ITP3,IWR,J,JROT,K,KV,N
      REAL*8  A1,A2,A13,A23,BFCT,
     1        C1A,C2A,DF,DSOC,E,EO,FBA,FIA,FJ,GB,GBA,GI,GIA,H,
     2        RMIN,RMINN,SB,SI,SL,V(N),VMX,VMXPR,XJ1
      DATA C1A/0.355028053887817D0/,C2A/0.258819403792807D0/
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IQTST=1
      RMINN=RMIN-H
c** Start by searching for third turning point.
      J=N
      IF(V(N).GT.E) GO TO 22
      DO  I=2,N
          J=J-1
          IF(V(J).GT.E) GO TO 10
          ENDDO
      GO TO 14
c** Check that there is a classically allowed region inside this point
c  and determine height of barrier maximum.
   10 II=J
      VMX=DSOC
      DO  I=2,J
          II=II-1
          IF(V(II).LE.E) GO TO 16
          IF(V(II).GT.VMX) VMX=V(II)
          ENDDO
c** Energy too high ... find no more than one turning point.
   14 XJ1=RMINN+J*H
c ... Search outward for barrier height to facilitate energy correction
      IF(J.EQ.1) J= 2
      K=J-1
      DO  I=J,N
          IF(V(I).GT.V(K)) GO TO 120
          K=I
          ENDDO
      VMX=V(K)
      GO TO 130
  120 K=K+2
      J=K-1
      DO  I=K,N
          IF(V(I).LT.V(J)) GO TO 126
          J=I
          ENDDO
  126 VMX=V(J)
  130 VMXPR=VMX/BFCT
      IF(IWR.NE.0) WRITE(6,608) JROT,EO,VMXPR,XJ1
      ITP3= J
      IQTST=-1
      GO TO 100
   16 ITP3= J+1
c** ITP3 is the first mesh point outside classically forbidden region
      GB=V(ITP3)-E
      GI=V(ITP3-1)-E
      FJ=GI/(GI-GB)
c** Treat quasibound levels as bound using outer boundary condition
c  of Airy function at third turning point ... as discussed by
c  R.J.Le Roy and R.B.Bernstein  in  J.Chem.Phys. 54,5114(1971).
      SL=(GI-GB)**(1.d0/3.d0)/H
      IF((SL*H).LT.1.d0) THEN
          A1=GI/(SL*H)**2
          A2=GB/(SL*H)**2
          A13=A1*A1*A1
          A23=A2*A2*A2
          FIA= 1.d0+ A13*(A13*(A13+72.D0)+2160.D0)/12960.D0
          GIA=A1+A1*A13*(A13*(A13+90.D0)+3780.D0)/45360.D0
          FBA= 1.d0+ A23*(A23*(A23+72.D0)+2160.D0)/12960.D0
          GBA=A2+A2*A23*(A23*(A23+90.D0)+3780.D0)/45360.D0
c** Airy function  Bi(X)  at points straddling 3-rd turning point
          SI=C1A*FIA+C2A*GIA
          SB=C1A*FBA+C2A*GBA
          GO TO 100
          ENDIF
c** If Airy function expansion unreliable, use zero slope at third
c  turning point as quasibound outer boundary condition.
      DF=GI-GB
      SI= 1.d0+ DF*FJ**3/6.d0
      SB= 1.d0 -DF*(1.d0- FJ)**3/6.d0
      IF(IWR.NE.0) WRITE(6,606) KV,JROT,EO,IT
      GO TO 100
c** If 3-rd turning point beyond range start with WKB wave function
c  at end of range.
   22 IF(IWR.NE.0) WRITE(6,607) JROT,EO
      ITP3= N
      IQTST=0
      GB=V(ITP3)-E
      GI=V(ITP3-1)-E
      VMX=V(ITP3)
      II=ITP3
      DO  I=2,ITP3
          II=II-1
          IF(V(II).LT.VMX) GO TO 100
          VMX=V(II)
          ENDDO
      IF(IWR.NE.0) WRITE(6,604)
c** End of quasibound level initialization schemes.
      IQTST=-9
  100 RETURN
  604 FORMAT(" **** QBOUND doesn't work ... no classically allowed regio
     1n accessible at this energy.")
  606 FORMAT(' *** CAUTION ***  v=',I3,'   J=',I3,'   E=',1PD13.6,
     1 '   IT=',I2/5x,'Airy initialization unstable so use  zero slope',
     2 'at  R(3-rd)' )
  607 FORMAT(' *** For  J=',I3,'  E=',F9.2,
     1  '  R(3-rd) > RMAX  & E < V(N)  so try WKB B.C. @ RMAX')
  608 FORMAT(' For J=',I3,'  ETRY=',F11.4,' > VMAX=',F11.4,
     1  '  find onee turn point:  R=',F6.2)
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
c** Subroutine to calculates quasibound level tunneling lifetime/width
c** For relevant theory see Le Roy & Liu [J.Chem.Phys.69,3622-31(1978)]
c  and Connor & Smith [Mol.Phys. 43, 397 (1981)] and Huang & Le Roy 
c  [J.Chem.Phys. 119, 7398 (2003); Erratum, ibid, 127, xxxx (2007)]
c** Final level width calculation from Eq.(4.5) of Connor & Smith.
c  Rearranged slightly for consistency with PotFit derivatives 9/05/02
c-----------------------------------------------------------------------
      SUBROUTINE WIDTH(KV,JROT,E,EO,DSOC,V,S,VMX,RMIN,H,BFCT,IWR,ITP1,
     1  ITP3,INNER,N,GAMA)
c++ "WIDTH" calls subroutine "LEVQAD" ++++++++++++++++++++++++++++++++++
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      INTEGER  I,IMM,INNER,IRM,ITP1,ITP1P,ITP1P1,ITP2,ITP2M,ITP2M2,
     1         ITP2P1,ITP2P2,ITP3,IWR,JROT,KV,KVI,KVO,
     2         M,M2,N,NN,NST
      REAL*8  ANS1,ANS2,ARG,BFCT,COR,
     1        D1,D2,D3,DFI,DSGB,DSGN,DSOC,DWEB,OMEGJC,
     2        E,EO,EMSC,EMV,G1,G2,G3,GA,GAMA,GAMALG,
     3        H,H2,HBW,HBWB,PI,PMX,RMIN,RMINN,RMX,RT,RT1,RT2,
     4        S(N),SM,TAU,TAULG,TI,TUN0,U1,U2,V(N),VMAX,VMX,
     7        XJ,XX
      CHARACTER*5 LWELL(2)
      DATA PI/3.141592653589793D0/
      DATA LWELL/'INNER','OUTER'/
      RMINN= RMIN- H
      H2= H*H
c** ITP1 is first mesh point to right of innermost turning point.
   40 ITP1P= ITP1+ 1
      ITP1P1= ITP1P+ 1
      IRM= ITP1- 1
c** Calculate JWKB tunneling probability from quadrature over barrier
c** First must locate 2-nd turning point.
      DO  I= ITP1P1,ITP3
          ITP2= I
          IF(V(I).GT.E) GO TO 202
          ENDDO
      GAMA= 0.d0
      GO TO 250
  202 ITP2P1= ITP2+ 1
      ITP2P2= ITP2+ 2
c** ITP2M is the last mesh point before the 2-nd turning point.
      ITP2M= ITP2- 1
      ITP2M2= ITP2- 2
      G1= V(ITP2M)- E
      G2= V(ITP2)- E
      GA= V(ITP2P1)- E
c** Quadrature over barrier starts here.
      CALL LEVQAD(G1,G2,GA,H,RT,ANS1,ANS2)
      SM= ANS2/H
      IF(GA.LT.0.d0) GO TO 218
      SM= SM+ 0.5d0*DSQRT(GA)
      PMX= VMX
      M2= ITP2P2
  204 DO  I=M2,ITP3
          M= I
          GA= V(I)- E
          IF(V(I).GT.PMX) PMX=V(I)
          IF(GA.LT.0.d0) GO TO 210
          SM= SM+ DSQRT(GA)
          ENDDO
      IF(V(M).GT.V(M-1)) THEN
          IF(IWR.NE.0) WRITE(6,602) KV,JROT
          GO TO 250
          ENDIF
      RMX= RMINN+ M*H
      U1= DSQRT(GA/(V(M)- DSOC))
      U2= DSQRT((E- DSOC)/(V(M)- DSOC))
      SM= SM- 0.5d0*DSQRT(GA)+ (DLOG((1.d0+U1)/U2)-U1)*RMX*
     1                                             DSQRT(V(M)- DSOC)/H
      XJ= (DSQRT(1.d0+ 4.d0*(V(M)-DSOC)*(RMX/H)**2)- 1.d0)*0.5d0
      IF(IWR.NE.0) WRITE(6,603) JROT,EO,XJ,RMX
      GO TO 218
  210 IF(M.LT.ITP3) THEN
c** If encounter a double-humped barrier, take care here.
          IF(IWR.NE.0) WRITE(6,609) KV,JROT,EO,M
          KVO= 0
          DSGN= DSIGN(1.d0,S(M-1))
c** Find the effective quantum number for the outer well
          DO  I= M,ITP3
              DSGB= DSGN
              DSGN= DSIGN(1.d0,S(I))
              IF((DSGN*DSGB).LT.0.d0) KVO=KVO+1
              ENDDO
          KVI= KV- KVO
          IF(INNER.EQ.0) THEN
c** For levels of outer well, get correct width by changing ITP1
              ITP1= M
              IF(IWR.GT.0) WRITE(6,610) KVO,LWELL(2)
              GO TO 40
              ENDIF
          IF(IWR.GT.0) WRITE(6,610) KVI,LWELL(1)
c** For "inner-well" levels, locate outer barrier
          DO  I= M,ITP3
              M2= I
              GA= V(I)- E
              IF(GA.GE.0.d0) GO TO 204
              ENDDO
          GO TO 218
          ENDIF 
      G3= V(M-2)- E
      G2= V(M-1)- E
      CALL LEVQAD(GA,G2,G3,H,RT,ANS1,ANS2)
      SM= SM- 0.5d0*DSQRT(G3)-DSQRT(G2) + ANS2/H
  218 EMSC= -SM/PI
      IF(INNER.GT.0) VMX= PMX
      VMAX= VMX/BFCT
c** Tunneling factors calculated here ** TUN0 is simple WKB result
c  as in Child's eqs.(57c) & (59).
c .....  EPSRJ= -2.* PI* EMSC 
      TUN0= 0.5d0*DEXP(2.d0*PI*EMSC)
c ... for permeability calculate Connor-Smith's Eq.(3.7) \omega=OMEGJC 
      OMEGJC= DSQRT(1.d0+ 2.d0*TUN0) - 1.d0
c ... alternate calculation to give better precision for small TUN0
      IF(TUN0.LT.1.d-5) OMEGJC= TUN0*(1.d0-0.5d0*TUN0*(1.d0-TUN0))
      OMEGJC= 4.d0*OMEGJC/(OMEGJC + 2.d0)
c** Quadrature for JWKB calculation of vibrational spacing in well HBW
      D1= E- V(IRM)
      D2= E- V(ITP1)
      D3= E- V(ITP1P)
      CALL LEVQAD(D1,D2,D3,H,RT,ANS1,ANS2)
      RT1= RT
      SM= ANS1/H
      IF(D3.LT.0.d0) GO TO 228
      SM= SM+ 0.5d0/DSQRT(D3)
      DO  I= ITP1P1,ITP2M2
          IMM= I
          EMV= E- V(I)
          IF(EMV.LT.0.d0) GO TO 222
          SM= SM+ 1.d0/DSQRT(EMV)
          ENDDO
      D3= E- V(ITP2M2)
      D2= E- V(ITP2M)
      D1= E- V(ITP2)
      GO TO 226
c** If encounter a double-minimum well, take care here.
  222 D1= EMV
      D2= E- V(IMM-1)
      D3= E- V(IMM-2)
      IF(IWR.NE.0) WRITE(6,605) KV,JROT,EO
  226 CALL LEVQAD(D1,D2,D3,H,RT,ANS1,ANS2)
      RT2=RT
      SM=SM-0.5d0/DSQRT(D3) + ANS1/H
c** Get HBW in same energy units (1/cm) associated with BFCT
  228 HBW=2.d0*PI/(BFCT*SM)
c** HBW fix up suggested by Child uses his eqs.(48)&(62) for HBW
c** Derivative of complex gamma function argument calculated as
c  per eq.(6.1.27) in Abramowitz and Stegun.
      NST= DABS(EMSC)*1.D2
      NST= MAX0(NST,4)
      ARG= -1.963510026021423d0
      DO  I= 0,NST
          NN= I
          XX= I + 0.5d0
          TI= 1.d0/(XX*((XX/EMSC)**2 + 1.d0))
          ARG= ARG+TI
          IF(DABS(TI).LT.1.D-10) GO TO 233
          ENDDO
c ... and use continuum approximation for tail of summation (???)
  233 COR= 0.5d0*(EMSC/(NN+1.d0))**2
      ARG= ARG+ COR- COR**2
c** Now use WKL's Weber fx. approx for (?) derivative of barrier integral ..
      DWEB= (EO-VMAX)*BFCT/(H2*EMSC)
      DFI= (DLOG(DABS(EMSC)) - ARG)*BFCT/(H2*DWEB)
      HBWB= 1.d0/(1.d0/HBW + DFI/(2.d0*PI))
c** Width from formula (4.5) of  Connor & Smith, Mol.Phys.43,397(1981)
c [neglect time delay integral past barrier in their Eq.(4.16)].
      IF(EMSC.GT.-25.D0) THEN
          GAMA= (HBWB/(2.d0*PI))* OMEGJC
          TAU= 0.D0
          IF(GAMA.GT.1.D-60) TAU= 5.308837457D-12/GAMA
c** GAM0 = TUN0*HBW/PI  is the simple WKB width GAMMA(0) discussed by
c  Le Roy & Liu in J.C.P.69,3622(1978).
          IF(IWR.GT.0) WRITE(6,601) TAU,GAMA,HBWB,VMAX
        ELSE
          GAMALG= DLOG10(HBWB/(2.d0*PI))+2.d0*PI*EMSC/2.302585093D0
          TAULG= DLOG10(5.308837457D-12)-GAMALG
          IF(IWR.GT.0) WRITE(6,611) TAULG,GAMALG,HBWB,VMAX
        ENDIF
  250 RETURN
  601 FORMAT('    Lifetime=',1PD10.3,'(s)   Width=',D10.3,'   dG/dv=',
     1 0PF7.2,'   V(max)=',F9.2)
  602 FORMAT(' *** WARNING ***  For   v =',I3,'   J =',I3,'   cannot cal
     1culate width since barrier maximum beyond range')
  603 FORMAT(' *** For  J=',I3,'  E=',F9.2,'  R(3-rd) beyond range so tu
     1nneling calculation uses'/8X,'pure centrifugal potential with  J(a
     2pp)=',F7.2,'  for  R > R(max)=',F7.2)
  605 FORMAT(' **** CAUTION *** Width estimate only qualitative, as have
     1 a double-minimum well for   E(v=',I3,', J=',I3,')=',F15.7/15X,
     2 'a more stable result may be obtained by searching for the quasib
     3ound levels using option: INNER > 0 .')
  609 FORMAT(' *** CAUTION - Permeability estimate not exact as have a d
     1ouble-humped barrier:  E(v=',I3,', J=',I3,') =',G15.8,I6)
  610 FORMAT(16X,'(NOTE: this has the node count of a   v=',I3,2X,A5,
     1 '-well level')
  611 FORMAT(12X,'Log10(lifetime/sec)=',F10.5,' ;   Log10(width/cm-1)=',
     1 F10.5,'   Spacing=',G12.5,'   V(max)=',G14.7,'(cm-1)')
      END
c**********************************************************************
      SUBROUTINE LEVQAD(Y1,Y2,Y3,H,RT,ANS1,ANS2)
c** Subroutine "LEVQAD" fits quadratic  Y = A + B*X + C*X**2  through
c  function values  Y1, Y2, Y3  at equally spaced points separated by
c  distance H, where  Y1 < 0  and (Y2,Y3 .ge.0), locates the function
c  zero (at RT, relative to  X1 < X2 = 0) between points X1 & X2, and
c  evaluates the integral from RT to R3 of   1/sqrt(Y)  , called
c  ANS1, and the integral (same range) of  sqrt(Y) , which is ANS2
c** Alternately, if Y1 & Y3 both  < 0  and only the middle point
c  Y2.ge.0 ,   fit the points to:  Y = A - B*(X-X0)**2 , locate the
c  turning points between which  Y(X) > 0  and evaluate these integrals
c  on this interval.  *************************************************
c----------------------------------------------------------------------
      REAL*8  A,ANS1,ANS2,B,C,CQ,H,HPI,R1,R2,RCQ,RR,RT,SL3,SLT,
     1        X0,Y1,Y2,Y3,ZT
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IF((Y1.GE.0).OR.(Y2.LT.0)) GO TO 99
      DATA HPI/1.570796326794896D0/
      IF(Y3.LT.0.d0) GO TO 50
c** Here treat case where both 'Y2' & 'Y3' are positive
      IF(DABS((Y2-Y1)/(Y3-Y2) -1.D0).LT.1.d-10) THEN
c ... special case of true (to 1/10^10) linearity ...
          RT= -H*Y2/(Y2-Y1)
          ANS1= 2.d0*(H-RT)/DSQRT(Y3)
          ANS2= ANS1*Y3/3.D0
          RETURN
          ENDIF
      C= (Y3-2.d0*Y2+Y1)/(2.d0*H*H)
      B= (Y3-Y2)/H-C*H
      A= Y2
      CQ= B**2- 4.d0*A*C
      RCQ= DSQRT(CQ)
      R1= (-B-RCQ)/(2.d0*C)
      R2= R1+ RCQ/C
      IF((R2.LE.0.d0).AND.(R2.GE.-H)) RT=R2
      IF((R1.LE.0.d0).AND.(R1.GE.-H)) RT=R1
      SL3= 2.d0*C*H+B
      SLT= 2.d0*C*RT+B
      IF(C.LT.0.d0) GO TO 10
      ANS1= DLOG((2.d0*DSQRT(C*Y3)+SL3)/SLT)/DSQRT(C)
      GO TO 20
   10 ANS1= -(DASIN(SL3/RCQ)- DSIGN(HPI,SLT))/DSQRT(-C)
   20 ANS2= (SL3*DSQRT(Y3)- CQ*ANS1/2.d0)/(4.d0*C)
      IF(RT.GE.H) WRITE(6,601) H,R1,R2
  601 FORMAT(' *** CAUTION *** in LEVQAD, turning point not between poin
     1ts 1 & 2.   H =',F9.6,'   R1 =',F9.6,'   R2 =',F9.6)
      RETURN
c** Here treat case when only 'Y2' is non-negative
   50 RR= (Y2-Y1)/(Y2-Y3)
      X0= H*(RR-1.d0)/((RR+1.d0)*2.d0)
      B= (Y2-Y1)/(H*(2.d0*X0+H))
      A= Y2+ B*X0**2
      ZT= DSQRT(A/B)
      RT= X0- ZT
      ANS1= 2.d0*HPI/DSQRT(B)
      ANS2= ANS1*A*0.5d0
      RETURN
   99 WRITE(6,602) Y1,Y2
  602 FORMAT(' *** ERROR in LEVQAD *** No turning point between 1-st two
     1 points as   Y1=',D10.3,'   Y2=',D10.3)
      ANS1= 0.d0
      ANS2= 0.d0
      RETURN
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE NLLSSRR(NDATA,NPTOT,NPMAX,IROUND,ROBUST,LPRINT,IFXP,
     1                           YO,YU,YD,PV,PU,PS,CM,TSTPS,TSTPU,DSE)
c**  Program for performing linear or non-linear least-squares fits and
c  (if desired) automatically using sequential rounding and refitting 
c  to minimize the numbers of parameter digits which must be quoted [see
c  R.J. Le Roy, J.Mol.Spectrosc. 191, 223-231 (1998)].         21/08/04
c
c  21/08/04 test version ... attempting to stablize non-linear fits by
c               scaling back parameter changes by factor of 4 or 16
c         [ & corrects  CM  for constrained-parameter cases!]
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c             COPYRIGHT 1998-2004  by  Robert J. Le Roy                +
c   Dept. of Chemistry, Univ. of Waterloo, Waterloo, Ontario, Canada   +
c    This software may not be sold or any other commercial use made    +
c      of it without the express written permission of the author.     +
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Program uses orthogonal decomposition of the "design" (partial 
c  derivative) matrix for the core locally linear (steepest descent) 
c  step, following a method introduced (to me) by Dr. Michael Dulick. 
c** If no parameters are free, simply return RMS(residuals) as
c  calculated from the input parameter values {PV(j)}.
c** A user MUST SUPPLY subroutine  DYIDPJ  to generate the predicted
c  value of each datum and the partial derivatives of each datum w.r.t.
c  each parameter (see below) from the current trial parameters.
c
c** On entry: 
c    NDATA  is the number of data to be fitted 
c    NPTOT  the total number of parameters in the model (.le.NPMAX).
c           If NPTOT.le.0 , assume  YD(i)=YO(i)  and calculate the (RMS 
c           dimensionless deviation)=DSE  from them & YU(i) 
c    NPMAX is the maximum number of model parameters allowed by current
c          external array sizes.  Should set internal NPINTMX = NPMAX 
c          (may be freely changed by the user).
c    IROUND .ne. 0  causes Sequential Rounding & Refitting to be 
c             performed, with each parameter being rounded at the 
c            |IROUND|'th sig. digit of its local incertainty.
c        > 0  rounding selects in turn remaining parameter with largest
c             relative uncertainy
c        < 0  round parameters sequentially from last to first
c        = 0  simply stops after full convergence (without rounding).
c    ROBUST > 0  causes fits to use Watson's ``robust'' weighting  
c        1/[u^2 +{(c-o)^2}/3].  ROBUST > 1 uses normal 1/u^2 on first
c        fit cycle and  'robust' on later cycles.
c    LPRINT  specifies the level of printing inside NLLSSRR
c          if: =  0, no print except for failed convergence.
c               < 0  only converged, unrounded parameters, PU & PS's
c              >= 1  print converged parameters, PU & PS's
c              >= 2  also print parameter change each rounding step
c              >= 3  also indicate nature of convergence
c              >= 4  also print convergence tests on each cycle
c              >= 5  also parameters changes & uncertainties, each cycle
c    IFXP(j)  specifies whether parameter  j  is to be held fixed
c           [IFXP > 0] or to be freely varied in the fit [IFXP= 0]
c    YO(i)  are the NDATA 'observed' data to be fitted  
c    YU(i)  are the uncertainties in these YO(i) values
c    PV(j)  are initial trial parameter values (for non-linear fits);  
c           should be set at zero for initially undefined parameters.
c
c** On Exit:   
c    YD(i)  is the array of differences  [Ycalc(i) - YO(i)]
c    PV(j)  are the final converged parameter values
c    PU(j)  are 95% confidence limit uncertainties in the PV(j)'s
c    PS(j)  are 'parameter sensitivities' for the PV(j)'s, defined such 
c           that the RMS displacement of predicted data  due to rounding
c           off parameter-j by PS(j) is .le. DSE/10*NPTOT
c    CM(j,k)  is the correlation matrix obtained by normalizing variance
c           /covariance matrix:  CM(j,k) = CM(j,k)/SQRT[CM(j,j)*CM(k,k)]
c    TSTPS = max{|delta[PV(j)]/PS(j)|}  is the parameter sensitivity 
c          convergence test:  delta[PV(j)] is last change in parameter-j
c    TSTPU = max{|delta[PV(j)]/PU(j)|}  is the parameter uncertainty 
c          convergence test:  delta[PV(j)] is last change in parameter-j
c    DSE    is the predicted (dimensionless) standard error of the fit
c
c  NOTE that the squared 95% confidence limit uncertainty in a property 
c  F({PV(j)}) defined in terms of the fitted parameters {PV(j)} (where
c  the L.H.S. involves  [row]*[matrix]*[column]  multiplication) is:
c  [D(F)]^2 = [PU(1)*dF/dPV(1), PU(2)*dF/dPV(2), ...]*[CM(j,k)]*
c                              [PU(2)*dF/dPV(1), PU(2)*dF/dPV(2), ...]
c
c** Externally dimension:  YO, YU and YD  .ge. NDATA 
c             PV, PU  and  PS  .ge.  NPTOT (say as NPMAX), 
c             CM   as a square matrix with column & row length  NPMAX
c***********************************************************************
      INTEGER NPINTMX
      PARAMETER (NPINTMX=2000)
      INTEGER I,J,K,L,IDF,ITER,NITER,IROUND,ISCAL,JROUND,LPRINT,NDATA,
     1 NPTOT,NPMAX,NPARM,NPFIT,JFIX,QUIT,ROBUST,
     2 IFXP(NPMAX),JFXP(NPINTMX)
      REAL*8  YO(NDATA), YU(NDATA), YD(NDATA), PV(NPTOT), PU(NPTOT), 
     1 PS(NPTOT),PSS(NPINTMX),PC(NPINTMX),PCS(NPINTMX),PX(NPINTMX),
     2 PY(NPINTMX),CM(NPMAX,NPMAX), F95(10),
     3 RMSR, RMSRB, DSE, TSTPS, TSTPSB, TSTPU, TFACT, S, UU, Zthrd
      DATA F95/12.7062D0,4.3027D0,3.1824D0,2.7764D0,2.5706D0,2.4469D0,
     1  2.3646D0,2.3060D0,2.2622D0,2.2281D0/
c
      IF((NPTOT.GT.NPMAX).OR.(NPTOT.GT.NPINTMX)
     1                    .OR.(NPTOT.GT.NDATA)) THEN
c** If array dimensioning inadequate, print warning & then STOP
          WRITE(6,602) NPTOT,NPINTMX,NPMAX,NDATA
          STOP
          ENDIF
      Zthrd= 0.d0
      IF(ROBUST.GE.2) Zthrd= 1.d0/3.d0
      TSTPS= 0.d0
      RMSR= 0.d0
      NITER= 0
      QUIT= 0
      NPARM= NPTOT
      DO J= 1, NPTOT
          PS(J)= 0.d0
          JFXP(J)= IFXP(J)
          IF(IFXP(J).GT.0) NPARM= NPARM- 1
          ENDDO
      NPFIT= NPARM
      JROUND= IABS(IROUND)
c=======================================================================
c** Beginning of loop to perform rounding (if desired).  NOTE that in 
c  sequential rounding, NPARM is the current (iteratively shrinking) 
c  number of free parameters. 
    6 IF(NPARM.GT.0) TSTPS= 9.d99
c** TFACT  is 95% student t-value for (NDATA-NPARM) degrees of freedom.
c [Approximate expression for (NDATA-NPARM).GT.10 accurate to ca. 0.002]
      TFACT= 0.D0
      IF(NDATA.GT.NPARM) THEN
          IDF= NDATA-NPARM
          IF(IDF.GT.10) TFACT= 1.960D0*DEXP(1.265D0/DFLOAT(IDF))
          IF(IDF.LE.10) TFACT= F95(IDF) 
        ELSE
          TFACT= 0.D0
        ENDIF
c======================================================================
c** Begin iterative convergence loop:  try for up to 30 cycles
      DO 50 ITER= 1, 30
          ISCAL= 0
          NITER= NITER+ 1
          DSE= 0.d0 
          TSTPSB= TSTPS
          RMSRB= RMSR
c** Zero out various arrays
          IF(NPARM.GT.0) THEN
              DO  I = 1,NPARM
c** PSS is the array of Saved Parameter Sensitivities from previous 
c   iteration to be carried into dyidpj subroutine - used in predicting
c   increment for derivatives by differences.
                  PSS(I)= PS(I)
c** PCS is the saved array of parameter changes from previous iteration
c   to be used (if necessary) to attempt to stablize fit
                  PCS(I)= PC(I)
                  PS(I) = 0.D0
                  PU(I) = 0.D0
                  PX(I) = 0.D0
                  PY(I) = 0.D0
                  DO  J = 1,NPARM
                      CM(I,J) = 0.D0
                      ENDDO
                  ENDDO
              ENDIF
c
c========Beginning of core linear least-squares step====================
c
c** Begin by forming the Jacobian Matrix from partial derivative matrix
          DO  I = 1,NDATA
c** User-supplied subroutine DYIDPJ uses current (trial) parameter 
c  values {PV} to generate predicted datum # I [y(calc;I)=UU] and its
c  partial derivatives w.r.t. each of the parameters, returning the 
c  latter in 1-D array PC.  See dummy sample version at end of listing.
c* NOTE 1: if more convenient, DYIDPJ could prepare the y(calc) values 
c     and derivatives for all data at the same time (when I=1), but only
c     returned the values here one datum at a time (for I > 1).]
c* NOTE 2: the partial derivative array PC returned by DYIDPJ must have
c     an entry for every parameter in the model, though for parameters 
c     which are held fixed [JFXP(j)=1], those PC(j) values are ignored.
              CALL DYIDPJ(I,NDATA,NPTOT,JFXP,YO(I),UU,PV,PC,PSS,RMSR)
              IF(NPARM.LT.NPTOT) THEN
c** For constrained parameter or sequential rounding, collapse partial 
c   derivative array here
                  DO  J= NPTOT,1,-1
                      IF(JFXP(J).GT.0) THEN
                          IF(J.LT.NPTOT) THEN
                              DO  K= J,NPTOT-1
                                  PC(K)= PC(K+1)
                                  ENDDO
                              ENDIF
                          PC(NPTOT)= 0.d0
                          ENDIF
                      ENDDO
                  ENDIF
              YD(I)= UU - YO(I)
              S = 1.D0 / YU(I)
              IF(Zthrd.GT.0.d0) S= 1.d0/DSQRT(YU(I)**2 + Zthrd*YD(I)**2)
              UU = - YD(I) * S
              DSE= DSE+ UU*UU
              IF(NPARM.GT.0) THEN
                  DO  J = 1,NPARM
                      PC(J) = PC(J)*S
                      PS(J) = PS(J)+ PC(J)**2
                      ENDDO
                  CALL QROD(NPARM,NPMAX,NPMAX,CM,PC,PU,UU,PX,PY)
                  ENDIF
              ENDDO
          RMSR= DSQRT(DSE/NDATA)
          IF(NPARM.LE.0) GO TO 60

c         IF((ITER.GT.1).AND.(RMSR.GT.RMSRB).AND.(ISCAL.LE.1)) THEN
c** LeRoy's Marquardt-like attempt to damp changes if RMSR increases ...
c             ISCAL= ISCAL+ 1
c             IF(LPRINT.GE.5) THEN
c                 WRITE(6,620) ITER,RMSR/RMSRB,ISCAL
c 620 FORMAT(' At Iteration',i3,'  with  RMSD/RMSDB=',1PD8.1,
c    1 "  Scale Param Changes by  (1/4)**",i1)
c                 WRITE(6,612) (J,PV(J),PU(J),PS(J),PC(J),J=1,NPTOT)
c                 ENDIF
c             DO  J= 1,NPTOT
c                 PC(J)= 0.25d0*PCS(J)
c                 PV(J)= PV(J)- 3.d0*PC(J)
c                 ENDDO
c             GOTO 10
c             ENDIF

c
c** Compute the inverse of  CM 
          CM(1,1) = 1.D0 / CM(1,1)
          DO  I = 2,NPARM
              L = I - 1
              DO  J = 1,L
                  S = 0.D0
                  DO  K = J,L
                      S = S + CM(K,I) * CM(J,K)
                      ENDDO
                  CM(J,I) = -S / CM(I,I)
                  ENDDO
              CM(I,I) = 1.D0 / CM(I,I)
              ENDDO
c
c** Solve for parameter changes  PC(j)
          DO  I = 1,NPARM
              J = NPARM - I + 1
              PC(J) = 0.D0
              DO  K = J,NPARM
                  PC(J) = PC(J) + CM(J,K) * PU(K)
                  ENDDO
              ENDDO
c
c** Get (upper triangular) "dispersion Matrix" [variance-covarience 
c  matrix  without the sigma^2 factor].
          DO  I = 1,NPARM
              DO  J = I,NPARM
                  UU = 0.D0
                  DO  K = J,NPARM
                      UU = UU + CM(I,K) * CM(J,K)
                      ENDDO
                  CM(I,J) = UU
                  ENDDO
              ENDDO
c** Generate core of Parameter Uncertainties  PU(j) and (symmetric)
c   correlation matrix  CM
          DO  J = 1,NPARM
              PU(J) = DSQRT(CM(J,J))
              DO  K= J,NPARM
                  CM(J,K)= CM(J,K)/PU(J)
                  ENDDO
              DO  K= 1,J
                  CM(K,J)= CM(K,J)/PU(J)
                  CM(J,K)= CM(K,J)
                  ENDDO
              ENDDO
c
c** Generate standard error  DSE = sigma^2,  and prepare to calculate 
c  Parameter Sensitivities PS
          IF(NDATA.GT.NPARM) THEN
              DSE= DSQRT(DSE/(NDATA-NPARM))
            ELSE
              DSE= 0.d0
            ENDIF
c** Use DSE to get final (95% confid. limit) parameter uncertainties PU
c** Calculate 'parameter sensitivities', changes in PV(j) which would 
c  change predictions of input data by an RMS average of  DSE*0.1/NPARM
          UU= DSE*0.1d0/DFLOAT(NPARM)
          S= DSE*TFACT
          DO  J = 1,NPARM
              PU(J)= S* PU(J)
              PS(J)= UU*DSQRT(NDATA/PS(J))
              ENDDO
c========End of core linear least-squares step==========================
c ... early exit if Rounding cycle finished ... 
          IF(QUIT.GT.0) GO TO 54
c
c** Next test for convergence 
          TSTPS= 0.D0
          TSTPU= 0.D0
          DO  J= 1, NPARM
              TSTPS= MAX(TSTPS,DABS(PC(J)/PS(J)))
              TSTPU= MAX(TSTPU,DABS(PC(J)/PU(J)))
              ENDDO
          IF(LPRINT.GE.4) WRITE(6,604) ITER,RMSR,TSTPS,TSTPU
c** Now ... update parameters (careful about rounding)
          DO  J= 1,NPTOT
              IF(JFXP(J).GT.0) THEN
c** If parameter held fixed (by input or rounding process), shift values
c   of change, sensitivity & uncertainty to correct label.
                  IF(J.LT.NPTOT) THEN
                      DO  I= NPTOT,J+1,-1
                          PC(I)= PC(I-1)
                          PS(I)= PS(I-1)
                          PU(I)= PU(I-1)
                          ENDDO
                      ENDIF
                  PC(J)= 0.d0
                  PS(J)= 0.d0
                  PU(J)= 0.d0
                ELSE
                  PV(J)= PV(J)+ PC(J)
                ENDIF
              ENDDO
          IF(LPRINT.GE.5) WRITE(6,612) (J,PV(J),PU(J),PS(J),PC(J),
     1                                                      J=1,NPTOT)
          IF(NITER.GT.1) THEN
c** Test for convergence:  for every parameter desire:
c  |parameter change| < |parameter sensitivity|,  but after iteration #5
c  STOP iterating if  Max{|change/sens.|} increases AND 
c  Max{|change/unc.|} < 0.01
              IF(TSTPS.GT.1.d0) THEN
                  IF((RMSR.GT.RMSRB).AND.(ITER.GT.5)) THEN
                      IF((TSTPU.LT.1.d-2).OR.((TSTPU.LT.0.5d0).AND.
     1                                             (ITER.GT.10))) THEN
                          IF(LPRINT.GE.3) WRITE(6,606) ITER,TSTPU,RMSR
                          GO TO 54
                          ENDIF
                      ENDIF
                ELSE
                  IF(LPRINT.GE.3) WRITE(6,608) ITER,TSTPS,RMSR
                  GO TO 54
                ENDIF
              ENDIF
cc        CALL FLUSH(6)
          IF(ROBUST.GT.0) Zthrd= 1.d0/3.d0
   50     CONTINUE
      WRITE(6,610) NPARM,NDATA,ITER,RMSR,TSTPS,TSTPU
c** End of iterative convergence loop for (in general) non-linear case.
c======================================================================
c
   54 IF(NPARM.LT.NPTOT) THEN
c** If necessary, redistribute correlation matrix elements to full 
c  NPTOT-element correlation matrix
          DO  J= 1,NPTOT
              IF(JFXP(J).GT.0) THEN  
c* If parameter J was held fixed
                  IF(J.LT.NPTOT) THEN
c ... then move every lower CM element down one row:
                      DO  I= NPTOT,J+1,-1
c ... For  K < J, just shift down or over to the right
                          IF(J.GT.1) THEN
                              DO  K= 1,J-1
                                  CM(I,K)= CM(I-1,K) 
                                  CM(K,I)= CM(I,K)
                                  ENDDO
                              ENDIF
c ... while for  K > J  also shift elements one column to the right
                          DO  K= NPTOT,J+1,-1
                              CM(I,K)= CM(I-1,K-1)
                              ENDDO
                          ENDDO
                      ENDIF
c ... and finally, insert appropriate row/column of zeros ....
                  DO  I= 1,NPTOT
                      CM(I,J)= 0.d0
                      CM(J,I)= 0.d0
                      ENDDO 
                  CM(J,J)= 1.d0
                  ENDIF
              ENDDO
          ENDIF
      IF(QUIT.GT.0) GOTO 60
      IF(NPARM.EQ.NPFIT) THEN
c** If desired, print unrounded parameters and fit properties
          IF(LPRINT.NE.0) THEN
              WRITE(6,616) NDATA,NPARM,RMSR,TSTPS
              WRITE(6,612) (J,PV(J),PU(J),PS(J),PC(J),J=1,NPTOT)
              ENDIF
          ENDIF
      IF(IROUND.EQ.0) RETURN
c** Automated 'Sequential Rounding and Refitting' section:  round 
c  selected parameter, fix it, and return (above) to repeat fit.
      IF(IROUND.LT.0) THEN
c ... if IROUND < 0, sequentially round off 'last' remaining parameter
          DO  J= 1, NPTOT
              IF(JFXP(J).LE.0) THEN
                  JFIX= J
                  ENDIF
              ENDDO
        ELSE
c ... if IROUND > 0, sequentially round off remaining parameter with
c                    largest relative uncertainty.
c ... First, select parameter JFIX with the largest relative uncertainty
          K= 0
          TSTPS= 0.d0
          DO  J= 1,NPTOT
              IF(JFXP(J).LE.0) THEN
                  K= K+1
                  TSTPSB= DABS(PU(J)/PV(J))
                  IF(TSTPSB.GT.TSTPS) THEN
                      JFIX= J
                      TSTPS= TSTPSB
                      ENDIF
                  ENDIF
              ENDDO 
        ENDIF
      UU= PV(JFIX)
      CALL ROUND(JROUND,NPMAX,NPTOT,NPTOT,JFIX,PV,PU,PS,CM)
      JFXP(JFIX)= 1
      IF(LPRINT.GE.2)
     1            WRITE(6,614) JFIX,UU,PU(JFIX),PS(JFIX),PV(JFIX),RMSR
      NPARM= NPARM-1
      IF(NPARM.EQ.0) THEN
c** After rounding complete, make one more pass with all non-fixed 
c  parameters set free to get full correct final correlation matrix, 
c  uncertainties & sensitivities
          NPARM= NPFIT
          QUIT= 1
          DO  J= 1,NPTOT
              JFXP(J)= IFXP(J)
              ENDDO
c ... reinitialize for derivative-by-differences calculation
          RMSR= 0.d0
          ENDIF
      GO TO 6
c
c** If no parameters varied or sequential rounding completed - simply 
c   calculate DSE from RMS residuals and return.
   60 DSE= 0.d0
      IF(NDATA.GT.NPFIT) THEN
          DSE= RMSR*DSQRT(DFLOAT(NDATA)/DFLOAT(NDATA-NPFIT))
        ELSE
          DSE= 0.d0
        ENDIF
      IF(NPFIT.GT.0) THEN
          IF(LPRINT.GT.0) THEN
c** Print final rounded parameters with original Uncert. & Sensitivities
              IF(QUIT.LT.1) WRITE(6,616) NDATA, NPFIT, RMSR, TSTPS
              IF(QUIT.EQ.1) WRITE(6,616) NDATA, NPFIT, RMSR
              DO  J= 1, NPTOT
                  IF(JFXP(J).GT.0) THEN
c** If parameter held fixed (by rounding process), shift values of
c   change, sensitivity & uncertainty to correct absolute number label.
                      DO  I= NPTOT,J+1,-1
                          PC(I)= PC(I-1)
                          PS(I)= PS(I-1)
                          PU(I)= PU(I-1)
                          ENDDO
                      PC(J)= 0.d0
                      PS(J)= 0.d0
                      PU(J)= 0.d0
                      ENDIF
                  ENDDO
              WRITE(6,612) (J,PV(J),PU(J),PS(J),PC(J),J=1,NPTOT)
              ENDIF
          ENDIF
      RETURN
c
  602 FORMAT(/' *** NLLSSRR problem:  [NPTOT=',i4,'] > min{NPINTMX=',
     1  i4,' NPMAX=',i4,', NDATA=',i6,'}')
  604 FORMAT(' After Cycle #',i2,':   DRMSD=',1PD10.3,'    test(PS)=',
     1  1PD8.1,'    test(PU)=',D8.1)
  606 FORMAT(' Effective',i3,'-cycle Cgce:  MAX{|change/unc.|}=',1PD8.1,
     1  ' < 0.01   DRMSD=',D10.3)
  608 FORMAT(' Full',i3,'-cycle convergence:  Max{|change/sens.|}=',
     1  1PD8.1,' < 1   DRMSD=',D10.2)
  610 FORMAT(' !! CAUTION !! fit of',i4,' parameters to',I6,' data not c
     1onverged after',i3,' Cycles'/5x,'DRMS(deviations)=',1PD10.3,
     2 '    test(PS) =',D9.2,'    test(PU) =',D9.2/1x,31('**'))
  612 FORMAT((4x,'PV(',i4,') =',1PD22.14,' (+/-',D8.1,')    PS=',d8.1,
     1  '   PC=',d8.1))
  614 FORMAT(' =',39('==')/' Round Off  PV(',i4,')=',1PD21.13,' (+/-',
     1 D9.2,')    PS=',d9.2/10x,'fix it as ',D21.13,'  & refit:  DRMS(de
     2viations)=',D10.3)
  616 FORMAT(/i6,' data fit to',i5,' param. yields   DRMS(devn)=',G11.4:
     1 '  tst(PS)=',1Pd8.1)
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
	SUBROUTINE QROD(N,NR,NC,A,R,F,B,GC,GS)
C** Performs ORTHOGONAL DECOMPOSITION OF THE LINEAR LEAST-SQUARES    
C            EQUATION J * X = F TO A * X = B(TRANSPOSE) * F WHERE   
C            J IS THE JACOBIAN IN WHICH THE FIRST N ROWS AND COLUMNS
C            ARE TRANSFORMED TO THE UPPER TRIANGULAR MATRIX A      
C            (J = B * A), X IS THE INDEPENDENT VARIABLE VECTOR, AND
C            F IS THE DEPENDENT VARIABLE VECTOR. THE TRANSFORMATION
C            IS APPLIED TO ONE ROW OF THE JACOBIAN MATRIX AT A TIME.
C  PARAMETERS :                                                   
C      N   -  (INTEGER) DIMENSION OF A TO BE TRANSFORMED.        
C      NR  -  (INTEGER) ROW DIMENSION OF A DECLARED IN CALLING PROGRAM.
C      NC  -  (INTEGER) Column DIMENSION OF F DECLARED IN CALLING PROGRAM.
C      A   -  (REAL*8 ARRAY OF DIMENSIONS .GE. N*N) UPPER TRIANGULAR
C             TRANSFORMATION MATRIX.                               
C      R   -  (REAL*8 LINEAR ARRAY OF DIMENSION .GE. N) ROW OF    
C             JACOBIAN TO BE ADDED.                             
C      F   -  (REAL*8 LINEAR ARRAY .GE. TO THE ROW DIMENSION OF THE
C             JACOBIAN) TRANSFORMED DEPENDENT VARIABLE MATRIX.    
C      B   -  (REAL*8) VALUE OF F THAT CORRESPONDS TO THE ADDED  
C             JACOBIAN ROW.                                     
C     GC   -  (REAL*8 LINEAR ARRAY .GE. N) GIVENS COSINE TRANSFORMATIONS.
C     GS   -  (REAL*8 LINEAR ARRAY .GE. N) GIVENS SINE TRANSFORMATIONS. 
C--------------------------------------------------------------------
C  AUTHOR : MICHAEL DULICK, Department of Chemistry,
C           UNIVERSITY OF WATERLOO, WATERLOO, ONTARIO N2L 3G1
C--------------------------------------------------------------------
      INTEGER  I,J,K,N,NC,NR
      REAL*8 A(NR,NC), R(N), F(NR), GC(N), GS(N), B, Z(2)
      DO 10 I = 1,N
          Z(1) = R(I)
          J = I - 1
          DO  K = 1,J
              Z(2) = GC(K) * A(K,I) + GS(K) * Z(1)
              Z(1) = GC(K) * Z(1) - GS(K) * A(K,I)
              A(K,I) = Z(2)
              ENDDO
          GC(I) = 1.D0
          GS(I) = 0.D0
          IF(Z(1) .EQ. 0.D0) GOTO 10
          IF(DABS(A(I,I)) .LT. DABS(Z(1))) THEN
              Z(2) = A(I,I) / Z(1)
              GS(I) = 1.D0 / DSQRT(1.D0 + Z(2) * Z(2))
              GC(I) = Z(2) * GS(I)
            ELSE
              Z(2) = Z(1) / A(I,I)
              GC(I) = 1.D0 / DSQRT(1.D0 + Z(2) * Z(2))
              GS(I) = Z(2) * GC(I)
            ENDIF
          A(I,I) = GC(I) * A(I,I) + GS(I) * Z(1)
          Z(2) = GC(I) * F(I) + GS(I) * B
          B = GC(I) * B - GS(I) * F(I)
          F(I) = Z(2)
   10     CONTINUE
      RETURN
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE ROUND(IROUND,NPMAX,NPARM,NPTOT,IPAR,PV,PU,PS,CM)
c** Subroutine to round off parameter # IPAR with value PV(IPAR) at the
c  |IROUND|'th significant digit of:  [its uncertainty  PU(IPAR)] . 
c** On return, the rounded value replaced the initial value  PV(IPAR).
c** Then ... use the correlation matrix CM and the uncertainties PU(I)
c  in the other (NPTOT-1) [or (NPARM-1) free] parameters to calculate 
c  the optimum compensating changes PV(I) in their values.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c                COPYRIGHT 1998  by  Robert J. Le Roy                  +
c   Dept. of Chemistry, Univ. of Waterloo, Waterloo, Ontario, Canada   +
c    This software may not be sold or any other commercial use made    +
c      of it without the express written permission of the author.     +
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      INTEGER    IROUND,NPMAX,NPARM,NPTOT,IPAR,I,IRND,KRND
      REAL*8  PU(NPMAX),PS(NPMAX),PV(NPMAX),CM(NPMAX,NPMAX),CNST,
     1        CRND,XRND,FCT,Z0
      DATA Z0/0.d0/
      CNST= PV(IPAR)
      XRND= DLOG10(PU(IPAR))
c** If appropriate, base last rounding step on sensitivity (not uncert.)
      IF((NPARM.EQ.1).AND.(PS(IPAR).LT.PU(IPAR))) XRND= DLOG10(PS(IPAR))
c** First ... fiddle with log's to perform the rounding
      IRND= INT(XRND)
      IF(XRND.GT.0) IRND=IRND+1
      IRND= IRND- IROUND
      FCT= 10.D0**IRND
      CRND= PV(IPAR)/FCT
      XRND= Z0
c ... if rounding goes past REAL*8 precision, retain unrounded constant
      IF(DABS(CRND).GE.1.D+16) THEN
          WRITE(6,601) IROUND,IPAR
           RETURN
           ENDIF
      IF(DABS(CRND).GE.1.D+8) THEN
c ... to avoid problems from overflow of I*4 integers ...
          KRND= NINT(CRND/1.D+8)
          XRND= KRND*1.D+8
          CRND= CRND-XRND
          XRND= XRND*FCT
          END IF
      IRND= NINT(CRND)
      CNST= IRND*FCT+ XRND
c????????????????
c** Zero parameters more aggressively ... if unc. > 2* value
        if(dabs(PU(IPAR)/PV(IPAR)).GT.2.d0) then
            cnst= 0.d0
            endif
c????????????????
c** Now ... combine rounding change in parameter # IPAR, together with
c  correlation matrix CM and parameter uncertainties PU to predict
c  changes in other parameters to optimally compensate for rounding off
c  of parameter-IPAR.  Method pointed out by Mary Thompson (Dept. of
c  Statistics, UW),
      IF(IPAR.GT.1) THEN
          XRND= (CNST-PV(IPAR))/PU(IPAR)
          DO  I= 1,NPTOT
              IF(I.NE.IPAR) THEN
                  PV(I)= PV(I)+ CM(IPAR,I)*PU(I)*XRND
                  ENDIF
              ENDDO
          ENDIF
      PV(IPAR)= CNST
      RETURN
  601 FORMAT(' =',39('==')/' Caution:',i3,'-digit rounding of parameter-
     1',i2,' would exceed (assumed) REAL*8'/' ********   precision overf
     2low at 1.D+16, so keep unrounded constant')
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
c     SUBROUTINE DYIDPJ(I,NDATA,NPTOT,IFXP,UU,PV,PD,PS,RMSR)
c** Illustrative dummy version of DYIDPJ for the case of a fit to a
c  power series of order (NPTOT-1) in X(i). ***  For datum number-i, 
c  calculate and return  PD(j)=[partial derivatives of datum-i] w.r.t. 
c  each of the free polynomial coefficients varied in the fit 
c  (for j=1 to NPTOT).  **  Elements of the integer array IFXP indicate
c  whether parameter j is being held fixed [IFXP(j) > 0] or varied in
c  the fit [IFXP(j).le.0].  If the former, the partial derivative 
c  for parameter j should be  PD(j)= 0.0. 
c* NOTE that  NDATA, PS and RMSR are useful for cases in which
c  derivatives-by-differences are generated (as for BCONT).
c=====================================================================
c** Use COMMON block(s) to bring in values of the independent variable 
c  [here XX(i)] and any other parameters or variables needeed to
c  calculate YC and the partial derivatives. 
c=====================================================================
c     INTEGER  I,J,NDATA,NPTOT,MXDATA,IFXP(NPTOT)
c     PARAMETER  (MXDATA= 501)
c     REAL*8  RMSR,YC,PV(NPTOT),PD(NPTOT),PS(NPTOT),POWER,XX(MXDATA)
c     COMMON /DATABLK/XX
c=====================================================================
c** NOTE BENE(!!) for non-linear fits, need to be sure that the
c  calculations of YC and PD(j) are based on the current UPDATED PV(j)
c  values.  If other (than PV) parameter labels are used internally
c  in the calculations, UPDATE them whenever (say)  I = 1 .
c=====================================================================
c     POWER= 1.D0
c     YC= PV(1)
c     PD(1)= POWER
c     DO 10 J= 2,NPTOT
c         POWER= POWER*XX(I)
c         YC= YC+ PV(J)*POWER
c         PD(J)= POWER
c  10     CONTINUE
c     RETURN
c     END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE FUNUNC(ISTATE,WRITFILE,OSEL,PU,CM)
c***********************************************************************
c** This subroutine will calculate the uncertainties in the radial
c  strength functions for V(r), phi(r), UA(r), UB(r), qA(r), qB(r) and 
c  wRAD(r) and print them out in `Tecplot' format.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** On entry:
c ... OSEL     print values at every OSEL'th mesh point
c ... ISTATE   electronic state counter
c ... PU(n)    uncertainties in the TOTPOTPAR parameters of the model
c ... CM(n,n)  (symmetric) correlation matrix from the fit
c=======================================================================
      INCLUDE 'arrsizes.h'
      INCLUDE 'BLKPOT.h'
      INCLUDE 'BLKDVDP.h'
      INCLUDE 'BLKBOB.h'
      INCLUDE 'BLKBOBRF.h'
      INCLUDE 'BLKCOUNT.h'
c-----------------------------------------------------------------------
      INTEGER I,J,JJ,ISTATE,LAM2,NPHII, OSEL(NSTATEMX)
      REAL*8 FU,FLAM,RDVAL,RDVAL2,RDVALLD, PU(NPARMX),PT(NPARMX),
     1  CM(NPARMX,NPARMX) 
      CHARACTER*20 WRITFILE
c
      REAL*8 Vsr(NPNTMX,NSTATEMX),Bsr(NPNTMX,NSTATEMX)
      INTEGER nPointSR(NSTATEMX)
      COMMON /VsrBLK/Vsr,Bsr,nPointSR
c
      REAL*8 Plr(NPNTMX,NSTATEMX),Blr(NPNTMX,NSTATEMX)
      INTEGER nPointLR(NSTATEMX)
      COMMON /PlrBLK/Plr,Blr,nPointLR
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      WRITE(10,900) 'V(r) ', ISTATE, 'V(r) ', WRITFILE
      DO I= 1,nPointSR(ISTATE),OSEL(ISTATE)
          RDVAL = RH(ISTATE)*DBLE(I-1) 
          WRITE(10,909) RDVAL,Vsr(I,ISTATE)
          END DO   
      WRITE(11,900) 'B(r) ', ISTATE, 'B(r) ', WRITFILE
      DO I= 1,nPointSR(ISTATE),OSEL(ISTATE)
          RDVAL = RH(ISTATE)*DBLE(I-1)
          WRITE(11,909) RDVAL,Bsr(I,ISTATE)
          END DO
      IF(NUA(ISTATE).GE.0) WRITE(12,900) 'UA(r)',ISTATE,'UA(r)',
     1                                                WRITFILE
      IF(NUB(ISTATE).GE.0) WRITE(13,900) 'UB(r)',ISTATE,'UB(r)',
     1                                                WRITFILE
      IF(NTA(ISTATE).GE.0) WRITE(14,900) 'qA(r)',ISTATE,'qA(r)',
     1                                                WRITFILE
      IF(NTB(ISTATE).GE.0) WRITE(15,900) 'qB(r)',ISTATE,'qB(r)',
     1                                                WRITFILE
      IF(NwCFT(ISTATE).GE.0) THEN
          WRITE(16,902) 'lambda(r)',ISTATE,'lambda(r)',WRITFILE
          LAM2= 2
          IF(IOMEG(ISTATE).GT.0) LAM2= 4*IOMEG(ISTATE)
          ENDIF
      NPHII= POTPARF(ISTATE) - MAX(NSphi(ISTATE),NLphi(ISTATE))
      DO  I= 1,NDATPT(ISTATE),OSEL(ISTATE)
          DO J= 1,TOTPOTPAR
              PT(J) = 0.0d0
              END DO
          RDVAL = RMIN(ISTATE) + RH(ISTATE) * DBLE(I-1)
          RDVAL2= RDVAL*RDVAL
          RDVALLD= RDVAL**LAM2
c ... first ... the potential function itself ...
          DO  J= POTPARI(ISTATE),POTPARF(ISTATE)
              PT(J)= PU(J)*DVtot(J,I)
              ENDDO
          CALL MMCALC(POTPARI(ISTATE),POTPARF(ISTATE),PT,CM,FU)
          WRITE(10,910) RDVAL,VPOT(I,ISTATE),FU
c ... then the exponent function \phi(i)
          IF(PSEL(ISTATE).NE.4) THEN
              DO  J= NPHII,POTPARF(ISTATE)
                  JJ= J-NPHII
                  PT(J)= PU(J)*DBDB(JJ,I,ISTATE)
                  ENDDO
              CALL MMCALC(NPHII,POTPARF(ISTATE),PT,CM,FU)
              WRITE(11,910) RDVAL,PHIFX(I,ISTATE),FU
              ENDIF
c ... adiabatic BOB correction function for atom-A
          IF(NUA(ISTATE).GE.0) THEN
              DO  J= UAPARI(ISTATE),UAPARF(ISTATE)
                  PT(J)= PU(J)*DVtot(J,I)
                  ENDDO
              CALL MMCALC(UAPARI(ISTATE),UAPARF(ISTATE),PT,CM,FU)
              WRITE(12,910) RDVAL,UAR(I,ISTATE),FU
              ENDIF
c ... adiabatic BOB correction function for atom-B
          IF(NUB(ISTATE).GE.0) THEN
              DO  J= UBPARI(ISTATE),UBPARF(ISTATE)
                  PT(J)= PU(J)*DVtot(J,I)
                  ENDDO
              CALL MMCALC(UBPARI(ISTATE),UBPARF(ISTATE),PT,CM,FU)
              WRITE(13,910) RDVAL,UBR(I,ISTATE),FU
              ENDIF
c ... centrifugal BOB correction function for atom-A
          IF(NTA(ISTATE).GE.0) THEN
              DO  J= TAPARI(ISTATE),TAPARF(ISTATE)
                  PT(J)= PU(J)*DVtot(J,I)*RDVAL2
                  ENDDO
              CALL MMCALC(TAPARI(ISTATE),TAPARF(ISTATE),PT,CM,FU)
              WRITE(14,910) RDVAL,TAR(I,ISTATE),FU
              ENDIF
c ... centrifugal BOB correction function for atom-B
          IF(NTB(ISTATE).GE.0) THEN
              DO  J= TBPARI(ISTATE),TBPARF(ISTATE)
                  PT(J)= PU(J)*DVtot(J,I)*RDVAL2
                  ENDDO
              CALL MMCALC(TBPARI(ISTATE),TBPARF(ISTATE),PT,CM,FU)
              WRITE(15,910) RDVAL,TBR(I,ISTATE),FU
              ENDIF
c ... Lambda/doublet-sigma doubling correction radial function 
          IF(NwCFT(ISTATE).GE.0) THEN
              DO  J= LDPARI(ISTATE),LDPARF(ISTATE)
                  PT(J)= PU(J)*DVtot(J,I)*RDVALLD
                  ENDDO
              CALL MMCALC(LDPARI(ISTATE),LDPARF(ISTATE),PT,CM,FU)
              FLAM= wRAD(I,ISTATE)*RDVALLD
              WRITE(16,910) RDVAL,FLAM,FU
              ENDIF
          END DO
      DO I= 10,nPointLR(ISTATE),10
          RDVAL= RMAX(ISTATE) + RH(ISTATE)*DBLE(I*OSEL(ISTATE))
          WRITE(10,909) RDVAL,Plr(I,ISTATE)
          WRITE(11,909) RDVAL,Blr(I,ISTATE)
          END DO
      RETURN
c-----------------------------------------------------------------------
  900 FORMAT(/'variables = "r", "',A5,'", "Uncertainty"'/
     1'zone T = "State',I2,1x,A5,2x,A20,'"',/)
  902 FORMAT(/'variables = "r", "',A9,'", "Uncertainty"'/
     1'zone T = "State',I2,2x,A9,2x,A20,'"',/)
  909 FORMAT(F12.6,1x,1PD22.13,5X,'0.0')
  910 FORMAT(F12.6,1x,1PD22.13,D11.3)
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE MMCALC(JI,JF,PT,CM,FU)
c***********************************************************************
c** For elements JI to JF of column vector PT(j) and rows/columns 
c  JI to JF of the square matrix CM, evaluate  FU**2 = PT^t CM PT
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Externally dimension  PT(NPARMX) and CM(NPARMX,NPARMX)
c***********************************************************************
      INCLUDE 'arrsizes.h'
      INTEGER I,J,JI,JF
      REAL*8  PT(NPARMX), CM(NPARMX,NPARMX), TVAL, FU
c=======================================================================
c** Initialize the variables.
      FU = 0.0d0
      DO  J= JI, JF
          TVAL= 0.d0
          DO  I= JI,JF
c** Multiply the vector PT with the correlation matrix (CM).
              TVAL= TVAL + PT(I)*CM(I,J)
              ENDDO
c** Now to multiply the vector (PT*CM) with the vector PT.
          FU= FU + TVAL * PT(J)
          ENDDO
c** Complete calculation of the uncertainty of the function.
      FU= DSQRT(FU)
      RETURN
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

C***********************************************************************
      SUBROUTINE GPROUND(IROUND,NPTOT,NPMAX,NPAR1,NPAR2,LPRINT,IFXP,
     1                                                          PV,PU)
c** Subroutine to round off parameters PV(i), i= NPAR1 to NPAR2, at the
c  |IROUND|'th significant digit of the smallest of their uncertainties
c  min{U(i)}.  This procedure does NOT attempt to correct the remaining
c  parameters to compensate for these changes (as ROUND does), so this
c  procedure is not appropriate for nonlinear parameters.
c** On return, the rounded values replaces the initial values of  PV(i).
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c                COPYRIGHT 2000-2004  by  Robert J. Le Roy             +
c   Dept. of Chemistry, Univ. of Waterloo, Waterloo, Ontario, Canada   +
c    This software may not be sold or any other commercial use made    +
c      of it without the express written permission of the author.     +
c                      Version of 27 January 2004                      +
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      INTEGER  IROUND,NPMAX,NPTOT,NPAR1,NPAR2,NPARM,IRND,KRND,LPRINT
      INTEGER  IFXP(NPTOT)
      REAL*8  PV(NPMAX),PU(NPMAX),CNST,CRND,XRND,FCT,YY,UNC
c
c** Loop over & round off the parameters # NPAR1 to NPAR2 
      IF(LPRINT.GE.2) WRITE(6,602)  NPAR2-NPAR1+1,NPTOT,NPAR1,NPAR2
      UNC= 99.d99
      DO  NPARM= NPAR1, NPAR2
          IF(PU(NPARM).LT.UNC) UNC= PU(NPARM)
          ENDDO
      DO  NPARM= NPAR1, NPAR2
c** First ... fiddle with log's to perform the rounding
          XRND= DLOG10(UNC)
          IRND= INT(XRND)
          IF(XRND.GT.0) IRND=IRND+1
          IRND= IRND- IROUND
          FCT= 10.D0**IRND
          CNST= PV(NPARM)
          YY= CNST
          CRND= PV(NPARM)/FCT
          XRND= 0.d0
c ... if rounding goes past REAL*8 precision, retain unrounded constant
          IF(DABS(CRND).GE.1.D+16) THEN
              WRITE(6,600) IROUND,NPARM
               RETURN
               ENDIF
          IF(DABS(CRND).GE.1.D+8) THEN
c ... to avoid problems from overflow of I*4 integers ...
              KRND= NINT(CRND/1.D+8)
              XRND= KRND*1.D+8
              CRND= CRND-XRND
              XRND= XRND*FCT
              END IF
          IRND= NINT(CRND)
          CNST= IRND*FCT+ XRND
          PV(NPARM) = CNST
          IFXP(NPARM)= 1
          IF(LPRINT.GE.2) WRITE(6,604) NPARM,YY,PV(NPARM)
  604 FORMAT(5x,'Round parameter #',i4,' from',G20.12,'  to',G20.12)
          ENDDO
          NPARM= NPARM- 1
      RETURN
  600 FORMAT(' =',39('==')/' Caution:',i3,'-digit rounding of parameter-
     1',i2,' would exceed (assumed) REAL*8'/' ********   precision overf
     2low at 1.D+16, so keep unrounded constant')
  602 FORMAT(' Rounding off ',i5,' of the ',i5,' parameters #:',i5,
     1 ' to',i5)
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

