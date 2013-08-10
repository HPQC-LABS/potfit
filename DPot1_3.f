
c***********************************************************************
        PROGRAM DPotFit
c***********************************************************************
c** Program "D(iatomic)Pot(ential)Fit" (DPotFit) for performing least-
c   squares fits of diatomic spectral data to molecular potential
c   energy functions for one or multiple electronic states.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++   COPYRIGHT 2006-2012  by R.J. Le Roy, J.Y. Seto and Y. Huang   +++
c   Dept. of Chemistry, Univ. of Waterloo, Waterloo, Ontario, Canada   +
c    This software may not be sold or any other commercial use made    +
c      of it without the express written permission of the authors.    +
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c                   Version of 18 November 2012
c             (just after removal of RREFad & RREFna & RREFw)
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
      INCLUDE 'BLKPARAM.h'
      INCLUDE 'BLKBOB.h'
      INCLUDE 'BLKCOUNT.h'
c-----------------------------------------------------------------------
      CHARACTER*40 DATAFILE,MAKEPRED
      CHARACTER*24 WRITFILE,TVNAME(NPARMX)
      CHARACTER*27 FN4,FN6,FN7,FN8,FN10,FN11,FN12,FN13,FN14,FN15,FN16,
     1   FN20,FN22,FN30
cc   1                                                     ,FN32
      INTEGER*4 lnblnk
      INTEGER I,J,ISTATE,ISOT,CHARGE,IPV,MKPRED,PRINP,
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
c** Set type statements for unused MASSES variables.
c
      CHARACTER*2  CATOM
      INTEGER GELGS(2,NISTPMX),GNS(2,NISTPMX)
      REAL*8 ABUND(2,NISTPMX)
c------------------------------------------------------------------------
      REAL*8 RDIST,VDIST,BETADIST
c
      REAL*8 Vsr(NPNTMX,NSTATEMX),Bsr(NPNTMX,NSTATEMX)
      INTEGER nPointSR(NSTATEMX)
      COMMON /VsrBLK/Vsr,Bsr,nPointSR
c
      REAL*8 Plr(NPNTMX,NSTATEMX),Blr(NPNTMX,NSTATEMX)
      INTEGER nPointLR(NSTATEMX)
      COMMON /PlrBLK/Plr,Blr,nPointLR
c**************************
      REAL*8 RMINT,RMAXT,RHT
      INTEGER NCNN
      DATA ZMASE /5.4857990945D-04/
      DATA MAKEPRED/'MAKEPRED'/
c=======================================================================
      SLABL(-5)='VS'
      SLABL(-4)='VR'
      SLABL(-3)='VV'
      SLABL(-2)='WI'
      SLABL(-1)='PA'
      SLABL(0)='FS'
c** uncertainties for data involving Quasibound level increased  
c   by  Fqb*width to DSQRT{u(_i;exp)**2 + (Fqb*width)**2}
      Fqb= 0.20d0
c=======================================================================
c** FSsame > 0  checks all FS to find those with a common (v',J',isot)
c   and the fit will use a single upper-state energy, instead of a
c   separate one for each series.
      FSsame= 0
c%%   FSsame= 1
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
      WRITE(FN6,'(2A)') WRITFILE(1:lnblnk(WRITFILE)),'.6'
      WRITE(FN7,'(2A)') WRITFILE(1:lnblnk(WRITFILE)),'.7'
      WRITE(FN8,'(2A)') WRITFILE(1:lnblnk(WRITFILE)),'.8'
      WRITE(FN20,'(2A)') WRITFILE(1:lnblnk(WRITFILE)),'.20'
      WRITE(FN22,'(2A)') WRITFILE(1:lnblnk(WRITFILE)),'.22'
      WRITE(FN30,'(2A)') WRITFILE(1:lnblnk(WRITFILE)),'.30'
      OPEN(UNIT=6,FILE=FN6)
      OPEN(UNIT=7,FILE=FN7)
      OPEN(UNIT=8,FILE=FN8)
      OPEN(UNIT=20,FILE=FN20)
      OPEN(UNIT=22,FILE=FN22)
      OPEN(UNIT=30,FILE=FN30)
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
cc        WRITE(FN32,'(2A)') WRITFILE(1:lnblnk(WRITFILE)),'.32'
cc        OPEN(UNIT=32, FILE= FN32)
cc        ENDIF
c!!!  
      I= 999
      WRITE(6,601) NISTP
      DO  ISOT= 1,NISTP
c
c** Read the mass numbers of the atoms in each of the isotopomers
c
c  MN(i,ISOT)  is the mass number for atom with atomic number AN(i)
c       [NOTE: be sure order of MN values consistent with that of AN's].
c       Choosing it .ne. value for some known isotope if that species
c       causes the average atomic mass to be used.
c=======================================================================
          READ(5,*) MN(1,ISOT), MN(2,ISOT)
c=======================================================================
          I= MIN(I,MN(1,ISOT),MN(2,ISOT))
          CALL MASSES(AN(1),MN(1,ISOT),CATOM,GELGS(1,ISOT),
     1                     GNS(1,ISOT),ZMASS(1,ISOT),ABUND(1,ISOT))
          IF(ISOT.EQ.1) NAME(1)= CATOM
          CALL MASSES(AN(2),MN(2,ISOT),CATOM,GELGS(2,ISOT),
     1                     GNS(2,ISOT),ZMASS(2,ISOT),ABUND(2,ISOT))
          IF(ISOT.EQ.1) NAME(2)= CATOM
          IF((AN(1).EQ.1).AND.(MN(1,ISOT).GT.3)) MN(1,ISOT)=MN(1,ISOT)-3
          IF((AN(2).EQ.1).AND.(MN(2,ISOT).GT.3)) MN(2,ISOT)=MN(2,ISOT)-3
          ZMASS(3,ISOT)= (ZMASS(1,ISOT)*ZMASS(2,ISOT))/
     1                    (ZMASS(1,ISOT)+ZMASS(2,ISOT)-CHARGE*ZMASE)
          WRITE(6,602) NAME(1),MN(1,ISOT),NAME(2),MN(2,ISOT),
     1                                          (ZMASS(J,ISOT),J=1,3)
          RSQMU(ISOT)= DSQRT(ZMASS(3,1)/ZMASS(3,ISOT))
          ENDDO
c... end of loop over isotopologues ....................................
      IF(I.EQ.0) WRITE(6,603)
      IF(CHARGE.NE.0) WRITE(6,597) CHARGE
      WRITE(6,599) DATAFILE,Fqb
cc    WRITE(6,*)
      IF(AN(1).EQ.AN(2)) WRITE(6,604)
  599 FORMAT(/' Use experimental data input file:  ',a30/' Uncertainties
     1 for transitions involving quasibound levels modified to):'/10x,
     2  'SQRT{(u(i;exp)**2 + (',f5.2,'*width)**2}')  
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
c  IOMEG(s) .EQ. -1  if it indicates a doublet SIGMA electronic state
c                 [other spin multiplets not yet coded]
c  IOMEG(s) .EQ. -2  indicated that the centrifugal potential strength 
c                    is to be scaled by [J(J+1) + 2]              
c  V(MIN/MAX)(s) Neglect data for electronic state vibrational levels
c                outside the range  VMIN  to  VMAX.
c  JTRUNC(s)     data with J > JTRUNC are not included in the fit.
c  EFSEL(s)  allows a user to consider data for:
c          * ONLY the e-parity levels of this state, if EFSEL > 0
c          * ONLY the f-parity levels of this state, if EFSEL < 0
c          * BOTH e- and f-parity levels of thsi state, if EFSEL = 0
c=======================================================================
          READ(5,*) SLABL(ISTATE), IOMEG(ISTATE), VMIN(ISTATE,1), 
     1                   VMAX(ISTATE,1), JTRUNC(ISTATE), EFSEL(ISTATE)
c======================================================================
          IF(NISTP.GT.1) THEN
              DO  ISOT= 2, NISTP
                  VMIN(ISTATE,ISOT)= VMIN(ISTATE,1)
                  VMAX(ISTATE,ISOT)= INT((VMAX(ISTATE,1)+1.0d0)/
     1                                              RSQMU(ISOT)-0.5d0)
                  ENDDO
              IF(VMAX(ISTATE,1).LT.0) 
c** If desired, read separate upper bound level for each isotopologue
c=======================================================================
     1            READ(5,*) (VMAX(ISTATE,ISOT), ISOT= 1, NISTP)
c=======================================================================
              ENDIF
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          CALL READPOT(ISTATE,SLABL,OSEL)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** These statements construct and define the names of output files
c   associated with WRITE's to channels 6-10 used by the program.
          IF(OSEL(ISTATE).GT.0) THEN
              WRITE(FN10,'(2A)') WRITFILE(1:lnblnk(WRITFILE)),'.10'
              WRITE(FN11,'(2A)') WRITFILE(1:lnblnk(WRITFILE)),'.11'
              OPEN(UNIT=10,FILE=FN10)
              OPEN(UNIT=11,FILE=FN11)
              IF(NUA(ISTATE).GE.0) THEN
                  WRITE(FN12,'(2A)') WRITFILE(1:lnblnk(WRITFILE)),'.12'
                  OPEN(UNIT=12, FILE=FN12)
                  ENDIF
              IF(NUB(ISTATE).GE.0) THEN
                  WRITE(FN13,'(2A)') WRITFILE(1:lnblnk(WRITFILE)),'.13'
                  OPEN(UNIT=13,FILE=FN13)
                  ENDIF
              IF(NTA(ISTATE).GE.0) THEN
                  WRITE(FN14,'(2A)') WRITFILE(1:lnblnk(WRITFILE)),'.14'
                  OPEN(UNIT=14,FILE=FN14)
                  ENDIF
              IF(NTB(ISTATE).GE.0) THEN
                  WRITE(FN15,'(2A)') WRITFILE(1:lnblnk(WRITFILE)),'.15'
                  OPEN(UNIT=15,FILE=FN15)
                  ENDIF
              IF(NwCFT(ISTATE).GE.0) THEN
                  WRITE(FN16,'(2A)') WRITFILE(1:lnblnk(WRITFILE)),'.16'
                  OPEN(UNIT=16,FILE=FN16)
                  ENDIF
              ENDIF
          PASok(ISTATE)= 1
          IF(PSEL(ISTATE).EQ.6) PASok(ISTATE)= 0
c** Call VGEN to prepare parameters for output in WRITEPOT
          IF(PSEL(ISTATE).EQ.2) CALL VGEN(ISTATE,1.0d0,VDIST,BETADIST,0)
          ENDDO
c** Now write summary of the initial potential parameters for each state
      CALL WRITEPOT(1,SLABL,NAME,DECM,PV,PU,PS)
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
          IF((PSEL(ISTATE).EQ.0).OR.(PSEL(ISTATE).EQ.-2)) GOTO 90
          IF(PSEL(ISTATE).EQ.-1) THEN
c... When using band constants for this state ... count them and label 
c    first and last for each level of each isotopologue ...
              DO  ISOT= 1, NISTP
                  DO  I= VMIN(ISTATE,ISOT),VMAX(ISTATE,ISOT)
                      IF(NBC(I,ISOT,ISTATE).GT.0) THEN
                          BCPARI(I,ISOT,ISTATE)= IPV+1
                          DO  J= 1,NBC(I,ISOT,ISTATE)
                              IPV= IPV+1
                              IFXPV(IPV)= 0
                              PV(IPV)= 0.d0
                              PU(IPV)= 0.d0
                              ENDDO
                          BCPARF(I,ISOT,ISTATE)= IPV
                          ENDIF
                      ENDDO
                  ENDDO
              GOTO 90
              ENDIF
c... For all types of fitted potential, count  Re
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
          IF(PSEL(ISTATE).LT.6) THEN
              IFXPV(IPV)= IFXDe(ISTATE)
              IPV= IPV+ 1
              POTPARF(ISTATE)= IPV
              IFXPV(IPV)= IFXRe(ISTATE)
              ENDIF
          IF((PSEL(ISTATE).GE.2).and.(PSEL(ISTATE).LE.5)) THEN
c... For MLJ/MLR family, count long-range parameters: first count Cm's
              DO  J= 1,NCMM(ISTATE)
c... additional Aubert-Frecon{3,6,6,8,8} parameters included in this count
                  IPV= IPV+ 1
                  POTPARF(ISTATE)= IPV
                  IFXPV(IPV)= IFXCm(J,ISTATE)
                  IF(IFXPV(IPV).GT.1) THEN
c!!! If constraining one or more Cm to be fixed at same value as for
c  an earlier (smaller ISTATE) state ....
                      WRITE(6,610) IPV,IFXPV(IPV)
  610 FORMAT('** Constrain PV(',i3,') = PV(',I3,')  in the fits')
                      ENDIF
                  ENDDO
              ENDIF
c... Now count [exponent] \beta_i expansion coefficients
          J=0
          IF(NSpow(ISTATE).LT.0) J=1
          DO I= J,MAX0(NSpow(ISTATE),NLpow(ISTATE))
              IPV= IPV+ 1
              POTPARF(ISTATE)= IPV
              IFXPV(IPV)= IFXBETA(I,ISTATE)
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
c..... end of parameter count/label loop!
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
c** Call subroutine to input experimental data in specified band-by-band
c   format, and do bookkeeping to characterize amounts of data of each
c   type.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IF(MKPRED.LE.0) OPEN(UNIT= 4, STATUS= 'OLD', FILE= DATAFILE)
c** when COMMON blocks check out ... introduce MKPRED option ......
      IF(MKPRED.GT.0) THEN
          WRITE(FN4,'(2A)') WRITFILE(1:lnblnk(WRITFILE)),'.4'
          OPEN(UNIT= 4, FILE= FN4)
          IF(UCUTOFF.LT.1.d0) UCUTOFF= 1.d0
          CALL MKPREDICT(NSTATES,NDAT)
          REWIND(4)
          ENDIF
      CALL READATA(NSTATES,PASok,UCUTOFF,JTRUNC,EFSEL,VMIN,VMAX,NDAT,
     1                                                 NOWIDTHS,PRINP)
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
c** Set the energy convergence criterion to be 1/100th of the smallest
c   experimental uncertainty. [UCUTOFF reset by READATA to that min. unc.]
      DO  ISTATE=1,NSTATES
          EPS(ISTATE)= DMIN1(UCUTOFF/100.0d0,1.d-06)
          WRITE(6,638) SLABL(ISTATE), EPS(ISTATE)
c** Initialize the dissociation energy ????
          DECM(ISTATE)= 0.0d0
          ENDDO
      IF(IROUND.NE.0) WRITE(6,685) IABS(IROUND)
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

      flush(6)

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Now Generate internal NLLSSRR variables {PV} from the external ones 
      CALL MAPPAR(NISTP,PV,0)
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
              IF(IFXFS(I-TOTPOTPAR).GT.0) IFXPV(I)= 1
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
              IF(IFXFS(I-TOTPOTPAR).GT.0) THEN

              write(6,888) I-TOTPOTPAR,TVALUE(I-TOTPOTPAR),
     1                 IFXFS(I-TOTPOTPAR),TVALUE(IFXFS(I-TOTPOTPAR))
  888 format(/' Following FS fit reset  T(',i3,')=',f12.4,
     1  '  equal   T(',I3,')=', F12.4)

                    TVALUE(I-TOTPOTPAR)= TVALUE(IFXFS(I-TOTPOTPAR))
                    IFXPV(I)= 1
                    ENDIF
              ENDDO
          NFPAR= NFPAR+ NFSTOT+ NTVALL(0)
          CALL MAPPAR(NISTP,PV,0)
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
c** Perform group rounding of all term values and/or fluorescence series
c   origins in single step
          IF((NFSTOT.GT.0).OR.(NTVALL(0).GT.0)) THEN
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
c ... renormalize DSE relative value for to ALL fitting parameters free
          DSE= DSE*DSQRT(DFLOAT(COUNTOT- (NFPAR- NFSTOT- NTVALL(0)))/
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
            ELSEIF(PSEL(ISTATE).GT.0) THEN
c** Calculation of the uncertainties for Te for each potential require 
c   elements from the correlation matrix.
              IF((IFXDE(1).LE.0).AND.(IFXDE(ISTATE).LE.0)) THEN
                  DECM(ISTATE)= CM(1,POTPARI(ISTATE))
                ELSE
                  DECM(ISTATE)= 0.0d0
                ENDIF
            ENDIF
          ENDDO
  690 FORMAT(/,1X,34('=='))
  691 FORMAT(' For fit of',I5,'  parameters to',I6,
     1  '  transitions,  DSE=',G15.8)
  692 FORMAT('  The following',I5,' Fluorescence Series Origins were det
     1ermined'/1x,30('--')/"  ( v', J', p'; ISTP)",4x,'T(value)',4x,
     2  'Uncertainty  Sensitivity'/1x,30('--'))
cc694 FORMAT(3X,'(',I3,',',I3,',',SP,I3,SS,';',I2,')',1X,1PD19.10,
  694 FORMAT(2X,'(',I4,',',I3,',',SP,I3,SS,';',I2,')',1X,F15.6,
     1  1PD11.1,D12.1)
  696 FORMAT('  State ',A2,' represented by the',I5,' individual term va
     1lues:'/1x,33('--')/" T(es: v', J', p';IS)  #dat",4x,'T(value)',4x,
     2  'Uncertainty  Sensitivity'/1x,33('--'))
  698 FORMAT(2X,A24,1PD19.10,D11.1,D12.1)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Call subroutine MAPPAR to convert internal NLLSSRR parameter array
c   back into external (logical) variable system.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      CALL MAPPAR(NISTP,PV,1)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Call subroutine VGEN to generate the potential function from the
c   final calculated converged parameters.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      DO  ISTATE= 1,NSTATES
          IF(PSEL(ISTATE).GT.0) THEN 
              RMINT= RD(1,ISTATE)
              RMAXT= RD(NDATPT(ISTATE),ISTATE)
              RHT= RD(2,ISTATE) - RD(1,ISTATE) 
              nPointSR(ISTATE)= RMINT/RHT 
              IF(OSEL(ISTATE) .NE. 0) THEN
                  IF(RMAXT .GT. 100.0) THEN
                      nPointLR(ISTATE)= 0
                    ELSE
                      nPointLR(ISTATE)= (100.0-RMAXT)
     &                                    /(RHT*OSEL(ISTATE))
                    ENDIF
                  ENDIF
              CALL VGEN(ISTATE,-1.0d0,VDIST,BETADIST,0)
              IF(OSEL(ISTATE).NE.0) THEN
                  DO  I= 1,nPointSR(ISTATE),OSEL(ISTATE)
c ... generate potential & exponent values in inner extrapolation region
                      RDIST= RHT*DBLE(I-1) 
                      CALL VGEN(ISTATE,RDIST,VDIST,BETADIST,0)
                      Vsr(I,ISTATE)= VDIST
                      Bsr(I,ISTATE)= BETADIST
                      ENDDO
                  DO  I= 1,nPointLR(ISTATE)
c ... generate potential & exponent values in outer extrapolation region
                      RDIST= RMAXT + RHT*DBLE(I*OSEL(ISTATE))
                      CALL VGEN(ISTATE,RDIST,VDIST,BETADIST,0)
                      Plr(I,ISTATE)= VDIST
                      Blr(I,ISTATE)= BETADIST
                      ENDDO
                  ENDIF
              ENDIF
          ENDDO
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Call subroutine to print out a summary of the converged and fixed
c   values to standard output (channel-6).
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      CALL WRITEPOT(2,SLABL,NAME,DECM,PV,PU,PS)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** If chosen, output file(s) will be created for the export of the
c   generated functions: V, BETAFX, UAR/UBR, or TAR/TBR and their
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
                  DO  I= 1,NDATPT(NSTATES),OSEL(ISTATE)
                      WRITE(18,900) RD(I,ISTATE),VPOT(I,NSTATES)
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
c                          COPYRIGHT 2005-2010
c** By R.J. Le Roy (with assistance from G.T. Kraemer & J.Y. Seto).
c           Last modified  25 July 2010 (added proton, d, t)
c***********************************************************************
      REAL*8 zm(123,0:10),mass,ab(123,10),abund
      INTEGER i,ian,imn,gel(123),nmn(123),mn(123,10),ns2(123,10),
     1        gelgs, gns
      CHARACTER*2 NAME,AT(123)
c
      DATA  at(1),gel(1),nmn(1),(mn(1,i),i=1,6)/' H',2,6,1,2,3,4,5,6/
      DATA  (zm(1,i),i=0,6)/1.00794d0, 1.00782503207d0, 2.0141017778d0,
     1                 3.0160492777d0,1.00727646677d0,2.013553212724d0,
     2                 3.0155007134d0/
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
      DATA (zm(64,i),i=0,7)/157.25d0, 151.9197910d0, 153.92086560d0,
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
     1    170.936323580d0, 171.9363815d0, 172.9382108d0, 173.9388621d0,
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
          IF(IMN.EQ.4) NAME=' p'
          IF(IMN.EQ.5) NAME=' d'
          IF(IMN.EQ.6) NAME=' t'
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
c             ********* Version of 14 July 2011 *********
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++  COPYRIGHT 1997-2011 by  Robert J. Le Roy & Dominique R.T. Appadoo +
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
c (ix) 2'nd virial coefficient data (also only for DPotFit applications)
c (x)  Potential function value from some other source (e.g., ab initio
c      energy high on the repusive wall.
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
c     VMIN(s,ISOT)  to  VMAX(s,ISOT),  AND  NEGLECT any data for which 
c     read-in uncertainty is  > UCUTOFF (cm-1).  EFSEL(s) > 0 causes 
c     f-parity levels to be neglected, EFSEL(s) < 0 omits e-parity levels
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
     1 VMIN(NSTATEMX,NISTPMX),VMAX(NSTATEMX,NISTPMX),
     2 NDAT(0:NVIBMX,NISTPMX,NSTATEMX),PASok(NSTATES)
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
              WRITE(6,607) SLABL(ISTATE),JTRUNC(ISTATE)
            ELSE
              WRITE(6,605) SLABL(ISTATE),-JTRUNC(ISTATE)
            ENDIF 
          WRITE(6,611) (VMIN(ISTATE,ISOT),VMAX(ISTATE,ISOT),ISOT,
     1                                                  ISOT= 1,NISTP)
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
              NVVPP(ISOT,ISTATE)= 0
              NWIDTH(ISOT,ISTATE)= 0
              NEBPAS(ISOT,ISTATE)= 0
              NVIRIAL(ISOT,ISTATE)= 0
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
            GOTO 40
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
c        energies [D-E(v,J,p)] for a given state.  Data Labels as for 'FS'
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
c   LABLP = 'VR' identifies this data group/band as a set of virial 
c           coefficients for electronic state LABLPP.  In this case, 
c           parameters VP, VPP are dummy variables.
c** STOP reading when run out of bands OR when read-in VPP is negative   
c-----------------------------------------------------------------------
      READ(4,*,END=40) VP(IBAND), VPP(IBAND), LABLP, LABLPP, MN1, MN2
c-----------------------------------------------------------------------
      IF(VP(IBAND).LT.0) GO TO 40
      IEP(IBAND)= -99
      IEPP(IBAND)= -99
      DO  I= -5, NSTATES
          IF(LABLP.EQ.SLABL(I)) IEP(IBAND)= I
          IF(LABLPP.EQ.SLABL(I)) IEPP(IBAND)= I
          ENDDO
c** Check that this isotopomer is one of those chosen to be fitted ...
      ISOT= 0
      DO  I= 1,NISTP
          IF((MN1.EQ.MN(1,I)).AND.(MN2.EQ.MN(2,I))) ISOT= I
          ENDDO
      ISTP(IBAND)= ISOT
      ESP= IEP(IBAND)
      ESPP= IEPP(IBAND)
      IF(IEP(IBAND).EQ.-3) THEN
c** For case in which the 'data' are potential function value(s)...
          COUNT= COUNT+ 1
          IFIRST(IBAND)= COUNT
c...  TEMP(i)= r(i) ;  FREQ(i)= V(r(i))  ; UFREQ= unc{V(r(i))}
c----------------------------------------------------------------------
   12     READ(4,*) TEMP(COUNT),FREQ(COUNT),UFREQ(COUNT)
c----------------------------------------------------------------------
          YUNC(COUNT)= UFREQ(COUNT)
          IF(TEMP(COUNT).GT.0.d0) THEN
c ... a negative input temperature implies end of virial data set
              IF(ISOT.LE.0) GOTO 12
c ... if this isotope not considered in the fit, ignore this datum
              IB(COUNT)= IBAND
              COUNT= COUNT+1
              GOTO 12
            ELSE
              GOTO 18
            ENDIF
          ENDIF
      IF(IEP(IBAND).EQ.-4) THEN
c** For case in which the data are virial coefficients
          COUNT= COUNT+ 1
          IFIRST(IBAND)= COUNT
c... TEMP(i)= temperature[K],  FREQ(i)= Bvir(i),  UFREQ(i)= unc{Bvir(i)}
c----------------------------------------------------------------------
   14     READ(4,*) TEMP(COUNT),FREQ(COUNT),UFREQ(COUNT)
c----------------------------------------------------------------------
          YUNC(COUNT)= UFREQ(COUNT)
          IF(TEMP(COUNT).GT.0.d0) THEN
c ... a negative input temperature implies end of virial data set
              IF(ISOT.LE.0) GOTO 14
c ... if this isotope not considered in the fit, ignore this datum
              IB(COUNT)= IBAND
              COUNT= COUNT+1
              GOTO 14
            ELSE
              GOTO 18
            ENDIF
          ENDIF
c... now ... for the case of spectroscopic data ...
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
      VMAXespp= 0
      VMINespp= 0
      VMAXesp= 0
      VMINesp= 0
      IF((ESPP.GT.0).AND.(ISOT.GT.0)) THEN
          VMAXespp= VMAX(ESPP,ISOT)
          VMINespp= VMIN(ESPP,ISOT)
          JTRUNCespp= JTRUNC(ESPP)
          IF(ISOT.GT.1) THEN
              JTRUNCespp= INT(JTRUNC(ESPP)/RSQMU(ISOT))
              ENDIF
          VMAXesp= VMAX(ESPP,ISOT)
          ENDIF
      IF((ESP.GT.0).AND.(ISOT.GT.0)) THEN
          VMAXesp= VMAX(ESP,ISOT)
          VMINesp= VMIN(ESP,ISOT)
          JTRUNCesp= JTRUNC(ESP)
          IF(ISOT.GT.1) THEN
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
      YUNC(COUNT)= UFREQ(COUNT)
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
c?? RJL What was the purpose of UMIN ?  a check on EPS ?
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
        ELSEIF(ESP.LE.-3) THEN
c... and for potential function values as data ...
          IF(TEMP(COUNT).LT.JMIN(IBAND)) JMIN(IBAND)= TEMP(COUNT) 
          IF(TEMP(COUNT).GT.JMAX(IBAND)) JMAX(IBAND)= TEMP(COUNT)
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
c** Define counters to label which f.s. is associated with band IBAND 
c ...  FSBAND(j)  is the absolute band number for the  j'th FS
c ...  NDF(IBAND)  if the FS number associated with band  IBAND
          FSBAND(NFSTOT)= IBAND
          NFS(IBAND)= NFSTOT
          NBANDFS(ISOT,ESPP)= NBANDFS(ISOT,ESPP)+ 1
          NBND= NBANDFS(ISOT,ESPP)
          NTRANSFS(ISOT,ESPP)= NTRANSFS(ISOT,ESPP)+NTRANS
c ... and then set up labels/ranges/properties for each band
          YPR(ISOT,ESPP,1,1,NBND)= VP(IBAND)
          YPR(ISOT,ESPP,1,2,NBND)= VPP(IBAND)
          YPR(ISOT,ESPP,1,3,NBND)= NTRANS
          YPR(ISOT,ESPP,1,4,NBND)= IBAND
          YPR(ISOT,ESPP,1,5,NBND)= JMIN(IBAND)
          YPR(ISOT,ESPP,1,6,NBND)= JMAX(IBAND)
          IFXFS(NFSTOT)= 0
          IF((NFSTOT.GT.1).AND.(FSsame.GT.0)) THEN
c** Finally - If desired (FSsame > 0) check to see if this band has the
c   same upper state as an FS for this isotopologue encountered earlier,
c   and if so (try) to relabel origin accordingly ...
              DO  I= 1, NFSTOT-1
                  IF((VP(IBAND).EQ.VP(FSBAND(I))).AND.
     1                   (VPP(IBAND).EQ.VPP(FSBAND(I))).AND.
     2                          (ISTP(IBAND).EQ.ISTP(FSBAND(I)))) THEN
c ... fix origin for this FS band to be the same as that for FS band I
                      IFXFS(NFSTOT)= I
                      WRITE(6,654) VP(IBAND),VPP(IBAND),ISTP(IBAND),
     1                                                        NFSTOT,I
  654 FORMAT(" NOTE that  FS(v'=",I4,", J'=",I3,", ISOT=",I2,") #",I4,
     1  " has same origin as  FS #",I4)
                      GOTO 20
                      ENDIF
                  ENDDO
   20         CONTINUE
              ENDIF
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
c  this isotopomer of this electronic state.  Expect to find all 
c  potential fx. values, virial coeficients, Tunneling Widths or PAS 
c  binding energies in a single group.
      IF(ESP.EQ.-4) THEN
c** Data are Virial Coefficients for electronic state IEPP= ESPP
          NVIRIAL(ISOT,ESPP)= NTRANS
          YPR(ISOT,ESPP,8,3,1)= NTRANS
          YPR(ISOT,ESPP,8,4,1)= IBAND
          YPR(ISOT,ESPP,8,5,1)= JMIN(IBAND)
          YPR(ISOT,ESPP,8,6,1)= JMAX(IBAND)
          ENDIF
c
      IF(ESP.EQ.-3) THEN
c** Data are not transition energies, but rather values of the potential
c  function at particular distances for electronic state s=IEPP  
          WRITE(6,612) ESPP,ISOT
          NVVPP(ISOT,ESPP)= NTRANS
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
c
c** Now return to read the next band
      GOTO 10
c========================================================================
c** Now, write a summary of the input data to the output file
   40 COUNTOT= COUNT
      NBANDTOT= 0
      DO  I= 1,NISTP
          NBANDTOT= NBANDTOT+ NBANDS(I)
          ENDDO
      ISOT= 1
      UCUTOFF= UMIN
      IF(FSOMIT.GT.0) WRITE(6,650) FSOMIT
      IF(PRINP.LE.0) RETURN
c** Write a summary of the data, one isotopomer at a time.
   26 WRITE(6,602) NBANDS(ISOT), (NAME(I),MN(I,ISOT),I=1,2)
c
      DO 50 ISTATE= 1,NSTATES
c ... For internal use, may wish to update VMAX(ISTATE,ISOT) to actual 
c  highest v in the data set for this state. ** Reactivate as needed.
c      VMAX(ISTATE,ISOT)= VMX(ISTATE)
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
      IF(NVVPP(ISOT,ISTATE).GT.0)THEN
c** Book-keeping for potential function values as data ....
          WRITE(6,618) NVVPP(ISOT,ISTATE),SLABL(ISTATE),(NAME(I),
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
      IF(NVIRIAL(ISOT,ISTATE).GT.0) THEN
c** Book-keeping for Virial data
          IBB= YPR(ISOT,ESPP,8,4,1)
          WRITE(6,642) NVIRIAL(ISOT,ISTATE), SLABL(ISTATE), 
     1                                      (NAME(I),MN(I,ISOT),I=1,2)
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
  603 FORMAT(/' Neglect data with:  Uncertainties > UCUTOFF=',1PD10.2,
     1  ' (cm-1)')
  605 FORMAT(7x,'and State ',A2,' data with  J < JTRUNC=',I4)
  607 FORMAT(7x,'and State ',A2,' data with  J > JTRUNC=',I4)
  611 FORMAT(29x,'or  v  outside range',i3,'  to',i4,'   for   ISOT=',
     1  i2:)
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
  612 FORMAT(/" NOTE that read-in potential fx. values for   ISTATE=",
     1   i2,'   ISOT=',i2/32x,' must be input as a single "band" or data
     1 group')
  614 FORMAT(1x,38('==')/I5,' Fluorescence transitions into State ',
     1 A2,2x,A2,'(',I3,')-',A2,'(',I3,')  in',i5,' series'/
     2 1x,38('==')/"   v'  j' p' ",'#data  v"min  v"max  Avge.Unc.  Max.
     3Unc.'/1x,51('-'))
  616 FORMAT(2I4,A3,I6,2I7,1x,1P2D10.1)
  618 FORMAT(1x,65('=')/1x,I3,' State ',A2,1x,A2,'(',I3,')-',A2,'(',I3,
     1 ') potential fx values treated as independent data'/1x,24('--')/
     2 '  #values   r(min)  r(max)  Avge.Unc.   Max.Unc.'/1x,24('--'))
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
  642 FORMAT(1x,70('=')/I4,' State ',A2,1x,A2,'(',I3,')-',A2,
     1 '(',I3,') Virial coefficients included in data set' )
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
      INTEGER I,J,P,IBAND,ISOT,ISTATE,NPTOT,LOWEST,
     1 VMAX(NSTATEMX,NISTPMX),NLV(0:NVIBMX,0:NVIBMX,-1:1),
     2 NTVS(NSTATEMX,NISTPMX),NTVALL(0:NSTATEMX)
      CHARACTER*24 TVNAME(NPARMX)
c=======================================================================
      WRITE(6,600) SLABL(ISTATE) 
      LOWEST= 1
      IF(ISTATE.GT.1) LOWEST= 0
      NTVALL(ISTATE)= 0
      DO  ISOT= 1, NISTP
c** First ... zero transition counter array for this isotopomer
          DO  I= 0, VMAX(ISTATE,ISOT)
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
                      IF(IEPP(IBAND).EQ.ISTATE) THEN
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
          DO  I= 0, VMAX(ISTATE,ISOT)
              DO  J= 0, NVIBMX
                  DO  P= -1,1
                      IF(NLV(I,J,P).GT.0) THEN
c!! For ParFit ONLY!!     IF(LOWEST.EQ.1) THEN
c!! If using term values for `lowest' state (defined as the first state
c!!considered), its lowest observed level for isotopologue-1 defines the
c!! absolute energy zero
c!!                           WRITE(6,606) I,J,P,ISOT,SLABL(ISTATE)
c!!                           LOWEST= 0
c!!                           NLV(I,J,P)= 0
c!!                           GOTO 20
c!!                           ENDIF
                          NPTOT= NPTOT+ 1
                          NTV(ISTATE,ISOT)= NTV(ISTATE,ISOT)+ 1
                          IF(NLV(I,J,P).EQ.1) NTVS(ISTATE,ISOT)=
     1                                           NTVS(ISTATE,ISOT) +1 
                          REWIND(30)
                          WRITE(30,700) SLABL(ISTATE),I,J,P,ISOT,
     1                                                      NLV(I,J,P)
ccc                       WRITE(31,700) SLABL(ISTATE),I,J,P,ISOT,
ccc  1                                                      NLV(I,J,P)
                          REWIND(30)
                          READ(30,*) TVNAME(NPTOT)
c ... reset NLV(v,J,p) as the parameter index for that term value
                          NLV(I,J,P)= NPTOT
                          ENDIF
   20                 CONTINUE
                      ENDDO
                  ENDDO
              ENDDO
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
     1                          TVLW(I)= NLV(VPP(IBAND),JPP(I),EFPP(I))
                      ENDDO
                  ENDIF
              ENDDO
          WRITE(6,608) NAME(1),MN(1,ISOT),NAME(2),MN(2,ISOT),
     1                              NTV(ISTATE,ISOT),NTVS(ISTATE,ISOT)
          NTVALL(ISTATE)= NTVALL(ISTATE)+ NTV(ISTATE,ISOT)
          ENDDO
c
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
  608 FORMAT(' For ',A2,'(',i3,')-',A2,'(',I3,')  fit to',i5,
     1 ' T(v,J,p) term values,'/20x,'of which',i5,' are involved in only
     2 one transition')
  700 FORMAT("'",'T(',A2,':',i3,',',i3,',',SP,i2,';',SS,i2,')',I5,"'")
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c**********************************************************************N(
      SUBROUTINE READPOT(ISTATE,SLABL,OSEL)
c**********************************************************************
c** This subroutine reads parameters that define the model potential or
c   parameter representation used for each state in the fit procedure
c   analytical molecular potentials for the direct Hamiltonian fitting
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c                   Version of 18 November 2012
c             (after removal of RREFns, RREFad & RREFw)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** On entry:
c   ISTATE  is the electronic state being fitted to
c   SLABL   is the two-character label identifying that state
c-----------------------------------------------------------------------
      INCLUDE 'arrsizes.h'
      INCLUDE 'BLKISOT.h'
      INCLUDE 'BLKPOT.h'
      INCLUDE 'BLKPARAM.h'
      INCLUDE 'BLKBOB.h'
c-----------------------------------------------------------------------
c** Type statements for input or local variables
      INTEGER I, ISTATE, IISTP, m, MMN, npow, VTST, OSEL(NSTATEMX)
      CHARACTER*2 SLABL(-5:NSTATEMX)
      CHARACTER*1 DASH
      REAL*8 ZMASE, RR(NPNTMX), VV(NPNTMX)
      DATA ZMASE /5.4857990945D-04/,DASH/'-'/
c
c** Set some defaults for parameters not common to all models ...
      IFXDe(ISTATE)= 1
      IFXRe(ISTATE)= 1
      DO  m= 1, NCMMax
          IFXCm(m,ISTATE)= 1
          ENDDO
c-----------------------------------------------------------------------
c** First choose potential model and select form of BOB representation
c  PSEL(s)  choses the type of analytical potential to be fitted to:
c         = -2 : represent each distinct observed level of this state as
c               an independent term value [an alternative to an 'FS' 
c               treatmentof transitions involving that state]
c         = -1 : represent the rotational sublevels for each v of each
c                isotopologue by Band Constant (!)
c         = 0 : Use a fixed potential defined by LEVEL's PREPOT routine
c         = 1 : Use an Expanded Morse Oscillator EMO(p) potential
c         = 2 : Use a Morse/Long-Range (MLR) Potential.
c         = 3 : Use a Double-Exponential Long-Range (DELR) Potential.
c         = 4 : Use a Tiemann/Hannover-polynomial-potential (HPP)
c         = 5 : Use a Tang-Toennies type potential 
c         = 6 : Use a Surkus Generalized Potential Energy Function (GPEF).
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
      IF((PSEL(ISTATE).EQ.-1).OR.(PSEL(ISTATE).EQ. -2)) THEN
          IF(PSEL(ISTATE).EQ.-2) WRITE(6,604) SLABL(ISTATE)
          IF(PSEL(ISTATE).EQ.-1) THEN
              WRITE(6,606) SLABL(ISTATE),(I,I=1,NISTP)
              WRITE(6,607) (DASH,I=1,NISTP)
c** If representing data for this state by fitted band constants, 
c   read in the number of band constants for each vibrational level
              DO  I= VMIN(ISTATE,1),VMAX(ISTATE,1)
c** For each isotopologue in each vibrational level, read the number of 
c  band constants to be used (fited to) to represent the data,
c=======================================================================
                  READ(5,*) VTST,(NBC(I,IISTP,ISTATE),IISTP= 1,NISTP)
c=======================================================================
                  IF(I.NE.VTST) THEN
c... Verify that band constant specification is for the correct vib level
                      WRITE(6,610) I,VTST
                      STOP
                      ENDIF
                  IF(NBC(I,IISTP,ISTATE).GT.NBCMX) 
     1                                      NBC(I,IISTP,ISTATE)= NBCMX
                  WRITE(6,608) I,(NBC(I,IISTP,ISTATE),IISTP= 1,NISTP)
                  ENDDO
              ENDIF
          NUA(ISTATE)= -1
          NUB(ISTATE)= -1
          NTA(ISTATE)= -1
          NTB(ISTATE)= -1
          RETURN
          ENDIF
c-----------------------------------------------------------------------
c** Now to read in the range and mesh for the numerical integration
c  RMIN/MAX(s)  define the range over which this potential is defined.
c  RH(s)  specifies radial mesh for numerical integration for this state
c=======================================================================
      READ(5,*) RMIN(ISTATE), RMAX(ISTATE), RH(ISTATE)
c=======================================================================
      NDATPT(ISTATE)= (RMAX(ISTATE)-RMIN(ISTATE))/RH(ISTATE)+1.0001d0
      NDATPT(ISTATE)= MIN(NPNTMX-1,NDATPT(ISTATE))
      RMAX(ISTATE)= RMIN(ISTATE) + RH(ISTATE)*DBLE(NDATPT(ISTATE)-1)
      DO  I= 1, NDATPT(ISTATE)
          RD(I,ISTATE)= RMIN(ISTATE)+ DBLE(I-1)*RH(ISTATE)
          ENDDO
      IF(PSEL(ISTATE).EQ.0) THEN
c-----------------------------------------------------------------------
c** For case of a fixed potential defined by read-in turning points, 
c  subroutine PREPOTT reads those points & generates potential array
c-----------------------------------------------------------------------
          DO  I= 1, NDATPT(ISTATE)
              RR(I)= RD(I,ISTATE)
              ENDDO
          WRITE(6,600)SLABL(ISTATE),RMIN(ISTATE),RMAX(ISTATE),RH(ISTATE)
          CALL PREPOTT(1,AN(1),AN(2),MN(1,1),MN(2,1),NDATPT(ISTATE),
     1                                             VLIM(ISTATE),RR,VV)
          DO  I= 1, NDATPT(ISTATE)
              VPOT(I,ISTATE)= VV(I)
              ENDDO
          NUA(ISTATE)= -1
          NUB(ISTATE)= -1
          NTA(ISTATE)= -1
          NTB(ISTATE)= -1
          NwCFT(ISTATE)= -1
          RETURN
          ENDIF
      IF((PSEL(ISTATE).GE.2).AND.(PSEL(ISTATE).LE.5)) THEN
c-----------------------------------------------------------------------
c** For MLR, DELR of HPP potential read number of terms NCMM in 
c   damped inverse-power long-range tail
c    uLR(R) = - SUM_{i=1}^{NCMM} Dm(R;MMLR(i) * CmVAL(i)/R**MMLR(i) 
c** If rhoAB .LE. 0.0  have NO damping functions: all  Dm(R)= 1.0
c   If rhpAB > 0.0  it is the molecule-dependent radial scaling factor
c                    of Douketis et al. [JCP 76, 3057 (1982)]
c     rhoAB =  2*rhoA*rhoB/(rhoA+rhoB)   where  rhoA  is the ionization
c                 potential ratio  (I_p^A/I_p^H)^0.66  for atom A
c
c  For rhoAB > 0.0,  IDF specifies damping s.th.  Dm(r)/r^m --> r^{IDF/2}
c                    IDSTT > 0  use Douketis et al. damping functions
c                    IDSTT .LE. 0  use Tang-Toennies damping functions
c   {NOTE:  IDF & IDSTT required only for testing - remove em later!}
c
c** IFXCm specifies whether this long-range coefficient is to be fitted
c  freely (when .LE.0), held fixed at the read-in value (when =1) or held
c  fixed at the value for another state, in which case the parameter value
c  'IFXCm(m,ISTATE)' is the no. of the parameter it is constrained to =
c** For alkali dimer A-state: Aubert-Frecon ULR(r) [PRA 55, 3458 (1997)]
c  is special case with MMLR= {3, 0, 6, 6, 8, 8, ...} 
c=======================================================================
          READ(5,*) NCMM(ISTATE), rhoAB(ISTATE), IDF(ISTATE),
     1                                                   IDSTT(ISTATE)
          DO  m= 1,NCMM(ISTATE)
              READ(5,*) MMLR(m,ISTATE), CmVAL(m,ISTATE), IFXCm(m,ISTATE)
              ENDDO
c=======================================================================
c** Note that HPP potentials have no damping
          IF(PSEL(ISTATE).EQ.4) rhoAB(ISTATE)= -1.d0
          ENDIF
c-----------------------------------------------------------------------
      IF(PSEL(ISTATE).EQ.6) THEN
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
      IF(PSEL(ISTATE).EQ.6) IFXDE(ISTATE)= 1
c=======================================================================
c** Read in parameters defining the exponent expansion \phi(r)
c  NSpow(s)  is the order of the phi_i exponent expansion for  R < Re
c  NLpow(s)  is the order of the phi_i exponent expansion for R > Re
c  nPB(s)  is the power  p  for  beta(r)  basic exponent variable
c  nQB(s)  is the power  q  for  beta(r)  exponent expansion variable
c  RREF(s)  defines the reference distance in the potential exponent 
c    expansion variable:  * for  RREF.le.0 , define parameter  RREF = Re
c      * for  RREF.gt.0 , fix parameter  RREF   at its read-in value
c=======================================================================
      READ(5,*) NSpow(ISTATE), NLpow(ISTATE), nPB(ISTATE), nQB(ISTATE),
     1                                                     RREF(ISTATE)
c=======================================================================
      IF((PSEL(ISTATE).EQ.2).OR.(PSEL(ISTATE).EQ.3)) THEN
          MMN= MMLR(NCMM(ISTATE),ISTATE)- MMLR(1,ISTATE)
          IF((NCMM(ISTATE).GT.1).AND.(nPB(ISTATE).LE.MMN)) 
     1                                    WRITE(6,628) nPB(ISTATE),MMN
          ENDIF
      IF(NSpow(ISTATE).GE.0) THEN
c** For NSpow .ge.0  the MLR exponent  \beta(yp)  is  'conventional'
c                    Huang-type constrained polynomial in  yp.
          npow= MAX(NSpow(ISTATE),NLpow(ISTATE))
          IF(PSEL(ISTATE).EQ.4) npow= npow+3
          DO  I= 0,npow
c** Read in trial initial trial parameters for exponent \phi(r)
c
c  BETA(i,s)   contains the expansion parameters defining the potential
c    for  PSEL.LE.3 : read-in values are the {npow+1} beta_i  exponent
c                   exponent expansion parameters defining the potential
c    for  PSEL = 4  :  read in the  {1+NLpow} expansion parameters plus
c                     b, `RINN, and ROUT  of the HPP form
c    for  PSEL = 6  :    read-in values are leading coefficients in 
c                Surkus' Generalized Potential Energy Function (GPEF).
c  IFXBETA(i,s) indicates whether each potential expansion coefficient
c             coefficient will be:    = 1: held fixed at read-in values.
c                                  .LE. 0: determined from fits.
c=======================================================================
              READ(5,*) BETA(I,ISTATE), IFXBETA(I,ISTATE)
c=======================================================================
              ENDDO
        ELSE
          npow= NLpow(ISTATE)
          DO  I= 1,npow
c-----------------------------------------------------------------------
c** For NSpow < 0  exponent is a natural spline function with values BETA
c       at the yp values ypBETA, and fixed to equal ypINF at  ypBETA=1
c=======================================================================
              READ(5,*)ypBETA(I,ISTATE),BETA(I,ISTATE),IFXBETA(I,ISTATE)
c=======================================================================
              ENDDO
          NLpow(ISTATE)= NLpow(ISTATE)+1
          ypBETA(NLpow(ISTATE),ISTATE)= 1.d0
          IFXBETA(NLpow(ISTATE),ISTATE)= 1
        ENDIF
      IF(PSEL(ISTATE).EQ.4) THEN
c** Constraints for Tiemann polynomial potential ....
          NSpow(ISTATE)= NLpow(ISTATE)
          nPB(ISTATE)= 1
          nQB(ISTATE)= 1
          RREF(ISTATE)= RE(ISTATE)
          ENDIF
c=======================================================================
c** Read parameters defining the BOB adiabatic radial functions
c  NUA/NUB(s)  specifies the order of the polynomial in  yp  defining 
c              the adiabatic BOB function for atom A/B
c         if < 0   do not read in any adiabatic BOB function parameters
c  pAD(s)/qAD(s)  are the powers defining the expansion variables
c  UA/UB(a,s)   are the adiabatic BOB function expansion coefficients
c  IFXU(A/B)(a,s) indicates whether each expansion coefficient is to be
c           > 0 :  held fixed at read-in value, or
c        .le. 0 :  varied in the fit
c  UAinf/UBinf  is the limiting asymptotic value of uA(r)/uB(r), as per
c               Theochem paper [internally stored as  UA(NUA+1), etc.]
c  IFXUAinf/IFXUBinf  specifies whether (>0) or not (.le.0)  UAinf/UBinf 
c               is to be held fixed at the read-in value
c=======================================================================
      READ(5,*) NUA(ISTATE),NUB(ISTATE),pAD(ISTATE),qAD(ISTATE)
      IF(NUA(ISTATE).GE.0) THEN
c... NOTE that  parameters  NUA(ISTATE)+1  are  UAinf & IFXUAinf ...
          NUA(ISTATE)= NUA(ISTATE)+ 1
          DO  I= 0, NUA(ISTATE)
              READ(5,*) UA(I,ISTATE), IFXUA(I,ISTATE)
              ENDDO
c=======================================================================
          IF(BOBCN(ISTATE).GT.0) THEN
              UA(NUA(ISTATE),ISTATE)= 0.d0
              IFXUA(NUA(ISTATE),ISTATE)= 1
              ENDIF
          ENDIF
c=======================================================================
      IF(NUB(ISTATE).GE.0) THEN
c... NOTE that  parameters  NUB(ISTATE)+1  are  UBinf & IFXUBinf ...
          NUB(ISTATE)= NUB(ISTATE)+ 1
          DO  I= 0, NUB(ISTATE)
              READ(5,*) UB(I,ISTATE), IFXUB(I,ISTATE)
              ENDDO
c=======================================================================
          IF(BOBCN(ISTATE).GT.0) THEN
              UB(NUB(ISTATE),ISTATE)= 0.d0
              IFXUB(NUB(ISTATE),ISTATE)= 1
              ENDIF
          ENDIF 
c***********************************************************************
c** Read parameters defining the BOB non-adiabatic centrifugal functions
c** If NISTP= 1 , read only one set of non-adiabatic parameters
c
c  NTA/NTB(s)  specifies the order of the polynomial in  yp  defining
c              the non-adiabatic centrifugal BOB functions for atom A/B
c         if < 0   do not read in any non-adiabatic BOB parameters
c  pNA(s)/qNA(s)  are the powers defining the expansion variables
c  TA/TB(a,s)   are the non-adiabatic centrifugal BOB expansion coeffts
c  IFXTA/IFXTB(a,s) indicates whether each expansion coefficient is to be
c           > 0 :  held fixed at read-in value, or
c        .le. 0 :  varied in the fit
c  TAinf/TBinf  is the limiting asymptotic value of qA(r)/qB(r), as per 
c               Theochem paper [internally stored as  TA(NTA+1), etc.]
c  IFXTAinf/IFXTBinf  specifies whether (>0) or not (.le.0)  TAinf/TBinf 
c               is to be held fixed at the read-in value
c=======================================================================
      READ(5,*) NTA(ISTATE),NTB(ISTATE),pNA(ISTATE),qNA(ISTATE)
      IF(NTA(ISTATE).GE.0) THEN
c... NOTE that  parameters  NTA(ISTATE)+1  are  TAinf & IFXTAinf ...
          NTA(ISTATE)= NTA(ISTATE)+ 1
          DO  I= 0, NTA(ISTATE)
              READ(5,*) TA(I,ISTATE), IFXTA(I,ISTATE)
              ENDDO
c=======================================================================
          IF(BOBCN(ISTATE).GT.0) THEN
              TA(NTA(ISTATE),ISTATE)= 0.d0
              IFXTA(NTA(ISTATE),ISTATE)= 1
              ENDIF
          ENDIF 
c=======================================================================
      IF(NTB(ISTATE).GE.0) THEN
c... NOTE that  parameters  NTB(ISTATE)+1  are  TBinf & IFXTBinf ...
          NTB(ISTATE)= NTB(ISTATE)+ 1
          DO  I= 0, NTB(ISTATE)
              READ(5,*) TB(I,ISTATE), IFXTB(I,ISTATE)
              ENDDO
c=======================================================================
          IF(BOBCN(ISTATE).GT.0) THEN
              TB(NTB(ISTATE),ISTATE)= 0.d0
              IFXTB(NTB(ISTATE),ISTATE)= 1
              ENDIF
          ENDIF
c
      NwCFT(ISTATE)= -1
      IF((IOMEG(ISTATE).GT.0).OR.(IOMEG(ISTATE).EQ.-1)) THEN
c-----------------------------------------------------------------------
c** If electronic angular momentum not zero for this state, read Lambda 
c       doubling or doublet Sigma radialfunction parameters.
c*  NwCFT(s)  is order of the polynomial representing the radial fx.
c*  qw(s)  defined nature of radial expansion variable:  
c         y_q= [R^{Pqw} - Re^{Pqw}]/[R^{Pqw} + Re^{Pqw}]  
c*  efREF(s) defines reference level for the Lambda doubling splitting
c            = -1  treats f level as the reference
c            =  0  treats the mid-point between e and f as reference
c            =  1  treats e level as the reference
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          READ(5,*) NwCFT(ISTATE), Pqw(ISTATE), efREF(ISTATE)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          IF(IABS(efREF(ISTATE)).GT.1) THEN
              WRITE(6,646) efREF(ISTATE)
              STOP
              ENDIF
          IF(NwCFT(ISTATE).GE.0) THEN
c... NOTE that  parameters  NwCFT(ISTATE)+1  are  wCFTinf & IFXwCFTinf 
ccc           NwCFT(ISTATE)= NwCFT(ISTATE)+ 1   !!! NOT ANY LONGER !!!
              DO  I= 0, NwCFT(ISTATE)
                  READ(5,*) wCFT(I,ISTATE), IFXwCFT(I,ISTATE)
                  ENDDO
c=======================================================================
              IF(IOMEG(ISTATE).GT.0) THEN
                  IF(efREF(ISTATE).EQ.-1) WRITE(6,640) SLABL(ISTATE)
                  IF(efREF(ISTATE).EQ.0) WRITE(6,642) SLABL(ISTATE)
                  IF(efREF(ISTATE).EQ.1) WRITE(6,644) SLABL(ISTATE)
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
  600 FORMAT(/'For state ',A2,' use a fixed potential defined by LEVEL s
     1ubroutine PREPOT'/4x,'Integrate from   RMIN=',f5.2,
     1  '   to   RMAX=',f6.2,'   with mesh   RH=',f8.5)
  604 FORMAT(/' For state ',A2,' represent level energies by independent
     1 term values')
  606 FORMAT(/' For state ',A2,' represent level energies by independent
     1 band constant for each'/40x,'vibrational level of each isotologue
     2'/3x,'No. band constants'/6x,'v','   isotop:  #',I1,
     3   11('   #',I1:) )
  607 FORMAT(4x,12('-'),12('----',A1:) )
  608 Format(3x,I4,9x,12I5)
  610 FORMAT(' *** Input ERROR *** band constant specification  v=',I3,
     1  ' .NE.', I3)
  626 FORMAT(/' For state   ',A2/4x,'integrate from   RMIN=',f5.2,
     1  '   to   RMAX=',f6.2,'   with mesh   RH=',f8.5)
  628 FORMAT(' ** Since  p=',i2,' .LE. [MMLR(NCMM)-MMLR(1)]=',i2,
     1  '  or   [MMLR(NCMM)-MMLR(1)].le.0,   STOP !!')
  640 FORMAT(/' ', A2,' state energies referenced to f-parity levels')
  642 FORMAT(/' ', A2,' state energies referenced to the mid-point betwe
     1en e and f-parity levels')
  644 FORMAT(/' ', A2,' state energies referenced to e-parity levels')
  646 FORMAT(/' *** INPUT ERROR ***  |efREF=',i3,'| > 1')
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

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
      INTEGER NPASS, NALL,ISTATE,IPV,I,ISOT,J,MMN,m,MM1, LSR,
     1                    JSTATE,IPVRe(NSTATEMX)
      REAL*8 DECM(NSTATEMX),BTEMP,UAT,UBT,BINF,RE3,RE6,RE8,T0,T1,ULRe,
     1 C3VAL,C6adj,C9adj,RET,RETSig,RETPi,RETp,RETm,DM(MMAX),DMP(MMAX),
     2 DMPP(MMAX),bTT(-2:2),cDS(-4:0),bDS(-4:0), PVSR
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
              WRITE(6,602) SLABL(ISTATE),nPB(ISTATE),nPB(ISTATE),
     1                                     NSpow(ISTATE),NLpow(ISTATE)
              IF(RREF(ISTATE).LE.0.d0) WRITE(6,553) 
              IF(RREF(ISTATE).GT.0.d0) WRITE(6,555) RREF(ISTATE),
     1                                                    RREF(ISTATE)
              ENDIF
          IF(PSEL(ISTATE).EQ.2) THEN
c** Header printout for MLR potential
              BINF= betaINF(ISTATE)
              WRITE(6,604) SLABL(ISTATE),nPB(ISTATE),nQB(ISTATE)
              IF(NSpow(ISTATE).GE.0) THEN
                  WRITE(6,605) nPB(ISTATE),nPB(ISTATE),nQB(ISTATE),
     1                                     NSpow(ISTATE),NLpow(ISTATE)
                ELSE
                  WRITE(6,680) NLpow(ISTATE)
                  BETA(NLpow(ISTATE),ISTATE)= BINF
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
              WRITE(6,682) BINF,MMLR(1,ISTATE),CmVAL(1,ISTATE),
     1                                                  MMLR(1,ISTATE)
              MM1= 2
              IF((NCMM(ISTATE).GE.4).AND.(MMLR(2,ISTATE).EQ.0)) THEN
                  WRITE(6,606) MMLR(2,ISTATE),CmVAL(2,ISTATE),0
                  MM1= 3
                  ENDIF
              IF(NCMM(ISTATE).GT.1) THEN
                  DO  m= MM1,NCMM(ISTATE)
                      IF(MMLR(m,ISTATE).LE.9) 
     1      WRITE(6,608) MMLR(m,ISTATE),CmVAL(m,ISTATE),MMLR(m,ISTATE)
                      IF(MMLR(m,ISTATE).GT.9) 
     1      WRITE(6,609) MMLR(m,ISTATE),CmVAL(m,ISTATE),MMLR(m,ISTATE)
                      ENDDO
                  ENDIF
              ENDIF
c
          IF(PSEL(ISTATE).EQ.3) THEN
c** Header printout for DELR potential form ...
              WRITE(6,612) SLABL(ISTATE),nPB(ISTATE),
     1                                     NSpow(ISTATE),NLpow(ISTATE)
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
              WRITE(6,623) SLABL(ISTATE), BETA(NLpow(ISTATE)+1,ISTATE),
     1        BETA(NLpow(ISTATE)+2,ISTATE),BETA(NLpow(ISTATE)+3,ISTATE),
     2                      (nPB(ISTATE)),BETA(NLpow(ISTATE)+1, ISTATE)
              MM1= 1
              DO  m= MM1,NCMM(ISTATE)
                  IF(MMLR(m,ISTATE).LE.9) 
     1      WRITE(6,608) MMLR(m,ISTATE),CmVAL(m,ISTATE),MMLR(m,ISTATE)
                  IF(MMLR(m,ISTATE).GT.9) 
     1      WRITE(6,609) MMLR(m,ISTATE),CmVAL(m,ISTATE),MMLR(m,ISTATE)
                  ENDDO
              ENDIF
c
          IF(PSEL(ISTATE).EQ.5) THEN
c** Header printout for Tang-Toennies type potential ...
              WRITE(6,626)
  626 FORMAT(/' State ',A2,' represented by a Tang-Toennies type potenti
     1al' )
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
                      WRITE(6,554) NAME(1)
                      ENDIF
                  IF(NUB(ISTATE).GE.0) THEN
                     WRITE(6,558) NAME(2),pAD(ISTATE),pAD(ISTATE),
     1                                       qAD(ISTATE),NUB(ISTATE)-1
                      WRITE(6,554) NAME(2)
                      ENDIF
                ELSE
                  WRITE(6,564) (pAD(ISTATE),i=1,5)
                ENDIF
              ENDIF
          IF((NTA(ISTATE).GE.0).OR.(NTB(ISTATE).GE.0)) THEN
c** Print description of centrifugal BOB functional forms ...
              IF(BOBCN(ISTATE).LE.0) THEN
                  WRITE(6,559) 
                  IF(NTA(ISTATE).GE.0) THEN
                      WRITE(6,560) NAME(1),pNA(ISTATE),pNA(ISTATE),
     1                                       qNA(ISTATE),NTA(ISTATE)-1
                      WRITE(6,554) NAME(1)
                      ENDIF
                  IF(NTB(ISTATE).GE.0) THEN
                      WRITE(6,560) NAME(2),pNA(ISTATE),pNA(ISTATE),
     1                                       qNA(ISTATE),NTB(ISTATE)-1
                      WRITE(6,554) NAME(2)
                      ENDIF
                ELSE
                  WRITE(6,562) (pNA(ISTATE),i=1,5)
                ENDIF
              ENDIF
          IF(NwCFT(ISTATE).GE.0) THEN
c** Print description of Lambda/2-Sigma doubling functional forms ...
              IF(IOMEG(ISTATE).GT.0) THEN
                  WRITE(6,618) 'Lambda',Pqw(ISTATE),NwCFT(ISTATE),
     1              Pqw(ISTATE),Pqw(ISTATE),Pqw(ISTATE),Pqw(ISTATE),
     2              Pqw(ISTATE)
                  IF(efREF(ISTATE).EQ.-1) WRITE(6,692) SLABL(ISTATE)
                  IF(efREF(ISTATE).EQ.0) WRITE(6,694) SLABL(ISTATE)
                  IF(efREF(ISTATE).EQ.1) WRITE(6,696) SLABL(ISTATE)
                  ENDIF
              IF(IOMEG(ISTATE).EQ.-1) THEN
                  WRITE(6,618) ' Gamma',Pqw(ISTATE),NwCFT(ISTATE),
     1              Pqw(ISTATE),Pqw(ISTATE),Pqw(ISTATE),Pqw(ISTATE),
     2              Pqw(ISTATE)
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
          IF(NPASS.GT.1) WRITE(20,671) NSpow(ISTATE), NLpow(ISTATE),
     1                          nPB(ISTATE), nQB(ISTATE), RREF(ISTATE)
          NALL= MAX(NSpow(ISTATE),NLpow(ISTATE))
          J=0
          IF(NSpow(ISTATE).LT.0) J=1
          DO  I=J, NALL
              IPV= IPV + 1
              IF(IFXBETA(I,ISTATE).LE.0) THEN
                  IF(DABS(BETA(I,ISTATE)).GT.PU(IPV)) THEN
                      IF(NSpow(ISTATE).GE.0) WRITE(6,640) 'be','ta',I,
     1                                   BETA(I,ISTATE),PU(IPV),PS(IPV)
                      IF(NSpow(ISTATE).LT.0) WRITE(6,640) 'be','ta',I,
     1                  BETA(I,ISTATE),PU(IPV),PS(IPV),ypBETA(I,ISTATE)
                    ELSE
                      IF(NSpow(ISTATE).GE.0) WRITE(6,641) 'be','ta',I,
     1                                   BETA(I,ISTATE),PU(IPV),PS(IPV)
                      IF(NSpow(ISTATE).LT.0) WRITE(6,641) 'be','ta',I,
     1                  BETA(I,ISTATE),PU(IPV),PS(IPV),ypBETA(I,ISTATE)
                    ENDIF
                ELSE
                  IF(NSpow(ISTATE).GE.0) WRITE(6,638) 'be','ta',I,
     1                                                   BETA(I,ISTATE)
                  IF(NSpow(ISTATE).LT.0) WRITE(6,638) 'be','ta',I,
     1                                  BETA(I,ISTATE),ypBETA(I,ISTATE)
                ENDIF
              IF(NPASS.GT.1) THEN
                  IF(NSpow(ISTATE).GE.0) THEN
                      WRITE(20,670) BETA(I,ISTATE),IFXBETA(I,ISTATE)
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
c                           IF(NSpow(ISTATE).GE.0) 
c     1        WRITE(6,640) 'pa','rm',I,BETA(I,ISTATE),PU(IPV),PS(IPV)
c                           IF(NSpow(ISTATE).LT.0) 
c     1        WRITE(6,640) 'pa','rm',I,BETA(I,ISTATE),PU(IPV),PS(IPV),
c     2                                               ypBETA(I,ISTATE)
c                          ELSE
c                           IF(NSpow(ISTATE).GE.0) 
c     1        WRITE(6,641) 'pa','rm',I,BETA(I,ISTATE),PU(IPV),PS(IPV)
c                           IF(NSpow(ISTATE).LT.0) 
c     1        WRITE(6,641) 'pa','rm',I,BETA(I,ISTATE),PU(IPV),PS(IPV),
c     2                                               ypBETA(I,ISTATE)
c                          ENDIF
c                      ELSE
c                        IF(NSpow(ISTATE).GE.0) 
c     1                        WRITE(6,638) 'pa','rm',I,BETA(I,ISTATE)
c                        IF(NSpow(ISTATE).LT.0) 
c     1       WRITE(6,638) 'pa','rm',I,BETA(I,ISTATE),ypBETA(I,ISTATE)
c                      ENDIF
c                   ENDDO
c               ENDIF

          IF(PSEL(ISTATE).EQ.1) BINF= BTEMP
c** Write out  phi_\infty  constant for the MLR/MLJ form
          IF(PSEL(ISTATE).EQ.2) THEN
               IF(NSpow(ISTATE).GE.0) WRITE(6,648) BINF
               IF(NSpow(ISTATE).LT.0) WRITE(6,649) BINF
               BTEMP= CmVAL(1,ISTATE)*2.d0*(2.d0*BINF - BTEMP)
     1                                        *RE(ISTATE)**nPB(ISTATE)
               WRITE(6,652) MMLR(1,ISTATE)+nPB(ISTATE),BTEMP
               ENDIF
c-----------------------------------------------------------------------
c** Writing out the adiabatic BOB radial function for atom A.
c-----------------------------------------------------------------------
          IF(NPASS.GT.1) WRITE(20,672) NUA(ISTATE)-1, NUB(ISTATE)-1,
     1                                        pAD(ISTATE), qAD(ISTATE)
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
                  IF(NPASS.GT.1) WRITE(20,670) UA(I,ISTATE),
     1                                                 IFXUA(I,ISTATE)
                  END DO
              IPV= IPV + 1
              IF(NPASS.GT.1) WRITE(20,674) UA(NUA(ISTATE),ISTATE),
     1                             IFXUA(NUA(ISTATE),ISTATE),'uAinf  '
              IF(IFXUA(NUA(ISTATE),ISTATE).LE.0) THEN
                  WRITE(6,644) ' u',NAME(1),UA(NUA(ISTATE),ISTATE),
     1                                                 PU(IPV),PS(IPV)
                ELSE
                  WRITE(6,646) ' u',NAME(1),UA(NUA(ISTATE),ISTATE)
                ENDIF
              ENDIF    
  668 FORMAT(1Pd20.12,d20.12,0PI3)
  670 FORMAT(1Pd20.12,0PI3,9x,'% ',A2,' IFX',A2)
  671 FORMAT(/2I3,I4,I3,1PD11.2,7x,'% NSpow NLpow nPB nQB RREF')
  672 FORMAT(/2I3,I4,I3,1PD11.2,D12.2,5x,
     1                            '% nuA nuB pAD qAD ')
  673 FORMAT(/2I3,I4,I3,1PD11.2,D11.2,5x,
     1                            '% nTA nTB pNA qNA ')
  674 FORMAT(1Pd20.12,0PI3,11x,'% ',a7)
  675 FORMAT(/3I3,1PD11.2,11x,'%  NwCFT Pqw efREF')
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
                  IF(NPASS.GT.1) WRITE(20,670) UB(I,ISTATE),
     1                                                 IFXUB(I,ISTATE)
                  END DO
              IPV= IPV + 1
              IF(NPASS.GT.1) WRITE(20,674) UB(NUB(ISTATE),ISTATE),
     1                             IFXUB(NUB(ISTATE),ISTATE),'uBinf  '
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
          IF(NPASS.GT.1) WRITE(20,673) NTA(ISTATE)-1, NTB(ISTATE)-1,
     1         pNA(ISTATE),qNA(ISTATE)
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
                  IF(NPASS.GT.1) WRITE(20,670) TA(I,ISTATE),
     1                                                 IFXTA(I,ISTATE)
                  END DO
              IPV= IPV + 1
              IF(NPASS.GT.1) WRITE(20,674) TA(NTA(ISTATE),ISTATE),
     1                             IFXTA(NTA(ISTATE),ISTATE),'tAinf  '
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
                  IF(NPASS.GT.1) WRITE(20,670) TB(I,ISTATE),
     1                                                 IFXTB(I,ISTATE)
                  END DO
              IPV= IPV + 1
              IF(NPASS.GT.1) WRITE(20,674) TB(NTB(ISTATE),ISTATE),
     1                             IFXTB(NTB(ISTATE),ISTATE),'tBinf  '
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
                  IF(NPASS.GT.1) WRITE(20,670) wCFT(I,ISTATE),
     1                                               IFXwCFT(I,ISTATE)
                  END DO
              ENDIF
   90     CONTINUE
      WRITE(6,600)
      RETURN
c-----------------------------------------------------------------------
  553 FORMAT(11x,'with radial expansion variable:   yp = (R^p - Re^p)/(R
     1^p + Re^p)')
  554 FORMAT(7x,'with ',A2,'-atom radial expansion variable:   yp = (R^p
     1 - Re^p)/(R^p + Re^p)')
  555 FORMAT(11x,'with radial variable:   yp = (R^p -',F9.6,'^p)/(R^p +'
     1  ,F9.6,'^p)')
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
  602 FORMAT(/' State ',A2,' represented by an EMO(p=q=',i2,') potential
     1 defined in terms of'/1x,4('=='),'  exponent coefficient:   beta(R
     2)= Sum{beta_i*y',i1,'^i}'/6x,'using exponent  power series order',
     4  i3,' for  R < Re','  and',i3,' for  R > Re')
  604 FORMAT(/' State ',A2,' represented by an MLR(p=',i2,', q=',i2,
     1  ') potential defined in terms of')
  605 FORMAT(1x,4('=='), '  exponent coefficient:  beta(R)= betaINF*y',
     1 i1,' +(1-y',I1,')*Sum{beta_i*y',i1,'^i}'/6x,'using exponent power
     2 series orders',i3,' for  R < Re  and',i3,' for  R > Re')
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
  682 FORMAT(4x,'and   betaINF=',F14.10,'   defined by:  C',I1,'=',
     1 1PD14.7,'[cm-1 Ang^',i1,']')
  608 FORMAT(48x,'C',I1,'=',1PD14.7,'[cm-1 Ang^',i1,']')
  609 FORMAT(47x,'C',I2,'=',1PD14.7,'[cm-1 Ang^{',i2,'}]')
  606 FORMAT(3x,'Use Aubert-Frecon model for uLR(r) with',
     1  4x,'C',I1,'=',1PD13.6,'[cm-1 Ang^',i1,']'/6x,'with retarded C3 a
     2nd  Delta(V)_{gu}^{(6,7)}')
  610 FORMAT(/' For state ',A2,"  use Surkus' Generalized Potential Ener
     1gy Function GPEF with"/1x,6('=='),'  expansion vble:   y_',i1, 
     2 '(r) = (r^',i1,' - re^',i1,')/(',F5.2,'*r^',i1,' +',F5.2,'*re^',
     3 i1,'p)')
  612 FORMAT(/' State ',A2,' represented by a DELR(p=',i2,') potential w
     1ith an exponent'/1x,4('=='),'  coefficient power series of order',
     2  i3,' for  R < Re','  and',i3,' for  R >= Re')
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
     1  '   C',0Pi2,'=',1Pd14.6)
  629 FORMAT((51x,'C',0Pi2,'=',1Pd14.6:)) 
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
  649 FORMAT(3X,'beta_INF',1PD21.12,7X,'--',12X,'--',5x,'at   yp=  1.000
     1000000')
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
c               ----- Version of  18 November 2012 -----
c           (after removal or RREFna, RREFad & RREFw variables!)
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
c* NSpow(s)  is the order of the  beta(r)  exponent expansion for R.le.Re
c* NLpow(s)  is the order of the  beta(r)  exponent expansion for R > Re
c!!!  The code (currently) assumes  NLpow .ge. NSpow !!!!
c*  MMLR(j,s)  are long-range inverse-powers for an MLR or DELR potential
c*  nPB(s)  the basic value of power p for the beta(r)  exponent function
c*  nQB(s)  the power p for the power series expansion variable in beta(r)
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
     9 B4,B5,xBETA(NbetaMX),rKL(NbetaMX,NbetaMX),C1LIM,BETA0,BETAN,
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
c** For the Expanded Morse Oscillator - with  NSpow/NLpow
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
              npow= NSpow(ISTATE)
              IF(RVAL.GT.RE(ISTATE)) npow= NLpow(ISTATE)
              DO  J= 1,npow
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
              DO  J= 0,npow
                  IPV= IPV+1
                  DVtot(IPV,I)= YPP
                  YPP= YPP*YP
                  ENDDO
              ENDDO
          ENDIF
ccccc Print for testing
cc      write(10,610) (RD(i,ISTATE),vpot(i,istate),BETAFX(i,istate),
cc   1                                        i= 1, NDATPT(ISTATE),50)
ccccc
c******** End preparation of Expanded Morse Potential Function *********

      IF(PSEL(ISTATE).EQ.2) THEN
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** For the  {Morse/Long-Range}_p  potential.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** First - define values & derivatives of uLR at Re
          IF((NCMM(ISTATE).GE.4).AND.(MMLR(2,ISTATE).LE.0)) THEN
c** First - define  C6adj & C9adj for the 2x2 and 3x3 Aubert-Frecon cases
              C3VAL= CmVAL(1,ISTATE)
              C3bar= C3VAL/DE(ISTATE)
              C6bar= CmVAL(3,ISTATE)/DE(ISTATE)
              C6adj= CmVAL(3,ISTATE) + 0.25D0*C3VAL*C3bar
              C6Pi= CmVAL(4,ISTATE)
              C9adj= 0.5d0*C3bar*C6adj
              dC9dC6= 0.5d0*C3bar
              dC9dC3= 0.5d0*C6bar + 1.5d0*dC9dC6**2
              dC9dDe= - dC9dC6*(C6bar + C3bar*dC9dC6)
              dC6dDe= - dC9dC6**2
              RE3= 1.d0/RE(ISTATE)**3
              RE6= RE3*RE3
              RE9= RE6*RE3
c** Include QED retardation function in C3 terms for Li2(A) !!
c  NOTE ... the numerical factor is  2\pi/\lambda  for this case
c  Ignore effect of retardation on the partial derivatives & C6adj & C9adj.
              RET= 9.36423830d-04*RE(ISTATE)
              RETSig= DCOS(RET) + (RET)*DSIN(RET)
              RETPi= RETSig - RET**2 *DCOS(RET)
              RETp= RETSig + 0.5d0*RETPi
              RETm= RETSig - 0.5d0*RETPi
              IF(MMLR(2,ISTATE).EQ.0) THEN
c ... for Aubert-Frecon 2x2 treatment of {C3,C6,C8} for Alkali A-state
                  T1= RE3*(C3VAL*RETm + RE3*(C6adj - C6Pi))/3.d0
                  IF(NCMM(ISTATE).GT.4) THEN
c ... extension for Aubert-Frecon Li2(A) {3,0,6,6,8,8} case ...
                      C8val= CmVAL(5,ISTATE)
                      C8Pi= CmVAL(6,ISTATE)
                      RE8= RE6/RE(ISTATE)**2
                      T1= T1+ (C8VAL - C8Pi)*RE8/3.d0
                      ENDIF
                  T0= DSQRT((T1- CmVAL(2,ISTATE))**2 + 8.d0*T1**2)
                  ULRe= 0.5d0*(-CmVAL(2,ISTATE) + RE3*(C3VAL*RETp
     1                   + RE3*(C6adj + C6Pi)) + T0) + RE9*C9adj 
                  IF(NCMM(ISTATE).GT.4) ULRe= ULRe
     1                                      + 0.5d0*(C8VAL + C8Pi)*RE8
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
                  IF(NCMM(ISTATE).GT.4) THEN
                      dLULRedCm(5)= RE8*T0P23/ULRe
                      dLULRedCm(6)= RE8*(1.d0 - T0P23)/ULRe
                      dLULRedRe= dLULRedRe - RE9*4.d0*(C8VAL + C8Pi
     1                            + 2.d0*T0P*(C8VAL - C8Pi)/3.d0)/ULRe
                      ENDIF
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
          IF(NSpow(ISTATE).LT.0) THEN
c*** For Pashov-natural-spline exponent coefficient ...
              DO  I= 1,NLpow(ISTATE)
                  xBETA(I)= ypBETA(I,ISTATE)
                  ENDDO
              BETA(NLpow(ISTATE),ISTATE)= BINF
              CALL Lkoef(NLpow(ISTATE),xBETA,rKL,NbetaMX)
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
              IF(NSpow(ISTATE).GE.0) THEN
c** For 'conventional' power series or A-F exponent function ...
                  VAL= BETA(0,ISTATE)
                  npow= NSpow(ISTATE)
                  IF(RVAL.GT.RE(ISTATE)) npow= NLpow(ISTATE)
                  DO  J= 1,npow
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
                  npow= NLpow(ISTATE)
                  DO  J= 1,npow
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
                  IF(MMLR(2,ISTATE).EQ.0) THEN
c ... extension for Li2(A) 2x2 {3,0,6,6,8,8} case ...
                      T1= (C3VAL*RETm + (C6adj - C6Pi)*RD3)*RD3/3.d0
                      IF(NCMM(ISTATE).GT.4) THEN
                          RD8= RD6/RVAL**2
                          T1= T1+ (C8VAL - C8Pi)*RD8/3.d0
                          ENDIF
                      T0= DSQRT((T1- CmVAL(2,ISTATE))**2 + 8.d0*T1**2)
                      ULR= 0.5d0*(-CmVAL(2,ISTATE) + RD3*(C3VAL*RETp
     1                   + RD3*(C6adj + C6Pi))) + 0.5d0*T0 + C9adj*RD9
                      IF(NCMM(ISTATE).GT.4) ULR= ULR
     1                                      + 0.5d0*(C8VAL + C8Pi)*RD8
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
                  IF(NCMM(ISTATE).GT.4) THEN
c... and then derivatives w.r.t. Auber-Frecon C8 values ...
                      IPV= IPV+1
                      DVtot(IPV,I)= YPP*((1.d0 - YP*YPE)*dLULRedCm(5)
     1                                               - dULRdCm(5)/ULR)
                      IPV= IPV+1
                      DVtot(IPV,I)= YPP*((1.d0 - YP*YPE)*dLULRedCm(6)
     1                                               - dULRdCm(6)/ULR)
                      ENDIF

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
c ... note that derivative w.r.t. atomic level splitting FIXED
                  DVtot(IPV,I)= 0.d0
                  IPV= IPV+1
                  DVtot(IPV,I)= YPP*((1.d0 - YP*YPE)*dLULRedCm(3)
     1                                               - dULRdCm(3)/ULR)
                  IPV= IPV+1
                  DVtot(IPV,I)= YPP*((1.d0 - YP*YPE)*dLULRedCm(4)
     1                                               - dULRdCm(4)/ULR)
                  IF(NCMM(ISTATE).GT.4) THEN
c ... and then derivatives w.r.t. Aubert-Frecon C8 values ...
                      IPV= IPV+1
                      DVtot(IPV,I)= YPP*((1.d0 - YP*YPE)*dLULRedCm(5)
     1                                               - dULRdCm(5)/ULR)
                      IPV= IPV+1
                      DVtot(IPV,I)= YPP*((1.d0 - YP*YPE)*dLULRedCm(6)
     1                                               - dULRdCm(6)/ULR)
                      ENDIF

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
              IF(NSpow(ISTATE).GE.0) THEN
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
                  DO  J= 0,npow
                      IPV= IPV+1
                      DVtot(IPV,I)= YPP
                      YPP= YPP*YQ
                      ENDDO
                ELSE
c... for Pashov-spline exponent cases...
                  YPP= YPP*YPE
                  DO  J= 1,NLpow(ISTATE)
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
     1                                        i= 1, NDATPT(ISTATE),50)
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
              npow= NSpow(ISTATE)
              IF(RVAL.GT.RE(ISTATE)) npow= NLpow(ISTATE)
              BETAFX(I,ISTATE)= BETA(0,ISTATE)
              DO  J= 1,npow
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
              IF(RVAL .LT. RE(ISTATE)) THEN
                  DO  J= 1,NSpow(ISTATE)
                      Btemp2= Btemp2 + BETA(J,ISTATE)*DBLE(J)
     1                                                   *YP**(J-1)
                      ENDDO
                  BTEMP= Btemp2
                ELSE
                  DO  J= 1,NLpow(ISTATE)
                      Btemp2= Btemp2 + BETA(J,ISTATE)*DBLE(J) 
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
c ... finally, derivatives of the potential w.r.t. the \beta_i
              IPV= IPV+ 1
              DVtot(IPV,I)= -PBtemp2 + PBTEMP + dVdBtemp*(Xtemp2-XTEMP)
              DO  J= 1,NSpow(ISTATE)
                  IPV= IPV+ 1
                  DVtot(IPV,I)= (-PBtemp2 + PBTEMP) * YP**J
                  ENDDO
              DO  J= NSpow(ISTATE)+1, NLpow(ISTATE)
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
              DO  J= 1, NLpow(ISTATE)
                  DBDB(J,I,ISTATE)= DBDB(J-1,I,ISTATE)*YP
                  ENDDO
              ENDDO
          ENDIF
******7* End Double-Exponential Long-Range Potential Function ********

      IF(PSEL(ISTATE).EQ.4) THEN
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** For the Tiemann Polynomial Potential Energy Function.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          BT= BETA(NLPow(ISTATE)+1, ISTATE)
          Rinn= BETA(NLPow(ISTATE)+2, ISTATE)
          Rout= BETA(NLPow(ISTATE)+3, ISTATE)
c** Determine analytic function attaching smoothly to inner wall of 
c  polynomial expansion at  R= Rinn < Rm
          YP= (Rinn - RE(ISTATE))/(Rinn+ BT*RE(ISTATE))
          YPP= 1.d0
          A1= BETA(0,ISTATE)
          A2= 0.d0
          DO  J= 1, NLpow(ISTATE)
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
          DO  J= 1, NLpow(ISTATE)
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
                  DO  J= 1,NLpow(ISTATE)
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
              DO  J= 1, NLpow(ISTATE)
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
c           DO  m= 1, NCMM(ISTATE)
c               MMLR1D(m)= MMLR(m,ISTATE)
c               ENDDO
c** Now, calculate A, D_e and their derivatives w.r.t. b and r_e
c           REQ= RE(ISTATE)
c           BTT= BETA(0,ISTATE)
c           KDER= 2
c           CALL dampF(REW,rhoAB(ISTATE),NCMM(ISTATE),
c    1              MMLR1D,IDF(ISTATE),IDSTT(ISTATE),KDER,Dm,Dmp,Dmpp)
c           ATT= 0.d0
c           VATT= 0.d0
c           DO  m= 1,NCMM(ISTATE)
c               TM= CMval(m,ISTATE)/REQ**MMLR1D(m)
c               VATT= VATT + Dm(m)*TM
c               ATT= ATT+ CMval(m,ISTATE)*(MMLR1D(m)*Dm(MMLR1D(m))/REQ
c    1                 - Dmp(MMLR1D(m)))/REQ**MMLR1D(m)
c               ENDDO
c           ATT= ATT*DEXP(BTT*RE(ISTATE))/BTT
c           DE(ISTATE)= ATT*DEXP(-BTT*REQ) - VATT
            ENDIF

      IF(PSEL(ISTATE).EQ.6) THEN
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** For the Surkus-type Generalized Potential Energy Function.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** First, we calculate the implied Dissociation Energy
        IF(AGPEF(ISTATE).NE.0.0d0) THEN
            YPP= 1.d0/AGPEF(ISTATE)**2
            VAL= YPP
            DO  I= 1, NLpow(ISTATE)
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
            DO  J= 1, NLpow(ISTATE)
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
            DO  J= 2, NLpow(ISTATE)
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
          BFCT= 16.85762920d0/(ZMASS(3,IISTP)*RDIST**2)
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

c***********************************************************************
      SUBROUTINE dampF(r,rhoAB,NCMM,MMLR,IDF,IDSTT,KDER,DM,DMP,DMPP)
c** Subroutine to generate values 'Dm' and its first `Dmp' and second
c   'Dmpp' derivatives w.r.t. R of the chosen version of the incomplete
c    gamma function damping function, for  m= 1 to MMAX.
c---------------------- RJL Version of 14 July 2011 --------------------
c-----------------------------------------------------------------------
c                 Upon Input
c* r - the radial distance in Angsroms (!) 
c* RHOab  'universal' scaling coefficient used for systems other than H_2
c       RHOab= 2*(RHOa*RHOb)/(RHOa+RHOb) where RHOa = (I_p^A/I_p^H)^0.66
c              where I_p^A is the ionization potential of atom A
c              and I_p^H is the ionization potential of atomic hydrogen
c* NCMM  the number of inverse-power terms to be considered
c* MMLR  are the powers of the NCMM inverse-power terms
c* IDF requires damping to be defined s.th.  Dm(r)/r^m --> r^{IDF/2}
c* IDSTT specifies damping function type:  > 0  use Douketis et al. form 
c                               if  IDSTT .LE. 0  use Tang-Toennies form
c* KDER:  if KDER.GT.0  the first derivative is also calculated 
c*        if KDER.GT.1  the second derivative is also calculated 
c-----------------------------------------------------------------------
c                 Upon Output
c  DM(m) - The value of the damping function for the long range term 
c          C_MMLR(m)/r^MMLR(m)    {m= 1, NCMM}
c  DMP(m) - The first derivative of the damping function  DM(m)
c  DMPP(m) - The second derivative of the damping function  DM(m)
c-----------------------------------------------------------------------
      INTEGER NCMM,NCMMax,MMLR(NCMM),IDF,IDSTT,KDER,IDFF,FIRST,
     1  Lsr,m,MM,MMAX
      REAL*8 r,rhoAB,bTT(-2:2),cDS(-4:0),bDS(-4:0),aTT,br,XP,YP,
     1  TK, DM(NCMM),DMP(NCMM),DMPP(NCMM),SM(-3:25),
     2  bpm(20,-2:0), cpm(20,-2:0),ZK
c------------------------------------------------------------------------
c  The following values for the numerical factors used in both TT and DS
c  were  normalized to the Hydrogen data presented
c  by Kreek and Meath in J.Chem.Phys. 50, 2289 (1969).
c  The ratio has been chosen such that  b= FACTOR*(I_p^X / I_p^H)^{2/3}
c  for the homoatomic diatomic species X_2, where I_p^A is the ionization
c------------------------------------------------------------------------
       DATA bTT/2.10d0,2.44d0,2.78d0,3.13d0,3.47d0/
       DATA bDS/2.50d0,2.90d0,3.3d0,3.69d0,3.95d0/
       DATA cDS/0.468d0,0.446d0,0.423d0,0.405d0,0.390d0/
       DATA FIRST/ 1/
       SAVE FIRST, bpm, cpm
c------------------------------------------------------------------------
      IF(RHOab.LE.0) THEN
          WRITE(6,602) RHOab
          STOP
          ENDIF
      IF(IDSTT.LE.0) THEN
c===========================================
c** For Tang-Toennies type damping functions
c===========================================
          IF((IDF.LT.-4).OR.(IDF.GT.4)) THEN
                WRITE(6,600) IDSTT,IDF
                STOP
                ENDIF
          Lsr= IDF/2
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
                  IF(KDER.GT.0) THEN
                      DMP(m)= aTT*XP*(SM(MM) - SM(MM-1))
                      IF(KDER.GT.1) DMPP(m)= -aTT*aTT*XP*(SM(MM) 
     1                                     - 2.d0*SM(MM-1) + SM(MM-2))
                      ENDIF
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
                  IF(KDER.GT.0) THEN
                      DMP(m)= aTT*XP*SM(m-1)
                      IF(KDER.GT.1)DMPP(m)= aTT*aTT*XP*(SM(m-2)-SM(m-1))
                      ENDIF
                  ENDDO
            ENDIF
          ENDIF
c
      IF(IDSTT.GT.0) THEN
c=======================================================================
c** For Douketis-Scoles-Marchetti-Zen-Thakkar type damping function ...
c=======================================================================
          IF((IDF.LT.-4).OR.(IDF.GT.0)) THEN
              WRITE(6,600) IDSTT,IDF
              STOP
              ENDIF
          IF(FIRST.EQ.1) THEN
              DO m= 1, 20
                  DO  IDFF= -2,0
                      bpm(m,IDFF)= bDS(IDFF)/DFLOAT(m)
                      cpm(m,IDFF)= cDS(IDFF)/DSQRT(DFLOAT(m))
                      ENDDO
                  ENDDO
              FIRST= 0 
              ENDIF
          br= rhoAB*r
          DO m= 1, NCMM
              MM= MMLR(m)
              XP= DEXP(-(bpm(MM,IDF) + cpm(MM,IDF)*br)*br)
              YP= 1.d0 - XP
              ZK= MM-1.d0
              DM(m)= YP**(MM-1)
c... Actually ...  DM(m)= YP**(MM + IDF/2)  :  set it up this way to 
c   avoid taking exponential of a logarithm for fractional powers (slow)
              IF(IDF.EQ.-4) THEN
                  ZK= ZK- 1.d0
                  DM(m)= DM(m)/YP
                  ENDIF
              IF(IDF.EQ.-3) THEN
                  ZK= ZK- 0.5d0
                  DM(m)= DM(m)/DSQRT(YP)
                  ENDIF
              IF(IDF.EQ.-1) THEN
                  ZK= ZK+ 0.5d0
                  DM(m)= DM(m)*DSQRT(YP)
                  ENDIF
              IF(IDF.EQ.0) THEN
                  ZK= MM
                  DM(m)= DM(m)*YP
                  ENDIF
              IF(KDER.GT.0) THEN
                  TK= bpm(MM,IDF) + 2.d0*cpm(MM,IDF)*br
                  DMP(m) = ZK*XP*rhoAB*TK*DM(m)/YP
                  IF(KDER.GT.1) THEN
c ... if desired ... calculate second derivative [for DELR case] {check this!}
                      DMPP(m)= (ZK-1.d0)*XP*TK*DMP(m)/YP
     1               - DMP(m)*TK + DMP(m)*2.d0*cpm(MM,IDF)*rhoAB**2/TK
                      ENDIF
                  ENDIF
              ENDDO   
          ENDIF  
      RETURN
  600 FORMAT(/,' *** ERROR ***  For  IDSTT=',i3,'   IDF=',i3,'  no dampi
     1ng function is defined')
  602 FORMAT( /,' ***ERROR ***  rhoAB=', F7.4,'  yields an invalid Dampi
     1ng Function definition')
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE AF3X3potRet(RDIST,DELTAE,C3val,C6val,
     1               C8val,De,ULR,DEIGM1,DEIGM3,DEIGM5,DEIGR,DEIGDe)
c***********************************************************************
c** Subroutine to generate the lowest eigenvalue of the 3x3 long-range
c  Li2 interaction matrix of Eq.(25) of J.Mol.Spectr. 268, 199 (2011), 
c  and its derivatives w.r.t. the potential function parameters.
c==> Input: distance  r= RDIST, DELTAE, C3, C6eff, C8 & De.
c==> Output:  ULR= -lambda(min) and iuts partial derivatives DEIGxx
c** Version from Nik in June 2011 for Real root case; cosmetically 
c    modified by RJL in Aug 2011
c-----------------------------------------------------------------------
      REAL*8  A(3,3),DM1(3,3),DM3(3,3),DM5(3,3),DR(3,3),
     1              DDe(3,3),Q(3,3)
      REAL*8  DEIGM1(1,1),DEIGM3(1,1),DEIGM5(1,1),DEIGR(1,1),
     1        DEIGDe(1,1), EIGVEC(3,1), RESID(3,1), W(3) 
      REAL*8  RDIST,RDIST2,RDIST3,DELTAE,C3val,C6val,C8val,De,ULR,
     1   RET,RETSig,RETPi,Modulus,M1,M3,M5,Z
      INTEGER          I,J,L,K
      M1= C3val
      M3= C6val
      M5= C8val
      RET= 9.36423830d-4*RDIST
      RETSig= DCOS(RET) + (RET)*DSIN(RET)
      RETPi= RETSig - RET**2 *DCOS(RET)
      RDIST2= RDIST**2
      RDIST3= RDIST*RDIST2
c      WRITE(25,*) 'Variables = "r", "U(r)","U(r)-U(r)^2/(4De)" ' 
c      WRITE(25,*) 'zone T = "U(r)"'
c  Initialize interaction matrix to 0.d0
      DO  I= 1,3
          A(I,I)=0.0D0
          ENDDO
ccccc Prepare interation matrix  A 
      A(1,1)= -(M1*RETSig+ M3/(RDIST3)+M5/(RDIST3*RDIST2))/(3.d0*RDIST3)
      A(1,2)= -(DSQRT(2.D0))*A(1,1)
      A(2,1)= A(1,2)
      A(1,3)= M1*RETPi/(DSQRT(6.D0)*RDIST3)
      A(3,1)= A(1,3)
      A(2,2)= 2*A(1,1) + DELTAE
      A(2,3)= A(1,3)/DSQRT(2.d0)
      A(3,2)= A(2,3)
      A(3,3)= DELTAE
cccccc Prepare radial derivative of interaction matrix (? is it needed ?)
      DR(1,1)= (3.d0*M1*RETSig + 6.d0*M3/RDIST3 
     1                  + 8.D0*M5/(RDIST3*RDIST2))/(3.d0*RDIST3*RDIST)
      DR(1,2)= -DSQRT(2.d0)*DR(1,1)
      DR(2,1)= DR(1,2)
      DR(2,2)= 2.d0*DR(1,1)
      DR(1,3)= -3.d0*A(1,3)/RDIST
      DR(3,1)= DR(1,3)
      DR(2,3)= -3.d0*A(2,3)/RDIST
      DR(3,2)= DR(2,3)
      DR(3,3)= 0.d0 
cccccc Partial derivative of interaction matric  H  w.r.t.  C3
      DM1(1,1)= -(RETSig + M1/(2.d0*De*RDIST3))/(3.d0*RDIST3)
      DM1(1,2)= -DSQRT(2.d0)*DM1(1,1)
      DM1(2,1)= DM1(1,2)
      DM1(2,2)= 2.d0*DM1(1,1)
      DM1(1,3)= RETPi/(DSQRT(6.d0)*RDIST3)
      DM1(3,1)= DM1(1,3)
      DM1(2,3)= DM1(1,3)/DSQRT(2.d0)
      DM1(3,2)= DM1(2,3)
      DM1(3,3)= 0.d0
cccccc Partial derivative of interaction matric  H  w.r.t.  C6
      DM3(1,1)= -1.d0/(3.d0*RDIST3**2)
      DM3(1,2)= -SQRT(2.d0)*DM3(1,1)
      DM3(1,3)= 0.D0
      DM3(2,1)= DM3(1,2)
      DM3(2,2)= 2.d0*DM3(1,1)
      DM3(2,3)= 0.D0
      DM3(3,1)= DM3(1,3)
      DM3(3,2)= DM3(2,3)
      DM3(3,3)= 0.D0
cccccc Partial derivative of interaction matric  H  w.r.t.  C8
      DM5(1,1)= DM3(1,1)/(RDIST2)
      DM5(1,2)= DM3(1,2)/(RDIST2)
      DM5(1,3)= 0.D0
      DM5(2,1)= DM3(1,2)
      DM5(2,2)= DM3(2,2)/(RDIST2)
      DM5(2,3)= 0.D0
      DM5(3,1)= DM5(1,3)
      DM5(3,2)= DM5(2,3)
      DM5(3,3)= 0.D0
cccccc Partial derivative of interaction matric  H  w.r.t.  De
      DDe(1,1)= M1**2/(12.D0*(RDIST3*De)**2)
      DDe(1,2)= -SQRT(2.D0)*DDe(1,1)
      DDe(1,3)= 0.D0
      DDe(2,1)= DDe(1,2)
      DDe(2,2)= 2.D0*DDe(1,1)
      DDe(2,3)= 0.d0
      DDe(3,1)= DDe(1,3)
      DDe(3,2)= DDe(2,3)
      DDe(3,3)= 0.D0
cccccc Call subroutine to prepare and invert interaction matrix  H
      CALL DSYEVJ3(A,Q,W)
      L=1
ccc Nor - identify the lowest eigenvalue of  H  and label it  L
      DO J=2,3
          IF (W(J) .LT. W(L)) THEN
              L=J
              ENDIF
          ENDDO  
      ULR= -W(L)
      DO I=1,3      
          EIGVEC(I,1) = Q(I,L)
          ENDDO  
      DEIGM1= -MATMUL(TRANSPOSE(EIGVEC),MATMUL(DM1,EIGVEC))
      DEIGM3= -MATMUL(TRANSPOSE(EIGVEC),MATMUL(DM3,EIGVEC))
      DEIGM5= -MATMUL(TRANSPOSE(EIGVEC),MATMUL(DM5,EIGVEC))           
      DEIGR = -MATMUL(TRANSPOSE(EIGVEC),MATMUL(DR,EIGVEC))
      DEIGDe= -MATMUL(TRANSPOSE(EIGVEC),MATMUL(DDe,EIGVEC))
c     WRITE(25,600) RDIST ,ULR 
c 600 FORMAT(2D16.7)
c     WRITE(26,601) RDIST , DEIGM1, DEIGR ,DEIGDe
c 601 FORMAT(4D16.7)  
      RETURN
      CONTAINS
c=======================================================================
      SUBROUTINE DSYEVJ3(A, Q, W)
c ----------------------------------------------------------------------------
c** Subroutine to setup and diagonalize the matrix  A  and return 
c   eigenvalues W and eigenvector matric  Q
      INTEGER N, I, X, Y, R
      PARAMETER (N=3)
      REAL*8 A(3,3), Q(3,3), W(3)
      REAL*8 SD, SO, S, C, T, G, H, Z, THETA, THRESH
c     Initialize Q to the identitity matrix
c     --- This loop can be omitted if only the eigenvalues are desired ---
      DO  X = 1, N
          Q(X,X) = 1.0D0
          DO  Y = 1, X-1
              Q(X, Y) = 0.0D0
              Q(Y, X) = 0.0D0
              ENDDO
          ENDDO
c Initialize W to diag(A)
      DO  X = 1, N
          W(X) = A(X, X)
          ENDDO
c Calculate SQR(tr(A))  
      SD = 0.0D0
      DO  X = 1, N
          SD = SD + ABS(W(X))
          ENDDO
      SD = SD**2
c Main iteration loop
      DO 40 I = 1, 50
c Test for convergence
          SO = 0.0D0
          DO  X = 1, N
              DO  Y = X+1, N
                  SO = SO + ABS(A(X, Y))
                  ENDDO
              ENDDO
          IF(SO .EQ. 0.0D0)  RETURN
          IF(I .LT. 4) THEN
              THRESH = 0.2D0 * SO / N**2
            ELSE
              THRESH = 0.0D0
            END IF
c Do sweep
          DO 60 X = 1, N
              DO 61 Y = X+1, N
                  G = 100.0D0 * ( ABS(A(X, Y)) )
                  IF ( I .GT. 4 .AND. ABS(W(X)) + G .EQ. ABS(W(X))
     $                          .AND. ABS(W(Y)) + G .EQ. ABS(W(Y))) THEN
                      A(X, Y) = 0.0D0
                    ELSE IF (ABS(A(X, Y)) .GT. THRESH) THEN
c Calculate Jacobi transformation
                      H = W(Y) - W(X)
                      IF ( ABS(H) + G .EQ. ABS(H) ) THEN
                          T = A(X, Y) / H
                        ELSE
                          THETA = 0.5D0 * H / A(X, Y)
                          IF (THETA .LT. 0.0D0) THEN
                              T= -1.0D0/(SQRT(1.0D0 + THETA**2)-THETA)
                            ELSE
                              T= 1.0D0/(SQRT(1.0D0 + THETA**2) + THETA)
                            END IF
                        END IF
                      C = 1.0D0 / SQRT( 1.0D0 + T**2 )
                      S = T * C
                      Z = T * A(X, Y)
c Apply Jacobi transformation
                      A(X, Y) = 0.0D0
                      W(X)    = W(X) - Z
                      W(Y)    = W(Y) + Z
                      DO  R = 1, X-1
                          T       = A(R, X)
                          A(R, X) = C * T - S * A(R, Y)
                          A(R, Y) = S * T + C * A(R, Y)
                          ENDDO
                      DO  R = X+1, Y-1
                          T       = A(X, R)
                          A(X, R) = C * T - S * A(R, Y)
                          A(R, Y) = S * T + C * A(R, Y)
                          ENDDO
                      DO  R = Y+1, N
                          T       = A(X, R)
                          A(X, R) = C * T - S * A(Y, R)
                          A(Y, R) = S * T + C * A(Y, R)
                          ENDDO
c Update eigenvectors
c --- This loop can be omitted if only the eigenvalues are desired ---
                      DO  R = 1, N
                          T       = Q(R, X)
                          Q(R, X) = C * T - S * Q(R, Y)
                          Q(R, Y) = S * T + C * Q(R, Y)
                          ENDDO
                    END IF
   61             CONTINUE
   60         CONTINUE
   40     CONTINUE
      PRINT *, "DSYEVJ3: No convergence."
      END SUBROUTINE DSYEVJ3
      END SUBROUTINE AF3X3potRet
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
              NVVPP(ISOT,ISTATE)= 0
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
      SUBROUTINE DYIDPJ(IDAT,NDATA,NPTOT,JFXP,YOBS,YC,PV,PD,PSS,
     1  RMSR)
c***********************************************************************
c** This program assumes that the upper state IEP is at a higher energy
c   than the lower state IEPP.
c** This subroutine returns the calculated value YC of datum IDAT, and
c  its partial derivatives PD(k) w.r.t. the NPTOT parameters PV.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++  COPYRIGHT 2007-11  by  R.J. Le Roy, Jenning Seto and Yiye Huang +
c   Dept. of Chemistry, Univ. of Waterloo, Waterloo, Ontario, Canada   +
c    This software may not be sold or any other commercial use made    +
c      of it without the express written permission of the authors.    +
c++++++++++++++++++++ (version of 24/02/2012) ++++++++++++++++++++++++++
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
      INCLUDE 'BLKPARAM.h'
      INCLUDE 'BLKBOB.h'
      INCLUDE 'BLKCOUNT.h'
c-----------------------------------------------------------------------
c** Common block for partial derivatives of potential at the one distance RDIST
      REAL*8 dVdPk(HPARMX)
      COMMON /dVdPkBLK/dVdPk
c=======================================================================
      INTEGER IDAT,NPTOT,NBAND,IISTP,ISTATE,I,J,JFXP(NPARMX),NDATA
c
c** Define parameters required locally and from NLLSSRR.
      REAL*8 RDIST,VDIST,BETADIST,EUP,ELW,YOBS,YC,EO,RMSR,width,
     1  PV(NPARMX),PD(NPARMX),PSS(NPARMX),UPPER(NPARMX),LOWER(NPARMX),
     2  DEDPK(HPARMX),BVIR,dBVIRdP(NPARMX)
c
c----------------------------------------------------------------------- c
c** Vibrational energies and rot. constants for generating trial energies
      INTEGER INNR(0:NVIBMX)
      REAL*8 ZK(0:NVIBMX,0:NROTMX,NISTPMX,NSTATEMX)
c
      IF(IDAT.EQ.1) THEN
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++ At beginning of each fit cycle (datum #1), re-map internal NLLSSRR
c  parameter array PV onto external (physical) variable set, and get
c  updated band constant array ZK for estimating trial eigenvalues.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          CALL MAPPAR(NISTP,PV,1)
          REWIND(22)
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
     1                                                      BETADIST,0)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Now generate band constants for each state and isotopomer for 
c  generating trial eigenvalues in fit calculations.  To take account of
c  the different vibrational ranges for different isotopomers, only 
c  generate values for levels to the input VMAX for isotopomer-1.
              IF(PSEL(ISTATE).GE.0) THEN
                  DO  IISTP= 1,NISTP
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Call subroutine INITDD to calculate band constants for initial
c   trial eigenvalue estimates as well as the partial derivatives.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                    CALL INITDD(ISTATE,IISTP,VMAX(ISTATE,IISTP),
     1                                                        ZK,INNR)
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
          CALL DEDP(IDAT,IEPP(NBAND),IISTP,ZMASS(3,IISTP),
     1               JP(IDAT),JPP(IDAT),EFPP(IDAT),ELW,width,ZK,LOWER)
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
        ELSEIF(IEP(NBAND).EQ.-1) THEN
c=======================================================================
c*** For PAS data ...
c=======================================================================
          CALL DEDP(IDAT,IEPP(NBAND),IISTP,ZMASS(3,IISTP),
     1               JP(IDAT),JPP(IDAT),EFPP(IDAT),ELW,width,ZK,LOWER)
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
        ELSEIF(IEP(NBAND).EQ.-2) THEN
c=======================================================================
c*** If datum is width of a tunneling-predissociation quasibound level
c ... for forward calculation, use widths from SCHRQ in DEDP
c=======================================================================
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
        ELSEIF(IEP(NBAND).EQ.-3) THEN
c=======================================================================
c*** If datum is the potential energy function at some specific distance
c=======================================================================
          RDIST= TEMP(IDAT)
ccc       ISTATE= IEPP(NBAND)
          CALL VGEN(IEPP(NBAND),RDIST,VDIST,BETADIST,IDAT)
          YC= VDIST
          DO  J= POTPARI(IEPP(NBAND)), POTPARF(IEPP(NBAND))
              UPPER(J)= dVdPk(J) 
              ENDDO
          RETURN
        ELSEIF(IEP(NBAND).EQ.-4) THEN
c=======================================================================
c*** For Virial coefficient data ...
c=======================================================================
          CALL DVIRDP(IDAT,IEPP(NBAND),ZMASS(3,IISTP),BVIR,PD)
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
              CALL DEDP(IDAT,IEP(NBAND),IISTP,ZMASS(3,IISTP),
     1              VP(NBAND),JP(IDAT),EFP(IDAT),EUP,width,ZK,UPPER)
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
     1           VPP(NBAND),JPP(IDAT),EFPP(IDAT),ELW,width,ZK,LOWER)
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
  600 FORMAT(' *** FAIL to find level(',A2,')   v=',I3,'   J=',I3,
     1  '  so ignore  YOBS(',i5,')=',f12.4)
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE INITDD(ISTATE,IISTP,VIMX,ZK,INNR)
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
      INCLUDE 'BLKBOB.h'
      INCLUDE 'BLKBOBRF.h'
      INCLUDE 'BLKDATA.h'
      INCLUDE 'BLKISOT.h'
      INCLUDE 'BLKPOT.h'
c-----------------------------------------------------------------------
c
      INTEGER ISTATE, IISTP, VIMX, AFLAG, 
     1  VIN, NBEG, NEND, WARN, IWR, LPRWF, I, J, INNODE, KV
c
      REAL*8 GV(0:NVIBMX),RCNST(NROTMX),RM2(NPNTMX),V(NPNTMX),
     1  SWF(NPNTMX),FWHM,PMAX,BFCT,BvWN,RHSQ,C3gu,T3
c
c** Vibrational energies and rot. constants for generating trial energies
      INTEGER INNR(0:NVIBMX)
      REAL*8 ZK(0:NVIBMX,0:NROTMX,NISTPMX,NSTATEMX),SWF2,qDBL,qFCT
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IF((ISTATE.EQ.1).AND.(IISTP.EQ.1)) THEN
          REWIND(7)
          REWIND(21)
          ENDIF
      VIN= VIMX
      AFLAG= 0
      WARN= 0
      IWR= -1
ccc   iwr= 5
      LPRWF= 0
      INNODE= 1
c
c** Calculate potential scaling factors for ALF and CDJOEL
c
      RHSQ= RH(ISTATE)*RH(ISTATE) 
      BFCT= (ZMASS(3,IISTP)/16.85762920D0)*RHSQ
      BvWN= 16.85762920D0/ZMASS(3,IISTP)
      qFCT= RH(ISTATE)*BvWN**(2*IOMEG(ISTATE))
c
c** Now generate the BFCT-scaled and adiabatically corrected potential
c   for the current isotope for use in SCHRQ, plus the 1/R**2 array for
c   CDC calculations in CDJOEL
c
      DO  I=1,NDATPT(ISTATE)
          V(I)= BFCT*(VPOT(I,ISTATE)+ ZMUA(IISTP,ISTATE)*UAR(I,ISTATE)
     1                             + ZMUB(IISTP,ISTATE)*UBR(I,ISTATE))
          RM2(I)= 1.d0/RD(I,ISTATE)**2
c** Special BOB correction for A-state Li2 and analogous cases. !!!!!!!!
          IF(IOMEG(ISTATE).EQ.-2) V(I)= V(I) + 2.d0*RHSQ*RM2(I)
          ENDDO 
      IF((NCMM(ISTATE).GE.3).AND.(MMLR(2,ISTATE).LE.0).AND.(AN(1).EQ.3)
     1        .AND.(AN(2).EQ.3).AND.(MN(1,IISTP).NE.MN(2,IISTP))) THEN
c** Add g/u symmetry breakdown correction for special case {6,7}Li2(A) !!
          C3gu= (2.d0/3.d0)*CmVAL(1,ISTATE)
          DO I= 1,NDATPT(ISTATE)
              T3= C3gu/RD(I,ISTATE)**3
              V(I)= V(I) + BFCT*(T3 - DSQRT(T3**2 + 0.03085959756d0))
              ENDDO
          ENDIF
c!!
      IF((NTA(ISTATE).GE.0).OR.(NTB(ISTATE).GE.0)) THEN
          DO  I= 1, NDATPT(ISTATE)
              RM2(I)= RM2(I)*(1.d0 + ZMTA(IISTP,ISTATE)*TAR(I,ISTATE) 
     1                             + ZMTB(IISTP,ISTATE)*TBR(I,ISTATE))
              ENDDO
          ENDIF
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Call subroutine ALF that will locate the needed vibrational levels 
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      CALL ALF(NDATPT(ISTATE),RMIN(ISTATE),RH(ISTATE),MMLR(1,ISTATE),V,
     1 SWF,VLIM(ISTATE),VIN,NVIBMX,AFLAG,ZMASS(3,IISTP),EPS(ISTATE),GV,
     1 INNODE,INNR,IWR)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** If a serious error occured during within ALF, then print out a 
c   warning, record the constants that we have and hope the program
c   doesn't call on the constants for which ALF could not calculate.
c
      IF (AFLAG.LT.0) THEN
          WRITE(6,600) ISTATE,IISTP
          IF(AFLAG.EQ.-1) THEN
              WRITE(6,601) VIN,VIMX, AFLAG
              STOP
              ENDIF
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
c   of this state for this isotopologue
      IF(NwCFT(ISTATE).LT.0) THEN
          WRITE(7,614)  NAME(1),MN(1,IISTP),NAME(2),MN(2,IISTP)
        ELSE
          WRITE(7,618)  NAME(1),MN(1,IISTP),NAME(2),MN(2,IISTP)
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
          CALL CDJOEL(GV(KV),NBEG,NEND,NPNTMX,BvWN,RH(ISTATE),WARN,V,
     1                                                  SWF,RM2,RCNST)
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
  601 FORMAT(4X,'ALF finds highest level   v=',i3,'  is below the desire
     1d (v=',I3,', J=',I3,')')
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
  614 FORMAT(/' For  ',A2,'(',I3,') - ',A2,'(',I3,')'/1x,11('--')/
     1  '   v       ','E',12x,'Bv',11x,'-Dv',13x,'Hv',13x,'Lv',
     2  12x,'Mv',13x,'Nv',13x,'Ov'/1x,58('=='))
  616 FORMAT(I4,f12.4,f14.10,7(1PD15.7))
  618 FORMAT(/' For  ',A2,'(',I3,') - ',A2,'(',I3,')'/1x,11('--')/
     1  '   v       ','E',12x,'Bv',11x,'-Dv',13x,'Hv',13x,'Lv',
     2  13x,'Mv',13x,'Nv',13x,'Ov'13x,'qB(v)'/1x,67('=='))
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
c          The first line identifies the level, gives the position of
c          1-st point and radial mesh, & states No. of  points.
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
      INTEGER IDAT, bandN, efPARITY, JRe
c
      INTEGER I,ICOR,ISTATE,IISTP,INNER,J,JROT,JIN,KVLEV,KV, 
     1  NBEG,NEND,INNODE, IWR, LPRWF
c
      REAL*8 ZMU,EO,BFCT,JFCT,JFCTP,FWHM,UMAX, RM2, ETRY, BMAX, DGDV2,
     1  SWF2,JFCTA,JFCTB, DUARe,DUBRe,DTARe,DTBRe,DLDRe,DEROT,DEROTB,
     2  JFCTDBL,JFCTD, muFCT,C3gu,T3, GV(0:NVIBMX),Vtotal(NPNTMX),
     3  SWF(NPNTMX), V(NPNTMX), DEDPK(HPARMX)
c
      COMMON /VBLIK/Vtotal
c
c** Vibrational energies and rotational constants.
c
      REAL*8 ZK(0:NVIBMX,0:NROTMX,NISTPMX,NSTATEMX)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c** Initializing the arrays and variables
      DATA IWR/0/,LPRWF/0/,INNODE/1/,BMAX/-9.d9/
c
c** To calculate the values for <vj|dV/dP(k)|vj>
c   we must first determine the wave equation for the vj state.
      muFCT= 16.85762920d0/ZMU
      BFCT= RH(ISTATE)*RH(ISTATE)/muFCT
      JFCT= DBLE(JROT*(JROT+1))
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
               JFCTP= JFCT
               IF(NBC(KVLEV,IISTP,ISTATE).GT.1) THEN
                   DO  J= 2, NBC(KVLEV,IISTP,ISTATE)
                       I= I+ 1
                       DEDPK(I)= JFCTP
                       EO= EO+ ZBC(KVLEV,J-1,IISTP,ISTATE)*JFCTP
                       JFCTP= JFCTP*JFCT 
                       ENDDO
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
              IF(IOMEG(ISTATE).EQ.-1) JFCTP= JFCTP - 2**I
              DEROT= ZK(KVLEV,I,IISTP,ISTATE) * JFCTP
c... if centrifugal term bigger than previous one - truncate summation
              IF(DABS(DEROT).GT.DABS(DEROTB)) GOTO 4
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
          IF(efPARITY.GT.0) JFCTDBL=  0.5d0*DBLE(JROT)*muFCT
          IF(efPARITY.EQ.0) JFCTDBL=  0.d0
          IF(efPARITY.LT.0) JFCTDBL= -0.5d0*DBLE(JROT+1)*muFCT
          ENDIF
c
c** Generating potential function including centrifugal, adiabatic BOB,
c   and rotational non-adiabatic BOB, and if appropriate, Lambda 
c   doubling or 2\Sigma doubling radial functions.
cc
c==> Model-A <==
c** First - for Li2(A), add BOB centrifugal shift  2*B(r)
cc    IF((NCMM(ISTATE).GE.3).AND.(MMLR(2,ISTATE).LE.0)
cc   1                      .AND.(IOMEG(ISTATE).NE.-2))
cc   2              JFCT= JFCT + 2.d0*ZMASS(3,1)/ZMASS(3,IISTP) - 2.d0
c==> Model-A <==
cc
c** Include ad BOB function for A-state Li2 & analogous cases
cc
      IF(IOMEG(ISTATE).LE.-2) JFCT= JFCT - IOMEG(ISTATE) 
      JFCT= JFCT* RH(ISTATE)* RH(ISTATE)
      JFCTD= JFCTDBL* RH(ISTATE)* RH(ISTATE)/muFCT
      bandN= IB(IDAT)
      DO I= 1,NDATPT(ISTATE)
          RM2= 1/RD(I,ISTATE)**2
          V(I)= BFCT*(VPOT(I,ISTATE) + ZMUA(IISTP,ISTATE)*UAR(I,ISTATE)
     1                              + ZMUB(IISTP,ISTATE)*UBR(I,ISTATE))
     2          + JFCT* (1.0d0 + ZMTA(IISTP,ISTATE)*TAR(I,ISTATE)
     1                         + ZMTB(IISTP,ISTATE)*TBR(I,ISTATE))*RM2
          IF((IOMEG(ISTATE).GT.0).OR.(IOMEG(ISTATE).EQ.-1)) THEN
c** Add radial potential for Lambda- or 2-Sigma splitting
c ... note that power of 1/r**2  included in  wRAD  array ...
              V(I)= V(I) + JFCTD*wRAD(I,ISTATE)
              ENDIF
          Vtotal(I)= V(I)/BFCT
          ENDDO
c!!
      IF((NCMM(ISTATE).GE.3).AND.(MMLR(2,ISTATE).LE.0)) THEN
          IF((AN(1).EQ.3).AND.(AN(2).EQ.3)
     1                         .AND.(MN(1,IISTP).NE.MN(2,IISTP))) THEN
c** Add g/u symmetry breakdown correction for special case {6,7}Li2(A) !!
              C3gu= (2.d0/3.d0)*CmVAL(1,ISTATE)
              DO I= 1,NDATPT(ISTATE)
                  T3= C3gu/RD(I,ISTATE)**3
                  Vtotal(I)= Vtotal(I)+ T3- DSQRT(T3**2+ 3.085959756d-2)
                  V(I)= Vtotal(I)*BFCT
                  ENDDO
              ENDIF
          ENDIF
      ICOR= 0
      INNER= 0
      EO= ETRY
   10 KV= KVLEV
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      CALL SCHRQ(KV,JROT,EO,FWHM,UMAX,VLIM(ISTATE),V,SWF,BFCT,
     1             EPS(ISTATE),RMIN(ISTATE),RH(ISTATE),NDATPT(ISTATE),
     2                               NBEG,NEND,INNODE,INNER,IWR,LPRWF)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IF(KV.NE.KVLEV) THEN
c** If SCHRQ found the wrong level ....
          write(6,666) KVLEV,JROT,ETRY,KV,EO
  666 FORMAT(' Search for v=',I3,'   J=',I3,'  starting from  E=',
     1  f9.2,' finds  E(v=',I3,')=', f9.2)
          ICOR= ICOR+1
          IF((ICOR.LE.10).AND.(KV.GE.0)) THEN
c... SCECOR uses semiclassical methods to estimate correct energy
              CALL SCECOR(KV,KVLEV,JROT,INNER,ICOR,IWR,EO,RH(ISTATE),
     1   BFCT,NDATPT(ISTATE),MMLR(1,ISTATE),V,BMAX,VLIM(ISTATE),DGDV2)
c***********************************************************************
              KV= KVLEV
              GOTO 10
              ENDIF 
c** If the calculated wavefunction is still for the wrong vibrational
c   level, then write out a warning and skip the calculation of the
c   partial derivatives (hence setting them to zero).
          WRITE(6,610) KVLEV,JROT,KV
c.. eigenvalue of -9.9d9 indicates that eigenvalue search failed completely
          EO= -9.9d9
          IWR= 0
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
      JRe= POTPARI(ISTATE)+ 1
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
      RETURN
c-----------------------------------------------------------------------
  610 FORMAT(' *** SCECOR fails. Seeking   v=',i3,', J=',i3,
     1   ';  Found  v=',I3)
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
c ... obtain potential function value at new turning point estimate
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
c
      REAL*8 dVdPk(HPARMX)
      COMMON /dVdPkBLK/dVdPk
c
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

c***********************************************************************
      SUBROUTINE DVIRDP(IDAT,ISTATE,ZMU,BVIR,dBVIRdP)
c=======================================================================
c  This subroutine calculates the second virial coefficient BVIR and its
c  partial derivatives with respect to the various parameters. It is 
c  used when virial data has been input into the program.  It performs a
c  classical calcultion with quantum corrections
c=======================================================================
      INCLUDE 'arrsizes.h'
      INCLUDE 'BLKPOT.h'
      INCLUDE 'BLKCOUNT.h'
      INCLUDE 'BLKDVDP.h'
      INCLUDE 'BLKISOT.h'
      INCLUDE 'BLKDATA.h'
c* Define types for local variables
       INTEGER check,j,counter,nParams,ISTATE,ISTART,
     1 IDAT,i,m,k,kk,jj,n
       REAL*8  XG(8),WG(8),x(4),w(4),BTemp,BTempInv,
     1 CONST(3),jump,a,b,YVAL(8),RVAL(8),
     2 VDIST(8),EXP_TERM,int_fact,Vsq,VPsq,Rsq,RINV,XTEMP,Bclass,Bq1,
     3 Bq2,INTEGRALS(NPARMX+3),error,BVIR,BETADISTa,dVdR(8),d2VdR(8),
     4 XTEMP1,dBcdP(NPARMX),dBq1dP(NPARMX),dBVIRdP(NPARMX),ZMU,dBVIR
c
       REAL*8, PARAMETER :: k_boltz=1.3806504D-23, NA=6.02214179D23,
     1  h = 6.62606896D-34, c = 2.99792458D10, k_cm = k_boltz/(h*c),
     2  h_cm = 1.d0/c, pi = 3.14159265358979323846, amu=1.660538782D-27,
     3  Cu= 16.8576292D0
c* Common block data for quadrature weights and points
      data x/0.960289856497536d0, 0.796666477413627d0,
     1       0.525532409916329d0, 0.183434642495650d0/,
     2     w/0.101228536290376d0, 0.222381034453374d0,
     3       0.313706645877887d0, 0.362683783378362d0/
c***********************************************************************
      ERROR= 0.001d0
c.. runs a loop to set all quadrature weights and points
      DO j= 1,4
          XG(j)= -x(j)
          WG(j)= w(j)
          XG(9-j)= x(j)
          WG(9-j)= w(j)
          ENDDO
c.. initializes the array dBVIRdP to zero.
      DO J= 1,HPARMX
          dBVIRdP(J)= 0.d0
          ENDDO
c.. takes the current temperature from the virial data and changes it to E cm
      BTemp = k_cm*Temp(IDAT)
c.. inverts BTemp
      BTempInv = 1.d0/BTemp
c.. sets the constants for the intergration of each correction term
      Const(1)= -2.d0*pi*0.602214179d0
      Const(2)= Cu*pi*BTempInv**3/(6.d0*ZMU)
      Const(3)= -(pi/6.d0)*(Cu*BTempInv**2/ZMU)**2

c.. sets the number of paramters used for the integration and initializes
c.  the integrals to zero
      nParams = 3 + NCMM(ISTATE) + NLpow(ISTATE)
      check = 1

c.. this first loop is here to repeat the calculations until the values
c** Outermost loop: repeatedly bisect overall [-1,+1] interval NBISMX
c  times until convergence to within absolute error  ERROR  is achieved
      n= 12
      DO  kk= 1,n
          DO J= 1,nParams+3
              INTEGRALS(J)= 0.d0
              ENDDO
          IF(check.EQ.-1) EXIT
          counter= 2.d0**kk

c.. futher subdivides the interval for gaussian integration formula
          jump= 2.d0/DBLE(counter)
          b= -1.d0
c. the first loop is over each subinterval within (-1,1)
          DO m= 1,counter
              a= b
              b= a+ jump
              IF(m.EQ.counter) b= 1.d0
              DO i= 1,8
c.. sets the y value for the gaussian formula
                  YVAL(i)= 0.5d0*((b - a)*XG(i) + (b + a))
c... mapping  r<->y=(r/re - 0.9999)/(r/re + 1.0001) 
c   where 'real' range is  0.01*Re  to infty - with  'p=2'
                  RVAL(i)= Re(ISTATE)*DSQRT((1.0001d0 
     1                           + 0.9999d0*YVAL(i))/(1.d0 - YVAL(i)))
                  VDIST(i)= 0.d0
                  ENDDO
c.. calls VGENp to find the neccesary derivatives
              CALL VGENp(ISTATE,RVAL,VDIST,BETADISTa,IDAT,dVdR,d2VdR)
c.. the next loop is over each Gaussian point within each subinterval
              DO i= 1,8
c.. some of the terms that will be used in later calculations of the
c.  virial coefficients are constructed here
                  EXP_TERM= DEXP(-VDIST(i)*BTempINV)
                  RINV= 1.d0/RVAL(i)
                  int_fact= (Re(ISTATE)/(1.d0 - YVAL(i)))**2
     1                                          *RINV*(b - a)*0.5d0
                  Vsq= VDIST(i)**2
                  VPsq= dVdR(i)**2
                  Rsq= RVAL(i)**2
                  XTEMP= d2VdR(i)**2/10.d0 + VPsq*RINV**2/5.d0
     1      + dVdR(i)**3*BTempInv*RINV/9.d0 - (VPsq*BTempInv)**2/72.d0
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c*     Finally the relevant partial derivatives of the virial 
c*     coefficients are calculated
                  DO J = 1,nParams
                      jj = ISTART + J
c* the derivative of the classical expression
                      dBcdP(J)= -BTempInv*EXP_TERM
     1                                   *DVtot(jj,i)*Rsq*int_fact
c* and of the first quantum correction 
                      XTEMP1 = -VPsq*BTempInv*DVtot(jj,i)
     1                                     + 2*dVdR(i)*dVpdP(jj,i)
                      dBq1dP(J)= EXP_TERM*XTEMP1*Rsq*int_fact
c..  As the final step these terms all added together in a weighted sum
                      INTEGRALS(J) = INTEGRALS(J) + Const(1)
     1                  *dBcdP(J)*WG(I) + Const(2)*dBq1dP(J)*WG(I)
                      ENDDO
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c. now the integrands are evaluated at each particular Gaussian point
c. and summed together with the proper weighting
                      Bclass= (EXP_TERM - 1.d0)*Rsq*int_fact
                      INTEGRALS(nParams+1)= INTEGRALS(NParams+1) 
     1                                                 + Bclass*WG(i)
                      Bq1= EXP_TERM*VPsq*Rsq*int_fact
                      INTEGRALS(nParams+2)= INTEGRALS(nParams+2) 
     1                                                    + Bq1*WG(i)
                      Bq2= EXP_TERM*XTEMP*Rsq*int_fact
                      INTEGRALS(nParams+3)= INTEGRALS(nParams+3)
     1                                                    + Bq2*WG(i)
                  ENDDO
              ENDDO
              
          DO j= 1,nParams
                  dBVIRdP(j)= INTEGRALS(j)
                  ENDDO
          dBVIR=BVIR
          BVIR= Const(1)*INTEGRALS(nParams+1) + 
     2    Const(2)*INTEGRALS(nParams+2) +
     3    Const(3)*Integrals(nParams+3)

c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c. Check to see if BVIR has converged
          IF(counter.GE.2) THEN
              error=0.000001
              IF(DABS(BVIR-dBVIR).LE.(abs(error*BVIR))) THEN
                  check= -1 
                  ENDIF
              ENDIF 
          ENDDO
      END 
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

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
c   RDIST:  at the 8 input RDIST(i) diustances, calculate potl & derivs
c     * return potential function at those points as VDIST, and the
c       first and second radial derivartives as  dVdR  &  d2VdR
c     * return potential function exponent coefficientw as BETADIS
c???  * skip partial derivative calculation if  IDAT.le.0
c     * If RDIST.le.0  calculate partial derivatives at distances
c             given by array RD(i,ISTATE) & return them in array DVtot
c** On entry via common blocks:
c  NSpow(s)  is the order of the  beta(r)  exponent expansion for R.le.Re
c  NLpow(s)  is the order of the  beta(r)  exponent expansion for R > Re
c!!!  The code (currently) assumes  NLpow .ge. NSpow !!!!
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
          npow= NSpow(ISTATE)
          IF(RVAL(I).GT.RE(ISTATE)) npow= NLpow(ISTATE)
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
          BFCT= 16.85762920d0/(ZMASS(3,IISTP)*RDIST(I)**2)
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

c**********************************************************************
      SUBROUTINE DIFFSTATS(NSTATES,ROBUST,MKPRED)
c** Subroutine to summarise dimensionless standard errors on a band-by-
c  band basis, and (if desired) print [obs.-calc.] values to channel-8.
c-----------------------------------------------------------------------
c                 Version of  28 November 2012
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
              IF( (DIV.GE.4.d0).AND.(DIV.LT.10.d0) ) marker='** '
              IF( (DIV.GE.8.d0) ) marker='***'
              WRITE(8,628) JP(I),JPP(I),NEF(EFPP(I)),FREQ(I),
     1                      UFREQ(I),DFREQ(I),DFREQ(I)/UFREQ(I),MARKER
              ENDDO
          WRITE(6,629)
          WRITE(8,629)
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
c
      IF(NVVPP(ISOT,ISTATE).GT.0) THEN
c** Book-keeping for potential function values as data .....
          IBB= YPR(ISOT,ISTATE,5,4,1)
          CALL BNDERR(IFIRST(IBB),ILAST(IBB),ROBUST,AVE,RMSR,SSQTOT,
     1                                                    DFREQ,UFREQ)
          WRITE(6,638) NVVPP(ISOT,ISTATE),SLABL(ISTATE),
     1              (NAME(I),MN(I,ISOT),I=1,2),YPR(ISOT,ISTATE,5,3,1),
     2              AVEUFREQ(IBB),MAXUFREQ(IBB),AVE,RMSR
          WRITE(8,638) NVVPP(ISOT,ISTATE),SLABL(ISTATE),
     1              (NAME(I),MN(I,ISOT),I=1,2),YPR(ISOT,ISTATE,5,3,1),
     2              AVEUFREQ(IBB),MAXUFREQ(IBB),AVE,RMSR
          DO  J= IFIRST(IBB),ILAST(IBB)
              WRITE(6,637) TEMP(J),FREQ(J),UFREQ(J),DFREQ(J),
     1                                               DFREQ(J)/UFREQ(J)
              WRITE(8,637) TEMP(J),FREQ(J),UFREQ(J),DFREQ(J),
     1                                               DFREQ(J)/UFREQ(J)
              ENDDO
          ENDIF
c
      IF(NVIRIAL(ISOT,ISTATE).GT.0) THEN
c** Book-keeping for  Virial Coefficient  data
          IBB= YPR(ISOT,ISTATE,8,4,1)
          CALL BNDERR(IFIRST(IBB),ILAST(IBB),ROBUST,AVE,RMSR,SSQTOT,
     1                                                    DFREQ,UFREQ)
          WRITE(6,634) NVIRIAL(ISOT,ISTATE),SLABL(ISTATE),
     1              (NAME(I),MN(I,ISOT),I=1,2),YPR(ISOT,ISTATE,8,3,1),
     2              AVEUFREQ(IBB),MAXUFREQ(IBB),AVE,RMSR
          WRITE(8,634) NVIRIAL(ISOT,ISTATE),SLABL(ISTATE),
     1              (NAME(I),MN(I,ISOT),I=1,2),YPR(ISOT,ISTATE,8,3,1),
     2              AVEUFREQ(IBB),MAXUFREQ(IBB),AVE,RMSR
          DO  J= IFIRST(IBB),ILAST(IBB)
              WRITE(6,636) TEMP(J),FREQ(J),UFREQ(J),DFREQ(J),
     1                                               DFREQ(J)/UFREQ(J)
              WRITE(8,636) TEMP(J),FREQ(J),UFREQ(J),DFREQ(J),
     1                                               DFREQ(J)/UFREQ(J)
c=======================================================
c   7 State X0 Ar( 40)-Xe(132)  Virial Coefficients
c======================================= Avge. =========
c         #data    Av.Unc.   Max.Unc.   Err/Unc   DRMSD
c-------------------------------------------------------
c            7     3.3D-08   3.3D-08   -0.15741   0.622
c==============================================  calc-obs
c    temp.     Bvir(obs)    u(Bvir)   calc-obs   /u(Bvir)
c--------------------------------------------------------
c   1234.00    -141.22       2.00      23.xxx     5.0000
c   1234.00    -141.22       2.00      23.xxx     5.0000
c   1234.00    -141.22       2.00      23.xxx     5.0000
c--------------------------------------------------------
              ENDDO
  634 FORMAT(/1x,55('=')/I5,'  State ',A2,2x,A2,'(',I3,')-',A2,'(',
     1 I3,')  Virial coefficients'/1x,10('===='),' Avge. ',('====')/5x,
     2 '#data    Av.Unc.    Max.Unc.    Err/Unc   DRMSD'/1x,10('-----')
     3 /I9,1PD12.2,D12.2,0PF11.5,F8.3/1x,23('=='),' calc-obs'/
     4  5x,'temp.',4x,'Bvir(obs)',4x,'u(Bvir)   calc-obs   /u(Bvir)'/
     5  1x,53('-'))
  636 FORMAT(F11.3,2F10.2,2F11.3)
  637 FORMAT(F11.6,F12.2,F10.2,2F11.3)
          ENDIF
  638 FORMAT(/1x,55('=')/I5,'  State ',A2,2x,A2,'(',I3,')-',A2,'(',
     1 I3,')  Potential fx. values'/1x,21('=='),' Avge. ',('====')/5x,
     2 '#data    Av.Unc.    Max.Unc.    Err/Unc   DRMSD'/1x,26('--')
     3 /I9,1PD12.2,D12.2,0PF11.5,F8.3/1x,24('=='),' calc-obs'/
     4  7x,'R',7x,'V(r)',8x,'u(V(r))   calc-obs   /u(V(r))'/
     5  1x,53('-'))
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
     1 I3,') potential fx. values treated as independent data'/
     2 1x,20('=='),'  Avge.  ',17('=')/' #data   v"min  v"max  AvgeUnc',
     1 '  Max.Unc.  Err/Unc  DRMSD'/1x,55('-')/I5,2I7,2x,1P2D9.1,0PF9.3,
     4 F8.3/1x,30('==')/'    v  p',8x,'Bv',7x,'u(Bv)',4x,
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
  627 FORMAT(1x,50('='),'  calc-obs'/'   v   j  p       PAS(Eb)        u
     1(Eb)       calc-obs   /u(FREQ)'/1x,30('--'))
  628 FORMAT(2I4,A3,F15.8,2F11.8,F11.4,1X,A3)
  629 FORMAT(1x,30('=='))
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
      SUBROUTINE PREPOTT(LNPT,IAN1,IAN2,IMN1,IMN2,NPP,VLIM,RR,VV)
c** Driver subroutine of package to generate a potential function VV(i) 
c  at the NPP input distances  RR(i)  by reading, interpolating over and
c  extrapolating beyond a set of up to NPTMX read-in points.  
c  Based subroutine PREPOT of program LEVEL.
c====================== Version of  8 June 2007 ========================
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
c        calculating adiabatic and/or non-adiabatic BOB correction fx.
c  NPP (integer) is the number of input distances  RR(i) (in Angstroms)
c        at which potential values  VV(i) (in cm-1) are to be generated
c  VLIM (cm-1) is the absolute energy at the potential asymptote
c  RR  (real array) is set of NPP distances where potential calculated
c----------------------
c**** Subroutine Output:
c----------------------
c  VV (real array) is the set of function values generated (in cm-1)
c  NCN is an integer power defining the asymptotically-dominant 
c      inverse-power long-range potential tail:  CNN/R**NCN 
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+ Calls GENINT (which calls PLYINTRP, SPLINT & SPLINE) 
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Set maximum array dimension for the input function values to be
c  interpolated over & extrapolated beyong
      INTEGER NTPMX
      PARAMETER (NTPMX= 1600) 
      INTEGER I,J,IAN1,IAN2,IMN1,IMN2,INPTS,ILR,IR2,JWR,LNPT,LWR,
     1  NCN,NLIN,NPP,NROW,NTP,NUSE
      REAL*8  RFACT,EFACT,RH,RMIN,VLIM,VSHIFT,VV(NPP),RR(NPP),RM2(NPP),
     1  XI(NTPMX),YI(NTPMX),RWR(20),VWR(20),VWRB(3),D1V(3),D1VB(3),
     2  D2V(3),CNN
c
c** Save variables needed for 'subsequent' LNPT.le.0 calls
      SAVE ILR,IR2,NTP,NUSE
      SAVE CNN,VSHIFT,XI,YI
c
      DATA VWRB/3*0.D0/,D1VB/3*0.D0/
c
      WRITE(6,600) VLIM
c** For a pointwise potential (NTP > 0), now read points & parameters
c  controlling how the interpolation/extrapolation is to be done.
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** NTP : define potential by interpolation over & extrapolation
c        beyond the NTP read-in turning points using subroutine GENINT.
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
c-------------------------------------------------------------------
      READ(5,*) NTP, NUSE, IR2, ILR, NCN, CNN
c-------------------------------------------------------------------
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
c-------------------------------------------------------------------
      READ(5,*) RFACT, EFACT, VSHIFT
      READ(5,*) (XI(I), YI(I), I= 1,NTP)
c-------------------------------------------------------------------
      WRITE(6,612) VSHIFT, RFACT, EFACT
      NROW= (NTP+2)/3
      DO  J= 1,NROW
          IF(EFACT.LE.10.D0) THEN
              WRITE(6,614) (XI(I),YI(I),I= J,NTP,NROW)
            ELSE
              WRITE(6,616) (XI(I),YI(I),I= J,NTP,NROW)
            ENDIF
          ENDDO
          WRITE(6,624)
      DO  I= 1,NTP
          YI(I)= YI(I)*EFACT+ VSHIFT
          XI(I)= XI(I)*RFACT
          ENDDO
      IF(IR2.GT.0) THEN
          DO  I= 1,NTP
              YI(I)= YI(I)*XI(I)**2
              ENDDO
          ENDIF
      IF((DABS(YI(NTP)-YI(NTP-1)).LE.0).AND.
     1                      (XI(NTP).LT.RR(NPP))) WRITE(6,618)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      CALL GENINT(LNPT,NPP,RR,VV,NUSE,IR2,NTP,XI,YI,VLIM,ILR,NCN,CNN)
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     IF((LNPT.GT.0).AND.(LPPOT.NE.0)) THEN
c** If desired, on the first pass (i.e. if LNPT > 0) print the potential
c         RH= RR(2)-RR(1)
c         INPTS= IABS(LPPOT)
c         IF(LPPOT.LT.0) THEN
c** Option to write resulting function compactly to channel-8. 
c             RMIN= RR(1)
c             NLIN= NPP/INPTS+ 1
c             WRITE(8,800) NLIN,VLIM
c             WRITE(8,802) (RR(I),VV(I),I= 1,NPP,INPTS)
c           ELSE
c** Option to print potential & its 1-st three derivatives, the latter
c  calculated by differences, assuming equally spaced RR(I) values.
c             WRITE(6,620)
c             NLIN= NPP/(2*INPTS)+1
c             RH= INPTS*RH
c             DO  I= 1,NLIN
c                 LWR= 1+ INPTS*(I-1)
c                 DO  J= 1,2
c                     JWR= LWR+(J-1)*NLIN*INPTS
c                     IF(JWR.LE.NPP) THEN
c                         RWR(J)= RR(JWR)
c                         VWR(J)= VV(JWR)
c                         D1V(J)= (VWR(J)-VWRB(J))/RH
c                         VWRB(J)= VWR(J)
c                         D2V(J)= (D1V(J)-D1VB(J))/RH
c                         D1VB(J)= D1V(J)
c                       ELSE
c                         RWR(J)= 0.d0
c                         VWR(J)= 0.d0
c                       ENDIF
c                     IF(I.LE.2) THEN
c                         D2V(J)= 0.d0
c                         IF(I.EQ.1) D1V(J)= 0.d0
c                         ENDIF
c                     ENDDO
c                 WRITE(6,622) (RWR(J),VWR(J),D1V(J),D2V(J),J= 1,2)
c                 ENDDO
c           ENDIF
c         ENDIF
      IF(LNPT.GT.0) WRITE(6,624)
      RETURN
  600 FORMAT(' State has energy asymptote:   Y(lim)=',F12.4,'[cm-1]')
  602 FORMAT(/' **** ERROR in dimensioning of arrays required'
     1 ,' by GENINT;   No. input points ',I5,' > NTPMX =',I4)
  604 FORMAT(' Perform',I3,'-point piecewise polynomial interpolation ov
     1er',I5,' input points' )
  606 FORMAT(' Perform cubic spline interpolation over the',I5,
     1  ' input points' )
  608 FORMAT(' Interpolation actually performed over modified input arra
     1y:   Y(I) * r(I)**2')
  610 FORMAT( ' Beyond read-in points extrapolate to limiting asymptotic
     1 behaviour:'/20x,'Y(r)  =  Y(lim) - (',D16.7,')/r**',I2)
  612 FORMAT(' To make input points Y(i) consistent with  Y(lim),  add'
     1 ,'  Y(shift)=',F12.4/' Scale input points:  (distance)*',
     2 1PD16.9,'  &  (energy)*',D16.9/13x,'to get required internal unit
     3s  [Angstroms & cm-1 for potentials]'/
     4  3('      r(i)         Y(i)  ')/3(3X,11('--')))
  614 FORMAT((3(F13.8,F12.4)))
  616 FORMAT((3(F12.6,F13.8)))
  618 FORMAT(/' !!! CAUTION !!! Last two mesh point  YI  values are equa
     1l'/17x,'so extrapolation to large  r  will be unreliable !!!'/)
  620 FORMAT(/'  Function and first 2 derivatives by differences'/
     1  2('     r       Y(r)     d1Y/dr1    d2Y/dr2')/2(2X,19('--')))
  622 FORMAT(2(0PF8.3,F11.3,1PD11.3,D10.2))
c 622 FORMAT(2(0PF7.2,F12.5,1PD11.3,D10.2))
  624 FORMAT(1x,38('--'))
  800 FORMAT(I7,' function values with asymptotic value:',F14.6)
  802 FORMAT((1X,3(F12.8,F14.6)))
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE GENINT(LNPT,NPP,XX,YY,NUSE,IR2,NTP,XI,YI,VLIM,ILR,
     1                                                        NCN,CNN)
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
     1  NN,NPP,NUSE,NUST,NORD,NCN,NCN2,NCN4,NTP,
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
  600 FORMAT('  ** CAUTION ** Exponential inward extrapolation fails'/
     1 16x,'since first 3 points not monotonic, ... so ...')
  602 FORMAT(' *** CAUTION ** inward extrapolation exponent coefficient
     1   C=',D12.4/10x,'could cause overflows, ... so ...')
  604 FORMAT(' *** CAUTION ** after 15 tries inward extrap. exponent coe
     1fft change is',1PD9.1)
  606 FORMAT(' Extrapolate to   X .le.',F7.4,'  with'/'   Y=',F13.3,
     1  SP,1PD15.6,' * exp(',SS,D13.6,'*X)')
  608 FORMAT(' Extrapolate to   X .le.',F8.4,'   with'/'   Y=',F13.3,
     1  SP,1PD16.7,' * [X - (',SS,F8.4,')]')
  610 FORMAT(' Extrapolate to  X .le.',F8.4,'   with   Y=',F12.3,
     1  SP,1PD15.6,')/X**1')
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
              NCN= NINT(BLR)
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
c         ENDIF
c 627 FORMAT('     Verify smoothness of outer join at   X=',F9.5/
c    1  (3X,3(F10.5,G15.7)))
      RETURN
  612 FORMAT('  *** BUT *** since exponent has positive coefficient, swi
     1tch form ...')
  614 FORMAT(' Function for  X .GE.',F8.4,'   generated as'/'   Y=',
     1  F12.4,' - (',1PD13.6,') * exp{-',0PF10.6,' * (r -',F9.6,')**2}')
  616 FORMAT(' Function for  X .GE.',F8.4,'   generated as'/'   Y=',
     1 F12.4,' - (',1PD13.6,') * r**',0PF10.6,'  * exp{-(',F11.6,'*r)}')
  618 FORMAT(' Extrapolate to  X .GE.',F8.4,'  using'/'   Y=',
     1  F12.4,SP,1PD15.6,'/X**(',SS,D13.6,')] ,  yielding   NCN=',I3)
  620 FORMAT(' Extrapolate to  X .GE.',F8.4,'  using'/'   Y=',
     1  F12.4,' - [',1PD13.6,'/X**',I1,SP,D14.6,'/X**',SS,I1,']')
  622 FORMAT(' Extrapolate to  X .GE.',F8.4,'  using'/
     1  '   Y=',F12.4,' - [',1PD13.6,'/X**',I1,SP,D14.6,'/X**',
     2  SS,I1,SP,D14.6,'/X**',SS,I2,']')
  624 FORMAT(' Function for  X .GE.',F7.3,'  generated by',I3,
     1 '-point inverse-power interpolation'/'   with leading term  1/r**
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
c+++++ Calls only subroutines SPLINE and PLYINTRP ++++++++++++++++++++++
c=======================================================================
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
ccc       CALL SPLINE(R1,V1,NTP,3,CSP,MAXSP,IER)
c... using Pashov 'natural spline' with zero 2'nd derivative @ end points
          CALL SPLINE(R1,V1,NTP,0,CSP,MAXSP,IER)
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
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

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
  250 IER= NER
      RETURN
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE MAPPAR(NISTP,PV,FORBAC)
c***********************************************************************
c** This subroutine will convert external logical physical parameters 
c  into the generic NLLSSRR parameter array PV or the reverse.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++ COPYRIGHT 1997-2011  by  J.Y. Seto & R.J. Le Roy  (ver. 14/07/2011)+
c   Dept. of Chemistry, Univ. of Waterloo, Waterloo, Ontario, Canada   +
c    This software may not be sold or any other commercial use made    +
c      of it without the express written permission of the authors.    +
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c    PV(i)   is the NLLSSRR parameter array.
c    FORBAC  is a flag to determine which way the parameters are mapped
c            FORBAC = 0 : Map internal PV to external variables.
c            FORBAC = 1 : Map external varuables to internal PV.
c*   NSTATES is the number of states being considered (in BLKPARAM)
c=======================================================================
      INCLUDE 'arrsizes.h'
      INCLUDE 'BLKPOT.h'
      INCLUDE 'BLKPARAM.h'
      INCLUDE 'BLKBOB.h'
c-----------------------------------------------------------------------
      INTEGER NISTP, m, FORBAC, ISTATE, IPV, I, J, ISOT
      REAL*8 PV(NPARMX)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Map external free parameters (De, Re, etc.) onto internal NLLSSRR
c  parameters PV(j)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IF(FORBAC.EQ.0) THEN
          IPV= 0
          DO  ISTATE= 1,NSTATES
              IF(PSEL(ISTATE).EQ.-1) THEN
c*** Manage parameters for term value mappings ...
                  DO  ISOT= 1, NISTP
                      DO  I= VMIN(ISTATE,ISOT),VMAX(ISTATE,ISOT)
                          IF(NBC(I,ISOT,ISTATE).GT.0) THEN
                              DO  J= 1,NBC(I,ISOT,ISTATE)
                                  IPV= IPV+1
                                  PV(IPV)= ZBC(I,J-1,ISOT,ISTATE)
                                  ENDDO
                              ENDIF
                          ENDDO
                      ENDDO
                  ENDIF
              IF(PSEL(ISTATE).GT.0) THEN
c*** Manage parameters for potential function mapping ...
                  IF(PSEL(ISTATE).LT.6) THEN
                      IPV= IPV+ 1
                      PV(IPV)= DE(ISTATE)
                      ENDIF
                  IPV= IPV+ 1
                  PV(IPV)= RE(ISTATE)
                  IF((PSEL(ISTATE).GE.2).AND.(PSEL(ISTATE).LE.5)) THEN
                      DO  m= 1,NCMM(ISTATE)
                          IPV= IPV+ 1
                          PV(IPV)= CmVAL(m,ISTATE)
                          ENDDO
                      ENDIF
                  J= 0
                  IF(NSpow(ISTATE).LT.0) J=1
                  DO  I= J,MAX(NSpow(ISTATE),NLpow(ISTATE))
                      IPV= IPV+ 1
                      PV(IPV)= BETA(I,ISTATE)
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
              IF(PSEL(ISTATE).EQ.-1) THEN
c*** Manage parameters for term value mappings ...
                  DO  ISOT= 1, NISTP
                      DO  I= VMIN(ISTATE,ISOT),VMAX(ISTATE,ISOT)
                          IF(NBC(I,ISOT,ISTATE).GT.0) THEN
                              DO  J= 1,NBC(I,ISOT,ISTATE)
                                  IPV= IPV+1
                                  ZBC(I,J-1,ISOT,ISTATE)= PV(IPV)
                                  ENDDO
                              ENDIF
                          ENDDO
                      ENDDO
                  ENDIF
              IF(PSEL(ISTATE).GT.0) THEN
c*** Manage parameters for potential function mappings ...
                  IF(PSEL(ISTATE).LT.6) THEN
                      IPV= IPV + 1
                      DE(ISTATE)= PV(IPV)
                      ENDIF
                  IPV= IPV + 1
                  RE(ISTATE) = PV(IPV)
                  IF((PSEL(ISTATE).GE.2).AND.(PSEL(ISTATE).LE.5)) THEN
                      DO  m= 1,NCMM(ISTATE)
                          IPV= IPV+ 1
                          CmVAL(m,ISTATE)= PV(IPV) 
                          ENDDO
                      ENDIF
                  J=0
                  IF(NSpow(ISTATE).LT.0) J=1
                  DO I= J, MAX(NSpow(ISTATE),NLpow(ISTATE))
                      IPV = IPV + 1
                      BETA(I,ISTATE) = PV(IPV)
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
      SUBROUTINE ALF(NDP,RMIN,RH,NCN,V,SWF,VLIM,KVMAX,NVIBMX,AFLAG,ZMU,
     1                                         EPS,GV,INNODE,INNR,IWR)
c***********************************************************************
c   Version 2.0 dated  June 21, 2011
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
c+++++++++++++   COPYRIGHT 2008-09  by  Robert J. Le Roy   +++++++++++++
c   Dept. of Chemistry, Univ. of Waterloo, Waterloo, Ontario, Canada   +
c    This software may not be sold or any other commercial use made    +
c     of it without the express written permission of the authors.     +
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++++++ Please inform me of any bugs, by phone at: (519)888-4051 +++++++
c+++++++++ by e-mail to: leroy@uwaterloo.ca , or by Post at: +++++++++++
c+++ Dept. of Chemistry, Univ. Waterloo, Waterloo, Ontario  N2L 3G1 ++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Uses the Schrodinger solver subroutine SCHRQ.
c
c** On entry:
c    NDP    is the number of datapoints used for the potential.
c    RMIN   is the inner radial distance of the potential (ang).
c    RH     is the meshvalue (ang).
c           NDP, RMIN, and RH define the radial range over which the
c           potential is defined.
c    NCN    is the (integer) inverse power defining the linmiting attractive
c           long-range behaviour of the potential.  For a barrier, set NCN=99
c    V(i)   is the scaled input potential (cm-1).
c           The scaling factor BFCT is (2*mu/hbar^2)*RH^2.
c    VLIM   is the potential asymptote (cm-1).
c    KVMAX  is v for the highest vibrational level we wish to find.
c    NVIBMX defines dimension of the external Gv array:  GV(0:NVIBMX)
c    AFLAG  is rot.quantum J for the (centrifugally distorted) potential
c    ZMU    is the reduced mass of the diatom (amu).
c    EPS    is the energy convergence criterion (cm-1).
c    INNODE specifies whether wave fx. initiation @ RMIN starts with a
c        note (normal case: INNODE > 0) or zero slope (when INNODE.le.0)
c    IWR    specifies the level of printing inside SCHRQ
c           <> 0 : print error & warning descriptions.
c           >= 1 : also print final eigenvalues & node count.
c           >= 2 : also show end-of-range wave function amplitudes.
c           >= 3 : print also intermediate trial eigenvalues, etc.
c
c** On exit:
c    KVMAX   is vib.quantum number for the highest vibrational level
c            found (may be less than the input value of KVMAX).
c    AFLAG   returns calculation outcome to calling program.
c            >=  0 : found all levels to v=KVMAX{input} & AFLAG= J 
c             = -1 : KVMAX larger than number of levels found.
c    GV(v)   contains the vibrational energy levels found for v=0-KVMAX
c    INNR(v) labels each level as belonging to the inner (INNR = 1) or
c            outer (INNR = 0) well.
c
c** Flags: Modify only when debugging.
c    AWO   specifies the level of printing inside ALF
c          <> 0 : print error & warning descriptions.
c          >  0 : also print intermediate ALF messages.
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
c** NF counts levels found in automatic search option
c
      INTEGER NDP,KVMAX,KV,KVB,KVBB,AFLAG,NF,NBEG,NEND,NVIBMX,
     1  NBEGG(0:NVIBMX),NENDD(0:NVIBMX),INNR(0:NVIBMX),ICOR,IWR,
     2  I,LTRY,AWO,INNODE,INNER,
     3  LPRWF,JROT,NCN,NPMIN, NPMAX, IPMIN(10),IPMINN
c
      REAL*8 RMIN,RMAX,RH,V(NDP),SWF(NDP),VLIM,EO,ZMU,EPS,BZ,BFCT,GAMA,
     1  VMIN,VMAX,PMAX, ESAV, ZPEHO, DGDV2, BMAX,
     2  GV(0:NVIBMX),VPMIN(10),RPMIN(10),VPMAX(10),RPMAX(10)
c
      DATA AWO/1/,LPRWF/0/,KVB/-1/,KVBB/-2/
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Check that the array dimensions are adequate.
      IF(KVMAX.GT.NVIBMX) THEN
          WRITE(6,602) KVMAX, NVIBMX
          STOP
          ENDIF
c
c** Initialize the remaining variables and flags.
      NF= 0                                ! NF is label of level being sought
      KV= 0
      INNER= 0
      LTRY= 0
c** Initialize level counters for each well.
      DO  I= 0,KVMAX
          INNR(I)= -1
          ENDDO
c** Store input rotational quantum number.
      JROT= AFLAG
      AFLAG= -1
c
c** Numerical factor  16.85762920 (+/- 0.00000011) based on Compton
c  wavelength of proton & proton mass (u) from 2002 physical constants.
      BZ= ZMU/16.85762920d0
      BFCT= BZ*RH*RH
c
c** RMAX is the outermost radial distance at which potential is defined
      RMAX= RMIN + DBLE(NDP-1)*RH
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Locate the potential minima.
      NPMIN= 0
      VMIN= 1.d99
      DO  I= 3,NDP-1
          IF((V(I).LT.V(I-1)).AND.(V(I).LT.V(I+1))) THEN
              NPMIN= NPMIN + 1
              IPMIN(NPMIN)= I
              RPMIN(NPMIN)= RMIN + DBLE(I-1)*RH
              VPMIN(NPMIN)= V(I)/BFCT
              IF(VPMIN(NPMIN).LT.VMIN) THEN
                  IPMINN= I
                  VMIN= VPMIN(NPMIN)
                  ENDIF
              IF(NPMIN.EQ.10) GOTO 10
              ENDIF
          END DO
   10 IF(NPMIN.EQ.0) THEN
          IF(V(2).LE.V(1)) THEN
c** If NO minimum & potential has negative slope, print a warning and stop
              WRITE(6,604) JROT
              KVMAX= -1
              RETURN
              ENDIF
c...  but if potl. alway has positive slope, mesh point 1 is minimum
          NPMIN= 1
          IPMIN(NPMIN)= 1
          VPMIN(NPMIN)= V(1)/BFCT
          RPMIN(NPMIN)= RMIN
          VMIN= RPMIN(NPMIN)
          WRITE(6,606) VPMIN(1),RMIN
          ENDIF
c
c** Locate any potential maxima past innermost minimum (if they exists).
      NPMAX= 0
      VMAX= -9.d99
      DO  I= IPMIN(1)+1,NDP-1
          IF((V(I).GT.V(I-1)).AND.(V(I).GT.V(I+1))) THEN
              NPMAX= NPMAX + 1
              RPMAX(NPMAX)= RMIN + DBLE(I-1)*RH
              VPMAX(NPMAX)= V(I)/BFCT
              IF(VPMAX(NPMAX).GT.VMAX) VMAX= VPMAX(NPMAX)
              IF(NPMAX.EQ.10) GOTO 20
              ENDIF
          END DO
   20 IF((NPMAX.EQ.0).OR.
     1         ((NPMAX.GT.0).AND.(RPMAX(NPMAX).LT.RPMIN(NPMIN)))) THEN
c** If no maxima found or there is no barrier past outermost minimum,
c   set an energy maximum to be the value at the end of the radial range.
          NPMAX= NPMAX+ 1
          RPMAX(NPMAX)= RMAX
c?? should this limit be set at  VLIM ??
          VPMAX(NPMAX)= V(NDP-1)/BFCT
          IF(VPMAX(NPMAX).GT.VMAX) VMAX= VPMAX(NPMAX)
          ENDIF
c
c** For multiple minima, print out potential extrema count
      IF(NPMIN.GT.1) THEN
          WRITE(6,614) NPMIN, (VPMIN(I),I= 1,NPMIN)
          WRITE(6,616) (RPMIN(I), I= 1,NPMIN)
          WRITE(6,618) NPMAX, (VPMAX(I),I= 1,NPMAX)
          WRITE(6,616) (RPMAX(I), I= 1,NPMAX)
          IF(NPMIN.GT.2) THEN
c** If potential has more than two minima - print warning & stop
              WRITE(6,620)
              STOP
              ENDIF
          ENDIF
c** Set BMAX as barrier height of double-minimum potential
      BMAX= -9.d+09
      IF(NPMIN.GT.1) THEN
          DO  I= 1,NPMAX
              IF((RPMAX(I).GT.RPMIN(1)).AND.(RPMAX(I).LT.RPMIN(2)))
     1            BMAX= VPMAX(I)
              ENDDO
          ENDIF
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c*** Use harmonic approximation to estimate zero point energy.
      ZPEHO= DSQRT((V(IPMINN+20)-V(IPMINN))/400.d0)/BFCT
      EO= VMIN + ZPEHO
c
c=========== Begin Actual Eigenvalue Calculation Loop Here =============
c** Compute eigenvalues ... etc. up to the KVMAX'th vibrational level.
c** When attempts to find the next eigenvalue fails, then perhaps the
c   next level is located in a second (inner) well. If so, then the
c   subroutine will set INNER = 1, and attempt to find that level.
c
      ICOR= 0
  100 KVBB= KVB
      KVB= KV
      KV= NF
  110 ESAV= EO
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Call subroutine SCHRQ to find eigenvalue EO and eigenfunction SWF(I).
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      CALL SCHRQ(KV,JROT,EO,GAMA,PMAX,VLIM,V,SWF,BFCT,EPS,RMIN,RH,NDP,
     1                               NBEG,NEND,INNODE,INNER,IWR,LPRWF)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IF(KV.LT.0) THEN
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** The SCHRQ error condition is KV < 0.  Allow for 3 cases:
c     EO > VMAX : energy from previous trial above potential maximum
c     NF = 0 : Looking for the first vibrational level (v = 0)
c     NF > 0 : Looking for the other vibrational levels (v > 0)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          IF(EO.GT.VMAX) THEN
c** For the case when the previous trial gave energy above the potential
c   maximum, make one last ditch attempt to find the highest bound level
c   (quasi or otherwise) in the potential.
              IF(LTRY.LT.1) THEN
                  LTRY= 1
                  KV= 999
                  EO= VMAX - 0.01d0
                  GOTO 110
c... if that was unsuccessful, then print out a warning and exit.
                ELSE
                  WRITE(6,622) NF, EO, VMAX
                  KV= NF-1
                  GOTO 200
                ENDIF
              ENDIF
          WRITE(6,624) NF,JROT,ESAV
c.. eigenvalue of -9.9d9 signifies that eigenvalue search failed completely
          KVMAX= NF-1
          EO= -9.9d9
          RETURN
          ENDIF
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** If calculated vibrational level is the desired level, NF, then ...
c   call SCECOR to calculate dG/dv and predict next higher level
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IF(KV.EQ.NF) THEN
          NBEGG(KV)= NBEG
          NENDD(KV)= NEND
          GV(NF)= EO
          INNR(NF)= INNER
  120     NF= NF + 1
          IF(NF.LE.KVMAX) THEN
              IF(INNR(NF).GT.0) GOTO 120
c... if the next level was found earlier in overshoot ... 
            ELSE
              IF((AWO.GT.0).AND.(IWR.GT.0)) WRITE(6,626) JROT,KVMAX
              AFLAG= JROT
              RETURN
            ENDIF
          ICOR= 0
          CALL SCECOR(KV,NF,JROT,INNER,ICOR,IWR,EO,RH,BFCT,NDP,NCN,V,
     1                                                BMAX,VLIM,DGDV2)
          IF(ICOR.GT.20) GOTO 200
          IF(EO.GT.VPMAX(NPMAX)) THEN
c... if estimated energy above highest barrier, set value below it
              EO=  VPMAX(NPMAX) - 0.10d0*DGDV2
              ICOR= ICOR+10
            ELSE
              IF(DGDV2.LT.0.d0) THEN
                  WRITE(6,628) JROT,EO
                  AFLAG= -1
                  GOTO 200
                  ENDIF
            ENDIF
          LTRY= 0
          KV= NF
          GOTO 100
          ENDIF
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IF(KV.NE.NF) THEN
c*** If last level found is not the desired one ...
          IF(INNR(KV).EQ.-1) THEN
c... Record vibrational level (if haven't already) for posterity.
              GV(KV)= EO
              INNR(KV)= INNER
              ENDIF
          ICOR= ICOR+1
          IF(ICOR.LE.20) THEN
c... Call subroutine using semiclassical methods to estimate correct energy
              CALL SCECOR(KV,NF,JROT,INNER,ICOR,IWR,EO,RH,BFCT,NDP,NCN,
     1                                              V,BMAX,VLIM,DGDV2)
              IF(EO.GT.VPMAX(NPMAX)) THEN
c... if estimated energy above highest barrier, set value below it
                  KV= 999
                  EO=  VPMAX(NPMAX) - 0.05d0*DGDV2
                  ENDIF
              GOTO 100
              ENDIF
c** If the calculated wavefunction is still for the wrong vibrational
c   level, then write out a warning return
          WRITE(6,630) NF,JROT
          KVMAX= NF-1
          ENDIF
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  200 IF(AFLAG.LT.0) THEN
c** If unable to find all KVMAX+1 levels requested, then return KVMAX as
c  v for the highest vibrational level actually found, and print out the
c  the energy of that level.
          KVMAX= KV
          IF(AWO.NE.0) WRITE(6,632) KVMAX, GV(KVMAX)
          ENDIF
      RETURN
c-----------------------------------------------------------------------
  602 FORMAT(/'  *** ALF ERROR ***'/4X,'Number of vib levels requested='
     1 ,i4,' exceeds internal ALF array dimension  NVIBMX=',i4)
  604 FORMAT(/' *** ALF ERROR ***   Find NO potential minima for   J=',
     1  i4)
  606 FORMAT(/'  ALF  finds onee potential minimum of',1PD15.7,
     1  '  at  R(1)=',0Pf9.6)
  608 FORMAT(/'  *** ALF WARNING ***'/4X,'There are',I3,'  potential ',
     1  A6,' in this potential. Stop searching after 10.')
  610 FORMAT(/'  *** ALF ERROR ***'/ 4X,'The potential turns over in the
     1 short range region at  R= ',G15.8)
  614 FORMAT(' Find',I3,'  potential minima:   Vmin=',8F11.3)
  616 FORMAT(15x,'at mesh points   R =',8f11.5)
  618 FORMAT(' Find',I3,'  potential maxima:   Vmax=',8F11.3)
  620 FORMAT(' *** So  STOP !!!!')
  622 FORMAT(/' ALF search finds next estimated trial energy  E(v=',I3,
     1 ')=',G15.8/8X,'lies above potential maximum or asymptote at  VMAX
     2=',G15.8)
  624 FORMAT(/' *** SCHRQ FAILS in ALF when searching for  v=',i3,
     1  ' J=',i3,'   with   EO=',f9.3/5x,'Check range and/or contact R.J
     2. Le Roy [leroy@uwaterloo.ca]')
  626 FORMAT(/' ALF successfully finds all (J=',i3,') vibrational levels
     1 up to   v= KVMAX=',I3)
  628 FORMAT(/' *** ERROR:   at   E(J=',i3,')=',f10.3,'  SCECOR  finds n
     1o Phase Integrals')
  630 FORMAT(4x,'ALF fails to find level   v=',i3,', J=',i3)
  632 FORMAT(' ALF finds the highest calculated level is  E(v=',I3,
     1  ')=',G15.8 /)
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE SCECOR(KV,KVLEV,JROT,INNER,ICOR,IWR,EO,RH,BFCT,NDP,
     1                                          NCN,V,BMAX,VLIM,DGDV2)
c** Subroutine calculates (approximate!) semiclassical estimate of 
c  dG/dv for level  v= KV  with energy  EO [cm-1]  on potential 
c  {V(i),i=1,NDP} (in 'internal BFCT units' {V[cm-1]*BFCT}), and uses
c  those results to estimate energy of level  KVLEV
c** If the 'clever' semiclassical procedure fails - try a brute force
c  step-by-step search, using alternately INNER & OUTER well starting
c** BMAX is internal barrier maximum energy for double-well case, 
c   and very large negative number for single-well potential
c** On return, negative DGDV2 signals error!  No phase integrals found
c=======================================================================
c    Version date:  12 June 2011
c***********************************************************************
      INTEGER I,II,I1,I2,I3,I4,IV1,IV2,INNER,ICOR,JROT,KV,KVB,KVLEV,
     1  KVDIF,NDP,NCN,IDIF,BRUTE,IB,IWR,NPMAX
      REAL*8 EO,DE0,RH,BFCT,ARG2,ARG3,EINT,VPH1,VPH2,DGDV1,DGDV2,DGDVM,
     1  DGDV2P,DGDVB,DGDVBP,EBRUTE,DEBRUTE,DE1,DE2,Y1,Y2,Y3,RT,ANS1,
     2  ANS2,XDIF,VLIM,BMAX,PNCN,PP1,V(NDP)
      SAVE BRUTE,EBRUTE,DEBRUTE,DGDVB
      DATA DGDVB/-1.d0/,KVB/-1/
c
      DGDV2= -1.d0
      EINT= EO*BFCT
      IF(KVLEV.EQ.0) DGDVB= -1.d0
      KVDIF= KVLEV- KV
      IF(ICOR.EQ.1) BRUTE= 0
      I3= NDP
      PNCN= DFLOAT(NCN-2)/DFLOAT(NCN+2)
      PP1= 1.d0/pNCN + 1.d0
c*** For Quasibound levels, first search inward to classically forbidden
      IF(EO.GT.VLIM) THEN
          PNCN=1.d0 
          PP1= 1.d0
          DO  I= NDP,1,-1
              I3= I
              IF(V(I).GT.EINT) GOTO 8
              ENDDO
          ENDIF
c*** First, search inward for outermost turning point
    8 DO  I= I3,1,-1
          I4= I
          IF(V(I).LT.EINT) GOTO 10
          ENDDO
c*** If never found an 'outer' turning point (e.g., above qbdd. barier)
c  then simply return with negative  DGDV2  as error flag
      RETURN
c... Now collect vibrational phase and its energy deriv. over outer well
   10 Y1= EINT- V(I4+1)
      Y2= EINT- V(I4)
      Y3= EINT- V(I4-1)
      CALL LEVQAD(Y1,Y2,Y3,RH,RT,ANS1,ANS2)
      ARG2= DSQRT(Y3)
      VPH2= 0.5d0*ARG2 + ANS2/RH
      DGDV2= 0.5d0/ARG2 + ANS1/RH
      DO  I= I4-2,1,-1
c... now, collect (v+1/2) and dv/dG integrals to next turning point ...
          II= I
          IF(V(I).GT.EINT) GO TO 12
          ARG3= ARG2
          ARG2= DSQRT(EINT - V(I))
          VPH2= VPH2+ ARG2
          DGDV2= DGDV2+ 1.d0/ARG2
          ENDDO
   12 I3= II+1
      Y1= EINT- V(I3-1)
      Y2= EINT- V(I3)
      Y3= EINT- V(I3+1)
      CALL LEVQAD(Y1,Y2,Y3,RH,RT,ANS1,ANS2)
      VPH2= (VPH2 - ARG2 - 0.5d0*ARG3 + ANS2/RH)/3.141592654d0
      DGDV2= DGDV2 -1.d0/ARG2 - 0.5d0/ARG3 + ANS1/RH
      DGDV2= 6.283185308d0/(BFCT*DGDV2)
c*** Next, search for innermost turning point
      DO  I= 1,NDP 
          I1= I
          IF(V(I).LT.EINT) GOTO 20
c... then collect vibrational phase and its energy deriv. over outer well
          ENDDO
   20 IF(I1.EQ.1) THEN
          WRITE(6,602) JROT,EO
          STOP
          ENDIF
      IF(I1.GE.I3) THEN
c*** For single-well potential or above barrier of double-well potential
          IF(IWR.GE.2) WRITE(6,600) ICOR,KV,JROT,EO,VPH2-0.5d0,DGDV2
          IF((KV.NE.(KVLEV-1)).AND.(DGDVB.GT.0.d0)) THEN
c... If got wrong level (KV not one below KVLEV) and NOT first call ...
              IF((EO-BMAX).GT.(2.d0*DGDV2)) THEN
c... 'Normal' case: use B-S plot area to estimate correct energy
                  DE0= KVDIF*(DGDV2- 0.5d0*(DGDV2-DGDVB)/DFLOAT(KV-KVB))
                  EO= EO+ DE0 
                  KV= KVB
c*** I don't recall what this is for, but it makes mischief in DPotFit!
cc                KVLEV= KV+1
                  RETURN
                ELSE
c... but close to barrier in double-well potential, switch to 'BRUTE'
                  BRUTE=BRUTE+ 1
                  DGDV1= DGDV2
                  XDIF= SIGN(1,KVDIF)  
                  GOTO 54
                ENDIF
              ENDIF
          IF(KV.EQ.0) THEN
c** Normally:  use B-S plot considerations to estimate next level energy
c... use harmonic estimate for v=1 energy
              EO= EO+ DGDV2
            ELSE
c... estimate Delta(G) based on linear Birge-Sponer
              DE0= 0.5d0*(3.d0*DGDV2 - DGDVB)
              IF((2.d0*DGDV2).GT.DGDVB) THEN
c... if linear Birge-Sponer predicts another level, then use it
                  EO= EO+ DE0
                ELSE
c... otherwise, use N-D theory extrapolation for next level...
                  DGDV2P= DGDV2**PNCN
                  DE0= (DGDV2P+DGDV2P-DGDVBP)
                  IF(DE0.GT.0.d0) THEN
                      DE0= (DE0**PP1- DGDV2P**PP1)/(PP1*(DGDV2P-DGDVBP))
                      EO= EO+ DE0
                    ELSE   
c... but if NDT(n=6) it predicts no more levels, go back to last DGDV
                      IF(IWR.GT.0) WRITE(6,604) KV,EO
                      ICOR= 100
  604 FORMAT(10x,'Find highest ound level is   E(v=',i3,')=',1PD18.10)
                    ENDIF
                ENDIF 
            ENDIF
          DGDVB= DGDV2
          DGDVBP= DGDVB**PNCN
          KVB= KV
          INNER= 0
          RETURN
          ENDIF
c
c*** For a double-well potential, collect vibrational phase and its 
c   energy derivative over the inner well
      Y1= EINT- V(I1-1)
      Y2= EINT- V(I1)
      Y3= EINT- V(I1+1)
      CALL LEVQAD(Y1,Y2,Y3,RH,RT,ANS1,ANS2)
      ARG2= DSQRT(Y3)
      VPH1= 0.5d0*ARG2 + ANS2/RH
      DGDV1= 0.5d0/ARG2 + ANS1/RH
      DO  I= I1+2,NDP
c... now, collect integral and count nodes outward to next turning point ...
          IF(V(I).GT.EINT) GO TO 22
          ARG3= ARG2
          ARG2= DSQRT(EINT - V(I))
          VPH1= VPH1+ ARG2
          DGDV1= DGDV1+ 1.d0/ARG2
          ENDDO
   22 I2= I-1
      Y1= EINT- V(I2+1)
      Y2= EINT- V(I2)
      Y3= EINT- V(I2-1)
      CALL LEVQAD(Y1,Y2,Y3,RH,RT,ANS1,ANS2)
      VPH1= (VPH1 - ARG2 - 0.5d0*ARG3 + ANS2/RH)/3.141592654d0
      DGDV1= DGDV1 -1.d0/ARG2 - 0.5d0/ARG3 + ANS1/RH
      DGDV1= 6.28318531d0/(BFCT*DGDV1)
      DGDVM= DGDV1*DGDV2/(DGDV1+DGDV2)
      IF(KVDIF.EQ.0) THEN
c** If already at level sought, return
          IF(IWR.GE.2) WRITE(6,610) KV,JROT,EO,VPH1-0.5d0,DGDV1,KVLEV,
     1                                           ICOR,VPH2-0.5d0,DGDV2
          RETURN
          ENDIF
c
c** Check whether looking for higher or lower level ...
      IDIF= SIGN(1,KVDIF)
      XDIF= IDIF
      IF((ICOR.GE.6).AND.((IABS(KVDIF).EQ.1).OR.(BRUTE.GT.0))) GOTO 50
c*** 'Conventional' semiclassical search for neared INNER or OUTER well level
      IF(INNER.LE.0) THEN
c... and current energy EO is for an outer-well level ...
          DE2= DGDV2*XDIF
          IV1= INT(VPH1+ 0.5d0)
          DE1= (DFLOAT(IV1) + 0.5d0 - VPH1)*DGDV1*XDIF
          IF(IWR.GE.2) WRITE(6,610) KV,JROT,EO,VPH1-0.5d0,DGDV1,KVLEV,
     1                                           ICOR,VPH2-0.5d0,DGDV2
   30     IF(DABS(DE1).LT.DABS(DE2)) THEN
              INNER= 1
              EO= EO+ DE1
              DE1= DGDV1*XDIF
            ELSE
              INNER= 0
              EO= EO+ DE2
            ENDIF
          KVDIF= KVDIF-IDIF
          IF(KVDIF.EQ.0) THEN
              RETURN
              ENDIF
          GOTO 30
          ENDIF
      IF(INNER.GT.0) THEN
c... and current energy EO is for an inner-well level ...
          DE1= DGDV1*XDIF
          IV2= INT(VPH2+ 0.5d0)
          DE2= (DFLOAT(IV2) + 0.5d0 - VPH2)*DGDV2*XDIF
          IF(IWR.GE.2) WRITE(6,610) KV,JROT,EO,VPH1-0.5d0,DGDV1,KVLEV,
     1                                           ICOR,VPH2-0.5d0,DGDV2
   40     IF(DABS(DE2).LT.DABS(DE1)) THEN
              INNER= 0
              EO= EO+ DE2
              DE2= DGDV2*XDIF
            ELSE
              INNER= 1
              EO= EO+ DE1
            ENDIF
          KVDIF= KVDIF-IDIF
          IF(KVDIF.EQ.0) THEN
              RETURN
              ENDIF
          GOTO 40
          ENDIF
   50 BRUTE= BRUTE+ 1
c*** Now .. Brute force search for desired level !
      IF(IWR.GE.2) WRITE(6,610) KV,JROT,EO,VPH1-0.5d0,DGDV1,KVLEV,
     1                                           ICOR,VPH2-0.5d0,DGDV2
   54 IF(BRUTE.EQ.1) THEN
c... in first brute-force step, use previous energy with opposite INNER
          EBRUTE= EO 
          IF(INNER.EQ.0) THEN
              INNER= 1
            ELSE
              INNER= 0
            ENDIF
          DEBRUTE= DMIN1(DGDV1,DGDV2)*XDIF*0.3d0
          RETURN
          ENDIF
      IB= BRUTE/2
c... in subsequent even steps, lower EO by DEBRUTE/10 for same INNER
      IF((IB+IB).EQ.BRUTE) THEN
          EBRUTE= EBRUTE+ DEBRUTE
          EO= EBRUTE
          RETURN
        ELSE
c... in subsequent odd steps, lower repeat previous EO with INNER changed
          IF(INNER.EQ.0) THEN
              INNER= 1
            ELSE
              INNER= 0
            ENDIF
          EO= EBRUTE
          RETURN
        ENDIF
c     RETURN
  600 FORMAT(' Single well  ICOR=',I2,':  E(v=',i3,',J=',I3,')=',f10.2,
     1 '  v(SC)=',F8.3,'  dGdv=',f8.3)
  602 FORMAT(/' *** ERROR ***  V(1) < E(J=',i3,')=',f10.2 )
  610 FORMAT(' Double well   E(v=',i3,', J=',I3,')=',f9.3,
     1 ':   v1(SC)=',F7.3,'   dGdv1=',f8.2/8x,'seeking  v=',I3,
     2 ' (ICOR=',I2,')',8x,':   v2(SC)=',F7.3,'   dGdv2=',f8.2 )
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE CDJOEL(EO,NBEG,NEND,NDIMR,BvWN,RH,WARN,V,WF0,RM2,RCNST)
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
c               COPYRIGHT 1994-2011  by  Robert J. Le Roy              +
c   Dept. of Chemistry, Univ. of Waterloo, Waterloo, Ontario, Canada   +
c    This software may not be sold or any other commercial use made    +
c      of it without the express written permission of the author.     +
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c  Authors: R.J. Le Roy & J. Tellinghuisen         Version of 05/02/2011
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Dimension:  potential arrays  and  vib. level arrays.
      INTEGER NDMINT
      INTEGER I,M,IPASS,M1,M2,NBEG,NEND,NDIMR,WARN
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
      DO  M= 2, 7
c** Kill nonsensical high-order CDCs (which can occur in double-well cases)
          IF(DABS(RCNST(M)).GT.DABS(RCNST(M-1))) THEN
              DO I= M, 7
                  RCNST(I)= 0.d0
                  ENDDO
              EXIT
              ENDIF
          ENDDO
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
c***** R.J. Le Roy  subroutine SCHRQ, last updated  12 June 2011 *******
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c                COPYRIGHT 2008-2011  by  Robert J. Le Roy             +
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
     2         KKV,KV,KVIN,LPRWF,M,MS,MSAVE,N,NBEG,NDN,NEND,NLINES,NPR
      REAL*8  BFCT,DE,DEP,DEPRN,DF,DOLD,DSOC,
     2        E,EEPS,EO,EPS,F,FX,GAMA,GI,GN,H,HT,PROD,PPROD,
     3        RATIN,RATOUT,RATST,RH,RINC,RMIN,RMINN,RR,RSTT,RWR(20),
     4        WF(N),SB,SI,SM,SN,SNEND,SPNEND,SRTGI,SRTGN,SWR(20),
     5        V(N),VLIM,VMAX,VMX,VPR,
     6        WKBTST,XEND,XPR,XPW,DXPW,Y1,Y2,Y3,YIN,YM,YOUT
      DATA RATST/1.D-9/,XPW/27.63d0/
      DATA NDN/15/
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      DXPW= XPW/NDN
      ICOR= 0
      KVIN= KV
      KV= -1
      RMINN= RMIN-RH
      GAMA= 0.d0
      VMAX= VLIM
      VMX= VMAX*BFCT
      H= RH
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
c** Start iterative loop; try to converge for up to 20 iterations.
      DO 90 IT= 1,20
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
                  IF(SM.GT.DXPW) EXIT
                  ENDDO
              IF(NEND.LT.N) NEND= NEND+ NDN
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
ccc               MS= M
                  Y2= Y2*SI
                  Y3= Y3*SI
                  SB= SB*SI
                  SI= 1.d0
                  ENDIF
              Y1= Y2
              Y2= Y3
c** Test for outermost maximum of wave function.
c... old S{max} matching condition - turning point works OK & is simpler.
ccc           IF((INNER.EQ.0).AND.(SI.LE.SB)) GO TO 32
c** Test for outermost well outer turning point 
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
          IF(INNER.GT.0) MSAVE= N
c** Actual outward integration loops start here
          DO  I= NBEG+2,MSAVE
              Y3= Y2+Y2-Y1+GI*SI
              GI= V(I)-E
              SI= Y3/(1.d0- HT*GI)
              WF(I)= SI
              IF(DABS(SI).GE.1.D+17) THEN
c** Renormalize to prevent overflow of  WF(I)  in classically forbidden
c  region where  V(I) .gt. E
                  SI= 1.d0/SI
                  DO  J= NBEG,I
                      WF(J)= WF(J)*SI
                      ENDDO
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
                  DO  J= NBEG,I
                      WF(J)= WF(J)*SI
                      ENDDO
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
c... note that by construction, at this point  WF(MSAVE)= 1.0
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
c** KV.ge.999  Option ... Search for highest bound level.  Adjust new
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
cc            IF(IT.EQ.20) THEN
cc                WRITE(6,617) IT,DEP
cc                WRITE(6,620) KVIN,JROT,ITER,DEPRN
cc                GO TO 100
cc                ENDIF
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
      IF(IWR.GE.2) WRITE(6,607) KV,JROT,EO,ITER,RR,NBEG,RATIN,INNER,
     1                                                     NEND,RATOUT
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
     1 6X,'M    R(M)   /WF(M)   /WF(M)   R(NEND) NBEG ITP1'/
     2  1X,99('-'))
  603 FORMAT(I4,1PD15.7,3D10.2,0PI7,F7.3,1P2D9.1,0PF8.2,I6,I5)
  604 FORMAT('   NOTE:  for  J=',I3,'   EO=',F12.4,' .ge. V(',i3,')=',
     1  F12.4)
  605 FORMAT(/' Solution of radial Schr. equation for   E(v=',I3,',J=',
     1  I3,') =',F15.7/2x,4('    R(I)   WF(I)   ')/2X,38('--') )
  606 FORMAT(2X,4(F8.3,F11.7))
  607 FORMAT('E(v=',I3,',J=',I3,')=',F11.4,1x,I3,' Iter   R(M)=',F6.3,
     1 '  WF(NBEG=',i6,')/WF(M)=',1PD8.1/ 37x,'INNER=',I2,6x,
     2 'WF(NEND=',i6,')/WF(M)=',D8.1)
  608 FORMAT(' *** SCHRQ Error:  E=',F9.2,' > V(',I6,')=',F9.2,
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
  619 FORMAT(1X,99('-'))
  620 FORMAT(' *** CAUTION for  v=',I3,'  J=',I3,"  SCHRQ doesn't conver
     1ge by  ITER=",I2,'  DE=',1PD9.2)
  701 FORMAT(/2x,'Level  v=',I3,'   J=',I3,'   E=',F12.4,' ,  wave funct
     1ion at',I6,' points.'/7x,'R(1-st)=',F12.8,'   mesh=',F12.8,
     2  '   NBEG=',I4,'   |LPRWF|=',I3)
  702 FORMAT((1X,4(0Pf9.4,1PD13.5)))
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE QBOUND(KV,JROT,E,EO,VMX,DSOC,V,RMIN,H,GB,GI,SB,SI,N,
     1  ITP3,IWR,IQTST,BFCT,IT)
c***********************************************************************
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
   10 II=J
c** Check that there is a classically allowed region inside this point
c  and determine height of barrier maximum.
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
c  Uses series expansions of Abramowitz & Stegun Eq.(10.4.3)
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
      NST= INT(DABS(EMSC)*1.D2)
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
  611 FORMAT(4X,'Log10(lifetime/sec)=',F10.5,';   Log10(width/cm-1)=',
     1 F10.5,'  dG/dv=',F7.2,'  V(max)=',f9.2,'(cm-1)')
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
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
c  on this interval.  **************************************************
c-----------------------------------------------------------------------
      REAL*8  A,ANS1,ANS2,B,C,CQ,H,HPI,R1,R2,RCQ,RR,RT,SL3,SLT,
     1        X0,Y1,Y2,Y3,ZT
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      DATA HPI/1.570796326794896D0/
      IF((Y1.GE.0).OR.(Y2.LT.0)) GO TO 99
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
      PARAMETER (NPINTMX=8000)
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
      DO 50 ITER= 1, 50
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
c!! First ... move derivative for constrained-parameter case
cc666 FORMAT(' For  IDAT=',I5,'  add PC(',I3,') =',1pD15.8,
cc   1  '  to PC(',0pI3,') =',1pD15.8)
                          IF(JFXP(J).GT.1) THEN
cc                            write(6,666) I,J,PC(J),JFXP(J),PC(JFXP(J))
                              PC(JFXP(J))= PC(JFXP(J))+ PC(J)
                              ENDIF
c  ... now continue collapsing partial derivative array
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
cc *** For 'Robust' fitting, adjust uncertainties here
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
                 IF(JFXP(J).GT.1) THEN
c** If this parameter constrained to equal some earlier parameter ....
                     PV(J)= PV(JFXP(J))
                     WRITE(6,668) J,JFXP(J),PV(J),ITER
                     ENDIF
  668 FORMAT(' Constrain  PV('i3,') = PV(',I3,') =',1pd15.8,
     1 '  on cycle',i3)

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
  604 FORMAT(' After Cycle #',i2,':   DRMSD=',1PD12.5,'    test(PS)=',
     1  1PD8.1,'   test(PU)=',D8.1)
  606 FORMAT(/' Effective',i3,'-cycle Cgce:  MAX{|change/unc.|}=',
     1  1PD8.1,' < 0.01   DRMSD=',D10.3)
  608 FORMAT(/' Full',i3,'-cycle convergence:  Max{|change/sens.|}=',
     1  1PD8.1,' < 1   DRMSD=',D10.2)
  610 FORMAT(/ ' !! CAUTION !! fit of',i4,' parameters to',I6,' data not
     1converged after',i3,' Cycles'/5x,'DRMS(deviations)=',1PD10.3,
     2 '    test(PS) =',D9.2,'    test(PU) =',D9.2/1x,31('**'))
  612 FORMAT((4x,'PV(',i4,') =',1PD22.14,' (+/-',D8.1,')    PS=',d8.1,
     1  '   PC=',d8.1))
  614 FORMAT(' =',39('==')/' Round Off  PV(',i4,')=',1PD21.13,' (+/-',
     1 D9.2,')    PS=',d9.2/10x,'fix it as ',D21.13,'  & refit:  DRMS(de
     2viations)=',D10.3)
  616 FORMAT(/i6,' data fit to',i5,' param. yields  DRMS(devn)=',
     1 1PD12.5:'  tst(PS)=',D8.1)
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
c  strength functions for V(r), phi(r), UA(r), UB(r), tA(r), tB(r) and 
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
      INTEGER I,J,JJ,ISTATE,LAM2,NBETAI, OSEL(NSTATEMX)
      REAL*8 FU,FLAM,RHT,RMAXT,RDVAL,RDVAL2,RDVALLD, PU(NPARMX),
     1  PT(NPARMX),CM(NPARMX,NPARMX) 
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
      RHT= RD(2,ISTATE)- RD(1,ISTATE)
      RMAXT= RD(NDATPT(ISTATE),ISTATE)
      DO I= 1,nPointSR(ISTATE),OSEL(ISTATE)
          RDVAL = RHT*DBLE(I-1) 
          WRITE(10,909) RDVAL,Vsr(I,ISTATE)
          END DO   
      WRITE(11,900) 'B(r) ', ISTATE, 'B(r) ', WRITFILE
      DO I= 1,nPointSR(ISTATE),OSEL(ISTATE)
          RDVAL = RHT*DBLE(I-1)
          WRITE(11,909) RDVAL,Bsr(I,ISTATE)
          END DO
      IF(NUA(ISTATE).GE.0) WRITE(12,900) 'UA(r)',ISTATE,'UA(r)',
     1                                                WRITFILE
      IF(NUB(ISTATE).GE.0) WRITE(13,900) 'UB(r)',ISTATE,'UB(r)',
     1                                                WRITFILE
      IF(NTA(ISTATE).GE.0) WRITE(14,900) 'tA(r)',ISTATE,'qA(r)',
     1                                                WRITFILE
      IF(NTB(ISTATE).GE.0) WRITE(15,900) 'tB(r)',ISTATE,'qB(r)',
     1                                                WRITFILE
      IF(NwCFT(ISTATE).GE.0) THEN
          WRITE(16,902) 'lambda(r)',ISTATE,'lambda(r)',WRITFILE
          LAM2= 2
          IF(IOMEG(ISTATE).GT.0) LAM2= 4*IOMEG(ISTATE)
          ENDIF
      NBETAI= POTPARF(ISTATE) - MAX(NSpow(ISTATE),NLpow(ISTATE))
      FU= 0.d0
      DO  I= 1,NDATPT(ISTATE),OSEL(ISTATE)
          DO J= 1,TOTPOTPAR
              PT(J) = 0.0d0
              END DO
          RDVAL = RD(I,ISTATE)
          RDVAL2= RDVAL*RDVAL
          RDVALLD= RDVAL**LAM2
c ... first ... the potential function itself ...
          DO  J= POTPARI(ISTATE),POTPARF(ISTATE)
              PT(J)= PU(J)*DVtot(J,I)
              ENDDO
          IF(PSEL(ISTATE).GT.0) 
     1           CALL MMCALC(POTPARI(ISTATE),POTPARF(ISTATE),PT,CM,FU)
          WRITE(10,910) RDVAL,VPOT(I,ISTATE),FU
c ... then the exponent function \phi(i)
          IF(PSEL(ISTATE).LE.5) THEN
              DO  J= NBETAI,POTPARF(ISTATE)
                  JJ= J-NBETAI
                  PT(J)= PU(J)*DBDB(JJ,I,ISTATE)
                  ENDDO
              CALL MMCALC(NBETAI,POTPARF(ISTATE),PT,CM,FU)
              WRITE(11,910) RDVAL,BETAFX(I,ISTATE),FU
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
          RDVAL= RMAXT + RHT*DBLE(I*OSEL(ISTATE))
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

c***********************************************************************
      double precision function Scalc(x,m,n,y,rKL,LMAX)
c** At the position 'x', evaluate the m'th Sm(x) function contributing 
c  the definition of the the natural cubic spline defined by 
c  function values at the  n  points  y(i) [i=1,n]
      INTEGER  LMAX,I,K,KK,M,N
      REAL*8  x,y1,y2,y(1:LMAX),rKL(1:LMAX,1:LMAX)
      k= 0
      kk= 0
      do i=2,n
c... select interval
          if ((x.gt.y(i-1)).and.(x.le.y(i)))  k=i
          end do
      if (x.lt.y(1)) then
          k=2
          kk=1
          end if
      if (x.gt.y(n)) then
          k=n
          kk=1
          end if
      if(x.eq.y(1)) k=2
      y1=y(k-1)
      y2=y(k)
      Scalc= 0.d0
      IF(kk.eq.0) 
     1    Scalc= rKL(m,k)*((y1-x)*(((y1-x)/(y1-y2))**2-1)/6)*(y1-y2)
     2         + rKL(m,k-1)*((x-y2)*(((x-y2)/(y1-y2))**2-1)/6)*(y1-y2)
      IF(k.EQ.m) Scalc= Scalc + (y1-x)/(y1-y2)
      IF(k-1.EQ.m) Scalc= Scalc + (x-y2)/(y1-y2)
c... Asen's original coding ...
cc       Scalc=ndirac(k,m)*A(x,y1,y2)+ndirac(k-1,m)*B(x,y1,y2)+
cc   +   C(x,y1,y2)*rKL(m,k)+D(x,y1,y2)*rKL(m,k-1)
cc       else
cc       Scalc=ndirac(k,m)*A(x,y1,y2)+ndirac(k-1,m)*B(x,y1,y2)
cc     A=(x1-z)/(x1-x2)
cc     B=(z-x2)/(x1-x2)
cc     C=((x1-z)*(((x1-z)/(x1-x2))**2-1)/6)*(x1-x2)
cc     D=((z-x2)*(((z-x2)/(x1-x2))**2-1)/6)*(x1-x2)
c... Asen's original coding ...
      end
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      double precision function Sprime(x,m,n,y,rKL,LMAX)
c** At the position 'x', evaluate the derivative w.r.t. x of the m'th 
c  Sm(x) function contributing the definition of the the natural cubic
c  spline defined by function values at the  n  points  y(i) [i=1,n]
      INTEGER i,k,kk,m,n,LMAX
      REAL*8 x,del,y1,y2,y(1:LMAX),rKL(1:LMAX,1:LMAX)
      k=0
      kk=0
      do i=2,n
          if((x.gt.y(i-1)).and.(x.le.y(i)))  k=i
          enddo
      if(x.lt.y(1)) then
          k=2
          kk=1
          end if
      if (x.gt.y(n)) then
          k=n
          kk=1
          end if
      if (x.eq.y(1)) k=2
      y1=y(k-1)
      y2=y(k)
      del=y1-y2
      Sprime= 0.d0
      if(kk.eq.0) Sprime= (del-3.d0*(y1-x)**2/del)*rKL(m,k)/6.d0 +
     1                        (3.d0*(x-y2)**2/del-del)*rKL(m,k-1)/6.d0
      IF(k-1.eq.m) Sprime= Sprime + 1.d0/del 
      IF(k.eq.m) Sprime= Sprime - 1.d0/del 
ccc     if(kk.eq.0) then
ccc         Sprim=ndirac(k-1,m)/del-ndirac(k,m)/del+
ccc  +                    (del-3*(y1-x)**2/del)*rKL(m,k)/6+
ccc  +                    (3*(x-y2)**2/del-del)*rKL(m,k-1)/6
ccc       else
ccc         Sprim=ndirac(k-1,m)/del-ndirac(k,m)/del
ccc       end if
      end

c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      subroutine Lkoef(n,x,A,LMAX)   
c*** Based on nespl subroutine          
      INTEGER LMAX
      INTEGER I,J,N,INDX(1:LMAX)
      REAL*8 X(1:LMAX),A(1:LMAX,1:LMAX),B(1:LMAX,1:LMAX), d
c
      DO  i= 1,LMAX
          DO  j= 1,LMAX
              A(i,j)= 0.d0
              B(i,j)= 0.d0
              ENDDO
          ENDDO
      A(1,1)= (x(3)-x(1))/3.d0
      A(1,2)= (x(3)-x(2))/6.d0
      do i= 2,n-3
          A(i,i-1)= (x(i+1)-x(i))/6.d0
          A(i,i)= (x(i+2)-x(i))/3.d0
          A(i,i+1)= (x(i+2)-x(i+1))/6.d0
          end do
      A(n-2,n-3)= (x(n-1)-x(n-2))/6.d0
      A(n-2,n-2)= (x(n)-x(n-2))/3.d0  
      do i= 1,n-2
          B(i,i)= 1.d0/(x(i+1)-x(i))
          B(i,i+1)= -1.d0/(x(i+2)-x(i+1))-1.d0/(x(i+1)-x(i))
          B(i,i+2)= 1.d0/(x(i+2)-x(i+1))
          end do  
      call ludcmp(A,n-2,LMAX,indx,d)
      do i= 1,n 
          call lubksb(A,n-2,LMAX,indx,B(1,i))
          end do 
      do i= 1,n-2
          do j= 1,n
              A(j,i+1)= B(i,j)
              end do
          end do 
      do i= 1,n
          A(i,1)= 0.0d0
          A(i,n)= 0.0d0
          end do
      end
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE ludcmp(a,n,np,indx,d)
      INTEGER n,np,indx(n),NMAX
      double precision d,a(np,np),TINY
      PARAMETER (NMAX= 500,TINY= 1.0e-20)
      INTEGER i,imax,j,k
      double precision aamax,dum,sum,vv(NMAX)
      d= 1.d0
      do  i= 1,n
          aamax= 0.d0
          do  j= 1,n
              if (abs(a(i,j)).gt.aamax) aamax= abs(a(i,j))
              enddo
          if (aamax.eq.0.) pause 'singular matrix in ludcmp'
          vv(i)= 1.d0/aamax
          enddo
      do  j= 1,n
          do  i= 1,j-1
              sum= a(i,j)
              do  k= 1,i-1
                  sum= sum-a(i,k)*a(k,j)
                  enddo
              a(i,j)= sum
              enddo
          aamax= 0.d0
          do  i= j,n
              sum= a(i,j)
              do  k= 1,j-1
                  sum= sum-a(i,k)*a(k,j)
                  enddo
              a(i,j)= sum
              dum= vv(i)*abs(sum)
              if (dum.ge.aamax) then
                  imax= i
                  aamax= dum
                  endif
              enddo
          if(j.ne.imax)then
              do  k= 1,n
                  dum= a(imax,k)
                  a(imax,k)= a(j,k)
                  a(j,k)= dum
                  enddo
              d= -d
              vv(imax)= vv(j)
              endif
          indx(j)= imax
          if(a(j,j).eq.0.)a(j,j)= TINY
              if(j.ne.n)then
                  dum= 1.d0/a(j,j)
                  do  i= j+1,n
                      a(i,j)= a(i,j)*dum
                      enddo
                  endif
          enddo
      return
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE lubksb(a,n,np,indx,b)
      INTEGER i,ii,j,ll, n,np,indx(n)
      double precision a(np,np),b(n), sum
      ii= 0
      do  i= 1,n
          ll= indx(i)
          sum= b(ll)
          b(ll)= b(i)
          if (ii.ne.0)then
              do  j= ii,i-1
                  sum= sum-a(i,j)*b(j)
                  enddo
            else if (sum.ne.0.) then
              ii= i
            endif
          b(i)= sum
          enddo
      do  i= n,1,-1
          sum= b(i)
          do  j= i+1,n
              sum= sum-a(i,j)*b(j)
              enddo
          b(i)= sum/a(i,i)
          enddo
      return
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
