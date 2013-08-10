c***********************************************************************
        PROGRAM DPotFit
c***********************************************************************
c** Program "D(iatomic)Pot(ential)Fit" (DPotFit) for performing least-
c   squares fits of diatomic spectral data to molecular potential
c   energy functions for one or multiple electronic states.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++   COPYRIGHT 2006-2013  by R.J. Le Roy, J.Y. Seto and Y. Huang   +++
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
      INCLUDE 'BLKPARAM.h'
      INCLUDE 'BLKBOB.h'
      INCLUDE 'BLKCOUNT.h'
c-----------------------------------------------------------------------
      CHARACTER*40 DATAFILE,MAKEPRED
      CHARACTER*24 WRITFILE,TVNAME(NPARMX)
      CHARACTER*27 FN4,FN6,FN7,FN8,FN10,FN11,FN12,FN13,FN14,FN15,FN16,
     1   FN17,FN20,FN22,FN30
cc   1                                                     ,FN32
      INTEGER*4 lnblnk
      INTEGER I,J,ISTATE,ISOT,CHARGE,IPV,MKPRED,PRINP,
     1 JTRUNC(NSTATEMX),EFSEL(NSTATEMX),PASok(NSTATEMX),
     2 NDAT(0:NVIBMX,NISTPMX,NSTATEMX)
      REAL*8 UCUTOFF,ZMASE,DECM(NSTATEMX)
c
      INTEGER NOWIDTHS
      COMMON /WIDTHBLK/NOWIDTHS
c                                                      
c** Parameters required for NLLSSRR.
c
      INTEGER NPTOT,CYCMAX,IROUND,ROBUST,LPRINT,SIROUND,NFPAR,uBv,
     1  IFXPV(NPARMX),SIFXPV(NPARMX)
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
c!!! REMOVE THIS OPTION - for such cases invoke a fake electronic state !
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
      READ(5,*)  AN(1), AN(2), CHARGE, NISTP, NSTATES, LPRINT, PRINP
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
      WRITE(6,601) NISTP
      DO  ISOT= 1,NISTP
c** Loop to read the mass numbers of the atoms in each of the isotopomers
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
          IF(I.EQ.0) WRITE(6,603)
          RSQMU(ISOT)= DSQRT(ZMASS(3,1)/ZMASS(3,ISOT))
          ENDDO
c... end of loop over isotopologues ....................................
      IF(CHARGE.NE.0) WRITE(6,597) CHARGE
      WRITE(6,599) DATAFILE,Fqb
      IF(AN(1).EQ.AN(2)) WRITE(6,604)
  599 FORMAT(/' Use experimental data input file:  ',a30/' Uncertainties
     1 for transitions involving quasibound levels modified to:'/20x,
     2  'SQRT{(u(i;exp)**2 + (',f5.2,'*width)**2}')
  597 FORMAT(1x,67('-')/' Since this is an ion with charge',SP,i3,
     1  ", use Watson's charge-modified reduced mass.")
  601 FORMAT(2X,'Input data for',I3,'  isotopomer(s)'/2X,16('**')/2X,
     1  '    Isotopomer       Mass of atom-1   Mass of atom-2    Reduced
     2 mass'/ 2X,'----------------- ',3('   --------------'))
  602 FORMAT(2X,A2,'(',I3,') - ',A2,'(',I3,')',3(3X,F14.9))
  603 FORMAT('  Note that   (Mass Number) = 0   causes the average atomi
     1c mass to be used.')
  604 FORMAT(' For electrically homonuclear molecules, BO correction fun
     1ctions are the same'/5x,'for both atoms, so only the first sets of
     2 correction coefficients'/5x,'UA(s) and TA(s) are used, and the ma
     3ss scaling factors are sums over'/5x,'the two individual atoms.')
c
      MKPRED= 0
      IF(DATAFILE.EQ.MAKEPRED) THEN
          MKPRED= 1
          ENDIF
c-----------------------------------------------------------------------
c  UCUTOFF   Neglect any input data with uncertainties > UCUTOFF (cm-1)
c
c  NOWIDTHS >  0  causes the program to ignore any tunneling widths in
c                 the data set and omit calculating partial derivatives
c                 of predissociation level widths w.r.t. potential param.
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
c
c  CYCMAX  sets an upper bound on the number of cycles to allowed in the 
c          least-squares fit subroutine NLLSSRR
c
c  uBv   defines whether (uBv > 0) or nor (uBv.LE.0) to compute the 
c    uncertainties in the calculated Gv & Bv values due to the fit 
c    uncertainties and write them to channel 17.
c=======================================================================
      READ(5,*)  UCUTOFF, NOWIDTHS, IROUND, ROBUST, CYCMAX, uBv 
c=======================================================================
      IF(IROUND.NE.0) WRITE(6,685) IABS(IROUND)
      IF(IROUND.GT.0) WRITE(6,686)
      IF(IROUND.LT.0) WRITE(6,687)
      IF(ROBUST.GT.0) THEN
          ROBUST= 2
          WRITE(6,596)
        ELSE
          WRITE(6,598)
        ENDIF
      WRITE(6,595) CYCMAX
  596 FORMAT( " Fit uses Watson's",' "Robust" data weighting [J.Mol/Spec
     1trosc. 219, 326 (2003)] '/20x,'1/[{unc(i)}^2 + {calc.-obs.}^2/3]')
  595 FORMAT(' Non-linear fits are allowed a maximum of  CYCMAX=', I4,' 
     1cycles')  
  598 FORMAT( ' Fit uses standard  1/[uncertainty(i)]**2  data weighting
     1')
  685 FORMAT(/' Apply "Sequential Rounding & Refitting" at digit-',
     1  i1,' of the (local) parameter')
  686 FORMAT(4x,'uncertainty, selecting remaining parameter with largest
     1 relative uncertainty')
  687 FORMAT(4x,'uncertainty, proceeding sequentially from the LAST para
     1meter to the FIRST.')

c!!!
cc    IF(NOWIDTHS.LE.0) THEN
cc        WRITE(FN32,'(2A)') WRITFILE(1:lnblnk(WRITFILE)),'.32'
cc        OPEN(UNIT=32, FILE= FN32)
cc        ENDIF
c!!!  
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
c             is to be scaled by [J(J+1) + 2] (special Li2 case)
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
          CALL READPOT(ISTATE,SLABL)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** These statements construct and define the names of output files
c   associated with WRITE's to channels 10-16 used by the program.
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
      IF(uBv.GE.0) THEN
c** If uBv > 0,  define the name of the output file for Bv & Gv uncert
          WRITE(FN17,'(2A)') WRITFILE(1:lnblnk(WRITFILE)),'.17'
          OPEN(UNIT=17,FILE=FN17)
          ENDIF
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
c** For Pashov-exponent MLR count starts with 1 for y_p = 1.0
          IF(APSE(ISTATE).GT.0) J=1
          DO I= J,Nbeta(ISTATE)
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
  626 FORMAT(/' *** Dimension Error *** [(total No. Hamiltonian parmaete
     1rs)=',i4,'] >  HPARMX=',i4)
  638 FORMAT(' State ',A2,'  Energy Convergence criterion EPS is',
     1  1PD8.1,' cm-1')
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
          CALL NLLSSRR(COUNTOT,NPTOT,NPARMX,CYCMAX,IROUND,ROBUST,LPRINT,
     1             IFXPV,FREQ,UFREQ,DFREQ,PV,PU,PS,CM,TSTPS,TSTPU,DSE)
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
      CALL NLLSSRR(COUNTOT,NPTOT,NPARMX,CYCMAX,IROUND,ROBUST,LPRINT,
     1             IFXPV,FREQ,UFREQ,DFREQ,PV,PU,PS,CM,TSTPS,TSTPU,DSE)
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
          CALL NLLSSRR(COUNTOT,NPTOT,NPARMX,CYCMAX,IROUND,ROBUST,LPRINT,
     1             IFXPV,FREQ,UFREQ,DFREQ,PV,PU,PS,CM,TSTPS,TSTPU,DSE)
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
                  CALL FUNUNC(ISTATE,WRITFILE,PU,CM)
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
      CALL DIFFSTATS(NSTATES,NFPAR,ROBUST,MKPRED)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IF(uBv.GT.0) THEN
c** If desired, calculate Bv uncertainties
          CALL UNCBV(NPTOT,PV,PU,CM)
          ENDIF
      STOP
  900 FORMAT(5X,G18.8,5X,G18.8)
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
