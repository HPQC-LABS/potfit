c======================================================================c
c|..optional...Register...optional...Register...optional...Register...|c
c|--------------------------------------------------------------------|c
c| You have choosen to download the following source code for my      |c
c| Fortran program dPotFit. I would appreciate it if you would please |c
c| go to the www address                                              |c
c|            http://scienide2.uwaterloo.ca/~rleroy/dPotFit16/   and  |c
c| fill in the registration form there if you wish to be accessible   |c
c| so that I can send you possible future updates and/or corrections  |c
c| for this code.  This address list will be held securely by me and  |c
c| used for no other purpose................. Robert J. Le Roy .......|c
c|..Register...optional....Register...optional...Register...optional..|c
c======================================================================c

c***********************************************************************
        PROGRAM dPotFit16
c***********************************************************************
c** Program "D(iatomic)Pot(ential)Fit" (dPotFit) for performing least-
c   squares fits of diatomic spectral data to molecular potential
c   energy functions for one or multiple electronic states.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++++++++++   COPYRIGHT 2006-2016  by R.J. Le Roy, +++++++++++++++++++++
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
c    one or more electronic states and one or more isotopologues, to
c    parameters defining the observed levels of each state.
c=======================================================================
c** Dimensioning parameters intrinsic to the program are input through
c          'arrsizes.h'
c** Parameters characterizing the problem and governing the fits are
c   read on  channel-5  while the raw data are read on  channel-4 .
c   Principle output goes to  channel-6  while higher channel numbers
c   are used for secondary or more detailed/voluminous output.
c***********************************************************************
cc    INCLUDE 'arrsizes.h'
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c** 'Block' Data Utility routine named: 'arrsizes.h' that governs 
c    array dimensioning in program  dPotFit
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER NISTPMX,NPARMX,NbetaMX,NBOBMX,HPARMX,NDATAMX,
     1  NVIBMX,NBCMX,NSTATEMX,NPNTMX,NROTMX,NCMMAX
c*  NISTPMX  is the maximum number of isotopomers allowed for fit
      PARAMETER (NISTPMX =  12)
c*  NSTATEMX  is maximum no. of electronic states which can be
c             simultaneously fitted to
      PARAMETER (NSTATEMX = 4)
c*  NPARMX  is the largest number of free parameters allowed for fit
c  Since FS origins may be parameters, this is also max. no, data bands
      PARAMETER (NPARMX  = 8000)
c*  NbetaMX  is the largest number of exponent parameters allowed for fit
      PARAMETER (NbetaMX  = 40)
c*  NBOBMX-1  is the highest-order polynomial expansion allowed for the
c      adiabatic or centrifugal Born-Oppenheimer breakdown functions, or 
c      the Lambda-doubling or 2\Sigma splitting radial strength functions
      PARAMETER (NBOBMX  = 15)
c*  HPARMX  is the largest number of Hamiltonian parameters of all types
c    (potential energy, BOB. etc.) for all states.
c           HPARMX >= NSTATEMX*[5 + (NbetaMX+1) + 5*(NBOBMX+1)]
      PARAMETER (HPARMX= NSTATEMX*(5 + (NbetaMX+1) + 5*(NBOBMX+1)))
cc    PARAMETER (HPARMX = 300)
c*  NDATAMX  is largest No. of individual data which may be considered
      PARAMETER (NDATAMX = 35000)
c*  NVIBMX  is the maximum number of vibrational levels of a single
c           state for which data are to be considered
      PARAMETER (NVIBMX    = 200)
** NBCMX  is the maximum number of band constants per vib level to be
c         allowed when doing band constant fits (PSEL= -1)
      PARAMETER (NBCMX = 8) 
c*  NPNTMX  is the largest number of potential data points that can be
c           stored in a single 1D array
      PARAMETER (NPNTMX = 90000)
c*  NROTMX  is the highest order of rotational constants calculated and
c            used for estimating level energies
      PARAMETER (NROTMX = 7)
c*  NCMMAX is the largest number of Cm terms in the MLR or DELR 
c            long-range potential
      PARAMETER (NCMMAX = 12)
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cc    INCLUDE 'BLKISOT.h'
c=======================================================================
c** Isotope/isotopologue numbers, masses & BOB mass scaling factors
c** Array ZK carries about the band constants for all levels of all ISOT
      INTEGER NISTP,NDUNMX,AN(2),MN(2,NISTPMX)
c** NDUNMX is a dummy parameter reqd. for portability of READATA
      PARAMETER (NDUNMX=0)
      REAL*8  ZMASS(3,NISTPMX),RSQMU(NISTPMX),RSQMUP(0:NDUNMX,NISTPMX),
     1 RMUP(0:9,NISTPMX),ZMUA(NISTPMX,NSTATEMX),ZMUB(NISTPMX,NSTATEMX),
     2 ZMTA(NISTPMX,NSTATEMX),ZMTB(NISTPMX,NSTATEMX),
     3  ZK(0:NVIBMX,0:NROTMX,NISTPMX,NSTATEMX)
c
      COMMON /BLKISOT/ZMASS,RSQMU,RSQMUP,RMUP,ZMUA,ZMUB,ZMTA,ZMTB,ZK,
     1  NISTP,AN,MN
c=======================================================================
cc    INCLUDE 'BLKDATA.h'
c=======================================================================
c** Type statements & common block for data
      REAL*8  FREQ(NDATAMX),UFREQ(NDATAMX),DFREQ(NDATAMX),TEMP(NDATAMX),
     1                                               YUNC(NDATAMX),Fqb
      INTEGER  COUNTOT,NFS1,NFSTOT,NBANDTOT,IB(NDATAMX),JP(NDATAMX),
     1 JPP(NDATAMX),VP(NPARMX),VPP(NPARMX),EFP(NDATAMX),EFPP(NDATAMX),
     2 TVUP(NDATAMX),TVLW(NDATAMX),FSBAND(NPARMX),IFXFS(NPARMX),
     3 NFS(NPARMX),IEP(NPARMX),IEPP(NPARMX),ISTP(NPARMX),
     4 IFIRST(NPARMX),ILAST(NPARMX),NTV(NSTATEMX,NISTPMX),FSsame,
     5 NTRANS(NPARMX),IBB(NISTPMX,NSTATEMX,9,NPARMX),JMIN(NPARMX),
     6 JMAX(NPARMX)
      CHARACTER*2 NAME(2) 
      CHARACTER*3 SLABL(-6:NSTATEMX)
      CHARACTER*30 BANDNAME(NPARMX)
      COMMON /DATABLK/Fqb,FREQ,UFREQ,YUNC,DFREQ,TEMP,COUNTOT,NFS1,
     1 NFSTOT,NBANDTOT,IB,JP,JPP,VP,VPP,EFP,EFPP,TVUP,TVLW,FSBAND,IFXFS,
     2 NFS,IEP,IEPP,ISTP,IFIRST,ILAST,NTV,FSsame,
     3 NTRANS,IBB,JMIN,JMAX,NAME,SLABL,BANDNAME
c=======================================================================
cc    INCLUDE 'BLKPOT.h'
c=======================================================================
c** Effective adiabatic radial potential variables.
      INTEGER BOBCN(NSTATEMX),PSEL(NSTATEMX),MAXMIN(NSTATEMX),
     1 IOMEG(NSTATEMX),Nbeta(NSTATEMX),APSE(NSTATEMX),IFXDE(NSTATEMX),
     2 IFXRE(NSTATEMX),IFXCm(NCMMax,NSTATEMX),
     3 IFXBETA(0:NbetaMX,NSTATEMX),NDATPT(NSTATEMX),NCMM(NSTATEMX),
     4 MMLR(NCMMax,NSTATEMX),nPB(NSTATEMX),nQB(NSTATEMX),pAD(NSTATEMX),
     5 qAD(NSTATEMX),LRad(NSTATEMX),pNA(NSTATEMX),qNA(NSTATEMX),
     6 Pqw(NSTATEMX),IVSR(NSTATEMX),IDSTT(NSTATEMX)
c
      REAL*8 DE(NSTATEMX),RE(NSTATEMX),BETA(0:NbetaMX,NSTATEMX),
     1 yqBETA(NbetaMX,NSTATEMX),BETAFX(NPNTMX,NSTATEMX),RH(NSTATEMX),
     2 RMIN(NSTATEMX),RMAX(NSTATEMX),VLIM(NSTATEMX),EPS(NSTATEMX),
     3 betaINF(NSTATEMX),AGPEF(NSTATEMX),BGPEF(NSTATEMX),
     4 CmVAL(NCMMax,NSTATEMX),CmEFF(NCMMax,NSTATEMX),rhoAB(NSTATEMX),
     5 AA(NSTATEMX),BB(NSTATEMX),RREF(NSTATEMX),ASO(NSTATEMX),
     6 R01(NSTATEMX),Q12(NSTATEMX),RD(NPNTMX,NSTATEMX),
     7 VPOT(NPNTMX,NSTATEMX),dCmA(NCMMax,NSTATEMX),dCmB(NCMMax,NSTATEMX)
c
      COMMON /BLKPOT/DE,RE,BETA,yqBETA,BETAFX,RH,RMIN,RMAX,VLIM,EPS,
     1 betaINF,AGPEF,BGPEF,CmVAL,CmEFF,rhoAB,AA,BB,RREF,ASO,R01,Q12,RD,
     2 VPOT,dCmA,dCmB, BOBCN,PSEL,MAXMIN,IOMEG,Nbeta,APSE,IFXDE,IFXRE,
     3 IFXCm,IFXBETA,NDATPT,NCMM,MMLR,nPB,nQB,pAD,qAD,LRad,pNA,qNA,Pqw,
     4 IVSR,IDSTT
c=======================================================================
cc    INCLUDE 'BLKPARAM.h'
c=======================================================================
c** Parameters and count-labels for band constant (PSEL=-1) or term
c   value (PSEL=-2) fits
      REAL*8 TVALUE(NPARMX),ZBC(0:NVIBMX,0:NROTMX,NISTPMX,NSTATEMX),
     1 ZQC(0:NVIBMX,0:NROTMX,NISTPMX,NSTATEMX)
c
      INTEGER NSTATES,NTVALL(0:NSTATEMX),NTVI(NSTATEMX),NTVF(NSTATEMX),
     1 VMIN(NSTATEMX,NISTPMX),VMAX(NSTATEMX,NISTPMX),JTRUNC(NSTATEMX),
     2 EFSEL(NSTATEMX),NBC(0:NVIBMX,NISTPMX,NSTATEMX),
     3 NQC(0:NVIBMX,NISTPMX,NSTATEMX),
     4 BCPARI(0:NVIBMX,NISTPMX,NSTATEMX),
     5 BCPARF(0:NVIBMX,NISTPMX,NSTATEMX),
     6 QCPARI(0:NVIBMX,NISTPMX,NSTATEMX),
     7 QCPARF(0:NVIBMX,NISTPMX,NSTATEMX)
      COMMON /BLKPARAM/TVALUE,ZBC,ZQC,NSTATES,NTVALL,NTVI,NTVF,VMIN,
     1      VMAX,JTRUNC,EFSEL,NBC,NQC,BCPARI,BCPARF,QCPARI,QCPARF
c=======================================================================
cc    INCLUDE 'BLKBOB.h'
c=======================================================================
c** Born-Oppenheimer Breakdown & doubling function parameters.
c**                       March 16 2012
c=======================================================================
      INTEGER NUA(NSTATEMX),NUB(NSTATEMX),NTA(NSTATEMX),NTB(NSTATEMX),
     1  IFXUA(0:NBOBMX,NSTATEMX),IFXUB(0:NBOBMX,NSTATEMX),
     2  IFXTA(0:NBOBMX,NSTATEMX),IFXTB(0:NBOBMX,NSTATEMX),
     3  NwCFT(NSTATEMX),IFXwCFT(0:NBOBMX,NSTATEMX),efREF(NSTATEMX)
c
      REAL*8 UA(0:NBOBMX,NSTATEMX),UB(0:NBOBMX,NSTATEMX),
     1  TA(0:NBOBMX,NSTATEMX),TB(0:NBOBMX,NSTATEMX),
     2   wCFT(0:NBOBMX,NSTATEMX)
c
      COMMON /BLKBOB/UA,UB,TA,TB,wCFT,NUA,NUB,NTA,NTB,NwCFT,
     1  IFXUA,IFXUB,IFXTA,IFXTB,IFXwCFT,efREF
c=======================================================================
cc    INCLUDE 'BLKCOUNT.h'
c=======================================================================
c    Block data file  BLKCOUNT.h
c=======================================================================
c** Counters for numbers of potential parameters of different types for 
c   each state
      INTEGER  TOTPOTPAR,POTPARI(NSTATEMX),POTPARF(NSTATEMX),
     1  UAPARI(NSTATEMX),UAPARF(NSTATEMX),UBPARI(NSTATEMX),
     2  UBPARF(NSTATEMX),TAPARI(NSTATEMX),TAPARF(NSTATEMX),
     3  TBPARI(NSTATEMX),TBPARF(NSTATEMX),LDPARI(NSTATEMX),
     4  LDPARF(NSTATEMX),HPARF(NSTATEMX),OSEL(NSTATEMX)
c
      COMMON /BLKCOUNT/TOTPOTPAR,POTPARI,POTPARF,UAPARI,UAPARF,UBPARI,
     1  UBPARF,TAPARI,TAPARF,TBPARI,TBPARF,LDPARI,LDPARF,HPARF,OSEL
c=======================================================================
c-----------------------------------------------------------------------
c** Common block for partial derivatives of potential at the one distance RDIST
c   and HPP derivatives for uncertainties
      REAL*8 dVdPk(HPARMX),dDe(0:NbetaMX),dDedRe
      COMMON /dVdPkBLK/dVdPk,dDe,dDedRe
c=======================================================================
      CHARACTER*40 DATAFILE,MAKEPRED
      CHARACTER*24 WRITFILE,TVNAME(NPARMX)
      CHARACTER*27 FN4,FN6,FN7,FN8,FN10,FN11,FN12,FN13,FN14,FN15,FN16,
     1   FN17,FN20,FN22,FN30
cc   1                                                     ,FN32
      INTEGER*4 lnblnk
      INTEGER I,J,ISTATE,ISOT,CHARGE,hCHARGE1,hCHARGE2,CHARGE3,IPV,
     1 MKPRED,PRINP,PASok(NSTATEMX),NDAT(0:NVIBMX,NISTPMX,NSTATEMX),
     2 NTVSTATE,NTVSSTAT,NTVSTOT,VMAXIN(NSTATEMX)
      REAL*8 UCUTOFF,ZMASE,HZMASE,DECM(NSTATEMX)
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
c** Set type statements for (unused) MASSES variables.
c
      CHARACTER*2  CATOM
      INTEGER GELGS(2,NISTPMX),GNS(2,NISTPMX)
      REAL*8 zIP,ABUND
c------------------------------------------------------------------------
c------------------------------------------------------------------------
c*** Common Block info for  fununc  calculations ***********************     
      REAL*8 Rsr(NPNTMX,NSTATEMX),Vsr(NPNTMX,NSTATEMX),
     1                                            Bsr(NPNTMX,NSTATEMX)
      INTEGER nPointSR(NSTATEMX)
      COMMON /VsrBLK/Rsr,Vsr,Bsr,nPointSR
c
      REAL*8 Rlr(NPNTMX,NSTATEMX),Plr(NPNTMX,NSTATEMX),
     1                                            Blr(NPNTMX,NSTATEMX)
      INTEGER nPointLR(NSTATEMX)
      COMMON /PlrBLK/Rlr,Plr,Blr,nPointLR
c-----------------------------------------------------------------------
c************************** misc. other variables **********************
      REAL*8 RDIST,VDIST,BETADIST,RMAXT,RHT,RHL
      INTEGER NCNN,NBCTOT
      DATA ZMASE /5.4857990946D-04/    !! 2010 physical constants d:mohr12
      DATA MAKEPRED/'MAKEPRED'/
c=======================================================================
      HZMASE= 0.5d0*ZMASE
      SLABL(-6)= '   '             !! data type not yet defined
      SLABL(-5)='VAC'              !! Accoustic Virial Coefficient
      SLABL(-4)='VIR'              !! Pressure Virial Coefficients  
      SLABL(-3)='VVV'              !! potential function values
      SLABL(-2)='WID'              !! tunneling level widths
      SLABL(-1)='PAS'              !! Photo-Association binding energies
      SLABL(0)='FLS'               !! fluorescence series
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
c  If(CHARGE.ne.0) use Watson's(JMS 1980) charge-modified reduced mass
c  (default case), OR assign gained/lost electrons masses (in units of
c   {m_e/2}) to one particle or the other using integers hCHARGE1 & hCHARGE2
c
c  NISTP   is the number of isotopologues to be simultaneously considered.
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
c**  IF |PRINP|=2 READ  title BANDNAME(IBAND) for each Band/Series on 1'st
c    line of input for that series and print it at the end of 'summary'
c    file for that series in the Channels 6 & 8 output.
c
c** For |CHARGE|.ne. 0: option to distribute missing/added e^- mass(es) 
c   Read # half-electron-masses to be added to/subtracted from standard
c   atomic masses to create standard 2-body reduced mass  m1*m2/(m1+m2)
c   For Watson's charge-adjusted reduced mass, set  hCHARGE1= hCHARGE2= 0
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
      IF(IABS(CHARGE).NE.0) READ(5,*) hCHARGE1, hCHARGE2
      READ(5,*)  DATAFILE
      READ(5,*)  WRITFILE
c=======================================================================
c** Now construct and define the names of output files associated with 
c   WRITE's to channels 6, 7, 8, 20, 22 & 30 used by the program.
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
c for a molecular ion, printout re. placement of +/- e^- mass(es)
      CHARGE3= 0
      IF(CHARGE.NE.0) THEN
          CHARGE3= hCHARGE1 + hCHARGE2
          IF((hCHARGE1.NE.0).OR.(hCHARGE2.NE.0)) THEN
c** If wish to add/subtract e- mass(es) to atomic mass of ions .....
              IF(CHARGE3.NE.2*CHARGE) THEN
c,,, if adding particle charges don't give total charge ... ERROR & STOP
                  WRITE(6,605) hCHARGE1,hCHARGE2,CHARGE
                  STOP
                  ENDIF
              WRITE(6,606) hCHARGE1,hCHARGE2,hCHARGE1,hCHARGE2
            ELSE
              WRITE(6,607) 
            ENDIF  
          ENDIF
      WRITE(6,601) NISTP
  606 FORMAT('  Reduced masses below are based on atoms 1 & 2 with charg
     1es (',SP,I2,'/2) and (',I2,'/2),'/8x,'respectively, with subtracti
     2on/addition of',SS,I2,' and',I2,' half-electron masses.'/)
  605 FORMAT(' *** ERROR *** atomic charges',SP,I3,'/2  and',I3,"/2   do
     1n't add up to total  CHARGE=",I3/10x,' !!! so STOP !!!!')
  607 FORMAT("  Reduced masses are Watson's charge-modified reduced mass
     1 for diatomic ions"/)
c
      DO  ISOT= 1,NISTP
c** Loop to read the mass numbers of the atoms in each of the isotopologues
c  MN(i,ISOT)  is the mass number for atom with atomic number AN(i)
c       [NOTE: be sure order of MN values consistent with that of AN's].
c       Choosing it .ne. value for some known isotope if that species
c       causes the average atomic mass to be used.
c=======================================================================
          READ(5,*) MN(1,ISOT), MN(2,ISOT)
c=======================================================================
          I= MIN(I,MN(1,ISOT),MN(2,ISOT))
          CALL MASSES(AN(1),MN(1,ISOT),CATOM,GELGS(1,ISOT),
     1                    GNS(1,ISOT),ZMASS(1,ISOT),ABUND)
          IF(ISOT.EQ.1) NAME(1)= CATOM
          CALL MASSES(AN(2),MN(2,ISOT),CATOM,GELGS(2,ISOT),
     1                    GNS(2,ISOT),ZMASS(2,ISOT),ABUND)
          IF(ISOT.EQ.1) NAME(2)= CATOM
          IF(CHARGE3.EQ.0) THEN     !! Watson charge modified mass
              ZMASS(3,ISOT)= (ZMASS(1,ISOT)*ZMASS(2,ISOT))/
     1                    (ZMASS(1,ISOT)+ZMASS(2,ISOT)-CHARGE*ZMASE)
            ELSE                          !! standard 2-body mass
              IF(CHARGE.NE.0) THEN             !! adjust masses for ion
                  ZMASS(1,ISOT)= ZMASS(1,ISOT) - hCHARGE1*HZMASE
                  ZMASS(2,ISOT)= ZMASS(2,ISOT) - hCHARGE2*HZMASE
                  ENDIF
              ZMASS(3,ISOT)= ZMASS(1,ISOT)*ZMASS(2,ISOT)/
     2                                 (ZMASS(1,ISOT) + ZMASS(2,ISOT))
            ENDIF
          WRITE(6,602) NAME(1),MN(1,ISOT),NAME(2),MN(2,ISOT),
     1                                          (ZMASS(J,ISOT),J=1,3)
          IF(I.EQ.0) WRITE(6,603)
          RSQMU(ISOT)= DSQRT(ZMASS(3,1)/ZMASS(3,ISOT))
          ENDDO
c... end of loop over isotopologues ....................................
ccc   IF(CHARGE.NE.0) WRITE(6,597) CHARGE
      WRITE(6,599) DATAFILE,Fqb
      IF(AN(1).EQ.AN(2)) WRITE(6,604)
  599 FORMAT(/' Use experimental data input file:  ',a30/' Uncertainties
     1 for transitions involving quasibound levels modified to:'/20x,
     2  'SQRT{(u(i;exp)**2 + (',f5.2,'*width)**2}')
cc597 FORMAT(1x,67('-')/' Since this is an ion with charge',SP,i3,
cc   1  ", use Watson's charge-modified reduced mass.")
  601 FORMAT(2X,'Input data for',I3,'  isotopologues(s)'/2X,16('**')/2X,
     1  '    Isotopologues    Mass of atom-1   Mass of atom-2    Reduced
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
c              Watson [J.Mol.Spectrosc. 219, 326 (2003)]) to be used
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
c
      DO  ISTATE= 1,NSTATES 
c-----------------------------------------------------------------------
c** Read parameters to characterize state & possibly restrict data used
c  SLABL(s)  is a 3-character alphameric label enclosed in single quotes
c            to identify the electronic state; e.g., 'XSG', 'A1P', ... etc.
c  IOMEG(s) .GE.0  is electronic angular momentum of singlet state with
c                    projection quantum number  Lambda= IOMEG
c  IOMEG(s) .EQ. -1  if it indicates a doublet SIGMA electronic state
c                 [other spin multiplets not yet coded]
c  IOMEG(s) .EQ. -2  indicated that the centrifugal potential strength 
c             factor is  [J(J+1) + 2]  (special Li2 case)
c  V(MIN/MAX)(s) Neglect data for electronic state vibrational levels
c                outside the range  VMIN  to  VMAX.
c  JTRUNC(s)     data with J > JTRUNC are not included in the fit.
c  EFSEL(s)  allows a user to consider data for:
c          * ONLY the e-parity levels of this state, if EFSEL > 0
c          * ONLY the f-parity levels of this state, if EFSEL < 0
c          * BOTH e- and f-parity levels of this state, if EFSEL = 0
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
              VMAXIN(ISTATE)= 1
              IF(VMAX(ISTATE,1).LT.0) THEN 
                  VMAXIN(ISTATE)= VMAX(ISTATE,1)
c** If desired, read separate upper bound level for each isotopologue
c=======================================================================
                  READ(5,*) (VMAX(ISTATE,ISOT), ISOT= 1, NISTP)
c=======================================================================
                  ENDIF
              ENDIF
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          CALL READPOT(ISTATE,SLABL)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** These statements construct and define the names of output files
c   associated with WRITE's to channels 10-16 used by the program.
          IF(OSEL(ISTATE).NE.0) THEN
              WRITE(FN10,'(2A)') WRITFILE(1:lnblnk(WRITFILE)),'.10'
              OPEN(UNIT=10,FILE=FN10)
cc            WRITE(FN11,'(2A)') WRITFILE(1:lnblnk(WRITFILE)),'.11'
cc            OPEN(UNIT=11,FILE=FN11)
              IF(OSEL(ISTATE).LT.0) THEN
                  IF(NUA(ISTATE).GE.0) THEN
                      WRITE(FN12,'(2A)') WRITFILE(1:lnblnk(WRITFILE)),
     1                                                           '.12'
                      OPEN(UNIT=12, FILE=FN12)
                      ENDIF
                  IF(NUB(ISTATE).GE.0) THEN
                      WRITE(FN13,'(2A)') WRITFILE(1:lnblnk(WRITFILE)),
     1                                                           '.13'
                      OPEN(UNIT=13,FILE=FN13)
                      ENDIF
                  IF(NTA(ISTATE).GE.0) THEN
                      WRITE(FN14,'(2A)') WRITFILE(1:lnblnk(WRITFILE)),
     1                                                           '.14'
                      OPEN(UNIT=14,FILE=FN14)
                      ENDIF
                  IF(NTB(ISTATE).GE.0) THEN
                      WRITE(FN15,'(2A)') WRITFILE(1:lnblnk(WRITFILE)),
     1                                                           '.15'
                      OPEN(UNIT=15,FILE=FN15)
                      ENDIF
                  IF(NwCFT(ISTATE).GE.0) THEN
                      WRITE(FN16,'(2A)') WRITFILE(1:lnblnk(WRITFILE)),
     1                                                           '.16'
                      OPEN(UNIT=16,FILE=FN16)
                      ENDIF
                  ENDIF
              ENDIF
          PASok(ISTATE)= 1
          IF(PSEL(ISTATE).EQ.6) PASok(ISTATE)= 0
c** Call VGEN to generate  betaINF  value for output in WRITEPOT
          IF(PSEL(ISTATE).EQ.2) THEN
              POTPARI(ISTATE)= 1
              CALL VGEN(ISTATE,1.0d0,VDIST,BETADIST,0)
              ENDIF
          ENDDO
      IF(uBv.GT.0) THEN
c** If uBv > 0,  define the name of the output file for Bv & Gv uncert
          WRITE(FN17,'(2A)') WRITFILE(1:lnblnk(WRITFILE)),'.17'
          OPEN(UNIT=17,FILE=FN17)
          ENDIF
c** Now write summary of the initial potential parameters for each state
      CALL WRITEPOT(1,SLABL,NAME,DECM,PV,PU,PS,CM,VMAXIN)
c
c** Now ... count potential parameters of various types for each state
c=======================================================================
c** Counters for numbers of potential parameters of different types for 
c   each state
c     COMMON /BLKCOUNT/TOTPOTPAR,POTPARI,POTPARF,UAPARI,UAPARF,
c    1  UBPARI,UBPARF,TAPARI,TAPARF,TBPARI,TBPARF,LDPARI,LDPARF
c=======================================================================
      TOTPOTPAR= 0
      NBCTOT= 0
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
                          NBCTOT= NBCTOT + NBC(I,ISOT,ISTATE)
                          BCPARF(I,ISOT,ISTATE)= IPV
                          ENDIF
                      IF(NQC(I,ISOT,ISTATE).GT.0) THEN
                          QCPARI(I,ISOT,ISTATE)= IPV+1
                          DO  J= 1,NQC(I,ISOT,ISTATE)
                              IPV= IPV+1
                              IFXPV(IPV)= 0
                              PV(IPV)= 0.d0
                              PU(IPV)= 0.d0
                              ENDDO
                          NBCTOT= NBCTOT + NQC(I,ISOT,ISTATE)
                          QCPARF(I,ISOT,ISTATE)= IPV
                          ENDIF
                      ENDDO
                  ENDDO
cc            TOTPOTPAR= IPV                  !! ?? save BC for GPROUND
              GOTO 90
              ENDIF
          POTPARI(ISTATE)= IPV+1
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
c... For all fitted potentials except  GPEF count De
          IF((PSEL(ISTATE).LT.4)) THEN
              IPV= IPV+ 1          
              IFXPV(IPV)= IFXDe(ISTATE)
              ENDIF
c... For all fitted potentials count Re
          IF((PSEL(ISTATE).LE.4)) THEN 
              IPV= IPV+1
              IF(PSEL(ISTATE).EQ.4) POTPARI(ISTATE)= IPV
              IFXPV(IPV)= IFXRe(ISTATE)
              ENDIF
          IF((PSEL(ISTATE).EQ.2).OR.(PSEL(ISTATE).EQ.3)) THEN
c... For MLR and  DELR, forms, count long-range parameters: count Cm's
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
c** For Pashov-exponent SE-MLR, or TT or HDF parameter count starts with 1
          IF((APSE(ISTATE).GT.0).OR.(PSEL(ISTATE).GE.6)) J=1
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
          TOTPOTPAR= IPV
c..... end of parameter count/label loop!
   90     CONTINUE
cc    IF(TOTPOTPAR.EQ.0) TOTPOTPAR= IPV     !! ?? spurious - unneeded ?
      IF(TOTPOTPAR.GT.HPARMX) THEN
          WRITE(6,626) TOTPOTPAR,HPARMX
          STOP
          ENDIF
      NPTOT= IPV
      NFPAR= 0
c** Count total free Hamiltonian fitting parameters
      DO  IPV= 1, TOTPOTPAR
          IF(IFXPV(IPV).LE.0) NFPAR= NFPAR+ 1
          ENDDO
c------------ Finished counting Hamiltonian Parameters------------------
  626 FORMAT(/' *** Dimension Error *** [(total No. Hamiltonian parmaete
     1rs)=',i4,'] >  HPARMX=',i4)
  638 FORMAT(' State ',A3,'  Energy Convergence criterion EPS is',
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
      CALL READATA(PASok,UCUTOFF,NDAT,NOWIDTHS,PRINP)
      NTVALL(0)= 0
      NTVSTOT= 0
      DO  ISTATE= 1,NSTATES
          IF(PSEL(ISTATE).EQ.-2) THEN
c... If this state to be represented by term values, determine the number
c    and add them to the parameter count
              NTVI(ISTATE)= NPTOT+ 1      !! note: TVSORT updates NPTOT
              CALL TVSORT(ISTATE,NPTOT,VMAX,NTVALL,NTVSSTAT,TVNAME)
              NTVALL(0)= NTVALL(0) + NTVALL(ISTATE)
              IF(NTVALL(ISTATE).GT.0) THEN
                  NTVF(ISTATE)= NPTOT
                  ENDIF
              IF(NTVSSTAT.GT.0) NTVSTOT= NTVSTOT+ NTVSSTAT
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
          EPS(ISTATE)= DMIN1(UCUTOFF/100.0d0,1.D-06)
cc        EPS(ISTATE)= MIN(UCUTOFF/10.0d0,1.d-06)
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
      IF((NFSTOT.GT.0).OR.(NTVALL(0).GT.0).OR.(NBCTOT.GT.0)) THEN
c** If HAVE fluorescence series and/or fitted term values and/or band 
c  constants ... first fix ALL potential parameters and fit to determine
c  estimates of the series origins and/or term values, and only THEN
c  free potential parameters too.  
c** Start by saving read-in values of  'IF(fix)' parameters
          DO  I= 1,TOTPOTPAR
              SIFXPV(I)= IFXPV(I)
              IFXPV(I)= 1
              ENDDO
          DO  I= TOTPOTPAR+1,NPTOT
              IFXPV(I)= 0
              PV(I)= 0.d0
cc            IF(IFXFS(I-TOTPOTPAR).GT.0) IFXPV(I)= 1  !!?????
cc            IF(IFXFS(I-TOTPOTPAR+NBCTOT).GT.0) IFXPV(I)= 1
              ENDDO
c** First, fit ONLY to Fluorescence series origins and/or free Term
c            Value and/or band constants
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
  888 format(/' Following FS fit, reset  T(',i3,')=',f12.4,
     1  '  equal   T(',I3,')=', F12.4)
                  TVALUE(I-TOTPOTPAR)= TVALUE(IFXFS(I-TOTPOTPAR))
                  IFXPV(I)= 1
                  ENDIF
              ENDDO
          NFPAR= NFPAR+ NFSTOT+ NTVALL(0)+ NBCTOT 
          CALL MAPPAR(NISTP,PV,0)
          ENDIF
c--- End of section to determine preliminary values of any fluorescence 
c    series origins, aterm Values or Band Constants
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
c** Perform group rounding of all band constants and/or term values, 
c   and/or fluorescence series origins in a single step
          IF((NFSTOT.GT.0).OR.(NTVALL(0).GT.0).OR.(NBCTOT.GT.0)) THEN
              IROUND= IABS(SIROUND) + 2
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
          ENDIF
c** Writing out the general information of the fit.
c-----------------------------------------------------------------------
      WRITE(6,691) NFPAR,COUNTOT,DSE
c-----------------------------------------------------------------------
c** Writing out the fluorescence band results.
c-----------------------------------------------------------------------
      IF(NFSTOT.GT.0) THEN
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
  691 FORMAT(/,1X,36('==')/' Fitting',I5,'  free parameters to',I6,
     1  '  transitions yields  DSE=',G15.8/1X,36('=='))
  692 FORMAT(/1X,33('==')/'  The following',I5,' Fluorescence Series Ori
     1gins were determined'/1x,30('--')/"  ( v', J', p'; ISTP)",4x,
     2 'T(value)',4x, 'Uncertainty  Sensitivity'/1x,30('--'))
cc694 FORMAT(3X,'(',I3,',',I3,',',SP,I3,SS,';',I2,')',1X,1PD19.10,
  694 FORMAT(2X,'(',I4,',',I3,',',SP,I3,SS,';',I2,')',1X,F15.6,
     1  1PD11.1,D12.1)
  696 FORMAT(/1X,33('==')/'  State ',A3,' represented by the',I5,' indiv
     1idual term values:'/1x,33('--')/" T(es: v', J', p';IS)  #dat",4x,
     2  'T(value)',4x,'Uncertainty  Sensitivity'/1x,33('--'))
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
              RMAXT= RD(NDATPT(ISTATE),ISTATE)
              RHT= RD(2,ISTATE) - RD(1,ISTATE) 
              nPointSR(ISTATE)= RD(1,ISTATE)/RHT 
              IF(nPointSR(ISTATE).GT.NPNTMX) nPointSR(ISTATE)= NPNTMX
              IF(OSEL(ISTATE).NE. 0) THEN
                  IF(RMAXT .GT. 100.0) THEN
                      nPointLR(ISTATE)= 0
                    ELSE
                      RHL= RHT*OSEL(ISTATE)
                      nPointLR(ISTATE)= (100.0-RMAXT)/RHL
                      IF(nPointLR(ISTATE).GT.NPNTMX) THEN
                          RHL= RHL*DFLOAT(nPointLR(ISTATE))/NPNTMX
                          nPointLR(ISTATE)= NPNTMX
                          ENDIF
                    ENDIF
                  ENDIF
              CALL VGEN(ISTATE,-1.0d0,VDIST,BETADIST,1)
              IB(NDATAMX)= NPARMX    ! omits centrifugal bits from  VGEN
              IF(OSEL(ISTATE).GT.0) THEN
                  J= MAX1(1.,OSEL(ISTATE)/10.)
                  DO  I= 1,nPointSR(ISTATE), J
c ... generate potential & exponent values in inner extrapolation region
                      RDIST= RHT*DBLE(I) 
                      Rsr(I,ISTATE)= RDIST
                      CALL VGEN(ISTATE,RDIST,VDIST,BETADIST,-1)
                      Vsr(I,ISTATE)= VDIST
                      Bsr(I,ISTATE)= BETADIST
                      ENDDO
                  DO  I= 1,nPointLR(ISTATE)
c ... generate potential & exponent values in outer extrapolation region
                      RDIST= RMAXT + RHL*DBLE(I)
                      Rlr(I,ISTATE)= RDIST
                      CALL VGEN(ISTATE,RDIST,VDIST,BETADIST,-1)
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
      CALL WRITEPOT(2,SLABL,NAME,DECM,PV,PU,PS,CM,VMAXIN)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** If chosen, output file(s) will be created for the export of the
c   generated functions: V, BETAFX, UAR/UBR, or TAR/TBR and their
c   respective uncertainties.
      DO  ISTATE= 1, NSTATES
          IF(OSEL(ISTATE).NE.0) THEN
              IF(PSEL(ISTATE).GT.0) THEN
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Call subroutine to print out the generated functions and their
c   respective uncertainties.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                  CALL FUNUNC(ISTATE,WRITFILE,PU,CM)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                ELSE 
                  DO  I= 1,NDATPT(NSTATES),IABS(OSEL(ISTATE))
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
      CALL DIFFSTATS(NSTATES,NFPAR,ROBUST,MKPRED,NPTOT,NTVSTOT,PRINP)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IF(uBv.GT.0) THEN
c** If desired, calculate Bv uncertainties
          CALL UNCBV(NPTOT,PV,PU,CM)
          ENDIF
      STOP
  900 FORMAT(5X,G18.8,5X,G18.8)
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE MASSES(IAN,IMN,NAME,GELGS,DGNS,MASS,ABUND)
c***********************************************************************
c** For isotope with (input) atomic number IAN and mass number IMN,
c  return (output):  (i) as the right-adjusted 2-character variable NAME
c  the alphabetic symbol for that element,  (ii) the ground state
c  electronic degeneracy GELGS, (iii) the nuclear spin degeneracy DGNS,
c  (iv) the atomic mass MASS [amu], and  (v) the natural isotopic
c  abundance ABUND [in percent].   GELGS values based on atomic states
c  in Moore's "Atomic Energy Level" tables, the isotope masses are taken
c  from the 2012 mass table [Wang, Audi, Wapstra, Kondev, MacCormick, Xu
c  & Pfeiffer, Chin.Phys.C 36, 1603-2014 (2012)] ,the proton, deuteron,
c  and triton masses are taken from the 2010 fundamental constants table 
c  [Mohr, Taylor, & Newell, Rev. Mod. Phys. 84, 1587-1591 (2012)] and other
c  quantities from Tables 6.2 and 6.3 of "Quantities, Units and Symbols in
c  Physical Chemistry", by Mills et al.(Blackwell,2'nd Edition, Oxford,1993).
c** If the input value of IMN does not equal one of the tabulated values
c  for atomic species IAN, return the abundance-averaged standard atomic
c  weight of that atom and set DGNS=-1 and ABUND=-1.
c** For Atomic number IAN=0 and isotope mass numbers IMN=1-3,  return the
c    masses of the proton, deuteron, and triton, p,d & t, respectively
c Masses and properties of selected Halo nuclei an unstable nuclei included
c                 COPYRIGHT 2005-2015  :  last  updated 10 January 2016
c** By R.J. Le Roy, with assistance from 
c                 G.T. Kraemer, J.Y. Seto and K.V. Slaughter.
c***********************************************************************
      REAL*8 zm(0:123,0:15),mass,ab(0:123,15),abund
      INTEGER i,ian,imn,gel(0:123),nmn(0:123),mn(0:123,15),
     1                                        gns(0:123,15),DGNS,gelgs
      CHARACTER*2 NAME,AT(0:123)
cc
      DATA  at(0),gel(0),nmn(0),(mn(0,i),i=1,3)/' p',1,3,1,2,3/
      DATA  (zm(0,i),i=0,3)/1.008d0,1.007276466812d0,2.013553212712d0,
     2                 3.0155007134d0/
      DATA  (gns(0,i),i=1,3)/2,3,2/
      DATA  (ab(0,i),i=1,3)/0.d0, 0.d0, 0.d0/
c
      DATA  at(1),gel(1),nmn(1),(mn(1,i),i=1,3)/' H',2,3,1,2,3/
      DATA  (zm(1,i),i=0,3)/1.00794d0, 1.00782503223d0, 2.01410177812d0,
     1                 3.0160492779d0/
      DATA  (gns(1,i),i=1,3)/2,3,2/
      DATA  (ab(1,i),i=1,3)/99.985d0,0.015d0,0.d0/
c
      DATA  at(2),gel(2),nmn(2),(mn(2,i),i=1,4)/'He',1,4,3,4,6,8/
      DATA  (zm(2,i),i=0,4)/4.002602d0, 3.0160293201d0, 4.00260325413d0,
     1                                        6.0188891d0, 8.033922d0/
      DATA  (gns(2,i),i=1,4)/2,1,1,1/
      DATA  (ab(2,i),i=1,4)/0.000137d0,99.999863d0, 2*0.d0/
c
      DATA  at(3),gel(3),nmn(3),(mn(3,i),i=1,6)/'Li',2,6,6,7,8,9,11,12/
      DATA  (zm(3,i),i=0,6)/6.941d0, 6.0151228874d0, 7.016003437d0,
     1     8.02248736d0,9.0267895d0,11.043798d0,12.05378d0/
      DATA  (gns(3,i),i=1,6)/3,4,5,4,4,1/
      DATA  (ab(3,i),i=1,6)/7.5d0, 92.5d0, 4*0.d0/
c
      DATA  at(4),gel(4),nmn(4),(mn(4,i),i=1,8)/'Be',1,8,7,9,10,11,12,
     1                                                       14,15,16/
      DATA  (zm(4,i),i=0,8)/9.012182d0, 7.01692983d0, 9.01218307d0,
     1 10.0135338d0, 11.021658d0, 12.026921d0, 14.04289d0, 15.05346d0,
     2 16.06192d0/
      DATA  (gns(4,i),i=1,8)/4,4,3,2,1,1,2,1/
      DATA  (ab(4,i),i=1,8)/0.d0, 100.d0, 6*0.d0/
c
      DATA at(5),gel(5),nmn(5),(mn(5,i),i=1,10)/' B',2,10,8,10,11,12,
     1                                              13,14,15,17,18,19/
      DATA (zm(5,i),i=0,10)/10.811d0, 8.0246072d0, 10.0129369d0, 
     1          11.0093054d0, 12.0143521d0, 13.0177802d0, 14.025404d0,
     2          15.031103d0, 17.04699d0, 18.05617d0,19.06373d0/
      DATA  (gns(5,i),i=1,10)/5,7,4,3,4,5,4,4,1,4/
      DATA  (ab(5,i),i=1,10)/0.d0, 19.9d0,80.1d0, 7*0.d0/
c
      DATA at(6),gel(6),nmn(6),(mn(6,i),i=1,14)/' C',1,14,9,10,11,12,13,
     1               14,15,16,17,18,19,20,21,22/
      DATA (zm(6,i),i=0,14)/12.011d0, 9.0310367d0, 10.0168532d0,
     1          11.0114336d0, 12.d0, 13.00335483507d0, 14.003241989d0, 
     1  15.0105993d0, 16.014701d0, 17.022586d0, 18.02676d0, 19.03481d0,
     2  20.04032d0, 21.04934d0, 22.05720d0/
      DATA  (gns(6,i),i=1,14)/4,1,4,1,2,1,2,1,4,1,2,1,2,1/
      DATA  (ab(6,i),i=1,14)/3*0.d0, 98.90d0,1.10d0, 9*0.d0/
c
      DATA at(7),gel(7),nmn(7),(mn(7,i),i=1,2)/' N',4,2,14,15/
      DATA (zm(7,i),i=0,2)/14.00674d0, 14.00307400443d0,15.0001088989d0/
      DATA (gns(7,i),i=1,2)/3,2/
      DATA (ab(7,i),i=1,2)/99.634d0,0.366d0/
c
      DATA at(8),gel(8),nmn(8),(mn(8,i),i=1,3)/' O',5,3,16,17,18/
      DATA (zm(8,i),i=0,3)/15.9994d0, 15.99491461957d0, 16.9991317565d0,
     1                      17.9991596129d0/
      DATA (gns(8,i),i=1,3)/1,6,1/
      DATA (ab(8,i),i=1,3)/99.762d0, 0.038d0, 0.200d0/
c
      DATA at(9),gel(9),nmn(9),(mn(9,i),i=1,1)/' F',4,1,19/
      DATA (zm(9,i),i=0,1)/18.9984032d0, 18.9984031627d0/
      DATA (gns(9,i),i=1,1)/2/
      DATA (ab(9,i),i=1,1)/100.d0/
c
      DATA at(10),gel(10),nmn(10),(mn(10,i),i=1,4)/'Ne',1,4,17,20,21,22/
      DATA (zm(10,i),i=0,4)/20.1797d0, 17.017672d0, 19.9924401762d0, 
     1                                   20.99384669d0,21.991385115d0/
      DATA (gns(10,i),i=1,4)/2,1,4,1/
      DATA (ab(10,i),i=1,4)/0.d0, 90.48d0, 0.27d0, 9.25d0/
c
      DATA at(11),gel(11),nmn(11),(mn(11,i),i=1,1)/'Na',2,1,23/
      DATA (zm(11,i),i=0,1)/22.989768d0, 22.9897692820d0/
      DATA (gns(11,i),i=1,1)/4/
      DATA (ab(11,i),i=1,1)/100.d0/
c
      DATA at(12),gel(12),nmn(12),(mn(12,i),i=1,3)/'Mg',1,3,24,25,26/
      DATA (zm(12,i),i=0,3)/24.3050d0, 23.985041698d0, 24.98583698d0,
     1                       25.98259297d0/
      DATA (gns(12,i),i=1,3)/1,6,1/
      DATA (ab(12,i),i=1,3)/78.99d0, 10.00d0, 11.01d0/
c
      DATA at(13),gel(13),nmn(13),(mn(13,i),i=1,1)/'Al',2,1,27/
      DATA (zm(13,i),i=0,1)/26.981539d0, 26.98153853d0/
      DATA (gns(13,i),i=1,1)/6/
      DATA (ab(13,i),i=1,1)/100.d0/
c
      DATA at(14),gel(14),nmn(14),(mn(14,i),i=1,3)/'Si',1,3,28,29,30/
      DATA (zm(14,i),i=0,3)/28.0855d0, 27.9769265346d0, 28.9764946649d0,
     1                       29.973770136d0/
      DATA (gns(14,i),i=1,3)/1,2,1/
      DATA (ab(14,i),i=1,3)/92.23d0, 4.67d0, 3.10d0/
 
      DATA at(15),gel(15),nmn(15),(mn(15,i),i=1,2)/' P',4,2,26,31/
      DATA (zm(15,i),i=0,2)/30.973762d0, 26.01178d0, 30.9737619984d0/
      DATA (gns(15,i),i=1,2)/15,2/
      DATA (ab(15,i),i=1,2)/0.d0, 100.d0/
c
      DATA at(16),gel(16),nmn(16),(mn(16,i),i=1,5)/' S',5,5,27,32,33,
     1                                                          34,36/
      DATA (zm(16,i),i=0,5)/32.066d0, 27.01883d0, 31.9720711744d0,
     1                   32.9714589098d0,33.96786700d0, 35.96708071d0/
      DATA (gns(16,i),i=1,5)/6,1,4,1,1/
      DATA (ab(16,i),i=1,5)/0.d0, 95.02d0, 0.75d0, 4.21d0, 0.02d0/
c
      DATA at(17),gel(17),nmn(17),(mn(17,i),i=1,2)/'Cl',4,2,35,37/
      DATA (zm(17,i),i=0,2)/35.4527d0, 34.96885268d0, 36.96590260d0/
      DATA (gns(17,i),i=1,2)/4,4/
      DATA (ab(17,i),i=1,2)/75.77d0, 24.23d0/
c
      DATA at(18),gel(18),nmn(18),(mn(18,i),i=1,3)/'Ar',1,3,36,38,40/
      DATA (zm(18,i),i=0,3)/39.948d0, 35.967545105d0, 37.96273211d0,
     1                       39.9623831237d0/
      DATA (gns(18,i),i=1,3)/1,1,1/
      DATA (ab(18,i),i=1,3)/0.337d0, 0.063d0, 99.600d0/
c
      DATA at(19),gel(19),nmn(19),(mn(19,i),i=1,3)/' K',2,3,39,40,41/
      DATA (zm(19,i),i=0,3)/39.0983d0, 38.963706486d0, 39.96399817d0,
     1                       40.961825258d0/
      DATA (gns(19,i),i=1,3)/4,9,4/
      DATA (ab(19,i),i=1,3)/93.2581d0, 0.0117d0, 6.7302d0/
 
      DATA at(20),gel(20),nmn(20),(mn(20,i),i=1,6)/'Ca',1,6,40,42,43,44,
     1                                              46,48/
      DATA (zm(20,i),i=0,6)/40.078d0, 39.962590864d0, 41.95861783d0,
     1         42.95876644d0, 43.9554816d0, 45.9536890d0, 47.95252277d0/
      DATA (gns(20,i),i=1,6)/1,1,8,1,1,1/
      DATA (ab(20,i),i=1,6)/96.941d0, 0.647d0, 0.135d0, 2.086d0,
     1                      0.004d0, 0.187d0/
c
      DATA at(21),gel(21),nmn(21),(mn(21,i),i=1,1)/'Sc',4,1,45/
      DATA (zm(21,i),i=0,1)/44.955910d0, 44.9559083d0/
      DATA (gns(21,i),i=1,1)/8/
      DATA (ab(21,i),i=1,1)/100.d0/
c
      DATA at(22),gel(22),nmn(22),(mn(22,i),i=1,5)/'Ti',5,5,46,47,48,49,
     1                                              50/
      DATA (zm(22,i),i=0,5)/47.88d0, 45.9526277d0, 46.9517588d0,
     1         47.9479420d0, 48.9478657d0, 49.9447869d0/
      DATA (gns(22,i),i=1,5)/1,6,1,8,1/
      DATA (ab(22,i),i=1,5)/8.0d0, 7.3d0, 73.8d0, 5.5d0, 5.4d0/
c
      DATA at(23),gel(23),nmn(23),(mn(23,i),i=1,2)/' V',4,2,50,51/
      DATA (zm(23,i),i=0,2)/50.9415d0, 49.9471560d0, 50.9439570d0/
      DATA (gns(23,i),i=1,2)/13,8/
      DATA (ab(23,i),i=1,2)/0.250d0, 99.750d0/
c
      DATA at(24),gel(24),nmn(24),(mn(24,i),i=1,4)/'Cr',7,4,50,52,53,54/
      DATA (zm(24,i),i=0,4)/51.9961d0, 49.9460418d0, 51.9405062d0,
     1                       52.9406481d0, 53.9388792d0/
      DATA (gns(24,i),i=1,4)/1,1,4,1/
      DATA (ab(24,i),i=1,4)/4.345d0, 83.789d0, 9.501d0, 2.365d0/
c
      DATA at(25),gel(25),nmn(25),(mn(25,i),i=1,1)/'Mn',6,1,55/
      DATA (zm(25,i),i=0,1)/54.93805d0, 54.938049d0/
      DATA (gns(25,i),i=1,1)/6/
      DATA (ab(25,i),i=1,1)/100.d0/
c
      DATA at(26),gel(26),nmn(26),(mn(26,i),i=1,4)/'Fe',9,4,54,56,57,58/
      DATA (zm(26,i),i=0,4)/55.847d0, 53.9396090d0, 55.9349363d0,
     1                       56.9353928d0, 57.9332744d0/
      DATA (gns(26,i),i=1,4)/1,1,2,1/
      DATA (ab(26,i),i=1,4)/5.8d0, 91.72d0, 2.2d0, 0.28d0/
c
      DATA at(27),gel(27),nmn(27),(mn(27,i),i=1,1)/'Co',10,1,59/
      DATA (zm(27,i),i=0,1)/58.93320d0, 58.9331943d0/
      DATA (gns(27,i),i=1,1)/8/
      DATA (ab(27,i),i=1,1)/100.d0/
c
      DATA at(28),gel(28),nmn(28),(mn(28,i),i=1,5)/'Ni',9,5,58,60,61,62,
     1                                              64/
      DATA (zm(28,i),i=0,5)/58.69d0, 57.9353424d0, 59.9307859d0,
     1         60.9310556d0, 61.9283454d0, 63.9279668d0/
      DATA (gns(28,i),i=1,5)/1,1,4,1,1/
      DATA (ab(28,i),i=1,5)/68.077d0,26.223d0,1.140d0,3.634d0,0.926d0/
c
      DATA at(29),gel(29),nmn(29),(mn(29,i),i=1,2)/'Cu',2,2,63,65/
      DATA (zm(29,i),i=0,2)/63.546d0, 62.9295977d0,64.9277897d0/
      DATA (gns(29,i),i=1,2)/4,4/
      DATA (ab(29,i),i=1,2)/69.17d0, 30.83d0/
c
      DATA at(30),gel(30),nmn(30),(mn(30,i),i=1,5)/'Zn',1,5,64,66,67,68,
     1                                              70/
      DATA (zm(30,i),i=0,5)/65.40d0, 63.9291420d0, 65.9260338d0,
     1         66.9271277d0, 67.9248446d0, 69.9253192d0/
      DATA (gns(30,i),i=1,5)/1,1,6,1,1/
      DATA (ab(30,i),i=1,5)/48.6d0, 27.9d0, 4.1d0, 18.8d0, 0.6d0/
c
      DATA at(31),gel(31),nmn(31),(mn(31,i),i=1,2)/'Ga',2,2,69,71/
      DATA (zm(31,i),i=0,2)/69.723d0, 68.9255735d0, 70.9247026d0/
      DATA (gns(31,i),i=1,2)/4,4/
      DATA (ab(31,i),i=1,2)/60.108d0, 39.892d0/
c
      DATA at(32),gel(32),nmn(32),(mn(32,i),i=1,5)/'Ge',1,5,70,72,73,74,
     1                                              76/
      DATA (zm(32,i),i=0,5)/72.61d0, 69.9242488d0, 71.92207583d0,
     1         72.92345896d0, 73.921177762d0, 75.921402726d0/
      DATA (gns(32,i),i=1,5)/1,1,10,1,1/
      DATA (ab(32,i),i=1,5)/21.23d0, 27.66d0, 7.73d0, 35.94d0, 7.44d0/
c
      DATA at(33),gel(33),nmn(33),(mn(33,i),i=1,1)/'As',4,1,75/
      DATA (zm(33,i),i=0,1)/74.92159d0, 74.9215946d0/
      DATA (gns(33,i),i=1,1)/4/
      DATA (ab(33,i),i=1,1)/100.d0/
c
      DATA at(34),gel(34),nmn(34),(mn(34,i),i=1,6)/'Se',5,6,74,76,77,78,
     1                                              80,82/
      DATA (zm(34,i),i=0,6)/78.96d0, 73.922475935d0, 75.919213704d0,
     1         76.91991415d0, 77.91730928d0, 79.9165218d0, 81.9166995d0/
      DATA (gns(34,i),i=1,6)/1,1,2,1,1,1/
      DATA (ab(34,i),i=1,6)/0.89d0, 9.36d0, 7.63d0, 23.78d0, 49.61d0,
     1                      8.73d0/
c
      DATA at(35),gel(35),nmn(35),(mn(35,i),i=1,2)/'Br',4,2,79,81/
      DATA (zm(35,i),i=0,2)/79.904d0, 78.9183376d0, 80.9162897d0/
      DATA (gns(35,i),i=1,2)/4,4/
      DATA (ab(35,i),i=1,2)/50.69d0, 49.31d0/
c
      DATA at(36),gel(36),nmn(36),(mn(36,i),i=1,6)/'Kr',1,6,78,80,82,83,
     1                                              84,86/
      DATA (zm(36,i),i=0,6)/83.80d0, 77.9203649d0, 79.9163781d0,
     1     81.9134827d0, 82.9141272d0, 83.911497728d0, 85.910610627d0/
      DATA (gns(36,i),i=1,6)/1,1,1,10,1,1/
      DATA (ab(36,i),i=1,6)/0.35d0, 2.25d0, 11.6d0, 11.5d0, 57.0d0,
     1                      17.3d0/
c
      DATA at(37),gel(37),nmn(37),(mn(37,i),i=1,2)/'Rb',2,2,85,87/
      DATA (zm(37,i),i=0,2)/85.4678d0, 84.911789738d0, 86.909180532d0/
      DATA (gns(37,i),i=1,2)/6,4/
      DATA (ab(37,i),i=1,2)/72.165d0, 27.835d0/
c
      DATA at(38),gel(38),nmn(38),(mn(38,i),i=1,4)/'Sr',1,4,84,86,87,88/
      DATA (zm(38,i),i=0,4)/87.62d0, 83.9134191d0, 85.9092606d0,
     1                      86.9088775d0, 87.9056125d0/
      DATA (gns(38,i),i=1,4)/1,1,10,1/
      DATA (ab(38,i),i=1,4)/0.56d0, 9.86d0, 7.00d0, 82.58d0/
c
      DATA at(39),gel(39),nmn(39),(mn(39,i),i=1,1)/' Y',4,1,89/
      DATA (zm(39,i),i=0,1)/88.90585d0, 88.9058403d0/
      DATA (gns(39,i),i=1,1)/2/
      DATA (ab(39,i),i=1,1)/100.d0/
c
      DATA at(40),gel(40),nmn(40),(mn(40,i),i=1,5)/'Zr',5,5,90,91,92,94,
     1                                              96/
      DATA (zm(40,i),i=0,5)/91.224d0, 89.9046977d0, 90.9056396d0,
     1                      91.9050347d0, 93.9063108d0, 95.9082714d0/
      DATA (gns(40,i),i=1,5)/1,6,1,1,1/
      DATA (ab(40,i),i=1,5)/51.45d0, 11.22d0, 17.15d0, 17.38d0, 2.80d0/
c
      DATA at(41),gel(41),nmn(41),(mn(41,i),i=1,1)/'Nb',2,1,93/
      DATA (zm(41,i),i=0,1)/92.90638d0, 92.9063730d0/
      DATA (gns(41,i),i=1,1)/10/
      DATA (ab(41,i),i=1,1)/100.d0/
c
      DATA at(42),gel(42),nmn(42),(mn(42,i),i=1,7)/'Mo',7,7,92,94,95,96,
     1                                              97,98,100/
      DATA (zm(42,i),i=0,7)/95.94d0, 91.9068080d0, 93.9050849d0,
     1        94.9058388d0, 95.9046761d0, 96.9060181d0, 97.9054048d0,
     2        99.9074718d0/
      DATA (gns(42,i),i=1,7)/1,1,6,1,6,1,1/
      DATA (ab(42,i),i=1,7)/14.84d0, 9.25d0, 15.92d0, 16.68d0, 9.55d0,
     1                      24.13d0, 9.63d0/
c
      DATA at(43),gel(43),nmn(43),(mn(43,i),i=1,1)/'Tc',6,1,98/
      DATA (zm(43,i),i=0,1)/97.907215d0, 97.907212d0/
      DATA (gns(43,i),i=1,1)/13/
      DATA (ab(43,i),i=1,1)/100.d0/
c
      DATA at(44),gel(44),nmn(44),(mn(44,i),i=1,7)/'Ru',11,7,96,98,99,
     1                                              100,101,102,104/
      DATA (zm(44,i),i=0,7)/101.07d0, 95.9075903d0, 97.905287d0,
     1     98.9059341d0, 99.9042143d0, 100.9055769d0, 101.9043441d0,
     2     103.9054275d0/
      DATA (gns(44,i),i=1,7)/1,1,6,1,6,1,1/
      DATA (ab(44,i),i=1,7)/5.52d0, 1.88d0, 12.7d0, 12.6d0, 17.0d0,
     1                      31.6d0, 18.7d0/
c
      DATA at(45),gel(45),nmn(45),(mn(45,i),i=1,1)/'Rh',10,1,103/
      DATA (zm(45,i),i=0,1)/102.90550d0, 102.9054980d0/
      DATA (gns(45,i),i=1,1)/2/
      DATA (ab(45,i),i=1,1)/100.d0/
c
      DATA at(46),gel(46),nmn(46),(mn(46,i),i=1,6)/'Pd',1,6,102,104,105,
     1                                              106,108,110/
      DATA (zm(46,i),i=0,6)/106.42d0, 101.9056022d0, 103.9040305d0,
     1       104.9050796d0, 105.9034804d0, 107.9038916d0, 109.9051722d0/
      DATA (gns(46,i),i=1,6)/1,1,6,1,1,1/
      DATA (ab(46,i),i=1,6)/1.02d0, 11.14d0, 22.33d0, 27.33d0, 26.46d0,
     1                      11.72d0/
c
      DATA at(47),gel(47),nmn(47),(mn(47,i),i=1,2)/'Ag',2,2,107,109/
      DATA (zm(47,i),i=0,2)/107.8682d0, 106.9050916d0, 108.9047553d0/
      DATA (gns(47,i),i=1,2)/2,2/
      DATA (ab(47,i),i=1,2)/51.839d0, 48.161d0/
c
      DATA at(48),gel(48),nmn(48),(mn(48,i),i=1,8)/'Cd',1,8,106,108,110,
     1                                             111,112,113,114,116/ 
      DATA (zm(48,i),i=0,8)/112.411d0, 105.9064599d0, 107.9041834d0, 
     1       109.9030066d0, 110.9041829d0, 111.9027629d0, 112.9044081d0,
     2       113.9033651d0, 115.90476315d0/
      DATA (gns(48,i),i=1,8)/1,1,1,2,1,2,1,1/
      DATA (ab(48,i),i=1,8)/1.25d0, 0.89d0, 12.49d0, 12.80d0, 24.13d0,
     1                      12.22d0, 28.73d0, 7.49d0/
c
      DATA at(49),gel(49),nmn(49),(mn(49,i),i=1,2)/'In',2,2,113,115/
      DATA (zm(49,i),i=0,2)/114.818d0, 112.9040618d0, 114.903878776d0/
      DATA  (gns(49,i),i=1,2)/10,10/
      DATA (ab(49,i),i=1,2)/4.3d0, 95.7d0/
c
      DATA at(50),gel(50),nmn(50),(mn(50,i),i=1,10)/'Sn',1,10,112,114,
     1                                 115,116,117,118,119,120,122,124/
      DATA (zm(50,i),i=0,10)/118.710d0, 111.9048239d0, 113.9027827d0,
     1    114.903344699d0, 115.90174280d0, 116.9029540d0, 117.9016066d0,
     2    118.9033112d0, 119.9022016d0, 121.9034438d0, 123.9052766d0/
      DATA (gns(50,i),i=1,10)/1,1,2,1,2,1,2,1,1,1/
      DATA (ab(50,i),i=1,10)/0.97d0, 0.65d0, 0.34d0, 14.53d0, 7.68d0,
     1                       24.23d0, 8.59d0, 32.59d0, 4.63d0, 5.79d0/
c
      DATA at(51),gel(51),nmn(51),(mn(51,i),i=1,2)/'Sb',4,2,121,123/
      DATA (zm(51,i),i=0,2)/121.757d0, 120.903812d0, 122.9042132d0/
      DATA (gns(51,i),i=1,2)/6,8/
      DATA (ab(51,i),i=1,2)/57.36d0, 42.64d0/
c
      DATA at(52),gel(52),nmn(52),(mn(52,i),i=1,8)/'Te',5,8,120,122,123,
     1                                             124,125,126,128,130/
      DATA (zm(52,i),i=0,8)/127.60d0, 119.904059d0, 121.9030435d0,
     1    122.9042698d0, 123.9028171d0, 124.9044299d0, 125.9033109d0,
     2    127.9044613d0, 129.906222749d0/
      DATA (gns(52,i),i=1,8)/1,1,2,1,2,1,1,1/
      DATA (ab(52,i),i=1,8)/0.096d0, 2.603d0, 0.908d0, 4.816d0,
     1                      7.139d0, 18.95d0, 31.69d0, 33.80d0/
c
      DATA at(53),gel(53),nmn(53),(mn(53,i),i=1,2)/' I',4,2,127,129/
      DATA (zm(53,i),i=0,2)/126.90447d0, 126.904472d0, 128.904984d0/
      DATA (gns(53,i),i=1,2)/6,8/
      DATA (ab(53,i),i=1,2)/100.d0,0.d0/
c
      DATA at(54),gel(54),nmn(54),(mn(54,i),i=1,9)/'Xe',1,9,124,126,128,
     1                                          129,130,131,132,134,136/
      DATA (zm(54,i),i=0,9)/131.29d0, 123.9058920d0, 125.904298d0,
     1    127.9035310d0, 128.904780861d0,129.903509350d0,130.90508406d0,
     2    131.904155086d0, 133.9053947d0, 135.907214484d0/
      DATA (gns(54,i),i=1,9)/1,1,1,2,1,4,1,1,1/
      DATA (ab(54,i),i=1,9)/0.10d0, 0.09d0, 1.91d0, 26.4d0, 4.1d0,
     1                      21.2d0, 26.9d0, 10.4d0, 8.9d0/
c
      DATA at(55),gel(55),nmn(55),(mn(55,i),i=1,1)/'Cs',2,1,133/
      DATA (zm(55,i),i=0,1)/132.90543d0, 132.905451961d0/
      DATA (gns(55,i),i=1,1)/8/
      DATA (ab(55,i),i=1,1)/100.d0/
c
      DATA at(56),gel(56),nmn(56),(mn(56,i),i=1,7)/'Ba',1,7,130,132,134,
     1                                             135,136,137,138/
      DATA (zm(56,i),i=0,7)/137.327d0, 129.9063207d0, 131.9050611d0,
     1    133.90450818d0, 134.90568838d0, 135.90457573d0, 136.9058271d0,
     2    137.9052470d0/
      DATA (gns(56,i),i=1,7)/1,1,1,4,1,4,1/
      DATA (ab(56,i),i=1,7)/0.106d0, 0.101d0, 2.417d0, 6.592d0, 
     1                      7.854d0, 11.23d0, 71.70d0/
c
      DATA at(57),gel(57),nmn(57),(mn(57,i),i=1,2)/'La',4,2,138,139/
      DATA (zm(57,i),i=0,2)/138.9055d0, 137.907115d0, 138.9063563d0/
      DATA (gns(57,i),i=1,2)/11,8/ 
      DATA (ab(57,i),i=1,2)/0.0902d0, 99.9098d0/
c
      DATA at(58),gel(58),nmn(58),(mn(58,i),i=1,4)/'Ce',9,4,136,138,140,
     1                                             142/
      DATA (zm(58,i),i=0,4)/140.115d0, 135.9071292d0, 137.905991d0,
     1    139.9054431d0, 141.9092504d0/
      DATA (gns(58,i),i=1,4)/1,1,1,1/
      DATA (ab(58,i),i=1,4)/0.19d0, 0.25d0, 88.48d0, 11.08d0/
c
      DATA at(59),gel(59),nmn(59),(mn(59,i),i=1,1)/'Pr',10,1,141/
      DATA (zm(59,i),i=0,1)/140.90765d0, 140.9076576d0/
      DATA (gns(59,i),i=1,1)/6/
      DATA (ab(59,i),i=1,1)/100.d0/
c
      DATA at(60),gel(60),nmn(60),(mn(60,i),i=1,7)/'Nd',9,7,142,143,144,
     1                                             145,146,148,150/
      DATA (zm(60,i),i=0,7)/144.24d0, 141.9077290d0, 142.9098200d0,
     1    143.9100930d0, 144.9125793d0, 145.9131226d0, 147.9168993d0,
     2    149.9209022d0/
      DATA (gns(60,i),i=1,7)/1,8,1,8,1,1,1/
      DATA (ab(60,i),i=1,7)/27.13d0, 12.18d0, 23.80d0, 8.30d0, 17.19d0,
     1                       5.76d0, 5.64d0/
c
      DATA at(61),gel(61),nmn(61),(mn(61,i),i=1,1)/'Pm',6,1,145/
      DATA (zm(61,i),i=0,1)/144.912743d0, 144.912756d0/
      DATA (gns(61,i),i=1,1)/6/
      DATA (ab(61,i),i=1,1)/100.d0/
c
      DATA at(62),gel(62),nmn(62),(mn(62,i),i=1,7)/'Sm',1,7,144,147,148,
     1                                             149,150,152,154/
      DATA (zm(62,i),i=0,7)/150.36d0, 143.9120065d0, 146.9149044d0,
     1    147.9148292d0, 148.9171921d0, 149.9172829d0, 151.9197397d0,
     2    153.9222169d0/
      DATA (gns(62,i),i=1,7)/1,8,1,8,1,1,1/
      DATA (ab(62,i),i=1,7)/3.1d0, 15.0d0, 11.3d0, 13.8d0, 7.4d0,
     1                      26.7d0, 22.7d0/
c
      DATA at(63),gel(63),nmn(63),(mn(63,i),i=1,2)/'Eu',8,2,151,153/
      DATA (zm(63,i),i=0,2)/151.965d0, 150.9198578d0, 152.9212380d0/
      DATA (gns(63,i),i=1,2)/6,6/
      DATA (ab(63,i),i=1,2)/47.8d0, 52.2d0/
c
      DATA at(64),gel(64),nmn(64),(mn(64,i),i=1,7)/'Gd',5,7,152,154,155,
     1                                              156,157,158,160/
      DATA (zm(64,i),i=0,7)/157.25d0, 151.9197995d0, 153.9208741d0,
     1    154.9226305d0, 155.9221312d0, 156.9239686d0, 157.9241123d0,
     2    159.9270624d0/
      DATA (gns(64,i),i=1,7)/1,1,4,1,4,1,1/
      DATA (ab(64,i),i=1,7)/0.20d0, 2.18d0, 14.80d0, 20.47d0, 15.65d0,
     1                      24.84d0, 21.86d0/
c
      DATA at(65),gel(65),nmn(65),(mn(65,i),i=1,1)/'Tb',16,1,159/
      DATA (zm(65,i),i=0,1)/158.92534d0, 158.9253547d0/
      DATA (gns(65,i),i=1,1)/4/
      DATA (ab(65,i),i=1,1)/100.d0/
c
      DATA at(66),gel(66),nmn(66),(mn(66,i),i=1,7)/'Dy',17,7,156,158,
     1                                           160,161,162,163,164/
      DATA (zm(66,i),i=0,7)/162.50d0, 155.9242847d0, 157.924416d0,
     1    159.9252046d0, 160.9269405d0, 161.9268056d0, 162.9287383d0,
     2    163.9291819d0/
      DATA (gns(66,i),i=1,7)/1,1,1,6,1,6,1/
      DATA (ab(66,i),i=1,7)/0.06d0, 0.10d0, 2.34d0, 18.9d0, 25.5d0,
     1                      24.9d0, 28.2d0/
c
      DATA at(67),gel(67),nmn(67),(mn(67,i),i=1,1)/'Ho',16,1,165/
      DATA (zm(67,i),i=0,1)/164.93032d0, 164.9303288d0/
      DATA (gns(67,i),i=1,1)/8/
      DATA (ab(67,i),i=1,1)/100.d0/
     
      DATA at(68),gel(68),nmn(68),(mn(68,i),i=1,6)/'Er',13,6,162,164,
     1                                            166,167,168,170/
      DATA (zm(68,i),i=0,6)/167.26d0, 161.9287884d0, 163.9292088d0,
     1    165.9302995d0, 166.9320546d0, 167.9323767d0, 169.9354702d0/
      DATA (gns(68,i),i=1,6)/1,1,1,8,1,1/
      DATA (ab(68,i),i=1,6)/0.14d0, 1.61d0, 33.6d0, 22.95d0, 26.8d0,
     1                      14.9d0/
c
      DATA at(69),gel(69),nmn(69),(mn(69,i),i=1,1)/'Tm',8,1,169/  
      DATA (zm(69,i),i=0,1)/168.93421d0, 168.9342179d0/
      DATA (gns(69,i),i=1,1)/2/
      DATA (ab(69,i),i=1,1)/100.d0/
c
      DATA at(70),gel(70),nmn(70),(mn(70,i),i=1,7)/'Yb',1,7,168,170,171,
     1                                            172,173,174,176/
      DATA (zm(70,i),i=0,7)/173.04d0, 167.9338896d0, 169.9347664d0,
     1    170.9363302d0, 171.9363859d0, 172.9382151d0, 173.9388664d0,
     2    175.9425764d0/
      DATA (gns(70,i),i=1,7)/1,1,2,1,6,1,1/
      DATA (ab(70,i),i=1,7)/0.13d0, 3.05d0, 14.3d0, 21.9d0, 16.12d0,
     1                      31.8d0, 12.7d0/
c
      DATA at(71),gel(71),nmn(71),(mn(71,i),i=1,2)/'Lu',4,2,175,176/
      DATA (zm(71,i),i=0,2)/174.967d0, 174.9407752d0, 175.9426897d0/
      DATA (gns(71,i),i=1,2)/6,15/
      DATA (ab(71,i),i=1,2)/97.41d0, 2.59d0/
c
      DATA at(72),gel(72),nmn(72),(mn(72,i),i=1,6)/'Hf',5,6,174,176,177,
     1                                             178,179,180/
      DATA (zm(72,i),i=0,6)/178.49d0, 173.9400461d0, 175.9414076d0,
     1    176.9432277d0, 177.9437058d0, 178.9458232d0, 179.9465570d0/
      DATA (gns(72,i),i=1,6)/1,1,8,1,10,1/
      DATA (ab(72,i),i=1,6)/0.162d0, 5.206d0, 18.606d0, 27.297d0,
     1                      13.629d0, 35.100d0/
c
      DATA at(73),gel(73),nmn(73),(mn(73,i),i=1,2)/'Ta',4,2,180,181/
      DATA (zm(73,i),i=0,2)/180.9479d0, 179.9474648d0, 180.9479958d0/
      DATA (gns(73,i),i=1,2)/17,8/
      DATA (ab(73,i),i=1,2)/0.012d0, 99.988d0/
c
      DATA at(74),gel(74),nmn(74),(mn(74,i),i=1,5)/' W',1,5,180,182,183,
     1                                             184,186/
      DATA (zm(74,i),i=0,5)/183.84d0, 179.9467108d0, 181.9482039d0,
     1    182.9502227d0, 183.9509309d0, 185.9543628d0/
      DATA (gns(74,i),i=1,5)/1,1,2,1,1/
      DATA (ab(74,i),i=1,5)/0.13d0, 26.3d0, 14.3d0, 30.67d0, 28.6d0/
c
      DATA at(75),gel(75),nmn(75),(mn(75,i),i=1,2)/'Re',6,2,185,187/
      DATA (zm(75,i),i=0,2)/186.207d0, 184.9529545d0, 186.9557501d0/
      DATA (gns(75,i),i=1,2)/6,6/
      DATA (ab(75,i),i=1,2)/37.40d0, 62.60d0/
c
      DATA at(76),gel(76),nmn(76),(mn(76,i),i=1,7)/'Os',9,7,184,186,187,
     1                                             188,189,190,192/
      DATA (zm(76,i),i=0,7)/190.23d0, 183.9524885d0, 185.9538350d0,
     1    186.9557474d0, 187.9558352d0, 188.9581442d0, 189.9584437d0,
     2    191.9614770d0/
      DATA (gns(76,i),i=1,7)/1,1,2,1,4,1,1/
      DATA (ab(76,i),i=1,7)/0.02d0, 1.58d0, 1.6d0, 13.3d0, 16.1d0,
     1                      26.4d0, 41.0d0/
c
      DATA at(77),gel(77),nmn(77),(mn(77,i),i=1,2)/'Ir',10,2,191,193/
      DATA (zm(77,i),i=0,2)/192.22d0, 190.9605893d0, 192.9629216d0/
      DATA (gns(77,i),i=1,2)/4,4/
      DATA (ab(77,i),i=1,2)/37.3d0, 62.7d0/
c
c
      DATA at(78),gel(78),nmn(78),(mn(78,i),i=1,6)/'Pt',7,6,190,192,194,
     1                                            195,196,198/
      DATA (zm(78,i),i=0,6)/195.08d0, 189.959930d0, 191.961039d0,
     1    193.9626809d0, 194.9647917d0, 195.9649521d0, 197.9678949d0/
      DATA (gns(78,i),i=1,6)/1,1,1,2,1,1/
      DATA (ab(78,i),i=1,6)/0.01d0,0.79d0,32.9d0,33.8d0,25.3d0,7.2d0/
c
      DATA at(79),gel(79),nmn(79),(mn(79,i),i=1,1)/'Au',2,1,197/
      DATA (zm(79,i),i=0,1)/196.96654d0, 196.9665688d0/
      DATA (gns(79,i),i=1,1)/4/
      DATA (ab(79,i),i=1,1)/100.d0/
c
      DATA at(80),gel(80),nmn(80),(mn(80,i),i=1,7)/'Hg',1,7,196,198,199,
     1                                            200,201,202,204/
      DATA (zm(80,i),i=0,7)/200.59d0, 195.965833d0, 197.9667686d0,
     1    198.9682806d0, 199.9683266d0, 200.9703028d0, 201.9706434d0,
     2    203.9734940d0/
      DATA (gns(80,i),i=1,7)/1,1,2,1,4,1,1/
      DATA (ab(80,i),i=1,7)/0.15d0, 9.97d0, 16.87d0, 23.10d0, 13.18d0,
     1                      29.86d0, 6.87d0/
c
      DATA at(81),gel(81),nmn(81),(mn(81,i),i=1,2)/'Tl',2,2,203,205/
      DATA (zm(81,i),i=0,2)/204.3833d0, 202.9723446d0, 204.9744278d0/
      DATA (gns(81,i),i=1,2)/2,2/
      DATA (ab(81,i),i=1,2)/29.524d0, 70.476d0/
c
      DATA at(82),gel(82),nmn(82),(mn(82,i),i=1,4)/'Pb',1,4,204,206,207,
     1                                             208/
      DATA (zm(82,i),i=0,4)/207.2d0, 203.9730440d0, 205.9744657d0,
     1    206.9758973d0, 207.9766525d0/
      DATA (gns(82,i),i=1,4)/1,1,2,1/
      DATA (ab(82,i),i=1,4)/1.4d0, 24.1d0, 22.1d0, 52.4d0/
c
      DATA at(83),gel(83),nmn(83),(mn(83,i),i=1,1)/'Bi',4,1,209/
      DATA (zm(83,i),i=0,1)/208.98037d0, 208.9803991d0/
      DATA (gns(83,i),i=1,1)/10/
      DATA (ab(83,i),i=1,1)/100.d0/
c
      DATA at(84),gel(84),nmn(84),(mn(84,i),i=1,1)/'Po',5,1,209/
      DATA (zm(84,i),i=0,1)/208.982404d0, 208.9824308d0/
      DATA (gns(84,i),i=1,1)/2/
      DATA (ab(84,i),i=1,1)/100.d0/
c
      DATA at(85),gel(85),nmn(85),(mn(85,i),i=1,1)/'At',-1,1,210/
      DATA (zm(85,i),i=0,1)/209.987126d0, 209.987148d0/
      DATA (gns(85,i),i=1,1)/11/
      DATA (ab(85,i),i=1,1)/100.d0/
c
      DATA at(86),gel(86),nmn(86),(mn(86,i),i=1,1)/'Rn',1,1,222/
      DATA (zm(86,i),i=0,1)/222.017571d0, 222.0175782d0/
      DATA (gns(86,i),i=1,1)/1/
      DATA (ab(86,i),i=1,1)/100.d0/
c
      DATA at(87),gel(87),nmn(87),(mn(87,i),i=1,1)/'Fr',-1,1,223/
      DATA (zm(87,i),i=0,1)/223.019733d0, 223.0197360d0/
      DATA (gns(87,i),i=1,1)/4/
      DATA (ab(87,i),i=1,1)/100.d0/
c
      DATA at(88),gel(88),nmn(88),(mn(88,i),i=1,1)/'Ra',1,1,226/
      DATA (zm(88,i),i=0,1)/226.025403d0, 226.0254103d0/
      DATA (gns(88,i),i=1,1)/1/
      DATA (ab(88,i),i=1,1)/100.d0/
c
      DATA at(89),gel(89),nmn(89),(mn(89,i),i=1,1)/'Ac',4,1,227/
      DATA (zm(89,i),i=0,1)/227.027750d0, 227.0277523d0/
      DATA (gns(89,i),i=1,1)/4/
      DATA (ab(89,i),i=1,1)/100.d0/
c
      DATA at(90),gel(90),nmn(90),(mn(90,i),i=1,1)/'Th',-1,1,232/
      DATA (zm(90,i),i=0,1)/232.038d0, 232.0380558d0/
      DATA (gns(90,i),i=1,1)/1/
      DATA (ab(90,i),i=1,1)/100.d0/
c
      DATA at(91),gel(91),nmn(91),(mn(91,i),i=1,1)/'Pa',-1,1,231/
      DATA (zm(91,i),i=0,1)/231.03588d0, 231.0358842d0/
      DATA (gns(91,i),i=1,1)/4/
      DATA (ab(91,i),i=1,1)/100.d0/
c
      DATA at(92),gel(92),nmn(92),(mn(92,i),i=1,4)/' U',-1,4,233,234,
     1                                             235,238/
      DATA (zm(92,i),i=0,4)/238.0289d0, 233.0396355d0, 234.0409523d0,
     1    235.0439301d0, 238.0507884d0/
      DATA (gns(92,i),i=1,4)/6,1,8,1/
      DATA (ab(92,i),i=1,4)/0.d0, 0.0055d0, 0.7200d0, 99.2745d0/
c
      DATA at(93),gel(93),nmn(93),(mn(93,i),i=1,1)/'Np',-1,1,237/
      DATA (zm(93,i),i=0,1)/237.0481678d0, 237.0481736d0/
      DATA (gns(93,i),i=1,1)/6/
      DATA (ab(93,i),i=1,1)/100.d0/
c
      DATA at(94),gel(94),nmn(94),(mn(94,i),i=1,1)/'Pu',-1,1,244/
      DATA (zm(94,i),i=0,1)/244.064199d0, 244.064205d0/
      DATA (gns(94,i),i=1,1)/1/
      DATA (ab(94,i),i=1,1)/100.d0/
c
      DATA at(95),gel(95),nmn(95),(mn(95,i),i=1,1)/'Am',-1,1,243/
      DATA (zm(95,i),i=0,1)/243.061375d0, 243.0613815d0/
      DATA (gns(95,i),i=1,1)/6/
      DATA (ab(95,i),i=1,1)/100.d0/
c
      DATA at(96),gel(96),nmn(96),(mn(96,i),i=1,1)/'Cm',-1,1,247/
      DATA (zm(96,i),i=0,1)/247.070347d0, 247.070354d0/
      DATA (gns(96,i),i=1,1)/10/
      DATA (ab(96,i),i=1,1)/100.d0/
c
      DATA at(97),gel(97),nmn(97),(mn(97,i),i=1,1)/'Bk',-1,1,247/
      DATA (zm(97,i),i=0,1)/247.070300d0, 247.070307d0/
      DATA (gns(97,i),i=1,1)/4/
      DATA (ab(97,i),i=1,1)/100.d0/
c
      DATA at(98),gel(98),nmn(98),(mn(98,i),i=1,1)/'Cf',-1,1,251/
      DATA (zm(98,i),i=0,1)/251.079580d0, 251.079589d0/
      DATA (gns(98,i),i=1,1)/2/
      DATA (ab(98,i),i=1,1)/100.d0/
c
      DATA at(99),gel(99),nmn(99),(mn(99,i),i=1,1)/'Es',-1,1,252/
      DATA (zm(99,i),i=0,1)/252.082944d0, 252.082980d0/
      DATA (gns(99,i),i=1,1)/11/
      DATA (ab(99,i),i=1,1)/100.d0/
c
      DATA at(100),gel(100),nmn(100),(mn(100,i),i=1,1)/'Fm',-1,1,257/
      DATA (zm(100,i),i=0,1)/257.095099d0, 257.095106d0/
      DATA (gns(100,i),i=1,1)/10/
      DATA (ab(100,i),i=1,1)/100.d0/
c
      DATA at(101),gel(101),nmn(101),(mn(101,i),i=1,1)/'Md',-1,1,258/
      DATA (zm(101,i),i=0,1)/258.09857d0, 258.098431d0/
      DATA (gns(101,i),i=1,1)/17/
      DATA (ab(101,i),i=1,1)/100.d0/
c
      DATA at(102),gel(102),nmn(102),(mn(102,i),i=1,1)/'No',-1,1,259/
      DATA (zm(102,i),i=0,1)/259.100931d0, 259.101030d0/
      DATA (gns(102,i),i=1,1)/10/
      DATA (ab(102,i),i=1,1)/100.d0/
c
      DATA at(103),gel(103),nmn(103),(mn(103,i),i=1,1)/'Lr',-1,1,260/
      DATA (zm(103,i),i=0,1)/260.105320d0, 260.105510d0/
      DATA (gns(103,i),i=1,1)/-1/
      DATA (ab(103,i),i=1,1)/100.d0/
c
      DATA at(104),gel(104),nmn(104),(mn(104,i),i=1,1)/'Rf',-1,1,261/
      DATA (zm(104,i),i=0,1)/261.10869d0, 261.108770d0/
      DATA (gns(104,i),i=1,1)/-1/
      DATA (ab(104,i),i=1,1)/100.d0/
c
      DATA at(105),gel(105),nmn(105),(mn(105,i),i=1,1)/'Db',-1,1,262/
      DATA (zm(105,i),i=0,1)/262.11376d0, 262.114070d0/
      DATA (gns(105,i),i=1,1)/-1/
      DATA (ab(105,i),i=1,1)/100.d0/
c
      DATA at(106),gel(106),nmn(106),(mn(106,i),i=1,1)/'Sg',-1,1,263/
      DATA (zm(106,i),i=0,1)/263.11822d0, 263.118290d0/
      DATA (gns(106,i),i=1,1)/-1/
      DATA (ab(106,i),i=1,1)/100.d0/
c
      DATA at(107),gel(107),nmn(107),(mn(107,i),i=1,1)/'Bh',-1,1,262/
      DATA (zm(107,i),i=0,1)/262.12293d0, 262.122970d0/
      DATA (gns(107,i),i=1,1)/-1/
      DATA (ab(107,i),i=1,1)/100.d0/
c
      DATA at(108),gel(108),nmn(108),(mn(108,i),i=1,1)/'Hs',-1,1,265/
      DATA (zm(108,i),i=0,1)/265.13016d0, 265.129793d0/
      DATA (gns(108,i),i=1,1)/-1/
      DATA (ab(108,i),i=1,1)/100.d0/
c
      DATA at(109),gel(109),nmn(109),(mn(109,i),i=1,1)/'Mt',-1,1,266/
      DATA (zm(109,i),i=0,1)/266.13764d0, 266.137370d0/
      DATA (gns(109,i),i=1,1)/-1/
      DATA (ab(109,i),i=1,1)/100.d0/
c
      IF((IAN.LT.0).OR.(IAN.GT.109)) THEN
          MASS= 0.d0
          NAME= 'XX'
          IMN= 0
          WRITE(6,601) IAN
          RETURN
        ELSE
          NAME= AT(IAN)
        ENDIF
      IF((IAN.EQ.1).AND.(IMN.GT.1)) THEN
c** Special case: insert common name for deuterium or tritium
          IF(IMN.EQ.2) NAME=' D'
          IF(IMN.EQ.3) NAME=' T'
          ENDIF
      IF((IAN.EQ.0).AND.(IMN.GT.1)) THEN
          IF(IMN.EQ.2) NAME=' d'
          IF(IMN.EQ.3) NAME=' t'
          ENDIF
      GELGS= GEL(IAN)
      MASS= -1.d0
      DGNS= -1
      ABUND = -1.d0
      DO  I= 1,NMN(IAN)
          if(i.gt.15)  write(6,606) ian,imn,nmn(ian)
          IF(IMN.EQ.MN(IAN,I)) THEN
              MASS= ZM(IAN,I)
              DGNS= gns(IAN,I)
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
  606  format(/' *** ERROR *** called MASSES for atom with  AN=',I4,
     1  '  MN=',I4,'n(MN)=',I4)
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE READATA(PASok,UCUTOFF,NDAT,NOWIDTHS,PRINP)
c***********************************************************************
c** Subroutine to read, do book-keeping for, and print summary of
c  experimental data used in fits to spectroscopic data for one or more
c  electronic states and one or more isotopologues. 
c             ********* Version of 11 July 2015 *********
c             last change ... add acoustic virial data type
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++  COPYRIGHT 1997-2015 by  Robert J. Le Roy & Dominique R.T. Appadoo +
c   Dept. of Chemistry, Univ. of Waterloo, Waterloo, Ontario, Canada   +
c    This software may not be sold or any other commercial use made    +
c      of it without the express written permission of the authors.    +
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** The present program version can treat seven types of experimental
c   experimental data, for up to NISTPMX isotopologues of a given species.
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
c (ix) 2'nd virial coefficient data (also only for dPotFit applications)
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
c    vibrational level-v of isotopologue-i of state-s [for NDEGB < 0 case]
c** This subroutine reads in the experimental data on channel-4
c-----------------------------------------------------------------------
cc    INCLUDE 'arrsizes.h'
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c** 'Block' Data Utility routine named: 'arrsizes.h' that governs 
c    array dimensioning in program  dPotFit
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER NISTPMX,NPARMX,NbetaMX,NBOBMX,HPARMX,NDATAMX,
     1  NVIBMX,NBCMX,NSTATEMX,NPNTMX,NROTMX,NCMMAX
c*  NISTPMX  is the maximum number of isotopomers allowed for fit
      PARAMETER (NISTPMX =  12)
c*  NSTATEMX  is maximum no. of electronic states which can be
c             simultaneously fitted to
      PARAMETER (NSTATEMX = 4)
c*  NPARMX  is the largest number of free parameters allowed for fit
c  Since FS origins may be parameters, this is also max. no, data bands
      PARAMETER (NPARMX  = 8000)
c*  NbetaMX  is the largest number of exponent parameters allowed for fit
      PARAMETER (NbetaMX  = 40)
c*  NBOBMX-1  is the highest-order polynomial expansion allowed for the
c      adiabatic or centrifugal Born-Oppenheimer breakdown functions, or 
c      the Lambda-doubling or 2\Sigma splitting radial strength functions
      PARAMETER (NBOBMX  = 15)
c*  HPARMX  is the largest number of Hamiltonian parameters of all types
c    (potential energy, BOB. etc.) for all states.
c           HPARMX >= NSTATEMX*[5 + (NbetaMX+1) + 5*(NBOBMX+1)]
      PARAMETER (HPARMX= NSTATEMX*(5 + (NbetaMX+1) + 5*(NBOBMX+1)))
cc    PARAMETER (HPARMX = 300)
c*  NDATAMX  is largest No. of individual data which may be considered
      PARAMETER (NDATAMX = 35000)
c*  NVIBMX  is the maximum number of vibrational levels of a single
c           state for which data are to be considered
      PARAMETER (NVIBMX    = 200)
** NBCMX  is the maximum number of band constants per vib level to be
c         allowed when doing band constant fits (PSEL= -1)
      PARAMETER (NBCMX = 8) 
c*  NPNTMX  is the largest number of potential data points that can be
c           stored in a single 1D array
      PARAMETER (NPNTMX = 90000)
c*  NROTMX  is the highest order of rotational constants calculated and
c            used for estimating level energies
      PARAMETER (NROTMX = 7)
c*  NCMMAX is the largest number of Cm terms in the MLR or DELR 
c            long-range potential
      PARAMETER (NCMMAX = 12)
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cc    INCLUDE 'BLKPOT.h'
c=======================================================================
c** Effective adiabatic radial potential variables.
      INTEGER BOBCN(NSTATEMX),PSEL(NSTATEMX),MAXMIN(NSTATEMX),
     1 IOMEG(NSTATEMX),Nbeta(NSTATEMX),APSE(NSTATEMX),IFXDE(NSTATEMX),
     2 IFXRE(NSTATEMX),IFXCm(NCMMax,NSTATEMX),
     3 IFXBETA(0:NbetaMX,NSTATEMX),NDATPT(NSTATEMX),NCMM(NSTATEMX),
     4 MMLR(NCMMax,NSTATEMX),nPB(NSTATEMX),nQB(NSTATEMX),pAD(NSTATEMX),
     5 qAD(NSTATEMX),LRad(NSTATEMX),pNA(NSTATEMX),qNA(NSTATEMX),
     6 Pqw(NSTATEMX),IVSR(NSTATEMX),IDSTT(NSTATEMX)
c
      REAL*8 DE(NSTATEMX),RE(NSTATEMX),BETA(0:NbetaMX,NSTATEMX),
     1 yqBETA(NbetaMX,NSTATEMX),BETAFX(NPNTMX,NSTATEMX),RH(NSTATEMX),
     2 RMIN(NSTATEMX),RMAX(NSTATEMX),VLIM(NSTATEMX),EPS(NSTATEMX),
     3 betaINF(NSTATEMX),AGPEF(NSTATEMX),BGPEF(NSTATEMX),
     4 CmVAL(NCMMax,NSTATEMX),CmEFF(NCMMax,NSTATEMX),rhoAB(NSTATEMX),
     5 AA(NSTATEMX),BB(NSTATEMX),RREF(NSTATEMX),ASO(NSTATEMX),
     6 R01(NSTATEMX),Q12(NSTATEMX),RD(NPNTMX,NSTATEMX),
     7 VPOT(NPNTMX,NSTATEMX),dCmA(NCMMax,NSTATEMX),dCmB(NCMMax,NSTATEMX)
c
      COMMON /BLKPOT/DE,RE,BETA,yqBETA,BETAFX,RH,RMIN,RMAX,VLIM,EPS,
     1 betaINF,AGPEF,BGPEF,CmVAL,CmEFF,rhoAB,AA,BB,RREF,ASO,R01,Q12,RD,
     2 VPOT,dCmA,dCmB, BOBCN,PSEL,MAXMIN,IOMEG,Nbeta,APSE,IFXDE,IFXRE,
     3 IFXCm,IFXBETA,NDATPT,NCMM,MMLR,nPB,nQB,pAD,qAD,LRad,pNA,qNA,Pqw,
     4 IVSR,IDSTT
c=======================================================================
cc    INCLUDE 'BLKISOT.h'
c=======================================================================
c** Isotope/isotopologue numbers, masses & BOB mass scaling factors
c** Array ZK carries about the band constants for all levels of all ISOT
      INTEGER NISTP,NDUNMX,AN(2),MN(2,NISTPMX)
c** NDUNMX is a dummy parameter reqd. for portability of READATA
      PARAMETER (NDUNMX=0)
      REAL*8  ZMASS(3,NISTPMX),RSQMU(NISTPMX),RSQMUP(0:NDUNMX,NISTPMX),
     1 RMUP(0:9,NISTPMX),ZMUA(NISTPMX,NSTATEMX),ZMUB(NISTPMX,NSTATEMX),
     2 ZMTA(NISTPMX,NSTATEMX),ZMTB(NISTPMX,NSTATEMX),
     3  ZK(0:NVIBMX,0:NROTMX,NISTPMX,NSTATEMX)
c
      COMMON /BLKISOT/ZMASS,RSQMU,RSQMUP,RMUP,ZMUA,ZMUB,ZMTA,ZMTB,ZK,
     1  NISTP,AN,MN
c=======================================================================
cc    INCLUDE 'BLKDATA.h'
c=======================================================================
c** Type statements & common block for data
      REAL*8  FREQ(NDATAMX),UFREQ(NDATAMX),DFREQ(NDATAMX),TEMP(NDATAMX),
     1                                               YUNC(NDATAMX),Fqb
      INTEGER  COUNTOT,NFS1,NFSTOT,NBANDTOT,IB(NDATAMX),JP(NDATAMX),
     1 JPP(NDATAMX),VP(NPARMX),VPP(NPARMX),EFP(NDATAMX),EFPP(NDATAMX),
     2 TVUP(NDATAMX),TVLW(NDATAMX),FSBAND(NPARMX),IFXFS(NPARMX),
     3 NFS(NPARMX),IEP(NPARMX),IEPP(NPARMX),ISTP(NPARMX),
     4 IFIRST(NPARMX),ILAST(NPARMX),NTV(NSTATEMX,NISTPMX),FSsame,
     5 NTRANS(NPARMX),IBB(NISTPMX,NSTATEMX,9,NPARMX),JMIN(NPARMX),
     6 JMAX(NPARMX)
      CHARACTER*2 NAME(2) 
      CHARACTER*3 SLABL(-6:NSTATEMX)
      CHARACTER*30 BANDNAME(NPARMX)
      COMMON /DATABLK/Fqb,FREQ,UFREQ,YUNC,DFREQ,TEMP,COUNTOT,NFS1,
     1 NFSTOT,NBANDTOT,IB,JP,JPP,VP,VPP,EFP,EFPP,TVUP,TVLW,FSBAND,IFXFS,
     2 NFS,IEP,IEPP,ISTP,IFIRST,ILAST,NTV,FSsame,
     3 NTRANS,IBB,JMIN,JMAX,NAME,SLABL,BANDNAME
c=======================================================================
cc    INCLUDE 'BLKPARAM.h'
c=======================================================================
c** Parameters and count-labels for band constant (PSEL=-1) or term
c   value (PSEL=-2) fits
      REAL*8 TVALUE(NPARMX),ZBC(0:NVIBMX,0:NROTMX,NISTPMX,NSTATEMX),
     1 ZQC(0:NVIBMX,0:NROTMX,NISTPMX,NSTATEMX)
c
      INTEGER NSTATES,NTVALL(0:NSTATEMX),NTVI(NSTATEMX),NTVF(NSTATEMX),
     1 VMIN(NSTATEMX,NISTPMX),VMAX(NSTATEMX,NISTPMX),JTRUNC(NSTATEMX),
     2 EFSEL(NSTATEMX),NBC(0:NVIBMX,NISTPMX,NSTATEMX),
     3 NQC(0:NVIBMX,NISTPMX,NSTATEMX),
     4 BCPARI(0:NVIBMX,NISTPMX,NSTATEMX),
     5 BCPARF(0:NVIBMX,NISTPMX,NSTATEMX),
     6 QCPARI(0:NVIBMX,NISTPMX,NSTATEMX),
     7 QCPARF(0:NVIBMX,NISTPMX,NSTATEMX)
      COMMON /BLKPARAM/TVALUE,ZBC,ZQC,NSTATES,NTVALL,NTVI,NTVF,VMIN,
     1      VMAX,JTRUNC,EFSEL,NBC,NQC,BCPARI,BCPARF,QCPARI,QCPARF
c=======================================================================
cc    INCLUDE 'BLKTYPE.h'
c=======================================================================
c** Type statements & common blocks for characterizing transitions
      REAL*8  AVEUFREQ(NPARMX),MAXUFREQ(NPARMX)
      INTEGER NTRANSFS(NISTPMX,NSTATEMX),
     1  NTRANSVIS(NISTPMX,NSTATEMX,NSTATEMX),
     1  NBANDEL(NISTPMX,NSTATEMX,NSTATEMX),
     2  NTRANSIR(NISTPMX,NSTATEMX),NTRANSMW(NISTPMX,NSTATEMX),
     3  NBANDFS(NISTPMX,NSTATEMX),NBANDVIS(NISTPMX,NSTATEMX),
     4  NBANDIR(NISTPMX,NSTATEMX),NBANDMW(NISTPMX,NSTATEMX),
     5  NVVPP(NISTPMX,NSTATEMX),NWIDTH(NISTPMX,NSTATEMX),
     6  NEBPAS(NISTPMX,NSTATEMX),NVIRIAL(NISTPMX,NSTATEMX),
     7  NAcVIR(NISTPMX,NSTATEMX),NBANDS(NISTPMX)
c
      COMMON /BLKTYPE/AVEUFREQ,MAXUFREQ,NTRANSFS,NTRANSVIS,NTRANSIR,
     1  NTRANSMW,NBANDFS,NBANDEL,NBANDVIS,NBANDIR,NBANDMW,NVVPP,NWIDTH,
     2  NEBPAS,NVIRIAL,NAcVIR,NBANDS
c=======================================================================
c
      INTEGER I,IBN,COUNT,IBAND,
     1  VMX(NSTATEMX),ISOT,NBND,ESP,ESPP,ISTATE,ISTATEE,MN1,MN2,PRINP,
     2  FSOMIT,VMAXesp,VMINesp,VMAXespp,VMINespp,JTRUNCesp,JTRUNCespp
      INTEGER NOWIDTHS,NDAT(0:NVIBMX,NISTPMX,NSTATEMX),PASok(NSTATEMX)
      REAL*8 UCUTOFF,UMIN,TOTUFREQ
      CHARACTER*3 NEF(-1:1)
      CHARACTER*3 LABLP,LABLPP
      CHARACTER*2 OLDSLABL(-6:0)
c-----------------------------------------------------------------------
      DATA NEF/'  f','   ','  e'/
c** maintain compatibility with old labeling method      
      OLDSLABL(-6)=' '                       !! awaiting new data type
      OLDSLABL(-5)='VA'
      OLDSLABL(-4)='VR'
      OLDSLABL(-3)='VV'
      OLDSLABL(-2)='WI'
      OLDSLABL(-1)='PA'
      OLDSLABL(0)='FS'
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
              NAcVIR(ISOT,ISTATE)= 0
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
c** Read all data for each isotopologue at one time.
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
c   AN(2)] identifying the particular isotopologue.  Note that LABLP also
c   identifies the type of data in the 'band' or data-group (see below).
c** LABLP = LABLPP  and  VP = VPP  for a microwave band
c   LABLP = LABLPP  and  VP.ne.VPP  for an infrared band 
c   LABLP = 'FLS'  identifies this data group/band as a fluorescence 
c           series from a single emitting level into vibrational levels
c           of electronic state LABLPP.  In this case: VP is the quantum
c           number v' for the emitting level, while VPP is actually the 
c           rotational quantum number J' for the emitting level and JP
c           [see below] the lower state vibrational quantum number v".
c   LABLP = 'PAS'  identifies this data group/band as a set of binding
c        energies [D-E(v,J,p)] for a given state.  Data Labels as for 'FS'
c   LABLP = 'BV'  identifies this data group/band as a set of Bv values
c           for electronic state LABLPP.  In this case, parameters  VP
c           & VPP are dummy variables, as are EFP, JPP and EFPP [see
c           below],  JP is actually the vibrational quantum number v",
c           FREQ the Bv value & UFREQ its uncertainty
c   LABLP = 'WID'  identifies this data group/band as a set of tunneling 
c           predissociation widths for electronic state LABLPP.  In this
c           case, parameters VP, VPP and EFP are dummy variables, while
c           the predissociating level is identified as: v"=JP, J"=JPP,
c           and parity p"=EFPP.
c   LABLP = 'VVV' to identify this as a set of potential fx. values 
c           e.g., ab initio values for the high repulsive wall. In this
c            case, parameters VP, VPP are dummy variables.
c   LABLP = 'VIR' identifies this data group/band as a set of virial 
c           coefficients for electronic state LABLPP.  In this case, 
c           parameters VP, VPP are dummy variables.
c   LABLP = 'VAC' identifies this data group/band as a set of virial 
c           coefficients for electronic state LABLPP.  In this case, 
c           parameters VP, VPP are dummy variables.
c** STOP reading when run out of bands OR when read-in VPP is negative   
c-----------------------------------------------------------------------
      IF((PRINP.EQ.2).OR.(PRINP.EQ.-2)) THEN
      READ(4,*,END=40) VP(IBAND), VPP(IBAND), LABLP, LABLPP, MN1, MN2,
     1                 BANDNAME(IBAND)     
      ELSE
      READ(4,*,END=40) VP(IBAND), VPP(IBAND), LABLP, LABLPP, MN1, MN2
      ENDIF
c-----------------------------------------------------------------------
      IF(VP(IBAND).LT.0) GO TO 40
      IEP(IBAND)= -99
      IEPP(IBAND)= -99
      DO  I= -6, 0
          IF(LABLP.EQ.OLDSLABL(I)) LABLP= SLABL(I)
          IF(LABLPP.EQ.OLDSLABL(I)) LABLPP= SLABL(I)
          ENDDO
      DO  I= -6, NSTATES
          IF(LABLP.EQ.SLABL(I)) IEP(IBAND)= I
          IF(LABLPP.EQ.SLABL(I)) IEPP(IBAND)= I
          ENDDO
c** Check that this isotopologue is one of those chosen to be fitted ...
      ISOT= 0
      DO  I= 1,NISTP
          IF((MN1.EQ.MN(1,I)).AND.(MN2.EQ.MN(2,I))) ISOT= I
          ENDDO
      ISTP(IBAND)= ISOT
      ESP= IEP(IBAND)
      ESPP= IEPP(IBAND)
      IF(IEP(IBAND).EQ.-3) THEN
c** For case in which the 'data' are potential function value(s) in cm-1 ...
          COUNT= COUNT+ 1
          IFIRST(IBAND)= COUNT
c...  TEMP(i)= r(i) ;  FREQ(i)= V(r(i))  ; UFREQ= unc{V(r(i))}
c----------------------------------------------------------------------
   12     READ(4,*) TEMP(COUNT),FREQ(COUNT),UFREQ(COUNT)
c----------------------------------------------------------------------
          YUNC(COUNT)= UFREQ(COUNT)
c ... a negative input distance implies end of potential energy data set
          IF(TEMP(COUNT).GT.0.d0) THEN
c ... if this isotope or state not considered, ignore this datum
              IF((ISOT.LE.0).OR.(ESPP.LT.-6)) GOTO 12
c ... if no potential used, ignore this datum  
              IF(PSEL(ESPP).LT.0) GOTO 12  
              IB(COUNT)= IBAND
              COUNT= COUNT+1
              GOTO 12
            ELSE
              GOTO 18
            ENDIF
          ENDIF
      IF((IEP(IBAND).EQ.-4).OR.(IEP(IBAND).EQ.-5)) THEN
c** For case in which the data are virial coefficients
          COUNT= COUNT+ 1
          IFIRST(IBAND)= COUNT
c... TEMP(i)= temperature[K],  FREQ(i)= Bvir(i),  UFREQ(i)= unc{Bvir(i)}
c----------------------------------------------------------------------
   14     READ(4,*) TEMP(COUNT),FREQ(COUNT),UFREQ(COUNT)
c----------------------------------------------------------------------
          YUNC(COUNT)= UFREQ(COUNT)
c ... negative input 'temperature' implies end of virial/PE data set
          IF(TEMP(COUNT).GT.0.d0) THEN
c ... if this isotope or state not considered, ignore this datum
              IF((ISOT.LE.0).OR.(ESPP.LT.-6)) GOTO 14
c ... if no potential used, ignore this datum  
              IF(PSEL(ESPP).LT.0) GOTO 14  
              IB(COUNT)= IBAND
              COUNT= COUNT+1
              GOTO 14
            ELSE       !! for 'TEMP'.LE.0
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
      NTRANS(IBAND)= 0
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
              VMAXespp= INT((VMAX(ESPP,ISOT)+0.5d0)/RSQMU(ISOT)-0.5d0) !! added
              VMINespp= INT((VMIN(ESPP,ISOT)+0.5d0)/RSQMU(ISOT)-0.5d0) !! added
              JTRUNCespp= INT(JTRUNC(ESPP)/RSQMU(ISOT))
              ENDIF
cc        VMAXesp= VMAX(ESPP,ISOT)         ?????? Huh ??????/
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
c** If this band is not for one of the chosen isotopologues or states
c  or 'property' datum w. no PEC, omit this datum from the fit
      IF((ISOT.EQ.0).OR.(ESPP.LT.-6)) GO TO 15
      IF((PSEL(ESPP).LT.0).AND.(ESP.LT.0)) GO TO 15
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
      NTRANS(IBAND)= ILAST(IBAND)-IFIRST(IBAND)+1
      IF(NTRANS(IBAND).GT.0) THEN
c** Treat PAS data as Fluorescence series unless  PASok > 0
          IF((IEP(IBAND).EQ.-1).AND.(PASok(IEPP(IBAND)).LE.0)) 
     1                                                    IEP(IBAND)=0
          IF((NTRANS(IBAND).EQ.1).AND.(LABLP.EQ.'FS')) THEN
c** Ignore any fluorescence series consisting of only one datum
              COUNT= COUNT-1
              IBAND= IBAND-1
              FSOMIT= FSOMIT+1
              GOTO 10
              ENDIF
          AVEUFREQ(IBAND)= TOTUFREQ/NTRANS(IBAND)
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
          NTRANSFS(ISOT,ESPP)= NTRANSFS(ISOT,ESPP)+NTRANS(IBAND)
c ... and then set up labels/ranges/properties for each band
          IBB(ISOT,ESPP,1,NBND)= IBAND
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
          NTRANSVIS(ISOT,ESP,ESPP)=
     1                     NTRANSVIS(ISOT,ESP,ESPP)+NTRANS(IBAND)
c ... and then set up labels/ranges/properties for each of them
          IBB(ISOT,ESPP,2,NBND)= IBAND
          ENDIF 
c
      IF((ESP.EQ.ESPP).AND.(VP(IBAND).NE.VPP(IBAND))) THEN
c** For an Infrared band of electronic state  s=ESPP=ESP
c** First cumulatively count the number of IR bands & transitions
          NBANDIR(ISOT,ESPP)= NBANDIR(ISOT,ESPP)+1
          NBND= NBANDIR(ISOT,ESPP)
          NTRANSIR(ISOT,ESPP)= NTRANSIR(ISOT,ESPP)+NTRANS(IBAND) 
c ... and then set up labels/ranges/properties for each of them
          IBB(ISOT,ESPP,3,NBND)= IBAND
          ENDIF
c
      IF((ESP.EQ.ESPP).AND.(VP(IBAND).EQ.VPP(IBAND))) THEN
c** For Microwave transitions in electronic state  s=ESPP=ESP
c** First cumulatively count the number of MW bands & transitions
          NBANDMW(ISOT,ESPP)= NBANDMW(ISOT,ESPP)+1
          NBND= NBANDMW(ISOT,ESPP)
          NTRANSMW(ISOT,ESPP)= NTRANSMW(ISOT,ESPP)+NTRANS(IBAND)
c ... and then set up labels/ranges/properties for each of them
          IBB(ISOT,ESPP,4,NBND)= IBAND
          ENDIF
c
c** NOTE ... in IBB array a last index counts bands of this type for 
c  this isotopologue of this electronic state.  Expect to find all 
c  potential fx. values, virial coeficients, Tunneling Widths, PAS 
c  binding energies, virial coeffts, ... etc. in a single group.
          IF(ESP.EQ.-5) THEN
c** Data are Accoustic Virial Coefficients for electronic state IEPP= ESPP
          NAcVIR(ISOT,ESPP)= NTRANS(IBAND)
          IBB(ISOT,ESPP,9,1)= IBAND
          ENDIF
c
          IF(ESP.EQ.-4) THEN
c** Data are pressure Virial Coefficients for electronic state IEPP= ESPP
          NVIRIAL(ISOT,ESPP)= NTRANS(IBAND)
          IBB(ISOT,ESPP,8,1)= IBAND
          ENDIF
c
          IF(ESP.EQ.-3) THEN
c** Data are not transition energies, but rather values of the potential
c  function at particular distances for electronic state s=IEPP  
              WRITE(6,612) LABLPP,ISOT
              NVVPP(ISOT,ESPP)= NTRANS(IBAND)
              IBB(ISOT,ESPP,5,1)= IBAND
              ENDIF
c
      IF(ESP.EQ.-2) THEN
c** Data are tunneling predissociation linewidths (in cm-1) for levels
c  of electronic state IEPP=ESPP
ccc       IF((NWIDTH(ISOT,ESPP).GT.0).AND.(NTRANS(IBAND).GT.0)) THEN
              WRITE(6,626) ESPP,ISOT
ccc           STOP
ccc           ENDIF
          NWIDTH(ISOT,ESPP)= NTRANS(IBAND)
          IBB(ISOT,ESPP,6,1)= IBAND
          ENDIF
c
      IF(ESP.EQ.-1) THEN
c** Data are PhotoAssociation Binding Energies (in cm-1) for levels
c  of electronic state IEPP=ESPP
          WRITE(6,636) LABLPP,ISOT
          NEBPAS(ISOT,ESPP)= NTRANS(IBAND)
          IBB(ISOT,ESPP,7,1)= IBAND
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
c** Write a summary of the data, one isotopologue at a time.
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
          IBN=IBB(ISOT,ISTATE,4,I)
              WRITE(6,606) VP(IBN),VPP(IBN),NTRANS(IBN),JMIN(IBN),
     1               JMAX(IBN),AVEUFREQ(IBN),MAXUFREQ(IBN)
              ENDDO
          ENDIF
c
      IF(NTRANSIR(ISOT,ISTATE).GT.0)THEN
c** Book-keeping for Infrared data
          WRITE(6,608) NTRANSIR(ISOT,ISTATE),SLABL(ISTATE),(NAME(I),
     1                          MN(I,ISOT),I=1,2),NBANDIR(ISOT,ISTATE)
          DO  I= 1,NBANDIR(ISOT,ISTATE)
          IBN=IBB(ISOT,ISTATE,3,I)
              WRITE(6,606) VP(IBN),VPP(IBN),NTRANS(IBN),
     1                  JMIN(IBN),JMAX(IBN),AVEUFREQ(IBN),MAXUFREQ(IBN)
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
              IBN=IBB(ISOT,ISTATE,2,I)
                  IF(IEP(IBN).EQ.ISTATEE) THEN
                      WRITE(6,606) VP(IBN),VPP(IBN),NTRANS(IBN),
     1                  JMIN(IBN),JMAX(IBN),AVEUFREQ(IBN),MAXUFREQ(IBN)
                      ENDIF
                  ENDDO
              ENDIF
          ENDDO
      IF(NTRANSFS(ISOT,ISTATE).GT.0)THEN
c** Book-keeping for Fluorescence data
          WRITE(6,614) NTRANSFS(ISOT,ISTATE),SLABL(ISTATE),(NAME(I),
     1                          MN(I,ISOT),I=1,2),NBANDFS(ISOT,ISTATE)
          DO  I= 1,NBANDFS(ISOT,ISTATE)
          IBN = IBB(ISOT,ISTATE,1,I)
              WRITE(6,616) VP(IBN),VPP(IBN),
     1                  NEF(EFP(IFIRST(IBB(ISOT,ISTATE,1,I)))),
     2     NTRANS(IBN),JMIN(IBN),JMAX(IBN),AVEUFREQ(IBN),MAXUFREQ(IBN)
              ENDDO
          ENDIF
      IF(NVVPP(ISOT,ISTATE).GT.0)THEN
c** Book-keeping for potential function values as data ....
          WRITE(6,618) NVVPP(ISOT,ISTATE),SLABL(ISTATE),(NAME(I),
     1                                               MN(I,ISOT),I=1,2)
          IBN=IBB(ISOT,ISTATE,5,1)
          WRITE(6,620) NTRANS(IBN),JMIN(IBN),JMAX(IBN),AVEUFREQ(IBN),
     2                               MAXUFREQ(IBN)
          ENDIF
      IF(NWIDTH(ISOT,ISTATE).GT.0) THEN
c** Book-keeping for  Tunneling Width  data
          WRITE(6,628) NWIDTH(ISOT,ISTATE),SLABL(ISTATE),(NAME(I),
     1                                               MN(I,ISOT),I=1,2)
          IBN=IBB(ISOT,ISTATE,6,1)
          WRITE(6,630) NTRANS(IBN),JMIN(IBN),JMAX(IBN),AVEUFREQ(IBN),
     3                                MAXUFREQ(IBN)
          ENDIF
      IF(NEBPAS(ISOT,ISTATE).GT.0) THEN
c** Book-keeping for  PAS Binding Energy  data
          WRITE(6,632) NEBPAS(ISOT,ISTATE),SLABL(ISTATE),(NAME(I),
     1                                               MN(I,ISOT),I=1,2)
          IBN=IBB(ISOT,ISTATE,7,1)
          WRITE(6,630) NTRANS(IBN),JMIN(IBN),JMAX(IBN),AVEUFREQ(IBN),
     3                                MAXUFREQ(IBN)
          ENDIF
      IF(NVIRIAL(ISOT,ISTATE).GT.0) THEN
c** Book-keeping for Virial data
          WRITE(6,642) NVIRIAL(ISOT,ISTATE), SLABL(ISTATE), 
     1                                      (NAME(I),MN(I,ISOT),I=1,2)
          ENDIF
      IF(NAcVIR(ISOT,ISTATE).GT.0) THEN
c** Book-keeping for Accoustic Virial data
          WRITE(6,644) NAcVIR(ISOT,ISTATE), SLABL(ISTATE), 
     1                                      (NAME(I),MN(I,ISOT),I=1,2)
          ENDIF
   50 CONTINUE
      IF(ISOT.LT.NISTP) THEN
c** If NISTP > 1, return to print data summaries for other isotopologues
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
  605 FORMAT(7x,'and State ',A3,' data with  J < JTRUNC=',I4)
  607 FORMAT(7x,'and State ',A3,' data with  J > JTRUNC=',I4)
  611 FORMAT(29x,'or  v  outside range',i3,'  to',i4,'   for   ISOT=',
     1  i2:)
  602 FORMAT(/1x,20('===')/'  *** Input data for',i5,' bands/series of '
     1  ,A2,'(',I3,')-',A2,'(',I3,') ***'/1x,20('==='))
  604 FORMAT(1x,28('--')/I5,' State ',A3,1x,A2,'(',I3,')-',A2,'(',I3,
     1 ') MW transitions in',i4,' sets'/1x,28('--')/"   v'  ",
     1 'v"  #data   Jmin   Jmax  Avge.Unc.  Max.Unc.'/1x,25('--'))
  606 FORMAT(I4,I4,3I7,1x,1P2D10.1)
  608 FORMAT(1x,32('--')/I6,' State ',A3,1x,A2,'(',I3,')-',A2,'(',I3,
     1 ') InfraRed transitions in',I4,' bands'/1x,32('--')/
     2 "   v'  ",'v"  #data   Jmin   Jmax  Avge.Unc.  Max.Unc.'/
     3 1x,25('--'))
  609 FORMAT(/' *** ERROR *** Dimension allocated for number of bands ex
     1ceeded:'/' (IBAND=',i4,') > (NBANDMX=',i4,')   so truncate input a
     2nd TRY to continue ...')
  610 FORMAT(/1x,35('==')/I6,1x,A2,'(',I3,')-',A2,'(',i3,')  {State ',
     1  A3,'}--{State ',A3,'} Transitions in',i4,' Bands'/1x,35('--')/
     2 "   v'",'  v"  #data   Jmin   Jmax  Avge.Unc.  Max.Unc.'/
     3 1x,25('--'))
  612 FORMAT(/" NOTE that read-in potential fx. values for   ISTATE= ",
     1   A3,'   ISOT=',i2/32x,' must be input as a single "band" or data
     1 group')
  614 FORMAT(1x,38('==')/I5,' Fluorescence transitions into State ',
     1 A3,2x,A2,'(',I3,')-',A2,'(',I3,')  in',i5,' series'/
     2 1x,38('==')/"   v'  j' p' ",'#data  v"min  v"max  Avge.Unc.  Max.
     3Unc.'/1x,51('-'))
  616 FORMAT(2I4,A3,I6,2I7,1x,1P2D10.1)
  618 FORMAT(1x,65('=')/1x,I3,' State ',A3,1x,A2,'(',I3,')-',A2,'(',I3,
     1 ') potential fx values treated as independent data'/1x,24('--')/
     2 '  #values   r(min)  r(max)  Avge.Unc.   Max.Unc.'/1x,24('--'))
  620 FORMAT(I7,I9,I8,3x,1P2D11.1)
  622 FORMAT(1x,25('===')/1x,25('==='))
  626 FORMAT(/" NOTE that all read-in Tunneling Widths for   ISTATE=",
     1 i2,'   ISOT=',i2/10x,' must be in a single "band" or data group')
cc626 FORMAT(/" *** STOP INPUT *** and put all read-in Tunneling Widths'
cc   1  '  for   ISTATE=",i2,'   ISOT=',i2/ 
cc   2  10x,'into one "band" or data group.')
  628 FORMAT(1x,61('=')/1x,I3,' State ',A3,1x,A2,'(',I3,')-',A2,'(',I3,
     1 ') Tunneling Widths included as data'/
     2 1x,61('-')/'  #values   v(min)  v(max)   Avge.Unc.   Max.Unc.'/
     3 1x,24('--'))
  630 FORMAT(I7,I9,I8,2x,1P2D11.1)
  632 FORMAT(1x,70('=')/I4,' State ',A3,1x,A2,'(',I3,')-',A2,'(',I3,
     1 ') PAS Binding Energies included in data set'/
     2 1x,70('-')/'  #values   v(min)  v(max)   Avge.Unc.   Max.Unc.'/
     3 1x,24('--'))
  636 FORMAT(/' NOTE that all read-in PAS Binding Energies for   ISTATE=
     1 ',a2,'  ISOT=',i2/10x,' must be in a single "band" or data group'
     2 )
  640 FORMAT(/' *** Input Data Count reaches',i6,' which EXCEEDS ARRAY L
     1IMIT of',i6)
  642 FORMAT(1x,70('=')/I4,' State ',A3,1x,A2,'(',I3,')-',A2,
     1 '(',I3,') Virial coefficients included in data set' )
  644 FORMAT(1x,70('=')/I4,' State ',A3,1x,A2,'(',I3,')-',A2,
     1 '(',I3,') Accoustic Virial coefficients included in data set' )
  650 FORMAT(/' Data input IGNORES',i4,' fluorescence series consisting'
     1 ,' of only  onee  line!')
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE TVSORT(ISTATE,NPTOT,VMAX,NTVALL,NTVSSTAT,TVNAME)
c***********************************************************************
c** Subroutine to sort through global data file, and for each isotopologue
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
c  NTVALL  is the total number of term value parameters for this state
c  NTVSSTAT  is the total number of term values this state associated 
c           with only a single transition
c  TVNAME(j)  is the alphameric name identifying term value parameter j
c
c** Internally
c-------------
c  NLV(v,J.p) * initially, counts transitions for level {v,J,p} of a 
c                          given isotopologue
c           * later reset it as the parameter index for that term value
c-----------------------------------------------------------------------
cc    INCLUDE 'arrsizes.h'
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c** 'Block' Data Utility routine named: 'arrsizes.h' that governs 
c    array dimensioning in program  dPotFit
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER NISTPMX,NPARMX,NbetaMX,NBOBMX,HPARMX,NDATAMX,
     1  NVIBMX,NBCMX,NSTATEMX,NPNTMX,NROTMX,NCMMAX
c*  NISTPMX  is the maximum number of isotopomers allowed for fit
      PARAMETER (NISTPMX =  12)
c*  NSTATEMX  is maximum no. of electronic states which can be
c             simultaneously fitted to
      PARAMETER (NSTATEMX = 4)
c*  NPARMX  is the largest number of free parameters allowed for fit
c  Since FS origins may be parameters, this is also max. no, data bands
      PARAMETER (NPARMX  = 8000)
c*  NbetaMX  is the largest number of exponent parameters allowed for fit
      PARAMETER (NbetaMX  = 40)
c*  NBOBMX-1  is the highest-order polynomial expansion allowed for the
c      adiabatic or centrifugal Born-Oppenheimer breakdown functions, or 
c      the Lambda-doubling or 2\Sigma splitting radial strength functions
      PARAMETER (NBOBMX  = 15)
c*  HPARMX  is the largest number of Hamiltonian parameters of all types
c    (potential energy, BOB. etc.) for all states.
c           HPARMX >= NSTATEMX*[5 + (NbetaMX+1) + 5*(NBOBMX+1)]
      PARAMETER (HPARMX= NSTATEMX*(5 + (NbetaMX+1) + 5*(NBOBMX+1)))
cc    PARAMETER (HPARMX = 300)
c*  NDATAMX  is largest No. of individual data which may be considered
      PARAMETER (NDATAMX = 35000)
c*  NVIBMX  is the maximum number of vibrational levels of a single
c           state for which data are to be considered
      PARAMETER (NVIBMX    = 200)
** NBCMX  is the maximum number of band constants per vib level to be
c         allowed when doing band constant fits (PSEL= -1)
      PARAMETER (NBCMX = 8) 
c*  NPNTMX  is the largest number of potential data points that can be
c           stored in a single 1D array
      PARAMETER (NPNTMX = 90000)
c*  NROTMX  is the highest order of rotational constants calculated and
c            used for estimating level energies
      PARAMETER (NROTMX = 7)
c*  NCMMAX is the largest number of Cm terms in the MLR or DELR 
c            long-range potential
      PARAMETER (NCMMAX = 12)
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cc    INCLUDE 'BLKISOT.h'
c=======================================================================
c** Isotope/isotopologue numbers, masses & BOB mass scaling factors
c** Array ZK carries about the band constants for all levels of all ISOT
      INTEGER NISTP,NDUNMX,AN(2),MN(2,NISTPMX)
c** NDUNMX is a dummy parameter reqd. for portability of READATA
      PARAMETER (NDUNMX=0)
      REAL*8  ZMASS(3,NISTPMX),RSQMU(NISTPMX),RSQMUP(0:NDUNMX,NISTPMX),
     1 RMUP(0:9,NISTPMX),ZMUA(NISTPMX,NSTATEMX),ZMUB(NISTPMX,NSTATEMX),
     2 ZMTA(NISTPMX,NSTATEMX),ZMTB(NISTPMX,NSTATEMX),
     3  ZK(0:NVIBMX,0:NROTMX,NISTPMX,NSTATEMX)
c
      COMMON /BLKISOT/ZMASS,RSQMU,RSQMUP,RMUP,ZMUA,ZMUB,ZMTA,ZMTB,ZK,
     1  NISTP,AN,MN
c=======================================================================
cc    INCLUDE 'BLKDATA.h'
c=======================================================================
c** Type statements & common block for data
      REAL*8  FREQ(NDATAMX),UFREQ(NDATAMX),DFREQ(NDATAMX),TEMP(NDATAMX),
     1                                               YUNC(NDATAMX),Fqb
      INTEGER  COUNTOT,NFS1,NFSTOT,NBANDTOT,IB(NDATAMX),JP(NDATAMX),
     1 JPP(NDATAMX),VP(NPARMX),VPP(NPARMX),EFP(NDATAMX),EFPP(NDATAMX),
     2 TVUP(NDATAMX),TVLW(NDATAMX),FSBAND(NPARMX),IFXFS(NPARMX),
     3 NFS(NPARMX),IEP(NPARMX),IEPP(NPARMX),ISTP(NPARMX),
     4 IFIRST(NPARMX),ILAST(NPARMX),NTV(NSTATEMX,NISTPMX),FSsame,
     5 NTRANS(NPARMX),IBB(NISTPMX,NSTATEMX,9,NPARMX),JMIN(NPARMX),
     6 JMAX(NPARMX)
      CHARACTER*2 NAME(2) 
      CHARACTER*3 SLABL(-6:NSTATEMX)
      CHARACTER*30 BANDNAME(NPARMX)
      COMMON /DATABLK/Fqb,FREQ,UFREQ,YUNC,DFREQ,TEMP,COUNTOT,NFS1,
     1 NFSTOT,NBANDTOT,IB,JP,JPP,VP,VPP,EFP,EFPP,TVUP,TVLW,FSBAND,IFXFS,
     2 NFS,IEP,IEPP,ISTP,IFIRST,ILAST,NTV,FSsame,
     3 NTRANS,IBB,JMIN,JMAX,NAME,SLABL,BANDNAME
c=======================================================================
c
      INTEGER I,J,P,IBAND,ISOT,ISTATE,NPTOT,LOWEST,
     1 VMAX(NSTATEMX,NISTPMX),NLV(0:NVIBMX,0:NVIBMX,-1:1),
     2 NTVS(NSTATEMX,NISTPMX),NTVALL(0:NSTATEMX),NTVSSTAT
      CHARACTER*24 TVNAME(NPARMX)
c=======================================================================
      WRITE(6,600) SLABL(ISTATE) 
      LOWEST= 1
      IF(ISTATE.GT.1) LOWEST= 0
      NTVALL(ISTATE)= 0
      NTVSSTAT= 0
      DO  ISOT= 1, NISTP
c** First ... zero transition counter array for this isotopologue
          DO  I= 0, VMAX(ISTATE,ISOT)
              DO  J= 0, NVIBMX 
                  DO  P= -1,1
                      NLV(I,J,P)= 0
                      ENDDO
                  ENDDO
              ENDDO
          DO  IBAND= 1, NBANDTOT
c** Then ... search for bands involving isotopologue ISOT in this state
              IF(((IEP(IBAND).EQ.ISTATE).OR.(IEPP(IBAND).EQ.ISTATE))
     1          .AND.(ISTP(IBAND).EQ.ISOT).AND.(IEP(IBAND).GE.0)) THEN
                  DO  I= IFIRST(IBAND), ILAST(IBAND)
c ... for each such band, loop over all transitions, and increment NLV 
c     for each {v,J,p} level encountered in a transision
                      IF(IEP(IBAND).EQ.ISTATE) THEN
                          IF(JP(I).GT.NVIBMX) THEN
c ... check for array dimension overruns
                              WRITE(6,602) ISTATE,ISOT,JP(I),NVIBMX
                              STOP
                              ENDIF
                          NLV(VP(IBAND),JP(I),EFP(I))= 
     1                                  NLV(VP(IBAND),JP(I),EFP(I))+ 1
                          ENDIF
                      IF(IEPP(IBAND).EQ.ISTATE) THEN
                          IF(JPP(I).GT.NVIBMX) THEN
                              WRITE(6,604) ISTATE,ISOT,JPP(I),NVIBMX
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
          NTVSSTAT= NTVSSTAT+ NTVS(ISTATE,ISOT)
          ENDDO
c
      RETURN
  600 FORMAT(/' For State ',A3,'  fit to individual term values for each
     1  {v,J,p,isot}'/1x,6('******'))
  602 FORMAT(/' *** ARRAY DIMENSION PROBLEM ***  JP(ISTATE)=',i2,
     1  ',ISOT=',I2,')=',i3,'  greater than  NVIBMX=',i4)
  604 FORMAT(/' *** ARRAY DIMENSION PROBLEM ***  JPP(ISTATE)=',i2,
     1  ',ISOT=',I2,')=',i3,'  greater than  NVIBMX=',i4)
  606 FORMAT(/'  Absolute zero of energy is fixed at level {v=',i3,
     1 ', J=',i3,', p=',i2,'}'/1x,12('**'),10x,'of isotopologue ',i2,
     2 ' of  State ',A3)
  608 FORMAT(' For ',A2,'(',i3,')-',A2,'(',I3,')  fit to',i5,
     1 ' T(v,J,p) term values,'/20x,'of which',i5,' are involved in only
     2 one transition')
  700 FORMAT("'",'T(',A3,':',i3,',',i3,',',SP,i2,';',SS,i2,')',I4,"'")
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c**********************************************************************N(
      SUBROUTINE READPOT(ISTATE,SLABL)
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
c   SLABL   is the three-character label identifying that state
c-----------------------------------------------------------------------
cc    INCLUDE 'arrsizes.h'
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c** 'Block' Data Utility routine named: 'arrsizes.h' that governs 
c    array dimensioning in program  dPotFit
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER NISTPMX,NPARMX,NbetaMX,NBOBMX,HPARMX,NDATAMX,
     1  NVIBMX,NBCMX,NSTATEMX,NPNTMX,NROTMX,NCMMAX
c*  NISTPMX  is the maximum number of isotopomers allowed for fit
      PARAMETER (NISTPMX =  12)
c*  NSTATEMX  is maximum no. of electronic states which can be
c             simultaneously fitted to
      PARAMETER (NSTATEMX = 4)
c*  NPARMX  is the largest number of free parameters allowed for fit
c  Since FS origins may be parameters, this is also max. no, data bands
      PARAMETER (NPARMX  = 8000)
c*  NbetaMX  is the largest number of exponent parameters allowed for fit
      PARAMETER (NbetaMX  = 40)
c*  NBOBMX-1  is the highest-order polynomial expansion allowed for the
c      adiabatic or centrifugal Born-Oppenheimer breakdown functions, or 
c      the Lambda-doubling or 2\Sigma splitting radial strength functions
      PARAMETER (NBOBMX  = 15)
c*  HPARMX  is the largest number of Hamiltonian parameters of all types
c    (potential energy, BOB. etc.) for all states.
c           HPARMX >= NSTATEMX*[5 + (NbetaMX+1) + 5*(NBOBMX+1)]
      PARAMETER (HPARMX= NSTATEMX*(5 + (NbetaMX+1) + 5*(NBOBMX+1)))
cc    PARAMETER (HPARMX = 300)
c*  NDATAMX  is largest No. of individual data which may be considered
      PARAMETER (NDATAMX = 35000)
c*  NVIBMX  is the maximum number of vibrational levels of a single
c           state for which data are to be considered
      PARAMETER (NVIBMX    = 200)
** NBCMX  is the maximum number of band constants per vib level to be
c         allowed when doing band constant fits (PSEL= -1)
      PARAMETER (NBCMX = 8) 
c*  NPNTMX  is the largest number of potential data points that can be
c           stored in a single 1D array
      PARAMETER (NPNTMX = 90000)
c*  NROTMX  is the highest order of rotational constants calculated and
c            used for estimating level energies
      PARAMETER (NROTMX = 7)
c*  NCMMAX is the largest number of Cm terms in the MLR or DELR 
c            long-range potential
      PARAMETER (NCMMAX = 12)
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cc    INCLUDE 'BLKISOT.h'
c=======================================================================
c** Isotope/isotopologue numbers, masses & BOB mass scaling factors
c** Array ZK carries about the band constants for all levels of all ISOT
      INTEGER NISTP,NDUNMX,AN(2),MN(2,NISTPMX)
c** NDUNMX is a dummy parameter reqd. for portability of READATA
      PARAMETER (NDUNMX=0)
      REAL*8  ZMASS(3,NISTPMX),RSQMU(NISTPMX),RSQMUP(0:NDUNMX,NISTPMX),
     1 RMUP(0:9,NISTPMX),ZMUA(NISTPMX,NSTATEMX),ZMUB(NISTPMX,NSTATEMX),
     2 ZMTA(NISTPMX,NSTATEMX),ZMTB(NISTPMX,NSTATEMX),
     3  ZK(0:NVIBMX,0:NROTMX,NISTPMX,NSTATEMX)
c
      COMMON /BLKISOT/ZMASS,RSQMU,RSQMUP,RMUP,ZMUA,ZMUB,ZMTA,ZMTB,ZK,
     1  NISTP,AN,MN
c=======================================================================
cc    INCLUDE 'BLKCOUNT.h'
c=======================================================================
c    Block data file  BLKCOUNT.h
c=======================================================================
c** Counters for numbers of potential parameters of different types for 
c   each state
      INTEGER  TOTPOTPAR,POTPARI(NSTATEMX),POTPARF(NSTATEMX),
     1  UAPARI(NSTATEMX),UAPARF(NSTATEMX),UBPARI(NSTATEMX),
     2  UBPARF(NSTATEMX),TAPARI(NSTATEMX),TAPARF(NSTATEMX),
     3  TBPARI(NSTATEMX),TBPARF(NSTATEMX),LDPARI(NSTATEMX),
     4  LDPARF(NSTATEMX),HPARF(NSTATEMX),OSEL(NSTATEMX)
c
      COMMON /BLKCOUNT/TOTPOTPAR,POTPARI,POTPARF,UAPARI,UAPARF,UBPARI,
     1  UBPARF,TAPARI,TAPARF,TBPARI,TBPARF,LDPARI,LDPARF,HPARF,OSEL
c=======================================================================
cc    INCLUDE 'BLKPOT.h'
c=======================================================================
c** Effective adiabatic radial potential variables.
      INTEGER BOBCN(NSTATEMX),PSEL(NSTATEMX),MAXMIN(NSTATEMX),
     1 IOMEG(NSTATEMX),Nbeta(NSTATEMX),APSE(NSTATEMX),IFXDE(NSTATEMX),
     2 IFXRE(NSTATEMX),IFXCm(NCMMax,NSTATEMX),
     3 IFXBETA(0:NbetaMX,NSTATEMX),NDATPT(NSTATEMX),NCMM(NSTATEMX),
     4 MMLR(NCMMax,NSTATEMX),nPB(NSTATEMX),nQB(NSTATEMX),pAD(NSTATEMX),
     5 qAD(NSTATEMX),LRad(NSTATEMX),pNA(NSTATEMX),qNA(NSTATEMX),
     6 Pqw(NSTATEMX),IVSR(NSTATEMX),IDSTT(NSTATEMX)
c
      REAL*8 DE(NSTATEMX),RE(NSTATEMX),BETA(0:NbetaMX,NSTATEMX),
     1 yqBETA(NbetaMX,NSTATEMX),BETAFX(NPNTMX,NSTATEMX),RH(NSTATEMX),
     2 RMIN(NSTATEMX),RMAX(NSTATEMX),VLIM(NSTATEMX),EPS(NSTATEMX),
     3 betaINF(NSTATEMX),AGPEF(NSTATEMX),BGPEF(NSTATEMX),
     4 CmVAL(NCMMax,NSTATEMX),CmEFF(NCMMax,NSTATEMX),rhoAB(NSTATEMX),
     5 AA(NSTATEMX),BB(NSTATEMX),RREF(NSTATEMX),ASO(NSTATEMX),
     6 R01(NSTATEMX),Q12(NSTATEMX),RD(NPNTMX,NSTATEMX),
     7 VPOT(NPNTMX,NSTATEMX),dCmA(NCMMax,NSTATEMX),dCmB(NCMMax,NSTATEMX)
c
      COMMON /BLKPOT/DE,RE,BETA,yqBETA,BETAFX,RH,RMIN,RMAX,VLIM,EPS,
     1 betaINF,AGPEF,BGPEF,CmVAL,CmEFF,rhoAB,AA,BB,RREF,ASO,R01,Q12,RD,
     2 VPOT,dCmA,dCmB, BOBCN,PSEL,MAXMIN,IOMEG,Nbeta,APSE,IFXDE,IFXRE,
     3 IFXCm,IFXBETA,NDATPT,NCMM,MMLR,nPB,nQB,pAD,qAD,LRad,pNA,qNA,Pqw,
     4 IVSR,IDSTT
c=======================================================================
cc    INCLUDE 'BLKPARAM.h'
c=======================================================================
c** Parameters and count-labels for band constant (PSEL=-1) or term
c   value (PSEL=-2) fits
      REAL*8 TVALUE(NPARMX),ZBC(0:NVIBMX,0:NROTMX,NISTPMX,NSTATEMX),
     1 ZQC(0:NVIBMX,0:NROTMX,NISTPMX,NSTATEMX)
c
      INTEGER NSTATES,NTVALL(0:NSTATEMX),NTVI(NSTATEMX),NTVF(NSTATEMX),
     1 VMIN(NSTATEMX,NISTPMX),VMAX(NSTATEMX,NISTPMX),JTRUNC(NSTATEMX),
     2 EFSEL(NSTATEMX),NBC(0:NVIBMX,NISTPMX,NSTATEMX),
     3 NQC(0:NVIBMX,NISTPMX,NSTATEMX),
     4 BCPARI(0:NVIBMX,NISTPMX,NSTATEMX),
     5 BCPARF(0:NVIBMX,NISTPMX,NSTATEMX),
     6 QCPARI(0:NVIBMX,NISTPMX,NSTATEMX),
     7 QCPARF(0:NVIBMX,NISTPMX,NSTATEMX)
      COMMON /BLKPARAM/TVALUE,ZBC,ZQC,NSTATES,NTVALL,NTVI,NTVF,VMIN,
     1      VMAX,JTRUNC,EFSEL,NBC,NQC,BCPARI,BCPARF,QCPARI,QCPARF
c=======================================================================
cc    INCLUDE 'BLKBOB.h'
c=======================================================================
c** Born-Oppenheimer Breakdown & doubling function parameters.
c**                       March 16 2012
c=======================================================================
      INTEGER NUA(NSTATEMX),NUB(NSTATEMX),NTA(NSTATEMX),NTB(NSTATEMX),
     1  IFXUA(0:NBOBMX,NSTATEMX),IFXUB(0:NBOBMX,NSTATEMX),
     2  IFXTA(0:NBOBMX,NSTATEMX),IFXTB(0:NBOBMX,NSTATEMX),
     3  NwCFT(NSTATEMX),IFXwCFT(0:NBOBMX,NSTATEMX),efREF(NSTATEMX)
c
      REAL*8 UA(0:NBOBMX,NSTATEMX),UB(0:NBOBMX,NSTATEMX),
     1  TA(0:NBOBMX,NSTATEMX),TB(0:NBOBMX,NSTATEMX),
     2   wCFT(0:NBOBMX,NSTATEMX)
c
      COMMON /BLKBOB/UA,UB,TA,TB,wCFT,NUA,NUB,NTA,NTB,NwCFT,
     1  IFXUA,IFXUB,IFXTA,IFXTB,IFXwCFT,efREF
c=======================================================================
c-----------------------------------------------------------------------
c** Type statements for input or local variables
      INTEGER I, I1, ISTATE, IISTP, m, MMN, VTST
      CHARACTER*3 SLABL(-6:NSTATEMX)
      REAL*8 ZMASE, RR(NPNTMX), VV(NPNTMX)
      DATA ZMASE /5.4857990945D-04/
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
c               treatment of transitions involving that state]
c         = -1 : represent the rotational sublevels for each v of each
c                isotopologue by Band Constants (!)
c         = 0 : Use a fixed potential defined by LEVEL's PREPOT routine
c         = 1 : Use an Expanded Morse Oscillator EMO(p) potential
c         = 2 : Use a Morse/Long-Range (MLR) Potential.
c         = 3 : Use a Double-Exponential Long-Range (DELR) Potential.
c         = 4 : Use a Surkus Generalized Potential Energy Function (GPEF).
c         = 5 : Use a Tiemann/Hannover-polynomial-potential (HPP)
c         = 6 : Use a Tang-Toennies type potential 
c         = 7 : Use an Aziz'ian HFD-C type potential 
c  MAXMIN(s)= 1  for a regular single-minimum potential, for which finding
c       more than one signals a bad model:  =2 for a double-minimum case
c  VLIM(s)  is the fixed absolute energy of the potential asymptote
c  BOBCN  is a flag to denote reference & scaling for BOB corrections 
c      = 0  using differences as per RJL [JMS 194,189(1999)]
c      = 1  use 'clamped nuclei' limit, m_e1/MASS scaling.
c  OSEL(s)  controls printout of radial function arrays to Ch. 10-16.
c       OSEL > 0: Export to file every OSEL'th point of final function
c=======================================================================
      READ(5,*) PSEL(ISTATE), VLIM(ISTATE), MAXMIN(ISTATE),
     1                                     BOBCN(ISTATE), OSEL(ISTATE)
c=======================================================================
      IF(OSEL(ISTATE).LE.0) OSEL(ISTATE)= 1
      IF((PSEL(ISTATE).EQ.-1).OR.(PSEL(ISTATE).EQ. -2)) THEN
          IF(PSEL(ISTATE).EQ.-2) THEN
c** For term value fits ... no further READs needed!
              WRITE(6,604) SLABL(ISTATE)
              RETURN
              ENDIF
          IF(PSEL(ISTATE).EQ.-1) THEN
c** If representing data for this state by fitted band constants, 
c   read in the number of band constants for each vibrational level
              DO  I= VMIN(ISTATE,1),VMAX(ISTATE,1)
c** For each isotopologue in each vibrational level, read the number of 
c  band constants to be used (fited to) to represent the data,
c=======================================================================
                  READ(5,*) VTST,(NBC(I,IISTP,ISTATE),IISTP= 1,NISTP)
                  IF(IOMEG(ISTATE).GT.0) 
     1                  READ(5,*) (NQC(I,IISTP,ISTATE),IISTP= 1,NISTP)
c=======================================================================
                  IF(I.NE.VTST) THEN
c... Verify that band constant specification is for the correct vib level
                      WRITE(6,610) I,VTST
                      STOP
                      ENDIF
                  DO IISTP= 1,NISTP       !! Check bounds on NBC & NQC
                      IF(NBC(I,IISTP,ISTATE).GT.NBCMX) 
     1                                      NBC(I,IISTP,ISTATE)= NBCMX
                      IF(IOMEG(ISTATE).LE.0) NQC(I,IISTP,ISTATE)= -1
                      IF(NQC(I,IISTP,ISTATE).GT.NBCMX) 
     1                                      NQC(I,IISTP,ISTATE)= NBCMX
                      ENDDO
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
      NDATPT(ISTATE)= MIN(NPNTMX,NDATPT(ISTATE))
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
      IF((PSEL(ISTATE).GE.2).AND.(PSEL(ISTATE).NE.4)) THEN
c-----------------------------------------------------------------------
c** For MLR, DELR and HPP, GTT or HFD potentials, read number of terms NCMM
c    in the {damped} inverse-power long-range tail
c    uLR(R) = - SUM_{i=1}^{NCMM} Dm(R;MMLR(i) * CmVAL(i)/R**MMLR(i) 
c** If rhoAB .LE. 0.0  have NO damping functions: all  Dm(R)= 1.0
c   If rhoAB > 0.0  recommend the molecule-dependent radial scaling 
c                    factor of Douketis et al. [JCP 76, 3057 (1982)]:
c     rhoAB =  2*rhoA*rhoB/(rhoA+rhoB)   where  rhoA  is the ionization
c                 potential ratio  (I_p^A/I_p^H)^0{2/3}  for atom A
c
c    IVSR specifies damping s.th.  Dm(r)/r^m --> r^{IVSR/2} as r->0.
c                    IDSTT > 0  use Douketis et al. damping functions
c                    IDSTT .LE. 0  use Tang-Toennies damping functions
c** IFXCm specifies whether this long-range coefficient is to be fitted
c  freely (when .LE.0), held fixed at the read-in value (when =1) or held
c  fixed at the value for another state, in which case the parameter value
c  'IFXCm(m,ISTATE)' is the no. of the parameter it is constrained to be = to
c
c** For Alkali dimer (nS + nP) states use Aubert-Frecon [PRA 55, 3458 (1997)] 
c  2x2 ULR(r) with NCMM= 7 & MMLR= {x, 3, 3, 6, 6, 8, 8} where x=0 for 
c  the A^1\Sigma_u^+ state and x=-1  for the b^3\Pi_u state, and the
c  read-in C_m's are, DELTAE, C3Sig, C3Pi,C6Sig, C6Pi, C8Sig and C8Pi .
c  FOR the 3x3 cases NCMM=10 and MMLR= {x, 3, 3, 3, 6, 6, 6, 8, 8, 8}
c  where   x= -2 for the c(1^3\Sigma_g^+) state (the lowest 3x3 root),
c  while  CnVAL= {DELTAE, C3Sig, C3Pi1, C3Pi3, C6Sig, C6Pi1, C6Pi3, 
c  C8Sig, C8Pi1, and C8Pi3 .
c  For all Cm's assume units units are cm-1*Angst^m   
c=======================================================================
          READ(5,*) NCMM(ISTATE), rhoAB(ISTATE), IVSR(ISTATE),
     1                                                   IDSTT(ISTATE)
          DO  m= 1,NCMM(ISTATE)
              READ(5,*) MMLR(m,ISTATE), CmVAL(m,ISTATE), IFXCm(m,ISTATE)
              ENDDO
c=======================================================================
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
      IF(PSEL(ISTATE).GE.4) IFXDE(ISTATE)= 1
      IF(PSEL(ISTATE).GE.5) IFXRE(ISTATE)= 1
c=======================================================================
c** Read parameters defining the exponent coefficient function \beta(r)
c* Nbeta(s) is order of the beta(r) exponent polynomial or # spline points
c  APSE(s).LE.0  to use {p,q}-type MLR exponent polynomial of order Nbeta(s)
c     if APSE(s) > 0, \beta(r) is Pashov spline defined by Nbeta(s) points
c  nQB(s)  is the power  q  for  beta(r)  exponent expansion variable
c  nPB(s)  is the power  p  radial and  beta(r)  switching fx. variables
c  RREF(s)  defines the reference distance in the potential exponent 
c    expansion variable:  * for  RREF.le.0 , define parameter  RREF = Re
c      * for  RREF.gt.0 , fix parameter  RREF   at its read-in value
c=======================================================================
      READ(5,*) Nbeta(ISTATE), APSE(ISTATE), nQB(ISTATE), nPB(ISTATE),
     1                                                     RREF(ISTATE)
c=======================================================================
      IF(Nbeta(ISTATE).GE.NbetaMX) THEN
          WRITE(6,648) ISTATE,Nbeta(ISTATE),NbetaMX
          STOP
          ENDIF
      IF((PSEL(ISTATE).EQ.2)) THEN    !! to test if nPB big enuf for MLR
          MMN= MMLR(NCMM(ISTATE),ISTATE)- MMLR(1,ISTATE)
          IF(MMLR(1,ISTATE).le.0)
     1                  MMN= MMLR(NCMM(ISTATE),ISTATE)- MMLR(2,ISTATE)
          IF((NCMM(ISTATE).GT.1).AND.(nPB(ISTATE).LE.MMN)) 
     1                                    WRITE(6,628) nPB(ISTATE),MMN
          ENDIF
      IF((PSEL(ISTATE).EQ.7).AND.((Nbeta(ISTATE).NE.5).AND.
     1                                     (Nbeta(ISTATE).NE.2))) THEN
          WRITE(6,629) Nbeta(ISTATE)
          STOP
          ENDIF
      IF(PSEL(ISTATE).NE.2) APSE(ISTATE)= 0
      IF(APSE(ISTATE).GT.0) THEN
          DO  I= 1, Nbeta(ISTATE)
c-----------------------------------------------------------------------
c** For SE-MLR  exponent is a natural spline function with values BETA
c       at the yq values yqBETA, and fixed to equal yqINF at  yqBETA=1
c=======================================================================
              READ(5,*)yqBETA(I,ISTATE),BETA(I,ISTATE),IFXBETA(I,ISTATE)
c=======================================================================
              ENDDO
          IF(yqBETA(Nbeta(ISTATE),ISTATE).LT.1.d0) THEN
c** Ensure outer endppoint is at yq= 1.d0
              Nbeta(ISTATE)= Nbeta(ISTATE)+1
              yqBETA(Nbeta(ISTATE),ISTATE)= 1.d0
              IFXBETA(Nbeta(ISTATE),ISTATE)= 1
              ENDIF
          ENDIF
c** For non-MLR or the PE-MLR, exponent  \beta(yp)  is  'conventional'
      IF((Nbeta(ISTATE).GE.0).AND.(APSE(ISTATE).LE.0)) THEN
          I1= 0
          IF(PSEL(ISTATE).GE.6) I1= 1   !! omit \beta(0) for TT or HFD
          IF(PSEL(ISTATE).EQ.6) THEN
              Nbeta(ISTATE)= 9
              IDSTT(ISTATE)= 0
              IVSR(ISTATE)= +2
              ENDIF
          IF(PSEL(ISTATE).EQ.5) Nbeta(ISTATE)= Nbeta(ISTATE) + 3
          DO  I= I1, Nbeta(ISTATE)
c** Read in trial initial trial parameters for exponent \beta(r)
c
c  BETA(i,s)   contains the expansion parameters defining the potential
c**  for  PSEL.LE.3 : read-in values are the {Nbeta+1} beta_i  exponent
c                   exponent expansion parameters defining the potential
c**  for  PSEL = 4  : read-in values are leading coefficients in 
c                Surkus' Generalized Potential Energy Function (GPEF).
c**  for  PSEL = 5  : read in the  {1+Nbeta} expansion parameters plus
c                         b, RINN, and ROUT  of the HPP form
c**  for  PSEL = 6,  Nbeta=9  Read in the \\beta_i of Eq.(32) for i=1-9
c** For PSEL=7: >> set Nbeta=4 to use the single global damping function for 
c   HFD-A,B, & C potentials:  f_1(x)= beta(0)*exp{-\beta(1)/x)^beta(2)}
c     while exponent is {\alpha*x + beta(3)*x^2} and  \gamma= beta(4).
c**             >> set Nbeta=3 to combine the overall damping function 
c      f_2(x)= [1 - r^{\beta(0)} exp{-\beta(1)(r}]  , and in this case
c     while exponent is {\alpha*x + beta(2)*x^2} and  \gamma= beta(3).
c  IFXBETA(i,s) indicates whether each potential expansion coefficient
c             coefficient will be:    = 1: held fixed at read-in values.
c                                  .LE. 0: determined from fits.
c=======================================================================
             READ(5,*) BETA(I,ISTATE), IFXBETA(I,ISTATE)
c=======================================================================
             IF(PSEL(ISTATE).GE.5) IFXBETA(I,ISTATE)= 1
             ENDDO
          ENDIF    
c** Note that HPP and (most) HFD potentials assume no damping
      IF((PSEL(ISTATE).EQ.7).and.(Nbeta(ISTATE).EQ.2)) THEN
          IVSR(ISTATE)= 0        !! for HFD-D potentials
          IDSTT(ISTATE)= 1
          ENDIF
      IF((PSEL(ISTATE).EQ.5).OR.((PSEL(ISTATE).EQ.7).and.
     1                    (Nbeta(ISTATE).EQ.4))) rhoAB(ISTATE)= -1.d0
      IF(PSEL(ISTATE).EQ.5) THEN
c** Constraints for Tiemann polynomial potential ....
          nPB(ISTATE)= 1
          nQB(ISTATE)= 1
          IFXDe(ISTATE)= 1
          RREF(ISTATE)= RE(ISTATE)
          ENDIF
c=======================================================================
c** Read parameters defining the BOB adiabatic radial functions
c*  NUA/NUB(s)  specifies the order of the polynomial in  yp  defining 
c              the adiabatic BOB function for atom A/B
c         if < 0   do not read in any adiabatic BOB function parameters
c*  pAD(s)/qAD(s)  are the powers defining the expansion variables
c*  LRad(s) determines whether (if > 0) or not (if .LE.0) isotope shift
c     C_m factors for atoms-A 'dCmA' and atoms B 'dCmB' are to be read in
c*  UA/UB(a,s)   are the adiabatic BOB function expansion coefficients
c*  IFXU(A/B)(a,s) indicates whether each expansion coefficient is to be
c           > 0 :  held fixed at read-in value, or
c        .le. 0 :  varied in the fit
c*  UAinf/UBinf  is the limiting asymptotic value of uA(r)/uB(r), as per
c               Theochem paper [internally stored as  UA(NUA+1), etc.]
c*  IFXUAinf/IFXUBinf  specifies whether (>0) or not (.le.0)  UAinf/UBinf 
c               is to be held fixed at the read-in value
c=======================================================================
      READ(5,*) NUA(ISTATE),NUB(ISTATE),qAD(ISTATE),pAD(ISTATE),
     1                                                    LRad(ISTATE)
      IF(((NUA(ISTATE).GE.0).OR.(NUB(ISTATE).GE.0))
     1              .AND.(PSEL(ISTATE).EQ.1)) pAD(ISTATE)= qAD(ISTATE)
c... NOTE never read delta Cm values unless PSEL = 1-3
          IF((PSEL(ISTATE).LT.1).OR.(PSEL(ISTATE).GT.3)) LRad(ISTATE)=0
          IF(LRad(ISTATE).GT.0) THEN
c. if desired, read \delta{Cm} values dCmA & dCmB for atoms-A & B one per line
              DO m=1,NCMM(ISTATE)
                  READ(5,*) dCmA(m,ISTATE)
                  ENDDO
              DO m=1,NCMM(ISTATE)
                  READ(5,*) dCmB(m,ISTATE)
                  ENDDO
          ENDIF
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
c  qNA(s)  is  the power defining the the form of the expansion variable
c  TA/TB(a,s)   are the non-adiabatic centrifugal BOB expansion coeffts
c  IFXTA/IFXTB(a,s) indicates whether each expansion coefficient is to be
c           > 0 :  held fixed at read-in value, or
c        .le. 0 :  varied in the fit
c  TAinf/TBinf  is the limiting asymptotic value of qA(r)/qB(r), as per 
c               Theochem paper [internally stored as  TA(NTA+1), etc.]
c  IFXTAinf/IFXTBinf  specifies whether (>0) or not (.le.0)  TAinf/TBinf 
c               is to be held fixed at the read-in value
c=======================================================================
      READ(5,*) NTA(ISTATE), NTB(ISTATE), qNA(ISTATE), pNA(ISTATE)
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
c*  Pqw(s)  defined nature of radial expansion variable:  
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
c   non-adiabatic centrifugal (ZMTA, ZMTB) BOB functions
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
  999 RETURN
  600 FORMAT(/'For state ',A3,' use a fixed potential defined by LEVEL s
     1ubroutine PREPOT'/4x,'Integrate from   RMIN=',f5.2,
     1  '   to   RMAX=',f6.2,'   with mesh   RH=',f8.5)
  604 FORMAT(/' For state ',A3,' represent level energies by independent
     1 term values')
  610 FORMAT(' *** Input ERROR *** band constant specification  v=',I3,
     1  ' .NE.', I3)
  626 FORMAT(/' For state   ',A3/4x,'integrate from   RMIN=',f5.2,
     1  '   to   RMAX=',f6.2,'   with mesh   RH=',f8.5)
  628 FORMAT(' ***** WARNING  p=',i2,' .LE.[MMLR(NCMM)-MMLR(1)]=',i2,
     1  ' *****'/"   so tail of MLR exponential will 'pollute' u_{LR}(r)
     2 behaviour"/(2x,19('****')))
  629 FORMAT(' *** ERROR *** For HFD potentials  Nbeta=',I3,'  should be
     1  5  or  2 !!')
  640 FORMAT(/' ', A3,' state energies referenced to f-parity levels')
  642 FORMAT(/' ', A3,' state energies referenced to the mid-point betwe
     1en e and f-parity levels')
  644 FORMAT(/' ', A3,' state energies referenced to e-parity levels')
  646 FORMAT(/' *** INPUT ERROR ***  |efREF=',i3,'| > 1')
  648 FORMAT(/'  For ISTATE=',I2,'  read-in  Nbeta=',I3,'  while NbetaMX
     1=',I3,'  so STOP!!' )
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

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
cc    INCLUDE 'arrsizes.h'
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c** 'Block' Data Utility routine named: 'arrsizes.h' that governs 
c    array dimensioning in program  dPotFit
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER NISTPMX,NPARMX,NbetaMX,NBOBMX,HPARMX,NDATAMX,
     1  NVIBMX,NBCMX,NSTATEMX,NPNTMX,NROTMX,NCMMAX
c*  NISTPMX  is the maximum number of isotopomers allowed for fit
      PARAMETER (NISTPMX =  12)
c*  NSTATEMX  is maximum no. of electronic states which can be
c             simultaneously fitted to
      PARAMETER (NSTATEMX = 4)
c*  NPARMX  is the largest number of free parameters allowed for fit
c  Since FS origins may be parameters, this is also max. no, data bands
      PARAMETER (NPARMX  = 8000)
c*  NbetaMX  is the largest number of exponent parameters allowed for fit
      PARAMETER (NbetaMX  = 40)
c*  NBOBMX-1  is the highest-order polynomial expansion allowed for the
c      adiabatic or centrifugal Born-Oppenheimer breakdown functions, or 
c      the Lambda-doubling or 2\Sigma splitting radial strength functions
      PARAMETER (NBOBMX  = 15)
c*  HPARMX  is the largest number of Hamiltonian parameters of all types
c    (potential energy, BOB. etc.) for all states.
c           HPARMX >= NSTATEMX*[5 + (NbetaMX+1) + 5*(NBOBMX+1)]
      PARAMETER (HPARMX= NSTATEMX*(5 + (NbetaMX+1) + 5*(NBOBMX+1)))
cc    PARAMETER (HPARMX = 300)
c*  NDATAMX  is largest No. of individual data which may be considered
      PARAMETER (NDATAMX = 35000)
c*  NVIBMX  is the maximum number of vibrational levels of a single
c           state for which data are to be considered
      PARAMETER (NVIBMX    = 200)
** NBCMX  is the maximum number of band constants per vib level to be
c         allowed when doing band constant fits (PSEL= -1)
      PARAMETER (NBCMX = 8) 
c*  NPNTMX  is the largest number of potential data points that can be
c           stored in a single 1D array
      PARAMETER (NPNTMX = 90000)
c*  NROTMX  is the highest order of rotational constants calculated and
c            used for estimating level energies
      PARAMETER (NROTMX = 7)
c*  NCMMAX is the largest number of Cm terms in the MLR or DELR 
c            long-range potential
      PARAMETER (NCMMAX = 12)
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cc    INCLUDE 'BLKPOT.h'
c=======================================================================
c** Effective adiabatic radial potential variables.
      INTEGER BOBCN(NSTATEMX),PSEL(NSTATEMX),MAXMIN(NSTATEMX),
     1 IOMEG(NSTATEMX),Nbeta(NSTATEMX),APSE(NSTATEMX),IFXDE(NSTATEMX),
     2 IFXRE(NSTATEMX),IFXCm(NCMMax,NSTATEMX),
     3 IFXBETA(0:NbetaMX,NSTATEMX),NDATPT(NSTATEMX),NCMM(NSTATEMX),
     4 MMLR(NCMMax,NSTATEMX),nPB(NSTATEMX),nQB(NSTATEMX),pAD(NSTATEMX),
     5 qAD(NSTATEMX),LRad(NSTATEMX),pNA(NSTATEMX),qNA(NSTATEMX),
     6 Pqw(NSTATEMX),IVSR(NSTATEMX),IDSTT(NSTATEMX)
c
      REAL*8 DE(NSTATEMX),RE(NSTATEMX),BETA(0:NbetaMX,NSTATEMX),
     1 yqBETA(NbetaMX,NSTATEMX),BETAFX(NPNTMX,NSTATEMX),RH(NSTATEMX),
     2 RMIN(NSTATEMX),RMAX(NSTATEMX),VLIM(NSTATEMX),EPS(NSTATEMX),
     3 betaINF(NSTATEMX),AGPEF(NSTATEMX),BGPEF(NSTATEMX),
     4 CmVAL(NCMMax,NSTATEMX),CmEFF(NCMMax,NSTATEMX),rhoAB(NSTATEMX),
     5 AA(NSTATEMX),BB(NSTATEMX),RREF(NSTATEMX),ASO(NSTATEMX),
     6 R01(NSTATEMX),Q12(NSTATEMX),RD(NPNTMX,NSTATEMX),
     7 VPOT(NPNTMX,NSTATEMX),dCmA(NCMMax,NSTATEMX),dCmB(NCMMax,NSTATEMX)
c
      COMMON /BLKPOT/DE,RE,BETA,yqBETA,BETAFX,RH,RMIN,RMAX,VLIM,EPS,
     1 betaINF,AGPEF,BGPEF,CmVAL,CmEFF,rhoAB,AA,BB,RREF,ASO,R01,Q12,RD,
     2 VPOT,dCmA,dCmB, BOBCN,PSEL,MAXMIN,IOMEG,Nbeta,APSE,IFXDE,IFXRE,
     3 IFXCm,IFXBETA,NDATPT,NCMM,MMLR,nPB,nQB,pAD,qAD,LRad,pNA,qNA,Pqw,
     4 IVSR,IDSTT
c=======================================================================
cc    INCLUDE 'BLKPARAM.h'
c=======================================================================
c** Parameters and count-labels for band constant (PSEL=-1) or term
c   value (PSEL=-2) fits
      REAL*8 TVALUE(NPARMX),ZBC(0:NVIBMX,0:NROTMX,NISTPMX,NSTATEMX),
     1 ZQC(0:NVIBMX,0:NROTMX,NISTPMX,NSTATEMX)
c
      INTEGER NSTATES,NTVALL(0:NSTATEMX),NTVI(NSTATEMX),NTVF(NSTATEMX),
     1 VMIN(NSTATEMX,NISTPMX),VMAX(NSTATEMX,NISTPMX),JTRUNC(NSTATEMX),
     2 EFSEL(NSTATEMX),NBC(0:NVIBMX,NISTPMX,NSTATEMX),
     3 NQC(0:NVIBMX,NISTPMX,NSTATEMX),
     4 BCPARI(0:NVIBMX,NISTPMX,NSTATEMX),
     5 BCPARF(0:NVIBMX,NISTPMX,NSTATEMX),
     6 QCPARI(0:NVIBMX,NISTPMX,NSTATEMX),
     7 QCPARF(0:NVIBMX,NISTPMX,NSTATEMX)
      COMMON /BLKPARAM/TVALUE,ZBC,ZQC,NSTATES,NTVALL,NTVI,NTVF,VMIN,
     1      VMAX,JTRUNC,EFSEL,NBC,NQC,BCPARI,BCPARF,QCPARI,QCPARF
c=======================================================================
cc    INCLUDE 'BLKISOT.h'
c=======================================================================
c** Isotope/isotopologue numbers, masses & BOB mass scaling factors
c** Array ZK carries about the band constants for all levels of all ISOT
      INTEGER NISTP,NDUNMX,AN(2),MN(2,NISTPMX)
c** NDUNMX is a dummy parameter reqd. for portability of READATA
      PARAMETER (NDUNMX=0)
      REAL*8  ZMASS(3,NISTPMX),RSQMU(NISTPMX),RSQMUP(0:NDUNMX,NISTPMX),
     1 RMUP(0:9,NISTPMX),ZMUA(NISTPMX,NSTATEMX),ZMUB(NISTPMX,NSTATEMX),
     2 ZMTA(NISTPMX,NSTATEMX),ZMTB(NISTPMX,NSTATEMX),
     3  ZK(0:NVIBMX,0:NROTMX,NISTPMX,NSTATEMX)
c
      COMMON /BLKISOT/ZMASS,RSQMU,RSQMUP,RMUP,ZMUA,ZMUB,ZMTA,ZMTB,ZK,
     1  NISTP,AN,MN
c=======================================================================
cc    INCLUDE 'BLKBOB.h'
c=======================================================================
c** Born-Oppenheimer Breakdown & doubling function parameters.
c**                       March 16 2012
c=======================================================================
      INTEGER NUA(NSTATEMX),NUB(NSTATEMX),NTA(NSTATEMX),NTB(NSTATEMX),
     1  IFXUA(0:NBOBMX,NSTATEMX),IFXUB(0:NBOBMX,NSTATEMX),
     2  IFXTA(0:NBOBMX,NSTATEMX),IFXTB(0:NBOBMX,NSTATEMX),
     3  NwCFT(NSTATEMX),IFXwCFT(0:NBOBMX,NSTATEMX),efREF(NSTATEMX)
c
      REAL*8 UA(0:NBOBMX,NSTATEMX),UB(0:NBOBMX,NSTATEMX),
     1  TA(0:NBOBMX,NSTATEMX),TB(0:NBOBMX,NSTATEMX),
     2   wCFT(0:NBOBMX,NSTATEMX)
c
      COMMON /BLKBOB/UA,UB,TA,TB,wCFT,NUA,NUB,NTA,NTB,NwCFT,
     1  IFXUA,IFXUB,IFXTA,IFXTB,IFXwCFT,efREF
c=======================================================================
cc    INCLUDE 'BLKCOUNT.h'
c=======================================================================
c    Block data file  BLKCOUNT.h
c=======================================================================
c** Counters for numbers of potential parameters of different types for 
c   each state
      INTEGER  TOTPOTPAR,POTPARI(NSTATEMX),POTPARF(NSTATEMX),
     1  UAPARI(NSTATEMX),UAPARF(NSTATEMX),UBPARI(NSTATEMX),
     2  UBPARF(NSTATEMX),TAPARI(NSTATEMX),TAPARF(NSTATEMX),
     3  TBPARI(NSTATEMX),TBPARF(NSTATEMX),LDPARI(NSTATEMX),
     4  LDPARF(NSTATEMX),HPARF(NSTATEMX),OSEL(NSTATEMX)
c
      COMMON /BLKCOUNT/TOTPOTPAR,POTPARI,POTPARF,UAPARI,UAPARF,UBPARI,
     1  UBPARF,TAPARI,TAPARF,TBPARI,TBPARF,LDPARI,LDPARF,HPARF,OSEL
c=======================================================================
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
cc     DATA cDS/0.468d0,0.446d0,0.423d0,0.405d0,0.390d0,0.d0,0.360d0,
       DATA cDS/0.468d0,0.446d0,0.423d0,0.400d0,0.390d0,0.d0,0.360d0,
     1            0.d0,0.340d0/
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
          AREF= RREF(ISTATE)
          IF(AREF.LE.0.d0) AREF= RE(ISTATE)
          IF(PSEL(ISTATE).EQ.1) THEN
c** Header printout for EMO potential
              WRITE(6,602) SLABL(ISTATE),Nbeta(ISTATE),nQB(ISTATE),
     1                                       nQB(ISTATE),Nbeta(ISTATE)
              IF(RREF(ISTATE).LE.0.d0) WRITE(6,552) (nQB(ISTATE),i=1,5)
              IF(RREF(ISTATE).GT.0.d0) WRITE(6,555) AREF,AREF
              ENDIF
          IF(PSEL(ISTATE).EQ.2) THEN
c** Header printout for MLR potential
              BINF= betaINF(ISTATE)
              WRITE(6,604) SLABL(ISTATE),nQB(ISTATE),nPB(ISTATE) 
              IF(APSE(ISTATE).LE.0) THEN
                  WRITE(6,605) nPB(ISTATE),nPB(ISTATE),nQB(ISTATE),
     1                                                   Nbeta(ISTATE)
                ELSE
                  WRITE(6,680) Nbeta(ISTATE)
                  BETA(Nbeta(ISTATE),ISTATE)= BINF
                ENDIF
              IF(RREF(ISTATE).LE.0.d0) WRITE(6,552) (nQB(ISTATE),i=1,5)
              IF(RREF(ISTATE).GT.0.d0) WRITE(6,555) AREF,AREF
              ENDIF
c
          IF(PSEL(ISTATE).EQ.3) THEN
c** Header printout for DELR potential form ...
              WRITE(6,612) SLABL(ISTATE),Nbeta(ISTATE),nQB(ISTATE),
     1                                       nQB(ISTATE),Nbeta(ISTATE)
              IF(RREF(ISTATE).LE.0.d0) WRITE(6,552) (nQB(ISTATE),i=1,5)
              IF(RREF(ISTATE).GT.0.d0) WRITE(6,555) AREF,AREF
              ENDIF
c
          IF(PSEL(ISTATE).EQ.4) THEN
c** Header printout for Surkus GPEF potential form ...
              WRITE(6,610) SLABL(ISTATE),(nPB(ISTATE),i=1,3),
     1             AGPEF(ISTATE),nPB(ISTATE),BGPEF(ISTATE),nPB(ISTATE)
              ENDIF
c
          IF(PSEL(ISTATE).EQ.5) THEN
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
     1s   gamma=',f9.6/'   beta1=',f12.8,'   beta2=',f9.6,5x,'A=',
     2 1PD16.9:"  &  Cm's:"/(3x,3('   C',I2,' =',D15.7:)))
  629 FORMAT('   De=',f10.4,'[cm-1]   Re=',f9.6,'[Angst.]':'   and   for
     1  r <',F9.6/'     Damping function  D(r)= exp[ -',f6.4,'*(',f9.6, 
     2  '/r -1.0)**',f5.2,']')
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
cc                WRITE(6,608) MMLR(1,ISTATE),CmVAL(1,ISTATE), 
cc                                 1 MMLR(1,ISTATE)
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
          IF((IOMEG(ISTATE).NE.0).AND.(PSEL(ISTATE).GE.0))
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
c** Writing out the De information.
c-----------------------------------------------------------------------
          IF((PSEL(ISTATE).GE.1).AND.(PSEL(ISTATE).NE.4)) THEN
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
          IF((PSEL(ISTATE).GE.2).AND.(PSEL(ISTATE).NE.4)) THEN
c-----------------------------------------------------------------------
c** For MLR or DELR, HPP, TT, or HFD write out the  Cm  information.
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
              yqRe=(RE(ISTATE)**nQB(ISTATE) - AREF**nQB(ISTATE))
     1                  /(RE(ISTATE)**nQB(ISTATE) + AREF**nQB(ISTATE))
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
                  CALL AFdiag(RE(ISTATE),NCMM(ISTATE),
     1             NCMMax,MMLR1D,PSEL(ISTATE),DECM(ISTATE),
     2             CmVAL(1:NCMMax,ISTATE),rhoAB(ISTATE),IVSR(ISTATE),
     3             IDSTT(ISTATE),VATTRe,dULRdCm,dVATTRE,dULRdDe)
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
          IF(NPASS.GT.1) WRITE(20,671) Nbeta(ISTATE), APSE(ISTATE),
     1                          nQB(ISTATE), nPB(ISTATE), RREF(ISTATE)
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
  671 FORMAT(/2I3,I4,I3,1PD11.2,8x,'% Nbeta APSE nQB nPB RREF') 
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
  564 FORMAT(4x,A10,A2,';R) = \sum\{u_i * [y',I1,']^i}   for  i= 0 to',
     1 i3)
  558 FORMAT(4x,A10,A2,';R) = u(inf)*y',i1,' + (1 - y',I1, ')*\Sum{u_i *
     1 [y',I1,']^i}   for  i= 0 to',i3)
  559 FORMAT(' Non-Adiabatic centrifugal BOB fx. with {p=',i2,', q=',i2,
     1  '} are scaled by  M(1)/M(A):')
  560 FORMAT(' Non-Adiabatic centrifugal BOB fx are power series in y_',
     1  I2,'(r) scaled by  m_e/M(A):')
  570 FORMAT(8x,'but replace  u(inf)  with  u(inf) + Sum{dCm/r**m}  wher
     1e:')
  571 FORMAT(48x,'dC',I1,'=',1PD14.7,'[cm-1 Ang^',i1,']') 
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
     1gy Function GPEF with"/1x,6('=='),'  expansion vble:   y_',i1, 
     2 '(r) = (r^',i1,' - re^',i1,')/(',F5.2,'*r^',i1,' +',F5.2,'*re^',
     3 i1,'p)')
 6102 FORMAT(' *** Input ERROR *** band constant specification  v=',I3,
     1  ' .NE.', I3)
  612 FORMAT(/' State ',A3,' represented by a DELR(N=',I2,'  q=',i2,')
     1potential with'/1x,4('=='),'   exponent coefficient:   beta(r)=
     Su 2m{beta_i*y',i1,'^i}'/48x,'with polynomial order   Nbeta=',I3)
  607 FORMAT(4x,'uLR inverse-power terms incorporate DS-type damping wit
     1h   rhoAB=',f10.7/11x,'defined to give very short-range  Dm(r)*Cm/
     2r^m  behaviour   r^{',SP,f4.1,'}'/8x,SS,'Dm(r)= [1 - exp(-',f5.2, 
     3 '(rhoAB*r)/m -',f6.3,'(rhoAB*r)^2/sqrt{m})]^{m',SP,F4.1,'}')
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
  614 FORMAT('   Parameter    Initial Value    Uncertainty   Sensitivity
     1')
  615 FORMAT('   Parameter      Final Value    Uncertainty   Sensitivity
     1')
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
  703 FORMAT(I3,F7.2,2I3,15x,' % NCMM RHOab IDF IDSTT') 
  704 FORMAT(I3,1P,D17.8,0P,I3,8x,' % MMLR(',I1,'), CmVAL(',I1,
     1  '), IFXCm(',I1,')')
  705 FORMAT(10I4)
  706 FORMAT(I3,1P,D17.8,0P,I3,8x,' % ',A10)
  708 FORMAT(39x,A10,'=',1PD14.7,'[cm-1 Ang^{',i2,'}]')
  720 FORMAT(2X,A6,3X,1PD21.12,3X,1PD8.1,6X,1PD8.1,10x)
  721 FORMAT(1X,'*',A6,3X,1PD21.12,3X,1PD8.1,6X,1PD8.1) 
  722 FORMAT(2X,A6,3X,1PD21.12,7X,'--',12X,'--')
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

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
c               ----- Version of  17 March 2016 -----
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
cc    INCLUDE 'arrsizes.h'
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c** 'Block' Data Utility routine named: 'arrsizes.h' that governs 
c    array dimensioning in program  dPotFit
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER NISTPMX,NPARMX,NbetaMX,NBOBMX,HPARMX,NDATAMX,
     1  NVIBMX,NBCMX,NSTATEMX,NPNTMX,NROTMX,NCMMAX
c*  NISTPMX  is the maximum number of isotopomers allowed for fit
      PARAMETER (NISTPMX =  12)
c*  NSTATEMX  is maximum no. of electronic states which can be
c             simultaneously fitted to
      PARAMETER (NSTATEMX = 4)
c*  NPARMX  is the largest number of free parameters allowed for fit
c  Since FS origins may be parameters, this is also max. no, data bands
      PARAMETER (NPARMX  = 8000)
c*  NbetaMX  is the largest number of exponent parameters allowed for fit
      PARAMETER (NbetaMX  = 40)
c*  NBOBMX-1  is the highest-order polynomial expansion allowed for the
c      adiabatic or centrifugal Born-Oppenheimer breakdown functions, or 
c      the Lambda-doubling or 2\Sigma splitting radial strength functions
      PARAMETER (NBOBMX  = 15)
c*  HPARMX  is the largest number of Hamiltonian parameters of all types
c    (potential energy, BOB. etc.) for all states.
c           HPARMX >= NSTATEMX*[5 + (NbetaMX+1) + 5*(NBOBMX+1)]
      PARAMETER (HPARMX= NSTATEMX*(5 + (NbetaMX+1) + 5*(NBOBMX+1)))
cc    PARAMETER (HPARMX = 300)
c*  NDATAMX  is largest No. of individual data which may be considered
      PARAMETER (NDATAMX = 35000)
c*  NVIBMX  is the maximum number of vibrational levels of a single
c           state for which data are to be considered
      PARAMETER (NVIBMX    = 200)
** NBCMX  is the maximum number of band constants per vib level to be
c         allowed when doing band constant fits (PSEL= -1)
      PARAMETER (NBCMX = 8) 
c*  NPNTMX  is the largest number of potential data points that can be
c           stored in a single 1D array
      PARAMETER (NPNTMX = 90000)
c*  NROTMX  is the highest order of rotational constants calculated and
c            used for estimating level energies
      PARAMETER (NROTMX = 7)
c*  NCMMAX is the largest number of Cm terms in the MLR or DELR 
c            long-range potential
      PARAMETER (NCMMAX = 12)
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cc    INCLUDE 'BLKISOT.h'
c=======================================================================
c** Isotope/isotopologue numbers, masses & BOB mass scaling factors
c** Array ZK carries about the band constants for all levels of all ISOT
      INTEGER NISTP,NDUNMX,AN(2),MN(2,NISTPMX)
c** NDUNMX is a dummy parameter reqd. for portability of READATA
      PARAMETER (NDUNMX=0)
      REAL*8  ZMASS(3,NISTPMX),RSQMU(NISTPMX),RSQMUP(0:NDUNMX,NISTPMX),
     1 RMUP(0:9,NISTPMX),ZMUA(NISTPMX,NSTATEMX),ZMUB(NISTPMX,NSTATEMX),
     2 ZMTA(NISTPMX,NSTATEMX),ZMTB(NISTPMX,NSTATEMX),
     3  ZK(0:NVIBMX,0:NROTMX,NISTPMX,NSTATEMX)
c
      COMMON /BLKISOT/ZMASS,RSQMU,RSQMUP,RMUP,ZMUA,ZMUB,ZMTA,ZMTB,ZK,
     1  NISTP,AN,MN
c=======================================================================
cc    INCLUDE 'BLKDATA.h'
c=======================================================================
c** Type statements & common block for data
      REAL*8  FREQ(NDATAMX),UFREQ(NDATAMX),DFREQ(NDATAMX),TEMP(NDATAMX),
     1                                               YUNC(NDATAMX),Fqb
      INTEGER  COUNTOT,NFS1,NFSTOT,NBANDTOT,IB(NDATAMX),JP(NDATAMX),
     1 JPP(NDATAMX),VP(NPARMX),VPP(NPARMX),EFP(NDATAMX),EFPP(NDATAMX),
     2 TVUP(NDATAMX),TVLW(NDATAMX),FSBAND(NPARMX),IFXFS(NPARMX),
     3 NFS(NPARMX),IEP(NPARMX),IEPP(NPARMX),ISTP(NPARMX),
     4 IFIRST(NPARMX),ILAST(NPARMX),NTV(NSTATEMX,NISTPMX),FSsame,
     5 NTRANS(NPARMX),IBB(NISTPMX,NSTATEMX,9,NPARMX),JMIN(NPARMX),
     6 JMAX(NPARMX)
      CHARACTER*2 NAME(2) 
      CHARACTER*3 SLABL(-6:NSTATEMX)
      CHARACTER*30 BANDNAME(NPARMX)
      COMMON /DATABLK/Fqb,FREQ,UFREQ,YUNC,DFREQ,TEMP,COUNTOT,NFS1,
     1 NFSTOT,NBANDTOT,IB,JP,JPP,VP,VPP,EFP,EFPP,TVUP,TVLW,FSBAND,IFXFS,
     2 NFS,IEP,IEPP,ISTP,IFIRST,ILAST,NTV,FSsame,
     3 NTRANS,IBB,JMIN,JMAX,NAME,SLABL,BANDNAME
c=======================================================================
cc    INCLUDE 'BLKPOT.h'
c=======================================================================
c** Effective adiabatic radial potential variables.
      INTEGER BOBCN(NSTATEMX),PSEL(NSTATEMX),MAXMIN(NSTATEMX),
     1 IOMEG(NSTATEMX),Nbeta(NSTATEMX),APSE(NSTATEMX),IFXDE(NSTATEMX),
     2 IFXRE(NSTATEMX),IFXCm(NCMMax,NSTATEMX),
     3 IFXBETA(0:NbetaMX,NSTATEMX),NDATPT(NSTATEMX),NCMM(NSTATEMX),
     4 MMLR(NCMMax,NSTATEMX),nPB(NSTATEMX),nQB(NSTATEMX),pAD(NSTATEMX),
     5 qAD(NSTATEMX),LRad(NSTATEMX),pNA(NSTATEMX),qNA(NSTATEMX),
     6 Pqw(NSTATEMX),IVSR(NSTATEMX),IDSTT(NSTATEMX)
c
      REAL*8 DE(NSTATEMX),RE(NSTATEMX),BETA(0:NbetaMX,NSTATEMX),
     1 yqBETA(NbetaMX,NSTATEMX),BETAFX(NPNTMX,NSTATEMX),RH(NSTATEMX),
     2 RMIN(NSTATEMX),RMAX(NSTATEMX),VLIM(NSTATEMX),EPS(NSTATEMX),
     3 betaINF(NSTATEMX),AGPEF(NSTATEMX),BGPEF(NSTATEMX),
     4 CmVAL(NCMMax,NSTATEMX),CmEFF(NCMMax,NSTATEMX),rhoAB(NSTATEMX),
     5 AA(NSTATEMX),BB(NSTATEMX),RREF(NSTATEMX),ASO(NSTATEMX),
     6 R01(NSTATEMX),Q12(NSTATEMX),RD(NPNTMX,NSTATEMX),
     7 VPOT(NPNTMX,NSTATEMX),dCmA(NCMMax,NSTATEMX),dCmB(NCMMax,NSTATEMX)
c
      COMMON /BLKPOT/DE,RE,BETA,yqBETA,BETAFX,RH,RMIN,RMAX,VLIM,EPS,
     1 betaINF,AGPEF,BGPEF,CmVAL,CmEFF,rhoAB,AA,BB,RREF,ASO,R01,Q12,RD,
     2 VPOT,dCmA,dCmB, BOBCN,PSEL,MAXMIN,IOMEG,Nbeta,APSE,IFXDE,IFXRE,
     3 IFXCm,IFXBETA,NDATPT,NCMM,MMLR,nPB,nQB,pAD,qAD,LRad,pNA,qNA,Pqw,
     4 IVSR,IDSTT
c=======================================================================
cc    INCLUDE 'BLKDVDP.h'
c=======================================================================
c** Partial derivative arrays for fits and uncertainties (fununc)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      REAL*8 DVtot(HPARMX,NPNTMX),DLDDRe(NPNTMX,NSTATEMX),
     1  DUADRe(NPNTMX,NSTATEMX),DUBDRe(NPNTMX,NSTATEMX),
     2  DTADRe(NPNTMX,NSTATEMX),DTBDRe(NPNTMX,NSTATEMX),
     3  DBDB(0:NbetaMX,NPNTMX,NSTATEMX),DBDRe(NPNTMX,NSTATEMX),
     4  dVpdP(HPARMX,NPNTMX)
      COMMON/BLKDVDP/DVtot,DUADRe,DUBDRe,DTADRe,DTBDRe,DLDDRe,DBDB,
     1 DBDRe,dVpdP
c=======================================================================
cc    INCLUDE 'BLKBOB.h'
c=======================================================================
c** Born-Oppenheimer Breakdown & doubling function parameters.
c**                       March 16 2012
c=======================================================================
      INTEGER NUA(NSTATEMX),NUB(NSTATEMX),NTA(NSTATEMX),NTB(NSTATEMX),
     1  IFXUA(0:NBOBMX,NSTATEMX),IFXUB(0:NBOBMX,NSTATEMX),
     2  IFXTA(0:NBOBMX,NSTATEMX),IFXTB(0:NBOBMX,NSTATEMX),
     3  NwCFT(NSTATEMX),IFXwCFT(0:NBOBMX,NSTATEMX),efREF(NSTATEMX)
c
      REAL*8 UA(0:NBOBMX,NSTATEMX),UB(0:NBOBMX,NSTATEMX),
     1  TA(0:NBOBMX,NSTATEMX),TB(0:NBOBMX,NSTATEMX),
     2   wCFT(0:NBOBMX,NSTATEMX)
c
      COMMON /BLKBOB/UA,UB,TA,TB,wCFT,NUA,NUB,NTA,NTB,NwCFT,
     1  IFXUA,IFXUB,IFXTA,IFXTB,IFXwCFT,efREF
c=======================================================================
cc    INCLUDE 'BLKBOBRF.h'
c=======================================================================
c** Born-Oppenheimer breakdown radial functions 
      REAL*8 UAR(NPNTMX,NSTATEMX),UBR(NPNTMX,NSTATEMX),
     1 TAR(NPNTMX,NSTATEMX),TBR(NPNTMX,NSTATEMX),wRAD(NPNTMX,NSTATEMX)
c
      COMMON /BLKBOBRF/UAR,UBR,TAR,TBR,wRAD
c=======================================================================
cc    INCLUDE 'BLKCOUNT.h'
c=======================================================================
c    Block data file  BLKCOUNT.h
c=======================================================================
c** Counters for numbers of potential parameters of different types for 
c   each state
      INTEGER  TOTPOTPAR,POTPARI(NSTATEMX),POTPARF(NSTATEMX),
     1  UAPARI(NSTATEMX),UAPARF(NSTATEMX),UBPARI(NSTATEMX),
     2  UBPARF(NSTATEMX),TAPARI(NSTATEMX),TAPARF(NSTATEMX),
     3  TBPARI(NSTATEMX),TBPARF(NSTATEMX),LDPARI(NSTATEMX),
     4  LDPARF(NSTATEMX),HPARF(NSTATEMX),OSEL(NSTATEMX)
c
      COMMON /BLKCOUNT/TOTPOTPAR,POTPARI,POTPARF,UAPARI,UAPARF,UBPARI,
     1  UBPARF,TAPARI,TAPARF,TBPARI,TBPARF,LDPARI,LDPARF,HPARF,OSEL
c=======================================================================
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
     3 YPA,YPB,YQA,YQB,YPE,YPM,YPMA,YPMB,YPP,YQP,REp,RDp,RDq,DYPDRE,
     4 DYQDRE,VAL,DVAL,HReP,HReQ,yqRe,dyqRedRe,betaRe,DbetaRe,yPOW,
     4 dAAdb0,dbetaFX,ULRe,dULRe,d2ULRe,SL,SLB,SLBB,AREF,AREFp,AREFq,T0,
     5 T1,T2,Scalc,dLULRedRe,dLULRedCm(NCMMax),dLULRedDe,dULRdDe,rhoINT,
     6 dULRdCm(NCMMax),dULRepdCm(NCMMax),dULRedCm(NCMMax),DVDD,RDIST,
     7 VDIST,BETADIST,X,BFCT,JFCT,JFCTLD,REadAp,REadBp,REadAq,REadBq,
     8 REnaAp,REnaBp,REnaAq,REnaBq,REwp,dC6dDe,dC9dC3,dC9dC6,dC9dDe,BT,
     9 Rinn,Rout,A0,A1,A2,A3,xBETA(NbetaMX),rKL(NbetaMX,NbetaMX),C1LIM,
     o B5,BETA0,BETAN,TM,VATT,dATTdRe,dATTdb,ATT,REq,VMIN,
     a REQQ,XRI,dXRI,fRO,XRIpw,XRO,dXRO,XROpw,ROmp2,dXROdRe,d2XROdRe,
     b DXRIdRe,d2XRIdRe,dCmp2dRe,EXPBI,BIrat,CMMp2,RMMp2,dAIdRe,dBIdRe,
     c VX,dVX,dVdRe,dDeROdRe,dDeRIdRe,dULRdR,dCmASUM, dCmBSUM,
     d dAI(0:NbetaMX),dBI(0:NbetaMX),dCmp2(0:NbetaMX),
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
      REp= RE(ISTATE)**nPB(ISTATE)
      REq= RE(ISTATE)**nQB(ISTATE)
      IF(RREF(ISTATE).LE.0) AREF= RE(ISTATE)
      IF(RREF(ISTATE).GT.0) AREF= RREF(ISTATE)
      AREFP= AREF**nPB(ISTATE)
      AREFQ= AREF**nQB(ISTATE)
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
c** First ... calculate the Extended Morse Oscillator exponent
          DO  I= ISTART,ISTOP
              RVAL= RD(I,ISTATE)
              IF(RDIST.GT.0.d0) RVAL= RDIST
              RDQ= RVAL**nQB(ISTATE)
              YQ= (RDQ - AREFQ)/(RDQ + AREFQ)
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
              IF(RREF(ISTATE).LE.0.d0) THEN
                  DBDRe(I,ISTATE)= -0.5d0*nPB(ISTATE)*(1.d0-YP**2)
     1                                                *DVAL/RE(ISTATE)
                  VAL= VAL - (RVAL- RE(ISTATE))*DBDRe(I,ISTATE) 
                  ENDIF
              IPV= IPV+1
              DVtot(IPV,I)= -YPP*VAL       !! derivative w.r.t Re 
              YQP= YPP*(RVAL - RE(ISTATE)) 
              DO  J= 0, Nbeta(ISTATE)
                  IPV= IPV+1
                  DVtot(IPV,I)= YQP        !! derivative w.r.t. \beta_i 
                  YQP= YQP*YQ 
                  ENDDO
              ENDDO
ccccc Print for testing
          rewind(10) 
          write(10,610) (RD(i,ISTATE),vpot(i,istate),BETAFX(i,istate),
     1                              i= 1, NDATPT(ISTATE),OSEL(ISTATE))
ccccc
          ENDIF
c........ End preparation of Expanded Morse Potential Function .........
      IF(PSEL(ISTATE).EQ.2) THEN
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** For the  {Morse/Long-Range}  potential.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** For 'normal' inverse-power sum MLR case, with or without damping,
c   set up and write 'Dattani-corrected' effective Cm values 
          DO m= 1, NCMM(ISTATE)
              CmEFF(m,ISTATE)= CmVAL(m,ISTATE)
              Cm1D(m)= CmVAL(m,ISTATE)
              MMLR1D(m)= MMLR(m,ISTATE)   !! powers in 0D array for dampF
              ENDDO
          MCMM= NCMM(ISTATE)
cc   ! Assuming {adj} corrections remembered from a previous call .....
cc        IF((RDIST.GT.0.d0).AND.((IDAT.GT.1).OR.(IDAT.LT.0))) GO TO 50
c*** ! ELSE - make (& write) 'Dattani' Cm{adj} correctons here *********
          CALL quadCORR(NCMM(ISTATE),MCMM,NCMMAX,MMLR1D,
     1                                        DE(ISTATE),Cm1D,CmEFF1D)
c??? try to devise sensible way to skip  quadCORR in late iterations
          IF(MCMM.GT.NCMM(ISTATE)) THEN
              DO m= NCMM(ISTATE),MCMM
                  CmEFF(m,ISTATE)= CmEFF1D(m)
                  ENDDO
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
              CALL AFdiag(RE(ISTATE),NCMM(ISTATE),
     1             NCMMax,MMLR1D,PSEL(ISTATE),DECM(ISTATE),
     2             Cm1D,rhoAB(ISTATE),IVSR(ISTATE),
     3             IDSTT(ISTATE),ULRe,dULRedCm,dLULRedRe,dULRdDe)
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
          DO  I= ISTART,ISTOP
c** Now - generate potential while looping over radial array
              RVAL= RD(I,ISTATE)
              IF(RDIST.GT.0.d0) RVAL= RDIST
              RDp= RVAL**nPB(ISTATE)
              RDp= RVAL**nPB(ISTATE)
              RDq= RVAL**nQB(ISTATE)
              YPE= (RDp-REp)/(RDp+REp)
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
     1                               Scalc(YQ,J,npow,xBETA,rKL,NbetaMX)
                      VAL= VAL+ DBDB(J,I,ISTATE)*BETA(J,ISTATE)
                      ENDDO
                  DBDRe(I,ISTATE)= -DBDB(npow,I,ISTATE)*dLULRedRe
                ENDIF
              BETAFX(I,ISTATE)= VAL
              XTEMP= DEXP(-VAL*YPE)
c** Now begin by generating  uLR(r)
c----- Special Aubert-Frecon cases ------------------------------------
              IF(MMLR(1,ISTATE).LE.0) THEN
c ... generate ULR for Aubert-Frecon type case ...
              CALL AFdiag(RE(ISTATE),NCMM(ISTATE),
     1             NCMMax,MMLR1D,PSEL(ISTATE),DECM(ISTATE),
     2             Cm1D,rhoAB(ISTATE),IVSR(ISTATE),
     3             IDSTT(ISTATE),ULR,dULRdCm,dULRdR,dULRdDe)
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
                  BETADIST= VAL
c... branch to skip derivatives and inclusion of centrifugal & BOB terms
                  IF(IDAT.LE.-1) GOTO 999  
                  ENDIF
              YPP= 2.d0*DE(ISTATE)*(1.0d0-XTEMP)*XTEMP   !! == DER
              IPV= IPVSTART+2
c... derivatives w.r.t long-range parameters ... 
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
     1                                       + VAL*DYPDRE + dLULRedRe)
              IF(APSE(ISTATE).LE.0) THEN
c... derivatives w.r.t. \beta_i for PE-MLR cases...
                  YPP= YPP*YPE*(1 - YP)
                  DO  J= 0, Nbeta(ISTATE)
                      IPV= IPV+1
                      DVtot(IPV,I)= YPP
                      YPP= YPP*YQ
                      ENDDO
c... derivative w.r.t. De  for 'conventional' power-series exponent
                  DVtot(IPVSTART+1,I)= DVDD + YPP*YP*YPE/DE(ISTATE)
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
              ENDDO                !! end of loop over MLR radial array
ccccc Print for testing
          rewind(10)
          write(10,610) (RD(i,ISTATE),vpot(i,istate),BETAFX(i,istate),
     1                              i= 1, NDATPT(ISTATE),OSEL(ISTATE))
  610 FORMAT(/(f10.4,f15.5,f12.6))
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
          yqRe= (REq - AREFQ)/(REq + AREFQ)       !! next - dyq/dr @ r_e
          dyqRedRe = 2.d0*nQB(ISTATE)*REq*AREFQ/(Re(ISTATE)*
     1                                               (REq + AREFQ)**2)
          IF(RREF(ISTATE).LE.0.d0) dyqRedRe= 0.d0
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
              CALL AFdiag(RE(ISTATE),NCMM(ISTATE),
     1             NCMMax,MMLR1D,PSEL(ISTATE),DECM(ISTATE),
     2             Cm1D,rhoAB(ISTATE),IVSR(ISTATE),
     3             IDSTT(ISTATE),ULRe,dULRedCm,dULRe,dULRdDe)
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
              RDQ= RVAL**NQB(ISTATE)
              YQ= (RDQ-AREFQ)/(RDQ+AREFQ)
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
              CALL AFdiag(RVAL,NCMM(ISTATE),
     1             NCMMax,MMLR1D,PSEL(ISTATE),DECM(ISTATE),
     2             Cm1D,rhoAB(ISTATE),IVSR(ISTATE),
     3             IDSTT(ISTATE),ULR,dULRdCm,dULRdR,dULRdDe)
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
              IF(RREF(ISTATE).LE.0.d0) 
     1             DVtot(IPV,I)= DVtot(IPV,I) - Btemp*DbetaRe*dyqRedRe
c... ** when calculating Cm derivatives, dULRe'/dCm has been excluded ** 
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
              IF(RREF(ISTATE).LE.0) DBDRe(I,ISTATE)= 1.d0
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
            RDp= RVAL**nPB(ISTATE)
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
     1                  *(nPB(ISTATE)/RE(ISTATE))/(AGPEF(ISTATE)*RDp +
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
        rewind(10)
        write(10,612) (RD(i,ISTATE),vpot(i,istate),
     1                              i= 1, NDATPT(ISTATE),OSEL(ISTATE))
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
      rewind(10)
      write(10,612) (RD(i,ISTATE),vpot(i,istate),
     1                              i= 1, NDATPT(ISTATE),OSEL(ISTATE))
          FLUSH(10)
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
          IF(ISTOP.GT.ISTART) WRITE(6,602) VMIN,REQQ
  602 FORMAT('    Extended TT potential has   VMIN=',f9.4,'   at  RMIN='
     1                                                          f8.5)
c...... Print for testing            
          rewind(10)
          write(10,612) (RD(i,ISTATE),vpot(i,istate),
     1                              i= 1, NDATPT(ISTATE),OSEL(ISTATE))
          FLUSH(10)
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
          rewind(10)
          write(10,612) (RD(i,ISTATE),vpot(i,istate),
     1                              i= 1, NDATPT(ISTATE),OSEL(ISTATE))
          FLUSH(10)
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
                          IPV= IPV+ 1
                          DVtot(IPV,I)= YPMA*YQP
                          ENDDO
                      ENDIF
                  IPV= IPV + 1
                  DVtot(IPV,I)= YPA
                  IF(LRad(ISTATE).EQ.0) THEN
                      UAR(I,ISTATE)= VAL*YPMA + 
     1                                      YPA*UA(NUA(ISTATE),ISTATE)
                    ELSE                !! Add up the \delta{Cm} terms
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

c===========================================================================
      SUBROUTINE quadCORR(NCMM,MCMM,NCMMAX,MMLR,De,CmVAL,CmEFF)
c===========================================================================
c** subroutine to generate and print MLR CmEFF values incorporating
c  quadratic 'Dattani' corrections to Cm values for both standard 'linear'
c  and A-F diagonalized uLR(r) functions for MLR potentials
c** Return MCMM= NCMM+1  for C9{adj} term for m_1= 3 potentials
c===========================================================================
      INTEGER NCMM,MCMM,NCMMAX,MMLR(NCMMAX)
      REAL*8 De,CmVAL(NCMMAX),CmEFF(NCMMAX)
c----------------------------------------------------------------------
      IF(MMLR(1).GT.0) THEN
c** For 'normal' inverse-power sum MLR case, with or without damping,
c   set up Dattani's 'Quadratic-corrected' effective Cm values 
          IF((MMLR(1).EQ.6).AND.(NCMM.GE.4)) THEN
c... First, consider C6/C12adj(C14adj) for MMLR(m)={6,8,10,(11),12,14} case
              IF(MMLR(4).EQ.12) THEN             ! explicitly MMLR(4)=12
                  CmEFF(4)= CmVAL(4)+ 0.25D0*CmVAL(1)**2/De
                  WRITE(6,710) MMLR(4),MMLR(4),CmEFF(4)
                  ENDIF
              IF(NCMM.GE.5) THEN
                 IF(MMLR(4).EQ.11) THEN         ! implicitly MMLR(5)=12
                     CmEFF(5)= CmVAL(5) + 0.25D0*CmVAL(1)**2/De
                     WRITE(6,710) MMLR(5),MMLR(5),CmEFF(5)
                     IF(NCMM.GE.6) THEN         ! implicitly MMLR(6)=14
                         CmEFF(6)= CmVAL(6)+ 0.5D0*CmVAL(1)*CmVAL(2)/De
                         WRITE(6,710) MMLR(6),MMLR(6),CmEFF(6)
                         ENDIF
                     ENDIF
                 IF(MMLR(4).EQ.12) THEN           ! assuming MMLR(5)=14
                     CmEFF(5)= CmVAL(5) + 0.5D0*CmVAL(1)*CmVAL(2)/De
                     WRITE(6,710) MMLR(5),MMLR(5),CmEFF(5)
                     ENDIF
                 ENDIF
              ENDIF
          IF((MMLR(1).EQ.5).AND.(NCMM.GE.4)) THEN
c... Then, consider C5/C10adj + C12adj for MMLR(m)={5,6,8,10,12,14} cases
              CmEFF(4)= CmVAL(4) + 0.25D0*CmVAL(1)**2/De
              WRITE(6,710) MMLR(4),MMLR(4),CmEFF(4)
              IF(NCMM.GE.5) THEN                 ! introduce C12^{adj}
                  CmEFF(5)= CmVAL(5) + 0.25D0*CmVAL(2)**2/De
                  WRITE(6,710) MMLR(5),MMLR(5),CmEFF(5)
                  IF(NCMM.GE.6) THEN             ! introduce C14^{adj}
                      CmEFF(6)= CmVAL(6) + 0.5D0*CmVAL(2)*CmVAL(3)/De
                      WRITE(6,710) MMLR(6),MMLR(6),CmEFF(6)
                      ENDIF
                  ENDIF
              ENDIF
          IF((MMLR(1).EQ.4).AND.(NCMM.GE.3).and.(MMLR(3).EQ.8)) THEN
c... Then, consider C4/C8adj + C12adj for MMLR(m)={4,6,8,10,12,14} cases
              CmEFF(3)= CmVAL(3) + 0.25D0*CmVAL(1)**2/De
              WRITE(6,712) MMLR(3),MMLR(3),CmEFF(3)
              IF(NCMM.GE.4) THEN                 ! implicitly MMLR(4)=10
                  CmEFF(4)= CmVAL(4) + 0.5D0*CmVAL(1)*CmVAL(2)/De
                  WRITE(6,710) MMLR(4),MMLR(4),CmEFF(4)
                  IF(NCMM.GE.5) THEN             ! implicitly MMLR(5)=12
                      CmEFF(5)= CmVAL(5) + 0.5D0*CmVAL(1)*CmVAL(3)/De
     1                                       + 0.25D0*CmVAL(2)**2/De
                      WRITE(6,710) MMLR(5),MMLR(5),CmEFF(5)
                      IF(NCMM.GE.6) THEN         ! implicitly MMLR(6)=14
                          CmEFF(6)= CmVAL(6)+ 0.5D0*CmVAL(2)*CmVAL(3)/De
     1                                      + 0.5D0*CmVAL(1)*CmVAL(4)/De
                          WRITE(6,710) MMLR(6),MMLR(6),CmEFF(6)
                          ENDIF
                      ENDIF
                  ENDIF
                ENDIF                   !! consider no further adjustment
          IF((MMLR(1).EQ.4).AND.(NCMM.GE.4).and.(MMLR(4).EQ.8)) THEN
c... Then, consider C4/C8adj + C12adj for MMLR(m)={4,6,7,8} cases
                  CmEFF(4)= CmVAL(4) + 0.25D0*CmVAL(1)**2/De
                  WRITE(6,712) MMLR(4),MMLR(4),CmEFF(4)
                  ENDIF
          IF((MMLR(1).EQ.3).AND.(NCMM.GE.2)) THEN
c... Then, consider C3/C6adj + C9adj for MMLR(m)={3,6,8,(9),10,12,14} cases
              CmEFF(2)= CmVAL(2) + 0.25D0*CmVAL(1)**2/De 
              WRITE(6,712) MMLR(2),MMLR(2),CmEFF(2)
              IF(NCMM.GE.3) THEN              ! introduce C9adj & MMLR=9
                  MCMM= NCMM + 1
                  MMLR(MCMM)= 9 
                  CmEFF(MCMM)= 0.5d0*CmVAL(1)*CmEFF(2)/De
                  WRITE(6,714) MMLR(MCMM),CmEFF(MCMM)
                  IF(NCMM.GE.5) THEN             ! implicitly MMLR(5)=12
                      CmEFF(5)= CmVAL(5) + 0.5D0*CmVAL(1)*CmEFF(MCMM)/De
     1                                         + 0.25D0*CmEFF(2)**2/De
                      WRITE(6,710) MMLR(5),MMLR(5),CmEFF(5)
                      IF(NCMM.GE.6) THEN         ! implicitly MMLR(6)=14
                          CmEFF(6)= CmVAL(6)+ 0.5D0*CmEFF(2)*CmVAL(3)/De
                          WRITE(6,710) MMLR(6),MMLR(6),CmEFF(6)
                          ENDIF
                      ENDIF
                  ENDIF
              ENDIF
          ENDIF
c======================================================================= c
c** End of  CmEFF= Cm + CmADJ  setup for non-AF case ===================
  710 Format("  'Quadratic correction' for   C",I2,'(MLR)   yields',
     1  6x,'C',I2,'{adj}=',1PD15.8)
  712 Format("  'Quadratic correction' for   C",I1,'(MLR)    yields',
     1  7x,'C'I1,'{adj}=',1PD15.8)
  714 Format("  'Quadratic corrn' for  MLR(m_1=3)  introduces   C",
     1    I1,'(',A4,',adj) =',1PD15.8)
  716 Format("  'Quadratic correction' for  C",I1,'(Sigma)  yields   C',
     1    I1,'(Sigma,adj)=',1PD15.8)
  718 Format("  'Quadratic correction' for  C",I1,'(^3Pi)   yields   C',
     1    I1,'(^3Pi,adj) =',1PD15.8)
  720 Format("  'Quadratic correction' for  C",I1,'(^1Pi)   yields  C',
     1    I1,'(^1Pi,adj) =',1PD15.8)
c=========================================================================      
      IF(MMLR(1).LE.0) THEN
c** implement Quadratic 'Dattani' MLR corrections for AF cases         
          IF(MMLR(1).GE.-1) THEN         !! first for the 2x2 cases ...
              CmEFF(4)= CmVAL(4) + 0.25*CmVAL(2)**2/De
              CmEFF(5)= CmVAL(5) + 0.25*CmVAL(3)**2/De
              WRITE(6,716) MMLR(4),MMLR(4),CmEFF(4)
              WRITE(6,718) MMLR(5),MMLR(5),CmEFF(5)
c*  prepare C9{adj} coefficients for addition to chosen root
              MMLR(8)= 9               !! These terms added just
              MMLR(9)= 9               !! before exit from  AFdiag
              Cmeff(8)= 0.5*CmVAL(2)*CmEFF(4)/De   
              WRITE(6,714) MMLR(8),'Sigm',CmEFF(8)
              Cmeff(9)= 0.5*CmVAL(3)*CmEFF(5)/De
              WRITE(6,714) MMLR(9),'^3Pi',CmEFF(9)
              ENDIF
          IF(MMLR(1).LE.-2) THEN         !! now for the 3x3 cases ...
              CmEFF(5)= CmVAL(5) + 0.25*CmVAL(2)**2/De
              WRITE(6,716) MMLR(5),MMLR(5),CmEFF(5)
              CmEFF(6)= CmVAL(6) + 0.25*CmVAL(3)**2/De
              WRITE(6,720) MMLR(6),MMLR(6),CmEFF(6)
              CmEFF(7)= CmVAL(7) + 0.25*CmVAL(4)**2/De
              WRITE(6,718) MMLR(7),MMLR(7),CmEFF(7)
c*  prepare C9{adj} coefficients for addition to chosen root
              MMLR(11)= 9               !! These terms added just
              MMLR(12)= 9               !! before exit from  AFdiag
              MMLR(13)= 9
              Cmeff(11)= 0.5*CmVAL(2)*CmEFF(5)/De   
              IF(MMLR(1).EQ.-2) WRITE(6,714) MMLR(11),'Sigm',CmEFF(11)
              Cmeff(12)= 0.5*CmVAL(3)*CmEFF(6)/De
              IF(MMLR(1).EQ.-3) WRITE(6,714) MMLR(12),'^3Pi',CmEFF(12)
              Cmeff(13)= 0.5*CmVAL(4)*CmEFF(7)/De
              IF(MMLR(1).EQ.-4) WRITE(6,714) MMLR(13),'^1Pi',CmEFF(13)
              ENDIF
          ENDIF
      RETURN
      END
c23456789012345678901234567890123456789012345678901234567890123456789012

c***********************************************************************
      SUBROUTINE dampF(r,rhoAB,NCMM,NCMMAX,MMLR,sVRS2,IDSTT,DM,DMP,DMPP)
c** Subroutine to generate values 'Dm' and its first `Dmp' and second
c   'Dmpp' derivatives w.r.t. r of the chosen form of the damping
c    function, for  m= 1 to MMAX.
c---------------------- RJL Version of 21 April 2016 -------------------
c-----------------------------------------------------------------------
c                 Upon Input
c* r - the radial distance in Angsroms (!) 
c* RHOab  'universal' scaling coefficient used for systems other than H_2
c       RHOab= 2*(RHOa*RHOb)/(RHOa+RHOb) where RHOa = (I_p^A/I_p^H)^0.66
c              where I_p^A is the ionization potential of atom A
c              and I_p^H is the ionization potential of atomic hydrogen
c* NCMM  the number of inverse-power terms to be considered
c* MMLR  are the powers of the NCMM inverse-power terms
c* sVRS2 defines damping s.th.  Dm(r)/r^m --> r^{sVRS2/2} as r --> 0
c* IDSTT specifies damping function type:  > 0  use Douketis et al. form 
c                               if  IDSTT .LE. 0  use Tang-Toennies form
c-----------------------------------------------------------------------
c                 Upon Output
c  DM(m) - The value of the damping function for the long range term 
c          C_MMLR(m)/r^MMLR(m)    {m= 1, NCMM}
c  DMP(m): 1'st derivative w.r.t. r of the damping function  DM(m)
c  DMPP(m): 2'nd derivative w.r.t. r of the damping function  DM(m)
c  IF(rhoAB.LE.0.0) return w. DM(m)= 1.0 & DMP(m)=DMPP(m)=0.0 for all m
c-----------------------------------------------------------------------
      INTEGER NCMM,NCMMAX,MMLR(NCMMAX),sVRS2,IDSTT,sVRS2F,FIRST, Lsr,m,
     1  MM,MMAX,MMTEMP
      REAL*8 r,rhoAB,bTT(-2:2),cDS(-4:4),bDS(-4:4),aTT,br,XP,YP,
     1  TK, DM(NCMMAX),DMP(NCMMAX),DMPP(NCMMAX),SM(-3:25),
     2  bpm(20,-4:0), cpm(20,-4:0),ZK
c------------------------------------------------------------------------
c  The following values for the numerical factors used in both TT and DS
c  were  normalized to the Hydrogen data presented
c  by Kreek and Meath in J.Chem.Phys. 50, 2289 (1969).
c  The ratio has been chosen such that  b= FACTOR*(I_p^X / I_p^H)^{2/3}
c  for the homoatomic diatomic species X_2, where I_p^A is the ionization
c------------------------------------------------------------------------
      DATA bTT/2.10d0,2.44d0,2.78d0,3.13d0,3.47d0/
      DATA bDS/2.50d0,2.90d0,3.30d0,3.69d0,3.95d0,0.d0,4.53d0,
     1         0.d0,4.99d0/
      DATA cDS/0.468d0,0.446d0,0.423d0,0.405d0,0.390d0,0.d0,
     1           0.360d0,0.d0,0.340d0/
c...For testing: precise Scolegian values of 'b' and 'c' for s=0 ......
cc    DATA bDS/2.50d0,2.90d0,3.30d0,3.69d0,3.968424883d0,0.d0,4.53d0,
cc    DATA cDS/0.468d0,0.446d0,0.423d0,0.405d0,0.3892460703d0,0.d0,
      DATA FIRST/ 1/
      SAVE FIRST, bpm, cpm
c-----------------------------------------------------------------------
      MMTEMP = MMLR(1)
c      IF(MMLR(1).LE.0) MMLR(1) = 1 !2025: No idea what's going on here.
      IF(RHOab.LE.0) THEN
          DO  m=1,NCMMax
              DM(m)=1.d0
              DMP(m)= 0.d0
              DMPP(m)= 0.d0
              ENDDO
          RETURN
          ENDIF
      IF(IDSTT.LE.0) THEN
c===========================================
c** For Tang-Toennies type damping functions
c===========================================
          Lsr= sVRS2/2
          IF((sVRS2.LT.-4).OR.(sVRS2.GT.4).OR.((2*LSR).NE.sVRS2)) THEN
                WRITE(6,600) 'TT',sVRS2
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
          IF((sVRS2.LT.-4).OR.(sVRS2.GT.4).OR.(sVRS2.EQ.1).OR.
     1                                              (sVRS2.EQ.3)) THEN
              WRITE(6,600) 'DS',sVRS2
              STOP
              ENDIF
          IF(FIRST.EQ.1) THEN
              DO m= 1, 20
                  DO  sVRS2F= -4,0
                      bpm(m,sVRS2F)= bDS(sVRS2F)/DFLOAT(m)
                      cpm(m,sVRS2F)= cDS(sVRS2F)/DSQRT(DFLOAT(m))
                      ENDDO
                  ENDDO
              FIRST= 0 
              ENDIF
          br= rhoAB*r
          DO m= 1, NCMM
              MM= MMLR(m)
              XP= DEXP(-(bpm(MM,sVRS2) + cpm(MM,sVRS2)*br)*br)
              YP= 1.d0 - XP
              ZK= MM + 0.5d0*sVRS2
              DM(m)= YP**ZK
              TK= (bpm(MM,sVRS2) + 2.d0*cpm(MM,sVRS2)*br)*rhoAB
              DMP(m) = ZK*XP*TK*DM(m)/YP
c ... calculate second derivative [for DELR case] {check this!}
              DMPP(m)= (ZK-1.d0)*DMP(m)*(XP*TK)/YP
     1               - DMP(m)*TK + DMP(m)*2.d0*cpm(MM,sVRS2)*rhoAB**2/TK
              ENDDO   
          ENDIF  
      MMLR(1) = MMTEMP
      RETURN
  600 FORMAT(/,' *** ERROR ***  For  ',A2,'-damping functions not yet de
     1fined for   sVRS2=',i3)
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
c***********************************************************************
      SUBROUTINE AFdiag(RDIST,NCMM,NCMMax,MMLR,PSEL,De,Cm,rhoAB,sVSR2,
     1                               IDSTT,ULR,dULRdCm,dULRdR,dULRdDe)
c***********************************************************************
c**   Aubert-Frecon Potential Model for u_{LR}(r)
c***********************************************************************
c** Subroutine to generate, at the onee distance RDIST, an eigenvalue 
c  of the 2x2 or 3x3 long-range interaction matrix described by Eqs.1
c and 10, resp., of J.Mol.Spec.188, 182 (1998) (Aubert-Frecon et al)
c** and its derivatives w.r.t. the C_m long-range parameters.
c***********************************************************************
c==> Input:  r= RDIST, NCMM, m=MMLR & Cm's, rhoAB, sVSR2, IDSTT
c==> Output: ULR, partial derivatives dULRdCm & radial derivative dULRdR
c-----------------------------------------------------------------------
c** Summer 2008 Original Version from Nike Dattani for 3x3 case
c** July 2014 incorporated 2x2 case, removed retardation terms and
c   incorporate damping  ...  by Kai Slaughter
c-----------------------------------------------------------------------
      INTEGER NCMMax
c-----------------------------------------------------------------------
      REAL*8 RDIST,Cm(NCMMax),ULR,dULRdCm(NCMMax),dULRdR,R2,R3,R5,
     1  R6,R8,R9,T1,T0,T2,T0P,T0P23,DDe1,DDe2,DDe3,DELTAE,Modulus,Z,
     2  Dm(NCMMax),Dmp(NCMMax),De,DDe(3,3),Dmpp(NCMMax),rhoAB,A(3,3),
     3  DR(3,3),Q(3,3),DMx(NCMMax,3,3),DMtemp(3,3),DEIGMx(NCMMax,1,1),
     4  DEIGMtemp(1,1),DEIGR(1,1),DEIGDe(1,1),EIGVEC(3,1),RESID(3,1),
     5  W(3),RPOW(NCMMax),dULRdDe
      INTEGER H,I,J,K,L,M,X,NCMM,MMLR(NCMMax),sVSR2,IDSTT,PSEL
c-----------------------------------------------------------------------
      DELTAE=Cm(1)
      R2= 1.d0/RDIST**2
      R3= R2/RDIST
      R5= R2*R3
      R6= R3*R3
      R8= R6*R2
c-----------------------------------------------------------------------
c....... for rhoAB.le.0.0   returns Dm(m)=1 & Dmp(m)=Dmpp(m)=0  
      CALL dampF(RDIST,rhoAB,NCMM,NCMMAX,MMLR,sVSR2,IDSTT,Dm,Dmp,Dmpp)
c-----------------------------------------------------------------------
      IF(MMLR(1).GE.-1) THEN           !!  For the A (0)  or b (-1) state
c***********************************************************************
c************* Aubert Frecon 2x2 case   NCMM= 7  and  ...
c***              Cm(1) = DELTAE
c***              Cm(2) = C3Sig
c***              Cm(3) = C3Pi
c***              Cm(4) = C6Sig
c***              Cm(5) = C6Pi
c***              Cm(6) = C8Sig
c***              Cm(7) = C8Pi
c***********************************************************************
          T1= R3*(Dm(2)*(Cm(2)-Cm(3)) + R3*Dm(4)*(Cm(4)-Cm(5)) + 
     1        R5*Dm(6)*(Cm(6)-Cm(7)))/3.d0
          T0= DSQRT((T1 - Cm(1))**2 + 8.d0*T1**2)
          ULR= 0.5d0*(-Cm(1) + R3*(Dm(2)*(Cm(2)+Cm(3)) + 
     1         R3*Dm(4)*(Cm(4)+Cm(5)) + R5*Dm(6)*(Cm(6)+Cm(7))) + T0)
c-----------------------------------------------------------------------
          IF(MMLR(1).EQ.0) THEN     
              ULR= ULR + Cm(8)*R3*R6          !! add C9{adj correction
              ENDIF
c...  adjustment for the b-state
          IF(MMLR(1).EQ.-1) THEN
              ULR=ULR-T0
              ULR= ULR + Cm(9)*R3*R6          !! add C9{adj correction
              ENDIF
c...  now get derivatives
          T0P= 0.5d0*(9.d0*T1 - Cm(1))/T0
          T0P23= 0.5d0 + T0P/3.d0
c...  another adjustment for the b-state
          IF(MMLR(1).EQ.-1) T0P23=T0P23-2.d0*T0P/3.d0
          dULRdCm(1)= 0.d0
          dULRdCm(2)= R3*(T0P23)
          dULRdCm(3)= R3*(1.d0-T0P23)
          dULRdCm(4)= R6*(T0P23)
          dULRdCm(5)= R6*(1.d0 - T0P23)
          dULRdCm(6)= R8*T0P23
          dULRdCm(7)= R8*(1.d0-T0P23)
          T2        =-T0P*R3*((Dm(2)*(Cm(2)-Cm(3))+R3*(Dm(4)*2.d0*(Cm(4)
     1                -Cm(5))+R2*Dm(6)*8.d0/3.d0*(Cm(6)-Cm(7))))/RDIST
     2                +(Dmp(2)*(Cm(2)-Cm(3))+R3*Dmp(4)*(Cm(4)-Cm(5))+
     3                R2*R3*Dmp(6)*(Cm(6)-Cm(7)))/3.d0)
          dULRdR    = -R3*((1.5d0*Dm(2)*(Cm(2)+Cm(3)) + R3*(Dm(4)*3.d0*
     1                (Cm(4)+Cm(5))+4.d0*Dm(6)*R2*(Cm(6)+Cm(7))))/RDIST
     2                + 0.5d0*(Dmp(2)*(Cm(2)+Cm(3)) + Dmp(4)*R3*(Cm(4)+
     3                Cm(5)) + Dmp(6)*R3*R2*(Cm(6)+Cm(7)))) + T2
c... and a final adjustment for the b-state
          IF(MMLR(1).EQ.-1) dULRdR= dULRdR- 2.d0*T2
c-----------------------------------------------------------------------
      ELSE
c***********************************************************************
c********* Aubert Frecon 3x3 case   NCMM= 10  and ...
c*********        Cm(1) = DELTAE
c*********        Cm(2) = C3Sig
c*********        Cm(3) = C3Pi1
c*********        Cm(4) = C3Pi3
c*********        Cm(5) = C6Sig
c*********        Cm(6) = C6Pi1
c*********        Cm(7) = C6Pi3
c*********        Cm(8) = C8Sig
c*********        Cm(9) = C8Pi1
c*********        Cm(10)= C8Pi3
c***********************************************************************      
c...      Initialize interaction matrix to 0.d0
          DO  I= 1,3
              DO J= 1,3
                  A(I,J)=0.0D0
                  DR(I,J)=0.d0
                  DO  K= 1,NCMMax
                      DMx(K,I,J)=0.d0
                      ENDDO
                  ENDDO
              ENDDO
c...      Prepare interaction matrix  A
          DO  I= 2,NCMM,3
              RPOW(I)= RDIST**MMLR(I)
              A(1,1)=A(1,1)-Dm(I)*(Cm(I)+Cm(I+1)+Cm(I+2))/(3.d0*RPOW(I))
              A(1,2)=A(1,2)-Dm(I)*(Cm(I+2)+Cm(I+1)-2.d0*Cm(I))/(RPOW(I))
              A(1,3)=A(1,3)-Dm(I)*(Cm(I+2)-Cm(I+1))/(RPOW(I))
              A(2,2)= A(2,2)-Dm(I)*(Cm(I+2)+Cm(I+1)+4.d0*Cm(I))
     1                                                  /(6.d0*RPOW(I))
              A(3,3)= A(3,3) - Dm(I)*(Cm(I+2)+Cm(I+1))/(2.d0*RPOW(I))
              ENDDO
          A(1,2) = A(1,2)/(3.d0*DSQRT(2.d0))
          A(2,1) = A(1,2)
          A(2,2) = A(2,2) + DELTAE
          A(2,3) = A(1,3)/(2.d0*DSQRT(3.d0))
          A(1,3) = A(1,3)/(DSQRT(6.d0))
          A(3,1) = A(1,3)
          A(3,2) = A(2,3)
          A(3,3) = A(3,3) + DELTAE
c...      Prepare radial derivative of interaction matrix (? is it needed ?)
          DO  I= 2,NCMM,3
              DR(1,1)= DR(1,1) + Dm(I)*MMLR(I)*(Cm(I)+Cm(I+1)+Cm(I+2))
     1                             /(3.d0*RPOW(I)*RDIST)
     2                    -Dmp(I)*(Cm(I)+Cm(I+1)+Cm(I+2))/(3.d0*RPOW(I))
              DR(1,2)= DR(1,2) + Dm(I)*MMLR(I)*(Cm(I+2)+Cm(I+1)-2.d0*
     1                           Cm(I))/(RPOW(I)*RDIST)
     2                    -Dmp(I)*(Cm(I+2)+Cm(I+1)-2.d0*Cm(I))/(RPOW(I))
              DR(1,3)= DR(1,3) + Dm(I)*MMLR(I)*(Cm(I+2)-Cm(I+1))
     1                            /(RPOW(I)*RDIST)
     2                        -Dmp(I)*(Cm(I+2)-Cm(I+1))/(RPOW(I))
              DR(2,2)= DR(2,2) + Dm(I)*MMLR(I)*(Cm(I+2)+Cm(I+1)+
     1                           4.d0*Cm(I))/(6.d0*RPOW(I)*RDIST)
     2                        -Dmp(I)*(Cm(I+2)+Cm(I+1)+4.d0*Cm(I))
     3                            /(6.d0*RPOW(I))
              DR(3,3)= DR(3,3) + Dm(I)*MMLR(I)*(Cm(I+2)+Cm(I+1))
     1                            /(2.d0*RPOW(I)*RDIST)
     2                        -Dmp(I)*(Cm(I+2)+Cm(I+1))/(2.d0*RPOW(I)) 
              ENDDO
          DR(1,2) = DR(1,2)/(3.d0*DSQRT(2.d0))
          DR(2,1) = DR(1,2)
          DR(2,3) = DR(1,3)/(2.d0*DSQRT(3.d0))
          DR(1,3) = DR(1,3)/(DSQRT(6.d0))
          DR(3,1) = DR(1,3)
          DR(3,2) = DR(2,3)
c...      Partial derivatives of interaction matrix A  w.r.t.  Cm's
          DO  I= 2,NCMM,3 
              DMx(I,1,1)= -Dm(I)/(3.d0*RPOW(I))     !! d{1,1}/dCm{Sig}
              DMx(I+1,1,1)= DMx(I,1,1)              !! d{1,1}/dCm{1Pi}
              DMx(I+2,1,1)= DMx(I,1,1)              !! d{1,1}/dCm{3Pi}
              DMx(I,1,2)= 2.d0*Dm(I)/(3.d0*DSQRT(2.d0)*RPOW(I))
              DMx(I+1,1,2)= -DMx(I,1,2)/2.d0        !! d{1,2}/dCm{1Pi}
              DMx(I+2,1,2)= DMx(I+1,1,2)            !! d{1,2}/dCm{3Pi}
              DMx(I,2,1)= DMx(I,1,2)
              DMx(I+1,2,1)= DMx(I+1,1,2)
              DMx(I+2,2,1)= DMx(I+2,1,2)
              DMx(I,1,3)= 0.d0                    !! no C3{sig} in {1,3}
              DMx(I,3,1)= 0.d0                    !! no C3{sig} in {3,1}
              DMx(I+1,1,3)= Dm(I)/(DSQRT(6.d0)*RPOW(I))
              DMx(I+2,1,3)= -DMx(I+1,1,3)
              DMx(I+1,3,1)= DMx(I+1,1,3)
              DMx(I+2,3,1)= DMx(I+2,1,3)
              DMx(I,2,2)= -2.d0*Dm(I)/(3.d0*RPOW(I))
              DMx(I+1,2,2)=  DMx(I,2,2)/4.d0
              DMx(I+2,2,2)=  DMx(I+1,2,2)
              DMx(I,2,3)= 0.d0
              DMx(I,3,2)= 0.d0
              DMx(I+1,2,3)= Dm(I)/(2.d0*DSQRT(3.d0)*RPOW(I))
              DMx(I+2,2,3)= -DMx(I+1,2,3)
              DMx(I+1,3,2)= DMx(I+1,2,3)                !! by symmetry
              DMx(I+2,3,2)= DMx(I+2,2,3)                !! by symmetry
              DMx(I,3,3)= 0.d0
              DMx(I+1,3,3)= -Dm(I)/(2.d0*RPOW(I))
              DMx(I+2,3,3)= DMx(I+1,3,3)
              IF((RPOW(I).EQ.6).AND.(PSEL.EQ.2)) THEN
c!! For an MLR PEF, adjust derivatives for d/dC3{C6^{adj}} term 
                  DO J= I-3,I-1                     !! for {1,1} terms
                     DMx(J,1,1)= DMx(J,1,1)*(1.d0 + Dm(J)*Cm(J)/
     1                                            (2.d0*De*RPOW(J)))
                     DMx(J,2,2)= DMx(J,2,2)*(1.d0 + Dm(J)*Cm(J)/
     1                                            (2.d0*De*RPOW(J)))
                     DMx(J,1,2)= DMx(J,1,2)*(1.d0 + Dm(J)*Cm(J)/
     1                                            (2.d0*De*RPOW(J)))
                     DMx(J,2,1)= DMx(J,2,1)
                     DMx(J,1,3)= DMx(J,1,3)*(1.d0 + Dm(J)*Cm(J)/
     1                                            (2.d0*De*RPOW(J)))
                     DMx(J,3,1)= DMx(J,3,1)
                     DMx(J,3,3)= DMx(J,3,3)*(1.d0 + Dm(J)*Cm(J)/
     1                                            (2.d0*De*RPOW(J)))
                     ENDDO
c!! and finally ... derivatives w.r.t. De 
                  DDE1= ((Dm(I-3)*Cm(I-3)/2.d0*De)**2/RPOW(I))
                  DDE2= ((Dm(I-2)*Cm(I-2)/2.d0*De)**2/RPOW(I))
                  DDE3= ((Dm(I-1)*Cm(I-1)/2.d0*De)**2/RPOW(I))
                  DDe(1,1)= (DDe1 + DDe2 + DDe3)/3.d0
                  DDe(1,2)= (-2.d0*DDe1 + DDe2 + DDe3)/(3.d0*SQRT(2.d0))
                  DDe(2,1)= DDe(1,2)
                  DDe(1,3)= (-DDe2 + Dde3)/SQRT(6.d0)
                  DDe(3,1)= DDe(1,3)
                  DDe(2,2)= (4.d0*DDe1 + DDe2 + DDe3)/6.d0
                  DDe(2,3)= DDe(1,3)/SQRT(2.d0)
                  DDe(3,2)= DDe(2,3)
                  DDe(3,3)= (DDe2 + DDe3)/2.d0
                  ENDIF   
              ENDDO
c...      Call subroutine to prepare and invert interaction matrix  A
          CALL DSYEVJ3(A,Q,W)
          L=1
c...      Now - identify the lowest eigenvalue of  A  and label it  L
          DO J=2,3
              IF (W(J) .LT. W(L)) THEN
                  L=J
                  ENDIF
              ENDDO
c...      Identifiy the highest eigenvalue of A and label it H
          H=1 
          DO J=2,3
              IF(W(J).GT.W(H)) THEN
                  H=J
                  ENDIF
              ENDDO
c...      Identify the middle eigenvalue of A and label it M
          M=1 
          DO J=2,3
              IF((J.NE.L).AND.(J.NE.H)) M= J
              ENDDO
c...      Select which eigenvalue to use based on user input
          IF(MMLR(1).EQ.-2) THEN 
              X = L
          ELSEIF(MMLR(1).EQ.-3) THEN
              X = M
          ELSE         
              X = H
              ENDIF
c...      determine ULR and eigenvectors
          ULR= -W(X)
          IF(MMLR(1).EQ.-2) ULR= ULR+ Cm(11)*R3*R6        !! C9adj term
          IF((MMLR(1).EQ.-3).OR.(MMLR(1).EQ.-4)) ULR = ULR + DELTAE
          IF(MMLR(1).EQ.-3) ULR= ULR+ Cm(12)*R3*R6        !! C9adj term
          IF(MMLR(1).EQ.-4) ULR= ULR+ Cm(13)*R3*R6        !! C9adj term
!!!!! print for testing   !! print for testing   !! print for testing
cc    WRITE(25,600) RDIST ,ULR, W(1),W(2),W(3)   !! print for testing
cc600 FORMAT(F12.4,1P,D16.7,2x,3D15.7)           !! print for testing
!!!!! print for testing   !! print for testing   !! print for testing
          DO I=1,3      
              EIGVEC(I,1) = Q(I,X)
              ENDDO 
cc  loop over values of m to determine partial derivatives w.r.t. each Cm
          DO I=2,NCMM
             DMtemp(1:3,1:3) = DMx(I,1:3,1:3) 
             DEIGMtemp= -MATMUL(TRANSPOSE(EIGVEC),MATMUL(DMtemp,EIGVEC))
             dULRdCm(I)= DEIGMtemp(1,1)
             ENDDO
          DEIGR = -MATMUL(TRANSPOSE(EIGVEC),MATMUL(DR,EIGVEC))
          dULRdR = DEIGR(1,1)
          DEIGDe = -MATMUL(TRANSPOSE(EIGVEC),MATMUL(DDe,EIGVEC))
          dULRdDe = DEIGDe(1,1)
c------------------------------------------------------------------------
          ENDIF
c------------------------------------------------------------------------
      RETURN
      CONTAINS
c=======================================================================
      SUBROUTINE DSYEVJ3(A, Q, W)
c ----------------------------------------------------------------------
c** Subroutine to setup and diagonalize the matrix  A  and return 
c   eigenvalues W and eigenvector matrix  Q
      INTEGER N, I, X, Y, R
      PARAMETER (N=3)
      REAL*8 A(3,3), Q(3,3), W(3)
      REAL*8 SD, SO, S, C, T, G, H, Z, THETA, THRESH
c     Initialize Q to the identitity matrix
c --- This loop can be omitted if only the eigenvalues are desired ---
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
      WRITE(6,'("DSYEVJ3: No convergence.")')
      END SUBROUTINE DSYEVJ3
      END SUBROUTINE AFdiag
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
c    vibrational level-v of isotopologue-i of state-s [for NDEGB < 0 case]
c** This subroutine reads in band specifications on Channel-5 and writes
c   the transition energy specifications to channel-4
c-----------------------------------------------------------------------
c                         Version of 1 September 2005
c-----------------------------------------------------------------------
cc    INCLUDE 'arrsizes.h'
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c** 'Block' Data Utility routine named: 'arrsizes.h' that governs 
c    array dimensioning in program  dPotFit
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER NISTPMX,NPARMX,NbetaMX,NBOBMX,HPARMX,NDATAMX,
     1  NVIBMX,NBCMX,NSTATEMX,NPNTMX,NROTMX,NCMMAX
c*  NISTPMX  is the maximum number of isotopomers allowed for fit
      PARAMETER (NISTPMX =  12)
c*  NSTATEMX  is maximum no. of electronic states which can be
c             simultaneously fitted to
      PARAMETER (NSTATEMX = 4)
c*  NPARMX  is the largest number of free parameters allowed for fit
c  Since FS origins may be parameters, this is also max. no, data bands
      PARAMETER (NPARMX  = 8000)
c*  NbetaMX  is the largest number of exponent parameters allowed for fit
      PARAMETER (NbetaMX  = 40)
c*  NBOBMX-1  is the highest-order polynomial expansion allowed for the
c      adiabatic or centrifugal Born-Oppenheimer breakdown functions, or 
c      the Lambda-doubling or 2\Sigma splitting radial strength functions
      PARAMETER (NBOBMX  = 15)
c*  HPARMX  is the largest number of Hamiltonian parameters of all types
c    (potential energy, BOB. etc.) for all states.
c           HPARMX >= NSTATEMX*[5 + (NbetaMX+1) + 5*(NBOBMX+1)]
      PARAMETER (HPARMX= NSTATEMX*(5 + (NbetaMX+1) + 5*(NBOBMX+1)))
cc    PARAMETER (HPARMX = 300)
c*  NDATAMX  is largest No. of individual data which may be considered
      PARAMETER (NDATAMX = 35000)
c*  NVIBMX  is the maximum number of vibrational levels of a single
c           state for which data are to be considered
      PARAMETER (NVIBMX    = 200)
** NBCMX  is the maximum number of band constants per vib level to be
c         allowed when doing band constant fits (PSEL= -1)
      PARAMETER (NBCMX = 8) 
c*  NPNTMX  is the largest number of potential data points that can be
c           stored in a single 1D array
      PARAMETER (NPNTMX = 90000)
c*  NROTMX  is the highest order of rotational constants calculated and
c            used for estimating level energies
      PARAMETER (NROTMX = 7)
c*  NCMMAX is the largest number of Cm terms in the MLR or DELR 
c            long-range potential
      PARAMETER (NCMMAX = 12)
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cc    INCLUDE 'BLKISOT.h'
c=======================================================================
c** Isotope/isotopologue numbers, masses & BOB mass scaling factors
c** Array ZK carries about the band constants for all levels of all ISOT
      INTEGER NISTP,NDUNMX,AN(2),MN(2,NISTPMX)
c** NDUNMX is a dummy parameter reqd. for portability of READATA
      PARAMETER (NDUNMX=0)
      REAL*8  ZMASS(3,NISTPMX),RSQMU(NISTPMX),RSQMUP(0:NDUNMX,NISTPMX),
     1 RMUP(0:9,NISTPMX),ZMUA(NISTPMX,NSTATEMX),ZMUB(NISTPMX,NSTATEMX),
     2 ZMTA(NISTPMX,NSTATEMX),ZMTB(NISTPMX,NSTATEMX),
     3  ZK(0:NVIBMX,0:NROTMX,NISTPMX,NSTATEMX)
c
      COMMON /BLKISOT/ZMASS,RSQMU,RSQMUP,RMUP,ZMUA,ZMUB,ZMTA,ZMTB,ZK,
     1  NISTP,AN,MN
c=======================================================================
cc    INCLUDE 'BLKDATA.h'
c=======================================================================
c** Type statements & common block for data
      REAL*8  FREQ(NDATAMX),UFREQ(NDATAMX),DFREQ(NDATAMX),TEMP(NDATAMX),
     1                                               YUNC(NDATAMX),Fqb
      INTEGER  COUNTOT,NFS1,NFSTOT,NBANDTOT,IB(NDATAMX),JP(NDATAMX),
     1 JPP(NDATAMX),VP(NPARMX),VPP(NPARMX),EFP(NDATAMX),EFPP(NDATAMX),
     2 TVUP(NDATAMX),TVLW(NDATAMX),FSBAND(NPARMX),IFXFS(NPARMX),
     3 NFS(NPARMX),IEP(NPARMX),IEPP(NPARMX),ISTP(NPARMX),
     4 IFIRST(NPARMX),ILAST(NPARMX),NTV(NSTATEMX,NISTPMX),FSsame,
     5 NTRANS(NPARMX),IBB(NISTPMX,NSTATEMX,9,NPARMX),JMIN(NPARMX),
     6 JMAX(NPARMX)
      CHARACTER*2 NAME(2) 
      CHARACTER*3 SLABL(-6:NSTATEMX)
      CHARACTER*30 BANDNAME(NPARMX)
      COMMON /DATABLK/Fqb,FREQ,UFREQ,YUNC,DFREQ,TEMP,COUNTOT,NFS1,
     1 NFSTOT,NBANDTOT,IB,JP,JPP,VP,VPP,EFP,EFPP,TVUP,TVLW,FSBAND,IFXFS,
     2 NFS,IEP,IEPP,ISTP,IFIRST,ILAST,NTV,FSsame,
     3 NTRANS,IBB,JMIN,JMAX,NAME,SLABL,BANDNAME
c=======================================================================
cc    INCLUDE 'BLKTYPE.h'
c=======================================================================
c** Type statements & common blocks for characterizing transitions
      REAL*8  AVEUFREQ(NPARMX),MAXUFREQ(NPARMX)
      INTEGER NTRANSFS(NISTPMX,NSTATEMX),
     1  NTRANSVIS(NISTPMX,NSTATEMX,NSTATEMX),
     1  NBANDEL(NISTPMX,NSTATEMX,NSTATEMX),
     2  NTRANSIR(NISTPMX,NSTATEMX),NTRANSMW(NISTPMX,NSTATEMX),
     3  NBANDFS(NISTPMX,NSTATEMX),NBANDVIS(NISTPMX,NSTATEMX),
     4  NBANDIR(NISTPMX,NSTATEMX),NBANDMW(NISTPMX,NSTATEMX),
     5  NVVPP(NISTPMX,NSTATEMX),NWIDTH(NISTPMX,NSTATEMX),
     6  NEBPAS(NISTPMX,NSTATEMX),NVIRIAL(NISTPMX,NSTATEMX),
     7  NAcVIR(NISTPMX,NSTATEMX),NBANDS(NISTPMX)
c
      COMMON /BLKTYPE/AVEUFREQ,MAXUFREQ,NTRANSFS,NTRANSVIS,NTRANSIR,
     1  NTRANSMW,NBANDFS,NBANDEL,NBANDVIS,NBANDIR,NBANDMW,NVVPP,NWIDTH,
     2  NEBPAS,NVIRIAL,NAcVIR,NBANDS
c=======================================================================
c-----------------------------------------------------------------------
c
      CHARACTER*3 LABLP,LABLPP
      INTEGER I,J,J2,JD,J2DL,J2DU,J2DD,JMAXX,PP,PPP,NTRANST,COUNT,
     1  IBAND,JMAXP(NPARMX),JMINP(NPARMX),
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
c** MN1 & MN2 identify the isotopologue
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
c** Determine the correct isotopologue-number for this band.
      DO  I= 1,NISTP
          IF((MN1.EQ.MN(1,I)).AND.(MN2.EQ.MN(2,I))) ISOT= I
          ENDDO
      ISTP(IBAND)= ISOT
      MAXUFREQ(IBAND)= 0
      JMAXP(IBAND)= JMAXX
      JMINP(IBAND)= 0
      NTRANST= 0
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
  400 FORMAT(2I4,"   '",A3,"'  '",A3,"'   ",2I4,"     'predictions' ")
  402 FORMAT(I4,I3,I5,I3,'    0.d0     1.0d-3')
  404 FORMAT('  -1 -1   -1 -1    -1.d0    -1.d-3'/)
      VMX(ESP)= MAX(VMX(ESP),VP(IBAND))
      VMX(ESPP)= MAX(VMX(ESPP),VPP(IBAND))
      ILAST(IBAND)= COUNT
      NTRANST= ILAST(IBAND)-IFIRST(IBAND)+1
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
cc    INCLUDE 'arrsizes.h'
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c** 'Block' Data Utility routine named: 'arrsizes.h' that governs 
c    array dimensioning in program  dPotFit
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER NISTPMX,NPARMX,NbetaMX,NBOBMX,HPARMX,NDATAMX,
     1  NVIBMX,NBCMX,NSTATEMX,NPNTMX,NROTMX,NCMMAX
c*  NISTPMX  is the maximum number of isotopomers allowed for fit
      PARAMETER (NISTPMX =  12)
c*  NSTATEMX  is maximum no. of electronic states which can be
c             simultaneously fitted to
      PARAMETER (NSTATEMX = 4)
c*  NPARMX  is the largest number of free parameters allowed for fit
c  Since FS origins may be parameters, this is also max. no, data bands
      PARAMETER (NPARMX  = 8000)
c*  NbetaMX  is the largest number of exponent parameters allowed for fit
      PARAMETER (NbetaMX  = 40)
c*  NBOBMX-1  is the highest-order polynomial expansion allowed for the
c      adiabatic or centrifugal Born-Oppenheimer breakdown functions, or 
c      the Lambda-doubling or 2\Sigma splitting radial strength functions
      PARAMETER (NBOBMX  = 15)
c*  HPARMX  is the largest number of Hamiltonian parameters of all types
c    (potential energy, BOB. etc.) for all states.
c           HPARMX >= NSTATEMX*[5 + (NbetaMX+1) + 5*(NBOBMX+1)]
      PARAMETER (HPARMX= NSTATEMX*(5 + (NbetaMX+1) + 5*(NBOBMX+1)))
cc    PARAMETER (HPARMX = 300)
c*  NDATAMX  is largest No. of individual data which may be considered
      PARAMETER (NDATAMX = 35000)
c*  NVIBMX  is the maximum number of vibrational levels of a single
c           state for which data are to be considered
      PARAMETER (NVIBMX    = 200)
** NBCMX  is the maximum number of band constants per vib level to be
c         allowed when doing band constant fits (PSEL= -1)
      PARAMETER (NBCMX = 8) 
c*  NPNTMX  is the largest number of potential data points that can be
c           stored in a single 1D array
      PARAMETER (NPNTMX = 90000)
c*  NROTMX  is the highest order of rotational constants calculated and
c            used for estimating level energies
      PARAMETER (NROTMX = 7)
c*  NCMMAX is the largest number of Cm terms in the MLR or DELR 
c            long-range potential
      PARAMETER (NCMMAX = 12)
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cc    INCLUDE 'BLKISOT.h'
c=======================================================================
c** Isotope/isotopologue numbers, masses & BOB mass scaling factors
c** Array ZK carries about the band constants for all levels of all ISOT
      INTEGER NISTP,NDUNMX,AN(2),MN(2,NISTPMX)
c** NDUNMX is a dummy parameter reqd. for portability of READATA
      PARAMETER (NDUNMX=0)
      REAL*8  ZMASS(3,NISTPMX),RSQMU(NISTPMX),RSQMUP(0:NDUNMX,NISTPMX),
     1 RMUP(0:9,NISTPMX),ZMUA(NISTPMX,NSTATEMX),ZMUB(NISTPMX,NSTATEMX),
     2 ZMTA(NISTPMX,NSTATEMX),ZMTB(NISTPMX,NSTATEMX),
     3  ZK(0:NVIBMX,0:NROTMX,NISTPMX,NSTATEMX)
c
      COMMON /BLKISOT/ZMASS,RSQMU,RSQMUP,RMUP,ZMUA,ZMUB,ZMTA,ZMTB,ZK,
     1  NISTP,AN,MN
c=======================================================================
cc    INCLUDE 'BLKDATA.h'
c=======================================================================
c** Type statements & common block for data
      REAL*8  FREQ(NDATAMX),UFREQ(NDATAMX),DFREQ(NDATAMX),TEMP(NDATAMX),
     1                                               YUNC(NDATAMX),Fqb
      INTEGER  COUNTOT,NFS1,NFSTOT,NBANDTOT,IB(NDATAMX),JP(NDATAMX),
     1 JPP(NDATAMX),VP(NPARMX),VPP(NPARMX),EFP(NDATAMX),EFPP(NDATAMX),
     2 TVUP(NDATAMX),TVLW(NDATAMX),FSBAND(NPARMX),IFXFS(NPARMX),
     3 NFS(NPARMX),IEP(NPARMX),IEPP(NPARMX),ISTP(NPARMX),
     4 IFIRST(NPARMX),ILAST(NPARMX),NTV(NSTATEMX,NISTPMX),FSsame,
     5 NTRANS(NPARMX),IBB(NISTPMX,NSTATEMX,9,NPARMX),JMIN(NPARMX),
     6 JMAX(NPARMX)
      CHARACTER*2 NAME(2) 
      CHARACTER*3 SLABL(-6:NSTATEMX)
      CHARACTER*30 BANDNAME(NPARMX)
      COMMON /DATABLK/Fqb,FREQ,UFREQ,YUNC,DFREQ,TEMP,COUNTOT,NFS1,
     1 NFSTOT,NBANDTOT,IB,JP,JPP,VP,VPP,EFP,EFPP,TVUP,TVLW,FSBAND,IFXFS,
     2 NFS,IEP,IEPP,ISTP,IFIRST,ILAST,NTV,FSsame,
     3 NTRANS,IBB,JMIN,JMAX,NAME,SLABL,BANDNAME
c=======================================================================
cc    INCLUDE 'BLKPOT.h'
c=======================================================================
c** Effective adiabatic radial potential variables.
      INTEGER BOBCN(NSTATEMX),PSEL(NSTATEMX),MAXMIN(NSTATEMX),
     1 IOMEG(NSTATEMX),Nbeta(NSTATEMX),APSE(NSTATEMX),IFXDE(NSTATEMX),
     2 IFXRE(NSTATEMX),IFXCm(NCMMax,NSTATEMX),
     3 IFXBETA(0:NbetaMX,NSTATEMX),NDATPT(NSTATEMX),NCMM(NSTATEMX),
     4 MMLR(NCMMax,NSTATEMX),nPB(NSTATEMX),nQB(NSTATEMX),pAD(NSTATEMX),
     5 qAD(NSTATEMX),LRad(NSTATEMX),pNA(NSTATEMX),qNA(NSTATEMX),
     6 Pqw(NSTATEMX),IVSR(NSTATEMX),IDSTT(NSTATEMX)
c
      REAL*8 DE(NSTATEMX),RE(NSTATEMX),BETA(0:NbetaMX,NSTATEMX),
     1 yqBETA(NbetaMX,NSTATEMX),BETAFX(NPNTMX,NSTATEMX),RH(NSTATEMX),
     2 RMIN(NSTATEMX),RMAX(NSTATEMX),VLIM(NSTATEMX),EPS(NSTATEMX),
     3 betaINF(NSTATEMX),AGPEF(NSTATEMX),BGPEF(NSTATEMX),
     4 CmVAL(NCMMax,NSTATEMX),CmEFF(NCMMax,NSTATEMX),rhoAB(NSTATEMX),
     5 AA(NSTATEMX),BB(NSTATEMX),RREF(NSTATEMX),ASO(NSTATEMX),
     6 R01(NSTATEMX),Q12(NSTATEMX),RD(NPNTMX,NSTATEMX),
     7 VPOT(NPNTMX,NSTATEMX),dCmA(NCMMax,NSTATEMX),dCmB(NCMMax,NSTATEMX)
c
      COMMON /BLKPOT/DE,RE,BETA,yqBETA,BETAFX,RH,RMIN,RMAX,VLIM,EPS,
     1 betaINF,AGPEF,BGPEF,CmVAL,CmEFF,rhoAB,AA,BB,RREF,ASO,R01,Q12,RD,
     2 VPOT,dCmA,dCmB, BOBCN,PSEL,MAXMIN,IOMEG,Nbeta,APSE,IFXDE,IFXRE,
     3 IFXCm,IFXBETA,NDATPT,NCMM,MMLR,nPB,nQB,pAD,qAD,LRad,pNA,qNA,Pqw,
     4 IVSR,IDSTT
c=======================================================================
cc    INCLUDE 'BLKPARAM.h'
c=======================================================================
c** Parameters and count-labels for band constant (PSEL=-1) or term
c   value (PSEL=-2) fits
      REAL*8 TVALUE(NPARMX),ZBC(0:NVIBMX,0:NROTMX,NISTPMX,NSTATEMX),
     1 ZQC(0:NVIBMX,0:NROTMX,NISTPMX,NSTATEMX)
c
      INTEGER NSTATES,NTVALL(0:NSTATEMX),NTVI(NSTATEMX),NTVF(NSTATEMX),
     1 VMIN(NSTATEMX,NISTPMX),VMAX(NSTATEMX,NISTPMX),JTRUNC(NSTATEMX),
     2 EFSEL(NSTATEMX),NBC(0:NVIBMX,NISTPMX,NSTATEMX),
     3 NQC(0:NVIBMX,NISTPMX,NSTATEMX),
     4 BCPARI(0:NVIBMX,NISTPMX,NSTATEMX),
     5 BCPARF(0:NVIBMX,NISTPMX,NSTATEMX),
     6 QCPARI(0:NVIBMX,NISTPMX,NSTATEMX),
     7 QCPARF(0:NVIBMX,NISTPMX,NSTATEMX)
      COMMON /BLKPARAM/TVALUE,ZBC,ZQC,NSTATES,NTVALL,NTVI,NTVF,VMIN,
     1      VMAX,JTRUNC,EFSEL,NBC,NQC,BCPARI,BCPARF,QCPARI,QCPARF
c=======================================================================
cc    INCLUDE 'BLKBOB.h'
c=======================================================================
c** Born-Oppenheimer Breakdown & doubling function parameters.
c**                       March 16 2012
c=======================================================================
      INTEGER NUA(NSTATEMX),NUB(NSTATEMX),NTA(NSTATEMX),NTB(NSTATEMX),
     1  IFXUA(0:NBOBMX,NSTATEMX),IFXUB(0:NBOBMX,NSTATEMX),
     2  IFXTA(0:NBOBMX,NSTATEMX),IFXTB(0:NBOBMX,NSTATEMX),
     3  NwCFT(NSTATEMX),IFXwCFT(0:NBOBMX,NSTATEMX),efREF(NSTATEMX)
c
      REAL*8 UA(0:NBOBMX,NSTATEMX),UB(0:NBOBMX,NSTATEMX),
     1  TA(0:NBOBMX,NSTATEMX),TB(0:NBOBMX,NSTATEMX),
     2   wCFT(0:NBOBMX,NSTATEMX)
c
      COMMON /BLKBOB/UA,UB,TA,TB,wCFT,NUA,NUB,NTA,NTB,NwCFT,
     1  IFXUA,IFXUB,IFXTA,IFXTB,IFXwCFT,efREF
c=======================================================================
cc    INCLUDE 'BLKCOUNT.h'
c=======================================================================
c    Block data file  BLKCOUNT.h
c=======================================================================
c** Counters for numbers of potential parameters of different types for 
c   each state
      INTEGER  TOTPOTPAR,POTPARI(NSTATEMX),POTPARF(NSTATEMX),
     1  UAPARI(NSTATEMX),UAPARF(NSTATEMX),UBPARI(NSTATEMX),
     2  UBPARF(NSTATEMX),TAPARI(NSTATEMX),TAPARF(NSTATEMX),
     3  TBPARI(NSTATEMX),TBPARF(NSTATEMX),LDPARI(NSTATEMX),
     4  LDPARF(NSTATEMX),HPARF(NSTATEMX),OSEL(NSTATEMX)
c
      COMMON /BLKCOUNT/TOTPOTPAR,POTPARI,POTPARF,UAPARI,UAPARF,UBPARI,
     1  UBPARF,TAPARI,TAPARF,TBPARI,TBPARF,LDPARI,LDPARF,HPARF,OSEL
c=======================================================================
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
          CALL DEDP(IDAT,IEPP(NBAND),IISTP,ZMASS(3,IISTP),
     1           JP(IDAT),JPP(IDAT),EFPP(IDAT),ELW,VMAXX(IEPP(NBAND)),
     2                                             width,LOWER,fcount)
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
          CALL DEDP(IDAT,IEPP(NBAND),IISTP,ZMASS(3,IISTP),
     1          JP(IDAT),JPP(IDAT),EFPP(IDAT),ELW,VMAXX(IEPP(NBAND)),
     2                                             width,LOWER,fcount)
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
          CALL DEDP(IDAT,IEPP(NBAND),IISTP,ZMASS(3,IISTP),
     1           JP(IDAT),JPP(IDAT),EFPP(IDAT),EO,VMAXX(IEPP(NBAND)),
     2                                             width,DEDPK,fcount)
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
              CALL DEDP(IDAT,IEP(NBAND),IISTP,ZMASS(3,IISTP),
     1             VP(NBAND),JP(IDAT),EFP(IDAT),EUP,VMAXX(IEP(NBAND)),
     2                                             width,UPPER,fcount)
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
cc    INCLUDE 'arrsizes.h'
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c** 'Block' Data Utility routine named: 'arrsizes.h' that governs 
c    array dimensioning in program  dPotFit
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER NISTPMX,NPARMX,NbetaMX,NBOBMX,HPARMX,NDATAMX,
     1  NVIBMX,NBCMX,NSTATEMX,NPNTMX,NROTMX,NCMMAX
c*  NISTPMX  is the maximum number of isotopomers allowed for fit
      PARAMETER (NISTPMX =  12)
c*  NSTATEMX  is maximum no. of electronic states which can be
c             simultaneously fitted to
      PARAMETER (NSTATEMX = 4)
c*  NPARMX  is the largest number of free parameters allowed for fit
c  Since FS origins may be parameters, this is also max. no, data bands
      PARAMETER (NPARMX  = 8000)
c*  NbetaMX  is the largest number of exponent parameters allowed for fit
      PARAMETER (NbetaMX  = 40)
c*  NBOBMX-1  is the highest-order polynomial expansion allowed for the
c      adiabatic or centrifugal Born-Oppenheimer breakdown functions, or 
c      the Lambda-doubling or 2\Sigma splitting radial strength functions
      PARAMETER (NBOBMX  = 15)
c*  HPARMX  is the largest number of Hamiltonian parameters of all types
c    (potential energy, BOB. etc.) for all states.
c           HPARMX >= NSTATEMX*[5 + (NbetaMX+1) + 5*(NBOBMX+1)]
      PARAMETER (HPARMX= NSTATEMX*(5 + (NbetaMX+1) + 5*(NBOBMX+1)))
cc    PARAMETER (HPARMX = 300)
c*  NDATAMX  is largest No. of individual data which may be considered
      PARAMETER (NDATAMX = 35000)
c*  NVIBMX  is the maximum number of vibrational levels of a single
c           state for which data are to be considered
      PARAMETER (NVIBMX    = 200)
** NBCMX  is the maximum number of band constants per vib level to be
c         allowed when doing band constant fits (PSEL= -1)
      PARAMETER (NBCMX = 8) 
c*  NPNTMX  is the largest number of potential data points that can be
c           stored in a single 1D array
      PARAMETER (NPNTMX = 90000)
c*  NROTMX  is the highest order of rotational constants calculated and
c            used for estimating level energies
      PARAMETER (NROTMX = 7)
c*  NCMMAX is the largest number of Cm terms in the MLR or DELR 
c            long-range potential
      PARAMETER (NCMMAX = 12)
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cc    INCLUDE 'BLKBOB.h'
c=======================================================================
c** Born-Oppenheimer Breakdown & doubling function parameters.
c**                       March 16 2012
c=======================================================================
      INTEGER NUA(NSTATEMX),NUB(NSTATEMX),NTA(NSTATEMX),NTB(NSTATEMX),
     1  IFXUA(0:NBOBMX,NSTATEMX),IFXUB(0:NBOBMX,NSTATEMX),
     2  IFXTA(0:NBOBMX,NSTATEMX),IFXTB(0:NBOBMX,NSTATEMX),
     3  NwCFT(NSTATEMX),IFXwCFT(0:NBOBMX,NSTATEMX),efREF(NSTATEMX)
c
      REAL*8 UA(0:NBOBMX,NSTATEMX),UB(0:NBOBMX,NSTATEMX),
     1  TA(0:NBOBMX,NSTATEMX),TB(0:NBOBMX,NSTATEMX),
     2   wCFT(0:NBOBMX,NSTATEMX)
c
      COMMON /BLKBOB/UA,UB,TA,TB,wCFT,NUA,NUB,NTA,NTB,NwCFT,
     1  IFXUA,IFXUB,IFXTA,IFXTB,IFXwCFT,efREF
c=======================================================================
cc    INCLUDE 'BLKBOBRF.h'
c=======================================================================
c** Born-Oppenheimer breakdown radial functions 
      REAL*8 UAR(NPNTMX,NSTATEMX),UBR(NPNTMX,NSTATEMX),
     1 TAR(NPNTMX,NSTATEMX),TBR(NPNTMX,NSTATEMX),wRAD(NPNTMX,NSTATEMX)
c
      COMMON /BLKBOBRF/UAR,UBR,TAR,TBR,wRAD
c=======================================================================
cc    INCLUDE 'BLKDATA.h'
c=======================================================================
c** Type statements & common block for data
      REAL*8  FREQ(NDATAMX),UFREQ(NDATAMX),DFREQ(NDATAMX),TEMP(NDATAMX),
     1                                               YUNC(NDATAMX),Fqb
      INTEGER  COUNTOT,NFS1,NFSTOT,NBANDTOT,IB(NDATAMX),JP(NDATAMX),
     1 JPP(NDATAMX),VP(NPARMX),VPP(NPARMX),EFP(NDATAMX),EFPP(NDATAMX),
     2 TVUP(NDATAMX),TVLW(NDATAMX),FSBAND(NPARMX),IFXFS(NPARMX),
     3 NFS(NPARMX),IEP(NPARMX),IEPP(NPARMX),ISTP(NPARMX),
     4 IFIRST(NPARMX),ILAST(NPARMX),NTV(NSTATEMX,NISTPMX),FSsame,
     5 NTRANS(NPARMX),IBB(NISTPMX,NSTATEMX,9,NPARMX),JMIN(NPARMX),
     6 JMAX(NPARMX)
      CHARACTER*2 NAME(2) 
      CHARACTER*3 SLABL(-6:NSTATEMX)
      CHARACTER*30 BANDNAME(NPARMX)
      COMMON /DATABLK/Fqb,FREQ,UFREQ,YUNC,DFREQ,TEMP,COUNTOT,NFS1,
     1 NFSTOT,NBANDTOT,IB,JP,JPP,VP,VPP,EFP,EFPP,TVUP,TVLW,FSBAND,IFXFS,
     2 NFS,IEP,IEPP,ISTP,IFIRST,ILAST,NTV,FSsame,
     3 NTRANS,IBB,JMIN,JMAX,NAME,SLABL,BANDNAME
c=======================================================================
cc    INCLUDE 'BLKISOT.h'
c=======================================================================
c** Isotope/isotopologue numbers, masses & BOB mass scaling factors
c** Array ZK carries about the band constants for all levels of all ISOT
      INTEGER NISTP,NDUNMX,AN(2),MN(2,NISTPMX)
c** NDUNMX is a dummy parameter reqd. for portability of READATA
      PARAMETER (NDUNMX=0)
      REAL*8  ZMASS(3,NISTPMX),RSQMU(NISTPMX),RSQMUP(0:NDUNMX,NISTPMX),
     1 RMUP(0:9,NISTPMX),ZMUA(NISTPMX,NSTATEMX),ZMUB(NISTPMX,NSTATEMX),
     2 ZMTA(NISTPMX,NSTATEMX),ZMTB(NISTPMX,NSTATEMX),
     3  ZK(0:NVIBMX,0:NROTMX,NISTPMX,NSTATEMX)
c
      COMMON /BLKISOT/ZMASS,RSQMU,RSQMUP,RMUP,ZMUA,ZMUB,ZMTA,ZMTB,ZK,
     1  NISTP,AN,MN
c=======================================================================
cc    INCLUDE 'BLKPOT.h'
c=======================================================================
c** Effective adiabatic radial potential variables.
      INTEGER BOBCN(NSTATEMX),PSEL(NSTATEMX),MAXMIN(NSTATEMX),
     1 IOMEG(NSTATEMX),Nbeta(NSTATEMX),APSE(NSTATEMX),IFXDE(NSTATEMX),
     2 IFXRE(NSTATEMX),IFXCm(NCMMax,NSTATEMX),
     3 IFXBETA(0:NbetaMX,NSTATEMX),NDATPT(NSTATEMX),NCMM(NSTATEMX),
     4 MMLR(NCMMax,NSTATEMX),nPB(NSTATEMX),nQB(NSTATEMX),pAD(NSTATEMX),
     5 qAD(NSTATEMX),LRad(NSTATEMX),pNA(NSTATEMX),qNA(NSTATEMX),
     6 Pqw(NSTATEMX),IVSR(NSTATEMX),IDSTT(NSTATEMX)
c
      REAL*8 DE(NSTATEMX),RE(NSTATEMX),BETA(0:NbetaMX,NSTATEMX),
     1 yqBETA(NbetaMX,NSTATEMX),BETAFX(NPNTMX,NSTATEMX),RH(NSTATEMX),
     2 RMIN(NSTATEMX),RMAX(NSTATEMX),VLIM(NSTATEMX),EPS(NSTATEMX),
     3 betaINF(NSTATEMX),AGPEF(NSTATEMX),BGPEF(NSTATEMX),
     4 CmVAL(NCMMax,NSTATEMX),CmEFF(NCMMax,NSTATEMX),rhoAB(NSTATEMX),
     5 AA(NSTATEMX),BB(NSTATEMX),RREF(NSTATEMX),ASO(NSTATEMX),
     6 R01(NSTATEMX),Q12(NSTATEMX),RD(NPNTMX,NSTATEMX),
     7 VPOT(NPNTMX,NSTATEMX),dCmA(NCMMax,NSTATEMX),dCmB(NCMMax,NSTATEMX)
c
      COMMON /BLKPOT/DE,RE,BETA,yqBETA,BETAFX,RH,RMIN,RMAX,VLIM,EPS,
     1 betaINF,AGPEF,BGPEF,CmVAL,CmEFF,rhoAB,AA,BB,RREF,ASO,R01,Q12,RD,
     2 VPOT,dCmA,dCmB, BOBCN,PSEL,MAXMIN,IOMEG,Nbeta,APSE,IFXDE,IFXRE,
     3 IFXCm,IFXBETA,NDATPT,NCMM,MMLR,nPB,nQB,pAD,qAD,LRad,pNA,qNA,Pqw,
     4 IVSR,IDSTT
c=======================================================================
c-----------------------------------------------------------------------
c
      INTEGER ISTATE,IISTP,VIMX,AFLAG,VIN,NBEG,NEND,WARN,IWR,LPRWF,I,J,
     1   INNODE,KV,NCN
c
      REAL*8 GV(0:NVIBMX),RCNST(NROTMX),RR(NPNTMX),RM2(NPNTMX),
     1  V(NPNTMX),SWF(NPNTMX),FWHM,PMAX,BFCT,BvWN,RHSQ,C3gu,T3,VMAXX
c
      INTEGER INNR(0:NVIBMX)
      REAL*8 SWF2,qDBL,qFCT,Cm1
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      VIN= VIMX
      AFLAG= 0
      WARN= 0
      IWR= -0           !!! should in general be -1
cc    IWR=  5           !!! use this when trouble shooting level search
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
      IF((NCMM(ISTATE).GE.3).AND.(MMLR(2,ISTATE).LE.0).AND.(AN(1).EQ.3)
     1        .AND.(AN(2).EQ.3).AND.(MN(1,IISTP).NE.MN(2,IISTP))) THEN
c** Add g/u symmetry breakdown correction for special case {6,7}Li2(A) !!
          C3gu= (2.d0/3.d0)*CmVAL(1,ISTATE)
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
      NCN= 999
      Cm1 = CmVAL(1,ISTATE)
      IF(MMLR(1,ISTATE).LE.0) Cm1 = CmVAL(2,ISTATE) 
      IF((PSEL(ISTATE).GE.2).AND.(Cm1.GT.0.d0)) 
     1    NCN= MMLR(1,ISTATE)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Call subroutine ALF that will locate the needed vibrational levels 
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      CALL ALF(NDATPT(ISTATE),RH(ISTATE),NCN,RR,V,SWF,VLIM(ISTATE),
     1  MAXMIN(ISTATE),VIN,NVIBMX,VMAXX,AFLAG,ZMASS(3,IISTP),
     2  EPS(ISTATE),GV,INNODE,INNR,IWR)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** If a serious error occured during within ALF, then print out a 
c   warning, record the constants that we have and hope the program
c   doesn't call on the constants for which ALF could not calculate.
c
      IF (AFLAG.LT.0) THEN
          WRITE(6,600) ISTATE,IISTP
          IF(AFLAG.EQ.-1) THEN
              WRITE(6,601) VIN,VIMX, AFLAG
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
cc    INCLUDE 'arrsizes.h'
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c** 'Block' Data Utility routine named: 'arrsizes.h' that governs 
c    array dimensioning in program  dPotFit
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER NISTPMX,NPARMX,NbetaMX,NBOBMX,HPARMX,NDATAMX,
     1  NVIBMX,NBCMX,NSTATEMX,NPNTMX,NROTMX,NCMMAX
c*  NISTPMX  is the maximum number of isotopomers allowed for fit
      PARAMETER (NISTPMX =  12)
c*  NSTATEMX  is maximum no. of electronic states which can be
c             simultaneously fitted to
      PARAMETER (NSTATEMX = 4)
c*  NPARMX  is the largest number of free parameters allowed for fit
c  Since FS origins may be parameters, this is also max. no, data bands
      PARAMETER (NPARMX  = 8000)
c*  NbetaMX  is the largest number of exponent parameters allowed for fit
      PARAMETER (NbetaMX  = 40)
c*  NBOBMX-1  is the highest-order polynomial expansion allowed for the
c      adiabatic or centrifugal Born-Oppenheimer breakdown functions, or 
c      the Lambda-doubling or 2\Sigma splitting radial strength functions
      PARAMETER (NBOBMX  = 15)
c*  HPARMX  is the largest number of Hamiltonian parameters of all types
c    (potential energy, BOB. etc.) for all states.
c           HPARMX >= NSTATEMX*[5 + (NbetaMX+1) + 5*(NBOBMX+1)]
      PARAMETER (HPARMX= NSTATEMX*(5 + (NbetaMX+1) + 5*(NBOBMX+1)))
cc    PARAMETER (HPARMX = 300)
c*  NDATAMX  is largest No. of individual data which may be considered
      PARAMETER (NDATAMX = 35000)
c*  NVIBMX  is the maximum number of vibrational levels of a single
c           state for which data are to be considered
      PARAMETER (NVIBMX    = 200)
** NBCMX  is the maximum number of band constants per vib level to be
c         allowed when doing band constant fits (PSEL= -1)
      PARAMETER (NBCMX = 8) 
c*  NPNTMX  is the largest number of potential data points that can be
c           stored in a single 1D array
      PARAMETER (NPNTMX = 90000)
c*  NROTMX  is the highest order of rotational constants calculated and
c            used for estimating level energies
      PARAMETER (NROTMX = 7)
c*  NCMMAX is the largest number of Cm terms in the MLR or DELR 
c            long-range potential
      PARAMETER (NCMMAX = 12)
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cc    INCLUDE 'BLKISOT.h'
c=======================================================================
c** Isotope/isotopologue numbers, masses & BOB mass scaling factors
c** Array ZK carries about the band constants for all levels of all ISOT
      INTEGER NISTP,NDUNMX,AN(2),MN(2,NISTPMX)
c** NDUNMX is a dummy parameter reqd. for portability of READATA
      PARAMETER (NDUNMX=0)
      REAL*8  ZMASS(3,NISTPMX),RSQMU(NISTPMX),RSQMUP(0:NDUNMX,NISTPMX),
     1 RMUP(0:9,NISTPMX),ZMUA(NISTPMX,NSTATEMX),ZMUB(NISTPMX,NSTATEMX),
     2 ZMTA(NISTPMX,NSTATEMX),ZMTB(NISTPMX,NSTATEMX),
     3  ZK(0:NVIBMX,0:NROTMX,NISTPMX,NSTATEMX)
c
      COMMON /BLKISOT/ZMASS,RSQMU,RSQMUP,RMUP,ZMUA,ZMUB,ZMTA,ZMTB,ZK,
     1  NISTP,AN,MN
c=======================================================================
cc    INCLUDE 'BLKDATA.h'
c=======================================================================
c** Type statements & common block for data
      REAL*8  FREQ(NDATAMX),UFREQ(NDATAMX),DFREQ(NDATAMX),TEMP(NDATAMX),
     1                                               YUNC(NDATAMX),Fqb
      INTEGER  COUNTOT,NFS1,NFSTOT,NBANDTOT,IB(NDATAMX),JP(NDATAMX),
     1 JPP(NDATAMX),VP(NPARMX),VPP(NPARMX),EFP(NDATAMX),EFPP(NDATAMX),
     2 TVUP(NDATAMX),TVLW(NDATAMX),FSBAND(NPARMX),IFXFS(NPARMX),
     3 NFS(NPARMX),IEP(NPARMX),IEPP(NPARMX),ISTP(NPARMX),
     4 IFIRST(NPARMX),ILAST(NPARMX),NTV(NSTATEMX,NISTPMX),FSsame,
     5 NTRANS(NPARMX),IBB(NISTPMX,NSTATEMX,9,NPARMX),JMIN(NPARMX),
     6 JMAX(NPARMX)
      CHARACTER*2 NAME(2) 
      CHARACTER*3 SLABL(-6:NSTATEMX)
      CHARACTER*30 BANDNAME(NPARMX)
      COMMON /DATABLK/Fqb,FREQ,UFREQ,YUNC,DFREQ,TEMP,COUNTOT,NFS1,
     1 NFSTOT,NBANDTOT,IB,JP,JPP,VP,VPP,EFP,EFPP,TVUP,TVLW,FSBAND,IFXFS,
     2 NFS,IEP,IEPP,ISTP,IFIRST,ILAST,NTV,FSsame,
     3 NTRANS,IBB,JMIN,JMAX,NAME,SLABL,BANDNAME
c=======================================================================
cc    INCLUDE 'BLKPOT.h'
c=======================================================================
c** Effective adiabatic radial potential variables.
      INTEGER BOBCN(NSTATEMX),PSEL(NSTATEMX),MAXMIN(NSTATEMX),
     1 IOMEG(NSTATEMX),Nbeta(NSTATEMX),APSE(NSTATEMX),IFXDE(NSTATEMX),
     2 IFXRE(NSTATEMX),IFXCm(NCMMax,NSTATEMX),
     3 IFXBETA(0:NbetaMX,NSTATEMX),NDATPT(NSTATEMX),NCMM(NSTATEMX),
     4 MMLR(NCMMax,NSTATEMX),nPB(NSTATEMX),nQB(NSTATEMX),pAD(NSTATEMX),
     5 qAD(NSTATEMX),LRad(NSTATEMX),pNA(NSTATEMX),qNA(NSTATEMX),
     6 Pqw(NSTATEMX),IVSR(NSTATEMX),IDSTT(NSTATEMX)
c
      REAL*8 DE(NSTATEMX),RE(NSTATEMX),BETA(0:NbetaMX,NSTATEMX),
     1 yqBETA(NbetaMX,NSTATEMX),BETAFX(NPNTMX,NSTATEMX),RH(NSTATEMX),
     2 RMIN(NSTATEMX),RMAX(NSTATEMX),VLIM(NSTATEMX),EPS(NSTATEMX),
     3 betaINF(NSTATEMX),AGPEF(NSTATEMX),BGPEF(NSTATEMX),
     4 CmVAL(NCMMax,NSTATEMX),CmEFF(NCMMax,NSTATEMX),rhoAB(NSTATEMX),
     5 AA(NSTATEMX),BB(NSTATEMX),RREF(NSTATEMX),ASO(NSTATEMX),
     6 R01(NSTATEMX),Q12(NSTATEMX),RD(NPNTMX,NSTATEMX),
     7 VPOT(NPNTMX,NSTATEMX),dCmA(NCMMax,NSTATEMX),dCmB(NCMMax,NSTATEMX)
c
      COMMON /BLKPOT/DE,RE,BETA,yqBETA,BETAFX,RH,RMIN,RMAX,VLIM,EPS,
     1 betaINF,AGPEF,BGPEF,CmVAL,CmEFF,rhoAB,AA,BB,RREF,ASO,R01,Q12,RD,
     2 VPOT,dCmA,dCmB, BOBCN,PSEL,MAXMIN,IOMEG,Nbeta,APSE,IFXDE,IFXRE,
     3 IFXCm,IFXBETA,NDATPT,NCMM,MMLR,nPB,nQB,pAD,qAD,LRad,pNA,qNA,Pqw,
     4 IVSR,IDSTT
c=======================================================================
cc    INCLUDE 'BLKPARAM.h'
c=======================================================================
c** Parameters and count-labels for band constant (PSEL=-1) or term
c   value (PSEL=-2) fits
      REAL*8 TVALUE(NPARMX),ZBC(0:NVIBMX,0:NROTMX,NISTPMX,NSTATEMX),
     1 ZQC(0:NVIBMX,0:NROTMX,NISTPMX,NSTATEMX)
c
      INTEGER NSTATES,NTVALL(0:NSTATEMX),NTVI(NSTATEMX),NTVF(NSTATEMX),
     1 VMIN(NSTATEMX,NISTPMX),VMAX(NSTATEMX,NISTPMX),JTRUNC(NSTATEMX),
     2 EFSEL(NSTATEMX),NBC(0:NVIBMX,NISTPMX,NSTATEMX),
     3 NQC(0:NVIBMX,NISTPMX,NSTATEMX),
     4 BCPARI(0:NVIBMX,NISTPMX,NSTATEMX),
     5 BCPARF(0:NVIBMX,NISTPMX,NSTATEMX),
     6 QCPARI(0:NVIBMX,NISTPMX,NSTATEMX),
     7 QCPARF(0:NVIBMX,NISTPMX,NSTATEMX)
      COMMON /BLKPARAM/TVALUE,ZBC,ZQC,NSTATES,NTVALL,NTVI,NTVF,VMIN,
     1      VMAX,JTRUNC,EFSEL,NBC,NQC,BCPARI,BCPARF,QCPARI,QCPARF
c=======================================================================
cc    INCLUDE 'BLKDVDP.h'
c=======================================================================
c** Partial derivative arrays for fits and uncertainties (fununc)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      REAL*8 DVtot(HPARMX,NPNTMX),DLDDRe(NPNTMX,NSTATEMX),
     1  DUADRe(NPNTMX,NSTATEMX),DUBDRe(NPNTMX,NSTATEMX),
     2  DTADRe(NPNTMX,NSTATEMX),DTBDRe(NPNTMX,NSTATEMX),
     3  DBDB(0:NbetaMX,NPNTMX,NSTATEMX),DBDRe(NPNTMX,NSTATEMX),
     4  dVpdP(HPARMX,NPNTMX)
      COMMON/BLKDVDP/DVtot,DUADRe,DUBDRe,DTADRe,DTBDRe,DLDDRe,DBDB,
     1 DBDRe,dVpdP
c=======================================================================
cc    INCLUDE 'BLKBOB.h'
c=======================================================================
c** Born-Oppenheimer Breakdown & doubling function parameters.
c**                       March 16 2012
c=======================================================================
      INTEGER NUA(NSTATEMX),NUB(NSTATEMX),NTA(NSTATEMX),NTB(NSTATEMX),
     1  IFXUA(0:NBOBMX,NSTATEMX),IFXUB(0:NBOBMX,NSTATEMX),
     2  IFXTA(0:NBOBMX,NSTATEMX),IFXTB(0:NBOBMX,NSTATEMX),
     3  NwCFT(NSTATEMX),IFXwCFT(0:NBOBMX,NSTATEMX),efREF(NSTATEMX)
c
      REAL*8 UA(0:NBOBMX,NSTATEMX),UB(0:NBOBMX,NSTATEMX),
     1  TA(0:NBOBMX,NSTATEMX),TB(0:NBOBMX,NSTATEMX),
     2   wCFT(0:NBOBMX,NSTATEMX)
c
      COMMON /BLKBOB/UA,UB,TA,TB,wCFT,NUA,NUB,NTA,NTB,NwCFT,
     1  IFXUA,IFXUB,IFXTA,IFXTB,IFXwCFT,efREF
c=======================================================================
cc    INCLUDE 'BLKBOBRF.h'
c=======================================================================
c** Born-Oppenheimer breakdown radial functions 
      REAL*8 UAR(NPNTMX,NSTATEMX),UBR(NPNTMX,NSTATEMX),
     1 TAR(NPNTMX,NSTATEMX),TBR(NPNTMX,NSTATEMX),wRAD(NPNTMX,NSTATEMX)
c
      COMMON /BLKBOBRF/UAR,UBR,TAR,TBR,wRAD
c=======================================================================
cc    INCLUDE 'BLKCOUNT.h'
c=======================================================================
c    Block data file  BLKCOUNT.h
c=======================================================================
c** Counters for numbers of potential parameters of different types for 
c   each state
      INTEGER  TOTPOTPAR,POTPARI(NSTATEMX),POTPARF(NSTATEMX),
     1  UAPARI(NSTATEMX),UAPARF(NSTATEMX),UBPARI(NSTATEMX),
     2  UBPARF(NSTATEMX),TAPARI(NSTATEMX),TAPARF(NSTATEMX),
     3  TBPARI(NSTATEMX),TBPARF(NSTATEMX),LDPARI(NSTATEMX),
     4  LDPARF(NSTATEMX),HPARF(NSTATEMX),OSEL(NSTATEMX)
c
      COMMON /BLKCOUNT/TOTPOTPAR,POTPARI,POTPARF,UAPARI,UAPARF,UBPARI,
     1  UBPARF,TAPARI,TAPARF,TBPARI,TBPARF,LDPARI,LDPARF,HPARF,OSEL
c=======================================================================
c=======================================================================
      INTEGER IDAT, bandN, efPARITY, JRe, fcount
c
      INTEGER I,ICOR,ISTATE,IISTP,INNER,J,JROT,JIN,KVLEV,KV, 
     1  NBEG,NEND,INNODE, IWR, LPRWF
      REAL*8 ZMU,EO,BFCT,JFCT,JFCTP,JFCTL,FWHM,UMAX, RM2, ETRY,VMAXX,
     1  DGDV2,SWF2,JFCTA,JFCTB,DUARe,DUBRe,DTARe,DTBRe,DLDRe,DEROT,
     2  DEROTB,JFCTDBL,JFCTD,muFCT,C3gu,T3,GV(0:NVIBMX),Vtotal(NPNTMX),
     3  SWF(NPNTMX), V(NPNTMX), DEDPK(HPARMX)
c
      COMMON /VBLIK/Vtotal
c
c** Vibrational Band Constants for generating trial energies.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c** Initializing the arrays and variables
      DATA IWR/0/,LPRWF/0/,INNODE/1/
c
c** To calculate the values for <vj|dV/dP(k)|vj>
c   we must first determine the wave equation for the vj state.
      muFCT= 16.857629206d0/ZMU
      BFCT= RH(ISTATE)*RH(ISTATE)/muFCT
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
              IF(IOMEG(ISTATE).EQ.-1) JFCTP= JFCTP - 2**I
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
     3                         + ZMTB(IISTP,ISTATE)*TBR(I,ISTATE))*RM2
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
            CALL SCECOR(KV,KVLEV,JROT,INNER,ICOR,IWR,EO,RH(ISTATE),BFCT,
     1       NDATPT(ISTATE),MMLR(1,ISTATE),V,VMAXX,VLIM(ISTATE),DGDV2)
c***********************************************************************
              KV= KVLEV
              GOTO 10
              ENDIF 
c** If the calculated wavefunction is still for the wrong vibrational
c   level, then write out a warning and skip the calculation of the
c   partial derivatives (hence setting them to zero).
          fcount= fcount+1
          WRITE(6,610) fcount,KVLEV,JROT,KV
          IF(fcount.ge.1000) THEN
              WRITE(6,612)
  612 FORMAT(/' *** Excessive SCECOR failures, so stop and figure out wh
     1y *** ' //)
              STOP
              ENDIF        
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
  610 FORMAT(' *** SCECOR failed',I4,' times,   Currently Seeking   v=',
     1   i3,', J=',i3,';  Found  v=',I3)
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
cc    INCLUDE 'arrsizes.h'
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c** 'Block' Data Utility routine named: 'arrsizes.h' that governs 
c    array dimensioning in program  dPotFit
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER NISTPMX,NPARMX,NbetaMX,NBOBMX,HPARMX,NDATAMX,
     1  NVIBMX,NBCMX,NSTATEMX,NPNTMX,NROTMX,NCMMAX
c*  NISTPMX  is the maximum number of isotopomers allowed for fit
      PARAMETER (NISTPMX =  12)
c*  NSTATEMX  is maximum no. of electronic states which can be
c             simultaneously fitted to
      PARAMETER (NSTATEMX = 4)
c*  NPARMX  is the largest number of free parameters allowed for fit
c  Since FS origins may be parameters, this is also max. no, data bands
      PARAMETER (NPARMX  = 8000)
c*  NbetaMX  is the largest number of exponent parameters allowed for fit
      PARAMETER (NbetaMX  = 40)
c*  NBOBMX-1  is the highest-order polynomial expansion allowed for the
c      adiabatic or centrifugal Born-Oppenheimer breakdown functions, or 
c      the Lambda-doubling or 2\Sigma splitting radial strength functions
      PARAMETER (NBOBMX  = 15)
c*  HPARMX  is the largest number of Hamiltonian parameters of all types
c    (potential energy, BOB. etc.) for all states.
c           HPARMX >= NSTATEMX*[5 + (NbetaMX+1) + 5*(NBOBMX+1)]
      PARAMETER (HPARMX= NSTATEMX*(5 + (NbetaMX+1) + 5*(NBOBMX+1)))
cc    PARAMETER (HPARMX = 300)
c*  NDATAMX  is largest No. of individual data which may be considered
      PARAMETER (NDATAMX = 35000)
c*  NVIBMX  is the maximum number of vibrational levels of a single
c           state for which data are to be considered
      PARAMETER (NVIBMX    = 200)
** NBCMX  is the maximum number of band constants per vib level to be
c         allowed when doing band constant fits (PSEL= -1)
      PARAMETER (NBCMX = 8) 
c*  NPNTMX  is the largest number of potential data points that can be
c           stored in a single 1D array
      PARAMETER (NPNTMX = 90000)
c*  NROTMX  is the highest order of rotational constants calculated and
c            used for estimating level energies
      PARAMETER (NROTMX = 7)
c*  NCMMAX is the largest number of Cm terms in the MLR or DELR 
c            long-range potential
      PARAMETER (NCMMAX = 12)
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cc    INCLUDE 'BLKCOUNT.h'
c=======================================================================
c    Block data file  BLKCOUNT.h
c=======================================================================
c** Counters for numbers of potential parameters of different types for 
c   each state
      INTEGER  TOTPOTPAR,POTPARI(NSTATEMX),POTPARF(NSTATEMX),
     1  UAPARI(NSTATEMX),UAPARF(NSTATEMX),UBPARI(NSTATEMX),
     2  UBPARF(NSTATEMX),TAPARI(NSTATEMX),TAPARF(NSTATEMX),
     3  TBPARI(NSTATEMX),TBPARF(NSTATEMX),LDPARI(NSTATEMX),
     4  LDPARF(NSTATEMX),HPARF(NSTATEMX),OSEL(NSTATEMX)
c
      COMMON /BLKCOUNT/TOTPOTPAR,POTPARI,POTPARF,UAPARI,UAPARF,UBPARI,
     1  UBPARF,TAPARI,TAPARF,TBPARI,TBPARF,LDPARI,LDPARF,HPARF,OSEL
c=======================================================================
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
cc    INCLUDE 'arrsizes.h'
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c** 'Block' Data Utility routine named: 'arrsizes.h' that governs 
c    array dimensioning in program  dPotFit
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER NISTPMX,NPARMX,NbetaMX,NBOBMX,HPARMX,NDATAMX,
     1  NVIBMX,NBCMX,NSTATEMX,NPNTMX,NROTMX,NCMMAX
c*  NISTPMX  is the maximum number of isotopomers allowed for fit
      PARAMETER (NISTPMX =  12)
c*  NSTATEMX  is maximum no. of electronic states which can be
c             simultaneously fitted to
      PARAMETER (NSTATEMX = 4)
c*  NPARMX  is the largest number of free parameters allowed for fit
c  Since FS origins may be parameters, this is also max. no, data bands
      PARAMETER (NPARMX  = 8000)
c*  NbetaMX  is the largest number of exponent parameters allowed for fit
      PARAMETER (NbetaMX  = 40)
c*  NBOBMX-1  is the highest-order polynomial expansion allowed for the
c      adiabatic or centrifugal Born-Oppenheimer breakdown functions, or 
c      the Lambda-doubling or 2\Sigma splitting radial strength functions
      PARAMETER (NBOBMX  = 15)
c*  HPARMX  is the largest number of Hamiltonian parameters of all types
c    (potential energy, BOB. etc.) for all states.
c           HPARMX >= NSTATEMX*[5 + (NbetaMX+1) + 5*(NBOBMX+1)]
      PARAMETER (HPARMX= NSTATEMX*(5 + (NbetaMX+1) + 5*(NBOBMX+1)))
cc    PARAMETER (HPARMX = 300)
c*  NDATAMX  is largest No. of individual data which may be considered
      PARAMETER (NDATAMX = 35000)
c*  NVIBMX  is the maximum number of vibrational levels of a single
c           state for which data are to be considered
      PARAMETER (NVIBMX    = 200)
** NBCMX  is the maximum number of band constants per vib level to be
c         allowed when doing band constant fits (PSEL= -1)
      PARAMETER (NBCMX = 8) 
c*  NPNTMX  is the largest number of potential data points that can be
c           stored in a single 1D array
      PARAMETER (NPNTMX = 90000)
c*  NROTMX  is the highest order of rotational constants calculated and
c            used for estimating level energies
      PARAMETER (NROTMX = 7)
c*  NCMMAX is the largest number of Cm terms in the MLR or DELR 
c            long-range potential
      PARAMETER (NCMMAX = 12)
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cc    INCLUDE 'BLKISOT.h'
c=======================================================================
c** Isotope/isotopologue numbers, masses & BOB mass scaling factors
c** Array ZK carries about the band constants for all levels of all ISOT
      INTEGER NISTP,NDUNMX,AN(2),MN(2,NISTPMX)
c** NDUNMX is a dummy parameter reqd. for portability of READATA
      PARAMETER (NDUNMX=0)
      REAL*8  ZMASS(3,NISTPMX),RSQMU(NISTPMX),RSQMUP(0:NDUNMX,NISTPMX),
     1 RMUP(0:9,NISTPMX),ZMUA(NISTPMX,NSTATEMX),ZMUB(NISTPMX,NSTATEMX),
     2 ZMTA(NISTPMX,NSTATEMX),ZMTB(NISTPMX,NSTATEMX),
     3  ZK(0:NVIBMX,0:NROTMX,NISTPMX,NSTATEMX)
c
      COMMON /BLKISOT/ZMASS,RSQMU,RSQMUP,RMUP,ZMUA,ZMUB,ZMTA,ZMTB,ZK,
     1  NISTP,AN,MN
c=======================================================================
cc    INCLUDE 'BLKDATA.h'
c=======================================================================
c** Type statements & common block for data
      REAL*8  FREQ(NDATAMX),UFREQ(NDATAMX),DFREQ(NDATAMX),TEMP(NDATAMX),
     1                                               YUNC(NDATAMX),Fqb
      INTEGER  COUNTOT,NFS1,NFSTOT,NBANDTOT,IB(NDATAMX),JP(NDATAMX),
     1 JPP(NDATAMX),VP(NPARMX),VPP(NPARMX),EFP(NDATAMX),EFPP(NDATAMX),
     2 TVUP(NDATAMX),TVLW(NDATAMX),FSBAND(NPARMX),IFXFS(NPARMX),
     3 NFS(NPARMX),IEP(NPARMX),IEPP(NPARMX),ISTP(NPARMX),
     4 IFIRST(NPARMX),ILAST(NPARMX),NTV(NSTATEMX,NISTPMX),FSsame,
     5 NTRANS(NPARMX),IBB(NISTPMX,NSTATEMX,9,NPARMX),JMIN(NPARMX),
     6 JMAX(NPARMX)
      CHARACTER*2 NAME(2) 
      CHARACTER*3 SLABL(-6:NSTATEMX)
      CHARACTER*30 BANDNAME(NPARMX)
      COMMON /DATABLK/Fqb,FREQ,UFREQ,YUNC,DFREQ,TEMP,COUNTOT,NFS1,
     1 NFSTOT,NBANDTOT,IB,JP,JPP,VP,VPP,EFP,EFPP,TVUP,TVLW,FSBAND,IFXFS,
     2 NFS,IEP,IEPP,ISTP,IFIRST,ILAST,NTV,FSsame,
     3 NTRANS,IBB,JMIN,JMAX,NAME,SLABL,BANDNAME
c=======================================================================
cc    INCLUDE 'BLKPOT.h'
c=======================================================================
c** Effective adiabatic radial potential variables.
      INTEGER BOBCN(NSTATEMX),PSEL(NSTATEMX),MAXMIN(NSTATEMX),
     1 IOMEG(NSTATEMX),Nbeta(NSTATEMX),APSE(NSTATEMX),IFXDE(NSTATEMX),
     2 IFXRE(NSTATEMX),IFXCm(NCMMax,NSTATEMX),
     3 IFXBETA(0:NbetaMX,NSTATEMX),NDATPT(NSTATEMX),NCMM(NSTATEMX),
     4 MMLR(NCMMax,NSTATEMX),nPB(NSTATEMX),nQB(NSTATEMX),pAD(NSTATEMX),
     5 qAD(NSTATEMX),LRad(NSTATEMX),pNA(NSTATEMX),qNA(NSTATEMX),
     6 Pqw(NSTATEMX),IVSR(NSTATEMX),IDSTT(NSTATEMX)
c
      REAL*8 DE(NSTATEMX),RE(NSTATEMX),BETA(0:NbetaMX,NSTATEMX),
     1 yqBETA(NbetaMX,NSTATEMX),BETAFX(NPNTMX,NSTATEMX),RH(NSTATEMX),
     2 RMIN(NSTATEMX),RMAX(NSTATEMX),VLIM(NSTATEMX),EPS(NSTATEMX),
     3 betaINF(NSTATEMX),AGPEF(NSTATEMX),BGPEF(NSTATEMX),
     4 CmVAL(NCMMax,NSTATEMX),CmEFF(NCMMax,NSTATEMX),rhoAB(NSTATEMX),
     5 AA(NSTATEMX),BB(NSTATEMX),RREF(NSTATEMX),ASO(NSTATEMX),
     6 R01(NSTATEMX),Q12(NSTATEMX),RD(NPNTMX,NSTATEMX),
     7 VPOT(NPNTMX,NSTATEMX),dCmA(NCMMax,NSTATEMX),dCmB(NCMMax,NSTATEMX)
c
      COMMON /BLKPOT/DE,RE,BETA,yqBETA,BETAFX,RH,RMIN,RMAX,VLIM,EPS,
     1 betaINF,AGPEF,BGPEF,CmVAL,CmEFF,rhoAB,AA,BB,RREF,ASO,R01,Q12,RD,
     2 VPOT,dCmA,dCmB, BOBCN,PSEL,MAXMIN,IOMEG,Nbeta,APSE,IFXDE,IFXRE,
     3 IFXCm,IFXBETA,NDATPT,NCMM,MMLR,nPB,nQB,pAD,qAD,LRad,pNA,qNA,Pqw,
     4 IVSR,IDSTT
c=======================================================================
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
cc    INCLUDE 'arrsizes.h'
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c** 'Block' Data Utility routine named: 'arrsizes.h' that governs 
c    array dimensioning in program  dPotFit
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER NISTPMX,NPARMX,NbetaMX,NBOBMX,HPARMX,NDATAMX,
     1  NVIBMX,NBCMX,NSTATEMX,NPNTMX,NROTMX,NCMMAX
c*  NISTPMX  is the maximum number of isotopomers allowed for fit
      PARAMETER (NISTPMX =  12)
c*  NSTATEMX  is maximum no. of electronic states which can be
c             simultaneously fitted to
      PARAMETER (NSTATEMX = 4)
c*  NPARMX  is the largest number of free parameters allowed for fit
c  Since FS origins may be parameters, this is also max. no, data bands
      PARAMETER (NPARMX  = 8000)
c*  NbetaMX  is the largest number of exponent parameters allowed for fit
      PARAMETER (NbetaMX  = 40)
c*  NBOBMX-1  is the highest-order polynomial expansion allowed for the
c      adiabatic or centrifugal Born-Oppenheimer breakdown functions, or 
c      the Lambda-doubling or 2\Sigma splitting radial strength functions
      PARAMETER (NBOBMX  = 15)
c*  HPARMX  is the largest number of Hamiltonian parameters of all types
c    (potential energy, BOB. etc.) for all states.
c           HPARMX >= NSTATEMX*[5 + (NbetaMX+1) + 5*(NBOBMX+1)]
      PARAMETER (HPARMX= NSTATEMX*(5 + (NbetaMX+1) + 5*(NBOBMX+1)))
cc    PARAMETER (HPARMX = 300)
c*  NDATAMX  is largest No. of individual data which may be considered
      PARAMETER (NDATAMX = 35000)
c*  NVIBMX  is the maximum number of vibrational levels of a single
c           state for which data are to be considered
      PARAMETER (NVIBMX    = 200)
** NBCMX  is the maximum number of band constants per vib level to be
c         allowed when doing band constant fits (PSEL= -1)
      PARAMETER (NBCMX = 8) 
c*  NPNTMX  is the largest number of potential data points that can be
c           stored in a single 1D array
      PARAMETER (NPNTMX = 90000)
c*  NROTMX  is the highest order of rotational constants calculated and
c            used for estimating level energies
      PARAMETER (NROTMX = 7)
c*  NCMMAX is the largest number of Cm terms in the MLR or DELR 
c            long-range potential
      PARAMETER (NCMMAX = 12)
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cc    INCLUDE 'BLKCOUNT.h'
c=======================================================================
c    Block data file  BLKCOUNT.h
c=======================================================================
c** Counters for numbers of potential parameters of different types for 
c   each state
      INTEGER  TOTPOTPAR,POTPARI(NSTATEMX),POTPARF(NSTATEMX),
     1  UAPARI(NSTATEMX),UAPARF(NSTATEMX),UBPARI(NSTATEMX),
     2  UBPARF(NSTATEMX),TAPARI(NSTATEMX),TAPARF(NSTATEMX),
     3  TBPARI(NSTATEMX),TBPARF(NSTATEMX),LDPARI(NSTATEMX),
     4  LDPARF(NSTATEMX),HPARF(NSTATEMX),OSEL(NSTATEMX)
c
      COMMON /BLKCOUNT/TOTPOTPAR,POTPARI,POTPARF,UAPARI,UAPARF,UBPARI,
     1  UBPARF,TAPARI,TAPARF,TBPARI,TBPARF,LDPARI,LDPARF,HPARF,OSEL
c=======================================================================
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

c***********************************************************************
      SUBROUTINE DVIRDP(IDAT,ISTATE,ZMU,BVIR,dBVIRdP)
c=======================================================================
c  This subroutine calculates the second virial coefficient BVIR and its
c  partial derivatives with respect to the various parameters. It is 
c  used when virial data has been input into the program.  It performs a
c  classical calcultion with two quantum corrections.   Update: 15/05/16
c=======================================================================
cc    INCLUDE 'arrsizes.h'
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c** 'Block' Data Utility routine named: 'arrsizes.h' that governs 
c    array dimensioning in program  dPotFit
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER NISTPMX,NPARMX,NbetaMX,NBOBMX,HPARMX,NDATAMX,
     1  NVIBMX,NBCMX,NSTATEMX,NPNTMX,NROTMX,NCMMAX
c*  NISTPMX  is the maximum number of isotopomers allowed for fit
      PARAMETER (NISTPMX =  12)
c*  NSTATEMX  is maximum no. of electronic states which can be
c             simultaneously fitted to
      PARAMETER (NSTATEMX = 4)
c*  NPARMX  is the largest number of free parameters allowed for fit
c  Since FS origins may be parameters, this is also max. no, data bands
      PARAMETER (NPARMX  = 8000)
c*  NbetaMX  is the largest number of exponent parameters allowed for fit
      PARAMETER (NbetaMX  = 40)
c*  NBOBMX-1  is the highest-order polynomial expansion allowed for the
c      adiabatic or centrifugal Born-Oppenheimer breakdown functions, or 
c      the Lambda-doubling or 2\Sigma splitting radial strength functions
      PARAMETER (NBOBMX  = 15)
c*  HPARMX  is the largest number of Hamiltonian parameters of all types
c    (potential energy, BOB. etc.) for all states.
c           HPARMX >= NSTATEMX*[5 + (NbetaMX+1) + 5*(NBOBMX+1)]
      PARAMETER (HPARMX= NSTATEMX*(5 + (NbetaMX+1) + 5*(NBOBMX+1)))
cc    PARAMETER (HPARMX = 300)
c*  NDATAMX  is largest No. of individual data which may be considered
      PARAMETER (NDATAMX = 35000)
c*  NVIBMX  is the maximum number of vibrational levels of a single
c           state for which data are to be considered
      PARAMETER (NVIBMX    = 200)
** NBCMX  is the maximum number of band constants per vib level to be
c         allowed when doing band constant fits (PSEL= -1)
      PARAMETER (NBCMX = 8) 
c*  NPNTMX  is the largest number of potential data points that can be
c           stored in a single 1D array
      PARAMETER (NPNTMX = 90000)
c*  NROTMX  is the highest order of rotational constants calculated and
c            used for estimating level energies
      PARAMETER (NROTMX = 7)
c*  NCMMAX is the largest number of Cm terms in the MLR or DELR 
c            long-range potential
      PARAMETER (NCMMAX = 12)
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cc    INCLUDE 'BLKPOT.h'
c=======================================================================
c** Effective adiabatic radial potential variables.
      INTEGER BOBCN(NSTATEMX),PSEL(NSTATEMX),MAXMIN(NSTATEMX),
     1 IOMEG(NSTATEMX),Nbeta(NSTATEMX),APSE(NSTATEMX),IFXDE(NSTATEMX),
     2 IFXRE(NSTATEMX),IFXCm(NCMMax,NSTATEMX),
     3 IFXBETA(0:NbetaMX,NSTATEMX),NDATPT(NSTATEMX),NCMM(NSTATEMX),
     4 MMLR(NCMMax,NSTATEMX),nPB(NSTATEMX),nQB(NSTATEMX),pAD(NSTATEMX),
     5 qAD(NSTATEMX),LRad(NSTATEMX),pNA(NSTATEMX),qNA(NSTATEMX),
     6 Pqw(NSTATEMX),IVSR(NSTATEMX),IDSTT(NSTATEMX)
c
      REAL*8 DE(NSTATEMX),RE(NSTATEMX),BETA(0:NbetaMX,NSTATEMX),
     1 yqBETA(NbetaMX,NSTATEMX),BETAFX(NPNTMX,NSTATEMX),RH(NSTATEMX),
     2 RMIN(NSTATEMX),RMAX(NSTATEMX),VLIM(NSTATEMX),EPS(NSTATEMX),
     3 betaINF(NSTATEMX),AGPEF(NSTATEMX),BGPEF(NSTATEMX),
     4 CmVAL(NCMMax,NSTATEMX),CmEFF(NCMMax,NSTATEMX),rhoAB(NSTATEMX),
     5 AA(NSTATEMX),BB(NSTATEMX),RREF(NSTATEMX),ASO(NSTATEMX),
     6 R01(NSTATEMX),Q12(NSTATEMX),RD(NPNTMX,NSTATEMX),
     7 VPOT(NPNTMX,NSTATEMX),dCmA(NCMMax,NSTATEMX),dCmB(NCMMax,NSTATEMX)
c
      COMMON /BLKPOT/DE,RE,BETA,yqBETA,BETAFX,RH,RMIN,RMAX,VLIM,EPS,
     1 betaINF,AGPEF,BGPEF,CmVAL,CmEFF,rhoAB,AA,BB,RREF,ASO,R01,Q12,RD,
     2 VPOT,dCmA,dCmB, BOBCN,PSEL,MAXMIN,IOMEG,Nbeta,APSE,IFXDE,IFXRE,
     3 IFXCm,IFXBETA,NDATPT,NCMM,MMLR,nPB,nQB,pAD,qAD,LRad,pNA,qNA,Pqw,
     4 IVSR,IDSTT
c=======================================================================
cc    INCLUDE 'BLKCOUNT.h'
c=======================================================================
c    Block data file  BLKCOUNT.h
c=======================================================================
c** Counters for numbers of potential parameters of different types for 
c   each state
      INTEGER  TOTPOTPAR,POTPARI(NSTATEMX),POTPARF(NSTATEMX),
     1  UAPARI(NSTATEMX),UAPARF(NSTATEMX),UBPARI(NSTATEMX),
     2  UBPARF(NSTATEMX),TAPARI(NSTATEMX),TAPARF(NSTATEMX),
     3  TBPARI(NSTATEMX),TBPARF(NSTATEMX),LDPARI(NSTATEMX),
     4  LDPARF(NSTATEMX),HPARF(NSTATEMX),OSEL(NSTATEMX)
c
      COMMON /BLKCOUNT/TOTPOTPAR,POTPARI,POTPARF,UAPARI,UAPARF,UBPARI,
     1  UBPARF,TAPARI,TAPARF,TBPARI,TBPARF,LDPARI,LDPARF,HPARF,OSEL
c=======================================================================
cc    INCLUDE 'BLKDVDP.h'
c=======================================================================
c** Partial derivative arrays for fits and uncertainties (fununc)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      REAL*8 DVtot(HPARMX,NPNTMX),DLDDRe(NPNTMX,NSTATEMX),
     1  DUADRe(NPNTMX,NSTATEMX),DUBDRe(NPNTMX,NSTATEMX),
     2  DTADRe(NPNTMX,NSTATEMX),DTBDRe(NPNTMX,NSTATEMX),
     3  DBDB(0:NbetaMX,NPNTMX,NSTATEMX),DBDRe(NPNTMX,NSTATEMX),
     4  dVpdP(HPARMX,NPNTMX)
      COMMON/BLKDVDP/DVtot,DUADRe,DUBDRe,DTADRe,DTBDRe,DLDDRe,DBDB,
     1 DBDRe,dVpdP
c=======================================================================
cc    INCLUDE 'BLKISOT.h'
c=======================================================================
c** Isotope/isotopologue numbers, masses & BOB mass scaling factors
c** Array ZK carries about the band constants for all levels of all ISOT
      INTEGER NISTP,NDUNMX,AN(2),MN(2,NISTPMX)
c** NDUNMX is a dummy parameter reqd. for portability of READATA
      PARAMETER (NDUNMX=0)
      REAL*8  ZMASS(3,NISTPMX),RSQMU(NISTPMX),RSQMUP(0:NDUNMX,NISTPMX),
     1 RMUP(0:9,NISTPMX),ZMUA(NISTPMX,NSTATEMX),ZMUB(NISTPMX,NSTATEMX),
     2 ZMTA(NISTPMX,NSTATEMX),ZMTB(NISTPMX,NSTATEMX),
     3  ZK(0:NVIBMX,0:NROTMX,NISTPMX,NSTATEMX)
c
      COMMON /BLKISOT/ZMASS,RSQMU,RSQMUP,RMUP,ZMUA,ZMUB,ZMTA,ZMTB,ZK,
     1  NISTP,AN,MN
c=======================================================================
cc    INCLUDE 'BLKDATA.h'
c=======================================================================
c** Type statements & common block for data
      REAL*8  FREQ(NDATAMX),UFREQ(NDATAMX),DFREQ(NDATAMX),TEMP(NDATAMX),
     1                                               YUNC(NDATAMX),Fqb
      INTEGER  COUNTOT,NFS1,NFSTOT,NBANDTOT,IB(NDATAMX),JP(NDATAMX),
     1 JPP(NDATAMX),VP(NPARMX),VPP(NPARMX),EFP(NDATAMX),EFPP(NDATAMX),
     2 TVUP(NDATAMX),TVLW(NDATAMX),FSBAND(NPARMX),IFXFS(NPARMX),
     3 NFS(NPARMX),IEP(NPARMX),IEPP(NPARMX),ISTP(NPARMX),
     4 IFIRST(NPARMX),ILAST(NPARMX),NTV(NSTATEMX,NISTPMX),FSsame,
     5 NTRANS(NPARMX),IBB(NISTPMX,NSTATEMX,9,NPARMX),JMIN(NPARMX),
     6 JMAX(NPARMX)
      CHARACTER*2 NAME(2) 
      CHARACTER*3 SLABL(-6:NSTATEMX)
      CHARACTER*30 BANDNAME(NPARMX)
      COMMON /DATABLK/Fqb,FREQ,UFREQ,YUNC,DFREQ,TEMP,COUNTOT,NFS1,
     1 NFSTOT,NBANDTOT,IB,JP,JPP,VP,VPP,EFP,EFPP,TVUP,TVLW,FSBAND,IFXFS,
     2 NFS,IEP,IEPP,ISTP,IFIRST,ILAST,NTV,FSsame,
     3 NTRANS,IBB,JMIN,JMAX,NAME,SLABL,BANDNAME
c=======================================================================
c* Define types for local variables
      INTEGER check,j,counter,nParams,ISTATE,ISTART,IDAT,i,m,k,kk,jj,n,
     1  VIRCNT(NSTATEMX) 
      REAL*8  XG(8),WG(8),x(4),w(4),BTemp,BTempInv,
     1 CONST(3),jump,a,b,YVAL(8),RDIST(8),Class,Q1corr,Q2corr,
     2 VDIST(8),EXP_TERM,int_fact,Vsq,VPsq,Rsq,RINV,XTEMP,Bclass,Bq1,
     3 Bq2,INTEGRALS(NPARMX+3),error,BVIR,dVdR(8),d2VdR(8),
     4 XTEMP1,dBcdP(NPARMX),dBq1dP(NPARMX),dBVIRdP(NPARMX),ZMU,dBVIR
c
       REAL*8, PARAMETER :: k_boltz=1.3806488D-23, NA=6.02214129D23,
     1  h = 6.62606957D-34, c = 2.99792458D10, k_cm = k_boltz/(h*c),
     2  h_cm = 1.d0/c, pi = 3.14159265358979323846, amu=1.660538921D-27,
     3  Cu= 16.857629206D0
c* Common block data for quadrature weights and points
      data x/0.960289856497536d0, 0.796666477413627d0,
     1       0.525532409916329d0, 0.183434642495650d0/,
     2     w/0.101228536290376d0, 0.222381034453374d0,
     3       0.313706645877887d0, 0.362683783378362d0/
c***********************************************************************
      ERROR= 0.001d0
      IF(VIRCNT(ISTATE).EQ.0) THEN
          VIRCNT(ISTATE)= VIRCNT(ISTATE) + 1
          WRITE(32,100)
          ENDIF
      ISTART = POTPARI(ISTATE)- 1
c.. runs a loop to set all quadrature weights and points
      DO j= 1,4
          XG(j)= -x(j)
          WG(j)= w(j)
          XG(9-j)= x(j)
          WG(9-j)= w(j)
          ENDDO
c.. initializes the array dBVIRdP to zero.
      DO J= 1,HPARMX+3
          dBVIRdP(J)= 0.d0
          ENDDO
c.. takes the current temperature from the virial data and changes it to E cm
      BTemp = k_cm*Temp(IDAT)
c.. inverts BTemp
      BTempInv = 1.d0/BTemp
c.. sets the constants for the intergration of each correction term
      Const(1)= -2.d0*pi*0.6022140857d0           !! updated Avogadro No.
      Const(2)= pi*0.6022140857d0*Cu*BTempInv**3/(6.d0*ZMU)
      Const(3)= -(pi*0.6022140857d0/6.d0)*(Cu*BTempInv**2/ZMU)**2
c.. sets the number of paramters used for the integration and initializes
c.  the integrals to zero
      nParams = 3 + NCMM(ISTATE) + Nbeta(ISTATE)
      check = 1
c.. this first loop is here to repeat the calculations until the values
c** Outermost loop: repeatedly bisect overall [-1,+1] interval NBISMX
c  times until convergence to within absolute error  ERROR  is achieved
      n= 12
      DO kk= 1,n
          DO J= 1,nParams+3
              INTEGRALS(J)= 0.d0
              ENDDO
          IF(check.EQ.-1) EXIT
          counter= 2**kk
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
c... q=1 mapping  r <-> y= (r/re - 1)/(r/re + 1) 
                  RDIST(i)= Re(ISTATE)*DSQRT((1.0d0 + YVAL(i))/
     1                                                (1.d0 - YVAL(i)))
                  VDIST(i)= 0.d0
                  ENDDO
c.. calls VGENP to find the neccesary PEC value and radial derivatives
              CALL VGENP(ISTATE,RDIST,VDIST,dVdR,d2VdR,IDAT)
c.. the next loop over each Gaussian point within each subinterval
              DO  i= 1,8
c.. some of the terms that will be used in later calculations of the
c.  virial coefficients are constructed here
                  IF(RDIST(I).LT.0.8d0*RE(ISTATE)) THEN
c** Check for and correct for potential turnover problems
                      IF((VDIST(I).LT.0.d0).OR.(DVDR(I).GT.0.d0)
c... or for potential becoming singular at very short dist (r < 0.1 Ang)
     1             .OR.(RDIST(I).LE.0.1d0).OR.(D2VDR(I).LT.0.d0)) THEN
                          VDIST(I)= 1.d6
                          dVdR(I)= -1.d6
                          d2VdR(I)= 1.d6
                          ENDIF
                      ENDIF
                  EXP_TERM= DEXP(-VDIST(i)*BTempINV)
                  RINV= 1.d0/RDIST(i)
                  int_fact= (Re(ISTATE)/(1.d0 - YVAL(i)))**2
     1                                     *RINV*(b - a)*0.5d0
                  Vsq= VDIST(i)**2
                  VPsq= dVdR(i)**2
                  Rsq= RDIST(i)**2
                  XTEMP= d2VdR(i)**2/10.d0 + VPsq*RINV**2/5.d0
     1      + dVdR(i)**3*BTempInv*RINV/9.d0 - (VPsq*BTempInv)**2/72.d0
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c*     Finally the relevant partial derivatives of the virial 
c*     coefficients w.r.t. potential param are calculated
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
                  ENDDO     !! end of 'i' loop over 8 gaussian points 
              ENDDO    !! end of 'm' loop over 'counter' radial subintervals
          IF(PSEL(ISTATE).LE.3) THEN
                  DO j= 1,nParams      ! partial derivatives for fit cases
                  dBVIRdP(j)= INTEGRALS(j)
                  ENDDO
              ENDIF
          dBVIR=BVIR
          Class= Const(1)*INTEGRALS(nParams+1)
          Q1corr= Const(2)*INTEGRALS(nParams+2)
          Q2corr= Const(3)*Integrals(nParams+3)
          BVIR= Class + Q1corr+ Q2corr
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c. Check to see if BVIR has converged
          IF(counter.GE.2) THEN
              error=0.000001
              IF(DABS(BVIR-dBVIR).LE.(abs(error*BVIR))) THEN
                  check= -1 
                  ENDIF
              ENDIF 
          ENDDO               !! end of 'kk' loop over # sets of subdivisions
      WRITE(32,102) Temp(IDAT),Class,Q1corr,Q2corr,BVIR,counter
  100 FORMAT('Temperature Classical   FirstQ   SecondQ   Total  counter
     1', /24('==='))
  102 FORMAT(2X,F7.2,2X,F8.3,2X,F8.3,2X,F8.3,2X,F8.3,I6)
      RETURN
      END 
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE DVACDP(IDAT,ISTATE,ZMU,BVIR,dBVIRdP)
c=======================================================================
c  This subroutine calculates the acoustic second virial coefficient BVIR 
c  and its partial derivatives with respect to the various parameters. It is 
c  used when virial data have been input into the program.  It performs a
c  classical calculation plus two quantum corrections
c=======================================================================
cc    INCLUDE 'arrsizes.h'
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c** 'Block' Data Utility routine named: 'arrsizes.h' that governs 
c    array dimensioning in program  dPotFit
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER NISTPMX,NPARMX,NbetaMX,NBOBMX,HPARMX,NDATAMX,
     1  NVIBMX,NBCMX,NSTATEMX,NPNTMX,NROTMX,NCMMAX
c*  NISTPMX  is the maximum number of isotopomers allowed for fit
      PARAMETER (NISTPMX =  12)
c*  NSTATEMX  is maximum no. of electronic states which can be
c             simultaneously fitted to
      PARAMETER (NSTATEMX = 4)
c*  NPARMX  is the largest number of free parameters allowed for fit
c  Since FS origins may be parameters, this is also max. no, data bands
      PARAMETER (NPARMX  = 8000)
c*  NbetaMX  is the largest number of exponent parameters allowed for fit
      PARAMETER (NbetaMX  = 40)
c*  NBOBMX-1  is the highest-order polynomial expansion allowed for the
c      adiabatic or centrifugal Born-Oppenheimer breakdown functions, or 
c      the Lambda-doubling or 2\Sigma splitting radial strength functions
      PARAMETER (NBOBMX  = 15)
c*  HPARMX  is the largest number of Hamiltonian parameters of all types
c    (potential energy, BOB. etc.) for all states.
c           HPARMX >= NSTATEMX*[5 + (NbetaMX+1) + 5*(NBOBMX+1)]
      PARAMETER (HPARMX= NSTATEMX*(5 + (NbetaMX+1) + 5*(NBOBMX+1)))
cc    PARAMETER (HPARMX = 300)
c*  NDATAMX  is largest No. of individual data which may be considered
      PARAMETER (NDATAMX = 35000)
c*  NVIBMX  is the maximum number of vibrational levels of a single
c           state for which data are to be considered
      PARAMETER (NVIBMX    = 200)
** NBCMX  is the maximum number of band constants per vib level to be
c         allowed when doing band constant fits (PSEL= -1)
      PARAMETER (NBCMX = 8) 
c*  NPNTMX  is the largest number of potential data points that can be
c           stored in a single 1D array
      PARAMETER (NPNTMX = 90000)
c*  NROTMX  is the highest order of rotational constants calculated and
c            used for estimating level energies
      PARAMETER (NROTMX = 7)
c*  NCMMAX is the largest number of Cm terms in the MLR or DELR 
c            long-range potential
      PARAMETER (NCMMAX = 12)
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cc    INCLUDE 'BLKPOT.h'
c=======================================================================
c** Effective adiabatic radial potential variables.
      INTEGER BOBCN(NSTATEMX),PSEL(NSTATEMX),MAXMIN(NSTATEMX),
     1 IOMEG(NSTATEMX),Nbeta(NSTATEMX),APSE(NSTATEMX),IFXDE(NSTATEMX),
     2 IFXRE(NSTATEMX),IFXCm(NCMMax,NSTATEMX),
     3 IFXBETA(0:NbetaMX,NSTATEMX),NDATPT(NSTATEMX),NCMM(NSTATEMX),
     4 MMLR(NCMMax,NSTATEMX),nPB(NSTATEMX),nQB(NSTATEMX),pAD(NSTATEMX),
     5 qAD(NSTATEMX),LRad(NSTATEMX),pNA(NSTATEMX),qNA(NSTATEMX),
     6 Pqw(NSTATEMX),IVSR(NSTATEMX),IDSTT(NSTATEMX)
c
      REAL*8 DE(NSTATEMX),RE(NSTATEMX),BETA(0:NbetaMX,NSTATEMX),
     1 yqBETA(NbetaMX,NSTATEMX),BETAFX(NPNTMX,NSTATEMX),RH(NSTATEMX),
     2 RMIN(NSTATEMX),RMAX(NSTATEMX),VLIM(NSTATEMX),EPS(NSTATEMX),
     3 betaINF(NSTATEMX),AGPEF(NSTATEMX),BGPEF(NSTATEMX),
     4 CmVAL(NCMMax,NSTATEMX),CmEFF(NCMMax,NSTATEMX),rhoAB(NSTATEMX),
     5 AA(NSTATEMX),BB(NSTATEMX),RREF(NSTATEMX),ASO(NSTATEMX),
     6 R01(NSTATEMX),Q12(NSTATEMX),RD(NPNTMX,NSTATEMX),
     7 VPOT(NPNTMX,NSTATEMX),dCmA(NCMMax,NSTATEMX),dCmB(NCMMax,NSTATEMX)
c
      COMMON /BLKPOT/DE,RE,BETA,yqBETA,BETAFX,RH,RMIN,RMAX,VLIM,EPS,
     1 betaINF,AGPEF,BGPEF,CmVAL,CmEFF,rhoAB,AA,BB,RREF,ASO,R01,Q12,RD,
     2 VPOT,dCmA,dCmB, BOBCN,PSEL,MAXMIN,IOMEG,Nbeta,APSE,IFXDE,IFXRE,
     3 IFXCm,IFXBETA,NDATPT,NCMM,MMLR,nPB,nQB,pAD,qAD,LRad,pNA,qNA,Pqw,
     4 IVSR,IDSTT
c=======================================================================
cc    INCLUDE 'BLKCOUNT.h'
c=======================================================================
c    Block data file  BLKCOUNT.h
c=======================================================================
c** Counters for numbers of potential parameters of different types for 
c   each state
      INTEGER  TOTPOTPAR,POTPARI(NSTATEMX),POTPARF(NSTATEMX),
     1  UAPARI(NSTATEMX),UAPARF(NSTATEMX),UBPARI(NSTATEMX),
     2  UBPARF(NSTATEMX),TAPARI(NSTATEMX),TAPARF(NSTATEMX),
     3  TBPARI(NSTATEMX),TBPARF(NSTATEMX),LDPARI(NSTATEMX),
     4  LDPARF(NSTATEMX),HPARF(NSTATEMX),OSEL(NSTATEMX)
c
      COMMON /BLKCOUNT/TOTPOTPAR,POTPARI,POTPARF,UAPARI,UAPARF,UBPARI,
     1  UBPARF,TAPARI,TAPARF,TBPARI,TBPARF,LDPARI,LDPARF,HPARF,OSEL
c=======================================================================
cc    INCLUDE 'BLKDVDP.h'
c=======================================================================
c** Partial derivative arrays for fits and uncertainties (fununc)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      REAL*8 DVtot(HPARMX,NPNTMX),DLDDRe(NPNTMX,NSTATEMX),
     1  DUADRe(NPNTMX,NSTATEMX),DUBDRe(NPNTMX,NSTATEMX),
     2  DTADRe(NPNTMX,NSTATEMX),DTBDRe(NPNTMX,NSTATEMX),
     3  DBDB(0:NbetaMX,NPNTMX,NSTATEMX),DBDRe(NPNTMX,NSTATEMX),
     4  dVpdP(HPARMX,NPNTMX)
      COMMON/BLKDVDP/DVtot,DUADRe,DUBDRe,DTADRe,DTBDRe,DLDDRe,DBDB,
     1 DBDRe,dVpdP
c=======================================================================
cc    INCLUDE 'BLKISOT.h'
c=======================================================================
c** Isotope/isotopologue numbers, masses & BOB mass scaling factors
c** Array ZK carries about the band constants for all levels of all ISOT
      INTEGER NISTP,NDUNMX,AN(2),MN(2,NISTPMX)
c** NDUNMX is a dummy parameter reqd. for portability of READATA
      PARAMETER (NDUNMX=0)
      REAL*8  ZMASS(3,NISTPMX),RSQMU(NISTPMX),RSQMUP(0:NDUNMX,NISTPMX),
     1 RMUP(0:9,NISTPMX),ZMUA(NISTPMX,NSTATEMX),ZMUB(NISTPMX,NSTATEMX),
     2 ZMTA(NISTPMX,NSTATEMX),ZMTB(NISTPMX,NSTATEMX),
     3  ZK(0:NVIBMX,0:NROTMX,NISTPMX,NSTATEMX)
c
      COMMON /BLKISOT/ZMASS,RSQMU,RSQMUP,RMUP,ZMUA,ZMUB,ZMTA,ZMTB,ZK,
     1  NISTP,AN,MN
c=======================================================================
cc    INCLUDE 'BLKDATA.h'
c=======================================================================
c** Type statements & common block for data
      REAL*8  FREQ(NDATAMX),UFREQ(NDATAMX),DFREQ(NDATAMX),TEMP(NDATAMX),
     1                                               YUNC(NDATAMX),Fqb
      INTEGER  COUNTOT,NFS1,NFSTOT,NBANDTOT,IB(NDATAMX),JP(NDATAMX),
     1 JPP(NDATAMX),VP(NPARMX),VPP(NPARMX),EFP(NDATAMX),EFPP(NDATAMX),
     2 TVUP(NDATAMX),TVLW(NDATAMX),FSBAND(NPARMX),IFXFS(NPARMX),
     3 NFS(NPARMX),IEP(NPARMX),IEPP(NPARMX),ISTP(NPARMX),
     4 IFIRST(NPARMX),ILAST(NPARMX),NTV(NSTATEMX,NISTPMX),FSsame,
     5 NTRANS(NPARMX),IBB(NISTPMX,NSTATEMX,9,NPARMX),JMIN(NPARMX),
     6 JMAX(NPARMX)
      CHARACTER*2 NAME(2) 
      CHARACTER*3 SLABL(-6:NSTATEMX)
      CHARACTER*30 BANDNAME(NPARMX)
      COMMON /DATABLK/Fqb,FREQ,UFREQ,YUNC,DFREQ,TEMP,COUNTOT,NFS1,
     1 NFSTOT,NBANDTOT,IB,JP,JPP,VP,VPP,EFP,EFPP,TVUP,TVLW,FSBAND,IFXFS,
     2 NFS,IEP,IEPP,ISTP,IFIRST,ILAST,NTV,FSsame,
     3 NTRANS,IBB,JMIN,JMAX,NAME,SLABL,BANDNAME
c=======================================================================
c* Define types for local variables
       INTEGER check,j,counter,nParams,ISTATE,ISTART,IDAT,i,m,k,kk,jj,
     1  n,IPASS,VACCNT(NSTATEMX)
       REAL*8  XG(8),WG(8),x(4),w(4),BTemp,BTempInv,BV,BVP,BVPP,BVsq,
     1 CONST(3),jump,a,b,YVAL(8),RDIST(8),BVPsq,BVPcu,BVPTtF,BVPPsq,
     2 VDIST(8),EXP_TERM,int_fact,Vsq,VPsq,Rsq,RINV,XTEMP,Bclass,Bq1,
     3 Bq2,INTEGRALS(NPARMX+3),error,BVIR,dVdR(8),d2VdR(8),Risq,
     4 XTEMP1,dBcdP(NPARMX),dBq1dP(NPARMX),dBVIRdP(NPARMX),ZMU,dBVIR
     5 ,Ccorr,Qcorr,Q2corr
c
       REAL*8, PARAMETER :: k_boltz=1.3806488D-23, NA=6.02214129D23,
     1  h = 6.62606957D-34, c = 2.99792458D10, k_cm = k_boltz/(h*c),
     2  h_cm = 1.d0/c, pi = 3.14159265358979323846, amu=1.660538921D-27,
     3  Cu= 16.857629206D0
c* Common block data for quadrature weights and points
      data x/0.960289856497536d0, 0.796666477413627d0,
     1       0.525532409916329d0, 0.183434642495650d0/,
     2     w/0.101228536290376d0, 0.222381034453374d0,
     3       0.313706645877887d0, 0.362683783378362d0/
c***********************************************************************
      ERROR= 0.001d0
      IF(VACCNT(ISTATE).EQ.0) THEN
          VACCNT(ISTATE)= VACCNT(ISTATE) + 1
          WRITE(34,100) 
          ENDIF
      ISTART = POTPARI(ISTATE)- 1
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
c.. takes the current temperature from the virial data and changes it to E cm-1
      BTemp = k_cm*Temp(IDAT)
c.. inverts BTemp
      BTempInv = 1.d0/BTemp
c.. sets the constants for the integration of each correction term
      Const(1)= 4.d0*pi*0.602214179d0
      Const(2)= Cu*pi*0.602214129d0*BTempInv/(3.d0*ZMU)
      Const(3)= (pi*0.602214129d0/36.d0)*(Cu*BTempInv/ZMU)**2
c.. sets the number of parameters used for the integration and initializes
c.  the integrals to zero
      nParams = 3 + NCMM(ISTATE) + Nbeta(ISTATE)
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
c.. further subdivides the interval for gaussian integration formula
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
                  RDIST(i)= Re(ISTATE)*DSQRT((1.0001d0 
     1                           + 0.9999d0*YVAL(i))/(1.d0 - YVAL(i)))
                  VDIST(i)= 0.d0
                  ENDDO
             IPASS= 0 
c.. calls VGENP to find the necessary derivatives
              CALL VGENP(ISTATE,RDIST,VDIST,dVdR,d2VdR,IDAT)
c.. the next loop is over each Gaussian point within each subinterval
              DO i= 1,8
c.. some of the terms that will be used in later calculations of the
c.  virial coefficients are constructed here
                  IF(RDIST(I).LT.0.8d0*RE(ISTATE)) THEN
c** Check for and correct for potential turnover problems
                      IF((VDIST(I).LT.0.d0).OR.(DVDR(I).GT.0.d0)
     1                                    .OR.(D2VDR(I).LT.0.d0)) THEN
                          VDIST(I)= 1.d6
                          dVdR(I)= -1.d6
                          d2VdR(I)= 1.d6
                          ENDIF
                      ENDIF
                  EXP_TERM= DEXP(-VDIST(i)*BTempINV)
                  RINV= 1.d0/RDIST(i)
                  int_fact= (Re(ISTATE)/(1.d0 - YVAL(i)))**2
     1                                          *RINV*(b - a)*0.5d0
                  Vsq= VDIST(i)**2
                  VPsq= dVdR(i)**2
                  Rsq= RDIST(i)**2
                  Risq= 1.d0/Rsq
                  BV= BTempInv*VDIST(I)
                  BVP= BTempInv*dVdR(I)
                  BVPP= BTempInv*d2VdR(I)
                  BVsq= BV**2
                  BVPsq= BVP**2
                  BVPcu= BVP**3
                  BVPTtF= BVP**4
                  BVPPsq= BVPP**2
                  XTEMP= ((-6.d0/5.d0)*BVPPsq - (12.d0/5.d0)*Risq*BVPsq
     1    - (20.d0/9.d0)*RINV*BVPcu + (13.d0/30.d0)*BVPTtF) + BV*((4.d0/
     2  5.d0)*BVPPsq + (8.d0/5.d0)*Risq*BVPsq + (56.d0/45.d0)*RINV*BVPcu
     3 - (1.d0/5.d0)*BVPTtF) + BVsq*((-4.d0/25.d0)*BVPPsq - (8.d0/25.d0)
     4      *Risq*BVPsq - (8.d0/45.d0)*RINV*BVPcu + (1.d0/45.d0)*BVPTtF)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c*     Finally the relevant partial derivatives of the virial 
c*     coefficients w.r.t. potential params are calculated
                  DO J = 1,nParams
                      jj = ISTART + J
c* the derivative of the classical expression
                      dBcdP(J)= EXP_TERM*Rsq*BTempInv*DVtot(jj,i)*(1.d0+
     1(2.d0/5.d0)*BV +(2.d0/15.d0)*BVsq - (2.d0/5.d0) - (4.d0/15.d0)*BV)
     2                                                         *int_fact
c* and of the first quantum correction 
                      XTEMP1= (3.d0/5.d0) - (2.d0/5.d0)*BV + 
     1                                                 (2.d0/15.d0)*BVsq
                      XTEMP1= -BVP*XTEMP1*DVtot(jj,i) + 2.d0*dVpdP(jj,i)
     1        *XTEMP1 + BVP*DVtot(jj,i)*((-2.d0/5.d0) + (4.d0/15.d0)*BV)
                      dBq1dP(J)= EXP_TERM*Rsq*BTempInv*BVP*XTEMP1*
     1                                                          int_fact
c..  As the final step these terms all added together in a weighted sum
                      INTEGRALS(J) = INTEGRALS(J) + Const(1)
     1                  *dBcdP(J)*WG(I) + Const(2)*dBq1dP(J)*WG(I)
                      ENDDO
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c. now the integrands are evaluated at each particular Gaussian point
c. and summed together with the proper weightings
                  Bclass= ( 1.d0 - EXP_TERM*(1.d0 + 2.d0*BV/5.d0 +
     1                                  (2.d0/15.d0)*BVsq))*Rsq*int_fact
                  INTEGRALS(nParams+1)= INTEGRALS(NParams+1) 
     1                                                 + Bclass*WG(i)
                  Bq1= EXP_TERM*BVPsq*Rsq*int_fact*((3.d0/5.d0) - 
     1                               (2.d0*BV/5.d0) + (2.d0/15.d0)*BVsq)
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
          Ccorr= Const(1)*INTEGRALS(nParams+1)
          Qcorr= Const(2)*INTEGRALS(nParams+2)
          Q2corr= Const(3)*Integrals(nParams+3)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c. Check to see if BVIR has converged
          IF(counter.GE.2) THEN
              error=0.000001
              IF(DABS(BVIR-dBVIR).LE.(abs(error*BVIR))) THEN
                  check= -1 
                  ENDIF
              ENDIF 
          ENDDO
          WRITE(34,102) Temp(IDAT),Ccorr,Qcorr,Q2corr,BVIR
  100 FORMAT('Temperature Classical   FirstQ   SecondQ   Total',
     1 /24('==='))
  102 FORMAT(2X,F7.2,2X,F8.3,2X,F8.3,2X,F8.3,2X,F8.3)
      RETURN
      END 
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE VGENP(ISTATE,RDIST,VDIST,dVdR,d2VdR2,IDAT)
c***********************************************************************
c** This subroutine will generate function values and derivatives
c   of Morse/Long-Range potentials as required for semiclassical 
c   calculation (with quantum corrections) of virial coefficients and
c   their analytical derivatives in direct hamiltonian fitting
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++  COPYRIGHT 2009-2016 by  R.J. Le Roy, Aleksander Cholewinski and ++
c                      Philip T. Myatt
c   Dept. of Chemistry, Univ. of Waterloo, Waterloo, Ontario, Canada   +
c    This software may not be sold or any other commercial use made    +
c      of it without the express written permission of the authors.    +
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c                   ----- Version of 17 March 2016 -----
c      (after PTW addition of G-TT and specialized HFD potentials)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** On entry:
c   ISTATE  is the electronic state being considered in this CALL.
c   RDIST:  at the 8 input RDIST(i) distances, calculate potl & derivs
c     * return potential function at those points as VDIST, and the
c       first and second radial derivartives as  dVdR  &  d2VdR
c???  * skip partial derivative calculation if  IDAT.le.0
c     * If RDIST.le.0  calculate partial derivatives at distances
c             given by array RD(i,ISTATE) & return them in array DVtot
c** On entry via common blocks:
c  APSE(s).le.0  to use {p,q}-type exponent polynomial of order Nbeta(s)
c     if APSE(s) > 0 \beta(r) is Pashov spline defined by Nbeta(s) points
c* Nbeta(s) is order of the beta(r) exponent polynomial or # spline points
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
cc    INCLUDE 'arrsizes.h'
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c** 'Block' Data Utility routine named: 'arrsizes.h' that governs 
c    array dimensioning in program  dPotFit
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER NISTPMX,NPARMX,NbetaMX,NBOBMX,HPARMX,NDATAMX,
     1  NVIBMX,NBCMX,NSTATEMX,NPNTMX,NROTMX,NCMMAX
c*  NISTPMX  is the maximum number of isotopomers allowed for fit
      PARAMETER (NISTPMX =  12)
c*  NSTATEMX  is maximum no. of electronic states which can be
c             simultaneously fitted to
      PARAMETER (NSTATEMX = 4)
c*  NPARMX  is the largest number of free parameters allowed for fit
c  Since FS origins may be parameters, this is also max. no, data bands
      PARAMETER (NPARMX  = 8000)
c*  NbetaMX  is the largest number of exponent parameters allowed for fit
      PARAMETER (NbetaMX  = 40)
c*  NBOBMX-1  is the highest-order polynomial expansion allowed for the
c      adiabatic or centrifugal Born-Oppenheimer breakdown functions, or 
c      the Lambda-doubling or 2\Sigma splitting radial strength functions
      PARAMETER (NBOBMX  = 15)
c*  HPARMX  is the largest number of Hamiltonian parameters of all types
c    (potential energy, BOB. etc.) for all states.
c           HPARMX >= NSTATEMX*[5 + (NbetaMX+1) + 5*(NBOBMX+1)]
      PARAMETER (HPARMX= NSTATEMX*(5 + (NbetaMX+1) + 5*(NBOBMX+1)))
cc    PARAMETER (HPARMX = 300)
c*  NDATAMX  is largest No. of individual data which may be considered
      PARAMETER (NDATAMX = 35000)
c*  NVIBMX  is the maximum number of vibrational levels of a single
c           state for which data are to be considered
      PARAMETER (NVIBMX    = 200)
** NBCMX  is the maximum number of band constants per vib level to be
c         allowed when doing band constant fits (PSEL= -1)
      PARAMETER (NBCMX = 8) 
c*  NPNTMX  is the largest number of potential data points that can be
c           stored in a single 1D array
      PARAMETER (NPNTMX = 90000)
c*  NROTMX  is the highest order of rotational constants calculated and
c            used for estimating level energies
      PARAMETER (NROTMX = 7)
c*  NCMMAX is the largest number of Cm terms in the MLR or DELR 
c            long-range potential
      PARAMETER (NCMMAX = 12)
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cc    INCLUDE 'BLKISOT.h'
c=======================================================================
c** Isotope/isotopologue numbers, masses & BOB mass scaling factors
c** Array ZK carries about the band constants for all levels of all ISOT
      INTEGER NISTP,NDUNMX,AN(2),MN(2,NISTPMX)
c** NDUNMX is a dummy parameter reqd. for portability of READATA
      PARAMETER (NDUNMX=0)
      REAL*8  ZMASS(3,NISTPMX),RSQMU(NISTPMX),RSQMUP(0:NDUNMX,NISTPMX),
     1 RMUP(0:9,NISTPMX),ZMUA(NISTPMX,NSTATEMX),ZMUB(NISTPMX,NSTATEMX),
     2 ZMTA(NISTPMX,NSTATEMX),ZMTB(NISTPMX,NSTATEMX),
     3  ZK(0:NVIBMX,0:NROTMX,NISTPMX,NSTATEMX)
c
      COMMON /BLKISOT/ZMASS,RSQMU,RSQMUP,RMUP,ZMUA,ZMUB,ZMTA,ZMTB,ZK,
     1  NISTP,AN,MN
c=======================================================================
cc    INCLUDE 'BLKDATA.h'
c=======================================================================
c** Type statements & common block for data
      REAL*8  FREQ(NDATAMX),UFREQ(NDATAMX),DFREQ(NDATAMX),TEMP(NDATAMX),
     1                                               YUNC(NDATAMX),Fqb
      INTEGER  COUNTOT,NFS1,NFSTOT,NBANDTOT,IB(NDATAMX),JP(NDATAMX),
     1 JPP(NDATAMX),VP(NPARMX),VPP(NPARMX),EFP(NDATAMX),EFPP(NDATAMX),
     2 TVUP(NDATAMX),TVLW(NDATAMX),FSBAND(NPARMX),IFXFS(NPARMX),
     3 NFS(NPARMX),IEP(NPARMX),IEPP(NPARMX),ISTP(NPARMX),
     4 IFIRST(NPARMX),ILAST(NPARMX),NTV(NSTATEMX,NISTPMX),FSsame,
     5 NTRANS(NPARMX),IBB(NISTPMX,NSTATEMX,9,NPARMX),JMIN(NPARMX),
     6 JMAX(NPARMX)
      CHARACTER*2 NAME(2) 
      CHARACTER*3 SLABL(-6:NSTATEMX)
      CHARACTER*30 BANDNAME(NPARMX)
      COMMON /DATABLK/Fqb,FREQ,UFREQ,YUNC,DFREQ,TEMP,COUNTOT,NFS1,
     1 NFSTOT,NBANDTOT,IB,JP,JPP,VP,VPP,EFP,EFPP,TVUP,TVLW,FSBAND,IFXFS,
     2 NFS,IEP,IEPP,ISTP,IFIRST,ILAST,NTV,FSsame,
     3 NTRANS,IBB,JMIN,JMAX,NAME,SLABL,BANDNAME
c=======================================================================
cc    INCLUDE 'BLKPOT.h'
c=======================================================================
c** Effective adiabatic radial potential variables.
      INTEGER BOBCN(NSTATEMX),PSEL(NSTATEMX),MAXMIN(NSTATEMX),
     1 IOMEG(NSTATEMX),Nbeta(NSTATEMX),APSE(NSTATEMX),IFXDE(NSTATEMX),
     2 IFXRE(NSTATEMX),IFXCm(NCMMax,NSTATEMX),
     3 IFXBETA(0:NbetaMX,NSTATEMX),NDATPT(NSTATEMX),NCMM(NSTATEMX),
     4 MMLR(NCMMax,NSTATEMX),nPB(NSTATEMX),nQB(NSTATEMX),pAD(NSTATEMX),
     5 qAD(NSTATEMX),LRad(NSTATEMX),pNA(NSTATEMX),qNA(NSTATEMX),
     6 Pqw(NSTATEMX),IVSR(NSTATEMX),IDSTT(NSTATEMX)
c
      REAL*8 DE(NSTATEMX),RE(NSTATEMX),BETA(0:NbetaMX,NSTATEMX),
     1 yqBETA(NbetaMX,NSTATEMX),BETAFX(NPNTMX,NSTATEMX),RH(NSTATEMX),
     2 RMIN(NSTATEMX),RMAX(NSTATEMX),VLIM(NSTATEMX),EPS(NSTATEMX),
     3 betaINF(NSTATEMX),AGPEF(NSTATEMX),BGPEF(NSTATEMX),
     4 CmVAL(NCMMax,NSTATEMX),CmEFF(NCMMax,NSTATEMX),rhoAB(NSTATEMX),
     5 AA(NSTATEMX),BB(NSTATEMX),RREF(NSTATEMX),ASO(NSTATEMX),
     6 R01(NSTATEMX),Q12(NSTATEMX),RD(NPNTMX,NSTATEMX),
     7 VPOT(NPNTMX,NSTATEMX),dCmA(NCMMax,NSTATEMX),dCmB(NCMMax,NSTATEMX)
c
      COMMON /BLKPOT/DE,RE,BETA,yqBETA,BETAFX,RH,RMIN,RMAX,VLIM,EPS,
     1 betaINF,AGPEF,BGPEF,CmVAL,CmEFF,rhoAB,AA,BB,RREF,ASO,R01,Q12,RD,
     2 VPOT,dCmA,dCmB, BOBCN,PSEL,MAXMIN,IOMEG,Nbeta,APSE,IFXDE,IFXRE,
     3 IFXCm,IFXBETA,NDATPT,NCMM,MMLR,nPB,nQB,pAD,qAD,LRad,pNA,qNA,Pqw,
     4 IVSR,IDSTT
c=======================================================================
cc    INCLUDE 'BLKDVDP.h'
c=======================================================================
c** Partial derivative arrays for fits and uncertainties (fununc)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      REAL*8 DVtot(HPARMX,NPNTMX),DLDDRe(NPNTMX,NSTATEMX),
     1  DUADRe(NPNTMX,NSTATEMX),DUBDRe(NPNTMX,NSTATEMX),
     2  DTADRe(NPNTMX,NSTATEMX),DTBDRe(NPNTMX,NSTATEMX),
     3  DBDB(0:NbetaMX,NPNTMX,NSTATEMX),DBDRe(NPNTMX,NSTATEMX),
     4  dVpdP(HPARMX,NPNTMX)
      COMMON/BLKDVDP/DVtot,DUADRe,DUBDRe,DTADRe,DTBDRe,DLDDRe,DBDB,
     1 DBDRe,dVpdP
c=======================================================================
cc    INCLUDE 'BLKBOB.h'
c=======================================================================
c** Born-Oppenheimer Breakdown & doubling function parameters.
c**                       March 16 2012
c=======================================================================
      INTEGER NUA(NSTATEMX),NUB(NSTATEMX),NTA(NSTATEMX),NTB(NSTATEMX),
     1  IFXUA(0:NBOBMX,NSTATEMX),IFXUB(0:NBOBMX,NSTATEMX),
     2  IFXTA(0:NBOBMX,NSTATEMX),IFXTB(0:NBOBMX,NSTATEMX),
     3  NwCFT(NSTATEMX),IFXwCFT(0:NBOBMX,NSTATEMX),efREF(NSTATEMX)
c
      REAL*8 UA(0:NBOBMX,NSTATEMX),UB(0:NBOBMX,NSTATEMX),
     1  TA(0:NBOBMX,NSTATEMX),TB(0:NBOBMX,NSTATEMX),
     2   wCFT(0:NBOBMX,NSTATEMX)
c
      COMMON /BLKBOB/UA,UB,TA,TB,wCFT,NUA,NUB,NTA,NTB,NwCFT,
     1  IFXUA,IFXUB,IFXTA,IFXTB,IFXwCFT,efREF
c=======================================================================
cc    INCLUDE 'BLKBOBRF.h'
c=======================================================================
c** Born-Oppenheimer breakdown radial functions 
      REAL*8 UAR(NPNTMX,NSTATEMX),UBR(NPNTMX,NSTATEMX),
     1 TAR(NPNTMX,NSTATEMX),TBR(NPNTMX,NSTATEMX),wRAD(NPNTMX,NSTATEMX)
c
      COMMON /BLKBOBRF/UAR,UBR,TAR,TBR,wRAD
c=======================================================================
cc    INCLUDE 'BLKCOUNT.h'
c=======================================================================
c    Block data file  BLKCOUNT.h
c=======================================================================
c** Counters for numbers of potential parameters of different types for 
c   each state
      INTEGER  TOTPOTPAR,POTPARI(NSTATEMX),POTPARF(NSTATEMX),
     1  UAPARI(NSTATEMX),UAPARF(NSTATEMX),UBPARI(NSTATEMX),
     2  UBPARF(NSTATEMX),TAPARI(NSTATEMX),TAPARF(NSTATEMX),
     3  TBPARI(NSTATEMX),TBPARF(NSTATEMX),LDPARI(NSTATEMX),
     4  LDPARF(NSTATEMX),HPARF(NSTATEMX),OSEL(NSTATEMX)
c
      COMMON /BLKCOUNT/TOTPOTPAR,POTPARI,POTPARF,UAPARI,UAPARF,UBPARI,
     1  UBPARF,TAPARI,TAPARF,TBPARI,TBPARF,LDPARI,LDPARF,HPARF,OSEL
c=======================================================================
c-----------------------------------------------------------------------
c** Define local variables ...
      INTEGER I,J,I1,ISTATE,IPV,IPVSTART,ISTART,ISTOP,LAMB2,m,npow,
     1  IDAT, NBAND, IISTP,MMLR1D(NCMMax)
      REAL*8 BTEMP,BINF,RVAL(8),RTEMP,RM2,XTEMP,PBTEMP,PETEMP,RET,
     1 FSW,Xtemp2,Btemp2,BMtemp,BMtemp2,RMF,PBtemp2,C3VAL,C3bar,C6bar,
     2 C6adj,C9adj,YP,YQ,YPA,YPB,YQA,YQB,YPE,YPM,YPMA,YPMB,YPP,YQP,YQPA,
     3 YQPB,REp,Req,RDp,RDq,DYPDRE,DYQDRE,VAL,DVAL,HReP,HReQ,SL,SLB,
     4 AREF,AREFp,AREFq, RE3,RE6,RE8,T0,T0P,T1,ULRe,Scalc,dLULRedCm(9),
     5 dLULRedRe,dLULRedDe,dULRdDe,dULRdCm(9),RD3,RD6,RD8,DVDD,RDIST(8),
     6 VDIST(8),BFCT,JFCT,JFCTLD,RETSig,RETPi,RETp,RETm,A0,A1,A2,T2,
     7 REpADA,REpADB,REqADA,REqADB,D2VAL,dYPdR,A3,X,VATT,dVATT,D2VATT,
     8 dYPEdR,dYQdR,d2YPdR,d2YQdR,d2YPEdR,RINV,dDULRdR,d2DULRdR,dULRdR,
     9 d2ULRdR,DXTEMP,D2XTEMP,dVdR(8),d2VdR2(8),dLULRdR,YPPP,dBdR,d2BdR,
     x DX,T1P,T1PP, dULRdRCm(9),dXdP(HPARMX),dXpdP(HPARMX),dLULRdCm(9),
     y DYPEDRE,dVALdRe,dYBdRe,dBpdRe,DYPpDRE,DYPEpdRE,DYQpDRE,dYBpdRe,
     z xBETA(NbetaMX),rKL(NbetaMX,NbetaMX),BR,r,bohr,rhoINT,f2,f2p,f2pp
c***********************************************************************
c** Common block for partial derivatives of potential at the one distance RDIST
c   and HPP derivatives for uncertainties
      REAL*8 dVdPk(HPARMX),dDe(0:NbetaMX),dDedRe
      COMMON /dVdPkBLK/dVdPk,dDe,dDedRe
c=======================================================================
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
      DATA bohr/0.52917721092d0/     !! 2010 physical constants d:mohr12
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
      IF((PSEL(ISTATE).GE.2).AND.(rhoAB(ISTATE).GT.0.d0)) THEN
c ... save uLR powers in a 1D array for calls to SUBROUTINE dampF
          DO  m= 1, NCMM(ISTATE)
              MMLR1D(m)= MMLR(m,ISTATE)
              ENDDO
          ENDIF
c** Initialize parameter counter for this state ...
      IPVSTART= POTPARI(ISTATE) - 1
c=======================================================================
c First ... for the case of an MLR potential ...
c-----------------------------------------------------------------------
          IF(PSEL(ISTATE).EQ.2) THEN
c** First - define values & derivatives of uLR at Re for MLR potential
          ULRe= 0.d0
          T1= 0.d0
          IF(rhoAB(ISTATE).GT.0.d0) THEN
              CALL dampF(RE(ISTATE),rhoAB(ISTATE),NCMM(ISTATE),NCMMAX,
     1                  MMLR1D,IVSR(ISTATE),IDSTT(ISTATE),Dm,Dmp,Dmpp)
              ENDIF
          DO  m= 1,NCMM(ISTATE)
              dLULRedCm(m)= 1.d0/RE(ISTATE)**MMLR(m,ISTATE)
              IF(rhoAB(ISTATE).GT.0.d0) dLULRedCm(m)= Dm(m)*dLULRedCm(m)
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
              npow= Nbeta(ISTATE)
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
                  CALL dampF(RVAL(I),rhoAB(ISTATE),NCMM(ISTATE),NCMMAX,
     1                  MMLR1D,IVSR(ISTATE),IDSTT(ISTATE),Dm,Dmp,Dmpp)
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
                      dULRdRCm(m)= -dULRdCm(m)*RINV*DBLE(MMLR(m,ISTATE))
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
c             BETADIST= VAL
              IF(IDAT.LE.0) GO TO 999
c-              ENDDO
              YPP= 2.d0*DE(ISTATE)*(1.0d0-XTEMP)*XTEMP
              IPV= IPVSTART+2
c... derivatives w.r.t R
              DXTEMP= XTEMP*(dLULRdR - dBdR)
              D2XTEMP= XTEMP*(dBdR**2 - d2BdR + (d2ULRdR 
     1                                           - 2*dBdR*dULRdR)/ULR)
              dVdR(I)= 2.d0*DE(ISTATE)*DXTEMP*(XTEMP - 1.d0)
              d2VdR2(I)= 2.d0*DE(ISTATE)*(DXTEMP**2 + D2XTEMP
     1                                                *(XTEMP - 1.d0))
c *** This is just to write the derivatives for testing
c             IF(RDIST.LT.0) WRITE (40,640) (RVAL,VVAL,dVdR,d2VdR2,
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
     1                                           + dULRdCm(m)/ULR)
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
          ENDIF
c-----------Finished calculations for MLR potential
 
C======================================================================
c For the Tang Toennies potential 
c----------------------------------------------------------------------
      IF(PSEL(ISTATE).EQ.6) THEN
          rhoINT= rhoAB(ISTATE)/3.13d0     !! remove btt(IVSR(ISTATE)/2)
          DO I= 1,8
              VATT= 0.d0
              dVATT= 0.d0
              d2VATT= 0.d0
              r= RDIST(I)
              IF(rhoAB(ISTATE).GT.0.d0) THEN
                  CALL dampF(r,rhoINT,NCMM(ISTATE),NCMMAX,
     1                    MMLR1D,IVSR(ISTATE),IDSTT(ISTATE),Dm,Dmp,Dmpp)
                  DO m= 1,NCMM(ISTATE)
                      T0= CMval(m,ISTATE)/r**MMLR1D(m)
                      VATT= VATT + T0*Dm(m)
                      dVATT= dVATT + T0*(Dmp(m) - Dm(m)*MMLR1D(m)/r)
                      d2VATT= d2VATT + T0*(Dmpp(m) - MMLR1D(m)*(2.d0*
     1                               Dmp(m)- Dm(m)*(MMLR1D(m)+1)/r)/r)
                      ENDDO
                ELSE
                  DO m= 1,NCMM(ISTATE)
                      T0= CMval(m,ISTATE)/r**MMLR1D(m)
                      VATT= VATT + T0
                      dVATT= dVATT + T0*MMLR1D(m)/r
                      d2VATT= d2VATT + T0*MMLR1D(m)*(MMLR1D(m)+1)/r**2
                      ENDDO
                  ENDIF
              T0= r*(BETA(1,ISTATE) + r*BETA(2,ISTATE))
     1                         + (BETA(3,ISTATE) + BETA(4,ISTATE)/r)/r
              T1= BETA(1,ISTATE) + 2.d0*r*BETA(2,ISTATE)
     1                 - (BETA(3,ISTATE) + 2.d0*BETA(4,ISTATE)/r)/r**2
              T2= 2*BETA(2,ISTATE) + (2.d0*BETA(3,ISTATE) 
     1                                   + 6.d0*BETA(4,ISTATE)/r)/r**3
              A0= BETA(5,ISTATE) + r*(BETA(6,ISTATE) 
     1                       + r*(BETA(8,ISTATE) + r*BETA(9,ISTATE)))
     2                                              + BETA(7,ISTATE)/r
              A1= BETA(6,ISTATE) + r*(2.d0*BETA(8,ISTATE)  
     1                  + 3.d0*r*BETA(9,ISTATE)) - BETA(7,ISTATE)/r**2
              A2= 2.d0*BETA(8,ISTATE)+ 6.d0*r*BETA(9,ISTATE)
     1                                      + 2.d0*BETA(7,ISTATE)/r**3
              DX= A0*EXP(-T0)
              VDIST(I)= DX - VATT
              dVdr(I)= DX*(A1/A0 - T1) - dVATT
              d2VdR2(I)= DX*((A2- 2.d0*T1*A1)/A0+ T1**2- T2) - d2VATT
              ENDDO
          ENDIF
c=======================================================================
c ....... for the case of an Aziz'ian HFD-ABC potential ...
c-----------------------------------------------------------------------
      IF((PSEL(ISTATE).EQ.7).AND.(Nbeta(ISTATE).EQ.5)) THEN
          A1= BETA(1,ISTATE)
          A2= BETA(2,ISTATE)
          A3= BETA(3,ISTATE)
          DO  I= 1,8
              r= RDIST(I)
              X= RDIST(I)/RE(ISTATE)
              VATT= 0.d0
              dVATT= 0.d0
              d2VATT= 0.d0
              T1= 1.d0
              T1P= 0.d0
              T1PP= 0.d0
              IF(r.LT.A2) THEN
                  T1= EXP(-A1*(A2/r - 1.d0)**A3)
                  T1P= (A1*A2*A3/(r**2))*((A2/r)-1.d0)**(A3-1.d0)
                  T1PP= T1P*T1P - (A1*A2*A3/r**3)*(A2*(A3-1.d0)/r)
     1       *((A2/r)-1.d0)**(A3-2.d0) - 2.d0*((A2/r)-1.d0)**(A3-1.d0)
                  T1P= T1*T1P
                  T1PP= T1*T1PP
                  ENDIF
              DO M= 1,NCMM(ISTATE)
                   T0= (CmVAL(m,ISTATE)/r**MMLR(m,ISTATE))
                   VATT= VATT+ T0*T1
                   dVATT= dVATT+ (T1P-MMLR(m,ISTATE)*T1/r)*T0
                   d2VATT= d2VATT + (T1PP-2.d0*(T1P*MMLR(m,ISTATE)/r) + 
     1                T1*((MMLR(m,ISTATE)**2)+MMLR(m,ISTATE))/(r**2))*T0
                   ENDDO 
              DX= AA(ISTATE)*(X**BETA(5,ISTATE))*EXP(-r*(BB(ISTATE)
     1                                            + r*BETA(4,ISTATE)))
              VDIST(I)= DX - VATT
              T0= BETA(5,ISTATE)/r - BB(ISTATE)- 2.d0*r*BETA(4,ISTATE)
              dVdR(I)= DX*T0 - DVATT
              d2VdR2(I)= DX*(T0**2 - BETA(5,ISTATE)/r**2 
     1                                 - 2.d0*BETA(4,ISTATE)) - d2VATT
              ENDDO
          ENDIF
c=======================================================================
c... Finally ...For the case of an Aziz'ian HFD-ID potential ...
C-----------------------------------------------------------------------
      IF((PSEL(ISTATE).EQ.7).AND.(Nbeta(ISTATE).EQ.2)) THEN
          A1= BETA(1,ISTATE)
          A2= BETA(2,ISTATE)
          DO  I= ISTART,ISTOP
              r= RDIST(I)
              CALL dampF(r,rhoAB(ISTATE),NCMM(ISTATE),NCMMAX,MMLR1D,
     1                         IVSR(ISTATE),IDSTT(ISTATE),DM,DMP,DMPP)
              X= r/RE(ISTATE)
              BR= RHOab(ISTATE)*r
              VATT= 0.d0
              dVATT= 0.d0
              D2VATT= 0.d0
              f2= (BR/bohr)**1.68d0 *EXP(-0.78d0*BR/bohr)
              f2p= 1.68d0/r - 0.78d0*RHOab(ISTATE)/bohr
              f2pp= - f2*(f2p**2 - 1.68d0/(r**2))
              f2p= -f2*f2p
              f2= 1.d0 - f2
              VATT= 0.d0
              dVATT= 0.d0
              d2VATT= 0.d0
              DO m= 1,NCMM(ISTATE)
                  T0= CmVAL(m,ISTATE)/r**MMLR1D(m)
                  VATT= VATT+ T0*DM(m)
                  dVATT= dVATT+ T0*(f2p*DM(m)+ f2*(DMP(m) -
     1                                       DM(m)*(MMLR1D(m)/r)))
                  d2VATT= d2VATT + T0*(f2pp*DM(m)+ f2*DMPP(m) + 
     1     2.d0*f2p*DMP(m) - 2.d0*(f2p*DM(m) + f2*DMP(m)*(MMLR1D(m)/r))
     2                     + f2*DM(m)*MMLR1D(m)*(MMLR1D(m)+ 1.d0)/r**2)
                  ENDDO
              DX= AA(ISTATE)*(X**A2)*EXP(-r*(BB(ISTATE) + r*A1))
              VDIST(I)= DX - f2*VATT
              T0= A2/r - BB(ISTATE) - 2.d0*r*A1
              dVdR(I)= DX*T0 - dVATT
              d2VdR2(I)= DX*(T0**2 - A2/r**2 - 2.d0*A1) - d2VATT
              ENDDO
          ENDIF
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
cc606 FORMAT(9('===')/'!!!! Extrapolate to correct ',A3,' inner-wall inf
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
          BFCT= 16.857629206d0/(ZMASS(3,IISTP)*RDIST(I)**2)
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
      SUBROUTINE DIFFSTATS(NSTATES,NFPAR,ROBUST,MKPRED,NPTOT,NTVSTOT,
     1                      PRINP)
c** Subroutine to summarise dimensionless standard errors on a band-by-
c  band basis, and (if desired) print [obs.-calc.] values to channel-8.
c-----------------------------------------------------------------------
c                 Version of  30 January  2013
c     last change - added NFPAR` to printout
c-----------------------------------------------------------------------
cc    INCLUDE 'arrsizes.h'
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c** 'Block' Data Utility routine named: 'arrsizes.h' that governs 
c    array dimensioning in program  dPotFit
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER NISTPMX,NPARMX,NbetaMX,NBOBMX,HPARMX,NDATAMX,
     1  NVIBMX,NBCMX,NSTATEMX,NPNTMX,NROTMX,NCMMAX
c*  NISTPMX  is the maximum number of isotopomers allowed for fit
      PARAMETER (NISTPMX =  12)
c*  NSTATEMX  is maximum no. of electronic states which can be
c             simultaneously fitted to
      PARAMETER (NSTATEMX = 4)
c*  NPARMX  is the largest number of free parameters allowed for fit
c  Since FS origins may be parameters, this is also max. no, data bands
      PARAMETER (NPARMX  = 8000)
c*  NbetaMX  is the largest number of exponent parameters allowed for fit
      PARAMETER (NbetaMX  = 40)
c*  NBOBMX-1  is the highest-order polynomial expansion allowed for the
c      adiabatic or centrifugal Born-Oppenheimer breakdown functions, or 
c      the Lambda-doubling or 2\Sigma splitting radial strength functions
      PARAMETER (NBOBMX  = 15)
c*  HPARMX  is the largest number of Hamiltonian parameters of all types
c    (potential energy, BOB. etc.) for all states.
c           HPARMX >= NSTATEMX*[5 + (NbetaMX+1) + 5*(NBOBMX+1)]
      PARAMETER (HPARMX= NSTATEMX*(5 + (NbetaMX+1) + 5*(NBOBMX+1)))
cc    PARAMETER (HPARMX = 300)
c*  NDATAMX  is largest No. of individual data which may be considered
      PARAMETER (NDATAMX = 35000)
c*  NVIBMX  is the maximum number of vibrational levels of a single
c           state for which data are to be considered
      PARAMETER (NVIBMX    = 200)
** NBCMX  is the maximum number of band constants per vib level to be
c         allowed when doing band constant fits (PSEL= -1)
      PARAMETER (NBCMX = 8) 
c*  NPNTMX  is the largest number of potential data points that can be
c           stored in a single 1D array
      PARAMETER (NPNTMX = 90000)
c*  NROTMX  is the highest order of rotational constants calculated and
c            used for estimating level energies
      PARAMETER (NROTMX = 7)
c*  NCMMAX is the largest number of Cm terms in the MLR or DELR 
c            long-range potential
      PARAMETER (NCMMAX = 12)
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cc    INCLUDE 'BLKISOT.h'
c=======================================================================
c** Isotope/isotopologue numbers, masses & BOB mass scaling factors
c** Array ZK carries about the band constants for all levels of all ISOT
      INTEGER NISTP,NDUNMX,AN(2),MN(2,NISTPMX)
c** NDUNMX is a dummy parameter reqd. for portability of READATA
      PARAMETER (NDUNMX=0)
      REAL*8  ZMASS(3,NISTPMX),RSQMU(NISTPMX),RSQMUP(0:NDUNMX,NISTPMX),
     1 RMUP(0:9,NISTPMX),ZMUA(NISTPMX,NSTATEMX),ZMUB(NISTPMX,NSTATEMX),
     2 ZMTA(NISTPMX,NSTATEMX),ZMTB(NISTPMX,NSTATEMX),
     3  ZK(0:NVIBMX,0:NROTMX,NISTPMX,NSTATEMX)
c
      COMMON /BLKISOT/ZMASS,RSQMU,RSQMUP,RMUP,ZMUA,ZMUB,ZMTA,ZMTB,ZK,
     1  NISTP,AN,MN
c=======================================================================
cc    INCLUDE 'BLKDATA.h'
c=======================================================================
c** Type statements & common block for data
      REAL*8  FREQ(NDATAMX),UFREQ(NDATAMX),DFREQ(NDATAMX),TEMP(NDATAMX),
     1                                               YUNC(NDATAMX),Fqb
      INTEGER  COUNTOT,NFS1,NFSTOT,NBANDTOT,IB(NDATAMX),JP(NDATAMX),
     1 JPP(NDATAMX),VP(NPARMX),VPP(NPARMX),EFP(NDATAMX),EFPP(NDATAMX),
     2 TVUP(NDATAMX),TVLW(NDATAMX),FSBAND(NPARMX),IFXFS(NPARMX),
     3 NFS(NPARMX),IEP(NPARMX),IEPP(NPARMX),ISTP(NPARMX),
     4 IFIRST(NPARMX),ILAST(NPARMX),NTV(NSTATEMX,NISTPMX),FSsame,
     5 NTRANS(NPARMX),IBB(NISTPMX,NSTATEMX,9,NPARMX),JMIN(NPARMX),
     6 JMAX(NPARMX)
      CHARACTER*2 NAME(2) 
      CHARACTER*3 SLABL(-6:NSTATEMX)
      CHARACTER*30 BANDNAME(NPARMX)
      COMMON /DATABLK/Fqb,FREQ,UFREQ,YUNC,DFREQ,TEMP,COUNTOT,NFS1,
     1 NFSTOT,NBANDTOT,IB,JP,JPP,VP,VPP,EFP,EFPP,TVUP,TVLW,FSBAND,IFXFS,
     2 NFS,IEP,IEPP,ISTP,IFIRST,ILAST,NTV,FSsame,
     3 NTRANS,IBB,JMIN,JMAX,NAME,SLABL,BANDNAME
c=======================================================================
cc    INCLUDE 'BLKTYPE.h'
c=======================================================================
c** Type statements & common blocks for characterizing transitions
      REAL*8  AVEUFREQ(NPARMX),MAXUFREQ(NPARMX)
      INTEGER NTRANSFS(NISTPMX,NSTATEMX),
     1  NTRANSVIS(NISTPMX,NSTATEMX,NSTATEMX),
     1  NBANDEL(NISTPMX,NSTATEMX,NSTATEMX),
     2  NTRANSIR(NISTPMX,NSTATEMX),NTRANSMW(NISTPMX,NSTATEMX),
     3  NBANDFS(NISTPMX,NSTATEMX),NBANDVIS(NISTPMX,NSTATEMX),
     4  NBANDIR(NISTPMX,NSTATEMX),NBANDMW(NISTPMX,NSTATEMX),
     5  NVVPP(NISTPMX,NSTATEMX),NWIDTH(NISTPMX,NSTATEMX),
     6  NEBPAS(NISTPMX,NSTATEMX),NVIRIAL(NISTPMX,NSTATEMX),
     7  NAcVIR(NISTPMX,NSTATEMX),NBANDS(NISTPMX)
c
      COMMON /BLKTYPE/AVEUFREQ,MAXUFREQ,NTRANSFS,NTRANSVIS,NTRANSIR,
     1  NTRANSMW,NBANDFS,NBANDEL,NBANDVIS,NBANDIR,NBANDMW,NVVPP,NWIDTH,
     2  NEBPAS,NVIRIAL,NAcVIR,NBANDS
c=======================================================================
c
      INTEGER I,IBN,ISOT,ISTATE,ISTATEE,J,K,NSTATES,MKPRED,ROBUST,NFPAR,
     1   NPTOT,NTVSTOT,PRINP
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
c** Summarize data discrepancies for one isotopologue at a time.
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
              IBN= IBB(ISOT,ISTATE,4,I)
              IF(MKPRED.LE.0) THEN
                  CALL BNDERR(IFIRST(IBN),ILAST(IBN),ROBUST,AVE,RMSR,
     1                                             SSQTOT,DFREQ,UFREQ)
                  RMSTOT= RMSTOT+ NTRANS(IBN)*RMSR**2
                  AVETOT= AVETOT+ NTRANS(IBN)*AVE
                  IF((PRINP.EQ.2).OR.(PRINP.EQ.-2)) THEN
                  WRITE(6,606) VP(IBN),VPP(IBN),NTRANS(IBN),JMIN(IBN),
     1                  JMAX(IBN),AVEUFREQ(IBN),MAXUFREQ(IBN),AVE,RMSR,
     2                  BANDNAME(IBN)
                  ELSE
                  WRITE(6,606) VP(IBN),VPP(IBN),NTRANS(IBN),JMIN(IBN),
     1                  JMAX(IBN),AVEUFREQ(IBN),MAXUFREQ(IBN),AVE,RMSR 
                  ENDIF
              ENDIF
              WRITE(8,605)
              IF(MKPRED.LE.0) THEN
                  IF((PRINP.EQ.2).OR.(PRINP.EQ.-2)) THEN    
                      WRITE(8,606) VP(IBN),VPP(IBN),NTRANS(IBN),
     1                  JMIN(IBN),JMAX(IBN),AVEUFREQ(IBN),MAXUFREQ(IBN),
     2                  AVE,RMSR,BANDNAME(IBN)
                     ELSE
                      WRITE(8,606) VP(IBN),VPP(IBN),NTRANS(IBN),
     1                  JMIN(IBN),JMAX(IBN),AVEUFREQ(IBN),MAXUFREQ(IBN),
     2                   AVE,RMSR
                     ENDIF
              ENDIF
              IF(MKPRED.GT.0) THEN
                  WRITE(8,606) VP(IBN),VPP(IBN),NTRANS(IBN),JMIN(IBN),
     1                                                        JMAX(IBN) 
                  WRITE(4,640) VP(IBN),VPP(IBN),SLABL(ISTATE),
     1                                SLABL(ISTATE),(MN(K,ISOT),K=1,2)
                  ENDIF
  640 FORMAT(/2I4,2(2x,"'",A3,"'"),2x,2I4)
              CALL PBNDERR(IBN,MKPRED,NEF)
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
              IBN= IBB(ISOT,ISTATE,3,I)
              IF(MKPRED.LE.0) THEN
                  CALL BNDERR(IFIRST(IBN),ILAST(IBN),ROBUST,AVE,RMSR,
     1                                             SSQTOT,DFREQ,UFREQ)
                  RMSTOT= RMSTOT+ NTRANS(IBN)*RMSR**2
                  AVETOT= AVETOT+ NTRANS(IBN)*AVE
                  IF((PRINP.EQ.2).OR.(PRINP.EQ.-2)) THEN
                  WRITE(6,606) VP(IBN),VPP(IBN),NTRANS(IBN),JMIN(IBN),
     1                  JMAX(IBN),AVEUFREQ(IBN),MAXUFREQ(IBN),AVE,RMSR,
     2                  BANDNAME(IBN)
                  ELSE
                  WRITE(6,606) VP(IBN),VPP(IBN),NTRANS(IBN),JMIN(IBN),
     1                  JMAX(IBN),AVEUFREQ(IBN),MAXUFREQ(IBN),AVE,RMSR
                  ENDIF
              ENDIF
              WRITE(8,605)
              IF(MKPRED.LE.0) THEN
                  IF((PRINP.EQ.2).OR.(PRINP.EQ.-2)) THEN 
                      WRITE(8,606) VP(IBN),VPP(IBN),NTRANS(IBN),
     1         JMIN(IBN),JMAX(IBN),AVEUFREQ(IBN),MAXUFREQ(IBN),AVE,RMSR,
     2                  BANDNAME(IBN)
                  ELSE
                       WRITE(8,606) VP(IBN),VPP(IBN),NTRANS(IBN),
     1         JMIN(IBN),JMAX(IBN),AVEUFREQ(IBN),MAXUFREQ(IBN),AVE,RMSR
                  ENDIF
              ENDIF
              IF(MKPRED.GT.0) THEN
                  WRITE(8,606) VP(IBN),VPP(IBN),NTRANS(IBN),JMIN(IBN),
     1                         JMAX(IBN)
                  WRITE(4,640) VP(IBN),VPP(IBN),
     1             SLABL(ISTATE),SLABL(ISTATE),
     2                                              (MN(K,ISOT),K=1,2)
                  ENDIF
              CALL PBNDERR(IBN,MKPRED,NEF)
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
                  IBN= IBB(ISOT,ISTATE,2,I)
                  IF(IEP(IBN).EQ.ISTATEE) THEN
                      IF(MKPRED.LE.0) THEN
                          CALL BNDERR(IFIRST(IBN),ILAST(IBN),ROBUST,AVE,
     1                                        RMSR,SSQTOT,DFREQ,UFREQ)
                          RMSTOT= RMSTOT+ NTRANS(IBN)*RMSR**2
                          AVETOT= AVETOT+ NTRANS(IBN)*AVE
                          IF((PRINP.EQ.2).OR.(PRINP.EQ.-2)) THEN
                          WRITE(6,606) VP(IBN),VPP(IBN),NTRANS(IBN),
     1         JMIN(IBN),JMAX(IBN),AVEUFREQ(IBN),MAXUFREQ(IBN),AVE,RMSR,
     2                            BANDNAME(IBN)
                          ELSE
                          WRITE(6,606) VP(IBN),VPP(IBN),NTRANS(IBN),
     1         JMIN(IBN),JMAX(IBN),AVEUFREQ(IBN),MAXUFREQ(IBN),AVE,RMSR
                          ENDIF
                          ENDIF
                      WRITE(8,605)
                      IF(MKPRED.LE.0) THEN
                          IF((PRINP.EQ.2).OR.(PRINP.EQ.-2)) THEN 
                              WRITE(8,606) 
     1                 VP(IBN),VPP(IBN),NTRANS(IBN),JMIN(IBN),JMAX(IBN),
     2                AVEUFREQ(IBN),MAXUFREQ(IBN),AVE,RMSR,BANDNAME(IBN)
                          ELSE
                                  WRITE(8,606)
     1                 VP(IBN),VPP(IBN),NTRANS(IBN),JMIN(IBN),JMAX(IBN),
     2                AVEUFREQ(IBN),MAXUFREQ(IBN),AVE,RMSR
                          ENDIF
                      ENDIF
                      IF(MKPRED.GT.0) THEN
                           WRITE(8,606) VP(IBN),VPP(IBN),
     1                  NTRANS(IBN),JMIN(IBN),JMAX(IBN) 
                           WRITE(4,640) VP(IBN),VPP(IBN),SLABL(ISTATE),
     1                                  SLABL(ISTATE),(MN(K,ISOT),K=1,2)
                  ENDIF
                      CALL PBNDERR(IBN,MKPRED,NEF)
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
              IBN= IBB(ISOT,ISTATE,1,I)
              CALL BNDERR(IFIRST(IBN),ILAST(IBN),ROBUST,AVE,RMSR,
     1                                             SSQTOT,DFREQ,UFREQ)
              RMSTOT= RMSTOT+ NTRANS(IBN)*RMSR**2
              AVETOT= AVETOT+ NTRANS(IBN)*AVE
              IF((PRINP.EQ.2).OR.(PRINP.EQ.-2)) THEN
                   WRITE(6,614) VP(IBN),VPP(IBN),NEF(EFP(IFIRST(IBN))),
     1                      NTRANS(IBN),JMIN(IBN),JMAX(IBN),
     2                      AVEUFREQ(IBN),MAXUFREQ(IBN),AVE,RMSR,
     3                      BANDNAME(IBN)
                  WRITE(8,617)
                  WRITE(8,614) VP(IBN),VPP(IBN),NEF(EFP(IFIRST(IBN))),
     1                      NTRANS(IBN),JMIN(IBN),JMAX(IBN),
     2                      AVEUFREQ(IBN),MAXUFREQ(IBN),AVE,RMSR,
     3                      BANDNAME(IBN)
              ELSE
                  WRITE(6,614) VP(IBN),VPP(IBN),NEF(EFP(IFIRST(IBN))),
     1                      NTRANS(IBN),JMIN(IBN),JMAX(IBN),
     2                      AVEUFREQ(IBN),MAXUFREQ(IBN),AVE,RMSR
                  WRITE(8,617)
                  WRITE(8,614) VP(IBN),VPP(IBN),NEF(EFP(IFIRST(IBN))),
     1                      NTRANS(IBN),JMIN(IBN),JMAX(IBN),
     2                      AVEUFREQ(IBN),MAXUFREQ(IBN),AVE,RMSR
              ENDIF
              CALL PBNDERR(IBN,MKPRED,NEF)
              ENDDO
              RMSTOT= DSQRT(RMSTOT/NTRANSFS(ISOT,ISTATE))
              AVETOT= AVETOT/NTRANSFS(ISOT,ISTATE)
              WRITE(6,632) NTRANSFS(ISOT,ISTATE),AVETOT,RMSTOT
          ENDIF
c
      IF(NEBPAS(ISOT,ISTATE).GT.0) THEN
c** Book-keeping for  PAS  data
          IBN= IBB(ISOT,ISTATE,7,1)
          CALL BNDERR(IFIRST(IBN),ILAST(IBN),ROBUST,AVE,RMSR,SSQTOT,
     1                                                    DFREQ,UFREQ)
          IF((PRINP.EQ.2).OR.(PRINP.EQ.-2)) THEN
          WRITE(6,626) NEBPAS(ISOT,ISTATE),SLABL(ISTATE),(NAME(I),
     1 MN(I,ISOT),I=1,2),NTRANS(IBN),JMIN(IBN),JMAX(IBN),
     2      AVEUFREQ(IBN),MAXUFREQ(IBN),AVE,RMSR,BANDNAME(IBN)
          WRITE(8,626) NEBPAS(ISOT,ISTATE),SLABL(ISTATE),(NAME(I),
     1 MN(I,ISOT),I=1,2),NTRANS(IBN),JMIN(IBN),JMAX(IBN),
     2      AVEUFREQ(IBN),MAXUFREQ(IBN),AVE,RMSR,BANDNAME(IBN)
          ELSE
          WRITE(6,626) NEBPAS(ISOT,ISTATE),SLABL(ISTATE),(NAME(I),
     1 MN(I,ISOT),I=1,2),NTRANS(IBN),JMIN(IBN),JMAX(IBN),
     2      AVEUFREQ(IBN),MAXUFREQ(IBN),AVE,RMSR
          WRITE(8,626) NEBPAS(ISOT,ISTATE),SLABL(ISTATE),(NAME(I),
     1 MN(I,ISOT),I=1,2),NTRANS(IBN),JMIN(IBN),JMAX(IBN),
     2 AVEUFREQ(IBN),MAXUFREQ(IBN),AVE,RMSR
          ENDIF
          WRITE(8,627)
          DO  I= IFIRST(IBN),ILAST(IBN)
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
      IF(NVVPP(ISOT,ISTATE).GT.0) THEN
c** Book-keeping for potential function values as data .....
          IBN= IBB(ISOT,ISTATE,5,1)
          CALL BNDERR(IFIRST(IBN),ILAST(IBN),ROBUST,AVE,RMSR,SSQTOT,
     1                                                    DFREQ,UFREQ)
          IF((PRINP.EQ.2).OR.(PRINP.EQ.-2)) THEN
          WRITE(6,638) NVVPP(ISOT,ISTATE),SLABL(ISTATE),
     1              (NAME(I),MN(I,ISOT),I=1,2),NTRANS(IBN),
     2              AVEUFREQ(IBN),MAXUFREQ(IBN),AVE,RMSR,
     3              BANDNAME(IBN)
          WRITE(8,638) NVVPP(ISOT,ISTATE),SLABL(ISTATE),
     1              (NAME(I),MN(I,ISOT),I=1,2),NTRANS(IBN),
     2              AVEUFREQ(IBN),MAXUFREQ(IBN),AVE,RMSR,
     3              BANDNAME(IBN)
          ELSE
          WRITE(6,639) NVVPP(ISOT,ISTATE),SLABL(ISTATE),
     1              (NAME(I),MN(I,ISOT),I=1,2),NTRANS(IBN),
     2              AVEUFREQ(IBN),MAXUFREQ(IBN),AVE,RMSR
          WRITE(8,639) NVVPP(ISOT,ISTATE),SLABL(ISTATE),
     1              (NAME(I),MN(I,ISOT),I=1,2),NTRANS(IBN),
     2              AVEUFREQ(IBN),MAXUFREQ(IBN),AVE,RMSR
          ENDIF
          DO  J= IFIRST(IBN),ILAST(IBN)
              WRITE(6,637) TEMP(J),FREQ(J),UFREQ(J),DFREQ(J),
     1                                               DFREQ(J)/UFREQ(J)
              WRITE(8,637) TEMP(J),FREQ(J),UFREQ(J),DFREQ(J),
     1                                               DFREQ(J)/UFREQ(J)
              ENDDO
          ENDIF
c
      IF(NWIDTH(ISOT,ISTATE).GT.0) THEN
c** Book-keeping for  Tunneling Width  data
          IBN= IBB(ISOT,ISTATE,6,1)    
          CALL BNDERR(IFIRST(IBN),ILAST(IBN),ROBUST,AVE,RMSR,SSQTOT,
     1                                                    DFREQ,UFREQ)
          IF((PRINP.EQ.2).OR.(PRINP.EQ.-2)) THEN
          WRITE(6,620) NWIDTH(ISOT,ISTATE),SLABL(ISTATE),
     1              (NAME(I),MN(I,ISOT),I=1,2),NTRANS(IBN),JMIN(IBN),
     2              JMAX(IBN),AVEUFREQ(IBN),MAXUFREQ(IBN),AVE,RMSR,
     3              BANDNAME(IBN)
          WRITE(8,620) NWIDTH(ISOT,ISTATE),SLABL(ISTATE),
     1              (NAME(I),MN(I,ISOT),I=1,2),NTRANS(IBN),JMIN(IBN),
     2              JMAX(IBN),AVEUFREQ(IBN),MAXUFREQ(IBN),AVE,RMSR,
     3              BANDNAME(IBN)
          ELSE
          WRITE(6,621) NWIDTH(ISOT,ISTATE),SLABL(ISTATE),
     1              (NAME(I),MN(I,ISOT),I=1,2),NTRANS(IBN),JMIN(IBN),
     2              JMAX(IBN),AVEUFREQ(IBN),MAXUFREQ(IBN),AVE,RMSR
          WRITE(8,621) NWIDTH(ISOT,ISTATE),SLABL(ISTATE),
     1              (NAME(I),MN(I,ISOT),I=1,2),NTRANS(IBN),JMIN(IBN),
     2              JMAX(IBN),AVEUFREQ(IBN),MAXUFREQ(IBN),AVE,RMSR
          ENDIF
          DO  J= IFIRST(IBN),ILAST(IBN)
              WRITE(6,622) JP(J),JPP(J),NEF(EFPP(J)),FREQ(J),UFREQ(J),
     1                                      DFREQ(J),DFREQ(J)/UFREQ(J)
              WRITE(8,622) JP(J),JPP(J),NEF(EFPP(J)),FREQ(J),
     1                             UFREQ(J),DFREQ(J),DFREQ(J)/UFREQ(J)
              ENDDO
          ENDIF
c
      IF(NVIRIAL(ISOT,ISTATE).GT.0) THEN
c** Book-keeping for  Virial Coefficient  data
          IBN= IBB(ISOT,ISTATE,8,1)
          CALL BNDERR(IFIRST(IBN),ILAST(IBN),ROBUST,AVE,RMSR,SSQTOT,
     1                                                    DFREQ,UFREQ)
          IF((PRINP.EQ.2).OR.(PRINP.EQ.-2)) THEN
              WRITE(6,634) NVIRIAL(ISOT,ISTATE),SLABL(ISTATE),
     1              (NAME(I),MN(I,ISOT),I=1,2),' Pressure',NTRANS(IBN),
     2              AVEUFREQ(IBN),MAXUFREQ(IBN),AVE,RMSR,
     3              BANDNAME(IBN)
              WRITE(8,634) NVIRIAL(ISOT,ISTATE),SLABL(ISTATE),
     1              (NAME(I),MN(I,ISOT),I=1,2),' Pressure', NTRANS(IBN),
     2              AVEUFREQ(IBN),MAXUFREQ(IBN),AVE,RMSR,
     3              BANDNAME(IBN)
            ELSE
              WRITE(6,635) NVIRIAL(ISOT,ISTATE),SLABL(ISTATE),
     1              (NAME(I),MN(I,ISOT),I=1,2),' Pressure',NTRANS(IBN),
     2              AVEUFREQ(IBN),MAXUFREQ(IBN),AVE,RMSR
              WRITE(8,635) NVIRIAL(ISOT,ISTATE),SLABL(ISTATE),
     1              (NAME(I),MN(I,ISOT),I=1,2),' Pressure',NTRANS(IBN),
     2              AVEUFREQ(IBN),MAXUFREQ(IBN),AVE,RMSR
            ENDIF
          DO  J= IFIRST(IBN),ILAST(IBN)
              WRITE(8,636) TEMP(J),FREQ(J),UFREQ(J),DFREQ(J),
     1                                               DFREQ(J)/UFREQ(J)
              ENDDO
          ENDIF
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
      IF(NAcVIR(ISOT,ISTATE).GT.0) THEN
c** Book-keeping for  Acoustic Virial Coefficient  data
          IBN= IBB(ISOT,ISTATE,9,1)
          CALL BNDERR(IFIRST(IBN),ILAST(IBN),ROBUST,AVE,RMSR,SSQTOT,
     1                                                    DFREQ,UFREQ)
          IF((PRINP.EQ.2).OR.(PRINP.EQ.-2)) THEN
              WRITE(6,634) NAcVIR(ISOT,ISTATE),SLABL(ISTATE),
     1              (NAME(I),MN(I,ISOT),I=1,2),'Accoustic',NTRANS(IBN),
     2              AVEUFREQ(IBN),MAXUFREQ(IBN),AVE,RMSR,
     3              BANDNAME(IBN)
              WRITE(8,634) NAcVIR(ISOT,ISTATE),SLABL(ISTATE),
     1              (NAME(I),MN(I,ISOT),I=1,2),'Accoustic', NTRANS(IBN),
     2              AVEUFREQ(IBN),MAXUFREQ(IBN),AVE,RMSR,
     3              BANDNAME(IBN)
            ELSE
              WRITE(6,635) NAcVIR(ISOT,ISTATE),SLABL(ISTATE),
     1              (NAME(I),MN(I,ISOT),I=1,2),'Accoustic',NTRANS(IBN),
     2              AVEUFREQ(IBN),MAXUFREQ(IBN),AVE,RMSR
              WRITE(8,635) NAcVIR(ISOT,ISTATE),SLABL(ISTATE),
     1              (NAME(I),MN(I,ISOT),I=1,2),'Accoustic',NTRANS(IBN),
     2              AVEUFREQ(IBN),MAXUFREQ(IBN),AVE,RMSR
            ENDIF
          DO  J= IFIRST(IBN),ILAST(IBN)
              WRITE(8,636) TEMP(J),FREQ(J),UFREQ(J),DFREQ(J),
     1                                               DFREQ(J)/UFREQ(J)
              ENDDO
          ENDIF
c** End of loop over the various (lower) electronic states
   90 CONTINUE
c=======================================================================
      IF(ISOT.LT.NISTP) THEN
c** If NISTP > 1, return to print data summaries for other isotopologues
          ISOT= ISOT+1
          GO TO 10
          ENDIF 
      RMSR= DSQRT(SSQTOT/COUNTOT)
      WRITE(6,624) NFPAR,COUNTOT,RMSR
      IF(NTVSTOT.GT.0) THEN
          RMSR= RMSR*SQRT(COUNTOT/DFLOAT(COUNTOT - NTVSTOT))
          WRITE(6,625) NTVSTOT,(NFPAR-NTVSTOT),(COUNTOT - NTVSTOT), RMSR
          ENDIF
      RETURN
  600 FORMAT(/1x,36('**')/'  Write to Channel-8 Predictions From Complet
     1e Set of Input Parameters!'/1x,36('**'))
  601 FORMAT(/1x,25('**')/'  Predictions From Complete Set of Input Para
     1meters!'/1x,25('**'))
  602 FORMAT(/1x,21('===')/'  *** Discrepancies for',I5,' bands/series o
     1f ',A2,'(',I3,')-',A2,'(',I3,') ***'/1x,21('==='))
  604 FORMAT(/1x,21('===')/I5,' State ',A3,1x,A2,'(',I3,')-',A2,'(',I3,
     1 ') MW transitions in',i4,' vib. levels')
  605 FORMAT(1x,16('==='),'== Avge. ========'/"   v' ",
     2  ' v" #data  J"min  J"max  Av.Unc.  Max.Unc.   Err/Unc   DRMSD'/
     1  1x,13('-----'))
  606 FORMAT(2I4,I6,3x,I4,3x,I4,1x,1P2D9.1,0PF11.5,F8.3,2x,A30)
  608 FORMAT(/1x,63('=')/I5,' State ',A3,1x,A2,'(',I3,')-',A2,'(',I3,
     1 ') InfraRed transitions in',I4,' bands')
  610 FORMAT(/1x,35('==')/I6,1x,A2,'(',I3,')-',A2,'(',i3,')  {State ',
     1  A3,'}--{State ',A3,'} Transitions in',i4,' bands')
  612 FORMAT(/1x,75('=')/I5,' Fluorescence transitions into State ',A3,
     1 2x,A2,'(',I3,')-',A2,'(',I3,') in',i5,' series')
  617 FORMAT(1x,52('='),'= Avge. ',15('=')/"   v'  j' p' ",
     2  '#data  v"min  v"max','  AvgeUnc  Max.Unc.   Err/Unc   DRMSD'/
     3  1x,25('---'))
  614 FORMAT(2I4,A3,I6,2I7,1x,1P2D9.1,0PF11.5,F8.3,A31)
  616 FORMAT(/1x,66('=')/1x,I3,' State ',A2,1x,A2,'(',I3,')-',A2,'(',
     1 I3,') potential fx. values treated as independent data'/
     2 1x,20('=='),'  Avge.  ',17('=')/' #data   v"min  v"max  AvgeUnc',
     1 '  Max.Unc.  Err/Unc  DRMSD'/1x,55('-')/I5,2I7,2x,1P2D9.1,0PF9.3,
     4 F8.3/1x,30('==')/'    v  p',8x,'Bv',7x,'u(Bv)',4x,
     5 '[calc-obs]  [calc-obs]/unc',/1x,30('--'))
  618 FORMAT(I5,A3,2x,F12.8,1PD9.1,0PF13.8,F12.4)
  620 FORMAT(/1x,73('=')/1x,I3,' State ',A3,1x,A2,'(',I3,')-',A2,'(',
     1 I3,') Tunneling Widths treated as independent data'/1x,20('=='),
     2 '  Avge.  ',24('=')/' #data   v"min  v"max  AvgeUnc  Max.Unc.  Er
     3r/Unc  DRMSD'/1x,55('-')/I5,2I7,2x,1P2D9.1,0PF9.3,F8.3,A32/
     4 1x,59('=')/'   v   J  p     Width',7x,'u(Width)  [calc-obs]  [cal
     5c-obs]/unc'/1x,59('-'))
  621 FORMAT(/1x,73('=')/1x,I3,' State ',A3,1x,A2,'(',I3,')-',A2,'(',
     1 I3,') Tunneling Widths treated as independent data'/1x,20('=='),
     2 '  Avge.  ',24('=')/' #data   v"min  v"max  AvgeUnc  Max.Unc.  Er
     3r/Unc  DRMSD'/1x,55('-')/I5,2I7,2x,1P2D9.1,0PF9.3,F8.3/
     4 1x,59('=')/'   v   J  p     Width',7x,'u(Width)  [calc-obs]  [cal
     5c-obs]/unc'/1x,59('-'))
  622 FORMAT(2I4,A3,1PD14.6,D10.1,D13.2,0PF10.3)
  624 FORMAT(/1x,39('==')/' Fit of ',I5,' total param to',i6,' data yiel
     1ds   DRMS(devn.)=',G15.8/1x,39('==')) 
  625 FORMAT('  ... & after correcting for the',I5,' term values with on
     1ly onee transition ...'/' Fit of ',I6,' final param to',i7, ' data
     2 yields   DRMS(devn.)=',G15.8/1x,39('=='))
  626 FORMAT(/1x,29('==')/I5,' PAS Binding Energies for State ',A3,2x,
     1 A2,'(',I3,')-',A2,'(',I3,')'/1x,50('='),' Avge. ',('=')/
     2 ' #data  v_min  v_max   AvgeUnc  Max.Unc.  Err/Unc  DRMSD'/
     3  1x,29('--')/I5,2I7,2x,1P2D9.1,0PF9.3,F8.3,A32)
  627 FORMAT(1x,51('='),'  calc-obs'/'   v   j  p       PAS(Eb)        u
     1(Eb)         calc-obs   /u(FREQ)'/1x,30('--'))
  628 FORMAT(2I4,A3,F15.8,F13.8,F13.8,F11.4,1X,A3)
  629 FORMAT(1x,30('=='))
  630 FORMAT(1x,7('--'),' For these',i6,' lines, overall:',F11.5,F8.3)
  632 FORMAT(1x,17('-'),' For these',i6,' lines, overall:',F11.5,F8.3)
  634 FORMAT(/1x,55('=')/I5,'  State ',A3,2x,A2,'(',I3,')-',A2,'(',
     1 I3,') ',A9,' Virial coefficients'/1x,10('===='),' Avge. ',
     2 ('====')/5x,'#data    Av.Unc.    Max.Unc.    Err/Unc   DRMSD'/
     3 1x,10('-----')/I9,1PD12.2,D12.2,0PF11.5,F8.3,A30/1x,23('=='),
     4 ' calc-obs'/ 5x,'temp.',4x,'Bvir(obs)',4x,
     5 'u(Bvir)   calc-obs   /u(Bvir)'/1x,53('-'))
  635 FORMAT(/1x,55('=')/I5,'  State ',A3,2x,A2,'(',I3,')-',A2,'(',
     1 I3,') ',A9,' Virial coefficients'/1x,10('===='),' Avge. ',
     2 ('====')/5x,'#data    Av.Unc.    Max.Unc.    Err/Unc   DRMSD'/
     3 1x,10('-----')/I9,1PD12.2,D12.2,0PF11.5,F8.3/1x,23('=='),
     4 ' calc-obs'/ 5x,'temp.',4x,'Bvir(obs)',4x,
     5 'u(Bvir)   calc-obs   /u(Bvir)'/1x,53('-'))
  636 FORMAT(F11.3,2F10.2,2F11.3)
c 637 FORMAT(F11.6,1P,D14.6,d10.1,D13.5,0P,F8.2)
  637 FORMAT(F11.6,F12.2,F10.2,2F11.3)
  638 FORMAT(/1x,55('=')/I5,'  State ',A3,2x,A2,'(',I3,')-',A2,'(',
     1 I3,')  Potential fx. values'/1x,21('=='),' Avge. ',('====')/5x,
     2 '#data    Av.Unc.    Max.Unc.    Err/Unc   DRMSD'/1x,26('--')
     3 /I9,1PD12.2,D12.2,0PF11.5,F8.3,A30/1x,24('=='),' calc-obs'/
     4  7x,'R',7x,'V(r)',8x,'u(V(r))   calc-obs   /u(V(r))'/
     5  1x,53('-'))
  639 FORMAT(/1x,55('=')/I5,'  State ',A3,2x,A2,'(',I3,')-',A2,'(',
     1 I3,')  Potential fx. values'/1x,21('=='),' Avge. ',('====')/5x,
     2 '#data    Av.Unc.    Max.Unc.    Err/Unc   DRMSD'/1x,26('--')
     3 /I9,1PD12.2,D12.2,0PF11.5,F8.3/1x,24('=='),' calc-obs'/
     4  7x,'R',7x,'V(r)',8x,'u(V(r))   calc-obs   /u(V(r))'/
     5  1x,53('-'))
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE BNDERR(FIRST,LAST,ROBUST,AVEDD,RMSDD,SSQTOT,DFREQ,
     1                                                          UFREQ)
c** Calculate the average (AVEDD) & the root mean square dimensionless 
c  deviation (RSMDD) for the band running from datum # FIRST to LAST.
cc    INCLUDE 'arrsizes.h'
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c** 'Block' Data Utility routine named: 'arrsizes.h' that governs 
c    array dimensioning in program  dPotFit
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER NISTPMX,NPARMX,NbetaMX,NBOBMX,HPARMX,NDATAMX,
     1  NVIBMX,NBCMX,NSTATEMX,NPNTMX,NROTMX,NCMMAX
c*  NISTPMX  is the maximum number of isotopomers allowed for fit
      PARAMETER (NISTPMX =  12)
c*  NSTATEMX  is maximum no. of electronic states which can be
c             simultaneously fitted to
      PARAMETER (NSTATEMX = 4)
c*  NPARMX  is the largest number of free parameters allowed for fit
c  Since FS origins may be parameters, this is also max. no, data bands
      PARAMETER (NPARMX  = 8000)
c*  NbetaMX  is the largest number of exponent parameters allowed for fit
      PARAMETER (NbetaMX  = 40)
c*  NBOBMX-1  is the highest-order polynomial expansion allowed for the
c      adiabatic or centrifugal Born-Oppenheimer breakdown functions, or 
c      the Lambda-doubling or 2\Sigma splitting radial strength functions
      PARAMETER (NBOBMX  = 15)
c*  HPARMX  is the largest number of Hamiltonian parameters of all types
c    (potential energy, BOB. etc.) for all states.
c           HPARMX >= NSTATEMX*[5 + (NbetaMX+1) + 5*(NBOBMX+1)]
      PARAMETER (HPARMX= NSTATEMX*(5 + (NbetaMX+1) + 5*(NBOBMX+1)))
cc    PARAMETER (HPARMX = 300)
c*  NDATAMX  is largest No. of individual data which may be considered
      PARAMETER (NDATAMX = 35000)
c*  NVIBMX  is the maximum number of vibrational levels of a single
c           state for which data are to be considered
      PARAMETER (NVIBMX    = 200)
** NBCMX  is the maximum number of band constants per vib level to be
c         allowed when doing band constant fits (PSEL= -1)
      PARAMETER (NBCMX = 8) 
c*  NPNTMX  is the largest number of potential data points that can be
c           stored in a single 1D array
      PARAMETER (NPNTMX = 90000)
c*  NROTMX  is the highest order of rotational constants calculated and
c            used for estimating level energies
      PARAMETER (NROTMX = 7)
c*  NCMMAX is the largest number of Cm terms in the MLR or DELR 
c            long-range potential
      PARAMETER (NCMMAX = 12)
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
      SUBROUTINE PBNDERR(IBN,MKPRED,NEF)
c** Print to channel-8 a listing of the [obs.-calc.] values for the band
c  running from datum # FIRST to LAST.                           
cc    INCLUDE 'arrsizes.h'             
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c** 'Block' Data Utility routine named: 'arrsizes.h' that governs 
c    array dimensioning in program  dPotFit
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER NISTPMX,NPARMX,NbetaMX,NBOBMX,HPARMX,NDATAMX,
     1  NVIBMX,NBCMX,NSTATEMX,NPNTMX,NROTMX,NCMMAX
c*  NISTPMX  is the maximum number of isotopomers allowed for fit
      PARAMETER (NISTPMX =  12)
c*  NSTATEMX  is maximum no. of electronic states which can be
c             simultaneously fitted to
      PARAMETER (NSTATEMX = 4)
c*  NPARMX  is the largest number of free parameters allowed for fit
c  Since FS origins may be parameters, this is also max. no, data bands
      PARAMETER (NPARMX  = 8000)
c*  NbetaMX  is the largest number of exponent parameters allowed for fit
      PARAMETER (NbetaMX  = 40)
c*  NBOBMX-1  is the highest-order polynomial expansion allowed for the
c      adiabatic or centrifugal Born-Oppenheimer breakdown functions, or 
c      the Lambda-doubling or 2\Sigma splitting radial strength functions
      PARAMETER (NBOBMX  = 15)
c*  HPARMX  is the largest number of Hamiltonian parameters of all types
c    (potential energy, BOB. etc.) for all states.
c           HPARMX >= NSTATEMX*[5 + (NbetaMX+1) + 5*(NBOBMX+1)]
      PARAMETER (HPARMX= NSTATEMX*(5 + (NbetaMX+1) + 5*(NBOBMX+1)))
cc    PARAMETER (HPARMX = 300)
c*  NDATAMX  is largest No. of individual data which may be considered
      PARAMETER (NDATAMX = 35000)
c*  NVIBMX  is the maximum number of vibrational levels of a single
c           state for which data are to be considered
      PARAMETER (NVIBMX    = 200)
** NBCMX  is the maximum number of band constants per vib level to be
c         allowed when doing band constant fits (PSEL= -1)
      PARAMETER (NBCMX = 8) 
c*  NPNTMX  is the largest number of potential data points that can be
c           stored in a single 1D array
      PARAMETER (NPNTMX = 90000)
c*  NROTMX  is the highest order of rotational constants calculated and
c            used for estimating level energies
      PARAMETER (NROTMX = 7)
c*  NCMMAX is the largest number of Cm terms in the MLR or DELR 
c            long-range potential
      PARAMETER (NCMMAX = 12)
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cc    INCLUDE 'BLKDATA.h'
c=======================================================================
c** Type statements & common block for data
      REAL*8  FREQ(NDATAMX),UFREQ(NDATAMX),DFREQ(NDATAMX),TEMP(NDATAMX),
     1                                               YUNC(NDATAMX),Fqb
      INTEGER  COUNTOT,NFS1,NFSTOT,NBANDTOT,IB(NDATAMX),JP(NDATAMX),
     1 JPP(NDATAMX),VP(NPARMX),VPP(NPARMX),EFP(NDATAMX),EFPP(NDATAMX),
     2 TVUP(NDATAMX),TVLW(NDATAMX),FSBAND(NPARMX),IFXFS(NPARMX),
     3 NFS(NPARMX),IEP(NPARMX),IEPP(NPARMX),ISTP(NPARMX),
     4 IFIRST(NPARMX),ILAST(NPARMX),NTV(NSTATEMX,NISTPMX),FSsame,
     5 NTRANS(NPARMX),IBB(NISTPMX,NSTATEMX,9,NPARMX),JMIN(NPARMX),
     6 JMAX(NPARMX)
      CHARACTER*2 NAME(2) 
      CHARACTER*3 SLABL(-6:NSTATEMX)
      CHARACTER*30 BANDNAME(NPARMX)
      COMMON /DATABLK/Fqb,FREQ,UFREQ,YUNC,DFREQ,TEMP,COUNTOT,NFS1,
     1 NFSTOT,NBANDTOT,IB,JP,JPP,VP,VPP,EFP,EFPP,TVUP,TVLW,FSBAND,IFXFS,
     2 NFS,IEP,IEPP,ISTP,IFIRST,ILAST,NTV,FSsame,
     3 NTRANS,IBB,JMIN,JMAX,NAME,SLABL,BANDNAME
c=======================================================================
      REAL*8 DIV
      INTEGER IBN,I,MKPRED
      CHARACTER*3 marker, NEF(-1:1)
c----------------------------------------------------------------------- 
      IF(MKPRED.LE.0) WRITE(8,600)
      IF(MKPRED.GT.0) WRITE(8,601)
      DO  I= IFIRST(IBN),ILAST(IBN)
          IF(MKPRED.LE.0) THEN
              DIV= DABS(DFREQ(I)/UFREQ(I))
              marker='   '
              IF( (DIV.GE.2.d0).AND.(DIV.LT.4.d0) ) marker='*  '
              IF( (DIV.GE.4.d0).AND.(DIV.LT.8.d0) ) marker='** '
              IF( (DIV.GE.8.d0) ) marker='***'
              IF(IEP(IBN).GT.0) WRITE(8,602) VP(IBN),JP(I),NEF(EFP(I)),
     1         VPP(IBN),JPP(I),NEF(EFPP(I)),FREQ(I),UFREQ(I),DFREQ(I),
     2                                      DFREQ(I)/UFREQ(I),marker
              IF(IEP(IBN).EQ.0) WRITE(8,602) VP(IBN),VPP(IBN),
     1                  NEF(EFP(I)),JP(I),JPP(I),NEF(EFPP(I)),FREQ(I),
     2                      UFREQ(I),DFREQ(I),DFREQ(I)/UFREQ(I),marker
            ELSE
              WRITE(8,602) VP(IBN),JP(I),NEF(EFP(I)),VPP(IBN),JPP(I),
     1                                           NEF(EFPP(I)),DFREQ(I)
              WRITE(4,608) JP(I),EFP(I),JPP(I),EFPP(I),DFREQ(I),UFREQ(I)
            ENDIF
          ENDDO
      WRITE(8,604)
      RETURN                       
  600 FORMAT(1x,60('='),'  calc-obs'/ "   v'  J' p'",
     1  '  v"  J" p"    FREQ(obs)     u(FREQ)     calc-obs  /u(FREQ)'/
     2  1x,69('-'))
  601 FORMAT(1x,36('=')/ "   v'  J' p'",'  v"  J" p"   FREQ(calc)'/
     1  1x,36('-'))
  602 FORMAT(2(2I4,A3),f14.6,2f13.7,f10.4:1x,A3)
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
cc    INCLUDE 'arrsizes.h'
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c** 'Block' Data Utility routine named: 'arrsizes.h' that governs 
c    array dimensioning in program  dPotFit
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER NISTPMX,NPARMX,NbetaMX,NBOBMX,HPARMX,NDATAMX,
     1  NVIBMX,NBCMX,NSTATEMX,NPNTMX,NROTMX,NCMMAX
c*  NISTPMX  is the maximum number of isotopomers allowed for fit
      PARAMETER (NISTPMX =  12)
c*  NSTATEMX  is maximum no. of electronic states which can be
c             simultaneously fitted to
      PARAMETER (NSTATEMX = 4)
c*  NPARMX  is the largest number of free parameters allowed for fit
c  Since FS origins may be parameters, this is also max. no, data bands
      PARAMETER (NPARMX  = 8000)
c*  NbetaMX  is the largest number of exponent parameters allowed for fit
      PARAMETER (NbetaMX  = 40)
c*  NBOBMX-1  is the highest-order polynomial expansion allowed for the
c      adiabatic or centrifugal Born-Oppenheimer breakdown functions, or 
c      the Lambda-doubling or 2\Sigma splitting radial strength functions
      PARAMETER (NBOBMX  = 15)
c*  HPARMX  is the largest number of Hamiltonian parameters of all types
c    (potential energy, BOB. etc.) for all states.
c           HPARMX >= NSTATEMX*[5 + (NbetaMX+1) + 5*(NBOBMX+1)]
      PARAMETER (HPARMX= NSTATEMX*(5 + (NbetaMX+1) + 5*(NBOBMX+1)))
cc    PARAMETER (HPARMX = 300)
c*  NDATAMX  is largest No. of individual data which may be considered
      PARAMETER (NDATAMX = 35000)
c*  NVIBMX  is the maximum number of vibrational levels of a single
c           state for which data are to be considered
      PARAMETER (NVIBMX    = 200)
** NBCMX  is the maximum number of band constants per vib level to be
c         allowed when doing band constant fits (PSEL= -1)
      PARAMETER (NBCMX = 8) 
c*  NPNTMX  is the largest number of potential data points that can be
c           stored in a single 1D array
      PARAMETER (NPNTMX = 90000)
c*  NROTMX  is the highest order of rotational constants calculated and
c            used for estimating level energies
      PARAMETER (NROTMX = 7)
c*  NCMMAX is the largest number of Cm terms in the MLR or DELR 
c            long-range potential
      PARAMETER (NCMMAX = 12)
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      INTEGER NTPMX,VMIN,ISTATE,IDAT,K
      PARAMETER (NTPMX= 1600) 
      INTEGER I,J,IAN1,IAN2,IMN1,IMN2,INPTS,ILR,IR2,JWR,LNPT,LWR,
     1  NCN,NLIN,NPP,NROW,NTP,NUSE
      REAL*8  RFACT,EFACT,RH,RMIN,VLIM,VSHIFT,VV(NPP),RR(NPP),RM2(NPP),
     1  XI(NTPMX),YI(NTPMX),RWR(20),VWR(20),VWRB(3),D1V(3),D1VB(3),
     2  D2V(3),CNN,RDIST(8),VDIST(8),DVDR(8),D2VDR2(8)
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
      VMIN= 1.0d9
      IF(NTP.LE.0) THEN 
          WRITE(6,601) !'comments to say PEC generated as analytic ...'
          DO I=1,NPP,8
              DO J=1,8
                  RDIST(J)= RR(J+I-1)
                  ENDDO
              CALL VGENP(ISTATE,RDIST,VDIST,DVDR,D2VdR2,IDAT)
              DO J=1,8
                  VV(J+I-1)= VDIST(J)
                  IF(VV(J+I-1).LT.VMIN) THEN    !! locate potential minimum
                      VMIN= VV(J+I-1)
                      ENDIF    
                  ENDDO
              ENDDO
              IF((I+J-1).LT.NPP) THEN
                  DO K= I+J-1,NPP
                      VV(K)= VLIM
                      ENDDO
                  ENDIF
c+++ Write for testing ++++++++++++++++++++++++++++++++++
cc            REWIND(10)
              WRITE(10,603)  (RR(I),VV(I),I= 1,NPP,20)
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++              
          RETURN
          ENDIF
c**** End of generation of non-standard PEC *****************************
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
c+++ Write for testing ++++++++++++++++++++++++++++++++++
      REWIND(10)
      WRITE(10,603)  (RR(I),VV(I),I= 1,NPP,20)
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++              
      IF(LNPT.GT.0) WRITE(6,624)
      RETURN
  600 FORMAT(' State has energy asymptote:   Y(lim)=',F12.4,'[cm-1]')
  601 FORMAT(/'NTP is set less than or equal to zero, so use some analyt
     1ic function')
  602 FORMAT(/' **** ERROR in dimensioning of arrays required'
     1 ,' by GENINT;   No. input points ',I5,' > NTPMX =',I4)
  603 FORMAT(5x, F12.4,f14.4)
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
c++ COPYRIGHT 1997-2016  by  J.Y. Seto & R.J. Le Roy  (ver. 27/03/2016)+
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
cc    INCLUDE 'arrsizes.h'
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c** 'Block' Data Utility routine named: 'arrsizes.h' that governs 
c    array dimensioning in program  dPotFit
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER NISTPMX,NPARMX,NbetaMX,NBOBMX,HPARMX,NDATAMX,
     1  NVIBMX,NBCMX,NSTATEMX,NPNTMX,NROTMX,NCMMAX
c*  NISTPMX  is the maximum number of isotopomers allowed for fit
      PARAMETER (NISTPMX =  12)
c*  NSTATEMX  is maximum no. of electronic states which can be
c             simultaneously fitted to
      PARAMETER (NSTATEMX = 4)
c*  NPARMX  is the largest number of free parameters allowed for fit
c  Since FS origins may be parameters, this is also max. no, data bands
      PARAMETER (NPARMX  = 8000)
c*  NbetaMX  is the largest number of exponent parameters allowed for fit
      PARAMETER (NbetaMX  = 40)
c*  NBOBMX-1  is the highest-order polynomial expansion allowed for the
c      adiabatic or centrifugal Born-Oppenheimer breakdown functions, or 
c      the Lambda-doubling or 2\Sigma splitting radial strength functions
      PARAMETER (NBOBMX  = 15)
c*  HPARMX  is the largest number of Hamiltonian parameters of all types
c    (potential energy, BOB. etc.) for all states.
c           HPARMX >= NSTATEMX*[5 + (NbetaMX+1) + 5*(NBOBMX+1)]
      PARAMETER (HPARMX= NSTATEMX*(5 + (NbetaMX+1) + 5*(NBOBMX+1)))
cc    PARAMETER (HPARMX = 300)
c*  NDATAMX  is largest No. of individual data which may be considered
      PARAMETER (NDATAMX = 35000)
c*  NVIBMX  is the maximum number of vibrational levels of a single
c           state for which data are to be considered
      PARAMETER (NVIBMX    = 200)
** NBCMX  is the maximum number of band constants per vib level to be
c         allowed when doing band constant fits (PSEL= -1)
      PARAMETER (NBCMX = 8) 
c*  NPNTMX  is the largest number of potential data points that can be
c           stored in a single 1D array
      PARAMETER (NPNTMX = 90000)
c*  NROTMX  is the highest order of rotational constants calculated and
c            used for estimating level energies
      PARAMETER (NROTMX = 7)
c*  NCMMAX is the largest number of Cm terms in the MLR or DELR 
c            long-range potential
      PARAMETER (NCMMAX = 12)
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cc    INCLUDE 'BLKPOT.h'
c=======================================================================
c** Effective adiabatic radial potential variables.
      INTEGER BOBCN(NSTATEMX),PSEL(NSTATEMX),MAXMIN(NSTATEMX),
     1 IOMEG(NSTATEMX),Nbeta(NSTATEMX),APSE(NSTATEMX),IFXDE(NSTATEMX),
     2 IFXRE(NSTATEMX),IFXCm(NCMMax,NSTATEMX),
     3 IFXBETA(0:NbetaMX,NSTATEMX),NDATPT(NSTATEMX),NCMM(NSTATEMX),
     4 MMLR(NCMMax,NSTATEMX),nPB(NSTATEMX),nQB(NSTATEMX),pAD(NSTATEMX),
     5 qAD(NSTATEMX),LRad(NSTATEMX),pNA(NSTATEMX),qNA(NSTATEMX),
     6 Pqw(NSTATEMX),IVSR(NSTATEMX),IDSTT(NSTATEMX)
c
      REAL*8 DE(NSTATEMX),RE(NSTATEMX),BETA(0:NbetaMX,NSTATEMX),
     1 yqBETA(NbetaMX,NSTATEMX),BETAFX(NPNTMX,NSTATEMX),RH(NSTATEMX),
     2 RMIN(NSTATEMX),RMAX(NSTATEMX),VLIM(NSTATEMX),EPS(NSTATEMX),
     3 betaINF(NSTATEMX),AGPEF(NSTATEMX),BGPEF(NSTATEMX),
     4 CmVAL(NCMMax,NSTATEMX),CmEFF(NCMMax,NSTATEMX),rhoAB(NSTATEMX),
     5 AA(NSTATEMX),BB(NSTATEMX),RREF(NSTATEMX),ASO(NSTATEMX),
     6 R01(NSTATEMX),Q12(NSTATEMX),RD(NPNTMX,NSTATEMX),
     7 VPOT(NPNTMX,NSTATEMX),dCmA(NCMMax,NSTATEMX),dCmB(NCMMax,NSTATEMX)
c
      COMMON /BLKPOT/DE,RE,BETA,yqBETA,BETAFX,RH,RMIN,RMAX,VLIM,EPS,
     1 betaINF,AGPEF,BGPEF,CmVAL,CmEFF,rhoAB,AA,BB,RREF,ASO,R01,Q12,RD,
     2 VPOT,dCmA,dCmB, BOBCN,PSEL,MAXMIN,IOMEG,Nbeta,APSE,IFXDE,IFXRE,
     3 IFXCm,IFXBETA,NDATPT,NCMM,MMLR,nPB,nQB,pAD,qAD,LRad,pNA,qNA,Pqw,
     4 IVSR,IDSTT
c=======================================================================
cc    INCLUDE 'BLKPARAM.h'
c=======================================================================
c** Parameters and count-labels for band constant (PSEL=-1) or term
c   value (PSEL=-2) fits
      REAL*8 TVALUE(NPARMX),ZBC(0:NVIBMX,0:NROTMX,NISTPMX,NSTATEMX),
     1 ZQC(0:NVIBMX,0:NROTMX,NISTPMX,NSTATEMX)
c
      INTEGER NSTATES,NTVALL(0:NSTATEMX),NTVI(NSTATEMX),NTVF(NSTATEMX),
     1 VMIN(NSTATEMX,NISTPMX),VMAX(NSTATEMX,NISTPMX),JTRUNC(NSTATEMX),
     2 EFSEL(NSTATEMX),NBC(0:NVIBMX,NISTPMX,NSTATEMX),
     3 NQC(0:NVIBMX,NISTPMX,NSTATEMX),
     4 BCPARI(0:NVIBMX,NISTPMX,NSTATEMX),
     5 BCPARF(0:NVIBMX,NISTPMX,NSTATEMX),
     6 QCPARI(0:NVIBMX,NISTPMX,NSTATEMX),
     7 QCPARF(0:NVIBMX,NISTPMX,NSTATEMX)
      COMMON /BLKPARAM/TVALUE,ZBC,ZQC,NSTATES,NTVALL,NTVI,NTVF,VMIN,
     1      VMAX,JTRUNC,EFSEL,NBC,NQC,BCPARI,BCPARF,QCPARI,QCPARF
c=======================================================================
cc    INCLUDE 'BLKBOB.h'
c=======================================================================
c** Born-Oppenheimer Breakdown & doubling function parameters.
c**                       March 16 2012
c=======================================================================
      INTEGER NUA(NSTATEMX),NUB(NSTATEMX),NTA(NSTATEMX),NTB(NSTATEMX),
     1  IFXUA(0:NBOBMX,NSTATEMX),IFXUB(0:NBOBMX,NSTATEMX),
     2  IFXTA(0:NBOBMX,NSTATEMX),IFXTB(0:NBOBMX,NSTATEMX),
     3  NwCFT(NSTATEMX),IFXwCFT(0:NBOBMX,NSTATEMX),efREF(NSTATEMX)
c
      REAL*8 UA(0:NBOBMX,NSTATEMX),UB(0:NBOBMX,NSTATEMX),
     1  TA(0:NBOBMX,NSTATEMX),TB(0:NBOBMX,NSTATEMX),
     2   wCFT(0:NBOBMX,NSTATEMX)
c
      COMMON /BLKBOB/UA,UB,TA,TB,wCFT,NUA,NUB,NTA,NTB,NwCFT,
     1  IFXUA,IFXUB,IFXTA,IFXTB,IFXwCFT,efREF
c=======================================================================
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
                              DO  J=1,NBC(I,ISOT,ISTATE)
                                  IPV= IPV+1
                                  PV(IPV)= ZBC(I,J-1,ISOT,ISTATE)
                                  ENDDO
                              IF(NQC(I,ISOT,ISTATE).GT.0) THEN  
                                  DO  J=1,NQC(I,ISOT,ISTATE)
                                      IPV= IPV+1
                                      PV(IPV)= ZQC(I,J-1,ISOT,ISTATE)
                                      ENDDO
                                   ENDIF
                              ENDIF  
                          ENDDO
                      DO  I= VMIN(ISTATE,ISOT),VMAX(ISTATE,ISOT)
                          ENDDO
                      ENDDO
                  ENDIF
              IF(PSEL(ISTATE).GT.0) THEN
c*** Manage parameters for potential function mapping ...
                  IF(PSEL(ISTATE).LT.4) THEN
                      IPV= IPV+ 1
                      PV(IPV)= DE(ISTATE)
                      ENDIF
                  IF(PSEL(ISTATE).LE.4) THEN
                      IPV= IPV+ 1
                      PV(IPV)= RE(ISTATE)
                      ENDIF
                  IF((PSEL(ISTATE).EQ.2).OR.(PSEL(ISTATE).EQ.3)) THEN
                      DO  m= 1,NCMM(ISTATE)
                          IPV= IPV+ 1
                          PV(IPV)= CmVAL(m,ISTATE)
                          ENDDO
                      ENDIF
                  J= 0         !! for all PECs except SE-MLR, TT or HDF 
                  IF((APSE(ISTATE).GT.0).OR.(PSEL(ISTATE).GE.6)) J=1
                  DO  I= J,Nbeta(ISTATE)
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
                              IF(NQC(I,ISOT,ISTATE).GT.0) THEN
                                  DO  J= 1,NQC(I,ISOT,ISTATE)
                                      IPV= IPV+1
                                      ZQC(I,J-1,ISOT,ISTATE)= PV(IPV)
                                      ENDDO
                                  ENDIF
                              ENDIF
                          ENDDO
                      ENDDO
                  ENDIF
              IF(PSEL(ISTATE).GT.0) THEN
c*** Manage parameters for potential function mappings ...
                  IF(PSEL(ISTATE).LT.4) THEN
                      IPV= IPV + 1
                      DE(ISTATE)= PV(IPV)
                      ENDIF
                  IF(PSEL(ISTATE).LE.4) THEN
                      IPV= IPV + 1
                      RE(ISTATE) = PV(IPV)
                      ENDIF
                  IF((PSEL(ISTATE).EQ.2).OR.(PSEL(ISTATE).EQ.3)) THEN
                      DO  m= 1,NCMM(ISTATE)
                          IPV= IPV+ 1
                          CmVAL(m,ISTATE)= PV(IPV) 
                          ENDDO
                      ENDIF
                  J=0         !! for all PECs except SE-MLR, TT or HDF 
                  IF((APSE(ISTATE).GT.0).OR.(PSEL(ISTATE).GE.6)) J=1
                  DO I= J, Nbeta(ISTATE)
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
      SUBROUTINE ALF(NDP,RH,NCN,RR,V,SWF,VLIM,MAXMIN,KVMAX,NVIBMX,VMAXX,
     1                               AFLAG,ZMU,EPS,GV,INNODE,INNR,IWR)
c***********************************************************************
c-----------------------------------------------------------------------
c** The subroutine ALF (Automatic vibrational Level Finder) will
c   automatically generate the eigenvalues from the first vibrational
c   level (v=0) to a user specified level (v=KVMAX) or the highest
c   allowed vibrational level of a given smooth single (or double)
c   minimum potential (V). These energies are stored and returned to the
c   calling program in the molecular constants array GV(v=0-KVMAX).
c** For any errors that cannot be resolved within the subroutine, ALF
c   returns AFLAG with a value that defines which error had occured.
c++++++++++   Version last updated  February 18, 2016 ++++++++++++++++++
c+++++++++++++ {removed BMAX, added VMAXX to CALL} +++++++++++++++++++++
c+++++++++++++   COPYRIGHT 2008-16  by  Robert J. Le Roy   +++++++++++++
c   Dept. of Chemistry, Univ. of Waterloo, Waterloo, Ontario, Canada   +
c    This software may not be sold or any other commercial use made    +
c     of it without the express written permission of the author.      +
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++++++ Please inform me of any bugs, by phone at: (519)888-4051 +++++++
c+++++++++ by e-mail to: leroy@uwaterloo.ca , or by Post at: +++++++++++
c+++ Dept. of Chemistry, Univ. Waterloo, Waterloo, Ontario  N2L 3G1 ++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Uses the Schrodinger solver subroutine SCHRQ.
c** On entry:
c    NDP    is the number of datapoints used for the potential.
c    RR(i)  is the array of radial distances (in Angst.), for i= 1, NDP
c    RH     is the radial mesh step size (in Angst).
c    NCN    is the (integer) inverse power defining the linmiting attractive
c           long-range behaviour of the potential.  For a barrier, set NCN=99
c    RR(i)  is the array of distances at which V(i) is defined
c    V(i)   is the scaled input potential (cm-1).
c           The scaling factor BFCT is (2*mu/hbar^2)*RH^2.
c    VLIM   is the potential asymptote (cm-1).
c    MAXMIN the code STOPS if a search finds more than MAXMIN potential minima
c    KVMAX  is v for the highest vibrational level we wish to find.
c    NVIBMX defines dimension of the external Gv array:  GV(0:NVIBMX)
c    AFLAG  is rot.quantum J for the (centrifugally distorted) potential
c    ZMU    is the reduced mass of the diatom (amu).
c    EPS    is the energy convergence criterion (cm-1).
c    INNODE specifies whether wave fx. initiation @ RMIN=RR(1) starts with
c        a node (normal case: INNODE > 0) or zero slope (when INNODE.le.0)
c    IWR    specifies the level of printing inside SCHRQ
c           <> 0 : print error & warning descriptions.
c           >= 1 : also print final eigenvalues & node count.
c           >= 2 : also show end-of-range wave function amplitudes.
c           >= 3 : print also intermediate trial eigenvalues, etc.
c
c** On exit:
c    KVMAX   is vib.quantum number for the highest vibrational level
c            found (may be less than the input value of KVMAX).
c    VMAXX   is MAX{energy at barrier maximim,asymptote}
c    AFLAG   returns calculation outcome to calling program.
c            >=  0 : found all levels to v=KVMAX{input} & AFLAG= J 
c             = -1 : KVMAX larger than number of levels found.
c    GV(v)   contains the vibrational energy levels found for v=0-KVMAX
c    INNR(v) labels each level as belonging to the inner (INNR = 1) or
c            outer (INNR = 0) well.
c
c** Flags: Modify only when debugging.
c    AWO   specifies the level of printing inside ALF
c          < or > 0 : print error & warning descriptions.
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
      IMPLICIT NONE
      INTEGER NDP,KVMAX,KV,KVB,KVBB,AFLAG,NF,NBEG,NEND,NVIBMX,
     1  NBEGG(0:NVIBMX),NENDD(0:NVIBMX),INNR(0:NVIBMX),ICOR,IWR,
     2  IPMIN(10),IPMINN,I,LTRY,AWO,INNODE,INNER,LPRWF,JROT,NCN,NPMIN,
     3  NPMAX,MAXMIN
c
      REAL*8 RMIN,RH,RBAR,RR(NDP),V(NDP),SWF(NDP),VLIM,EO,ZMU,EPS,
     1  BZ,BFCT,GAMA,VMIN,VMAX,VMAXX,PMAX, ESAV, ZPEHO, DGDV2, 
     2  GV(0:NVIBMX),VPMIN(10),RPMIN(10),VPMAX(10),RPMAX(10)
c
      DATA AWO/1/,LPRWF/0/,KVB/-1/,KVBB/-2/
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Check that the array dimensions are adequate.
      RMIN= RR(1)
      IF(KVMAX.GT.NVIBMX) THEN
          WRITE(6,602) KVMAX, NVIBMX
          STOP
          ENDIF
c
c** Initialize the remaining variables and flags.
      NF= 0                          ! NF is label of level being sought
      LTRY= 0
c** Initialize level counters for each well.
      DO  I= 0,KVMAX
          INNR(I)= -2
          ENDDO
c** Store input rotational quantum number.
      JROT= AFLAG
      AFLAG= -1
c
c** Numerical factor  16.857629206 (+/- 0.000,000,013) based on Compton
c  wavelength of proton & proton mass (u) from 2011 physical constants.
      BZ= ZMU/16.857629206d0
      BFCT= BZ*RH*RH
c
cc    IWR=5
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Locate the potential minima.
      NPMIN= 0
      VMIN= 1.d99
      DO  I= 2,NDP-1
          IF((V(I).LT.V(I-1)).AND.(V(I).LT.V(I+1))) THEN
c.... at each minimum located ...
              NPMIN= NPMIN + 1
              IPMIN(NPMIN)= I
              RPMIN(NPMIN)= RR(I)
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
          RPMIN(NPMIN)= RR(1)
          VMIN= RPMIN(NPMIN)
          WRITE(6,606) VPMIN(1),RR(1)
          ENDIF
c
c** Locate any potential maxima past innermost minimum (if they exists).
      NPMAX= 0
      VMAX= -9.d99
      DO  I= IPMIN(1)+1,NDP-1
          IF((V(I).GT.V(I-1)).AND.(V(I).GT.V(I+1))) THEN
              NPMAX= NPMAX + 1
              RPMAX(NPMAX)= RR(I) 
              VPMAX(NPMAX)= V(I)/BFCT
              IF(VPMAX(NPMAX).GT.VMAX) VMAX= VPMAX(NPMAX)
              IF(NPMAX.EQ.9) GOTO 20               !! array bound stop
              ENDIF
          ENDDO
c** Whether or not internal maxima found, add end-of-range as maximum
   20 NPMAX= NPMAX+ 1
      RPMAX(NPMAX)= RR(NDP)
c?? should this limit be set at  VLIM ??  ... naaahhh
      VPMAX(NPMAX)= V(NDP)/BFCT
      IF(VPMAX(NPMAX).GT.VMAX) VMAX= VPMAX(NPMAX)
      VMAXX= VPMAX(NPMAX)    
      IF(VMAXX.LT.VLIM) VMAXX= VLIM
c
c** For multiple minima, print out potential extrema count
      IF(NPMIN.GT.1) THEN
          WRITE(6,614) NPMIN, (VPMIN(I),I= 1,NPMIN)
          WRITE(6,616) (RPMIN(I), I= 1,NPMIN)
          WRITE(6,618) NPMAX, (VPMAX(I),I= 1,NPMAX)
          WRITE(6,616) (RPMAX(I), I= 1,NPMAX)
          IF(NPMIN.GT.MAXMIN) THEN
c** If PEF has more than MAXMIN minima - print warning & stop
              WRITE(6,620)
              STOP
              ENDIF
          ENDIF
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c*** Use harmonic approximation to estimate zero point energy.
      ZPEHO= DSQRT((V(IPMINN+20)-V(IPMINN))/400.d0)/BFCT
      EO= VMIN + ZPEHO
      EO= VMIN + ZPEHO
      IF(EO.GT.VLIM) THEN
          WRITE(6,612) EO,VLIM
          EO= VLIM - 2.d0
          ENDIF
c
c=========== Begin Actual Eigenvalue Calculation Loop Here =============
c** Compute eigenvalues ... etc. up to the KVMAX'th vibrational level.
c** When attempts to find the next eigenvalue fails, then perhaps the
c   next level is located in a second (inner) well. If so, then the
c   subroutine will set INNER = 1, and attempt to find that level.
c
      ICOR= 0
      INNER= 0
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
c   maximum/asymptote, make one last ditch attempt to find the highest 
c   bound level (quasi or otherwise) in the potential.
              IF(LTRY.LT.1) THEN
                  LTRY= 1
                  KV= 999
                  EO= VMAX - 0.0001d0
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
      IF((NPMIN.GT.1).AND.(EO.LT.VPMAX(1))) THEN    
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Begin by asking if the current level is in a double minimum potential
c   and if so, whether it lies below the barrier maximim and if so, 
c   calculate RBAR = <v,J|r|v,J> to see which well it lies in
          RBAR= 0.d0 
          DO I= NBEG,NEND
          RBAR= RBAR+ RR(I)*SWF(I)**2
          ENDDO
          RBAR= RBAR*RH
          INNER= 0
          IF(RBAR.LT.RPMAX(1)) INNER= 1
          IF(IWR.GT.0) write(6,777) RBAR,RPMAX(1),INNER
  777 FORMAT('  Since   RBAR=',F8.3,'   and  RPMAX=',F8.3,'   set INNER
     1=',I2)         
          ENDIF
      IF(KV.EQ.NF) THEN
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** If calculated vibrational level is the desired level, NF, then increase
c   NF by one and call SCECOR to calculate dG/dv and predict next higher level
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          NBEGG(KV)= NBEG
          NENDD(KV)= NEND
          GV(NF)= EO
          INNR(NF)= INNER
  120     NF= NF + 1
          IF(NF.GT.KVMAX) THEN
c** If we have found all desired levels, then RETURN
              IF((AWO.GT.0).AND.(IWR.GT.0)) WRITE(6,626) JROT,KVMAX
              AFLAG= JROT
              RETURN
              ENDIF
c... Check whether the next level had been found earlier in overshoot.
c    If so, count it in and skip on to the next one
          IF(INNR(NF).GE.0) THEN
              EO= GV(NF)
              INNER= INNR(NF)
              KV= NF
              GOTO 120
              ENDIF
          ICOR= 0
c*** NOW, call SCECOR to calculate dG/dv and predict next higher level
c** EO enters as G(KV) & exits as predicted G(NF=KV+1) w. predicted INNER
          CALL SCECOR(KV,NF,JROT,INNER,ICOR,IWR,EO,RH,BFCT,NDP,NCN,V,
     1                                               VMAXX,VLIM,DGDV2)
          IF(ICOR.GE.11) THEN
              KVMAX= KV             !! for case when vD-v < 1 for v=KV
              GOTO 200
              ENDIF
          IF(EO.GT.VPMAX(NPMAX)) THEN
c... if estimated energy above highest barrier, set value slightly below it
              EO=  VPMAX(NPMAX) - 0.10d0*DGDV2
              ICOR= ICOR+10
            ELSE
              IF(DGDV2.LT.0.d0) THEN
c... SCECOR returned negative phase integral, so quit loop & RETURN
                  WRITE(6,628) JROT,EO
                  AFLAG= -1
                  GOTO 200
                  ENDIF
            ENDIF
          LTRY= 0
          GOTO 100
          ENDIF
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IF(KV.NE.NF) THEN
c*** If last level found was not the desired one ...
         IF(INNR(KV).LT.-1) THEN
c... Record vibrational level (if haven't already) for posterity.
              GV(KV)= EO
              INNR(KV)= INNER
              ENDIF
          ICOR= ICOR+1
          IF(ICOR.LE.10) THEN
c... Call subroutine using semiclassical methods to estimate correct energy
              CALL SCECOR(KV,NF,JROT,INNER,ICOR,IWR,EO,RH,BFCT,NDP,NCN,
     1                                             V,VMAXX,VLIM,DGDV2)
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
          KVMAX= KV         !! modified 10/03/15 !! changed back 9/05/15
          IF(AWO.NE.0) WRITE(6,632) KV, GV(KVMAX)
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
  612 FORMAT('  *** Caution ***  H.Osc.ZPE places   E=',F10.2, '   above
     1   VLIM=',F12.2)
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
     1  ')=',1PD15.8 /)
        END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE SCECOR(KV,KVLEV,JROT,INNER,ICOR,IWR,EO,RH,BFCT,NDP,
     1                                    NCN,V,VMAXX,VLIM,DGDV2)
c** Subroutine calculates (approximate!) semiclassical estimate of 
c  dG/dv for level  v= KV  with energy  EO [cm-1]  on potential 
c  {V(i),i=1,NDP} (in 'internal BFCT units' {V[cm-1]*BFCT}), and uses
c  those results to estimate energy of level  KVLEV (usually = KV+1)
c** If the 'clever' semiclassical procedure fails - try a brute force
c  step-by-step search, using alternately INNER & OUTER well starting
c** VMAXX is height of outermost maximum, or VLIM for barrierless case
c** On return, negative DGDV2 signals error!  No phase integrals found
c=======================================================================
c             Version date:  18 February 2016  {removed BMAX}
c***********************************************************************
      INTEGER I,II,I1,I2,I3,I4,IV1,IV2,INNER,ICOR,JROT,KV,KVB,KVLEV,
     1  KVDIF,NDP,NCN,IDIF,BRUTE,IB,IWR,NPMAX
      REAL*8 EO,DE0,RH,BFCT,ARG2,ARG3,EINT,VPH1,VPH2,DGDV1,DGDV2,DGDVM,
     1  DGDV2P,DGDVB,DGDVBP,EBRUTE,DEBRUTE,DE1,DE2,Y1,Y2,Y3,RT,ANS1,dv1,
     2  dv2,ANS2,XDIF,VLIM,VMAXX,PNCN,PWCN,PP1,VDMV,ENEXT,V(NDP)
      SAVE BRUTE,EBRUTE,DEBRUTE,DGDVB
      DATA DGDVB/-1.d0/,KVB/-1/
c
      DGDV2= -1.d0
      EINT= EO*BFCT
      IF(KVLEV.EQ.0) DGDVB= -1.d0
      KVDIF= KVLEV- KV
      IF(ICOR.EQ.1) BRUTE= 0
      PWCN= 2.d0*NCN/DABS(NCN- 2.d0)
      PNCN= ABS(NCN-2)/DFLOAT(NCN+2)
      DGDVBP= DGDVB**PNCN
      PP1= 1.d0/pNCN + 1.d0
      I3= NDP
      IF(EO.GT.VLIM) THEN
c*** For Quasibound levels, first search inward to classically forbidden
          PWCN= 2.d0
          PNCN= 1.d0 
          PP1= 1.d0
          DO  I= NDP,1,-1
              I3= I
              IF(V(I).GT.EINT) GOTO 8
              ENDDO
          ENDIF
c*** Now, search inward for outermost well turning point
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
c... here collect (v+1/2) and dv/dG integrals over outer well ....
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
c*** Next, search outward from RMIN for innermost turning point
      DO  I= 1,NDP 
          I1= I
          IF(V(I).LT.EINT) GOTO 20
          ENDDO
   20 IF(I1.EQ.1) THEN
c... but if RMIN is in the classically allowed region ... STOP here          
          WRITE(6,602) JROT,EO
          STOP
          ENDIF
      IF(I1.GE.I3) THEN
c*** For single-well potential or above barrier of double-well potential
c   use N-D theory estimate based on 'vD-v' from ratio  of Eb to dG/dv
          VDMV= PWCN*(VMAXX-EO)/DGDV2
          ENEXT= VMAXX - (VMAXX-EO)*((VDMV- KVDIF)/VDMV)**PWCN
          IF(IWR.GE.2) THEN
              IF(ABS(EO).GT.1.d0) WRITE(6,600) ICOR,KV,JROT,EO, 
     1                                                VPH2-0.5d0,DGDV2
              IF(ABS(EO).LE.1.d0) WRITE(6,601) ICOR,KV,JROT,EO,
     1                                                VPH2-0.5d0,DGDV2
              WRITE(6,606) VDMV,ENEXT
              ENDIF    
          IF(VDMV.LT.1.d0) THEN
              ICOR= 100
              IF(IWR.GT.0) WRITE(6,604) KV,EO
            ELSE
               EO= ENEXT
            ENDIF   
          DGDVB= DGDV2
          DGDVBP= DGDVB**PNCN
          KVB= KV
          INNER= 0
          RETURN
          ENDIF
c
c*** For a double-well potential, now collect vibrational phase and its 
c   energy derivative over the inner well
      Y1= EINT- V(I1-1)
      Y2= EINT- V(I1)
      Y3= EINT- V(I1+1)
      CALL LEVQAD(Y1,Y2,Y3,RH,RT,ANS1,ANS2)
      ARG2= DSQRT(Y3)
      VPH1= 0.5d0*ARG2 + ANS2/RH
      DGDV1= 0.5d0/ARG2 + ANS1/RH
      DO  I= I1+2,NDP
c... now, collect integral and count nodes outward to second turning point ...
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
c** Not at right level - Check whether looking for higher or lower level ...
      IDIF= SIGN(1,KVDIF)
      XDIF= IDIF
      IF((ICOR.GE.3).AND.((IABS(KVDIF).EQ.1).OR.(BRUTE.GT.0))) GOTO 50
c*** 'Conventional' semiclassical search for nearest INNER or OUTER well level
c... first, determine whether starting level KV was really INNER or OUTER
      dv1= (VPH1-0.5d0) - NINT(VPH1-0.5d0)
      dv2=(VPH2-0.5d0) - NINT(VPH2-0.5d0) 
      IF((DABS(dv2).GT.0.1).AND.(DABS(dv1).LT.0.1)) THEN
          INNER=1
          ENDIF
      IF(INNER.EQ.0) THEN
c... and if current energy EO is for an outer-well level ...
          DE2= DGDV2*XDIF
          IF(IDIF.GT.0) DE1= (Ceiling(VPH1-0.5d0) - (VPH1-0.5d0))*DGDV1
          IF(IDIF.LE.0) DE1= -((VPH1-0.5d0)- Floor(VPH1-0.5d0))*DGDV1
          IF(IWR.GE.2) WRITE(6,610) KV,JROT,EO,VPH1-0.5d0,DGDV1,KVLEV,
     1                                           ICOR,VPH2-0.5d0,DGDV2
        ELSE
c... and if current energy EO is for an inner-well level ...
          DE1= DGDV1*XDIF
          IF(IDIF.GT.0) DE2= (Ceiling(VPH2-0.5d0) - (VPH2-0.5d0))*DGDV2
          IF(IDIF.LE.0) DE2= -(1.d0 - dv2)*DGDV2
          IF(IWR.GE.2) WRITE(6,610) KV,JROT,EO,VPH1-0.5d0,DGDV1,KVLEV,
     1                                           ICOR,VPH2-0.5d0,DGDV2
        ENDIF  
      IF(DABS(DE2).LT.DABS(DE1)) THEN
c... for case in which predict that next level will be OUTER
          INNER= 0
          EO= EO+ DE2
        ELSE
c... for case in which predict that next level will be INNER
          INNER= 1
          EO= EO+ DE1
        ENDIF
      RETURN
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
c... in subsequent EVEN steps, lower EO by DEBRUTE/10 for same INNER
      IF((IB+IB).EQ.BRUTE) THEN
          EBRUTE= EBRUTE+ DEBRUTE
          EO= EBRUTE
          RETURN
        ELSE
c... in subsequent ODD steps, lower repeat previous EO with INNER changed
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
  601 FORMAT(' Single well  ICOR=',I2,':  E(v=',i3,',J=',I3,')=',
     1  1PD12.4,'  v(SC)=',0PF8.3, /63x,'dGdv=',1PD12.4)
  602 FORMAT(/' *** ERROR ***  V(1) < E(J=',i3,')=',f10.2 )
  604 FORMAT(10x,'Find highest level of this potential is   E(v=',i3,
     1                                                   ')=',1PD18.10)
  606 FORMAT(40x,'(vD-v)=',f10.4,'   E(next)=',1PD12.4)
  610 FORMAT(' Double well   E(v=',i3,', J=',I3,')=',f9.3,
     1 ':   v1(SC)=',F7.3,'   dGdv1=',f8.2/8x,'seeking  v=',I3,
     2 ' (ICOR=',I2,')',8x,':   v2(SC)=',F7.3,'   dGdv2=',f8.2 )
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
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
c               COPYRIGHT 1994-2016  by  Robert J. Le Roy              +
c   Dept. of Chemistry, Univ. of Waterloo, Waterloo, Ontario, Canada   +
c    This software may not be sold or any other commercial use made    +
c      of it without the express written permission of the author.     +
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c  Authors: R.J. Le Roy & J. Tellinghuisen         Version of 23/05/2016
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
cc    INCLUDE 'arrsizes.h'             !! bring in array size parameters
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c** 'Block' Data Utility routine named: 'arrsizes.h' that governs 
c    array dimensioning in program  dPotFit
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER NISTPMX,NPARMX,NbetaMX,NBOBMX,HPARMX,NDATAMX,
     1  NVIBMX,NBCMX,NSTATEMX,NPNTMX,NROTMX,NCMMAX
c*  NISTPMX  is the maximum number of isotopomers allowed for fit
      PARAMETER (NISTPMX =  12)
c*  NSTATEMX  is maximum no. of electronic states which can be
c             simultaneously fitted to
      PARAMETER (NSTATEMX = 4)
c*  NPARMX  is the largest number of free parameters allowed for fit
c  Since FS origins may be parameters, this is also max. no, data bands
      PARAMETER (NPARMX  = 8000)
c*  NbetaMX  is the largest number of exponent parameters allowed for fit
      PARAMETER (NbetaMX  = 40)
c*  NBOBMX-1  is the highest-order polynomial expansion allowed for the
c      adiabatic or centrifugal Born-Oppenheimer breakdown functions, or 
c      the Lambda-doubling or 2\Sigma splitting radial strength functions
      PARAMETER (NBOBMX  = 15)
c*  HPARMX  is the largest number of Hamiltonian parameters of all types
c    (potential energy, BOB. etc.) for all states.
c           HPARMX >= NSTATEMX*[5 + (NbetaMX+1) + 5*(NBOBMX+1)]
      PARAMETER (HPARMX= NSTATEMX*(5 + (NbetaMX+1) + 5*(NBOBMX+1)))
cc    PARAMETER (HPARMX = 300)
c*  NDATAMX  is largest No. of individual data which may be considered
      PARAMETER (NDATAMX = 35000)
c*  NVIBMX  is the maximum number of vibrational levels of a single
c           state for which data are to be considered
      PARAMETER (NVIBMX    = 200)
** NBCMX  is the maximum number of band constants per vib level to be
c         allowed when doing band constant fits (PSEL= -1)
      PARAMETER (NBCMX = 8) 
c*  NPNTMX  is the largest number of potential data points that can be
c           stored in a single 1D array
      PARAMETER (NPNTMX = 90000)
c*  NROTMX  is the highest order of rotational constants calculated and
c            used for estimating level energies
      PARAMETER (NROTMX = 7)
c*  NCMMAX is the largest number of Cm terms in the MLR or DELR 
c            long-range potential
      PARAMETER (NCMMAX = 12)
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c** Dimension:  potential arrays  and  vib. level arrays.
c===============================================================      
      INTEGER I,M,IPASS,M1,M2,NBEG,NEND,WARN
      REAL*8 V(NPNTMX),WF0(NPNTMX),RM2(NPNTMX),P(NPNTMX),WF1(NPNTMX),
     1                                       WF2(NPNTMX),RCNST(NROTMX)
      REAL*8 BvWN,DV,DVV,HVV,HV2,LVV,LV2,MVV,MV2,NVV,OVV,EO,E,RH,RHSQ,
     1  ZTW,AR,R2IN,G2,G3,P0,P1,P2,P3,PI,PIF,PRS,PRT,V1,V2,V3,Y1,Y2,Y3,
     2  TSTHv,TSTLv,TSTMv,AMB,AMB1,AMB2,
     3  OV,OV01,OV02,OV03,OV11,OV12,OV13,OV22,OV23,OV33,
     4  PER01,PER02,PER03,PER11,PER12,PER13,PER22,PER23,PER33
c
      IF(NEND.GT.NPNTMX) THEN
          WRITE(6,602) NEND,NPNTMX
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
     1  ' > NPNTMX=',i6)
  603 FORMAT(' ** CAUTION ** Comparison tests for Hv, Lv & Mv give:',
     1 3(1Pd9.1))
  604 FORMAT(' ** CAUTION ** CDJOEL orthogonality tests OV01,OV02 & OV03
     1:',3(1Pd9.1))
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
c***** R.J. Le Roy  subroutine SCHRQ, last modified  9 May 2015 ********
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c                COPYRIGHT 2008-2014  by  Robert J. Le Roy             +
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
c** For energy in (cm-1), BFCT=ZMU(u)*H(Angst)**2/16.857629206 (1/cm-1)
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
c** Start iterative loop; try to converge for up to 30 iterations.
      DO 90 IT= 1,30
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
cc           WRITE(6,603) IT,EO,F,DF,DEPRN,MSAVE,RR,RATIN,RATOUT,
cc   1                                                  XEND,NBEG,ITP1
             WRITE(6,603) IT,EO,DEPRN,MSAVE,RR,RATIN,RATOUT,
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
              IF(IWR.NE.0) WRITE(6,617) JROT,EO,IT,DEP
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
     1            .AND.(RMIN.GT.0.d0)) WRITE(6,614) KVIN,JROT,EO,RATIN
          IF((E.LT.DSOC).AND.(DABS(RATOUT).GT.RATST)) THEN
              WKBTST=0.5d0*DABS(V(NEND)-V(NEND-1))/DSQRT((V(NEND)-E)**3)
              IF(WKBTST.GT.1.d-3) WRITE(6,615) KVIN,JROT,EO,RATOUT,
     1                                                     RATST,WKBTST
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
     1   '  INNER=',i2)
cc602 FORMAT(' ITER    ETRIAL',8X,'F(E)      DF(E)     D(E)',
cc   1 6X,'M    R(M)   /WF(M)   /WF(M)   R(NEND) NBEG ITP1'/
cc   2  1X,99('-'))
  602 FORMAT(' ITER    ETRIAL',7X,'D(E)      M    r(M) wf(1)/wf(M) wf(NE
     1ND)/wf(M) R(NEND) NBEG ITP1'/1X,85('-'))
  603 FORMAT(I4,1PD15.7,D10.2,0P,I7,F7.2,1P2D9.1,0PF8.2,I5,I5)
  604 FORMAT('   NOTE:  for  J=',I3,'   EO=',F12.4,' .ge. V(',i3,')=',
     1  F12.4)
  605 FORMAT(/' Solution of radial Schr. equation for   E(v=',I3,',J=',
     1  I3,') =',F15.7/2x,4('    R(I)   WF(I)   ')/2X,38('--') )
  606 FORMAT(2X,4(F8.3,F11.7))
  607 FORMAT('E(v=',I3,',J=',I3,')=',F11.4,I4,' Iter  R(M)=',F6.2,
     1 '  WF(NBEG=',i6,')/WF(M)=',1PD8.1/36x,'INNER=',I2,6x,
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
  614 FORMAT(' ****** For  v=',I3,', J=',I3,'  E=',G15.8/16x,
     1 'WF(first)/WF(Max)=',D9.2,'  suggests  RMIN  may be too large')
  615 FORMAT(' ****** For  v=',I3,',J=',I3,'  E=',1PD13.6,
     1 '  WF(NEND)/WF(Max)=',D8.1,' >',D8.1/4X,'& initialization ',
     2 'quality test ',1PD8.1,' > 1.D-3   so RMAX may be too small')
  616 FORMAT(' ** WARNING *** For  v=',I2,', J=',I3,' at  E=',G14.7,
     1  ':  inward propagation finds no turning point ... Energy too low
     2 or potential too weak' )
  617 FORMAT(' ** @ J=',I3,'  E=',1PD9.2,' SCHRQ has cgce prob at  IT=',
     1 0P,I3,', so halve  DE=',1PD10.2 )
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
c  [J.Chem.Phys. 119, 7398 (2003); Erratum, ibid, 126, 169904 (2007)]
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
  611 FORMAT(12X,'Log10(lifetime/sec)=',F10.5,' ;   Log10(width/cm-1)=',
     1 F10.5,'   dG/dv=',G12.5,'   V(max)=',G14.7,'(cm-1)')
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
      SUBROUTINE NLLSSRR(NDATA,NPTOT,NPMAX,CYCMAX,IROUND,ROBUST,LPRINT,
     1                      IFXP,YO,YU,YD,PV,PU,PS,CM,TSTPS,TSTPU,DSE)
c**  Program for performing linear or non-linear least-squares fits and
c  (if desired) automatically using sequential rounding and refitting 
c  to minimize the numbers of parameter digits which must be quoted [see
c  R.J. Le Roy, J.Mol.Spectrosc. 191, 223-231 (1998)].         23/03/16
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c             COPYRIGHT 1998-2016  by  Robert J. Le Roy                +
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
c    CYCMAX is the upper bound on the allowed number of iterative cycles
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
c             .NE.0  also print DRMSD and convergence tests on each cycle
c                    and indicate nature of convergence
c              >= 1  also parameters changes & uncertainties, each cycle
c              >= 2  also print parameter change each rounding step
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
      INTEGER I,J,K,L,IDF,ITER,NITER,CYCMAX,IROUND,ISCAL,JROUND,LPRINT,
     1 NDATA,NPTOT,NPMAX,NPARM,NPFIT,JFIX,QUIT,ROBUST,
     2 IFXP(NPMAX),JFXP(NPINTMX)
      REAL*8  YO(NDATA), YU(NDATA), YD(NDATA), PV(NPTOT), PU(NPTOT), 
     1 PS(NPTOT),PSS(NPINTMX),PC(NPINTMX),PX(NPINTMX),
     2 PY(NPINTMX),CM(NPMAX,NPMAX), F95(10),
     3 RMSR, RMSRB, DSE, TSTPS, TSTPSB, TSTPU, TFACT, S, YC, Zthrd
      DATA F95/12.7062D0,4.3027D0,3.1824D0,2.7764D0,2.5706D0,2.4469D0,
     1  2.3646D0,2.3060D0,2.2622D0,2.2281D0/
      IF((NPTOT.GT.NPMAX).OR.(NPTOT.GT.NPINTMX).OR.(NPINTMX.NE.NPMAX)
     1                                      .OR.(NPTOT.GT.NDATA)) THEN
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
      DO 50 ITER= 1, CYCMAX
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
c  values {PV} to generate predicted datum # I [y(calc;I)=YC] and its
c  partial derivatives w.r.t. each of the parameters, returning the 
c  latter in 1-D array PC.  See dummy sample version at end of listing.
c* NOTE 1: if more convenient, DYIDPJ could prepare the y(calc) values 
c     and derivatives for all data at the same time (when I=1), but only
c     returned the values here one datum at a time (for I > 1).]
c* NOTE 2: the partial derivative array PC returned by DYIDPJ must have
c     an entry for every parameter in the model, though for parameters 
c     which are held fixed [JFXP(j)=1], those PC(j) values are ignored.
              CALL DYIDPJ(I,NDATA,NPTOT,YO(I),YC,PV,PC)
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
              YD(I)= YC - YO(I)
              S = 1.D0/YU(I)
cc *** For 'Robust' fitting, adjust uncertainties here
              IF(Zthrd.GT.0.d0) S= 1.d0/DSQRT(YU(I)**2 + Zthrd*YD(I)**2)
              YC= -YD(I)*S
              DSE= DSE+ YC*YC
              IF(NPARM.GT.0) THEN
                  DO  J = 1,NPARM
                      PC(J) = PC(J)*S
                      PS(J) = PS(J)+ PC(J)**2
                      ENDDO
                  CALL QROD(NPARM,NPMAX,NPMAX,CM,PC,PU,YC,PX,PY)
                  ENDIF
              ENDDO
          RMSR= DSQRT(DSE/NDATA)
          IF(NPARM.LE.0) GO TO 60
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
                  YC = 0.D0
                  DO  K = J,NPARM
                      YC = YC + CM(I,K) * CM(J,K)
                      ENDDO
                  CM(I,J) = YC
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
          YC= DSE*0.1d0/DFLOAT(NPARM)
          S= DSE*TFACT
          DO  J = 1,NPARM
              PU(J)= S* PU(J)
              PS(J)= YC*DSQRT(NDATA/PS(J))
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
          IF(LPRINT.NE.0) WRITE(6,604) ITER,RMSR,TSTPS,TSTPU
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
          IF(LPRINT.GE.1) WRITE(6,612) (J,PV(J),PU(J),PS(J),PC(J),
     1                                                      J=1,NPTOT)
          IF(ITER.GT.1) THEN
c** New Convergence test is to require  RMSD to be constant to 1 part in
c   10^7 in adjacent cycles (unlikely to occur by accident)
              IF(ABS((RMSR/RMSRB)-1.d0).LT.1.d-07) THEN
                  IF(LPRINT.NE.0) WRITE(6,607) ITER,
     1                                      ABS(RMSR/RMSRB-1.d0),TSTPS
                  GO TO 54
                  ENDIF
              ENDIF
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
          JFIX= NPTOT
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
      YC= PV(JFIX)
      CALL ROUND(JROUND,NPMAX,NPTOT,NPTOT,JFIX,PV,PU,PS,CM)
      JFXP(JFIX)= 1
      IF(LPRINT.GE.2)
     1       WRITE(6,614) JFIX,YC,PU(JFIX),PS(JFIX),JFIX,PV(JFIX),RMSR
      NPARM= NPARM-1
      IF(NPARM.EQ.0) THEN
c** After rounding complete, make one more pass with all non-fixed 
c  parameters set free to get full correct final correlation matrix, 
c  uncertainties & sensitivities.  Don't update parameters on this pass!
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
  604 FORMAT(' After Cycle #',i2,':  DRMSD=',1PD14.7,'    test(PS)=',
     1  1PD8.1,'   test(PU)=',D8.1)
  606 FORMAT(/' Effective',i3,'-cycle Cgce:  MAX{|change/unc.|}=',
     1  1PD8.1,' < 0.01   DRMSD=',D10.3)
  607 FORMAT(/' Full',i3,'-cycle convergence:  {ABS(RMSR/RMSRB)-1}=',
     1  1PD9.2,'  TSTPS=',D8.1)
  610 FORMAT(/ ' !! CAUTION !! fit of',i5,' parameters to',I6,' data not
     1 converged after',i3,' Cycles'/5x,'DRMS(deviations)=',1PD10.3,
     2 '    test(PS) =',D9.2,'    test(PU) =',D9.2/1x,31('**'))
  612 FORMAT((3x,'PV(',i4,') =',1PD22.14,' (+/-',D8.1,')    PS=',d8.1,
     1  '   PC=',d9.1))
  614 FORMAT(' =',39('==')/' Round Off  PV(',i4,')=',1PD21.13,' (+/-',
     1 D9.2,')    PS=',d9.2/4x,'fix PV(',I4,') as ',D19.11,
     2 '  & refit:  DRMS(deviations)=',D12.5)
  616 FORMAT(/i6,' data fit to',i5,' param. yields  DRMS(devn)=',
     1 1PD14.7:'  tst(PS)=',D8.1)
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
          IF(DABS(Z(1)).LE.0.D0) GOTO 10
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
      CNST= IRND*FCT+ XRND                      !! rounded constant !!
c** Zero parameters more aggressively ... if unc. > 2* value
        if(dabs(PU(IPAR)/PV(IPAR)).GT.2.d0) then
            cnst= 0.d0
            endif
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
c     SUBROUTINE DYIDPJ(I,NDATA,NPTOT,YC,PV,PD)
c** Illustrative dummy version of DYIDPJ for the case of a fit to a
c  power series of order (NPTOT-1) in X(i). ***  For datum number-i, 
c  calculate and return  PD(j)=[partial derivatives of datum-i] w.r.t. 
c  each of the free polynomial coefficients varied in the fit 
c  (for j=1 to NPTOT).  **  Elements of the integer array IFXP indicate
c  whether parameter j is being held fixed [IFXP(j) > 0] or varied in
c  the fit [IFXP(j).le.0].  If the former, the partial derivative 
c  for parameter j should be  PD(j)= 0.0. 
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
      SUBROUTINE FUNUNC(ISTATE,WRITFILE,PU,CM)
c***********************************************************************
c** This subroutine will calculate the uncertainties in the radial
c  strength functions for V(r), phi(r), UA(r), UB(r), tA(r), tB(r) and 
c  wRAD(r) and print them out in `Tecplot' format.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** On entry:
c ... OSEL     print PEF values at every |OSEL|'th mesh point
c           and if OSEL < 0, also the BOB fx an every |OSEL|'th mesh point
c ... ISTATE   electronic state counter
c ... PU(n)    uncertainties in the TOTPOTPAR parameters of the model
c ... CM(n,n)  (symmetric) correlation matrix from the fit
c=======================================================================
cc    INCLUDE 'arrsizes.h'
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c** 'Block' Data Utility routine named: 'arrsizes.h' that governs 
c    array dimensioning in program  dPotFit
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER NISTPMX,NPARMX,NbetaMX,NBOBMX,HPARMX,NDATAMX,
     1  NVIBMX,NBCMX,NSTATEMX,NPNTMX,NROTMX,NCMMAX
c*  NISTPMX  is the maximum number of isotopomers allowed for fit
      PARAMETER (NISTPMX =  12)
c*  NSTATEMX  is maximum no. of electronic states which can be
c             simultaneously fitted to
      PARAMETER (NSTATEMX = 4)
c*  NPARMX  is the largest number of free parameters allowed for fit
c  Since FS origins may be parameters, this is also max. no, data bands
      PARAMETER (NPARMX  = 8000)
c*  NbetaMX  is the largest number of exponent parameters allowed for fit
      PARAMETER (NbetaMX  = 40)
c*  NBOBMX-1  is the highest-order polynomial expansion allowed for the
c      adiabatic or centrifugal Born-Oppenheimer breakdown functions, or 
c      the Lambda-doubling or 2\Sigma splitting radial strength functions
      PARAMETER (NBOBMX  = 15)
c*  HPARMX  is the largest number of Hamiltonian parameters of all types
c    (potential energy, BOB. etc.) for all states.
c           HPARMX >= NSTATEMX*[5 + (NbetaMX+1) + 5*(NBOBMX+1)]
      PARAMETER (HPARMX= NSTATEMX*(5 + (NbetaMX+1) + 5*(NBOBMX+1)))
cc    PARAMETER (HPARMX = 300)
c*  NDATAMX  is largest No. of individual data which may be considered
      PARAMETER (NDATAMX = 35000)
c*  NVIBMX  is the maximum number of vibrational levels of a single
c           state for which data are to be considered
      PARAMETER (NVIBMX    = 200)
** NBCMX  is the maximum number of band constants per vib level to be
c         allowed when doing band constant fits (PSEL= -1)
      PARAMETER (NBCMX = 8) 
c*  NPNTMX  is the largest number of potential data points that can be
c           stored in a single 1D array
      PARAMETER (NPNTMX = 90000)
c*  NROTMX  is the highest order of rotational constants calculated and
c            used for estimating level energies
      PARAMETER (NROTMX = 7)
c*  NCMMAX is the largest number of Cm terms in the MLR or DELR 
c            long-range potential
      PARAMETER (NCMMAX = 12)
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cc    INCLUDE 'BLKPOT.h'
c=======================================================================
c** Effective adiabatic radial potential variables.
      INTEGER BOBCN(NSTATEMX),PSEL(NSTATEMX),MAXMIN(NSTATEMX),
     1 IOMEG(NSTATEMX),Nbeta(NSTATEMX),APSE(NSTATEMX),IFXDE(NSTATEMX),
     2 IFXRE(NSTATEMX),IFXCm(NCMMax,NSTATEMX),
     3 IFXBETA(0:NbetaMX,NSTATEMX),NDATPT(NSTATEMX),NCMM(NSTATEMX),
     4 MMLR(NCMMax,NSTATEMX),nPB(NSTATEMX),nQB(NSTATEMX),pAD(NSTATEMX),
     5 qAD(NSTATEMX),LRad(NSTATEMX),pNA(NSTATEMX),qNA(NSTATEMX),
     6 Pqw(NSTATEMX),IVSR(NSTATEMX),IDSTT(NSTATEMX)
c
      REAL*8 DE(NSTATEMX),RE(NSTATEMX),BETA(0:NbetaMX,NSTATEMX),
     1 yqBETA(NbetaMX,NSTATEMX),BETAFX(NPNTMX,NSTATEMX),RH(NSTATEMX),
     2 RMIN(NSTATEMX),RMAX(NSTATEMX),VLIM(NSTATEMX),EPS(NSTATEMX),
     3 betaINF(NSTATEMX),AGPEF(NSTATEMX),BGPEF(NSTATEMX),
     4 CmVAL(NCMMax,NSTATEMX),CmEFF(NCMMax,NSTATEMX),rhoAB(NSTATEMX),
     5 AA(NSTATEMX),BB(NSTATEMX),RREF(NSTATEMX),ASO(NSTATEMX),
     6 R01(NSTATEMX),Q12(NSTATEMX),RD(NPNTMX,NSTATEMX),
     7 VPOT(NPNTMX,NSTATEMX),dCmA(NCMMax,NSTATEMX),dCmB(NCMMax,NSTATEMX)
c
      COMMON /BLKPOT/DE,RE,BETA,yqBETA,BETAFX,RH,RMIN,RMAX,VLIM,EPS,
     1 betaINF,AGPEF,BGPEF,CmVAL,CmEFF,rhoAB,AA,BB,RREF,ASO,R01,Q12,RD,
     2 VPOT,dCmA,dCmB, BOBCN,PSEL,MAXMIN,IOMEG,Nbeta,APSE,IFXDE,IFXRE,
     3 IFXCm,IFXBETA,NDATPT,NCMM,MMLR,nPB,nQB,pAD,qAD,LRad,pNA,qNA,Pqw,
     4 IVSR,IDSTT
c=======================================================================
cc    INCLUDE 'BLKDVDP.h'
c=======================================================================
c** Partial derivative arrays for fits and uncertainties (fununc)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      REAL*8 DVtot(HPARMX,NPNTMX),DLDDRe(NPNTMX,NSTATEMX),
     1  DUADRe(NPNTMX,NSTATEMX),DUBDRe(NPNTMX,NSTATEMX),
     2  DTADRe(NPNTMX,NSTATEMX),DTBDRe(NPNTMX,NSTATEMX),
     3  DBDB(0:NbetaMX,NPNTMX,NSTATEMX),DBDRe(NPNTMX,NSTATEMX),
     4  dVpdP(HPARMX,NPNTMX)
      COMMON/BLKDVDP/DVtot,DUADRe,DUBDRe,DTADRe,DTBDRe,DLDDRe,DBDB,
     1 DBDRe,dVpdP
c=======================================================================
cc    INCLUDE 'BLKBOB.h'
c=======================================================================
c** Born-Oppenheimer Breakdown & doubling function parameters.
c**                       March 16 2012
c=======================================================================
      INTEGER NUA(NSTATEMX),NUB(NSTATEMX),NTA(NSTATEMX),NTB(NSTATEMX),
     1  IFXUA(0:NBOBMX,NSTATEMX),IFXUB(0:NBOBMX,NSTATEMX),
     2  IFXTA(0:NBOBMX,NSTATEMX),IFXTB(0:NBOBMX,NSTATEMX),
     3  NwCFT(NSTATEMX),IFXwCFT(0:NBOBMX,NSTATEMX),efREF(NSTATEMX)
c
      REAL*8 UA(0:NBOBMX,NSTATEMX),UB(0:NBOBMX,NSTATEMX),
     1  TA(0:NBOBMX,NSTATEMX),TB(0:NBOBMX,NSTATEMX),
     2   wCFT(0:NBOBMX,NSTATEMX)
c
      COMMON /BLKBOB/UA,UB,TA,TB,wCFT,NUA,NUB,NTA,NTB,NwCFT,
     1  IFXUA,IFXUB,IFXTA,IFXTB,IFXwCFT,efREF
c=======================================================================
cc    INCLUDE 'BLKBOBRF.h'
c=======================================================================
c** Born-Oppenheimer breakdown radial functions 
      REAL*8 UAR(NPNTMX,NSTATEMX),UBR(NPNTMX,NSTATEMX),
     1 TAR(NPNTMX,NSTATEMX),TBR(NPNTMX,NSTATEMX),wRAD(NPNTMX,NSTATEMX)
c
      COMMON /BLKBOBRF/UAR,UBR,TAR,TBR,wRAD
c=======================================================================
cc    INCLUDE 'BLKCOUNT.h'
c=======================================================================
c    Block data file  BLKCOUNT.h
c=======================================================================
c** Counters for numbers of potential parameters of different types for 
c   each state
      INTEGER  TOTPOTPAR,POTPARI(NSTATEMX),POTPARF(NSTATEMX),
     1  UAPARI(NSTATEMX),UAPARF(NSTATEMX),UBPARI(NSTATEMX),
     2  UBPARF(NSTATEMX),TAPARI(NSTATEMX),TAPARF(NSTATEMX),
     3  TBPARI(NSTATEMX),TBPARF(NSTATEMX),LDPARI(NSTATEMX),
     4  LDPARF(NSTATEMX),HPARF(NSTATEMX),OSEL(NSTATEMX)
c
      COMMON /BLKCOUNT/TOTPOTPAR,POTPARI,POTPARF,UAPARI,UAPARF,UBPARI,
     1  UBPARF,TAPARI,TAPARF,TBPARI,TBPARF,LDPARI,LDPARF,HPARF,OSEL
c=======================================================================
c-----------------------------------------------------------------------
      INTEGER I,J,JJ,ISTATE,LAM2,NBETAI
      REAL*8 FU,FLAM,RHT,RMAXT,RDVAL,RDVAL2,RDVALLD, PU(NPARMX),
     1  PT(NPARMX),CM(NPARMX,NPARMX) 
      CHARACTER*20 WRITFILE
c
c------------------------------------------------------------------------
c*** Common Block info for  fununc  calculations ***********************  
      REAL*8 Rsr(NPNTMX,NSTATEMX),Vsr(NPNTMX,NSTATEMX),
     1                                            Bsr(NPNTMX,NSTATEMX)
      INTEGER nPointSR(NSTATEMX)
      COMMON /VsrBLK/Rsr,Vsr,Bsr,nPointSR
c
      REAL*8 Rlr(NPNTMX,NSTATEMX),Plr(NPNTMX,NSTATEMX),
     1                                            Blr(NPNTMX,NSTATEMX)
      INTEGER nPointLR(NSTATEMX)
      COMMON /PlrBLK/Rlr,Plr,Blr,nPointLR
c------------------------------------------------------------------------
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      RHT= RD(2,ISTATE)- RD(1,ISTATE)
      RMAXT= RD(NDATPT(ISTATE),ISTATE)
      WRITE(10,900) 'V(r) ', ISTATE, 'V(r) ', WRITFILE
ccc   WRITE(11,900) 'B(r) ', ISTATE, 'B(r) ', WRITFILE
      DO I= 1,nPointSR(ISTATE),MAX(1,IABS(OSEL(ISTATE))/10)
          WRITE(10,909) Rsr(I,ISTATE),Vsr(I,ISTATE)
ccc       WRITE(11,909) Rsr(I,ISTATE),Bsr(I,ISTATE)
          END DO
      IF(OSEL(ISTATE).LT.0) THEN 
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
          ENDIF
      NBETAI= POTPARF(ISTATE) - Nbeta(ISTATE)
      FU= 0.d0
      DO  I= 1,NDATPT(ISTATE),IABS(OSEL(ISTATE))
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
c ... then the exponent coefficient function \betai(i)
ccc       IF(PSEL(ISTATE).LE.5) THEN
ccc           DO  J= NBETAI,POTPARF(ISTATE)
ccc               JJ= J-NBETAI
ccc               PT(J)= PU(J)*DBDB(JJ,I,ISTATE)
ccc               ENDDO
ccc           CALL MMCALC(NBETAI,POTPARF(ISTATE),PT,CM,FU)
ccc           WRITE(11,910) RDVAL,BETAFX(I,ISTATE),FU
ccc           ENDIF
c ... adiabatic BOB correction function for atom-A
          IF(OSEL(1).LT.0) THEN
              IF(NUA(ISTATE).GE.1) THEN
                  DO  J= UAPARI(ISTATE),UAPARF(ISTATE)
                      PT(J)= PU(J)*DVtot(J,I)
                      ENDDO
                  CALL MMCALC(UAPARI(ISTATE),UAPARF(ISTATE),PT,CM,FU)
                  WRITE(12,910) RDVAL,UAR(I,ISTATE),FU
                  ENDIF
c ... adiabatic BOB correction function for atom-B
              IF(NUB(ISTATE).GE.1) THEN
                  DO  J= UBPARI(ISTATE),UBPARF(ISTATE)
                      PT(J)= PU(J)*DVtot(J,I)
                      ENDDO
                  CALL MMCALC(UBPARI(ISTATE),UBPARF(ISTATE),PT,CM,FU)
                  WRITE(13,910) RDVAL,UBR(I,ISTATE),FU
                  ENDIF
c ... centrifugal BOB correction function for atom-A
              IF(NTA(ISTATE).GE.1) THEN
                  DO  J= TAPARI(ISTATE),TAPARF(ISTATE)
                      PT(J)= PU(J)*DVtot(J,I)*RDVAL2
                      ENDDO
                  CALL MMCALC(TAPARI(ISTATE),TAPARF(ISTATE),PT,CM,FU)
                  WRITE(14,910) RDVAL,TAR(I,ISTATE),FU
                  ENDIF
c ... centrifugal BOB correction function for atom-B
              IF(NTB(ISTATE).GE.1) THEN
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
              ENDIF
          END DO
      DO I= 1,nPointLR(ISTATE)
          WRITE(10,909) Rlr(I,ISTATE),Plr(I,ISTATE)
ccc       WRITE(11,909) Rlr(I,ISTATE),Blr(I,ISTATE)
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
cc    INCLUDE 'arrsizes.h'
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c** 'Block' Data Utility routine named: 'arrsizes.h' that governs 
c    array dimensioning in program  dPotFit
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER NISTPMX,NPARMX,NbetaMX,NBOBMX,HPARMX,NDATAMX,
     1  NVIBMX,NBCMX,NSTATEMX,NPNTMX,NROTMX,NCMMAX
c*  NISTPMX  is the maximum number of isotopomers allowed for fit
      PARAMETER (NISTPMX =  12)
c*  NSTATEMX  is maximum no. of electronic states which can be
c             simultaneously fitted to
      PARAMETER (NSTATEMX = 4)
c*  NPARMX  is the largest number of free parameters allowed for fit
c  Since FS origins may be parameters, this is also max. no, data bands
      PARAMETER (NPARMX  = 8000)
c*  NbetaMX  is the largest number of exponent parameters allowed for fit
      PARAMETER (NbetaMX  = 40)
c*  NBOBMX-1  is the highest-order polynomial expansion allowed for the
c      adiabatic or centrifugal Born-Oppenheimer breakdown functions, or 
c      the Lambda-doubling or 2\Sigma splitting radial strength functions
      PARAMETER (NBOBMX  = 15)
c*  HPARMX  is the largest number of Hamiltonian parameters of all types
c    (potential energy, BOB. etc.) for all states.
c           HPARMX >= NSTATEMX*[5 + (NbetaMX+1) + 5*(NBOBMX+1)]
      PARAMETER (HPARMX= NSTATEMX*(5 + (NbetaMX+1) + 5*(NBOBMX+1)))
cc    PARAMETER (HPARMX = 300)
c*  NDATAMX  is largest No. of individual data which may be considered
      PARAMETER (NDATAMX = 35000)
c*  NVIBMX  is the maximum number of vibrational levels of a single
c           state for which data are to be considered
      PARAMETER (NVIBMX    = 200)
** NBCMX  is the maximum number of band constants per vib level to be
c         allowed when doing band constant fits (PSEL= -1)
      PARAMETER (NBCMX = 8) 
c*  NPNTMX  is the largest number of potential data points that can be
c           stored in a single 1D array
      PARAMETER (NPNTMX = 90000)
c*  NROTMX  is the highest order of rotational constants calculated and
c            used for estimating level energies
      PARAMETER (NROTMX = 7)
c*  NCMMAX is the largest number of Cm terms in the MLR or DELR 
c            long-range potential
      PARAMETER (NCMMAX = 12)
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
      INTEGER  IROUND,NPMAX,NPTOT,NPAR1,NPAR2,IPAR,IRND,KRND,LPRINT
      INTEGER  IFXP(NPTOT)
      REAL*8  PV(NPMAX),PU(NPMAX),CNST,CRND,XRND,FCT,YY,UNC
c !!  This only makes sense if ALL param have same magnitude (e.g. Tvj's)
c** Loop over & round off the parameters # NPAR1 to NPAR2 
cc    IF(LPRINT.GE.2) WRITE(6,602)  NPAR2-NPAR1+1,NPTOT,NPAR1,NPAR2
cc    UNC= 99.d99
cc    DO  IPAR= NPAR1, NPAR2      !! search for smallest uncertainty
cc        IF(PU(IPAR).LT.UNC) UNC= PU(IPAR)   !! which is/was used
cc        ENDDO                    !! to round ALL parameters!
      DO  IPAR= NPAR1, NPAR2
c** First ... fiddle with log's to perform the rounding
          XRND= DLOG10(PU(IPAR))
          IRND= INT(XRND)
          IF(XRND.GT.0) IRND=IRND+1
          IRND= IRND- IROUND
          FCT= 10.D0**IRND
          CNST= PV(IPAR)
          YY= CNST
          CRND= PV(IPAR)/FCT
          XRND= 0.d0
c ... if rounding goes past REAL*8 precision, retain unrounded constant
          IF(DABS(CRND).GE.1.D+16) THEN
              WRITE(6,600) IROUND,IPAR
               GO TO 20
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
          PV(IPAR) = CNST
          IFXP(IPAR)= 1
          IF(LPRINT.GE.2) WRITE(6,604) IPAR,YY,PV(IPAR)
  604 FORMAT(5x,'Round parameter #',i4,' from',G20.12,'  to',G20.12)
   20     CONTINUE
          ENDDO
      IPAR= IPAR- 1
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
c** At the position 'x', Scalc is returned as the value of the m'th 
c  of the 'n' Sm(x) function defining a natural cubic spline through the
c  mesh points located at  x= y(x_i), for i=1,n.  LMAX specifies the 
c  maximum number of mesh points x= y(x_i) allowed by the calling program
c---------------------------------------------------------------------
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
c** Call this subroutine with list of the 'n' spline x_i values in array 
c   'x' with maximum dimension 'LMAX' and it will return the LMAX x LMAX
c   array of 'rKL' coefficients used for generating the 'n' S_n(x) 
c   spline coefficient functions
c-----------------------------------------------------------------------
c   
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
          if (aamax.eq.0.) WRITE(6,*) 'singular matrix in ludcmp'
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

c***********************************************************************
      SUBROUTINE UNCBV(NPTOT,PV,PU,CM)
c***********************************************************************
c**  Subroutine to compute the uncertainties in calcuated Bv values, 
c  1) by the finite difference R(0) approximation
c  2) from the LIDE expression
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** On entry:
c ... OSEL     print values at every OSEL'th mesh point
c ... ISTATE   electronic state counter
c ... PU(n)    uncertainties in the TOTPOTPAR parameters of the model
c ... CM(n,n)  (symmetric) correlation matrix from the fit
c=======================================================================
cc    INCLUDE 'arrsizes.h'
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c** 'Block' Data Utility routine named: 'arrsizes.h' that governs 
c    array dimensioning in program  dPotFit
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER NISTPMX,NPARMX,NbetaMX,NBOBMX,HPARMX,NDATAMX,
     1  NVIBMX,NBCMX,NSTATEMX,NPNTMX,NROTMX,NCMMAX
c*  NISTPMX  is the maximum number of isotopomers allowed for fit
      PARAMETER (NISTPMX =  12)
c*  NSTATEMX  is maximum no. of electronic states which can be
c             simultaneously fitted to
      PARAMETER (NSTATEMX = 4)
c*  NPARMX  is the largest number of free parameters allowed for fit
c  Since FS origins may be parameters, this is also max. no, data bands
      PARAMETER (NPARMX  = 8000)
c*  NbetaMX  is the largest number of exponent parameters allowed for fit
      PARAMETER (NbetaMX  = 40)
c*  NBOBMX-1  is the highest-order polynomial expansion allowed for the
c      adiabatic or centrifugal Born-Oppenheimer breakdown functions, or 
c      the Lambda-doubling or 2\Sigma splitting radial strength functions
      PARAMETER (NBOBMX  = 15)
c*  HPARMX  is the largest number of Hamiltonian parameters of all types
c    (potential energy, BOB. etc.) for all states.
c           HPARMX >= NSTATEMX*[5 + (NbetaMX+1) + 5*(NBOBMX+1)]
      PARAMETER (HPARMX= NSTATEMX*(5 + (NbetaMX+1) + 5*(NBOBMX+1)))
cc    PARAMETER (HPARMX = 300)
c*  NDATAMX  is largest No. of individual data which may be considered
      PARAMETER (NDATAMX = 35000)
c*  NVIBMX  is the maximum number of vibrational levels of a single
c           state for which data are to be considered
      PARAMETER (NVIBMX    = 200)
** NBCMX  is the maximum number of band constants per vib level to be
c         allowed when doing band constant fits (PSEL= -1)
      PARAMETER (NBCMX = 8) 
c*  NPNTMX  is the largest number of potential data points that can be
c           stored in a single 1D array
      PARAMETER (NPNTMX = 90000)
c*  NROTMX  is the highest order of rotational constants calculated and
c            used for estimating level energies
      PARAMETER (NROTMX = 7)
c*  NCMMAX is the largest number of Cm terms in the MLR or DELR 
c            long-range potential
      PARAMETER (NCMMAX = 12)
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cc    INCLUDE 'BLKPOT.h'
c=======================================================================
c** Effective adiabatic radial potential variables.
      INTEGER BOBCN(NSTATEMX),PSEL(NSTATEMX),MAXMIN(NSTATEMX),
     1 IOMEG(NSTATEMX),Nbeta(NSTATEMX),APSE(NSTATEMX),IFXDE(NSTATEMX),
     2 IFXRE(NSTATEMX),IFXCm(NCMMax,NSTATEMX),
     3 IFXBETA(0:NbetaMX,NSTATEMX),NDATPT(NSTATEMX),NCMM(NSTATEMX),
     4 MMLR(NCMMax,NSTATEMX),nPB(NSTATEMX),nQB(NSTATEMX),pAD(NSTATEMX),
     5 qAD(NSTATEMX),LRad(NSTATEMX),pNA(NSTATEMX),qNA(NSTATEMX),
     6 Pqw(NSTATEMX),IVSR(NSTATEMX),IDSTT(NSTATEMX)
c
      REAL*8 DE(NSTATEMX),RE(NSTATEMX),BETA(0:NbetaMX,NSTATEMX),
     1 yqBETA(NbetaMX,NSTATEMX),BETAFX(NPNTMX,NSTATEMX),RH(NSTATEMX),
     2 RMIN(NSTATEMX),RMAX(NSTATEMX),VLIM(NSTATEMX),EPS(NSTATEMX),
     3 betaINF(NSTATEMX),AGPEF(NSTATEMX),BGPEF(NSTATEMX),
     4 CmVAL(NCMMax,NSTATEMX),CmEFF(NCMMax,NSTATEMX),rhoAB(NSTATEMX),
     5 AA(NSTATEMX),BB(NSTATEMX),RREF(NSTATEMX),ASO(NSTATEMX),
     6 R01(NSTATEMX),Q12(NSTATEMX),RD(NPNTMX,NSTATEMX),
     7 VPOT(NPNTMX,NSTATEMX),dCmA(NCMMax,NSTATEMX),dCmB(NCMMax,NSTATEMX)
c
      COMMON /BLKPOT/DE,RE,BETA,yqBETA,BETAFX,RH,RMIN,RMAX,VLIM,EPS,
     1 betaINF,AGPEF,BGPEF,CmVAL,CmEFF,rhoAB,AA,BB,RREF,ASO,R01,Q12,RD,
     2 VPOT,dCmA,dCmB, BOBCN,PSEL,MAXMIN,IOMEG,Nbeta,APSE,IFXDE,IFXRE,
     3 IFXCm,IFXBETA,NDATPT,NCMM,MMLR,nPB,nQB,pAD,qAD,LRad,pNA,qNA,Pqw,
     4 IVSR,IDSTT
c=======================================================================
cc    INCLUDE 'BLKDVDP.h'
c=======================================================================
c** Partial derivative arrays for fits and uncertainties (fununc)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      REAL*8 DVtot(HPARMX,NPNTMX),DLDDRe(NPNTMX,NSTATEMX),
     1  DUADRe(NPNTMX,NSTATEMX),DUBDRe(NPNTMX,NSTATEMX),
     2  DTADRe(NPNTMX,NSTATEMX),DTBDRe(NPNTMX,NSTATEMX),
     3  DBDB(0:NbetaMX,NPNTMX,NSTATEMX),DBDRe(NPNTMX,NSTATEMX),
     4  dVpdP(HPARMX,NPNTMX)
      COMMON/BLKDVDP/DVtot,DUADRe,DUBDRe,DTADRe,DTBDRe,DLDDRe,DBDB,
     1 DBDRe,dVpdP
c=======================================================================
cc    INCLUDE 'BLKBOB.h'
c=======================================================================
c** Born-Oppenheimer Breakdown & doubling function parameters.
c**                       March 16 2012
c=======================================================================
      INTEGER NUA(NSTATEMX),NUB(NSTATEMX),NTA(NSTATEMX),NTB(NSTATEMX),
     1  IFXUA(0:NBOBMX,NSTATEMX),IFXUB(0:NBOBMX,NSTATEMX),
     2  IFXTA(0:NBOBMX,NSTATEMX),IFXTB(0:NBOBMX,NSTATEMX),
     3  NwCFT(NSTATEMX),IFXwCFT(0:NBOBMX,NSTATEMX),efREF(NSTATEMX)
c
      REAL*8 UA(0:NBOBMX,NSTATEMX),UB(0:NBOBMX,NSTATEMX),
     1  TA(0:NBOBMX,NSTATEMX),TB(0:NBOBMX,NSTATEMX),
     2   wCFT(0:NBOBMX,NSTATEMX)
c
      COMMON /BLKBOB/UA,UB,TA,TB,wCFT,NUA,NUB,NTA,NTB,NwCFT,
     1  IFXUA,IFXUB,IFXTA,IFXTB,IFXwCFT,efREF
c=======================================================================
cc    INCLUDE 'BLKPARAM.h'
c=======================================================================
c** Parameters and count-labels for band constant (PSEL=-1) or term
c   value (PSEL=-2) fits
      REAL*8 TVALUE(NPARMX),ZBC(0:NVIBMX,0:NROTMX,NISTPMX,NSTATEMX),
     1 ZQC(0:NVIBMX,0:NROTMX,NISTPMX,NSTATEMX)
c
      INTEGER NSTATES,NTVALL(0:NSTATEMX),NTVI(NSTATEMX),NTVF(NSTATEMX),
     1 VMIN(NSTATEMX,NISTPMX),VMAX(NSTATEMX,NISTPMX),JTRUNC(NSTATEMX),
     2 EFSEL(NSTATEMX),NBC(0:NVIBMX,NISTPMX,NSTATEMX),
     3 NQC(0:NVIBMX,NISTPMX,NSTATEMX),
     4 BCPARI(0:NVIBMX,NISTPMX,NSTATEMX),
     5 BCPARF(0:NVIBMX,NISTPMX,NSTATEMX),
     6 QCPARI(0:NVIBMX,NISTPMX,NSTATEMX),
     7 QCPARF(0:NVIBMX,NISTPMX,NSTATEMX)
      COMMON /BLKPARAM/TVALUE,ZBC,ZQC,NSTATES,NTVALL,NTVI,NTVF,VMIN,
     1      VMAX,JTRUNC,EFSEL,NBC,NQC,BCPARI,BCPARF,QCPARI,QCPARF
c=======================================================================
cc    INCLUDE 'BLKISOT.h'
c=======================================================================
c** Isotope/isotopologue numbers, masses & BOB mass scaling factors
c** Array ZK carries about the band constants for all levels of all ISOT
      INTEGER NISTP,NDUNMX,AN(2),MN(2,NISTPMX)
c** NDUNMX is a dummy parameter reqd. for portability of READATA
      PARAMETER (NDUNMX=0)
      REAL*8  ZMASS(3,NISTPMX),RSQMU(NISTPMX),RSQMUP(0:NDUNMX,NISTPMX),
     1 RMUP(0:9,NISTPMX),ZMUA(NISTPMX,NSTATEMX),ZMUB(NISTPMX,NSTATEMX),
     2 ZMTA(NISTPMX,NSTATEMX),ZMTB(NISTPMX,NSTATEMX),
     3  ZK(0:NVIBMX,0:NROTMX,NISTPMX,NSTATEMX)
c
      COMMON /BLKISOT/ZMASS,RSQMU,RSQMUP,RMUP,ZMUA,ZMUB,ZMTA,ZMTB,ZK,
     1  NISTP,AN,MN
c=======================================================================
cc    INCLUDE 'BLKBOBRF.h'
c=======================================================================
c** Born-Oppenheimer breakdown radial functions 
      REAL*8 UAR(NPNTMX,NSTATEMX),UBR(NPNTMX,NSTATEMX),
     1 TAR(NPNTMX,NSTATEMX),TBR(NPNTMX,NSTATEMX),wRAD(NPNTMX,NSTATEMX)
c
      COMMON /BLKBOBRF/UAR,UBR,TAR,TBR,wRAD
c=======================================================================
cc    INCLUDE 'BLKCOUNT.h'
c=======================================================================
c    Block data file  BLKCOUNT.h
c=======================================================================
c** Counters for numbers of potential parameters of different types for 
c   each state
      INTEGER  TOTPOTPAR,POTPARI(NSTATEMX),POTPARF(NSTATEMX),
     1  UAPARI(NSTATEMX),UAPARF(NSTATEMX),UBPARI(NSTATEMX),
     2  UBPARF(NSTATEMX),TAPARI(NSTATEMX),TAPARF(NSTATEMX),
     3  TBPARI(NSTATEMX),TBPARF(NSTATEMX),LDPARI(NSTATEMX),
     4  LDPARF(NSTATEMX),HPARF(NSTATEMX),OSEL(NSTATEMX)
c
      COMMON /BLKCOUNT/TOTPOTPAR,POTPARI,POTPARF,UAPARI,UAPARF,UBPARI,
     1  UBPARF,TAPARI,TAPARF,TBPARI,TBPARF,LDPARI,LDPARF,HPARF,OSEL
c=======================================================================
cc    INCLUDE 'BLKDATA.h'
c=======================================================================
c** Type statements & common block for data
      REAL*8  FREQ(NDATAMX),UFREQ(NDATAMX),DFREQ(NDATAMX),TEMP(NDATAMX),
     1                                               YUNC(NDATAMX),Fqb
      INTEGER  COUNTOT,NFS1,NFSTOT,NBANDTOT,IB(NDATAMX),JP(NDATAMX),
     1 JPP(NDATAMX),VP(NPARMX),VPP(NPARMX),EFP(NDATAMX),EFPP(NDATAMX),
     2 TVUP(NDATAMX),TVLW(NDATAMX),FSBAND(NPARMX),IFXFS(NPARMX),
     3 NFS(NPARMX),IEP(NPARMX),IEPP(NPARMX),ISTP(NPARMX),
     4 IFIRST(NPARMX),ILAST(NPARMX),NTV(NSTATEMX,NISTPMX),FSsame,
     5 NTRANS(NPARMX),IBB(NISTPMX,NSTATEMX,9,NPARMX),JMIN(NPARMX),
     6 JMAX(NPARMX)
      CHARACTER*2 NAME(2) 
      CHARACTER*3 SLABL(-6:NSTATEMX)
      CHARACTER*30 BANDNAME(NPARMX)
      COMMON /DATABLK/Fqb,FREQ,UFREQ,YUNC,DFREQ,TEMP,COUNTOT,NFS1,
     1 NFSTOT,NBANDTOT,IB,JP,JPP,VP,VPP,EFP,EFPP,TVUP,TVLW,FSBAND,IFXFS,
     2 NFS,IEP,IEPP,ISTP,IFIRST,ILAST,NTV,FSsame,
     3 NTRANS,IBB,JMIN,JMAX,NAME,SLABL,BANDNAME
c=======================================================================
c-----------------------------------------------------------------------
      INTEGER IISTP, IV, I, J, JROT,KVLEV,efPARITY,NPTOT,NB1,COUNT1,
     1    ISTATE, NBEG,NEND,INNODE,INNER,IWR,LPRWF,WARN, fcount
      REAL*8 PV(NPARMX), PU(NPARMX), PD(NPARMX), CM(NPARMX,NPARMX), 
     1  PQ(NPARMX),PT(NPARMX),DEDPK(HPARMX),V1D(NPNTMX),SWF(NPNTMX),
     2  DVDPK(NPNTMX),Bunc1, Bunc2, Eunc2, CALC, EIV, FWHM, BFCT,BvWN,
     3  EO,Bv,UMAX,dBdPk
c** First - set up fake R(0) datum
      WARN= 1
      COUNT1= COUNTOT+1
      NB1= NBANDTOT+1
      IB(COUNT1)= NB1
      JPP(COUNT1)= 0
      JP(COUNT1)= 1
      FREQ(COUNT1)= 0.d0
      UFREQ(COUNT1)= 0.d0
      DO ISTATE= 1, NSTATES
          IF(PSEL(ISTATE).LE.0) CYCLE
          IEP(NB1)= ISTATE
          IEPP(NB1)= ISTATE
          DO IISTP= 1,NISTP
              write(7,600) 
              ISTP(NB1)= IISTP
cc            WRITE(17,600) IISTP, SLABL(ISTATE)
c** Generate 1-D potential for this state/isotopologue to prepare ...
              BFCT=RH(ISTATE)*RH(ISTATE)*ZMASS(3,IISTP)/16.857629206d0
              BvWN= 16.857629206D0/ZMASS(3,IISTP)
              DO I= 1,NDATPT(ISTATE)
                  V1D(I)= BFCT*(VPOT(I,ISTATE) 
     1                              + ZMUA(IISTP,ISTATE)*UAR(I,ISTATE)
     2                             + ZMUB(IISTP,ISTATE)*UBR(I,ISTATE))
                  ENDDO
c** for each vibrational level of each isotopologue ....
              DO  IV= VMIN(ISTATE,IISTP), VMAX(ISTATE,IISTP)
                  VP(NB1)= IV
                  VPP(NB1)= IV
c** Next get partial derivatives from Bv as half of a pure R(0) energy
                  CALL DYIDPJ(COUNT1,COUNT1,NPTOT,FREQ(COUNT1),CALC,PV,
     1                                                             PD)
c ... set up column vector for uncertainty calculation
                  DO  J= POTPARI(ISTATE), POTPARF(ISTATE)
                      PT(J)= 0.5*PD(J)*PU(J)
                      ENDDO
c... Uncertainty calculation for this level via "R(0)" approvimation 
                  CALL MMCALC(1,NPTOT,PT,CM,Bunc1)
c   
cc                WRITE(7,602)  IV,0.5d0*CALC/ZK(IV,1,IISTP,ISTATE),
cc   1                                                         Bunc1
c... Call SCHRQ to get energy and wavefunction for exact calculation
                  JROT= 0
                  KVLEV= IV
                  EO= ZK(IV,0,IISTP,ISTATE)
                  Bv= ZK(IV,1,IISTP,ISTATE)
                  CALL SCHRQ(IV,JROT,EO,FWHM,UMAX,VLIM(ISTATE),V1D,SWF,
     1                       BFCT,EPS(ISTATE),RMIN(ISTATE),RH(ISTATE),
     2                NDATPT(ISTATE),NBEG,NEND,INNODE,INNER,IWR,LPRWF)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                  efPARITY= 0
c ... call DEDP to get eigenvalue derivatives DEDPK(j) for E & Bv unc.
                  CALL DEDP(COUNT1,ISTATE,IISTP,ZMASS(3,IISTP),KVLEV,
     1                JROT,efPARITY,EO,VLIM(ISTATE),FWHM,DEDPK,fcount)
c... Now Loop over potential parameters for this state, calculating
                  DO J= POTPARI(ISTATE), HPARF(ISTATE)
                      PQ(J)= DEDPK(J)*PU(J) 
c... first - define 1D partial derivative array
                      DO I= 1,NDATPT(ISTATE)
                          DVDPK(I)= BFCT*dVtot(J,I)
                          ENDDO
cc                    CALL dPSIdp(ISTATE,IISTP,EO,NBEG,NEND,
cc                ?? problem to solve in my 'spare time'
cc   1               NDATPT(ISTATE),BvWN,V1D,SWF,DEDPK(J),dBdPk,dVdPk)
cc
cc                    WRITE(17,604)  J,0.5*PD(J),dBdPk
cc
                      PT(J)= dBdPk*PU(J)
                      CONTINUE
                      ENDDO
c... get eigenvalue Gv uncertainties
                  CALL MMCALC(1,NPTOT,PQ,CM,Eunc2)
c... get Bv uncertainties
                  CALL MMCALC(1,NPTOT,PT,CM,Bunc2)
ccc
cc                WRITE(17,610) IV,EO,Eunc2,Bv,Bunc1
                  WRITE(7,610) IV,EO,Eunc2,Bv,Bunc1
cc                WRITE(17,608) IV,Bunc1,Bunc2
ccc               WRITE(7,608) IV,Bunc1,Bunc2
                  ENDDO
              ENDDO
          ENDDO
      RETURN
  600 FORMAT(:/'  Predict Bv uncertainties for isotopologue-',I1,' in st
     1ate ',A3)
  602 FORMAT(' For  v=',I3,'    Bv{R(0)/2}/Bv(exact)=',F13.10,
     1   '       unc(Bv)=',1PD13.6)
  604 FORMAT('   for param(',I2,')    dBvdp(approx)=',1P1D13.6,
     1  '   Bvdp(numerical)=',D13.6)
  608 FORMAT(' For  v=',I3,'    u(Bv{R(0)/2})=',F13.10,
     1   '    u(Bv;numerical)=',1PD13.6)
  610 FORMAT(' v=',I3,'    E=', F12.4,'  u(E)=',F10.6,'  Bv=',f11.8,
     1  '  u(Bv)=',1PD13.6)
      END
c234567890 234567890 234567890 234567890 234567890234567890 234567890
c234567890 234567890 234567890 234567890 234567890234567890 234567890

