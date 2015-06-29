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
      INCLUDE 'arrsizes.h'
      INCLUDE 'BLKISOT.h'
      INCLUDE 'BLKCOUNT.h'
      INCLUDE 'BLKPOT.h'
      INCLUDE 'BLKPARAM.h'
      INCLUDE 'BLKBOB.h'
c-----------------------------------------------------------------------
c** Type statements for input or local variables
      INTEGER I, ISTATE, IISTP, m, MMN, npow, VTST
      CHARACTER*3 SLABL(-5:NSTATEMX)
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
c         = 4 : Use a Tiemann/Hannover-polynomial-potential (HPP)
c         = 5 : Use a Tang-Toennies type potential 
c         = 6 : Use an Aziz'ian HFD-C type potential 
c         = 7 : Use a Surkus Generalized Potential Energy Function (GPEF).
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
      IF((PSEL(ISTATE).GE.2).AND.(PSEL(ISTATE).LE.6)) THEN
c-----------------------------------------------------------------------
c** For MLR, DELR, HPP, TT or HFD-C potential read number of terms NCMM in 
c   the {damped} inverse-power long-range tail
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
c** For Li2(2S + 2P) states use Aubert-Frecon [PRA 55, 3458 (1997)] 
c  ULR(r) with NCMM= 6 & MMLR= {3, x, 6, 6, 8, 8} where x=0 for 
c  the A^1\Sigma_u^+ state, x=-1  for the c(1^3\Sigma_g^+) state & x=-2 for
c  the b^3\Pi_u state.
c**  For Aziz'ian HDF form (PSEL=6) input Cm's are dimensionless
c   'reduced' values.  Otherwise (PSEL< 6) units are cm-1*Angst^m
c=======================================================================
          READ(5,*) NCMM(ISTATE), rhoAB(ISTATE), IVSR(ISTATE),
     1                                                   IDSTT(ISTATE)
          DO  m= 1,NCMM(ISTATE)
              READ(5,*) MMLR(m,ISTATE), CmVAL(m,ISTATE), IFXCm(m,ISTATE)
              ENDDO
c=======================================================================
c** Note that HPP potentials have no damping
          IF(PSEL(ISTATE).EQ.4) rhoAB(ISTATE)= -1.d0
          ENDIF
c-----------------------------------------------------------------------
      IF(PSEL(ISTATE).EQ.7) THEN
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
      IF((PSEL(ISTATE).EQ.4).OR.(PSEL(ISTATE).EQ.6)) IFXDE(ISTATE)= 1
c=======================================================================
c** Read parameters defining the exponent coefficient function \beta(r)
c  NSR(s).ge.0  to use {p,q}-type exponent polynomial of order Nbeta(s)
c     if NSR(s) < 0, \beta(r) is Pashov spline defined by Nbeta(s) points
c* Nbeta(s) is order of the beta(r) exponent polynomial or # spline points
c  nPB(s)  is the power  p  for  beta(r)  basic exponent variable
c  nQB(s)  is the power  q  for  beta(r)  exponent expansion variable
c  RREF(s)  defines the reference distance in the potential exponent 
c    expansion variable:  * for  RREF.le.0 , define parameter  RREF = Re
c      * for  RREF.gt.0 , fix parameter  RREF   at its read-in value
c=======================================================================
      READ(5,*) NSR(ISTATE), Nbeta(ISTATE), nPB(ISTATE), nQB(ISTATE),
     1                                                     RREF(ISTATE)
c=======================================================================
      IF(Nbeta(ISTATE).GE.NbetaMX) THEN
          WRITE(6,755) ISTATE,Nbeta(ISTATE),NbetaMX
  755 FORMAT(/'  For ISTATE=',I2,'  read-in  Nbeta=',I3,'  while NbetaMX
     1=',I3,'  so STOP!!' )
          STOP
          ENDIF
      IF((PSEL(ISTATE).EQ.2)) THEN
          MMN= MMLR(NCMM(ISTATE),ISTATE)- MMLR(1,ISTATE)
          IF((NCMM(ISTATE).GT.1).AND.(nPB(ISTATE).LE.MMN)) THEN
              IF(MMN.LE.0) THEN
                 WRITE(6,629)
                 STOP
              ELSE
                 WRITE(6,628) nPB(ISTATE),MMN
                 ENDIF
              ENDIF
          ENDIF
      IF(((PSEL(ISTATE).eq.3).AND.(NSR(ISTATE).lt.0))) 
     1                                      NSR(ISTATE)= Nbeta(ISTATE)
      IF(PSEL(ISTATE).EQ.5) Nbeta(ISTATE)= 0
      IF(PSEL(ISTATE).EQ.6) Nbeta(ISTATE)= 4
      IF((NSR(ISTATE).GE.0).OR.(PSEL(ISTATE).NE.5)) THEN
c** For NSR .ge.0  the MLR exponent  \beta(yp)  is  'conventional'
c  Huang-type constrained polynomial in  yp with polynomial order  NSR
c   for r < r_e and  NBETA for r.ge.r_e
          npow= MAX(NSR(ISTATE),Nbeta(ISTATE))
          IF(PSEL(ISTATE).EQ.4) npow= npow+3
          IF(Nbeta(ISTATE).GE.0) THEN
              DO  I= 0,npow
c** Read in trial initial trial parameters for exponent \beta(r)
c
c  BETA(i,s)   contains the expansion parameters defining the potential
c    for  PSEL.LE.3 : read-in values are the {npow+1} beta_i  exponent
c                   exponent expansion parameters defining the potential
c    for  PSEL = 4  : read in the  {1+Nbeta} expansion parameters plus
c                         b, RINN, and ROUT  of the HPP form
c** For PSEL=6: >> set Nbeta=4 to use the single global damping function for 
c   HDFD-A,B, & C potentials:  f_1(x)= beta(0)*exp{-\beta(1)/x)^beta(2)}
c     while exponent is {\alpha*x + beta(3)*x^2} and  \gamma= beta(4).
c**             >> set Nbeta=3 to combine the overall damping function 
c      f_2(x)= [1 - r^{\beta(0)} exp{-\beta(1)(r}]  , and in this case
c     while exponent is {\alpha*x + beta(2)*x^2} and  \gamma= beta(3).
c**  for  PSEL = 7  : read-in values are leading coefficients in 
c                Surkus' Generalized Potential Energy Function (GPEF).
c  IFXBETA(i,s) indicates whether each potential expansion coefficient
c             coefficient will be:    = 1: held fixed at read-in values.
c                                  .LE. 0: determined from fits.
c=======================================================================
                  READ(5,*) BETA(I,ISTATE), IFXBETA(I,ISTATE)
c=======================================================================
                  ENDDO
              ENDIF    
        ELSE
          npow= Nbeta(ISTATE)
          DO  I= 1,npow
c-----------------------------------------------------------------------
c** For NSR .lt 0  exponent is a natural spline function with values BETA
c       at the yp values ypBETA, and fixed to equal ypINF at  ypBETA=1
c=======================================================================
              READ(5,*)ypBETA(I,ISTATE),BETA(I,ISTATE),IFXBETA(I,ISTATE)
c=======================================================================
              ENDDO
          Nbeta(ISTATE)= Nbeta(ISTATE)+1
          ypBETA(Nbeta(ISTATE),ISTATE)= 1.d0
          IFXBETA(Nbeta(ISTATE),ISTATE)= 1
        ENDIF
      IF(PSEL(ISTATE).EQ.4) THEN
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
      READ(5,*) NUA(ISTATE),NUB(ISTATE),pAD(ISTATE),qAD(ISTATE),
     1                                                    LRad(ISTATE)
      IF(((NUA(ISTATE).GE.0).OR.(NUB(ISTATE).GE.0))
     1              .AND.(PSEL(ISTATE).EQ.1)) pAD(ISTATE)= qAD(ISTATE)
c... NOTE never read delta Cm values unless PSEL = 2-6
          IF((PSEL(ISTATE).LT.2).OR.(PSEL(ISTATE).GT.6)) LRad(ISTATE)=0
          IF(LRad(ISTATE).EQ.1) THEN
c... Now read \delta{Cm} values dCmA & dCmB for atoms-A & B, desired
              READ(5,*) (dCmA(m,ISTATE), m=1,NCMM(ISTATE))     
              READ(5,*) (dCmB(m,ISTATE), m=1,NCMM(ISTATE))     
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
      READ(5,*) NTA(ISTATE), NTB(ISTATE), pNA(ISTATE), qNA(ISTATE)
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
     2 behaviour"/(1x,20('****')))
  629 FORMAT(' ** Since  [MMLR(NCMM)-MMLR(1)].le.0,   STOP !!')
  640 FORMAT(/' ', A3,' state energies referenced to f-parity levels')
  642 FORMAT(/' ', A3,' state energies referenced to the mid-point betwe
     1en e and f-parity levels')
  644 FORMAT(/' ', A3,' state energies referenced to e-parity levels')
  646 FORMAT(/' *** INPUT ERROR ***  |efREF=',i3,'| > 1')
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
