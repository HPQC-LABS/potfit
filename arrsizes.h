c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c** 'Block' Data Utility routine named: 'arrsizes.h' that governs 
c    array dimensioning in program  DPotFit 2.0
c-----------------------------------------------------------------------
      INTEGER NISTPMX,NPARMX,NbetaMX,NBOBMX,HPARMX,NDATAMX,
     1  NVIBMX,NBCMX,NSTATEMX,NPNTMX,NROTMX,NCMMax

c*  NISTPMX  is the maximum number of isotopomers allowed for fit
      PARAMETER (NISTPMX =  12)

c*  NSTATEMX  is maximum no. of electronic states which can be
c             simultaneously fitted to
      PARAMETER (NSTATEMX = 4)

c*  NPARMX  is the largest number of free parameters allowed for fit
c  Since FS origins may be parameters, this is also max. no, data bands
      PARAMETER (NPARMX  = 6000)

c*  NbetaMX  is the largest number of exponent parameters allowed for fit
      PARAMETER (NbetaMX  = 30)

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
      PARAMETER (NDATAMX = 33000)

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

c*  NCMMax is the largest number of Cm terms in the MLR or DELR 
c            long-range potential
      PARAMETER (NCMMax = 8)

