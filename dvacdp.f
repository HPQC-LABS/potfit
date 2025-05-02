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

