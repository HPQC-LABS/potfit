      INTEGER MXPARM,MXDATA,NTP
      PARAMETER (MXPARM=12,MXDATA=500,NTP=30000)
      real*8  yl0(0:20),yl1(0:20),Evj(mxdata),unc(mxdata),
     1  YO(mxdata),dY(mxdata),dydp(mxdata,mxparm),pv(mxparm),
     2  pu(mxparm),ps(mxparm),cm(mxparm,mxparm),DCMPAR(mxparm),
     3  RPT,VT(NTP),VJ(NPT),dVdP(NPT,mxparm)
      integer i,j,idat,ndata,npar,max(0:4),iv(500),ij(500)
c** ndata   id the number of E(v,J) data to read in and fit to
c** lmax(0)  is the power of the Dunham expansion for G_v used to
c        generate trial energues:  Yl0(i), starting with  Te
c** lmax(1)  is the power of the Dunham expansion for B_v used to
c        generate trial energues:  Yl1(i)
c** Evj(iv(i),ij(i))  are the input (binding) ennergies to be 
c       fitted, and   unc(i)  their individual uncertainties.
c-----------------------------------------------------------------
      read(5,*) RH,RMIN
      read(5,*) ndata,lmax(1),lmax(2)
      read(5,*) (Yl0(i),i=1,lmax(0))
      read(5,*) (Yl1(i),i=1,lmax(1))
      read(5,*) (iv(i),ij(i),Evj(i),unc(i),i=1,ndata)
c-----------------------------------------------------------------
      write(6,600) (Yl0(i),i=1,lmax(0))
      write(6,602) (Yl1(i),i=1,lmax(1))
  600 format(' Generate trial energies from',i3,' Dunham vibrationa
     1l coefficients'/(4x,1pD15.7,4d15.7))
  602 format(5x,' and Dunham rotational coefficients'/
     2  (4x,1pD15.7,4d15.7))
      do  i= 1,npt
          rpt(i)= rmin+ (i-1)*RH
          enddo
      write(6,604) RMIN,RPT(npt),RH
  604 format(/' Integration range is from',f6.2,'  to',f7.2,
     1  '  with mesh   RH=',f7.4)
      write(6,606) 
  606 format(/' Fit to the',I4,' input binding energies (uncertaint
     1ies)'/(2('   E(v=',I3,' J=',i2,')=',f8.4,' (+/-',f6.3,')'))) 
c
c** End of data input ... READ initial trial potential parameters
c-----------------------------------------------------------------
      read(5,*)  (DCMPAR(j),j=1,10)
      write(6,610)  (DCMPAR(j),j=1,10)
  610 format(/' Initial trial potential parameters are:'/
     1  (1p5d16.7))
c==================================================================
c** End of input ... begin non-linear fitting loop
c==================================================================
      do  itry= 1, 50
c
c** CALL subroutine with current parameter set DCMPAR to return 
c   trial potential function array  VT(i)  and the  10  partial 
c   derivative columns of the array  dVdP(i,j) (for  j=1,10) at the 
c   NPT  radial mesh points  RPT(i).
          CALL VPOT(NPT,RPT,DCMPAR,VT,dVdP)
c ... now convert potential to SCHRQ's correct internal units ...
          do  i= 1, NPT
              VT(i)= VT(i)*BFCT
              enddo
c** Perform loop over all data to generate current calculated
c  energy  ETRY  and partial derivatives required by  LLSQF
c
          do  idat= 1, ndata
c ... first, calvulate trial energy from the Dunham coefficients
c
              etry= xxxxxxxx   ... etc.

c ... now call SCHRQ to get eigenvalue and eigenfunction
c ...... after first generating centrifugally distorted potential
              do  i= 1, NPT
                  VJ(i)= VT(i)+ EJJ/RPT(i)**2
                  enddo
              CALL SCHRQ( ;;;;; )
c** Now ... zero partial derivative array ...
              do  j= 1,npar
                  dydp(idat,j)= 0.d0
                  enddo
c ... and then calculate partial derivatives required for  LLSQF
              do  i= 1, NPT
                  ps2= psi(i)**2
                  do  j= 1, npar
                      dydp(idat,j)= dydp(idat,j)+ ps2*dVdP(i,j)
                      enddo
                  enddo
              enddo
          call LLSQF(ndata,npar,mxdata,mxparm,Evj,unc,dydp,dy,
     1                                         PV,PU,PS,CM,SERR)
c ... now update potential parameters & test for convergence
          tstps= 0.d0
          tstpu= 0.d0
          do  j= 1, npar
              tstpu= dmax1(tstpu,dabs(pv(j)/pu(j))
              tstps= dmax1(tstps,dabs(pv(j)/ps(j))
              DCMPAR(j)= DCMPAR(j)+ pv(j)
              enddo
          if(tstps.lt.1.d0) goto 88
          enddo
c=================end of iterative fitting loop===================
      write(6,620) itry,tstps,tstpu
  620 format(' *** NOT CONVERGED *** After',i3,' cycles  tstPS)=',
     1  1pd8.1,'   tst(PU)=',d8.1)
   88 










c***********************************************************************
c**** R.J. Le Roy  subroutine SCHRQ, last updated  6 February 1999 *****
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c                COPYRIGHT 1998  by  Robert J. Le Roy                  +
c   Dept. of Chemistry, Univ. of Waterloo, Waterloo, Ontario, Canada   +
c    This software may not be sold or any other commercial use made    +
c      of it without the express written permission of the author.     +
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** SCHRQ solves radial Schrodinger equation in dimensionless form
c  d2S/dR2 = - (E-V(R))*S(R) ,  where S(I) is the wave function.
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
c** For energy in (cm-1), BFCT=ZMU(u)*H(Angst)**2/16.8576314 (1/cm-1)
c** INNER specifies wave function matching (& initiation) condition *
c*  For INNER = 0 , match inward & outward solutions at outermost
c  wave function maximum;  otherwise match at inner edge of classically
c  allowed region. ** INNER<0  uses zero slope inner boundary condition.
c** For most normal cases set INNER=0 ,  but ......
c* to find "inner-well-dominated" solutions of an asymmetric double
c  minimum potential, set  INNER > 0 .
c* To find symmetric eigenfunctions of a symmetric potential, set
c  INNER < 0  & start integration (set RMIN) at potential mid point.
c**********************************************************************

      SUBROUTINE SCHRQ(KV,JROT,EO,GAMA,VMAX,VLIM,V,S,BFCT,EEPS,RMIN,RH,
     1  N,NBEG,NEND,INNER,IWR,LPRWF)

c** Output vibrational quantum number KV, eigenvalue EO, normalized
c  wave function S(I), and range, NBEG .le. I .le. NEND  over
c  which S(I) is defined. *** Have set  S(I)=0  outside this range.

c* (NBEG,NEND), defined by requiring  abs(S(I)) < RATST=1.D-9  outside.

c** If(LPRWF.gt.0) print wavefunction S(I) every LPRWF-th point.
c* If(LPRWF.lt.0) "punch" (i.e., WRITE(8,XXX)) every |LPRWF|-th point
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
c  mesh point |KV| , based on incoming value of wave function S(|KV|)
c  and of the wavefunction derivative at that point, SPNEND, which is
c  brought in as S(|KV|-1).  For a hard wall condition at mesh point
c  |KV|, set S(|KV|)=0 and S(|KV|-1)= -1 before entry.
c**********************************************************************
c++ "SCHRQ" calls subroutineas "QBOUND" and "WIDTH", and the latter
c++ calls "LEVQAD" .
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION  RWR(20),SWR(20),S(N),V(N)
      DATA Z0/0.D0/,ZH/0.5D0/,Z1/1.D0/,Z2/2.D0/,Z48/48.D0/

c** To minimize effort, choose inward integration starting point for
c  truly bound state (at mesh point NEND) such that normalized wave
c  function there is  .lt. ca.  RATST = exp(-XPW)

      DATA RATST/1.D-9/,XPW/20.7D0/
      IR4= 0
      NDN= 20
      KVIN= KV
      KV= -1
      RMINN= RMIN-RH
      DXPW= (XPW+2.3D0)/NDN
      GAMA= Z0
      VMAX= VLIM
      VMX= VMAX*BFCT
      H= RH
      H2= H*H
      HT= Z1/12.D+0
      E= EO*BFCT
      EPS= EEPS*BFCT
      DSOC= VLIM*BFCT
      DE= Z0
      RATIN= Z0
      RATOUT= Z0
      IF(IWR.GT.2) THEN
          IF(KVIN.GE.998) then
              WRITE(6,610) EO
            ELSE
              WRITE(6,601) KVIN,JROT,EO
            ENDIF
          WRITE(6,602)
        ENDIF
      NEND= N
      IF(KVIN.LT.-10) THEN
          NEND= -KVIN
          SNEND= S(NEND)
          SPNEND= S(NEND-1)
          ENDIF
      JQTST = 0
c** Start iterative loop; try to converge for up to 15 iterations.
      DO 90 IT= 1,15
      ITER = IT
      IF(INNER.NE.0) GO TO 38
    8 IF(KVIN.LT.-10) THEN

c** If desired, (KVIN < -10) outer boundary set at NEND=|KVIN| and 
c  initialize wavefunction with log-derivative condition based on value
c  SNEND and derivative SPNEND at that mesh point.

                GN= V(NEND)-E
                GI= V(NEND-1)-E
                SB= SNEND
                SI= SB*(Z1+ ZH*GN)- RH*SPNEND
                GO TO 24
            END IF
      IF(E.GE.DSOC) GO TO 20

c** For truly bound state attempt to minimize computation by starting
c  outer wave function initialization at smallest distance for which
c  still have    RATOUT .le. RATST = exp(-XPW)
c** First estimate JWKB exponent for start at end of range.

      IF(ITER.GT.3) GO TO 18
      NEND= N
      M= NEND
   10 GI= V(M)-E

c** First do rough search for outer turning point
      IF(GI.LE.Z0) GO TO 12
      GN= GI
      IR4= M
      M= M-NDN
      IF(M) 18,18,10

c** Linear estimate for quadrature near 2-nd turning point.
   12 IF(M.GE.N) GO TO 998
      FX= GN/(GI-GN)
      SM= ZH*(Z1+FX)*DSQRT(GN)
      MS= M+2*NDN
      DO 14 M= MS,N,NDN
      NEND= M
      SM= SM+DSQRT(V(M)-E)
      IF(SM.GE.DXPW) GO TO 16
   14 CONTINUE
   16 IF(NEND.LT.N) NEND= NEND+NDN
   18 GN= V(NEND)-E
      GI= V(NEND-1)-E
      M= NEND-1
      IF(GI.GE.Z0) GO TO 22
      IF(NEND.GE.N) GO TO 998
      NEND= N
      GO TO 18

c** For quasibound levels, initialize wave function in "QBOUND"
   20 CALL QBOUND(KVIN,JROT,E,EO,VMX,DSOC,V,RMIN,H,GN,GI,SB,SI,N,ITP3,
     1  IWR,IQTST,BFCT,IT)
      NEND= ITP3
      VMAX= VMX/BFCT
      IR4= ITP3
      IF(IQTST) 21,22,24
   21 JQTST = JQTST+IQTST
      IF((JQTST.LE.-2).OR.(VMAX.LT.VLIM)) GO TO 999

c** Try up to once to find level using trial value just below maximum
      EO = VMAX-0.5D0
      E = EO*BFCT
      GO TO 90

c* Initialize wave function as 1-st order WKB solution increasing inward
   22 SRTGN= DSQRT(GN)
      SRTGI= DSQRT(GI)
      SB= Z1
      SI= SB*DSQRT(SRTGN/SRTGI)*DEXP((SRTGN+SRTGI)/Z2)
      IF(SB.LE.SI) GO TO 24

c** WOOPS ... something slipped up so initialize with a node.
      IF(IWR.NE.0) WRITE(6,618) JROT,EO,SB,SI
      SI= Z1
      SB= Z0
   24 M= NEND-1
      Y1= (Z1-HT*GN)*SB
      Y2= (Z1-HT*GI)*SI
      S(NEND)= SB
      S(NEND-1)= SI
      M2= NEND
      NENDCH= NEND
      IBEGIN= 3
      IF(INNER.NE.0) IBEGIN= ITP1+2

c** Actual inward integration loop starts here
      DO 30 I= IBEGIN,NEND
      M= M-1
      Y3= Y2+Y2-Y1+GI*SI
      GI= V(M)-E
      SB= SI
      SI= Y3/(Z1-HT*GI)
      S(M)= SI
      IF(DABS(SI).GE.1.D+17) THEN

c** Renormalize to prevent overflow of  S(I)  in classically
c  forbidden region where  (V(I) .gt. E)
          SI= Z1/SI
          DO 26 J= M,M2
   26     S(J)= S(J)*SI
          NENDCH= M2
          M2= M
          Y2= Y2*SI
          Y3= Y3*SI
          SB= SB*SI
          SI= Z1
        ENDIF

c** Test for outermost maximum of wave function.
      Y1= Y2
      Y2= Y3
      IF(INNER.NE.0) GO TO 30
      IF(SI.LE.SB) GO TO 32
   30 CONTINUE
      IF(INNER.NE.0) GO TO 32

c** Error mode ... find no wave function maximum.
      KV= -2
      IF(IWR.NE.0) WRITE(6,616) KV,JROT,EO
      GO TO 999

c** Scale outer part of wave function before proceding
   32 SI= Z1/SI
      MSAVE= M
      RR= RMINN+MSAVE*H
      YIN= Y1*SI
      RATOUT= S(NEND)*SI
      NEND= NENDCH
      DO 34 J= MSAVE,NEND
   34 S(J)= S(J)*SI
      IF(INNER.NE.0) GO TO 70
c-------------------------------------------------------------------
c** Set up to prepare for outward integration **********************
   38 NBEG= 1
      IF(INNER.LT.0) THEN

c** Option to initialize with zero slope at beginning of the range
          SB= Z1
          GN= V(1)-E
          Y1= SB*(Z1-HT*GN)
          Y2= Y1+GN*SB/Z2
          GI= V(2)-E
          SI= Y2/(Z1-HT*GI)
        ELSE

c** Initialize outward integration with a node at beginning of range
   40     GN= V(NBEG)-E
          IF(GN.LT.10.D0) GO TO 42

c** If potential has (V(1)-E) so high that H is (locally) much too
c  large, then shift inner starting point outward.
          NBEG= NBEG+1
          IF(NBEG.LT.N) GO TO 40
          IF(IWR.NE.0) WRITE(6,613)
          GO TO 999
   42     IF((ITER.GT.1).OR.(IWR.EQ.0)) GO TO 44
          IF(NBEG.GT.1) WRITE(6,609) JROT,EO,NBEG
          IF(GN.LE.Z0) WRITE(6,604) JROT,EO,E,V(NBEG),NBEG

c** Initialize outward wave function with a node:  S(NBEG) = 0.
   44     SB= Z0
          SI= Z1
          GI= V(NBEG+1)-E
          Y1= SB*(Z1-HT*GN)
          Y2= SI*(Z1-HT*GI)
      END IF
c
      S(NBEG)= SB
      S(NBEG+1)= SI
      NBEGB= NBEG
      NBEG2= NBEG+2
      IF(INNER.NE.0) MSAVE= N

c** Actual outward integration loops start here
      DO 50 I= NBEG2,MSAVE
      Y3= Y2+Y2-Y1+GI*SI
      GI= V(I)-E
      SI= Y3/(Z1-HT*GI)
      S(I)= SI
      IF(DABS(SI).GE.1.D+17) THEN

c** Renormalize to prevent overflow of  S(I)  in classically forbidden
c  region where  V(I) .gt. E
          SI= Z1/SI
          NBEG= NBEGB
          DO 46 J= NBEG,I
   46     S(J)= S(J)*SI
          NBEGB= I
          Y2= Y2*SI
          Y3= Y3*SI
          SI= Z1
        ENDIF
      Y1= Y2
      Y2= Y3
      ITP1= I
      IF(GI.LE.Z0) GO TO 52

c** Exit from this loop at onset of classically allowed region
   50 CONTINUE
      MM= MSAVE
      IF((INNER.EQ.0).AND.(GN.LE.Z0)) GO TO 60
      IF(IWR.NE.0) WRITE(6,612) KVIN,JROT,EO,MSAVE
      GO TO 999
   52 ITP1P= ITP1+1
      MM= ITP1
      IF(INNER.NE.0) GO TO 60
      DO 58 I= ITP1P,MSAVE
      Y3= Y2+Y2-Y1+GI*SI
      GI= V(I)-E
      SI= Y3/(Z1-HT*GI)
      S(I)= SI
      IF(DABS(SI).LE.1.D+17) GO TO 56

c** Renormalize to prevent overflow of  S(I) , as needed.
      SI= Z1/SI
      NBEG= NBEGB
      DO 54 J= NBEG,I
   54 S(J)= S(J)*SI
      NBEGB= I
      Y2= Y2*SI
      Y3= Y3*SI
      SI= Z1
   56 Y1= Y2
   58 Y2= Y3
      MM= MSAVE

c** Finished outward integration.  Normalize w.r.t. S(MSAVE)
   60 SI= Z1/SI
      YOUT= Y1*SI
      YM= Y2*SI
      RATIN= S(NBEG+1)*SI
      DO 62 I= NBEG,MM
   62 S(I)= S(I)*SI
      IF(INNER.NE.0) GO TO 8

c----- Finished numerical integration ... now correct trial energy
c** DF*H  is the integral of  (S(I))**2 dR
   70 DF= Z0
      DO 72 J= NBEG,NEND
   72 DF= DF+S(J)**2

c** Add edge correction to DF assuming wave function dies off as simple
c  exponential past R(NEND);  matters only if S(NEND) unusually large.
      IF((E.GT.DSOC).OR.(S(NEND).EQ.0)) GO TO 73
      IF((KVIN.GE.-10).AND.(S(NEND-1)/S(NEND).GT.Z1))
     1    DF= DF+ S(NEND)**2/(Z2*DLOG(S(NEND-1)/S(NEND)))
   73 F= (-YOUT-YIN+Z2*YM+GI)
      DOLD= DE
      IF(DABS(F).LE.1.D+30) GO TO 74
      F= 9.9D+30
      DF= F
      DE= DABS(0.01D+0 *(DSOC-E))
      GO TO 76
   74 DE= F/DF
   76 IF(IWR.LE.2) GO TO 80
      DEPRN = DE/BFCT
      XEND= RMINN+NEND*H

c** RATIN & RATOUT  are wave fx. amplitude at inner/outer ends of range
c  relative to its value at outermost extremum.
      WRITE(6,603) IT,EO,F,DF,DEPRN,MSAVE,RR,RATIN,RATOUT,XEND,NBEG,ITP1
c** KV.ge.998  Option ... Search for highest bound level.  Adjust new
c  trial energy downward if it would have been above dissociation.
   80 IF((KVIN.GE.998).AND.(E.GT.DSOC)) E= E/Z2
      EO= E/BFCT

c** Test trial eigenvalue for convergence
      IF(DABS(DE).LE.DABS(EPS)) GO TO 100
      E= E+DE
      EO= E/BFCT
      IF(IT.LE.4) GO TO 90
      IF(DABS(DE).LT.DABS(DOLD)) GO TO 90

c** Adjust energy increment if having convergence difficulties.  Not
c  usually needed except for some quasibounds extremely near  VMAX .
      IF((DOLD*DE).GT.Z0) GO TO 90
      ICOR= ICOR+1
      DEP= DE/BFCT
      IF(IWR.NE.0) WRITE(6,617) IT,DEP
      DE= ZH*DE
      E= E-DE
      EO= E/BFCT
   90 CONTINUE
      E= E-DE
      EO= E/BFCT
      DEPRN= DE/BFCT
      IF(IWR.NE.0) WRITE(6,620) KVIN,JROT,ITER,EO,DEPRN

c** End of iterative loop which searches for eigenvalue ************
c-------------------------------------------------------------------*
  100 IF(IWR.EQ.0) GO TO 104
      IF(IWR.GE.3) WRITE(6,619)
      IF((DABS(RATIN).GT.RATST).AND.(INNER.GE.0))
     1         WRITE(6,614) JROT,EO,RATIN
      IF((E.LT.DSOC).AND.(DABS(RATOUT).GT.RATST)) THEN
          WKBTST= ZH*DABS(V(NEND)-V(NEND-1))/DSQRT((V(NEND)-E)**3)
          IF(WKBTST.GT.1.d-3) WRITE(6,615) JROT,EO,RATOUT,RATST,WKBTST
          ENDIF
  104 KKV = 0

c** Perform node count on converged solution
      PROD= S(ITP1)*S(ITP1-1)
      J1= ITP1+1
      J2= NEND-1
      DO 110 J= J1, J2
      PPROD= PROD
      PROD= S(J)*S(J-1)
      IF((PPROD.LE.Z0).AND.(PROD.GT.Z0)) KKV= KKV+1
  110 CONTINUE
      KV = KKV

c** Normalize & find interval (NBEG,NEND) where S(I) is non-negligible
      SN= Z1/DSQRT(H*DF)
      DO 116 I= NBEG,NEND
  116 S(I)= S(I)*SN
      IF(ITP1.LE.1) GO TO 122
      J= ITP1P
      DO 118 I= 1,ITP1
      J= J-1
      IF(DABS(S(J)).LT.RATST) GO TO 119
  118 CONTINUE
  119 NBEG= J
      IF(NBEG.LE.1) GO TO 122
      J= J-1
      DO 120 I= 1,J
  120 S(I)= Z0
  122 IF(KVIN.GE.-10)THEN

c** For "non-wall" cases, move NEND inward to where wavefunction 
c  "non-negligible"
          J= NEND-1
          DO 124 I= NBEG,NEND
          IF(DABS(S(J)).GT.RATST) GO TO 126
  124     J= J-1
  126     NEND= J+1
          END IF
      IF(NEND.LT.N) THEN

c** Zero out wavefunction array at distances past NEND
                NENDP= NEND+1
                DO 128 I= NENDP,N
  128           S(I)= Z0
             END IF
      IF(LPRWF) 133,140,135

c** If desired, "punch" (i.e., WRITE(8,XXX)) every |LPRWF|-th point 
c  wave function value onto disk, starting at the NBEG-th mesh point.
  133 JPSIQ= -LPRWF
      NPR= 1+(NEND-NBEG)/JPSIQ
      RINC= RH*JPSIQ
      RSTT= RMINN+NBEG*RH

c** Write every JPSIQ-th point of the wave function for level  v=KV
c  J=JROT , beginning at mesh point NBEG & distance RSTT where
c  the NPR values written separated by mesh step RINC=JPSIQ*RH
      WRITE(8,701) KV,JROT,EO,NPR,RSTT,RINC,NBEG,JPSIQ
      WRITE(8,702) (RMINN+I*RH,S(I),I=NBEG,NEND,JPSIQ)
      GO TO 140

c** Print solutions every  LPRWF-th  point, 6 to a line, in columns.
  135 NLINES= ((1+(NEND-NBEG)/LPRWF)+3)/4
      IPSID= LPRWF*NLINES
      WRITE(6,605) KV,JROT,EO
      DO 138 J= 1,NLINES
      JJ= NBEG+(J-1)*LPRWF
      IJK= 0
      DO 136 IJ= JJ,NEND,IPSID
      IJK= IJK+1
      RWR(IJK)= RMINN+IJ*H
  136 SWR(IJK)= S(IJ)
  138 WRITE(6,606) (RWR(I),SWR(I),I= 1,IJK)
  140 IF(IWR.EQ.1) WRITE(6,607) KV,JROT,EO
      IF(IWR.GE.2) WRITE(6,607) KV,JROT,EO,ITER,RR,RATIN,RATOUT

c** For quasibound levels, calculate width in subroutine "WIDTH"
      IF((E.GT.DSOC).AND.(KVIN.GT.-10)) CALL WIDTH(KV,JROT,E,EO,DSOC,
     1  V,S,VMX,RMIN,H,BFCT,IWR,ITP1,MSAVE,ITP3,INNER,N,GAMA)
      RETURN

c** ERROR condition if  E.gt.V(R)  at outer end of integration range.
  998 XPR= RMINN+M*H
      VPR= V(M)/BFCT
      IF(IWR.NE.0) WRITE(6,608) EO,M,VPR,XPR,IT

c** Return in error mode
  999 KV= -1
      RETURN

  601 FORMAT(/' Solve for  v=',I3,'   J=',I3,'   ETRIAL=', G14.7,
     1   24x,'S(1st)   S(NEND)' )
  602 FORMAT(' ITER    ETRIAL',10X,'F(E)',8X,'DF(E)',7X,'D(E)',
     1 7X,'M    R(M)     /S(M)     /S(M)  R(NEND) NBEG ITP1'/
     2  1X,53('--'))
  603 FORMAT(I4,G16.8,1P3D12.4,I6,0PF8.4,1P2D10.1,0PF8.4,2I5)
  604 FORMAT('      NOTE:  for   J =',I3,'   EO =',F12.4,'   E=',D13.6,
     1 ' .ge. V(R)=',D13.6,'   at initial mesh point',I6)
  605 FORMAT(/' Solution of radial Schr. equation for   E(v=',I3,',J=',
     1  I3,') =',G14.7/2x,4('    R(I)    S(I)   ')/2X,38('--') )
  606 FORMAT(2X,4(F8.3,F11.7))
  607 FORMAT(' E(v=',I3,', J=',I3,')=',F10.3,2X,I4,' Iterations'
     1 ,'   R(M)=',F7.4,'   S(NBEG)/S(Max)=',1PD8.1,'   S(NEND)/S(Max)='
     2  ,D8.1)
  608 FORMAT(' *** ERROR *** At outer limit ...  E =',F12.4,' .GT. V(',
     1  I5,')=',F12.4,'   at   R =',F9.4,'   for   IT =',I2)
  609 FORMAT(' *** For   J =',I3,'    E =',F12.4,',  integration cannot'
     1 ,' start till past mesh point #',I4,',  so RMIN smaller than need
     2ed')
  610 FORMAT(/' Attempt to find the highest bound level starting from',
     1 '   ETRIAL =',1PD9.2)
  612 FORMAT(/' *** ERROR *** for   v =',I3,'   J =',I3,'   E =',
     1  G13.6,'  Innermost turning point not found by   M = MSAVE =',I5)
  613 FORMAT(/' *** ERROR in potential array ... V(I) everywhere',
     1 ' too big to integrate with given  increment')
  614 FORMAT(' *** CAUTION *** For  J=',I3,'  E=',G15.8/18x,
     1 'S(first)/S(Max)=',D9.2,'  suggests  RMIN  may be too large')
  615 FORMAT(' ** CAUTION ** For  J=',I3,'  E=',1PD13.6,
     1 '   S(NEND)/S(Max)=',D8.1,' >',D8.1/4X,'& initialization ',
     2 'quality test ',1PD8.1,' > 1.D-3   so RMAX may be too small')
  616 FORMAT(' ** WARNING *** For  v=',I2,', J=',I3,' at  E=',
     1 G14.7,', S always has negative slope ... Energy too low or potent
     2ial too weak' )
  617 FORMAT(5X,'*** SCHRQ convergence problem so for   IT =',I2,
     1 '   cut   DE =',G11.4,'  in HALF' )
  618 FORMAT(' *** WARNING *** For   J =',I3,'    E =',F9.3,'    JWKB in
     1itialization gives   SB/SI =',D10.3,'/',D10.3,',  so set   SB=0.')
  619 FORMAT(1X,54('--'))
  620 FORMAT('  *** CAUTION ***   v =',I2,'   J =',I3,'    no convergenc
     1e after',I4,' tries.   E =',G14.7,'   DE =',G14.7)
  701 FORMAT(/2x,'Level  v=',I3,'   J=',I3,'   E=',F12.4,' ,  wave funct
     1ion at',I6,' points.'/7x,'R(1-st)=',F12.8,'   mesh=',F12.8,
     2  '   NBEG=',I4,'   |LPRWF|=',I3)
  702 FORMAT((1X,4(f9.4,f10.6)))
      END












c***********************************************************************
      SUBROUTINE LLSQF(NDATA,NPARM,MXDATA,MXPARM,YO,UY,DYDP,DY,P,UP,PS,
     1                 CM,SERR)
c
c
c**  Program for performing linear least squares fits using orthogonal 
c  decomposition of the Design (partial derivative) matrix.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c                COPYRIGHT 1993  by  Robert J. Le Roy                  +
c   Dept. of Chemistry, Univ. of Waterloo, Waterloo, Ontario, Canada   +
c    This software may not be sold or any other commercial use made    +
c      of it without the express written permission of the author.     +
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c** This version of the program is designed for the data sets of modest
c  size where it is convenient to generate and store the complete 
c  partial derivative matrix prior to calling LLSQF.  If this is not the
c  case, subroutine version LLSQFVL, which generates this partial 
c  derivative array one row at a time through calls to a user-supplied
c  subroutine, should be used.
c*************************************************************************
c
c
c** On entry: NDATA  is the number of data to be fitted (.le.MXDATA)
c
c             NPARM  the number of parameters to be varied (.le.MXPARM)
c                 If NPARM.le.0 , set  DY(i)=YO(i)  and calculate the 
c                 (weighted RMS deviation)=SERR  from them & the UY(i)'s
c
c             MXDATA & MXPARM are array dimension parameters (see below)
c                 Internal array sizes currently assume  MXPARM .le. 60
c
c             YO(i)  are the NDATA 'observed' data;  for iterative 
c                  non-linear fits these are:  [Y(obs,i) - Y(trial,i)]
c
c             UY(i)  are the uncertainties in these YO(i) values
c
c             DYDP(i,j)  is the partial derivative array  dYO(i)/dP(j)
c
c** On Exit: P(j)  are the fitted parameter values;  for iterative
c                  non-linear fits these are the parameter changes
c
c            UP(j) are 95% confidence limit uncertainties in the P(j)'s
c
c            PS(j) are 'parameter sensitivities' for the P(j)'s, defined
c               such that the RMS displacement of predicted data  due to
c               rounding off parameter-j by PS(j) is .le. SERR/10*NPARM
c
c            SERR  is the standard error of the fit
c
c            DY(i) is the array of differences  [YO(i) - Ycalc(i)]
c
c            CM(j,k)  is the correlation matrix obtained by normalizing
c   variance/covariance matrix:  CM(j,k) = CM(j,k)/SQRT[CM(j,j)*CM(k,k)]
c
c** The squared 95% confidence limit uncertainty in a property F({P(j)})
c  defined in terms of the fitted parameters {P(j)} is (where the
c  L.H.S. involves  [row]*[matrix]*[column]  multiplication):
c  [D(F)]^2 = [UP(1)*dF/dP(1), UP(2)*dF/dP(2), ...]*[CM(j,k)]*
c                              [UP(2)*dF/dP(1), UP(2)*dF/dP(2), ...]
c
c** Externally dimension:  YO, UY and DY  .ge. NDATA (say as MXDATA),
c             P, UP  and  PS  .ge.  NPARM (say as MXPARM), 
c             DYDP  with column length MXDATA and row length .ge. NPARM
c             CM   as square matrix with column length  MXPARM
c
c  Authors: Michael Dulick  &  Robert J. Le Roy, Department of Chemistry
c    U. of Waterloo, Waterloo, Ontario  N2L 3G1.    Version of:  1/09/93
c***********************************************************************
c
c
      INTEGER I,J,K,L,M,IDF,NDATA,MXDATA,NPARM,MXPARM
      REAL*8  YO(NDATA), UY(NDATA), DY(NDATA), P(NPARM), UP(NPARM), 
     1   PS(NPARM), DYDP(MXDATA,NPARM), CM(MXPARM,MXPARM), SERR,
     2   PX(60), F95(10), TFACT, S, U
      DATA F95/12.7062D0,4.3027D0,3.1824D0,2.7764D0,2.5706D0,2.4469D0,
     1  2.3646D0,2.3060D0,2.2622D0,2.2281D0/
c
      IF((NDATA.GT.MXDATA).OR.(NPARM.GT.MXPARM).OR.(NPARM.GT.60)
     1                    .OR.(NPARM.GT.NDATA)) GOTO 90 
      IF(NPARM.LE.0) GO TO 80
c
c** TFACT  is 95% student t-value for (NDATA-NPARM) degrees of freedom.
c [Approximate expression for (NDATA-NPARM).GT.10 accurate to ca. 0.002]
c
      TFACT= 0.D0
      IF(NDATA.GT.NPARM) THEN
          IDF= NDATA-NPARM
          IF(IDF.GT.10) TFACT= 1.960D0*DEXP(1.265D0/DFLOAT(IDF))
          IF(IDF.LE.10) TFACT= F95(IDF) 
        ELSE
          TFACT= 0.D0
        ENDIF
      DO 5 I = 1,NDATA
    5 DY(I) = 0.D0
      DO 7 I = 1,NPARM
        PS(I) = 0.D0
        PX(I) = 0.D0
        DO 8 J = 1,NPARM
    8     CM(I,J) = 0.D0
    7 CONTINUE
c
c** Begin by forming the Jacobian Matrix from the input partial 
c  derivative matrix DYDP.  For VERY large data sets, these partial 
c  derivatives may be generated inside this loop (see version LLSQFVL).
c
      DO 11 I = 1,NDATA
        S = 1.D0 / UY(I)
        U = YO(I) * S
        DO 12 J = 1,NPARM
   12     P(J) = DYDP(I,J) * S
        CALL QROD(NPARM,MXPARM,MXPARM,CM,P,DY,U,PS,PX)
        IF(I .GT. NPARM) DY(I) = U
   11 CONTINUE
c
c** Compute the inverse of  CM 
c
      CM(1,1) = 1.D0 / CM(1,1)
      DO 15 I = 2,NPARM
        L = I - 1
        DO 16 J = 1,L
          S = 0.D0
          DO 17 K = J,L
   17       S = S + CM(K,I) * CM(J,K)
   16     CM(J,I) = -S / CM(I,I)
   15   CM(I,I) = 1.D0 / CM(I,I)
c
c** Solve for parameter values  P(j)
c
      DO 18 I = 1,NPARM
        J = NPARM - I + 1
        P(J) = 0.D0
        DO 19 K = J,NPARM
   19     P(J) = P(J) + CM(J,K) * DY(K)
   18 CONTINUE
c
c** Calculate standard deviation of the fit
c
      S = 0.D0
      M = NPARM + 1
      DO 20 I = M,NDATA
   20   S = S + DY(I) * DY(I)
      IF(NDATA.GT.NPARM) THEN
          S = S / DFLOAT(NDATA - NPARM)
        else
          S = 0.d0
        endif
      SERR= DSQRT(S)
c
c** Get (upper triangular) Covarience Matrix  CM(j,k)
c  
      DO 24 I = 1,NPARM
        DO 24 J = I,NPARM
          U = 0.D0
          DO 22 K = J,NPARM
   22       U = U + CM(I,K) * CM(J,K)
   24     CM(I,J) = S * U
c
c** Generate Parameter Uncertainties  UP(j) and (symmetric) 
c                    covariance matrix
c
      DO 30 J = 1,NPARM
        UP(J) = DSQRT(CM(J,J))
        DO 26 K= J,NPARM
   26     CM(J,K)= CM(J,K)/UP(J)
        DO 28 K= 1,J
          CM(K,J)= CM(K,J)/UP(J)
   28     CM(J,K)= CM(K,J)
        UP(J)= TFACT*UP(J)
   30   PX(J)= 0.d0
c
c** Generate differences:   DY(i) = [YO(i) - Ycalc(i)]  and prepare to
c             calculate parameter truncation allowances
c
      DO 40 I = 1,NDATA
        S = 1.D0 / UY(I)
        U = 0.D0
        DO 38 J = 1,NPARM
          PX(J)= PX(J)+ (DYDP(I,J)*S)**2
   38     U = U + DYDP(I,J) * P(J)
        DY(I) = YO(I) - U
   40   CONTINUE
c
c** Calculate the 'parameter sensitivities', changes in P(j) which would 
c  change predictions of input data by an RMS average of  SERR*0.1/NPARM
c
      U= SERR*0.1d0/DFLOAT(NPARM)
      DO 44 J = 1,NPARM
   44 PS(J)= U*DSQRT(NDATA/PX(J))
c
      RETURN
c
c** If no parameters varied - simply calculate RMS deviation = SERR
c
   80 SERR= 0.D0
      DO 84 I= 1,NDATA
      DY(I)= YO(I)
   84 SERR= SERR+ (DY(I)/UY(I))**2
      SERR= DSQRT(SERR/DFLOAT(NDATA))
      RETURN
c
c** If array dimensioning inadequate, print warning & then STOP
c
   90 WRITE(6,601) NDATA,MXDATA,NPARM,MXPARM
  601 FORMAT(/' *** Dimensioning problems in LLSQF *** (NDATA, MXDATA, N
     1PARM, MXPARM)  =  (',I5,4(' ,',I5),' )')
      STOP
      END
C
C
c_______________________________________________________________________


