c***********************************************************************
      SUBROUTINE ALF(NDP,RH,NCN,RR,V,SWF,VLIM,MAXMIN,KVMAX,KVtop,NVIBMX,
     1                         VMAXX,AFLAG,ZMU,EPS,GV,INNODE,INNR,IWR)
c***********************************************************************
c-----------------------------------------------------------------------
c** The subroutine ALF (Automatic vibrational Level Finder) will
c   automatically generate the eigenvalues from vibrational level v=0
c   to a user-specified level (v=KVMAX) or the highest vibrational
c   level it can find in a given smooth single (or double)
c   minimum potential (V). These energies are stored and returned to the
c   calling program in the molecular constants array GV(v=0-KVMAX).
c** For any errors that cannot be resolved within the subroutine, ALF
c   returns AFLAG with a value that defines which error had occured.
c++++++++++   Version last updated  by Shirin March 24, 2017 +++++++++++
c+++++++ {tinkered with ICOR cases aimed at double well PEFs +++++++++++
c+++++++++++++   COPYRIGHT 2008-17  by  Robert J. Le Roy   +++++++++++++
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
c         .NE. 0 : print error & warning descriptions.
c           >= 1 : also print final eigenvalues & node count.
c           >= 2 : also show end-of-range wave function amplitudes.
c           >= 3 : print also intermediate trial eigenvalues, etc.
c     SWF not needed for in/out {I think} but allows use of same long array
c
c** On exit:
c    KVtop   is vib.quantum number for the highest vibrational level
c            found (may be less than the input desired value of KVMAX).
c    VMAXX   is MAX{energy at barrier maximim,asymptote}
c    AFLAG   returns calculation outcome to calling program.
c            >=  0 : found all levels to v=KVMAX{input} & AFLAG= J 
c             = -1 : highest level found  KVtop < KVMAX
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
      INTEGER NDP,KVMAX,KVtop,KV,KVB,KVBB,AFLAG,NF,NBEG,NEND,NVIBMX,
     1  NBEGG(0:NVIBMX),NENDD(0:NVIBMX),INNR(0:NVIBMX),ICOR,IWR,INNSAV,
     2  IPMIN(10),IPMINN,I,LTRY,AWO,INNODE,INNER,LPRWF,JROT,NCN,NPMIN,
     3  NPMAX,MAXMIN
c
      REAL*8 RMIN,RH,RBAR,RR(NDP),V(NDP),SWF(NDP),VLIM,EO,EX,ZMU,EPS,
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
      INNER= 0
      KV= 0
      KVB= -1
      KVtop= KVMAX
c** Initialize level labels and energies for each well
      DO  I= 0,KVMAX
          INNR(I)= -2
          GV(I)= 0.d0
          ENDDO
c** Store input rotational quantum number.
      JROT= AFLAG
      AFLAG= -1
c
c** Numerical factor  16.857629206 (+/- 0.000,000,013) based on Compton
c  wavelength of proton & proton mass (u) from 2011 physical constants.
      BZ= ZMU/16.857629206d0
      BFCT= BZ*RH*RH
!!!   IWR=5             !!! for testing
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
              IF(NPMAX.EQ.10) GOTO 20               !! array bound stop
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
      IF((NPMIN.GT.1).AND.(IWR.NE.0)) THEN
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
      LTRY= 0
      KVB= KV
      KV= NF
  110 ESAV= EO
      INNSAV= INNER
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
          KVtop= NF-1
          EO= -9.9d9
          AFLAG= -1
cc        IF(ICOR.EQ.0) GOTO 115
          IF(JROT.EQ.0) STOP
          GOTO 115 
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
          IF(IWR.GT.1) write(6,777) RBAR,RPMAX(1),INNER
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
              KVtop= KV             !! for case when vD-v < 1 for v=KV
              GOTO 200
              ENDIF
c*** For double well potential if the new energy calculated in within 1
c*** wavenumber of the previous one and the energy for the previous 2
c** levels is also within 1 wavenumber of each other, increase this new
c** energy calculated by DGDV2 (the outer level energy)

         IF (NPMIN.EQ.2) THEN
          IF (((EO- GV(NF-1)).LE.1).AND.(GV(NF-1)- GV(NF-2)).LE.1) THEN
               EO= EO+ DGDV2
            ELSE
               EO= EO
          ENDIF
         ENDIF
!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!        PWCN= 2.d0*NCN/FLOAT(NCN-2)
!!        PWCNinv= 1.d0/PWCN
!!        IF(KV.GE.1) THEN   !! alternate predictor for next levl: see LEVEL 
!!            VDMV= ((VLIM-GV(KV))/(VLIM- !GV(KV-1)))**PWCNinv      !!  diary
!!            VDMV= VDMV/(1.d0 - VDMV)
!!            EO= VLIM - (VLIM-GV(KV))*((VDMV-
!1.d0)/VDMV)**PWCN
!!            IF(IWR.GE.2) WRITE(6,6122) VDMV,EO
!6122 FORMAT(" Applying NDT to past 2 Eb's predicts
!(vD-v)=",f9.4,
!!   1  '   E(next)=',1PD13.5)
!!            ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          IF(EO.GT.VPMAX(NPMAX)) THEN
c... if estimated energy above highest barrier, set value slightly below it
              EO=  VPMAX(NPMAX) - 0.10d0*DGDV2
              ICOR= ICOR+10
            ELSE
              IF(DGDV2.LT.0.d0) THEN
c... SCECOR returned negative phase integral, so quit loop & RETURN
                  WRITE(6,628) JROT,EO
                  AFLAG= -1
                  GOTO 110
                  ENDIF
            ENDIF
          LTRY= 0
          GOTO 100
          ENDIF
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  115    IF(KV.NE.NF) THEN
c*** If last level found was not the desired one ...
         IF(KV.GE.0) THEN  
             IF(INNR(KV).LT.-1) THEN
c... Record vib. level actually found (if haven't already) for posterity.
                  GV(KV)= EO
                  INNR(KV)= INNER
                  ENDIF
              ENDIF
              IF(KV.EQ.NF-1) THEN         !! if re-found previous level
cc               EO= EO- DGDV2   
                  EO= ESAV + (ESAV-EO)  !!Previous energy - new energy (lower)
                  KV= NF      !! try doubling interval from initial ESAV
                  GOTO 110
                  ENDIF
          ICOR= ICOR+1
          IF(ICOR.EQ.1) THEN
c!! on first failure,switch the value of INNER and try again w. same EO
              INNER= INNSAV
              IF(INNER.GT.0) THEN
                  INNER= 0
                ELSE
                  INNER= 1
                ENDIF  
              EO= ESAV
              KV= NF
              GOTO 110
              ENDIF
          WRITE(6,*)  "   ICOR=",ICOR,"   v=",KV,"   J=",JROT
          IF((ICOR.EQ.2).AND.(NPMIN.EQ.1)) THEN
c!! after second failure for single minimum PEF, assume looking for last
c    Q'bdd level, so check for proximity to barrier max and if so,
c   set value slightly below it
              EO=  VPMAX(1) - 0.10d0*DGDV2
              ICOR= ICOR+ 10
              GO TO 110
              ENDIF
          IF((ICOR.GE.2).AND.(NPMIN.EQ.2).AND.(ICOR.LE.100)) THEN 
c!! after second failure for double minimum PEF step down from bad uppr
c   level in steps of 0.05*spacing
              IF (ICOR.EQ.2) EX= EO     !! Remember the original trial energy
                IF (GV(NF-1).GE.EX) THEN
                  IF((ICOR/2*2).EQ.ICOR) THEN    !! ICOR even
                  EO= GV(NF-1)+ (GV(NF-1)- EX)*0.02*(ICOR/2)
                  KV= NF
                  INNER=0
                  GO TO 110
                ELSE                          !! ICOR odd
                  INNER= INNSAV
                  IF(INNER.GT.0) INNER= 0
                  IF(INNER.EQ.0) INNER= 1
                  EO= ESAV
                  KV= NF
                  GOTO 110
                  ENDIF
              ELSE        
                 IF((ICOR/2*2).EQ.ICOR) THEN    !! ICOR even
                  EO= EX- (EX - GV(NF-1))*0.02*(ICOR/2)
                  KV= NF
                  INNER=0
                  GO TO 110
                ELSE                          !! ICOR odd
                  INNER= INNSAV
                  IF(INNER.GT.0) INNER= 0
                  IF(INNER.EQ.0) INNER= 1
                  EO= ESAV
                  KV= NF
                  GOTO 110
                  ENDIF
                ENDIF
                ENDIF
         IF(ICOR.EQ.100) THEN     !!Small intervals so go until Icor 95
c... Call subroutine using semiclassical methods to estimate correct energy
c...        CALL SCECOR(KV,NF,JROT,INNER,ICOR,IWR,EO,RH,BFCT,NDP,NCN,
c...    1                                             V,VMAXX,VLIM,DGDV2)
c...             IF(EO.GT.VPMAX(NPMAX)) THEN
c... if estimated energy above highest barrier, set value below it
c...               KV= 999
c...                EO=  VPMAX(NPMAX) - 0.05d0*DGDV2
c...                 ENDIF
c...             GOTO 100
              STOP
              ENDIF
c** If the calculated wavefunction is still for the wrong vibrational
c   level, then write out a warning return
            WRITE(6,630) NF,JROT
            KVtop= NF-1                   !! from  KVMAX= NF-1  17/02/17
            ENDIF
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  200 IF(AFLAG.LT.0) THEN
c** If unable to find all KVMAX+1 levels requested, then return KVtop as
c  v for the highest vibrational level actually found, and print out the
c  the energy of that level.
          IF(KV.GE.0) KVtop= KV
          IF(AWO.NE.0) WRITE(6,632) KVtop, GV(KVtop)
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

