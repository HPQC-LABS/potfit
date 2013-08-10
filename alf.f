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
c** Numerical factor  16.857629205 (+/- 0.000,000,02) based on Compton
c  wavelength of proton & proton mass (u) from 2011 physical constants.
      BZ= ZMU/16.857629205d0
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

