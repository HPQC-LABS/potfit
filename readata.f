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
      INCLUDE 'arrsizes.h'
      INCLUDE 'BLKPOT.h'
      INCLUDE 'BLKISOT.h'
      INCLUDE 'BLKDATA.h'
      INCLUDE 'BLKPARAM.h'
      INCLUDE 'BLKTYPE.h'
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
!!        IF(ISOT.GT.1) THEN       !! not needed - done in calling program
!!            VMAXespp= INT((VMAX(ESPP,ISOT)+0.5d0)/RSQMU(ISOT)-0.5d0) !! added
!!            VMINespp= INT((VMIN(ESPP,ISOT)+0.5d0)/RSQMU(ISOT)-0.5d0) !! added
!!            JTRUNCespp= INT(JTRUNC(ESPP)/RSQMU(ISOT))
!!            ENDIF
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
c** Data are Acoustic Virial Coefficients for electronic state IEPP= ESPP
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
c** Book-keeping for Acoustic Virial data
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
     1 '(',I3,') Acoustic Virial coefficients included in data set' )
  650 FORMAT(/' Data input IGNORES',i4,' fluorescence series consisting'
     1 ,' of only  onee  line!')
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

