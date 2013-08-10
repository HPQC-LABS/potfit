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

