c***********************************************************************
      SUBROUTINE TVSORT(ISTATE,NPTOT,VMAX,NTVALL,TVNAME)
c***********************************************************************
c** Subroutine to sort through global data file, and for each isotopomer
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
c  NTVALL  is the number of term value parameters for this state
c  TVNAME(j)  is the alphameric name identifying term value parameter j
c
c** Internally
c-------------
c  NLV(v,J.p) * initially, counts transitions for level {v,J,p} of a 
c                          given isotopologue
c           * later reset it as the parameter index for that term value
c-----------------------------------------------------------------------
      INCLUDE 'arrsizes.h'
      INCLUDE 'BLKISOT.h'
      INCLUDE 'BLKDATA.h'
c
      INTEGER I,J,P,IBAND,ISOT,ISTATE,NPTOT,LOWEST,
     1 VMAX(NSTATEMX,NISTPMX),NLV(0:NVIBMX,0:NVIBMX,-1:1),
     2 NTVS(NSTATEMX,NISTPMX),NTVALL(0:NSTATEMX)
      CHARACTER*24 TVNAME(NPARMX)
c=======================================================================
      WRITE(6,600) SLABL(ISTATE) 
      LOWEST= 1
      IF(ISTATE.GT.1) LOWEST= 0
      NTVALL(ISTATE)= 0
      DO  ISOT= 1, NISTP
c** First ... zero transition counter array for this isotopomer
          DO  I= 0, VMAX(ISTATE,ISOT)
              DO  J= 0, NVIBMX 
                  DO  P= -1,1
                      NLV(I,J,P)= 0
                      ENDDO
                  ENDDO
              ENDDO
          DO  IBAND= 1, NBANDTOT
c** Then ... search for bands involving isotopomer ISOT in this state
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
ccc                       WRITE(31,700) SLABL(ISTATE),I,J,P,ISOT,
ccc  1                                                      NLV(I,J,P)
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
          ENDDO
c
      RETURN
  600 FORMAT(/' For State ',A2,'  fit to individual term values for each
     1  {v,J,p,isot}'/1x,6('******'))
  602 FORMAT(/' *** ARRAY DIMENSION PROBLEM ***  JP(ISTATE)=',i2,
     1  ',ISOT=',I2,')=',i3,'  greater than  NVIBMX=',i4)
  604 FORMAT(/' *** ARRAY DIMENSION PROBLEM ***  JPP(ISTATE)=',i2,
     1  ',ISOT=',I2,')=',i3,'  greater than  NVIBMX=',i4)
  606 FORMAT(/'  Absolute zero of energy is fixed at level {v=',i3,
     1 ', J=',i3,', p=',i2,'}'/1x,12('**'),10x,'of isotopomer ',i2,
     2 ' of  State ',A2)
  608 FORMAT(' For ',A2,'(',i3,')-',A2,'(',I3,')  fit to',i5,
     1 ' T(v,J,p) term values,'/20x,'of which',i5,' are involved in only
     2 one transition')
  700 FORMAT("'",'T(',A2,':',i3,',',i3,',',SP,i2,';',SS,i2,')',I5,"'")
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

