c**********************************************************************
      SUBROUTINE DIFFSTATS(NSTATES,NPTOT,ROBUST,MKPRED)
c** Subroutine to summarise dimensionless standard errors on a band-by-
c  band basis, and (if desired) print [obs.-calc.] values to channel-8.
c-----------------------------------------------------------------------
c                 Version of  30 January  2013
c     last change - added NPTOT to printout
c-----------------------------------------------------------------------
      INCLUDE 'arrsizes.h'
      INCLUDE 'BLKISOT.h'
      INCLUDE 'BLKDATA.h'
      INCLUDE 'BLKTYPE.h'
c
      INTEGER I,IBB,ISOT,ISTATE,ISTATEE,J,K,NSTATES,MKPRED,ROBUST,NPTOT
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
c** Summarize data discrepancies for one isotopomer at a time.
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
              IBB= YPR(ISOT,ISTATE,4,4,I)
              IF(MKPRED.LE.0) THEN
                  CALL BNDERR(IFIRST(IBB),ILAST(IBB),ROBUST,AVE,RMSR,
     1                                             SSQTOT,DFREQ,UFREQ)
                  RMSTOT= RMSTOT+ YPR(ISOT,ISTATE,4,3,I)*RMSR**2
                  AVETOT= AVETOT+ YPR(ISOT,ISTATE,4,3,I)*AVE
                  WRITE(6,606)YPR(ISOT,ISTATE,4,2,I),
     1                  YPR(ISOT,ISTATE,4,1,I),YPR(ISOT,ISTATE,4,3,I),
     2                  YPR(ISOT,ISTATE,4,5,I),YPR(ISOT,ISTATE,4,6,I),
     3                  AVEUFREQ(IBB),MAXUFREQ(IBB),AVE,RMSR
                  ENDIF
              WRITE(8,605)
              IF(MKPRED.LE.0) WRITE(8,606) YPR(ISOT,ISTATE,4,2,I),
     1                  YPR(ISOT,ISTATE,4,1,I),YPR(ISOT,ISTATE,4,3,I),
     2                  YPR(ISOT,ISTATE,4,5,I),YPR(ISOT,ISTATE,4,6,I),
     4                  AVEUFREQ(IBB),MAXUFREQ(IBB),AVE,RMSR
              IF(MKPRED.GT.0) THEN
                  WRITE(8,606) YPR(ISOT,ISTATE,4,2,I),
     1                  YPR(ISOT,ISTATE,4,1,I),YPR(ISOT,ISTATE,4,3,I),
     2                  YPR(ISOT,ISTATE,4,5,I),YPR(ISOT,ISTATE,4,6,I)
                  WRITE(4,640) YPR(ISOT,ISTATE,4,2,I),
     1             YPR(ISOT,ISTATE,4,1,I),SLABL(ISTATE),SLABL(ISTATE),
     2                                              (MN(K,ISOT),K=1,2)
                  ENDIF
  640 FORMAT(/2I4,2(2x,"'",A2,"'"),2x,2I4)
              CALL PBNDERR(IBB,MKPRED,NEF)
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
              IBB= YPR(ISOT,ISTATE,3,4,I)
              IF(MKPRED.LE.0) THEN
                  CALL BNDERR(IFIRST(IBB),ILAST(IBB),ROBUST,AVE,RMSR,
     1                                             SSQTOT,DFREQ,UFREQ)
                  RMSTOT= RMSTOT+ YPR(ISOT,ISTATE,3,3,I)*RMSR**2
                  AVETOT= AVETOT+ YPR(ISOT,ISTATE,3,3,I)*AVE
                  WRITE(6,606) YPR(ISOT,ISTATE,3,2,I),
     1                  YPR(ISOT,ISTATE,3,1,I),YPR(ISOT,ISTATE,3,3,I),
     2                  YPR(ISOT,ISTATE,3,5,I),YPR(ISOT,ISTATE,3,6,I),
     3                  AVEUFREQ(IBB),MAXUFREQ(IBB),AVE,RMSR
                  ENDIF
              WRITE(8,605)
              IF(MKPRED.LE.0) WRITE(8,606) YPR(ISOT,ISTATE,3,2,I),
     1                  YPR(ISOT,ISTATE,3,1,I),YPR(ISOT,ISTATE,3,3,I),
     2                  YPR(ISOT,ISTATE,3,5,I),YPR(ISOT,ISTATE,3,6,I),
     3                  AVEUFREQ(IBB),MAXUFREQ(IBB),AVE,RMSR
              IF(MKPRED.GT.0) THEN
                  WRITE(8,606) YPR(ISOT,ISTATE,3,2,I),
     1                  YPR(ISOT,ISTATE,3,1,I),YPR(ISOT,ISTATE,3,3,I),
     2                  YPR(ISOT,ISTATE,3,5,I),YPR(ISOT,ISTATE,3,6,I)
                  WRITE(4,640) YPR(ISOT,ISTATE,3,2,I),
     1             YPR(ISOT,ISTATE,3,1,I),SLABL(ISTATE),SLABL(ISTATE),
     2                                              (MN(K,ISOT),K=1,2)
                  ENDIF
              CALL PBNDERR(IBB,MKPRED,NEF)
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
                  IBB= YPR(ISOT,ISTATE,2,4,I)
                  IF(IEP(IBB).EQ.ISTATEE) THEN
                      IF(MKPRED.LE.0) THEN
                          CALL BNDERR(IFIRST(IBB),ILAST(IBB),ROBUST,AVE,
     1                                        RMSR,SSQTOT,DFREQ,UFREQ)
                          RMSTOT= RMSTOT+ YPR(ISOT,ISTATE,2,3,I)*RMSR**2
                          AVETOT= AVETOT+ YPR(ISOT,ISTATE,2,3,I)*AVE
                          WRITE(6,606) YPR(ISOT,ISTATE,2,2,I),
     1                  YPR(ISOT,ISTATE,2,1,I),YPR(ISOT,ISTATE,2,3,I),
     2                  YPR(ISOT,ISTATE,2,5,I),YPR(ISOT,ISTATE,2,6,I),
     3                            AVEUFREQ(IBB),MAXUFREQ(IBB),AVE,RMSR
                          ENDIF
                      WRITE(8,605)
                      IF(MKPRED.LE.0) WRITE(8,606) 
     1                  YPR(ISOT,ISTATE,2,2,I),YPR(ISOT,ISTATE,2,1,I),
     2                  YPR(ISOT,ISTATE,2,3,I),YPR(ISOT,ISTATE,2,5,I),
     3             YPR(ISOT,ISTATE,2,6,I),AVEUFREQ(IBB),MAXUFREQ(IBB),
     4                                                        AVE,RMSR
                      IF(MKPRED.GT.0) THEN
                           WRITE(8,606) YPR(ISOT,ISTATE,2,2,I),
     1                  YPR(ISOT,ISTATE,2,1,I),YPR(ISOT,ISTATE,2,3,I),
     2                  YPR(ISOT,ISTATE,2,5,I),YPR(ISOT,ISTATE,2,6,I)
                           WRITE(4,640) YPR(ISOT,ISTATE,2,2,I),
     1             YPR(ISOT,ISTATE,2,1,I),SLABL(ISTATE),SLABL(ISTATE),
     2                                              (MN(K,ISOT),K=1,2)
                  ENDIF
                      CALL PBNDERR(IBB,MKPRED,NEF)
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
              IBB= YPR(ISOT,ISTATE,1,4,I)
              CALL BNDERR(IFIRST(IBB),ILAST(IBB),ROBUST,AVE,RMSR,
     1                                             SSQTOT,DFREQ,UFREQ)
              RMSTOT= RMSTOT+ YPR(ISOT,ISTATE,1,3,I)*RMSR**2
              AVETOT= AVETOT+ YPR(ISOT,ISTATE,1,3,I)*AVE
              WRITE(6,614) YPR(ISOT,ISTATE,1,1,I),
     1                   YPR(ISOT,ISTATE,1,2,I),NEF(EFP(IFIRST(IBB))),
     2                  YPR(ISOT,ISTATE,1,3,I),YPR(ISOT,ISTATE,1,5,I),
     3                  YPR(ISOT,ISTATE,1,6,I),
     4                  AVEUFREQ(IBB),MAXUFREQ(IBB),AVE,RMSR
              WRITE(8,617)
              WRITE(8,614) YPR(ISOT,ISTATE,1,1,I),
     1                  YPR(ISOT,ISTATE,1,2,I),NEF(EFP(IFIRST(IBB))),
     2                  YPR(ISOT,ISTATE,1,3,I),YPR(ISOT,ISTATE,1,5,I),
     3                  YPR(ISOT,ISTATE,1,6,I),
     4                  AVEUFREQ(IBB),MAXUFREQ(IBB),AVE,RMSR
              CALL PBNDERR(IBB,MKPRED,NEF)
              ENDDO
              RMSTOT= DSQRT(RMSTOT/NTRANSFS(ISOT,ISTATE))
              AVETOT= AVETOT/NTRANSFS(ISOT,ISTATE)
              WRITE(6,632) NTRANSFS(ISOT,ISTATE),AVETOT,RMSTOT
          ENDIF
c
      IF(NEBPAS(ISOT,ISTATE).GT.0) THEN
c** Book-keeping for  PAS  data
          IBB= YPR(ISOT,ISTATE,7,4,1)
          CALL BNDERR(IFIRST(IBB),ILAST(IBB),ROBUST,AVE,RMSR,SSQTOT,
     1                                                    DFREQ,UFREQ)
          WRITE(6,626) NEBPAS(ISOT,ISTATE),SLABL(ISTATE),(NAME(I),
     1 MN(I,ISOT),I=1,2),YPR(ISOT,ISTATE,7,3,1),YPR(ISOT,ISTATE,7,5,1),
     2      YPR(ISOT,ISTATE,7,6,1),AVEUFREQ(IBB),MAXUFREQ(IBB),AVE,RMSR
          WRITE(8,626) NEBPAS(ISOT,ISTATE),SLABL(ISTATE),(NAME(I),
     1 MN(I,ISOT),I=1,2),YPR(ISOT,ISTATE,7,3,1),YPR(ISOT,ISTATE,7,5,1),
     2      YPR(ISOT,ISTATE,7,6,1),AVEUFREQ(IBB),MAXUFREQ(IBB),AVE,RMSR
          WRITE(8,627)
          DO  I= IFIRST(IBB),ILAST(IBB)
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
      IF(NWIDTH(ISOT,ISTATE).GT.0) THEN
c** Book-keeping for  Tunneling Width  data
          IBB= YPR(ISOT,ISTATE,6,4,1)    
          CALL BNDERR(IFIRST(IBB),ILAST(IBB),ROBUST,AVE,RMSR,SSQTOT,
     1                                                    DFREQ,UFREQ)
          WRITE(6,620) NWIDTH(ISOT,ISTATE),SLABL(ISTATE),
     1              (NAME(I),MN(I,ISOT),I=1,2),YPR(ISOT,ISTATE,6,3,1),
     2              YPR(ISOT,ISTATE,6,5,1),YPR(ISOT,ISTATE,6,6,1),
     3              AVEUFREQ(IBB),MAXUFREQ(IBB),AVE,RMSR
          WRITE(8,620) NWIDTH(ISOT,ISTATE),SLABL(ISTATE),
     1              (NAME(I),MN(I,ISOT),I=1,2),YPR(ISOT,ISTATE,6,3,1),
     2              YPR(ISOT,ISTATE,6,5,1),YPR(ISOT,ISTATE,6,6,1),
     3              AVEUFREQ(IBB),MAXUFREQ(IBB),AVE,RMSR
          DO  J= IFIRST(IBB),ILAST(IBB)
              WRITE(6,622) JP(J),JPP(J),NEF(EFPP(J)),FREQ(J),UFREQ(J),
     1                                      DFREQ(J),DFREQ(J)/UFREQ(J)
              WRITE(8,622) JP(J),JPP(J),NEF(EFPP(J)),FREQ(J),
     1                             UFREQ(J),DFREQ(J),DFREQ(J)/UFREQ(J)
              ENDDO
          ENDIF
c
      IF(NVVPP(ISOT,ISTATE).GT.0) THEN
c** Book-keeping for potential function values as data .....
          IBB= YPR(ISOT,ISTATE,5,4,1)
          CALL BNDERR(IFIRST(IBB),ILAST(IBB),ROBUST,AVE,RMSR,SSQTOT,
     1                                                    DFREQ,UFREQ)
          WRITE(6,638) NVVPP(ISOT,ISTATE),SLABL(ISTATE),
     1              (NAME(I),MN(I,ISOT),I=1,2),YPR(ISOT,ISTATE,5,3,1),
     2              AVEUFREQ(IBB),MAXUFREQ(IBB),AVE,RMSR
          WRITE(8,638) NVVPP(ISOT,ISTATE),SLABL(ISTATE),
     1              (NAME(I),MN(I,ISOT),I=1,2),YPR(ISOT,ISTATE,5,3,1),
     2              AVEUFREQ(IBB),MAXUFREQ(IBB),AVE,RMSR
          DO  J= IFIRST(IBB),ILAST(IBB)
              WRITE(6,637) TEMP(J),FREQ(J),UFREQ(J),DFREQ(J),
     1                                               DFREQ(J)/UFREQ(J)
              WRITE(8,637) TEMP(J),FREQ(J),UFREQ(J),DFREQ(J),
     1                                               DFREQ(J)/UFREQ(J)
              ENDDO
          ENDIF
c
      IF(NVIRIAL(ISOT,ISTATE).GT.0) THEN
c** Book-keeping for  Virial Coefficient  data
          IBB= YPR(ISOT,ISTATE,8,4,1)
          CALL BNDERR(IFIRST(IBB),ILAST(IBB),ROBUST,AVE,RMSR,SSQTOT,
     1                                                    DFREQ,UFREQ)
          WRITE(6,634) NVIRIAL(ISOT,ISTATE),SLABL(ISTATE),
     1              (NAME(I),MN(I,ISOT),I=1,2),YPR(ISOT,ISTATE,8,3,1),
     2              AVEUFREQ(IBB),MAXUFREQ(IBB),AVE,RMSR
          WRITE(8,634) NVIRIAL(ISOT,ISTATE),SLABL(ISTATE),
     1              (NAME(I),MN(I,ISOT),I=1,2),YPR(ISOT,ISTATE,8,3,1),
     2              AVEUFREQ(IBB),MAXUFREQ(IBB),AVE,RMSR
          DO  J= IFIRST(IBB),ILAST(IBB)
              WRITE(6,636) TEMP(J),FREQ(J),UFREQ(J),DFREQ(J),
     1                                               DFREQ(J)/UFREQ(J)
              WRITE(8,636) TEMP(J),FREQ(J),UFREQ(J),DFREQ(J),
     1                                               DFREQ(J)/UFREQ(J)
c=======================================================
c   7 State X0 Ar( 40)-Xe(132)  Virial Coefficients
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
              ENDDO
  634 FORMAT(/1x,55('=')/I5,'  State ',A2,2x,A2,'(',I3,')-',A2,'(',
     1 I3,')  Virial coefficients'/1x,10('===='),' Avge. ',('====')/5x,
     2 '#data    Av.Unc.    Max.Unc.    Err/Unc   DRMSD'/1x,10('-----')
     3 /I9,1PD12.2,D12.2,0PF11.5,F8.3/1x,23('=='),' calc-obs'/
     4  5x,'temp.',4x,'Bvir(obs)',4x,'u(Bvir)   calc-obs   /u(Bvir)'/
     5  1x,53('-'))
  636 FORMAT(F11.3,2F10.2,2F11.3)
  637 FORMAT(F11.6,F12.2,F10.2,2F11.3)
          ENDIF
  638 FORMAT(/1x,55('=')/I5,'  State ',A2,2x,A2,'(',I3,')-',A2,'(',
     1 I3,')  Potential fx. values'/1x,21('=='),' Avge. ',('====')/5x,
     2 '#data    Av.Unc.    Max.Unc.    Err/Unc   DRMSD'/1x,26('--')
     3 /I9,1PD12.2,D12.2,0PF11.5,F8.3/1x,24('=='),' calc-obs'/
     4  7x,'R',7x,'V(r)',8x,'u(V(r))   calc-obs   /u(V(r))'/
     5  1x,53('-'))
c** End of loop over the various (lower) electronic states
   90 CONTINUE
c=======================================================================
      IF(ISOT.LT.NISTP) THEN
c** If NISTP > 1, return to print data summaries for other isotopomers
          ISOT= ISOT+1
          GO TO 10
          ENDIF 
      RMSR= DSQRT(SSQTOT/COUNTOT)
      WRITE(6,624) NPTOT,COUNTOT,RMSR
      RETURN
  600 FORMAT(/1x,36('**')/'  Write to Channel-8 Predictions From Complet
     1e Set of Input Parameters!'/1x,36('**'))
  601 FORMAT(/1x,25('**')/'  Predictions From Complete Set of Input Para
     1meters!'/1x,25('**'))
  602 FORMAT(/1x,21('===')/'  *** Discrepancies for',I5,' bands/series o
     1f ',A2,'(',I3,')-',A2,'(',I3,') ***'/1x,21('==='))
  604 FORMAT(/1x,21('===')/I5,' State ',A2,1x,A2,'(',I3,')-',A2,'(',I3,
     1 ') MW transitions in',i4,' vib. levels')
  605 FORMAT(1x,16('==='),'== Avge. ========'/"   v' ",
     2  ' v" #data  J"min  J"max  Av.Unc.  Max.Unc.   Err/Unc   DRMSD'/
     1  1x,13('-----'))
  606 FORMAT(2I4,I6,3x,I4,3x,I4,1x,1P2D9.1,0PF11.5,F8.3)
  608 FORMAT(/1x,63('=')/I5,' State ',A2,1x,A2,'(',I3,')-',A2,'(',I3,
     1 ') InfraRed transitions in',I4,' bands')
  610 FORMAT(/1x,35('==')/I6,1x,A2,'(',I3,')-',A2,'(',i3,')  {State ',
     1  A2,'}--{State ',A2,'} Transitions in',i4,' bands')
  612 FORMAT(/1x,75('=')/I5,' Fluorescence transitions into State ',A2,
     1 2x,A2,'(',I3,')-',A2,'(',I3,') in',i5,' series')
  617 FORMAT(1x,52('='),'= Avge. ',15('=')/"   v'  j' p' ",
     2  '#data  v"min  v"max','  AvgeUnc  Max.Unc.   Err/Unc   DRMSD'/
     3  1x,25('---'))
  614 FORMAT(2I4,A3,I6,2I7,1x,1P2D9.1,0PF11.5,F8.3)
  616 FORMAT(/1x,66('=')/1x,I3,' State ',A2,1x,A2,'(',I3,')-',A2,'(',
     1 I3,') potential fx. values treated as independent data'/
     2 1x,20('=='),'  Avge.  ',17('=')/' #data   v"min  v"max  AvgeUnc',
     1 '  Max.Unc.  Err/Unc  DRMSD'/1x,55('-')/I5,2I7,2x,1P2D9.1,0PF9.3,
     4 F8.3/1x,30('==')/'    v  p',8x,'Bv',7x,'u(Bv)',4x,
     5 '[calc-obs]  [calc-obs]/unc',/1x,30('--'))
  618 FORMAT(I5,A3,2x,F12.8,1PD9.1,0PF13.8,F12.4)
  620 FORMAT(/1x,73('=')/1x,I3,' State ',A2,1x,A2,'(',I3,')-',A2,'(',
     1 I3,') Tunneling Widths treated as independent data'/1x,20('=='),
     2 '  Avge.  ',24('=')/' #data   v"min  v"max  AvgeUnc  Max.Unc.  Er
     3r/Unc  DRMSD'/1x,55('-')/I5,2I7,2x,1P2D9.1,0PF9.3,F8.3/
     4 1x,59('=')/'   v   J  p     Width',7x,'u(Width)  [calc-obs]  [cal
     5c-obs]/unc'/1x,59('-'))
  622 FORMAT(2I4,A3,1PD14.6,D10.1,D13.2,0PF10.3)
  624 FORMAT(/1x,33('==')/' Fit of ',I5,' param to',i6,' data:  DRMS(dev
     1n)=',G13.6/1x,33('==')) 
  626 FORMAT(/1x,29('==')/I5,' PAS Binding Energies for State ',A2,2x,
     1 A2,'(',I3,')-',A2,'(',I3,')'/1x,50('='),' Avge. ',('=')/
     2 ' #data  v_min  v_max   AvgeUnc  Max.Unc.  Err/Unc  DRMSD'/
     3  1x,29('--')/I5,2I7,2x,1P2D9.1,0PF9.3,F8.3)
  627 FORMAT(1x,50('='),'  calc-obs'/'   v   j  p       PAS(Eb)        u
     1(Eb)       calc-obs   /u(FREQ)'/1x,30('--'))
  628 FORMAT(2I4,A3,F15.8,2F11.8,F11.4,1X,A3)
  629 FORMAT(1x,30('=='))
  630 FORMAT(1x,7('--'),' For these',i6,' lines, overall:',F11.5,F8.3)
  632 FORMAT(1x,17('-'),' For these',i6,' lines, overall:',F11.5,F8.3)
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE BNDERR(FIRST,LAST,ROBUST,AVEDD,RMSDD,SSQTOT,DFREQ,
     1                                                          UFREQ)
c** Calculate the average (AVEDD) & the root mean square dimensionless 
c  deviation (RSMDD) for the band running from datum # FIRST to LAST.
      INCLUDE 'arrsizes.h'
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
      SUBROUTINE PBNDERR(IBB,MKPRED,NEF)
c** Print to channel-8 a listing of the [obs.-calc.] values for the band
c  running from datum # FIRST to LAST.                           
      INCLUDE 'arrsizes.h'             
      INCLUDE 'BLKDATA.h'
      REAL*8 DIV
      INTEGER IBB,I,MKPRED
      CHARACTER*3 marker, NEF(-1:1)
c----------------------------------------------------------------------- 
      IF(MKPRED.LE.0) WRITE(8,600)
      IF(MKPRED.GT.0) WRITE(8,601)
      DO  I= IFIRST(IBB),ILAST(IBB)
          IF(MKPRED.LE.0) THEN
              DIV= DABS(DFREQ(I)/UFREQ(I))
              marker='   '
              IF( (DIV.GE.2.d0).AND.(DIV.LT.4.d0) ) marker='*  '
              IF( (DIV.GE.4.d0).AND.(DIV.LT.8.d0) ) marker='** '
              IF( (DIV.GE.8.d0) ) marker='***'
              IF(IEP(IBB).GT.0) WRITE(8,602) VP(IBB),JP(I),NEF(EFP(I)),
     1         VPP(IBB),JPP(I),NEF(EFPP(I)),FREQ(I),UFREQ(I),DFREQ(I),
     2                                      DFREQ(I)/UFREQ(I),marker
              IF(IEP(IBB).EQ.0) WRITE(8,602) VP(IBB),VPP(IBB),
     1                  NEF(EFP(I)),JP(I),JPP(I),NEF(EFPP(I)),FREQ(I),
     2                      UFREQ(I),DFREQ(I),DFREQ(I)/UFREQ(I),marker
            ELSE
              WRITE(8,602) VP(IBB),JP(I),NEF(EFP(I)),VPP(IBB),JPP(I),
     1                                           NEF(EFPP(I)),DFREQ(I)
              WRITE(4,608) JP(I),EFP(I),JPP(I),EFPP(I),DFREQ(I),UFREQ(I)
c* Print predictions in alternate (Lyon) format
c             WRITE(11,606)VP(IBB),VPP(IBB),JPP(I),JP(I)-JPP(I),DFREQ(I)
c 606 FORMAT(2I4,I5,I4,f13.4)
            ENDIF
          ENDDO
      WRITE(8,604)
      RETURN                       
  600 FORMAT(1x,59('='),'  calc-obs'/ "   v'  J' p'",
     1  '  v"  J" p"    FREQ(obs)     u(FREQ)    calc-obs  /u(FREQ)'/
     2  1x,69('-'))
  601 FORMAT(1x,36('=')/ "   v'  J' p'",'  v"  J" p"   FREQ(calc)'/
     1  1x,36('-'))
  602 FORMAT(2(2I4,A3),f14.6,2f12.6,f10.4,1x,A3)
  604 FORMAT(1x,69('-'))
  608 FORMAT(I5,I3,I5,I3,F13.4,F9.4)
      END   
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

