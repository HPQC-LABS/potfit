C***********************************************************************
      SUBROUTINE GPROUND(IROUND,NPTOT,NPMAX,NPAR1,NPAR2,LPRINT,IFXP,
     1                                                          PV,PU)
c** Subroutine to round off parameters PV(i), i= NPAR1 to NPAR2, at the
c  |IROUND|'th significant digit of the smallest of their uncertainties
c  min{U(i)}.  This procedure does NOT attempt to correct the remaining
c  parameters to compensate for these changes (as ROUND does), so this
c  procedure is not appropriate for nonlinear parameters.
c** On return, the rounded values replaces the initial values of  PV(i).
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c                COPYRIGHT 2000-2004  by  Robert J. Le Roy             +
c   Dept. of Chemistry, Univ. of Waterloo, Waterloo, Ontario, Canada   +
c    This software may not be sold or any other commercial use made    +
c      of it without the express written permission of the author.     +
c                      Version of 27 January 2004                      +
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      INTEGER  IROUND,NPMAX,NPTOT,NPAR1,NPAR2,NPARM,IRND,KRND,LPRINT
      INTEGER  IFXP(NPTOT)
      REAL*8  PV(NPMAX),PU(NPMAX),CNST,CRND,XRND,FCT,YY,UNC
c
c** Loop over & round off the parameters # NPAR1 to NPAR2 
      IF(LPRINT.GE.2) WRITE(6,602)  NPAR2-NPAR1+1,NPTOT,NPAR1,NPAR2
      UNC= 99.d99
      DO  NPARM= NPAR1, NPAR2
          IF(PU(NPARM).LT.UNC) UNC= PU(NPARM)
          ENDDO
      DO  NPARM= NPAR1, NPAR2
c** First ... fiddle with log's to perform the rounding
          XRND= DLOG10(UNC)
          IRND= INT(XRND)
          IF(XRND.GT.0) IRND=IRND+1
          IRND= IRND- IROUND
          FCT= 10.D0**IRND
          CNST= PV(NPARM)
          YY= CNST
          CRND= PV(NPARM)/FCT
          XRND= 0.d0
c ... if rounding goes past REAL*8 precision, retain unrounded constant
          IF(DABS(CRND).GE.1.D+16) THEN
              WRITE(6,600) IROUND,NPARM
               RETURN
               ENDIF
          IF(DABS(CRND).GE.1.D+8) THEN
c ... to avoid problems from overflow of I*4 integers ...
              KRND= NINT(CRND/1.D+8)
              XRND= KRND*1.D+8
              CRND= CRND-XRND
              XRND= XRND*FCT
              END IF
          IRND= NINT(CRND)
          CNST= IRND*FCT+ XRND
          PV(NPARM) = CNST
          IFXP(NPARM)= 1
          IF(LPRINT.GE.2) WRITE(6,604) NPARM,YY,PV(NPARM)
  604 FORMAT(5x,'Round parameter #',i4,' from',G20.12,'  to',G20.12)
          ENDDO
      NPARM= NPARM- 1
      RETURN
  600 FORMAT(' =',39('==')/' Caution:',i3,'-digit rounding of parameter-
     1',i2,' would exceed (assumed) REAL*8'/' ********   precision overf
     2low at 1.D+16, so keep unrounded constant')
  602 FORMAT(' Rounding off ',i5,' of the ',i5,' parameters #:',i5,
     1 ' to',i5)
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

