c***********************************************************************
      SUBROUTINE FUNUNC(ISTATE,WRITFILE,PU,CM)
c***********************************************************************
c** This subroutine will calculate the uncertainties in the radial
c  strength functions for V(r), phi(r), UA(r), UB(r), tA(r), tB(r) and 
c  wRAD(r) and print them out in `Tecplot' format.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** On entry:
c ... OSEL     print PEF values at every |OSEL|'th mesh point
c           and if OSEL < 0, also the BOB fx an every |OSEL|'th mesh point
c ... ISTATE   electronic state counter
c ... PU(n)    uncertainties in the TOTPOTPAR parameters of the model
c ... CM(n,n)  (symmetric) correlation matrix from the fit
c=======================================================================
      INCLUDE 'arrsizes.h'
      INCLUDE 'BLKPOT.h'
      INCLUDE 'BLKDVDP.h'
      INCLUDE 'BLKBOB.h'
      INCLUDE 'BLKBOBRF.h'
      INCLUDE 'BLKCOUNT.h'
c-----------------------------------------------------------------------
      INTEGER I,J,JJ,ISTATE,LAM2,NBETAI
      REAL*8 FU,FLAM,RHT,RMAXT,RDVAL,RDVAL2,RDVALLD, PU(NPARMX),
     1  PT(NPARMX),CM(NPARMX,NPARMX) 
      CHARACTER*20 WRITFILE
c
c------------------------------------------------------------------------
c*** Common Block info for  fununc  calculations ***********************  
      REAL*8 Rsr(NPNTMX,NSTATEMX),Vsr(NPNTMX,NSTATEMX),
     1                                            Bsr(NPNTMX,NSTATEMX)
      INTEGER nPointSR(NSTATEMX)
      COMMON /VsrBLK/Rsr,Vsr,Bsr,nPointSR
c
      REAL*8 Rlr(NPNTMX,NSTATEMX),Plr(NPNTMX,NSTATEMX),
     1                                            Blr(NPNTMX,NSTATEMX)
      INTEGER nPointLR(NSTATEMX)
      COMMON /PlrBLK/Rlr,Plr,Blr,nPointLR
c------------------------------------------------------------------------
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      RHT= RD(2,ISTATE)- RD(1,ISTATE)
      RMAXT= RD(NDATPT(ISTATE),ISTATE)
      WRITE(10,900) 'V(r) ', ISTATE, 'V(r) ', WRITFILE
      WRITE(11,900) 'B(r) ', ISTATE, 'B(r) ', WRITFILE
      DO I= 1,nPointSR(ISTATE),MAX(1,IABS(OSEL(ISTATE))/10)
          WRITE(10,909) Rsr(I,ISTATE),Vsr(I,ISTATE)
ccc       WRITE(11,909) Rsr(I,ISTATE),Bsr(I,ISTATE)
          END DO
      IF(OSEL(ISTATE).LT.0) THEN 
          IF(NUA(ISTATE).GE.0) WRITE(12,900) 'UA(r)',ISTATE,'UA(r)',
     1                                                WRITFILE
          IF(NUB(ISTATE).GE.0) WRITE(13,900) 'UB(r)',ISTATE,'UB(r)',
     1                                                WRITFILE
          IF(NTA(ISTATE).GE.0) WRITE(14,900) 'tA(r)',ISTATE,'qA(r)',
     1                                                WRITFILE
          IF(NTB(ISTATE).GE.0) WRITE(15,900) 'tB(r)',ISTATE,'qB(r)',
     1                                                WRITFILE
          IF(NwCFT(ISTATE).GE.0) THEN
              WRITE(16,902) 'lambda(r)',ISTATE,'lambda(r)',WRITFILE
              LAM2= 2
              IF(IOMEG(ISTATE).GT.0) LAM2= 4*IOMEG(ISTATE)
              ENDIF
          ENDIF
      NBETAI= POTPARF(ISTATE) - Nbeta(ISTATE)
      FU= 0.d0
      DO  I= 1,NDATPT(ISTATE),IABS(OSEL(ISTATE))
          DO J= 1,TOTPOTPAR
              PT(J) = 0.0d0
              END DO
          RDVAL = RD(I,ISTATE)
          RDVAL2= RDVAL*RDVAL
          RDVALLD= RDVAL**LAM2
c ... first ... the potential function itself ...
          DO  J= POTPARI(ISTATE),POTPARF(ISTATE)
              PT(J)= PU(J)*DVtot(J,I)
              ENDDO
          IF(PSEL(ISTATE).GT.0) 
     1           CALL MMCALC(POTPARI(ISTATE),POTPARF(ISTATE),PT,CM,FU)
          WRITE(10,910) RDVAL,VPOT(I,ISTATE),FU
c ... then the exponent coefficient function \betai(i)
ccc       IF(PSEL(ISTATE).LE.5) THEN
ccc           DO  J= NBETAI,POTPARF(ISTATE)
ccc               JJ= J-NBETAI
ccc               PT(J)= PU(J)*DBDB(JJ,I,ISTATE)
ccc               ENDDO
ccc           CALL MMCALC(NBETAI,POTPARF(ISTATE),PT,CM,FU)
ccc           WRITE(11,910) RDVAL,BETAFX(I,ISTATE),FU
ccc           ENDIF
c ... adiabatic BOB correction function for atom-A
          IF(OSEL(ISTATE).LT.0) THEN
              IF(NUA(ISTATE).GE.1) THEN
                  DO  J= UAPARI(ISTATE),UAPARF(ISTATE)
                      PT(J)= PU(J)*DVtot(J,I)
                      ENDDO
                  CALL MMCALC(UAPARI(ISTATE),UAPARF(ISTATE),PT,CM,FU)
                  WRITE(12,910) RDVAL,UAR(I,ISTATE),FU
                  ENDIF
c ... adiabatic BOB correction function for atom-B
              IF(NUB(ISTATE).GE.1) THEN
                  DO  J= UBPARI(ISTATE),UBPARF(ISTATE)
                      PT(J)= PU(J)*DVtot(J,I)
                      ENDDO
                  CALL MMCALC(UBPARI(ISTATE),UBPARF(ISTATE),PT,CM,FU)
                  WRITE(13,910) RDVAL,UBR(I,ISTATE),FU
                  ENDIF
c ... centrifugal BOB correction function for atom-A
              IF(NTA(ISTATE).GE.1) THEN
                  DO  J= TAPARI(ISTATE),TAPARF(ISTATE)
                      PT(J)= PU(J)*DVtot(J,I)*RDVAL2
                      ENDDO
                  CALL MMCALC(TAPARI(ISTATE),TAPARF(ISTATE),PT,CM,FU)
                  WRITE(14,910) RDVAL,TAR(I,ISTATE),FU
                  ENDIF
c ... centrifugal BOB correction function for atom-B
              IF(NTB(ISTATE).GE.1) THEN
                  DO  J= TBPARI(ISTATE),TBPARF(ISTATE)
                      PT(J)= PU(J)*DVtot(J,I)*RDVAL2
                      ENDDO
                  CALL MMCALC(TBPARI(ISTATE),TBPARF(ISTATE),PT,CM,FU)
                  WRITE(15,910) RDVAL,TBR(I,ISTATE),FU
                  ENDIF
c ... Lambda/doublet-sigma doubling correction radial function 
              IF(NwCFT(ISTATE).GE.0) THEN
                  DO  J= LDPARI(ISTATE),LDPARF(ISTATE)
                      PT(J)= PU(J)*DVtot(J,I)*RDVALLD
                      ENDDO
                  CALL MMCALC(LDPARI(ISTATE),LDPARF(ISTATE),PT,CM,FU)
                  FLAM= wRAD(I,ISTATE)*RDVALLD
                  WRITE(16,910) RDVAL,FLAM,FU
                  ENDIF
              ENDIF
          END DO
      DO I= 1,nPointLR(ISTATE)
          WRITE(10,909) Rlr(I,ISTATE),Plr(I,ISTATE)
          WRITE(11,909) Rlr(I,ISTATE),Blr(I,ISTATE)
          END DO
      RETURN
c-----------------------------------------------------------------------
  900 FORMAT(/'variables = "r", "',A5,'", "Uncertainty"'/
     1'zone T = "State',I2,1x,A5,2x,A20,'"',/)
  902 FORMAT(/'variables = "r", "',A9,'", "Uncertainty"'/
     1'zone T = "State',I2,2x,A9,2x,A20,'"',/)
  909 FORMAT(F12.6,1x,1PD22.13,5X,'0.0')
  910 FORMAT(F12.6,1x,1PD22.13,D11.3)
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE MMCALC(JI,JF,PT,CM,FU)
c***********************************************************************
c** For elements JI to JF of column vector PT(j) and rows/columns 
c  JI to JF of the square matrix CM, evaluate  FU**2 = PT^t CM PT
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Externally dimension  PT(NPARMX) and CM(NPARMX,NPARMX)
c***********************************************************************
      INCLUDE 'arrsizes.h'
      INTEGER I,J,JI,JF
      REAL*8  PT(NPARMX), CM(NPARMX,NPARMX), TVAL, FU
c=======================================================================
c** Initialize the variables.
      FU = 0.0d0
      DO  J= JI, JF
          TVAL= 0.d0
          DO  I= JI,JF
c** Multiply the vector PT with the correlation matrix (CM).
              TVAL= TVAL + PT(I)*CM(I,J)
              ENDDO
c** Now to multiply the vector (PT*CM) with the vector PT.
          FU= FU + TVAL * PT(J)
          ENDDO
c** Complete calculation of the uncertainty of the function.
      FU= DSQRT(FU)
      RETURN
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
