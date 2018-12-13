c***********************************************************************
      SUBROUTINE MAPPAR(NISTP,PV,FORBAC)
c***********************************************************************
c** This subroutine will convert external logical physical parameters 
c  into the generic NLLSSRR parameter array PV or the reverse.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++ COPYRIGHT 1997-2011  by  J.Y. Seto & R.J. Le Roy  (ver. 14/07/2011)+
c   Dept. of Chemistry, Univ. of Waterloo, Waterloo, Ontario, Canada   +
c    This software may not be sold or any other commercial use made    +
c      of it without the express written permission of the authors.    +
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c    PV(i)   is the NLLSSRR parameter array.
c    FORBAC  is a flag to determine which way the parameters are mapped
c            FORBAC = 0 : Map internal PV to external variables.
c            FORBAC = 1 : Map external varuables to internal PV.
c*   NSTATES is the number of states being considered (in BLKPARAM)
c=======================================================================
      INCLUDE 'arrsizes.h'
      INCLUDE 'BLKPOT.h'
      INCLUDE 'BLKPARAM.h'
      INCLUDE 'BLKBOB.h'
c-----------------------------------------------------------------------
      INTEGER NISTP, m, FORBAC, ISTATE, IPV, I, J, ISOT
      REAL*8 PV(NPARMX)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Map external free parameters (De, Re, etc.) onto internal NLLSSRR
c  parameters PV(j)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IF(FORBAC.EQ.0) THEN
          IPV= 0
          DO  ISTATE= 1,NSTATES
              IF(PSEL(ISTATE).EQ.-1) THEN
c*** Manage parameters for term value mappings ...
                  DO  ISOT= 1, NISTP
                      DO  I= VMIN(ISTATE,ISOT),VMAX(ISTATE,ISOT)
                          IF(NBC(I,ISOT,ISTATE).GT.0) THEN
                              DO  J=1,NBC(I,ISOT,ISTATE)
                                  IPV= IPV+1
                                  PV(IPV)= ZBC(I,J-1,ISOT,ISTATE)
                                  ENDDO
                              IF(NQC(I,ISOT,ISTATE).GT.0) THEN  
                                  DO  J=1,NQC(I,ISOT,ISTATE)
                                      IPV= IPV+1
                                      PV(IPV)= ZQC(I,J-1,ISOT,ISTATE)
                                      ENDDO
                                   ENDIF
                              ENDIF  
                          ENDDO
                      DO  I= VMIN(ISTATE,ISOT),VMAX(ISTATE,ISOT)
                          ENDDO
                      ENDDO
                  ENDIF
              IF(PSEL(ISTATE).GT.0) THEN
c*** Manage parameters for potential function mapping ...
                  IF((PSEL(ISTATE).LE.6).and.(PSEL(ISTATE).NE.4)) THEN
                      IPV= IPV+ 1
                      PV(IPV)= DE(ISTATE)
                      ENDIF
                  IPV= IPV+ 1
                  PV(IPV)= RE(ISTATE)
                  IF((PSEL(ISTATE).GE.2).AND.(PSEL(ISTATE).LE.6)) THEN
                      DO  m= 1,NCMM(ISTATE)
                          IPV= IPV+ 1
                          PV(IPV)= CmVAL(m,ISTATE)
                          ENDDO
                      ENDIF
                  J= 0
                  IF(NSR(ISTATE).LT.0) J=1   !! for Pashov spline exponent
                  DO  I= J,Nbeta(ISTATE)
                      IPV= IPV+ 1
                      PV(IPV)= BETA(I,ISTATE)
                      ENDDO
                  IF(NUA(ISTATE).GE.0) THEN
                      DO  I= 0,NUA(ISTATE)
                          IPV= IPV+ 1
                          PV(IPV) = UA(I,ISTATE)
                          ENDDO
                      ENDIF
                  IF(NUB(ISTATE).GE.0) THEN
                      DO  I= 0,NUB(ISTATE)
                          IPV= IPV+ 1
                          PV(IPV) = UB(I,ISTATE)
                          ENDDO
                      ENDIF
                  IF(NTA(ISTATE).GE.0) THEN
                      DO  I= 0,NTA(ISTATE)
                          IPV= IPV+ 1
                          PV(IPV) = TA(I,ISTATE)
                          ENDDO
                      ENDIF
                  IF(NTB(ISTATE).GE.0) THEN
                      DO  I= 0,NTB(ISTATE)
                          IPV= IPV+ 1
                          PV(IPV) = TB(I,ISTATE)
                          ENDDO
                      ENDIF
                  IF(NwCFT(ISTATE).GE.0) THEN
                      DO  I= 0, NwCFT(ISTATE)
                          IPV= IPV+ 1
                          PV(IPV) = wCFT(I,ISTATE)
                          ENDDO
                      ENDIF
                  ENDIF
              ENDDO
        ELSEIF(FORBAC.EQ.1) THEN
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Convert internal NLLSSRR parameter array back into external
c   (logical) variable system (De, Re, etc.).
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          IPV = 0
          DO  ISTATE=1,NSTATES
              IF(PSEL(ISTATE).EQ.-1) THEN
c*** Manage parameters for term value mappings ...
                  DO  ISOT= 1, NISTP
                      DO  I= VMIN(ISTATE,ISOT),VMAX(ISTATE,ISOT)
                          IF(NBC(I,ISOT,ISTATE).GT.0) THEN
                              DO  J= 1,NBC(I,ISOT,ISTATE)
                                  IPV= IPV+1
                                  ZBC(I,J-1,ISOT,ISTATE)= PV(IPV)
                                  ENDDO
                              IF(NQC(I,ISOT,ISTATE).GT.0) THEN
                                  DO  J= 1,NQC(I,ISOT,ISTATE)
                                      IPV= IPV+1
                                      ZQC(I,J-1,ISOT,ISTATE)= PV(IPV)
                                      ENDDO
                                  ENDIF
                              ENDIF
                          ENDDO
                      ENDDO
                  ENDIF
              IF(PSEL(ISTATE).GT.0) THEN
c*** Manage parameters for potential function mappings ...
                  IF((PSEL(ISTATE).LE.6).and.(PSEL(ISTATE).NE.4)) THEN
                      IPV= IPV + 1
                      DE(ISTATE)= PV(IPV)
                      ENDIF
                  IPV= IPV + 1
                  RE(ISTATE) = PV(IPV)
                  IF((PSEL(ISTATE).GE.2).AND.(PSEL(ISTATE).LE.6)) THEN
                      DO  m= 1,NCMM(ISTATE)
                          IPV= IPV+ 1
                          CmVAL(m,ISTATE)= PV(IPV) 
                          ENDDO
                      ENDIF
                  J=0
                  IF(NSR(ISTATE).LT.0) J=1   !! for Pashov spline exponent
                  DO I= J, Nbeta(ISTATE)
                      IPV = IPV + 1
                      BETA(I,ISTATE) = PV(IPV)
                      ENDDO
                  IF(NUA(ISTATE).GE.0) THEN
                      DO I= 0,NUA(ISTATE)
                          IPV = IPV + 1
                          UA(I,ISTATE) = PV(IPV)
                          ENDDO
                      ENDIF
                  IF(NUB(ISTATE).GE.0) THEN
                      DO I= 0,NUB(ISTATE)
                          IPV = IPV + 1
                          UB(I,ISTATE) = PV(IPV)
                          ENDDO
                      ENDIF
                  IF(NTA(ISTATE).GE.0) THEN
                      DO I= 0,NTA(ISTATE)
                          IPV = IPV + 1
                          TA(I,ISTATE) = PV(IPV)
                          ENDDO
                      ENDIF
                  IF(NTB(ISTATE).GE.0) THEN
                      DO I= 0,NTB(ISTATE)
                          IPV = IPV + 1
                          TB(I,ISTATE) = PV(IPV)
                          ENDDO
                      ENDIF
                  IF(NwCFT(ISTATE).GE.0) THEN
                      DO I=0,NwCFT(ISTATE)
                          IPV = IPV + 1
                          wCFT(I,ISTATE) = PV(IPV)
                          ENDDO
                      ENDIF
                  ENDIF
              ENDDO
        ENDIF
      RETURN
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
