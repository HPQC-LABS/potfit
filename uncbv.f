c***********************************************************************
      SUBROUTINE UNCBV(NPTOT,PV,PU,CM)
c***********************************************************************
c**  Subroutine to compute the uncertainties in calcuated Bv values, 
c  1) by the finite difference R(0) approximation
c  2) from the LIDE expression
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** On entry:
c ... OSEL     print values at every OSEL'th mesh point
c ... ISTATE   electronic state counter
c ... PU(n)    uncertainties in the TOTPOTPAR parameters of the model
c ... CM(n,n)  (symmetric) correlation matrix from the fit
c=======================================================================
      INCLUDE 'arrsizes.h'
      INCLUDE 'BLKPOT.h'
      INCLUDE 'BLKDVDP.h'
      INCLUDE 'BLKBOB.h'
      INCLUDE 'BLKPARAM.h'
      INCLUDE 'BLKISOT.h'
      INCLUDE 'BLKBOBRF.h'
      INCLUDE 'BLKCOUNT.h'
      INCLUDE 'BLKDATA.h'
c-----------------------------------------------------------------------
      INTEGER IISTP, IV, I, J, JROT,KVLEV,efPARITY,NPTOT,NB1,COUNT1,
     1    ISTATE, NBEG,NEND,INNODE,INNER,IWR,LPRWF,WARN, fcount
      REAL*8 PV(NPARMX), PU(NPARMX), PD(NPARMX), CM(NPARMX,NPARMX), 
     1  PQ(NPARMX),PT(NPARMX),DEDPK(HPARMX),V1D(NPNTMX),SWF(NPNTMX),
     2  DVDPK(NPNTMX),Bunc1, Bunc2, Eunc2, CALC, EIV, FWHM, BFCT,BvWN,
     3  EO,Bv,UMAX,dBdPk
c** First - set up fake R(0) datum
      WARN= 1
      COUNT1= COUNTOT+1
      NB1= NBANDTOT+1
      IB(COUNT1)= NB1
      JPP(COUNT1)= 0
      JP(COUNT1)= 1
      FREQ(COUNT1)= 0.d0
      UFREQ(COUNT1)= 0.d0
      DO ISTATE= 1, NSTATES
          IF(PSEL(ISTATE).LE.0) CYCLE
          IEP(NB1)= ISTATE
          IEPP(NB1)= ISTATE
          DO IISTP= 1,NISTP
              write(7,600) 
              ISTP(NB1)= IISTP
cc            WRITE(17,600) IISTP, SLABL(ISTATE)
c** Generate 1-D potential for this state/isotopologue to prepare ...
              BFCT=RH(ISTATE)*RH(ISTATE)*ZMASS(3,IISTP)/16.857629206d0
              BvWN= 16.857629206D0/ZMASS(3,IISTP)
              DO I= 1,NDATPT(ISTATE)
                  V1D(I)= BFCT*(VPOT(I,ISTATE) 
     1                              + ZMUA(IISTP,ISTATE)*UAR(I,ISTATE)
     2                             + ZMUB(IISTP,ISTATE)*UBR(I,ISTATE))
                  ENDDO
c** for each vibrational level of each isotopologue ....
              DO  IV= VMIN(ISTATE,IISTP), VMAX(ISTATE,IISTP)
                  VP(NB1)= IV
                  VPP(NB1)= IV
c** Next get partial derivatives from Bv as half of a pure R(0) energy
                  CALL DYIDPJ(COUNT1,COUNT1,NPTOT,FREQ(COUNT1),CALC,PV,
     1                                                             PD)
c ... set up column vector for uncertainty calculation
                  DO  J= POTPARI(ISTATE), POTPARF(ISTATE)
                      PT(J)= 0.5*PD(J)*PU(J)
                      ENDDO
c... Uncertainty calculation for this level via "R(0)" approvimation 
                  CALL MMCALC(1,NPTOT,PT,CM,Bunc1)
c   
cc                WRITE(7,602)  IV,0.5d0*CALC/ZK(IV,1,IISTP,ISTATE),
cc   1                                                         Bunc1
c... Call SCHRQ to get energy and wavefunction for exact calculation
                  JROT= 0
                  KVLEV= IV
                  EO= ZK(IV,0,IISTP,ISTATE)
                  Bv= ZK(IV,1,IISTP,ISTATE)
                  CALL SCHRQ(IV,JROT,EO,FWHM,UMAX,VLIM(ISTATE),V1D,SWF,
     1                       BFCT,EPS(ISTATE),RMIN(ISTATE),RH(ISTATE),
     2                NDATPT(ISTATE),NBEG,NEND,INNODE,INNER,IWR,LPRWF)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                  efPARITY= 0
c ... call DEDP to get eigenvalue derivatives DEDPK(j) for E & Bv unc.
                  CALL DEDP(COUNT1,ISTATE,IISTP,ZMASS(3,IISTP),KVLEV,
     1                JROT,efPARITY,EO,VLIM(ISTATE),FWHM,DEDPK,fcount)
c... Now Loop over potential parameters for this state, calculating
                  DO J= POTPARI(ISTATE), HPARF(ISTATE)
                      PQ(J)= DEDPK(J)*PU(J) 
c... first - define 1D partial derivative array
                      DO I= 1,NDATPT(ISTATE)
                          DVDPK(I)= BFCT*dVtot(J,I)
                          ENDDO
cc                    CALL dPSIdp(ISTATE,IISTP,EO,NBEG,NEND,
cc                ?? problem to solve in my 'spare time'
cc   1               NDATPT(ISTATE),BvWN,V1D,SWF,DEDPK(J),dBdPk,dVdPk)
cc
cc                    WRITE(17,604)  J,0.5*PD(J),dBdPk
cc
                      PT(J)= dBdPk*PU(J)
                      CONTINUE
                      ENDDO
c... get eigenvalue Gv uncertainties
                  CALL MMCALC(1,NPTOT,PQ,CM,Eunc2)
c... get Bv uncertainties
                  CALL MMCALC(1,NPTOT,PT,CM,Bunc2)
ccc
cc                WRITE(17,610) IV,EO,Eunc2,Bv,Bunc1
                  WRITE(7,610) IV,EO,Eunc2,Bv,Bunc1
cc                WRITE(17,608) IV,Bunc1,Bunc2
ccc               WRITE(7,608) IV,Bunc1,Bunc2
                  ENDDO
              ENDDO
          ENDDO
      RETURN
  600 FORMAT(:/'  Predict Bv uncertainties for isotopologue-',I1,' in st
     1ate ',A3)
  602 FORMAT(' For  v=',I3,'    Bv{R(0)/2}/Bv(exact)=',F13.10,
     1   '       unc(Bv)=',1PD13.6)
  604 FORMAT('   for param(',I2,')    dBvdp(approx)=',1P1D13.6,
     1  '   Bvdp(numerical)=',D13.6)
  608 FORMAT(' For  v=',I3,'    u(Bv{R(0)/2})=',F13.10,
     1   '    u(Bv;numerical)=',1PD13.6)
  610 FORMAT(' v=',I3,'    E=', F12.4,'  u(E)=',F10.6,'  Bv=',f11.8,
     1  '  u(Bv)=',1PD13.6)
      END
c234567890 234567890 234567890 234567890 234567890234567890 234567890

