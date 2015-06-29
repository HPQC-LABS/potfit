c***********************************************************************
      SUBROUTINE SCECOR(KV,KVLEV,JROT,INNER,ICOR,IWR,EO,RH,BFCT,NDP,
     1                                    NCN,V,BMAX,VMAXX,VLIM,DGDV2)
c** Subroutine calculates (approximate!) semiclassical estimate of 
c  dG/dv for level  v= KV  with energy  EO [cm-1]  on potential 
c  {V(i),i=1,NDP} (in 'internal BFCT units' {V[cm-1]*BFCT}), and uses
c  those results to estimate energy of level  KVLEV (usually = KV+1)
c** If the 'clever' semiclassical procedure fails - try a brute force
c  step-by-step search, using alternately INNER & OUTER well starting
c** BMAX is internal barrier maximum energy for double-well case, 
c   and very large negative number for single-well potential
c** VMAXX is height of outermost maximum, or VLIM for barrierless case
c** On return, negative DGDV2 signals error!  No phase integrals found
c=======================================================================
c                   Version date:  14 Sept 2014
c***********************************************************************
      INTEGER I,II,I1,I2,I3,I4,IV1,IV2,INNER,ICOR,JROT,KV,KVB,KVLEV,
     1  KVDIF,NDP,NCN,IDIF,BRUTE,IB,IWR,NPMAX
      REAL*8 EO,DE0,RH,BFCT,ARG2,ARG3,EINT,VPH1,VPH2,DGDV1,DGDV2,DGDVM,
     1  DGDV2P,DGDVB,DGDVBP,EBRUTE,DEBRUTE,DE1,DE2,Y1,Y2,Y3,RT,ANS1,dv1,
     2  dv2,ANS2,XDIF,VLIM,VMAXX,BMAX,PNCN,PWCN,PP1,VDMV,ENEXT,V(NDP)
      SAVE BRUTE,EBRUTE,DEBRUTE,DGDVB
      DATA DGDVB/-1.d0/,KVB/-1/
c
      DGDV2= -1.d0
      EINT= EO*BFCT
      IF(KVLEV.EQ.0) DGDVB= -1.d0
      KVDIF= KVLEV- KV
      IF(ICOR.EQ.1) BRUTE= 0
      PWCN= 2.d0*NCN/DABS(NCN- 2.d0)
      PNCN= ABS(NCN-2)/DFLOAT(NCN+2)
      DGDVBP= DGDVB**PNCN
      PP1= 1.d0/pNCN + 1.d0
      I3= NDP
      IF(EO.GT.VLIM) THEN
c*** For Quasibound levels, first search inward to classically forbidden
          PWCN= 2.d0
          PNCN= 1.d0 
          PP1= 1.d0
          DO  I= NDP,1,-1
              I3= I
              IF(V(I).GT.EINT) GOTO 8
              ENDDO
          ENDIF
c*** Now, search inward for outermost well turning point
    8 DO  I= I3,1,-1
          I4= I
          IF(V(I).LT.EINT) GOTO 10
          ENDDO
c*** If never found an 'outer' turning point (e.g., above qbdd. barier)
c  then simply return with negative  DGDV2  as error flag
      RETURN
c... Now collect vibrational phase and its energy deriv. over outer well
   10 Y1= EINT- V(I4+1)
      Y2= EINT- V(I4)
      Y3= EINT- V(I4-1)
      CALL LEVQAD(Y1,Y2,Y3,RH,RT,ANS1,ANS2)
      ARG2= DSQRT(Y3)
      VPH2= 0.5d0*ARG2 + ANS2/RH
      DGDV2= 0.5d0/ARG2 + ANS1/RH
      DO  I= I4-2,1,-1
c... here collect (v+1/2) and dv/dG integrals over outer well ....
          II= I
          IF(V(I).GT.EINT) GO TO 12
          ARG3= ARG2
          ARG2= DSQRT(EINT - V(I))
          VPH2= VPH2+ ARG2
          DGDV2= DGDV2+ 1.d0/ARG2
          ENDDO
   12 I3= II+1
      Y1= EINT- V(I3-1)
      Y2= EINT- V(I3)
      Y3= EINT- V(I3+1)
      CALL LEVQAD(Y1,Y2,Y3,RH,RT,ANS1,ANS2)
      VPH2= (VPH2 - ARG2 - 0.5d0*ARG3 + ANS2/RH)/3.141592654d0
      DGDV2= DGDV2 -1.d0/ARG2 - 0.5d0/ARG3 + ANS1/RH
      DGDV2= 6.283185308d0/(BFCT*DGDV2)
c*** Next, search outward from RMIN for innermost turning point
      DO  I= 1,NDP 
          I1= I
          IF(V(I).LT.EINT) GOTO 20
          ENDDO
   20 IF(I1.EQ.1) THEN
c... but if RMIN is in the classically allowed region ... STOP here          
          WRITE(6,602) JROT,EO
          STOP
          ENDIF
      IF(I1.GE.I3) THEN
c*** For single-well potential or above barrier of double-well potential
c   use N-D theory estimate based on 'vD-v' from ratio  of Eb to dG/dv
          VDMV= PWCN*(VMAXX-EO)/DGDV2
          ENEXT= VMAXX - (VMAXX-EO)*((VDMV- KVDIF)/VDMV)**PWCN
          IF(IWR.GE.2) THEN
              IF(ABS(EO).GT.1.d0) WRITE(6,600) ICOR,KV,JROT,EO, 
     1                                                VPH2-0.5d0,DGDV2
              IF(ABS(EO).LE.1.d0) WRITE(6,601) ICOR,KV,JROT,EO,
     1                                                VPH2-0.5d0,DGDV2
              WRITE(6,606) VDMV,ENEXT
              ENDIF    
ccccccc????? Redundant stuff now ???????????????????????????????????????
cc        IF((KV.LT.(KVLEV-1)).AND.(DGDVB.GT.0.d0)) THEN
c... If got wrong level (KV not one below KVLEV) and NOT first call ...
cc            IF((EO-BMAX).GT.(2.d0*DGDV2)) THEN
c  For eneries well above the barrier of a double minimum potrnti
c... 'Normal' case: use B-S plot area to estimate correct energy
cc                DE0= KVDIF*(DGDV2- 0.5d0*(DGDV2-DGDVB)/DFLOAT(KV-KVB))
cc                EO= EO+ DE0 
cc                KV= KVB
cc                KVLEV= KV+1
cc                RETURN
cc              ELSE
c... but close to barrier in double-well potential, switch to 'BRUTE'
cc                BRUTE=BRUTE+ 1
cc                DGDV1= DGDV2
cc                XDIF= SIGN(1,KVDIF)  
cc                GOTO 54
cc              ENDIF
cc            ENDIF
ccccccc????? Redundant stuff now ???????????????????????????????????????
          IF(VDMV.LT.1.d0) THEN
              ICOR= 100
              IF(IWR.GT.0) WRITE(6,604) KV,EO
            ELSE
               EO= ENEXT
            ENDIF   
          DGDVB= DGDV2
          DGDVBP= DGDVB**PNCN
          KVB= KV
          INNER= 0
          RETURN
          ENDIF
c
c*** For a double-well potential, now collect vibrational phase and its 
c   energy derivative over the inner well
      Y1= EINT- V(I1-1)
      Y2= EINT- V(I1)
      Y3= EINT- V(I1+1)
      CALL LEVQAD(Y1,Y2,Y3,RH,RT,ANS1,ANS2)
      ARG2= DSQRT(Y3)
      VPH1= 0.5d0*ARG2 + ANS2/RH
      DGDV1= 0.5d0/ARG2 + ANS1/RH
      DO  I= I1+2,NDP
c... now, collect integral and count nodes outward to second turning point ...
          IF(V(I).GT.EINT) GO TO 22
          ARG3= ARG2
          ARG2= DSQRT(EINT - V(I))
          VPH1= VPH1+ ARG2
          DGDV1= DGDV1+ 1.d0/ARG2
          ENDDO
   22 I2= I-1
      Y1= EINT- V(I2+1)
      Y2= EINT- V(I2)
      Y3= EINT- V(I2-1)
      CALL LEVQAD(Y1,Y2,Y3,RH,RT,ANS1,ANS2)
      VPH1= (VPH1 - ARG2 - 0.5d0*ARG3 + ANS2/RH)/3.141592654d0
      DGDV1= DGDV1 -1.d0/ARG2 - 0.5d0/ARG3 + ANS1/RH
      DGDV1= 6.28318531d0/(BFCT*DGDV1)
      DGDVM= DGDV1*DGDV2/(DGDV1+DGDV2)
      IF(KVDIF.EQ.0) THEN
c** If already at level sought, return
          IF(IWR.GE.2) WRITE(6,610) KV,JROT,EO,VPH1-0.5d0,DGDV1,KVLEV,
     1                                           ICOR,VPH2-0.5d0,DGDV2
          RETURN
          ENDIF
c
c** Not at right level - Check whether looking for higher or lower level ...
      IDIF= SIGN(1,KVDIF)
      XDIF= IDIF
      IF((ICOR.GE.3).AND.((IABS(KVDIF).EQ.1).OR.(BRUTE.GT.0))) GOTO 50
c*** 'Conventional' semiclassical search for nearest INNER or OUTER well level
c... first, determine whether starting level KV was really INNER or OUTER
      dv1= (VPH1-0.5d0) - NINT(VPH1-0.5d0)
      dv2=(VPH2-0.5d0) - NINT(VPH2-0.5d0) 
      IF((DABS(dv2).GT.0.1).AND.(DABS(dv1).LT.0.1)) THEN
          INNER=1
          ENDIF
      IF(INNER.EQ.0) THEN
c... and if current energy EO is for an outer-well level ...
          DE2= DGDV2*XDIF
          IF(IDIF.GT.0) DE1= (Ceiling(VPH1-0.5d0) - (VPH1-0.5d0))*DGDV1
          IF(IDIF.LE.0) DE1= -((VPH1-0.5d0)- Floor(VPH1-0.5d0))*DGDV1
          IF(IWR.GE.2) WRITE(6,610) KV,JROT,EO,VPH1-0.5d0,DGDV1,KVLEV,
     1                                           ICOR,VPH2-0.5d0,DGDV2
        ELSE
c... and if current energy EO is for an inner-well level ...
          DE1= DGDV1*XDIF
          IF(IDIF.GT.0) DE2= (Ceiling(VPH2-0.5d0) - (VPH2-0.5d0))*DGDV2
          IF(IDIF.LE.0) DE2= -(1.d0 - dv2)*DGDV2
          IF(IWR.GE.2) WRITE(6,610) KV,JROT,EO,VPH1-0.5d0,DGDV1,KVLEV,
     1                                           ICOR,VPH2-0.5d0,DGDV2
        ENDIF  
      IF(DABS(DE2).LT.DABS(DE1)) THEN
c... for case in which predict that next level will be OUTER
          INNER= 0
          EO= EO+ DE2
        ELSE
c... for case in which predict that next level will be INNER
          INNER= 1
          EO= EO+ DE1
        ENDIF
      RETURN
   50 BRUTE= BRUTE+ 1
c*** Now .. Brute force search for desired level !
      IF(IWR.GE.2) WRITE(6,610) KV,JROT,EO,VPH1-0.5d0,DGDV1,KVLEV,
     1                                           ICOR,VPH2-0.5d0,DGDV2
   54 IF(BRUTE.EQ.1) THEN
c... in first brute-force step, use previous energy with opposite INNER
          EBRUTE= EO 
          IF(INNER.EQ.0) THEN
              INNER= 1
            ELSE
              INNER= 0
            ENDIF
          DEBRUTE= DMIN1(DGDV1,DGDV2)*XDIF*0.3d0
          RETURN
          ENDIF
      IB= BRUTE/2
c... in subsequent EVEN steps, lower EO by DEBRUTE/10 for same INNER
      IF((IB+IB).EQ.BRUTE) THEN
          EBRUTE= EBRUTE+ DEBRUTE
          EO= EBRUTE
          RETURN
        ELSE
c... in subsequent ODD steps, lower repeat previous EO with INNER changed
          IF(INNER.EQ.0) THEN
              INNER= 1
            ELSE
              INNER= 0
            ENDIF
          EO= EBRUTE
          RETURN
        ENDIF
c     RETURN
  600 FORMAT(' Single well  ICOR=',I2,':  E(v=',i3,',J=',I3,')=',f10.2,
     1 '  v(SC)=',F8.3,'  dGdv=',f8.3)
  601 FORMAT(' Single well  ICOR=',I2,':  E(v=',i3,',J=',I3,')=',
     1  1PD12.4,'  v(SC)=',0PF8.3, /63x,'dGdv=',1PD12.4)
  602 FORMAT(/' *** ERROR ***  V(1) < E(J=',i3,')=',f10.2 )
  604 FORMAT(10x,'Find highest level of this potential is   E(v=',i3,
     1                                                   ')=',1PD18.10)
  606 FORMAT(40x,'(vD-v)=',f10.4,'   E(next)=',1PD12.4)
  610 FORMAT(' Double well   E(v=',i3,', J=',I3,')=',f9.3,
     1 ':   v1(SC)=',F7.3,'   dGdv1=',f8.2/8x,'seeking  v=',I3,
     2 ' (ICOR=',I2,')',8x,':   v2(SC)=',F7.3,'   dGdv2=',f8.2 )
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

