c***********************************************************************
      SUBROUTINE dampF(r,rhoAB,NCMM,MMLR,IVSR,IDSTT,DM,DMP,DMPP)
c** Subroutine to generate values 'Dm' and its first `Dmp' and second
c   'Dmpp' derivatives w.r.t. r of the chosen form of the damping
c    function, for  m= 1 to MMAX.
c---------------------- RJL Version of 05 August 2014 ------------------
c-----------------------------------------------------------------------
c                 Upon Input
c* r - the radial distance in Angsroms (!) 
c* RHOab  'universal' scaling coefficient used for systems other than H_2
c       RHOab= 2*(RHOa*RHOb)/(RHOa+RHOb) where RHOa = (I_p^A/I_p^H)^0.66
c              where I_p^A is the ionization potential of atom A
c              and I_p^H is the ionization potential of atomic hydrogen
c* NCMM  the number of inverse-power terms to be considered
c* MMLR  are the powers of the NCMM inverse-power terms
c* IVSR defines damping s.th.  Dm(r)/r^m --> r^{IVSR/2} as r --> 0
c* IDSTT specifies damping function type:  > 0  use Douketis et al. form 
c                               if  IDSTT .LE. 0  use Tang-Toennies form
c-----------------------------------------------------------------------
c                 Upon Output
c  DM(m) - The value of the damping function for the long range term 
c          C_MMLR(m)/r^MMLR(m)    {m= 1, NCMM}
c  DMP(m): first derivative of the damping function  DM(m) w.r.t. r
c  DMPP(m): second derivative of the damping function  DM(m) w.r.t. r
c-----------------------------------------------------------------------
      INTEGER NCMM,NCMMax,MMLR(NCMM),IVSR,IDSTT,IVSRF,FIRST, Lsr,m,
     1  MM,MMAX,MMTEMP
      REAL*8 r,rhoAB,bTT(-2:2),cDS(-4:4),bDS(-4:4),aTT,br,XP,YP,
     1  TK, DM(NCMM),DMP(NCMM),DMPP(NCMM),SM(-3:25),
     2  bpm(20,-4:0), cpm(20,-4:0),ZK
c------------------------------------------------------------------------
c  The following values for the numerical factors used in both TT and DS
c  were  normalized to the Hydrogen data presented
c  by Kreek and Meath in J.Chem.Phys. 50, 2289 (1969).
c  The ratio has been chosen such that  b= FACTOR*(I_p^X / I_p^H)^{2/3}
c  for the homoatomic diatomic species X_2, where I_p^A is the ionization
c------------------------------------------------------------------------
       DATA bTT/2.10d0,2.44d0,2.78d0,3.13d0,3.47d0/
       DATA bDS/2.50d0,2.90d0,3.30d0,3.69d0,3.95d0,0.d0,4.53d0,0.d0,
     1          4.99d0/
       DATA cDS/0.468d0,0.446d0,0.423d0,0.405d0,0.390d0,0.d0,0.360d0,
cc     DATA cDS/0.468d0,0.446d0,0.423d0,0.400d0,0.390d0,0.d0,0.360d0,
     1            0.d0,0.340d0/
       DATA FIRST/ 1/
       SAVE FIRST, bpm, cpm
c------------------------------------------------------------------------
       MMTEMP = MMLR(1)              ! temporary substn.: corrected@end
       IF(MMLR(1).LE.0) MMLR(1) = 1  ! enforce stability for A-F cases
      IF(RHOab.LE.0) THEN
          DO  m=1,NCMMax
              DM(m)=1.d0
              DMP(m)= 0.d0
              DMPP(m)= 0.d0
              ENDDO
          WRITE(6,602) RHOab
          RETURN
          ENDIF
      IF(IDSTT.LE.0) THEN
c===========================================
c** For Tang-Toennies type damping functions
c===========================================
          Lsr= IVSR/2
          IF((IVSR.LT.-4).OR.(IVSR.GT.4).OR.((2*LSR).NE.IVSR)) THEN
                WRITE(6,600) 'TT',IVSR
                STOP
                ENDIF
          MMAX= MMLR(NCMM) + Lsr - 1
          aTT= RHOab*bTT(Lsr)
          br= aTT*r
          XP= DEXP(-br)
          SM(-3)= 0.d0
          SM(-2)= 0.d0
          SM(-1)= 0.d0
          SM(0)=  1.d0
          TK= 1.d0
          IF(br.GT.0.5d0) THEN
              DO  m= 1,MMAX
                  TK= TK*br/DFLOAT(m)
                  SM(m)= SM(m-1)+ TK
                  ENDDO
              DO m= 1, NCMM
                  MM= MMLR(m) - 1 + Lsr
                  DM(m)= 1.d0 - XP*SM(MM)
                  DMP(m)= aTT*XP*(SM(MM) - SM(MM-1))
                  DMPP(m)= -aTT*aTT*XP*(SM(MM) 
     1                                     - 2.d0*SM(MM-1) + SM(MM-2))
                  ENDDO
c-----------------------------------------------------------------------
c  The above section handles the calculation of the value of the damping
c  function for most values of r.  However, at very small r that algorithm
c  becomes unstable due to numerical noise.  To avoid this, if the 
c  argument is very small it is re-evaluated as a finite sum ...
c-----------------------------------------------------------------------
            ELSE
              MMAX= MMAX+5
              DO  m= 1, MMAX
c... NOTE that here SM(m) is the m'th term  (b*r)^m/m!  [not a sum]
                  SM(m)= SM(m-1)*br/DFLOAT(m)
                  ENDDO
              DO  m= 1, NCMM
                  MM= MMLR(m) + Lsr
                  DM(m)= XP*(SM(MM)+ SM(MM+1)+ SM(MM+2)+ SM(MM+3) 
     1                                                     + SM(MM+4))
                  DMP(m)= aTT*XP*SM(m-1)
                  DMPP(m)= aTT*aTT*XP*(SM(m-2)-SM(m-1))
                  ENDDO
            ENDIF
          ENDIF
c
      IF(IDSTT.GT.0) THEN
c=======================================================================
c** For Douketis-Scoles-Marchetti-Zen-Thakkar type damping function ...
c=======================================================================
          IF((IVSR.LT.-4).OR.(IVSR.GT.4).OR.(IVSR.EQ.1).OR.(IVSR.EQ.3)) 
     1    THEN
              WRITE(6,600) 'DS',IVSR
              STOP
              ENDIF
          IF(FIRST.EQ.1) THEN
              DO m= 1, 20
                  DO  IVSRF= -4,0
                      bpm(m,IVSRF)= bDS(IVSRF)/DFLOAT(m)
                      cpm(m,IVSRF)= cDS(IVSRF)/DSQRT(DFLOAT(m))
                      ENDDO
                  ENDDO
              FIRST= 0 
              ENDIF
          br= rhoAB*r
          DO m= 1, NCMM
              MM= MMLR(m)
              XP= DEXP(-(bpm(MM,IVSR) + cpm(MM,IVSR)*br)*br)
              YP= 1.d0 - XP
              ZK= MM + 0.5d0*IVSR
              DM(m)= YP**ZK
              TK= (bpm(MM,IVSR) + 2.d0*cpm(MM,IVSR)*br)*rhoAB
              DMP(m) = ZK*XP*TK*DM(m)/YP
c ... calculate second derivative [for DELR case] {check this!}
              DMPP(m)= (ZK-1.d0)*DMP(m)*(XP*TK)/YP
     1               - DMP(m)*TK + DMP(m)*2.d0*cpm(MM,IVSR)*rhoAB**2/TK
              ENDDO   
          ENDIF  
      MMLR(1) = MMTEMP
      RETURN
  600 FORMAT(/,' *** ERROR ***  For  ',A2,'-damping functions not yet de
     1fined for   IVSR=',i3)
  602 FORMAT( /,' ***ERROR ***  should not call dampF when rhoAB=',
     1  F7.4)
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
