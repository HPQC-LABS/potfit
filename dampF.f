c***********************************************************************
      SUBROUTINE dampF(r,rhoAB,NCMM,MMLR,IDF,IDSTT,KDER,DM,DMP,DMPP)
c** Subroutine to generate values 'Dm' and its first `Dmp' and second
c   'Dmpp' derivatives w.r.t. R of the chosen form of the damping
c    function, for  m= 1 to MMAX.
c---------------------- RJL Version of 13 January 2013 -----------------
c-----------------------------------------------------------------------
c                 Upon Input
c* r - the radial distance in Angsroms (!) 
c* RHOab  'universal' scaling coefficient used for systems other than H_2
c       RHOab= 2*(RHOa*RHOb)/(RHOa+RHOb) where RHOa = (I_p^A/I_p^H)^0.66
c              where I_p^A is the ionization potential of atom A
c              and I_p^H is the ionization potential of atomic hydrogen
c* NCMM  the number of inverse-power terms to be considered
c* MMLR  are the powers of the NCMM inverse-power terms
c* IDF requires damping to be defined s.th.  Dm(r)/r^m --> r^{IDF/2}
c* IDSTT specifies damping function type:  > 0  use Douketis et al. form 
c                               if  IDSTT .LE. 0  use Tang-Toennies form
c* KDER:  if KDER.GT.0  the first derivative is also calculated 
c*        if KDER.GT.1  the second derivative is also calculated 
c-----------------------------------------------------------------------
c                 Upon Output
c  DM(m) - The value of the damping function for the long range term 
c          C_MMLR(m)/r^MMLR(m)    {m= 1, NCMM}
c  DMP(m) - The first derivative of the damping function  DM(m)
c  DMPP(m) - The second derivative of the damping function  DM(m)
c-----------------------------------------------------------------------
      INTEGER NCMM,NCMMax,MMLR(NCMM),IDF,IDSTT,KDER,IDFF,FIRST,
     1  Lsr,m,MM,MMAX
      REAL*8 r,rhoAB,bTT(-2:2),cDS(-4:0),bDS(-4:0),aTT,br,XP,YP,
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
       DATA bDS/2.50d0,2.90d0,3.3d0,3.69d0,3.95d0/
       DATA cDS/0.468d0,0.446d0,0.423d0,0.405d0,0.390d0/
       DATA FIRST/ 1/
       SAVE FIRST, bpm, cpm
c------------------------------------------------------------------------
      IF(RHOab.LE.0) THEN
          WRITE(6,602) RHOab
          STOP
          ENDIF
      IF(IDSTT.LE.0) THEN
c===========================================
c** For Tang-Toennies type damping functions
c===========================================
          Lsr= IDF/2
          IF((IDF.LT.-4).OR.(IDF.GT.4).OR.((2*LSR).NE.IDF)) THEN
                WRITE(6,600) IDSTT,IDF
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
                  IF(KDER.GT.0) THEN
                      DMP(m)= aTT*XP*(SM(MM) - SM(MM-1))
                      IF(KDER.GT.1) DMPP(m)= -aTT*aTT*XP*(SM(MM) 
     1                                     - 2.d0*SM(MM-1) + SM(MM-2))
                      ENDIF
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
                  IF(KDER.GT.0) THEN
                      DMP(m)= aTT*XP*SM(m-1)
                      IF(KDER.GT.1)DMPP(m)= aTT*aTT*XP*(SM(m-2)-SM(m-1))
                      ENDIF
                  ENDDO
            ENDIF
          ENDIF
c
      IF(IDSTT.GT.0) THEN
c=======================================================================
c** For Douketis-Scoles-Marchetti-Zen-Thakkar type damping function ...
c=======================================================================
          IF((IDF.LT.-4).OR.(IDF.GT.0)) THEN
              WRITE(6,600) IDSTT,IDF
              STOP
              ENDIF
          IF(FIRST.EQ.1) THEN
              DO m= 1, 20
                  DO  IDFF= -4,0
                      bpm(m,IDFF)= bDS(IDFF)/DFLOAT(m)
                      cpm(m,IDFF)= cDS(IDFF)/DSQRT(DFLOAT(m))
                      ENDDO
                  ENDDO
              FIRST= 0 
              ENDIF
          br= rhoAB*r
          DO m= 1, NCMM
              MM= MMLR(m)
              XP= DEXP(-(bpm(MM,IDF) + cpm(MM,IDF)*br)*br)
              YP= 1.d0 - XP
              ZK= MM-1.d0
              DM(m)= YP**(MM-1)
c... Actually ...  DM(m)= YP**(MM + IDF/2)  :  set it up this way to 
c   avoid taking exponential of a logarithm for fractional powers (slow)
              IF(IDF.EQ.-4) THEN
                  ZK= ZK- 1.d0
                  DM(m)= DM(m)/YP
                  ENDIF
              IF(IDF.EQ.-3) THEN
                  ZK= ZK- 0.5d0
                  DM(m)= DM(m)/DSQRT(YP)
                  ENDIF
              IF(IDF.EQ.-1) THEN
                  ZK= ZK+ 0.5d0
                  DM(m)= DM(m)*DSQRT(YP)
                  ENDIF
              IF(IDF.EQ.0) THEN
                  ZK= MM
                  DM(m)= DM(m)*YP
                  ENDIF
              IF(KDER.GT.0) THEN
                  TK= bpm(MM,IDF) + 2.d0*cpm(MM,IDF)*br
                  DMP(m) = ZK*XP*rhoAB*TK*DM(m)/YP
                  IF(KDER.GT.1) THEN
c ... if desired ... calculate second derivative [for DELR case] {check this!}
                      DMPP(m)= (ZK-1.d0)*XP*TK*DMP(m)/YP
     1               - DMP(m)*TK + DMP(m)*2.d0*cpm(MM,IDF)*rhoAB**2/TK
                      ENDIF
                  ENDIF
              ENDDO   
          ENDIF  
      RETURN
  600 FORMAT(/,' *** ERROR ***  For  IDSTT=',i3,'   IDF=',i3,'  no dampi
     1ng function is defined')
  602 FORMAT( /,' ***ERROR ***  rhoAB=', F7.4,'  yields an invalid Dampi
     1ng Function definition')
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
