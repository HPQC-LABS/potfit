c=======================================================================
c** Born-Oppenheimer Breakdown & doubling function parameters.
c**                       March 16 2012
c=======================================================================
      INTEGER NUA(NSTATEMX),NUB(NSTATEMX),NTA(NSTATEMX),NTB(NSTATEMX),
     1  IFXUA(0:NBOBMX,NSTATEMX),IFXUB(0:NBOBMX,NSTATEMX),
     2  IFXTA(0:NBOBMX,NSTATEMX),IFXTB(0:NBOBMX,NSTATEMX),
     3  NwCFT(NSTATEMX),IFXwCFT(0:NBOBMX,NSTATEMX),efREF(NSTATEMX)
c
      REAL*8 UA(0:NBOBMX,NSTATEMX),UB(0:NBOBMX,NSTATEMX),
     1  TA(0:NBOBMX,NSTATEMX),TB(0:NBOBMX,NSTATEMX),
     2   wCFT(0:NBOBMX,NSTATEMX)
c
      COMMON /BLKBOB/UA,UB,TA,TB,wCFT,NUA,NUB,NTA,NTB,NwCFT,
     1  IFXUA,IFXUB,IFXTA,IFXTB,IFXwCFT,efREF
c=======================================================================
