c=======================================================================
c    Block data file  BLKCOUNT.h
c=======================================================================
c** Counters for numbers of potential parameters of different types for 
c   each state
      INTEGER  TOTPOTPAR,POTPARI(NSTATEMX),POTPARF(NSTATEMX),
     1  UAPARI(NSTATEMX),UAPARF(NSTATEMX),UBPARI(NSTATEMX),
     2  UBPARF(NSTATEMX),TAPARI(NSTATEMX),TAPARF(NSTATEMX),
     3  TBPARI(NSTATEMX),TBPARF(NSTATEMX),LDPARI(NSTATEMX),
     4  LDPARF(NSTATEMX),HPARF(NSTATEMX),OSEL(NSTATEMX)
c
      COMMON /BLKCOUNT/TOTPOTPAR,POTPARI,POTPARF,UAPARI,UAPARF,UBPARI,
     1  UBPARF,TAPARI,TAPARF,TBPARI,TBPARF,LDPARI,LDPARF,HPARF,OSEL
c=======================================================================
