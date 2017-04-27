c=======================================================================
c** Partial derivative arrays for fits and uncertainties (fununc)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      REAL*8 DVtot(HPARMX,NPNTMX),DLDDRe(NPNTMX,NSTATEMX),
     1  DUADRe(NPNTMX,NSTATEMX),DUBDRe(NPNTMX,NSTATEMX),
     2  DTADRe(NPNTMX,NSTATEMX),DTBDRe(NPNTMX,NSTATEMX),
     3  DBDB(0:NbetaMX,NPNTMX,NSTATEMX),DBDRe(NPNTMX,NSTATEMX),
     4  dVpdP(HPARMX,NPNTMX)
      COMMON/BLKDVDP/DVtot,DUADRe,DUBDRe,DTADRe,DTBDRe,DLDDRe,DBDB,
     1 DBDRe,dVpdP
c=======================================================================
