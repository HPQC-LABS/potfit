c=======================================================================
c** Parameters and count-labels for band constant (PSEL=-1) or term
c   value (PSEL=-2) fits
c
      REAL*8 TVALUE(NPARMX),ZBC(0:NVIBMX,0:NROTMX,NISTPMX,NSTATEMX)
c
      INTEGER NSTATES,NTVALL(0:NSTATEMX),NTVI(NSTATEMX),NTVF(NSTATEMX),
     1 VMIN(NSTATEMX,NISTPMX),VMAX(NSTATEMX,NISTPMX),
     2 NBC(0:NVIBMX,NISTPMX,NSTATEMX),BCPARI(0:NVIBMX,NISTPMX,NSTATEMX),
     3 BCPARF(0:NVIBMX,NISTPMX,NSTATEMX)
      COMMON /BLKPARAM/TVALUE,ZBC,NSTATES,NTVALL,NTVI,NTVF,VMIN,VMAX,
     1       NBC,BCPARI,BCPARF
c=======================================================================
