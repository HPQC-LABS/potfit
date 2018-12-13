c=======================================================================
c** Parameters and count-labels for band constant (PSEL=-1) or term
c   value (PSEL=-2) fits
c
      REAL*8 TVALUE(NPARMX),ZBC(0:NVIBMX,0:NROTMX,NISTPMX,NSTATEMX),
     1 ZQC(0:NVIBMX,0:NROTMX,NISTPMX,NSTATEMX)
c
      INTEGER NSTATES,NTVALL(0:NSTATEMX),NTVI(NSTATEMX),NTVF(NSTATEMX),
     1 VMIN(NSTATEMX,NISTPMX),VMAX(NSTATEMX,NISTPMX),JTRUNC(NSTATEMX),
     2 EFSEL(NSTATEMX),NBC(0:NVIBMX,NISTPMX,NSTATEMX),
     3 NQC(0:NVIBMX,NISTPMX,NSTATEMX),
     4 BCPARI(0:NVIBMX,NISTPMX,NSTATEMX),
     5 BCPARF(0:NVIBMX,NISTPMX,NSTATEMX),
     6 QCPARI(0:NVIBMX,NISTPMX,NSTATEMX),
     7 QCPARF(0:NVIBMX,NISTPMX,NSTATEMX)
      COMMON /BLKPARAM/TVALUE,ZBC,ZQC,NSTATES,NTVALL,NTVI,NTVF,VMIN,
     1      VMAX,JTRUNC,EFSEL,NBC,NQC,BCPARI,BCPARF,QCPARI,QCPARF
c=======================================================================
