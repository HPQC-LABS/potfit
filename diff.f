=============================
diff dvirdp.f tmp.f > diff.f
=============================

15,17c15,17
<        INTEGER check(0:NPARMX+3),j,counter,nParams,ISTATE,ISTART,
<      1 IDAT,i,m,k,kk,jj,ii,f
<        REAL*8  XG(8),WG(8),x(4),w(4),BTemp,BTempInv,BTOT(3),
---
>        INTEGER check,j,counter,nParams,ISTATE,ISTART,
>      1 IDAT,i,m,k,kk,jj,n
>        REAL*8  XG(8),WG(8),x(4),w(4),BTemp,BTempInv,
21,22c21
<      4 XTEMP1,dBcdP(NPARMX),dBq1dP(NPARMX),dBVIRdP(NPARMX),ZMU,etime,
<      5 elapsed(2),timestart,timestop
---
>      4 XTEMP1,dBcdP(NPARMX),dBq1dP(NPARMX),dBVIRdP(NPARMX),ZMU,dBVIR
29d27
<       data f/0/
34d31
<       save f
55,56d51
< c.. initializes the total correction term and its derivative to zero
<       BTOT = 0.d0
67,72c62
<       DO ii=0,nParams+3
<           check(ii) = 1
<           ENDDO
<       DO J= 1,nParams+3
<           INTEGRALS(J)= 0.d0
<           ENDDO
---
>       check = 1
75c65
< c.  converge to within error, or until they are repeated 7 times without
---
> c.  converge to within error, or until they are repeated n times without
77,78c67,72
<       DO  kk= 0,7
<           IF(check(0).EQ.-1) EXIT
---
>       n= 22
>       DO  kk= 1,n
>           DO J= 1,nParams+3
>               INTEGRALS(J)= 0.d0
>               ENDDO
>           IF(check.EQ.-1) EXIT
120d113
<                       IF(check(J).eq.1) THEN
122,123c115,116
<                           dBcdP(J)= -BTempInv*EXP_TERM
<      1                                       *DVtot(jj,i)*Rsq*int_fact
---
>                       dBcdP(J)= -BTempInv*EXP_TERM
>      1                                   *DVtot(jj,i)*Rsq*int_fact
125,127c118,120
<                           XTEMP1 = -VPsq*BTempInv*DVtot(jj,i)
<      1                                         + 2*dVdR(i)*dVpdP(jj,i)
<                           dBq1dP(J)= EXP_TERM*XTEMP1*Rsq*int_fact
---
>                       XTEMP1 = -VPsq*BTempInv*DVtot(jj,i)
>      1                                     + 2*dVdR(i)*dVpdP(jj,i)
>                       dBq1dP(J)= EXP_TERM*XTEMP1*Rsq*int_fact
131,133c124,125
<                           INTEGRALS(J) = INTEGRALS(J) + Const(1)
<      1                      *dBcdP(J)*WG(I) + Const(2)*dBq1dP(J)*WG(I)
<                           ENDIF
---
>                       INTEGRALS(J) = INTEGRALS(J) + Const(1)
>      1                  *dBcdP(J)*WG(I) + Const(2)*dBq1dP(J)*WG(I)
138d129
<                   IF(check(nParams+1).EQ.1) THEN
142,143d132
<                       ENDIF
<                   IF(check(nParams+2).EQ.1) THEN
147,148d135
<                       ENDIF
<                   IF(check(nParams+3).EQ.1) THEN
152d138
<                       ENDIF
155,196c141,142
< c. this section checks each integral to see if it has converged,
< c. assigning a check value of -1 to those that have, which will stop
< c. further individual calculations for them
<           IF(counter.GE.2) THEN
<               error= 0.000001
<               DO i= 1,nParams
<                   IF(integrals(i).NE.0.d0) THEN
<                       print*,i,temp(idat), integrals(i), check(i)
<                       ENDIF
<                   IF(check(i).EQ.1) THEN
<                       IF(abs(INTEGRALS(i) - dBTOTdp(i)).LE.(abs(error
<      1                                  *INTEGRALS(i)))) THEN
<                          check(i)= -1
<                          ENDIF
<                       dBTOTdp(i)= INTEGRALS(i)
<                       INTEGRALS(i)= 0.d0
<                       ENDIF
<                   ENDDO
<               DO i= 1,3
<                   IF(check(nParams+i).EQ.1) THEN
<                      IF(abs(INTEGRALS(nParams+i)-BTOT(i)).LE.(abs(error
<      1                                   *INTEGRALS(nParams+i)))) THEN
<                          check(nParams+i)= -1
<                          ENDIF
<                       BTOT(i)= INTEGRALS(nParams+i)
<                       INTEGRALS(nParams+i)= 0.d0
<                       ENDIF
<                   ENDDO
< c. this next section checks if all the integrals have converged. If they
< c. have, then check(0) is set to -1 to stop the calculations there
<               check(0)= -1
<               DO i= 1,nParams+3
<                   IF(check(i).EQ.1) THEN
<                       check(0)= 1
<                       EXIT
<                       ENDIF
<                   ENDDO
<             ELSE
<               DO i= 1,3
<                   BTOT(i)= INTEGRALS(nParams+i)
<                   ENDDO
<               DO j= 1,nParams
---
>               
>           DO j= 1,nParams
199,203c145,160
<               ENDIF
<           ENDDO
<       BVIR= Const(1)*BTOT(1) + Const(2)*BTOT(2) + Const(3)*BTOT(3)
<       DO J= 1,nParams
<           dBVIRdP(J)= dBTOTdp(J)
---
>           dBVIR=BVIR
>           BVIR= Const(1)*INTEGRALS(nParams+1) + 
>      2    Const(2)*INTEGRALS(nParams+2) +
>      3    Const(3)*Integrals(nParams+3)
> 
> c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
> c. Check to see if BVIR has converged
>           IF(counter.GE.2) THEN
>               error=0.000001
>               IF(abs(BVIR-dBVIR).LE.(abs(error*BVIR))) THEN
>                   check= -1 
>                   ENDIF
>               ENDIF 
>           DO J= 1,nParams
>               dBVIRdP(J)= dBTOTdp(J)
>               ENDDO
205d161
<       RETURN
