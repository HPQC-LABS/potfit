      PROGRAM testretardation
      INTEGER I
      REAL*8 r,V,C6,C6ulr,dVdR,dVdalphaD,A(1000000),B(1000000)
      C6= 6.7185d6
      WRITE(*,'(/"zone V")')
      DO I=1,1000000
          r = I*0.0001d0
          CALL RETARDATION(r,C6,V,dVdR,dVdalphaD)
          IF(ceiling(r*100).eq.r*100)
     1    WRITE(*,'(F6.2,1P,E18.8,0P)') r,V
          A(I) = dVdR
          B(I) = V
          ENDDO
      WRITE(*,'(/"zone dV/dr")')
      DO I = 1,1000000
          r = I*0.0001d0
          IF(ceiling(r*100).eq.r*100)
     1    WRITE(*,'(F6.2,1P,E18.8,0P)') r,A(I)
          ENDDO
      WRITE(*,'(/"zone approx. dV/dr")')
      DO I = 1,999999
          r = I*0.0001d0
          IF(ceiling(r*100).eq.r*100)
     1    WRITE(*,'(F6.2,1P,E18.8,0P)') r,(B(I+1)-B(I))/(0.0001d0)
          ENDDO
      WRITE(*,'(/"zone no retardation")')
      DO I=1,1000000
          r = I*0.0001d0
          C6ulr = -C6/r**6
          IF(ceiling(r*100).eq.r*100)
     1    WRITE(*,'(F6.2,1P,E18.8,0P)') r,C6ulr
          ENDDO 
      ENDPROGRAM

c*********************************************************************** 
      SUBROUTINE RETARDATION(rval,C6,V,dVdR,dVdC6)
c*********************************************************************** 
c...  Subroutine to determine retardation replacement to C6/r^6 term
c...  in ULR(r) function. 
c...  INPUT= r(angstroms),alphaD(a.u.)
c...  OUTPUT= V(r)(cm-1) (replacment for C6 term),dVdR,dVdalphaD
c----------------------------------------------------------------------- 
      INTEGER K,first 
      REAL*8 alphaFS,alphaDomega(50),omega(50),w(50),r,rval,E,
     1       V,dVdR,dVdC6,numerator,denominator,dnum,C6,OaF,OaF2,
     2       R1,R2,R3,R4,wAdO2 
c----------------------------------------------------------------------- 
      alphaFS= 7.29735257d-3
      E= 43486.5564d0*4.5563d-6  !convert from cm-1 to hartree
      r=rval*1.88972613392d0     !convert from angstroms to bohr radii
c-----------------------------------------------------------------------
      OPEN(UNIT= 3,STATUS='OLD',FILE='retardation.3')
      IF(FIRST.EQ.0) THEN
          DO  K=1,50
              READ(3,*) w(K),omega(K),alphaDomega(K)
              ENDDO
              FIRST=1
          ENDIF
c-----------------------------------------------------------------------
      numerator = 0.d0
      denominator = 0.d0
      dnum = 0.d0
      R1 = 1.d0/r
      R2 = R1*R1 
      R3 = R1*R2
      R4 = R2*R2
      DO  K=1,50
          OaF = omega(k)*alphaFS
          OaF2 = OaF*OaF
          wAdO2 = w(K)*alphaDomega(K)**2
          numerator = numerator + wAdO2
     1    *dexp(-2.d0*OaF*r)*((OaF*r)**4 + 2.d0*(OaF*r)**3 +
     2    5.d0*(OaF*r)**2 + 6.d0*(OaF*r)+3.d0)
          dnum = dnum + wAdO2*dexp(-2.d0*OaF*r)*(2.d0*OaF*R2*(OaF2*OaF2 
     3           +2.d0*(OaF2*OaF)*R1+5.d0*(OaF2)*R2+6.d0*(OaF)*R3
     6           +3.d0*R4)+R3*(2.d0*(OaF2*OaF2)+6.d0*(OaF2*OaF)*R1
     9           +20.d0*(OaF2)*R2+30.d0*(OaF)*R3+18.d0*R4))
          denominator = denominator + wAdO2
          ENDDO
      dVdC6 = -numerator/(denominator*9.d0*rval**6) ! num/den is dimensionless
      V = C6*dVdC6   
      dVdR = C6*dnum/(denominator*9.d0)*(1.88972613392d0)**7
c-----------------------------------------------------------------------
c*********************************************************************** 
      END
c*********************************************************************** 
