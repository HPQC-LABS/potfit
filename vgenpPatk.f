c***********************************************************************
      SUBROUTINE VGENP(ISTATE,RDIST,VDIST,dVdR,d2VdR2,IDAT)
c***********************************************************************
C Fit of the Ar-Ar potential
C K. Patkowski and K. Szalewicz, to be submitted (2010).
C
C Input:
C r - distance in Angstrom
C
C Output:
C value  - interaction potential in cm-1
C deriv  - analytic first derivative of potential in cm-1/bohr
c***********************************************************************
      INCLUDE 'arrsizes.h'
      INTEGER I,IDAT,ISTATE
      REAL*8 RDIST(8),VDIST(8),BETADIST,dVdR(8),d2VdR2(8),ar2pot,ar2der
c
      BETADIST= 0.d0
      DO I= 1,8
          VDIST(I)= ar2pot(RDIST(I))
          d2VdR2(I)= 0.0
          dVdR(I)=  ar2der(RDIST(I))
          ENDDO
cc    write(6,*) RDIST,'    ',VDIST,'     ',dRdV
      return
      end
c***********************************************************************
c***********************************************************************

c*********************************************************************
      real*8 function ar2pot(r)
      implicit real*8 (a-h,o-z)
      dimension cn(6)
      data cn /64.288984d0,1514.86211d0,50240.0d0,1898195.d0,
     $  86445426.d0,4619452502.d0/
      data ifirst /1/
      data coef /127641878.945519894d0/
      data coef2/-26138949.621478189d0/
      data coef3/-115672346.174201056d0/
      data coef4/2064381.526204719d0/
      data coef5/-58371.409016267d0/
      data alpha /1.553386357296d0/
      data ddd /2.393847610341d0/
      data conv /219474.63d0/
      data rjoin /1.3d0/
      data e0/324.d0/
      data a0/-84134425.1811d0/
      data b0/33397004.2848d0/
      data c0/-4443000.0d0/
      data alpha0/2.3577d0/
      data beta0/-1.2756d0/
      save
c
      nc=6
      if(ifirst.eq.1) then
      do nn=1,nc
        cn(nn) = cn(nn)*conv
      enddo
      e0=e0*conv
      ifirst = 0
      endif
      if (r.ge.rjoin) then
        rpt = r/0.529177209d0
        term = (coef+coef2*rpt+coef3/rpt+coef4*rpt*rpt
     $         +coef5*rpt*rpt*rpt)*dexp(-alpha*rpt)
        sumc=0.0d0
        do nn=1,nc
          sumc=sumc+d(2*nn+4,ddd,rpt)*cn(nn)/rpt**(2*nn+4)
        enddo
        asy = sumc
        term = term - asy 
        ar2pot = term
      else
        rpt = r/0.529177209d0
        ar2pot=dexp(-alpha0*rpt-beta0*rpt*rpt)
        ar2pot=ar2pot*(e0/rpt+a0+b0*rpt+c0*rpt*rpt)
      end if
      return
      end
cc
      real*8 function ar2der(r)
      implicit real*8 (a-h,o-z)
      dimension cn(6)
      data cn /64.288984d0,1514.86211d0,50240.0d0,1898195.d0,
     $  86445426.d0,4619452502.d0/
      data ifirst /1/
      data coef /127641878.945519894d0/
      data coef2/-26138949.621478189d0/
      data coef3/-115672346.174201056d0/
      data coef4/2064381.526204719d0/
      data coef5/-58371.409016267d0/
      data alpha /1.553386357296d0/
      data ddd /2.393847610341d0/
      data conv /219474.63d0/
      data rjoin /1.3d0/
      data e0/324.d0/
      data a0/-84134425.1811d0/
      data b0/33397004.2848d0/
      data c0/-4443000.0d0/
      data alpha0/2.3577d0/
      data beta0/-1.2756d0/
      save
c
      nc=6
      if(ifirst.eq.1) then
      do nn=1,nc
        cn(nn) = cn(nn)*conv
      enddo
      e0=e0*conv
      ifirst = 0
      endif
      if (r.ge.rjoin) then
        rpt = r/0.529177209d0
        expfac = exp(-alpha*rpt)
        dercoe = -alpha
        der= coef * dercoe * expfac
        der = der + coef2 * (1.0 + rpt*dercoe) * expfac
        der = der + coef3 * (-1.0/rpt**2 + dercoe/rpt)*expfac
        der = der + coef4 * (2.d0*rpt + rpt*rpt*dercoe) * expfac
        der = der + coef5 * (3.d0*rpt*rpt + rpt*rpt*rpt*dercoe) * expfac
        sumcd=0.0d0
        do nn=1,nc
            nnn=2*nn+4
            sumcd = sumcd + cn(nn)*(-nnn/rpt**(nnn+1)
     &                    + ddd*exp(-ddd*rpt)*d1(nnn,ddd,rpt)/rpt**nnn
     &                    - exp(-ddd*rpt)*da(nnn,ddd,rpt) )
        enddo
        der = der - sumcd
        ar2der = der
      else
        rpt = r/0.529177209d0
        e1=dexp(-alpha0*rpt-beta0*rpt*rpt)
        de1=-alpha0-2.d0*beta0*rpt
        e2=(e0/rpt+a0+b0*rpt+c0*rpt*rpt)
        de2=-e0/(rpt*rpt)+b0+2.0*c0*rpt
        ar2der=e1*de2+e2*de1*e1
      end if
      return
      end
cc
      double precision function d(n,beta,r)
c
c     calculate the damping factor
c
      implicit real*8 (a-h,o-z)
      br=beta*r
      sum=1.0d0
      term=1.0d0
      ncn=n
      do i=1,ncn
        term=term*br/i
        sum=sum+term
      enddo
      d=1.0d0 - dexp(-br)*sum
c     write(6,*) n,beta,r,d
      return
      end
cc
      double precision function d1(n,beta,r)
c
c     calculate part of the derivative of the damping factor
c
      implicit real*8 (a-h,o-z)
      br=beta*r
      sum=1.0d0
      term=1.0d0
      ncn=n
      do i=1,ncn
        term=term*br/i
        sum=sum+term
      enddo
      d1= sum
c     write(6,*) n,beta,r,d
      return
      end
cc
      double precision function da(n,beta,r)
c
c     calculate part of the derivative of the damping factor
c
      implicit real*8 (a-h,o-z)
      br=beta*r
      sum=-n/r**(n+1)
      term=1.0d0
      ncn=n
      do i=1,ncn-1
        term=term*br/i
        sum=sum+term*(i-n)/r**(n+1)
      enddo
      da = sum
c     write(6,*) n,beta,r,d
      return
      end

