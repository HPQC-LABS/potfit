c***********************************************************************
      double precision function Scalc(x,m,n,y,rKL,LMAX)
c** At the position 'x', Scalc is returned as the value of the m'th 
c  of the 'n' Sm(x) function defining a natural cubic spline through the
c  mesh points located at  x= y(x_i), for i=1,n.  LMAX specifies the 
c  maximum number of mesh points x= y(x_i) allowed by the calling program
c---------------------------------------------------------------------
      INTEGER  LMAX,I,K,KK,M,N
      REAL*8  x,y1,y2,y(1:LMAX),rKL(1:LMAX,1:LMAX)
      k= 0
      kk= 0
      do i=2,n
c... select interval
          if ((x.gt.y(i-1)).and.(x.le.y(i)))  k=i
          end do
      if (x.lt.y(1)) then
          k=2
          kk=1
          end if
      if (x.gt.y(n)) then
          k=n
          kk=1
          end if
      if(x.eq.y(1)) k=2
      y1=y(k-1)
      y2=y(k)
      Scalc= 0.d0
      IF(kk.eq.0) 
     1    Scalc= rKL(m,k)*((y1-x)*(((y1-x)/(y1-y2))**2-1)/6)*(y1-y2)
     2         + rKL(m,k-1)*((x-y2)*(((x-y2)/(y1-y2))**2-1)/6)*(y1-y2)
      IF(k.EQ.m) Scalc= Scalc + (y1-x)/(y1-y2)
      IF(k-1.EQ.m) Scalc= Scalc + (x-y2)/(y1-y2)
c... Asen's original coding ...
cc       Scalc=ndirac(k,m)*A(x,y1,y2)+ndirac(k-1,m)*B(x,y1,y2)+
cc   +   C(x,y1,y2)*rKL(m,k)+D(x,y1,y2)*rKL(m,k-1)
cc       else
cc       Scalc=ndirac(k,m)*A(x,y1,y2)+ndirac(k-1,m)*B(x,y1,y2)
cc     A=(x1-z)/(x1-x2)
cc     B=(z-x2)/(x1-x2)
cc     C=((x1-z)*(((x1-z)/(x1-x2))**2-1)/6)*(x1-x2)
cc     D=((z-x2)*(((z-x2)/(x1-x2))**2-1)/6)*(x1-x2)
c... Asen's original coding ...
      end
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      double precision function Sprime(x,m,n,y,rKL,LMAX)
c** At the position 'x', evaluate the derivative w.r.t. x of the m'th 
c  Sm(x) function contributing the definition of the the natural cubic
c  spline defined by function values at the  n  points  y(i) [i=1,n]
      INTEGER i,k,kk,m,n,LMAX
      REAL*8 x,del,y1,y2,y(1:LMAX),rKL(1:LMAX,1:LMAX)
      k=0
      kk=0
      do i=2,n
          if((x.gt.y(i-1)).and.(x.le.y(i)))  k=i
          enddo
      if(x.lt.y(1)) then
          k=2
          kk=1
          end if
      if (x.gt.y(n)) then
          k=n
          kk=1
          end if
      if (x.eq.y(1)) k=2
      y1=y(k-1)
      y2=y(k)
      del=y1-y2
      Sprime= 0.d0
      if(kk.eq.0) Sprime= (del-3.d0*(y1-x)**2/del)*rKL(m,k)/6.d0 +
     1                        (3.d0*(x-y2)**2/del-del)*rKL(m,k-1)/6.d0
      IF(k-1.eq.m) Sprime= Sprime + 1.d0/del 
      IF(k.eq.m) Sprime= Sprime - 1.d0/del 
ccc     if(kk.eq.0) then
ccc         Sprim=ndirac(k-1,m)/del-ndirac(k,m)/del+
ccc  +                    (del-3*(y1-x)**2/del)*rKL(m,k)/6+
ccc  +                    (3*(x-y2)**2/del-del)*rKL(m,k-1)/6
ccc       else
ccc         Sprim=ndirac(k-1,m)/del-ndirac(k,m)/del
ccc       end if
      end
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      subroutine Lkoef(n,x,A,LMAX)   
c** Call this subroutine with list of the 'n' spline x_i values in array 
c   'x' with maximum dimension 'LMAX' and it will return the LMAX x LMAX
c   array of 'rKL' coefficients used for generating the 'n' S_n(x) 
c   spline coefficient functions
c-----------------------------------------------------------------------
c   
c*** Based on nespl subroutine          
      INTEGER LMAX
      INTEGER I,J,N,INDX(1:LMAX)
      REAL*8 X(1:LMAX),A(1:LMAX,1:LMAX),B(1:LMAX,1:LMAX), d
c
      DO  i= 1,LMAX
          DO  j= 1,LMAX
              A(i,j)= 0.d0
              B(i,j)= 0.d0
              ENDDO
          ENDDO
      A(1,1)= (x(3)-x(1))/3.d0
      A(1,2)= (x(3)-x(2))/6.d0
      do i= 2,n-3
          A(i,i-1)= (x(i+1)-x(i))/6.d0
          A(i,i)= (x(i+2)-x(i))/3.d0
          A(i,i+1)= (x(i+2)-x(i+1))/6.d0
          end do
      A(n-2,n-3)= (x(n-1)-x(n-2))/6.d0
      A(n-2,n-2)= (x(n)-x(n-2))/3.d0  
      do i= 1,n-2
          B(i,i)= 1.d0/(x(i+1)-x(i))
          B(i,i+1)= -1.d0/(x(i+2)-x(i+1))-1.d0/(x(i+1)-x(i))
          B(i,i+2)= 1.d0/(x(i+2)-x(i+1))
          end do  
      call ludcmp(A,n-2,LMAX,indx,d)
      do i= 1,n 
          call lubksb(A,n-2,LMAX,indx,B(1,i))
          end do 
      do i= 1,n-2
          do j= 1,n
              A(j,i+1)= B(i,j)
              end do
          end do 
      do i= 1,n
          A(i,1)= 0.0d0
          A(i,n)= 0.0d0
          end do
      end
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE ludcmp(a,n,np,indx,d)
      INTEGER n,np,indx(n),NMAX
      double precision d,a(np,np),TINY
      PARAMETER (NMAX= 500,TINY= 1.0e-20)
      INTEGER i,imax,j,k
      double precision aamax,dum,sum,vv(NMAX)
      d= 1.d0
      do  i= 1,n
          aamax= 0.d0
          do  j= 1,n
              if (abs(a(i,j)).gt.aamax) aamax= abs(a(i,j))
              enddo
          if (aamax.eq.0.) WRITE(6,*) 'singular matrix in ludcmp'
          vv(i)= 1.d0/aamax
          enddo
      do  j= 1,n
          do  i= 1,j-1
              sum= a(i,j)
              do  k= 1,i-1
                  sum= sum-a(i,k)*a(k,j)
                  enddo
              a(i,j)= sum
              enddo
          aamax= 0.d0
          do  i= j,n
              sum= a(i,j)
              do  k= 1,j-1
                  sum= sum-a(i,k)*a(k,j)
                  enddo
              a(i,j)= sum
              dum= vv(i)*abs(sum)
              if (dum.ge.aamax) then
                  imax= i
                  aamax= dum
                  endif
              enddo
          if(j.ne.imax)then
              do  k= 1,n
                  dum= a(imax,k)
                  a(imax,k)= a(j,k)
                  a(j,k)= dum
                  enddo
              d= -d
              vv(imax)= vv(j)
              endif
          indx(j)= imax
          if(a(j,j).eq.0.)a(j,j)= TINY
              if(j.ne.n)then
                  dum= 1.d0/a(j,j)
                  do  i= j+1,n
                      a(i,j)= a(i,j)*dum
                      enddo
                  endif
          enddo
      return
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE lubksb(a,n,np,indx,b)
      INTEGER i,ii,j,ll, n,np,indx(n)
      double precision a(np,np),b(n), sum
      ii= 0
      do  i= 1,n
          ll= indx(i)
          sum= b(ll)
          b(ll)= b(i)
          if (ii.ne.0)then
              do  j= ii,i-1
                  sum= sum-a(i,j)*b(j)
                  enddo
            else if (sum.ne.0.) then
              ii= i
            endif
          b(i)= sum
          enddo
      do  i= n,1,-1
          sum= b(i)
          do  j= i+1,n
              sum= sum-a(i,j)*b(j)
              enddo
          b(i)= sum/a(i,i)
          enddo
      return
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

