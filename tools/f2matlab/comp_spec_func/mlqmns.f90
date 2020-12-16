PROGRAM mlqmns
 
! Code converted using TO_F90 by Alan Miller
! Date: 2004-06-28  Time: 12:58:12

!       =========================================================
!       Purpose: This program computes the associated Legendre
!                functions Qmn(x) and their derivatives Qmn'(x)
!                for a given order using subroutine LQMNS
!       Input :  x --- Argument of Qmn(x)
!                m --- Order of Qmn(x),  m = 0,1,2,...
!                n --- Degree of Qmn(x), n = 0,1,2,...
!       Output:  QM(n) --- Qmn(x)
!                QD(n) --- Qmn'(x)
!       Examples:
!                m = 1,  N = 5,  x = .5
!                n        Qmn(x)           Qmn'(x)
!               -------------------------------------
!                0    .11547005D+01    .76980036D+00
!                1    .10530633D+01    .23771592D+01
!                2   -.72980606D+00    .51853281D+01
!                3   -.24918526D+01    .10914062D+01
!                4   -.19340866D+01   -.11454786D+02
!                5    .93896830D+00   -.18602587D+02

!                m = 2,  N = 5,  x = 2.5
!                n        Qmn(x)           Qmn'(x)
!               -------------------------------------
!                0    .95238095D+00   -.52607710D+00
!                1    .38095238D+00   -.36281179D+00
!                2    .12485160D+00   -.17134314D+00
!                3    .36835513D-01   -.66284127D-01
!                4    .10181730D-01   -.22703958D-01
!                5    .26919481D-02   -.71662396D-02
!       =========================================================

IMPLICIT DOUBLE PRECISION (q,x,y)
DIMENSION qm(0:200),qd(0:200)
WRITE(*,*)'Please enter m, N, and x '
!        READ(*,*)M,N,X
m=1
n=5
x=.5
WRITE(*,30)m,n,x
CALL lqmns(m,n,x,qm,qd)
WRITE(*,*)
WRITE(*,*)'  n        Qmn(x)           Qmn''(X)'
WRITE(*,*)' -------------------------------------'
DO  j=0,n
  WRITE(*,20)j,qm(j),qd(j)
END DO
20      FORMAT(1X,i3,2D17.8)
30      FORMAT(1X,'m =',i2,',  ','n =',i2,',  ','x =',f5.1)
END PROGRAM mlqmns


SUBROUTINE lqmns(m,n,x,qm,qd)

!       ========================================================
!       Purpose: Compute associated Legendre functions Qmn(x)
!                and Qmn'(x) for a given order
!       Input :  x --- Argument of Qmn(x)
!                m --- Order of Qmn(x),  m = 0,1,2,...
!                n --- Degree of Qmn(x), n = 0,1,2,...
!       Output:  QM(n) --- Qmn(x)
!                QD(n) --- Qmn'(x)
!       ========================================================


INTEGER, INTENT(IN)                      :: m
INTEGER, INTENT(IN)                      :: n
DOUBLE PRECISION, INTENT(IN)             :: x
DOUBLE PRECISION, INTENT(OUT)            :: qm(0:n)
DOUBLE PRECISION, INTENT(OUT)            :: qd(0:n)
IMPLICIT DOUBLE PRECISION (a-h,o-z)


DO  k=0,n
  qm(k)=0.0D0
  qd(k)=0.0D0
END DO
IF (DABS(x) == 1.0D0) THEN
  DO  k=0,n
    qm(k)=1.0D+300
    qd(k)=1.0D+300
  END DO
  RETURN
END IF
ls=1
IF (DABS(x) > 1.0D0) ls=-1
xq=DSQRT(ls*(1.0D0-x*x))
q0=0.5D0*DLOG(DABS((x+1.0)/(x-1.0)))
q00=q0
q10=-1.0D0/xq
q01=x*q0-1.0D0
q11=-ls*xq*(q0+x/(1.0D0-x*x))
qf0=q00
qf1=q10
DO  k=2,m
  qm0=-2.0D0*(k-1.0)/xq*x*qf1-ls*(k-1.0)*(2.0-k)*qf0
  qf0=qf1
  qf1=qm0
END DO
IF (m == 0) qm0=q00
IF (m == 1) qm0=q10
qm(0)=qm0
IF (DABS(x) < 1.0001D0) THEN
  IF (m == 0.AND.n > 0) THEN
    qf0=q00
    qf1=q01
    DO  k=2,n
      qf2=((2.0*k-1.0D0)*x*qf1-(k-1.0)*qf0)/k
      qm(k)=qf2
      qf0=qf1
      qf1=qf2
    END DO
  END IF
  qg0=q01
  qg1=q11
  DO  k=2,m
    qm1=-2.0D0*(k-1.0)/xq*x*qg1-ls*k*(3.0-k)*qg0
    qg0=qg1
    qg1=qm1
  END DO
  IF (m == 0) qm1=q01
  IF (m == 1) qm1=q11
  qm(1)=qm1
  IF (m == 1.AND.n > 1) THEN
    qh0=q10
    qh1=q11
    DO  k=2,n
      qh2=((2.0*k-1.0D0)*x*qh1-k*qh0)/(k-1.0)
      qm(k)=qh2
      qh0=qh1
      qh1=qh2
    END DO
  ELSE IF (m >= 2) THEN
    qg0=q00
    qg1=q01
    qh0=q10
    qh1=q11
    DO  l=2,n
      q0l=((2.0D0*l-1.0D0)*x*qg1-(l-1.0D0)*qg0)/l
      q1l=((2.0*l-1.0D0)*x*qh1-l*qh0)/(l-1.0D0)
      qf0=q0l
      qf1=q1l
      DO  k=2,m
        qmk=-2.0D0*(k-1.0)/xq*x*qf1-ls*(k+l-1.0)* (l+2.0-k)*qf0
        qf0=qf1
        qf1=qmk
      END DO
      qm(l)=qmk
      qg0=qg1
      qg1=q0l
      qh0=qh1
      qh1=q1l
    END DO
  END IF
ELSE
  IF (DABS(x) > 1.1) THEN
    km=40+m+n
  ELSE
    km=(40+m+n)*INT(-1.0-1.8*LOG(x-1.0))
  END IF
  qf2=0.0D0
  qf1=1.0D0
  DO  k=km,0,-1
    qf0=((2.0*k+3.0D0)*x*qf1-(k+2.0-m)*qf2)/(k+m+1.0)
    IF (k <= n) qm(k)=qf0
    qf2=qf1
    qf1=qf0
  END DO
  DO  k=0,n
    qm(k)=qm(k)*qm0/qf0
  END DO
END IF
IF (DABS(x) < 1.0D0) THEN
  DO  k=0,n
    qm(k)=(-1)**m*qm(k)
  END DO
END IF
qd(0)=((1.0D0-m)*qm(1)-x*qm(0))/(x*x-1.0)
DO  k=1,n
  qd(k)=(k*x*qm(k)-(k+m)*qm(k-1))/(x*x-1.0)
END DO
RETURN
END SUBROUTINE lqmns
