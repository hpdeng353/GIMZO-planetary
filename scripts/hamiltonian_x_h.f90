 program hamiltonian_x_h

! Solves the forced equation for nonlinear vertical oscillations in a
! warped disc:

!   H_tt + H = J^(-gamma)

!   J = H - X (cos t)^2

! Looks for a 2 pi-periodic symmetric solution as a function of the
! real variable X.

! Calculates the dimensionless Hamiltonian (negative of the averaged
! Lagrangian of the vertical oscillator) and its derivative wrt X.

! This version assumes that Z = X is real.

! The solutions are labelled by H(0) and the value of X is found by
! Newton iteration.

! x: t
! y1: H
! y2: H_t
! y3: integral for averaged Lagrangian
! y4: dy1/dlambda1
! y5: dy2/dlambda1
! y6: dy3/dlambda1

! lambda1: X

! mu1: gamma
! mu2: H(0)

use, intrinsic :: iso_c_binding

implicit none

integer(C_INT) :: ny,nlambda,nmu,nr  ! nr number of rbins
parameter (ny=6,nlambda=1,nmu=2,nr=300)

logical :: disp
integer(C_INT) :: stepmax,loop,it,itmax,loop1
real(C_DOUBLE) :: x1,x2,mu(nmu),h,hmax,tol,mm,dmm,delta,acc,twopi,&
  xx,hh,dhhdxx,offset,approx1,approx2,approx3
real(C_DOUBLE) :: y1(ny),y2(ny),lambda(nlambda),H0(nr)

integer :: i
!open (unit=11,file='H0.txt')
!read (11,*) (H0(i),i=1,8)
!close (11)
twopi=dble(8)*atan(dble(1))

mu(1)=(5.0d0/3.0d0)
!mu(1)=1.0d0
mu(2)=1.0d0

lambda(1)=0.001d0

open (unit=11,file='Hphi.txt')
do loop1=1,nr
   
do loop=1,loop1*10

mu(2)=1.0d0-dble(loop)*0.0001d0

itmax=100
acc=1.0d-10

do it=1,itmax

x1=0.0d0
x2=twopi
y1(1)=mu(2)
y1(2)=0.0d0
y1(3)=0.0d0
y1(4)=0.0d0
y1(5)=0.0d0
y1(6)=0.0d0
disp=.false.
h=0.01d0*(x2-x1)
hmax=h
stepmax=1000000
tol=1.0d-12

call odeint(ny,nlambda,nmu,x1,x2,y1,y2,lambda,mu,disp,h,hmax,stepmax,tol)

mm=y2(2)
dmm=y2(5)

delta=-(mm/dmm)

lambda(1)=lambda(1)+delta

!write (*,'(2d20.12)') lambda(1),mm

if (abs(delta).lt.acc) goto 1

enddo

write (*,'(a7)') 'Failed'
stop

1 continue

if (mu(1).eq.dble(1)) then
  offset=-0.5d0
else
  offset=-0.5d0-(dble(1)/(mu(1)-dble(1)))
endif

xx=lambda(1)
hh=(y2(3)/twopi)-offset
dhhdxx=(y2(6)/twopi)

approx1=offset+0.5d0*xx
approx2=approx1+(dble(1)/dble(8))*(mu(1)/(mu(1)+dble(1)))*xx*xx&
  +(dble(3)/dble(4))*(mu(1)/(dble(3)-mu(1)))*xx*xx
approx3=approx2+(dble(1)/dble(48))*(mu(1)/((mu(1)+dble(1))*(mu(1)+dble(1))))&
  *xx*xx*xx&
  +(dble(9)/dble(8))*(mu(1)/((dble(3)-mu(1))*(dble(3)-mu(1))))*xx*xx*xx

!write (*,'(8d20.12)') mu(1),mu(2),xx,hh,dhhdxx,&
!     approx1,approx2,approx3
enddo


!stop

x1=0.0d0
x2=twopi
y1(1)=mu(2)
y1(2)=0.0d0
y1(3)=0.0d0
y1(4)=0.0d0
y1(5)=0.0d0
y1(6)=0.0d0
disp=.true.
h=0.01d0*(x2-x1)
hmax=h
stepmax=1000000
tol=1.0d-12
write (*,'(d20.12)') mu(2)
call odeint(ny,nlambda,nmu,x1,x2,y1,y2,lambda,mu,disp,h,hmax,stepmax,tol)
enddo
close (11)
end program hamiltonian_x_h

!

subroutine odeint(ny,nlambda,nmu,x1,x2,y1,y2,lambda,mu,disp,h,hmax,stepmax,tol)

! Integrates the set of real ordinary differential equations of a
! real variable, y'(x)=f(x,y), using a 4th/5th-order embedded Runge-
! Kutta method with adaptive stepsize control (Numerical Recipes,
! Section 16.2).  The equations may contain real eigenvalues and
! real parameters.  The function f(x,y) is defined in the subroutine
! "derivs".  The real function "error" supplies a dimensionless local
! error estimate based on the local error in y.  The subroutine
! "display" can be used to display the solution during the integration.

! ny (integer): order of the system (number of functions y)

! nlambda (integer): number of real eigenvalues lambda

! nmu (integer): number of real parameters nmu

! x1 (real): initial value of x

! x2 (real): final value of x

! y1 (real array): initial values of y

! y2 (real array): final values of y (to be returned)

! lambda (real array): real eigenvalues

! mu (real array): real parameters

! disp (logical): set if the subroutine "display" is to be called
! during the integration

! h (real): initial stepsize (negative if x2<x1)

! hmax (real): maximum permitted stepsize (negative if x2<x1)

! stepmax (integer): maximum number of (successful) steps permitted

! tol (real): tolerance to be achieved

!

use, intrinsic :: iso_c_binding

implicit none

integer(C_INT) :: nymax
parameter (nymax=100)

logical :: disp
integer(C_INT) :: ny,nlambda,nmu,stepmax
real(C_DOUBLE) :: x1,x2,mu(nmu),h,hmax,tol
real(C_DOUBLE) :: y1(ny),y2(ny),lambda(nlambda)

logical :: backwards,last,good
integer(C_INT) :: i,j
real(C_DOUBLE) :: x
real(C_DOUBLE) :: y(nymax)

if (x1.gt.x2) then
  backwards=.true.
else if (x1.lt.x2) then
  backwards=.false.
else
  write (*,'(a13)') 'odeint: x1=x2'
  stop
endif

x=x1
do i=1,ny
  y(i)=y1(i)
enddo

do i=1,stepmax
  last=.false.
  if (disp) call display(ny,nlambda,nmu,x,y,lambda,mu)
  if (((.not.backwards).and.((x+h).ge.x2)).or.(backwards.and.((x+h).le.x2)))&
    then
    h=x2-x
    last=.true.
  endif
  call goodstep(ny,nlambda,nmu,x,y,lambda,mu,h,hmax,tol,good)
  if (last.and.good) then
    do j=1,ny
      y2(j)=y(j)
    enddo
    if (disp) call display(ny,nlambda,nmu,x,y,lambda,mu)
    return
  endif
enddo

write (*,'(a22)') 'odeint: too many steps'
stop

end subroutine odeint

!

subroutine goodstep(ny,nlambda,nmu,x,y,lambda,mu,h,hmax,tol,good)

!

use, intrinsic :: iso_c_binding

implicit none

integer(C_INT) :: nymax
parameter (nymax=100)

logical :: good
integer(C_INT) :: ny,nlambda,nmu
real(C_DOUBLE) :: x,mu(nmu),h,hmax,tol
real(C_DOUBLE) :: y(ny),lambda(nlambda)

integer(C_INT) :: i
real(C_DOUBLE) :: e,error,xnext
real(C_DOUBLE) :: ynext(nymax),yerror(nymax)

good=.true.
1 call rkstep(ny,nlambda,nmu,x,y,lambda,mu,h,ynext,yerror)
e=(error(ny,nlambda,nmu,x,y,lambda,mu,h,ynext,yerror)/tol)
if (e.lt.1.0d0) then
  x=x+h
  do i=1,ny
    y(i)=ynext(i)
  enddo
  h=h*dmin1(5.0d0,dmax1(0.1d0,0.9d0*(e**(-0.2d0))))
  if (dabs(h).gt.dabs(hmax)) h=hmax
  return
else
  good=.false.
  h=h*dmin1(5.0d0,dmax1(0.1d0,0.9d0*(e**(-0.2d0))))
  xnext=x+h
  if (xnext.eq.x) then
    write (*,'(a25)') 'qstep: stepsize underflow'
    stop
  endif
  goto 1
endif

end subroutine goodstep

!

subroutine rkstep(ny,nlambda,nmu,x,y,lambda,mu,h,ynext,yerror)

!

use, intrinsic :: iso_c_binding

implicit none

integer(C_INT) :: nymax
parameter (nymax=100)

integer(C_INT) :: ny,nlambda,nmu
real(C_DOUBLE) :: x,mu(nmu),h
real(C_DOUBLE) :: y(ny),lambda(nlambda),ynext(ny),yerror(ny)

integer(C_INT) :: i
real(C_DOUBLE) :: a1,a2,a3,a4,a5,a6,b21,b31,b32,b41,b42,b43,b51,b52,b53,b54,&
  b61,b62,b63,b64,b65,c1,c2,c3,c4,c5,c6,d1,d2,d3,d4,d5,d6
real(C_DOUBLE) :: ytemp(nymax),f1(nymax),f2(nymax),f3(nymax),f4(nymax),&
  f5(nymax),f6(nymax)

parameter (a1=0.0d0,&
  a2=(1.0d0/5.0d0),&
  a3=(3.0d0/10.0d0),&
  a4=(3.0d0/5.0d0),&
  a5=1.0d0,&
  a6=(7.0d0/8.0d0),&
  b21=(1.0d0/5.0d0),&
  b31=(3.0d0/40.0d0),&
  b32=(9.0d0/40.0d0),&
  b41=(3.0d0/10.0d0),&
  b42=-(9.0d0/10.0d0),&
  b43=(6.0d0/5.0d0),&
  b51=-(11.0d0/54.0d0),&
  b52=(5.0d0/2.0d0),&
  b53=-(70.0d0/27.0d0),&
  b54=(35.0d0/27.0d0),&
  b61=(1631.0d0/55296.0d0),&
  b62=(175.0d0/512.0d0),&
  b63=(575.0d0/13824.0d0),&
  b64=(44275.0d0/110592.0d0),&
  b65=(253.0d0/4096.0d0),&
  c1=(37.0d0/378.0d0),&
  c2=0.0d0,&
  c3=(250.0d0/621.0d0),&
  c4=(125.0d0/594.0d0),&
  c5=0.0d0,&
  c6=(512.0d0/1771.0d0),&
  d1=c1-(2825.0d0/27648.0d0),&
  d2=c2-0.0d0,&
  d3=c3-(18575.0d0/48384.0d0),&
  d4=c4-(13525.0d0/55296.0d0),&
  d5=c5-(277.0d0/14336.0d0),&
  d6=c6-(1.0d0/4.0d0))

call derivs(ny,nlambda,nmu,x+a1*h,y,lambda,mu,f1)
do i=1,ny
  ytemp(i)=y(i)+h*b21*f1(i)
enddo
call derivs(ny,nlambda,nmu,x+a2*h,ytemp,lambda,mu,f2)
do i=1,ny
  ytemp(i)=y(i)+h*(b31*f1(i)+b32*f2(i))
enddo
call derivs(ny,nlambda,nmu,x+a3*h,ytemp,lambda,mu,f3)
do i=1,ny
  ytemp(i)=y(i)+h*(b41*f1(i)+b42*f2(i)+b43*f3(i))
enddo
call derivs(ny,nlambda,nmu,x+a4*h,ytemp,lambda,mu,f4)
do i=1,ny
  ytemp(i)=y(i)+h*(b51*f1(i)+b52*f2(i)+b53*f3(i)+b54*f4(i))
enddo
call derivs(ny,nlambda,nmu,x+a5*h,ytemp,lambda,mu,f5)
do i=1,ny
  ytemp(i)=y(i)+h*(b61*f1(i)+b62*f2(i)+b63*f3(i)+b64*f4(i)+b65*f5(i))
enddo
call derivs(ny,nlambda,nmu,x+a6*h,ytemp,lambda,mu,f6)
do i=1,ny
  ynext(i)=y(i)+h*(c1*f1(i)+c2*f2(i)+c3*f3(i)+c4*f4(i)+c5*f5(i)+c6*f6(i))
  yerror(i)=h*(d1*f1(i)+d2*f2(i)+d3*f3(i)+d4*f4(i)+d5*f5(i)+d6*f6(i))
enddo

end subroutine rkstep

!

subroutine derivs(ny,nlambda,nmu,x,y,lambda,mu,f)

!

use, intrinsic :: iso_c_binding

implicit none

integer(C_INT) :: ny,nlambda,nmu
real(C_DOUBLE) :: x,mu(nmu)
real(C_DOUBLE) :: y(ny),lambda(nlambda),f(ny)

real(C_DOUBLE) :: j,jmg,jmgm,jmgp,dj

j=y(1)-lambda(1)*cos(x)*cos(x)
jmg=j**(-mu(1))
jmgm=j**(-mu(1)-dble(1))
jmgp=j**(-mu(1)+dble(1))
dj=y(4)-cos(x)*cos(x)

f(1)=y(2)
f(2)=-y(1)+jmg

f(4)=y(5)
f(5)=-y(4)-mu(1)*dj*jmgm

if (mu(1).eq.dble(1)) then
  f(3)=0.5d0*(y(2)*y(2)-y(1)*y(1))+log(j)
  f(6)=y(5)*y(2)-y(4)*y(1)+(dj/j)
else
  f(3)=0.5d0*(y(2)*y(2)-y(1)*y(1))-(jmgp/(mu(1)-dble(1)))
  f(6)=y(5)*y(2)-y(4)*y(1)-dj*jmg
endif

end subroutine derivs

!

function error(ny,nlambda,nmu,x,y,lambda,mu,h,ynext,yerror)

!

use, intrinsic :: iso_c_binding

implicit none

integer(C_INT) :: ny,nlambda,nmu
real(C_DOUBLE) :: error,x,mu(nmu),h
real(C_DOUBLE) :: y(ny),lambda(nlambda),ynext(ny),yerror(ny)

integer(C_INT) :: i
real(C_DOUBLE) :: ymax,yemax

! Define error estimate appropriate to application...

ymax=0.0d0
yemax=0.0d0
do i=1,3
  ymax=dmax1(ymax,dabs(y(i)))
  yemax=dmax1(yemax,dabs(yerror(i)))
enddo
!error=yemax/ymax
error=yemax

end function error

!

subroutine display(ny,nlambda,nmu,x,y,lambda,mu)

!

use, intrinsic :: iso_c_binding

implicit none

integer(C_INT) :: ny,nlambda,nmu
real(C_DOUBLE) :: x,mu(nmu)
real(C_DOUBLE) :: y(ny),lambda(nlambda)

write (11,'(2d20.12)') x, y(1)

end subroutine display

!
