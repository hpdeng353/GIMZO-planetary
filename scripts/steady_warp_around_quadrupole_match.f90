program steady_warp_around_quadrupole_match

! Solves for a steady warp around an inner axisymmetric quadrupole.
! Inviscid disc, untwisted warp.
! Nonlinear theory, but assumes f^2 << 1.
! Assumes reference disc has constant H/r = epsilon and constant Sigma H.
! Matches on to a multiple of the linear solution at the inner radius.
! Uses a logarithmic radial coordinate.

! x: ln(r)

! y1: beta (tilt angle)
! y2: g (torque-related variable)

! mu1: R (characteristic radius related to inner quadrupole coefficient)
! mu2: epsilon (H/r)
! mu3: r_in (inner radius for matching to linear solution)
! mu4: r_out (outer radius)
! mu5: beta_in (quantifies the amplitude of the linear solution
!      used for matching at the inner radius in terms of the tilt angle
!      in radians that would occur at the outer radius according to
!      linear theory)

use, intrinsic :: iso_c_binding

implicit none

integer(C_INT) :: ny,nlambda,nmu
parameter (ny=2,nlambda=1,nmu=5)

integer(C_INT) :: n_spline
parameter (n_spline=6781)
real(C_DOUBLE) :: x_spline(n_spline),&
  y1_spline(n_spline),y1pp_spline(n_spline),&
  y2_spline(n_spline),y2pp_spline(n_spline),&
  y3_spline(n_spline),y3pp_spline(n_spline)
common /spline/ x_spline,y1_spline,y1pp_spline,y2_spline,y2pp_spline,&
  y3_spline,y3pp_spline

logical :: disp
integer(C_INT) :: stepmax,loop,it,itmax,n_spline_verify
real(C_DOUBLE) :: x1,x2,mu(nmu),h,hmax,tol,acc
real(C_DOUBLE) :: y1(ny),y2(ny),lambda(nlambda),mm,dmm,delta,xxx

integer :: i
open (unit=11,file='spline_gamma53.txt')
read (11, *) n_spline_verify
!open (unit=11,file='spline_gamma53',form='unformatted')
!read (11) n_spline_verify
if (n_spline_verify.ne.n_spline) then
  write (*,'(a27)') 'Incorrect spline parameters'
  stop
endif

!read (11) x_spline,y1_spline,y1pp_spline,y2_spline,y2pp_spline,y3_spline,&
!  y3pp_spline

read (11,*) (x_spline(i),i=1,n_spline)
read (11,*) (y1_spline(i),i=1,n_spline)
read (11,*) (y1pp_spline(i),i=1,n_spline)
read (11,*) (y2_spline(i),i=1,n_spline)
read (11,*) (y2pp_spline(i),i=1,n_spline)
read (11,*) (y3_spline(i),i=1,n_spline)
read (11,*) (y3pp_spline(i),i=1,n_spline)

close (11)


mu(1)=1.0d0
mu(2)=0.02d0
mu(3)=1.0d0
mu(4)=1.0d10

mu(5)=1.0d0

acc=1.0d-12
itmax=100

x1=log(mu(3))
x2=log(mu(4))
y1(1)=mu(5)*exp(-((mu(1)*mu(1))/(mu(2)*mu(3)*mu(3))))
y1(2)=(y1(1)/(dble(2)*mu(2)))

disp=.true.
!disp=.false.
h=0.001d0*(x2-x1)
hmax=h
stepmax=1000000
tol=1.0d-12

call odeint(ny,nlambda,nmu,x1,x2,y1,y2,lambda,mu,disp,h,hmax,stepmax,tol)

end program steady_warp_around_quadrupole_match

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

open (unit=11,file='results.txt')

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

close(11)

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

integer(C_INT) :: n_spline
parameter (n_spline=6781)
real(C_DOUBLE) :: x_spline(n_spline),&
  y1_spline(n_spline),y1pp_spline(n_spline),&
  y2_spline(n_spline),y2pp_spline(n_spline),&
  y3_spline(n_spline),y3pp_spline(n_spline)
common /spline/ x_spline,y1_spline,y1pp_spline,y2_spline,y2pp_spline,&
  y3_spline,y3pp_spline

real(C_DOUBLE) :: r,rr,epsilon,xxx,y1,dxxx,dy1,y1p,y1pp

r=exp(x)
rr=mu(1)
epsilon=mu(2)

! xxx means sqrt(xx*dhhdxx*dhhdxx*dhhdxx)
! y1 means sqrt(xx/dhhdxx)
! found from cubic spline interpolation

xxx=(rr/r)*sqrt(1.0d0-1.5d0*sin(y(1))*sin(y(1)))*y(2)
dxxx=-(rr/r)*1.5d0*sin(y(1))*cos(y(1))*&
  (y(3)/sqrt(1.0d0-1.5d0*sin(y(1))*sin(y(1))))*y(2)+&
  (rr/r)*sqrt(1.0d0-1.5d0*sin(y(1))*sin(y(1)))*y(4)

call splint(n_spline,x_spline,y1_spline,y1pp_spline,xxx,y1,y1p,y1pp)
dy1=y1p*dxxx

f(1)=(rr/r)*sqrt(1.0d0-1.5d0*sin(y(1))*sin(y(1)))*y1
f(2)=((rr*rr)/(epsilon*epsilon*r*r))*cos(y(1))*sin(y(1))

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
do i=1,2
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

integer(C_INT) :: n_spline
parameter (n_spline=6781)
real(C_DOUBLE) :: x_spline(n_spline),&
  y1_spline(n_spline),y1pp_spline(n_spline),&
  y2_spline(n_spline),y2pp_spline(n_spline),&
  y3_spline(n_spline),y3pp_spline(n_spline)
common /spline/ x_spline,y1_spline,y1pp_spline,y2_spline,y2pp_spline,&
  y3_spline,y3pp_spline

real(C_DOUBLE) :: r,xxx,y1,y1p,y1pp,y2,y2p,y2pp,y3,y3p,y3pp,dhhdxx,f

r=exp(x)
xxx=(mu(1)/r)*sqrt(1.0d0-1.5d0*sin(y(1))*sin(y(1)))*y(2)
call splint(n_spline,x_spline,y1_spline,y1pp_spline,xxx,y1,y1p,y1pp)
call splint(n_spline,x_spline,y2_spline,y2pp_spline,xxx,y2,y2p,y2pp)
call splint(n_spline,x_spline,y3_spline,y3pp_spline,xxx,y3,y3p,y3pp)
dhhdxx=sqrt(xxx/y1)
f=mu(2)*(y(2)/dhhdxx)
!write (*,'(4d20.12)') r,y(1),y(2),xx
write (11,'(8d20.12)') r,y(1),mu(5)*exp(-((mu(1)*mu(1))/(mu(2)*r*r))),&
  y(2),f,xxx,y2,y3


end subroutine display

!

subroutine splint(n,xa,ya,y2a,x,y,yp,ypp)

use, intrinsic :: iso_c_binding

implicit none

integer(C_INT) :: n
real(C_DOUBLE) :: xa(n),ya(n),y2a(n)
real(C_DOUBLE) :: x,y,yp,ypp

integer(C_INT) :: k,khi,klo
real(C_DOUBLE) :: a,b,h

klo=1
khi=n
1 if ((khi-klo).gt.1) then
  k=(khi+klo)/2
  if (xa(k).gt.x) then
    khi=k
  else
    klo=k
  endif
  goto 1
endif
h=xa(khi)-xa(klo)
if (h.eq.dble(0)) then
  write (*,'(a20)') 'splint: bad xa input'
  stop
endif
a=((xa(khi)-x)/h)
b=((x-xa(klo))/h)
y=a*ya(klo)+b*ya(khi)+&
  ((a*a*a-a)*y2a(klo)+(b*b*b-b)*y2a(khi))*((h*h)/6.0d0)
yp=((ya(khi)-ya(klo))/h)+&
  ((3.0d0*b*b-1.0d0)*y2a(khi)-(3.0d0*a*a-1.0d0)*y2a(klo))*(h/6.0d0)
ypp=a*y2a(klo)+b*y2a(khi)

end subroutine splint

!
