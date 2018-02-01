subroutine sweepx1(xexpand)

! This subroutine performs 1D hydro sweeps in the X direction, 
! looping over js by ks rows.  At the end, data is packed into
! a SEND buffer for data transpose in sweepy.
!-----------------------------------------------------------------------

! GLOBALS
use zone
use global
use sweeps

IMPLICIT NONE

! LOCALS
INTEGER :: i, j, k, n, m
REAL :: xexpand, rad

INTEGER :: nleftalt	! 0 = basic reflection

!-----------------------------------------------------------------------

sweep  = 'x'
ngeom  = ngeomx
nleft  = nleftx
nright = nrightx
nmin   = 7
nmax   = imax + 6

nleftalt = 7	! Define alternate boundary conditition here

! xexpand - expand grid (and contract data)
do i = 1, imax
  n = i + 6
  xa0(n) = zxa(i)
  dx0(n) = zdx(i)
enddo
xa0(nmax+1) = xa0(nmax) + dx0(nmax)
if (xexpand .gt. 0.0) then
  do n = nmin, nmax+1
    xa0(n) = xa0(n) * (1. + xexpand)
  enddo
  do n = nmin, nmax
    dx0(n) = xa0(n+1) - xa0(n)
  enddo
endif

! Now Loop over each row...

do k = 1, ks
 do j = 1, js
   theta = zyc(mypey*js+j)

   ! Put state variables into 1D arrays, padding with 6 ghost zones
   nleft = nleftx
   do i = 1,imax
     n = i + 6
     r  (n) = zro(i,j,k)
     p  (n) = zpr(i,j,k)
     u  (n) = zux(i,j,k)
     v  (n) = zuy(i,j,k)
     w  (n) = zuz(i,j,k)
     f  (n) = zfl(i,j,k)
     g  (n) = zgm(i,j,k)
     ga (n) = zga(i,j,k) + dt
     c  (n) = zcl(i,j,k)

!     if ( (i < 3) .and. (u(n) .le. 0.0) ) nleft = nleftalt

!     xa0(n) = zxa(i)
!     dx0(n) = zdx(i)
     xa (n) = zxa(i)
     dx (n) = zdx(i)
     p  (n) = max(smallp,p(n))
     e  (n) = p(n)/(r(n)*(g(n)-1.))+0.5*(u(n)**2+v(n)**2+w(n)**2)

     rad = xa(n) + 0.5*dx(n)
     if ( (i < 3) .and. (u(n) .le. - dinflo / r(n) / (2.*pi*(rad*uinflo)**2)) ) nleft = nleftalt

   enddo

   ! Do 1D hydro update using PPMLR
   call ppmlr 

   ! Put updated values into send buffer for sweepy, dropping ghost zones
   do i = 1, imax
     n = i + 6
     send1(1,k,j,i) = r(n)
     send1(2,k,j,i) = p(n)
     send1(3,k,j,i) = u(n)
     send1(4,k,j,i) = v(n)
     send1(5,k,j,i) = w(n)
     send1(6,k,j,i) = f(n)
     send1(7,k,j,i) = g(n)
     send1(8,k,j,i) = ga(n)
     send1(9,k,j,i) = c(n)
   enddo

 enddo
enddo

do i = 1, imax
  n = i + 6
  zxa(i) = xa0(n)
  zdx(i) = dx0(n)
  zxc(i) = xa0(n) + 0.5*dx0(n)
enddo

return
end

