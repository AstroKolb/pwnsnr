subroutine init

! version: 020315

! Sod shock tube problem (a whimpy test) in 2 or 3 dimensions
! 24jan92 blondin
!=======================================================================

! GLOBALS
use global
use zone
use mpi

IMPLICIT NONE

! LOCALS


! version
CHARACTER(6), PARAMETER :: version = '020118'

! counters
INTEGER :: i, j, k, mpierr
INTEGER :: n, nm, ii

! ICs
REAL :: vt, Rsh, d_plateau, r_plateau, Lpsr, capA, alpha
REAL :: dscale, rscale, uscale, pscale, rad_max, den_max, pre_max
REAL :: rad_pwn_max, rscale_pwn, uscale_pwn, dscale_pwn, pscale_pwn
REAL, DIMENSION(2048) :: rss, dss, pss, uss
REAL, DIMENSION(400)  :: rpwn, dpwn, ppwn, upwn, gpwn

! grid values
REAL :: psr_vel, radSN, sinalpha, cosalpha, zed
REAL :: velocity, density, pressure, frac, gamma, color

! grid geometry
REAL :: xmin, xmax, ymin, ymax, zmin, zmax, zoom

! timestep calculation
REAL :: rodt, ridt, xvel, yvel, zvel, width, widthz, widthy
REAL :: rn, ph, sh

!--------------------------------------------------------------------------------
! Set up geometry and boundary conditions of grid
!
! Boundary condition flags : nleft, nright
!   = 0  :  reflecting boundary condition
!   = 1  :  inflow/outflow boundary condition
!   = 2  :  fixed inflow boundary condition
!   = 3  :  periodic
!   = 6  :  fixed outflow boundary condition
! Geometry flag : ngeom                         |  Cartesian:
!   = 0  :  planar                              |    gx = 0, gy = 0, gz = 0
!   = 1  :  cylindrical radial                  |  Cylindrical:
!   = 2  :  spherical   radial             3D= {     gx = 1, gy = 3, gz = 0
!   = 3  :  cylindrical angle                   |
!   = 4  :  spherical polar angle (theta)       |  Spherical:
!   = 5  :  spherical azimu angle (phi)         |    gx = 2, gy = 4, gz = 5

! Define the problem...

ngeomx = 2
ngeomy = 3 
ngeomz = 5

nleftx = 6
nrightx= 6
nlefty = 3
nrighty= 3
nleftz = 0
nrightz= 0

ymin   = 0.
ymax   = 2.
zmin   = 0.
zmax   = 1.

! If any dimension is angular, multiply coordinates by pi...
if(ngeomy > 2) then
   ymin = ymin * pi
   ymax = ymax * pi
endif
if(ngeomz > 2) then
   zmin = zmin * pi
   zmax = zmax * pi
endif

!======================================================================
! Set up parameters from the problem (Translational Self-Similar Driven
!				      wave with Stationary Pulsar       )

if (mype == 0) then
open (unit=15, form='formatted', file='output/' // trim(prefix) // '.icons')
   write(15,*) 'Grid    = ', imax, jmax, kmax
   write(15,*) 'version = ', version
   write(15,*) 'profile = ', trim(profile)
   write(15,*) 'Esn     = ', Esn,     ' * 1.00e+51 erg'
   write(15,*) 'Mej     = ', Mej,     ' * 2.00e+33 g'
   write(15,*) 'd0      = ', d0,      ' * 1.00e-24 g/cm**3'
   write(15,*) 'age     = ', age,     ' * 3.15e+07 s'
   write(15,*) 'Lpsr0   = ', Lpsr0,   ' * 1.00e+38 erg/s'
   write(15,*) 'v_psr   = ', v_psr,   ' * 1.00e+05 cm/s'
   write(15,*) 'pbi     = ', pbi
   write(15,*) 'psrLt   = ', psrLt,   ' * 3.15e+07 s'
   write(15,*) 'x_amb   = ', x_amb
   write(15,*) 'H_amb   = ', H_amb,   ' * 1.00e+18 cm'
   write(15,*) 'z_amb   = ', z_amb,   ' * 1.00e+18 cm'
   write(15,*) 'psi_amb = ', psi_amb, ' * pi / 4.0 rad'
   write(15,*) 'seed    = ', seed
   close(15)
end if

alpha  = 1.048  ! 1.256 for n=6?

!!! RUN SPECIFIC PARAMETER !!!

! convert input data to correct units
Esn     = Esn     * 1.e51		! E51  to  erg
Mej     = Mej     * 2.0e33		! SM   to  g
d0      = d0      * 1.e-24		!      to  g/cm**3
age     = age     * 3.15e7		! yrs  to  s
Lpsr0   = Lpsr0   * 1.e38		!      to  erg/s
v_psr   = v_psr   * 1.e5		! km/s to  cm/s
psrLT   = psrLT   * 3.15e7		! yrs  to  s
H_amb   = H_amb   * 1.e18		!      to  cm
z_amb   = z_amb   * 1.e18       !      to  cm
psi_amb = psi_amb * pi/4.		!      to  rad

! Calculate the size of the SNR and the break point in the ejecta profile (assuming nssdw = 6)
vt   = sqrt((10.*nssdw-50.)/(3.*nssdw-9.)*Esn/Mej)
capA = (5.*nssdw-25.)/(2.*pi*nssdw)*Esn/vt**5.0
Rsh  = alpha * (capA * vt**nssdw/d0)**(1.0/nssdw)*age**((nssdw-3.0)/nssdw)

d_plateau = capA / age**3
r_plateau = vt * age

!----------------------------------------------------------------------!
! read in 1D data for SSDW; scaling to above parameters
! input data is scaled to a radius of 1.0, age of 0.5, ambient density of 1.0

open(27,file='src/assets/' // trim(profile))
 do n = 1, 2048 
  read(27,*) rss(n), dss(n), pss(n), uss(n)
 enddo
close(27)

n = 2048
do while (dss(n) < 1.1)
 n = n - 1
enddo
dscale = d0
rscale = Rsh / rss(n)
uscale = (rss(1)*Rsh/rss(n)) / age / uss(1)  ! force uss(1) = rss(1)/age
pscale = dscale * uscale**2
time   = age

do n = 1, 2048
 rss(n) = rss(n) * rscale
 dss(n) = dss(n) * dscale
 pss(n) = pss(n) * pscale
 uss(n) = uss(n) * uscale
enddo

rad_max = rss(2048)
den_max = dss(2048)
pre_max = pss(2048)
assdw   = dss(1) * uss(1)**nssdw

write(*,*) assdw, dss(1), uss(1), dss(1)*uss(1)**nssdw

!----------------------------------------------------------------------!
! read in 1D data for PWN inside expanding homogeneous ejecta
! default vel1e3.dat .. various other PWN ICs available in /initpack

open(27,file='src/assets/pwn00.dat')
do n = 1, 400
 read(27,*) rpwn(n), dpwn(n), ppwn(n), upwn(n), gpwn(n)
enddo
close(27)

!----------------------------------------------------------------------!
! Scale the Input Files

! pick a radius to match up the PWN solution with the SSDW solution
Lpsr = Lpsr0 * (1. + age / psrLT) ** ( (pbi + 1.) / (1. - pbi) )   ! eqn (6) in BCF
rad_pwn_max = 1.24 * ((Esn/Mej)**3*(Lpsr/Mej)**2)**0.1 * age**1.2  ! eqn (5) in BCF


! determine scale factors
rscale_pwn = rad_pwn_max / rpwn(399)
uscale_pwn = (rad_pwn_max/age) / upwn(399)
dscale_pwn = d_plateau / dpwn(399)
pscale_pwn = dscale_pwn * uscale_pwn**2

! scale
do n = 1, 400
 rpwn(n) = rpwn(n) * rscale_pwn
 dpwn(n) = dpwn(n) * dscale_pwn
 ppwn(n) = ppwn(n) * pscale_pwn
 upwn(n) = upwn(n) * uscale_pwn
enddo


! check that simulation times are consistent in two input files
if (abs(rpwn(400)/upwn(400)/time - 1.0) > 0.1) then
  write(8,*) 'DID NOT START SIMULATION:'
  write(8,*) '   ssdw time = ', time
  write(8,*) '   pwn  time = ', rpwn(400)/upwn(400)
  stop
endif

! find radius of reverse wind shock
n = 5
do while (dpwn(n)<dpwn(n-1))
 n = n+1
enddo

psr_vel = upwn(1)
z_psr   = v_psr * time 

! grid scale
xmax   = 1.10 * rad_max + z_psr
xmin   = rpwn(n-4)

!=======================================================================
! set time and cycle counters

timep  = age
timed  = age
timem  = 0.0
ncycle = 0
ncycp  = 0
ncycd  = 0
ncycm  = 0
nfile  = 1000

! Set up grid coordinates
!zoom = exp(log(xmax/xmin)/float(imax-1))	! for logarithmic grid (slow; accurate) 
zoom = 1.00					! for linear grid

call grid(imax,xmin,xmax,zxa,zxc,zdx,zoom)
call grid(jmax,ymin,ymax,zya,zyc,zdy,1.00)
call grid(kmax,zmin,zmax,zza,zzc,zdz,1.00)

if (ndim <= 2) zzc(1) = 0.0
if (ndim == 1) zyc(1) = 0.0

!=========================================================================
! Log parameters of problem in history file

if (mype == 0) then
  write (8,*) 'Self-Similar Driven Wave'
  write (8,"(' Grid dimensions: ',i4,' x ',i4)") imax,jmax
  write (8,"(' Grid aspect ratio dy/dx = ',f4.2)") zdy(jmax/2)*zxa(imax/2)/zdx(imax/2)
  write (8,*) 
  write (8,"(' Ejecta power-law index = ',i2)") nssdw
  write (8,*) 'Constant density CSM'
  write (8,*) 
  write (8,*) 
endif

!=========================================================================
! initialize grid:

do i = 1, imax

  ! simple random number generation
  call srand(seed)			! "		"
  ph = rand()*20. + 10.			! and set these to some value
  sh = rand()*10. + 10.			!	
  ph = rand()*20. + 10.			!

  do j = 1, js

    ! mapping function for the offset SNR
    radSN    = sqrt(zxc(i)**2 + z_psr**2 - 2.0*zxc(i)*z_psr*cos(zyc(mypey*js+j)))
    sinalpha = sin(zyc(mypey*js+j)) * z_psr/radSN
    cosalpha = sqrt(1.0-sinalpha**2)
    zed      = zxc(i) * cos(zyc(mypey*js+j) - psi_amb) - z_psr*cos(psi_amb) - z_amb

    ! flip cosalpha if necessary
    if ( zyc(mypey*js+j) > pi/2. .and. zyc(mypey*js+j) < 3*pi/2. ) then
!	do nothing
    elseif ( zxc(i) < z_psr*cos(zyc(mypey*js+j))) then
	cosalpha = - cosalpha
    endif

    if (radSN < rss(1)) then		! plateau
	velocity = radSN / time
	density  = min(d_plateau, assdw / (velocity**nssdw))
	pressure = pss(1)
	color = 0.0
    elseif (radSN > rad_max) then	! ISM
	velocity = 0.0
	density  = d0 * (1.0-(2.0-x_amb)/x_amb * tanh(zed/H_amb))
	pressure = pre_max
	color = 1.0
    else				! SNR
	n = 2
	do while (rss(n) < radSN)
	  n = n + 1
	enddo
	nm = n - 1
	frac = (radSN-rss(nm))/(rss(n)-rss(nm))
	density  = frac*dss(n) + (1.0-frac)*dss(nm)
	pressure = frac*pss(n) + (1.0-frac)*pss(nm)
	velocity = frac*uss(n) + (1.0-frac)*uss(nm)
	color = 0.0

	! density perturbation
	if (radSN < 0.7*rad_max .and. radSN > 0.5*rad_max) then
	  density = density + (rand()*0.2+0.1)*density*sin(ph * zyc(mypey*js+j))
	  density = density + (rand()*0.1+0.1)*density*sin(sh * zyc(mypey*js+j))
	end if

    endif

    do k = 1, ks
	zro(i,j,k) = density
	zpr(i,j,k) = pressure
	zux(i,j,k) = velocity * cosalpha + v_psr * cos(zyc(mypey*js+j))
	zuy(i,j,k) = velocity * sinalpha - v_psr * sin(zyc(mypey*js+j))
	zuz(i,j,k) = 0.0
	zfl(i,j,k) = 0.0
	zgm(i,j,k) = 5./3.
	zga(i,j,k) = age
	zcl(i,j,k) = color
    enddo
  enddo
enddo


! Add the PWN:
i = 1
do while(zxc(i) < rad_pwn_max)
  n = 2
  do while(rpwn(n) < zxc(i))
    n = n + 1
  enddo
  nm = n - 1
  frac = (zxc(i)-rpwn(nm))/(rpwn(n)-rpwn(nm))
  density  = frac*dpwn(n) + (1.0-frac)*dpwn(nm)
  pressure = frac*ppwn(n) + (1.0-frac)*ppwn(nm)
  velocity = frac*upwn(n) + (1.0-frac)*upwn(nm)
  gamma    = frac*gpwn(n) + (1.0-frac)*gpwn(nm)

  do j = 1, js
    do k = 1, ks
	zro(i,j,k) = density
	zpr(i,j,k) = pressure
	zux(i,j,k) = velocity
	zuy(i,j,k) = 0.0
	zuz(i,j,k) = 0.0
	zfl(i,j,k) = 0.0
	zgm(i,j,k) = gamma
	zga(i,j,k) = age
	zcl(i,j,k) = 0.0
    enddo
  enddo
  i = i + 1
enddo


! Boundary Condition
dinflo = Lpsr
pinflo = (0.01*psr_vel)**2
uinflo = psr_vel
vinflo = 0.0
winflo = 0.0
ginflo = 4./3.
cinflo = 0.0

potflo = zpr(imax,js/2,1)
gotflo = 5./3.
cotflo = 1.0

!########################################################################
! Compute Courant-limited timestep

ridt = 0.

if (ndim == 2) then

 do j = 1, js
  do i = 1, imax
    widthy = zdy(j+mypey*js)
    if(ngeomy > 2) widthy = widthy*zxc(i)
    width  = min(zdx(i),widthy)
    svel = sqrt(gam*zpr(i,j,1)/zro(i,j,1))/width
    xvel = abs(zux(i,j,1)) / zdx(i)
    yvel = abs(zuy(i,j,1)) / widthy
    ridt = max(xvel,yvel,svel,ridt)
  enddo
 enddo

else

 do k = 1, ks
  do j = 1, js
   do i = 1, imax
     widthy = zdy(j+mypey*js)
     widthz = zdz(k+mypez*ks)
     if(ngeomy.gt.2) widthy = widthy*zxc(i)
     if(ngeomz.gt.2) widthz = widthz*zxc(i)
     if(ngeomz.eq.5) widthz = widthz*sin(zyc(j+mypey*js))
     width  = min(zdx(i),widthy,widthz)
     svel = sqrt(gam*zpr(i,j,k)/zro(i,j,k))/width
     xvel = abs(zux(i,j,k)) / zdx(i)
     yvel = abs(zuy(i,j,k)) / widthy
     zvel = abs(zuz(i,j,k)) / widthz
     ridt = max(xvel,yvel,zvel,svel,ridt)
   enddo
  enddo
 enddo

endif

call MPI_ALLREDUCE( ridt, rodt, 1, VH1_DATATYPE, MPI_MAX, MPI_COMM_WORLD, mpierr )
dt = courant / rodt

return
end

!#####################################################################

subroutine grid( nzones, xmin, xmax, xa, xc, dx, zoom )

! Create grid to cover physical size from xmin to xmax
! number of physical grid zones is nzones
!
! xa(1) is left boundary location - at xmin
! xa(nzones+1) is right boundary location - at xmax
!----------------------------------------------------------------------
IMPLICIT NONE

! LOCALS
integer :: nzones, n
real, dimension(nzones) :: xa, dx, xc
real :: dxfac, xmin, xmax, zoom

!=======================================================================

if (zoom==1.0d0) then
  dxfac = (xmax - xmin) / float(nzones)
  do n = 1, nzones
    xa(n) = xmin + (n-1)*dxfac
    dx(n) = dxfac
    xc(n) = xa(n) + 0.5*dx(n)
  enddo
else
  dx(1) = (xmax-xmin)*(zoom-1.0)/(zoom**nzones - 1.0)
  xa(1) = xmin
  xc(1) = xa(1) + 0.5d0*dx(1)
  do n = 2, nzones
    xa(n) = xa(n-1) + dx(n-1)
    dx(n) = dx(n-1) * zoom
    xc(n) = xa(n) + 0.5d0*dx(n)
  enddo
endif

return
end
