program post
! post.f90 pwnsnr post-processing script
! cekolb - 12:15:14

use NETCDF

! LOCALS

INTEGER :: n, i, j, k, ndim, nvar, nv
INTEGER :: imax, jmax, kmax


CHARACTER(LEN=1) :: tmp1
CHARACTER(LEN=4) :: tmp4, tmp2
CHARACTER(LEN=80) :: varname
CHARACTER(LEN=50) :: prefix
CHARACTER(LEN=50) :: filename
CHARACTER(LEN=50) :: ghostfile
CHARACTER(LEN=50) :: geofile 
CHARACTER(LEN=50) :: casefile
CHARACTER(LEN=50) :: varfile

CHARACTER(LEN=80) :: label1a, label1b, label1c, label1d, label1e, label1f, label1g, label1h, label1i
CHARACTER(LEN=80) :: label2, label3, label4, coord1, coord2, coord3

INTEGER :: ipart, nsteps
INTEGER :: ncid, ncstat, ncid_slabs, rank

INTEGER :: n, i, j, k, ndim, nvar, nv
INTEGER :: imax, jmax, kmax
INTEGER, DIMENSION(4) :: start
INTEGER, DIMENSION(:,:), ALLOCATABLE :: ghostcells
INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: gs

REAL, DIMENSION(:,:,:), ALLOCATABLE :: var, slab
REAL, DIMENSION(:), ALLOCATABLE :: zxc, zyc, zzc
REAL, DIMENSION(0:500) :: time

