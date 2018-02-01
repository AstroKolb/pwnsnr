subroutine prin_ensight

! collect all rows of data onto mypey=0, and then 
! write out pez netcdf data files containing all variables (plus time as an attribute).
! These binary data files can then be merged later using mergeslabs.f90.

! output now prints extra fluid variable 'color' 
!---------------------------------------------------------------------------

use global
use zone
use mpi

IMPLICIT NONE

! LOCALS
LOGICAL :: periodic_in_y, periodic_in_z

CHARACTER(LEN=1) :: char
CHARACTER(LEN=4) :: tmp1, tmp2
CHARACTER(LEN=50) :: filename, casefile, geomfile, viewvh1, testfile, zrofile, zprfile, zuxfile, zuyfile, zuzfile
CHARACTER(LEN=80) :: label1a, label1b, label1c, label1d, label1e, label1f, label1g, label1h
CHARACTER(LEN=80) :: label2, label3, label4 

INTEGER, PARAMETER :: nvarout = 8   ! # of arrays to write into netcdf file
INTEGER :: gathbuffer_size, m, jsk, nv, jjmax, kkmax, ipart, i, j, k, mpierr

REAL(4) :: sx(imax), sy(jmax), sz(kmax)
REAL(4), DIMENSION(imax,nvarout,jmax/pey,kmax/pez) :: send_buff
REAL(4), DIMENSION(imax,nvarout,jmax/pey,kmax/pez,pey) :: recv_buff

REAL(4), ALLOCATABLE, DIMENSION(:,:,:) :: var

!------------------------------------------------------------------------------

! everybody loads up a send buffer; data is gathered on mypey=0 procs.

do k = 1, ks
 do j = 1, js
  do i = 1, imax
    send_buff(i,1,j,k) = zro(i,j,k)
    send_buff(i,2,j,k) = zpr(i,j,k)
    send_buff(i,3,j,k) = zgm(i,j,k)
    send_buff(i,4,j,k) = zga(i,j,k)
    send_buff(i,5,j,k) = zcl(i,j,k)
    send_buff(i,6,j,k) = zux(i,j,k)
    send_buff(i,7,j,k) = zuy(i,j,k)
    send_buff(i,8,j,k) = zuz(i,j,k)
  enddo
 enddo
enddo
    
gathbuffer_size = imax * js * ks * nvarout
call MPI_GATHER(send_buff, gathbuffer_size, MPI_REAL, recv_buff, gathbuffer_size, MPI_REAL, 0, MPI_COMM_ROW, mpierr)

 casefile = 'output/' // trim(prefix) // '.case'
 ipart  = 1
 label2 = 'part'
 label3 = 'block'
 label1c = 'junk line 2'

! adjust for periodic boundaries in phi if necessary
! geometry options are   0,0,0; 1,3,0; 2,4,5

  kkmax = kmax
  jjmax = jmax

  periodic_in_z = .false.
  if (ngeomz == 5) then
   if (((zzc(kmax) - zzc(1)) > 6.0).and.(kmax > 1)) then
     kkmax = kmax + 1
     periodic_in_z = .true.
     label1c = 'z periodic'
   endif
  endif

  periodic_in_y = .false.
  if (ngeomy > 2) then
   if ((zyc(jmax) - zyc(1)) > 6.0) then
     jjmax = jmax + 1
     periodic_in_y = .true.
     label1c = 'y periodic'
   endif
  endif

!--------------------------------------------------------------------------------
if (mype == 0) then   ! only mype=0 writes to case and geometry files

 write(8,6001) trim(prefix)//tmp1, time, ncycle

 if (nfile == 1000) then  ! create geometry file and open/start case file

!######################  CASEFILE   ########################

  open(unit=14,file=casefile)
  write(14,*) 'FORMAT'
  write(14,*) 'type: ensight gold'
  write(14,*) 'GEOMETRY'
  write(14,*) 'model: ', trim(prefix)//'.geo' !, '     change_coords_only'
  write(14,*) 'VARIABLE'
  write(14,*) 'scalar per node:    density   ', trim(prefix) // '_****.den'
  write(14,*) 'scalar per node:    pressure  ', trim(prefix) // '_****.pre'
  write(14,*) 'scalar per node:    gamma     ', trim(prefix) // '_****.gam'
  write(14,*) 'scalar per node:    age       ', trim(prefix) // '_****.age'
  write(14,*) 'scalar per node:    color     ', trim(prefix) // '_****.clr'
  write(14,*) 'vector per node:    velocity  ', trim(prefix) // '_****.vel'
  write(14,*) 'TIME'
  write(14,*) 'time set:    1'
  write(14,*) 'number of steps:    ', n_pframes   ! number of data files calculated in vhone.f90
  write(14,*) 'filename start number:  ', nfile + 1000
  write(14,*) 'filename increment:    1'
  write(14,fmt="('time values: ')",advance="no") 
  print *, "writing case file"

!######################  GEOFILE   ########################


  ALLOCATE(var(imax,jjmax,kkmax) )  ! allocate buffer for curvilinear coordinates and data arrays

! set up labels for EnSight files
  if (ngeomy==3) then
    label1b = 'cylindrical geometry file from VH-1'
    label1h = 'block curvilinear'
  else if (ngeomy==4) then
    label1b = 'spherical geometry file from VH-1'
    label1h = 'block curvilinear'
  else
    label1b = 'cartesian geometry file from VH-1'
    label1h = 'block rectilinear'
  endif
  label1a = 'Fortran Binary'
  label1d = 'node id off'
  label1e = 'element id off'
  label1f = 'part '
  label1g = 'cartesian coordinates'

  geomfile = 'output/' // trim(prefix) // '.geo'
  open(unit=15,file=geomfile,form='unformatted')
  print *, "writing geometry file"
  write(15) label1a
  write(15) label1b
  write(15) label1c
  write(15) label1d
  write(15) label1e
  write(15) label1f
  write(15) ipart
  write(15) label1g
  write(15) label1h
  write(15) imax, jjmax, kkmax 

  if (ngeomy < 3) then

    sx = zxc(:imax)
    sy = zyc(:jmax)
    sz = zzc(:kmax)
    write(15) sx
    write(15) sy
    write(15) sz

  else if (ngeomy == 3) then

    do k = 1, kmax
     do j = 1, jmax
      do i = 1, imax
       var(i,j,k) = zxc(i)*cos(zyc(j))/zxc(imax)
      enddo
     enddo
    enddo

    if (periodic_in_y) then
     do k = 1, kmax
      do i = 1, imax
        var(i,jjmax,k) = var(i,1,k)
      enddo
     enddo
    endif

    write(15) var

    do k = 1, kmax
     do j = 1, jmax
      do i = 1, imax
       var(i,j,k) = zxc(i)*sin(zyc(j))/zxc(imax)
      enddo
     enddo
    enddo

    if (periodic_in_y) then
     do k = 1, kmax
      do i = 1, imax
        var(i,jjmax,k) = var(i,1,k)
      enddo
     enddo
    endif

    write(15) var

    do k = 1, kmax
     do j = 1, jmax
      do i = 1, imax
       var(i,j,k) = zzc(k)
      enddo
     enddo
    enddo

    if (periodic_in_y) then
     do k = 1, kmax
      do i = 1, imax
        var(i,jjmax,k) = var(i,1,k)
      enddo
     enddo
    endif

    write(15) var

  else

    do k = 1, kmax
     do j = 1, jmax
      do i = 1, imax
       var(i,j,k) = zxc(i)*sin(zyc(j))*cos(zzc(k))
      enddo
     enddo
    enddo

    if (periodic_in_z) then
     do j = 1, jmax
      do i = 1, imax
        var(i,j,kkmax) = var(i,j,1)
      enddo
     enddo
    endif

    write(15) var

    do k = 1, kmax
     do j = 1, jmax
      do i = 1, imax
       var(i,j,k) = zxc(i)*sin(zyc(j))*sin(zzc(k))
      enddo
     enddo
    enddo

    if (periodic_in_z) then
     do j = 1, jmax
      do i = 1, imax
        var(i,j,kkmax) = var(i,j,1)
      enddo
     enddo
    endif

    write(15) var

    do k = 1, kmax
     do j = 1, jmax
      do i = 1, imax
       var(i,j,k) = zxc(i)*cos(zyc(j))
      enddo
     enddo
    enddo

    if (periodic_in_z) then
     do j = 1, jmax
      do i = 1, imax
        var(i,j,kkmax) = var(i,j,1)
      enddo
     enddo
    endif

    write(15) var

  endif
  close(15)

  DEALLOCATE( var )

 else

  open(unit=14,file=casefile,form='formatted',position='append')

 endif

 write(14,fmt="(1pe13.5)",advance="no") time
 if ((nfile - 10*(nfile/10)) == 0) write(14,*)

 close(14)


endif

!####################  DATA FILES  ##############################

! only the mypey=0 procs unload receive buffer and write data to disk

if (mypey == 0) then

 if (ndim == 2) then

  ALLOCATE(var(imax,jjmax,kkmax) )  ! allocate buffer for data arrays (here kkmax = 1)

! For 2D create filename from integer nfile and prefix such that filename looks like prefix_1000.den
  write(tmp2,"(i4)") nfile + 1000
  nfile = nfile + encount


  !!!! DENSITY !!!!!
  label4 = 'density'
  filename = 'output/' // trim(prefix) // '_' // tmp2 // '.den'
  do m = 1, npey
   do j = 1, js
    jsk = (m-1)*js + j
    do i = 1, imax
      var(i,jsk,1) = recv_buff(i,1,j,1,m)
    enddo
   enddo
  enddo
  if (periodic_in_y) then
    do i = 1, imax
      var(i,jjmax,1) = var(i,1,1)
    enddo
  endif
  open(unit=15,file=filename,form='unformatted')
    write(15) label4
    write(15) label2
    write(15) ipart
    write(15) label3
    write(15) var
  close(15)


  !!! PRESSURE !!!
  label4 = 'pressure'
  filename = 'output/' // trim(prefix) // '_' // tmp2 // '.pre'
  do m = 1, npey
   do j = 1, js
    jsk = (m-1)*js + j
    do i = 1, imax
      var(i,jsk,1) = recv_buff(i,2,j,1,m)
    enddo
   enddo
  enddo
  if (periodic_in_y) then
    do i = 1, imax
      var(i,jjmax,1) = var(i,1,1)
    enddo
  endif
  open(unit=15,file=filename,form='unformatted')
    write(15) label4
    write(15) label2
    write(15) ipart
    write(15) label3
    write(15) var
  close(15)


  !!! GAMMA !!!
  label4 = 'gamma'
  filename = 'output/' // trim(prefix) // '_' // tmp2 // '.gam'
  do m = 1, npey
   do j = 1, js
    jsk = (m-1)*js + j
    do i = 1, imax
      var(i,jsk,1) = recv_buff(i,3,j,1,m)
    enddo
   enddo
  enddo
  if (periodic_in_y) then
    do i = 1, imax
      var(i,jjmax,1) = var(i,1,1)
    enddo
  endif
  open(unit=15,file=filename,form='unformatted')
    write(15) label4
    write(15) label2
    write(15) ipart
    write(15) label3
    write(15) var
  close(15)


  !!! AGE !!!
  label4 = 'age'
  filename = 'output/' // trim(prefix) // '_' // tmp2 // '.age'
  do m = 1, npey
   do j = 1, js
    jsk = (m-1)*js + j
    do i = 1, imax
      var(i,jsk,1) = recv_buff(i,4,j,1,m)
    enddo
   enddo
  enddo
  if (periodic_in_y) then
    do i = 1, imax
      var(i,jjmax,1) = var(i,1,1)
    enddo
  endif
  open(unit=15,file=filename,form='unformatted')
    write(15) label4
    write(15) label2
    write(15) ipart
    write(15) label3
    write(15) var
  close(15)


  !!! COLOR !!!
  label4 = 'color'
  filename = 'output/' // trim(prefix) // '_' // tmp2 // '.clr'
  do m = 1, npey
   do j = 1, js
    jsk = (m-1)*js + j
    do i = 1, imax
      var(i,jsk,1) = recv_buff(i,5,j,1,m)
    enddo
   enddo
  enddo
  if (periodic_in_y) then
    do i = 1, imax
      var(i,jjmax,1) = var(i,1,1)
    enddo
  endif
  open(unit=15,file=filename,form='unformatted')
    write(15) label4
    write(15) label2
    write(15) ipart
    write(15) label3
    write(15) var
  close(15)


  !!! VELOCITY !!!
  label4 = 'velocity'
  filename = 'output/' // trim(prefix) // '_' // tmp2 // '.vel'
  open(unit=15,file=filename,form='unformatted')
   write(15) label4
   write(15) label2
   write(15) ipart
   write(15) label3

   do m = 1, npey
    do j = 1, js
     jsk = (m-1)*js + j
     do i = 1, imax
       var(i,jsk,1) = recv_buff(i,3,j,1,m) !# 1 was k. Possibly needs to be for 3D
     enddo
    enddo
   enddo
   if (periodic_in_y) then
    do i = 1, imax
      var(i,jjmax,1) = var(i,1,1)
    enddo
   endif
   write(15) var


   do m = 1, npey
    do j = 1, js
     jsk = (m-1)*js + j
     do i = 1, imax
       var(i,jsk,1) = recv_buff(i,4,j,1,m) !# 1 was k. Possibly needs to be for 3D
     enddo
    enddo
   enddo
   if (periodic_in_y) then
    do i = 1, imax
      var(i,jjmax,1) = var(i,1,1)
    enddo
   endif
   write(15) var

   var = 0.0
   write(15) var
   close(15)

   DEALLOCATE( var )

 else ! ndim = 3

  ALLOCATE(var(imax,jjmax,kmax/pez) )  ! allocate buffer for data arrays 

! Create filename from integers nfile and mypez and prefix such that filename
! looks like prefx_1000.0000.bin where 1000 is the value of nfile and 0000 is mypez

  write(tmp1,"(i4)") nfile + 1000
  nfile = nfile + 1
  write(tmp2,"(i4)") mypez
  do i = 1, 4
   if ((tmp2(i:i)) == ' ') tmp2(i:i) = '0'
  enddo
  filename = 'output/' // trim(prefix) // '_' // tmp1 // '.' // tmp2 // '.bin'
  open(unit=15,file=filename,form='unformatted')
  do nv = 1, nvarout

   do m = 1, npey
    do k = 1, ks
     do j = 1, js
      jsk = (m-1)*js + j
      do i = 1, imax
        var(i,jsk,k) = recv_buff(i,nv,j,k,m)
      enddo
     enddo
    enddo
   enddo
   if (periodic_in_y) then
    do k = 1, ks
     do i = 1, imax
       var(i,jjmax,k) = var(i,1,k)
     enddo
    enddo
   endif
   write(15) var

  enddo
  close(15)

  DEALLOCATE( var )

 endif

endif

6001 format('Wrote files for ',a12,' to disk at time =',1pe12.5,' (ncycle =', i8,')')

return
end



