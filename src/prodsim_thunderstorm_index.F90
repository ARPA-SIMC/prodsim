MODULE misc_computations
USE missing_values
USE phys_const
IMPLICIT NONE

CONTAINS

FUNCTION grid_ddx(field, dxm1)
REAL,INTENT(in) :: field(:,:), dxm1(:,:)
REAL :: grid_ddx(SIZE(field,1),SIZE(field,2))

INTEGER :: i, j

grid_ddx = 0. ! set to 0 on frame
DO j = 2, SIZE(dxm1, 2) - 1
  DO i = 2, SIZE(dxm1, 1) - 1
    grid_ddx(i,j) = (field(i+1,j)-field(i-1,j))*dxm1(i,j)
  ENDDO
ENDDO

END FUNCTION grid_ddx

FUNCTION grid_ddy(field, dym1)
REAL,INTENT(in) :: field(:,:), dym1(:,:)
REAL :: grid_ddy(SIZE(field,1),SIZE(field,2))

INTEGER :: i, j

grid_ddy = 0. ! set to 0 on frame
DO j = 2, SIZE(dym1, 2) - 1
  DO i = 2, SIZE(dym1, 1) - 1
    grid_ddy(i,j) = (field(i,j+1)-field(i,j-1))*dym1(i,j)
  ENDDO
ENDDO

END FUNCTION grid_ddy

SUBROUTINE grid_metric_terms(lon, lat, dxm1, dym1)
DOUBLE PRECISION,INTENT(in) :: lon(:,:), lat(:,:)
REAL,INTENT(out),ALLOCATABLE :: dxm1(:,:), dym1(:,:)

INTEGER :: i, j

! autoallocate does not work with -fcheck-bounds
ALLOCATE(dxm1(SIZE(lon, 1),SIZE(lon, 2)), dym1(SIZE(lon, 1),SIZE(lon, 2)))
!dxm1 = lon
!dym1 = lat
dxm1 = rmiss
dym1 = rmiss

DO j = 2, SIZE(lon, 2) - 1
  DO i = 2, SIZE(lon, 1) - 1
    dxm1(i,j) = 1./dist(lon(i-1,j), lat(i-1,j), lon(i+1,j), lat(i+1,j))
    dym1(i,j) = 1./dist(lon(i,j-1), lat(i,j-1), lon(i,j+1), lat(i,j+1))
  ENDDO
ENDDO

END SUBROUTINE grid_metric_terms

FUNCTION dist(lon1, lat1, lon2, lat2)
DOUBLE PRECISION,INTENT(in) :: lon1, lat1, lon2, lat2
REAL :: dist

REAL :: x, y

x = REAL((lon2 - lon1)*COS(((lat1 + lat2)/2.)*degrad))
y = REAL(lat2 - lat1)
dist = SQRT(x**2 + y**2)*degrad*rearth

END FUNCTION dist

FUNCTION grid_vorticity_z(u, v, dxm1, dym1)
REAL,INTENT(in) :: u(:,:), v(:,:), dxm1(:,:), dym1(:,:)
REAL :: grid_vorticity_z(SIZE(u,1),SIZE(u,2))

grid_vorticity_z = grid_ddx(v, dxm1) - grid_ddy(u, dym1)

END FUNCTION grid_vorticity_z

FUNCTION grid_horiz_divergence(u, v, dxm1, dym1)
REAL,INTENT(in) :: u(:,:), v(:,:), dxm1(:,:), dym1(:,:)
REAL :: grid_horiz_divergence(SIZE(u,1),SIZE(u,2))

grid_horiz_divergence = grid_ddx(u, dxm1) + grid_ddy(v, dym1)

END FUNCTION grid_horiz_divergence

FUNCTION grid_horiz_advection(u, v, field, dxm1, dym1)
REAL,INTENT(in) :: u(:,:), v(:,:), field(:,:), dxm1(:,:), dym1(:,:)
REAL :: grid_horiz_advection(SIZE(u,1),SIZE(u,2))

grid_horiz_advection = u*grid_ddx(field, dxm1) + v*grid_ddy(field, dym1)

END FUNCTION grid_horiz_advection

!FUNCTION compute_vert_integral(var, p, t, z)
!REAL,INTENT(in) :: var(:,:,:), p(:,:,:), t(:,:,:), z(:,:,:)
!REAL,ALLOCATABLE :: compute_vert_integral(:,:)
!
!!compute_vert_integral = compute_vert_integral + var/gearth*dp
!compute_vert_integral = compute_vert_integral + var*rd*t/p*dz
!
!END FUNCTION compute_vert_integral

FUNCTION mask_average(field, mask, nzones)
REAL,INTENT(in) :: field(:,:)
INTEGER,INTENT(in) :: mask(:,:)
INTEGER,INTENT(in) :: nzones
REAL :: mask_average(nzones)

INTEGER :: i, nval

DO i = 1, nzones
  nval = COUNT(mask == i)
  IF (nval > 0) THEN
    mask_average(i) = SUM(field, mask=(mask == i))/nval
  ELSE
    mask_average(i) = rmiss
  ENDIF
ENDDO

END FUNCTION mask_average

FUNCTION mask_gt_threshold(field, mask, nzones, thr)
REAL,INTENT(in) :: field(:,:)
INTEGER,INTENT(in) :: mask(:,:)
INTEGER,INTENT(in) :: nzones
REAL,INTENT(in) :: thr
REAL :: mask_gt_threshold(nzones)

INTEGER :: i, nval

DO i = 1, nzones
  nval = COUNT(mask == i)
  IF (nval > 0) THEN
    mask_gt_threshold(i) = REAL(COUNT(field > thr .AND. mask == i))/nval
  ELSE
    mask_gt_threshold(i) = rmiss
  ENDIF
ENDDO

END FUNCTION mask_gt_threshold

END MODULE misc_computations

! export LOG4C_PRIORITY=info
! export LOG4C_APPENDER=stderr
PROGRAM prodsim_thunderstorm_index
USE log4fortran
USE err_handling
USE missing_values
USE char_utilities
!USE phys_const
USE optionparser_class
USE vol7d_var_class
USE grid_class
USE volgrid6d_var_class
USE vol7d_level_class
USE volgrid6d_class
USE misc_computations
IMPLICIT NONE

INTEGER :: category, ier
CHARACTER(len=512) :: a_name, input_file, output_file, mask_file
TYPE(optionparser) :: opt
INTEGER :: optind, optstatus
LOGICAL :: version, ldisplay
TYPE(volgrid6d),POINTER :: volgrid_tmp(:), volgrid_tmpr(:)
TYPE(volgrid6d) :: volgridz, volgridsurf, volgridua, volgridmask
REAL,ALLOCATABLE :: dxm1(:,:), dym1(:,:)
INTEGER,ALLOCATABLE :: intmask(:,:)
INTEGER :: nzones
! variable indices
INTEGER :: ip, iu, iv, iw, it, iqv, iqc, iqi, it2, itd2, iu10, iv10, itp
! level indices
INTEGER :: l500 !...
! fields to be computed
REAL,ALLOCATABLE :: vorticity(:,:), tadvection(:,:)
! spatialised fields
REAL,ALLOCATABLE :: example_index1(:), example_index2(:)

! initialise logging
CALL l4f_launcher(a_name,a_name_force='prodsim_vg6d_tcorr')
ier=l4f_init()
! set a_name
category=l4f_category_get(a_name//'.main')

! define the option parser
opt = optionparser_new(description_msg= &
 'Preprocess data for computing thunderstorm index from 3d volumes of grib data.', &
 usage_msg='Usage: prodsim_thunderstorm_index [options] inputz inputsurf inputua outputfile')

! for generating mask file:
! vg6d_transform --trans-type=maskgen --sub-type=poly \
!  --coord-file=/usr/local/share/libsim/macroaree_er --coord-FORMAT=shp \
!  orog.grib mask.grib
CALL optionparser_add(opt, ' ', 'mask-file', mask_file, default='', &
 help='mask file for spatial averaging of the results, if empty computations &
 &will be performed point by point')

! display option
CALL optionparser_add(opt, 'd', 'display', ldisplay, help= &
 'briefly display the data volumes imported')

! help options
CALL optionparser_add_help(opt, 'h', 'help', help='show an help message and exit')
CALL optionparser_add(opt, ' ', 'version', version, help='show version and exit')

! parse options and check for errors
CALL optionparser_parse(opt, optind, optstatus)

IF (optstatus == optionparser_help) THEN
  CALL exit(0) ! generate a clean manpage
ELSE IF (optstatus == optionparser_err) THEN
  CALL l4f_category_log(category,L4F_ERROR,'in command-line arguments')
  CALL raise_fatal_error()
ENDIF
IF (version) THEN
  WRITE(*,'(A,1X,A)')'prodsim_thunderstorm_index','0.1'
  CALL exit(0)
ENDIF

IF (optind+3 /= iargc()) THEN
  CALL optionparser_printhelp(opt)
  CALL l4f_category_log(category,L4F_ERROR,'input and/or output file(s) missing')
  CALL raise_fatal_error()
  CALL EXIT(1)
ENDIF

! inputz
CALL getarg(optind, input_file)
CALL l4f_category_log(category,L4F_INFO,'inputz file: '//TRIM(input_file))

CALL IMPORT(volgrid_tmp, filename=input_file, decode=.TRUE., dup_mode=0, &
 time_definition=0, categoryappend='inputz')
IF (SIZE(volgrid_tmp) > 1) THEN
  CALL l4f_category_log(category, L4F_ERROR, &
   t2c(SIZE(volgrid_tmp))//' grids found in '//TRIM(input_file))
  CALL raise_fatal_error()
ENDIF
volgridz = volgrid_tmp(1)
DEALLOCATE(volgrid_tmp)

! inputsurf
CALL getarg(optind+1, input_file)
CALL l4f_category_log(category,L4F_INFO,'inputsurf file: '//TRIM(input_file))

CALL IMPORT(volgrid_tmp, filename=input_file, decode=.TRUE., dup_mode=0, &
 time_definition=0, categoryappend='inputsurf')
IF (SIZE(volgrid_tmp) > 1) THEN
  CALL l4f_category_log(category, L4F_ERROR, &
   t2c(SIZE(volgrid_tmp))//' grids found in '//TRIM(input_file))
  CALL raise_fatal_error()
ENDIF
! round the volume to flatten similar level and timeranges
CALL rounding(volgrid_tmp, volgrid_tmpr, level=almost_equal_levels, nostatproc=.TRUE.)
CALL delete(volgrid_tmp)
volgridsurf = volgrid_tmpr(1)
DEALLOCATE(volgrid_tmpr)
IF (volgridsurf%griddim /= volgridz%griddim) THEN
  CALL display(volgridsurf%griddim)
  CALL display(volgridz%griddim)
  CALL l4f_category_log(category, L4F_ERROR, &
   'grid in '//TRIM(input_file)//' differs from grid in inputz file')
  CALL raise_fatal_error()
ENDIF

! inputua
CALL getarg(optind+2, input_file)
CALL l4f_category_log(category,L4F_INFO,'inputua file: '//TRIM(input_file))

CALL IMPORT(volgrid_tmp, filename=input_file, decode=.TRUE., dup_mode=0, &
 time_definition=0, categoryappend='inputua')
IF (SIZE(volgrid_tmp) > 1) THEN
  CALL l4f_category_log(category, L4F_ERROR, &
   t2c(SIZE(volgrid_tmp))//' grids found in '//TRIM(input_file))
  CALL raise_fatal_error()
ENDIF
volgridua = volgrid_tmp(1)
DEALLOCATE(volgrid_tmp)
IF (volgridua%griddim /= volgridz%griddim) THEN
  CALL display(volgridua%griddim)
  CALL display(volgridz%griddim)
  CALL l4f_category_log(category, L4F_ERROR, &
   'grid in '//TRIM(input_file)//' differs from grid in inputz file')
  CALL raise_fatal_error()
ENDIF
IF (SIZE(volgridua%level) /= SIZE(volgridz%level)) THEN
  CALL l4f_category_log(category, L4F_ERROR, &
   t2c(SIZE(volgridua%level))//'/'//t2c(SIZE(volgridz%level)))
  CALL raise_fatal_error()
ENDIF

IF (mask_file /= '') THEN
  CALL IMPORT(volgrid_tmp, filename=mask_file, decode=.TRUE., dup_mode=0, &
 time_definition=0, categoryappend='inputmask')
  IF (SIZE(volgrid_tmp) > 1) THEN
    CALL l4f_category_log(category, L4F_ERROR, &
     t2c(SIZE(volgrid_tmp))//' grids found in '//TRIM(mask_file))
    CALL raise_fatal_error()
  ENDIF
  volgridmask = volgrid_tmp(1)
  DEALLOCATE(volgrid_tmp)
  IF (volgridmask%griddim /= volgridz%griddim) THEN
    call display(volgridmask%griddim)
    call display(volgridz%griddim)
    CALL l4f_category_log(category, L4F_ERROR, &
     'grid in '//TRIM(mask_file)//' differs from grid in inputz file')
    CALL raise_fatal_error()
  ENDIF
! is this necessary for IMPLICIT allocation?
  ALLOCATE(intmask(SIZE(volgridmask%voldati,1),SIZE(volgridmask%voldati,2)))
  WHERE(c_e(volgridmask%voldati(:,:,1,1,1,1)))
    intmask = NINT(volgridmask%voldati(:,:,1,1,1,1))
  ELSEWHERE
    intmask = imiss
  END WHERE
  nzones = MAXVAL(intmask, mask=c_e(intmask))
ENDIF
  

IF (ldisplay) THEN
  PRINT*,'input volume >>>>>>>>>>>>>>>>>>>>'
  CALL display(volgridz)
  CALL display(volgridsurf)
  CALL display(volgridua)
ENDIF

!j = volgrid(1)%griddim%dim%nx*volgrid(1)%griddim%dim%ny
!gridsize = SIZE(volgrid(1)%voldati,1)*SIZE(volgrid(1)%voldati,2)
!IF (j /= gridsize .OR. gridsize == 0) THEN
!  CALL l4f_category_log(category, L4F_ERROR, &
!   'error inconsistent grid sizes: '//t2c(j)//','//t2c(gridsize)//' in '//TRIM(input_file))
!  CALL raise_fatal_error()
!ENDIF

! if coordinates of input grid are needed, do the following, then
! coordinates will be allocated in arrays volgrid(1)%griddim%dim%lon
! volgrid(1)%griddim%dim%lat
CALL unproj(volgridz%griddim)
CALL grid_metric_terms(volgridz%griddim%dim%lon, volgridz%griddim%dim%lat, &
 dxm1, dym1)

! compute variable indices
! P B10004, U B11003, V B11004, W B11006, T B12101, QV B13001, QC B13192, QI B13193
! T2 B12101, TD2 B12102, U10 B11003, V10 B11004, TP B13011
! OMEGA B11005
ip = vartable_index(volgridua%var, 'B10004')
iu = vartable_index(volgridua%var, 'B11003')
iv = vartable_index(volgridua%var, 'B11004')
iw = vartable_index(volgridua%var, 'B11006')
it = vartable_index(volgridua%var, 'B12101')
iqv = vartable_index(volgridua%var, 'B13001')
iqc = vartable_index(volgridua%var, 'B13192')
iqi = vartable_index(volgridua%var, 'B13193')

it2 = vartable_index(volgridsurf%var, 'B12101')
itd2 = vartable_index(volgridsurf%var, 'B12102')
iu10 = vartable_index(volgridsurf%var, 'B11003')
iv10 = vartable_index(volgridsurf%var, 'B11004')
itp = vartable_index(volgridsurf%var, 'B13011')

! associate vertical level indices to pseudo pressure levels
l500 = 30 ! example

! examples opf grid computations
! remember voldati(x,y,level,time,timerange,var)
IF (c_e(iu) .AND. c_e(iv) .AND. c_e(it)) THEN

  vorticity = grid_vorticity_z( &
   volgridua%voldati(:,:,l500, 1, 1, iu), &
   volgridua%voldati(:,:,l500, 1, 1, iv), &
   dxm1, dym1)

  tadvection = grid_horiz_advection( &
   volgridua%voldati(:,:,l500, 1, 1, iu), &
   volgridua%voldati(:,:,l500, 1, 1, iv), &
   volgridua%voldati(:,:,l500, 1, 1, it), &
   dxm1, dym1)

! example of spatialization with mask
  IF (ALLOCATED(intmask)) THEN
! implicit allocation
    example_index1 = mask_average(vorticity, intmask, nzones)
    example_index2 = mask_gt_threshold(tadvection, intmask, nzones, 0.001)
  ENDIF
ELSE
  CALL l4f_category_log(category, L4F_ERROR, &
   'u, v or t missing in upper air data')
  CALL raise_fatal_error()
ENDIF

! output
CALL getarg(optind+3, output_file)
CALL l4f_category_log(category,L4F_INFO,'output file: '//TRIM(output_file))
! output done here...
IF (ALLOCATED(example_index1)) THEN
  PRINT*,'Average vorticity'
  PRINT*,example_index1
ENDIF
IF (ALLOCATED(example_index2)) THEN
  PRINT*,'Advection over threshold'
  PRINT*,example_index2
ENDIF

CALL delete(volgridz)
CALL delete(volgridsurf)
CALL delete(volgridua)

CONTAINS

FUNCTION vartable_index(varlist, btable)
TYPE(volgrid6d_var),INTENT(in) :: varlist(:)
CHARACTER(len=*) :: btable
INTEGER :: vartable_index

INTEGER :: i
TYPE(vol7d_var) :: varbufr

vartable_index = imiss
DO i = 1, SIZE(varlist)
  varbufr = convert(varlist(i))
  IF (varbufr%btable == btable) THEN
    vartable_index = i
    EXIT
  ENDIF
ENDDO

END FUNCTION vartable_index

END PROGRAM prodsim_thunderstorm_index
