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



FUNCTION vertical_shear(u1,v1,u2,v2) ! First two terms refer to higher level
  REAL,INTENT(in) :: u1(:,:), v1(:,:), u2(:,:), v2(:,:)
  REAL :: vertical_shear(SIZE(u1,1),SIZE(u1,2))

  vertical_shear = (SQRT(u1**2+v1**2)-SQRT(u2**2+v2**2))
  
END FUNCTION vertical_shear



! is this the right definition of jet?

FUNCTION jet(u, v)
  REAL,INTENT(in) :: u(:,:), v(:,:)
  REAL :: jet(SIZE(u,1),SIZE(u,2))

  jet = SQRT(u**2+v**2)
  
END FUNCTION jet



FUNCTION kindex(t850,t500,td850,t700,td700)
  REAL,INTENT(in) :: t850(:,:), t500(:,:), td850(:,:), t700(:,:), td700(:,:)
  REAL :: kindex(SIZE(t850,1),SIZE(t850,2)) !!! ???? Che size ha?
  
  kindex = ((t850 -t500) + td850 - (t700 - td700) -273.15)
  
END FUNCTION kindex



FUNCTION mcs(LI,Sh03,AvvT700)
  REAL,INTENT(in) :: LI(:,:), Sh03(:,:), AvvT700(:,:)
  REAL :: mcs(SIZE(LI,1),SIZE(LI,2))

  mcs = ((-(LI + 4.4))/3.3 + (Sh03 - 11.5)/5 + (AvvT700 - 0.000045)/0.000073)

END FUNCTION mcs


FUNCTION TWC(p,tv,q,dz)
  REAL,INTENT(in) :: p,tv,q,dz
 ! INTEGER, INTENT(in) :: np 
 ! REAL :: dz
  REAL :: TWC

  TWC = ((p/(287.05*tv))*q*dz)
  

END FUNCTION TWC


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

FUNCTION mask_average_miss(field, mask, nzones)
REAL,INTENT(in) :: field(:,:)
INTEGER,INTENT(in) :: mask(:,:)
INTEGER,INTENT(in) :: nzones
REAL :: mask_average_miss(nzones)

INTEGER :: i, nval

DO i = 1, nzones
  nval = COUNT(mask == i .AND. c_e(field))
  IF (nval > 0) THEN
    mask_average_miss(i) = SUM(field, mask=(mask == i .AND. c_e(field)))/nval
  ELSE
    mask_average_miss(i) = rmiss
  ENDIF
ENDDO

END FUNCTION mask_average_miss

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
    mask_gt_threshold(i) = REAL(COUNT(field > thr .AND. mask == i))*100./nval
  ELSE
    mask_gt_threshold(i) = rmiss
  ENDIF
ENDDO

END FUNCTION mask_gt_threshold

FUNCTION mask_gt_threshold_less(field, mask, nzones, thr)
REAL,INTENT(in) :: field(:,:)
INTEGER,INTENT(in) :: mask(:,:)
INTEGER,INTENT(in) :: nzones
REAL,INTENT(in) :: thr
REAL :: mask_gt_threshold_less(nzones)

INTEGER :: i, nval

DO i = 1, nzones
  nval = COUNT(mask == i)
  IF (nval > 0) THEN
    mask_gt_threshold_less(i) = REAL(COUNT(field < thr .AND. mask == i))*100./nval
  ELSE
    mask_gt_threshold_less(i) = rmiss
  ENDIF
ENDDO

END FUNCTION mask_gt_threshold_less

END MODULE misc_computations

! export LOG4C_PRIORITY=info
! export LOG4C_APPENDER=stderr
PROGRAM prodsim_thunderstorm_index_isobar
! input
! * lifted index già calcolato (SLI tab 203? param 147)
! * omega a 700hPa, u,v a 500hPa
! * z a 500hPa, t,td(q) a 500,700,850hPa
! output
! Lifted index, Omega 700hPa, adv geop 500hPa, K index
USE log4fortran
USE err_handling
USE missing_values
USE char_utilities
USE file_utilities
!USE phys_const
USE optionparser_class
USE vol7d_var_class
USE grid_class
USE volgrid6d_var_class
USE vol7d_level_class
USE gridinfo_class
USE volgrid6d_class
USE volgrid6d_class_compute
USE alchimia
USE volgrid6d_alchimia_class
USE termo
USE misc_computations
USE termolib
USE datetime_class
IMPLICIT NONE

INTEGER :: category, ier, k, f, j, i, nt
CHARACTER(len=512) :: a_name, input_file, output_file, mask_file, output_csv
TYPE(optionparser) :: opt
INTEGER :: optind, optstatus
LOGICAL :: version, ldisplay, lappend
TYPE(volgrid6d) :: volgridisobar, volgridpress, volgridmask, volgrid_tmp
TYPE(arrayof_gridinfo):: gridinfoisobar
TYPE(fndsv) :: vfn, vfnoracle
CHARACTER(len=10), ALLOCATABLE :: vl(:)
REAL,ALLOCATABLE :: dxm1(:,:), dym1(:,:)
INTEGER,ALLOCATABLE :: intmask(:,:)
INTEGER :: nzones
! variable indices
INTEGER :: ip, iu, iv, iw, it, iqv, iqc, iqi, it2, itd2, iu10, iv10, itp, itd, irh, iomega, ih, ili,  icape, ih0
! level indices
INTEGER :: l500, l700, l850, lfullatm, l10m, lmu, lml, l0deg
! fields to be computed
REAL,ALLOCATABLE :: geopadvection12h(:,:), kind(:,:)

! spatialised fields
REAL,ALLOCATABLE :: ind_omega700(:), ind_geopadv(:), ind_kindex(:), ind_sli(:), &
 ind_deepshear(:), ind_h0(:), ind_cape_mu(:), ind_cape_ml(:)

CHARACTER(len=12) :: filetimename

! output file
TYPE(csv_record) :: csv_writer

! initialise logging
CALL l4f_launcher(a_name,a_name_force='prodsim_thunderstorm_index')
ier=l4f_init()
! set a_name
category=l4f_category_get(TRIM(a_name)//'.main')

! define the option parser
opt = optionparser_new(description_msg= &
 'Preprocess isobaric gridded data volumes in grib format for computing thunderstorm indices.', &
 usage_msg='Usage: prodsim_thunderstorm_index_isobar [options] inputisobar outputfile')

CALL optionparser_add(opt, ' ', 'mask-file', mask_file, default='', &
 help='mask file for spatial averaging of the results, the file can to be generated &
 &from a shapefile using `vg6d_transform --trans-type=maskgen --sub-type=poly &
 &--coord-file=<shapefile> --coord-format=shp fieldin maskout`')

CALL optionparser_add(opt, ' ', 'append', lappend, &
 help='append output to csv file instead of truncating')

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
  WRITE(*,'(A,1X,A)')'prodsim_thunderstorm_index_isobar','0.1'
  CALL exit(0)
ENDIF

IF (optind + 1 /= iargc()) THEN
  CALL optionparser_printhelp(opt)
  CALL l4f_category_log(category,L4F_ERROR,'input and/or output file(s) missing')
  CALL raise_fatal_error()
  CALL EXIT(1)
ENDIF

IF (mask_file /= '') THEN
  CALL read_input_volume(mask_file, volgridmask)
! is this necessary for IMPLICIT allocation?
  ALLOCATE(intmask(SIZE(volgridmask%voldati,1),SIZE(volgridmask%voldati,2)))
  WHERE(c_e(volgridmask%voldati(:,:,1,1,1,1)))
    intmask = NINT(volgridmask%voldati(:,:,1,1,1,1))
  ELSEWHERE
    intmask = imiss
  END WHERE
  nzones = MAXVAL(intmask, mask=c_e(intmask))
ENDIF

! inputsurf
!CALL getarg(optind, input_file)
!CALL read_input_volume(input_file, volgrid_tmp, volgridlev)
! round the volume to flatten similar levels and timeranges
!CALL rounding(volgrid_tmp, volgridsurf, level=almost_equal_levels, nostatproc=.TRUE.)
!CALL delete(volgrid_tmp)

! inputisobar
CALL getarg(optind, input_file)
! import in 2 steps in order to keep a copy of data in a gridinfo
CALL IMPORT(gridinfoisobar, input_file)
CALL read_input_gridinfo(gridinfoisobar, volgridisobar, input_file)
! compute pressure as a variable
CALL volgrid6d_compute_vert_coord_var(volgridisobar, vol7d_level_new(100), &
 volgridpress)
! append pressure to the gridinfo and reimport in a volume
CALL export(volgridpress, gridinfoisobar, clone=.TRUE.)
CALL delete(volgridisobar)
CALL delete(volgridpress)
CALL read_input_gridinfo(gridinfoisobar, volgridisobar, TRIM(input_file)//'+pressure')
CALL delete(gridinfoisobar)
! output file
CALL getarg(optind+1, output_csv)

IF (ldisplay) THEN
  PRINT*,'input volume >>>>>>>>>>>>>>>>>>>>'
  CALL display(volgridisobar)
ENDIF

CALL register_termo(vfn)
IF (ldisplay) CALL display(vfn)
ALLOCATE(vl(8))
vl = (/'B11003', 'B11004', 'B12101', 'B12103', 'B11005', 'B10008', 'B13234', 'B13230'/)

! Fill low qv values to avoid problems in thermodynamic functions
iqv = vartable_index(volgridisobar%var, 'B13001')
IF (c_e(iqv)) THEN
  WHERE (.NOT.c_e(volgridisobar%voldati(:,:,:,1,1,iqv)) .OR. volgridisobar%voldati(:,:,:,1,1,iqv) < 1.E-7)
    volgridisobar%voldati(:,:,:,1,1,iqv) = 1.E-5
  END WHERE
ENDIF 

IF (alchemy(volgridisobar, vfn,vl,volgrid_tmp,copy=.TRUE.,vfnoracle=vfnoracle) == 0) THEN
  CALL register_termo(vfn)
  IF (ldisplay) CALL display(vfn)
  CALL delete(volgridisobar)
  volgridisobar = volgrid_tmp
ELSE
  CALL l4f_category_log(category,L4F_ERROR,'computing output variables')
  CALL raise_fatal_error()
ENDIF

! if coordinates of input grid are needed, do the following, then
! coordinates will be allocated in arrays volgrid(1)%griddim%dim%lon
! volgrid(1)%griddim%dim%lat
CALL unproj(volgridisobar%griddim)
CALL grid_metric_terms(volgridisobar%griddim%dim%lon, &
 volgridisobar%griddim%dim%lat, dxm1, dym1)

! locate variables in volume
! P B10004, U B11003, V B11004, W B11006, T B12101, QV B13001, QC B13192, QI B13193
! T2 B12101, TD2 B12102, U10 B11003, V10 B11004, TP B13011
! OMEGA B11005
!ip = vartable_index(volgridua%var, 'B10004')
iu = vartable_index(volgridisobar%var, 'B11003')
iv = vartable_index(volgridisobar%var, 'B11004')
it = vartable_index(volgridisobar%var, 'B12101')
itd = vartable_index(volgridisobar%var, 'B12103')
iqv = vartable_index(volgridisobar%var, 'B13001')
iomega = vartable_index(volgridisobar%var, 'B11005')
ih = vartable_index(volgridisobar%var, 'B10008') ! geopotential
ili = vartable_index(volgridisobar%var, 'B13234') ! surface lifted index

! Need to add also cape_ml, cape_mu and hzerocl
icape = vargrib_index(volgridisobar%var, (/imiss, 0, 7, 6/))
ih0 = vargrib_index(volgridisobar%var, (/imiss, 0, 3, 6/))

IF (ldisplay) THEN
  PRINT*,'processed volume >>>>>>>>>>>>>>>>>>>>'
  CALL display(volgridisobar)
ENDIF

!IF (.NOT.ALL(c_e((/iu,iv,it,itd,iqv,iomega,ih,ili/)))) THEN
!  CALL l4f_category_log(category,L4F_ERROR,'some output variables missing')
!  CALL raise_fatal_error()
!ENDIF

! locate pressure levels in volume
l500 = index(volgridisobar%level, vol7d_level_new(100, 50000))
l700 = index(volgridisobar%level, vol7d_level_new(100, 70000))
l850 = index(volgridisobar%level, vol7d_level_new(100, 85000))
lfullatm = INDEX(volgridisobar%level, vol7d_level_new(10))
l10m = index(volgridisobar%level, vol7d_level_new(103, 10000))
lml = index(volgridisobar%level, vol7d_level_new(192))
lmu = index(volgridisobar%level, vol7d_level_new(193))
l0deg = index(volgridisobar%level, vol7d_level_new(4, imiss, 101, imiss))

!IF (.NOT.ALL(c_e((/l500,l700,l850,lfullatm/)))) THEN
!  CALL l4f_category_log(category,L4F_ERROR,'some pressure levels missing')
!  CALL raise_fatal_error()
!ENDIF

! Omega 
IF (ALLOCATED(intmask).AND.c_e(iomega).AND.c_e(l700)) THEN
  ind_omega700 = mask_gt_threshold_less(volgridisobar%voldati(:,:,l700,1,1,iomega), intmask, nzones, -0.5)
ENDIF

! geopotential advection in 12h
IF (c_e(l500).AND.c_e(iu).AND.c_e(iv).AND.c_e(ih)) THEN
  geopadvection12h = (grid_horiz_advection( &
   volgridisobar%voldati(:,:,l500, 1, 1, iu), &
   volgridisobar%voldati(:,:,l500, 1, 1, iv), &
   volgridisobar%voldati(:,:,l500, 1, 1, ih), dxm1, dym1)*4320)

  IF (ALLOCATED(intmask)) THEN
    ind_geopadv = mask_average(geopadvection12h, intmask, nzones)
  ENDIF
ENDIF

! K index
IF (c_e(l500).AND.c_e(l700).AND.c_e(l850).AND.c_e(it).AND.c_e(itd)) THEN
  kind = kindex(volgridisobar%voldati(:,:,l850, 1, 1, it), &
   volgridisobar%voldati(:,:,l500, 1, 1, it), &
   volgridisobar%voldati(:,:,l850, 1, 1, itd), &
   volgridisobar%voldati(:,:,l700, 1, 1, it), &
   volgridisobar%voldati(:,:,l700, 1, 1, itd))
  IF (ALLOCATED(intmask)) THEN
    ind_kindex = mask_average(kind, intmask, nzones)
  ENDIF
ENDIF

! surface lifted index
IF (ALLOCATED(intmask).AND.c_e(ili)) THEN
  ind_sli = mask_average(volgridisobar%voldati(:,:,lfullatm, 1, 1, ili), intmask, nzones)
ENDIF

! deep shear index
IF (ALLOCATED(intmask).AND.c_e(l500).AND.c_e(l10m).AND.c_e(iu).AND.c_e(iv)) THEN
  ind_deepshear = mask_average(vertical_shear( &
   volgridisobar%voldati(:,:,l500,1,1,iu), &
   volgridisobar%voldati(:,:,l500,1,1,iv), &
   volgridisobar%voldati(:,:,l10m,1,1,iu), &
   volgridisobar%voldati(:,:,l10m,1,1,iv)), intmask, nzones)
ENDIF

! index h0
IF (ALLOCATED(intmask).AND.c_e(l0deg).AND.c_e(ih0)) THEN
! set points with 0c below surface to a true missing value
  WHERE (volgridisobar%voldati(:,:,l0deg,1,1,ih0) < -990.)
    volgridisobar%voldati(:,:,l0deg,1,1,ih0) = rmiss
  END WHERE
  ind_h0 = mask_average_miss(volgridisobar%voldati(:,:,l0deg,1,1,ih0), &
   intmask, nzones)
ENDIF

! cape mu
IF (ALLOCATED(intmask).AND.c_e(lmu).AND.c_e(icape)) THEN
  ind_cape_mu = mask_average(volgridisobar%voldati(:,:,lmu,1,1,icape), &
   intmask, nzones)
ENDIF

! cape ml
IF (ALLOCATED(intmask).AND.c_e(lml).AND.c_e(icape)) THEN
  ind_cape_ml = mask_average(volgridisobar%voldati(:,:,lml,1,1,icape), &
   intmask, nzones)
ENDIF


CALL getval(volgridisobar%time(1), simpledate=filetimename)
IF (lappend) THEN
  OPEN(unit=2,file=output_csv, position='append')
ELSE
  OPEN(unit=2,file=output_csv)
  WRITE(2,'(6(A,'',''))')'Data','Macroarea','%VV700', &
   'AvvGeop500','Kindex','LI','DeepShear','H0','CAPE_MU','CAPE_ML'
ENDIF


DO i = 1, nzones
  CALL init(csv_writer)
  CALL csv_record_addfield(csv_writer, filetimename(1:10))
  CALL csv_record_addfield(csv_writer, i)
  CALL addfield_if_assoc(csv_writer, ind_omega700, i)
  CALL addfield_if_assoc(csv_writer, ind_geopadv, i)
  CALL addfield_if_assoc(csv_writer, ind_kindex, i)
  CALL addfield_if_assoc(csv_writer, ind_sli, i)
  CALL addfield_if_assoc(csv_writer, ind_deepshear, i)
  CALL addfield_if_assoc(csv_writer, ind_h0, i)
  CALL addfield_if_assoc(csv_writer, ind_cape_mu, i)
  CALL addfield_if_assoc(csv_writer, ind_cape_ml, i)
  WRITE(2,'(A)')csv_record_getrecord(csv_writer)
  CALL delete(csv_writer)
ENDDO

CLOSE(2)

CALL delete(volgridisobar)

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

FUNCTION vargrib_index(varlist, vargrib)
TYPE(volgrid6d_var),INTENT(in) :: varlist(:)
INTEGER,INTENT(in) :: vargrib(4)
INTEGER :: vargrib_index

INTEGER :: i

vargrib_index = imiss
DO i = 1, SIZE(varlist)
  IF ((varlist(i)%centre == vargrib(1) .OR. .NOT.c_e(vargrib(1))) .AND. &
   (varlist(i)%discipline == vargrib(2) .OR. .NOT.c_e(vargrib(2))) .AND. &
   (varlist(i)%category == vargrib(3) .OR. .NOT.c_e(vargrib(3))) .AND. &
   (varlist(i)%number == vargrib(4) .OR. .NOT.c_e(vargrib(4)))) THEN
    vargrib_index = i
    EXIT
  ENDIF
ENDDO

END FUNCTION vargrib_index

SUBROUTINE read_input_volume(infile, outvol, comparevol)
CHARACTER(len=*),INTENT(in) :: infile
TYPE(volgrid6d),INTENT(inout) :: outvol
TYPE(volgrid6d),INTENT(in),OPTIONAL :: comparevol

TYPE(volgrid6d),POINTER :: volgrid_tmp(:)=>NULL()

CALL l4f_category_log(category,L4F_INFO,'importing volume from file: '//TRIM(infile))

CALL import(volgrid_tmp, filename=infile, decode=.TRUE., dup_mode=0, &
 time_definition=0, categoryappend='inputvol')

IF (.NOT.ASSOCIATED(volgrid_tmp)) THEN
  CALL l4f_category_log(category, L4F_ERROR, &
   'error importing volume from '//TRIM(infile))
  CALL raise_fatal_error()
ENDIF

IF (SIZE(volgrid_tmp) > 1) THEN
  CALL l4f_category_log(category, L4F_ERROR, &
   t2c(SIZE(volgrid_tmp))//' grids found in '//TRIM(infile)//', 1 expected')
  CALL raise_fatal_error()
ENDIF

IF (PRESENT(comparevol)) THEN
  IF (volgrid_tmp(1)%griddim /= comparevol%griddim) THEN
    CALL display(volgrid_tmp(1)%griddim)
    CALL display(comparevol%griddim)
    CALL l4f_category_log(category, L4F_ERROR, &
     'grid in '//TRIM(infile)//' differs from reference grid')
    CALL raise_fatal_error()
  ENDIF
ENDIF

outvol = volgrid_tmp(1)
DEALLOCATE(volgrid_tmp)

END SUBROUTINE read_input_volume


SUBROUTINE read_input_gridinfo(ingridinfo, outvol, infile, comparevol)
TYPE(arrayof_gridinfo),INTENT(in) :: ingridinfo
TYPE(volgrid6d),INTENT(inout) :: outvol
CHARACTER(len=*),INTENT(in) :: infile
TYPE(volgrid6d),INTENT(in),OPTIONAL :: comparevol

TYPE(volgrid6d),POINTER :: volgrid_tmp(:)=>NULL()

CALL l4f_category_log(category,L4F_INFO,'importing volume from gridinfo: '//TRIM(infile))

CALL import(volgrid_tmp, ingridinfo, decode=.TRUE., dup_mode=0, &
 time_definition=0, clone=.TRUE., categoryappend='inputvol')

IF (.NOT.ASSOCIATED(volgrid_tmp)) THEN
  CALL l4f_category_log(category, L4F_ERROR, &
   'error importing volume from '//TRIM(infile))
  CALL raise_fatal_error()
ENDIF

IF (SIZE(volgrid_tmp) > 1) THEN
  CALL l4f_category_log(category, L4F_ERROR, &
   t2c(SIZE(volgrid_tmp))//' grids found in '//TRIM(infile)//', 1 expected')
  CALL raise_fatal_error()
ENDIF

IF (PRESENT(comparevol)) THEN
  IF (volgrid_tmp(1)%griddim /= comparevol%griddim) THEN
    CALL display(volgrid_tmp(1)%griddim)
    CALL display(comparevol%griddim)
    CALL l4f_category_log(category, L4F_ERROR, &
     'grid in '//TRIM(infile)//' differs from reference grid')
    CALL raise_fatal_error()
  ENDIF
ENDIF

outvol = volgrid_tmp(1)
DEALLOCATE(volgrid_tmp)

END SUBROUTINE read_input_gridinfo


SUBROUTINE addfield_if_assoc(writer, arr, i)
TYPE(csv_record),INTENT(inout) :: writer
REAL,ALLOCATABLE,INTENT(in) :: arr(:)
INTEGER,intent(in) :: i

IF (ALLOCATED(arr)) THEN
  CALL csv_record_addfield(writer, arr(i))
ELSE
  CALL csv_record_addfield(writer, rmiss)
ENDIF

END SUBROUTINE addfield_if_assoc

END PROGRAM prodsim_thunderstorm_index_isobar
