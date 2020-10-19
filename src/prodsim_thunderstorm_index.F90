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
PROGRAM prodsim_thunderstorm_index
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
USE volgrid6d_class
USE misc_computations
USE termolib
USE datetime_class
IMPLICIT NONE

INTEGER :: category, ier, k, f, j, i, nt
CHARACTER(len=512) :: a_name, input_file, output_file, mask_file, output_csv
TYPE(optionparser) :: opt
INTEGER :: optind, optstatus
LOGICAL :: version, ldisplay
TYPE(volgrid6d) :: volgridz, volgridsurf, volgridua, volgridmask, volgridlev, volgrid_tmp
REAL,ALLOCATABLE :: dxm1(:,:), dym1(:,:)
INTEGER,ALLOCATABLE :: intmask(:,:)
INTEGER :: nzones
! variable indices
INTEGER :: ip, iu, iv, iw, it, iqv, iqc, iqi, it2, itd2, iu10, iv10, itp, itd, irh, iomega, ih
! level indices
INTEGER :: l250, l300, l500, l700, l850, l925, l950, l1000 !...
! fields to be computed
REAL,ALLOCATABLE :: vorticity(:,:), omega(:,:), vertwindshear(:,:), SLI(:,:), jets(:,:)
REAL,ALLOCATABLE :: inibition(:,:), kind(:,:), tadvection700(:,:), mcsind(:,:), relhum(:,:), lfc(:,:), lcl(:,:)
REAL,ALLOCATABLE :: templcl(:,:), tdewadvection12h(:,:), equilibrium(:,:), capeindex(:,:), wclayer(:,:), virtualt(:,:)
REAL,ALLOCATABLE :: totalwc(:,:), tadvection12h(:,:), hgeopotential500(:,:), geopadvection12h(:,:), teta(:,:), tamblcl(:,:)
REAL,ALLOCATABLE :: tsatlcl(:,:), samix(:,:), pccl(:,:), tambccl(:,:), tsatccl(:,:), equilibriumccl(:,:)
REAL,ALLOCATABLE :: volumeinv(:,:,:,:,:,:)

! spatialised fields
REAL,ALLOCATABLE :: example_index1(:), example_omega(:), example_vws(:), example_sli(:)
REAL,ALLOCATABLE :: example_jet(:), example_cin(:), example_kindex(:), example_mcsindex(:), example_relhum(:), example_lfc(:)
REAL,ALLOCATABLE :: example_lcl(:), example_equilibrium(:), example_cape(:), example_twc(:), example_tadv12h(:)
REAL,ALLOCATABLE :: example_tdewadv12h(:), example_hgeopot(:), example_geopavv(:), example_ccl(:), example_omega300(:)
REAL,ALLOCATABLE :: example_vws500(:), example_jet925(:)

! pressure and temperature values at levels
REAL, DIMENSION(32) ::  aw ! allocate correctly!!!

CHARACTER(len=12) :: filetimename

! output file
TYPE(csv_record) :: csv_writer

! initialise logging
CALL l4f_launcher(a_name,a_name_force='prodsim_thunderstorm_index')
ier=l4f_init()
! set a_name
category=l4f_category_get(a_name//'.main')

! define the option parser
opt = optionparser_new(description_msg= &
     'Preprocess data for computing thunderstorm indices from 3d volumes of grib data.', &
!!!!! SE USO OPERAZIONIDATA_PARALLEL DECOMMENTARE E COMMENTARE RIGA SUCCESSIVA     
! usage_msg='Usage: prodsim_thunderstorm_index [options] inputz inputsurf inputua inputlev inputcsv outputfile')
 usage_msg='Usage: prodsim_thunderstorm_index [options] inputz inputsurf inputua inputlev outputfile')

! for generating mask file:
! vg6d_transform --trans-type=maskgen --sub-type=poly \
!  --coord-file=/usr/local/share/libsim/macroaree_er --coord-FORMAT=shp \
!  orog.grib mask.grib
CALL optionparser_add(opt, ' ', 'mask-file', mask_file, default='', &
 help='mask file for spatial averaging of the results, if empty computations &
 &will be performed point by point')

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
  WRITE(*,'(A,1X,A)')'prodsim_thunderstorm_index','0.2'
  CALL exit(0)
ENDIF

IF (optind+5 /= iargc()) THEN
  CALL optionparser_printhelp(opt)
  CALL l4f_category_log(category,L4F_ERROR,'input and/or output file(s) missing')
  CALL raise_fatal_error()
  CALL EXIT(1)
ENDIF

IF (optind+4 /= iargc()) THEN
  CALL optionparser_printhelp(opt)
  CALL l4f_category_log(category,L4F_ERROR,'input and/or output file(s) missing')
  CALL raise_fatal_error()
  CALL EXIT(1)
ENDIF

! inputz
CALL getarg(optind, input_file)
CALL read_input_volume(input_file, volgridz)

! inputsurf
CALL getarg(optind+1, input_file)
CALL read_input_volume(input_file, volgrid_tmp, volgridz)
! round the volume to flatten similar level and timeranges
CALL rounding(volgrid_tmp, volgridsurf, level=almost_equal_levels, nostatproc=.TRUE.)
CALL delete(volgrid_tmp)

! inputua
CALL getarg(optind+2, input_file)
CALL read_input_volume(input_file, volgridua, volgridz)
IF (SIZE(volgridua%level) /= SIZE(volgridz%level)) THEN
  CALL l4f_category_log(category, L4F_ERROR, &
   'n. of levels in '//input_file//' differs from reference value '//&
   t2c(SIZE(volgridua%level))//'/'//t2c(SIZE(volgridz%level)))
  CALL raise_fatal_error()
ENDIF

IF (mask_file /= '') THEN
CALL read_input_volume(mask_file, volgridmask, volgridz)
! is this necessary for IMPLICIT allocation?
  ALLOCATE(intmask(SIZE(volgridmask%voldati,1),SIZE(volgridmask%voldati,2)))
  WHERE(c_e(volgridmask%voldati(:,:,1,1,1,1)))
    intmask = NINT(volgridmask%voldati(:,:,1,1,1,1))
  ELSEWHERE
    intmask = imiss
  END WHERE
  nzones = MAXVAL(intmask, mask=c_e(intmask))
ENDIF

! inputlev
CALL getarg(optind+3, input_file)
CALL read_input_volume(input_file, volgridlev, volgridz)
CALL l4f_category_log(category,L4F_INFO,'inputlev file: '//TRIM(input_file))

CALL getarg(optind+5, output_csv) ! +4?

IF (ldisplay) THEN
  PRINT*,'input volume >>>>>>>>>>>>>>>>>>>>'
  CALL display(volgridz)
  CALL display(volgridsurf)
  CALL display(volgridua)
  CALL display(volgridlev)
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
iqv = vartable_index(volgridua%var, 'B13001') !Specific Humidity
iqc = vartable_index(volgridua%var, 'B13192') !Cloud liquid water content
iqi = vartable_index(volgridua%var, 'B13193') !Cloud ice content
itd = vartable_index(volgridua%var, 'B12103') !Dew-Point Temperature
irh = vartable_index(volgridua%var, 'B13003') !Relative Humidity
iomega = vartable_index(volgridua%var, 'B11005')

ih = vartable_index(volgridlev%var, 'B10007') !Level high

it2 = vartable_index(volgridsurf%var, 'B12101') !Temperature/Dry-Bulb Temperature
itd2 = vartable_index(volgridsurf%var, 'B12102') !Wet-Bulb Temperature
iu10 = vartable_index(volgridsurf%var, 'B11003') !U-component
iv10 = vartable_index(volgridsurf%var, 'B11004') !V-component
itp = vartable_index(volgridsurf%var, 'B13011') !Total Precipitation/Total Water Equivalent


! associate vertical level indices to pseudo pressure levels

!l250 = 5
!l300 = 6
!l500 = 10
!l700 = 15
!l850 = 21
!l925 = 25
!l950 = 27
!l1000 = 32

l250 = 11
l300 = 13
l500 = 19
l700 = 26
l850 = 32
l925 = 37
l950 = 39
l1000 = 45


! layer depth

!IF (c_e(ih)) THEN
!   DO k = 1, (SIZE(volgridlev%level)-1)
!      depth(k) = volgridlev%voldati(1,1,k,1,1,ih)-volgridlev%voldati(1,1,k+1,1,1,ih)
!   ENDDO
!ELSE
!   CALL l4f_category_log(category, L4F_ERROR, &
!        'h missing')
!   CALL raise_fatal_error()
!ENDIF



! Omega 

IF (c_e(iw) .AND. c_e(ip) .AND. c_e(it)) THEN
   omega=omega_simple(volgridua%voldati(:,:,l300,1,1,it), &
        volgridua%voldati(:,:,l300,1,1,ip), &
        volgridua%voldati(:,:,l300,1,1,iw))
   IF (ALLOCATED(intmask)) THEN
     ! example_omega300=mask_average(omega, intmask,nzones)
      example_omega300 = mask_gt_threshold_less(omega, intmask, nzones, -0.5)
     ! IF (ALLOCATED(example_omega300)) THEN
     !    PRINT*,'% Vertical velocity at 300hPa (Pa/s)'
     !    PRINT*, example_omega300
     ! ENDIF
   ENDIF
   omega=omega_simple(volgridua%voldati(:,:,l700,1,1,it), &
        volgridua%voldati(:,:,l700,1,1,ip), &
        volgridua%voldati(:,:,l700,1,1,iw))
   IF (ALLOCATED(intmask)) THEN
      !example_omega=mask_average(omega, intmask,nzones)
      example_omega = mask_gt_threshold_less(omega, intmask, nzones, -0.5)
     ! IF (ALLOCATED(example_omega)) THEN
     !    PRINT*,'% Vertical velocity at 700hPa (Pa/s)'
     !    PRINT*, example_omega
     ! ENDIF
   ENDIF
ELSE
   CALL l4f_category_log(category, L4F_ERROR, &
        'w, p or t missing in upper air data')
   CALL raise_fatal_error()
ENDIF

! Metto un valore minimo di qv dove manca
WHERE (.NOT.c_e(volgridua%voldati(:,:,:,1,1,iqv)) .OR. volgridua%voldati(:,:,:,1,1,iqv) < 1.E-7)
   volgridua%voldati(:,:,:,1,1,iqv) = 1.E-5
END WHERE

!ALLOCATE (dew(SIZE(volgridua%voldati,1),SIZE(volgridua%voldati,2)))

DO f = 1, SIZE(volgridua%level)
   volgridua%voldati(:,:,f,1,1,itd) = td_pq((volgridua%voldati(:,:,f,1,1,ip)/100), &
        volgridua%voldati(:,:,f,1,1,iqv))
ENDDO


!Vertical Wind Shear 

IF (c_e(iu) .AND. c_e(iv)) THEN
   vertwindshear = vertical_shear(volgridua%voldati(:,:,l500,1,1,iu), &
        volgridua%voldati(:,:,l500,1,1,iv),volgridua%voldati(:,:,l950,1,1,iu), &
        volgridua%voldati(:,:,l950,1,1,iv))
   IF (ALLOCATED(intmask)) THEN
      !example_vws500 = mask_average(vertwindshear, intmask, nzones)
      example_vws500 = mask_gt_threshold(vertwindshear, intmask, nzones, 15.)
     ! IF (ALLOCATED(example_vws500)) THEN
     !    PRINT*,'% Vertical wind shear 500-950 hPa (m/s)'
     !    PRINT*, example_vws500
     ! ENDIF
   ENDIF
   vertwindshear = vertical_shear(volgridua%voldati(:,:,l700,1,1,iu), &
        volgridua%voldati(:,:,l700,1,1,iv),volgridua%voldati(:,:,l1000,1,1,iu), &
        volgridua%voldati(:,:,l1000,1,1,iv))
   IF (ALLOCATED(intmask)) THEN
      !example_vws = mask_average(vertwindshear, intmask, nzones)
      example_vws = mask_gt_threshold(vertwindshear, intmask, nzones, 10.)
     ! IF (ALLOCATED(example_vws)) THEN
     !    PRINT*,'% Vertical wind shear 700-1000 hPa (m/s)'
     !    PRINT*, example_vws
     ! ENDIF
   ENDIF
ELSE
   CALL l4f_category_log(category, L4F_ERROR, &
        'u or v missing in upper air data')
   CALL raise_fatal_error()
ENDIF


! Jet

IF (c_e(iu) .AND. c_e(iv)) THEN
   jets = jet(volgridua%voldati(:,:,l925,1,1,iu),volgridua%voldati(:,:,l925,1,1,iv))
   IF (ALLOCATED(intmask)) THEN
      example_jet925 = mask_average(jets, intmask, nzones)
     ! IF (ALLOCATED(example_jet925)) THEN
     !    PRINT*, 'Jet at 925 hPa (m/s)'
     !    PRINT*, example_jet925
     ! ENDIF
   ENDIF
   jets = jet(volgridua%voldati(:,:,l250,1,1,iu),volgridua%voldati(:,:,l250,1,1,iv))
   IF (ALLOCATED(intmask)) THEN
      example_jet = mask_average(jets, intmask, nzones)
     ! IF (ALLOCATED(example_jet)) THEN
     !    PRINT*, 'Jet at 250 hPa (m/s)'
     !    PRINT*, example_jet
     ! ENDIF
   ENDIF
ELSE
   CALL l4f_category_log(category, L4F_ERROR, &
        'u or v missing in upper air data')
   CALL raise_fatal_error()
ENDIF

!volumeinv=volgridua%voldati(:,:,size(volgridua%voldati,3):1:-1,:,:,:)
ALLOCATE(volumeinv(size(volgridua%voldati,3),size(volgridua%voldati,1), &
     size(volgridua%voldati,2), size(volgridua%voldati,4), &
     size(volgridua%voldati,5),size(volgridua%voldati,6)))
DO k = 1, size(volgridua%voldati,3)
   volumeinv(k,:,:,:,:,:) = volgridua%voldati(:,:,size(volgridua%voldati,3)-k+1,:,:,:)
ENDDO

volumeinv(:,:,:,:,:,ip)=volumeinv(:,:,:,:,:,ip)/100


! Surface Lifted Index

NT=size(volgridua%voldati,3)

ALLOCATE (SLI(SIZE(volgridua%voldati,1),SIZE(volgridua%voldati,2)))

IF (c_e(ip) .AND. c_e(it) .AND. c_e(itd)) THEN
   
   DO j = 1 , size(volgridua%voldati,2)
      DO i = 1 , size(volgridua%voldati,1)
         SLI(i,j) = SI_LI(volumeinv(:,i,j,1,1,ip),volumeinv(:,i,j,1,1,it), &
              volumeinv(:,i,j,1,1,itd), NT)     
      ENDDO
   ENDDO
   IF (ALLOCATED(intmask)) THEN
      example_sli = mask_average(SLI, intmask, nzones)
   ENDIF
ELSE
   CALL l4f_category_log(category, L4F_ERROR, &
        'p, t or td missing in upper air data')
   CALL raise_fatal_error()
ENDIF

!IF (ALLOCATED(example_sli)) THEN
!   PRINT*,'Lifted Index (SLI)'
!   PRINT*,example_sli
!ENDIF

! CIN

ALLOCATE (inibition(SIZE(volgridua%voldati,1),SIZE(volgridua%voldati,2)))

IF (c_e(ip) .AND. c_e(it) .AND. c_e(itd)) THEN

   DO j = 1 , size(volgridua%voldati,2)
      DO i = 1 , size(volgridua%voldati,1)
         inibition(i,j) = cin((volumeinv(:,i,j,1,1,ip)), &
              volumeinv(:,i,j,1,1,it), &
              volumeinv(:,i,j,1,1,itd), NT)
      ENDDO
   ENDDO
   IF (ALLOCATED(intmask)) THEN
      example_cin = mask_average(inibition, intmask, nzones)
   ENDIF
ELSE
   CALL l4f_category_log(category, L4F_ERROR, &
        'p, t or td missing in upper air data')
   CALL raise_fatal_error()
ENDIF

! K index

IF (c_e(it) .AND. c_e(itd)) THEN

   kind = kindex(volgridua%voldati(:,:,l850, 1, 1, it), &
        volgridua%voldati(:,:,l500, 1, 1, it), &
        volgridua%voldati(:,:,l850, 1, 1, itd), &
        volgridua%voldati(:,:,l700, 1, 1, it), &
        volgridua%voldati(:,:,l700, 1, 1, itd))
   IF (ALLOCATED(intmask)) THEN
      example_kindex = mask_average(kind, intmask, nzones)
   ENDIF
ELSE
   CALL l4f_category_log(category, L4F_ERROR, &
        't or td missing in upper air data')
   CALL raise_fatal_error()
ENDIF


! MCS index

IF (c_e(iu) .AND. c_e(iv) .AND. c_e(it)) THEN
   tadvection700 = grid_horiz_advection( &
        volgridua%voldati(:,:,l700, 1, 1, iu), &
        volgridua%voldati(:,:,l700, 1, 1, iv), &
        volgridua%voldati(:,:,l700, 1, 1, it), &
        dxm1, dym1)

   mcsind = mcs(SLI,vertwindshear,tadvection700)
   IF (ALLOCATED(intmask)) THEN
    !  example_mcsindex = mask_average(mcsind, intmask, nzones)
      example_mcsindex = mask_gt_threshold(mcsind, intmask, nzones, 0.)
   ENDIF
ELSE
   CALL l4f_category_log(category, L4F_ERROR, &
        'u, v or t missing in upper air data')
   CALL raise_fatal_error()
ENDIF


! Total Water Content

ALLOCATE(virtualt(SIZE(volgridua%voldati,1),SIZE(volgridua%voldati,2)))
ALLOCATE(wclayer(SIZE(volgridua%voldati,1),SIZE(volgridua%voldati,2)))
ALLOCATE(totalwc(SIZE(volgridua%voldati,1),SIZE(volgridua%voldati,2)))

!k = (SIZE(volgridua%voldati,3)-SIZE(volgridlev%level)+2) ! k=9
!d = 1
!DO WHILE (k<(SIZE(volgridua%voldati,3)+1)) ! k<33
!   DO j = 1 , size(volgridua%voldati,2)
!      DO i = 1 , size(volgridua%voldati,1)
!         virtualt(i,j) = TVIR(volgridua%voldati(i,j,k,1,1,itd), &
!              volgridua%voldati(i,j,k,1,1,it), &
!              volgridua%voldati(i,j,k,1,1,ip))
!         wclayer(i,j) = TWC(volgridua%voldati(i,j,k,1,1,ip),virtualt(i,j), &
!              volgridua%voldati(i,j,k,1,1,iqv), &
!              (volgridlev%voldati(1,1,d,1,1,ih)-volgridlev%voldati(1,1,d+1,1,1,ih)))
!         IF (K>(SIZE(volgridua%voldati,3)-SIZE(volgridlev%level)+2)) THEN
!            totalwc(i,j) = totalwc(i,j) + wclayer(i,j)
!         ELSE
!            totalwc(i,j) = wclayer(i,j)
!         ENDIF
!      ENDDO
!   ENDDO
!   k=k+1
!   d=d+1
!ENDDO

k = 1
DO WHILE (k<(SIZE(volgridua%voldati,3)+1)) ! k<46
   DO j = 1 , size(volgridua%voldati,2)
      DO i = 1 , size(volgridua%voldati,1)
         virtualt(i,j) = TVIR(volgridua%voldati(i,j,k,1,1,itd), &
              volgridua%voldati(i,j,k,1,1,it), &
              volgridua%voldati(i,j,k,1,1,ip))
         wclayer(i,j) = TWC(volgridua%voldati(i,j,k,1,1,ip),virtualt(i,j), &
              volgridua%voldati(i,j,k,1,1,iqv), &
              (volgridlev%voldati(1,1,k,1,1,ih)-volgridlev%voldati(1,1,k+1,1,1,ih)))
         IF (K>1) THEN
            totalwc(i,j) = totalwc(i,j) + wclayer(i,j)
         ELSE
            totalwc(i,j) = wclayer(i,j)
         ENDIF
      ENDDO
   ENDDO
   k=k+1
ENDDO


IF (ALLOCATED(intmask)) THEN
   example_twc = mask_average(totalwc, intmask, nzones)
ENDIF

! Relative humidity
   
relhum = volgridua%voldati(:,:,l500, 1, 1, irh)
IF (ALLOCATED(intmask)) THEN
   example_relhum = mask_average(relhum, intmask, nzones)
ENDIF

! LFC

!W(TD,PT) mixing ratio [g/Kg]
!!!!!!!!!!!!!!!!!!!!!!!!!!! SPESSORE STRATO IMPOSTATO MANUALMENTE IN alfc
ALLOCATE(lfc(SIZE(volgridua%voldati,1),SIZE(volgridua%voldati,2)))
DO j = 1 , size(volgridua%voldati,2)
   DO i = 1 , size(volgridua%voldati,1)
      
      k=size(volgridua%voldati,3)
      
      DO WHILE (k .GT. 0)
         aw(size(volgridua%voldati,3)-k+1) = W(volgridua%voldati(i,j,k,1,1,itd), &
              (volgridua%voldati(i,j,k,1,1,ip)/100))
         k=k-1
      ENDDO
      lfc(i,j)=alfc(volumeinv(:,i,j,1,1,ip),volumeinv(:,i,j,1,1,it),aw,NT,50.)
   ENDDO
ENDDO

IF (ALLOCATED(intmask)) THEN
   example_lfc = mask_average(lfc, intmask, nzones)
ENDIF

! LCL
! Equilibrium Level

! AEQL(pt,tt,nt,plow,ths)

ALLOCATE(templcl(SIZE(volgridua%voldati,1),SIZE(volgridua%voldati,2)))
ALLOCATE(teta(SIZE(volgridua%voldati,1),SIZE(volgridua%voldati,2)))
ALLOCATE(lcl(SIZE(volgridua%voldati,1),SIZE(volgridua%voldati,2)))
ALLOCATE(tamblcl(SIZE(volgridua%voldati,1),SIZE(volgridua%voldati,2)))
ALLOCATE(tsatlcl(SIZE(volgridua%voldati,1),SIZE(volgridua%voldati,2)))
ALLOCATE(equilibrium(SIZE(volgridua%voldati,1),SIZE(volgridua%voldati,2)))
ALLOCATE(samix(SIZE(volgridua%voldati,1),SIZE(volgridua%voldati,2)))




DO j = 1 , size(volgridua%voldati,2)
   DO i = 1 , size(volgridua%voldati,1)
      
      k=size(volgridua%voldati,3)
      
      DO WHILE (k .GT. 0)
         aw(size(volgridua%voldati,3)-k+1) = W(volgridua%voldati(i,j,k,1,1,itd), &
              (volgridua%voldati(i,j,k,1,1,ip)/100))
         k=k-1
      ENDDO
      
      templcl(i,j) = TL(volgridua%voldati(i,j,size(volgridua%voldati,3),1,1,it), &
           volgridua%voldati(i,j,size(volgridua%voldati,3),1,1,itd))
      
      teta(i,j) = O(volgridua%voldati(i,j,size(volgridua%voldati,3),1,1,it), &
           (volgridua%voldati(i,j,size(volgridua%voldati,3),1,1,ip)/100))
      
      samix(i,j) = PKAVR(volumeinv(:,i,j,1,1,ip), &
           aw,NT,volgridua%voldati(i,j,size(volgridua%voldati,3),1,1,ip)/100, &
           ((volgridua%voldati(i,j,size(volgridua%voldati,3),1,1,ip)/100)-50))
      
      lcl(i,j) = ALCLM(samix(i,j),teta(i,j), &
           (volgridua%voldati(i,j,size(volgridua%voldati,3),1,1,ip)/100))
      
   ENDDO
ENDDO


IF (ALLOCATED(intmask)) THEN
   example_lcl = mask_average(lcl, intmask, nzones)
ENDIF

!calcolo ccl che mi permette di avere il livello di equilibrio

ALLOCATE(pccl(SIZE(volgridua%voldati,1),SIZE(volgridua%voldati,2)))
ALLOCATE(tambccl(SIZE(volgridua%voldati,1),SIZE(volgridua%voldati,2)))
ALLOCATE(tsatccl(SIZE(volgridua%voldati,1),SIZE(volgridua%voldati,2)))
ALLOCATE(equilibriumccl(SIZE(volgridua%voldati,1),SIZE(volgridua%voldati,2)))

DO j = 1 , size(volgridua%voldati,2)
   DO i = 1 , size(volgridua%voldati,1)
      
      k=size(volgridua%voldati,3)
      
      DO WHILE (k .GT. 0)
         aw(size(volgridua%voldati,3)-k+1) = W(volgridua%voldati(i,j,k,1,1,itd), &
              (volgridua%voldati(i,j,k,1,1,ip)/100))
         k=k-1
      ENDDO
      samix(i,j) = PKAVR(volumeinv(:,i,j,1,1,ip), &
           aw,NT,volgridua%voldati(i,j,size(volgridua%voldati,3),1,1,ip)/100, &
           ((volgridua%voldati(i,j,size(volgridua%voldati,3),1,1,ip)/100)-50))
      
      pccl(i,j)= CCL((volgridua%voldati(i,j,size(volgridua%voldati,3),1,1,ip)/100)-50, &
           volumeinv(:,i,j,1,1,ip),volumeinv(:,i,j,1,1,it),volumeinv(:,i,j,1,1,itd), &
           samix(i,j), NT)
      
      tambccl(i,j) = tmr(samix(i,j),pccl(i,j))
      tsatccl(i,j) = OS(tambccl(i,j),pccl(i,j))
      IF (lfc(i,j) == 0) THEN
         equilibriumccl(i,j) = 0.
      ELSE
         equilibriumccl(i,j)=AEQL(volumeinv(:,i,j,1,1,ip),volumeinv(:,i,j,1,1,it), &
              NT,pccl(i,j),tsatccl(i,j))  
      ENDIF
   ENDDO
ENDDO


IF (ALLOCATED(intmask)) THEN
   example_equilibrium = mask_average(equilibriumccl, intmask, nzones)
   example_ccl = mask_average(pccl, intmask, nzones)
ENDIF

! cape(pt,tt,np,lfc,eql)

ALLOCATE(capeindex(SIZE(volgridua%voldati,1),SIZE(volgridua%voldati,2)))

DO j = 1 , size(volgridua%voldati,2)
   DO i = 1 , size(volgridua%voldati,1)
      IF (lfc(i,j) == 0) THEN
         capeindex(i,j) = 0.
      ELSE
         capeindex(i,j)=cape(volumeinv(:,i,j,1,1,ip),volumeinv(:,i,j,1,1,it), &
              NT,lfc(i,j),equilibriumccl(i,j))
      ENDIF
   ENDDO
ENDDO
IF (ALLOCATED(intmask)) THEN
   example_cape = mask_average(capeindex, intmask, nzones)
ENDIF

! Temperature and Temperature Dew Point Advection at 500hPa in 12h

tadvection12h = (grid_horiz_advection( &
     volgridua%voldati(:,:,l500, 1, 1, iu), &
     volgridua%voldati(:,:,l500, 1, 1, iv), &
     volgridua%voldati(:,:,l500, 1, 1, it), &
     dxm1, dym1)*43200)

tdewadvection12h = (grid_horiz_advection( &
     volgridua%voldati(:,:,l850, 1, 1, iu), &
     volgridua%voldati(:,:,l850, 1, 1, iv), &
     volgridua%voldati(:,:,l850, 1, 1, itd), &
     dxm1, dym1)*43200)


ALLOCATE(hgeopotential500(SIZE(volgridua%voldati,1),SIZE(volgridua%voldati,2)))

DO j = 1 , size(volgridua%voldati,2)
   DO i = 1 , size(volgridua%voldati,1)
      hgeopotential500(i,j) = Z(500., &
           volumeinv(:,i,j,1,1,ip),volumeinv(:,i,j,1,1,it),volumeinv(:,i,j,1,1,itd),NT)       
   ENDDO
ENDDO
      
      geopadvection12h = (grid_horiz_advection( &
           volgridua%voldati(:,:,l500, 1, 1, iu), &
           volgridua%voldati(:,:,l500, 1, 1, iv), &
           hgeopotential500, dxm1, dym1)*4320)
! hgeopotential500 è stato calcolato in conrrispondenza della 500hPa,
! mentre u e v fanno riferimento al livello più vicino alla 500hPa 

IF (ALLOCATED(intmask)) THEN
  ! example_tadv12h = mask_average(tadvection12h, intmask, nzones)
   example_tadv12h = mask_gt_threshold_less(tadvection12h, intmask, nzones, -5.)
  ! example_tdewadv12h = mask_average(tdewadvection12h, intmask, nzones)
   example_tdewadv12h = mask_gt_threshold(tdewadvection12h, intmask, nzones, 10.)
   example_hgeopot = mask_average(hgeopotential500, intmask, nzones)
   example_geopavv = mask_average(geopadvection12h, intmask, nzones)
ENDIF


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! per indici che hanno bisogno dei contact point tra i parametri di input procedere così
!DO j=,size(volgridua%voldati,2)
!   DO i=,size(volgridua%voldati,1)
!      res(i,j) = (volgridua%voldati(:,i,j,1,1,ip),volgridua%voldati(:,i,j,1,1,it),  ..., size(volgridua%voldati,3), 
!   ENDDO
!ENDDO


! examples opf grid computations
! remember voldati(x,y,level,time,timerange,var)
IF (c_e(iu) .AND. c_e(iv) .AND. c_e(it)) THEN
   
   vorticity = grid_vorticity_z( &
        volgridua%voldati(:,:,l500, 1, 1, iu), &
        volgridua%voldati(:,:,l500, 1, 1, iv), &
        dxm1, dym1)*1.E5
   
   
! example of spatialization with mask
   IF (ALLOCATED(intmask)) THEN
! implicit allocation
      ! example_index1 = mask_average(vorticity, intmask, nzones)
      !VORTICITY as percentage over threshold
      example_index1 = mask_gt_threshold(vorticity, intmask, nzones, 0.)
   ENDIF
ELSE
   CALL l4f_category_log(category, L4F_ERROR, &
        'u, v or t missing in upper air data')
   CALL raise_fatal_error()
ENDIF


CALL getval(volgridua%time(1), simpledate=filetimename)

OPEN(unit=2,file=output_csv,position="append")
WRITE(2,'(19(A,'',''))')"Data","Macroarea","%VV300","%VV700", &
 "%VWS500.950","%VWS700.1000","Jet925","Jet250", &
 "LI","CAPE","CIN","Kindex","%MCSindex","TWC","R.H.500", &
 "%AvvT500","%AvvTd850","%VortRel500","AvvGeop500"

DO i = 1, nzones
  CALL init(csv_writer)
  CALL csv_record_addfield(csv_writer, filetimename(1:10))
  CALL csv_record_addfield(csv_writer, i)
  CALL csv_record_addfield(csv_writer, example_omega300(i))
  CALL csv_record_addfield(csv_writer, example_omega(i))
  CALL csv_record_addfield(csv_writer, example_vws500(i))
  CALL csv_record_addfield(csv_writer, example_vws(i))
  CALL csv_record_addfield(csv_writer, example_jet925(i))
  CALL csv_record_addfield(csv_writer, example_jet(i))
  CALL csv_record_addfield(csv_writer, example_sli(i))
  CALL csv_record_addfield(csv_writer, example_cape(i))
  CALL csv_record_addfield(csv_writer, example_cin(i))
  CALL csv_record_addfield(csv_writer, example_kindex(i))
  CALL csv_record_addfield(csv_writer, example_mcsindex(i))
  CALL csv_record_addfield(csv_writer, example_twc(i))
  CALL csv_record_addfield(csv_writer, example_relhum(i))
  CALL csv_record_addfield(csv_writer, example_tadv12h(i))
  CALL csv_record_addfield(csv_writer, example_tdewadv12h(i))
  CALL csv_record_addfield(csv_writer, example_index1(i))
  CALL csv_record_addfield(csv_writer, example_geopavv(i))
  WRITE(2,'(A)')csv_record_getrecord(csv_writer)
  CALL delete(csv_writer)
ENDDO

CLOSE(2)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CALL delete(volgridz)
CALL delete(volgridsurf)
CALL delete(volgridua)
CALL delete(volgridlev)

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



END PROGRAM prodsim_thunderstorm_index
