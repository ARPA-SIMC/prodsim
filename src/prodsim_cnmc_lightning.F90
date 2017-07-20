! file in:
! /autofs/webarchive/prodotti_infomet/lampi
!
! unpacking suggestion:
! tar -x --transform='s?.*/??g' -v \
!  -f /autofs/webarchive/prodotti_infomet/lampi/LAMPI_20170712_ITALIA.tar.gz
!
! run example:
! prodsim_cnmc_lightning --comp-frac-valid=0.9 \
!  --comp-start=2017-07-12 --comp-STOP=2017-07-13 --comp-step='0 01' \
!  --file-TEMPLATE='CNMC_LAM_%d_ITALIA_SF@@@@@@_@@@@@@@@@@@@_@@@_005_@@@@.TXT' \
!  --coord-file=/scratch/dcesari/gis/vector/macroaree_italia/aree_allertamento_ll \
!  out.csv
PROGRAM prodsim_cnmc_lightning
USE log4fortran
USE err_handling
USE missing_values
USE optionparser_class
USE georef_coord_class
USE datetime_class
USE cnmc_lightning_mo
IMPLICIT NONE

INTEGER :: category, ier
CHARACTER(len=512) :: a_name, coord_file, file_template, input_file, output_file, line
CHARACTER(len=12) :: filetimename
CHARACTER(len=8) :: comp_algo
TYPE(optionparser) :: opt
INTEGER :: optind, optstatus
LOGICAL :: version
CHARACTER(len=23) :: comp_start, comp_stop, comp_step, file_step
TYPE(datetime):: c_sta, c_sto, comptime, filetime
TYPE(timedelta) :: c_ste, f_ste, validity_interval
TYPE(arrayof_georef_coord_array) :: poly
REAL :: comp_frac_valid
LOGICAL :: center_time
INTEGER :: i, ds, ierr, t_sec, v_sec
TYPE(cnmc_lightning_t) :: lightning
REAL,ALLOCATABLE :: val(:)


! initialise logging
CALL l4f_launcher(a_name,a_name_force='prodsim_vg6d_tcorr')
ier=l4f_init()
! set a_name
category=l4f_category_get(a_name//'.main')

! define the option parser
opt = optionparser_new(description_msg= &
 'Aggregate pointwise lightning data in discrete time intervals and in space (polygons).', &
 usage_msg='Usage: prodsim_cnmc_lightning [options] output')

CALL optionparser_add(opt, ' ', 'coord-file', coord_file, default='', &
 help='file in shp format with polygons for spatial aggregation of the results')

CALL optionparser_add(opt, ' ', 'comp-algo', comp_algo, 'count', &
 help='type of aggregation algorithm, ''count'' for counting, &
 &''sum'' for summing ''max'' for taking the maximum of all values in &
 &space and time, ''logsum''??')

CALL optionparser_add(opt, ' ', 'comp-step', comp_step, '0000000000 01:00:00.000', help= &
 'length of time aggregation step in the format &
 &''YYYYMMDD hh:mm:ss.msc'', it can be simplified up to the form ''D hh''')

CALL optionparser_add(opt, ' ', 'comp-start', comp_start, '', help= &
 'start of processing time interval &
 &in the format iso ''YYYY-MM-DD hh:mm:ss.msc'' where characters on the right are optional')

CALL optionparser_add(opt, ' ', 'comp-stop', comp_stop, '', help= &
 'end of processing time interval &
 &in the format iso ''YYYY-MM-DD hh:mm:ss.msc'' where characters on the right are optional')

CALL optionparser_add(opt, ' ', 'comp-frac-valid', comp_frac_valid, 1., help= &
 '(from 0. to 1.) fraction of input data that has to be valid in order to consider a &
 &statistically processed value acceptable')

CALL optionparser_add(opt, ' ', 'file-template', file_template, '', help= &
 'template for determining input file names, &
 &the special string "%d" will be translated into YYYYMMDDHHMM, &
 &indicating end of time interval of the searched file')

CALL optionparser_add(opt, ' ', 'file-step', file_step, '0000000000 00:05:00.000', help= &
 'length of time aggregation step in each file &
 &in the format ''YYYYMMDD hh:mm:ss.msc'', it can be simplified up to the form ''D hh''')

CALL optionparser_add(opt, ' ', 'center-time', center_time, help= &
 'set the time of output data at the center of the processing interval &
 &rather than at the end of the interval as is usually done for observations')

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

ds = INDEX(file_template, '%d')
IF (ds <= 0) THEN
  CALL l4f_category_log(category,L4F_ERROR,'file-template does not contain date format ''%d''')
  CALL raise_fatal_error()
ENDIF
c_sta = datetime_new(isodate=comp_start)
IF (.NOT.c_e(c_sta)) THEN
  CALL l4f_category_log(category,L4F_ERROR,'argument comp-start '//TRIM(comp_start)//' not provided or not valid')
  CALL raise_fatal_error()
ENDIF
c_sto = datetime_new(isodate=comp_stop)
IF (.NOT.c_e(c_sto)) THEN
  CALL l4f_category_log(category,L4F_ERROR,'argument comp-stop '//TRIM(comp_stop)//' not provided or not valid')
  CALL raise_fatal_error()
ENDIF
c_ste = timedelta_new(isodate=comp_step)
IF (.NOT.c_e(c_ste) .OR. c_ste <= timedelta_0) THEN
  CALL l4f_category_log(category,L4F_ERROR,'argument comp-step '//TRIM(comp_step)//' not provided or not valid')
  CALL raise_fatal_error()
ENDIF
f_ste = timedelta_new(isodate=file_step)
IF (.NOT.c_e(f_ste) .OR. f_ste <= timedelta_0) THEN
  CALL l4f_category_log(category,L4F_ERROR,'argument file-step '//TRIM(file_step)//' not provided or not valid')
  CALL raise_fatal_error()
ENDIF

CALL import(poly, coord_file)
IF (poly%arraysize <= 0) THEN
  CALL l4f_category_log(category, L4F_ERROR, &
   'error importing shapefile '//TRIM(coord_file))
  CALL raise_fatal_error()
ENDIF
IF (optind > iargc()) THEN
  CALL optionparser_printhelp(opt)
  CALL l4f_category_log(category,L4F_ERROR,'output file missing')
  CALL raise_fatal_error()
  CALL EXIT(1)
ENDIF

CALL getarg(optind, output_file)
IF (output_file == '-') THEN
  CALL l4f_category_log(category, L4F_INFO, 'trying /dev/stdout as stdout unit')
  output_file = '/dev/stdout'
ELSE
  CALL l4f_category_log(category,L4F_INFO,'output file: '//TRIM(output_file))
ENDIF
OPEN(10, file=output_file)

CALL getval(c_ste, asec=t_sec)
ALLOCATE(val(poly%arraysize))
comptime = c_sta + c_ste
loop_comptime: DO WHILE(comptime <= c_sto)
  val(:) = 0.
  validity_interval = timedelta_0
  filetime = comptime - c_ste + f_ste
  loop_filetime: DO WHILE(filetime <= comptime)
    CALL getval(filetime, simpledate=filetimename)
    input_file = file_template(:ds-1)//filetimename//file_template(ds+2:)
    OPEN(11, file=input_file, status='OLD', iostat=ierr)
    IF (ierr /= 0) THEN
! a file is missing, skip
      filetime = filetime + f_ste
      CYCLE loop_filetime
    ENDIF

    loop_fileline: DO WHILE(.TRUE.)
      READ(11,'(A)',iostat=ierr)line
      IF (ierr /= 0) EXIT loop_fileline
      CALL lightning%readline(line)
      IF (lightning%c_e()) THEN ! process the datum
        DO i = 1, poly%arraysize
          IF (inside(lightning%coord, poly%array(i))) THEN
            val(i) = val(i) + 1.
            EXIT
          ENDIF
        ENDDO
      ENDIF
    ENDDO loop_fileline
    CLOSE(11)
! increment counters
    validity_interval = validity_interval + f_ste
    filetime = filetime + f_ste
  ENDDO loop_filetime
! check whether we have enough data
  CALL getval(validity_interval, asec=v_sec)
  IF (REAL(v_sec)/REAL(t_sec) >= comp_frac_valid) THEN
    IF (center_time) THEN
      CALL getval(comptime-c_ste/2, simpledate=filetimename)
    ELSE
      CALL getval(comptime, simpledate=filetimename)
    ENDIF
    
    DO i = 1, poly%arraysize
      WRITE(10,'(A,'','',I0,'','',I0)')TRIM(filetimename),i,INT(val(i))
    ENDDO
  ENDIF
! increment counter
  comptime = comptime + c_ste
ENDDO loop_comptime

CLOSE(10)

END PROGRAM prodsim_cnmc_lightning
