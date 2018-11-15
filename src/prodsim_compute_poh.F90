! gfortran -I/usr/lib64/gfortran/modules -o etop2grib etop2grib.f90 -lnetcdff -lnetcdf -llog4fortran -lsim_base -lsim_vol7d -lsim_volgrid6d

! Per ottenere i dati di zero termico la query consigliata per un DS
! di Cosmo grib1 e`:

! arki-query --data 'reftime: >= 2018-09-01 00:00:00, <= 2018-09-02 00:00:00; product:GRIB1,200,201,84; level:GRIB1,4; timerange:GRIB1,0,1h,0h or GRIB1,0,2h,0h or GRIB1,0,3h,0h or GRIB1,0,4h,0h or GRIB1,0,5h,0h or GRIB1,0,6h,0h or GRIB1,0,7h,0h or GRIB1,0,8h,0h or GRIB1,0,9h,0h or GRIB1,0,10h,0h or GRIB1,0,11h,0h or GRIB1,0,12h,0h' http://arkimet.metarpa:8090/dataset/lmsmr4x52 | \
! vg6d_transform --trans-mode=s --trans-type=zoom --sub-type=coordbb --ilon=8.4 --ilat=43.3 --flon=13.3 --flat=46.1 - iso0.grib

! aggiustare reftime secondo necessita`, la vg6d_trasform riduce i
! dati per permettere l'elaborazione di intervalli piu` lunghi in
! un'esecuzione sola, usare coordbb e non coord per assicurarsi di
! includere tutta l'area radar.

PROGRAM prodsim_compute_poh
USE log4fortran
USE missing_values
USE err_handling
USE optionparser_class
USE datetime_class
USE geo_coord_class
USE volgrid6d_class
USE grid_transform_class
USE simc_radar_mo
IMPLICIT NONE

! logging
INTEGER :: category
INTEGER :: ier
CHARACTER(len=512):: a_name
! radar
TYPE(radar_desc) :: radarlist(2)
! options
TYPE(geo_coord) :: tmpcoord
TYPE(optionparser) :: opt
INTEGER :: optind, optstatus
INTEGER :: iargc
CHARACTER(len=23) :: start_date, end_date, comp_step
TYPE(timedelta) :: c_i
LOGICAL :: ldisplay, version
CHARACTER(len=512):: input_file, output_template, output_file

INTEGER :: i
TYPE(volgrid6d),POINTER :: voltmp(:)
TYPE(volgrid6d) :: vol_iso0, vol_etop, vol_iso0_rad
TYPE(transform_def) :: mod2rad

! get launcher name
CALL l4f_launcher(a_name, a_name_force="etop2grib")
! log4fortran init
ier = l4f_init()
! set a_name
category = l4f_category_get(a_name//".main")

! initialise radarlist
CALL init(tmpcoord, 10.508333D0, 44.791667D0)
radarlist(1) = radar_desc(coord=tmpcoord, radius=125000.D0, name='RAD:IYai,PLC:itgat')
CALL init(tmpcoord, 11.623611D0, 44.654722D0)
radarlist(2) = radar_desc(coord=tmpcoord, radius=125000.D0, name='RAD:IY46,PLC:itspc')
CALL init(vol_iso0)
CALL init(vol_etop)

! define the option parser
opt = optionparser_new(description_msg= &
 'Compute Probability Of Hail from radar data in netcdf and 0deg isotherm height &
& model data for a time interval.', &
 usage_msg='Usage: prodsim_compute_poh [options] inputgrib inputnetcdf1 [inputnetcdf2...] outputfile')

CALL optionparser_add(opt, 's', 'start-date', start_date, '', help= &
 'initial date for processing data &
 &with format iso YYYY-MM-DD hh:mm:ss.msc where characters on the right are optional, &
 &the default is the first time level available in the 0deg isotherm height grib file')
CALL optionparser_add(opt, 'e', 'end-date', end_date, '', help= &
 'final date for processing data &
 &with format iso YYYY-MM-DD hh:mm:ss.msc where characters on the right are optional, &
 &the default is the last time level available in the 0deg isotherm height grib file')
CALL optionparser_add(opt, ' ', 'comp-step', comp_step, '0000000000 01:00:00.000', help= &
 'length of computation step in the format &
 &''YYYYMMDD hh:mm:ss.msc'', it can be simplified up to the form ''D hh''')
CALL optionparser_add(opt, ' ', 'output-template', output_template, '', help= &
 'template for defining output grid, default is the grid of the radar netcdf files')
CALL optionparser_add(opt, 'd', 'display', ldisplay, help= &
 'briefly display the data volume imported, warning: this option is incompatible &
 &with output on stdout.')

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
  WRITE(*,'(A,1X,A)')'prodsim_compute_poh','0.1' !VERSION
  CALL exit(0)
ENDIF

! check input/output files
i = iargc() - optind
IF (i < 0) THEN
  CALL optionparser_printhelp(opt)
  CALL l4f_category_log(category,L4F_ERROR,'input grib file missing')
  CALL raise_fatal_error()
ELSE IF (i < 1) THEN
  CALL optionparser_printhelp(opt)
  CALL l4f_category_log(category,L4F_ERROR,'input netcdf file(s) missing')
  CALL raise_fatal_error()
ELSE IF (i < 2) THEN
  CALL optionparser_printhelp(opt)
  CALL l4f_category_log(category,L4F_ERROR,'output file missing')
  CALL raise_fatal_error()
ENDIF

CALL getarg(optind, input_file)

! read 0deg isotherm grib file
CALL l4f_category_log(category,L4F_INFO,'Importing grib file '//TRIM(input_file))
! error checking
CALL IMPORT(voltmp, input_file, time_definition=1) ! verification time
IF (.NOT.ASSOCIATED(voltmp)) THEN
  CALL l4f_category_log(category,L4F_ERROR,'Error importing grib file '//TRIM(input_file))
  CALL raise_fatal_error()
ENDIF
IF (SIZE(voltmp) /= 1) THEN
  CALL l4f_category_log(category,L4F_ERROR,'Error, grib file has '//t2c(SIZE(voltmp))//' different grids, this is not allowed')
  CALL raise_fatal_error()
ENDIF
! flatten timeranges
CALL rounding(voltmp(1), vol_iso0, timerange=voltmp(1)%timerange)
CALL delete(voltmp)
! define time interval
IF (SIZE(vol_iso0%time) > 1) THEN
  c_i = vol_iso0%time(2) - vol_iso0%time(1)
ELSE
  c_i = timedelta_new(hour=1)
ENDIF
DO i = 2, SIZE(vol_iso0%time)
  IF (vol_iso0%time(i) - vol_iso0%time(i-1) /= c_i) THEN
    CALL l4f_category_log(category,L4F_ERROR,'Error, grib file has non regular time steps, this is not allowed')
    CALL raise_fatal_error()
  ENDIF
ENDDO

IF (ldisplay) CALL display(vol_iso0)

DO i = optind + 1, iargc() - 1
  CALL getarg(i, input_file)
  CALL read_netcdf(input_file, vol_etop, vol_iso0%time, c_i, radarlist, category)
ENDDO

CALL init(mod2rad, 'inter', 'bilin')
CALL transform(mod2rad, vol_etop%griddim, vol_iso0, vol_iso0_rad, &
 categoryappend='mod2rad', clone=.TRUE., decode=.TRUE.)
! free memory
CALL delete(vol_iso0)
IF (ldisplay) CALL display(vol_iso0_rad)
PRINT*,SHAPE(vol_iso0_rad%voldati)
PRINT*,SHAPE(vol_etop%voldati)
vol_iso0_rad%voldati = compute_poh_2(vol_etop%voldati, vol_iso0_rad%voldati)
PRINT*,'dati validi:',COUNT(c_e(vol_iso0_rad%voldati))
PRINT*,MAXVAL(vol_iso0_rad%voldati, mask=c_e(vol_iso0_rad%voldati))
PRINT*,MINVAL(vol_iso0_rad%voldati, mask=c_e(vol_iso0_rad%voldati))
!PRINT*,'dati fuori griglia',COUNT(.NOT.c_e(vol_iso0_rad%voldati)),SIZE(vol_iso0_rad%voldati)

CALL getarg(iargc(), output_file)
ALLOCATE(voltmp(1))
voltmp(1) = vol_iso0_rad
CALL export(voltmp, filename=output_file)
! compute POH
! interpolate once more to template
! export to file

END PROGRAM prodsim_compute_poh
