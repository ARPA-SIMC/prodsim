PROGRAM prodsim_compute_anomaly
USE log4fortran
USE grid_class
USE volgrid6d_class
USE vol7d_level_class
USE array_utilities
USE err_handling
USE char_utilities
USE optionparser_class
USE datetime_class
IMPLICIT NONE

INTEGER :: category, ier, i, j, n, gridsize, climmonth, runhour
TYPE(arrayof_integer) :: timerange_list
CHARACTER(len=512) :: a_name, input_file1, input_file2, output_file
TYPE(optionparser) :: opt
INTEGER :: optind, optstatus
LOGICAL :: version, ldisplay

TYPE(volgrid6d),POINTER  :: volgrid(:),  volgrid_tmp(:),  volgridclimt(:)

!questa chiamata prende dal launcher il nome univoco
CALL l4f_launcher(a_name,a_name_force="anomalie_prova")

!init di log4fortran
ier=l4f_init()

!imposta a_name
category=l4f_category_get(TRIM(a_name)//".main")

! define the option parser
opt = optionparser_new(description_msg= &
 'Compute anomalies with respect to a monthly climate dataset.', &
 usage_msg='Usage: prodsim_compute_anomaly [options] climatefile inputfile outputfile')

! display option
CALL optionparser_add(opt, ' ', 'display', ldisplay, help= &
 'briefly display the data volumes imported.')

CALL optionparser_add(opt, ' ', 'timerange-list', timerange_list, help= &
 'comma-separated list of timeranges (end of period) in hour to keep in output')

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
  WRITE(*,'(A,1X,A)')'interpolamare','0.1'
  CALL exit(0)
ENDIF

IF (optind+2 > iargc()) THEN
  CALL optionparser_printhelp(opt)
  CALL l4f_category_log(category,L4F_ERROR,'input and/or output file missing')
  CALL raise_fatal_error()
  CALL EXIT(1)
ENDIF

! last argument is output file
CALL getarg(iargc(), output_file)

CALL l4f_category_log(category,L4F_INFO,'output file: '//TRIM(output_file))


! climate file
CALL getarg(optind, input_file1)
CALL l4f_category_log(category,L4F_INFO,'importing file: '//TRIM(input_file1))
CALL IMPORT(volgridclimt, filename=input_file1, decode=.TRUE., dup_mode=0, &
 time_definition=0)!, categoryappend="input_volume")
IF (.NOT.ASSOCIATED(volgridclimt)) THEN
  CALL l4f_category_log(category, L4F_ERROR, &
   'error importing input volume from file '//TRIM(input_file1))
  CALL raise_fatal_error()
ENDIF

n = volgridclimt(1)%griddim%dim%nx*volgridclimt(1)%griddim%dim%ny

gridsize = SIZE(volgridclimt(1)%voldati,1)*SIZE(volgridclimt(1)%voldati,2)
IF (n /= gridsize .OR. gridsize == 0) THEN
  CALL l4f_category_log(category, L4F_ERROR, &
   'error inconsistent grid sizes: '//t2c(j)//','//t2c(gridsize)//' in '//TRIM(input_file1))
  CALL raise_fatal_error()
ENDIF
CALL unproj(volgridclimt(1)%griddim)  
!call getval(volgridclimt(1)%time(1), month=climmonth, hour=runhour) !inserisce nella variabile climmonth il mese (month) contenuto in questo dato
!CALL display(volgridclimt(1)%time(1))

! actual file
optind = optind + 1
CALL getarg(optind, input_file2)
CALL l4f_category_log(category,L4F_INFO,'importing file: '//TRIM(input_file2))
CALL IMPORT(volgrid, filename=input_file2, decode=.TRUE., dup_mode=0, &
 time_definition=0)
IF (.NOT.ASSOCIATED(volgrid)) THEN
  CALL l4f_category_log(category, L4F_ERROR, &
   'error importing input volume from file '//TRIM(input_file2))
  CALL raise_fatal_error()
ENDIF

CALL rounding(volgrid, volgrid_tmp, level=almost_equal_levels, nostatproc=.TRUE.)
CALL delete(volgrid)
volgrid => volgrid_tmp
NULLIFY(volgrid_tmp)

j = volgrid(1)%griddim%dim%nx*volgrid(1)%griddim%dim%ny
gridsize = SIZE(volgrid(1)%voldati,1)*SIZE(volgrid(1)%voldati,2)
IF (j /= gridsize .OR. gridsize == 0) THEN
  CALL l4f_category_log(category, L4F_ERROR, &
   'error inconsistent grid sizes: '//t2c(j)//','//t2c(gridsize)//' in '//TRIM(input_file2))
  CALL raise_fatal_error()
ENDIF

CALL unproj(volgrid(1)%griddim)

CALL getval(volgrid(1)%time(1), month=climmonth, hour=runhour)
!CALL display(volgrid(1)%timerange(5))
!
!PRINT*,'quinto timerange:',volgrid(1)%timerange(5)%timerange, volgrid(1)%timerange(5)%p1, &
! volgrid(1)%timerange(5)%p2

! filter timeranges if requested (e.g. for locating maximum or minimum)
! better to filter before import on the gridinfo object
CALL packarray(timerange_list)
IF (timerange_list%arraysize > 0) THEN
! convert to seconds
  timerange_list%array(:) = timerange_list%array(:)*3600
  DO i = 1, SIZE(volgrid(1)%timerange)
    IF (ALL(volgrid(1)%timerange(i)%p1 /= timerange_list%array(:))) THEN ! to be eliminated
      volgrid(1)%voldati(:,:,1,1,i,:) = rmiss
      DO j = 1, SIZE(volgrid(1)%var)
        CALL delete(volgrid(1)%gaid(1,1,i,j))
      ENDDO
    ENDIF
  ENDDO
ENDIF

! compute anomalies
DO j = 1, SIZE(volgrid(1)%var)
  CALL display(volgrid(1)%var(j))
  PRINT*,'variabile:',volgrid(1)%var(j)%category,volgrid(1)%var(j)%number
  DO i = 1, SIZE(volgrid(1)%timerange)
    WHERE (c_e(volgrid(1)%voldati(:,:,1,1,i,j)) .AND. c_e(volgridclimt(1)%voldati(:,:,1,climmonth,1,1)))
      volgrid(1)%voldati(:,:,1,1,i,j)=(volgrid(1)%voldati(:,:,1,1,i,j)-volgridclimt(1)%voldati(:,:,1,climmonth,1,1)) ! (x,y,level,time,timerange,var)
    ELSEWHERE
      volgrid(1)%voldati(:,:,1,1,i,j)=rmiss
    END WHERE
  ! PRINT*, SIZE(volgrid(1)%voldati(:,:,1,1,i,j)-volgridclimt(1)%voldati(:,:,1,1,1,j))
  ENDDO
ENDDO
CALL EXPORT(volgrid, filename=output_file)

END PROGRAM prodsim_compute_anomaly
