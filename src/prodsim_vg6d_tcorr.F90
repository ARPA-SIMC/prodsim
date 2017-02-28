! export LOG4C_PRIORITY=info
! export LOG4C_APPENDER=stderr
PROGRAM prodsim_vg6d_tcorr
USE log4fortran
USE err_handling
USE missing_values
USE char_utilities
USE phys_const
USE optionparser_class
USE vol7d_var_class
USE grid_class
USE volgrid6d_var_class
USE vol7d_level_class
USE volgrid6d_class
IMPLICIT NONE

INTEGER :: category, ier, i, j, k, gridsize, tindex
CHARACTER(len=512) :: a_name, input_orography, output_orography, &
 input_file, output_file
CHARACTER(len=12) :: tcorr_method
REAL :: tgrad
TYPE(optionparser) :: opt
INTEGER :: optind, optstatus
LOGICAL :: version, ldisplay
TYPE(vol7d_var) :: varbufr
TYPE(volgrid6d),POINTER  :: volgrid(:),  volgrid_tmp(:), volgrid_io(:), volgrid_oo(:)

! innitialise logging
CALL l4f_launcher(a_name,a_name_force='prodsim_vg6d_tcorr')
ier=l4f_init()
! set a_name
category=l4f_category_get(a_name//'.main')

! define the option parser
opt = optionparser_new(description_msg= &
 'Tool for correcting near-surface temperature according to height difference &
 &between model orography and detailed orography.', &
 usage_msg='Usage: prodsim_vg6d_tcorr [options] inputfile outputfile')

CALL optionparser_add(opt, ' ', 'tcorr-method', tcorr_method, default='dry', &
 help='method for determining the vertical temperature gradient, &
 &''dry'' for dry adiabatic gradient, &
 &''user'' for constant gradient &
 &provided by the user with the --tgrad argument')
tgrad = rmiss
CALL optionparser_add(opt, ' ', 'tgrad', tgrad, &
 help='constant vertical temperature gradient in K/m, to be used with &
 &--tcorr-method=user, it should be <0 for temperature decreasing with height')
input_orography = ''
CALL optionparser_add(opt, ' ', 'input-orograhy', input_orography, &
 help='name of file containing the orography associated to the input &
 &temperature data, it should be on the same grid as the input data')
output_orography = ''
CALL optionparser_add(opt, ' ', 'output-orograhy', output_orography, &
 help='name of file containing the target orography to which temperature &
 &should be corrected in output, it should be on the same grid as the input data')

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
  WRITE(*,'(A,1X,A)')'prodsim_vg6d_tcorr','0.1'
  CALL exit(0)
ENDIF

IF (tcorr_method == 'dry') THEN
  tgrad = - gearth/cpd
ELSE IF (tcorr_method == 'user') THEN
  IF (.NOT.c_e(tgrad)) THEN
    CALL l4f_category_log(category,L4F_ERROR, &
     'argument --tcorr-method=user requires specification of --tgrad')
    CALL raise_fatal_error()
  ENDIF
ELSE
  CALL l4f_category_log(category,L4F_ERROR, &
   'value '//TRIM(tcorr_method)//' not valid for argument --tcorr-method')
  CALL raise_fatal_error()
ENDIF

IF (optind+1 /= iargc()) THEN
  CALL optionparser_printhelp(opt)
  CALL l4f_category_log(category,L4F_ERROR,'input and/or output file missing')
  CALL raise_fatal_error()
  CALL EXIT(1)
ENDIF

! last argument is output file
CALL getarg(iargc(), output_file)

CALL l4f_category_log(category,L4F_INFO,'output file: '//TRIM(output_file))

IF (input_orography /= '') THEN
  CALL IMPORT(volgrid_io, filename=input_orography, decode=.TRUE., dup_mode=0, &
   time_definition=0, categoryappend='input_oro')
  IF (SIZE(volgrid_io) > 1) THEN
    CALL l4f_category_log(category, L4F_ERROR, &
     'error '//t2c(SIZE(volgrid_io))//' grids found in '//TRIM(input_orography))
    CALL raise_fatal_error()
  ENDIF
ELSE
  CALL l4f_category_log(category, L4F_ERROR, &
   'error input orography not provided')
  CALL raise_fatal_error()
ENDIF

IF (output_orography /= '') THEN
  CALL IMPORT(volgrid_oo, filename=output_orography, decode=.TRUE., dup_mode=0, &
   time_definition=0, categoryappend='output_oro')
  IF (SIZE(volgrid_oo) > 1) THEN
    CALL l4f_category_log(category, L4F_ERROR, &
     'error '//t2c(SIZE(volgrid_oo))//' grids found in '//TRIM(output_orography))
    CALL raise_fatal_error()
  ENDIF
ELSE
  CALL l4f_category_log(category, L4F_ERROR, &
   'error output orography not provided')
  CALL raise_fatal_error()
ENDIF
! check for orography variables, etc.?

! loop on input file(s)
!DO WHILE(optind <= iargc()-1)
CALL getarg(optind, input_file)
CALL l4f_category_log(category,L4F_INFO,'importing file: '//TRIM(input_file))
CALL IMPORT(volgrid, filename=input_file, decode=.TRUE., dup_mode=0, &
 time_definition=0, categoryappend='input_volume')
IF (.NOT.ASSOCIATED(volgrid)) THEN
  CALL l4f_category_log(category, L4F_ERROR, &
   'error importing input volume from file '//TRIM(input_file))
  CALL raise_fatal_error()
ENDIF

IF (SIZE(volgrid) > 1) THEN
  CALL l4f_category_log(category, L4F_ERROR, &
   'error '//t2c(SIZE(volgrid))//' grids found in '//TRIM(input_file))
  CALL raise_fatal_error()
ENDIF

! round the volume to flatten similar level and timeranges
!CALL rounding(volgrid, volgrid_tmp, level=almost_equal_levels, &
! nostatproc=.TRUE.)
!CALL delete(volgrid)
!volgrid => volgrid_tmp
!NULLIFY(volgrid_tmp)

IF (ldisplay) THEN
  PRINT*,'input volume >>>>>>>>>>>>>>>>>>>>'
  CALL display(volgrid)
ENDIF

! check for consistency
IF (SIZE(volgrid(1)%level) /= 1) THEN
  CALL l4f_category_log(category, L4F_ERROR, &
   'error '//t2c(SIZE(volgrid(1)%level))//' levels found in '//TRIM(input_file))
  CALL raise_fatal_error()
ENDIF

j = volgrid(1)%griddim%dim%nx*volgrid(1)%griddim%dim%ny
gridsize = SIZE(volgrid(1)%voldati,1)*SIZE(volgrid(1)%voldati,2)
IF (j /= gridsize .OR. gridsize == 0) THEN
  CALL l4f_category_log(category, L4F_ERROR, &
   'error inconsistent grid sizes: '//t2c(j)//','//t2c(gridsize)//' in '//TRIM(input_file))
  CALL raise_fatal_error()
ENDIF

!PRINT*,ASSOCIATED(volgrid(1)%griddim%dim%lon),ASSOCIATED(volgrid(1)%griddim%dim%lat)
! if coordinates of input grid are needed, do the following, then
! coordinates will be allocated in arrays volgrid(1)%griddim%dim%lon
! volgrid(1)%griddim%dim%lat
!  CALL unproj(volgrid(1)%griddim)
!  PRINT*,ASSOCIATED(volgrid(1)%griddim%dim%lon),ASSOCIATED(volgrid(1)%griddim%dim%lat)

tindex = imiss
DO j = 1, SIZE(volgrid(1)%var)
  varbufr = convert(volgrid(1)%var(j))
  IF (varbufr%btable == 'B12101') THEN
    tindex = j
    EXIT
  ENDIF
ENDDO
! avoid check on variable for now
tindex = 1

IF (c_e(tindex)) THEN
  DO k = 1, SIZE(volgrid(1)%timerange)
    DO j = 1, SIZE(volgrid(1)%time)
      DO i = 1, SIZE(volgrid(1)%level)
        IF (volgrid(1)%level(i)%level1 == 1 .OR. &
         volgrid(1)%level(i)%level1 == 103) THEN ! only fixed height over surface
! (x,y,level,time,timerange,var)
          WHERE(c_e(volgrid(1)%voldati(:,:,i,j,k,tindex)))
            volgrid(1)%voldati(:,:,i,j,k,tindex) = &
             volgrid(1)%voldati(:,:,i,j,k,tindex) + &
             tgrad*(volgrid_oo(1)%voldati(:,:,1,1,1,1) - &
             volgrid_io(1)%voldati(:,:,1,1,1,1))
          END WHERE
        ENDIF
      ENDDO
    ENDDO
  ENDDO
ENDIF

CALL export(volgrid, output_file)
CALL delete(volgrid)
!  optind = optind + 1
!ENDDO

END PROGRAM prodsim_vg6d_tcorr
