PROGRAM prodsim_expand_field
USE log4fortran
USE optionparser_class
USE err_handling
USE missing_values
USE array_utilities
USE gridinfo_class
USE grid_id_class
IMPLICIT NONE

INTEGER :: category, ier
CHARACTER(len=512) :: a_name, input_file, output_file, mask_file
TYPE(optionparser) :: opt
INTEGER :: optind, optstatus
LOGICAL :: version
TYPE(arrayof_real) :: maskbounds
INTEGER :: num_iter, i, fieldshape(2)
TYPE(arrayof_gridinfo) :: maskgrid
REAL,ALLOCATABLE :: maskfield(:,:), field(:,:)
TYPE(gridinfo_def) :: gridinfo
TYPE(grid_file_id) :: ifile, ofile
TYPE(grid_id) :: input_grid_id


CALL l4f_launcher(a_name, a_name_force="expand_field")
ier = l4f_init()
category=l4f_category_get(TRIM(a_name)//".main")

! define the option parser
opt = optionparser_new(description_msg= &
 'Iteratively expand the values of all the fields within a grib file towards &
 &the points having missing value using the average of the neighbouring valid points. &
 &The points with missing values are by default deduced from each input field, &
 &but they can be optionally specified by means of a separate mask field &
 &on the same grid and of an optional interval of values to be considered &
 &as valid in the mask field', &
 usage_msg='Usage: prodsim_expand_field [options] inputfile outputfile')

! define command-line options
mask_file = cmiss
CALL optionparser_add(opt, ' ', 'mask-file', mask_file, help= &
 'a file in grib format specifying an optional mask of missing values, &
 &only the first message is considered, if this argument is not provided &
 &the missing values are taken from the input file itself, &
 &see also --mask-values')
CALL optionparser_add(opt, ' ', 'maskbounds', maskbounds, help= &
 'two comma-separated values indicating the interval of values in &
 &--mask-file that indicate valid points, including the extremes, &
 &if not present, all and only valid values in mask-file &
 &indicate valid points')
CALL optionparser_add(opt, ' ', 'num-iter', num_iter, 1, help= &
 'number of iterations on the expansion process')

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
  WRITE(*,'(A,1X,A)')'prodsim_expand_field',VERSION
  CALL exit(0)
ENDIF

IF (optind + 1 <= iargc()) THEN
  CALL getarg(optind, input_file)
  IF (input_file == '-') THEN
    CALL l4f_category_log(category, L4F_INFO, 'trying /dev/stdin as stdin unit')
    input_file = '/dev/stdin'
  ENDIF

  optind = optind+1
  CALL getarg(optind, output_file)
  IF (output_file == '-') THEN
    CALL l4f_category_log(category, L4F_INFO, 'trying /dev/stdout as stdout unit')
    output_file = '/dev/stdout'
  ENDIF

ELSE
  CALL optionparser_printhelp(opt)
  CALL l4f_category_log(category, L4F_FATAL, 'input or output file missing')
  CALL raise_fatal_error()
ENDIF


IF (c_e(mask_file)) THEN
! import all the file into an array of gridinfo
  CALL IMPORT(maskgrid, mask_file, categoryappend='maskgrid')
  IF (maskgrid%arraysize < 1) THEN
    CALL l4f_category_log(category, L4F_ERROR, &
     'error importing mask grid file '//TRIM(mask_file))
    CALL raise_fatal_error()
  ENDIF
! get first field from gridinfo array
  CALL IMPORT(maskgrid%array(1))
  ALLOCATE(maskfield(maskgrid%array(1)%griddim%dim%nx, maskgrid%array(1)%griddim%dim%ny))
  maskfield(:,:) = decode_gridinfo(maskgrid%array(1))
  CALL delete(maskgrid)

  IF (maskbounds%arraysize >= 2) THEN ! maskbounds provided, set mask accordingly
    WHERE(maskfield < maskbounds%array(1) .OR. maskfield > maskbounds%array(2))
      maskfield = rmiss
    END WHERE
  ENDIF
  fieldshape = SHAPE(maskfield)
ENDIF

ifile = grid_file_id_new(input_file,'r')
ofile = grid_file_id_new(output_file,'w')

DO WHILE (.TRUE.)
  input_grid_id = grid_id_new(ifile)
  IF (.NOT.c_e(input_grid_id)) THEN ! THEN because of a bug in gfortran?!
    EXIT
  ENDIF

  CALL l4f_category_log(category,L4F_INFO,"importing gridinfo")
  CALL init(gridinfo, gaid=input_grid_id, categoryappend="imported")
  CALL IMPORT(gridinfo)

  CALL l4f_category_log(category,L4F_INFO,"gridinfo imported")

  IF (ALLOCATED(maskfield)) THEN ! preliminary check
    IF (fieldshape(1) /= gridinfo%griddim%dim%nx .OR. &
     fieldshape(2) /= gridinfo%griddim%dim%ny) THEN
      CALL l4f_category_log(category, L4F_ERROR, &
       'mask and input grid do not match: '// &
       t2c(fieldshape(1))//' '//t2c(gridinfo%griddim%dim%nx)//', '// &
       t2c(fieldshape(2))//' '//t2c(gridinfo%griddim%dim%nx))
      CALL raise_error()
      CALL delete(gridinfo)
      CYCLE
    ENDIF
  ENDIF
  
  ALLOCATE(field(gridinfo%griddim%dim%nx,gridinfo%griddim%dim%ny))
  field(:,:) = decode_gridinfo(gridinfo)

  IF (ALLOCATED(maskfield)) THEN ! apply mask if required
    WHERE(.NOT.c_e(maskfield))
      field = rmiss
    END WHERE
  ENDIF
  DO i = 1, num_iter
    CALL expand_field(field)
  ENDDO

  CALL encode_gridinfo(gridinfo, field(:,:))
  CALL export(gridinfo)

  CALL export(gridinfo%gaid, ofile)

  CALL delete(gridinfo)
  DEALLOCATE(field)

ENDDO

CALL delete(ifile)
CALL delete(ofile)

CALL l4f_category_log(category,L4F_INFO,"end")

CALL l4f_category_delete(category)
ier=l4f_fini()


CONTAINS


SUBROUTINE expand_field(field)
REAL,INTENT(inout) :: field(:,:)

INTEGER :: nx, ny, i, j, im1, ip1, jm1, jp1
REAL,ALLOCATABLE :: cfield(:,:)

! make a copy to avoid using updated data
cfield = field

nx = SIZE(field, 1)
ny = SIZE(field, 2)

DO j = 1, ny
  jm1 = MAX(j-1, 1)
  jp1 = MIN(j+1, ny)

  DO i = 2, nx - 1
    IF (.NOT.c_e(field(i,j))) THEN
      im1 = MAX(i-1, 1)
      ip1 = MIN(i+1, nx)
      field(i,j) = avg_of_4(cfield(im1,j), cfield(ip1,j), &
       cfield(i,jm1), cfield(i,jp1))
    ENDIF
  ENDDO
ENDDO

END SUBROUTINE expand_field


FUNCTION avg_of_4(v1, v2, v3, v4)
REAL,INTENT(in) :: v1, v2, v3, v4
REAL :: avg_of_4

INTEGER :: nv

avg_of_4 = 0.
nv = 0
! shameful way
IF (c_e(v1)) THEN
  avg_of_4 = avg_of_4 + v1
  nv = nv + 1
ENDIF
IF (c_e(v2)) THEN
  avg_of_4 = avg_of_4 + v2
  nv = nv + 1
ENDIF
IF (c_e(v3)) THEN
  avg_of_4 = avg_of_4 + v3
  nv = nv + 1
ENDIF
IF (c_e(v4)) THEN
  avg_of_4 = avg_of_4 + v4
  nv = nv + 1
ENDIF
IF (nv == 0) THEN
  avg_of_4 = rmiss
ELSE
  avg_of_4 = avg_of_4/nv
ENDIF

END FUNCTION avg_of_4

END PROGRAM prodsim_expand_field
