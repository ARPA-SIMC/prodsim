PROGRAM prodsim_vg6d_gen_obs_mask
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
CHARACTER(len=12) :: coord_format, radius_unit
CHARACTER(len=512) :: a_name, coord_file, grid_file, box_file, output_file
REAL :: radius
TYPE(optionparser) :: opt
INTEGER :: optind, optstatus
LOGICAL :: version, ldisplay
TYPE(vol7d) :: v7d_coord
TYPE(vol7d_dballe) :: v7d_ana
TYPE(arrayof_georef_coord_array) :: poly
INTEGER :: polytopo
TYPE(volgrid6d),POINTER  :: volgrid(:),  volgrid_tmp(:), volgrid_io(:), volgrid_oo(:)

TYPE(arrayof_doubleprecision) :: lon, lat
DOUBLE PRECISION,ALLOCATABLE :: lon_array(:), lat_array(:)

! innitialise logging
CALL l4f_launcher(a_name,a_name_force='prodsim_vg6d_gen_obs_mask')
ier=l4f_init()
! set a_name
category=l4f_category_get(a_name//'.main')

! define the option parser
opt = optionparser_new(description_msg= &
 'Tool for generating a mask of grid points having a &
 &minimum distance from a set of sparse points (stations) &
 &higher than a specified threshold. &
 &For each grid point, the check can be performed on set of sparse points &
 &lying in a (coarser) grid box in which also the grid point is contained &
 &or on all the set of sparse points provided.', &
 usage_msg='Usage: prodsim_vg6d_gen_obs_mask [options] gridfile outputfile')

CALL optionparser_add(opt, 'a', 'lon', lon, (/10.D0/), help= &
 'comma-separated list of longitudes of sparse points, &
 &alternative to --coord-file')
CALL optionparser_add(opt, 'b', 'lat', lat, (/45.D0/), help= &
 'comma-separated list of latitudes of sparse points, &
 &alternative to --coord-file')
CALL optionparser_add(opt, 'c', 'coord-file', coord_file, help= &
 'file with coordinates of sparse points, alternative to --lon, --lat')
coord_file=cmiss
CALL optionparser_add(opt, ' ', 'coord-format', coord_format, 'BUFR', &
 & help='format of input file with coordinates, ''native'' for vol7d native binary file &
 & ''BUFR'' for BUFR file, ''CREX'' for CREX file ''shp'' for shapefile')

CALL optionparser_add(opt, ' ', 'radius', radius, 20000., help= &
 'radius of search for sparse points around each grid point in m')
!CALL optionparser_add(opt, ' ', 'radius-unit', radius_unit, 'm', help= &
! 'unit of radius specified by --radius, acceptable values ''m'' (meters) &
! &or ''grid'' (fractionary units of grid boxes)')

!CALL optionparser_add(opt, ' ', 'grid-file', grid_file, '', help= &
! 'grib file defining the grid, output file will be on this grid')
CALL optionparser_add(opt, ' ', 'box-file', box_file, '', help= &
 'grib file defining the boxes over which computations are performed, if empty &
 & the computations will be performed on the whole grid domain, &
 &if provided it should have a coarser resolution than grid-file, &
 &but it may be on a different geographical projection')

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
  WRITE(*,'(A,1X,A)')'prodsim_vg6d_gen_obs_mask','0.1'
  CALL exit(0)
ENDIF

IF (optind+1 /= iargc()) THEN
  CALL optionparser_printhelp(opt)
  CALL l4f_category_log(category,L4F_ERROR,'grid and/or output file missing')
  CALL raise_fatal_error()
  CALL EXIT(1)
ENDIF

! last argument is output file
CALL getarg(optind, grid_file)
CALL getarg(optind+1, output_file)


IF (c_e(coord_file)) THEN
  CALL l4f_category_log(category,L4F_DEBUG,'start import coord')
  IF (coord_format == 'native') THEN
    CALL l4f_category_log(category,L4F_DEBUG,'execute import coord native')
    CALL import(v7d_coord, filename=coord_file)

  ELSE IF (coord_format == 'BUFR' .OR. coord_format == 'CREX') THEN
    CALL l4f_category_log(category,L4F_DEBUG,'execute import coord v7d')
    CALL init(v7d_ana, filename=coord_file, format=coord_format, file=.TRUE., &
     write=.FALSE., categoryappend="anagrafica")
    CALL import(v7d_ana, anaonly=.TRUE.)
    v7d_coord = v7d_ana%vol7d
! destroy v7d_ana without deallocating the contents passed to v7d
    CALL init(v7d_ana%vol7d)
    CALL delete(v7d_ana)

  ELSE IF (coord_format == 'shp') THEN
    CALL l4f_category_log(category,L4F_DEBUG,'execute import coord shape')
    CALL import(poly, coord_file)
    IF (poly%arraysize <= 0) THEN
      CALL l4f_category_log(category, L4F_ERROR, &
       'error importing shapefile '//TRIM(coord_file))
      CALL raise_fatal_error()
    ENDIF
    CALL getval(poly%array(1), topo=polytopo)
    IF (polytopo == georef_coord_array_point) THEN ! topology suitable for sparse points
      CALL vol7d_alloc(v7d_coord, nana=poly%arraysize)
      CALL vol7d_alloc_vol(v7d_coord)
      DO i = 1, poly%arraysize
        CALL getval(poly%array(i), x=lon_array,y=lat_array) ! improve!!!!
        CALL init(v7d_coord%ana(i), lon=lon_array(1), lat=lat_array(1))
      ENDDO
      CALL delete(poly)
      CALL l4f_category_log(category,L4F_INFO, &
       'shapefile '//TRIM(coord_file)//' interpreted as sparse point list')
    ELSE ! topology suitable for polygons, nothing to do
      CALL l4f_category_log(category,L4F_ERROR,'shapefile '//TRIM(coord_file)//&
       ' contains polygons instead of points')
      CALL raise_fatal_error()
      CALL EXIT(1)
    ENDIF

  ELSE IF (coord_format == 'grib_api') THEN
    CALL l4f_category_log(category,L4F_DEBUG,'execute import coord grib_api')
    CALL import(maskgrid, coord_file, categoryappend='maskgrid')
    IF (maskgrid%arraysize < 1) THEN
      CALL l4f_category_log(category, L4F_ERROR, &
       'error importing mask grid file '//TRIM(coord_file))
      CALL raise_fatal_error()
    ENDIF
    CALL import(maskgrid%array(1))
    ALLOCATE(maskfield(maskgrid%array(1)%griddim%dim%nx, maskgrid%array(1)%griddim%dim%ny))
    maskfield(:,:) = decode_gridinfo(maskgrid%array(1))
    CALL delete(maskgrid)

  ELSE
    CALL l4f_category_log(category, L4F_ERROR, &
     'error in command-line arguments, format '// &
     TRIM(coord_format)//' in --coord-format not valid or not supported.')
    CALL raise_fatal_error()
  ENDIF

  CALL l4f_category_log(category,L4F_DEBUG,'end import coord')
ENDIF




! from grid_transform_class, code for mapping points to box:
! for grid
! in e out TYPE(griddim_def)
CALL get_val(in, nx=innx, ny=inny)
CALL griddim_setsteps(out)
CALL get_val(out, nx=outnx, ny=outny, &
 xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)
! unlike before, here index arrays must have the shape of input grid
ALLOCATE(inter_index_x(innx,inny), &
 inter_index_y(innx,inny))

! compute coordinates of input grid in geo system
CALL unproj(in) ! TODO costringe a dichiarare in INTENT(inout), si puo` evitare?
! use find_index in the opposite way, here extrap does not make sense
CALL find_index(out, 'near', &
 outnx, outny, xmin, xmax, ymin, ymax, &
 in%dim%lon, in%dim%lat, .FALSE., &
 inter_index_x, inter_index_y)

! for sparse points
innx=SIZE(v7d_in%ana)
inny=1
!?!?CALL griddim_setsteps(out)
CALL get_val(out, nx=outnx, ny=outny, &
 xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)
! index arrays must have the shape of input grid
ALLOCATE(lon(innx),lat(inny))
ALLOCATE(inter_index_x(innx,inny), &
 inter_index_y(innx,inny))

! get coordinates of input grid in geo system
CALL getval(v7d_in%ana(:)%coord,lon=lon,lat=lat)
! use find_index in the opposite way, here extrap does not make sense
CALL find_index(out,'near',&
 outnx, outny, xmin, xmax, ymin, ymax, &
 lon, lat, .FALSE., &
 inter_index_x(:,1), inter_index_y(:,1))

END PROGRAM prodsim_vg6d_gen_obs_mask
