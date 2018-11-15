MODULE simc_radar_mo
USE log4fortran
USE missing_values
USE err_handling
USE char_utilities
USE array_utilities
USE datetime_class
USE geo_coord_class
USE vol7d_level_class
USE vol7d_timerange_class
USE volgrid6d_var_class
USE grid_class
USE volgrid6d_class
USE netcdf
IMPLICIT NONE

PRIVATE

TYPE radar_desc
  TYPE(geo_coord) :: coord
  DOUBLE PRECISION :: radius
  REAL,ALLOCATABLE :: mask(:,:)
  CHARACTER(len=100) :: name
END TYPE radar_desc

PUBLIC radar_desc, radar_desc_make_mask, read_netcdf, &
 compute_poh_1, compute_poh_2

CONTAINS

SUBROUTINE radar_desc_make_mask(this, grid)
TYPE(radar_desc),INTENT(inout) :: this
TYPE(griddim_def),INTENT(inout) :: grid ! inout necessario per unproj

INTEGER :: i, j
TYPE(geo_coord) :: gp

CALL unproj(grid)
ALLOCATE(this%mask(grid%dim%nx, grid%dim%ny))
this%mask(:,:) = rmiss
DO j = 1, grid%dim%ny
  DO i = 1, grid%dim%nx
    CALL init(gp, grid%dim%lon(i,j) , grid%dim%lat(i,j))
    IF (geo_coord_dist(gp, this%coord) < this%radius) THEN
      this%mask(i,j) = 1.0
    ENDIF
  ENDDO
ENDDO

END SUBROUTINE radar_desc_make_mask


SUBROUTINE read_netcdf(name_nc, vol_nc, timelist, timestep, radarlist, category)
CHARACTER(len=*),INTENT(in) :: name_nc
TYPE(volgrid6d),INTENT(inout) :: vol_nc
TYPE(datetime),INTENT(in) :: timelist(:)
TYPE(timedelta),INTENT(in) :: timestep
TYPE(radar_desc),INTENT(inout) :: radarlist(:)
INTEGER,INTENT(in) :: category

INTEGER :: i, j, il, ts, nradar
INTEGER :: ncid, ndim, nvar, natt, nunlimdimid
INTEGER :: id_lon, id_lat, id_geo, id_mesh, id_time
INTEGER :: dim_lon, dim_lat, dim_geo, dim_mesh
INTEGER :: varid, varid_time, dim_time
INTEGER :: delta_t(1),dimid(7)
TYPE(datetime) :: rtime, vtime
INTEGER :: varid_h, varid_geo, varid_mesh, varid_lon, varid_lat 

CHARACTER(len=100) :: time, varname, chmiss, radname
REAL, ALLOCATABLE :: databuff(:,:,:)
REAL :: hmiss
DOUBLE PRECISION, ALLOCATABLE :: geo_lim(:), mesh_xy(:), lon(:), lat(:)
TYPE(griddim_def) :: griddim_nc

! Apertura del file netcdf 
CALL check(nf90_open(name_nc, NF90_NOWRITE, ncid), "Apertura")
! Estraggo le dimensioni delle variabili immagazzinate
CALL check(nf90_inq_dimid(ncid, "lon", id_lon))
CALL check(nf90_inq_dimid(ncid, "lat", id_lat))
CALL check(nf90_inq_dimid(ncid, "geo_dim", id_geo))
CALL check(nf90_inq_dimid(ncid, "mesh_dim", id_mesh))
CALL check(nf90_inquire_dimension(ncid, id_lon, len=dim_lon))
CALL check(nf90_inquire_dimension(ncid, id_lat, len=dim_lat))
CALL check(nf90_inquire_dimension(ncid, id_geo, len=dim_geo))
CALL check(nf90_inquire_dimension(ncid, id_mesh, len=dim_mesh))
! Estraggo le variabili
!call check(nf90_inq_varid(ncid, "lon", varid_lon))
!call check(nf90_inq_varid(ncid, "lat", varid_lat))
CALL check(nf90_inq_varid(ncid, "geo_dim", varid_geo))
CALL check(nf90_inq_varid(ncid, "mesh_dim", varid_mesh))
!CALL check(nf90_inquire_attribute(ncid, NF90_GLOBAL, "RADARS_NAME", attnum=radname_id))

! get radar name as attribute andidentify it
CALL check(nf90_get_att(ncid, NF90_GLOBAL, "RADARS_NAME", radname))
DO i = 1, SIZE(radarlist)
  IF (TRIM(radarlist(i)%name) == TRIM(radname)) THEN
    nradar = i
    PRINT*,'Found radar',TRIM(radname)
    EXIT
  ENDIF
ENDDO
IF (i > SIZE(radarlist)) THEN
  CALL l4f_category_log(category,L4F_ERROR,'Unknown radar '// &
   TRIM(radname))
  CALL raise_fatal_error()
ENDIF

! Estraggo l'istante di emissione del dato
CALL check(nf90_inq_dimid(ncid, "time", id_time) )
CALL check(nf90_inq_varid(ncid, "time", varid_time) )
CALL check(nf90_get_att(ncid, varid_time, "units", time) )
CALL check(nf90_inquire_dimension(ncid, id_time, len=dim_time) )

IF (time(1:12) == 'hours since ') THEN
  il = 12
  ts = 1
ELSE IF (time(1:11) == 'hour since ') THEN
  il = 11
  ts = 1
ELSE IF (time(1:11) == 'hours before ') THEN
  il = 13
  ts = -1
ELSE IF (time(1:11) == 'hour before ') THEN
  il = 12
  ts = -1
ELSE
  CALL l4f_category_log(category,L4F_ERROR,'time unit: '//TRIM(time)//' non gestito')
  CALL raise_fatal_error()
ENDIF
CALL init(rtime, isodate=time(il+1:))
PRINT*,'Detected time:',to_char(rtime)

CALL check(nf90_inquire(ncid, ndim, nvar, natt, nunlimdimid))
varid_h = imiss
DO varid = 1, nvar
  CALL check(nf90_inquire_variable(ncid, varid, name=varname, ndims=j, dimids=dimid))
  IF (varname == "HGHT") THEN
    IF (.NOT.ALL(dimid(1:3) == (/id_lon,id_lat,id_time/))) THEN
      CALL l4f_category_log(category,L4F_ERROR,'dimensions for variable ' &
       //TRIM(varname)//' not in the expected order')
      CALL raise_fatal_error()
    ENDIF
    varid_h = varid
  ENDIF
ENDDO
IF (.NOT.c_e(varid_h)) THEN
  CALL l4f_category_log(category,L4F_ERROR,'variabile HGHT non trovata')
  CALL raise_fatal_error()
ENDIF
CALL check(nf90_get_att(ncid, varid_h, "var_missing", chmiss))
j = f_nblnk(chmiss)
hmiss = c2r(TRIM(chmiss(j:))) ! perche' non funziona c2r (hmiss = rmiss)?
!PRINT*,'missingc',TRIM(chmiss(j:)),'missingr',hmiss
hmiss = -9999.90 ! per ora faccio a mano

!ALLOCATE (lon(dim_lon), lat(dim_lat), geo_lim(dim_geo), mesh_xy(dim_mesh))
ALLOCATE (geo_lim(dim_geo), mesh_xy(dim_mesh))

!CALL check(nf90_get_var(ncid, varid_lat, lat))
!CALL check(nf90_get_var(ncid, varid_lon, lon))
CALL check(nf90_get_var(ncid, varid_geo, geo_lim))
CALL check(nf90_get_var(ncid, varid_mesh, mesh_xy))

ALLOCATE(databuff(dim_lon,dim_lat,1))

IF (.NOT.ASSOCIATED(vol_nc%voldati)) THEN
! first time allocate volume
  CALL init(griddim_nc, proj_type='regular_ll', nx=dim_lon, ny=dim_lat, &
   xmin=geo_lim(2), xmax=geo_lim(4), ymin=geo_lim(1), ymax=geo_lim(3), &
   component_flag=0, categoryappend="generated")
  CALL init(vol_nc, griddim_nc, 1)

  CALL volgrid6d_alloc(vol_nc, ntime=SIZE(timelist), nlevel=1, ntimerange=1, &
   nvar=1, ini=.TRUE.)

  CALL volgrid6d_alloc_vol(vol_nc, ini=.TRUE., inivol=.TRUE., decode=.TRUE.)
  vol_nc%time(:) = timelist(:)
  vol_nc%timerange(1) = vol7d_timerange_new(254, 0, 0)
  vol_nc%level(1) = vol7d_level_new(1)
  vol_nc%var(1) = volgrid6d_var_new(80, 2, 6)
! create also radar mask
  DO i = 1, SIZE(radarlist)
    CALL radar_desc_make_mask(radarlist(i), vol_nc%griddim)
  ENDDO

ELSE
  ! volume already allocated, just check grid
  IF (dim_lon /= griddim_nc%dim%nx .OR. dim_lat /= griddim_nc%dim%ny) THEN
    CALL l4f_category_log(category,L4F_ERROR,'netcdf file '// &
     TRIM(name_nc)//' has a grid different from that of previous file(s)')
    CALL raise_fatal_error()
  ENDIF
ENDIF

DO i = 1, dim_time
  CALL check(nf90_get_var(ncid, varid_time, delta_t, (/i/),(/i/)))
  vtime = rtime + ts*timedelta_new(hour=delta_t(1)) + timestep/2
  vtime = vtime + MOD(vtime, timestep) ! round to step (tipically 1h)
  j = firsttrue(vtime == vol_nc%time)
  PRINT*,'Assigned time:',to_char(vtime),' found at pos.',j
  IF (j > 0) THEN
! get data array
    CALL check(nf90_get_var(ncid, varid_h, databuff, (/1,1,i/),(/dim_lon,dim_lat,i/)))
! apply radar mask
    WHERE(.NOT.c_e(radarlist(nradar)%mask(:,:)))
      databuff(:,:,1) = rmiss
    END WHERE
! merge with existing data
    WHERE(c_e(databuff(:,:,1)))
      WHERE(c_e(vol_nc%voldati(:,:,1,j,1,1)))
        vol_nc%voldati(:,:,1,j,1,1) = MAX(vol_nc%voldati(:,:,1,j,1,1), databuff(:,:,1))
      ELSEWHERE
        vol_nc%voldati(:,:,1,j,1,1) = databuff(:,:,1)
      ENDWHERE
    ENDWHERE
  ELSE
    CALL l4f_category_log(category,L4F_WARN,'time in netcdf file does not fit in available interval')
  ENDIF
ENDDO

PRINT*,'n. dati mancanti in volgrid',COUNT(.NOT.c_e(vol_nc%voldati))
PRINT*,'n. dati validi in volgrid',COUNT(c_e(vol_nc%voldati))
PRINT*,'n. dati validi e con 45dbZ in volgrid',COUNT(c_e(vol_nc%voldati) .AND. vol_nc%voldati /= hmiss)
PRINT*,'n. dati validi e senza 45dbZ in volgrid',COUNT(c_e(vol_nc%voldati) .AND. vol_nc%voldati == hmiss)

CALL check(nf90_close(ncid))

CONTAINS

SUBROUTINE check(status, message)
INTEGER,INTENT(in) :: status
CHARACTER(len=*),OPTIONAL :: message    

IF (status /= nf90_noerr) THEN
  CALL l4f_category_log(category,L4F_ERROR,&
   nf90_strerror(status)//optio_c(message,40))
  CALL raise_fatal_error()
ENDIF
END SUBROUTINE check

END SUBROUTINE read_netcdf

ELEMENTAL FUNCTION compute_poh_1(ecotop, iso0) RESULT(poh)
REAL,INTENT(in) :: ecotop
REAL,INTENT(in) :: iso0

REAL :: poh
REAL :: delta

IF (c_e(ecotop) .AND. c_e(iso0)) THEN
  IF (iso0 > -100. .AND. ecotop > -1.) THEN
    delta = ecotop - iso0/1000.
    poh = 0.319 + 0.133*delta
    poh = MAX(0.,MIN(1.,poh))
  ELSE
    poh = 0.
  ENDIF
ELSE
  poh = rmiss
ENDIF
    
END FUNCTION compute_poh_1

ELEMENTAL FUNCTION compute_poh_2(ecotop, iso0) RESULT(poh)
REAL,INTENT(in) :: ecotop
REAL,INTENT(in) :: iso0

REAL :: poh
REAL :: delta

IF (c_e(ecotop) .AND. c_e(iso0)) THEN
  IF (iso0 > -100. .AND. ecotop > -1.) THEN
    delta = ecotop - iso0/1000.
    poh = -1.20231 + 1.00184*Delta - 0.17018*delta**2 + 0.01086*delta**3
    poh = MAX(0.,MIN(1.,poh))
  ELSE
    poh = 0.
  ENDIF
ELSE
  poh = rmiss
ENDIF
    
END FUNCTION compute_poh_2

END MODULE simc_radar_mo
