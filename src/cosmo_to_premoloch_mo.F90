MODULE cosmo_to_premoloch_mo
USE grib_api
IMPLICIT NONE

INTEGER,PARAMETER :: wp=SELECTED_REAL_KIND(6)
! double  precision SELECTED_REAL_KIND(13)
! single  precision SELECTED_REAL_KIND(6)

INTEGER,PARAMETER :: nsoiltyp=10

REAL  (KIND=wp)     ::  &
!   a) parameters describing the soil water budget
 cporv (nsoiltyp), &  !  pore volume (fraction of volume)
 cfcap (nsoiltyp), &  !  field capacity (fraction of volume)
 cpwp  (nsoiltyp), &  !  plant wilting point (fraction of volume)
 cadp  (nsoiltyp)     !  air dryness point (fraction of volume)

! soil type:   ice    rock    sand    sandy   loam   clay      clay    peat    sea     sea  
! (by index)                          loam           loam                     water    ice

DATA  cporv / 1.E-10_wp, 1.E-10_wp, 0.364_wp  , 0.445_wp  , 0.455_wp  , 0.475_wp  , 0.507_wp  , 0.863_wp  , 1.E-10_wp, 1.E-10_wp /
DATA  cfcap / 1.E-10_wp, 1.E-10_wp, 0.196_wp  , 0.260_wp  , 0.340_wp  , 0.370_wp  , 0.463_wp  , 0.763_wp  , 1.E-10_wp, 1.E-10_wp /
DATA  cpwp  / 0.0_wp   , 0.0_wp   , 0.042_wp  , 0.100_wp  , 0.110_wp  , 0.185_wp  , 0.257_wp  , 0.265_wp  , 0.0_wp   ,  0.0_wp   /
DATA  cadp  / 0.0_wp   , 0.0_wp   , 0.012_wp  , 0.030_wp  , 0.035_wp  , 0.060_wp  , 0.065_wp  , 0.098_wp  , 0.0_wp   ,  0.0_wp   /

INTEGER :: ni=-1, nj=-1

REAL(kind=wp),ALLOCATABLE :: soiltype(:), landfrac(:)


CONTAINS


SUBROUTINE write_msg(gid, fid, hyblayer, hyblevel)
INTEGER,INTENT(inout) :: gid, fid
LOGICAL,INTENT(in),OPTIONAL :: hyblayer, hyblevel

IF (PRESENT(hyblayer)) THEN
  IF (hyblayer) THEN
    CALL grib_set(gid, 'typeOfFirstFixedSurface', 105)
    CALL grib_set(gid, 'typeOfSecondFixedSurface', 105)
  ENDIF
ENDIF
IF (PRESENT(hyblevel)) THEN
  IF (hyblevel) THEN
    CALL grib_set(gid, 'typeOfFirstFixedSurface', 105)
  ENDIF
ENDIF
CALL grib_set(gid, 'subCentre', 103)
CALL grib_write(gid, fid)

END SUBROUTINE write_msg


FUNCTION check_grid(gid)
INTEGER,INTENT(inout) :: gid
LOGICAL :: check_grid

INTEGER :: lni, lnj, ier

check_grid = .TRUE.
IF (ni == -1 .OR. nj == -1) THEN
  CALL grib_get(gid, 'Ni', ni, ier)
  IF (ier /= GRIB_SUCCESS) check_grid = .FALSE.
  CALL grib_get(gid, 'Nj', nj, ier)
  IF (ier /= GRIB_SUCCESS) check_grid = .FALSE.
ELSE
  CALL grib_get(gid, 'Ni', lni, ier)
  IF (ier /= GRIB_SUCCESS) check_grid = .FALSE.
  CALL grib_get(gid, 'Nj', lnj, ier)
  IF (ier /= GRIB_SUCCESS) check_grid = .FALSE.
  check_grid = check_grid .AND. (lni == ni .AND. lnj == nj)
ENDIF

END FUNCTION check_grid


! allocate array field and store into it the values from the grib
! provided
SUBROUTINE store_field(gid, field)
INTEGER,INTENT(inout) :: gid
REAL(kind=wp),ALLOCATABLE,INTENT(inout) :: field(:)

ALLOCATE(field(ni*nj))
CALL grib_get(gid, 'values', field)

END SUBROUTINE store_field

! useful grib_ls
!grib_ls -p count,discipline,parameterCategory,parameterNumber,shortName,typeOfLevel
!grib_ls -p count,discipline,parameterCategory,parameterNumber,shortName,typeOfFirstFixedSurface,scaleFactorOfFirstFixedSurface,scaledValueOfFirstFixedSurface,typeOfSecondFixedSurface,scaleFactorOfSecondFixedSurface,scaledValueOfSecondFixedSurface

! convert soil moisture
! From column-integrated soil moisture (kg m-2)
! To soil moisture relative index
SUBROUTINE compute_rel_moist(gid)
INTEGER,INTENT(inout) :: gid

REAL(kind=wp),ALLOCATABLE :: rel_moist(:)
INTEGER :: f1, v1, f2, v2, ij, st
REAL(kind=wp) :: tl, bl, levthick

CALL grib_get(gid, 'scaleFactorOfFirstFixedSurface', f1)
CALL grib_get(gid, 'scaledValueOfFirstFixedSurface', v1)
tl = REAL(v1,kind=wp)*10.0_wp**(-f1)
CALL grib_get(gid, 'scaleFactorOfSecondFixedSurface', f2)
CALL grib_get(gid, 'scaledValueOfSecondFixedSurface', v2)
bl = REAL(v2,kind=wp)*10.0_wp**(-f2)
levthick = (bl-tl)/1000.0_wp ! kg m-2 => kg m-3 => m3 m-3

CALL store_field(gid, rel_moist)

DO ij = 1, ni*nj
  st = NINT(soiltype(ij))
  IF (st > 0 .AND. st <= nsoiltyp) THEN
    rel_moist(ij) = (rel_moist(ij)/levthick-cadp(st))/(cporv(st)-cadp(st))
    rel_moist(ij) = MIN(MAX(rel_moist(ij), 0.0_wp), 1.0_wp)
  ELSE
    rel_moist(ij) = 1.0 ! sea?
  ENDIF
ENDDO

! from layer to level
! because of grib_api behavior, second surface has to be set to
! missing before setting the values of the first surface
CALL grib_set(gid, 'typeOfSecondFixedSurface', 255)
CALL grib_set_missing(gid, 'scaleFactorOfSecondFixedSurface')
CALL grib_set_missing(gid, 'scaledValueOfSecondFixedSurface')
CALL grib_set(gid, 'typeOfFirstFixedSurface', 106)
IF (f1 == f2) THEN ! simple case
  IF (MOD(v1+v2, 2) == 0) THEN
    CALL grib_set(gid, 'scaleFactorOfFirstFixedSurface', f1)
    CALL grib_set(gid, 'scaledValueOfFirstFixedSurface', (v1+v2)/2)
  ELSE
    CALL grib_set(gid, 'scaleFactorOfFirstFixedSurface', f1+1)
    CALL grib_set(gid, 'scaledValueOfFirstFixedSurface', (v1+v2)*5)
  ENDIF
ELSE
  CALL grib_set(gid, 'scaleFactorOfFirstFixedSurface', MAX(f1,f2)+1)
  CALL grib_set(gid, 'scaledValueOfFirstFixedSurface', &
   NINT((bl+tl)/2.0_wp*10.0_wp**(MAX(f1,f2)+1)))
ENDIF

CALL grib_set(gid, 'values', rel_moist)

END SUBROUTINE compute_rel_moist


SUBROUTINE compute_ice_cover(gid)
INTEGER,INTENT(inout) :: gid

REAL(kind=wp),ALLOCATABLE :: ice_cover(:)

CALL store_field(gid, ice_cover)
WHERE (ice_cover > 0.0_wp)
  ice_cover = 1.0_wp
ELSEWHERE
  ice_cover = 0.0_wp
END WHERE
CALL grib_set(gid, 'values', ice_cover)

END SUBROUTINE compute_ice_cover


SUBROUTINE process_const_fields(fid, ofid)
INTEGER,INTENT(inout) :: fid, ofid

INTEGER :: gid, ier, d, c, n, l1, l2

DO WHILE(.TRUE.)
  CALL grib_new_from_file(fid, gid, ier)
  IF (ier /= GRIB_SUCCESS) EXIT
  IF (.NOT. check_grid(gid)) THEN
    PRINT*,'Griglia non coincide'
    STOP
  ENDIF

  CALL grib_get(gid, 'discipline', d)
  CALL grib_get(gid, 'parameterCategory', c)
  CALL grib_get(gid, 'parameterNumber', n)

  CALL grib_get(gid, 'typeOfFirstFixedSurface', l1)
  CALL grib_get(gid, 'typeOfSecondFixedSurface', l2)

  IF (l1 == 1) THEN ! surface
    SELECT CASE(d)
    CASE(0) ! atmosphere
      SELECT CASE(c)
      CASE(3)
        SELECT CASE(n)
        CASE(6) ! hsurf
          CALL grib_set(gid, 'discipline', 2)
          CALL grib_set(gid, 'parameterCategory', 0)
          CALL grib_set(gid, 'parameterNumber', 7)
          CALL grib_set(gid, 'typeOfSecondFixedSurface', 255)
          CALL grib_set_missing(gid, 'scaleFactorOfSecondFixedSurface')
          CALL grib_set_missing(gid, 'scaledValueOfSecondFixedSurface')
          CALL write_msg(gid, ofid)
!        CASE(4) ! geopotential
        END SELECT
      END SELECT
    CASE(2) ! surface
      SELECT CASE(c)
      CASE(0)
        SELECT CASE(n)
        CASE(0) ! land fraction
          CALL store_field(gid, landfrac)
!        CALL grib_set(gid, 'parameterNumber', 218)
          CALL write_msg(gid, ofid)
        END SELECT
      CASE(3)
        SELECT CASE(n)
        CASE(196) ! soil type
          CALL store_field(gid, soiltype)
        END SELECT
      END SELECT
    END SELECT
  ENDIF

  CALL grib_release(gid)
ENDDO

END SUBROUTINE process_const_fields


SUBROUTINE process_progn_fields(fid, ofid)
INTEGER,INTENT(inout) :: fid, ofid

INTEGER :: gid, clgid, ier, d, c, n, l1, l2, v1

DO WHILE(.TRUE.)
  CALL grib_new_from_file(fid, gid, ier)
  IF (ier /= GRIB_SUCCESS) EXIT
  IF (.NOT. check_grid(gid)) THEN
    PRINT*,'Griglia non coincide'
    STOP
  ENDIF

  CALL grib_get(gid, 'discipline', d)
  CALL grib_get(gid, 'parameterCategory', c)
  CALL grib_get(gid, 'parameterNumber', n)

  CALL grib_get(gid, 'typeOfFirstFixedSurface', l1)
  CALL grib_get(gid, 'typeOfSecondFixedSurface', l2)

  IF (l1 == 150 .AND. l2 == 150) THEN ! upper air full levels
    SELECT CASE(d)
    CASE(0) ! atmosphere
      SELECT CASE(c)
      CASE(0) ! temperature
        SELECT CASE(n)
        CASE(0) ! temperature
          CALL write_msg(gid, ofid, hyblayer=.TRUE.)
        END SELECT
      CASE(1) ! moisture
        SELECT CASE(n)
        CASE(0) ! q
          CALL write_msg(gid, ofid, hyblayer=.TRUE.)
        CASE(22) ! qc?
          CALL grib_set(gid, 'parameterNumber', 83)
          CALL write_msg(gid, ofid, hyblayer=.TRUE.)
        CASE(82) ! qi?
          CALL grib_set(gid, 'parameterNumber', 84)
          CALL write_msg(gid, ofid, hyblayer=.TRUE.)
        END SELECT
      CASE(2) ! momentum
        SELECT CASE(n)
        CASE(2,3) ! u, v
          CALL write_msg(gid, ofid, hyblayer=.TRUE.)
        END SELECT
      CASE(3) ! mass
        SELECT CASE(n)
        CASE(0) ! p
          CALL write_msg(gid, ofid, hyblayer=.TRUE.)
        END SELECT
      END SELECT
    END SELECT

  ELSE IF (l1 == 150 .AND. l2 == 255) THEN ! upper air half levels
    SELECT CASE(d)
    CASE(0) ! atmosphere
      SELECT CASE(c)
      CASE(2) ! momentum
        SELECT CASE(n)
        CASE(9) ! w
          CALL write_msg(gid, ofid, hyblevel=.TRUE.)
        END SELECT
      END SELECT
    END SELECT

  ELSE IF (l1 == 106) THEN ! soil levels
    SELECT CASE(d)
    CASE(2) ! surface
      SELECT CASE(c)
      CASE(3)
        SELECT CASE(n)
        CASE(18) ! soil temperature

          CALL grib_get(gid, 'scaledValueOfFirstFixedSurface', v1)
          IF (v1 == 0) THEN ! assume l2 == missing, surface
! duplicate to sea temperature
            CALL grib_clone(gid, clgid)
!            CALL grib_set(clgid, 'discipline', 10)
!            CALL grib_set(clgid, 'parameterCategory', 3)
!            CALL grib_set(clgid, 'parameterNumber', 0)
! pretend to be first deep soil level (warning depth may change!)
!            CALL grib_set(clgid, 'scaleFactorOfFirstFixedSurface', 3)
!            CALL grib_set(clgid, 'scaledValueOfFirstFixedSurface', 5)
!            CALL write_msg(clgid, ofid)
! duplicate to skin temperature
            CALL grib_set(clgid, 'discipline', 0)
            CALL grib_set(clgid, 'parameterCategory', 0)
            CALL grib_set(clgid, 'parameterNumber', 0)
            CALL grib_set(clgid, 'typeOfFirstFixedSurface', 1)
            CALL grib_set_missing(clgid, 'scaleFactorOfFirstFixedSurface')
            CALL grib_set_missing(clgid, 'scaledValueOfFirstFixedSurface')
            CALL write_msg(clgid, ofid)
            CALL grib_release(clgid)
          ELSE ! write only deeper levels
            CALL grib_set(gid, 'discipline', 2)
            CALL grib_set(gid, 'parameterCategory', 0)
            CALL grib_set(gid, 'parameterNumber', 2)
            CALL write_msg(gid, ofid)
          ENDIF
        CASE(20) ! soil moisture
          CALL grib_set(gid, 'parameterCategory', 0)
          CALL grib_set(gid, 'parameterNumber', 198) ! local
          CALL  compute_rel_moist(gid)
!        CALL grib_set(gid, 'discipline', 2)
          CALL write_msg(gid, ofid)
        CASE(22) ! soil ice
          CALL  compute_rel_moist(gid)
!        CALL grib_set(gid, 'discipline', 2)
          CALL grib_set(gid, 'parameterCategory', 0)
          CALL grib_set(gid, 'parameterNumber', 199) ! local
          CALL write_msg(gid, ofid)
        END SELECT
      END SELECT
    END SELECT

  ELSE IF (l1 == 1) THEN ! surface
    SELECT CASE(d)
    CASE(0) ! atmosphere
      SELECT CASE(c)
      CASE(1) ! moisture
        SELECT CASE(n)
        CASE(60) ! snow depth
          CALL grib_set(gid, 'parameterNumber', 13) ! ~
          CALL write_msg(gid, ofid)
        END SELECT
      CASE(3) ! mass
        SELECT CASE(n)
        CASE(0) ! pressure
          CALL write_msg(gid, ofid)
        END SELECT
      END SELECT
    CASE(10) ! ocean
      SELECT CASE(c)
      CASE(2) ! ice
        SELECT CASE(n)
        CASE(1) ! ice thickness
          CALL write_msg(gid, ofid)
          CALL compute_ice_cover(gid)
          CALL grib_set(gid, 'parameterNumber', 0)
          CALL write_msg(gid, ofid)
        END SELECT
      END SELECT
    END SELECT
  ENDIF
  CALL grib_release(gid)
ENDDO

END SUBROUTINE process_progn_fields

END MODULE cosmo_to_premoloch_mo
