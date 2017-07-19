MODULE cnmc_lightning_mo
USE kinds
USE missing_values
USE datetime_class
USE georef_coord_class
IMPLICIT NONE

TYPE cnmc_lightning_t
  TYPE(datetime) :: time=datetime_miss
  TYPE(georef_coord) :: coord
  REAL :: current=rmiss
  CHARACTER(len=2) :: unit=''
  CHARACTER(len=2) :: lighttype=''
  CONTAINS
  PROCEDURE :: readline => lightning_readline
  PROCEDURE :: c_e => lightning_c_e
END TYPE cnmc_lightning_t

PRIVATE
PUBLIC cnmc_lightning_t

CONTAINS

SUBROUTINE lightning_readline(this, buffer)
CLASS(cnmc_lightning_t),INTENT(out) :: this
CHARACTER(len=*) :: buffer

TYPE(cnmc_lightning_t) :: lbuff
INTEGER :: ierr
REAL(kind=fp_d) :: lon, lat

! initialise local object
lbuff%time = datetime_new(isodate=buffer(:19))
READ(buffer, '(20X,F7.0,1X,F8.0,1X,F7.0,1X,A,1X,A)', iostat=ierr) lat, lon, &
 lbuff%current, lbuff%unit, lbuff%lighttype

IF (ierr == 0 .AND. c_e(lbuff%time)) THEN
  lbuff%coord = georef_coord_new(x=lon, y=lat)
  SELECT TYPE(this)
  TYPE is (cnmc_lightning_t)
    this = lbuff
  END SELECT
ENDIF

END SUBROUTINE lightning_readline


FUNCTION lightning_c_e(this)
CLASS(cnmc_lightning_t),INTENT(in) :: this
LOGICAL :: lightning_c_e

lightning_c_e = c_e(this%time) .AND. c_e(this%coord) .AND. c_e(this%current)

END FUNCTION lightning_c_e

END MODULE cnmc_lightning_mo
