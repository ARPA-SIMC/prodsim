! Program to recode cosmo grib2 output to fit premoloch input
! conventions. It should be called with the cosmo constant file
! (e.g. lfff00000000c) as first argument and with all the desired
! cosmo prognostic files (e.g.lfff????0000) as following arguments
PROGRAM cosmo_to_premoloch
USE cosmo_to_premoloch_mo
IMPLICIT NONE

CHARACTER(len=512) :: input_file, output_file
INTEGER :: ifid, ofid, i, ier

DO i = 2, iargc()
! open output file grib_002
  WRITE(output_file,'(''grib_'',I3.3)')i-1
  CALL grib_open_file(ofid, output_file, 'w', ier)
  IF (ier /= GRIB_SUCCESS) THEN
    PRINT*,'Error opening output file ',TRIM(output_file)
    STOP
  ENDIF

  IF (i == 2) THEN ! digression for processing constant data
    CALL getarg(1, input_file)
    CALL grib_open_file(ifid, input_file, 'r', ier)
    IF (ier /= GRIB_SUCCESS) THEN
      PRINT*,'Error opening constant file ',TRIM(input_file)
      STOP
    ENDIF
    CALL process_const_fields(ifid, ofid)
    CALL grib_close_file(ifid)
  ENDIF

! process prognostic data
  CALL getarg(i, input_file)
  CALL grib_open_file(ifid, input_file, 'r', ier)
  IF (ier /= GRIB_SUCCESS) THEN
    PRINT*,'Error opening input file ',TRIM(input_file)
    STOP
  ENDIF
  CALL process_progn_fields(ifid, ofid)

  CALL grib_close_file(ifid)
  CALL grib_close_file(ofid)
ENDDO

END PROGRAM cosmo_to_premoloch
