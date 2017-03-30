PROGRAM prodsim_v7d_fast_read
! This program reads a csv file in the dball-e format with the
! assumption that it contains observations all having the same
! metadata except coordinate and writes a v7d native file, this allows
! a much faster read than traditional bufr format for big datasets. It
! accepts only 2 parameters: input and output file.
USE file_utilities
USE vol7d_class
IMPLICIT NONE

INTEGER :: nl, ierr
CHARACTER(len=512) :: input_file, output_file, line
CHARACTER(len=20) :: date, btable, network
DOUBLE PRECISION :: lon, lat
TYPE(csv_record) :: csv_line
TYPE(vol7d) :: v7d

CALL getarg(1, input_file)
CALL getarg(2, output_file)

OPEN(10, file=input_file)

!READ(10,*)
nl = 0
DO WHILE(.TRUE.)
  READ(10,*, iostat=ierr)
  IF (ierr /= 0) EXIT
  nl = nl + 1
ENDDO
PRINT*,'Scorso ',TRIM(input_file),nl,'record'

CALL init(v7d)
CALL vol7d_alloc(v7d, nana=nl, nnetwork=1, ntimerange=1, nlevel=1, &
 ndativarr=1)
CALL vol7d_alloc_vol(v7d)
REWIND(10)
!READ(10,*)
nl = 0
DO WHILE(.TRUE.)
  READ(10,'(A)', iostat=ierr) line
  IF (ierr /= 0) EXIT
  nl = nl + 1
  CALL init(csv_line, line)
  CALL csv_record_getfield(csv_line, lon, ierr)
  CALL csv_record_getfield(csv_line, lat, ierr)
  CALL init(v7d%ana(nl),lon=lon, lat=lat)
  IF (nl == 1) THEN
! network
    CALL csv_record_getfield(csv_line, network)
    CALL init(v7d%network(1), network)
! datetime
    CALL csv_record_getfield(csv_line, date)
    CALL init(v7d%time(1), isodate=date)
! level
    CALL csv_record_getfield(csv_line, v7d%level(1)%level1)
    CALL csv_record_getfield(csv_line, v7d%level(1)%l1)
    CALL csv_record_getfield(csv_line, v7d%level(1)%level2)
    CALL csv_record_getfield(csv_line, v7d%level(1)%l2 )
! timerange
    CALL csv_record_getfield(csv_line, v7d%timerange(1)%timerange)
    CALL csv_record_getfield(csv_line, v7d%timerange(1)%p1)
    CALL csv_record_getfield(csv_line, v7d%timerange(1)%p2)
! variable
    CALL csv_record_getfield(csv_line, btable)
    CALL init(v7d%dativar%r(1), btable=btable)
  ELSE
! network
    CALL csv_record_getfield(csv_line)
! datetime
    CALL csv_record_getfield(csv_line)
! level
    CALL csv_record_getfield(csv_line)
    CALL csv_record_getfield(csv_line)
    CALL csv_record_getfield(csv_line)
    CALL csv_record_getfield(csv_line)
! timerange
    CALL csv_record_getfield(csv_line)
    CALL csv_record_getfield(csv_line)
    CALL csv_record_getfield(csv_line)
! variable
    CALL csv_record_getfield(csv_line)
  ENDIF
  CALL csv_record_getfield(csv_line, v7d%voldatir(nl,1,1,1,1,1))
  CALL delete(csv_line)
ENDDO
PRINT*,'Letto ',TRIM(input_file),nl,'record'

CALL export(v7d, filename=output_file)
PRINT*,'Scritto ',TRIM(output_file)


END PROGRAM prodsim_v7d_fast_read
