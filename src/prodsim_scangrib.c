/*
 * Program template for scanning a single grib file and doing
 * something useful with information such as reference time, forecast
 * time, parameter, latitude, longitude, value.
 *
 * Notice: latitude and longitude are correct only with regular_ll
 * projection.
 *
 * Must be linked with GRIB-API library (-lgrib_api)
 * https://software.ecmwf.int/wiki/display/GRIB/What+is+GRIB-API
 *
 * Based on examples:
 * https://software.ecmwf.int/wiki/display/GRIB/get.c
 * https://software.ecmwf.int/wiki/display/GRIB/iterator_bitmap.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "grib_api.h"
 
int main(int argc, char** argv)
{
  FILE* in = NULL;
  char* filename;
  grib_handle *h = NULL;
  int err = 0;

  /* double *values = NULL; */
  char dataDate[9], dataTime[5];
  long stepRange, table2Version, indicatorOfParameter;
  size_t len = 0, values_len = 0, n = 0;

  size_t bmp_len = 0;
  long bitmapPresent = 0;
  long *bitmap = NULL;
  grib_iterator* iter=NULL;
  double lat,lon,value;

  filename = argv[1];
  
  in = fopen(filename,"r");
  if(!in) {
    fprintf(stderr, "ERROR: unable to open file %s\n",filename);
    return 1;
  }
 
  /* loop on messages in file */
  while ((h = grib_handle_new_from_file(0,in,&err)) != NULL ) {
    if (err != GRIB_SUCCESS) GRIB_CHECK(err,0);
  
    /* get reference date and time (dataDate, dataTime) as a char,
       computed keys from yearOfCentury, month, day, hour, minute,
       centuryOfReferenceTimeOfData */

    /* GRIB_CHECK(grib_get_length(h, "dataDate", &len), 0); */
    /* if (len > sizeof(dataDate)) { */
    /*   fprintf(stderr, "ERROR: unexpected length of dataDate %d\n",len); */
    /*   return 1; */
    /* } */
    len = sizeof(dataDate);
    grib_get_string(h,"dataDate",dataDate,&len);

    /* GRIB_CHECK(grib_get_length(h, "dataTime", &len), 0); */
    /* if (len > sizeof(dataTime)) { */
    /*   fprintf(stderr, "ERROR: unexpected length of dataTime %d\n",len); */
    /*   return 1; */
    /* } */
    len = sizeof(dataTime);
    grib_get_string(h,"dataTime",dataTime,&len);

    GRIB_CHECK(grib_get_long(h,"stepRange",&stepRange),0);
    GRIB_CHECK(grib_get_long(h,"table2Version", &table2Version), 0);
    GRIB_CHECK(grib_get_long(h,"indicatorOfParameter", &indicatorOfParameter), 0);

    GRIB_CHECK(grib_get_long(h,"bitmapPresent",&bitmapPresent),0);
    if (bitmapPresent) {
      GRIB_CHECK(grib_get_size(h,"bitmap",&bmp_len),0);
      bitmap = malloc(bmp_len*sizeof(long));
      GRIB_CHECK(grib_get_long_array(h,"bitmap",bitmap,&bmp_len),0);
      /* printf("Bitmap is present. Num = %ld\n", bmp_len); */
    }
    /* Sanity check. Number of values must match number in bitmap */
    GRIB_CHECK(grib_get_size(h,"values",&values_len),0);
    /* explicit values not needed with following iterator */
    /* values = malloc(values_len*sizeof(double)); */
    /* GRIB_CHECK(grib_get_double_array(h,"values",values,&values_len),0); */
    if (bitmapPresent) {
      assert(values_len == bmp_len);
    }
 
    /* A new iterator on lat/lon/values is created from the message handle h */
    iter = grib_iterator_new(h,0,&err);
    if (err != GRIB_SUCCESS) GRIB_CHECK(err,0);

    n = 0;
    /* Loop on all the lat/lon/values. Only print non-missing values */
    while(grib_iterator_next(iter,&lat,&lon,&value))
      {
	/* Consult bitmap to see if the n'th value is missing */
	int is_missing_val = (bitmapPresent && bitmap[n] == 0);
	if (!is_missing_val) {
	  /* do something here with
	     dataDate,dataTime,stepRange,table2Version,indicatorOfParameter */
	  if (n < 8) /* for demo */
	    printf("%s %s,%d,%d,%d,%f,%f,%f\n",dataDate,dataTime,stepRange,
		   table2Version,indicatorOfParameter,lat,lon,value);
	}
	n++;
      }
    /* Check number of elements in iterator matches value count */
    assert(n == values_len);
    /* free(values); */
 
    grib_iterator_delete(iter);
    grib_handle_delete(h);
  }

  fclose(in);
  return 0;
}


