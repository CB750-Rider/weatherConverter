/* conversionVectorHelpers.c
 *
 * Copyright (C) 2016 Dr. Lee Burchett
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE. 
 *
 *  Created on: Aug 11, 2016
 *      Author: Dr. Lee Burchett 
 *
 *  Last Modified: 6 June, 2017
 * 
 * Thes file contains tools for interacting with routines that are useful 
 * for converting the myriad ways that researchers use to represent humidity in 
 * the atmosphere. The goal is to provide a tool to save time on writing and 
 * debugging these conversions. If any conversion data appear suspect. Please 
 * contact me with suggested changes at lee.r.burchett@gmail.com.
 *

 */

#include "weatherConversion.h"
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <errno.h>

char *_weather_converter_error_messages[N_WEATHER_CONVERSION_ERRORS] = {
	"SUCCESS (NO ERROR REPORTED)",
	"TEMPERATURE DATA ARE REQUIRED FOR CONVERSION BUT ARE NOT AVAILABLE",
	"PRESSURE DATA ARE REQUIRED FOR CONVERSION BUT ARE NOT AVAILABLE",
	"HUMIDITY DATA ARE REQUIRED FOR CONVERSION BUT ARE NOT AVAILABLE",
	"WIND DATA ARE REQUIRED FOR CONVERSION BUT ARE NOT AVAILABLE",
	"THE REQUESTED FIELD HAS NOT BEEN ALLOCATED",
	"A FILE INPUT/OUTPUT ERROR WAS DETECTED",
	"A VECTOR THAT HAS ALREADY BEEN INITIALIZED WAS INITIALIZED AGAIN, AND MAY CAUSE A MEMORY LEAK.",
	"AN ERROR WAS ENCOUNTERED WHEN ALLOCATING MEMORY",
	"AN UNKNOWN OR UNSPECIFIED ERROR WAS ENCOUNTERED",
    "AN UNRECOVERABLE ERROR WAS ENCOUNTERED DURING PROGRAM INITIALIZATION"};
/* Keep track of any initialized vectors. */
static unsigned int N_vector = 0;
static WEATHER_CONVERSION_VECTOR **initialized_vectors=NULL;
FILE *weather_converter_error_stream = NULL;

/* Initialization tracking functions */
static BOOLEAN isAlreadyInitialized(WEATHER_CONVERSION_VECTOR *WX){
	unsigned int i;
	for(i=0;i<N_vector;i++){
		if(WX == initialized_vectors[i]) return TRUE;
	}
	return FALSE;
}
static WEATHER_CONVERSION_ERROR add_vector_to_init_list(WEATHER_CONVERSION_VECTOR *WX){
	/* create space for copying to */
	initialized_vectors = realloc(initialized_vectors,sizeof(WEATHER_CONVERSION_VECTOR*)*(N_vector+1));
	if (initialized_vectors == NULL)return ALLOCATION_ERROR;
	/* add the new entry */
	initialized_vectors[N_vector] = WX;
	N_vector++;
	return WEATHER_CONVERSION_SUCCESS;
}
static void remove_vector_from_init_list(WEATHER_CONVERSION_VECTOR *WX){
	/* Search through the list for this vector, and remove it from the list. */
	unsigned int i,j;
	for(i=0;i<N_vector;i++){
		if(initialized_vectors[i]==WX){
			for(j=i+1;j<N_vector;j++)
				initialized_vectors[j-1] = initialized_vectors[j];
			break;
		}
	}
	N_vector--;
	if(N_vector>0)
		initialized_vectors = realloc(initialized_vectors,sizeof(WEATHER_CONVERSION_VECTOR*)*(N_vector));
	else{
		free(initialized_vectors);
		initialized_vectors=NULL;
	}
}
WEATHER_CONVERTER_FIELD getEnumFromFlag(const char *flag) {
    WEATHER_CONVERTER_FIELD i;
    for(i=0;i<_N_WEATHER_FIELDS;i++){
        if (strcmp(flag, _weather_converter_field_flags[i]) == 0)
            return i;
    }
    return -1;
}
const char *getFieldFlag(WEATHER_CONVERTER_FIELD wcf){
    if((wcf < 0) || (wcf >= _N_WEATHER_FIELDS))
        return "<UNKNOWN FIELD ENUM>";
    return _weather_converter_field_flags[wcf];
}
const char *getFieldName(WEATHER_CONVERTER_FIELD wcf){
    if((wcf < 0) || (wcf >= _N_WEATHER_FIELDS))
        return "<UNKNOWN FIELD ENUM>";
    return _weather_converter_field_names[wcf];
}
const char *getFieldUnitsAbbr(WEATHER_CONVERTER_FIELD wcf){
    if((wcf < 0) || (wcf >= _N_WEATHER_FIELDS))
        return "<UNKNOWN FIELD ENUM>";
    return _weather_converter_field_units[wcf];
}
int get_N_WEATHER_FIELDS() {
    return _N_WEATHER_FIELDS;
}
void reportWeatherConversionError(char *msg, ...){
	/* This function prints errors to stdout by default, but this will also
	work with other streams. */

	va_list ap;

	va_start(ap,msg);

	if (weather_converter_error_stream==NULL)
		weather_converter_error_stream = stdout;
	else {
		if (ftell(weather_converter_error_stream) == -1) {
			/* There may be a problem with the output stream. Check errno to be sure. */
			if (errno == EBADF) { /* Ignore errors from piping output and the like. */
				printf("\n***ERROR*** Cannot write error message:\n");
				vprintf("msg", ap);
				printf("\nto the error output stream because the file is not open.\n"
					"Additional errors will be redirected to STDOUT.\n***END ERROR***.\n");
				weather_converter_error_stream = stdout;
				return;
			}
		}
	}
		

	vfprintf(weather_converter_error_stream,msg,ap);
	va_end(ap);
}
WEATHER_CONVERSION_ERROR openWeatherConversionVector(WEATHER_CONVERSION_VECTOR *WX, unsigned int N){
	/* Allocate space for the weatherConversionVector. */
	WEATHER_CONVERTER_FIELD fi;
	_wcCheck(initializeVector(WX));
	WX->N = N;
	WX->standardPressure = _weather_converter_site_defaults[_STANDARD_PRESSURE];
	WX->xCO2 = _weather_converter_site_defaults[_XCO2];
	WX->surfaceHeight = _weather_converter_site_defaults[_SURFACE_HEIGHT];
	WX->surfacePressure = _weather_converter_site_defaults[_SURFACE_PRESSURE];
	WX->latitude = _weather_converter_site_defaults[_SITE_LATITUDE];
	WX->surfaceTemperature = _weather_converter_site_defaults[_SURFACE_TEMPERATURE];
	for(fi=0;fi<_N_WEATHER_FIELDS;fi++){
		WX->val[fi] = (double *)malloc(N*sizeof(double));
		WX->allocated[fi] = TRUE;
		WX->populated[fi] = FALSE;
	}
	WX->quiet = FALSE;
	return WEATHER_CONVERSION_SUCCESS;
}
WEATHER_CONVERSION_ERROR freeWeatherConversionVector(WEATHER_CONVERSION_VECTOR *WX){
	WEATHER_CONVERTER_FIELD fi;
	for(fi=0;fi<_N_WEATHER_FIELDS;fi++){
		if(WX->allocated[fi])
			free(WX->val[fi]);
		WX->allocated[fi] = FALSE;
	}
	remove_vector_from_init_list(WX);
	return WEATHER_CONVERSION_SUCCESS;
}
WEATHER_CONVERSION_ERROR setWeatherField(WEATHER_CONVERSION_VECTOR *WX, double *x, WEATHER_CONVERTER_FIELD fi){
	unsigned int i;
	if(WX->allocated[fi]){
		for(i=0;i<WX->N;i++)
			WX->val[fi][i] = x[i];
		WX->populated[fi] = TRUE;
		return WEATHER_CONVERSION_SUCCESS;
	}
	printf("weatherConversion.h:setWeatherField() Field has not been allocated.\n");
	return FIELD_NOT_ALLOCATED;
}
WEATHER_CONVERSION_ERROR initializeVector(WEATHER_CONVERSION_VECTOR *WX){
	/* Set up an empty and open weather conversion vector. This should only be
	 * called once for each vector. */
	WEATHER_CONVERSION_ERROR rv = WEATHER_CONVERSION_SUCCESS;
	WEATHER_CONVERTER_FIELD i;
	if(isAlreadyInitialized(WX)) rv = DOUBLE_INITIALIZATION_ERROR;
	WX->N=0;
	WX->standardPressure=0.0;
	WX->xCO2=0.0;
	for(i=0;i<_N_WEATHER_FIELDS;i++){
		WX->val[i]=NULL;
		WX->allocated[i]=FALSE;
		WX->populated[i]=FALSE;
	}
	/*
	WX->f.alloc = (WEATHER_CONVERSION_ERROR (*)(void *, unsigned int ))&openWeatherConversionVector;
	WX->f.free = (WEATHER_CONVERSION_ERROR (*)(void *))&freeWeatherConversionVector;
	WX->f.clear_all = (WEATHER_CONVERSION_ERROR (*)(void *))&clearAllWeatherConversionFields;
	WX->f.change_temp =  (WEATHER_CONVERSION_ERROR (*) (void *,double *,WEATHER_CONVERTER_FIELD))&changeTemperatures;
	WX->f.change_pressure = (WEATHER_CONVERSION_ERROR (*) (void *,double *,WEATHER_CONVERTER_FIELD))&changePressures;
	WX->f.change_humidity = (WEATHER_CONVERSION_ERROR (*) (void *,double *,WEATHER_CONVERTER_FIELD))&changeHumidity;
	WX->f.set_field = (WEATHER_CONVERSION_ERROR (*) (void *,double *,WEATHER_CONVERTER_FIELD))&setWeatherField;
	WX->f.run_conversion = (WEATHER_CONVERSION_ERROR (*)(void *))&setAllFields;
	WX->f.integrate_column_water_density = (double (*)(void *, double, double))&integrateColumnWaterDensity;
	WX->f.integrate_column_water_number_density = (double(*)(void *, double, double))&integrateColumnWaterNumberDensity;
	WX->f.integrate_column_moist_air_density = (double(*)(void *, double, double))&integrateColumnMoistAirDensity;
	WX->f.integrate_column_moist_air_number_density = (double(*)(void *, double, double))&integrateColumnMoistAirNumberDensity;
	TODO Bring this back or take it out */
	if(rv==DOUBLE_INITIALIZATION_ERROR)return rv;
	_wcCheck(add_vector_to_init_list(WX));
	return WEATHER_CONVERSION_SUCCESS;
}
WEATHER_CONVERSION_ERROR clearAllWeatherConversionFields(WEATHER_CONVERSION_VECTOR *WX){
	WEATHER_CONVERTER_FIELD fi;
	/* If this vector is uninitialized, there is little that can be done.*/
	if(!isAlreadyInitialized(WX)){
		_wcCheck(initializeVector(WX));
		return FIELD_NOT_ALLOCATED;
	}
	/* Make sure the allocations look okay. If so, then set the populated field
	 * to FALSE so that the data will be ignored. */
	for(fi=0;fi<_N_WEATHER_FIELDS;fi++){
		if(WX->allocated[fi] && (WX->val[fi] != NULL))
			WX->populated[fi]=FALSE;
		else
			return FIELD_NOT_ALLOCATED;
	}
	return WEATHER_CONVERSION_SUCCESS;
}
WEATHER_CONVERSION_ERROR changeTemperatures(WEATHER_CONVERSION_VECTOR *WX, double *t, WEATHER_CONVERTER_FIELD fi){
	/* Clear all temperatures, set the field and reprocess. */
	WX->populated[_TEMPERATURE_C] = WX->populated[_TEMPERATURE_K] =
			WX->populated[_TEMPERATURE_F] = WX->populated[_POTENTIAL_TEMPERATURE] =
					FALSE;

	_wcCheck(setWeatherField(WX,t,fi));

	_wcCheck(setAllFields(WX));

	return WEATHER_CONVERSION_SUCCESS;
}
WEATHER_CONVERSION_ERROR changePressures(WEATHER_CONVERSION_VECTOR *WX, double *p,WEATHER_CONVERTER_FIELD fi){
	/* Clear all fields used to get the pressure, reset, and reprocess. */
	WX->populated[_PRESSURE] = WX->populated[_POTENTIAL_TEMPERATURE] = FALSE;

	_wcCheck(setWeatherField(WX,p,fi));

	_wcCheck(setAllFields(WX));

	return WEATHER_CONVERSION_SUCCESS;
}
WEATHER_CONVERSION_ERROR changeHumidity(WEATHER_CONVERSION_VECTOR *WX, double *h,WEATHER_CONVERTER_FIELD fi){
	/* Clear out all possible humidity fields */
	WEATHER_CONVERTER_FIELD i, hflds[] = {_RELATIVE_HUMIDITY, _VAPOR_PRESSURE,
			_POTENTIAL_VAPOR_PRESSURE, _MOLE_MIXING_RATIO, _MASS_MIXING_RATIO,
			_DEW_POINT_C,_DEW_POINT_K,_DEW_POINT_F,_SPECIFIC_HUMIDITY,
			_ABSOLUTE_HUMIDITY, _MOIST_AIR_DENSITY};
	for(i=0;i<(sizeof(hflds)/sizeof(WEATHER_CONVERTER_FIELD));i++){
		WX->populated[hflds[i]] = FALSE;
	}

	_wcCheck(setWeatherField(WX,h,fi));

	_wcCheck(setAllFields(WX));

	return WEATHER_CONVERSION_SUCCESS;
}
WEATHER_CONVERSION_ERROR changeHeight(WEATHER_CONVERSION_VECTOR *WX, double *h,WEATHER_CONVERTER_FIELD fi){
	/* Clear out all possible humidity fields */
	WEATHER_CONVERTER_FIELD i, hflds[] = {_GEOPOTENTIAL_HEIGHT, _HEIGHT_AGL,
	_HEIGHT_AMSL};
	for(i=0;i<(sizeof(hflds)/sizeof(WEATHER_CONVERTER_FIELD));i++){
		WX->populated[hflds[i]] = FALSE;
	}

	_wcCheck(setWeatherField(WX,h,fi));

	_wcCheck(setAllFields(WX));

	return WEATHER_CONVERSION_SUCCESS;
}
WEATHER_CONVERSION_ERROR parseSiteSettingsLine(char *line, WEATHER_CONVERSION_VECTOR *WX){
	/* Convenience function to read a line and look for any setting flags for the site location. */
	SITE_SPECIFIC_SETTINGS si;
	char *sp,*dmp;
	double new_value;
	static BOOLEAN surface_height_set = FALSE;
	static BOOLEAN surface_pressure_set = FALSE;

	if (line == NULL) {
		/* Called when no more data are being read in. */
		if (surface_height_set) {
			if (surface_pressure_set) {
				/* We are done. */
				return WEATHER_CONVERSION_SUCCESS;
			} /* If both are set, then we are done. */
			else{
				/* Figure out the pressure from the height. */
				WX->surfacePressure = standardAtmosPressureAtAltitude(WX->surfaceHeight);
				return WEATHER_CONVERSION_SUCCESS;
			} /* If the surface pressure wasn't set. */
		}
		else {
			if (surface_pressure_set) {
				WX->surfaceHeight = standardAtmosAltitudeAtPressure(WX->surfacePressure);
				return WEATHER_CONVERSION_SUCCESS;
			}
			else {
				/* Do our best guess */
				if (WX->populated[_PRESSURE])
					WX->surfacePressure = WX->val[_PRESSURE][0];
				else
					WX->surfacePressure = 1013.25;
				WX->surfaceHeight = standardAtmosAltitudeAtPressure(WX->surfacePressure);
				return WEATHER_CONVERSION_SUCCESS;
			} /* If the surface pressure wasn't set either. */
		} /* If the surface heights were not set. */
	} /* If there was nothing in the line, then figure out the pressures and temperatures. */ 

	for(si=0;si<_N_WEATHER_SITE_SPECIFIC_SETTINGS; si++){
		if((sp=strstr(line,_weather_converter_site_setting_flags[si])) != NULL){
			new_value = strtod(
					sp+strlen(_weather_converter_site_setting_flags[si]),
					&dmp);
			if( (errno==ERANGE || errno==EINVAL) &&
									new_value == 0.0){
								printf("Error %s encountered when attempting to set\n"
										" %s. The default value will be used instead. \n",
										strerror(errno),
										_weather_converter_field_flags[si]);
							} /* If there was a problem */
			else{
				switch(si){
				case _STANDARD_PRESSURE:
					WX->standardPressure = new_value;
					break;
				case _XCO2:
					WX->xCO2 = new_value;
					break;
				case _SITE_LATITUDE:
					if (new_value < -90.0 && new_value > 90.0)
						printf("Error in Site Latitude Value, %lf. The value should be between -90 and 90.\n"
								"The default value will be used instead.\n",new_value);
					else
						WX->latitude = new_value;
					break;
				case _SURFACE_HEIGHT:
					WX->surfaceHeight = new_value;
					surface_height_set = TRUE;
					break;
				case _SURFACE_PRESSURE:
					WX->surfacePressure = new_value;
					surface_height_set = TRUE;
					break;
				default:
					printf("Error encountered when parsing site specific settings.\n"
							"A setting for %s was requested, but the parsing function\n"
							"does not know how to handle these values yet.\n"
							"Please see conversionVectorHelpers:parseSiteSettingLine().\n",
							_weather_converter_field_flags[si]);
					break;
				} /* End the switch*/
			} /* End if there were no errors in reading the value */
		} /* End if we found a flag in the line  */
	} /* End for each possible site setting */

	return WEATHER_CONVERSION_SUCCESS;
} 
void setQuiet(WEATHER_CONVERSION_VECTOR *WX){
    WX->quiet = TRUE;
}
void setNotQuiet(WEATHER_CONVERSION_VECTOR *WX){
    WX->quiet = FALSE;
}
void printWeatherConversionVectorMetadata(WEATHER_CONVERSION_VECTOR *WX){
    int i;
    printf("Weather Conversion Vector Metadata:\n");
    printf(" N: %u\n", WX->N);
    printf(" allocated: ");
    for(i=0;i<_N_WEATHER_FIELDS;i++) printf(", %u", WX->allocated[i]);
    printf("\n populated: ");
    for(i=0;i<_N_WEATHER_FIELDS;i++) printf(", %u", WX->populated[i]);
    printf("\n standandPressure: %lf mb\n", WX->standardPressure);
    printf(" carbon dioxide concentration: %lf ppm\n", WX->xCO2);
    printf(" latitude: %lf degrees North\n", WX->latitude);
    printf(" surfaceHeight: %lf m MSL\n", WX->surfaceHeight);
    printf(" surfacePressure: %lf mb\n", WX->surfacePressure);
    printf(" surfaceTemperature: %lf K\n", WX->surfaceTemperature);
}/* End function */
