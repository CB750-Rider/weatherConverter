/*
 * conversionVectorHelpers.c
 *
 *  Created on: Aug 11, 2016
 *      Author: Dr. Lee Burchett Booz | Allen Hamilton
 *
 *  Last Modified: 6 June, 2017
 * 
 * Thes file contains tools for interacting with routines that are useful 
 * for converting the myriad ways that researchers use to represent humidity in 
 * the atmosphere. The goal is to provide a tool to save time on writing and 
 * debugging these conversions. If any conversion data appear suspect. Please 
 * contact me with suggested changes at lee.r.burchett@gmail.com.
 *
 * Copyright (c) 2016 Dr. Lee Burchett
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
 */

#include "weatherConversion.h"
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
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
	"AN UNKNOWN OR UNSPECIFIED ERROR WAS ENCOUNTERED"};
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
	initialized_vectors = realloc(initialized_vectors,sizeof(WEATHER_CONVERSION_VECTOR*)*(N_vector));
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
	WX->standardPressure = 1000.0;
	WX->xCO2 = 390.0;
	for(fi=0;fi<_N_WEATHER_FIELDS;fi++){
		WX->val[fi] = (double *)malloc(N*sizeof(double));
		WX->allocated[fi] = TRUE;
		WX->populated[fi] = FALSE;
	}
	return WEATHER_CONVERSION_SUCCESS;
}

WEATHER_CONVERSION_ERROR freeWeatherConversionVector(WEATHER_CONVERSION_VECTOR *WX){
	WEATHER_CONVERTER_FIELD fi;
	for(fi=0;fi<_N_WEATHER_FIELDS;fi++){
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
	WX->f.alloc = (WEATHER_CONVERSION_ERROR (*)(void *, unsigned int ))&openWeatherConversionVector;
	WX->f.free = (WEATHER_CONVERSION_ERROR (*)(void *))&freeWeatherConversionVector;
	WX->f.clear_all = (WEATHER_CONVERSION_ERROR (*)(void *))&clearAllWeatherConversionFields;
	WX->f.change_temp =  (WEATHER_CONVERSION_ERROR (*) (void *,double *,WEATHER_CONVERTER_FIELD))&changeTemperatures;
	WX->f.change_pressure = (WEATHER_CONVERSION_ERROR (*) (void *,double *,WEATHER_CONVERTER_FIELD))&changePressures;
	WX->f.change_humidity = (WEATHER_CONVERSION_ERROR (*) (void *,double *,WEATHER_CONVERTER_FIELD))&changeHumidity;
	WX->f.set_field = (WEATHER_CONVERSION_ERROR (*) (void *,double *,WEATHER_CONVERTER_FIELD))&setWeatherField;
	WX->f.run_conversion = (WEATHER_CONVERSION_ERROR (*)(void *))&setAllFields;
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
			_DEW_POINT,_SPECIFIC_HUMIDITY,_ABSOLUTE_HUMIDITY, _MOIST_AIR_DENSITY};
	for(i=0;i<(sizeof(hflds)/sizeof(WEATHER_CONVERTER_FIELD));i++){
		WX->populated[hflds[i]] = FALSE;
	}

	_wcCheck(setWeatherField(WX,h,fi));

	_wcCheck(setAllFields(WX));

	return WEATHER_CONVERSION_SUCCESS;
}
