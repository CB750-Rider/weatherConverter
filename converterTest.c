/* converterTest.c
Copyright (C) 2016 Dr. Lee R. Burchett

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE. 
 *
 *  Created on: Aug 11, 2016
 *      Author: Dr. Lee R. Burchett
 
 converter Test is a short command line utility that is used to verify that the weather 
 conversion library functions as intended.
 
 */

#include "weatherConversion.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

void importFile(const char *fname, WEATHER_CONVERSION_VECTOR *OUT);
double compare(WEATHER_CONVERSION_VECTOR *A, WEATHER_CONVERSION_VECTOR *B);
double relError(double *a, double *b, unsigned int N);
double _relError(double a, double b);
void setTestVector(WEATHER_CONVERSION_VECTOR *TST, WEATHER_CONVERSION_VECTOR *STD, WEATHER_CONVERTER_FIELD field);

int main(){
	WEATHER_CONVERSION_VECTOR STANDARD;
	WEATHER_CONVERSION_VECTOR TEST;

	importFile("humidityTest.csv",&STANDARD);
	openWeatherConversionVector(&TEST, STANDARD.N);

	WEATHER_CONVERTER_FIELD ri;

	for(ri=_RELATIVE_HUMIDITY;ri<_MOIST_AIR_DENSITY;ri++){
		printf("Now starting with %s.\n",_weather_converter_field_names[ri]);
		setTestVector(&TEST,&STANDARD,ri);
		printf("MEAN ABSOLUTE RELATIVE ERROR = %g.\n",compare(&TEST,&STANDARD));
	}

	return 0;
}

void importFile(const char *fname, WEATHER_CONVERSION_VECTOR *OUT){
	FILE *fp;
	fp = fopen(fname,"r");
	uint N,i;
	char line[2000],*dmp;

	for(N=0;feof(fp)==0;N++)
		if(fgets(line,2000,fp)==NULL){
			printf("converterTest.c:importFile() Successfully read %d lines in %s.\n",N,fname);
			break;
		}

	rewind(fp);

	openWeatherConversionVector(OUT, N);
	for(i=0;feof(fp)==0;i++){
		if (i==N) break;
		if(fgets(line,2000,fp)==NULL)
			printf("converterTest.c:importFile() Error reading line %d in %s.\n",i,fname);
		OUT->val[_TEMPERATURE_K][i]=strtod(line,&dmp);
		OUT->val[_PRESSURE][i]=strtod(dmp,&dmp);
		OUT->val[_RELATIVE_HUMIDITY][i]=strtod(dmp,&dmp);
		OUT->val[_ABSOLUTE_HUMIDITY][i]=strtod(dmp,&dmp);
		OUT->val[_DEW_POINT_K][i]=strtod(dmp,&dmp);
		OUT->val[_ENHANCEMENT_FACTOR][i]=strtod(dmp,&dmp);
		OUT->val[_MASS_MIXING_RATIO][i]=strtod(dmp,&dmp);
		OUT->val[_MOIST_AIR_DENSITY][i]=strtod(dmp,&dmp);
		OUT->val[_MOLE_MIXING_RATIO][i]=strtod(dmp,&dmp);
		OUT->val[_POTENTIAL_TEMPERATURE][i]=strtod(dmp,&dmp);
		OUT->val[_SATURATION_MIXING_RATIO][i]=strtod(dmp,&dmp);
		OUT->val[_SATURATION_VAPOR_PRESSURE][i]=strtod(dmp,&dmp);
		OUT->val[_SPECIFIC_HUMIDITY][i]=strtod(dmp,&dmp);
		OUT->val[_VAPOR_PRESSURE][i]=strtod(dmp,&dmp);
		OUT->val[_VIRTUAL_POTENTIAL_TEMPERATURE][i]=strtod(dmp,&dmp);
		OUT->val[_VIRTUAL_TEMPERATURE][i]=strtod(dmp,&dmp);
		OUT->val[_V_WIND][i] = 1.0;
		OUT->val[_U_WIND][i] = 1.0;
	}
	OUT->populated[_TEMPERATURE_K]=TRUE;
	OUT->populated[_PRESSURE]=TRUE;
	OUT->populated[_RELATIVE_HUMIDITY]=TRUE;
	OUT->populated[_ABSOLUTE_HUMIDITY]=TRUE;
	OUT->populated[_DEW_POINT_K]=TRUE;
	OUT->populated[_ENHANCEMENT_FACTOR]=TRUE;
	OUT->populated[_MASS_MIXING_RATIO]=TRUE;
	OUT->populated[_MOIST_AIR_DENSITY]=TRUE;
	OUT->populated[_MOLE_MIXING_RATIO]=TRUE;
	OUT->populated[_POTENTIAL_TEMPERATURE]=TRUE;
	OUT->populated[_SATURATION_MIXING_RATIO]=TRUE;
	OUT->populated[_SATURATION_VAPOR_PRESSURE]=TRUE;
	OUT->populated[_SPECIFIC_HUMIDITY]=TRUE;
	OUT->populated[_VAPOR_PRESSURE]=TRUE;
	OUT->populated[_VIRTUAL_POTENTIAL_TEMPERATURE]=TRUE;
	OUT->populated[_VIRTUAL_TEMPERATURE]=TRUE;
	OUT->populated[_U_WIND]=TRUE;
	OUT->populated[_V_WIND]=TRUE;

	setAllFields(OUT);
	/* Fill in the remaining fields */
	fclose(fp);
}

double compare(WEATHER_CONVERSION_VECTOR *TST, WEATHER_CONVERSION_VECTOR *STD){
	double totalErr = 0.0,tmpError;
	WEATHER_CONVERTER_FIELD ri;

	for(ri=0;ri<_N_WEATHER_FIELDS;ri++){
		if((tmpError=relError(TST->val[ri],STD->val[ri],TST->N)) > 1.0e-6)
			printf("A high relative error, %lf, was seen for field # %d, %s.\n",tmpError,ri,_weather_converter_field_names[ri]);
		totalErr += tmpError;
	}

	return totalErr/(double)TST->N/(double)_N_WEATHER_FIELDS;
}

double relError(double *a, double *b, unsigned int N){
	/* Recursive add relative error.*/
	int brk;

	switch(N){
	case 0:
		return 0.0;
		break;
	case 1:
		return _relError(a[0],b[0]);
	break;
	case 2:
		return _relError(a[0],b[0]) + _relError(a[1],b[1]);
		break;
	case 3:
		return _relError(a[0],b[0]) + _relError(a[1],b[1]) + _relError(a[2],b[2]);
		break;

	case 4:
		return _relError(a[0],b[0]) + _relError(a[1],b[1]) + _relError(a[2],b[2]) + _relError(a[3],b[3]);
		break;
	default:
		brk = N/2;
		return relError(a,b,brk) + relError(a+brk,b+brk,N-brk);
		break;
	}
}

double _relError(double a, double b){
	if((a==0.0)&&(b==0.0)) return 0.0;
	return pow((a-b) / (fabs(a)>fabs(b)?a:b),2.0);
}
void setTestVector(WEATHER_CONVERSION_VECTOR *TST, WEATHER_CONVERSION_VECTOR *STD, WEATHER_CONVERTER_FIELD field){
	WEATHER_CONVERTER_FIELD ri;
	for(ri=0;ri<_N_WEATHER_FIELDS;ri++)
		TST->populated[ri]=FALSE;
	setWeatherField(TST, STD->val[_U_WIND], _U_WIND);
	setWeatherField(TST, STD->val[_V_WIND], _V_WIND);
	setWeatherField(TST, STD->val[_TEMPERATURE_K], _TEMPERATURE_K);
	setWeatherField(TST, STD->val[_PRESSURE], _PRESSURE);
	setWeatherField(TST, STD->val[field], field);
	setAllFields(TST); /* Attempt to set fields*/
}
