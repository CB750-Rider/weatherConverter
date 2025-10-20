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
#include <errno.h>
#include <string.h>

void importFile(const char *fname, WEATHER_CONVERSION_VECTOR *OUT);
void saveToFile(WEATHER_CONVERSION_VECTOR *V, const char *fname);
void setUpTest(WEATHER_CONVERSION_VECTOR *V);
double compare(WEATHER_CONVERSION_VECTOR *A, WEATHER_CONVERSION_VECTOR *B);
double relError(double *a, double *b, unsigned int N);
void setTestVector(WEATHER_CONVERSION_VECTOR *TST, WEATHER_CONVERSION_VECTOR *STD, WEATHER_CONVERTER_FIELD field);

#define uint size_t

static 	WEATHER_CONVERTER_FIELD WC_field_list[] = {
		_TEMPERATURE_K,
		_PRESSURE,
		_RELATIVE_HUMIDITY,
		_ABSOLUTE_HUMIDITY,
		_DEW_POINT_K,
		_ENHANCEMENT_FACTOR,
		_MASS_MIXING_RATIO,
		_MOIST_AIR_DENSITY,
		_MOLE_MIXING_RATIO,
		_POTENTIAL_TEMPERATURE,
		_SATURATION_MIXING_RATIO,
		_SATURATION_VAPOR_PRESSURE,
		_SPECIFIC_HUMIDITY,
		_VAPOR_PRESSURE,
		_VIRTUAL_POTENTIAL_TEMPERATURE,
		_VIRTUAL_TEMPERATURE};

int main(){
	WEATHER_CONVERSION_VECTOR STANDARD;
	WEATHER_CONVERSION_VECTOR TEST;
	WEATHER_CONVERTER_FIELD ri;

	printf("Running tests to see how well the conversions work.\n");

	printf("Performing a self-consistency test. (Can we go forward and back).\n");
	setUpTest(&STANDARD);
	printf("Test is now set up. Opening the conversion vector.\n");
	openWeatherConversionVector(&TEST, STANDARD.N);
	printf("Vector open was successful.\n");
	for(ri=_RELATIVE_HUMIDITY;ri<_MOIST_AIR_DENSITY;ri++){
		printf("Now performing the conversion starting with temperature, pressure, and %s.\n",_weather_converter_field_names[ri]);
		setTestVector(&TEST,&STANDARD,ri);
		printf("MEAN SQUARED RELATIVE ERROR = %g.\n",compare(&TEST,&STANDARD));
	}
	freeWeatherConversionVector(&STANDARD);
	freeWeatherConversionVector(&TEST);

	printf("Performing a comparison to standard data.\n");
	importFile("humidityTest.csv",&STANDARD);
	openWeatherConversionVector(&TEST, STANDARD.N);
	for(ri=_RELATIVE_HUMIDITY;ri<_MOIST_AIR_DENSITY;ri++){
		printf("Now performing the conversion starting with temperature, pressure, and %s.\n",_weather_converter_field_names[ri]);
		setTestVector(&TEST,&STANDARD,ri);
		printf("MEAN SQUARED RELATIVE ERROR = %g.\n",compare(&TEST,&STANDARD));
	}
	freeWeatherConversionVector(&STANDARD);
	freeWeatherConversionVector(&TEST);

	return 0;
}

void setUpTest(WEATHER_CONVERSION_VECTOR *V){
	double P[] = {1.,5.,10.,15.,20.,40.,80.,100.,200.,300.,400.,500.,600.,700.,800.,900.,950.,1000.,1050.,1100.,1150.,1200.,1250.,1300.}; /* in mb */
	double T=250.,_T[] = {250.,5.,320.}; /* Temperature range (min,step,max) */
	double RH=0.,_RH[] = {0.1,4.99,100.}; /* Relative Humidity Range */

	size_t iP,iT,iRH,idx,NP=sizeof(P)/sizeof(double);
	size_t NT = (size_t) ( (_T[2]-_T[0])/_T[1]);
	size_t NRH = (size_t) ( (_RH[2]-_RH[0])/_RH[1]);

	printf("Opening the vector.");
	openWeatherConversionVector(V, (unsigned int)NP*NT*NRH);
	printf(" Done.")
	for(iP=0;iP<NP;iP++){
		for(iT=0;iT<NT;iT++){
			T = _T[0] + (double)iT*_T[1];
			for(iRH=0;iRH<NRH;iRH++){
				RH = _RH[0] + (double)iRH*_RH[1];
				idx = iP*NT*NRH + iT*NRH + iRH;
				V->val[_TEMPERATURE_K][idx] = T;
				V->val[_PRESSURE][idx] = P[iP];
				V->val[_RELATIVE_HUMIDITY][idx] = RH;
				V->val[_U_WIND][idx] = 1.0;
				V->val[_V_WIND][idx] = 1.0;
			}
		}
	}

	V->populated[_TEMPERATURE_K]=TRUE;
	V->populated[_PRESSURE]=TRUE;
	V->populated[_RELATIVE_HUMIDITY]=TRUE;
	V->populated[_U_WIND]=TRUE;
	V->populated[_V_WIND]=TRUE;

	setAllFields(V);

	saveToFile(V,"temporaryHumidityTest.csv");
}
void saveToFile(WEATHER_CONVERSION_VECTOR *V, const char *fname){
	FILE *fp;
	fp = fopen(fname,"w");
	uint i,j,N=(uint)sizeof( WC_field_list)/sizeof(WEATHER_CONVERTER_FIELD);
	WEATHER_CONVERTER_FIELD fi;

	if(fp==NULL){
		printf("Unable to open %s for writing.\n",fname);
		return;
	}

	for(i=0;i<V->N;i++){
		fi = WC_field_list[0];
		fprintf(fp,"%18.16lf",V->val[fi][i]);
		for(j=1;j<N;j++){
			fi = WC_field_list[j];
			fprintf(fp,",%18.16lf",V->val[fi][i]); /* <- added a comma */
		}
		fprintf(fp,"\n");
	}
	fclose(fp);
}
void importFile(const char *fname, WEATHER_CONVERSION_VECTOR *OUT){
	FILE *fp;
	fp = fopen(fname,"r");
	uint N, i;
	char line[2000],*dmp;
	uint j,NF=(uint)sizeof( WC_field_list)/sizeof(WEATHER_CONVERTER_FIELD);
	WEATHER_CONVERTER_FIELD fi;

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
		dmp = line;
		for(j=0;j<NF;j++){/* Iterate over all the standard input fields. */
			fi = WC_field_list[j];
			if((OUT->val[fi][i]=strtod(dmp,&dmp))==0.0){
				if(errno==EINVAL){
					fprintf(stderr,"Error %s reported when attempting to convert %s.\n",strerror(errno),dmp);
					exit(-3);
				}
			}
			if(dmp[0]==',')dmp++;
			if(dmp[0]=='\n' || dmp[0]=='\0')break;
		}
		OUT->val[_V_WIND][i] = 1.0;
		OUT->val[_U_WIND][i] = 1.0;
	}
	for(j=0;j<NF;j++){
		fi = WC_field_list[j];
		OUT->populated[fi]=TRUE;
	}
	OUT->populated[_U_WIND] = OUT->populated[_V_WIND]=TRUE;

	setAllFields(OUT);
	/* Fill in the remaining fields */
	fclose(fp);
}

double compare(WEATHER_CONVERSION_VECTOR *TST, WEATHER_CONVERSION_VECTOR *STD){
	double totalErr = 0.0,tmpError;
	WEATHER_CONVERTER_FIELD ri;

	for(ri=0;ri<_N_WEATHER_FIELDS;ri++){
		if(ri==_OTHER_INPUT)continue;
		if((tmpError=relError(TST->val[ri],STD->val[ri],TST->N)) > ((double)TST->N)*1.0e-6)
			printf("A high relative error, %lf, was seen for field # %d, %s.\n",tmpError/((double)TST->N),ri,_weather_converter_field_names[ri]);
		totalErr += tmpError;
	}

	return totalErr/(double)_N_WEATHER_FIELDS/(double)TST->N;
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
