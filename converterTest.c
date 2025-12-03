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

#ifndef uint
	#define uint size_t
#endif 
int importFile(const char *fname, WEATHER_CONVERSION_VECTOR *OUT);
void saveSpeedOfSound(M_BY_N_ARRAY *A, const char *fname);
void saveToFile(WEATHER_CONVERSION_VECTOR *V, const char *fname, const WEATHER_CONVERTER_FIELD *WC_field_list, uint N);
void setUpTest(WEATHER_CONVERSION_VECTOR *V);
double compare(WEATHER_CONVERSION_VECTOR *A, WEATHER_CONVERSION_VECTOR *B, BOOLEAN *C);
void recordErrorProneConversions(WEATHER_CONVERSION_VECTOR *A, WEATHER_CONVERSION_VECTOR *B, BOOLEAN *C, WEATHER_CONVERTER_FIELD si, const char *f);
void addField(FILE *, WEATHER_CONVERSION_VECTOR *, WEATHER_CONVERTER_FIELD , const char *);
void recordSetTestVector(WEATHER_CONVERSION_VECTOR *, 
	WEATHER_CONVERSION_VECTOR *, BOOLEAN *, 
	const WEATHER_CONVERTER_FIELD, const char *);
double relError(double *a, double *b, unsigned int N);
void setTestVector(WEATHER_CONVERSION_VECTOR *TST, WEATHER_CONVERSION_VECTOR *STD, WEATHER_CONVERTER_FIELD field);

static 	WEATHER_CONVERTER_FIELD WC_field_list_part[] = {
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
#define N_PART (uint)sizeof( WC_field_list_part)/sizeof(WEATHER_CONVERTER_FIELD)
static 	WEATHER_CONVERTER_FIELD WC_field_list_full[] = {
	_TEMPERATURE_C, /* 0 Degrees C */
	_TEMPERATURE_K, /* 1 Kelvin */
	_TEMPERATURE_F, /* 2 Degrees F */
	_U_WIND, /* 3 meters / second*/
	_V_WIND,  /* 4 meters / second*/
	_WIND_SPEED,  /* 5 meters / second*/
	_WIND_DIRECTION, /* 6 degrees from N */
	_PRESSURE, /* 7 millibar */
	_POTENTIAL_TEMPERATURE, /* 8 Kelvin */
	_VIRTUAL_TEMPERATURE, /* 9 Kelvin */
	_VIRTUAL_POTENTIAL_TEMPERATURE, /* 10 Kelvin */
	_SATURATION_VAPOR_PRESSURE, /* 11 millibar */
	_SATURATION_MIXING_RATIO, /* 12 grams water vapor / kilogram dry air */
	_ENHANCEMENT_FACTOR, /* 13 unitless enhancement factor (non-ideal behavior of moist air) */
	_RELATIVE_HUMIDITY, /* 14 percent */
	_VAPOR_PRESSURE, /* 15 millibar */
	_POTENTIAL_VAPOR_PRESSURE, /* 16 millibar */
	_MOLE_MIXING_RATIO, /* 17 mole water vapor / mole moist air */
	_MASS_MIXING_RATIO, /* 18 grams water vapor / kilogram dry air */
	_DEW_POINT_C, /* 19 Degrees C */
	_DEW_POINT_K, /* 20 Kelvin */
	_DEW_POINT_F, /* 21 Degrees F */
	_SPECIFIC_HUMIDITY, /* 22 grams water vapor / kilogram moist air */
	_ABSOLUTE_HUMIDITY, /* 23 grams water vapor / meter^3 */
	_MOIST_AIR_DENSITY, /* 24 grams / meter^3 */
	_MOIST_AIR_NUMBER_DENSITY, /* 25 mole / meter^3 */
	_WATER_VAPOR_NUMBER_DENSITY, /* 26 mole / meter^3 */
	_GEOPOTENTIAL_HEIGHT, /* 27 in meters */
	_HEIGHT_AGL, /* 28 in meters */
	_HEIGHT_AMSL, /* 29 in meters */
	_DRY_AIR_NUMBER_DENSITY, /* 30 moles / meter^3 */
	_DRY_AIR_DENSITY, /* 31 grams / meter^3 */
	_SPEED_OF_SOUND, /* 32 meters / second */
	_MOIST_AIR_MOLAR_MASS, /* 33 g / mol */
	_OTHER_INPUT /* 34 Anything ... */};
#define N_FULL (uint)sizeof( WC_field_list_full)/sizeof(WEATHER_CONVERTER_FIELD)

int main(){
	WEATHER_CONVERSION_VECTOR STANDARD;
	WEATHER_CONVERSION_VECTOR TEST;
	WEATHER_CONVERTER_FIELD ri;
	BOOLEAN record_conversion[_N_WEATHER_FIELDS] = {0};

	printf("Running tests to see how well the conversions work.\n");

	printf("Performing a self-consistency test. (Can we go forward and back).\n");
	setUpTest(&STANDARD);
	printf("Test is now set up. Opening the conversion vector.\n");
	
	openWeatherConversionVector(&TEST, STANDARD.N);
	printf("Vector open was successful.\n");
	setQuiet(&STANDARD);
	setNotQuiet(&STANDARD);
	for(ri=_RELATIVE_HUMIDITY;ri<_MOIST_AIR_DENSITY;ri++){
		printf("Now performing the conversion starting with temperature, pressure, and %s.\n",_weather_converter_field_names[ri]);
		setTestVector(&TEST, &STANDARD, ri);
		printf("MEAN SQUARED RELATIVE ERROR = %g.\n",compare(&TEST, &STANDARD, record_conversion));
		recordErrorProneConversions(&TEST, &STANDARD, record_conversion, ri, "self-consistency");
	}
	freeWeatherConversionVector(&STANDARD);
	freeWeatherConversionVector(&TEST);

	printf("Performing a comparison to standard data.\n");
	if(importFile("humidityTest.csv", &STANDARD)){
		openWeatherConversionVector(&TEST, STANDARD.N);
		for(ri=_RELATIVE_HUMIDITY;ri<_MOIST_AIR_DENSITY;ri++){
			printf("Now performing the conversion starting with temperature, pressure, and %s.\n",_weather_converter_field_names[ri]);
			setTestVector(&TEST,&STANDARD,ri);
			printf("MEAN SQUARED RELATIVE ERROR = %g.\n",compare(&TEST,&STANDARD, record_conversion));
			recordErrorProneConversions(&TEST, &STANDARD, record_conversion, ri, "standard-data");
		}
		freeWeatherConversionVector(&STANDARD);
		freeWeatherConversionVector(&TEST);
	}
	else{
		printf("Error opening the standard data file 'humidityTest.csv'\n");
		printf("Skipping the standard comparison.\n");
	}
	return 0;
}

void setUpTest(WEATHER_CONVERSION_VECTOR *V){
	double P[] = {1.,5.,10.,15.,20.,40.,80.,100.,200.,300.,400.,500.,600.,700.,800.,900.,950.,1000.,1050.,1100.,1150.,1200.,1250.,1300.}; /* in mb */
	double T=250.,_T[] = {250.,5.,320.}; /* Temperature range (min,step,max) */
	double RH=0.,_RH[] = {0.0,4.99,100.}; /* Relative Humidity Range */

	size_t iP,iT,iRH,idx,NP=sizeof(P)/sizeof(double);
	size_t NT = (size_t) ( (_T[2]-_T[0])/_T[1]);
	size_t NRH = (size_t) ( (_RH[2]-_RH[0])/_RH[1]);

	printf("Opening the vector.");
	openWeatherConversionVector(V, (unsigned int)NP*NT*NRH);
	printf(" Done.");
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

	saveToFile(V, "temporaryHumidityTest.csv", WC_field_list_part, N_PART);
	saveToFile(V, "fullTestVector.csv", WC_field_list_full, N_FULL);
	saveSpeedOfSound(&(V->speedOfSoundDispersion), "speedOfSoundDispersion.csv");
}

void saveSpeedOfSound(M_BY_N_ARRAY *A, const char *fname){
	FILE *fp = fopen(fname, "w");
	unsigned int i, j;
	if(fp==NULL){
		printf("Unable to open %s for writing.\n",fname);
		return;
	}

	/* Header */
	fprintf(fp,"Index");
	for(i=0;i<A->M;i++){
		fprintf(fp,",%18.16lf",A->X[i]);
	}
	fprintf(fp,"\n");

	for(i=0;i<A->N;i++){
		fprintf(fp,"%u",i);
		for(j=0;j<A->M;j++){
			fprintf(fp,",%18.16lf",A->val[j][i]);
		}
		fprintf(fp,"\n");
	}
	fclose(fp);

}
void saveToFile(WEATHER_CONVERSION_VECTOR *V, const char *fname,
				const WEATHER_CONVERTER_FIELD *WC_field_list,
				uint N){
	FILE *fp;
	fp = fopen(fname,"w");
	uint i,j;
	WEATHER_CONVERTER_FIELD fi;

	if(fp==NULL){
		printf("Unable to open %s for writing.\n",fname);
		return;
	}

	/* Create a header */
	fi = WC_field_list[0];
	fprintf(fp, "%s (%s)", _weather_converter_field_names[fi], _weather_converter_field_units[fi]);
	for(j=1;j<N;j++){
		fi = WC_field_list[j];
		fprintf(fp, ",%s (%s)", _weather_converter_field_names[fi], _weather_converter_field_units[fi]);
	}
	fprintf(fp,"\n");

	for(i=0;i<V->N;i++){
		fi = WC_field_list[0];
		if(V->populated[fi]==FALSE) continue;
		fprintf(fp,"%18.16lf",V->val[fi][i]);
		for(j=1;j<N;j++){
			fi = WC_field_list[j];
			if(V->populated[fi]==FALSE) continue;
			fprintf(fp,",%18.16lf",V->val[fi][i]); /* <- added a comma */
		}
		fprintf(fp,"\n");
	}
	fclose(fp);
}

int importFile(const char *fname, WEATHER_CONVERSION_VECTOR *OUT){
	FILE *fp;
	fp = fopen(fname,"r");
	if (fp==NULL) return 0;
	uint N, i;
	char line[2000],*dmp;
	uint j,NF=(uint)sizeof( WC_field_list_part)/sizeof(WEATHER_CONVERTER_FIELD);
	WEATHER_CONVERTER_FIELD fi;

	for(N=0;feof(fp)==0;N++)
		if(fgets(line,2000,fp)==NULL){
			printf("converterTest.c:importFile() Successfully read %lu lines in %s.\n",N, fname);
			break;
		}

	rewind(fp);

	openWeatherConversionVector(OUT, N);
	for(i=0;feof(fp)==0;i++){
		if (i==N) break;
		if(fgets(line,2000,fp)==NULL)
			printf("converterTest.c:importFile() Error reading line %lu in %s.\n",i,fname);
		dmp = line;
		for(j=0;j<NF;j++){/* Iterate over all the standard input fields. */
			fi = WC_field_list_part[j];
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
		fi = WC_field_list_part[j];
		OUT->populated[fi]=TRUE;
	}
	OUT->populated[_U_WIND] = OUT->populated[_V_WIND]=TRUE;

	setAllFields(OUT);
	/* Fill in the remaining fields */
	fclose(fp);
	return 1;
}

double compare(WEATHER_CONVERSION_VECTOR *TST, WEATHER_CONVERSION_VECTOR *STD,
			   BOOLEAN *record_conversion){
	double totalErr = 0.0,tmpError;
	WEATHER_CONVERTER_FIELD ri;

	for(ri=0;ri<_N_WEATHER_FIELDS;ri++){
		if(ri==_OTHER_INPUT)continue;
		if((TST->populated[ri]==FALSE) || (STD->populated[ri]==FALSE))continue;
		if((tmpError=relError(TST->val[ri],STD->val[ri],TST->N)) > ((double)TST->N)*1.0e-6){
			printf("A high relative error, %lf, was seen for field # %d, %s.\n",tmpError/((double)TST->N),ri,_weather_converter_field_names[ri]);	
			record_conversion[ri] = TRUE;
		}
		else
			record_conversion[ri] = FALSE;
		totalErr += tmpError;
	}

	return totalErr/(double)_N_WEATHER_FIELDS/(double)TST->N;
}

void recordErrorProneConversions(WEATHER_CONVERSION_VECTOR *TEST, 
	WEATHER_CONVERSION_VECTOR *STANDARD, BOOLEAN *record_conversion, 
	const WEATHER_CONVERTER_FIELD si, const char *outfile_base){
	/* If there were any fields that did not compare well, then re-run the conversion,
but compare the processes. */
WEATHER_CONVERTER_FIELD i;
for(i=0; i<_N_WEATHER_FIELDS; i++){
	if(record_conversion[i]){
		recordSetTestVector(TEST, STANDARD, record_conversion, si, outfile_base);
		return;
	}
}
return;
}

void addField(FILE *fp, WEATHER_CONVERSION_VECTOR *WX, WEATHER_CONVERTER_FIELD fi, const char *prefix){
	fprintf(fp, "%s\"%s\":[",prefix, _weather_converter_field_flags[fi]);
	size_t i;
	if(WX->N > 0){
		fprintf(fp,"%g", WX->val[fi][0]);
		for(i=1;i<WX->N;i++)
			fprintf(fp,",%g",WX->val[fi][i]);
	}	
		fprintf(fp, "]");
}

void recordSetTestVector(WEATHER_CONVERSION_VECTOR *TST, 
	WEATHER_CONVERSION_VECTOR *STD, BOOLEAN *record_conversion, 
	const WEATHER_CONVERTER_FIELD field, const char *outfile_base){
	WEATHER_CONVERTER_FIELD ri;
	FILE *fp;
	char fname[1000 + 1];

	snprintf(fname, 1000, "%s_from_%s_TEST.csv", outfile_base, _weather_converter_field_flags[field]);
	saveToFile(TST, fname, WC_field_list_part, N_PART);

	snprintf(fname, 1000, "%s_from_%s.json", outfile_base, _weather_converter_field_flags[field]);

	if ((fp = fopen(fname, "w")) == NULL){
		printf("Unable to open %s for writing.", fname);
		return;
	}
	
	fprintf(fp, "{\n");

    fprintf(fp, " \"standardVariables\": {\n");
		addField(fp, STD, 0, "  ");
		for(ri=1;ri<_N_WEATHER_FIELDS;ri++){
			if(STD->populated[ri])
				addField(fp, STD, ri, ",\n  ");
		}

	fprintf(fp, " },\n");

	fprintf(fp, " \"setVariables\": {\n");
		addField(fp, STD, _U_WIND, "  ");		
		addField(fp, STD, _V_WIND, ",\n  ");		
		addField(fp, STD, _TEMPERATURE_K, ",\n  ");		
		addField(fp, STD, _PRESSURE, ",\n  ");		
		addField(fp, STD, field, ",\n  ");		
	fprintf(fp, " },\n");

	fprintf(fp, " \"calculatedVariables\": {\n");
		for(ri=0;ri<_N_WEATHER_FIELDS;ri++){
			if(record_conversion[ri] && STD->populated[ri]){
				addField(fp, TST, ri, "  ");
				break;
			}
		}
		for(;ri<_N_WEATHER_FIELDS;ri++){
			if(record_conversion[ri] && STD->populated[ri]){
				addField(fp, TST, ri, ",\n  ");
			}
		}

	fprintf(fp, " }\n");

	fprintf(fp, "}");

	fclose(fp);
		
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
