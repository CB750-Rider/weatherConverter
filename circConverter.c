/*
 * converterTest.c
 *
 *  Created on: Aug 11, 2016
 *      Author: lee
 */

#include "weatherConversion.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

int importFile(const char *fname, WEATHER_CONVERSION_VECTOR *OUT);
void exportData(const char *fname, WEATHER_CONVERSION_VECTOR *OUT);
void fileConvert(const char *fname);

int main(int argc, char *argv[]){

	int i;

	if (argc < 2){
		printf("Please provide a file name to open and process.\n");
		return 0;
	}

	for(i=1;i<argc;i++){
		printf("Now converting %s.\n",argv[i]);
		fileConvert(argv[i]);
	}

	printf("Conversions Complete.\n");

	return 0;
}

void fileConvert(const char *fname){
	WEATHER_CONVERSION_VECTOR work;
	char outFileName[500];

	snprintf(outFileName,500,"%s_converted.csv",fname);

	if(importFile(fname, &work) == 0){
		setAllFields(&work); /* performs conversion */
		exportData(outFileName,&work);
		freeWeatherConversionVector(&work);
	}
}

int importFile(const char *fname, WEATHER_CONVERSION_VECTOR *OUT){
	uint N,i,dmpi;
	long data_start;
	char line[2000],dmp[2000];
	FILE *fp= fopen(fname,"r");

	if(fp==NULL){
		printf("Error opening %s for reading.\n",fname);
		return -1;
	}

	/* The first three lines are boiler-plate */
	for(N=0;N<3;N++){
		if(fgets(line,2000,fp)==NULL){
			printf("circConverter.c:importFile() Error reading line %d of the header in in %s.\n",N,fname);
			break;
		}
	}
	data_start = ftell(fp);
	/* Count the rest of the lines */
	for(N=0;feof(fp)==0;N++)
		if(fgets(line,2000,fp)==NULL){
			printf("circConverter.c:importFile() Successfully read %d lines in %s.\n",N,fname);
			break;
		}

	fseek(fp,data_start,SEEK_SET);

	openWeatherConversionVector(OUT, N);

	for(i=0;feof(fp)==0;i++){
		if (i>=N) break;
		if(fgets(line,2000,fp)==NULL)
			printf("converterTest.c:importFile() Error reading line %d in %s.\n",i,fname);
		if(sscanf(line,"%u %lf %lf %lg %s",&dmpi,
										  OUT->val[_PRESSURE]+i,
										  OUT->val[_TEMPERATURE_K]+i,
										  OUT->val[_MOLE_MIXING_RATIO]+i,
										  dmp)!=5) printf("Error converting values in line %d.\n",i);
		OUT->val[_V_WIND][i] = 1.0;
		OUT->val[_U_WIND][i] = 1.0;
	}
	OUT->populated[_TEMPERATURE_K]=TRUE;
	OUT->populated[_PRESSURE]=TRUE;
	OUT->populated[_RELATIVE_HUMIDITY]=FALSE;
	OUT->populated[_ABSOLUTE_HUMIDITY]=FALSE;
	OUT->populated[_DEW_POINT_C]=FALSE;
	OUT->populated[_DEW_POINT_F]=FALSE;
	OUT->populated[_DEW_POINT_K]=FALSE;
	OUT->populated[_ENHANCEMENT_FACTOR]=FALSE;
	OUT->populated[_MASS_MIXING_RATIO]=FALSE;
	OUT->populated[_MOIST_AIR_DENSITY]=FALSE;
	OUT->populated[_MOLE_MIXING_RATIO]=TRUE;
	OUT->populated[_POTENTIAL_TEMPERATURE]=FALSE;
	OUT->populated[_SATURATION_MIXING_RATIO]=FALSE;
	OUT->populated[_SATURATION_VAPOR_PRESSURE]=FALSE;
	OUT->populated[_SPECIFIC_HUMIDITY]=FALSE;
	OUT->populated[_VAPOR_PRESSURE]=FALSE;
	OUT->populated[_VIRTUAL_POTENTIAL_TEMPERATURE]=FALSE;
	OUT->populated[_VIRTUAL_TEMPERATURE]=FALSE;
	OUT->populated[_U_WIND]=TRUE;
	OUT->populated[_V_WIND]=TRUE;

	fclose(fp);
	return 0;
}

void exportData(const char *fname, WEATHER_CONVERSION_VECTOR *OUT){
	uint i;
	FILE *fp= fopen(fname,"w");

	if(fp==NULL){
		printf("Error opening %s for writing output.\n",fname);
		return;
	}
	fprintf(fp,"Level Num,P(MB),T(K),RH,Mol Mixing Ratio\n");
	for(i=0;i<OUT->N;i++){
		fprintf(fp,"%u,%lf,%lf,%lg,%lg\n",i+1,
									 OUT->val[_TEMPERATURE_K][i],
									 OUT->val[_PRESSURE][i],
									 OUT->val[_RELATIVE_HUMIDITY][i],
									 OUT->val[_MOLE_MIXING_RATIO][i]);
	}

	fclose(fp);
}

