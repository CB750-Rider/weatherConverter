/* weatherConversion.h
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
 * Thes file contains header prototypes and tools for routines that are useful 
 * for converting the myriad ways that researchers use to represent humidity in 
 * the atmosphere. The goal is to provide a tool to save time on writing and 
 * debugging these conversions. If any conversion data appear suspect. Please 
 * contact me with suggested changes at lee.r.burchett@gmail.com.
 *

 */

#ifdef _MSC_BUILD
#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif
#endif

#ifndef WEATHERCONVERSION_H_
#define WEATHERCONVERSION_H_

#include <math.h>
#include <stdio.h>

typedef enum{
	BACKWARD,
	FOREWARD
} CALCULATION_DIRECTION;

#ifndef BOOLEAN_DEFINED
#define BOOLEAN_DEFINED 1
typedef enum{
	FALSE,
	TRUE
} BOOLEAN;
#endif

#define WATER_MOLAR_MASS 18.015261322024042

/* Check for successful function completion. If there is a problem, then send an
 * error to the user. */
#define _wcCheck(fun) do{ \
	if((_weather_converter_global_return = fun) != WEATHER_CONVERSION_SUCCESS){ \
		reportWeatherConversionError("\
***ERROR IN %s at %d***\n\
  Error reported when calling:\n\
  " #fun ".\n\
  The reported error was %s.\n\
***END ERROR***\n", \
		__FILE__,__LINE__, _weather_converter_error_messages[_weather_converter_global_return]); \
		return _weather_converter_global_return; \
}}while(0)

/* TODO Add potential temp, virt temp, virt potnl temp in F and C,
Add Wind speed in knots 
Add other pressure units 
*/
typedef enum {
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
	_SPEED_OF_SOUND, /* 30 in meters / second */
	_OTHER_INPUT, /* 31 Anything else. This should be the last enum before _N_WEATHER_FIELDS. */
	_N_WEATHER_FIELDS
} WEATHER_CONVERTER_FIELD;

typedef enum{
	WEATHER_CONVERSION_SUCCESS,
	NO_TEMPERATURE_PRESENT,
	NO_PRESSURE_PRESENT,
	NO_HUMIDITY_PRESENT,
	NO_WIND_PRESENT,
	FIELD_NOT_ALLOCATED,
	WC_FILE_IO_ERROR,
	DOUBLE_INITIALIZATION_ERROR,
	ALLOCATION_ERROR,
	UNKNOWN_ERROR,
	INITIALIZATION_ERROR,
	N_WEATHER_CONVERSION_ERRORS
} WEATHER_CONVERSION_ERROR;

typedef enum{
	/* Verify default values in WeatherConversion.c*/
	_STANDARD_PRESSURE, /* Standard pressure to use (in mb) default is 1000 */
	_XCO2, /* CO2 mole mixing ratio in ppm default is 390.0*/
	_SITE_LATITUDE,  /* In degrees, default is 35 */
	_SURFACE_HEIGHT, /* Elevation in meters above mean sea level, default is 0 */
	_SURFACE_PRESSURE, /* Surface pressure at 0 m AMSL in mb, default is 1013.25 */
	_SURFACE_TEMPERATURE, /* Surface temperature, default is 273.15 K (0° C)*/
	_N_WEATHER_SITE_SPECIFIC_SETTINGS
} SITE_SPECIFIC_SETTINGS;

/*
typedef struct{
	WEATHER_CONVERSION_ERROR (*alloc)(void *WEATHER_CONVERSION_VECTOR, unsigned int size);
	WEATHER_CONVERSION_ERROR (*free) (void *WEATHER_CONVERSION_VECTOR);
	WEATHER_CONVERSION_ERROR (*clear_all) (void *WEATHER_CONVERSION_VECTOR);
    WEATHER_CONVERSION_ERROR (*change_temp) (void *WEATHER_CONVERSION_VECTOR,
                                             double *new_temps,
											 WEATHER_CONVERTER_FIELD temp_type);
    WEATHER_CONVERSION_ERROR (*change_pressure) (void *WEATHER_CONVERSION_VECTOR,
                                                 double *new_pressure,
												 WEATHER_CONVERTER_FIELD pressure_type);
    WEATHER_CONVERSION_ERROR (*change_humidity) (void *WEATHER_CONVERSION_VECTOR,
                                             double *new_humidity,
                                             WEATHER_CONVERTER_FIELD humidity_type);
    WEATHER_CONVERSION_ERROR (*change_height) (void *WEATHER_CONVERSION_VECTOR,
                                             double *new_height,
                                             WEATHER_CONVERTER_FIELD height_type);
	WEATHER_CONVERSION_ERROR (*set_field) (void *WEATHER_CONVERSION_VECTOR,
										   double *value,
                                           WEATHER_CONVERTER_FIELD field);
    WEATHER_CONVERSION_ERROR (*run_conversion) (void *WEATHER_CONVERSION_VECTOR);
    double (*integrate_column_water_density)(void *WEATHER_CONVERSION_VECTOR,double z0, double z1);
	double(*integrate_column_water_number_density)(void *WEATHER_CONVERSION_VECTOR, double z0, double z1);
	double(*integrate_column_moist_air_density)(void *WEATHER_CONVERSION_VECTOR, double z0, double z1);
	double(*integrate_column_moist_air_number_density)(void *WEATHER_CONVERSION_VECTOR, double z0, double z1);
} WEATHER_CONVERSION_FUNCTIONS;
 */

typedef struct{
	unsigned int N;
	double *val[_N_WEATHER_FIELDS];
	unsigned int allocated[_N_WEATHER_FIELDS];
	unsigned int populated[_N_WEATHER_FIELDS];
	double standardPressure; /* in millibars */
	double xCO2; /* in ppm */
	double vaporArealDensity; /* For the column in g / m^2 */
	double vaporArealNumberDensity; /* For the column in mol / m^2 */
	double moistAirArealDensity; /* For the column in g / m^2 */
	double moistAirArealNumberDensity; /* For the column in mol / m^2 */
	double latitude; /* In degrees */
	double surfaceHeight; /* Above mean sea level, in meters */
	double surfacePressure; /* in mb */
	double surfaceTemperature; /* in Kelvin */
	int quiet; /* Verbosity flag */
} WEATHER_CONVERSION_VECTOR;

/* Strings defined in weatherConversion.c */
extern char *_weather_converter_field_names[_N_WEATHER_FIELDS];
extern char *_weather_converter_field_flags[_N_WEATHER_FIELDS];
extern char *_weather_converter_field_units_full[_N_WEATHER_FIELDS];
extern char *_weather_converter_field_units[_N_WEATHER_FIELDS];
extern char *_weather_converter_site_setting_names[_N_WEATHER_SITE_SPECIFIC_SETTINGS] ;
extern char *_weather_converter_site_setting_flags[_N_WEATHER_SITE_SPECIFIC_SETTINGS];
extern char *_weather_converter_site_units_full[_N_WEATHER_SITE_SPECIFIC_SETTINGS] ;
extern char *_weather_converter_site_units[_N_WEATHER_SITE_SPECIFIC_SETTINGS];
extern double _weather_converter_site_defaults[_N_WEATHER_SITE_SPECIFIC_SETTINGS];
/* Functions defined in weatherConversion.c */
double latitude_gravity(double lat);
double latitude_earth_radius(double lat);
double free_air_gravity(double lat, double h);
double m_sToKnots(double ms);
double knotsTom_s(double kt);
double KtoF(double K);
double FtoK(double F);
double CtoF(double Cent);
double FtoC(double F);
double CtoK(double C);
double KtoC(double K);
double calcDewpoint(double T, double RH, double SVP);
double goffGratch(double T);
double DgoffGratch_dT(double T);
double goffGratchIce(double T);
double DgoffGratchIce_dT(double T);
double enhancementFactor(double T, double P);
double calcAbsolute(double P, double T, double moleMixingRatio);
double calcCompressibility(double P, double T, double xv);
double calcPotentialTemperature(double T, double P, double P0, CALCULATION_DIRECTION D);
double calcVirtualTemperature(double T, double mr, CALCULATION_DIRECTION D);
double calcMoistAirDensity(double T, double P, double xv, double xCO2);
double moistAirNumberDensity(double P, double T, double xv);
double dZ_dxv(double T, double p, double x);
double calcSpecificHumidity(double mr);
double calcMassMixingRatio(double P, double vp);
double calcMMRfromAbsoluteHumidity(double P, double T, double a);
double calcSHfromRH(double T, double P, double xCO2, double SVP, double RH, double ef);
double estRHfromSH(double SH, double T, double P, double xCO2, double SVP, double ef, double RH);
int findStartIndex(double *x, double x0,int N);
int findEndIndex(double *x,double x0,int N);
double logarithmicRule(double *x, double *y, int i);
double trapezoidalRule(double *x, double *y, int i);
double logarithmicPartial(double *x, double *y, int i, double a, double b);
double trapezoidalPartial(double *x, double *y, int i, double a, double b);
double integrateColumnWaterDensity(WEATHER_CONVERSION_VECTOR *WX,double z0, double z1);
double integrateColumnWaterNumberDensity(WEATHER_CONVERSION_VECTOR *WX, double z0, double z1);
double integrateColumnMoistAirDensity(WEATHER_CONVERSION_VECTOR *WX, double z0, double z1);
double integrateColumnMoistAirNumberDensity(WEATHER_CONVERSION_VECTOR *WX, double z0, double z1);
double standardAtmosAltitudeAtPressure(double P);
double standardAtmosPressureAtAltitude(double Z);
double speedOfSound(double T, double mass_mix_ratio, double P);
double _relError(double a, double b);

WEATHER_CONVERTER_FIELD presentHumidity(WEATHER_CONVERSION_VECTOR *WX);
WEATHER_CONVERSION_ERROR setAllFields(WEATHER_CONVERSION_VECTOR *WX);
WEATHER_CONVERSION_ERROR setTemperatures(WEATHER_CONVERSION_VECTOR *WX);
WEATHER_CONVERSION_ERROR setPressures(WEATHER_CONVERSION_VECTOR *WX);
WEATHER_CONVERSION_ERROR setHeights(WEATHER_CONVERSION_VECTOR *WX);
WEATHER_CONVERSION_ERROR setWinds(WEATHER_CONVERSION_VECTOR *WX);
WEATHER_CONVERSION_ERROR humidityConversion(WEATHER_CONVERSION_VECTOR *WX);
/* Variables required for error reporting. Definitions are in conversionVectorHelpers.c */
extern char *_weather_converter_error_messages[N_WEATHER_CONVERSION_ERRORS];
extern WEATHER_CONVERSION_ERROR _weather_converter_global_return;
extern FILE *weather_converter_error_stream;
/* Functions in conversionVectorHelpers.c */
void reportWeatherConversionError(char *msg, ...);
WEATHER_CONVERSION_ERROR openWeatherConversionVector(WEATHER_CONVERSION_VECTOR *WX, unsigned int N);
WEATHER_CONVERSION_ERROR freeWeatherConversionVector(WEATHER_CONVERSION_VECTOR *WX);
WEATHER_CONVERSION_ERROR setWeatherField(WEATHER_CONVERSION_VECTOR *WX, double *x, WEATHER_CONVERTER_FIELD fi);
WEATHER_CONVERSION_ERROR initializeVector(WEATHER_CONVERSION_VECTOR *WX);
WEATHER_CONVERSION_ERROR clearAllWeatherConversionFields(WEATHER_CONVERSION_VECTOR *WX);
WEATHER_CONVERSION_ERROR changeTemperatures(WEATHER_CONVERSION_VECTOR *WX, double *t, WEATHER_CONVERTER_FIELD fi);
WEATHER_CONVERSION_ERROR changePressures(WEATHER_CONVERSION_VECTOR *WX, double *p,  WEATHER_CONVERTER_FIELD fi);
WEATHER_CONVERSION_ERROR changeHumidity(WEATHER_CONVERSION_VECTOR *WX, double *h,  WEATHER_CONVERTER_FIELD fi);
WEATHER_CONVERSION_ERROR changeHeight(WEATHER_CONVERSION_VECTOR *WX, double *h,WEATHER_CONVERTER_FIELD fi);
WEATHER_CONVERSION_ERROR parseSiteSettingnsLine(char *line, WEATHER_CONVERSION_VECTOR *WX);

#endif // WEATHERCONVERSION_H_ 
