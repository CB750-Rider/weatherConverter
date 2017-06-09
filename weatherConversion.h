/*
 * weatherConversion.h
 *
 *  Created on: Aug 11, 2016
 *      Author: Dr. Lee Burchett 
 *
 *      C code for converting/populating weather data.
 */

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
	_DEW_POINT, /* 19 Kelvin */
	_SPECIFIC_HUMIDITY, /* 20 grams water vapor / kilogram moist air */
	_ABSOLUTE_HUMIDITY, /* 21 grams water vapor / meter^3 */
	_MOIST_AIR_DENSITY, /* 22 grams / meter^3 */
	_OTHER_INPUT, /* Anything else */
	_N_WEATHER_FIELDS
} WEATHER_CONVERTER_FIELD;

typedef enum{
	WEATHER_CONVERSION_SUCCESS,
	NO_TEMPERATURE_PRESENT,
	NO_PRESSURE_PRESENT,
	NO_HUMIDITY_PRESENT,
	NO_WIND_PRESENT,
	FIELD_NOT_ALLOCATED,
	FILE_IO_ERROR,
	DOUBLE_INITIALIZATION_ERROR,
	ALLOCATION_ERROR,
	UNKNOWN_ERROR,
	N_WEATHER_CONVERSION_ERRORS
} WEATHER_CONVERSION_ERROR;

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
	WEATHER_CONVERSION_ERROR (*set_field) (void *WEATHER_CONVERSION_VECTOR,
										   double *value,
                                           WEATHER_CONVERTER_FIELD field);
    WEATHER_CONVERSION_ERROR (*run_conversion) (void *WEATHER_CONVERSION_VECTOR);                                       
} WEATHER_CONVERSION_FUNCTIONS;

typedef struct{
	unsigned int N;
	double *val[_N_WEATHER_FIELDS];
	unsigned int allocated[_N_WEATHER_FIELDS];
	unsigned int populated[_N_WEATHER_FIELDS];
	double standardPressure; /* in millibars */
	double xCO2; /* in ppm*/
    WEATHER_CONVERSION_FUNCTIONS f; /* Reference available functions */
} WEATHER_CONVERSION_VECTOR;

/* Strings defined in weatherConversion.c */
char *_weather_converter_field_names[_N_WEATHER_FIELDS];
char *_weather_converter_field_flags[_N_WEATHER_FIELDS];
char *_weather_converter_field_units_full[_N_WEATHER_FIELDS];
char *_weather_converter_field_units[_N_WEATHER_FIELDS];
/* Functions defined in weatherConversion.c */
double m_sToKnots(double ms);
double knotsTom_s(double kt);
double KtoF(double K);
double FtoK(double F);
double CtoF(double Cent);
double FtoC(double F);
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
double calcMixingRatio(double P, double vp);
double calcMMRfromAbsoluteHumidity(double P, double T, double a);
double calcSHfromRH(double T, double P, double xCO2, double SVP, double RH, double ef);
double estRHfromSH(double SH, double T, double P, double xCO2, double SVP, double ef, double RH);
WEATHER_CONVERTER_FIELD presentHumidity(WEATHER_CONVERSION_VECTOR *WX);
WEATHER_CONVERSION_ERROR setAllFields(WEATHER_CONVERSION_VECTOR *WX);
WEATHER_CONVERSION_ERROR setTemperatures(WEATHER_CONVERSION_VECTOR *WX);
WEATHER_CONVERSION_ERROR setPressures(WEATHER_CONVERSION_VECTOR *WX);
WEATHER_CONVERSION_ERROR setWinds(WEATHER_CONVERSION_VECTOR *WX);
WEATHER_CONVERSION_ERROR humidityConversion(WEATHER_CONVERSION_VECTOR *WX);
WEATHER_CONVERSION_ERROR altitudeCalcualtion(WEATHER_CONVERSION_VECTOR *WX);
/* Variables required for error reporting. Definitions are in conversionVectorHelpers.c */
char *_weather_converter_error_messages[N_WEATHER_CONVERSION_ERRORS];
WEATHER_CONVERSION_ERROR _weather_converter_global_return;
FILE *weather_converter_error_stream;
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

#endif /* WEATHERCONVERSION_H_ */
