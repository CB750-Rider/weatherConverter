#include "weatherConversion.h"

/*
Store all the enum-dependent strings we use.
*/

const char *_weather_converter_field_names[_N_WEATHER_FIELDS] = {
													"Temperature", /*0*/
													"Temperature", /*1*/
													"Temperature", /*2*/
													"U wind", /*3*/
													"V wind", /*4*/
													"Wind Speed", /*5*/
													"Wind Direction", /*6*/
													"Pressure", /*7*/
													"Potential Temperature", /*8*/
													"Virtual Temperature", /*9*/
													"Virtual Potential Temperature", /*10*/
													"Saturation Vapor Pressure", /*11*/
													"Saturation Mass Mixing Ratio", /*12*/
													"Enhancement Factor", /*13*/
													"Relative Humidity", /*14*/
													"Vapor Pressure", /*15*/
													"Potential Vapor Pressure", /*16*/
													"Mole Mixing Ratio", /*17*/
													"Mass Mixing Ratio", /*18*/
													"Dew Point C", /*19*/
													"Dew Point K", /*20*/
													"Dew Point F", /*21*/
													"Specific Humidity", /*22*/
													"Absolute Humidity", /*23*/
													"Moist Air Density", /*24*/
													"Moist Air Number Density", /*25*/
													"Water Vapor Number Density", /*26*/
													"Geopotential Height", /*27*/
													"Height Above Ground Level", /*28*/
													"Height Above Mean Sea Level", /*29*/
													"Dry Air Number Density", /*30*/
													"Dry Air Density", /*31*/
													"Speed of Sound (dispersionless)", /*32*/
													"Moist Air Molar Mass", /*33*/
													"Other Input Field"/*34*/};

const char *_weather_converter_field_flags[_N_WEATHER_FIELDS] = {
	"_TEMPERATURE_C", /* 0 Degrees C */
	"_TEMPERATURE_K", /* 1 Kelvin */
	"_TEMPERATURE_F", /* 2 Degrees F */
	"_U_WIND", /* 3 meters / second*/
	"_V_WIND",  /* 4 meters / second*/
	"_WIND_SPEED",  /* 5 meters / second*/
	"_WIND_DIRECTION", /* 6 degrees from N */
	"_PRESSURE", /* 7 millibar */
	"_POTENTIAL_TEMPERATURE", /* 8 Kelvin */
	"_VIRTUAL_TEMPERATURE", /* 9 Kelvin */
	"_VIRTUAL_POTENTIAL_TEMPERATURE", /* 10 Kelvin */
	"_SATURATION_VAPOR_PRESSURE", /* 11 millibar */
	"_SATURATION_MIXING_RATIO", /* 12 grams water vapor / kilogram dry air */
	"_ENHANCEMENT_FACTOR", /* 13 unitless enhancement factor (non-ideal behavior of moist air) */
	"_RELATIVE_HUMIDITY", /* 14 percent */
	"_VAPOR_PRESSURE", /* 15 millibar */
	"_POTENTIAL_VAPOR_PRESSURE", /* 16 millibar */
	"_MOLE_MIXING_RATIO", /* 17 mole water vapor / mole moist air */
	"_MASS_MIXING_RATIO", /* 18 grams water vapor / kilogram dry air */
	"_DEW_POINT_C", /* 19 Kelvin */
	"_DEW_POINT_K", /* 20 Kelvin */
	"_DEW_POINT_F", /* 21 Kelvin */
	"_SPECIFIC_HUMIDITY", /* 22 grams water vapor / kilogram moist air */
	"_ABSOLUTE_HUMIDITY", /* 23 grams water vapor / meter^3 */
	"_MOIST_AIR_DENSITY", /* 24 grams / meter^3 */
	"_MOIST_AIR_NUMBER_DENSITY", /* 25 mole / meter^3 */
	"_WATER_VAPOR_NUMBER_DENSITY", /* 26 mole / meter^3 */
	"_GEOPOTENTIAL_HEIGHT", /* 27 meters */
	"_HEIGHT_AGL", /* 28 meters */
	"_HEIGHT_AMSL", /* 29 meters*/
	"_DRY_AIR_NUMBER_DENSITY", /* 30 moles / meter^3 */
	"_DRY_AIR_DENSITY", /* 31 grams / meter^3 */
	"_SPEED_OF_SOUND", /* 32 meters / second */
	"_MOIST_AIR_MOLAR_MASS", /* 33 grams / mole */
    "_OTHER_INPUT" /* 34 For any other field */};

const char *_weather_converter_field_units_full[_N_WEATHER_FIELDS] = {
													"degrees Celsius", /*0*/
													"Kelvin", /*1*/
													"degrees Fahrenheit", /*2*/
													"meters / second", /*3*/
													"meters / second", /*4*/
													"meters / second", /*5*/
													"degrees from North", /*6*/
													"millibar", /*7*/
													"Kelvin", /*8*/
													"Kelvin", /*9*/
													"Kelvin", /*10*/
													"millibar", /*11*/
													"grams water vapor / kilogram moist air", /*12*/
													"unitless", /*13*/
													"percent", /*14*/
													"millibar", /*15*/
													"millibar", /*16*/
													"moles water vapor / mole moist air", /*17*/
													"grams water vapor / kilogram moist air", /*18*/
													"degrees Celsius", /*19*/
													"Kelvin", /*20*/
													"degrees Fahrenheit", /*21*/
													"grams water vapor / kilogram moist air", /*22*/
													"grams water vapor / meter^3", /*23*/
													"grams / meter^3", /*24*/
													"mole / meter^3",/*25*/
													"mole / meter^3",/*26*/
													"meters", /*27*/
													"meters", /*28*/
													"meters", /*29*/
													"moles / meter^3", /*30*/
													"grams / meter^3", /*31*/
													"meters / second", /*32*/
													"grams / mole", /*33*/
													"undefined"/*34*/};

const char *_weather_converter_field_units[_N_WEATHER_FIELDS] = {
													"°C", /*0*/
													"K", /*1*/
													"°F", /*2*/
													"m/s", /*3*/
													"m/s", /*4*/
													"m/s", /*5*/
													"°", /*6*/
													"mb", /*7*/
													"K", /*8*/
													"K", /*9*/
													"K", /*10*/
													"mb", /*11*/
													"g/kg", /*12*/
													"<unitless>", /*13*/
													"%%", /*14*/
													"mb", /*15*/
													"mb", /*16*/
													"mol/mol", /*17*/
													"g/kg", /*18*/
													"°C", /*19*/
													"K", /*20*/
													"°F", /*21*/
													"g/kg", /*22*/
													"g/m^3", /*23*/
													"g/m^3", /*24*/
													"mol/m^3",/*25*/
													"mol/m^3",/*26*/
													"m", /*27*/
													"m", /*28*/
													"m", /*29*/
													"mol/m^3", /*30*/
													"g/m^3", /*31*/
													"m/s", /*32*/
													"g/mol", /*33*/
													"----"/*34*/};

const char *_weather_converter_site_setting_names[_N_WEATHER_SITE_SPECIFIC_SETTINGS]  = {
		"Standard Pressure",
		"CO2 Parts Per Million",
		"Latitude",
		"Surface Height AMSL",
		"Surface Pressure",
		"Surface Temperature",
		"CH4 Parts Per Billion",
		"N2O Parts Per Billion",
		"SF6 Parts Per Billion"
};
const char *_weather_converter_site_setting_flags[_N_WEATHER_SITE_SPECIFIC_SETTINGS] = {
		"_STANDARD_PRESSURE", /* Standard pressure to use (in mb) default is 1000 */
		"_XCO2", /* CO2 mole mixing ratio in ppm default is 390.0*/
		"_SITE_LATITUDE",  /* In degrees, default is 35 */
		"_SURFACE_HEIGHT", /* Elevation in meters above mean sea level, default is 0 */
		"_SURFACE_PRESSURE", /* Surface pressure at 0 m AMSL in mb, default is 1013.25 */
		"_SURFACE_TEMPERATURE", /* Surface pressure at 0 m AMSL in mb, default is 1013.25 */
		"_XCH4",
		"_XN2O",
		"_XSF6"
};
const char *_weather_converter_site_units_full[_N_WEATHER_SITE_SPECIFIC_SETTINGS] = {
		"millibars",
		"parts per million",
		"degrees latitutude",
		"meters above mean sea level",
		"millibars",
		"Kelvin",
		"parts per billion",
		"parts per billion",
		"parts per trillion"
};
const char *_weather_converter_site_units[_N_WEATHER_SITE_SPECIFIC_SETTINGS] = {
		"mb",
		"ppm",
		"deg",
		"m",
		"mb",
		"K",
		"ppb",
		"ppb",
		"ppt"
};
