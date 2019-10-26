/* weatherConversion.c
 
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
 *      Author: Dr. Lee Burchett  
 *
 *  Last Modified: 20 August, 2019
 * 
 * Thes file contains routines useful for converting the myriad ways that 
 * researchers use to represent humidity in the atmosphere. The goal is to 
 * provide a tool to save time on writing and debugging these conversions. 
 * If any conversion data appear suspect. Please contact me with suggested
 * changes at lee.r.burchett@gmail.com.

 */

#include "weatherConversion.h"
#ifdef _MSC_BUILD
#define _USE_MATH_DEFINES
#include <math.h>
#define _MY_INFINITY HUGE_VAL
#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif
#else
#define _MY_INFINITY 1.0/0.0
#endif
#include <stdlib.h>
#include <stdio.h>

static int is_unstable(double num,double denom, double X0);

char *_weather_converter_field_names[_N_WEATHER_FIELDS] = {
													"Temperature C", /*0*/
													"Temperature K", /*1*/
													"Temperature F", /*2*/
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
													"Other Input Field"/*30*/};

char *_weather_converter_field_flags[_N_WEATHER_FIELDS] = {
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
    "_OTHER_INPUT" /* 30 For any other field */};

char *_weather_converter_field_units_full[_N_WEATHER_FIELDS] = {
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
													"undefined"/*30*/};

char *_weather_converter_field_units[_N_WEATHER_FIELDS] = {
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
													" ", /*13*/
													"%", /*14*/
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
													"----"/*30*/};

char *_weather_converter_site_setting_names[_N_WEATHER_SITE_SPECIFIC_SETTINGS]  = {
		"Standard Pressure",
		"CO2 Parts Per Million",
		"Latitude",
		"Surface Height AMSL",
		"Surface Pressure"
};
char *_weather_converter_site_setting_flags[_N_WEATHER_SITE_SPECIFIC_SETTINGS] = {
		"_STANDARD_PRESSURE", /* Standard pressure to use (in mb) default is 1000 */
		"_XCO2", /* CO2 mole mixing ratio in ppm default is 390.0*/
		"_SITE_LATITUDE",  /* In degrees, default is 35 */
		"_SURFACE_HEIGHT", /* Elevation in meters above mean sea level, default is 0 */
		"_SURFACE_PRESSURE"
};
char *_weather_converter_site_units_full[_N_WEATHER_SITE_SPECIFIC_SETTINGS] = {
		"millibars",
		"parts per million",
		"degrees latitutude",
		"meters above mean sea level",
		"millibars"
};
char *_weather_converter_site_units[_N_WEATHER_SITE_SPECIFIC_SETTINGS] = {
		"mb",
		"ppm",
		"deg",
		"m",
		"mb"
};
double _weather_converter_site_defaults[_N_WEATHER_SITE_SPECIFIC_SETTINGS] = {
		1000.0, /* Standard pressure */
		390.0, /* CO2 ppm */
		35.0, /* Site Latitutde */
		0.0, /* Site surface elevation AMSL*/
		1013.25 /* Surface pressure */
};
int is_unstable(double num,double denom,double X0){
	/* Check for numerical stability of division */
	double lnum,lden,lx;
	if(num==0.0 || denom==0.0)return 1;
	lnum = log10(fabs(num));
	lden = log10(fabs(denom));
	lx = log10(fabs(X0));
	/* If the ratio of lnum to lden is much greater than X0, then
	 * the results are probably ill-conditioned */
	if( (lnum - lden) > (lx - 2)) return 1;
	return (fabs(lnum-lden) > 12);
}
double latitude_gravity(double lat){
	/* Gravity formula 1980 from The Geodesist's Handbook 2000, p132
	 * input is latitude in degrees
	 * output is the acceleration of gravity at sea level at latitude  */
	double sterm = pow(sin(lat*M_PI/180.0),2.0);
	double gamma_e = 9.7803267715; /* Standard equatorial gravity */
	double k = 0.001931851353; /* See handbook p 131 */
	double e2 = 0.00669438002290; /* first Earth eccentricity squared */
	return gamma_e*(1 + k*sterm)/sqrt(1-e2*sterm);
}
double latitude_earth_radius(double lat){
	double a = 6378137.0; /* Earth semi-major axis */
	double b = 6356752.3; /* Earth semi-minor axis */
	double cl2  = pow(cos(lat*M_PI/180.0),1.0);
	double sl2  = pow(sin(lat*M_PI/180.0),1.0);
	return sqrt((a*a*a*a*cl2 + b*b*b*b*sl2)/(a*a*cl2 + b*b*sl2));
}
double free_air_gravity(double lat, double h){
	/* Using the universal law of gravitation and local radius and height
	 * inputs are latitude in degrees and height in meters */
	double R = latitude_earth_radius(lat);
	return latitude_gravity(lat)*pow(1.0/(1.0 + h/R),2.0);
}
double m_sToKnots(double ms){
	return ms*1.9438444924574;}
double knotsTom_s(double kt){
	return kt* 0.51444444444;}
double KtoF(double K){
	return 9.0/5.0*KtoC(K)+32;}
double FtoK(double F){
	return CtoK(5.0/9.0*(F-32.0));}
double CtoF(double C){
	return (9.0/5.0*C + 32.0);}
double FtoC(double F){
	return 5.0/9.0*(F-32.0);}
double CtoK(double C){
	return C + 273.15;
}
double KtoC(double K){
	return K - 273.15;
}
double calcDewpoint(double T, double RH, double SVP){
	/* Based on the method presented by NOAA,
	 * http://www.srh.noaa.gov/images/epz/wxcalc/wetBulbTdFromRh.pdf
	 *  finds the dewpoint in Kelvin using:
	 *  T	the temperature (K),
	 *  RH	relative humidity (%), and
	 *  SVP	saturation vapor pressure (mb)
	 *  Since this function only approximately inverts the Goff-Gratch equation,
	 *  use Newton's method to improve the inversion to machine precision. */
	double term1, dewPoint, alpha, num, denom;

	term1 = log(SVP*RH/611.2);
	dewPoint = 243.75*term1/(17.67-term1)+273.15;
	/* Now improve the estimate of the dew point using
	 * f(dewpoint) = SVP*RH/100 - goffGratch(dewPoint)
	 * and find where f=0. */
	alpha = SVP*RH/100;
	/* We start off close so 4 iterations of Newton's method is plenty. */
	num = (alpha-goffGratch(dewPoint)); /* One */
	denom = DgoffGratch_dT(dewPoint);
	if(is_unstable(num,denom,dewPoint))
		return dewPoint;
	dewPoint += num/denom;

	num = (alpha-goffGratch(dewPoint)); /* Two */
	denom = DgoffGratch_dT(dewPoint);
	if(is_unstable(num,denom,dewPoint))
		return dewPoint;
	dewPoint += num/denom;

	num = (alpha-goffGratch(dewPoint)); /* Three */
	denom = DgoffGratch_dT(dewPoint);
	if(is_unstable(num,denom,dewPoint))
		return dewPoint;
	dewPoint += num/denom;

	num = (alpha-goffGratch(dewPoint)); /* Four */
	denom = DgoffGratch_dT(dewPoint);
	if(is_unstable(num,denom,dewPoint))
		return dewPoint;
	dewPoint += num/denom;

	/*for(i=0;i<4;i++)
		dewPoint += (alpha-goffGratch(dewPoint))/DgoffGratch_dT(dewPoint);*/

	return dewPoint;
}
double goffGratch(double T){
	/* Goff Gatch equation returns vapor pressure (in millibars) from temperature in Kelvin.
	 * References:
     *http://en.wikipedia.org/wiki/Goff%E2%80%93Gratch_equation
     *http://cires.colorado.edu/~voemel/vp.html*/
    double TS = 373.16,loge;
    loge =  -7.90298 * ((TS / T) - 1.0) +
             5.02808 * log10(TS / T) -
             1.3816e-7 * (pow(10.0 , 11.344 * (1.0 - (T / TS))) - 1.0) +
             8.1328e-3 * (pow(10.0 , (-3.49149 * (TS / T - 1.0))) - 1.0) +
             log10(101324.6);
    return pow(10.0,loge-2.0);
}
double DgoffGratch_dT(double T){
	/* Calculate the derivative of the Goff Gratch equation with respect to
	 * temperature. */
	double TS = 373.16, out, a, b, c, d;
	out = log(10.0)*goffGratch(T);
	a = 7.90298*TS/T/T;
	b = -5.02808/T/log(10.0);
	c = -1.3816e-7*pow(10.0,11.344*(1.0 - T/TS))*log(10.0)*-11.344/TS;
	d = 8.1328e-3*pow(10.0,-3.49149*(TS/T - 1.0)) * log(10.0) * 3.49149*TS/T/T;
	out *= (a+b+c+d);
	return out;
}
double goffGratchIce(double T){
	/* Goff Gatch equation returns vapor pressure (in millibars) from temperature in Kelvin.
	 * References:
     *http://en.wikipedia.org/wiki/Goff%E2%80%93Gratch_equation
     *http://cires.colorado.edu/~voemel/vp.html*/
    double TS = 273.16,loge;
    loge =  -9.09718 * ((TS / T) - 1.0) -
             3.56654 * log10(TS / T) +
             0.876793 * (1.0 - (T/TS)) +
             log10(6.1071);
    return pow(10.0,loge-2.0);
}
double DgoffGratchIce_dT(double T){
	/* Calculate the derivative of the Ice Goff Gratch equation with respect to
	 * temperature. */
	double TS = 273.16, out, a, b, c;
	out = log(10.0)*goffGratch(T);
	a = 9.09718*TS/T/T;
	b = 3.56654/T/log(10.0);
	c = -0.876793/TS;
	out *= (a+b+c);
	return out;
}
double enhancementFactor(double T, double P){
	/* Find the enhancement factor from temperature (Kelvin) and Pressure (mb).
	 * Taken from Phillip E. Ciddor,
	 *	"Refractive index of air: new equations for the visible and
	 *	near infrared" Applied Optics, Vol 35, No 9. 20 March 1996.
	 *	beta is modified from Ciddor's paper to accomidate pressure
	 *	in mb.
	 *
	 *	 Aslo at NIST http://www.nist.gov/calibrations/upload/metv29i1p67-2.pdf
	 *
	 *	Inputs :
	 *	T temperature in Kelvin
	 *	P pressure in mb
	 *
	 *
	 *	 Outputs:
	 *	 f     - Enhancement factor for water vapor in moist air. (unitless)
	 */

	double alpha = 1.00062,
		   beta = 3.14e-6,
		   gamma_= 5.60e-7;

	return alpha + beta*P + gamma_*pow(T-273.15,2.0);

}
double calcAbsolute(double T, double P, double moleMixingRatio){
   /* %CALCABSOLUTEHUMIDITY Calculate absolute humidity
    % Calculate the absolute humidity for the given pressure, mole mixing ratio and
    % temperature. Not that many people use absolute humidity, but here it
    % is.
    %
    % Inputs:
    %   P	pressure (mB or hPa)
    %   T	temperature (kelvin)
    	moleMixingRatio
    %
    % Outputs:
    %   absoluteHumidity (g/m^3)
    %
    % References:
    %   See NIST http://www.nist.gov/calibrations/upload/metv29i1p67-2.pdf */
	double beta,Z; /* Number density of moist air, Compressibility of moist air */

	Z = calcCompressibility(T,P,moleMixingRatio);
	beta = P*100/Z/T/8.314598;

	return beta*moleMixingRatio*18.01528;
}
double calcCompressibility(double T, double P, double xv){
	/* Find the compressibility of moist air.
	 * Taken from Phillip E. Ciddor,
	 *	"Refractive index of air: new equations for the visible and
	 *	near infrared" Applied Optics, Vol 35, No 9. 20 March 1996.
	 *
	 * Inputs are
	 * T 	Temperature (K)
	 * P	Pressure (mb)
	 * xv	mole mixing ratio of water vapor to moist air. */
	double t,p,a[3],b[2],c[2],d,e;

	/* Set dimensions and convert units*/
	t = T - 273.15;
	p = P*100;

	/* Ciddor Equation Constasts */
	a[0]=1.58123e-6;  a[1]=-2.9331e-8;  a[2]=1.1043e-10;
	b[0]=5.707e-6;    b[1]=-2.051e-8;
	c[0]=1.9898e-4;   c[1]=-2.376e-6;
	d=1.83e-11;     e=-0.765e-8;

	return 1.0 + pow((p/T),2.0)*(d+e*xv*xv) - p/T*(a[0] + a[1]*t +
	    a[2]*t*t + (b[0] + b[1]*t)*xv + (c[0] + c[1]*t)*xv*xv);
}
double calcPotentialTemperature(double T, double P, double P0, CALCULATION_DIRECTION D){
	return T*pow(P0/P,2.0/7.0*(double)D);
}
double calcVirtualTemperature(double T, double mr, CALCULATION_DIRECTION D){
	/* See http://glossary.ametsoc.org/wiki/Virtual_temperature */
	if (D)
		return T * (1.0 + 0.61 * mr/1000.0);
	else
		return T / (1.0 + 0.61 * mr/1000.0);
}
double calcMoistAirDensity(double T, double P, double xv, double xCO2){
	/* Calculate the moist air density using the method from Ciddor's paper.
	 * Inputs:
	 * T	Temperatuer (K)
	 * P	Pressure (mb)
	 * xv 	Mole mixing ratio of water vapor to moist air
	 * xCO2	Mole mixing ratio of CO2 to dry air (ppm) */
	double Ma,Z,R=8.314510; /* Mass of dry air, compressibility of moist air, gas constant */

	Ma = 28.9635 + 12.011e-6*(xCO2 - 400.0); /* Molar mass of dry air */
	Z = calcCompressibility(T,P,xv);

	return (P*100.0*Ma/Z/T/R);
}
double moistAirNumberDensity(double T, double P, double xv){
	/* Moist air number density in mol / m^3. Also called "beta"*/
	double Z;

	Z = calcCompressibility(T,P,xv);

	return ((100*P)/(Z*T*8.314598));
}
double dZ_dxv(double T, double p, double x){
	/* The first derivative of compressibility with respect to molar mixing ratio */
	double t,b0,b1,c0,c1,e;
	t = T-273.15;
	p *= 100.0; /*convert to pascals*/

	/* Ciddor constants */
	/*a0=1.58123e-6;  a1=-2.9331e-8;  a2=1.1043e-10;*/
	b0=5.707e-6;    b1=-2.051e-8;
	c0=1.9898e-4;   c1=-2.376e-6;
	/*d=1.83e-11;*/	e=-0.765e-8;

	return (p*pow(0.27315e3+t,-2.0) * (2*e*p*x - (0.27315e3+t)*(b0+b1*t + 2.0*(c0 + c1*t)*x)));
}
double dBeta_dxv(double T, double P, double xv){
	/* The partial of the moist air number density with respect to molar mixing ratio */
	double Z;
	Z = calcCompressibility(T,P,xv);

	return (-moistAirNumberDensity(T,P,xv)/Z*dZ_dxv(T,P,xv));
}
double calcSpecificHumidity(double mr){
	/*
	 *
    %CALCSPECIFICHUMIDITY Calculate specific humidity
    % Calculate the specific humidity based upon a given mixing ratio.
    %
    % Inputs:
    %   mixingRatio      [vector] = Mixing Ratio (g / kg)
    %
    % Outputs:
    %   specificHumidity [vector] = Specific Humidity (g/kg)
    %
    % References:
    %   Wallace, John M., Atmsopheric Science: An Introductor Survey 2nd Edition, 2006 - 3.5.1 Unlabeled Equation
    %   http://hurri.kean.edu/~yoh/calculations/moisture/Equations/moist.html */

	double sh;
	/* convert to g/g*/
	sh = mr/1000.0;
	sh = sh/(1.0+sh);
	/* convert back to g/kg*/
	return sh*1000.0;
}
double calcMassMixingRatio(double P, double vp){
	/*     %CALCMIXINGRATIO Calculate mixing ratio
    % Calculate mixing ratio for a specified pressure and temperature. For
    % saturated mixing ratios use temperature. For standard mixing ratios
    % use dewpoint. This may become unstable in rare cases where the
    % pressure equals the vapor pressure.
    %
    % Inputs:
    %   P    Atmospheric pressure (mb or hPa)
    %   vp		vapor pressure (mb or hPa)
    %
    % Outputs:
    %   ratio       [vector] = The mixing ratio (g/kg)
    %
    % References:
    %   http://www.wrh.noaa.gov/slc/projects/wxcalc/formulas/mixingRatio.pdf*/
	double epsilon = 0.62197;

	/* Pressure should always be greater than vapor pressure */
	P = P > 1.00001*vp ? P:vp*1.00001;

	return 1000.0*epsilon*vp/(P-vp);
}
double calcMMRfromAbsoluteHumidity(double P, double T, double a){
	double R = 8.314598; /* Gas constant */
	double Mv = 18.01528; /* Molar mass of H20 */
	double beta0,p,xv,num,denom;
	int i;

	p = P*100; /* Convert to pascals */
	beta0 = p/R/T; /* Initial guess of density from ideal gas law. */
	xv = a/Mv/beta0; /* Initial guess of molar mixing ratio */
	for(i=0;i<6;i++){ /* Six iterations of newton's method */
		num = a- moistAirNumberDensity(T,P,xv)*xv*Mv;
		denom = dBeta_dxv(T,P,xv)*xv*Mv + moistAirNumberDensity(T,P,xv)*Mv;
		xv = xv + num/denom;
	}
	return xv;
}
double calcSHfromRH(double T, double P, double xCO2, double SVP, double RH, double ef){
	/* calculate Specific Humidity from Relative Humidity */
	double MMR,AH,MAD;
	/* First we need the mole mixing ratio */
	MMR =  ef*RH*SVP*P/100.0;
	/* Then the absolute humidity */
	AH = calcAbsolute(T,P,MMR);
	/* Next, the moist air density */
	MAD = calcMoistAirDensity(T,P,MMR,xCO2);

	return 1000.0*AH/MAD;

}
double estRHfromSH(double SH, double T, double P, double xCO2, double SVP, double ef, double RH){
	/* Primarily for cases where the SH is near 1 so, finding RH via Mass Mixing
	 * Ratio is numerically unstable, this method starts with an estimate of the
	 * RH, and refines the estimate using Newton's Method. */
	int step;
	double RHp,RHm,del = 1.0,dS_dR,dRH,dS,num;

	for(step=0;step<6;step++){
		/* use centered differencing to get the slope */
		RHm = (RH - del) >= 0.0 ? RH-del : RH;
		RHp = (RH + del) <= 1.0 ? RH+del : RH;

		/* Find the dSH / dRH. */
		dRH = (RHp - RHm);
		dS = calcSHfromRH(T,P,xCO2,SVP,RHp,ef) - calcSHfromRH(T,P,xCO2,SVP,RHm,ef);

		if (is_unstable(dS,dRH,RH)) break;
		else
			dS_dR =  dS/dRH;

		/* Get the update in RH */
		num = (calcSHfromRH(T,P,xCO2,SVP,RH,ef) - SH);

		if(is_unstable(num,dS_dR,RH)) break;
		else
			del = num/dS_dR;

		/* update RH */
		RH -= del;

		if(is_unstable(del,RH,RH)) break;

		/* Decide to continue */
		if (fabs(del / RH) < 1e-12) break;
	}

	return RH;

}
int findStartIndex(double *x, double x0,int N){
	int start_i;
	/* Find the first entry above z0 */
	for(start_i=0;start_i<N;start_i++){
		if(x0 <= x[start_i])
			break;
	}
	return start_i;
}
int findEndIndex(double *x,double x0,int N){
	int end_i;
	/* Find the last entry below x0 */
	for(end_i=N-1;end_i>0;end_i--){
		if(x0 >= x[end_i])
			break;
	}
	return end_i;
}
double logarithmicRule(double *x, double *y, int i){
	/* Integrate using a logarithmic interpolant. This expression can be found
	 * by integrating a logarithmic function of the form exp(a*x + b) from
	 * x = x0 to x = x1 after finding a & b so that the function passes through
	 * y(x0) and y(x1). */
	double raty = y[i+1]/y[i];
	double dx = x[i+1] - x[i];
	if ((raty == 1.0) || isinf(raty) || isnan(raty))
		/* If something went wrong, use linear interpolation. */
		return dx*(y[i]);
	else
		return y[i]*dx/log(raty)*(raty-1.0);
}
double logarithmicPartial(double *x, double *y, int i,double a, double b){
	/* Integrate using a logarithmic interpolant. This expression can be found
	 * by integrating a logarithmic function of the form exp(a*x + b) from
	 * x = x0 to x = x1 after finding a & b so that the function passes through
	 * y(x0) and y(x1). */
	double raty = y[i+1]/y[i];
	double dx = x[i+1] - x[i];
	if (raty == 1.0)
		return dx*y[i];
	else
		return y[i]*dx/log(raty)*(pow(raty,(b - x[i])/dx)-pow(raty,(a-x[i])/dx));
}
double trapezoidalRule(double *x, double *y, int i){
	/* Integrate using the trapezoidal rule . */
	double dx = x[i+1] - x[i];
	return 0.5*(dx*y[i] + dx*y[i+1]);
}
double trapezoidalPartial(double *x, double *y, int i,double a, double b){
	/* Integrate using the trapezoidal rule . */
	double alpha = (y[i+1] - y[i])/(x[i+1] - x[i]);
	double xa = a - x[i];
	double xb = b - x[i];
	return 0.5*alpha*(xb*xb - xa*xa) + y[i]*(xb-xa);
}
double integrateColumnLogarithmic(double *x, double *y,double x0, double x1, unsigned int N) {
	/* Integrate y over the domain of x from x0 to x1 using a log-linear interpolant */
	unsigned int start_i, end_i, i;
	double out = 0.0;

	start_i = findStartIndex(x, x0, N);
	if (start_i >= N - 1) return out;

	end_i = findEndIndex(x, x1, N);
	if (end_i <= 0) return out;

	if (start_i > 0)
		out += logarithmicPartial(x,
			y,
			start_i - 1,
			x0,
			x[start_i]);
	if (end_i < N - 1)
		out += logarithmicPartial(x,
			y,
			end_i,
			x[end_i],
			x1);
	for (i = start_i; i<end_i; i++)
		out += logarithmicRule(x,
			y,
			i);
	return out;
}
double integrateColumnTrapezoid(double *x, double *y, double x0, double x1, unsigned int N) {
	/* Integrate y over the domain of x from x0 to x1 using a linear interpolant */
	unsigned int start_i, end_i, i;
	double out = 0.0;

	start_i = findStartIndex(x, x0, N);
	if (start_i >= N - 1) return out;

	end_i = findEndIndex(x, x1, N);
	if (end_i <= 0) return out;

	if (start_i > 0)
		out += trapezoidalPartial(x,
			y,
			start_i - 1,
			x0,
			x[start_i]);
	if (end_i < N - 1)
		out += trapezoidalPartial(x,
			y,
			end_i,
			x[end_i],
			x1);
	for (i = start_i; i<end_i; i++)
		out += trapezoidalRule(x,
			y,
			i);
	return out;
}
double integrateColumnWaterDensity(WEATHER_CONVERSION_VECTOR *WX,double z0, double z1){
	/* Integrate the Column Water Vapor Density */
	if(WX->populated[_ABSOLUTE_HUMIDITY] && WX->populated[_HEIGHT_AGL])
		return integrateColumnLogarithmic(WX->val[_HEIGHT_AGL], WX->val[_ABSOLUTE_HUMIDITY],z0, z1, WX->N);
	else {
		if (WX->populated[_ABSOLUTE_HUMIDITY]) {
			if(setHeights(WX)!=WEATHER_CONVERSION_SUCCESS)
				printf("weatherConversion.c:integrateColumnWaterDenisty(): WARNING: Unable to calculate the column water vapor density because _HEIGHT_AGL is not populated.\n");
			else
				return integrateColumnLogarithmic(WX->val[_HEIGHT_AGL], WX->val[_ABSOLUTE_HUMIDITY], z0, z1, WX->N);
		}			
		else if(WX->populated[_HEIGHT_AGL])
			printf("weatherConversion.c:integrateColumnWaterDenisty(): WARNING: Unable to calculate the column water vapor density because _ABSOLUTE_HUMIDITY is not populated.\n");
		else
			printf("weatherConversion.c:integrateColumnWaterDenisty(): WARNING: Unable to calculate the column water vapor density because HEIGHT_AGL and _ABSOLUTE_HUMIDITY are not populated.\n");
		return 0.0;
	}
}
double integrateColumnWaterNumberDensity(WEATHER_CONVERSION_VECTOR *WX, double z0, double z1) {
	return integrateColumnWaterDensity(WX, z0, z1) / WATER_MOLAR_MASS;
}
double integrateColumnMoistAirDensity(WEATHER_CONVERSION_VECTOR *WX, double z0, double z1) {
	/* Integrate the Column Moist Air Density */
	if (WX->populated[_MOIST_AIR_DENSITY] && WX->populated[_HEIGHT_AGL])
		return integrateColumnLogarithmic(WX->val[_HEIGHT_AGL], WX->val[_MOIST_AIR_DENSITY], z0, z1, WX->N);
	else {
		if (WX->populated[_MOIST_AIR_DENSITY]) {
			if (setHeights(WX) != WEATHER_CONVERSION_SUCCESS)
				printf("weatherConversion.c:integrateColumnWaterDenisty(): WARNING: Unable to calculate the column water vapor density because _HEIGHT_AGL is not populated.\n");
			else
				return integrateColumnLogarithmic(WX->val[_HEIGHT_AGL], WX->val[_MOIST_AIR_DENSITY], z0, z1, WX->N);
		}
		else if (WX->populated[_HEIGHT_AGL])
			printf("weatherConversion.c:integrateColumnWaterDenisty(): WARNING: Unable to calculate the column water vapor density because _MOIST_AIR_DENSITY is not populated.\n");
		else
			printf("weatherConversion.c:integrateColumnWaterDenisty(): WARNING: Unable to calculate the column water vapor density because HEIGHT_AGL and _MOIST_AIR_DENSITY are not populated.\n");
		return 0.0;
	}
}
double integrateColumnMoistAirNumberDensity(WEATHER_CONVERSION_VECTOR *WX, double z0, double z1) {
	/* Integrate the Column Moist Air Number Density */
	if (WX->populated[_MOIST_AIR_NUMBER_DENSITY] && WX->populated[_HEIGHT_AGL])
		return integrateColumnLogarithmic(WX->val[_HEIGHT_AGL], WX->val[_MOIST_AIR_NUMBER_DENSITY], z0, z1, WX->N);
	else {
		if (WX->populated[_MOIST_AIR_NUMBER_DENSITY]) {
			if (setHeights(WX) != WEATHER_CONVERSION_SUCCESS)
				printf("weatherConversion.c:integrateColumnWaterDenisty(): WARNING: Unable to calculate the column water vapor density because _HEIGHT_AGL is not populated.\n");
			else
				return integrateColumnLogarithmic(WX->val[_HEIGHT_AGL], WX->val[_MOIST_AIR_NUMBER_DENSITY], z0, z1, WX->N);
		}
		else if (WX->populated[_HEIGHT_AGL])
			printf("weatherConversion.c:integrateColumnWaterDenisty(): WARNING: Unable to calculate the column water vapor density because _MOIST_AIR_NUMBER_DENSITY is not populated.\n");
		else
			printf("weatherConversion.c:integrateColumnWaterDenisty(): WARNING: Unable to calculate the column water vapor density because HEIGHT_AGL and _MOIST_AIR_NUMBER_DENSITY are not populated.\n");
		return 0.0;
	}
}
WEATHER_CONVERTER_FIELD presentHumidity(WEATHER_CONVERSION_VECTOR *WX){
	WEATHER_CONVERTER_FIELD fi;
	for(fi=_RELATIVE_HUMIDITY;fi<_N_WEATHER_FIELDS;fi++)
		if(WX->populated[fi]) return fi;
	return _N_WEATHER_FIELDS;
}
WEATHER_CONVERSION_ERROR setAllFields(WEATHER_CONVERSION_VECTOR *WX){
	WEATHER_CONVERSION_ERROR retVal;
	unsigned int i;
	/* Winds do not depend on other variables. Attempt to set them
	 * first */
	if((retVal=setWinds(WX))!=WEATHER_CONVERSION_SUCCESS)
		printf("weatherConversion.h:setAllFields() unable to set wind fields.\nContinuing on to other conversions.\n");

	if((retVal=setTemperatures(WX))!=WEATHER_CONVERSION_SUCCESS){
		printf("weatherConversion.h:setAllFields() unable to set temperature fields.\nAborting further conversions.\n");
		return retVal;
	}
	if((retVal=setPressures(WX))!=WEATHER_CONVERSION_SUCCESS){
		printf("weatherConversion.h:setAllFields() unable to set pressure fields.\nAborting further conversions.\n");
		return retVal;
	}
	switch (presentHumidity(WX)){
	case _RELATIVE_HUMIDITY:
		/* This is the value we need. Stop here. */
		break;
	case _VAPOR_PRESSURE:
		if(WX->populated[_RELATIVE_HUMIDITY]==FALSE)
			for(i=0;i<WX->N;i++)
				WX->val[_RELATIVE_HUMIDITY][i] = 100.0*WX->val[_VAPOR_PRESSURE][i]/WX->val[_SATURATION_VAPOR_PRESSURE][i];
		WX->populated[_RELATIVE_HUMIDITY] = TRUE;
		break;
	case _POTENTIAL_VAPOR_PRESSURE:
		if(WX->populated[_VAPOR_PRESSURE]==FALSE)
			for(i=0;i<WX->N;i++)
				WX->val[_VAPOR_PRESSURE][i] = WX->val[_POTENTIAL_VAPOR_PRESSURE][i]*WX->val[_PRESSURE][i]/WX->standardPressure;
		if(WX->populated[_RELATIVE_HUMIDITY]==FALSE)
			for(i=0;i<WX->N;i++)
				WX->val[_RELATIVE_HUMIDITY][i] = 100.0*WX->val[_VAPOR_PRESSURE][i]/WX->val[_SATURATION_VAPOR_PRESSURE][i];
		WX->populated[_VAPOR_PRESSURE]=TRUE;
		WX->populated[_RELATIVE_HUMIDITY] = TRUE;
		break;
	case _MOLE_MIXING_RATIO:
		if(WX->populated[_RELATIVE_HUMIDITY]==FALSE)
			for(i=0;i<WX->N;i++)
				WX->val[_RELATIVE_HUMIDITY][i] = 100.0*WX->val[_MOLE_MIXING_RATIO][i]*WX->val[_PRESSURE][i] /
											WX->val[_SATURATION_VAPOR_PRESSURE][i]/WX->val[_ENHANCEMENT_FACTOR][i];
		WX->populated[_RELATIVE_HUMIDITY] = TRUE;
		break;
	case _MASS_MIXING_RATIO:
		if(WX->populated[_RELATIVE_HUMIDITY]==FALSE){
			for(i=0;i<WX->N;i++)
					WX->val[_RELATIVE_HUMIDITY][i] =
							100.0*WX->val[_MASS_MIXING_RATIO][i]/WX->val[_SATURATION_MIXING_RATIO][i];
		}
		WX->populated[_RELATIVE_HUMIDITY] = TRUE;
		break;
	case _DEW_POINT_F:
		for(i=0;i<WX->N;i++){
			WX->val[_DEW_POINT_K][i] = FtoK(WX->val[_DEW_POINT_F][i]);
			WX->val[_DEW_POINT_C][i] = FtoC(WX->val[_DEW_POINT_F][i]);
		}
		if(WX->populated[_VAPOR_PRESSURE]==FALSE)
			for(i=0;i<WX->N;i++)
				WX->val[_VAPOR_PRESSURE][i] = goffGratch(WX->val[_DEW_POINT_K][i] );
		if(WX->populated[_RELATIVE_HUMIDITY]==FALSE)
			for(i=0;i<WX->N;i++)
				WX->val[_RELATIVE_HUMIDITY][i] = 100.0*WX->val[_VAPOR_PRESSURE][i]/WX->val[_SATURATION_VAPOR_PRESSURE][i];
		WX->populated[_DEW_POINT_K] =
				WX->populated[_DEW_POINT_C] =
						WX->populated[_VAPOR_PRESSURE] =
								WX->populated[_RELATIVE_HUMIDITY] = TRUE;
		break;
	case _DEW_POINT_C:
		for(i=0;i<WX->N;i++){
			WX->val[_DEW_POINT_K][i] = CtoK(WX->val[_DEW_POINT_C][i]);
			WX->val[_DEW_POINT_F][i] = CtoF(WX->val[_DEW_POINT_C][i]);
		}
		if(WX->populated[_VAPOR_PRESSURE]==FALSE)
			for(i=0;i<WX->N;i++)
				WX->val[_VAPOR_PRESSURE][i] = goffGratch(WX->val[_DEW_POINT_K][i] );
		if(WX->populated[_RELATIVE_HUMIDITY]==FALSE)
			for(i=0;i<WX->N;i++)
				WX->val[_RELATIVE_HUMIDITY][i] = 100.0*WX->val[_VAPOR_PRESSURE][i]/WX->val[_SATURATION_VAPOR_PRESSURE][i];
		WX->populated[_DEW_POINT_K] =
				WX->populated[_DEW_POINT_F] =
						WX->populated[_VAPOR_PRESSURE] =
								WX->populated[_RELATIVE_HUMIDITY] = TRUE;
		break;
	case _DEW_POINT_K:
		if(WX->populated[_VAPOR_PRESSURE]==FALSE)
			for(i=0;i<WX->N;i++)
				WX->val[_VAPOR_PRESSURE][i] = goffGratch(WX->val[_DEW_POINT_K][i] );
		if(WX->populated[_RELATIVE_HUMIDITY]==FALSE)
			for(i=0;i<WX->N;i++)
				WX->val[_RELATIVE_HUMIDITY][i] = 100.0*WX->val[_VAPOR_PRESSURE][i]/WX->val[_SATURATION_VAPOR_PRESSURE][i];
		WX->populated[_VAPOR_PRESSURE] = TRUE;
		WX->populated[_RELATIVE_HUMIDITY] = TRUE;
		break;
	case _SPECIFIC_HUMIDITY:
		if(WX->populated[_MASS_MIXING_RATIO]==FALSE)
			for(i=0;i<WX->N;i++)
				WX->val[_MASS_MIXING_RATIO][i] = WX->val[_SPECIFIC_HUMIDITY][i]/(1.0-WX->val[_SPECIFIC_HUMIDITY][i]/1000.0);
		if(WX->populated[_RELATIVE_HUMIDITY]==FALSE)
			for(i=0;i<WX->N;i++){
				/*if (WX->val[_SPECIFIC_HUMIDITY][i] > 0.99)
					WX->val[_RELATIVE_HUMIDITY][i] = estRHfromSH(
							WX->val[_SPECIFIC_HUMIDITY][i],
							WX->val[_TEMPERATURE_K][i],
							WX->val[_PRESSURE][i],
							WX->xCO2,
							WX->val[_SATURATION_VAPOR_PRESSURE][i],
							WX->val[_ENHANCEMENT_FACTOR][i],
							100);
				else*/
					WX->val[_RELATIVE_HUMIDITY][i] =
						100.0*WX->val[_MASS_MIXING_RATIO][i]/WX->val[_SATURATION_MIXING_RATIO][i];

			}
		WX->populated[_MASS_MIXING_RATIO] = TRUE;
		WX->populated[_RELATIVE_HUMIDITY] = TRUE;
		break;
	case _ABSOLUTE_HUMIDITY:
		if(WX->populated[_MOLE_MIXING_RATIO]==FALSE)
			for(i=0;i<WX->N;i++)
				WX->val[_MOLE_MIXING_RATIO][i] = calcMMRfromAbsoluteHumidity(WX->val[_PRESSURE][i], WX->val[_TEMPERATURE_K][i],
					WX->val[_ABSOLUTE_HUMIDITY][i]);
		if(WX->populated[_RELATIVE_HUMIDITY]==FALSE)
			for(i=0;i<WX->N;i++)
				WX->val[_RELATIVE_HUMIDITY][i] = 100.0*(WX->val[_MOLE_MIXING_RATIO][i]*WX->val[_PRESSURE][i]) /
														(WX->val[_SATURATION_VAPOR_PRESSURE][i]*WX->val[_ENHANCEMENT_FACTOR][i]);

		WX->populated[_MOLE_MIXING_RATIO] = TRUE;
		WX->populated[_RELATIVE_HUMIDITY] = TRUE;
		break;
	case _MOIST_AIR_DENSITY:
		printf("weatherConversion.h:humidityConversion() Conversion from Moist air density is NOT SUPPORTED YET!! But it is possible.\n");
		for(i=0;i<WX->N;i++){
		}
		return NO_HUMIDITY_PRESENT;
		break;
	case _MOIST_AIR_NUMBER_DENSITY:
		printf("weatherConversion.h:humidityConversion() Conversion from Moist air number density is NOT SUPPORTED YET!! But it is possible.\n");
		for (i = 0; i<WX->N; i++) {
		}
		return NO_HUMIDITY_PRESENT;
		break;
	case _WATER_VAPOR_NUMBER_DENSITY:
		printf("weatherConversion.h:humidityConversion() Conversion from Water vapor number density is NOT SUPPORTED YET!! But it is possible.\n");
		for (i = 0; i<WX->N; i++) {
		}
		return NO_HUMIDITY_PRESENT;
		break;
	default:
		printf("weatherConversion.h:setAllFields() unable to set humidity fields.\nAborting conversions.\n");
		return NO_HUMIDITY_PRESENT;
		break;
	}
	/* Starting from relative Humidity, perform all conversions. */
	humidityConversion(WX);
	/* Density and pressure may be required to get the heights, so do this now. */
	if((retVal=setHeights(WX))!=WEATHER_CONVERSION_SUCCESS){
		printf("weatherConversion.h:setAllFields() unable to set height fields.\nIntegrating column values is not possible.\n");
	}
	else{
	/* Find the integrated column densities from the ground level up to the Top of the atmosphere. */
	WX->vaporArealDensity = WX->f.integrate_column_water_density(WX,
			                0.0,
							WX->val[_HEIGHT_AGL][WX->N-1]);
	WX->vaporArealNumberDensity =	WX->f.integrate_column_water_number_density(WX,
									0.0,
									WX->val[_HEIGHT_AGL][WX->N - 1]);
	WX->moistAirArealDensity =	WX->f.integrate_column_moist_air_density(WX,
								0.0,
								WX->val[_HEIGHT_AGL][WX->N - 1]);
	WX->moistAirArealNumberDensity =	WX->f.integrate_column_moist_air_number_density(WX,
										0.0,
										WX->val[_HEIGHT_AGL][WX->N - 1]);
	}
	return retVal;
}
WEATHER_CONVERSION_ERROR setTemperatures(WEATHER_CONVERSION_VECTOR *WX){
	unsigned int i;
	/* Not that its likely, but maybe we got potential Temp and pressure only. */
	if(WX->populated[_PRESSURE] && WX->populated[_POTENTIAL_TEMPERATURE]){
		if(WX->populated[_TEMPERATURE_K]==FALSE)
			for(i=0;i<WX->N;i++)
				WX->val[_TEMPERATURE_K][i] = calcPotentialTemperature(WX->val[_POTENTIAL_TEMPERATURE][i], WX->val[_PRESSURE][i], WX->standardPressure, BACKWARD);
		WX->populated[_TEMPERATURE_K] = TRUE;
	}

	if(WX->populated[_TEMPERATURE_C]){
		if(WX->populated[_TEMPERATURE_F]==FALSE)
			for(i=0;i<WX->N;i++)
				WX->val[_TEMPERATURE_F][i] = CtoF(WX->val[_TEMPERATURE_C][i]);
		if(WX->populated[_TEMPERATURE_K]==FALSE)
			for(i=0;i<WX->N;i++)
				WX->val[_TEMPERATURE_K][i] = WX->val[_TEMPERATURE_C][i] + 273.15;

		WX->populated[_TEMPERATURE_F] = TRUE;
		WX->populated[_TEMPERATURE_K] = TRUE;
	}
	else if(WX->populated[_TEMPERATURE_F]){
		if(WX->populated[_TEMPERATURE_C]==FALSE)
			for(i=0;i<WX->N;i++)
				WX->val[_TEMPERATURE_C][i] = FtoC(WX->val[_TEMPERATURE_F][i]);
		if(WX->populated[_TEMPERATURE_K]==FALSE)
			for(i=0;i<WX->N;i++)
				WX->val[_TEMPERATURE_K][i] = WX->val[_TEMPERATURE_C][i] + 273.15;
		WX->populated[_TEMPERATURE_C] = TRUE;
		WX->populated[_TEMPERATURE_K] = TRUE;
	}
	else if(WX->populated[_TEMPERATURE_K]){
		if(WX->populated[_TEMPERATURE_C]==FALSE)
			for(i=0;i<WX->N;i++)
				WX->val[_TEMPERATURE_C][i] = WX->val[_TEMPERATURE_K][i] - 273.15;
		if(WX->populated[_TEMPERATURE_F]==FALSE)
			for(i=0;i<WX->N;i++)
				WX->val[_TEMPERATURE_F][i] = CtoF(WX->val[_TEMPERATURE_C][i]);

		WX->populated[_TEMPERATURE_C] = TRUE;
		WX->populated[_TEMPERATURE_F] = TRUE;
	}
	else
		return NO_TEMPERATURE_PRESENT;

	if(WX->populated[_SATURATION_VAPOR_PRESSURE]==FALSE)
		for(i=0;i<WX->N;i++)
			WX->val[_SATURATION_VAPOR_PRESSURE][i] = goffGratch(WX->val[_TEMPERATURE_K][i]);
	WX->populated[_SATURATION_VAPOR_PRESSURE] = TRUE;
	/* Saturation mixing ratio is set with the pressures. */
	return WEATHER_CONVERSION_SUCCESS;
}
WEATHER_CONVERSION_ERROR setPressures(WEATHER_CONVERSION_VECTOR *WX){
	unsigned int i;
	if(WX->populated[_PRESSURE]){
		if(WX->populated[_POTENTIAL_TEMPERATURE]==FALSE)
			for(i=0;i<WX->N;i++)
				WX->val[_POTENTIAL_TEMPERATURE][i] = calcPotentialTemperature(WX->val[_TEMPERATURE_K][i], WX->val[_PRESSURE][i],WX->standardPressure,FOREWARD);
		WX->populated[_POTENTIAL_TEMPERATURE]=TRUE;
	}
	else if(WX->populated[_POTENTIAL_TEMPERATURE] && WX->populated[_TEMPERATURE_K]){
		if(WX->populated[_PRESSURE]==FALSE)
			for(i=0;i<WX->N;i++)
				WX->val[_PRESSURE][i] = calcPotentialTemperature(WX->standardPressure, WX->val[_TEMPERATURE_K][i],WX->val[_POTENTIAL_TEMPERATURE][i],FOREWARD);
		WX->populated[_PRESSURE]=TRUE;
		}
	else
		return NO_PRESSURE_PRESENT;
	if(WX->populated[_ENHANCEMENT_FACTOR]==FALSE)
		for(i=0;i<WX->N;i++)
			WX->val[_ENHANCEMENT_FACTOR][i] = enhancementFactor( WX->val[_TEMPERATURE_K][i],WX->val[_PRESSURE][i]);
	if(WX->populated[_SATURATION_MIXING_RATIO]==FALSE){
		for(i=0;i<WX->N;i++){
			/* If the actual pressure is less than the saturation pressure, then assume that there is no dry air. * /
			if(WX->val[_PRESSURE][i] < WX->val[_SATURATION_VAPOR_PRESSURE][i])
				WX->val[_SATURATION_MIXING_RATIO][i] = _MY_INFINITY;
			else*/
				WX->val[_SATURATION_MIXING_RATIO][i] = calcMassMixingRatio(WX->val[_PRESSURE][i],WX->val[_SATURATION_VAPOR_PRESSURE][i]);
						/*622.0*WX->val[_SATURATION_VAPOR_PRESSURE][i]/(WX->val[_PRESSURE][i]-WX->val[_SATURATION_VAPOR_PRESSURE][i]);*/
		}
	}

	WX->populated[_ENHANCEMENT_FACTOR]=TRUE;
	WX->populated[_SATURATION_MIXING_RATIO]=TRUE;
	return WEATHER_CONVERSION_SUCCESS;
}
double hydrostaticDZ(WEATHER_CONVERSION_VECTOR *WX, int i){
	double dP,rho,g;
	if((i<1) || (i>=(int)WX->N))return 0.0;
	dP = (WX->val[_PRESSURE][i-1] - WX->val[_PRESSURE][i] )*100.0; /* Convert to Pascals */
	rho =( WX->val[_MOIST_AIR_DENSITY][i] + WX->val[_MOIST_AIR_DENSITY][i-1])/2000.0; /* Convert to kg/m3 */
	g = free_air_gravity(WX->latitude,WX->val[_HEIGHT_AMSL][i-1]);
	return dP/rho/g;
}
WEATHER_CONVERSION_ERROR setHeights(WEATHER_CONVERSION_VECTOR *WX){
	unsigned int i;
	double R = latitude_earth_radius(WX->latitude);
	double max_ratio = 10.0*(1.0 - nextafter(1.0,0.0));
	if(WX->populated[_GEOPOTENTIAL_HEIGHT]){
		/* Convert from Geopotential Height to Above Mean Sea Level and Above Ground level
		 * Theoretically, GPH cannot be the same as the Earth's radius expect where the AMSL height is infinite. */
		for(i=0;i<WX->N;i++){
			if(WX->val[_GEOPOTENTIAL_HEIGHT][i] / R > max_ratio )
				WX->val[_HEIGHT_AMSL][i]  = R*1e16;
			else
				WX->val[_HEIGHT_AMSL][i] = WX->val[_GEOPOTENTIAL_HEIGHT][i]/(1.0 - WX->val[_GEOPOTENTIAL_HEIGHT][i]/R);
			WX->val[_HEIGHT_AGL][i] = WX->val[_HEIGHT_AMSL][i] - WX->surfaceHeight;
		}
		WX->populated[_HEIGHT_AMSL]=WX->populated[_HEIGHT_AGL] = TRUE;
	}
	else if(WX->populated[_HEIGHT_AMSL]){
		/* Convert from height Above Mean Sea Level to Geopotential Height and AGL*/
		for(i=0;i<WX->N;i++){
			WX->val[_GEOPOTENTIAL_HEIGHT][i] = WX->val[_HEIGHT_AMSL][i]/(1.0 + WX->val[_HEIGHT_AMSL][i]/R);
			WX->val[_HEIGHT_AGL][i] = WX->val[_HEIGHT_AMSL][i] - WX->surfaceHeight;
		}
		WX->populated[_GEOPOTENTIAL_HEIGHT]=WX->populated[_HEIGHT_AGL] = TRUE;
	}
	else if(WX->populated[_HEIGHT_AGL]){
		/* Convert from height Above Ground Level to Geopotential Height and AMSL */
		for(i=0;i<WX->N;i++){
			WX->val[_HEIGHT_AMSL][i] = WX->val[_HEIGHT_AGL][i] + WX->surfaceHeight;
			WX->val[_GEOPOTENTIAL_HEIGHT][i] = WX->val[_HEIGHT_AMSL][i]/(1.0 + WX->val[_HEIGHT_AMSL][i]/R);
		}
		WX->populated[_GEOPOTENTIAL_HEIGHT]=WX->populated[_HEIGHT_AMSL] = TRUE;
	}
	else{
		/* Estimate height from the hydrostatic equation */
		if(WX->populated[_PRESSURE] == FALSE) 	return NO_PRESSURE_PRESENT;
		if(WX->populated[_MOIST_AIR_DENSITY]==FALSE) return NO_HUMIDITY_PRESENT;
		WX->val[_HEIGHT_AGL][0] = 0.0;
		WX->val[_HEIGHT_AMSL][0] = WX->surfaceHeight;
		WX->val[_GEOPOTENTIAL_HEIGHT][0] =WX->val[_HEIGHT_AMSL][0]/(1.0 + WX->val[_HEIGHT_AMSL][0]/R);
		for(i=1;i<WX->N;i++){
			WX->val[_HEIGHT_AMSL][i]  = WX->val[_HEIGHT_AMSL][i-1]  + hydrostaticDZ(WX,i);
			WX->val[_GEOPOTENTIAL_HEIGHT][i]= WX->val[_HEIGHT_AMSL][i]/(1.0 + WX->val[_HEIGHT_AMSL][i]/R);
			WX->val[_HEIGHT_AGL][i] = WX->val[_HEIGHT_AMSL][i] - WX->surfaceHeight;
		}
		WX->populated[_GEOPOTENTIAL_HEIGHT]=
				WX->populated[_HEIGHT_AMSL] =
						WX->populated[_HEIGHT_AGL] = TRUE;
	}
	return WEATHER_CONVERSION_SUCCESS;
}
WEATHER_CONVERSION_ERROR setWinds(WEATHER_CONVERSION_VECTOR *WX){
	unsigned int i;
	if(WX->populated[_U_WIND] && WX->populated[_V_WIND]){
		if(WX->populated[_WIND_SPEED]==FALSE)
			for(i=0;i<WX->N;i++)
				WX->val[_WIND_SPEED][i] = sqrt(WX->val[_U_WIND][i]*WX->val[_U_WIND][i] +
										   WX->val[_V_WIND][i]*WX->val[_V_WIND][i]);
		if(WX->populated[_WIND_DIRECTION]==FALSE)
			for(i=0;i<WX->N;i++)
				WX->val[_WIND_DIRECTION][i] = 180.0*atan2(-WX->val[_V_WIND][i],WX->val[_U_WIND][i])/M_PI;

		WX->populated[_WIND_SPEED] = TRUE;
		WX->populated[_WIND_DIRECTION] = TRUE;
	}
	else if (WX->populated[_WIND_SPEED] && WX->populated[_WIND_DIRECTION]){
		if(WX->populated[_U_WIND]==FALSE)
			for(i=0;i<WX->N;i++)
				WX->val[_U_WIND][i] = cos(WX->val[_WIND_DIRECTION][i]*M_PI/180.0)*WX->val[_WIND_SPEED][i];
		if(WX->populated[_V_WIND]==FALSE)
			for(i=0;i<WX->N;i++)
				WX->val[_V_WIND][i] = -sin(WX->val[_WIND_DIRECTION][i]*M_PI/180.0)*WX->val[_WIND_SPEED][i];
		WX->populated[_U_WIND] = TRUE;
		WX->populated[_V_WIND] = TRUE;
	}
	else
		return NO_WIND_PRESENT;

	return WEATHER_CONVERSION_SUCCESS;
}
WEATHER_CONVERSION_ERROR humidityConversion(WEATHER_CONVERSION_VECTOR *WX){
	unsigned int i;
	if(WX->populated[_DEW_POINT_K]==FALSE){
		for(i=0;i<WX->N;i++)
			WX->val[_DEW_POINT_K][i] = calcDewpoint(WX->val[_TEMPERATURE_K][i], WX->val[_RELATIVE_HUMIDITY][i],
					WX->val[_SATURATION_VAPOR_PRESSURE][i]);
		WX->populated[_DEW_POINT_K]=TRUE;
	}
	if(WX->populated[_DEW_POINT_C]==FALSE){
		for(i=0;i<WX->N;i++)
			WX->val[_DEW_POINT_C][i] = KtoC(WX->val[_DEW_POINT_K][i]);
		WX->populated[_DEW_POINT_C]=TRUE;
	}
	if(WX->populated[_DEW_POINT_F]==FALSE){
		for(i=0;i<WX->N;i++)
			WX->val[_DEW_POINT_F][i] = KtoF(WX->val[_DEW_POINT_K][i]);
		WX->populated[_DEW_POINT_F]=TRUE;
	}
	if(WX->populated[_VAPOR_PRESSURE]==FALSE){
		for(i=0;i<WX->N;i++)
			WX->val[_VAPOR_PRESSURE][i] = WX->val[_RELATIVE_HUMIDITY][i]*WX->val[_SATURATION_VAPOR_PRESSURE][i]/100.0;
		WX->populated[_VAPOR_PRESSURE]=TRUE;
	}
	if(WX->populated[_MOLE_MIXING_RATIO]==FALSE){
		for(i=0;i<WX->N;i++)
			WX->val[_MOLE_MIXING_RATIO][i] = WX->val[_ENHANCEMENT_FACTOR][i]*WX->val[_RELATIVE_HUMIDITY][i]*
												WX->val[_SATURATION_VAPOR_PRESSURE][i]/WX->val[_PRESSURE][i]/100.0;
		WX->populated[_MOLE_MIXING_RATIO]=TRUE;
	}
	if(WX->populated[_ABSOLUTE_HUMIDITY]==FALSE){
		for(i=0;i<WX->N;i++)
			WX->val[_ABSOLUTE_HUMIDITY][i] = calcAbsolute(WX->val[_TEMPERATURE_K][i],WX->val[_PRESSURE][i],
															WX->val[_MOLE_MIXING_RATIO][i]);
		WX->populated[_ABSOLUTE_HUMIDITY]=TRUE;
	}
	if(WX->populated[_MOIST_AIR_DENSITY]==FALSE){
		for(i=0;i<WX->N;i++)
			WX->val[_MOIST_AIR_DENSITY][i] = calcMoistAirDensity(WX->val[_TEMPERATURE_K][i], WX->val[_PRESSURE][i],
												WX->val[_MOLE_MIXING_RATIO][i], WX->xCO2);
		WX->populated[_MOIST_AIR_DENSITY]=TRUE;
	}
	if(WX->populated[_MASS_MIXING_RATIO]==FALSE){
		for(i=0;i<WX->N;i++)
				WX->val[_MASS_MIXING_RATIO][i] = WX->val[_RELATIVE_HUMIDITY][i]*WX->val[_SATURATION_MIXING_RATIO][i]/100.0;
		WX->populated[_MASS_MIXING_RATIO]=TRUE;
	}
	if(WX->populated[_VIRTUAL_TEMPERATURE]==FALSE){
		for(i=0;i<WX->N;i++)
			WX->val[_VIRTUAL_TEMPERATURE][i] = calcVirtualTemperature(WX->val[_TEMPERATURE_K][i],
														WX->val[_MASS_MIXING_RATIO][i],FOREWARD);
		WX->populated[_VIRTUAL_TEMPERATURE]=TRUE;
	}
	if(WX->populated[_SPECIFIC_HUMIDITY]==FALSE){
		for(i=0;i<WX->N;i++)
			WX->val[_SPECIFIC_HUMIDITY][i] = calcSpecificHumidity(WX->val[_MASS_MIXING_RATIO][i]) ;
					/* TODO This section was included at some point, but I don't know why.
					 * If it does'n appear to mess up anything else, we can remove. LRB Oct 25, 2019.
					 * 1000.0*WX->val[_ABSOLUTE_HUMIDITY][i]/WX->val[_MOIST_AIR_DENSITY][i];*/
		WX->populated[_SPECIFIC_HUMIDITY]=TRUE;
	}
	if(WX->populated[_VIRTUAL_POTENTIAL_TEMPERATURE]==FALSE){
		for(i=0;i<WX->N;i++)
			WX->val[_VIRTUAL_POTENTIAL_TEMPERATURE][i] = calcVirtualTemperature(WX->val[_POTENTIAL_TEMPERATURE][i],
															WX->val[_MASS_MIXING_RATIO][i],FOREWARD);
		WX->populated[_VIRTUAL_POTENTIAL_TEMPERATURE]=TRUE;
	}
	if(WX->populated[_POTENTIAL_VAPOR_PRESSURE]==FALSE){
		for(i=0;i<WX->N;i++)
			WX->val[_POTENTIAL_VAPOR_PRESSURE][i] = WX->val[_VAPOR_PRESSURE][i]*WX->standardPressure/WX->val[_PRESSURE][i];
		WX->populated[_POTENTIAL_VAPOR_PRESSURE]=TRUE;
	}
	if (WX->populated[_WATER_VAPOR_NUMBER_DENSITY] == FALSE) {
		for (i = 0; i<WX->N; i++)
			WX->val[_WATER_VAPOR_NUMBER_DENSITY][i] = WX->val[_ABSOLUTE_HUMIDITY][i] / WATER_MOLAR_MASS;
		WX->populated[_WATER_VAPOR_NUMBER_DENSITY] = TRUE;
	}
	if (WX->populated[_MOIST_AIR_NUMBER_DENSITY] == FALSE) {
		for (i = 0; i<WX->N; i++)
			WX->val[_MOIST_AIR_NUMBER_DENSITY][i] = WX->val[_WATER_VAPOR_NUMBER_DENSITY][i] / WX->val[_MOLE_MIXING_RATIO][i];
		WX->populated[_MOIST_AIR_NUMBER_DENSITY] = TRUE;
	}
	return WEATHER_CONVERSION_SUCCESS;
}
double standardAtmosAltitudeAtPressure(double P){
	/* NOAA 1976 Standard atmosphere table from Jacobson Fundamentals of
	 * atmospheric modeling */
static double alt[] = {
				  0.0,   0.1,   0.2,   0.3,   0.4,   0.5,   0.6,   0.7,   0.8,   0.9,
				  1.0,   1.5,   2.0,   2.5,   3.0,   3.5,   4.0,   4.5,   5.0,   5.5,
				  6.0,   6.5,   7.0,   7.5,   8.0,   8.5,   9.0,   9.5, 10.0, 11.0,
				12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 21.0,
				22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0, 30.0, 31.0,
				32.0, 33.0, 34.0, 35.0, 36.0, 37.0, 38.0, 39.0, 40.0, 41.0,
				42.0, 43.0, 44.0, 45.0, 46.0, 47.0, 48.0, 49.0, 50.0, 55.0,
				60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0, 95.0, 100.0};
static double prs[] = {
				1013.25, 1001.20, 989.45, 977.72, 966.11, 954.61, 943.22,   931.94,  920.77, 909.71,
				  898.80,   845.59,   795.0,   746.9,   701.2,   657.8,   616.6,     577.5,    540.5,   505.4,
				    472.2,     440.7,   411.1,   383.0,   356.5,   331.5,   308.0,     285.8,    265.0, 227.0,
				    194.0,     165.8,   141.7,   121.1,   103.5,     88.5,     75.7,       64.7,      55.3, 47.3,
				      40.5,       34.7,     29.7,    25.5,      21.9,     18.8,     16.2,       13.9,      12.0, 10.3,
				      8.89,       7.67,     6.63,    5.75,      4.99,     4.33,     3.77,       3.29,      2.87, 2.51,
				      2.20,       1.93,     1.69,    1.49,      1.31,     1.16,     1.02,     0.903,    0.798, 0.425,
				    0.220,     0.109, 0.0522,0.0239,  0.0105, 0.0045, 0.0018, 0.00076,0.00032};
int N = sizeof(prs)/sizeof(double),i;
double gamma;

/* Find the interval that we are in. */
for(i=1;i<N-1;i++){
	if(P>prs[i]) break;
}

/* Get the scaling term */
gamma = log(P/prs[i-1])/log(prs[i]/prs[i-1]);

/* Return the value in meters instead of km. */
return 1000.0*(alt[i]*gamma + alt[i-1]*(1.0 - gamma));
}
double standardAtmosPressureAtAltitude(double Z) {
	/* NOAA 1976 Standard atmosphere table from Jacobson Fundamentals of
	* atmospheric modeling */
	static double alt[] = {
		0.0,   0.1,   0.2,   0.3,   0.4,   0.5,   0.6,   0.7,   0.8,   0.9,
		1.0,   1.5,   2.0,   2.5,   3.0,   3.5,   4.0,   4.5,   5.0,   5.5,
		6.0,   6.5,   7.0,   7.5,   8.0,   8.5,   9.0,   9.5, 10.0, 11.0,
		12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 21.0,
		22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0, 30.0, 31.0,
		32.0, 33.0, 34.0, 35.0, 36.0, 37.0, 38.0, 39.0, 40.0, 41.0,
		42.0, 43.0, 44.0, 45.0, 46.0, 47.0, 48.0, 49.0, 50.0, 55.0,
		60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0, 95.0, 100.0 };
	static double prs[] = {
		1013.25, 1001.20, 989.45, 977.72, 966.11, 954.61, 943.22,   931.94,  920.77, 909.71,
		898.80,   845.59,   795.0,   746.9,   701.2,   657.8,   616.6,     577.5,    540.5,   505.4,
		472.2,     440.7,   411.1,   383.0,   356.5,   331.5,   308.0,     285.8,    265.0, 227.0,
		194.0,     165.8,   141.7,   121.1,   103.5,     88.5,     75.7,       64.7,      55.3, 47.3,
		40.5,       34.7,     29.7,    25.5,      21.9,     18.8,     16.2,       13.9,      12.0, 10.3,
		8.89,       7.67,     6.63,    5.75,      4.99,     4.33,     3.77,       3.29,      2.87, 2.51,
		2.20,       1.93,     1.69,    1.49,      1.31,     1.16,     1.02,     0.903,    0.798, 0.425,
		0.220,     0.109, 0.0522,0.0239,  0.0105, 0.0045, 0.0018, 0.00076,0.00032 };
	int N = sizeof(prs) / sizeof(double), i;
	double gamma;
	Z /= 1000.0; /* Convert to km*/

	/* Find the interval that we are in. */
	for (i = 1; i<N - 1; i++) {
		if (Z<alt[i]) break;
	}

	/* Get the scaling term */
	gamma = log(prs[i] / prs[i - 1]) / (alt[i] - alt[i - 1]);

	/* Return the value in meters instead of km. */
	return prs[i-1]*exp(gamma*(Z - alt[i-1]));
}
double _relError(double a, double b){
	if((a==0.0)&&(b==0.0)) return 0.0;
	return pow((a-b) / (fabs(a)>fabs(b)?a:b),2.0);
};
