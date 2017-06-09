/*
 * weatherConversion.c
 *
 *  Created on: Aug 11, 2016
 *      Author: lee
 */

#include "weatherConversion.h"
#ifdef _MSC_BUILD
#define _USE_MATH_DEFINES
#include <math.h>
#define _MY_INFINITY HUGE_VAL
#else
#define _MY_INFINITY 1.0/0.0
#endif
#include <stdlib.h>
#include <stdio.h>

char *_weather_converter_field_names[_N_WEATHER_FIELDS] = {
													"Temperature C",
													"Temperature K",
													"Temperature F",
													"U wind",
													"V wind",
													"Wind Speed",
													"Wind Direction",
													"Pressure",
													"Potential Temperature",
													"Virtual Temperature",
													"Virtual Potential Temperature",
													"Saturation Vapor Pressure",
													"Saturation Mass Mixing Ratio",
													"Enhancement Factor",
													"Relative Humidity",
													"Vapor Pressure",
													"Potential Vapor Pressure",
													"Mole Mixing Ratio",
													"Mass Mixing Ratio",
													"Dew Point",
													"Specific Humidity",
													"Absolute Humidity",
													"Moist Air Density",
													"Other Input Field"};

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
	"_DEW_POINT", /* 19 Kelvin */
	"_SPECIFIC_HUMIDITY", /* 20 grams water vapor / kilogram moist air */
	"_ABSOLUTE_HUMIDITY", /* 21 grams water vapor / meter^3 */
	"_MOIST_AIR_DENSITY", /* 22 grams / meter^3 */ 
    "_OTHER_INPUT" /* For any other field */};

char *_weather_converter_field_units_full[_N_WEATHER_FIELDS] = {
													"degrees Celsius",
													"Kelvin",
													"degrees Fahrenheit",
													"meters / second",
													"meters / second",
													"meters / second",
													"degrees from North",
													"millibar",
													"Kelvin",
													"Kelvin",
													"Kelvin",
													"millibar",
													"grams water vapor / kilogram moist air",
													"unitless",
													"percent",
													"millibar",
													"millibar",
													"moles water vapor / mole moist air",
													"grams water vapor / kilogram moist air",
													"Kelvin",
													"grams water vapor / kilogram moist air",
													"grams water vapor / meter^3",
													"grams / meter^3",
													"undefined"};

char *_weather_converter_field_units[_N_WEATHER_FIELDS] = {
													"°C",
													"K",
													"°F",
													"m/s",
													"m/s",
													"m/s",
													"°",
													"mb",
													"K",
													"K",
													"K",
													"mb",
													"g/kg",
													" ",
													"%%",
													"mb",
													"mb",
													"mol/mol",
													"g/kg",
													"K",
													"g/kg",
													"g/m^3",
													"g/m^3",
													"----"};

double m_sToKnots(double ms){
	return ms*1.9438444924574;}

double knotsTom_s(double kt){
	return kt* 0.51444444444;}

double KtoF(double K){
	return 9.0/5.0*(K-273.15)+32;}

double FtoK(double F){
	return (5.0/9.0*(F-32.0) + 273.15);}

double CtoF(double C){
	return (9.0/5.0*C + 32.0);}

double FtoC(double F){
	return 5.0/9.0*(F-32.0);}

double calcDewpoint(double T, double RH, double SVP){
	/* Based on the method presented by NOAA,
	 * http://www.srh.noaa.gov/images/epz/wxcalc/wetBulbTdFromRh.pdf
	 *  finds the dewpoint in Kelvin using:
	 *  T	the temperature (K),
	 *  RH	relative humidity (%), and
	 *  SVP	saturation vapor pressure (mb)
	 *  Since this function only approximately inverts the Goff-Gratch equation,
	 *  use Newton's method to improve the inversion to machine precision. */
	double term1, dewPoint, alpha;
	int i;

	term1 = log(SVP*RH/611.2);
	dewPoint = 243.75*term1/(17.67-term1)+273.15;
	/* Now improve the estimate of the dew point using
	 * f(dewpoint) = es*RH/100 - goffGratch(dewPoint)
	 * and find where f=0. */
	alpha = SVP*RH/100;
	for(i=0;i<4;i++) /* We start off close so 4 iterations of Newton's method is plenty. */
		dewPoint += (alpha-goffGratch(dewPoint))/DgoffGratch_dT(dewPoint);

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
    %   See NIST http://www.nist.gov/calibrations/upload/metv29i1p67-2.pdf*/
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

double calcMixingRatio(double P, double vp){
	/*     %CALCMIXINGRATIO Calculate mixing ratio
    % Calculate mixing ratio for a specified pressure and temperature. For
    % saturated mixing ratios use temperature. For standard mixing ratios
    % use dewpoint.
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
	for(i=0;i<6;i++){ /* Six iteratinos of newton's method */
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
	double RHp,RHm,del = 1.0,dS_dR;

	for(step=0;step<6;step++){
		/* use centered differencing to get the slope */
		RHm = (RH - del) >= 0.0 ? RH-del : RH;
		RHp = (RH + del) <= 1.0 ? RH+del : RH;

		/* Find the dSH / dRH. */
		dS_dR = (calcSHfromRH(T,P,xCO2,SVP,RHp,ef) - calcSHfromRH(T,P,xCO2,SVP,RHm,ef) ) /
				(RHp - RHm);

		/* Get the update in RH */
		del = (calcSHfromRH(T,P,xCO2,SVP,RH,ef) - SH)/dS_dR;

		/* update RH */
		RH -= del;

		/* Decide to continue */
		if (fabs(del / RH) < 1e-12) break;
	}

	return RH;

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
							100.0*WX->val[_MASS_MIXING_RATIO][i]*WX->val[_PRESSURE][i]/
							(622.0*WX->val[_SATURATION_VAPOR_PRESSURE][i]+
							WX->val[_MASS_MIXING_RATIO][i]*WX->val[_SATURATION_VAPOR_PRESSURE][i]);
		}
		WX->populated[_RELATIVE_HUMIDITY] = TRUE;
		break;
	case _DEW_POINT:
		if(WX->populated[_VAPOR_PRESSURE]==FALSE)
			for(i=0;i<WX->N;i++)
				WX->val[_VAPOR_PRESSURE][i] = goffGratch(WX->val[_DEW_POINT][i] );
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
				/* TODO The MMR may blow up for SH => 1.0, so find a stable way to get RH. */
				if (WX->val[_SPECIFIC_HUMIDITY][i] > 0.99)
					WX->val[_RELATIVE_HUMIDITY][i] = estRHfromSH(
							WX->val[_SPECIFIC_HUMIDITY][i],
							WX->val[_TEMPERATURE_K][i],
							WX->val[_PRESSURE][i],
							WX->xCO2,
							WX->val[_SATURATION_VAPOR_PRESSURE][i],
							WX->val[_ENHANCEMENT_FACTOR][i],
							100);
				else
					WX->val[_RELATIVE_HUMIDITY][i] =
						100.0*WX->val[_MASS_MIXING_RATIO][i]*WX->val[_PRESSURE][i]/
						(622.0 + WX->val[_MASS_MIXING_RATIO][i])/WX->val[_SATURATION_VAPOR_PRESSURE][i];

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
		printf("weatherConversion.h:humidityConversion() Conversion from Moist air density NOT SUPPORTED YET!! But it is possible.\n");
		for(i=0;i<WX->N;i++){
			}
		//WX->populated[_RELATIVE_HUMIDITY] = TRUE;
		return NO_HUMIDITY_PRESENT;
		break;
	default:
		printf("weatherConversion.h:setAllFields() unable to set humidity fields.\nAborting conversions.\n");
		return NO_HUMIDITY_PRESENT;
		break;
	}
	/* Starting from relative Humidity, perform all conversions */
	humidityConversion(WX);

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
			/* If the actual pressure is less than the saturation pressure, then assume that there is no dry air. */
			if(WX->val[_PRESSURE][i] < WX->val[_SATURATION_VAPOR_PRESSURE][i])
				WX->val[_SATURATION_MIXING_RATIO][i] = _MY_INFINITY;
			else
				WX->val[_SATURATION_MIXING_RATIO][i] = 622.0*WX->val[_SATURATION_VAPOR_PRESSURE][i]/(WX->val[_PRESSURE][i]-WX->val[_SATURATION_VAPOR_PRESSURE][i]);
		}
	}

	WX->populated[_ENHANCEMENT_FACTOR]=TRUE;
	WX->populated[_SATURATION_MIXING_RATIO]=TRUE;
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
	if(WX->populated[_DEW_POINT]==FALSE){
		for(i=0;i<WX->N;i++)
			WX->val[_DEW_POINT][i] = calcDewpoint(WX->val[_TEMPERATURE_K][i], WX->val[_RELATIVE_HUMIDITY][i],
					WX->val[_SATURATION_VAPOR_PRESSURE][i]);
		WX->populated[_DEW_POINT]=TRUE;
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
	if(WX->populated[_VIRTUAL_TEMPERATURE]==FALSE){
		for(i=0;i<WX->N;i++)
			WX->val[_VIRTUAL_TEMPERATURE][i] = calcVirtualTemperature(WX->val[_TEMPERATURE_K][i],
														WX->val[_MASS_MIXING_RATIO][i],FOREWARD);
		WX->populated[_VIRTUAL_TEMPERATURE]=TRUE;
	}
	if(WX->populated[_MOIST_AIR_DENSITY]==FALSE){
		for(i=0;i<WX->N;i++)
			WX->val[_MOIST_AIR_DENSITY][i] = calcMoistAirDensity(WX->val[_TEMPERATURE_K][i], WX->val[_PRESSURE][i],
												WX->val[_MOLE_MIXING_RATIO][i], WX->xCO2);
		WX->populated[_MOIST_AIR_DENSITY]=TRUE;
	}
	if(WX->populated[_MASS_MIXING_RATIO]==FALSE){
		for(i=0;i<WX->N;i++){
			if(WX->val[_PRESSURE][i] > WX->val[_SATURATION_MIXING_RATIO][i])
				WX->val[_MASS_MIXING_RATIO][i] = 622.0*WX->val[_VAPOR_PRESSURE][i]/(WX->val[_PRESSURE][i]-WX->val[_VAPOR_PRESSURE][i]);
			else
				WX->val[_MASS_MIXING_RATIO][i] = _MY_INFINITY;
		}

		WX->populated[_MASS_MIXING_RATIO]=TRUE;
	}
	if(WX->populated[_SPECIFIC_HUMIDITY]==FALSE){
		for(i=0;i<WX->N;i++)
			WX->val[_SPECIFIC_HUMIDITY][i] = 1000.0*WX->val[_ABSOLUTE_HUMIDITY][i]/WX->val[_MOIST_AIR_DENSITY][i];
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
	return WEATHER_CONVERSION_SUCCESS;
}

