#include "weatherConversion.h"
#include "table_8.h"
#include <math.h>

#define len(A) sizeof(A)/sizeof(A[0])

int binary_find(double *X, double x, int N){
    if (N==1) return 0;
    else if(x > X[N/2])  // Left 1/2
        return binary_find(X, x, N/2);
    else // Right 1/2
        return N/2 + binary_find(X+N/2, x, N-N/2);
}

double tree_search(double *Y, double x) {
    int N = len(pressures);
    int idx = binary_find(pressures, x, N);
    idx = idx > N - 2 ? N -2 : idx;
    double dp = pressures[idx+1] - pressures[idx];
    double dy = log10(Y[idx+1]) - log10(Y[idx]);
    double p = x - pressures[idx];
    return pow(10.0, log10(Y[idx]) + p*dy/dy);
}

double v_tree_search(double *Y, double x) {
    int N = len(pressures);
    int idx = binary_find(pressures, x, N);
    idx = idx > N - 2 ? N -2 : idx;
    double dp = pressures[idx+1] - pressures[idx];
    double dy = log10(Y[idx+1]) - log10(Y[idx]);
    double p = x - pressures[idx];
    return pow(10.0, log10(Y[idx]) + p*dy/dy);
}

double moles_N2_P(double P) { return tree_search(N2, P); }


double moles_O_P(double P) { return tree_search(O, P); }


double moles_O2_P(double P) { return tree_search(O2, P); }


double moles_Ar_P(double P) { return tree_search(Ar, P); }


double moles_He_P(double P) { return tree_search(He, P); }


double moles_H_P(double P) { 
    if (P < 4.35307e-6)
        return tree_search(H, P);
    else
        return 0.0;
}


// conversionVectorHelpers.c 
// @ 93

WEATHER_CONVERTER_FIELD getEnumFromFlag(const char *flag) {
    WEATHER_CONVERTER_FIELD i;
    for(i=0;i<_N_WEATHER_FIELDS;i++){
        if (strcmp(flag, _weather_converter_field_flags[i]) == 0)
            return i;
    }
    return -1;
}
const char *getFieldFlag(WEATHER_CONVERTER_FIELD wcf){
    if((wcf < 0) || (wcf >= _N_WEATHER_FIELDS))
        return "<UNKNOWN FIELD ENUM>";
    return _weather_converter_field_flags[wcf];
}
const char *getFieldName(WEATHER_CONVERTER_FIELD wcf){
    if((wcf < 0) || (wcf >= _N_WEATHER_FIELDS))
        return "<UNKNOWN FIELD ENUM>";
    return _weather_converter_field_names[wcf];
}
const char *getFieldUnitsAbbr(WEATHER_CONVERTER_FIELD wcf){
    if((wcf < 0) || (wcf >= _N_WEATHER_FIELDS))
        return "<UNKNOWN FIELD ENUM>";
    return _weather_converter_field_units[wcf];
}
int get_N_WEATHER_FIELDS() {
    return _N_WEATHER_FIELDS;
}

// @122, remove blank lines

// @134
WX->latitude = ...
WX->surfaceTemperature = _weather_converter_site_defualts[_SURFACE_TEMPERATURE];
for(fi=0;...

...}
WX->quiet = FALSE;
return WEATHER_CONVERSION_SUCCESS; ...
...CONVERTER_FIELD fi;
// printf("\n  *** Freeqing weather converter at %p\n", WX);
for(fi=0...
if(WX->allocated[fi])
    free(WX->val[fi]);

// @ 175 - 187 Block comment WX->f.alloc = (WEATH... ...WX->f.integrate_column_moist...NumberDensity;
// TODO Bring this back or take it out */
if(rv==DOUBLE_IN....

//@260 Spelling fix in "Settings"
WEATHER_CONVERSION_ERROR parseSiteSettingsLine(char *line, WEATHER_CONVERSION_VECTOR *WX)

// @346 remove "/* End function */"
{...
}
void setQuiet(WEATHER_CONVERSION_VECTOR *WX){
    WX->quiet = TRUE;
}
void setNotQuiet(WEATHER_CONVERSION_VECTOR *WX){
    WX->quiet = FALSE;
}
void printWeatherConversionVectorMetadata(WEATHER_CONVERSION_VECTOR *WX){
    int i;
    printf("Weather Conversion Vector Metadata:\n");
    printf(" N: %u\n", WX->N);
    printf(" allocated: ");
    for(i=0;i<_N_WEATHER_FIELDS;i++) printf(", %u", WX->allocated[i]);
    printf(" populated: ");
    for(i=0;i<_N_WEATHER_FIELDS;i++) printf(", %u", WX->populated[i]);
    printf("\n standandPressure: %lf mb\n", WX->standandPressure);
    printf(" carbon dioxide concentration: %lf ppm\n", WX->xCO2);
    prnitf(" latitude: %lf degrees North\n", WX->latitude);
    printf(" surfaceHeight: %lf m MSL\n", WX->surfaceHeight);
    printf(" surfacePressure: %lf mb\n", WX->surfacePressure);
    printf(" surfaceTemperature: %lf K\n", WX->surfaceTemperature);
}/* End function */

/* converterTest.c */
//@34
int importFile(const char *fname, WEATHER_CONVERSION_VECTOR *OUT);

//@69
openWeatherConversionVector(...
setQuiet(&STDARDAND);
setNotQuiet(&STANDARD);

// @78-80
if (importFile("humidityTest.csv", &STANDARD)){
    openWeatherConversionVector(&TEST, STANDARD.N);
    for(ri=_RELATIVE_HUMIDITY;ri<_MOIST_AIR_DENSITY;ri++) {
        printf("Now performing the conversion starting with temperature, pressure, and %s.\n",_weather_converter_field_names[ri]);
        setTestVector(&TEST,&STANDARD,ri);
        printf("MEAN SQUARED RELATIVE ERROR = %g.\n", compare(&TEST,&STANDARD));
    }
    freeWeatherConversionVector(&STANDARD);
    freeWeatherConversionVector(&TEST);
}
else{
    printf("Error opening the stdandard data file 'humidityTest.csv'.\n");
    printf("Skipping the standard data comparison.\n");
}


//@149
    ...fclose(fp);
}
int importFile(const char *fname, WEATHER_CONVERSION_VECTOR *OUT) {
    FILE *fp;
    fp = fopen(fname,"r");
    if (fp==NULL) return 0;
    ...
    fclose(fp);
    return 1;
}
//deleted line

/* weatherConversion.c */

//@50 (note we dropped units from field names)

const char *_weather_converter_field_names[_N_WEATHER_FIELDS] = {
    "Temperature", /*0*/
    "Temperature", /*1*/
    "Temperature", /*2*/
    ...
    "Dry Air Number Density", /*30*/
    "Dry Air Density", /*31*/
    "Other Input Field" /*32*/};

const char *_weather_converter_field_flags[_N_WEATHER_FIELDS] = {
    ...
    "_DRY_AIR_NUMBER_DENSITY", /* 30 */
    "_DRY_AIR_DENSITY", /* 31 */
    "_OTHER_INPUT", /* 32 */};

const char *_weather_converter_field_units_full[_N_WEATHER_FIELDS] = {
    ...
    "moles / meter^3", /* 30 */
    "grams / meter^3", /* 31 */
    "undefined", /* 32 */};

const char *_weather_converter_field_units[_N_WEATHER_FIELDS] = {
    ...
    "%%", /*14*/
    ...
    "mol/m^3", /*30*/
    "g/m^3", /*31*/
    "----" /*32*/};

char *_weather_converter_site_setting...

    "Surface Pressure",
    "Surface Temperature"
};

char *_weather_converter_site_setting_flags[_N_WEATHER_SITE_SPECIFIC_SETTINGS]...
    "_SURFACE_PRESSURE",
    "_SURFACE_TEMPERATURE"
};

char *weather_converter_site_units_full...
    ...
    "Kelvin"
};

char *weather_coverter_site_units[...
    ...
    "mb",
    "K"
};

double _weather_converter_site_defaults[...
    1013.25, /* Surface pressure */
    273.15, /* Surface temperature */
};

//@410
beta = P*100/Z/T/GAS_CONSTANT

//@458
    if(P < 3.7338e-3) return calcDryAirDensity(T, P, xCO2);

    MA = ....

// @467 Consider using GAS_CONSTANT?

if(P < 3.7338e-3) return dryAirNumberDensity(T, P);

Z = calcCompressibility(T,P,xv); 

return ((100*P)/Z*T*GAS_CONSTANT);
}
double calcDryAirDensity(double T, double P, double xCO2){
    /* Calculate the dry air density using the method from Ciddor's paper. Except above
     * the US1976 86 km pressure level, where the number densities change.
     * Inputs:
     * T    Temperature (K)
     * P    Pressure (mb)
     * xCO2 Mole mixing ratio of CO2 to dry air (ppm)  */
     double MA, R=8.314510; /* Mass of dry air, compressibility of moist air, gas constant */
     double molar_mass[] = {28.0134, 15.999, 15.999*2.0, 39.948, 4.002602, 1.008};

     if (P > 3.7338e-3){
        MA = 28.9635 + 12.011e-6*(xCO2 - 400); /* Molar mass of dry air */
        return (MA*(100*P)/T/R/AVAGADRO_CONSTANT)
     }
     else{
        return (molar_mass[0] * moles_N2_P(P)
              + molar_mass[1] * moles_O_P(P)
              + molar_mass[2] * moles_O2_P(P)
              + molar_mass[3] * moles_Ar_P(P)
              + molar_mass[4] * moles_He_P(P)
              + molar_mass[5] * moles_H_P(P));
     }
}

double dryAirNumberDensity(double T, double P) {
    /* Calculate the dry air density using the method from Ciddor's paper. Except above
     * the US1976 86 km pressure level, where the number densities change.
     * Inputs:
     * T    Temperature (K)
     * P    Pressure (mb)*/
     double R=8.314510; /* gas constant */

     if (P > 3.7338e-3){
        return ((100*P)/T/R);
     }
     else
        return moles_N2_P(P)
             + moles_O_P(P)
             + moles_O2_P(P)
             + moles_Ar_P(P)
             + moles_He_P(P)
             + moles_H_P(P);
}
double calcPressure(double T, double rho){
    return rho*GAS_CONSTANT*T;
}

// @533 not sure why we have mr = mr;
/* Pressure should always be greater than vapor pressure for this to work. So, we force the matter. */
// P = P > 1.00001*zp ? P : vp*1.00001;

return 1000.0*epsilon*vp/(P-vp);
}
double calcVaporPressureFromMassMixingRatio(double P, double mr){
    /*  %CALCVAPORPRESSUREFROMMASSMIXINGRATIO Calculate vapor pressure from mass mixing ratio. 
    
    TODO Fill this in, reference above. */
    if (isinf(mr)) {return P;}
        double epsilon = 0.62197;

    mr = mr; 

    return mr*P/(epsilon*1000.0 + mr)
}
double calcMMRfromAbsoluteHumid...
    double R = GAS_CONSTANT; /* Gas constant */
    ...
    //int in_roi = ((P>0.01) && (P<50));
    ...
    // if (in_roi) printf("Convergence @ P=%lg mk: %lg, ",P, xv);
    for...
        num = a - moistAirNumberDensity(T,P,xv)*xv*Mv;
        ...
        //if (in_roi) printf("%lg, ",xv);
    }
    //if (in_roi)printf("\n");
    return xv;
}

//@778
...
WEATHER_CONVERTER_FIELD fi;
if(WX->populated[_ABSOLUTE_HUMIDITY]) return _ABSOLUTE_HUMIDITY
for(fi=_RELATIVE_HUMIDITY...

// @788
if (WX->quiet != TRUE)
printf(...


// @825
//- case _MASS_MIXING_RATIO
case _ABSOLUTE_HUMIDITY:
    if(WX->populated[_MOLE_MIXING_RATIO]==FALSE){
        for(i=0;i<WX->N;i++)
            WX->val[_MOLE_MIXING_RATIO][i] = calcMMRfromAbsoluteHumidity(WX->val[_PRESSURE][i], WX->val[_TEMPERATURE_K][i],
            WX->val[_ABSOLUTE_HUMIDITY][i]);
        WX->populated[_MOLE_MIXING_RATIO] = TRUE;
    }
    if(WX->populated[_RELATIVE_HUMIDITY]==FALSE){
        for(i=0;i<WX->N;i++)
            WX->val[_RELATIVE_HUMIDITY][i] = 100.0*(WX->val[_MOLE_MIXING_RATIO][i]*WX->val[_PRESSURE][i]) /
        (WX->val[_SATURATION_VAPOR_PRESSURE][i]*WX->val[_ENHANCEMENT_FACTOR][i]);
        WX->populated[_RELATIVE_HUMIDITY] = TRUE;
    }
    break; ...


// @879
...WX->val[_MASS_MIXING_RATIO][i] ...
if(WX->populated[_VAPOR_PRESSURE] == FALSE)
        for(i=0;i<WX->N;i++)
            WX->val[_VAPOR_PRESSURE][i] = goffGratch(WX->val[_DEW_POINT][i]);
if(WX->populated[_RELATIVE_HUMIDITY] == FALSE)
        for(i=0;i<WX->N;i++)
            WX->val[_RELATIVE_HUMIDITY][i] = 100.0*WX->val[_VAPOR_PRESSURE][i] = WX->val[_SATURATION_VAPOR_PRESSURE][i];
}
WX->populated[_MASS_MIXING_RATIO] = TRUE;
WX->populated[_VAPOR_PRESSURE] = TRUE; 
WX->populated[_RELATIVE_HUMIDITY] = TRUE;
break;
case _MASS_MIXING_RATIO:
if(WX->populated[_VAPOR_PRESSURE]==FALSE){
        for(i=0;i<WX->N;i++)
            WX->val[_VAPOR_PRESSURE][i] = calcVaporPressureFromMassMixingRatio(
                WX->val[_PRESSURE][i], WX->val[_MASS_MIXING_RATIO][i]);
}
if(WX->populated[_RELATIVE_HUMIDITY]==FALSE){
        for(i=0;i<WX->N;i++)
            WX->val[_RELATIVE_HUMIDITY][i] = 100.0*WX->val[_VAPOR_PRESSURE][i] = WX->val[_SATURATION_VAPOR_PRESSURE][i];
}
WX->populated[_RELATIVE_HUMIDITY] = TRUE;
WX->populated[_VAPOR_PRESSURE] = TRUE;
break;
case _MOIST_AIR_DENSITY: ...


// @935 This was just cleaning up spacing over the next few lines (< 10)

// @1023
    WX->populated[PRESSURE]=TRUE;
    }
    else if(WX->populated[_MOIST_AIR_DENSITY] && WX->populated[_TEMPERATURE_K]){
        if(WX->populated[_PRESSURE]==FALSE)
            for(i=0;i<WX->N;i++)
                WX->val[_PRESSURE][i] = calcPressure(WX->val[_TEMPERATURE_K][i], WX->val[_MOIST_AIR_DENSITY][i]);
        WX->populated[_PRESSURE] = TRUE;
    }
    else if(WX->populated[_DRY_AIR_DENSITY] && WX->populated[_ABSOLUTE_HUMIDITY] && WX->populated[_TEMPERATURE_K]){
        if(WX->populated[_MOIST_AIR_DENSITY]==FALSE)
            for(i=0;i<WX->N;i++)
                WX->val[_MOIST_AIR_DENSITY][i] = WX->val[_DRY_AIR_DENSITY][i] + WX->val[_ABSOLUTE_HUMIDITY][i];
        if(WX->populated[_PRESSURE]==FALSE)
            for(i=0;i<WX->N;i++)
                WX->val[_PRESSURE][i] = calcPotentialTemperature(WX->standandPressure, WX->val[_TEMPERATURE_K][i],WX->val[_POTENTIAL_TEMPERATURE][i],FOREWARD);
        WX->populated[_MOIST_AIR_DENSITY] = WX->populated[_PRESSURE] = TRUE;

    }
    else    
        return NO_PRESSURE_PRESENT;
    ...
    if(WX->populated[_SATURATION_MIXING_RATIO]==FALSE){
        for(i=0;i<WX->N;i++)
            WX->val[_SATURATION_MIXING_RATIO][i] = calcMassMixingRatio(WX->val[_PRESSURE][i], WX->val[_SATURATION_VAPOR_PRESSURE][i]);
    }...

// @1054
double max_ratio = 0.9; // 10.0 * (1.0 - nextafter(1.0, 0.0));
if(WX->populated[_GEOPOTENTIAL_HEIGHT] && WX->populated[_HEIGHT_AMSL] && WX->populated[_HEIGHT_AGL])
    return WEATHER_CONVERSION_SUCCESS;
if(WX->populated[_GEOPOTENTIAL_HEIGHT]){
    if(WX->populated[_HEIGHT_AMSL]){
        for(i=0;i<WX->N;i++)
            WX->val[_HEIGHT_AGL][i] = WX->val[_HEIGHT_AMSL] - WX->surfaceHeight;
    }
    else{
        if(WX->populated[_HEIGHT_AGL]){
            for(i=0;i<WX->N;i++)
                WX->val[_HEIGHT_AMSL] = WX->val[_HEIGHT_AGL] + WX->surfaceHeight;
        }
        else{
            /* Convert from Geopotential Height to Above Mean Sea Level and Above Ground Level.
            * Theoretically, GPH cannot be tha same as Earth's radius except where the AMSL height is infinite. */
            for(i=0;i<WX->N;i++){
                if(WK->val[_GEOPOTENTIAL_HEIGHT][i] / R > max_ratio ){
                    WX->val[_HEIGHT_AMSL][i] = R*1e16;
                    printf("Ratio issue in %s at %d! Ratio = %g, Z = %g, Max Ratio = %g, Radius  = %g\n". __FILE__, __LINE__, WX->val[_GEOPOTENTIAL_HEIGHT][i] / R, WX->val[_GEOPOTENTIAL_HEIGHT][i], max_ratio, R );
                }
                else
                    WX->val[_HEIGHT_AMSL][i] = WX->val[_GEOPOTENTIAL_HEIGHT][i] / (1.0 - WX->val[_GEOPOTENTIAL_HEIGHT][i]/R);
                WX->val[_HEIGHT_AGL][i] = WX->val[_HEIGHT_AMSL] - WX->surfaceHeight;
            } /* for each height */
        } /* If neither Height MSL or Height AGL are provided */
        WX->populated[_HEIGHT_AMSL]=WX->populated[_HEIGHT_AGL] = TRUE;
    } /* If Hight AMSL is not provided*/
} /* If Geopotential Height is provided*/
...
//printf("%lf ", WX->val[_HEIGHT_AGL][i])


//@1154

if(WX->populated[_MOLE_MIXING_RATIO]==FALSE){
    for(i=0;i<WX->N;i++){
        if (isnan(WX->val[_ENHANCEMENT_FACTOR][i]) || isnan(WX->val[_RELATIVE_HUMIDITY][i])){
            WX->val[_MOLE_MIXING_RATIO][i] = 0.0;
        }
        else {
            WX->val[_MOLE_MIXING_RATIO][i] = WX->val[_ENHANCEMENT_FACTOR][i] * WX->val[_RELATIVE_HUMIDITY][i] * WX->val[_SATURATION_VAPOR_PRESSURE][i] / WX->val[PRESSURE][i] / 100.0;
        }
    }
    WX->populated[_MOLE_MIXING_RATIO] = TRUE;
}

//@1166
if(WX->populated[_MOIST_AIR_DENSITY]==FALSE){
    for(i=0;i<WX->N;i++){
       WX->val[_MOIST_AIR_DENSITY][i] = calcMoistAirDensity(WX->val[_TEMPERATURE_K][i], WX->val[_PRESSURE][i],
            WX->val[_MOLE_MIXING_RATIO][i], WX->val[_XCO2][i]); 
    }  
    WX->populated[_MOIST_AIR_DENSITY]=TRUE;
}
if(WX->populated[_MASS_MIXING_RATIO]==FALSE){
    for(i=0;i<WX->N;i++){
        WX->val[_MASS_MIXING_RATIO][i] = calcMassMixingRatio(WX->val[_PRESSURE][i],
        WX->val[_VAPOR_PRESSURE][i]);    
    }  
    WX->populated[_MASS_MIXING_RATIO]=TRUE;
}

//@1185 Removed TODO comment on 3 lines

//@1206
...
WX->populated[_WATER_VAPOR_NUMBER_DENSITY] = TRUE;
}
if(WX->populated[_DRY_AIR_NUMBER_DENSITY] ==  FALSE) {
    for(i=0;i<WX->N;i++){
        WX->val[_DRY_AIR_NUMBER_DENSITY][i] = dryAirNumberDensity(WX->val[_PRESSURE][i],
        WX->val[_TEMPERATURE_K][i]);    
    }  
    WX->populated[_DRY_AIR_NUMBER_DENSITY] = TRUE;
}
if(WX->populated[_MOIST_AIR_NUMBER_DENSITY]==FALSE) {
    for(i=0;i<WX->N;i++){
        if(isnan(WX->val[_MOLE_MIXING_RATIO][i]) || WX->val[_WATER_VAPOR_NUMBER_DENSITY][i] == 0){
            WX->val[_MOIST_AIR_NUMBER_DENSITY][i] = WX->val[_DRY_AIR_NUMBER_DENSITY][i];
        }    
        else {
            WX->val[_MOIST_AIR_NUMBER_DENSITY][i] = WX->val[_WATER_VAPOR_NUMBER_DENSITY][i] / WX->val[_MOLE_MIXING_RATIO][i];
        }
    }  
    WX->populated[_MOIST_AIR_NUMBER_DENSITY] = TRUE;
}
if(WX->populated[_DRY_AIR_DENSITY]==FALSE){
    for(i=0;i<WX->N;i++){
       WX->val[_DRY_AIR_DENSITY][i] = calcDryAirDensity(WX->val[_PRESSURE][i],
        WX->val[_TEMPERATURE_K][i], WX->val[_XCO2][i]); 
    }  
    WX->populated[_DRY_AIR_DENSITY] = TRUE;
}
...


/* weatherConversion.h */

// @57
#define AVAGADRO_CONSTANT 6.02214076e23
#define GAS_CONSTANT 8.314598

// @72 spelling
/* TODO Add potential temp, virt, temp, virt potential temp in F and C, */

// @107
_DRY_AIR_NUMBER_DENSITY, /* 30 mole / meter^3 */
_DRY_AIR_DENSITY, /* 31 grams / meter^3 */
_OTHER_INPUT, /* 32 Anything else. This should be the enum before _N_WEATHER_FIELDS. */

// @133
_SURFACE_PRESSURE ...
_SURFACE_TEMPERATURE, /* Surface temperature in Kelvin, default is 273.15 */
_N_WEATHER_SITE_SPECIFIC_SETTINGS

// @177
double surfaceTemperature; /* in Kelvin */
BOOLEAN quiet; /* Suprress warning-level messages */
/*WEATHER_CONVERSION_FUNCITONS f; / * Reference available functions. TODO Bring this back or remove it. */

/* Strings defined in Wetaher Conversion.c  */
extern const char *_weather_converter_field_names[_N_WEATHER_FIELDS];
extern const char *_weather_converter_field_flags[_N_WEATHER_FIELDS];
extern const char *_weather_converter_field_units_full[_N_WEATHER_FIELDS];
extern const char *_weather_converter_field_units[_N_WEATHER_FIELDS];

//@214
double moistAirNumberDensity...
double calcDryAirDensity(double T, double P, double xCO2);
double calcPressure(double T, double rho);
double dryAirNumberDensity(double P, double T);
...
double calcMassMixingRatio...
double calcVaporPressureFromMassMixingRatio(double P, double mr);
...

// @245
/* Functions in conversionVectorHelpers.c*/
WEATHER_CONVERTER_FIELD getEnumFromFlag(const char *flag);
const char *getFieldFlag(WEATHER_CONVERTER_FIELD wcf);
const char *getFieldName(WEATHER_CONVERTER_FIELD wcf);
const char *getFieldUnitsAbbr(WEATHER_CONVERTER_FIELD wcf);
int get_N_WEATHER_FIELDS();


//@256
WEATHER_CONVERSION_ERROR changeHeight...
WEATHER_CONVERSION_ERROR parseSiteSettingsLine(char *line, WEATHER_CONVERSION_VECTOR *WX);
void setQuiet(WEATHER_CONVERSION_VECTOR *WX);
void setNotQuiet(WEATHER_CONVERSION_VECTOR *WX);
void printWeatherConversionVectorMetadata(WEATHER_CONVERSION_VECTOR *WX);

/* Functions defined in US1976_Standard_Atmos_Table_8.c */
double moles_N2_P(double P);
double moles_O_P(double P);
double moles_O2_P(double P);
double moles_Ar_P(double P);
double moles_He_P(double P);
double moles_H_P(double P);





