#include "weatherConversion.h"
#include "table_8.h"
#include <math.h>
#include <string.h>

#define len(A) sizeof(A)/sizeof(A[0])

#define MOLAR_MASS_M0 28.9522 // See US1976 Standard Atmosphere, Sec 1.2.4, 3rd paragraph
/* See US1976 Standard Atmosphere Table 3 */
#define MOLAR_MASS_CO2 44.009995
#define MOLAR_MASS_Ne 20.183
#define MOLAR_MASS_Kr 83.8
#define MOLAR_MASS_Xe 131.30
#define MOLAR_MASS_CH4 16.04303
#define MOLAR_MASS_H2 2.01594
#define MOLAR_MASS_N2O 44.0128
#define MOLAR_MASS_SF6 146.055
#define FRAC_VOL_N2 0.78084
#define FRAC_VOL_O2 0.209476
#define FRAC_VOL_Ar 9.34e-3
#define FRAC_VOL_Ne 1.818e-5
#define FRAC_VOL_He 5.24e-6
#define FRAC_VOL_Kr 1.14e-6
#define FRAC_VOL_Xe 8.7e-8
#define FRAC_VOL_H2 5.00e-7
// TODO Make these Dynamic
#define FRAC_VOL_CH4 2.00e-6
#define FRAC_VOL_N2O 3.39e-7
#define FRAC_VOL_SF6 1.225e-11


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
    return pow(10.0, log10(Y[idx]) + p*dy/dp);
}

double v_tree_search(double *Y, double x) {
    int N = len(pressures);
    int idx = binary_find(pressures, x, N);
    idx = idx > N - 2 ? N -2 : idx;
    double dp = pressures[idx+1] - pressures[idx];
    double dy = log10(Y[idx+1]) - log10(Y[idx]);
    double p = x - pressures[idx];
    return pow(10.0, log10(Y[idx]) + p*dy/dp);
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


double molar_density_N2_P(double P) {return MOLAR_MASS_N2 * moles_N2_P(P);}
double molar_density_O_P(double P) {return MOLAR_MASS_O * moles_O_P(P);}
double molar_density_O2_P(double P) {return MOLAR_MASS_O2 * moles_O2_P(P);}
double molar_density_Ar_P(double P) {return MOLAR_MASS_Ar * moles_Ar_P(P);}
double molar_density_He_P(double P) {return MOLAR_MASS_He * moles_He_P(P);}
double molar_density_H_P(double P) {return MOLAR_MASS_H * moles_H_P(P);}

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