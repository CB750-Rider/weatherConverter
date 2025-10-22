#include "weatherConversion.h"
#include "table_8.h"
#include <math.h>
#include <string.h>

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

