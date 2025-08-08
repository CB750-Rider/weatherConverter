#include <stdio.h>
#include <stdlib.h>
#include "weatherConversion.h"

int main(int argc, char *argv[]){
    char *cp;
    double P = strtod(argv[1], &cp);
    double vp = strtod(argv[2], &cp);
    double mr = calcMassMixingRatio(P, vp);

    printf("Starting with P=%lf, vp=%lf, we get mr=%lf.\n", P, vp, mr);

    vp = calcVaporPressureFromMassMixingRatio(P ,mr);

    printf("Starting with P=%lf, mr=%lf, we get vp=%lf.\n", P, mr, vp);

    return 0;
}