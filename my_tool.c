#include <stdio.h>
#include "table_8.h"

int main(int argc, char *argv[]){
    int i;
    for(i=0;H[i]<1e-13;i++)
        printf("Pressure[%d]:H[%d] = %lg : %lg.\n", i, i, presures[i], H[i]);
    printf("Pressure[%d]:H[%d] = %lg : %lg.\n", i, i, presures[i], H[i]);

    return 0;
}