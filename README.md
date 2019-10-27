Weather Converter
=================

Weather Converter is intended to provide libraries for conversion of common weather quantities with a focus on humidity values. The converterTest.c function gives an example of how to use the library. It is intended to be used for development, to verify that a change has not caused errors in the output. 

To use the library in your program, create a `WEATHER_CONVERSION_VECTOR`, open it up, populated what values you have, and run `setAllFields`. After that you will be able to read all the converted fields. 

```C

void setInputData(WEATHER_CONVERSION_VECTOR *V,uint N){
	uint i;
	double drh = 40.0/(double)N;
	for(i=0;i<N;i++){
		V->val[_TEMPERATURE_K][i] = 295.0;
		V->val[_PRESSURE][i] = 1100; /* in millibars */
		V->val[_RELATIVE_HUMIDITY][i] = 40.0 + (double)i * drh;
	}
	V->populated[_TEMPERATURE_K]=TRUE;
	V->populated[_PRESSURE]=TRUE;
	V->populated[_RELATIVE_HUMIDITY]=TRUE;
}

int main(){
WEATHER_CONVERSION_VECTOR V;
uint i;

openWeatherConversionVector(&V,10);

setInputData(&V,10);

setAllFields(&V);

/* Do something with the converted data...*/

freeWeatherConversionVector(&V);

return 0;
}
```

If you happen to have your data in a .csv format, you can use my wxFileConverter application to perform the conversion. 
