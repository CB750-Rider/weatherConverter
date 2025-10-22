CFLAGS += -Wall -fmessage-length=0 -fPIC
LFLAGS += -lm
SRCS = conversionVectorHelpers.c weatherConversion.c US1976_Standard_Atmos_Table_8.c
TST_SRC = converterTest.c
DEPS = weatherConversion.h
OBJS = $(SRCS:.c=.o)
TST_OBJ = $(TST_SRC:.c=.o)

all: release 

debug: CFLAGS += -g3 -O0 -D__BASE_FILE_NAME__=\"$(notdir $<)\"
release: CFLAGS += -O3 -D__BASE_FILE_NAME__=\"$(notdir $<)\"

.PHONY= test clean all libs

debug: libweatherconverter.so libweatherconverter.a testWeatherConversion

release: libweatherconverter.so libweatherconverter.a testWeatherConversion

%.o: %.c $(DEPS)
	$(CC) -c $(CFLAGS) -o $@ $<

libweatherconverter.so: $(OBJS)
	$(CC) -shared -o $@ $^

libweatherconverter.a: $(OBJS)
	ar r $@ $^
	ranlib $@

testWeatherConversion: $(OBJS) $(TST_OBJ)
	$(CC) -o $@ $^ $(LFLAGS)
	@echo "Run testWeatherConversion to check and see that everything is working as expected."

clean:
	rm -f $(OBJS)
