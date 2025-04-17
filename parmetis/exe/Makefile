# module purge
# module load compilers/gcc/10.3.0 generators/make/4.3/gcc-10.3.0 tools/cJSON/1.7.18 tpl/gklib/gklib tpl/metis/metis
CC = gcc
SRC = ../src

CFLAGS = -c
CMD = -fPIC -Wall -Wwrite-strings \
      -Wno-unknown-pragmas -Wno-lto-type-mismatch \
      -fstack-protector -fvisibility=hidden -g -O3 \
      -march=native -mtune=native

INCL = -I${CJSON_INCLUDE_PATH} \
	-I${GKLIB_INCLUDE} \
	-I${METIS_INCLUDE}


LIB_PATH = -L${CJSON_LIB_PATH} \
	   -Wl,-rpath,${CJSON_LIB_PATH} \
	   -L${METIS_LIB} \
	   -Wl,-rpath,${METIS_LIB} \
	   -L${GKLIB_LIB} \
	   -Wl,-rpath,${GKLIB_LIB}

LIB_FLAG = -lcjson -lGKlib -lmetis -lGKlib -lm

LIB = $(LIB_PATH) $(LIB_FLAG)
OBJECT := $(patsubst %.c,%.o,$(wildcard $(SRC)/*.c))

.PHONY: all clean

all: app_metis_exe

app_metis_exe: $(OBJECT)
	$(CC) $(CMD) $(INCL) -o app_metis_exe $(OBJECT) $(LIB)

%.o: %.c
	$(CC) $(CFLAGS) $(CMD) $(INCL) "$<" -o "$@"

clean:
	rm -f app_metis_exe $(OBJECT)