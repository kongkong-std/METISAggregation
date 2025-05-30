# module purge
# module load generators/make/4.3/gcc-10.3.0 mpi/mpich/3.2.1/gcc-10.3.0 compilers/gcc/10.3.0
# module load tpl/gklib/gklib tpl/metis/metis tpl/parmetis/parmetis
CC = mpicc
SRC = ../src

CFLAGS = -c
CMD = -fPIC -Wall -Wwrite-strings \
      -Wno-unknown-pragmas -Wno-lto-type-mismatch \
      -fstack-protector -fvisibility=hidden -g -O3 \
      -march=native -mtune=native

INCL = -I${MPICH_HOME}/include \
	-I${GKLIB_INCLUDE} \
	-I${METIS_INCLUDE} \
	-I${PARMETIS_INCLUDE}


LIB_PATH = -Wl,-rpath,${MPICH_HOME}/lib \
	   -L${MPICH_HOME}/lib \
	   -Wl,-rpath,/opt/apps/compilers/gcc-10.3.0/lib/gcc/x86_64-pc-linux-gnu/10.3.0 \
	   -L/opt/apps/compilers/gcc-10.3.0/lib/gcc/x86_64-pc-linux-gnu/10.3.0 \
	   -Wl,-rpath,/opt/apps/compilers/gcc-10.3.0/lib64 \
	   -L/opt/apps/compilers/gcc-10.3.0/lib64 \
	   -Wl,-rpath,/opt/apps/compilers/gcc-10.3.0/lib \
	   -L/opt/apps/compilers/gcc-10.3.0/lib \
	   -L${METIS_LIB} \
	   -Wl,-rpath,${METIS_LIB} \
	   -L${GKLIB_LIB} \
	   -Wl,-rpath,${GKLIB_LIB} \
	   -L${PARMETIS_LIB} \
	   -Wl,-rpath,${PARMETIS_LIB}

LIB_FLAG = -lpthread -lparmetis -lmetis -lGKlib -lmetis -lGKlib \
	   -lm -lX11 -ldl -lmpifort -lmpi -lgfortran -lm -lgfortran \
	   -lm -lgcc_s -lquadmath -lstdc++ -lquadmath -ldl

LIB = $(LIB_PATH) $(LIB_FLAG)
OBJECT := $(patsubst %.c,%.o,$(wildcard $(SRC)/*.c))

.PHONY: all clean

all: app_parmetis_exe

app_parmetis_exe: $(OBJECT)
	$(CC) $(CMD) $(INCL) -o app_parmetis_exe $(OBJECT) $(LIB)

%.o: %.c
	$(CC) $(CFLAGS) $(CMD) $(INCL) "$<" -o "$@"

clean:
	rm -f app_parmetis_exe $(OBJECT)