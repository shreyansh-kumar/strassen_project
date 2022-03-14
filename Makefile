CC = gcc-11
LD = gcc-11

MATRIX_SOURCE  =multiplication.c
MATRIX_SOURCE +=main.c
MATRIX_SOURCE +=matrix.c

MATRIX_OBJECTS = $(addsuffix .o,$(basename $(MATRIX_SOURCE)))

MATRIX_TARGET = MatrixMultiplication

MATRIX_CFLAGS  =-Wall -std=c17 -g -O0 $(CFLAGS)#-DDEBUG_TRANSFORMATIONS #-fsanitize=unsigned-integer-overflow -fno-sanitize-recover=unsigned-integer-overflow 
MATRIX_LDFLAGS =-g -O0 $(LDFLAGS)#-fsanitize=unsigned-integer-overflow -fno-sanitize-recover=unsigned-integer-overflow 

all : MatrixImplementation

clean:
	rm -r *.o MatrixImplementation

%.o : %.c
	${CC} ${MATRIX_CFLAGS} -c $<

MatrixImplementation : $(MATRIX_OBJECTS)
	$(LD) $(MATRIX_LDFLAGS) -o $@ $(MATRIX_OBJECTS)