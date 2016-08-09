CC            = cc
CFLAGS        = -g
OUTPUT_OPTION = -MMD -MP

SOURCE = $(wildcard *10.c)
OBJS   = $(SOURCE:.c=.o)
DEPS   = $(SOURCE:.c=.d)
-include ${DEPS}

ex10: ${OBJS}
	${CC} -o ex10 ${OBJS}

etags:
	etags *.c *.h
