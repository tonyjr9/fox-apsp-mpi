CC      = gcc
CFLAGS  = -O2 -Wall -Wextra -std=c11 -Isrc/common
COMMON  = src/common/matrix.c src/common/matrix_io.c

test_io: tests/test_io.c $(COMMON)
	$(CC) $(CFLAGS) -o $@ $^
