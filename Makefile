CC       ?= gcc
MPICC    ?= mpicc
CFLAGS   ?= -O3 -Wall -lm -Wextra -std=c11 -march=native
CPPFLAGS += -Ilib/common -Ilib/seq -Ilib/mpi -MMD -MP

BUILD := build

COMMON_SRC := $(wildcard lib/common/*.c)
COMMON_OBJ := $(COMMON_SRC:%.c=$(BUILD)/%.o)

SEQ_OBJ := $(BUILD)/seq/main.o $(BUILD)/lib/seq/apsp_seq.o
FOX_OBJ := $(BUILD)/fox/main.o $(BUILD)/lib/mpi/fox.o

TESTS    := test_matrix test_matrix_io test_minplus
TEST_BIN := $(addprefix $(BUILD)/tests/,$(TESTS))

.PHONY: all test check check-seq time-seq clean

# The MPI target will be enabled once Fox's implementation is complete.
all: seq

fox: $(FOX_OBJ) $(COMMON_OBJ)
	$(MPICC) $(CFLAGS) -o $@ $^ -lm

seq: $(SEQ_OBJ) $(COMMON_OBJ)
	$(CC) $(CFLAGS) -o $@ $^

# src/main.c is compiled twice: with MPI (-DUSE_MPI) and without.
$(BUILD)/fox/main.o: src/main.c
	@mkdir -p $(dir $@)
	$(MPICC) $(CPPFLAGS) $(CFLAGS) -DUSE_MPI -c $< -o $@

$(BUILD)/seq/main.o: src/main.c
	@mkdir -p $(dir $@)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(BUILD)/lib/mpi/%.o: lib/mpi/%.c
	@mkdir -p $(dir $@)
	$(MPICC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(BUILD)/%.o: %.c
	@mkdir -p $(dir $@)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

# Unit tests (no MPI needed)
$(BUILD)/tests/%: tests/%.c tests/test_util.h $(COMMON_OBJ)
	@mkdir -p $(dir $@)
	$(CC) $(CPPFLAGS) -Itests $(CFLAGS) -o $@ $< $(COMMON_OBJ)

test: $(TEST_BIN)
	@for t in $(TEST_BIN); do ./$$t || exit 1; done

check-seq: seq
	./tests/run_tests.sh

check: test check-seq

INPUT ?= tests/examples/input/input6

time-seq: seq
	start=$$(date +%s%N); ./seq < "$(INPUT)" > /dev/null; end=$$(date +%s%N); \
	awk -v start="$$start" -v end="$$end" 'BEGIN { printf "elapsed_ms=%.3f\n", (end - start) / 1000000 }'

clean:
	rm -rf $(BUILD) seq fox

-include $(shell find $(BUILD) -name '*.d' 2>/dev/null)
