CC       ?= gcc
MPICC    ?= mpicc
MPIRUN   ?= mpirun --oversubscribe
CFLAGS   ?= -O3 -Wall -Wextra -std=c11 -march=native
CPPFLAGS += -Ilib/common -Ilib/seq -Ilib/mpi -MMD -MP
LDLIBS   += -lm

BUILD := build

COMMON_SRC := $(wildcard lib/common/*.c)
COMMON_OBJ := $(COMMON_SRC:%.c=$(BUILD)/%.o)

SEQ_LIB_OBJ := $(BUILD)/lib/seq/apsp_seq.o
SEQ_OBJ     := $(BUILD)/seq/main.o $(SEQ_LIB_OBJ)
FOX_OBJ     := $(BUILD)/fox/main.o $(BUILD)/lib/mpi/fox.o $(BUILD)/lib/mpi/foxgrid.o $(BUILD)/lib/mpi/apsp_mpi.o

TESTS    := test_matrix test_matrix_io test_minplus
TEST_BIN := $(addprefix $(BUILD)/tests/,$(TESTS))

# Process counts used by check-fox
NP ?= 1 4 9
MATRICES_2_TO_THE = 6 7 8 9 10 11 12 13 14
GEN_SEED ?= 42

.PHONY: all test check check-seq check-fox bench time-seq gen examples-bin clean

all: seq fox

fox: $(FOX_OBJ) $(COMMON_OBJ)
	$(MPICC) $(CFLAGS) -o $@ $^ $(LDLIBS)

seq: $(SEQ_OBJ) $(COMMON_OBJ)
	$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS)

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
$(BUILD)/tests/%: tests/%.c tests/test_util.h $(COMMON_OBJ) $(SEQ_LIB_OBJ)
	@mkdir -p $(dir $@)
	$(CC) $(CPPFLAGS) -Itests $(CFLAGS) -o $@ $< $(COMMON_OBJ) $(SEQ_LIB_OBJ)

# Random matrix generator for benchmark inputs (writes rand%u_input%u.bin to the cwd)
$(BUILD)/tests/matgen: tests/matgen.c $(COMMON_OBJ)
	@mkdir -p $(dir $@)
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ $< $(COMMON_OBJ)

# Text -> binary matrix converter (stdin -> stdout)
$(BUILD)/tests/binconv: tests/binconv.c $(COMMON_OBJ)
	@mkdir -p $(dir $@)
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ $< $(COMMON_OBJ)

# Convert the text examples to the binary format next to the originals
examples-bin: $(BUILD)/tests/binconv
	sh ./tests/examples/to_bin.sh

test: $(TEST_BIN)
	@for t in $(TEST_BIN); do ./$$t || exit 1; done

check-seq: seq examples-bin
	./tests/run_tests.sh

check-fox: fox examples-bin
	@for np in $(NP); do \
	  BIN=./fox NP=$$np MPIRUN="$(MPIRUN)" ./tests/run_tests.sh || exit 1; \
	done

check: test check-seq check-fox

bench: seq fox examples-bin
	MPIRUN="$(MPIRUN)" ./bench/run_bench.sh

INPUT ?= tests/examples/input/input6.bin

time-seq: seq
	start=$$(date +%s%N); ./seq "$(INPUT)" /dev/null; end=$$(date +%s%N); \
	awk -v start="$$start" -v end="$$end" 'BEGIN { printf "elapsed_ms=%.3f\n", (end - start) / 1000000 }'

gen: $(BUILD)/tests/matgen
	@for power in $(MATRICES_2_TO_THE); do \
	  n=$$(awk -v p="$$power" 'BEGIN { print 2 ^ p }'); \
	  echo "gen: N = $$n (2^$$power)"; \
	  $(BUILD)/tests/matgen $$n 1 $(GEN_SEED) || exit 1; \
	done
	@echo "gen: zipping rand0_input*.bin into rand0_input.zip"; \
	zip rand0_input.zip rand0_input*.bin || exit 1

clean:
	rm -rf $(BUILD) seq fox

-include $(shell find $(BUILD) -name '*.d' 2>/dev/null)
