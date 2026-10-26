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
FOX_LIB_OBJ := $(BUILD)/lib/mpi/fox.o $(BUILD)/lib/mpi/foxgrid.o
# fox: the assignment's interface (text matrix on stdin, result on stdout).
# fox_bin, fox_shared: binary input and output files, used by tests and benchmarks.
FOX_OBJ        := $(BUILD)/fox/main_text.o $(FOX_LIB_OBJ) $(BUILD)/lib/mpi/apsp_mpi.o
FOX_BIN_OBJ    := $(BUILD)/fox/main.o $(FOX_LIB_OBJ) $(BUILD)/lib/mpi/apsp_mpi.o
FOX_SHARED_OBJ := $(BUILD)/fox/main.o $(FOX_LIB_OBJ) $(BUILD)/fox_shared/apsp_mpi.o

TESTS    := test_matrix test_matrix_io test_minplus
TEST_BIN := $(addprefix $(BUILD)/tests/,$(TESTS))

# Process counts used by check-fox
NP ?= 1 4 9
MATRICES_2_TO_THE = 6 7 8 9 10 11 12 13 14
GEN_SEED ?= 42

.PHONY: all test check check-seq check-fox bench time-seq gen examples-bin report dist clean

all: seq fox_bin fox_shared fox  # fox last, so its note is the last thing printed

fox: $(FOX_OBJ) $(COMMON_OBJ)
	$(MPICC) $(CFLAGS) -o $@ $^ $(LDLIBS)
	@echo "Note: ./fox reads the text matrix from stdin and prints the result to stdout,"
	@echo "      as the assignment asks. We also implemented reading and writing a binary"
	@echo "      format (./fox_bin input.bin output.bin, see README), which is faster and"
	@echo "      uses less memory; the tests and benchmarks use it."

fox_bin: $(FOX_BIN_OBJ) $(COMMON_OBJ)
	$(MPICC) $(CFLAGS) -o $@ $^ $(LDLIBS)

# Shared variant: same Fox steps, blocks read from a full local copy instead
# of being communicated (lib/mpi/fox.h). Only apsp_mpi.c differs.
fox_shared: $(FOX_SHARED_OBJ) $(COMMON_OBJ)
	$(MPICC) $(CFLAGS) -o $@ $^ $(LDLIBS)

seq: $(SEQ_OBJ) $(COMMON_OBJ)
	$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS)

# src/main.c is compiled three times: with MPI and text I/O (fox), with MPI and
# binary I/O (fox_bin, fox_shared), and without MPI (seq).
$(BUILD)/fox/main_text.o: src/main.c
	@mkdir -p $(dir $@)
	$(MPICC) $(CPPFLAGS) $(CFLAGS) -DUSE_MPI -DTEXT_IO -c $< -o $@

$(BUILD)/fox/main.o: src/main.c
	@mkdir -p $(dir $@)
	$(MPICC) $(CPPFLAGS) $(CFLAGS) -DUSE_MPI -c $< -o $@

$(BUILD)/fox_shared/apsp_mpi.o: lib/mpi/apsp_mpi.c
	@mkdir -p $(dir $@)
	$(MPICC) $(CPPFLAGS) $(CFLAGS) -DFOX_SHARED -c $< -o $@

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

check-fox: fox_bin examples-bin
	@for np in $(NP); do \
	  BIN=./fox_bin NP=$$np MPIRUN="$(MPIRUN)" ./tests/run_tests.sh || exit 1; \
	done

check: test check-seq check-fox

bench: seq fox_bin examples-bin
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

# Report: graphs from the CSVs in bench/report/, then the PDF (its tables read
# the same CSVs). Needs python3 with matplotlib, and typst.
report:
	python3 bench/plot_lab.py bench/report/lab_examples.csv report/lab_examples_speedup.png "Fox on LabCC, assignment inputs"
	python3 bench/plot_lab.py bench/report/mac_examples.csv report/mac_examples_speedup.png "Fox on an Apple M1 Max, assignment inputs"
	python3 bench/plot_lab.py bench/report/lab_rand0.csv report/lab_rand0_speedup.png "Fox on LabCC, random matrices"
	python3 bench/plot_power2.py bench/report/results_power2.csv report/power2_fox_vs_shared.png
	typst compile --root . report/report.typ report/report.pdf

# Delivery zip: sources, build and test scripts, the small examples, the
# report PDF (in report/) and the CSVs it uses. No build output, large inputs or other results.
DIST     := fox-apsp-mpi
DIST_SRC := Makefile README.md setup.md hostfile.txt src lib \
            $(wildcard tests/*.c tests/*.h tests/*.sh) tests/examples/to_bin.sh \
            tests/examples/input/input5 tests/examples/input/input6 \
            tests/examples/output/output5 tests/examples/output/output6 \
            $(wildcard bench/*.sh bench/*.py) bench/report report/report.pdf

dist: report
	rm -rf $(BUILD)/$(DIST) $(DIST).zip
	mkdir -p $(BUILD)/$(DIST)
	cp -r --parents $(DIST_SRC) $(BUILD)/$(DIST)
	cd $(BUILD) && zip -qr ../$(DIST).zip $(DIST)
	@echo "dist: $(DIST).zip"; unzip -l $(DIST).zip | tail -1

clean:
	rm -rf $(BUILD) seq fox fox_bin fox_shared

-include $(shell find $(BUILD) -name '*.d' 2>/dev/null)
