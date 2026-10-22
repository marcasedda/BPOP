# Independent build of the local Population III variant.
CXX ?= g++
CXXFLAGS ?= -O3 -std=c++17

PROJECT_ROOT := $(abspath ../..)
BUILD_DIR := .build
RESULT := BPOP.x
OBJECTS := $(BUILD_DIR)/Main_Bpop_dev_p3.o $(BUILD_DIR)/funcGR.o
BUILD_INCLUDES := -I$(BUILD_DIR) -I$(PROJECT_ROOT)/include
TABLES := $(PROJECT_ROOT)/include/redshift_time.txt \
	$(wildcard $(PROJECT_ROOT)/include/ClustersRhM/scottrule_*.ttt)
DATAFILES := data_1e-11.dat spectrum1e-11.txt spectrum_cleaned00000000001.txt

.PHONY: all clean run
all: $(RESULT)

$(RESULT): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(OBJECTS) $(LDLIBS) -o $@

$(BUILD_DIR):
	mkdir -p $@

# Absolute repository location plus the data directory relative to that root.
$(BUILD_DIR)/build_paths.h: Makefile | $(BUILD_DIR)
	printf '%s\n' '#define BPOP_PROJECT_DIR "$(PROJECT_ROOT)/"' \
		'#define BPOP_DATA_SUBDIR "$(patsubst $(PROJECT_ROOT)/%,%,$(CURDIR))/"' > $@

$(BUILD_DIR)/Main_Bpop_dev_p3.o: Main_Bpop_dev_p3.cpp input_params_p3.h $(BUILD_DIR)/build_paths.h Makefile
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(BUILD_INCLUDES) -MMD -MP -c $< -o $@

$(BUILD_DIR)/funcGR.o: $(PROJECT_ROOT)/src/funcGR.cpp Makefile | $(BUILD_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(BUILD_INCLUDES) -MMD -MP -c $< -o $@

run: $(RESULT)
	@set -eu; \
	for input in $(DATAFILES) $(TABLES); do \
		test -s "$$input" || { echo "Missing or empty input: $$input" >&2; exit 1; }; \
	done; \
	for env in YC GC NC; do \
		for kind in mass radi; do \
			test -s "$(PROJECT_ROOT)/include/ClustersRhM/scottrule_$$kind$$env.ttt" || \
				{ echo "Missing cluster table: $$kind$$env" >&2; exit 1; }; \
		done; \
	done; \
	mkdir -p "$(PROJECT_ROOT)/tmp"; \
	run_dir=$$(mktemp -d "$(PROJECT_ROOT)/tmp/pop3-run.XXXXXX"); \
	cp $(RESULT) input_params_p3.h $(BUILD_DIR)/build_paths.h "$$run_dir/"; \
	echo "Run directory: $$run_dir"; \
	cd "$$run_dir"; \
	./$(RESULT) > Output.log 2>&1 || \
		{ mv Output.log Output.txt; echo "Run failed; see $$run_dir/Output.txt" >&2; exit 1; }; \
	mv Output.log Output.txt; \
	grep -q '^Elapsed time = ' Output.txt || \
		{ echo "Run exited before completion; see $$run_dir/Output.txt" >&2; exit 1; }; \
	echo "Completed. Log: $$run_dir/Output.txt"

# Preserve input datasets and all completed run directories.
clean:
	rm -f $(OBJECTS) $(OBJECTS:.o=.d) $(BUILD_DIR)/build_paths.h $(RESULT)
	@rmdir $(BUILD_DIR) 2>/dev/null || true

-include $(OBJECTS:.o=.d)
