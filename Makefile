# Seam carving - build, test, benchmark (GNU make; on Windows use mingw32-make, or double-click build_all.bat).
# Nothing is left behind: each target compiles its program into the system temp folder, runs it and deletes it.
# Flags: keep in sync with tools/mutation_test.py.

ifeq ($(origin CXX),default)
  CXX := g++
endif
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra -pthread -Iinclude -Ithird_party

ifeq ($(OS),Windows_NT)
  EXE := .exe
  # Static linking avoids the Git Bash libstdc++ DLL clash.
  LDFLAGS += -static
else
  EXE :=
endif

# First working Python: python3, python or py.
ifeq ($(OS),Windows_NT)
  PYTHON ?= $(firstword $(foreach p,python3 python py,$(if $(filter 1,$(shell $(p) -c "print(1)" 2>&1)),$(p))))
else
  PYTHON ?= $(firstword $(shell command -v python3 python 2>/dev/null))
endif

TMP     := $(shell $(PYTHON) -c "import tempfile, pathlib; print(pathlib.Path(tempfile.gettempdir()).as_posix())")
# A fresh name per run, so two runs at once (or a leftover file) can never block each other.
UID     := $(shell $(PYTHON) -c "import uuid; print(uuid.uuid4().hex[:8])")
TESTBIN := $(TMP)/sc_run_tests_$(UID)$(EXE)
TOOLBIN := $(TMP)/sc_seamcarve_$(UID)$(EXE)

# RUN <program> <args>: run the program, delete it (pass or fail), then exit with its exit code.
RUN := $(PYTHON) -c "import subprocess, os, sys; r = subprocess.call(sys.argv[1:]); os.remove(sys.argv[1]); sys.exit(r)"

.PHONY: all test bench bench-quick demo viz results mutation clean

# Plain "make" compiles and runs every test, printing each result in the terminal.
# One module only: make test ARGS="--filter DP"
all: test

test:
	$(CXX) $(CXXFLAGS) tests/run_tests.cpp -o $(TESTBIN) $(LDFLAGS)
	$(RUN) $(TESTBIN) $(ARGS)

bench:
	$(CXX) $(CXXFLAGS) src/main.cpp -o $(TOOLBIN) $(LDFLAGS)
	$(RUN) $(TOOLBIN) bench --suite all --out results

bench-quick:
	$(CXX) $(CXXFLAGS) src/main.cpp -o $(TOOLBIN) $(LDFLAGS)
	$(RUN) $(TOOLBIN) bench --suite all --quick --out results

demo:
	$(CXX) $(CXXFLAGS) src/main.cpp -o $(TOOLBIN) $(LDFLAGS)
	$(RUN) $(TOOLBIN) demo --out results/demo

viz:
	$(PYTHON) tools/make_viz.py

results: test demo bench viz

mutation:
	$(PYTHON) tools/mutation_test.py

clean:
	$(PYTHON) -c "import shutil; shutil.rmtree('tests_tmp', True)"
