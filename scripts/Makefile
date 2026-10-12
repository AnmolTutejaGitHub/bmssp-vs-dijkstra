# Builds three binaries from the same sources:
#   bin/sssp        timing build (no instrumentation)
#   bin/sssp_count  operation-count build (-DCOUNT_OPS)
#   bin/sssp_debug  assertions + sanitizers, used by `make check`
CXX      ?= c++
CXXFLAGS ?= -std=c++17 -O3 -DNDEBUG -Wall -Wextra
# On macOS, point the compiler at the libc++ headers of the active SDK (some
# Command Line Tools installations ship an incomplete copy next to the compiler).
SDKCXX   := $(shell xcrun --show-sdk-path 2>/dev/null)/usr/include/c++/v1
ifneq ($(wildcard $(SDKCXX)/vector),)
CXXFLAGS += -isystem $(SDKCXX)
DBGINC   := -isystem $(SDKCXX)
endif
SRC      := src/main.cpp
HDR      := $(wildcard src/*.hpp)

all: bin/sssp bin/sssp_count

bin/sssp: $(SRC) $(HDR)
	@mkdir -p bin
	$(CXX) $(CXXFLAGS) -o $@ $(SRC)

bin/sssp_count: $(SRC) $(HDR)
	@mkdir -p bin
	$(CXX) $(CXXFLAGS) -DCOUNT_OPS -o $@ $(SRC)

bin/sssp_debug: $(SRC) $(HDR)
	@mkdir -p bin
	$(CXX) -std=c++17 -O1 -g -DCOUNT_OPS -fsanitize=address,undefined -fno-omit-frame-pointer -Wall -Wextra $(DBGINC) -o $@ $(SRC)

check: bin/sssp_debug
	./bin/sssp_debug validate --trials 300 --seed 7 --maxn 40

clean:
	rm -rf bin

.PHONY: all check clean
