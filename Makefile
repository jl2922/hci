# Default options.
CXX := mpic++
CXX_WARNING_OPTIONS := -Wall -Wextra -Wno-expansion-to-defined -Wno-int-in-bool-context
CXX_MALLOC_OPTIONS := -fno-builtin-malloc -fno-builtin-calloc -fno-builtin-realloc -fno-builtin-free
CXX_ARCH_TUNE := -march='native' -mtune='native' # Needs recompile on different CPU models.
CXXFLAGS := -std=c++17 -O3 -fopenmp $(CXX_ARCH_TUNE) $(CXX_WARNING_OPTIONS) $(CXX_MALLOC_OPTIONS)
LDLIBS := -pthread -lboost_mpi -lboost_serialization -lprotobuf -lpthread
SRC_DIR := src
OBJ_DIR := build
EXE := hci.x
TEST_EXE := $(OBJ_DIR)/hci_test.x

# Libraries.
UNAME := $(shell uname)
HOSTNAME := $(shell hostname)
ifeq ($(UNAME), Linux)
	TOOLS_DIR := $(HOME)/tools
	EIGEN_DIR := $(TOOLS_DIR)/eigen
	BOOST_DIR := $(TOOLS_DIR)/boost
	PROTOBUF_DIR := $(TOOLS_DIR)/protobuf
	GPERFTOOLS_DIR := $(TOOLS_DIR)/gperftools
	CXXFLAGS := $(CXXFLAGS) -I $(EIGEN_DIR)/include -I $(BOOST_DIR)/include -I $(PROTOBUF_DIR)/include
	LDLIBS := -L $(BOOST_DIR)/lib -L $(GPERFTOOLS_DIR)/lib -L $(PROTOBUF_DIR)/lib $(LDLIBS) -ltcmalloc
endif
ifeq ($(UNAME), Darwin)
	LDLIBS := $(LDLIBS) -ltcmalloc_minimal
endif

# Load Makefile.config if exists.
LOCAL_MAKEFILE := local.mk
ifneq ($(wildcard $(LOCAL_MAKEFILE)),)
	include $(LOCAL_MAKEFILE)
endif

# Sources and intermediate objects.
PROTO_SRC := $(SRC_DIR)/data.proto
SRCS := $(shell find $(SRC_DIR) \
		! -name "main.cc" ! -name "*_test.cc" -name "*.cc")
TESTS := $(shell find $(SRC_DIR) -name "*_test.cc")
HEADERS := $(shell find $(SRC_DIR) -name "*.h")
MAIN := $(SRC_DIR)/main.cc
OBJS := $(SRCS:$(SRC_DIR)/%.cc=$(OBJ_DIR)/%.o)
TEST_OBJS := $(TESTS:$(SRC_DIR)/%.cc=$(OBJ_DIR)/%.o)

# Test related.
GTEST_DIR := gtest/googletest
GMOCK_DIR := gtest/googlemock
GTEST_ALL_SRC := ${GTEST_DIR}/src/gtest-all.cc
GMOCK_ALL_SRC := ${GMOCK_DIR}/src/gmock-all.cc
TEST_MAIN_SRC := ${GMOCK_DIR}/src/gmock_main.cc
TEST_MAIN := $(OBJ_DIR)/gtest_main.o
TEST_CXXFLAGS := $(CXXFLAGS) -isystem $(GTEST_DIR)/include -isystem $(GMOCK_DIR)/include -pthread
TEST_LIB := $(OBJ_DIR)/libgtest.a

.PHONY: all test proto clean

all: $(EXE)

test: $(TEST_EXE)
	./$(TEST_EXE) --gtest_filter=-*LargeTest.*

all_tests: $(TEST_EXE)
	./$(TEST_EXE)

proto: $(PROTO_SRC)
	protoc -I=$(SRC_DIR) --cpp_out=$(SRC_DIR) --python_out=. $(PROTO_SRC)

clean:
	rm -rf $(OBJ_DIR)
	rm -f ./$(EXE)
	rm -f ./$(TEST_EXE)

# Main program.

$(EXE): $(OBJS) $(MAIN) $(HEADERS)
	$(CXX) $(CXXFLAGS) $(MAIN) $(OBJS) -o $(EXE) $(LDLIBS)
	
$(OBJS): $(OBJ_DIR)/%.o: $(SRC_DIR)/%.cc $(HEADERS)
	mkdir -p $(@D) && $(CXX) $(CXXFLAGS) -c $< -o $@

# Test specific.
$(OBJ_DIR)/gtest-all.o: $(GTEST_ALL_SRC)
	mkdir -p $(@D) && $(CXX) $(TEST_CXXFLAGS) -I$(GTEST_DIR) -I$(GMOCK_DIR) -c $(GTEST_ALL_SRC) -o $@

$(OBJ_DIR)/gmock-all.o: $(GMOCK_ALL_SRC)
	mkdir -p $(@D) && $(CXX) $(TEST_CXXFLAGS) -I$(GTEST_DIR) -I$(GMOCK_DIR) -c $(GMOCK_ALL_SRC) -o $@

$(TEST_LIB): $(OBJ_DIR)/gtest-all.o $(OBJ_DIR)/gmock-all.o
	$(AR) $(ARFLAGS) $@ $(OBJ_DIR)/gtest-all.o $(OBJ_DIR)/gmock-all.o

$(TEST_OBJS): $(OBJ_DIR)/%.o: $(SRC_DIR)/%.cc $(HEADERS)
	mkdir -p $(@D) && $(CXX) $(TEST_CXXFLAGS) -c $< -o $@

$(TEST_MAIN): $(TEST_MAIN_SRC)
	mkdir -p $(@D) && $(CXX) -I$(GTEST_DIR) -I$(GMOCK_DIR) $(TEST_CXXFLAGS) -c $< -o $@

$(TEST_EXE): $(TEST_OBJS) $(OBJS) $(TEST_MAIN) $(TEST_LIB)
	$(CXX) $(CXXFLAGS) $(TEST_OBJS) $(OBJS) $(TEST_MAIN) $(TEST_LIB) -o $(TEST_EXE) $(LDLIBS) -lpthread
