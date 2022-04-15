
# for Intel compiler on Hawk
#AIMS_DIR=/home/c.sacps4/aims/build.so
#AIMS_LIB=-laims.210513.scalapack.mpi
#DFTBP_HOME=/home/c.sacps4/dftbp/lib64 
#DFTBP_LIB_DIR=${DFTBP_HOME}/lib64 
#LIBS='-lmkl_rt'


BUILD_PATH ?= ${PWD}/build
ASI_HOME ?= ${HOME}/opt/asi
ASI_LIB_DIR ?= ${ASI_HOME}/lib/
ASI_INCLUDE_DIR ?= ${ASI_HOME}/include/


AIMS_DIR=${HOME}/prg/aims/build-so-3
AIMS_LIB=-laims.220309.scalapack.mpi

DFTBP_LIB_DIR ?= ${HOME}/opt/dftbp-mpi/lib
DFTBP_LIBS=-lscalapack-openmpi -ldftbplus -lmudpack -lscalapackfx -lmpifx


TEST_SOURCES = $(wildcard src/test_*.cpp)
TEST_OBJECTS = $(patsubst src/%.cpp,$(BUILD_PATH)/%.o,$(TEST_SOURCES))
AIMS_TESTS := $(patsubst src/test_%.cpp,$(BUILD_PATH)/test_%-aims.x,$(TEST_SOURCES))
DFTBP_TESTS := $(patsubst src/test_%.cpp,$(BUILD_PATH)/test_%-dftbp.x, $(filter-out $(wildcard src/*initdm* src/*dm_min*), $(TEST_SOURCES)))

COMPILE_CMD := mpicxx -c -g -std=c++14 -fPIC -I./include/
COMPILE_TEST_CMD := $(COMPILE_CMD) -I$(ASI_INCLUDE_DIR)
COMPILE_AIMS_CMD  := $(COMPILE_TEST_CMD)
COMPILE_DFTBP_CMD := $(COMPILE_TEST_CMD)
LINK_AIMS_CMD  := mpicxx -Wl,-copy-dt-needed-entries -L$(AIMS_DIR) $(AIMS_LIB) $(BUILD_PATH)/make_config_aims.o
LINK_DFTBP_CMD := mpicxx -Wl,-copy-dt-needed-entries -L$(DFTBP_LIB_DIR) -L$(ASI_LIB_DIR) $(DFTBP_LIBS) -lasidftbp $(BUILD_PATH)/make_config_dftbp.o


all : build_dir $(BUILD_PATH)/make_config_aims.o $(BUILD_PATH)/make_config_dftbp.o $(TEST_OBJECTS) $(AIMS_TESTS) $(DFTBP_TESTS)
	@echo "AIMS_TESTS: ${AIMS_TESTS}"
	@echo "DFTBP_TESTS: ${DFTBP_TESTS}"

build_dir:
	mkdir -p $(BUILD_PATH)


$(BUILD_PATH)/make_config_aims.o : src/make_config_aims.cpp include/codespec.hpp
	$(COMPILE_CMD) src/make_config_aims.cpp  -o $(BUILD_PATH)/make_config_aims.o

$(BUILD_PATH)/make_config_dftbp.o : src/make_config_dftbp.cpp include/codespec.hpp
	$(COMPILE_CMD) src/make_config_dftbp.cpp -o $(BUILD_PATH)/make_config_dftbp.o

$(TEST_OBJECTS): $(BUILD_PATH)/%.o : src/%.cpp
	$(COMPILE_TEST_CMD) $< -o $@

$(AIMS_TESTS): $(BUILD_PATH)/test_%-aims.x : $(BUILD_PATH)/test_%.o
	$(LINK_AIMS_CMD) $< -o $@

$(DFTBP_TESTS): $(BUILD_PATH)/test_%-dftbp.x : $(BUILD_PATH)/test_%.o
	$(LINK_DFTBP_CMD) $< -o $@
	
clean : 
	rm ${BUILD_PATH}/*.o ${BUILD_PATH}/*.x

