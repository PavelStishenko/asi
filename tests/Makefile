BUILD_PATH ?= ${PWD}/build
ASI_HOME ?= ${HOME}/opt/asi
ASI_LIB_DIR ?= ${ASI_HOME}/lib/
ASI_INCLUDE_DIR ?= ${ASI_HOME}/include/


AIMS_DIR ?= ${HOME}/prg/aims/build-so-3
AIMS_LIB ?= -laims.220309.scalapack.mpi

DFTBP_LIB_DIR ?= ${HOME}/opt/dftbp-mpi/lib
DFTBP_LIBS=-lscalapack-openmpi -ldftbplus -lmudpack -lscalapackfx -lmpifx


TEST_SOURCES = $(wildcard src/test_*.cpp)
TEST_OBJECTS = $(patsubst src/%.cpp,$(BUILD_PATH)/%.o,$(TEST_SOURCES))
AIMS_TESTS := $(patsubst src/test_%.cpp,$(BUILD_PATH)/test_%-aims.x,$(TEST_SOURCES))
DFTBP_TESTS := $(patsubst src/test_%.cpp,$(BUILD_PATH)/test_%-dftbp.x, $(filter-out $(wildcard src/*initdm* src/*dm_min*), $(TEST_SOURCES)))

CC ?= mpicxx
COMPILE_CMD := $(CC) -c -g -std=c++14 -fPIC -I./include/ -I$(BOOST_DIR)/include
COMPILE_TEST_CMD := $(COMPILE_CMD) -I$(ASI_INCLUDE_DIR)
COMPILE_AIMS_CMD  := $(COMPILE_TEST_CMD)
COMPILE_DFTBP_CMD := $(COMPILE_TEST_CMD)
LINK_AIMS_CMD  := $(CC) -Wl,-copy-dt-needed-entries -L$(AIMS_DIR) -Wl,-no-as-needed -Wl,-start-group $(AIMS_LIB) $(BUILD_PATH)/make_config_aims.o
LINK_DFTBP_CMD := $(CC) -Wl,-copy-dt-needed-entries -L$(DFTBP_LIB_DIR) -L$(ASI_LIB_DIR) -Wl,-start-group $(DFTBP_LIBS) -lasidftbp $(BUILD_PATH)/make_config_dftbp.o


all : $(AIMS_TESTS) $(DFTBP_TESTS)
	@echo "AIMS_TESTS: ${AIMS_TESTS}"
	@echo "DFTBP_TESTS: ${DFTBP_TESTS}"

$(BUILD_PATH) :
	mkdir -p $(BUILD_PATH)

dftbp_tests : $(DFTBP_TESTS)

aims_tests : $(AIMS_TESTS)

$(BUILD_PATH)/make_config_aims.o : $(BUILD_PATH) src/make_config_aims.cpp include/codespec.hpp
	$(COMPILE_CMD) src/make_config_aims.cpp  -o $(BUILD_PATH)/make_config_aims.o

$(BUILD_PATH)/make_config_dftbp.o : $(BUILD_PATH) src/make_config_dftbp.cpp include/codespec.hpp
	$(COMPILE_CMD) src/make_config_dftbp.cpp -o $(BUILD_PATH)/make_config_dftbp.o

$(TEST_OBJECTS): $(BUILD_PATH)/%.o : src/%.cpp $(BUILD_PATH)
	$(COMPILE_TEST_CMD) $< -o $@

$(AIMS_TESTS): $(BUILD_PATH)/test_%-aims.x : $(BUILD_PATH)/test_%.o $(BUILD_PATH)/make_config_aims.o
	$(LINK_AIMS_CMD) $< -Wl,-end-group -o $@

$(DFTBP_TESTS): $(BUILD_PATH)/test_%-dftbp.x : $(BUILD_PATH)/test_%.o $(BUILD_PATH)/make_config_dftbp.o
	$(LINK_DFTBP_CMD) $< -Wl,-end-group -o $@
	
clean : 
	rm ${BUILD_PATH}/*.o ${BUILD_PATH}/*.x

