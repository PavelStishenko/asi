DFTBP_INCLUDE ?= ${HOME}/opt/dftbp-mpi/include
DFTBP_LIB_DIR ?= ${HOME}/opt/dftbp-mpi/lib
INSTALL_PREFIX ?= ${HOME}/opt/asi
BUILD_PATH ?= ${PWD}/build

all : ${BUILD_PATH}/libasidftbp.so 

${BUILD_PATH}/asidftbp.o : src/dftbp/asi.cpp include/*.h
	mkdir -p $(BUILD_PATH)
	mpicxx -c -std=c++11 -fPIC -I./include -I${DFTBP_INCLUDE}  src/dftbp/asi.cpp -o ${BUILD_PATH}/asidftbp.o

${BUILD_PATH}/libasidftbp.so : ${BUILD_PATH}/asidftbp.o
	mpicxx -shared -Wl,--no-undefined -L${DFTBP_LIB_DIR} -Wl,-start-group -ldftbplus ${BUILD_PATH}/asidftbp.o  -Wl,-end-group -o ${BUILD_PATH}/libasidftbp.so

clean : 
	rm ${BUILD_PATH}/*.o ${BUILD_PATH}/lib*.so

install : ${BUILD_PATH}/libasidftbp.so include/asi.h
	install -d ${INSTALL_PREFIX}/lib
	install -d ${INSTALL_PREFIX}/include
	cp ${BUILD_PATH}/lib*.so ${INSTALL_PREFIX}/lib
	cp include/asi.h ${INSTALL_PREFIX}/include

