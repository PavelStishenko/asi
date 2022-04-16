DFTBP_INCLUDE ?= ${HOME}/opt/dftbp/include
INSTALL_PREFIX ?= ${HOME}/opt/asi
BUILD_PATH ?= ${PWD}/build

all : build_dir ${BUILD_PATH}/asidftbp.so

build_dir:
	mkdir -p $(BUILD_PATH)

${BUILD_PATH}/asidftbp.o : src/dftbp/asi.cpp include/*.h
	mpicxx -c -std=c++11 -fPIC -I./include -I${DFTBP_INCLUDE}  src/dftbp/asi.cpp -o ${BUILD_PATH}/asidftbp.o

${BUILD_PATH}/asidftbp.so : ${BUILD_PATH}/asidftbp.o
	mpicxx -shared ${BUILD_PATH}/asidftbp.o -o ${BUILD_PATH}/libasidftbp.so

clean : 
	rm ${BUILD_PATH}/*.o ${BUILD_PATH}/lib*.so

install : ${BUILD_PATH}/asidftbp.so include/asi.h
	install -D ${BUILD_PATH}/lib*.so ${INSTALL_PREFIX}/lib
	install -d ${INSTALL_PREFIX}/include
	cp include/asi.h ${INSTALL_PREFIX}/include

