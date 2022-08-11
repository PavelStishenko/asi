export TESTING_PYTHON=~/prg/regr/ve/bin/python3
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:~/prg/aims/build-so-GOO/install/lib/
export ASI_LIB_PATH=`ls -tr ~/prg/aims/build-so-GOO/install/lib/libaims.*so | tail`
export AIMS_SPECIES_DIR=~/prg/aims/species_defaults/defaults_2010/tight/
# not used export ASI_HOME=~/prg/asi/install 
export ASI_INCLUDE_DIR=~/prg/asi/include
export AIMS_DIR=~/prg/aims/build-so-GOO/install/lib
export AIMS_LIB="-laims"
export CC=mpicxx

export PS1="(aims) $PS1"
