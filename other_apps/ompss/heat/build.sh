#!/bin/bash

VERSION=$1
ADD_COMPILE_FLAGS=$2

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
DIR=$DIR/../../

export CFLAGS=" ${ADD_COMPILE_FLAGS}"
export CXXFLAGS=" ${ADD_COMPILE_FLAGS}"

# Check if parsec hooks should be used.
# Note that you will need to compile them if you want to use them (utils/hooks)
if [[ $ADD_COMPILE_FLAGS == *"ENABLE_PARSEC_HOOKS"* ]]
then
    export CFLAGS="$CFLAGS -I${DIR}/utils/hooks/include"
    export CXXFLAGS="$CXXFLAGS -I${DIR}/utils/hooks/include"
    export LDFLAGS="-L${DIR}/utils/hooks/lib"
    export LIBS="-lhooks"
fi

if [ ${VERSION} = "ompss" ]; then
	echo -e "\033[32mCleaning directory\033[m"
	make clean
	echo -e "\033[32mCompiling ${VERSION} version\033[m"
	make heat
	echo -e "\033[32mInstalling ${VERSION} version\033[m"
	make install
	echo -e "\033[32mFinal clean ${VERSION} version\033[m"
	make distclean
	echo -e "\033[32mDone!\033[m"
elif [ ${VERSION} = "serial" ]; then
	echo -e "\033[32mCleaning directory\033[m"
	make clean
	echo -e "\033[32mCompiling ${VERSION} version\033[m"
	make heat-seq
	echo -e "\033[32mInstalling ${VERSION} version\033[m"
	make install
	echo -e "\033[32mFinal clean ${VERSION} version\033[m"
	make distclean
	echo -e "\033[32mDone!\033[m"
else
	echo -e "\033[31m${VERSION} version not supported!\033[m"
fi
