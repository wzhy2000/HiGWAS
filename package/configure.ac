## gwas.lasso package,  Zhong Wang<zw355@cornell.edu>
##
## This scipt is a modification of the configure.in file
## from the magma package (version 0.2.1)
## http://cran.r-project.org/web/packages/magma/
## The copyright information for the magma package is as follows:
## Copyright (C) 2010 Brian J Smith
## Licensed under GNU GPL 3 or later

AC_INIT([gwas.lasso], 0.0.1)
#AC_CONFIG_AUX_DIR([tools])

#################### GCC Compiler ####################

AC_PROG_CC(gcc)

COMPILE_CUDA=""
LINKER=`"${R_HOME}/bin/R" CMD config CC`
CUDA_LIBS=""

#################### CUDA Toolkit ####################

AC_ARG_WITH([cuda-home],
	[AC_HELP_STRING([--with-cuda-home=DIR],
		[full path to CUDA installation @<:@/usr/local/cuda@:>@])],
	[CUDA_HOME="$withval"], [CUDA_HOME=/usr/local/cuda])
AC_MSG_RESULT([setting CUDA home directory to ${CUDA_HOME}])

NVCC=${CUDA_HOME}/bin/nvcc
AC_CHECK_FILE([${NVCC}],, [AC_MSG_WARN([nvcc not found!])])
AC_CHECK_FILE([${NVCC}], [use_cuda=YES])

NCFLAGS="-arch=sm_50 -O2"
CUDA_INC=${CUDA_HOME}/include
AC_MSG_CHECKING([for CUDA headers])
AC_CHECK_FILE([${CUDA_INC}/cublas.h],
	[NCFLAGS="${NCFLAGS} -I${CUDA_INC} -DUSECUDA"],
	[AC_MSG_WARN([CUDA headers not found!])]
)

if test "${use_cuda}" == "YES"; then

AC_MSG_CHECKING([whether this is a 64-bit version of CUDA])
AC_CHECK_FILE([${CUDA_HOME}/lib64/libcublas.so], [check_bit=64-bit])
if test "${check_bit}" == "64-bit"; then
	BITS="64"
	LDFLAGS="-I${CUDA_INC} -L${CUDA_HOME}/lib64 -DUSECUDA"
else
	BITS="32"
	LDFLAGS="-I${CUDA_INC} -L${CUDA_HOME}/lib -DUSECUDA"
fi

COMPILE_CUDA="-I${CUDA_INC} -DUSECUDA -m${BITS}"
LINKER="${NVCC} $(CUDA_INC)"
CUDA_LIBS="-lcuda -lcublas -lcublas_device -lcudadevrt -arch=sm_50"

fi



#################### R HOME LIBRARY ####################
: ${R_HOME=`R RHOME`}
AC_ARG_WITH([r-home],
	[AC_HELP_STRING([--with-r-home=DIR],
		[full path to R installation @<:@@:>@])],
	[R_HOME=$withval]
)
AC_MSG_RESULT([setting R home directory to ${R_HOME}])
AC_CHECK_PROG([check_R], [R], [yes],, [${R_HOME}/bin])
if test -z "${check_R}"; then
	AC_MSG_ERROR([R not found!])
fi 


AC_MSG_RESULT([setting R header files to ${R_INCLUDE_DIR}])
R_CPICFLAGS=`"${R_HOME}/bin/R" CMD config CPICFLAGS`
NCFLAGS="${NCFLAGS} -Xcompiler \"-O2 ${R_CPICFLAGS} -I${R_INCLUDE_DIR}\" -dlink -lcublas"


#################### SET FLAGS ####################

COMPILER_FLAGS="${COMPILE_CUDA} -I${R_INCLUDE_DIR} -Wall -fPIC -D_GNU_SOURCE "

COMPILER_OPTIMIZE_FLAGS="-funroll-loops -fomit-frame-pointer -ffast-math"
COMPILER_FLAGS="-O2 ${COMPILER_OPTIMIZE_FLAGS} ${COMPILER_FLAGS}"

CFLAGS="-std=gnu99 ${COMPILER_FLAGS}"
CXXFLAGS="${COMPILER_FLAGS}"

AC_SUBST([NVCC])
AC_SUBST([LINKER])
AC_SUBST([CFLAGS])
AC_SUBST([CXXFLAGS])
AC_SUBST([NCFLAGS])
AC_SUBST([LDFLAGS])
AC_SUBST([CUDA_LIBS])
AC_SUBST([COMPILE_CUDA])

AC_CONFIG_FILES([src/Makefile])
AC_OUTPUT
