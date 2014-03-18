#!/usr/bin/env bash

# Takes fortran files that call blas routines and makes
#   new fortran files that call slightly different blas routines.
# This is because the fortran file expected 32-bit blas,
#   but most modern blas is 64-bit for integer types.
# When using Matlab, use -lmwblascompat32 to use 32-bit blas,
#   but this uses names like ddot32 instead of ddot

FILES='lbfgsb.f linpack.f'
DIR='Lbfgsb.3.0'

for file in $FILES
do
    OUTPUT=${file%.f}32.f
    echo Converting $DIR/$file to $DIR/$OUTPUT

    sed 's/ddot/ddot32/' $DIR/$file > $DIR/$OUTPUT
    sed -i 's/dcopy/dcopy32/' $DIR/$OUTPUT
    sed -i 's/dscal/dscal32/' $DIR/$OUTPUT
    sed -i 's/daxpy/daxpy32/' $DIR/$OUTPUT
done
