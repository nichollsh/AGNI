# Script to generate Mk_cmd for socrates

genpath="Mk_cmd_gen"

rm -f $genpath

echo "# Generated automatically" >> $genpath 
echo "# System: $(uname -a) " >> $genpath
echo "# Date: $(date) " >> $genpath
echo " " >> $genpath

fc="$(which gfortran)"
echo "FORTCOMP        = $fc -O3 -fPIC -c " >> $genpath
echo "LINK            = $fc -O3 -fPIC " >> $genpath

ar="$(which ar)"
echo "LIBLINK         = $ar rvu " >> $genpath

cdfinc="$(nf-config --includedir)"
echo "INCCDF_PATH     = $cdfinc " >> $genpath

cdflib="$(nf-config --flibs | cut -d " " -f 1 | tail -c +3)"
echo "LIBCDF_PATH     = $cdflib " >> $genpath
echo "LIBCDF_NAME     = netcdff " >> $genpath
echo "OMPARG          = -fopenmp " >> $genpath

echo " " >> $genpath 
echo "LIBSUFFIX       = a " >> $genpath

echo " " >> $genpath
echo ".SUFFIXES: \$(SUFFIXES) .f90 " >> $genpath 

echo " " >> $genpath
