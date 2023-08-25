
# linux / gcc options
LIB = .so
LDFLAGS = -ldl
SHLIBFLAGS = -shared
FC=gfortran
# FFLAGS= -g
FFLAGS = -O2 -Wall -fPIC

# SOCRATES library code must be built with -fPIC (eg edit Mk_cmd to add this option)
# library location RAD_BIN is set by set_rad_env script in the top-level socrates working folder

# FDEPS = paleo_fortran.mod
FOBJ = Utilities_CF.o StrCtrl_C.o StrDim_C.o StrSpecData_C.o StrAtm_C.o StrCld_C.o StrAer_C.o StrBound_C.o StrOut_C.o SOCRATES_C.o \
       StrSpecDataDim_C.o  StrSpecDataBasic_C.o StrSpecDataSolar_C.o StrSpecDataRayleigh_C.o StrSpecDataGas_C.o StrSpecDataPlanck_C.o \
	   StrSpecDataCont_C.o StrSpecDataContGen_C.o StrSpecDataDrop_C.o \
       StrSpecDataAerosol_C.o StrSpecDataIce_C.o StrSpecDataVar_C.o StrSpecDataPhotol_C.o StrSpecDataMap_C.o

#=============================================================================

# use -I flag to tell gfortran where to look for SOCRATES .mod files
%.o: ../src/%.f90 
	$(FC) -c -o $@ $< $(FFLAGS) -I$(RAD_BIN)

%.o: ../gen/%.f90 
	$(FC) -c -o $@ $< $(FFLAGS) -I$(RAD_BIN)

libSOCRATES_C$(LIB): $(FOBJ)
	$(FC) $(SHLIBFLAGS) $(FFLAGS) -o $@ $^ $(RAD_BIN)/radlib.a

clean:
	rm -f *.o *.mod libSOCRATES_C$(LIB)

.PHONY: clean
