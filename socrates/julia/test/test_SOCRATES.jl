
# Path to SOCRATES svn top level: Edit for path on local computer
SOCRATES_DIR = "../.."

import SOCRATES


##################################
# Test minimal argument passing
##################################

din = 42.0

dout = SOCRATES.test_double_val(din)
println("SOCRATES.test_double_val ", din, " --> ", dout)

dout = SOCRATES.test_double_ref()
println("SOCRATES.test_double_ref C_NULL --> ", dout)
dout = SOCRATES.test_double_ref(din)
println("SOCRATES.test_double_ref ", din, " --> ", dout)

##############################################
# Test create and delete of SOCRATES structs
# (handles containing Ptr{Cvoid} to Fortran types)
##############################################

control = SOCRATES.StrCtrl()
println("StrCtrl: ", control)
SOCRATES.dump_properties(control)
println()

spectrum = SOCRATES.StrSpecData()
println("StrSpecData: ", spectrum)


#################################
# Test reading a spectral file
#################################

# SOCRATES data and example folders relative to SOCRATES_DIR
RAD_DATA = joinpath(SOCRATES_DIR, "data")

spectral_file = joinpath(RAD_DATA, "spectra/ga7/sp_lw_ga7")

spectrum = SOCRATES.StrSpecData()
println("StrSpecData: ", spectrum)

SOCRATES.read_spectrum(spectral_file, spectrum)
# SOCRATES.set_spectrum(Ptr_StrSpecData=spectrum, spectral_file=spectral_file, l_all_gasses=true)

SOCRATES.dump_spectrum(spectrum)
println()

#################################
# Test creating a StrDim
#################################

dimen = SOCRATES.StrDim()
println("StrDim: ")
SOCRATES.dump_properties(dimen)
println()

# test modifying a field
dimen.nd_tile = 2
println("StrDim: (after modifying nd_tile = 2) ")
SOCRATES.dump_properties(dimen)


